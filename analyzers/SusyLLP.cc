#include "SusyLLP.h"
#include "RazorHelper.h"
#include "SusyLLPTree.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"

//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"

#define _debug 0
#define _debug_npu 0
#define _debug_sync 0
#define _debug_met 0
#define _debug_jet 0
#define _debug_match 0
#define _debug_avgH 0
#define _debug_trg 0
#define _debug_calojet 0
#define _run_calojet_ 0

#define N_MAX_LLP_DAUGHTERS 4
#define N_MAX_LLP_GRAND_DAUGHTERS 4
#define N_MAX_LLPS 2
#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define NTriggersMAX 602 //Number of trigger in the .dat file
using namespace std;

struct leptons
{
  TLorentzVector lepton;
  int pdgId;
  // float dZ;
  // bool passLooseId;
  // bool passMediumId;
  // bool passTightId;
  // bool passId;
};


struct jets
{
  TLorentzVector jet;
  float time;
  int jetNeutralHadronMultiplicity;
  int jetChargedHadronMultiplicity;
  int jetMuonMultiplicity;
  int jetElectronMultiplicity;
  int jetPhotonMultiplicity;
  float jetNeutralHadronEnergyFraction;
  float jetChargedHadronEnergyFraction;
  float jetMuonEnergyFraction;
  float jetElectronEnergyFraction;
  float jetPhotonEnergyFraction;
  float jetCSV;
  int ecalNRechits;
  float ecalRechitE;

  float jetGammaMax_ET;
  float jetMinDeltaRPVTracks;

  float jetChargedEMEnergyFraction;
  float jetNeutralEMEnergyFraction;
};

struct calojets
{
  TLorentzVector calojet;
  float time;
  bool passId;
  bool isCSVT;
  bool matched;
  int ecalNRechits;
  float ecalRechitE;
  float calojet_EMEnergyFraction;
  float calojet_HadronicEnergyFraction;
  float calojetGammaMax_ET;
  float calojetMinDeltaRPVTracks;
  float calojetPtAllPVTracks;

};

//pt comparison
//not used so far
struct greater_than_pt
{
  inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
};

//lepton highest pt comparator
struct largest_pt_lep
{
  inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt_lep;

//jet highest pt comparator
struct largest_pt_jet
{
  inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;

//calojet highest pt comparator
struct largest_pt_calojet
{
  inline bool operator() (const calojets& p1, const calojets& p2){return p1.calojet.Pt() > p2.calojet.Pt();}
} my_largest_pt_calojet;

//Analyze
//void SusyLLP::Analyze(bool isData, int options, string outputfilename, string analysisTag)
void SusyLLP::Analyze(bool isData, int options, string outputfilename, string analysisTag, string process)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  //-----------option----------
  //---------------------------
  int option;
  std::string label;
  bool pf;

  //HUNDRED'S DIGIT
  //option of run condor or locally
  if (options < 200){
    option = 1; // used when running condor
  }
  else{
    option = 0;// used when running locally
    //options need to be larger than 200, if do local test run
    cout << "option = 0, running locally, load aux locally \n";
  }

  //TEN'S DIGIT
  // label of signal / bkg
  if ((options/10)%10 == 1){
    label = "wH";
    cout << "process label: " << label << "\n";
  }
  else if ((options/10) % 10 == 2){
    label = "zH";
    cout << "process label: " << label << "\n";
  }
  else if ((options/10) % 10 == 3){
    label = "bkg_wH";
    cout << "process label: " << label << "\n";
  }
  else if ((options/10) % 10 == 4){
    label = "bkg_zH";
    cout << "process label: " << label << "\n";
  }
  else if ((options/10) % 10 == 5){
    label = "HH";
    cout << "process label: " << label << "\n";
  }
  else if ((options/10) % 10 == 6){
    label = "bkg_HH";
    cout << "process label: " << label << "\n";
  }
  else{
    cout << "What process it is? Label not defined. \n";
  }

  //UNIT'S DIGIT
  // pf option
  if(options%10==1){
    pf = true;
  }
  else{
    pf = false;
  }

  // DATA or MC
  if( isData )
  {
    std::cout << "[INFO]: running on data with label: " << label << " and option: " << option << " and pfjet is " << pf << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with label: " << label << " and option: " << option << " and pfjet is " << pf << std::endl;
  }

  const float ELE_MASS = 0.000511;
  const float MU_MASS  = 0.105658;
  const float Z_MASS   = 91.2;
  const float H_MASS   = 125.0;
   
  string dataset = "94X";
  //Analysis Tag
  //reference in RazorHelper
  if (analysisTag == ""){
    analysisTag = "Razor2016_80X";
  }

  int wzId;
  int NTrigger;//Number of trigger in trigger paths
  int elePt_cut = 0;
  int muonPt_cut = 0;
  uint nLepton_cut = 0;



  if (label == "HH" || label == "bkg_HH" ){
    NTrigger = 1;
    muonPt_cut = 15;
    elePt_cut = 15;
    nLepton_cut = 0;
    }
  if (label == "zH" || label == "bkg_zH" ){
    NTrigger = 4;
    muonPt_cut = 15;
    elePt_cut = 15;
    nLepton_cut = 2;
    }
  //else{}
  if (label == "wH" || label == "bkg_wH" ){
    NTrigger = 2;
    muonPt_cut = 27;
    elePt_cut = 32;
    nLepton_cut = 1;
  }

  int trigger_paths[NTrigger];
  if (label == "wH" || label == "bkg_wH"){
    wzId = 24;
    trigger_paths[0] = 87;
    trigger_paths[1] = 135;
    // trigger_paths[2] = 310;
  }
  else if (label == "zH" || label == "bkg_zH"){
    wzId = 23;
    trigger_paths[0] = 177;
    trigger_paths[1] = 362;
    // trigger_paths[2] = 310;
    trigger_paths[2] = 87;
    trigger_paths[3] = 135;
  }
  else if (label == "HH" || label == "bkg_HH"){
     wzId = 25;
     trigger_paths[0] = 310;
  }
  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "SusyLLPTree.root";
  TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
  //RazorLiteTree *llp_tree = new RazorLiteTree;
  SusyLLPTree *llp_tree = new SusyLLPTree;
  llp_tree->CreateTree();
  llp_tree->tree_->SetAutoFlush(0);
  llp_tree->InitTree();
  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *generatedEvents = new TH1F("generatedEvents", "generatedEvents", 1, 1, 2);
  TH1F *trig = new TH1F("trig", "trig", 1, 1, 2);
  TH1F *trig_lepId = new TH1F("trig_lepId", "trig_lepId", 1, 1, 2);
  TH1F *trig_lepId_dijet = new TH1F("trig_lepId_dijet", "trig_lepId_dijet", 1, 1, 2);


  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/cms_lpc_llp/llp_analyzer/data/JEC/";
  if(cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1

  cout << "Getting JEC parameters from " << pathname << endl;

  std::vector<JetCorrectorParameters> correctionParameters;
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));

  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);

  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  if (process == ""){
    process = "ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Summer16_2016";
  }

  if (analysisTag == "Razor2015_76X") helper = new RazorHelper("Razor2015_76X", isData, false);
  else if (analysisTag == "Razor2016_MoriondRereco") helper = new RazorHelper("Razor2016_MoriondRereco", isData, false);
  //else if (analysisTag == "CT2016_07Aug2017Rereco") helper = new RazorHelper("CT2016_07Aug2017Rereco", isData, false);
  else if (analysisTag == "CT2016_07Aug2017Rereco") helper = new RazorHelper("CT2016_07Aug2017Rereco", isData, false, process.c_str());
  //else if (analysisTag == "CT2016_07Aug2017Rereco") helper = new RazorHelper("CT2016_07Aug2017Rereco", isData, false, "QCD_HT50toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16_2016");
  //else if (analysisTag == "CT2016_07Aug2017Rereco") helper = new RazorHelper("CT2016_07Aug2017Rereco", isData, false, "ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Summer16_2016");
  else helper = new RazorHelper(analysisTag, isData, false);

  //----------
  //pu histo
  //----------
  //TH1D* puhisto = new TH1D("pileup", "", 50, 0, 50);
  //histogram containing total number of processed events (for normalization)
  //TH1F *histNPV = new TH1F("NPV", "NPV", 2, -0.5, 1.5);
  //TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  //TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
  //TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
  //TH1F *SumPdfWeights = new TH1F("SumPdfWeights", "SumPdfWeights", NUM_PDF_WEIGHTS, -0.5, NUM_PDF_WEIGHTS-0.5);



  //*************************************************************************
  //Look over Input File Events
  //*************************************************************************
  if (fChain == 0) return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {

    //begin event
    if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //GetEntry(ientry);
    nb = fChain->GetEntry(jentry); nbytes += nb;

    //fill normalization histogram
    if(_debug) std::cout << "deb0 " << jentry << std::endl;
    llp_tree->InitVariables();
    if(_debug) std::cout << "deb1 " << jentry << std::endl;
    if (label =="bkg_wH"|| label == "bkg_zH" || label == "bkg_HH"){
      if (isData)
      {
        NEvents->Fill(1);
        llp_tree->weight = 1;
      }
      else
      {
        //NEvents->Fill(genWeight);
        //llp_tree->weight = genWeight;
        NEvents->Fill(1);
        llp_tree->weight = 1;
      }

    }
    else{
      generatedEvents->Fill(1);
      llp_tree->weight = 1;
    }
    if(_debug) std::cout << "deb2 " << jentry << std::endl;
    //event info
    llp_tree->runNum = runNum;
    llp_tree->lumiSec = lumiNum;
    llp_tree->evtNum = eventNum;
    if(_debug) std::cout << "deb3 " << jentry << std::endl;
    if(_debug) std::cout << "nBunchXing " << nBunchXing << std::endl;
    if (label == "zH" || label == "wH" ){
      NEvents->Fill(1);
      bool wzFlag = false;
      for (int i=0; i < nGenParticle; ++i)
      {
        // if (abs(gParticleId[i]) == wzId && gParticleStatus[i] == 22)
        //if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == wzId)
        if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] && abs(gParticleMotherId[i]) == wzId)
        {
          wzFlag = true;
        }

      }
      if ( wzFlag == false ) continue;

    }
    else if (label == "HH"){
      NEvents->Fill(1);
      bool wzFlag = false;
      for (int i=0; i < nGenParticle; ++i)
      {
        if ((abs(gParticleId[i]) == 5) && gParticleStatus[i] && abs(gParticleMotherId[i]) == wzId)
        {
          wzFlag = true;
        }

      }
    //if(_debug) std::cout << "wzFlag " << wzFlag << std::endl;
      //if ( wzFlag == false ) continue;

    }

    if(!isData)
    {
    if(_debug) std::cout << "nBunchXing " << nBunchXing << std::endl;
      for (int i=0; i < nBunchXing; i++)
      {
        if (BunchXing[i] == 0)
        {
          llp_tree->npu = int(nPUmean[i]);
      if(_debug_npu) 
      {
	std::cout << "npu " << llp_tree->npu << std::endl;
        std::cout << "nPUmean[i] " << nPUmean[i] << std::endl;
      }
        }
      }

      llp_tree->pileupWeight = 1;
      llp_tree->pileupWeight = helper->getPileupWeight(llp_tree->npu);
      llp_tree->pileupWeightUp = helper->getPileupWeightUp(llp_tree->npu) / llp_tree->pileupWeight;
      llp_tree->pileupWeightDown = helper->getPileupWeightDown(llp_tree->npu) / llp_tree->pileupWeight;
      if(_debug_npu) 
      {
        std::cout << "pileupWeightUp " << llp_tree->pileupWeightUp << std::endl;
	std::cout << "pileupWeightDown " << llp_tree->pileupWeightDown << std::endl;
      }
    }

    if(_debug_met) std::cout << "npu " << llp_tree->npu << std::endl;
    if(_debug && llp_tree->npu != 0 ) std::cout << "npu " << llp_tree->npu << std::endl;
    //get NPU
    llp_tree->npv = nPV;
    llp_tree->rho = fixedGridRhoFastjetAll;
    llp_tree->met = metType1Pt;
    if(_debug_met) std::cout << "met " << llp_tree->met << std::endl;
    if( llp_tree->met < 120. ) continue;
    //if( llp_tree->met < 150. ) continue;
    if(_debug_met) std::cout << "metType1Pt passed" << metType1Pt << std::endl;
    llp_tree->metPhi = metType1Phi;
    if(_debug) std::cout << "npv " << llp_tree->npv << std::endl;
    TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM( metType1Pt, 0, metType1Phi, 0 );

   //met filters
   llp_tree->Flag2_globalSuperTightHalo2016Filter          = Flag2_globalSuperTightHalo2016Filter;
   llp_tree->Flag2_globalTightHalo2016Filter               = Flag2_globalTightHalo2016Filter;
   llp_tree->Flag2_goodVertices                            = Flag2_goodVertices;
   llp_tree->Flag2_BadChargedCandidateFilter               = Flag2_BadChargedCandidateFilter;
   llp_tree->Flag2_BadPFMuonFilter                         = Flag2_BadPFMuonFilter;
   llp_tree->Flag2_EcalDeadCellTriggerPrimitiveFilter      = Flag2_EcalDeadCellTriggerPrimitiveFilter;
   llp_tree->Flag2_HBHENoiseFilter                         = Flag2_HBHENoiseFilter;
   llp_tree->Flag2_HBHEIsoNoiseFilter                      = Flag2_HBHEIsoNoiseFilter;
   llp_tree->Flag2_ecalBadCalibFilter                      = Flag2_ecalBadCalibFilter;
   llp_tree->Flag2_eeBadScFilter                           = Flag2_eeBadScFilter;

    //Triggers
    for(int i = 0; i < NTriggersMAX; i++){
      llp_tree->HLTDecision[i] = HLTDecision[i];
    }
    if(_debug_trg) std::cout << "begin: 310 " << HLTDecision[310] << std::endl;
    if(_debug_trg) std::cout << "begin: 310 " << llp_tree->HLTDecision[310] << std::endl;
    bool triggered = false;
    for(int i = 0; i < NTrigger; i++)
    {
    if(_debug_trg) std::cout << "i " << i << ", NTrigger "<< NTrigger << std::endl;
   
      int trigger_temp = trigger_paths[i];
    if(_debug_trg) std::cout << "temp  " << trigger_paths[i] << ", triggered "<< triggered << std::endl;

      triggered = triggered || HLTDecision[trigger_temp];
    if(_debug_trg) std::cout << "is triggered ?"<< triggered << std::endl;

    }
    if (triggered) trig->Fill(1);
    if(_debug) std::cout << "triggered " << triggered << std::endl;
    //*************************************************************************
    //Start Object Selection
    //*************************************************************************
    //sync 
    if(_debug_sync)
    {
    std::cout << "nMuons " << nMuons << std::endl;
    std::cout << "nElectron " << nElectrons << std::endl;
    std::cout << "nPhoton " << nPhotons << std::endl;
    std::cout << "nTau " << nTaus << std::endl;
    }


    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
    //twiki (https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation)
    for( int i = 0; i < nMuons; i++ )
    {
       if(_debug_sync)
       {
    	   std::cout << "nMuons " << nMuons << std::endl;
      	   std::cout << "iMuon " << i << ", Pt " << muonPt[i] << ", Eta " << muonEta[i]<< ", Phi " << muonPhi[i]<< ", E " << muonE[i] << std::endl;
      	   std::cout << "iMuon " << i << ", muon_chargedIso[i] " << muon_chargedIso[i] << ", muon_neutralHadIso[i] " << muon_neutralHadIso[i]<< ", muon_photonIso[i] " << muon_photonIso[i]<< ", muon_pileupIso[i] " << muon_pileupIso[i] << std::endl;
      	   std::cout << "iMuon " << i << ", muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i] " << muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i] << ", std::max( stuff, 0)  " << std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) << ", muon_chargedIso[i] + std::max( stuff, 0) " << (muon_chargedIso[i] + std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) ) << ", (muon_chargedIso[i] + std::max( stuff, 0))/muonPt " << (muon_chargedIso[i] + std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) )/muonPt[i] << std::endl;
       }

	if(muonPt[i] <= 10 ) continue;
        if(fabs(muonEta[i]) > 2.4) continue;

        if(!muonIsLoose[i]) continue;
	if( (muon_chargedIso[i] + std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) )/muonPt[i] >= 0.25) continue;


        //remove overlaps
        bool overlap = false;
        for(auto& lep : Leptons)
        {
          if (RazorAnalyzerLLP::deltaR(muonEta[i],muonPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
        //if(overlap) continue;

        leptons tmpMuon;
        tmpMuon.lepton.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i], MU_MASS);
        tmpMuon.pdgId = 13 * -1 * muonCharge[i];

        Leptons.push_back(tmpMuon);

        llp_tree->muon_pileupIso[llp_tree->nMuons] = muon_pileupIso[i];
        llp_tree->muon_chargedIso[llp_tree->nMuons] = muon_chargedIso[i];
        llp_tree->muon_photonIso[llp_tree->nMuons] = muon_photonIso[i];
        llp_tree->muon_neutralHadIso[llp_tree->nMuons] = muon_neutralHadIso[i];

        llp_tree->muonIsLoose[llp_tree->nMuons] = muonIsLoose[i];

        llp_tree->muonPt[llp_tree->nMuons] = muonPt[i];
        llp_tree->muonEta[llp_tree->nMuons] = muonEta[i];
        llp_tree->muonE[llp_tree->nMuons] = muonE[i];
        llp_tree->muonPhi[llp_tree->nMuons] = muonPhi[i];

        llp_tree->nMuons++;
    }
    
    //-------------------------------
    //Electrons
    //-------------------------------
    for( int i = 0; i < nElectrons; i++ )
    {
       if(_debug_sync)
       {
    	   std::cout << "nElectrons " << nElectrons << std::endl;
      	   std::cout << "iElectron " << i << ", Pt " << elePt[i] << ", Eta " << eleEta[i]<< ", Phi " << elePhi[i]<< ", E " << eleE[i] << std::endl;
      	   std::cout << "iElectron " << i << ", ele_passCutBasedIDVeto[i] " << ele_passCutBasedIDVeto[i]<< std::endl;
       }

	if(elePt[i] <= 10 ) continue;
        if(fabs(eleEta[i]) > 2.5) continue;

        if(!ele_passCutBasedIDVeto[i]) continue;

        //remove overlaps
        bool overlap = false;
        for(auto& lep : Leptons)
        {
          if (RazorAnalyzerLLP::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
        if(overlap) continue;
        leptons tmpElectron;
        tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);

        Leptons.push_back(tmpElectron);

        llp_tree->ele_passCutBasedIDVeto[llp_tree->nElectrons] = ele_passCutBasedIDVeto[i];

        llp_tree->elePt[llp_tree->nElectrons] = elePt[i];
        llp_tree->eleEta[llp_tree->nElectrons] = eleEta[i];
        llp_tree->eleE[llp_tree->nElectrons] = eleE[i];
        llp_tree->elePhi[llp_tree->nElectrons] = elePhi[i];

        llp_tree->nElectrons++;
    }

    //-------------------------------
    //Taus
    //-------------------------------
    for( int i = 0; i < nTaus; i++ )
    {
       if(_debug_sync)
       {
    	   std::cout << "nTaus " << nTaus << std::endl;
      	   std::cout << "iTau " << i << ", Pt " << tauPt[i] << ", Eta " << tauEta[i]<< ", Phi " << tauPhi[i]<< ", E " << tauE[i] << std::endl;
      	   std::cout << "iTau " << i << ", tau_IsLoose[i] " << tau_IsLoose[i]<< std::endl;
       }

	if(tauPt[i] <= 18 ) continue;
        if(fabs(tauEta[i]) > 2.3) continue;

	if(!tau_IsLoose[i]) continue;

        llp_tree->tau_IsLoose[llp_tree->nTaus] = tau_IsLoose[i];

        llp_tree->tauPt[llp_tree->nTaus] = tauPt[i];
        llp_tree->tauEta[llp_tree->nTaus] = tauEta[i];
        llp_tree->tauE[llp_tree->nTaus] = tauE[i];
        llp_tree->tauPhi[llp_tree->nTaus] = tauPhi[i];

        llp_tree->nTaus++;
    }
    
    //-------------------------------
    //Photons
    //-------------------------------
    for( int i = 0; i < nPhotons; i++ )
    {
       if(_debug_sync)
       {
    	   std::cout << "nPhotons " << nPhotons << std::endl;
      	   std::cout << "iPhoton " << i << ", Pt " << phoPt[i] << ", Eta " << phoEta[i]<< ", Phi " << phoPhi[i]<< ", E " << phoE[i] << std::endl;
      	   std::cout << "iPhoton " << i << ", pho_passCutBasedIDLoose[i] " << pho_passCutBasedIDLoose[i]<< std::endl;
       }

	if(phoPt[i] <= 15 ) continue;
        if(fabs(phoEta[i]) > 2.5) continue;

	if(!pho_passCutBasedIDLoose[i]) continue;

        llp_tree->pho_passCutBasedIDLoose[llp_tree->nPhotons] = pho_passCutBasedIDLoose[i];

        llp_tree->phoPt[llp_tree->nPhotons] = phoPt[i];
        llp_tree->phoEta[llp_tree->nPhotons] = phoEta[i];
        llp_tree->phoE[llp_tree->nPhotons] = phoE[i];
        llp_tree->phoPhi[llp_tree->nPhotons] = phoPhi[i];

        llp_tree->nPhotons++;
    }

    //-------------------------------
    //Leptons
    //-------------------------------
    sort(Leptons.begin(), Leptons.end(), my_largest_pt_lep);
    for ( auto &tmp : Leptons )
    {
      llp_tree->lepE[llp_tree->nLeptons]      = tmp.lepton.E();
      llp_tree->lepPt[llp_tree->nLeptons]     = tmp.lepton.Pt();
      llp_tree->lepEta[llp_tree->nLeptons]    = tmp.lepton.Eta();
      llp_tree->lepPhi[llp_tree->nLeptons]    = tmp.lepton.Phi();
      llp_tree->lepPdgId[llp_tree->nLeptons]  = tmp.pdgId;
      //llp_tree->lepDZ[llp_tree->nLeptons]     = tmp.dZ;
      //llp_tree->lepPassId[llp_tree->nLeptons] = tmp.passId;
      if(_debug) std::cout << "lepE " << tmp.lepton.E() << std::endl;


      // std::cout << "lepton pdg " << llp_tree->lepPdgId[llp_tree->nLeptons] << std::endl;
      llp_tree->nLeptons++;
    }
    if (triggered) trig_lepId->Fill(1);


  //-----------------------------------------------
  //Select Jets
  //-----------------------------------------------
  //std::vector<double> jetPtVector;
  //std::vector<double> jetCISVVector;
  std::vector<jets> Jets;
  //auto highest = [](auto a, auto b) { return a > b; };
  //cout <<"nJets :" << nJets << std::endl;

  if(_debug_jet) std::cout << "nJets " << nJets << std::endl;
  if(_debug_jet) std::cout << "jetE 0 " << jetE[0] << std::endl;
  //if(_debug_jet) std::cout << "jetChargedEMEnergyFraction 0 " << jetChargedEMEnergyFraction[0] << std::endl;
  //if(_debug_jet) std::cout << "jetNeutralEMEnergyFraction 0 " << jetNeutralEMEnergyFraction[0] << std::endl;
  if(_debug_jet) std::cout << "jetGammaMax_ET 0 " << jetGammaMax_ET[0] << std::endl;

  float ht = 0.;

  for(int i = 0; i < nJets; i++)
  {

    ht += jetPt[i];

    //------------------------------------------------------------
    //exclude selected muons and electrons from the jet collection
    //------------------------------------------------------------
    double deltaR = -1;
    for(auto& lep : Leptons){
      double thisDR = RazorAnalyzerLLP::deltaR(jetEta[i],jetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
      if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
    }
    if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

    //------------------------------------------------------------
    //Apply Jet Energy and Resolution Corrections
    //------------------------------------------------------------
    /*
    double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
       fixedGridRhoFastjetAll, jetJetArea[i] , JetCorrector);
    //cout <<"JEC :" << JEC << std::endl;

      TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );
*/
      TLorentzVector thisJet = makeTLorentzVector( jetPt[i], jetEta[i], jetPhi[i], jetE[i] );
      if( thisJet.Pt() < 30 ) continue;//According to the April 1st 2015 AN
      //if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
      if( fabs( thisJet.Eta() ) >= 2.4 ) continue;
      //if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
      // if ( !jetPassIDLoose[i] ) continue;
      // if (!(jetRechitE[i] > 0.0)) continue;
      // if(jetNRechits[i]<10) continue;
      // if(jetRechitT[i] < 0.0) continue;
      // if ((jetChargedHadronEnergyFraction[i]+jetChargedEMEnergyFraction[i]) > 0.4) continue;
      // if ((jetChargedHadronEnergyFraction[i]+jetNeutralHadronEnergyFraction[i])/(jetChargedEMEnergyFraction[i]+jetNeutralEMEnergyFraction[i]) < 0.2) continue;

      // std::cout <<jetRechitT[i] << "," << jetRechitE[i] <<  "," << jetNRechits[i] << std::endl;

      jets tmpJet;
      tmpJet.jet    = thisJet;
      tmpJet.time   = jetRechitT[i];
      tmpJet.ecalNRechits = jetNRechits[i];
      tmpJet.ecalRechitE = jetRechitE[i];
      tmpJet.jetNeutralHadronMultiplicity = jetNeutralHadronMultiplicity[i];
      tmpJet.jetChargedHadronMultiplicity = jetChargedHadronMultiplicity[i];
      tmpJet.jetMuonMultiplicity = jetMuonMultiplicity[i];
      tmpJet.jetElectronMultiplicity = jetElectronMultiplicity[i];
      tmpJet.jetPhotonMultiplicity = jetPhotonMultiplicity[i];
      tmpJet.jetNeutralHadronEnergyFraction = jetNeutralHadronEnergyFraction[i];
      tmpJet.jetChargedHadronEnergyFraction = jetChargedHadronEnergyFraction[i];
      tmpJet.jetMuonEnergyFraction = jetMuonEnergyFraction[i];
      tmpJet.jetElectronEnergyFraction = jetElectronEnergyFraction[i];
      tmpJet.jetPhotonEnergyFraction = jetPhotonEnergyFraction[i];
      tmpJet.jetCSV = jetCISV[i];

      tmpJet.jetGammaMax_ET = jetGammaMax_ET[i];
      tmpJet.jetMinDeltaRPVTracks = jetMinDeltaRPVTracks[i];

      tmpJet.jetChargedEMEnergyFraction = jetChargedEMEnergyFraction[i];
      tmpJet.jetNeutralEMEnergyFraction = jetNeutralEMEnergyFraction[i];
      Jets.push_back(tmpJet);

    }
    llp_tree->HT = ht;
    //-----------------------------
    //Require at least 2 jets
    //-----------------------------
/*
    if(pf)
    {
      if( Jets.size() < 1 ) continue;

    }
    else
    {
      if( caloJets.size() < 1 ) continue;

    }
*/
    if (triggered) trig_lepId_dijet->Fill(1);
    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

    if (Jets.size()>0)
    {
      llp_tree->jetMet_dPhi = RazorAnalyzerLLP::deltaPhi(jetPhi[0],metType1Phi);
      //TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM( metType1Pt, 0, metType1Phi, 0 );
      TLorentzVector jet0 = makeTLorentzVectorPtEtaPhiM( jetPt[0], 0, jetPhi[0], 0 );
      llp_tree->jetMet_dPhiStar = RazorAnalyzerLLP::deltaPhi(jetPhi[0],  (t1PFMET+jet0).Phi() );
    }
    else{
      llp_tree->jetMet_dPhi = -999.;
      llp_tree->jetMet_dPhiStar = -999.;
    }
    float jetMet_dPhiMin_temp = 999 ; 
    float jetMet_dPhiStarMin_temp = 999 ; 
    float jetMet_dPhiMin4_temp = 999 ; 

    for ( auto &tmp : Jets )
    {
      llp_tree->jetNeutralHadronMultiplicity[llp_tree->nJets] = tmp.jetNeutralHadronMultiplicity;
      llp_tree->jetChargedHadronMultiplicity[llp_tree->nJets] = tmp.jetChargedHadronMultiplicity;
      llp_tree->jetMuonMultiplicity[llp_tree->nJets] = tmp.jetMuonMultiplicity;
      llp_tree->jetElectronMultiplicity[llp_tree->nJets] = tmp.jetElectronMultiplicity;
      llp_tree->jetPhotonMultiplicity[llp_tree->nJets] = tmp.jetPhotonMultiplicity;
      llp_tree->jetNeutralHadronEnergyFraction[llp_tree->nJets] = tmp.jetNeutralHadronEnergyFraction;
      llp_tree->jetChargedHadronEnergyFraction[llp_tree->nJets] = tmp.jetChargedHadronEnergyFraction;
      llp_tree->jetMuonEnergyFraction[llp_tree->nJets] = tmp.jetMuonEnergyFraction;
      llp_tree->jetElectronEnergyFraction[llp_tree->nJets] = tmp.jetElectronEnergyFraction;
      llp_tree->jetPhotonEnergyFraction[llp_tree->nJets] = tmp.jetPhotonEnergyFraction;
      llp_tree->jetCSV[llp_tree->nJets] = tmp.jetCSV;

      llp_tree->jetPt[llp_tree->nJets] = tmp.jet.Pt();
      llp_tree->jetEta[llp_tree->nJets] = tmp.jet.Eta();
      llp_tree->jetE[llp_tree->nJets] = tmp.jet.E();
      llp_tree->jetPhi[llp_tree->nJets] = tmp.jet.Phi();
      llp_tree->jetTime[llp_tree->nJets] = tmp.time;
      llp_tree->ecalNRechits[llp_tree->nJets] = tmp.ecalNRechits;
      llp_tree->ecalRechitE[llp_tree->nJets] = tmp.ecalRechitE;

      llp_tree->jetGammaMax_ET[llp_tree->nJets] = tmp.jetGammaMax_ET;
      llp_tree->jetMinDeltaRPVTracks[llp_tree->nJets] = tmp.jetMinDeltaRPVTracks;

      llp_tree->jetChargedEMEnergyFraction[llp_tree->nJets] = tmp.jetChargedEMEnergyFraction;
      llp_tree->jetNeutralEMEnergyFraction[llp_tree->nJets] = tmp.jetNeutralEMEnergyFraction;

      //std::cout << "jetEta " << tmp.jet.Eta() << std::endl;
      //std::cout << "jetEta " << llp_tree->jetEta[llp_tree->nJets] << std::endl;
      
      if(jetMet_dPhiMin4_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi)) && llp_tree->nJets < 4)
      {
        jetMet_dPhiMin4_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi));
      }
      if (jetMet_dPhiMin_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi)))
      {
        jetMet_dPhiMin_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi));
      }     
      TLorentzVector jet_temp = makeTLorentzVectorPtEtaPhiM( tmp.jet.Pt(), 0, tmp.jet.Phi(), 0 );
      if (jetMet_dPhiStarMin_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(), (t1PFMET+jet_temp).Phi() )))
      {
        jetMet_dPhiStarMin_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(), (t1PFMET+jet_temp).Phi() ));
      }     
 
      llp_tree->nJets++;
    }
    if(_debug) std::cout << "nJets in tree " << llp_tree->nJets << std::endl;
    llp_tree->jetMet_dPhiMin = jetMet_dPhiMin_temp;
    llp_tree->jetMet_dPhiStarMin = jetMet_dPhiStarMin_temp;
    llp_tree->jetMet_dPhiMin4 = jetMet_dPhiMin4_temp;


    
    //gLLP
    for(int i = 0; i <N_MAX_LLPS; i++){
      llp_tree->gLLP_travel_time[i] = gLLP_travel_time[i];
      llp_tree->gLLP_e[i] = gLLP_e[i];
      llp_tree->gLLP_pt[i] = gLLP_pt[i];
      llp_tree->gLLP_eta[i] = gLLP_eta[i];
      llp_tree->gLLP_beta[i] = gLLP_beta[i];
      llp_tree->gLLP_phi[i] = gLLP_phi[i];
      llp_tree->gLLP_decay_vertex_x[i] = gLLP_decay_vertex_x[i];
      llp_tree->gLLP_decay_vertex_y[i] = gLLP_decay_vertex_y[i];
      llp_tree->gLLP_decay_vertex_z[i] = gLLP_decay_vertex_z[i];
      llp_tree->gLLP_prod_vertex_x[i] = gLLP_prod_vertex_x[i];
      llp_tree->gLLP_prod_vertex_y[i] = gLLP_prod_vertex_y[i];
      llp_tree->gLLP_prod_vertex_z[i] = gLLP_prod_vertex_z[i];

      //acceptance
      double decay_radius = sqrt( pow(gLLP_decay_vertex_x[i],2) + pow(gLLP_decay_vertex_y[i],2) );
      llp_tree->gLLP_decay_vertex_r[i] = decay_radius;
      double decay_z = abs( gLLP_decay_vertex_z[i] );
      double decay_eta = abs( gLLP_eta[i] );
      //if(decay_eta != 666) std::cout << "deb gLLP decay radius " << decay_radius << " z " << decay_z  << " eta " << decay_eta << endl;

      // barrel, 30 cm < r < 1.84 m, |z| < 3.76m  ( ecal outter edge )
      // EB_z = 268.36447217; // cm 129*sinh(1.479)
      bool inEB = false;
      if(decay_radius > 30. && decay_radius < 184. && decay_z < 376. && decay_eta != 666)
      {
	if(i==0)
	{
	  inEB = true;
          llp_tree->gLLP0_EB = inEB;
	}
	else
	{
          llp_tree->gLLP1_EB = true;
	}
      }
      //if(decay_eta != 666 && i==0) std::cout << "deb gLLP EB: " << jentry << " evtNum " << llp_tree->evtNum  << " index of llp " << i << " ; 0  " << llp_tree->gLLP0_EB  << " 1 " << llp_tree->gLLP1_EB << endl;

      // endcap, 1.5 < |eta| < 2.6, 1 m < |z| <  3.9 m
      bool inEE = false;
      if(decay_eta > 1.5 && decay_eta < 2.6 && decay_z > 100. && decay_z < 390. && decay_eta != 666)
      {
	if(i==0)
	{
          inEE = true;
          llp_tree->gLLP0_EE = inEE;
	}
	else
	{
          llp_tree->gLLP1_EE = true;
	}
      }


    }

    llp_tree->gLLP_dr = deltaR(gLLP_eta[0], gLLP_phi[0], gLLP_eta[1], gLLP_phi[1]);

    //gLLP daughters
    for(int i = 0; i <N_MAX_LLP_DAUGHTERS; i++){
/*
      llp_tree->gen_time[i] = gen_time[i];
      llp_tree->photon_travel_time[i] = photon_travel_time[i];
      llp_tree->gLLP_daughter_travel_time[i] = gLLP_daughter_travel_time[i];
      llp_tree->gLLP_daughter_e[i] = gLLP_daughter_e[i];
      llp_tree->gLLP_daughter_pt[i] = gLLP_daughter_pt[i];
      llp_tree->gLLP_daughter_eta[i] = gLLP_daughter_eta[i];
      llp_tree->gLLP_daughter_phi[i] = gLLP_daughter_phi[i];
      llp_tree->gLLP_daughter_eta_ecalcorr[i] = gLLP_daughter_eta_ecalcorr[i];
      llp_tree->gLLP_min_delta_r_match_jet[i] = gLLP_min_delta_r_match_jet[i];
      llp_tree->gLLP_daughter_match_jet_index[i] = gLLP_daughter_match_jet_index[i];
*/
      llp_tree->gLLP_daughter_EB[i] = gLLP_daughter_EB[i];
      llp_tree->gLLP_daughter_ETL[i] = gLLP_daughter_ETL[i];

      llp_tree->gLLP_daughter_photon_travel_time_EB[i] = gLLP_daughter_photon_travel_time_EB[i];
      llp_tree->gLLP_daughter_photon_travel_time_ETL[i] = gLLP_daughter_photon_travel_time_ETL[i];

      llp_tree->gLLP_daughter_travel_time_EB[i] = gLLP_daughter_travel_time_EB[i];
      llp_tree->gLLP_daughter_travel_time_ETL[i] = gLLP_daughter_travel_time_ETL[i];

      llp_tree->gen_time_daughter_EB[i] = gen_time_daughter_EB[i];
      llp_tree->gen_time_daughter_ETL[i] = gen_time_daughter_ETL[i];

      llp_tree->gLLP_daughter_match_jet_index[i] = gLLP_daughter_match_jet_index[i];
      llp_tree->gLLP_daughter_match_calojet_index[i] = gLLP_daughter_match_calojet_index[i];

      llp_tree->gLLP_daughter_min_delta_r_match_jet[i] = gLLP_daughter_min_delta_r_match_jet[i];
      llp_tree->gLLP_daughter_min_delta_r_match_calojet[i] = gLLP_daughter_min_delta_r_match_calojet[i];

      llp_tree->gLLP_daughter_id[i] = gLLP_daughter_id[i];
      llp_tree->gLLP_daughter_mass[i] = gLLP_daughter_mass[i];
      llp_tree->gLLP_daughter_e[i] = gLLP_daughter_e[i];
      llp_tree->gLLP_daughter_pt[i] = gLLP_daughter_pt[i];
      llp_tree->gLLP_daughter_eta[i] = gLLP_daughter_eta[i];
      llp_tree->gLLP_daughter_phi[i] = gLLP_daughter_phi[i];
      llp_tree->gLLP_daughter_eta_ecalcorr[i] = gLLP_daughter_eta_ecalcorr[i];
      llp_tree->gLLP_daughter_phi_ecalcorr[i] = gLLP_daughter_phi_ecalcorr[i];
    }

    //gLLP grandaughters
    for(int i = 0; i <N_MAX_LLP_GRAND_DAUGHTERS; i++){
      llp_tree->gLLP_grandaughter_EB[i] = gLLP_grandaughter_EB[i];
      llp_tree->gLLP_grandaughter_ETL[i] = gLLP_grandaughter_ETL[i];

      llp_tree->gLLP_grandaughter_photon_travel_time_EB[i] = gLLP_grandaughter_photon_travel_time_EB[i];
      llp_tree->gLLP_grandaughter_photon_travel_time_ETL[i] = gLLP_grandaughter_photon_travel_time_ETL[i];

      llp_tree->gLLP_grandaughter_travel_time_EB[i] = gLLP_grandaughter_travel_time_EB[i];
      llp_tree->gLLP_grandaughter_travel_time_ETL[i] = gLLP_grandaughter_travel_time_ETL[i];

      llp_tree->gen_time_grandaughter_EB[i] = gen_time_grandaughter_EB[i];
      llp_tree->gen_time_grandaughter_ETL[i] = gen_time_grandaughter_ETL[i];

      //llp_tree->gLLP_grandaughter_match_jet_index[i] = gLLP_grandaughter_match_jet_index[i];
      if(_debug_match) std::cout << "evt: "<<llp_tree->evtNum <<" n_jet " << llp_tree->nJets << std::endl; 
      double min_delta_r = 666.;
      int match_jet_index = -666;
      for ( int i_jet = 0; i_jet < llp_tree->nJets; i_jet++ )
      {
       double current_delta_r = deltaR(gLLP_grandaughter_eta_ecalcorr[i], gLLP_grandaughter_phi_ecalcorr[i], llp_tree->jetEta[i_jet], llp_tree->jetPhi[i_jet]);
      if(_debug_match) std::cout << " i_jet " << i_jet << ", match_jet_index " << match_jet_index << ", min_delta_r " << min_delta_r <<", current_delta_r "<< current_delta_r << std::endl;
       if ( current_delta_r < min_delta_r )
       {
         min_delta_r = current_delta_r;
         match_jet_index = i_jet;
       }
      if(_debug_match) std::cout << " i_jet " << i_jet << ", match_jet_index " << match_jet_index << ", min_delta_r " << min_delta_r << std::endl;
      }//end matching to jets 
      
      if ( min_delta_r < 0.4 )
      //if ( min_delta_r < 20 )
      {
      llp_tree->gLLP_grandaughter_match_jet_index[i] = match_jet_index;
      llp_tree->gLLP_grandaughter_min_delta_r_match_jet[i] = min_delta_r;
      //llp_tree->matched[match_jet_index] = true;
      if(i<N_MAX_LLP_GRAND_DAUGHTERS/2) llp_tree->jet_matched_gLLP0_grandaughter[match_jet_index] = true;
      else llp_tree->jet_matched_gLLP1_grandaughter[match_jet_index] = true;
      if(_debug_match) std::cout << " i " << i << ", match_jet_index " << match_jet_index << ", min_delta_r " << min_delta_r << std::endl;
      }
      
      //llp_tree->gLLP_grandaughter_min_delta_r_match_jet[i] = gLLP_grandaughter_min_delta_r_match_jet[i];

      if(_debug_calojet)
      {
      llp_tree->gLLP_grandaughter_match_calojet_index[i] = gLLP_grandaughter_match_calojet_index[i];
      llp_tree->gLLP_grandaughter_min_delta_r_match_calojet[i] = gLLP_grandaughter_min_delta_r_match_calojet[i];
      }

      llp_tree->gLLP_grandaughter_id[i] = gLLP_grandaughter_id[i];
      llp_tree->gLLP_grandaughter_mass[i] = gLLP_grandaughter_mass[i];
      llp_tree->gLLP_grandaughter_e[i] = gLLP_grandaughter_e[i];
      llp_tree->gLLP_grandaughter_pt[i] = gLLP_grandaughter_pt[i];
      llp_tree->gLLP_grandaughter_eta[i] = gLLP_grandaughter_eta[i];
      llp_tree->gLLP_grandaughter_phi[i] = gLLP_grandaughter_phi[i];
      llp_tree->gLLP_grandaughter_eta_ecalcorr[i] = gLLP_grandaughter_eta_ecalcorr[i];
      llp_tree->gLLP_grandaughter_phi_ecalcorr[i] = gLLP_grandaughter_phi_ecalcorr[i];
    }
    //dr
    llp_tree->gLLP_grandaughter_dr1 = deltaR(gLLP_grandaughter_eta_ecalcorr[0], gLLP_grandaughter_phi_ecalcorr[0], gLLP_grandaughter_eta_ecalcorr[1], gLLP_grandaughter_phi_ecalcorr[1]);
    llp_tree->gLLP_grandaughter_dr2 = deltaR(gLLP_grandaughter_eta_ecalcorr[2], gLLP_grandaughter_phi_ecalcorr[2], gLLP_grandaughter_eta_ecalcorr[3], gLLP_grandaughter_phi_ecalcorr[3]);
    if(_debug_match) std::cout << " dr1(b1, b2) " << llp_tree->gLLP_grandaughter_dr1 << ", dr2(b1, b2) " << llp_tree->gLLP_grandaughter_dr2 << std::endl;
    if(_debug_match) std::cout << " gLLP_grandaughter_match_jet_index " << gLLP_grandaughter_match_jet_index[0] << ", " << gLLP_grandaughter_match_jet_index[1] << ", " << gLLP_grandaughter_match_jet_index[2] << ", "  << gLLP_grandaughter_match_jet_index[3] << std::endl;

    if(llp_tree->gLLP_grandaughter_match_jet_index[0]!=-666 && llp_tree->gLLP_grandaughter_match_jet_index[1]!=-666)
    {
	llp_tree->gLLP_grandaughter_matched_jet_dr1 = deltaR(llp_tree->jetEta[llp_tree->gLLP_grandaughter_match_jet_index[0]], llp_tree->jetPhi[llp_tree->gLLP_grandaughter_match_jet_index[0]], llp_tree->jetEta[llp_tree->gLLP_grandaughter_match_jet_index[1]], llp_tree->jetPhi[llp_tree->gLLP_grandaughter_match_jet_index[1]] ); 
    }

    if(llp_tree->gLLP_grandaughter_match_jet_index[2]!=-666 && llp_tree->gLLP_grandaughter_match_jet_index[3]!=-666)
    {
	llp_tree->gLLP_grandaughter_matched_jet_dr2 = deltaR(llp_tree->jetEta[llp_tree->gLLP_grandaughter_match_jet_index[2]], llp_tree->jetPhi[llp_tree->gLLP_grandaughter_match_jet_index[2]], llp_tree->jetEta[llp_tree->gLLP_grandaughter_match_jet_index[3]], llp_tree->jetPhi[llp_tree->gLLP_grandaughter_match_jet_index[3]] ); 
    }
    if(_debug_match) std::cout << " dr1(b1j, b2j) " << llp_tree->gLLP_grandaughter_matched_jet_dr1 << ", dr2(b1j, b2j) " << llp_tree->gLLP_grandaughter_matched_jet_dr2 << std::endl;

    if(_debug_trg) std::cout << "end: 310 " << HLTDecision[310] << std::endl;
    if(_debug_trg) std::cout << "end: 310 " << llp_tree->HLTDecision[310] << std::endl;

    llp_tree->tree_->Fill();
  }

    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    llp_tree->tree_->Write();
    outFile->Write();
    outFile->Close();
}
