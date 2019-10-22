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
#define _debug_jet 0
#define _debug_calojet 0

#define N_MAX_LLP_DAUGHTERS 4
#define N_MAX_LLP_GRAND_DAUGHTERS 4
#define N_MAX_LLPS 2
#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define NTriggersMAX 601 //Number of trigger in the .dat file
using namespace std;

struct leptons
{
  TLorentzVector lepton;
  int pdgId;
  float dZ;
  // bool passLooseId;
  // bool passMediumId;
  // bool passTightId;
  bool passId;
};


struct jets
{
  TLorentzVector jet;
  float time;
  bool passId;
  // bool passLooseId;
  // bool passMediumId;
  // bool passTightId;
  bool isCSVL;
  //bool matched;
  bool jet_matched_gLLP_daughter;
  bool jet_matched_gLLP_grandaughter;
  int ecalNRechits;
  float ecalRechitE;
  float jetChargedEMEnergyFraction;
  float jetNeutralEMEnergyFraction;
  float jetChargedHadronEnergyFraction;
  float jetNeutralHadronEnergyFraction;
  float jetGammaMax_ET;
  float jetMinDeltaRPVTracks;
  float jetPtAllPVTracks;
  float energy_frac;
  float sig_et1;
  float sig_et2;

};

struct calojets
{
  TLorentzVector calojet;
  float time;
  bool passId;
  bool isCSVL;
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
void SusyLLP::Analyze(bool isData, int options, string outputfilename, string analysisTag)
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
    cout << "signal / bkg label: " << label << "\n";
  }
  else if ((options/10) % 10 == 2){
    label = "zH";
    cout << "signal / bkg label: " << label << "\n";
  }
  else if ((options/10) % 10 == 3){
    label = "bkg_wH";
    cout << "signal / bkg label: " << label << "\n";
  }
  else if ((options/10) % 10 == 4){
    label = "bkg_zH";
    cout << "signal / bkg label: " << label << "\n";
  }
  else{
    cout << "What signal / bkg it is? Label not defined. \n";
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
  if (analysisTag == "Razor2015_76X") helper = new RazorHelper("Razor2015_76X", isData, false);
  else if (analysisTag == "Razor2016_MoriondRereco") helper = new RazorHelper("Razor2016_MoriondRereco", isData, false);
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
    if (label =="bkg_wH"|| label == "bkg_zH"){
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
    if (label == "zH" || label == "wH"){
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

    for (int i=0; i < nBunchXing; ++i)
    {
      if (BunchXing[i] == 0)
      {
        llp_tree->npu = nPUmean[i];
      }
    }
    if(_debug && llp_tree->npu != 0 ) std::cout << "npu " << llp_tree->npu << std::endl;
    //get NPU
    llp_tree->npv = nPV;
    llp_tree->rho = fixedGridRhoFastjetAll;
    llp_tree->met = metType1Pt;
    llp_tree->metPhi = metType1Phi;
    if(_debug) std::cout << "npv " << llp_tree->npv << std::endl;

    //Triggers
    for(int i = 0; i < NTriggersMAX; i++){
      llp_tree->HLTDecision[i] = HLTDecision[i];
    }
    bool triggered = false;
    for(int i = 0; i < NTrigger; i++)
    {
      int trigger_temp = trigger_paths[i];

      triggered = triggered || HLTDecision[trigger_temp];

    }
    if (triggered) trig->Fill(1);
    if(_debug) std::cout << "triggered " << triggered << std::endl;
    //*************************************************************************
    //Start Object Selection
    //*************************************************************************
    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
    if(_debug) std::cout << "nMuons " << nMuons << std::endl;
    for( int i = 0; i < nMuons; i++ )
    {
      if(!isMuonPOGLooseMuon(i)) continue;
      if(muonPt[i] < muonPt_cut) continue;
      if(fabs(muonEta[i]) > 2.4) continue;

      //remove overlaps
      bool overlap = false;
      for(auto& lep : Leptons)
      {
        if (RazorAnalyzer::deltaR(muonEta[i],muonPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
      }
      if(overlap) continue;

      leptons tmpMuon;
      tmpMuon.lepton.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i], MU_MASS);
      tmpMuon.pdgId = 13 * -1 * muonCharge[i];
      tmpMuon.dZ = muon_dZ[i];
      tmpMuon.passId = isMuonPOGTightMuon(i);

      Leptons.push_back(tmpMuon);
    }
    //std::cout << "deb6 " << jentry << std::endl;
    //-------------------------------
    //Electrons
    //-------------------------------
    if(_debug) std::cout << "nElectrons " << nElectrons << std::endl;
    for( int i = 0; i < nElectrons; i++ )
    {


      if (!isEGammaPOGLooseElectron(i, true, true, true, "Summer16")) continue;



      if(elePt[i] < elePt_cut) continue;

      if(fabs(eleEta[i]) > 2.5) continue;

      //remove overlaps
      bool overlap = false;
      for(auto& lep : Leptons)
      {
        if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
      }
      if(overlap) continue;
      // std::cout << "here" << std::endl;
      leptons tmpElectron;
      tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);
      tmpElectron.pdgId = 11 * -1 * eleCharge[i];
      tmpElectron.dZ = ele_dZ[i];
      // tmpElectron.passId = passMVALooseElectronID(i) && passEGammaPOGLooseElectronIso(i);
      Leptons.push_back(tmpElectron);
    }

    sort(Leptons.begin(), Leptons.end(), my_largest_pt_lep);
    //std::cout << "deb7 " << jentry << std::endl;
    for ( auto &tmp : Leptons )
    {
      llp_tree->lepE[llp_tree->nLeptons]      = tmp.lepton.E();
      llp_tree->lepPt[llp_tree->nLeptons]     = tmp.lepton.Pt();
      llp_tree->lepEta[llp_tree->nLeptons]    = tmp.lepton.Eta();
      llp_tree->lepPhi[llp_tree->nLeptons]    = tmp.lepton.Phi();
      llp_tree->lepPdgId[llp_tree->nLeptons]  = tmp.pdgId;
      llp_tree->lepDZ[llp_tree->nLeptons]     = tmp.dZ;
      llp_tree->lepPassId[llp_tree->nLeptons] = tmp.passId;
      if(_debug) std::cout << "lepE " << tmp.lepton.E() << std::endl;


      // std::cout << "lepton pdg " << llp_tree->lepPdgId[llp_tree->nLeptons] << std::endl;
      llp_tree->nLeptons++;
    }

    //----------------
    //Find Z Candidate
    //----------------


    double ZMass = -999;
    double ZPt = -999;
    double tmpDistToZPole = 9999;
    pair<uint,uint> ZCandidateLeptonIndex;
    bool foundZ = false;
    TLorentzVector ZCandidate;
    if(_debug) std::cout << "Leptons.size() " << Leptons.size() << std::endl;
    for( uint i = 0; i < Leptons.size(); i++ )
    {
      for( uint j = 0; j < Leptons.size(); j++ )
      {
        if (!( Leptons[i].pdgId == -1*Leptons[j].pdgId )) continue;// same flavor opposite charge
        double tmpMass = (Leptons[i].lepton+Leptons[j].lepton).M();

        //select the pair closest to Z pole mass
        if ( fabs( tmpMass - Z_MASS) < tmpDistToZPole)
        {
          tmpDistToZPole = tmpMass;
          if (Leptons[i].pdgId > 0)
          {
            ZCandidateLeptonIndex = pair<int,int>(i,j);
          }
          else
          {
            ZCandidateLeptonIndex = pair<int,int>(j,i);
          }
          ZMass = tmpMass;
          ZPt = (Leptons[i].lepton+Leptons[j].lepton).Pt();
          ZCandidate = Leptons[i].lepton+Leptons[j].lepton;
          foundZ = true;
        }
      }
    }

    if (foundZ && fabs(ZMass-Z_MASS) < 30.0)
    {
      llp_tree->ZMass = ZMass;
      llp_tree->ZPt   = ZPt;
      llp_tree->ZEta  = ZCandidate.Eta();
      llp_tree->ZPhi  = ZCandidate.Phi();
      llp_tree->ZleptonIndex1 = ZCandidateLeptonIndex.first;
      llp_tree->ZleptonIndex2 = ZCandidateLeptonIndex.second;

      //match to gen leptons
      //if (abs(lep1Id) == 11) lep1IsPrompt = matchesGenElectron(lep1Eta,lep1Phi);
      //else lep1IsPrompt = matchesGenMuon(lep1Eta,lep1Phi);
      //if (abs(lep2Id) == 11) lep2IsPrompt = matchesGenElectron(lep2Eta,lep2Phi);
      //else lep2IsPrompt = matchesGenMuon(lep2Eta,lep2Phi);
    } // endif foundZ
    //------------------------
    //require 1 lepton
    //------------------------
    // if (nMuons == 0 && !(nElectrons == 0)){
    //   std::cout <<nMuons << "," << nElectrons <<  "," << llp_tree->nLeptons <<  "," << llp_tree->met << std::endl;
    // }

    //if ( Leptons.size() < nLepton_cut ) continue;
    //TLorentzVector met;
    //TLorentzVector visible = Leptons[0].lepton;
    //met.SetPtEtaPhiE(metType1Pt,0,metType1Phi,metType1Pt);
    //llp_tree->MT = GetMT(visible,met);

    // else{
    //   if ( Leptons.size() < 2 ) continue;
    //   if (!(foundZ && fabs(ZMass-Z_MASS) < 15.0 )) continue;
    // }
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
  if(_debug_jet) std::cout << "jetChargedEMEnergyFraction 0 " << jetChargedEMEnergyFraction[0] << std::endl;
  if(_debug_jet) std::cout << "jetNeutralEMEnergyFraction 0 " << jetNeutralEMEnergyFraction[0] << std::endl;
  if(_debug_jet) std::cout << "jetGammaMax_ET 0 " << jetGammaMax_ET[0] << std::endl;
  for(int i = 0; i < nJets; i++)
  {

    //------------------------------------------------------------
    //exclude selected muons and electrons from the jet collection
    //------------------------------------------------------------
    double deltaR = -1;
    for(auto& lep : Leptons){
      double thisDR = RazorAnalyzer::deltaR(jetEta[i],jetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
      if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
    }
    if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

    //------------------------------------------------------------
    //Apply Jet Energy and Resolution Corrections
    //------------------------------------------------------------
    double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
       fixedGridRhoFastjetAll, jetJetArea[i] , JetCorrector);
    //cout <<"JEC :" << JEC << std::endl;

      TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );

      if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
      if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
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
      tmpJet.passId = jetPassIDLoose[i];
      //tmpJet.matched = jet_matched[i];
      tmpJet.jet_matched_gLLP_daughter = jet_matched_gLLP_daughter[i];
      tmpJet.jet_matched_gLLP_grandaughter = jet_matched_gLLP_grandaughter[i];
      tmpJet.energy_frac = jet_energy_frac[i];
      tmpJet.sig_et1 = jet_sig_et1[i];
      tmpJet.sig_et2 = jet_sig_et2[i];
      // std::cout<<tmpJet.sig_et1<<","<<jet_sig_et1[i]<<std::endl;
      tmpJet.isCSVL = isCSVL(i);
      //if (isCSVL(i)) NBJet20++;
      //if (isCSVL(i) && thisJet.Pt() > 30) NBJet30++;
      tmpJet.ecalNRechits = jetNRechits[i];
      tmpJet.ecalRechitE = jetRechitE[i];
      tmpJet.jetChargedEMEnergyFraction = jetChargedEMEnergyFraction[i];
      tmpJet.jetNeutralEMEnergyFraction = jetNeutralEMEnergyFraction[i];
      tmpJet.jetChargedHadronEnergyFraction = jetChargedHadronEnergyFraction[i];
      tmpJet.jetNeutralHadronEnergyFraction = jetNeutralHadronEnergyFraction[i];
      tmpJet.jetGammaMax_ET = jetGammaMax_ET[i];
      tmpJet.jetMinDeltaRPVTracks = jetMinDeltaRPVTracks[i];
      tmpJet.jetPtAllPVTracks = jetPtAllPVTracks[i];

      if(_debug_jet) std::cout << "jetE " << jetE[i] << std::endl;
      if(_debug_jet) std::cout << "jetChargedEMEnergyFraction " << jetChargedEMEnergyFraction[i] << std::endl;
      if(_debug_jet) std::cout << "jetNeutralEMEnergyFraction " << jetNeutralEMEnergyFraction[i] << std::endl;
      if(_debug_jet) std::cout << "jetGammaMax_ET " << jetGammaMax_ET[i] << std::endl;
      if(_debug_jet) std::cout << "jetMinDeltaRPVTracks " << jetMinDeltaRPVTracks[i] << std::endl;
      if(_debug_jet) std::cout << "jetPtAllPVTracks " << jetPtAllPVTracks[i] << std::endl;

      Jets.push_back(tmpJet);

    }
    std::vector<calojets> caloJets;
    //auto highest = [](auto a, auto b) { return a > b; };

    if(_debug) std::cout << "nCaloJets " << nCaloJets << std::endl;
    for(int i = 0; i < nCaloJets; i++)
    {

      //------------------------------------------------------------
      //exclude selected muons and electrons from the jet collection
      //------------------------------------------------------------
      double deltaR = -1;
      for(auto& lep : Leptons){
        double thisDR = RazorAnalyzer::deltaR(calojetEta[i],calojetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
        if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      }
      if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

      //------------------------------------------------------------
      //Apply Jet Energy and Resolution Corrections
      //------------------------------------------------------------
      double JEC = JetEnergyCorrectionFactor(calojetPt[i], calojetEta[i], calojetPhi[i], calojetE[i],
         fixedGridRhoFastjetAll, calojetJetArea[i] , JetCorrector);

        TLorentzVector thisJet = makeTLorentzVector( calojetPt[i]*JEC,calojetEta[i], calojetPhi[i], calojetE[i]*JEC );

        if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
        if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
        // if ( !jetPassIDLoose[i] ) continue;
        // if (!(jetRechitE[i] > 0.0)) continue;
        // if(jetNRechits[i]<10) continue;
        // if(jetRechitT[i] < 0.0) continue;
        // if ((jetChargedHadronEnergyFraction[i]+jetChargedEMEnergyFraction[i]) > 0.4) continue;
        // if ((jetChargedHadronEnergyFraction[i]+jetNeutralHadronEnergyFraction[i])/(jetChargedEMEnergyFraction[i]+jetNeutralEMEnergyFraction[i]) < 0.2) continue;

        // std::cout <<jetRechitT[i] << "," << jetRechitE[i] <<  "," << jetNRechits[i] << std::endl;


        calojets tmpJet;
        tmpJet.calojet    = thisJet;
        tmpJet.time   = calojetRechitT[i];
        tmpJet.passId = calojetPassIDLoose[i];
        tmpJet.isCSVL = isCSVL(i);
        //if (isCSVL(i)) NBJet20++;
        //if (isCSVL(i) && thisJet.Pt() > 30) NBJet30++;
        tmpJet.ecalNRechits = calojetNRechits[i];
        tmpJet.ecalRechitE = calojetRechitE[i];

        tmpJet.calojet_EMEnergyFraction = calojet_EMEnergyFraction[i];
        tmpJet.calojet_HadronicEnergyFraction = calojet_HadronicEnergyFraction[i];
        //tmpJet.jetChargedEMEnergyFraction = calojetChargedEMEnergyFraction[i];
        //tmpJet.jetNeutralEMEnergyFraction = calojetNeutralEMEnergyFraction[i];
        //tmpJet.jetChargedHadronEnergyFraction = calojetChargedHadronEnergyFraction[i];
        //tmpJet.jetNeutralHadronEnergyFraction = calojetNeutralHadronEnergyFraction[i];
        tmpJet.calojetGammaMax_ET = calojetGammaMax_ET[i];
        tmpJet.calojetMinDeltaRPVTracks = calojetMinDeltaRPVTracks[i];
        tmpJet.calojetPtAllPVTracks = calojetPtAllPVTracks[i];

        if(_debug_calojet) std::cout << "calojetE " << calojetE[i] << std::endl;
        if(_debug_calojet) std::cout << "calojetChargedEMEnergyFraction " << calojetChargedEMEnergyFraction[i] << std::endl;
        if(_debug_calojet) std::cout << "calojetNeutralEMEnergyFraction " << calojetNeutralEMEnergyFraction[i] << std::endl;
        if(_debug_calojet) std::cout << "calojetGammaMax_ET " << calojetGammaMax_ET[i] << std::endl;
        if(_debug_calojet) std::cout << "calojetMinDeltaRPVTracks " << calojetMinDeltaRPVTracks[i] << std::endl;
        if(_debug_calojet) std::cout << "calojetPtAllPVTracks " << calojetPtAllPVTracks[i] << std::endl;

        caloJets.push_back(tmpJet);

      }
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

    for ( auto &tmp : Jets )
    {
      llp_tree->jetE[llp_tree->nJets] = tmp.jet.E();
      llp_tree->jetEt[llp_tree->nJets] = tmp.jet.Et();
      llp_tree->jetPt[llp_tree->nJets] = tmp.jet.Pt();
      llp_tree->jetEta[llp_tree->nJets] = tmp.jet.Eta();
      llp_tree->jetPhi[llp_tree->nJets] = tmp.jet.Phi();
      llp_tree->jetTime[llp_tree->nJets] = tmp.time;
      llp_tree->jetPassId[llp_tree->nJets] = tmp.passId;
      //llp_tree->matched[llp_tree->nJets] = tmp.matched;
      llp_tree->jet_matched_gLLP_daughter[llp_tree->nJets] = tmp.jet_matched_gLLP_daughter;
      llp_tree->jet_matched_gLLP_grandaughter[llp_tree->nJets] = tmp.jet_matched_gLLP_grandaughter;
      llp_tree->jet_sig_et1[llp_tree->nJets] = tmp.sig_et1;
      llp_tree->jet_sig_et2[llp_tree->nJets] = tmp.sig_et2;
      llp_tree->jet_energy_frac[llp_tree->nJets] = tmp.energy_frac;
      llp_tree->ecalNRechits[llp_tree->nJets] = tmp.ecalNRechits;
      llp_tree->ecalRechitE[llp_tree->nJets] = tmp.ecalRechitE;
      llp_tree->jetChargedEMEnergyFraction[llp_tree->nJets] = tmp.jetChargedEMEnergyFraction;
      llp_tree->jetNeutralEMEnergyFraction[llp_tree->nJets] = tmp.jetNeutralEMEnergyFraction;
      llp_tree->jetChargedHadronEnergyFraction[llp_tree->nJets] = tmp.jetChargedHadronEnergyFraction;
      llp_tree->jetNeutralHadronEnergyFraction[llp_tree->nJets] = tmp.jetNeutralHadronEnergyFraction;
      llp_tree->jetGammaMax_ET[llp_tree->nJets] = tmp.jetGammaMax_ET;
      llp_tree->jetMinDeltaRPVTracks[llp_tree->nJets] = tmp.jetMinDeltaRPVTracks;
      llp_tree->jetPtAllPVTracks[llp_tree->nJets] = tmp.jetPtAllPVTracks;

      if(_debug_jet) std::cout << "jetE " << tmp.jet.E() << std::endl;
      if(_debug_jet) std::cout << "jetChargedEMEnergyFraction " << tmp.jetChargedEMEnergyFraction << std::endl;
      if(_debug_jet) std::cout << "jetNeutralEMEnergyFraction " << tmp.jetNeutralEMEnergyFraction << std::endl;
      if(_debug_jet) std::cout << "jetGammaMax_ET " << tmp.jetGammaMax_ET << std::endl;
      if(_debug_jet) std::cout << "jetMinDeltaRPVTracks " << tmp.jetMinDeltaRPVTracks << std::endl;
      if(_debug_jet) std::cout << "jetPtAllPVTracks " << tmp.jetPtAllPVTracks << std::endl;

      // std::cout <<tmp.time << "," <<tmp.ecalRechitE <<  "," << tmp.ecalNRechits << llp_tree->nJets<<std::endl;

      llp_tree->nJets++;
    }
    if(_debug) std::cout << "nJets in tree " << llp_tree->nJets << std::endl;
    sort(caloJets.begin(), caloJets.end(), my_largest_pt_calojet);

    for ( auto &tmp : caloJets )
    {
      llp_tree->calojetE[llp_tree->nCaloJets] = tmp.calojet.E();
      llp_tree->calojetEt[llp_tree->nCaloJets] = tmp.calojet.Et();
      llp_tree->calojetPt[llp_tree->nCaloJets] = tmp.calojet.Pt();
      llp_tree->calojetEta[llp_tree->nCaloJets] = tmp.calojet.Eta();
      llp_tree->calojetPhi[llp_tree->nCaloJets] = tmp.calojet.Phi();
      llp_tree->calojetTime[llp_tree->nCaloJets] = tmp.time;
      llp_tree->calojetPassId[llp_tree->nCaloJets] = tmp.passId;
      llp_tree->calojetNRechits[llp_tree->nCaloJets] = tmp.ecalNRechits;
      llp_tree->calojetRechitE[llp_tree->nCaloJets] = tmp.ecalRechitE;

      llp_tree->calojet_EMEnergyFraction[llp_tree->nCaloJets] = tmp.calojet_EMEnergyFraction;
      llp_tree->calojet_HadronicEnergyFraction[llp_tree->nCaloJets] = tmp.calojet_HadronicEnergyFraction;
      llp_tree->calojetGammaMax_ET[llp_tree->nCaloJets] = tmp.calojetGammaMax_ET;
      llp_tree->calojetMinDeltaRPVTracks[llp_tree->nCaloJets] = tmp.calojetMinDeltaRPVTracks;
      llp_tree->calojetPtAllPVTracks[llp_tree->nCaloJets] = tmp.calojetPtAllPVTracks;

      // std::cout <<tmp.time << "," <<tmp.ecalRechitE <<  "," << tmp.ecalNRechits << llp_tree->nJets<<std::endl;
      if(_debug_calojet) std::cout << "calojetE " << tmp.calojet.E() << std::endl;
      if(_debug_calojet) std::cout << "calojet_EMEnergyFraction " << tmp.calojet_EMEnergyFraction << std::endl;
      if(_debug_calojet) std::cout << "calojet_HadronicEnergyFraction " << tmp.calojet_HadronicEnergyFraction << std::endl;
      if(_debug_calojet) std::cout << "calojetGammaMax_ET " << tmp.calojetGammaMax_ET << std::endl;
      if(_debug_calojet) std::cout << "calojetMinDeltaRPVTracks " << tmp.calojetMinDeltaRPVTracks << std::endl;
      if(_debug_calojet) std::cout << "calojetPtAllPVTracks " << tmp.calojetPtAllPVTracks << std::endl;


      llp_tree->nCaloJets++;
    }
    if(_debug) std::cout << "nCaloJets in tree " << llp_tree->nCaloJets << std::endl;
    //std::cout << "deb fill: " << llp_tree->nLeptons << " " << jentry << endl;
    
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

    }

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

      llp_tree->gLLP_grandaughter_match_jet_index[i] = gLLP_grandaughter_match_jet_index[i];
      llp_tree->gLLP_grandaughter_match_calojet_index[i] = gLLP_grandaughter_match_calojet_index[i];

      llp_tree->gLLP_grandaughter_min_delta_r_match_jet[i] = gLLP_grandaughter_min_delta_r_match_jet[i];
      llp_tree->gLLP_grandaughter_min_delta_r_match_calojet[i] = gLLP_grandaughter_min_delta_r_match_calojet[i];

      llp_tree->gLLP_grandaughter_id[i] = gLLP_grandaughter_id[i];
      llp_tree->gLLP_grandaughter_mass[i] = gLLP_grandaughter_mass[i];
      llp_tree->gLLP_grandaughter_e[i] = gLLP_grandaughter_e[i];
      llp_tree->gLLP_grandaughter_pt[i] = gLLP_grandaughter_pt[i];
      llp_tree->gLLP_grandaughter_eta[i] = gLLP_grandaughter_eta[i];
      llp_tree->gLLP_grandaughter_phi[i] = gLLP_grandaughter_phi[i];
      llp_tree->gLLP_grandaughter_eta_ecalcorr[i] = gLLP_grandaughter_eta_ecalcorr[i];
      llp_tree->gLLP_grandaughter_phi_ecalcorr[i] = gLLP_grandaughter_phi_ecalcorr[i];
    }

    llp_tree->tree_->Fill();
  }

    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();
}
