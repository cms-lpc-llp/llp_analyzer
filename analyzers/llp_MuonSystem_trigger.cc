#include "llp_MuonSystem_trigger.h"
#include "RazorHelper.h"
#include "LiteLiteTreeMuonSystem.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
#include "DBSCAN.h"
//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"


using namespace std::chrono;
using namespace std;

struct greater_than_pt
{
  inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
};

struct leptons
{
  TLorentzVector lepton;
  int pdgId;
  float dZ;
  // bool passLooseId;
  // bool passMediumId;
  bool passId;
  bool passVetoId;
  bool passLooseIso;
  bool passTightIso;
  bool passVTightIso;
  bool passVVTightIso;

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
  int ecalNRechits;
  float ecalRechitE;
  float jetChargedEMEnergyFraction;
  float jetNeutralEMEnergyFraction;
  float jetChargedHadronEnergyFraction;
  float jetNeutralHadronEnergyFraction;
  bool jetPassMuFrac;
  float jetPtJESUp;
  float jetPtJESDown;
  float jetEJESUp;
  float jetEJESDown;
  float JecUnc;

};

//lepton highest pt comparator
struct largest_pt
{
  inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt;


struct largest_pt_jet
{
  inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;

void llp_MuonSystem_trigger::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  //options format: MH/MX/ctau/condor: 1000/300/0/1
  // mh can be 3-4 digits, mx is always 3 digits, ctau is one digit(number of zeros), last digit is condor option
  // mh can be 3-4 digits, mx is always 3 digits, ctau is 2 digit(number of zeros), last digit is condor option
  //
  //
  // int mx = int(options/1000)%1000;
  // int mh = options/1000000;
  // int ctau = pow(10, int(options/10)%10) * int(int(options/100)%10);
  //
  // cout<<"mh "<<mh<<", mx "<<mx<<", ctau "<<ctau<<endl;

  bool eMuFlag = int(options/10) == 1;
  int option = options%10;
  // if (options % 1){
  //   option = 1; // used when running condor
  // }
  // else{
  //   option = 0;// used when running locally
  // }

  if( isData )
  {
    std::cout << "[INFO]: running on data with option: " << option << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with option: " << option << std::endl;
  }
  if( eMuFlag )
  {
    std::cout << "[INFO]: running with muon" << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running with electrons " << option << std::endl;
  }



  const float ELE_MASS = 0.000511;
  const float MU_MASS  = 0.105658;
  const float Z_MASS   = 91.2;

  if (analysisTag == ""){
    analysisTag = "Razor2016_80X";

  }
  int wzId;


  const int zh_lepton0_cut = 15;
  const int zh_lepton1_cut = 15;

  const int wh_muonPt_cut = 25;
  const int wh_elePt_cut = 35;


  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile;
  outFile = new TFile(outfilename.c_str(), "RECREATE");


  LiteLiteTreeMuonSystem *MuonSystem = new LiteLiteTreeMuonSystem;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  // for signals, need one output file for each signal point
  map<pair<int,int>, TFile*> Files2D;
  map<pair<int,int>, TTree*> Trees2D;
  map<pair<int,int>, TH1F*> NEvents2D;
  map<pair<int,int>, TH1F*> accep2D;
  map<pair<int,int>, TH1F*> accep_met2D;
  map<pair<int,int>, TH1F*> Total2D;


  // map<pair<int,int>, TH1F*> smsSumScaleWeights2D;
  // map<pair<int,int>, TH1F*> smsSumPdfWeights2D;
  // map<pair<int,int>, TH1F*> smsNISRJets2D;
  // map<pair<int,int>, TH1F*> smsPtISR2D;
  // map<pair<int,int>, TH1F*> smsNPV2D;




  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  // TH1F *NEvents_genweight = new TH1F("NEvents_genweight", "NEvents_genweight", 1, 1, 2);

  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
  if(cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1



  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  helper = new RazorHelper(analysisTag, isData, false);



  std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
  std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();



  //*************************************************************************
  //Look over Input File Events
  //*************************************************************************
  if (fChain == 0) return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nbytes = 0, nb = 0;
  clock_t start, end;
  start = clock();
  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {

    //begin event

    if(jentry % 10000 == 0)
    {
      end = clock();
      double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
      cout << "Processing entry " << jentry << endl;
      cout << "Time taken by program is : " << time_taken << endl;
      start = clock();
    }
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //GetEntry(ientry);
    nb = fChain->GetEntry(jentry); nbytes += nb;

    //fill normalization histogram
    //std::cout << "deb0 " << jentry << std::endl;
    MuonSystem->InitVariables();
    //std::cout << "deb1 " << jentry << std::endl;



    //event info
    if (isData)
    {
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    }
    else
    {
      // cout<<*lheComments<<endl;

      MuonSystem->weight = genWeight;
      NEvents->Fill(1, genWeight);
      // NEvents_genweight->Fill(1);
    }
    MuonSystem->runNum = runNum;
    MuonSystem->lumiSec = lumiNum;
    MuonSystem->evtNum = eventNum;

    bool wzFlag = false;
    // cout<<"ngenparticles: "<<nGenParticle<<endl;
    if (!isData)
    {


      for (int i=0; i < nBunchXing; i++)
      {
        if (BunchXing[i] == 0)
        {
          MuonSystem->npu = nPUmean[i];
        }
      }
      MuonSystem->pileupWeight = helper->getPileupWeight(MuonSystem->npu);

    }


    //get NPU
    MuonSystem->npv = nPV;
    // MuonSystem->rho = fixedGridRhoFastjetAll;
    MuonSystem->met = metType1Pt;
    // MuonSystem->metPhi = metType1Phi;



    std::pair<double,double> corrected_met;
    if (analysisTag=="Razor2016_07Aug2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2016, !isData, nPV);
    else if (analysisTag=="Razor2017_17Nov2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2017, !isData, nPV);
    else if (analysisTag=="Razor2018_17SeptEarlyReReco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2018, !isData, nPV);



    //Triggers
    // for(int i = 0; i < NTriggersMAX; i++){
    //   MuonSystem->HLTDecision[i] = HLTDecision[i];
    // }

    if (analysisTag=="Razor2016_07Aug2017Rereco")
    {
      MuonSystem->METTrigger = HLTDecision[310] || HLTDecision[467];
    }
    else
    {
      // MuonSystem->METTrigger = HLTDecision[310] || HLTDecision[467] || HLTDecision[724] || HLTDecision[729] || HLTDecision[730] || HLTDecision[733];
      MuonSystem->METTrigger = HLTDecision[310] || HLTDecision[467] || HLTDecision[703] || HLTDecision[717] || HLTDecision[710] || HLTDecision[709];

    }
    MuonSystem->Flag2_all = Flag2_HBHENoiseFilter && Flag2_HBHEIsoNoiseFilter && Flag2_BadPFMuonFilter && Flag2_globalSuperTightHalo2016Filter && Flag2_EcalDeadCellTriggerPrimitiveFilter;
    if (isData) MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_eeBadScFilter;

    if (analysisTag!="Razor2016_07Aug2017Rereco")
    {
      MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_ecalBadCalibFilter;
    }



  MuonSystem->HLT_PFMET120_PFMHT120_IDTight = HLTDecision[310];
  MuonSystem->HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = HLTDecision[467];
  MuonSystem->HLT_PFMET140_PFMHT140_IDTight = HLTDecision[703];
  MuonSystem->HLT_PFMET120_PFMHT120_IDTight_PFHT60 = HLTDecision[709];
  MuonSystem->HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = HLTDecision[710];
  MuonSystem->HLT_PFMETNoMu140_PFMHTNoMu140_IDTight = HLTDecision[717];
  MuonSystem->HLT_IsoMu27 = HLTDecision[136];
  MuonSystem->HLT_IsoMu27 = HLTDecision[136];
  MuonSystem->HLT_IsoMu27 = HLTDecision[136];
  if (analysisTag=="Razor2016_07Aug2017Rereco")
  {
    MuonSystem->HLT_IsoEle = HLTDecision[79];
  }
  else if(analysisTag=="Razor2017_17Nov2017Rereco")
  {
    MuonSystem->HLT_IsoEle = HLTDecision[625];
  }
  else{
    MuonSystem->HLT_IsoEle = HLTDecision[87];
  }

  //
  // if (eMuFlag)
  // {
  //   if (!HLTDecision[136])continue;
  //   if (nMuons!=1)continue;
  // }
  // else{
  //   if (!MuonSystem->HLT_IsoEle)continue;
  //   if (nElectrons!=1)continue;
  // }


  if (!isData)
  {
    for (int i=0; i < nGenParticle; i++)
    {

      if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == 24)
      { // choosing only the W->munu events
        wzFlag = true;
        MuonSystem->gLepId = gParticleId[i];
        MuonSystem->gLepPt = gParticlePt[i];
        MuonSystem->gLepEta = gParticleEta[i];
        MuonSystem->gLepE = gParticleE[i];
        MuonSystem->gLepPhi = gParticlePhi[i];
      }
      else if (abs(gParticleId[i]) == 15 && gParticleStatus[i] == 2 && abs(gParticleMotherId[i]) == 24){
        wzFlag = true;
        MuonSystem->gLepId = gParticleId[i];
        MuonSystem->gLepPt = gParticlePt[i];
        MuonSystem->gLepEta = gParticleEta[i];
        MuonSystem->gLepE = gParticleE[i];
        MuonSystem->gLepPhi = gParticlePhi[i];
      }
    }
  }

    // if (!MuonSystem->Flag2_all) continue;
    //*************************************************************************
    //Start Object Selection
    //*************************************************************************

    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
    float MET_NoMuX = 0.0;
    float MET_NoMuY = 0.0;

    for( int i = 0; i < nMuons; i++ )
    {
      float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];



      if(!isMuonPOGTightMuon(i, true, false)) continue;
      if(muonPt[i] < 29) continue;
      // if(muonPt[i] > 100) continue;

      if (muonIso > 0.15) continue;
      // if(fabs(muonEta[i]) > 2.4) continue;

      TLorentzVector thisMuon;
      thisMuon.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i], MU_MASS);


      MET_NoMuX += thisMuon.Px();
      MET_NoMuY += thisMuon.Py();
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
      tmpMuon.passId = isMuonPOGTightMuon(i, true, false);
      tmpMuon.passLooseIso = muonIso<0.25;
      tmpMuon.passTightIso = muonIso<0.15;
      tmpMuon.passVTightIso = muonIso<0.10;
      tmpMuon.passVVTightIso = muonIso<0.05;
      tmpMuon.passVetoId = false;

      Leptons.push_back(tmpMuon);
    }
    // if (eMuFlag && Leptons.size() != 1)continue;
    // if (!eMuFlag && Leptons.size() != 0)continue;

    //-------------------------------
    //Electrons
    //-------------------------------
    for( int i = 0; i < nElectrons; i++ )
    {
      // if(isEGammaPOGTightElectron(i, true, true, true, "2017_94X")!=isEGammaPOGTightElectron(i, true, false, true, "vid"))
      // {
      //   cout<<isEGammaPOGTightElectron(i, true, true, true, "2017_94X")<<","<<isEGammaPOGTightElectron(i, true, false, true, "vid")<<endl;
      //
      //
      // }
      if (!isEGammaPOGTightElectron(i, true, false, true, "vid")) continue;
      if(elePt[i] < 40) continue;
      if(fabs(eleEta[i]) > 2.4) continue;

      //remove overlaps
      bool overlap = false;
      for(auto& lep : Leptons)
      {
        if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
      }
      if(overlap) continue;
      leptons tmpElectron;
      tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);
      tmpElectron.pdgId = 11 * -1 * eleCharge[i];
      tmpElectron.dZ = ele_dZ[i];
      tmpElectron.passId = isEGammaPOGLooseElectron(i, true, true, true, "2017_94X");

      tmpElectron.passVetoId = isEGammaPOGVetoElectron(i, true, true, true, "vid");
      Leptons.push_back(tmpElectron);
    }
    // if (Leptons.size() != 1)continue;


    sort(Leptons.begin(), Leptons.end(), my_largest_pt);

    for ( auto &tmp : Leptons )
    {
      MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
      MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
      MuonSystem->nLeptons++;
    }


    for (int i=0; i < nGenParticle; i++)
    {

      if (abs(gParticleId[i]) == 24)
      { // choosing only the W->munu events
        MuonSystem->gWPt = gParticlePt[i];

      }
    }

    // std::vector<jets> Jets;

    // for(int i = 0; i < nJets; i++)
    // {
    //
    //   //------------------------------------------------------------
    //   //exclude selected muons and electrons from the jet collection
    //   //------------------------------------------------------------
    //   double deltaR = -1;
    //   for(auto& lep : Leptons){
    //     double thisDR = RazorAnalyzer::deltaR(jetEta[i],jetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
    //     if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
    //   }
    //   if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton
    //
    //   //------------------------------------------------------------
    //   //Apply Jet Energy and Resolution Corrections
    //   //------------------------------------------------------------
    //   // double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
    //      // fixedGridRhoFastjetAll, jetJetArea[i] , JetCorrector);
    //  // cout<<"before JEC"<<endl;
    //  double JEC = 1.0;
    //   double jetCorrPt = jetPt[i]*JEC;
    //   double jetCorrE = jetE[i]*JEC;
    //   TLorentzVector thisJet = makeTLorentzVector( jetCorrPt, jetEta[i], jetPhi[i], jetCorrE );
    //
    //   if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
    //
    //   if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
    //   if ( !jetPassIDLoose[i] ) continue;
    //
    //
    //   jets tmpJet;
    //   tmpJet.jet    = thisJet;
    //   tmpJet.passId = jetPassIDTight[i];
    //
    //   Jets.push_back(tmpJet);
    //
    //   }
    //
    //   sort(Jets.begin(), Jets.end(), my_largest_pt_jet);
    //
    //   for ( auto &tmp : Jets )
    //   {
    //     MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
    //     MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
    //     MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
    //     MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
    //
    //     MuonSystem->nJets++;
    //   }

    TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);
  	MuonSystem->metNoMu    = sqrt( pow(PFMET.Px() + MET_NoMuX,2) + pow(PFMET.Py() + MET_NoMuY,2) );



    MuonSystem->tree_->Fill();




  }


  cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
  cout << "Writing output trees..." << endl;
  outFile->cd();
  MuonSystem->tree_->Write();
  NEvents->Write();
  // outFile->Write();
  outFile->Close();


}
