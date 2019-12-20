#include "llp_MuonSystem_cluster.h"
#include "RazorHelper.h"
#include "LiteTreeMuonSystem.h"
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
};

//lepton highest pt comparator
struct largest_pt
{
  inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt;

//jet highest pt comparator
struct largest_pt_jet
{
  inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;

int cscChamber(double x, double y, double z)
{
  double r = sqrt(x*x+y*y);
  int sign_z = TMath::Sign(1.0, z);
  // if (r > 80 && r < 283 && abs(z) > 568  && abs(z) < 632) return sign_z*11;
  // if (r > 255 && r < 465.0 && abs(z) > 668.3 && abs(z) < 724) return sign_z*12;
  // if (r > 485.5 && r < 695.5 && abs(z) > 686 && abs(z) < 724) return sign_z*13;
  // if (r > 118.5 && r < 345 && abs(z) > 791 && abs(z) < 849.5) return sign_z*21;
  // if (r > 337.5 && r < 695.5 && abs(z) > 791 && abs(z) < 849.5) return sign_z*22;
  // if (r > 140.5 && r < 345 && abs(z) > 911.5 && abs(z) < 970) return sign_z*31;
  // if (r > 337.5 && r < 695.5 && abs(z) > 911.5 && abs(z) < 970) return sign_z*32;
  // if (r > 157.5 && r < 345 && abs(z) > 1002 && abs(z) < 1060.5) return sign_z*41;
  // if (r > 337.5 && r < 695.5 && abs(z) > 1002 && abs(z) < 1060.5) return sign_z*42;
  if (r < 283 && abs(z) > 568  && abs(z) < 632) return sign_z*11;
  if (r < 470.0 && abs(z) > 668.3 && abs(z) < 724) return sign_z*12;
  if (r > 480.0 && abs(z) > 686 && abs(z) < 724) return sign_z*13;
  if (r < 345 && abs(z) > 791 && abs(z) < 849.5) return sign_z*21;
  if (r > 337.5 && abs(z) > 791 && abs(z) < 849.5) return sign_z*22;
  if (r < 345 && abs(z) > 911.5 && abs(z) < 970) return sign_z*31;
  if (r > 337.5 && abs(z) > 911.5 && abs(z) < 970) return sign_z*32;
  if (r < 345 && abs(z) > 1002 && abs(z) < 1060.5) return sign_z*41;
  if (r > 337.5 && abs(z) > 1002 && abs(z) < 1060.5) return sign_z*42;
  return -999;
};
int cscStation(double x, double y, double z)
{
  double r = sqrt(x*x+y*y);
  // z = abs(z);
  int sign_z = TMath::Sign(1.0, z);
  if (r < 283 && abs(z) > 568  && abs(z) < 632) return sign_z*1;
  if (r < 470.0 && abs(z) > 668.3 && abs(z) < 724) return sign_z*1;
  if (r > 480.0 && abs(z) > 686 && abs(z) < 724) return sign_z*1;
  if (r < 345 && abs(z) > 791 && abs(z) < 849.5) return sign_z*2;
  if (r > 337.5 && abs(z) > 791 && abs(z) < 849.5) return sign_z*2;
  if (r < 345 && abs(z) > 911.5 && abs(z) < 970) return sign_z*3;
  if (r > 337.5 && abs(z) > 911.5 && abs(z) < 970) return sign_z*3;
  if (r < 345 && abs(z) > 1002 && abs(z) < 1060.5) return sign_z*4;
  if (r > 337.5 && abs(z) > 1002 && abs(z) < 1060.5) return sign_z*4;
  return -999;
};


void llp_MuonSystem_cluster::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  int option;
  std::string label;
  if (options < 20){
    option = 1; // used when running condor
  }
  else{
    option = 0;// used when running locally
  }
  if (options%10 == 1){
    label = "wH";
  }
  else if (options % 10 == 2){
    label = "zH";
  }
  else if (options % 10 == 3){
    label = "bkg_wH";
  }
  else{
    label = "bkg_zH";
  }
  if( isData )
  {
    std::cout << "[INFO]: running on data with label: " << label << " and option: " << option << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with label: " << label << " and option: " << option << std::endl;
  }

  const float ELE_MASS = 0.000511;
  const float MU_MASS  = 0.105658;
  const float Z_MASS   = 91.2;

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
  else{
    NTrigger = 2;
    muonPt_cut = 25;
    elePt_cut = 35;
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
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
  LiteTreeMuonSystem *MuonSystem = new LiteTreeMuonSystem;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();
  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *NEvents_genweight = new TH1F("NEvents_genweight", "NEvents_genweight", 1, 1, 2);

  // TH1F *generatedEvents = new TH1F("generatedEvents", "generatedEvents", 1, 1, 2);
  // TH1F *trig = new TH1F("trig", "trig", 1, 1, 2);
  // TH1F *trig_lepId = new TH1F("trig_lepId", "trig_lepId", 1, 1, 2);
  // TH1F *trig_lepId_dijet = new TH1F("trig_lepId_dijet", "trig_lepId_dijet", 1, 1, 2);

  cout << "here"  << endl;
  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
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
  clock_t start, end;
  start = clock();
  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {

    //begin event
    // clock_t start, end;
    // double time_taken;
    // start = clock();
    if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;

    if(jentry % 10000 == 0)
    {
      end = clock();
      double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
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

    if (isData)
    {
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    }
    else
    {
      //NEvents->Fill(genWeight);
      MuonSystem->weight = genWeight;
      NEvents->Fill(1);
      NEvents_genweight->Fill(1, genWeight);
      // cout<<genWeight<<endl;
    }


    //event info
    MuonSystem->runNum = runNum;
    MuonSystem->lumiSec = lumiNum;
    MuonSystem->evtNum = eventNum;
    //std::cout << "deb3 " << jentry << std::endl;

    bool wzFlag = false;
    for (int i=0; i < nGenParticle; i++)
    {

      if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == wzId)
      { // choosing only the W->munu events
        wzFlag = true;
        /*MuonSystem->gLepId = gParticleId[i];
        MuonSystem->gLepPt = gParticlePt[i];
        MuonSystem->gLepEta = gParticleEta[i];
        MuonSystem->gLepE = gParticleE[i];
        MuonSystem->gLepPhi = gParticlePhi[i];*/
      }
      else if (abs(gParticleId[i]) == 15 && gParticleStatus[i] == 2 && abs(gParticleMotherId[i]) == wzId){
        wzFlag = true;
        /*MuonSystem->gLepId = gParticleId[i];
        MuonSystem->gLepPt = gParticlePt[i];
        MuonSystem->gLepEta = gParticleEta[i];
        MuonSystem->gLepE = gParticleE[i];
        MuonSystem->gLepPhi = gParticlePhi[i];*/
      }

    }
    for(int i = 0; i < 2;i++)
    {
      MuonSystem->gLLP_eta[i] = gLLP_eta[i];
      MuonSystem->gLLP_phi[i] = gLLP_phi[i];
      MuonSystem->gLLP_decay_vertex_r[i] = sqrt(gLLP_decay_vertex_x[i]*gLLP_decay_vertex_x[i]+gLLP_decay_vertex_y[i]*gLLP_decay_vertex_y[i]);
      MuonSystem->gLLP_decay_vertex_x[i] = gLLP_decay_vertex_x[i];
      MuonSystem->gLLP_decay_vertex_y[i] = gLLP_decay_vertex_y[i];
      MuonSystem->gLLP_decay_vertex_z[i] = gLLP_decay_vertex_z[i];
      float beta = gLLP_beta[i];
      float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[i], 2) + pow(MuonSystem->gLLP_decay_vertex_z[i],2));
      float gamma = 1.0/sqrt(1-beta*beta);
      MuonSystem->gLLP_ctau[i] = gLLP_decay_vertex/(beta * gamma);
      MuonSystem->gLLP_beta[i] = gLLP_beta[i];

      if (abs(MuonSystem->gLLP_eta[i]) < 2.4 && abs(MuonSystem->gLLP_eta[i]) > 0.9
        && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>568
        && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) MuonSystem->gLLP_csc[i] = true;


    }
    // if ( wzFlag == true ) generatedEvents->Fill(1);;
    for (int i=0; i < nBunchXing; i++)
    {
      if (BunchXing[i] == 0)
      {
        MuonSystem->npu = nPUmean[i];
      }
    }
    //get NPU
    MuonSystem->npv = nPV;
    // MuonSystem->rho = fixedGridRhoFastjetAll;
    MuonSystem->met = metType1Pt;
    MuonSystem->metPhi = metType1Phi;

    //Triggers
    for(int i = 0; i < 150; i++){
      MuonSystem->HLTDecision[i] = HLTDecision[i];
    }
    bool triggered = false;
    for(int i = 0; i < NTrigger; i++)
    {
      int trigger_temp = trigger_paths[i];

      triggered = triggered || HLTDecision[trigger_temp];

    }
    // if (triggered) trig->Fill(1);


    //*************************************************************************
    //Start Object Selection
    //*************************************************************************

    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
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
      // tmpMuon.passId = isLooseMuon(i);
      tmpMuon.passVetoId = false;
      // tmpMuon.passVetoPOGId = false;

      Leptons.push_back(tmpMuon);
    }
    //-------------------------------
    //Electrons
    //-------------------------------
    for( int i = 0; i < nElectrons; i++ )
    {

      if (!isEGammaPOGLooseElectron(i, true, true, true, "Summer16")) continue;
      if(elePt[i] < elePt_cut) continue;
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
      tmpElectron.passId = isEGammaPOGTightElectron(i, true, true, true, "Summer16");
      tmpElectron.passVetoId = isEGammaPOGVetoElectron(i, true, true, true, "Summer16");
      Leptons.push_back(tmpElectron);
    }

    sort(Leptons.begin(), Leptons.end(), my_largest_pt);
    //std::cout << "deb7 " << jentry << std::endl;
    for ( auto &tmp : Leptons )
    {
      MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
      // MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
      // MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
      MuonSystem->lepPassId[MuonSystem->nLeptons] = tmp.passId;
      //MuonSystem->lepPassVetoId[MuonSystem->nLeptons] = tmp.passVetoId;

      // std::cout << "lepton pdg " << MuonSystem->lepPdgId[MuonSystem->nLeptons] << std::endl;
      MuonSystem->nLeptons++;
    }

    //----------------
    //Find Z Candidate
    //----------------
    /*
    double ZMass = -999;
    double ZPt = -999;
    double tmpDistToZPole = 9999;
    pair<uint,uint> ZCandidateLeptonIndex;
    bool foundZ = false;
    TLorentzVector ZCandidate;
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
      MuonSystem->ZMass = ZMass;
      MuonSystem->ZPt   = ZPt;
      MuonSystem->ZEta  = ZCandidate.Eta();
      MuonSystem->ZPhi  = ZCandidate.Phi();
      MuonSystem->ZleptonIndex1 = ZCandidateLeptonIndex.first;
      MuonSystem->ZleptonIndex2 = ZCandidateLeptonIndex.second;

    } // endif foundZ
    */
    //------------------------
    //require 1 lepton
    //------------------------

    if ( !(Leptons.size() == nLepton_cut )) continue;
    // end = clock();
    // time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    // if (nCsc >= 30) cout << "Time taken by program is after lepton : " << time_taken << endl;
    TLorentzVector met;
    met.SetPtEtaPhiE(metType1Pt,0,metType1Phi,metType1Pt);
    if ( Leptons.size() > 0 )
    {
      TLorentzVector visible = Leptons[0].lepton;
      MuonSystem->MT = GetMT(visible,met);
    }
    // if (triggered) trig_lepId->Fill(1);


  //-----------------------------------------------
  //Select Jets
  //-----------------------------------------------

  std::vector<jets> Jets;
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

    TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );

    if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
    if( fabs( thisJet.Eta() ) >= 2.4 ) continue;
    if ( !jetPassIDLoose[i] ) continue;
    // if (!(jetRechitE[i] > 0.0)) continue;

      // std::cout <<jetRechitT[i] << "," << jetRechitE[i] <<  "," << jetNRechits[i] << std::endl;


    jets tmpJet;
    tmpJet.jet    = thisJet;
    // tmpJet.time   = jetRechitT[i];
    tmpJet.passId = jetPassIDLoose[i];
    // tmpJet.isCSVL = isCSVL(i);
    //if (isCSVL(i)) NBJet20++;
    //if (isCSVL(i) && thisJet.Pt() > 30) NBJet30++;
    // tmpJet.ecalNRechits = jetNRechits[i];
    // tmpJet.ecalRechitE = jetRechitE[i];
    // tmpJet.jetChargedEMEnergyFraction = jetChargedEMEnergyFraction[i];
    // tmpJet.jetNeutralEMEnergyFraction = jetNeutralEMEnergyFraction[i];
    // tmpJet.jetChargedHadronEnergyFraction = jetChargedHadronEnergyFraction[i];
    // tmpJet.jetNeutralHadronEnergyFraction = jetNeutralHadronEnergyFraction[i];
    Jets.push_back(tmpJet);

    }

    //-----------------------------
    //Require at least 2 jets
    //-----------------------------
    //if( Jets.size() < 2 ) continue;
    // if (triggered) trig_lepId_dijet->Fill(1);
    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

    for ( auto &tmp : Jets )
    {
      // MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
      MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
      // MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
      // MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
      // MuonSystem->jetTime[MuonSystem->nJets] = tmp.time;
      // MuonSystem->jetPassId[MuonSystem->nJets] = tmp.passId;
      // MuonSystem->ecalNRechits[MuonSystem->nJets] = tmp.ecalNRechits;
      // MuonSystem->ecalNRechits[MuonSystem->nJets] = tmp.ecalRechitE;
      // MuonSystem->jetChargedEMEnergyFraction[MuonSystem->nJets] = tmp.jetChargedEMEnergyFraction;
      // MuonSystem->jetNeutralEMEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralEMEnergyFraction;
      // MuonSystem->jetChargedHadronEnergyFraction[MuonSystem->nJets] = tmp.jetChargedHadronEnergyFraction;
      // MuonSystem->jetNeutralHadronEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralHadronEnergyFraction;
      MuonSystem->nJets++;
    }



    //-----------------------------
    // CSC INFO
    //-----------------------------
    //
    // if( nCsc < 10 ) continue;//require at least 30 segments in the CSCs
    // end = clock();
    // time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    // cout << "Time taken by program is after jets : " << time_taken << endl;

    MuonSystem->nCsc = 0;
    vector<Point> points;
    for(int i = 0; i < nCsc; i++)
    {
      // if (cscT[i] > 22.0) continue;
      // if (cscT[i] < -12.5) continue;

      // MuonSystem->cscPhi[MuonSystem->nCsc]           = cscPhi[i];   //[nCsc]
      // MuonSystem->cscEta[MuonSystem->nCsc]           = cscEta[i];   //[nCsc]
      // MuonSystem->cscX[MuonSystem->nCsc]             = cscX[i];   //[nCsc]
      // MuonSystem->cscY[MuonSystem->nCsc]             = cscY[i];   //[nCsc]
      // MuonSystem->cscZ[MuonSystem->nCsc]             = cscZ[i];   //[nCsc]
      // MuonSystem->cscDirectionX[MuonSystem->nCsc]             = cscDirectionX[i];   //
      // MuonSystem->cscDirectionY[MuonSystem->nCsc]             = cscDirectionY[i];   //
      // MuonSystem->cscDirectionZ[MuonSystem->nCsc]             = cscDirectionZ[i];   //
      // MuonSystem->cscStation[MuonSystem->nCsc] = cscStation(cscX[i],cscY[i],cscZ[i]);
      // MuonSystem->cscChamber[MuonSystem->nCsc] = cscChamber(cscX[i],cscY[i],cscZ[i]);

      // for dbscan
      Point p;
      p.phi = cscPhi[i];
      p.eta = cscEta[i];
      p.x = cscX[i];
      p.y = cscY[i];
      p.z = cscZ[i];
      p.t = cscT[i];
      p.dirX = cscDirectionX[i];
      p.dirY = cscDirectionY[i];
      p.dirZ = cscDirectionZ[i];
      // p.station = MuonSystem->cscStation[MuonSystem->nCsc];
      // p.chamber = MuonSystem->cscChamber[MuonSystem->nCsc];
      p.station = cscStation(cscX[i],cscY[i],cscZ[i]);
      p.chamber = cscChamber(cscX[i],cscY[i],cscZ[i]);
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);

      // MuonSystem->cscNRecHits[MuonSystem->nCsc]      = cscNRecHits[i];   //[nCsc]
      // MuonSystem->cscNRecHits_flag[MuonSystem->nCsc] = cscNRecHits_flag[i];   //[nCsc]
      // MuonSystem->cscT[MuonSystem->nCsc]             = cscT[i];   //[nCsc]
      // MuonSystem->cscChi2[MuonSystem->nCsc]          = cscChi2[i];   //[nCsc]

      //page 141 of tdr: https://cds.cern.ch/record/343814/files/LHCC-97-032.pdf

      MuonSystem->nCsc++;
    }



    //****************************************
    // CLUSTERING ALGORITHM WITHOUT TIME CUT
    //****************************************

    int min_point = 10;//6
    float epsilon = 0.2;//100
    //run db scan only with points in the station
    DBSCAN ds(min_point, epsilon, points);
    ds.run();
    ds.result();
    ds.clusterMoments();
    ds.vertexing();
    ds.sort_clusters();

    /*for ( auto &tmp : ds.CscCluster )
    {
      MuonSystem->cscClusterX[MuonSystem->nCscClusters] =tmp.x;
      MuonSystem->cscClusterY[MuonSystem->nCscClusters] =tmp.y;
      MuonSystem->cscClusterZ[MuonSystem->nCscClusters] =tmp.z;
      MuonSystem->cscClusterTime[MuonSystem->nCscClusters] = tmp.t;
      MuonSystem->cscClusterEta[MuonSystem->nCscClusters] =tmp.eta;
      MuonSystem->cscClusterPhi[MuonSystem->nCscClusters] = tmp.phi;
      MuonSystem->cscClusterMajorAxis[MuonSystem->nCscClusters] =tmp.MajorAxis;
      MuonSystem->cscClusterMinorAxis[MuonSystem->nCscClusters] =tmp.MinorAxis;
      MuonSystem->cscClusterXSpread[MuonSystem->nCscClusters] =tmp.XSpread;
      MuonSystem->cscClusterYSpread[MuonSystem->nCscClusters] =tmp.YSpread;
      MuonSystem->cscClusterZSpread[MuonSystem->nCscClusters] =tmp.ZSpread;
      MuonSystem->cscClusterEtaPhiSpread[MuonSystem->nCscClusters] =tmp.EtaPhiSpread;
      MuonSystem->cscClusterEtaSpread[MuonSystem->nCscClusters] =tmp.EtaSpread;
      MuonSystem->cscClusterPhiSpread[MuonSystem->nCscClusters] = tmp.PhiSpread;
      MuonSystem->cscClusterTimeSpread[MuonSystem->nCscClusters] = tmp.TSpread;
      MuonSystem->cscClusterSize[MuonSystem->nCscClusters] = tmp.nCscSegments;

      MuonSystem->cscClusterMaxChamber[MuonSystem->nCscClusters] = tmp.maxChamber;
      MuonSystem->cscClusterMaxChamberRatio[MuonSystem->nCscClusters] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
      MuonSystem->cscClusterNChamber[MuonSystem->nCscClusters] = tmp.nChamber;
      MuonSystem->cscClusterMaxStation[MuonSystem->nCscClusters] = tmp.maxStation;
      MuonSystem->cscClusterMaxStationRatio[MuonSystem->nCscClusters] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
      MuonSystem->cscClusterNStation[MuonSystem->nCscClusters] = tmp.nStation;
      MuonSystem->cscClusterMe11Ratio[MuonSystem->nCscClusters] = tmp.Me11Ratio;
      MuonSystem->cscClusterMe12Ratio[MuonSystem->nCscClusters] = tmp.Me12Ratio;
      MuonSystem->cscClusterVertexR[MuonSystem->nCscClusters] = tmp.vertex_r;
      MuonSystem->cscClusterVertexZ[MuonSystem->nCscClusters] = tmp.vertex_z;
      MuonSystem->cscClusterVertexChi2[MuonSystem->nCscClusters] = tmp.vertex_chi2;
      MuonSystem->cscClusterVertexDis[MuonSystem->nCscClusters] = tmp.vertex_dis;
      MuonSystem->cscClusterVertexN[MuonSystem->nCscClusters] = tmp.vertex_n;
      MuonSystem->cscClusterVertexN1[MuonSystem->nCscClusters] = tmp.vertex_n1;
      MuonSystem->cscClusterVertexN5[MuonSystem->nCscClusters] = tmp.vertex_n5;
      MuonSystem->cscClusterVertexN15[MuonSystem->nCscClusters] = tmp.vertex_n15;
      MuonSystem->cscClusterVertexN20[MuonSystem->nCscClusters] = tmp.vertex_n20;
      MuonSystem->cscClusterVertexN10[MuonSystem->nCscClusters] = tmp.vertex_n10;
      // for (unsigned int j = 0; j < tmp.segment_id.size(); j++)
      // {
      //   MuonSystem->cscLabels[j] = MuonSystem->nCscClusters;
      // } it should be cscLabels[tmp.segment_id[j]]

      //Jet veto/ muon veto
      MuonSystem->cscClusterJetVetoPt[MuonSystem->nCscClusters] = 0.0;
      // MuonSystem->cscClusterCaloJetVeto[MuonSystem->nCscClusters] = 0.0;
      MuonSystem->cscClusterMuonVetoPt[MuonSystem->nCscClusters] = 0.0;

      for (int j = 0; j < nJets; j++)
      {
        // if (jetPt[j]<JET_PT_CUT) continue;
        if (abs(jetEta[j])>3) continue;
        if (RazorAnalyzer::deltaR(tmp.eta, tmp.phi, jetEta[j],jetPhi[j]) < 0.4 && jetPt[j] > MuonSystem->cscClusterJetVetoPt[MuonSystem->nCscClusters]){
          MuonSystem->cscClusterJetVetoPt[MuonSystem->nCscClusters]  = jetPt[j];
          MuonSystem->cscClusterJetVetoE[MuonSystem->nCscClusters]  = jetE[j];

        }
      }
      // for (int j = 0; j < nCaloJets; j++)
      // {
      //   // if (calojetPt[j]<JET_PT_CUT) continue;
      //   if (abs(calojetEta[j])>3) continue;
      //   if (RazorAnalyzer::deltaR(tmp.eta, tmp.phi, calojetEta[j],calojetPhi[j]) < 0.4 && calojetPt[j] > MuonSystem->cscClusterCaloJetVeto[MuonSystem->nCscClusters])
      //   {
      //     MuonSystem->cscClusterCaloJetVeto[MuonSystem->nCscClusters] = calojetPt[j];
      //     MuonSystem->cscClusterCaloJetVetoE[MuonSystem->nCscClusters] = calojetE[j];
      //
      //   }
      // }
      for (int j = 0; j < nMuons; j++)
      {
        // if (muonPt[j]<MUON_PT_CUT) continue;
        if (abs(muonEta[j])>3) continue;
        if (RazorAnalyzer::deltaR(tmp.eta, tmp.phi, muonEta[j],muonPhi[j]) < 0.4 && muonPt[j] > MuonSystem->cscClusterMuonVetoPt[MuonSystem->nCscClusters])
        {
          MuonSystem->cscClusterMuonVetoPt[MuonSystem->nCscClusters] = muonPt[j];
          MuonSystem->cscClusterMuonVetoE[MuonSystem->nCscClusters] = muonE[j];

        }
      }
      bool me1112_veto = MuonSystem->cscClusterMe11Ratio[MuonSystem->nCscClusters] == 0.0 && MuonSystem->cscClusterMe12Ratio[MuonSystem->nCscClusters] == 0.0;
      // if (MuonSystem->cscClusterJetVetoPt[MuonSystem->nCscClusters] < JET_PT_CUT) MuonSystem->nCsc_JetVetoCluster0p4 += tmp.nCscSegments;
      // if (MuonSystem->cscClusterJetVetoPt[MuonSystem->nCscClusters] < JET_PT_CUT && MuonSystem->cscClusterMuonVetoPt[MuonSystem->nCscClusters] < MUON_PT_CUT) MuonSystem->nCsc_JetMuonVetoCluster0p4 += tmp.nCscSegments;
      // if (MuonSystem->cscClusterJetVetoPt[MuonSystem->nCscClusters] < JET_PT_CUT && me1112_veto) MuonSystem->nCsc_JetVetoCluster0p4_Me1112Veto+= tmp.nCscSegments;
      // if (MuonSystem->cscClusterJetVetoPt[MuonSystem->nCscClusters] < JET_PT_CUT && MuonSystem->cscClusterMuonVetoPt[MuonSystem->nCscClusters] < MUON_PT_CUT && me1112_veto) MuonSystem->nCsc_JetMuonVetoCluster0p4_Me1112Veto+= tmp.nCscSegments;



      MuonSystem->nCscClusters++;
    }
    if(MuonSystem->nCscClusters == 0) continue;
    */
    //****************************************
    // CLUSTERING ALGORITHM WITH TIME CUT
    //****************************************
    /*
    points.clear();
    for(int i = 0; i < nCsc; i++)
    {
      if(isData){
        if (cscT[i] > -28.0) continue;
        if (cscT[i] < -62.5) continue;
      }
      else{
        if (cscT[i] > 22.0) continue;
        if (cscT[i] < -12.5) continue;
      }

      // for dbscan
      // Point p;
      // p.phi = MuonSystem->cscPhi[i];
      // p.eta = MuonSystem->cscEta[i];
      // p.x = MuonSystem->cscX[i];
      // p.y = MuonSystem->cscY[i];
      // p.z = MuonSystem->cscZ[i];
      // p.t = MuonSystem->cscT[i];
      // p.dirX = MuonSystem->cscDirectionX[i];
      // p.dirY = MuonSystem->cscDirectionY[i];
      // p.dirZ = MuonSystem->cscDirectionZ[i];
      // p.station = MuonSystem->cscStation[i];
      // p.chamber = MuonSystem->cscChamber[i];
      // p.clusterID = UNCLASSIFIED;
      // points.push_back(p);
      Point p;
      p.phi = cscPhi[i];
      p.eta = cscEta[i];
      p.x = cscX[i];
      p.y = cscY[i];
      p.z = cscZ[i];
      p.t = cscT[i];
      p.dirX = cscDirectionX[i];
      p.dirY = cscDirectionY[i];
      p.dirZ = cscDirectionZ[i];
      p.station = cscStation(cscX[i],cscY[i],cscZ[i]);
      p.chamber = cscChamber(cscX[i],cscY[i],cscZ[i]);
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);

    }

    //run db scan only with points in time
    DBSCAN ds_it(min_point, epsilon, points);
    ds_it.run();
    ds_it.result();
    ds_it.clusterMoments();
    ds_it.vertexing();
    ds_it.sort_clusters();
    for ( auto &tmp : ds_it.CscCluster )
    {
      MuonSystem->cscITClusterX[MuonSystem->nCscITClusters] =tmp.x;
      MuonSystem->cscITClusterY[MuonSystem->nCscITClusters] =tmp.y;
      MuonSystem->cscITClusterZ[MuonSystem->nCscITClusters] =tmp.z;
      MuonSystem->cscITClusterTime[MuonSystem->nCscITClusters] = tmp.t;
      MuonSystem->cscITClusterEta[MuonSystem->nCscITClusters] =tmp.eta;
      MuonSystem->cscITClusterPhi[MuonSystem->nCscITClusters] = tmp.phi;
      MuonSystem->cscITClusterMajorAxis[MuonSystem->nCscITClusters] =tmp.MajorAxis;
      MuonSystem->cscITClusterMinorAxis[MuonSystem->nCscITClusters] =tmp.MinorAxis;
      MuonSystem->cscITClusterXSpread[MuonSystem->nCscITClusters] =tmp.XSpread;
      MuonSystem->cscITClusterYSpread[MuonSystem->nCscITClusters] =tmp.YSpread;
      MuonSystem->cscITClusterZSpread[MuonSystem->nCscITClusters] =tmp.ZSpread;
      MuonSystem->cscITClusterEtaPhiSpread[MuonSystem->nCscITClusters] =tmp.EtaPhiSpread;
      MuonSystem->cscITClusterEtaSpread[MuonSystem->nCscITClusters] =tmp.EtaSpread;
      MuonSystem->cscITClusterPhiSpread[MuonSystem->nCscITClusters] = tmp.PhiSpread;
      MuonSystem->cscITClusterTimeSpread[MuonSystem->nCscITClusters] = tmp.TSpread;
      MuonSystem->cscITClusterSize[MuonSystem->nCscITClusters] = tmp.nCscSegments;

      MuonSystem->cscITClusterMaxChamber[MuonSystem->nCscITClusters] = tmp.maxChamber;
      MuonSystem->cscITClusterMaxChamberRatio[MuonSystem->nCscITClusters] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
      MuonSystem->cscITClusterNChamber[MuonSystem->nCscITClusters] = tmp.nChamber;
      MuonSystem->cscITClusterMaxStation[MuonSystem->nCscITClusters] = tmp.maxStation;
      MuonSystem->cscITClusterMaxStationRatio[MuonSystem->nCscITClusters] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
      MuonSystem->cscITClusterNStation[MuonSystem->nCscITClusters] = tmp.nStation;

      MuonSystem->cscITClusterMe11Ratio[MuonSystem->nCscITClusters] = tmp.Me11Ratio;
      MuonSystem->cscITClusterMe12Ratio[MuonSystem->nCscITClusters] = tmp.Me12Ratio;
      MuonSystem->cscITClusterVertexR[MuonSystem->nCscITClusters] = tmp.vertex_r;
      MuonSystem->cscITClusterVertexZ[MuonSystem->nCscITClusters] = tmp.vertex_z;
      MuonSystem->cscITClusterVertexChi2[MuonSystem->nCscITClusters] = tmp.vertex_chi2;
      MuonSystem->cscITClusterVertexDis[MuonSystem->nCscITClusters] = tmp.vertex_dis;
      MuonSystem->cscITClusterVertexN[MuonSystem->nCscITClusters] = tmp.vertex_n;
      MuonSystem->cscITClusterVertexN1[MuonSystem->nCscITClusters] = tmp.vertex_n1;
      MuonSystem->cscITClusterVertexN5[MuonSystem->nCscITClusters] = tmp.vertex_n5;
      MuonSystem->cscITClusterVertexN15[MuonSystem->nCscITClusters] = tmp.vertex_n15;
      MuonSystem->cscITClusterVertexN20[MuonSystem->nCscITClusters] = tmp.vertex_n20;
      MuonSystem->cscITClusterVertexN10[MuonSystem->nCscITClusters] = tmp.vertex_n10;
      // for (unsigned int j = 0; j < tmp.segment_id.size(); j++)
      // {
      //   MuonSystem->cscITLabels[j] = MuonSystem->nCscClusters;
      // }


      //Jet veto / muon veto
      MuonSystem->cscITClusterJetVeto[MuonSystem->nCscITClusters] = 0.0;
      MuonSystem->cscITClusterCaloJetVeto[MuonSystem->nCscITClusters] = 0.0;
      MuonSystem->cscITClusterMuonVeto[MuonSystem->nCscITClusters] = 0.0;

      for (int j = 0; j < nJets; j++)
      {
        // if (jetPt[j]<JET_PT_CUT) continue;
        if (abs(jetEta[j])>3) continue;
        if (RazorAnalyzer::deltaR(tmp.eta, tmp.phi, jetEta[j],jetPhi[j]) < 0.4 && jetPt[j] > MuonSystem->cscITClusterJetVeto[MuonSystem->nCscITClusters]){
          MuonSystem->cscITClusterJetVeto[MuonSystem->nCscITClusters]  = jetPt[j];
          MuonSystem->cscITClusterJetVetoE[MuonSystem->nCscITClusters]  = jetE[j];

        }
      }
      for (int j = 0; j < nCaloJets; j++)
      {
        // if (calojetPt[j]<JET_PT_CUT) continue;
        if (abs(calojetEta[j])>3) continue;
        if (RazorAnalyzer::deltaR(tmp.eta, tmp.phi, calojetEta[j],calojetPhi[j]) < 0.4 && calojetPt[j] > MuonSystem->cscITClusterCaloJetVeto[MuonSystem->nCscITClusters])
        {
          MuonSystem->cscITClusterCaloJetVeto[MuonSystem->nCscITClusters] = calojetPt[j];
          MuonSystem->cscITClusterCaloJetVetoE[MuonSystem->nCscITClusters] = calojetE[j];

        }
      }
      for (int j = 0; j < nMuons; j++)
      {
        // if (muonPt[j]<MUON_PT_CUT) continue;
        if (abs(muonEta[j])>3) continue;
        if (RazorAnalyzer::deltaR(tmp.eta, tmp.phi, muonEta[j],muonPhi[j]) < 0.4 && muonPt[j] > MuonSystem->cscITClusterMuonVeto[MuonSystem->nCscITClusters])
        {
          MuonSystem->cscITClusterMuonVeto[MuonSystem->nCscITClusters] = muonPt[j];
          MuonSystem->cscITClusterMuonVetoE[MuonSystem->nCscITClusters] = muonE[j];

        }
      }
      bool me1112_veto = MuonSystem->cscITClusterMe11Ratio[MuonSystem->nCscITClusters] == 0.0 && MuonSystem->cscITClusterMe12Ratio[MuonSystem->nCscITClusters] == 0.0;
      if (MuonSystem->cscITClusterJetVeto[MuonSystem->nCscITClusters] < JET_PT_CUT) MuonSystem->nCsc_JetVetoITCluster0p4 += tmp.nCscSegments;
      if (MuonSystem->cscITClusterJetVeto[MuonSystem->nCscITClusters] < JET_PT_CUT && MuonSystem->cscITClusterMuonVeto[MuonSystem->nCscITClusters] < MUON_PT_CUT) MuonSystem->nCsc_JetMuonVetoITCluster0p4 += tmp.nCscSegments;
      if (MuonSystem->cscITClusterJetVeto[MuonSystem->nCscITClusters] < JET_PT_CUT && me1112_veto) MuonSystem->nCsc_JetVetoITCluster0p4_Me1112Veto+= tmp.nCscSegments;
      if (MuonSystem->cscITClusterJetVeto[MuonSystem->nCscITClusters] < JET_PT_CUT && MuonSystem->cscITClusterMuonVeto[MuonSystem->nCscITClusters] < MUON_PT_CUT && me1112_veto) MuonSystem->nCsc_JetMuonVetoITCluster0p4_Me1112Veto+= tmp.nCscSegments;

      MuonSystem->nCscITClusters++;
    }
    //match cscITClusters with cscClusters
    for(int i = 0; i < MuonSystem->nCscITClusters; i++)
    {
     float minDeltaR = 100000;
     for(int j = 0;j < MuonSystem->nCscClusters;j++)
     {
       float deltaR_temp = RazorAnalyzer::deltaR(MuonSystem->cscITClusterEta[i],MuonSystem->cscITClusterPhi[i],MuonSystem->cscClusterEta[j],MuonSystem->cscClusterPhi[j]);
       if ( deltaR_temp < minDeltaR)
       {
         minDeltaR = deltaR_temp;
         MuonSystem->cscITCluster_match_cscCluster_index[i] = j;
         MuonSystem->cscITCluster_cscCluster_SizeRatio[i] = 1.0* MuonSystem->cscITClusterSize[i]/MuonSystem->cscClusterSize[j];
       }
     }
   }*/
    MuonSystem->tree_->Fill();
  }

    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();
}
