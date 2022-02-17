#include "llp_MuonSystem_bdt.h"
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
  bool passTightIso;
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
  bool jetPassMuFrac;

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


void llp_MuonSystem_bdt::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  //---------------------------
  //options format: MH/MX/ctau/condor: 1000/300/0/1
  // mh can be 3-4 digits, mx is always 3 digits, ctau is one digit(number of zeros), last digit is condor option
  // mh can be 3-4 digits, mx is always 3 digits, ctau is 6 digit last digit is condor option
  bool brem = int(options/100)==1;//1 is muon brem, 2 is OOT data
  bool signalScan = int(options/10)%10 == 1; //middle digit
  int option = options%10;//last digit



  if( isData )
  {
    std::cout << "[INFO]: running on data with option: " << option << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with option: " << option << std::endl;
  }
  if( signalScan )
  {
    std::cout << "[INFO]: running with Signal scan" << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running without Signal scan " << option << std::endl;
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
  //
  // if (label == "zH" || label == "bkg_zH" ){
  //   NTrigger = 4;
  //   muonPt_cut = 15;
  //   elePt_cut = 15;
  //   nLepton_cut = 2;
  // }
  // else{
  //   NTrigger = 2;
  //   muonPt_cut = 25;
  //   elePt_cut = 35;
  //   nLepton_cut = 1;
  // }


  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  map<pair<int,int>, TFile*> Files2D;
  map<pair<int,int>, TTree*> Trees2D;
  map<pair<int,int>, TH1F*> NEvents2D;
  // map<int, TFile*> Files2D;
  // map<int, TTree*> Trees2D;
  // map<int, TH1F*> NEvents2D;

  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile;
  if (isData || !signalScan) outFile = new TFile(outfilename.c_str(), "RECREATE");
  LiteTreeMuonSystem *MuonSystem = new LiteTreeMuonSystem;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();
  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *NEvents_genweight = new TH1F("NEvents_genweight", "NEvents_genweight", 1, 1, 2);

  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
  if(cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1

  // cout << "Getting JEC parameters from " << pathname << endl;
  //
  // std::vector<JetCorrectorParameters> correctionParameters;
  // correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
  // correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
  // correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));
  //
  // FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);




  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  // if (analysisTag == "Razor2018_17SeptEarlyReReco") helper = new RazorHelper("Razor2018_17SeptEarlyReReco", isData, false);
  helper = new RazorHelper(analysisTag, isData, false);



  std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
  std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();
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
    if (!isData && signalScan)
    {

      string mh_substring = lheComments->substr(lheComments->find("MH-")+3);
      int mh = stoi(mh_substring.substr(0,mh_substring.find('_')));
      string mx_substring = lheComments->substr(lheComments->find("MS-")+3);
      int mx = stoi(mx_substring.substr(0,mx_substring.find('_')));
      string ctau_substring = lheComments->substr(lheComments->find("ctauS-")+6);
      int ctau = stoi(ctau_substring.substr(0,ctau_substring.find('_')));
      MuonSystem->mH = mh;
      MuonSystem->mX = mx;
      MuonSystem->ctau = ctau;

      // if (mh2 != mh || mx2!=mx || ctau2!=ctau) continue;
      // cout<<*lheComments<<endl;

      pair<int,int> signalPair = make_pair(mx, ctau);

      // if (Files2D.count(mx) == 0){ //create file and tree
      if (Files2D.count(signalPair) == 0){
        //format file name
        string thisFileName = outfilename;
        thisFileName.erase(thisFileName.end()-5, thisFileName.end());
        // thisFileName += "_" + to_string(mx)  + ".root";
        // Files2D[mx] = new TFile(thisFileName.c_str(), "recreate");
        // Trees2D[mx] =  MuonSystem->tree_->CloneTree(0);
        // NEvents2D[mx] = new TH1F(Form("NEvents%d", mx), "NEvents", 1,0.5,1.5);
        thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";

        Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
        Trees2D[signalPair] =  MuonSystem->tree_->CloneTree(0);
        NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1,0.5,1.5);


        cout << "Created new output file " << thisFileName << endl;
      }
          //Fill NEvents hist
      // NEvents2D[mx]->Fill(1.0, genWeight);

      NEvents2D[signalPair]->Fill(1.0, genWeight);


    }
    if (isData)
    {
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    }
    else
    {
      //NEvents->Fill(genWeight);
      MuonSystem->weight = genWeight;
      NEvents->Fill(1, genWeight);
      // NEvents_genweight->Fill(1);
    }


    //event info
    MuonSystem->runNum = runNum;
    MuonSystem->lumiSec = lumiNum;
    MuonSystem->evtNum = eventNum;

    bool wzFlag = false;
    // cout<<"ngenparticles: "<<nGenParticle<<endl;
    if (!isData)
    {
      for (int i=0; i < nGenParticle; i++)
      {

        if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == wzId)
        { // choosing only the W->munu events
          wzFlag = true;
          MuonSystem->gLepId = gParticleId[i];
          MuonSystem->gLepPt = gParticlePt[i];
          MuonSystem->gLepEta = gParticleEta[i];
          MuonSystem->gLepE = gParticleE[i];
          MuonSystem->gLepPhi = gParticlePhi[i];
        }
        else if (abs(gParticleId[i]) == 15 && gParticleStatus[i] == 2 && abs(gParticleMotherId[i]) == wzId){
          wzFlag = true;
          MuonSystem->gLepId = gParticleId[i];
          MuonSystem->gLepPt = gParticlePt[i];
          MuonSystem->gLepEta = gParticleEta[i];
          MuonSystem->gLepE = gParticleE[i];
          MuonSystem->gLepPhi = gParticlePhi[i];
        }
        if (abs(gParticleId[i])== 25)
        {
          MuonSystem->gHiggsPt = gParticlePt[i];
          MuonSystem->gHiggsEta = gParticleEta[i];
          MuonSystem->gHiggsPhi = gParticlePhi[i];
          MuonSystem->gHiggsE = gParticleE[i];

        }
        MuonSystem->gParticleStatus[MuonSystem->nGenParticle] = gParticleStatus[i];
        MuonSystem->gParticleId[MuonSystem->nGenParticle]  = gParticleId[i];
        MuonSystem->gParticleMotherId[MuonSystem->nGenParticle]  = gParticleMotherId[i];
        MuonSystem->gParticlePt[MuonSystem->nGenParticle]  = gParticlePt[i];
        MuonSystem->gParticleEta[MuonSystem->nGenParticle]  = gParticleEta[i];
        MuonSystem->gParticlePhi[MuonSystem->nGenParticle]  = gParticlePhi[i];
        MuonSystem->gParticleE[MuonSystem->nGenParticle]  = gParticleE[i];
        MuonSystem->nGenParticle++;
        // cout<<"genparticles: "<<MuonSystem->nGenParticle<<endl;


      }
      // cout<<nGenJets<<endl;
      // MuonSystem->nGenJets = 0;
      for(int i=0; i < nGenJets; i++)
      {
        // cout<<genJetE[i]<<","<<MuonSystem->nGenJets<<","<<nGenJets<<endl;
        MuonSystem->genJetE[MuonSystem->nGenJets] = genJetE[i];
        MuonSystem->genJetPt[MuonSystem->nGenJets] = genJetPt[i];
        MuonSystem->genJetEta[MuonSystem->nGenJets] = genJetEta[i];
        MuonSystem->genJetPhi[MuonSystem->nGenJets] = genJetPhi[i];
        // MuonSystem->genJetMET[MuonSystem->nGenJets] = genJetMET[i];
        MuonSystem->nGenJets++;
      }
      MuonSystem->genMetPtTrue = genMetPtTrue;
      MuonSystem->genMetPhiTrue = genMetPhiTrue;
      MuonSystem->genMetPtCalo = genMetPtCalo;
      MuonSystem->genMetPhiCalo = genMetPhiCalo;
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
      for (int i=0; i < nBunchXing; i++)
      {
        if (BunchXing[i] == 0)
        {
          MuonSystem->npu = nPUmean[i];
        }
      }
      MuonSystem->pileupWeight = helper->getPileupWeight(MuonSystem->npu);
      MuonSystem->pileupWeightUp = helper->getPileupWeightUp(MuonSystem->npu) / MuonSystem->pileupWeight;
      MuonSystem->pileupWeightDown = helper->getPileupWeightDown(MuonSystem->npu) / MuonSystem->pileupWeight;
    }

    //get NPU
    MuonSystem->npv = nPV;
    MuonSystem->rho = fixedGridRhoFastjetAll;
    MuonSystem->met = metType1Pt;
    MuonSystem->metPhi = metType1Phi;
    MuonSystem->metJESUp = MuonSystem->met;
    MuonSystem->metJESDown = MuonSystem->met;
    // if (MuonSystem->met < 200) continue;
    //Triggers
    for(int i = 0; i < NTriggersMAX; i++){
      MuonSystem->HLTDecision[i] = HLTDecision[i];
    }
    // flags
    // MuonSystem->Flag_HBHENoiseFilter = Flag_HBHENoiseFilter;
    // MuonSystem->Flag_HBHEIsoNoiseFilter = Flag_HBHEIsoNoiseFilter;
    // MuonSystem->Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
    // MuonSystem->Flag_CSCTightHaloFilter = Flag_CSCTightHaloFilter;
    // MuonSystem->Flag_goodVertices = Flag_goodVertices;
    // MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;
    // MuonSystem->Flag_all = Flag_goodVertices && Flag_HBHEIsoNoiseFilter && Flag_BadPFMuonFilter && Flag_CSCTightHaloFilter && Flag_goodVertices && Flag_ecalBadCalibFilter;


    MuonSystem->Flag2_HBHENoiseFilter = Flag2_HBHENoiseFilter;
    MuonSystem->Flag2_HBHEIsoNoiseFilter = Flag2_HBHEIsoNoiseFilter;
    MuonSystem->Flag2_BadPFMuonFilter = Flag2_BadPFMuonFilter;
    MuonSystem->Flag2_globalSuperTightHalo2016Filter = Flag2_globalSuperTightHalo2016Filter;
    MuonSystem->Flag2_globalTightHalo2016Filter = Flag2_globalTightHalo2016Filter;
    MuonSystem->Flag2_BadChargedCandidateFilter =Flag2_BadChargedCandidateFilter;
    // Flag2_goodVertices = Flag2_goodVertices;
    MuonSystem->Flag2_EcalDeadCellTriggerPrimitiveFilter = Flag2_EcalDeadCellTriggerPrimitiveFilter;
    MuonSystem->Flag2_ecalBadCalibFilter = Flag2_ecalBadCalibFilter;
    MuonSystem->Flag2_eeBadScFilter = Flag2_eeBadScFilter;
    MuonSystem->Flag2_all = Flag2_HBHENoiseFilter && Flag2_HBHEIsoNoiseFilter && Flag2_BadPFMuonFilter && Flag2_globalSuperTightHalo2016Filter && Flag2_EcalDeadCellTriggerPrimitiveFilter;
    if (isData) MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag_eeBadScFilter;

    if (analysisTag!="Razor2016_07Aug2017Rereco")
    {
      MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_ecalBadCalibFilter;
    }


    if (!MuonSystem->Flag2_all) continue;

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
      if(muonPt[i] < zh_lepton1_cut) continue;
      if(fabs(muonEta[i]) > 2.4) continue;
      float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
      //if (muonIso>=0.15) continue;

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
      tmpMuon.passId = isMuonPOGTightMuon(i, true,false);

      tmpMuon.passTightIso = muonIso<0.15;
      // tmpMuon.passId = isLooseMuon(i);
      tmpMuon.passVetoId = false;
      // tmpMuon.passVetoPOGId = false;

      Leptons.push_back(tmpMuon);
    }
    // cout<<"jere"<<endl;

    //-------------------------------
    //Electrons
    //-------------------------------
    for( int i = 0; i < nElectrons; i++ )
    {

      if (!isEGammaPOGLooseElectron(i, true, true, true, "Summer16")) continue;
      if(elePt[i] < zh_lepton1_cut) continue;
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


    //----------------
    //Find Z Candidate
    //----------------

    double ZMass = -999;
    double ZPt = -999;
    double tmpDistToZPole = 9999;
    pair<uint,uint> ZCandidateLeptonIndex;
    bool foundZ = false;
    TLorentzVector ZCandidate;
    double leadingLepPt = 0.0;
    for( uint i = 0; i < Leptons.size(); i++ )
    {
      for( uint j = 0; j < Leptons.size(); j++ )
      {
        if (!( Leptons[i].pdgId == -1*Leptons[j].pdgId )) continue;// same flavor opposite charge
        if (abs(Leptons[i].pdgId)!=13)continue; //only select Z->mumu events
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
          leadingLepPt = max(Leptons[i].lepton.Pt(),Leptons[j].lepton.Pt());
          foundZ = true;
        }
      }
    }

    if (foundZ && Leptons.size() == 2 && leadingLepPt > zh_lepton0_cut)
    {
      MuonSystem->ZMass = ZMass;
      MuonSystem->ZPt   = ZPt;
      MuonSystem->ZEta  = ZCandidate.Eta();
      MuonSystem->ZPhi  = ZCandidate.Phi();
      MuonSystem->ZleptonIndex1 = ZCandidateLeptonIndex.first;
      MuonSystem->ZleptonIndex2 = ZCandidateLeptonIndex.second;
      MuonSystem->category = 2;
    } // endif foundZ
    else{
      for ( unsigned int i = Leptons.size(); i>0; --i )
      {
        int index = i-1;
        if (abs(Leptons[index].pdgId) == 13 && Leptons[index].lepton.Pt() < wh_muonPt_cut)
        {
           Leptons.erase(Leptons.begin() + index);
        }
        else if (abs(Leptons[index].pdgId) == 11 && Leptons[index].lepton.Pt() < wh_elePt_cut)
        {
          Leptons.erase(Leptons.begin() + index);
        }
      }
      if (Leptons.size() == 1) MuonSystem->category = 1;
      else MuonSystem->category = 0;
    }
    for ( auto &tmp : Leptons )
    {
      MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
      MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
      MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
      MuonSystem->lepPassId[MuonSystem->nLeptons] = tmp.passId;
      MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
      if(!isData)
      {
        for (int i=0; i < nGenParticle; i++)
        {
          if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && abs(gParticleMotherId[i]) == 23) MuonSystem->lepFromZ[MuonSystem->nLeptons] = true;
        }
      }
    MuonSystem->nLeptons++;
   }

    TLorentzVector met;
    met.SetPtEtaPhiE(metType1Pt,0,metType1Phi,metType1Pt);
    if ( Leptons.size() > 0 )
    {
      TLorentzVector visible = Leptons[0].lepton;
      MuonSystem->MT = GetMT(visible,met);
    }
    if(brem && (MuonSystem->category!=2))continue; //only select Z->mumu events if selecting muon brem
    if(brem && (MuonSystem->lepPt[MuonSystem->ZleptonIndex1] < 50 || MuonSystem->lepPt[MuonSystem->ZleptonIndex2] < 50)) continue;
    if(brem && !isData && !MuonSystem->lepFromZ[MuonSystem->ZleptonIndex1] && !MuonSystem->lepFromZ[MuonSystem->ZleptonIndex2])continue;
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
    // double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
       // fixedGridRhoFastjetAll, jetJetArea[i] , JetCorrector);
   // cout<<"before JEC"<<endl;
   double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
 						   fixedGridRhoFastjetAll, jetJetArea[i],
 						   runNum,
 						   JetCorrectorIOV,JetCorrector);
    double jetCorrPt = jetPt[i]*JEC;
    double jetCorrE = jetE[i]*JEC;
    // cout<<JEC<<endl;
    // cout<<"corrected pt: "<<jetCorrPt<<", "<<"eta: "<<jetEta[i]<<endl;
    TLorentzVector thisJet = makeTLorentzVector( jetCorrPt, jetEta[i], jetPhi[i], jetCorrE );

    if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
    if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
    if ( !jetPassIDLoose[i] ) continue;

    jets tmpJet;
    tmpJet.jet    = thisJet;
    tmpJet.passId = jetPassIDTight[i];

    //tmpJet.jetPassMuFrac = jetPassMuFrac[i];
    //tmpJet.jetNeutralHadronEnergyFraction = jetNeutralHadronEnergyFraction[i];
    //tmpJet.jetNeutralEMEnergyFraction = jetNeutralEMEnergyFraction[i];
    //tmpJet.jetChargedEMEnergyFraction = jetChargedEMEnergyFraction[i];
    //tmpJet.jetChargedHadronEnergyFraction = jetChargedHadronEnergyFraction[i];

    // calculate jet energy scale uncertainty
    double unc = helper->getJecUnc( jetCorrPt, jetEta[i], runNum ); //use run=999 as default
    double jetPtJESUp = jetCorrPt*(1+unc);
    double jetPtJESDown = jetCorrPt/(1+unc);
    double jetEJESUp = jetCorrE*(1+unc);
    double jetEJESDown = jetCorrE/(1+unc);
    TLorentzVector thisJetJESUp = makeTLorentzVector(jetPtJESUp, jetEta[i], jetPhi[i], jetEJESUp);
    TLorentzVector thisJetJESDown = makeTLorentzVector(jetPtJESDown, jetEta[i], jetPhi[i], jetEJESDown);
    if (jetPtJESUp > 20) MuonSystem->metJESUp +=  -1 * (thisJetJESUp.Pt() - thisJet.Pt());
    if (jetPtJESDown > 20) MuonSystem->metJESDown +=  -1 * (thisJetJESDown.Pt() - thisJet.Pt());


    Jets.push_back(tmpJet);

    }

    //-----------------------------
    //Require at least 2 jets
    //-----------------------------
    // if( Jets.size() < 2 ) continue;
    // if (triggered) trig_lepId_dijet->Fill(1);
    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);
    if (Jets.size()>0)
    {
      MuonSystem->jetMet_dPhi = RazorAnalyzer::deltaPhi(jetPhi[0],metType1Phi);
    }
    else{
      MuonSystem->jetMet_dPhi = -999.;
    }
    double jetMet_dPhiMin_temp = 999.;
    double jetMet_dPhiMin4_temp = 999.;
    for ( auto &tmp : Jets )
    {
      MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
      MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
      MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
      MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
      MuonSystem->jetTime[MuonSystem->nJets] = tmp.time;
      if (jetMet_dPhiMin4_temp > abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi)) && MuonSystem->nJets < 4)
      {
        jetMet_dPhiMin4_temp = abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi));

      }
      if (jetMet_dPhiMin_temp > abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi)))
      {
        jetMet_dPhiMin_temp = abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi));

      }
      MuonSystem->jetPassId[MuonSystem->nJets] = tmp.passId;
      // MuonSystem->ecalNRechits[MuonSystem->nJets] = tmp.ecalNRechits;
      // MuonSystem->ecalNRechits[MuonSystem->nJets] = tmp.ecalRechitE;
      // MuonSystem->jetChargedEMEnergyFraction[MuonSystem->nJets] = tmp.jetChargedEMEnergyFraction;
      // MuonSystem->jetNeutralEMEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralEMEnergyFraction;
      // MuonSystem->jetChargedHadronEnergyFraction[MuonSystem->nJets] = tmp.jetChargedHadronEnergyFraction;
      // MuonSystem->jetNeutralHadronEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralHadronEnergyFraction;
      // MuonSystem->jetPassMuFrac[MuonSystem->nJets] = tmp.jetPassMuFrac;
      float min_deltaR = 15.;
      int index = 999;
      // cout<<"nGenJets"<<nGenJets<<endl;
      for(int i=0; i < nGenJets; i++)
      {

        double current_delta_r = RazorAnalyzer::deltaPhi(genJetPhi[i],jetPhi[MuonSystem->nJets]);
        if (current_delta_r < min_deltaR)
        {
          min_deltaR = current_delta_r;
          index = i;
        }
      }
      if (min_deltaR < 0.4)
      {
        MuonSystem->jet_match_genJet_minDeltaR[MuonSystem->nJets] = min_deltaR;
        MuonSystem->jet_match_genJet_index[MuonSystem->nJets] = index;
        MuonSystem->jet_match_genJet_pt[MuonSystem->nJets] = jetPt[index];
      }
      MuonSystem->nJets++;
    }


    MuonSystem-> jetMet_dPhiMin = jetMet_dPhiMin_temp;
    MuonSystem-> jetMet_dPhiMin4 = jetMet_dPhiMin4_temp;

    //-----------------------------
    // CSC INFO
    //-----------------------------
    vector<Point> points;
    points.clear();
    MuonSystem->nCscRechits  = 0;
    int nCscRechitsChamberPlus11(0), nCscRechitsChamberPlus12(0), nCscRechitsChamberPlus13(0), nCscRechitsChamberPlus21(0), nCscRechitsChamberPlus22(0),
    nCscRechitsChamberPlus31(0), nCscRechitsChamberPlus32(0), nCscRechitsChamberPlus41(0), nCscRechitsChamberPlus42(0),
    nCscRechitsChamberMinus11(0), nCscRechitsChamberMinus12(0), nCscRechitsChamberMinus13(0), nCscRechitsChamberMinus21(0), nCscRechitsChamberMinus22(0),
    nCscRechitsChamberMinus31(0), nCscRechitsChamberMinus32(0), nCscRechitsChamberMinus41(0), nCscRechitsChamberMinus42(0);
    for (int i = 0; i < ncscRechits; i++) {

      Point p;
      p.phi = cscRechitsPhi[i];
      p.eta = cscRechitsEta[i];
      p.x = cscRechitsX[i];
      p.y = cscRechitsY[i];
      p.z = cscRechitsZ[i];
      p.t = cscRechitsTpeak[i];
      p.twire = cscRechitsTwire[i];
      p.station = cscRechitsStation[i];
      p.chamber = cscRechitsChamber[i];
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);
      if (cscRechitsChamber[i] == 11) nCscRechitsChamberPlus11++;
      if (cscRechitsChamber[i] == 12) nCscRechitsChamberPlus12++;
      if (cscRechitsChamber[i] == 13) nCscRechitsChamberPlus13++;
      if (cscRechitsChamber[i] == 21) nCscRechitsChamberPlus21++;
      if (cscRechitsChamber[i] == 22) nCscRechitsChamberPlus22++;
      if (cscRechitsChamber[i] == 31) nCscRechitsChamberPlus31++;
      if (cscRechitsChamber[i] == 32) nCscRechitsChamberPlus32++;
      if (cscRechitsChamber[i] == 41) nCscRechitsChamberPlus41++;
      if (cscRechitsChamber[i] == 42) nCscRechitsChamberPlus42++;
      if (cscRechitsChamber[i] == -11) nCscRechitsChamberMinus11++;
      if (cscRechitsChamber[i] == -12) nCscRechitsChamberMinus12++;
      if (cscRechitsChamber[i] == -13) nCscRechitsChamberMinus13++;
      if (cscRechitsChamber[i] == -21) nCscRechitsChamberMinus21++;
      if (cscRechitsChamber[i] == -22) nCscRechitsChamberMinus22++;
      if (cscRechitsChamber[i] == -31) nCscRechitsChamberMinus31++;
      if (cscRechitsChamber[i] == -32) nCscRechitsChamberMinus32++;
      if (cscRechitsChamber[i] == -41) nCscRechitsChamberMinus41++;
      if (cscRechitsChamber[i] == -42) nCscRechitsChamberMinus42++;




    }
    if (nCscRechitsChamberPlus11 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus12 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus13 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus21 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus22 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus31 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus32 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus41 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberPlus42 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus11 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus12 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus13 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus21 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus22 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus31 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus32 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus41 > 50) MuonSystem->nCscRings++;
    if (nCscRechitsChamberMinus42 > 50) MuonSystem->nCscRings++;
    int nDTRechitsChamberMinus12(0),nDTRechitsChamberMinus11(0),nDTRechitsChamber10(0),nDTRechitsChamberPlus11(0),nDTRechitsChamberPlus12(0),nDTRechitsChamberMinus22(0),
    nDTRechitsChamberMinus21(0),nDTRechitsChamber20(0),nDTRechitsChamberPlus21(0),nDTRechitsChamberPlus22(0),nDTRechitsChamberMinus32(0),nDTRechitsChamberMinus31(0),
    nDTRechitsChamber30(0),nDTRechitsChamberPlus31(0),nDTRechitsChamberPlus32(0),nDTRechitsChamberMinus42(0),nDTRechitsChamberMinus41(0),
    nDTRechitsChamber40(0),nDTRechitsChamberPlus41(0),nDTRechitsChamberPlus42(0);
    for (int i = 0; i < nDtRechits; i++) {


      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus12++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus11++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 0) nDTRechitsChamber10++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus11++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus12++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus22++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus21++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 0) nDTRechitsChamber20++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus21++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus22++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus32++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus31++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 0) nDTRechitsChamber30++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus31++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus32++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus42++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus41++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 0) nDTRechitsChamber40++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus41++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus42++;
    }
    if (nDTRechitsChamberMinus12 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus11 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamber10 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus11 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus12 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus22 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus21 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamber20 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus21 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus22 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus32 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus31 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamber30 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus31 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus32 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus42 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberMinus41 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamber40 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus41 > 50) MuonSystem->nDtRings++;
    if (nDTRechitsChamberPlus42 > 50) MuonSystem->nDtRings++;

    if (MuonSystem->nDtRings+MuonSystem->nCscRings > 10) continue;
    //Do DBSCAN Clustering
    int min_point = 50;  //minimum number of segments to call it a cluster
    float epsilon = 0.2; //cluster radius parameter
    DBSCAN ds(min_point, epsilon, points);
    ds.run();
    ds.result();
    ds.clusterMoments();//but need this step to merge
    ds.sort_clusters();//might not need to sort

    ds.merge_clusters();
    ds.result();
    ds.clusterMoments();
    ds.sort_clusters();

    //Save cluster information
    MuonSystem->nCscRechitClusters = 0;
    for ( auto &tmp : ds.CscCluster ) {
        if (isData && !brem && tmp.tTotal > -12.5)continue;//OOT data
        if (brem && abs(tmp.tTotal) > 12.5)continue;//ZMu data or DY samples
        if (!isData && !brem && abs(tmp.tTotal) > 12.5)continue;//actual signal
        if(abs(tmp.eta)>=2.0) continue;
        if(tmp.TSpread > 20)continue;

        if(!brem  && tmp.nCscSegmentChamberPlus11 != 0)continue;
        if(!brem  && tmp.nCscSegmentChamberPlus12 != 0)continue;
        if(!brem && tmp.nCscSegmentChamberMinus11 != 0)continue;
        if(!brem  && tmp.nCscSegmentChamberMinus12 != 0)continue;
        if (brem && abs(tmp.maxChamber)<=12)continue;



        MuonSystem->cscRechitClusterX[MuonSystem->nCscRechitClusters] =tmp.x;
        MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] =tmp.y;
        MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] =tmp.z;
        MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] = tmp.t;
        MuonSystem->cscRechitClusterTimeTotal[MuonSystem->nCscRechitClusters] = tmp.tTotal;
        MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] =tmp.eta;
        MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters] = tmp.phi;
        MuonSystem->cscRechitClusterMajorAxis[MuonSystem->nCscRechitClusters] =tmp.MajorAxis;
        MuonSystem->cscRechitClusterMinorAxis[MuonSystem->nCscRechitClusters] =tmp.MinorAxis;
        MuonSystem->cscRechitClusterXSpread[MuonSystem->nCscRechitClusters] =tmp.XSpread;
        MuonSystem->cscRechitClusterYSpread[MuonSystem->nCscRechitClusters] =tmp.YSpread;
        MuonSystem->cscRechitClusterZSpread[MuonSystem->nCscRechitClusters] =tmp.ZSpread;
        MuonSystem->cscRechitClusterEtaPhiSpread[MuonSystem->nCscRechitClusters] =tmp.EtaPhiSpread;
        MuonSystem->cscRechitClusterXYSpread[MuonSystem->nCscRechitClusters] =tmp.XYSpread;

        MuonSystem->cscRechitClusterEtaSpread[MuonSystem->nCscRechitClusters] =tmp.EtaSpread;
        MuonSystem->cscRechitClusterPhiSpread[MuonSystem->nCscRechitClusters] = tmp.PhiSpread;
        MuonSystem->cscRechitClusterTimeSpread[MuonSystem->nCscRechitClusters] = tmp.TSpread;
        MuonSystem->cscRechitClusterSize[MuonSystem->nCscRechitClusters] = tmp.nCscSegments;
        //if(!isData)
        //{
        //  int phi_corr_index = 0;
        //  int r_corr_index = 0;
        //  if (abs(tmp.eta) < 1)
        //  {
        //    phi_corr_index = 5;
        //    r_corr_index = 3;
        //  }
        //  else if (abs(tmp.eta) < 1.2)
        //  {
        //    phi_corr_index = 4;
        //    r_corr_index = 2;
        //  }
        //  else if (abs(tmp.eta) < 1.4)
        //  {
        //    phi_corr_index = 2;
        //    r_corr_index = 1;
        //  }
        //  else if (abs(tmp.eta) < 1.6)
        //  {
        //    phi_corr_index = 1;
        //    r_corr_index = 1;
        //  }
        //  else if (abs(tmp.eta) < 1.8)
        //  {
        //    phi_corr_index = 5;
        //    r_corr_index = 1;
        //  }
        //  else if (abs(tmp.eta) < 2)
        //  {
        //    phi_corr_index = 5;
        //    r_corr_index = 0;
        //  }
        //  MuonSystem->cscRechitClusterRSpread_corr[MuonSystem->nCscRechitClusters] = tmp.RSpread_corr[phi_corr_index][r_corr_index];
        //  MuonSystem->cscRechitClusterXYSpread_corr[MuonSystem->nCscRechitClusters] = tmp.XYSpread_corr[phi_corr_index][r_corr_index];
        //  MuonSystem->cscRechitClusterXSpread_corr[MuonSystem->nCscRechitClusters] = tmp.XSpread_corr[phi_corr_index][r_corr_index];
        //  MuonSystem->cscRechitClusterYSpread_corr[MuonSystem->nCscRechitClusters] = tmp.YSpread_corr[phi_corr_index][r_corr_index];
        //  MuonSystem->cscRechitClusterPhiSpread_corr[MuonSystem->nCscRechitClusters] = tmp.PhiSpread_corr[phi_corr_index][r_corr_index];
        //  MuonSystem->cscRechitClusterEtaSpread_corr[MuonSystem->nCscRechitClusters] = tmp.EtaSpread_corr[phi_corr_index][r_corr_index];
        //  MuonSystem->cscRechitClusterEtaPhiSpread_corr[MuonSystem->nCscRechitClusters] = tmp.EtaPhiSpread_corr[phi_corr_index][r_corr_index];
        //  MuonSystem->cscRechitClusterRSpread_corr[MuonSystem->nCscRechitClusters] = tmp.RSpread_corr[phi_corr_index][r_corr_index];

        //}


        MuonSystem->cscRechitClusterNRechitChamberPlus11[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus11;
        MuonSystem->cscRechitClusterNRechitChamberPlus12[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus12;
        MuonSystem->cscRechitClusterNRechitChamberPlus13[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus13;
        MuonSystem->cscRechitClusterNRechitChamberPlus21[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus21;
        MuonSystem->cscRechitClusterNRechitChamberPlus22[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus22;
        MuonSystem->cscRechitClusterNRechitChamberPlus31[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus31;
        MuonSystem->cscRechitClusterNRechitChamberPlus32[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus32;
        MuonSystem->cscRechitClusterNRechitChamberPlus41[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus41;
        MuonSystem->cscRechitClusterNRechitChamberPlus42[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus42;
        MuonSystem->cscRechitClusterNRechitChamberMinus11[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus11;
        MuonSystem->cscRechitClusterNRechitChamberMinus12[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus12;
        MuonSystem->cscRechitClusterNRechitChamberMinus13[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus13;
        MuonSystem->cscRechitClusterNRechitChamberMinus21[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus21;
        MuonSystem->cscRechitClusterNRechitChamberMinus22[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus22;
        MuonSystem->cscRechitClusterNRechitChamberMinus31[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus31;
        MuonSystem->cscRechitClusterNRechitChamberMinus32[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus32;
        MuonSystem->cscRechitClusterNRechitChamberMinus41[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus41;
        MuonSystem->cscRechitClusterNRechitChamberMinus42[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus42;
        MuonSystem->cscRechitClusterMaxChamber[MuonSystem->nCscRechitClusters] = tmp.maxChamber;
        MuonSystem->cscRechitClusterMaxChamberRatio[MuonSystem->nCscRechitClusters] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
        MuonSystem->cscRechitClusterNChamber[MuonSystem->nCscRechitClusters] = tmp.nChamber;
        MuonSystem->cscRechitClusterMaxStation[MuonSystem->nCscRechitClusters] = tmp.maxStation;
        MuonSystem->cscRechitClusterMaxStationRatio[MuonSystem->nCscRechitClusters] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
        MuonSystem->cscRechitClusterNStation[MuonSystem->nCscRechitClusters] = tmp.nStation;
        MuonSystem->cscRechitClusterNStation5[MuonSystem->nCscRechitClusters] = tmp.nStation5;
        MuonSystem->cscRechitClusterNStation10[MuonSystem->nCscRechitClusters] = tmp.nStation10;
        MuonSystem->cscRechitClusterNStation10perc[MuonSystem->nCscRechitClusters] = tmp.nStation10perc;
        MuonSystem->cscRechitClusterAvgStation[MuonSystem->nCscRechitClusters] = tmp.avgStation;
        MuonSystem->cscRechitClusterAvgStation5[MuonSystem->nCscRechitClusters] = tmp.avgStation5;
        MuonSystem->cscRechitClusterAvgStation10[MuonSystem->nCscRechitClusters] = tmp.avgStation10;
        MuonSystem->cscRechitClusterAvgStation10perc[MuonSystem->nCscRechitClusters] = tmp.avgStation10perc;
        MuonSystem->cscRechitClusterMe11Ratio[MuonSystem->nCscRechitClusters] = tmp.Me11Ratio;
        MuonSystem->cscRechitClusterMe12Ratio[MuonSystem->nCscRechitClusters] = tmp.Me12Ratio;
       //  for (unsigned int j = 0; j < tmp.segment_id.size(); j++)
       // {
       //   MuonSystem->cscRechitsCluster2Id[tmp.segment_id[j]] = MuonSystem->nCscRechitClusters;
       // }
       //match to MB1 DT segments
       for (int i = 0; i < nDtSeg; i++) {
         if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
         {
           if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters] ++;
         }

       }
       //match to RPC hits in RE1/2
       for (int i = 0; i < nRpc; i++) {
         float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
         if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
         {
           if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)
           {
             MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters] ++;
           }
           if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
           {
             MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters] ++;
           }

         }

       }
       if ( !brem && MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters]>0) continue;
       if (!brem && MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters]>0) continue;
       if (!brem && MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters]>0) continue;

       if(MuonSystem->category == 2 && brem)
       {
	 if (abs(MuonSystem->lepPdgId[MuonSystem->ZleptonIndex1])!=13)continue;
         bool match1 = RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex1], MuonSystem->lepPhi[MuonSystem->ZleptonIndex1], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4;
         bool match2 = RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex2], MuonSystem->lepPhi[MuonSystem->ZleptonIndex2], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4;
         //if(match1==false && match2 == false) continue;
	 // tag and probe
	 bool sel1 = match1 && MuonSystem->lepPassId[MuonSystem->ZleptonIndex2] && MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex2];
	bool sel2 = match2 && MuonSystem->lepPassId[MuonSystem->ZleptonIndex1] && MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
	if (!(sel1 || sel2))continue;
       }


        //Jet veto/ muon veto
        MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterGenMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterGenMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;

        // jet veto
        for(int i = 0; i < nJets; i++)
        {
          if (fabs(jetEta[i]>3.0)) continue;
          if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] ) {
            MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]  = jetPt[i];
          }
          if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetE[i] > MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] ) {
            MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters]  = jetE[i];
          }
        }

        for(int i = 0; i < nCscRechitClusters; i++)
        {
          if (fabs(jetEta[i]>3.0)) continue;
          if (RazorAnalyzer::deltaR(cscRechitClusterEta[i], cscRechitClusterPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && cscRechitClusterJetVetoPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] ) {
            MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]  = cscRechitClusterJetVetoPt[i];
          }
          if (RazorAnalyzer::deltaR(cscRechitClusterEta[i], cscRechitClusterPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && cscRechitClusterJetVetoE[i] > MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] ) {
            MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters]  = cscRechitClusterJetVetoE[i];
          }
        }
        for(int i = 0; i < nMuons; i++)
        {
          if (fabs(muonEta[i]>3.0)) continue;
          float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
          if (muonIso>0.15)continue;
          if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] ) {
            MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]  = muonPt[i];

          }
          if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonE[i] > MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] ) {
            MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters]  = muonE[i];
          }
        }
        if (brem && MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] < 50)continue;//match to muon the pass tight iso, with pt>50
        if (!brem && MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]>10)continue;
        if (!brem && MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]>20)continue;

        // match to gen-level muon
        if(!isData)
        {
          for(int i = 0; i < nGenParticle; i++)
          {
            if (abs(gParticleId[i])!=13) continue;

            // if (fabs(muonEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->cscRechitClusterGenMuonVetoPt[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterGenMuonVetoPt[MuonSystem->nCscRechitClusters]  = gParticlePt[i];
            }
            if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonE[i] > MuonSystem->cscRechitClusterGenMuonVetoE[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterGenMuonVetoE[MuonSystem->nCscRechitClusters]  = gParticleE[i];
            }
          }
          if (brem && MuonSystem->cscRechitClusterGenMuonVetoPt[MuonSystem->nCscRechitClusters] < 50)continue;
        }

        // match to gen level LLP
        float min_deltaR = 15.;
        int index = 999;
        for(int j = 0; j < 2;j++)
        {

          double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], gLLP_eta[j], gLLP_phi[j]);
          if (current_delta_r < min_deltaR)
          {
            min_deltaR = current_delta_r;
            index = j;
          }
        }
        if (min_deltaR < 0.4)
        {
          MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = true;
          MuonSystem->cscRechitCluster_match_gLLP_minDeltaR[MuonSystem->nCscRechitClusters] = min_deltaR;
          MuonSystem->cscRechitCluster_match_gLLP_index[MuonSystem->nCscRechitClusters] = index;
          MuonSystem->cscRechitCluster_match_gLLP_eta[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_eta[index];
          MuonSystem->cscRechitCluster_match_gLLP_phi[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_phi[index];
          MuonSystem->cscRechitCluster_match_gLLP_decay_r[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_r[index];
          MuonSystem->cscRechitCluster_match_gLLP_decay_x[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_x[index];
          MuonSystem->cscRechitCluster_match_gLLP_decay_y[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_y[index];
          MuonSystem->cscRechitCluster_match_gLLP_decay_z[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
          MuonSystem->cscRechitCluster_match_gLLP_ctau[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_ctau[index];
          MuonSystem->cscRechitCluster_match_gLLP_beta[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_beta[index];
          MuonSystem->cscRechitCluster_match_gLLP_csc[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_csc[index];
          if(!isData && !brem && !MuonSystem->gLLP_csc[index]) continue;
        }
        else{
          if (!isData && !brem) continue;
        }
        MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhi);

        MuonSystem->nCscRechitClusters++;
    }

    if (MuonSystem->nCscRechitClusters == 0) continue;
    /*if(nCscRechitClusters == 0) continue;

    for(int i = 0; i < nCscRechitClusters; i++)
    {
      // if (isData && cscRechitClusterTime[i] > -12.5)continue;
      if (cscRechitClusterTime[i] > -12.5) continue;

      // if (abs(cscRechitClusterMaxChamber[i]) <= 12) continue;
      // if(!cscRechitClusterNRechitChamberPlus11[i] == 0)continue;
      // if(!cscRechitClusterNRechitChamberPlus12[i] == 0)continue;
      // if(!cscRechitClusterNRechitChamberMinus11[i] == 0)continue;
      // if(!cscRechitClusterNRechitChamberMinus12[i] == 0)continue;
      // if(abs(cscRechitClusterEta[i])>=2.1) continue;
      MuonSystem->cscRechitClusterX[MuonSystem->nCscRechitClusters] =cscRechitClusterX[i];
      MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] =cscRechitClusterY[i];
      MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] =cscRechitClusterZ[i];
      MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] =cscRechitClusterTime[i];
      MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] =cscRechitClusterEta[i];
      MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters] =cscRechitClusterPhi[i];
      MuonSystem->cscRechitClusterMajorAxis[MuonSystem->nCscRechitClusters] =cscRechitClusterMajorAxis[i];
      MuonSystem->cscRechitClusterMinorAxis[MuonSystem->nCscRechitClusters] =cscRechitClusterMinorAxis[i];
      MuonSystem->cscRechitClusterXSpread[MuonSystem->nCscRechitClusters] =cscRechitClusterXSpread[i];
      MuonSystem->cscRechitClusterYSpread[MuonSystem->nCscRechitClusters] =cscRechitClusterYSpread[i];
      MuonSystem->cscRechitClusterZSpread[MuonSystem->nCscRechitClusters] =cscRechitClusterZSpread[i];
      MuonSystem->cscRechitClusterEtaPhiSpread[MuonSystem->nCscRechitClusters] =cscRechitClusterEtaPhiSpread[i];
      MuonSystem->cscRechitClusterEtaSpread[MuonSystem->nCscRechitClusters] =cscRechitClusterEtaSpread[i];
      MuonSystem->cscRechitClusterPhiSpread[MuonSystem->nCscRechitClusters] = cscRechitClusterPhiSpread[i];
      MuonSystem->cscRechitClusterTimeSpread[MuonSystem->nCscRechitClusters] =cscRechitClusterTimeSpread[i];

      MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] =cscRechitClusterJetVetoPt[i];
      MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] =cscRechitClusterJetVetoE[i];
      MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] =cscRechitClusterMuonVetoPt[i];
      MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] =cscRechitClusterMuonVetoE[i];

      MuonSystem->cscRechitClusterSize[MuonSystem->nCscRechitClusters] = cscRechitClusterSize[i];
      MuonSystem->cscRechitClusterMaxChamber[MuonSystem->nCscRechitClusters] =cscRechitClusterMaxChamber[i];
      MuonSystem->cscRechitClusterMaxChamberRatio[MuonSystem->nCscRechitClusters] = cscRechitClusterMaxChamberRatio[i];
      MuonSystem->cscRechitClusterNChamber[MuonSystem->nCscRechitClusters] =cscRechitClusterNChamber[i];
      MuonSystem->cscRechitClusterMaxStation[MuonSystem->nCscRechitClusters] = cscRechitClusterMaxStation[i];
      MuonSystem->cscRechitClusterMaxStationRatio[MuonSystem->nCscRechitClusters] = cscRechitClusterMaxStationRatio[i];
      MuonSystem->cscRechitClusterNStation[MuonSystem->nCscRechitClusters] = cscRechitClusterNStation[i];
      MuonSystem->cscRechitClusterMe11Ratio[MuonSystem->nCscRechitClusters] =cscRechitClusterMe11Ratio[i];
      MuonSystem->cscRechitClusterMe12Ratio[MuonSystem->nCscRechitClusters] = cscRechitClusterMe12Ratio[i];

      MuonSystem->cscRechitClusterNRechitChamberPlus11[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus11[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus12[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus12[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus13[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus13[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus21[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus21[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus22[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus22[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus31[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus31[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus32[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus32[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus41[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus41[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus42[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus42[i];

      MuonSystem->cscRechitClusterNRechitChamberMinus11[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus11[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus12[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus12[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus13[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus13[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus21[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus21[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus22[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus22[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus31[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus31[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus32[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus32[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus41[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus41[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus42[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus42[i];
      bool me1112_veto = abs(MuonSystem->cscRechitClusterMaxChamber[MuonSystem->nCscRechitClusters]) > 12;
      // if (MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] < JET_PT_CUT) MuonSystem->nCsc_JetVetoCluster0p4 += cscRechitClusterSize[i];
      // if (MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] < JET_PT_CUT && MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] < MUON_PT_CUT) MuonSystem->nCsc_JetMuonVetoCluster0p4 += cscRechitClusterSize[i];
      // if (MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] < JET_PT_CUT && me1112_veto) MuonSystem->nCsc_JetVetoCluster0p4_Me1112Veto+= cscRechitClusterSize[i];
      if (MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] < JET_PT_CUT && MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] < MUON_PT_CUT && me1112_veto) MuonSystem->nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto+= cscRechitClusterSize[i];
      //
      float min_deltaR = 15.;
      int index = 999;
      for(int j = 0; j < 2;j++)
      {

        double current_delta_r = RazorAnalyzer::deltaR(cscRechitClusterEta[MuonSystem->nCscRechitClusters], cscRechitClusterPhi[MuonSystem->nCscRechitClusters], gLLP_eta[j], gLLP_phi[j]);
        if (current_delta_r < min_deltaR)
        {
          min_deltaR = current_delta_r;
          index = j;
        }
      }
      if (min_deltaR < 0.4)
      {
        MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = true;
        MuonSystem->cscRechitCluster_match_gLLP_minDeltaR[MuonSystem->nCscRechitClusters] = min_deltaR;
        MuonSystem->cscRechitCluster_match_gLLP_index[MuonSystem->nCscRechitClusters] = index;
        // MuonSystem->cscRechitCluster_match_gLLP_pt[MuonSystem->nCscRechitClusters] = index;
        // MuonSystem->cscRechitCluster_match_gLLP_pt[MuonSystem->nCscRechitClusters] = index;

        MuonSystem->cscRechitCluster_match_gLLP_eta[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_eta[index];
        MuonSystem->cscRechitCluster_match_gLLP_phi[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_phi[index];
        MuonSystem->cscRechitCluster_match_gLLP_decay_r[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_r[index];
        MuonSystem->cscRechitCluster_match_gLLP_decay_x[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_x[index];
        MuonSystem->cscRechitCluster_match_gLLP_decay_y[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_y[index];
        MuonSystem->cscRechitCluster_match_gLLP_decay_z[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
        MuonSystem->cscRechitCluster_match_gLLP_ctau[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_ctau[index];
        MuonSystem->cscRechitCluster_match_gLLP_beta[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_beta[index];
        MuonSystem->cscRechitCluster_match_gLLP_csc[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_csc[index];

        // if(!isData && !MuonSystem->gLLP_csc[index]) continue;

      }
      // else{
      //   if (!isData) continue;
      // }
      MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhi);

      MuonSystem->nCscRechitClusters++;
    }
    if (MuonSystem->nCscRechitClusters == 0) continue;
    */

    if(!isData && signalScan)
    {
      // Trees2D[MuonSystem->mX]->Fill();
      pair<int,int> smsPair = make_pair(MuonSystem->mX, MuonSystem->ctau);
      Trees2D[smsPair]->Fill();
    }
    else
    {
      MuonSystem->tree_->Fill();
    }
  }

  if(!isData && signalScan)
  {
    for(auto &filePtr : Files2D)
     {
       cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
       filePtr.second->cd();
       Trees2D[filePtr.first]->Write();
       NEvents2D[filePtr.first]->Write("NEvents");
       filePtr.second->Close();

     }
  }
  else
  {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();

    NEvents->Write();

    outFile->Close();
  }
}
