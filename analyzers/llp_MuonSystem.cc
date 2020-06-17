#include "llp_MuonSystem.h"
#include "RazorHelper.h"
#include "LiteTreeMuonSystem.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
// #include "DBSCAN.h"
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
int dtRing(double x, double y, double z)
{
  double r = sqrt(x*x+y*y);
  int sign_z = TMath::Sign(1.0, z);
  if ((r > 402 && r < 449) || (r > 490.5 && r < 533.5) || (r > 597.5 && r < 636) || (r > 700 && r < 738)){
    if (abs(z) < 126.8) return 0;
    if (abs(z) > 126.8 && abs(z) < 395.4) return sign_z*1;
    if (abs(z) > 395.5 && abs(z) < 661) return sign_z*2;
  }
  return -999;
}
int dtStation(double x, double y, double z)
{
  double r = sqrt(x*x+y*y);
  int sign_z = TMath::Sign(1.0, z);
  if (abs(z) < 661){
    if (r > 402 && r < 449) return sign_z*1;
    if (r > 490.5 && r < 533.5) return sign_z*2;
    if (r > 597.5 && r < 636) return sign_z*3;
    if (r > 700 && r < 738) return sign_z*4;
  }
  return -999;
}

void llp_MuonSystem::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  bool signalScan = int(options/10) == 1;
  int option = options%10;
  if (options == 1){
    option = 1; // used when running condor
  }
  else{
    option = 0;// used when running locally
  }


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
    std::cout << "[INFO]: running without Signal scan" << std::endl;
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
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile;
  if (isData || !signalScan) outFile = new TFile(outfilename.c_str(), "RECREATE");
  
  LiteTreeMuonSystem *MuonSystem = new LiteTreeMuonSystem;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  // for signals, need one output file for each signal point
  map<pair<int,int>, TFile*> Files2D;
  map<pair<int,int>, TTree*> Trees2D;
  map<pair<int,int>, TH1F*> NEvents2D;
  // map<pair<int,int>, TH1F*> smsSumWeights2D;
  // map<pair<int,int>, TH1F*> smsSumScaleWeights2D;
  // map<pair<int,int>, TH1F*> smsSumPdfWeights2D;
  // map<pair<int,int>, TH1F*> smsNISRJets2D;
  // map<pair<int,int>, TH1F*> smsPtISR2D;
  // map<pair<int,int>, TH1F*> smsNPV2D;

  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  //TH1F *NEvents_genweight = new TH1F("NEvents_genweight", "NEvents_genweight", 1, 1, 2);

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
  //if (analysisTag == "Razor2018_17SeptEarlyReReco") helper = new RazorHelper("Razor2018_17SeptEarlyReReco", isData, false);
  //else helper = new RazorHelper(analysisTag, isData, false);
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
      
      if (Files2D.count(signalPair) == 0){ //create file and tree
	//format file name
	string thisFileName = outfilename;
	thisFileName.erase(thisFileName.end()-5, thisFileName.end());
	thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";
	
	Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
	Trees2D[signalPair] =  MuonSystem->tree_->CloneTree(0);
	NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1,0.5,1.5);
	
	
	cout << "Created new output file " << thisFileName << endl;
      }
      //Fill NEvents hist
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
      //NEvents_genweight->Fill(1, genWeight);
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
        if (abs(gParticleId[i])== 25 || abs(gParticleId[i]) == 35)
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

      MuonSystem->higgsPtWeight = helper->getHiggsPtWeight(MuonSystem->gHiggsPt);
      for (int i = 0; i < 9; i++)
      {
	MuonSystem->higgsPtWeightSys[i] = helper->getHiggsPtWeightSys(MuonSystem->gHiggsPt, i) / MuonSystem->higgsPtWeight;

      }
      MuonSystem->sf_facScaleUp = MuonSystem->higgsPtWeightSys[5];
      MuonSystem->sf_facScaleDown = MuonSystem->higgsPtWeightSys[3];
      MuonSystem->sf_renScaleUp = MuonSystem->higgsPtWeightSys[7];
      MuonSystem->sf_renScaleDown = MuonSystem->higgsPtWeightSys[1];
      MuonSystem->sf_facRenScaleUp = MuonSystem->higgsPtWeightSys[8];
      MuonSystem->sf_facRenScaleDown = MuonSystem->higgsPtWeightSys[0];
      // cout<<nGenJets<<endl;
      // MuonSystem->nGenJets = 0;
      for(int i=0; i < nGenJets; i++)
      {
        // cout<<genJetE[i]<<","<<MuonSystem->nGenJets<<","<<nGenJets<<endl;
        MuonSystem->genJetE[MuonSystem->nGenJets] = genJetE[i];
        MuonSystem->genJetPt[MuonSystem->nGenJets] = genJetPt[i];
        MuonSystem->genJetEta[MuonSystem->nGenJets] = genJetEta[i];
        MuonSystem->genJetPhi[MuonSystem->nGenJets] = genJetPhi[i];
        //MuonSystem->genJetMET[MuonSystem->nGenJets] = genJetMET[i];
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
      if (abs(MuonSystem->gLLP_decay_vertex_z[i])< 661.0
        && MuonSystem->gLLP_decay_vertex_r[i] < 738.0
	&& MuonSystem->gLLP_decay_vertex_r[i] > 380.0) MuonSystem->gLLP_dt[i] = true;
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
    if (MuonSystem->met < 200) continue;
    //Triggers
    for(int i = 0; i < NTriggersMAX; i++){
      MuonSystem->HLTDecision[i] = HLTDecision[i];
    }
    // flags
    MuonSystem->Flag_HBHENoiseFilter = Flag_HBHENoiseFilter;
    MuonSystem->Flag_HBHEIsoNoiseFilter = Flag_HBHEIsoNoiseFilter;
    MuonSystem->Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
    //MuonSystem->Flag_CSCTightHaloFilter = Flag_CSCTightHaloFilter;
    MuonSystem->Flag_globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
    MuonSystem->Flag_goodVertices = Flag_goodVertices;
    MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;
    MuonSystem->Flag_BadChargedCandidateFilter = Flag_BadChargedCandidateFilter;
    MuonSystem->Flag_eeBadScFilter = Flag_eeBadScFilter;
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
    if (isData) MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_eeBadScFilter;

    if (analysisTag!="Razor2016_07Aug2017Rereco")
    {
      MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_ecalBadCalibFilter;
    }

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

    if (foundZ && fabs(ZMass-Z_MASS) < 30.0 && Leptons.size() == 2 && leadingLepPt > zh_lepton0_cut)
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
      MuonSystem->nLeptons++;
    }

    TLorentzVector met;
    met.SetPtEtaPhiE(metType1Pt,0,metType1Phi,metType1Pt);
    if ( Leptons.size() > 0 )
    {
      TLorentzVector visible = Leptons[0].lepton;
      MuonSystem->MT = GetMT(visible,met);
    }


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

    tmpJet.jetPassMuFrac = jetPassMuFrac[i];
    tmpJet.jetNeutralHadronEnergyFraction = jetNeutralHadronEnergyFraction[i];
    tmpJet.jetNeutralEMEnergyFraction = jetNeutralEMEnergyFraction[i];
    tmpJet.jetChargedEMEnergyFraction = jetChargedEMEnergyFraction[i];
    tmpJet.jetChargedHadronEnergyFraction = jetChargedHadronEnergyFraction[i];

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
      MuonSystem->jetChargedEMEnergyFraction[MuonSystem->nJets] = tmp.jetChargedEMEnergyFraction;
      MuonSystem->jetNeutralEMEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralEMEnergyFraction;
      MuonSystem->jetChargedHadronEnergyFraction[MuonSystem->nJets] = tmp.jetChargedHadronEnergyFraction;
      MuonSystem->jetNeutralHadronEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralHadronEnergyFraction;
      MuonSystem->jetPassMuFrac[MuonSystem->nJets] = tmp.jetPassMuFrac;
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
    //if(nCscClusters == 0) continue;

    /*for(int i = 0; i < nCscClusters; i++)
    {
      MuonSystem->cscClusterX[MuonSystem->nCscClusters] =cscClusterX[i];
      MuonSystem->cscClusterY[MuonSystem->nCscClusters] =cscClusterY[i];
      MuonSystem->cscClusterZ[MuonSystem->nCscClusters] =cscClusterZ[i];
      MuonSystem->cscClusterTime[MuonSystem->nCscClusters] =cscClusterTime[i];
      MuonSystem->cscClusterEta[MuonSystem->nCscClusters] =cscClusterEta[i];
      MuonSystem->cscClusterPhi[MuonSystem->nCscClusters] =cscClusterPhi[i];
      MuonSystem->cscClusterMajorAxis[MuonSystem->nCscClusters] =cscClusterMajorAxis[i];
      MuonSystem->cscClusterMinorAxis[MuonSystem->nCscClusters] =cscClusterMinorAxis[i];
      MuonSystem->cscClusterXSpread[MuonSystem->nCscClusters] =cscClusterXSpread[i];
      MuonSystem->cscClusterYSpread[MuonSystem->nCscClusters] =cscClusterYSpread[i];
      MuonSystem->cscClusterZSpread[MuonSystem->nCscClusters] =cscClusterZSpread[i];
      MuonSystem->cscClusterEtaPhiSpread[MuonSystem->nCscClusters] =cscClusterEtaPhiSpread[i];
      MuonSystem->cscClusterEtaSpread[MuonSystem->nCscClusters] =cscClusterEtaSpread[i];
      MuonSystem->cscClusterPhiSpread[MuonSystem->nCscClusters] = cscClusterPhiSpread[i];
      MuonSystem->cscClusterTimeSpread[MuonSystem->nCscClusters] =cscClusterTimeSpread[i];

      MuonSystem->cscClusterJetVetoPt[MuonSystem->nCscClusters] =cscClusterJetVetoPt[i];
      MuonSystem->cscClusterJetVetoE[MuonSystem->nCscClusters] =cscClusterJetVetoE[i];
      MuonSystem->cscClusterMuonVetoPt[MuonSystem->nCscClusters] =cscClusterMuonVetoPt[i];
      MuonSystem->cscClusterMuonVetoE[MuonSystem->nCscClusters] =cscClusterMuonVetoE[i];

      MuonSystem->cscClusterSize[MuonSystem->nCscClusters] = cscClusterSize[i];
      MuonSystem->cscClusterMaxChamber[MuonSystem->nCscClusters] =cscClusterMaxChamber[i];
      MuonSystem->cscClusterMaxChamberRatio[MuonSystem->nCscClusters] = cscClusterMaxChamberRatio[i];
      MuonSystem->cscClusterNChamber[MuonSystem->nCscClusters] =cscClusterNChamber[i];
      MuonSystem->cscClusterMaxStation[MuonSystem->nCscClusters] = cscClusterMaxStation[i];
      MuonSystem->cscClusterMaxStationRatio[MuonSystem->nCscClusters] = cscClusterMaxStationRatio[i];
      MuonSystem->cscClusterNStation[MuonSystem->nCscClusters] = cscClusterNStation[i];
      MuonSystem->cscClusterMe11Ratio[MuonSystem->nCscClusters] =cscClusterMe11Ratio[i];
      MuonSystem->cscClusterMe12Ratio[MuonSystem->nCscClusters] = cscClusterMe12Ratio[i];

      MuonSystem->cscClusterNSegmentChamberPlus11[MuonSystem->nCscClusters] = cscClusterNSegmentChamberPlus11[i];
      MuonSystem->cscClusterNSegmentChamberPlus12[MuonSystem->nCscClusters] = cscClusterNSegmentChamberPlus12[i];
      MuonSystem->cscClusterNSegmentChamberPlus13[MuonSystem->nCscClusters] = cscClusterNSegmentChamberPlus13[i];
      MuonSystem->cscClusterNSegmentChamberPlus21[MuonSystem->nCscClusters] = cscClusterNSegmentChamberPlus21[i];
      MuonSystem->cscClusterNSegmentChamberPlus22[MuonSystem->nCscClusters] = cscClusterNSegmentChamberPlus22[i];
      MuonSystem->cscClusterNSegmentChamberPlus31[MuonSystem->nCscClusters] = cscClusterNSegmentChamberPlus31[i];
      MuonSystem->cscClusterNSegmentChamberPlus32[MuonSystem->nCscClusters] = cscClusterNSegmentChamberPlus32[i];
      MuonSystem->cscClusterNSegmentChamberPlus41[MuonSystem->nCscClusters] = cscClusterNSegmentChamberPlus41[i];
      MuonSystem->cscClusterNSegmentChamberPlus42[MuonSystem->nCscClusters] = cscClusterNSegmentChamberPlus42[i];

      MuonSystem->cscClusterNSegmentChamberMinus11[MuonSystem->nCscClusters] = cscClusterNSegmentChamberMinus11[i];
      MuonSystem->cscClusterNSegmentChamberMinus12[MuonSystem->nCscClusters] = cscClusterNSegmentChamberMinus12[i];
      MuonSystem->cscClusterNSegmentChamberMinus13[MuonSystem->nCscClusters] = cscClusterNSegmentChamberMinus13[i];
      MuonSystem->cscClusterNSegmentChamberMinus21[MuonSystem->nCscClusters] = cscClusterNSegmentChamberMinus21[i];
      MuonSystem->cscClusterNSegmentChamberMinus22[MuonSystem->nCscClusters] = cscClusterNSegmentChamberMinus22[i];
      MuonSystem->cscClusterNSegmentChamberMinus31[MuonSystem->nCscClusters] = cscClusterNSegmentChamberMinus31[i];
      MuonSystem->cscClusterNSegmentChamberMinus32[MuonSystem->nCscClusters] = cscClusterNSegmentChamberMinus32[i];
      MuonSystem->cscClusterNSegmentChamberMinus41[MuonSystem->nCscClusters] = cscClusterNSegmentChamberMinus41[i];
      MuonSystem->cscClusterNSegmentChamberMinus42[MuonSystem->nCscClusters] = cscClusterNSegmentChamberMinus42[i];
      bool me1112_veto = MuonSystem->cscClusterMe11Ratio[MuonSystem->nCscClusters] == 0.0 && MuonSystem->cscClusterMe12Ratio[MuonSystem->nCscClusters] == 0.0;
      // if (MuonSystem->cscClusterJetVetoPt[MuonSystem->nCscClusters] < JET_PT_CUT) MuonSystem->nCsc_JetVetoCluster0p4 += cscClusterSize[i];
      // if (MuonSystem->cscClusterJetVetoPt[MuonSystem->nCscClusters] < JET_PT_CUT && MuonSystem->cscClusterMuonVetoPt[MuonSystem->nCscClusters] < MUON_PT_CUT) MuonSystem->nCsc_JetMuonVetoCluster0p4 += cscClusterSize[i];
      // if (MuonSystem->cscClusterJetVetoPt[MuonSystem->nCscClusters] < JET_PT_CUT && me1112_veto) MuonSystem->nCsc_JetVetoCluster0p4_Me1112Veto+= cscClusterSize[i];
      // if (MuonSystem->cscClusterJetVetoPt[MuonSystem->nCscClusters] < JET_PT_CUT && MuonSystem->cscClusterMuonVetoPt[MuonSystem->nCscClusters] < MUON_PT_CUT && me1112_veto) MuonSystem->nCsc_JetMuonVetoCluster0p4_Me1112Veto+= cscClusterSize[i];
      float min_deltaR = 15.;
      int index = 999;
      for(int j = 0; j < 2;j++)
      {

        double current_delta_r = RazorAnalyzer::deltaR(cscClusterEta[i], cscClusterPhi[i], gLLP_eta[j], gLLP_phi[j]);
        if (current_delta_r < min_deltaR)
        {
          min_deltaR = current_delta_r;
          index = j;
        }
      }
      if (min_deltaR < 0.4)
      {
        MuonSystem->cscCluster_match_gLLP[i] = true;
        MuonSystem->cscCluster_match_gLLP_minDeltaR[i] = min_deltaR;
        MuonSystem->cscCluster_match_gLLP_index[i] = index;
      }

      //match cluster to met
      MuonSystem->cscClusterMet_dPhi[MuonSystem->nCscClusters] =  RazorAnalyzer::deltaPhi(cscClusterPhi[i],MuonSystem->metPhi);




      MuonSystem->nCscClusters++;
    }*/


    /*if (isData)
    {
      for(int i = 0; i < nCscSegClusters; i++)
      {
        MuonSystem->cscSegClusterX[MuonSystem->nCscSegClusters] =cscSegClusterX[i];
        MuonSystem->cscSegClusterY[MuonSystem->nCscSegClusters] =cscSegClusterY[i];
        MuonSystem->cscSegClusterZ[MuonSystem->nCscSegClusters] =cscSegClusterZ[i];
        MuonSystem->cscSegClusterTime[MuonSystem->nCscSegClusters] =cscSegClusterTime[i];
        MuonSystem->cscSegClusterEta[MuonSystem->nCscSegClusters] =cscSegClusterEta[i];
        MuonSystem->cscSegClusterPhi[MuonSystem->nCscSegClusters] =cscSegClusterPhi[i];
        MuonSystem->cscSegClusterMajorAxis[MuonSystem->nCscSegClusters] =cscSegClusterMajorAxis[i];
        MuonSystem->cscSegClusterMinorAxis[MuonSystem->nCscSegClusters] =cscSegClusterMinorAxis[i];
        MuonSystem->cscSegClusterXSpread[MuonSystem->nCscSegClusters] =cscSegClusterXSpread[i];
        MuonSystem->cscSegClusterYSpread[MuonSystem->nCscSegClusters] =cscSegClusterYSpread[i];
        MuonSystem->cscSegClusterZSpread[MuonSystem->nCscSegClusters] =cscSegClusterZSpread[i];
        MuonSystem->cscSegClusterEtaPhiSpread[MuonSystem->nCscSegClusters] =cscSegClusterEtaPhiSpread[i];
        MuonSystem->cscSegClusterEtaSpread[MuonSystem->nCscSegClusters] =cscSegClusterEtaSpread[i];
        MuonSystem->cscSegClusterPhiSpread[MuonSystem->nCscSegClusters] = cscSegClusterPhiSpread[i];
        MuonSystem->cscSegClusterTimeSpread[MuonSystem->nCscSegClusters] =cscSegClusterTimeSpread[i];

        MuonSystem->cscSegClusterJetVetoPt[MuonSystem->nCscSegClusters] =cscSegClusterJetVetoPt[i];
        MuonSystem->cscSegClusterJetVetoE[MuonSystem->nCscSegClusters] =cscSegClusterJetVetoE[i];
        MuonSystem->cscSegClusterMuonVetoPt[MuonSystem->nCscSegClusters] =cscSegClusterMuonVetoPt[i];
        MuonSystem->cscSegClusterMuonVetoE[MuonSystem->nCscSegClusters] =cscSegClusterMuonVetoE[i];

        MuonSystem->cscSegClusterSize[MuonSystem->nCscSegClusters] = cscSegClusterSize[i];
        MuonSystem->cscSegClusterMaxChamber[MuonSystem->nCscSegClusters] =cscSegClusterMaxChamber[i];
        MuonSystem->cscSegClusterMaxChamberRatio[MuonSystem->nCscSegClusters] = cscSegClusterMaxChamberRatio[i];
        MuonSystem->cscSegClusterNChamber[MuonSystem->nCscSegClusters] =cscSegClusterNChamber[i];
        MuonSystem->cscSegClusterMaxStation[MuonSystem->nCscSegClusters] = cscSegClusterMaxStation[i];
        MuonSystem->cscSegClusterMaxStationRatio[MuonSystem->nCscSegClusters] = cscSegClusterMaxStationRatio[i];
        MuonSystem->cscSegClusterNStation[MuonSystem->nCscSegClusters] = cscSegClusterNStation[i];
        MuonSystem->cscSegClusterMe11Ratio[MuonSystem->nCscSegClusters] =cscSegClusterMe11Ratio[i];
        MuonSystem->cscSegClusterMe12Ratio[MuonSystem->nCscSegClusters] = cscSegClusterMe12Ratio[i];

        MuonSystem->cscSegClusterNSegmentChamberPlus11[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberPlus11[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus12[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberPlus12[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus13[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberPlus13[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus21[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberPlus21[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus22[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberPlus22[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus31[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberPlus31[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus32[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberPlus32[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus41[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberPlus41[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus42[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberPlus42[i];

        MuonSystem->cscSegClusterNSegmentChamberMinus11[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberMinus11[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus12[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberMinus12[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus13[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberMinus13[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus21[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberMinus21[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus22[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberMinus22[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus31[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberMinus31[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus32[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberMinus32[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus41[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberMinus41[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus42[MuonSystem->nCscSegClusters] = cscSegClusterNSegmentChamberMinus42[i];
        bool me1112_veto = MuonSystem->cscSegClusterMe11Ratio[MuonSystem->nCscSegClusters] == 0.0 && MuonSystem->cscSegClusterMe12Ratio[MuonSystem->nCscSegClusters] == 0.0;
        if (MuonSystem->cscSegClusterJetVetoPt[MuonSystem->nCscSegClusters] < JET_PT_CUT) MuonSystem->nCsc_JetVetoCluster0p4 += cscSegClusterSize[i];
        if (MuonSystem->cscSegClusterJetVetoPt[MuonSystem->nCscSegClusters] < JET_PT_CUT && MuonSystem->cscSegClusterMuonVetoPt[MuonSystem->nCscSegClusters] < MUON_PT_CUT) MuonSystem->nCsc_JetMuonVetoCluster0p4 += cscSegClusterSize[i];
        if (MuonSystem->cscSegClusterJetVetoPt[MuonSystem->nCscSegClusters] < JET_PT_CUT && me1112_veto) MuonSystem->nCsc_JetVetoCluster0p4_Me1112Veto+= cscSegClusterSize[i];
        if (MuonSystem->cscSegClusterJetVetoPt[MuonSystem->nCscSegClusters] < JET_PT_CUT && MuonSystem->cscSegClusterMuonVetoPt[MuonSystem->nCscSegClusters] < MUON_PT_CUT && me1112_veto) MuonSystem->nCsc_JetMuonVetoCluster0p4_Me1112Veto+= cscSegClusterSize[i];

        MuonSystem->cscSegClusterMet_dPhi[MuonSystem->nCscSegClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscSegClusterPhi[MuonSystem->nCscSegClusters],MuonSystem->metPhi);

        MuonSystem->nCscSegClusters++;
      }
    }
    else
    {
      for(int i = 0; i < nCscClusters; i++)
      {
        MuonSystem->cscSegClusterX[MuonSystem->nCscSegClusters] =cscClusterX[i];
        MuonSystem->cscSegClusterY[MuonSystem->nCscSegClusters] =cscClusterY[i];
        MuonSystem->cscSegClusterZ[MuonSystem->nCscSegClusters] =cscClusterZ[i];
        MuonSystem->cscSegClusterTime[MuonSystem->nCscSegClusters] =cscClusterTime[i];
        MuonSystem->cscSegClusterEta[MuonSystem->nCscSegClusters] =cscClusterEta[i];
        MuonSystem->cscSegClusterPhi[MuonSystem->nCscSegClusters] =cscClusterPhi[i];
        MuonSystem->cscSegClusterMajorAxis[MuonSystem->nCscSegClusters] =cscClusterMajorAxis[i];
        MuonSystem->cscSegClusterMinorAxis[MuonSystem->nCscSegClusters] =cscClusterMinorAxis[i];
        MuonSystem->cscSegClusterXSpread[MuonSystem->nCscSegClusters] =cscClusterXSpread[i];
        MuonSystem->cscSegClusterYSpread[MuonSystem->nCscSegClusters] =cscClusterYSpread[i];
        MuonSystem->cscSegClusterZSpread[MuonSystem->nCscSegClusters] =cscClusterZSpread[i];
        MuonSystem->cscSegClusterEtaPhiSpread[MuonSystem->nCscSegClusters] =cscClusterEtaPhiSpread[i];
        MuonSystem->cscSegClusterEtaSpread[MuonSystem->nCscSegClusters] =cscClusterEtaSpread[i];
        MuonSystem->cscSegClusterPhiSpread[MuonSystem->nCscSegClusters] = cscClusterPhiSpread[i];
        MuonSystem->cscSegClusterTimeSpread[MuonSystem->nCscSegClusters] =cscClusterTimeSpread[i];

        MuonSystem->cscSegClusterJetVetoPt[MuonSystem->nCscSegClusters] =cscClusterJetVetoPt[i];
        MuonSystem->cscSegClusterJetVetoE[MuonSystem->nCscSegClusters] =cscClusterJetVetoE[i];
        MuonSystem->cscSegClusterMuonVetoPt[MuonSystem->nCscSegClusters] =cscClusterMuonVetoPt[i];
        MuonSystem->cscSegClusterMuonVetoE[MuonSystem->nCscSegClusters] =cscClusterMuonVetoE[i];

        MuonSystem->cscSegClusterSize[MuonSystem->nCscSegClusters] = cscClusterSize[i];
        MuonSystem->cscSegClusterMaxChamber[MuonSystem->nCscSegClusters] =cscClusterMaxChamber[i];
        MuonSystem->cscSegClusterMaxChamberRatio[MuonSystem->nCscSegClusters] = cscClusterMaxChamberRatio[i];
        MuonSystem->cscSegClusterNChamber[MuonSystem->nCscSegClusters] =cscClusterNChamber[i];
        MuonSystem->cscSegClusterMaxStation[MuonSystem->nCscSegClusters] = cscClusterMaxStation[i];
        MuonSystem->cscSegClusterMaxStationRatio[MuonSystem->nCscSegClusters] = cscClusterMaxStationRatio[i];
        MuonSystem->cscSegClusterNStation[MuonSystem->nCscSegClusters] = cscClusterNStation[i];
        MuonSystem->cscSegClusterMe11Ratio[MuonSystem->nCscSegClusters] =cscClusterMe11Ratio[i];
        MuonSystem->cscSegClusterMe12Ratio[MuonSystem->nCscSegClusters] = cscClusterMe12Ratio[i];

        MuonSystem->cscSegClusterNSegmentChamberPlus11[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberPlus11[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus12[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberPlus12[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus13[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberPlus13[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus21[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberPlus21[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus22[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberPlus22[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus31[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberPlus31[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus32[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberPlus32[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus41[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberPlus41[i];
        MuonSystem->cscSegClusterNSegmentChamberPlus42[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberPlus42[i];

        MuonSystem->cscSegClusterNSegmentChamberMinus11[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberMinus11[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus12[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberMinus12[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus13[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberMinus13[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus21[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberMinus21[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus22[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberMinus22[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus31[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberMinus31[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus32[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberMinus32[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus41[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberMinus41[i];
        MuonSystem->cscSegClusterNSegmentChamberMinus42[MuonSystem->nCscSegClusters] = cscClusterNSegmentChamberMinus42[i];

        MuonSystem->cscSegClusterMet_dPhi[MuonSystem->nCscSegClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscSegClusterPhi[MuonSystem->nCscSegClusters],MuonSystem->metPhi);
        if (!isData) {
          //match to genparticles
          float min_deltaR = 15.;
          int index = 999;
          for(int j = 0; j < nGenParticle; j++)
          {

              if (abs(gParticleId[j]) >= 100 && abs(gParticleId[j]) <=350) continue;

              double current_delta_r = deltaR(cscClusterEta[i], cscClusterPhi[i], gParticleEta[j], gParticlePhi[j]);
              // cout<<"min_deltaR: "<<min_deltaR<<endl;

              if (current_delta_r < min_deltaR)
              {
                min_deltaR = current_delta_r;
                index = j;
              }
          }
          if (min_deltaR < 0.4)
          {
              MuonSystem->cscSegCluster_match_gParticle_minDeltaR[MuonSystem->nCscSegClusters] = min_deltaR;
              MuonSystem->cscSegCluster_match_gParticle_index[MuonSystem->nCscSegClusters] = index;
              MuonSystem->cscSegCluster_match_gParticle_id[MuonSystem->nCscSegClusters] = gParticleId[index];
          }
        }
        MuonSystem->nCscSegClusters++;
      }
    }*/

    /*
   for(int i = 0; i < nCscRechitClusters; i++)
    {
      // if (cscRechitClusterTime[i] > -12.5) continue;
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

      MuonSystem->cscRechitCluster_match_gParticle_id[MuonSystem->nCscRechitClusters] = cscRechitCluster_match_gParticle_id[i];

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
      }
      // else{
      //   continue;
      // }
      MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhi);

      MuonSystem->nCscRechitClusters++;
    }
    */

   //----------------------------
   // RPCs
   //----------------------------

   for(int i = 0; i < nRpc; i++)
   {
     MuonSystem->rpcX[MuonSystem->nRpc] = rpcX[i];
     MuonSystem->rpcY[MuonSystem->nRpc] = rpcY[i];
     MuonSystem->rpcZ[MuonSystem->nRpc] = rpcZ[i];
     MuonSystem->rpcEta[MuonSystem->nRpc] = rpcEta[i];
     MuonSystem->rpcPhi[MuonSystem->nRpc] = rpcPhi[i];
     MuonSystem->rpcBx[MuonSystem->nRpc] = rpcBx[i];

     MuonSystem->nRpc++;
   }

   //----------------------------
   // HO
   //----------------------------

   for(int i = 0; i < nHORechits; i++)
   {
     MuonSystem->hoRechit_X[MuonSystem->nHORechits] = hoRechit_X[i];
     MuonSystem->hoRechit_Y[MuonSystem->nHORechits] = hoRechit_Y[i];
     MuonSystem->hoRechit_Z[MuonSystem->nHORechits] = hoRechit_Z[i];
     MuonSystem->hoRechit_Eta[MuonSystem->nHORechits] = hoRechit_Eta[i];
     MuonSystem->hoRechit_Phi[MuonSystem->nHORechits] = hoRechit_Phi[i];
     MuonSystem->hoRechit_T[MuonSystem->nHORechits] = hoRechit_T[i];
     MuonSystem->hoRechit_E[MuonSystem->nHORechits] = hoRechit_E[i];

     MuonSystem->nHORechits++;
   }

    //-----------------------------
    // DT INFO
    //-----------------------------
    // if( nDt < 10 ) continue;
   //if(MuonSystem->category !=1) continue;
   //if(nDtRechitClusters == 0) continue;
   for(int i = 0; i < nDtRechits; i++)
   {
     MuonSystem->dtRechitX[MuonSystem->nDtRechits] =dtRechitX[i];
     MuonSystem->dtRechitY[MuonSystem->nDtRechits] =dtRechitY[i];
     MuonSystem->dtRechitZ[MuonSystem->nDtRechits] =dtRechitZ[i];
     MuonSystem->dtRechitEta[MuonSystem->nDtRechits] =dtRechitEta[i];
     MuonSystem->dtRechitPhi[MuonSystem->nDtRechits] =dtRechitPhi[i];
     MuonSystem->dtRechitStation[MuonSystem->nDtRechits] =dtRechitStation[i];
     MuonSystem->dtRechitWheel[MuonSystem->nDtRechits] =dtRechitWheel[i];

     MuonSystem->nDtRechits++;
   }

    for(int i = 0; i < nDtRechitClusters; i++)
    {
      MuonSystem->dtRechitClusterX[MuonSystem->nDtRechitClusters] =dtRechitClusterX[i];
      MuonSystem->dtRechitClusterY[MuonSystem->nDtRechitClusters] =dtRechitClusterY[i];
      MuonSystem->dtRechitClusterZ[MuonSystem->nDtRechitClusters] =dtRechitClusterZ[i];
      MuonSystem->dtRechitClusterTime[MuonSystem->nDtRechitClusters] =dtRechitClusterTime[i];
      MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters] =dtRechitClusterEta[i];
      MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters] =dtRechitClusterPhi[i];
      MuonSystem->dtRechitClusterMajorAxis[MuonSystem->nDtRechitClusters] =dtRechitClusterMajorAxis[i];
      MuonSystem->dtRechitClusterMinorAxis[MuonSystem->nDtRechitClusters] =dtRechitClusterMinorAxis[i];
      MuonSystem->dtRechitClusterXSpread[MuonSystem->nDtRechitClusters] =dtRechitClusterXSpread[i];
      MuonSystem->dtRechitClusterYSpread[MuonSystem->nDtRechitClusters] =dtRechitClusterYSpread[i];
      MuonSystem->dtRechitClusterZSpread[MuonSystem->nDtRechitClusters] =dtRechitClusterZSpread[i];
      MuonSystem->dtRechitClusterEtaPhiSpread[MuonSystem->nDtRechitClusters] =dtRechitClusterEtaPhiSpread[i];
      MuonSystem->dtRechitClusterEtaSpread[MuonSystem->nDtRechitClusters] =dtRechitClusterEtaSpread[i];
      MuonSystem->dtRechitClusterPhiSpread[MuonSystem->nDtRechitClusters] = dtRechitClusterPhiSpread[i];
      MuonSystem->dtRechitClusterTimeSpread[MuonSystem->nDtRechitClusters] =dtRechitClusterTimeSpread[i];

      MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] =dtRechitClusterJetVetoPt[i];
      MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] =dtRechitClusterJetVetoE[i];
      MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] =dtRechitClusterMuonVetoPt[i];
      MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] =dtRechitClusterMuonVetoE[i];

      MuonSystem->dtRechitClusterSize[MuonSystem->nDtRechitClusters] = dtRechitClusterSize[i];
      MuonSystem->dtRechitClusterMaxChamber[MuonSystem->nDtRechitClusters] =dtRechitClusterMaxChamber[i];
      MuonSystem->dtRechitClusterMaxChamberRatio[MuonSystem->nDtRechitClusters] = dtRechitClusterMaxChamberRatio[i];
      MuonSystem->dtRechitClusterNChamber[MuonSystem->nDtRechitClusters] =dtRechitClusterNChamber[i];
      MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters] = dtRechitClusterMaxStation[i];
      MuonSystem->dtRechitClusterMaxStationRatio[MuonSystem->nDtRechitClusters] = dtRechitClusterMaxStationRatio[i];
      MuonSystem->dtRechitClusterNStation[MuonSystem->nDtRechitClusters] = dtRechitClusterNStation[i];
      //MuonSystem->dtRechitClusterMe11Ratio[MuonSystem->nDtRechitClusters] =dtRechitClusterMe11Ratio[i];
      //MuonSystem->dtRechitClusterMe12Ratio[MuonSystem->nDtRechitClusters] = dtRechitClusterMe12Ratio[i];
      MuonSystem->dtRechitClusterNSegmentStation1[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentStation1[i];
      MuonSystem->dtRechitClusterNSegmentStation2[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentStation2[i];
      MuonSystem->dtRechitClusterNSegmentStation3[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentStation3[i];
      MuonSystem->dtRechitClusterNSegmentStation4[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentStation4[i];
      /*MuonSystem->dtRechitClusterNSegmentChamberPlus11[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberPlus11[i];
      MuonSystem->dtRechitClusterNSegmentChamberPlus12[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberPlus12[i];
      MuonSystem->dtRechitClusterNSegmentChamberPlus13[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberPlus13[i];
      MuonSystem->dtRechitClusterNSegmentChamberPlus21[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberPlus21[i];
      MuonSystem->dtRechitClusterNSegmentChamberPlus22[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberPlus22[i];
      MuonSystem->dtRechitClusterNSegmentChamberPlus31[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberPlus31[i];
      MuonSystem->dtRechitClusterNSegmentChamberPlus32[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberPlus32[i];
      MuonSystem->dtRechitClusterNSegmentChamberPlus41[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberPlus41[i];
      MuonSystem->dtRechitClusterNSegmentChamberPlus42[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberPlus42[i];

      MuonSystem->dtRechitClusterNSegmentChamberMinus11[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberMinus11[i];
      MuonSystem->dtRechitClusterNSegmentChamberMinus12[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberMinus12[i];
      MuonSystem->dtRechitClusterNSegmentChamberMinus13[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberMinus13[i];
      MuonSystem->dtRechitClusterNSegmentChamberMinus21[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberMinus21[i];
      MuonSystem->dtRechitClusterNSegmentChamberMinus22[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberMinus22[i];
      MuonSystem->dtRechitClusterNSegmentChamberMinus31[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberMinus31[i];
      MuonSystem->dtRechitClusterNSegmentChamberMinus32[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberMinus32[i];
      MuonSystem->dtRechitClusterNSegmentChamberMinus41[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberMinus41[i];
      MuonSystem->dtRechitClusterNSegmentChamberMinus42[MuonSystem->nDtRechitClusters] = dtRechitClusterNSegmentChamberMinus42[i];
      */
      //bool me1112_veto = MuonSystem->dtRechitClusterMe11Ratio[MuonSystem->nDtRechitClusters] == 0.0 && MuonSystem->dtRechitClusterMe12Ratio[MuonSystem->nDtRechitClusters] == 0.0;
      // if (MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] < JET_PT_CUT) MuonSystem->nDt_JetVetoCluster0p4 += dtRechitClusterSize[i];
      // if (MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] < JET_PT_CUT && MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] < MUON_PT_CUT) MuonSystem->nDt_JetMuonVetoCluster0p4 += dtRechitClusterSize[i];
      // if (MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] < JET_PT_CUT && me1112_veto) MuonSystem->nDt_JetVetoCluster0p4_Me1112Veto+= dtRechitClusterSize[i];
      // if (MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] < JET_PT_CUT && MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] < MUON_PT_CUT && me1112_veto) MuonSystem->nDt_JetMuonVetoCluster0p4_Me1112Veto+= dtRechitClusterSize[i];
      float min_deltaR = 15.;
      int index = 999;
      /*for(int j = 0; j < 2;j++)
      {

        double current_delta_r = RazorAnalyzer::deltaR(dtRechitClusterEta[i], dtRechitClusterPhi[i], gLLP_eta[j], gLLP_phi[j]);
        if (current_delta_r < min_deltaR)
        {
          min_deltaR = current_delta_r;
          index = j;
        }
      }
      if (min_deltaR < 0.4)
      {
        MuonSystem->dtRechitCluster_match_gLLP[i] = true;
        MuonSystem->dtRechitCluster_match_gLLP_minDeltaR[i] = min_deltaR;
        MuonSystem->dtRechitCluster_match_gLLP_index[i] = index;
	}*/

      MuonSystem->nDtRechitClusters++;
    }

    if(!isData && signalScan)
    {
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
    // outFile->Write();
    outFile->Close();
  }
}
