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

  bool signalScan = int(options/10) == 1;
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
  // TH1F *NEvents_genweight = new TH1F("NEvents_genweight", "NEvents_genweight", 1, 1, 2);

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


    // cout<<*lheComments<<endl;
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
        if (abs(gParticleId[i])== 25 || abs(gParticleId[i] == 35))
        {
	        MuonSystem->gHiggsPt = gParticlePt[i];
          MuonSystem->gHiggsEta = gParticleEta[i];
          MuonSystem->gHiggsPhi = gParticlePhi[i];
          MuonSystem->gHiggsE = gParticleE[i];

        }
        // MuonSystem->gParticleStatus[MuonSystem->nGenParticle] = gParticleStatus[i];
        // MuonSystem->gParticleId[MuonSystem->nGenParticle]  = gParticleId[i];
        // MuonSystem->gParticleMotherId[MuonSystem->nGenParticle]  = gParticleMotherId[i];
        // MuonSystem->gParticlePt[MuonSystem->nGenParticle]  = gParticlePt[i];
        // MuonSystem->gParticleEta[MuonSystem->nGenParticle]  = gParticleEta[i];
        // MuonSystem->gParticlePhi[MuonSystem->nGenParticle]  = gParticlePhi[i];
        // MuonSystem->gParticleE[MuonSystem->nGenParticle]  = gParticleE[i];
        // MuonSystem->nGenParticle++;
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
      /*for(int i=0; i < nGenJets; i++)
      {
        // cout<<genJetE[i]<<","<<MuonSystem->nGenJets<<","<<nGenJets<<endl;
        MuonSystem->genJetE[MuonSystem->nGenJets] = genJetE[i];
        MuonSystem->genJetPt[MuonSystem->nGenJets] = genJetPt[i];
        MuonSystem->genJetEta[MuonSystem->nGenJets] = genJetEta[i];
        MuonSystem->genJetPhi[MuonSystem->nGenJets] = genJetPhi[i];
        // MuonSystem->genJetMET[MuonSystem->nGenJets] = genJetMET[i];
        // MuonSystem->nGenJets++;
      }*/
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
    std::pair<double,double> corrected_met;
    if (analysisTag=="Razor2016_07Aug2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2016, !isData, nPV);
    else if (analysisTag=="Razor2017_17Nov2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2017, !isData, nPV);
    else if (analysisTag=="Razor2018_17SeptEarlyReReco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2018, !isData, nPV);

    MuonSystem->metXYCorr = corrected_met.first;
    MuonSystem->metPhiXYCorr = corrected_met.second;


    // if (MuonSystem->met < 200) continue;
    // if (nCscRechitClusters==0) continue;


    //Triggers
    for(int i = 0; i < NTriggersMAX; i++){
      MuonSystem->HLTDecision[i] = HLTDecision[i];
    }
    // flags
    MuonSystem->Flag_HBHENoiseFilter = Flag_HBHENoiseFilter;
    MuonSystem->Flag_HBHEIsoNoiseFilter = Flag_HBHEIsoNoiseFilter;
    MuonSystem->Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
    // MuonSystem->Flag_CSCTightHaloFilter = Flag_CSCTightHaloFilter;
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

    // if (!MuonSystem->Flag2_all) continue;
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
      float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];

      tmpMuon.passLooseIso = muonIso<0.25;
      tmpMuon.passTightIso = muonIso<0.15;
      tmpMuon.passVTightIso = muonIso<0.10;
      tmpMuon.passVVTightIso = muonIso<0.05;


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

    if (foundZ  && Leptons.size() == 2 && leadingLepPt > zh_lepton0_cut)
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
      MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
      MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
      MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
      MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;


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
    // JEC = 1.0;
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
        if (tmp.jet.Pt()>30 && abs(tmp.jet.Eta())<2.4)
        {
          jetMet_dPhiMin_temp = abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi));

        }
      }
      // if (tmp.jet.Pt()>20 )MuonSystem->HT = MuonSystem->HT + tmp.jet.Pt();
      MuonSystem->jetTightPassId[MuonSystem->nJets] = tmp.passId;
      MuonSystem->jetChargedEMEnergyFraction[MuonSystem->nJets] = tmp.jetChargedEMEnergyFraction;
      MuonSystem->jetNeutralEMEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralEMEnergyFraction;
      MuonSystem->jetChargedHadronEnergyFraction[MuonSystem->nJets] = tmp.jetChargedHadronEnergyFraction;
      MuonSystem->jetNeutralHadronEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralHadronEnergyFraction;
      MuonSystem->jetPassMuFrac[MuonSystem->nJets] = tmp.jetPassMuFrac;

      float min_deltaR = 15.;
      int index = 999;
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
    // cout<<"number of clusters: "<<nCscRechitClusters<<endl;
   /*for(int i = 0; i < nCscRechitClusters; i++)
    {
      // if (cscRechitClusterTime[i] > -12.5) continue;
      // if(cscRechitClusterNRechitChamberPlus11[i] != 0)continue;
      // if(cscRechitClusterNRechitChamberPlus12[i] != 0)continue;
      // if(cscRechitClusterNRechitChamberMinus11[i] != 0)continue;
      // if(cscRechitClusterNRechitChamberMinus12[i] != 0)continue;

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

        MuonSystem->cscRechitCluster_match_gLLP_eta[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_eta[index];
        MuonSystem->cscRechitCluster_match_gLLP_phi[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_phi[index];
        MuonSystem->cscRechitCluster_match_gLLP_decay_r[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_r[index];
        MuonSystem->cscRechitCluster_match_gLLP_decay_x[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_x[index];
        MuonSystem->cscRechitCluster_match_gLLP_decay_y[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_y[index];
        MuonSystem->cscRechitCluster_match_gLLP_decay_z[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
        MuonSystem->cscRechitCluster_match_gLLP_ctau[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_ctau[index];
        MuonSystem->cscRechitCluster_match_gLLP_beta[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_beta[index];
        MuonSystem->cscRechitCluster_match_gLLP_csc[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_csc[index];
      }

      MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhi);

      MuonSystem->nCscRechitClusters++;
    }
    */
    // count dt rechits
    MuonSystem->nDTRechits  = 0;
    for (int i = 0; i < nDtRechits; i++) {

      if (dtRechitY[i]>=0.0) MuonSystem->nDTPositiveYRechits++;
      else MuonSystem->nDTNegativeYRechits++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -2) MuonSystem->nDTRechitsChamberMinus12++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -1) MuonSystem->nDTRechitsChamberMinus11++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 0) MuonSystem->nDTRechitsChamber10++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 1) MuonSystem->nDTRechitsChamberPlus11++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 2) MuonSystem->nDTRechitsChamberPlus12++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -2) MuonSystem->nDTRechitsChamberMinus22++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -1) MuonSystem->nDTRechitsChamberMinus21++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 0) MuonSystem->nDTRechitsChamber20++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 1) MuonSystem->nDTRechitsChamberPlus21++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 2) MuonSystem->nDTRechitsChamberPlus22++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -2) MuonSystem->nDTRechitsChamberMinus32++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -1) MuonSystem->nDTRechitsChamberMinus31++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 0) MuonSystem->nDTRechitsChamber30++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 1) MuonSystem->nDTRechitsChamberPlus31++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 2) MuonSystem->nDTRechitsChamberPlus32++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -2) MuonSystem->nDTRechitsChamberMinus42++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -1) MuonSystem->nDTRechitsChamberMinus41++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 0) MuonSystem->nDTRechitsChamber40++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 1) MuonSystem->nDTRechitsChamberPlus41++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 2) MuonSystem->nDTRechitsChamberPlus42++;


      MuonSystem->nDTRechits++;
    }
    if ( MuonSystem->nDTRechitsChamberMinus12 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberMinus11 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamber10 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberPlus11 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberPlus12 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberMinus22 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberMinus21 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamber20 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberPlus21 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberPlus22 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberMinus32 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberMinus31 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamber30 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberPlus31 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberPlus32 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberMinus42 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberMinus41 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamber40 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberPlus41 > 50) MuonSystem->nDtRings++;
    if ( MuonSystem->nDTRechitsChamberPlus42 > 50) MuonSystem->nDtRings++;







    // cout << "Number of rec hits: "<<ncscRechits<<endl;
    vector<Point> points;
    vector<int> cscRechitsClusterId;
    points.clear();
    MuonSystem->nCscRechits  = 0;

    for (int i = 0; i < ncscRechits; i++) {
      // if (cscRechitsQuality[i]!=1) continue;
      // if (cscRechitsQuality[i]>2) cout<<cscRechitsQuality[i]<<endl;
      // MuonSystem->cscRechitsPhi[MuonSystem->nCscRechits]           = cscRechitsPhi[i];   //[nCsc]
      // MuonSystem->cscRechitsEta[MuonSystem->nCscRechits]           = cscRechitsEta[i];   //[nCsc]
      // MuonSystem->cscRechitsX[MuonSystem->nCscRechits]             = cscRechitsX[i];   //[nCsc]
      // MuonSystem->cscRechitsY[MuonSystem->nCscRechits]             = cscRechitsY[i];   //[nCsc]
      // MuonSystem->cscRechitsZ[MuonSystem->nCscRechits]             = cscRechitsZ[i];   //[nCsc]
      // MuonSystem->cscRechitsTpeak[MuonSystem->nCscRechits] = cscRechitsTpeak[i];
      // MuonSystem->cscRechitsTwire[MuonSystem->nCscRechits] = cscRechitsTwire[i];
      // MuonSystem->cscRechitsQuality[MuonSystem->nCscRechits] = cscRechitsQuality[i];
      // MuonSystem->cscRechitsStation[MuonSystem->nCscRechits] = cscRechitsStation[i];
      // MuonSystem->cscRechitsChamber[MuonSystem->nCscRechits] = cscRechitsChamber[i];
      // static int station(int index) { return ((index >> START_STATION) & MASK_STATION); }
      // cout<<cscRechitsStation[i]<<", " << ((cscRechitsDetId[i] >> 12) & 07)<<","<<((cscRechitsDetId[i] >> 3) & 077)<<endl;

      // MuonSystem->cscRechitsChamber[MuonSystem->nCscRechits] = ((index >> START_CHAMBER) & MASK_CHAMBER);
      //pick out the right bits for chamber
      int chamber = ((cscRechitsDetId[i] >> 3) & 077); //https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonDetId/interface/CSCDetId.h#L147
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
      cscRechitsClusterId.push_back(-1);
      if (cscRechitsY[i]>=0.0)
      {
        MuonSystem->nCscPositiveYRechits++;
        MuonSystem->cscPosTpeak = MuonSystem->cscPosTpeak + cscRechitsTpeak[i];
      }
      else
      {
        MuonSystem->nCscNegativeYRechits++;
        MuonSystem->cscNegTpeak = MuonSystem->cscNegTpeak + cscRechitsTpeak[i];
      }
      if (cscRechitsTpeak[i]<-12.5)MuonSystem->nEarlyCscRechits++;
      if (cscRechitsTpeak[i]>12.5)MuonSystem->nLateCscRechits++;
      if (cscRechitsTpeak[i]<-25)MuonSystem->nEarly2CscRechits++;
      if (cscRechitsTpeak[i]>25)MuonSystem->nLate2CscRechits++;
      if (cscRechitsChamber[i] == 11) MuonSystem->nCscRechitsChamberPlus11++;
      if (cscRechitsChamber[i] == 12) MuonSystem->nCscRechitsChamberPlus12++;
      if (cscRechitsChamber[i] == 13) MuonSystem->nCscRechitsChamberPlus13++;
      if (cscRechitsChamber[i] == 21) MuonSystem->nCscRechitsChamberPlus21++;
      if (cscRechitsChamber[i] == 22) MuonSystem->nCscRechitsChamberPlus22++;
      if (cscRechitsChamber[i] == 31) MuonSystem->nCscRechitsChamberPlus31++;
      if (cscRechitsChamber[i] == 32) MuonSystem->nCscRechitsChamberPlus32++;
      if (cscRechitsChamber[i] == 41) MuonSystem->nCscRechitsChamberPlus41++;
      if (cscRechitsChamber[i] == 42) MuonSystem->nCscRechitsChamberPlus42++;
      if (cscRechitsChamber[i] == -11) MuonSystem->nCscRechitsChamberMinus11++;
      if (cscRechitsChamber[i] == -12) MuonSystem->nCscRechitsChamberMinus12++;
      if (cscRechitsChamber[i] == -13) MuonSystem->nCscRechitsChamberMinus13++;
      if (cscRechitsChamber[i] == -21) MuonSystem->nCscRechitsChamberMinus21++;
      if (cscRechitsChamber[i] == -22) MuonSystem->nCscRechitsChamberMinus22++;
      if (cscRechitsChamber[i] == -31) MuonSystem->nCscRechitsChamberMinus31++;
      if (cscRechitsChamber[i] == -32) MuonSystem->nCscRechitsChamberMinus32++;
      if (cscRechitsChamber[i] == -41) MuonSystem->nCscRechitsChamberMinus41++;
      if (cscRechitsChamber[i] == -42) MuonSystem->nCscRechitsChamberMinus42++;
      // MuonSystem->nCscRechits++;
    }
    MuonSystem->nCscRings = 0;
    if ( MuonSystem->nCscRechitsChamberPlus11 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus12 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus13 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus21 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus22 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus31 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus32 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus41 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberPlus42 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus11 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus12 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus13 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus21 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus22 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus31 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus32 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus41 > 50) MuonSystem->nCscRings++;
    if ( MuonSystem->nCscRechitsChamberMinus42 > 50) MuonSystem->nCscRings++;
    //Do DBSCAN Clustering
    int min_point = 50;  //minimum number of segments to call it a cluster
    float epsilon = 0.2; //cluster radius parameter
    DBSCAN ds(min_point, epsilon, points);
    ds.run();
    ds.result();
    ds.clusterMoments();
    ds.sort_clusters();

    //Save cluster information
    // cout<<"done clustering: "<<ds.CscCluster.size()<<","<<MuonSystem->nCscRechitClusters<<endl;
    MuonSystem->nCscRechitClusters2 = 0;
    for ( auto &tmp : ds.CscCluster ) {
        // if(tmp.nCscSegmentChamberPlus11 != 0)continue;
        // if(tmp.nCscSegmentChamberPlus12 != 0)continue;
        // if(tmp.nCscSegmentChamberMinus11 != 0)continue;
        // if(tmp.nCscSegmentChamberMinus12 != 0)continue;
        // if(abs(tmp.eta)>=2.1) continue;
        MuonSystem->cscRechitCluster2X[MuonSystem->nCscRechitClusters2] =tmp.x;
        MuonSystem->cscRechitCluster2Y[MuonSystem->nCscRechitClusters2] =tmp.y;
        MuonSystem->cscRechitCluster2Z[MuonSystem->nCscRechitClusters2] =tmp.z;
        MuonSystem->cscRechitCluster2Time[MuonSystem->nCscRechitClusters2] = tmp.t;
        MuonSystem->cscRechitCluster2TimeTotal[MuonSystem->nCscRechitClusters2] = tmp.tTotal;
        MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2] =tmp.eta;
        MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2] = tmp.phi;
        MuonSystem->cscRechitCluster2MajorAxis[MuonSystem->nCscRechitClusters2] =tmp.MajorAxis;
        MuonSystem->cscRechitCluster2MinorAxis[MuonSystem->nCscRechitClusters2] =tmp.MinorAxis;
        MuonSystem->cscRechitCluster2XSpread[MuonSystem->nCscRechitClusters2] =tmp.XSpread;
        MuonSystem->cscRechitCluster2YSpread[MuonSystem->nCscRechitClusters2] =tmp.YSpread;
        MuonSystem->cscRechitCluster2ZSpread[MuonSystem->nCscRechitClusters2] =tmp.ZSpread;
        MuonSystem->cscRechitCluster2EtaPhiSpread[MuonSystem->nCscRechitClusters2] =tmp.EtaPhiSpread;
        MuonSystem->cscRechitCluster2XYSpread[MuonSystem->nCscRechitClusters2] =tmp.XYSpread;
        MuonSystem->cscRechitCluster2RSpread[MuonSystem->nCscRechitClusters2] =tmp.RSpread;
        if (!isData)
        {
          // MuonSystem->cscRechitCluster2XYSpread_phi0p5[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_phi0p5;
          // MuonSystem->cscRechitCluster2XSpread_phi0p5[MuonSystem->nCscRechitClusters2] = tmp.XSpread_phi0p5;
          // MuonSystem->cscRechitCluster2YSpread_phi0p5[MuonSystem->nCscRechitClusters2] = tmp.YSpread_phi0p5;
          // MuonSystem->cscRechitCluster2PhiSpread_phi0p5[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_phi0p5;
          // MuonSystem->cscRechitCluster2EtaPhiSpread_phi0p5[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_phi0p5;
          //
          // MuonSystem->cscRechitCluster2XYSpread_phi0p55[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_phi0p55;
          // MuonSystem->cscRechitCluster2XSpread_phi0p55[MuonSystem->nCscRechitClusters2] = tmp.XSpread_phi0p55;
          // MuonSystem->cscRechitCluster2YSpread_phi0p55[MuonSystem->nCscRechitClusters2] = tmp.YSpread_phi0p55;
          // MuonSystem->cscRechitCluster2PhiSpread_phi0p55[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_phi0p55;
          // MuonSystem->cscRechitCluster2EtaPhiSpread_phi0p55[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_phi0p55;
          //
          // MuonSystem->cscRechitCluster2XYSpread_phi0p6[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_phi0p6;
          // MuonSystem->cscRechitCluster2XSpread_phi0p6[MuonSystem->nCscRechitClusters2] = tmp.XSpread_phi0p6;
          // MuonSystem->cscRechitCluster2YSpread_phi0p6[MuonSystem->nCscRechitClusters2] = tmp.YSpread_phi0p6;
          // MuonSystem->cscRechitCluster2PhiSpread_phi0p6[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_phi0p6;
          // MuonSystem->cscRechitCluster2EtaPhiSpread_phi0p6[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_phi0p6;
          //
          // MuonSystem->cscRechitCluster2XYSpread_phi0p65[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_phi0p65;
          // MuonSystem->cscRechitCluster2XSpread_phi0p65[MuonSystem->nCscRechitClusters2] = tmp.XSpread_phi0p65;
          // MuonSystem->cscRechitCluster2YSpread_phi0p65[MuonSystem->nCscRechitClusters2] = tmp.YSpread_phi0p65;
          // MuonSystem->cscRechitCluster2PhiSpread_phi0p65[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_phi0p65;
          // MuonSystem->cscRechitCluster2EtaPhiSpread_phi0p65[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_phi0p65;
          //
          // MuonSystem->cscRechitCluster2XYSpread_phi0p7[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_phi0p7;
          // MuonSystem->cscRechitCluster2XSpread_phi0p7[MuonSystem->nCscRechitClusters2] = tmp.XSpread_phi0p7;
          // MuonSystem->cscRechitCluster2YSpread_phi0p7[MuonSystem->nCscRechitClusters2] = tmp.YSpread_phi0p7;
          // MuonSystem->cscRechitCluster2PhiSpread_phi0p7[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_phi0p7;
          // MuonSystem->cscRechitCluster2EtaPhiSpread_phi0p7[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_phi0p7;
          //
          // MuonSystem->cscRechitCluster2XYSpread_phi0p75[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_phi0p75;
          // MuonSystem->cscRechitCluster2XSpread_phi0p75[MuonSystem->nCscRechitClusters2] = tmp.XSpread_phi0p75;
          // MuonSystem->cscRechitCluster2YSpread_phi0p75[MuonSystem->nCscRechitClusters2] = tmp.YSpread_phi0p75;
          // MuonSystem->cscRechitCluster2PhiSpread_phi0p75[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_phi0p75;
          // MuonSystem->cscRechitCluster2EtaPhiSpread_phi0p75[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_phi0p75;
          //
          // MuonSystem->cscRechitCluster2XYSpread_r1p2[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_r1p2;
          // MuonSystem->cscRechitCluster2XSpread_r1p2[MuonSystem->nCscRechitClusters2] = tmp.XSpread_r1p2;
          // MuonSystem->cscRechitCluster2YSpread_r1p2[MuonSystem->nCscRechitClusters2] = tmp.YSpread_r1p2;
          // MuonSystem->cscRechitCluster2PhiSpread_r1p2[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_r1p2;
          // MuonSystem->cscRechitCluster2EtaSpread_r1p2[MuonSystem->nCscRechitClusters2] = tmp.EtaSpread_r1p2;
          // MuonSystem->cscRechitCluster2EtaPhiSpread_r1p2[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_r1p2;
          // MuonSystem->cscRechitCluster2RSpread_r1p2[MuonSystem->nCscRechitClusters2] = tmp.RSpread_r1p2;
          //
          //
          // MuonSystem->cscRechitCluster2XYSpread_phi0p7_r1p2[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_phi0p7_r1p2;
          // MuonSystem->cscRechitCluster2XSpread_phi0p7_r1p2[MuonSystem->nCscRechitClusters2] = tmp.XSpread_phi0p7_r1p2;
          // MuonSystem->cscRechitCluster2YSpread_phi0p7_r1p2[MuonSystem->nCscRechitClusters2] = tmp.YSpread_phi0p7_r1p2;
          // MuonSystem->cscRechitCluster2PhiSpread_phi0p7_r1p2[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_phi0p7_r1p2;
          // MuonSystem->cscRechitCluster2EtaSpread_phi0p7_r1p2[MuonSystem->nCscRechitClusters2] = tmp.EtaSpread_phi0p7_r1p2;
          // MuonSystem->cscRechitCluster2EtaPhiSpread_phi0p7_r1p2[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_phi0p7_r1p2;
          // MuonSystem->cscRechitCluster2RSpread_phi0p7_r1p2[MuonSystem->nCscRechitClusters2] = tmp.RSpread_phi0p7_r1p2;
          //
          // MuonSystem->cscRechitCluster2XYSpread_phi0p7_r1p3[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_phi0p7_r1p3;
          // MuonSystem->cscRechitCluster2XSpread_phi0p7_r1p3[MuonSystem->nCscRechitClusters2] = tmp.XSpread_phi0p7_r1p3;
          // MuonSystem->cscRechitCluster2YSpread_phi0p7_r1p3[MuonSystem->nCscRechitClusters2] = tmp.YSpread_phi0p7_r1p3;
          // MuonSystem->cscRechitCluster2PhiSpread_phi0p7_r1p3[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_phi0p7_r1p3;
          // MuonSystem->cscRechitCluster2EtaSpread_phi0p7_r1p3[MuonSystem->nCscRechitClusters2] = tmp.EtaSpread_phi0p7_r1p3;
          // MuonSystem->cscRechitCluster2EtaPhiSpread_phi0p7_r1p3[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_phi0p7_r1p3;
          // MuonSystem->cscRechitCluster2RSpread_phi0p7_r1p3[MuonSystem->nCscRechitClusters2] = tmp.RSpread_phi0p7_r1p3;
          //
          // MuonSystem->cscRechitCluster2RSpread_phi0p7_r1p1[MuonSystem->nCscRechitClusters2] = tmp.RSpread_phi0p7_r1p1;
          // MuonSystem->cscRechitCluster2XYSpread_phi0p7_r1p1[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_phi0p7_r1p1;
          // MuonSystem->cscRechitCluster2XSpread_phi0p7_r1p1[MuonSystem->nCscRechitClusters2] = tmp.XSpread_phi0p7_r1p1;
          // MuonSystem->cscRechitCluster2YSpread_phi0p7_r1p1[MuonSystem->nCscRechitClusters2] = tmp.YSpread_phi0p7_r1p1;
          // MuonSystem->cscRechitCluster2PhiSpread_phi0p7_r1p1[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_phi0p7_r1p1;
          // MuonSystem->cscRechitCluster2EtaSpread_phi0p7_r1p1[MuonSystem->nCscRechitClusters2] = tmp.EtaSpread_phi0p7_r1p1;
          // MuonSystem->cscRechitCluster2EtaPhiSpread_phi0p7_r1p1[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_phi0p7_r1p1;
          // MuonSystem->cscRechitCluster2RSpread_phi0p7_r1p1[MuonSystem->nCscRechitClusters2] = tmp.RSpread_phi0p7_r1p1;
          //
          // MuonSystem->cscRechitCluster2RSpread_phi0p7_r1p15[MuonSystem->nCscRechitClusters2] = tmp.RSpread_phi0p7_r1p15;
          // MuonSystem->cscRechitCluster2XYSpread_phi0p7_r1p15[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_phi0p7_r1p15;
          // MuonSystem->cscRechitCluster2XSpread_phi0p7_r1p15[MuonSystem->nCscRechitClusters2] = tmp.XSpread_phi0p7_r1p15;
          // MuonSystem->cscRechitCluster2YSpread_phi0p7_r1p15[MuonSystem->nCscRechitClusters2] = tmp.YSpread_phi0p7_r1p15;
          // MuonSystem->cscRechitCluster2PhiSpread_phi0p7_r1p15[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_phi0p7_r1p15;
          // MuonSystem->cscRechitCluster2EtaSpread_phi0p7_r1p15[MuonSystem->nCscRechitClusters2] = tmp.EtaSpread_phi0p7_r1p15;
          // MuonSystem->cscRechitCluster2EtaPhiSpread_phi0p7_r1p15[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_phi0p7_r1p15;
          // MuonSystem->cscRechitCluster2RSpread_phi0p7_r1p15[MuonSystem->nCscRechitClusters2] = tmp.RSpread_phi0p7_r1p15;
          //
          // MuonSystem->cscRechitCluster2RSpread_phi0p7_r1p25[MuonSystem->nCscRechitClusters2] = tmp.RSpread_phi0p7_r1p25;
          // MuonSystem->cscRechitCluster2XYSpread_phi0p7_r1p25[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_phi0p7_r1p25;
          // MuonSystem->cscRechitCluster2XSpread_phi0p7_r1p25[MuonSystem->nCscRechitClusters2] = tmp.XSpread_phi0p7_r1p25;
          // MuonSystem->cscRechitCluster2YSpread_phi0p7_r1p25[MuonSystem->nCscRechitClusters2] = tmp.YSpread_phi0p7_r1p25;
          // MuonSystem->cscRechitCluster2PhiSpread_phi0p7_r1p25[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_phi0p7_r1p25;
          // MuonSystem->cscRechitCluster2EtaSpread_phi0p7_r1p25[MuonSystem->nCscRechitClusters2] = tmp.EtaSpread_phi0p7_r1p25;
          // MuonSystem->cscRechitCluster2EtaPhiSpread_phi0p7_r1p25[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_phi0p7_r1p25;
          // MuonSystem->cscRechitCluster2RSpread_phi0p7_r1p25[MuonSystem->nCscRechitClusters2] = tmp.RSpread_phi0p7_r1p25;

          for (int j = 0; j < N_phicorr; j++)
          {
            for (int k = 0; k < N_rcorr; k++)
            {
              MuonSystem->cscRechitCluster2RSpread_corr[MuonSystem->nCscRechitClusters2][j][k] = tmp.RSpread_corr[j][k];
              MuonSystem->cscRechitCluster2XYSpread_corr[MuonSystem->nCscRechitClusters2][j][k] = tmp.XYSpread_corr[j][k];
              MuonSystem->cscRechitCluster2XSpread_corr[MuonSystem->nCscRechitClusters2][j][k] = tmp.XSpread_corr[j][k];
              MuonSystem->cscRechitCluster2YSpread_corr[MuonSystem->nCscRechitClusters2][j][k] = tmp.YSpread_corr[j][k];
              MuonSystem->cscRechitCluster2PhiSpread_corr[MuonSystem->nCscRechitClusters2][j][k] = tmp.PhiSpread_corr[j][k];
              MuonSystem->cscRechitCluster2EtaSpread_corr[MuonSystem->nCscRechitClusters2][j][k] = tmp.EtaSpread_corr[j][k];
              MuonSystem->cscRechitCluster2EtaPhiSpread_corr[MuonSystem->nCscRechitClusters2][j][k] = tmp.EtaPhiSpread_corr[j][k];
            }
          }

        }


        MuonSystem->cscRechitCluster2EtaSpread[MuonSystem->nCscRechitClusters2] =tmp.EtaSpread;
        MuonSystem->cscRechitCluster2PhiSpread[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread;
        MuonSystem->cscRechitCluster2TimeSpread[MuonSystem->nCscRechitClusters2] = tmp.TSpread;
        MuonSystem->cscRechitCluster2Size[MuonSystem->nCscRechitClusters2] = tmp.nCscSegments;

        MuonSystem->cscRechitCluster2NRechitChamberPlus11[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberPlus11;
        MuonSystem->cscRechitCluster2NRechitChamberPlus12[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberPlus12;
        MuonSystem->cscRechitCluster2NRechitChamberPlus13[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberPlus13;
        MuonSystem->cscRechitCluster2NRechitChamberPlus21[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberPlus21;
        MuonSystem->cscRechitCluster2NRechitChamberPlus22[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberPlus22;
        MuonSystem->cscRechitCluster2NRechitChamberPlus31[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberPlus31;
        MuonSystem->cscRechitCluster2NRechitChamberPlus32[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberPlus32;
        MuonSystem->cscRechitCluster2NRechitChamberPlus41[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberPlus41;
        MuonSystem->cscRechitCluster2NRechitChamberPlus42[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberPlus42;
        MuonSystem->cscRechitCluster2NRechitChamberMinus11[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberMinus11;
        MuonSystem->cscRechitCluster2NRechitChamberMinus12[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberMinus12;
        MuonSystem->cscRechitCluster2NRechitChamberMinus13[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberMinus13;
        MuonSystem->cscRechitCluster2NRechitChamberMinus21[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberMinus21;
        MuonSystem->cscRechitCluster2NRechitChamberMinus22[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberMinus22;
        MuonSystem->cscRechitCluster2NRechitChamberMinus31[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberMinus31;
        MuonSystem->cscRechitCluster2NRechitChamberMinus32[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberMinus32;
        MuonSystem->cscRechitCluster2NRechitChamberMinus41[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberMinus41;
        MuonSystem->cscRechitCluster2NRechitChamberMinus42[MuonSystem->nCscRechitClusters2] = tmp.nCscSegmentChamberMinus42;
        MuonSystem->cscRechitCluster2MaxChamber[MuonSystem->nCscRechitClusters2] = tmp.maxChamber;
        MuonSystem->cscRechitCluster2MaxChamberRatio[MuonSystem->nCscRechitClusters2] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
        MuonSystem->cscRechitCluster2NChamber[MuonSystem->nCscRechitClusters2] = tmp.nChamber;
        MuonSystem->cscRechitCluster2MaxStation[MuonSystem->nCscRechitClusters2] = tmp.maxStation;
        MuonSystem->cscRechitCluster2MaxStationRatio[MuonSystem->nCscRechitClusters2] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
        MuonSystem->cscRechitCluster2NStation[MuonSystem->nCscRechitClusters2] = tmp.nStation;
        MuonSystem->cscRechitCluster2NStation5[MuonSystem->nCscRechitClusters2] = tmp.nStation5;
        MuonSystem->cscRechitCluster2NStation10perc[MuonSystem->nCscRechitClusters2] = tmp.nStation10perc;
        MuonSystem->cscRechitCluster2AvgStation[MuonSystem->nCscRechitClusters2] = tmp.avgStation;
        MuonSystem->cscRechitCluster2AvgStation5[MuonSystem->nCscRechitClusters2] = tmp.avgStation5;
        MuonSystem->cscRechitCluster2AvgStation10perc[MuonSystem->nCscRechitClusters2] = tmp.avgStation10perc;

        MuonSystem->cscRechitCluster2Me11Ratio[MuonSystem->nCscRechitClusters2] = tmp.Me11Ratio;
        MuonSystem->cscRechitCluster2Me12Ratio[MuonSystem->nCscRechitClusters2] = tmp.Me12Ratio;
       //  for (unsigned int j = 0; j < tmp.segment_id.size(); j++)
       // {
       //   // MuonSystem->cscRechitsCluster2Id[tmp.segment_id[j]] = MuonSystem->nCscRechitClusters2;
       //   cscRechitsClusterId[tmp.segment_id[j]] = MuonSystem->nCscRechitClusters2+1;
       // }
       if(MuonSystem->category == 2)
       {
         if (MuonSystem->nLeptons<2) cout<<"less than 2"<<endl;
         if (RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex1], MuonSystem->lepPhi[MuonSystem->ZleptonIndex1], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4) {
           MuonSystem->cscRechitCluster2ZLep1[MuonSystem->nCscRechitClusters2] = true;
           MuonSystem->cscRechitCluster2ZLep1Id[MuonSystem->nCscRechitClusters2] = MuonSystem->lepPdgId[MuonSystem->ZleptonIndex1];
           MuonSystem->cscRechitCluster2ZLep1LooseIso[MuonSystem->nCscRechitClusters2]  = MuonSystem->lepPassLooseIso[MuonSystem->ZleptonIndex1];
           MuonSystem->cscRechitCluster2ZLep1TightIso[MuonSystem->nCscRechitClusters2]  = MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
           MuonSystem->cscRechitCluster2ZLep1VTightIso[MuonSystem->nCscRechitClusters2]  = MuonSystem->lepPassVTightIso[MuonSystem->ZleptonIndex1];
           MuonSystem->cscRechitCluster2ZLep1VVTightIso[MuonSystem->nCscRechitClusters2]  = MuonSystem->lepPassVVTightIso[MuonSystem->ZleptonIndex1];
           MuonSystem->cscRechitCluster2ZLep1TightId[MuonSystem->nCscRechitClusters2]  = MuonSystem->lepPassId[MuonSystem->ZleptonIndex1];
         }
         if (RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex2], MuonSystem->lepPhi[MuonSystem->ZleptonIndex2], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4) {
           MuonSystem->cscRechitCluster2ZLep2[MuonSystem->nCscRechitClusters2] = true;
           MuonSystem->cscRechitCluster2ZLep2Id[MuonSystem->nCscRechitClusters2] = MuonSystem->lepPdgId[MuonSystem->ZleptonIndex2];
           MuonSystem->cscRechitCluster2ZLep2LooseIso[MuonSystem->nCscRechitClusters2]  = MuonSystem->lepPassLooseIso[MuonSystem->ZleptonIndex1];
           MuonSystem->cscRechitCluster2ZLep2TightIso[MuonSystem->nCscRechitClusters2]  = MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
           MuonSystem->cscRechitCluster2ZLep2VTightIso[MuonSystem->nCscRechitClusters2]  = MuonSystem->lepPassVTightIso[MuonSystem->ZleptonIndex1];
           MuonSystem->cscRechitCluster2ZLep2VVTightIso[MuonSystem->nCscRechitClusters2]  = MuonSystem->lepPassVVTightIso[MuonSystem->ZleptonIndex1];
           MuonSystem->cscRechitCluster2ZLep2TightId[MuonSystem->nCscRechitClusters2]  = MuonSystem->lepPassId[MuonSystem->ZleptonIndex1];
         }
       }

        //Jet veto/ muon veto
        MuonSystem->cscRechitCluster2JetVetoPt[MuonSystem->nCscRechitClusters2] = 0.0;
        MuonSystem->cscRechitCluster2JetVetoE[MuonSystem->nCscRechitClusters2] = 0.0;
        MuonSystem->cscRechitCluster2MuonVetoPt[MuonSystem->nCscRechitClusters2] = 0.0;
        MuonSystem->cscRechitCluster2MuonVetoE[MuonSystem->nCscRechitClusters2] = 0.0;
        MuonSystem->cscRechitCluster2GenMuonVetoPt[MuonSystem->nCscRechitClusters2] = 0.0;
        MuonSystem->cscRechitCluster2GenMuonVetoE[MuonSystem->nCscRechitClusters2] = 0.0;
        MuonSystem->cscRechitCluster2IsoMuonVetoPt[MuonSystem->nCscRechitClusters2] = 0.0;

        // jet veto
        for(int i = 0; i < nJets; i++)
        {
          if (fabs(jetEta[i]>3.0)) continue;
          if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 && jetPt[i] > MuonSystem->cscRechitCluster2JetVetoPt[MuonSystem->nCscRechitClusters2] ) {
            MuonSystem->cscRechitCluster2JetVetoPt[MuonSystem->nCscRechitClusters2]  = jetPt[i];
          }
          if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2], MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 && jetE[i] > MuonSystem->cscRechitCluster2JetVetoE[MuonSystem->nCscRechitClusters2] ) {
            MuonSystem->cscRechitCluster2JetVetoE[MuonSystem->nCscRechitClusters2]  = jetE[i];
          }
        }
        float min_deltaR = 15.;
        int index = 999;
        // for(int i = 0; i < nCscRechitClusters; i++)
        // {
        //   double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2], MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2], cscRechitClusterEta[i], cscRechitClusterPhi[i]);
        //   if (current_delta_r < min_deltaR)
        //   {
        //     min_deltaR = current_delta_r;
        //     index = i;
        //   }
        //
        //
        // }
        // if (min_deltaR < 0.4)
        // {
        //   MuonSystem->cscRechitCluster2JetVetoPt[MuonSystem->nCscRechitClusters2]  = cscRechitClusterJetVetoPt[index];
        //   MuonSystem->cscRechitCluster2JetVetoE[MuonSystem->nCscRechitClusters2]  = cscRechitClusterJetVetoE[index];
        //   MuonSystem->cscRechitCluster2JetVetoPt_0p6[MuonSystem->nCscRechitClusters2]  = cscRechitClusterJetVetoPt[index];
        //   MuonSystem->cscRechitCluster2JetVetoE_0p6[MuonSystem->nCscRechitClusters2]  = cscRechitClusterJetVetoE[index];
        //   MuonSystem->cscRechitCluster2JetVetoPt_0p8[MuonSystem->nCscRechitClusters2]  = cscRechitClusterJetVetoPt[index];
        //   MuonSystem->cscRechitCluster2JetVetoE_0p8[MuonSystem->nCscRechitClusters2]  = cscRechitClusterJetVetoE[index];
        //
        // }
        // for(int i = 0; i < nJets; i++)
        // {
        //   if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.6 && jetPt[i] > MuonSystem->cscRechitCluster2JetVetoPt_0p6[MuonSystem->nCscRechitClusters2] ) {
        //     MuonSystem->cscRechitCluster2JetVetoPt_0p6[MuonSystem->nCscRechitClusters2]  = jetPt[i];
        //     MuonSystem->cscRechitCluster2JetVetoE_0p6[MuonSystem->nCscRechitClusters2]  = jetE[i];
        //   }
        //   if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.8 && jetPt[i] > MuonSystem->cscRechitCluster2JetVetoPt_0p8[MuonSystem->nCscRechitClusters2] ) {
        //     MuonSystem->cscRechitCluster2JetVetoPt_0p8[MuonSystem->nCscRechitClusters2]  = jetPt[i];
        //     MuonSystem->cscRechitCluster2JetVetoE_0p8[MuonSystem->nCscRechitClusters2]  = jetE[i];
        //   }
        // }

        for(int i = 0; i < nMuons; i++)
        {
          if (fabs(muonEta[i]>3.0)) continue;
          float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
          if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 && muonPt[i] > MuonSystem->cscRechitCluster2MuonVetoPt[MuonSystem->nCscRechitClusters2] ) {
            MuonSystem->cscRechitCluster2MuonVetoPt[MuonSystem->nCscRechitClusters2]  = muonPt[i];
            MuonSystem->cscRechitCluster2MuonVetoE[MuonSystem->nCscRechitClusters2]  = muonE[i];
            MuonSystem->cscRechitCluster2MuonVetoPhi[MuonSystem->nCscRechitClusters2]  = muonPhi[i];
            MuonSystem->cscRechitCluster2MuonVetoEta[MuonSystem->nCscRechitClusters2]  = muonEta[i];
            MuonSystem->cscRechitCluster2MuonVetoEta[MuonSystem->nCscRechitClusters2]  = muonEta[i];
            MuonSystem->cscRechitCluster2MuonVetoLooseIso[MuonSystem->nCscRechitClusters2]  = muonIso<0.25;
            MuonSystem->cscRechitCluster2MuonVetoTightIso[MuonSystem->nCscRechitClusters2]  = muonIso<0.15;
            MuonSystem->cscRechitCluster2MuonVetoVTightIso[MuonSystem->nCscRechitClusters2]  = muonIso<0.10;
            MuonSystem->cscRechitCluster2MuonVetoVVTightIso[MuonSystem->nCscRechitClusters2]  = muonIso<0.05;
            MuonSystem->cscRechitCluster2MuonVetoTightId[MuonSystem->nCscRechitClusters2]  = isMuonPOGTightMuon(i);


          }


        }

        // match to gen-level muon
        if(!isData)
        {
          for(int i = 0; i < nGenParticle; i++)
          {
            if (abs(gParticleId[i])!=13) continue;
            if (abs(gParticleMotherId[i])>24 || abs(gParticleMotherId[i])<23) continue;
            if (abs(gParticleStatus[i])!=1)continue;
            // if (fabs(muonEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 ) {
              MuonSystem->cscRechitCluster2GenMuonVetoPt[MuonSystem->nCscRechitClusters2]  = gParticlePt[i];
              MuonSystem->cscRechitCluster2GenMuonVetoE[MuonSystem->nCscRechitClusters2]  = gParticleE[i];

            }

          }
          min_deltaR = 15.;
          index = 999;
          for(int j = 0; j < nGenParticle; j++)
          {

            double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2], MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2], gParticleEta[j], gParticlePhi[j]);

            if (current_delta_r < min_deltaR)
            {
              min_deltaR = current_delta_r;
              index = j;
            }
          }
          if (min_deltaR < 0.4)
          {

            MuonSystem->cscRechitCluster2_match_gParticle[MuonSystem->nCscRechitClusters2] = true;
            MuonSystem->cscRechitCluster2_match_gParticle_minDeltaR[MuonSystem->nCscRechitClusters2] = min_deltaR;
            MuonSystem->cscRechitCluster2_match_gParticle_index[MuonSystem->nCscRechitClusters2] = index;
            MuonSystem->cscRechitCluster2_match_gParticle_id[MuonSystem->nCscRechitClusters2] = gParticleId[index];

            MuonSystem->cscRechitCluster2_match_gParticle_eta[MuonSystem->nCscRechitClusters2] = gParticleEta[index];
            MuonSystem->cscRechitCluster2_match_gParticle_phi[MuonSystem->nCscRechitClusters2] = gParticlePhi[index];
            MuonSystem->cscRechitCluster2_match_gParticle_E[MuonSystem->nCscRechitClusters2] = gParticleE[index];
            MuonSystem->cscRechitCluster2_match_gParticle_pt[MuonSystem->nCscRechitClusters2] = gParticlePt[index];
            MuonSystem->cscRechitCluster2_match_gParticle_MotherId[MuonSystem->nCscRechitClusters2] = gParticleMotherId[index];


          }
          // match to gen level LLP
          min_deltaR = 15.;
          index = 999;
          for(int j = 0; j < 2;j++)
          {

            double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2], MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2], gLLP_eta[j], gLLP_phi[j]);
            if (current_delta_r < min_deltaR)
            {
              min_deltaR = current_delta_r;
              index = j;
            }
          }
          if (min_deltaR < 0.4)
          {
            MuonSystem->cscRechitCluster2_match_gLLP[MuonSystem->nCscRechitClusters2] = true;
            MuonSystem->cscRechitCluster2_match_gLLP_minDeltaR[MuonSystem->nCscRechitClusters2] = min_deltaR;
            MuonSystem->cscRechitCluster2_match_gLLP_index[MuonSystem->nCscRechitClusters2] = index;
            MuonSystem->cscRechitCluster2_match_gLLP_eta[MuonSystem->nCscRechitClusters2] = MuonSystem->gLLP_eta[index];
            MuonSystem->cscRechitCluster2_match_gLLP_phi[MuonSystem->nCscRechitClusters2] = MuonSystem->gLLP_phi[index];
            MuonSystem->cscRechitCluster2_match_gLLP_decay_r[MuonSystem->nCscRechitClusters2] = MuonSystem->gLLP_decay_vertex_r[index];
            MuonSystem->cscRechitCluster2_match_gLLP_decay_x[MuonSystem->nCscRechitClusters2] = MuonSystem->gLLP_decay_vertex_x[index];
            MuonSystem->cscRechitCluster2_match_gLLP_decay_y[MuonSystem->nCscRechitClusters2] = MuonSystem->gLLP_decay_vertex_y[index];
            MuonSystem->cscRechitCluster2_match_gLLP_decay_z[MuonSystem->nCscRechitClusters2] = MuonSystem->gLLP_decay_vertex_z[index];
            MuonSystem->cscRechitCluster2_match_gLLP_ctau[MuonSystem->nCscRechitClusters2] = MuonSystem->gLLP_ctau[index];
            MuonSystem->cscRechitCluster2_match_gLLP_beta[MuonSystem->nCscRechitClusters2] = MuonSystem->gLLP_beta[index];
            MuonSystem->cscRechitCluster2_match_gLLP_csc[MuonSystem->nCscRechitClusters2] = MuonSystem->gLLP_csc[index];
          }
        }
        for (int i = 0; i < ncscRechits; i++) {
          if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 )
          {
            MuonSystem->cscRechitCluster2_match_cscRechits_0p4[MuonSystem->nCscRechitClusters2] ++;
          }
          // if(abs(cscRechitsEta[i])>2.1)
          // {
          //   if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 )
          //   {
          //     MuonSystem->cscRechitCluster2_match_highEta_0p4[MuonSystem->nCscRechitClusters2] ++;
          //   }
          //   if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.6 )
          //   {
          //     MuonSystem->cscRechitCluster2_match_highEta_0p6[MuonSystem->nCscRechitClusters2] ++;
          //   }
          //   if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.8 )
          //   {
          //     MuonSystem->cscRechitCluster2_match_highEta_0p8[MuonSystem->nCscRechitClusters2] ++;
          //   }
          // }
          if (!(abs(cscRechitsChamber[i]) == 11 || abs(cscRechitsChamber[i])==12))continue;
          if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 )
          {
            MuonSystem->cscRechitCluster2_match_Me1112_0p4[MuonSystem->nCscRechitClusters2] ++;
            if (abs(cscRechitsChamber[i]) == 11) MuonSystem->cscRechitCluster2_match_Me11_0p4[MuonSystem->nCscRechitClusters2] ++;
            if (abs(cscRechitsChamber[i]) == 12) MuonSystem->cscRechitCluster2_match_Me12_0p4[MuonSystem->nCscRechitClusters2] ++;

          }
          if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.6 )
          {
            MuonSystem->cscRechitCluster2_match_Me1112_0p6[MuonSystem->nCscRechitClusters2] ++;
            if (abs(cscRechitsChamber[i]) == 11) MuonSystem->cscRechitCluster2_match_Me11_0p6[MuonSystem->nCscRechitClusters2] ++;
            if (abs(cscRechitsChamber[i]) == 12) MuonSystem->cscRechitCluster2_match_Me12_0p6[MuonSystem->nCscRechitClusters2] ++;
          }
          if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.8 )
          {
            MuonSystem->cscRechitCluster2_match_Me1112_0p8[MuonSystem->nCscRechitClusters2] ++;
            if (abs(cscRechitsChamber[i]) == 11) MuonSystem->cscRechitCluster2_match_Me11_0p8[MuonSystem->nCscRechitClusters2] ++;
            if (abs(cscRechitsChamber[i]) == 12) MuonSystem->cscRechitCluster2_match_Me12_0p8[MuonSystem->nCscRechitClusters2] ++;
          }

        }
        for (int i = 0; i < nCscSeg; i++) {
          if (RazorAnalyzer::deltaR(cscSegEta[i], cscSegPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 )
          {
            MuonSystem->cscRechitCluster2_match_cscSeg_0p4[MuonSystem->nCscRechitClusters2] ++;
            if (abs(cscSegChamber[i]) == 11) MuonSystem->cscRechitCluster2_match_ME11Seg_0p4[MuonSystem->nCscRechitClusters2] ++;
            if (abs(cscSegChamber[i]) == 12) MuonSystem->cscRechitCluster2_match_ME12Seg_0p4[MuonSystem->nCscRechitClusters2] ++;

          }

          if (RazorAnalyzer::deltaR(cscSegEta[i], cscSegPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.6 )
          {
            MuonSystem->cscRechitCluster2_match_cscSeg_0p6[MuonSystem->nCscRechitClusters2] ++;
            if (abs(cscSegChamber[i]) == 11) MuonSystem->cscRechitCluster2_match_ME11Seg_0p6[MuonSystem->nCscRechitClusters2] ++;
            if (abs(cscSegChamber[i]) == 12) MuonSystem->cscRechitCluster2_match_ME12Seg_0p6[MuonSystem->nCscRechitClusters2] ++;
          }
        }
        //match to MB1 DT hits
        for (int i = 0; i < nDtRechits; i++) {
          if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 )
          {
            MuonSystem->cscRechitCluster2_match_dtRechits_0p4[MuonSystem->nCscRechitClusters2] ++;
            if (dtRechitStation[i] == 1) MuonSystem->cscRechitCluster2_match_MB1_0p4[MuonSystem->nCscRechitClusters2] ++;
          }

          if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.6 )
          {
            MuonSystem->cscRechitCluster2_match_dtRechits_0p6[MuonSystem->nCscRechitClusters2] ++;
            if (dtRechitStation[i] == 1) MuonSystem->cscRechitCluster2_match_MB1_0p6[MuonSystem->nCscRechitClusters2] ++;
          }

        }
        //match to MB1 DT segments
        for (int i = 0; i < nDtSeg; i++) {
          if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 )
          {
            MuonSystem->cscRechitCluster2_match_dtSeg_0p4[MuonSystem->nCscRechitClusters2] ++;
            if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster2_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters2] ++;
          }

          if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.6 )
          {
            MuonSystem->cscRechitCluster2_match_dtSeg_0p6[MuonSystem->nCscRechitClusters2] ++;
            if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster2_match_MB1Seg_0p6[MuonSystem->nCscRechitClusters2] ++;
          }

        }
        //match to RPC hits in RE1/2
        for (int i = 0; i < nRpc; i++) {
          float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
          if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 )
          {
            if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)
            {
              MuonSystem->cscRechitCluster2_match_RE12_0p4[MuonSystem->nCscRechitClusters2] ++;
            }
            if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
            {
              MuonSystem->cscRechitCluster2_match_RB1_0p4[MuonSystem->nCscRechitClusters2] ++;
            }

          }
          // if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.6 )
          // {
          //   if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)
          //   {
          //     MuonSystem->cscRechitCluster2_match_RE12_0p6[MuonSystem->nCscRechitClusters2] ++;
          //   }
          //   if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
          //   {
          //     MuonSystem->cscRechitCluster2_match_RB1_0p6[MuonSystem->nCscRechitClusters2] ++;
          //   }
          //
          // }

        }





        MuonSystem->cscRechitCluster2Met_dPhi[MuonSystem->nCscRechitClusters2] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2],MuonSystem->metPhi);
        MuonSystem->cscRechitCluster2MetXYCorr_dPhi[MuonSystem->nCscRechitClusters2] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2],MuonSystem->metPhiXYCorr);

        // cout<<MuonSystem->nCscRechitClusters2<<endl;
        MuonSystem->nCscRechitClusters2++;
    }
    //isolation

    for(int i = 0; i < MuonSystem->nCscRechitClusters2;i++)
    {
      float closest_cluster_dR = 15;
      int index = 999;
      for(int j = 0; j < MuonSystem->nCscRechitClusters2;j++)
      {
        float dR_temp = RazorAnalyzer::deltaR(MuonSystem->cscRechitCluster2Eta[i],MuonSystem->cscRechitCluster2Phi[i], MuonSystem->cscRechitCluster2Eta[j],MuonSystem->cscRechitCluster2Phi[j]);
        if ( dR_temp < closest_cluster_dR )
        {
          closest_cluster_dR = dR_temp;
          index = j;
        }
      }
      MuonSystem->cscRechitCluster2_match_cluster_dR[i] = closest_cluster_dR;
      MuonSystem->cscRechitCluster2_match_cluster_index[i] = index;
    }
    // cluster merging

    ds.merge_clusters();
    ds.result();
    ds.clusterMoments();
    ds.sort_clusters();




    MuonSystem->nCscRechitClusters3 = 0;
    for ( auto &tmp : ds.CscCluster ) {

        MuonSystem->cscRechitCluster3X[MuonSystem->nCscRechitClusters3] =tmp.x;
        MuonSystem->cscRechitCluster3Y[MuonSystem->nCscRechitClusters3] =tmp.y;
        MuonSystem->cscRechitCluster3Z[MuonSystem->nCscRechitClusters3] =tmp.z;
        MuonSystem->cscRechitCluster3Time[MuonSystem->nCscRechitClusters3] = tmp.t;
        MuonSystem->cscRechitCluster3TimeTotal[MuonSystem->nCscRechitClusters3] = tmp.tTotal;
        MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3] =tmp.eta;
        MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3] = tmp.phi;
        MuonSystem->cscRechitCluster3MajorAxis[MuonSystem->nCscRechitClusters3] =tmp.MajorAxis;
        MuonSystem->cscRechitCluster3MinorAxis[MuonSystem->nCscRechitClusters3] =tmp.MinorAxis;
        MuonSystem->cscRechitCluster3XSpread[MuonSystem->nCscRechitClusters3] =tmp.XSpread;
        MuonSystem->cscRechitCluster3YSpread[MuonSystem->nCscRechitClusters3] =tmp.YSpread;
        MuonSystem->cscRechitCluster3ZSpread[MuonSystem->nCscRechitClusters3] =tmp.ZSpread;
        MuonSystem->cscRechitCluster3EtaPhiSpread[MuonSystem->nCscRechitClusters3] =tmp.EtaPhiSpread;
        MuonSystem->cscRechitCluster3XYSpread[MuonSystem->nCscRechitClusters3] =tmp.XYSpread;
        MuonSystem->cscRechitCluster3RSpread[MuonSystem->nCscRechitClusters3] =tmp.RSpread;
        if(!isData)
        {
          for (int j = 0; j < N_phicorr; j++)
          {
            for (int k = 0; k < N_rcorr; k++)
            {
              MuonSystem->cscRechitCluster3RSpread_corr[MuonSystem->nCscRechitClusters3][j][k] = tmp.RSpread_corr[j][k];
              MuonSystem->cscRechitCluster3XYSpread_corr[MuonSystem->nCscRechitClusters3][j][k] = tmp.XYSpread_corr[j][k];
              MuonSystem->cscRechitCluster3XSpread_corr[MuonSystem->nCscRechitClusters3][j][k] = tmp.XSpread_corr[j][k];
              MuonSystem->cscRechitCluster3YSpread_corr[MuonSystem->nCscRechitClusters3][j][k] = tmp.YSpread_corr[j][k];
              MuonSystem->cscRechitCluster3PhiSpread_corr[MuonSystem->nCscRechitClusters3][j][k] = tmp.PhiSpread_corr[j][k];
              MuonSystem->cscRechitCluster3EtaSpread_corr[MuonSystem->nCscRechitClusters3][j][k] = tmp.EtaSpread_corr[j][k];
              MuonSystem->cscRechitCluster3EtaPhiSpread_corr[MuonSystem->nCscRechitClusters3][j][k] = tmp.EtaPhiSpread_corr[j][k];
              MuonSystem->cscRechitCluster3RSpread_corr[MuonSystem->nCscRechitClusters3][j][k] = tmp.RSpread_corr[j][k];
            }
          }

          int phi_corr_index = 0;
          int r_corr_index = 0;
          if (abs(tmp.eta) < 1)
          {
            phi_corr_index = 5;
            r_corr_index = 3;
          }
          else if (abs(tmp.eta) < 1.2)
          {
            phi_corr_index = 4;
            r_corr_index = 2;
          }
          else if (abs(tmp.eta) < 1.4)
          {
            phi_corr_index = 2;
            r_corr_index = 1;
          }
          else if (abs(tmp.eta) < 1.6)
          {
            phi_corr_index = 1;
            r_corr_index = 1;
          }
          else if (abs(tmp.eta) < 1.8)
          {
            phi_corr_index = 5;
            r_corr_index = 1;
          }
          else if (abs(tmp.eta) < 2)
          {
            phi_corr_index = 5;
            r_corr_index = 0;
          }
          MuonSystem->cscRechitCluster3RSpread_corr2[MuonSystem->nCscRechitClusters3] = tmp.RSpread_corr[phi_corr_index][r_corr_index];
          MuonSystem->cscRechitCluster3XYSpread_corr2[MuonSystem->nCscRechitClusters3] = tmp.XYSpread_corr[phi_corr_index][r_corr_index];
          MuonSystem->cscRechitCluster3XSpread_corr2[MuonSystem->nCscRechitClusters3] = tmp.XSpread_corr[phi_corr_index][r_corr_index];
          MuonSystem->cscRechitCluster3YSpread_corr2[MuonSystem->nCscRechitClusters3] = tmp.YSpread_corr[phi_corr_index][r_corr_index];
          MuonSystem->cscRechitCluster3PhiSpread_corr2[MuonSystem->nCscRechitClusters3] = tmp.PhiSpread_corr[phi_corr_index][r_corr_index];
          MuonSystem->cscRechitCluster3EtaSpread_corr2[MuonSystem->nCscRechitClusters3] = tmp.EtaSpread_corr[phi_corr_index][r_corr_index];
          MuonSystem->cscRechitCluster3EtaPhiSpread_corr2[MuonSystem->nCscRechitClusters3] = tmp.EtaPhiSpread_corr[phi_corr_index][r_corr_index];
          MuonSystem->cscRechitCluster3RSpread_corr2[MuonSystem->nCscRechitClusters3] = tmp.RSpread_corr[phi_corr_index][r_corr_index];


        }


        MuonSystem->cscRechitCluster3EtaSpread[MuonSystem->nCscRechitClusters3] =tmp.EtaSpread;
        MuonSystem->cscRechitCluster3PhiSpread[MuonSystem->nCscRechitClusters3] = tmp.PhiSpread;
        MuonSystem->cscRechitCluster3TimeSpread[MuonSystem->nCscRechitClusters3] = tmp.TSpread;
        MuonSystem->cscRechitCluster3Size[MuonSystem->nCscRechitClusters3] = tmp.nCscSegments;

        MuonSystem->cscRechitCluster3NRechitChamberPlus11[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus11;
        MuonSystem->cscRechitCluster3NRechitChamberPlus12[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus12;
        MuonSystem->cscRechitCluster3NRechitChamberPlus13[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus13;
        MuonSystem->cscRechitCluster3NRechitChamberPlus21[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus21;
        MuonSystem->cscRechitCluster3NRechitChamberPlus22[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus22;
        MuonSystem->cscRechitCluster3NRechitChamberPlus31[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus31;
        MuonSystem->cscRechitCluster3NRechitChamberPlus32[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus32;
        MuonSystem->cscRechitCluster3NRechitChamberPlus41[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus41;
        MuonSystem->cscRechitCluster3NRechitChamberPlus42[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus42;
        MuonSystem->cscRechitCluster3NRechitChamberMinus11[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus11;
        MuonSystem->cscRechitCluster3NRechitChamberMinus12[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus12;
        MuonSystem->cscRechitCluster3NRechitChamberMinus13[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus13;
        MuonSystem->cscRechitCluster3NRechitChamberMinus21[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus21;
        MuonSystem->cscRechitCluster3NRechitChamberMinus22[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus22;
        MuonSystem->cscRechitCluster3NRechitChamberMinus31[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus31;
        MuonSystem->cscRechitCluster3NRechitChamberMinus32[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus32;
        MuonSystem->cscRechitCluster3NRechitChamberMinus41[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus41;
        MuonSystem->cscRechitCluster3NRechitChamberMinus42[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus42;
        MuonSystem->cscRechitCluster3MaxChamber[MuonSystem->nCscRechitClusters3] = tmp.maxChamber;
        MuonSystem->cscRechitCluster3MaxChamberRatio[MuonSystem->nCscRechitClusters3] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
        MuonSystem->cscRechitCluster3NChamber[MuonSystem->nCscRechitClusters3] = tmp.nChamber;
        MuonSystem->cscRechitCluster3MaxStation[MuonSystem->nCscRechitClusters3] = tmp.maxStation;
        MuonSystem->cscRechitCluster3MaxStationRatio[MuonSystem->nCscRechitClusters3] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
        MuonSystem->cscRechitCluster3NStation[MuonSystem->nCscRechitClusters2] = tmp.nStation;
        MuonSystem->cscRechitCluster3NStation5[MuonSystem->nCscRechitClusters2] = tmp.nStation5;
        MuonSystem->cscRechitCluster3NStation10perc[MuonSystem->nCscRechitClusters2] = tmp.nStation10perc;
        MuonSystem->cscRechitCluster3AvgStation[MuonSystem->nCscRechitClusters2] = tmp.avgStation;
        MuonSystem->cscRechitCluster3AvgStation5[MuonSystem->nCscRechitClusters2] = tmp.avgStation5;
        MuonSystem->cscRechitCluster3AvgStation10perc[MuonSystem->nCscRechitClusters2] = tmp.avgStation10perc;
        MuonSystem->cscRechitCluster3Me11Ratio[MuonSystem->nCscRechitClusters3] = tmp.Me11Ratio;
        MuonSystem->cscRechitCluster3Me12Ratio[MuonSystem->nCscRechitClusters3] = tmp.Me12Ratio;
        if(MuonSystem->category == 2)
        {
          if (RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex1], MuonSystem->lepPhi[MuonSystem->ZleptonIndex1], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4) {
            MuonSystem->cscRechitCluster3ZLep1[MuonSystem->nCscRechitClusters3] = true;
            MuonSystem->cscRechitCluster3ZLep1Id[MuonSystem->nCscRechitClusters3] = MuonSystem->lepPdgId[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep1LooseIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassLooseIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep1TightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep1VTightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep1VVTightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassVVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep1TightId[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassId[MuonSystem->ZleptonIndex1];

          }
          if (RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex2], MuonSystem->lepPhi[MuonSystem->ZleptonIndex2], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4) {
            MuonSystem->cscRechitCluster3ZLep2[MuonSystem->nCscRechitClusters3] = true;
            MuonSystem->cscRechitCluster3ZLep2Id[MuonSystem->nCscRechitClusters3] = MuonSystem->lepPdgId[MuonSystem->ZleptonIndex2];
            MuonSystem->cscRechitCluster3ZLep2LooseIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassLooseIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep2TightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep2VTightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep2VVTightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassVVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep2TightId[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassId[MuonSystem->ZleptonIndex1];
          }
        }
        //Jet veto/ muon veto
        MuonSystem->cscRechitCluster3JetVetoPt[MuonSystem->nCscRechitClusters3] = 0.0;
        MuonSystem->cscRechitCluster3JetVetoE[MuonSystem->nCscRechitClusters3] = 0.0;
        MuonSystem->cscRechitCluster3MuonVetoPt[MuonSystem->nCscRechitClusters3] = 0.0;
        MuonSystem->cscRechitCluster3MuonVetoE[MuonSystem->nCscRechitClusters3] = 0.0;
        MuonSystem->cscRechitCluster3GenMuonVetoPt[MuonSystem->nCscRechitClusters3] = 0.0;
        MuonSystem->cscRechitCluster3GenMuonVetoE[MuonSystem->nCscRechitClusters3] = 0.0;
        MuonSystem->cscRechitCluster3IsoMuonVetoPt[MuonSystem->nCscRechitClusters3] = 0.0;

        // jet veto
        for(int i = 0; i < nJets; i++)
        {
          if (fabs(jetEta[i]>3.0)) continue;
          if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 && jetPt[i] > MuonSystem->cscRechitCluster3JetVetoPt[MuonSystem->nCscRechitClusters3] ) {
            MuonSystem->cscRechitCluster3JetVetoPt[MuonSystem->nCscRechitClusters3]  = jetPt[i];
          }
          if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 && jetE[i] > MuonSystem->cscRechitCluster3JetVetoE[MuonSystem->nCscRechitClusters3] ) {
            MuonSystem->cscRechitCluster3JetVetoE[MuonSystem->nCscRechitClusters3]  = jetE[i];
          }
        }
        float min_deltaR = 15.;
        int index = 999;
        // for(int i = 0; i < nCscRechitClusters; i++)
        // {
        //   double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3], cscRechitClusterEta[i], cscRechitClusterPhi[i]);
        //   if (current_delta_r < 0.4 && MuonSystem->cscRechitCluster3JetVetoPt[MuonSystem->nCscRechitClusters3] < cscRechitClusterJetVetoPt[i])
        //   {
        //     MuonSystem->cscRechitCluster3JetVetoPt[MuonSystem->nCscRechitClusters3]  = cscRechitClusterJetVetoPt[i];
        //     MuonSystem->cscRechitCluster3JetVetoE[MuonSystem->nCscRechitClusters3]  = cscRechitClusterJetVetoE[i];
        //   }
        //   if (current_delta_r < 0.6 && MuonSystem->cscRechitCluster3JetVetoPt[MuonSystem->nCscRechitClusters3] < cscRechitClusterJetVetoPt[i])
        //   {
        //     MuonSystem->cscRechitCluster3JetVetoPt_0p6[MuonSystem->nCscRechitClusters3]  = cscRechitClusterJetVetoPt[i];
        //     MuonSystem->cscRechitCluster3JetVetoE_0p6[MuonSystem->nCscRechitClusters3]  = cscRechitClusterJetVetoE[i];
        //   }
        //   if (current_delta_r < 0.8 && MuonSystem->cscRechitCluster3JetVetoPt[MuonSystem->nCscRechitClusters3] < cscRechitClusterJetVetoPt[i])
        //   {
        //     MuonSystem->cscRechitCluster3JetVetoPt_0p8[MuonSystem->nCscRechitClusters3]  = cscRechitClusterJetVetoPt[i];
        //     MuonSystem->cscRechitCluster3JetVetoE_0p8[MuonSystem->nCscRechitClusters3]  = cscRechitClusterJetVetoE[i];
        //   }
        //
        //
        // }

        // for(int i = 0; i < nJets; i++)
        // {
        //   if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.6 && jetPt[i] > MuonSystem->cscRechitCluster3JetVetoPt_0p6[MuonSystem->nCscRechitClusters3] ) {
        //     MuonSystem->cscRechitCluster3JetVetoPt_0p6[MuonSystem->nCscRechitClusters3]  = jetPt[i];
        //     MuonSystem->cscRechitCluster3JetVetoE_0p6[MuonSystem->nCscRechitClusters3]  = jetE[i];
        //   }
        //   if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.8 && jetPt[i] > MuonSystem->cscRechitCluster3JetVetoPt_0p8[MuonSystem->nCscRechitClusters3] ) {
        //     MuonSystem->cscRechitCluster3JetVetoPt_0p8[MuonSystem->nCscRechitClusters3]  = jetPt[i];
        //     MuonSystem->cscRechitCluster3JetVetoE_0p8[MuonSystem->nCscRechitClusters3]  = jetE[i];
        //   }
        // }

        for(int i = 0; i < nMuons; i++)
        {
          if (fabs(muonEta[i]>3.0)) continue;
          float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
          if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 && muonPt[i] > MuonSystem->cscRechitCluster3MuonVetoPt[MuonSystem->nCscRechitClusters3] ) {
            MuonSystem->cscRechitCluster3MuonVetoPt[MuonSystem->nCscRechitClusters3]  = muonPt[i];
            MuonSystem->cscRechitCluster3MuonVetoE[MuonSystem->nCscRechitClusters3]  = muonE[i];
            MuonSystem->cscRechitCluster3MuonVetoPhi[MuonSystem->nCscRechitClusters3]  = muonPhi[i];
            MuonSystem->cscRechitCluster3MuonVetoEta[MuonSystem->nCscRechitClusters3]  = muonEta[i];
            MuonSystem->cscRechitCluster3MuonVetoLooseIso[MuonSystem->nCscRechitClusters3]  = muonIso<0.25;
            MuonSystem->cscRechitCluster3MuonVetoTightIso[MuonSystem->nCscRechitClusters3]  = muonIso<0.15;
            MuonSystem->cscRechitCluster3MuonVetoVTightIso[MuonSystem->nCscRechitClusters3]  = muonIso<0.10;
            MuonSystem->cscRechitCluster3MuonVetoVVTightIso[MuonSystem->nCscRechitClusters3]  = muonIso<0.05;
            MuonSystem->cscRechitCluster3MuonVetoTightId[MuonSystem->nCscRechitClusters3]  = isMuonPOGTightMuon(i);

          }
          // if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.6 && muonPt[i] > MuonSystem->cscRechitCluster3MuonVetoPt_0p6[MuonSystem->nCscRechitClusters3] ) {
          //   MuonSystem->cscRechitCluster3MuonVetoPt_0p6[MuonSystem->nCscRechitClusters3]  = muonPt[i];
          //   MuonSystem->cscRechitCluster3MuonVetoE_0p6[MuonSystem->nCscRechitClusters3]  = muonE[i];
          // }
          // if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.8 && muonPt[i] > MuonSystem->cscRechitCluster3MuonVetoPt_0p8[MuonSystem->nCscRechitClusters3] ) {
          //   MuonSystem->cscRechitCluster3MuonVetoPt_0p8[MuonSystem->nCscRechitClusters3]  = muonPt[i];
          //   MuonSystem->cscRechitCluster3MuonVetoE_0p8[MuonSystem->nCscRechitClusters3]  = muonE[i];
          // }
          //check if muon is isolated


        }

        // match to gen-level muon
        if(!isData)
        {
          for(int i = 0; i < nGenParticle; i++)
          {
            if (abs(gParticleId[i])!=13) continue;
            if (abs(gParticleMotherId[i])>24 || abs(gParticleMotherId[i])<23) continue;
            if (abs(gParticleStatus[i])!=1)continue;
            // if (fabs(muonEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 ) {
              MuonSystem->cscRechitCluster3GenMuonVetoPt[MuonSystem->nCscRechitClusters3]  = gParticlePt[i];
              MuonSystem->cscRechitCluster3GenMuonVetoE[MuonSystem->nCscRechitClusters3]  = gParticleE[i];

            }

          }
          min_deltaR = 15.;
          index = 999;
          for(int j = 0; j < nGenParticle; j++)
          {

            double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3], gParticleEta[j], gParticlePhi[j]);

            if (current_delta_r < min_deltaR)
            {
              min_deltaR = current_delta_r;
              index = j;
            }
          }
          if (min_deltaR < 0.4)
          {

            MuonSystem->cscRechitCluster3_match_gParticle[MuonSystem->nCscRechitClusters3] = true;
            MuonSystem->cscRechitCluster3_match_gParticle_minDeltaR[MuonSystem->nCscRechitClusters3] = min_deltaR;
            MuonSystem->cscRechitCluster3_match_gParticle_index[MuonSystem->nCscRechitClusters3] = index;
            MuonSystem->cscRechitCluster3_match_gParticle_id[MuonSystem->nCscRechitClusters3] = gParticleId[index];

            MuonSystem->cscRechitCluster3_match_gParticle_eta[MuonSystem->nCscRechitClusters3] = gParticleEta[index];
            MuonSystem->cscRechitCluster3_match_gParticle_phi[MuonSystem->nCscRechitClusters3] = gParticlePhi[index];
            MuonSystem->cscRechitCluster3_match_gParticle_E[MuonSystem->nCscRechitClusters3] = gParticleE[index];
            MuonSystem->cscRechitCluster3_match_gParticle_pt[MuonSystem->nCscRechitClusters3] = gParticlePt[index];
            MuonSystem->cscRechitCluster3_match_gParticle_MotherId[MuonSystem->nCscRechitClusters3] = gParticleMotherId[index];


          }
          // match to gen level LLP
          min_deltaR = 15.;
          index = 999;
          for(int j = 0; j < 2;j++)
          {

            double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3], gLLP_eta[j], gLLP_phi[j]);
            if (current_delta_r < min_deltaR)
            {
              min_deltaR = current_delta_r;
              index = j;
            }
          }
          if (min_deltaR < 0.4)
          {
            MuonSystem->cscRechitCluster3_match_gLLP[MuonSystem->nCscRechitClusters3] = true;
            MuonSystem->cscRechitCluster3_match_gLLP_minDeltaR[MuonSystem->nCscRechitClusters3] = min_deltaR;
            MuonSystem->cscRechitCluster3_match_gLLP_index[MuonSystem->nCscRechitClusters3] = index;
            MuonSystem->cscRechitCluster3_match_gLLP_eta[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_eta[index];
            MuonSystem->cscRechitCluster3_match_gLLP_phi[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_phi[index];
            MuonSystem->cscRechitCluster3_match_gLLP_decay_r[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_decay_vertex_r[index];
            MuonSystem->cscRechitCluster3_match_gLLP_decay_x[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_decay_vertex_x[index];
            MuonSystem->cscRechitCluster3_match_gLLP_decay_y[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_decay_vertex_y[index];
            MuonSystem->cscRechitCluster3_match_gLLP_decay_z[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_decay_vertex_z[index];
            MuonSystem->cscRechitCluster3_match_gLLP_ctau[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_ctau[index];
            MuonSystem->cscRechitCluster3_match_gLLP_beta[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_beta[index];
            MuonSystem->cscRechitCluster3_match_gLLP_csc[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_csc[index];
          }
        }
        for (int i = 0; i < ncscRechits; i++) {
          if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 )
          {
            MuonSystem->cscRechitCluster3_match_cscRechits_0p4[MuonSystem->nCscRechitClusters3] ++;
          }
          // if(abs(cscRechitsEta[i])>2.1)
          // {
          //   if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 )
          //   {
          //     MuonSystem->cscRechitCluster3_match_highEta_0p4[MuonSystem->nCscRechitClusters3] ++;
          //   }
          //   if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.6 )
          //   {
          //     MuonSystem->cscRechitCluster3_match_highEta_0p6[MuonSystem->nCscRechitClusters3] ++;
          //   }
          //   if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.8 )
          //   {
          //     MuonSystem->cscRechitCluster3_match_highEta_0p8[MuonSystem->nCscRechitClusters3] ++;
          //   }
          // }
          if (!(abs(cscRechitsChamber[i]) == 11 || abs(cscRechitsChamber[i])==12))continue;
          if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 )
          {
            MuonSystem->cscRechitCluster3_match_Me1112_0p4[MuonSystem->nCscRechitClusters3] ++;
            if (abs(cscRechitsChamber[i]) == 11) MuonSystem->cscRechitCluster3_match_Me11_0p4[MuonSystem->nCscRechitClusters3] ++;
            if (abs(cscRechitsChamber[i]) == 12) MuonSystem->cscRechitCluster3_match_Me12_0p4[MuonSystem->nCscRechitClusters3] ++;

          }
          if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.6 )
          {
            MuonSystem->cscRechitCluster3_match_Me1112_0p6[MuonSystem->nCscRechitClusters3] ++;
            if (abs(cscRechitsChamber[i]) == 11) MuonSystem->cscRechitCluster3_match_Me11_0p6[MuonSystem->nCscRechitClusters3] ++;
            if (abs(cscRechitsChamber[i]) == 12) MuonSystem->cscRechitCluster3_match_Me12_0p6[MuonSystem->nCscRechitClusters3] ++;
          }
          if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.8 )
          {
            MuonSystem->cscRechitCluster3_match_Me1112_0p8[MuonSystem->nCscRechitClusters3] ++;
            if (abs(cscRechitsChamber[i]) == 11) MuonSystem->cscRechitCluster3_match_Me11_0p8[MuonSystem->nCscRechitClusters3] ++;
            if (abs(cscRechitsChamber[i]) == 12) MuonSystem->cscRechitCluster3_match_Me12_0p8[MuonSystem->nCscRechitClusters3] ++;
          }

        }
        for (int i = 0; i < nCscSeg; i++) {
          if (RazorAnalyzer::deltaR(cscSegEta[i], cscSegPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 )
          {
            MuonSystem->cscRechitCluster3_match_cscSeg_0p4[MuonSystem->nCscRechitClusters3] ++;
            if (abs(cscSegChamber[i]) == 11) MuonSystem->cscRechitCluster3_match_ME11Seg_0p4[MuonSystem->nCscRechitClusters3] ++;
            if (abs(cscSegChamber[i]) == 12) MuonSystem->cscRechitCluster3_match_ME12Seg_0p4[MuonSystem->nCscRechitClusters3] ++;

          }

          if (RazorAnalyzer::deltaR(cscSegEta[i], cscSegPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.6 )
          {
            MuonSystem->cscRechitCluster3_match_cscSeg_0p6[MuonSystem->nCscRechitClusters3] ++;
            if (abs(cscSegChamber[i]) == 11) MuonSystem->cscRechitCluster3_match_ME11Seg_0p6[MuonSystem->nCscRechitClusters3] ++;
            if (abs(cscSegChamber[i]) == 12) MuonSystem->cscRechitCluster3_match_ME12Seg_0p6[MuonSystem->nCscRechitClusters3] ++;
          }
        }
        //match to MB1 DT hits
        for (int i = 0; i < nDtRechits; i++) {
          if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 )
          {
            MuonSystem->cscRechitCluster3_match_dtRechits_0p4[MuonSystem->nCscRechitClusters3] ++;
            if (dtRechitStation[i] == 1) MuonSystem->cscRechitCluster3_match_MB1_0p4[MuonSystem->nCscRechitClusters3] ++;
          }

          if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.6 )
          {
            MuonSystem->cscRechitCluster3_match_dtRechits_0p6[MuonSystem->nCscRechitClusters3] ++;
            if (dtRechitStation[i] == 1) MuonSystem->cscRechitCluster3_match_MB1_0p6[MuonSystem->nCscRechitClusters3] ++;
          }

        }
        //match to MB1 DT segments
        for (int i = 0; i < nDtSeg; i++) {
          if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 )
          {
            MuonSystem->cscRechitCluster3_match_dtSeg_0p4[MuonSystem->nCscRechitClusters3] ++;
            if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster3_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters3] ++;
          }

          if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.6 )
          {
            MuonSystem->cscRechitCluster3_match_dtSeg_0p6[MuonSystem->nCscRechitClusters3] ++;
            if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster3_match_MB1Seg_0p6[MuonSystem->nCscRechitClusters3] ++;
          }

        }
        //match to RPC hits in RE1/2
        for (int i = 0; i < nRpc; i++) {
          float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
          if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 )
          {
            if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)
            {
              MuonSystem->cscRechitCluster3_match_RE12_0p4[MuonSystem->nCscRechitClusters3] ++;
            }
            if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
            {
              MuonSystem->cscRechitCluster3_match_RB1_0p4[MuonSystem->nCscRechitClusters3] ++;
            }

          }
          if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.6 )
          {
            if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)
            {
              MuonSystem->cscRechitCluster3_match_RE12_0p6[MuonSystem->nCscRechitClusters3] ++;
            }
            if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
            {
              MuonSystem->cscRechitCluster3_match_RB1_0p6[MuonSystem->nCscRechitClusters3] ++;
            }

          }

        }





        MuonSystem->cscRechitCluster3Met_dPhi[MuonSystem->nCscRechitClusters3] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3],MuonSystem->metPhi);
        MuonSystem->cscRechitCluster3MetXYCorr_dPhi[MuonSystem->nCscRechitClusters3] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3],MuonSystem->metPhiXYCorr);

        // cout<<MuonSystem->nCscRechitClusters3<<endl;
        MuonSystem->nCscRechitClusters3++;
    }


    //  for (unsigned int j = 0; j < tmp.segment_id.size(); j++)
    // {
    //   MuonSystem->cscRechitsCluster2Id[tmp.segment_id[j]] = MuonSystem->nCscRechitClusters2;
    // }
    // if(MuonSystem->nCscRechitClusters3==0 && MuonSystem->nCscRechitClusters2==0 && MuonSystem->nCscRechitClusters==0) continue;
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
