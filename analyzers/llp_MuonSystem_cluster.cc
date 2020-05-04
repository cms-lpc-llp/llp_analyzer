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
        // MuonSystem->genJetMET[MuonSystem->nGenJets] = genJetMET[i];
        // MuonSystem->nGenJets++;
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
    std::pair<double,double> corrected_met;
    if (analysisTag=="Razor2016_07Aug2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2016, !isData, nPV);
    else if (analysisTag=="Razor2017_17Nov2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2017, !isData, nPV);
    else if (analysisTag=="Razor2018_17SeptEarlyReReco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2018, !isData, nPV);

    MuonSystem->metXYCorr = corrected_met.first;
    MuonSystem->metPhiXYCorr = corrected_met.second;




    // if (MuonSystem->met < 200) continue;
    // if (nCscRechitClusters==0) continue;

    // bool choose = eventNum == 1243610572 || eventNum == 572324158 || eventNum == 1451893272 || eventNum == 44602077 || eventNum == 1670944854 || eventNum == 631780136 ||
    // eventNum == 83056980 || eventNum == 4325526 || eventNum == 593372604 || eventNum == 704776501 || eventNum == 2105752624 || eventNum == 2165046124 ||
    // eventNum == 255458933 || eventNum == 354792156 || eventNum == 1224149306 || eventNum == 64618353 || eventNum == 1069480720 || eventNum == 453016592 ||
    // eventNum == 1270945645 || eventNum == 96810985 || eventNum == 1362289194 || eventNum == 133669729;



    bool choose = eventNum == 100543182 || eventNum == 952719763 || eventNum == 1451893272 || eventNum == 816983307 || eventNum == 490853685
    || eventNum == 225725803 || eventNum == 476534769 || eventNum == 426365350 || eventNum == 1670944854 || eventNum == 631780136
    || eventNum == 83056980 || eventNum == 304638627 || eventNum == 289348625 || eventNum == 593372604 || eventNum == 704776501
    || eventNum == 602262988 || eventNum == 1667434810 || eventNum == 352616181 || eventNum == 182247249 || eventNum == 198364713
    || eventNum == 354792156 || eventNum == 1224149306|| eventNum == 1069480720 || eventNum == 96810985 || eventNum == 1362289194;
    if (!choose)continue;



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
      tmpMuon.passId = isMuonPOGLooseMuon(i);
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
    // if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
    // if ( !jetPassIDLoose[i] ) continue;

    jets tmpJet;
    tmpJet.jet    = thisJet;
    tmpJet.passId = jetPassIDLoose[i];

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
      if (tmp.jet.Pt()>20 )MuonSystem->HT = MuonSystem->HT + tmp.jet.Pt();
      MuonSystem->jetPassId[MuonSystem->nJets] = tmp.passId;
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

      if( tmp.jet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
      if( fabs( tmp.jet.Eta() ) >= 3.0 ) continue;
      if ( !tmp.passId ) continue;
      MuonSystem->nJets++;
    }


    MuonSystem-> jetMet_dPhiMin = jetMet_dPhiMin_temp;
    MuonSystem-> jetMet_dPhiMin4 = jetMet_dPhiMin4_temp;

    //-----------------------------
    // CSC INFO
    //-----------------------------
    // cout<<"number of clusters: "<<nCscRechitClusters<<endl;
   for(int i = 0; i < nCscRechitClusters; i++)
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

    // cout << "Number of rec hits: "<<ncscRechits<<endl;
    vector<Point> points;
    points.clear();
    MuonSystem->nCscRechits  = 0;

    for (int i = 0; i < ncscRechits; i++) {
      // if (cscRechitsQuality[i]!=1) continue;
      // if (cscRechitsQuality[i]>2) cout<<cscRechitsQuality[i]<<endl;
      MuonSystem->cscRechitsPhi[MuonSystem->nCscRechits]           = cscRechitsPhi[i];   //[nCsc]
      MuonSystem->cscRechitsEta[MuonSystem->nCscRechits]           = cscRechitsEta[i];   //[nCsc]
      MuonSystem->cscRechitsX[MuonSystem->nCscRechits]             = cscRechitsX[i];   //[nCsc]
      MuonSystem->cscRechitsY[MuonSystem->nCscRechits]             = cscRechitsY[i];   //[nCsc]
      MuonSystem->cscRechitsZ[MuonSystem->nCscRechits]             = cscRechitsZ[i];   //[nCsc]
      MuonSystem->cscRechitsTpeak[MuonSystem->nCscRechits] = cscRechitsTpeak[i];
      MuonSystem->cscRechitsTwire[MuonSystem->nCscRechits] = cscRechitsTwire[i];
      MuonSystem->cscRechitsQuality[MuonSystem->nCscRechits] = cscRechitsQuality[i];
      MuonSystem->cscRechitsStation[MuonSystem->nCscRechits] = cscRechitsStation[i];
      MuonSystem->cscRechitsChamber[MuonSystem->nCscRechits] = cscRechitsChamber[i];
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
      if (cscRechitsY[i]>=0.0) MuonSystem->nCscPositiveYRechits++;
      else MuonSystem->nCscNegativeYRechits++;
      if (cscRechitsTpeak[i]<-12.5)MuonSystem->nEarlyCscRechits++;
      if (cscRechitsTpeak[i]>12.5)MuonSystem->nLateCscRechits++;
      if (cscRechitsTpeak[i]<-25)MuonSystem->nEarly2CscRechits++;
      if (cscRechitsTpeak[i]>25)MuonSystem->nLate2CscRechits++;
      if (cscRechitsChamber[i] == 11) MuonSystem->nCscRechitsChamberPlus11[chamber-1]++;
      if (cscRechitsChamber[i] == 12) MuonSystem->nCscRechitsChamberPlus12[chamber-1]++;
      if (cscRechitsChamber[i] == 13) MuonSystem->nCscRechitsChamberPlus13[chamber-1]++;
      if (cscRechitsChamber[i] == 21) MuonSystem->nCscRechitsChamberPlus21[chamber-1]++;
      if (cscRechitsChamber[i] == 22) MuonSystem->nCscRechitsChamberPlus22[chamber-1]++;
      if (cscRechitsChamber[i] == 31) MuonSystem->nCscRechitsChamberPlus31[chamber-1]++;
      if (cscRechitsChamber[i] == 32) MuonSystem->nCscRechitsChamberPlus32[chamber-1]++;
      if (cscRechitsChamber[i] == 41) MuonSystem->nCscRechitsChamberPlus41[chamber-1]++;
      if (cscRechitsChamber[i] == 42) MuonSystem->nCscRechitsChamberPlus42[chamber-1]++;
      if (cscRechitsChamber[i] == -11) MuonSystem->nCscRechitsChamberMinus11[chamber-1]++;
      if (cscRechitsChamber[i] == -12) MuonSystem->nCscRechitsChamberMinus12[chamber-1]++;
      if (cscRechitsChamber[i] == -13) MuonSystem->nCscRechitsChamberMinus13[chamber-1]++;
      if (cscRechitsChamber[i] == -21) MuonSystem->nCscRechitsChamberMinus21[chamber-1]++;
      if (cscRechitsChamber[i] == -22) MuonSystem->nCscRechitsChamberMinus22[chamber-1]++;
      if (cscRechitsChamber[i] == -31) MuonSystem->nCscRechitsChamberMinus31[chamber-1]++;
      if (cscRechitsChamber[i] == -32) MuonSystem->nCscRechitsChamberMinus32[chamber-1]++;
      if (cscRechitsChamber[i] == -41) MuonSystem->nCscRechitsChamberMinus41[chamber-1]++;
      if (cscRechitsChamber[i] == -42) MuonSystem->nCscRechitsChamberMinus42[chamber-1]++;
      MuonSystem->nCscRechits++;
    }
    //Do DBSCAN Clustering
    int min_point = 50;  //minimum number of segments to call it a cluster
    float epsilon = 0.2; //cluster radius parameter
    DBSCAN ds(min_point, epsilon, points);
    ds.run();
    ds.result();
    ds.clusterMoments();
    ds.sort_clusters();
    //Save cluster information
    cout<<"done clustering: "<<ds.CscCluster.size()<<endl;
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

        MuonSystem->cscRechitCluster2XYSpread_phi0p5[MuonSystem->nCscRechitClusters2] = tmp.XYSpread_phi0p5;
        MuonSystem->cscRechitCluster2XSpread_phi0p5[MuonSystem->nCscRechitClusters2] = tmp.XSpread_phi0p5;
        MuonSystem->cscRechitCluster2YSpread_phi0p5[MuonSystem->nCscRechitClusters2] = tmp.YSpread_phi0p5;
        MuonSystem->cscRechitCluster2PhiSpread_phi0p5[MuonSystem->nCscRechitClusters2] = tmp.PhiSpread_phi0p5;
        MuonSystem->cscRechitCluster2EtaPhiSpread_phi0p5[MuonSystem->nCscRechitClusters2] = tmp.EtaPhiSpread_phi0p5;
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

        MuonSystem->cscRechitCluster2Me11Ratio[MuonSystem->nCscRechitClusters2] = tmp.Me11Ratio;
        MuonSystem->cscRechitCluster2Me12Ratio[MuonSystem->nCscRechitClusters2] = tmp.Me12Ratio;
        for (unsigned int j = 0; j < tmp.segment_id.size(); j++)
       {
         MuonSystem->cscRechitsCluster2Id[tmp.segment_id[j]] = MuonSystem->nCscRechitClusters2;
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
        // for(int i = 0; i < nJets; i++)
        // {
        //   if (fabs(jetEta[i]>3.0)) continue;
        //   if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 && jetPt[i] > MuonSystem->cscRechitCluster2JetVetoPt[MuonSystem->nCscRechitClusters2] ) {
        //     MuonSystem->cscRechitCluster2JetVetoPt[MuonSystem->nCscRechitClusters2]  = jetPt[i];
        //   }
        //   if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2], MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 && jetE[i] > MuonSystem->cscRechitCluster2JetVetoE[MuonSystem->nCscRechitClusters2] ) {
        //     MuonSystem->cscRechitCluster2JetVetoE[MuonSystem->nCscRechitClusters2]  = jetE[i];
        //   }
        // }
        float min_deltaR = 15.;
        int index = 999;
        for(int i = 0; i < nCscRechitClusters; i++)
        {
          double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2], MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2], cscRechitClusterEta[i], cscRechitClusterPhi[i]);
          if (current_delta_r < min_deltaR)
          {
            min_deltaR = current_delta_r;
            index = i;
          }


        }
        if (min_deltaR < 0.4)
        {
          MuonSystem->cscRechitCluster2JetVetoPt[MuonSystem->nCscRechitClusters2]  = cscRechitClusterJetVetoPt[index];
          MuonSystem->cscRechitCluster2JetVetoE[MuonSystem->nCscRechitClusters2]  = cscRechitClusterJetVetoE[index];

        }


        for(int i = 0; i < nMuons; i++)
        {
          if (fabs(muonEta[i]>3.0)) continue;
          bool muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.25;
          if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 && muonPt[i] > MuonSystem->cscRechitCluster2MuonVetoPt[MuonSystem->nCscRechitClusters2] ) {
            MuonSystem->cscRechitCluster2MuonVetoPt[MuonSystem->nCscRechitClusters2]  = muonPt[i];
            MuonSystem->cscRechitCluster2MuonVetoE[MuonSystem->nCscRechitClusters2]  = muonE[i];
            MuonSystem->cscRechitCluster2MuonVetoPhi[MuonSystem->nCscRechitClusters2]  = muonPhi[i];
            MuonSystem->cscRechitCluster2MuonVetoEta[MuonSystem->nCscRechitClusters2]  = muonEta[i];
            if (muonIso) MuonSystem->cscRechitCluster2MuonVetoIso[MuonSystem->nCscRechitClusters2]  = true;

          }
          //check if muon is isolated
          if (!muonIso) continue;
          if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 && muonPt[i] > MuonSystem->cscRechitCluster2IsoMuonVetoPt[MuonSystem->nCscRechitClusters2] ) {
            MuonSystem->cscRechitCluster2IsoMuonVetoPt[MuonSystem->nCscRechitClusters2]  = muonPt[i];
            MuonSystem->cscRechitCluster2IsoMuonVetoE[MuonSystem->nCscRechitClusters2]  = muonE[i];
            MuonSystem->cscRechitCluster2IsoMuonVetoPhi[MuonSystem->nCscRechitClusters2]  = muonPhi[i];
            MuonSystem->cscRechitCluster2IsoMuonVetoEta[MuonSystem->nCscRechitClusters2]  = muonEta[i];
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
          if (!(abs(cscRechitsChamber[i]) == 11 || abs(cscRechitsChamber[i])==12))continue;
          if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 )
          {
            MuonSystem->cscRechitCluster2_match_Me1112_0p4[MuonSystem->nCscRechitClusters2] ++;
          }
          if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.8 )
          {
            MuonSystem->cscRechitCluster2_match_Me1112_0p8[MuonSystem->nCscRechitClusters2] ++;
          }

        }
        for (int i = 0; i < nDtRechits; i++) {
          if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->cscRechitCluster2Eta[MuonSystem->nCscRechitClusters2],MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2]) < 0.4 )
          {
            MuonSystem->cscRechitCluster2_match_dtRechits_0p4[MuonSystem->nCscRechitClusters2] ++;
            if (dtRechitStation[i] == 1) MuonSystem->cscRechitCluster2_match_MB1_0p4[MuonSystem->nCscRechitClusters2] ++;
          }

        }
        MuonSystem->cscRechitCluster2Met_dPhi[MuonSystem->nCscRechitClusters2] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2],MuonSystem->metPhi);
        MuonSystem->cscRechitCluster2MetXYCorr_dPhi[MuonSystem->nCscRechitClusters2] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitCluster2Phi[MuonSystem->nCscRechitClusters2],MuonSystem->metPhiXYCorr);

        cout<<MuonSystem->nCscRechitClusters2<<endl;
        MuonSystem->nCscRechitClusters2++;
    }


    // if(MuonSystem->nCscRechitClusters2==0 && MuonSystem->nCscRechitClusters==0) continue;
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
