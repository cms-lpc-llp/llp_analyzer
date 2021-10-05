#include "llp_MuonSystem_combine.h"
#include "RazorHelper.h"
#include "TreeMuonSystemCombination.h"
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

  float electronEnergyFraction;
  float neutralEmEnergyFraction;
  float chargedHadronEnergyFraction;
  float neutralHadronEnergyFraction;
  float muonEnergyFraction;

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

void llp_MuonSystem_combine::Analyze(bool isData, int options, string outputfilename, string analysisTag)
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


  TreeMuonSystemCombination *MuonSystem = new TreeMuonSystemCombination;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  // for signals, need one output file for each signal point
  map<pair<int,int>, TFile*> Files2D;
  map<pair<int,int>, TTree*>Trees2D;
  map<pair<int,int>, TH1F*> NEvents2D;
  map<pair<int,int>, TH1F*> accep2D;
  map<pair<int,int>, TH1F*> accep_met2D;
  map<pair<int,int>, TH1F*> Total2D;



  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *Total = new TH1F("Total", "Total", 1, 1, 2);

  TH1F *accep = new TH1F("accep", "acceptance", 1, 1, 2);
  TH1F *accep_met = new TH1F("accep_met", "acceptance_met", 1, 1, 2);

  TH1F *Nmet200 = new TH1F("Nmet200", "Nmet200", 1, 1, 2);
  TH1F *NmetFilter = new TH1F("NmetFilter", "NmetFilter", 1, 1, 2);
  TH1F *Nlep0 = new TH1F("Nlep0", "Nlep0", 1, 1, 2);
  TH1F *Njet1 = new TH1F("Njet1", "Njet1", 1, 1, 2);
  TH1F *NcosmicVeto = new TH1F("NcosmicVeto", "NcosmicVeto", 1, 1, 2);



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
    // std::cout << "deb1 " << jentry << std::endl;

    // std::cout  << ncscRechits << std::endl;




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
        Total2D[signalPair] = new TH1F(Form("Total%d%d", mx, ctau), "Total", 1,0.5,1.5);
        accep2D[signalPair] = new TH1F(Form("accep2D%d%d", mx, ctau), "acceptance", 1,0.5,1.5);
        accep_met2D[signalPair] = new TH1F(Form("accep_met2D%d%d", mx, ctau), "acceptance_met", 1,0.5,1.5);



        cout << "Created new output file " << thisFileName << endl;
      }
          //Fill NEvents hist
      // cout<<"here"<<endl;
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
    if (!isData)
    {
        if (analysisTag=="Razor2016_07Aug2017Rereco") MuonSystem->MC_condition = 2016;
        else if (analysisTag=="Razor2017_17Nov2017Rereco") MuonSystem->MC_condition = 2017;
        else if (analysisTag=="Razor2018_17SeptEarlyReReco") MuonSystem->MC_condition = 2018;

    }
    bool wzFlag = false;
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
        if (abs(gParticleId[i]) == 24 && gParticleStatus[i]==62)
        {
          MuonSystem->gWPt = gParticlePt[i];

        }
        if (abs(gParticleId[i])== 25 || abs(gParticleId[i] == 35))
        {
	        MuonSystem->gHiggsPt = gParticlePt[i];
          MuonSystem->gHiggsEta = gParticleEta[i];
          MuonSystem->gHiggsPhi = gParticlePhi[i];
          MuonSystem->gHiggsE = gParticleE[i];

        }
        // if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == 23)
        // { //  Z->mumu/Z->ee
        //   MuonSystem->ZCategory  = 0;
        //
        // }
        // else if (abs(gParticleId[i]) == 15 && gParticleStatus[i] == 2 && abs(gParticleMotherId[i]) == 23){
        //   //  Z->tautau
        //   MuonSystem->ZCategory  = 0;
        // }
        // else if ((abs(gParticleId[i]) == 12 || abs(gParticleId[i]) == 14 || abs(gParticleId[i]) == 16) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == 23){
        //   //Z->nunu
        //   MuonSystem->ZCategory  = 1;
        // }
        // else if ((abs(gParticleId[i]) < 6) && gParticleStatus[i] == 23 && abs(gParticleMotherId[i]) == 23){
        //   //Z->qq
        //   MuonSystem->ZCategory  = 2;
        // }



        MuonSystem->gParticleStatus[MuonSystem->nGenParticle] = gParticleStatus[i];
        MuonSystem->gParticleId[MuonSystem->nGenParticle]  = gParticleId[i];
        MuonSystem->gParticleMotherId[MuonSystem->nGenParticle]  = gParticleMotherId[i];
        MuonSystem->gParticlePt[MuonSystem->nGenParticle]  = gParticlePt[i];
        MuonSystem->gParticleEta[MuonSystem->nGenParticle]  = gParticleEta[i];
        MuonSystem->gParticlePhi[MuonSystem->nGenParticle]  = gParticlePhi[i];
        MuonSystem->gParticleE[MuonSystem->nGenParticle]  = gParticleE[i];
        // MuonSystem->nGenParticle++;
        // cout<<"genparticles: "<<MuonSystem->nGenParticle<<endl;


      }

      MuonSystem->higgsPtWeight = helper->getHiggsPtWeight(MuonSystem->gHiggsPt);
      //for (unsigned int i = 0; i < 9; i++)
      //{
      //  MuonSystem->higgsPtWeightSys[i] = helper->getHiggsPtWeightSys(MuonSystem->gHiggsPt, i) / MuonSystem->higgsPtWeight;
      //  MuonSystem->scaleWeights[i]= (*scaleWeights)[i]/genWeight;
      //}
      //MuonSystem->sf_facScaleUp = MuonSystem->higgsPtWeightSys[5];
      //MuonSystem->sf_facScaleDown = MuonSystem->higgsPtWeightSys[3];
      //MuonSystem->sf_renScaleUp = MuonSystem->higgsPtWeightSys[7];
      //MuonSystem->sf_renScaleDown = MuonSystem->higgsPtWeightSys[1];
      //MuonSystem->sf_facRenScaleUp = MuonSystem->higgsPtWeightSys[8];
      //MuonSystem->sf_facRenScaleDown = MuonSystem->higgsPtWeightSys[0];

      MuonSystem->genMetPtTrue = genMetPtTrue;
      MuonSystem->genMetPhiTrue = genMetPhiTrue;
      MuonSystem->genMetPtCalo = genMetPtCalo;
      MuonSystem->genMetPhiCalo = genMetPhiCalo;
      for(int i = 0; i < 2;i++)
      {
        MuonSystem->gLLP_eta[i] = gLLP_eta[i];
        MuonSystem->gLLP_phi[i] = gLLP_phi[i];
        MuonSystem->gLLP_e[i] = gLLP_e[i];
        MuonSystem->gLLP_pt[i] = gLLP_pt[i];


        MuonSystem->gLLP_decay_vertex_r[i] = sqrt(gLLP_decay_vertex_x[i]*gLLP_decay_vertex_x[i]+gLLP_decay_vertex_y[i]*gLLP_decay_vertex_y[i]);
        MuonSystem->gLLP_decay_vertex_x[i] = gLLP_decay_vertex_x[i];
        MuonSystem->gLLP_decay_vertex_y[i] = gLLP_decay_vertex_y[i];
        MuonSystem->gLLP_decay_vertex_z[i] = gLLP_decay_vertex_z[i];
        float beta = gLLP_beta[i];
        float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[i], 2) + pow(MuonSystem->gLLP_decay_vertex_z[i],2));
        float gamma = 1.0/sqrt(1-beta*beta);
        MuonSystem->gLLP_ctau[i] = gLLP_decay_vertex/(beta * gamma);
        MuonSystem->gLLP_beta[i] = gLLP_beta[i];

          // if (abs(MuonSystem->gLLP_eta[i]) < 2.4 && abs(MuonSystem->gLLP_eta[i]) > 0.9
          //   && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>568
          //   && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) MuonSystem->gLLP_csc[i] = true;

            if (abs(MuonSystem->gLLP_eta[i]) < 2.4
              && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>400
              && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) MuonSystem->gLLP_csc[i] = true;
            if (abs(MuonSystem->gLLP_decay_vertex_z[i])< 661.0
              && MuonSystem->gLLP_decay_vertex_r[i] < 738.0
               && MuonSystem->gLLP_decay_vertex_r[i] > 380.0) MuonSystem->gLLP_dt[i] = true;


      }
      for(int i = 0; i < 4;i++)
      {
        MuonSystem->gLLP_daughter_pt[i] = gLLP_daughter_pt[i];
        MuonSystem->gLLP_daughter_e[i] = gLLP_daughter_e[i];

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

    }//end of isData

      //get NPU
      MuonSystem->npv = nPV;
      MuonSystem->rho = fixedGridRhoFastjetAll;
      MuonSystem->met = metType1Pt;
      MuonSystem->metPhi = metType1Phi;
      MuonSystem->metJESUp = MuonSystem->met;
      MuonSystem->metJESDown = MuonSystem->met;

      MuonSystem->gLLP_deltaR = deltaR(gLLP_eta[0], gLLP_phi[0], gLLP_eta[1], gLLP_phi[1]);




      MuonSystem->metSF = helper->getMetTriggerSF(MuonSystem->met);

      if(signalScan && !isData)Total2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->higgsPtWeight*MuonSystem->pileupWeight);
      // if(!isData)
      // {
      //   if (MuonSystem->gLLP_csc[0] == false && MuonSystem->gLLP_csc[1] == false)continue;
      // }

      if(signalScan && !isData)accep2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->higgsPtWeight*MuonSystem->pileupWeight);
      else if (!isData) accep->Fill(1.0, genWeight*MuonSystem->higgsPtWeight*MuonSystem->pileupWeight);


      std::pair<double,double> corrected_met;
      if (analysisTag=="Razor2016_07Aug2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2016, !isData, nPV);
      else if (analysisTag=="Razor2017_17Nov2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2017, !isData, nPV);
      else if (analysisTag=="Razor2018_17SeptEarlyReReco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2018, !isData, nPV);

      MuonSystem->metXYCorr = corrected_met.first;
      MuonSystem->metPhiXYCorr = corrected_met.second;

      MuonSystem->metEENoise = MuonSystem->metXYCorr;
      MuonSystem->metHEM = MuonSystem->metXYCorr;
      MuonSystem->metPhiEENoise = MuonSystem->metPhiXYCorr;
      MuonSystem->metPhiHEM = MuonSystem->metPhiXYCorr;
      // if (MuonSystem->met < 200) continue;
      // if (nCscRechitClusters==0) continue;


      //Triggers
      for(int i = 0; i < NTriggersMAX; i++){
        MuonSystem->HLTDecision[i] = HLTDecision[i];
      }

      if (analysisTag=="Razor2016_07Aug2017Rereco" || analysisTag == "Razor2016_Source2018")
      {
        MuonSystem->METTrigger = HLTDecision[310] || HLTDecision[467];
        MuonSystem->METNoMuTrigger = HLTDecision[467];
      }
      else
      {
        MuonSystem->METTrigger = HLTDecision[310] || HLTDecision[467] || HLTDecision[703] || HLTDecision[717] || HLTDecision[710] || HLTDecision[709];
        MuonSystem->METNoMuTrigger = HLTDecision[467] ||  HLTDecision[717] || HLTDecision[710];

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

      if (analysisTag!="Razor2016_07Aug2017Rereco")MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_ecalBadCalibFilter;

      // if (MuonSystem->evtNum!=931158463)continue;

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
        if(muonPt[i] < 25) continue;
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

        tmpMuon.passVetoId = false;
        Leptons.push_back(tmpMuon);
      }

      //-------------------------------
      //Electrons
      //-------------------------------
      for( int i = 0; i < nElectrons; i++ )
      {

        // if (!isEGammaPOGLooseElectron(i, true, true, true, "Summer16")) continue;
        // if (!isEGammaPOGLooseElectron(i, true, false, true, "vid")) continue;
        // if(!isEGammaPOGLooseElectron(i, true, true, true, "2017_94X"))continue;
        if (!isEGammaPOGLooseElectron(i, true, false, true, "vid")) continue;

        if(elePt[i] < 35) continue;
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
      MuonSystem->category = MuonSystem->nLeptons;
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
    float MetXCorr_JESUp = 0.;
    float MetYCorr_JESUp = 0.;
    float MetXCorr_JESDown = 0.;
    float MetYCorr_JESDown = 0.;
    float MetXCorr_HEM = 0.;
    float MetYCorr_HEM = 0.;
    float MetXCorr_EENoise = 0.;
    float MetYCorr_EENoise = 0.;

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


      double jetCorrPt = jetPt[i];
      double jetCorrE = jetE[i];
      TLorentzVector thisJet = makeTLorentzVector( jetCorrPt, jetEta[i], jetPhi[i], jetCorrE );



      if (thisJet.Eta()>-3.0 && thisJet.Eta()<-1.3 && thisJet.Phi() >-1.57 && thisJet.Phi() <-0.87 && analysisTag == "Razor2018_17SeptEarlyReReco")
      {
        MetXCorr_HEM += thisJet.Px();
        MetYCorr_HEM += thisJet.Py();
      }
      if (fabs(thisJet.Eta())> 2.65 && fabs(thisJet.Eta())<3.139 && thisJet.Pt() < 50  && analysisTag == "Razor2017_17Nov2017Rereco")
      {
        MetXCorr_EENoise += thisJet.Px();
        MetYCorr_EENoise += thisJet.Py();
        if (eventNum==123969624)
        {
          cout<<i<<", "<<thisJet.Eta()<<","<<jetEta[i]<<endl;
          cout<<MetXCorr_EENoise<<", "<<MetYCorr_EENoise<<", "<<thisJet.Px()<<", "<<thisJet.Py()<<","<<thisJet.Pt()<<","<<jetPt[i]<<endl;
        }
      }
      if (fabs(thisJet.Eta())> 2.25 && fabs(thisJet.Eta())<3.0 && thisJet.Pt() > 100 && (analysisTag == "Razor2016_07Aug2017Rereco" || analysisTag == "Razor2017_17Nov2017Rereco"))
      {
        MuonSystem->EE_prefiring = false;
      }


      if (fabs(thisJet.Eta()) >= 3.0)continue;

      jets tmpJet;
      tmpJet.jet    = thisJet;
      tmpJet.passId = isPFTightJet(i, true,analysisTag);

      tmpJet.jetPassMuFrac = jetPassMuFrac[i];
      tmpJet.jetNeutralHadronEnergyFraction = jetNeutralHadronEnergyFraction[i];
      // tmpJet.jetNeutralEMEnergyFraction = jetNeutralEMEnergyFraction[i];
      // tmpJet.jetChargedEMEnergyFraction = jetChargedEMEnergyFraction[i];
      // tmpJet.jetChargedHadronEnergyFraction = jetChargedHadronEnergyFraction[i];


      // calculate jet energy scale uncertainty
      double unc = helper->getJecUnc( jetCorrPt, jetEta[i], runNum ); //use run=999 as default
      tmpJet.jetPtJESUp = jetCorrPt*(1+unc);
      tmpJet.jetPtJESDown = jetCorrPt*(1-unc);
      tmpJet.jetEJESUp = jetCorrE*(1+unc);
      tmpJet.jetEJESDown = jetCorrE*(1-unc);
      tmpJet.JecUnc = unc;
      TLorentzVector thisJetJESUp = makeTLorentzVector(tmpJet.jetPtJESUp, jetEta[i], jetPhi[i], tmpJet.jetEJESUp);
      TLorentzVector thisJetJESDown = makeTLorentzVector(tmpJet.jetPtJESDown, jetEta[i], jetPhi[i], tmpJet.jetEJESDown);
      if (tmpJet.jetPtJESUp > 10)
      {
        MetXCorr_JESUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
        MetYCorr_JESUp += -1 * (thisJetJESUp.Py() - thisJet.Py());

      }
      if (tmpJet.jetPtJESDown > 10)
      {
        MetXCorr_JESDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
        MetYCorr_JESDown += -1 * (thisJetJESDown.Py() - thisJet.Py());

      }
      // if ( !jetPassIDLoose[i] ) continue;
      // if(!isPFTightJet(i, true,analysisTag))continue;
      if( thisJet.Pt() < 30 ) continue;//According to the April 1st 2015 AN
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
        if(tmp.jet.Pt()<30)continue;
        if(abs(tmp.jet.Eta())>2.4)continue;
        MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
        MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
        MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
        MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
        MuonSystem->jetTime[MuonSystem->nJets] = tmp.time;
        MuonSystem->jetPtJESUp[MuonSystem->nJets] = tmp.jetPtJESUp;
        MuonSystem->jetPtJESDown[MuonSystem->nJets] = tmp.jetPtJESDown;
        MuonSystem->jetEJESUp[MuonSystem->nJets] = tmp.jetEJESUp;
        MuonSystem->jetEJESDown[MuonSystem->nJets] = tmp.jetEJESDown;
        MuonSystem->JecUnc[MuonSystem->nJets] = tmp.JecUnc;
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
        // MuonSystem->jetChargedEMEnergyFraction[MuonSystem->nJets] = tmp.jetChargedEMEnergyFraction;
        // MuonSystem->jetNeutralEMEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralEMEnergyFraction;
        // MuonSystem->jetChargedHadronEnergyFraction[MuonSystem->nJets] = tmp.jetChargedHadronEnergyFraction;
        // MuonSystem->jetNeutralHadronEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralHadronEnergyFraction;



        // MuonSystem->jetPassMuFrac[MuonSystem->nJets] = tmp.jetPassMuFrac;


        MuonSystem->nJets++;
      }


      MuonSystem-> jetMet_dPhiMin = jetMet_dPhiMin_temp;
      MuonSystem-> jetMet_dPhiMin4 = jetMet_dPhiMin4_temp;
      // TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);
      TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(MuonSystem->metXYCorr, 0, MuonSystem->metPhiXYCorr, 0);

      //JES up
      float PFMetXJESUp   = PFMET.Px() + MetXCorr_JESUp;
      float PFMetYJESUp   = PFMET.Py() + MetYCorr_JESUp;
      MuonSystem->metJESUp    = sqrt( pow(PFMetXJESUp,2) + pow(PFMetYJESUp,2) );
      MuonSystem->metPhiJESUp    = atan(PFMetYJESUp/PFMetXJESUp);
      if  (PFMetXJESUp < 0.0) MuonSystem->metPhiJESUp = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiJESUp,0.0);
      MuonSystem->metJESUpSF = helper->getMetTriggerSF(MuonSystem->metJESUp)/MuonSystem->metSF;

      //JES down
      float PFMetXJESDown   = PFMET.Px() + MetXCorr_JESDown;
      float PFMetYJESDown   = PFMET.Py() + MetYCorr_JESDown;
      MuonSystem->metJESDown    =  sqrt( pow(PFMetXJESDown,2) + pow(PFMetYJESDown,2) );
      MuonSystem->metPhiJESDown    = atan(PFMetYJESDown/PFMetXJESDown);
      if  (PFMetXJESUp < 0.0) MuonSystem->metPhiJESDown = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiJESDown,0.0);
      MuonSystem->metJESDownSF = helper->getMetTriggerSF(MuonSystem->metJESDown)/MuonSystem->metSF;

      //HEM
      float PFMetXHEM   = PFMET.Px() + MetXCorr_HEM;
      float PFMetYHEM   = PFMET.Py() + MetYCorr_HEM;
      MuonSystem->metHEM    = sqrt( pow(PFMetXHEM,2) + pow(PFMetYHEM,2) );
      MuonSystem->metPhiHEM    = atan(PFMetYHEM/PFMetXHEM);
      if  (PFMetXHEM < 0.0) MuonSystem->metPhiHEM = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiHEM,0.0);


      //EENoise
      float PFMetXEENoise   = PFMET.Px() + MetXCorr_EENoise;
      float PFMetYEENoise   = PFMET.Py() + MetYCorr_EENoise;
      MuonSystem->metEENoise    = sqrt( pow(PFMetXEENoise,2) + pow(PFMetYEENoise,2) );
      MuonSystem->metPhiEENoise    = atan(PFMetYEENoise/PFMetXEENoise);
      if  (PFMetXEENoise < 0.0) MuonSystem->metPhiEENoise = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiEENoise,0.0);

     //if (!MuonSystem->METNoMuTrigger) continue;
     //if (MuonSystem->metEENoise < 200) continue;

      if(signalScan && !isData)accep_met2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->higgsPtWeight*MuonSystem->pileupWeight*MuonSystem->metSF);
      else if(!isData) accep_met->Fill(1.0, genWeight*MuonSystem->higgsPtWeight*MuonSystem->pileupWeight*MuonSystem->metSF);
      else Nmet200->Fill(1.0);


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



        if (dtRechitStation[i] == 1) MuonSystem->nDTRechitsStation1++;
        if (dtRechitStation[i] == 2) MuonSystem->nDTRechitsStation2++;
        if (dtRechitStation[i] == 3) MuonSystem->nDTRechitsStation3++;
        if (dtRechitStation[i] == 4) MuonSystem->nDTRechitsStation4++;

        if (dtRechitWheel[i] == -2) MuonSystem->nDTRechitsWheelMinus2++;
        if (dtRechitWheel[i] == -1) MuonSystem->nDTRechitsWheelMinus1++;
        if (dtRechitWheel[i] == 0) MuonSystem->nDTRechitsWheel0++;
        if (dtRechitWheel[i] == 1) MuonSystem->nDTRechitsWheelPlus1++;
        if (dtRechitWheel[i] == 2) MuonSystem->nDTRechitsWheelPlus2++;

        // MuonSystem->dtRechitsX[MuonSystem->nDTRechits] =  dtRechitX[i];
        // MuonSystem->dtRechitsY[MuonSystem->nDTRechits] =  dtRechitY[i];
        // MuonSystem->dtRechitsZ[MuonSystem->nDTRechits] =  dtRechitZ[i];
        // MuonSystem->dtRechitsEta[MuonSystem->nDTRechits] =  dtRechitEta[i];
        // MuonSystem->dtRechitsPhi[MuonSystem->nDTRechits] = dtRechitPhi[i];
        // MuonSystem->dtRechitsWheel[MuonSystem->nDTRechits] = dtRechitWheel[i];
        // MuonSystem->dtRechitsStation[MuonSystem->nDTRechits] = dtRechitStation[i];
        // MuonSystem->nDTRechits++;
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


      if (MuonSystem->nDTRechitsStation1 > 25) MuonSystem->nDtStations25++;
      if (MuonSystem->nDTRechitsStation2 > 25) MuonSystem->nDtStations25++;
      if (MuonSystem->nDTRechitsStation3 > 25) MuonSystem->nDtStations25++;
      if (MuonSystem->nDTRechitsStation4 > 25) MuonSystem->nDtStations25++;

      if (MuonSystem->nDTRechitsWheelMinus2 > 25) MuonSystem->nDtWheels25++;
      if (MuonSystem->nDTRechitsWheelMinus1 > 25) MuonSystem->nDtWheels25++;
      if (MuonSystem->nDTRechitsWheel0 > 25) MuonSystem->nDtWheels25++;
      if (MuonSystem->nDTRechitsWheelPlus1 > 25) MuonSystem->nDtWheels25++;
      if (MuonSystem->nDTRechitsWheelPlus2 > 25) MuonSystem->nDtWheels25++;


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

        int layer = (cscRechitsDetId[i] & 07);

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
        p.layer = layer;
        p.superlayer = 999;
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


      //ds.merge_clusters();
      //ds.result();
      //ds.clusterMoments();
      //ds.sort_clusters();




      MuonSystem->nCscRechitClusters = 0;
      for ( auto &tmp : ds.CscCluster ) {

          MuonSystem->cscRechitClusterX[MuonSystem->nCscRechitClusters] =tmp.x;
          MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] =tmp.y;
          MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] =tmp.z;
          MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] = tmp.t;
          MuonSystem->cscRechitClusterTimeTotal[MuonSystem->nCscRechitClusters] = tmp.tTotal;
          MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] =tmp.eta;
          MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters] = tmp.phi;

          MuonSystem->cscRechitClusterTimeSpread[MuonSystem->nCscRechitClusters] = tmp.TSpread;
          MuonSystem->cscRechitClusterSize[MuonSystem->nCscRechitClusters] = tmp.nCscSegments;

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

          MuonSystem->cscRechitClusterNStation10[MuonSystem->nCscRechitClusters] = tmp.nStation10;
          MuonSystem->cscRechitClusterAvgStation10[MuonSystem->nCscRechitClusters] = tmp.avgStation10;


          MuonSystem->cscRechitClusterMe11Ratio[MuonSystem->nCscRechitClusters] = tmp.Me11Ratio;
          MuonSystem->cscRechitClusterMe12Ratio[MuonSystem->nCscRechitClusters] = tmp.Me12Ratio;


          MuonSystem->cscRechitClusterNLayersChamberPlus11[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus11;
          MuonSystem->cscRechitClusterNLayersChamberPlus12[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus12;
          MuonSystem->cscRechitClusterNLayersChamberPlus13[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus13;
          MuonSystem->cscRechitClusterNLayersChamberPlus21[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus21;
          MuonSystem->cscRechitClusterNLayersChamberPlus22[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus22;
          MuonSystem->cscRechitClusterNLayersChamberPlus31[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus31;
          MuonSystem->cscRechitClusterNLayersChamberPlus32[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus32;
          MuonSystem->cscRechitClusterNLayersChamberPlus41[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus41;
          MuonSystem->cscRechitClusterNLayersChamberPlus42[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus42;
          MuonSystem->cscRechitClusterNLayersChamberMinus11[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus11;
          MuonSystem->cscRechitClusterNLayersChamberMinus12[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus12;
          MuonSystem->cscRechitClusterNLayersChamberMinus13[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus13;
          MuonSystem->cscRechitClusterNLayersChamberMinus21[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus21;
          MuonSystem->cscRechitClusterNLayersChamberMinus22[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus22;
          MuonSystem->cscRechitClusterNLayersChamberMinus31[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus31;
          MuonSystem->cscRechitClusterNLayersChamberMinus32[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus32;
          MuonSystem->cscRechitClusterNLayersChamberMinus41[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus41;
          MuonSystem->cscRechitClusterNLayersChamberMinus42[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus42;
          //  for (unsigned int j = 0; j < tmp.segment_id.size(); j++)
          // {
          //   MuonSystem->cscRechitsCluster2Id[tmp.segment_id[j]] = MuonSystem->nCscRechitClusters2;
          // }
          // if(isData && MuonSystem->nCscRechitClusters==0) continue;







          //Jet veto/ muon veto
          MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterGenMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterGenMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterIsoMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;

          // jet veto
          for(int i = 0; i < nJets; i++)
          {
            if (fabs(jetEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]  = jetPt[i];
              MuonSystem->cscRechitClusterJetVetoEta[MuonSystem->nCscRechitClusters]  = jetEta[i];
              MuonSystem->cscRechitClusterJetVetoPhi[MuonSystem->nCscRechitClusters]  = jetPhi[i];

            }
            if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetE[i] > MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters]  = jetE[i];
            }
          }
          float min_deltaR = 15.;
          int index = 999;


          for(int i = 0; i < nMuons; i++)
          {
            if (fabs(muonEta[i]>3.0)) continue;
            float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
            if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]  = muonPt[i];
              MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters]  = muonE[i];
              MuonSystem->cscRechitClusterMuonVetoPhi[MuonSystem->nCscRechitClusters]  = muonPhi[i];
              MuonSystem->cscRechitClusterMuonVetoEta[MuonSystem->nCscRechitClusters]  = muonEta[i];
              MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters]  = muon_isGlobal[i];



              // MuonSystem->cscRechitClusterMuonVetoLooseIso[MuonSystem->nCscRechitClusters]  = muonIso<0.25;
              // MuonSystem->cscRechitClusterMuonVetoTightIso[MuonSystem->nCscRechitClusters]  = muonIso<0.15;
              // MuonSystem->cscRechitClusterMuonVetoVTightIso[MuonSystem->nCscRechitClusters]  = muonIso<0.10;
              // MuonSystem->cscRechitClusterMuonVetoVVTightIso[MuonSystem->nCscRechitClusters]  = muonIso<0.05;
              // MuonSystem->cscRechitClusterMuonVetoTightId[MuonSystem->nCscRechitClusters]  = isMuonPOGTightMuon(i);
              MuonSystem->cscRechitClusterMuonVetoLooseId[MuonSystem->nCscRechitClusters]  = isMuonPOGLooseMuon(i);


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
              if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 ) {
                MuonSystem->cscRechitClusterGenMuonVetoPt[MuonSystem->nCscRechitClusters]  = gParticlePt[i];
                MuonSystem->cscRechitClusterGenMuonVetoE[MuonSystem->nCscRechitClusters]  = gParticleE[i];

              }

            }
            min_deltaR = 15.;
            index = 999;
            for(int j = 0; j < nGenParticle; j++)
            {

              double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], gParticleEta[j], gParticlePhi[j]);

              if (current_delta_r < min_deltaR)
              {
                min_deltaR = current_delta_r;
                index = j;
              }
            }

            // match to gen level LLP
            min_deltaR = 15.;
            index = 999;
            for(int j = 0; j < 2;j++)
            {

              double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], gLLP_eta[j], gLLP_phi[j]);
              if (current_delta_r < min_deltaR)
              {
                min_deltaR = current_delta_r;
                index = j;
              }
            }
            if (min_deltaR < 0.4)MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = true;
            else MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = false;


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
            MuonSystem->cscRechitCluster_match_gLLP_dt[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_dt[index];
            MuonSystem->cscRechitCluster_match_gLLP_e[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_e[index];



          }


          //match to MB1 DT segments
          for (int i = 0; i < nDtSeg; i++) {
            if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
            {
              MuonSystem->cscRechitCluster_match_dtSeg_0p4[MuonSystem->nCscRechitClusters] ++;
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



          // MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhi);
          // MuonSystem->cscRechitClusterMetXYCorr_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhiXYCorr);
          MuonSystem->cscRechitClusterMetEENoise_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhiEENoise);
          // MuonSystem->cscRechitClusterMetHEM_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhiHEM);

          MuonSystem->nCscRechitClusters++;
      }



      // DT cluster

      points.clear();

      for (int i = 0; i < nDtRechits; i++) {
        Point p;
        p.phi = dtRechitPhi[i];
        p.eta = dtRechitEta[i];
        p.x = dtRechitX[i];
        p.y = dtRechitY[i];
        p.z = dtRechitZ[i];
        p.t = dtRechitTime[i];
        p.twire = dtRechitTime[i];
        p.station = dtRechitStation[i];
        p.chamber = dtRechitWheel[i];
        // p.superlayer = 999;
        // p.superlayer = dtRechitSuperLayer[i];
        p.superlayer = 0;
        p.clusterID = UNCLASSIFIED;
        points.push_back(p);

      }
      //Do DBSCAN Clustering
      int min_point_dt = 50;  //minimum number of segments to call it a cluster
      float epsilon_dt = 0.2; //cluster radius parameter
      DBSCAN ds_dtRechit(min_point_dt, epsilon_dt, points);
      ds_dtRechit.run();
      ds_dtRechit.result();
      ds_dtRechit.clusterMoments();
      ds_dtRechit.sort_clusters();
      ds_dtRechit.merge_clusters();
      ds_dtRechit.result();
      ds_dtRechit.clusterMoments();
      ds_dtRechit.sort_clusters();



      MuonSystem->nDtRechitClusters = 0;

      for ( auto &tmp : ds_dtRechit.CscCluster ) {
        for (unsigned int j = 0; j < tmp.segment_id.size(); j++)
        {
         MuonSystem->dtRechitsClusterId[tmp.segment_id[j]] = MuonSystem->nDtRechitClusters;
        }
          MuonSystem->dtRechitClusterX[MuonSystem->nDtRechitClusters] =tmp.x;
          MuonSystem->dtRechitClusterY[MuonSystem->nDtRechitClusters] =tmp.y;
          MuonSystem->dtRechitClusterZ[MuonSystem->nDtRechitClusters] =tmp.z;
          if (abs(tmp.z) < 126.8) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 0;
          else if (tmp.z > 126.8 && tmp.z < 395.4) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 1;
          else if (tmp.z < -126.8 && tmp.z > -395.4)MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -1;
          else if (tmp.z<0) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -2;
          else MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 2;


          MuonSystem->dtRechitClusterTime[MuonSystem->nDtRechitClusters] =tmp.t;
          MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters] =tmp.eta;
          MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters] =tmp.phi;
          MuonSystem->dtRechitClusterTimeSpread[MuonSystem->nDtRechitClusters] =tmp.TSpread;


          MuonSystem->dtRechitClusterSize[MuonSystem->nDtRechitClusters] = tmp.nCscSegments;


          MuonSystem->dtRechitClusterNSegmentStation1[MuonSystem->nDtRechitClusters] = tmp.nDtSegmentStation1;
        	MuonSystem->dtRechitClusterNSegmentStation2[MuonSystem->nDtRechitClusters] = tmp.nDtSegmentStation2;
        	MuonSystem->dtRechitClusterNSegmentStation3[MuonSystem->nDtRechitClusters] = tmp.nDtSegmentStation3;
        	MuonSystem->dtRechitClusterNSegmentStation4[MuonSystem->nDtRechitClusters] = tmp.nDtSegmentStation4;

        	MuonSystem->dtRechitClusterMaxChamber[MuonSystem->nDtRechitClusters] = tmp.maxChamber;
        	MuonSystem->dtRechitClusterMaxChamberRatio[MuonSystem->nDtRechitClusters] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
        	MuonSystem->dtRechitClusterNChamber[MuonSystem->nDtRechitClusters] = tmp.nChamber;
        	MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters] = tmp.maxStation;
        	MuonSystem->dtRechitClusterMaxStationRatio[MuonSystem->nDtRechitClusters] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
        	MuonSystem->dtRechitClusterNStation[MuonSystem->nDtRechitClusters] = tmp.nStation;
	  MuonSystem->dtRechitClusterNStation10[MuonSystem->nDtRechitClusters] = tmp.nStation10;
          MuonSystem->dtRechitClusterAvgStation10[MuonSystem->nDtRechitClusters] = tmp.avgStation10;

          //Jet veto/ muon veto
          MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = 0.0;


          // jet veto
          for(int i = 0; i < nJets; i++)
          {
            if (fabs(jetEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters]  = jetPt[i];
              MuonSystem->dtRechitClusterJetVetoEta[MuonSystem->nDtRechitClusters]  = jetEta[i];
              MuonSystem->dtRechitClusterJetVetoPhi[MuonSystem->nDtRechitClusters]  = jetPhi[i];
            }
            if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetE[i] > MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters]  = jetE[i];
            }
          }


          for(int i = 0; i < nMuons; i++)
          {
            if (fabs(muonEta[i]>3.0)) continue;
            float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
            if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters]  = muonPt[i];
              MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters]  = muonE[i];
              MuonSystem->dtRechitClusterMuonVetoPhi[MuonSystem->nDtRechitClusters]  = muonPhi[i];
              MuonSystem->dtRechitClusterMuonVetoEta[MuonSystem->nDtRechitClusters]  = muonEta[i];
              MuonSystem->dtRechitClusterMuonVetoGlobal[MuonSystem->nDtRechitClusters]  = muon_isGlobal[i];
              MuonSystem->dtRechitClusterMuonVetoLooseId[MuonSystem->nDtRechitClusters]  = isMuonPOGLooseMuon(i);


            }
          }
          // match to gen-level muon
          float min_deltaR = 15.;
          int index = 999;
          if(!isData)
          {
            for(int i = 0; i < nGenParticle; i++)
            {
              if (abs(gParticleStatus[i])>2)continue;
              double current_delta_r = RazorAnalyzer::deltaR(tmp.eta, tmp.phi, gParticleEta[i], gParticlePhi[i]);
              if (current_delta_r < min_deltaR)
              {
                min_deltaR = current_delta_r;
                index = i;
              }
            }

            MuonSystem->dtRechitCluster_match_gParticle_deltaR[MuonSystem->nDtRechitClusters]  = min_deltaR;
            MuonSystem->dtRechitCluster_match_gParticle_Id[MuonSystem->nDtRechitClusters]  = gParticleId[index];
            MuonSystem->dtRechitCluster_match_gParticle_Pt[MuonSystem->nDtRechitClusters]  = gParticlePt[index];
            MuonSystem->dtRechitCluster_match_gParticle_Eta[MuonSystem->nDtRechitClusters]  = gParticleEta[index];
            MuonSystem->dtRechitCluster_match_gParticle_Phi[MuonSystem->nDtRechitClusters]  = gParticlePhi[index];
            MuonSystem->dtRechitCluster_match_gParticle_E[MuonSystem->nDtRechitClusters]  = gParticleE[index];
            MuonSystem->dtRechitCluster_match_gParticle_Status[MuonSystem->nDtRechitClusters]  = gParticleStatus[index];
            MuonSystem->dtRechitCluster_match_gParticle_MotherId[MuonSystem->nDtRechitClusters]  = gParticleMotherId[index];

          }



          min_deltaR = 15.;
          index = 999;
          for(int j = 0; j < 2;j++)
          {
            double current_delta_r = RazorAnalyzer::deltaR(tmp.eta, tmp.phi, gLLP_eta[j], gLLP_phi[j]);
            if (current_delta_r < min_deltaR)
            {
              min_deltaR = current_delta_r;
              index = j;
            }
          }
          if (min_deltaR < 0.4)MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = true;
          else MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = false;

           MuonSystem->dtRechitCluster_match_gLLP_minDeltaR[MuonSystem->nDtRechitClusters] = min_deltaR;
           MuonSystem->dtRechitCluster_match_gLLP_index[MuonSystem->nDtRechitClusters] = index;
           MuonSystem->dtRechitCluster_match_gLLP_eta[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_eta[index];
           MuonSystem->dtRechitCluster_match_gLLP_phi[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_phi[index];
           MuonSystem->dtRechitCluster_match_gLLP_decay_r[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_r[index];
           MuonSystem->dtRechitCluster_match_gLLP_decay_x[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_x[index];
           MuonSystem->dtRechitCluster_match_gLLP_decay_y[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_y[index];
           MuonSystem->dtRechitCluster_match_gLLP_decay_z[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
           MuonSystem->dtRechitCluster_match_gLLP_ctau[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_ctau[index];
           MuonSystem->dtRechitCluster_match_gLLP_beta[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_beta[index];
           MuonSystem->dtRechitCluster_match_gLLP_csc[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_csc[index];
           MuonSystem->dtRechitCluster_match_gLLP_dt[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_dt[index];
           MuonSystem->dtRechitCluster_match_gLLP_e[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_e[index];


          //match to MB1 DT segments
          MuonSystem->nCscRechits = ncscRechits;
          MuonSystem->nDtSeg=nDtSeg;

          // MuonSystem->dtRechitCluster_match_dtSegT_dR0p4.push_back({});
          //
          // MuonSystem->dtRechitCluster_match_dtSegX_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegY_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegZ_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegPhi_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegEta_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegWheel_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegStation_dR0p4.push_back({});
          //
          //
          // MuonSystem->dtRechitCluster_match_dtSegT_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegX_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegY_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegZ_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegPhi_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegEta_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegWheel_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_dtSegStation_sameStation_dR0p4.push_back({});
          //
          //
          // MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCX_dPhi0p5.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCY_dPhi0p5.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCZ_dPhi0p5.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCPhi_dPhi0p5.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCEta_dPhi0p5.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCRing_dPhi0p5.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCLayer_dPhi0p5.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCSector_dPhi0p5.push_back({});
          //
          //
          //
          // MuonSystem->dtRechitCluster_match_RPCBx_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCX_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCY_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCZ_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCPhi_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCEta_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCRing_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCLayer_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCSector_dR0p4.push_back({});
          //
          // MuonSystem->dtRechitCluster_match_RPCBx_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCX_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCY_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCZ_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCPhi_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCEta_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCRing_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCLayer_sameStation_dR0p4.push_back({});
          // MuonSystem->dtRechitCluster_match_RPCSector_sameStation_dR0p4.push_back({});
         //
         //  for (int i = 0; i < nDtSeg; i++) {
         //    if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5 )
         //    {
         //      if (dtSegStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1Seg_0p5[MuonSystem->nDtRechitClusters] ++;
         //      if (dtSegTime[i] > -9999.0)
         //      {
         //        MuonSystem->dtRechitCluster_match_dtSeg_0p5[MuonSystem->nDtRechitClusters] ++;
         //        MuonSystem->dtRechitCluster_match_dtSegTime_0p5[MuonSystem->nDtRechitClusters] += dtSegTime[i];
         //        if (dtSegStation[i] == MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters])
         //        {
         //          MuonSystem->dtRechitCluster_match_dtSeg_sameStation_0p5[MuonSystem->nDtRechitClusters] ++;
         //          MuonSystem->dtRechitCluster_match_dtSegTime_sameStation_0p5[MuonSystem->nDtRechitClusters] += dtSegTime[i];
         //        }
         //      }
         //
         //    }
         //    if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
         //    {
         //      if (dtSegStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1Seg_0p4[MuonSystem->nDtRechitClusters] ++;
         //      if (dtSegTime[i] > -9999.0)
         //      {
         //        // temp.push_back(dtSegTime[i]);
         //        MuonSystem->dtRechitCluster_match_dtSegT_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegTime[i]);
         //        MuonSystem->dtRechitCluster_match_dtSegX_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegX[i]);
         //        MuonSystem->dtRechitCluster_match_dtSegY_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegY[i]);
         //        MuonSystem->dtRechitCluster_match_dtSegZ_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegZ[i]);
         //        MuonSystem->dtRechitCluster_match_dtSegPhi_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegPhi[i]);
         //        MuonSystem->dtRechitCluster_match_dtSegEta_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegEta[i]);
         //        MuonSystem->dtRechitCluster_match_dtSegWheel_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegWheel[i]);
         //        MuonSystem->dtRechitCluster_match_dtSegStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegStation[i]);
         //
         //        MuonSystem->dtRechitCluster_match_dtSegTime_0p4[MuonSystem->nDtRechitClusters] += dtSegTime[i];
         //        MuonSystem->dtRechitCluster_match_dtSeg_0p4[MuonSystem->nDtRechitClusters] ++;
         //        if (dtSegStation[i] == MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters])
         //        {
         //          MuonSystem->dtRechitCluster_match_dtSeg_sameStation_0p4[MuonSystem->nDtRechitClusters] ++;
         //          MuonSystem->dtRechitCluster_match_dtSegTime_sameStation_0p4[MuonSystem->nDtRechitClusters] += dtSegTime[i];
         //          MuonSystem->dtRechitCluster_match_dtSegT_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegTime[i]);
         //          MuonSystem->dtRechitCluster_match_dtSegX_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegX[i]);
         //          MuonSystem->dtRechitCluster_match_dtSegY_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegY[i]);
         //          MuonSystem->dtRechitCluster_match_dtSegZ_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegZ[i]);
         //          MuonSystem->dtRechitCluster_match_dtSegPhi_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegPhi[i]);
         //          MuonSystem->dtRechitCluster_match_dtSegEta_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegEta[i]);
         //          MuonSystem->dtRechitCluster_match_dtSegWheel_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegWheel[i]);
         //          MuonSystem->dtRechitCluster_match_dtSegStation_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(dtSegStation[i]);
         //        }
         //      }
         //    }
         //
         //  }
         //
         //  MuonSystem->dtRechitCluster_match_dtSegTime_0p5[MuonSystem->nDtRechitClusters] /= MuonSystem->dtRechitCluster_match_dtSeg_0p5[MuonSystem->nDtRechitClusters];
         //  MuonSystem->dtRechitCluster_match_dtSegTime_0p4[MuonSystem->nDtRechitClusters] /= MuonSystem->dtRechitCluster_match_dtSeg_0p4[MuonSystem->nDtRechitClusters];
         //  MuonSystem->dtRechitCluster_match_dtSegTime_sameStation_0p4[MuonSystem->nDtRechitClusters] /= MuonSystem->dtRechitCluster_match_dtSeg_sameStation_0p4[MuonSystem->nDtRechitClusters];
         //  MuonSystem->dtRechitCluster_match_dtSegTime_sameStation_0p5[MuonSystem->nDtRechitClusters] /= MuonSystem->dtRechitCluster_match_dtSeg_sameStation_0p5[MuonSystem->nDtRechitClusters];
         //
         //
         //  //calculate time spread
         //  for (int i = 0; i < nDtSeg; i++) {
         //    if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5 )
         //    {
         //      if (dtSegTime[i] > -9999.0)
         //      {
         //        MuonSystem->dtRechitCluster_match_dtSegTimeSpread_0p5[MuonSystem->nDtRechitClusters] += pow(dtSegTime[i] - MuonSystem->dtRechitCluster_match_dtSegTime_0p5[MuonSystem->nDtRechitClusters],2) ;
         //        if (dtSegStation[i] == MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters])
         //        {
         //          MuonSystem->dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5[MuonSystem->nDtRechitClusters] += pow(dtSegTime[i] - MuonSystem->dtRechitCluster_match_dtSegTime_0p5[MuonSystem->nDtRechitClusters],2) ;
         //        }
         //      }
         //
         //    }
         //    if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
         //    {
         //      if (dtSegTime[i] > -9999.0)
         //      {
         //        MuonSystem->dtRechitCluster_match_dtSegTimeSpread_0p4[MuonSystem->nDtRechitClusters] += pow(dtSegTime[i] - MuonSystem->dtRechitCluster_match_dtSegTime_0p4[MuonSystem->nDtRechitClusters],2) ;
         //        if (dtSegStation[i] == MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters])
         //        {
         //          MuonSystem->dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4[MuonSystem->nDtRechitClusters] += pow(dtSegTime[i] - MuonSystem->dtRechitCluster_match_dtSegTime_0p5[MuonSystem->nDtRechitClusters],2) ;
         //        }
         //      }
         //    }
         //
         //  }
         //
         //  MuonSystem->dtRechitCluster_match_dtSegTimeSpread_0p5[MuonSystem->nDtRechitClusters] = sqrt(MuonSystem->dtRechitCluster_match_dtSegTimeSpread_0p5[MuonSystem->nDtRechitClusters]/MuonSystem->dtRechitCluster_match_dtSeg_0p5[MuonSystem->nDtRechitClusters]);
         //  MuonSystem->dtRechitCluster_match_dtSegTimeSpread_0p4[MuonSystem->nDtRechitClusters] = sqrt(MuonSystem->dtRechitCluster_match_dtSegTimeSpread_0p4[MuonSystem->nDtRechitClusters]/MuonSystem->dtRechitCluster_match_dtSeg_0p4[MuonSystem->nDtRechitClusters]);
         //
         // MuonSystem->dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5[MuonSystem->nDtRechitClusters] = sqrt(MuonSystem->dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5[MuonSystem->nDtRechitClusters]/MuonSystem->dtRechitCluster_match_dtSeg_sameStation_0p5[MuonSystem->nDtRechitClusters]);
         // MuonSystem->dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4[MuonSystem->nDtRechitClusters] = sqrt(MuonSystem->dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4[MuonSystem->nDtRechitClusters]/MuonSystem->dtRechitCluster_match_dtSeg_sameStation_0p4[MuonSystem->nDtRechitClusters]);


          for (int i = 0; i < nDtRechits; i++) {
            if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5 )
            {
              if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters] ++;
            }
            if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
            {
              if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p4[MuonSystem->nDtRechitClusters] ++;
            }
            if(abs(dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])==1 && dtRechitStation[i] == 1)
            {
              if (abs(RazorAnalyzer::deltaPhi(dtRechitPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < TMath::Pi()/4.0 )
              {
                if (dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] == 1)
                {
                  MuonSystem->dtRechitCluster_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters] ++;
                }
                else
                {
                  MuonSystem->dtRechitCluster_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters] ++;
                }
              }
            }



          }

         //  MuonSystem->nRpc = nRpc;
         //  //match to RPC hits with dPhi<0.5 and same wheel in DT
         //  for (int i = 0; i < nRpc; i++) {
         //    float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
         //    if (rpcRegion[i]!=0) continue;
         //    if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5 )
         //    {
         //      if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
         //      {
         //        MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters] ++;
         //        MuonSystem->dtRechitCluster_match_RPCTime_dPhi0p5[MuonSystem->nDtRechitClusters] += rpcBx[i];
         //
         //        MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcBx[i]);
         //        MuonSystem->dtRechitCluster_match_RPCX_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcX[i]);
         //        MuonSystem->dtRechitCluster_match_RPCY_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcY[i]);
         //        MuonSystem->dtRechitCluster_match_RPCZ_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcZ[i]);
         //        MuonSystem->dtRechitCluster_match_RPCPhi_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcPhi[i]);
         //        MuonSystem->dtRechitCluster_match_RPCEta_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcEta[i]);
         //        MuonSystem->dtRechitCluster_match_RPCRing_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcRing[i]);
         //        MuonSystem->dtRechitCluster_match_RPCLayer_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcLayer[i]);
         //        MuonSystem->dtRechitCluster_match_RPCSector_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcSector[i]);
         //
         //      }
         //
         //    }
         //    if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
         //    {
         //        MuonSystem->dtRechitCluster_match_RPChits_dR0p4[MuonSystem->nDtRechitClusters] ++;
         //        MuonSystem->dtRechitCluster_match_RPCTime_dR0p4[MuonSystem->nDtRechitClusters]  += rpcBx[i];
         //        MuonSystem->dtRechitCluster_match_RPCBx_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcBx[i]);
         //        MuonSystem->dtRechitCluster_match_RPCX_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcX[i]);
         //        MuonSystem->dtRechitCluster_match_RPCY_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcY[i]);
         //        MuonSystem->dtRechitCluster_match_RPCZ_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcZ[i]);
         //        MuonSystem->dtRechitCluster_match_RPCPhi_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcPhi[i]);
         //        MuonSystem->dtRechitCluster_match_RPCEta_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcEta[i]);
         //        MuonSystem->dtRechitCluster_match_RPCRing_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcRing[i]);
         //        MuonSystem->dtRechitCluster_match_RPCLayer_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcLayer[i]);
         //        MuonSystem->dtRechitCluster_match_RPCSector_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcSector[i]);
         //        if (rpcStation[i] == MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters])
         //        {
         //          MuonSystem->dtRechitCluster_match_RPChits_sameStation_dR0p4[MuonSystem->nDtRechitClusters] ++;
         //          MuonSystem->dtRechitCluster_match_RPCTime_sameStation_dR0p4[MuonSystem->nDtRechitClusters]  += rpcBx[i];
         //          MuonSystem->dtRechitCluster_match_RPCBx_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcBx[i]);
         //          MuonSystem->dtRechitCluster_match_RPCX_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcX[i]);
         //          MuonSystem->dtRechitCluster_match_RPCY_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcY[i]);
         //          MuonSystem->dtRechitCluster_match_RPCZ_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcZ[i]);
         //          MuonSystem->dtRechitCluster_match_RPCPhi_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcPhi[i]);
         //          MuonSystem->dtRechitCluster_match_RPCEta_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcEta[i]);
         //          MuonSystem->dtRechitCluster_match_RPCRing_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcRing[i]);
         //          MuonSystem->dtRechitCluster_match_RPCLayer_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcLayer[i]);
         //          MuonSystem->dtRechitCluster_match_RPCSector_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcSector[i]);
         //        }
         //    }
         //
         //  }
         // MuonSystem->dtRechitCluster_match_RPCTime_dPhi0p5[MuonSystem->nDtRechitClusters] /= 1.0* MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters];
         // MuonSystem->dtRechitCluster_match_RPCTime_dR0p4[MuonSystem->nDtRechitClusters] /= 1.0* MuonSystem->dtRechitCluster_match_RPChits_dR0p4[MuonSystem->nDtRechitClusters];
         // MuonSystem->dtRechitCluster_match_RPCTime_sameStation_dR0p4[MuonSystem->nDtRechitClusters] /= 1.0* MuonSystem->dtRechitCluster_match_RPChits_sameStation_dR0p4[MuonSystem->nDtRechitClusters];
         // for (int i = 0; i < nRpc; i++) {
         //   if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5 )
         //   {
         //     if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
         //     {
         //       MuonSystem->dtRechitCluster_match_RPCTimeSpread_dPhi0p5[MuonSystem->nDtRechitClusters] += pow(rpcBx[i] - MuonSystem->dtRechitCluster_match_RPCTime_dPhi0p5[MuonSystem->nDtRechitClusters],2);
         //     }
         //
         //   }
         //   if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
         //   {
         //       MuonSystem->dtRechitCluster_match_RPCTimeSpread_dR0p4[MuonSystem->nDtRechitClusters]  += pow(rpcBx[i] - MuonSystem->dtRechitCluster_match_RPCTime_dR0p4[MuonSystem->nDtRechitClusters],2);
         //       if (rpcStation[i] == MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters])
         //       {
         //         MuonSystem->dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[MuonSystem->nDtRechitClusters]  += pow(rpcBx[i] - MuonSystem->dtRechitCluster_match_RPCTime_sameStation_dR0p4[MuonSystem->nDtRechitClusters],2);
         //       }
         //   }
         //
         // }
         //
         // MuonSystem->dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[MuonSystem->nDtRechitClusters] = sqrt(MuonSystem->dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[MuonSystem->nDtRechitClusters]/MuonSystem->dtRechitCluster_match_RPChits_sameStation_dR0p4[MuonSystem->nDtRechitClusters]);
         // MuonSystem->dtRechitCluster_match_RPCTimeSpread_dR0p4[MuonSystem->nDtRechitClusters] = sqrt(MuonSystem->dtRechitCluster_match_RPCTimeSpread_dR0p4[MuonSystem->nDtRechitClusters]/MuonSystem->dtRechitCluster_match_RPChits_dR0p4[MuonSystem->nDtRechitClusters]);
         // MuonSystem->dtRechitCluster_match_RPCTimeSpread_dPhi0p5[MuonSystem->nDtRechitClusters] = sqrt(MuonSystem->dtRechitCluster_match_RPCTimeSpread_dPhi0p5[MuonSystem->nDtRechitClusters]/MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]);
         std::vector<int> dtRechitCluster_match_rpcBx;

         //match to RPC hits with dPhi<0.5 and same wheel in DT
         for (int i = 0; i < nRpc; i++) {
           float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
           if (rpcRegion[i]!=0) continue;
           if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5 )
           {
             if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
             {
               dtRechitCluster_match_rpcBx.push_back(rpcBx[i]);
               MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]++;
               if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters] ++;

             }
           }
           if(RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
           {
             if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_0p4[MuonSystem->nDtRechitClusters] ++;
           }
         }
         int max_occurence = 0;
         int max_bx = -999;
         for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++)
         {
           int counter = 0;
           for(unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j ++)
           {
             if (dtRechitCluster_match_rpcBx[j] == dtRechitCluster_match_rpcBx[l]) counter++;
           }
           if (counter>max_occurence)
           {
             max_occurence = counter;
             max_bx = dtRechitCluster_match_rpcBx[l];
           }
           // cout<<dtRechitCluster_match_rpcBx[l]<<endl;
         }
         // cout<<max_occurence<<","<<max_bx<<endl;
         // cout<<"###################################################"<<endl;
           MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters] = max_bx;

          MuonSystem->dtRechitClusterMetEENoise_dPhi[MuonSystem->nDtRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhiEENoise);





          MuonSystem->nDtRechitClusters++;
        }


        // correct DT clusters
        // points.clear();
        //
        // for (int i = 0; i < nDtRechits; i++) {
        //   Point p;
        //   p.phi = dtRechitCorrectPhi[i];
        //   p.eta = dtRechitCorrectEta[i];
        //   p.x = dtRechitCorrectX[i];
        //   p.y = dtRechitCorrectY[i];
        //   p.z = dtRechitCorrectZ[i];
        //   p.t = dtRechitTime[i];
        //   p.twire = dtRechitTime[i];
        //   p.station = dtRechitStation[i];
        //   p.chamber = dtRechitWheel[i];
        //   p.superlayer = dtRechitSuperLayer[i];
        //   p.clusterID = UNCLASSIFIED;
        //   points.push_back(p);
        //
        // }
        // //Do DBSCAN Clustering
        // // int min_point_dt = 50;  //minimum number of segments to call it a cluster
        // // float epsilon_dt = 0.2; //cluster radius parameter
        // DBSCAN ds_dtRechit_correct(min_point_dt, epsilon_dt, points);
        // ds_dtRechit_correct.run();
        // ds_dtRechit_correct.result();
        // ds_dtRechit_correct.clusterMoments();
        // ds_dtRechit_correct.sort_clusters();
        // ds_dtRechit_correct.merge_clusters();
        // ds_dtRechit_correct.result();
        // ds_dtRechit_correct.clusterMoments();
        // ds_dtRechit_correct.sort_clusters();
        //
        //
        //
        // MuonSystem->nDtRechitClusters2 = 0;
        //
        // for ( auto &tmp : ds_dtRechit_correct.CscCluster ) {
        //   for (unsigned int j = 0; j < tmp.segment_id.size(); j++)
        //   {
        //    MuonSystem->dtRechitsClusterId[tmp.segment_id[j]] = MuonSystem->nDtRechitClusters2;
        //   }
        //     MuonSystem->dtRechitCluster2X[MuonSystem->nDtRechitClusters2] =tmp.x;
        //     MuonSystem->dtRechitCluster2Y[MuonSystem->nDtRechitClusters2] =tmp.y;
        //     MuonSystem->dtRechitCluster2Z[MuonSystem->nDtRechitClusters2] =tmp.z;
        //     if (abs(tmp.z) < 126.8) MuonSystem->dtRechitCluster2Wheel[MuonSystem->nDtRechitClusters2] = 0;
        //     else if (tmp.z > 126.8 && tmp.z < 395.4) MuonSystem->dtRechitCluster2Wheel[MuonSystem->nDtRechitClusters2] = 1;
        //     else if (tmp.z < -126.8 && tmp.z > -395.4)MuonSystem->dtRechitCluster2Wheel[MuonSystem->nDtRechitClusters2] = -1;
        //     else if (tmp.z<0) MuonSystem->dtRechitCluster2Wheel[MuonSystem->nDtRechitClusters2] = -2;
        //     else MuonSystem->dtRechitCluster2Wheel[MuonSystem->nDtRechitClusters2] = 2;
        //
        //
        //     MuonSystem->dtRechitCluster2Time[MuonSystem->nDtRechitClusters2] =tmp.t;
        //     MuonSystem->dtRechitCluster2Eta[MuonSystem->nDtRechitClusters2] =tmp.eta;
        //     MuonSystem->dtRechitCluster2Phi[MuonSystem->nDtRechitClusters2] =tmp.phi;
        //     MuonSystem->dtRechitCluster2TimeSpread[MuonSystem->nDtRechitClusters2] =tmp.TSpread;
        //
        //
        //     MuonSystem->dtRechitCluster2Size[MuonSystem->nDtRechitClusters2] = tmp.nCscSegments;
        //
        //
        //     MuonSystem->dtRechitCluster2NSegmentStation1[MuonSystem->nDtRechitClusters2] = tmp.nDtSegmentStation1;
        //     MuonSystem->dtRechitCluster2NSegmentStation2[MuonSystem->nDtRechitClusters2] = tmp.nDtSegmentStation2;
        //     MuonSystem->dtRechitCluster2NSegmentStation3[MuonSystem->nDtRechitClusters2] = tmp.nDtSegmentStation3;
        //     MuonSystem->dtRechitCluster2NSegmentStation4[MuonSystem->nDtRechitClusters2] = tmp.nDtSegmentStation4;
        //
        //     MuonSystem->dtRechitCluster2MaxChamber[MuonSystem->nDtRechitClusters2] = tmp.maxChamber;
        //     MuonSystem->dtRechitCluster2MaxChamberRatio[MuonSystem->nDtRechitClusters2] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
        //     MuonSystem->dtRechitCluster2NChamber[MuonSystem->nDtRechitClusters2] = tmp.nChamber;
        //     MuonSystem->dtRechitCluster2MaxStation[MuonSystem->nDtRechitClusters2] = tmp.maxStation;
        //     MuonSystem->dtRechitCluster2MaxStationRatio[MuonSystem->nDtRechitClusters2] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
        //     MuonSystem->dtRechitCluster2NStation[MuonSystem->nDtRechitClusters2] = tmp.nStation;
        //     MuonSystem->dtRechitCluster2NStation10[MuonSystem->nDtRechitClusters2] = tmp.nStation10;
        //     MuonSystem->dtRechitCluster2AvgStation10[MuonSystem->nDtRechitClusters2] = tmp.avgStation10;
        //
        //     //Jet veto/ muon veto
        //     MuonSystem->dtRechitCluster2JetVetoPt[MuonSystem->nDtRechitClusters2] = 0.0;
        //     MuonSystem->dtRechitCluster2JetVetoE[MuonSystem->nDtRechitClusters2] = 0.0;
        //     MuonSystem->dtRechitCluster2MuonVetoPt[MuonSystem->nDtRechitClusters2] = 0.0;
        //     MuonSystem->dtRechitCluster2MuonVetoE[MuonSystem->nDtRechitClusters2] = 0.0;
        //
        //
        //     // jet veto
        //     for(int i = 0; i < nJets; i++)
        //     {
        //       if (fabs(jetEta[i]>3.0)) continue;
        //       if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitCluster2Eta[MuonSystem->nDtRechitClusters2],MuonSystem->dtRechitCluster2Phi[MuonSystem->nDtRechitClusters2]) < 0.4 && jetPt[i] > MuonSystem->dtRechitCluster2JetVetoPt[MuonSystem->nDtRechitClusters2] ) {
        //         MuonSystem->dtRechitCluster2JetVetoPt[MuonSystem->nDtRechitClusters2]  = jetPt[i];
        //         MuonSystem->dtRechitCluster2JetVetoEta[MuonSystem->nDtRechitClusters2]  = jetEta[i];
        //         MuonSystem->dtRechitCluster2JetVetoPhi[MuonSystem->nDtRechitClusters2]  = jetPhi[i];
        //       }
        //       if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitCluster2Eta[MuonSystem->nDtRechitClusters2], MuonSystem->dtRechitCluster2Phi[MuonSystem->nDtRechitClusters2]) < 0.4 && jetE[i] > MuonSystem->dtRechitCluster2JetVetoE[MuonSystem->nDtRechitClusters2] ) {
        //         MuonSystem->dtRechitCluster2JetVetoE[MuonSystem->nDtRechitClusters2]  = jetE[i];
        //       }
        //     }
        //
        //
        //     for(int i = 0; i < nMuons; i++)
        //     {
        //       if (fabs(muonEta[i]>3.0)) continue;
        //       float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
        //       if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->dtRechitCluster2Eta[MuonSystem->nDtRechitClusters2],MuonSystem->dtRechitCluster2Phi[MuonSystem->nDtRechitClusters2]) < 0.4 && muonPt[i] > MuonSystem->dtRechitCluster2MuonVetoPt[MuonSystem->nDtRechitClusters2] ) {
        //         MuonSystem->dtRechitCluster2MuonVetoPt[MuonSystem->nDtRechitClusters2]  = muonPt[i];
        //         MuonSystem->dtRechitCluster2MuonVetoE[MuonSystem->nDtRechitClusters2]  = muonE[i];
        //         MuonSystem->dtRechitCluster2MuonVetoPhi[MuonSystem->nDtRechitClusters2]  = muonPhi[i];
        //         MuonSystem->dtRechitCluster2MuonVetoEta[MuonSystem->nDtRechitClusters2]  = muonEta[i];
        //         MuonSystem->dtRechitCluster2MuonVetoGlobal[MuonSystem->nDtRechitClusters2]  = muon_isGlobal[i];
        //         MuonSystem->dtRechitCluster2MuonVetoLooseId[MuonSystem->nDtRechitClusters2]  = isMuonPOGLooseMuon(i);
        //
        //
        //       }
        //     }
        //     // match to gen-level muon
        //     float min_deltaR = 15.;
        //     int index = 999;
        //     if(!isData)
        //     {
        //       for(int i = 0; i < nGenParticle; i++)
        //       {
        //         if (abs(gParticleStatus[i])>2)continue;
        //         double current_delta_r = RazorAnalyzer::deltaR(tmp.eta, tmp.phi, gParticleEta[i], gParticlePhi[i]);
        //         if (current_delta_r < min_deltaR)
        //         {
        //           min_deltaR = current_delta_r;
        //           index = i;
        //         }
        //       }
        //
        //       MuonSystem->dtRechitCluster2_match_gParticle_deltaR[MuonSystem->nDtRechitClusters2]  = min_deltaR;
        //       MuonSystem->dtRechitCluster2_match_gParticle_Id[MuonSystem->nDtRechitClusters2]  = gParticleId[index];
        //       MuonSystem->dtRechitCluster2_match_gParticle_Pt[MuonSystem->nDtRechitClusters2]  = gParticlePt[index];
        //       MuonSystem->dtRechitCluster2_match_gParticle_Eta[MuonSystem->nDtRechitClusters2]  = gParticleEta[index];
        //       MuonSystem->dtRechitCluster2_match_gParticle_Phi[MuonSystem->nDtRechitClusters2]  = gParticlePhi[index];
        //       MuonSystem->dtRechitCluster2_match_gParticle_E[MuonSystem->nDtRechitClusters2]  = gParticleE[index];
        //       MuonSystem->dtRechitCluster2_match_gParticle_Status[MuonSystem->nDtRechitClusters2]  = gParticleStatus[index];
        //       MuonSystem->dtRechitCluster2_match_gParticle_MotherId[MuonSystem->nDtRechitClusters2]  = gParticleMotherId[index];
        //
        //     }
        //
        //
        //
        //     min_deltaR = 15.;
        //     index = 999;
        //     for(int j = 0; j < 2;j++)
        //     {
        //       double current_delta_r = RazorAnalyzer::deltaR(tmp.eta, tmp.phi, gLLP_eta[j], gLLP_phi[j]);
        //       if (current_delta_r < min_deltaR)
        //       {
        //         min_deltaR = current_delta_r;
        //         index = j;
        //       }
        //     }
        //     if (min_deltaR < 0.4)MuonSystem->dtRechitCluster2_match_gLLP[MuonSystem->nDtRechitClusters2] = true;
        //     else MuonSystem->dtRechitCluster2_match_gLLP[MuonSystem->nDtRechitClusters2] = false;
        //
        //      MuonSystem->dtRechitCluster2_match_gLLP_minDeltaR[MuonSystem->nDtRechitClusters2] = min_deltaR;
        //      MuonSystem->dtRechitCluster2_match_gLLP_index[MuonSystem->nDtRechitClusters2] = index;
        //      MuonSystem->dtRechitCluster2_match_gLLP_eta[MuonSystem->nDtRechitClusters2] = MuonSystem->gLLP_eta[index];
        //      MuonSystem->dtRechitCluster2_match_gLLP_phi[MuonSystem->nDtRechitClusters2] = MuonSystem->gLLP_phi[index];
        //      MuonSystem->dtRechitCluster2_match_gLLP_decay_r[MuonSystem->nDtRechitClusters2] = MuonSystem->gLLP_decay_vertex_r[index];
        //      MuonSystem->dtRechitCluster2_match_gLLP_decay_x[MuonSystem->nDtRechitClusters2] = MuonSystem->gLLP_decay_vertex_x[index];
        //      MuonSystem->dtRechitCluster2_match_gLLP_decay_y[MuonSystem->nDtRechitClusters2] = MuonSystem->gLLP_decay_vertex_y[index];
        //      MuonSystem->dtRechitCluster2_match_gLLP_decay_z[MuonSystem->nDtRechitClusters2] = MuonSystem->gLLP_decay_vertex_z[index];
        //      MuonSystem->dtRechitCluster2_match_gLLP_ctau[MuonSystem->nDtRechitClusters2] = MuonSystem->gLLP_ctau[index];
        //      MuonSystem->dtRechitCluster2_match_gLLP_beta[MuonSystem->nDtRechitClusters2] = MuonSystem->gLLP_beta[index];
        //      MuonSystem->dtRechitCluster2_match_gLLP_csc[MuonSystem->nDtRechitClusters2] = MuonSystem->gLLP_csc[index];
        //      MuonSystem->dtRechitCluster2_match_gLLP_dt[MuonSystem->nDtRechitClusters2] = MuonSystem->gLLP_dt[index];
        //      MuonSystem->dtRechitCluster2_match_gLLP_e[MuonSystem->nDtRechitClusters2] = MuonSystem->gLLP_e[index];
        //
        //
        //     //match to MB1 DT segments
        //     MuonSystem->nCscRechits = ncscRechits;
        //     MuonSystem->nDtSeg=nDtSeg;
        //
        //
        //
        //     for (int i = 0; i < nDtRechits; i++) {
        //       if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->dtRechitCluster2Eta[MuonSystem->nDtRechitClusters2],MuonSystem->dtRechitCluster2Phi[MuonSystem->nDtRechitClusters2]) < 0.5 )
        //       {
        //         if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster2_match_MB1hits_0p5[MuonSystem->nDtRechitClusters2] ++;
        //       }
        //       if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->dtRechitCluster2Eta[MuonSystem->nDtRechitClusters2],MuonSystem->dtRechitCluster2Phi[MuonSystem->nDtRechitClusters2]) < 0.4 )
        //       {
        //         if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster2_match_MB1hits_0p4[MuonSystem->nDtRechitClusters2] ++;
        //       }
        //       if(abs(dtRechitWheel[i]-MuonSystem->dtRechitCluster2Wheel[MuonSystem->nDtRechitClusters2])==1 && dtRechitStation[i] == 1)
        //       {
        //         if (abs(RazorAnalyzer::deltaPhi(dtRechitPhi[i], MuonSystem->dtRechitCluster2Phi[MuonSystem->nDtRechitClusters2])) < TMath::Pi()/4.0 )
        //         {
        //           if (dtRechitWheel[i]-MuonSystem->dtRechitCluster2Wheel[MuonSystem->nDtRechitClusters2] == 1)
        //           {
        //             MuonSystem->dtRechitCluster2_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters2] ++;
        //           }
        //           else
        //           {
        //             MuonSystem->dtRechitCluster2_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters2] ++;
        //           }
        //         }
        //       }
        //
        //
        //
        //     }
        //
        //    std::vector<int> dtRechitCluster_match_rpcBx;
        //
        //    //match to RPC hits with dPhi<0.5 and same wheel in DT
        //    for (int i = 0; i < nRpc; i++) {
        //      float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
        //      if (rpcRegion[i]!=0) continue;
        //      if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], MuonSystem->dtRechitCluster2Phi[MuonSystem->nDtRechitClusters2])) < 0.5 )
        //      {
        //        if (rpcRing[i] == MuonSystem->dtRechitCluster2Wheel[MuonSystem->nDtRechitClusters2])
        //        {
        //          dtRechitCluster_match_rpcBx.push_back(rpcBx[i]);
        //          MuonSystem->dtRechitCluster2_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters2]++;
        //          if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster2_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters2] ++;
        //
        //        }
        //      }
        //      if(RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitCluster2Eta[MuonSystem->nDtRechitClusters2], MuonSystem->dtRechitCluster2Phi[MuonSystem->nDtRechitClusters2]) < 0.4 )
        //      {
        //        if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster2_match_RB1_0p4[MuonSystem->nDtRechitClusters2] ++;
        //      }
        //    }
        //    int max_occurence = 0;
        //    int max_bx = -999;
        //    for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++)
        //    {
        //      int counter = 0;
        //      for(unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j ++)
        //      {
        //        if (dtRechitCluster_match_rpcBx[j] == dtRechitCluster_match_rpcBx[l]) counter++;
        //      }
        //      if (counter>max_occurence)
        //      {
        //        max_occurence = counter;
        //        max_bx = dtRechitCluster_match_rpcBx[l];
        //      }
        //      // cout<<dtRechitCluster_match_rpcBx[l]<<endl;
        //    }
        //    // cout<<max_occurence<<","<<max_bx<<endl;
        //    // cout<<"###################################################"<<endl;
        //      MuonSystem->dtRechitCluster2_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters2] = max_bx;
        //
        //     MuonSystem->dtRechitCluster2MetEENoise_dPhi[MuonSystem->nDtRechitClusters2] =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitCluster2Phi[MuonSystem->nDtRechitClusters2],MuonSystem->metPhiEENoise);
        //
        //
        //
        //
        //
        //     MuonSystem->nDtRechitClusters2++;
        //   }




      //if (  MuonSystem->nDtRechitClusters == 0 &&   MuonSystem->nCscRechitClusters == 0)continue;
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
           Total2D[filePtr.first]->Write("Total");
           accep2D[filePtr.first]->Write("acceptance");
           accep_met2D[filePtr.first]->Write("acceptance_met");
           filePtr.second->Close();

         }
      }
      else if (!isData)
      {
         cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
         cout << "Writing output trees..." << endl;
         outFile->cd();
         MuonSystem->tree_->Write();
         NEvents->Write();
         accep->Write("acceptance");
         accep_met->Write("acceptance_met");
         outFile->Close();
      }


      else
      {
        cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
        cout << "Writing output trees..." << endl;
        outFile->cd();
        MuonSystem->tree_->Write();
        Nmet200->Write();
        NmetFilter->Write();
        Nlep0->Write();
        Njet1->Write();
        NcosmicVeto->Write();
        NEvents->Write();
        // outFile->Write();
        outFile->Close();
      }
}
