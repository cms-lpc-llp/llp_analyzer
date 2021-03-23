#include "llp_MuonSystem_hnl.h"
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

//jet highest pt comparator
struct largest_pt_jet
{
  inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;

void llp_MuonSystem_hnl::Analyze(bool isData, int options, string outputfilename, string analysisTag)
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
  map<pair<int,int>, TTree*>Trees2D;
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
    bool hnl_model = false;
    bool wzFlag = false;
    // cout<<"ngenparticles: "<<nGenParticle<<endl;
    if (!isData)
    {
      for (int i=0; i < nGenParticle; i++)
      {
        if(abs(gParticleId[i])==9900012 || abs(gParticleId[i])==9900014 || abs(gParticleId[i])==9900016)hnl_model = true;
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
      for (unsigned int i = 0; i < 9; i++)
      {
        MuonSystem->higgsPtWeightSys[i] = helper->getHiggsPtWeightSys(MuonSystem->gHiggsPt, i) / MuonSystem->higgsPtWeight;
        MuonSystem->scaleWeights[i]= (*scaleWeights)[i]/genWeight;
      }
      MuonSystem->sf_facScaleUp = MuonSystem->higgsPtWeightSys[5];
      MuonSystem->sf_facScaleDown = MuonSystem->higgsPtWeightSys[3];
      MuonSystem->sf_renScaleUp = MuonSystem->higgsPtWeightSys[7];
      MuonSystem->sf_renScaleDown = MuonSystem->higgsPtWeightSys[1];
      MuonSystem->sf_facRenScaleUp = MuonSystem->higgsPtWeightSys[8];
      MuonSystem->sf_facRenScaleDown = MuonSystem->higgsPtWeightSys[0];


      // cout<<nGenJets<<endl;
      // MuonSystem->nGenJets = 0;
      // for(int i=0; i < nGenJets; i++)
      // {
      //   // cout<<genJetE[i]<<","<<MuonSystem->nGenJets<<","<<nGenJets<<endl;
      //   MuonSystem->genJetE[MuonSystem->nGenJets] = genJetE[i];
      //   MuonSystem->genJetPt[MuonSystem->nGenJets] = genJetPt[i];
      //   MuonSystem->genJetEta[MuonSystem->nGenJets] = genJetEta[i];
      //   MuonSystem->genJetPhi[MuonSystem->nGenJets] = genJetPhi[i];
      //   // MuonSystem->genJetMET[MuonSystem->nGenJets] = genJetMET[i];
      //   MuonSystem->nGenJets++;
      // }
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
        MuonSystem->gLLP_e[i] = gLLP_e[i];
        MuonSystem->gLLP_pt[i] = gLLP_pt[i];
        MuonSystem->gLLP_lepdPhi[i] = deltaPhi(MuonSystem->gLepPhi,MuonSystem->gLLP_phi[i]);
        if (hnl_model)
        {
          MuonSystem->gLLP_EMFracP[i] = gLLP_daughter_pt[0]*cosh(gLLP_daughter_eta[0])/(gLLP_pt[i]*cosh(gLLP_eta[i]));
          MuonSystem->gLLP_EMFracPz[i] = gLLP_daughter_pt[0]*sinh(gLLP_daughter_eta[0])/(gLLP_pt[i]*sinh(gLLP_eta[i]));
          MuonSystem->gLLP_EMFracE[i] = gLLP_daughter_e[0]/(gLLP_e[i]);
          MuonSystem->gLLP_EMFracEz[i] = gLLP_daughter_e[0]/cosh(gLLP_daughter_eta[0])*sinh(gLLP_daughter_eta[0])/(gLLP_e[i]/cosh(gLLP_eta[i])*sinh(gLLP_eta[i]));

          MuonSystem->gLLP_visE[i] = gLLP_e[i];
          MuonSystem->gLLP_visEz[i] = gLLP_e[i]/cosh(gLLP_eta[i])*sinh(gLLP_eta[i]);
          MuonSystem->gLLP_visP[i] = gLLP_pt[i]*cosh(gLLP_eta[i]);
          MuonSystem->gLLP_visPz[i] = gLLP_pt[i]*sinh(gLLP_eta[i]);
        }


        // cout<<MuonSystem->gLLP_e[i]/gamma<<endl;
        // if (abs(MuonSystem->gLLP_eta[i]) < 2.4 && abs(MuonSystem->gLLP_eta[i]) > 0.9
        //   && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>568
        //   && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) MuonSystem->gLLP_csc[i] = true;
          if (abs(MuonSystem->gLLP_eta[i]) < 2.4 
            && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>400
            && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) MuonSystem->gLLP_csc[i] = true;


      }
      MuonSystem->gLLP_deltaR = deltaR(gLLP_eta[0],gLLP_phi[0],gLLP_eta[1],gLLP_phi[1]);

      for(int i = 0; i < 4;i++)
      {
        MuonSystem->gLLP_daughter_id[i] = gLLP_daughter_id[i];
        MuonSystem->gLLP_daughter_pt[i] = gLLP_daughter_pt[i];
        MuonSystem->gLLP_daughter_eta[i] = gLLP_daughter_eta[i];
        MuonSystem->gLLP_daughter_phi[i] = gLLP_daughter_phi[i];
        MuonSystem->gLLP_daughter_e[i] = gLLP_daughter_e[i];
        MuonSystem->gLLP_daughter_mass[i] = gLLP_daughter_mass[i];


      }
      if(!isData && !hnl_model)//twin higgs model
      {
        // cout<<"hnl"<<MuonSystem->gLLP_daughter_id[0]<<endl;
        int temp_id = 999;
        float temp_energy = 999.;
        int temp_status = 999;
        bool has_kaon = false;
        // if (abs(MuonSystem->gLLP_daughter_id[0])==15) //4tau model
        // {
        //   for (int i=0; i < nGenParticle; i++)
        //   {
        //     if (abs(gParticleMotherId[i])!=15)continue;//only do this for tautau mode
        //     if (abs(gParticleStatus[gParticleMotherIndex[i]])!=2)continue;
        //     if (abs(gParticleMotherId[gParticleMotherIndex[i]])!=9000006)continue;//grandparent is LLP
        //     if (abs(gParticleId[i])==310 || abs(gParticleId[i])==311 || abs(gParticleId[i])==130)has_kaon = true;
        //     if (gParticleId[i] == temp_id && temp_status == gParticleStatus[i] && abs(gParticleE[i]-temp_energy)<0.01)continue;
        //
        //     int llp_index = 1;
        //     if (gParticleMotherId[gParticleMotherIndex[i]]>0) llp_index -= 1;
        //     if (abs(gParticleId[i])==12 || abs(gParticleId[i])==14 || abs(gParticleId[i])==16 || abs(gParticleId[i])==13)continue;
        //     MuonSystem->gLLP_visE[llp_index] += gParticleE[i];
        //     MuonSystem->gLLP_visEz[llp_index] += gParticleE[i]/cosh(gParticleEta[i])*sinh(gParticleEta[i]);
        //     MuonSystem->gLLP_visP[llp_index] += gParticlePt[i]*cosh(gParticleEta[i]);
        //     MuonSystem->gLLP_visPz[llp_index] += gParticlePt[i]*sinh(gParticleEta[i]);
        //
        //     if (abs(gParticleId[i])==11 || abs(gParticleId[i])==22 || abs(gParticleId[i])==111) //EM don't count muon energy abs(gParticleId[i])==13 ||
        //     {
        //       MuonSystem->gLLP_EMFracE[llp_index] += gParticleE[i];
        //       MuonSystem->gLLP_EMFracEz[llp_index] += gParticleE[i]/cosh(gParticleEta[i])*sinh(gParticleEta[i]);
        //       MuonSystem->gLLP_EMFracP[llp_index] += gParticlePt[i]*cosh(gParticleEta[i]);
        //       MuonSystem->gLLP_EMFracPz[llp_index] += gParticlePt[i]*sinh(gParticleEta[i]);
        //     }
        //     temp_id = gParticleId[i];
        //     temp_energy = gParticleE[i];
        //     temp_status = gParticleStatus[i];
        //   }
        //
        //   if(has_kaon)continue;
        // }
        if (abs(MuonSystem->gLLP_daughter_id[0])==5 || abs(MuonSystem->gLLP_daughter_id[0])==1 || abs(MuonSystem->gLLP_daughter_id[0])==15) //4b and 4d
        {
          for (int i=0; i < nGenParticle; i++)
          {
            if (abs(MuonSystem->gLLP_daughter_id[0])==15 && abs(gParticleMotherId[i])!=15)continue;

            if (abs(gParticleStatus[i])!=1 && abs(MuonSystem->gLLP_daughter_id[0])!=15)continue;
            if (abs(gParticleId[i])==12 || abs(gParticleId[i])==14 || abs(gParticleId[i])==16 || abs(gParticleId[i])==13)continue;
            if (gParticleId[i] == temp_id && temp_status == gParticleStatus[i] && abs(gParticleE[i]-temp_energy)<0.01)continue;

            // cout<<i<<"here"<<endl;
            int llp_index = 1;
            bool matched_llp0 = false;
            bool matched_llp1 = false;
            float llp0_dis = sqrt(pow(gParticleProdVertexX[i]-gLLP_decay_vertex_x[0],2)+pow(gParticleProdVertexY[i]-gLLP_decay_vertex_y[0],2)+pow(gParticleProdVertexZ[i]-gLLP_decay_vertex_z[0],2));
            float llp1_dis = sqrt(pow(gParticleProdVertexX[i]-gLLP_decay_vertex_x[1],2)+pow(gParticleProdVertexY[i]-gLLP_decay_vertex_y[1],2)+pow(gParticleProdVertexZ[i]-gLLP_decay_vertex_z[1],2));
            // if (RazorAnalyzer::deltaR(gParticleEta[i],gParticlePhi[i], gLLP_eta[0], gLLP_phi[0]) < 0.4 && llp0_dis<100) matched_llp0 = true;
            // if (RazorAnalyzer::deltaR(gParticleEta[i],gParticlePhi[i], gLLP_eta[1], gLLP_phi[1]) < 0.4 && llp1_dis<100) matched_llp1 = true;

            if (llp0_dis<100) matched_llp0 = true;
            if (llp1_dis<100) matched_llp1 = true;

            if (matched_llp0 && matched_llp1) //both matched
            {
              if (llp0_dis<llp1_dis)llp_index = 0;
              else llp_index = 1;
            }
            else if( matched_llp0)
            {
              llp_index = 0;
            }
            else if(matched_llp1)
            {
              llp_index = 1;
            }
            else continue; //both not mached

            if (abs(gParticleId[i])==310 || abs(gParticleId[i])==311 || abs(gParticleId[i])==130)has_kaon = true; //remove matched kaons

            MuonSystem->gLLP_multiplicity[llp_index] += 1;

            MuonSystem->gLLP_visE[llp_index] += gParticleE[i];
            MuonSystem->gLLP_visEz[llp_index] += gParticleE[i]/cosh(gParticleEta[i])*sinh(gParticleEta[i]);
            MuonSystem->gLLP_visP[llp_index] += gParticlePt[i]*cosh(gParticleEta[i]);
            MuonSystem->gLLP_visPz[llp_index] += gParticlePt[i]*sinh(gParticleEta[i]);

            if (abs(gParticleId[i])==11 || abs(gParticleId[i])==22 || abs(gParticleId[i])==111) //EM don't count muon energy abs(gParticleId[i])==13 ||
            {
              MuonSystem->gLLP_EMFracE[llp_index] += gParticleE[i];
              MuonSystem->gLLP_EMFracEz[llp_index] += gParticleE[i]/cosh(gParticleEta[i])*sinh(gParticleEta[i]);
              MuonSystem->gLLP_EMFracP[llp_index] += gParticlePt[i]*cosh(gParticleEta[i]);
              MuonSystem->gLLP_EMFracPz[llp_index] += gParticlePt[i]*sinh(gParticleEta[i]);
              MuonSystem->gLLP_EM_multiplicity[llp_index] += 1;
            }
            temp_id = gParticleId[i];
            temp_energy = gParticleE[i];
            temp_status = gParticleStatus[i];
          }
          if(has_kaon)continue;
        }
        else //ee, pipi，kk

        {
          for(int i = 0; i < 4;i++)
          {
            int llp_index = 1;
            if (i<2) llp_index = 0;

            MuonSystem->gLLP_visE[llp_index] += gLLP_daughter_e[i];
            MuonSystem->gLLP_visEz[llp_index] += gLLP_daughter_e[i]/cosh(gLLP_daughter_eta[i])*sinh(gLLP_daughter_eta[i]);
            MuonSystem->gLLP_visP[llp_index] += gLLP_daughter_pt[i]*cosh(gLLP_daughter_eta[i]);
            MuonSystem->gLLP_visPz[llp_index] += gLLP_daughter_pt[i]*sinh(gLLP_daughter_eta[i]);
            if (abs(gLLP_daughter_id[i])==11 || abs(gLLP_daughter_id[i])==22 || abs(gLLP_daughter_id[i])==111) //EM don't count muon energy abs(gParticleId[i])==13 ||
            {
              MuonSystem->gLLP_EMFracE[llp_index] += gLLP_daughter_e[i];
              MuonSystem->gLLP_EMFracEz[llp_index] += gLLP_daughter_e[i]/cosh(gLLP_daughter_eta[i])*sinh(gLLP_daughter_eta[i]);
              MuonSystem->gLLP_EMFracP[llp_index] += gLLP_daughter_pt[i]*cosh(gLLP_daughter_eta[i]);
              MuonSystem->gLLP_EMFracPz[llp_index] += gLLP_daughter_pt[i]*sinh(gLLP_daughter_eta[i]);
            }

          }
        }
        for (int i=0; i < 2; i++)
        {
          if (MuonSystem->gLLP_visE[i] > 0)
          {
            MuonSystem->gLLP_EMFracE[i] /= MuonSystem->gLLP_visE[i];
            MuonSystem->gLLP_EMFracEz[i] /= MuonSystem->gLLP_visEz[i];
            MuonSystem->gLLP_EMFracP[i] /= MuonSystem->gLLP_visP[i];
            MuonSystem->gLLP_EMFracPz[i] /= MuonSystem->gLLP_visPz[i];
          }
        }
        // else MuonSystem->gLLP_daughterKaon[llp_index] = false;
      }

      // MuonSystem->gLLP_daughter_deltaR[0] = deltaR(gLLP_eta[0],gLLP_phi[0],gLLP_eta[1],gLLP_phi[1]);
      // MuonSystem->gLLP_daughter_deltaR[1] = deltaR(gLLP_eta[2],gLLP_phi[2],gLLP_eta[3],gLLP_phi[3]);
      // cout<<"LLP"<<endl;
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

    MuonSystem->metSF = helper->getMetTriggerSF(MuonSystem->met);

    if(signalScan && !isData)Total2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->higgsPtWeight*MuonSystem->pileupWeight);





    std::pair<double,double> corrected_met;
    if (analysisTag=="Razor2016_07Aug2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2016, !isData, nPV);
    else if (analysisTag=="Razor2017_17Nov2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2017, !isData, nPV);
    else if (analysisTag=="Razor2018_17SeptEarlyReReco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2018, !isData, nPV);

    MuonSystem->metXYCorr = corrected_met.first;
    MuonSystem->metPhiXYCorr = corrected_met.second;

    if(!isData && !hnl_model)
    {
      if (MuonSystem->gLLP_csc[0] == false && MuonSystem->gLLP_csc[1] == false)continue;
    }

    if(signalScan && !isData)accep2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->higgsPtWeight*MuonSystem->pileupWeight);
   else if (!isData) accep->Fill(1.0, genWeight*MuonSystem->higgsPtWeight*MuonSystem->pileupWeight);


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

    // if (MuonSystem->met < 200 || !MuonSystem->METNoMuTrigger) continue;

    if(signalScan && !isData)accep_met2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->higgsPtWeight*MuonSystem->pileupWeight*MuonSystem->metSF);
    else if(!isData) accep_met->Fill(1.0, genWeight*MuonSystem->higgsPtWeight*MuonSystem->pileupWeight);
    else Nmet200->Fill(1.0);

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

    if (analysisTag!="Razor2016_07Aug2017Rereco" && analysisTag!="Razor2016_Source2018")
    {
      MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_ecalBadCalibFilter;
    }

    // if (isData){
    //   if (!MuonSystem->Flag2_all) continue;
    //   NmetFilter->Fill(1.0);
    // }

    //*************************************************************************
    //Start Object Selection
    //*************************************************************************

    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
    for( int i = 0; i < nMuons; i++ )
    {
      float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];


      if(!isMuonPOGLooseMuon(i)) continue;
      if(muonPt[i] < 20) continue;


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
      tmpMuon.passId = isMuonPOGTightMuon(i, true, false);
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

      // if (!isEGammaPOGLooseElectron(i, true, true, true, "Summer16")) continue;
      // if (!isEGammaPOGLooseElectron(i, true, true, true, "2017_94X")) continue;


      if (!isEGammaPOGLooseElectron(i, true, false, true, "vid")) continue;


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
      tmpElectron.passId = isEGammaPOGLooseElectron(i, true, true, true, "2017_94X");

      tmpElectron.passVetoId = isEGammaPOGLooseElectron(i, true, true, true, "vid");
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
      for( uint j = i+1; j < Leptons.size(); j++ )
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

    // if (foundZ  && Leptons.size() == 2 && leadingLepPt > zh_lepton0_cut)
    if (foundZ  && Leptons.size() == 2 )
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
    bool tag = false;
    for ( auto &tmp : Leptons )
    {
      MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
      MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
      MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
      MuonSystem->lepPassId[MuonSystem->nLeptons] = tmp.passId;
      MuonSystem->lepPassVetoId[MuonSystem->nLeptons] = tmp.passVetoId;

      MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
      MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
      MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
      MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
      if (!isData)
      {
        MuonSystem->lepTriggerSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), true, true, true);
        MuonSystem->lepTightIdSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, false);
        MuonSystem->lepLooseIdSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, true);
        MuonSystem->lepTightIsoSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, false);
        MuonSystem->lepLooseIsoSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, true);
        MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), true, true, true);
        MuonSystem->lepTightIdMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, false);
        MuonSystem->lepLooseIdMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, true);
        MuonSystem->lepTightIsoMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, false);
        MuonSystem->lepLooseIsoMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, true);
      }



      if (MuonSystem->lepPassTightIso[MuonSystem->nLeptons] && MuonSystem->lepPassId[MuonSystem->nLeptons] && !tag)
      {
        MuonSystem->lepTag[MuonSystem->nLeptons] = true;
        if(!isData)
        {
          MuonSystem->lepSF[MuonSystem->nLeptons] = MuonSystem->lepTriggerSF[MuonSystem->nLeptons] * MuonSystem->lepTightIdSF[MuonSystem->nLeptons] * MuonSystem->lepTightIsoSF[MuonSystem->nLeptons];
          MuonSystem->lepEff[MuonSystem->nLeptons] = MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepTightIdMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepTightIsoMCEfficiency[MuonSystem->nLeptons];
        }
        else{
          MuonSystem->lepSF[MuonSystem->nLeptons] = 1.0;
          MuonSystem->lepEff[MuonSystem->nLeptons] = 1.0;
        }
        tag = true;
      }
      else
      {
        MuonSystem->lepTag[MuonSystem->nLeptons] = false;
        if(!isData)
        {
          MuonSystem->lepSF[MuonSystem->nLeptons] = MuonSystem->lepTriggerSF[MuonSystem->nLeptons] * MuonSystem->lepLooseIdSF[MuonSystem->nLeptons] * MuonSystem->lepLooseIsoSF[MuonSystem->nLeptons];
          MuonSystem->lepEff[MuonSystem->nLeptons] = MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepLooseIdMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepLooseIsoMCEfficiency[MuonSystem->nLeptons];
        }
        else{
          MuonSystem->lepSF[MuonSystem->nLeptons] = 1.0;
          MuonSystem->lepEff[MuonSystem->nLeptons] = 1.0;
        }
      }



      if(!isData)
      {
        for (int i=0; i < nGenParticle; i++)
        {
          if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && abs(gParticleMotherId[i]) == 23) MuonSystem->lepFromZ[MuonSystem->nLeptons] = true;
        }
      }
      MuonSystem->nLeptons++;
    }
    // if(isData)
    // {
    //   if(MuonSystem->nLeptons>0)continue;
    //   Nlep0->Fill(1.0);
    // }

    // if(MuonSystem->nLeptons!=2)continue;
    // if(MuonSystem->category!=2)continue;
    // if (abs(MuonSystem->lepPdgId[0])!=13)continue;
    // if (abs(MuonSystem->ZMass)<50)continue;
    // if (abs(MuonSystem->lepPt[0])<30)continue;
    // if (abs(MuonSystem->lepPt[1])<30)continue;
    // if(MuonSystem->nLeptons == 2)
    // {
    //   if(MuonSystem->lepTag[0] == MuonSystem->lepTag[1]) continue; //require one tag one probe
    //   if(!isData)
    //   {
    //     MuonSystem->lepOverallSF = 1.0 - (1.0 - MuonSystem->lepSF[0] * MuonSystem->lepEff[0]) * (1 - MuonSystem->lepSF[1] * MuonSystem->lepEff[1]);
    //     MuonSystem->lepOverallSF = MuonSystem->lepOverallSF / (1.0 - (1.0 - MuonSystem->lepEff[0]) * (1 - MuonSystem->lepEff[1]));
    //   }
    //   else{
    //     MuonSystem->lepOverallSF = 1.0;
    //   }
    //
    //
    // }
    TLorentzVector met;
    met.SetPtEtaPhiE(metType1Pt,0,metType1Phi,metType1Pt);
    if ( Leptons.size() > 0 )
    {
      TLorentzVector visible = Leptons[0].lepton;
      MuonSystem->MT = GetMT(visible,met);
    }

    // cout<<"njets:"<<nJets<<endl;
  //-----------------------------------------------
  //Select Jets
  //-----------------------------------------------
  std::vector<jets> Jets;
  float MetXCorr_JESUp = 0;
  float MetYCorr_JESUp = 0;
  float MetXCorr_JESDown = 0;
  float MetYCorr_JESDown = 0;
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
    JEC = 1.0;
    double jetCorrPt = jetPt[i]*JEC;
    double jetCorrE = jetE[i]*JEC;
    // cout<<JEC<<endl;
    // cout<<"corrected pt: "<<jetCorrPt<<", "<<"eta: "<<jetEta[i]<<endl;
    TLorentzVector thisJet = makeTLorentzVector( jetCorrPt, jetEta[i], jetPhi[i], jetCorrE );

    // if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
    if( fabs( thisJet.Eta() ) >= 3.0 ) continue;


    jets tmpJet;
    tmpJet.jet    = thisJet;
    tmpJet.passId = jetPassIDTight[i];

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
    if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
    if ( !jetPassIDLoose[i] ) continue;

    Jets.push_back(tmpJet);

    }
    // cout<<"here"<<nJets<<endl;

    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);
    // cout<<"here"<<nJets<<endl;
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
      if(tmp.jet.Pt()<50)continue;
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
      MuonSystem->jetTightPassId[MuonSystem->nJets] = tmp.passId;


      // float min_deltaR = 15.;
      // int index = 999;
      // for(int i=0; i < nGenJets; i++)
      // {
      //
      //   double current_delta_r = RazorAnalyzer::deltaPhi(genJetPhi[i],jetPhi[MuonSystem->nJets]);
      //   if (current_delta_r < min_deltaR)
      //   {
      //     min_deltaR = current_delta_r;
      //     index = i;
      //   }
      // }
      // if (min_deltaR < 0.4)
      // {
      //   MuonSystem->jet_match_genJet_minDeltaR[MuonSystem->nJets] = min_deltaR;
      //   MuonSystem->jet_match_genJet_index[MuonSystem->nJets] = index;
      //   MuonSystem->jet_match_genJet_pt[MuonSystem->nJets] = jetPt[index];
      // }

      MuonSystem->nJets++;
    }
    // if (isData && MuonSystem->nJets==0)continue;
    // cout<<"here"<<MuonSystem->nJets<<endl;

    Njet1->Fill(1.0);
    MuonSystem-> jetMet_dPhiMin = jetMet_dPhiMin_temp;
    MuonSystem-> jetMet_dPhiMin4 = jetMet_dPhiMin4_temp;

    TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);

    //JES up
  	float PFMetXJESUp   = PFMET.Px() + MetXCorr_JESUp;
  	float PFMetYJESUp   = PFMET.Py() + MetYCorr_JESUp;
  	MuonSystem->metJESUp    = sqrt( pow(PFMetXJESUp,2) + pow(PFMetYJESUp,2) );
    MuonSystem->metJESUpSF = helper->getMetTriggerSF(MuonSystem->metJESUp)/MuonSystem->metSF;

  	//JES down
  	float PFMetXJESDown   = PFMET.Px() + MetXCorr_JESDown;
  	float PFMetYJESDown   = PFMET.Py() + MetYCorr_JESDown;
  	MuonSystem->metJESDown    =  sqrt( pow(PFMetXJESDown,2) + pow(PFMetYJESDown,2) );
    MuonSystem->metJESDownSF = helper->getMetTriggerSF(MuonSystem->metJESDown)/MuonSystem->metSF;

    // cout<<"met"<<endl;

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



    // cout<<"jets"<<endl;




    // cout << "Number of rec hits: "<<ncscRechits<<endl;
    vector<Point> points;
    vector<int> cscRechitsClusterId;
    points.clear();
    MuonSystem->nCscRechits  = 0;

    for (int i = 0; i < ncscRechits; i++) {
      // if(abs(cscRechitsStation[i])!=2)continue;
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
    // if(isData && (MuonSystem->nDtRings+MuonSystem->nCscRings)>=10)continue;
    NcosmicVeto->Fill(1.0);

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
    /*MuonSystem->nCscRechitClusters2 = 0;
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
        MuonSystem->cscRechitCluster2NStation10[MuonSystem->nCscRechitClusters2] = tmp.nStation10;
        MuonSystem->cscRechitCluster2NStation10perc[MuonSystem->nCscRechitClusters2] = tmp.nStation10perc;
        MuonSystem->cscRechitCluster2AvgStation[MuonSystem->nCscRechitClusters2] = tmp.avgStation;
        MuonSystem->cscRechitCluster2AvgStation5[MuonSystem->nCscRechitClusters2] = tmp.avgStation5;
        MuonSystem->cscRechitCluster2AvgStation10[MuonSystem->nCscRechitClusters2] = tmp.avgStation10;
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
            MuonSystem->cscRechitCluster2_match_gLLP_e[MuonSystem->nCscRechitClusters2] = MuonSystem->gLLP_e[index];
            MuonSystem->cscRechitCluster2_match_gLLP_pt[MuonSystem->nCscRechitClusters2] = MuonSystem->gLLP_pt[index];

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
    }*/
    // cluster merging

    ds.merge_clusters();
    ds.result();
    ds.clusterMoments();
    ds.sort_clusters();




    MuonSystem->nCscRechitClusters3 = 0;
    for ( auto &tmp : ds.CscCluster ) {
        // if (tmp.tTotal>-12.5)continue;
        // if( tmp.nStation10 <2)continue;
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
        MuonSystem->cscRechitCluster3DeltaRSpread[MuonSystem->nCscRechitClusters3] =tmp.DeltaRSpread;


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
        MuonSystem->cscRechitCluster3NStation[MuonSystem->nCscRechitClusters3] = tmp.nStation;
        MuonSystem->cscRechitCluster3NStation5[MuonSystem->nCscRechitClusters3] = tmp.nStation5;
        MuonSystem->cscRechitCluster3NStation10[MuonSystem->nCscRechitClusters3] = tmp.nStation10;
        MuonSystem->cscRechitCluster3NStation10perc[MuonSystem->nCscRechitClusters3] = tmp.nStation10perc;
        MuonSystem->cscRechitCluster3AvgStation[MuonSystem->nCscRechitClusters3] = tmp.avgStation;
        MuonSystem->cscRechitCluster3AvgStation5[MuonSystem->nCscRechitClusters3] = tmp.avgStation5;
        MuonSystem->cscRechitCluster3AvgStation10[MuonSystem->nCscRechitClusters3] = tmp.avgStation10;
        MuonSystem->cscRechitCluster3AvgStation10perc[MuonSystem->nCscRechitClusters3] = tmp.avgStation10perc;
        MuonSystem->cscRechitCluster3Me11Ratio[MuonSystem->nCscRechitClusters3] = tmp.Me11Ratio;
        MuonSystem->cscRechitCluster3Me12Ratio[MuonSystem->nCscRechitClusters3] = tmp.Me12Ratio;
        // if(MuonSystem->category == 2)
        // {
        //   if (RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex1], MuonSystem->lepPhi[MuonSystem->ZleptonIndex1], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4) {
        //     MuonSystem->cscRechitCluster3ZLep1[MuonSystem->nCscRechitClusters3] = true;
        //     MuonSystem->cscRechitCluster3ZLep1Id[MuonSystem->nCscRechitClusters3] = MuonSystem->lepPdgId[MuonSystem->ZleptonIndex1];
        //     MuonSystem->cscRechitCluster3ZLep1LooseIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassLooseIso[MuonSystem->ZleptonIndex1];
        //     MuonSystem->cscRechitCluster3ZLep1TightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
        //     MuonSystem->cscRechitCluster3ZLep1VTightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassVTightIso[MuonSystem->ZleptonIndex1];
        //     MuonSystem->cscRechitCluster3ZLep1VVTightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassVVTightIso[MuonSystem->ZleptonIndex1];
        //     MuonSystem->cscRechitCluster3ZLep1TightId[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassId[MuonSystem->ZleptonIndex1];
        //
        //   }
        //   if (RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex2], MuonSystem->lepPhi[MuonSystem->ZleptonIndex2], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4) {
        //     MuonSystem->cscRechitCluster3ZLep2[MuonSystem->nCscRechitClusters3] = true;
        //     MuonSystem->cscRechitCluster3ZLep2Id[MuonSystem->nCscRechitClusters3] = MuonSystem->lepPdgId[MuonSystem->ZleptonIndex2];
        //     MuonSystem->cscRechitCluster3ZLep2LooseIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassLooseIso[MuonSystem->ZleptonIndex1];
        //     MuonSystem->cscRechitCluster3ZLep2TightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
        //     MuonSystem->cscRechitCluster3ZLep2VTightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassVTightIso[MuonSystem->ZleptonIndex1];
        //     MuonSystem->cscRechitCluster3ZLep2VVTightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassVVTightIso[MuonSystem->ZleptonIndex1];
        //     MuonSystem->cscRechitCluster3ZLep2TightId[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassId[MuonSystem->ZleptonIndex1];
        //   }
        // }
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
        // genjet veto
        for(int i = 0; i < nGenJets; i++)
        {
          if (fabs(genJetEta[i]>3.0)) continue;
          if (RazorAnalyzer::deltaR(genJetEta[i], genJetPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 && genJetPt[i] > MuonSystem->cscRechitCluster3GenJetVetoPt[MuonSystem->nCscRechitClusters3] ) {
            MuonSystem->cscRechitCluster3GenJetVetoPt[MuonSystem->nCscRechitClusters3]  = genJetPt[i];
          }
          if (RazorAnalyzer::deltaR(genJetEta[i], genJetPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 && genJetE[i] > MuonSystem->cscRechitCluster3GenJetVetoE[MuonSystem->nCscRechitClusters3] ) {
            MuonSystem->cscRechitCluster3GenJetVetoE[MuonSystem->nCscRechitClusters3]  = genJetE[i];
          }
        }
        float min_deltaR = 15.;
        int index = 999;


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
          // min_deltaR = 15.;
          // index = 999;
          // for(int j = 0; j < nGenParticle; j++)
          // {
          //
          //   double current_delta_r = RazorAnalyzer::deltaR(MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3], gParticleEta[j], gParticlePhi[j]);
          //
          //   if (current_delta_r < min_deltaR)
          //   {
          //     min_deltaR = current_delta_r;
          //     index = j;
          //   }
          // }
          // if (min_deltaR < 0.4)
          // {
          //
          //   MuonSystem->cscRechitCluster3_match_gParticle[MuonSystem->nCscRechitClusters3] = true;
          //   MuonSystem->cscRechitCluster3_match_gParticle_minDeltaR[MuonSystem->nCscRechitClusters3] = min_deltaR;
          //   MuonSystem->cscRechitCluster3_match_gParticle_index[MuonSystem->nCscRechitClusters3] = index;
          //   MuonSystem->cscRechitCluster3_match_gParticle_id[MuonSystem->nCscRechitClusters3] = gParticleId[index];
          //
          //   MuonSystem->cscRechitCluster3_match_gParticle_eta[MuonSystem->nCscRechitClusters3] = gParticleEta[index];
          //   MuonSystem->cscRechitCluster3_match_gParticle_phi[MuonSystem->nCscRechitClusters3] = gParticlePhi[index];
          //   MuonSystem->cscRechitCluster3_match_gParticle_E[MuonSystem->nCscRechitClusters3] = gParticleE[index];
          //   MuonSystem->cscRechitCluster3_match_gParticle_pt[MuonSystem->nCscRechitClusters3] = gParticlePt[index];
          //   MuonSystem->cscRechitCluster3_match_gParticle_MotherId[MuonSystem->nCscRechitClusters3] = gParticleMotherId[index];
          //
          //
          // }
          // for(int i = 0; i<4;i++)
          // {
          //   MuonSystem->cscRechitCluster3_match_gLLP_daughter0_deltaR[MuonSystem->nCscRechitClusters3] = deltaR(gLLP_daughter_eta[0],gLLP_daughter_phi[0],MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]);
          //   MuonSystem->cscRechitCluster3_match_gLLP_daughter1_deltaR[MuonSystem->nCscRechitClusters3] = deltaR(gLLP_daughter_eta[1],gLLP_daughter_phi[1],MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]);
          //   MuonSystem->cscRechitCluster3_match_gLLP_daughter2_deltaR[MuonSystem->nCscRechitClusters3] = deltaR(gLLP_daughter_eta[2],gLLP_daughter_phi[2],MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]);
          //   MuonSystem->cscRechitCluster3_match_gLLP_daughter3_deltaR[MuonSystem->nCscRechitClusters3] = deltaR(gLLP_daughter_eta[3],gLLP_daughter_phi[3],MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]);
          //
          //
          // }


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
            MuonSystem->cscRechitCluster3_match_gLLP_e[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_e[index];
            MuonSystem->cscRechitCluster3_match_gLLP_pt[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_pt[index];
            MuonSystem->cscRechitCluster3_match_gLLP_lepdPhi[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_lepdPhi[index];
            MuonSystem->cscRechitCluster3_match_gLLP_EMFracP[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_EMFracP[index];
            MuonSystem->cscRechitCluster3_match_gLLP_EMFracPz[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_EMFracPz[index];
            MuonSystem->cscRechitCluster3_match_gLLP_EMFracE[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_EMFracE[index];
            MuonSystem->cscRechitCluster3_match_gLLP_EMFracEz[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_EMFracEz[index];
            MuonSystem->cscRechitCluster3_match_gLLP_visP[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_visP[index];
            MuonSystem->cscRechitCluster3_match_gLLP_visPz[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_visPz[index];
            MuonSystem->cscRechitCluster3_match_gLLP_visE[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_visE[index];
            MuonSystem->cscRechitCluster3_match_gLLP_visEz[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_visEz[index];
            MuonSystem->cscRechitCluster3_match_gLLP_daughterKaon[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_daughterKaon[index];
            MuonSystem->cscRechitCluster3_match_gLLP_multiplicity[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_multiplicity[index];
            MuonSystem->cscRechitCluster3_match_gLLP_EM_multiplicity[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_EM_multiplicity[index];


            // MuonSystem->cscRechitCluster3_match_gLLP_other_eta[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_eta[1-index];
            // MuonSystem->cscRechitCluster3_match_gLLP_other_phi[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_phi[1-index];
            // MuonSystem->cscRechitCluster3_match_gLLP_other_decay_r[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_decay_vertex_r[1-index];
            // MuonSystem->cscRechitCluster3_match_gLLP_other_decay_x[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_decay_vertex_x[1-index];
            // MuonSystem->cscRechitCluster3_match_gLLP_other_decay_y[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_decay_vertex_y[1-index];
            // MuonSystem->cscRechitCluster3_match_gLLP_other_decay_z[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_decay_vertex_z[1-index];
            // MuonSystem->cscRechitCluster3_match_gLLP_other_ctau[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_ctau[1-index];
            // MuonSystem->cscRechitCluster3_match_gLLP_other_beta[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_beta[1-index];
            // MuonSystem->cscRechitCluster3_match_gLLP_other_csc[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_csc[1-index];
            // MuonSystem->cscRechitCluster3_match_gLLP_other_e[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_e[1-index];
            // MuonSystem->cscRechitCluster3_match_gLLP_other_pt[MuonSystem->nCscRechitClusters3] = MuonSystem->gLLP_pt[1-index];

            // if (index == 0)
            // {
            //
            //   if(MuonSystem->cscRechitCluster3_match_gLLP_daughter2_deltaR[MuonSystem->nCscRechitClusters3]<MuonSystem->cscRechitCluster3_match_gLLP_daughter3_deltaR[MuonSystem->nCscRechitClusters3])
            //   {
            //     MuonSystem->cscRechitCluster3_match_gLLP_other_daughter_deltaR[MuonSystem->nCscRechitClusters3] = MuonSystem->cscRechitCluster3_match_gLLP_daughter2_deltaR[MuonSystem->nCscRechitClusters3];
            //     MuonSystem->cscRechitCluster3_match_gLLP_other_daughter_index[MuonSystem->nCscRechitClusters3] = 2;
            //   }
            //   else
            //   {
            //     MuonSystem->cscRechitCluster3_match_gLLP_other_daughter_deltaR[MuonSystem->nCscRechitClusters3] = MuonSystem->cscRechitCluster3_match_gLLP_daughter3_deltaR[MuonSystem->nCscRechitClusters3];
            //     MuonSystem->cscRechitCluster3_match_gLLP_other_daughter_index[MuonSystem->nCscRechitClusters3] = 3;
            //   }
            // }
            // else
            // {
            //   if(MuonSystem->cscRechitCluster3_match_gLLP_daughter0_deltaR[MuonSystem->nCscRechitClusters3]<MuonSystem->cscRechitCluster3_match_gLLP_daughter1_deltaR[MuonSystem->nCscRechitClusters3])
            //   {
            //     MuonSystem->cscRechitCluster3_match_gLLP_other_daughter_deltaR[MuonSystem->nCscRechitClusters3] = MuonSystem->cscRechitCluster3_match_gLLP_daughter0_deltaR[MuonSystem->nCscRechitClusters3];
            //     MuonSystem->cscRechitCluster3_match_gLLP_other_daughter_index[MuonSystem->nCscRechitClusters3] = 0;
            //   }
            //   else
            //   {
            //     MuonSystem->cscRechitCluster3_match_gLLP_other_daughter_deltaR[MuonSystem->nCscRechitClusters3] = MuonSystem->cscRechitCluster3_match_gLLP_daughter1_deltaR[MuonSystem->nCscRechitClusters3];
            //     MuonSystem->cscRechitCluster3_match_gLLP_other_daughter_index[MuonSystem->nCscRechitClusters3] = 1;
            //   }
            // }


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
          //match to check for beam halo
          if (RazorAnalyzer::deltaPhi(dtRechitPhi[i], MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.2 )
          {
            MuonSystem->cscRechitCluster3_match_dtRechits_phi0p2[MuonSystem->nCscRechitClusters3] ++;
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
        MuonSystem->nCscRechitClusters3++;
    }



    //  for (unsigned int j = 0; j < tmp.segment_id.size(); j++)
    // {
    //   MuonSystem->cscRechitsCluster2Id[tmp.segment_id[j]] = MuonSystem->nCscRechitClusters2;
    // }
    // if(isData&& MuonSystem->nCscRechitClusters3==0 && MuonSystem->nCscRechitClusters2==0 && MuonSystem->nCscRechitClusters==0) continue;
    // cout<<"fill"<<endl;
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
