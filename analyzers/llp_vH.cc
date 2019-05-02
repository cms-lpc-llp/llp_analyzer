#include "llp_vH.h"
#include "RazorHelper.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"

//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define NTriggersMAX 601 //Number of trigger in the .dat file
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


class RazorLiteTree
{

public:
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t rho, weight;
  Float_t met, metPhi;
  //leptons
  int nLeptons;
  float lepE[N_MAX_LEPTONS];
  float lepPt[N_MAX_LEPTONS];
  float lepEta[N_MAX_LEPTONS];
  float lepPhi[N_MAX_LEPTONS];
  int  lepPdgId[N_MAX_LEPTONS];
  float lepDZ[N_MAX_LEPTONS];
  // bool lepLoosePassId[N_MAX_LEPTONS];
  // bool lepMediumPassId[N_MAX_LEPTONS];
  // bool lepTightPassId[N_MAX_LEPTONS];
  bool lepPassId[N_MAX_LEPTONS];

  //Z-candidate
  float MT;
  float ZMass;
  float ZPt;
  float ZEta;
  float ZPhi;
  int ZleptonIndex1;
  int ZleptonIndex2;
  //jets
  int nJets;
  float jetE[N_MAX_JETS];
  float jetPt[N_MAX_JETS];
  float jetEta[N_MAX_JETS];
  float jetPhi[N_MAX_JETS];
  float jetTime[N_MAX_JETS];
  // bool jetLoosePassId[N_MAX_JETS];
  bool jetPassId[N_MAX_JETS];
  // bool jetTightPassId[N_MAX_JETS];
  bool HLTDecision[NTriggersMAX];

  UInt_t wzevtNum,trig, trig_lepId, trig_lepId_dijet; //number of events that pass each criteria



  TTree *tree_;
  TFile *f_;

  RazorLiteTree()
  {
    InitVariables();
  };

  ~RazorLiteTree()
  {
    if (f_) f_->Close();
  };

  void InitVariables()
  {
    runNum=0; lumiSec=0; evtNum=0; category=0;
    npv=0; npu=0; rho=-1; weight=-1;
    met=-1; metPhi=-1;

    //leptons
    nLeptons = 0;
    for( int i = 0; i < N_MAX_LEPTONS; i++ )
    {
      lepE[i]      = -999.;
      lepPt[i]     = -999.;
      lepEta[i]    = -999.;
      lepPhi[i]    = -999.;
      lepPdgId[i]  = -999;
      lepDZ[i]     = -999.;
      // lepLoosePassId[i] = false;
      // lepMediumPassId[i] = false;
      // lepTightPassId[i] = false;
      lepPassId[i] = false;
    }
    //Z-candidate
    ZMass = -999.; ZPt = -999.; ZEta = -999.; ZPhi = -999.;
    MT = -999.;
    ZleptonIndex1 = -999; ZleptonIndex2 = -999;
    //jets
    nJets = 0;
    for( int i = 0; i < N_MAX_JETS; i++ )
    {
      jetE[i]      = -999.;
      jetPt[i]     = -999.;
      jetEta[i]    = -999.;
      jetPhi[i]    = -999.;
      jetTime[i]   = -999.;
      // jetLoosePassId[i] = false;
      jetPassId[i] = false;
      // jetTightPassId[i] = false;
    }

    for(int i = 0; i <NTriggersMAX; i++){
      HLTDecision[i] = false;
    }

  };

  void LoadTree(const char* file)
  {
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->Get("vH"));
    InitTree();
    assert(tree_);
  };

  void CreateTree()
  {
    tree_ = new TTree("vH","vH");
    f_ = 0;

    tree_->Branch("runNum",      &runNum,     "runNum/i");      // event run number
    tree_->Branch("lumiSec",     &lumiSec,    "lumiSec/i");     // event lumi section
    tree_->Branch("evtNum",      &evtNum,     "evtNum/i");      // event number
    tree_->Branch("category",    &category,   "category/i");    // dilepton category
    tree_->Branch("npv",         &npv,        "npv/i");         // number of primary vertices
    tree_->Branch("npu",         &npu,        "npu/i");         // number of in-time PU events (MC)
    tree_->Branch("weight",      &weight,     "weight/F");
    tree_->Branch("rho",         &rho,        "rho/F");
    tree_->Branch("met",         &met,        "met/F");         // MET
    tree_->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)
    //leptons
    tree_->Branch("nLeptons",  &nLeptons, "nLeptons/I");
    tree_->Branch("lepE",      lepE,      "lepE[nLeptons]/F");
    tree_->Branch("lepPt",     lepPt,     "lepPt[nLeptons]/F");
    tree_->Branch("lepEta",    lepEta,    "lepEta[nLeptons]/F");
    tree_->Branch("lepPhi",    lepPhi,    "lepPhi[nLeptons]/F");
    tree_->Branch("lepPdgId",  lepPdgId,  "lepPdgId[nLeptons]/I");
    tree_->Branch("lepDZ",     lepDZ,     "lepDZ[nLeptons]/F");
    tree_->Branch("lepPassId", lepPassId, "lepPassId[nLeptons]/O");
    // tree_->Branch("lepLoosePassId", lepLoosePassId, "lepLoosePassId[nLeptons]/O");
    // tree_->Branch("lepMediumPassId", lepMediumPassId, "lepMediumPassId[nLeptons]/O");
    // tree_->Branch("lepTightPassId", lepTightPassId, "lepTightPassId[nLeptons]/O");

    //Z-candidate
    tree_->Branch("MT",      &MT,  "MT/F");
    tree_->Branch("ZMass",      &ZMass,  "ZMass/F");
    tree_->Branch("ZPt",        &ZPt,    "ZPt/F");
    tree_->Branch("ZEta",       &ZEta,   "ZEta/F");
    tree_->Branch("ZPhi",       &ZPhi,   "ZPhi/F");
    tree_->Branch("ZleptonIndex1", &ZleptonIndex1, "ZleptonIndex1/I");
    tree_->Branch("ZleptonIndex2", &ZleptonIndex2, "ZleptonIndex2/I");
    //jets
    tree_->Branch("nJets",     &nJets,    "nJets/I");
    tree_->Branch("jetE",      jetE,      "jetE[nJets]/F");
    tree_->Branch("jetPt",     jetPt,     "jetPt[nJets]/F");
    tree_->Branch("jetEta",    jetEta,    "jetEta[nJets]/F");
    tree_->Branch("jetPhi",    jetPhi,    "jetPhi[nJets]/F");
    tree_->Branch("jetTime",   jetTime,   "jetTime[nJets]/F");
    tree_->Branch("jetPassId", jetPassId, "jetPassId[nJets]/O");
    // tree_->Branch("jetLoosePassId", jetLoosePassId, "jetLoosePassId[nJets]/O");
    // tree_->Branch("jetTightPassId", jetTightPassId, "jetTightPassId[nJets]/O");
    tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[601]/O"); //hardcoded



  };

  void InitTree()
  {
    assert(tree_);
    InitVariables();

    tree_->SetBranchAddress("runNum",      &runNum);
    tree_->SetBranchAddress("lumiSec",     &lumiSec);
    tree_->SetBranchAddress("evtNum",      &evtNum);
    tree_->SetBranchAddress("category",    &category);
    tree_->SetBranchAddress("npv",         &npv);
    tree_->SetBranchAddress("npu",         &npu);
    tree_->SetBranchAddress("weight",      &weight);
    tree_->SetBranchAddress("rho",         &rho);
    tree_->SetBranchAddress("met",         &met);
    tree_->SetBranchAddress("metPhi",      &metPhi);
    //Leptons
    tree_->SetBranchAddress("nLeptons",    &nLeptons);
    tree_->SetBranchAddress("lepE",        lepE);
    tree_->SetBranchAddress("lepPt",       lepPt);
    tree_->SetBranchAddress("lepEta",      lepEta);
    tree_->SetBranchAddress("lepPhi",      lepPhi);
    tree_->SetBranchAddress("lepPdgId",  lepPdgId);
    tree_->SetBranchAddress("lepDZ",     lepDZ);
    // tree_->SetBranchAddress("lepLoosePassId", lepLoosePassId);
    // tree_->SetBranchAddress("lepMediumPassId", lepMediumPassId);
    // tree_->SetBranchAddress("lepTightPassId", lepTightPassId);
    tree_->SetBranchAddress("lepPassId", lepPassId);

    //Z-candidate
    tree_->SetBranchAddress("ZMass",       &ZMass);
    tree_->SetBranchAddress("ZPt",         &ZPt);
    tree_->SetBranchAddress("ZEta",        &ZEta);
    tree_->SetBranchAddress("ZPhi",        &ZPhi);
    tree_->SetBranchAddress("ZleptonIndex1", &ZleptonIndex1);
    tree_->SetBranchAddress("ZleptonIndex2", &ZleptonIndex2);
    tree_->SetBranchAddress("MT", &MT);

    //jets
    tree_->SetBranchAddress("nJets",     &nJets);
    tree_->SetBranchAddress("jetE",      jetE);
    tree_->SetBranchAddress("jetPt",     jetPt);
    tree_->SetBranchAddress("jetEta",    jetEta);
    tree_->SetBranchAddress("jetPhi",    jetPhi);
    tree_->SetBranchAddress("jetTime",   jetTime);
    tree_->SetBranchAddress("jetPassId", jetPassId);
    // tree_->SetBranchAddress("jetLoosePassId", jetLoosePassId);
    // tree_->SetBranchAddress("jetTightPassId", jetTightPassId);
    // triggers
    tree_->SetBranchAddress("HLTDecision",   HLTDecision);

  };

};

void llp_vH::Analyze(bool isData, int options, string outputfilename, string analysisTag)
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
  else{
    label = "bkg";
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
  int NTrigger;;//Number of trigger in trigger paths





  if (label == "zH"){
    NTrigger = 4;
    }
  else{
    NTrigger = 2;
  }

  int trigger_paths[NTrigger];
  if (label == "wH" || label == "bkg"){
    wzId = 24;
    trigger_paths[0] = 87;
    trigger_paths[1] = 135;
    // trigger_paths[2] = 310;
  }
  else if (label == "zH"){
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
  if (outfilename == "") outfilename = "vH_Tree.root";
  TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
  RazorLiteTree *vH = new RazorLiteTree;
  vH->CreateTree();
  vH->tree_->SetAutoFlush(0);
  vH->InitTree();
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
    //std::cout << "deb0 " << jentry << std::endl;
    vH->InitVariables();
    //std::cout << "deb1 " << jentry << std::endl;
    if (label =="bkg"){
      if (isData)
      {
        NEvents->Fill(1);
        vH->weight = 1;
      }
      else
      {
        //NEvents->Fill(genWeight);
        //vH->weight = genWeight;
        NEvents->Fill(1);
        vH->weight = 1;
      }

    }
    else{
      generatedEvents->Fill(1);
      vH->weight = 1;
    }
    //std::cout << "deb2 " << jentry << std::endl;
    //event info
    vH->runNum = runNum;
    vH->lumiSec = lumiNum;
    vH->evtNum = eventNum;
    //std::cout << "deb3 " << jentry << std::endl;
    if (label == "zH" || label == "wH"){
      bool wzFlag = false;
      for (int i=0; i < nGenParticle; ++i)
      {
        // if (abs(gParticleId[i]) == wzId && gParticleStatus[i] == 22)
        // {
        if (abs(gParticleId[i]) == 13 && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == wzId)
        { // choosing only the W->munu events
          wzFlag = true;
        }

      }
      if ( wzFlag == false ) continue;
      NEvents->Fill(1);
    }

    for (int i=0; i < nBunchXing; ++i)
    {
      if (BunchXing[i] == 0)
      {
        vH->npu = nPUmean[i];
      }
    }
    //get NPU
    vH->npv = nPV;
    vH->rho = fixedGridRhoFastjetAll;
    vH->met = metType1Pt;
    vH->metPhi = metType1Phi;

    //Triggers
    for(int i = 0; i < NTriggersMAX; i++){
      vH->HLTDecision[i] = HLTDecision[i];
    }
    bool triggered = false;
    for(int i = 0; i < NTrigger; i++)
    {
      int trigger_temp = trigger_paths[i];

      triggered = triggered || HLTDecision[trigger_temp];

    }
    if (triggered) trig->Fill(1);
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
      if(muonPt[i] < 15) continue;
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
    for( int i = 0; i < nElectrons; i++ )
    {


      // if(!(passMVAVetoElectronID(i) &&
      // ( (fabs(eleEta[i]) < 1.5 && fabs(ele_d0[i]) < 0.0564) ||
      // (fabs(eleEta[i]) >= 1.5 && fabs(ele_d0[i]) < 0.222))
      // && passEGammaPOGVetoElectronIso(i))) continue;
      if (!isEGammaPOGVetoElectron(i, true, true, true, "Summer16")) continue;

      if (!( (fabs(eleEta[i]) < 1.5 && fabs(ele_d0[i]) < 0.0564) ||
      (fabs(eleEta[i]) >= 1.5 && fabs(ele_d0[i]) < 0.222))) continue;

      if(elePt[i] < 15) continue;

      if(fabs(eleEta[i]) > 2.4) continue;

      //remove overlaps
      bool overlap = false;
      for(auto& lep : Leptons)
      {
        if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
      }
      if(overlap) continue;
      std::cout << "here" << std::endl;
      leptons tmpElectron;
      tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);
      tmpElectron.pdgId = 11 * -1 * eleCharge[i];
      tmpElectron.dZ = ele_dZ[i];
      // tmpElectron.passId = passMVALooseElectronID(i) && passEGammaPOGLooseElectronIso(i);
      Leptons.push_back(tmpElectron);
    }

    sort(Leptons.begin(), Leptons.end(), my_largest_pt);
    //std::cout << "deb7 " << jentry << std::endl;
    for ( auto &tmp : Leptons )
    {
      vH->lepE[vH->nLeptons]      = tmp.lepton.E();
      vH->lepPt[vH->nLeptons]     = tmp.lepton.Pt();
      vH->lepEta[vH->nLeptons]    = tmp.lepton.Eta();
      vH->lepPhi[vH->nLeptons]    = tmp.lepton.Phi();
      vH->lepPdgId[vH->nLeptons]  = tmp.pdgId;
      vH->lepDZ[vH->nLeptons]     = tmp.dZ;
      vH->lepPassId[vH->nLeptons] = tmp.passId;


      // std::cout << "lepton pdg " << vH->lepPdgId[vH->nLeptons] << std::endl;
      vH->nLeptons++;
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
      vH->ZMass = ZMass;
      vH->ZPt   = ZPt;
      vH->ZEta  = ZCandidate.Eta();
      vH->ZPhi  = ZCandidate.Phi();
      vH->ZleptonIndex1 = ZCandidateLeptonIndex.first;
      vH->ZleptonIndex2 = ZCandidateLeptonIndex.second;

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
    //   std::cout <<nMuons << "," << nElectrons <<  "," << vH->nLeptons <<  "," << vH->met << std::endl;
    // }

    if ( Leptons.size() < 1 ) continue;
    TLorentzVector met;
    TLorentzVector visible = Leptons[0].lepton;
    met.SetPtEtaPhiE(metType1Pt,0,metType1Phi,metType1Pt);
    vH->MT = GetMT(visible,met);

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
      if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
      if ( !jetPassIDLoose[i] ) continue;

      jets tmpJet;
      tmpJet.jet    = thisJet;
      tmpJet.time   = jetRechitT[i];
      tmpJet.passId = jetPassIDTight[i];
      tmpJet.isCSVL = isCSVL(i);
      Jets.push_back(tmpJet);
      //if (isCSVL(i)) NBJet20++;
      //if (isCSVL(i) && thisJet.Pt() > 30) NBJet30++;

    }

    //-----------------------------
    //Require at least 2 jets
    //-----------------------------
    // if( Jets.size() < 2 ) continue;
    if (triggered) trig_lepId_dijet->Fill(1);
    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

    for ( auto &tmp : Jets )
    {
      vH->jetE[vH->nJets] = tmp.jet.E();
      vH->jetPt[vH->nJets] = tmp.jet.Pt();
      vH->jetEta[vH->nJets] = tmp.jet.Eta();
      vH->jetPhi[vH->nJets] = tmp.jet.Phi();
      vH->jetTime[vH->nJets] = tmp.time;
      vH->jetPassId[vH->nJets] = tmp.passId;

      vH->nJets++;
    }
    //std::cout << "deb fill: " << vH->nLeptons << " " << jentry << endl;
    vH->tree_->Fill();
  }

    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();
}
