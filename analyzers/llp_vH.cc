#include "llp_vH.h"
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
  bool matched;
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

/*
class RazorLiteTree
{

public:
  UInt_t  runNum,
//pt comparison
//not used so fa
////pt comparison
//not used so farr lumiSec, evtNum;
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
  float jetEt[N_MAX_JETS];
  float jetPt[N_MAX_JETS];
  float jetEta[N_MAX_JETS];
  float jetPhi[N_MAX_JETS];
  float jetTime[N_MAX_JETS];
  float ecalNRechits[N_MAX_JETS];
  float ecalRechitE[N_MAX_JETS];
  float jetChargedEMEnergyFraction[N_MAX_JETS];
  float jetNeutralEMEnergyFraction[N_MAX_JETS];
  float jetChargedHadronEnergyFraction[N_MAX_JETS];
  float jetNeutralHadronEnergyFraction[N_MAX_JETS];
  float jetGammaMax_ET[N_MAX_JETS];
  float jetMinDeltaRPVTracks[N_MAX_JETS];
  float jetPtAllPVTracks[N_MAX_JETS];
  // bool jetLoosePassId[N_MAX_JETS];
  bool jetPassId[N_MAX_JETS];
  bool matched[N_MAX_JETS];
  // bool jetTightPassId[N_MAX_JETS];
  float jet_energy_frac[N_MAX_JETS];
  float jet_sig_et1[N_MAX_JETS];
  float jet_sig_et2[N_MAX_JETS];

  //calojet
  int nCaloJets;
  float calojetE[N_MAX_JETS];
  float calojetEt[N_MAX_JETS];
  float calojetPt[N_MAX_JETS];
  float calojetEta[N_MAX_JETS];
  float calojetPhi[N_MAX_JETS];
  float calojetTime[N_MAX_JETS];
  float calojetNRechits[N_MAX_JETS];
  float calojetRechitE[N_MAX_JETS];
  float calojetChargedEMEnergyFraction[N_MAX_JETS];
  float calojetNeutralEMEnergyFraction[N_MAX_JETS];
  float calojetChargedHadronEnergyFraction[N_MAX_JETS];
  float calojetNeutralHadronEnergyFraction[N_MAX_JETS];
  float calojetGammaMax_ET[N_MAX_JETS];
  float calojetMinDeltaRPVTracks[N_MAX_JETS];
  float calojetPtAllPVTracks[N_MAX_JETS];
  // bool jetLoosePassId[N_MAX_JETS];
  bool calojetPassId[N_MAX_JETS];
  // // bool jetTightPassId[N_MAX_JETS];

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
    nCaloJets = 0;
    for( int i = 0; i < N_MAX_JETS; i++ )
    {
      jetE[i]      = -999.;
      jetEt[i] = -999.;
      jetPt[i]     = -999.;
      jetEta[i]    = -999.;
      jetPhi[i]    = -999.;
      jetTime[i]   = -999.;
      // jetLoosePassId[i] = false;
      jetPassId[i] = false;
      matched[i] = false;
      jet_energy_frac[i] = -999.;
      jet_sig_et1[i] = -999.;
      jet_sig_et2[i] = -999.;
      ecalNRechits[i] = -999;
      ecalRechitE[i] = -999.;
      jetChargedEMEnergyFraction[i] = -999.;
      jetNeutralEMEnergyFraction[i] = -999.;
      jetChargedHadronEnergyFraction[i] = -999.;
      jetNeutralHadronEnergyFraction[i] = -999.;
      jetGammaMax_ET[i] = -999.;
      jetMinDeltaRPVTracks[i] = -999.;
      jetPtAllPVTracks[i] = -999.;

      calojetE[i] = -999.;
      calojetEt[i] = -999.;
      calojetPt[i] = -999.;
      calojetEta[i] = -999.;
      calojetPhi[i] = -999.;
      calojetTime[i] = -999.;
      calojetPassId[i] = -999.;
      calojetGammaMax_ET[i] = -999.;
      calojetMinDeltaRPVTracks[i] = -999.;
      calojetPtAllPVTracks[i] = -999.;
      calojetNRechits[i] = -999.;
      calojetRechitE[i] = -999.;
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
    tree_->Branch("jetEt",      jetEt,      "jetEt[nJets]/F");
    tree_->Branch("jetPt",     jetPt,     "jetPt[nJets]/F");
    tree_->Branch("jetEta",    jetEta,    "jetEta[nJets]/F");
    tree_->Branch("jetPhi",    jetPhi,    "jetPhi[nJets]/F");
    tree_->Branch("jetTime",   jetTime,   "jetTime[nJets]/F");
    tree_->Branch("jetPassId", jetPassId, "jetPassId[nJets]/O");
    tree_->Branch("matched", matched, "matched[nJets]/O");
    tree_->Branch("jet_energy_frac", jet_energy_frac, "jet_energy_frac[nJets]/F");
    tree_->Branch("jet_sig_et1", jet_sig_et1, "jet_sig_et1[nJets]/F");
    tree_->Branch("jet_sig_et2", jet_sig_et2, "jet_sig_et2[nJets]/F");


    tree_->Branch("jetGammaMax_ET", jetGammaMax_ET, "jetGammaMax_ET[nJets]/F");
    tree_->Branch("jetMinDeltaRPVTracks", jetMinDeltaRPVTracks, "jetMinDeltaRPVTracks[nJets]/F");
    tree_->Branch("jetPtAllPVTracks", jetPtAllPVTracks, "jetPtAllPVTracks[nJets]/F");
    tree_->Branch("ecalNRechits",   ecalNRechits,   "ecalNRechits[nJets]/F");
    tree_->Branch("ecalRechitE", ecalRechitE, "ecalRechitE[nJets]/F");
    tree_->Branch("jetChargedEMEnergyFraction",   jetChargedEMEnergyFraction,   "jetChargedEMEnergyFraction[nJets]/F");
    tree_->Branch("jetNeutralEMEnergyFraction",   jetNeutralEMEnergyFraction,   "jetNeutralEMEnergyFraction[nJets]/F");
    tree_->Branch("jetChargedHadronEnergyFraction",   jetChargedHadronEnergyFraction,   "jetChargedHadronEnergyFraction[nJets]/F");
    tree_->Branch("jetNeutralHadronEnergyFraction",   jetNeutralHadronEnergyFraction,   "jetNeutralHadronEnergyFraction[nJets]/F");


      //calojet
    tree_->Branch("nCaloJets",     &nCaloJets,    "nCaloJets/I");
    tree_->Branch("calojetE",      calojetE,      "calojetE[nCaloJets]/F");
    tree_->Branch("calojetEt",      calojetEt,      "calojetEt[nCaloJets]/F");
    tree_->Branch("calojetPt",     calojetPt,     "calojetPt[nCaloJets]/F");
    tree_->Branch("calojetEta",    calojetEta,    "calojetEta[nCaloJets]/F");
    tree_->Branch("calojetPhi",    calojetPhi,    "calojetPhi[nCaloJets]/F");
    tree_->Branch("calojetTime",   calojetTime,   "calojetTime[nCaloJets]/F");
    tree_->Branch("calojetPassId", calojetPassId, "calojetPassId[nCaloJets]/O");
    tree_->Branch("calojetGammaMax_ET", calojetGammaMax_ET, "calojetGammaMax_ET[nCaloJets]/F");
    tree_->Branch("calojetMinDeltaRPVTracks", calojetMinDeltaRPVTracks, "calojetMinDeltaRPVTracks[nCaloJets]/F");
    tree_->Branch("calojetPtAllPVTracks", calojetPtAllPVTracks, "calojetPtAllPVTracks[nCaloJets]/F");
    tree_->Branch("calojetNRechits",   calojetNRechits,   "calojetNRechits[nCaloJets]/F");
    tree_->Branch("calojetRechitE", calojetRechitE, "calojetRechitE[nCaloJets]/F");




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
    tree_->SetBranchAddress("jetEt",      jetEt);

    tree_->SetBranchAddress("jetPt",     jetPt);
    tree_->SetBranchAddress("jetEta",    jetEta);
    tree_->SetBranchAddress("jetPhi",    jetPhi);
    tree_->SetBranchAddress("jetTime",   jetTime);
    tree_->SetBranchAddress("jetPassId", jetPassId);
    tree_->SetBranchAddress("matched", matched);
    tree_->SetBranchAddress("jet_energy_frac", jet_energy_frac);
    tree_->SetBranchAddress("jet_sig_et1", jet_sig_et1);
    tree_->SetBranchAddress("jet_sig_et2", jet_sig_et2);

    tree_->SetBranchAddress("ecalNRechits",   ecalNRechits);
    tree_->SetBranchAddress("ecalRechitE", ecalRechitE);
    tree_->SetBranchAddress("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction);
    tree_->SetBranchAddress("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction);
    tree_->SetBranchAddress("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction);
    tree_->SetBranchAddress("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction);
    tree_->SetBranchAddress("jetGammaMax_ET", jetGammaMax_ET);
    tree_->SetBranchAddress("jetMinDeltaRPVTracks", jetMinDeltaRPVTracks);
    tree_->SetBranchAddress("jetPtAllPVTracks", jetPtAllPVTracks);
    //calojets
    tree_->SetBranchAddress("nCaloJets",     &nCaloJets);
    tree_->SetBranchAddress("calojetE",      calojetE);
    tree_->SetBranchAddress("calojetEt",      calojetEt);
    tree_->SetBranchAddress("calojetPt",     calojetPt);
    tree_->SetBranchAddress("calojetEta",    calojetEta);
    tree_->SetBranchAddress("calojetPhi",    calojetPhi);
    tree_->SetBranchAddress("calojetTime",   calojetTime);
    tree_->SetBranchAddress("calojetPassId", calojetPassId);
    tree_->SetBranchAddress("calojetRechitE",   calojetRechitE);
    // tree_->SetBranchAddress("calojetRechitT_rms", calojetRechitT_rms);
    // tree_->SetBranchAddress("calojet_EMEnergyFraction", calojet_EMEnergyFraction);
    // tree_->SetBranchAddress("calojet_HadronicEnergyFraction", calojet_HadronicEnergyFraction);

    tree_->SetBranchAddress("calojetGammaMax_ET", calojetGammaMax_ET);
    tree_->SetBranchAddress("calojetMinDeltaRPVTracks", calojetMinDeltaRPVTracks);
    tree_->SetBranchAddress("calojetPtAllPVTracks", calojetPtAllPVTracks);
    // tree_->SetBranchAddress("jetLoosePassId", jetLoosePassId);
    // tree_->SetBranchAddress("jetTightPassId", jetTightPassId);
    // triggers
    tree_->SetBranchAddress("HLTDecision",   HLTDecision);

  };

};
*/
void llp_vH::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  int option;
  std::string label;
  bool pf;
  if (options < 200){
    option = 1; // used when running condor
  }
  else{
    option = 0;// used when running locally
  }
  if ((options/10)%10 == 1){
    label = "wH";
  }
  else if ((options/10) % 10 == 2){
    label = "zH";
  }
  else if ((options/10) % 10 == 3){
    label = "bkg_wH";
  }
  else{
    label = "bkg_zH";
  }
  if(options%10==1){
    pf = true;
  }
  else{
    pf = false;
  }
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
  if (outfilename == "") outfilename = "vH_Tree.root";
  TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
  //RazorLiteTree *vH = new RazorLiteTree;
  SusyLLPTree *vH = new SusyLLPTree;
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
    if (label =="bkg_wH"|| label == "bkg_zH"){
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
      NEvents->Fill(1);
      bool wzFlag = false;
      for (int i=0; i < nGenParticle; ++i)
      {
        // if (abs(gParticleId[i]) == wzId && gParticleStatus[i] == 22)
        // {
        if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == wzId)
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

    if ( Leptons.size() < nLepton_cut ) continue;
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
      tmpJet.matched = jet_matched[i];
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

      Jets.push_back(tmpJet);

    }
    std::vector<jets> caloJets;
    //auto highest = [](auto a, auto b) { return a > b; };

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


        jets tmpJet;
        tmpJet.jet    = thisJet;
        tmpJet.time   = calojetRechitT[i];
        tmpJet.passId = calojetPassIDLoose[i];
        tmpJet.isCSVL = isCSVL(i);
        //if (isCSVL(i)) NBJet20++;
        //if (isCSVL(i) && thisJet.Pt() > 30) NBJet30++;
        tmpJet.ecalNRechits = calojetNRechits[i];
        tmpJet.ecalRechitE = calojetRechitE[i];

        tmpJet.jetGammaMax_ET = calojetGammaMax_ET[i];
        tmpJet.jetMinDeltaRPVTracks = calojetMinDeltaRPVTracks[i];
        tmpJet.jetPtAllPVTracks = calojetPtAllPVTracks[i];

        caloJets.push_back(tmpJet);

      }
    //-----------------------------
    //Require at least 2 jets
    //-----------------------------
    if(pf)
    {
      if( Jets.size() < 1 ) continue;

    }
    else
    {
      if( caloJets.size() < 1 ) continue;

    }
    if (triggered) trig_lepId_dijet->Fill(1);
    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

    for ( auto &tmp : Jets )
    {
      vH->jetE[vH->nJets] = tmp.jet.E();
      vH->jetEt[vH->nJets] = tmp.jet.Et();
      vH->jetPt[vH->nJets] = tmp.jet.Pt();
      vH->jetEta[vH->nJets] = tmp.jet.Eta();
      vH->jetPhi[vH->nJets] = tmp.jet.Phi();
      vH->jetTime[vH->nJets] = tmp.time;
      vH->jetPassId[vH->nJets] = tmp.passId;
      vH->matched[vH->nJets] = tmp.matched;
      vH->jet_sig_et1[vH->nJets] = tmp.sig_et1;
      vH->jet_sig_et2[vH->nJets] = tmp.sig_et2;
      vH->jet_energy_frac[vH->nJets] = tmp.energy_frac;
      vH->ecalNRechits[vH->nJets] = tmp.ecalNRechits;
      vH->ecalRechitE[vH->nJets] = tmp.ecalRechitE;
      vH->jetChargedEMEnergyFraction[vH->nJets] = tmp.jetChargedEMEnergyFraction;
      vH->jetNeutralEMEnergyFraction[vH->nJets] = tmp.jetNeutralEMEnergyFraction;
      vH->jetChargedHadronEnergyFraction[vH->nJets] = tmp.jetChargedHadronEnergyFraction;
      vH->jetNeutralHadronEnergyFraction[vH->nJets] = tmp.jetNeutralHadronEnergyFraction;
      vH->jetGammaMax_ET[vH->nJets] = tmp.jetGammaMax_ET;
      vH->jetMinDeltaRPVTracks[vH->nJets] = tmp.jetMinDeltaRPVTracks;
      vH->jetPtAllPVTracks[vH->nJets] = tmp.jetPtAllPVTracks;

      // std::cout <<tmp.time << "," <<tmp.ecalRechitE <<  "," << tmp.ecalNRechits << vH->nJets<<std::endl;

      vH->nJets++;
    }
    sort(caloJets.begin(), caloJets.end(), my_largest_pt_jet);

    for ( auto &tmp : caloJets )
    {
      vH->calojetE[vH->nCaloJets] = tmp.jet.E();
      vH->calojetEt[vH->nCaloJets] = tmp.jet.Et();
      vH->calojetPt[vH->nCaloJets] = tmp.jet.Pt();
      vH->calojetEta[vH->nCaloJets] = tmp.jet.Eta();
      vH->calojetPhi[vH->nCaloJets] = tmp.jet.Phi();
      vH->calojetTime[vH->nCaloJets] = tmp.time;
      vH->calojetPassId[vH->nCaloJets] = tmp.passId;
      vH->calojetNRechits[vH->nCaloJets] = tmp.ecalNRechits;
      vH->calojetRechitE[vH->nCaloJets] = tmp.ecalRechitE;

      vH->calojetGammaMax_ET[vH->nCaloJets] = tmp.jetGammaMax_ET;
      vH->calojetMinDeltaRPVTracks[vH->nCaloJets] = tmp.jetMinDeltaRPVTracks;
      vH->calojetPtAllPVTracks[vH->nCaloJets] = tmp.jetPtAllPVTracks;

      // std::cout <<tmp.time << "," <<tmp.ecalRechitE <<  "," << tmp.ecalNRechits << vH->nJets<<std::endl;

      vH->nCaloJets++;
    }
    //std::cout << "deb fill: " << vH->nLeptons << " " << jentry << endl;
    vH->tree_->Fill();
  }

    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();
}
