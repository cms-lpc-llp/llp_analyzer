#include "llp_MuonSystem.h"
#include "RazorHelper.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
#include "DBSCAN.h"
//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 2000
#define NTriggersMAX 601 //Number of trigger in the .dat file
#define N_CSC_CUT 20
#define JET_PT_CUT 10
using namespace std;

struct greater_than_pt
{
  inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
};
struct cscCluster
{
  float x;
  float y;
  float z;
  float eta;
  float phi;
  int nCscSegments;
  bool jetVeto;
  bool muonVeto;
  int station;
  float radius;
  float MajorAxis;
  float MinorAxis;
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

struct largest_nCsc
{
  inline bool operator() (const cscCluster& c1, const cscCluster& c2){return c1.nCscSegments > c2.nCscSegments;}
} my_largest_nCsc;

int cscStation(double x, double y, double z)
{
  double r = sqrt(x*x+y*y);
  z = abs(z);
  if (r > 100 && r < 263 && z > 568  && z < 632) return 11;
  if (r > 275 && r < 465.0 && z > 668.3 && z < 724) return 12;
  if (r > 505.5 && r < 695.5 && z > 686 && z < 724) return 13;
  if (r > 138.5 && r < 345 && z > 791 && z < 849.5) return 21;
  if (r > 357.5 && r < 695.5 && z > 791 && z < 849.5) return 22;
  if (r > 160.5 && r < 345 && z > 911.5 && z < 970) return 31;
  if (r > 357.5 && r < 695.5 && z > 911.5 && z < 970) return 32;
  if (r > 177.5 && r < 345 && z > 1002 && z < 1060.5) return 41;
  if (r > 357.5 && r < 695.5 && z > 1002 && z < 1060.5) return 42;
  return -999;
};
class LiteTreeMuonSystem
{

public:
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  category;
  UInt_t  npv, npu;
  float rho, weight;
  float met, metPhi;
  int lepGenId;

  //csc
  int           nCsc;
  int           cscLabels[N_MAX_CSC];
  float         cscPhi[N_MAX_CSC];   //[nCsc]
  float         cscEta[N_MAX_CSC];   //[nCsc]
  float         cscX[N_MAX_CSC];   //[nCsc]
  float         cscY[N_MAX_CSC];   //[nCsc]
  float         cscZ[N_MAX_CSC];   //[nCsc]
  float         cscNRecHits[N_MAX_CSC];   //[nCsc]
  float         cscNRecHits_flag[N_MAX_CSC];   //[nCsc]
  float         cscNRecHits_jetveto0p4[N_MAX_CSC];   //[nCsc]
  float         cscNRecHits_jetveto0p8[N_MAX_CSC];   //[nCsc]
  float         cscT[N_MAX_CSC];   //[nCsc]
  float         cscChi2[N_MAX_CSC];   //[nCsc]

  int           nCscClusters;
  bool          cscClusterJetVeto[N_MAX_CSC];   //[nCsc]
  bool          cscClusterMuonVeto[N_MAX_CSC];   //[nCsc]

  float         cscClusterX[N_MAX_CSC];   //[nCsc]
  float         cscClusterY[N_MAX_CSC];   //[nCsc]
  float         cscClusterZ[N_MAX_CSC];   //[nCsc]
  float         cscClusterRadius[N_MAX_CSC];   //[nCsc]
  float         cscClusterMajorAxis[N_MAX_CSC];
  float         cscClusterMinorAxis[N_MAX_CSC];


  float         cscClusterEta[N_MAX_CSC];   //[nCsc]
  float         cscClusterPhi[N_MAX_CSC];   //[nCsc]
  int           cscClusterSize[N_MAX_CSC];
  int           cscClusterStation[N_MAX_CSC];


  int           nCsc_Me11Veto;
  int           nCsc_Me12Veto;
  int           nCsc_Me1112Veto;

  int           nCsc_JetVetoCluster0p4;
  int           nCsc_JetMuonVetoCluster0p4;
  int           nCsc_JetVetoCluster0p4_Me11Veto;
  int           nCsc_JetVetoCluster0p4_Me1112Veto;
  int           nCsc_JetMuonVetoCluster0p4_Me1112Veto;
  int           nCsc0_JetMuonVetoCluster0p4;
  int           nCsc_recoJetVeto0p4;
  int           nCsc_recoJetVeto0p8;
  int           nCsc_JetLepVeto0p4;
  int           nCsc_JetLepVeto0p8;
  int           nCsc_recoJetVeto0p4_Me11Veto;
  int           nCsc_recoJetVeto0p8_Me11Veto;
  int           nCsc_JetLepVeto0p4_Me11Veto;
  int           nCsc_JetLepVeto0p8_Me11Veto;
  int           nCsc_recoJetVeto0p4_Me1112Veto;
  int           nCsc_recoJetVeto0p8_Me1112Veto;
  int           nCsc_JetLepVeto0p4_Me1112Veto;
  int           nCsc_JetLepVeto0p8_Me1112Veto;
  bool          event_recoJetVeto0p4;
  bool          event_recoJetVeto0p8;
  bool          event_JetLepVeto0p4;
  bool          event_JetLepVeto0p8;
  bool          event_Me1112Veto;
  bool          event_Me11Veto;
  bool          event_Me12Veto;
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
  float ecalNRechits[N_MAX_JETS];
  float ecalRechitE[N_MAX_JETS];
  float jetChargedEMEnergyFraction[N_MAX_JETS];
  float jetNeutralEMEnergyFraction[N_MAX_JETS];
  float jetChargedHadronEnergyFraction[N_MAX_JETS];
  float jetNeutralHadronEnergyFraction[N_MAX_JETS];


  // bool jetLoosePassId[N_MAX_JETS];
  bool jetPassId[N_MAX_JETS];
  // bool jetTightPassId[N_MAX_JETS];
  bool HLTDecision[NTriggersMAX];

  UInt_t wzevtNum,trig, trig_lepId, trig_lepId_dijet; //number of events that pass each criteria



  TTree *tree_;
  TFile *f_;

  LiteTreeMuonSystem()
  {
    InitVariables();
  };

  ~LiteTreeMuonSystem()
  {
    if (f_) f_->Close();
  };

  void InitVariables()
  {
    runNum=0; lumiSec=0; evtNum=0; category=0;
    npv=0; npu=0; rho=-1; weight=-1;
    met=-1; metPhi=-1;
    lepGenId = 0;
    //CSC
    nCsc = 0;
    nCscClusters = 0;
    nCsc_JetVetoCluster0p4 = 0;
    nCsc_JetMuonVetoCluster0p4 = 0;
    nCsc_JetVetoCluster0p4_Me11Veto = 0;
    nCsc_JetVetoCluster0p4_Me1112Veto = 0;
    nCsc_JetMuonVetoCluster0p4_Me1112Veto = 0;
    nCsc0_JetMuonVetoCluster0p4 = 0;
    // nCsc_JetVetoCluster0p4_Me1112Veto = 0;
    // nCsc0_JetVetoCluster0p4_Me1112Veto = 0;
    nCsc_Me11Veto = 0;
    nCsc_Me12Veto = 0;
    nCsc_Me1112Veto = 0;
    nCsc_recoJetVeto0p4 = 0;
    nCsc_recoJetVeto0p8 = 0;
    nCsc_JetLepVeto0p4 = 0;
    nCsc_JetLepVeto0p8 = 0;
    nCsc_recoJetVeto0p4_Me11Veto = 0;
    nCsc_recoJetVeto0p8_Me11Veto = 0;
    nCsc_JetLepVeto0p4_Me11Veto = 0;
    nCsc_JetLepVeto0p8_Me11Veto = 0;
    nCsc_recoJetVeto0p4_Me1112Veto = 0;
    nCsc_recoJetVeto0p8_Me1112Veto = 0;
    nCsc_JetLepVeto0p4_Me1112Veto = 0;
    nCsc_JetLepVeto0p8_Me1112Veto = 0;

    event_recoJetVeto0p4 = true;
    event_recoJetVeto0p8 = true;
    event_JetLepVeto0p4 = true;
    event_JetLepVeto0p8 = true;
    event_Me11Veto = true;
    event_Me12Veto = true;
    event_Me1112Veto = true;
    for( int i = 0; i < N_MAX_CSC; i++ )
    {
      cscLabels[i] = -999;
      cscPhi[i] = -999;   //[nCsc]
      cscEta[i] = -999;   //[nCsc]
      cscX[i] = -999;   //[nCsc]
      cscY[i] = -999;   //[nCsc]
      cscZ[i] = -999;   //[nCsc]
      cscNRecHits[i] = -999;   //[nCsc]
      cscNRecHits_flag[i] = -999;   //[nCsc]
      cscT[i] = -999;   //[nCsc]
      cscChi2[i] = -999;   //[nCsc]
      cscClusterSize[i] = -999;
      cscClusterX[i] = -999.;
      cscClusterY[i] = -999.;
      cscClusterZ[i] = -999.;
      cscClusterRadius[i] = -999.;
      cscClusterMajorAxis[i] = -999.;
      cscClusterMinorAxis[i] = -999.;

      cscClusterEta[i] = -999.;
      cscClusterPhi[i] = -999.;
      cscClusterJetVeto[i] = true;
      cscClusterMuonVeto[i] = true;
      cscClusterStation[i] = -999;
    }

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
      ecalNRechits[i] = -999.;
      ecalRechitE[i] = -999.;
      jetChargedEMEnergyFraction[i] = -999.;
      jetNeutralEMEnergyFraction[i] = -999.;
      jetChargedHadronEnergyFraction[i] = -999.;
      jetNeutralHadronEnergyFraction[i] = -999.;
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
    tree_ = dynamic_cast<TTree*>(f_->Get("MuonSystem"));
    InitTree();
    assert(tree_);
  };

  void CreateTree()
  {
    tree_ = new TTree("MuonSystem","MuonSystem");
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
    tree_->Branch("lepGenId",      &lepGenId,     "lepGenId/I");      // phi(MET)

    //CSC
    tree_->Branch("nCsc",             &nCsc, "nCsc/I");
    tree_->Branch("nCscClusters",             &nCscClusters, "nCscClusters/I");
    tree_->Branch("cscClusterX",             cscClusterX,             "cscClusterX[nCscClusters]/F");
    tree_->Branch("cscClusterY",             cscClusterY,             "cscClusterY[nCscClusters]/F");
    tree_->Branch("cscClusterZ",             cscClusterZ,             "cscClusterZ[nCscClusters]/F");
    tree_->Branch("cscClusterRadius",             cscClusterRadius,             "cscClusterRadius[nCscClusters]/F");
    tree_->Branch("cscClusterMajorAxis",             cscClusterMajorAxis,             "cscClusterMajorAxis[nCscClusters]/F");
    tree_->Branch("cscClusterMinorAxis",             cscClusterMinorAxis,             "cscClusterMinorAxis[nCscClusters]/F");

    tree_->Branch("cscClusterPhi",             cscClusterPhi,             "cscClusterPhi[nCscClusters]/F");
    tree_->Branch("cscClusterEta",             cscClusterEta,             "cscClusterEta[nCscClusters]/F");
    tree_->Branch("cscClusterJetVeto",             cscClusterJetVeto,             "cscClusterJetVeto[nCscClusters]/O");
    tree_->Branch("cscClusterMuonVeto",             cscClusterMuonVeto,             "cscClusterMuonVeto[nCscClusters]/O");

    tree_->Branch("cscClusterSize",             cscClusterSize,             "cscClusterSize[nCscClusters]/I");
    tree_->Branch("cscClusterStation",             cscClusterStation,             "cscClusterStation[nCscClusters]/I");

    tree_->Branch("nCsc_Me11Veto",             &nCsc_Me11Veto,"nCsc_Me11Veto/I");
    tree_->Branch("nCsc_Me12Veto",             &nCsc_Me12Veto,"nCsc_Me12Veto/I");
    tree_->Branch("nCsc_Me1112Veto",             &nCsc_Me1112Veto,"nCsc_Me1112Veto/I");
    tree_->Branch("nCsc_JetVetoCluster0p4",             &nCsc_JetVetoCluster0p4,"nCsc_JetVetoCluster0p4/I");
    tree_->Branch("nCsc_JetMuonVetoCluster0p4",             &nCsc_JetMuonVetoCluster0p4,"nCsc_JetMuonVetoCluster0p4/I");

    tree_->Branch("nCsc_JetVetoCluster0p4_Me11Veto",             &nCsc_JetVetoCluster0p4_Me11Veto,"nCsc_JetVetoCluster0p4_Me11Veto/I");
    tree_->Branch("nCsc_JetVetoCluster0p4_Me1112Veto",             &nCsc_JetVetoCluster0p4_Me1112Veto,"nCsc_JetVetoCluster0p4_Me1112Veto/I");
    tree_->Branch("nCsc_JetMuonVetoCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoCluster0p4_Me1112Veto,"nCsc_JetMuonVetoCluster0p4_Me1112Veto/I");

    tree_->Branch("nCsc0_JetMuonVetoCluster0p4",             &nCsc0_JetMuonVetoCluster0p4,"nCsc0_JetMuonVetoCluster0p4/I");

    tree_->Branch("nCsc_recoJetVeto0p4",             &nCsc_recoJetVeto0p4,"nCsc_recoJetVeto0p4/I");
    tree_->Branch("nCsc_recoJetVeto0p8",             &nCsc_recoJetVeto0p8, "nCsc_recoJetVeto0p8/I");
    tree_->Branch("nCsc_JetLepVeto0p4",             &nCsc_JetLepVeto0p4, "nCsc_JetLepVeto0p4/I");
    tree_->Branch("nCsc_JetLepVeto0p8",             &nCsc_JetLepVeto0p8, "nCsc_JetLepVeto0p8/I");

    tree_->Branch("nCsc_recoJetVeto0p4_Me11Veto",             &nCsc_recoJetVeto0p4_Me11Veto, "nCsc_recoJetVeto0p4_Me11Veto/I");
    tree_->Branch("nCsc_recoJetVeto0p8_Me11Veto",             &nCsc_recoJetVeto0p8_Me11Veto, "nCsc_recoJetVeto0p8_Me11Veto/I");
    tree_->Branch("nCsc_JetLepVeto0p4_Me11Veto",             &nCsc_JetLepVeto0p4_Me11Veto, "nCsc_JetLepVeto0p4_Me11Veto/I");
    tree_->Branch("nCsc_JetLepVeto0p8_Me11Veto",             &nCsc_JetLepVeto0p8_Me11Veto, "nCsc_JetLepVeto0p8_Me11Veto/I");

    tree_->Branch("nCsc_recoJetVeto0p4_Me1112Veto",             &nCsc_recoJetVeto0p4_Me1112Veto, "nCsc_recoJetVeto0p4_Me1112Veto/I");
    tree_->Branch("nCsc_recoJetVeto0p8_Me1112Veto",             &nCsc_recoJetVeto0p8_Me1112Veto, "nCsc_recoJetVeto0p8_Me1112Veto/I");
    tree_->Branch("nCsc_JetLepVeto0p4_Me1112Veto",             &nCsc_JetLepVeto0p4_Me1112Veto, "nCsc_JetLepVeto0p4_Me1112Veto/I");
    tree_->Branch("nCsc_JetLepVeto0p8_Me1112Veto",             &nCsc_JetLepVeto0p8_Me1112Veto, "nCsc_JetLepVeto0p8_Me1112Veto/I");


    tree_->Branch("event_recoJetVeto0p4",             &event_recoJetVeto0p4, "event_recoJetVeto0p4/O");
    tree_->Branch("event_recoJetVeto0p8",             &event_recoJetVeto0p8, "event_recoJetVeto0p8/O");
    tree_->Branch("event_JetLepVeto0p4",             &event_JetLepVeto0p4, "event_JetLepVeto0p4/O");
    tree_->Branch("event_JetLepVeto0p8",             &event_JetLepVeto0p8, "event_JetLepVeto0p8/O");

    tree_->Branch("event_Me11Veto",             &event_Me11Veto,"event_Me11Veto/O");
    tree_->Branch("event_Me12Veto",             &event_Me12Veto,"event_Me12Veto/O");
    tree_->Branch("event_Me1112Veto",             &event_Me1112Veto,"event_Me1112Veto/O");

    tree_->Branch("cscLabels",             cscLabels,             "cscLabels[nCsc]/I");

    tree_->Branch("cscPhi",           cscPhi,           "cscPhi[nCsc]/F");
    tree_->Branch("cscEta",           cscEta,           "cscEta[nCsc]/F");
    tree_->Branch("cscX",             cscX,             "cscX[nCsc]/F");
    tree_->Branch("cscY",             cscY,             "cscY[nCsc]/F");
    tree_->Branch("cscZ",             cscZ,             "cscZ[nCsc]/F");
    tree_->Branch("cscNRecHits",      cscNRecHits,      "cscNRecHits[nCsc]/F");
    tree_->Branch("cscNRecHits_flag", cscNRecHits_flag, "cscNRecHits_flag[nCsc]/F");
    tree_->Branch("cscT",             cscT,             "cscT[nCsc]/F");
    tree_->Branch("cscChi2",          cscChi2,          "cscChi2[nCsc]/F");

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
    /*
    //Z-candidate
    tree_->Branch("MT",      &MT,  "MT/F");
    tree_->Branch("ZMass",      &ZMass,  "ZMass/F");
    tree_->Branch("ZPt",        &ZPt,    "ZPt/F");
    tree_->Branch("ZEta",       &ZEta,   "ZEta/F");
    tree_->Branch("ZPhi",       &ZPhi,   "ZPhi/F");
    tree_->Branch("ZleptonIndex1", &ZleptonIndex1, "ZleptonIndex1/I");
    tree_->Branch("ZleptonIndex2", &ZleptonIndex2, "ZleptonIndex2/I");
    */
    //jets
    tree_->Branch("nJets",     &nJets,    "nJets/I");
    tree_->Branch("jetE",      jetE,      "jetE[nJets]/F");
    tree_->Branch("jetPt",     jetPt,     "jetPt[nJets]/F");
    tree_->Branch("jetEta",    jetEta,    "jetEta[nJets]/F");
    tree_->Branch("jetPhi",    jetPhi,    "jetPhi[nJets]/F");
    // tree_->Branch("jetTime",   jetTime,   "jetTime[nJets]/F");
    tree_->Branch("jetPassId", jetPassId, "jetPassId[nJets]/O");
    // tree_->Branch("ecalNRechits",   ecalNRechits,   "ecalNRechits[nJets]/F");
    // tree_->Branch("ecalRechitE", ecalRechitE, "ecalRechitE[nJets]/F");
    // tree_->Branch("jetLoosePassId", jetLoosePassId, "jetLoosePassId[nJets]/O");
    // tree_->Branch("jetTightPassId", jetTightPassId, "jetTightPassId[nJets]/O");
    tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[601]/O"); //hardcoded
    // tree_->Branch("jetChargedEMEnergyFraction",   jetChargedEMEnergyFraction,   "jetChargedEMEnergyFraction[nJets]/F");
    // tree_->Branch("jetNeutralEMEnergyFraction",   jetNeutralEMEnergyFraction,   "jetNeutralEMEnergyFraction[nJets]/F");
    // tree_->Branch("jetChargedHadronEnergyFraction",   jetChargedHadronEnergyFraction,   "jetChargedHadronEnergyFraction[nJets]/F");
    // tree_->Branch("jetNeutralHadronEnergyFraction",   jetNeutralHadronEnergyFraction,   "jetNeutralHadronEnergyFraction[nJets]/F");
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
    tree_->SetBranchAddress("lepGenId",      &lepGenId);
    //CSC
    tree_->SetBranchAddress("nCsc",             &nCsc);
    tree_->SetBranchAddress("cscClusterStation",             &cscClusterStation);

    tree_->SetBranchAddress("nCscClusters",             &nCscClusters);
    tree_->SetBranchAddress("cscClusterX",             cscClusterX);
    tree_->SetBranchAddress("cscClusterY",             cscClusterY);
    tree_->SetBranchAddress("cscClusterZ",             cscClusterZ);
    tree_->SetBranchAddress("cscClusterRadius",             cscClusterRadius);
    tree_->SetBranchAddress("cscClusterMajorAxis",             cscClusterMajorAxis);
    tree_->SetBranchAddress("cscClusterMinorAxis",             cscClusterMinorAxis);

    tree_->SetBranchAddress("cscClusterEta",             cscClusterEta);
    tree_->SetBranchAddress("cscClusterPhi",             cscClusterPhi);
    tree_->SetBranchAddress("cscClusterJetVeto",             cscClusterJetVeto);
    tree_->SetBranchAddress("cscClusterMuonVeto",             cscClusterMuonVeto);

    tree_->SetBranchAddress("cscClusterSize",             cscClusterSize);

    tree_->SetBranchAddress("nCsc_Me11Veto",             &nCsc_Me11Veto);
    tree_->SetBranchAddress("nCsc_Me12Veto",             &nCsc_Me12Veto);
    tree_->SetBranchAddress("nCsc_Me1112Veto",             &nCsc_Me1112Veto);

    tree_->SetBranchAddress("nCsc_recoJetVeto0p4",             &nCsc_recoJetVeto0p4);
    tree_->SetBranchAddress("nCsc_recoJetVeto0p8",             &nCsc_recoJetVeto0p8);
    tree_->SetBranchAddress("nCsc_JetLepVeto0p4",             &nCsc_JetLepVeto0p4);
    tree_->SetBranchAddress("nCsc_JetLepVeto0p8",             &nCsc_JetLepVeto0p8);

    tree_->SetBranchAddress("nCsc_recoJetVeto0p4_Me11Veto",             &nCsc_recoJetVeto0p4_Me11Veto);
    tree_->SetBranchAddress("nCsc_recoJetVeto0p8_Me11Veto",             &nCsc_recoJetVeto0p8_Me11Veto);
    tree_->SetBranchAddress("nCsc_JetLepVeto0p4_Me11Veto",             &nCsc_JetLepVeto0p4_Me11Veto);
    tree_->SetBranchAddress("nCsc_JetLepVeto0p8_Me11Veto",             &nCsc_JetLepVeto0p8_Me11Veto);

    tree_->SetBranchAddress("nCsc_recoJetVeto0p4_Me1112Veto",             &nCsc_recoJetVeto0p4_Me1112Veto);
    tree_->SetBranchAddress("nCsc_recoJetVeto0p8_Me1112Veto",             &nCsc_recoJetVeto0p8_Me1112Veto);
    tree_->SetBranchAddress("nCsc_JetLepVeto0p4_Me1112Veto",             &nCsc_JetLepVeto0p4_Me1112Veto);
    tree_->SetBranchAddress("nCsc_JetLepVeto0p8_Me1112Veto",             &nCsc_JetLepVeto0p8_Me1112Veto);


    tree_->SetBranchAddress("event_recoJetVeto0p4",             &event_recoJetVeto0p4);
    tree_->SetBranchAddress("event_recoJetVeto0p8",             &event_recoJetVeto0p8);
    tree_->SetBranchAddress("event_JetLepVeto0p4",             &event_JetLepVeto0p4);
    tree_->SetBranchAddress("event_JetLepVeto0p8",             &event_JetLepVeto0p8);

    tree_->SetBranchAddress("event_Me11Veto",             &event_Me11Veto);
    tree_->SetBranchAddress("event_Me12Veto",             &event_Me12Veto);
    tree_->SetBranchAddress("event_Me1112Veto",             &event_Me1112Veto);


    tree_->SetBranchAddress("cscLabels",             cscLabels);

    tree_->SetBranchAddress("cscPhi",           cscPhi);
    tree_->SetBranchAddress("cscEta",           cscEta);
    tree_->SetBranchAddress("cscX",             cscX);
    tree_->SetBranchAddress("cscY",             cscY);
    tree_->SetBranchAddress("cscZ",             cscZ);
    tree_->SetBranchAddress("cscNRecHits",      cscNRecHits);
    tree_->SetBranchAddress("cscNRecHits_flag", cscNRecHits_flag);
    tree_->SetBranchAddress("cscT",             cscT);
    tree_->SetBranchAddress("cscChi2",          cscChi2);
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
    /*
    //Z-candidate
    tree_->SetBranchAddress("ZMass",       &ZMass);
    tree_->SetBranchAddress("ZPt",         &ZPt);
    tree_->SetBranchAddress("ZEta",        &ZEta);
    tree_->SetBranchAddress("ZPhi",        &ZPhi);
    tree_->SetBranchAddress("ZleptonIndex1", &ZleptonIndex1);
    tree_->SetBranchAddress("ZleptonIndex2", &ZleptonIndex2);
    tree_->SetBranchAddress("MT", &MT);
    */
    //jets
    tree_->SetBranchAddress("nJets",     &nJets);
    tree_->SetBranchAddress("jetE",      jetE);
    tree_->SetBranchAddress("jetPt",     jetPt);
    tree_->SetBranchAddress("jetEta",    jetEta);
    tree_->SetBranchAddress("jetPhi",    jetPhi);
    // tree_->SetBranchAddress("jetTime",   jetTime);
    tree_->SetBranchAddress("jetPassId", jetPassId);
    // tree_->SetBranchAddress("ecalNRechits",   ecalNRechits);
    // tree_->SetBranchAddress("ecalRechitE", ecalRechitE);
    // tree_->SetBranchAddress("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction);
    // tree_->SetBranchAddress("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction);
    // tree_->SetBranchAddress("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction);
    // tree_->SetBranchAddress("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction);

    // tree_->SetBranchAddress("jetLoosePassId", jetLoosePassId);
    // tree_->SetBranchAddress("jetTightPassId", jetTightPassId);
    // triggers
    tree_->SetBranchAddress("HLTDecision",   HLTDecision);
  };

};

void llp_MuonSystem::Analyze(bool isData, int options, string outputfilename, string analysisTag)
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
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
  LiteTreeMuonSystem *MuonSystem = new LiteTreeMuonSystem;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();
  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *generatedEvents = new TH1F("generatedEvents", "generatedEvents", 1, 1, 2);
  TH1F *trig = new TH1F("trig", "trig", 1, 1, 2);
  TH1F *trig_lepId = new TH1F("trig_lepId", "trig_lepId", 1, 1, 2);
  TH1F *trig_lepId_dijet = new TH1F("trig_lepId_dijet", "trig_lepId_dijet", 1, 1, 2);


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

  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {

    //begin event
    if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //GetEntry(ientry);
    nb = fChain->GetEntry(jentry); nbytes += nb;

    //fill normalization histogram
    //std::cout << "deb0 " << jentry << std::endl;
    MuonSystem->InitVariables();
    //std::cout << "deb1 " << jentry << std::endl;
    if (label =="bkg_wH"|| label == "bkg_zH"){
      if (isData)
      {
        NEvents->Fill(1);
        MuonSystem->weight = 1;
      }
      else
      {
        //NEvents->Fill(genWeight);
        //MuonSystem->weight = genWeight;
        NEvents->Fill(1);
        MuonSystem->weight = 1;
      }

    }
    else{
      NEvents->Fill(1);

      MuonSystem->weight = 1;
    }
    //std::cout << "deb2 " << jentry << std::endl;
    //event info
    MuonSystem->runNum = runNum;
    MuonSystem->lumiSec = lumiNum;
    MuonSystem->evtNum = eventNum;
    //std::cout << "deb3 " << jentry << std::endl;

    bool wzFlag = false;
    for (int i=0; i < nGenParticle; ++i)
    {

      if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == wzId)
      { // choosing only the W->munu events
        wzFlag = true;
        MuonSystem->lepGenId = gParticleId[i];
      }

    }
    if ( wzFlag == true ) generatedEvents->Fill(1);;
    // NEvents->Fill(1);

    for (int i=0; i < nBunchXing; ++i)
    {
      if (BunchXing[i] == 0)
      {
        MuonSystem->npu = nPUmean[i];
      }
    }
    //get NPU
    MuonSystem->npv = nPV;
    MuonSystem->rho = fixedGridRhoFastjetAll;
    MuonSystem->met = metType1Pt;
    MuonSystem->metPhi = metType1Phi;

    //Triggers
    for(int i = 0; i < NTriggersMAX; i++){
      MuonSystem->HLTDecision[i] = HLTDecision[i];
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
    //CSC INFO
    if( nCsc < 30 ) continue;//require at least 30 segments in the CSCs
    MuonSystem->nCsc = nCsc;
    vector<Point> points;

    for(int i = 0; i < nCsc; i++)
    {
      bool me11 = false;
      bool me12 = false;
      MuonSystem->cscPhi[i]           = cscPhi[i];   //[nCsc]
      MuonSystem->cscEta[i]           = cscEta[i];   //[nCsc]
      MuonSystem->cscX[i]             = cscX[i];   //[nCsc]
      MuonSystem->cscY[i]             = cscY[i];   //[nCsc]
      MuonSystem->cscZ[i]             = cscZ[i];   //[nCsc]
      // for dbscan
      Point p;
      p.x = cscX[i];
      p.y = cscY[i];
      p.z = cscZ[i];
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);

      MuonSystem->cscNRecHits[i]      = cscNRecHits[i];   //[nCsc]
      MuonSystem->cscNRecHits_flag[i] = cscNRecHits_flag[i];   //[nCsc]
      MuonSystem->cscT[i]             = cscT[i];   //[nCsc]
      MuonSystem->cscChi2[i]          = cscChi2[i];   //[nCsc]
      float cscR = sqrt(cscX[i]*cscX[i]+cscY[i]*cscY[i]);

      //page 141 of tdr: https://cds.cern.ch/record/343814/files/LHCC-97-032.pdf
      if (abs(cscZ[i]) > 568 && abs(cscZ[i]) < 632 ) me11 = true;
      if (abs(cscZ[i]) > 663 && abs(cscZ[i]) < 724 && abs(cscR) > 275 && abs(cscR) < 465) me12 = true;
      if (me11) MuonSystem->nCsc_Me11Veto++;
      if (me12) MuonSystem->nCsc_Me12Veto++;
    }
    if (MuonSystem->nCsc_Me11Veto > N_CSC_CUT)MuonSystem->event_Me11Veto = false;
    if (MuonSystem->nCsc_Me12Veto > N_CSC_CUT)MuonSystem->event_Me12Veto = false;
    if ((MuonSystem->nCsc_Me11Veto + MuonSystem->nCsc_Me12Veto)> N_CSC_CUT)MuonSystem->event_Me1112Veto = false;
    MuonSystem->nCsc_Me1112Veto = MuonSystem->nCsc - (MuonSystem->nCsc_Me11Veto + MuonSystem->nCsc_Me12Veto);
    MuonSystem->nCsc_Me11Veto = MuonSystem->nCsc - MuonSystem->nCsc_Me11Veto;
    MuonSystem->nCsc_Me12Veto = MuonSystem->nCsc - MuonSystem->nCsc_Me12Veto;

    //dbscan of csc segments
    // constructor
    //v2:10,200
    //v3: 10, 100
    //v4:
    int min_point = 10;//6
    float epsilon = 100;//100
    DBSCAN ds(min_point, epsilon, points);
    int nClusters = ds.run();
    int clusterSize[nClusters] = {0};
    int cscLabels_temp[nCsc] = {-999};
    float clusterEta[nClusters] = {-999.};
    float clusterPhi[nClusters] = {-999.};
    float clusterX[nClusters] = {-999.};
    float clusterY[nClusters] = {-999.};
    float clusterZ[nClusters] = {-999.};
    float clusterRadius[nClusters] = {-999.};
    float clusterMajorAxis[nClusters] = {-999.};
    float clusterMinorAxis[nClusters] = {-999.};

    ds.result(nClusters, cscLabels_temp,clusterSize, clusterX, clusterY, clusterZ, clusterEta, clusterPhi, clusterRadius);
    ds.clusterMoments(nClusters, clusterMajorAxis, clusterMinorAxis, clusterEta, clusterPhi, clusterSize);
    for(int i = 0;i<nCsc;i++)
    {
      MuonSystem->cscLabels[i] = cscLabels_temp[i];
    }

    std::vector<cscCluster> CscCluster;
    for(int i = 0; i < nClusters; i++)
    {
      cscCluster tmpCluster;
      tmpCluster.x = clusterX[i];
      tmpCluster.y = clusterY[i];
      tmpCluster.z = clusterZ[i];
      tmpCluster.radius = clusterRadius[i];
      tmpCluster.MajorAxis = clusterMajorAxis[i];
      tmpCluster.MinorAxis = clusterMinorAxis[i];

      tmpCluster.eta = clusterEta[i];
      tmpCluster.phi = clusterPhi[i];
      tmpCluster.nCscSegments = clusterSize[i];
      tmpCluster.station = cscStation(clusterX[i],clusterY[i],clusterZ[i]);

      bool nCscFlag_recoJetVeto0p4 = true;
      bool nCscFlag_recoMuonVeto0p4 = true;

      bool me11 = false;
      bool me12 = false;
      float clusterR = sqrt(clusterX[i]*clusterX[i]+clusterY[i]*clusterY[i]);
      if (abs(clusterZ[i]) > 568 && abs(clusterZ[i]) < 632 ) me11 = true;
      if (abs(clusterZ[i]) > 663 && abs(clusterZ[i]) < 724 && abs(clusterR) > 275 && abs(clusterR) < 465) me12 = true;
      for (int j = 0; j < nJets; j++)
      {
        if (jetPt[j]<JET_PT_CUT) continue;
        if (abs(jetEta[j])>3) continue;

        if (RazorAnalyzer::deltaR(clusterEta[i],clusterPhi[i],jetEta[j],jetPhi[j]) < 0.4) nCscFlag_recoJetVeto0p4 = false;
      }
      for (int j = 0; j < nMuons; j++)
      {
        if (muonPt[j]<20) continue;
        if (abs(muonEta[j])>3) continue;

        if (RazorAnalyzer::deltaR(clusterEta[i],clusterPhi[i],muonEta[j],muonPhi[j]) < 0.4) nCscFlag_recoMuonVeto0p4 = false;
      }
      tmpCluster.jetVeto = nCscFlag_recoJetVeto0p4;
      tmpCluster.muonVeto = nCscFlag_recoMuonVeto0p4;

      if (nCscFlag_recoJetVeto0p4) MuonSystem->nCsc_JetVetoCluster0p4 += clusterSize[i];
      if (nCscFlag_recoJetVeto0p4 && nCscFlag_recoMuonVeto0p4) MuonSystem->nCsc_JetMuonVetoCluster0p4+= clusterSize[i];

      if (nCscFlag_recoJetVeto0p4 && (!me11)) MuonSystem->nCsc_JetVetoCluster0p4_Me11Veto+= clusterSize[i];
      if (nCscFlag_recoJetVeto0p4 && (!me11) && (!me12)) MuonSystem->nCsc_JetVetoCluster0p4_Me1112Veto+= clusterSize[i];
      if (nCscFlag_recoJetVeto0p4 && nCscFlag_recoMuonVeto0p4 && (!me11) && (!me12)) MuonSystem->nCsc_JetMuonVetoCluster0p4_Me1112Veto+= clusterSize[i];

      CscCluster.push_back(tmpCluster);
    }
    sort(CscCluster.begin(), CscCluster.end(), my_largest_nCsc);
    //most populated cluster with jet veto
    bool first = true;
    for ( auto &tmp : CscCluster )
    {
      MuonSystem->cscClusterX[MuonSystem->nCscClusters] =tmp.x;
      MuonSystem->cscClusterY[MuonSystem->nCscClusters] =tmp.y;
      MuonSystem->cscClusterZ[MuonSystem->nCscClusters] =tmp.z;
      MuonSystem->cscClusterEta[MuonSystem->nCscClusters] =tmp.eta;
      MuonSystem->cscClusterPhi[MuonSystem->nCscClusters] = tmp.phi;
      MuonSystem->cscClusterRadius[MuonSystem->nCscClusters] =tmp.radius;
      MuonSystem->cscClusterMajorAxis[MuonSystem->nCscClusters] =tmp.MajorAxis;
      MuonSystem->cscClusterMinorAxis[MuonSystem->nCscClusters] =tmp.MinorAxis;

      MuonSystem->cscClusterSize[MuonSystem->nCscClusters] = tmp.nCscSegments;
      MuonSystem->cscClusterJetVeto[MuonSystem->nCscClusters] = tmp.jetVeto;
      MuonSystem->cscClusterMuonVeto[MuonSystem->nCscClusters] = tmp.muonVeto;
      MuonSystem->cscClusterStation[MuonSystem->nCscClusters] = tmp.station;
      if (tmp.jetVeto && tmp.muonVeto && first){
        first = false;
        MuonSystem->nCsc0_JetMuonVetoCluster0p4 = tmp.nCscSegments;
      }

      // std::cout << "lepton pdg " << MuonSystem->lepPdgId[MuonSystem->nLeptons] << std::endl;
      MuonSystem->nCscClusters++;
    }





    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
    for( int i = 0; i < nMuons; i++ )
    {
      // if(!isMuonPOGLooseMuon(i)) continue;
      if(muonPt[i] < 20) continue;
      if(fabs(muonEta[i]) > 3) continue;

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
      // tmpElectron.passId = passMVALooseElectronID(i) && passEGammaPOGLooseElectronIso(i);
      Leptons.push_back(tmpElectron);
    }

    sort(Leptons.begin(), Leptons.end(), my_largest_pt);
    //std::cout << "deb7 " << jentry << std::endl;
    for ( auto &tmp : Leptons )
    {
      MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
      MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
      MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
      MuonSystem->lepPassId[MuonSystem->nLeptons] = tmp.passId;


      // std::cout << "lepton pdg " << MuonSystem->lepPdgId[MuonSystem->nLeptons] << std::endl;
      MuonSystem->nLeptons++;
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
      MuonSystem->ZMass = ZMass;
      MuonSystem->ZPt   = ZPt;
      MuonSystem->ZEta  = ZCandidate.Eta();
      MuonSystem->ZPhi  = ZCandidate.Phi();
      MuonSystem->ZleptonIndex1 = ZCandidateLeptonIndex.first;
      MuonSystem->ZleptonIndex2 = ZCandidateLeptonIndex.second;

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
    //   std::cout <<nMuons << "," << nElectrons <<  "," << MuonSystem->nLeptons <<  "," << MuonSystem->met << std::endl;
    // }

    //if ( Leptons.size() < nLepton_cut ) continue;
    TLorentzVector met;

    met.SetPtEtaPhiE(metType1Pt,0,metType1Phi,metType1Pt);
    if ( Leptons.size() > 0 )
    {
      TLorentzVector visible = Leptons[0].lepton;
      MuonSystem->MT = GetMT(visible,met);
    }
    // else{
    //   if ( Leptons.size() < 2 ) continue;
    //   if (!(foundZ && fabs(ZMass-Z_MASS) < 15.0 )) continue;
    // }
    if (triggered) trig_lepId->Fill(1);
  //Jet vetoing
  for (int i = 0; i < nCsc; i++)
  {
    bool nCscFlag_recoJetVeto0p4 = true;
    bool nCscFlag_recoJetVeto0p8 = true;
    bool me11 = false;
    bool me12 = false;
    float cscR = sqrt(cscX[i]*cscX[i]+cscY[i]*cscY[i]);
    if (abs(cscZ[i]) > 568 && abs(cscZ[i]) < 632 ) me11 = true;
    if (abs(cscZ[i]) > 663 && abs(cscZ[i]) < 724 && abs(cscR) > 275 && abs(cscR) < 465) me12 = true;
    for (int j = 0; j < nJets; j++)
    {
      if (jetPt[j]<JET_PT_CUT) continue;
      if (abs(jetEta[j])>3) continue;
      if (RazorAnalyzer::deltaR(cscEta[i],cscPhi[i],jetEta[j],jetPhi[j]) < 0.4) {

        nCscFlag_recoJetVeto0p4 = false;
      }
      if (RazorAnalyzer::deltaR(cscEta[i],cscPhi[i],jetEta[j],jetPhi[j]) < 0.8){
        nCscFlag_recoJetVeto0p8 = false;
      }

    }
    if (nCscFlag_recoJetVeto0p4) MuonSystem->nCsc_recoJetVeto0p4++;
    if (nCscFlag_recoJetVeto0p8) MuonSystem->nCsc_recoJetVeto0p8++;
    if (nCscFlag_recoJetVeto0p4 && (!me11)) MuonSystem->nCsc_recoJetVeto0p4_Me11Veto++;
    if (nCscFlag_recoJetVeto0p8 && (!me11)) MuonSystem->nCsc_recoJetVeto0p8_Me11Veto++;
    if (nCscFlag_recoJetVeto0p4 && (!me11) && (!me12)) MuonSystem->nCsc_recoJetVeto0p4_Me1112Veto++;
    if (nCscFlag_recoJetVeto0p8 && (!me11) && (!me12)) MuonSystem->nCsc_recoJetVeto0p8_Me1112Veto++;

  }
  if ((MuonSystem->nCsc - MuonSystem->nCsc_recoJetVeto0p4)  > N_CSC_CUT)  MuonSystem->event_recoJetVeto0p4 = false;
  if ((MuonSystem->nCsc -MuonSystem->nCsc_recoJetVeto0p8)  > N_CSC_CUT)  MuonSystem->event_recoJetVeto0p8 = false;

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

      if( thisJet.Pt() < 5 ) continue;//According to the April 1st 2015 AN
      if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
      // if ( !jetPassIDLoose[i] ) continue;
      // if (!(jetRechitE[i] > 0.0)) continue;

      // std::cout <<jetRechitT[i] << "," << jetRechitE[i] <<  "," << jetNRechits[i] << std::endl;


      jets tmpJet;
      tmpJet.jet    = thisJet;
      tmpJet.time   = jetRechitT[i];
      tmpJet.passId = jetPassIDTight[i];
      tmpJet.isCSVL = isCSVL(i);
      //if (isCSVL(i)) NBJet20++;
      //if (isCSVL(i) && thisJet.Pt() > 30) NBJet30++;
      tmpJet.ecalNRechits = jetNRechits[i];
      tmpJet.ecalRechitE = jetRechitE[i];
      tmpJet.jetChargedEMEnergyFraction = jetChargedEMEnergyFraction[i];
      tmpJet.jetNeutralEMEnergyFraction = jetNeutralEMEnergyFraction[i];
      tmpJet.jetChargedHadronEnergyFraction = jetChargedHadronEnergyFraction[i];
      tmpJet.jetNeutralHadronEnergyFraction = jetNeutralHadronEnergyFraction[i];
      Jets.push_back(tmpJet);

    }

    //-----------------------------
    //Require at least 2 jets
    //-----------------------------
    //if( Jets.size() < 2 ) continue;
    if (triggered) trig_lepId_dijet->Fill(1);
    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

    for ( auto &tmp : Jets )
    {
      MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
      MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
      MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
      MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
      MuonSystem->jetTime[MuonSystem->nJets] = tmp.time;
      MuonSystem->jetPassId[MuonSystem->nJets] = tmp.passId;
      MuonSystem->ecalNRechits[MuonSystem->nJets] = tmp.ecalNRechits;
      MuonSystem->ecalNRechits[MuonSystem->nJets] = tmp.ecalRechitE;
      MuonSystem->jetChargedEMEnergyFraction[MuonSystem->nJets] = tmp.jetChargedEMEnergyFraction;
      MuonSystem->jetNeutralEMEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralEMEnergyFraction;
      MuonSystem->jetChargedHadronEnergyFraction[MuonSystem->nJets] = tmp.jetChargedHadronEnergyFraction;
      MuonSystem->jetNeutralHadronEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralHadronEnergyFraction;
      // std::cout <<tmp.time << "," <<tmp.ecalRechitE <<  "," << tmp.ecalNRechits << MuonSystem->nJets<<std::endl;

      MuonSystem->nJets++;
    }

    //Jet vetoing
    for (int i = 0; i < nCsc; i++)
    {
      bool nCscFlag_JetLepVeto0p4 = true;
      bool nCscFlag_JetLepVeto0p8 = true;
      bool me11 = false;
      bool me12 = false;
      float cscR = sqrt(cscX[i]*cscX[i]+cscY[i]*cscY[i]);
      if (abs(cscZ[i]) > 568 && abs(cscZ[i]) < 632 ) me11 = true;
      if (abs(cscZ[i]) > 663 && abs(cscZ[i]) < 724 && abs(cscR) > 275 && abs(cscR) < 465) me12 = true;
      for (int j = 0; j < MuonSystem->nJets; j++)
      {
        if (jetPt[j]<5) continue;
        if (abs(jetEta[j])>3) continue;
        if (RazorAnalyzer::deltaR(cscEta[i],cscPhi[i],MuonSystem->jetEta[j],MuonSystem->jetPhi[j]) < 0.4) {
          nCscFlag_JetLepVeto0p4 = false;
        }
        if (RazorAnalyzer::deltaR(cscEta[i],cscPhi[i],MuonSystem->jetEta[j],MuonSystem->jetPhi[j]) < 0.8) {
          nCscFlag_JetLepVeto0p8 = false;
        }

      }
      for (int j = 0; j < MuonSystem->nLeptons; j++)
      {
        if (RazorAnalyzer::deltaR(cscEta[i],cscPhi[i],MuonSystem->lepEta[j],MuonSystem->lepPhi[j]) < 0.4) {
          nCscFlag_JetLepVeto0p4 = false;
        }
        if (RazorAnalyzer::deltaR(cscEta[i],cscPhi[i],MuonSystem->lepEta[j],MuonSystem->lepPhi[j]) < 0.8) {
          nCscFlag_JetLepVeto0p8 = false;
        }
      }

      if (nCscFlag_JetLepVeto0p4) MuonSystem->nCsc_JetLepVeto0p4++;
      if (nCscFlag_JetLepVeto0p8) MuonSystem->nCsc_JetLepVeto0p8++;
      if (nCscFlag_JetLepVeto0p4 && (!me11)) MuonSystem->nCsc_JetLepVeto0p4_Me11Veto++;
      if (nCscFlag_JetLepVeto0p8 && (!me11)) MuonSystem->nCsc_JetLepVeto0p8_Me11Veto++;
      if (nCscFlag_JetLepVeto0p4 && (!me11) && (!me12)) MuonSystem->nCsc_JetLepVeto0p4_Me1112Veto++;
      if (nCscFlag_JetLepVeto0p8 && (!me11) && (!me12)) MuonSystem->nCsc_JetLepVeto0p8_Me1112Veto++;


    }
    if ((MuonSystem->nCsc-MuonSystem->nCsc_JetLepVeto0p4)  > N_CSC_CUT)  MuonSystem->event_JetLepVeto0p4 = false;
    if ((MuonSystem->nCsc-MuonSystem->nCsc_JetLepVeto0p8)  > N_CSC_CUT)  MuonSystem->event_JetLepVeto0p8 = false;

    //std::cout << "deb fill: " << MuonSystem->nLeptons << " " << jentry << endl;
    MuonSystem->tree_->Fill();
  }

    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();
}
