#include "RazorHelper.h"
#include "LiteLiteTreeMuonSystem.h"
#include "assert.h"
#include "TTree.h"
#include "DBSCAN.h"

// Constructor
LiteLiteTreeMuonSystem::LiteLiteTreeMuonSystem()
{
  InitVariables();
};
LiteLiteTreeMuonSystem::~LiteLiteTreeMuonSystem()
{
  if (f_) f_->Close();
};
void LiteLiteTreeMuonSystem::InitVariables()
{
  runNum=0; lumiSec=0; evtNum=0; category=0;
  npv=0; npu=0;
  pileupWeight = 0;

  weight=-1.0;rho=-1;
  met=-1; metPhi=-1;
  metNoMu = -1.;

  mH = 0; mX = 0; ctau = 0;
  Flag2_all = false;
  gWPt = 0.0;
  gLepId = 0;
  gLepPt = 0.; gLepPhi = 0.; gLepEta = 0.; gLepE = 0.;
  //leptons
  nLeptons = 0;
  for( int i = 0; i < N_MAX_LEPTONS; i++ )
  {
    lepE[i]      = -999.;
    lepPt[i]     = -999.;
    lepEta[i]    = -999.;
    lepPhi[i]    = -999.;
    lepPdgId[i]  = -999;


  }
  //jets
  nJets = 0;
  for( int i = 0; i < N_MAX_JETS; i++ )
  {
    jetE[i]      = -999.;
    jetPt[i]     = -999.;
    jetEta[i]    = -999.;
    jetPhi[i]    = -999.;
  }
  //Z-candidate
  MT = -999.;

  //
  // for(int i = 0; i <NTriggersMAX; i++){
  //   HLTDecision[i] = false;
  // }
  HLT_PFMET120_PFMHT120_IDTight = false;
  HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = false;
  HLT_PFMET140_PFMHT140_IDTight = false;
  HLT_PFMET120_PFMHT120_IDTight_PFHT60 = false;
  HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 =false;
  HLT_PFMETNoMu140_PFMHTNoMu140_IDTight = false;
  HLT_IsoMu27 = false;
  HLT_IsoEle = false;

  METTrigger = false;
};

void LiteLiteTreeMuonSystem::InitTree()
{
  assert(tree_);
  InitVariables();

  tree_->SetBranchAddress("runNum",      &runNum);
  tree_->SetBranchAddress("lumiSec",     &lumiSec);
  tree_->SetBranchAddress("evtNum",      &evtNum);
  tree_->SetBranchAddress("category",    &category);
  tree_->SetBranchAddress("mX",      &mX);
  tree_->SetBranchAddress("mH",      &mH);
  tree_->SetBranchAddress("ctau",      &ctau);

  tree_->SetBranchAddress("npv",         &npv);
  tree_->SetBranchAddress("npu",         &npu);
  tree_->SetBranchAddress("weight",      &weight);



    tree_->SetBranchAddress("gLepId",      &gLepId);
    tree_->SetBranchAddress("gLepPt",      &gLepPt);
    tree_->SetBranchAddress("gLepPhi",      &gLepPhi);
    tree_->SetBranchAddress("gLepE",      &gLepE);
    tree_->SetBranchAddress("gLepEta",      &gLepEta);
  tree_->SetBranchAddress("pileupWeight",      &pileupWeight);
  tree_->SetBranchAddress("Flag2_all",      &Flag2_all);

  tree_->SetBranchAddress("gWPt",      &gWPt);

  tree_->SetBranchAddress("rho",         &rho);
  tree_->SetBranchAddress("met",         &met);
  tree_->SetBranchAddress("metNoMu",         &metNoMu);


  tree_->SetBranchAddress("metPhi",      &metPhi);


  tree_->SetBranchAddress("nLeptons",    &nLeptons);
  tree_->SetBranchAddress("lepE",        lepE);
  tree_->SetBranchAddress("lepPt",       lepPt);
  tree_->SetBranchAddress("lepEta",      lepEta);
  tree_->SetBranchAddress("lepPhi",      lepPhi);
  tree_->SetBranchAddress("lepPdgId",  lepPdgId);

  //jets
  tree_->SetBranchAddress("nJets",     &nJets);
  tree_->SetBranchAddress("jetE",      jetE);
  tree_->SetBranchAddress("jetPt",     jetPt);
  tree_->SetBranchAddress("jetEta",    jetEta);
  tree_->SetBranchAddress("jetPhi",    jetPhi);



  //Z-candidate

  tree_->SetBranchAddress("MT", &MT);
  // triggers
  // tree_->SetBranchAddress("HLTDecision",   HLTDecision);
  tree_->SetBranchAddress("METTrigger",   &METTrigger);
  tree_->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight",   &HLT_PFMET120_PFMHT120_IDTight);
  tree_->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",   &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
  tree_->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight",   &HLT_PFMET140_PFMHT140_IDTight);
  tree_->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60",   &HLT_PFMET120_PFMHT120_IDTight_PFHT60);
  tree_->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",   &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
  tree_->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",   &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
  tree_->SetBranchAddress("HLT_IsoMu27",   &HLT_IsoMu27);
  tree_->SetBranchAddress("HLT_IsoEle",   &HLT_IsoEle);


};

void LiteLiteTreeMuonSystem::LoadTree(const char* file)
{
  f_ = TFile::Open(file);
  assert(f_);
  tree_ = dynamic_cast<TTree*>(f_->Get("MuonSystem"));
  InitTree();
  assert(tree_);
};

void LiteLiteTreeMuonSystem::CreateTree()
{
  tree_ = new TTree("MuonSystem","MuonSystem");
  f_ = 0;

  tree_->Branch("runNum",      &runNum,     "runNum/i");      // event run number
  tree_->Branch("lumiSec",     &lumiSec,    "lumiSec/i");     // event lumi section
  tree_->Branch("evtNum",      &evtNum,     "evtNum/i");      // event number
  tree_->Branch("mH",      &mH,     "mH/I");      // event number
  tree_->Branch("mX",      &mX,     "mX/I");      // event number
  tree_->Branch("ctau",      &ctau,     "ctau/I");      // event number

  tree_->Branch("category",    &category,   "category/i");    // dilepton category
  tree_->Branch("npv",         &npv,        "npv/i");         // number of primary vertices
  tree_->Branch("npu",         &npu,        "npu/i");         // number of in-time PU events (MC)
  tree_->Branch("weight",      &weight,     "weight/F");

  tree_->Branch("pileupWeight",      &pileupWeight,     "pileupWeight/F");

    tree_->Branch("gLepId",      &gLepId,     "gLepId/I");      // phi(MET)
    tree_->Branch("gLepPt",      &gLepPt,     "gLepPt/F");      // phi(MET)
    tree_->Branch("gLepE",      &gLepE,     "gLepE/F");      // phi(MET)
    tree_->Branch("gLepEta",      &gLepEta,     "gLepEta/F");      // phi(MET)
    tree_->Branch("gLepPhi",      &gLepPhi,     "gLepPhi/F");      // phi(MET)

  tree_->Branch("gWPt",         &gWPt,        "gWPt/F");


  tree_->Branch("rho",         &rho,        "rho/F");
  tree_->Branch("met",         &met,        "met/F");         // MET
  tree_->Branch("metNoMu",         &metNoMu,        "metNoMu/F");         // MET


  tree_->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)

  // //leptons
  tree_->Branch("nLeptons",  &nLeptons, "nLeptons/I");
  tree_->Branch("lepE",      lepE,      "lepE[nLeptons]/F");
  tree_->Branch("lepPt",     lepPt,     "lepPt[nLeptons]/F");
  tree_->Branch("lepEta",    lepEta,    "lepEta[nLeptons]/F");
  tree_->Branch("lepPhi",    lepPhi,    "lepPhi[nLeptons]/F");
  tree_->Branch("lepPdgId",  lepPdgId,  "lepPdgId[nLeptons]/I");

    //jets
    tree_->Branch("nJets",     &nJets,    "nJets/I");
    tree_->Branch("jetE",      jetE,      "jetE[nJets]/F");
    tree_->Branch("jetPt",     jetPt,     "jetPt[nJets]/F");
    tree_->Branch("jetEta",    jetEta,    "jetEta[nJets]/F");
    tree_->Branch("jetPhi",    jetPhi,    "jetPhi[nJets]/F");

  tree_->Branch("Flag2_all",      &Flag2_all,  "Flag2_all/O");

  tree_->Branch("MT",      &MT,  "MT/F");
  //Z-candidate

  tree_->Branch("HLT_PFMET120_PFMHT120_IDTight",   HLT_PFMET120_PFMHT120_IDTight, "HLT_PFMET120_PFMHT120_IDTight/O");
  tree_->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",   HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight/O");
  tree_->Branch("HLT_PFMET140_PFMHT140_IDTight",   HLT_PFMET140_PFMHT140_IDTight, "HLT_PFMET140_PFMHT140_IDTight/O");
  tree_->Branch("HLT_PFMET120_PFMHT120_IDTight_PFHT60",   HLT_PFMET120_PFMHT120_IDTight_PFHT60, "HLT_PFMET120_PFMHT120_IDTight_PFHT60/O");
  tree_->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",   HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60/O");
  tree_->Branch("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",   HLT_PFMETNoMu140_PFMHTNoMu140_IDTight, "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight/O");
  tree_->Branch("HLT_IsoMu27",   HLT_IsoMu27, "HLT_IsoMu27/O");
  tree_->Branch("HLT_IsoEle",   HLT_IsoEle, "HLT_IsoEle/O");



  // tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[982]/O"); //hardcoded
  tree_->Branch("METTrigger", METTrigger, "METTrigger/O"); //hardcoded

};
