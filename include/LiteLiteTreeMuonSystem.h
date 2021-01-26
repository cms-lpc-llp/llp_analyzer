// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef LiteLiteTreeMuonSystem_H
#define LiteLiteTreeMuonSystem_H

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 200
#define N_MAX_CSCRECHITS 5000
#define NTriggersMAX 982 // Number of trigger in the .dat file
#define N_CSC_CUT 20
#define JET_PT_CUT 10
#define MUON_PT_CUT 20
#define N_MAX_GPARTICLES 2000

#include <iostream>
#include <string>
#include <sys/stat.h>
#include "assert.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TRandom.h"
#include "FactorizedJetCorrector.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "SimpleJetResolution.h"
#include "BTagCalibrationStandalone.h"
#include "TTree.h"
#include "DBSCAN.h"

#include "RazorAnalyzer.h"

#include "RazorHelper.h"

class LiteLiteTreeMuonSystem
{

public:
  LiteLiteTreeMuonSystem();
  ~LiteLiteTreeMuonSystem();
  // LiteLiteTreeMuonSystem::LiteLiteTreeMuonSystem()
  // {
  //   InitVariables();
  // };
  // LiteLiteTreeMuonSystem::~LiteLiteTreeMuonSystem()
  // {
  //   if (f_) f_->Close();
  // };
  TTree *tree_;
  TFile *f_;

  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  category;
  UInt_t  npv, npu;
  float rho, weight;
  float pileupWeight;


  float met, metPhi, metNoMu;
  bool Flag_HBHENoiseFilter, Flag_HBHEIsoNoiseFilter, Flag_BadPFMuonFilter, Flag_globalSuperTightHalo2016Filter,
  Flag_CSCTightHaloFilter, Flag_BadChargedCandidateFilter, Flag_eeBadScFilter, Flag_goodVertices, Flag_ecalBadCalibFilter, Flag_all;
  int mH, mX, ctau;
  bool Flag2_all;
  bool  METTrigger;
  float  gWPt;


  int gLepId;
  float gLepPt, gLepPhi, gLepEta, gLepE;
  // bool HLTDecision[NTriggersMAX];
bool HLT_PFMET120_PFMHT120_IDTight;
bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
bool HLT_PFMET140_PFMHT140_IDTight;
bool HLT_PFMET120_PFMHT120_IDTight_PFHT60;
bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
bool HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;
bool HLT_IsoMu27;
bool HLT_IsoEle;
int nJets;
float jetE[N_MAX_JETS];
float jetPt[N_MAX_JETS];
float jetEta[N_MAX_JETS];
float jetPhi[N_MAX_JETS];

  //leptons
  int nLeptons;
  float lepE[N_MAX_LEPTONS];
  float lepPt[N_MAX_LEPTONS];
  float lepEta[N_MAX_LEPTONS];
  float lepPhi[N_MAX_LEPTONS];
  int  lepPdgId[N_MAX_LEPTONS];
  float MT;

  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();

};
#endif
