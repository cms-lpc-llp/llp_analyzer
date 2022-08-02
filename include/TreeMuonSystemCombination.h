// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef TreeMuonSystemCombination_H
#define TreeMuonSystemCombination_H

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 200
#define N_MAX_CSCRECHITS 5000
#define N_MAX_DTRECHITS 20000
#define NTriggersMAX 982 // Number of trigger in the .dat file
#define N_CSC_CUT 20
#define JET_PT_CUT 10
#define MUON_PT_CUT 20
#define N_MAX_GPARTICLES 5000

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

class TreeMuonSystemCombination
{

public:
  TreeMuonSystemCombination();
  ~TreeMuonSystemCombination();
  // TreeMuonSystemCombination::TreeMuonSystemCombination()
  // {
  //   InitVariables();
  // };
  // TreeMuonSystemCombination::~TreeMuonSystemCombination()
  // {
  //   if (f_) f_->Close();
  // };
  TTree *tree_;
  TFile *f_;

  UInt_t  runNum, lumiSec, evtNum, MC_condition;
  UInt_t  category;
  UInt_t  npv, npu;
  float rho, weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  float higgsPtWeight;
  float higgsPtWeightSys[9];
  float scaleWeights[9];
  int ZCategory;
  float metSF;
  float lepOverallSF;
  float sf_facScaleUp, sf_facScaleDown, sf_renScaleUp, sf_renScaleDown, sf_facRenScaleUp, sf_facRenScaleDown;
  bool EE_prefiring;
  float met, metPhi, HT, jetMet_dPhi, jetMet_dPhiMin,jetMet_dPhiMin4, metJESUp, metPhiJESUp, metJESDown, metPhiJESDown, metJESDownSF, metJESUpSF, metNoMu, metEENoise, metPhiEENoise, metHEM, metPhiHEM, metHEMXYCorr, metPhiHEMXYCorr;
  float metXYCorr, metPhiXYCorr, metEENoiseXYCorr, metPhiEENoiseXYCorr;
  bool Flag_HBHENoiseFilter, Flag_HBHEIsoNoiseFilter, Flag_BadPFMuonFilter, Flag_globalSuperTightHalo2016Filter,
  Flag_CSCTightHaloFilter, Flag_BadChargedCandidateFilter, Flag_eeBadScFilter, Flag_goodVertices, Flag_ecalBadCalibFilter, Flag_all;
  int mH, mX, ctau;

  bool Flag2_HBHENoiseFilter, Flag2_HBHEIsoNoiseFilter, Flag2_BadPFMuonFilter, Flag2_globalSuperTightHalo2016Filter,
  Flag2_globalTightHalo2016Filter, Flag2_BadChargedCandidateFilter, Flag2_EcalDeadCellTriggerPrimitiveFilter,
  Flag2_ecalBadCalibFilter, Flag2_eeBadScFilter, Flag2_all;

  float gWPt;
  int gLepId;
  float gLepPt, gLepPhi, gLepEta, gLepE;
  int nGenParticle;
  int gParticleId[N_MAX_GPARTICLES];
  int gParticleStatus[N_MAX_GPARTICLES];
  int gParticleMotherId[N_MAX_GPARTICLES];
  int gParticleMotherIndex[N_MAX_GPARTICLES];
  float gParticlePt[N_MAX_GPARTICLES];
  float gParticleEta[N_MAX_GPARTICLES];
  float gParticlePhi[N_MAX_GPARTICLES];
  float gParticleE[N_MAX_GPARTICLES];
  float gParticleProdVertexX[N_MAX_GPARTICLES];
  float gParticleProdVertexY[N_MAX_GPARTICLES];
  float gParticleProdVertexZ[N_MAX_GPARTICLES];


  float genMetPtTrue;
  float genMetPhiTrue;
  float genMetPtCalo;
  float genMetPhiCalo;


  int nGenJets;
  float genJetE[N_MAX_GPARTICLES];
  float genJetPt[N_MAX_GPARTICLES];
  float genJetEta[N_MAX_GPARTICLES];
  float genJetPhi[N_MAX_GPARTICLES];
  float genJetMET[N_MAX_GPARTICLES];

  // int         nCosmicTwoLegClusters;
  // float         cosmicTwoLegClusterChi2Reduced[N_MAX_CSC];
  // int         cosmicTwoLegClusterSize[N_MAX_CSC];
  // float         cosmicTwoLegClusterChi2[N_MAX_CSC];
  // int         cosmicTwoLegCluster1NSegmentStation1[N_MAX_CSC];
  // int         cosmicTwoLegCluster1NSegmentStation2[N_MAX_CSC];
  // int         cosmicTwoLegCluster1NSegmentStation3[N_MAX_CSC];
  // int         cosmicTwoLegCluster1NSegmentStation4[N_MAX_CSC];
  // int         cosmicTwoLegCluster2NSegmentStation1[N_MAX_CSC];
  // int         cosmicTwoLegCluster2NSegmentStation2[N_MAX_CSC];
  // int         cosmicTwoLegCluster2NSegmentStation3[N_MAX_CSC];
  // int         cosmicTwoLegCluster2NSegmentStation4[N_MAX_CSC];
  // int         cosmicTwoLegCluster2NStation[N_MAX_CSC];
  // int         cosmicTwoLegCluster1NStation[N_MAX_CSC];
  // int         cosmicTwoLegCluster2Size[N_MAX_CSC];
  // int         cosmicTwoLegCluster1Size[N_MAX_CSC];
  // int         cosmicTwoLegCluster2Index[N_MAX_CSC];
  // int         cosmicTwoLegCluster1Index[N_MAX_CSC];
  // float       cosmicTwoLegCluster_m_xz[N_MAX_CSC];
  // float       cosmicTwoLegCluster_c_xz[N_MAX_CSC];
  // float       cosmicTwoLegCluster_m_yz[N_MAX_CSC];
  // float       cosmicTwoLegCluster_c_yz[N_MAX_CSC];
  //
  // int         nCosmicOneLegClusters;
  // float         cosmicOneLegClusterChi2Reduced[N_MAX_CSC];
  // int         cosmicOneLegClusterSize[N_MAX_CSC];
  // float         cosmicOneLegClusterChi2[N_MAX_CSC];
  // int         cosmicOneLegClusterNSegmentStation1[N_MAX_CSC];
  // int         cosmicOneLegClusterNSegmentStation2[N_MAX_CSC];
  // int         cosmicOneLegClusterNSegmentStation3[N_MAX_CSC];
  // int         cosmicOneLegClusterNSegmentStation4[N_MAX_CSC];
  // int         cosmicOneLegClusterNStation[N_MAX_CSC];
  // float       cosmicOneLegCluster_m_xz[N_MAX_CSC];
  // float       cosmicOneLegCluster_c_xz[N_MAX_CSC];
  // float       cosmicOneLegCluster_m_yz[N_MAX_CSC];
  // float       cosmicOneLegCluster_c_yz[N_MAX_CSC];

  // vector<float>  cosmicOneLegClusterdtSegX[N_MAX_CSC];
  // vector<float>  cosmicOneLegClusterdtSegY[N_MAX_CSC];
  // vector<float>  cosmicOneLegClusterdtSegZ[N_MAX_CSC];

  //csc
  int           nCscRechits;
  int           nEarlyCscRechits;
  int           nLateCscRechits;
  int           nEarly2CscRechits;
  int           nLate2CscRechits;
  int           nCscRings;
  int           nCscPositiveYRechits;
  int           nCscNegativeYRechits;
  float         cscPosTpeak;
  float         cscNegTpeak;
  int           nCscRechitsChamberPlus11;
  int           nCscRechitsChamberPlus12;
  int           nCscRechitsChamberPlus13;
  int           nCscRechitsChamberPlus21;
  int           nCscRechitsChamberPlus22;
  int           nCscRechitsChamberPlus31;
  int           nCscRechitsChamberPlus32;
  int           nCscRechitsChamberPlus41;
  int           nCscRechitsChamberPlus42;

  int           nCscRechitsChamberMinus11;
  int           nCscRechitsChamberMinus12;
  int           nCscRechitsChamberMinus13;
  int           nCscRechitsChamberMinus21;
  int           nCscRechitsChamberMinus22;
  int           nCscRechitsChamberMinus31;
  int           nCscRechitsChamberMinus32;
  int           nCscRechitsChamberMinus41;
  int           nCscRechitsChamberMinus42;
  //
  //
  // int           nCscRechitsChamberPlus11[36];
  // int           nCscRechitsChamberPlus12[36];
  // int           nCscRechitsChamberPlus13[36];
  // int           nCscRechitsChamberPlus21[36];
  // int           nCscRechitsChamberPlus22[36];
  // int           nCscRechitsChamberPlus31[36];
  // int           nCscRechitsChamberPlus32[36];
  // int           nCscRechitsChamberPlus41[36];
  // int           nCscRechitsChamberPlus42[36];
  //
  // int           nCscRechitsChamberMinus11[36];
  // int           nCscRechitsChamberMinus12[36];
  // int           nCscRechitsChamberMinus13[36];
  // int           nCscRechitsChamberMinus21[36];
  // int           nCscRechitsChamberMinus22[36];
  // int           nCscRechitsChamberMinus31[36];
  // int           nCscRechitsChamberMinus32[36];
  // int           nCscRechitsChamberMinus41[36];
  // int           nCscRechitsChamberMinus42[36];

  int           nDTRechits;
  int           nDtRings;
  int           nDtStations25;
  int           nDtWheels25;
  int           nDTRechitsStation1;
  int           nDTRechitsStation2;
  int           nDTRechitsStation3;
  int           nDTRechitsStation4;
  int 		nDTRechitsSector[4][5][12];//[NStattion][NWheels][NSector]
  // int 		nDTSegSector[4][5][12];//[NStattion][NWheels][NSector]
  int           nDTRechitsWheelMinus2;
  int           nDTRechitsWheelMinus1;
  int           nDTRechitsWheel0;
  int           nDTRechitsWheelPlus1;
  int           nDTRechitsWheelPlus2;
  int           nDTPositiveYRechits;
  int           nDTNegativeYRechits;
  int           nDTRechitsChamberMinus12;
  int           nDTRechitsChamberMinus11;
  int           nDTRechitsChamber10;
  int           nDTRechitsChamberPlus11;
  int           nDTRechitsChamberPlus12;
  int           nDTRechitsChamberMinus22;
  int           nDTRechitsChamberMinus21;
  int           nDTRechitsChamber20;
  int           nDTRechitsChamberPlus21;
  int           nDTRechitsChamberPlus22;
  int           nDTRechitsChamberMinus32;
  int           nDTRechitsChamberMinus31;
  int           nDTRechitsChamber30;
  int           nDTRechitsChamberPlus31;
  int           nDTRechitsChamberPlus32;
  int           nDTRechitsChamberMinus42;
  int           nDTRechitsChamberMinus41;
  int           nDTRechitsChamber40;
  int           nDTRechitsChamberPlus41;
  int           nDTRechitsChamberPlus42;
  float         dtRechitsX[N_MAX_DTRECHITS];   //[nCsc]
  float         dtRechitsY[N_MAX_DTRECHITS];   //[nCsc]
  float         dtRechitsZ[N_MAX_DTRECHITS];   //[nCsc]
  float         dtRechitsPhi[N_MAX_DTRECHITS];   //[nCsc]
  float         dtRechitsEta[N_MAX_DTRECHITS];   //[nCsc]
  int         dtRechitsStation[N_MAX_DTRECHITS];   //[nCsc]
  int         dtRechitsWheel[N_MAX_DTRECHITS];   //[nCsc]
  int         dtRechitsClusterId[N_MAX_DTRECHITS];
  int           nRpc;
  float         rpcPhi[N_MAX_CSCRECHITS];   //[nCsc]
  float         rpcEta[N_MAX_CSCRECHITS];   //[nCsc]
  bool         rpc_RE12[N_MAX_CSCRECHITS];   //[nCsc]
  bool         rpc_RB1[N_MAX_CSCRECHITS];   //[nCsc]

  int           nDtSeg;
  float         dtSegX[N_MAX_CSCRECHITS];   //[nCsc]
  float         dtSegY[N_MAX_CSCRECHITS];   //[nCsc]
  float         dtSegZ[N_MAX_CSCRECHITS];   //[nCsc]

  float         dtSegPhi[N_MAX_CSCRECHITS];   //[nCsc]
  float         dtSegEta[N_MAX_CSCRECHITS];   //[nCsc]
  int         dtSegStation[N_MAX_CSCRECHITS];   //[nCsc]
  int         dtSegWheel[N_MAX_CSCRECHITS];   //[nCsc]

  float         cscRechitsPhi[N_MAX_CSCRECHITS];   //[nCsc]
  float         cscRechitsEta[N_MAX_CSCRECHITS];   //[nCsc]
  float         cscRechitsX[N_MAX_CSCRECHITS];   //[nCsc]
  float         cscRechitsY[N_MAX_CSCRECHITS];   //[nCsc]
  float         cscRechitsZ[N_MAX_CSCRECHITS];   //[nCsc]
  int         cscRechitsQuality[N_MAX_CSCRECHITS];   //[nCsc]

  float         cscRechitsTpeak[N_MAX_CSCRECHITS];   //[nCsc]
  float         cscRechitsTwire[N_MAX_CSCRECHITS];   //[nCsc]
  int         cscRechitsChamber[N_MAX_CSCRECHITS];   //[nCsc]
  int         cscRechitsStation[N_MAX_CSCRECHITS];   //[nCsc]
  int           cscRechitsClusterId[N_MAX_CSCRECHITS];
  int           cscRechitsCluster2Id[N_MAX_CSCRECHITS];

  int           nDtRechitClusters;
  float           dtRechitClusterMaxDPhi[N_MAX_CSC];
  bool           dtRechitClusterOverlap[N_MAX_CSC];
  int           dtRechitClusterMaxDPhi_index[N_MAX_CSC];
  bool          dtRechitClusterNoiseVeto[N_MAX_CSC];
  int           dtRechitClusterNSegStation1[N_MAX_CSC];
  int           dtRechitClusterNSegStation2[N_MAX_CSC];
  int           dtRechitClusterNSegStation3[N_MAX_CSC];
  int           dtRechitClusterNSegStation4[N_MAX_CSC];
  int           dtRechitClusterNOppositeSegStation1[N_MAX_CSC];
  int           dtRechitClusterNOppositeSegStation2[N_MAX_CSC];
  int           dtRechitClusterNOppositeSegStation3[N_MAX_CSC];
  int           dtRechitClusterNOppositeSegStation4[N_MAX_CSC];
  int           dtRechitClusterNSegmentStation1[N_MAX_CSC];
  int           dtRechitClusterNSegmentStation2[N_MAX_CSC];
  int           dtRechitClusterNSegmentStation3[N_MAX_CSC];
  int           dtRechitClusterNSegmentStation4[N_MAX_CSC];
  bool         dtRechitClusterZLep1[N_MAX_CSC];
  bool         dtRechitClusterZLep2[N_MAX_CSC];
   int         dtRechitClusterZLep1Id[N_MAX_CSC];
   int         dtRechitClusterZLep2Id[N_MAX_CSC];
  bool         dtRechitClusterZLep1TightId[N_MAX_CSC];
  bool         dtRechitClusterZLep1LooseIso[N_MAX_CSC];
  bool         dtRechitClusterZLep1TightIso[N_MAX_CSC];
  bool         dtRechitClusterZLep1VTightIso[N_MAX_CSC];
  bool         dtRechitClusterZLep1VVTightIso[N_MAX_CSC];
  bool         dtRechitClusterZLep1Tag[N_MAX_CSC];

  bool         dtRechitClusterZLep2TightId[N_MAX_CSC];
  bool         dtRechitClusterZLep2LooseIso[N_MAX_CSC];
  bool         dtRechitClusterZLep2TightIso[N_MAX_CSC];
  bool         dtRechitClusterZLep2VTightIso[N_MAX_CSC];
  bool         dtRechitClusterZLep2VVTightIso[N_MAX_CSC];
  bool         dtRechitClusterZLep2Tag[N_MAX_CSC];

  int         dtRechitCluster_match_gParticle_Id[N_MAX_CSC];
  float       dtRechitCluster_match_gParticle_Pt[N_MAX_CSC];
  float       dtRechitCluster_match_gParticle_Eta[N_MAX_CSC];
  float       dtRechitCluster_match_gParticle_Phi[N_MAX_CSC];
  float       dtRechitCluster_match_gParticle_E[N_MAX_CSC];
  int         dtRechitCluster_match_gParticle_Status[N_MAX_CSC];
  int         dtRechitCluster_match_gParticle_MotherId[N_MAX_CSC];
  float       dtRechitCluster_match_gParticle_deltaR[N_MAX_CSC];

  int           dtRechitCluster_match_MB1hits_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_MB1hits_0p5[N_MAX_CSC];
  int           dtRechitCluster_match_MB1hits_cosmics_plus[N_MAX_CSC];
  int           dtRechitCluster_match_MB1hits_cosmics_minus[N_MAX_CSC];
  int           dtRechitCluster_match_MB1Seg_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_MB1Seg_0p5[N_MAX_CSC];
  int           dtRechitCluster_match_dtSeg_0p5[N_MAX_CSC];
  float         dtRechitCluster_match_dtSegTime_0p5[N_MAX_CSC];
  int           dtRechitCluster_match_dtSeg_0p4[N_MAX_CSC];
  float         dtRechitCluster_match_dtSegTime_0p4[N_MAX_CSC];
  float         dtRechitCluster_match_dtSegTimeSpread_0p5[N_MAX_CSC];
  float         dtRechitCluster_match_dtSegTimeSpread_0p4[N_MAX_CSC];

  int           dtRechitCluster_match_dtSeg_sameStation_0p5[N_MAX_CSC];
  float         dtRechitCluster_match_dtSegTime_sameStation_0p5[N_MAX_CSC];
  int           dtRechitCluster_match_dtSeg_sameStation_0p4[N_MAX_CSC];
  float         dtRechitCluster_match_dtSegTime_sameStation_0p4[N_MAX_CSC];
  float         dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5[N_MAX_CSC];
  float         dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_RPCBx_dPhi0p5[N_MAX_CSC];
  int           dtRechitCluster_match_RB1_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_RB1_dPhi0p5[N_MAX_CSC];
  //
  // vector<vector<float> > dtRechitCluster_match_dtSegT_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_dtSegX_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_dtSegY_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_dtSegZ_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_dtSegPhi_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_dtSegEta_dR0p4;
  // vector<vector<int> > dtRechitCluster_match_dtSegWheel_dR0p4;
  // vector<vector<int> > dtRechitCluster_match_dtSegStation_dR0p4;
  //
  //
  // vector<vector<float> > dtRechitCluster_match_dtSegT_sameStation_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_dtSegX_sameStation_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_dtSegY_sameStation_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_dtSegZ_sameStation_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_dtSegPhi_sameStation_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_dtSegEta_sameStation_dR0p4;
  // vector<vector<int> > dtRechitCluster_match_dtSegWheel_sameStation_dR0p4;
  // vector<vector<int> > dtRechitCluster_match_dtSegStation_sameStation_dR0p4;
  //
  //
  // vector<vector<float> > dtRechitCluster_match_RPCBx_dPhi0p5;
  // vector<vector<float> > dtRechitCluster_match_RPCX_dPhi0p5;
  // vector<vector<float> > dtRechitCluster_match_RPCY_dPhi0p5;
  // vector<vector<float> > dtRechitCluster_match_RPCZ_dPhi0p5;
  // vector<vector<float> > dtRechitCluster_match_RPCPhi_dPhi0p5;
  // vector<vector<float> > dtRechitCluster_match_RPCEta_dPhi0p5;
  // vector<vector<int> > dtRechitCluster_match_RPCRing_dPhi0p5;
  // vector<vector<int> > dtRechitCluster_match_RPCLayer_dPhi0p5;
  // vector<vector<int> > dtRechitCluster_match_RPCSector_dPhi0p5;
  //
  //
  //
  // vector<vector<float> > dtRechitCluster_match_RPCBx_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_RPCX_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_RPCY_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_RPCZ_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_RPCPhi_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_RPCEta_dR0p4;
  // vector<vector<int> > dtRechitCluster_match_RPCRing_dR0p4;
  // vector<vector<int> > dtRechitCluster_match_RPCLayer_dR0p4;
  // vector<vector<int> > dtRechitCluster_match_RPCSector_dR0p4;
  //
  // vector<vector<float> > dtRechitCluster_match_RPCBx_sameStation_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_RPCX_sameStation_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_RPCY_sameStation_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_RPCZ_sameStation_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_RPCPhi_sameStation_dR0p4;
  // vector<vector<float> > dtRechitCluster_match_RPCEta_sameStation_dR0p4;
  // vector<vector<int> > dtRechitCluster_match_RPCRing_sameStation_dR0p4;
  // vector<vector<int> > dtRechitCluster_match_RPCLayer_sameStation_dR0p4;
  // vector<vector<int> > dtRechitCluster_match_RPCSector_sameStation_dR0p4;


  float         dtRechitCluster_match_RPCTime_dPhi0p5[N_MAX_CSC];
  float         dtRechitCluster_match_RPCTimeSpread_dPhi0p5[N_MAX_CSC];
  int         dtRechitCluster_match_RPChits_dPhi0p5[N_MAX_CSC];
  float         dtRechitCluster_match_RPCTime_dR0p4[N_MAX_CSC];
  float         dtRechitCluster_match_RPCTimeSpread_dR0p4[N_MAX_CSC];
  int         dtRechitCluster_match_RPChits_dR0p4[N_MAX_CSC];
  float         dtRechitCluster_match_RPCTime_sameStation_dR0p4[N_MAX_CSC];
  float         dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[N_MAX_CSC];
  int         dtRechitCluster_match_RPChits_sameStation_dR0p4[N_MAX_CSC];

  bool          dtRechitCluster_match_gLLP[N_MAX_CSC];
  int           dtRechitCluster_match_gLLP_index[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_minDeltaR[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_eta[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_phi[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_decay_r[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_decay_x[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_decay_y[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_decay_z[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_ctau[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_beta[N_MAX_CSC];
  bool         dtRechitCluster_match_gLLP_csc[N_MAX_CSC];
  bool         dtRechitCluster_match_gLLP_dt[N_MAX_CSC];
  int         dtRechitCluster_match_gLLP_multiplicity[N_MAX_CSC];
  int         dtRechitCluster_match_gLLP_EM_multiplicity[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_e[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_pt[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_EMFracE[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_EMFracEz[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_EMFracP[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_EMFracPz[N_MAX_CSC];
  bool          dtRechitCluster_match_gLLP_daughterKaon[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_visE[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_visEz[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_visP[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_visPz[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_lepdPhi[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_daughter0_deltaR[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_daughter1_deltaR[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_daughter2_deltaR[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_daughter3_deltaR[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_other_daughter_deltaR[N_MAX_CSC];
  int         dtRechitCluster_match_gLLP_other_daughter_index[N_MAX_CSC];

  float         dtRechitCluster_match_gLLP_other_eta[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_other_phi[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_other_decay_r[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_other_decay_x[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_other_decay_y[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_other_decay_z[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_other_ctau[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_other_beta[N_MAX_CSC];
  bool         dtRechitCluster_match_gLLP_other_csc[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_other_e[N_MAX_CSC];
  float         dtRechitCluster_match_gLLP_other_pt[N_MAX_CSC];
  float         dtRechitClusterX[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterY[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterZ[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterTime[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterTimeTotal[N_MAX_CSC];
  float         dtRechitClusterTimeWire[N_MAX_CSC];
  float         dtRechitClusterTimeWirePruned[N_MAX_CSC];
  int         dtRechitClusterWheel[N_MAX_CSC];

  float         dtRechitClusterTimeSpread[N_MAX_CSC];
  float         dtRechitClusterTimeWireSpread[N_MAX_CSC];
  float         dtRechitClusterTimeTotalSpread[N_MAX_CSC];
  float         dtRechitClusterTimeTotalSpreadPruned[N_MAX_CSC];

  float         dtRechitClusterGenMuonDeltaR[N_MAX_CSC];
  float         dtRechitClusterMajorAxis[N_MAX_CSC];
  float         dtRechitClusterMinorAxis[N_MAX_CSC];
  float         dtRechitClusterXSpread[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterYSpread[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterZSpread[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterXYSpread[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterRSpread[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterDeltaRSpread[N_MAX_CSC];


  float         dtRechitClusterEtaPhiSpread[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterEtaSpread[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterPhiSpread[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterEta[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterPhi[N_MAX_CSC];   //[nCsc]
  int           dtRechitClusterSize[N_MAX_CSC];
  int           dtRechitClusterNoiseHit[N_MAX_CSC];
  int           dtRechitClusterNoiseHitStation1[N_MAX_CSC];
  int           dtRechitClusterNoiseHitStation2[N_MAX_CSC];
  int           dtRechitClusterNoiseHitStation3[N_MAX_CSC];
  int           dtRechitClusterNoiseHitStation4[N_MAX_CSC];
  float         dtRechitClusterMe11Ratio[N_MAX_CSC];
  float         dtRechitClusterMe12Ratio[N_MAX_CSC];


  float         dtRechitClusterMaxStationRatio[N_MAX_CSC];   //[nCsc]
  int           dtRechitClusterMaxStation[N_MAX_CSC];   //[nCsc]
  int           dtRechitClusterNStation[N_MAX_CSC];
  int           dtRechitClusterNStation5[N_MAX_CSC];
  int           dtRechitClusterNStation10[N_MAX_CSC];
  int           dtRechitClusterNStation10perc[N_MAX_CSC];
  float          dtRechitClusterAvgStation[N_MAX_CSC];
  float          dtRechitClusterAvgStation5[N_MAX_CSC];
  float          dtRechitClusterAvgStation10[N_MAX_CSC];
  float          dtRechitClusterAvgStation10perc[N_MAX_CSC];
  float         dtRechitClusterMaxChamberRatio[N_MAX_CSC];   //[nCsc]
  int           dtRechitClusterMaxChamber[N_MAX_CSC];   //[nCsc]
  int           dtRechitClusterNChamber[N_MAX_CSC];

  float         dtRechitClusterJetVetoElectronEnergyFraction[N_MAX_CSC];
  float         dtRechitClusterJetVetoPhotonEnergyFraction[N_MAX_CSC];
  float         dtRechitClusterJetVetoChargedHadronEnergyFraction[N_MAX_CSC];
  float         dtRechitClusterJetVetoNeutralHadronEnergyFraction[N_MAX_CSC];
  float         dtRechitClusterJetVetoMuonEnergyFraction[N_MAX_CSC];
  float         dtRechitClusterJetVetoEta[N_MAX_CSC];
  float         dtRechitClusterJetVetoPhi[N_MAX_CSC];
  float         dtRechitClusterJetVetoPt[N_MAX_CSC];
  float         dtRechitClusterJetVetoPtJESUp[N_MAX_CSC];
  float         dtRechitClusterJetVetoPtJESDown[N_MAX_CSC];
  float         dtRechitClusterJetVetoE[N_MAX_CSC];
  float         dtRechitClusterTightJetVetoEta[N_MAX_CSC];
  float         dtRechitClusterTightJetVetoPhi[N_MAX_CSC];
  float         dtRechitClusterTightJetVetoPt[N_MAX_CSC];
  float         dtRechitClusterTightJetVetoPtJESUp[N_MAX_CSC];
  float         dtRechitClusterTightJetVetoPtJESDown[N_MAX_CSC];

  float         dtRechitClusterGenJetVetoPt[N_MAX_CSC];
  float         dtRechitClusterGenJetVetoE[N_MAX_CSC];
  float         dtRechitClusterMuonVetoPt[N_MAX_CSC];
  float         dtRechitClusterMuonVetoE[N_MAX_CSC];
  float         dtRechitClusterMuonVetoPhi[N_MAX_CSC];
  float         dtRechitClusterMuonVetoEta[N_MAX_CSC];
  float         dtRechitClusterJetVetoPt_0p6[N_MAX_CSC];
  float         dtRechitClusterJetVetoPt_0p8[N_MAX_CSC];
  float         dtRechitClusterJetVetoE_0p6[N_MAX_CSC];
  float         dtRechitClusterJetVetoE_0p8[N_MAX_CSC];
  float         dtRechitClusterMuonVetoPt_0p6[N_MAX_CSC];
  float         dtRechitClusterMuonVetoPt_0p8[N_MAX_CSC];
  float         dtRechitClusterMuonVetoE_0p6[N_MAX_CSC];
  float         dtRechitClusterMuonVetoE_0p8[N_MAX_CSC];
  bool          dtRechitClusterMuonVetoLooseIso[N_MAX_CSC];
  bool          dtRechitClusterMuonVetoTightIso[N_MAX_CSC];
  bool          dtRechitClusterMuonVetoVTightIso[N_MAX_CSC];
  bool          dtRechitClusterMuonVetoVVTightIso[N_MAX_CSC];
  bool          dtRechitClusterMuonVetoTightId[N_MAX_CSC];
  bool          dtRechitClusterMuonVetoLooseId[N_MAX_CSC];
  bool          dtRechitClusterMuonVetoGlobal[N_MAX_CSC];

  bool          dtRechitClusterMuonVetoIso[N_MAX_CSC];
  float         dtRechitClusterIsoMuonVetoPt[N_MAX_CSC];
  float         dtRechitClusterIsoMuonVetoE[N_MAX_CSC];
  float         dtRechitClusterIsoMuonVetoPhi[N_MAX_CSC];
  float         dtRechitClusterIsoMuonVetoEta[N_MAX_CSC];
  float         dtRechitClusterGenMuonVetoPt[N_MAX_CSC];
  float         dtRechitClusterGenMuonVetoE[N_MAX_CSC];
  float         dtRechitClusterGenMuonVetoProdX[N_MAX_CSC];
  float         dtRechitClusterGenMuonVetoProdY[N_MAX_CSC];
  float         dtRechitClusterGenMuonVetoProdZ[N_MAX_CSC];
  int         dtRechitClusterGenMuonVetoLLPIndex[N_MAX_CSC];
  float         dtRechitClusterGenMuonVetoLLPDist[N_MAX_CSC];
  int           dtRechitClusterMuonVetoType[N_MAX_CSC];


  int           dtRechitClusterNRechitChamberPlus11[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberPlus12[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberPlus13[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberPlus21[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberPlus22[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberPlus31[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberPlus32[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberPlus41[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberPlus42[N_MAX_CSC];

  int           dtRechitClusterNRechitChamberMinus11[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberMinus12[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberMinus13[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberMinus21[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberMinus22[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberMinus31[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberMinus32[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberMinus41[N_MAX_CSC];
  int           dtRechitClusterNRechitChamberMinus42[N_MAX_CSC];
  float         dtRechitClusterMet_dPhi[N_MAX_CSC];
  float         dtRechitClusterMetJESUp_dPhi[N_MAX_CSC];
  float         dtRechitClusterMetJESDown_dPhi[N_MAX_CSC];
  float         dtRechitClusterMetXYCorr_dPhi[N_MAX_CSC];

  float     dtRechitClusterMetHEM_dPhi[N_MAX_CSC];
  float     dtRechitClusterMetHEMXYCorr_dPhi[N_MAX_CSC];
  float     dtRechitClusterMetEENoise_dPhi[N_MAX_CSC];
  float     dtRechitClusterMetEENoiseXYCorr_dPhi[N_MAX_CSC];
  float     dtRechitClusterMetJesUp_dPhi[N_MAX_CSC];
  float     dtRechitClusterMetJesDown_dPhi[N_MAX_CSC];

  int           dtRechitClusterNLayersChamberPlus11[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberPlus12[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberPlus13[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberPlus21[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberPlus22[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberPlus31[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberPlus32[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberPlus41[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberPlus42[N_MAX_CSC];

  int           dtRechitClusterNLayersChamberMinus11[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberMinus12[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberMinus13[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberMinus21[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberMinus22[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberMinus31[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberMinus32[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberMinus41[N_MAX_CSC];
  int           dtRechitClusterNLayersChamberMinus42[N_MAX_CSC];



  int           nCscRechitClusters;

  int           cscRechitCluster_match_Me1112_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_Me1112_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_Me1112_0p8[N_MAX_CSC];
  int           cscRechitCluster_match_Me11_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_Me11_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_Me11_0p8[N_MAX_CSC];
  int           cscRechitCluster_match_Me12_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_Me12_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_Me12_0p8[N_MAX_CSC];
  int           cscRechitCluster_match_cscRechits_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_cscSeg_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_ME11Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_ME12Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_cscSeg_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_ME11Seg_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_ME12Seg_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_dtRechits_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_dtRechits_phi0p2[N_MAX_CSC];
  int           cscRechitCluster_match_MB1_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_dtRechits_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_MB1_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_dtSeg_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_MB1Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_dtSeg_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_MB1Seg_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_RB1_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_RE12_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_RB1_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_RE12_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_highEta_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_highEta_0p6[N_MAX_CSC];
  int           cscRechitCluster_match_highEta_0p8[N_MAX_CSC];
  float           cscRechitCluster_match_cluster_dR[N_MAX_CSC];
  int           cscRechitCluster_match_cluster_index[N_MAX_CSC];
  bool          cscRechitCluster_match_gParticle[N_MAX_CSC];
  float         cscRechitCluster_match_gParticle_minDeltaR[N_MAX_CSC];

  int           cscRechitCluster_match_gParticleMotherId[N_MAX_CSC];

    float           cscRechitCluster_match_gParticle_deltaR[N_MAX_CSC];
    float           cscRechitCluster_match_gParticlePt[N_MAX_CSC];
    float           cscRechitCluster_match_gParticleEta[N_MAX_CSC];
    float           cscRechitCluster_match_gParticlePhi[N_MAX_CSC];
    float           cscRechitCluster_match_gParticleE[N_MAX_CSC];
    int           cscRechitCluster_match_gParticleStatus[N_MAX_CSC];
    int         cscRechitCluster_match_gParticleId[N_MAX_CSC];
    int         cscRechitCluster_match_gParticleIndex[N_MAX_CSC];
    float         cscRechitCluster_match_gParticleProdVertexX[N_MAX_CSC];
    float         cscRechitCluster_match_gParticleProdVertexY[N_MAX_CSC];
    float         cscRechitCluster_match_gParticleProdVertexZ[N_MAX_CSC];
  bool          cscRechitCluster_match_gLLP[N_MAX_CSC];
  int           cscRechitCluster_match_gLLP_index[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_minDeltaR[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_eta[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_phi[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_decay_r[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_decay_x[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_decay_y[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_decay_z[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_ctau[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_beta[N_MAX_CSC];
  bool         cscRechitCluster_match_gLLP_csc[N_MAX_CSC];
  bool         cscRechitCluster_match_gLLP_dt[N_MAX_CSC];
  int         cscRechitCluster_match_gLLP_multiplicity[N_MAX_CSC];
  int         cscRechitCluster_match_gLLP_EM_multiplicity[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_e[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_pt[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_EMFracE[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_EMFracEz[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_EMFracP[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_EMFracPz[N_MAX_CSC];
  bool          cscRechitCluster_match_gLLP_daughterKaon[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_visE[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_visEz[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_visP[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_visPz[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_lepdPhi[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_daughter0_deltaR[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_daughter1_deltaR[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_daughter2_deltaR[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_daughter3_deltaR[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_other_daughter_deltaR[N_MAX_CSC];
  int         cscRechitCluster_match_gLLP_other_daughter_index[N_MAX_CSC];

  float         cscRechitCluster_match_gLLP_other_eta[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_other_phi[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_other_decay_r[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_other_decay_x[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_other_decay_y[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_other_decay_z[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_other_ctau[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_other_beta[N_MAX_CSC];
  bool         cscRechitCluster_match_gLLP_other_csc[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_other_e[N_MAX_CSC];
  float         cscRechitCluster_match_gLLP_other_pt[N_MAX_CSC];
  float         cscRechitClusterX[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterY[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterZ[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterTime[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterTimeTotal[N_MAX_CSC];
  float         cscRechitClusterTimeWeighted[N_MAX_CSC];

  float         cscRechitClusterTimeSpread[N_MAX_CSC];
  float         cscRechitClusterTimeSpreadWeighted[N_MAX_CSC];
  float         cscRechitClusterTimeSpreadWeightedAll[N_MAX_CSC];

  float         cscRechitClusterGenMuonDeltaR[N_MAX_CSC];
  float         cscRechitClusterMajorAxis[N_MAX_CSC];
  float         cscRechitClusterMinorAxis[N_MAX_CSC];
  float         cscRechitClusterXSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterYSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterZSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterXYSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterRSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterDeltaRSpread[N_MAX_CSC];


  float         cscRechitClusterXYSpread_r1p2[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterXSpread_r1p2[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterYSpread_r1p2[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterPhiSpread_r1p2[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterEtaSpread_r1p2[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterEtaPhiSpread_r1p2[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterRSpread_r1p2[N_MAX_CSC];   //[nCsc]

  float         cscRechitClusterEtaPhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterEtaSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterPhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterEta[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterPhi[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterSize[N_MAX_CSC];
  float         cscRechitClusterMe11Ratio[N_MAX_CSC];
  float         cscRechitClusterMe12Ratio[N_MAX_CSC];


  float         cscRechitClusterMaxStationRatio[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterMaxStation[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterNStation[N_MAX_CSC];
  int           cscRechitClusterNStation5[N_MAX_CSC];
  int           cscRechitClusterNStation10[N_MAX_CSC];
  int           cscRechitClusterNStation10perc[N_MAX_CSC];
  float          cscRechitClusterAvgStation[N_MAX_CSC];
  float          cscRechitClusterAvgStation5[N_MAX_CSC];
  float          cscRechitClusterAvgStation10[N_MAX_CSC];
  float          cscRechitClusterAvgStation10perc[N_MAX_CSC];
  float         cscRechitClusterMaxChamberRatio[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterMaxChamber[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterNChamber[N_MAX_CSC];

  float         cscRechitClusterJetVetoElectronEnergyFraction[N_MAX_CSC];
  float         cscRechitClusterJetVetoPhotonEnergyFraction[N_MAX_CSC];
  float         cscRechitClusterJetVetoChargedHadronEnergyFraction[N_MAX_CSC];
  float         cscRechitClusterJetVetoNeutralHadronEnergyFraction[N_MAX_CSC];
  float         cscRechitClusterJetVetoMuonEnergyFraction[N_MAX_CSC];
  float         cscRechitClusterJetVetoEta[N_MAX_CSC];
  float         cscRechitClusterJetVetoPhi[N_MAX_CSC];
  float         cscRechitClusterJetVetoPt[N_MAX_CSC];
  float         cscRechitClusterJetVetoPtJESUp[N_MAX_CSC];
  float         cscRechitClusterJetVetoPtJESDown[N_MAX_CSC];
  float         cscRechitClusterTightJetVetoEta[N_MAX_CSC];
  float         cscRechitClusterTightJetVetoPhi[N_MAX_CSC];
  float         cscRechitClusterTightJetVetoPt[N_MAX_CSC];
  float         cscRechitClusterTightJetVetoPtJESUp[N_MAX_CSC];
  float         cscRechitClusterTightJetVetoPtJESDown[N_MAX_CSC];

  float         cscRechitClusterJetVetoE[N_MAX_CSC];
  float         cscRechitClusterGenJetVetoPt[N_MAX_CSC];
  float         cscRechitClusterGenJetVetoE[N_MAX_CSC];
  float         cscRechitClusterMuonVetoPt[N_MAX_CSC];
  float         cscRechitClusterMuonVetoE[N_MAX_CSC];
  float         cscRechitClusterMuonVetoPhi[N_MAX_CSC];
  float         cscRechitClusterMuonVetoEta[N_MAX_CSC];
  float         cscRechitClusterJetVetoPt_0p6[N_MAX_CSC];
  float         cscRechitClusterJetVetoPt_0p8[N_MAX_CSC];
  float         cscRechitClusterJetVetoE_0p6[N_MAX_CSC];
  float         cscRechitClusterJetVetoE_0p8[N_MAX_CSC];
  float         cscRechitClusterMuonVetoPt_0p6[N_MAX_CSC];
  float         cscRechitClusterMuonVetoPt_0p8[N_MAX_CSC];
  float         cscRechitClusterMuonVetoE_0p6[N_MAX_CSC];
  float         cscRechitClusterMuonVetoE_0p8[N_MAX_CSC];
  bool         cscRechitClusterZLep1[N_MAX_CSC];
  bool         cscRechitClusterZLep2[N_MAX_CSC];
  int         cscRechitClusterZLep1Id[N_MAX_CSC];
  int         cscRechitClusterZLep2Id[N_MAX_CSC];
  bool         cscRechitClusterZLep1TightId[N_MAX_CSC];
  bool         cscRechitClusterZLep1LooseIso[N_MAX_CSC];
  bool         cscRechitClusterZLep1TightIso[N_MAX_CSC];
  bool         cscRechitClusterZLep1VTightIso[N_MAX_CSC];
  bool         cscRechitClusterZLep1VVTightIso[N_MAX_CSC];
  bool         cscRechitClusterZLep1Tag[N_MAX_CSC];

  bool         cscRechitClusterZLep2TightId[N_MAX_CSC];
  bool         cscRechitClusterZLep2LooseIso[N_MAX_CSC];
  bool         cscRechitClusterZLep2TightIso[N_MAX_CSC];
  bool         cscRechitClusterZLep2VTightIso[N_MAX_CSC];
  bool         cscRechitClusterZLep2VVTightIso[N_MAX_CSC];
  bool         cscRechitClusterZLep2Tag[N_MAX_CSC];
  bool          cscRechitClusterMuonVetoLooseIso[N_MAX_CSC];
  bool          cscRechitClusterMuonVetoTightIso[N_MAX_CSC];
  bool          cscRechitClusterMuonVetoVTightIso[N_MAX_CSC];
  bool          cscRechitClusterMuonVetoVVTightIso[N_MAX_CSC];
  bool          cscRechitClusterMuonVetoTightId[N_MAX_CSC];
  bool          cscRechitClusterMuonVetoLooseId[N_MAX_CSC];
  bool          cscRechitClusterMuonVetoGlobal[N_MAX_CSC];

  bool          cscRechitClusterMuonVetoIso[N_MAX_CSC];
  float         cscRechitClusterIsoMuonVetoPt[N_MAX_CSC];
  float         cscRechitClusterIsoMuonVetoE[N_MAX_CSC];
  float         cscRechitClusterIsoMuonVetoPhi[N_MAX_CSC];
  float         cscRechitClusterIsoMuonVetoEta[N_MAX_CSC];
  float         cscRechitClusterGenMuonVetoPt[N_MAX_CSC];
  float         cscRechitClusterGenMuonVetoE[N_MAX_CSC];
  float         cscRechitClusterGenMuonVetoProdX[N_MAX_CSC];
  float         cscRechitClusterGenMuonVetoProdY[N_MAX_CSC];
  float         cscRechitClusterGenMuonVetoProdZ[N_MAX_CSC];
  int         cscRechitClusterGenMuonVetoLLPIndex[N_MAX_CSC];
  float         cscRechitClusterGenMuonVetoLLPDist[N_MAX_CSC];
  int           cscRechitClusterMuonVetoType[N_MAX_CSC];


  int           cscRechitClusterNRechitChamberPlus11[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus12[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus13[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus21[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus22[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus31[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus32[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus41[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberPlus42[N_MAX_CSC];

  int           cscRechitClusterNRechitChamberMinus11[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus12[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus13[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus21[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus22[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus31[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus32[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus41[N_MAX_CSC];
  int           cscRechitClusterNRechitChamberMinus42[N_MAX_CSC];
  float         cscRechitClusterMet_dPhi[N_MAX_CSC];
  float         cscRechitClusterMetJESUp_dPhi[N_MAX_CSC];
  float         cscRechitClusterMetJESDown_dPhi[N_MAX_CSC];
  float         cscRechitClusterMetXYCorr_dPhi[N_MAX_CSC];

  float     cscRechitClusterMetHEM_dPhi[N_MAX_CSC];
  float     cscRechitClusterMetHEMXYCorr_dPhi[N_MAX_CSC];
  float     cscRechitClusterMetEENoise_dPhi[N_MAX_CSC];
  float     cscRechitClusterMetEENoiseXYCorr_dPhi[N_MAX_CSC];
  float     cscRechitClusterMetJesUp_dPhi[N_MAX_CSC];
  float     cscRechitClusterMetJesDown_dPhi[N_MAX_CSC];

  int           cscRechitClusterNLayersChamberPlus11[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberPlus12[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberPlus13[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberPlus21[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberPlus22[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberPlus31[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberPlus32[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberPlus41[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberPlus42[N_MAX_CSC];

  int           cscRechitClusterNLayersChamberMinus11[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberMinus12[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberMinus13[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberMinus21[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberMinus22[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberMinus31[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberMinus32[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberMinus41[N_MAX_CSC];
  int           cscRechitClusterNLayersChamberMinus42[N_MAX_CSC];

  //gLLP
  float gLLP_eta[2];
  float gLLP_phi[2];
  float gLLP_csc[2];
  float gLLP_dt[2];
  float gLLP_beta[2];
  float gLLP_maxMatchedDis[2];
  int gLLP_match_dtRechits[2];
  int gLLP_match_cscRechits[2];


  float gLLP_e[2];
  bool gLLP_daughterKaon[2];
  float gLLP_pt[2];
  float gLLP_lepdPhi[2];
  int gLLP_multiplicity[2];
  int gLLP_multiplicity20[2];

  int gLLP_EM_multiplicity[2];



  float gLLP_ctau[2];
  float gLLP_decay_vertex_r[2];
  float gLLP_decay_vertex_x[2];
  float gLLP_decay_vertex_y[2];
  float gLLP_decay_vertex_z[2];
  float gLLP_EMFracE[2];
  float gLLP_EMFracEz[2];
  float gLLP_EMFracP[2];
  float gLLP_EMFracPz[2];
  float gLLP_visE[2];
  float gLLP_visE20[2];
  float gLLP_visEz[2];
  float gLLP_visP[2];
  float gLLP_visPz[2];
  float gLLP_deltaR;
  float gLLP_daughter_deltaR[2];
  float gLLP_daughter_pt[4];
  int   gLLP_daughter_id[4];
  float gLLP_daughter_eta[4];
  float gLLP_daughter_phi[4];
  float gLLP_daughter_e[4];
  float gLLP_daughter_mass[4];
  float gHiggsPt;
  float gHiggsEta;
  float gHiggsPhi;
  float gHiggsE;

  float gLLP_match_jet_minDeltaR[2];
  int gLLP_match_jet_index[2];
  float gLLP_match_jet_pt[2];
  //leptons

  int nLeptons;
  float lepE[N_MAX_LEPTONS];
  float lepPt[N_MAX_LEPTONS];
  float lepEta[N_MAX_LEPTONS];
  float lepPhi[N_MAX_LEPTONS];
  int  lepPdgId[N_MAX_LEPTONS];
  float lepDZ[N_MAX_LEPTONS];
  float lepTriggerSF[N_MAX_LEPTONS];
  float lepTightIdSF[N_MAX_LEPTONS];
  float lepLooseIdSF[N_MAX_LEPTONS];
  float lepTightIsoSF[N_MAX_LEPTONS];
  float lepLooseIsoSF[N_MAX_LEPTONS];
  float lepTriggerMCEfficiency[N_MAX_LEPTONS];
  float lepTightIdMCEfficiency[N_MAX_LEPTONS];
  float lepLooseIdMCEfficiency[N_MAX_LEPTONS];
  float lepTightIsoMCEfficiency[N_MAX_LEPTONS];
  float lepLooseIsoMCEfficiency[N_MAX_LEPTONS];
  float lepEff[N_MAX_LEPTONS];
  float lepSF[N_MAX_LEPTONS];
  bool lepTag[N_MAX_LEPTONS];


  int nGlobalMuons;
  float GlobalMuonPt[N_MAX_LEPTONS];
  float GlobalMuonEta[N_MAX_LEPTONS];
  float GlobalMuonPhi[N_MAX_LEPTONS];
  bool GlobalMuonLooseId[N_MAX_LEPTONS];

  // bool lepLoosePassId[N_MAX_LEPTONS];
  // bool lepMediumPassId[N_MAX_LEPTONS];
  // bool lepTightPassId[N_MAX_LEPTONS];
  bool lepPassVetoId[N_MAX_LEPTONS];
  bool lepFromZ[N_MAX_LEPTONS];
  bool lepLooseId[N_MAX_LEPTONS];
  bool lepTightId[N_MAX_LEPTONS];
  bool lepGlobal[N_MAX_LEPTONS];
  bool lepPassLooseIso[N_MAX_LEPTONS];
  bool lepPassTightIso[N_MAX_LEPTONS];
  bool lepPassVTightIso[N_MAX_LEPTONS];
  bool lepPassVVTightIso[N_MAX_LEPTONS];

  //Z-candidate
  float MT;
  float ZMass1;
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
  float jetPtJESUp[N_MAX_JETS];
  float jetPtJESDown[N_MAX_JETS];
  float jetEJESUp[N_MAX_JETS];
  float jetEJESDown[N_MAX_JETS];
  float JecUnc[N_MAX_JETS];

  float jetChargedEMEnergyFraction[N_MAX_JETS];
  float jetNeutralEMEnergyFraction[N_MAX_JETS];
  float jetChargedHadronEnergyFraction[N_MAX_JETS];
  float jetNeutralHadronEnergyFraction[N_MAX_JETS];



  bool jetPassMuFrac[N_MAX_JETS];
  float jet_match_llp_minDeltaR[N_MAX_JETS];
  int jet_match_llp_index[N_MAX_JETS];
  float jet_match_llp_pt[N_MAX_JETS];

  float jet_match_genJet_minDeltaR[N_MAX_JETS];
  int jet_match_genJet_index[N_MAX_JETS];
  float jet_match_genJet_pt[N_MAX_JETS];
  // bool jetLoosePassId[N_MAX_JETS];
  bool jetPassId[N_MAX_JETS];
  bool jetTightPassId[N_MAX_JETS];
  bool HLTDecision[NTriggersMAX];
  bool METTrigger;
  bool METNoMuTrigger;
  UInt_t wzevtNum,trig, trig_lepId, trig_lepId_dijet; //number of events that pass each criteria




  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();

};
#endif
