// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef HNLMuonSystemTree_H
#define HNLMuonSystemTree_H

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 200
#define N_MAX_CSCRECHITS 5000
#define NTriggersMAX 1201 // Number of trigger in the .dat file
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

class HNLMuonSystemTree
{

public:
  HNLMuonSystemTree();
  ~HNLMuonSystemTree();

  TTree *tree_;
  TFile *f_;

  UInt_t  runNum, lumiSec, evtNum;
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
  float gParticlePt[N_MAX_GPARTICLES];
  float gParticleEta[N_MAX_GPARTICLES];
  float gParticlePhi[N_MAX_GPARTICLES];
  float gParticleE[N_MAX_GPARTICLES];

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
  float         dtRechitsPhi[N_MAX_CSCRECHITS];   //[nCsc]
  float         dtRechitsEta[N_MAX_CSCRECHITS];   //[nCsc]
  int         dtRechitsStation[N_MAX_CSCRECHITS];   //[nCsc]
  int         dtRechitsWheel[N_MAX_CSCRECHITS];   //[nCsc]
  int           nRpc;
  float         rpcPhi[N_MAX_CSCRECHITS];   //[nCsc]
  float         rpcEta[N_MAX_CSCRECHITS];   //[nCsc]
  bool         rpc_RE12[N_MAX_CSCRECHITS];   //[nCsc]
  bool         rpc_RB1[N_MAX_CSCRECHITS];   //[nCsc]

  int           nDtSeg;
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

    int           nCscClusters;
    bool          cscCluster_match_gLLP[N_MAX_CSC];
    int           cscCluster_match_gLLP_index[N_MAX_CSC];
    int           cscCluster_match_gLLP_minDeltaR[N_MAX_CSC];
    float         cscClusterX[N_MAX_CSC];   //[nCsc]
    float         cscClusterY[N_MAX_CSC];   //[nCsc]
    float         cscClusterZ[N_MAX_CSC];   //[nCsc]
    float         cscClusterTime[N_MAX_CSC];   //[nCsc]
    float         cscClusterTimeSpread[N_MAX_CSC];
    float         cscClusterGenMuonDeltaR[N_MAX_CSC];
    float         cscClusterMajorAxis[N_MAX_CSC];
    float         cscClusterMinorAxis[N_MAX_CSC];
    float         cscClusterXSpread[N_MAX_CSC];   //[nCsc]
    float         cscClusterYSpread[N_MAX_CSC];   //[nCsc]
    float         cscClusterZSpread[N_MAX_CSC];   //[nCsc]
    float         cscClusterEtaPhiSpread[N_MAX_CSC];   //[nCsc]
    float         cscClusterEtaSpread[N_MAX_CSC];   //[nCsc]
    float         cscClusterPhiSpread[N_MAX_CSC];   //[nCsc]
    float         cscClusterEta[N_MAX_CSC];   //[nCsc]
    float         cscClusterPhi[N_MAX_CSC];   //[nCsc]
    int           cscClusterSize[N_MAX_CSC];
    float         cscClusterMe11Ratio[N_MAX_CSC];
    float         cscClusterMe12Ratio[N_MAX_CSC];

    float         cscClusterVertexR[N_MAX_CSC];   //[nCsc]
    float         cscClusterVertexZ[N_MAX_CSC];   //[nCsc]
    int           cscClusterVertexN[N_MAX_CSC];   //[nCsc]
    int           cscClusterVertexN1[N_MAX_CSC];   //[nCsc]
    int           cscClusterVertexN5[N_MAX_CSC];   //[nCsc]
    int           cscClusterVertexN10[N_MAX_CSC];   //[nCsc]
    int           cscClusterVertexN15[N_MAX_CSC];   //[nCsc]
    int           cscClusterVertexN20[N_MAX_CSC];   //[nCsc]
    float         cscClusterVertexChi2[N_MAX_CSC];   //[nCsc]
    float         cscClusterVertexDis[N_MAX_CSC];   //[nCsc]
    float         cscClusterMaxStationRatio[N_MAX_CSC];   //[nCsc]
    int           cscClusterMaxStation[N_MAX_CSC];   //[nCsc]
    int           cscClusterNStation[N_MAX_CSC];
    int           cscClusterNStation10[N_MAX_CSC];
    float           cscClusterAvgStation[N_MAX_CSC];
    float           cscClusterAvgStation10[N_MAX_CSC];

    float         cscClusterMaxChamberRatio[N_MAX_CSC];   //[nCsc]
    int           cscClusterMaxChamber[N_MAX_CSC];   //[nCsc]
    int           cscClusterNChamber[N_MAX_CSC];
    float         cscClusterJetVetoPt[N_MAX_CSC];
    float         cscClusterJetVetoE[N_MAX_CSC];
    float         cscClusterMuonVetoPt[N_MAX_CSC];
    float         cscClusterMuonVetoE[N_MAX_CSC];
    int           cscClusterNSegmentChamberPlus11[N_MAX_CSC];
    int           cscClusterNSegmentChamberPlus12[N_MAX_CSC];
    int           cscClusterNSegmentChamberPlus13[N_MAX_CSC];
    int           cscClusterNSegmentChamberPlus21[N_MAX_CSC];
    int           cscClusterNSegmentChamberPlus22[N_MAX_CSC];
    int           cscClusterNSegmentChamberPlus31[N_MAX_CSC];
    int           cscClusterNSegmentChamberPlus32[N_MAX_CSC];
    int           cscClusterNSegmentChamberPlus41[N_MAX_CSC];
    int           cscClusterNSegmentChamberPlus42[N_MAX_CSC];

    int           cscClusterNSegmentChamberMinus11[N_MAX_CSC];
    int           cscClusterNSegmentChamberMinus12[N_MAX_CSC];
    int           cscClusterNSegmentChamberMinus13[N_MAX_CSC];
    int           cscClusterNSegmentChamberMinus21[N_MAX_CSC];
    int           cscClusterNSegmentChamberMinus22[N_MAX_CSC];
    int           cscClusterNSegmentChamberMinus31[N_MAX_CSC];
    int           cscClusterNSegmentChamberMinus32[N_MAX_CSC];
    int           cscClusterNSegmentChamberMinus41[N_MAX_CSC];
    int           cscClusterNSegmentChamberMinus42[N_MAX_CSC];
    float         cscClusterMet_dPhi[N_MAX_CSC];

  int           nCscSegClusters;

  int           cscSegCluster_match_Me1112_0p4[N_MAX_CSC];
  int           cscSegCluster_match_Me1112_0p6[N_MAX_CSC];
  int           cscSegCluster_match_Me1112_0p8[N_MAX_CSC];
  int           cscSegCluster_match_Me11_0p4[N_MAX_CSC];
  int           cscSegCluster_match_Me11_0p6[N_MAX_CSC];
  int           cscSegCluster_match_Me11_0p8[N_MAX_CSC];
  int           cscSegCluster_match_Me12_0p4[N_MAX_CSC];
  int           cscSegCluster_match_Me12_0p6[N_MAX_CSC];
  int           cscSegCluster_match_Me12_0p8[N_MAX_CSC];
  int           cscSegCluster_match_cscRechits_0p4[N_MAX_CSC];
  int           cscSegCluster_match_cscSeg_0p4[N_MAX_CSC];
  int           cscSegCluster_match_ME11Seg_0p4[N_MAX_CSC];
  int           cscSegCluster_match_ME12Seg_0p4[N_MAX_CSC];
  int           cscSegCluster_match_cscSeg_0p6[N_MAX_CSC];
  int           cscSegCluster_match_ME11Seg_0p6[N_MAX_CSC];
  int           cscSegCluster_match_ME12Seg_0p6[N_MAX_CSC];
  int           cscSegCluster_match_dtRechits_0p4[N_MAX_CSC];
  int           cscSegCluster_match_MB1_0p4[N_MAX_CSC];
  int           cscSegCluster_match_dtRechits_0p6[N_MAX_CSC];
  int           cscSegCluster_match_MB1_0p6[N_MAX_CSC];
  int           cscSegCluster_match_dtSeg_0p4[N_MAX_CSC];
  int           cscSegCluster_match_MB1Seg_0p4[N_MAX_CSC];
  int           cscSegCluster_match_dtSeg_0p6[N_MAX_CSC];
  int           cscSegCluster_match_MB1Seg_0p6[N_MAX_CSC];
  int           cscSegCluster_match_RB1_0p4[N_MAX_CSC];
  int           cscSegCluster_match_RE12_0p4[N_MAX_CSC];
  int           cscSegCluster_match_RB1_0p6[N_MAX_CSC];
  int           cscSegCluster_match_RE12_0p6[N_MAX_CSC];
  int           cscSegCluster_match_highEta_0p4[N_MAX_CSC];
  int           cscSegCluster_match_highEta_0p6[N_MAX_CSC];
  int           cscSegCluster_match_highEta_0p8[N_MAX_CSC];
  float           cscSegCluster_match_cluster_dR[N_MAX_CSC];
  int           cscSegCluster_match_cluster_index[N_MAX_CSC];
  bool          cscSegCluster_match_gParticle[N_MAX_CSC];
  float         cscSegCluster_match_gParticle_minDeltaR[N_MAX_CSC];
  int           cscSegCluster_match_gParticle_index[N_MAX_CSC];
  int           cscSegCluster_match_gParticle_id[N_MAX_CSC];
  float         cscSegCluster_match_gParticle_eta[N_MAX_CSC];
  float         cscSegCluster_match_gParticle_phi[N_MAX_CSC];
  float         cscSegCluster_match_gParticle_E[N_MAX_CSC];
  float         cscSegCluster_match_gParticle_pt[N_MAX_CSC];
  int           cscSegCluster_match_gParticle_MotherId[N_MAX_CSC];

  bool          cscSegCluster_match_gLLP[N_MAX_CSC];
  int           cscSegCluster_match_gLLP_index[N_MAX_CSC];
  float         cscSegCluster_match_gLLP_minDeltaR[N_MAX_CSC];
  float         cscSegCluster_match_gLLP_eta[N_MAX_CSC];
  float         cscSegCluster_match_gLLP_phi[N_MAX_CSC];
  float         cscSegCluster_match_gLLP_decay_r[N_MAX_CSC];
  float         cscSegCluster_match_gLLP_decay_x[N_MAX_CSC];
  float         cscSegCluster_match_gLLP_decay_y[N_MAX_CSC];
  float         cscSegCluster_match_gLLP_decay_z[N_MAX_CSC];
  float         cscSegCluster_match_gLLP_ctau[N_MAX_CSC];
  float         cscSegCluster_match_gLLP_beta[N_MAX_CSC];
  bool         cscSegCluster_match_gLLP_csc[N_MAX_CSC];
  float         cscSegClusterX[N_MAX_CSC];   //[nCsc]
  float         cscSegClusterY[N_MAX_CSC];   //[nCsc]
  float         cscSegClusterZ[N_MAX_CSC];   //[nCsc]
  float         cscSegClusterTime[N_MAX_CSC];   //[nCsc]
  float         cscSegClusterTimeTotal[N_MAX_CSC];

  float         cscSegClusterTimeSpread[N_MAX_CSC];
  float         cscSegClusterGenMuonDeltaR[N_MAX_CSC];
  float         cscSegClusterMajorAxis[N_MAX_CSC];
  float         cscSegClusterMinorAxis[N_MAX_CSC];
  float         cscSegClusterXSpread[N_MAX_CSC];   //[nCsc]
  float         cscSegClusterYSpread[N_MAX_CSC];   //[nCsc]
  float         cscSegClusterZSpread[N_MAX_CSC];   //[nCsc]
  float         cscSegClusterXYSpread[N_MAX_CSC];   //[nCsc]
  float         cscSegClusterRSpread[N_MAX_CSC];   //[nCsc]


  float         cscSegClusterXYSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]
  float         cscSegClusterXSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]
  float         cscSegClusterYSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]
  float         cscSegClusterPhiSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]
  float         cscSegClusterEtaSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]
  float         cscSegClusterEtaPhiSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]
  float         cscSegClusterRSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]

    float         cscSegClusterXYSpread_corr2[N_MAX_CSC];   //[nCsc]
    float         cscSegClusterXSpread_corr2[N_MAX_CSC];   //[nCsc],
    float         cscSegClusterYSpread_corr2[N_MAX_CSC];   //[nCsc]
    float         cscSegClusterPhiSpread_corr2[N_MAX_CSC];   //[nCsc]
    float         cscSegClusterEtaSpread_corr2[N_MAX_CSC];   //[nCsc]
    float         cscSegClusterEtaPhiSpread_corr2[N_MAX_CSC];   //[nCsc]
    float         cscSegClusterRSpread_corr2[N_MAX_CSC];   //[nCsc]


  float         cscSegClusterEtaPhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscSegClusterEtaSpread[N_MAX_CSC];   //[nCsc]
  float         cscSegClusterPhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscSegClusterEta[N_MAX_CSC];   //[nCsc]
  float         cscSegClusterPhi[N_MAX_CSC];   //[nCsc]
  int           cscSegClusterSize[N_MAX_CSC];
  float         cscSegClusterMe11Ratio[N_MAX_CSC];
  float         cscSegClusterMe12Ratio[N_MAX_CSC];


  float         cscSegClusterMaxStationRatio[N_MAX_CSC];   //[nCsc]
  int           cscSegClusterMaxStation[N_MAX_CSC];   //[nCsc]
  int           cscSegClusterNStation[N_MAX_CSC];
  int           cscSegClusterNStation5[N_MAX_CSC];
  int           cscSegClusterNStation10perc[N_MAX_CSC];
  float          cscSegClusterAvgStation[N_MAX_CSC];
  float          cscSegClusterAvgStation5[N_MAX_CSC];
  float           cscSegClusterAvgStation10perc[N_MAX_CSC];
  float         cscSegClusterMaxChamberRatio[N_MAX_CSC];   //[nCsc]
  int           cscSegClusterMaxChamber[N_MAX_CSC];   //[nCsc]
  int           cscSegClusterNChamber[N_MAX_CSC];
  float         cscSegClusterJetVetoPt[N_MAX_CSC];
  float         cscSegClusterJetVetoE[N_MAX_CSC];
  float         cscSegClusterMuonVetoPt[N_MAX_CSC];
  float         cscSegClusterMuonVetoE[N_MAX_CSC];
  float         cscSegClusterMuonVetoPhi[N_MAX_CSC];
  float         cscSegClusterMuonVetoEta[N_MAX_CSC];
  float         cscSegClusterJetVetoPt_0p6[N_MAX_CSC];
  float         cscSegClusterJetVetoPt_0p8[N_MAX_CSC];
  float         cscSegClusterJetVetoE_0p6[N_MAX_CSC];
  float         cscSegClusterJetVetoE_0p8[N_MAX_CSC];
  float         cscSegClusterMuonVetoPt_0p6[N_MAX_CSC];
  float         cscSegClusterMuonVetoPt_0p8[N_MAX_CSC];
  float         cscSegClusterMuonVetoE_0p6[N_MAX_CSC];
  float         cscSegClusterMuonVetoE_0p8[N_MAX_CSC];
  bool         cscSegClusterZLep1[N_MAX_CSC];
  bool         cscSegClusterZLep2[N_MAX_CSC];
  int         cscSegClusterZLep1Id[N_MAX_CSC];
  int         cscSegClusterZLep2Id[N_MAX_CSC];
  bool         cscSegClusterZLep1TightId[N_MAX_CSC];
  bool         cscSegClusterZLep1LooseIso[N_MAX_CSC];
  bool         cscSegClusterZLep1TightIso[N_MAX_CSC];
  bool         cscSegClusterZLep1VTightIso[N_MAX_CSC];
  bool         cscSegClusterZLep1VVTightIso[N_MAX_CSC];
  bool         cscSegClusterZLep2TightId[N_MAX_CSC];
  bool         cscSegClusterZLep2LooseIso[N_MAX_CSC];
  bool         cscSegClusterZLep2TightIso[N_MAX_CSC];
  bool         cscSegClusterZLep2VTightIso[N_MAX_CSC];
  bool         cscSegClusterZLep2VVTightIso[N_MAX_CSC];
  bool          cscSegClusterMuonVetoLooseIso[N_MAX_CSC];
  bool          cscSegClusterMuonVetoTightIso[N_MAX_CSC];
  bool          cscSegClusterMuonVetoVTightIso[N_MAX_CSC];
  bool          cscSegClusterMuonVetoVVTightIso[N_MAX_CSC];
  bool          cscSegClusterMuonVetoTightId[N_MAX_CSC];


  bool          cscSegClusterMuonVetoIso[N_MAX_CSC];
  float         cscSegClusterIsoMuonVetoPt[N_MAX_CSC];
  float         cscSegClusterIsoMuonVetoE[N_MAX_CSC];
  float         cscSegClusterIsoMuonVetoPhi[N_MAX_CSC];
  float         cscSegClusterIsoMuonVetoEta[N_MAX_CSC];
  float         cscSegClusterGenMuonVetoPt[N_MAX_CSC];
  float         cscSegClusterGenMuonVetoE[N_MAX_CSC];
  int           cscSegClusterNRechitChamberPlus11[N_MAX_CSC];
  int           cscSegClusterNRechitChamberPlus12[N_MAX_CSC];
  int           cscSegClusterNRechitChamberPlus13[N_MAX_CSC];
  int           cscSegClusterNRechitChamberPlus21[N_MAX_CSC];
  int           cscSegClusterNRechitChamberPlus22[N_MAX_CSC];
  int           cscSegClusterNRechitChamberPlus31[N_MAX_CSC];
  int           cscSegClusterNRechitChamberPlus32[N_MAX_CSC];
  int           cscSegClusterNRechitChamberPlus41[N_MAX_CSC];
  int           cscSegClusterNRechitChamberPlus42[N_MAX_CSC];

  int           cscSegClusterNRechitChamberMinus11[N_MAX_CSC];
  int           cscSegClusterNRechitChamberMinus12[N_MAX_CSC];
  int           cscSegClusterNRechitChamberMinus13[N_MAX_CSC];
  int           cscSegClusterNRechitChamberMinus21[N_MAX_CSC];
  int           cscSegClusterNRechitChamberMinus22[N_MAX_CSC];
  int           cscSegClusterNRechitChamberMinus31[N_MAX_CSC];
  int           cscSegClusterNRechitChamberMinus32[N_MAX_CSC];
  int           cscSegClusterNRechitChamberMinus41[N_MAX_CSC];
  int           cscSegClusterNRechitChamberMinus42[N_MAX_CSC];
  float         cscSegClusterMet_dPhi[N_MAX_CSC];
  float         cscSegClusterMetXYCorr_dPhi[N_MAX_CSC];

  int           nCsc_JetVetoCluster0p4;
  int           nCsc_JetMuonVetoCluster0p4;
  int           nCsc_JetVetoCluster0p4_Me1112Veto;
  int           nCsc_JetMuonVetoCluster0p4_Me1112Veto;

  int           nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto;

  int           nCscRechitClusters;



  int           cscRechitCluster_match_gParticle_id[N_MAX_CSC];
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
  float         cscRechitClusterX[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterY[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterZ[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterTime[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterTimeTotal[N_MAX_CSC];

  float         cscRechitClusterTimeSpread[N_MAX_CSC];
  float         cscRechitClusterGenMuonDeltaR[N_MAX_CSC];
  float         cscRechitClusterMajorAxis[N_MAX_CSC];
  float         cscRechitClusterMinorAxis[N_MAX_CSC];
  float         cscRechitClusterXYSpread[N_MAX_CSC];   //[nCsc]

  float         cscRechitClusterXSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterYSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterZSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterEtaPhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterEtaSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterPhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterEta[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterPhi[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterSize[N_MAX_CSC];
  float         cscRechitClusterMe11Ratio[N_MAX_CSC];
  float         cscRechitClusterMe12Ratio[N_MAX_CSC];


    float         cscRechitClusterXYSpread_corr[N_MAX_CSC];   //[nCsc]
    float         cscRechitClusterXSpread_corr[N_MAX_CSC];   //[nCsc]
    float         cscRechitClusterYSpread_corr[N_MAX_CSC];   //[nCsc]
    float         cscRechitClusterPhiSpread_corr[N_MAX_CSC];   //[nCsc]
    float         cscRechitClusterEtaSpread_corr[N_MAX_CSC];   //[nCsc]
    float         cscRechitClusterEtaPhiSpread_corr[N_MAX_CSC];   //[nCsc]
    float         cscRechitClusterRSpread_corr[N_MAX_CSC];   //[nCsc]




  float         cscRechitClusterMaxStationRatio[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterMaxStation[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterNStation[N_MAX_CSC];
  int           cscRechitClusterNStation5[N_MAX_CSC];
  int           cscRechitClusterNStation10[N_MAX_CSC];
  int           cscRechitClusterNStation10perc[N_MAX_CSC];

  float          cscRechitClusterAvgStation[N_MAX_CSC];
  float           cscRechitClusterAvgStation5[N_MAX_CSC];
  float           cscRechitClusterAvgStation10[N_MAX_CSC];
  float          cscRechitClusterAvgStation10perc[N_MAX_CSC];

  float         cscRechitClusterMaxChamberRatio[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterMaxChamber[N_MAX_CSC];   //[nCsc]
  int           cscRechitClusterNChamber[N_MAX_CSC];
  float         cscRechitClusterJetVetoPt[N_MAX_CSC];
  float         cscRechitClusterJetVetoE[N_MAX_CSC];
  float         cscRechitClusterMuonVetoPt[N_MAX_CSC];
  float         cscRechitClusterMuonVetoE[N_MAX_CSC];
  int           cscRechitCluster_match_MB1Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_RB1_0p4[N_MAX_CSC];
  int           cscRechitCluster_match_RE12_0p4[N_MAX_CSC];

  float         cscRechitClusterGenMuonVetoPt[N_MAX_CSC];
  float         cscRechitClusterGenMuonVetoE[N_MAX_CSC];
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


  int           nCscRechitClusters2;

  int           cscRechitCluster2_match_Me1112_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_Me1112_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_Me1112_0p8[N_MAX_CSC];
  int           cscRechitCluster2_match_Me11_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_Me11_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_Me11_0p8[N_MAX_CSC];
  int           cscRechitCluster2_match_Me12_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_Me12_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_Me12_0p8[N_MAX_CSC];
  int           cscRechitCluster2_match_cscRechits_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_cscSeg_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_ME11Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_ME12Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_cscSeg_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_ME11Seg_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_ME12Seg_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_dtRechits_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_MB1_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_dtRechits_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_MB1_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_dtSeg_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_MB1Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_dtSeg_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_MB1Seg_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_RB1_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_RE12_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_RB1_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_RE12_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_highEta_0p4[N_MAX_CSC];
  int           cscRechitCluster2_match_highEta_0p6[N_MAX_CSC];
  int           cscRechitCluster2_match_highEta_0p8[N_MAX_CSC];
  float           cscRechitCluster2_match_cluster_dR[N_MAX_CSC];
  int           cscRechitCluster2_match_cluster_index[N_MAX_CSC];
  bool          cscRechitCluster2_match_gParticle[N_MAX_CSC];
  float         cscRechitCluster2_match_gParticle_minDeltaR[N_MAX_CSC];
  int           cscRechitCluster2_match_gParticle_index[N_MAX_CSC];
  int           cscRechitCluster2_match_gParticle_id[N_MAX_CSC];
  float         cscRechitCluster2_match_gParticle_eta[N_MAX_CSC];
  float         cscRechitCluster2_match_gParticle_phi[N_MAX_CSC];
  float         cscRechitCluster2_match_gParticle_E[N_MAX_CSC];
  float         cscRechitCluster2_match_gParticle_pt[N_MAX_CSC];
  int           cscRechitCluster2_match_gParticle_MotherId[N_MAX_CSC];

  bool          cscRechitCluster2_match_gLLP[N_MAX_CSC];
  int           cscRechitCluster2_match_gLLP_index[N_MAX_CSC];
  float         cscRechitCluster2_match_gLLP_minDeltaR[N_MAX_CSC];
  float         cscRechitCluster2_match_gLLP_eta[N_MAX_CSC];
  float         cscRechitCluster2_match_gLLP_phi[N_MAX_CSC];
  float         cscRechitCluster2_match_gLLP_decay_r[N_MAX_CSC];
  float         cscRechitCluster2_match_gLLP_decay_x[N_MAX_CSC];
  float         cscRechitCluster2_match_gLLP_decay_y[N_MAX_CSC];
  float         cscRechitCluster2_match_gLLP_decay_z[N_MAX_CSC];
  float         cscRechitCluster2_match_gLLP_ctau[N_MAX_CSC];
  float         cscRechitCluster2_match_gLLP_beta[N_MAX_CSC];
  bool         cscRechitCluster2_match_gLLP_csc[N_MAX_CSC];
  float         cscRechitCluster2_match_gLLP_e[N_MAX_CSC];
  float         cscRechitCluster2_match_gLLP_pt[N_MAX_CSC];

  float         cscRechitCluster2X[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster2Y[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster2Z[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster2Time[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster2TimeTotal[N_MAX_CSC];

  float         cscRechitCluster2TimeSpread[N_MAX_CSC];
  float         cscRechitCluster2GenMuonDeltaR[N_MAX_CSC];
  float         cscRechitCluster2MajorAxis[N_MAX_CSC];
  float         cscRechitCluster2MinorAxis[N_MAX_CSC];
  float         cscRechitCluster2XSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster2YSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster2ZSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster2XYSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster2RSpread[N_MAX_CSC];   //[nCsc]


  float         cscRechitCluster2XYSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]
  float         cscRechitCluster2XSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]
  float         cscRechitCluster2YSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]
  float         cscRechitCluster2PhiSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]
  float         cscRechitCluster2EtaSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]
  float         cscRechitCluster2EtaPhiSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]
  float         cscRechitCluster2RSpread_corr[N_MAX_CSC][N_phicorr][N_rcorr];   //[nCsc]


  float         cscRechitCluster2EtaPhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster2EtaSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster2PhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster2Eta[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster2Phi[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster2Size[N_MAX_CSC];
  float         cscRechitCluster2Me11Ratio[N_MAX_CSC];
  float         cscRechitCluster2Me12Ratio[N_MAX_CSC];


  float         cscRechitCluster2MaxStationRatio[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster2MaxStation[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster2NStation[N_MAX_CSC];
  int           cscRechitCluster2NStation5[N_MAX_CSC];
  int           cscRechitCluster2NStation10[N_MAX_CSC];
  int           cscRechitCluster2NStation10perc[N_MAX_CSC];
  float          cscRechitCluster2AvgStation[N_MAX_CSC];
  float          cscRechitCluster2AvgStation5[N_MAX_CSC];
  float          cscRechitCluster2AvgStation10[N_MAX_CSC];
  float           cscRechitCluster2AvgStation10perc[N_MAX_CSC];
  float         cscRechitCluster2MaxChamberRatio[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster2MaxChamber[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster2NChamber[N_MAX_CSC];
  float         cscRechitCluster2JetVetoPt[N_MAX_CSC];
  float         cscRechitCluster2JetVetoE[N_MAX_CSC];
  float         cscRechitCluster2MuonVetoPt[N_MAX_CSC];
  float         cscRechitCluster2MuonVetoE[N_MAX_CSC];
  float         cscRechitCluster2MuonVetoPhi[N_MAX_CSC];
  float         cscRechitCluster2MuonVetoEta[N_MAX_CSC];
  float         cscRechitCluster2JetVetoPt_0p6[N_MAX_CSC];
  float         cscRechitCluster2JetVetoPt_0p8[N_MAX_CSC];
  float         cscRechitCluster2JetVetoE_0p6[N_MAX_CSC];
  float         cscRechitCluster2JetVetoE_0p8[N_MAX_CSC];
  float         cscRechitCluster2MuonVetoPt_0p6[N_MAX_CSC];
  float         cscRechitCluster2MuonVetoPt_0p8[N_MAX_CSC];
  float         cscRechitCluster2MuonVetoE_0p6[N_MAX_CSC];
  float         cscRechitCluster2MuonVetoE_0p8[N_MAX_CSC];
  bool         cscRechitCluster2ZLep1[N_MAX_CSC];
  bool         cscRechitCluster2ZLep2[N_MAX_CSC];
  int         cscRechitCluster2ZLep1Id[N_MAX_CSC];
  int         cscRechitCluster2ZLep2Id[N_MAX_CSC];
  bool         cscRechitCluster2ZLep1TightId[N_MAX_CSC];
  bool         cscRechitCluster2ZLep1LooseIso[N_MAX_CSC];
  bool         cscRechitCluster2ZLep1TightIso[N_MAX_CSC];
  bool         cscRechitCluster2ZLep1VTightIso[N_MAX_CSC];
  bool         cscRechitCluster2ZLep1VVTightIso[N_MAX_CSC];
  bool         cscRechitCluster2ZLep2TightId[N_MAX_CSC];
  bool         cscRechitCluster2ZLep2LooseIso[N_MAX_CSC];
  bool         cscRechitCluster2ZLep2TightIso[N_MAX_CSC];
  bool         cscRechitCluster2ZLep2VTightIso[N_MAX_CSC];
  bool         cscRechitCluster2ZLep2VVTightIso[N_MAX_CSC];
  bool          cscRechitCluster2MuonVetoLooseIso[N_MAX_CSC];
  bool          cscRechitCluster2MuonVetoTightIso[N_MAX_CSC];
  bool          cscRechitCluster2MuonVetoVTightIso[N_MAX_CSC];
  bool          cscRechitCluster2MuonVetoVVTightIso[N_MAX_CSC];
  bool          cscRechitCluster2MuonVetoTightId[N_MAX_CSC];


  bool          cscRechitCluster2MuonVetoIso[N_MAX_CSC];
  float         cscRechitCluster2IsoMuonVetoPt[N_MAX_CSC];
  float         cscRechitCluster2IsoMuonVetoE[N_MAX_CSC];
  float         cscRechitCluster2IsoMuonVetoPhi[N_MAX_CSC];
  float         cscRechitCluster2IsoMuonVetoEta[N_MAX_CSC];
  float         cscRechitCluster2GenMuonVetoPt[N_MAX_CSC];
  float         cscRechitCluster2GenMuonVetoE[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberPlus11[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberPlus12[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberPlus13[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberPlus21[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberPlus22[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberPlus31[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberPlus32[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberPlus41[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberPlus42[N_MAX_CSC];

  int           cscRechitCluster2NRechitChamberMinus11[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberMinus12[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberMinus13[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberMinus21[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberMinus22[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberMinus31[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberMinus32[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberMinus41[N_MAX_CSC];
  int           cscRechitCluster2NRechitChamberMinus42[N_MAX_CSC];
  float         cscRechitCluster2Met_dPhi[N_MAX_CSC];
  float         cscRechitCluster2MetXYCorr_dPhi[N_MAX_CSC];

  int           nCscRechitClusters3;

  int           cscRechitCluster3_match_Me1112_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_Me1112_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_Me1112_0p8[N_MAX_CSC];
  int           cscRechitCluster3_match_Me11_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_Me11_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_Me11_0p8[N_MAX_CSC];
  int           cscRechitCluster3_match_Me12_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_Me12_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_Me12_0p8[N_MAX_CSC];
  int           cscRechitCluster3_match_cscRechits_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_cscSeg_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_ME11Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_ME12Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_cscSeg_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_ME11Seg_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_ME12Seg_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_dtRechits_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_dtRechits_phi0p2[N_MAX_CSC];
  int           cscRechitCluster3_match_MB1_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_dtRechits_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_MB1_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_dtSeg_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_MB1Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_dtSeg_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_MB1Seg_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_RB1_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_RE12_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_RB1_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_RE12_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_highEta_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_highEta_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_highEta_0p8[N_MAX_CSC];
  float           cscRechitCluster3_match_cluster_dR[N_MAX_CSC];
  int           cscRechitCluster3_match_cluster_index[N_MAX_CSC];
  bool          cscRechitCluster3_match_gParticle[N_MAX_CSC];
  float         cscRechitCluster3_match_gParticle_minDeltaR[N_MAX_CSC];
  int           cscRechitCluster3_match_gParticle_index[N_MAX_CSC];
  int           cscRechitCluster3_match_gParticle_id[N_MAX_CSC];
  float         cscRechitCluster3_match_gParticle_eta[N_MAX_CSC];
  float         cscRechitCluster3_match_gParticle_phi[N_MAX_CSC];
  float         cscRechitCluster3_match_gParticle_E[N_MAX_CSC];
  float         cscRechitCluster3_match_gParticle_pt[N_MAX_CSC];
  int           cscRechitCluster3_match_gParticle_MotherId[N_MAX_CSC];

  bool          cscRechitCluster3_match_gLLP[N_MAX_CSC];
  int           cscRechitCluster3_match_gLLP_index[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_minDeltaR[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_eta[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_phi[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_decay_r[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_decay_x[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_decay_y[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_decay_z[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_ctau[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_beta[N_MAX_CSC];
  bool         cscRechitCluster3_match_gLLP_csc[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_e[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_pt[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_EMFracE[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_EMFracEz[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_EMFracP[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_EMFracPz[N_MAX_CSC];
  bool          cscRechitCluster3_match_gLLP_daughterKaon[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_visE[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_visEz[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_visP[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_visPz[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_lepdPhi[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_daughter0_deltaR[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_daughter1_deltaR[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_daughter2_deltaR[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_daughter3_deltaR[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_other_daughter_deltaR[N_MAX_CSC];
  int         cscRechitCluster3_match_gLLP_other_daughter_index[N_MAX_CSC];

  float         cscRechitCluster3_match_gLLP_other_eta[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_other_phi[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_other_decay_r[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_other_decay_x[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_other_decay_y[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_other_decay_z[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_other_ctau[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_other_beta[N_MAX_CSC];
  bool         cscRechitCluster3_match_gLLP_other_csc[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_other_e[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_other_pt[N_MAX_CSC];
  float         cscRechitCluster3X[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3Y[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3Z[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3Time[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3TimeTotal[N_MAX_CSC];
  float         cscRechitCluster3TimeWire[N_MAX_CSC];
  float         cscRechitCluster3TimeWirePruned[N_MAX_CSC];

  float         cscRechitCluster3TimeSpread[N_MAX_CSC];
  float         cscRechitCluster3TimeWireSpread[N_MAX_CSC];
  float         cscRechitCluster3TimeTotalSpread[N_MAX_CSC];
  float         cscRechitCluster3TimeTotalSpreadPruned[N_MAX_CSC];

  float         cscRechitCluster3GenMuonDeltaR[N_MAX_CSC];
  float         cscRechitCluster3MajorAxis[N_MAX_CSC];
  float         cscRechitCluster3MinorAxis[N_MAX_CSC];
  float         cscRechitCluster3XSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3YSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3ZSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3XYSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3RSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3DeltaRSpread[N_MAX_CSC];


  float         cscRechitCluster3XYSpread_r1p2[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3XSpread_r1p2[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3YSpread_r1p2[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3PhiSpread_r1p2[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3EtaSpread_r1p2[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3EtaPhiSpread_r1p2[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3RSpread_r1p2[N_MAX_CSC];   //[nCsc]

  float         cscRechitCluster3EtaPhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3EtaSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3PhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3Eta[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3Phi[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster3Size[N_MAX_CSC];
  float         cscRechitCluster3Me11Ratio[N_MAX_CSC];
  float         cscRechitCluster3Me12Ratio[N_MAX_CSC];


  float         cscRechitCluster3MaxStationRatio[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster3MaxStation[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster3NStation[N_MAX_CSC];
  int           cscRechitCluster3NStation5[N_MAX_CSC];
  int           cscRechitCluster3NStation10[N_MAX_CSC];
  int           cscRechitCluster3NStation10perc[N_MAX_CSC];
  float          cscRechitCluster3AvgStation[N_MAX_CSC];
  float          cscRechitCluster3AvgStation5[N_MAX_CSC];
  float          cscRechitCluster3AvgStation10[N_MAX_CSC];
  float          cscRechitCluster3AvgStation10perc[N_MAX_CSC];
  float         cscRechitCluster3MaxChamberRatio[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster3MaxChamber[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster3NChamber[N_MAX_CSC];

  float         cscRechitCluster3JetVetoElectronEnergyFraction[N_MAX_CSC];
  float         cscRechitCluster3JetVetoPhotonEnergyFraction[N_MAX_CSC];
  float         cscRechitCluster3JetVetoChargedHadronEnergyFraction[N_MAX_CSC];
  float         cscRechitCluster3JetVetoNeutralHadronEnergyFraction[N_MAX_CSC];
  float         cscRechitCluster3JetVetoMuonEnergyFraction[N_MAX_CSC];
  float         cscRechitCluster3JetVetoEta[N_MAX_CSC];
  float         cscRechitCluster3JetVetoPhi[N_MAX_CSC];
float         cscRechitCluster3JetVetoPt[N_MAX_CSC];
  float         cscRechitCluster3JetVetoE[N_MAX_CSC];
  float         cscRechitCluster3GenJetVetoPt[N_MAX_CSC];
  float         cscRechitCluster3GenJetVetoE[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoPt[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoE[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoPhi[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoEta[N_MAX_CSC];
  float         cscRechitCluster3JetVetoPt_0p6[N_MAX_CSC];
  float         cscRechitCluster3JetVetoPt_0p8[N_MAX_CSC];
  float         cscRechitCluster3JetVetoE_0p6[N_MAX_CSC];
  float         cscRechitCluster3JetVetoE_0p8[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoPt_0p6[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoPt_0p8[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoE_0p6[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoE_0p8[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1[N_MAX_CSC];
  bool         cscRechitCluster3ZLep2[N_MAX_CSC];
  int         cscRechitCluster3ZLep1Id[N_MAX_CSC];
  int         cscRechitCluster3ZLep2Id[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1TightId[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1LooseIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1TightIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1VTightIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1VVTightIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1Tag[N_MAX_CSC];

  bool         cscRechitCluster3ZLep2TightId[N_MAX_CSC];
  bool         cscRechitCluster3ZLep2LooseIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep2TightIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep2VTightIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep2VVTightIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep2Tag[N_MAX_CSC];
  bool          cscRechitCluster3MuonVetoLooseIso[N_MAX_CSC];
  bool          cscRechitCluster3MuonVetoTightIso[N_MAX_CSC];
  bool          cscRechitCluster3MuonVetoVTightIso[N_MAX_CSC];
  bool          cscRechitCluster3MuonVetoVVTightIso[N_MAX_CSC];
  bool          cscRechitCluster3MuonVetoTightId[N_MAX_CSC];
  bool          cscRechitCluster3MuonVetoLooseId[N_MAX_CSC];

  bool          cscRechitCluster3MuonVetoIso[N_MAX_CSC];
  float         cscRechitCluster3IsoMuonVetoPt[N_MAX_CSC];
  float         cscRechitCluster3IsoMuonVetoE[N_MAX_CSC];
  float         cscRechitCluster3IsoMuonVetoPhi[N_MAX_CSC];
  float         cscRechitCluster3IsoMuonVetoEta[N_MAX_CSC];
  float         cscRechitCluster3GenMuonVetoPt[N_MAX_CSC];
  float         cscRechitCluster3GenMuonVetoE[N_MAX_CSC];
  float         cscRechitCluster3GenMuonVetoProdX[N_MAX_CSC];
  float         cscRechitCluster3GenMuonVetoProdY[N_MAX_CSC];
  float         cscRechitCluster3GenMuonVetoProdZ[N_MAX_CSC];
  int         cscRechitCluster3GenMuonVetoLLPIndex[N_MAX_CSC];
  float         cscRechitCluster3GenMuonVetoLLPDist[N_MAX_CSC];
  int           cscRechitCluster3MuonVetoType[N_MAX_CSC];


  int           cscRechitCluster3NRechitChamberPlus11[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus12[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus13[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus21[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus22[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus31[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus32[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus41[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus42[N_MAX_CSC];

  int           cscRechitCluster3NRechitChamberMinus11[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus12[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus13[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus21[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus22[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus31[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus32[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus41[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus42[N_MAX_CSC];
  float         cscRechitCluster3Met_dPhi[N_MAX_CSC];
  float         cscRechitCluster3MetXYCorr_dPhi[N_MAX_CSC];

  float     cscRechitCluster3MetHEM_dPhi[N_MAX_CSC];
  float     cscRechitCluster3MetHEMXYCorr_dPhi[N_MAX_CSC];
  float     cscRechitCluster3MetEENoise_dPhi[N_MAX_CSC];
  float     cscRechitCluster3MetEENoiseXYCorr_dPhi[N_MAX_CSC];
  float     cscRechitCluster3MetJesUp_dPhi[N_MAX_CSC];
  float     cscRechitCluster3MetJesDown_dPhi[N_MAX_CSC];


  // //csc intime cluster
  // int           nCscITClusters;
  // float         cscITClusterJetVeto[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterCaloJetVeto[N_MAX_CSC];
  // float         cscITClusterMuonVeto[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterJetVetoE[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterCaloJetVetoE[N_MAX_CSC];
  // float         cscITClusterMuonVetoE[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterX[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterY[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterZ[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterRadius[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterTime[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterTimeSpread[N_MAX_CSC];
  // float         cscITClusterTimeRMS[N_MAX_CSC];
  // float         cscITClusterGenMuonDeltaR[N_MAX_CSC];
  // float         cscITClusterMajorAxis[N_MAX_CSC];
  // float         cscITClusterMinorAxis[N_MAX_CSC];
  // float         cscITClusterXSpread[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterYSpread[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterZSpread[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterEtaPhiSpread[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterEtaSpread[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterPhiSpread[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterEta[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterPhi[N_MAX_CSC];   //[nCsc]
  // int           cscITClusterSize[N_MAX_CSC];
  // float         cscITClusterMe11Ratio[N_MAX_CSC];
  // float         cscITClusterMe12Ratio[N_MAX_CSC];
  //
  // float         cscITClusterVertexR[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterVertexZ[N_MAX_CSC];   //[nCsc]
  // int           cscITClusterVertexN[N_MAX_CSC];   //[nCsc]
  // int           cscITClusterVertexN1[N_MAX_CSC];   //[nCsc]
  // int           cscITClusterVertexN5[N_MAX_CSC];   //[nCsc]
  // int           cscITClusterVertexN10[N_MAX_CSC];   //[nCsc]
  // int           cscITClusterVertexN15[N_MAX_CSC];   //[nCsc]
  // int           cscITClusterVertexN20[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterVertexChi2[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterVertexDis[N_MAX_CSC];   //[nCsc]
  // float         cscITClusterMaxStationRatio[N_MAX_CSC];   //[nCsc]
  // int           cscITClusterMaxStation[N_MAX_CSC];   //[nCsc]
  // int           cscITClusterNStation[N_MAX_CSC];
  // float         cscITClusterMaxChamberRatio[N_MAX_CSC];   //[nCsc]
  // int           cscITClusterMaxChamber[N_MAX_CSC];   //[nCsc]
  // int           cscITClusterNChamber[N_MAX_CSC];
  // int           cscITCluster_match_cscCluster_index[N_MAX_CSC];
  // float         cscITCluster_cscCluster_SizeRatio[N_MAX_CSC];
  // int           nCsc_JetVetoITCluster0p4;
  // int           nCsc_JetMuonVetoITCluster0p4;
  // int           nCsc_JetVetoITCluster0p4_Me1112Veto;
  // int           nCsc_JetMuonVetoITCluster0p4_Me1112Veto;
  //gLLP
  float gLLP_eta[2];
  float gLLP_phi[2];
  float gLLP_csc[2];
  float gLLP_beta[2];
  float gLLP_e[2];
  bool gLLP_daughterKaon[2];
  float gLLP_pt[2];
  float gLLP_lepdPhi[2];
  int gLLP_multiplicity[2];



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

  //leptons
  int nMuons;
  float muonPt[N_MAX_LEPTONS];
  float muonEta[N_MAX_LEPTONS];
  float muonPhi[N_MAX_LEPTONS];
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

  // bool lepLoosePassId[N_MAX_LEPTONS];
  // bool lepMediumPassId[N_MAX_LEPTONS];
  // bool lepTightPassId[N_MAX_LEPTONS];
  bool lepPassVetoId[N_MAX_LEPTONS];
  bool lepFromZ[N_MAX_LEPTONS];
  bool lepPassId[N_MAX_LEPTONS];
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


  float ecalNRechits[N_MAX_JETS];
  float ecalRechitE[N_MAX_JETS];
  float jetElectronEnergyFraction[N_MAX_JETS];
  float jetPhotonEnergyFraction[N_MAX_JETS];
  float jetChargedHadronEnergyFraction[N_MAX_JETS];
  float jetNeutralHadronEnergyFraction[N_MAX_JETS];

  float jetMuonEnergyFraction[N_MAX_JETS];


  bool jetPassMuFrac[N_MAX_JETS];
  float jet_match_genJet_minDeltaR[N_MAX_JETS];
  int jet_match_genJet_index[N_MAX_JETS];
  float jet_match_genJet_pt[N_MAX_JETS];

  // bool jetLoosePassId[N_MAX_JETS];
  bool jetPassId[N_MAX_JETS];
  bool jetTightPassId[N_MAX_JETS];
  bool HLTDecision[NTriggersMAX];
  bool SingleMuonTrigger;
  bool SingleEleTrigger;
  bool SingleLepTrigger;
  UInt_t wzevtNum,trig, trig_lepId, trig_lepId_dijet; //number of events that pass each criteria




  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();

};
#endif
