// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef TreeMuonSystemParticleGun_H
#define TreeMuonSystemParticleGun_H

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 200
#define N_MAX_CSCRECHITS 5000
#define N_MAX_DTRECHITS 20000
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

class TreeMuonSystemParticleGun
{

public:
  TreeMuonSystemParticleGun();
  ~TreeMuonSystemParticleGun();
  // TreeMuonSystemParticleGun::TreeMuonSystemParticleGun()
  // {
  //   InitVariables();
  // };
  // TreeMuonSystemParticleGun::~TreeMuonSystemParticleGun()
  // {
  //   if (f_) f_->Close();
  // };
  TTree *tree_;
  TFile *f_;

  UInt_t  runNum, lumiSec, evtNum, MC_condition,npv;


  //particle
  int particle1_id;
  float particle1_pt;
  float particle1_eta;
  float particle1_phi;
  float particle1_e;
  int particle2_id;
  float particle2_pt;
  float particle2_eta;
  float particle2_phi;
  float particle2_e;

  int nJets;
  float jetE[N_MAX_JETS];
  float jetPt[N_MAX_JETS];
  float jetEta[N_MAX_JETS];
  float jetPhi[N_MAX_JETS];
  bool jetTightPassId[N_MAX_JETS];

  int           nDtRechitClusters;
  float         dtRechitClusterMaxDPhi[N_MAX_CSC];
  bool          dtRechitClusterOverlap[N_MAX_CSC];
  int           dtRechitClusterMaxDPhi_index[N_MAX_CSC];
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
  int           dtRechitCluster_match_MB1hits_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_MB1hits_0p5[N_MAX_CSC];
  int           dtRechitCluster_match_MB1hits_cosmics_plus[N_MAX_CSC];
  int           dtRechitCluster_match_MB1hits_cosmics_minus[N_MAX_CSC];
  int           dtRechitCluster_match_MB1Seg_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_MB1Seg_0p5[N_MAX_CSC];
  int           dtRechitCluster_match_RPCBx_dPhi0p5[N_MAX_CSC];
  int           dtRechitCluster_match_RB1_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_RB1_dPhi0p5[N_MAX_CSC];

  float         dtRechitClusterJetVetoEta[N_MAX_CSC];
  float         dtRechitClusterJetVetoPhi[N_MAX_CSC];
  float         dtRechitClusterJetVetoPt[N_MAX_CSC];
  float         dtRechitClusterJetVetoE[N_MAX_CSC];
  

  float         dtRechitClusterMuonVetoPt[N_MAX_CSC];
  float         dtRechitClusterMuonVetoE[N_MAX_CSC];
  float         dtRechitClusterMuonVetoPhi[N_MAX_CSC];
  float         dtRechitClusterMuonVetoEta[N_MAX_CSC];
  bool          dtRechitClusterMuonVetoLooseId[N_MAX_CSC];
  bool          dtRechitClusterMuonVetoGlobal[N_MAX_CSC];

  float         dtRechitCluster_match_RPCTime_dPhi0p5[N_MAX_CSC];
  float         dtRechitCluster_match_RPCTimeSpread_dPhi0p5[N_MAX_CSC];
  int         dtRechitCluster_match_RPChits_dPhi0p5[N_MAX_CSC];
  float         dtRechitCluster_match_RPCTime_dR0p4[N_MAX_CSC];
  float         dtRechitCluster_match_RPCTimeSpread_dR0p4[N_MAX_CSC];
  int         dtRechitCluster_match_RPChits_dR0p4[N_MAX_CSC];
  float         dtRechitCluster_match_RPCTime_sameStation_dR0p4[N_MAX_CSC];
  float         dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[N_MAX_CSC];
  int         dtRechitCluster_match_RPChits_sameStation_dR0p4[N_MAX_CSC];

  float         dtRechitClusterX[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterY[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterZ[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterTime[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterTimeTotal[N_MAX_CSC];
  float         dtRechitClusterTimeWire[N_MAX_CSC];
  float         dtRechitClusterTimeWirePruned[N_MAX_CSC];
  int           dtRechitClusterWheel[N_MAX_CSC];

  float         dtRechitClusterTimeSpread[N_MAX_CSC];
  float         dtRechitClusterTimeWireSpread[N_MAX_CSC];
  float         dtRechitClusterTimeTotalSpread[N_MAX_CSC];
  float         dtRechitClusterTimeTotalSpreadPruned[N_MAX_CSC];

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
  float         cscRechitCluster_match_cluster_dR[N_MAX_CSC];
  int           cscRechitCluster_match_cluster_index[N_MAX_CSC];


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




  float         cscRechitClusterJetVetoEta[N_MAX_CSC];
  float         cscRechitClusterJetVetoPhi[N_MAX_CSC];
  float         cscRechitClusterJetVetoPt[N_MAX_CSC];
  float         cscRechitClusterJetVetoE[N_MAX_CSC];


  float         cscRechitClusterMuonVetoPt[N_MAX_CSC];
  float         cscRechitClusterMuonVetoE[N_MAX_CSC];
  float         cscRechitClusterMuonVetoPhi[N_MAX_CSC];
  float         cscRechitClusterMuonVetoEta[N_MAX_CSC];
  bool          cscRechitClusterMuonVetoLooseId[N_MAX_CSC];
  bool          cscRechitClusterMuonVetoGlobal[N_MAX_CSC];
  float           cscRechitClusterMetEENoise_dPhi[N_MAX_CSC];
float           cscRechitCluster_match_gLLP_deltaR[N_MAX_CSC];

  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();

};
#endif
