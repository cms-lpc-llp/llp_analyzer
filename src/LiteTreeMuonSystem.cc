#include "RazorHelper.h"
#include "LiteTreeMuonSystem.h"
#include "assert.h"
#include "TTree.h"
#include "DBSCAN.h"

// Constructor
LiteTreeMuonSystem::LiteTreeMuonSystem()
{
  InitVariables();
};
LiteTreeMuonSystem::~LiteTreeMuonSystem()
{
  if (f_) f_->Close();
};
void LiteTreeMuonSystem::InitVariables()
{
  runNum=0; lumiSec=0; evtNum=0; category=0;
  npv=0; npu=0;
  pileupWeight = 0; pileupWeightUp = 0; pileupWeightDown = 0;
  higgsPtWeight = 0;
  sf_facScaleUp = 0; sf_facScaleDown = 0; sf_renScaleUp = 0; sf_renScaleDown = 0; sf_facRenScaleUp = 0; sf_facRenScaleDown = 0;

  for (int i = 0; i < 9; i++)
  {
    higgsPtWeightSys[i] = 0;
  }
  weight=-1.0;//rho=-1;
  met=-1; metPhi=-1; metXYCorr=-1; metPhiXYCorr=-1; HT = 0.0; jetMet_dPhi = -999.;jetMet_dPhiMin = 999.;jetMet_dPhiMin4 = 999.;
  Flag_HBHENoiseFilter = false; Flag_HBHEIsoNoiseFilter = false; Flag_BadPFMuonFilter = false; Flag_CSCTightHaloFilter = false; Flag_goodVertices = false;
  Flag_ecalBadCalibFilter = false; Flag_all = false; Flag_globalSuperTightHalo2016Filter = false; Flag_BadChargedCandidateFilter = false; Flag_eeBadScFilter = false;

  Flag2_HBHENoiseFilter = false; Flag2_HBHEIsoNoiseFilter = false; Flag2_BadPFMuonFilter = false; Flag2_globalSuperTightHalo2016Filter = false;
  Flag2_globalTightHalo2016Filter = false; Flag2_BadChargedCandidateFilter = false; Flag2_EcalDeadCellTriggerPrimitiveFilter = false; Flag2_ecalBadCalibFilter = false;
  Flag2_eeBadScFilter = false;
  Flag2_all = false;
  mH = 0; mX = 0; ctau = 0;



  metJESUp = -999.;metJESDown = -999.;
  gLepId = 0;
  gLepPt = 0.; gLepPhi = 0.; gLepEta = 0.; gLepE = 0.;
  gHiggsPt = 0.; gHiggsPhi = 0.; gHiggsEta = 0.; gHiggsE = 0.;
  //CSC
  // nCscSegClusters = 0;
  nCscRechitClusters = 0;
  nCscRechitClusters2 = 0;
  nCscRechitClusters3 = 0;
  nCscClusters = 0;
  nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto = 0;
  nCscRechits = 0;
  nEarlyCscRechits = 0;
  nLateCscRechits = 0;
  nEarly2CscRechits = 0;
  nLate2CscRechits = 0;
  nCscPositiveYRechits = 0;
  nCscNegativeYRechits = 0;
  cscPosTpeak = 0.0;
  cscNegTpeak = 0.0;
  nDtRings = 0;
  nCscRings = 0;


  nCscRechitsChamberPlus11 = 0;
  nCscRechitsChamberPlus12 = 0;
  nCscRechitsChamberPlus13 = 0;
  nCscRechitsChamberPlus21 = 0;
  nCscRechitsChamberPlus22 = 0;
  nCscRechitsChamberPlus31 = 0;
  nCscRechitsChamberPlus32 = 0;
  nCscRechitsChamberPlus41 = 0;
  nCscRechitsChamberPlus42 = 0;
  nCscRechitsChamberMinus11 = 0;
  nCscRechitsChamberMinus12 = 0;
  nCscRechitsChamberMinus13 = 0;
  nCscRechitsChamberMinus21 = 0;
  nCscRechitsChamberMinus22 = 0;
  nCscRechitsChamberMinus31 = 0;
  nCscRechitsChamberMinus32 = 0;
  nCscRechitsChamberMinus41 = 0;
  nCscRechitsChamberMinus42 = 0;  // nCscClusters[i] = 0;



  // for (int i = 0; i<36; i++)
  // {
  //   nCscRechitsChamberPlus11[i] = 0;
  //   nCscRechitsChamberPlus12[i] = 0;
  //   nCscRechitsChamberPlus13[i] = 0;
  //   nCscRechitsChamberPlus21[i] = 0;
  //   nCscRechitsChamberPlus22[i] = 0;
  //   nCscRechitsChamberPlus31[i] = 0;
  //   nCscRechitsChamberPlus32[i] = 0;
  //   nCscRechitsChamberPlus41[i] = 0;
  //   nCscRechitsChamberPlus42[i] = 0;
  //   nCscRechitsChamberMinus11[i] = 0;
  //   nCscRechitsChamberMinus12[i] = 0;
  //   nCscRechitsChamberMinus13[i] = 0;
  //   nCscRechitsChamberMinus21[i] = 0;
  //   nCscRechitsChamberMinus22[i] = 0;
  //   nCscRechitsChamberMinus31[i] = 0;
  //   nCscRechitsChamberMinus32[i] = 0;
  //   nCscRechitsChamberMinus41[i] = 0;
  //   nCscRechitsChamberMinus42[i] = 0;  // nCscClusters[i] = 0;
  //
  // }

  nDTRechits = 0;
  nDTNegativeYRechits  = 0;
  nDTPositiveYRechits = 0;
  nDTRechitsChamberMinus12 = 0;
  nDTRechitsChamberMinus11 = 0;
  nDTRechitsChamber10 = 0;
  nDTRechitsChamberPlus11 = 0;
  nDTRechitsChamberPlus12 = 0;
  nDTRechitsChamberMinus22 = 0;
  nDTRechitsChamberMinus21 = 0;
  nDTRechitsChamber20 = 0;
  nDTRechitsChamberPlus21 = 0;
  nDTRechitsChamberPlus22 = 0;
  nDTRechitsChamberMinus32 = 0;
  nDTRechitsChamberMinus31 = 0;
  nDTRechitsChamber30 = 0;
  nDTRechitsChamberPlus31 = 0;
  nDTRechitsChamberPlus32 = 0;
  nDTRechitsChamberMinus42 = 0;
  nDTRechitsChamberMinus41 = 0;
  nDTRechitsChamber40 = 0;
  nDTRechitsChamberPlus41 = 0;
  nDTRechitsChamberPlus42 = 0;
  // nCscITClusters = 0;
  // nCsc_JetVetoITCluster0p4 = 0;
  // nCsc_JetMuonVetoITCluster0p4 = 0;
  // nCsc_JetVetoITCluster0p4_Me1112Veto = 0;
  // nCsc_JetMuonVetoITCluster0p4_Me1112Veto = 0;
  // nCsc_JetVetoCluster0p4 = 0;
  // nCsc_JetMuonVetoCluster0p4 = 0;
  // nCsc_JetVetoCluster0p4_Me1112Veto = 0;
  // nCsc_JetMuonVetoCluster0p4_Me1112Veto = 0;

  for( int i = 0; i < N_MAX_CSCRECHITS; i++ )
  {
    cscRechitsStation[i] = -999;
    cscRechitsChamber[i] = -999;

    cscRechitsPhi[i] = -999;   //[nCsc]
    cscRechitsEta[i] = -999;   //[nCsc]
    cscRechitsQuality[i] = 999;
    cscRechitsX[i] = -999;   //[nCsc]
    cscRechitsY[i] = -999;   //[nCsc]
    cscRechitsZ[i] = -999;   //[nCsc]
    cscRechitsTpeak[i] = -999;   //[nCsc]
    cscRechitsTwire[i] = -999;   //[nCsc]
    cscRechitsClusterId[i] = -999;   //[nCsc]
    cscRechitsCluster2Id[i] = -999;   //[nCsc]

  }
  for( int i = 0; i < N_MAX_CSC; i++ )
  {





    cscCluster_match_gLLP[i] = false;
    cscCluster_match_gLLP_minDeltaR[i] = 999;
    cscCluster_match_gLLP_index[i] = 999;
    cscClusterSize[i] = -999;
    cscClusterX[i] = -999.;
    cscClusterY[i] = -999.;
    cscClusterZ[i] = -999.;
    cscClusterTime[i] = -999.;
    cscClusterGenMuonDeltaR[i] = 999.;
    cscClusterTimeSpread[i] = -999.;
    cscClusterMajorAxis[i] = -999.;
    cscClusterMinorAxis[i] = -999.;
    cscClusterXSpread[i] = -999.;
    cscClusterYSpread[i] = -999.;
    cscClusterZSpread[i] = -999.;

    cscClusterEtaPhiSpread[i] = -999.;
    cscClusterEtaSpread[i] = -999.;
    cscClusterPhiSpread[i] = -999.;
    cscClusterEta[i] = -999.;
    cscClusterPhi[i] = -999.;
    cscClusterJetVetoPt[i] = 0.0;
    // cscClusterCaloJetVeto[i] = 0.0;
    cscClusterMuonVetoPt[i] = 0.0;
    cscClusterJetVetoE[i] = 0.0;
    // cscClusterCaloJetVetoE[i] = 0.0;
    cscClusterMuonVetoE[i] = 0.0;
    cscClusterNChamber[i] = -999;
    cscClusterMaxChamberRatio[i] = -999.;
    cscClusterMaxChamber[i] = -999;
    cscClusterNStation[i] = -999;
    cscClusterMaxStationRatio[i] = -999.;
    cscClusterMaxStation[i] = -999;
    cscClusterMe11Ratio[i] = -999.;
    cscClusterMe12Ratio[i] = -999.;
    cscClusterNSegmentChamberPlus11[i] = -999;
    cscClusterNSegmentChamberPlus12[i] = -999;
    cscClusterNSegmentChamberPlus13[i] = -999;
    cscClusterNSegmentChamberPlus21[i] = -999;
    cscClusterNSegmentChamberPlus22[i] = -999;
    cscClusterNSegmentChamberPlus31[i] = -999;
    cscClusterNSegmentChamberPlus32[i] = -999;
    cscClusterNSegmentChamberPlus41[i] = -999;
    cscClusterNSegmentChamberPlus42[i] = -999;
    cscClusterNSegmentChamberMinus11[i] = -999;
    cscClusterNSegmentChamberMinus12[i] = -999;
    cscClusterNSegmentChamberMinus13[i] = -999;
    cscClusterNSegmentChamberMinus21[i] = -999;
    cscClusterNSegmentChamberMinus22[i] = -999;
    cscClusterNSegmentChamberMinus31[i] = -999;
    cscClusterNSegmentChamberMinus32[i] = -999;
    cscClusterNSegmentChamberMinus41[i] = -999;
    cscClusterNSegmentChamberMinus42[i] = -999;
    // cscClusterVertexR[i] = 0.0;
    // cscClusterVertexZ[i] = 0.0;
    // cscClusterVertexDis[i] = 0.0;
    // cscClusterVertexChi2[i] = 0.0;
    // cscClusterVertexN[i] = 0;
    // cscClusterVertexN1[i] = 0;
    // cscClusterVertexN5[i] = 0;
    // cscClusterVertexN10[i] = 0;
    // cscClusterVertexN15[i] = 0;
    // cscClusterVertexN20[i] = 0;
    /*cscSegClusterSize[i] = -999;
    cscSegClusterX[i] = -999.;
    cscSegClusterY[i] = -999.;
    cscSegClusterZ[i] = -999.;
    cscSegClusterTime[i] = -999.;
    cscSegClusterGenMuonDeltaR[i] = 999.;
    cscSegClusterTimeSpread[i] = -999.;
    cscSegClusterMajorAxis[i] = -999.;
    cscSegClusterMinorAxis[i] = -999.;
    cscSegClusterXSpread[i] = -999.;
    cscSegClusterYSpread[i] = -999.;
    cscSegClusterZSpread[i] = -999.;
    cscSegClusterEtaPhiSpread[i] = -999.;
    cscSegClusterEtaSpread[i] = -999.;
    cscSegClusterPhiSpread[i] = -999.;
    cscSegClusterEta[i] = -999.;
    cscSegClusterPhi[i] = -999.;
    cscSegClusterJetVetoPt[i] = 0.0;
    // cscSegClusterCaloJetVeto[i] = 0.0;
    cscSegClusterMuonVetoPt[i] = 0.0;
    cscSegClusterJetVetoE[i] = 0.0;
    // cscSegClusterCaloJetVetoE[i] = 0.0;
    cscSegClusterMuonVetoE[i] = 0.0;
    cscSegClusterNChamber[i] = -999;
    cscSegClusterMaxChamberRatio[i] = -999.;
    cscSegClusterMaxChamber[i] = -999;
    cscSegClusterNStation[i] = -999;
    cscSegClusterMaxStationRatio[i] = -999.;
    cscSegClusterMaxStation[i] = -999;
    cscSegClusterMe11Ratio[i] = -999.;
    cscSegClusterMe12Ratio[i] = -999.;
    cscSegClusterNSegmentChamberPlus11[i] = -999;
    cscSegClusterNSegmentChamberPlus12[i] = -999;
    cscSegClusterNSegmentChamberPlus13[i] = -999;
    cscSegClusterNSegmentChamberPlus21[i] = -999;
    cscSegClusterNSegmentChamberPlus22[i] = -999;
    cscSegClusterNSegmentChamberPlus31[i] = -999;
    cscSegClusterNSegmentChamberPlus32[i] = -999;
    cscSegClusterNSegmentChamberPlus41[i] = -999;
    cscSegClusterNSegmentChamberPlus42[i] = -999;
    cscSegClusterNSegmentChamberMinus11[i] = -999;
    cscSegClusterNSegmentChamberMinus12[i] = -999;
    cscSegClusterNSegmentChamberMinus13[i] = -999;
    cscSegClusterNSegmentChamberMinus21[i] = -999;
    cscSegClusterNSegmentChamberMinus22[i] = -999;
    cscSegClusterNSegmentChamberMinus31[i] = -999;
    cscSegClusterNSegmentChamberMinus32[i] = -999;
    cscSegClusterNSegmentChamberMinus41[i] = -999;
    cscSegClusterNSegmentChamberMinus42[i] = -999;
    cscSegCluster_match_gParticle_id[i] = -999;
    cscSegCluster_match_gParticle_index[i] = -999;
    cscSegCluster_match_gParticle_minDeltaR[i] = -999.;
    cscSegClusterMet_dPhi[i] = -999.;*/


    cscRechitCluster2_match_gParticle[i] = false;
    cscRechitCluster2_match_gParticle_minDeltaR[i] = -999.;
    cscRechitCluster2_match_gParticle_index[i] = -999;
    cscRechitCluster2_match_gParticle_id[i] = -999;
    cscRechitCluster2_match_gParticle_eta[i] = -999.;
    cscRechitCluster2_match_gParticle_phi[i] = -999.;
    cscRechitCluster2_match_gParticle_E[i] = -999.;
    cscRechitCluster2_match_gParticle_pt[i] = -999.;
    cscRechitCluster2_match_gParticle_MotherId[i]  = -999;
    cscRechitCluster_match_gLLP[i] = false;
    cscRechitCluster_match_gLLP_minDeltaR[i] = 999;
    cscRechitCluster_match_gLLP_index[i] = 999;
    cscRechitCluster_match_gLLP_eta[i] = 999.;
    cscRechitCluster_match_gLLP_phi[i] = 999.;
    cscRechitCluster_match_gLLP_decay_r[i] = 999.;
    cscRechitCluster_match_gLLP_decay_x[i] = 999.;
    cscRechitCluster_match_gLLP_decay_y[i] = 999.;
    cscRechitCluster_match_gLLP_decay_z[i] = 999.;
    cscRechitCluster_match_gLLP_ctau[i] = 999.;
    cscRechitCluster_match_gLLP_beta[i] = 999.;
    cscRechitCluster_match_gLLP_csc[i] = false;



    cscRechitClusterSize[i] = -999;
    cscRechitClusterX[i] = -999.;
    cscRechitClusterY[i] = -999.;
    cscRechitClusterZ[i] = -999.;
    cscRechitClusterTime[i] = -999.;
    cscRechitClusterTimeTotal[i] = -999.;
    cscRechitClusterGenMuonDeltaR[i] = 999.;
    cscRechitClusterTimeSpread[i] = -999.;
    cscRechitClusterMajorAxis[i] = -999.;
    cscRechitClusterMinorAxis[i] = -999.;
    cscRechitClusterXYSpread[i] = -999.;

    cscRechitClusterXSpread[i] = -999.;
    cscRechitClusterYSpread[i] = -999.;
    cscRechitClusterZSpread[i] = -999.;

    cscRechitClusterEtaPhiSpread[i] = -999.;
    cscRechitClusterEtaSpread[i] = -999.;
    cscRechitClusterPhiSpread[i] = -999.;
    cscRechitClusterEta[i] = -999.;
    cscRechitClusterPhi[i] = -999.;
    cscRechitClusterJetVetoPt[i] = 0.0;
    // cscRechitClusterCaloJetVeto[i] = 0.0;
    cscRechitClusterMuonVetoPt[i] = 0.0;
    cscRechitClusterJetVetoE[i] = 0.0;
    // cscRechitClusterCaloJetVetoE[i] = 0.0;
    cscRechitClusterMuonVetoE[i] = 0.0;



    cscRechitClusterXYSpread_corr[i] = -999.;
    cscRechitClusterXSpread_corr[i] = -999.;
    cscRechitClusterYSpread_corr[i] = -999.;
    cscRechitClusterPhiSpread_corr[i] = -999.;
    cscRechitClusterEtaSpread_corr[i] = -999.;
    cscRechitClusterEtaPhiSpread_corr[i] = -999.;
    cscRechitClusterRSpread_corr[i] = -999.;

    cscRechitClusterGenMuonVetoPt[i] = 0.0;
    cscRechitClusterGenMuonVetoE[i] = 0.0;
    cscRechitCluster_match_MB1Seg_0p4[i] = 0;
    cscRechitCluster_match_RB1_0p4[i] = 0;
    cscRechitCluster_match_RE12_0p4[i] = 0;

    cscRechitClusterNChamber[i] = -999;
    cscRechitClusterMaxChamberRatio[i] = -999.;
    cscRechitClusterMaxChamber[i] = -999;
    cscRechitClusterNStation[i] = -999;
    cscRechitClusterNStation5[i] = -999;
    cscRechitClusterNStation10perc[i] = -999;
    cscRechitClusterAvgStation[i] = -999.;
    cscRechitClusterAvgStation5[i] = -999.;
    cscRechitClusterAvgStation10perc[i] = -999.;

    cscRechitClusterMaxStationRatio[i] = -999.;
    cscRechitClusterMaxStation[i] = -999;
    cscRechitClusterMe11Ratio[i] = -999.;
    cscRechitClusterMe12Ratio[i] = -999.;
    cscRechitClusterNRechitChamberPlus11[i] = -999;
    cscRechitClusterNRechitChamberPlus12[i] = -999;
    cscRechitClusterNRechitChamberPlus13[i] = -999;
    cscRechitClusterNRechitChamberPlus21[i] = -999;
    cscRechitClusterNRechitChamberPlus22[i] = -999;
    cscRechitClusterNRechitChamberPlus31[i] = -999;
    cscRechitClusterNRechitChamberPlus32[i] = -999;
    cscRechitClusterNRechitChamberPlus41[i] = -999;
    cscRechitClusterNRechitChamberPlus42[i] = -999;
    cscRechitClusterNRechitChamberMinus11[i] = -999;
    cscRechitClusterNRechitChamberMinus12[i] = -999;
    cscRechitClusterNRechitChamberMinus13[i] = -999;
    cscRechitClusterNRechitChamberMinus21[i] = -999;
    cscRechitClusterNRechitChamberMinus22[i] = -999;
    cscRechitClusterNRechitChamberMinus31[i] = -999;
    cscRechitClusterNRechitChamberMinus32[i] = -999;
    cscRechitClusterNRechitChamberMinus41[i] = -999;
    cscRechitClusterNRechitChamberMinus42[i] = -999;
    cscRechitClusterMet_dPhi[i] = 999.;

    cscRechitCluster2_match_Me1112_0p4[i] = 0;
    cscRechitCluster2_match_Me1112_0p6[i] = 0;
    cscRechitCluster2_match_Me1112_0p8[i] = 0;
    cscRechitCluster2_match_Me11_0p4[i] = 0;
    cscRechitCluster2_match_Me11_0p6[i] = 0;
    cscRechitCluster2_match_Me11_0p8[i] = 0;
    cscRechitCluster2_match_Me12_0p4[i] = 0;
    cscRechitCluster2_match_Me12_0p6[i] = 0;
    cscRechitCluster2_match_Me12_0p8[i] = 0;

    cscRechitCluster2_match_cscRechits_0p4[i] = 0;

    cscRechitCluster2_match_cscSeg_0p4[i] = 0;
    cscRechitCluster2_match_ME11Seg_0p4[i] = 0;
    cscRechitCluster2_match_ME12Seg_0p4[i] = 0;
    cscRechitCluster2_match_cscSeg_0p6[i] = 0;
    cscRechitCluster2_match_ME11Seg_0p6[i] = 0;
    cscRechitCluster2_match_ME12Seg_0p6[i] = 0;

    cscRechitCluster2_match_dtRechits_0p4[i] = 0;
    cscRechitCluster2_match_dtRechits_0p6[i] = 0;
    cscRechitCluster2_match_MB1_0p4[i] = 0;
    cscRechitCluster2_match_MB1_0p6[i] = 0;
    cscRechitCluster2_match_dtSeg_0p4[i] = 0;
    cscRechitCluster2_match_dtSeg_0p6[i] = 0;
    cscRechitCluster2_match_MB1Seg_0p4[i] = 0;
    cscRechitCluster2_match_MB1Seg_0p6[i] = 0;
    cscRechitCluster2_match_RB1_0p4[i] = 0;
    cscRechitCluster2_match_RE12_0p4[i] = 0;
    cscRechitCluster2_match_RB1_0p6[i] = 0;
    cscRechitCluster2_match_RE12_0p6[i] = 0;
    cscRechitCluster2_match_highEta_0p4[i] = 0;
    cscRechitCluster2_match_highEta_0p6[i] = 0;
    cscRechitCluster2_match_highEta_0p8[i] = 0;
    cscRechitCluster2_match_cluster_dR[i] = 999.;
    cscRechitCluster2_match_cluster_index[i] = 999;

      cscRechitCluster2_match_gLLP[i] = false;
      cscRechitCluster2_match_gLLP_minDeltaR[i] = 999;
      cscRechitCluster2_match_gLLP_index[i] = 999;
      cscRechitCluster2_match_gLLP_eta[i] = 999.;
      cscRechitCluster2_match_gLLP_phi[i] = 999.;
      cscRechitCluster2_match_gLLP_decay_r[i] = 999.;
      cscRechitCluster2_match_gLLP_decay_x[i] = 999.;
      cscRechitCluster2_match_gLLP_decay_y[i] = 999.;
      cscRechitCluster2_match_gLLP_decay_z[i] = 999.;
      cscRechitCluster2_match_gLLP_ctau[i] = 999.;
      cscRechitCluster2_match_gLLP_beta[i] = 999.;
      cscRechitCluster2_match_gLLP_csc[i] = false;



      cscRechitCluster2Size[i] = -999;
      cscRechitCluster2X[i] = -999.;
      cscRechitCluster2Y[i] = -999.;
      cscRechitCluster2Z[i] = -999.;
      cscRechitCluster2Time[i] = -999.;
      cscRechitCluster2TimeTotal[i] = -999.;
      cscRechitCluster2GenMuonDeltaR[i] = 999.;
      cscRechitCluster2TimeSpread[i] = -999.;
      cscRechitCluster2MajorAxis[i] = -999.;
      cscRechitCluster2MinorAxis[i] = -999.;
      cscRechitCluster2XSpread[i] = -999.;
      cscRechitCluster2YSpread[i] = -999.;
      cscRechitCluster2ZSpread[i] = -999.;
      cscRechitCluster2EtaPhiSpread[i] = -999.;
      cscRechitCluster2XYSpread[i] = -999.;
      cscRechitCluster2RSpread[i] = -999.;

      for (int j = 0; j < 20; j++)
      {
        for(int k = 0; k < 20; k++)
        {
          cscRechitCluster2XYSpread_corr[i][j][k] = -999.;
          cscRechitCluster2XSpread_corr[i][j][k] = -999.;
          cscRechitCluster2YSpread_corr[i][j][k] = -999.;
          cscRechitCluster2PhiSpread_corr[i][j][k] = -999.;
          cscRechitCluster2EtaSpread_corr[i][j][k] = -999.;
          cscRechitCluster2EtaPhiSpread_corr[i][j][k] = -999.;
          cscRechitCluster2RSpread_corr[i][j][k] = -999.;
        }
      }

      cscRechitCluster2EtaSpread[i] = -999.;
      cscRechitCluster2PhiSpread[i] = -999.;
      cscRechitCluster2Eta[i] = -999.;
      cscRechitCluster2Phi[i] = -999.;
      cscRechitCluster2JetVetoPt[i] = 0.0;
      cscRechitCluster2JetVetoE[i] = 0.0;
      cscRechitCluster2MuonVetoPt[i] = 0.0;
      cscRechitCluster2MuonVetoE[i] = 0.0;
      cscRechitCluster2MuonVetoPhi[i] = 0.0;
      cscRechitCluster2MuonVetoEta[i] = 0.0;
      cscRechitCluster2MuonVetoIso[i] = false;
      cscRechitCluster2MuonVetoLooseIso[i] = false;
      cscRechitCluster2MuonVetoTightIso[i] = false;
      cscRechitCluster2MuonVetoVTightIso[i] = false;
      cscRechitCluster2MuonVetoVVTightIso[i] = false;
      cscRechitCluster2MuonVetoTightId[i] = false;

      cscRechitCluster2IsoMuonVetoPt[i] = 0.0;
      cscRechitCluster2IsoMuonVetoE[i] = 0.0;
      cscRechitCluster2IsoMuonVetoPhi[i] = 0.0;
      cscRechitCluster2IsoMuonVetoEta[i] = 0.0;
      cscRechitCluster2GenMuonVetoE[i] = 0.0;
      cscRechitCluster2GenMuonVetoPt[i] = 0.0;
      cscRechitCluster2JetVetoPt_0p6[i] = 0.0;
      cscRechitCluster2JetVetoPt_0p8[i] = 0.0;
      cscRechitCluster2JetVetoE_0p6[i] = 0.0;
      cscRechitCluster2JetVetoE_0p8[i] = 0.0;
      cscRechitCluster2MuonVetoPt_0p6[i] = 0.0;
      cscRechitCluster2MuonVetoPt_0p8[i] = 0.0;
      cscRechitCluster2MuonVetoE_0p6[i] = 0.0;
      cscRechitCluster2MuonVetoE_0p8[i] = 0.0;
      cscRechitCluster2ZLep1[i] = false;
      cscRechitCluster2ZLep2[i] = false;
      cscRechitCluster2ZLep1Id[i] = -999;
      cscRechitCluster2ZLep2Id[i] = -999;

      cscRechitCluster2ZLep1LooseIso[i] = false;
      cscRechitCluster2ZLep1TightIso[i] = false;
      cscRechitCluster2ZLep1VTightIso[i] = false;
      cscRechitCluster2ZLep1VVTightIso[i] = false;
      cscRechitCluster2ZLep1TightId[i] = false;
      cscRechitCluster2ZLep2LooseIso[i] = false;
      cscRechitCluster2ZLep2TightIso[i] = false;
      cscRechitCluster2ZLep2VTightIso[i] = false;
      cscRechitCluster2ZLep2VVTightIso[i] = false;
      cscRechitCluster2ZLep2TightId[i] = false;

      cscRechitCluster2NChamber[i] = -999;
      cscRechitCluster2MaxChamberRatio[i] = -999.;
      cscRechitCluster2MaxChamber[i] = -999;
      cscRechitCluster2NStation[i] = -999;
      cscRechitCluster2NStation5[i] = -999;
      cscRechitCluster2NStation10perc[i] = -999;
      cscRechitCluster2AvgStation[i] = -999.;
      cscRechitCluster2AvgStation5[i] = -999.;
      cscRechitCluster2AvgStation10perc[i] = -999.;
      cscRechitCluster2MaxStationRatio[i] = -999.;
      cscRechitCluster2MaxStation[i] = -999;
      cscRechitCluster2Me11Ratio[i] = -999.;
      cscRechitCluster2Me12Ratio[i] = -999.;
      cscRechitCluster2NRechitChamberPlus11[i] = -999;
      cscRechitCluster2NRechitChamberPlus12[i] = -999;
      cscRechitCluster2NRechitChamberPlus13[i] = -999;
      cscRechitCluster2NRechitChamberPlus21[i] = -999;
      cscRechitCluster2NRechitChamberPlus22[i] = -999;
      cscRechitCluster2NRechitChamberPlus31[i] = -999;
      cscRechitCluster2NRechitChamberPlus32[i] = -999;
      cscRechitCluster2NRechitChamberPlus41[i] = -999;
      cscRechitCluster2NRechitChamberPlus42[i] = -999;
      cscRechitCluster2NRechitChamberMinus11[i] = -999;
      cscRechitCluster2NRechitChamberMinus12[i] = -999;
      cscRechitCluster2NRechitChamberMinus13[i] = -999;
      cscRechitCluster2NRechitChamberMinus21[i] = -999;
      cscRechitCluster2NRechitChamberMinus22[i] = -999;
      cscRechitCluster2NRechitChamberMinus31[i] = -999;
      cscRechitCluster2NRechitChamberMinus32[i] = -999;
      cscRechitCluster2NRechitChamberMinus41[i] = -999;
      cscRechitCluster2NRechitChamberMinus42[i] = -999;
      cscRechitCluster2Met_dPhi[i] = 999.;
      cscRechitCluster2MetXYCorr_dPhi[i] = 999.;

      cscRechitCluster3_match_Me1112_0p4[i] = 0;
      cscRechitCluster3_match_Me1112_0p6[i] = 0;
      cscRechitCluster3_match_Me1112_0p8[i] = 0;
      cscRechitCluster3_match_Me11_0p4[i] = 0;
      cscRechitCluster3_match_Me11_0p6[i] = 0;
      cscRechitCluster3_match_Me11_0p8[i] = 0;
      cscRechitCluster3_match_Me12_0p4[i] = 0;
      cscRechitCluster3_match_Me12_0p6[i] = 0;
      cscRechitCluster3_match_Me12_0p8[i] = 0;

      cscRechitCluster3_match_cscRechits_0p4[i] = 0;

      cscRechitCluster3_match_cscSeg_0p4[i] = 0;
      cscRechitCluster3_match_ME11Seg_0p4[i] = 0;
      cscRechitCluster3_match_ME12Seg_0p4[i] = 0;
      cscRechitCluster3_match_cscSeg_0p6[i] = 0;
      cscRechitCluster3_match_ME11Seg_0p6[i] = 0;
      cscRechitCluster3_match_ME12Seg_0p6[i] = 0;

      cscRechitCluster3_match_dtRechits_0p4[i] = 0;
      cscRechitCluster3_match_dtRechits_0p6[i] = 0;
      cscRechitCluster3_match_MB1_0p4[i] = 0;
      cscRechitCluster3_match_MB1_0p6[i] = 0;
      cscRechitCluster3_match_dtSeg_0p4[i] = 0;
      cscRechitCluster3_match_dtSeg_0p6[i] = 0;
      cscRechitCluster3_match_MB1Seg_0p4[i] = 0;
      cscRechitCluster3_match_MB1Seg_0p6[i] = 0;
      cscRechitCluster3_match_RB1_0p4[i] = 0;
      cscRechitCluster3_match_RE12_0p4[i] = 0;
      cscRechitCluster3_match_RB1_0p6[i] = 0;
      cscRechitCluster3_match_RE12_0p6[i] = 0;
      cscRechitCluster3_match_highEta_0p4[i] = 0;
      cscRechitCluster3_match_highEta_0p6[i] = 0;
      cscRechitCluster3_match_highEta_0p8[i] = 0;
      cscRechitCluster3_match_cluster_dR[i] = 999.;
      cscRechitCluster3_match_cluster_index[i] = 999;


      cscRechitCluster3_match_gParticle[i] = false;
      cscRechitCluster3_match_gParticle_minDeltaR[i] = -999.;
      cscRechitCluster3_match_gParticle_index[i] = -999;
      cscRechitCluster3_match_gParticle_id[i] = -999;
      cscRechitCluster3_match_gParticle_eta[i] = -999.;
      cscRechitCluster3_match_gParticle_phi[i] = -999.;
      cscRechitCluster3_match_gParticle_E[i] = -999.;
      cscRechitCluster3_match_gParticle_pt[i] = -999.;
      cscRechitCluster3_match_gParticle_MotherId[i]  = -999;

        cscRechitCluster3_match_gLLP[i] = false;
        cscRechitCluster3_match_gLLP_minDeltaR[i] = 999;
        cscRechitCluster3_match_gLLP_index[i] = 999;
        cscRechitCluster3_match_gLLP_eta[i] = 999.;
        cscRechitCluster3_match_gLLP_phi[i] = 999.;
        cscRechitCluster3_match_gLLP_decay_r[i] = 999.;
        cscRechitCluster3_match_gLLP_decay_x[i] = 999.;
        cscRechitCluster3_match_gLLP_decay_y[i] = 999.;
        cscRechitCluster3_match_gLLP_decay_z[i] = 999.;
        cscRechitCluster3_match_gLLP_ctau[i] = 999.;
        cscRechitCluster3_match_gLLP_beta[i] = 999.;
        cscRechitCluster3_match_gLLP_csc[i] = false;



        cscRechitCluster3Size[i] = -999;
        cscRechitCluster3X[i] = -999.;
        cscRechitCluster3Y[i] = -999.;
        cscRechitCluster3Z[i] = -999.;
        cscRechitCluster3Time[i] = -999.;
        cscRechitCluster3TimeTotal[i] = -999.;
        cscRechitCluster3GenMuonDeltaR[i] = 999.;
        cscRechitCluster3TimeSpread[i] = -999.;
        cscRechitCluster3MajorAxis[i] = -999.;
        cscRechitCluster3MinorAxis[i] = -999.;
        cscRechitCluster3XSpread[i] = -999.;
        cscRechitCluster3YSpread[i] = -999.;
        cscRechitCluster3ZSpread[i] = -999.;
        cscRechitCluster3EtaPhiSpread[i] = -999.;
        cscRechitCluster3XYSpread[i] = -999.;
        cscRechitCluster3RSpread[i] = -999.;

        cscRechitCluster3XYSpread_corr2[i] = -999.;
        cscRechitCluster3XSpread_corr2[i] = -999.;
        cscRechitCluster3YSpread_corr2[i] = -999.;
        cscRechitCluster3PhiSpread_corr2[i] = -999.;
        cscRechitCluster3EtaSpread_corr2[i] = -999.;
        cscRechitCluster3EtaPhiSpread_corr2[i] = -999.;
        cscRechitCluster3RSpread_corr2[i] = -999.;
        for (int j = 0; j < N_phicorr; j++)
        {
          for(int k = 0; k < N_rcorr; k++)
          {
            cscRechitCluster3XYSpread_corr[i][j][k] = -999.;
            cscRechitCluster3XSpread_corr[i][j][k] = -999.;
            cscRechitCluster3YSpread_corr[i][j][k] = -999.;
            cscRechitCluster3PhiSpread_corr[i][j][k] = -999.;
            cscRechitCluster3EtaSpread_corr[i][j][k] = -999.;
            cscRechitCluster3EtaPhiSpread_corr[i][j][k] = -999.;
            cscRechitCluster3RSpread_corr[i][j][k] = -999.;
          }
        }

        cscRechitCluster3EtaSpread[i] = -999.;
        cscRechitCluster3PhiSpread[i] = -999.;
        cscRechitCluster3Eta[i] = -999.;
        cscRechitCluster3Phi[i] = -999.;
        cscRechitCluster3JetVetoPt[i] = 0.0;
        cscRechitCluster3JetVetoE[i] = 0.0;
        cscRechitCluster3MuonVetoPt[i] = 0.0;
        cscRechitCluster3MuonVetoE[i] = 0.0;
        cscRechitCluster3MuonVetoPhi[i] = 0.0;
        cscRechitCluster3MuonVetoEta[i] = 0.0;
        cscRechitCluster3MuonVetoLooseIso[i] = false;
        cscRechitCluster3MuonVetoTightIso[i] = false;
        cscRechitCluster3MuonVetoVTightIso[i] = false;
        cscRechitCluster3MuonVetoVVTightIso[i] = false;
        cscRechitCluster3MuonVetoTightId[i] = false;
        cscRechitCluster3MuonVetoIso[i] = false;
        cscRechitCluster3IsoMuonVetoPt[i] = 0.0;
        cscRechitCluster3IsoMuonVetoE[i] = 0.0;
        cscRechitCluster3IsoMuonVetoPhi[i] = 0.0;
        cscRechitCluster3IsoMuonVetoEta[i] = 0.0;
        cscRechitCluster3GenMuonVetoE[i] = 0.0;
        cscRechitCluster3GenMuonVetoPt[i] = 0.0;
        cscRechitCluster3JetVetoPt_0p6[i] = 0.0;
        cscRechitCluster3JetVetoPt_0p8[i] = 0.0;
        cscRechitCluster3JetVetoE_0p6[i] = 0.0;
        cscRechitCluster3JetVetoE_0p8[i] = 0.0;
        cscRechitCluster3MuonVetoPt_0p6[i] = 0.0;
        cscRechitCluster3MuonVetoPt_0p8[i] = 0.0;
        cscRechitCluster3MuonVetoE_0p6[i] = 0.0;
        cscRechitCluster3MuonVetoE_0p8[i] = 0.0;

        cscRechitCluster3ZLep1[i] = false;
        cscRechitCluster3ZLep2[i] = false;
        cscRechitCluster3ZLep1Id[i] = -999;
        cscRechitCluster3ZLep2Id[i] = -999;
        cscRechitCluster3ZLep1LooseIso[i] = false;
        cscRechitCluster3ZLep1TightIso[i] = false;
        cscRechitCluster3ZLep1VTightIso[i] = false;
        cscRechitCluster3ZLep1VVTightIso[i] = false;
        cscRechitCluster3ZLep1TightId[i] = false;
        cscRechitCluster3ZLep2LooseIso[i] = false;
        cscRechitCluster3ZLep2TightIso[i] = false;
        cscRechitCluster3ZLep2VTightIso[i] = false;
        cscRechitCluster3ZLep2VVTightIso[i] = false;
        cscRechitCluster3ZLep2TightId[i] = false;
        cscRechitCluster3NChamber[i] = -999;
        cscRechitCluster3MaxChamberRatio[i] = -999.;
        cscRechitCluster3MaxChamber[i] = -999;
        cscRechitCluster3NStation[i] = -999;
        cscRechitCluster3NStation5[i] = -999;
        cscRechitCluster3NStation10perc[i] = -999;
        cscRechitCluster3AvgStation[i] = -999.;
        cscRechitCluster3AvgStation5[i] = -999.;
        cscRechitCluster3AvgStation10perc[i] = -999.;
        cscRechitCluster3MaxStationRatio[i] = -999.;
        cscRechitCluster3MaxStation[i] = -999;
        cscRechitCluster3Me11Ratio[i] = -999.;
        cscRechitCluster3Me12Ratio[i] = -999.;
        cscRechitCluster3NRechitChamberPlus11[i] = -999;
        cscRechitCluster3NRechitChamberPlus12[i] = -999;
        cscRechitCluster3NRechitChamberPlus13[i] = -999;
        cscRechitCluster3NRechitChamberPlus21[i] = -999;
        cscRechitCluster3NRechitChamberPlus22[i] = -999;
        cscRechitCluster3NRechitChamberPlus31[i] = -999;
        cscRechitCluster3NRechitChamberPlus32[i] = -999;
        cscRechitCluster3NRechitChamberPlus41[i] = -999;
        cscRechitCluster3NRechitChamberPlus42[i] = -999;
        cscRechitCluster3NRechitChamberMinus11[i] = -999;
        cscRechitCluster3NRechitChamberMinus12[i] = -999;
        cscRechitCluster3NRechitChamberMinus13[i] = -999;
        cscRechitCluster3NRechitChamberMinus21[i] = -999;
        cscRechitCluster3NRechitChamberMinus22[i] = -999;
        cscRechitCluster3NRechitChamberMinus31[i] = -999;
        cscRechitCluster3NRechitChamberMinus32[i] = -999;
        cscRechitCluster3NRechitChamberMinus41[i] = -999;
        cscRechitCluster3NRechitChamberMinus42[i] = -999;
        cscRechitCluster3Met_dPhi[i] = 999.;
        cscRechitCluster3MetXYCorr_dPhi[i] = 999.;
  }

  for(int i = 0;i<2;i++)
  {
    gLLP_eta[i] = 0.0;
    gLLP_phi[i] = 0.0;
    gLLP_beta[i] = 0.0;
    gLLP_csc[i] = 0.0;
    gLLP_ctau[i] = 0.0;
    gLLP_decay_vertex_r[i] = 0.0;
    gLLP_decay_vertex_x[i] = 0.0;
    gLLP_decay_vertex_y[i] = 0.0;
    gLLP_decay_vertex_z[i] = 0.0;



  }
  genMetPtTrue = -999.;
  genMetPhiTrue = -999.;
  genMetPtCalo = -999.;
  genMetPhiCalo = -999.;
  nGenParticle = 0;
  nGenJets = 0;
  for( int i = 0; i < N_MAX_GPARTICLES; i++ )
  {
    gParticleId[i] = 0;
    gParticleStatus[i] = 999;
    gParticleMotherId[i] = 0;
    gParticlePt[i] = -999.;
    gParticleEta[i] = -999.;
    gParticlePhi[i] = -999.;
    gParticleE[i] = -999.;
    genJetE[i] = -999.;
    genJetPt[i] = -999.;
    genJetEta[i] = -999.;
    genJetPhi[i] = -999.;
    genJetMET[i] = -999.;

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
    lepPassVetoId[i] = false;
    lepPassId[i] = false;
    lepPassLooseIso[i] = false;
    lepPassTightIso[i] = false;
    lepPassVTightIso[i] = false;
    lepPassVTightIso[i] = false;

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
    jet_match_genJet_minDeltaR[i] = -999.;
    jet_match_genJet_index[i] = -999;
    jet_match_genJet_pt[i] = -999.;
    jetTightPassId[i] = false;
  }

  for(int i = 0; i <NTriggersMAX; i++){
    HLTDecision[i] = false;
  }

};

void LiteTreeMuonSystem::InitTree()
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
  tree_->SetBranchAddress("higgsPtWeight",      &higgsPtWeight);
  tree_->SetBranchAddress("higgsPtWeightSys",      higgsPtWeightSys);
  tree_->SetBranchAddress("sf_facScaleUp",      &sf_facScaleUp);
  tree_->SetBranchAddress("sf_facScaleDown",      &sf_facScaleDown);
  tree_->SetBranchAddress("sf_renScaleUp",      &sf_renScaleUp);
  tree_->SetBranchAddress("sf_renScaleDown",      &sf_renScaleDown);
  tree_->SetBranchAddress("sf_facRenScaleUp",      &sf_facRenScaleUp);
  tree_->SetBranchAddress("sf_facRenScaleDown",      &sf_facRenScaleDown);

  tree_->SetBranchAddress("pileupWeight",      &pileupWeight);
  tree_->SetBranchAddress("pileupWeightUp",      &pileupWeightUp);
  tree_->SetBranchAddress("pileupWeightDown",      &pileupWeightDown);
  tree_->SetBranchAddress("Flag_HBHENoiseFilter",      &Flag_HBHENoiseFilter);
  tree_->SetBranchAddress("Flag_HBHEIsoNoiseFilter",      &Flag_HBHEIsoNoiseFilter);
  tree_->SetBranchAddress("Flag_BadPFMuonFilter",      &Flag_BadPFMuonFilter);
  tree_->SetBranchAddress("Flag_CSCTightHaloFilter",      &Flag_CSCTightHaloFilter);
  tree_->SetBranchAddress("Flag_BadChargedCandidateFilter",      &Flag_BadChargedCandidateFilter);
  tree_->SetBranchAddress("Flag_eeBadScFilter",      &Flag_eeBadScFilter);
  tree_->SetBranchAddress("Flag_globalSuperTightHalo2016Filter",      &Flag_globalSuperTightHalo2016Filter);
  tree_->SetBranchAddress("Flag_goodVertices",      &Flag_goodVertices);
  tree_->SetBranchAddress("Flag_ecalBadCalibFilter",      &Flag_ecalBadCalibFilter);
  // tree_->SetBranchAddress("Flag_all",      &Flag_all);

  tree_->SetBranchAddress("Flag2_HBHENoiseFilter",      &Flag2_HBHENoiseFilter);
  tree_->SetBranchAddress("Flag2_HBHEIsoNoiseFilter",      &Flag2_HBHEIsoNoiseFilter);
  tree_->SetBranchAddress("Flag2_BadPFMuonFilter",      &Flag2_BadPFMuonFilter);
  tree_->SetBranchAddress("Flag2_globalSuperTightHalo2016Filter",      &Flag2_globalSuperTightHalo2016Filter);
  tree_->SetBranchAddress("Flag2_globalTightHalo2016Filter",      &Flag2_globalTightHalo2016Filter);
  tree_->SetBranchAddress("Flag2_BadChargedCandidateFilter",      &Flag2_BadChargedCandidateFilter);
  tree_->SetBranchAddress("Flag2_EcalDeadCellTriggerPrimitiveFilter",      &Flag2_EcalDeadCellTriggerPrimitiveFilter);
  tree_->SetBranchAddress("Flag2_ecalBadCalibFilter",      &Flag2_ecalBadCalibFilter);
  tree_->SetBranchAddress("Flag2_eeBadScFilter",      &Flag2_eeBadScFilter);
  tree_->SetBranchAddress("Flag2_all",      &Flag2_all);


  // tree_->SetBranchAddress("rho",         &rho);
  tree_->SetBranchAddress("met",         &met);
  tree_->SetBranchAddress("HT",         &HT);

  tree_->SetBranchAddress("metPhi",      &metPhi);
  tree_->SetBranchAddress("metXYCorr",      &metXYCorr);
  tree_->SetBranchAddress("metPhiXYCorr",      &metPhiXYCorr);

  tree_->SetBranchAddress("jetMet_dPhi",      &jetMet_dPhi);
  tree_->SetBranchAddress("jetMet_dPhiMin",      &jetMet_dPhiMin);
  tree_->SetBranchAddress("jetMet_dPhiMin4",      &jetMet_dPhiMin4);

  tree_->SetBranchAddress("metJESUp",      &metJESUp);
  tree_->SetBranchAddress("metJESDown",      &metJESDown);

  tree_->SetBranchAddress("genMetPtTrue",         &genMetPtTrue);
  tree_->SetBranchAddress("genMetPhiTrue",      &genMetPhiTrue);
  tree_->SetBranchAddress("genMetPtCalo",      &genMetPtCalo);
  tree_->SetBranchAddress("genMetPhiCalo",      &genMetPhiCalo);

  tree_->SetBranchAddress("nGenParticle",      &nGenParticle);
  tree_->SetBranchAddress("gParticleId",      &gParticleId);
  tree_->SetBranchAddress("gParticleStatus",      &gParticleStatus);
  tree_->SetBranchAddress("gParticleMotherId",      &gParticleMotherId);
  tree_->SetBranchAddress("gParticleE",      &gParticleE);
  tree_->SetBranchAddress("gParticlePt",      &gParticlePt);
  tree_->SetBranchAddress("gParticleEta",      &gParticleEta);
  tree_->SetBranchAddress("gParticlePhi",      &gParticlePhi);

  tree_->SetBranchAddress("nGenJets",      &nGenJets);
  tree_->SetBranchAddress("genJetE",      &genJetE);
  tree_->SetBranchAddress("genJetPt",      &genJetPt);
  tree_->SetBranchAddress("genJetEta",      &genJetEta);
  tree_->SetBranchAddress("genJetPhi",      &genJetPhi);
  tree_->SetBranchAddress("genJetMET",      &genJetMET);


  tree_->SetBranchAddress("gLepId",      &gLepId);
  tree_->SetBranchAddress("gLepPt",      &gLepPt);
  tree_->SetBranchAddress("gLepPhi",      &gLepPhi);
  tree_->SetBranchAddress("gLepE",      &gLepE);
  tree_->SetBranchAddress("gLepEta",      &gLepEta);
  tree_->SetBranchAddress("gHiggsPt",      &gHiggsPt);
  tree_->SetBranchAddress("gHiggsPhi",      &gHiggsPhi);
  tree_->SetBranchAddress("gHiggsE",      &gHiggsE);
  tree_->SetBranchAddress("gHiggsEta",      &gHiggsEta);
  //CSC
  tree_->SetBranchAddress("nCscRechits",             &nCscRechits);
  tree_->SetBranchAddress("nCscPositiveYRechits",             &nCscPositiveYRechits);
  tree_->SetBranchAddress("nCscNegativeYRechits",             &nCscNegativeYRechits);
  tree_->SetBranchAddress("cscPosTpeak",             &cscPosTpeak);
  tree_->SetBranchAddress("cscNegTpeak",             &cscNegTpeak);


  tree_->SetBranchAddress("nEarlyCscRechits",             &nEarlyCscRechits);
  tree_->SetBranchAddress("nLateCscRechits",             &nLateCscRechits);
  tree_->SetBranchAddress("nEarly2CscRechits",             &nEarly2CscRechits);
  tree_->SetBranchAddress("nLate2CscRechits",             &nLate2CscRechits);
  tree_->SetBranchAddress("nLate2CscRechits",             &nLate2CscRechits);
  tree_->SetBranchAddress("nCscRings",             &nCscRings);

  tree_->SetBranchAddress("nCscRechitsChamberPlus11",           &nCscRechitsChamberPlus11);
  tree_->SetBranchAddress("nCscRechitsChamberPlus12",           &nCscRechitsChamberPlus12);
  tree_->SetBranchAddress("nCscRechitsChamberPlus13",           &nCscRechitsChamberPlus13);
  tree_->SetBranchAddress("nCscRechitsChamberPlus21",           &nCscRechitsChamberPlus21);
  tree_->SetBranchAddress("nCscRechitsChamberPlus22",           &nCscRechitsChamberPlus22);
  tree_->SetBranchAddress("nCscRechitsChamberPlus31",           &nCscRechitsChamberPlus31);
  tree_->SetBranchAddress("nCscRechitsChamberPlus32",           &nCscRechitsChamberPlus32);
  tree_->SetBranchAddress("nCscRechitsChamberPlus41",           &nCscRechitsChamberPlus41);
  tree_->SetBranchAddress("nCscRechitsChamberPlus42",           &nCscRechitsChamberPlus42);

  tree_->SetBranchAddress("nCscRechitsChamberMinus11",            &nCscRechitsChamberMinus11);
  tree_->SetBranchAddress("nCscRechitsChamberMinus12",            &nCscRechitsChamberMinus12);
  tree_->SetBranchAddress("nCscRechitsChamberMinus13",            &nCscRechitsChamberMinus13);
  tree_->SetBranchAddress("nCscRechitsChamberMinus21",            &nCscRechitsChamberMinus21);
  tree_->SetBranchAddress("nCscRechitsChamberMinus22",            &nCscRechitsChamberMinus22);
  tree_->SetBranchAddress("nCscRechitsChamberMinus31",            &nCscRechitsChamberMinus31);
  tree_->SetBranchAddress("nCscRechitsChamberMinus32",            &nCscRechitsChamberMinus32);
  tree_->SetBranchAddress("nCscRechitsChamberMinus41",            &nCscRechitsChamberMinus41);
  tree_->SetBranchAddress("nCscRechitsChamberMinus42",            &nCscRechitsChamberMinus42);


  tree_->SetBranchAddress("nDTRechits",            &nDTRechits);
  tree_->SetBranchAddress("nDTPositiveYRechits",            &nDTPositiveYRechits);
  tree_->SetBranchAddress("nDTNegativeYRechits",            &nDTNegativeYRechits);
  tree_->SetBranchAddress("nDtRings",             &nDtRings);

  tree_->SetBranchAddress("nDTRechitsChamberMinus12",            &nDTRechitsChamberMinus12);
  tree_->SetBranchAddress("nDTRechitsChamberMinus11",            &nDTRechitsChamberMinus11);
  tree_->SetBranchAddress("nDTRechitsChamber10",            &nDTRechitsChamber10);
  tree_->SetBranchAddress("nDTRechitsChamberPlus11",            &nDTRechitsChamberPlus11);
  tree_->SetBranchAddress("nDTRechitsChamberPlus12",            &nDTRechitsChamberPlus12);
  tree_->SetBranchAddress("nDTRechitsChamberMinus22",            &nDTRechitsChamberMinus22);
  tree_->SetBranchAddress("nDTRechitsChamberMinus21",            &nDTRechitsChamberMinus21);
  tree_->SetBranchAddress("nDTRechitsChamber20",            &nDTRechitsChamber20);
  tree_->SetBranchAddress("nDTRechitsChamberPlus21",            &nDTRechitsChamberPlus21);
  tree_->SetBranchAddress("nDTRechitsChamberPlus22",            &nDTRechitsChamberPlus22);
  tree_->SetBranchAddress("nDTRechitsChamberMinus32",            &nDTRechitsChamberMinus32);
  tree_->SetBranchAddress("nDTRechitsChamberMinus31",            &nDTRechitsChamberMinus31);
  tree_->SetBranchAddress("nDTRechitsChamber30",            &nDTRechitsChamber30);

  tree_->SetBranchAddress("nDTRechitsChamberPlus31",            &nDTRechitsChamberPlus31);
  tree_->SetBranchAddress("nDTRechitsChamberPlus32",            &nDTRechitsChamberPlus32);
  tree_->SetBranchAddress("nDTRechitsChamberMinus42",            &nDTRechitsChamberMinus42);
  tree_->SetBranchAddress("nDTRechitsChamberMinus41",            &nDTRechitsChamberMinus41);
  tree_->SetBranchAddress("nDTRechitsChamber40",            &nDTRechitsChamber40);
  tree_->SetBranchAddress("nDTRechitsChamberPlus41",            &nDTRechitsChamberPlus41);
  tree_->SetBranchAddress("nDTRechitsChamberPlus42",            &nDTRechitsChamberPlus42);

  tree_->SetBranchAddress("cscRechitsStation",             cscRechitsStation);
  tree_->SetBranchAddress("cscRechitsChamber",             cscRechitsChamber);

  tree_->SetBranchAddress("cscRechitsPhi",           cscRechitsPhi);
  tree_->SetBranchAddress("cscRechitsEta",           cscRechitsEta);
  tree_->SetBranchAddress("cscRechitsQuality",           cscRechitsQuality);
  tree_->SetBranchAddress("cscRechitsX",             cscRechitsX);
  tree_->SetBranchAddress("cscRechitsY",             cscRechitsY);
  tree_->SetBranchAddress("cscRechitsZ",             cscRechitsZ);
  tree_->SetBranchAddress("cscRechitsTpeak",             cscRechitsTpeak);
  tree_->SetBranchAddress("cscRechitsTwire",             cscRechitsTwire);
  tree_->SetBranchAddress("cscRechitsClusterId",             cscRechitsClusterId);
  tree_->SetBranchAddress("cscRechitsCluster2Id",             cscRechitsCluster2Id);

  tree_->SetBranchAddress("nCscClusters",             &nCscClusters);
  tree_->SetBranchAddress("cscCluster_match_gLLP",             &cscCluster_match_gLLP);
  tree_->SetBranchAddress("cscCluster_match_gLLP_minDeltaR",             &cscCluster_match_gLLP_minDeltaR);
  tree_->SetBranchAddress("cscCluster_match_gLLP_index",             &cscCluster_match_gLLP_index);

  // tree_->SetBranchAddress("nCsc_JetVetoCluster0p4",             &nCsc_JetVetoCluster0p4);
  // tree_->SetBranchAddress("nCsc_JetMuonVetoCluster0p4",             &nCsc_JetMuonVetoCluster0p4);
  // tree_->SetBranchAddress("nCsc_JetVetoCluster0p4_Me1112Veto",             &nCsc_JetVetoCluster0p4_Me1112Veto);
  // tree_->SetBranchAddress("nCsc_JetMuonVetoCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoCluster0p4_Me1112Veto);
  tree_->SetBranchAddress("cscClusterMe11Ratio",             &cscClusterMe11Ratio);
  tree_->SetBranchAddress("cscClusterMe12Ratio",             &cscClusterMe12Ratio);
  tree_->SetBranchAddress("cscClusterMaxStation",             &cscClusterMaxStation);
  tree_->SetBranchAddress("cscClusterMaxStationRatio",             &cscClusterMaxStationRatio);
  tree_->SetBranchAddress("cscClusterNStation",             &cscClusterNStation);
  tree_->SetBranchAddress("cscClusterMaxChamber",             &cscClusterMaxChamber);
  tree_->SetBranchAddress("cscClusterMaxChamberRatio",             &cscClusterMaxChamberRatio);
  tree_->SetBranchAddress("cscClusterNChamber",             &cscClusterNChamber);
  tree_->SetBranchAddress("cscClusterX",             cscClusterX);
  tree_->SetBranchAddress("cscClusterY",             cscClusterY);
  tree_->SetBranchAddress("cscClusterZ",             cscClusterZ);
  tree_->SetBranchAddress("cscClusterTime",             cscClusterTime);
  tree_->SetBranchAddress("cscClusterGenMuonDeltaR",             cscClusterGenMuonDeltaR);

  tree_->SetBranchAddress("cscClusterTimeSpread",             cscClusterTimeSpread);
  tree_->SetBranchAddress("cscClusterMajorAxis",             cscClusterMajorAxis);
  tree_->SetBranchAddress("cscClusterMinorAxis",             cscClusterMinorAxis);
  tree_->SetBranchAddress("cscClusterXSpread",             cscClusterXSpread);
  tree_->SetBranchAddress("cscClusterYSpread",             cscClusterYSpread);
  tree_->SetBranchAddress("cscClusterZSpread",             cscClusterZSpread);
  tree_->SetBranchAddress("cscClusterEtaPhiSpread",             cscClusterEtaPhiSpread);
  tree_->SetBranchAddress("cscClusterEtaSpread",             cscClusterEtaSpread);
  tree_->SetBranchAddress("cscClusterPhiSpread",             cscClusterPhiSpread);
  tree_->SetBranchAddress("cscClusterEta",             cscClusterEta);
  tree_->SetBranchAddress("cscClusterPhi",             cscClusterPhi);
  tree_->SetBranchAddress("cscClusterJetVetoPt",             cscClusterJetVetoPt);
  tree_->SetBranchAddress("cscClusterMuonVetoPt",             cscClusterMuonVetoPt);
  tree_->SetBranchAddress("cscClusterJetVetoE",             cscClusterJetVetoE);
  tree_->SetBranchAddress("cscClusterMuonVetoE",             cscClusterMuonVetoE);
  tree_->SetBranchAddress("cscClusterSize",             cscClusterSize);

  tree_->SetBranchAddress("cscClusterNSegmentChamberPlus11",             cscClusterNSegmentChamberPlus11);
  tree_->SetBranchAddress("cscClusterNSegmentChamberPlus12",             cscClusterNSegmentChamberPlus12);
  tree_->SetBranchAddress("cscClusterNSegmentChamberPlus13",             cscClusterNSegmentChamberPlus13);
  tree_->SetBranchAddress("cscClusterNSegmentChamberPlus21",             cscClusterNSegmentChamberPlus21);
  tree_->SetBranchAddress("cscClusterNSegmentChamberPlus22",             cscClusterNSegmentChamberPlus22);
  tree_->SetBranchAddress("cscClusterNSegmentChamberPlus31",             cscClusterNSegmentChamberPlus31);
  tree_->SetBranchAddress("cscClusterNSegmentChamberPlus32",             cscClusterNSegmentChamberPlus32);
  tree_->SetBranchAddress("cscClusterNSegmentChamberPlus41",             cscClusterNSegmentChamberPlus41);
  tree_->SetBranchAddress("cscClusterNSegmentChamberPlus42",             cscClusterNSegmentChamberPlus42);

  tree_->SetBranchAddress("cscClusterNSegmentChamberMinus11",             cscClusterNSegmentChamberMinus11);
  tree_->SetBranchAddress("cscClusterNSegmentChamberMinus12",             cscClusterNSegmentChamberMinus12);
  tree_->SetBranchAddress("cscClusterNSegmentChamberMinus13",             cscClusterNSegmentChamberMinus13);
  tree_->SetBranchAddress("cscClusterNSegmentChamberMinus21",             cscClusterNSegmentChamberMinus21);
  tree_->SetBranchAddress("cscClusterNSegmentChamberMinus22",             cscClusterNSegmentChamberMinus22);
  tree_->SetBranchAddress("cscClusterNSegmentChamberMinus31",             cscClusterNSegmentChamberMinus31);
  tree_->SetBranchAddress("cscClusterNSegmentChamberMinus32",             cscClusterNSegmentChamberMinus32);
  tree_->SetBranchAddress("cscClusterNSegmentChamberMinus41",             cscClusterNSegmentChamberMinus41);
  tree_->SetBranchAddress("cscClusterNSegmentChamberMinus42",             cscClusterNSegmentChamberMinus42);


  // CSC CLUSTER

  /*tree_->SetBranchAddress("nCscSegClusters",             &nCscSegClusters);
  tree_->SetBranchAddress("cscSegClusterMet_dPhi",             &cscSegClusterMet_dPhi);


  // tree_->SetBranchAddress("nCsc_JetVetoSegCluster0p4",             &nCsc_JetVetoSegCluster0p4);
  // tree_->SetBranchAddress("nCsc_JetMuonVetoSegCluster0p4",             &nCsc_JetMuonVetoSegCluster0p4);
  // tree_->SetBranchAddress("nCsc_JetVetoSegCluster0p4_Me1112Veto",             &nCsc_JetVetoSegCluster0p4_Me1112Veto);
  // tree_->SetBranchAddress("nCsc_JetMuonVetoSegCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoCluster0p4_Me1112Veto);
  tree_->SetBranchAddress("cscSegClusterMe11Ratio",             &cscSegClusterMe11Ratio);
  tree_->SetBranchAddress("cscSegClusterMe12Ratio",             &cscSegClusterMe12Ratio);
  tree_->SetBranchAddress("cscSegClusterMaxStation",             &cscSegClusterMaxStation);
  tree_->SetBranchAddress("cscSegClusterMaxStationRatio",             &cscSegClusterMaxStationRatio);
  tree_->SetBranchAddress("cscSegClusterNStation",             &cscSegClusterNStation);
  tree_->SetBranchAddress("cscSegClusterMaxChamber",             &cscSegClusterMaxChamber);
  tree_->SetBranchAddress("cscSegClusterMaxChamberRatio",             &cscSegClusterMaxChamberRatio);
  tree_->SetBranchAddress("cscSegClusterNChamber",             &cscSegClusterNChamber);
  tree_->SetBranchAddress("cscSegClusterX",             cscSegClusterX);
  tree_->SetBranchAddress("cscSegClusterY",             cscSegClusterY);
  tree_->SetBranchAddress("cscSegClusterZ",             cscSegClusterZ);
  tree_->SetBranchAddress("cscSegClusterTime",             cscSegClusterTime);
  tree_->SetBranchAddress("cscSegClusterGenMuonDeltaR",             cscSegClusterGenMuonDeltaR);

  tree_->SetBranchAddress("cscSegClusterTimeSpread",             cscSegClusterTimeSpread);
  tree_->SetBranchAddress("cscSegClusterMajorAxis",             cscSegClusterMajorAxis);
  tree_->SetBranchAddress("cscSegClusterMinorAxis",             cscSegClusterMinorAxis);
  tree_->SetBranchAddress("cscSegClusterXSpread",             cscSegClusterXSpread);
  tree_->SetBranchAddress("cscSegClusterYSpread",             cscSegClusterYSpread);
  tree_->SetBranchAddress("cscSegClusterZSpread",             cscSegClusterZSpread);
  tree_->SetBranchAddress("cscSegClusterEtaPhiSpread",             cscSegClusterEtaPhiSpread);
  tree_->SetBranchAddress("cscSegClusterEtaSpread",             cscSegClusterEtaSpread);
  tree_->SetBranchAddress("cscSegClusterPhiSpread",             cscSegClusterPhiSpread);
  tree_->SetBranchAddress("cscSegClusterEta",             cscSegClusterEta);
  tree_->SetBranchAddress("cscSegClusterPhi",             cscSegClusterPhi);
  tree_->SetBranchAddress("cscSegClusterJetVetoPt",             cscSegClusterJetVetoPt);
  tree_->SetBranchAddress("cscSegClusterMuonVetoPt",             cscSegClusterMuonVetoPt);
  tree_->SetBranchAddress("cscSegClusterJetVetoE",             cscSegClusterJetVetoE);
  tree_->SetBranchAddress("cscSegClusterMuonVetoE",             cscSegClusterMuonVetoE);
  tree_->SetBranchAddress("cscSegClusterSize",             cscSegClusterSize);

  tree_->SetBranchAddress("cscSegClusterNSegmentChamberPlus11",             cscSegClusterNSegmentChamberPlus11);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberPlus12",             cscSegClusterNSegmentChamberPlus12);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberPlus13",             cscSegClusterNSegmentChamberPlus13);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberPlus21",             cscSegClusterNSegmentChamberPlus21);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberPlus22",             cscSegClusterNSegmentChamberPlus22);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberPlus31",             cscSegClusterNSegmentChamberPlus31);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberPlus32",             cscSegClusterNSegmentChamberPlus32);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberPlus41",             cscSegClusterNSegmentChamberPlus41);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberPlus42",             cscSegClusterNSegmentChamberPlus42);

  tree_->SetBranchAddress("cscSegClusterNSegmentChamberMinus11",             cscSegClusterNSegmentChamberMinus11);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberMinus12",             cscSegClusterNSegmentChamberMinus12);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberMinus13",             cscSegClusterNSegmentChamberMinus13);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberMinus21",             cscSegClusterNSegmentChamberMinus21);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberMinus22",             cscSegClusterNSegmentChamberMinus22);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberMinus31",             cscSegClusterNSegmentChamberMinus31);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberMinus32",             cscSegClusterNSegmentChamberMinus32);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberMinus41",             cscSegClusterNSegmentChamberMinus41);
  tree_->SetBranchAddress("cscSegClusterNSegmentChamberMinus42",             cscSegClusterNSegmentChamberMinus42);
  tree_->SetBranchAddress("cscSegCluster_match_gParticle_id",             &cscSegCluster_match_gParticle_id);
  tree_->SetBranchAddress("cscSegCluster_match_gParticle_minDeltaR",             &cscSegCluster_match_gParticle_minDeltaR);
  tree_->SetBranchAddress("cscSegCluster_match_gParticle_index",             &cscSegCluster_match_gParticle_index);
  */
  // CSC CLUSTER

  tree_->SetBranchAddress("nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto);
  tree_->SetBranchAddress("nCscRechitClusters",             &nCscRechitClusters);
  tree_->SetBranchAddress("cscRechitCluster_match_gParticle_id",             &cscRechitCluster_match_gParticle_id);

  tree_->SetBranchAddress("cscRechitCluster_match_gLLP",             &cscRechitCluster_match_gLLP);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_index",             &cscRechitCluster_match_gLLP_index);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_minDeltaR",             &cscRechitCluster_match_gLLP_minDeltaR);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_eta",             &cscRechitCluster_match_gLLP_eta);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_phi",             &cscRechitCluster_match_gLLP_phi);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_decay_r",             &cscRechitCluster_match_gLLP_decay_r);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_decay_x",             &cscRechitCluster_match_gLLP_decay_x);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_decay_y",             &cscRechitCluster_match_gLLP_decay_y);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_decay_z",             &cscRechitCluster_match_gLLP_decay_z);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_ctau",             &cscRechitCluster_match_gLLP_ctau);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_beta",             &cscRechitCluster_match_gLLP_beta);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_csc",             &cscRechitCluster_match_gLLP_csc);

  tree_->SetBranchAddress("cscRechitClusterMe11Ratio",             &cscRechitClusterMe11Ratio);
  tree_->SetBranchAddress("cscRechitClusterMe12Ratio",             &cscRechitClusterMe12Ratio);
  tree_->SetBranchAddress("cscRechitClusterMaxStation",             &cscRechitClusterMaxStation);
  tree_->SetBranchAddress("cscRechitClusterMaxStationRatio",             &cscRechitClusterMaxStationRatio);
  tree_->SetBranchAddress("cscRechitClusterNStation",             &cscRechitClusterNStation);
  tree_->SetBranchAddress("cscRechitClusterNStation5",             &cscRechitClusterNStation5);
  tree_->SetBranchAddress("cscRechitClusterNStation10perc",             &cscRechitClusterNStation10perc);
  tree_->SetBranchAddress("cscRechitClusterAvgStation",             &cscRechitClusterAvgStation);
  tree_->SetBranchAddress("cscRechitClusterAvgStation5",             &cscRechitClusterAvgStation5);
  tree_->SetBranchAddress("cscRechitClusterAvgStation10perc",             &cscRechitClusterAvgStation10perc);

  tree_->SetBranchAddress("cscRechitClusterMaxChamber",             &cscRechitClusterMaxChamber);
  tree_->SetBranchAddress("cscRechitClusterMaxChamberRatio",             &cscRechitClusterMaxChamberRatio);
  tree_->SetBranchAddress("cscRechitClusterNChamber",             &cscRechitClusterNChamber);
  tree_->SetBranchAddress("cscRechitClusterX",             cscRechitClusterX);
  tree_->SetBranchAddress("cscRechitClusterY",             cscRechitClusterY);
  tree_->SetBranchAddress("cscRechitClusterZ",             cscRechitClusterZ);
  tree_->SetBranchAddress("cscRechitClusterTime",             cscRechitClusterTime);
  tree_->SetBranchAddress("cscRechitClusterTimeTotal",             cscRechitClusterTimeTotal);

  tree_->SetBranchAddress("cscRechitClusterGenMuonDeltaR",             cscRechitClusterGenMuonDeltaR);

  tree_->SetBranchAddress("cscRechitClusterTimeSpread",             cscRechitClusterTimeSpread);
  tree_->SetBranchAddress("cscRechitClusterMajorAxis",             cscRechitClusterMajorAxis);
  tree_->SetBranchAddress("cscRechitClusterMinorAxis",             cscRechitClusterMinorAxis);
  tree_->SetBranchAddress("cscRechitClusterXYSpread",             cscRechitClusterXYSpread);

  tree_->SetBranchAddress("cscRechitClusterXSpread",             cscRechitClusterXSpread);
  tree_->SetBranchAddress("cscRechitClusterYSpread",             cscRechitClusterYSpread);
  tree_->SetBranchAddress("cscRechitClusterZSpread",             cscRechitClusterZSpread);
  tree_->SetBranchAddress("cscRechitClusterEtaPhiSpread",             cscRechitClusterEtaPhiSpread);
  tree_->SetBranchAddress("cscRechitClusterEtaSpread",             cscRechitClusterEtaSpread);
  tree_->SetBranchAddress("cscRechitClusterPhiSpread",             cscRechitClusterPhiSpread);


    tree_->SetBranchAddress("cscRechitClusterXYSpread_corr",             cscRechitClusterXYSpread_corr);
    tree_->SetBranchAddress("cscRechitClusterXSpread_corr",             cscRechitClusterXSpread_corr);
    tree_->SetBranchAddress("cscRechitClusterYSpread_corr",             cscRechitClusterYSpread_corr);
    tree_->SetBranchAddress("cscRechitClusterPhiSpread_corr",             cscRechitClusterPhiSpread_corr);
    tree_->SetBranchAddress("cscRechitClusterEtaSpread_corr",             cscRechitClusterEtaSpread_corr);
    tree_->SetBranchAddress("cscRechitClusterEtaPhiSpread_corr",             cscRechitClusterEtaPhiSpread_corr);
    tree_->SetBranchAddress("cscRechitClusterRSpread_corr",             cscRechitClusterRSpread_corr);

  tree_->SetBranchAddress("cscRechitClusterEta",             cscRechitClusterEta);
  tree_->SetBranchAddress("cscRechitClusterPhi",             cscRechitClusterPhi);
  tree_->SetBranchAddress("cscRechitClusterJetVetoPt",             cscRechitClusterJetVetoPt);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoPt",             cscRechitClusterMuonVetoPt);
  tree_->SetBranchAddress("cscRechitClusterJetVetoE",             cscRechitClusterJetVetoE);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoE",             cscRechitClusterMuonVetoE);
  tree_->SetBranchAddress("cscRechitClusterGenMuonVetoE",             cscRechitClusterGenMuonVetoE);
  tree_->SetBranchAddress("cscRechitClusterGenMuonVetoPt",             cscRechitClusterGenMuonVetoPt);

  tree_->SetBranchAddress("cscRechitClusterSize",             cscRechitClusterSize);

  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus11",             cscRechitClusterNRechitChamberPlus11);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus12",             cscRechitClusterNRechitChamberPlus12);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus13",             cscRechitClusterNRechitChamberPlus13);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus21",             cscRechitClusterNRechitChamberPlus21);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus22",             cscRechitClusterNRechitChamberPlus22);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus31",             cscRechitClusterNRechitChamberPlus31);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus32",             cscRechitClusterNRechitChamberPlus32);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus41",             cscRechitClusterNRechitChamberPlus41);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberPlus42",             cscRechitClusterNRechitChamberPlus42);

  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus11",             cscRechitClusterNRechitChamberMinus11);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus12",             cscRechitClusterNRechitChamberMinus12);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus13",             cscRechitClusterNRechitChamberMinus13);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus21",             cscRechitClusterNRechitChamberMinus21);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus22",             cscRechitClusterNRechitChamberMinus22);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus31",             cscRechitClusterNRechitChamberMinus31);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus32",             cscRechitClusterNRechitChamberMinus32);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus41",             cscRechitClusterNRechitChamberMinus41);
  tree_->SetBranchAddress("cscRechitClusterNRechitChamberMinus42",             cscRechitClusterNRechitChamberMinus42);

  tree_->SetBranchAddress("cscRechitCluster_match_MB1Seg_0p4",             cscRechitCluster_match_MB1Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_RB1_0p4",             cscRechitCluster_match_RB1_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_RE12_0p4",             cscRechitCluster_match_RE12_0p4);


  tree_->SetBranchAddress("cscRechitClusterMet_dPhi",             cscRechitClusterMet_dPhi);

  tree_->SetBranchAddress("nCscRechitClusters2",             &nCscRechitClusters2);
  tree_->SetBranchAddress("cscRechitCluster2_match_Me1112_0p4",             &cscRechitCluster2_match_Me1112_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_Me1112_0p6",             &cscRechitCluster2_match_Me1112_0p6);
  tree_->SetBranchAddress("cscRechitCluster2_match_Me1112_0p8",             &cscRechitCluster2_match_Me1112_0p8);
  tree_->SetBranchAddress("cscRechitCluster2_match_Me11_0p4",             &cscRechitCluster2_match_Me11_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_Me11_0p6",             &cscRechitCluster2_match_Me11_0p6);
  tree_->SetBranchAddress("cscRechitCluster2_match_Me11_0p8",             &cscRechitCluster2_match_Me11_0p8);
  tree_->SetBranchAddress("cscRechitCluster2_match_Me12_0p4",             &cscRechitCluster2_match_Me12_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_Me12_0p6",             &cscRechitCluster2_match_Me12_0p6);
  tree_->SetBranchAddress("cscRechitCluster2_match_Me12_0p8",             &cscRechitCluster2_match_Me12_0p8);
  tree_->SetBranchAddress("cscRechitCluster2_match_gParticle_id",             &cscRechitCluster2_match_gParticle_id);
  tree_->SetBranchAddress("cscRechitCluster2_match_gParticle",             &cscRechitCluster2_match_gParticle);
  tree_->SetBranchAddress("cscRechitCluster2_match_gParticle_minDeltaR",             &cscRechitCluster2_match_gParticle_minDeltaR);
  tree_->SetBranchAddress("cscRechitCluster2_match_gParticle_index",             &cscRechitCluster2_match_gParticle_index);
  tree_->SetBranchAddress("cscRechitCluster2_match_gParticle_eta",             &cscRechitCluster2_match_gParticle_eta);
  tree_->SetBranchAddress("cscRechitCluster2_match_gParticle_phi",             &cscRechitCluster2_match_gParticle_phi);
  tree_->SetBranchAddress("cscRechitCluster2_match_gParticle_E",             &cscRechitCluster2_match_gParticle_E);
  tree_->SetBranchAddress("cscRechitCluster2_match_gParticle_pt",             &cscRechitCluster2_match_gParticle_pt);
  tree_->SetBranchAddress("cscRechitCluster2_match_gParticle_MotherId",             &cscRechitCluster2_match_gParticle_MotherId);

  tree_->SetBranchAddress("cscRechitCluster2_match_gLLP",             &cscRechitCluster2_match_gLLP);
  tree_->SetBranchAddress("cscRechitCluster2_match_gLLP_index",             &cscRechitCluster2_match_gLLP_index);
  tree_->SetBranchAddress("cscRechitCluster2_match_gLLP_minDeltaR",             &cscRechitCluster2_match_gLLP_minDeltaR);
  tree_->SetBranchAddress("cscRechitCluster2_match_gLLP_eta",             &cscRechitCluster2_match_gLLP_eta);
  tree_->SetBranchAddress("cscRechitCluster2_match_gLLP_phi",             &cscRechitCluster2_match_gLLP_phi);
  tree_->SetBranchAddress("cscRechitCluster2_match_gLLP_decay_r",             &cscRechitCluster2_match_gLLP_decay_r);
  tree_->SetBranchAddress("cscRechitCluster2_match_gLLP_decay_x",             &cscRechitCluster2_match_gLLP_decay_x);
  tree_->SetBranchAddress("cscRechitCluster2_match_gLLP_decay_y",             &cscRechitCluster2_match_gLLP_decay_y);
  tree_->SetBranchAddress("cscRechitCluster2_match_gLLP_decay_z",             &cscRechitCluster2_match_gLLP_decay_z);
  tree_->SetBranchAddress("cscRechitCluster2_match_gLLP_ctau",             &cscRechitCluster2_match_gLLP_ctau);
  tree_->SetBranchAddress("cscRechitCluster2_match_gLLP_beta",             &cscRechitCluster2_match_gLLP_beta);
  tree_->SetBranchAddress("cscRechitCluster2_match_gLLP_csc",             &cscRechitCluster2_match_gLLP_csc);

  tree_->SetBranchAddress("cscRechitCluster2Me11Ratio",             &cscRechitCluster2Me11Ratio);
  tree_->SetBranchAddress("cscRechitCluster2Me12Ratio",             &cscRechitCluster2Me12Ratio);
  tree_->SetBranchAddress("cscRechitCluster2MaxStation",             &cscRechitCluster2MaxStation);
  tree_->SetBranchAddress("cscRechitCluster2MaxStationRatio",             &cscRechitCluster2MaxStationRatio);
  tree_->SetBranchAddress("cscRechitCluster2NStation",             &cscRechitCluster2NStation);
  tree_->SetBranchAddress("cscRechitCluster2NStation5",             &cscRechitCluster2NStation5);
  tree_->SetBranchAddress("cscRechitCluster2NStation10perc",             &cscRechitCluster2NStation10perc);
  tree_->SetBranchAddress("cscRechitCluster2AvgStation",             &cscRechitCluster2AvgStation);
  tree_->SetBranchAddress("cscRechitCluster2AvgStation5",             &cscRechitCluster2AvgStation5);
  tree_->SetBranchAddress("cscRechitCluster2AvgStation10perc",             &cscRechitCluster2AvgStation10perc);
  tree_->SetBranchAddress("cscRechitCluster2MaxChamber",             &cscRechitCluster2MaxChamber);
  tree_->SetBranchAddress("cscRechitCluster2MaxChamberRatio",             &cscRechitCluster2MaxChamberRatio);
  tree_->SetBranchAddress("cscRechitCluster2NChamber",             &cscRechitCluster2NChamber);
  tree_->SetBranchAddress("cscRechitCluster2X",             cscRechitCluster2X);
  tree_->SetBranchAddress("cscRechitCluster2Y",             cscRechitCluster2Y);
  tree_->SetBranchAddress("cscRechitCluster2Z",             cscRechitCluster2Z);
  tree_->SetBranchAddress("cscRechitCluster2Time",             cscRechitCluster2Time);
  tree_->SetBranchAddress("cscRechitCluster2TimeTotal",             cscRechitCluster2TimeTotal);

  tree_->SetBranchAddress("cscRechitCluster2GenMuonDeltaR",             cscRechitCluster2GenMuonDeltaR);

  tree_->SetBranchAddress("cscRechitCluster2TimeSpread",             cscRechitCluster2TimeSpread);
  tree_->SetBranchAddress("cscRechitCluster2MajorAxis",             cscRechitCluster2MajorAxis);
  tree_->SetBranchAddress("cscRechitCluster2MinorAxis",             cscRechitCluster2MinorAxis);
  tree_->SetBranchAddress("cscRechitCluster2RSpread",             cscRechitCluster2RSpread);

  tree_->SetBranchAddress("cscRechitCluster2XSpread",             cscRechitCluster2XSpread);
  tree_->SetBranchAddress("cscRechitCluster2YSpread",             cscRechitCluster2YSpread);
  tree_->SetBranchAddress("cscRechitCluster2ZSpread",             cscRechitCluster2ZSpread);
  tree_->SetBranchAddress("cscRechitCluster2EtaPhiSpread",             cscRechitCluster2EtaPhiSpread);
  tree_->SetBranchAddress("cscRechitCluster2XYSpread",             cscRechitCluster2XYSpread);
  tree_->SetBranchAddress("cscRechitCluster2EtaSpread",             cscRechitCluster2EtaSpread);
  tree_->SetBranchAddress("cscRechitCluster2PhiSpread",             cscRechitCluster2PhiSpread);
  tree_->SetBranchAddress("cscRechitCluster2Eta",             cscRechitCluster2Eta);
  tree_->SetBranchAddress("cscRechitCluster2Phi",             cscRechitCluster2Phi);

  tree_->SetBranchAddress("cscRechitCluster2XYSpread_corr",             cscRechitCluster2XYSpread_corr);
  tree_->SetBranchAddress("cscRechitCluster2XSpread_corr",             cscRechitCluster2XSpread_corr);
  tree_->SetBranchAddress("cscRechitCluster2YSpread_corr",             cscRechitCluster2YSpread_corr);
  tree_->SetBranchAddress("cscRechitCluster2PhiSpread_corr",             cscRechitCluster2PhiSpread_corr);
  tree_->SetBranchAddress("cscRechitCluster2EtaSpread_corr",             cscRechitCluster2EtaSpread_corr);
  tree_->SetBranchAddress("cscRechitCluster2EtaPhiSpread_corr",             cscRechitCluster2EtaPhiSpread_corr);
  tree_->SetBranchAddress("cscRechitCluster2RSpread_corr",             cscRechitCluster2RSpread_corr);


  tree_->SetBranchAddress("cscRechitCluster2JetVetoPt",             cscRechitCluster2JetVetoPt);
  tree_->SetBranchAddress("cscRechitCluster2JetVetoE",             cscRechitCluster2JetVetoE);
  tree_->SetBranchAddress("cscRechitCluster2MuonVetoPt",             cscRechitCluster2MuonVetoPt);
  tree_->SetBranchAddress("cscRechitCluster2MuonVetoE",             cscRechitCluster2MuonVetoE);

  tree_->SetBranchAddress("cscRechitCluster2JetVetoPt_0p6",             cscRechitCluster2JetVetoPt_0p6);
  tree_->SetBranchAddress("cscRechitCluster2JetVetoPt_0p8",             cscRechitCluster2JetVetoPt_0p8);
  tree_->SetBranchAddress("cscRechitCluster2JetVetoE_0p6",             cscRechitCluster2JetVetoE_0p6);
  tree_->SetBranchAddress("cscRechitCluster2JetVetoE_0p8",             cscRechitCluster2JetVetoE_0p8);
  tree_->SetBranchAddress("cscRechitCluster2MuonVetoPt_0p6",             cscRechitCluster2MuonVetoPt_0p6);
  tree_->SetBranchAddress("cscRechitCluster2MuonVetoPt_0p8",             cscRechitCluster2MuonVetoPt_0p8);
  tree_->SetBranchAddress("cscRechitCluster2MuonVetoE_0p6",             cscRechitCluster2MuonVetoE_0p6);
  tree_->SetBranchAddress("cscRechitCluster2MuonVetoE_0p8",             cscRechitCluster2MuonVetoE_0p8);

  tree_->SetBranchAddress("cscRechitCluster2ZLep1",             cscRechitCluster2ZLep1);
  tree_->SetBranchAddress("cscRechitCluster2ZLep2",             cscRechitCluster2ZLep2);
  tree_->SetBranchAddress("cscRechitCluster2ZLep1Id",             cscRechitCluster2ZLep1Id);
  tree_->SetBranchAddress("cscRechitCluster2ZLep2Id",             cscRechitCluster2ZLep2Id);

  tree_->SetBranchAddress("cscRechitCluster2ZLep1LooseIso",             cscRechitCluster2ZLep1LooseIso);
  tree_->SetBranchAddress("cscRechitCluster2ZLep1TightIso",             cscRechitCluster2ZLep1TightIso);
  tree_->SetBranchAddress("cscRechitCluster2ZLep1VTightIso",             cscRechitCluster2ZLep1VTightIso);
  tree_->SetBranchAddress("cscRechitCluster2ZLep1VVTightIso",             cscRechitCluster2ZLep1VVTightIso);
  tree_->SetBranchAddress("cscRechitCluster2ZLep1TightId",             cscRechitCluster2ZLep1TightId);
  tree_->SetBranchAddress("cscRechitCluster2ZLep2LooseIso",             cscRechitCluster2ZLep2LooseIso);
  tree_->SetBranchAddress("cscRechitCluster2ZLep2TightIso",             cscRechitCluster2ZLep2TightIso);
  tree_->SetBranchAddress("cscRechitCluster2ZLep2VTightIso",             cscRechitCluster2ZLep2VTightIso);
  tree_->SetBranchAddress("cscRechitCluster2ZLep2VVTightIso",             cscRechitCluster2ZLep2VVTightIso);
  tree_->SetBranchAddress("cscRechitCluster2ZLep2TightId",             cscRechitCluster2ZLep2TightId);

  tree_->SetBranchAddress("cscRechitCluster2MuonVetoPhi",             cscRechitCluster2MuonVetoPhi);
  tree_->SetBranchAddress("cscRechitCluster2MuonVetoEta",             cscRechitCluster2MuonVetoEta);
  tree_->SetBranchAddress("cscRechitCluster2MuonVetoIso",             cscRechitCluster2MuonVetoIso);

  tree_->SetBranchAddress("cscRechitCluster2MuonVetoLooseIso",             cscRechitCluster2MuonVetoLooseIso);
  tree_->SetBranchAddress("cscRechitCluster2MuonVetoTightIso",             cscRechitCluster2MuonVetoTightIso);
  tree_->SetBranchAddress("cscRechitCluster2MuonVetoVTightIso",             cscRechitCluster2MuonVetoVTightIso);
  tree_->SetBranchAddress("cscRechitCluster2MuonVetoVVTightIso",             cscRechitCluster2MuonVetoVVTightIso);
  tree_->SetBranchAddress("cscRechitCluster2MuonVetoTightId",             cscRechitCluster2MuonVetoTightId);

  tree_->SetBranchAddress("cscRechitCluster2IsoMuonVetoPt",             cscRechitCluster2IsoMuonVetoPt);
  tree_->SetBranchAddress("cscRechitCluster2IsoMuonVetoE",             cscRechitCluster2IsoMuonVetoE);
  tree_->SetBranchAddress("cscRechitCluster2IsoMuonVetoPhi",             cscRechitCluster2IsoMuonVetoPhi);
  tree_->SetBranchAddress("cscRechitCluster2IsoMuonVetoEta",             cscRechitCluster2IsoMuonVetoEta);
  tree_->SetBranchAddress("cscRechitCluster2GenMuonVetoPt",             cscRechitCluster2GenMuonVetoPt);
  tree_->SetBranchAddress("cscRechitCluster2GenMuonVetoE",             cscRechitCluster2GenMuonVetoE);
  tree_->SetBranchAddress("cscRechitCluster2Size",             cscRechitCluster2Size);

  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberPlus11",             cscRechitCluster2NRechitChamberPlus11);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberPlus12",             cscRechitCluster2NRechitChamberPlus12);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberPlus13",             cscRechitCluster2NRechitChamberPlus13);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberPlus21",             cscRechitCluster2NRechitChamberPlus21);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberPlus22",             cscRechitCluster2NRechitChamberPlus22);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberPlus31",             cscRechitCluster2NRechitChamberPlus31);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberPlus32",             cscRechitCluster2NRechitChamberPlus32);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberPlus41",             cscRechitCluster2NRechitChamberPlus41);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberPlus42",             cscRechitCluster2NRechitChamberPlus42);

  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberMinus11",             cscRechitCluster2NRechitChamberMinus11);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberMinus12",             cscRechitCluster2NRechitChamberMinus12);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberMinus13",             cscRechitCluster2NRechitChamberMinus13);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberMinus21",             cscRechitCluster2NRechitChamberMinus21);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberMinus22",             cscRechitCluster2NRechitChamberMinus22);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberMinus31",             cscRechitCluster2NRechitChamberMinus31);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberMinus32",             cscRechitCluster2NRechitChamberMinus32);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberMinus41",             cscRechitCluster2NRechitChamberMinus41);
  tree_->SetBranchAddress("cscRechitCluster2NRechitChamberMinus42",             cscRechitCluster2NRechitChamberMinus42);

  tree_->SetBranchAddress("cscRechitCluster2_match_cscRechits_0p4",             cscRechitCluster2_match_cscRechits_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_cscSeg_0p4",             cscRechitCluster2_match_cscSeg_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_ME11Seg_0p4",             cscRechitCluster2_match_ME11Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_ME12Seg_0p4",             cscRechitCluster2_match_ME12Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_cscSeg_0p6",             cscRechitCluster2_match_cscSeg_0p6);
  tree_->SetBranchAddress("cscRechitCluster2_match_ME11Seg_0p6",             cscRechitCluster2_match_ME11Seg_0p6);
  tree_->SetBranchAddress("cscRechitCluster2_match_ME12Seg_0p6",             cscRechitCluster2_match_ME12Seg_0p6);

  tree_->SetBranchAddress("cscRechitCluster2_match_dtRechits_0p4",             cscRechitCluster2_match_dtRechits_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_MB1_0p4",             cscRechitCluster2_match_MB1_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_dtRechits_0p6",             cscRechitCluster2_match_dtRechits_0p6);
  tree_->SetBranchAddress("cscRechitCluster2_match_MB1_0p6",             cscRechitCluster2_match_MB1_0p6);
  tree_->SetBranchAddress("cscRechitCluster2_match_dtSeg_0p4",             cscRechitCluster2_match_dtSeg_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_MB1Seg_0p4",             cscRechitCluster2_match_MB1Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_dtSeg_0p6",             cscRechitCluster2_match_dtSeg_0p6);
  tree_->SetBranchAddress("cscRechitCluster2_match_MB1Seg_0p6",             cscRechitCluster2_match_MB1Seg_0p6);
  tree_->SetBranchAddress("cscRechitCluster2_match_RB1_0p4",             cscRechitCluster2_match_RB1_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_RE12_0p4",             cscRechitCluster2_match_RE12_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_RB1_0p6",             cscRechitCluster2_match_RB1_0p6);
  tree_->SetBranchAddress("cscRechitCluster2_match_RE12_0p6",             cscRechitCluster2_match_RE12_0p6);
  tree_->SetBranchAddress("cscRechitCluster2_match_highEta_0p4",             cscRechitCluster2_match_highEta_0p4);
  tree_->SetBranchAddress("cscRechitCluster2_match_highEta_0p6",             cscRechitCluster2_match_highEta_0p6);
  tree_->SetBranchAddress("cscRechitCluster2_match_highEta_0p8",             cscRechitCluster2_match_highEta_0p8);
  tree_->SetBranchAddress("cscRechitCluster2_match_cluster_dR",             cscRechitCluster2_match_cluster_dR);
  tree_->SetBranchAddress("cscRechitCluster2_match_cluster_index",             cscRechitCluster2_match_cluster_index);

  tree_->SetBranchAddress("cscRechitCluster2Met_dPhi",             cscRechitCluster2Met_dPhi);
  tree_->SetBranchAddress("cscRechitCluster2MetXYCorr_dPhi",             cscRechitCluster2MetXYCorr_dPhi);

  tree_->SetBranchAddress("nCscRechitClusters3",             &nCscRechitClusters3);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me1112_0p4",             &cscRechitCluster3_match_Me1112_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me1112_0p6",             &cscRechitCluster3_match_Me1112_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me1112_0p8",             &cscRechitCluster3_match_Me1112_0p8);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me11_0p4",             &cscRechitCluster3_match_Me11_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me11_0p6",             &cscRechitCluster3_match_Me11_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me11_0p8",             &cscRechitCluster3_match_Me11_0p8);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me12_0p4",             &cscRechitCluster3_match_Me12_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me12_0p6",             &cscRechitCluster3_match_Me12_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me12_0p8",             &cscRechitCluster3_match_Me12_0p8);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_id",             &cscRechitCluster3_match_gParticle_id);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle",             &cscRechitCluster3_match_gParticle);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_minDeltaR",             &cscRechitCluster3_match_gParticle_minDeltaR);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_index",             &cscRechitCluster3_match_gParticle_index);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_eta",             &cscRechitCluster3_match_gParticle_eta);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_phi",             &cscRechitCluster3_match_gParticle_phi);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_E",             &cscRechitCluster3_match_gParticle_E);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_pt",             &cscRechitCluster3_match_gParticle_pt);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_MotherId",             &cscRechitCluster3_match_gParticle_MotherId);

  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP",             &cscRechitCluster3_match_gLLP);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_index",             &cscRechitCluster3_match_gLLP_index);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_minDeltaR",             &cscRechitCluster3_match_gLLP_minDeltaR);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_eta",             &cscRechitCluster3_match_gLLP_eta);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_phi",             &cscRechitCluster3_match_gLLP_phi);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_decay_r",             &cscRechitCluster3_match_gLLP_decay_r);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_decay_x",             &cscRechitCluster3_match_gLLP_decay_x);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_decay_y",             &cscRechitCluster3_match_gLLP_decay_y);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_decay_z",             &cscRechitCluster3_match_gLLP_decay_z);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_ctau",             &cscRechitCluster3_match_gLLP_ctau);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_beta",             &cscRechitCluster3_match_gLLP_beta);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_csc",             &cscRechitCluster3_match_gLLP_csc);

  tree_->SetBranchAddress("cscRechitCluster3Me11Ratio",             &cscRechitCluster3Me11Ratio);
  tree_->SetBranchAddress("cscRechitCluster3Me12Ratio",             &cscRechitCluster3Me12Ratio);
  tree_->SetBranchAddress("cscRechitCluster3MaxStation",             &cscRechitCluster3MaxStation);
  tree_->SetBranchAddress("cscRechitCluster3MaxStationRatio",             &cscRechitCluster3MaxStationRatio);
  tree_->SetBranchAddress("cscRechitCluster3NStation",             &cscRechitCluster3NStation);
  tree_->SetBranchAddress("cscRechitCluster3NStation5",             &cscRechitCluster3NStation5);
  tree_->SetBranchAddress("cscRechitCluster3NStation10perc",             &cscRechitCluster3NStation10perc);
  tree_->SetBranchAddress("cscRechitCluster3AvgStation",             &cscRechitCluster3AvgStation);
  tree_->SetBranchAddress("cscRechitCluster3AvgStation5",             &cscRechitCluster3AvgStation5);
  tree_->SetBranchAddress("cscRechitCluster3AvgStation10perc",             &cscRechitCluster3AvgStation10perc);
  tree_->SetBranchAddress("cscRechitCluster3MaxChamber",             &cscRechitCluster3MaxChamber);
  tree_->SetBranchAddress("cscRechitCluster3MaxChamberRatio",             &cscRechitCluster3MaxChamberRatio);
  tree_->SetBranchAddress("cscRechitCluster3NChamber",             &cscRechitCluster3NChamber);
  tree_->SetBranchAddress("cscRechitCluster3X",             cscRechitCluster3X);
  tree_->SetBranchAddress("cscRechitCluster3Y",             cscRechitCluster3Y);
  tree_->SetBranchAddress("cscRechitCluster3Z",             cscRechitCluster3Z);
  tree_->SetBranchAddress("cscRechitCluster3Time",             cscRechitCluster3Time);
  tree_->SetBranchAddress("cscRechitCluster3TimeTotal",             cscRechitCluster3TimeTotal);

  tree_->SetBranchAddress("cscRechitCluster3GenMuonDeltaR",             cscRechitCluster3GenMuonDeltaR);

  tree_->SetBranchAddress("cscRechitCluster3TimeSpread",             cscRechitCluster3TimeSpread);
  tree_->SetBranchAddress("cscRechitCluster3MajorAxis",             cscRechitCluster3MajorAxis);
  tree_->SetBranchAddress("cscRechitCluster3MinorAxis",             cscRechitCluster3MinorAxis);
  tree_->SetBranchAddress("cscRechitCluster3RSpread",             cscRechitCluster3RSpread);

  tree_->SetBranchAddress("cscRechitCluster3XSpread",             cscRechitCluster3XSpread);
  tree_->SetBranchAddress("cscRechitCluster3YSpread",             cscRechitCluster3YSpread);
  tree_->SetBranchAddress("cscRechitCluster3ZSpread",             cscRechitCluster3ZSpread);
  tree_->SetBranchAddress("cscRechitCluster3EtaPhiSpread",             cscRechitCluster3EtaPhiSpread);
  tree_->SetBranchAddress("cscRechitCluster3XYSpread",             cscRechitCluster3XYSpread);
  tree_->SetBranchAddress("cscRechitCluster3EtaSpread",             cscRechitCluster3EtaSpread);
  tree_->SetBranchAddress("cscRechitCluster3PhiSpread",             cscRechitCluster3PhiSpread);
  tree_->SetBranchAddress("cscRechitCluster3Eta",             cscRechitCluster3Eta);
  tree_->SetBranchAddress("cscRechitCluster3Phi",             cscRechitCluster3Phi);


  tree_->SetBranchAddress("cscRechitCluster3XYSpread_corr",             cscRechitCluster3XYSpread_corr);
  tree_->SetBranchAddress("cscRechitCluster3XSpread_corr",             cscRechitCluster3XSpread_corr);
  tree_->SetBranchAddress("cscRechitCluster3YSpread_corr",             cscRechitCluster3YSpread_corr);
  tree_->SetBranchAddress("cscRechitCluster3PhiSpread_corr",             cscRechitCluster3PhiSpread_corr);
  tree_->SetBranchAddress("cscRechitCluster3EtaSpread_corr",             cscRechitCluster3EtaSpread_corr);
  tree_->SetBranchAddress("cscRechitCluster3EtaPhiSpread_corr",             cscRechitCluster3EtaPhiSpread_corr);
  tree_->SetBranchAddress("cscRechitCluster3RSpread_corr",             cscRechitCluster3RSpread_corr);

    tree_->SetBranchAddress("cscRechitCluster3XYSpread_corr2",             cscRechitCluster3XYSpread_corr2);
    tree_->SetBranchAddress("cscRechitCluster3XSpread_corr2",             cscRechitCluster3XSpread_corr2);
    tree_->SetBranchAddress("cscRechitCluster3YSpread_corr2",             cscRechitCluster3YSpread_corr2);
    tree_->SetBranchAddress("cscRechitCluster3PhiSpread_corr2",             cscRechitCluster3PhiSpread_corr2);
    tree_->SetBranchAddress("cscRechitCluster3EtaSpread_corr2",             cscRechitCluster3EtaSpread_corr2);
    tree_->SetBranchAddress("cscRechitCluster3EtaPhiSpread_corr2",             cscRechitCluster3EtaPhiSpread_corr2);
    tree_->SetBranchAddress("cscRechitCluster3RSpread_corr2",             cscRechitCluster3RSpread_corr2);

  tree_->SetBranchAddress("cscRechitCluster3JetVetoPt",             cscRechitCluster3JetVetoPt);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoE",             cscRechitCluster3JetVetoE);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoPt",             cscRechitCluster3MuonVetoPt);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoE",             cscRechitCluster3MuonVetoE);

  tree_->SetBranchAddress("cscRechitCluster3JetVetoPt_0p6",             cscRechitCluster3JetVetoPt_0p6);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoPt_0p8",             cscRechitCluster3JetVetoPt_0p8);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoE_0p6",             cscRechitCluster3JetVetoE_0p6);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoE_0p8",             cscRechitCluster3JetVetoE_0p8);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoPt_0p6",             cscRechitCluster3MuonVetoPt_0p6);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoPt_0p8",             cscRechitCluster3MuonVetoPt_0p8);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoE_0p6",             cscRechitCluster3MuonVetoE_0p6);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoE_0p8",             cscRechitCluster3MuonVetoE_0p8);

  tree_->SetBranchAddress("cscRechitCluster3ZLep1",             cscRechitCluster3ZLep1);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2",             cscRechitCluster3ZLep2);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1Id",             cscRechitCluster3ZLep1Id);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2Id",             cscRechitCluster3ZLep2Id);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1LooseIso",             cscRechitCluster3ZLep1LooseIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1TightIso",             cscRechitCluster3ZLep1TightIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1VTightIso",             cscRechitCluster3ZLep1VTightIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1VVTightIso",             cscRechitCluster3ZLep1VVTightIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1TightId",             cscRechitCluster3ZLep1TightId);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2LooseIso",             cscRechitCluster3ZLep2LooseIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2TightIso",             cscRechitCluster3ZLep2TightIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2VTightIso",             cscRechitCluster3ZLep2VTightIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2VVTightIso",             cscRechitCluster3ZLep2VVTightIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2TightId",             cscRechitCluster3ZLep2TightId);

  tree_->SetBranchAddress("cscRechitCluster3MuonVetoPhi",             cscRechitCluster3MuonVetoPhi);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoEta",             cscRechitCluster3MuonVetoEta);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoIso",             cscRechitCluster3MuonVetoIso);

    tree_->SetBranchAddress("cscRechitCluster3MuonVetoLooseIso",             cscRechitCluster3MuonVetoLooseIso);
    tree_->SetBranchAddress("cscRechitCluster3MuonVetoTightIso",             cscRechitCluster3MuonVetoTightIso);
    tree_->SetBranchAddress("cscRechitCluster3MuonVetoVTightIso",             cscRechitCluster3MuonVetoVTightIso);
    tree_->SetBranchAddress("cscRechitCluster3MuonVetoVVTightIso",             cscRechitCluster3MuonVetoVVTightIso);
    tree_->SetBranchAddress("cscRechitCluster3MuonVetoTightId",             cscRechitCluster3MuonVetoTightId);


  tree_->SetBranchAddress("cscRechitCluster3IsoMuonVetoPt",             cscRechitCluster3IsoMuonVetoPt);
  tree_->SetBranchAddress("cscRechitCluster3IsoMuonVetoE",             cscRechitCluster3IsoMuonVetoE);
  tree_->SetBranchAddress("cscRechitCluster3IsoMuonVetoPhi",             cscRechitCluster3IsoMuonVetoPhi);
  tree_->SetBranchAddress("cscRechitCluster3IsoMuonVetoEta",             cscRechitCluster3IsoMuonVetoEta);
  tree_->SetBranchAddress("cscRechitCluster3GenMuonVetoPt",             cscRechitCluster3GenMuonVetoPt);
  tree_->SetBranchAddress("cscRechitCluster3GenMuonVetoE",             cscRechitCluster3GenMuonVetoE);
  tree_->SetBranchAddress("cscRechitCluster3Size",             cscRechitCluster3Size);

  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus11",             cscRechitCluster3NRechitChamberPlus11);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus12",             cscRechitCluster3NRechitChamberPlus12);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus13",             cscRechitCluster3NRechitChamberPlus13);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus21",             cscRechitCluster3NRechitChamberPlus21);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus22",             cscRechitCluster3NRechitChamberPlus22);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus31",             cscRechitCluster3NRechitChamberPlus31);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus32",             cscRechitCluster3NRechitChamberPlus32);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus41",             cscRechitCluster3NRechitChamberPlus41);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus42",             cscRechitCluster3NRechitChamberPlus42);

  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus11",             cscRechitCluster3NRechitChamberMinus11);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus12",             cscRechitCluster3NRechitChamberMinus12);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus13",             cscRechitCluster3NRechitChamberMinus13);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus21",             cscRechitCluster3NRechitChamberMinus21);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus22",             cscRechitCluster3NRechitChamberMinus22);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus31",             cscRechitCluster3NRechitChamberMinus31);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus32",             cscRechitCluster3NRechitChamberMinus32);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus41",             cscRechitCluster3NRechitChamberMinus41);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus42",             cscRechitCluster3NRechitChamberMinus42);

  tree_->SetBranchAddress("cscRechitCluster3_match_cscRechits_0p4",             cscRechitCluster3_match_cscRechits_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_cscSeg_0p4",             cscRechitCluster3_match_cscSeg_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_ME11Seg_0p4",             cscRechitCluster3_match_ME11Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_ME12Seg_0p4",             cscRechitCluster3_match_ME12Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_cscSeg_0p6",             cscRechitCluster3_match_cscSeg_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_ME11Seg_0p6",             cscRechitCluster3_match_ME11Seg_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_ME12Seg_0p6",             cscRechitCluster3_match_ME12Seg_0p6);

  tree_->SetBranchAddress("cscRechitCluster3_match_dtRechits_0p4",             cscRechitCluster3_match_dtRechits_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_MB1_0p4",             cscRechitCluster3_match_MB1_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_dtRechits_0p6",             cscRechitCluster3_match_dtRechits_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_MB1_0p6",             cscRechitCluster3_match_MB1_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_dtSeg_0p4",             cscRechitCluster3_match_dtSeg_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_MB1Seg_0p4",             cscRechitCluster3_match_MB1Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_dtSeg_0p6",             cscRechitCluster3_match_dtSeg_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_MB1Seg_0p6",             cscRechitCluster3_match_MB1Seg_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_RB1_0p4",             cscRechitCluster3_match_RB1_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_RE12_0p4",             cscRechitCluster3_match_RE12_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_RB1_0p6",             cscRechitCluster3_match_RB1_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_RE12_0p6",             cscRechitCluster3_match_RE12_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_highEta_0p4",             cscRechitCluster3_match_highEta_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_highEta_0p6",             cscRechitCluster3_match_highEta_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_highEta_0p8",             cscRechitCluster3_match_highEta_0p8);
  tree_->SetBranchAddress("cscRechitCluster3_match_cluster_dR",             cscRechitCluster3_match_cluster_dR);
  tree_->SetBranchAddress("cscRechitCluster3_match_cluster_index",             cscRechitCluster3_match_cluster_index);

  tree_->SetBranchAddress("cscRechitCluster3Met_dPhi",             cscRechitCluster3Met_dPhi);
  tree_->SetBranchAddress("cscRechitCluster3MetXYCorr_dPhi",             cscRechitCluster3MetXYCorr_dPhi);



  // CSC IT CLUSTER
  /*tree_->SetBranchAddress("nCscITClusters",             &nCscITClusters);
  tree_->SetBranchAddress("nCsc_JetVetoITCluster0p4",             &nCsc_JetVetoITCluster0p4);
  tree_->SetBranchAddress("nCsc_JetMuonVetoITCluster0p4",             &nCsc_JetMuonVetoITCluster0p4);
  tree_->SetBranchAddress("nCsc_JetVetoITCluster0p4_Me1112Veto",             &nCsc_JetVetoITCluster0p4_Me1112Veto);
  tree_->SetBranchAddress("nCsc_JetMuonVetoITCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoITCluster0p4_Me1112Veto);

  tree_->SetBranchAddress("cscITClusterMe11Ratio",             &cscITClusterMe11Ratio);
  tree_->SetBranchAddress("cscITClusterMe12Ratio",             &cscITClusterMe12Ratio);
  tree_->SetBranchAddress("cscITClusterMaxStation",             &cscITClusterMaxStation);
  tree_->SetBranchAddress("cscITClusterMaxStationRatio",             &cscITClusterMaxStationRatio);
  tree_->SetBranchAddress("cscITClusterNStation",             &cscITClusterNStation);
  tree_->SetBranchAddress("cscITClusterMaxChamber",             &cscITClusterMaxChamber);
  tree_->SetBranchAddress("cscITClusterMaxChamberRatio",             &cscITClusterMaxChamberRatio);
  tree_->SetBranchAddress("cscITClusterNChamber",             &cscITClusterNChamber);
  tree_->SetBranchAddress("cscITClusterVertexR",             cscITClusterVertexR);
  tree_->SetBranchAddress("cscITClusterVertexZ",             cscITClusterVertexZ);
  tree_->SetBranchAddress("cscITClusterVertexDis",             cscITClusterVertexDis);
  tree_->SetBranchAddress("cscITClusterVertexChi2",             cscITClusterVertexChi2);
  tree_->SetBranchAddress("cscITClusterVertexN",             cscITClusterVertexN);
  tree_->SetBranchAddress("cscITClusterVertexN1",             cscITClusterVertexN1);
  tree_->SetBranchAddress("cscITClusterVertexN5",             cscITClusterVertexN5);
  tree_->SetBranchAddress("cscITClusterVertexN10",             cscITClusterVertexN10);
  tree_->SetBranchAddress("cscITClusterVertexN15",             cscITClusterVertexN15);
  tree_->SetBranchAddress("cscITClusterVertexN20",             cscITClusterVertexN20);
  tree_->SetBranchAddress("cscITClusterX",             cscITClusterX);
  tree_->SetBranchAddress("cscITClusterY",             cscITClusterY);
  tree_->SetBranchAddress("cscITClusterZ",             cscITClusterZ);
  tree_->SetBranchAddress("cscITClusterTime",             cscITClusterTime);
  tree_->SetBranchAddress("cscITClusterTimeRMS",             cscITClusterTimeRMS);
  tree_->SetBranchAddress("cscITClusterGenMuonDeltaR",             cscITClusterGenMuonDeltaR);

  tree_->SetBranchAddress("cscITClusterTimeSpread",             cscITClusterTimeSpread);
  tree_->SetBranchAddress("cscITClusterRadius",             cscITClusterRadius);
  tree_->SetBranchAddress("cscITClusterMajorAxis",             cscITClusterMajorAxis);
  tree_->SetBranchAddress("cscITClusterMinorAxis",             cscITClusterMinorAxis);
  tree_->SetBranchAddress("cscITClusterXSpread",             cscITClusterXSpread);
  tree_->SetBranchAddress("cscITClusterYSpread",             cscITClusterYSpread);
  tree_->SetBranchAddress("cscITClusterZSpread",             cscITClusterZSpread);
  tree_->SetBranchAddress("cscITClusterEtaPhiSpread",             cscITClusterEtaPhiSpread);
  tree_->SetBranchAddress("cscITClusterEtaSpread",             cscITClusterEtaSpread);
  tree_->SetBranchAddress("cscITClusterPhiSpread",             cscITClusterPhiSpread);
  tree_->SetBranchAddress("cscITClusterEta",             cscITClusterEta);
  tree_->SetBranchAddress("cscITClusterPhi",             cscITClusterPhi);
  tree_->SetBranchAddress("cscITClusterJetVeto",             cscITClusterJetVeto);
  tree_->SetBranchAddress("cscITClusterCaloJetVeto",             cscITClusterCaloJetVeto);
  tree_->SetBranchAddress("cscITClusterMuonVeto",             cscITClusterMuonVeto);
  tree_->SetBranchAddress("cscITClusterJetVetoE",             cscITClusterJetVetoE);
  tree_->SetBranchAddress("cscITClusterCaloJetVetoE",             cscITClusterCaloJetVetoE);
  tree_->SetBranchAddress("cscITClusterMuonVetoE",             cscITClusterMuonVetoE);
  tree_->SetBranchAddress("cscITClusterSize",             cscITClusterSize);


  tree_->SetBranchAddress("cscITCluster_match_cscCluster_index",             cscITCluster_match_cscCluster_index);
  tree_->SetBranchAddress("cscITCluster_cscCluster_SizeRatio",             cscITCluster_cscCluster_SizeRatio);
  */


  tree_->SetBranchAddress("gLLP_eta",    gLLP_eta);
  tree_->SetBranchAddress("gLLP_phi",    gLLP_phi);
  tree_->SetBranchAddress("gLLP_csc",    gLLP_csc);
  tree_->SetBranchAddress("gLLP_ctau",    gLLP_ctau);
  tree_->SetBranchAddress("gLLP_beta",    gLLP_beta);

  tree_->SetBranchAddress("gLLP_decay_vertex_r",    gLLP_decay_vertex_r);
  tree_->SetBranchAddress("gLLP_decay_vertex_x",    gLLP_decay_vertex_x);
  tree_->SetBranchAddress("gLLP_decay_vertex_y",    gLLP_decay_vertex_y);
  tree_->SetBranchAddress("gLLP_decay_vertex_z",    gLLP_decay_vertex_z);

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
  tree_->SetBranchAddress("lepPassVetoId", lepPassVetoId);
  tree_->SetBranchAddress("lepPassLooseIso", lepPassLooseIso);
  tree_->SetBranchAddress("lepPassTightIso", lepPassTightIso);
  tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);
  tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);


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
  tree_->SetBranchAddress("jet_match_genJet_pt", jet_match_genJet_pt);
  tree_->SetBranchAddress("jet_match_genJet_index", jet_match_genJet_index);
  tree_->SetBranchAddress("jet_match_genJet_minDeltaR", jet_match_genJet_minDeltaR);

  // tree_->SetBranchAddress("ecalNRechits",   ecalNRechits);*/
  // tree_->SetBranchAddress("ecalRechitE", ecalRechitE);
  tree_->SetBranchAddress("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction);
  tree_->SetBranchAddress("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction);
  tree_->SetBranchAddress("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction);
  tree_->SetBranchAddress("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction);

  // tree_->SetBranchAddress("jetLoosePassId", jetLoosePassId);
  tree_->SetBranchAddress("jetTightPassId", jetTightPassId);
  // triggers
  tree_->SetBranchAddress("HLTDecision",   HLTDecision);
};

void LiteTreeMuonSystem::LoadTree(const char* file)
{
  f_ = TFile::Open(file);
  assert(f_);
  tree_ = dynamic_cast<TTree*>(f_->Get("MuonSystem"));
  InitTree();
  assert(tree_);
};

void LiteTreeMuonSystem::CreateTree()
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
  tree_->Branch("higgsPtWeight",      &higgsPtWeight,     "higgsPtWeight/F");
  tree_->Branch("higgsPtWeightSys",      higgsPtWeightSys,     "higgsPtWeightSys[9]/F");
  tree_->Branch("sf_facScaleUp",      &sf_facScaleUp,     "sf_facScaleUp/F");
  tree_->Branch("sf_facScaleDown",      &sf_facScaleDown,     "sf_facScaleDown/F");
  tree_->Branch("sf_renScaleUp",      &sf_renScaleUp,     "sf_renScaleUp/F");
  tree_->Branch("sf_renScaleDown",      &sf_renScaleDown,     "sf_renScaleDown/F");
  tree_->Branch("sf_facRenScaleUp",      &sf_facRenScaleUp,     "sf_facRenScaleUp/F");
  tree_->Branch("sf_facRenScaleDown",      &sf_facRenScaleDown,     "sf_facRenScaleDown/F");

  tree_->Branch("pileupWeight",      &pileupWeight,     "pileupWeight/F");
  tree_->Branch("pileupWeightUp",      &pileupWeightUp,     "pileupWeightUp/F");
  tree_->Branch("pileupWeightDown",      &pileupWeightDown,     "pileupWeightDown/F");
  tree_->Branch("Flag_HBHENoiseFilter",      &Flag_HBHENoiseFilter,     "Flag_HBHENoiseFilter/O");
  tree_->Branch("Flag_BadPFMuonFilter",      &Flag_BadPFMuonFilter,     "Flag_BadPFMuonFilter/O");
  tree_->Branch("Flag_HBHEIsoNoiseFilter",      &Flag_HBHEIsoNoiseFilter,     "Flag_HBHEIsoNoiseFilter/O");
  tree_->Branch("Flag_CSCTightHaloFilter",      &Flag_CSCTightHaloFilter,     "Flag_CSCTightHaloFilter/O");
  tree_->Branch("Flag_globalSuperTightHalo2016Filter",      &Flag_globalSuperTightHalo2016Filter,     "Flag_globalSuperTightHalo2016Filter/O");
  tree_->Branch("Flag_goodVertices",      &Flag_goodVertices,     "Flag_goodVertices/O");
  tree_->Branch("Flag_ecalBadCalibFilter",      &Flag_ecalBadCalibFilter,     "Flag_ecalBadCalibFilter/O");
  tree_->Branch("Flag_BadChargedCandidateFilter",      &Flag_BadChargedCandidateFilter,     "Flag_BadChargedCandidateFilter/O");
  tree_->Branch("Flag_eeBadScFilter",      &Flag_eeBadScFilter,     "Flag_eeBadScFilter/O");
  tree_->Branch("Flag_all",      &Flag_all,     "Flag_all/O");

  tree_->Branch("Flag2_HBHENoiseFilter",      &Flag2_HBHENoiseFilter,     "Flag2_HBHENoiseFilter/O");
  tree_->Branch("Flag2_HBHEIsoNoiseFilter",      &Flag2_HBHEIsoNoiseFilter,     "Flag2_HBHEIsoNoiseFilter/O");
  tree_->Branch("Flag2_BadPFMuonFilter",      &Flag2_BadPFMuonFilter,     "Flag2_BadPFMuonFilter/O");
  tree_->Branch("Flag2_globalSuperTightHalo2016Filter",      &Flag2_globalSuperTightHalo2016Filter,     "Flag2_globalSuperTightHalo2016Filter/O");
  tree_->Branch("Flag2_globalTightHalo2016Filter",      &Flag2_globalTightHalo2016Filter,     "Flag2_globalTightHalo2016Filter/O");
  tree_->Branch("Flag2_BadChargedCandidateFilter",      &Flag2_BadChargedCandidateFilter,     "Flag2_BadChargedCandidateFilter/O");
  tree_->Branch("Flag2_EcalDeadCellTriggerPrimitiveFilter",      &Flag2_EcalDeadCellTriggerPrimitiveFilter,     "Flag2_EcalDeadCellTriggerPrimitiveFilter/O");
  tree_->Branch("Flag2_ecalBadCalibFilter",      &Flag2_ecalBadCalibFilter,     "Flag2_ecalBadCalibFilter/O");
  tree_->Branch("Flag2_eeBadScFilter",      &Flag2_eeBadScFilter,     "Flag2_eeBadScFilter/O");
  tree_->Branch("Flag2_all",      &Flag2_all,     "Flag2_all/O");



  // tree_->Branch("rho",         &rho,        "rho/F");
  tree_->Branch("met",         &met,        "met/F");         // MET
  tree_->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)
  tree_->Branch("metXYCorr",      &metXYCorr,     "metXYCorr/F");      // phi(MET)
  tree_->Branch("metPhiXYCorr",      &metPhiXYCorr,     "metPhiXYCorr/F");      // phi(MET)

  tree_->Branch("HT",      &HT,     "HT/F");      // phi(MET)

  tree_->Branch("jetMet_dPhi",      &jetMet_dPhi,     "jetMet_dPhi/F");      // phi(MET)
  tree_->Branch("jetMet_dPhiMin",      &jetMet_dPhiMin,     "jetMet_dPhiMin/F");      // phi(MET)
  tree_->Branch("jetMet_dPhiMin4",      &jetMet_dPhiMin4,     "jetMet_dPhiMin4/F");      // phi(MET)

  tree_->Branch("metJESUp",      &metJESUp,     "metJESUp/F");      // phi(MET)
  tree_->Branch("metJESDown",      &metJESDown,     "metJESDown/F");      // phi(MET)

  tree_->Branch("genMetPtTrue",         &genMetPtTrue,        "genMetPtTrue/F");         // MET
  tree_->Branch("genMetPhiTrue",      &genMetPhiTrue,     "genMetPhiTrue/F");      // phi(MET)
  tree_->Branch("genMetPtCalo",         &genMetPtCalo,        "genMetPtCalo/F");         // MET
  tree_->Branch("genMetPhiCalo",      &genMetPhiCalo,     "genMetPhiCalo/F");      // phi(MET)

  tree_->Branch("nGenParticle",      &nGenParticle,   "nGenParticle/I");
  tree_->Branch("gParticleId",      gParticleId,  "gParticleId[nGenParticle]/I");
  tree_->Branch("gParticleStatus",      gParticleStatus,  "gParticleStatus[nGenParticle]/I");
  tree_->Branch("gParticleMotherId",      gParticleMotherId,  "gParticleMotherId[nGenParticle]/I");
  tree_->Branch("gParticleE",      gParticleE,  "gParticleE[nGenParticle]/F");
  tree_->Branch("gParticlePt",      gParticlePt,  "gParticlePt[nGenParticle]/F");
  tree_->Branch("gParticleEta",      gParticleEta,  "gParticleEta[nGenParticle]/F");
  tree_->Branch("gParticlePhi",      gParticlePhi,  "gParticlePhi[nGenParticle]/F");

  tree_->Branch("nGenJets",      &nGenJets,  "nGenJets/I");
  tree_->Branch("genJetE",      genJetE,  "genJetE[nGenJets]/F");
  tree_->Branch("genJetPt",      genJetPt,  "genJetPt[nGenJets]/F");
  tree_->Branch("genJetEta",      genJetEta,  "genJetEta[nGenJets]/F");
  tree_->Branch("genJetPhi",      genJetPhi,  "genJetPhi[nGenJets]/F");
  tree_->Branch("genJetMET",      genJetMET,  "genJetMET[nGenJets]/F");


  // tree_->Branch("gLepId",      &gLepId,     "gLepId/I");      // phi(MET)
  // tree_->Branch("gLepPt",      &gLepPt,     "gLepPt/F");      // phi(MET)
  // tree_->Branch("gLepE",      &gLepE,     "gLepE/F");      // phi(MET)
  // tree_->Branch("gLepEta",      &gLepEta,     "gLepEta/F");      // phi(MET)
  // tree_->Branch("gLepPhi",      &gLepPhi,     "gLepPhi/F");      // phi(MET)
  tree_->Branch("gHiggsPt",      &gHiggsPt,     "gHiggsPt/F");      // phi(MET)
  tree_->Branch("gHiggsE",      &gHiggsE,     "gHiggsE/F");      // phi(MET)
  tree_->Branch("gHiggsEta",      &gHiggsEta,     "gHiggsEta/F");      // phi(MET)
  tree_->Branch("gHiggsPhi",      &gHiggsPhi,     "gHiggsPhi/F");      // phi(MET)
  //CSC

  tree_->Branch("nCscRechits",             &nCscRechits, "nCscRechits/I");
  tree_->Branch("nEarlyCscRechits",             &nEarlyCscRechits, "nEarlyCscRechits/I");
  tree_->Branch("nLateCscRechits",             &nLateCscRechits, "nLateCscRechits/I");
  tree_->Branch("nEarly2CscRechits",             &nEarly2CscRechits, "nEarly2CscRechits/I");
  tree_->Branch("nLate2CscRechits",             &nLate2CscRechits, "nLate2CscRechits/I");
  tree_->Branch("nCscRings",             &nCscRings, "nCscRings/I");
  tree_->Branch("nCscPositiveYRechits",             &nCscPositiveYRechits, "nCscPositiveYRechits/I");
  tree_->Branch("nCscNegativeYRechits",             &nCscNegativeYRechits, "nCscNegativeYRechits/I");
  tree_->Branch("cscPosTpeak",             &cscPosTpeak, "cscPosTpeak/F");
  tree_->Branch("cscNegTpeak",             &cscNegTpeak, "cscNegTpeak/F");

  tree_->Branch("nCscRechitsChamberPlus11",            &nCscRechitsChamberPlus11,             "nCscRechitsChamberPlus11/I");
  tree_->Branch("nCscRechitsChamberPlus12",            &nCscRechitsChamberPlus12,             "nCscRechitsChamberPlus12/I");
  tree_->Branch("nCscRechitsChamberPlus13",            &nCscRechitsChamberPlus13,             "nCscRechitsChamberPlus13/I");
  tree_->Branch("nCscRechitsChamberPlus21",            &nCscRechitsChamberPlus21,             "nCscRechitsChamberPlus21/I");
  tree_->Branch("nCscRechitsChamberPlus22",            &nCscRechitsChamberPlus22,             "nCscRechitsChamberPlus22/I");
  tree_->Branch("nCscRechitsChamberPlus31",            &nCscRechitsChamberPlus31,             "nCscRechitsChamberPlus31/I");
  tree_->Branch("nCscRechitsChamberPlus32",            &nCscRechitsChamberPlus32,             "nCscRechitsChamberPlus32/I");
  tree_->Branch("nCscRechitsChamberPlus41",            &nCscRechitsChamberPlus41,             "nCscRechitsChamberPlus41/I");
  tree_->Branch("nCscRechitsChamberPlus42",            &nCscRechitsChamberPlus42,             "nCscRechitsChamberPlus42/I");
  tree_->Branch("nCscRechitsChamberMinus11",            &nCscRechitsChamberMinus11,             "nCscRechitsChamberMinus11/I");
  tree_->Branch("nCscRechitsChamberMinus12",            &nCscRechitsChamberMinus12,             "nCscRechitsChamberMinus12/I");
  tree_->Branch("nCscRechitsChamberMinus13",            &nCscRechitsChamberMinus13,             "nCscRechitsChamberMinus13/I");
  tree_->Branch("nCscRechitsChamberMinus21",            &nCscRechitsChamberMinus21,             "nCscRechitsChamberMinus21/I");
  tree_->Branch("nCscRechitsChamberMinus22",            &nCscRechitsChamberMinus22,             "nCscRechitsChamberMinus22/I");
  tree_->Branch("nCscRechitsChamberMinus31",            &nCscRechitsChamberMinus31,             "nCscRechitsChamberMinus31/I");
  tree_->Branch("nCscRechitsChamberMinus32",            &nCscRechitsChamberMinus32,             "nCscRechitsChamberMinus32/I");
  tree_->Branch("nCscRechitsChamberMinus41",            &nCscRechitsChamberMinus41,             "nCscRechitsChamberMinus41/I");
  tree_->Branch("nCscRechitsChamberMinus42",            &nCscRechitsChamberMinus42,             "nCscRechitsChamberMinus42/I");


  tree_->Branch("nDTRechits",            &nDTRechits,             "nDTRechits/I");
  tree_->Branch("nDtRings",             &nDtRings, "nDtRings/I");
  tree_->Branch("nDTPositiveYRechits",             &nDTPositiveYRechits, "nDTPositiveYRechits/I");
  tree_->Branch("nDTNegativeYRechits",             &nDTNegativeYRechits, "nDTNegativeYRechits/I");

  tree_->Branch("nDTRechitsChamberMinus12",            &nDTRechitsChamberMinus12,             "nDTRechitsChamberMinus12/I");
  tree_->Branch("nDTRechitsChamberMinus11",            &nDTRechitsChamberMinus11,             "nDTRechitsChamberMinus11/I");
  tree_->Branch("nDTRechitsChamber10",            &nDTRechitsChamber10,             "nDTRechitsChamber10/I");
  tree_->Branch("nDTRechitsChamberPlus11",            &nDTRechitsChamberPlus11,             "nDTRechitsChamberPlus11/I");
  tree_->Branch("nDTRechitsChamberPlus12",            &nDTRechitsChamberPlus12,             "nDTRechitsChamberPlus12/I");
  tree_->Branch("nDTRechitsChamberMinus22",            &nDTRechitsChamberMinus22,             "nDTRechitsChamberMinus22/I");
  tree_->Branch("nDTRechitsChamberMinus21",            &nDTRechitsChamberMinus21,             "nDTRechitsChamberMinus21/I");
  tree_->Branch("nDTRechitsChamber20",            &nDTRechitsChamber20,             "nDTRechitsChamber20/I");
  tree_->Branch("nDTRechitsChamberPlus21",            &nDTRechitsChamberPlus21,             "nDTRechitsChamberPlus21/I");
  tree_->Branch("nDTRechitsChamberPlus22",            &nDTRechitsChamberPlus22,             "nDTRechitsChamberPlus22/I");
  tree_->Branch("nDTRechitsChamberMinus32",            &nDTRechitsChamberMinus32,             "nDTRechitsChamberMinus32/I");
  tree_->Branch("nDTRechitsChamberMinus31",            &nDTRechitsChamberMinus31,             "nDTRechitsChamberMinus31/I");
  tree_->Branch("nDTRechitsChamber30",            &nDTRechitsChamber30,             "nDTRechitsChamber30/I");

  tree_->Branch("nDTRechitsChamberPlus31",            &nDTRechitsChamberPlus31,             "nDTRechitsChamberPlus31/I");
  tree_->Branch("nDTRechitsChamberPlus32",            &nDTRechitsChamberPlus32,             "nDTRechitsChamberPlus32/I");
  tree_->Branch("nDTRechitsChamberMinus42",            &nDTRechitsChamberMinus42,             "nDTRechitsChamberMinus42/I");
  tree_->Branch("nDTRechitsChamberMinus41",            &nDTRechitsChamberMinus41,             "nDTRechitsChamberMinus41/I");
  tree_->Branch("nDTRechitsChamber40",            &nDTRechitsChamber40,             "nDTRechitsChamber40/I");
  tree_->Branch("nDTRechitsChamberPlus41",            &nDTRechitsChamberPlus41,             "nDTRechitsChamberPlus41/I");
  tree_->Branch("nDTRechitsChamberPlus42",            &nDTRechitsChamberPlus42,             "nDTRechitsChamberPlus42/I");
  tree_->Branch("cscRechitsStation",             cscRechitsStation,             "cscRechitsStation[nCscRechits]/I");
  tree_->Branch("cscRechitsChamber",             cscRechitsChamber,             "cscRechitsChamber[nCscRechits]/I");

  tree_->Branch("cscRechitsPhi",           cscRechitsPhi,           "cscRechitsPhi[nCscRechits]/F");
  tree_->Branch("cscRechitsEta",           cscRechitsEta,           "cscRechitsEta[nCscRechits]/F");
  tree_->Branch("cscRechitsQuality",             cscRechitsQuality,             "cscRechitsQuality[nCscRechits]/I");

  tree_->Branch("cscRechitsX",             cscRechitsX,             "cscRechitsX[nCscRechits]/F");

  tree_->Branch("cscRechitsY",             cscRechitsY,             "cscRechitsY[nCscRechits]/F");
  tree_->Branch("cscRechitsZ",             cscRechitsZ,             "cscRechitsZ[nCscRechits]/F");
  tree_->Branch("cscRechitsTwire",             cscRechitsTwire,             "cscRechitsTwire[nCscRechits]/F");
  tree_->Branch("cscRechitsTpeak",             cscRechitsTpeak,             "cscRechitsTpeak[nCscRechits]/F");
  tree_->Branch("cscRechitsClusterId",             cscRechitsClusterId,             "cscRechitsClusterId[nCscRechits]/I");
  tree_->Branch("cscRechitsCluster2Id",             cscRechitsCluster2Id,             "cscRechitsCluster2Id[nCscRechits]/I");


/*  tree_->Branch("nCscClusters",             &nCscClusters, "nCscClusters/I");
  tree_->Branch("cscCluster_match_gLLP",             cscCluster_match_gLLP, "cscCluster_match_gLLP[nCscClusters]/O");
  tree_->Branch("cscCluster_match_gLLP_index",             cscCluster_match_gLLP_index, "cscCluster_match_gLLP_index[nCscClusters]/I");
  tree_->Branch("cscCluster_match_gLLP_minDeltaR",             cscCluster_match_gLLP_minDeltaR, "cscCluster_match_gLLP_minDeltaR[nCscClusters]/I");

  // tree_->Branch("nCsc_JetVetoCluster0p4",             &nCsc_JetVetoCluster0p4, "nCsc_JetVetoCluster0p4/I");
  // tree_->Branch("nCsc_JetMuonVetoCluster0p4",             &nCsc_JetMuonVetoCluster0p4, "nCsc_JetMuonVetoCluster0p4/I");
  // tree_->Branch("nCsc_JetVetoCluster0p4_Me1112Veto",             &nCsc_JetVetoCluster0p4_Me1112Veto, "nCsc_JetVetoCluster0p4_Me1112Veto/I");
  // tree_->Branch("nCsc_JetMuonVetoCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoCluster0p4_Me1112Veto, "nCsc_JetMuonVetoCluster0p4_Me1112Veto/I");
  tree_->Branch("cscClusterX",             cscClusterX,             "cscClusterX[nCscClusters]/F");
  tree_->Branch("cscClusterY",             cscClusterY,             "cscClusterY[nCscClusters]/F");
  tree_->Branch("cscClusterZ",             cscClusterZ,             "cscClusterZ[nCscClusters]/F");
  tree_->Branch("cscClusterTime",             cscClusterTime,             "cscClusterTime[nCscClusters]/F");
  tree_->Branch("cscClusterTimeSpread",             cscClusterTimeSpread,             "cscClusterTimeSpread[nCscClusters]/F");
  tree_->Branch("cscClusterGenMuonDeltaR",             cscClusterGenMuonDeltaR,             "cscClusterGenMuonDeltaR[nCscClusters]/F");

  tree_->Branch("cscClusterMajorAxis",             cscClusterMajorAxis,             "cscClusterMajorAxis[nCscClusters]/F");
  tree_->Branch("cscClusterMinorAxis",             cscClusterMinorAxis,             "cscClusterMinorAxis[nCscClusters]/F");
  tree_->Branch("cscClusterEtaPhiSpread",             cscClusterEtaPhiSpread,             "cscClusterEtaPhiSpread[nCscClusters]/F");
  tree_->Branch("cscClusterPhiSpread",             cscClusterPhiSpread,             "cscClusterPhiSpread[nCscClusters]/F");
  tree_->Branch("cscClusterEtaSpread",             cscClusterEtaSpread,             "cscClusterEtaSpread[nCscClusters]/F");
  tree_->Branch("cscClusterXSpread",             cscClusterXSpread,             "cscClusterXSpread[nCscClusters]/F");
  tree_->Branch("cscClusterYSpread",             cscClusterYSpread,             "cscClusterYSpread[nCscClusters]/F");
  tree_->Branch("cscClusterZSpread",             cscClusterZSpread,             "cscClusterZSpread[nCscClusters]/F");
  tree_->Branch("cscClusterPhi",             cscClusterPhi,             "cscClusterPhi[nCscClusters]/F");
  tree_->Branch("cscClusterEta",             cscClusterEta,             "cscClusterEta[nCscClusters]/F");
  tree_->Branch("cscClusterJetVetoPt",             cscClusterJetVetoPt,             "cscClusterJetVetoPt[nCscClusters]/F");
  tree_->Branch("cscClusterMuonVetoPt",             cscClusterMuonVetoPt,             "cscClusterMuonVetoPt[nCscClusters]/F");
  // tree_->Branch("cscClusterCaloJetVeto",             cscClusterCaloJetVeto,             "cscClusterCaloJetVeto[nCscClusters]/F");
  tree_->Branch("cscClusterJetVetoE",             cscClusterJetVetoE,             "cscClusterJetVetoE[nCscClusters]/F");
  tree_->Branch("cscClusterMuonVetoE",             cscClusterMuonVetoE,             "cscClusterMuonVetoE[nCscClusters]/F");
  tree_->Branch("cscClusterSize",             cscClusterSize,             "cscClusterSize[nCscClusters]/I");
  tree_->Branch("cscClusterMe11Ratio",             cscClusterMe11Ratio,             "cscClusterMe11Ratio[nCscClusters]/F");
  tree_->Branch("cscClusterMe12Ratio",             cscClusterMe12Ratio,             "cscClusterMe12Ratio[nCscClusters]/F");
  tree_->Branch("cscClusterNStation",             cscClusterNStation,             "cscClusterNStation[nCscClusters]/I");
  tree_->Branch("cscClusterMaxStation",             cscClusterMaxStation,             "cscClusterMaxStation[nCscClusters]/I");
  tree_->Branch("cscClusterMaxStationRatio",             cscClusterMaxStationRatio,             "cscClusterMaxStationRatio[nCscClusters]/F");
  tree_->Branch("cscClusterNChamber",             cscClusterNChamber,             "cscClusterNChamber[nCscClusters]/I");
  tree_->Branch("cscClusterMaxChamber",             cscClusterMaxChamber,             "cscClusterMaxChamber[nCscClusters]/I");
  tree_->Branch("cscClusterMaxChamberRatio",             cscClusterMaxChamberRatio,             "cscClusterMaxChamberRatio[nCscClusters]/F");
  tree_->Branch("cscClusterNSegmentChamberPlus11",             cscClusterNSegmentChamberPlus11,             "cscClusterNSegmentChamberPlus11[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberPlus12",             cscClusterNSegmentChamberPlus12,             "cscClusterNSegmentChamberPlus12[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberPlus13",             cscClusterNSegmentChamberPlus13,             "cscClusterNSegmentChamberPlus13[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberPlus21",             cscClusterNSegmentChamberPlus21,             "cscClusterNSegmentChamberPlus21[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberPlus22",             cscClusterNSegmentChamberPlus22,             "cscClusterNSegmentChamberPlus22[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberPlus31",             cscClusterNSegmentChamberPlus31,             "cscClusterNSegmentChamberPlus31[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberPlus32",             cscClusterNSegmentChamberPlus32,             "cscClusterNSegmentChamberPlus32[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberPlus41",             cscClusterNSegmentChamberPlus41,             "cscClusterNSegmentChamberPlus41[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberPlus42",             cscClusterNSegmentChamberPlus42,             "cscClusterNSegmentChamberPlus42[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberMinus11",             cscClusterNSegmentChamberMinus11,             "cscClusterNSegmentChamberMinus11[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberMinus12",             cscClusterNSegmentChamberMinus12,             "cscClusterNSegmentChamberMinus12[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberMinus13",             cscClusterNSegmentChamberMinus13,             "cscClusterNSegmentChamberMinus13[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberMinus21",             cscClusterNSegmentChamberMinus21,             "cscClusterNSegmentChamberMinus21[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberMinus22",             cscClusterNSegmentChamberMinus22,             "cscClusterNSegmentChamberMinus22[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberMinus31",             cscClusterNSegmentChamberMinus31,             "cscClusterNSegmentChamberMinus31[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberMinus32",             cscClusterNSegmentChamberMinus32,             "cscClusterNSegmentChamberMinus32[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberMinus41",             cscClusterNSegmentChamberMinus41,             "cscClusterNSegmentChamberMinus41[nCscClusters]/I");
  tree_->Branch("cscClusterNSegmentChamberMinus42",             cscClusterNSegmentChamberMinus42,             "cscClusterNSegmentChamberMinus42[nCscClusters]/I");
*/
  // all csc SegClusters
  /*
  tree_->Branch("nCscSegClusters",             &nCscSegClusters, "nCscSegClusters/I");
  // tree_->Branch("nCsc_JetVetoSegCluster0p4",             &nCsc_JetVetoSegCluster0p4, "nCsc_JetVetoSegCluster0p4/I");
  // tree_->Branch("nCsc_JetMuonVetoSegCluster0p4",             &nCsc_JetMuonVetoSegCluster0p4, "nCsc_JetMuonVetoSegCluster0p4/I");
  // tree_->Branch("nCsc_JetVetoSegCluster0p4_Me1112Veto",             &nCsc_JetVetoSegCluster0p4_Me1112Veto, "nCsc_JetVetoSegCluster0p4_Me1112Veto/I");
  // tree_->Branch("nCsc_JetMuonVetoSegCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoSegCluster0p4_Me1112Veto, "nCsc_JetMuonVetoSegCluster0p4_Me1112Veto/I");
  tree_->Branch("cscSegClusterX",             cscSegClusterX,             "cscSegClusterX[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterY",             cscSegClusterY,             "cscSegClusterY[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterZ",             cscSegClusterZ,             "cscSegClusterZ[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterTime",             cscSegClusterTime,             "cscSegClusterTime[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterTimeSpread",             cscSegClusterTimeSpread,             "cscSegClusterTimeSpread[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterGenMuonDeltaR",             cscSegClusterGenMuonDeltaR,             "cscSegClusterGenMuonDeltaR[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterMet_dPhi",             cscSegClusterMet_dPhi,             "cscSegClusterMet_dPhi[nCscSegClusters]/F");

  tree_->Branch("cscSegClusterMajorAxis",             cscSegClusterMajorAxis,             "cscSegClusterMajorAxis[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterMinorAxis",             cscSegClusterMinorAxis,             "cscSegClusterMinorAxis[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterEtaPhiSpread",             cscSegClusterEtaPhiSpread,             "cscSegClusterEtaPhiSpread[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterPhiSpread",             cscSegClusterPhiSpread,             "cscSegClusterPhiSpread[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterEtaSpread",             cscSegClusterEtaSpread,             "cscSegClusterEtaSpread[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterXSpread",             cscSegClusterXSpread,             "cscSegClusterXSpread[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterYSpread",             cscSegClusterYSpread,             "cscSegClusterYSpread[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterZSpread",             cscSegClusterZSpread,             "cscSegClusterZSpread[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterPhi",             cscSegClusterPhi,             "cscSegClusterPhi[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterEta",             cscSegClusterEta,             "cscSegClusterEta[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterJetVetoPt",             cscSegClusterJetVetoPt,             "cscSegClusterJetVetoPt[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterMuonVetoPt",             cscSegClusterMuonVetoPt,             "cscSegClusterMuonVetoPt[nCscSegClusters]/F");
  // tree_->Branch("cscSegClusterCaloJetVeto",             cscSegClusterCaloJetVeto,             "cscSegClusterCaloJetVeto[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterJetVetoE",             cscSegClusterJetVetoE,             "cscSegClusterJetVetoE[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterMuonVetoE",             cscSegClusterMuonVetoE,             "cscSegClusterMuonVetoE[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterSize",             cscSegClusterSize,             "cscSegClusterSize[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterMe11Ratio",             cscSegClusterMe11Ratio,             "cscSegClusterMe11Ratio[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterMe12Ratio",             cscSegClusterMe12Ratio,             "cscSegClusterMe12Ratio[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterNStation",             cscSegClusterNStation,             "cscSegClusterNStation[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterMaxStation",             cscSegClusterMaxStation,             "cscSegClusterMaxStation[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterMaxStationRatio",             cscSegClusterMaxStationRatio,             "cscSegClusterMaxStationRatio[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterNChamber",             cscSegClusterNChamber,             "cscSegClusterNChamber[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterMaxChamber",             cscSegClusterMaxChamber,             "cscSegClusterMaxChamber[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterMaxChamberRatio",             cscSegClusterMaxChamberRatio,             "cscSegClusterMaxChamberRatio[nCscSegClusters]/F");
  tree_->Branch("cscSegClusterNSegmentChamberPlus11",             cscSegClusterNSegmentChamberPlus11,             "cscSegClusterNSegmentChamberPlus11[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberPlus12",             cscSegClusterNSegmentChamberPlus12,             "cscSegClusterNSegmentChamberPlus12[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberPlus13",             cscSegClusterNSegmentChamberPlus13,             "cscSegClusterNSegmentChamberPlus13[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberPlus21",             cscSegClusterNSegmentChamberPlus21,             "cscSegClusterNSegmentChamberPlus21[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberPlus22",             cscSegClusterNSegmentChamberPlus22,             "cscSegClusterNSegmentChamberPlus22[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberPlus31",             cscSegClusterNSegmentChamberPlus31,             "cscSegClusterNSegmentChamberPlus31[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberPlus32",             cscSegClusterNSegmentChamberPlus32,             "cscSegClusterNSegmentChamberPlus32[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberPlus41",             cscSegClusterNSegmentChamberPlus41,             "cscSegClusterNSegmentChamberPlus41[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberPlus42",             cscSegClusterNSegmentChamberPlus42,             "cscSegClusterNSegmentChamberPlus42[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberMinus11",             cscSegClusterNSegmentChamberMinus11,             "cscSegClusterNSegmentChamberMinus11[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberMinus12",             cscSegClusterNSegmentChamberMinus12,             "cscSegClusterNSegmentChamberMinus12[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberMinus13",             cscSegClusterNSegmentChamberMinus13,             "cscSegClusterNSegmentChamberMinus13[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberMinus21",             cscSegClusterNSegmentChamberMinus21,             "cscSegClusterNSegmentChamberMinus21[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberMinus22",             cscSegClusterNSegmentChamberMinus22,             "cscSegClusterNSegmentChamberMinus22[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberMinus31",             cscSegClusterNSegmentChamberMinus31,             "cscSegClusterNSegmentChamberMinus31[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberMinus32",             cscSegClusterNSegmentChamberMinus32,             "cscSegClusterNSegmentChamberMinus32[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberMinus41",             cscSegClusterNSegmentChamberMinus41,             "cscSegClusterNSegmentChamberMinus41[nCscSegClusters]/I");
  tree_->Branch("cscSegClusterNSegmentChamberMinus42",             cscSegClusterNSegmentChamberMinus42,             "cscSegClusterNSegmentChamberMinus42[nCscSegClusters]/I");
  tree_->Branch("cscSegCluster_match_gParticle_id",             cscSegCluster_match_gParticle_id,             "cscSegCluster_match_gParticle_id[nCscSegClusters]/I");
  tree_->Branch("cscSegCluster_match_gParticle_minDeltaR",             cscSegCluster_match_gParticle_minDeltaR,             "cscSegCluster_match_gParticle_minDeltaR[nCscSegClusters]/F");
  tree_->Branch("cscSegCluster_match_gParticle_index",             cscSegCluster_match_gParticle_index,             "cscSegCluster_match_gParticle_index[nCscSegClusters]/I");
  */

    // all csc RechitClusters

    tree_->Branch("nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto, "nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto/I");
    tree_->Branch("nCscRechitClusters",             &nCscRechitClusters, "nCscRechitClusters/I");
    tree_->Branch("cscRechitCluster_match_gLLP",             cscRechitCluster_match_gLLP,             "cscRechitCluster_match_gLLP[nCscRechitClusters]/O");
    tree_->Branch("cscRechitCluster_match_gLLP_minDeltaR",             cscRechitCluster_match_gLLP_minDeltaR,             "cscRechitCluster_match_gLLP_minDeltaR[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_index",             cscRechitCluster_match_gLLP_index,             "cscRechitCluster_match_gLLP_index[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_gParticle_id",             cscRechitCluster_match_gParticle_id,             "cscRechitCluster_match_gParticle_id[nCscRechitClusters]/O");
    tree_->Branch("cscRechitCluster_match_gLLP_eta",             cscRechitCluster_match_gLLP_eta, "cscRechitCluster_match_gLLP_eta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_phi",             cscRechitCluster_match_gLLP_phi, "cscRechitCluster_match_gLLP_phi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_decay_r",             cscRechitCluster_match_gLLP_decay_r, "cscRechitCluster_match_gLLP_decay_r[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_decay_x",             cscRechitCluster_match_gLLP_decay_x, "cscRechitCluster_match_gLLP_decay_x[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_decay_y",             cscRechitCluster_match_gLLP_decay_y, "cscRechitCluster_match_gLLP_decay_y[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_decay_z",             cscRechitCluster_match_gLLP_decay_z, "cscRechitCluster_match_gLLP_decay_z[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_ctau",             cscRechitCluster_match_gLLP_ctau, "cscRechitCluster_match_gLLP_ctau[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_beta",             cscRechitCluster_match_gLLP_beta, "cscRechitCluster_match_gLLP_beta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_csc",             cscRechitCluster_match_gLLP_csc, "cscRechitCluster_match_gLLP_csc[nCscRechitClusters]/O");

    tree_->Branch("cscRechitClusterX",             cscRechitClusterX,             "cscRechitClusterX[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterY",             cscRechitClusterY,             "cscRechitClusterY[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterZ",             cscRechitClusterZ,             "cscRechitClusterZ[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTime",             cscRechitClusterTime,             "cscRechitClusterTime[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeTotal",             cscRechitClusterTimeTotal,             "cscRechitClusterTimeTotal[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterTimeSpread",             cscRechitClusterTimeSpread,             "cscRechitClusterTimeSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterGenMuonDeltaR",             cscRechitClusterGenMuonDeltaR,             "cscRechitClusterGenMuonDeltaR[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterMajorAxis",             cscRechitClusterMajorAxis,             "cscRechitClusterMajorAxis[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMinorAxis",             cscRechitClusterMinorAxis,             "cscRechitClusterMinorAxis[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEtaPhiSpread",             cscRechitClusterEtaPhiSpread,             "cscRechitClusterEtaPhiSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterPhiSpread",             cscRechitClusterPhiSpread,             "cscRechitClusterPhiSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEtaSpread",             cscRechitClusterEtaSpread,             "cscRechitClusterEtaSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterXYSpread",             cscRechitClusterXYSpread,             "cscRechitClusterXYSpread[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterXSpread",             cscRechitClusterXSpread,             "cscRechitClusterXSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterYSpread",             cscRechitClusterYSpread,             "cscRechitClusterYSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterZSpread",             cscRechitClusterZSpread,             "cscRechitClusterZSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterPhi",             cscRechitClusterPhi,             "cscRechitClusterPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEta",             cscRechitClusterEta,             "cscRechitClusterEta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoPt",             cscRechitClusterJetVetoPt,             "cscRechitClusterJetVetoPt[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoPt",             cscRechitClusterMuonVetoPt,             "cscRechitClusterMuonVetoPt[nCscRechitClusters]/F");
    // tree_->Branch("cscRechitClusterCaloJetVeto",             cscRechitClusterCaloJetVeto,             "cscRechitClusterCaloJetVeto[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterXYSpread_corr",             cscRechitClusterXYSpread_corr, "cscRechitClusterXYSpread_corr[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterXSpread_corr",             cscRechitClusterXSpread_corr, "cscRechitClusterXSpread_corr[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterYSpread_corr",            cscRechitClusterYSpread_corr, "cscRechitClusterYSpread_corr[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterPhiSpread_corr",         cscRechitClusterPhiSpread_corr, "cscRechitClusterPhiSpread_corr[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEtaSpread_corr",         cscRechitClusterEtaSpread_corr, "cscRechitClusterEtaSpread_corr[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEtaPhiSpread_corr",      cscRechitClusterEtaPhiSpread_corr, "cscRechitClusterEtaPhiSpread_corr[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterRSpread_corr",           cscRechitClusterRSpread_corr, "cscRechitClusterRSpread_corr[nCscRechitClusters]/F");


    tree_->Branch("cscRechitCluster_match_MB1Seg_0p4",             cscRechitCluster_match_MB1Seg_0p4, "cscRechitCluster_match_MB1Seg_0p4[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_RB1_0p4",             cscRechitCluster_match_RB1_0p4, "cscRechitCluster_match_RB1_0p4[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_RE12_0p4",             cscRechitCluster_match_RE12_0p4, "cscRechitCluster_match_RE12_0p4[nCscRechitClusters]/I");

    tree_->Branch("cscRechitClusterJetVetoE",             cscRechitClusterJetVetoE,             "cscRechitClusterJetVetoE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoE",             cscRechitClusterMuonVetoE,             "cscRechitClusterMuonVetoE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterGenMuonVetoE",             cscRechitClusterGenMuonVetoE,             "cscRechitClusterGenMuonVetoE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterGenMuonVetoPt",             cscRechitClusterGenMuonVetoPt,             "cscRechitClusterGenMuonVetoPt[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterSize",             cscRechitClusterSize,             "cscRechitClusterSize[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterMe11Ratio",             cscRechitClusterMe11Ratio,             "cscRechitClusterMe11Ratio[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMe12Ratio",             cscRechitClusterMe12Ratio,             "cscRechitClusterMe12Ratio[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterNStation",             cscRechitClusterNStation,             "cscRechitClusterNStation[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNStation5",             cscRechitClusterNStation5,             "cscRechitClusterNStation5[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNStation10perc",             cscRechitClusterNStation10perc,             "cscRechitClusterNStation10perc[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterAvgStation",             cscRechitClusterAvgStation,             "cscRechitClusterAvgStation[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterAvgStation5",             cscRechitClusterAvgStation5,             "cscRechitClusterAvgStation5[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterAvgStation10perc",             cscRechitClusterAvgStation10perc,             "cscRechitClusterAvgStation10perc[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterMaxStation",             cscRechitClusterMaxStation,             "cscRechitClusterMaxStation[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterMaxStationRatio",             cscRechitClusterMaxStationRatio,             "cscRechitClusterMaxStationRatio[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterNChamber",             cscRechitClusterNChamber,             "cscRechitClusterNChamber[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterMaxChamber",             cscRechitClusterMaxChamber,             "cscRechitClusterMaxChamber[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterMaxChamberRatio",             cscRechitClusterMaxChamberRatio,             "cscRechitClusterMaxChamberRatio[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterNRechitChamberPlus11",             cscRechitClusterNRechitChamberPlus11,             "cscRechitClusterNRechitChamberPlus11[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus12",             cscRechitClusterNRechitChamberPlus12,             "cscRechitClusterNRechitChamberPlus12[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus13",             cscRechitClusterNRechitChamberPlus13,             "cscRechitClusterNRechitChamberPlus13[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus21",             cscRechitClusterNRechitChamberPlus21,             "cscRechitClusterNRechitChamberPlus21[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus22",             cscRechitClusterNRechitChamberPlus22,             "cscRechitClusterNRechitChamberPlus22[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus31",             cscRechitClusterNRechitChamberPlus31,             "cscRechitClusterNRechitChamberPlus31[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus32",             cscRechitClusterNRechitChamberPlus32,             "cscRechitClusterNRechitChamberPlus32[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus41",             cscRechitClusterNRechitChamberPlus41,             "cscRechitClusterNRechitChamberPlus41[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberPlus42",             cscRechitClusterNRechitChamberPlus42,             "cscRechitClusterNRechitChamberPlus42[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus11",             cscRechitClusterNRechitChamberMinus11,             "cscRechitClusterNRechitChamberMinus11[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus12",             cscRechitClusterNRechitChamberMinus12,             "cscRechitClusterNRechitChamberMinus12[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus13",             cscRechitClusterNRechitChamberMinus13,             "cscRechitClusterNRechitChamberMinus13[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus21",             cscRechitClusterNRechitChamberMinus21,             "cscRechitClusterNRechitChamberMinus21[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus22",             cscRechitClusterNRechitChamberMinus22,             "cscRechitClusterNRechitChamberMinus22[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus31",             cscRechitClusterNRechitChamberMinus31,             "cscRechitClusterNRechitChamberMinus31[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus32",             cscRechitClusterNRechitChamberMinus32,             "cscRechitClusterNRechitChamberMinus32[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus41",             cscRechitClusterNRechitChamberMinus41,             "cscRechitClusterNRechitChamberMinus41[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNRechitChamberMinus42",             cscRechitClusterNRechitChamberMinus42,             "cscRechitClusterNRechitChamberMinus42[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterMet_dPhi",             cscRechitClusterMet_dPhi,             "cscRechitClusterMet_dPhi[nCscRechitClusters]/F");

    tree_->Branch("nCscRechitClusters2",             &nCscRechitClusters2, "nCscRechitClusters2/I");
    tree_->Branch("cscRechitCluster2_match_Me1112_0p4",             cscRechitCluster2_match_Me1112_0p4,             "cscRechitCluster2_match_Me1112_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_Me1112_0p6",             cscRechitCluster2_match_Me1112_0p6,             "cscRechitCluster2_match_Me1112_0p6[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_Me1112_0p8",             cscRechitCluster2_match_Me1112_0p8,             "cscRechitCluster2_match_Me1112_0p8[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_Me11_0p4",             cscRechitCluster2_match_Me11_0p4,             "cscRechitCluster2_match_Me11_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_Me11_0p6",             cscRechitCluster2_match_Me11_0p6,             "cscRechitCluster2_match_Me11_0p6[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_Me11_0p8",             cscRechitCluster2_match_Me11_0p8,             "cscRechitCluster2_match_Me11_0p8[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_Me12_0p4",             cscRechitCluster2_match_Me12_0p4,             "cscRechitCluster2_match_Me12_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_Me12_0p6",             cscRechitCluster2_match_Me12_0p6,             "cscRechitCluster2_match_Me12_0p6[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_Me12_0p8",             cscRechitCluster2_match_Me12_0p8,             "cscRechitCluster2_match_Me12_0p8[nCscRechitClusters2]/I");

    tree_->Branch("cscRechitCluster2_match_gParticle_id",             cscRechitCluster2_match_gParticle_id,             "cscRechitCluster2_match_gParticle_id[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_gParticle",             cscRechitCluster2_match_gParticle,             "cscRechitCluster2_match_gParticle[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2_match_gParticle_minDeltaR",             cscRechitCluster2_match_gParticle_minDeltaR,             "cscRechitCluster2_match_gParticle_minDeltaR[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gParticle_index",             cscRechitCluster2_match_gParticle_index,             "cscRechitCluster2_match_gParticle_index[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_gParticle_eta",             cscRechitCluster2_match_gParticle_eta,             "cscRechitCluster2_match_gParticle_eta[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gParticle_phi",             cscRechitCluster2_match_gParticle_phi,             "cscRechitCluster2_match_gParticle_phi[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gParticle_E",             cscRechitCluster2_match_gParticle_E,             "cscRechitCluster2_match_gParticle_E[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gParticle_pt",             cscRechitCluster2_match_gParticle_pt,             "cscRechitCluster2_match_gParticle_pt[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gParticle_MotherId",             cscRechitCluster2_match_gParticle_MotherId,             "cscRechitCluster2_match_gParticle_MotherId[nCscRechitClusters2]/I");

    tree_->Branch("cscRechitCluster2_match_gLLP",             cscRechitCluster2_match_gLLP,             "cscRechitCluster2_match_gLLP[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2_match_gLLP_minDeltaR",             cscRechitCluster2_match_gLLP_minDeltaR,             "cscRechitCluster2_match_gLLP_minDeltaR[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gLLP_index",             cscRechitCluster2_match_gLLP_index,             "cscRechitCluster2_match_gLLP_index[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_gLLP_eta",             cscRechitCluster2_match_gLLP_eta, "cscRechitCluster2_match_gLLP_eta[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gLLP_phi",             cscRechitCluster2_match_gLLP_phi, "cscRechitCluster2_match_gLLP_phi[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gLLP_decay_r",             cscRechitCluster2_match_gLLP_decay_r, "cscRechitCluster2_match_gLLP_decay_r[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gLLP_decay_x",             cscRechitCluster2_match_gLLP_decay_x, "cscRechitCluster2_match_gLLP_decay_x[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gLLP_decay_y",             cscRechitCluster2_match_gLLP_decay_y, "cscRechitCluster2_match_gLLP_decay_y[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gLLP_decay_z",             cscRechitCluster2_match_gLLP_decay_z, "cscRechitCluster2_match_gLLP_decay_z[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gLLP_ctau",             cscRechitCluster2_match_gLLP_ctau, "cscRechitCluster2_match_gLLP_ctau[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gLLP_beta",             cscRechitCluster2_match_gLLP_beta, "cscRechitCluster2_match_gLLP_beta[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_gLLP_csc",             cscRechitCluster2_match_gLLP_csc, "cscRechitCluster2_match_gLLP_csc[nCscRechitClusters2]/O");

    tree_->Branch("cscRechitCluster2X",             cscRechitCluster2X,             "cscRechitCluster2X[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2Y",             cscRechitCluster2Y,             "cscRechitCluster2Y[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2Z",             cscRechitCluster2Z,             "cscRechitCluster2Z[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2Time",             cscRechitCluster2Time,             "cscRechitCluster2Time[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2TimeTotal",             cscRechitCluster2TimeTotal,             "cscRechitCluster2TimeTotal[nCscRechitClusters2]/F");

    tree_->Branch("cscRechitCluster2TimeSpread",             cscRechitCluster2TimeSpread,             "cscRechitCluster2TimeSpread[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2GenMuonDeltaR",             cscRechitCluster2GenMuonDeltaR,             "cscRechitCluster2GenMuonDeltaR[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2XYSpread",             cscRechitCluster2XYSpread,             "cscRechitCluster2XYSpread[nCscRechitClusters2]/F");

    tree_->Branch("cscRechitCluster2MajorAxis",             cscRechitCluster2MajorAxis,             "cscRechitCluster2MajorAxis[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2MinorAxis",             cscRechitCluster2MinorAxis,             "cscRechitCluster2MinorAxis[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2EtaPhiSpread",             cscRechitCluster2EtaPhiSpread,             "cscRechitCluster2EtaPhiSpread[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2PhiSpread",             cscRechitCluster2PhiSpread,             "cscRechitCluster2PhiSpread[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2EtaSpread",             cscRechitCluster2EtaSpread,             "cscRechitCluster2EtaSpread[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2XSpread",             cscRechitCluster2XSpread,             "cscRechitCluster2XSpread[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2RSpread",             cscRechitCluster2RSpread,             "cscRechitCluster2RSpread[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2YSpread",             cscRechitCluster2YSpread,             "cscRechitCluster2YSpread[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2ZSpread",             cscRechitCluster2ZSpread,             "cscRechitCluster2ZSpread[nCscRechitClusters2]/F");

    tree_->Branch("cscRechitCluster2XYSpread_corr",             &cscRechitCluster2XYSpread_corr,  Form("cscRechitCluster2XYSpread_corr[nCscRechitClusters2][%d][%d]/F",N_phicorr, N_rcorr));
    tree_->Branch("cscRechitCluster2XSpread_corr",             &cscRechitCluster2XSpread_corr, Form("cscRechitCluster2XSpread_corr[nCscRechitClusters2][%d][%d]/F",N_phicorr, N_rcorr));
    tree_->Branch("cscRechitCluster2YSpread_corr",             &cscRechitCluster2YSpread_corr,  Form("cscRechitCluster2YSpread_corr[nCscRechitClusters2][%d][%d]/F",N_phicorr, N_rcorr));
    tree_->Branch("cscRechitCluster2PhiSpread_corr",             &cscRechitCluster2PhiSpread_corr,  Form("cscRechitCluster2PhiSpread_corr[nCscRechitClusters2][%d][%d]/F",N_phicorr, N_rcorr));
    tree_->Branch("cscRechitCluster2EtaSpread_corr",             &cscRechitCluster2EtaSpread_corr,  Form("cscRechitCluster2EtaSpread_corr[nCscRechitClusters2][%d][%d]/F",N_phicorr, N_rcorr));
    tree_->Branch("cscRechitCluster2EtaPhiSpread_corr",             &cscRechitCluster2EtaPhiSpread_corr,  Form("cscRechitCluster2EtaPhiSpread_corr[nCscRechitClusters2][%d][%d]/F",N_phicorr, N_rcorr));
    tree_->Branch("cscRechitCluster2RSpread_corr",             &cscRechitCluster2RSpread_corr,  Form("cscRechitCluster2RSpread_corr[nCscRechitClusters2][%d][%d]/F",N_phicorr, N_rcorr));

      tree_->Branch("cscRechitCluster2Phi",             cscRechitCluster2Phi,             "cscRechitCluster2Phi[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2Eta",             cscRechitCluster2Eta,             "cscRechitCluster2Eta[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2JetVetoPt",             cscRechitCluster2JetVetoPt,             "cscRechitCluster2JetVetoPt[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2JetVetoE",             cscRechitCluster2JetVetoE,             "cscRechitCluster2JetVetoE[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2MuonVetoPt",             cscRechitCluster2MuonVetoPt,             "cscRechitCluster2MuonVetoPt[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2MuonVetoE",             cscRechitCluster2MuonVetoE,             "cscRechitCluster2MuonVetoE[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2MuonVetoPhi",             cscRechitCluster2MuonVetoPhi,             "cscRechitCluster2MuonVetoPhi[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2MuonVetoEta",             cscRechitCluster2MuonVetoEta,             "cscRechitCluster2MuonVetoEta[nCscRechitClusters2]/F");

    tree_->Branch("cscRechitCluster2ZLep1",             cscRechitCluster2ZLep1,             "cscRechitCluster2ZLep1[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2ZLep2",             cscRechitCluster2ZLep2,             "cscRechitCluster2ZLep2[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2ZLep1Id",             cscRechitCluster2ZLep1Id,             "cscRechitCluster2ZLep1Id[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2ZLep2Id",             cscRechitCluster2ZLep2Id,             "cscRechitCluster2ZLep2Id[nCscRechitClusters2]/I");

    tree_->Branch("cscRechitCluster2ZLep1LooseIso",             cscRechitCluster2ZLep1LooseIso,             "cscRechitCluster2ZLep1LooseIso[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2ZLep1TightIso",             cscRechitCluster2ZLep1TightIso,             "cscRechitCluster2ZLep1TightIso[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2ZLep1VTightIso",             cscRechitCluster2ZLep1VTightIso,             "cscRechitCluster2ZLep1VTightIso[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2ZLep1VVTightIso",             cscRechitCluster2ZLep1VVTightIso,             "cscRechitCluster2ZLep1VVTightIso[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2ZLep1TightId",             cscRechitCluster2ZLep1TightId,             "cscRechitCluster2ZLep1TightId[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2ZLep2LooseIso",             cscRechitCluster2ZLep2LooseIso,             "cscRechitCluster2ZLep2LooseIso[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2ZLep2TightIso",             cscRechitCluster2ZLep2TightIso,             "cscRechitCluster2ZLep2TightIso[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2ZLep2VTightIso",             cscRechitCluster2ZLep2VTightIso,             "cscRechitCluster2ZLep2VTightIso[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2ZLep2VVTightIso",             cscRechitCluster2ZLep2VVTightIso,             "cscRechitCluster2ZLep2VVTightIso[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2ZLep2TightId",             cscRechitCluster2ZLep2TightId,             "cscRechitCluster2ZLep2TightId[nCscRechitClusters2]/O");




    tree_->Branch("cscRechitCluster2JetVetoPt_0p6",             cscRechitCluster2JetVetoPt_0p6,        "cscRechitCluster2JetVetoPt_0p6[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2JetVetoPt_0p8",             cscRechitCluster2JetVetoPt_0p8,        "cscRechitCluster2JetVetoPt_0p8[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2JetVetoE_0p6",             cscRechitCluster2JetVetoE_0p6,        "cscRechitCluster2JetVetoE_0p6[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2JetVetoE_0p8",             cscRechitCluster2JetVetoE_0p8,        "cscRechitCluster2JetVetoE_0p8[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2MuonVetoPt_0p6",             cscRechitCluster2MuonVetoPt_0p6,        "cscRechitCluster2MuonVetoPt_0p6[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2MuonVetoPt_0p8",             cscRechitCluster2MuonVetoPt_0p8,        "cscRechitCluster2MuonVetoPt_0p8[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2MuonVetoE_0p6",             cscRechitCluster2MuonVetoE_0p6,        "cscRechitCluster2MuonVetoE_0p6[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2MuonVetoE_0p8",             cscRechitCluster2MuonVetoE_0p8,        "cscRechitCluster2MuonVetoE_0p8[nCscRechitClusters2]/F");

      tree_->Branch("cscRechitCluster2MuonVetoLooseIso",             cscRechitCluster2MuonVetoLooseIso,             "cscRechitCluster2MuonVetoLooseIso[nCscRechitClusters2]/O");
      tree_->Branch("cscRechitCluster2MuonVetoTightIso",             cscRechitCluster2MuonVetoTightIso,             "cscRechitCluster2MuonVetoTightIso[nCscRechitClusters2]/O");
      tree_->Branch("cscRechitCluster2MuonVetoVTightIso",             cscRechitCluster2MuonVetoVTightIso,             "cscRechitCluster2MuonVetoVTightIso[nCscRechitClusters2]/O");
      tree_->Branch("cscRechitCluster2MuonVetoVVTightIso",             cscRechitCluster2MuonVetoVVTightIso,             "cscRechitCluster2MuonVetoVVTightIso[nCscRechitClusters2]/O");
      tree_->Branch("cscRechitCluster2MuonVetoTightId",             cscRechitCluster2MuonVetoTightId,             "cscRechitCluster2MuonVetoTightId[nCscRechitClusters2]/O");


    tree_->Branch("cscRechitCluster2MuonVetoIso",             cscRechitCluster2MuonVetoIso,             "cscRechitCluster2MuonVetoIso[nCscRechitClusters2]/O");
    tree_->Branch("cscRechitCluster2IsoMuonVetoPt",             cscRechitCluster2IsoMuonVetoPt,             "cscRechitCluster2IsoMuonVetoPt[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2IsoMuonVetoE",             cscRechitCluster2IsoMuonVetoE,             "cscRechitCluster2IsoMuonVetoE[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2IsoMuonVetoPhi",             cscRechitCluster2IsoMuonVetoPhi,             "cscRechitCluster2IsoMuonVetoPhi[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2IsoMuonVetoEta",             cscRechitCluster2IsoMuonVetoEta,             "cscRechitCluster2IsoMuonVetoEta[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2GenMuonVetoPt",             cscRechitCluster2GenMuonVetoPt,             "cscRechitCluster2GenMuonVetoPt[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2GenMuonVetoE",             cscRechitCluster2GenMuonVetoE,             "cscRechitCluster2GenMuonVetoE[nCscRechitClusters2]/F");

    tree_->Branch("cscRechitCluster2Size",             cscRechitCluster2Size,             "cscRechitCluster2Size[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2Me11Ratio",             cscRechitCluster2Me11Ratio,             "cscRechitCluster2Me11Ratio[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2Me12Ratio",             cscRechitCluster2Me12Ratio,             "cscRechitCluster2Me12Ratio[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2NStation",             cscRechitCluster2NStation,             "cscRechitCluster2NStation[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NStation5",             cscRechitCluster2NStation5,             "cscRechitCluster2NStation5[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NStation10perc",             cscRechitCluster2NStation10perc,             "cscRechitCluster2NStation10perc[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2AvgStation",             cscRechitCluster2AvgStation,             "cscRechitCluster2AvgStation[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2AvgStation5",             cscRechitCluster2AvgStation5,             "cscRechitCluster2AvgStation5[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2AvgStation10perc",             cscRechitCluster2AvgStation10perc,             "cscRechitCluster2AvgStation10perc[nCscRechitClusters2]/F");


    tree_->Branch("cscRechitCluster2MaxStation",             cscRechitCluster2MaxStation,             "cscRechitCluster2MaxStation[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2MaxStationRatio",             cscRechitCluster2MaxStationRatio,             "cscRechitCluster2MaxStationRatio[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2NChamber",             cscRechitCluster2NChamber,             "cscRechitCluster2NChamber[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2MaxChamber",             cscRechitCluster2MaxChamber,             "cscRechitCluster2MaxChamber[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2MaxChamberRatio",             cscRechitCluster2MaxChamberRatio,             "cscRechitCluster2MaxChamberRatio[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2NRechitChamberPlus11",             cscRechitCluster2NRechitChamberPlus11,             "cscRechitCluster2NRechitChamberPlus11[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberPlus12",             cscRechitCluster2NRechitChamberPlus12,             "cscRechitCluster2NRechitChamberPlus12[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberPlus13",             cscRechitCluster2NRechitChamberPlus13,             "cscRechitCluster2NRechitChamberPlus13[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberPlus21",             cscRechitCluster2NRechitChamberPlus21,             "cscRechitCluster2NRechitChamberPlus21[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberPlus22",             cscRechitCluster2NRechitChamberPlus22,             "cscRechitCluster2NRechitChamberPlus22[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberPlus31",             cscRechitCluster2NRechitChamberPlus31,             "cscRechitCluster2NRechitChamberPlus31[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberPlus32",             cscRechitCluster2NRechitChamberPlus32,             "cscRechitCluster2NRechitChamberPlus32[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberPlus41",             cscRechitCluster2NRechitChamberPlus41,             "cscRechitCluster2NRechitChamberPlus41[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberPlus42",             cscRechitCluster2NRechitChamberPlus42,             "cscRechitCluster2NRechitChamberPlus42[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberMinus11",             cscRechitCluster2NRechitChamberMinus11,             "cscRechitCluster2NRechitChamberMinus11[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberMinus12",             cscRechitCluster2NRechitChamberMinus12,             "cscRechitCluster2NRechitChamberMinus12[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberMinus13",             cscRechitCluster2NRechitChamberMinus13,             "cscRechitCluster2NRechitChamberMinus13[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberMinus21",             cscRechitCluster2NRechitChamberMinus21,             "cscRechitCluster2NRechitChamberMinus21[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberMinus22",             cscRechitCluster2NRechitChamberMinus22,             "cscRechitCluster2NRechitChamberMinus22[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberMinus31",             cscRechitCluster2NRechitChamberMinus31,             "cscRechitCluster2NRechitChamberMinus31[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberMinus32",             cscRechitCluster2NRechitChamberMinus32,             "cscRechitCluster2NRechitChamberMinus32[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberMinus41",             cscRechitCluster2NRechitChamberMinus41,             "cscRechitCluster2NRechitChamberMinus41[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2NRechitChamberMinus42",             cscRechitCluster2NRechitChamberMinus42,             "cscRechitCluster2NRechitChamberMinus42[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2Met_dPhi",             cscRechitCluster2Met_dPhi,             "cscRechitCluster2Met_dPhi[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2MetXYCorr_dPhi",             cscRechitCluster2MetXYCorr_dPhi,             "cscRechitCluster2MetXYCorr_dPhi[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_cscRechits_0p4",             cscRechitCluster2_match_cscRechits_0p4,             "cscRechitCluster2_match_cscRechits_0p4[nCscRechitClusters2]/I");

    tree_->Branch("cscRechitCluster2_match_cscSeg_0p4",             cscRechitCluster2_match_cscSeg_0p4,             "cscRechitCluster2_match_cscSeg_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_ME11Seg_0p4",             cscRechitCluster2_match_ME11Seg_0p4,             "cscRechitCluster2_match_ME11Seg_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_ME12Seg_0p4",             cscRechitCluster2_match_ME12Seg_0p4,             "cscRechitCluster2_match_ME12Seg_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_cscSeg_0p6",             cscRechitCluster2_match_cscSeg_0p6,             "cscRechitCluster2_match_cscSeg_0p6[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_ME11Seg_0p6",             cscRechitCluster2_match_ME11Seg_0p6,             "cscRechitCluster2_match_ME11Seg_0p6[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_ME12Seg_0p6",             cscRechitCluster2_match_ME12Seg_0p6,             "cscRechitCluster2_match_ME12Seg_0p6[nCscRechitClusters2]/I");


    tree_->Branch("cscRechitCluster2_match_dtRechits_0p4",             cscRechitCluster2_match_dtRechits_0p4,             "cscRechitCluster2_match_dtRechits_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_MB1_0p4",             cscRechitCluster2_match_MB1_0p4,             "cscRechitCluster2_match_MB1_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_dtRechits_0p6",             cscRechitCluster2_match_dtRechits_0p6,             "cscRechitCluster2_match_dtRechits_0p6[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_MB1_0p6",             cscRechitCluster2_match_MB1_0p6,             "cscRechitCluster2_match_MB1_0p6[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_dtSeg_0p4",             cscRechitCluster2_match_dtSeg_0p4,             "cscRechitCluster2_match_dtSeg_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_MB1Seg_0p4",             cscRechitCluster2_match_MB1Seg_0p4,             "cscRechitCluster2_match_MB1Seg_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_dtSeg_0p6",             cscRechitCluster2_match_dtSeg_0p6,             "cscRechitCluster2_match_dtSeg_0p6[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_MB1Seg_0p6",             cscRechitCluster2_match_MB1Seg_0p6,             "cscRechitCluster2_match_MB1Seg_0p6[nCscRechitClusters2]/I");

    tree_->Branch("cscRechitCluster2_match_RB1_0p4",             cscRechitCluster2_match_RB1_0p4,             "cscRechitCluster2_match_RB1_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_RE12_0p4",             cscRechitCluster2_match_RE12_0p4,             "cscRechitCluster2_match_RE12_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_RB1_0p6",             cscRechitCluster2_match_RB1_0p6,             "cscRechitCluster2_match_RB1_0p6[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_RE12_0p6",             cscRechitCluster2_match_RE12_0p6,             "cscRechitCluster2_match_RE12_0p6[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_highEta_0p4",             cscRechitCluster2_match_highEta_0p4,             "cscRechitCluster2_match_highEta_0p4[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_highEta_0p6",             cscRechitCluster2_match_highEta_0p6,             "cscRechitCluster2_match_highEta_0p6[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_highEta_0p8",             cscRechitCluster2_match_highEta_0p8,             "cscRechitCluster2_match_highEta_0p8[nCscRechitClusters2]/I");
    tree_->Branch("cscRechitCluster2_match_cluster_dR",             cscRechitCluster2_match_cluster_dR,             "cscRechitCluster2_match_cluster_dR[nCscRechitClusters2]/F");
    tree_->Branch("cscRechitCluster2_match_cluster_index",             cscRechitCluster2_match_cluster_index,             "cscRechitCluster2_match_cluster_index[nCscRechitClusters2]/I");


    tree_->Branch("nCscRechitClusters3",             &nCscRechitClusters3, "nCscRechitClusters3/I");
    tree_->Branch("cscRechitCluster3_match_Me1112_0p4",             cscRechitCluster3_match_Me1112_0p4,             "cscRechitCluster3_match_Me1112_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me1112_0p6",             cscRechitCluster3_match_Me1112_0p6,             "cscRechitCluster3_match_Me1112_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me1112_0p8",             cscRechitCluster3_match_Me1112_0p8,             "cscRechitCluster3_match_Me1112_0p8[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me11_0p4",             cscRechitCluster3_match_Me11_0p4,             "cscRechitCluster3_match_Me11_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me11_0p6",             cscRechitCluster3_match_Me11_0p6,             "cscRechitCluster3_match_Me11_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me11_0p8",             cscRechitCluster3_match_Me11_0p8,             "cscRechitCluster3_match_Me11_0p8[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me12_0p4",             cscRechitCluster3_match_Me12_0p4,             "cscRechitCluster3_match_Me12_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me12_0p6",             cscRechitCluster3_match_Me12_0p6,             "cscRechitCluster3_match_Me12_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me12_0p8",             cscRechitCluster3_match_Me12_0p8,             "cscRechitCluster3_match_Me12_0p8[nCscRechitClusters3]/I");

    tree_->Branch("cscRechitCluster3_match_gParticle_id",             cscRechitCluster3_match_gParticle_id,             "cscRechitCluster3_match_gParticle_id[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_gParticle",             cscRechitCluster3_match_gParticle,             "cscRechitCluster3_match_gParticle[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3_match_gParticle_minDeltaR",             cscRechitCluster3_match_gParticle_minDeltaR,             "cscRechitCluster3_match_gParticle_minDeltaR[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gParticle_index",             cscRechitCluster3_match_gParticle_index,             "cscRechitCluster3_match_gParticle_index[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_gParticle_eta",             cscRechitCluster3_match_gParticle_eta,             "cscRechitCluster3_match_gParticle_eta[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gParticle_phi",             cscRechitCluster3_match_gParticle_phi,             "cscRechitCluster3_match_gParticle_phi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gParticle_E",             cscRechitCluster3_match_gParticle_E,             "cscRechitCluster3_match_gParticle_E[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gParticle_pt",             cscRechitCluster3_match_gParticle_pt,             "cscRechitCluster3_match_gParticle_pt[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gParticle_MotherId",             cscRechitCluster3_match_gParticle_MotherId,             "cscRechitCluster3_match_gParticle_MotherId[nCscRechitClusters3]/I");

    tree_->Branch("cscRechitCluster3_match_gLLP",             cscRechitCluster3_match_gLLP,             "cscRechitCluster3_match_gLLP[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3_match_gLLP_minDeltaR",             cscRechitCluster3_match_gLLP_minDeltaR,             "cscRechitCluster3_match_gLLP_minDeltaR[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_index",             cscRechitCluster3_match_gLLP_index,             "cscRechitCluster3_match_gLLP_index[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_gLLP_eta",             cscRechitCluster3_match_gLLP_eta, "cscRechitCluster3_match_gLLP_eta[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_phi",             cscRechitCluster3_match_gLLP_phi, "cscRechitCluster3_match_gLLP_phi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_decay_r",             cscRechitCluster3_match_gLLP_decay_r, "cscRechitCluster3_match_gLLP_decay_r[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_decay_x",             cscRechitCluster3_match_gLLP_decay_x, "cscRechitCluster3_match_gLLP_decay_x[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_decay_y",             cscRechitCluster3_match_gLLP_decay_y, "cscRechitCluster3_match_gLLP_decay_y[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_decay_z",             cscRechitCluster3_match_gLLP_decay_z, "cscRechitCluster3_match_gLLP_decay_z[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_ctau",             cscRechitCluster3_match_gLLP_ctau, "cscRechitCluster3_match_gLLP_ctau[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_beta",             cscRechitCluster3_match_gLLP_beta, "cscRechitCluster3_match_gLLP_beta[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_csc",             cscRechitCluster3_match_gLLP_csc, "cscRechitCluster3_match_gLLP_csc[nCscRechitClusters3]/O");

    tree_->Branch("cscRechitCluster3X",             cscRechitCluster3X,             "cscRechitCluster3X[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3Y",             cscRechitCluster3Y,             "cscRechitCluster3Y[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3Z",             cscRechitCluster3Z,             "cscRechitCluster3Z[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3Time",             cscRechitCluster3Time,             "cscRechitCluster3Time[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3TimeTotal",             cscRechitCluster3TimeTotal,             "cscRechitCluster3TimeTotal[nCscRechitClusters3]/F");

    tree_->Branch("cscRechitCluster3TimeSpread",             cscRechitCluster3TimeSpread,             "cscRechitCluster3TimeSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3GenMuonDeltaR",             cscRechitCluster3GenMuonDeltaR,             "cscRechitCluster3GenMuonDeltaR[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3XYSpread",             cscRechitCluster3XYSpread,             "cscRechitCluster3XYSpread[nCscRechitClusters3]/F");

    tree_->Branch("cscRechitCluster3MajorAxis",             cscRechitCluster3MajorAxis,             "cscRechitCluster3MajorAxis[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MinorAxis",             cscRechitCluster3MinorAxis,             "cscRechitCluster3MinorAxis[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3EtaPhiSpread",             cscRechitCluster3EtaPhiSpread,             "cscRechitCluster3EtaPhiSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3PhiSpread",             cscRechitCluster3PhiSpread,             "cscRechitCluster3PhiSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3EtaSpread",             cscRechitCluster3EtaSpread,             "cscRechitCluster3EtaSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3XSpread",             cscRechitCluster3XSpread,             "cscRechitCluster3XSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3RSpread",             cscRechitCluster3RSpread,             "cscRechitCluster3RSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3YSpread",             cscRechitCluster3YSpread,             "cscRechitCluster3YSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3ZSpread",             cscRechitCluster3ZSpread,             "cscRechitCluster3ZSpread[nCscRechitClusters3]/F");


    tree_->Branch("cscRechitCluster3XYSpread_corr",             &cscRechitCluster3XYSpread_corr, Form("cscRechitCluster3XYSpread_corr[nCscRechitClusters3][%d][%d]/F",N_phicorr, N_rcorr));
    tree_->Branch("cscRechitCluster3XSpread_corr",             &cscRechitCluster3XSpread_corr, Form("cscRechitCluster3XSpread_corr[nCscRechitClusters3][%d][%d]/F",N_phicorr, N_rcorr));
    tree_->Branch("cscRechitCluster3YSpread_corr",             &cscRechitCluster3YSpread_corr, Form("cscRechitCluster3YSpread_corr[nCscRechitClusters3][%d][%d]/F",N_phicorr, N_rcorr));
    tree_->Branch("cscRechitCluster3PhiSpread_corr",             &cscRechitCluster3PhiSpread_corr, Form("cscRechitCluster3PhiSpread_corr[nCscRechitClusters3][%d][%d]/F",N_phicorr, N_rcorr));
    tree_->Branch("cscRechitCluster3EtaSpread_corr",             &cscRechitCluster3EtaSpread_corr, Form("cscRechitCluster3EtaSpread_corr[nCscRechitClusters3][%d][%d]/F",N_phicorr, N_rcorr));
    tree_->Branch("cscRechitCluster3EtaPhiSpread_corr",             &cscRechitCluster3EtaPhiSpread_corr, Form("cscRechitCluster3EtaPhiSpread_corr[nCscRechitClusters3][%d][%d]/F",N_phicorr, N_rcorr));
    tree_->Branch("cscRechitCluster3RSpread_corr",             &cscRechitCluster3RSpread_corr, Form("cscRechitCluster3RSpread_corr[nCscRechitClusters3][%d][%d]/F",N_phicorr, N_rcorr));

    tree_->Branch("cscRechitCluster3XYSpread_corr2",             cscRechitCluster3XYSpread_corr2, "cscRechitCluster3XYSpread_corr2[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3XSpread_corr2",             cscRechitCluster3XSpread_corr2, "cscRechitCluster3XSpread_corr2[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3YSpread_corr2",            cscRechitCluster3YSpread_corr2, "cscRechitCluster3YSpread_corr2[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3PhiSpread_corr2",         cscRechitCluster3PhiSpread_corr2, "cscRechitCluster3PhiSpread_corr2[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3EtaSpread_corr2",         cscRechitCluster3EtaSpread_corr2, "cscRechitCluster3EtaSpread_corr2[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3EtaPhiSpread_corr2",      cscRechitCluster3EtaPhiSpread_corr2, "cscRechitCluster3EtaPhiSpread_corr2[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3RSpread_corr2",           cscRechitCluster3RSpread_corr2, "cscRechitCluster3RSpread_corr2[nCscRechitClusters3]/F");

    tree_->Branch("cscRechitCluster3Phi",             cscRechitCluster3Phi,             "cscRechitCluster3Phi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3Eta",             cscRechitCluster3Eta,             "cscRechitCluster3Eta[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoPt",             cscRechitCluster3JetVetoPt,             "cscRechitCluster3JetVetoPt[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoE",             cscRechitCluster3JetVetoE,             "cscRechitCluster3JetVetoE[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoPt",             cscRechitCluster3MuonVetoPt,             "cscRechitCluster3MuonVetoPt[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoE",             cscRechitCluster3MuonVetoE,             "cscRechitCluster3MuonVetoE[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoPhi",             cscRechitCluster3MuonVetoPhi,             "cscRechitCluster3MuonVetoPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoEta",             cscRechitCluster3MuonVetoEta,             "cscRechitCluster3MuonVetoEta[nCscRechitClusters3]/F");




    tree_->Branch("cscRechitCluster3JetVetoPt_0p6",             cscRechitCluster3JetVetoPt_0p6,        "cscRechitCluster3JetVetoPt_0p6[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoPt_0p8",             cscRechitCluster3JetVetoPt_0p8,        "cscRechitCluster3JetVetoPt_0p8[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoE_0p6",             cscRechitCluster3JetVetoE_0p6,        "cscRechitCluster3JetVetoE_0p6[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoE_0p8",             cscRechitCluster3JetVetoE_0p8,        "cscRechitCluster3JetVetoE_0p8[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoPt_0p6",             cscRechitCluster3MuonVetoPt_0p6,        "cscRechitCluster3MuonVetoPt_0p6[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoPt_0p8",             cscRechitCluster3MuonVetoPt_0p8,        "cscRechitCluster3MuonVetoPt_0p8[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoE_0p6",             cscRechitCluster3MuonVetoE_0p6,        "cscRechitCluster3MuonVetoE_0p6[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoE_0p8",             cscRechitCluster3MuonVetoE_0p8,        "cscRechitCluster3MuonVetoE_0p8[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoLooseIso",             cscRechitCluster3MuonVetoLooseIso,             "cscRechitCluster3MuonVetoLooseIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3MuonVetoTightIso",             cscRechitCluster3MuonVetoTightIso,             "cscRechitCluster3MuonVetoTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3MuonVetoVTightIso",             cscRechitCluster3MuonVetoVTightIso,             "cscRechitCluster3MuonVetoVTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3MuonVetoVVTightIso",             cscRechitCluster3MuonVetoVVTightIso,             "cscRechitCluster3MuonVetoVVTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3MuonVetoTightId",             cscRechitCluster3MuonVetoTightId,             "cscRechitCluster3MuonVetoTightId[nCscRechitClusters3]/O");


    tree_->Branch("cscRechitCluster3MuonVetoIso",             cscRechitCluster3MuonVetoIso,             "cscRechitCluster3MuonVetoIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3IsoMuonVetoPt",             cscRechitCluster3IsoMuonVetoPt,             "cscRechitCluster3IsoMuonVetoPt[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3IsoMuonVetoE",             cscRechitCluster3IsoMuonVetoE,             "cscRechitCluster3IsoMuonVetoE[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3IsoMuonVetoPhi",             cscRechitCluster3IsoMuonVetoPhi,             "cscRechitCluster3IsoMuonVetoPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3IsoMuonVetoEta",             cscRechitCluster3IsoMuonVetoEta,             "cscRechitCluster3IsoMuonVetoEta[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3GenMuonVetoPt",             cscRechitCluster3GenMuonVetoPt,             "cscRechitCluster3GenMuonVetoPt[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3GenMuonVetoE",             cscRechitCluster3GenMuonVetoE,             "cscRechitCluster3GenMuonVetoE[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3ZLep1",             cscRechitCluster3ZLep1,             "cscRechitCluster3ZLep1[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2",             cscRechitCluster3ZLep2,             "cscRechitCluster3ZLep2[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep1Id",             cscRechitCluster3ZLep1Id,             "cscRechitCluster3ZLep1Id[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3ZLep2Id",             cscRechitCluster3ZLep2Id,             "cscRechitCluster3ZLep2Id[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster2ZLep1LooseIso",             cscRechitCluster3ZLep1LooseIso,             "cscRechitCluster3ZLep1LooseIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep1TightIso",             cscRechitCluster3ZLep1TightIso,             "cscRechitCluster3ZLep1TightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep1VTightIso",             cscRechitCluster3ZLep1VTightIso,             "cscRechitCluster3ZLep1VTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep1VVTightIso",             cscRechitCluster3ZLep1VVTightIso,             "cscRechitCluster3ZLep1VVTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep1TightId",             cscRechitCluster3ZLep1TightId,             "cscRechitCluster3ZLep1TightId[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2LooseIso",             cscRechitCluster3ZLep2LooseIso,             "cscRechitCluster3ZLep2LooseIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2TightIso",             cscRechitCluster3ZLep2TightIso,             "cscRechitCluster3ZLep2TightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2VTightIso",             cscRechitCluster3ZLep2VTightIso,             "cscRechitCluster3ZLep2VTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2VVTightIso",             cscRechitCluster3ZLep2VVTightIso,             "cscRechitCluster3ZLep2VVTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2TightId",             cscRechitCluster3ZLep2TightId,             "cscRechitCluster3ZLep2TightId[nCscRechitClusters3]/O");


    tree_->Branch("cscRechitCluster3Size",             cscRechitCluster3Size,             "cscRechitCluster3Size[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3Me11Ratio",             cscRechitCluster3Me11Ratio,             "cscRechitCluster3Me11Ratio[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3Me12Ratio",             cscRechitCluster3Me12Ratio,             "cscRechitCluster3Me12Ratio[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3NStation",             cscRechitCluster3NStation,             "cscRechitCluster3NStation[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NStation5",             cscRechitCluster3NStation5,             "cscRechitCluster3NStation5[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NStation10perc",             cscRechitCluster3NStation10perc,             "cscRechitCluster3NStation10perc[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3AvgStation",             cscRechitCluster3AvgStation,             "cscRechitCluster3AvgStation[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3AvgStation5",             cscRechitCluster3AvgStation5,             "cscRechitCluster3AvgStation5[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3AvgStation10perc",             cscRechitCluster3AvgStation10perc,             "cscRechitCluster3AvgStation10perc[nCscRechitClusters3]/F");


    tree_->Branch("cscRechitCluster3MaxStation",             cscRechitCluster3MaxStation,             "cscRechitCluster3MaxStation[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3MaxStationRatio",             cscRechitCluster3MaxStationRatio,             "cscRechitCluster3MaxStationRatio[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3NChamber",             cscRechitCluster3NChamber,             "cscRechitCluster3NChamber[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3MaxChamber",             cscRechitCluster3MaxChamber,             "cscRechitCluster3MaxChamber[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3MaxChamberRatio",             cscRechitCluster3MaxChamberRatio,             "cscRechitCluster3MaxChamberRatio[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus11",             cscRechitCluster3NRechitChamberPlus11,             "cscRechitCluster3NRechitChamberPlus11[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus12",             cscRechitCluster3NRechitChamberPlus12,             "cscRechitCluster3NRechitChamberPlus12[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus13",             cscRechitCluster3NRechitChamberPlus13,             "cscRechitCluster3NRechitChamberPlus13[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus21",             cscRechitCluster3NRechitChamberPlus21,             "cscRechitCluster3NRechitChamberPlus21[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus22",             cscRechitCluster3NRechitChamberPlus22,             "cscRechitCluster3NRechitChamberPlus22[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus31",             cscRechitCluster3NRechitChamberPlus31,             "cscRechitCluster3NRechitChamberPlus31[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus32",             cscRechitCluster3NRechitChamberPlus32,             "cscRechitCluster3NRechitChamberPlus32[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus41",             cscRechitCluster3NRechitChamberPlus41,             "cscRechitCluster3NRechitChamberPlus41[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus42",             cscRechitCluster3NRechitChamberPlus42,             "cscRechitCluster3NRechitChamberPlus42[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus11",             cscRechitCluster3NRechitChamberMinus11,             "cscRechitCluster3NRechitChamberMinus11[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus12",             cscRechitCluster3NRechitChamberMinus12,             "cscRechitCluster3NRechitChamberMinus12[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus13",             cscRechitCluster3NRechitChamberMinus13,             "cscRechitCluster3NRechitChamberMinus13[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus21",             cscRechitCluster3NRechitChamberMinus21,             "cscRechitCluster3NRechitChamberMinus21[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus22",             cscRechitCluster3NRechitChamberMinus22,             "cscRechitCluster3NRechitChamberMinus22[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus31",             cscRechitCluster3NRechitChamberMinus31,             "cscRechitCluster3NRechitChamberMinus31[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus32",             cscRechitCluster3NRechitChamberMinus32,             "cscRechitCluster3NRechitChamberMinus32[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus41",             cscRechitCluster3NRechitChamberMinus41,             "cscRechitCluster3NRechitChamberMinus41[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus42",             cscRechitCluster3NRechitChamberMinus42,             "cscRechitCluster3NRechitChamberMinus42[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3Met_dPhi",             cscRechitCluster3Met_dPhi,             "cscRechitCluster3Met_dPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MetXYCorr_dPhi",             cscRechitCluster3MetXYCorr_dPhi,             "cscRechitCluster3MetXYCorr_dPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_cscRechits_0p4",             cscRechitCluster3_match_cscRechits_0p4,             "cscRechitCluster3_match_cscRechits_0p4[nCscRechitClusters3]/I");

    tree_->Branch("cscRechitCluster3_match_cscSeg_0p4",             cscRechitCluster3_match_cscSeg_0p4,             "cscRechitCluster3_match_cscSeg_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_ME11Seg_0p4",             cscRechitCluster3_match_ME11Seg_0p4,             "cscRechitCluster3_match_ME11Seg_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_ME12Seg_0p4",             cscRechitCluster3_match_ME12Seg_0p4,             "cscRechitCluster3_match_ME12Seg_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_cscSeg_0p6",             cscRechitCluster3_match_cscSeg_0p6,             "cscRechitCluster3_match_cscSeg_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_ME11Seg_0p6",             cscRechitCluster3_match_ME11Seg_0p6,             "cscRechitCluster3_match_ME11Seg_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_ME12Seg_0p6",             cscRechitCluster3_match_ME12Seg_0p6,             "cscRechitCluster3_match_ME12Seg_0p6[nCscRechitClusters3]/I");


    tree_->Branch("cscRechitCluster3_match_dtRechits_0p4",             cscRechitCluster3_match_dtRechits_0p4,             "cscRechitCluster3_match_dtRechits_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_MB1_0p4",             cscRechitCluster3_match_MB1_0p4,             "cscRechitCluster3_match_MB1_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_dtRechits_0p6",             cscRechitCluster3_match_dtRechits_0p6,             "cscRechitCluster3_match_dtRechits_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_MB1_0p6",             cscRechitCluster3_match_MB1_0p6,             "cscRechitCluster3_match_MB1_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_dtSeg_0p4",             cscRechitCluster3_match_dtSeg_0p4,             "cscRechitCluster3_match_dtSeg_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_MB1Seg_0p4",             cscRechitCluster3_match_MB1Seg_0p4,             "cscRechitCluster3_match_MB1Seg_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_dtSeg_0p6",             cscRechitCluster3_match_dtSeg_0p6,             "cscRechitCluster3_match_dtSeg_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_MB1Seg_0p6",             cscRechitCluster3_match_MB1Seg_0p6,             "cscRechitCluster3_match_MB1Seg_0p6[nCscRechitClusters3]/I");

    tree_->Branch("cscRechitCluster3_match_RB1_0p4",             cscRechitCluster3_match_RB1_0p4,             "cscRechitCluster3_match_RB1_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_RE12_0p4",             cscRechitCluster3_match_RE12_0p4,             "cscRechitCluster3_match_RE12_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_RB1_0p6",             cscRechitCluster3_match_RB1_0p6,             "cscRechitCluster3_match_RB1_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_RE12_0p6",             cscRechitCluster3_match_RE12_0p6,             "cscRechitCluster3_match_RE12_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_highEta_0p4",             cscRechitCluster3_match_highEta_0p4,             "cscRechitCluster3_match_highEta_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_highEta_0p6",             cscRechitCluster3_match_highEta_0p6,             "cscRechitCluster3_match_highEta_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_highEta_0p8",             cscRechitCluster3_match_highEta_0p8,             "cscRechitCluster3_match_highEta_0p8[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_cluster_dR",             cscRechitCluster3_match_cluster_dR,             "cscRechitCluster3_match_cluster_dR[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_cluster_index",             cscRechitCluster3_match_cluster_index,             "cscRechitCluster3_match_cluster_index[nCscRechitClusters3]/I");


      // INTIME CSC cluster
  /*
  tree_->Branch("nCscITClusters",             &nCscITClusters, "nCscITClusters/I");
  tree_->Branch("nCsc_JetVetoITCluster0p4",             &nCsc_JetVetoITCluster0p4, "nCsc_JetVetoITCluster0p4/I");
  tree_->Branch("nCsc_JetMuonVetoITCluster0p4",             &nCsc_JetMuonVetoITCluster0p4, "nCsc_JetMuonVetoITCluster0p4/I");
  tree_->Branch("nCsc_JetVetoITCluster0p4_Me1112Veto",             &nCsc_JetVetoITCluster0p4_Me1112Veto, "nCsc_JetVetoITCluster0p4_Me1112Veto/I");
  tree_->Branch("nCsc_JetMuonVetoITCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoITCluster0p4_Me1112Veto, "nCsc_JetMuonVetoITCluster0p4_Me1112Veto/I");
  tree_->Branch("cscITClusterX",             cscITClusterX,             "cscITClusterX[nCscITClusters]/F");
  tree_->Branch("cscITClusterY",             cscITClusterY,             "cscITClusterY[nCscITClusters]/F");
  tree_->Branch("cscITClusterZ",             cscITClusterZ,             "cscITClusterZ[nCscITClusters]/F");
  tree_->Branch("cscITClusterTime",             cscITClusterTime,             "cscITClusterTime[nCscITClusters]/F");
  tree_->Branch("cscITClusterTimeSpread",             cscITClusterTimeSpread,             "cscITClusterTimeSpread[nCscITClusters]/F");
  tree_->Branch("cscITClusterTimeRMS",             cscITClusterTimeRMS,             "cscITClusterTimeRMS[nCscITClusters]/F");
  tree_->Branch("cscITClusterGenMuonDeltaR",             cscITClusterGenMuonDeltaR,             "cscITClusterGenMuonDeltaR[nCscITClusters]/F");

  tree_->Branch("cscITClusterRadius",             cscITClusterRadius,             "cscITClusterRadius[nCscITClusters]/F");
  tree_->Branch("cscITClusterMajorAxis",             cscITClusterMajorAxis,             "cscITClusterMajorAxis[nCscITClusters]/F");
  tree_->Branch("cscITClusterMinorAxis",             cscITClusterMinorAxis,             "cscITClusterMinorAxis[nCscITClusters]/F");
  tree_->Branch("cscITClusterEtaPhiSpread",             cscITClusterEtaPhiSpread,             "cscITClusterEtaPhiSpread[nCscITClusters]/F");
  tree_->Branch("cscITClusterPhiSpread",             cscITClusterPhiSpread,             "cscITClusterPhiSpread[nCscITClusters]/F");
  tree_->Branch("cscITClusterEtaSpread",             cscITClusterEtaSpread,             "cscITClusterEtaSpread[nCscITClusters]/F");
  tree_->Branch("cscITClusterXSpread",             cscITClusterXSpread,             "cscITClusterXSpread[nCscITClusters]/F");
  tree_->Branch("cscITClusterYSpread",             cscITClusterYSpread,             "cscITClusterYSpread[nCscITClusters]/F");
  tree_->Branch("cscITClusterZSpread",             cscITClusterZSpread,             "cscITClusterZSpread[nCscITClusters]/F");
  tree_->Branch("cscITClusterPhi",             cscITClusterPhi,             "cscITClusterPhi[nCscITClusters]/F");
  tree_->Branch("cscITClusterEta",             cscITClusterEta,             "cscITClusterEta[nCscITClusters]/F");
  tree_->Branch("cscITClusterJetVeto",             cscITClusterJetVeto,             "cscITClusterJetVeto[nCscITClusters]/F");
  tree_->Branch("cscITClusterMuonVeto",             cscITClusterMuonVeto,             "cscITClusterMuonVeto[nCscITClusters]/F");
  tree_->Branch("cscITClusterCaloJetVeto",             cscITClusterCaloJetVeto,             "cscITClusterCaloJetVeto[nCscITClusters]/F");
  tree_->Branch("cscITClusterJetVetoE",             cscITClusterJetVetoE,             "cscITClusterJetVetoE[nCscITClusters]/F");
  tree_->Branch("cscITClusterMuonVetoE",             cscITClusterMuonVetoE,             "cscITClusterMuonVetoE[nCscITClusters]/F");
  tree_->Branch("cscITClusterCaloJetVetoE",             cscITClusterCaloJetVetoE,             "cscITClusterCaloJetVetoE[nCscITClusters]/F");
  tree_->Branch("cscITClusterSize",             cscITClusterSize,             "cscITClusterSize[nCscITClusters]/I");
  tree_->Branch("cscITClusterMe11Ratio",             cscITClusterMe11Ratio,             "cscITClusterMe11Ratio[nCscITClusters]/F");
  tree_->Branch("cscITClusterMe12Ratio",             cscITClusterMe12Ratio,             "cscITClusterMe12Ratio[nCscITClusters]/F");

  tree_->Branch("cscITClusterNStation",             cscITClusterNStation,             "cscITClusterNStation[nCscITClusters]/I");
  tree_->Branch("cscITClusterMaxStation",             cscITClusterMaxStation,             "cscITClusterMaxStation[nCscITClusters]/I");
  tree_->Branch("cscITClusterMaxStationRatio",             cscITClusterMaxStationRatio,             "cscITClusterMaxStationRatio[nCscITClusters]/F");

  tree_->Branch("cscITClusterNChamber",             cscITClusterNChamber,             "cscITClusterNChamber[nCscITClusters]/I");
  tree_->Branch("cscITClusterMaxChamber",             cscITClusterMaxChamber,             "cscITClusterMaxChamber[nCscITClusters]/I");
  tree_->Branch("cscITClusterMaxChamberRatio",             cscITClusterMaxChamberRatio,             "cscITClusterMaxChamberRatio[nCscITClusters]/F");
  tree_->Branch("cscITClusterVertexR",             cscITClusterVertexR,             "cscITClusterVertexR[nCscITClusters]/F");
  tree_->Branch("cscITClusterVertexZ",             cscITClusterVertexZ,             "cscITClusterVertexZ[nCscITClusters]/F");
  tree_->Branch("cscITClusterVertexDis",             cscITClusterVertexDis,             "cscITClusterVertexDis[nCscITClusters]/F");
  tree_->Branch("cscITClusterVertexChi2",             cscITClusterVertexChi2,             "cscITClusterVertexChi2[nCscITClusters]/F");
  tree_->Branch("cscITClusterVertexN1",             cscITClusterVertexN1,             "cscITClusterVertexN1[nCscITClusters]/I");
  tree_->Branch("cscITClusterVertexN5",             cscITClusterVertexN5,             "cscITClusterVertexN5[nCscITClusters]/I");
  tree_->Branch("cscITClusterVertexN10",             cscITClusterVertexN10,             "cscITClusterVertexN10[nCscITClusters]/I");
  tree_->Branch("cscITClusterVertexN15",             cscITClusterVertexN15,             "cscITClusterVertexN15[nCscITClusters]/I");
  tree_->Branch("cscITClusterVertexN20",             cscITClusterVertexN20,             "cscITClusterVertexN20[nCscITClusters]/I");
  tree_->Branch("cscITClusterVertexN",             cscITClusterVertexN,             "cscITClusterVertexN[nCscITClusters]/I");
  tree_->Branch("cscITCluster_match_cscCluster_index",             cscITCluster_match_cscCluster_index,             "cscITCluster_match_cscCluster_index[nCscITClusters]/I");
  tree_->Branch("cscITCluster_cscCluster_SizeRatio",             cscITCluster_cscCluster_SizeRatio,             "cscITCluster_cscCluster_SizeRatio[nCscITClusters]/F");
  */
  //gLLP branches
  tree_->Branch("gLLP_eta",          gLLP_eta,          "gLLP_eta[2]/F");
  tree_->Branch("gLLP_phi",          gLLP_phi,          "gLLP_phi[2]/F");
  tree_->Branch("gLLP_csc",          gLLP_csc,          "gLLP_csc[2]/F");
  tree_->Branch("gLLP_beta",          gLLP_beta,          "gLLP_beta[2]/F");

  tree_->Branch("gLLP_ctau",          gLLP_ctau,          "gLLP_ctau[2]/F");

  tree_->Branch("gLLP_decay_vertex_r",          gLLP_decay_vertex_r,          "gLLP_decay_vertex_r[2]/F");
  tree_->Branch("gLLP_decay_vertex_x",          gLLP_decay_vertex_x,          "gLLP_decay_vertex_x[2]/F");
  tree_->Branch("gLLP_decay_vertex_y",          gLLP_decay_vertex_y,          "gLLP_decay_vertex_y[2]/F");
  tree_->Branch("gLLP_decay_vertex_z",          gLLP_decay_vertex_z,          "gLLP_decay_vertex_z[2]/F");

  //leptons
  tree_->Branch("nLeptons",  &nLeptons, "nLeptons/I");
  tree_->Branch("lepE",      lepE,      "lepE[nLeptons]/F");
  tree_->Branch("lepPt",     lepPt,     "lepPt[nLeptons]/F");
  tree_->Branch("lepEta",    lepEta,    "lepEta[nLeptons]/F");
  tree_->Branch("lepPhi",    lepPhi,    "lepPhi[nLeptons]/F");
  tree_->Branch("lepPdgId",  lepPdgId,  "lepPdgId[nLeptons]/I");
  tree_->Branch("lepDZ",     lepDZ,     "lepDZ[nLeptons]/F");
  tree_->Branch("lepPassId", lepPassId, "lepPassId[nLeptons]/O");
  tree_->Branch("lepPassLooseIso", lepPassLooseIso, "lepPassLooseIso[nLeptons]/O");
  tree_->Branch("lepPassTightIso", lepPassTightIso, "lepPassTightIso[nLeptons]/O");
  tree_->Branch("lepPassVTightIso", lepPassVTightIso, "lepPassVTightIso[nLeptons]/O");
  tree_->Branch("lepPassVVTightIso", lepPassVVTightIso, "lepPassVVTightIso[nLeptons]/O");

  // tree_->Branch("lepPassVetoId", lepPassVetoId, "lepPassVetoId[nLeptons]/O");

  // tree_->Branch("lepLoosePassId", lepLoosePassId, "lepLoosePassId[nLeptons]/O");
  // tree_->Branch("lepMediumPassId", lepMediumPassId, "lepMediumPassId[nLeptons]/O");
  // tree_->Branch("lepTightPassId", lepTightPassId, "lepTightPassId[nLeptons]/O");


  tree_->Branch("MT",      &MT,  "MT/F");
  //Z-candidate

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
  tree_->Branch("jet_match_genJet_pt", jet_match_genJet_pt, "jet_match_genJet_pt[nJets]/F");
  tree_->Branch("jet_match_genJet_index", jet_match_genJet_index, "jet_match_genJet_index[nJets]/I");
  tree_->Branch("jet_match_genJet_minDeltaR", jet_match_genJet_minDeltaR, "jet_match_genJet_minDeltaR[nJets]/F");

  // tree_->Branch("ecalNRechits",   ecalNRechits,   "ecalNRechits[nJets]/F");
  // tree_->Branch("ecalRechitE", ecalRechitE, "ecalRechitE[nJets]/F");
  // tree_->Branch("jetLoosePassId", jetLoosePassId, "jetLoosePassId[nJets]/O");
  tree_->Branch("jetTightPassId", jetTightPassId, "jetTightPassId[nJets]/O");
  tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[601]/O"); //hardcoded
  tree_->Branch("jetChargedEMEnergyFraction",   jetChargedEMEnergyFraction,   "jetChargedEMEnergyFraction[nJets]/F");
  tree_->Branch("jetNeutralEMEnergyFraction",   jetNeutralEMEnergyFraction,   "jetNeutralEMEnergyFraction[nJets]/F");
  tree_->Branch("jetChargedHadronEnergyFraction",   jetChargedHadronEnergyFraction,   "jetChargedHadronEnergyFraction[nJets]/F");
  tree_->Branch("jetNeutralHadronEnergyFraction",   jetNeutralHadronEnergyFraction,   "jetNeutralHadronEnergyFraction[nJets]/F");
};
