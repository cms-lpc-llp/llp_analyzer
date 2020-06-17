#include "RazorHelper.h"
#include "LiteTreeMuonSystem.h"
#include "assert.h"
#include "TTree.h"

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
  nCsc = 0;
  nCscSegClusters = 0;
  nCscRechitClusters = 0;
  nCscClusters = 0;
  nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto = 0;
  nDt = 0;
  nRpc = 0;
  nHORechits = 0;
  // nCscSegClusters = 0;
  // nCscRechitClusters = 0;
  nCscClusters = 0;
  nDtClusters = 0;
  nDtRechitClusters = 0;
  nDtRechits = 0;
  // nCscClusters = 0;
  // nCscITClusters = 0;
  // nCsc_JetVetoITCluster0p4 = 0;
  // nCsc_JetMuonVetoITCluster0p4 = 0;
  // nCsc_JetVetoITCluster0p4_Me1112Veto = 0;
  // nCsc_JetMuonVetoITCluster0p4_Me1112Veto = 0;
  // nCsc_JetVetoCluster0p4 = 0;
  // nCsc_JetMuonVetoCluster0p4 = 0;
  // nCsc_JetVetoCluster0p4_Me1112Veto = 0;
  // nCsc_JetMuonVetoCluster0p4_Me1112Veto = 0;
  for( int i = 0; i < N_MAX_CSC; i++ )
  {
    cscLabels[i] = -999;
    cscITLabels[i] = -999;
    cscStation[i] = -999;
    cscChamber[i] = -999;

    cscPhi[i] = -999;   //[nCsc]
    cscEta[i] = -999;   //[nCsc]
    cscX[i] = -999;   //[nCsc]
    cscY[i] = -999;   //[nCsc]
    cscZ[i] = -999;   //[nCsc]
    cscDirectionX[i] = -999;   //[nCsc]
    cscDirectionY[i] = -999;   //[nCsc]
    cscDirectionZ[i] = -999;   //[nCsc]
    cscNRecHits[i] = -999;   //[nCsc]
    cscNRecHits_flag[i] = -999;   //[nCsc]
    cscT[i] = -999;   //[nCsc]
    cscChi2[i] = -999;   //[nCsc]

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
    cscSegClusterSize[i] = -999;
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
    cscSegClusterMet_dPhi[i] = -999.;



    cscRechitCluster_match_gLLP[i] = false;
    cscRechitCluster_match_gLLP_minDeltaR[i] = 999;
    cscRechitCluster_match_gLLP_index[i] = 999;
    cscRechitClusterSize[i] = -999;
    cscRechitClusterX[i] = -999.;
    cscRechitClusterY[i] = -999.;
    cscRechitClusterZ[i] = -999.;
    cscRechitClusterTime[i] = -999.;
    cscRechitClusterGenMuonDeltaR[i] = 999.;
    cscRechitClusterTimeSpread[i] = -999.;
    cscRechitClusterMajorAxis[i] = -999.;
    cscRechitClusterMinorAxis[i] = -999.;
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
    cscRechitClusterNChamber[i] = -999;
    cscRechitClusterMaxChamberRatio[i] = -999.;
    cscRechitClusterMaxChamber[i] = -999;
    cscRechitClusterNStation[i] = -999;
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
  }
  for( int i = 0; i < N_MAX_RPC; i++ )
  {
    rpcX[i] = -999;
    rpcY[i] = -999;
    rpcZ[i] = -999;
    rpcEta[i] = -999;
    rpcPhi[i] = -999;
    rpcBx[i] = -999;
  }
  for( int i = 0; i < N_MAX_HO; i++ )
  {
    hoRechit_X[i] = -999;
    hoRechit_Y[i] = -999;
    hoRechit_Z[i] = -999;
    hoRechit_Eta[i] = -999;
    hoRechit_Phi[i] = -999;
    hoRechit_T[i] = -999;
    hoRechit_E[i] = -999;
  }
  for( int i = 0; i < N_MAX_DT; i++ )
  {
    dtLabels[i] = -999;
    dtITLabels[i] = -999;
    dtStation[i] = -999;
    dtChamber[i] = -999;

    dtPhi[i] = -999;   //[nDt]
    dtEta[i] = -999;   //[nDt]
    dtX[i] = -999;   //[nDt]
    dtY[i] = -999;   //[nDt]
    dtZ[i] = -999;   //[nDt]
    dtDirectionX[i] = -999;   //[nDt]
    dtDirectionY[i] = -999;   //[nDt]
    dtDirectionZ[i] = -999;   //[nDt]
    dtNRecHits[i] = -999;   //[nDt]
    dtNRecHits_flag[i] = -999;   //[nDt]
    dtT[i] = -999;   //[nDt]
    dtChi2[i] = -999;   //[nDt]

    dtRechitStation[i] = -999;
    dtRechitWheel[i] = -999;

    dtRechitPhi[i] = -999;   //[nDtRechits]
    dtRechitEta[i] = -999;   //[nDtRechits]
    dtRechitX[i] = -999;   //[nDtRechits]
    dtRechitY[i] = -999;   //[nDtRechits]
    dtRechitZ[i] = -999;   //[nDtRechits]

    dtCluster_match_gLLP[i] = false;
    dtCluster_match_gLLP_minDeltaR[i] = 999;
    dtCluster_match_gLLP_index[i] = 999;
    dtClusterSize[i] = -999;
    dtClusterX[i] = -999.;
    dtClusterY[i] = -999.;
    dtClusterZ[i] = -999.;
    dtClusterTime[i] = -999.;
    dtClusterGenMuonDeltaR[i] = 999.;
    dtClusterTimeSpread[i] = -999.;
    dtClusterMajorAxis[i] = -999.;
    dtClusterMinorAxis[i] = -999.;
    dtClusterXSpread[i] = -999.;
    dtClusterYSpread[i] = -999.;
    dtClusterZSpread[i] = -999.;
    dtClusterEtaPhiSpread[i] = -999.;
    dtClusterEtaSpread[i] = -999.;
    dtClusterPhiSpread[i] = -999.;
    dtClusterEta[i] = -999.;
    dtClusterPhi[i] = -999.;
    dtClusterJetVetoPt[i] = 0.0;
    // dtClusterCaloJetVeto[i] = 0.0;
    dtClusterMuonVetoPt[i] = 0.0;
    dtClusterJetVetoE[i] = 0.0;
    // dtClusterCaloJetVetoE[i] = 0.0;
    dtClusterMuonVetoE[i] = 0.0;
    dtClusterNChamber[i] = -999;
    dtClusterMaxChamberRatio[i] = -999.;
    dtClusterMaxChamber[i] = -999;
    dtClusterNStation[i] = -999;
    dtClusterMaxStationRatio[i] = -999.;
    dtClusterMaxStation[i] = -999;
    dtClusterMe11Ratio[i] = -999.;
    dtClusterMe12Ratio[i] = -999.;
    dtClusterNSegmentChamberPlus11[i] = -999;
    dtClusterNSegmentChamberPlus12[i] = -999;
    dtClusterNSegmentChamberPlus13[i] = -999;
    dtClusterNSegmentChamberPlus21[i] = -999;
    dtClusterNSegmentChamberPlus22[i] = -999;
    dtClusterNSegmentChamberPlus31[i] = -999;
    dtClusterNSegmentChamberPlus32[i] = -999;
    dtClusterNSegmentChamberPlus41[i] = -999;
    dtClusterNSegmentChamberPlus42[i] = -999;
    dtClusterNSegmentChamberMinus11[i] = -999;
    dtClusterNSegmentChamberMinus12[i] = -999;
    dtClusterNSegmentChamberMinus13[i] = -999;
    dtClusterNSegmentChamberMinus21[i] = -999;
    dtClusterNSegmentChamberMinus22[i] = -999;
    dtClusterNSegmentChamberMinus31[i] = -999;
    dtClusterNSegmentChamberMinus32[i] = -999;
    dtClusterNSegmentChamberMinus41[i] = -999;
    dtClusterNSegmentChamberMinus42[i] = -999;
    /*dtSegClusterSize[i] = -999;
    dtSegClusterX[i] = -999.;
    dtSegClusterY[i] = -999.;
    dtSegClusterZ[i] = -999.;
    dtSegClusterTime[i] = -999.;
    dtSegClusterGenMuonDeltaR[i] = 999.;
    dtSegClusterTimeSpread[i] = -999.;
    dtSegClusterMajorAxis[i] = -999.;
    dtSegClusterMinorAxis[i] = -999.;
    dtSegClusterXSpread[i] = -999.;
    dtSegClusterYSpread[i] = -999.;
    dtSegClusterZSpread[i] = -999.;
    dtSegClusterEtaPhiSpread[i] = -999.;
    dtSegClusterEtaSpread[i] = -999.;
    dtSegClusterPhiSpread[i] = -999.;
    dtSegClusterEta[i] = -999.;
    dtSegClusterPhi[i] = -999.;
    dtSegClusterJetVetoPt[i] = 0.0;
    // dtSegClusterCaloJetVeto[i] = 0.0;
    dtSegClusterMuonVetoPt[i] = 0.0;
    dtSegClusterJetVetoE[i] = 0.0;
    // dtSegClusterCaloJetVetoE[i] = 0.0;
    dtSegClusterMuonVetoE[i] = 0.0;
    dtSegClusterNChamber[i] = -999;
    dtSegClusterMaxChamberRatio[i] = -999.;
    dtSegClusterMaxChamber[i] = -999;
    dtSegClusterNStation[i] = -999;
    dtSegClusterMaxStationRatio[i] = -999.;
    dtSegClusterMaxStation[i] = -999;
    dtSegClusterMe11Ratio[i] = -999.;
    dtSegClusterMe12Ratio[i] = -999.;
    dtSegClusterNSegmentChamberPlus11[i] = -999;
    dtSegClusterNSegmentChamberPlus12[i] = -999;
    dtSegClusterNSegmentChamberPlus13[i] = -999;
    dtSegClusterNSegmentChamberPlus21[i] = -999;
    dtSegClusterNSegmentChamberPlus22[i] = -999;
    dtSegClusterNSegmentChamberPlus31[i] = -999;
    dtSegClusterNSegmentChamberPlus32[i] = -999;
    dtSegClusterNSegmentChamberPlus41[i] = -999;
    dtSegClusterNSegmentChamberPlus42[i] = -999;
    dtSegClusterNSegmentChamberMinus11[i] = -999;
    dtSegClusterNSegmentChamberMinus12[i] = -999;
    dtSegClusterNSegmentChamberMinus13[i] = -999;
    dtSegClusterNSegmentChamberMinus21[i] = -999;
    dtSegClusterNSegmentChamberMinus22[i] = -999;
    dtSegClusterNSegmentChamberMinus31[i] = -999;
    dtSegClusterNSegmentChamberMinus32[i] = -999;
    dtSegClusterNSegmentChamberMinus41[i] = -999;
    dtSegClusterNSegmentChamberMinus42[i] = -999;
    // dtClusterVertexR[i] = 0.0;
    // dtClusterVertexZ[i] = 0.0;
    // dtClusterVertexDis[i] = 0.0;
    // dtClusterVertexChi2[i] = 0.0;
    // dtClusterVertexN[i] = 0;
    // dtClusterVertexN1[i] = 0;
    // dtClusterVertexN5[i] = 0;
    // dtClusterVertexN10[i] = 0;
    // dtClusterVertexN15[i] = 0;
    // dtClusterVertexN20[i] = 0;
    */
    dtRechitCluster_match_gParticle_id[i] = -999;
    dtRechitCluster_match_gParticle_index[i] = -999;
    dtRechitCluster_match_gParticle_minDeltaR[i] = -999.;
    dtRechitClusterSize[i] = -999;
    dtRechitClusterX[i] = -999.;
    dtRechitClusterY[i] = -999.;
    dtRechitClusterZ[i] = -999.;
    dtRechitClusterTime[i] = -999.;
    dtRechitClusterGenMuonDeltaR[i] = 999.;
    dtRechitClusterTimeSpread[i] = -999.;
    dtRechitClusterMajorAxis[i] = -999.;
    dtRechitClusterMinorAxis[i] = -999.;
    dtRechitClusterXSpread[i] = -999.;
    dtRechitClusterYSpread[i] = -999.;
    dtRechitClusterZSpread[i] = -999.;
    dtRechitClusterEtaPhiSpread[i] = -999.;
    dtRechitClusterEtaSpread[i] = -999.;
    dtRechitClusterPhiSpread[i] = -999.;
    dtRechitClusterEta[i] = -999.;
    dtRechitClusterPhi[i] = -999.;
    dtRechitClusterJetVetoPt[i] = 0.0;
    dtRechitClusterCaloJetVeto[i] = 0.0;
    dtRechitClusterMuonVetoPt[i] = 0.0;
    dtRechitClusterJetVetoE[i] = 0.0;
    // dtRechitClusterCaloJetVetoE[i] = 0.0;
    dtRechitClusterMuonVetoE[i] = 0.0;
    dtRechitClusterNChamber[i] = -999;
    dtRechitClusterMaxChamberRatio[i] = -999.;
    dtRechitClusterMaxChamber[i] = -999;
    dtRechitClusterNStation[i] = -999;
    dtRechitClusterMaxStationRatio[i] = -999.;
    dtRechitClusterMaxStation[i] = -999;
    /*dtRechitClusterMe11Ratio[i] = -999.;
    dtRechitClusterMe12Ratio[i] = -999.;
    dtRechitClusterNRechitChamberPlus11[i] = -999;
    dtRechitClusterNRechitChamberPlus12[i] = -999;
    dtRechitClusterNRechitChamberPlus13[i] = -999;
    dtRechitClusterNRechitChamberPlus21[i] = -999;
    dtRechitClusterNRechitChamberPlus22[i] = -999;
    dtRechitClusterNRechitChamberPlus31[i] = -999;
    dtRechitClusterNRechitChamberPlus32[i] = -999;
    dtRechitClusterNRechitChamberPlus41[i] = -999;
    dtRechitClusterNRechitChamberPlus42[i] = -999;
    dtRechitClusterNRechitChamberMinus11[i] = -999;
    dtRechitClusterNRechitChamberMinus12[i] = -999;
    dtRechitClusterNRechitChamberMinus13[i] = -999;
    dtRechitClusterNRechitChamberMinus21[i] = -999;
    dtRechitClusterNRechitChamberMinus22[i] = -999;
    dtRechitClusterNRechitChamberMinus31[i] = -999;
    dtRechitClusterNRechitChamberMinus32[i] = -999;
    dtRechitClusterNRechitChamberMinus41[i] = -999;
    dtRechitClusterNRechitChamberMinus42[i] = -999;
    */
    dtRechitClusterNSegmentStation1[i] = -999;
    dtRechitClusterNSegmentStation2[i] = -999;
    dtRechitClusterNSegmentStation3[i] = -999;
    dtRechitClusterNSegmentStation4[i] = -999;
  }
  for(int i = 0;i<2;i++)
  {
    gLLP_eta[i] = 0.0;
    gLLP_phi[i] = 0.0;
    gLLP_beta[i] = 0.0;
    gLLP_csc[i] = 0.0;
    gLLP_dt[i] = 0.0;
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
  }
  //muons
  nMuons = 0;
  for( int i = 0; i < N_MAX_MUONS; i++ )
  {
    muonE[i]      = -999.;
    muonPt[i]     = -999.;
    muonEta[i]    = -999.;
    muonPhi[i]    = -999.;
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
    // jetTightPassId[i] = false;
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
  tree_->SetBranchAddress("nCsc",             &nCsc);
  tree_->SetBranchAddress("cscLabels",             cscLabels);
  tree_->SetBranchAddress("cscITLabels",             cscITLabels);
  tree_->SetBranchAddress("cscStation",             cscStation);
  tree_->SetBranchAddress("cscChamber",             cscChamber);

  tree_->SetBranchAddress("cscPhi",           cscPhi);
  tree_->SetBranchAddress("cscEta",           cscEta);
  tree_->SetBranchAddress("cscX",             cscX);
  tree_->SetBranchAddress("cscY",             cscY);
  tree_->SetBranchAddress("cscZ",             cscZ);
  tree_->SetBranchAddress("cscDirectionX",             cscDirectionX);
  tree_->SetBranchAddress("cscDirectionY",             cscDirectionY);
  tree_->SetBranchAddress("cscDirectionZ",             cscDirectionZ);
  tree_->SetBranchAddress("cscNRecHits",      cscNRecHits);
  tree_->SetBranchAddress("cscNRecHits_flag", cscNRecHits_flag);
  tree_->SetBranchAddress("cscT",             cscT);
  tree_->SetBranchAddress("cscChi2",          cscChi2);

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

  tree_->SetBranchAddress("nCscSegClusters",             &nCscSegClusters);
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

  //RPC
  tree_->SetBranchAddress("nRpc",             &nRpc);
  tree_->SetBranchAddress("rpcPhi",           rpcPhi);
  tree_->SetBranchAddress("rpcEta",           rpcEta);
  tree_->SetBranchAddress("rpcX",             rpcX);
  tree_->SetBranchAddress("rpcY",             rpcY);
  tree_->SetBranchAddress("rpcZ",             rpcZ);
  tree_->SetBranchAddress("rpcBx",            rpcBx);

  //HO
  tree_->SetBranchAddress("nHORechits",             &nHORechits);
  tree_->SetBranchAddress("hoRechit_Phi",           hoRechit_Phi);
  tree_->SetBranchAddress("hoRechit_Eta",           hoRechit_Eta);
  tree_->SetBranchAddress("hoRechit_X",             hoRechit_X);
  tree_->SetBranchAddress("hoRechit_Y",             hoRechit_Y);
  tree_->SetBranchAddress("hoRechit_Z",             hoRechit_Z);
  tree_->SetBranchAddress("hoRechit_T",             hoRechit_T);
  tree_->SetBranchAddress("hoRechit_E",             hoRechit_E);

  //DT
  tree_->SetBranchAddress("nDt",             &nDt);
  tree_->SetBranchAddress("dtLabels",             dtLabels);
  tree_->SetBranchAddress("dtITLabels",             dtITLabels);
  tree_->SetBranchAddress("dtStation",             dtStation);
  tree_->SetBranchAddress("dtChamber",             dtChamber);

  tree_->SetBranchAddress("dtPhi",           dtPhi);
  tree_->SetBranchAddress("dtEta",           dtEta);
  tree_->SetBranchAddress("dtX",             dtX);
  tree_->SetBranchAddress("dtY",             dtY);
  tree_->SetBranchAddress("dtZ",             dtZ);
  tree_->SetBranchAddress("dtDirectionX",             dtDirectionX);
  tree_->SetBranchAddress("dtDirectionY",             dtDirectionY);
  tree_->SetBranchAddress("dtDirectionZ",             dtDirectionZ);
  tree_->SetBranchAddress("dtNRecHits",      dtNRecHits);
  tree_->SetBranchAddress("dtNRecHits_flag", dtNRecHits_flag);
  tree_->SetBranchAddress("dtT",             dtT);
  tree_->SetBranchAddress("dtChi2",          dtChi2);

  tree_->SetBranchAddress("nDtRechits",            &nDtRechits);   
  tree_->SetBranchAddress("dtRechitStation",             dtRechitStation);
  tree_->SetBranchAddress("dtRechitWheel",             dtRechitWheel);
  tree_->SetBranchAddress("dtRechitPhi",           dtRechitPhi);
  tree_->SetBranchAddress("dtRechitEta",           dtRechitEta);
  tree_->SetBranchAddress("dtRechitX",             dtRechitX);
  tree_->SetBranchAddress("dtRechitY",             dtRechitY);
  tree_->SetBranchAddress("dtRechitZ",             dtRechitZ);

  tree_->SetBranchAddress("nDtClusters",             &nDtClusters);
  tree_->SetBranchAddress("dtCluster_match_gLLP",             &dtCluster_match_gLLP);
  tree_->SetBranchAddress("dtCluster_match_gLLP_minDeltaR",             &dtCluster_match_gLLP_minDeltaR);
  tree_->SetBranchAddress("dtCluster_match_gLLP_index",             &dtCluster_match_gLLP_index);

  // tree_->SetBranchAddress("nDt_JetVetoCluster0p4",             &nDt_JetVetoCluster0p4);
  // tree_->SetBranchAddress("nDt_JetMuonVetoCluster0p4",             &nDt_JetMuonVetoCluster0p4);
  // tree_->SetBranchAddress("nDt_JetVetoCluster0p4_Me1112Veto",             &nDt_JetVetoCluster0p4_Me1112Veto);
  // tree_->SetBranchAddress("nDt_JetMuonVetoCluster0p4_Me1112Veto",             &nDt_JetMuonVetoCluster0p4_Me1112Veto);
  tree_->SetBranchAddress("dtClusterMe11Ratio",             &dtClusterMe11Ratio);
  tree_->SetBranchAddress("dtClusterMe12Ratio",             &dtClusterMe12Ratio);
  tree_->SetBranchAddress("dtClusterMaxStation",             &dtClusterMaxStation);
  tree_->SetBranchAddress("dtClusterMaxStationRatio",             &dtClusterMaxStationRatio);
  tree_->SetBranchAddress("dtClusterNStation",             &dtClusterNStation);
  tree_->SetBranchAddress("dtClusterMaxChamber",             &dtClusterMaxChamber);
  tree_->SetBranchAddress("dtClusterMaxChamberRatio",             &dtClusterMaxChamberRatio);
  tree_->SetBranchAddress("dtClusterNChamber",             &dtClusterNChamber);
  tree_->SetBranchAddress("dtClusterX",             dtClusterX);
  tree_->SetBranchAddress("dtClusterY",             dtClusterY);
  tree_->SetBranchAddress("dtClusterZ",             dtClusterZ);
  tree_->SetBranchAddress("dtClusterTime",             dtClusterTime);
  tree_->SetBranchAddress("dtClusterGenMuonDeltaR",             dtClusterGenMuonDeltaR);

  tree_->SetBranchAddress("dtClusterTimeSpread",             dtClusterTimeSpread);
  tree_->SetBranchAddress("dtClusterMajorAxis",             dtClusterMajorAxis);
  tree_->SetBranchAddress("dtClusterMinorAxis",             dtClusterMinorAxis);
  tree_->SetBranchAddress("dtClusterXSpread",             dtClusterXSpread);
  tree_->SetBranchAddress("dtClusterYSpread",             dtClusterYSpread);
  tree_->SetBranchAddress("dtClusterZSpread",             dtClusterZSpread);
  tree_->SetBranchAddress("dtClusterEtaPhiSpread",             dtClusterEtaPhiSpread);
  tree_->SetBranchAddress("dtClusterEtaSpread",             dtClusterEtaSpread);
  tree_->SetBranchAddress("dtClusterPhiSpread",             dtClusterPhiSpread);
  tree_->SetBranchAddress("dtClusterEta",             dtClusterEta);
  tree_->SetBranchAddress("dtClusterPhi",             dtClusterPhi);
  tree_->SetBranchAddress("dtClusterJetVetoPt",             dtClusterJetVetoPt);
  tree_->SetBranchAddress("dtClusterMuonVetoPt",             dtClusterMuonVetoPt);
  tree_->SetBranchAddress("dtClusterJetVetoE",             dtClusterJetVetoE);
  tree_->SetBranchAddress("dtClusterMuonVetoE",             dtClusterMuonVetoE);
  tree_->SetBranchAddress("dtClusterSize",             dtClusterSize);

  tree_->SetBranchAddress("dtClusterNSegmentChamberPlus11",             dtClusterNSegmentChamberPlus11);
  tree_->SetBranchAddress("dtClusterNSegmentChamberPlus12",             dtClusterNSegmentChamberPlus12);
  tree_->SetBranchAddress("dtClusterNSegmentChamberPlus13",             dtClusterNSegmentChamberPlus13);
  tree_->SetBranchAddress("dtClusterNSegmentChamberPlus21",             dtClusterNSegmentChamberPlus21);
  tree_->SetBranchAddress("dtClusterNSegmentChamberPlus22",             dtClusterNSegmentChamberPlus22);
  tree_->SetBranchAddress("dtClusterNSegmentChamberPlus31",             dtClusterNSegmentChamberPlus31);
  tree_->SetBranchAddress("dtClusterNSegmentChamberPlus32",             dtClusterNSegmentChamberPlus32);
  tree_->SetBranchAddress("dtClusterNSegmentChamberPlus41",             dtClusterNSegmentChamberPlus41);
  tree_->SetBranchAddress("dtClusterNSegmentChamberPlus42",             dtClusterNSegmentChamberPlus42);

  tree_->SetBranchAddress("dtClusterNSegmentChamberMinus11",             dtClusterNSegmentChamberMinus11);
  tree_->SetBranchAddress("dtClusterNSegmentChamberMinus12",             dtClusterNSegmentChamberMinus12);
  tree_->SetBranchAddress("dtClusterNSegmentChamberMinus13",             dtClusterNSegmentChamberMinus13);
  tree_->SetBranchAddress("dtClusterNSegmentChamberMinus21",             dtClusterNSegmentChamberMinus21);
  tree_->SetBranchAddress("dtClusterNSegmentChamberMinus22",             dtClusterNSegmentChamberMinus22);
  tree_->SetBranchAddress("dtClusterNSegmentChamberMinus31",             dtClusterNSegmentChamberMinus31);
  tree_->SetBranchAddress("dtClusterNSegmentChamberMinus32",             dtClusterNSegmentChamberMinus32);
  tree_->SetBranchAddress("dtClusterNSegmentChamberMinus41",             dtClusterNSegmentChamberMinus41);
  tree_->SetBranchAddress("dtClusterNSegmentChamberMinus42",             dtClusterNSegmentChamberMinus42);


  // DT CLUSTER
  /*tree_->SetBranchAddress("nDtSegClusters",             &nDtSegClusters);
  // tree_->SetBranchAddress("nDt_JetVetoSegCluster0p4",             &nDt_JetVetoSegCluster0p4);
  // tree_->SetBranchAddress("nDt_JetMuonVetoSegCluster0p4",             &nDt_JetMuonVetoSegCluster0p4);
  // tree_->SetBranchAddress("nDt_JetVetoSegCluster0p4_Me1112Veto",             &nDt_JetVetoSegCluster0p4_Me1112Veto);
  // tree_->SetBranchAddress("nDt_JetMuonVetoSegCluster0p4_Me1112Veto",             &nDt_JetMuonVetoCluster0p4_Me1112Veto);
  tree_->SetBranchAddress("dtSegClusterMe11Ratio",             &dtSegClusterMe11Ratio);
  tree_->SetBranchAddress("dtSegClusterMe12Ratio",             &dtSegClusterMe12Ratio);
  tree_->SetBranchAddress("dtSegClusterMaxStation",             &dtSegClusterMaxStation);
  tree_->SetBranchAddress("dtSegClusterMaxStationRatio",             &dtSegClusterMaxStationRatio);
  tree_->SetBranchAddress("dtSegClusterNStation",             &dtSegClusterNStation);
  tree_->SetBranchAddress("dtSegClusterMaxChamber",             &dtSegClusterMaxChamber);
  tree_->SetBranchAddress("dtSegClusterMaxChamberRatio",             &dtSegClusterMaxChamberRatio);
  tree_->SetBranchAddress("dtSegClusterNChamber",             &dtSegClusterNChamber);
  tree_->SetBranchAddress("dtSegClusterX",             dtSegClusterX);
  tree_->SetBranchAddress("dtSegClusterY",             dtSegClusterY);
  tree_->SetBranchAddress("dtSegClusterZ",             dtSegClusterZ);
  tree_->SetBranchAddress("dtSegClusterTime",             dtSegClusterTime);
  tree_->SetBranchAddress("dtSegClusterGenMuonDeltaR",             dtSegClusterGenMuonDeltaR);

  tree_->SetBranchAddress("dtSegClusterTimeSpread",             dtSegClusterTimeSpread);
  tree_->SetBranchAddress("dtSegClusterMajorAxis",             dtSegClusterMajorAxis);
  tree_->SetBranchAddress("dtSegClusterMinorAxis",             dtSegClusterMinorAxis);
  tree_->SetBranchAddress("dtSegClusterXSpread",             dtSegClusterXSpread);
  tree_->SetBranchAddress("dtSegClusterYSpread",             dtSegClusterYSpread);
  tree_->SetBranchAddress("dtSegClusterZSpread",             dtSegClusterZSpread);
  tree_->SetBranchAddress("dtSegClusterEtaPhiSpread",             dtSegClusterEtaPhiSpread);
  tree_->SetBranchAddress("dtSegClusterEtaSpread",             dtSegClusterEtaSpread);
  tree_->SetBranchAddress("dtSegClusterPhiSpread",             dtSegClusterPhiSpread);
  tree_->SetBranchAddress("dtSegClusterEta",             dtSegClusterEta);
  tree_->SetBranchAddress("dtSegClusterPhi",             dtSegClusterPhi);
  tree_->SetBranchAddress("dtSegClusterJetVetoPt",             dtSegClusterJetVetoPt);
  tree_->SetBranchAddress("dtSegClusterMuonVetoPt",             dtSegClusterMuonVetoPt);
  tree_->SetBranchAddress("dtSegClusterJetVetoE",             dtSegClusterJetVetoE);
  tree_->SetBranchAddress("dtSegClusterMuonVetoE",             dtSegClusterMuonVetoE);
  tree_->SetBranchAddress("dtSegClusterSize",             dtSegClusterSize);

  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus11",             dtSegClusterNSegmentChamberPlus11);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus12",             dtSegClusterNSegmentChamberPlus12);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus13",             dtSegClusterNSegmentChamberPlus13);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus21",             dtSegClusterNSegmentChamberPlus21);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus22",             dtSegClusterNSegmentChamberPlus22);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus31",             dtSegClusterNSegmentChamberPlus31);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus32",             dtSegClusterNSegmentChamberPlus32);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus41",             dtSegClusterNSegmentChamberPlus41);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus42",             dtSegClusterNSegmentChamberPlus42);

  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus11",             dtSegClusterNSegmentChamberMinus11);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus12",             dtSegClusterNSegmentChamberMinus12);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus13",             dtSegClusterNSegmentChamberMinus13);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus21",             dtSegClusterNSegmentChamberMinus21);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus22",             dtSegClusterNSegmentChamberMinus22);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus31",             dtSegClusterNSegmentChamberMinus31);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus32",             dtSegClusterNSegmentChamberMinus32);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus41",             dtSegClusterNSegmentChamberMinus41);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus42",             dtSegClusterNSegmentChamberMinus42);
  */

  // DT CLUSTER
  tree_->SetBranchAddress("nDtRechitClusters",             &nDtRechitClusters);
  tree_->SetBranchAddress("dtRechitCluster_match_gParticle_id",             &dtRechitCluster_match_gParticle_id);
  tree_->SetBranchAddress("dtRechitCluster_match_gParticle_index",             &dtRechitCluster_match_gParticle_index);
  tree_->SetBranchAddress("dtRechitCluster_match_gParticle_minDeltaR",             &dtRechitCluster_match_gParticle_minDeltaR);
  //tree_->SetBranchAddress("nDt_JetVetoRechitCluster0p4",             &nDt_JetVetoRechitCluster0p4);
  //tree_->SetBranchAddress("nDt_JetMuonVetoRechitCluster0p4",             &nDt_JetMuonVetoRechitCluster0p4);
  //tree_->SetBranchAddress("nDt_JetVetoRechitCluster0p4_Me1112Veto",             &nDt_JetVetoRechitCluster0p4_Me1112Veto);
  //tree_->SetBranchAddress("nDt_JetMuonVetoRechitCluster0p4_Me1112Veto",             &nDt_JetMuonVetoRechitCluster0p4_Me1112Veto);
  //tree_->SetBranchAddress("dtRechitClusterMe11Ratio",             &dtRechitClusterMe11Ratio);
  //tree_->SetBranchAddress("dtRechitClusterMe12Ratio",             &dtRechitClusterMe12Ratio);
  tree_->SetBranchAddress("dtRechitClusterMaxStation",             &dtRechitClusterMaxStation);
  tree_->SetBranchAddress("dtRechitClusterMaxStationRatio",             &dtRechitClusterMaxStationRatio);
  tree_->SetBranchAddress("dtRechitClusterNStation",             &dtRechitClusterNStation);
  tree_->SetBranchAddress("dtRechitClusterMaxChamber",             &dtRechitClusterMaxChamber);
  tree_->SetBranchAddress("dtRechitClusterMaxChamberRatio",             &dtRechitClusterMaxChamberRatio);
  tree_->SetBranchAddress("dtRechitClusterNChamber",             &dtRechitClusterNChamber);
  tree_->SetBranchAddress("dtRechitClusterX",             dtRechitClusterX);
  tree_->SetBranchAddress("dtRechitClusterY",             dtRechitClusterY);
  tree_->SetBranchAddress("dtRechitClusterZ",             dtRechitClusterZ);
  tree_->SetBranchAddress("dtRechitClusterTime",             dtRechitClusterTime);
  tree_->SetBranchAddress("dtRechitClusterGenMuonDeltaR",             dtRechitClusterGenMuonDeltaR);

  tree_->SetBranchAddress("dtRechitClusterTimeSpread",             dtRechitClusterTimeSpread);
  tree_->SetBranchAddress("dtRechitClusterMajorAxis",             dtRechitClusterMajorAxis);
  tree_->SetBranchAddress("dtRechitClusterMinorAxis",             dtRechitClusterMinorAxis);
  tree_->SetBranchAddress("dtRechitClusterXSpread",             dtRechitClusterXSpread);
  tree_->SetBranchAddress("dtRechitClusterYSpread",             dtRechitClusterYSpread);
  tree_->SetBranchAddress("dtRechitClusterZSpread",             dtRechitClusterZSpread);
  tree_->SetBranchAddress("dtRechitClusterEtaPhiSpread",             dtRechitClusterEtaPhiSpread);
  tree_->SetBranchAddress("dtRechitClusterEtaSpread",             dtRechitClusterEtaSpread);
  tree_->SetBranchAddress("dtRechitClusterPhiSpread",             dtRechitClusterPhiSpread);
  tree_->SetBranchAddress("dtRechitClusterEta",             dtRechitClusterEta);
  tree_->SetBranchAddress("dtRechitClusterPhi",             dtRechitClusterPhi);
  tree_->SetBranchAddress("dtRechitClusterJetVetoPt",             dtRechitClusterJetVetoPt);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoPt",             dtRechitClusterMuonVetoPt);
  tree_->SetBranchAddress("dtRechitClusterJetVetoE",             dtRechitClusterJetVetoE);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoE",             dtRechitClusterMuonVetoE);
  tree_->SetBranchAddress("dtRechitClusterSize",             dtRechitClusterSize);
  tree_->SetBranchAddress("dtRechitClusterNSegmentStation1",             dtRechitClusterNSegmentStation1); 
  tree_->SetBranchAddress("dtRechitClusterNSegmentStation2",             dtRechitClusterNSegmentStation2); 
  tree_->SetBranchAddress("dtRechitClusterNSegmentStation3",             dtRechitClusterNSegmentStation3); 
  tree_->SetBranchAddress("dtRechitClusterNSegmentStation4",             dtRechitClusterNSegmentStation4); 
  /*tree_->SetBranchAddress("dtRechitClusterNRechitChamberPlus11",             dtRechitClusterNRechitChamberPlus11);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberPlus12",             dtRechitClusterNRechitChamberPlus12);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberPlus13",             dtRechitClusterNRechitChamberPlus13);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberPlus21",             dtRechitClusterNRechitChamberPlus21);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberPlus22",             dtRechitClusterNRechitChamberPlus22);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberPlus31",             dtRechitClusterNRechitChamberPlus31);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberPlus32",             dtRechitClusterNRechitChamberPlus32);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberPlus41",             dtRechitClusterNRechitChamberPlus41);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberPlus42",             dtRechitClusterNRechitChamberPlus42);

  tree_->SetBranchAddress("dtRechitClusterNRechitChamberMinus11",             dtRechitClusterNRechitChamberMinus11);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberMinus12",             dtRechitClusterNRechitChamberMinus12);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberMinus13",             dtRechitClusterNRechitChamberMinus13);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberMinus21",             dtRechitClusterNRechitChamberMinus21);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberMinus22",             dtRechitClusterNRechitChamberMinus22);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberMinus31",             dtRechitClusterNRechitChamberMinus31);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberMinus32",             dtRechitClusterNRechitChamberMinus32);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberMinus41",             dtRechitClusterNRechitChamberMinus41);
  tree_->SetBranchAddress("dtRechitClusterNRechitChamberMinus42",             dtRechitClusterNRechitChamberMinus42);
  */
  // DT CLUSTER
  /*tree_->SetBranchAddress("nDtSegClusters",             &nDtSegClusters);
  // tree_->SetBranchAddress("nDt_JetVetoSegCluster0p4",             &nDt_JetVetoSegCluster0p4);
  // tree_->SetBranchAddress("nDt_JetMuonVetoSegCluster0p4",             &nDt_JetMuonVetoSegCluster0p4);
  // tree_->SetBranchAddress("nDt_JetVetoSegCluster0p4_Me1112Veto",             &nDt_JetVetoSegCluster0p4_Me1112Veto);
  // tree_->SetBranchAddress("nDt_JetMuonVetoSegCluster0p4_Me1112Veto",             &nDt_JetMuonVetoCluster0p4_Me1112Veto);
  tree_->SetBranchAddress("dtSegClusterMe11Ratio",             &dtSegClusterMe11Ratio);
  tree_->SetBranchAddress("dtSegClusterMe12Ratio",             &dtSegClusterMe12Ratio);
  tree_->SetBranchAddress("dtSegClusterMaxStation",             &dtSegClusterMaxStation);
  tree_->SetBranchAddress("dtSegClusterMaxStationRatio",             &dtSegClusterMaxStationRatio);
  tree_->SetBranchAddress("dtSegClusterNStation",             &dtSegClusterNStation);
  tree_->SetBranchAddress("dtSegClusterMaxChamber",             &dtSegClusterMaxChamber);
  tree_->SetBranchAddress("dtSegClusterMaxChamberRatio",             &dtSegClusterMaxChamberRatio);
  tree_->SetBranchAddress("dtSegClusterNChamber",             &dtSegClusterNChamber);
  tree_->SetBranchAddress("dtSegClusterX",             dtSegClusterX);
  tree_->SetBranchAddress("dtSegClusterY",             dtSegClusterY);
  tree_->SetBranchAddress("dtSegClusterZ",             dtSegClusterZ);
  tree_->SetBranchAddress("dtSegClusterTime",             dtSegClusterTime);
  tree_->SetBranchAddress("dtSegClusterGenMuonDeltaR",             dtSegClusterGenMuonDeltaR);

  tree_->SetBranchAddress("dtSegClusterTimeSpread",             dtSegClusterTimeSpread);
  tree_->SetBranchAddress("dtSegClusterMajorAxis",             dtSegClusterMajorAxis);
  tree_->SetBranchAddress("dtSegClusterMinorAxis",             dtSegClusterMinorAxis);
  tree_->SetBranchAddress("dtSegClusterXSpread",             dtSegClusterXSpread);
  tree_->SetBranchAddress("dtSegClusterYSpread",             dtSegClusterYSpread);
  tree_->SetBranchAddress("dtSegClusterZSpread",             dtSegClusterZSpread);
  tree_->SetBranchAddress("dtSegClusterEtaPhiSpread",             dtSegClusterEtaPhiSpread);
  tree_->SetBranchAddress("dtSegClusterEtaSpread",             dtSegClusterEtaSpread);
  tree_->SetBranchAddress("dtSegClusterPhiSpread",             dtSegClusterPhiSpread);
  tree_->SetBranchAddress("dtSegClusterEta",             dtSegClusterEta);
  tree_->SetBranchAddress("dtSegClusterPhi",             dtSegClusterPhi);
  tree_->SetBranchAddress("dtSegClusterJetVetoPt",             dtSegClusterJetVetoPt);
  tree_->SetBranchAddress("dtSegClusterMuonVetoPt",             dtSegClusterMuonVetoPt);
  tree_->SetBranchAddress("dtSegClusterJetVetoE",             dtSegClusterJetVetoE);
  tree_->SetBranchAddress("dtSegClusterMuonVetoE",             dtSegClusterMuonVetoE);
  tree_->SetBranchAddress("dtSegClusterSize",             dtSegClusterSize);

  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus11",             dtSegClusterNSegmentChamberPlus11);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus12",             dtSegClusterNSegmentChamberPlus12);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus13",             dtSegClusterNSegmentChamberPlus13);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus21",             dtSegClusterNSegmentChamberPlus21);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus22",             dtSegClusterNSegmentChamberPlus22);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus31",             dtSegClusterNSegmentChamberPlus31);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus32",             dtSegClusterNSegmentChamberPlus32);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus41",             dtSegClusterNSegmentChamberPlus41);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberPlus42",             dtSegClusterNSegmentChamberPlus42);

  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus11",             dtSegClusterNSegmentChamberMinus11);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus12",             dtSegClusterNSegmentChamberMinus12);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus13",             dtSegClusterNSegmentChamberMinus13);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus21",             dtSegClusterNSegmentChamberMinus21);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus22",             dtSegClusterNSegmentChamberMinus22);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus31",             dtSegClusterNSegmentChamberMinus31);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus32",             dtSegClusterNSegmentChamberMinus32);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus41",             dtSegClusterNSegmentChamberMinus41);
  tree_->SetBranchAddress("dtSegClusterNSegmentChamberMinus42",             dtSegClusterNSegmentChamberMinus42);

  // CSC CLUSTER

  tree_->SetBranchAddress("nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto);
  tree_->SetBranchAddress("nCscRechitClusters",             &nCscRechitClusters);
  tree_->SetBranchAddress("cscRechitCluster_match_gParticle_id",             &cscRechitCluster_match_gParticle_id);

  tree_->SetBranchAddress("cscRechitCluster_match_gLLP",             &cscRechitCluster_match_gLLP);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_index",             &cscRechitCluster_match_gLLP_index);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_minDeltaR",             &cscRechitCluster_match_gLLP_minDeltaR);
  tree_->SetBranchAddress("cscRechitClusterMe11Ratio",             &cscRechitClusterMe11Ratio);
  tree_->SetBranchAddress("cscRechitClusterMe12Ratio",             &cscRechitClusterMe12Ratio);
  tree_->SetBranchAddress("cscRechitClusterMaxStation",             &cscRechitClusterMaxStation);
  tree_->SetBranchAddress("cscRechitClusterMaxStationRatio",             &cscRechitClusterMaxStationRatio);
  tree_->SetBranchAddress("cscRechitClusterNStation",             &cscRechitClusterNStation);
  tree_->SetBranchAddress("cscRechitClusterMaxChamber",             &cscRechitClusterMaxChamber);
  tree_->SetBranchAddress("cscRechitClusterMaxChamberRatio",             &cscRechitClusterMaxChamberRatio);
  tree_->SetBranchAddress("cscRechitClusterNChamber",             &cscRechitClusterNChamber);
  tree_->SetBranchAddress("cscRechitClusterX",             cscRechitClusterX);
  tree_->SetBranchAddress("cscRechitClusterY",             cscRechitClusterY);
  tree_->SetBranchAddress("cscRechitClusterZ",             cscRechitClusterZ);
  tree_->SetBranchAddress("cscRechitClusterTime",             cscRechitClusterTime);
  tree_->SetBranchAddress("cscRechitClusterGenMuonDeltaR",             cscRechitClusterGenMuonDeltaR);

  tree_->SetBranchAddress("cscRechitClusterTimeSpread",             cscRechitClusterTimeSpread);
  tree_->SetBranchAddress("cscRechitClusterMajorAxis",             cscRechitClusterMajorAxis);
  tree_->SetBranchAddress("cscRechitClusterMinorAxis",             cscRechitClusterMinorAxis);
  tree_->SetBranchAddress("cscRechitClusterXSpread",             cscRechitClusterXSpread);
  tree_->SetBranchAddress("cscRechitClusterYSpread",             cscRechitClusterYSpread);
  tree_->SetBranchAddress("cscRechitClusterZSpread",             cscRechitClusterZSpread);
  tree_->SetBranchAddress("cscRechitClusterEtaPhiSpread",             cscRechitClusterEtaPhiSpread);
  tree_->SetBranchAddress("cscRechitClusterEtaSpread",             cscRechitClusterEtaSpread);
  tree_->SetBranchAddress("cscRechitClusterPhiSpread",             cscRechitClusterPhiSpread);
  tree_->SetBranchAddress("cscRechitClusterEta",             cscRechitClusterEta);
  tree_->SetBranchAddress("cscRechitClusterPhi",             cscRechitClusterPhi);
  tree_->SetBranchAddress("cscRechitClusterJetVetoPt",             cscRechitClusterJetVetoPt);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoPt",             cscRechitClusterMuonVetoPt);
  tree_->SetBranchAddress("cscRechitClusterJetVetoE",             cscRechitClusterJetVetoE);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoE",             cscRechitClusterMuonVetoE);
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

  // CSC IT CLUSTER
  tree_->SetBranchAddress("nCscITClusters",             &nCscITClusters);
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
  tree_->SetBranchAddress("gLLP_dt",     gLLP_dt);
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

  //Muons
  tree_->SetBranchAddress("nMuons",    &nMuons);
  tree_->SetBranchAddress("muonE",        muonE);
  tree_->SetBranchAddress("muonPt",       muonPt);
  tree_->SetBranchAddress("muonEta",      muonEta);
  tree_->SetBranchAddress("muonPhi",      muonPhi);

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
  // tree_->SetBranchAddress("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction);
  // tree_->SetBranchAddress("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction);
  // tree_->SetBranchAddress("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction);
  // tree_->SetBranchAddress("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction);

  // tree_->SetBranchAddress("jetLoosePassId", jetLoosePassId);
  // tree_->SetBranchAddress("jetTightPassId", jetTightPassId);
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
  tree_->Branch("mH",      &mH,     "mH/I");
  tree_->Branch("mX",      &mX,     "mX/I");
  tree_->Branch("ctau",      &ctau,     "ctau/I");
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
  tree_->Branch("metXYCorr",      &metXYCorr,     "metXYCorr/F");
  tree_->Branch("metPhiXYCorr",      &metPhiXYCorr,     "metPhiXYCorr/F");

  tree_->Branch("HT",      &HT,     "HT/F");

  tree_->Branch("jetMet_dPhi",      &jetMet_dPhi,     "jetMet_dPhi/F");
  tree_->Branch("jetMet_dPhiMin",      &jetMet_dPhiMin,     "jetMet_dPhiMin/F");
  tree_->Branch("jetMet_dPhiMin4",      &jetMet_dPhiMin4,     "jetMet_dPhiMin4/F");

  tree_->Branch("metJESUp",      &metJESUp,     "metJESUp/F");
  tree_->Branch("metJESDown",      &metJESDown,     "metJESDown/F");

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
  /*
  tree_->Branch("nCsc",             &nCsc, "nCsc/I");
  tree_->Branch("cscITLabels",             cscITLabels,             "cscITLabels[nCsc]/I");

  tree_->Branch("cscLabels",             cscLabels,             "cscLabels[nCsc]/I");
  tree_->Branch("cscStation",             cscStation,             "cscStation[nCsc]/I");
  tree_->Branch("cscChamber",             cscChamber,             "cscChamber[nCsc]/I");

  tree_->Branch("cscPhi",           cscPhi,           "cscPhi[nCsc]/F");
  tree_->Branch("cscEta",           cscEta,           "cscEta[nCsc]/F");
  tree_->Branch("cscX",             cscX,             "cscX[nCsc]/F");
  tree_->Branch("cscY",             cscY,             "cscY[nCsc]/F");
  tree_->Branch("cscZ",             cscZ,             "cscZ[nCsc]/F");
  tree_->Branch("cscDirectionX",             cscDirectionX,             "cscDirectionX[nCsc]/F");
  tree_->Branch("cscDirectionY",             cscDirectionY,             "cscDirectionY[nCsc]/F");
  tree_->Branch("cscDirectionZ",             cscDirectionZ,             "cscDirectionZ[nCsc]/F");
  tree_->Branch("cscNRecHits",      cscNRecHits,      "cscNRecHits[nCsc]/F");
  tree_->Branch("cscNRecHits_flag", cscNRecHits_flag, "cscNRecHits_flag[nCsc]/F");
  tree_->Branch("cscT",             cscT,             "cscT[nCsc]/F");
  tree_->Branch("cscChi2",          cscChi2,          "cscChi2[nCsc]/F");

  tree_->Branch("nCscClusters",             &nCscClusters, "nCscClusters/I");
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

  // all csc SegClusters

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


    // all csc RechitClusters

    tree_->Branch("nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto, "nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto/I");
    tree_->Branch("nCscRechitClusters",             &nCscRechitClusters, "nCscRechitClusters/I");
    tree_->Branch("cscRechitCluster_match_gLLP",             cscRechitCluster_match_gLLP,             "cscRechitCluster_match_gLLP[nCscRechitClusters]/O");
    tree_->Branch("cscRechitCluster_match_gLLP_minDeltaR",             cscRechitCluster_match_gLLP_minDeltaR,             "cscRechitCluster_match_gLLP_minDeltaR[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_index",             cscRechitCluster_match_gLLP_index,             "cscRechitCluster_match_gLLP_index[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_gParticle_id",             cscRechitCluster_match_gParticle_id,             "cscRechitCluster_match_gParticle_id[nCscRechitClusters]/O");

    tree_->Branch("cscRechitClusterX",             cscRechitClusterX,             "cscRechitClusterX[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterY",             cscRechitClusterY,             "cscRechitClusterY[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterZ",             cscRechitClusterZ,             "cscRechitClusterZ[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTime",             cscRechitClusterTime,             "cscRechitClusterTime[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeSpread",             cscRechitClusterTimeSpread,             "cscRechitClusterTimeSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterGenMuonDeltaR",             cscRechitClusterGenMuonDeltaR,             "cscRechitClusterGenMuonDeltaR[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterMajorAxis",             cscRechitClusterMajorAxis,             "cscRechitClusterMajorAxis[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMinorAxis",             cscRechitClusterMinorAxis,             "cscRechitClusterMinorAxis[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEtaPhiSpread",             cscRechitClusterEtaPhiSpread,             "cscRechitClusterEtaPhiSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterPhiSpread",             cscRechitClusterPhiSpread,             "cscRechitClusterPhiSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEtaSpread",             cscRechitClusterEtaSpread,             "cscRechitClusterEtaSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterXSpread",             cscRechitClusterXSpread,             "cscRechitClusterXSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterYSpread",             cscRechitClusterYSpread,             "cscRechitClusterYSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterZSpread",             cscRechitClusterZSpread,             "cscRechitClusterZSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterPhi",             cscRechitClusterPhi,             "cscRechitClusterPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEta",             cscRechitClusterEta,             "cscRechitClusterEta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoPt",             cscRechitClusterJetVetoPt,             "cscRechitClusterJetVetoPt[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoPt",             cscRechitClusterMuonVetoPt,             "cscRechitClusterMuonVetoPt[nCscRechitClusters]/F");
    // tree_->Branch("cscRechitClusterCaloJetVeto",             cscRechitClusterCaloJetVeto,             "cscRechitClusterCaloJetVeto[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoE",             cscRechitClusterJetVetoE,             "cscRechitClusterJetVetoE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoE",             cscRechitClusterMuonVetoE,             "cscRechitClusterMuonVetoE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterSize",             cscRechitClusterSize,             "cscRechitClusterSize[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterMe11Ratio",             cscRechitClusterMe11Ratio,             "cscRechitClusterMe11Ratio[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMe12Ratio",             cscRechitClusterMe12Ratio,             "cscRechitClusterMe12Ratio[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterNStation",             cscRechitClusterNStation,             "cscRechitClusterNStation[nCscRechitClusters]/I");
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
  */
    //RPC
    tree_->Branch("nRpc",             &nRpc,           "nRpc/I");
    tree_->Branch("rpcPhi",           rpcPhi,           "rpcPhi[nRpc]/F");
    tree_->Branch("rpcEta",           rpcEta,           "rpcEta[nRpc]/F");
    tree_->Branch("rpcX",             rpcX,             "rpcX[nRpc]/F");
    tree_->Branch("rpcY",             rpcY,             "rpcY[nRpc]/F");
    tree_->Branch("rpcZ",             rpcZ,             "rpcZ[nRpc]/F");
    tree_->Branch("rpcBx",             rpcBx,             "rpcBx[nRpc]/I");

    //HO
    tree_->Branch("nHORechits",             &nHORechits,           "nHORechits/I");
    tree_->Branch("hoRechit_Phi",           hoRechit_Phi,           "hoRechit_Phi[nHORechits]/F");
    tree_->Branch("hoRechit_Eta",           hoRechit_Eta,           "hoRechit_Eta[nHORechits]/F");
    tree_->Branch("hoRechit_X",             hoRechit_X,             "hoRechit_X[nHORechits]/F");
    tree_->Branch("hoRechit_Y",             hoRechit_Y,             "hoRechit_Y[nHORechits]/F");
    tree_->Branch("hoRechit_Z",             hoRechit_Z,             "hoRechit_Z[nHORechits]/F");
    tree_->Branch("hoRechit_T",             hoRechit_T,             "hoRechit_T[nHORechits]/F");
    tree_->Branch("hoRechit_E",             hoRechit_E,             "hoRechit_E[nHORechits]/F");
  
    //DT

  tree_->Branch("nDt",             &nDt, "nDt/I");
  tree_->Branch("dtITLabels",             dtITLabels,             "dtITLabels[nDt]/I");

  tree_->Branch("dtLabels",             dtLabels,             "dtLabels[nDt]/I");
  tree_->Branch("dtStation",             dtStation,             "dtStation[nDt]/I");
  tree_->Branch("dtChamber",             dtChamber,             "dtChamber[nDt]/I");

  //tree_->Branch("dtPhi",           dtPhi,           "dtPhi[nDt]/F");
  //tree_->Branch("dtEta",           dtEta,           "dtEta[nDt]/F");
  //tree_->Branch("dtX",             dtX,             "dtX[nDt]/F");
  //tree_->Branch("dtY",             dtY,             "dtY[nDt]/F");
  //tree_->Branch("dtZ",             dtZ,             "dtZ[nDt]/F");
  //tree_->Branch("dtDirectionX",             dtDirectionX,             "dtDirectionX[nDt]/F");
  //tree_->Branch("dtDirectionY",             dtDirectionY,             "dtDirectionY[nDt]/F");
  //tree_->Branch("dtDirectionZ",             dtDirectionZ,             "dtDirectionZ[nDt]/F");
  //tree_->Branch("dtNRecHits",      dtNRecHits,      "dtNRecHits[nDt]/F");
  //tree_->Branch("dtNRecHits_flag", dtNRecHits_flag, "dtNRecHits_flag[nDt]/F");
  //tree_->Branch("dtT",             dtT,             "dtT[nDt]/F");
  //tree_->Branch("dtChi2",          dtChi2,          "dtChi2[nDt]/F");

  tree_->Branch("nDtRechits",            &nDtRechits,           "nDtRechits/I");
  tree_->Branch("dtRechitStation",             dtRechitStation,             "dtRechitStation[nDtRechits]/I");
  tree_->Branch("dtRechitWheel",             dtRechitWheel,             "dtRechitWheel[nDtRechits]/I");
  tree_->Branch("dtRechitPhi",           dtRechitPhi,           "dtRechitPhi[nDtRechits]/F");
  tree_->Branch("dtRechitEta",           dtRechitEta,           "dtRechitEta[nDtRechits]/F");
  tree_->Branch("dtRechitX",             dtRechitX,             "dtRechitX[nDtRechits]/F");
  tree_->Branch("dtRechitY",             dtRechitY,             "dtRechitY[nDtRechits]/F");
  tree_->Branch("dtRechitZ",             dtRechitZ,             "dtRechitZ[nDtRechits]/F");

  /*tree_->Branch("nDtClusters",             &nDtClusters, "nDtClusters/I");
  tree_->Branch("dtCluster_match_gLLP",             dtCluster_match_gLLP, "dtCluster_match_gLLP[nDtClusters]/O");
  tree_->Branch("dtCluster_match_gLLP_index",             dtCluster_match_gLLP_index, "dtCluster_match_gLLP_index[nDtClusters]/I");
  tree_->Branch("dtCluster_match_gLLP_minDeltaR",             dtCluster_match_gLLP_minDeltaR, "dtCluster_match_gLLP_minDeltaR[nDtClusters]/I");

  // tree_->Branch("nDt_JetVetoCluster0p4",             &nDt_JetVetoCluster0p4, "nDt_JetVetoCluster0p4/I");
  // tree_->Branch("nDt_JetMuonVetoCluster0p4",             &nDt_JetMuonVetoCluster0p4, "nDt_JetMuonVetoCluster0p4/I");
  // tree_->Branch("nDt_JetVetoCluster0p4_Me1112Veto",             &nDt_JetVetoCluster0p4_Me1112Veto, "nDt_JetVetoCluster0p4_Me1112Veto/I");
  // tree_->Branch("nDt_JetMuonVetoCluster0p4_Me1112Veto",             &nDt_JetMuonVetoCluster0p4_Me1112Veto, "nDt_JetMuonVetoCluster0p4_Me1112Veto/I");
  tree_->Branch("dtClusterX",             dtClusterX,             "dtClusterX[nDtClusters]/F");
  tree_->Branch("dtClusterY",             dtClusterY,             "dtClusterY[nDtClusters]/F");
  tree_->Branch("dtClusterZ",             dtClusterZ,             "dtClusterZ[nDtClusters]/F");
  tree_->Branch("dtClusterTime",             dtClusterTime,             "dtClusterTime[nDtClusters]/F");
  tree_->Branch("dtClusterTimeSpread",             dtClusterTimeSpread,             "dtClusterTimeSpread[nDtClusters]/F");
  tree_->Branch("dtClusterGenMuonDeltaR",             dtClusterGenMuonDeltaR,             "dtClusterGenMuonDeltaR[nDtClusters]/F");

  tree_->Branch("dtClusterMajorAxis",             dtClusterMajorAxis,             "dtClusterMajorAxis[nDtClusters]/F");
  tree_->Branch("dtClusterMinorAxis",             dtClusterMinorAxis,             "dtClusterMinorAxis[nDtClusters]/F");
  tree_->Branch("dtClusterEtaPhiSpread",             dtClusterEtaPhiSpread,             "dtClusterEtaPhiSpread[nDtClusters]/F");
  tree_->Branch("dtClusterPhiSpread",             dtClusterPhiSpread,             "dtClusterPhiSpread[nDtClusters]/F");
  tree_->Branch("dtClusterEtaSpread",             dtClusterEtaSpread,             "dtClusterEtaSpread[nDtClusters]/F");
  tree_->Branch("dtClusterXSpread",             dtClusterXSpread,             "dtClusterXSpread[nDtClusters]/F");
  tree_->Branch("dtClusterYSpread",             dtClusterYSpread,             "dtClusterYSpread[nDtClusters]/F");
  tree_->Branch("dtClusterZSpread",             dtClusterZSpread,             "dtClusterZSpread[nDtClusters]/F");
  tree_->Branch("dtClusterPhi",             dtClusterPhi,             "dtClusterPhi[nDtClusters]/F");
  tree_->Branch("dtClusterEta",             dtClusterEta,             "dtClusterEta[nDtClusters]/F");
  tree_->Branch("dtClusterJetVetoPt",             dtClusterJetVetoPt,             "dtClusterJetVetoPt[nDtClusters]/F");
  tree_->Branch("dtClusterMuonVetoPt",             dtClusterMuonVetoPt,             "dtClusterMuonVetoPt[nDtClusters]/F");
  // tree_->Branch("dtClusterCaloJetVeto",             dtClusterCaloJetVeto,             "dtClusterCaloJetVeto[nDtClusters]/F");
  tree_->Branch("dtClusterJetVetoE",             dtClusterJetVetoE,             "dtClusterJetVetoE[nDtClusters]/F");
  tree_->Branch("dtClusterMuonVetoE",             dtClusterMuonVetoE,             "dtClusterMuonVetoE[nDtClusters]/F");
  tree_->Branch("dtClusterSize",             dtClusterSize,             "dtClusterSize[nDtClusters]/I");
  tree_->Branch("dtClusterMe11Ratio",             dtClusterMe11Ratio,             "dtClusterMe11Ratio[nDtClusters]/F");
  tree_->Branch("dtClusterMe12Ratio",             dtClusterMe12Ratio,             "dtClusterMe12Ratio[nDtClusters]/F");
  tree_->Branch("dtClusterNStation",             dtClusterNStation,             "dtClusterNStation[nDtClusters]/I");
  tree_->Branch("dtClusterMaxStation",             dtClusterMaxStation,             "dtClusterMaxStation[nDtClusters]/I");
  tree_->Branch("dtClusterMaxStationRatio",             dtClusterMaxStationRatio,             "dtClusterMaxStationRatio[nDtClusters]/F");
  tree_->Branch("dtClusterNChamber",             dtClusterNChamber,             "dtClusterNChamber[nDtClusters]/I");
  tree_->Branch("dtClusterMaxChamber",             dtClusterMaxChamber,             "dtClusterMaxChamber[nDtClusters]/I");
  tree_->Branch("dtClusterMaxChamberRatio",             dtClusterMaxChamberRatio,             "dtClusterMaxChamberRatio[nDtClusters]/F");
  tree_->Branch("dtClusterNSegmentChamberPlus11",             dtClusterNSegmentChamberPlus11,             "dtClusterNSegmentChamberPlus11[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberPlus12",             dtClusterNSegmentChamberPlus12,             "dtClusterNSegmentChamberPlus12[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberPlus13",             dtClusterNSegmentChamberPlus13,             "dtClusterNSegmentChamberPlus13[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberPlus21",             dtClusterNSegmentChamberPlus21,             "dtClusterNSegmentChamberPlus21[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberPlus22",             dtClusterNSegmentChamberPlus22,             "dtClusterNSegmentChamberPlus22[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberPlus31",             dtClusterNSegmentChamberPlus31,             "dtClusterNSegmentChamberPlus31[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberPlus32",             dtClusterNSegmentChamberPlus32,             "dtClusterNSegmentChamberPlus32[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberPlus41",             dtClusterNSegmentChamberPlus41,             "dtClusterNSegmentChamberPlus41[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberPlus42",             dtClusterNSegmentChamberPlus42,             "dtClusterNSegmentChamberPlus42[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberMinus11",             dtClusterNSegmentChamberMinus11,             "dtClusterNSegmentChamberMinus11[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberMinus12",             dtClusterNSegmentChamberMinus12,             "dtClusterNSegmentChamberMinus12[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberMinus13",             dtClusterNSegmentChamberMinus13,             "dtClusterNSegmentChamberMinus13[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberMinus21",             dtClusterNSegmentChamberMinus21,             "dtClusterNSegmentChamberMinus21[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberMinus22",             dtClusterNSegmentChamberMinus22,             "dtClusterNSegmentChamberMinus22[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberMinus31",             dtClusterNSegmentChamberMinus31,             "dtClusterNSegmentChamberMinus31[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberMinus32",             dtClusterNSegmentChamberMinus32,             "dtClusterNSegmentChamberMinus32[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberMinus41",             dtClusterNSegmentChamberMinus41,             "dtClusterNSegmentChamberMinus41[nDtClusters]/I");
  tree_->Branch("dtClusterNSegmentChamberMinus42",             dtClusterNSegmentChamberMinus42,             "dtClusterNSegmentChamberMinus42[nDtClusters]/I");

  // all dt SegClusters
  tree_->Branch("nDtSegClusters",             &nDtSegClusters, "nDtSegClusters/I");
  // tree_->Branch("nDt_JetVetoSegCluster0p4",             &nDt_JetVetoSegCluster0p4, "nDt_JetVetoSegCluster0p4/I");
  // tree_->Branch("nDt_JetMuonVetoSegCluster0p4",             &nDt_JetMuonVetoSegCluster0p4, "nDt_JetMuonVetoSegCluster0p4/I");
  // tree_->Branch("nDt_JetVetoSegCluster0p4_Me1112Veto",             &nDt_JetVetoSegCluster0p4_Me1112Veto, "nDt_JetVetoSegCluster0p4_Me1112Veto/I");
  // tree_->Branch("nDt_JetMuonVetoSegCluster0p4_Me1112Veto",             &nDt_JetMuonVetoSegCluster0p4_Me1112Veto, "nDt_JetMuonVetoSegCluster0p4_Me1112Veto/I");
  tree_->Branch("dtSegClusterX",             dtSegClusterX,             "dtSegClusterX[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterY",             dtSegClusterY,             "dtSegClusterY[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterZ",             dtSegClusterZ,             "dtSegClusterZ[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterTime",             dtSegClusterTime,             "dtSegClusterTime[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterTimeSpread",             dtSegClusterTimeSpread,             "dtSegClusterTimeSpread[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterGenMuonDeltaR",             dtSegClusterGenMuonDeltaR,             "dtSegClusterGenMuonDeltaR[nDtSegClusters]/F");

  tree_->Branch("dtSegClusterMajorAxis",             dtSegClusterMajorAxis,             "dtSegClusterMajorAxis[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterMinorAxis",             dtSegClusterMinorAxis,             "dtSegClusterMinorAxis[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterEtaPhiSpread",             dtSegClusterEtaPhiSpread,             "dtSegClusterEtaPhiSpread[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterPhiSpread",             dtSegClusterPhiSpread,             "dtSegClusterPhiSpread[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterEtaSpread",             dtSegClusterEtaSpread,             "dtSegClusterEtaSpread[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterXSpread",             dtSegClusterXSpread,             "dtSegClusterXSpread[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterYSpread",             dtSegClusterYSpread,             "dtSegClusterYSpread[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterZSpread",             dtSegClusterZSpread,             "dtSegClusterZSpread[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterPhi",             dtSegClusterPhi,             "dtSegClusterPhi[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterEta",             dtSegClusterEta,             "dtSegClusterEta[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterJetVetoPt",             dtSegClusterJetVetoPt,             "dtSegClusterJetVetoPt[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterMuonVetoPt",             dtSegClusterMuonVetoPt,             "dtSegClusterMuonVetoPt[nDtSegClusters]/F");
  // tree_->Branch("dtSegClusterCaloJetVeto",             dtSegClusterCaloJetVeto,             "dtSegClusterCaloJetVeto[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterJetVetoE",             dtSegClusterJetVetoE,             "dtSegClusterJetVetoE[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterMuonVetoE",             dtSegClusterMuonVetoE,             "dtSegClusterMuonVetoE[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterSize",             dtSegClusterSize,             "dtSegClusterSize[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterMe11Ratio",             dtSegClusterMe11Ratio,             "dtSegClusterMe11Ratio[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterMe12Ratio",             dtSegClusterMe12Ratio,             "dtSegClusterMe12Ratio[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterNStation",             dtSegClusterNStation,             "dtSegClusterNStation[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterMaxStation",             dtSegClusterMaxStation,             "dtSegClusterMaxStation[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterMaxStationRatio",             dtSegClusterMaxStationRatio,             "dtSegClusterMaxStationRatio[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterNChamber",             dtSegClusterNChamber,             "dtSegClusterNChamber[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterMaxChamber",             dtSegClusterMaxChamber,             "dtSegClusterMaxChamber[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterMaxChamberRatio",             dtSegClusterMaxChamberRatio,             "dtSegClusterMaxChamberRatio[nDtSegClusters]/F");
  tree_->Branch("dtSegClusterNSegmentChamberPlus11",             dtSegClusterNSegmentChamberPlus11,             "dtSegClusterNSegmentChamberPlus11[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberPlus12",             dtSegClusterNSegmentChamberPlus12,             "dtSegClusterNSegmentChamberPlus12[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberPlus13",             dtSegClusterNSegmentChamberPlus13,             "dtSegClusterNSegmentChamberPlus13[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberPlus21",             dtSegClusterNSegmentChamberPlus21,             "dtSegClusterNSegmentChamberPlus21[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberPlus22",             dtSegClusterNSegmentChamberPlus22,             "dtSegClusterNSegmentChamberPlus22[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberPlus31",             dtSegClusterNSegmentChamberPlus31,             "dtSegClusterNSegmentChamberPlus31[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberPlus32",             dtSegClusterNSegmentChamberPlus32,             "dtSegClusterNSegmentChamberPlus32[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberPlus41",             dtSegClusterNSegmentChamberPlus41,             "dtSegClusterNSegmentChamberPlus41[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberPlus42",             dtSegClusterNSegmentChamberPlus42,             "dtSegClusterNSegmentChamberPlus42[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberMinus11",             dtSegClusterNSegmentChamberMinus11,             "dtSegClusterNSegmentChamberMinus11[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberMinus12",             dtSegClusterNSegmentChamberMinus12,             "dtSegClusterNSegmentChamberMinus12[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberMinus13",             dtSegClusterNSegmentChamberMinus13,             "dtSegClusterNSegmentChamberMinus13[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberMinus21",             dtSegClusterNSegmentChamberMinus21,             "dtSegClusterNSegmentChamberMinus21[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberMinus22",             dtSegClusterNSegmentChamberMinus22,             "dtSegClusterNSegmentChamberMinus22[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberMinus31",             dtSegClusterNSegmentChamberMinus31,             "dtSegClusterNSegmentChamberMinus31[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberMinus32",             dtSegClusterNSegmentChamberMinus32,             "dtSegClusterNSegmentChamberMinus32[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberMinus41",             dtSegClusterNSegmentChamberMinus41,             "dtSegClusterNSegmentChamberMinus41[nDtSegClusters]/I");
  tree_->Branch("dtSegClusterNSegmentChamberMinus42",             dtSegClusterNSegmentChamberMinus42,             "dtSegClusterNSegmentChamberMinus42[nDtSegClusters]/I");
  */
    // all dt RechitClusters
    tree_->Branch("nDtRechitClusters",             &nDtRechitClusters, "nDtRechitClusters/I");
    tree_->Branch("dtRechitCluster_match_gParticle_id",             dtRechitCluster_match_gParticle_id,             "dtRechitCluster_match_gParticle_id[nCscSegClusters]/I");
    tree_->Branch("dtRechitCluster_match_gParticle_minDeltaR",             dtRechitCluster_match_gParticle_minDeltaR,             "dtRechitCluster_match_gParticle_minDeltaR[nCscSegClusters]/F");
    tree_->Branch("dtRechitCluster_match_gParticle_index",             dtRechitCluster_match_gParticle_index,             "dtRechitCluster_match_gParticle_index[nCscSegClusters]/I");
    // tree_->Branch("nDt_JetVetoRechitCluster0p4",             &nDt_JetVetoRechitCluster0p4, "nDt_JetVetoRechitCluster0p4/I");
    // tree_->Branch("nDt_JetMuonVetoRechitCluster0p4",             &nDt_JetMuonVetoRechitCluster0p4, "nDt_JetMuonVetoRechitCluster0p4/I");
    // tree_->Branch("nDt_JetVetoRechitCluster0p4_Me1112Veto",             &nDt_JetVetoRechitCluster0p4_Me1112Veto, "nDt_JetVetoRechitCluster0p4_Me1112Veto/I");
    // tree_->Branch("nDt_JetMuonVetoRechitCluster0p4_Me1112Veto",             &nDt_JetMuonVetoRechitCluster0p4_Me1112Veto, "nDt_JetMuonVetoRechitCluster0p4_Me1112Veto/I");
    tree_->Branch("dtRechitClusterX",             dtRechitClusterX,             "dtRechitClusterX[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterY",             dtRechitClusterY,             "dtRechitClusterY[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterZ",             dtRechitClusterZ,             "dtRechitClusterZ[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterTime",             dtRechitClusterTime,             "dtRechitClusterTime[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterTimeSpread",             dtRechitClusterTimeSpread,             "dtRechitClusterTimeSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterGenMuonDeltaR",             dtRechitClusterGenMuonDeltaR,             "dtRechitClusterGenMuonDeltaR[nDtRechitClusters]/F");

    tree_->Branch("dtRechitClusterMajorAxis",             dtRechitClusterMajorAxis,             "dtRechitClusterMajorAxis[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterMinorAxis",             dtRechitClusterMinorAxis,             "dtRechitClusterMinorAxis[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterEtaPhiSpread",             dtRechitClusterEtaPhiSpread,             "dtRechitClusterEtaPhiSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterPhiSpread",             dtRechitClusterPhiSpread,             "dtRechitClusterPhiSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterEtaSpread",             dtRechitClusterEtaSpread,             "dtRechitClusterEtaSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterXSpread",             dtRechitClusterXSpread,             "dtRechitClusterXSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterYSpread",             dtRechitClusterYSpread,             "dtRechitClusterYSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterZSpread",             dtRechitClusterZSpread,             "dtRechitClusterZSpread[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterPhi",             dtRechitClusterPhi,             "dtRechitClusterPhi[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterEta",             dtRechitClusterEta,             "dtRechitClusterEta[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterJetVetoPt",             dtRechitClusterJetVetoPt,             "dtRechitClusterJetVetoPt[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterMuonVetoPt",             dtRechitClusterMuonVetoPt,             "dtRechitClusterMuonVetoPt[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterCaloJetVeto",             dtRechitClusterCaloJetVeto,             "dtRechitClusterCaloJetVeto[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterJetVetoE",             dtRechitClusterJetVetoE,             "dtRechitClusterJetVetoE[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterMuonVetoE",             dtRechitClusterMuonVetoE,             "dtRechitClusterMuonVetoE[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterSize",             dtRechitClusterSize,             "dtRechitClusterSize[nDtRechitClusters]/I");
    //tree_->Branch("dtRechitClusterMe11Ratio",             dtRechitClusterMe11Ratio,             "dtRechitClusterMe11Ratio[nDtRechitClusters]/F");
    //tree_->Branch("dtRechitClusterMe12Ratio",             dtRechitClusterMe12Ratio,             "dtRechitClusterMe12Ratio[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterNStation",             dtRechitClusterNStation,             "dtRechitClusterNStation[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterMaxStation",             dtRechitClusterMaxStation,             "dtRechitClusterMaxStation[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterMaxStationRatio",             dtRechitClusterMaxStationRatio,             "dtRechitClusterMaxStationRatio[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterNChamber",             dtRechitClusterNChamber,             "dtRechitClusterNChamber[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterMaxChamber",             dtRechitClusterMaxChamber,             "dtRechitClusterMaxChamber[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterMaxChamberRatio",             dtRechitClusterMaxChamberRatio,             "dtRechitClusterMaxChamberRatio[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterNSegmentStation1",             dtRechitClusterNSegmentStation1,             "dtRechitClusterNSegmentStation1[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNSegmentStation2",             dtRechitClusterNSegmentStation2,             "dtRechitClusterNSegmentStation2[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNSegmentStation3",             dtRechitClusterNSegmentStation3,             "dtRechitClusterNSegmentStation3[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNSegmentStation4",             dtRechitClusterNSegmentStation4,             "dtRechitClusterNSegmentStation4[nDtRechitClusters]/I");
    /*    tree_->Branch("dtRechitClusterNRechitChamberPlus11",             dtRechitClusterNRechitChamberPlus11,             "dtRechitClusterNRechitChamberPlus11[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberPlus12",             dtRechitClusterNRechitChamberPlus12,             "dtRechitClusterNRechitChamberPlus12[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberPlus13",             dtRechitClusterNRechitChamberPlus13,             "dtRechitClusterNRechitChamberPlus13[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberPlus21",             dtRechitClusterNRechitChamberPlus21,             "dtRechitClusterNRechitChamberPlus21[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberPlus22",             dtRechitClusterNRechitChamberPlus22,             "dtRechitClusterNRechitChamberPlus22[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberPlus31",             dtRechitClusterNRechitChamberPlus31,             "dtRechitClusterNRechitChamberPlus31[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberPlus32",             dtRechitClusterNRechitChamberPlus32,             "dtRechitClusterNRechitChamberPlus32[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberPlus41",             dtRechitClusterNRechitChamberPlus41,             "dtRechitClusterNRechitChamberPlus41[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberPlus42",             dtRechitClusterNRechitChamberPlus42,             "dtRechitClusterNRechitChamberPlus42[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberMinus11",             dtRechitClusterNRechitChamberMinus11,             "dtRechitClusterNRechitChamberMinus11[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberMinus12",             dtRechitClusterNRechitChamberMinus12,             "dtRechitClusterNRechitChamberMinus12[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberMinus13",             dtRechitClusterNRechitChamberMinus13,             "dtRechitClusterNRechitChamberMinus13[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberMinus21",             dtRechitClusterNRechitChamberMinus21,             "dtRechitClusterNRechitChamberMinus21[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberMinus22",             dtRechitClusterNRechitChamberMinus22,             "dtRechitClusterNRechitChamberMinus22[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberMinus31",             dtRechitClusterNRechitChamberMinus31,             "dtRechitClusterNRechitChamberMinus31[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberMinus32",             dtRechitClusterNRechitChamberMinus32,             "dtRechitClusterNRechitChamberMinus32[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberMinus41",             dtRechitClusterNRechitChamberMinus41,             "dtRechitClusterNRechitChamberMinus41[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNRechitChamberMinus42",             dtRechitClusterNRechitChamberMinus42,             "dtRechitClusterNRechitChamberMinus42[nDtRechitClusters]/I");
    */
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
  
  tree_->Branch("gLLP_beta",          gLLP_beta,          "gLLP_beta[2]/F");
  tree_->Branch("gLLP_csc",           gLLP_csc,           "gLLP_csc[2]/F");
  tree_->Branch("gLLP_dt",           gLLP_dt,           "gLLP_dt[2]/F");
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
  // tree_->Branch("lepPassVetoId", lepPassVetoId, "lepPassVetoId[nLeptons]/O");

  // tree_->Branch("lepLoosePassId", lepLoosePassId, "lepLoosePassId[nLeptons]/O");
  // tree_->Branch("lepMediumPassId", lepMediumPassId, "lepMediumPassId[nLeptons]/O");
  // tree_->Branch("lepTightPassId", lepTightPassId, "lepTightPassId[nLeptons]/O");

  //muons
  tree_->Branch("nMuons",  &nMuons, "nMuons/I");
  tree_->Branch("muonE",      muonE,      "muonE[nMuons]/F");
  tree_->Branch("muonPt",     muonPt,     "muonPt[nMuons]/F");
  tree_->Branch("muonEta",    muonEta,    "muonEta[nMuons]/F");
  tree_->Branch("muonPhi",    muonPhi,    "muonPhi[nMuons]/F");

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
  // tree_->Branch("jetTightPassId", jetTightPassId, "jetTightPassId[nJets]/O");
  tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[601]/O"); //hardcoded
  // tree_->Branch("jetChargedEMEnergyFraction",   jetChargedEMEnergyFraction,   "jetChargedEMEnergyFraction[nJets]/F");
  // tree_->Branch("jetNeutralEMEnergyFraction",   jetNeutralEMEnergyFraction,   "jetNeutralEMEnergyFraction[nJets]/F");
  // tree_->Branch("jetChargedHadronEnergyFraction",   jetChargedHadronEnergyFraction,   "jetChargedHadronEnergyFraction[nJets]/F");
  // tree_->Branch("jetNeutralHadronEnergyFraction",   jetNeutralHadronEnergyFraction,   "jetNeutralHadronEnergyFraction[nJets]/F");
};
