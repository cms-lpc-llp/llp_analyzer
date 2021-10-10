#include "RazorHelper.h"
#include "TreeMuonSystemCombination.h"
#include "assert.h"
#include "TTree.h"
#include "DBSCAN.h"

// Constructor
TreeMuonSystemCombination::TreeMuonSystemCombination()
{
  InitVariables();
};
TreeMuonSystemCombination::~TreeMuonSystemCombination()
{
  if (f_) f_->Close();
};
void TreeMuonSystemCombination::InitVariables()
{
  runNum=0; lumiSec=0; evtNum=0; category=0; MC_condition = 0;
  npv=0; npu=0;
  pileupWeight = 0; pileupWeightUp = 0; pileupWeightDown = 0;
  higgsPtWeight = 0;
  lepOverallSF = 1.0;
  metSF = 1.0;

  ZCategory = 3;
  sf_facScaleUp = 0; sf_facScaleDown = 0; sf_renScaleUp = 0; sf_renScaleDown = 0; sf_facRenScaleUp = 0; sf_facRenScaleDown = 0;

  for (int i = 0; i < 9; i++)
  {
    higgsPtWeightSys[i] = 0;
    scaleWeights[i] = 0.0;
  }
  weight=-1.0;rho=-1;
  met=-1; metPhi=-1; metXYCorr=-1; metPhiXYCorr=-1; HT = 0.0; jetMet_dPhi = -999.;jetMet_dPhiMin = 999.;jetMet_dPhiMin4 = 999.;
  metNoMu = -1.;metPhiJESUp = -999.; metPhiJESDown = -999.;  metHEM=-999.; metPhiHEM = -999.;metHEMXYCorr=-999.; metPhiHEMXYCorr = -999.;
  EE_prefiring = true;
  Flag_HBHENoiseFilter = false; Flag_HBHEIsoNoiseFilter = false; Flag_BadPFMuonFilter = false; Flag_CSCTightHaloFilter = false; Flag_goodVertices = false;
  Flag_ecalBadCalibFilter = false; Flag_all = false; Flag_globalSuperTightHalo2016Filter = false; Flag_BadChargedCandidateFilter = false; Flag_eeBadScFilter = false;

  Flag2_HBHENoiseFilter = false; Flag2_HBHEIsoNoiseFilter = false; Flag2_BadPFMuonFilter = false; Flag2_globalSuperTightHalo2016Filter = false;
  Flag2_globalTightHalo2016Filter = false; Flag2_BadChargedCandidateFilter = false; Flag2_EcalDeadCellTriggerPrimitiveFilter = false; Flag2_ecalBadCalibFilter = false;
  Flag2_eeBadScFilter = false;
  Flag2_all = false;
  mH = 0; mX = 0; ctau = 0;



  metJESUp = -999.;metJESDown = -999.; metEENoise = -999.;metPhiEENoise = -999.;metEENoiseXYCorr = -999.;metPhiEENoiseXYCorr = -999.;
  gWPt = 0.0;
  gLepId = 0;
  gLepPt = 0.; gLepPhi = 0.; gLepEta = 0.; gLepE = 0.;
  gHiggsPt = 0.; gHiggsPhi = 0.; gHiggsEta = 0.; gHiggsE = 0.;
  //CSC

  nCscRechitClusters = 0;
  nDtRechitClusters = 0;
  nDtRechitClusters2 = 0;

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
  nDtStations25 = 0;
  nDtWheels25 = 0;
  nDTRechitsStation1 = 0;
  nDTRechitsStation2 = 0;
  nDTRechitsStation3 = 0;
  nDTRechitsStation4 = 0;

  nDTRechitsWheelMinus2 = 0;
  nDTRechitsWheelMinus1 = 0;
  nDTRechitsWheel0 = 0;
  nDTRechitsWheelPlus1 = 0;
  nDTRechitsWheelPlus2 = 0;
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
  nRpc = 0;
  nDTRechits = 0;
  nDtSeg = 0;
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

    rpcPhi[i] = -999.;   //[nCsc]
    rpcEta[i] = -999.;   //[nCsc]
    rpc_RE12[i] = false;   //[nCsc]
    rpc_RB1[i] = false;   //[nCsc]


     dtSegPhi[i] = -999.;   //[nCsc]
     dtSegEta[i] = -999.;   //[nCsc]
     dtSegStation[i] = -999;   //[nCsc]
     dtSegWheel[i] = -999;   //[nCsc]

  }

  for( int i = 0; i < N_MAX_DTRECHITS; i++ )
  {
    dtRechitsX[i] = -999.;   //[nCsc]
    dtRechitsY[i] = -999.;   //[nCsc]
    dtRechitsZ[i] = -999.;   //[nCsc]
    dtRechitsPhi[i] = -999.;   //[nCsc]
    dtRechitsEta[i] = -999.;   //[nCsc]
    dtRechitsStation[i] = -999;   //[nCsc]
    dtRechitsWheel[i] = -999;   //[nCsc]
    dtRechitsClusterId[i] = -999;   //[nCsc]

  }


  for( int i = 0; i < N_MAX_CSC; i++ )
  {




      cscRechitCluster_match_Me1112_0p4[i] = 0;
      cscRechitCluster_match_Me1112_0p6[i] = 0;
      cscRechitCluster_match_Me1112_0p8[i] = 0;
      cscRechitCluster_match_Me11_0p4[i] = 0;
      cscRechitCluster_match_Me11_0p6[i] = 0;
      cscRechitCluster_match_Me11_0p8[i] = 0;
      cscRechitCluster_match_Me12_0p4[i] = 0;
      cscRechitCluster_match_Me12_0p6[i] = 0;
      cscRechitCluster_match_Me12_0p8[i] = 0;

      cscRechitCluster_match_cscRechits_0p4[i] = 0;

      cscRechitCluster_match_cscSeg_0p4[i] = 0;
      cscRechitCluster_match_ME11Seg_0p4[i] = 0;
      cscRechitCluster_match_ME12Seg_0p4[i] = 0;
      cscRechitCluster_match_cscSeg_0p6[i] = 0;
      cscRechitCluster_match_ME11Seg_0p6[i] = 0;
      cscRechitCluster_match_ME12Seg_0p6[i] = 0;

      cscRechitCluster_match_dtRechits_0p4[i] = 0;
      cscRechitCluster_match_dtRechits_0p6[i] = 0;
      cscRechitCluster_match_dtRechits_phi0p2[i] = 0;
      cscRechitCluster_match_MB1_0p4[i] = 0;
      cscRechitCluster_match_MB1_0p6[i] = 0;
      cscRechitCluster_match_dtSeg_0p4[i] = 0;
      cscRechitCluster_match_dtSeg_0p6[i] = 0;
      cscRechitCluster_match_MB1Seg_0p4[i] = 0;
      cscRechitCluster_match_MB1Seg_0p6[i] = 0;
      cscRechitCluster_match_RB1_0p4[i] = 0;
      cscRechitCluster_match_RE12_0p4[i] = 0;
      cscRechitCluster_match_RB1_0p6[i] = 0;
      cscRechitCluster_match_RE12_0p6[i] = 0;

      cscRechitCluster_match_cluster_dR[i] = 999.;
      cscRechitCluster_match_cluster_index[i] = 999;


      cscRechitCluster_match_gParticle[i] = false;
      cscRechitCluster_match_gParticle_minDeltaR[i] = -999.;
      cscRechitCluster_match_gParticle_index[i] = -999;
      cscRechitCluster_match_gParticle_id[i] = -999;
      cscRechitCluster_match_gParticle_eta[i] = -999.;
      cscRechitCluster_match_gParticle_phi[i] = -999.;
      cscRechitCluster_match_gParticle_E[i] = -999.;
      cscRechitCluster_match_gParticle_pt[i] = -999.;
      cscRechitCluster_match_gParticle_MotherId[i]  = -999;

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
        cscRechitCluster_match_gLLP_dt[i] = false;
        cscRechitCluster_match_gLLP_multiplicity[i] = 999;
        cscRechitCluster_match_gLLP_EM_multiplicity[i] = 999;



        cscRechitCluster_match_gLLP_e[i] = 999.;
        cscRechitCluster_match_gLLP_pt[i] = 999.;
        cscRechitCluster_match_gLLP_EMFracE[i] = 999.;
        cscRechitCluster_match_gLLP_EMFracEz[i] = 999.;
        cscRechitCluster_match_gLLP_EMFracP[i] = 999.;
        cscRechitCluster_match_gLLP_EMFracPz[i] = 999.;
        cscRechitCluster_match_gLLP_lepdPhi[i] = 999.;
        cscRechitCluster_match_gLLP_daughter0_deltaR[i] = 999.0;
        cscRechitCluster_match_gLLP_daughter1_deltaR[i] = 999.0;
        cscRechitCluster_match_gLLP_daughter2_deltaR[i] = 999.0;
        cscRechitCluster_match_gLLP_daughter3_deltaR[i] = 999.0;
        cscRechitCluster_match_gLLP_other_daughter_deltaR[i] = 999.0;
        cscRechitCluster_match_gLLP_other_daughter_index[i] = 999;
        cscRechitCluster_match_gLLP_daughterKaon[i] = false;


        cscRechitCluster_match_gLLP_other_eta[i] = 999.0;
        cscRechitCluster_match_gLLP_other_phi[i] = 999.0;
        cscRechitCluster_match_gLLP_other_decay_r[i] = 999.0;
        cscRechitCluster_match_gLLP_other_decay_x[i] = 999.0;
        cscRechitCluster_match_gLLP_other_decay_y[i] = 999.0;
        cscRechitCluster_match_gLLP_other_decay_z[i] = 999.0;
        cscRechitCluster_match_gLLP_other_ctau[i] = 999.0;
        cscRechitCluster_match_gLLP_other_beta[i] = 999.0;
        cscRechitCluster_match_gLLP_other_csc[i] = false;
        cscRechitCluster_match_gLLP_other_e[i] = 999.0;
        cscRechitCluster_match_gLLP_other_pt[i] = 999.0;
        cscRechitClusterSize[i] = -999;
        cscRechitClusterX[i] = -999.;
        cscRechitClusterY[i] = -999.;
        cscRechitClusterZ[i] = -999.;
        cscRechitClusterTime[i] = -999.;
        cscRechitClusterTimeTotal[i] = -999.;
        cscRechitClusterTimeWire[i] = -999.;
        cscRechitClusterTimeWirePruned[i] = -999.;


        cscRechitClusterGenMuonDeltaR[i] = 999.;
        cscRechitClusterTimeSpread[i] = -999.;
        cscRechitClusterTimeTotalSpread[i] = -999.;
        cscRechitClusterTimeTotalSpreadPruned[i] = -999.;
        cscRechitClusterTimeWireSpread[i] = -999.;

        cscRechitClusterMajorAxis[i] = -999.;
        cscRechitClusterMinorAxis[i] = -999.;
        cscRechitClusterXSpread[i] = -999.;
        cscRechitClusterYSpread[i] = -999.;
        cscRechitClusterZSpread[i] = -999.;
        cscRechitClusterEtaPhiSpread[i] = -999.;
        cscRechitClusterXYSpread[i] = -999.;
        cscRechitClusterRSpread[i] = -999.;
        cscRechitClusterDeltaRSpread[i] =999.;

        cscRechitClusterEtaSpread[i] = -999.;
        cscRechitClusterPhiSpread[i] = -999.;
        cscRechitClusterEta[i] = -999.;
        cscRechitClusterPhi[i] = -999.;
        cscRechitClusterJetVetoPt[i] = 0.0;
        cscRechitClusterJetVetoEta[i] = 0.0;
        cscRechitClusterJetVetoPhi[i] = 0.0;

        cscRechitClusterJetVetoElectronEnergyFraction[i] = 0.0;
        cscRechitClusterJetVetoPhotonEnergyFraction[i] = 0.0;
        cscRechitClusterJetVetoChargedHadronEnergyFraction[i] = 0.0;
        cscRechitClusterJetVetoNeutralHadronEnergyFraction[i] = 0.0;
        cscRechitClusterJetVetoMuonEnergyFraction[i] = 0.0;
        cscRechitClusterJetVetoE[i] = 0.0;
        cscRechitClusterGenJetVetoPt[i] = 0.0;
        cscRechitClusterGenJetVetoE[i] = 0.0;
        cscRechitClusterMuonVetoPt[i] = 0.0;
        cscRechitClusterMuonVetoE[i] = 0.0;
        cscRechitClusterMuonVetoPhi[i] = 0.0;
        cscRechitClusterMuonVetoEta[i] = 0.0;
        cscRechitClusterMuonVetoLooseIso[i] = false;
        cscRechitClusterMuonVetoTightIso[i] = false;
        cscRechitClusterMuonVetoVTightIso[i] = false;
        cscRechitClusterMuonVetoVVTightIso[i] = false;
        cscRechitClusterMuonVetoTightId[i] = false;
        cscRechitClusterMuonVetoLooseId[i] = false;
        cscRechitClusterMuonVetoGlobal[i] = false;
        cscRechitClusterMuonVetoIso[i] = false;
        cscRechitClusterIsoMuonVetoPt[i] = 0.0;
        cscRechitClusterIsoMuonVetoE[i] = 0.0;
        cscRechitClusterIsoMuonVetoPhi[i] = 0.0;
        cscRechitClusterIsoMuonVetoEta[i] = 0.0;
        cscRechitClusterGenMuonVetoE[i] = 0.0;
        cscRechitClusterGenMuonVetoPt[i] = 0.0;
        cscRechitClusterMuonVetoType[i] = 999;

        // cscRechitClusterGenMuonVetoProdX[i] = 0.0;
        // cscRechitClusterGenMuonVetoProdY[i] = 0.0;
        // cscRechitClusterGenMuonVetoProdZ[i] = 0.0;
        // cscRechitClusterGenMuonVetoLLPDist[i] = 999.;
        // cscRechitClusterGenMuonVetoLLPIndex[i] = 999;

        cscRechitClusterJetVetoPt_0p6[i] = 0.0;
        cscRechitClusterJetVetoPt_0p8[i] = 0.0;
        cscRechitClusterJetVetoE_0p6[i] = 0.0;
        cscRechitClusterJetVetoE_0p8[i] = 0.0;
        cscRechitClusterMuonVetoPt_0p6[i] = 0.0;
        cscRechitClusterMuonVetoPt_0p8[i] = 0.0;
        cscRechitClusterMuonVetoE_0p6[i] = 0.0;
        cscRechitClusterMuonVetoE_0p8[i] = 0.0;

        cscRechitClusterZLep1[i] = false;
        cscRechitClusterZLep2[i] = false;
        cscRechitClusterZLep2Tag[i] = false;
        cscRechitClusterZLep1Tag[i] = false;
        cscRechitClusterZLep1Id[i] = -999;
        cscRechitClusterZLep2Id[i] = -999;
        cscRechitClusterZLep1LooseIso[i] = false;
        cscRechitClusterZLep1TightIso[i] = false;
        cscRechitClusterZLep1VTightIso[i] = false;
        cscRechitClusterZLep1VVTightIso[i] = false;
        cscRechitClusterZLep1TightId[i] = false;
        cscRechitClusterZLep2LooseIso[i] = false;
        cscRechitClusterZLep2TightIso[i] = false;
        cscRechitClusterZLep2VTightIso[i] = false;
        cscRechitClusterZLep2VVTightIso[i] = false;
        cscRechitClusterZLep2TightId[i] = false;
        cscRechitClusterNChamber[i] = -999;
        cscRechitClusterMaxChamberRatio[i] = -999.;
        cscRechitClusterMaxChamber[i] = -999;
        cscRechitClusterNStation[i] = -999;
        cscRechitClusterNStation5[i] = -999;
        cscRechitClusterNStation10[i] = -999;
        cscRechitClusterNStation10perc[i] = -999;
        cscRechitClusterAvgStation[i] = -999.;
        cscRechitClusterAvgStation5[i] = -999.;
        cscRechitClusterAvgStation10[i] = -999.;
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
        cscRechitClusterMetXYCorr_dPhi[i] = 999.;

        cscRechitClusterMetHEM_dPhi[i] = 999.;
        cscRechitClusterMetHEMXYCorr_dPhi[i] = 999.;
        cscRechitClusterMetEENoise_dPhi[i] = 999.;
        cscRechitClusterMetEENoiseXYCorr_dPhi[i] = 999.;
        cscRechitClusterMetJesUp_dPhi[i] = 999.;
        cscRechitClusterMetJesDown_dPhi[i] = 999.;


        cscRechitClusterNLayersChamberPlus11[i] = -999;
        cscRechitClusterNLayersChamberPlus12[i] = -999;
        cscRechitClusterNLayersChamberPlus13[i] = -999;
        cscRechitClusterNLayersChamberPlus21[i] = -999;
        cscRechitClusterNLayersChamberPlus22[i] = -999;
        cscRechitClusterNLayersChamberPlus31[i] = -999;
        cscRechitClusterNLayersChamberPlus32[i] = -999;
        cscRechitClusterNLayersChamberPlus41[i] = -999;
        cscRechitClusterNLayersChamberPlus42[i] = -999;
        cscRechitClusterNLayersChamberMinus11[i] = -999;
        cscRechitClusterNLayersChamberMinus12[i] = -999;
        cscRechitClusterNLayersChamberMinus13[i] = -999;
        cscRechitClusterNLayersChamberMinus21[i] = -999;
        cscRechitClusterNLayersChamberMinus22[i] = -999;
        cscRechitClusterNLayersChamberMinus31[i] = -999;
        cscRechitClusterNLayersChamberMinus32[i] = -999;
        cscRechitClusterNLayersChamberMinus41[i] = -999;
        cscRechitClusterNLayersChamberMinus42[i] = -999;

	dtRechitClusterZLep1[i] = false;
        dtRechitClusterZLep2[i] = false;
        dtRechitClusterZLep2Tag[i] = false;
        dtRechitClusterZLep1Tag[i] = false;
        dtRechitClusterZLep1Id[i] = -999;
        dtRechitClusterZLep2Id[i] = -999;
        dtRechitClusterZLep1LooseIso[i] = false;
        dtRechitClusterZLep1TightIso[i] = false;
        dtRechitClusterZLep1VTightIso[i] = false;
        dtRechitClusterZLep1VVTightIso[i] = false;
        dtRechitClusterZLep1TightId[i] = false;
        dtRechitClusterZLep2LooseIso[i] = false;
        dtRechitClusterZLep2TightIso[i] = false;
        dtRechitClusterZLep2VTightIso[i] = false;
        dtRechitClusterZLep2VVTightIso[i] = false;
        dtRechitClusterZLep2TightId[i] = false;


        dtRechitCluster_match_gParticle_Id[i] = -999;
        dtRechitCluster_match_gParticle_Pt[i] = -999.;
        dtRechitCluster_match_gParticle_Eta[i] = -999.;
        dtRechitCluster_match_gParticle_Phi[i] = -999.;
        dtRechitCluster_match_gParticle_E[i] = -999.;
        dtRechitCluster_match_gParticle_Status[i] = -999;
        dtRechitCluster_match_gParticle_MotherId[i] = -999;
        dtRechitCluster_match_gParticle_deltaR[i] = -999.;

        dtRechitCluster_match_MB1hits_0p4[i] = 0;
        dtRechitCluster_match_MB1hits_0p5[i] = 0;
        dtRechitCluster_match_MB1hits_cosmics_plus[i] = 0;
        dtRechitCluster_match_MB1hits_cosmics_minus[i] = 0;
        dtRechitCluster_match_MB1Seg_0p4[i] = 0;
        dtRechitCluster_match_MB1Seg_0p5[i] = 0;
        dtRechitCluster_match_RPChits_dPhi0p5[i] = 0;
        dtRechitCluster_match_RPCBx_dPhi0p5[i] = 0;
        dtRechitCluster_match_RB1_0p4[i] = 0;
        dtRechitCluster_match_RB1_dPhi0p5[i] = 0;
        dtRechitCluster_match_dtSeg_0p5[i] = 0;
        dtRechitCluster_match_dtSegTime_0p5[i] = 0;
        dtRechitCluster_match_dtSeg_0p4[i] = 0;
        dtRechitCluster_match_dtSegTime_0p4[i] = 0;

        dtRechitCluster_match_dtSegTimeSpread_0p5[i] = 0;
        dtRechitCluster_match_dtSegTimeSpread_0p4[i] = 0;
        dtRechitCluster_match_dtSeg_sameStation_0p5[i] = 0;
        dtRechitCluster_match_dtSegTime_sameStation_0p5[i] = 0;
        dtRechitCluster_match_dtSeg_sameStation_0p4[i] = 0;
        dtRechitCluster_match_dtSegTime_sameStation_0p4[i] = 0;

        dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5[i] = 0;
        dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4[i] = 0;



        dtRechitCluster_match_RPCTime_dPhi0p5[i] = 0;
        dtRechitCluster_match_RPCTimeSpread_dPhi0p5[i] = 0;
        dtRechitCluster_match_RPCTime_dR0p4[i] = 0;
        dtRechitCluster_match_RPCTimeSpread_dR0p4[i] = 0;
        dtRechitCluster_match_RPChits_dR0p4[i] = 0;
        dtRechitCluster_match_RPCTime_sameStation_dR0p4[i] = 0;
        dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[i] = 0;
        dtRechitCluster_match_RPChits_sameStation_dR0p4[i] = 0;



          dtRechitCluster_match_gLLP[i] = false;
          dtRechitCluster_match_gLLP_minDeltaR[i] = 999;
          dtRechitCluster_match_gLLP_index[i] = 999;
          dtRechitCluster_match_gLLP_eta[i] = 999.;
          dtRechitCluster_match_gLLP_phi[i] = 999.;
          dtRechitCluster_match_gLLP_decay_r[i] = 999.;
          dtRechitCluster_match_gLLP_decay_x[i] = 999.;
          dtRechitCluster_match_gLLP_decay_y[i] = 999.;
          dtRechitCluster_match_gLLP_decay_z[i] = 999.;
          dtRechitCluster_match_gLLP_ctau[i] = 999.;
          dtRechitCluster_match_gLLP_beta[i] = 999.;
          dtRechitCluster_match_gLLP_csc[i] = false;
          dtRechitCluster_match_gLLP_dt[i] = false;
          dtRechitCluster_match_gLLP_multiplicity[i] = 999;
          dtRechitCluster_match_gLLP_EM_multiplicity[i] = 999;



          dtRechitCluster_match_gLLP_e[i] = 999.;
          dtRechitCluster_match_gLLP_pt[i] = 999.;
          dtRechitCluster_match_gLLP_EMFracE[i] = 999.;
          dtRechitCluster_match_gLLP_EMFracEz[i] = 999.;
          dtRechitCluster_match_gLLP_EMFracP[i] = 999.;
          dtRechitCluster_match_gLLP_EMFracPz[i] = 999.;
          dtRechitCluster_match_gLLP_lepdPhi[i] = 999.;
          dtRechitCluster_match_gLLP_daughter0_deltaR[i] = 999.0;
          dtRechitCluster_match_gLLP_daughter1_deltaR[i] = 999.0;
          dtRechitCluster_match_gLLP_daughter2_deltaR[i] = 999.0;
          dtRechitCluster_match_gLLP_daughter3_deltaR[i] = 999.0;
          dtRechitCluster_match_gLLP_other_daughter_deltaR[i] = 999.0;
          dtRechitCluster_match_gLLP_other_daughter_index[i] = 999;
          dtRechitCluster_match_gLLP_daughterKaon[i] = false;


          dtRechitCluster_match_gLLP_other_eta[i] = 999.0;
          dtRechitCluster_match_gLLP_other_phi[i] = 999.0;
          dtRechitCluster_match_gLLP_other_decay_r[i] = 999.0;
          dtRechitCluster_match_gLLP_other_decay_x[i] = 999.0;
          dtRechitCluster_match_gLLP_other_decay_y[i] = 999.0;
          dtRechitCluster_match_gLLP_other_decay_z[i] = 999.0;
          dtRechitCluster_match_gLLP_other_ctau[i] = 999.0;
          dtRechitCluster_match_gLLP_other_beta[i] = 999.0;
          dtRechitCluster_match_gLLP_other_csc[i] = false;
          dtRechitCluster_match_gLLP_other_e[i] = 999.0;
          dtRechitCluster_match_gLLP_other_pt[i] = 999.0;
          dtRechitClusterSize[i] = -999;
          dtRechitClusterX[i] = -999.;
          dtRechitClusterY[i] = -999.;
          dtRechitClusterZ[i] = -999.;
          dtRechitClusterTime[i] = -999.;
          dtRechitClusterTimeTotal[i] = -999.;
          dtRechitClusterTimeWire[i] = -999.;
          dtRechitClusterTimeWirePruned[i] = -999.;
          dtRechitClusterWheel[i] = -999;


          dtRechitClusterGenMuonDeltaR[i] = 999.;
          dtRechitClusterTimeSpread[i] = -999.;
          dtRechitClusterTimeTotalSpread[i] = -999.;
          dtRechitClusterTimeTotalSpreadPruned[i] = -999.;
          dtRechitClusterTimeWireSpread[i] = -999.;

          dtRechitClusterMajorAxis[i] = -999.;
          dtRechitClusterMinorAxis[i] = -999.;
          dtRechitClusterXSpread[i] = -999.;
          dtRechitClusterYSpread[i] = -999.;
          dtRechitClusterZSpread[i] = -999.;
          dtRechitClusterEtaPhiSpread[i] = -999.;
          dtRechitClusterXYSpread[i] = -999.;
          dtRechitClusterRSpread[i] = -999.;
          dtRechitClusterDeltaRSpread[i] =999.;

          dtRechitClusterEtaSpread[i] = -999.;
          dtRechitClusterPhiSpread[i] = -999.;
          dtRechitClusterEta[i] = -999.;
          dtRechitClusterPhi[i] = -999.;
          dtRechitClusterJetVetoPt[i] = 0.0;
          dtRechitClusterJetVetoEta[i] = 0.0;
          dtRechitClusterJetVetoPhi[i] = 0.0;

          dtRechitClusterJetVetoElectronEnergyFraction[i] = 0.0;
          dtRechitClusterJetVetoPhotonEnergyFraction[i] = 0.0;
          dtRechitClusterJetVetoChargedHadronEnergyFraction[i] = 0.0;
          dtRechitClusterJetVetoNeutralHadronEnergyFraction[i] = 0.0;
          dtRechitClusterJetVetoMuonEnergyFraction[i] = 0.0;
          dtRechitClusterJetVetoE[i] = 0.0;
          dtRechitClusterGenJetVetoPt[i] = 0.0;
          dtRechitClusterGenJetVetoE[i] = 0.0;
          dtRechitClusterMuonVetoPt[i] = 0.0;
          dtRechitClusterMuonVetoE[i] = 0.0;
          dtRechitClusterMuonVetoPhi[i] = 0.0;
          dtRechitClusterMuonVetoEta[i] = 0.0;
          dtRechitClusterMuonVetoLooseIso[i] = false;
          dtRechitClusterMuonVetoTightIso[i] = false;
          dtRechitClusterMuonVetoVTightIso[i] = false;
          dtRechitClusterMuonVetoVVTightIso[i] = false;
          dtRechitClusterMuonVetoTightId[i] = false;
          dtRechitClusterMuonVetoLooseId[i] = false;
          dtRechitClusterMuonVetoGlobal[i] = false;
          dtRechitClusterMuonVetoIso[i] = false;
          dtRechitClusterIsoMuonVetoPt[i] = 0.0;
          dtRechitClusterIsoMuonVetoE[i] = 0.0;
          dtRechitClusterIsoMuonVetoPhi[i] = 0.0;
          dtRechitClusterIsoMuonVetoEta[i] = 0.0;
          dtRechitClusterGenMuonVetoE[i] = 0.0;
          dtRechitClusterGenMuonVetoPt[i] = 0.0;
          dtRechitClusterMuonVetoType[i] = 999;

          // dtRechitClusterGenMuonVetoProdX[i] = 0.0;
          // dtRechitClusterGenMuonVetoProdY[i] = 0.0;
          // dtRechitClusterGenMuonVetoProdZ[i] = 0.0;
          // dtRechitClusterGenMuonVetoLLPDist[i] = 999.;
          // dtRechitClusterGenMuonVetoLLPIndex[i] = 999;

          dtRechitClusterJetVetoPt_0p6[i] = 0.0;
          dtRechitClusterJetVetoPt_0p8[i] = 0.0;
          dtRechitClusterJetVetoE_0p6[i] = 0.0;
          dtRechitClusterJetVetoE_0p8[i] = 0.0;
          dtRechitClusterMuonVetoPt_0p6[i] = 0.0;
          dtRechitClusterMuonVetoPt_0p8[i] = 0.0;
          dtRechitClusterMuonVetoE_0p6[i] = 0.0;
          dtRechitClusterMuonVetoE_0p8[i] = 0.0;

          // dtRechitClusterZLep1[i] = false;
          // dtRechitClusterZLep2[i] = false;
          // dtRechitClusterZLep2Tag[i] = false;
          // dtRechitClusterZLep1Tag[i] = false;
          // dtRechitClusterZLep1Id[i] = -999;
          // dtRechitClusterZLep2Id[i] = -999;
          // dtRechitClusterZLep1LooseIso[i] = false;
          // dtRechitClusterZLep1TightIso[i] = false;
          // dtRechitClusterZLep1VTightIso[i] = false;
          // dtRechitClusterZLep1VVTightIso[i] = false;
          // dtRechitClusterZLep1TightId[i] = false;
          // dtRechitClusterZLep2LooseIso[i] = false;
          // dtRechitClusterZLep2TightIso[i] = false;
          // dtRechitClusterZLep2VTightIso[i] = false;
          // dtRechitClusterZLep2VVTightIso[i] = false;
          // dtRechitClusterZLep2TightId[i] = false;
          dtRechitClusterNChamber[i] = -999;
          dtRechitClusterMaxChamberRatio[i] = -999.;
          dtRechitClusterMaxChamber[i] = -999;
          dtRechitClusterNStation[i] = -999;
          dtRechitClusterNStation5[i] = -999;
          dtRechitClusterNStation10[i] = -999;
          dtRechitClusterNStation10perc[i] = -999;
          dtRechitClusterAvgStation[i] = -999.;
          dtRechitClusterAvgStation5[i] = -999.;
          dtRechitClusterAvgStation10[i] = -999.;
          dtRechitClusterAvgStation10perc[i] = -999.;
          dtRechitClusterMaxStationRatio[i] = -999.;
          dtRechitClusterMaxStation[i] = -999;
          dtRechitClusterMe11Ratio[i] = -999.;
          dtRechitClusterMe12Ratio[i] = -999.;

          dtRechitClusterNSegmentStation1[i] = -999;
          dtRechitClusterNSegmentStation2[i] = -999;
          dtRechitClusterNSegmentStation3[i] = -999;
          dtRechitClusterNSegmentStation4[i] = -999;


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
          dtRechitClusterMet_dPhi[i] = 999.;
          dtRechitClusterMetXYCorr_dPhi[i] = 999.;

          dtRechitClusterMetHEM_dPhi[i] = 999.;
          dtRechitClusterMetHEMXYCorr_dPhi[i] = 999.;
          dtRechitClusterMetEENoise_dPhi[i] = 999.;
          dtRechitClusterMetEENoiseXYCorr_dPhi[i] = 999.;
          dtRechitClusterMetJesUp_dPhi[i] = 999.;
          dtRechitClusterMetJesDown_dPhi[i] = 999.;


          dtRechitClusterNLayersChamberPlus11[i] = -999;
          dtRechitClusterNLayersChamberPlus12[i] = -999;
          dtRechitClusterNLayersChamberPlus13[i] = -999;
          dtRechitClusterNLayersChamberPlus21[i] = -999;
          dtRechitClusterNLayersChamberPlus22[i] = -999;
          dtRechitClusterNLayersChamberPlus31[i] = -999;
          dtRechitClusterNLayersChamberPlus32[i] = -999;
          dtRechitClusterNLayersChamberPlus41[i] = -999;
          dtRechitClusterNLayersChamberPlus42[i] = -999;
          dtRechitClusterNLayersChamberMinus11[i] = -999;
          dtRechitClusterNLayersChamberMinus12[i] = -999;
          dtRechitClusterNLayersChamberMinus13[i] = -999;
          dtRechitClusterNLayersChamberMinus21[i] = -999;
          dtRechitClusterNLayersChamberMinus22[i] = -999;
          dtRechitClusterNLayersChamberMinus31[i] = -999;
          dtRechitClusterNLayersChamberMinus32[i] = -999;
          dtRechitClusterNLayersChamberMinus41[i] = -999;
          dtRechitClusterNLayersChamberMinus42[i] = -999;


          //correct clusters

  	dtRechitCluster2ZLep1[i] = false;
          dtRechitCluster2ZLep2[i] = false;
          dtRechitCluster2ZLep2Tag[i] = false;
          dtRechitCluster2ZLep1Tag[i] = false;
          dtRechitCluster2ZLep1Id[i] = -999;
          dtRechitCluster2ZLep2Id[i] = -999;
          dtRechitCluster2ZLep1LooseIso[i] = false;
          dtRechitCluster2ZLep1TightIso[i] = false;
          dtRechitCluster2ZLep1VTightIso[i] = false;
          dtRechitCluster2ZLep1VVTightIso[i] = false;
          dtRechitCluster2ZLep1TightId[i] = false;
          dtRechitCluster2ZLep2LooseIso[i] = false;
          dtRechitCluster2ZLep2TightIso[i] = false;
          dtRechitCluster2ZLep2VTightIso[i] = false;
          dtRechitCluster2ZLep2VVTightIso[i] = false;
          dtRechitCluster2ZLep2TightId[i] = false;


          dtRechitCluster2_match_gParticle_Id[i] = -999;
          dtRechitCluster2_match_gParticle_Pt[i] = -999.;
          dtRechitCluster2_match_gParticle_Eta[i] = -999.;
          dtRechitCluster2_match_gParticle_Phi[i] = -999.;
          dtRechitCluster2_match_gParticle_E[i] = -999.;
          dtRechitCluster2_match_gParticle_Status[i] = -999;
          dtRechitCluster2_match_gParticle_MotherId[i] = -999;
          dtRechitCluster2_match_gParticle_deltaR[i] = -999.;

          dtRechitCluster2_match_MB1hits_0p4[i] = 0;
          dtRechitCluster2_match_MB1hits_0p5[i] = 0;
          dtRechitCluster2_match_MB1hits_cosmics_plus[i] = 0;
          dtRechitCluster2_match_MB1hits_cosmics_minus[i] = 0;
          dtRechitCluster2_match_MB1Seg_0p4[i] = 0;
          dtRechitCluster2_match_MB1Seg_0p5[i] = 0;
          dtRechitCluster2_match_RPChits_dPhi0p5[i] = 0;
          dtRechitCluster2_match_RPCBx_dPhi0p5[i] = 0;
          dtRechitCluster2_match_RB1_0p4[i] = 0;
          dtRechitCluster2_match_RB1_dPhi0p5[i] = 0;
          dtRechitCluster2_match_dtSeg_0p5[i] = 0;
          dtRechitCluster2_match_dtSegTime_0p5[i] = 0;
          dtRechitCluster2_match_dtSeg_0p4[i] = 0;
          dtRechitCluster2_match_dtSegTime_0p4[i] = 0;

          dtRechitCluster2_match_dtSegTimeSpread_0p5[i] = 0;
          dtRechitCluster2_match_dtSegTimeSpread_0p4[i] = 0;
          dtRechitCluster2_match_dtSeg_sameStation_0p5[i] = 0;
          dtRechitCluster2_match_dtSegTime_sameStation_0p5[i] = 0;
          dtRechitCluster2_match_dtSeg_sameStation_0p4[i] = 0;
          dtRechitCluster2_match_dtSegTime_sameStation_0p4[i] = 0;

          dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p5[i] = 0;
          dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p4[i] = 0;



          dtRechitCluster2_match_RPCTime_dPhi0p5[i] = 0;
          dtRechitCluster2_match_RPCTimeSpread_dPhi0p5[i] = 0;
          dtRechitCluster2_match_RPCTime_dR0p4[i] = 0;
          dtRechitCluster2_match_RPCTimeSpread_dR0p4[i] = 0;
          dtRechitCluster2_match_RPChits_dR0p4[i] = 0;
          dtRechitCluster2_match_RPCTime_sameStation_dR0p4[i] = 0;
          dtRechitCluster2_match_RPCTimeSpread_sameStation_dR0p4[i] = 0;
          dtRechitCluster2_match_RPChits_sameStation_dR0p4[i] = 0;



            dtRechitCluster2_match_gLLP[i] = false;
            dtRechitCluster2_match_gLLP_minDeltaR[i] = 999;
            dtRechitCluster2_match_gLLP_index[i] = 999;
            dtRechitCluster2_match_gLLP_eta[i] = 999.;
            dtRechitCluster2_match_gLLP_phi[i] = 999.;
            dtRechitCluster2_match_gLLP_decay_r[i] = 999.;
            dtRechitCluster2_match_gLLP_decay_x[i] = 999.;
            dtRechitCluster2_match_gLLP_decay_y[i] = 999.;
            dtRechitCluster2_match_gLLP_decay_z[i] = 999.;
            dtRechitCluster2_match_gLLP_ctau[i] = 999.;
            dtRechitCluster2_match_gLLP_beta[i] = 999.;
            dtRechitCluster2_match_gLLP_csc[i] = false;
            dtRechitCluster2_match_gLLP_dt[i] = false;
            dtRechitCluster2_match_gLLP_multiplicity[i] = 999;
            dtRechitCluster2_match_gLLP_EM_multiplicity[i] = 999;



            dtRechitCluster2_match_gLLP_e[i] = 999.;
            dtRechitCluster2_match_gLLP_pt[i] = 999.;
            dtRechitCluster2_match_gLLP_EMFracE[i] = 999.;
            dtRechitCluster2_match_gLLP_EMFracEz[i] = 999.;
            dtRechitCluster2_match_gLLP_EMFracP[i] = 999.;
            dtRechitCluster2_match_gLLP_EMFracPz[i] = 999.;
            dtRechitCluster2_match_gLLP_lepdPhi[i] = 999.;
            dtRechitCluster2_match_gLLP_daughter0_deltaR[i] = 999.0;
            dtRechitCluster2_match_gLLP_daughter1_deltaR[i] = 999.0;
            dtRechitCluster2_match_gLLP_daughter2_deltaR[i] = 999.0;
            dtRechitCluster2_match_gLLP_daughter3_deltaR[i] = 999.0;
            dtRechitCluster2_match_gLLP_other_daughter_deltaR[i] = 999.0;
            dtRechitCluster2_match_gLLP_other_daughter_index[i] = 999;
            dtRechitCluster2_match_gLLP_daughterKaon[i] = false;


            dtRechitCluster2_match_gLLP_other_eta[i] = 999.0;
            dtRechitCluster2_match_gLLP_other_phi[i] = 999.0;
            dtRechitCluster2_match_gLLP_other_decay_r[i] = 999.0;
            dtRechitCluster2_match_gLLP_other_decay_x[i] = 999.0;
            dtRechitCluster2_match_gLLP_other_decay_y[i] = 999.0;
            dtRechitCluster2_match_gLLP_other_decay_z[i] = 999.0;
            dtRechitCluster2_match_gLLP_other_ctau[i] = 999.0;
            dtRechitCluster2_match_gLLP_other_beta[i] = 999.0;
            dtRechitCluster2_match_gLLP_other_csc[i] = false;
            dtRechitCluster2_match_gLLP_other_e[i] = 999.0;
            dtRechitCluster2_match_gLLP_other_pt[i] = 999.0;
            dtRechitCluster2Size[i] = -999;
            dtRechitCluster2X[i] = -999.;
            dtRechitCluster2Y[i] = -999.;
            dtRechitCluster2Z[i] = -999.;
            dtRechitCluster2Time[i] = -999.;
            dtRechitCluster2TimeTotal[i] = -999.;
            dtRechitCluster2TimeWire[i] = -999.;
            dtRechitCluster2TimeWirePruned[i] = -999.;
            dtRechitCluster2Wheel[i] = -999;


            dtRechitCluster2GenMuonDeltaR[i] = 999.;
            dtRechitCluster2TimeSpread[i] = -999.;
            dtRechitCluster2TimeTotalSpread[i] = -999.;
            dtRechitCluster2TimeTotalSpreadPruned[i] = -999.;
            dtRechitCluster2TimeWireSpread[i] = -999.;

            dtRechitCluster2MajorAxis[i] = -999.;
            dtRechitCluster2MinorAxis[i] = -999.;
            dtRechitCluster2XSpread[i] = -999.;
            dtRechitCluster2YSpread[i] = -999.;
            dtRechitCluster2ZSpread[i] = -999.;
            dtRechitCluster2EtaPhiSpread[i] = -999.;
            dtRechitCluster2XYSpread[i] = -999.;
            dtRechitCluster2RSpread[i] = -999.;
            dtRechitCluster2DeltaRSpread[i] =999.;

            dtRechitCluster2EtaSpread[i] = -999.;
            dtRechitCluster2PhiSpread[i] = -999.;
            dtRechitCluster2Eta[i] = -999.;
            dtRechitCluster2Phi[i] = -999.;
            dtRechitCluster2JetVetoPt[i] = 0.0;
            dtRechitCluster2JetVetoEta[i] = 0.0;
            dtRechitCluster2JetVetoPhi[i] = 0.0;

            dtRechitCluster2JetVetoElectronEnergyFraction[i] = 0.0;
            dtRechitCluster2JetVetoPhotonEnergyFraction[i] = 0.0;
            dtRechitCluster2JetVetoChargedHadronEnergyFraction[i] = 0.0;
            dtRechitCluster2JetVetoNeutralHadronEnergyFraction[i] = 0.0;
            dtRechitCluster2JetVetoMuonEnergyFraction[i] = 0.0;
            dtRechitCluster2JetVetoE[i] = 0.0;
            dtRechitCluster2GenJetVetoPt[i] = 0.0;
            dtRechitCluster2GenJetVetoE[i] = 0.0;
            dtRechitCluster2MuonVetoPt[i] = 0.0;
            dtRechitCluster2MuonVetoE[i] = 0.0;
            dtRechitCluster2MuonVetoPhi[i] = 0.0;
            dtRechitCluster2MuonVetoEta[i] = 0.0;
            dtRechitCluster2MuonVetoLooseIso[i] = false;
            dtRechitCluster2MuonVetoTightIso[i] = false;
            dtRechitCluster2MuonVetoVTightIso[i] = false;
            dtRechitCluster2MuonVetoVVTightIso[i] = false;
            dtRechitCluster2MuonVetoTightId[i] = false;
            dtRechitCluster2MuonVetoLooseId[i] = false;
            dtRechitCluster2MuonVetoGlobal[i] = false;
            dtRechitCluster2MuonVetoIso[i] = false;
            dtRechitCluster2IsoMuonVetoPt[i] = 0.0;
            dtRechitCluster2IsoMuonVetoE[i] = 0.0;
            dtRechitCluster2IsoMuonVetoPhi[i] = 0.0;
            dtRechitCluster2IsoMuonVetoEta[i] = 0.0;
            dtRechitCluster2GenMuonVetoE[i] = 0.0;
            dtRechitCluster2GenMuonVetoPt[i] = 0.0;
            dtRechitCluster2MuonVetoType[i] = 999;

            // dtRechitCluster2GenMuonVetoProdX[i] = 0.0;
            // dtRechitCluster2GenMuonVetoProdY[i] = 0.0;
            // dtRechitCluster2GenMuonVetoProdZ[i] = 0.0;
            // dtRechitCluster2GenMuonVetoLLPDist[i] = 999.;
            // dtRechitCluster2GenMuonVetoLLPIndex[i] = 999;

            dtRechitCluster2JetVetoPt_0p6[i] = 0.0;
            dtRechitCluster2JetVetoPt_0p8[i] = 0.0;
            dtRechitCluster2JetVetoE_0p6[i] = 0.0;
            dtRechitCluster2JetVetoE_0p8[i] = 0.0;
            dtRechitCluster2MuonVetoPt_0p6[i] = 0.0;
            dtRechitCluster2MuonVetoPt_0p8[i] = 0.0;
            dtRechitCluster2MuonVetoE_0p6[i] = 0.0;
            dtRechitCluster2MuonVetoE_0p8[i] = 0.0;

            // dtRechitCluster2ZLep1[i] = false;
            // dtRechitCluster2ZLep2[i] = false;
            // dtRechitCluster2ZLep2Tag[i] = false;
            // dtRechitCluster2ZLep1Tag[i] = false;
            // dtRechitCluster2ZLep1Id[i] = -999;
            // dtRechitCluster2ZLep2Id[i] = -999;
            // dtRechitCluster2ZLep1LooseIso[i] = false;
            // dtRechitCluster2ZLep1TightIso[i] = false;
            // dtRechitCluster2ZLep1VTightIso[i] = false;
            // dtRechitCluster2ZLep1VVTightIso[i] = false;
            // dtRechitCluster2ZLep1TightId[i] = false;
            // dtRechitCluster2ZLep2LooseIso[i] = false;
            // dtRechitCluster2ZLep2TightIso[i] = false;
            // dtRechitCluster2ZLep2VTightIso[i] = false;
            // dtRechitCluster2ZLep2VVTightIso[i] = false;
            // dtRechitCluster2ZLep2TightId[i] = false;
            dtRechitCluster2NChamber[i] = -999;
            dtRechitCluster2MaxChamberRatio[i] = -999.;
            dtRechitCluster2MaxChamber[i] = -999;
            dtRechitCluster2NStation[i] = -999;
            dtRechitCluster2NStation5[i] = -999;
            dtRechitCluster2NStation10[i] = -999;
            dtRechitCluster2NStation10perc[i] = -999;
            dtRechitCluster2AvgStation[i] = -999.;
            dtRechitCluster2AvgStation5[i] = -999.;
            dtRechitCluster2AvgStation10[i] = -999.;
            dtRechitCluster2AvgStation10perc[i] = -999.;
            dtRechitCluster2MaxStationRatio[i] = -999.;
            dtRechitCluster2MaxStation[i] = -999;
            dtRechitCluster2Me11Ratio[i] = -999.;
            dtRechitCluster2Me12Ratio[i] = -999.;

            dtRechitCluster2NSegmentStation1[i] = -999;
            dtRechitCluster2NSegmentStation2[i] = -999;
            dtRechitCluster2NSegmentStation3[i] = -999;
            dtRechitCluster2NSegmentStation4[i] = -999;


            dtRechitCluster2NRechitChamberPlus11[i] = -999;
            dtRechitCluster2NRechitChamberPlus12[i] = -999;
            dtRechitCluster2NRechitChamberPlus13[i] = -999;
            dtRechitCluster2NRechitChamberPlus21[i] = -999;
            dtRechitCluster2NRechitChamberPlus22[i] = -999;
            dtRechitCluster2NRechitChamberPlus31[i] = -999;
            dtRechitCluster2NRechitChamberPlus32[i] = -999;
            dtRechitCluster2NRechitChamberPlus41[i] = -999;
            dtRechitCluster2NRechitChamberPlus42[i] = -999;
            dtRechitCluster2NRechitChamberMinus11[i] = -999;
            dtRechitCluster2NRechitChamberMinus12[i] = -999;
            dtRechitCluster2NRechitChamberMinus13[i] = -999;
            dtRechitCluster2NRechitChamberMinus21[i] = -999;
            dtRechitCluster2NRechitChamberMinus22[i] = -999;
            dtRechitCluster2NRechitChamberMinus31[i] = -999;
            dtRechitCluster2NRechitChamberMinus32[i] = -999;
            dtRechitCluster2NRechitChamberMinus41[i] = -999;
            dtRechitCluster2NRechitChamberMinus42[i] = -999;
            dtRechitCluster2Met_dPhi[i] = 999.;
            dtRechitCluster2MetXYCorr_dPhi[i] = 999.;

            dtRechitCluster2MetHEM_dPhi[i] = 999.;
            dtRechitCluster2MetHEMXYCorr_dPhi[i] = 999.;
            dtRechitCluster2MetEENoise_dPhi[i] = 999.;
            dtRechitCluster2MetEENoiseXYCorr_dPhi[i] = 999.;
            dtRechitCluster2MetJesUp_dPhi[i] = 999.;
            dtRechitCluster2MetJesDown_dPhi[i] = 999.;


            dtRechitCluster2NLayersChamberPlus11[i] = -999;
            dtRechitCluster2NLayersChamberPlus12[i] = -999;
            dtRechitCluster2NLayersChamberPlus13[i] = -999;
            dtRechitCluster2NLayersChamberPlus21[i] = -999;
            dtRechitCluster2NLayersChamberPlus22[i] = -999;
            dtRechitCluster2NLayersChamberPlus31[i] = -999;
            dtRechitCluster2NLayersChamberPlus32[i] = -999;
            dtRechitCluster2NLayersChamberPlus41[i] = -999;
            dtRechitCluster2NLayersChamberPlus42[i] = -999;
            dtRechitCluster2NLayersChamberMinus11[i] = -999;
            dtRechitCluster2NLayersChamberMinus12[i] = -999;
            dtRechitCluster2NLayersChamberMinus13[i] = -999;
            dtRechitCluster2NLayersChamberMinus21[i] = -999;
            dtRechitCluster2NLayersChamberMinus22[i] = -999;
            dtRechitCluster2NLayersChamberMinus31[i] = -999;
            dtRechitCluster2NLayersChamberMinus32[i] = -999;
            dtRechitCluster2NLayersChamberMinus41[i] = -999;
            dtRechitCluster2NLayersChamberMinus42[i] = -999;

      // for(int j = 0;j<N_MAX_CSC;j++)
      // {
      //   dtRechitCluster_match_dtSegT_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegX_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegY_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegZ_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegPhi_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegEta_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegWheel_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegStation_dR0p4[i][j] = -999;
      //
      //
      //   dtRechitCluster_match_dtSegT_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegX_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegY_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegZ_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegPhi_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegEta_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegWheel_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_dtSegStation_sameStation_dR0p4[i][j] = -999;
      //
      //
      //   dtRechitCluster_match_RPCBx_dPhi0p5[i][j] = -999;
      //   dtRechitCluster_match_RPCX_dPhi0p5[i][j] = -999;
      //   dtRechitCluster_match_RPCY_dPhi0p5[i][j] = -999;
      //   dtRechitCluster_match_RPCZ_dPhi0p5[i][j] = -999;
      //   dtRechitCluster_match_RPCPhi_dPhi0p5[i][j] = -999;
      //   dtRechitCluster_match_RPCEta_dPhi0p5[i][j] = -999;
      //   dtRechitCluster_match_RPCRing_dPhi0p5[i][j] = -999;
      //   dtRechitCluster_match_RPCLayer_dPhi0p5[i][j] = -999;
      //   dtRechitCluster_match_RPCSector_dPhi0p5[i][j] = -999;
      //
      //
      //
      //   dtRechitCluster_match_RPCBx_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCX_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCY_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCZ_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCPhi_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCEta_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCRing_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCLayer_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCSector_dR0p4[i][j] = -999;
      //
      //   dtRechitCluster_match_RPCBx_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCX_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCY_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCZ_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCPhi_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCEta_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCRing_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCLayer_sameStation_dR0p4[i][j] = -999;
      //   dtRechitCluster_match_RPCSector_sameStation_dR0p4[i][j] = -999;
      // }

  }

// dtRechitCluster_match_dtSegT_dR0p4.clear();
//
//   dtRechitCluster_match_dtSegX_dR0p4.clear();
//   dtRechitCluster_match_dtSegY_dR0p4.clear();
//   dtRechitCluster_match_dtSegZ_dR0p4.clear();
//   dtRechitCluster_match_dtSegPhi_dR0p4.clear();
//   dtRechitCluster_match_dtSegEta_dR0p4.clear();
//   dtRechitCluster_match_dtSegWheel_dR0p4.clear();
//   dtRechitCluster_match_dtSegStation_dR0p4.clear();
//
//
//   dtRechitCluster_match_dtSegT_sameStation_dR0p4.clear();
//   dtRechitCluster_match_dtSegX_sameStation_dR0p4.clear();
//   dtRechitCluster_match_dtSegY_sameStation_dR0p4.clear();
//   dtRechitCluster_match_dtSegZ_sameStation_dR0p4.clear();
//   dtRechitCluster_match_dtSegPhi_sameStation_dR0p4.clear();
//   dtRechitCluster_match_dtSegEta_sameStation_dR0p4.clear();
//   dtRechitCluster_match_dtSegWheel_sameStation_dR0p4.clear();
//   dtRechitCluster_match_dtSegStation_sameStation_dR0p4.clear();
//
//
//   dtRechitCluster_match_RPCBx_dPhi0p5.clear();
//   dtRechitCluster_match_RPCX_dPhi0p5.clear();
//   dtRechitCluster_match_RPCY_dPhi0p5.clear();
//   dtRechitCluster_match_RPCZ_dPhi0p5.clear();
//   dtRechitCluster_match_RPCPhi_dPhi0p5.clear();
//   dtRechitCluster_match_RPCEta_dPhi0p5.clear();
//   dtRechitCluster_match_RPCRing_dPhi0p5.clear();
//   dtRechitCluster_match_RPCLayer_dPhi0p5.clear();
//   dtRechitCluster_match_RPCSector_dPhi0p5.clear();
//
//
//
//   dtRechitCluster_match_RPCBx_dR0p4.clear();
//   dtRechitCluster_match_RPCX_dR0p4.clear();
//   dtRechitCluster_match_RPCY_dR0p4.clear();
//   dtRechitCluster_match_RPCZ_dR0p4.clear();
//   dtRechitCluster_match_RPCPhi_dR0p4.clear();
//   dtRechitCluster_match_RPCEta_dR0p4.clear();
//   dtRechitCluster_match_RPCRing_dR0p4.clear();
//   dtRechitCluster_match_RPCLayer_dR0p4.clear();
//   dtRechitCluster_match_RPCSector_dR0p4.clear();
//
//   dtRechitCluster_match_RPCBx_sameStation_dR0p4.clear();
//   dtRechitCluster_match_RPCX_sameStation_dR0p4.clear();
//   dtRechitCluster_match_RPCY_sameStation_dR0p4.clear();
//   dtRechitCluster_match_RPCZ_sameStation_dR0p4.clear();
//   dtRechitCluster_match_RPCPhi_sameStation_dR0p4.clear();
//   dtRechitCluster_match_RPCEta_sameStation_dR0p4.clear();
//   dtRechitCluster_match_RPCRing_sameStation_dR0p4.clear();
//   dtRechitCluster_match_RPCLayer_sameStation_dR0p4.clear();
//   dtRechitCluster_match_RPCSector_sameStation_dR0p4.clear();

  for(int i = 0;i<2;i++)
  {
    gLLP_multiplicity[i]= 0;
    gLLP_multiplicity20[i]= 0;
    gLLP_EM_multiplicity[i]= 0;
    gLLP_maxMatchedDis[i] = -999.;
    gLLP_eta[i] = 0.0;
    gLLP_phi[i] = 0.0;
    gLLP_beta[i] = 0.0;
    gLLP_e[i] = 0.0;
    gLLP_pt[i] = 0.0;
    gLLP_lepdPhi[i] = 0.0;
    gLLP_csc[i] = 0.0;
    gLLP_dt[i] = 0.0;
    gLLP_ctau[i] = 0.0;
    gLLP_decay_vertex_r[i] = 0.0;
    gLLP_decay_vertex_x[i] = 0.0;
    gLLP_decay_vertex_y[i] = 0.0;
    gLLP_decay_vertex_z[i] = 0.0;
    gLLP_daughterKaon[i] = false;
    gLLP_EMFracE[i] = 0.;
    gLLP_EMFracEz[i] = 0.;
    gLLP_EMFracP[i] = 0.;
    gLLP_EMFracPz[i] = 0.;
    gLLP_visE[i] = 0.;
    gLLP_visE20[i] = 0.;
    gLLP_visEz[i] = 0.;
    gLLP_visP[i] = 0.;
    gLLP_visPz[i] = 0.;
    gLLP_daughter_deltaR[i] = -999.0;
    gLLP_match_jet_minDeltaR[i] = -999.;
    gLLP_match_jet_index[i] = -999;
    gLLP_match_jet_pt[i] = -999.;

  }
  gLLP_deltaR  = -999.0;
  for(int i = 0;i<4;i++)
  {
    gLLP_daughter_pt[i] = -999.0;
    gLLP_daughter_eta[i] = -999.0;
    gLLP_daughter_phi[i] = -999.0;
    gLLP_daughter_e[i] = -999.0;
    gLLP_daughter_mass[i] = -999.0;
    gLLP_daughter_id[i] = 999;

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
    lepFromZ[i] = false;
    lepPassLooseIso[i] = false;
    lepPassTightIso[i] = false;
    lepPassVTightIso[i] = false;
    lepPassVTightIso[i] = false;
    lepEff[i] = 1.0;
    lepSF[i] = 1.0;
    lepTriggerSF[i] = 1.0;
    lepTightIdSF[i] = 1.0;
    lepLooseIdSF[i] = 1.0;
    lepTightIsoSF[i] = 1.0;
    lepLooseIsoSF[i] = 1.0;
    // lepTriggerMCEfficiency[i] = 1.0;
    // lepTightIdMCEfficiency[i] = 1.0;
    // lepLooseIdMCEfficiency[i] = 1.0;
    // lepTightIsoMCEfficiency[i] = 1.0;
    // lepLooseIsoMCEfficiency[i] = 1.0;
    lepTag[i] = false;

  }
  //Z-candidate
  ZMass1 = -999.; ZMass = -999.; ZPt = -999.; ZEta = -999.; ZPhi = -999.;
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
    jetPtJESUp[i] = -999.;
    jetPtJESDown[i] = -999.;
    jetEJESUp[i] = -999.;
    jetEJESDown[i] = -999.;
    JecUnc[i] = -999.;


    jetChargedEMEnergyFraction[i] = -999.;
    jetNeutralEMEnergyFraction[i] = -999.;
    jetChargedHadronEnergyFraction[i] = -999.;
    jetNeutralHadronEnergyFraction[i] = -999.;

    // jetMuonEnergyFraction[i] = -999.;
    jet_match_llp_minDeltaR[i] = -999.;
    jet_match_llp_index[i] = -999;
    jet_match_llp_pt[i] = -999.;



    jet_match_genJet_minDeltaR[i] = -999.;
    jet_match_genJet_index[i] = -999;
    jet_match_genJet_pt[i] = -999.;
    jetTightPassId[i] = false;
  }

  for(int i = 0; i <NTriggersMAX; i++){
    HLTDecision[i] = false;
  }
  METTrigger = false;
  METNoMuTrigger = false;
};

void TreeMuonSystemCombination::InitTree()
{
  assert(tree_);
  InitVariables();

  tree_->SetBranchAddress("runNum",      &runNum);
  tree_->SetBranchAddress("MC_condition",      &MC_condition);
  tree_->SetBranchAddress("lumiSec",     &lumiSec);
  tree_->SetBranchAddress("evtNum",      &evtNum);
  tree_->SetBranchAddress("category",    &category);
  tree_->SetBranchAddress("mX",      &mX);
  tree_->SetBranchAddress("mH",      &mH);
  tree_->SetBranchAddress("ctau",      &ctau);
  tree_->SetBranchAddress("ZCategory",      &ZCategory);

  tree_->SetBranchAddress("npv",         &npv);
  tree_->SetBranchAddress("npu",         &npu);
  tree_->SetBranchAddress("weight",      &weight);
  tree_->SetBranchAddress("higgsPtWeight",      &higgsPtWeight);
  tree_->SetBranchAddress("higgsPtWeightSys",      &higgsPtWeightSys);
  tree_->SetBranchAddress("scaleWeights",      &scaleWeights);


  tree_->SetBranchAddress("lepOverallSF",      &lepOverallSF);


  tree_->SetBranchAddress("sf_facScaleUp",      &sf_facScaleUp);
  tree_->SetBranchAddress("sf_facScaleDown",      &sf_facScaleDown);
  tree_->SetBranchAddress("sf_renScaleUp",      &sf_renScaleUp);
  tree_->SetBranchAddress("sf_renScaleDown",      &sf_renScaleDown);
  tree_->SetBranchAddress("sf_facRenScaleUp",      &sf_facRenScaleUp);
  tree_->SetBranchAddress("sf_facRenScaleDown",      &sf_facRenScaleDown);
  tree_->SetBranchAddress("metSF",      &metSF);

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
  tree_->SetBranchAddress("EE_prefiring",      &EE_prefiring);


  tree_->SetBranchAddress("rho",         &rho);
  tree_->SetBranchAddress("met",         &met);
  tree_->SetBranchAddress("HT",         &HT);
  tree_->SetBranchAddress("metNoMu",         &metNoMu);


  tree_->SetBranchAddress("metPhi",      &metPhi);
  tree_->SetBranchAddress("metXYCorr",      &metXYCorr);
  tree_->SetBranchAddress("metPhiXYCorr",      &metPhiXYCorr);

  tree_->SetBranchAddress("jetMet_dPhi",      &jetMet_dPhi);
  tree_->SetBranchAddress("jetMet_dPhiMin",      &jetMet_dPhiMin);
  tree_->SetBranchAddress("jetMet_dPhiMin4",      &jetMet_dPhiMin4);

  tree_->SetBranchAddress("metJESUp",      &metJESUp);
  tree_->SetBranchAddress("metJESDown",      &metJESDown);
  tree_->SetBranchAddress("metPhiJESUp",      &metPhiJESUp);
  tree_->SetBranchAddress("metPhiJESDown",      &metPhiJESDown);
  tree_->SetBranchAddress("metJESUpSF",      &metJESUpSF);
  tree_->SetBranchAddress("metJESDownSF",      &metJESDownSF);
  tree_->SetBranchAddress("metEENoise",      &metEENoise);
  tree_->SetBranchAddress("metPhiEENoise",      &metPhiEENoise);
  tree_->SetBranchAddress("metEENoiseXYCorr",      &metEENoiseXYCorr);
tree_->SetBranchAddress("metPhiEENoiseXYCorr",      &metPhiEENoiseXYCorr);

  tree_->SetBranchAddress("metHEM",      &metHEM);
  tree_->SetBranchAddress("metPhiHEM",      &metPhiHEM);

    tree_->SetBranchAddress("metHEMXYCorr",      &metHEMXYCorr);
    tree_->SetBranchAddress("metPhiHEMXYCorr",      &metPhiHEMXYCorr);


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

  tree_->SetBranchAddress("gWPt",      &gWPt);

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

  tree_->SetBranchAddress("nRpc",            &nRpc);
  tree_->SetBranchAddress("nDtSeg",            &nDtSeg);

  tree_->SetBranchAddress("nDTRechits",            &nDTRechits);
  tree_->SetBranchAddress("nDTPositiveYRechits",            &nDTPositiveYRechits);
  tree_->SetBranchAddress("nDTNegativeYRechits",            &nDTNegativeYRechits);
  tree_->SetBranchAddress("nDtRings",             &nDtRings);
  tree_->SetBranchAddress("nDtWheels25",             &nDtWheels25);
  tree_->SetBranchAddress("nDtStations25",             &nDtStations25);

  tree_->SetBranchAddress("nDTRechitsWheelMinus2",             &nDTRechitsWheelMinus2);
  tree_->SetBranchAddress("nDTRechitsWheelMinus1",             &nDTRechitsWheelMinus1);
  tree_->SetBranchAddress("nDTRechitsWheel0",             &nDTRechitsWheel0);
  tree_->SetBranchAddress("nDTRechitsWheelPlus1",             &nDTRechitsWheelPlus1);
  tree_->SetBranchAddress("nDTRechitsWheelPlus2",             &nDTRechitsWheelPlus2);

  tree_->SetBranchAddress("nDTRechitsStation1",             &nDTRechitsStation1);
  tree_->SetBranchAddress("nDTRechitsStation2",             &nDTRechitsStation2);
  tree_->SetBranchAddress("nDTRechitsStation3",             &nDTRechitsStation3);
  tree_->SetBranchAddress("nDTRechitsStation4",             &nDTRechitsStation4);



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

  tree_->SetBranchAddress("rpcPhi",             rpcPhi);
  tree_->SetBranchAddress("rpcEta",             rpcEta);
  tree_->SetBranchAddress("rpc_RE12",             rpc_RE12);
  tree_->SetBranchAddress("rpc_RB1",             rpc_RB1);
  tree_->SetBranchAddress("dtRechitsX",             dtRechitsX);
  tree_->SetBranchAddress("dtRechitsY",             dtRechitsY);
  tree_->SetBranchAddress("dtRechitsZ",             dtRechitsZ);
  tree_->SetBranchAddress("dtRechitsPhi",             dtRechitsPhi);
  tree_->SetBranchAddress("dtRechitsEta",             dtRechitsEta);
  tree_->SetBranchAddress("dtRechitsStation",             dtRechitsStation);
  tree_->SetBranchAddress("dtRechitsWheel",             dtRechitsWheel);
  tree_->SetBranchAddress("dtRechitsClusterId",             dtRechitsClusterId);
  tree_->SetBranchAddress("dtSegPhi",             dtSegPhi);
  tree_->SetBranchAddress("dtSegEta",             dtSegEta);
  tree_->SetBranchAddress("dtSegStation",             dtSegStation);
  tree_->SetBranchAddress("dtSegWheel",             dtSegWheel);


  tree_->SetBranchAddress("nDtRechitClusters",             &nDtRechitClusters);
  tree_->SetBranchAddress("nDtRechitClusters2",             &nDtRechitClusters2);


  tree_->SetBranchAddress("dtRechitCluster_match_gLLP",             &dtRechitCluster_match_gLLP);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_index",             &dtRechitCluster_match_gLLP_index);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_minDeltaR",             &dtRechitCluster_match_gLLP_minDeltaR);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_eta",             &dtRechitCluster_match_gLLP_eta);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_phi",             &dtRechitCluster_match_gLLP_phi);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_decay_r",             &dtRechitCluster_match_gLLP_decay_r);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_decay_x",             &dtRechitCluster_match_gLLP_decay_x);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_decay_y",             &dtRechitCluster_match_gLLP_decay_y);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_decay_z",             &dtRechitCluster_match_gLLP_decay_z);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_ctau",             &dtRechitCluster_match_gLLP_ctau);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_beta",             &dtRechitCluster_match_gLLP_beta);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_csc",             &dtRechitCluster_match_gLLP_csc);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_dt",             &dtRechitCluster_match_gLLP_dt);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_multiplicity",             &dtRechitCluster_match_gLLP_multiplicity);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_EM_multiplicity",             &dtRechitCluster_match_gLLP_EM_multiplicity);

  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_e",             &dtRechitCluster_match_gLLP_e);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_pt",             &dtRechitCluster_match_gLLP_pt);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_EMFracE",             &dtRechitCluster_match_gLLP_EMFracE);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_EMFracEz",             &dtRechitCluster_match_gLLP_EMFracEz);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_EMFracP",             &dtRechitCluster_match_gLLP_EMFracP);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_EMFracPz",             &dtRechitCluster_match_gLLP_EMFracPz);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_lepdPhi",             &dtRechitCluster_match_gLLP_lepdPhi);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_daughterKaon",             &dtRechitCluster_match_gLLP_daughterKaon);

  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_daughter0_deltaR",             &dtRechitCluster_match_gLLP_daughter0_deltaR);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_daughter1_deltaR",             &dtRechitCluster_match_gLLP_daughter1_deltaR);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_daughter2_deltaR",             &dtRechitCluster_match_gLLP_daughter2_deltaR);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_daughter3_deltaR",             &dtRechitCluster_match_gLLP_daughter3_deltaR);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_daughter_deltaR",             &dtRechitCluster_match_gLLP_other_daughter_deltaR);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_daughter_index",             &dtRechitCluster_match_gLLP_other_daughter_index);

  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_daughter_deltaR",             &dtRechitCluster_match_gLLP_other_daughter_deltaR);

 tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_eta",             &dtRechitCluster_match_gLLP_other_eta);
 tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_phi",             &dtRechitCluster_match_gLLP_other_phi);
 tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_decay_r",             &dtRechitCluster_match_gLLP_other_decay_r);
 tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_decay_x",             &dtRechitCluster_match_gLLP_other_decay_x);
 tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_decay_y",             &dtRechitCluster_match_gLLP_other_decay_y);
 tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_decay_z",             &dtRechitCluster_match_gLLP_other_decay_z);
 tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_ctau",             &dtRechitCluster_match_gLLP_other_ctau);
 tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_beta",             &dtRechitCluster_match_gLLP_other_beta);
 tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_csc",             &dtRechitCluster_match_gLLP_other_csc);
 tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_e",             &dtRechitCluster_match_gLLP_other_e);
 tree_->SetBranchAddress("dtRechitCluster_match_gLLP_other_pt",             &dtRechitCluster_match_gLLP_other_pt);

  tree_->SetBranchAddress("dtRechitClusterMe11Ratio",             &dtRechitClusterMe11Ratio);
  tree_->SetBranchAddress("dtRechitClusterMe12Ratio",             &dtRechitClusterMe12Ratio);
  tree_->SetBranchAddress("dtRechitClusterMaxStation",             &dtRechitClusterMaxStation);
  tree_->SetBranchAddress("dtRechitClusterMaxStationRatio",             &dtRechitClusterMaxStationRatio);
  tree_->SetBranchAddress("dtRechitClusterNStation",             &dtRechitClusterNStation);
  tree_->SetBranchAddress("dtRechitClusterNStation5",             &dtRechitClusterNStation5);
  tree_->SetBranchAddress("dtRechitClusterNStation10",             &dtRechitClusterNStation10);
  tree_->SetBranchAddress("dtRechitClusterNStation10perc",             &dtRechitClusterNStation10perc);
  tree_->SetBranchAddress("dtRechitClusterAvgStation",             &dtRechitClusterAvgStation);
  tree_->SetBranchAddress("dtRechitClusterAvgStation5",             &dtRechitClusterAvgStation5);
  tree_->SetBranchAddress("dtRechitClusterAvgStation10",             &dtRechitClusterAvgStation10);
  tree_->SetBranchAddress("dtRechitClusterAvgStation10perc",             &dtRechitClusterAvgStation10perc);
  tree_->SetBranchAddress("dtRechitClusterMaxChamber",             &dtRechitClusterMaxChamber);
  tree_->SetBranchAddress("dtRechitClusterMaxChamberRatio",             &dtRechitClusterMaxChamberRatio);
  tree_->SetBranchAddress("dtRechitClusterNChamber",             &dtRechitClusterNChamber);
  tree_->SetBranchAddress("dtRechitClusterX",             dtRechitClusterX);
  tree_->SetBranchAddress("dtRechitClusterY",             dtRechitClusterY);
  tree_->SetBranchAddress("dtRechitClusterZ",             dtRechitClusterZ);
  tree_->SetBranchAddress("dtRechitClusterTime",             dtRechitClusterTime);
  tree_->SetBranchAddress("dtRechitClusterTimeWire",             dtRechitClusterTimeWire);
  tree_->SetBranchAddress("dtRechitClusterTimeWirePruned",             dtRechitClusterTimeWirePruned);
  tree_->SetBranchAddress("dtRechitClusterTimeTotal",             dtRechitClusterTimeTotal);
  tree_->SetBranchAddress("dtRechitClusterWheel",             dtRechitClusterWheel);

  tree_->SetBranchAddress("dtRechitClusterGenMuonDeltaR",             dtRechitClusterGenMuonDeltaR);

  tree_->SetBranchAddress("dtRechitClusterTimeSpread",             dtRechitClusterTimeSpread);
  tree_->SetBranchAddress("dtRechitClusterTimeTotalSpread",             dtRechitClusterTimeTotalSpread);
  tree_->SetBranchAddress("dtRechitClusterTimeTotalSpreadPruned",             dtRechitClusterTimeTotalSpreadPruned);
  tree_->SetBranchAddress("dtRechitClusterTimeWireSpread",             dtRechitClusterTimeWireSpread);
  tree_->SetBranchAddress("dtRechitClusterMajorAxis",             dtRechitClusterMajorAxis);
  tree_->SetBranchAddress("dtRechitClusterMinorAxis",             dtRechitClusterMinorAxis);
  tree_->SetBranchAddress("dtRechitClusterRSpread",             dtRechitClusterRSpread);

  tree_->SetBranchAddress("dtRechitClusterXSpread",             dtRechitClusterXSpread);
  tree_->SetBranchAddress("dtRechitClusterYSpread",             dtRechitClusterYSpread);
  tree_->SetBranchAddress("dtRechitClusterZSpread",             dtRechitClusterZSpread);
  tree_->SetBranchAddress("dtRechitClusterEtaPhiSpread",             dtRechitClusterEtaPhiSpread);
  tree_->SetBranchAddress("dtRechitClusterXYSpread",             dtRechitClusterXYSpread);
  tree_->SetBranchAddress("dtRechitClusterEtaSpread",             dtRechitClusterEtaSpread);
  tree_->SetBranchAddress("dtRechitClusterPhiSpread",             dtRechitClusterPhiSpread);
  tree_->SetBranchAddress("dtRechitClusterDeltaRSpread",             dtRechitClusterDeltaRSpread);
  tree_->SetBranchAddress("dtRechitClusterEta",             dtRechitClusterEta);
  tree_->SetBranchAddress("dtRechitClusterPhi",             dtRechitClusterPhi);


  tree_->SetBranchAddress("dtRechitClusterGenJetVetoPt",             dtRechitClusterGenJetVetoPt);
  tree_->SetBranchAddress("dtRechitClusterGenJetVetoE",             dtRechitClusterGenJetVetoE);
  tree_->SetBranchAddress("dtRechitClusterJetVetoPt",             dtRechitClusterJetVetoPt);
  tree_->SetBranchAddress("dtRechitClusterJetVetoEta",             dtRechitClusterJetVetoEta);
  tree_->SetBranchAddress("dtRechitClusterJetVetoPhi",             dtRechitClusterJetVetoPhi);

  tree_->SetBranchAddress("dtRechitClusterJetVetoElectronEnergyFraction",             dtRechitClusterJetVetoElectronEnergyFraction);
  tree_->SetBranchAddress("dtRechitClusterJetVetoPhotonEnergyFraction",             dtRechitClusterJetVetoPhotonEnergyFraction);
  tree_->SetBranchAddress("dtRechitClusterJetVetoNeutralHadronEnergyFraction",             dtRechitClusterJetVetoNeutralHadronEnergyFraction);
  tree_->SetBranchAddress("dtRechitClusterJetVetoChargedHadronEnergyFraction",             dtRechitClusterJetVetoChargedHadronEnergyFraction);
  tree_->SetBranchAddress("dtRechitClusterJetVetoMuonEnergyFraction",             dtRechitClusterJetVetoMuonEnergyFraction);


  tree_->SetBranchAddress("dtRechitClusterJetVetoE",             dtRechitClusterJetVetoE);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoPt",             dtRechitClusterMuonVetoPt);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoE",             dtRechitClusterMuonVetoE);

  tree_->SetBranchAddress("dtRechitClusterJetVetoPt_0p6",             dtRechitClusterJetVetoPt_0p6);
  tree_->SetBranchAddress("dtRechitClusterJetVetoPt_0p8",             dtRechitClusterJetVetoPt_0p8);
  tree_->SetBranchAddress("dtRechitClusterJetVetoE_0p6",             dtRechitClusterJetVetoE_0p6);
  tree_->SetBranchAddress("dtRechitClusterJetVetoE_0p8",             dtRechitClusterJetVetoE_0p8);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoPt_0p6",             dtRechitClusterMuonVetoPt_0p6);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoPt_0p8",             dtRechitClusterMuonVetoPt_0p8);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoE_0p6",             dtRechitClusterMuonVetoE_0p6);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoE_0p8",             dtRechitClusterMuonVetoE_0p8);

  // tree_->SetBranchAddress("dtRechitClusterZLep1",             dtRechitClusterZLep1);
  // tree_->SetBranchAddress("dtRechitClusterZLep2",             dtRechitClusterZLep2);
  // tree_->SetBranchAddress("dtRechitClusterZLep1Tag",             dtRechitClusterZLep1Tag);
  // tree_->SetBranchAddress("dtRechitClusterZLep2Tag",             dtRechitClusterZLep2Tag);
  // tree_->SetBranchAddress("dtRechitClusterZLep1Id",             dtRechitClusterZLep1Id);
  // tree_->SetBranchAddress("dtRechitClusterZLep2Id",             dtRechitClusterZLep2Id);
  // tree_->SetBranchAddress("dtRechitClusterZLep1LooseIso",             dtRechitClusterZLep1LooseIso);
  // tree_->SetBranchAddress("dtRechitClusterZLep1TightIso",             dtRechitClusterZLep1TightIso);
  // tree_->SetBranchAddress("dtRechitClusterZLep1VTightIso",             dtRechitClusterZLep1VTightIso);
  // tree_->SetBranchAddress("dtRechitClusterZLep1VVTightIso",             dtRechitClusterZLep1VVTightIso);
  // tree_->SetBranchAddress("dtRechitClusterZLep1TightId",             dtRechitClusterZLep1TightId);
  // tree_->SetBranchAddress("dtRechitClusterZLep2LooseIso",             dtRechitClusterZLep2LooseIso);
  // tree_->SetBranchAddress("dtRechitClusterZLep2TightIso",             dtRechitClusterZLep2TightIso);
  // tree_->SetBranchAddress("dtRechitClusterZLep2VTightIso",             dtRechitClusterZLep2VTightIso);
  // tree_->SetBranchAddress("dtRechitClusterZLep2VVTightIso",             dtRechitClusterZLep2VVTightIso);
  // tree_->SetBranchAddress("dtRechitClusterZLep2TightId",             dtRechitClusterZLep2TightId);

  tree_->SetBranchAddress("dtRechitClusterMuonVetoPhi",             dtRechitClusterMuonVetoPhi);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoEta",             dtRechitClusterMuonVetoEta);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoIso",             dtRechitClusterMuonVetoIso);

    tree_->SetBranchAddress("dtRechitClusterMuonVetoLooseIso",             dtRechitClusterMuonVetoLooseIso);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoTightIso",             dtRechitClusterMuonVetoTightIso);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoVTightIso",             dtRechitClusterMuonVetoVTightIso);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoVVTightIso",             dtRechitClusterMuonVetoVVTightIso);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoTightId",             dtRechitClusterMuonVetoTightId);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoLooseId",             dtRechitClusterMuonVetoLooseId);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoGlobal",             dtRechitClusterMuonVetoGlobal);


  tree_->SetBranchAddress("dtRechitClusterIsoMuonVetoPt",             dtRechitClusterIsoMuonVetoPt);
  tree_->SetBranchAddress("dtRechitClusterIsoMuonVetoE",             dtRechitClusterIsoMuonVetoE);
  tree_->SetBranchAddress("dtRechitClusterIsoMuonVetoPhi",             dtRechitClusterIsoMuonVetoPhi);
  tree_->SetBranchAddress("dtRechitClusterIsoMuonVetoEta",             dtRechitClusterIsoMuonVetoEta);
  tree_->SetBranchAddress("dtRechitClusterGenMuonVetoPt",             dtRechitClusterGenMuonVetoPt);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoType",             dtRechitClusterMuonVetoType);
  // tree_->SetBranchAddress("dtRechitClusterGenMuonVetoProdX",             dtRechitClusterGenMuonVetoProdX);
  // tree_->SetBranchAddress("dtRechitClusterGenMuonVetoProdY",             dtRechitClusterGenMuonVetoProdY);
  // tree_->SetBranchAddress("dtRechitClusterGenMuonVetoProdZ",             dtRechitClusterGenMuonVetoProdZ);
  // tree_->SetBranchAddress("dtRechitClusterGenMuonVetoLLPDist",             dtRechitClusterGenMuonVetoLLPDist);
  // tree_->SetBranchAddress("dtRechitClusterGenMuonVetoLLPIndex",             dtRechitClusterGenMuonVetoLLPIndex);

  tree_->SetBranchAddress("dtRechitClusterGenMuonVetoE",             dtRechitClusterGenMuonVetoE);
  tree_->SetBranchAddress("dtRechitClusterSize",             dtRechitClusterSize);

  tree_->SetBranchAddress("dtRechitClusterNSegmentStation1",             dtRechitClusterNSegmentStation1);
  tree_->SetBranchAddress("dtRechitClusterNSegmentStation2",             dtRechitClusterNSegmentStation2);
  tree_->SetBranchAddress("dtRechitClusterNSegmentStation3",             dtRechitClusterNSegmentStation3);
  tree_->SetBranchAddress("dtRechitClusterNSegmentStation4",             dtRechitClusterNSegmentStation4);



  tree_->SetBranchAddress("dtRechitClusterNRechitChamberPlus11",             dtRechitClusterNRechitChamberPlus11);
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

  tree_->SetBranchAddress("dtRechitClusterNLayersChamberPlus11",             dtRechitClusterNLayersChamberPlus11);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberPlus12",             dtRechitClusterNLayersChamberPlus12);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberPlus13",             dtRechitClusterNLayersChamberPlus13);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberPlus21",             dtRechitClusterNLayersChamberPlus21);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberPlus22",             dtRechitClusterNLayersChamberPlus22);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberPlus31",             dtRechitClusterNLayersChamberPlus31);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberPlus32",             dtRechitClusterNLayersChamberPlus32);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberPlus41",             dtRechitClusterNLayersChamberPlus41);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberPlus42",             dtRechitClusterNLayersChamberPlus42);

  tree_->SetBranchAddress("dtRechitClusterNLayersChamberMinus11",             dtRechitClusterNLayersChamberMinus11);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberMinus12",             dtRechitClusterNLayersChamberMinus12);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberMinus13",             dtRechitClusterNLayersChamberMinus13);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberMinus21",             dtRechitClusterNLayersChamberMinus21);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberMinus22",             dtRechitClusterNLayersChamberMinus22);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberMinus31",             dtRechitClusterNLayersChamberMinus31);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberMinus32",             dtRechitClusterNLayersChamberMinus32);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberMinus41",             dtRechitClusterNLayersChamberMinus41);
  tree_->SetBranchAddress("dtRechitClusterNLayersChamberMinus42",             dtRechitClusterNLayersChamberMinus42);




  tree_->SetBranchAddress("dtRechitCluster_match_gParticle_Id",             dtRechitCluster_match_gParticle_Id);
  tree_->SetBranchAddress("dtRechitCluster_match_gParticle_Pt",             dtRechitCluster_match_gParticle_Pt);
  tree_->SetBranchAddress("dtRechitCluster_match_gParticle_Eta",             dtRechitCluster_match_gParticle_Eta);
  tree_->SetBranchAddress("dtRechitCluster_match_gParticle_Phi",             dtRechitCluster_match_gParticle_Phi);
  tree_->SetBranchAddress("dtRechitCluster_match_gParticle_E",             dtRechitCluster_match_gParticle_E);
  tree_->SetBranchAddress("dtRechitCluster_match_gParticle_Status",             dtRechitCluster_match_gParticle_Status);
  tree_->SetBranchAddress("dtRechitCluster_match_gParticle_MotherId",             dtRechitCluster_match_gParticle_MotherId);
  tree_->SetBranchAddress("dtRechitCluster_match_gParticle_deltaR",             dtRechitCluster_match_gParticle_deltaR);


  tree_->SetBranchAddress("dtRechitCluster_match_MB1hits_0p4",             dtRechitCluster_match_MB1hits_0p4);
  tree_->SetBranchAddress("dtRechitCluster_match_MB1hits_0p5",             dtRechitCluster_match_MB1hits_0p5);
  tree_->SetBranchAddress("dtRechitCluster_match_MB1hits_cosmics_plus",             dtRechitCluster_match_MB1hits_cosmics_plus);
  tree_->SetBranchAddress("dtRechitCluster_match_MB1hits_cosmics_minus",             dtRechitCluster_match_MB1hits_cosmics_minus);
  tree_->SetBranchAddress("dtRechitCluster_match_MB1Seg_0p4",             dtRechitCluster_match_MB1Seg_0p4);
  tree_->SetBranchAddress("dtRechitCluster_match_MB1Seg_0p5",             dtRechitCluster_match_MB1Seg_0p5);
  tree_->SetBranchAddress("dtRechitCluster_match_RPChits_dPhi0p5",             dtRechitCluster_match_RPChits_dPhi0p5);
  tree_->SetBranchAddress("dtRechitCluster_match_RPCBx_dPhi0p5",             dtRechitCluster_match_RPCBx_dPhi0p5);
  tree_->SetBranchAddress("dtRechitCluster_match_RB1_0p4",             dtRechitCluster_match_RB1_0p4);
  tree_->SetBranchAddress("dtRechitCluster_match_RB1_dPhi0p5",             dtRechitCluster_match_RB1_dPhi0p5);

  tree_->SetBranchAddress("dtRechitCluster_match_dtSeg_0p5",             dtRechitCluster_match_dtSeg_0p5);
  tree_->SetBranchAddress("dtRechitCluster_match_dtSegTime_0p5",             dtRechitCluster_match_dtSegTime_0p5);
  tree_->SetBranchAddress("dtRechitCluster_match_dtSeg_0p4",             dtRechitCluster_match_dtSeg_0p4);
  tree_->SetBranchAddress("dtRechitCluster_match_dtSegTime_0p4",             dtRechitCluster_match_dtSegTime_0p4);
  tree_->SetBranchAddress("dtRechitCluster_match_dtSegTimeSpread_0p5",             dtRechitCluster_match_dtSegTimeSpread_0p5);
  tree_->SetBranchAddress("dtRechitCluster_match_dtSegTimeSpread_0p4",             dtRechitCluster_match_dtSegTimeSpread_0p4);


    tree_->SetBranchAddress("dtRechitCluster_match_dtSeg_sameStation_0p5",             dtRechitCluster_match_dtSeg_sameStation_0p5);
    tree_->SetBranchAddress("dtRechitCluster_match_dtSegTime_sameStation_0p5",             dtRechitCluster_match_dtSegTime_sameStation_0p5);
    tree_->SetBranchAddress("dtRechitCluster_match_dtSeg_sameStation_0p4",             dtRechitCluster_match_dtSeg_sameStation_0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_dtSegTime_sameStation_0p4",             dtRechitCluster_match_dtSegTime_sameStation_0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5",             dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5);
    tree_->SetBranchAddress("dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4",             dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4);


    tree_->SetBranchAddress("dtRechitCluster_match_RPCTime_dPhi0p5",             dtRechitCluster_match_RPCTime_dPhi0p5);
    tree_->SetBranchAddress("dtRechitCluster_match_RPCTimeSpread_dPhi0p5",             dtRechitCluster_match_RPCTimeSpread_dPhi0p5);
    tree_->SetBranchAddress("dtRechitCluster_match_RPCTime_dR0p4",             dtRechitCluster_match_RPCTime_dR0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_RPCTimeSpread_dR0p4",             dtRechitCluster_match_RPCTimeSpread_dR0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_RPChits_dR0p4",             dtRechitCluster_match_RPChits_dR0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_RPCTime_sameStation_dR0p4",             dtRechitCluster_match_RPCTime_sameStation_dR0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4",             dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_RPChits_sameStation_dR0p4",             dtRechitCluster_match_RPChits_sameStation_dR0p4);






  tree_->SetBranchAddress("dtRechitClusterMet_dPhi",             dtRechitClusterMet_dPhi);
  tree_->SetBranchAddress("dtRechitClusterMetXYCorr_dPhi",             dtRechitClusterMetXYCorr_dPhi);


  tree_->SetBranchAddress("dtRechitClusterMetHEM_dPhi",             dtRechitClusterMetHEM_dPhi);
  tree_->SetBranchAddress("dtRechitClusterMetHEMXYCorr_dPhi",             dtRechitClusterMetHEMXYCorr_dPhi);
  tree_->SetBranchAddress("dtRechitClusterMetEENoise_dPhi",             dtRechitClusterMetEENoise_dPhi);
  tree_->SetBranchAddress("dtRechitClusterMetEENoiseXYCorr_dPhi",             dtRechitClusterMetEENoiseXYCorr_dPhi);
  tree_->SetBranchAddress("dtRechitClusterMetJesUp_dPhi",             dtRechitClusterMetJesUp_dPhi);
  tree_->SetBranchAddress("dtRechitClusterMetJesDown_dPhi",             dtRechitClusterMetJesDown_dPhi);



  tree_->SetBranchAddress("dtRechitClusterMetJesDown_dPhi",             dtRechitClusterMetJesDown_dPhi);


          tree_->SetBranchAddress("dtRechitClusterZLep1",             dtRechitClusterZLep1);
          tree_->SetBranchAddress("dtRechitClusterZLep2",             dtRechitClusterZLep2);
          tree_->SetBranchAddress("dtRechitClusterZLep1Tag",          dtRechitClusterZLep1Tag);
          tree_->SetBranchAddress("dtRechitClusterZLep2Tag",          dtRechitClusterZLep2Tag);
          tree_->SetBranchAddress("dtRechitClusterZLep1Id",           dtRechitClusterZLep1Id);
          tree_->SetBranchAddress("dtRechitClusterZLep2Id",           dtRechitClusterZLep2Id);
          tree_->SetBranchAddress("dtRechitClusterZLep1LooseIso",     dtRechitClusterZLep1LooseIso);
          tree_->SetBranchAddress("dtRechitClusterZLep1TightIso",     dtRechitClusterZLep1TightIso);
          tree_->SetBranchAddress("dtRechitClusterZLep1VTightIso",    dtRechitClusterZLep1VTightIso);
          tree_->SetBranchAddress("dtRechitClusterZLep1VVTightIso",   dtRechitClusterZLep1VVTightIso);
          tree_->SetBranchAddress("dtRechitClusterZLep1TightId",      dtRechitClusterZLep1TightId);
          tree_->SetBranchAddress("dtRechitClusterZLep2LooseIso",     dtRechitClusterZLep2LooseIso);
          tree_->SetBranchAddress("dtRechitClusterZLep2TightIso",     dtRechitClusterZLep2TightIso);
          tree_->SetBranchAddress("dtRechitClusterZLep2VTightIso",    dtRechitClusterZLep2VTightIso);
          tree_->SetBranchAddress("dtRechitClusterZLep2VVTightIso",   dtRechitClusterZLep2VVTightIso);
          tree_->SetBranchAddress("dtRechitClusterZLep2TightId",      dtRechitClusterZLep2TightId);
    // vector<vector<float> > * dtRechitCluster_match_dtSegT_dR0p4 = 0;
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegT_dR0p4",             &dtRechitCluster_match_dtSegT_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegX_dR0p4",             dtRechitCluster_match_dtSegX_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegY_dR0p4",             dtRechitCluster_match_dtSegY_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegZ_dR0p4",             dtRechitCluster_match_dtSegZ_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegPhi_dR0p4",             dtRechitCluster_match_dtSegPhi_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegEta_dR0p4",             dtRechitCluster_match_dtSegEta_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegWheel_dR0p4",             dtRechitCluster_match_dtSegWheel_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegStation_dR0p4",             dtRechitCluster_match_dtSegStation_dR0p4);
      //
      //
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegT_sameStation_dR0p4",             dtRechitCluster_match_dtSegT_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegX_sameStation_dR0p4",             dtRechitCluster_match_dtSegX_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegY_sameStation_dR0p4",             dtRechitCluster_match_dtSegY_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegZ_sameStation_dR0p4",             dtRechitCluster_match_dtSegZ_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegPhi_sameStation_dR0p4",             dtRechitCluster_match_dtSegPhi_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegEta_sameStation_dR0p4",             dtRechitCluster_match_dtSegEta_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegWheel_sameStation_dR0p4",             dtRechitCluster_match_dtSegWheel_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_dtSegStation_sameStation_dR0p4",             dtRechitCluster_match_dtSegStation_sameStation_dR0p4);
      //
      //
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCBx_dPhi0p5",             dtRechitCluster_match_RPCBx_dPhi0p5);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCX_dPhi0p5",             dtRechitCluster_match_RPCX_dPhi0p5);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCY_dPhi0p5",             dtRechitCluster_match_RPCY_dPhi0p5);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCZ_dPhi0p5",             dtRechitCluster_match_RPCZ_dPhi0p5);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCPhi_dPhi0p5",             dtRechitCluster_match_RPCPhi_dPhi0p5);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCEta_dPhi0p5",             dtRechitCluster_match_RPCEta_dPhi0p5);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCRing_dPhi0p5",             dtRechitCluster_match_RPCRing_dPhi0p5);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCLayer_dPhi0p5",             dtRechitCluster_match_RPCLayer_dPhi0p5);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCSector_dPhi0p5",             dtRechitCluster_match_RPCSector_dPhi0p5);
      //
      //
      //
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCBx_dR0p4",             dtRechitCluster_match_RPCBx_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCX_dR0p4",             dtRechitCluster_match_RPCX_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCY_dR0p4",             dtRechitCluster_match_RPCY_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCZ_dR0p4",             dtRechitCluster_match_RPCZ_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCPhi_dR0p4",             dtRechitCluster_match_RPCPhi_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCEta_dR0p4",             dtRechitCluster_match_RPCEta_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCRing_dR0p4",             dtRechitCluster_match_RPCRing_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCLayer_dR0p4",             dtRechitCluster_match_RPCLayer_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCSector_dR0p4",             dtRechitCluster_match_RPCSector_dR0p4);
      //
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCBx_sameStation_dR0p4",             dtRechitCluster_match_RPCBx_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCX_sameStation_dR0p4",             dtRechitCluster_match_RPCX_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCY_sameStation_dR0p4",             dtRechitCluster_match_RPCY_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCZ_sameStation_dR0p4",             dtRechitCluster_match_RPCZ_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCPhi_sameStation_dR0p4",             dtRechitCluster_match_RPCPhi_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCEta_sameStation_dR0p4",             dtRechitCluster_match_RPCEta_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCRing_sameStation_dR0p4",             dtRechitCluster_match_RPCRing_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCLayer_sameStation_dR0p4",             dtRechitCluster_match_RPCLayer_sameStation_dR0p4);
      // tree_->SetBranchAddress("dtRechitCluster_match_RPCSector_sameStation_dR0p4",             dtRechitCluster_match_RPCSector_sameStation_dR0p4);


    //correct clusters



      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP",             &dtRechitCluster2_match_gLLP);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_index",             &dtRechitCluster2_match_gLLP_index);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_minDeltaR",             &dtRechitCluster2_match_gLLP_minDeltaR);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_eta",             &dtRechitCluster2_match_gLLP_eta);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_phi",             &dtRechitCluster2_match_gLLP_phi);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_decay_r",             &dtRechitCluster2_match_gLLP_decay_r);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_decay_x",             &dtRechitCluster2_match_gLLP_decay_x);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_decay_y",             &dtRechitCluster2_match_gLLP_decay_y);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_decay_z",             &dtRechitCluster2_match_gLLP_decay_z);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_ctau",             &dtRechitCluster2_match_gLLP_ctau);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_beta",             &dtRechitCluster2_match_gLLP_beta);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_csc",             &dtRechitCluster2_match_gLLP_csc);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_dt",             &dtRechitCluster2_match_gLLP_dt);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_multiplicity",             &dtRechitCluster2_match_gLLP_multiplicity);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_EM_multiplicity",             &dtRechitCluster2_match_gLLP_EM_multiplicity);

      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_e",             &dtRechitCluster2_match_gLLP_e);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_pt",             &dtRechitCluster2_match_gLLP_pt);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_EMFracE",             &dtRechitCluster2_match_gLLP_EMFracE);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_EMFracEz",             &dtRechitCluster2_match_gLLP_EMFracEz);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_EMFracP",             &dtRechitCluster2_match_gLLP_EMFracP);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_EMFracPz",             &dtRechitCluster2_match_gLLP_EMFracPz);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_lepdPhi",             &dtRechitCluster2_match_gLLP_lepdPhi);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_daughterKaon",             &dtRechitCluster2_match_gLLP_daughterKaon);

      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_daughter0_deltaR",             &dtRechitCluster2_match_gLLP_daughter0_deltaR);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_daughter1_deltaR",             &dtRechitCluster2_match_gLLP_daughter1_deltaR);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_daughter2_deltaR",             &dtRechitCluster2_match_gLLP_daughter2_deltaR);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_daughter3_deltaR",             &dtRechitCluster2_match_gLLP_daughter3_deltaR);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_daughter_deltaR",             &dtRechitCluster2_match_gLLP_other_daughter_deltaR);
      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_daughter_index",             &dtRechitCluster2_match_gLLP_other_daughter_index);

      tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_daughter_deltaR",             &dtRechitCluster2_match_gLLP_other_daughter_deltaR);

     tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_eta",             &dtRechitCluster2_match_gLLP_other_eta);
     tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_phi",             &dtRechitCluster2_match_gLLP_other_phi);
     tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_decay_r",             &dtRechitCluster2_match_gLLP_other_decay_r);
     tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_decay_x",             &dtRechitCluster2_match_gLLP_other_decay_x);
     tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_decay_y",             &dtRechitCluster2_match_gLLP_other_decay_y);
     tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_decay_z",             &dtRechitCluster2_match_gLLP_other_decay_z);
     tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_ctau",             &dtRechitCluster2_match_gLLP_other_ctau);
     tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_beta",             &dtRechitCluster2_match_gLLP_other_beta);
     tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_csc",             &dtRechitCluster2_match_gLLP_other_csc);
     tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_e",             &dtRechitCluster2_match_gLLP_other_e);
     tree_->SetBranchAddress("dtRechitCluster2_match_gLLP_other_pt",             &dtRechitCluster2_match_gLLP_other_pt);

      tree_->SetBranchAddress("dtRechitCluster2Me11Ratio",             &dtRechitCluster2Me11Ratio);
      tree_->SetBranchAddress("dtRechitCluster2Me12Ratio",             &dtRechitCluster2Me12Ratio);
      tree_->SetBranchAddress("dtRechitCluster2MaxStation",             &dtRechitCluster2MaxStation);
      tree_->SetBranchAddress("dtRechitCluster2MaxStationRatio",             &dtRechitCluster2MaxStationRatio);
      tree_->SetBranchAddress("dtRechitCluster2NStation",             &dtRechitCluster2NStation);
      tree_->SetBranchAddress("dtRechitCluster2NStation5",             &dtRechitCluster2NStation5);
      tree_->SetBranchAddress("dtRechitCluster2NStation10",             &dtRechitCluster2NStation10);
      tree_->SetBranchAddress("dtRechitCluster2NStation10perc",             &dtRechitCluster2NStation10perc);
      tree_->SetBranchAddress("dtRechitCluster2AvgStation",             &dtRechitCluster2AvgStation);
      tree_->SetBranchAddress("dtRechitCluster2AvgStation5",             &dtRechitCluster2AvgStation5);
      tree_->SetBranchAddress("dtRechitCluster2AvgStation10",             &dtRechitCluster2AvgStation10);
      tree_->SetBranchAddress("dtRechitCluster2AvgStation10perc",             &dtRechitCluster2AvgStation10perc);
      tree_->SetBranchAddress("dtRechitCluster2MaxChamber",             &dtRechitCluster2MaxChamber);
      tree_->SetBranchAddress("dtRechitCluster2MaxChamberRatio",             &dtRechitCluster2MaxChamberRatio);
      tree_->SetBranchAddress("dtRechitCluster2NChamber",             &dtRechitCluster2NChamber);
      tree_->SetBranchAddress("dtRechitCluster2X",             dtRechitCluster2X);
      tree_->SetBranchAddress("dtRechitCluster2Y",             dtRechitCluster2Y);
      tree_->SetBranchAddress("dtRechitCluster2Z",             dtRechitCluster2Z);
      tree_->SetBranchAddress("dtRechitCluster2Time",             dtRechitCluster2Time);
      tree_->SetBranchAddress("dtRechitCluster2TimeWire",             dtRechitCluster2TimeWire);
      tree_->SetBranchAddress("dtRechitCluster2TimeWirePruned",             dtRechitCluster2TimeWirePruned);
      tree_->SetBranchAddress("dtRechitCluster2TimeTotal",             dtRechitCluster2TimeTotal);
      tree_->SetBranchAddress("dtRechitCluster2Wheel",             dtRechitCluster2Wheel);

      tree_->SetBranchAddress("dtRechitCluster2GenMuonDeltaR",             dtRechitCluster2GenMuonDeltaR);

      tree_->SetBranchAddress("dtRechitCluster2TimeSpread",             dtRechitCluster2TimeSpread);
      tree_->SetBranchAddress("dtRechitCluster2TimeTotalSpread",             dtRechitCluster2TimeTotalSpread);
      tree_->SetBranchAddress("dtRechitCluster2TimeTotalSpreadPruned",             dtRechitCluster2TimeTotalSpreadPruned);
      tree_->SetBranchAddress("dtRechitCluster2TimeWireSpread",             dtRechitCluster2TimeWireSpread);
      tree_->SetBranchAddress("dtRechitCluster2MajorAxis",             dtRechitCluster2MajorAxis);
      tree_->SetBranchAddress("dtRechitCluster2MinorAxis",             dtRechitCluster2MinorAxis);
      tree_->SetBranchAddress("dtRechitCluster2RSpread",             dtRechitCluster2RSpread);

      tree_->SetBranchAddress("dtRechitCluster2XSpread",             dtRechitCluster2XSpread);
      tree_->SetBranchAddress("dtRechitCluster2YSpread",             dtRechitCluster2YSpread);
      tree_->SetBranchAddress("dtRechitCluster2ZSpread",             dtRechitCluster2ZSpread);
      tree_->SetBranchAddress("dtRechitCluster2EtaPhiSpread",             dtRechitCluster2EtaPhiSpread);
      tree_->SetBranchAddress("dtRechitCluster2XYSpread",             dtRechitCluster2XYSpread);
      tree_->SetBranchAddress("dtRechitCluster2EtaSpread",             dtRechitCluster2EtaSpread);
      tree_->SetBranchAddress("dtRechitCluster2PhiSpread",             dtRechitCluster2PhiSpread);
      tree_->SetBranchAddress("dtRechitCluster2DeltaRSpread",             dtRechitCluster2DeltaRSpread);
      tree_->SetBranchAddress("dtRechitCluster2Eta",             dtRechitCluster2Eta);
      tree_->SetBranchAddress("dtRechitCluster2Phi",             dtRechitCluster2Phi);


      tree_->SetBranchAddress("dtRechitCluster2GenJetVetoPt",             dtRechitCluster2GenJetVetoPt);
      tree_->SetBranchAddress("dtRechitCluster2GenJetVetoE",             dtRechitCluster2GenJetVetoE);
      tree_->SetBranchAddress("dtRechitCluster2JetVetoPt",             dtRechitCluster2JetVetoPt);
      tree_->SetBranchAddress("dtRechitCluster2JetVetoEta",             dtRechitCluster2JetVetoEta);
      tree_->SetBranchAddress("dtRechitCluster2JetVetoPhi",             dtRechitCluster2JetVetoPhi);

      tree_->SetBranchAddress("dtRechitCluster2JetVetoElectronEnergyFraction",             dtRechitCluster2JetVetoElectronEnergyFraction);
      tree_->SetBranchAddress("dtRechitCluster2JetVetoPhotonEnergyFraction",             dtRechitCluster2JetVetoPhotonEnergyFraction);
      tree_->SetBranchAddress("dtRechitCluster2JetVetoNeutralHadronEnergyFraction",             dtRechitCluster2JetVetoNeutralHadronEnergyFraction);
      tree_->SetBranchAddress("dtRechitCluster2JetVetoChargedHadronEnergyFraction",             dtRechitCluster2JetVetoChargedHadronEnergyFraction);
      tree_->SetBranchAddress("dtRechitCluster2JetVetoMuonEnergyFraction",             dtRechitCluster2JetVetoMuonEnergyFraction);


      tree_->SetBranchAddress("dtRechitCluster2JetVetoE",             dtRechitCluster2JetVetoE);
      tree_->SetBranchAddress("dtRechitCluster2MuonVetoPt",             dtRechitCluster2MuonVetoPt);
      tree_->SetBranchAddress("dtRechitCluster2MuonVetoE",             dtRechitCluster2MuonVetoE);

      tree_->SetBranchAddress("dtRechitCluster2JetVetoPt_0p6",             dtRechitCluster2JetVetoPt_0p6);
      tree_->SetBranchAddress("dtRechitCluster2JetVetoPt_0p8",             dtRechitCluster2JetVetoPt_0p8);
      tree_->SetBranchAddress("dtRechitCluster2JetVetoE_0p6",             dtRechitCluster2JetVetoE_0p6);
      tree_->SetBranchAddress("dtRechitCluster2JetVetoE_0p8",             dtRechitCluster2JetVetoE_0p8);
      tree_->SetBranchAddress("dtRechitCluster2MuonVetoPt_0p6",             dtRechitCluster2MuonVetoPt_0p6);
      tree_->SetBranchAddress("dtRechitCluster2MuonVetoPt_0p8",             dtRechitCluster2MuonVetoPt_0p8);
      tree_->SetBranchAddress("dtRechitCluster2MuonVetoE_0p6",             dtRechitCluster2MuonVetoE_0p6);
      tree_->SetBranchAddress("dtRechitCluster2MuonVetoE_0p8",             dtRechitCluster2MuonVetoE_0p8);

      // tree_->SetBranchAddress("dtRechitCluster2ZLep1",             dtRechitCluster2ZLep1);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep2",             dtRechitCluster2ZLep2);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep1Tag",             dtRechitCluster2ZLep1Tag);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep2Tag",             dtRechitCluster2ZLep2Tag);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep1Id",             dtRechitCluster2ZLep1Id);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep2Id",             dtRechitCluster2ZLep2Id);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep1LooseIso",             dtRechitCluster2ZLep1LooseIso);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep1TightIso",             dtRechitCluster2ZLep1TightIso);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep1VTightIso",             dtRechitCluster2ZLep1VTightIso);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep1VVTightIso",             dtRechitCluster2ZLep1VVTightIso);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep1TightId",             dtRechitCluster2ZLep1TightId);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep2LooseIso",             dtRechitCluster2ZLep2LooseIso);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep2TightIso",             dtRechitCluster2ZLep2TightIso);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep2VTightIso",             dtRechitCluster2ZLep2VTightIso);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep2VVTightIso",             dtRechitCluster2ZLep2VVTightIso);
      // tree_->SetBranchAddress("dtRechitCluster2ZLep2TightId",             dtRechitCluster2ZLep2TightId);

      tree_->SetBranchAddress("dtRechitCluster2MuonVetoPhi",             dtRechitCluster2MuonVetoPhi);
      tree_->SetBranchAddress("dtRechitCluster2MuonVetoEta",             dtRechitCluster2MuonVetoEta);
      tree_->SetBranchAddress("dtRechitCluster2MuonVetoIso",             dtRechitCluster2MuonVetoIso);

        tree_->SetBranchAddress("dtRechitCluster2MuonVetoLooseIso",             dtRechitCluster2MuonVetoLooseIso);
        tree_->SetBranchAddress("dtRechitCluster2MuonVetoTightIso",             dtRechitCluster2MuonVetoTightIso);
        tree_->SetBranchAddress("dtRechitCluster2MuonVetoVTightIso",             dtRechitCluster2MuonVetoVTightIso);
        tree_->SetBranchAddress("dtRechitCluster2MuonVetoVVTightIso",             dtRechitCluster2MuonVetoVVTightIso);
        tree_->SetBranchAddress("dtRechitCluster2MuonVetoTightId",             dtRechitCluster2MuonVetoTightId);
        tree_->SetBranchAddress("dtRechitCluster2MuonVetoLooseId",             dtRechitCluster2MuonVetoLooseId);
        tree_->SetBranchAddress("dtRechitCluster2MuonVetoGlobal",             dtRechitCluster2MuonVetoGlobal);


      tree_->SetBranchAddress("dtRechitCluster2IsoMuonVetoPt",             dtRechitCluster2IsoMuonVetoPt);
      tree_->SetBranchAddress("dtRechitCluster2IsoMuonVetoE",             dtRechitCluster2IsoMuonVetoE);
      tree_->SetBranchAddress("dtRechitCluster2IsoMuonVetoPhi",             dtRechitCluster2IsoMuonVetoPhi);
      tree_->SetBranchAddress("dtRechitCluster2IsoMuonVetoEta",             dtRechitCluster2IsoMuonVetoEta);
      tree_->SetBranchAddress("dtRechitCluster2GenMuonVetoPt",             dtRechitCluster2GenMuonVetoPt);
      tree_->SetBranchAddress("dtRechitCluster2MuonVetoType",             dtRechitCluster2MuonVetoType);

      // tree_->SetBranchAddress("dtRechitCluster2GenMuonVetoProdX",             dtRechitCluster2GenMuonVetoProdX);
      // tree_->SetBranchAddress("dtRechitCluster2GenMuonVetoProdY",             dtRechitCluster2GenMuonVetoProdY);
      // tree_->SetBranchAddress("dtRechitCluster2GenMuonVetoProdZ",             dtRechitCluster2GenMuonVetoProdZ);
      // tree_->SetBranchAddress("dtRechitCluster2GenMuonVetoLLPDist",             dtRechitCluster2GenMuonVetoLLPDist);
      // tree_->SetBranchAddress("dtRechitCluster2GenMuonVetoLLPIndex",             dtRechitCluster2GenMuonVetoLLPIndex);

      tree_->SetBranchAddress("dtRechitCluster2GenMuonVetoE",             dtRechitCluster2GenMuonVetoE);
      tree_->SetBranchAddress("dtRechitCluster2Size",             dtRechitCluster2Size);

      tree_->SetBranchAddress("dtRechitCluster2NSegmentStation1",             dtRechitCluster2NSegmentStation1);
      tree_->SetBranchAddress("dtRechitCluster2NSegmentStation2",             dtRechitCluster2NSegmentStation2);
      tree_->SetBranchAddress("dtRechitCluster2NSegmentStation3",             dtRechitCluster2NSegmentStation3);
      tree_->SetBranchAddress("dtRechitCluster2NSegmentStation4",             dtRechitCluster2NSegmentStation4);



      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberPlus11",             dtRechitCluster2NRechitChamberPlus11);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberPlus12",             dtRechitCluster2NRechitChamberPlus12);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberPlus13",             dtRechitCluster2NRechitChamberPlus13);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberPlus21",             dtRechitCluster2NRechitChamberPlus21);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberPlus22",             dtRechitCluster2NRechitChamberPlus22);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberPlus31",             dtRechitCluster2NRechitChamberPlus31);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberPlus32",             dtRechitCluster2NRechitChamberPlus32);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberPlus41",             dtRechitCluster2NRechitChamberPlus41);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberPlus42",             dtRechitCluster2NRechitChamberPlus42);

      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberMinus11",             dtRechitCluster2NRechitChamberMinus11);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberMinus12",             dtRechitCluster2NRechitChamberMinus12);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberMinus13",             dtRechitCluster2NRechitChamberMinus13);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberMinus21",             dtRechitCluster2NRechitChamberMinus21);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberMinus22",             dtRechitCluster2NRechitChamberMinus22);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberMinus31",             dtRechitCluster2NRechitChamberMinus31);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberMinus32",             dtRechitCluster2NRechitChamberMinus32);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberMinus41",             dtRechitCluster2NRechitChamberMinus41);
      tree_->SetBranchAddress("dtRechitCluster2NRechitChamberMinus42",             dtRechitCluster2NRechitChamberMinus42);

      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberPlus11",             dtRechitCluster2NLayersChamberPlus11);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberPlus12",             dtRechitCluster2NLayersChamberPlus12);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberPlus13",             dtRechitCluster2NLayersChamberPlus13);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberPlus21",             dtRechitCluster2NLayersChamberPlus21);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberPlus22",             dtRechitCluster2NLayersChamberPlus22);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberPlus31",             dtRechitCluster2NLayersChamberPlus31);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberPlus32",             dtRechitCluster2NLayersChamberPlus32);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberPlus41",             dtRechitCluster2NLayersChamberPlus41);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberPlus42",             dtRechitCluster2NLayersChamberPlus42);

      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberMinus11",             dtRechitCluster2NLayersChamberMinus11);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberMinus12",             dtRechitCluster2NLayersChamberMinus12);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberMinus13",             dtRechitCluster2NLayersChamberMinus13);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberMinus21",             dtRechitCluster2NLayersChamberMinus21);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberMinus22",             dtRechitCluster2NLayersChamberMinus22);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberMinus31",             dtRechitCluster2NLayersChamberMinus31);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberMinus32",             dtRechitCluster2NLayersChamberMinus32);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberMinus41",             dtRechitCluster2NLayersChamberMinus41);
      tree_->SetBranchAddress("dtRechitCluster2NLayersChamberMinus42",             dtRechitCluster2NLayersChamberMinus42);




      tree_->SetBranchAddress("dtRechitCluster2_match_gParticle_Id",             dtRechitCluster2_match_gParticle_Id);
      tree_->SetBranchAddress("dtRechitCluster2_match_gParticle_Pt",             dtRechitCluster2_match_gParticle_Pt);
      tree_->SetBranchAddress("dtRechitCluster2_match_gParticle_Eta",             dtRechitCluster2_match_gParticle_Eta);
      tree_->SetBranchAddress("dtRechitCluster2_match_gParticle_Phi",             dtRechitCluster2_match_gParticle_Phi);
      tree_->SetBranchAddress("dtRechitCluster2_match_gParticle_E",             dtRechitCluster2_match_gParticle_E);
      tree_->SetBranchAddress("dtRechitCluster2_match_gParticle_Status",             dtRechitCluster2_match_gParticle_Status);
      tree_->SetBranchAddress("dtRechitCluster2_match_gParticle_MotherId",             dtRechitCluster2_match_gParticle_MotherId);
      tree_->SetBranchAddress("dtRechitCluster2_match_gParticle_deltaR",             dtRechitCluster2_match_gParticle_deltaR);


      tree_->SetBranchAddress("dtRechitCluster2_match_MB1hits_0p4",             dtRechitCluster2_match_MB1hits_0p4);
      tree_->SetBranchAddress("dtRechitCluster2_match_MB1hits_0p5",             dtRechitCluster2_match_MB1hits_0p5);
      tree_->SetBranchAddress("dtRechitCluster2_match_MB1hits_cosmics_plus",             dtRechitCluster2_match_MB1hits_cosmics_plus);
      tree_->SetBranchAddress("dtRechitCluster2_match_MB1hits_cosmics_minus",             dtRechitCluster2_match_MB1hits_cosmics_minus);
      tree_->SetBranchAddress("dtRechitCluster2_match_MB1Seg_0p4",             dtRechitCluster2_match_MB1Seg_0p4);
      tree_->SetBranchAddress("dtRechitCluster2_match_MB1Seg_0p5",             dtRechitCluster2_match_MB1Seg_0p5);
      tree_->SetBranchAddress("dtRechitCluster2_match_RPChits_dPhi0p5",             dtRechitCluster2_match_RPChits_dPhi0p5);
      tree_->SetBranchAddress("dtRechitCluster2_match_RPCBx_dPhi0p5",             dtRechitCluster2_match_RPCBx_dPhi0p5);
      tree_->SetBranchAddress("dtRechitCluster2_match_RB1_0p4",             dtRechitCluster2_match_RB1_0p4);
      tree_->SetBranchAddress("dtRechitCluster2_match_RB1_dPhi0p5",             dtRechitCluster2_match_RB1_dPhi0p5);

      tree_->SetBranchAddress("dtRechitCluster2_match_dtSeg_0p5",             dtRechitCluster2_match_dtSeg_0p5);
      tree_->SetBranchAddress("dtRechitCluster2_match_dtSegTime_0p5",             dtRechitCluster2_match_dtSegTime_0p5);
      tree_->SetBranchAddress("dtRechitCluster2_match_dtSeg_0p4",             dtRechitCluster2_match_dtSeg_0p4);
      tree_->SetBranchAddress("dtRechitCluster2_match_dtSegTime_0p4",             dtRechitCluster2_match_dtSegTime_0p4);
      tree_->SetBranchAddress("dtRechitCluster2_match_dtSegTimeSpread_0p5",             dtRechitCluster2_match_dtSegTimeSpread_0p5);
      tree_->SetBranchAddress("dtRechitCluster2_match_dtSegTimeSpread_0p4",             dtRechitCluster2_match_dtSegTimeSpread_0p4);


        tree_->SetBranchAddress("dtRechitCluster2_match_dtSeg_sameStation_0p5",             dtRechitCluster2_match_dtSeg_sameStation_0p5);
        tree_->SetBranchAddress("dtRechitCluster2_match_dtSegTime_sameStation_0p5",             dtRechitCluster2_match_dtSegTime_sameStation_0p5);
        tree_->SetBranchAddress("dtRechitCluster2_match_dtSeg_sameStation_0p4",             dtRechitCluster2_match_dtSeg_sameStation_0p4);
        tree_->SetBranchAddress("dtRechitCluster2_match_dtSegTime_sameStation_0p4",             dtRechitCluster2_match_dtSegTime_sameStation_0p4);
        tree_->SetBranchAddress("dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p5",             dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p5);
        tree_->SetBranchAddress("dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p4",             dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p4);


        tree_->SetBranchAddress("dtRechitCluster2_match_RPCTime_dPhi0p5",             dtRechitCluster2_match_RPCTime_dPhi0p5);
        tree_->SetBranchAddress("dtRechitCluster2_match_RPCTimeSpread_dPhi0p5",             dtRechitCluster2_match_RPCTimeSpread_dPhi0p5);
        tree_->SetBranchAddress("dtRechitCluster2_match_RPCTime_dR0p4",             dtRechitCluster2_match_RPCTime_dR0p4);
        tree_->SetBranchAddress("dtRechitCluster2_match_RPCTimeSpread_dR0p4",             dtRechitCluster2_match_RPCTimeSpread_dR0p4);
        tree_->SetBranchAddress("dtRechitCluster2_match_RPChits_dR0p4",             dtRechitCluster2_match_RPChits_dR0p4);
        tree_->SetBranchAddress("dtRechitCluster2_match_RPCTime_sameStation_dR0p4",             dtRechitCluster2_match_RPCTime_sameStation_dR0p4);
        tree_->SetBranchAddress("dtRechitCluster2_match_RPCTimeSpread_sameStation_dR0p4",             dtRechitCluster2_match_RPCTimeSpread_sameStation_dR0p4);
        tree_->SetBranchAddress("dtRechitCluster2_match_RPChits_sameStation_dR0p4",             dtRechitCluster2_match_RPChits_sameStation_dR0p4);






      tree_->SetBranchAddress("dtRechitCluster2Met_dPhi",             dtRechitCluster2Met_dPhi);
      tree_->SetBranchAddress("dtRechitCluster2MetXYCorr_dPhi",             dtRechitCluster2MetXYCorr_dPhi);


      tree_->SetBranchAddress("dtRechitCluster2MetHEM_dPhi",             dtRechitCluster2MetHEM_dPhi);
      tree_->SetBranchAddress("dtRechitCluster2MetHEMXYCorr_dPhi",             dtRechitCluster2MetHEMXYCorr_dPhi);
      tree_->SetBranchAddress("dtRechitCluster2MetEENoise_dPhi",             dtRechitCluster2MetEENoise_dPhi);
      tree_->SetBranchAddress("dtRechitCluster2MetEENoiseXYCorr_dPhi",             dtRechitCluster2MetEENoiseXYCorr_dPhi);
      tree_->SetBranchAddress("dtRechitCluster2MetJesUp_dPhi",             dtRechitCluster2MetJesUp_dPhi);
      tree_->SetBranchAddress("dtRechitCluster2MetJesDown_dPhi",             dtRechitCluster2MetJesDown_dPhi);



      tree_->SetBranchAddress("dtRechitCluster2MetJesDown_dPhi",             dtRechitCluster2MetJesDown_dPhi);


  tree_->SetBranchAddress("nCscRechitClusters",             &nCscRechitClusters);
  tree_->SetBranchAddress("cscRechitCluster_match_Me1112_0p4",             &cscRechitCluster_match_Me1112_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_Me1112_0p6",             &cscRechitCluster_match_Me1112_0p6);
  tree_->SetBranchAddress("cscRechitCluster_match_Me1112_0p8",             &cscRechitCluster_match_Me1112_0p8);
  tree_->SetBranchAddress("cscRechitCluster_match_Me11_0p4",             &cscRechitCluster_match_Me11_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_Me11_0p6",             &cscRechitCluster_match_Me11_0p6);
  tree_->SetBranchAddress("cscRechitCluster_match_Me11_0p8",             &cscRechitCluster_match_Me11_0p8);
  tree_->SetBranchAddress("cscRechitCluster_match_Me12_0p4",             &cscRechitCluster_match_Me12_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_Me12_0p6",             &cscRechitCluster_match_Me12_0p6);
  tree_->SetBranchAddress("cscRechitCluster_match_Me12_0p8",             &cscRechitCluster_match_Me12_0p8);
  // tree_->SetBranchAddress("cscRechitCluster_match_gParticle_id",             &cscRechitCluster_match_gParticle_id);
  // tree_->SetBranchAddress("cscRechitCluster_match_gParticle",             &cscRechitCluster_match_gParticle);
  // tree_->SetBranchAddress("cscRechitCluster_match_gParticle_minDeltaR",             &cscRechitCluster_match_gParticle_minDeltaR);
  // tree_->SetBranchAddress("cscRechitCluster_match_gParticle_index",             &cscRechitCluster_match_gParticle_index);
  // tree_->SetBranchAddress("cscRechitCluster_match_gParticle_eta",             &cscRechitCluster_match_gParticle_eta);
  // tree_->SetBranchAddress("cscRechitCluster_match_gParticle_phi",             &cscRechitCluster_match_gParticle_phi);
  // tree_->SetBranchAddress("cscRechitCluster_match_gParticle_E",             &cscRechitCluster_match_gParticle_E);
  // tree_->SetBranchAddress("cscRechitCluster_match_gParticle_pt",             &cscRechitCluster_match_gParticle_pt);
  // tree_->SetBranchAddress("cscRechitCluster_match_gParticle_MotherId",             &cscRechitCluster_match_gParticle_MotherId);

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
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_dt",             &cscRechitCluster_match_gLLP_dt);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_multiplicity",             &cscRechitCluster_match_gLLP_multiplicity);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_EM_multiplicity",             &cscRechitCluster_match_gLLP_EM_multiplicity);

  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_e",             &cscRechitCluster_match_gLLP_e);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_pt",             &cscRechitCluster_match_gLLP_pt);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_EMFracE",             &cscRechitCluster_match_gLLP_EMFracE);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_EMFracEz",             &cscRechitCluster_match_gLLP_EMFracEz);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_EMFracP",             &cscRechitCluster_match_gLLP_EMFracP);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_EMFracPz",             &cscRechitCluster_match_gLLP_EMFracPz);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_lepdPhi",             &cscRechitCluster_match_gLLP_lepdPhi);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_daughterKaon",             &cscRechitCluster_match_gLLP_daughterKaon);

  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_daughter0_deltaR",             &cscRechitCluster_match_gLLP_daughter0_deltaR);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_daughter1_deltaR",             &cscRechitCluster_match_gLLP_daughter1_deltaR);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_daughter2_deltaR",             &cscRechitCluster_match_gLLP_daughter2_deltaR);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_daughter3_deltaR",             &cscRechitCluster_match_gLLP_daughter3_deltaR);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_daughter_deltaR",             &cscRechitCluster_match_gLLP_other_daughter_deltaR);
  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_daughter_index",             &cscRechitCluster_match_gLLP_other_daughter_index);

  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_daughter_deltaR",             &cscRechitCluster_match_gLLP_other_daughter_deltaR);

 tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_eta",             &cscRechitCluster_match_gLLP_other_eta);
 tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_phi",             &cscRechitCluster_match_gLLP_other_phi);
 tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_decay_r",             &cscRechitCluster_match_gLLP_other_decay_r);
 tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_decay_x",             &cscRechitCluster_match_gLLP_other_decay_x);
 tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_decay_y",             &cscRechitCluster_match_gLLP_other_decay_y);
 tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_decay_z",             &cscRechitCluster_match_gLLP_other_decay_z);
 tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_ctau",             &cscRechitCluster_match_gLLP_other_ctau);
 tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_beta",             &cscRechitCluster_match_gLLP_other_beta);
 tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_csc",             &cscRechitCluster_match_gLLP_other_csc);
 tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_e",             &cscRechitCluster_match_gLLP_other_e);
 tree_->SetBranchAddress("cscRechitCluster_match_gLLP_other_pt",             &cscRechitCluster_match_gLLP_other_pt);

  tree_->SetBranchAddress("cscRechitClusterMe11Ratio",             &cscRechitClusterMe11Ratio);
  tree_->SetBranchAddress("cscRechitClusterMe12Ratio",             &cscRechitClusterMe12Ratio);
  tree_->SetBranchAddress("cscRechitClusterMaxStation",             &cscRechitClusterMaxStation);
  tree_->SetBranchAddress("cscRechitClusterMaxStationRatio",             &cscRechitClusterMaxStationRatio);
  tree_->SetBranchAddress("cscRechitClusterNStation",             &cscRechitClusterNStation);
  tree_->SetBranchAddress("cscRechitClusterNStation5",             &cscRechitClusterNStation5);
  tree_->SetBranchAddress("cscRechitClusterNStation10",             &cscRechitClusterNStation10);
  tree_->SetBranchAddress("cscRechitClusterNStation10perc",             &cscRechitClusterNStation10perc);
  tree_->SetBranchAddress("cscRechitClusterAvgStation",             &cscRechitClusterAvgStation);
  tree_->SetBranchAddress("cscRechitClusterAvgStation5",             &cscRechitClusterAvgStation5);
  tree_->SetBranchAddress("cscRechitClusterAvgStation10",             &cscRechitClusterAvgStation10);
  tree_->SetBranchAddress("cscRechitClusterAvgStation10perc",             &cscRechitClusterAvgStation10perc);
  tree_->SetBranchAddress("cscRechitClusterMaxChamber",             &cscRechitClusterMaxChamber);
  tree_->SetBranchAddress("cscRechitClusterMaxChamberRatio",             &cscRechitClusterMaxChamberRatio);
  tree_->SetBranchAddress("cscRechitClusterNChamber",             &cscRechitClusterNChamber);
  tree_->SetBranchAddress("cscRechitClusterX",             cscRechitClusterX);
  tree_->SetBranchAddress("cscRechitClusterY",             cscRechitClusterY);
  tree_->SetBranchAddress("cscRechitClusterZ",             cscRechitClusterZ);
  tree_->SetBranchAddress("cscRechitClusterTime",             cscRechitClusterTime);
  tree_->SetBranchAddress("cscRechitClusterTimeWire",             cscRechitClusterTimeWire);
  tree_->SetBranchAddress("cscRechitClusterTimeWirePruned",             cscRechitClusterTimeWirePruned);
  tree_->SetBranchAddress("cscRechitClusterTimeTotal",             cscRechitClusterTimeTotal);

  tree_->SetBranchAddress("cscRechitClusterGenMuonDeltaR",             cscRechitClusterGenMuonDeltaR);

  tree_->SetBranchAddress("cscRechitClusterTimeSpread",             cscRechitClusterTimeSpread);
  tree_->SetBranchAddress("cscRechitClusterTimeTotalSpread",             cscRechitClusterTimeTotalSpread);
  tree_->SetBranchAddress("cscRechitClusterTimeTotalSpreadPruned",             cscRechitClusterTimeTotalSpreadPruned);
  tree_->SetBranchAddress("cscRechitClusterTimeWireSpread",             cscRechitClusterTimeWireSpread);
  tree_->SetBranchAddress("cscRechitClusterMajorAxis",             cscRechitClusterMajorAxis);
  tree_->SetBranchAddress("cscRechitClusterMinorAxis",             cscRechitClusterMinorAxis);
  tree_->SetBranchAddress("cscRechitClusterRSpread",             cscRechitClusterRSpread);

  tree_->SetBranchAddress("cscRechitClusterXSpread",             cscRechitClusterXSpread);
  tree_->SetBranchAddress("cscRechitClusterYSpread",             cscRechitClusterYSpread);
  tree_->SetBranchAddress("cscRechitClusterZSpread",             cscRechitClusterZSpread);
  tree_->SetBranchAddress("cscRechitClusterEtaPhiSpread",             cscRechitClusterEtaPhiSpread);
  tree_->SetBranchAddress("cscRechitClusterXYSpread",             cscRechitClusterXYSpread);
  tree_->SetBranchAddress("cscRechitClusterEtaSpread",             cscRechitClusterEtaSpread);
  tree_->SetBranchAddress("cscRechitClusterPhiSpread",             cscRechitClusterPhiSpread);
  tree_->SetBranchAddress("cscRechitClusterDeltaRSpread",             cscRechitClusterDeltaRSpread);
  tree_->SetBranchAddress("cscRechitClusterEta",             cscRechitClusterEta);
  tree_->SetBranchAddress("cscRechitClusterPhi",             cscRechitClusterPhi);


  tree_->SetBranchAddress("cscRechitClusterGenJetVetoPt",             cscRechitClusterGenJetVetoPt);
  tree_->SetBranchAddress("cscRechitClusterGenJetVetoE",             cscRechitClusterGenJetVetoE);
  tree_->SetBranchAddress("cscRechitClusterJetVetoPt",             cscRechitClusterJetVetoPt);
  tree_->SetBranchAddress("cscRechitClusterJetVetoEta",             cscRechitClusterJetVetoEta);
  tree_->SetBranchAddress("cscRechitClusterJetVetoPhi",             cscRechitClusterJetVetoPhi);

  tree_->SetBranchAddress("cscRechitClusterJetVetoElectronEnergyFraction",             cscRechitClusterJetVetoElectronEnergyFraction);
  tree_->SetBranchAddress("cscRechitClusterJetVetoPhotonEnergyFraction",             cscRechitClusterJetVetoPhotonEnergyFraction);
  tree_->SetBranchAddress("cscRechitClusterJetVetoNeutralHadronEnergyFraction",             cscRechitClusterJetVetoNeutralHadronEnergyFraction);
  tree_->SetBranchAddress("cscRechitClusterJetVetoChargedHadronEnergyFraction",             cscRechitClusterJetVetoChargedHadronEnergyFraction);
  tree_->SetBranchAddress("cscRechitClusterJetVetoMuonEnergyFraction",             cscRechitClusterJetVetoMuonEnergyFraction);


  tree_->SetBranchAddress("cscRechitClusterJetVetoE",             cscRechitClusterJetVetoE);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoPt",             cscRechitClusterMuonVetoPt);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoE",             cscRechitClusterMuonVetoE);

  tree_->SetBranchAddress("cscRechitClusterJetVetoPt_0p6",             cscRechitClusterJetVetoPt_0p6);
  tree_->SetBranchAddress("cscRechitClusterJetVetoPt_0p8",             cscRechitClusterJetVetoPt_0p8);
  tree_->SetBranchAddress("cscRechitClusterJetVetoE_0p6",             cscRechitClusterJetVetoE_0p6);
  tree_->SetBranchAddress("cscRechitClusterJetVetoE_0p8",             cscRechitClusterJetVetoE_0p8);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoPt_0p6",             cscRechitClusterMuonVetoPt_0p6);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoPt_0p8",             cscRechitClusterMuonVetoPt_0p8);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoE_0p6",             cscRechitClusterMuonVetoE_0p6);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoE_0p8",             cscRechitClusterMuonVetoE_0p8);

  tree_->SetBranchAddress("cscRechitClusterZLep1",             cscRechitClusterZLep1);
  tree_->SetBranchAddress("cscRechitClusterZLep2",             cscRechitClusterZLep2);
  tree_->SetBranchAddress("cscRechitClusterZLep1Tag",             cscRechitClusterZLep1Tag);
  tree_->SetBranchAddress("cscRechitClusterZLep2Tag",             cscRechitClusterZLep2Tag);
  tree_->SetBranchAddress("cscRechitClusterZLep1Id",             cscRechitClusterZLep1Id);
  tree_->SetBranchAddress("cscRechitClusterZLep2Id",             cscRechitClusterZLep2Id);
  tree_->SetBranchAddress("cscRechitClusterZLep1LooseIso",             cscRechitClusterZLep1LooseIso);
  tree_->SetBranchAddress("cscRechitClusterZLep1TightIso",             cscRechitClusterZLep1TightIso);
  tree_->SetBranchAddress("cscRechitClusterZLep1VTightIso",             cscRechitClusterZLep1VTightIso);
  tree_->SetBranchAddress("cscRechitClusterZLep1VVTightIso",             cscRechitClusterZLep1VVTightIso);
  tree_->SetBranchAddress("cscRechitClusterZLep1TightId",             cscRechitClusterZLep1TightId);
  tree_->SetBranchAddress("cscRechitClusterZLep2LooseIso",             cscRechitClusterZLep2LooseIso);
  tree_->SetBranchAddress("cscRechitClusterZLep2TightIso",             cscRechitClusterZLep2TightIso);
  tree_->SetBranchAddress("cscRechitClusterZLep2VTightIso",             cscRechitClusterZLep2VTightIso);
  tree_->SetBranchAddress("cscRechitClusterZLep2VVTightIso",             cscRechitClusterZLep2VVTightIso);
  tree_->SetBranchAddress("cscRechitClusterZLep2TightId",             cscRechitClusterZLep2TightId);

  tree_->SetBranchAddress("cscRechitClusterMuonVetoPhi",             cscRechitClusterMuonVetoPhi);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoEta",             cscRechitClusterMuonVetoEta);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoIso",             cscRechitClusterMuonVetoIso);

    tree_->SetBranchAddress("cscRechitClusterMuonVetoLooseIso",             cscRechitClusterMuonVetoLooseIso);
    tree_->SetBranchAddress("cscRechitClusterMuonVetoTightIso",             cscRechitClusterMuonVetoTightIso);
    tree_->SetBranchAddress("cscRechitClusterMuonVetoVTightIso",             cscRechitClusterMuonVetoVTightIso);
    tree_->SetBranchAddress("cscRechitClusterMuonVetoVVTightIso",             cscRechitClusterMuonVetoVVTightIso);
    tree_->SetBranchAddress("cscRechitClusterMuonVetoTightId",             cscRechitClusterMuonVetoTightId);
    tree_->SetBranchAddress("cscRechitClusterMuonVetoLooseId",             cscRechitClusterMuonVetoLooseId);
    tree_->SetBranchAddress("cscRechitClusterMuonVetoGlobal",             cscRechitClusterMuonVetoGlobal);


  tree_->SetBranchAddress("cscRechitClusterIsoMuonVetoPt",             cscRechitClusterIsoMuonVetoPt);
  tree_->SetBranchAddress("cscRechitClusterIsoMuonVetoE",             cscRechitClusterIsoMuonVetoE);
  tree_->SetBranchAddress("cscRechitClusterIsoMuonVetoPhi",             cscRechitClusterIsoMuonVetoPhi);
  tree_->SetBranchAddress("cscRechitClusterIsoMuonVetoEta",             cscRechitClusterIsoMuonVetoEta);
  tree_->SetBranchAddress("cscRechitClusterGenMuonVetoPt",             cscRechitClusterGenMuonVetoPt);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoType",             cscRechitClusterMuonVetoType);
  // tree_->SetBranchAddress("cscRechitClusterGenMuonVetoProdX",             cscRechitClusterGenMuonVetoProdX);
  // tree_->SetBranchAddress("cscRechitClusterGenMuonVetoProdY",             cscRechitClusterGenMuonVetoProdY);
  // tree_->SetBranchAddress("cscRechitClusterGenMuonVetoProdZ",             cscRechitClusterGenMuonVetoProdZ);
  // tree_->SetBranchAddress("cscRechitClusterGenMuonVetoLLPDist",             cscRechitClusterGenMuonVetoLLPDist);
  // tree_->SetBranchAddress("cscRechitClusterGenMuonVetoLLPIndex",             cscRechitClusterGenMuonVetoLLPIndex);

  tree_->SetBranchAddress("cscRechitClusterGenMuonVetoE",             cscRechitClusterGenMuonVetoE);
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

  tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus11",             cscRechitClusterNLayersChamberPlus11);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus12",             cscRechitClusterNLayersChamberPlus12);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus13",             cscRechitClusterNLayersChamberPlus13);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus21",             cscRechitClusterNLayersChamberPlus21);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus22",             cscRechitClusterNLayersChamberPlus22);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus31",             cscRechitClusterNLayersChamberPlus31);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus32",             cscRechitClusterNLayersChamberPlus32);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus41",             cscRechitClusterNLayersChamberPlus41);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus42",             cscRechitClusterNLayersChamberPlus42);

  tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus11",             cscRechitClusterNLayersChamberMinus11);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus12",             cscRechitClusterNLayersChamberMinus12);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus13",             cscRechitClusterNLayersChamberMinus13);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus21",             cscRechitClusterNLayersChamberMinus21);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus22",             cscRechitClusterNLayersChamberMinus22);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus31",             cscRechitClusterNLayersChamberMinus31);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus32",             cscRechitClusterNLayersChamberMinus32);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus41",             cscRechitClusterNLayersChamberMinus41);
  tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus42",             cscRechitClusterNLayersChamberMinus42);

  tree_->SetBranchAddress("cscRechitCluster_match_dtRechits_0p4",             cscRechitCluster_match_dtRechits_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_MB1_0p4",             cscRechitCluster_match_MB1_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_dtRechits_phi0p2",             cscRechitCluster_match_dtRechits_phi0p2);


  tree_->SetBranchAddress("cscRechitCluster_match_dtSeg_0p4",             cscRechitCluster_match_dtSeg_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_MB1Seg_0p4",             cscRechitCluster_match_MB1Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_RB1_0p4",             cscRechitCluster_match_RB1_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_RE12_0p4",             cscRechitCluster_match_RE12_0p4);

  tree_->SetBranchAddress("cscRechitCluster_match_cluster_dR",             cscRechitCluster_match_cluster_dR);
  tree_->SetBranchAddress("cscRechitCluster_match_cluster_index",             cscRechitCluster_match_cluster_index);

  tree_->SetBranchAddress("cscRechitClusterMet_dPhi",             cscRechitClusterMet_dPhi);
  tree_->SetBranchAddress("cscRechitClusterMetXYCorr_dPhi",             cscRechitClusterMetXYCorr_dPhi);


  tree_->SetBranchAddress("cscRechitClusterMetHEM_dPhi",             cscRechitClusterMetHEM_dPhi);
  tree_->SetBranchAddress("cscRechitClusterMetHEMXYCorr_dPhi",             cscRechitClusterMetHEMXYCorr_dPhi);
  tree_->SetBranchAddress("cscRechitClusterMetEENoise_dPhi",             cscRechitClusterMetEENoise_dPhi);
  tree_->SetBranchAddress("cscRechitClusterMetEENoiseXYCorr_dPhi",             cscRechitClusterMetEENoiseXYCorr_dPhi);
  tree_->SetBranchAddress("cscRechitClusterMetJesUp_dPhi",             cscRechitClusterMetJesUp_dPhi);
  tree_->SetBranchAddress("cscRechitClusterMetJesDown_dPhi",             cscRechitClusterMetJesDown_dPhi);


  tree_->SetBranchAddress("gLLP_multiplicity",    gLLP_multiplicity);
  tree_->SetBranchAddress("gLLP_multiplicity20",    gLLP_multiplicity20);
  tree_->SetBranchAddress("gLLP_EM_multiplicity",    gLLP_EM_multiplicity);

  tree_->SetBranchAddress("gLLP_eta",    gLLP_eta);
  tree_->SetBranchAddress("gLLP_phi",    gLLP_phi);
  tree_->SetBranchAddress("gLLP_csc",    gLLP_csc);
  tree_->SetBranchAddress("gLLP_dt",    gLLP_dt);
  tree_->SetBranchAddress("gLLP_ctau",    gLLP_ctau);
  tree_->SetBranchAddress("gLLP_beta",    gLLP_beta);
  tree_->SetBranchAddress("gLLP_maxMatchedDis",    gLLP_maxMatchedDis);


  tree_->SetBranchAddress("gLLP_e",    gLLP_e);
  tree_->SetBranchAddress("gLLP_pt",    gLLP_pt);
  tree_->SetBranchAddress("gLLP_lepdPhi",    gLLP_lepdPhi);
  tree_->SetBranchAddress("gLLP_EMFracE",    gLLP_EMFracE);
  tree_->SetBranchAddress("gLLP_EMFracEz",    gLLP_EMFracEz);
  tree_->SetBranchAddress("gLLP_EMFracP",    gLLP_EMFracP);
  tree_->SetBranchAddress("gLLP_EMFracPz",    gLLP_EMFracPz);
  tree_->SetBranchAddress("gLLP_daughterKaon",    gLLP_daughterKaon);

  tree_->SetBranchAddress("gLLP_decay_vertex_r",    gLLP_decay_vertex_r);
  tree_->SetBranchAddress("gLLP_decay_vertex_x",    gLLP_decay_vertex_x);
  tree_->SetBranchAddress("gLLP_decay_vertex_y",    gLLP_decay_vertex_y);
  tree_->SetBranchAddress("gLLP_decay_vertex_z",    gLLP_decay_vertex_z);

  tree_->SetBranchAddress("gLLP_deltaR",    &gLLP_deltaR);
  tree_->SetBranchAddress("gLLP_daughter_deltaR",    gLLP_daughter_deltaR);



tree_->SetBranchAddress("gLLP_daughter_id",          gLLP_daughter_id);
tree_->SetBranchAddress("gLLP_daughter_pt",          gLLP_daughter_pt);
    tree_->SetBranchAddress("gLLP_daughter_eta",          gLLP_daughter_eta);
    tree_->SetBranchAddress("gLLP_daughter_phi",          gLLP_daughter_phi);
    tree_->SetBranchAddress("gLLP_daughter_e",          gLLP_daughter_e);
    tree_->SetBranchAddress("gLLP_daughter_mass",          gLLP_daughter_mass);

  //Leptons
  tree_->SetBranchAddress("nLeptons",    &nLeptons);
  tree_->SetBranchAddress("lepE",        lepE);
  tree_->SetBranchAddress("lepPt",       lepPt);
  tree_->SetBranchAddress("lepEta",      lepEta);
  tree_->SetBranchAddress("lepPhi",      lepPhi);
  tree_->SetBranchAddress("lepPdgId",  lepPdgId);
  tree_->SetBranchAddress("lepDZ",     lepDZ);
  tree_->SetBranchAddress("lepEff", lepEff);
  tree_->SetBranchAddress("lepSF", lepSF);

  tree_->SetBranchAddress("lepTriggerSF", lepTriggerSF);
  tree_->SetBranchAddress("lepTightIdSF", lepTightIdSF);
  tree_->SetBranchAddress("lepLooseIdSF", lepLooseIdSF);
  tree_->SetBranchAddress("lepTightIsoSF", lepTightIsoSF);
  tree_->SetBranchAddress("lepLooseIsoSF", lepLooseIsoSF);
  // tree_->SetBranchAddress("lepTriggerMCEfficiency", lepTriggerMCEfficiency);
  // tree_->SetBranchAddress("lepTightIdMCEfficiency", lepTightIdMCEfficiency);
  // tree_->SetBranchAddress("lepLooseIdMCEfficiency", lepLooseIdMCEfficiency);
  // tree_->SetBranchAddress("lepTightIsoMCEfficiency", lepTightIsoMCEfficiency);
  // tree_->SetBranchAddress("lepLooseIsoMCEfficiency", lepLooseIsoMCEfficiency);
  tree_->SetBranchAddress("lepTag", lepTag);



  // tree_->SetBranchAddress("lepLoosePassId", lepLoosePassId);
  // tree_->SetBranchAddress("lepMediumPassId", lepMediumPassId);
  // tree_->SetBranchAddress("lepTightPassId", lepTightPassId);
  tree_->SetBranchAddress("lepPassId", lepPassId);
  tree_->SetBranchAddress("lepFromZ", lepFromZ);

  tree_->SetBranchAddress("lepPassVetoId", lepPassVetoId);
  tree_->SetBranchAddress("lepPassLooseIso", lepPassLooseIso);
  tree_->SetBranchAddress("lepPassTightIso", lepPassTightIso);
  tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);
  tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);


  //Z-candidate

  tree_->SetBranchAddress("ZMass1",       &ZMass1);

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

  tree_->SetBranchAddress("jetPtJESUp", jetPtJESUp);
  tree_->SetBranchAddress("jetPtJESDown", jetPtJESDown);
  tree_->SetBranchAddress("jetEJESUp", jetEJESUp);
  tree_->SetBranchAddress("jetEJESDown", jetEJESDown);
  tree_->SetBranchAddress("JecUnc", JecUnc);
  tree_->SetBranchAddress("gLLP_match_jet_pt", gLLP_match_jet_pt);
  tree_->SetBranchAddress("gLLP_match_jet_index", gLLP_match_jet_index);
  tree_->SetBranchAddress("gLLP_match_jet_minDeltaR", gLLP_match_jet_minDeltaR);
  tree_->SetBranchAddress("jet_match_llp_pt", jet_match_llp_pt);
  tree_->SetBranchAddress("jet_match_llp_index", jet_match_llp_index);
  tree_->SetBranchAddress("jet_match_llp_minDeltaR", jet_match_llp_minDeltaR);

  tree_->SetBranchAddress("jet_match_genJet_pt", jet_match_genJet_pt);
  tree_->SetBranchAddress("jet_match_genJet_index", jet_match_genJet_index);
  tree_->SetBranchAddress("jet_match_genJet_minDeltaR", jet_match_genJet_minDeltaR);

  tree_->SetBranchAddress("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction);
  tree_->SetBranchAddress("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction);
  tree_->SetBranchAddress("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction);
  tree_->SetBranchAddress("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction);
  // tree_->SetBranchAddress("jetMuonEnergyFraction", jetMuonEnergyFraction);

  // tree_->SetBranchAddress("jetLoosePassId", jetLoosePassId);
  tree_->SetBranchAddress("jetTightPassId", jetTightPassId);
  // triggers
  tree_->SetBranchAddress("HLTDecision",   HLTDecision);
  tree_->SetBranchAddress("METTrigger",   &METTrigger);
  tree_->SetBranchAddress("METNoMuTrigger",   &METNoMuTrigger);


};

void TreeMuonSystemCombination::LoadTree(const char* file)
{
  f_ = TFile::Open(file);
  assert(f_);
  tree_ = dynamic_cast<TTree*>(f_->Get("MuonSystem"));
  InitTree();
  assert(tree_);
};

void TreeMuonSystemCombination::CreateTree()
{
  tree_ = new TTree("MuonSystem","MuonSystem");
  f_ = 0;

  tree_->Branch("runNum",      &runNum,     "runNum/i");      // event run number
  tree_->Branch("MC_condition",      &MC_condition,     "MC_condition/i");      // event run number
  tree_->Branch("lumiSec",     &lumiSec,    "lumiSec/i");     // event lumi section
  tree_->Branch("evtNum",      &evtNum,     "evtNum/i");      // event number
  tree_->Branch("mH",      &mH,     "mH/I");      // event number
  tree_->Branch("mX",      &mX,     "mX/I");      // event number
  tree_->Branch("ctau",      &ctau,     "ctau/I");      // event number
  tree_->Branch("ZCategory",    &ZCategory,   "ZCategory/i");    // dilepton category

  tree_->Branch("category",    &category,   "category/i");    // dilepton category
  tree_->Branch("npv",         &npv,        "npv/i");         // number of primary vertices
  tree_->Branch("npu",         &npu,        "npu/i");         // number of in-time PU events (MC)
  tree_->Branch("weight",      &weight,     "weight/F");
  tree_->Branch("higgsPtWeight",      &higgsPtWeight,     "higgsPtWeight/F");
  tree_->Branch("higgsPtWeightSys",      higgsPtWeightSys,     "higgsPtWeightSys[9]/F");
  tree_->Branch("scaleWeights",      scaleWeights,     "scaleWeights[9]/F");
  tree_->Branch("lepOverallSF",      &lepOverallSF,     "lepOverallSF/F");


  tree_->Branch("sf_facScaleUp",      &sf_facScaleUp,     "sf_facScaleUp/F");
  tree_->Branch("sf_facScaleDown",      &sf_facScaleDown,     "sf_facScaleDown/F");
  tree_->Branch("sf_renScaleUp",      &sf_renScaleUp,     "sf_renScaleUp/F");
  tree_->Branch("sf_renScaleDown",      &sf_renScaleDown,     "sf_renScaleDown/F");
  tree_->Branch("sf_facRenScaleUp",      &sf_facRenScaleUp,     "sf_facRenScaleUp/F");
  tree_->Branch("sf_facRenScaleDown",      &sf_facRenScaleDown,     "sf_facRenScaleDown/F");
  tree_->Branch("metSF",      &metSF,     "metSF/F");

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
  tree_->Branch("EE_prefiring",      &EE_prefiring,     "EE_prefiring/O");



  tree_->Branch("rho",         &rho,        "rho/F");
  tree_->Branch("met",         &met,        "met/F");         // MET
  tree_->Branch("metNoMu",         &metNoMu,        "metNoMu/F");         // MET


  tree_->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)
  tree_->Branch("metXYCorr",      &metXYCorr,     "metXYCorr/F");      // phi(MET)
  tree_->Branch("metPhiXYCorr",      &metPhiXYCorr,     "metPhiXYCorr/F");      // phi(MET)

  tree_->Branch("HT",      &HT,     "HT/F");      // phi(MET)

  tree_->Branch("jetMet_dPhi",      &jetMet_dPhi,     "jetMet_dPhi/F");      // phi(MET)
  tree_->Branch("jetMet_dPhiMin",      &jetMet_dPhiMin,     "jetMet_dPhiMin/F");      // phi(MET)
  tree_->Branch("jetMet_dPhiMin4",      &jetMet_dPhiMin4,     "jetMet_dPhiMin4/F");      // phi(MET)

  tree_->Branch("metJESUp",      &metJESUp,     "metJESUp/F");      // phi(MET)
  tree_->Branch("metJESDown",      &metJESDown,     "metJESDown/F");      // phi(MET)
  tree_->Branch("metPhiJESUp",      &metPhiJESUp,     "metPhiJESUp/F");      // phi(metPhi)
  tree_->Branch("metPhiJESDown",      &metPhiJESDown,     "metPhiJESDown/F");      // phi(metPhi)
  tree_->Branch("metJESUpSF",      &metJESUpSF,     "metJESUpSF/F");      // phi(MET)
  tree_->Branch("metJESDownSF",      &metJESDownSF,     "metJESDownSF/F");      // phi(MET)
  tree_->Branch("metEENoise",      &metEENoise,     "metEENoise/F");      // phi(MET)
  tree_->Branch("metPhiEENoise",      &metPhiEENoise,     "metPhiEENoise/F");      // phi(MET)
  tree_->Branch("metHEM",      &metHEM,     "metHEM/F");      // phi(MET)
  tree_->Branch("metPhiHEM",      &metPhiHEM,     "metPhiHEM/F");      // phi(MET)
  tree_->Branch("metEENoiseXYCorr",      &metEENoiseXYCorr,     "metEENoiseXYCorr/F");      // phi(MET)
  tree_->Branch("metPhiEENoiseXYCorr",      &metPhiEENoiseXYCorr,     "metPhiEENoiseXYCorr/F");      // phi(MET)
  tree_->Branch("metHEMXYCorr",      &metHEMXYCorr,     "metHEMXYCorr/F");      // phi(MET)
  tree_->Branch("metPhiHEMXYCorr",      &metPhiHEMXYCorr,     "metPhiHEMXYCorr/F");      // phi(MET)
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

  tree_->Branch("gWPt",         &gWPt,        "gWPt/F");

  tree_->Branch("gLepId",      &gLepId,     "gLepId/I");      // phi(MET)
  tree_->Branch("gLepPt",      &gLepPt,     "gLepPt/F");      // phi(MET)
  tree_->Branch("gLepE",      &gLepE,     "gLepE/F");      // phi(MET)
  tree_->Branch("gLepEta",      &gLepEta,     "gLepEta/F");      // phi(MET)
  tree_->Branch("gLepPhi",      &gLepPhi,     "gLepPhi/F");      // phi(MET)
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


  tree_->Branch("nRpc",            &nRpc,             "nRpc/I");
  tree_->Branch("nDtSeg",            &nDtSeg,             "nDtSeg/I");

  tree_->Branch("nDTRechits",            &nDTRechits,             "nDTRechits/I");
  tree_->Branch("nDtRings",             &nDtRings, "nDtRings/I");
  tree_->Branch("nDtWheels25",             &nDtWheels25, "nDtWheels25/I");
  tree_->Branch("nDtStations25",             &nDtStations25, "nDtStations25/I");
  tree_->Branch("nDTPositiveYRechits",             &nDTPositiveYRechits, "nDTPositiveYRechits/I");
  tree_->Branch("nDTNegativeYRechits",             &nDTNegativeYRechits, "nDTNegativeYRechits/I");

  tree_->Branch("nDTRechitsWheelMinus2",             &nDTRechitsWheelMinus2, "nDTRechitsWheelMinus2/I");
  tree_->Branch("nDTRechitsWheelMinus1",             &nDTRechitsWheelMinus1, "nDTRechitsWheelMinus1/I");
  tree_->Branch("nDTRechitsWheel0",             &nDTRechitsWheel0, "nDTRechitsWheel0/I");
  tree_->Branch("nDTRechitsWheelPlus1",             &nDTRechitsWheelPlus1, "nDTRechitsWheelPlus1/I");
  tree_->Branch("nDTRechitsWheelPlus2",             &nDTRechitsWheelPlus2, "nDTRechitsWheelPlus2/I");

  tree_->Branch("nDTRechitsStation1",             &nDTRechitsStation1, "nDTRechitsStation1/I");
  tree_->Branch("nDTRechitsStation2",             &nDTRechitsStation2, "nDTRechitsStation2/I");
  tree_->Branch("nDTRechitsStation3",             &nDTRechitsStation3, "nDTRechitsStation3/I");
  tree_->Branch("nDTRechitsStation4",             &nDTRechitsStation4, "nDTRechitsStation4/I");


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
  // tree_->Branch("cscRechitsQuality",             cscRechitsQuality,             "cscRechitsQuality[nCscRechits]/I");

  // tree_->Branch("cscRechitsX",             cscRechitsX,             "cscRechitsX[nCscRechits]/F");

  // tree_->Branch("cscRechitsY",             cscRechitsY,             "cscRechitsY[nCscRechits]/F");
  // tree_->Branch("cscRechitsZ",             cscRechitsZ,             "cscRechitsZ[nCscRechits]/F");
  // tree_->Branch("cscRechitsTwire",             cscRechitsTwire,             "cscRechitsTwire[nCscRechits]/F");
  // tree_->Branch("cscRechitsTpeak",             cscRechitsTpeak,             "cscRechitsTpeak[nCscRechits]/F");
  // tree_->Branch("cscRechitsClusterId",             cscRechitsClusterId,             "cscRechitsClusterId[nCscRechits]/I");
  // tree_->Branch("cscRechitsCluster2Id",             cscRechitsCluster2Id,             "cscRechitsCluster2Id[nCscRechits]/I");
  tree_->Branch("dtRechitsX",             dtRechitsX,             "dtRechitsX[nDTRechits]/F");
  tree_->Branch("dtRechitsY",             dtRechitsY,             "dtRechitsY[nDTRechits]/F");
  tree_->Branch("dtRechitsZ",             dtRechitsZ,             "dtRechitsZ[nDTRechits]/F");

  tree_->Branch("dtRechitsEta",             dtRechitsEta,             "dtRechitsEta[nDTRechits]/F");
  tree_->Branch("dtRechitsPhi",             dtRechitsPhi,             "dtRechitsPhi[nDTRechits]/F");
  tree_->Branch("dtRechitsStation",             dtRechitsStation,             "dtRechitsStation[nDTRechits]/I");
  tree_->Branch("dtRechitsWheel",             dtRechitsWheel,             "dtRechitsWheel[nDTRechits]/I");
  tree_->Branch("dtRechitsClusterId",             dtRechitsClusterId,             "dtRechitsClusterId[nDTRechits]/I");



  tree_->Branch("rpcEta",             rpcEta,             "rpcEta[nRpc]/F");
  tree_->Branch("rpcPhi",             rpcPhi,             "rpcPhi[nRpc]/F");
  tree_->Branch("rpc_RE12",             rpc_RE12,             "rpc_RE12[nRpc]/O");
  tree_->Branch("rpc_RB1",             rpc_RB1,             "rpc_RB1[nRpc]/O");




    tree_->Branch("dtSegEta",             dtSegEta,             "dtSegEta[nDtSeg]/F");
    tree_->Branch("dtSegPhi",             dtSegPhi,             "dtSegPhi[nDtSeg]/F");
    tree_->Branch("dtSegWheel",             dtSegWheel,             "dtSegWheel[nDtSeg]/I");
    tree_->Branch("dtSegStation",             dtSegStation,             "dtSegStation[nDtSeg]/I");

    tree_->Branch("nCscRechitClusters",             &nCscRechitClusters, "nCscRechitClusters/I");

    // tree_->Branch("cscRechitCluster_match_gParticle_id",             cscRechitCluster_match_gParticle_id,             "cscRechitCluster_match_gParticle_id[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitCluster_match_gParticle",             cscRechitCluster_match_gParticle,             "cscRechitCluster_match_gParticle[nCscRechitClusters]/O");
    // tree_->Branch("cscRechitCluster_match_gParticle_minDeltaR",             cscRechitCluster_match_gParticle_minDeltaR,             "cscRechitCluster_match_gParticle_minDeltaR[nCscRechitClusters]/F");
    // tree_->Branch("cscRechitCluster_match_gParticle_index",             cscRechitCluster_match_gParticle_index,             "cscRechitCluster_match_gParticle_index[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitCluster_match_gParticle_eta",             cscRechitCluster_match_gParticle_eta,             "cscRechitCluster_match_gParticle_eta[nCscRechitClusters]/F");
    // tree_->Branch("cscRechitCluster_match_gParticle_phi",             cscRechitCluster_match_gParticle_phi,             "cscRechitCluster_match_gParticle_phi[nCscRechitClusters]/F");
    // tree_->Branch("cscRechitCluster_match_gParticle_E",             cscRechitCluster_match_gParticle_E,             "cscRechitCluster_match_gParticle_E[nCscRechitClusters]/F");
    // tree_->Branch("cscRechitCluster_match_gParticle_pt",             cscRechitCluster_match_gParticle_pt,             "cscRechitCluster_match_gParticle_pt[nCscRechitClusters]/F");
    // tree_->Branch("cscRechitCluster_match_gParticle_MotherId",             cscRechitCluster_match_gParticle_MotherId,             "cscRechitCluster_match_gParticle_MotherId[nCscRechitClusters]/I");

    tree_->Branch("cscRechitCluster_match_gLLP",             cscRechitCluster_match_gLLP,             "cscRechitCluster_match_gLLP[nCscRechitClusters]/O");
    tree_->Branch("cscRechitCluster_match_gLLP_minDeltaR",             cscRechitCluster_match_gLLP_minDeltaR,             "cscRechitCluster_match_gLLP_minDeltaR[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_index",             cscRechitCluster_match_gLLP_index,             "cscRechitCluster_match_gLLP_index[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_gLLP_eta",             cscRechitCluster_match_gLLP_eta, "cscRechitCluster_match_gLLP_eta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_phi",             cscRechitCluster_match_gLLP_phi, "cscRechitCluster_match_gLLP_phi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_decay_r",             cscRechitCluster_match_gLLP_decay_r, "cscRechitCluster_match_gLLP_decay_r[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_decay_x",             cscRechitCluster_match_gLLP_decay_x, "cscRechitCluster_match_gLLP_decay_x[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_decay_y",             cscRechitCluster_match_gLLP_decay_y, "cscRechitCluster_match_gLLP_decay_y[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_decay_z",             cscRechitCluster_match_gLLP_decay_z, "cscRechitCluster_match_gLLP_decay_z[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_ctau",             cscRechitCluster_match_gLLP_ctau, "cscRechitCluster_match_gLLP_ctau[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_beta",             cscRechitCluster_match_gLLP_beta, "cscRechitCluster_match_gLLP_beta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_csc",             cscRechitCluster_match_gLLP_csc, "cscRechitCluster_match_gLLP_csc[nCscRechitClusters]/O");
    tree_->Branch("cscRechitCluster_match_gLLP_dt",             cscRechitCluster_match_gLLP_dt, "cscRechitCluster_match_gLLP_dt[nCscRechitClusters]/O");
    tree_->Branch("cscRechitCluster_match_gLLP_multiplicity",             cscRechitCluster_match_gLLP_multiplicity, "cscRechitCluster_match_gLLP_multiplicity[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_gLLP_EM_multiplicity",             cscRechitCluster_match_gLLP_EM_multiplicity, "cscRechitCluster_match_gLLP_EM_multiplicity[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_gLLP_daughterKaon",             cscRechitCluster_match_gLLP_daughterKaon, "cscRechitCluster_match_gLLP_daughterKaon[nCscRechitClusters]/O");
    tree_->Branch("cscRechitCluster_match_gLLP_e",             cscRechitCluster_match_gLLP_e, "cscRechitCluster_match_gLLP_e[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_pt",             cscRechitCluster_match_gLLP_pt, "cscRechitCluster_match_gLLP_pt[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_EMFracE",             cscRechitCluster_match_gLLP_EMFracE, "cscRechitCluster_match_gLLP_EMFracE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_EMFracEz",             cscRechitCluster_match_gLLP_EMFracEz, "cscRechitCluster_match_gLLP_EMFracEz[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_EMFracP",             cscRechitCluster_match_gLLP_EMFracP, "cscRechitCluster_match_gLLP_EMFracP[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_EMFracPz",             cscRechitCluster_match_gLLP_EMFracPz, "cscRechitCluster_match_gLLP_EMFracPz[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_visE",             cscRechitCluster_match_gLLP_visE, "cscRechitCluster_match_gLLP_visE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_visEz",             cscRechitCluster_match_gLLP_visEz, "cscRechitCluster_match_gLLP_visEz[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_visP",             cscRechitCluster_match_gLLP_visP, "cscRechitCluster_match_gLLP_visP[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_visPz",             cscRechitCluster_match_gLLP_visPz, "cscRechitCluster_match_gLLP_visPz[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_lepdPhi",             cscRechitCluster_match_gLLP_lepdPhi, "cscRechitCluster_match_gLLP_lepdPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_daughter0_deltaR",             cscRechitCluster_match_gLLP_daughter0_deltaR, "cscRechitCluster_match_gLLP_daughter0_deltaR[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_daughter1_deltaR",             cscRechitCluster_match_gLLP_daughter1_deltaR, "cscRechitCluster_match_gLLP_daughter1_deltaR[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_daughter2_deltaR",             cscRechitCluster_match_gLLP_daughter2_deltaR, "cscRechitCluster_match_gLLP_daughter2_deltaR[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_daughter3_deltaR",             cscRechitCluster_match_gLLP_daughter3_deltaR, "cscRechitCluster_match_gLLP_daughter3_deltaR[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_other_daughter_deltaR",             cscRechitCluster_match_gLLP_other_daughter_deltaR, "cscRechitCluster_match_gLLP_other_daughter_deltaR[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_other_daughter_index",             cscRechitCluster_match_gLLP_other_daughter_index, "cscRechitCluster_match_gLLP_other_daughter_index[nCscRechitClusters]/I");


    tree_->Branch("cscRechitCluster_match_gLLP_other_eta",             cscRechitCluster_match_gLLP_other_eta, "cscRechitCluster_match_gLLP_other_eta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_other_phi",             cscRechitCluster_match_gLLP_other_phi, "cscRechitCluster_match_gLLP_other_phi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_other_decay_r",             cscRechitCluster_match_gLLP_other_decay_r, "cscRechitCluster_match_gLLP_other_decay_r[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_other_decay_x",             cscRechitCluster_match_gLLP_other_decay_x, "cscRechitCluster_match_gLLP_other_decay_x[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_other_decay_y",             cscRechitCluster_match_gLLP_other_decay_y, "cscRechitCluster_match_gLLP_other_decay_y[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_other_decay_z",             cscRechitCluster_match_gLLP_other_decay_z, "cscRechitCluster_match_gLLP_other_decay_z[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_other_ctau",             cscRechitCluster_match_gLLP_other_ctau, "cscRechitCluster_match_gLLP_other_ctau[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_other_beta",             cscRechitCluster_match_gLLP_other_beta, "cscRechitCluster_match_gLLP_other_beta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_other_csc",             cscRechitCluster_match_gLLP_other_csc, "cscRechitCluster_match_gLLP_other_csc[nCscRechitClusters]/O");
    tree_->Branch("cscRechitCluster_match_gLLP_other_e",             cscRechitCluster_match_gLLP_other_e, "cscRechitCluster_match_gLLP_other_e[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_gLLP_other_pt",             cscRechitCluster_match_gLLP_other_pt, "cscRechitCluster_match_gLLP_other_pt[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterX",             cscRechitClusterX,             "cscRechitClusterX[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterY",             cscRechitClusterY,             "cscRechitClusterY[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterZ",             cscRechitClusterZ,             "cscRechitClusterZ[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTime",             cscRechitClusterTime,             "cscRechitClusterTime[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeWire",             cscRechitClusterTimeWire,             "cscRechitClusterTimeWire[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeWirePruned",             cscRechitClusterTimeWirePruned,             "cscRechitClusterTimeWirePruned[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeTotal",             cscRechitClusterTimeTotal,             "cscRechitClusterTimeTotal[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterTimeSpread",             cscRechitClusterTimeSpread,             "cscRechitClusterTimeSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeTotalSpread",             cscRechitClusterTimeTotalSpread,             "cscRechitClusterTimeTotalSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeTotalSpreadPruned",             cscRechitClusterTimeTotalSpreadPruned,             "cscRechitClusterTimeTotalSpreadPruned[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeWireSpread",             cscRechitClusterTimeWireSpread,             "cscRechitClusterTimeWireSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterGenMuonDeltaR",             cscRechitClusterGenMuonDeltaR,             "cscRechitClusterGenMuonDeltaR[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterXYSpread",             cscRechitClusterXYSpread,             "cscRechitClusterXYSpread[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterMajorAxis",             cscRechitClusterMajorAxis,             "cscRechitClusterMajorAxis[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMinorAxis",             cscRechitClusterMinorAxis,             "cscRechitClusterMinorAxis[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEtaPhiSpread",             cscRechitClusterEtaPhiSpread,             "cscRechitClusterEtaPhiSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterPhiSpread",             cscRechitClusterPhiSpread,             "cscRechitClusterPhiSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEtaSpread",             cscRechitClusterEtaSpread,             "cscRechitClusterEtaSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterDeltaRSpread",             cscRechitClusterDeltaRSpread,             "cscRechitClusterDeltaRSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterXSpread",             cscRechitClusterXSpread,             "cscRechitClusterXSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterRSpread",             cscRechitClusterRSpread,             "cscRechitClusterRSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterYSpread",             cscRechitClusterYSpread,             "cscRechitClusterYSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterZSpread",             cscRechitClusterZSpread,             "cscRechitClusterZSpread[nCscRechitClusters]/F");


    tree_->Branch("cscRechitClusterPhi",             cscRechitClusterPhi,             "cscRechitClusterPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterEta",             cscRechitClusterEta,             "cscRechitClusterEta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoPt",             cscRechitClusterJetVetoPt,             "cscRechitClusterJetVetoPt[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoEta",             cscRechitClusterJetVetoEta,             "cscRechitClusterJetVetoEta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoPhi",             cscRechitClusterJetVetoPhi,             "cscRechitClusterJetVetoPhi[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterJetVetoE",             cscRechitClusterJetVetoE,             "cscRechitClusterJetVetoE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterGenJetVetoPt",             cscRechitClusterGenJetVetoPt,             "cscRechitClusterGenJetVetoPt[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterGenJetVetoE",             cscRechitClusterGenJetVetoE,             "cscRechitClusterGenJetVetoE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoPt",             cscRechitClusterMuonVetoPt,             "cscRechitClusterMuonVetoPt[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoE",             cscRechitClusterMuonVetoE,             "cscRechitClusterMuonVetoE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoPhi",             cscRechitClusterMuonVetoPhi,             "cscRechitClusterMuonVetoPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoEta",             cscRechitClusterMuonVetoEta,             "cscRechitClusterMuonVetoEta[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterJetVetoElectronEnergyFraction",             cscRechitClusterJetVetoElectronEnergyFraction,             "cscRechitClusterJetVetoElectronEnergyFraction[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoPhotonEnergyFraction",             cscRechitClusterJetVetoPhotonEnergyFraction,             "cscRechitClusterJetVetoPhotonEnergyFraction[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoChargedHadronEnergyFraction",             cscRechitClusterJetVetoChargedHadronEnergyFraction,             "cscRechitClusterJetVetoChargedHadronEnergyFraction[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoNeutralHadronEnergyFraction",             cscRechitClusterJetVetoNeutralHadronEnergyFraction,             "cscRechitClusterJetVetoNeutralHadronEnergyFraction[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoMuonEnergyFraction",             cscRechitClusterJetVetoMuonEnergyFraction,             "cscRechitClusterJetVetoMuonEnergyFraction[nCscRechitClusters]/F");


    tree_->Branch("cscRechitClusterJetVetoPt_0p6",             cscRechitClusterJetVetoPt_0p6,        "cscRechitClusterJetVetoPt_0p6[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoPt_0p8",             cscRechitClusterJetVetoPt_0p8,        "cscRechitClusterJetVetoPt_0p8[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoE_0p6",             cscRechitClusterJetVetoE_0p6,        "cscRechitClusterJetVetoE_0p6[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoE_0p8",             cscRechitClusterJetVetoE_0p8,        "cscRechitClusterJetVetoE_0p8[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoPt_0p6",             cscRechitClusterMuonVetoPt_0p6,        "cscRechitClusterMuonVetoPt_0p6[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoPt_0p8",             cscRechitClusterMuonVetoPt_0p8,        "cscRechitClusterMuonVetoPt_0p8[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoE_0p6",             cscRechitClusterMuonVetoE_0p6,        "cscRechitClusterMuonVetoE_0p6[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoE_0p8",             cscRechitClusterMuonVetoE_0p8,        "cscRechitClusterMuonVetoE_0p8[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoLooseIso",             cscRechitClusterMuonVetoLooseIso,             "cscRechitClusterMuonVetoLooseIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterMuonVetoTightIso",             cscRechitClusterMuonVetoTightIso,             "cscRechitClusterMuonVetoTightIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterMuonVetoVTightIso",             cscRechitClusterMuonVetoVTightIso,             "cscRechitClusterMuonVetoVTightIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterMuonVetoVVTightIso",             cscRechitClusterMuonVetoVVTightIso,             "cscRechitClusterMuonVetoVVTightIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterMuonVetoTightId",             cscRechitClusterMuonVetoTightId,             "cscRechitClusterMuonVetoTightId[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterMuonVetoLooseId",             cscRechitClusterMuonVetoLooseId,             "cscRechitClusterMuonVetoLooseId[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterMuonVetoGlobal",             cscRechitClusterMuonVetoGlobal,             "cscRechitClusterMuonVetoGlobal[nCscRechitClusters]/O");


    tree_->Branch("cscRechitClusterMuonVetoIso",             cscRechitClusterMuonVetoIso,             "cscRechitClusterMuonVetoIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterIsoMuonVetoPt",             cscRechitClusterIsoMuonVetoPt,             "cscRechitClusterIsoMuonVetoPt[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterIsoMuonVetoE",             cscRechitClusterIsoMuonVetoE,             "cscRechitClusterIsoMuonVetoE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterIsoMuonVetoPhi",             cscRechitClusterIsoMuonVetoPhi,             "cscRechitClusterIsoMuonVetoPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterIsoMuonVetoEta",             cscRechitClusterIsoMuonVetoEta,             "cscRechitClusterIsoMuonVetoEta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterGenMuonVetoPt",             cscRechitClusterGenMuonVetoPt,             "cscRechitClusterGenMuonVetoPt[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoType",             cscRechitClusterMuonVetoType,             "cscRechitClusterMuonVetoType[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterGenMuonVetoE",             cscRechitClusterGenMuonVetoE,             "cscRechitClusterGenMuonVetoE[nCscRechitClusters]/F");
    // tree_->Branch("cscRechitClusterGenMuonVetoProdX",             cscRechitClusterGenMuonVetoProdX,             "cscRechitClusterGenMuonVetoProdX[nCscRechitClusters]/F");
    // tree_->Branch("cscRechitClusterGenMuonVetoProdY",             cscRechitClusterGenMuonVetoProdY,             "cscRechitClusterGenMuonVetoProdY[nCscRechitClusters]/F");
    // tree_->Branch("cscRechitClusterGenMuonVetoProdZ",             cscRechitClusterGenMuonVetoProdZ,             "cscRechitClusterGenMuonVetoProdZ[nCscRechitClusters]/F");
    // tree_->Branch("cscRechitClusterGenMuonVetoLLPDist",             cscRechitClusterGenMuonVetoLLPDist,             "cscRechitClusterGenMuonVetoLLPDist[nCscRechitClusters]/F");
    // tree_->Branch("cscRechitClusterGenMuonVetoLLPIndex",             cscRechitClusterGenMuonVetoLLPIndex,             "cscRechitClusterGenMuonVetoLLPIndex[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterZLep1",             cscRechitClusterZLep1,             "cscRechitClusterZLep1[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep2",             cscRechitClusterZLep2,             "cscRechitClusterZLep2[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep1Tag",             cscRechitClusterZLep1Tag,             "cscRechitClusterZLep1Tag[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep2Tag",             cscRechitClusterZLep2Tag,             "cscRechitClusterZLep2Tag[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep1Id",             cscRechitClusterZLep1Id,             "cscRechitClusterZLep1Id[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterZLep2Id",             cscRechitClusterZLep2Id,             "cscRechitClusterZLep2Id[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster2ZLep1LooseIso",             cscRechitClusterZLep1LooseIso,             "cscRechitClusterZLep1LooseIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep1TightIso",             cscRechitClusterZLep1TightIso,             "cscRechitClusterZLep1TightIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep1VTightIso",             cscRechitClusterZLep1VTightIso,             "cscRechitClusterZLep1VTightIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep1VVTightIso",             cscRechitClusterZLep1VVTightIso,             "cscRechitClusterZLep1VVTightIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep1TightId",             cscRechitClusterZLep1TightId,             "cscRechitClusterZLep1TightId[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep2LooseIso",             cscRechitClusterZLep2LooseIso,             "cscRechitClusterZLep2LooseIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep2TightIso",             cscRechitClusterZLep2TightIso,             "cscRechitClusterZLep2TightIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep2VTightIso",             cscRechitClusterZLep2VTightIso,             "cscRechitClusterZLep2VTightIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep2VVTightIso",             cscRechitClusterZLep2VVTightIso,             "cscRechitClusterZLep2VVTightIso[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterZLep2TightId",             cscRechitClusterZLep2TightId,             "cscRechitClusterZLep2TightId[nCscRechitClusters]/O");


    tree_->Branch("cscRechitClusterSize",             cscRechitClusterSize,             "cscRechitClusterSize[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterMe11Ratio",             cscRechitClusterMe11Ratio,             "cscRechitClusterMe11Ratio[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMe12Ratio",             cscRechitClusterMe12Ratio,             "cscRechitClusterMe12Ratio[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterNStation",             cscRechitClusterNStation,             "cscRechitClusterNStation[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNStation5",             cscRechitClusterNStation5,             "cscRechitClusterNStation5[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNStation10",             cscRechitClusterNStation10,             "cscRechitClusterNStation10[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNStation10perc",             cscRechitClusterNStation10perc,             "cscRechitClusterNStation10perc[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterAvgStation",             cscRechitClusterAvgStation,             "cscRechitClusterAvgStation[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterAvgStation5",             cscRechitClusterAvgStation5,             "cscRechitClusterAvgStation5[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterAvgStation10",             cscRechitClusterAvgStation10,             "cscRechitClusterAvgStation10[nCscRechitClusters]/F");
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
    tree_->Branch("cscRechitClusterMetXYCorr_dPhi",             cscRechitClusterMetXYCorr_dPhi,             "cscRechitClusterMetXYCorr_dPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_cscRechits_0p4",             cscRechitCluster_match_cscRechits_0p4,             "cscRechitCluster_match_cscRechits_0p4[nCscRechitClusters]/I");

    tree_->Branch("cscRechitClusterNLayersChamberPlus11",             cscRechitClusterNLayersChamberPlus11,             "cscRechitClusterNLayersChamberPlus11[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberPlus12",             cscRechitClusterNLayersChamberPlus12,             "cscRechitClusterNLayersChamberPlus12[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberPlus13",             cscRechitClusterNLayersChamberPlus13,             "cscRechitClusterNLayersChamberPlus13[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberPlus21",             cscRechitClusterNLayersChamberPlus21,             "cscRechitClusterNLayersChamberPlus21[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberPlus22",             cscRechitClusterNLayersChamberPlus22,             "cscRechitClusterNLayersChamberPlus22[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberPlus31",             cscRechitClusterNLayersChamberPlus31,             "cscRechitClusterNLayersChamberPlus31[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberPlus32",             cscRechitClusterNLayersChamberPlus32,             "cscRechitClusterNLayersChamberPlus32[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberPlus41",             cscRechitClusterNLayersChamberPlus41,             "cscRechitClusterNLayersChamberPlus41[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberPlus42",             cscRechitClusterNLayersChamberPlus42,             "cscRechitClusterNLayersChamberPlus42[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberMinus11",             cscRechitClusterNLayersChamberMinus11,             "cscRechitClusterNLayersChamberMinus11[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberMinus12",             cscRechitClusterNLayersChamberMinus12,             "cscRechitClusterNLayersChamberMinus12[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberMinus13",             cscRechitClusterNLayersChamberMinus13,             "cscRechitClusterNLayersChamberMinus13[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberMinus21",             cscRechitClusterNLayersChamberMinus21,             "cscRechitClusterNLayersChamberMinus21[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberMinus22",             cscRechitClusterNLayersChamberMinus22,             "cscRechitClusterNLayersChamberMinus22[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberMinus31",             cscRechitClusterNLayersChamberMinus31,             "cscRechitClusterNLayersChamberMinus31[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberMinus32",             cscRechitClusterNLayersChamberMinus32,             "cscRechitClusterNLayersChamberMinus32[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberMinus41",             cscRechitClusterNLayersChamberMinus41,             "cscRechitClusterNLayersChamberMinus41[nCscRechitClusters]/I");
    tree_->Branch("cscRechitClusterNLayersChamberMinus42",             cscRechitClusterNLayersChamberMinus42,             "cscRechitClusterNLayersChamberMinus42[nCscRechitClusters]/I");


    tree_->Branch("cscRechitClusterMetHEM_dPhi",             cscRechitClusterMetHEM_dPhi,             "cscRechitClusterMetHEM_dPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMetHEMXYCorr_dPhi",             cscRechitClusterMetHEMXYCorr_dPhi,             "cscRechitClusterMetHEMXYCorr_dPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMetEENoise_dPhi",             cscRechitClusterMetEENoise_dPhi,             "cscRechitClusterMetEENoise_dPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMetEENoiseXYCorr_dPhi",             cscRechitClusterMetEENoiseXYCorr_dPhi,             "cscRechitClusterMetEENoiseXYCorr_dPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMetJesUp_dPhi",             cscRechitClusterMetJesUp_dPhi,             "cscRechitClusterMetJesUp_dPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMetJesDown_dPhi",             cscRechitClusterMetJesDown_dPhi,             "cscRechitClusterMetJesDown_dPhi[nCscRechitClusters]/F");



    tree_->Branch("cscRechitCluster_match_dtRechits_phi0p2",             cscRechitCluster_match_dtRechits_phi0p2,             "cscRechitCluster_match_dtRechits_phi0p2[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_dtRechits_0p4",             cscRechitCluster_match_dtRechits_0p4,             "cscRechitCluster_match_dtRechits_0p4[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_MB1_0p4",             cscRechitCluster_match_MB1_0p4,             "cscRechitCluster_match_MB1_0p4[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_dtSeg_0p4",             cscRechitCluster_match_dtSeg_0p4,             "cscRechitCluster_match_dtSeg_0p4[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_MB1Seg_0p4",             cscRechitCluster_match_MB1Seg_0p4,             "cscRechitCluster_match_MB1Seg_0p4[nCscRechitClusters]/I");

    tree_->Branch("cscRechitCluster_match_RB1_0p4",             cscRechitCluster_match_RB1_0p4,             "cscRechitCluster_match_RB1_0p4[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_RE12_0p4",             cscRechitCluster_match_RE12_0p4,             "cscRechitCluster_match_RE12_0p4[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_cluster_dR",             cscRechitCluster_match_cluster_dR,             "cscRechitCluster_match_cluster_dR[nCscRechitClusters]/F");
    tree_->Branch("cscRechitCluster_match_cluster_index",             cscRechitCluster_match_cluster_index,             "cscRechitCluster_match_cluster_index[nCscRechitClusters]/I");




    tree_->Branch("nDtRechitClusters",             &nDtRechitClusters, "nDtRechitClusters/I");
    tree_->Branch("nDtRechitClusters2",             &nDtRechitClusters2, "nDtRechitClusters2/I");



        // tree_->Branch("dtRechitCluster_match_gParticle_id",             dtRechitCluster_match_gParticle_id,             "dtRechitCluster_match_gParticle_id[nDtRechitClusters]/I");
        // tree_->Branch("dtRechitCluster_match_gParticle",             dtRechitCluster_match_gParticle,             "dtRechitCluster_match_gParticle[nDtRechitClusters]/O");
        // tree_->Branch("dtRechitCluster_match_gParticle_minDeltaR",             dtRechitCluster_match_gParticle_minDeltaR,             "dtRechitCluster_match_gParticle_minDeltaR[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitCluster_match_gParticle_index",             dtRechitCluster_match_gParticle_index,             "dtRechitCluster_match_gParticle_index[nDtRechitClusters]/I");
        // tree_->Branch("dtRechitCluster_match_gParticle_eta",             dtRechitCluster_match_gParticle_eta,             "dtRechitCluster_match_gParticle_eta[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitCluster_match_gParticle_phi",             dtRechitCluster_match_gParticle_phi,             "dtRechitCluster_match_gParticle_phi[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitCluster_match_gParticle_E",             dtRechitCluster_match_gParticle_E,             "dtRechitCluster_match_gParticle_E[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitCluster_match_gParticle_pt",             dtRechitCluster_match_gParticle_pt,             "dtRechitCluster_match_gParticle_pt[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitCluster_match_gParticle_MotherId",             dtRechitCluster_match_gParticle_MotherId,             "dtRechitCluster_match_gParticle_MotherId[nDtRechitClusters]/I");

        tree_->Branch("dtRechitCluster_match_gLLP",             dtRechitCluster_match_gLLP,             "dtRechitCluster_match_gLLP[nDtRechitClusters]/O");
        tree_->Branch("dtRechitCluster_match_gLLP_minDeltaR",             dtRechitCluster_match_gLLP_minDeltaR,             "dtRechitCluster_match_gLLP_minDeltaR[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_index",             dtRechitCluster_match_gLLP_index,             "dtRechitCluster_match_gLLP_index[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_gLLP_eta",             dtRechitCluster_match_gLLP_eta, "dtRechitCluster_match_gLLP_eta[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_phi",             dtRechitCluster_match_gLLP_phi, "dtRechitCluster_match_gLLP_phi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_decay_r",             dtRechitCluster_match_gLLP_decay_r, "dtRechitCluster_match_gLLP_decay_r[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_decay_x",             dtRechitCluster_match_gLLP_decay_x, "dtRechitCluster_match_gLLP_decay_x[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_decay_y",             dtRechitCluster_match_gLLP_decay_y, "dtRechitCluster_match_gLLP_decay_y[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_decay_z",             dtRechitCluster_match_gLLP_decay_z, "dtRechitCluster_match_gLLP_decay_z[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_ctau",             dtRechitCluster_match_gLLP_ctau, "dtRechitCluster_match_gLLP_ctau[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_beta",             dtRechitCluster_match_gLLP_beta, "dtRechitCluster_match_gLLP_beta[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_csc",             dtRechitCluster_match_gLLP_csc, "dtRechitCluster_match_gLLP_csc[nDtRechitClusters]/O");
        tree_->Branch("dtRechitCluster_match_gLLP_dt",             dtRechitCluster_match_gLLP_dt, "dtRechitCluster_match_gLLP_dt[nDtRechitClusters]/O");
        tree_->Branch("dtRechitCluster_match_gLLP_multiplicity",             dtRechitCluster_match_gLLP_multiplicity, "dtRechitCluster_match_gLLP_multiplicity[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_gLLP_EM_multiplicity",             dtRechitCluster_match_gLLP_EM_multiplicity, "dtRechitCluster_match_gLLP_EM_multiplicity[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_gLLP_daughterKaon",             dtRechitCluster_match_gLLP_daughterKaon, "dtRechitCluster_match_gLLP_daughterKaon[nDtRechitClusters]/O");
        tree_->Branch("dtRechitCluster_match_gLLP_e",             dtRechitCluster_match_gLLP_e, "dtRechitCluster_match_gLLP_e[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_pt",             dtRechitCluster_match_gLLP_pt, "dtRechitCluster_match_gLLP_pt[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_EMFracE",             dtRechitCluster_match_gLLP_EMFracE, "dtRechitCluster_match_gLLP_EMFracE[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_EMFracEz",             dtRechitCluster_match_gLLP_EMFracEz, "dtRechitCluster_match_gLLP_EMFracEz[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_EMFracP",             dtRechitCluster_match_gLLP_EMFracP, "dtRechitCluster_match_gLLP_EMFracP[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_EMFracPz",             dtRechitCluster_match_gLLP_EMFracPz, "dtRechitCluster_match_gLLP_EMFracPz[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_visE",             dtRechitCluster_match_gLLP_visE, "dtRechitCluster_match_gLLP_visE[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_visEz",             dtRechitCluster_match_gLLP_visEz, "dtRechitCluster_match_gLLP_visEz[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_visP",             dtRechitCluster_match_gLLP_visP, "dtRechitCluster_match_gLLP_visP[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_visPz",             dtRechitCluster_match_gLLP_visPz, "dtRechitCluster_match_gLLP_visPz[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_lepdPhi",             dtRechitCluster_match_gLLP_lepdPhi, "dtRechitCluster_match_gLLP_lepdPhi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_daughter0_deltaR",             dtRechitCluster_match_gLLP_daughter0_deltaR, "dtRechitCluster_match_gLLP_daughter0_deltaR[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_daughter1_deltaR",             dtRechitCluster_match_gLLP_daughter1_deltaR, "dtRechitCluster_match_gLLP_daughter1_deltaR[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_daughter2_deltaR",             dtRechitCluster_match_gLLP_daughter2_deltaR, "dtRechitCluster_match_gLLP_daughter2_deltaR[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_daughter3_deltaR",             dtRechitCluster_match_gLLP_daughter3_deltaR, "dtRechitCluster_match_gLLP_daughter3_deltaR[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_other_daughter_deltaR",             dtRechitCluster_match_gLLP_other_daughter_deltaR, "dtRechitCluster_match_gLLP_other_daughter_deltaR[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_other_daughter_index",             dtRechitCluster_match_gLLP_other_daughter_index, "dtRechitCluster_match_gLLP_other_daughter_index[nDtRechitClusters]/I");


        tree_->Branch("dtRechitCluster_match_gLLP_other_eta",             dtRechitCluster_match_gLLP_other_eta, "dtRechitCluster_match_gLLP_other_eta[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_other_phi",             dtRechitCluster_match_gLLP_other_phi, "dtRechitCluster_match_gLLP_other_phi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_other_decay_r",             dtRechitCluster_match_gLLP_other_decay_r, "dtRechitCluster_match_gLLP_other_decay_r[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_other_decay_x",             dtRechitCluster_match_gLLP_other_decay_x, "dtRechitCluster_match_gLLP_other_decay_x[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_other_decay_y",             dtRechitCluster_match_gLLP_other_decay_y, "dtRechitCluster_match_gLLP_other_decay_y[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_other_decay_z",             dtRechitCluster_match_gLLP_other_decay_z, "dtRechitCluster_match_gLLP_other_decay_z[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_other_ctau",             dtRechitCluster_match_gLLP_other_ctau, "dtRechitCluster_match_gLLP_other_ctau[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_other_beta",             dtRechitCluster_match_gLLP_other_beta, "dtRechitCluster_match_gLLP_other_beta[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_other_csc",             dtRechitCluster_match_gLLP_other_csc, "dtRechitCluster_match_gLLP_other_csc[nDtRechitClusters]/O");
        tree_->Branch("dtRechitCluster_match_gLLP_other_e",             dtRechitCluster_match_gLLP_other_e, "dtRechitCluster_match_gLLP_other_e[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gLLP_other_pt",             dtRechitCluster_match_gLLP_other_pt, "dtRechitCluster_match_gLLP_other_pt[nDtRechitClusters]/F");

        tree_->Branch("dtRechitClusterX",             dtRechitClusterX,             "dtRechitClusterX[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterY",             dtRechitClusterY,             "dtRechitClusterY[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterZ",             dtRechitClusterZ,             "dtRechitClusterZ[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterTime",             dtRechitClusterTime,             "dtRechitClusterTime[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterTimeWire",             dtRechitClusterTimeWire,             "dtRechitClusterTimeWire[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterTimeWirePruned",             dtRechitClusterTimeWirePruned,             "dtRechitClusterTimeWirePruned[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterTimeTotal",             dtRechitClusterTimeTotal,             "dtRechitClusterTimeTotal[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterWheel",             dtRechitClusterWheel,             "dtRechitClusterWheel[nDtRechitClusters]/I");

        tree_->Branch("dtRechitClusterTimeSpread",             dtRechitClusterTimeSpread,             "dtRechitClusterTimeSpread[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterTimeTotalSpread",             dtRechitClusterTimeTotalSpread,             "dtRechitClusterTimeTotalSpread[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterTimeTotalSpreadPruned",             dtRechitClusterTimeTotalSpreadPruned,             "dtRechitClusterTimeTotalSpreadPruned[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterTimeWireSpread",             dtRechitClusterTimeWireSpread,             "dtRechitClusterTimeWireSpread[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterGenMuonDeltaR",             dtRechitClusterGenMuonDeltaR,             "dtRechitClusterGenMuonDeltaR[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterXYSpread",             dtRechitClusterXYSpread,             "dtRechitClusterXYSpread[nDtRechitClusters]/F");

        tree_->Branch("dtRechitClusterMajorAxis",             dtRechitClusterMajorAxis,             "dtRechitClusterMajorAxis[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMinorAxis",             dtRechitClusterMinorAxis,             "dtRechitClusterMinorAxis[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterEtaPhiSpread",             dtRechitClusterEtaPhiSpread,             "dtRechitClusterEtaPhiSpread[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterPhiSpread",             dtRechitClusterPhiSpread,             "dtRechitClusterPhiSpread[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterEtaSpread",             dtRechitClusterEtaSpread,             "dtRechitClusterEtaSpread[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterDeltaRSpread",             dtRechitClusterDeltaRSpread,             "dtRechitClusterDeltaRSpread[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterXSpread",             dtRechitClusterXSpread,             "dtRechitClusterXSpread[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterRSpread",             dtRechitClusterRSpread,             "dtRechitClusterRSpread[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterYSpread",             dtRechitClusterYSpread,             "dtRechitClusterYSpread[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterZSpread",             dtRechitClusterZSpread,             "dtRechitClusterZSpread[nDtRechitClusters]/F");


        tree_->Branch("dtRechitClusterPhi",             dtRechitClusterPhi,             "dtRechitClusterPhi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterEta",             dtRechitClusterEta,             "dtRechitClusterEta[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoPt",             dtRechitClusterJetVetoPt,             "dtRechitClusterJetVetoPt[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoEta",             dtRechitClusterJetVetoEta,             "dtRechitClusterJetVetoEta[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoPhi",             dtRechitClusterJetVetoPhi,             "dtRechitClusterJetVetoPhi[nDtRechitClusters]/F");

        tree_->Branch("dtRechitClusterJetVetoE",             dtRechitClusterJetVetoE,             "dtRechitClusterJetVetoE[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterGenJetVetoPt",             dtRechitClusterGenJetVetoPt,             "dtRechitClusterGenJetVetoPt[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterGenJetVetoE",             dtRechitClusterGenJetVetoE,             "dtRechitClusterGenJetVetoE[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoPt",             dtRechitClusterMuonVetoPt,             "dtRechitClusterMuonVetoPt[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoE",             dtRechitClusterMuonVetoE,             "dtRechitClusterMuonVetoE[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoPhi",             dtRechitClusterMuonVetoPhi,             "dtRechitClusterMuonVetoPhi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoEta",             dtRechitClusterMuonVetoEta,             "dtRechitClusterMuonVetoEta[nDtRechitClusters]/F");

        tree_->Branch("dtRechitClusterJetVetoElectronEnergyFraction",             dtRechitClusterJetVetoElectronEnergyFraction,             "dtRechitClusterJetVetoElectronEnergyFraction[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoPhotonEnergyFraction",             dtRechitClusterJetVetoPhotonEnergyFraction,             "dtRechitClusterJetVetoPhotonEnergyFraction[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoChargedHadronEnergyFraction",             dtRechitClusterJetVetoChargedHadronEnergyFraction,             "dtRechitClusterJetVetoChargedHadronEnergyFraction[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoNeutralHadronEnergyFraction",             dtRechitClusterJetVetoNeutralHadronEnergyFraction,             "dtRechitClusterJetVetoNeutralHadronEnergyFraction[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoMuonEnergyFraction",             dtRechitClusterJetVetoMuonEnergyFraction,             "dtRechitClusterJetVetoMuonEnergyFraction[nDtRechitClusters]/F");


        tree_->Branch("dtRechitClusterJetVetoPt_0p6",             dtRechitClusterJetVetoPt_0p6,        "dtRechitClusterJetVetoPt_0p6[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoPt_0p8",             dtRechitClusterJetVetoPt_0p8,        "dtRechitClusterJetVetoPt_0p8[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoE_0p6",             dtRechitClusterJetVetoE_0p6,        "dtRechitClusterJetVetoE_0p6[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoE_0p8",             dtRechitClusterJetVetoE_0p8,        "dtRechitClusterJetVetoE_0p8[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoPt_0p6",             dtRechitClusterMuonVetoPt_0p6,        "dtRechitClusterMuonVetoPt_0p6[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoPt_0p8",             dtRechitClusterMuonVetoPt_0p8,        "dtRechitClusterMuonVetoPt_0p8[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoE_0p6",             dtRechitClusterMuonVetoE_0p6,        "dtRechitClusterMuonVetoE_0p6[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoE_0p8",             dtRechitClusterMuonVetoE_0p8,        "dtRechitClusterMuonVetoE_0p8[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoLooseIso",             dtRechitClusterMuonVetoLooseIso,             "dtRechitClusterMuonVetoLooseIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterMuonVetoTightIso",             dtRechitClusterMuonVetoTightIso,             "dtRechitClusterMuonVetoTightIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterMuonVetoVTightIso",             dtRechitClusterMuonVetoVTightIso,             "dtRechitClusterMuonVetoVTightIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterMuonVetoVVTightIso",             dtRechitClusterMuonVetoVVTightIso,             "dtRechitClusterMuonVetoVVTightIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterMuonVetoTightId",             dtRechitClusterMuonVetoTightId,             "dtRechitClusterMuonVetoTightId[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterMuonVetoLooseId",             dtRechitClusterMuonVetoLooseId,             "dtRechitClusterMuonVetoLooseId[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterMuonVetoGlobal",             dtRechitClusterMuonVetoGlobal,             "dtRechitClusterMuonVetoGlobal[nDtRechitClusters]/O");


        tree_->Branch("dtRechitClusterMuonVetoIso",             dtRechitClusterMuonVetoIso,             "dtRechitClusterMuonVetoIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterIsoMuonVetoPt",             dtRechitClusterIsoMuonVetoPt,             "dtRechitClusterIsoMuonVetoPt[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterIsoMuonVetoE",             dtRechitClusterIsoMuonVetoE,             "dtRechitClusterIsoMuonVetoE[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterIsoMuonVetoPhi",             dtRechitClusterIsoMuonVetoPhi,             "dtRechitClusterIsoMuonVetoPhi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterIsoMuonVetoEta",             dtRechitClusterIsoMuonVetoEta,             "dtRechitClusterIsoMuonVetoEta[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterGenMuonVetoPt",             dtRechitClusterGenMuonVetoPt,             "dtRechitClusterGenMuonVetoPt[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoType",             dtRechitClusterMuonVetoType,             "dtRechitClusterMuonVetoType[nDtRechitClusters]/I");
        // tree_->Branch("dtRechitClusterGenMuonVetoE",             dtRechitClusterGenMuonVetoE,             "dtRechitClusterGenMuonVetoE[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitClusterGenMuonVetoProdX",             dtRechitClusterGenMuonVetoProdX,             "dtRechitClusterGenMuonVetoProdX[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitClusterGenMuonVetoProdY",             dtRechitClusterGenMuonVetoProdY,             "dtRechitClusterGenMuonVetoProdY[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitClusterGenMuonVetoProdZ",             dtRechitClusterGenMuonVetoProdZ,             "dtRechitClusterGenMuonVetoProdZ[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitClusterGenMuonVetoLLPDist",             dtRechitClusterGenMuonVetoLLPDist,             "dtRechitClusterGenMuonVetoLLPDist[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitClusterGenMuonVetoLLPIndex",             dtRechitClusterGenMuonVetoLLPIndex,             "dtRechitClusterGenMuonVetoLLPIndex[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterZLep1",             dtRechitClusterZLep1,             "dtRechitClusterZLep1[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep2",             dtRechitClusterZLep2,             "dtRechitClusterZLep2[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep1Tag",             dtRechitClusterZLep1Tag,             "dtRechitClusterZLep1Tag[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep2Tag",             dtRechitClusterZLep2Tag,             "dtRechitClusterZLep2Tag[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep1Id",             dtRechitClusterZLep1Id,             "dtRechitClusterZLep1Id[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterZLep2Id",             dtRechitClusterZLep2Id,             "dtRechitClusterZLep2Id[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster2ZLep1LooseIso",             dtRechitClusterZLep1LooseIso,             "dtRechitClusterZLep1LooseIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep1TightIso",             dtRechitClusterZLep1TightIso,             "dtRechitClusterZLep1TightIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep1VTightIso",             dtRechitClusterZLep1VTightIso,             "dtRechitClusterZLep1VTightIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep1VVTightIso",             dtRechitClusterZLep1VVTightIso,             "dtRechitClusterZLep1VVTightIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep1TightId",             dtRechitClusterZLep1TightId,             "dtRechitClusterZLep1TightId[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep2LooseIso",             dtRechitClusterZLep2LooseIso,             "dtRechitClusterZLep2LooseIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep2TightIso",             dtRechitClusterZLep2TightIso,             "dtRechitClusterZLep2TightIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep2VTightIso",             dtRechitClusterZLep2VTightIso,             "dtRechitClusterZLep2VTightIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep2VVTightIso",             dtRechitClusterZLep2VVTightIso,             "dtRechitClusterZLep2VVTightIso[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterZLep2TightId",             dtRechitClusterZLep2TightId,             "dtRechitClusterZLep2TightId[nDtRechitClusters]/O");


        tree_->Branch("dtRechitClusterSize",             dtRechitClusterSize,             "dtRechitClusterSize[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterMe11Ratio",             dtRechitClusterMe11Ratio,             "dtRechitClusterMe11Ratio[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMe12Ratio",             dtRechitClusterMe12Ratio,             "dtRechitClusterMe12Ratio[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterNStation",             dtRechitClusterNStation,             "dtRechitClusterNStation[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNStation5",             dtRechitClusterNStation5,             "dtRechitClusterNStation5[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNStation10",             dtRechitClusterNStation10,             "dtRechitClusterNStation10[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNStation10perc",             dtRechitClusterNStation10perc,             "dtRechitClusterNStation10perc[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterAvgStation",             dtRechitClusterAvgStation,             "dtRechitClusterAvgStation[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterAvgStation5",             dtRechitClusterAvgStation5,             "dtRechitClusterAvgStation5[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterAvgStation10",             dtRechitClusterAvgStation10,             "dtRechitClusterAvgStation10[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterAvgStation10perc",             dtRechitClusterAvgStation10perc,             "dtRechitClusterAvgStation10perc[nDtRechitClusters]/F");


        tree_->Branch("dtRechitClusterMaxStation",             dtRechitClusterMaxStation,             "dtRechitClusterMaxStation[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterMaxStationRatio",             dtRechitClusterMaxStationRatio,             "dtRechitClusterMaxStationRatio[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterNChamber",             dtRechitClusterNChamber,             "dtRechitClusterNChamber[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterMaxChamber",             dtRechitClusterMaxChamber,             "dtRechitClusterMaxChamber[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterMaxChamberRatio",             dtRechitClusterMaxChamberRatio,             "dtRechitClusterMaxChamberRatio[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterNSegmentStation1",             dtRechitClusterNSegmentStation1,             "dtRechitClusterNSegmentStation1[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNSegmentStation2",             dtRechitClusterNSegmentStation2,             "dtRechitClusterNSegmentStation2[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNSegmentStation3",             dtRechitClusterNSegmentStation3,             "dtRechitClusterNSegmentStation3[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNSegmentStation4",             dtRechitClusterNSegmentStation4,             "dtRechitClusterNSegmentStation4[nDtRechitClusters]/I");

        tree_->Branch("dtRechitClusterNRechitChamberPlus11",             dtRechitClusterNRechitChamberPlus11,             "dtRechitClusterNRechitChamberPlus11[nDtRechitClusters]/I");
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
        tree_->Branch("dtRechitClusterMet_dPhi",             dtRechitClusterMet_dPhi,             "dtRechitClusterMet_dPhi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMetXYCorr_dPhi",             dtRechitClusterMetXYCorr_dPhi,             "dtRechitClusterMetXYCorr_dPhi[nDtRechitClusters]/F");

        tree_->Branch("dtRechitClusterNLayersChamberPlus11",             dtRechitClusterNLayersChamberPlus11,             "dtRechitClusterNLayersChamberPlus11[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberPlus12",             dtRechitClusterNLayersChamberPlus12,             "dtRechitClusterNLayersChamberPlus12[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberPlus13",             dtRechitClusterNLayersChamberPlus13,             "dtRechitClusterNLayersChamberPlus13[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberPlus21",             dtRechitClusterNLayersChamberPlus21,             "dtRechitClusterNLayersChamberPlus21[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberPlus22",             dtRechitClusterNLayersChamberPlus22,             "dtRechitClusterNLayersChamberPlus22[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberPlus31",             dtRechitClusterNLayersChamberPlus31,             "dtRechitClusterNLayersChamberPlus31[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberPlus32",             dtRechitClusterNLayersChamberPlus32,             "dtRechitClusterNLayersChamberPlus32[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberPlus41",             dtRechitClusterNLayersChamberPlus41,             "dtRechitClusterNLayersChamberPlus41[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberPlus42",             dtRechitClusterNLayersChamberPlus42,             "dtRechitClusterNLayersChamberPlus42[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberMinus11",             dtRechitClusterNLayersChamberMinus11,             "dtRechitClusterNLayersChamberMinus11[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberMinus12",             dtRechitClusterNLayersChamberMinus12,             "dtRechitClusterNLayersChamberMinus12[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberMinus13",             dtRechitClusterNLayersChamberMinus13,             "dtRechitClusterNLayersChamberMinus13[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberMinus21",             dtRechitClusterNLayersChamberMinus21,             "dtRechitClusterNLayersChamberMinus21[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberMinus22",             dtRechitClusterNLayersChamberMinus22,             "dtRechitClusterNLayersChamberMinus22[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberMinus31",             dtRechitClusterNLayersChamberMinus31,             "dtRechitClusterNLayersChamberMinus31[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberMinus32",             dtRechitClusterNLayersChamberMinus32,             "dtRechitClusterNLayersChamberMinus32[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberMinus41",             dtRechitClusterNLayersChamberMinus41,             "dtRechitClusterNLayersChamberMinus41[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNLayersChamberMinus42",             dtRechitClusterNLayersChamberMinus42,             "dtRechitClusterNLayersChamberMinus42[nDtRechitClusters]/I");


        tree_->Branch("dtRechitClusterMetHEM_dPhi",             dtRechitClusterMetHEM_dPhi,             "dtRechitClusterMetHEM_dPhi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMetHEMXYCorr_dPhi",             dtRechitClusterMetHEMXYCorr_dPhi,             "dtRechitClusterMetHEMXYCorr_dPhi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMetEENoise_dPhi",             dtRechitClusterMetEENoise_dPhi,             "dtRechitClusterMetEENoise_dPhi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMetEENoiseXYCorr_dPhi",             dtRechitClusterMetEENoiseXYCorr_dPhi,             "dtRechitClusterMetEENoiseXYCorr_dPhi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMetJesUp_dPhi",             dtRechitClusterMetJesUp_dPhi,             "dtRechitClusterMetJesUp_dPhi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMetJesDown_dPhi",             dtRechitClusterMetJesDown_dPhi,             "dtRechitClusterMetJesDown_dPhi[nDtRechitClusters]/F");


        tree_->Branch("dtRechitCluster_match_dtSeg_0p5",             dtRechitCluster_match_dtSeg_0p5,             "dtRechitCluster_match_dtSeg_0p5[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_dtSegTime_0p5",             dtRechitCluster_match_dtSegTime_0p5,             "dtRechitCluster_match_dtSegTime_0p5[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_dtSeg_0p4",             dtRechitCluster_match_dtSeg_0p4,             "dtRechitCluster_match_dtSeg_0p4[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_dtSegTime_0p4",             dtRechitCluster_match_dtSegTime_0p4,             "dtRechitCluster_match_dtSegTime_0p4[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_dtSegTimeSpread_0p5",             dtRechitCluster_match_dtSegTimeSpread_0p5,             "dtRechitCluster_match_dtSegTimeSpread_0p5[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_dtSegTimeSpread_0p4",             dtRechitCluster_match_dtSegTimeSpread_0p4,             "dtRechitCluster_match_dtSegTimeSpread_0p4[nDtRechitClusters]/F");


          tree_->Branch("dtRechitCluster_match_dtSeg_sameStation_0p5",             dtRechitCluster_match_dtSeg_sameStation_0p5,             "dtRechitCluster_match_dtSeg_sameStation_0p5[nDtRechitClusters]/I");
          tree_->Branch("dtRechitCluster_match_dtSegTime_sameStation_0p5",             dtRechitCluster_match_dtSegTime_sameStation_0p5,             "dtRechitCluster_match_dtSegTime_sameStation_0p5[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_dtSeg_sameStation_0p4",             dtRechitCluster_match_dtSeg_sameStation_0p4,             "dtRechitCluster_match_dtSeg_sameStation_0p4[nDtRechitClusters]/I");
          tree_->Branch("dtRechitCluster_match_dtSegTime_sameStation_0p4",             dtRechitCluster_match_dtSegTime_sameStation_0p4,             "dtRechitCluster_match_dtSegTime_sameStation_0p4[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5",             dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5,             "dtRechitCluster_match_dtSegTimeSpread_sameStation_0p5[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4",             dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4,             "dtRechitCluster_match_dtSegTimeSpread_sameStation_0p4[nDtRechitClusters]/F");


          tree_->Branch("dtRechitCluster_match_RPCTime_dPhi0p5",             dtRechitCluster_match_RPCTime_dPhi0p5,             "dtRechitCluster_match_RPCTime_dPhi0p5[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_RPCTimeSpread_dPhi0p5",             dtRechitCluster_match_RPCTimeSpread_dPhi0p5,             "dtRechitCluster_match_RPCTimeSpread_dPhi0p5[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_RPCTime_dR0p4",             dtRechitCluster_match_RPCTime_dR0p4,             "dtRechitCluster_match_RPCTime_dR0p4[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_RPCTimeSpread_dR0p4",             dtRechitCluster_match_RPCTimeSpread_dR0p4,             "dtRechitCluster_match_RPCTimeSpread_dR0p4[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_RPChits_dR0p4",             dtRechitCluster_match_RPChits_dR0p4,             "dtRechitCluster_match_RPChits_dR0p4[nDtRechitClusters]/I");
          tree_->Branch("dtRechitCluster_match_RPCTime_sameStation_dR0p4",             dtRechitCluster_match_RPCTime_sameStation_dR0p4,             "dtRechitCluster_match_RPCTime_sameStation_dR0p4[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4",             dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4,             "dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_RPChits_sameStation_dR0p4",             dtRechitCluster_match_RPChits_sameStation_dR0p4,             "dtRechitCluster_match_RPChits_sameStation_dR0p4[nDtRechitClusters]/I");



        tree_->Branch("dtRechitCluster_match_gParticle_Id",             dtRechitCluster_match_gParticle_Id,             "dtRechitCluster_match_gParticle_Id[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_gParticle_Pt",             dtRechitCluster_match_gParticle_Pt,             "dtRechitCluster_match_gParticle_Pt[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gParticle_Eta",             dtRechitCluster_match_gParticle_Eta,             "dtRechitCluster_match_gParticle_Eta[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gParticle_Phi",             dtRechitCluster_match_gParticle_Phi,             "dtRechitCluster_match_gParticle_Phi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gParticle_E",             dtRechitCluster_match_gParticle_E,             "dtRechitCluster_match_gParticle_E[nDtRechitClusters]/F");
        tree_->Branch("dtRechitCluster_match_gParticle_Status",             dtRechitCluster_match_gParticle_Status,             "dtRechitCluster_match_gParticle_Status[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_gParticle_MotherId",             dtRechitCluster_match_gParticle_MotherId,             "dtRechitCluster_match_gParticle_MotherId[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_gParticle_deltaR",             dtRechitCluster_match_gParticle_deltaR,             "dtRechitCluster_match_gParticle_deltaR[nDtRechitClusters]/F");


        tree_->Branch("dtRechitCluster_match_RPChits_dPhi0p5",             dtRechitCluster_match_RPChits_dPhi0p5,             "dtRechitCluster_match_RPChits_dPhi0p5[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_RPCBx_dPhi0p5",             dtRechitCluster_match_RPCBx_dPhi0p5,             "dtRechitCluster_match_RPCBx_dPhi0p5[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_RB1_0p4",             dtRechitCluster_match_RB1_0p4,             "dtRechitCluster_match_RB1_0p4[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_RB1_dPhi0p5",             dtRechitCluster_match_RB1_dPhi0p5,             "dtRechitCluster_match_RB1_dPhi0p5[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_MB1Seg_0p4",             dtRechitCluster_match_MB1Seg_0p4,             "dtRechitCluster_match_MB1Seg_0p4[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_MB1Seg_0p5",             dtRechitCluster_match_MB1Seg_0p5,             "dtRechitCluster_match_MB1Seg_0p5[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_MB1hits_0p4",             dtRechitCluster_match_MB1hits_0p4,             "dtRechitCluster_match_MB1hits_0p4[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_MB1hits_0p5",             dtRechitCluster_match_MB1hits_0p5,             "dtRechitCluster_match_MB1hits_0p5[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_MB1hits_cosmics_plus",             dtRechitCluster_match_MB1hits_cosmics_plus,             "dtRechitCluster_match_MB1hits_cosmics_plus[nDtRechitClusters]/I");
        tree_->Branch("dtRechitCluster_match_MB1hits_cosmics_minus",             dtRechitCluster_match_MB1hits_cosmics_minus,             "dtRechitCluster_match_MB1hits_cosmics_minus[nDtRechitClusters]/I");




        // tree_->Branch("dtRechitCluster_match_dtSegT_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegT_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegX_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegX_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegY_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegY_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegZ_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegZ_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegEta_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegEta_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegPhi_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegPhi_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegWheel_dR0p4", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegWheel_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegWheel_dR0p4", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegWheel_dR0p4);
        //
        //
        // tree_->Branch("dtRechitCluster_match_dtSegT_sameStation_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegT_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegX_sameStation_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegX_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegY_sameStation_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegY_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegZ_sameStation_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegZ_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegEta_sameStation_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegEta_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegPhi_sameStation_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegPhi_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegWheel_sameStation_dR0p4", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegWheel_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_dtSegWheel_sameStation_dR0p4", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_dtSegWheel_sameStation_dR0p4);
        //
        //
        // tree_->Branch("dtRechitCluster_match_RPCBx_dPhi0p5", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCBx_dPhi0p5);
        // tree_->Branch("dtRechitCluster_match_RPCX_dPhi0p5", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCX_dPhi0p5);
        // tree_->Branch("dtRechitCluster_match_RPCY_dPhi0p5", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCY_dPhi0p5);
        // tree_->Branch("dtRechitCluster_match_RPCZ_dPhi0p5", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCZ_dPhi0p5);
        // tree_->Branch("dtRechitCluster_match_RPCPhi_dPhi0p5", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCPhi_dPhi0p5);
        // tree_->Branch("dtRechitCluster_match_RPCEta_dPhi0p5", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCEta_dPhi0p5);
        // tree_->Branch("dtRechitCluster_match_RPCRing_dPhi0p5", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_RPCRing_dPhi0p5);
        // tree_->Branch("dtRechitCluster_match_RPCLayer_dPhi0p5", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_RPCLayer_dPhi0p5);
        // tree_->Branch("dtRechitCluster_match_RPCSector_dPhi0p5", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_RPCSector_dPhi0p5);
        //
        //
        //
        // tree_->Branch("dtRechitCluster_match_RPCBx_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCBx_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCX_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCX_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCY_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCY_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCZ_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCZ_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCPhi_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCPhi_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCEta_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCEta_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCRing_dR0p4", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_RPCRing_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCLayer_dR0p4", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_RPCLayer_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCSector_dR0p4", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_RPCSector_dR0p4);
        //
        // tree_->Branch("dtRechitCluster_match_RPCBx_sameStation_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCBx_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCX_sameStation_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCX_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCY_sameStation_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCY_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCZ_sameStation_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCZ_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCPhi_sameStation_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCPhi_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCEta_sameStation_dR0p4", "std::vector<vector<float>>(nCscRechitClusters)",&dtRechitCluster_match_RPCEta_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCRing_sameStation_dR0p4", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_RPCRing_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCLayer_sameStation_dR0p4", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_RPCLayer_sameStation_dR0p4);
        // tree_->Branch("dtRechitCluster_match_RPCSector_sameStation_dR0p4", "std::vector<vector<int>>(nCscRechitClusters)",&dtRechitCluster_match_RPCSector_sameStation_dR0p4);


// correct DT clusters


tree_->Branch("dtRechitCluster2_match_gLLP",             dtRechitCluster2_match_gLLP,             "dtRechitCluster2_match_gLLP[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2_match_gLLP_minDeltaR",             dtRechitCluster2_match_gLLP_minDeltaR,             "dtRechitCluster2_match_gLLP_minDeltaR[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_index",             dtRechitCluster2_match_gLLP_index,             "dtRechitCluster2_match_gLLP_index[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_gLLP_eta",             dtRechitCluster2_match_gLLP_eta, "dtRechitCluster2_match_gLLP_eta[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_phi",             dtRechitCluster2_match_gLLP_phi, "dtRechitCluster2_match_gLLP_phi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_decay_r",             dtRechitCluster2_match_gLLP_decay_r, "dtRechitCluster2_match_gLLP_decay_r[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_decay_x",             dtRechitCluster2_match_gLLP_decay_x, "dtRechitCluster2_match_gLLP_decay_x[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_decay_y",             dtRechitCluster2_match_gLLP_decay_y, "dtRechitCluster2_match_gLLP_decay_y[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_decay_z",             dtRechitCluster2_match_gLLP_decay_z, "dtRechitCluster2_match_gLLP_decay_z[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_ctau",             dtRechitCluster2_match_gLLP_ctau, "dtRechitCluster2_match_gLLP_ctau[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_beta",             dtRechitCluster2_match_gLLP_beta, "dtRechitCluster2_match_gLLP_beta[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_csc",             dtRechitCluster2_match_gLLP_csc, "dtRechitCluster2_match_gLLP_csc[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2_match_gLLP_dt",             dtRechitCluster2_match_gLLP_dt, "dtRechitCluster2_match_gLLP_dt[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2_match_gLLP_multiplicity",             dtRechitCluster2_match_gLLP_multiplicity, "dtRechitCluster2_match_gLLP_multiplicity[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_gLLP_EM_multiplicity",             dtRechitCluster2_match_gLLP_EM_multiplicity, "dtRechitCluster2_match_gLLP_EM_multiplicity[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_gLLP_daughterKaon",             dtRechitCluster2_match_gLLP_daughterKaon, "dtRechitCluster2_match_gLLP_daughterKaon[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2_match_gLLP_e",             dtRechitCluster2_match_gLLP_e, "dtRechitCluster2_match_gLLP_e[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_pt",             dtRechitCluster2_match_gLLP_pt, "dtRechitCluster2_match_gLLP_pt[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_EMFracE",             dtRechitCluster2_match_gLLP_EMFracE, "dtRechitCluster2_match_gLLP_EMFracE[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_EMFracEz",             dtRechitCluster2_match_gLLP_EMFracEz, "dtRechitCluster2_match_gLLP_EMFracEz[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_EMFracP",             dtRechitCluster2_match_gLLP_EMFracP, "dtRechitCluster2_match_gLLP_EMFracP[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_EMFracPz",             dtRechitCluster2_match_gLLP_EMFracPz, "dtRechitCluster2_match_gLLP_EMFracPz[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_visE",             dtRechitCluster2_match_gLLP_visE, "dtRechitCluster2_match_gLLP_visE[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_visEz",             dtRechitCluster2_match_gLLP_visEz, "dtRechitCluster2_match_gLLP_visEz[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_visP",             dtRechitCluster2_match_gLLP_visP, "dtRechitCluster2_match_gLLP_visP[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_visPz",             dtRechitCluster2_match_gLLP_visPz, "dtRechitCluster2_match_gLLP_visPz[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_lepdPhi",             dtRechitCluster2_match_gLLP_lepdPhi, "dtRechitCluster2_match_gLLP_lepdPhi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_daughter0_deltaR",             dtRechitCluster2_match_gLLP_daughter0_deltaR, "dtRechitCluster2_match_gLLP_daughter0_deltaR[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_daughter1_deltaR",             dtRechitCluster2_match_gLLP_daughter1_deltaR, "dtRechitCluster2_match_gLLP_daughter1_deltaR[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_daughter2_deltaR",             dtRechitCluster2_match_gLLP_daughter2_deltaR, "dtRechitCluster2_match_gLLP_daughter2_deltaR[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_daughter3_deltaR",             dtRechitCluster2_match_gLLP_daughter3_deltaR, "dtRechitCluster2_match_gLLP_daughter3_deltaR[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_other_daughter_deltaR",             dtRechitCluster2_match_gLLP_other_daughter_deltaR, "dtRechitCluster2_match_gLLP_other_daughter_deltaR[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_other_daughter_index",             dtRechitCluster2_match_gLLP_other_daughter_index, "dtRechitCluster2_match_gLLP_other_daughter_index[nDtRechitClusters2]/I");


tree_->Branch("dtRechitCluster2_match_gLLP_other_eta",             dtRechitCluster2_match_gLLP_other_eta, "dtRechitCluster2_match_gLLP_other_eta[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_other_phi",             dtRechitCluster2_match_gLLP_other_phi, "dtRechitCluster2_match_gLLP_other_phi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_other_decay_r",             dtRechitCluster2_match_gLLP_other_decay_r, "dtRechitCluster2_match_gLLP_other_decay_r[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_other_decay_x",             dtRechitCluster2_match_gLLP_other_decay_x, "dtRechitCluster2_match_gLLP_other_decay_x[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_other_decay_y",             dtRechitCluster2_match_gLLP_other_decay_y, "dtRechitCluster2_match_gLLP_other_decay_y[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_other_decay_z",             dtRechitCluster2_match_gLLP_other_decay_z, "dtRechitCluster2_match_gLLP_other_decay_z[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_other_ctau",             dtRechitCluster2_match_gLLP_other_ctau, "dtRechitCluster2_match_gLLP_other_ctau[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_other_beta",             dtRechitCluster2_match_gLLP_other_beta, "dtRechitCluster2_match_gLLP_other_beta[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_other_csc",             dtRechitCluster2_match_gLLP_other_csc, "dtRechitCluster2_match_gLLP_other_csc[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2_match_gLLP_other_e",             dtRechitCluster2_match_gLLP_other_e, "dtRechitCluster2_match_gLLP_other_e[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gLLP_other_pt",             dtRechitCluster2_match_gLLP_other_pt, "dtRechitCluster2_match_gLLP_other_pt[nDtRechitClusters2]/F");

tree_->Branch("dtRechitCluster2X",             dtRechitCluster2X,             "dtRechitCluster2X[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2Y",             dtRechitCluster2Y,             "dtRechitCluster2Y[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2Z",             dtRechitCluster2Z,             "dtRechitCluster2Z[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2Time",             dtRechitCluster2Time,             "dtRechitCluster2Time[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2TimeWire",             dtRechitCluster2TimeWire,             "dtRechitCluster2TimeWire[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2TimeWirePruned",             dtRechitCluster2TimeWirePruned,             "dtRechitCluster2TimeWirePruned[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2TimeTotal",             dtRechitCluster2TimeTotal,             "dtRechitCluster2TimeTotal[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2Wheel",             dtRechitCluster2Wheel,             "dtRechitCluster2Wheel[nDtRechitClusters2]/I");

tree_->Branch("dtRechitCluster2TimeSpread",             dtRechitCluster2TimeSpread,             "dtRechitCluster2TimeSpread[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2TimeTotalSpread",             dtRechitCluster2TimeTotalSpread,             "dtRechitCluster2TimeTotalSpread[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2TimeTotalSpreadPruned",             dtRechitCluster2TimeTotalSpreadPruned,             "dtRechitCluster2TimeTotalSpreadPruned[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2TimeWireSpread",             dtRechitCluster2TimeWireSpread,             "dtRechitCluster2TimeWireSpread[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2GenMuonDeltaR",             dtRechitCluster2GenMuonDeltaR,             "dtRechitCluster2GenMuonDeltaR[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2XYSpread",             dtRechitCluster2XYSpread,             "dtRechitCluster2XYSpread[nDtRechitClusters2]/F");

tree_->Branch("dtRechitCluster2MajorAxis",             dtRechitCluster2MajorAxis,             "dtRechitCluster2MajorAxis[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MinorAxis",             dtRechitCluster2MinorAxis,             "dtRechitCluster2MinorAxis[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2EtaPhiSpread",             dtRechitCluster2EtaPhiSpread,             "dtRechitCluster2EtaPhiSpread[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2PhiSpread",             dtRechitCluster2PhiSpread,             "dtRechitCluster2PhiSpread[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2EtaSpread",             dtRechitCluster2EtaSpread,             "dtRechitCluster2EtaSpread[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2DeltaRSpread",             dtRechitCluster2DeltaRSpread,             "dtRechitCluster2DeltaRSpread[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2XSpread",             dtRechitCluster2XSpread,             "dtRechitCluster2XSpread[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2RSpread",             dtRechitCluster2RSpread,             "dtRechitCluster2RSpread[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2YSpread",             dtRechitCluster2YSpread,             "dtRechitCluster2YSpread[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2ZSpread",             dtRechitCluster2ZSpread,             "dtRechitCluster2ZSpread[nDtRechitClusters2]/F");


tree_->Branch("dtRechitCluster2Phi",             dtRechitCluster2Phi,             "dtRechitCluster2Phi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2Eta",             dtRechitCluster2Eta,             "dtRechitCluster2Eta[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2JetVetoPt",             dtRechitCluster2JetVetoPt,             "dtRechitCluster2JetVetoPt[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2JetVetoEta",             dtRechitCluster2JetVetoEta,             "dtRechitCluster2JetVetoEta[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2JetVetoPhi",             dtRechitCluster2JetVetoPhi,             "dtRechitCluster2JetVetoPhi[nDtRechitClusters2]/F");

tree_->Branch("dtRechitCluster2JetVetoE",             dtRechitCluster2JetVetoE,             "dtRechitCluster2JetVetoE[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2GenJetVetoPt",             dtRechitCluster2GenJetVetoPt,             "dtRechitCluster2GenJetVetoPt[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2GenJetVetoE",             dtRechitCluster2GenJetVetoE,             "dtRechitCluster2GenJetVetoE[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MuonVetoPt",             dtRechitCluster2MuonVetoPt,             "dtRechitCluster2MuonVetoPt[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MuonVetoE",             dtRechitCluster2MuonVetoE,             "dtRechitCluster2MuonVetoE[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MuonVetoPhi",             dtRechitCluster2MuonVetoPhi,             "dtRechitCluster2MuonVetoPhi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MuonVetoEta",             dtRechitCluster2MuonVetoEta,             "dtRechitCluster2MuonVetoEta[nDtRechitClusters2]/F");

tree_->Branch("dtRechitCluster2JetVetoElectronEnergyFraction",             dtRechitCluster2JetVetoElectronEnergyFraction,             "dtRechitCluster2JetVetoElectronEnergyFraction[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2JetVetoPhotonEnergyFraction",             dtRechitCluster2JetVetoPhotonEnergyFraction,             "dtRechitCluster2JetVetoPhotonEnergyFraction[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2JetVetoChargedHadronEnergyFraction",             dtRechitCluster2JetVetoChargedHadronEnergyFraction,             "dtRechitCluster2JetVetoChargedHadronEnergyFraction[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2JetVetoNeutralHadronEnergyFraction",             dtRechitCluster2JetVetoNeutralHadronEnergyFraction,             "dtRechitCluster2JetVetoNeutralHadronEnergyFraction[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2JetVetoMuonEnergyFraction",             dtRechitCluster2JetVetoMuonEnergyFraction,             "dtRechitCluster2JetVetoMuonEnergyFraction[nDtRechitClusters2]/F");


tree_->Branch("dtRechitCluster2JetVetoPt_0p6",             dtRechitCluster2JetVetoPt_0p6,        "dtRechitCluster2JetVetoPt_0p6[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2JetVetoPt_0p8",             dtRechitCluster2JetVetoPt_0p8,        "dtRechitCluster2JetVetoPt_0p8[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2JetVetoE_0p6",             dtRechitCluster2JetVetoE_0p6,        "dtRechitCluster2JetVetoE_0p6[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2JetVetoE_0p8",             dtRechitCluster2JetVetoE_0p8,        "dtRechitCluster2JetVetoE_0p8[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MuonVetoPt_0p6",             dtRechitCluster2MuonVetoPt_0p6,        "dtRechitCluster2MuonVetoPt_0p6[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MuonVetoPt_0p8",             dtRechitCluster2MuonVetoPt_0p8,        "dtRechitCluster2MuonVetoPt_0p8[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MuonVetoE_0p6",             dtRechitCluster2MuonVetoE_0p6,        "dtRechitCluster2MuonVetoE_0p6[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MuonVetoE_0p8",             dtRechitCluster2MuonVetoE_0p8,        "dtRechitCluster2MuonVetoE_0p8[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MuonVetoLooseIso",             dtRechitCluster2MuonVetoLooseIso,             "dtRechitCluster2MuonVetoLooseIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2MuonVetoTightIso",             dtRechitCluster2MuonVetoTightIso,             "dtRechitCluster2MuonVetoTightIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2MuonVetoVTightIso",             dtRechitCluster2MuonVetoVTightIso,             "dtRechitCluster2MuonVetoVTightIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2MuonVetoVVTightIso",             dtRechitCluster2MuonVetoVVTightIso,             "dtRechitCluster2MuonVetoVVTightIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2MuonVetoTightId",             dtRechitCluster2MuonVetoTightId,             "dtRechitCluster2MuonVetoTightId[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2MuonVetoLooseId",             dtRechitCluster2MuonVetoLooseId,             "dtRechitCluster2MuonVetoLooseId[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2MuonVetoGlobal",             dtRechitCluster2MuonVetoGlobal,             "dtRechitCluster2MuonVetoGlobal[nDtRechitClusters2]/O");


tree_->Branch("dtRechitCluster2MuonVetoIso",             dtRechitCluster2MuonVetoIso,             "dtRechitCluster2MuonVetoIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2IsoMuonVetoPt",             dtRechitCluster2IsoMuonVetoPt,             "dtRechitCluster2IsoMuonVetoPt[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2IsoMuonVetoE",             dtRechitCluster2IsoMuonVetoE,             "dtRechitCluster2IsoMuonVetoE[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2IsoMuonVetoPhi",             dtRechitCluster2IsoMuonVetoPhi,             "dtRechitCluster2IsoMuonVetoPhi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2IsoMuonVetoEta",             dtRechitCluster2IsoMuonVetoEta,             "dtRechitCluster2IsoMuonVetoEta[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2GenMuonVetoPt",             dtRechitCluster2GenMuonVetoPt,             "dtRechitCluster2GenMuonVetoPt[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MuonVetoType",             dtRechitCluster2MuonVetoType,             "dtRechitCluster2MuonVetoType[nDtRechitClusters2]/I");
// tree_->Branch("dtRechitCluster2GenMuonVetoE",             dtRechitCluster2GenMuonVetoE,             "dtRechitCluster2GenMuonVetoE[nDtRechitClusters2]/F");
// tree_->Branch("dtRechitCluster2GenMuonVetoProdX",             dtRechitCluster2GenMuonVetoProdX,             "dtRechitCluster2GenMuonVetoProdX[nDtRechitClusters2]/F");
// tree_->Branch("dtRechitCluster2GenMuonVetoProdY",             dtRechitCluster2GenMuonVetoProdY,             "dtRechitCluster2GenMuonVetoProdY[nDtRechitClusters2]/F");
// tree_->Branch("dtRechitCluster2GenMuonVetoProdZ",             dtRechitCluster2GenMuonVetoProdZ,             "dtRechitCluster2GenMuonVetoProdZ[nDtRechitClusters2]/F");
// tree_->Branch("dtRechitCluster2GenMuonVetoLLPDist",             dtRechitCluster2GenMuonVetoLLPDist,             "dtRechitCluster2GenMuonVetoLLPDist[nDtRechitClusters2]/F");
// tree_->Branch("dtRechitCluster2GenMuonVetoLLPIndex",             dtRechitCluster2GenMuonVetoLLPIndex,             "dtRechitCluster2GenMuonVetoLLPIndex[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2ZLep1",             dtRechitCluster2ZLep1,             "dtRechitCluster2ZLep1[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep2",             dtRechitCluster2ZLep2,             "dtRechitCluster2ZLep2[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep1Tag",             dtRechitCluster2ZLep1Tag,             "dtRechitCluster2ZLep1Tag[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep2Tag",             dtRechitCluster2ZLep2Tag,             "dtRechitCluster2ZLep2Tag[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep1Id",             dtRechitCluster2ZLep1Id,             "dtRechitCluster2ZLep1Id[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2ZLep2Id",             dtRechitCluster2ZLep2Id,             "dtRechitCluster2ZLep2Id[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster22ZLep1LooseIso",             dtRechitCluster2ZLep1LooseIso,             "dtRechitCluster2ZLep1LooseIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep1TightIso",             dtRechitCluster2ZLep1TightIso,             "dtRechitCluster2ZLep1TightIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep1VTightIso",             dtRechitCluster2ZLep1VTightIso,             "dtRechitCluster2ZLep1VTightIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep1VVTightIso",             dtRechitCluster2ZLep1VVTightIso,             "dtRechitCluster2ZLep1VVTightIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep1TightId",             dtRechitCluster2ZLep1TightId,             "dtRechitCluster2ZLep1TightId[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep2LooseIso",             dtRechitCluster2ZLep2LooseIso,             "dtRechitCluster2ZLep2LooseIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep2TightIso",             dtRechitCluster2ZLep2TightIso,             "dtRechitCluster2ZLep2TightIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep2VTightIso",             dtRechitCluster2ZLep2VTightIso,             "dtRechitCluster2ZLep2VTightIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep2VVTightIso",             dtRechitCluster2ZLep2VVTightIso,             "dtRechitCluster2ZLep2VVTightIso[nDtRechitClusters2]/O");
tree_->Branch("dtRechitCluster2ZLep2TightId",             dtRechitCluster2ZLep2TightId,             "dtRechitCluster2ZLep2TightId[nDtRechitClusters2]/O");


tree_->Branch("dtRechitCluster2Size",             dtRechitCluster2Size,             "dtRechitCluster2Size[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2Me11Ratio",             dtRechitCluster2Me11Ratio,             "dtRechitCluster2Me11Ratio[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2Me12Ratio",             dtRechitCluster2Me12Ratio,             "dtRechitCluster2Me12Ratio[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2NStation",             dtRechitCluster2NStation,             "dtRechitCluster2NStation[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NStation5",             dtRechitCluster2NStation5,             "dtRechitCluster2NStation5[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NStation10",             dtRechitCluster2NStation10,             "dtRechitCluster2NStation10[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NStation10perc",             dtRechitCluster2NStation10perc,             "dtRechitCluster2NStation10perc[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2AvgStation",             dtRechitCluster2AvgStation,             "dtRechitCluster2AvgStation[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2AvgStation5",             dtRechitCluster2AvgStation5,             "dtRechitCluster2AvgStation5[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2AvgStation10",             dtRechitCluster2AvgStation10,             "dtRechitCluster2AvgStation10[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2AvgStation10perc",             dtRechitCluster2AvgStation10perc,             "dtRechitCluster2AvgStation10perc[nDtRechitClusters2]/F");


tree_->Branch("dtRechitCluster2MaxStation",             dtRechitCluster2MaxStation,             "dtRechitCluster2MaxStation[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2MaxStationRatio",             dtRechitCluster2MaxStationRatio,             "dtRechitCluster2MaxStationRatio[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2NChamber",             dtRechitCluster2NChamber,             "dtRechitCluster2NChamber[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2MaxChamber",             dtRechitCluster2MaxChamber,             "dtRechitCluster2MaxChamber[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2MaxChamberRatio",             dtRechitCluster2MaxChamberRatio,             "dtRechitCluster2MaxChamberRatio[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2NSegmentStation1",             dtRechitCluster2NSegmentStation1,             "dtRechitCluster2NSegmentStation1[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NSegmentStation2",             dtRechitCluster2NSegmentStation2,             "dtRechitCluster2NSegmentStation2[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NSegmentStation3",             dtRechitCluster2NSegmentStation3,             "dtRechitCluster2NSegmentStation3[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NSegmentStation4",             dtRechitCluster2NSegmentStation4,             "dtRechitCluster2NSegmentStation4[nDtRechitClusters2]/I");

tree_->Branch("dtRechitCluster2NRechitChamberPlus11",             dtRechitCluster2NRechitChamberPlus11,             "dtRechitCluster2NRechitChamberPlus11[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberPlus12",             dtRechitCluster2NRechitChamberPlus12,             "dtRechitCluster2NRechitChamberPlus12[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberPlus13",             dtRechitCluster2NRechitChamberPlus13,             "dtRechitCluster2NRechitChamberPlus13[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberPlus21",             dtRechitCluster2NRechitChamberPlus21,             "dtRechitCluster2NRechitChamberPlus21[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberPlus22",             dtRechitCluster2NRechitChamberPlus22,             "dtRechitCluster2NRechitChamberPlus22[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberPlus31",             dtRechitCluster2NRechitChamberPlus31,             "dtRechitCluster2NRechitChamberPlus31[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberPlus32",             dtRechitCluster2NRechitChamberPlus32,             "dtRechitCluster2NRechitChamberPlus32[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberPlus41",             dtRechitCluster2NRechitChamberPlus41,             "dtRechitCluster2NRechitChamberPlus41[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberPlus42",             dtRechitCluster2NRechitChamberPlus42,             "dtRechitCluster2NRechitChamberPlus42[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberMinus11",             dtRechitCluster2NRechitChamberMinus11,             "dtRechitCluster2NRechitChamberMinus11[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberMinus12",             dtRechitCluster2NRechitChamberMinus12,             "dtRechitCluster2NRechitChamberMinus12[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberMinus13",             dtRechitCluster2NRechitChamberMinus13,             "dtRechitCluster2NRechitChamberMinus13[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberMinus21",             dtRechitCluster2NRechitChamberMinus21,             "dtRechitCluster2NRechitChamberMinus21[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberMinus22",             dtRechitCluster2NRechitChamberMinus22,             "dtRechitCluster2NRechitChamberMinus22[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberMinus31",             dtRechitCluster2NRechitChamberMinus31,             "dtRechitCluster2NRechitChamberMinus31[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberMinus32",             dtRechitCluster2NRechitChamberMinus32,             "dtRechitCluster2NRechitChamberMinus32[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberMinus41",             dtRechitCluster2NRechitChamberMinus41,             "dtRechitCluster2NRechitChamberMinus41[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NRechitChamberMinus42",             dtRechitCluster2NRechitChamberMinus42,             "dtRechitCluster2NRechitChamberMinus42[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2Met_dPhi",             dtRechitCluster2Met_dPhi,             "dtRechitCluster2Met_dPhi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MetXYCorr_dPhi",             dtRechitCluster2MetXYCorr_dPhi,             "dtRechitCluster2MetXYCorr_dPhi[nDtRechitClusters2]/F");

tree_->Branch("dtRechitCluster2NLayersChamberPlus11",             dtRechitCluster2NLayersChamberPlus11,             "dtRechitCluster2NLayersChamberPlus11[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberPlus12",             dtRechitCluster2NLayersChamberPlus12,             "dtRechitCluster2NLayersChamberPlus12[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberPlus13",             dtRechitCluster2NLayersChamberPlus13,             "dtRechitCluster2NLayersChamberPlus13[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberPlus21",             dtRechitCluster2NLayersChamberPlus21,             "dtRechitCluster2NLayersChamberPlus21[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberPlus22",             dtRechitCluster2NLayersChamberPlus22,             "dtRechitCluster2NLayersChamberPlus22[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberPlus31",             dtRechitCluster2NLayersChamberPlus31,             "dtRechitCluster2NLayersChamberPlus31[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberPlus32",             dtRechitCluster2NLayersChamberPlus32,             "dtRechitCluster2NLayersChamberPlus32[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberPlus41",             dtRechitCluster2NLayersChamberPlus41,             "dtRechitCluster2NLayersChamberPlus41[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberPlus42",             dtRechitCluster2NLayersChamberPlus42,             "dtRechitCluster2NLayersChamberPlus42[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberMinus11",             dtRechitCluster2NLayersChamberMinus11,             "dtRechitCluster2NLayersChamberMinus11[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberMinus12",             dtRechitCluster2NLayersChamberMinus12,             "dtRechitCluster2NLayersChamberMinus12[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberMinus13",             dtRechitCluster2NLayersChamberMinus13,             "dtRechitCluster2NLayersChamberMinus13[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberMinus21",             dtRechitCluster2NLayersChamberMinus21,             "dtRechitCluster2NLayersChamberMinus21[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberMinus22",             dtRechitCluster2NLayersChamberMinus22,             "dtRechitCluster2NLayersChamberMinus22[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberMinus31",             dtRechitCluster2NLayersChamberMinus31,             "dtRechitCluster2NLayersChamberMinus31[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberMinus32",             dtRechitCluster2NLayersChamberMinus32,             "dtRechitCluster2NLayersChamberMinus32[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberMinus41",             dtRechitCluster2NLayersChamberMinus41,             "dtRechitCluster2NLayersChamberMinus41[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2NLayersChamberMinus42",             dtRechitCluster2NLayersChamberMinus42,             "dtRechitCluster2NLayersChamberMinus42[nDtRechitClusters2]/I");


tree_->Branch("dtRechitCluster2MetHEM_dPhi",             dtRechitCluster2MetHEM_dPhi,             "dtRechitCluster2MetHEM_dPhi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MetHEMXYCorr_dPhi",             dtRechitCluster2MetHEMXYCorr_dPhi,             "dtRechitCluster2MetHEMXYCorr_dPhi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MetEENoise_dPhi",             dtRechitCluster2MetEENoise_dPhi,             "dtRechitCluster2MetEENoise_dPhi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MetEENoiseXYCorr_dPhi",             dtRechitCluster2MetEENoiseXYCorr_dPhi,             "dtRechitCluster2MetEENoiseXYCorr_dPhi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MetJesUp_dPhi",             dtRechitCluster2MetJesUp_dPhi,             "dtRechitCluster2MetJesUp_dPhi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2MetJesDown_dPhi",             dtRechitCluster2MetJesDown_dPhi,             "dtRechitCluster2MetJesDown_dPhi[nDtRechitClusters2]/F");


tree_->Branch("dtRechitCluster2_match_dtSeg_0p5",             dtRechitCluster2_match_dtSeg_0p5,             "dtRechitCluster2_match_dtSeg_0p5[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_dtSegTime_0p5",             dtRechitCluster2_match_dtSegTime_0p5,             "dtRechitCluster2_match_dtSegTime_0p5[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_dtSeg_0p4",             dtRechitCluster2_match_dtSeg_0p4,             "dtRechitCluster2_match_dtSeg_0p4[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_dtSegTime_0p4",             dtRechitCluster2_match_dtSegTime_0p4,             "dtRechitCluster2_match_dtSegTime_0p4[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_dtSegTimeSpread_0p5",             dtRechitCluster2_match_dtSegTimeSpread_0p5,             "dtRechitCluster2_match_dtSegTimeSpread_0p5[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_dtSegTimeSpread_0p4",             dtRechitCluster2_match_dtSegTimeSpread_0p4,             "dtRechitCluster2_match_dtSegTimeSpread_0p4[nDtRechitClusters2]/F");


  tree_->Branch("dtRechitCluster2_match_dtSeg_sameStation_0p5",             dtRechitCluster2_match_dtSeg_sameStation_0p5,             "dtRechitCluster2_match_dtSeg_sameStation_0p5[nDtRechitClusters2]/I");
  tree_->Branch("dtRechitCluster2_match_dtSegTime_sameStation_0p5",             dtRechitCluster2_match_dtSegTime_sameStation_0p5,             "dtRechitCluster2_match_dtSegTime_sameStation_0p5[nDtRechitClusters2]/F");
  tree_->Branch("dtRechitCluster2_match_dtSeg_sameStation_0p4",             dtRechitCluster2_match_dtSeg_sameStation_0p4,             "dtRechitCluster2_match_dtSeg_sameStation_0p4[nDtRechitClusters2]/I");
  tree_->Branch("dtRechitCluster2_match_dtSegTime_sameStation_0p4",             dtRechitCluster2_match_dtSegTime_sameStation_0p4,             "dtRechitCluster2_match_dtSegTime_sameStation_0p4[nDtRechitClusters2]/F");
  tree_->Branch("dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p5",             dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p5,             "dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p5[nDtRechitClusters2]/F");
  tree_->Branch("dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p4",             dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p4,             "dtRechitCluster2_match_dtSegTimeSpread_sameStation_0p4[nDtRechitClusters2]/F");


  tree_->Branch("dtRechitCluster2_match_RPCTime_dPhi0p5",             dtRechitCluster2_match_RPCTime_dPhi0p5,             "dtRechitCluster2_match_RPCTime_dPhi0p5[nDtRechitClusters2]/F");
  tree_->Branch("dtRechitCluster2_match_RPCTimeSpread_dPhi0p5",             dtRechitCluster2_match_RPCTimeSpread_dPhi0p5,             "dtRechitCluster2_match_RPCTimeSpread_dPhi0p5[nDtRechitClusters2]/F");
  tree_->Branch("dtRechitCluster2_match_RPCTime_dR0p4",             dtRechitCluster2_match_RPCTime_dR0p4,             "dtRechitCluster2_match_RPCTime_dR0p4[nDtRechitClusters2]/F");
  tree_->Branch("dtRechitCluster2_match_RPCTimeSpread_dR0p4",             dtRechitCluster2_match_RPCTimeSpread_dR0p4,             "dtRechitCluster2_match_RPCTimeSpread_dR0p4[nDtRechitClusters2]/F");
  tree_->Branch("dtRechitCluster2_match_RPChits_dR0p4",             dtRechitCluster2_match_RPChits_dR0p4,             "dtRechitCluster2_match_RPChits_dR0p4[nDtRechitClusters2]/I");
  tree_->Branch("dtRechitCluster2_match_RPCTime_sameStation_dR0p4",             dtRechitCluster2_match_RPCTime_sameStation_dR0p4,             "dtRechitCluster2_match_RPCTime_sameStation_dR0p4[nDtRechitClusters2]/F");
  tree_->Branch("dtRechitCluster2_match_RPCTimeSpread_sameStation_dR0p4",             dtRechitCluster2_match_RPCTimeSpread_sameStation_dR0p4,             "dtRechitCluster2_match_RPCTimeSpread_sameStation_dR0p4[nDtRechitClusters2]/F");
  tree_->Branch("dtRechitCluster2_match_RPChits_sameStation_dR0p4",             dtRechitCluster2_match_RPChits_sameStation_dR0p4,             "dtRechitCluster2_match_RPChits_sameStation_dR0p4[nDtRechitClusters2]/I");



tree_->Branch("dtRechitCluster2_match_gParticle_Id",             dtRechitCluster2_match_gParticle_Id,             "dtRechitCluster2_match_gParticle_Id[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_gParticle_Pt",             dtRechitCluster2_match_gParticle_Pt,             "dtRechitCluster2_match_gParticle_Pt[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gParticle_Eta",             dtRechitCluster2_match_gParticle_Eta,             "dtRechitCluster2_match_gParticle_Eta[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gParticle_Phi",             dtRechitCluster2_match_gParticle_Phi,             "dtRechitCluster2_match_gParticle_Phi[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gParticle_E",             dtRechitCluster2_match_gParticle_E,             "dtRechitCluster2_match_gParticle_E[nDtRechitClusters2]/F");
tree_->Branch("dtRechitCluster2_match_gParticle_Status",             dtRechitCluster2_match_gParticle_Status,             "dtRechitCluster2_match_gParticle_Status[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_gParticle_MotherId",             dtRechitCluster2_match_gParticle_MotherId,             "dtRechitCluster2_match_gParticle_MotherId[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_gParticle_deltaR",             dtRechitCluster2_match_gParticle_deltaR,             "dtRechitCluster2_match_gParticle_deltaR[nDtRechitClusters2]/F");


tree_->Branch("dtRechitCluster2_match_RPChits_dPhi0p5",             dtRechitCluster2_match_RPChits_dPhi0p5,             "dtRechitCluster2_match_RPChits_dPhi0p5[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_RPCBx_dPhi0p5",             dtRechitCluster2_match_RPCBx_dPhi0p5,             "dtRechitCluster2_match_RPCBx_dPhi0p5[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_RB1_0p4",             dtRechitCluster2_match_RB1_0p4,             "dtRechitCluster2_match_RB1_0p4[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_RB1_dPhi0p5",             dtRechitCluster2_match_RB1_dPhi0p5,             "dtRechitCluster2_match_RB1_dPhi0p5[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_MB1Seg_0p4",             dtRechitCluster2_match_MB1Seg_0p4,             "dtRechitCluster2_match_MB1Seg_0p4[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_MB1Seg_0p5",             dtRechitCluster2_match_MB1Seg_0p5,             "dtRechitCluster2_match_MB1Seg_0p5[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_MB1hits_0p4",             dtRechitCluster2_match_MB1hits_0p4,             "dtRechitCluster2_match_MB1hits_0p4[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_MB1hits_0p5",             dtRechitCluster2_match_MB1hits_0p5,             "dtRechitCluster2_match_MB1hits_0p5[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_MB1hits_cosmics_plus",             dtRechitCluster2_match_MB1hits_cosmics_plus,             "dtRechitCluster2_match_MB1hits_cosmics_plus[nDtRechitClusters2]/I");
tree_->Branch("dtRechitCluster2_match_MB1hits_cosmics_minus",             dtRechitCluster2_match_MB1hits_cosmics_minus,             "dtRechitCluster2_match_MB1hits_cosmics_minus[nDtRechitClusters2]/I");



  //gLLP branches
  tree_->Branch("gLLP_multiplicity",          gLLP_multiplicity,          "gLLP_multiplicity[2]/I");
  tree_->Branch("gLLP_multiplicity20",          gLLP_multiplicity20,          "gLLP_multiplicity20[2]/I");
  tree_->Branch("gLLP_EM_multiplicity",          gLLP_EM_multiplicity,          "gLLP_EM_multiplicity[2]/I");
  tree_->Branch("gLLP_eta",          gLLP_eta,          "gLLP_eta[2]/F");
  tree_->Branch("gLLP_phi",          gLLP_phi,          "gLLP_phi[2]/F");
  tree_->Branch("gLLP_csc",          gLLP_csc,          "gLLP_csc[2]/F");
  tree_->Branch("gLLP_dt",          gLLP_dt,          "gLLP_dt[2]/F");
  tree_->Branch("gLLP_beta",          gLLP_beta,          "gLLP_beta[2]/F");
  tree_->Branch("gLLP_maxMatchedDis",          gLLP_maxMatchedDis,          "gLLP_maxMatchedDis[2]/F");



  tree_->Branch("gLLP_e",          gLLP_e,          "gLLP_e[2]/F");
  tree_->Branch("gLLP_pt",          gLLP_pt,          "gLLP_pt[2]/F");
  tree_->Branch("gLLP_lepdPhi",          gLLP_lepdPhi,          "gLLP_lepdPhi[2]/F");
  tree_->Branch("gLLP_daughterKaon",          gLLP_daughterKaon,          "gLLP_daughterKaon[2]/O");

  tree_->Branch("gLLP_ctau",          gLLP_ctau,          "gLLP_ctau[2]/F");
  tree_->Branch("gLLP_EMFracE",          gLLP_EMFracE,          "gLLP_EMFracE[2]/F");
  tree_->Branch("gLLP_EMFracEz",          gLLP_EMFracEz,          "gLLP_EMFracEz[2]/F");
  tree_->Branch("gLLP_EMFracP",          gLLP_EMFracP,          "gLLP_EMFracP[2]/F");
  tree_->Branch("gLLP_EMFracPz",          gLLP_EMFracPz,          "gLLP_EMFracPz[2]/F");
  tree_->Branch("gLLP_visE",          gLLP_visE,          "gLLP_visE[2]/F");
  tree_->Branch("gLLP_visE20",          gLLP_visE20,          "gLLP_visE20[2]/F");
  tree_->Branch("gLLP_visEz",          gLLP_visEz,          "gLLP_visEz[2]/F");
  tree_->Branch("gLLP_visP",          gLLP_visP,          "gLLP_visP[2]/F");
  tree_->Branch("gLLP_visPz",          gLLP_visPz,          "gLLP_visPz[2]/F");
  tree_->Branch("gLLP_match_jet_pt", gLLP_match_jet_pt, "gLLP_match_jet_pt[2]/F");
  tree_->Branch("gLLP_match_jet_index", gLLP_match_jet_index, "gLLP_match_jet_index[2]/I");
  tree_->Branch("gLLP_match_jet_minDeltaR", gLLP_match_jet_minDeltaR, "gLLP_match_jet_minDeltaR[2]/F");
  tree_->Branch("gLLP_decay_vertex_r",          gLLP_decay_vertex_r,          "gLLP_decay_vertex_r[2]/F");
  tree_->Branch("gLLP_decay_vertex_x",          gLLP_decay_vertex_x,          "gLLP_decay_vertex_x[2]/F");
  tree_->Branch("gLLP_decay_vertex_y",          gLLP_decay_vertex_y,          "gLLP_decay_vertex_y[2]/F");
  tree_->Branch("gLLP_decay_vertex_z",          gLLP_decay_vertex_z,          "gLLP_decay_vertex_z[2]/F");
  tree_->Branch("gLLP_deltaR",          &gLLP_deltaR,          "gLLP_deltaR/F");
  tree_->Branch("gLLP_daughter_deltaR",          gLLP_daughter_deltaR,          "gLLP_daughter_deltaR[2]/F");

  tree_->Branch("gLLP_daughter_pt",          gLLP_daughter_pt,          "gLLP_daughter_pt[4]/F");
  tree_->Branch("gLLP_daughter_id",          gLLP_daughter_id,          "gLLP_daughter_id[4]/I");
  tree_->Branch("gLLP_daughter_eta",          gLLP_daughter_eta,          "gLLP_daughter_eta[4]/F");
  tree_->Branch("gLLP_daughter_phi",          gLLP_daughter_phi,          "gLLP_daughter_phi[4]/F");
  tree_->Branch("gLLP_daughter_e",          gLLP_daughter_e,          "gLLP_daughter_e[4]/F");
  tree_->Branch("gLLP_daughter_mass",          gLLP_daughter_mass,          "gLLP_daughter_mass[4]/F");



  //leptons
  tree_->Branch("nLeptons",  &nLeptons, "nLeptons/I");
  tree_->Branch("lepE",      lepE,      "lepE[nLeptons]/F");
  tree_->Branch("lepPt",     lepPt,     "lepPt[nLeptons]/F");
  tree_->Branch("lepEta",    lepEta,    "lepEta[nLeptons]/F");
  tree_->Branch("lepPhi",    lepPhi,    "lepPhi[nLeptons]/F");
  tree_->Branch("lepPdgId",  lepPdgId,  "lepPdgId[nLeptons]/I");
  tree_->Branch("lepDZ",     lepDZ,     "lepDZ[nLeptons]/F");
  tree_->Branch("lepPassId", lepPassId, "lepPassId[nLeptons]/O");
  tree_->Branch("lepFromZ", lepFromZ, "lepFromZ[nLeptons]/O");
  tree_->Branch("lepEff", lepEff, "lepEff[nLeptons]/F");
  tree_->Branch("lepSF", lepSF, "lepSF[nLeptons]/F");
  tree_->Branch("lepTriggerSF", lepTriggerSF, "lepTriggerSF[nLeptons]/F");
  tree_->Branch("lepTightIdSF", lepTightIdSF, "lepTightIdSF[nLeptons]/F");
  tree_->Branch("lepLooseIdSF", lepLooseIdSF, "lepLooseIdSF[nLeptons]/F");
  tree_->Branch("lepTightIsoSF", lepTightIsoSF, "lepTightIsoSF[nLeptons]/F");
  tree_->Branch("lepLooseIsoSF", lepLooseIsoSF, "lepLooseIsoSF[nLeptons]/F");
  // tree_->Branch("lepTriggerMCEfficiency", lepTriggerMCEfficiency, "lepTriggerMCEfficiency[nLeptons]/F");
  // tree_->Branch("lepTightIdMCEfficiency", lepTightIdMCEfficiency, "lepTightIdMCEfficiency[nLeptons]/F");
  // tree_->Branch("lepLooseIdMCEfficiency", lepLooseIdMCEfficiency, "lepLooseIdMCEfficiency[nLeptons]/F");
  // tree_->Branch("lepTightIsoMCEfficiency", lepTightIsoMCEfficiency, "lepTightIsoMCEfficiency[nLeptons]/F");
  // tree_->Branch("lepLooseIsoMCEfficiency", lepLooseIsoMCEfficiency, "lepLooseIsoMCEfficiency[nLeptons]/F");
  tree_->Branch("lepTag", lepTag, "lepTag[nLeptons]/O");



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
  tree_->Branch("ZMass1",      &ZMass1,  "ZMass1/F");

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

  tree_->Branch("jetPtJESUp",     jetPtJESUp,     "jetPtJESUp[nJets]/F");
  tree_->Branch("jetPtJESDown",     jetPtJESDown,     "jetPtJESDown[nJets]/F");
  tree_->Branch("jetEJESUp",      jetEJESUp,      "jetEJESUp[nJets]/F");
  tree_->Branch("jetEJESDown",      jetEJESDown,      "jetEJESDown[nJets]/F");
  tree_->Branch("JecUnc",      JecUnc,      "JecUnc[nJets]/F");

  tree_->Branch("jet_match_llp_pt", jet_match_llp_pt, "jet_match_llp_pt[nJets]/F");
  tree_->Branch("jet_match_llp_index", jet_match_llp_index, "jet_match_llp_index[nJets]/I");
  tree_->Branch("jet_match_llp_minDeltaR", jet_match_llp_minDeltaR, "jet_match_llp_minDeltaR[nJets]/F");
  tree_->Branch("jet_match_genJet_pt", jet_match_genJet_pt, "jet_match_genJet_pt[nJets]/F");
  tree_->Branch("jet_match_genJet_index", jet_match_genJet_index, "jet_match_genJet_index[nJets]/I");
  tree_->Branch("jet_match_genJet_minDeltaR", jet_match_genJet_minDeltaR, "jet_match_genJet_minDeltaR[nJets]/F");
  // tree_->Branch("jetLoosePassId", jetLoosePassId, "jetLoosePassId[nJets]/O");
  tree_->Branch("jetTightPassId", jetTightPassId, "jetTightPassId[nJets]/O");
  tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[982]/O"); //hardcoded
  tree_->Branch("METTrigger", METTrigger, "METTrigger/O"); //hardcoded
  tree_->Branch("METNoMuTrigger", METNoMuTrigger, "METNoMuTrigger/O"); //hardcoded

  // tree_->Branch("jetMuonEnergyFraction",   jetMuonEnergyFraction,   "jetMuonEnergyFraction[nJets]/F");
  tree_->Branch("jetChargedEMEnergyFraction",   jetChargedEMEnergyFraction,   "jetChargedEMEnergyFraction[nJets]/F");
  tree_->Branch("jetChargedHadronEnergyFraction",   jetChargedHadronEnergyFraction,   "jetChargedHadronEnergyFraction[nJets]/F");
  tree_->Branch("jetNeutralEMEnergyFraction",   jetNeutralEMEnergyFraction,   "jetNeutralEMEnergyFraction[nJets]/F");
  tree_->Branch("jetNeutralHadronEnergyFraction",   jetNeutralHadronEnergyFraction,   "jetNeutralHadronEnergyFraction[nJets]/F");
};
