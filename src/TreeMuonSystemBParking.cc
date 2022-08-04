#include "RazorHelper.h"
#include "TreeMuonSystemBParking.h"
#include "assert.h"
#include "TTree.h"
#include "DBSCAN.h"

// Constructor
TreeMuonSystemBParking::TreeMuonSystemBParking()
{
  InitVariables();
};
TreeMuonSystemBParking::~TreeMuonSystemBParking()
{
  if (f_) f_->Close();
};
void TreeMuonSystemBParking::InitVariables()
{
  runNum=0; lumiSec=0; evtNum=0;  MC_condition = 0;npv=0;
  rho=-1;

  Flag2_HBHENoiseFilter = false; Flag2_HBHEIsoNoiseFilter = false; Flag2_BadPFMuonFilter = false; Flag2_globalSuperTightHalo2016Filter = false;
  Flag2_globalTightHalo2016Filter = false; Flag2_BadChargedCandidateFilter = false; Flag2_EcalDeadCellTriggerPrimitiveFilter = false; Flag2_ecalBadCalibFilter = false;
  Flag2_eeBadScFilter = false;
  Flag2_all = false;

 metEENoise = -999.;metPhiEENoise = -999.;


   gLLP_eta = 0.0;
   gLLP_phi = 0.0;
   gLLP_beta = 0.0;
   gLLP_e = 0.0;
   gLLP_pt = 0.0;
   gLLP_csc = 0.0;
   gLLP_dt = 0.0;
   gLLP_ctau = 0.0;
   gLLP_decay_vertex_r = 0.0;
   gLLP_decay_vertex_x = 0.0;
   gLLP_decay_vertex_y = 0.0;
   gLLP_decay_vertex_z = 0.0;

   nLeptons = 0;
   for( int i = 0; i < N_MAX_LEPTONS; i++ )
   {

     lepE[i]      = -999.;
     lepPt[i]     = -999.;
     lepEta[i]    = -999.;
     lepPhi[i]    = -999.;
     lepPdgId[i]  = -999;
     lepDZ[i]     = -999.;

     lepLooseId[i] = false;
     lepTightId[i] = false;

     lepPassLooseIso[i] = false;
     lepPassTightIso[i] = false;
     lepPassVTightIso[i] = false;
     lepPassVTightIso[i] = false;


   }

   nJets = 0;
   for( int i = 0; i < N_MAX_JETS; i++ )
   {
     jetE[i]      = -999.;
     jetPt[i]     = -999.;
     jetEta[i]    = -999.;
     jetPhi[i]    = -999.;
     jetTightPassId[i] = false;
   }

   for(int i = 0; i <NTriggersMAX; i++){
     HLTDecision[i] = false;
   }
    //CSC

  nCscRechitClusters = 0;
  // nDtRechitClusters = 0;

  nCscRechits = 0;
  nEarlyCscRechits = 0;
  nLateCscRechits = 0;
  nEarly2CscRechits = 0;
  nLate2CscRechits = 0;
  nCscPositiveYRechits = 0;
  nCscNegativeYRechits = 0;
  cscPosTpeak = 0.0;
  cscNegTpeak = 0.0;
  // nDtRings = 0;
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


  /*nDTNegativeYRechits  = 0;
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
  */


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

      cscRechitCluster_match_gLLP_deltaR[i] = 999.;


        cscRechitClusterSize[i] = -999;
        cscRechitClusterX[i] = -999.;
        cscRechitClusterY[i] = -999.;
        cscRechitClusterZ[i] = -999.;
        cscRechitClusterTime[i] = -999.;
        cscRechitClusterTimeWeighted[i] = -999.;
        cscRechitClusterTimeTotal[i] = -999.;


        cscRechitClusterGenMuonDeltaR[i] = 999.;
        cscRechitClusterTimeSpread[i] = -999.;
        cscRechitClusterTimeSpreadWeighted[i] = -999.;
        cscRechitClusterTimeSpreadWeightedAll[i] = -999.;


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


        cscRechitClusterJetVetoPt[i] = 0.0;
        cscRechitClusterJetVetoEta[i] = 0.0;
        cscRechitClusterJetVetoPhi[i] = 0.0;
        cscRechitClusterJetVetoE[i] = 0.0;
        cscRechitClusterMuonVetoPt[i] = 0.0;
        cscRechitClusterMuonVetoE[i] = 0.0;
        cscRechitClusterMuonVetoPhi[i] = 0.0;
        cscRechitClusterMuonVetoEta[i] = 0.0;
        cscRechitClusterMuonVetoLooseId[i] = false;
        cscRechitClusterMuonVetoGlobal[i] = false;
        cscRechitClusterMetEENoise_dPhi[i] = 999.;
        // cscRechitClusterNLayersChamberPlus11[i] = -999;
        // cscRechitClusterNLayersChamberPlus12[i] = -999;
        // cscRechitClusterNLayersChamberPlus13[i] = -999;
        // cscRechitClusterNLayersChamberPlus21[i] = -999;
        // cscRechitClusterNLayersChamberPlus22[i] = -999;
        // cscRechitClusterNLayersChamberPlus31[i] = -999;
        // cscRechitClusterNLayersChamberPlus32[i] = -999;
        // cscRechitClusterNLayersChamberPlus41[i] = -999;
        // cscRechitClusterNLayersChamberPlus42[i] = -999;
        // cscRechitClusterNLayersChamberMinus11[i] = -999;
        // cscRechitClusterNLayersChamberMinus12[i] = -999;
        // cscRechitClusterNLayersChamberMinus13[i] = -999;
        // cscRechitClusterNLayersChamberMinus21[i] = -999;
        // cscRechitClusterNLayersChamberMinus22[i] = -999;
        // cscRechitClusterNLayersChamberMinus31[i] = -999;
        // cscRechitClusterNLayersChamberMinus32[i] = -999;
        // cscRechitClusterNLayersChamberMinus41[i] = -999;
        // cscRechitClusterNLayersChamberMinus42[i] = -999;

        dtRechitClusterMaxDPhi[i] = 0.;
        dtRechitClusterMaxDPhi_index[i] = 999;
        dtRechitClusterNSegStation1[i] = 0;
        dtRechitClusterNSegStation2[i] = 0;
        dtRechitClusterNSegStation3[i] = 0;
        dtRechitClusterNSegStation4[i] = 0;

        dtRechitClusterNOppositeSegStation1[i] = 0;
        dtRechitClusterNOppositeSegStation2[i] = 0;
        dtRechitClusterNOppositeSegStation3[i] = 0;
        dtRechitClusterNOppositeSegStation4[i] = 0;


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



        dtRechitCluster_match_RPCTime_dPhi0p5[i] = 0;
        dtRechitCluster_match_RPCTimeSpread_dPhi0p5[i] = 0;
        dtRechitCluster_match_RPCTime_dR0p4[i] = 0;
        dtRechitCluster_match_RPCTimeSpread_dR0p4[i] = 0;
        dtRechitCluster_match_RPChits_dR0p4[i] = 0;
        dtRechitCluster_match_RPCTime_sameStation_dR0p4[i] = 0;
        dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[i] = 0;
        dtRechitCluster_match_RPChits_sameStation_dR0p4[i] = 0;
        dtRechitCluster_match_gLLP_deltaR[i] = 999.;
        dtRechitClusterJetVetoPt[i] = 0.0;
        dtRechitClusterJetVetoEta[i] = 0.0;
        dtRechitClusterJetVetoPhi[i] = 0.0;
        dtRechitClusterJetVetoE[i] = 0.0;
        dtRechitClusterMuonVetoPt[i] = 0.0;
        dtRechitClusterMuonVetoE[i] = 0.0;
        dtRechitClusterMuonVetoPhi[i] = 0.0;
        dtRechitClusterMuonVetoEta[i] = 0.0;
        dtRechitClusterMuonVetoLooseId[i] = false;
        dtRechitClusterMuonVetoGlobal[i] = false;
        dtRechitClusterMetEENoise_dPhi[i] = 999.;


          dtRechitClusterSize[i] = -999;
          dtRechitClusterOverlap[i] = false;
          dtRechitClusterNoiseHit[i] = 0;
          dtRechitClusterNoiseHitStation1[i] = 0;
          dtRechitClusterNoiseHitStation2[i] = 0;
          dtRechitClusterNoiseHitStation3[i] = 0;
          dtRechitClusterNoiseHitStation4[i] = 0;
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

      // for(int j = 0;j<N_MAX_CSC;j++)
      // {


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


};

void TreeMuonSystemBParking::InitTree()
{
  assert(tree_);
  InitVariables();

  tree_->SetBranchAddress("runNum",      &runNum);
  tree_->SetBranchAddress("MC_condition",      &MC_condition);
  tree_->SetBranchAddress("lumiSec",     &lumiSec);
  tree_->SetBranchAddress("evtNum",      &evtNum);



  tree_->SetBranchAddress("npv",         &npv);
  tree_->SetBranchAddress("rho",         &rho);

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


    tree_->SetBranchAddress("metEENoise",      &metEENoise);
    tree_->SetBranchAddress("metPhiEENoise",      &metPhiEENoise);


      tree_->SetBranchAddress("gLLP_eta",    &gLLP_eta);
      tree_->SetBranchAddress("gLLP_phi",    &gLLP_phi);
      tree_->SetBranchAddress("gLLP_csc",    &gLLP_csc);
      tree_->SetBranchAddress("gLLP_dt",    &gLLP_dt);
      tree_->SetBranchAddress("gLLP_ctau",    &gLLP_ctau);
      tree_->SetBranchAddress("gLLP_beta",    &gLLP_beta);
      tree_->SetBranchAddress("gLLP_e",    &gLLP_e);
      tree_->SetBranchAddress("gLLP_pt",    &gLLP_pt);
      tree_->SetBranchAddress("gLLP_decay_vertex_r",    &gLLP_decay_vertex_r);
      tree_->SetBranchAddress("gLLP_decay_vertex_x",    &gLLP_decay_vertex_x);
      tree_->SetBranchAddress("gLLP_decay_vertex_y",    &gLLP_decay_vertex_y);
      tree_->SetBranchAddress("gLLP_decay_vertex_z",    &gLLP_decay_vertex_z);

      tree_->SetBranchAddress("nLeptons",    &nLeptons);
      tree_->SetBranchAddress("lepE",        lepE);
      tree_->SetBranchAddress("lepPt",       lepPt);
      tree_->SetBranchAddress("lepEta",      lepEta);
      tree_->SetBranchAddress("lepPhi",      lepPhi);
      tree_->SetBranchAddress("lepPdgId",  lepPdgId);
      tree_->SetBranchAddress("lepDZ",     lepDZ);

      tree_->SetBranchAddress("lepLooseId", lepLooseId);
      tree_->SetBranchAddress("lepTightId", lepTightId);

      tree_->SetBranchAddress("lepPassLooseIso", lepPassLooseIso);
      tree_->SetBranchAddress("lepPassTightIso", lepPassTightIso);
      tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);
      tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);

      //jets
      tree_->SetBranchAddress("nJets",     &nJets);
      tree_->SetBranchAddress("jetE",      jetE);
      tree_->SetBranchAddress("jetPt",     jetPt);
      tree_->SetBranchAddress("jetEta",    jetEta);
      tree_->SetBranchAddress("jetPhi",    jetPhi);
      tree_->SetBranchAddress("jetTightPassId", jetTightPassId);
      // triggers
      tree_->SetBranchAddress("HLTDecision",   HLTDecision);

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



  tree_->SetBranchAddress("nDtRechitClusters",             &nDtRechitClusters);

tree_->SetBranchAddress("dtRechitClusterMaxDPhi",             &dtRechitClusterMaxDPhi);
tree_->SetBranchAddress("dtRechitClusterMaxDPhi_index",             &dtRechitClusterMaxDPhi_index);

  tree_->SetBranchAddress("dtRechitClusterNSegStation1",             &dtRechitClusterNSegStation1);
  tree_->SetBranchAddress("dtRechitClusterNSegStation2",             &dtRechitClusterNSegStation2);
  tree_->SetBranchAddress("dtRechitClusterNSegStation3",             &dtRechitClusterNSegStation3);
  tree_->SetBranchAddress("dtRechitClusterNSegStation4",             &dtRechitClusterNSegStation4);

  tree_->SetBranchAddress("dtRechitClusterNOppositeSegStation1",             &dtRechitClusterNOppositeSegStation1);
  tree_->SetBranchAddress("dtRechitClusterNOppositeSegStation2",             &dtRechitClusterNOppositeSegStation2);
  tree_->SetBranchAddress("dtRechitClusterNOppositeSegStation3",             &dtRechitClusterNOppositeSegStation3);
  tree_->SetBranchAddress("dtRechitClusterNOppositeSegStation4",             &dtRechitClusterNOppositeSegStation4);

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


  tree_->SetBranchAddress("dtRechitClusterJetVetoPt",             dtRechitClusterJetVetoPt);
  tree_->SetBranchAddress("dtRechitClusterJetVetoEta",             dtRechitClusterJetVetoEta);
  tree_->SetBranchAddress("dtRechitClusterJetVetoE",             dtRechitClusterJetVetoE);
  tree_->SetBranchAddress("dtRechitClusterJetVetoPhi",             dtRechitClusterJetVetoPhi);

  tree_->SetBranchAddress("dtRechitClusterMuonVetoPt",             dtRechitClusterMuonVetoPt);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoE",             dtRechitClusterMuonVetoE);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoPhi",             dtRechitClusterMuonVetoPhi);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoEta",             dtRechitClusterMuonVetoEta);

  tree_->SetBranchAddress("dtRechitClusterMuonVetoLooseId",             dtRechitClusterMuonVetoLooseId);
  tree_->SetBranchAddress("dtRechitClusterMuonVetoGlobal",             dtRechitClusterMuonVetoGlobal);
  tree_->SetBranchAddress("dtRechitClusterMetEENoise_dPhi",             dtRechitClusterMetEENoise_dPhi);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_deltaR",             dtRechitCluster_match_gLLP_deltaR);


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

  tree_->SetBranchAddress("dtRechitClusterSize",             dtRechitClusterSize);
  tree_->SetBranchAddress("dtRechitClusterOverlap",             dtRechitClusterOverlap);
  tree_->SetBranchAddress("dtRechitClusterNoiseHit",             dtRechitClusterNoiseHit);
  tree_->SetBranchAddress("dtRechitClusterNoiseHitStation1",             dtRechitClusterNoiseHitStation1);
  tree_->SetBranchAddress("dtRechitClusterNoiseHitStation2",             dtRechitClusterNoiseHitStation2);
  tree_->SetBranchAddress("dtRechitClusterNoiseHitStation3",             dtRechitClusterNoiseHitStation3);
  tree_->SetBranchAddress("dtRechitClusterNoiseHitStation4",             dtRechitClusterNoiseHitStation4);

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


    tree_->SetBranchAddress("dtRechitCluster_match_RPCTime_dPhi0p5",             dtRechitCluster_match_RPCTime_dPhi0p5);
    tree_->SetBranchAddress("dtRechitCluster_match_RPCTimeSpread_dPhi0p5",             dtRechitCluster_match_RPCTimeSpread_dPhi0p5);
    tree_->SetBranchAddress("dtRechitCluster_match_RPCTime_dR0p4",             dtRechitCluster_match_RPCTime_dR0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_RPCTimeSpread_dR0p4",             dtRechitCluster_match_RPCTimeSpread_dR0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_RPChits_dR0p4",             dtRechitCluster_match_RPChits_dR0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_RPCTime_sameStation_dR0p4",             dtRechitCluster_match_RPCTime_sameStation_dR0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4",             dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_RPChits_sameStation_dR0p4",             dtRechitCluster_match_RPChits_sameStation_dR0p4);







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
  tree_->SetBranchAddress("cscRechitClusterTimeWeighted",             cscRechitClusterTimeWeighted);
  tree_->SetBranchAddress("cscRechitClusterTimeTotal",             cscRechitClusterTimeTotal);

  tree_->SetBranchAddress("cscRechitClusterGenMuonDeltaR",             cscRechitClusterGenMuonDeltaR);

  tree_->SetBranchAddress("cscRechitClusterTimeSpread",             cscRechitClusterTimeSpread);
  tree_->SetBranchAddress("cscRechitClusterTimeSpreadWeighted",             cscRechitClusterTimeSpreadWeighted);
  tree_->SetBranchAddress("cscRechitClusterTimeSpreadWeightedAll",             cscRechitClusterTimeSpreadWeightedAll);

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


  tree_->SetBranchAddress("cscRechitClusterJetVetoPt",             cscRechitClusterJetVetoPt);
  tree_->SetBranchAddress("cscRechitClusterJetVetoEta",             cscRechitClusterJetVetoEta);
  tree_->SetBranchAddress("cscRechitClusterJetVetoE",             cscRechitClusterJetVetoE);
  tree_->SetBranchAddress("cscRechitClusterJetVetoPhi",             cscRechitClusterJetVetoPhi);

  tree_->SetBranchAddress("cscRechitClusterMuonVetoPt",             cscRechitClusterMuonVetoPt);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoE",             cscRechitClusterMuonVetoE);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoPhi",             cscRechitClusterMuonVetoPhi);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoEta",             cscRechitClusterMuonVetoEta);

  tree_->SetBranchAddress("cscRechitClusterMuonVetoLooseId",             cscRechitClusterMuonVetoLooseId);
  tree_->SetBranchAddress("cscRechitClusterMuonVetoGlobal",             cscRechitClusterMuonVetoGlobal);
  tree_->SetBranchAddress("cscRechitClusterMetEENoise_dPhi",             cscRechitClusterMetEENoise_dPhi);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus11",             cscRechitClusterNLayersChamberPlus11);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus12",             cscRechitClusterNLayersChamberPlus12);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus13",             cscRechitClusterNLayersChamberPlus13);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus21",             cscRechitClusterNLayersChamberPlus21);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus22",             cscRechitClusterNLayersChamberPlus22);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus31",             cscRechitClusterNLayersChamberPlus31);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus32",             cscRechitClusterNLayersChamberPlus32);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus41",             cscRechitClusterNLayersChamberPlus41);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberPlus42",             cscRechitClusterNLayersChamberPlus42);
  //
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus11",             cscRechitClusterNLayersChamberMinus11);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus12",             cscRechitClusterNLayersChamberMinus12);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus13",             cscRechitClusterNLayersChamberMinus13);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus21",             cscRechitClusterNLayersChamberMinus21);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus22",             cscRechitClusterNLayersChamberMinus22);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus31",             cscRechitClusterNLayersChamberMinus31);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus32",             cscRechitClusterNLayersChamberMinus32);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus41",             cscRechitClusterNLayersChamberMinus41);
  // tree_->SetBranchAddress("cscRechitClusterNLayersChamberMinus42",             cscRechitClusterNLayersChamberMinus42);

  tree_->SetBranchAddress("cscRechitCluster_match_dtRechits_0p4",             cscRechitCluster_match_dtRechits_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_MB1_0p4",             cscRechitCluster_match_MB1_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_dtRechits_phi0p2",             cscRechitCluster_match_dtRechits_phi0p2);


  tree_->SetBranchAddress("cscRechitCluster_match_dtSeg_0p4",             cscRechitCluster_match_dtSeg_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_MB1Seg_0p4",             cscRechitCluster_match_MB1Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_RB1_0p4",             cscRechitCluster_match_RB1_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_RE12_0p4",             cscRechitCluster_match_RE12_0p4);

  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_deltaR",             cscRechitCluster_match_gLLP_deltaR);

};

void TreeMuonSystemBParking::LoadTree(const char* file)
{
  f_ = TFile::Open(file);
  assert(f_);
  tree_ = dynamic_cast<TTree*>(f_->Get("MuonSystem"));
  InitTree();
  assert(tree_);
};

void TreeMuonSystemBParking::CreateTree()
{
  tree_ = new TTree("MuonSystem","MuonSystem");
  f_ = 0;

  tree_->Branch("runNum",      &runNum,     "runNum/i");      // event run number
  tree_->Branch("MC_condition",      &MC_condition,     "MC_condition/i");      // event run number
  tree_->Branch("lumiSec",     &lumiSec,    "lumiSec/i");     // event lumi section
  tree_->Branch("evtNum",      &evtNum,     "evtNum/i");      // event number
  tree_->Branch("rho",         &rho,        "rho/F");

    tree_->Branch("npv",         &npv,        "npv/i");         // number of primary vertices


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
    tree_->Branch("metEENoise",      &metEENoise,     "metEENoise/F");      // phi(MET)
    tree_->Branch("metPhiEENoise",      &metPhiEENoise,     "metPhiEENoise/F");      // phi(MET)


      tree_->Branch("gLLP_eta",          &gLLP_eta,          "gLLP_eta/F");
      tree_->Branch("gLLP_phi",          &gLLP_phi,          "gLLP_phi/F");
      tree_->Branch("gLLP_csc",          &gLLP_csc,          "gLLP_csc/F");
      tree_->Branch("gLLP_dt",          &gLLP_dt,          "gLLP_dt/F");
      tree_->Branch("gLLP_beta",          &gLLP_beta,          "gLLP_beta/F");
      tree_->Branch("gLLP_e",          &gLLP_e,          "gLLP_e/F");
      tree_->Branch("gLLP_pt",          &gLLP_pt,          "gLLP_pt/F");
      tree_->Branch("gLLP_ctau",          &gLLP_ctau,          "gLLP_ctau/F");
      tree_->Branch("gLLP_decay_vertex_r",          &gLLP_decay_vertex_r,          "gLLP_decay_vertex_r/F");
      tree_->Branch("gLLP_decay_vertex_x",          &gLLP_decay_vertex_x,          "gLLP_decay_vertex_x/F");
      tree_->Branch("gLLP_decay_vertex_y",          &gLLP_decay_vertex_y,          "gLLP_decay_vertex_y/F");
      tree_->Branch("gLLP_decay_vertex_z",          &gLLP_decay_vertex_z,          "gLLP_decay_vertex_z/F");

      //leptons
      tree_->Branch("nLeptons",  &nLeptons, "nLeptons/I");
      tree_->Branch("lepE",      lepE,      "lepE[nLeptons]/F");
      tree_->Branch("lepPt",     lepPt,     "lepPt[nLeptons]/F");
      tree_->Branch("lepEta",    lepEta,    "lepEta[nLeptons]/F");
      tree_->Branch("lepPhi",    lepPhi,    "lepPhi[nLeptons]/F");
      tree_->Branch("lepPdgId",  lepPdgId,  "lepPdgId[nLeptons]/I");
      tree_->Branch("lepDZ",     lepDZ,     "lepDZ[nLeptons]/F");
      tree_->Branch("lepLooseId", lepLooseId, "lepLooseId[nLeptons]/O");
      tree_->Branch("lepTightId", lepTightId, "lepTightId[nLeptons]/O");

      //jets
      tree_->Branch("nJets",     &nJets,    "nJets/I");
      tree_->Branch("jetE",      jetE,      "jetE[nJets]/F");
      tree_->Branch("jetPt",     jetPt,     "jetPt[nJets]/F");
      tree_->Branch("jetEta",    jetEta,    "jetEta[nJets]/F");
      tree_->Branch("jetPhi",    jetPhi,    "jetPhi[nJets]/F");
      tree_->Branch("jetTightPassId", jetTightPassId, "jetTightPassId[nJets]/O");

      tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[1201]/O"); //hardcoded

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



    tree_->Branch("nCscRechitClusters",             &nCscRechitClusters, "nCscRechitClusters/I");

    tree_->Branch("cscRechitClusterX",             cscRechitClusterX,             "cscRechitClusterX[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterY",             cscRechitClusterY,             "cscRechitClusterY[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterZ",             cscRechitClusterZ,             "cscRechitClusterZ[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTime",             cscRechitClusterTime,             "cscRechitClusterTime[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeWeighted",             cscRechitClusterTimeWeighted,             "cscRechitClusterTimeWeighted[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeTotal",             cscRechitClusterTimeTotal,             "cscRechitClusterTimeTotal[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterTimeSpread",             cscRechitClusterTimeSpread,             "cscRechitClusterTimeSpread[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeSpreadWeighted",             cscRechitClusterTimeSpreadWeighted,             "cscRechitClusterTimeSpreadWeighted[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterTimeSpreadWeightedAll",             cscRechitClusterTimeSpreadWeightedAll,             "cscRechitClusterTimeSpreadWeightedAll[nCscRechitClusters]/F");
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
    tree_->Branch("cscRechitCluster_match_cscRechits_0p4",             cscRechitCluster_match_cscRechits_0p4,             "cscRechitCluster_match_cscRechits_0p4[nCscRechitClusters]/I");

    // tree_->Branch("cscRechitClusterNLayersChamberPlus11",             cscRechitClusterNLayersChamberPlus11,             "cscRechitClusterNLayersChamberPlus11[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberPlus12",             cscRechitClusterNLayersChamberPlus12,             "cscRechitClusterNLayersChamberPlus12[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberPlus13",             cscRechitClusterNLayersChamberPlus13,             "cscRechitClusterNLayersChamberPlus13[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberPlus21",             cscRechitClusterNLayersChamberPlus21,             "cscRechitClusterNLayersChamberPlus21[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberPlus22",             cscRechitClusterNLayersChamberPlus22,             "cscRechitClusterNLayersChamberPlus22[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberPlus31",             cscRechitClusterNLayersChamberPlus31,             "cscRechitClusterNLayersChamberPlus31[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberPlus32",             cscRechitClusterNLayersChamberPlus32,             "cscRechitClusterNLayersChamberPlus32[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberPlus41",             cscRechitClusterNLayersChamberPlus41,             "cscRechitClusterNLayersChamberPlus41[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberPlus42",             cscRechitClusterNLayersChamberPlus42,             "cscRechitClusterNLayersChamberPlus42[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberMinus11",             cscRechitClusterNLayersChamberMinus11,             "cscRechitClusterNLayersChamberMinus11[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberMinus12",             cscRechitClusterNLayersChamberMinus12,             "cscRechitClusterNLayersChamberMinus12[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberMinus13",             cscRechitClusterNLayersChamberMinus13,             "cscRechitClusterNLayersChamberMinus13[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberMinus21",             cscRechitClusterNLayersChamberMinus21,             "cscRechitClusterNLayersChamberMinus21[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberMinus22",             cscRechitClusterNLayersChamberMinus22,             "cscRechitClusterNLayersChamberMinus22[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberMinus31",             cscRechitClusterNLayersChamberMinus31,             "cscRechitClusterNLayersChamberMinus31[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberMinus32",             cscRechitClusterNLayersChamberMinus32,             "cscRechitClusterNLayersChamberMinus32[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberMinus41",             cscRechitClusterNLayersChamberMinus41,             "cscRechitClusterNLayersChamberMinus41[nCscRechitClusters]/I");
    // tree_->Branch("cscRechitClusterNLayersChamberMinus42",             cscRechitClusterNLayersChamberMinus42,             "cscRechitClusterNLayersChamberMinus42[nCscRechitClusters]/I");



    tree_->Branch("cscRechitCluster_match_dtRechits_phi0p2",             cscRechitCluster_match_dtRechits_phi0p2,             "cscRechitCluster_match_dtRechits_phi0p2[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_dtRechits_0p4",             cscRechitCluster_match_dtRechits_0p4,             "cscRechitCluster_match_dtRechits_0p4[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_MB1_0p4",             cscRechitCluster_match_MB1_0p4,             "cscRechitCluster_match_MB1_0p4[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_dtSeg_0p4",             cscRechitCluster_match_dtSeg_0p4,             "cscRechitCluster_match_dtSeg_0p4[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_MB1Seg_0p4",             cscRechitCluster_match_MB1Seg_0p4,             "cscRechitCluster_match_MB1Seg_0p4[nCscRechitClusters]/I");

    tree_->Branch("cscRechitCluster_match_RB1_0p4",             cscRechitCluster_match_RB1_0p4,             "cscRechitCluster_match_RB1_0p4[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_RE12_0p4",             cscRechitCluster_match_RE12_0p4,             "cscRechitCluster_match_RE12_0p4[nCscRechitClusters]/I");
    tree_->Branch("cscRechitCluster_match_gLLP_deltaR",             cscRechitCluster_match_gLLP_deltaR,             "cscRechitCluster_match_gLLP_deltaR[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterJetVetoPt",             cscRechitClusterJetVetoPt,             "cscRechitClusterJetVetoPt[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoEta",             cscRechitClusterJetVetoEta,             "cscRechitClusterJetVetoEta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterJetVetoPhi",             cscRechitClusterJetVetoPhi,             "cscRechitClusterJetVetoPhi[nCscRechitClusters]/F");

    tree_->Branch("cscRechitClusterJetVetoE",             cscRechitClusterJetVetoE,             "cscRechitClusterJetVetoE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoPt",             cscRechitClusterMuonVetoPt,             "cscRechitClusterMuonVetoPt[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoE",             cscRechitClusterMuonVetoE,             "cscRechitClusterMuonVetoE[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoPhi",             cscRechitClusterMuonVetoPhi,             "cscRechitClusterMuonVetoPhi[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoEta",             cscRechitClusterMuonVetoEta,             "cscRechitClusterMuonVetoEta[nCscRechitClusters]/F");
    tree_->Branch("cscRechitClusterMuonVetoLooseId",             cscRechitClusterMuonVetoLooseId,             "cscRechitClusterMuonVetoLooseId[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterMuonVetoGlobal",             cscRechitClusterMuonVetoGlobal,             "cscRechitClusterMuonVetoGlobal[nCscRechitClusters]/O");
    tree_->Branch("cscRechitClusterMetEENoise_dPhi",             cscRechitClusterMetEENoise_dPhi,             "cscRechitClusterMetEENoise_dPhi[nCscRechitClusters]/F");



    tree_->Branch("nDtRechitClusters",             &nDtRechitClusters, "nDtRechitClusters/I");
    tree_->Branch("dtRechitClusterMaxDPhi",             dtRechitClusterMaxDPhi,             "dtRechitClusterMaxDPhi[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterMaxDPhi_index",             dtRechitClusterMaxDPhi_index,             "dtRechitClusterMaxDPhi_index[nDtRechitClusters]/I");

    tree_->Branch("dtRechitClusterNSegStation1",             dtRechitClusterNSegStation1,             "dtRechitClusterNSegStation1[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNSegStation2",             dtRechitClusterNSegStation2,             "dtRechitClusterNSegStation2[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNSegStation3",             dtRechitClusterNSegStation3,             "dtRechitClusterNSegStation3[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNSegStation4",             dtRechitClusterNSegStation4,             "dtRechitClusterNSegStation4[nDtRechitClusters]/I");

    tree_->Branch("dtRechitClusterNOppositeSegStation1",             dtRechitClusterNOppositeSegStation1,             "dtRechitClusterNOppositeSegStation1[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNOppositeSegStation2",             dtRechitClusterNOppositeSegStation2,             "dtRechitClusterNOppositeSegStation2[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNOppositeSegStation3",             dtRechitClusterNOppositeSegStation3,             "dtRechitClusterNOppositeSegStation3[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNOppositeSegStation4",             dtRechitClusterNOppositeSegStation4,             "dtRechitClusterNOppositeSegStation4[nDtRechitClusters]/I");

        // tree_->Branch("dtRechitCluster_match_gParticle_id",             dtRechitCluster_match_gParticle_id,             "dtRechitCluster_match_gParticle_id[nDtRechitClusters]/I");
        // tree_->Branch("dtRechitCluster_match_gParticle",             dtRechitCluster_match_gParticle,             "dtRechitCluster_match_gParticle[nDtRechitClusters]/O");
        // tree_->Branch("dtRechitCluster_match_gParticle_minDeltaR",             dtRechitCluster_match_gParticle_minDeltaR,             "dtRechitCluster_match_gParticle_minDeltaR[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitCluster_match_gParticle_index",             dtRechitCluster_match_gParticle_index,             "dtRechitCluster_match_gParticle_index[nDtRechitClusters]/I");
        // tree_->Branch("dtRechitCluster_match_gParticle_eta",             dtRechitCluster_match_gParticle_eta,             "dtRechitCluster_match_gParticle_eta[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitCluster_match_gParticle_phi",             dtRechitCluster_match_gParticle_phi,             "dtRechitCluster_match_gParticle_phi[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitCluster_match_gParticle_E",             dtRechitCluster_match_gParticle_E,             "dtRechitCluster_match_gParticle_E[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitCluster_match_gParticle_pt",             dtRechitCluster_match_gParticle_pt,             "dtRechitCluster_match_gParticle_pt[nDtRechitClusters]/F");
        // tree_->Branch("dtRechitCluster_match_gParticle_MotherId",             dtRechitCluster_match_gParticle_MotherId,             "dtRechitCluster_match_gParticle_MotherId[nDtRechitClusters]/I");

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

        tree_->Branch("dtRechitClusterOverlap",             dtRechitClusterOverlap,             "dtRechitClusterOverlap[nDtRechitClusters]/O");

        tree_->Branch("dtRechitClusterSize",             dtRechitClusterSize,             "dtRechitClusterSize[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNoiseHit",             dtRechitClusterNoiseHit,             "dtRechitClusterNoiseHit[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNoiseHitStation1",             dtRechitClusterNoiseHitStation1,             "dtRechitClusterNoiseHitStation1[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNoiseHitStation2",             dtRechitClusterNoiseHitStation2,             "dtRechitClusterNoiseHitStation2[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNoiseHitStation3",             dtRechitClusterNoiseHitStation3,             "dtRechitClusterNoiseHitStation3[nDtRechitClusters]/I");
        tree_->Branch("dtRechitClusterNoiseHitStation4",             dtRechitClusterNoiseHitStation4,             "dtRechitClusterNoiseHitStation4[nDtRechitClusters]/I");
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



          tree_->Branch("dtRechitCluster_match_RPCTime_dPhi0p5",             dtRechitCluster_match_RPCTime_dPhi0p5,             "dtRechitCluster_match_RPCTime_dPhi0p5[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_RPCTimeSpread_dPhi0p5",             dtRechitCluster_match_RPCTimeSpread_dPhi0p5,             "dtRechitCluster_match_RPCTimeSpread_dPhi0p5[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_RPCTime_dR0p4",             dtRechitCluster_match_RPCTime_dR0p4,             "dtRechitCluster_match_RPCTime_dR0p4[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_RPCTimeSpread_dR0p4",             dtRechitCluster_match_RPCTimeSpread_dR0p4,             "dtRechitCluster_match_RPCTimeSpread_dR0p4[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_RPChits_dR0p4",             dtRechitCluster_match_RPChits_dR0p4,             "dtRechitCluster_match_RPChits_dR0p4[nDtRechitClusters]/I");
          tree_->Branch("dtRechitCluster_match_RPCTime_sameStation_dR0p4",             dtRechitCluster_match_RPCTime_sameStation_dR0p4,             "dtRechitCluster_match_RPCTime_sameStation_dR0p4[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4",             dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4,             "dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[nDtRechitClusters]/F");
          tree_->Branch("dtRechitCluster_match_RPChits_sameStation_dR0p4",             dtRechitCluster_match_RPChits_sameStation_dR0p4,             "dtRechitCluster_match_RPChits_sameStation_dR0p4[nDtRechitClusters]/I");


          tree_->Branch("dtRechitCluster_match_gLLP_deltaR",             dtRechitCluster_match_gLLP_deltaR,             "dtRechitCluster_match_gLLP_deltaR[nDtRechitClusters]/F");

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

        tree_->Branch("dtRechitClusterJetVetoPt",             dtRechitClusterJetVetoPt,             "dtRechitClusterJetVetoPt[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoEta",             dtRechitClusterJetVetoEta,             "dtRechitClusterJetVetoEta[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterJetVetoPhi",             dtRechitClusterJetVetoPhi,             "dtRechitClusterJetVetoPhi[nDtRechitClusters]/F");

        tree_->Branch("dtRechitClusterJetVetoE",             dtRechitClusterJetVetoE,             "dtRechitClusterJetVetoE[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoPt",             dtRechitClusterMuonVetoPt,             "dtRechitClusterMuonVetoPt[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoE",             dtRechitClusterMuonVetoE,             "dtRechitClusterMuonVetoE[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoPhi",             dtRechitClusterMuonVetoPhi,             "dtRechitClusterMuonVetoPhi[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoEta",             dtRechitClusterMuonVetoEta,             "dtRechitClusterMuonVetoEta[nDtRechitClusters]/F");
        tree_->Branch("dtRechitClusterMuonVetoLooseId",             dtRechitClusterMuonVetoLooseId,             "dtRechitClusterMuonVetoLooseId[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterMuonVetoGlobal",             dtRechitClusterMuonVetoGlobal,             "dtRechitClusterMuonVetoGlobal[nDtRechitClusters]/O");
        tree_->Branch("dtRechitClusterMetEENoise_dPhi",             dtRechitClusterMetEENoise_dPhi,             "dtRechitClusterMetEENoise_dPhi[nDtRechitClusters]/F");

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



};
