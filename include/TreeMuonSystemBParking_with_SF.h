// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef TreeMuonSystemBParking_H
#define TreeMuonSystemBParking_H

#define N_MAX_LEPTONS 100
#define N_MAX_TRACKS 2000
#define N_MAX_JETS 100
#define N_MAX_CSC 200
#define N_MAX_CSCRECHITS 5000
#define N_MAX_DTRECHITS 20000
#define NTriggersMAX 1201 // Number of trigger in the .dat file
#define N_CSC_CUT 20
#define JET_PT_CUT 10
#define MUON_PT_CUT 20
#define N_MAX_GPARTICLES 5000
#define MAX_MuonHLTFilters 100

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

class TreeMuonSystemBParking
{

public:
  TreeMuonSystemBParking();
  ~TreeMuonSystemBParking();
  // TreeMuonSystemBParking::TreeMuonSystemBParking()
  // {
  //   InitVariables();
  // };
  // TreeMuonSystemBParking::~TreeMuonSystemBParking()
  // {
  //   if (f_) f_->Close();
  // };
  TTree *tree_;
  TFile *f_;

  UInt_t  runNum, lumiSec, evtNum, MC_condition,npv;

  float rho, met, metPhi, metEENoise, metPhiEENoise;


  bool Flag2_all;

  //gLLP
  float gLLP_csc;
  float gLLP_dt;
  float gLLP_beta;
  float gLLP_eta;
  float gLLP_phi;
  float gLLP_e;
  float gLLP_pt;
  float gLLP_ctau;
  float gLLP_decay_vertex_r;
  float gLLP_decay_vertex_x;
  float gLLP_decay_vertex_y;
  float gLLP_decay_vertex_z;

  int nGenParticles;
  float gParticleE[N_MAX_GPARTICLES];
  float gParticlePt[N_MAX_GPARTICLES];
  float gParticleEta[N_MAX_GPARTICLES];
  float gParticlePhi[N_MAX_GPARTICLES];
  int gParticleMotherId[N_MAX_GPARTICLES];
  int gParticleMotherIndex[N_MAX_GPARTICLES];
  int gParticleId[N_MAX_GPARTICLES];
    
    int nTracks;
    float track_Pt[N_MAX_TRACKS];
    float track_Eta[N_MAX_TRACKS];
    float track_Phi[N_MAX_TRACKS];
    
    /*
    float track_pt_sum_0p2_DT;
    float track_pt_sum_0p3_DT;
    float track_pt_sum_0p4_DT;
    float track_pt_sum_0p5_DT;
    float leading_track_pt_0p2_DT;
    float leading_track_pt_0p3_DT;
    float leading_track_pt_0p4_DT;
    float leading_track_pt_0p5_DT;
    int matched_track_size_0p2_DT;
    int matched_track_size_0p3_DT;
    int matched_track_size_0p4_DT;
    int matched_track_size_0p5_DT;

    float track_pt_sum_0p2_CSC;
    float track_pt_sum_0p3_CSC;
    float track_pt_sum_0p4_CSC;
    float track_pt_sum_0p5_CSC;
    float leading_track_pt_0p2_CSC;
    float leading_track_pt_0p3_CSC;
    float leading_track_pt_0p4_CSC;
    float leading_track_pt_0p5_CSC;
    int matched_track_size_0p2_CSC;
    int matched_track_size_0p3_CSC;
    int matched_track_size_0p4_CSC;
    int matched_track_size_0p5_CSC;
    */
    
    int nLeptons;
    float lepE[N_MAX_LEPTONS];
    float lepPt[N_MAX_LEPTONS];
    float lepEta[N_MAX_LEPTONS];
    float lepPhi[N_MAX_LEPTONS];
    int  lepPdgId[N_MAX_LEPTONS];
    float lepDZ[N_MAX_LEPTONS];
    float lepDXY[N_MAX_LEPTONS];
    float lepDXYErr[N_MAX_LEPTONS];
    float lepSF[N_MAX_LEPTONS];
    bool lepLooseId[N_MAX_LEPTONS];
    bool lepTightId[N_MAX_LEPTONS];
    bool lepPassLooseIso[N_MAX_LEPTONS];
    bool lepPassTightIso[N_MAX_LEPTONS];
    bool lepPassVTightIso[N_MAX_LEPTONS];
    bool lepPassVVTightIso[N_MAX_LEPTONS];
    unsigned int lepMuonType[N_MAX_LEPTONS];//only assigned for muons
    unsigned int lepMuonQuality[N_MAX_LEPTONS];//only assigned for muons
    //bool lep_passSingleMuTagFilter[N_MAX_LEPTONS];//only assigned for muons
    bool  lepMuon_passHLTFilter[N_MAX_LEPTONS][MAX_MuonHLTFilters];

    int nJets;
    float jetE[N_MAX_JETS];
    float jetPt[N_MAX_JETS];
    float jetEta[N_MAX_JETS];
    float jetPhi[N_MAX_JETS];
    bool jetTightPassId[N_MAX_JETS];
    bool HLTDecision[NTriggersMAX];



  //csc
  int           nCscRechits;

  int           nCscRings;

  // int           nCscRechitsChamberPlus11;
  // int           nCscRechitsChamberPlus12;
  // int           nCscRechitsChamberPlus13;
  // int           nCscRechitsChamberPlus21;
  // int           nCscRechitsChamberPlus22;
  // int           nCscRechitsChamberPlus31;
  // int           nCscRechitsChamberPlus32;
  // int           nCscRechitsChamberPlus41;
  // int           nCscRechitsChamberPlus42;
  //
  // int           nCscRechitsChamberMinus11;
  // int           nCscRechitsChamberMinus12;
  // int           nCscRechitsChamberMinus13;
  // int           nCscRechitsChamberMinus21;
  // int           nCscRechitsChamberMinus22;
  // int           nCscRechitsChamberMinus31;
  // int           nCscRechitsChamberMinus32;
  // int           nCscRechitsChamberMinus41;
  // int           nCscRechitsChamberMinus42;


  int           nDtRechits;
  int           nDtRings;
  int           nDtStations25;
  int           nDtWheels25;
  // int           nDTRechitsStation1;
  // int           nDTRechitsStation2;
  // int           nDTRechitsStation3;
  // int           nDTRechitsStation4;
  // int           nDTRechitsWheelMinus2;
  // int           nDTRechitsWheelMinus1;
  // int           nDTRechitsWheel0;
  // int           nDTRechitsWheelPlus1;
  // int           nDTRechitsWheelPlus2;
  // int           nDTRechitsChamberMinus12;
  // int           nDTRechitsChamberMinus11;
  // int           nDTRechitsChamber10;
  // int           nDTRechitsChamberPlus11;
  // int           nDTRechitsChamberPlus12;
  // int           nDTRechitsChamberMinus22;
  // int           nDTRechitsChamberMinus21;
  // int           nDTRechitsChamber20;
  // int           nDTRechitsChamberPlus21;
  // int           nDTRechitsChamberPlus22;
  // int           nDTRechitsChamberMinus32;
  // int           nDTRechitsChamberMinus31;
  // int           nDTRechitsChamber30;
  // int           nDTRechitsChamberPlus31;
  // int           nDTRechitsChamberPlus32;
  // int           nDTRechitsChamberMinus42;
  // int           nDTRechitsChamberMinus41;
  // int           nDTRechitsChamber40;
  // int           nDTRechitsChamberPlus41;
  // int           nDTRechitsChamberPlus42;


  int           nDtRechitClusters;
  float           dtRechitClusterMaxDPhi[N_MAX_CSC];
  bool           dtRechitClusterOverlap[N_MAX_CSC];
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

  float cscRechitClusterMatchedTrackSumPt_0p2[N_MAX_CSC];
  float cscRechitClusterMatchedTrackLeadPt_0p2[N_MAX_CSC];
  int   cscRechitClusterMatchedTrackSize_0p2[N_MAX_CSC];

  float cscRechitClusterMatchedTrackSumPt_0p3[N_MAX_CSC];
  float cscRechitClusterMatchedTrackLeadPt_0p3[N_MAX_CSC];
  int   cscRechitClusterMatchedTrackSize_0p3[N_MAX_CSC];

  float cscRechitClusterMatchedTrackSumPt_0p4[N_MAX_CSC];
  float cscRechitClusterMatchedTrackLeadPt_0p4[N_MAX_CSC];
  int   cscRechitClusterMatchedTrackSize_0p4[N_MAX_CSC];

  float cscRechitClusterMatchedTrackSumPt_0p5[N_MAX_CSC];
  float cscRechitClusterMatchedTrackLeadPt_0p5[N_MAX_CSC];
  int   cscRechitClusterMatchedTrackSize_0p5[N_MAX_CSC];

  float dtRechitClusterMatchedTrackSumPt_0p2[N_MAX_CSC];
  float dtRechitClusterMatchedTrackLeadPt_0p2[N_MAX_CSC];
  int   dtRechitClusterMatchedTrackSize_0p2[N_MAX_CSC];

  float dtRechitClusterMatchedTrackSumPt_0p3[N_MAX_CSC];
  float dtRechitClusterMatchedTrackLeadPt_0p3[N_MAX_CSC];
  int   dtRechitClusterMatchedTrackSize_0p3[N_MAX_CSC];

  float dtRechitClusterMatchedTrackSumPt_0p4[N_MAX_CSC];
  float dtRechitClusterMatchedTrackLeadPt_0p4[N_MAX_CSC];
  int   dtRechitClusterMatchedTrackSize_0p4[N_MAX_CSC];

  float dtRechitClusterMatchedTrackSumPt_0p5[N_MAX_CSC];
  float dtRechitClusterMatchedTrackLeadPt_0p5[N_MAX_CSC];
  int   dtRechitClusterMatchedTrackSize_0p5[N_MAX_CSC];

  float cscRechitClusterMatchedTrackSumPt_trk_pos_0p2[N_MAX_CSC];
  float cscRechitClusterMatchedTrackLeadPt_trk_pos_0p2[N_MAX_CSC];
  int   cscRechitClusterMatchedTrackSize_trk_pos_0p2[N_MAX_CSC];

  float cscRechitClusterMatchedTrackSumPt_trk_pos_0p3[N_MAX_CSC];
  float cscRechitClusterMatchedTrackLeadPt_trk_pos_0p3[N_MAX_CSC];
  int   cscRechitClusterMatchedTrackSize_trk_pos_0p3[N_MAX_CSC];

  float cscRechitClusterMatchedTrackSumPt_trk_pos_0p4[N_MAX_CSC];
  float cscRechitClusterMatchedTrackLeadPt_trk_pos_0p4[N_MAX_CSC];
  int   cscRechitClusterMatchedTrackSize_trk_pos_0p4[N_MAX_CSC];

  float cscRechitClusterMatchedTrackSumPt_trk_pos_0p5[N_MAX_CSC];
  float cscRechitClusterMatchedTrackLeadPt_trk_pos_0p5[N_MAX_CSC];
  int   cscRechitClusterMatchedTrackSize_trk_pos_0p5[N_MAX_CSC];

  float dtRechitClusterMatchedTrackSumPt_trk_pos_0p2[N_MAX_CSC];
  float dtRechitClusterMatchedTrackLeadPt_trk_pos_0p2[N_MAX_CSC];
  int   dtRechitClusterMatchedTrackSize_trk_pos_0p2[N_MAX_CSC];

  float dtRechitClusterMatchedTrackSumPt_trk_pos_0p3[N_MAX_CSC];
  float dtRechitClusterMatchedTrackLeadPt_trk_pos_0p3[N_MAX_CSC];
  int   dtRechitClusterMatchedTrackSize_trk_pos_0p3[N_MAX_CSC];

  float dtRechitClusterMatchedTrackSumPt_trk_pos_0p4[N_MAX_CSC];
  float dtRechitClusterMatchedTrackLeadPt_trk_pos_0p4[N_MAX_CSC];
  int   dtRechitClusterMatchedTrackSize_trk_pos_0p4[N_MAX_CSC];

  float dtRechitClusterMatchedTrackSumPt_trk_pos_0p5[N_MAX_CSC];
  float dtRechitClusterMatchedTrackLeadPt_trk_pos_0p5[N_MAX_CSC];
  int   dtRechitClusterMatchedTrackSize_trk_pos_0p5[N_MAX_CSC];

  int           dtRechitCluster_match_MB1hits_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_MB1hits_0p5[N_MAX_CSC];
  int           dtRechitCluster_match_MB1hits_cosmics_plus[N_MAX_CSC];
  int           dtRechitCluster_match_MB1hits_cosmics_minus[N_MAX_CSC];
  int           dtRechitCluster_match_MB1Seg_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_MB1Seg_0p5[N_MAX_CSC];

  int           dtRechitCluster_match_RPCBx_dPhi0p5[N_MAX_CSC];
  int           dtRechitCluster_match_RB1_0p4[N_MAX_CSC];
  int           dtRechitCluster_match_RB1_dPhi0p5[N_MAX_CSC];
  //

    float         dtRechitClusterJetVetoEta[N_MAX_CSC];
    float         dtRechitClusterJetVetoPhi[N_MAX_CSC];
    float         dtRechitClusterJetVetoPt[N_MAX_CSC];
    float         dtRechitClusterJetVetoE[N_MAX_CSC];

    float         dtRechitClusterGenMuonVetoPt[N_MAX_CSC];
    float         dtRechitClusterGenMuonVetoPt_dR0p8[N_MAX_CSC];
    float         dtRechitClusterMuonVetoPt[N_MAX_CSC];
    float         dtRechitClusterMuonVetoE[N_MAX_CSC];
    float         dtRechitClusterMuonVetoPhi[N_MAX_CSC];
    float         dtRechitClusterMuonVetoEta[N_MAX_CSC];
    bool          dtRechitClusterMuonVetoLooseId[N_MAX_CSC];
    bool          dtRechitClusterMuonVetoGlobal[N_MAX_CSC];
float           dtRechitClusterMetEENoise_dPhi[N_MAX_CSC];
float         dtRechitCluster_match_gLLP_deltaR[N_MAX_CSC];
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

  float         dtRechitClusterX[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterY[N_MAX_CSC];   //[nCsc]
  float         dtRechitClusterZ[N_MAX_CSC];   //[nCsc]

  int         dtRechitClusterWheel[N_MAX_CSC];


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
  float           cscRechitCluster_match_cluster_dR[N_MAX_CSC];
  int           cscRechitCluster_match_cluster_index[N_MAX_CSC];


  float         cscRechitClusterX[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterY[N_MAX_CSC];   //[nCsc]
  float         cscRechitClusterZ[N_MAX_CSC];   //[nCsc]
  // float         cscRechitClusterTime[N_MAX_CSC];   //[nCsc]
  // float         cscRechitClusterTimeTotal[N_MAX_CSC];
  float         cscRechitClusterTimeWeighted[N_MAX_CSC];

  // float         cscRechitClusterTimeSpread[N_MAX_CSC];
  // float         cscRechitClusterTimeSpreadWeighted[N_MAX_CSC];
  float         cscRechitClusterTimeSpreadWeightedAll[N_MAX_CSC];

  float         cscRechitClusterGenMuonDeltaR[N_MAX_CSC];


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

  // int           cscRechitClusterNLayersChamberPlus11[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberPlus12[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberPlus13[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberPlus21[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberPlus22[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberPlus31[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberPlus32[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberPlus41[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberPlus42[N_MAX_CSC];
  //
  // int           cscRechitClusterNLayersChamberMinus11[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberMinus12[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberMinus13[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberMinus21[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberMinus22[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberMinus31[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberMinus32[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberMinus41[N_MAX_CSC];
  // int           cscRechitClusterNLayersChamberMinus42[N_MAX_CSC];




  // float         cscRechitClusterJetVetoEta[N_MAX_CSC];
  // float         cscRechitClusterJetVetoPhi[N_MAX_CSC];
  float         cscRechitClusterJetVetoPt[N_MAX_CSC];
  // float         cscRechitClusterJetVetoE[N_MAX_CSC];

  float         cscRechitClusterGenMuonVetoPt[N_MAX_CSC];
  float         cscRechitClusterGenMuonVetoPt_dR0p8[N_MAX_CSC];

  float         cscRechitClusterMuonVetoPt[N_MAX_CSC];
  // float         cscRechitClusterMuonVetoE[N_MAX_CSC];
  // float         cscRechitClusterMuonVetoPhi[N_MAX_CSC];
  // float         cscRechitClusterMuonVetoEta[N_MAX_CSC];
  // bool          cscRechitClusterMuonVetoLooseId[N_MAX_CSC];
  // bool          cscRechitClusterMuonVetoGlobal[N_MAX_CSC];
  float           cscRechitClusterMetEENoise_dPhi[N_MAX_CSC];
float           cscRechitCluster_match_gLLP_deltaR[N_MAX_CSC];

  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();

};
#endif
