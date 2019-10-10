// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef LiteTreeMuonSystem_H
#define LiteTreeMuonSystem_H

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 2000
#define NTriggersMAX 601 //Number of trigger in the .dat file
#define N_CSC_CUT 20
#define JET_PT_CUT 10
#define MUON_PT_CUT 20

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

#include "RazorAnalyzer.h"

#include "RazorHelper.h"

class LiteTreeMuonSystem
{

public:
  LiteTreeMuonSystem();
  ~LiteTreeMuonSystem();
  // LiteTreeMuonSystem::LiteTreeMuonSystem()
  // {
  //   InitVariables();
  // };
  // LiteTreeMuonSystem::~LiteTreeMuonSystem()
  // {
  //   if (f_) f_->Close();
  // };
  TTree *tree_;
  TFile *f_;

  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  category;
  UInt_t  npv, npu;
  float rho, weight;
  float met, metPhi;
  int gLepId;
  float gLepPt, gLepPhi, gLepEta, gLepE;


  //csc
  int           nCsc;
  // int           cscLabels[N_MAX_CSC];
  // int           cscITLabels[N_MAX_CSC];
  // int           cscStation[N_MAX_CSC];
  // int           cscChamber[N_MAX_CSC];
  //
  // float         cscPhi[N_MAX_CSC];   //[nCsc]
  // float         cscEta[N_MAX_CSC];   //[nCsc]
  // float         cscX[N_MAX_CSC];   //[nCsc]
  // float         cscY[N_MAX_CSC];   //[nCsc]
  // float         cscZ[N_MAX_CSC];   //[nCsc]
  // float         cscDirectionX[N_MAX_CSC];   //[nCsc]
  // float         cscDirectionY[N_MAX_CSC];   //[nCsc]
  // float         cscDirectionZ[N_MAX_CSC];   //[nCsc]
  // float         cscNRecHits[N_MAX_CSC];   //[nCsc]
  // float         cscNRecHits_flag[N_MAX_CSC];   //[nCsc]
  // float         cscNRecHits_jetveto0p4[N_MAX_CSC];   //[nCsc]
  // float         cscNRecHits_jetveto0p8[N_MAX_CSC];   //[nCsc]
  // float         cscT[N_MAX_CSC];   //[nCsc]
  // float         cscChi2[N_MAX_CSC];   //[nCsc]

  int           nCscClusters;
  float         cscClusterJetVeto[N_MAX_CSC];   //[nCsc]
  float         cscClusterCaloJetVeto[N_MAX_CSC];
  float         cscClusterMuonVeto[N_MAX_CSC];   //[nCsc]
  float         cscClusterJetVetoE[N_MAX_CSC];   //[nCsc]
  float         cscClusterCaloJetVetoE[N_MAX_CSC];
  float         cscClusterMuonVetoE[N_MAX_CSC];   //[nCsc]

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
  float         cscClusterMaxChamberRatio[N_MAX_CSC];   //[nCsc]
  int           cscClusterMaxChamber[N_MAX_CSC];   //[nCsc]
  int           cscClusterNChamber[N_MAX_CSC];
  int           nCsc_JetVetoCluster0p4;
  int           nCsc_JetMuonVetoCluster0p4;
  int           nCsc_JetVetoCluster0p4_Me1112Veto;
  int           nCsc_JetMuonVetoCluster0p4_Me1112Veto;
  //csc intime cluster
  int           nCscITClusters;
  float         cscITClusterJetVeto[N_MAX_CSC];   //[nCsc]
  float         cscITClusterCaloJetVeto[N_MAX_CSC];
  float         cscITClusterMuonVeto[N_MAX_CSC];   //[nCsc]
  float         cscITClusterJetVetoE[N_MAX_CSC];   //[nCsc]
  float         cscITClusterCaloJetVetoE[N_MAX_CSC];
  float         cscITClusterMuonVetoE[N_MAX_CSC];   //[nCsc]
  float         cscITClusterX[N_MAX_CSC];   //[nCsc]
  float         cscITClusterY[N_MAX_CSC];   //[nCsc]
  float         cscITClusterZ[N_MAX_CSC];   //[nCsc]
  float         cscITClusterRadius[N_MAX_CSC];   //[nCsc]
  float         cscITClusterTime[N_MAX_CSC];   //[nCsc]
  float         cscITClusterTimeSpread[N_MAX_CSC];
  float         cscITClusterTimeRMS[N_MAX_CSC];
  float         cscITClusterGenMuonDeltaR[N_MAX_CSC];
  float         cscITClusterMajorAxis[N_MAX_CSC];
  float         cscITClusterMinorAxis[N_MAX_CSC];
  float         cscITClusterXSpread[N_MAX_CSC];   //[nCsc]
  float         cscITClusterYSpread[N_MAX_CSC];   //[nCsc]
  float         cscITClusterZSpread[N_MAX_CSC];   //[nCsc]
  float         cscITClusterEtaPhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscITClusterEtaSpread[N_MAX_CSC];   //[nCsc]
  float         cscITClusterPhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscITClusterEta[N_MAX_CSC];   //[nCsc]
  float         cscITClusterPhi[N_MAX_CSC];   //[nCsc]
  int           cscITClusterSize[N_MAX_CSC];
  float         cscITClusterMe11Ratio[N_MAX_CSC];
  float         cscITClusterMe12Ratio[N_MAX_CSC];

  float         cscITClusterVertexR[N_MAX_CSC];   //[nCsc]
  float         cscITClusterVertexZ[N_MAX_CSC];   //[nCsc]
  int           cscITClusterVertexN[N_MAX_CSC];   //[nCsc]
  int           cscITClusterVertexN1[N_MAX_CSC];   //[nCsc]
  int           cscITClusterVertexN5[N_MAX_CSC];   //[nCsc]
  int           cscITClusterVertexN10[N_MAX_CSC];   //[nCsc]
  int           cscITClusterVertexN15[N_MAX_CSC];   //[nCsc]
  int           cscITClusterVertexN20[N_MAX_CSC];   //[nCsc]
  float         cscITClusterVertexChi2[N_MAX_CSC];   //[nCsc]
  float         cscITClusterVertexDis[N_MAX_CSC];   //[nCsc]
  float         cscITClusterMaxStationRatio[N_MAX_CSC];   //[nCsc]
  int           cscITClusterMaxStation[N_MAX_CSC];   //[nCsc]
  int           cscITClusterNStation[N_MAX_CSC];
  float         cscITClusterMaxChamberRatio[N_MAX_CSC];   //[nCsc]
  int           cscITClusterMaxChamber[N_MAX_CSC];   //[nCsc]
  int           cscITClusterNChamber[N_MAX_CSC];
  int           cscITCluster_match_cscCluster_index[N_MAX_CSC];
  float         cscITCluster_cscCluster_SizeRatio[N_MAX_CSC];
  int           nCsc_JetVetoITCluster0p4;
  int           nCsc_JetMuonVetoITCluster0p4;
  int           nCsc_JetVetoITCluster0p4_Me1112Veto;
  int           nCsc_JetMuonVetoITCluster0p4_Me1112Veto;
  //gLLP
  float gLLP_eta[2];
  float gLLP_phi[2];
  float gLLP_csc[2];
  float gLLP_beta[2];
  float gLLP_ctau[2];
  float gLLP_decay_vertex_r[2];
  float gLLP_decay_vertex_x[2];
  float gLLP_decay_vertex_y[2];
  float gLLP_decay_vertex_z[2];

  //leptons
  int nLeptons;
  float lepE[N_MAX_LEPTONS];
  float lepPt[N_MAX_LEPTONS];
  float lepEta[N_MAX_LEPTONS];
  float lepPhi[N_MAX_LEPTONS];
  int  lepPdgId[N_MAX_LEPTONS];
  float lepDZ[N_MAX_LEPTONS];
  // bool lepLoosePassId[N_MAX_LEPTONS];
  // bool lepMediumPassId[N_MAX_LEPTONS];
  // bool lepTightPassId[N_MAX_LEPTONS];
  bool lepPassVetoId[N_MAX_LEPTONS];

  bool lepPassId[N_MAX_LEPTONS];

  //Z-candidate
  float MT;
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
  float ecalNRechits[N_MAX_JETS];
  float ecalRechitE[N_MAX_JETS];
  float jetChargedEMEnergyFraction[N_MAX_JETS];
  float jetNeutralEMEnergyFraction[N_MAX_JETS];
  float jetChargedHadronEnergyFraction[N_MAX_JETS];
  float jetNeutralHadronEnergyFraction[N_MAX_JETS];
  // bool jetLoosePassId[N_MAX_JETS];
  bool jetPassId[N_MAX_JETS];
  // bool jetTightPassId[N_MAX_JETS];
  bool HLTDecision[NTriggersMAX];
  UInt_t wzevtNum,trig, trig_lepId, trig_lepId_dijet; //number of events that pass each criteria




  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();

};
#endif
