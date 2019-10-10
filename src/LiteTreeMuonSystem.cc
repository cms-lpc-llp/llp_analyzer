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
  npv=0; npu=0; rho=-1; weight=-1.0;
  met=-1; metPhi=-1;
  gLepId = 0;
  gLepPt = 0.; gLepPhi = 0.; gLepEta = 0.; gLepE = 0.;
  //CSC
  nCsc = 0;
  nCscClusters = 0;
  nCscITClusters = 0;
  nCsc_JetVetoITCluster0p4 = 0;
  nCsc_JetMuonVetoITCluster0p4 = 0;
  nCsc_JetVetoITCluster0p4_Me1112Veto = 0;
  nCsc_JetMuonVetoITCluster0p4_Me1112Veto = 0;
  nCsc_JetVetoCluster0p4 = 0;
  nCsc_JetMuonVetoCluster0p4 = 0;
  nCsc_JetVetoCluster0p4_Me1112Veto = 0;
  nCsc_JetMuonVetoCluster0p4_Me1112Veto = 0;
  for( int i = 0; i < N_MAX_CSC; i++ )
  {
    // cscLabels[i] = -999;
    // cscITLabels[i] = -999;
    // cscStation[i] = -999;
    // cscChamber[i] = -999;
    //
    // cscPhi[i] = -999;   //[nCsc]
    // cscEta[i] = -999;   //[nCsc]
    // cscX[i] = -999;   //[nCsc]
    // cscY[i] = -999;   //[nCsc]
    // cscZ[i] = -999;   //[nCsc]
    // cscDirectionX[i] = -999;   //[nCsc]
    // cscDirectionY[i] = -999;   //[nCsc]
    // cscDirectionZ[i] = -999;   //[nCsc]
    // cscNRecHits[i] = -999;   //[nCsc]
    // cscNRecHits_flag[i] = -999;   //[nCsc]
    // cscT[i] = -999;   //[nCsc]
    // cscChi2[i] = -999;   //[nCsc]

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
    cscClusterJetVeto[i] = 0.0;
    cscClusterCaloJetVeto[i] = 0.0;
    cscClusterMuonVeto[i] = 0.0;
    cscClusterJetVetoE[i] = 0.0;
    cscClusterCaloJetVetoE[i] = 0.0;
    cscClusterMuonVetoE[i] = 0.0;
    cscClusterNChamber[i] = -999;
    cscClusterMaxChamberRatio[i] = -999.;
    cscClusterMaxChamber[i] = -999;
    cscClusterNStation[i] = -999;
    cscClusterMaxStationRatio[i] = -999.;
    cscClusterMaxStation[i] = -999;
    cscClusterMe11Ratio[i] = -999.;
    cscClusterMe12Ratio[i] = -999.;
    cscClusterVertexR[i] = 0.0;
    cscClusterVertexZ[i] = 0.0;
    cscClusterVertexDis[i] = 0.0;
    cscClusterVertexChi2[i] = 0.0;
    cscClusterVertexN[i] = 0;
    cscClusterVertexN1[i] = 0;
    cscClusterVertexN5[i] = 0;
    cscClusterVertexN10[i] = 0;
    cscClusterVertexN15[i] = 0;
    cscClusterVertexN20[i] = 0;

    cscITClusterSize[i] = -999;
    cscITClusterX[i] = -999.;
    cscITClusterY[i] = -999.;
    cscITClusterZ[i] = -999.;
    cscITClusterTime[i] = -999.;
    cscITClusterTimeRMS[i] = -999.;
    cscITClusterGenMuonDeltaR[i] = 999.;
    cscITClusterTimeSpread[i] = -999.;
    cscITClusterRadius[i] = -999.;
    cscITClusterMajorAxis[i] = -999.;
    cscITClusterMinorAxis[i] = -999.;
    cscITClusterXSpread[i] = -999.;
    cscITClusterYSpread[i] = -999.;
    cscITClusterZSpread[i] = -999.;
    cscITClusterEtaPhiSpread[i] = -999.;
    cscITClusterEtaSpread[i] = -999.;
    cscITClusterPhiSpread[i] = -999.;
    cscITClusterEta[i] = -999.;
    cscITClusterPhi[i] = -999.;
    cscITClusterJetVeto[i] = 0.0;
    cscITClusterCaloJetVeto[i] = 0.0;
    cscITClusterMuonVeto[i] = 0.0;
    cscITClusterJetVetoE[i] = 0.0;
    cscITClusterCaloJetVetoE[i] = 0.0;
    cscITClusterMuonVetoE[i] = 0.0;
    cscITClusterNChamber[i] = -999;
    cscITClusterMaxChamberRatio[i] = -999.;
    cscITClusterMaxChamber[i] = -999;
    cscITClusterNStation[i] = -999;
    cscITClusterMaxStationRatio[i] = -999.;
    cscITClusterMaxStation[i] = -999;
    cscITClusterMe11Ratio[i] = -999.;
    cscITClusterMe12Ratio[i] = -999.;
    cscITClusterVertexR[i] = 0.0;
    cscITClusterVertexZ[i] = 0.0;
    cscITClusterVertexDis[i] = 0.0;
    cscITClusterVertexChi2[i] = 0.0;
    cscITClusterVertexN[i] = 0;
    cscITClusterVertexN1[i] = 0;
    cscITClusterVertexN5[i] = 0;
    cscITClusterVertexN10[i] = 0;
    cscITClusterVertexN15[i] = 0;
    cscITClusterVertexN20[i] = 0;
    cscITCluster_match_cscCluster_index[i] = -999;
    cscITCluster_cscCluster_SizeRatio[i] = -999.;
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
  tree_->SetBranchAddress("npv",         &npv);
  tree_->SetBranchAddress("npu",         &npu);
  tree_->SetBranchAddress("weight",      &weight);
  tree_->SetBranchAddress("rho",         &rho);
  tree_->SetBranchAddress("met",         &met);
  tree_->SetBranchAddress("metPhi",      &metPhi);
  tree_->SetBranchAddress("gLepId",      &gLepId);
  tree_->SetBranchAddress("gLepPt",      &gLepPt);
  tree_->SetBranchAddress("gLepPhi",      &gLepPhi);
  tree_->SetBranchAddress("gLepE",      &gLepE);
  tree_->SetBranchAddress("gLepEta",      &gLepEta);

  //CSC
  tree_->SetBranchAddress("nCsc",             &nCsc);
  // tree_->SetBranchAddress("cscLabels",             cscLabels);
  // tree_->SetBranchAddress("cscITLabels",             cscITLabels);
  // tree_->SetBranchAddress("cscStation",             cscStation);
  // tree_->SetBranchAddress("cscChamber",             cscChamber);
  //
  // tree_->SetBranchAddress("cscPhi",           cscPhi);
  // tree_->SetBranchAddress("cscEta",           cscEta);
  // tree_->SetBranchAddress("cscX",             cscX);
  // tree_->SetBranchAddress("cscY",             cscY);
  // tree_->SetBranchAddress("cscZ",             cscZ);
  // tree_->SetBranchAddress("cscDirectionX",             cscDirectionX);
  // tree_->SetBranchAddress("cscDirectionY",             cscDirectionY);
  // tree_->SetBranchAddress("cscDirectionZ",             cscDirectionZ);
  // tree_->SetBranchAddress("cscNRecHits",      cscNRecHits);
  // tree_->SetBranchAddress("cscNRecHits_flag", cscNRecHits_flag);
  // tree_->SetBranchAddress("cscT",             cscT);
  // tree_->SetBranchAddress("cscChi2",          cscChi2);

  // CSC CLUSTER
  tree_->SetBranchAddress("nCscClusters",             &nCscClusters);
  tree_->SetBranchAddress("nCsc_JetVetoCluster0p4",             &nCsc_JetVetoCluster0p4);
  tree_->SetBranchAddress("nCsc_JetMuonVetoCluster0p4",             &nCsc_JetMuonVetoCluster0p4);
  tree_->SetBranchAddress("nCsc_JetVetoCluster0p4_Me1112Veto",             &nCsc_JetVetoCluster0p4_Me1112Veto);
  tree_->SetBranchAddress("nCsc_JetMuonVetoCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoCluster0p4_Me1112Veto);


  tree_->SetBranchAddress("cscClusterMe11Ratio",             &cscClusterMe11Ratio);
  tree_->SetBranchAddress("cscClusterMe12Ratio",             &cscClusterMe12Ratio);
  tree_->SetBranchAddress("cscClusterMaxStation",             &cscClusterMaxStation);
  tree_->SetBranchAddress("cscClusterMaxStationRatio",             &cscClusterMaxStationRatio);
  tree_->SetBranchAddress("cscClusterNStation",             &cscClusterNStation);
  tree_->SetBranchAddress("cscClusterMaxChamber",             &cscClusterMaxChamber);
  tree_->SetBranchAddress("cscClusterMaxChamberRatio",             &cscClusterMaxChamberRatio);
  tree_->SetBranchAddress("cscClusterNChamber",             &cscClusterNChamber);
  tree_->SetBranchAddress("cscClusterVertexR",             cscClusterVertexR);
  tree_->SetBranchAddress("cscClusterVertexZ",             cscClusterVertexZ);
  tree_->SetBranchAddress("cscClusterVertexDis",             cscClusterVertexDis);
  tree_->SetBranchAddress("cscClusterVertexChi2",             cscClusterVertexChi2);
  tree_->SetBranchAddress("cscClusterVertexN",             cscClusterVertexN);
  tree_->SetBranchAddress("cscClusterVertexN1",             cscClusterVertexN1);
  tree_->SetBranchAddress("cscClusterVertexN5",             cscClusterVertexN5);
  tree_->SetBranchAddress("cscClusterVertexN10",             cscClusterVertexN10);
  tree_->SetBranchAddress("cscClusterVertexN15",             cscClusterVertexN15);
  tree_->SetBranchAddress("cscClusterVertexN20",             cscClusterVertexN20);
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
  tree_->SetBranchAddress("cscClusterJetVeto",             cscClusterJetVeto);
  tree_->SetBranchAddress("cscClusterCaloJetVeto",             cscClusterCaloJetVeto);
  tree_->SetBranchAddress("cscClusterMuonVeto",             cscClusterMuonVeto);
  tree_->SetBranchAddress("cscClusterJetVetoE",             cscClusterJetVetoE);
  tree_->SetBranchAddress("cscClusterCaloJetVetoE",             cscClusterCaloJetVetoE);
  tree_->SetBranchAddress("cscClusterMuonVetoE",             cscClusterMuonVetoE);

  tree_->SetBranchAddress("cscClusterSize",             cscClusterSize);

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

  /*
  //Z-candidate
  tree_->SetBranchAddress("ZMass",       &ZMass);
  tree_->SetBranchAddress("ZPt",         &ZPt);
  tree_->SetBranchAddress("ZEta",        &ZEta);
  tree_->SetBranchAddress("ZPhi",        &ZPhi);
  tree_->SetBranchAddress("ZleptonIndex1", &ZleptonIndex1);
  tree_->SetBranchAddress("ZleptonIndex2", &ZleptonIndex2);
  tree_->SetBranchAddress("MT", &MT);
  */
  //jets
  tree_->SetBranchAddress("nJets",     &nJets);
  tree_->SetBranchAddress("jetE",      jetE);
  tree_->SetBranchAddress("jetPt",     jetPt);
  tree_->SetBranchAddress("jetEta",    jetEta);
  tree_->SetBranchAddress("jetPhi",    jetPhi);
  // tree_->SetBranchAddress("jetTime",   jetTime);
  tree_->SetBranchAddress("jetPassId", jetPassId);
  // tree_->SetBranchAddress("ecalNRechits",   ecalNRechits);
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
  tree_->Branch("category",    &category,   "category/i");    // dilepton category
  tree_->Branch("npv",         &npv,        "npv/i");         // number of primary vertices
  tree_->Branch("npu",         &npu,        "npu/i");         // number of in-time PU events (MC)
  tree_->Branch("weight",      &weight,     "weight/F");
  tree_->Branch("rho",         &rho,        "rho/F");
  tree_->Branch("met",         &met,        "met/F");         // MET
  tree_->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)
  tree_->Branch("gLepId",      &gLepId,     "gLepId/I");      // phi(MET)
  tree_->Branch("gLepPt",      &gLepPt,     "gLepPt/F");      // phi(MET)
  tree_->Branch("gLepE",      &gLepE,     "gLepE/F");      // phi(MET)
  tree_->Branch("gLepEta",      &gLepEta,     "gLepEta/F");      // phi(MET)
  tree_->Branch("gLepPhi",      &gLepPhi,     "gLepPhi/F");      // phi(MET)

  //CSC
  tree_->Branch("nCsc",             &nCsc, "nCsc/I");
  // tree_->Branch("cscITLabels",             cscITLabels,             "cscITLabels[nCsc]/I");
  //
  // tree_->Branch("cscLabels",             cscLabels,             "cscLabels[nCsc]/I");
  // tree_->Branch("cscStation",             cscStation,             "cscStation[nCsc]/I");
  // tree_->Branch("cscChamber",             cscChamber,             "cscChamber[nCsc]/I");
  //
  // tree_->Branch("cscPhi",           cscPhi,           "cscPhi[nCsc]/F");
  // tree_->Branch("cscEta",           cscEta,           "cscEta[nCsc]/F");
  // tree_->Branch("cscX",             cscX,             "cscX[nCsc]/F");
  // tree_->Branch("cscY",             cscY,             "cscY[nCsc]/F");
  // tree_->Branch("cscZ",             cscZ,             "cscZ[nCsc]/F");
  // tree_->Branch("cscDirectionX",             cscDirectionX,             "cscDirectionX[nCsc]/F");
  // tree_->Branch("cscDirectionY",             cscDirectionY,             "cscDirectionY[nCsc]/F");
  // tree_->Branch("cscDirectionZ",             cscDirectionZ,             "cscDirectionZ[nCsc]/F");
  // tree_->Branch("cscNRecHits",      cscNRecHits,      "cscNRecHits[nCsc]/F");
  // tree_->Branch("cscNRecHits_flag", cscNRecHits_flag, "cscNRecHits_flag[nCsc]/F");
  // tree_->Branch("cscT",             cscT,             "cscT[nCsc]/F");
  // tree_->Branch("cscChi2",          cscChi2,          "cscChi2[nCsc]/F");

  // all csc clusters
  tree_->Branch("nCscClusters",             &nCscClusters, "nCscClusters/I");
  tree_->Branch("nCsc_JetVetoCluster0p4",             &nCsc_JetVetoCluster0p4, "nCsc_JetVetoCluster0p4/I");
  tree_->Branch("nCsc_JetMuonVetoCluster0p4",             &nCsc_JetMuonVetoCluster0p4, "nCsc_JetMuonVetoCluster0p4/I");
  tree_->Branch("nCsc_JetVetoCluster0p4_Me1112Veto",             &nCsc_JetVetoCluster0p4_Me1112Veto, "nCsc_JetVetoCluster0p4_Me1112Veto/I");
  tree_->Branch("nCsc_JetMuonVetoCluster0p4_Me1112Veto",             &nCsc_JetMuonVetoCluster0p4_Me1112Veto, "nCsc_JetMuonVetoCluster0p4_Me1112Veto/I");
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
  tree_->Branch("cscClusterJetVeto",             cscClusterJetVeto,             "cscClusterJetVeto[nCscClusters]/F");
  tree_->Branch("cscClusterMuonVeto",             cscClusterMuonVeto,             "cscClusterMuonVeto[nCscClusters]/F");
  tree_->Branch("cscClusterCaloJetVeto",             cscClusterCaloJetVeto,             "cscClusterCaloJetVeto[nCscClusters]/F");
  tree_->Branch("cscClusterJetVetoE",             cscClusterJetVetoE,             "cscClusterJetVetoE[nCscClusters]/F");
  tree_->Branch("cscClusterMuonVetoE",             cscClusterMuonVetoE,             "cscClusterMuonVetoE[nCscClusters]/F");
  tree_->Branch("cscClusterCaloJetVetoE",             cscClusterCaloJetVetoE,             "cscClusterCaloJetVetoE[nCscClusters]/F");
  tree_->Branch("cscClusterSize",             cscClusterSize,             "cscClusterSize[nCscClusters]/I");
  tree_->Branch("cscClusterMe11Ratio",             cscClusterMe11Ratio,             "cscClusterMe11Ratio[nCscClusters]/F");
  tree_->Branch("cscClusterMe12Ratio",             cscClusterMe12Ratio,             "cscClusterMe12Ratio[nCscClusters]/F");

  tree_->Branch("cscClusterNStation",             cscClusterNStation,             "cscClusterNStation[nCscClusters]/I");
  tree_->Branch("cscClusterMaxStation",             cscClusterMaxStation,             "cscClusterMaxStation[nCscClusters]/I");
  tree_->Branch("cscClusterMaxStationRatio",             cscClusterMaxStationRatio,             "cscClusterMaxStationRatio[nCscClusters]/F");

  tree_->Branch("cscClusterNChamber",             cscClusterNChamber,             "cscClusterNChamber[nCscClusters]/I");
  tree_->Branch("cscClusterMaxChamber",             cscClusterMaxChamber,             "cscClusterMaxChamber[nCscClusters]/I");
  tree_->Branch("cscClusterMaxChamberRatio",             cscClusterMaxChamberRatio,             "cscClusterMaxChamberRatio[nCscClusters]/F");
  tree_->Branch("cscClusterVertexR",             cscClusterVertexR,             "cscClusterVertexR[nCscClusters]/F");
  tree_->Branch("cscClusterVertexZ",             cscClusterVertexZ,             "cscClusterVertexZ[nCscClusters]/F");
  tree_->Branch("cscClusterVertexDis",             cscClusterVertexDis,             "cscClusterVertexDis[nCscClusters]/F");
  tree_->Branch("cscClusterVertexChi2",             cscClusterVertexChi2,             "cscClusterVertexChi2[nCscClusters]/F");
  tree_->Branch("cscClusterVertexN1",             cscClusterVertexN1,             "cscClusterVertexN1[nCscClusters]/I");
  tree_->Branch("cscClusterVertexN5",             cscClusterVertexN5,             "cscClusterVertexN5[nCscClusters]/I");
  tree_->Branch("cscClusterVertexN10",             cscClusterVertexN10,             "cscClusterVertexN10[nCscClusters]/I");
  tree_->Branch("cscClusterVertexN15",             cscClusterVertexN15,             "cscClusterVertexN15[nCscClusters]/I");
  tree_->Branch("cscClusterVertexN20",             cscClusterVertexN20,             "cscClusterVertexN20[nCscClusters]/I");
  tree_->Branch("cscClusterVertexN",             cscClusterVertexN,             "cscClusterVertexN[nCscClusters]/I");

  // INTIME CSC cluster
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
  tree_->Branch("lepPassVetoId", lepPassVetoId, "lepPassVetoId[nLeptons]/O");

  // tree_->Branch("lepLoosePassId", lepLoosePassId, "lepLoosePassId[nLeptons]/O");
  // tree_->Branch("lepMediumPassId", lepMediumPassId, "lepMediumPassId[nLeptons]/O");
  // tree_->Branch("lepTightPassId", lepTightPassId, "lepTightPassId[nLeptons]/O");
  /*
  //Z-candidate
  tree_->Branch("MT",      &MT,  "MT/F");
  tree_->Branch("ZMass",      &ZMass,  "ZMass/F");
  tree_->Branch("ZPt",        &ZPt,    "ZPt/F");
  tree_->Branch("ZEta",       &ZEta,   "ZEta/F");
  tree_->Branch("ZPhi",       &ZPhi,   "ZPhi/F");
  tree_->Branch("ZleptonIndex1", &ZleptonIndex1, "ZleptonIndex1/I");
  tree_->Branch("ZleptonIndex2", &ZleptonIndex2, "ZleptonIndex2/I");
  */
  //jets
  tree_->Branch("nJets",     &nJets,    "nJets/I");
  tree_->Branch("jetE",      jetE,      "jetE[nJets]/F");
  tree_->Branch("jetPt",     jetPt,     "jetPt[nJets]/F");
  tree_->Branch("jetEta",    jetEta,    "jetEta[nJets]/F");
  tree_->Branch("jetPhi",    jetPhi,    "jetPhi[nJets]/F");
  // tree_->Branch("jetTime",   jetTime,   "jetTime[nJets]/F");
  tree_->Branch("jetPassId", jetPassId, "jetPassId[nJets]/O");
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
