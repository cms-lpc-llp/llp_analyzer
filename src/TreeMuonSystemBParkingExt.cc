#include "RazorHelper.h"
#include "TreeMuonSystemBParkingExt.h"
#include "assert.h"
#include "TTree.h"
#include "DBSCAN.h"

// Constructor
TreeMuonSystemBParkingExt::TreeMuonSystemBParkingExt()
{
  InitVariables();
};
TreeMuonSystemBParkingExt::~TreeMuonSystemBParkingExt()
{
  if (f_) f_->Close();
};
void TreeMuonSystemBParkingExt::InitVariables()
{
  runNum=0; lumiSec=0; evtNum=0;  MC_condition = 0;npv=0; npu=0; pileupWeight=0.0; pileupWeightUp=0.0; pileupWeightDown=0.0;
  rho=-1;

  Flag2_all = false;
  metPhi = -999.;
  met= -999.;
  metEENoise = -999.;metPhiEENoise = -999.; // metSF = -777;

   genWeight=-999.;

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

   nGenParticles = 0;
   for( int i = 0; i < N_MAX_GPARTICLES; i++ )
   {
      gParticleE[i] =0.0;
      gParticlePt[i] = 0.0;
      gParticleEta[i] = 0.0;
      gParticlePhi[i] = 0.0;
      gParticleMotherId[i] = 999999;
      gParticleId[i] = 999999;
      gParticleMotherIndex[i] = 999999;
   }

   nLeptons = 0;
   for( int i = 0; i < N_MAX_LEPTONS; i++ )
   {

     lepE[i]      = -999.;
     lepPt[i]     = -999.;
     lepEta[i]    = -999.;
     lepPhi[i]    = -999.;
     lepPdgId[i]  = -999;
     lepDZ[i]     = -999.;
     lepDXY[i]    = -999.;
     lepDXYErr[i] = -999.;
     lepSF[i]     = -999.;
     lepSFup[i]   = -999.;
     lepSFdn[i]   = -999.;

     lepLooseId[i] = false;
     lepTightId[i] = false;

     lepPassLooseIso[i] = false;
     lepPassTightIso[i] = false;
     lepPassVTightIso[i] = false;
     lepMuonType[i] = 0;
     lepMuonQuality[i] = 0;
     for (int q=0;q<MAX_MuonHLTFilters;q++) lepMuon_passHLTFilter[i][q] = false;

   }

   nTracks = 0;
   for( int i = 0; i < N_MAX_TRACKS; i++ )
   {
     track_Pt [i] = -999.;
     track_Eta[i] = -999.;
     track_Phi[i] = -999.;
   }

   nJets = 0;
   for( int i = 0; i < N_MAX_JETS; i++ )
   {
     jetJecUnc[i]    = -999.;
     jetE[i]         = -999.;
     jetEJESDown[i]  = -999.;
     jetEJESUp[i]    = -999.;
     jetPt[i]        = -999.;
     jetPtJESDown[i] = -999.;
     jetPtJESUp[i]   = -999.;
     jetEta[i]    = -999.;
     jetPhi[i]    = -999.;
     jetCISV[i]    = -999.;
     jetCMVA[i]    = -999.;
     jetTightPassId[i] = false;
   }

   for(int i = 0; i <NTriggersMAX; i++){
     HLTDecision[i] = false;
   }
    //CSC

  nCscRechitClusters = 0;
  // nDtRechitClusters = 0;

  nCscRechits = 0;
  nDtRechits = 0;

  nDtRings = 0;
  nCscRings = 0;
  nDtStations25 = 0;
  nDtWheels25 = 0;

  for( int i = 0; i < N_MAX_CSC; i++ )
  {
      cscRechitCluster_match_dtSeg_0p4[i] = 0;
      cscRechitCluster_match_MB1Seg_0p4[i] = 0;
      cscRechitCluster_match_RB1_0p4[i] = 0;
      cscRechitCluster_match_RE12_0p4[i] = 0;

      cscRechitCluster_match_gLLP_deltaR[i] = 999.;

      cscRechitClusterSize[i] = -999;
      cscRechitClusterX[i] = -999.;
      cscRechitClusterY[i] = -999.;
      cscRechitClusterZ[i] = -999.;
      cscRechitClusterTime[i]         = -999.;
      cscRechitClusterTimeWeighted[i] = -999.;
      cscRechitClusterTimeTotal[i]    = -999.;

      cscRechitClusterGenMuonDeltaR[i] = 999.;
      cscRechitClusterTimeSpread[i]            = -999.;
      cscRechitClusterTimeSpreadWeighted[i]    = -999.;
      cscRechitClusterTimeSpreadWeightedAll[i] = -999.;


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
      cscRechitClusterMuonVetoPt[i] = 0.0;
      cscRechitClusterGenMuonVetoPt[i] = 0.0;
      cscRechitClusterGenMuonVetoPt_dR0p8[i] = 0.0;
      cscRechitClusterMetEENoise_dPhi[i] = 999.;

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
      dtRechitClusterMuonVetoPt[i] = 0.0;
      dtRechitClusterGenMuonVetoPt[i] = 0.0;
      dtRechitClusterGenMuonVetoPt_dR0p8[i] = 0.0;
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

      dtRechitClusterWheel[i] = -999;


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
      }

};

void TreeMuonSystemBParkingExt::InitTree()
{
  assert(tree_);
  InitVariables();

  tree_->SetBranchAddress("runNum",      &runNum);
  tree_->SetBranchAddress("MC_condition",      &MC_condition);
  tree_->SetBranchAddress("lumiSec",     &lumiSec);
  tree_->SetBranchAddress("evtNum",      &evtNum);

  tree_->SetBranchAddress("npv",         &npv);
  tree_->SetBranchAddress("rho",         &rho);

  tree_->SetBranchAddress("npu",         &npu);
  tree_->SetBranchAddress("pileupWeight",     &pileupWeight);
  tree_->SetBranchAddress("pileupWeightUp",   &pileupWeightUp);
  tree_->SetBranchAddress("pileupWeightDown", &pileupWeightDown);
  tree_->SetBranchAddress("genWeight", &genWeight);

  tree_->SetBranchAddress("Flag2_all",      &Flag2_all);

  tree_->SetBranchAddress("metEENoise",      &metEENoise);
  tree_->SetBranchAddress("metPhiEENoise",      &metPhiEENoise);
  tree_->SetBranchAddress("met",      &met);
  tree_->SetBranchAddress("metPhi",      &metPhi);

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

  tree_->SetBranchAddress("nGenParticles",    &nGenParticles);
  tree_->SetBranchAddress("gParticleE",        gParticleE);
  tree_->SetBranchAddress("gParticlePt",       gParticlePt);
  tree_->SetBranchAddress("gParticleEta",      gParticleEta);
  tree_->SetBranchAddress("gParticlePhi",      gParticlePhi);
  tree_->SetBranchAddress("gParticleMotherId",      gParticleMotherId);
  tree_->SetBranchAddress("gParticleId",      gParticleId);
  tree_->SetBranchAddress("gParticleMotherIndex",      gParticleMotherIndex);

  tree_->SetBranchAddress("nLeptons",    &nLeptons);
  tree_->SetBranchAddress("lepE",        lepE);
  tree_->SetBranchAddress("lepPt",       lepPt);
  tree_->SetBranchAddress("lepEta",      lepEta);
  tree_->SetBranchAddress("lepPhi",      lepPhi);
  tree_->SetBranchAddress("lepPdgId",  lepPdgId);
  tree_->SetBranchAddress("lepDZ",     lepDZ);
  tree_->SetBranchAddress("lepDXY",     lepDXY);
  tree_->SetBranchAddress("lepDXYErr",     lepDXYErr);
  tree_->SetBranchAddress("lepSF",     lepSF);
  tree_->SetBranchAddress("lepSFup",     lepSFup);
  tree_->SetBranchAddress("lepSFdn",     lepSFdn);
  tree_->SetBranchAddress("lepMuonType",     lepMuonType);
  tree_->SetBranchAddress("lepMuonQuality",     lepMuonQuality);
  tree_->SetBranchAddress("lepMuon_passHLTFilter",     lepMuon_passHLTFilter);

  tree_->SetBranchAddress("lepLooseId", lepLooseId);
  tree_->SetBranchAddress("lepTightId", lepTightId);

  tree_->SetBranchAddress("lepPassLooseIso", lepPassLooseIso);
  tree_->SetBranchAddress("lepPassTightIso", lepPassTightIso);
  tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);

  //jets
  tree_->SetBranchAddress("nJets",     &nJets);
  tree_->SetBranchAddress("jetEJESUp",   jetEJESUp);
  tree_->SetBranchAddress("jetEJESDown", jetEJESDown);
  tree_->SetBranchAddress("jetE",      jetE);
  tree_->SetBranchAddress("jetJecUnc",      jetJecUnc);
  tree_->SetBranchAddress("jetPtJESUp",     jetPtJESUp);
  tree_->SetBranchAddress("jetPtJESDown",     jetPtJESDown);
  tree_->SetBranchAddress("jetPt",     jetPt);
  tree_->SetBranchAddress("jetEta",    jetEta);
  tree_->SetBranchAddress("jetPhi",    jetPhi);
  tree_->SetBranchAddress("jetCISV",    jetCISV);
  tree_->SetBranchAddress("jetCMVA",    jetCMVA);
  tree_->SetBranchAddress("jetTightPassId", jetTightPassId);
  // triggers
  tree_->SetBranchAddress("HLTDecision",   HLTDecision);
  tree_->SetBranchAddress("nCscRechits",   &nCscRechits);

  tree_->SetBranchAddress("nCscRings",             &nCscRings);


  tree_->SetBranchAddress("nDtRechits",            &nDtRechits);
  tree_->SetBranchAddress("nDtRings",             &nDtRings);

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

  tree_->SetBranchAddress("dtRechitClusterWheel",             dtRechitClusterWheel);


  tree_->SetBranchAddress("dtRechitClusterEta",             dtRechitClusterEta);
  tree_->SetBranchAddress("dtRechitClusterPhi",             dtRechitClusterPhi);


  tree_->SetBranchAddress("dtRechitClusterJetVetoPt",             dtRechitClusterJetVetoPt);

  tree_->SetBranchAddress("dtRechitClusterMuonVetoPt",             dtRechitClusterMuonVetoPt);
  tree_->SetBranchAddress("dtRechitClusterGenMuonVetoPt",             dtRechitClusterGenMuonVetoPt);
  tree_->SetBranchAddress("dtRechitClusterGenMuonVetoPt_dR0p8",             dtRechitClusterGenMuonVetoPt_dR0p8);

  tree_->SetBranchAddress("dtRechitClusterMetEENoise_dPhi",             dtRechitClusterMetEENoise_dPhi);
  tree_->SetBranchAddress("dtRechitCluster_match_gLLP_deltaR",             dtRechitCluster_match_gLLP_deltaR);


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

  tree_->SetBranchAddress("nCscRechitClusters",             &nCscRechitClusters);
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
  tree_->SetBranchAddress("cscRechitClusterTimeWeighted",     cscRechitClusterTimeWeighted);
  tree_->SetBranchAddress("cscRechitClusterTimeTotal",        cscRechitClusterTimeTotal);

  tree_->SetBranchAddress("cscRechitClusterGenMuonDeltaR",             cscRechitClusterGenMuonDeltaR);

  tree_->SetBranchAddress("cscRechitClusterTimeSpread",             cscRechitClusterTimeSpread);
  tree_->SetBranchAddress("cscRechitClusterTimeSpreadWeighted",     cscRechitClusterTimeSpreadWeighted);
  tree_->SetBranchAddress("cscRechitClusterTimeSpreadWeightedAll",  cscRechitClusterTimeSpreadWeightedAll);

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

  tree_->SetBranchAddress("cscRechitClusterMuonVetoPt",             cscRechitClusterMuonVetoPt);
  tree_->SetBranchAddress("cscRechitClusterGenMuonVetoPt_dR0p8",             cscRechitClusterGenMuonVetoPt_dR0p8);
  tree_->SetBranchAddress("cscRechitClusterGenMuonVetoPt",             cscRechitClusterGenMuonVetoPt);
  tree_->SetBranchAddress("cscRechitClusterMetEENoise_dPhi",            cscRechitClusterMetEENoise_dPhi);
  tree_->SetBranchAddress("cscRechitCluster_match_dtSeg_0p4",           cscRechitCluster_match_dtSeg_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_MB1Seg_0p4",          cscRechitCluster_match_MB1Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_RB1_0p4",             cscRechitCluster_match_RB1_0p4);
  tree_->SetBranchAddress("cscRechitCluster_match_RE12_0p4",            cscRechitCluster_match_RE12_0p4);

  tree_->SetBranchAddress("cscRechitCluster_match_gLLP_deltaR",         cscRechitCluster_match_gLLP_deltaR);

};

void TreeMuonSystemBParkingExt::LoadTree(const char* file)
{
  f_ = TFile::Open(file);
  assert(f_);
  tree_ = dynamic_cast<TTree*>(f_->Get("MuonSystem"));
  InitTree();
  assert(tree_);
};

void TreeMuonSystemBParkingExt::CreateTree()
{
  tree_ = new TTree("MuonSystem","MuonSystem");
  f_ = 0;

  tree_->Branch("runNum",      &runNum,     "runNum/i");      // event run number
  tree_->Branch("MC_condition",      &MC_condition,     "MC_condition/i");      // event run number
  tree_->Branch("lumiSec",     &lumiSec,    "lumiSec/i");     // event lumi section
  tree_->Branch("evtNum",      &evtNum,     "evtNum/i");      // event number
  tree_->Branch("rho",         &rho,        "rho/F");
  tree_->Branch("npv",         &npv,        "npv/I");         // number of primary vertices

  tree_->Branch("npu",         &npu,        "npu/I");         // number of pileup
  tree_->Branch("pileupWeight",         &pileupWeight,        "pileupWeight/F"); // pileup weight
  tree_->Branch("pileupWeightUp",       &pileupWeightUp,      "pileupWeightUp/F"); // pileup weight up
  tree_->Branch("pileupWeightDown",     &pileupWeightDown,    "pileupWeightDown/F"); // pileup weight down
  tree_->Branch("genWeight",     &genWeight,    "genWeight/F"); 

  tree_->Branch("Flag2_all",      &Flag2_all,     "Flag2_all/O");
  tree_->Branch("metEENoise",      &metEENoise,     "metEENoise/F");      // phi(MET)
  // tree_->Branch("metSF",      &metSF,     "metSF/F");      // SF
  tree_->Branch("metPhiEENoise",      &metPhiEENoise,     "metPhiEENoise/F");      // phi(MET)
  tree_->Branch("met",      &met,     "met/F");      // phi(MET)
  tree_->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)

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

  tree_->Branch("nGenParticles",    &nGenParticles, "nGenParticles/I");
  tree_->Branch("gParticleE",        gParticleE,      "gParticleE[nGenParticles]/F");
  tree_->Branch("gParticlePt",        gParticlePt,      "gParticlePt[nGenParticles]/F");
  tree_->Branch("gParticleEta",        gParticleEta,      "gParticleEta[nGenParticles]/F");
  tree_->Branch("gParticlePhi",        gParticlePhi,      "gParticlePhi[nGenParticles]/F");
  tree_->Branch("gParticleId",        gParticleId,      "gParticleId[nGenParticles]/I");
  tree_->Branch("gParticleMotherId",        gParticleMotherId,      "gParticleMotherId[nGenParticles]/I");
  tree_->Branch("gParticleMotherIndex",        gParticleMotherIndex,      "gParticleMotherIndex[nGenParticles]/I");

  //leptons
  tree_->Branch("nLeptons",  &nLeptons, "nLeptons/I");
  tree_->Branch("lepE",      lepE,      "lepE[nLeptons]/F");
  tree_->Branch("lepPt",     lepPt,     "lepPt[nLeptons]/F");
  tree_->Branch("lepMuonType",   lepMuonType,     "lepMuonType[nLeptons]/i");
  tree_->Branch("lepMuonQuality",     lepMuonQuality,     "lepMuonQuality[nLeptons]/i");
  //tree_->Branch("lep_passSingleMuTagFilter",     lep_passSingleMuTagFilter,     "lep_passSingleMuTagFilter[nLeptons]/O");
  tree_->Branch("lepMuon_passHLTFilter", &lepMuon_passHLTFilter, Form("lepMuon_passHLTFilter[nLeptons][%d]/O",MAX_MuonHLTFilters));

  tree_->Branch("lepEta",    lepEta,    "lepEta[nLeptons]/F");
  tree_->Branch("lepPhi",    lepPhi,    "lepPhi[nLeptons]/F");
  tree_->Branch("lepPdgId",  lepPdgId,  "lepPdgId[nLeptons]/I");
  tree_->Branch("lepDZ",     lepDZ,     "lepDZ[nLeptons]/F");
  tree_->Branch("lepDXY",     lepDXY,     "lepDXY[nLeptons]/F");
  tree_->Branch("lepDXYErr",     lepDXYErr,     "lepDXYErr[nLeptons]/F");
  tree_->Branch("lepSF",     lepSF,     "lepSF[nLeptons]/F");
  tree_->Branch("lepSFup",     lepSFup,     "lepSFup[nLeptons]/F");
  tree_->Branch("lepSFdn",     lepSFdn,     "lepSFdn[nLeptons]/F");
  tree_->Branch("lepLooseId", lepLooseId, "lepLooseId[nLeptons]/O");
  tree_->Branch("lepTightId", lepTightId, "lepTightId[nLeptons]/O");
  tree_->Branch("lepPassLooseIso", lepPassLooseIso, "lepPassLooseIso[nLeptons]/O");
  tree_->Branch("lepPassTightIso", lepPassTightIso, "lepPassTightIso[nLeptons]/O");
  tree_->Branch("lepPassVTightIso", lepPassVTightIso, "lepPassVTightIso[nLeptons]/O");

  //jets
  tree_->Branch("nJets",     &nJets,    "nJets/I");
  tree_->Branch("jetEJESDown",  jetEJESDown,    "jetEJESDown[nJets]/F");
  tree_->Branch("jetEJESUp",    jetEJESUp,      "jetEJESUp[nJets]/F");
  tree_->Branch("jetE",      jetE,      "jetE[nJets]/F");
  tree_->Branch("jetJecUnc",      jetJecUnc,      "jetJecUnc[nJets]/F");
  tree_->Branch("jetPtJESDown", jetPtJESDown,   "jetPtJESDown[nJets]/F");
  tree_->Branch("jetPtJESUp",   jetPtJESUp,     "jetPtJESUp[nJets]/F");
  tree_->Branch("jetPt",     jetPt,     "jetPt[nJets]/F");
  tree_->Branch("jetEta",    jetEta,    "jetEta[nJets]/F");
  tree_->Branch("jetPhi",    jetPhi,    "jetPhi[nJets]/F");
  tree_->Branch("jetCISV",    jetCISV,    "jetCISV[nJets]/F");
  tree_->Branch("jetCMVA",    jetCMVA,    "jetCMVA[nJets]/F");
  tree_->Branch("jetTightPassId", jetTightPassId, "jetTightPassId[nJets]/O");

  tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[1201]/O"); //hardcoded

  tree_->Branch("nCscRechits",             &nCscRechits, "nCscRechits/I");
  tree_->Branch("nCscRings",             &nCscRings, "nCscRings/I");


  tree_->Branch("nDtRechits",            &nDtRechits,             "nDtRechits/I");
  tree_->Branch("nDtRings",             &nDtRings, "nDtRings/I");

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

  tree_->Branch("cscRechitCluster_match_dtSeg_0p4",             cscRechitCluster_match_dtSeg_0p4,           "cscRechitCluster_match_dtSeg_0p4[nCscRechitClusters]/I");
  tree_->Branch("cscRechitCluster_match_MB1Seg_0p4",            cscRechitCluster_match_MB1Seg_0p4,          "cscRechitCluster_match_MB1Seg_0p4[nCscRechitClusters]/I");
  tree_->Branch("cscRechitCluster_match_RB1_0p4",               cscRechitCluster_match_RB1_0p4,             "cscRechitCluster_match_RB1_0p4[nCscRechitClusters]/I");
  tree_->Branch("cscRechitCluster_match_RE12_0p4",              cscRechitCluster_match_RE12_0p4,            "cscRechitCluster_match_RE12_0p4[nCscRechitClusters]/I");
  tree_->Branch("cscRechitCluster_match_gLLP_deltaR",           cscRechitCluster_match_gLLP_deltaR,         "cscRechitCluster_match_gLLP_deltaR[nCscRechitClusters]/F");

  tree_->Branch("cscRechitClusterJetVetoPt",             cscRechitClusterJetVetoPt,             "cscRechitClusterJetVetoPt[nCscRechitClusters]/F");
  tree_->Branch("cscRechitClusterMuonVetoPt",             cscRechitClusterMuonVetoPt,             "cscRechitClusterMuonVetoPt[nCscRechitClusters]/F");
  tree_->Branch("cscRechitClusterGenMuonVetoPt_dR0p8",             cscRechitClusterGenMuonVetoPt_dR0p8,             "cscRechitClusterGenMuonVetoPt_dR0p8[nCscRechitClusters]/F");
  tree_->Branch("cscRechitClusterGenMuonVetoPt",             cscRechitClusterGenMuonVetoPt,             "cscRechitClusterGenMuonVetoPt[nCscRechitClusters]/F");
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

  tree_->Branch("dtRechitClusterX",             dtRechitClusterX,             "dtRechitClusterX[nDtRechitClusters]/F");
  tree_->Branch("dtRechitClusterY",             dtRechitClusterY,             "dtRechitClusterY[nDtRechitClusters]/F");
  tree_->Branch("dtRechitClusterZ",             dtRechitClusterZ,             "dtRechitClusterZ[nDtRechitClusters]/F");
  tree_->Branch("dtRechitClusterWheel",             dtRechitClusterWheel,             "dtRechitClusterWheel[nDtRechitClusters]/I");

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

  tree_->Branch("dtRechitClusterJetVetoE",             dtRechitClusterJetVetoE,             "dtRechitClusterJetVetoE[nDtRechitClusters]/F");
  tree_->Branch("dtRechitClusterGenMuonVetoPt",             dtRechitClusterGenMuonVetoPt,             "dtRechitClusterGenMuonVetoPt[nDtRechitClusters]/F");
  tree_->Branch("dtRechitClusterGenMuonVetoPt_dR0p8",             dtRechitClusterGenMuonVetoPt_dR0p8,             "dtRechitClusterGenMuonVetoPt_dR0p8[nDtRechitClusters]/F");
  tree_->Branch("dtRechitClusterMuonVetoPt",             dtRechitClusterMuonVetoPt,             "dtRechitClusterMuonVetoPt[nDtRechitClusters]/F");
  tree_->Branch("dtRechitClusterMetEENoise_dPhi",             dtRechitClusterMetEENoise_dPhi,             "dtRechitClusterMetEENoise_dPhi[nDtRechitClusters]/F");

};
