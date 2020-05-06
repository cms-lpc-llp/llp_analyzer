// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef SusyLLPTree_H
#define SusyLLPTree_H

#define LLP_ARRAY_SIZE 2
#define LLP_DAUGHTER_ARRAY_SIZE 4
#define LLP_GRAND_DAUGHTER_ARRAY_SIZE 4
#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 2000
#define NTriggersMAX 602 //Number of trigger in the .dat file
#define N_CSC_CUT 20
#define JET_PT_CUT 10
#define MUON_PT_CUT 20

#define OBJECTARRAYSIZE 100
#define CSCRECHITARRAYSIZE 100000
#define RECHITARRAYSIZE 20000
#define HORECHITARRAYSIZE 2000
#define GENPARTICLEARRAYSIZE 2000
#define MAX_NPV 1000
#define MAX_NPFCAND 5000
#define MAX_NPU 1000
#define MAX_NBX 1000
#define LLP_ARRAY_SIZE 2
#define LLP_DAUGHTER_ARRAY_SIZE 4
#define LLP_GRAND_DAUGHTER_ARRAY_SIZE 4
#define STRIP_DIGI_THRESHOLD 13.3

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

class SusyLLPTree
{

public:
  SusyLLPTree();
  ~SusyLLPTree();

  //tree
  TTree *tree_;
  TFile *f_;

  //event info
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t pileupWeight, pileupWeightUp, pileupWeightDown;
  Float_t rho, weight;
  Float_t met, metPhi;
  Float_t HT;
  Float_t jetMet_dPhi;
  Float_t jetMet_dPhiStar;
  Float_t jetMet_dPhiMin;
  Float_t jetMet_dPhiStarMin;
  Float_t jetMet_dPhiMin4;

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
  bool lepPassId[N_MAX_LEPTONS];

  //Z-candidate
  float MT;
  float ZMass;
  float ZPt;
  float ZEta;
  float ZPhi;
  int ZleptonIndex1;
  int ZleptonIndex2;

  //h-candidate
  float dr_max;
  float dm_hh; 
  float mh1; 
  float mh2; 
  float avg_mh; 

/*
 int nTaus;
 int nMuons;
 int nPhotons;
 int nElectrons;
*/

 //muons
 int nMuons;
 float muonE[OBJECTARRAYSIZE];
 float muonPt[OBJECTARRAYSIZE];
 float muonEta[OBJECTARRAYSIZE];
 float muonPhi[OBJECTARRAYSIZE];
 int muonCharge[OBJECTARRAYSIZE];//muon charge
 bool muonIsLoose[OBJECTARRAYSIZE];
 bool muonIsMedium[OBJECTARRAYSIZE];
 bool muonIsTight[OBJECTARRAYSIZE];
/*
 float muon_d0[OBJECTARRAYSIZE];//transverse impact paramenter
 float muon_dZ[OBJECTARRAYSIZE];//impact parameter
 float muon_ip3d[OBJECTARRAYSIZE];//3d impact paramenter
 float muon_ip3dSignificance[OBJECTARRAYSIZE];//3d impact paramenter/error
 unsigned int muonType[OBJECTARRAYSIZE];//muonTypeBit: global, tracker, standalone
 unsigned int muonQuality[OBJECTARRAYSIZE];//muonID Quality Bits
*/
 float muon_pileupIso[OBJECTARRAYSIZE];
 float muon_chargedIso[OBJECTARRAYSIZE];
 float muon_photonIso[OBJECTARRAYSIZE];
 float muon_neutralHadIso[OBJECTARRAYSIZE];
/*
 float muon_ptrel[OBJECTARRAYSIZE];
 float muon_chargedMiniIso[OBJECTARRAYSIZE];
 float muon_photonAndNeutralHadronMiniIso[OBJECTARRAYSIZE];
 float muon_chargedPileupMiniIso[OBJECTARRAYSIZE];
 float muon_activityMiniIsoAnnulus[OBJECTARRAYSIZE];
 bool  muon_passSingleMuTagFilter[OBJECTARRAYSIZE];
 bool  muon_passHLTFilter[OBJECTARRAYSIZE][MAX_MuonHLTFilters];
 float muon_validFractionTrackerHits[OBJECTARRAYSIZE];
 bool  muon_isGlobal[OBJECTARRAYSIZE];
 float muon_normChi2[OBJECTARRAYSIZE];
 float muon_chi2LocalPosition[OBJECTARRAYSIZE];
 float muon_kinkFinder[OBJECTARRAYSIZE];
 float muon_segmentCompatability[OBJECTARRAYSIZE];
 bool muonIsICHEPMedium[OBJECTARRAYSIZE];
*/

 //electrons
 int nElectrons;
 float eleE[OBJECTARRAYSIZE];
 float elePt[OBJECTARRAYSIZE];
 float eleEta[OBJECTARRAYSIZE];
 float elePhi[OBJECTARRAYSIZE];
 float eleCharge[OBJECTARRAYSIZE];
/*
 float eleE_SC[OBJECTARRAYSIZE];
 float eleEta_SC[OBJECTARRAYSIZE];
 float elePhi_SC[OBJECTARRAYSIZE];
 float eleSigmaIetaIeta[OBJECTARRAYSIZE];
 float eleFull5x5SigmaIetaIeta[OBJECTARRAYSIZE];
 float eleR9[OBJECTARRAYSIZE];
 float ele_dEta[OBJECTARRAYSIZE];
 float ele_dPhi[OBJECTARRAYSIZE];
 float ele_HoverE[OBJECTARRAYSIZE];
 float ele_d0[OBJECTARRAYSIZE];
 float ele_dZ[OBJECTARRAYSIZE];
 float ele_ip3d[OBJECTARRAYSIZE];
 float ele_ip3dSignificance[OBJECTARRAYSIZE];
 float ele_pileupIso[OBJECTARRAYSIZE];
 float ele_chargedIso[OBJECTARRAYSIZE];
 float ele_photonIso[OBJECTARRAYSIZE];
 float ele_neutralHadIso[OBJECTARRAYSIZE];
 int   ele_MissHits[OBJECTARRAYSIZE];
 bool  ele_PassConvVeto[OBJECTARRAYSIZE];
 float ele_OneOverEminusOneOverP[OBJECTARRAYSIZE];
 float ele_IDMVAGeneralPurpose[OBJECTARRAYSIZE];
 int   ele_IDMVACategoryGeneralPurpose[OBJECTARRAYSIZE];
 float ele_IDMVAHZZ[OBJECTARRAYSIZE];
 int   ele_IDMVACategoryHZZ[OBJECTARRAYSIZE];
 float ele_RegressionE[OBJECTARRAYSIZE];
 float ele_CombineP4[OBJECTARRAYSIZE];
 float ele_ptrel[OBJECTARRAYSIZE];
 float ele_chargedMiniIso[OBJECTARRAYSIZE];
 float ele_photonAndNeutralHadronMiniIso[OBJECTARRAYSIZE];
 float ele_chargedPileupMiniIso[OBJECTARRAYSIZE];
 float ele_activityMiniIsoAnnulus[OBJECTARRAYSIZE];
*/
 bool ele_passCutBasedIDVeto[OBJECTARRAYSIZE];
 bool ele_passCutBasedIDLoose[OBJECTARRAYSIZE];
 bool ele_passCutBasedIDMedium[OBJECTARRAYSIZE];
 bool ele_passCutBasedIDTight[OBJECTARRAYSIZE];
/*
 bool ele_passMVAIsoIDWP80[OBJECTARRAYSIZE];
 bool ele_passMVAIsoIDWP90[OBJECTARRAYSIZE];
 bool ele_passMVAIsoIDWPHZZ[OBJECTARRAYSIZE];
 bool ele_passMVAIsoIDWPLoose[OBJECTARRAYSIZE];
 bool ele_passMVANoIsoIDWP80[OBJECTARRAYSIZE];
 bool ele_passMVANoIsoIDWP90[OBJECTARRAYSIZE];
 bool ele_passMVANoIsoIDWPLoose[OBJECTARRAYSIZE];
 bool ele_passSingleEleTagFilter[OBJECTARRAYSIZE];
 bool ele_passTPOneTagFilter[OBJECTARRAYSIZE];
 bool ele_passTPTwoTagFilter[OBJECTARRAYSIZE];
 bool ele_passTPOneProbeFilter[OBJECTARRAYSIZE];
 bool ele_passTPTwoProbeFilter[OBJECTARRAYSIZE];
 bool ele_passHLTFilter[OBJECTARRAYSIZE][MAX_ElectronHLTFilters];
 vector<vector<uint> > ele_EcalRechitID;
 vector<vector<uint> > *ele_EcalRechitIndex;
 vector<uint> ele_SeedRechitID;
 vector<uint> *ele_SeedRechitIndex;
*/

 //taus
 int nTaus;
 float tauE[OBJECTARRAYSIZE];
 float tauPt[OBJECTARRAYSIZE];
 float tauEta[OBJECTARRAYSIZE];
 float tauPhi[OBJECTARRAYSIZE];
 bool tau_IsLoose[OBJECTARRAYSIZE];
 bool tau_IsMedium[OBJECTARRAYSIZE];
 bool tau_IsTight[OBJECTARRAYSIZE];
 bool tau_passEleVetoLoose[OBJECTARRAYSIZE];
 bool tau_passEleVetoMedium[OBJECTARRAYSIZE];
 bool tau_passEleVetoTight[OBJECTARRAYSIZE];
 bool tau_passMuVetoLoose[OBJECTARRAYSIZE];
 bool tau_passMuVetoMedium[OBJECTARRAYSIZE];
 bool tau_passMuVetoTight[OBJECTARRAYSIZE];
/*
 UInt_t tau_ID[OBJECTARRAYSIZE];//tauID Bits
 float tau_combinedIsoDeltaBetaCorr3Hits[OBJECTARRAYSIZE];
 float tau_chargedIsoPtSum[OBJECTARRAYSIZE];
 float tau_neutralIsoPtSum[OBJECTARRAYSIZE];
 float tau_puCorrPtSum[OBJECTARRAYSIZE];
 float tau_eleVetoMVA[OBJECTARRAYSIZE];
 int tau_eleVetoCategory[OBJECTARRAYSIZE];
 float tau_muonVetoMVA[OBJECTARRAYSIZE];
 float tau_isoMVAnewDMwLT[OBJECTARRAYSIZE];
 float tau_isoMVAnewDMwoLT[OBJECTARRAYSIZE];
 float tau_leadCandPt[OBJECTARRAYSIZE];
 int tau_leadCandID[OBJECTARRAYSIZE];
 float tau_leadChargedHadrCandPt[OBJECTARRAYSIZE];
 int tau_leadChargedHadrCandID[OBJECTARRAYSIZE];
*/

 //photons
 int nPhotons;
 //int nPhotons_overlap;
 float phoE[OBJECTARRAYSIZE];
 float phoPt[OBJECTARRAYSIZE];
 float phoEta[OBJECTARRAYSIZE];
 float phoPhi[OBJECTARRAYSIZE];
/*
 float phoSigmaIetaIeta[OBJECTARRAYSIZE];
 float phoFull5x5SigmaIetaIeta[OBJECTARRAYSIZE];
 float phoR9[OBJECTARRAYSIZE];
 float pho_sminor[OBJECTARRAYSIZE];
 float pho_smajor[OBJECTARRAYSIZE];
 float pho_HoverE[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPtAllVertices[OBJECTARRAYSIZE][MAX_NPV];
 float pho_sumChargedHadronPt[OBJECTARRAYSIZE];
 float pho_sumNeutralHadronEt[OBJECTARRAYSIZE];
 float pho_sumPhotonEt[OBJECTARRAYSIZE];
 float pho_ecalPFClusterIso[OBJECTARRAYSIZE];
 float pho_hcalPFClusterIso[OBJECTARRAYSIZE];
 float pho_trkSumPtHollowConeDR03[OBJECTARRAYSIZE];
 float pho_sumWorstVertexChargedHadronPt[OBJECTARRAYSIZE];
 float pho_pfIsoChargedHadronIso[OBJECTARRAYSIZE];
 float pho_pfIsoChargedHadronIsoWrongVtx[OBJECTARRAYSIZE];
 float pho_pfIsoNeutralHadronIso[OBJECTARRAYSIZE];
 float pho_pfIsoPhotonIso[OBJECTARRAYSIZE];
 float pho_pfIsoModFrixione[OBJECTARRAYSIZE];
 float pho_pfIsoSumPUPt[OBJECTARRAYSIZE];
 bool  pho_isConversion[OBJECTARRAYSIZE];
 bool  pho_passEleVeto[OBJECTARRAYSIZE];
 float pho_RegressionE[OBJECTARRAYSIZE];
 float pho_RegressionEUncertainty[OBJECTARRAYSIZE];
 float pho_IDMVA[OBJECTARRAYSIZE];
 float pho_superClusterEnergy[OBJECTARRAYSIZE];
 float pho_superClusterRawEnergy[OBJECTARRAYSIZE];
 float pho_superClusterEta[OBJECTARRAYSIZE];
 float pho_superClusterPhi[OBJECTARRAYSIZE];
 float pho_superClusterX[OBJECTARRAYSIZE];
 float pho_superClusterY[OBJECTARRAYSIZE];
 float pho_superClusterZ[OBJECTARRAYSIZE];
 bool pho_hasPixelSeed[OBJECTARRAYSIZE];
 bool pho_passHLTFilter[OBJECTARRAYSIZE][MAX_PhotonHLTFilters];
*/
 bool pho_passCutBasedIDLoose[OBJECTARRAYSIZE];
 bool pho_passCutBasedIDMedium[OBJECTARRAYSIZE];
 bool pho_passCutBasedIDTight[OBJECTARRAYSIZE];
 bool pho_passMVAIDWP80[OBJECTARRAYSIZE];
 bool pho_passMVAIDWP90[OBJECTARRAYSIZE];
/*
 int pho_convType[OBJECTARRAYSIZE];
 float pho_convTrkZ[OBJECTARRAYSIZE];
 float pho_convTrkClusZ[OBJECTARRAYSIZE];
 float pho_vtxSumPx[OBJECTARRAYSIZE][MAX_NPV];
 float pho_vtxSumPy[OBJECTARRAYSIZE][MAX_NPV];
 bool  pho_isStandardPhoton[OBJECTARRAYSIZE];
 bool  pho_seedRecHitSwitchToGain6[OBJECTARRAYSIZE];
 bool  pho_seedRecHitSwitchToGain1[OBJECTARRAYSIZE];
 bool  pho_anyRecHitSwitchToGain6[OBJECTARRAYSIZE];
 bool  pho_anyRecHitSwitchToGain1[OBJECTARRAYSIZE];
 vector<vector<uint> > pho_EcalRechitID;
 vector<vector<uint> > *pho_EcalRechitIndex;
 vector<uint>  pho_SeedRechitID;
 vector<uint>  *pho_SeedRechitIndex;
//extra stuff
 float pho_sumChargedHadronPt_NewPV_NoTiming[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing120_TrkVtx[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing50_TrkPho[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing80_TrkPho[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing100_TrkPho[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing120_TrkPho[OBJECTARRAYSIZE];
 float pho_superClusterSeedX[OBJECTARRAYSIZE];
 float pho_superClusterSeedY[OBJECTARRAYSIZE];
 float pho_superClusterSeedZ[OBJECTARRAYSIZE];
 float pho_superClusterSeedT[OBJECTARRAYSIZE];
 float pho_superClusterSeedE[OBJECTARRAYSIZE];
 float pho_pfClusterSeedE[OBJECTARRAYSIZE];
*/

  //jets
 int nJets;
 float jetE[N_MAX_JETS];
 float jetPt[N_MAX_JETS];
 float jetEta[N_MAX_JETS];
 float jetPhi[N_MAX_JETS];
 float jetCSV[N_MAX_JETS];
 float jetCISV[N_MAX_JETS];
 float jetProbb[N_MAX_JETS];
 float jetProbc[N_MAX_JETS];
 float jetProbudsg[N_MAX_JETS];
 float jetProbbb[N_MAX_JETS];
 float jetMass[N_MAX_JETS];
 float jetJetArea[N_MAX_JETS];
 float jetPileupE[N_MAX_JETS];
 float jetPileupId[N_MAX_JETS];
 int   jetPileupIdFlag[N_MAX_JETS];
 bool  jetPassIDLoose[N_MAX_JETS];
 bool  jetPassIDTight[N_MAX_JETS];
 bool  jetPassMuFrac[N_MAX_JETS];
 bool  jetPassEleFrac[N_MAX_JETS];
 int   jetPartonFlavor[N_MAX_JETS];
 int   jetHadronFlavor[N_MAX_JETS];
 float jetElectronEnergyFraction[N_MAX_JETS];
 float jetPhotonEnergyFraction[N_MAX_JETS];
 float jetChargedHadronEnergyFraction[N_MAX_JETS];
 float jetNeutralHadronEnergyFraction[N_MAX_JETS];
 float jetMuonEnergyFraction[N_MAX_JETS];
 float jetHOEnergyFraction[N_MAX_JETS];
 float jetHFHadronEnergyFraction[N_MAX_JETS];
 float jetHFEMEnergyFraction[N_MAX_JETS];
 int   jetChargedHadronMultiplicity[N_MAX_JETS];
 int   jetNeutralHadronMultiplicity[N_MAX_JETS];
 int   jetElectronMultiplicity[N_MAX_JETS];
 int   jetPhotonMultiplicity[N_MAX_JETS];
 int   jetMuonMultiplicity[N_MAX_JETS];
 float jetAllMuonPt[N_MAX_JETS];
 float jetAllMuonEta[N_MAX_JETS];
 float jetAllMuonPhi[N_MAX_JETS];
 float jetAllMuonM[N_MAX_JETS];
 float jetPtWeightedDZ[N_MAX_JETS];
 int   jetNRechits[N_MAX_JETS];
 float jetRechitE[N_MAX_JETS];
 float jetRechitT[N_MAX_JETS];
 float jetRechitT_rms[N_MAX_JETS];
 float jetRechitE_Error[N_MAX_JETS];
 float jetRechitT_Error[N_MAX_JETS];

 float jetGammaMax[N_MAX_JETS];
 float jetGammaMax_ET[N_MAX_JETS];
 float jetGammaMax_EM[N_MAX_JETS];
 float jetGammaMax_Hadronic[N_MAX_JETS];
 float jetAlphaMax[N_MAX_JETS];
 float jetBetaMax[N_MAX_JETS];

 float jetPtAllTracks[N_MAX_JETS];
 float jetPtAllPVTracks[N_MAX_JETS];
 float jetMedianTheta2D[N_MAX_JETS];
 float jetMedianIP[N_MAX_JETS];
 float jetMinDeltaRAllTracks[N_MAX_JETS];
 float jetMinDeltaRPVTracks[N_MAX_JETS];

 float jet_energy_frac[N_MAX_JETS];
 float jet_sig_et1[N_MAX_JETS];
 float jet_sig_et2[N_MAX_JETS];
 bool jet_matched[N_MAX_JETS];
 bool jet_matched_gLLP0_grandaughter[N_MAX_JETS];
 bool jet_matched_gLLP1_grandaughter[N_MAX_JETS];

 float jetGammaMax_wp[N_MAX_JETS];
 float jetGammaMax_ET_wp[N_MAX_JETS];
 float jetGammaMax_EM_wp[N_MAX_JETS];
 float jetGammaMax_Hadronic_wp[N_MAX_JETS];
 float jetAlphaMax_wp[N_MAX_JETS];
 float jetBetaMax_wp[N_MAX_JETS];

 float jetPtAllTracks_wp[N_MAX_JETS];
 float jetPtAllPVTracks_wp[N_MAX_JETS];
 float jetMedianTheta2D_wp[N_MAX_JETS];
 float jetMedianIP_wp[N_MAX_JETS];
 float jetMinDeltaRAllTracks_wp[N_MAX_JETS];
 float jetMinDeltaRPVTracks_wp[N_MAX_JETS];

  float jetTime[N_MAX_JETS];
  float ecalNRechits[N_MAX_JETS];
  float ecalRechitE[N_MAX_JETS];

  float jetChargedEMEnergyFraction[N_MAX_JETS];
  float jetNeutralEMEnergyFraction[N_MAX_JETS];
  /*
  int nJets;
  int nBJets;
  float jetE[N_MAX_JETS];
  float jetEt[N_MAX_JETS];
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
  float jetHoverE[N_MAX_JETS];
  float jetGammaMax_ET[N_MAX_JETS];
  float jetMinDeltaRPVTracks[N_MAX_JETS];
  float jetPtAllPVTracks[N_MAX_JETS];
  float jetMinDeltaRAllTracks[N_MAX_JETS];
  float jetPtAllTracks[N_MAX_JETS];
  float jetGammaMax_ET_wp[N_MAX_JETS];
  float jetMinDeltaRPVTracks_wp[N_MAX_JETS];
  float jetPtAllPVTracks_wp[N_MAX_JETS];
  float jetMinDeltaRAllTracks_wp[N_MAX_JETS];
  float jetPtAllTracks_wp[N_MAX_JETS];
  // bool jetLoosePassId[N_MAX_JETS];
  bool jetPassId[N_MAX_JETS];
  bool jetCSVT[N_MAX_JETS];
  float jetCISV[N_MAX_JETS];
  bool matched[N_MAX_JETS];
  bool jet_matched_gLLP0_daughter[N_MAX_JETS];
  bool jet_matched_gLLP1_daughter[N_MAX_JETS];
  bool jet_matched_gLLP0_grandaughter[N_MAX_JETS];
  bool jet_matched_gLLP1_grandaughter[N_MAX_JETS];
  // bool jetTightPassId[N_MAX_JETS];
  float jet_energy_frac[N_MAX_JETS];
  float jet_sig_et1[N_MAX_JETS];
  float jet_sig_et2[N_MAX_JETS];
*/

  // met filters
   bool Flag2_globalSuperTightHalo2016Filter;
   bool Flag2_globalTightHalo2016Filter;
   bool Flag2_goodVertices;
   bool Flag2_BadChargedCandidateFilter;
   bool Flag2_BadPFMuonFilter;
   bool Flag2_EcalDeadCellTriggerPrimitiveFilter;
   bool Flag2_HBHENoiseFilter;
   bool Flag2_HBHEIsoNoiseFilter;
   bool Flag2_ecalBadCalibFilter;
   bool Flag2_eeBadScFilter;

   //MC
   int nGenJets;
   float genJetE[OBJECTARRAYSIZE];
   float genJetPt[OBJECTARRAYSIZE];
   float genJetEta[OBJECTARRAYSIZE];
   float genJetPhi[OBJECTARRAYSIZE];
   float genJetMET[OBJECTARRAYSIZE];
   float genMetPtCalo;
   float genMetPhiCalo;
   float genMetPtTrue;
   float genMetPhiTrue;
   float genVertexX;
   float genVertexY;
   float genVertexZ;
   float genVertexT;
   float genWeight;
/*
   unsigned int genSignalProcessID;
   float genQScale;
   float genAlphaQCD;
   float genAlphaQED;
   string lheComments;
   vector<float> *scaleWeights;
   vector<float> *pdfWeights;
   vector<float> *alphasWeights;
  
   int firstPdfWeight;
   int lastPdfWeight;
   int firstAlphasWeight;
   int lastAlphasWeight;
*/
   //gen info
   int nGenParticle;
   int gParticleMotherId[GENPARTICLEARRAYSIZE];
   int gParticleMotherIndex[GENPARTICLEARRAYSIZE];
   int gParticleId[GENPARTICLEARRAYSIZE];
   int gParticleStatus[GENPARTICLEARRAYSIZE];
   float gParticleE[GENPARTICLEARRAYSIZE];
   float gParticlePt[GENPARTICLEARRAYSIZE];
   float gParticlePx[GENPARTICLEARRAYSIZE];
   float gParticlePy[GENPARTICLEARRAYSIZE];
   float gParticlePz[GENPARTICLEARRAYSIZE];
   float gParticleEta[GENPARTICLEARRAYSIZE];
   float gParticlePhi[GENPARTICLEARRAYSIZE];
  
   float gParticleProdVertexX[GENPARTICLEARRAYSIZE];
   float gParticleProdVertexY[GENPARTICLEARRAYSIZE];
   float gParticleProdVertexZ[GENPARTICLEARRAYSIZE];
  
   float gParticleDecayVertexX[GENPARTICLEARRAYSIZE];
   float gParticleDecayVertexY[GENPARTICLEARRAYSIZE];
   float gParticleDecayVertexZ[GENPARTICLEARRAYSIZE];
   float gLLP_prod_vertex_x[LLP_ARRAY_SIZE];
   float gLLP_prod_vertex_y[LLP_ARRAY_SIZE];
   float gLLP_prod_vertex_z[LLP_ARRAY_SIZE];
   float gLLP_decay_vertex_x[LLP_ARRAY_SIZE];
   float gLLP_decay_vertex_y[LLP_ARRAY_SIZE];
   float gLLP_decay_vertex_z[LLP_ARRAY_SIZE];
   float gLLP_beta[LLP_ARRAY_SIZE];
   float gLLP_travel_time[LLP_ARRAY_SIZE];
   float gLLP_pt[LLP_ARRAY_SIZE];
   float gLLP_e[LLP_ARRAY_SIZE];
   float gLLP_eta[LLP_ARRAY_SIZE];
   float gLLP_phi[LLP_ARRAY_SIZE];
   bool gLLP_csc[LLP_ARRAY_SIZE];
   bool gLLP_dt[LLP_ARRAY_SIZE];
  
   float photon_travel_time[LLP_DAUGHTER_ARRAY_SIZE];
   float photon_travel_time_pv[LLP_DAUGHTER_ARRAY_SIZE];
  
   float gen_time[LLP_DAUGHTER_ARRAY_SIZE];
   float gen_time_pv[LLP_DAUGHTER_ARRAY_SIZE];
   float gLLP_daughter_travel_time[LLP_DAUGHTER_ARRAY_SIZE];
   int   gLLP_daughter_id[LLP_DAUGHTER_ARRAY_SIZE];
   float gLLP_daughter_pt[LLP_DAUGHTER_ARRAY_SIZE];
   float gLLP_daughter_eta[LLP_DAUGHTER_ARRAY_SIZE];
   float gLLP_daughter_phi[LLP_DAUGHTER_ARRAY_SIZE];
   float gLLP_daughter_eta_ecalcorr[LLP_DAUGHTER_ARRAY_SIZE];
   float gLLP_daughter_phi_ecalcorr[LLP_DAUGHTER_ARRAY_SIZE];
   float gLLP_daughter_e[LLP_DAUGHTER_ARRAY_SIZE];
   float gLLP_daughter_mass[LLP_DAUGHTER_ARRAY_SIZE];
  
    //grandaughters
   bool gLLP_grandaughter_EB[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
   bool gLLP_grandaughter_ETL[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  
   float gLLP_grandaughter_photon_travel_time_EB[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
   float gLLP_grandaughter_photon_travel_time_ETL[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  
   float gLLP_grandaughter_travel_time_EB[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
   float gLLP_grandaughter_travel_time_ETL[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  
   float gen_time_grandaughter_EB[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
   float gen_time_grandaughter_ETL[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  
   int   gLLP_grandaughter_id[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
   float gLLP_grandaughter_pt[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
   float gLLP_grandaughter_eta[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
   float gLLP_grandaughter_phi[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
   float gLLP_grandaughter_eta_ecalcorr[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
   float gLLP_grandaughter_phi_ecalcorr[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
   float gLLP_grandaughter_e[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
   float gLLP_grandaughter_mass[LLP_GRAND_DAUGHTER_ARRAY_SIZE];

  //MET
 float metType1Pt;
 float metType1Px;
 float metType1Py;
 float metType1Eta;
 float metType1Phi; 

  //HLT 
  bool HLTDecision[NTriggersMAX];

  UInt_t wzevtNum,trig, trig_lepId, trig_lepId_dijet; //number of events that pass each criteria

  //functions
  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();

};
#endif
