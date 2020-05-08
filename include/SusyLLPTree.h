// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef SusyLLPTree_H
#define SusyLLPTree_H

#define LLP_ARRAY_SIZE 2
#define LLP_DAUGHTER_ARRAY_SIZE 4
#define LLP_GRAND_DAUGHTER_ARRAY_SIZE 4
#define N_MAX_LEPTONS 20
#define N_MAX_JETS 20
#define N_MAX_CSC 2000
#define NTriggersMAX 602 //Number of trigger in the .dat file
#define N_CSC_CUT 20
#define JET_PT_CUT 10
#define MUON_PT_CUT 20

#define N_MAX_MUONS 20
#define N_MAX_ELECTRONS 20
#define N_MAX_TAUS 20
#define N_MAX_PHOTONS 20

#define NTriggersMAX 602 //Number of trigger in the .dat file
#define OBJSIZE 20
#define CSCRECHITARRAYSIZE 100000
#define RECHITARRAYSIZE 20000
#define HORECHITARRAYSIZE 2000
#define GENPARTICLEARRAYSIZE 2000
#define MAX_NPV 20
#define MAX_NPFCAND 5000
#define MAX_NPU 20
#define MAX_NBX 20
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
		float muonE[OBJSIZE];
		float muonPt[OBJSIZE];
		float muonEta[OBJSIZE];
		float muonPhi[OBJSIZE];
		int muonCharge[OBJSIZE];//muon charge
		bool muonIsLoose[OBJSIZE];
		bool muonIsMedium[OBJSIZE];
		bool muonIsTight[OBJSIZE];
		/*
		   float muon_d0[OBJSIZE];//transverse impact paramenter
		   float muon_dZ[OBJSIZE];//impact parameter
		   float muon_ip3d[OBJSIZE];//3d impact paramenter
		   float muon_ip3dSignificance[OBJSIZE];//3d impact paramenter/error
		   unsigned int muonType[OBJSIZE];//muonTypeBit: global, tracker, standalone
		   unsigned int muonQuality[OBJSIZE];//muonID Quality Bits
		   */
		float muon_pileupIso[OBJSIZE];
		float muon_chargedIso[OBJSIZE];
		float muon_photonIso[OBJSIZE];
		float muon_neutralHadIso[OBJSIZE];
		/*
		   float muon_ptrel[OBJSIZE];
		   float muon_chargedMiniIso[OBJSIZE];
		   float muon_photonAndNeutralHadronMiniIso[OBJSIZE];
		   float muon_chargedPileupMiniIso[OBJSIZE];
		   float muon_activityMiniIsoAnnulus[OBJSIZE];
		   bool  muon_passSingleMuTagFilter[OBJSIZE];
		   bool  muon_passHLTFilter[OBJSIZE][MAX_MuonHLTFilters];
		   float muon_validFractionTrackerHits[OBJSIZE];
		   bool  muon_isGlobal[OBJSIZE];
		   float muon_normChi2[OBJSIZE];
		   float muon_chi2LocalPosition[OBJSIZE];
		   float muon_kinkFinder[OBJSIZE];
		   float muon_segmentCompatability[OBJSIZE];
		   bool muonIsICHEPMedium[OBJSIZE];
		   */

		//electrons
		int nElectrons;
		float eleE[OBJSIZE];
		float elePt[OBJSIZE];
		float eleEta[OBJSIZE];
		float elePhi[OBJSIZE];
		float eleCharge[OBJSIZE];
		/*
		   float eleE_SC[OBJSIZE];
		   float eleEta_SC[OBJSIZE];
		   float elePhi_SC[OBJSIZE];
		   float eleSigmaIetaIeta[OBJSIZE];
		   float eleFull5x5SigmaIetaIeta[OBJSIZE];
		   float eleR9[OBJSIZE];
		   float ele_dEta[OBJSIZE];
		   float ele_dPhi[OBJSIZE];
		   float ele_HoverE[OBJSIZE];
		   float ele_d0[OBJSIZE];
		   float ele_dZ[OBJSIZE];
		   float ele_ip3d[OBJSIZE];
		   float ele_ip3dSignificance[OBJSIZE];
		   float ele_pileupIso[OBJSIZE];
		   float ele_chargedIso[OBJSIZE];
		   float ele_photonIso[OBJSIZE];
		   float ele_neutralHadIso[OBJSIZE];
		   int   ele_MissHits[OBJSIZE];
		   bool  ele_PassConvVeto[OBJSIZE];
		   float ele_OneOverEminusOneOverP[OBJSIZE];
		   float ele_IDMVAGeneralPurpose[OBJSIZE];
		   int   ele_IDMVACategoryGeneralPurpose[OBJSIZE];
		   float ele_IDMVAHZZ[OBJSIZE];
		   int   ele_IDMVACategoryHZZ[OBJSIZE];
		   float ele_RegressionE[OBJSIZE];
		   float ele_CombineP4[OBJSIZE];
		   float ele_ptrel[OBJSIZE];
		   float ele_chargedMiniIso[OBJSIZE];
		   float ele_photonAndNeutralHadronMiniIso[OBJSIZE];
		   float ele_chargedPileupMiniIso[OBJSIZE];
		   float ele_activityMiniIsoAnnulus[OBJSIZE];
		   */
		bool ele_passCutBasedIDVeto[OBJSIZE];
		bool ele_passCutBasedIDLoose[OBJSIZE];
		bool ele_passCutBasedIDMedium[OBJSIZE];
		bool ele_passCutBasedIDTight[OBJSIZE];
		/*
		   bool ele_passMVAIsoIDWP80[OBJSIZE];
		   bool ele_passMVAIsoIDWP90[OBJSIZE];
		   bool ele_passMVAIsoIDWPHZZ[OBJSIZE];
		   bool ele_passMVAIsoIDWPLoose[OBJSIZE];
		   bool ele_passMVANoIsoIDWP80[OBJSIZE];
		   bool ele_passMVANoIsoIDWP90[OBJSIZE];
		   bool ele_passMVANoIsoIDWPLoose[OBJSIZE];
		   bool ele_passSingleEleTagFilter[OBJSIZE];
		   bool ele_passTPOneTagFilter[OBJSIZE];
		   bool ele_passTPTwoTagFilter[OBJSIZE];
		   bool ele_passTPOneProbeFilter[OBJSIZE];
		   bool ele_passTPTwoProbeFilter[OBJSIZE];
		   bool ele_passHLTFilter[OBJSIZE][MAX_ElectronHLTFilters];
		   vector<vector<uint> > ele_EcalRechitID;
		   vector<vector<uint> > *ele_EcalRechitIndex;
		   vector<uint> ele_SeedRechitID;
		   vector<uint> *ele_SeedRechitIndex;
		   */

		//taus
		int nTaus;
		float tauE[OBJSIZE];
		float tauPt[OBJSIZE];
		float tauEta[OBJSIZE];
		float tauPhi[OBJSIZE];
		bool tau_IsLoose[OBJSIZE];
		bool tau_IsMedium[OBJSIZE];
		bool tau_IsTight[OBJSIZE];
		bool tau_passEleVetoLoose[OBJSIZE];
		bool tau_passEleVetoMedium[OBJSIZE];
		bool tau_passEleVetoTight[OBJSIZE];
		bool tau_passMuVetoLoose[OBJSIZE];
		bool tau_passMuVetoMedium[OBJSIZE];
		bool tau_passMuVetoTight[OBJSIZE];
		/*
		   UInt_t tau_ID[OBJSIZE];//tauID Bits
		   float tau_combinedIsoDeltaBetaCorr3Hits[OBJSIZE];
		   float tau_chargedIsoPtSum[OBJSIZE];
		   float tau_neutralIsoPtSum[OBJSIZE];
		   float tau_puCorrPtSum[OBJSIZE];
		   float tau_eleVetoMVA[OBJSIZE];
		   int tau_eleVetoCategory[OBJSIZE];
		   float tau_muonVetoMVA[OBJSIZE];
		   float tau_isoMVAnewDMwLT[OBJSIZE];
		   float tau_isoMVAnewDMwoLT[OBJSIZE];
		   float tau_leadCandPt[OBJSIZE];
		   int tau_leadCandID[OBJSIZE];
		   float tau_leadChargedHadrCandPt[OBJSIZE];
		   int tau_leadChargedHadrCandID[OBJSIZE];
		   */

		//photons
		int nPhotons;
		//int nPhotons_overlap;
		float phoE[OBJSIZE];
		float phoPt[OBJSIZE];
		float phoEta[OBJSIZE];
		float phoPhi[OBJSIZE];
		/*
		   float phoSigmaIetaIeta[OBJSIZE];
		   float phoFull5x5SigmaIetaIeta[OBJSIZE];
		   float phoR9[OBJSIZE];
		   float pho_sminor[OBJSIZE];
		   float pho_smajor[OBJSIZE];
		   float pho_HoverE[OBJSIZE];
		   float pho_sumChargedHadronPtAllVertices[OBJSIZE][MAX_NPV];
		   float pho_sumChargedHadronPt[OBJSIZE];
		   float pho_sumNeutralHadronEt[OBJSIZE];
		   float pho_sumPhotonEt[OBJSIZE];
		   float pho_ecalPFClusterIso[OBJSIZE];
		   float pho_hcalPFClusterIso[OBJSIZE];
		   float pho_trkSumPtHollowConeDR03[OBJSIZE];
		   float pho_sumWorstVertexChargedHadronPt[OBJSIZE];
		   float pho_pfIsoChargedHadronIso[OBJSIZE];
		   float pho_pfIsoChargedHadronIsoWrongVtx[OBJSIZE];
		   float pho_pfIsoNeutralHadronIso[OBJSIZE];
		   float pho_pfIsoPhotonIso[OBJSIZE];
		   float pho_pfIsoModFrixione[OBJSIZE];
		   float pho_pfIsoSumPUPt[OBJSIZE];
		   bool  pho_isConversion[OBJSIZE];
		   bool  pho_passEleVeto[OBJSIZE];
		   float pho_RegressionE[OBJSIZE];
		   float pho_RegressionEUncertainty[OBJSIZE];
		   float pho_IDMVA[OBJSIZE];
		   float pho_superClusterEnergy[OBJSIZE];
		   float pho_superClusterRawEnergy[OBJSIZE];
		   float pho_superClusterEta[OBJSIZE];
		   float pho_superClusterPhi[OBJSIZE];
		   float pho_superClusterX[OBJSIZE];
		   float pho_superClusterY[OBJSIZE];
		   float pho_superClusterZ[OBJSIZE];
		   bool pho_hasPixelSeed[OBJSIZE];
		   bool pho_passHLTFilter[OBJSIZE][MAX_PhotonHLTFilters];
		   */
		bool pho_passCutBasedIDLoose[OBJSIZE];
		bool pho_passCutBasedIDMedium[OBJSIZE];
		bool pho_passCutBasedIDTight[OBJSIZE];
		bool pho_passMVAIDWP80[OBJSIZE];
		bool pho_passMVAIDWP90[OBJSIZE];
		/*
		   int pho_convType[OBJSIZE];
		   float pho_convTrkZ[OBJSIZE];
		   float pho_convTrkClusZ[OBJSIZE];
		   float pho_vtxSumPx[OBJSIZE][MAX_NPV];
		   float pho_vtxSumPy[OBJSIZE][MAX_NPV];
		   bool  pho_isStandardPhoton[OBJSIZE];
		   bool  pho_seedRecHitSwitchToGain6[OBJSIZE];
		   bool  pho_seedRecHitSwitchToGain1[OBJSIZE];
		   bool  pho_anyRecHitSwitchToGain6[OBJSIZE];
		   bool  pho_anyRecHitSwitchToGain1[OBJSIZE];
		   vector<vector<uint> > pho_EcalRechitID;
		   vector<vector<uint> > *pho_EcalRechitIndex;
		   vector<uint>  pho_SeedRechitID;
		   vector<uint>  *pho_SeedRechitIndex;
		//extra stuff
		float pho_sumChargedHadronPt_NewPV_NoTiming[OBJSIZE];
		float pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx[OBJSIZE];
		float pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx[OBJSIZE];
		float pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx[OBJSIZE];
		float pho_sumChargedHadronPt_NewPV_Timing120_TrkVtx[OBJSIZE];
		float pho_sumChargedHadronPt_NewPV_Timing50_TrkPho[OBJSIZE];
		float pho_sumChargedHadronPt_NewPV_Timing80_TrkPho[OBJSIZE];
		float pho_sumChargedHadronPt_NewPV_Timing100_TrkPho[OBJSIZE];
		float pho_sumChargedHadronPt_NewPV_Timing120_TrkPho[OBJSIZE];
		float pho_superClusterSeedX[OBJSIZE];
		float pho_superClusterSeedY[OBJSIZE];
		float pho_superClusterSeedZ[OBJSIZE];
		float pho_superClusterSeedT[OBJSIZE];
		float pho_superClusterSeedE[OBJSIZE];
		float pho_pfClusterSeedE[OBJSIZE];
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

		float jetEcalE[N_MAX_JETS];
		float jetHcalE[N_MAX_JETS];
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
		float genJetE[OBJSIZE];
		float genJetPt[OBJSIZE];
		float genJetEta[OBJSIZE];
		float genJetPhi[OBJSIZE];
		float genJetMET[OBJSIZE];
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
