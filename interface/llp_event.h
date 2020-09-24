//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 10 19:28:06 2019 by ROOT version 6.10/08
// from TTree llp_event/selected AOD information for llp_event analyses
// found on file: llp_ntupler_126.root
//////////////////////////////////////////////////////////

#ifndef llp_event_h
#define llp_event_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

#define OBJECTARRAYSIZE 100
#define MAX_NPFCAND 1000
#define RECHITARRAYSIZE 1000

class llp_event {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Fixed size dimensions of array or collections stored in the TTree if any.

		// Declaration of leaf types
		Bool_t          isData;
		Int_t           nPV;
		UInt_t          runNum;
		UInt_t          nSlimmedSecondV;
		UInt_t          lumiNum;
		UInt_t          eventNum;
		UInt_t          eventTime;
		Float_t         pvX;
		Float_t         pvY;
		Float_t         pvZ;
		Float_t         fixedGridRhoAll;
		Float_t         fixedGridRhoFastjetAll;
		Float_t         fixedGridRhoFastjetAllCalo;
		Float_t         fixedGridRhoFastjetCentralCalo;
		Float_t         fixedGridRhoFastjetCentralChargedPileUp;
		Float_t         fixedGridRhoFastjetCentralNeutral;
		Int_t           nPVAll;
		Float_t         pvAllX[200];   //[nPVAll]
		Float_t         pvAllY[200];   //[nPVAll]
		Float_t         pvAllZ[200];   //[nPVAll]
		Float_t         pvAllLogSumPtSq[200];   //[nPVAll]
		Float_t         pvAllSumPx[200];   //[nPVAll]
		Float_t         pvAllSumPy[200];   //[nPVAll]
		Int_t           nBunchXing;
		Int_t           BunchXing[50];   //[nBunchXing]
		Int_t           nPU[50];   //[nBunchXing]
		Float_t         nPUmean[50];   //[nBunchXing]
		Int_t           nMuons;
		Float_t         muonE[2000];   //[nMuons]
		Float_t         muonPt[2000];   //[nMuons]
		Float_t         muonEta[2000];   //[nMuons]
		Float_t         muonPhi[2000];   //[nMuons]
		Int_t           muonCharge[2000];   //[nMuons]
		Bool_t          muonIsLoose[2000];   //[nMuons]
		Bool_t          muonIsMedium[2000];   //[nMuons]
		Bool_t          muonIsTight[2000];   //[nMuons]
		Float_t         muon_d0[2000];   //[nMuons]
		Float_t         muon_dZ[2000];   //[nMuons]
		Float_t         muon_ip3d[2000];   //[nMuons]
		Float_t         muon_ip3dSignificance[2000];   //[nMuons]
		UInt_t          muonType[2000];   //[nMuons]
		UInt_t          muonQuality[2000];   //[nMuons]
		Float_t         muon_pileupIso[2000];   //[nMuons]
		Float_t         muon_chargedIso[2000];   //[nMuons]
		Float_t         muon_photonIso[2000];   //[nMuons]
		Float_t         muon_neutralHadIso[2000];   //[nMuons]
		Float_t         muon_ptrel[2000];   //[nMuons]
		Float_t         muon_chargedMiniIso[2000];   //[nMuons]
		Float_t         muon_photonAndNeutralHadronMiniIso[2000];   //[nMuons]
		Float_t         muon_chargedPileupMiniIso[2000];   //[nMuons]
		Float_t         muon_activityMiniIsoAnnulus[2000];   //[nMuons]
		Bool_t          muon_passSingleMuTagFilter[2000];   //[nMuons]
		Bool_t          muon_passHLTFilter[2000][100];   //[nMuons]
		Float_t         muon_validFractionTrackerHits[2000];  //[nMuons]
		Bool_t          muon_isGlobal[2000];   //[nMuons]
		Float_t         muon_normChi2[2000];   //[nMuons]
		Float_t         muon_chi2LocalPosition[2000];   //[nMuons]
		Float_t         muon_kinkFinder[2000];   //[nMuons]
		Float_t         muon_segmentCompatability[2000];   //[nMuons]
		Bool_t          muonIsICHEPMedium[2000];   //[nMuons]
		Int_t           nElectrons;
		Float_t         eleE[2000];   //[nElectrons]
		Float_t         elePt[2000];   //[nElectrons]
		Float_t         eleEta[2000];   //[nElectrons]
		Float_t         elePhi[2000];   //[nElectrons]
		Float_t         eleCharge[2000];   //[nElectrons]
		Float_t         eleEta_SC[2000];   //[nElectrons]
		Float_t         eleSigmaIetaIeta[2000];   //[nElectrons]
		Float_t         eleFull5x5SigmaIetaIeta[2000];   //[nElectrons]
		Float_t         eleR9[2000];   //[nElectrons]
		Float_t         ele_dEta[2000];   //[nElectrons]
		Float_t         ele_dPhi[2000];   //[nElectrons]
		Float_t         ele_HoverE[2000];   //[nElectrons]
		Float_t         ele_d0[2000];   //[nElectrons]
		Float_t         ele_dZ[2000];   //[nElectrons]
		Float_t         ele_ip3d[2000];   //[nElectrons]
		Float_t         ele_ip3dSignificance[2000];   //[nElectrons]
		Float_t         ele_pileupIso[2000];   //[nElectrons]
		Float_t         ele_chargedIso[2000];   //[nElectrons]
		Float_t         ele_photonIso[2000];   //[nElectrons]
		Float_t         ele_neutralHadIso[2000];   //[nElectrons]
		Int_t           ele_MissHits[2000];   //[nElectrons]
		Bool_t          ele_PassConvVeto[2000];   //[nElectrons]
		Float_t         ele_OneOverEminusOneOverP[2000];  //[nElectrons]
		Float_t         ele_IDMVAGeneralPurpose[2000];   //[nElectrons]
		Int_t           ele_IDMVACategoryGeneralPurpose[2000];   //[nElectrons]
		Float_t         ele_IDMVAHZZ[2000];   //[nElectrons]
		Int_t           ele_IDMVACategoryHZZ[2000];   //[nElectrons]
		Float_t         ele_RegressionE[2000];   //[nElectrons]
		Float_t         ele_CombineP4[2000];   //[nElectrons]
		Float_t         ele_ptrel[2000];   //[nElectrons]
		Float_t         ele_chargedMiniIso[2000];   //[nElectrons]
		Float_t         ele_photonAndNeutralHadronMiniIso[2000];   //[nElectrons]
		Float_t         ele_chargedPileupMiniIso[2000];   //[nElectrons]
		Float_t         ele_activityMiniIsoAnnulus[2000];   //[nElectrons]
		Bool_t 	   ele_passCutBasedIDVeto[2000];
		Bool_t 	   ele_passCutBasedIDLoose[2000];
		Bool_t 	   ele_passCutBasedIDMedium[2000];
		Bool_t 	   ele_passCutBasedIDTight[2000];
		Bool_t 	   ele_passMVAIsoIDWP80[2000];
		Bool_t 	   ele_passMVAIsoIDWP90[2000];
		Bool_t 	   ele_passMVAIsoIDWPHZZ[2000];
		Bool_t 	   ele_passMVAIsoIDWPLoose[2000];
		Bool_t 	   ele_passMVANoIsoIDWP80[2000];
		Bool_t 	   ele_passMVANoIsoIDWP90[2000];
		Bool_t 	   ele_passMVANoIsoIDWPLoose[2000];
		Bool_t          ele_passSingleEleTagFilter[2000];   //[nElectrons]
		Bool_t          ele_passTPOneTagFilter[2000];   //[nElectrons]
		Bool_t          ele_passTPTwoTagFilter[2000];   //[nElectrons]
		Bool_t          ele_passTPOneProbeFilter[2000];   //[nElectrons]
		Bool_t          ele_passTPTwoProbeFilter[2000];   //[nElectrons]
		Bool_t          ele_passHLTFilter[2000][100];   //[nElectrons]
		Int_t           nTaus;
		Float_t         tauE[2000];   //[nTaus]
		Float_t         tauPt[2000];   //[nTaus]
		Float_t         tauEta[2000];   //[nTaus]
		Float_t         tauPhi[2000];   //[nTaus]
		Bool_t          tau_IsLoose[2000];   //[nTaus]
		Bool_t          tau_IsMedium[2000];   //[nTaus]
		Bool_t          tau_IsTight[2000];   //[nTaus]
		Bool_t          tau_passEleVetoLoose[2000];   //[nTaus]
		Bool_t          tau_passEleVetoMedium[2000];   //[nTaus]
		Bool_t          tau_passEleVetoTight[2000];   //[nTaus]
		Bool_t          tau_passMuVetoLoose[2000];   //[nTaus]
		Bool_t          tau_passMuVetoMedium[2000];   //[nTaus]
		Bool_t          tau_passMuVetoTight[2000];   //[nTaus]
		UInt_t          tau_ID[2000];   //[nTaus]
		Float_t         tau_combinedIsoDeltaBetaCorr3Hits[2000];   //[nTaus]
		Float_t         tau_chargedIsoPtSum[2000];   //[nTaus]
		Float_t         tau_neutralIsoPtSum[2000];   //[nTaus]
		Float_t         tau_puCorrPtSum[2000];   //[nTaus]
		Float_t         tau_eleVetoMVA[2000];   //[nTaus]
		Int_t           tau_eleVetoCategory[2000];   //[nTaus]
		Float_t         tau_muonVetoMVA[2000];   //[nTaus]
		Float_t         tau_isoMVAnewDMwLT[2000];   //[nTaus]
		Float_t         tau_isoMVAnewDMwoLT[2000];   //[nTaus]
		Float_t         tau_leadCandPt[2000];   //[nTaus]
		Int_t           tau_leadCandID[2000];   //[nTaus]
		Float_t         tau_leadChargedHadrCandPt[2000];   //[nTaus]
		Int_t           tau_leadChargedHadrCandID[2000];   //[nTaus]
		UInt_t          nIsoPFCandidates;
		Float_t         isoPFCandidatePt[2000];   //[nIsoPFCandidates]
		Float_t         isoPFCandidateEta[2000];   //[nIsoPFCandidates]
		Float_t         isoPFCandidatePhi[2000];   //[nIsoPFCandidates]
		Float_t         isoPFCandidateIso04[2000];  //[nIsoPFCandidates]
		Float_t         isoPFCandidateD0[2000];  //[nIsoPFCandidates]
		Int_t           isoPFCandidatePdgId[2000];  //[nIsoPFCandidates]
		Int_t           nPhotons;
		Int_t           nPhotons_overlap;
		Float_t         phoE[2000];   //[nPhotons]
		Float_t         phoPt[2000];   //[nPhotons]
		Float_t         phoEta[2000];   //[nPhotons]
		Float_t         phoPhi[2000];   //[nPhotons]
		Float_t         phoSigmaIetaIeta[2000];   //[nPhotons]
		Float_t         phoFull5x5SigmaIetaIeta[2000];   //[nPhotons]
		Float_t         phoR9[2000];   //[nPhotons]
		Float_t         pho_sminor[2000];   //[nPhotons]
		Float_t         pho_smajor[2000];   //[nPhotons]
		Float_t         pho_HoverE[2000];   //[nPhotons]
		Float_t         pho_sumChargedHadronPtAllVertices[2000][200];   //[nPhotons]
		Float_t         pho_sumChargedHadronPt[2000];   //[nPhotons]
		Float_t         pho_sumNeutralHadronEt[2000];   //[nPhotons]
		Float_t         pho_sumPhotonEt[2000];   //[nPhotons]
		Float_t         pho_ecalPFClusterIso[2000];   //[nPhotons]
		Float_t         pho_hcalPFClusterIso[2000];   //[nPhotons]
		Float_t         pho_trkSumPtHollowConeDR03[2000];   //[nPhotons]
		Float_t         pho_sumWorstVertexChargedHadronPt[2000];   //[nPhotons]
		Float_t         pho_pfIsoChargedHadronIso[2000];   //[nPhotons]
		Float_t         pho_pfIsoChargedHadronIsoWrongVtx[2000];   //[nPhotons]
		Float_t         pho_pfIsoNeutralHadronIso[2000];   //[nPhotons]
		Float_t         pho_pfIsoPhotonIso[2000];   //[nPhotons]
		Float_t         pho_pfIsoModFrixione[2000];   //[nPhotons]
		Float_t         pho_pfIsoSumPUPt[2000];   //[nPhotons]
		Bool_t          pho_isConversion[2000];   //[nPhotons]
		Bool_t          pho_passEleVeto[2000];   //[nPhotons]
		Float_t         pho_RegressionE[2000];   //[nPhotons]
		Float_t         pho_RegressionEUncertainty[2000];   //[nPhotons]
		Float_t         pho_IDMVA[2000];   //[nPhotons]
		Float_t         pho_superClusterEnergy[2000];   //[nPhotons]
		Float_t         pho_superClusterRawEnergy[2000];   //[nPhotons]
		Float_t         pho_superClusterEta[2000];   //[nPhotons]
		Float_t         pho_superClusterPhi[2000];   //[nPhotons]
		Float_t         pho_superClusterX[2000];   //[nPhotons]
		Float_t         pho_superClusterY[2000];   //[nPhotons]
		Float_t         pho_superClusterZ[2000];   //[nPhotons]
		Bool_t          pho_hasPixelSeed[2000];   //[nPhotons]
		Bool_t          pho_passHLTFilter[2000][100];   //[nPhotons]
		Bool_t	   pho_passCutBasedIDLoose[2000];
		Bool_t	   pho_passCutBasedIDMedium[2000];
		Bool_t	   pho_passCutBasedIDTight[2000];
		Bool_t	   pho_passMVAIDWP80[2000];
		Bool_t	   pho_passMVAIDWP90[2000];
		Int_t           pho_convType[2000];   //[nPhotons]
		Float_t         pho_convTrkZ[2000];   //[nPhotons]
		Float_t         pho_convTrkClusZ[2000];   //[nPhotons]
		Float_t         pho_vtxSumPx[2000][200];   //[nPhotons]
		Float_t         pho_vtxSumPy[2000][200];   //[nPhotons]
		Bool_t          pho_isStandardPhoton[2000];   //[nPhotons]
		Float_t         pho_seedRecHitSwitchToGain6[2000];  //[nPhotons]
		Float_t         pho_seedRecHitSwitchToGain1[2000];  //[nPhotons]
		Float_t         pho_anyRecHitSwitchToGain6[2000];  //[nPhotons]
		Float_t         pho_anyRecHitSwitchToGain1[2000];  //[nPhotons]
		std::vector<float>   *ecalRechit_Eta;
		std::vector<float>   *ecalRechit_Phi;
		std::vector<float>   *ecalRechit_X;
		std::vector<float>   *ecalRechit_Y;
		std::vector<float>   *ecalRechit_Z;
		std::vector<float>   *ecalRechit_E;
		std::vector<float>   *ecalRechit_T;
		std::vector<unsigned int> *ecalRechit_ID;
		std::vector<bool>    *ecalRechit_FlagOOT;
		std::vector<bool>    *ecalRechit_GainSwitch1;
		std::vector<bool>    *ecalRechit_GainSwitch6;
		std::vector<float>   *ecalRechit_transpCorr;
		Int_t           nCsc;
		Float_t         cscPhi[2000];   //[nCsc]
		Float_t         cscEta[2000];   //[nCsc]
		Float_t         cscX[2000];   //[nCsc]
		Float_t         cscY[2000];   //[nCsc]
		Float_t         cscZ[2000];   //[nCsc]
		Float_t         cscDirectionX[2000];   //[nCsc]
		Float_t         cscDirectionY[2000];   //[nCsc]
		Float_t         cscDirectionZ[2000];   //[nCsc]
		Float_t         cscNRecHits[2000];   //[nCsc]
		Float_t         cscNRecHits_flag[2000];   //[nCsc]
		Float_t         cscT[2000];   //[nCsc]
		Float_t         cscChi2[2000];   //[nCsc]
		Int_t           nRpc;
		Float_t         rpcPhi[2000];   //[nRpc]
		Float_t         rpcEta[2000];   //[nRpc]
		Float_t         rpcX[2000];   //[nRpc]
		Float_t         rpcY[2000];   //[nRpc]
		Float_t         rpcZ[2000];   //[nRpc]
		Float_t         rpcT[2000];   //[nRpc]
		Float_t         rpcTError[2000];   //[nRpc]
		Int_t           nDt;
		Float_t         dtPhi[2000];   //[nDt]
		Float_t         dtEta[2000];   //[nDt]
		Float_t         dtX[2000];   //[nDt]
		Float_t         dtY[2000];   //[nDt]
		Float_t         dtZ[2000];   //[nDt]
		Float_t         dtDirX[2000];   //[nDt]
		Float_t         dtDirY[2000];   //[nDt]
		Float_t         dtDirZ[2000];   //[nDt]
		Float_t         dtT[2000];   //[nDt]
		Float_t         dtTError[2000];   //[nDt]
		Int_t           nCaloJets;
		//AK4 Jets
		int nJets;
		float jetE[OBJECTARRAYSIZE];
		float jetPt[OBJECTARRAYSIZE];
		float jetEta[OBJECTARRAYSIZE];
		float jetPhi[OBJECTARRAYSIZE];
		float jetCSV[OBJECTARRAYSIZE];
		float jetCISV[OBJECTARRAYSIZE];
		float jetProbb[OBJECTARRAYSIZE];
		float jetProbc[OBJECTARRAYSIZE];
		float jetProbudsg[OBJECTARRAYSIZE];
		float jetProbbb[OBJECTARRAYSIZE];
		float jetMass[OBJECTARRAYSIZE];
		float jetJetArea[OBJECTARRAYSIZE];
		float jetPileupE[OBJECTARRAYSIZE];
		float jetPileupId[OBJECTARRAYSIZE];
		int   jetPileupIdFlag[OBJECTARRAYSIZE];
		bool  jetPassIDLoose[OBJECTARRAYSIZE];
		bool  jetPassIDTight[OBJECTARRAYSIZE];
		bool  jetPassMuFrac[OBJECTARRAYSIZE];
		bool  jetPassEleFrac[OBJECTARRAYSIZE];
		int   jetPartonFlavor[OBJECTARRAYSIZE];
		int   jetHadronFlavor[OBJECTARRAYSIZE];
		float jetElectronEnergyFraction[OBJECTARRAYSIZE];
		float jetPhotonEnergyFraction[OBJECTARRAYSIZE];
		float jetChargedHadronEnergyFraction[OBJECTARRAYSIZE];
		float jetNeutralHadronEnergyFraction[OBJECTARRAYSIZE];
		float jetMuonEnergyFraction[OBJECTARRAYSIZE];
		float jetHOEnergyFraction[OBJECTARRAYSIZE];
		float jetHFHadronEnergyFraction[OBJECTARRAYSIZE];
		float jetHFEMEnergyFraction[OBJECTARRAYSIZE];
		int   jetChargedHadronMultiplicity[OBJECTARRAYSIZE];
		int   jetNeutralHadronMultiplicity[OBJECTARRAYSIZE];
		int   jetElectronMultiplicity[OBJECTARRAYSIZE];
		int   jetPhotonMultiplicity[OBJECTARRAYSIZE];
		int   jetMuonMultiplicity[OBJECTARRAYSIZE];
		float jetQGLikelihood[OBJECTARRAYSIZE];
		int   jetNSV[OBJECTARRAYSIZE];
		int   jetNSVCand[OBJECTARRAYSIZE];
		int   jetNVertexTracks[OBJECTARRAYSIZE];
		int   jetNSelectedTracks[OBJECTARRAYSIZE];
		float jetDRSVJet[OBJECTARRAYSIZE];
		float jetFlightDist2D[OBJECTARRAYSIZE];
		float jetFlightDist2DError[OBJECTARRAYSIZE];
		float jetFlightDist3D[OBJECTARRAYSIZE];
		float jetFlightDist3DError[OBJECTARRAYSIZE];
		float jetSV_x[OBJECTARRAYSIZE];
		float jetSV_y[OBJECTARRAYSIZE];
		float jetSV_z[OBJECTARRAYSIZE];
		int   jetSVNTracks[OBJECTARRAYSIZE];
		float jetSVMass[OBJECTARRAYSIZE];
		float jetAllMuonPt[OBJECTARRAYSIZE];
		float jetAllMuonEta[OBJECTARRAYSIZE];
		float jetAllMuonPhi[OBJECTARRAYSIZE];
		float jetAllMuonM[OBJECTARRAYSIZE];
		float jetPtWeightedDZ[OBJECTARRAYSIZE];
		int   jetNRechits[OBJECTARRAYSIZE];
		float jetRechitE[OBJECTARRAYSIZE];
		float jetRechitT[OBJECTARRAYSIZE];
		float jetRechitT_rms[OBJECTARRAYSIZE];
		float jetRechitE_Error[OBJECTARRAYSIZE];
		float jetRechitT_Error[OBJECTARRAYSIZE];

		float jetGammaMax[OBJECTARRAYSIZE];
		float jetGammaMax_ET[OBJECTARRAYSIZE];
		float jetGammaMax_EM[OBJECTARRAYSIZE];
		float jetGammaMax_Hadronic[OBJECTARRAYSIZE];
		float jetAlphaMax[OBJECTARRAYSIZE];
		float jetBetaMax[OBJECTARRAYSIZE];

		float jetPtAllTracks[OBJECTARRAYSIZE];
		float jetPtAllPVTracks[OBJECTARRAYSIZE];
		float jetMedianTheta2D[OBJECTARRAYSIZE];
		float jetMedianIP[OBJECTARRAYSIZE];
		float jetMinDeltaRAllTracks[OBJECTARRAYSIZE];
		float jetMinDeltaRPVTracks[OBJECTARRAYSIZE];

		float jet_energy_frac[OBJECTARRAYSIZE];
		float jet_sig_et1[OBJECTARRAYSIZE];
		float jet_sig_et2[OBJECTARRAYSIZE];
		bool jet_matched[OBJECTARRAYSIZE];
		bool jet_matched_gLLP0_grandaughter[OBJECTARRAYSIZE];
		bool jet_matched_gLLP1_grandaughter[OBJECTARRAYSIZE];

		float jetGammaMax_wp[OBJECTARRAYSIZE];
		float jetGammaMax_ET_wp[OBJECTARRAYSIZE];
		float jetGammaMax_EM_wp[OBJECTARRAYSIZE];
		float jetGammaMax_Hadronic_wp[OBJECTARRAYSIZE];
		float jetAlphaMax_wp[OBJECTARRAYSIZE];
		float jetBetaMax_wp[OBJECTARRAYSIZE];

		float jetPtAllTracks_wp[OBJECTARRAYSIZE];
		float jetPtAllPVTracks_wp[OBJECTARRAYSIZE];
		float jetMedianTheta2D_wp[OBJECTARRAYSIZE];
		float jetMedianIP_wp[OBJECTARRAYSIZE];
		float jetMinDeltaRAllTracks_wp[OBJECTARRAYSIZE];
		float jetMinDeltaRPVTracks_wp[OBJECTARRAYSIZE];

		float jetChargedEMEnergyFraction[OBJECTARRAYSIZE];
		float jetNeutralEMEnergyFraction[OBJECTARRAYSIZE];

		int  jetNPFCands[OBJECTARRAYSIZE];
		int  jetPFCandIndex[OBJECTARRAYSIZE][MAX_NPFCAND];
		int  jetAllPFCandIndex[OBJECTARRAYSIZE][MAX_NPFCAND];

		UInt_t          nFatJets;
		Float_t         fatJetE[400];   //[nFatJets]
		Float_t         fatJetPt[400];   //[nFatJets]
		Float_t         fatJetEta[400];   //[nFatJets]
		Float_t         fatJetPhi[400];   //[nFatJets]
		Float_t         fatJetCorrectedPt[400];   //[nFatJets]
		Float_t         fatJetPrunedM[400];   //[nFatJets]
		Float_t         fatJetTrimmedM[400];   //[nFatJets]
		Float_t         fatJetFilteredM[400];   //[nFatJets]
		Float_t         fatJetSoftDropM[400];   //[nFatJets]
		Float_t         fatJetCorrectedSoftDropM[400];   //[nFatJets]
		Float_t         fatJetUncorrectedSoftDropM[400];   //[nFatJets]
		Float_t         fatJetTau1[400];   //[nFatJets]
		Float_t         fatJetTau2[400];   //[nFatJets]
		Float_t         fatJetTau3[400];   //[nFatJets]
		Float_t         fatJetMaxSubjetCSV[400];   //[nFatJets]
		Bool_t          fatJetPassIDLoose[400];   //[nFatJets]
		Bool_t          fatJetPassIDTight[400];   //[nFatJets]
		Float_t         metPt;
		Float_t         metPhi;
		Float_t         sumMET;
		Float_t         metType0Pt;
		Float_t         metType0Phi;
		Float_t         metType1Pt_raw;
		Float_t         metType1Pt;
		Float_t         metType1Px;
		Float_t         metType1Py;
		Float_t         metType1Eta;
		Float_t         metType1Phi;
		Float_t         metType1Phi_raw;
		Float_t         metType0Plus1Pt;
		Float_t         metType0Plus1Phi;
		Float_t         metNoHFPt;
		Float_t         metNoHFPhi;
		Float_t         metPuppiPt;
		Float_t         metPuppiPhi;
		Float_t         metCaloPt;
		Float_t         metCaloPhi;
		Bool_t          Flag_HBHENoiseFilter;
		Bool_t          Flag_HBHETightNoiseFilter;
		Bool_t          Flag_HBHEIsoNoiseFilter;
		Bool_t          Flag_badChargedCandidateFilter;
		Bool_t          Flag_badMuonFilter;
		Bool_t          Flag_badGlobalMuonFilter;
		Bool_t          Flag_duplicateMuonFilter;
		Bool_t          Flag_CSCTightHaloFilter;
		Bool_t          Flag_hcalLaserEventFilter;
		Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
		Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
		Bool_t          Flag_goodVertices;
		Bool_t          Flag_trackingFailureFilter;
		Bool_t          Flag_eeBadScFilter;
		Bool_t          Flag_ecalLaserCorrFilter;
		Bool_t          Flag_trkPOGFilters;
		Bool_t          Flag_trkPOG_manystripclus53X;
		Bool_t          Flag_trkPOG_toomanystripclus53X;
		Bool_t          Flag_trkPOG_logErrorTooManyClusters;
		Bool_t          Flag_BadPFMuonFilter;
		Bool_t          Flag_BadChargedCandidateFilter;
		Bool_t          Flag_ecalBadCalibFilter;
		Bool_t          Flag_METFilters;
		Bool_t	   Flag2_globalSuperTightHalo2016Filter;
		Bool_t	   Flag2_globalTightHalo2016Filter;
		Bool_t	   Flag2_goodVertices;
		Bool_t	   Flag2_BadChargedCandidateFilter;
		Bool_t	   Flag2_BadPFMuonFilter;
		Bool_t	   Flag2_EcalDeadCellTriggerPrimitiveFilter;
		Bool_t	   Flag2_HBHENoiseFilter;
		Bool_t	   Flag2_HBHEIsoNoiseFilter;
		Bool_t	   Flag2_ecalBadCalibFilter;
		Bool_t	   Flag2_eeBadScFilter;
		Bool_t          HLTDecision[1201];
		Int_t           HLTPrescale[1201];
		Int_t           nGenJets;
		Float_t         genJetE[500];   //[nGenJets]
		Float_t         genJetPt[500];   //[nGenJets]
		Float_t         genJetEta[500];   //[nGenJets]
		Float_t         genJetPhi[500];   //[nGenJets]
		Float_t         genMetPtCalo;
		Float_t         genMetPhiCalo;
		Float_t         genMetPtTrue;
		Float_t         genMetPhiTrue;
		Float_t         genVertexX;
		Float_t         genVertexY;
		Float_t         genVertexZ;
		Float_t         genVertexT;
		Float_t         genWeight;
		UInt_t          genSignalProcessID;
		Float_t         genQScale;
		Float_t         genAlphaQCD;
		Float_t         genAlphaQED;
		std::string		*lheComments;
		std::vector<float>   *scaleWeights;
		std::vector<float>   *pdfWeights;
		std::vector<float>   *alphasWeights;
		Int_t           nGenParticle;
		Int_t           gParticleMotherId[4000];   //[nGenParticle]
		Int_t           gParticleMotherIndex[4000];   //[nGenParticle]
		Int_t           gParticleId[4000];   //[nGenParticle]
		Int_t           gParticleStatus[4000];   //[nGenParticle]
		Float_t         gParticleE[4000];   //[nGenParticle]
		Float_t         gParticlePt[4000];   //[nGenParticle]
		Float_t         gParticleEta[4000];   //[nGenParticle]
		Float_t         gParticlePhi[4000];   //[nGenParticle]

		//gLLP
		Float_t         gLLP_travel_time[2];
		Float_t         gLLP_e[2];
		Float_t         gLLP_pt[2];
		Float_t         gLLP_eta[2];
		Float_t         gLLP_beta[2];
		Float_t         gLLP_phi[2];
		Float_t         gLLP_decay_vertex_x[2];
		Float_t         gLLP_decay_vertex_y[2];
		Float_t         gLLP_decay_vertex_z[2];
		Float_t         gLLP_prod_vertex_x[2];
		Float_t         gLLP_prod_vertex_y[2];
		Float_t         gLLP_prod_vertex_z[2];
		/*
		//gLLP_daughter
		Float_t         gen_time[4];
		Float_t         photon_travel_time[4];
		Float_t         gLLP_daughter_travel_time[4];
		Float_t         gLLP_daughter_e[4];
		Float_t         gLLP_daughter_pt[4];
		Float_t         gLLP_daughter_eta[4];
		Float_t         gLLP_daughter_phi[4];
		Float_t         gLLP_daughter_eta_ecalcorr[4];
		Float_t         gLLP_daughter_phi_ecalcorr[4];
		Float_t         gLLP_min_delta_r_match_jet[4];
		UInt_t          gLLP_daughter_match_jet_index[4];
		*/
		//daughters
		bool gLLP_daughter_EB[4]; 
		bool gLLP_daughter_ETL[4];

		float gLLP_daughter_photon_travel_time_EB[4];
		float gLLP_daughter_photon_travel_time_ETL[4];

		float gLLP_daughter_travel_time_EB[4];
		float gLLP_daughter_travel_time_ETL[4];

		float gen_time_daughter_EB[4];
		float gen_time_daughter_ETL[4];

		int   gLLP_daughter_id[4];
		float gLLP_daughter_pt[4];
		float gLLP_daughter_eta[4];
		float gLLP_daughter_phi[4];
		float gLLP_daughter_eta_ecalcorr[4];
		float gLLP_daughter_phi_ecalcorr[4];
		float gLLP_daughter_e[4];
		float gLLP_daughter_mass[4];

		unsigned int gLLP_daughter_match_jet_index[4];
		float gLLP_daughter_min_delta_r_match_jet[4];

		//grandaughters
		bool gLLP_grandaughter_EB[4]; 
		bool gLLP_grandaughter_ETL[4];

		float gLLP_grandaughter_photon_travel_time_EB[4];
		float gLLP_grandaughter_photon_travel_time_ETL[4];

		float gLLP_grandaughter_travel_time_EB[4];
		float gLLP_grandaughter_travel_time_ETL[4];

		float gen_time_grandaughter_EB[4];
		float gen_time_grandaughter_ETL[4];

		int   gLLP_grandaughter_id[4];
		float gLLP_grandaughter_pt[4];
		float gLLP_grandaughter_eta[4];
		float gLLP_grandaughter_phi[4];
		float gLLP_grandaughter_eta_ecalcorr[4];
		float gLLP_grandaughter_phi_ecalcorr[4];
		float gLLP_grandaughter_e[4];
		float gLLP_grandaughter_mass[4];

		unsigned int gLLP_grandaughter_match_jet_index[4];
		float gLLP_grandaughter_min_delta_r_match_jet[4];

		//Tracks
		int nTracks;
		float track_Pt[RECHITARRAYSIZE];
		float track_Eta[RECHITARRAYSIZE];
		float track_Phi[RECHITARRAYSIZE];
		int   track_charge[RECHITARRAYSIZE];
		int   track_bestVertexIndex[RECHITARRAYSIZE];
		int   track_nMissingInnerHits[RECHITARRAYSIZE];
		int   track_nMissingOuterHits[RECHITARRAYSIZE];
		int   track_nPixelHits[RECHITARRAYSIZE];
		int   track_nHits[RECHITARRAYSIZE]; 
		float track_angle[RECHITARRAYSIZE];
		float track_dxyToBS[RECHITARRAYSIZE];
		float track_dxyErr[RECHITARRAYSIZE];
		float track_dzToPV[RECHITARRAYSIZE];
		float track_dzErr[RECHITARRAYSIZE];
		float track_chi2[RECHITARRAYSIZE];
		int   track_ndof[RECHITARRAYSIZE];

		// List of branches
		TBranch        *b_isData;   //!
		TBranch        *b_nPV;   //!
		TBranch        *b_runNum;   //!
		TBranch        *b_nSlimmedSecondV;   //!
		TBranch        *b_lumiNum;   //!
		TBranch        *b_eventNum;   //!
		TBranch        *b_eventTime;   //!
		TBranch        *b_pvX;   //!
		TBranch        *b_pvY;   //!
		TBranch        *b_pvZ;   //!
		TBranch        *b_fixedGridRhoAll;   //!
		TBranch        *b_fixedGridRhoFastjetAll;   //!
		TBranch        *b_fixedGridRhoFastjetAllCalo;   //!
		TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
		TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!
		TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
		TBranch        *b_nPVAll;   //!
		TBranch        *b_pvAllX;   //!
		TBranch        *b_pvAllY;   //!
		TBranch        *b_pvAllZ;   //!
		TBranch        *b_pvAllLogSumPtSq;   //!
		TBranch        *b_pvAllSumPx;   //!
		TBranch        *b_pvAllSumPy;   //!
		TBranch        *b_nBunchXing;   //!
		TBranch        *b_BunchXing;   //!
		TBranch        *b_nPU;   //!
		TBranch        *b_nPUmean;   //!
		TBranch        *b_nMuons;   //!
		TBranch        *b_muonE;   //!
		TBranch        *b_muonPt;   //!
		TBranch        *b_muonEta;   //!
		TBranch        *b_muonPhi;   //!
		TBranch        *b_muonCharge;   //!
		TBranch        *b_muonIsLoose;   //!
		TBranch        *b_muonIsMedium;   //!
		TBranch        *b_muonIsTight;   //!
		TBranch        *b_muon_d0;   //!
		TBranch        *b_muon_dZ;   //!
		TBranch        *b_muon_ip3d;   //!
		TBranch        *b_muon_ip3dSignificance;   //!
		TBranch        *b_muonType;   //!
		TBranch        *b_muonQuality;   //!
		TBranch        *b_muon_pileupIso;   //!
		TBranch        *b_muon_chargedIso;   //!
		TBranch        *b_muon_photonIso;   //!
		TBranch        *b_muon_neutralHadIso;   //!
		TBranch        *b_muon_ptrel;   //!
		TBranch        *b_muon_chargedMiniIso;   //!
		TBranch        *b_muon_photonAndNeutralHadronMiniIso;   //!
		TBranch        *b_muon_chargedPileupMiniIso;   //!
		TBranch        *b_muon_activityMiniIsoAnnulus;   //!
		TBranch        *b_muon_passSingleMuTagFilter;   //!
		TBranch        *b_muon_passHLTFilter;   //!
		TBranch        *b_muon_validFractionTrackerHits;   //!
		TBranch        *b_muon_isGlobal;   //!
		TBranch        *b_muon_normChi2;   //!
		TBranch        *b_muon_chi2LocalPosition;   //!
		TBranch        *b_muon_kinkFinder;   //!
		TBranch        *b_muon_segmentCompatability;   //!
		TBranch        *b_muonIsICHEPMedium;   //!
		TBranch        *b_nElectrons;   //!
		TBranch        *b_eleE;   //!
		TBranch        *b_elePt;   //!
		TBranch        *b_eleEta;   //!
		TBranch        *b_elePhi;   //!
		TBranch        *b_eleCharge;   //!
		TBranch        *b_eleEta_SC;   //!
		TBranch        *b_eleSigmaIetaIeta;   //!
		TBranch        *b_eleFull5x5SigmaIetaIeta;   //!
		TBranch        *b_eleR9;   //!
		TBranch        *b_ele_dEta;   //!
		TBranch        *b_ele_dPhi;   //!
		TBranch        *b_ele_HoverE;   //!
		TBranch        *b_ele_d0;   //!
		TBranch        *b_ele_dZ;   //!
		TBranch        *b_ele_ip3d;   //!
		TBranch        *b_ele_ip3dSignificance;   //!
		TBranch        *b_ele_pileupIso;   //!
		TBranch        *b_ele_chargedIso;   //!
		TBranch        *b_ele_photonIso;   //!
		TBranch        *b_ele_neutralHadIso;   //!
		TBranch        *b_ele_MissHits;   //!
		TBranch        *b_ele_PassConvVeto;   //!
		TBranch        *b_ele_OneOverEminusOneOverP;   //!
		TBranch        *b_ele_IDMVAGeneralPurpose;   //!
		TBranch        *b_ele_IDMVACategoryGeneralPurpose;   //!
		TBranch        *b_ele_IDMVAHZZ;   //!
		TBranch        *b_ele_IDMVACategoryHZZ;   //!
		TBranch        *b_ele_RegressionE;   //!
		TBranch        *b_ele_CombineP4;   //!
		TBranch        *b_ele_ptrel;   //!
		TBranch        *b_ele_chargedMiniIso;   //!
		TBranch        *b_ele_photonAndNeutralHadronMiniIso;   //!
		TBranch        *b_ele_chargedPileupMiniIso;   //!
		TBranch        *b_ele_activityMiniIsoAnnulus;   //!
		TBranch 	  *b_ele_passCutBasedIDVetog;
		TBranch 	  *b_ele_passCutBasedIDLooseg;
		TBranch 	  *b_ele_passCutBasedIDMediumg;
		TBranch 	  *b_ele_passCutBasedIDTightg;
		TBranch 	  *b_ele_passMVAIsoIDWP80g;
		TBranch 	  *b_ele_passMVAIsoIDWP90g;
		TBranch 	  *b_ele_passMVAIsoIDWPHZZg;
		TBranch 	  *b_ele_passMVAIsoIDWPLooseg;
		TBranch 	  *b_ele_passMVANoIsoIDWP80g;
		TBranch 	  *b_ele_passMVANoIsoIDWP90g;
		TBranch 	  *b_ele_passMVANoIsoIDWPLooseg;
		TBranch        *b_ele_passSingleEleTagFilter;   //!
		TBranch        *b_ele_passTPOneTagFilter;   //!
		TBranch        *b_ele_passTPTwoTagFilter;   //!
		TBranch        *b_ele_passTPOneProbeFilter;   //!
		TBranch        *b_ele_passTPTwoProbeFilter;   //!
		TBranch        *b_ele_passHLTFilter;   //!
		TBranch        *b_nTaus;   //!
		TBranch        *b_tauE;   //!
		TBranch        *b_tauPt;   //!
		TBranch        *b_tauEta;   //!
		TBranch        *b_tauPhi;   //!
		TBranch        *b_tau_IsLoose;   //!
		TBranch        *b_tau_IsMedium;   //!
		TBranch        *b_tau_IsTight;   //!
		TBranch        *b_tau_passEleVetoLoose;   //!
		TBranch        *b_tau_passEleVetoMedium;   //!
		TBranch        *b_tau_passEleVetoTight;   //!
		TBranch        *b_tau_passMuVetoLoose;   //!
		TBranch        *b_tau_passMuVetoMedium;   //!
		TBranch        *b_tau_passMuVetoTight;   //!
		TBranch        *b_tau_ID;   //!
		TBranch        *b_tau_combinedIsoDeltaBetaCorr3Hits;   //!
		TBranch        *b_tau_chargedIsoPtSum;   //!
		TBranch        *b_tau_neutralIsoPtSum;   //!
		TBranch        *b_tau_puCorrPtSum;   //!
		TBranch        *b_tau_eleVetoMVA;   //!
		TBranch        *b_tau_eleVetoCategory;   //!
		TBranch        *b_tau_muonVetoMVA;   //!
		TBranch        *b_tau_isoMVAnewDMwLT;   //!
		TBranch        *b_tau_isoMVAnewDMwoLT;   //!
		TBranch        *b_tau_leadCandPt;   //!
		TBranch        *b_tau_leadCandID;   //!
		TBranch        *b_tau_leadChargedHadrCandPt;   //!
		TBranch        *b_tau_leadChargedHadrCandID;   //!
		TBranch        *b_nIsoPFCandidates;   //!
		TBranch        *b_isoPFCandidatePt;   //!
		TBranch        *b_isoPFCandidateEta;   //!
		TBranch        *b_isoPFCandidatePhi;   //!
		TBranch        *b_isoPFCandidateIso04;   //!
		TBranch        *b_isoPFCandidateD0;   //!
		TBranch        *b_isoPFCandidatePdgId;   //!
		TBranch        *b_nPhotons;   //!
		TBranch        *b_nPhotons_overlap;   //!
		TBranch        *b_phoE;   //!
		TBranch        *b_phoPt;   //!
		TBranch        *b_phoEta;   //!
		TBranch        *b_phoPhi;   //!
		TBranch        *b_phoSigmaIetaIeta;   //!
		TBranch        *b_phoFull5x5SigmaIetaIeta;   //!
		TBranch        *b_phoR9;   //!
		TBranch        *b_pho_sminor;   //!
		TBranch        *b_pho_smajor;   //!
		TBranch        *b_pho_HoverE;   //!
		TBranch        *b_pho_sumChargedHadronPtAllVertices;   //!
		TBranch        *b_pho_sumChargedHadronPt;   //!
		TBranch        *b_pho_sumNeutralHadronEt;   //!
		TBranch        *b_pho_sumPhotonEt;   //!
		TBranch        *b_pho_ecalPFClusterIso;   //!
		TBranch        *b_pho_hcalPFClusterIso;   //!
		TBranch        *b_pho_trkSumPtHollowConeDR03;   //!
		TBranch        *b_pho_sumWorstVertexChargedHadronPt;   //!
		TBranch        *b_pho_pfIsoChargedHadronIso;   //!
		TBranch        *b_pho_pfIsoChargedHadronIsoWrongVtx;   //!
		TBranch        *b_pho_pfIsoNeutralHadronIso;   //!
		TBranch        *b_pho_pfIsoPhotonIso;   //!
		TBranch        *b_pho_pfIsoModFrixione;   //!
		TBranch        *b_pho_pfIsoSumPUPt;   //!
		TBranch        *b_pho_isConversion;   //!
		TBranch        *b_pho_passEleVeto;   //!
		TBranch        *b_pho_RegressionE;   //!
		TBranch        *b_pho_RegressionEUncertainty;   //!
		TBranch        *b_pho_IDMVA;   //!
		TBranch        *b_pho_superClusterEnergy;   //!
		TBranch        *b_pho_superClusterRawEnergy;   //!
		TBranch        *b_pho_superClusterEta;   //!
		TBranch        *b_pho_superClusterPhi;   //!
		TBranch        *b_pho_superClusterX;   //!
		TBranch        *b_pho_superClusterY;   //!
		TBranch        *b_pho_superClusterZ;   //!
		TBranch        *b_pho_hasPixelSeed;   //!
		TBranch        *b_pho_passHLTFilter;   //!
		TBranch	  *b_pho_passCutBasedIDLoose;
		TBranch	  *b_pho_passCutBasedIDMedium;
		TBranch	  *b_pho_passCutBasedIDTight;
		TBranch	  *b_pho_passMVAIDWP80;
		TBranch	  *b_pho_passMVAIDWP90;
		TBranch        *b_pho_convType;   //!
		TBranch        *b_pho_convTrkZ;   //!
		TBranch        *b_pho_convTrkClusZ;   //!
		TBranch        *b_pho_vtxSumPx;   //!
		TBranch        *b_pho_vtxSumPy;   //!
		TBranch        *b_pho_isStandardPhoton;   //!
		TBranch        *b_pho_seedRecHitSwitchToGain6;   //!
		TBranch        *b_pho_seedRecHitSwitchToGain1;   //!
		TBranch        *b_pho_anyRecHitSwitchToGain6;   //!
		TBranch        *b_pho_anyRecHitSwitchToGain1;   //!
		TBranch        *b_ecalRechit_Eta;   //!
		TBranch        *b_ecalRechit_Phi;   //!
		TBranch        *b_ecalRechit_X;   //!
		TBranch        *b_ecalRechit_Y;   //!
		TBranch        *b_ecalRechit_Z;   //!
		TBranch        *b_ecalRechit_E;   //!
		TBranch        *b_ecalRechit_T;   //!
		TBranch        *b_ecalRechit_ID;   //!
		TBranch        *b_ecalRechit_FlagOOT;   //!
		TBranch        *b_ecalRechit_GainSwitch1;   //!
		TBranch        *b_ecalRechit_GainSwitch6;   //!
		TBranch        *b_ecalRechit_transpCorr;   //!
		TBranch        *b_nCsc;   //!
		TBranch        *b_cscPhi;   //!
		TBranch        *b_cscEta;   //!
		TBranch        *b_cscX;   //!
		TBranch        *b_cscY;   //!
		TBranch        *b_cscZ;   //!
		TBranch        *b_cscDirectionX;   //!
		TBranch        *b_cscDirectionY;   //!
		TBranch        *b_cscDirectionZ;   //!
		TBranch        *b_cscNRecHits;   //!
		TBranch        *b_cscNRecHits_flag;   //!
		TBranch        *b_cscT;   //!
		TBranch        *b_cscChi2;   //!
		TBranch        *b_nRpc;   //!
		TBranch        *b_rpcPhi;   //!
		TBranch        *b_rpcEta;   //!
		TBranch        *b_rpcX;   //!
		TBranch        *b_rpcY;   //!
		TBranch        *b_rpcZ;   //!
		TBranch        *b_rpcT;   //!
		TBranch        *b_rpcTError;   //!
		TBranch        *b_nDt;   //!
		TBranch        *b_dtPhi;   //!
		TBranch        *b_dtEta;   //!
		TBranch        *b_dtX;   //!
		TBranch        *b_dtY;   //!
		TBranch        *b_dtZ;   //!
		TBranch        *b_dtDirX;   //!
		TBranch        *b_dtDirY;   //!
		TBranch        *b_dtDirZ;   //!
		TBranch        *b_dtT;   //!
		TBranch        *b_dtTError;   //!
		TBranch        *b_nCaloJets;   //!
		TBranch        *b_nJets;
		TBranch        *b_jetE;
		TBranch        *b_jetPt;
		TBranch        *b_jetEta;
		TBranch        *b_jetPhi;
		TBranch        *b_jetCSV;
		TBranch        *b_jetCISV;
		TBranch        *b_jetProbb;
		TBranch        *b_jetProbc;
		TBranch        *b_jetProbudsg;
		TBranch        *b_jetProbbb;
		TBranch        *b_jetMass;
		TBranch        *b_jetJetArea;
		TBranch        *b_jetPileupE;
		TBranch        *b_jetPileupId;
		TBranch        *b_jetPileupIdFlag;
		TBranch        *b_jetPassIDLoose;
		TBranch        *b_jetPassIDTight;
		TBranch        *b_jetPassMuFrac;
		TBranch        *b_jetPassEleFrac;
		TBranch        *b_jetPartonFlavor;
		TBranch        *b_jetHadronFlavor;
		TBranch        *b_jetElectronEnergyFraction;
		TBranch        *b_jetPhotonEnergyFraction;
		TBranch        *b_jetChargedHadronEnergyFraction;
		TBranch        *b_jetNeutralHadronEnergyFraction;
		TBranch        *b_jetMuonEnergyFraction;
		TBranch        *b_jetHOEnergyFraction;
		TBranch        *b_jetHFHadronEnergyFraction;
		TBranch        *b_jetHFEMEnergyFraction;
		TBranch        *b_jetChargedHadronMultiplicity;
		TBranch        *b_jetNeutralHadronMultiplicity;
		TBranch        *b_jetElectronMultiplicity;
		TBranch        *b_jetPhotonMultiplicity;
		TBranch        *b_jetMuonMultiplicity;
		TBranch        *b_jetAllMuonPt;
		TBranch        *b_jetAllMuonEta;
		TBranch        *b_jetAllMuonPhi;
		TBranch        *b_jetAllMuonM;
		TBranch        *b_jetPtWeightedDZ;
		TBranch        *b_jetNRechits;
		TBranch        *b_jetRechitE;
		TBranch        *b_jetRechitT;
		TBranch        *b_jetRechitT_rms;
		TBranch        *b_jetRechitE_Error;
		TBranch        *b_jetRechitT_Error;

		TBranch        *b_jetGammaMax;
		TBranch        *b_jetGammaMax_ET;
		TBranch        *b_jetGammaMax_EM;
		TBranch        *b_jetGammaMax_Hadronic;
		TBranch        *b_jetAlphaMax;
		TBranch        *b_jetBetaMax;

		TBranch        *b_jetPtAllTracks;
		TBranch        *b_jetPtAllPVTracks;
		TBranch        *b_jetMedianTheta2D;
		TBranch        *b_jetMedianIP;
		TBranch        *b_jetMinDeltaRAllTracks;
		TBranch        *b_jetMinDeltaRPVTracks;

		TBranch        *b_jet_energy_frac;
		TBranch        *b_jet_sig_et1;
		TBranch        *b_jet_sig_et2;
		TBranch        *b_jet_matched;
		TBranch        *b_jet_matched_gLLP0_grandaughter;
		TBranch        *b_jet_matched_gLLP1_grandaughter;

		TBranch        *b_jetGammaMax_wp;
		TBranch        *b_jetGammaMax_ET_wp;
		TBranch        *b_jetGammaMax_EM_wp;
		TBranch        *b_jetGammaMax_Hadronic_wp;
		TBranch        *b_jetAlphaMax_wp;
		TBranch        *b_jetBetaMax_wp;

		TBranch        *b_jetPtAllTracks_wp;
		TBranch        *b_jetPtAllPVTracks_wp;
		TBranch        *b_jetMedianTheta2D_wp;
		TBranch        *b_jetMedianIP_wp;
		TBranch        *b_jetMinDeltaRAllTracks_wp;
		TBranch        *b_jetMinDeltaRPVTracks_wp;

		TBranch        *b_jetChargedEMEnergyFraction;   //!
		TBranch        *b_jetNeutralEMEnergyFraction;   //!
		TBranch        *b_nFatJets;   //!
		TBranch        *b_fatJetE;   //!
		TBranch        *b_fatJetPt;   //!
		TBranch        *b_fatJetEta;   //!
		TBranch        *b_fatJetPhi;   //!
		TBranch        *b_fatJetCorrectedPt;   //!
		TBranch        *b_fatJetPrunedM;   //!
		TBranch        *b_fatJetTrimmedM;   //!
		TBranch        *b_fatJetFilteredM;   //!
		TBranch        *b_fatJetSoftDropM;   //!
		TBranch        *b_fatJetCorrectedSoftDropM;   //!
		TBranch        *b_fatJetUncorrectedSoftDropM;   //!
		TBranch        *b_fatJetTau1;   //!
		TBranch        *b_fatJetTau2;   //!
		TBranch        *b_fatJetTau3;   //!
		TBranch        *b_fatJetMaxSubjetCSV;   //!
		TBranch        *b_fatJetPassIDLoose;   //!
		TBranch        *b_fatJetPassIDTight;   //!
		TBranch        *b_metPt;   //!
		TBranch        *b_metPhi;   //!
		TBranch        *b_sumMET;   //!
		TBranch        *b_metType0Pt;   //!
		TBranch        *b_metType0Phi;   //!
		TBranch        *b_metType1Pt_raw;   //!
		TBranch        *b_metType1Pt;   //!
		TBranch        *b_metType1Px;   //!
		TBranch        *b_metType1Py;   //!
		TBranch        *b_metType1Eta;   //!
		TBranch        *b_metType1Phi;   //!
		TBranch        *b_metType1Phi_raw;   //!
		TBranch        *b_metType0Plus1Pt;   //!
		TBranch        *b_metType0Plus1Phi;   //!
		TBranch        *b_metNoHFPt;   //!
		TBranch        *b_metNoHFPhi;   //!
		TBranch        *b_metPuppiPt;   //!
		TBranch        *b_metPuppiPhi;   //!
		TBranch        *b_metCaloPt;   //!
		TBranch        *b_metCaloPhi;   //!
		TBranch        *b_Flag_HBHENoiseFilter;   //!
		TBranch        *b_Flag_HBHETightNoiseFilter;   //!
		TBranch        *b_Flag_HBHEIsoNoiseFilter;   //!
		TBranch        *b_Flag_badChargedCandidateFilter;   //!
		TBranch        *b_Flag_badMuonFilter;   //!
		TBranch        *b_Flag_badGlobalMuonFilter;   //!
		TBranch        *b_Flag_duplicateMuonFilter;   //!
		TBranch        *b_Flag_CSCTightHaloFilter;   //!
		TBranch        *b_Flag_hcalLaserEventFilter;   //!
		TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
		TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
		TBranch        *b_Flag_goodVertices;   //!
		TBranch        *b_Flag_trackingFailureFilter;   //!
		TBranch        *b_Flag_eeBadScFilter;   //!
		TBranch        *b_Flag_ecalLaserCorrFilter;   //!
		TBranch        *b_Flag_trkPOGFilters;   //!
		TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
		TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
		TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
		TBranch        *b_Flag_BadPFMuonFilter;   //!
		TBranch        *b_Flag_BadChargedCandidateFilter;   //!
		TBranch        *b_Flag_ecalBadCalibFilter;   //!
		TBranch        *b_Flag_METFilters;   //!
		TBranch        *b_Flag2_globalSuperTightHalo2016Filter;
		TBranch        *b_Flag2_globalTightHalo2016Filter;
		TBranch        *b_Flag2_goodVertices;
		TBranch        *b_Flag2_BadChargedCandidateFilter;
		TBranch        *b_Flag2_BadPFMuonFilter;
		TBranch        *b_Flag2_EcalDeadCellTriggerPrimitiveFilter;
		TBranch        *b_Flag2_HBHENoiseFilter;
		TBranch        *b_Flag2_HBHEIsoNoiseFilter;
		TBranch        *b_Flag2_ecalBadCalibFilter;
		TBranch        *b_Flag2_eeBadScFilter;
		TBranch        *b_HLTDecision;   //!
		TBranch        *b_HLTPrescale;   //!
		TBranch        *b_nGenJets;   //!
		TBranch        *b_genJetE;   //!
		TBranch        *b_genJetPt;   //!
		TBranch        *b_genJetEta;   //!
		TBranch        *b_genJetPhi;   //!
		TBranch        *b_genMetPtCalo;   //!
		TBranch        *b_genMetPhiCalo;   //!
		TBranch        *b_genMetPtTrue;   //!
		TBranch        *b_genMetPhiTrue;   //!
		TBranch        *b_genVertexX;   //!
		TBranch        *b_genVertexY;   //!
		TBranch        *b_genVertexZ;   //!
		TBranch        *b_genVertexT;   //!
		TBranch        *b_genWeight;   //!
		TBranch        *b_genSignalProcessID;   //!
		TBranch        *b_genQScale;   //!
		TBranch        *b_genAlphaQCD;   //!
		TBranch        *b_genAlphaQED;   //!
		TBranch        *b_lheComments;   //!
		TBranch        *b_scaleWeights;   //!
		TBranch        *b_pdfWeights;   //!
		TBranch        *b_alphasWeights;   //!
		TBranch        *b_nGenParticle;   //!
		TBranch        *b_gParticleMotherId;   //!
		TBranch        *b_gParticleMotherIndex;   //!
		TBranch        *b_gParticleId;   //!
		TBranch        *b_gParticleStatus;   //!
		TBranch        *b_gParticleE;   //!
		TBranch        *b_gParticlePt;   //!
		TBranch        *b_gParticleEta;   //!
		TBranch        *b_gParticlePhi;   //!
		TBranch        *b_gLLP_travel_time;
		TBranch        *b_gLLP_e;
		TBranch        *b_gLLP_pt;
		TBranch        *b_gLLP_eta;
		TBranch        *b_gLLP_beta;
		TBranch        *b_gLLP_phi;
		TBranch        *b_gLLP_decay_vertex_x;
		TBranch        *b_gLLP_decay_vertex_y;
		TBranch        *b_gLLP_decay_vertex_z;
		TBranch        *b_gLLP_prod_vertex_x;
		TBranch        *b_gLLP_prod_vertex_y;
		TBranch        *b_gLLP_prod_vertex_z;
		/*
		   TBranch        *b_gen_time;
		   TBranch        *b_photon_travel_time;
		   TBranch        *b_gLLP_daughter_travel_time;
		   TBranch        *b_gLLP_daughter_e;
		   TBranch        *b_gLLP_daughter_pt;
		   TBranch        *b_gLLP_daughter_eta;
		   TBranch        *b_gLLP_daughter_phi;
		   TBranch        *b_gLLP_daughter_eta_ecalcorr;
		   TBranch        *b_gLLP_daughter_phi_ecalcorr;
		   TBranch        *b_gLLP_min_delta_r_match_jet;
		   TBranch        *b_gLLP_daughter_match_jet_index;
		   */
		//daughters
		TBranch *b_gLLP_daughter_EB; 
		TBranch *b_gLLP_daughter_ETL;

		TBranch *b_gLLP_daughter_photon_travel_time_EB;
		TBranch *b_gLLP_daughter_photon_travel_time_ETL;

		TBranch *b_gLLP_daughter_travel_time_EB;
		TBranch *b_gLLP_daughter_travel_time_ETL;

		TBranch *b_gen_time_daughter_EB;
		TBranch *b_gen_time_daughter_ETL;

		TBranch *b_gLLP_daughter_id;
		TBranch *b_gLLP_daughter_pt;
		TBranch *b_gLLP_daughter_eta;
		TBranch *b_gLLP_daughter_phi;
		TBranch *b_gLLP_daughter_eta_ecalcorr;
		TBranch *b_gLLP_daughter_phi_ecalcorr;
		TBranch *b_gLLP_daughter_e;
		TBranch *b_gLLP_daughter_mass;

		TBranch *b_gLLP_daughter_match_jet_index;
		TBranch *b_gLLP_daughter_min_delta_r_match_jet;

		//grandaughters
		TBranch *b_gLLP_grandaughter_EB; 
		TBranch *b_gLLP_grandaughter_ETL;

		TBranch *b_gLLP_grandaughter_photon_travel_time_EB;
		TBranch *b_gLLP_grandaughter_photon_travel_time_ETL;

		TBranch *b_gLLP_grandaughter_travel_time_EB;
		TBranch *b_gLLP_grandaughter_travel_time_ETL;

		TBranch *b_gen_time_grandaughter_EB;
		TBranch *b_gen_time_grandaughter_ETL;

		TBranch *b_gLLP_grandaughter_id;
		TBranch *b_gLLP_grandaughter_pt;
		TBranch *b_gLLP_grandaughter_eta;
		TBranch *b_gLLP_grandaughter_phi;
		TBranch *b_gLLP_grandaughter_eta_ecalcorr;
		TBranch *b_gLLP_grandaughter_phi_ecalcorr;
		TBranch *b_gLLP_grandaughter_e;
		TBranch *b_gLLP_grandaughter_mass;

		TBranch *b_gLLP_grandaughter_match_jet_index;
		TBranch *b_gLLP_grandaughter_min_delta_r_match_jet;

		TBranch *b_nTracks;
		TBranch *b_track_Pt;
		TBranch *b_track_Eta;
		TBranch *b_track_Phi;
		TBranch *b_track_charge;
		TBranch *b_track_bestVertexIndex;
		TBranch *b_track_nMissingInnerHits;
		TBranch *b_track_nMissingOuterHits;
		TBranch *b_track_nPixelHits;
		TBranch *b_track_nHits; 
		TBranch *b_track_angle;
		TBranch *b_track_dxyToBS;
		TBranch *b_track_dxyErr;
		TBranch *b_track_dzToPV;
		TBranch *b_track_dzErr;
		TBranch *b_track_chi2;
		TBranch *b_track_ndof;

		llp_event(TTree *tree=0);
		virtual ~llp_event();
		//virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		//virtual void     Loop();
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef llp_event_cxx
/* llp_event::llp_event(TTree *tree) : fChain(0) */
/* { */
/* 	// if parameter tree is not specified (or zero), connect the file */
/* 	// used to generate this class and read the Tree. */
/* 	if (tree == 0) { */
/* 		TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("llp_ntupler_126.root"); */
/* 		if (!f || !f->IsOpen()) { */
/* 			f = new TFile("llp_ntupler_126.root"); */
/* 		} */
/* 		TDirectory * dir = (TDirectory*)f->Get("llp_ntupler_126.root:/ntuples"); */
/* 		dir->GetObject("llp",tree); */

/* 	} */
/* 	Init(tree); */
/* } */


/*
   Int_t llp_event::Cut(L)
   {
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
return 1;
}
*/
#endif // #ifdef llp_event_cxx
