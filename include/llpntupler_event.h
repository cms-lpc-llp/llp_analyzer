//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 14 14:17:17 2019 by ROOT version 6.10/09
// from TTree llp/selected AOD information for llp analyses
// found on file: /mnt/hadoop/store/group/phys_exotica/delayedjets/llpntuple/V1p8/MC_Summer16/v2/sixie/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/Run2_LLPNtupler_V1p8_MC_Summer16_RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v2_v2_v2/191011_040441/0000/llp_ntupler_1.root
//////////////////////////////////////////////////////////

#ifndef llpntupler_event_h
#define llpntupler_event_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class llpntupler_event {
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
   Int_t           nBunchXing;
   Int_t           BunchXing[16];   //[nBunchXing]
   Int_t           nPU[16];   //[nBunchXing]
   Float_t         nPUmean[16];   //[nBunchXing]
   Int_t           nMuons;
   Float_t         muonE[9];   //[nMuons]
   Float_t         muonPt[9];   //[nMuons]
   Float_t         muonEta[9];   //[nMuons]
   Float_t         muonPhi[9];   //[nMuons]
   Int_t           muonCharge[9];   //[nMuons]
   Bool_t          muonIsLoose[9];   //[nMuons]
   Bool_t          muonIsMedium[9];   //[nMuons]
   Bool_t          muonIsTight[9];   //[nMuons]
   Float_t         muon_d0[9];   //[nMuons]
   Float_t         muon_dZ[9];   //[nMuons]
   Float_t         muon_ip3d[9];   //[nMuons]
   Float_t         muon_ip3dSignificance[9];   //[nMuons]
   UInt_t          muonType[9];   //[nMuons]
   UInt_t          muonQuality[9];   //[nMuons]
   Float_t         muon_pileupIso[9];   //[nMuons]
   Float_t         muon_chargedIso[9];   //[nMuons]
   Float_t         muon_photonIso[9];   //[nMuons]
   Float_t         muon_neutralHadIso[9];   //[nMuons]
   Float_t         muon_ptrel[9];   //[nMuons]
   Bool_t          muon_passSingleMuTagFilter[9];   //[nMuons]
   Bool_t          muon_passHLTFilter[9][100];   //[nMuons]
   Float_t         muon_validFractionTrackerHits[9];   //[nMuons]
   Bool_t          muon_isGlobal[9];   //[nMuons]
   Float_t         muon_normChi2[9];   //[nMuons]
   Float_t         muon_chi2LocalPosition[9];   //[nMuons]
   Float_t         muon_kinkFinder[9];   //[nMuons]
   Float_t         muon_segmentCompatability[9];   //[nMuons]
   Bool_t          muonIsICHEPMedium[9];   //[nMuons]
   Int_t           nElectrons;
   Float_t         eleE[5];   //[nElectrons]
   Float_t         elePt[5];   //[nElectrons]
   Float_t         eleEta[5];   //[nElectrons]
   Float_t         elePhi[5];   //[nElectrons]
   Float_t         eleCharge[5];   //[nElectrons]
   Float_t         eleEta_SC[5];   //[nElectrons]
   Float_t         eleSigmaIetaIeta[5];   //[nElectrons]
   Float_t         eleFull5x5SigmaIetaIeta[5];   //[nElectrons]
   Float_t         eleR9[5];   //[nElectrons]
   Float_t         ele_dEta[5];   //[nElectrons]
   Float_t         ele_dPhi[5];   //[nElectrons]
   Float_t         ele_HoverE[5];   //[nElectrons]
   Float_t         ele_d0[5];   //[nElectrons]
   Float_t         ele_dZ[5];   //[nElectrons]
   Float_t         ele_ip3d[5];   //[nElectrons]
   Float_t         ele_ip3dSignificance[5];   //[nElectrons]
   Float_t         ele_pileupIso[5];   //[nElectrons]
   Float_t         ele_chargedIso[5];   //[nElectrons]
   Float_t         ele_photonIso[5];   //[nElectrons]
   Float_t         ele_neutralHadIso[5];   //[nElectrons]
   Int_t           ele_MissHits[5];   //[nElectrons]
   Bool_t          ele_PassConvVeto[5];   //[nElectrons]
   Float_t         ele_OneOverEminusOneOverP[5];   //[nElectrons]
   Float_t         ele_IDMVAGeneralPurpose[5];   //[nElectrons]
   Int_t           ele_IDMVACategoryGeneralPurpose[5];   //[nElectrons]
   Float_t         ele_IDMVAHZZ[5];   //[nElectrons]
   Int_t           ele_IDMVACategoryHZZ[5];   //[nElectrons]
   Bool_t          ele_passSingleEleTagFilter[5];   //[nElectrons]
   Bool_t          ele_passTPOneTagFilter[5];   //[nElectrons]
   Bool_t          ele_passTPTwoTagFilter[5];   //[nElectrons]
   Bool_t          ele_passTPOneProbeFilter[5];   //[nElectrons]
   Bool_t          ele_passTPTwoProbeFilter[5];   //[nElectrons]
   Bool_t          ele_passHLTFilter[5][100];   //[nElectrons]
   Int_t           nPhotons;
   Int_t           nPhotons_overlap;
   Float_t         phoE[5];   //[nPhotons]
   Float_t         phoPt[5];   //[nPhotons]
   Float_t         phoEta[5];   //[nPhotons]
   Float_t         phoPhi[5];   //[nPhotons]
   Float_t         phoSigmaIetaIeta[5];   //[nPhotons]
   Float_t         phoFull5x5SigmaIetaIeta[5];   //[nPhotons]
   Float_t         phoR9[5];   //[nPhotons]
   Float_t         pho_sminor[5];   //[nPhotons]
   Float_t         pho_smajor[5];   //[nPhotons]
   Float_t         pho_HoverE[5];   //[nPhotons]
   Float_t         pho_sumChargedHadronPtAllVertices[5][1000];   //[nPhotons]
   Float_t         pho_sumChargedHadronPt[5];   //[nPhotons]
   Float_t         pho_sumNeutralHadronEt[5];   //[nPhotons]
   Float_t         pho_sumPhotonEt[5];   //[nPhotons]
   Float_t         pho_ecalPFClusterIso[5];   //[nPhotons]
   Float_t         pho_hcalPFClusterIso[5];   //[nPhotons]
   Float_t         pho_trkSumPtHollowConeDR03[5];   //[nPhotons]
   Float_t         pho_sumWorstVertexChargedHadronPt[5];   //[nPhotons]
   Float_t         pho_pfIsoChargedHadronIso[5];   //[nPhotons]
   Float_t         pho_pfIsoChargedHadronIsoWrongVtx[5];   //[nPhotons]
   Float_t         pho_pfIsoNeutralHadronIso[5];   //[nPhotons]
   Float_t         pho_pfIsoPhotonIso[5];   //[nPhotons]
   Float_t         pho_pfIsoModFrixione[5];   //[nPhotons]
   Float_t         pho_pfIsoSumPUPt[5];   //[nPhotons]
   Bool_t          pho_isConversion[5];   //[nPhotons]
   Bool_t          pho_passEleVeto[5];   //[nPhotons]
   Float_t         pho_RegressionE[5];   //[nPhotons]
   Float_t         pho_RegressionEUncertainty[5];   //[nPhotons]
   Float_t         pho_IDMVA[5];   //[nPhotons]
   Float_t         pho_superClusterEnergy[5];   //[nPhotons]
   Float_t         pho_superClusterRawEnergy[5];   //[nPhotons]
   Float_t         pho_superClusterEta[5];   //[nPhotons]
   Float_t         pho_superClusterPhi[5];   //[nPhotons]
   Float_t         pho_superClusterX[5];   //[nPhotons]
   Float_t         pho_superClusterY[5];   //[nPhotons]
   Float_t         pho_superClusterZ[5];   //[nPhotons]
   Bool_t          pho_hasPixelSeed[5];   //[nPhotons]
   Bool_t          pho_passHLTFilter[5][100];   //[nPhotons]
   Int_t           pho_convType[5];   //[nPhotons]
   Float_t         pho_convTrkZ[5];   //[nPhotons]
   Float_t         pho_convTrkClusZ[5];   //[nPhotons]
   Float_t         pho_vtxSumPx[5][1000];   //[nPhotons]
   Float_t         pho_vtxSumPy[5][1000];   //[nPhotons]
   Bool_t          pho_isStandardPhoton[5];   //[nPhotons]
   Float_t         pho_seedRecHitSwitchToGain6[5];   //[nPhotons]
   Float_t         pho_seedRecHitSwitchToGain1[5];   //[nPhotons]
   Float_t         pho_anyRecHitSwitchToGain6[5];   //[nPhotons]
   Float_t         pho_anyRecHitSwitchToGain1[5];   //[nPhotons]
   Int_t           nCscClusters;
   Float_t         cscClusterX[2];   //[nCscClusters]
   Float_t         cscClusterY[2];   //[nCscClusters]
   Float_t         cscClusterZ[2];   //[nCscClusters]
   Float_t         cscClusterTime[2];   //[nCscClusters]
   Float_t         cscClusterTimeSpread[2];   //[nCscClusters]
   Float_t         cscClusterGenMuonDeltaR[2];   //[nCscClusters]
   Float_t         cscClusterMajorAxis[2];   //[nCscClusters]
   Float_t         cscClusterMinorAxis[2];   //[nCscClusters]
   Float_t         cscClusterEtaPhiSpread[2];   //[nCscClusters]
   Float_t         cscClusterPhiSpread[2];   //[nCscClusters]
   Float_t         cscClusterEtaSpread[2];   //[nCscClusters]
   Float_t         cscClusterXSpread[2];   //[nCscClusters]
   Float_t         cscClusterYSpread[2];   //[nCscClusters]
   Float_t         cscClusterZSpread[2];   //[nCscClusters]
   Float_t         cscClusterPhi[2];   //[nCscClusters]
   Float_t         cscClusterEta[2];   //[nCscClusters]
   Float_t         cscClusterJetVetoPt[2];   //[nCscClusters]
   Float_t         cscClusterJetVetoE[2];   //[nCscClusters]
   Float_t         cscClusterMuonVetoPt[2];   //[nCscClusters]
   Float_t         cscClusterMuonVetoE[2];   //[nCscClusters]
   Float_t         cscClusterCaloJetVeto[2];   //[nCscClusters]
   Int_t           cscClusterSize[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberPlus11[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberPlus12[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberPlus13[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberPlus21[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberPlus22[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberPlus31[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberPlus32[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberPlus41[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberPlus42[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberMinus11[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberMinus12[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberMinus13[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberMinus21[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberMinus22[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberMinus31[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberMinus32[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberMinus41[2];   //[nCscClusters]
   Int_t           cscClusterNSegmentChamberMinus42[2];   //[nCscClusters]
   Float_t         cscClusterMe11Ratio[2];   //[nCscClusters]
   Float_t         cscClusterMe12Ratio[2];   //[nCscClusters]
   Int_t           cscClusterNStation[2];   //[nCscClusters]
   Int_t           cscClusterMaxStation[2];   //[nCscClusters]
   Float_t         cscClusterMaxStationRatio[2];   //[nCscClusters]
   Int_t           cscClusterNChamber[2];   //[nCscClusters]
   Int_t           cscClusterMaxChamber[2];   //[nCscClusters]
   Float_t         cscClusterMaxChamberRatio[2];   //[nCscClusters]
   Float_t         cscClusterVertexR[2];   //[nCscClusters]
   Float_t         cscClusterVertexZ[2];   //[nCscClusters]
   Float_t         cscClusterVertexDis[2];   //[nCscClusters]
   Float_t         cscClusterVertexChi2[2];   //[nCscClusters]
   Int_t           cscClusterVertexN1[2];   //[nCscClusters]
   Int_t           cscClusterVertexN5[2];   //[nCscClusters]
   Int_t           cscClusterVertexN10[2];   //[nCscClusters]
   Int_t           cscClusterVertexN15[2];   //[nCscClusters]
   Int_t           cscClusterVertexN20[2];   //[nCscClusters]
   Int_t           cscClusterVertexN[2];   //[nCscClusters]
   Int_t           nJets;
   Float_t         jetE[55];   //[nJets]
   Float_t         jetPt[55];   //[nJets]
   Float_t         jetEta[55];   //[nJets]
   Float_t         jetEt[55];   //[nJets]
   Float_t         jetPhi[55];   //[nJets]
   Float_t         jetCSV[55];   //[nJets]
   Float_t         jetCISV[55];   //[nJets]
   Float_t         jetProbb[55];   //[nJets]
   Float_t         jetProbc[55];   //[nJets]
   Float_t         jetProbudsg[55];   //[nJets]
   Float_t         jetProbbb[55];   //[nJets]
   Float_t         jetMass[55];   //[nJets]
   Float_t         jetJetArea[55];   //[nJets]
   Float_t         jetPileupE[55];   //[nJets]
   Float_t         jetPileupId[55];   //[nJets]
   Int_t           jetPileupIdFlag[55];   //[nJets]
   Bool_t          jetPassIDLoose[55];   //[nJets]
   Bool_t          jetPassIDTight[55];   //[nJets]
   Bool_t          jetPassMuFrac[55];   //[nJets]
   Bool_t          jetPassEleFrac[55];   //[nJets]
   Int_t           jetPartonFlavor[55];   //[nJets]
   Int_t           jetHadronFlavor[55];   //[nJets]
   Float_t         jetChargedEMEnergyFraction[55];   //[nJets]
   Float_t         jetNeutralEMEnergyFraction[55];   //[nJets]
   Float_t         jetChargedHadronEnergyFraction[55];   //[nJets]
   Float_t         jetNeutralHadronEnergyFraction[55];   //[nJets]
   Float_t         jetMuonEnergyFraction[55];   //[nJets]
   Float_t         jetHOEnergyFraction[55];   //[nJets]
   Float_t         jetHFHadronEnergyFraction[55];   //[nJets]
   Float_t         jetHFEMEnergyFraction[55];   //[nJets]
   Float_t         jetAllMuonPt[55];   //[nJets]
   Float_t         jetAllMuonEta[55];   //[nJets]
   Float_t         jetAllMuonPhi[55];   //[nJets]
   Float_t         jetAllMuonM[55];   //[nJets]
   Float_t         jetPtWeightedDZ[55];   //[nJets]
   Int_t           jetNRechits[55];   //[nJets]
   Float_t         jetRechitE[55];   //[nJets]
   Float_t         jetRechitT[55];   //[nJets]
   Float_t         jetRechitT_rms[55];   //[nJets]
   Float_t         jetRechitE_Error[55];   //[nJets]
   Float_t         jetRechitT_Error[55];   //[nJets]
   Float_t         jetAlphaMax[55];   //[nJets]
   Float_t         jetBetaMax[55];   //[nJets]
   Float_t         jetGammaMax_ET[55];   //[nJets]
   Float_t         jetGammaMax_EM[55];   //[nJets]
   Float_t         jetGammaMax_Hadronic[55];   //[nJets]
   Float_t         jetGammaMax[55];   //[nJets]
   Float_t         jetPtAllTracks[55];   //[nJets]
   Float_t         jetPtAllPVTracks[55];   //[nJets]
   Float_t         jetMedianTheta2D[55];   //[nJets]
   Float_t         jetMedianIP[55];   //[nJets]
   Float_t         jetMinDeltaRAllTracks[55];   //[nJets]
   Float_t         jetMinDeltaRPVTracks[55];   //[nJets]
   Float_t         jet_sig_et1[55];   //[nJets]
   Float_t         jet_sig_et2[55];   //[nJets]
   Float_t         jet_energy_frac[55];   //[nJets]
   Bool_t          jet_matched[55];   //[nJets]
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
   Bool_t          HLTDecision[602];
   Int_t           HLTPrescale[602];
   Int_t           nGenJets;
   Float_t         genJetE[73];   //[nGenJets]
   Float_t         genJetPt[73];   //[nGenJets]
   Float_t         genJetEta[73];   //[nGenJets]
   Float_t         genJetPhi[73];   //[nGenJets]
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
   vector<float>   *scaleWeights;
   vector<float>   *pdfWeights;
   vector<float>   *alphasWeights;
   Int_t           nGenParticle;
   Int_t           gParticleMotherId[39];   //[nGenParticle]
   Int_t           gParticleMotherIndex[39];   //[nGenParticle]
   Int_t           gParticleId[39];   //[nGenParticle]
   Int_t           gParticleStatus[39];   //[nGenParticle]
   Float_t         gParticleE[39];   //[nGenParticle]
   Float_t         gParticlePt[39];   //[nGenParticle]
   Float_t         gParticlePx[39];   //[nGenParticle]
   Float_t         gParticlePy[39];   //[nGenParticle]
   Float_t         gParticlePz[39];   //[nGenParticle]
   Float_t         gParticleEta[39];   //[nGenParticle]
   Float_t         gParticlePhi[39];   //[nGenParticle]
   Float_t         gParticleProdVertexX[39];   //[nGenParticle]
   Float_t         gParticleProdVertexY[39];   //[nGenParticle]
   Float_t         gParticleProdVertexZ[39];   //[nGenParticle]
   Float_t         gLLP_prod_vertex_x[2];
   Float_t         gLLP_prod_vertex_y[2];
   Float_t         gLLP_prod_vertex_z[2];
   Float_t         gLLP_decay_vertex_x[2];
   Float_t         gLLP_decay_vertex_y[2];
   Float_t         gLLP_decay_vertex_z[2];
   Float_t         gLLP_beta[2];
   Float_t         gLLP_e[2];
   Float_t         gLLP_pt[2];
   Float_t         gLLP_eta[2];
   Float_t         gLLP_phi[2];
   Float_t         gLLP_travel_time[2];
   Float_t         gLLP_daughter_travel_time[4];
   Int_t           gLLP_daughter_id[4];
   Float_t         gLLP_daughter_pt[4];
   Float_t         gLLP_daughter_eta[4];
   Float_t         gLLP_daughter_phi[4];
   Float_t         gLLP_daughter_eta_ecalcorr[4];
   Float_t         gLLP_daughter_phi_ecalcorr[4];
   Float_t         gLLP_daughter_e[4];
   Float_t         photon_travel_time[4];
   Float_t         gen_time[4];
   Float_t         gen_time_pv[4];
   Float_t         gLLP_min_delta_r_match_calojet[4];
   UInt_t          gLLP_daughter_match_calojet_index[4];
   UInt_t          gLLP_daughter_match_jet_index[4];
   Float_t         gLLP_min_delta_r_match_jet[4];

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
   TBranch        *b_ele_passSingleEleTagFilter;   //!
   TBranch        *b_ele_passTPOneTagFilter;   //!
   TBranch        *b_ele_passTPTwoTagFilter;   //!
   TBranch        *b_ele_passTPOneProbeFilter;   //!
   TBranch        *b_ele_passTPTwoProbeFilter;   //!
   TBranch        *b_ele_passHLTFilter;   //!
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
   TBranch        *b_nCscClusters;   //!
   TBranch        *b_cscClusterX;   //!
   TBranch        *b_cscClusterY;   //!
   TBranch        *b_cscClusterZ;   //!
   TBranch        *b_cscClusterTime;   //!
   TBranch        *b_cscClusterTimeSpread;   //!
   TBranch        *b_cscClusterGenMuonDeltaR;   //!
   TBranch        *b_cscClusterMajorAxis;   //!
   TBranch        *b_cscClusterMinorAxis;   //!
   TBranch        *b_cscClusterEtaPhiSpread;   //!
   TBranch        *b_cscClusterPhiSpread;   //!
   TBranch        *b_cscClusterEtaSpread;   //!
   TBranch        *b_cscClusterXSpread;   //!
   TBranch        *b_cscClusterYSpread;   //!
   TBranch        *b_cscClusterZSpread;   //!
   TBranch        *b_cscClusterPhi;   //!
   TBranch        *b_cscClusterEta;   //!
   TBranch        *b_cscClusterJetVetoPt;   //!
   TBranch        *b_cscClusterJetVetoE;   //!
   TBranch        *b_cscClusterMuonVetoPt;   //!
   TBranch        *b_cscClusterMuonVetoE;   //!
   TBranch        *b_cscClusterCaloJetVeto;   //!
   TBranch        *b_cscClusterSize;   //!
   TBranch        *b_cscClusterNSegmentChamberPlus11;   //!
   TBranch        *b_cscClusterNSegmentChamberPlus12;   //!
   TBranch        *b_cscClusterNSegmentChamberPlus13;   //!
   TBranch        *b_cscClusterNSegmentChamberPlus21;   //!
   TBranch        *b_cscClusterNSegmentChamberPlus22;   //!
   TBranch        *b_cscClusterNSegmentChamberPlus31;   //!
   TBranch        *b_cscClusterNSegmentChamberPlus32;   //!
   TBranch        *b_cscClusterNSegmentChamberPlus41;   //!
   TBranch        *b_cscClusterNSegmentChamberPlus42;   //!
   TBranch        *b_cscClusterNSegmentChamberMinus11;   //!
   TBranch        *b_cscClusterNSegmentChamberMinus12;   //!
   TBranch        *b_cscClusterNSegmentChamberMinus13;   //!
   TBranch        *b_cscClusterNSegmentChamberMinus21;   //!
   TBranch        *b_cscClusterNSegmentChamberMinus22;   //!
   TBranch        *b_cscClusterNSegmentChamberMinus31;   //!
   TBranch        *b_cscClusterNSegmentChamberMinus32;   //!
   TBranch        *b_cscClusterNSegmentChamberMinus41;   //!
   TBranch        *b_cscClusterNSegmentChamberMinus42;   //!
   TBranch        *b_cscClusterMe11Ratio;   //!
   TBranch        *b_cscClusterMe12Ratio;   //!
   TBranch        *b_cscClusterNStation;   //!
   TBranch        *b_cscClusterMaxStation;   //!
   TBranch        *b_cscClusterMaxStationRatio;   //!
   TBranch        *b_cscClusterNChamber;   //!
   TBranch        *b_cscClusterMaxChamber;   //!
   TBranch        *b_cscClusterMaxChamberRatio;   //!
   TBranch        *b_cscClusterVertexR;   //!
   TBranch        *b_cscClusterVertexZ;   //!
   TBranch        *b_cscClusterVertexDis;   //!
   TBranch        *b_cscClusterVertexChi2;   //!
   TBranch        *b_cscClusterVertexN1;   //!
   TBranch        *b_cscClusterVertexN5;   //!
   TBranch        *b_cscClusterVertexN10;   //!
   TBranch        *b_cscClusterVertexN15;   //!
   TBranch        *b_cscClusterVertexN20;   //!
   TBranch        *b_cscClusterVertexN;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetEt;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetCSV;   //!
   TBranch        *b_jetCISV;   //!
   TBranch        *b_jetProbb;   //!
   TBranch        *b_jetProbc;   //!
   TBranch        *b_jetProbudsg;   //!
   TBranch        *b_jetProbbb;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetJetArea;   //!
   TBranch        *b_jetPileupE;   //!
   TBranch        *b_jetPileupId;   //!
   TBranch        *b_jetPileupIdFlag;   //!
   TBranch        *b_jetPassIDLoose;   //!
   TBranch        *b_jetPassIDTight;   //!
   TBranch        *b_jetPassMuFrac;   //!
   TBranch        *b_jetPassEleFrac;   //!
   TBranch        *b_jetPartonFlavor;   //!
   TBranch        *b_jetHadronFlavor;   //!
   TBranch        *b_jetChargedEMEnergyFraction;   //!
   TBranch        *b_jetNeutralEMEnergyFraction;   //!
   TBranch        *b_jetChargedHadronEnergyFraction;   //!
   TBranch        *b_jetNeutralHadronEnergyFraction;   //!
   TBranch        *b_jetMuonEnergyFraction;   //!
   TBranch        *b_jetHOEnergyFraction;   //!
   TBranch        *b_jetHFHadronEnergyFraction;   //!
   TBranch        *b_jetHFEMEnergyFraction;   //!
   TBranch        *b_jetAllMuonPt;   //!
   TBranch        *b_jetAllMuonEta;   //!
   TBranch        *b_jetAllMuonPhi;   //!
   TBranch        *b_jetAllMuonM;   //!
   TBranch        *b_jetPtWeightedDZ;   //!
   TBranch        *b_jetNRechits;   //!
   TBranch        *b_jetRechitE;   //!
   TBranch        *b_jetRechitT;   //!
   TBranch        *b_jetRechitT_rms;   //!
   TBranch        *b_jetRechitE_Error;   //!
   TBranch        *b_jetRechitT_Error;   //!
   TBranch        *b_jetAlphaMax;   //!
   TBranch        *b_jetBetaMax;   //!
   TBranch        *b_jetGammaMax_ET;   //!
   TBranch        *b_jetGammaMax_EM;   //!
   TBranch        *b_jetGammaMax_Hadronic;   //!
   TBranch        *b_jetGammaMax;   //!
   TBranch        *b_jetPtAllTracks;   //!
   TBranch        *b_jetPtAllPVTracks;   //!
   TBranch        *b_jetMedianTheta2D;   //!
   TBranch        *b_jetMedianIP;   //!
   TBranch        *b_jetMinDeltaRAllTracks;   //!
   TBranch        *b_jetMinDeltaRPVTracks;   //!
   TBranch        *b_jet_sig_et1;   //!
   TBranch        *b_jet_sig_et2;   //!
   TBranch        *b_jet_energy_frac;   //!
   TBranch        *b_jet_matched;   //!
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
   TBranch        *b_gParticlePx;   //!
   TBranch        *b_gParticlePy;   //!
   TBranch        *b_gParticlePz;   //!
   TBranch        *b_gParticleEta;   //!
   TBranch        *b_gParticlePhi;   //!
   TBranch        *b_gParticleProdVertexX;   //!
   TBranch        *b_gParticleProdVertexY;   //!
   TBranch        *b_gParticleProdVertexZ;   //!
   TBranch        *b_gLLP_prod_vertex_x;   //!
   TBranch        *b_gLLP_prod_vertex_y;   //!
   TBranch        *b_gLLP_prod_vertex_z;   //!
   TBranch        *b_gLLP_decay_vertex_x;   //!
   TBranch        *b_gLLP_decay_vertex_y;   //!
   TBranch        *b_gLLP_decay_vertex_z;   //!
   TBranch        *b_gLLP_beta;   //!
   TBranch        *b_gLLP_e;   //!
   TBranch        *b_gLLP_pt;   //!
   TBranch        *b_gLLP_eta;   //!
   TBranch        *b_gLLP_phi;   //!
   TBranch        *b_gLLP_travel_time;   //!
   TBranch        *b_gLLP_daughter_travel_time;   //!
   TBranch        *b_gLLP_daughter_id;   //!
   TBranch        *b_gLLP_daughter_pt;   //!
   TBranch        *b_gLLP_daughter_eta;   //!
   TBranch        *b_gLLP_daughter_phi;   //!
   TBranch        *b_gLLP_daughter_eta_ecalcorr;   //!
   TBranch        *b_gLLP_daughter_phi_ecalcorr;   //!
   TBranch        *b_gLLP_daughter_e;   //!
   TBranch        *b_photon_travel_time;   //!
   TBranch        *b_gen_time;   //!
   TBranch        *b_gen_time_pv;   //!
   TBranch        *b_gLLP_min_delta_r_match_calojet;   //!
   TBranch        *b_gLLP_daughter_match_calojet_index;   //!
   TBranch        *b_gLLP_daughter_match_jet_index;   //!
   TBranch        *b_gLLP_min_delta_r_match_jet;   //!

   llpntupler_event(TTree *tree=0);
   virtual ~llpntupler_event();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef llpntupler_event_cxx
llpntupler_event::llpntupler_event(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/mnt/hadoop/store/group/phys_exotica/delayedjets/llpntuple/V1p8/MC_Summer16/v2/sixie/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/Run2_LLPNtupler_V1p8_MC_Summer16_RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v2_v2_v2/191011_040441/0000/llp_ntupler_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/mnt/hadoop/store/group/phys_exotica/delayedjets/llpntuple/V1p8/MC_Summer16/v2/sixie/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/Run2_LLPNtupler_V1p8_MC_Summer16_RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v2_v2_v2/191011_040441/0000/llp_ntupler_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/mnt/hadoop/store/group/phys_exotica/delayedjets/llpntuple/V1p8/MC_Summer16/v2/sixie/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/Run2_LLPNtupler_V1p8_MC_Summer16_RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v2_v2_v2/191011_040441/0000/llp_ntupler_1.root:/ntuples");
      dir->GetObject("llp",tree);

   }
   Init(tree);
}

llpntupler_event::~llpntupler_event()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t llpntupler_event::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t llpntupler_event::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void llpntupler_event::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   scaleWeights = 0;
   pdfWeights = 0;
   alphasWeights = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("runNum", &runNum, &b_runNum);
   fChain->SetBranchAddress("nSlimmedSecondV", &nSlimmedSecondV, &b_nSlimmedSecondV);
   fChain->SetBranchAddress("lumiNum", &lumiNum, &b_lumiNum);
   fChain->SetBranchAddress("eventNum", &eventNum, &b_eventNum);
   fChain->SetBranchAddress("eventTime", &eventTime, &b_eventTime);
   fChain->SetBranchAddress("pvX", &pvX, &b_pvX);
   fChain->SetBranchAddress("pvY", &pvY, &b_pvY);
   fChain->SetBranchAddress("pvZ", &pvZ, &b_pvZ);
   fChain->SetBranchAddress("fixedGridRhoAll", &fixedGridRhoAll, &b_fixedGridRhoAll);
   fChain->SetBranchAddress("nBunchXing", &nBunchXing, &b_nBunchXing);
   fChain->SetBranchAddress("BunchXing", BunchXing, &b_BunchXing);
   fChain->SetBranchAddress("nPU", nPU, &b_nPU);
   fChain->SetBranchAddress("nPUmean", nPUmean, &b_nPUmean);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("muonE", muonE, &b_muonE);
   fChain->SetBranchAddress("muonPt", muonPt, &b_muonPt);
   fChain->SetBranchAddress("muonEta", muonEta, &b_muonEta);
   fChain->SetBranchAddress("muonPhi", muonPhi, &b_muonPhi);
   fChain->SetBranchAddress("muonCharge", muonCharge, &b_muonCharge);
   fChain->SetBranchAddress("muonIsLoose", muonIsLoose, &b_muonIsLoose);
   fChain->SetBranchAddress("muonIsMedium", muonIsMedium, &b_muonIsMedium);
   fChain->SetBranchAddress("muonIsTight", muonIsTight, &b_muonIsTight);
   fChain->SetBranchAddress("muon_d0", muon_d0, &b_muon_d0);
   fChain->SetBranchAddress("muon_dZ", muon_dZ, &b_muon_dZ);
   fChain->SetBranchAddress("muon_ip3d", muon_ip3d, &b_muon_ip3d);
   fChain->SetBranchAddress("muon_ip3dSignificance", muon_ip3dSignificance, &b_muon_ip3dSignificance);
   fChain->SetBranchAddress("muonType", muonType, &b_muonType);
   fChain->SetBranchAddress("muonQuality", muonQuality, &b_muonQuality);
   fChain->SetBranchAddress("muon_pileupIso", muon_pileupIso, &b_muon_pileupIso);
   fChain->SetBranchAddress("muon_chargedIso", muon_chargedIso, &b_muon_chargedIso);
   fChain->SetBranchAddress("muon_photonIso", muon_photonIso, &b_muon_photonIso);
   fChain->SetBranchAddress("muon_neutralHadIso", muon_neutralHadIso, &b_muon_neutralHadIso);
   fChain->SetBranchAddress("muon_ptrel", muon_ptrel, &b_muon_ptrel);
   fChain->SetBranchAddress("muon_passSingleMuTagFilter", muon_passSingleMuTagFilter, &b_muon_passSingleMuTagFilter);
   fChain->SetBranchAddress("muon_passHLTFilter", muon_passHLTFilter, &b_muon_passHLTFilter);
   fChain->SetBranchAddress("muon_validFractionTrackerHits", muon_validFractionTrackerHits, &b_muon_validFractionTrackerHits);
   fChain->SetBranchAddress("muon_isGlobal", muon_isGlobal, &b_muon_isGlobal);
   fChain->SetBranchAddress("muon_normChi2", muon_normChi2, &b_muon_normChi2);
   fChain->SetBranchAddress("muon_chi2LocalPosition", muon_chi2LocalPosition, &b_muon_chi2LocalPosition);
   fChain->SetBranchAddress("muon_kinkFinder", muon_kinkFinder, &b_muon_kinkFinder);
   fChain->SetBranchAddress("muon_segmentCompatability", muon_segmentCompatability, &b_muon_segmentCompatability);
   fChain->SetBranchAddress("muonIsICHEPMedium", muonIsICHEPMedium, &b_muonIsICHEPMedium);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("eleE", eleE, &b_eleE);
   fChain->SetBranchAddress("elePt", elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleCharge", eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleEta_SC", eleEta_SC, &b_eleEta_SC);
   fChain->SetBranchAddress("eleSigmaIetaIeta", eleSigmaIetaIeta, &b_eleSigmaIetaIeta);
   fChain->SetBranchAddress("eleFull5x5SigmaIetaIeta", eleFull5x5SigmaIetaIeta, &b_eleFull5x5SigmaIetaIeta);
   fChain->SetBranchAddress("eleR9", eleR9, &b_eleR9);
   fChain->SetBranchAddress("ele_dEta", ele_dEta, &b_ele_dEta);
   fChain->SetBranchAddress("ele_dPhi", ele_dPhi, &b_ele_dPhi);
   fChain->SetBranchAddress("ele_HoverE", ele_HoverE, &b_ele_HoverE);
   fChain->SetBranchAddress("ele_d0", ele_d0, &b_ele_d0);
   fChain->SetBranchAddress("ele_dZ", ele_dZ, &b_ele_dZ);
   fChain->SetBranchAddress("ele_ip3d", ele_ip3d, &b_ele_ip3d);
   fChain->SetBranchAddress("ele_ip3dSignificance", ele_ip3dSignificance, &b_ele_ip3dSignificance);
   fChain->SetBranchAddress("ele_pileupIso", ele_pileupIso, &b_ele_pileupIso);
   fChain->SetBranchAddress("ele_chargedIso", ele_chargedIso, &b_ele_chargedIso);
   fChain->SetBranchAddress("ele_photonIso", ele_photonIso, &b_ele_photonIso);
   fChain->SetBranchAddress("ele_neutralHadIso", ele_neutralHadIso, &b_ele_neutralHadIso);
   fChain->SetBranchAddress("ele_MissHits", ele_MissHits, &b_ele_MissHits);
   fChain->SetBranchAddress("ele_PassConvVeto", ele_PassConvVeto, &b_ele_PassConvVeto);
   fChain->SetBranchAddress("ele_OneOverEminusOneOverP", ele_OneOverEminusOneOverP, &b_ele_OneOverEminusOneOverP);
   fChain->SetBranchAddress("ele_IDMVAGeneralPurpose", ele_IDMVAGeneralPurpose, &b_ele_IDMVAGeneralPurpose);
   fChain->SetBranchAddress("ele_IDMVACategoryGeneralPurpose", ele_IDMVACategoryGeneralPurpose, &b_ele_IDMVACategoryGeneralPurpose);
   fChain->SetBranchAddress("ele_IDMVAHZZ", ele_IDMVAHZZ, &b_ele_IDMVAHZZ);
   fChain->SetBranchAddress("ele_IDMVACategoryHZZ", ele_IDMVACategoryHZZ, &b_ele_IDMVACategoryHZZ);
   fChain->SetBranchAddress("ele_passSingleEleTagFilter", ele_passSingleEleTagFilter, &b_ele_passSingleEleTagFilter);
   fChain->SetBranchAddress("ele_passTPOneTagFilter", ele_passTPOneTagFilter, &b_ele_passTPOneTagFilter);
   fChain->SetBranchAddress("ele_passTPTwoTagFilter", ele_passTPTwoTagFilter, &b_ele_passTPTwoTagFilter);
   fChain->SetBranchAddress("ele_passTPOneProbeFilter", ele_passTPOneProbeFilter, &b_ele_passTPOneProbeFilter);
   fChain->SetBranchAddress("ele_passTPTwoProbeFilter", ele_passTPTwoProbeFilter, &b_ele_passTPTwoProbeFilter);
   fChain->SetBranchAddress("ele_passHLTFilter", ele_passHLTFilter, &b_ele_passHLTFilter);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("nPhotons_overlap", &nPhotons_overlap, &b_nPhotons_overlap);
   fChain->SetBranchAddress("phoE", phoE, &b_phoE);
   fChain->SetBranchAddress("phoPt", phoPt, &b_phoPt);
   fChain->SetBranchAddress("phoEta", phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoSigmaIetaIeta", phoSigmaIetaIeta, &b_phoSigmaIetaIeta);
   fChain->SetBranchAddress("phoFull5x5SigmaIetaIeta", phoFull5x5SigmaIetaIeta, &b_phoFull5x5SigmaIetaIeta);
   fChain->SetBranchAddress("phoR9", phoR9, &b_phoR9);
   fChain->SetBranchAddress("pho_sminor", pho_sminor, &b_pho_sminor);
   fChain->SetBranchAddress("pho_smajor", pho_smajor, &b_pho_smajor);
   fChain->SetBranchAddress("pho_HoverE", pho_HoverE, &b_pho_HoverE);
   fChain->SetBranchAddress("pho_sumChargedHadronPtAllVertices", pho_sumChargedHadronPtAllVertices, &b_pho_sumChargedHadronPtAllVertices);
   fChain->SetBranchAddress("pho_sumChargedHadronPt", pho_sumChargedHadronPt, &b_pho_sumChargedHadronPt);
   fChain->SetBranchAddress("pho_sumNeutralHadronEt", pho_sumNeutralHadronEt, &b_pho_sumNeutralHadronEt);
   fChain->SetBranchAddress("pho_sumPhotonEt", pho_sumPhotonEt, &b_pho_sumPhotonEt);
   fChain->SetBranchAddress("pho_ecalPFClusterIso", pho_ecalPFClusterIso, &b_pho_ecalPFClusterIso);
   fChain->SetBranchAddress("pho_hcalPFClusterIso", pho_hcalPFClusterIso, &b_pho_hcalPFClusterIso);
   fChain->SetBranchAddress("pho_trkSumPtHollowConeDR03", pho_trkSumPtHollowConeDR03, &b_pho_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("pho_sumWorstVertexChargedHadronPt", pho_sumWorstVertexChargedHadronPt, &b_pho_sumWorstVertexChargedHadronPt);
   fChain->SetBranchAddress("pho_pfIsoChargedHadronIso", pho_pfIsoChargedHadronIso, &b_pho_pfIsoChargedHadronIso);
   fChain->SetBranchAddress("pho_pfIsoChargedHadronIsoWrongVtx", pho_pfIsoChargedHadronIsoWrongVtx, &b_pho_pfIsoChargedHadronIsoWrongVtx);
   fChain->SetBranchAddress("pho_pfIsoNeutralHadronIso", pho_pfIsoNeutralHadronIso, &b_pho_pfIsoNeutralHadronIso);
   fChain->SetBranchAddress("pho_pfIsoPhotonIso", pho_pfIsoPhotonIso, &b_pho_pfIsoPhotonIso);
   fChain->SetBranchAddress("pho_pfIsoModFrixione", pho_pfIsoModFrixione, &b_pho_pfIsoModFrixione);
   fChain->SetBranchAddress("pho_pfIsoSumPUPt", pho_pfIsoSumPUPt, &b_pho_pfIsoSumPUPt);
   fChain->SetBranchAddress("pho_isConversion", pho_isConversion, &b_pho_isConversion);
   fChain->SetBranchAddress("pho_passEleVeto", pho_passEleVeto, &b_pho_passEleVeto);
   fChain->SetBranchAddress("pho_RegressionE", pho_RegressionE, &b_pho_RegressionE);
   fChain->SetBranchAddress("pho_RegressionEUncertainty", pho_RegressionEUncertainty, &b_pho_RegressionEUncertainty);
   fChain->SetBranchAddress("pho_IDMVA", pho_IDMVA, &b_pho_IDMVA);
   fChain->SetBranchAddress("pho_superClusterEnergy", pho_superClusterEnergy, &b_pho_superClusterEnergy);
   fChain->SetBranchAddress("pho_superClusterRawEnergy", pho_superClusterRawEnergy, &b_pho_superClusterRawEnergy);
   fChain->SetBranchAddress("pho_superClusterEta", pho_superClusterEta, &b_pho_superClusterEta);
   fChain->SetBranchAddress("pho_superClusterPhi", pho_superClusterPhi, &b_pho_superClusterPhi);
   fChain->SetBranchAddress("pho_superClusterX", pho_superClusterX, &b_pho_superClusterX);
   fChain->SetBranchAddress("pho_superClusterY", pho_superClusterY, &b_pho_superClusterY);
   fChain->SetBranchAddress("pho_superClusterZ", pho_superClusterZ, &b_pho_superClusterZ);
   fChain->SetBranchAddress("pho_hasPixelSeed", pho_hasPixelSeed, &b_pho_hasPixelSeed);
   fChain->SetBranchAddress("pho_passHLTFilter", pho_passHLTFilter, &b_pho_passHLTFilter);
   fChain->SetBranchAddress("pho_convType", pho_convType, &b_pho_convType);
   fChain->SetBranchAddress("pho_convTrkZ", pho_convTrkZ, &b_pho_convTrkZ);
   fChain->SetBranchAddress("pho_convTrkClusZ", pho_convTrkClusZ, &b_pho_convTrkClusZ);
   fChain->SetBranchAddress("pho_vtxSumPx", pho_vtxSumPx, &b_pho_vtxSumPx);
   fChain->SetBranchAddress("pho_vtxSumPy", pho_vtxSumPy, &b_pho_vtxSumPy);
   fChain->SetBranchAddress("pho_isStandardPhoton", pho_isStandardPhoton, &b_pho_isStandardPhoton);
   fChain->SetBranchAddress("pho_seedRecHitSwitchToGain6", pho_seedRecHitSwitchToGain6, &b_pho_seedRecHitSwitchToGain6);
   fChain->SetBranchAddress("pho_seedRecHitSwitchToGain1", pho_seedRecHitSwitchToGain1, &b_pho_seedRecHitSwitchToGain1);
   fChain->SetBranchAddress("pho_anyRecHitSwitchToGain6", pho_anyRecHitSwitchToGain6, &b_pho_anyRecHitSwitchToGain6);
   fChain->SetBranchAddress("pho_anyRecHitSwitchToGain1", pho_anyRecHitSwitchToGain1, &b_pho_anyRecHitSwitchToGain1);
   fChain->SetBranchAddress("nCscClusters", &nCscClusters, &b_nCscClusters);
   fChain->SetBranchAddress("cscClusterX", cscClusterX, &b_cscClusterX);
   fChain->SetBranchAddress("cscClusterY", cscClusterY, &b_cscClusterY);
   fChain->SetBranchAddress("cscClusterZ", cscClusterZ, &b_cscClusterZ);
   fChain->SetBranchAddress("cscClusterTime", cscClusterTime, &b_cscClusterTime);
   fChain->SetBranchAddress("cscClusterTimeSpread", cscClusterTimeSpread, &b_cscClusterTimeSpread);
   fChain->SetBranchAddress("cscClusterGenMuonDeltaR", cscClusterGenMuonDeltaR, &b_cscClusterGenMuonDeltaR);
   fChain->SetBranchAddress("cscClusterMajorAxis", cscClusterMajorAxis, &b_cscClusterMajorAxis);
   fChain->SetBranchAddress("cscClusterMinorAxis", cscClusterMinorAxis, &b_cscClusterMinorAxis);
   fChain->SetBranchAddress("cscClusterEtaPhiSpread", cscClusterEtaPhiSpread, &b_cscClusterEtaPhiSpread);
   fChain->SetBranchAddress("cscClusterPhiSpread", cscClusterPhiSpread, &b_cscClusterPhiSpread);
   fChain->SetBranchAddress("cscClusterEtaSpread", cscClusterEtaSpread, &b_cscClusterEtaSpread);
   fChain->SetBranchAddress("cscClusterXSpread", cscClusterXSpread, &b_cscClusterXSpread);
   fChain->SetBranchAddress("cscClusterYSpread", cscClusterYSpread, &b_cscClusterYSpread);
   fChain->SetBranchAddress("cscClusterZSpread", cscClusterZSpread, &b_cscClusterZSpread);
   fChain->SetBranchAddress("cscClusterPhi", cscClusterPhi, &b_cscClusterPhi);
   fChain->SetBranchAddress("cscClusterEta", cscClusterEta, &b_cscClusterEta);
   fChain->SetBranchAddress("cscClusterJetVetoPt", cscClusterJetVetoPt, &b_cscClusterJetVetoPt);
   fChain->SetBranchAddress("cscClusterJetVetoE", cscClusterJetVetoE, &b_cscClusterJetVetoE);
   fChain->SetBranchAddress("cscClusterMuonVetoPt", cscClusterMuonVetoPt, &b_cscClusterMuonVetoPt);
   fChain->SetBranchAddress("cscClusterMuonVetoE", cscClusterMuonVetoE, &b_cscClusterMuonVetoE);
   fChain->SetBranchAddress("cscClusterCaloJetVeto", cscClusterCaloJetVeto, &b_cscClusterCaloJetVeto);
   fChain->SetBranchAddress("cscClusterSize", cscClusterSize, &b_cscClusterSize);
   fChain->SetBranchAddress("cscClusterNSegmentChamberPlus11", cscClusterNSegmentChamberPlus11, &b_cscClusterNSegmentChamberPlus11);
   fChain->SetBranchAddress("cscClusterNSegmentChamberPlus12", cscClusterNSegmentChamberPlus12, &b_cscClusterNSegmentChamberPlus12);
   fChain->SetBranchAddress("cscClusterNSegmentChamberPlus13", cscClusterNSegmentChamberPlus13, &b_cscClusterNSegmentChamberPlus13);
   fChain->SetBranchAddress("cscClusterNSegmentChamberPlus21", cscClusterNSegmentChamberPlus21, &b_cscClusterNSegmentChamberPlus21);
   fChain->SetBranchAddress("cscClusterNSegmentChamberPlus22", cscClusterNSegmentChamberPlus22, &b_cscClusterNSegmentChamberPlus22);
   fChain->SetBranchAddress("cscClusterNSegmentChamberPlus31", cscClusterNSegmentChamberPlus31, &b_cscClusterNSegmentChamberPlus31);
   fChain->SetBranchAddress("cscClusterNSegmentChamberPlus32", cscClusterNSegmentChamberPlus32, &b_cscClusterNSegmentChamberPlus32);
   fChain->SetBranchAddress("cscClusterNSegmentChamberPlus41", cscClusterNSegmentChamberPlus41, &b_cscClusterNSegmentChamberPlus41);
   fChain->SetBranchAddress("cscClusterNSegmentChamberPlus42", cscClusterNSegmentChamberPlus42, &b_cscClusterNSegmentChamberPlus42);
   fChain->SetBranchAddress("cscClusterNSegmentChamberMinus11", cscClusterNSegmentChamberMinus11, &b_cscClusterNSegmentChamberMinus11);
   fChain->SetBranchAddress("cscClusterNSegmentChamberMinus12", cscClusterNSegmentChamberMinus12, &b_cscClusterNSegmentChamberMinus12);
   fChain->SetBranchAddress("cscClusterNSegmentChamberMinus13", cscClusterNSegmentChamberMinus13, &b_cscClusterNSegmentChamberMinus13);
   fChain->SetBranchAddress("cscClusterNSegmentChamberMinus21", cscClusterNSegmentChamberMinus21, &b_cscClusterNSegmentChamberMinus21);
   fChain->SetBranchAddress("cscClusterNSegmentChamberMinus22", cscClusterNSegmentChamberMinus22, &b_cscClusterNSegmentChamberMinus22);
   fChain->SetBranchAddress("cscClusterNSegmentChamberMinus31", cscClusterNSegmentChamberMinus31, &b_cscClusterNSegmentChamberMinus31);
   fChain->SetBranchAddress("cscClusterNSegmentChamberMinus32", cscClusterNSegmentChamberMinus32, &b_cscClusterNSegmentChamberMinus32);
   fChain->SetBranchAddress("cscClusterNSegmentChamberMinus41", cscClusterNSegmentChamberMinus41, &b_cscClusterNSegmentChamberMinus41);
   fChain->SetBranchAddress("cscClusterNSegmentChamberMinus42", cscClusterNSegmentChamberMinus42, &b_cscClusterNSegmentChamberMinus42);
   fChain->SetBranchAddress("cscClusterMe11Ratio", cscClusterMe11Ratio, &b_cscClusterMe11Ratio);
   fChain->SetBranchAddress("cscClusterMe12Ratio", cscClusterMe12Ratio, &b_cscClusterMe12Ratio);
   fChain->SetBranchAddress("cscClusterNStation", cscClusterNStation, &b_cscClusterNStation);
   fChain->SetBranchAddress("cscClusterMaxStation", cscClusterMaxStation, &b_cscClusterMaxStation);
   fChain->SetBranchAddress("cscClusterMaxStationRatio", cscClusterMaxStationRatio, &b_cscClusterMaxStationRatio);
   fChain->SetBranchAddress("cscClusterNChamber", cscClusterNChamber, &b_cscClusterNChamber);
   fChain->SetBranchAddress("cscClusterMaxChamber", cscClusterMaxChamber, &b_cscClusterMaxChamber);
   fChain->SetBranchAddress("cscClusterMaxChamberRatio", cscClusterMaxChamberRatio, &b_cscClusterMaxChamberRatio);
   fChain->SetBranchAddress("cscClusterVertexR", cscClusterVertexR, &b_cscClusterVertexR);
   fChain->SetBranchAddress("cscClusterVertexZ", cscClusterVertexZ, &b_cscClusterVertexZ);
   fChain->SetBranchAddress("cscClusterVertexDis", cscClusterVertexDis, &b_cscClusterVertexDis);
   fChain->SetBranchAddress("cscClusterVertexChi2", cscClusterVertexChi2, &b_cscClusterVertexChi2);
   fChain->SetBranchAddress("cscClusterVertexN1", cscClusterVertexN1, &b_cscClusterVertexN1);
   fChain->SetBranchAddress("cscClusterVertexN5", cscClusterVertexN5, &b_cscClusterVertexN5);
   fChain->SetBranchAddress("cscClusterVertexN10", cscClusterVertexN10, &b_cscClusterVertexN10);
   fChain->SetBranchAddress("cscClusterVertexN15", cscClusterVertexN15, &b_cscClusterVertexN15);
   fChain->SetBranchAddress("cscClusterVertexN20", cscClusterVertexN20, &b_cscClusterVertexN20);
   fChain->SetBranchAddress("cscClusterVertexN", cscClusterVertexN, &b_cscClusterVertexN);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("jetE", jetE, &b_jetE);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetEt", jetEt, &b_jetEt);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetCSV", jetCSV, &b_jetCSV);
   fChain->SetBranchAddress("jetCISV", jetCISV, &b_jetCISV);
   fChain->SetBranchAddress("jetProbb", jetProbb, &b_jetProbb);
   fChain->SetBranchAddress("jetProbc", jetProbc, &b_jetProbc);
   fChain->SetBranchAddress("jetProbudsg", jetProbudsg, &b_jetProbudsg);
   fChain->SetBranchAddress("jetProbbb", jetProbbb, &b_jetProbbb);
   fChain->SetBranchAddress("jetMass", jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetJetArea", jetJetArea, &b_jetJetArea);
   fChain->SetBranchAddress("jetPileupE", jetPileupE, &b_jetPileupE);
   fChain->SetBranchAddress("jetPileupId", jetPileupId, &b_jetPileupId);
   fChain->SetBranchAddress("jetPileupIdFlag", jetPileupIdFlag, &b_jetPileupIdFlag);
   fChain->SetBranchAddress("jetPassIDLoose", jetPassIDLoose, &b_jetPassIDLoose);
   fChain->SetBranchAddress("jetPassIDTight", jetPassIDTight, &b_jetPassIDTight);
   fChain->SetBranchAddress("jetPassMuFrac", jetPassMuFrac, &b_jetPassMuFrac);
   fChain->SetBranchAddress("jetPassEleFrac", jetPassEleFrac, &b_jetPassEleFrac);
   fChain->SetBranchAddress("jetPartonFlavor", jetPartonFlavor, &b_jetPartonFlavor);
   fChain->SetBranchAddress("jetHadronFlavor", jetHadronFlavor, &b_jetHadronFlavor);
   fChain->SetBranchAddress("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction, &b_jetChargedEMEnergyFraction);
   fChain->SetBranchAddress("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction, &b_jetNeutralEMEnergyFraction);
   fChain->SetBranchAddress("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction, &b_jetChargedHadronEnergyFraction);
   fChain->SetBranchAddress("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction, &b_jetNeutralHadronEnergyFraction);
   fChain->SetBranchAddress("jetMuonEnergyFraction", jetMuonEnergyFraction, &b_jetMuonEnergyFraction);
   fChain->SetBranchAddress("jetHOEnergyFraction", jetHOEnergyFraction, &b_jetHOEnergyFraction);
   fChain->SetBranchAddress("jetHFHadronEnergyFraction", jetHFHadronEnergyFraction, &b_jetHFHadronEnergyFraction);
   fChain->SetBranchAddress("jetHFEMEnergyFraction", jetHFEMEnergyFraction, &b_jetHFEMEnergyFraction);
   fChain->SetBranchAddress("jetAllMuonPt", jetAllMuonPt, &b_jetAllMuonPt);
   fChain->SetBranchAddress("jetAllMuonEta", jetAllMuonEta, &b_jetAllMuonEta);
   fChain->SetBranchAddress("jetAllMuonPhi", jetAllMuonPhi, &b_jetAllMuonPhi);
   fChain->SetBranchAddress("jetAllMuonM", jetAllMuonM, &b_jetAllMuonM);
   fChain->SetBranchAddress("jetPtWeightedDZ", jetPtWeightedDZ, &b_jetPtWeightedDZ);
   fChain->SetBranchAddress("jetNRechits", jetNRechits, &b_jetNRechits);
   fChain->SetBranchAddress("jetRechitE", jetRechitE, &b_jetRechitE);
   fChain->SetBranchAddress("jetRechitT", jetRechitT, &b_jetRechitT);
   fChain->SetBranchAddress("jetRechitT_rms", jetRechitT_rms, &b_jetRechitT_rms);
   fChain->SetBranchAddress("jetRechitE_Error", jetRechitE_Error, &b_jetRechitE_Error);
   fChain->SetBranchAddress("jetRechitT_Error", jetRechitT_Error, &b_jetRechitT_Error);
   fChain->SetBranchAddress("jetAlphaMax", jetAlphaMax, &b_jetAlphaMax);
   fChain->SetBranchAddress("jetBetaMax", jetBetaMax, &b_jetBetaMax);
   fChain->SetBranchAddress("jetGammaMax_ET", jetGammaMax_ET, &b_jetGammaMax_ET);
   fChain->SetBranchAddress("jetGammaMax_EM", jetGammaMax_EM, &b_jetGammaMax_EM);
   fChain->SetBranchAddress("jetGammaMax_Hadronic", jetGammaMax_Hadronic, &b_jetGammaMax_Hadronic);
   fChain->SetBranchAddress("jetGammaMax", jetGammaMax, &b_jetGammaMax);
   fChain->SetBranchAddress("jetPtAllTracks", jetPtAllTracks, &b_jetPtAllTracks);
   fChain->SetBranchAddress("jetPtAllPVTracks", jetPtAllPVTracks, &b_jetPtAllPVTracks);
   fChain->SetBranchAddress("jetMedianTheta2D", jetMedianTheta2D, &b_jetMedianTheta2D);
   fChain->SetBranchAddress("jetMedianIP", jetMedianIP, &b_jetMedianIP);
   fChain->SetBranchAddress("jetMinDeltaRAllTracks", jetMinDeltaRAllTracks, &b_jetMinDeltaRAllTracks);
   fChain->SetBranchAddress("jetMinDeltaRPVTracks", jetMinDeltaRPVTracks, &b_jetMinDeltaRPVTracks);
   fChain->SetBranchAddress("jet_sig_et1", jet_sig_et1, &b_jet_sig_et1);
   fChain->SetBranchAddress("jet_sig_et2", jet_sig_et2, &b_jet_sig_et2);
   fChain->SetBranchAddress("jet_energy_frac", jet_energy_frac, &b_jet_energy_frac);
   fChain->SetBranchAddress("jet_matched", jet_matched, &b_jet_matched);
   fChain->SetBranchAddress("metPt", &metPt, &b_metPt);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("sumMET", &sumMET, &b_sumMET);
   fChain->SetBranchAddress("metType0Pt", &metType0Pt, &b_metType0Pt);
   fChain->SetBranchAddress("metType0Phi", &metType0Phi, &b_metType0Phi);
   fChain->SetBranchAddress("metType1Pt_raw", &metType1Pt_raw, &b_metType1Pt_raw);
   fChain->SetBranchAddress("metType1Pt", &metType1Pt, &b_metType1Pt);
   fChain->SetBranchAddress("metType1Px", &metType1Px, &b_metType1Px);
   fChain->SetBranchAddress("metType1Py", &metType1Py, &b_metType1Py);
   fChain->SetBranchAddress("metType1Eta", &metType1Eta, &b_metType1Eta);
   fChain->SetBranchAddress("metType1Phi", &metType1Phi, &b_metType1Phi);
   fChain->SetBranchAddress("metType1Phi_raw", &metType1Phi_raw, &b_metType1Phi_raw);
   fChain->SetBranchAddress("metType0Plus1Pt", &metType0Plus1Pt, &b_metType0Plus1Pt);
   fChain->SetBranchAddress("metType0Plus1Phi", &metType0Plus1Phi, &b_metType0Plus1Phi);
   fChain->SetBranchAddress("metNoHFPt", &metNoHFPt, &b_metNoHFPt);
   fChain->SetBranchAddress("metNoHFPhi", &metNoHFPhi, &b_metNoHFPhi);
   fChain->SetBranchAddress("metPuppiPt", &metPuppiPt, &b_metPuppiPt);
   fChain->SetBranchAddress("metPuppiPhi", &metPuppiPhi, &b_metPuppiPhi);
   fChain->SetBranchAddress("metCaloPt", &metCaloPt, &b_metCaloPt);
   fChain->SetBranchAddress("metCaloPhi", &metCaloPhi, &b_metCaloPhi);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHETightNoiseFilter", &Flag_HBHETightNoiseFilter, &b_Flag_HBHETightNoiseFilter);
   fChain->SetBranchAddress("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, &b_Flag_HBHEIsoNoiseFilter);
   fChain->SetBranchAddress("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter, &b_Flag_badChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_badMuonFilter", &Flag_badMuonFilter, &b_Flag_badMuonFilter);
   fChain->SetBranchAddress("Flag_badGlobalMuonFilter", &Flag_badGlobalMuonFilter, &b_Flag_badGlobalMuonFilter);
   fChain->SetBranchAddress("Flag_duplicateMuonFilter", &Flag_duplicateMuonFilter, &b_Flag_duplicateMuonFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, &b_Flag_trackingFailureFilter);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   fChain->SetBranchAddress("HLTDecision", HLTDecision, &b_HLTDecision);
   fChain->SetBranchAddress("HLTPrescale", HLTPrescale, &b_HLTPrescale);
   fChain->SetBranchAddress("nGenJets", &nGenJets, &b_nGenJets);
   fChain->SetBranchAddress("genJetE", genJetE, &b_genJetE);
   fChain->SetBranchAddress("genJetPt", genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetEta", genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPhi", genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genMetPtCalo", &genMetPtCalo, &b_genMetPtCalo);
   fChain->SetBranchAddress("genMetPhiCalo", &genMetPhiCalo, &b_genMetPhiCalo);
   fChain->SetBranchAddress("genMetPtTrue", &genMetPtTrue, &b_genMetPtTrue);
   fChain->SetBranchAddress("genMetPhiTrue", &genMetPhiTrue, &b_genMetPhiTrue);
   fChain->SetBranchAddress("genVertexX", &genVertexX, &b_genVertexX);
   fChain->SetBranchAddress("genVertexY", &genVertexY, &b_genVertexY);
   fChain->SetBranchAddress("genVertexZ", &genVertexZ, &b_genVertexZ);
   fChain->SetBranchAddress("genVertexT", &genVertexT, &b_genVertexT);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genSignalProcessID", &genSignalProcessID, &b_genSignalProcessID);
   fChain->SetBranchAddress("genQScale", &genQScale, &b_genQScale);
   fChain->SetBranchAddress("genAlphaQCD", &genAlphaQCD, &b_genAlphaQCD);
   fChain->SetBranchAddress("genAlphaQED", &genAlphaQED, &b_genAlphaQED);
   fChain->SetBranchAddress("scaleWeights", &scaleWeights, &b_scaleWeights);
   fChain->SetBranchAddress("pdfWeights", &pdfWeights, &b_pdfWeights);
   fChain->SetBranchAddress("alphasWeights", &alphasWeights, &b_alphasWeights);
   fChain->SetBranchAddress("nGenParticle", &nGenParticle, &b_nGenParticle);
   fChain->SetBranchAddress("gParticleMotherId", gParticleMotherId, &b_gParticleMotherId);
   fChain->SetBranchAddress("gParticleMotherIndex", gParticleMotherIndex, &b_gParticleMotherIndex);
   fChain->SetBranchAddress("gParticleId", gParticleId, &b_gParticleId);
   fChain->SetBranchAddress("gParticleStatus", gParticleStatus, &b_gParticleStatus);
   fChain->SetBranchAddress("gParticleE", gParticleE, &b_gParticleE);
   fChain->SetBranchAddress("gParticlePt", gParticlePt, &b_gParticlePt);
   fChain->SetBranchAddress("gParticlePx", gParticlePx, &b_gParticlePx);
   fChain->SetBranchAddress("gParticlePy", gParticlePy, &b_gParticlePy);
   fChain->SetBranchAddress("gParticlePz", gParticlePz, &b_gParticlePz);
   fChain->SetBranchAddress("gParticleEta", gParticleEta, &b_gParticleEta);
   fChain->SetBranchAddress("gParticlePhi", gParticlePhi, &b_gParticlePhi);
   fChain->SetBranchAddress("gParticleProdVertexX", gParticleProdVertexX, &b_gParticleProdVertexX);
   fChain->SetBranchAddress("gParticleProdVertexY", gParticleProdVertexY, &b_gParticleProdVertexY);
   fChain->SetBranchAddress("gParticleProdVertexZ", gParticleProdVertexZ, &b_gParticleProdVertexZ);
   fChain->SetBranchAddress("gLLP_prod_vertex_x", gLLP_prod_vertex_x, &b_gLLP_prod_vertex_x);
   fChain->SetBranchAddress("gLLP_prod_vertex_y", gLLP_prod_vertex_y, &b_gLLP_prod_vertex_y);
   fChain->SetBranchAddress("gLLP_prod_vertex_z", gLLP_prod_vertex_z, &b_gLLP_prod_vertex_z);
   fChain->SetBranchAddress("gLLP_decay_vertex_x", gLLP_decay_vertex_x, &b_gLLP_decay_vertex_x);
   fChain->SetBranchAddress("gLLP_decay_vertex_y", gLLP_decay_vertex_y, &b_gLLP_decay_vertex_y);
   fChain->SetBranchAddress("gLLP_decay_vertex_z", gLLP_decay_vertex_z, &b_gLLP_decay_vertex_z);
   fChain->SetBranchAddress("gLLP_beta", gLLP_beta, &b_gLLP_beta);
   fChain->SetBranchAddress("gLLP_e", gLLP_e, &b_gLLP_e);
   fChain->SetBranchAddress("gLLP_pt", gLLP_pt, &b_gLLP_pt);
   fChain->SetBranchAddress("gLLP_eta", gLLP_eta, &b_gLLP_eta);
   fChain->SetBranchAddress("gLLP_phi", gLLP_phi, &b_gLLP_phi);
   fChain->SetBranchAddress("gLLP_travel_time", gLLP_travel_time, &b_gLLP_travel_time);
   fChain->SetBranchAddress("gLLP_daughter_travel_time", gLLP_daughter_travel_time, &b_gLLP_daughter_travel_time);
   fChain->SetBranchAddress("gLLP_daughter_id", gLLP_daughter_id, &b_gLLP_daughter_id);
   fChain->SetBranchAddress("gLLP_daughter_pt", gLLP_daughter_pt, &b_gLLP_daughter_pt);
   fChain->SetBranchAddress("gLLP_daughter_eta", gLLP_daughter_eta, &b_gLLP_daughter_eta);
   fChain->SetBranchAddress("gLLP_daughter_phi", gLLP_daughter_phi, &b_gLLP_daughter_phi);
   fChain->SetBranchAddress("gLLP_daughter_eta_ecalcorr", gLLP_daughter_eta_ecalcorr, &b_gLLP_daughter_eta_ecalcorr);
   fChain->SetBranchAddress("gLLP_daughter_phi_ecalcorr", gLLP_daughter_phi_ecalcorr, &b_gLLP_daughter_phi_ecalcorr);
   fChain->SetBranchAddress("gLLP_daughter_e", gLLP_daughter_e, &b_gLLP_daughter_e);
   fChain->SetBranchAddress("photon_travel_time", photon_travel_time, &b_photon_travel_time);
   fChain->SetBranchAddress("gen_time", gen_time, &b_gen_time);
   fChain->SetBranchAddress("gen_time_pv", gen_time_pv, &b_gen_time_pv);
   fChain->SetBranchAddress("gLLP_min_delta_r_match_calojet", gLLP_min_delta_r_match_calojet, &b_gLLP_min_delta_r_match_calojet);
   fChain->SetBranchAddress("gLLP_daughter_match_calojet_index", gLLP_daughter_match_calojet_index, &b_gLLP_daughter_match_calojet_index);
   fChain->SetBranchAddress("gLLP_daughter_match_jet_index", gLLP_daughter_match_jet_index, &b_gLLP_daughter_match_jet_index);
   fChain->SetBranchAddress("gLLP_min_delta_r_match_jet", gLLP_min_delta_r_match_jet, &b_gLLP_min_delta_r_match_jet);
   Notify();
}

Bool_t llpntupler_event::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void llpntupler_event::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t llpntupler_event::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef llpntupler_event_cxx
