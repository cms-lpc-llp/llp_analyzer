//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun  7 21:27:17 2017 by ROOT version 6.02/05
// from TTree RazorEvents/selected miniAOD information
// found on file: root://eoscms//eos/cms/store/group/phys_susy/razor/run2/Run2RazorNtupleV3.14/MC_Summer16_EcalRechits/RunIISpring16/v2/sixie/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Run2RazorNtuplerV3p14_ToCERN_MC_Summer16_EcalRechits_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_v2_v1/170607_051412/0000/razorNtuple_1.root
//////////////////////////////////////////////////////////

#ifndef RazorEvents_h
#define RazorEvents_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class RazorEvents {
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
   Int_t           BunchXing[20];   //[nBunchXing]
   Int_t           nPU[20];   //[nBunchXing]
   Float_t         nPUmean[20];   //[nBunchXing]
   Int_t           nMuons;
   Float_t         muonE[700];   //[nMuons]
   Float_t         muonPt[700];   //[nMuons]
   Float_t         muonEta[700];   //[nMuons]
   Float_t         muonPhi[700];   //[nMuons]
   Int_t           muonCharge[700];   //[nMuons]
   Bool_t          muonIsLoose[700];   //[nMuons]
   Bool_t          muonIsMedium[700];   //[nMuons]
   Bool_t          muonIsTight[700];   //[nMuons]
   Float_t         muon_d0[700];   //[nMuons]
   Float_t         muon_dZ[700];   //[nMuons]
   Float_t         muon_ip3d[700];   //[nMuons]
   Float_t         muon_ip3dSignificance[700];   //[nMuons]
   UInt_t          muonType[700];   //[nMuons]
   UInt_t          muonQuality[700];   //[nMuons]
   Float_t         muon_pileupIso[700];   //[nMuons]
   Float_t         muon_chargedIso[700];   //[nMuons]
   Float_t         muon_photonIso[700];   //[nMuons]
   Float_t         muon_neutralHadIso[700];   //[nMuons]
   Float_t         muon_ptrel[700];   //[nMuons]
   Float_t         muon_chargedMiniIso[700];   //[nMuons]
   Float_t         muon_photonAndNeutralHadronMiniIso[700];   //[nMuons]
   Float_t         muon_chargedPileupMiniIso[700];   //[nMuons]
   Float_t         muon_activityMiniIsoAnnulus[700];   //[nMuons]
   Bool_t          muon_passSingleMuTagFilter[700];   //[nMuons]
   Bool_t          muon_passHLTFilter[700][100];   //[nMuons]
   Float_t         muon_validFractionTrackerHits[700];   //[nMuons]
   Bool_t          muon_isGlobal[700];   //[nMuons]
   Float_t         muon_normChi2[700];   //[nMuons]
   Float_t         muon_chi2LocalPosition[700];   //[nMuons]
   Float_t         muon_kinkFinder[700];   //[nMuons]
   Float_t         muon_segmentCompatability[700];   //[nMuons]
   Bool_t          muonIsICHEPMedium[700];   //[nMuons]
   Int_t           nElectrons;
   Float_t         eleE[700];   //[nElectrons]
   Float_t         elePt[700];   //[nElectrons]
   Float_t         eleEta[700];   //[nElectrons]
   Float_t         elePhi[700];   //[nElectrons]
   Float_t         eleCharge[700];   //[nElectrons]
   Float_t         eleEta_SC[700];   //[nElectrons]
   Float_t         eleSigmaIetaIeta[700];   //[nElectrons]
   Float_t         eleFull5x5SigmaIetaIeta[700];   //[nElectrons]
   Float_t         eleR9[700];   //[nElectrons]
   Float_t         ele_dEta[700];   //[nElectrons]
   Float_t         ele_dPhi[700];   //[nElectrons]
   Float_t         ele_HoverE[700];   //[nElectrons]
   Float_t         ele_d0[700];   //[nElectrons]
   Float_t         ele_dZ[700];   //[nElectrons]
   Float_t         ele_ip3d[700];   //[nElectrons]
   Float_t         ele_ip3dSignificance[700];   //[nElectrons]
   Float_t         ele_pileupIso[700];   //[nElectrons]
   Float_t         ele_chargedIso[700];   //[nElectrons]
   Float_t         ele_photonIso[700];   //[nElectrons]
   Float_t         ele_neutralHadIso[700];   //[nElectrons]
   Int_t           ele_MissHits[700];   //[nElectrons]
   Bool_t          ele_PassConvVeto[700];   //[nElectrons]
   Float_t         ele_OneOverEminusOneOverP[700];   //[nElectrons]
   Float_t         ele_IDMVAGeneralPurpose[700];   //[nElectrons]
   Int_t           ele_IDMVACategoryGeneralPurpose[700];   //[nElectrons]
   Float_t         ele_IDMVAHZZ[700];   //[nElectrons]
   Int_t           ele_IDMVACategoryHZZ[700];   //[nElectrons]
   Float_t         ele_RegressionE[700];   //[nElectrons]
   Float_t         ele_CombineP4[700];   //[nElectrons]
   Float_t         ele_ptrel[700];   //[nElectrons]
   Float_t         ele_chargedMiniIso[700];   //[nElectrons]
   Float_t         ele_photonAndNeutralHadronMiniIso[700];   //[nElectrons]
   Float_t         ele_chargedPileupMiniIso[700];   //[nElectrons]
   Float_t         ele_activityMiniIsoAnnulus[700];   //[nElectrons]
   Bool_t          ele_passSingleEleTagFilter[700];   //[nElectrons]
   Bool_t          ele_passTPOneTagFilter[700];   //[nElectrons]
   Bool_t          ele_passTPTwoTagFilter[700];   //[nElectrons]
   Bool_t          ele_passTPOneProbeFilter[700];   //[nElectrons]
   Bool_t          ele_passTPTwoProbeFilter[700];   //[nElectrons]
   Bool_t          ele_passHLTFilter[700][100];   //[nElectrons]
   std::vector<std::vector<unsigned int> > *ele_EcalRechitIndex;
   std::vector<unsigned int> *ele_SeedRechitIndex;
   Int_t           nTaus;
   Float_t         tauE[700];   //[nTaus]
   Float_t         tauPt[700];   //[nTaus]
   Float_t         tauEta[700];   //[nTaus]
   Float_t         tauPhi[700];   //[nTaus]
   Bool_t          tau_IsLoose[700];   //[nTaus]
   Bool_t          tau_IsMedium[700];   //[nTaus]
   Bool_t          tau_IsTight[700];   //[nTaus]
   Bool_t          tau_passEleVetoLoose[700];   //[nTaus]
   Bool_t          tau_passEleVetoMedium[700];   //[nTaus]
   Bool_t          tau_passEleVetoTight[700];   //[nTaus]
   Bool_t          tau_passMuVetoLoose[700];   //[nTaus]
   Bool_t          tau_passMuVetoMedium[700];   //[nTaus]
   Bool_t          tau_passMuVetoTight[700];   //[nTaus]
   UInt_t          tau_ID[700];   //[nTaus]
   Float_t         tau_combinedIsoDeltaBetaCorr3Hits[700];   //[nTaus]
   Float_t         tau_chargedIsoPtSum[700];   //[nTaus]
   Float_t         tau_neutralIsoPtSum[700];   //[nTaus]
   Float_t         tau_puCorrPtSum[700];   //[nTaus]
   Float_t         tau_eleVetoMVA[700];   //[nTaus]
   Int_t           tau_eleVetoCategory[700];   //[nTaus]
   Float_t         tau_muonVetoMVA[700];   //[nTaus]
   Float_t         tau_isoMVAnewDMwLT[700];   //[nTaus]
   Float_t         tau_isoMVAnewDMwoLT[700];   //[nTaus]
   Float_t         tau_leadCandPt[700];   //[nTaus]
   Int_t           tau_leadCandID[700];   //[nTaus]
   Float_t         tau_leadChargedHadrCandPt[700];   //[nTaus]
   Int_t           tau_leadChargedHadrCandID[700];   //[nTaus]
   Int_t           nPhotons;
   Int_t           nPhotons_overlap;
   Float_t         phoE[700];   //[nPhotons]
   Float_t         phoPt[700];   //[nPhotons]
   Float_t         phoEta[700];   //[nPhotons]
   Float_t         phoPhi[700];   //[nPhotons]
   Float_t         phoSigmaIetaIeta[700];   //[nPhotons]
   Float_t         phoFull5x5SigmaIetaIeta[700];   //[nPhotons]
   Float_t         phoR9[700];   //[nPhotons]
   Float_t         pho_HoverE[700];   //[nPhotons]
   Float_t         pho_sumChargedHadronPtAllVertices[700][200];   //[nPhotons]
   Float_t         pho_sumChargedHadronPt[700];   //[nPhotons]
   Float_t         pho_sumNeutralHadronEt[700];   //[nPhotons]
   Float_t         pho_sumPhotonEt[700];   //[nPhotons]
   Float_t         pho_ecalPFClusterIso[700];   //[nPhotons]
   Float_t         pho_hcalPFClusterIso[700];   //[nPhotons]
   Float_t         pho_trkSumPtHollowConeDR03[700];   //[nPhotons]
   Float_t         pho_sumWorstVertexChargedHadronPt[700];   //[nPhotons]
   Float_t         pho_pfIsoChargedHadronIso[700];   //[nPhotons]
   Float_t         pho_pfIsoChargedHadronIsoWrongVtx[700];   //[nPhotons]
   Float_t         pho_pfIsoNeutralHadronIso[700];   //[nPhotons]
   Float_t         pho_pfIsoPhotonIso[700];   //[nPhotons]
   Float_t         pho_pfIsoModFrixione[700];   //[nPhotons]
   Float_t         pho_pfIsoSumPUPt[700];   //[nPhotons]
   Bool_t          pho_isConversion[700];   //[nPhotons]
   Bool_t          pho_passEleVeto[700];   //[nPhotons]
   Float_t         pho_RegressionE[700];   //[nPhotons]
   Float_t         pho_RegressionEUncertainty[700];   //[nPhotons]
   Float_t         pho_IDMVA[700];   //[nPhotons]
   Float_t         pho_superClusterEnergy[700];   //[nPhotons]
   Float_t         pho_superClusterRawEnergy[700];   //[nPhotons]
   Float_t         pho_superClusterEta[700];   //[nPhotons]
   Float_t         pho_superClusterPhi[700];   //[nPhotons]
   Float_t         pho_superClusterX[700];   //[nPhotons]
   Float_t         pho_superClusterY[700];   //[nPhotons]
   Float_t         pho_superClusterZ[700];   //[nPhotons]
   Bool_t          pho_hasPixelSeed[700];   //[nPhotons]
   Bool_t          pho_isStandardPhoton[700];   //[nPhotons]
   Bool_t          pho_passHLTFilter[700][100];   //[nPhotons]
   Int_t           pho_convType[700];   //[nPhotons]
   Float_t         pho_convTrkZ[700];   //[nPhotons]
   Float_t         pho_convTrkClusZ[700];   //[nPhotons]
   Float_t         pho_vtxSumPx[700][200];   //[nPhotons]
   Float_t         pho_vtxSumPy[700][200];   //[nPhotons]
   Float_t         pho_seedRecHitSwitchToGain6[700];   //[nPhotons]
   Float_t         pho_seedRecHitSwitchToGain1[700];   //[nPhotons]
   Float_t         pho_anyRecHitSwitchToGain6[700];   //[nPhotons]
   Float_t         pho_anyRecHitSwitchToGain1[700];   //[nPhotons]
   std::vector<std::vector<unsigned int> > *pho_EcalRechitIndex;
   std::vector<unsigned int> *pho_SeedRechitIndex;
   Int_t           nJets;
   Float_t         jetE[900];   //[nJets]
   Float_t         jetPt[900];   //[nJets]
   Float_t         jetEta[900];   //[nJets]
   Float_t         jetPhi[900];   //[nJets]
   Float_t         jetCSV[900];   //[nJets]
   Float_t         jetCISV[900];   //[nJets]
   Float_t         jetMass[900];   //[nJets]
   Float_t         jetJetArea[900];   //[nJets]
   Float_t         jetPileupE[900];   //[nJets]
   Float_t         jetPileupId[900];   //[nJets]
   Int_t           jetPileupIdFlag[900];   //[nJets]
   Bool_t          jetPassIDLoose[900];   //[nJets]
   Bool_t          jetPassIDTight[900];   //[nJets]
   Bool_t          jetPassMuFrac[900];   //[nJets]
   Bool_t          jetPassEleFrac[900];   //[nJets]
   Int_t           jetPartonFlavor[900];   //[nJets]
   Int_t           jetHadronFlavor[900];   //[nJets]
   Float_t         jetChargedEMEnergyFraction[900];   //[nJets]
   Float_t         jetNeutralEMEnergyFraction[900];   //[nJets]
   Float_t         jetChargedHadronEnergyFraction[900];   //[nJets]
   Float_t         jetNeutralHadronEnergyFraction[900];   //[nJets]
   Float_t         jetMuonEnergyFraction[900];   //[nJets]
   Float_t         jetHOEnergyFraction[900];   //[nJets]
   Float_t         jetHFHadronEnergyFraction[900];   //[nJets]
   Float_t         jetHFEMEnergyFraction[900];   //[nJets]
   Float_t         jetAllMuonPt[900];   //[nJets]
   Float_t         jetAllMuonEta[900];   //[nJets]
   Float_t         jetAllMuonPhi[900];   //[nJets]
   Float_t         jetAllMuonM[900];   //[nJets]
   Float_t         jetPtWeightedDZ[900];   //[nJets]
   UInt_t          nFatJets;
   Float_t         fatJetE[900];
   Float_t         fatJetPt[900];
   Float_t         fatJetEta[900];
   Float_t         fatJetPhi[900];
   Float_t         fatJetCorrectedPt[900];
   Float_t         fatJetPrunedM[900];
   Float_t         fatJetTrimmedM[900];
   Float_t         fatJetFilteredM[900];
   Float_t         fatJetSoftDropM[900];
   Float_t         fatJetCorrectedSoftDropM[900];
   Float_t         fatJetUncorrectedSoftDropM[900];
   Float_t         fatJetTau1[900];
   Float_t         fatJetTau2[900];
   Float_t         fatJetTau3[900];
   Float_t         fatJetMaxSubjetCSV[900];
   Bool_t          fatJetPassIDLoose[900];
   Bool_t          fatJetPassIDTight[900];
   Float_t         metPt;
   Float_t         metPhi;
   Float_t         sumMET;
   Float_t         metType0Pt;
   Float_t         metType0Phi;
   Float_t         metType1Pt;
   Float_t         metType1Pt_raw;
   Float_t         metType1Phi;
   Float_t         metType1Phi_raw;
   Float_t         metType0Plus1Pt;
   Float_t         metType0Plus1Phi;
   Float_t         metEGCleanPt;
   Float_t         metEGCleanPhi;
   Float_t         metMuEGCleanPt;
   Float_t         metMuEGCleanPhi;
   Float_t         metMuEGCleanCorrPt;
   Float_t         metMuEGCleanCorrPhi;
   Float_t         metUncorrectedPt;
   Float_t         metUncorrectedPhi;
   Float_t         metNoHFPt;
   Float_t         metNoHFPhi;
   Float_t         metPuppiPt;
   Float_t         metPuppiPhi;
   Float_t         metCaloPt;
   Float_t         metCaloPhi;
   Float_t         metType1PtJetResUp;
   Float_t         metType1PtJetResDown;
   Float_t         metType1PtJetEnUp;
   Float_t         metType1PtJetEnDown;
   Float_t         metType1PtMuonEnUp;
   Float_t         metType1PtMuonEnDown;
   Float_t         metType1PtElectronEnUp;
   Float_t         metType1PtElectronEnDown;
   Float_t         metType1PtTauEnUp;
   Float_t         metType1PtTauEnDown;
   Float_t         metType1PtUnclusteredEnUp;
   Float_t         metType1PtUnclusteredEnDown;
   Float_t         metType1PtPhotonEnUp;
   Float_t         metType1PtPhotonEnDown;
   Float_t         metType1PhiJetResUp;
   Float_t         metType1PhiJetResDown;
   Float_t         metType1PhiJetEnUp;
   Float_t         metType1PhiJetEnDown;
   Float_t         metType1PhiMuonEnUp;
   Float_t         metType1PhiMuonEnDown;
   Float_t         metType1PhiElectronEnUp;
   Float_t         metType1PhiElectronEnDown;
   Float_t         metType1PhiTauEnUp;
   Float_t         metType1PhiTauEnDown;
   Float_t         metType1PhiUnclusteredEnUp;
   Float_t         metType1PhiUnclusteredEnDown;
   Float_t         metType1PhiPhotonEnUp;
   Float_t         metType1PhiPhotonEnDown;
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHETightNoiseFilter;
   Bool_t          Flag_HBHEIsoNoiseFilter;
   Bool_t          Flag_badChargedCandidateFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_badMuonFilter;
   Bool_t          Flag_BadPFMuonFilter;
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
   Bool_t          Flag_METFilters;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          HLTDecision[300];
   Int_t           HLTPrescale[300];
   Float_t         HLTMR;
   Float_t         HLTRSQ;
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
   Int_t           nGenJets;
   Float_t         genJetE[50];   //[nGenJets]
   Float_t         genJetPt[50];   //[nGenJets]
   Float_t         genJetEta[50];   //[nGenJets]
   Float_t         genJetPhi[50];   //[nGenJets]
   Float_t         genMetPt;
   Float_t         genMetPhi;
   Float_t         genVertexX;
   Float_t         genVertexY;
   Float_t         genVertexZ;
   Float_t	   genVertexT;
   Float_t         genWeight;
   UInt_t          genSignalProcessID;
   Float_t         genQScale;
   Float_t         genAlphaQCD;
   Float_t         genAlphaQED;
   std::string     *lheComments;
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
   Float_t         gParticlePx[4000];   //[nGenParticle]
   Float_t         gParticlePy[4000];   //[nGenParticle]
   Float_t         gParticlePz[4000];   //[nGenParticle]
   Float_t         gParticleEta[4000];   //[nGenParticle]
   Float_t         gParticlePhi[4000];   //[nGenParticle]
   Float_t         gParticleDecayVertexX[4000];   //[nGenParticle]
   Float_t         gParticleDecayVertexY[4000];   //[nGenParticle]
   Float_t         gParticleDecayVertexZ[4000];   //[nGenParticle]

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
   TBranch        *b_ele_passSingleEleTagFilter;   //!
   TBranch        *b_ele_passTPOneTagFilter;   //!
   TBranch        *b_ele_passTPTwoTagFilter;   //!
   TBranch        *b_ele_passTPOneProbeFilter;   //!
   TBranch        *b_ele_passTPTwoProbeFilter;   //!
   TBranch        *b_ele_passHLTFilter;   //!
   TBranch        *b_ele_EcalRechitIndex;   //!
   TBranch        *b_ele_SeedRechitIndex;   //!
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
   TBranch        *b_nPhotons;   //!
   TBranch        *b_nPhotons_overlap;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoPt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoSigmaIetaIeta;   //!
   TBranch        *b_phoFull5x5SigmaIetaIeta;   //!
   TBranch        *b_phoR9;   //!
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
   TBranch        *b_pho_isStandardPhoton;   //!
   TBranch        *b_pho_passHLTFilter;   //!
   TBranch        *b_pho_convType;   //!
   TBranch        *b_pho_convTrkZ;   //!
   TBranch        *b_pho_convTrkClusZ;   //!
   TBranch        *b_pho_vtxSumPx;   //!
   TBranch        *b_pho_vtxSumPy;   //!
   TBranch        *b_pho_seedRecHitSwitchToGain6;   //!
   TBranch        *b_pho_seedRecHitSwitchToGain1;   //!
   TBranch        *b_pho_anyRecHitSwitchToGain6;   //!
   TBranch        *b_pho_anyRecHitSwitchToGain1;   //!
   TBranch        *b_pho_EcalRechitIndex;   //!
   TBranch        *b_pho_SeedRechitIndex;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetCSV;   //!
   TBranch        *b_jetCISV;   //!
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
   TBranch        *b_nFatJets;
   TBranch        *b_fatJetE;
   TBranch        *b_fatJetPt;
   TBranch        *b_fatJetEta;
   TBranch        *b_fatJetPhi;
   TBranch        *b_fatJetCorrectedPt;
   TBranch        *b_fatJetPrunedM;
   TBranch        *b_fatJetTrimmedM;
   TBranch        *b_fatJetFilteredM;
   TBranch        *b_fatJetSoftDropM;
   TBranch        *b_fatJetCorrectedSoftDropM;
   TBranch        *b_fatJetUncorrectedSoftDropM;
   TBranch        *b_fatJetTau1;
   TBranch        *b_fatJetTau2;
   TBranch        *b_fatJetTau3;
   TBranch        *b_fatJetMaxSubjetCSV;
   TBranch        *b_fatJetPassIDLoose;
   TBranch        *b_fatJetPassIDTight;
   TBranch        *b_metPt;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_sumMET;   //!
   TBranch        *b_metType0Pt;   //!
   TBranch        *b_metType0Phi;   //!
   TBranch        *b_metType1Pt;   //!
   TBranch        *b_metType1Pt_raw;   //!
   TBranch        *b_metType1Phi;   //!
   TBranch        *b_metType1Phi_raw;   //!
   TBranch        *b_metType0Plus1Pt;   //!
   TBranch        *b_metType0Plus1Phi;   //!
   TBranch        *b_metEGCleanPt;   //!
   TBranch        *b_metEGCleanPhi;   //!
   TBranch        *b_metMuEGCleanPt;   //!
   TBranch        *b_metMuEGCleanPhi;   //!
   TBranch        *b_metMuEGCleanCorrPt;   //!
   TBranch        *b_metMuEGCleanCorrPhi;   //!
   TBranch        *b_metUncorrectedPt;   //!
   TBranch        *b_metUncorrectedPhi;   //!
   TBranch        *b_metNoHFPt;   //!
   TBranch        *b_metNoHFPhi;   //!
   TBranch        *b_metPuppiPt;   //!
   TBranch        *b_metPuppiPhi;   //!
   TBranch        *b_metCaloPt;   //!
   TBranch        *b_metCaloPhi;   //!
   TBranch        *b_metType1PtJetResUp;   //!
   TBranch        *b_metType1PtJetResDown;   //!
   TBranch        *b_metType1PtJetEnUp;   //!
   TBranch        *b_metType1PtJetEnDown;   //!
   TBranch        *b_metType1PtMuonEnUp;   //!
   TBranch        *b_metType1PtMuonEnDown;   //!
   TBranch        *b_metType1PtElectronEnUp;   //!
   TBranch        *b_metType1PtElectronEnDown;   //!
   TBranch        *b_metType1PtTauEnUp;   //!
   TBranch        *b_metType1PtTauEnDown;   //!
   TBranch        *b_metType1PtUnclusteredEnUp;   //!
   TBranch        *b_metType1PtUnclusteredEnDown;   //!
   TBranch        *b_metType1PtPhotonEnUp;   //!
   TBranch        *b_metType1PtPhotonEnDown;   //!
   TBranch        *b_metType1PhiJetResUp;   //!
   TBranch        *b_metType1PhiJetResDown;   //!
   TBranch        *b_metType1PhiJetEnUp;   //!
   TBranch        *b_metType1PhiJetEnDown;   //!
   TBranch        *b_metType1PhiMuonEnUp;   //!
   TBranch        *b_metType1PhiMuonEnDown;   //!
   TBranch        *b_metType1PhiElectronEnUp;   //!
   TBranch        *b_metType1PhiElectronEnDown;   //!
   TBranch        *b_metType1PhiTauEnUp;   //!
   TBranch        *b_metType1PhiTauEnDown;   //!
   TBranch        *b_metType1PhiUnclusteredEnUp;   //!
   TBranch        *b_metType1PhiUnclusteredEnDown;   //!
   TBranch        *b_metType1PhiPhotonEnUp;   //!
   TBranch        *b_metType1PhiPhotonEnDown;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHETightNoiseFilter;   //!
   TBranch        *b_Flag_HBHEIsoNoiseFilter;   //!
   TBranch        *b_Flag_badChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_badMuonFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
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
   TBranch        *b_Flag_METFilters;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_HLTDecision;   //!
   TBranch        *b_HLTPrescale;   //!
   TBranch        *b_HLTMR;   //!
   TBranch        *b_HLTRSQ;   //!
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
   TBranch        *b_nGenJets;   //!
   TBranch        *b_genJetE;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genMetPt;   //!
   TBranch        *b_genMetPhi;   //!
   TBranch        *b_genVertexX;   //!
   TBranch        *b_genVertexY;   //!
   TBranch        *b_genVertexZ;   //!
   TBranch	  *b_genVertexT;
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
   TBranch        *b_gParticlePx;   //!
   TBranch        *b_gParticlePy;   //!
   TBranch        *b_gParticlePz;   //!
   TBranch        *b_gParticleEta;   //!
   TBranch        *b_gParticlePhi;   //!
   TBranch        *b_gParticleDecayVertexX;   //!
   TBranch        *b_gParticleDecayVertexY;   //!
   TBranch        *b_gParticleDecayVertexZ;   //!

   RazorEvents(TTree *tree=0);
   virtual ~RazorEvents();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RazorEvents_cxx
RazorEvents::RazorEvents(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/phys_susy/razor/run2/Run2RazorNtupleV3.14/MC_Summer16_EcalRechits/RunIISpring16/v2/sixie/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Run2RazorNtuplerV3p14_ToCERN_MC_Summer16_EcalRechits_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_v2_v1/170607_051412/0000/razorNtuple_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//eos/cms/store/group/phys_susy/razor/run2/Run2RazorNtupleV3.14/MC_Summer16_EcalRechits/RunIISpring16/v2/sixie/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Run2RazorNtuplerV3p14_ToCERN_MC_Summer16_EcalRechits_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_v2_v1/170607_051412/0000/razorNtuple_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://eoscms//eos/cms/store/group/phys_susy/razor/run2/Run2RazorNtupleV3.14/MC_Summer16_EcalRechits/RunIISpring16/v2/sixie/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Run2RazorNtuplerV3p14_ToCERN_MC_Summer16_EcalRechits_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_v2_v1/170607_051412/0000/razorNtuple_1.root:/ntuples");
      dir->GetObject("RazorEvents",tree);

   }
   Init(tree);
}

RazorEvents::~RazorEvents()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RazorEvents::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RazorEvents::LoadTree(Long64_t entry)
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

void RazorEvents::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   lheComments = 0;
   ele_EcalRechitIndex = 0;
   ele_SeedRechitIndex = 0;
   pho_EcalRechitIndex = 0;
   pho_SeedRechitIndex = 0;
   ecalRechit_Eta = 0;
   ecalRechit_Phi = 0;
   ecalRechit_X = 0;
   ecalRechit_Y = 0;
   ecalRechit_Z = 0;
   ecalRechit_E = 0;
   ecalRechit_T = 0;
   ecalRechit_ID = 0;
   ecalRechit_FlagOOT = 0;
   ecalRechit_GainSwitch1 = 0;
   ecalRechit_GainSwitch6 = 0;
   ecalRechit_transpCorr = 0;
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
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("fixedGridRhoFastjetAllCalo", &fixedGridRhoFastjetAllCalo, &b_fixedGridRhoFastjetAllCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, &b_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("nPVAll", &nPVAll, &b_nPVAll);
   fChain->SetBranchAddress("pvAllX", pvAllX, &b_pvAllX);
   fChain->SetBranchAddress("pvAllY", pvAllY, &b_pvAllY);
   fChain->SetBranchAddress("pvAllZ", pvAllZ, &b_pvAllZ);
   fChain->SetBranchAddress("pvAllLogSumPtSq", pvAllLogSumPtSq, &b_pvAllLogSumPtSq);
   fChain->SetBranchAddress("pvAllSumPx", pvAllSumPx, &b_pvAllSumPx);
   fChain->SetBranchAddress("pvAllSumPy", pvAllSumPy, &b_pvAllSumPy);
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
   fChain->SetBranchAddress("muon_chargedMiniIso", muon_chargedMiniIso, &b_muon_chargedMiniIso);
   fChain->SetBranchAddress("muon_photonAndNeutralHadronMiniIso", muon_photonAndNeutralHadronMiniIso, &b_muon_photonAndNeutralHadronMiniIso);
   fChain->SetBranchAddress("muon_chargedPileupMiniIso", muon_chargedPileupMiniIso, &b_muon_chargedPileupMiniIso);
   fChain->SetBranchAddress("muon_activityMiniIsoAnnulus", muon_activityMiniIsoAnnulus, &b_muon_activityMiniIsoAnnulus);
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
   fChain->SetBranchAddress("ele_RegressionE", ele_RegressionE, &b_ele_RegressionE);
   fChain->SetBranchAddress("ele_CombineP4", ele_CombineP4, &b_ele_CombineP4);
   fChain->SetBranchAddress("ele_ptrel", ele_ptrel, &b_ele_ptrel);
   fChain->SetBranchAddress("ele_chargedMiniIso", ele_chargedMiniIso, &b_ele_chargedMiniIso);
   fChain->SetBranchAddress("ele_photonAndNeutralHadronMiniIso", ele_photonAndNeutralHadronMiniIso, &b_ele_photonAndNeutralHadronMiniIso);
   fChain->SetBranchAddress("ele_chargedPileupMiniIso", ele_chargedPileupMiniIso, &b_ele_chargedPileupMiniIso);
   fChain->SetBranchAddress("ele_activityMiniIsoAnnulus", ele_activityMiniIsoAnnulus, &b_ele_activityMiniIsoAnnulus);
   fChain->SetBranchAddress("ele_passSingleEleTagFilter", ele_passSingleEleTagFilter, &b_ele_passSingleEleTagFilter);
   fChain->SetBranchAddress("ele_passTPOneTagFilter", ele_passTPOneTagFilter, &b_ele_passTPOneTagFilter);
   fChain->SetBranchAddress("ele_passTPTwoTagFilter", ele_passTPTwoTagFilter, &b_ele_passTPTwoTagFilter);
   fChain->SetBranchAddress("ele_passTPOneProbeFilter", ele_passTPOneProbeFilter, &b_ele_passTPOneProbeFilter);
   fChain->SetBranchAddress("ele_passTPTwoProbeFilter", ele_passTPTwoProbeFilter, &b_ele_passTPTwoProbeFilter);
   fChain->SetBranchAddress("ele_passHLTFilter", ele_passHLTFilter, &b_ele_passHLTFilter);
   fChain->SetBranchAddress("ele_EcalRechitIndex", &ele_EcalRechitIndex, &b_ele_EcalRechitIndex);
   fChain->SetBranchAddress("ele_SeedRechitIndex", &ele_SeedRechitIndex, &b_ele_SeedRechitIndex);
   fChain->SetBranchAddress("nTaus", &nTaus, &b_nTaus);
   fChain->SetBranchAddress("tauE", tauE, &b_tauE);
   fChain->SetBranchAddress("tauPt", tauPt, &b_tauPt);
   fChain->SetBranchAddress("tauEta", tauEta, &b_tauEta);
   fChain->SetBranchAddress("tauPhi", tauPhi, &b_tauPhi);
   fChain->SetBranchAddress("tau_IsLoose", tau_IsLoose, &b_tau_IsLoose);
   fChain->SetBranchAddress("tau_IsMedium", tau_IsMedium, &b_tau_IsMedium);
   fChain->SetBranchAddress("tau_IsTight", tau_IsTight, &b_tau_IsTight);
   fChain->SetBranchAddress("tau_passEleVetoLoose", tau_passEleVetoLoose, &b_tau_passEleVetoLoose);
   fChain->SetBranchAddress("tau_passEleVetoMedium", tau_passEleVetoMedium, &b_tau_passEleVetoMedium);
   fChain->SetBranchAddress("tau_passEleVetoTight", tau_passEleVetoTight, &b_tau_passEleVetoTight);
   fChain->SetBranchAddress("tau_passMuVetoLoose", tau_passMuVetoLoose, &b_tau_passMuVetoLoose);
   fChain->SetBranchAddress("tau_passMuVetoMedium", tau_passMuVetoMedium, &b_tau_passMuVetoMedium);
   fChain->SetBranchAddress("tau_passMuVetoTight", tau_passMuVetoTight, &b_tau_passMuVetoTight);
   fChain->SetBranchAddress("tau_ID", tau_ID, &b_tau_ID);
   fChain->SetBranchAddress("tau_combinedIsoDeltaBetaCorr3Hits", tau_combinedIsoDeltaBetaCorr3Hits, &b_tau_combinedIsoDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_chargedIsoPtSum", tau_chargedIsoPtSum, &b_tau_chargedIsoPtSum);
   fChain->SetBranchAddress("tau_neutralIsoPtSum", tau_neutralIsoPtSum, &b_tau_neutralIsoPtSum);
   fChain->SetBranchAddress("tau_puCorrPtSum", tau_puCorrPtSum, &b_tau_puCorrPtSum);
   fChain->SetBranchAddress("tau_eleVetoMVA", tau_eleVetoMVA, &b_tau_eleVetoMVA);
   fChain->SetBranchAddress("tau_eleVetoCategory", tau_eleVetoCategory, &b_tau_eleVetoCategory);
   fChain->SetBranchAddress("tau_muonVetoMVA", tau_muonVetoMVA, &b_tau_muonVetoMVA);
   fChain->SetBranchAddress("tau_isoMVAnewDMwLT", tau_isoMVAnewDMwLT, &b_tau_isoMVAnewDMwLT);
   fChain->SetBranchAddress("tau_isoMVAnewDMwoLT", tau_isoMVAnewDMwoLT, &b_tau_isoMVAnewDMwoLT);
   fChain->SetBranchAddress("tau_leadCandPt", tau_leadCandPt, &b_tau_leadCandPt);
   fChain->SetBranchAddress("tau_leadCandID", tau_leadCandID, &b_tau_leadCandID);
   fChain->SetBranchAddress("tau_leadChargedHadrCandPt", tau_leadChargedHadrCandPt, &b_tau_leadChargedHadrCandPt);
   fChain->SetBranchAddress("tau_leadChargedHadrCandID", tau_leadChargedHadrCandID, &b_tau_leadChargedHadrCandID);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("nPhotons_overlap", &nPhotons_overlap, &b_nPhotons_overlap);
   fChain->SetBranchAddress("phoE", phoE, &b_phoE);
   fChain->SetBranchAddress("phoPt", phoPt, &b_phoPt);
   fChain->SetBranchAddress("phoEta", phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoSigmaIetaIeta", phoSigmaIetaIeta, &b_phoSigmaIetaIeta);
   fChain->SetBranchAddress("phoFull5x5SigmaIetaIeta", phoFull5x5SigmaIetaIeta, &b_phoFull5x5SigmaIetaIeta);
   fChain->SetBranchAddress("phoR9", phoR9, &b_phoR9);
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
   fChain->SetBranchAddress("pho_isStandardPhoton", pho_isStandardPhoton, &b_pho_isStandardPhoton);
   fChain->SetBranchAddress("pho_passHLTFilter", pho_passHLTFilter, &b_pho_passHLTFilter);
   fChain->SetBranchAddress("pho_convType", pho_convType, &b_pho_convType);
   fChain->SetBranchAddress("pho_convTrkZ", pho_convTrkZ, &b_pho_convTrkZ);
   fChain->SetBranchAddress("pho_convTrkClusZ", pho_convTrkClusZ, &b_pho_convTrkClusZ);
   fChain->SetBranchAddress("pho_vtxSumPx", pho_vtxSumPx, &b_pho_vtxSumPx);
   fChain->SetBranchAddress("pho_vtxSumPy", pho_vtxSumPy, &b_pho_vtxSumPy);
   fChain->SetBranchAddress("pho_seedRecHitSwitchToGain6", pho_seedRecHitSwitchToGain6, &b_pho_seedRecHitSwitchToGain6);
   fChain->SetBranchAddress("pho_seedRecHitSwitchToGain1", pho_seedRecHitSwitchToGain1, &b_pho_seedRecHitSwitchToGain1);
   fChain->SetBranchAddress("pho_anyRecHitSwitchToGain6", pho_anyRecHitSwitchToGain6, &b_pho_anyRecHitSwitchToGain6);
   fChain->SetBranchAddress("pho_anyRecHitSwitchToGain1", pho_anyRecHitSwitchToGain1, &b_pho_anyRecHitSwitchToGain1);
   fChain->SetBranchAddress("pho_EcalRechitIndex", &pho_EcalRechitIndex, &b_pho_EcalRechitIndex);
   fChain->SetBranchAddress("pho_SeedRechitIndex", &pho_SeedRechitIndex, &b_pho_SeedRechitIndex);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("jetE", jetE, &b_jetE);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetCSV", jetCSV, &b_jetCSV);
   fChain->SetBranchAddress("jetCISV", jetCISV, &b_jetCISV);
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
   fChain->SetBranchAddress("nFatJets", &nFatJets, &b_nFatJets);
   fChain->SetBranchAddress("fatJetE", fatJetE, &b_fatJetE);
   fChain->SetBranchAddress("fatJetPt", fatJetPt, &b_fatJetPt);
   fChain->SetBranchAddress("fatJetEta", fatJetEta, &b_fatJetEta);
   fChain->SetBranchAddress("fatJetPhi", fatJetPhi, &b_fatJetPhi);
   fChain->SetBranchAddress("fatJetCorrectedPt", fatJetCorrectedPt, &b_fatJetCorrectedPt);
   fChain->SetBranchAddress("fatJetPrunedM", fatJetPrunedM, &b_fatJetPrunedM);
   fChain->SetBranchAddress("fatJetTrimmedM", fatJetTrimmedM, &b_fatJetTrimmedM);
   fChain->SetBranchAddress("fatJetFilteredM", fatJetFilteredM, &b_fatJetFilteredM);
   fChain->SetBranchAddress("fatJetSoftDropM", fatJetSoftDropM, &b_fatJetSoftDropM);
   fChain->SetBranchAddress("fatJetCorrectedSoftDropM", fatJetCorrectedSoftDropM, &b_fatJetCorrectedSoftDropM);
   fChain->SetBranchAddress("fatJetUncorrectedSoftDropM", fatJetUncorrectedSoftDropM, &b_fatJetUncorrectedSoftDropM);
   fChain->SetBranchAddress("fatJetTau1", fatJetTau1, &b_fatJetTau1);
   fChain->SetBranchAddress("fatJetTau2", fatJetTau2, &b_fatJetTau2);
   fChain->SetBranchAddress("fatJetTau3", fatJetTau3, &b_fatJetTau3);
   fChain->SetBranchAddress("fatJetMaxSubjetCSV", fatJetMaxSubjetCSV, &b_fatJetMaxSubjetCSV);
   fChain->SetBranchAddress("fatJetPassIDLoose", fatJetPassIDLoose, &b_fatJetPassIDLoose);
   fChain->SetBranchAddress("fatJetPassIDTight", fatJetPassIDTight, &b_fatJetPassIDTight);
   fChain->SetBranchAddress("metPt", &metPt, &b_metPt);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("sumMET", &sumMET, &b_sumMET);
   fChain->SetBranchAddress("metType0Pt", &metType0Pt, &b_metType0Pt);
   fChain->SetBranchAddress("metType0Phi", &metType0Phi, &b_metType0Phi);
   fChain->SetBranchAddress("metType1Pt", &metType1Pt, &b_metType1Pt);
   fChain->SetBranchAddress("metType1Pt_raw", &metType1Pt_raw, &b_metType1Pt_raw);
   fChain->SetBranchAddress("metType1Phi", &metType1Phi, &b_metType1Phi);
   fChain->SetBranchAddress("metType1Phi_raw", &metType1Phi_raw, &b_metType1Phi_raw);
   fChain->SetBranchAddress("metType0Plus1Pt", &metType0Plus1Pt, &b_metType0Plus1Pt);
   fChain->SetBranchAddress("metType0Plus1Phi", &metType0Plus1Phi, &b_metType0Plus1Phi);
   fChain->SetBranchAddress("metEGCleanPt", &metEGCleanPt, &b_metEGCleanPt);
   fChain->SetBranchAddress("metEGCleanPhi", &metEGCleanPhi, &b_metEGCleanPhi);
   fChain->SetBranchAddress("metMuEGCleanPt", &metMuEGCleanPt, &b_metMuEGCleanPt);
   fChain->SetBranchAddress("metMuEGCleanPhi", &metMuEGCleanPhi, &b_metMuEGCleanPhi);
   fChain->SetBranchAddress("metMuEGCleanCorrPt", &metMuEGCleanCorrPt, &b_metMuEGCleanCorrPt);
   fChain->SetBranchAddress("metMuEGCleanCorrPhi", &metMuEGCleanCorrPhi, &b_metMuEGCleanCorrPhi);
   fChain->SetBranchAddress("metUncorrectedPt", &metUncorrectedPt, &b_metUncorrectedPt);
   fChain->SetBranchAddress("metUncorrectedPhi", &metUncorrectedPhi, &b_metUncorrectedPhi);
   fChain->SetBranchAddress("metNoHFPt", &metNoHFPt, &b_metNoHFPt);
   fChain->SetBranchAddress("metNoHFPhi", &metNoHFPhi, &b_metNoHFPhi);
   fChain->SetBranchAddress("metPuppiPt", &metPuppiPt, &b_metPuppiPt);
   fChain->SetBranchAddress("metPuppiPhi", &metPuppiPhi, &b_metPuppiPhi);
   fChain->SetBranchAddress("metCaloPt", &metCaloPt, &b_metCaloPt);
   fChain->SetBranchAddress("metCaloPhi", &metCaloPhi, &b_metCaloPhi);
   fChain->SetBranchAddress("metType1PtJetResUp", &metType1PtJetResUp, &b_metType1PtJetResUp);
   fChain->SetBranchAddress("metType1PtJetResDown", &metType1PtJetResDown, &b_metType1PtJetResDown);
   fChain->SetBranchAddress("metType1PtJetEnUp", &metType1PtJetEnUp, &b_metType1PtJetEnUp);
   fChain->SetBranchAddress("metType1PtJetEnDown", &metType1PtJetEnDown, &b_metType1PtJetEnDown);
   fChain->SetBranchAddress("metType1PtMuonEnUp", &metType1PtMuonEnUp, &b_metType1PtMuonEnUp);
   fChain->SetBranchAddress("metType1PtMuonEnDown", &metType1PtMuonEnDown, &b_metType1PtMuonEnDown);
   fChain->SetBranchAddress("metType1PtElectronEnUp", &metType1PtElectronEnUp, &b_metType1PtElectronEnUp);
   fChain->SetBranchAddress("metType1PtElectronEnDown", &metType1PtElectronEnDown, &b_metType1PtElectronEnDown);
   fChain->SetBranchAddress("metType1PtTauEnUp", &metType1PtTauEnUp, &b_metType1PtTauEnUp);
   fChain->SetBranchAddress("metType1PtTauEnDown", &metType1PtTauEnDown, &b_metType1PtTauEnDown);
   fChain->SetBranchAddress("metType1PtUnclusteredEnUp", &metType1PtUnclusteredEnUp, &b_metType1PtUnclusteredEnUp);
   fChain->SetBranchAddress("metType1PtUnclusteredEnDown", &metType1PtUnclusteredEnDown, &b_metType1PtUnclusteredEnDown);
   fChain->SetBranchAddress("metType1PtPhotonEnUp", &metType1PtPhotonEnUp, &b_metType1PtPhotonEnUp);
   fChain->SetBranchAddress("metType1PtPhotonEnDown", &metType1PtPhotonEnDown, &b_metType1PtPhotonEnDown);
   fChain->SetBranchAddress("metType1PhiJetResUp", &metType1PhiJetResUp, &b_metType1PhiJetResUp);
   fChain->SetBranchAddress("metType1PhiJetResDown", &metType1PhiJetResDown, &b_metType1PhiJetResDown);
   fChain->SetBranchAddress("metType1PhiJetEnUp", &metType1PhiJetEnUp, &b_metType1PhiJetEnUp);
   fChain->SetBranchAddress("metType1PhiJetEnDown", &metType1PhiJetEnDown, &b_metType1PhiJetEnDown);
   fChain->SetBranchAddress("metType1PhiMuonEnUp", &metType1PhiMuonEnUp, &b_metType1PhiMuonEnUp);
   fChain->SetBranchAddress("metType1PhiMuonEnDown", &metType1PhiMuonEnDown, &b_metType1PhiMuonEnDown);
   fChain->SetBranchAddress("metType1PhiElectronEnUp", &metType1PhiElectronEnUp, &b_metType1PhiElectronEnUp);
   fChain->SetBranchAddress("metType1PhiElectronEnDown", &metType1PhiElectronEnDown, &b_metType1PhiElectronEnDown);
   fChain->SetBranchAddress("metType1PhiTauEnUp", &metType1PhiTauEnUp, &b_metType1PhiTauEnUp);
   fChain->SetBranchAddress("metType1PhiTauEnDown", &metType1PhiTauEnDown, &b_metType1PhiTauEnDown);
   fChain->SetBranchAddress("metType1PhiUnclusteredEnUp", &metType1PhiUnclusteredEnUp, &b_metType1PhiUnclusteredEnUp);
   fChain->SetBranchAddress("metType1PhiUnclusteredEnDown", &metType1PhiUnclusteredEnDown, &b_metType1PhiUnclusteredEnDown);
   fChain->SetBranchAddress("metType1PhiPhotonEnUp", &metType1PhiPhotonEnUp, &b_metType1PhiPhotonEnUp);
   fChain->SetBranchAddress("metType1PhiPhotonEnDown", &metType1PhiPhotonEnDown, &b_metType1PhiPhotonEnDown);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHETightNoiseFilter", &Flag_HBHETightNoiseFilter, &b_Flag_HBHETightNoiseFilter);
   fChain->SetBranchAddress("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, &b_Flag_HBHEIsoNoiseFilter);
   fChain->SetBranchAddress("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter, &b_Flag_badChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_badMuonFilter", &Flag_badMuonFilter, &b_Flag_badMuonFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
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
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("HLTDecision", HLTDecision, &b_HLTDecision);
   fChain->SetBranchAddress("HLTPrescale", HLTPrescale, &b_HLTPrescale);
   fChain->SetBranchAddress("HLTMR", &HLTMR, &b_HLTMR);
   fChain->SetBranchAddress("HLTRSQ", &HLTRSQ, &b_HLTRSQ);
   fChain->SetBranchAddress("ecalRechit_Eta", &ecalRechit_Eta, &b_ecalRechit_Eta);
   fChain->SetBranchAddress("ecalRechit_Phi", &ecalRechit_Phi, &b_ecalRechit_Phi);
   fChain->SetBranchAddress("ecalRechit_X", &ecalRechit_X, &b_ecalRechit_X);
   fChain->SetBranchAddress("ecalRechit_Y", &ecalRechit_Y, &b_ecalRechit_Y);
   fChain->SetBranchAddress("ecalRechit_Z", &ecalRechit_Z, &b_ecalRechit_Z);
   fChain->SetBranchAddress("ecalRechit_E", &ecalRechit_E, &b_ecalRechit_E);
   fChain->SetBranchAddress("ecalRechit_T", &ecalRechit_T, &b_ecalRechit_T);
   fChain->SetBranchAddress("ecalRechit_ID", &ecalRechit_ID, &b_ecalRechit_ID);
   fChain->SetBranchAddress("ecalRechit_FlagOOT", &ecalRechit_FlagOOT, &b_ecalRechit_FlagOOT);
   fChain->SetBranchAddress("ecalRechit_GainSwitch1", &ecalRechit_GainSwitch1, &b_ecalRechit_GainSwitch1);
   fChain->SetBranchAddress("ecalRechit_GainSwitch6", &ecalRechit_GainSwitch6, &b_ecalRechit_GainSwitch6);
   fChain->SetBranchAddress("ecalRechit_transpCorr", &ecalRechit_transpCorr, &b_ecalRechit_transpCorr);
   fChain->SetBranchAddress("nGenJets", &nGenJets, &b_nGenJets);
   fChain->SetBranchAddress("genJetE", genJetE, &b_genJetE);
   fChain->SetBranchAddress("genJetPt", genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetEta", genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPhi", genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genMetPt", &genMetPt, &b_genMetPt);
   fChain->SetBranchAddress("genMetPhi", &genMetPhi, &b_genMetPhi);
   fChain->SetBranchAddress("genVertexX", &genVertexX, &b_genVertexX);
   fChain->SetBranchAddress("genVertexY", &genVertexY, &b_genVertexY);
   fChain->SetBranchAddress("genVertexZ", &genVertexZ, &b_genVertexZ);
   fChain->SetBranchAddress("genVertexT", &genVertexT, &b_genVertexT);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genSignalProcessID", &genSignalProcessID, &b_genSignalProcessID);
   fChain->SetBranchAddress("genQScale", &genQScale, &b_genQScale);
   fChain->SetBranchAddress("genAlphaQCD", &genAlphaQCD, &b_genAlphaQCD);
   fChain->SetBranchAddress("genAlphaQED", &genAlphaQED, &b_genAlphaQED);
   fChain->SetBranchAddress("lheComments", &lheComments, &b_lheComments);   
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
   fChain->SetBranchAddress("gParticleDecayVertexX", gParticleDecayVertexX, &b_gParticleDecayVertexX);
   fChain->SetBranchAddress("gParticleDecayVertexY", gParticleDecayVertexY, &b_gParticleDecayVertexY);
   fChain->SetBranchAddress("gParticleDecayVertexZ", gParticleDecayVertexZ, &b_gParticleDecayVertexZ);
   Notify();
}

Bool_t RazorEvents::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RazorEvents::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RazorEvents::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RazorEvents_cxx
