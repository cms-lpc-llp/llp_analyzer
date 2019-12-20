//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 29 22:19:14 2016 by ROOT version 6.06/01
// from TTree RazorEventsUpgradeTiming/selected miniAOD information
// found on file: /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/run2/Run2RazorNtupleV3.6/MC/RunIISpring16/v1/sixie/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Run2RazorNtuplerV3p6_ToCERN_MC_25ns_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1_v1_v1/160925_172232/0000/razorNtuple_756.root
//////////////////////////////////////////////////////////

#ifndef RazorEventsUpgradeTiming_h
#define RazorEventsUpgradeTiming_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class RazorEventsUpgradeTiming {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Bool_t          isData;
   Int_t           nPV;
   UInt_t          runNum;
   UInt_t          lumiNum;
   UInt_t          eventNum;
   Float_t         pvX;
   Float_t         pvY;
   Float_t         pvZ;
   Float_t         fixedGridRhoAll;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetAllCalo;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   Float_t         beamSpotX=0.0;
   Float_t         beamSpotY=0.0;
   Float_t         beamSpotZ=0.0;
   Int_t           nPVAll;
   Int_t           pvIndex[500];   //[nPVAll]
   Float_t         pvAllX[500];   //[nPVAll]
   Float_t         pvAllY[500];   //[nPVAll]
   Float_t         pvAllZ[500];   //[nPVAll]
   Float_t         pvAllT[500];   //[nPVAll]
   Float_t         pvAllLogSumPtSq[500];   //[nPVAll]
   Float_t         pvAllSumPt[500];   //[nPVAll]
   Float_t         pvAllSumPx[500];   //[nPVAll]
   Int_t         pvNtrack_orig[500];   //[nPVAll]
   Float_t         pvAllSumPy[500];   //[nPVAll]
   Float_t         pvAllLogSumPtSq_dt[500];   //[nPVAll]
   Float_t         pvAllSumPx_dt[500];   //[nPVAll]
   Int_t         pvNtrack_orig_dt[500];   //[nPVAll]
   Float_t         pvAllSumPt_dt[500];   //[nPVAll]
   Float_t         pvAllSumPy_dt[500];   //[nPVAll]

   Int_t           nBunchXing;
   Int_t           BunchXing[20];   //[nBunchXing]
   Int_t           nPU[20];   //[nBunchXing]
   Float_t         nPUmean[20];   //[nBunchXing]
   Int_t           nMuons;
   Float_t         muonE[40];   //[nMuons]
   Float_t         muonPt[40];   //[nMuons]
   Float_t         muonEta[40];   //[nMuons]
   Float_t         muonPhi[40];   //[nMuons]
   Int_t           muonCharge[40];   //[nMuons]
   Bool_t          muonIsLoose[40];   //[nMuons]
   Bool_t          muonIsMedium[40];   //[nMuons]
   Bool_t          muonIsTight[40];   //[nMuons]
   Float_t         muon_d0[40];   //[nMuons]
   Float_t         muon_dZ[40];   //[nMuons]
   Float_t         muon_ip3d[40];   //[nMuons]
   Float_t         muon_ip3dSignificance[40];   //[nMuons]
   UInt_t          muonType[40];   //[nMuons]
   UInt_t          muonQuality[40];   //[nMuons]
   Float_t         muon_pileupIso[40];   //[nMuons]
   Float_t         muon_chargedIso[40];   //[nMuons]
   Float_t         muon_photonIso[40];   //[nMuons]
   Float_t         muon_neutralHadIso[40];   //[nMuons]
   Float_t         muon_ptrel[40];   //[nMuons]
   Float_t         muon_chargedMiniIso[40];   //[nMuons]
   Float_t         muon_photonAndNeutralHadronMiniIso[40];   //[nMuons]
   Float_t         muon_chargedPileupMiniIso[40];   //[nMuons]
   Float_t         muon_activityMiniIsoAnnulus[40];   //[nMuons]
   Bool_t          muon_passSingleMuTagFilter[40];   //[nMuons]
   Bool_t          muon_passHLTFilter[40][100];   //[nMuons]
   Float_t         muon_validFractionTrackerHits[40];   //[nMuons]
   Bool_t          muon_isGlobal[40];   //[nMuons]
   Float_t         muon_normChi2[40];   //[nMuons]
   Float_t         muon_chi2LocalPosition[40];   //[nMuons]
   Float_t         muon_kinkFinder[40];   //[nMuons]
   Float_t         muon_segmentCompatability[40];   //[nMuons]
   Bool_t          muonIsICHEPMedium[40];   //[nMuons]
   Int_t           nElectrons;
   Float_t         eleE[40];   //[nElectrons]
   Float_t         elePt[40];   //[nElectrons]
   Float_t         eleEta[40];   //[nElectrons]
   Float_t         elePhi[40];   //[nElectrons]
   Float_t         eleCharge[40];   //[nElectrons]
   Float_t         eleEta_SC[40];   //[nElectrons]
   Float_t         eleSigmaIetaIeta[40];   //[nElectrons]
   Float_t         eleFull5x5SigmaIetaIeta[40];   //[nElectrons]
   Float_t         eleR9[40];   //[nElectrons]
   Float_t         ele_dEta[40];   //[nElectrons]
   Float_t         ele_dPhi[40];   //[nElectrons]
   Float_t         ele_HoverE[40];   //[nElectrons]
   Float_t         ele_d0[40];   //[nElectrons]
   Float_t         ele_dZ[40];   //[nElectrons]
   Float_t         ele_ip3d[40];   //[nElectrons]
   Float_t         ele_ip3dSignificance[40];   //[nElectrons]
   Float_t         ele_pileupIso[40];   //[nElectrons]
   Float_t         ele_chargedIso[40];   //[nElectrons]
   Float_t         ele_photonIso[40];   //[nElectrons]
   Float_t         ele_neutralHadIso[40];   //[nElectrons]
   Int_t           ele_MissHits[40];   //[nElectrons]
   Bool_t          ele_PassConvVeto[40];   //[nElectrons]
   Float_t         ele_OneOverEminusOneOverP[40];   //[nElectrons]
   Float_t         ele_IDMVATrig[40];   //[nElectrons]
   Float_t         ele_IDMVANonTrig[40];   //[nElectrons]
   Float_t         ele_RegressionE[40];   //[nElectrons]
   Float_t         ele_CombineP4[40];   //[nElectrons]
   Float_t         ele_ptrel[40];   //[nElectrons]
   Float_t         ele_chargedMiniIso[40];   //[nElectrons]
   Float_t         ele_photonAndNeutralHadronMiniIso[40];   //[nElectrons]
   Float_t         ele_chargedPileupMiniIso[40];   //[nElectrons]
   Float_t         ele_activityMiniIsoAnnulus[40];   //[nElectrons]
   Bool_t          ele_passSingleEleTagFilter[40];   //[nElectrons]
   Bool_t          ele_passTPOneTagFilter[40];   //[nElectrons]
   Bool_t          ele_passTPTwoTagFilter[40];   //[nElectrons]
   Bool_t          ele_passTPOneProbeFilter[40];   //[nElectrons]
   Bool_t          ele_passTPTwoProbeFilter[40];   //[nElectrons]
   Bool_t          ele_passHLTFilter[40][100];   //[nElectrons]
   Int_t           nTaus;
   Float_t         tauE[40];   //[nTaus]
   Float_t         tauPt[40];   //[nTaus]
   Float_t         tauEta[40];   //[nTaus]
   Float_t         tauPhi[40];   //[nTaus]
   Bool_t          tau_IsLoose[40];   //[nTaus]
   Bool_t          tau_IsMedium[40];   //[nTaus]
   Bool_t          tau_IsTight[40];   //[nTaus]
   Bool_t          tau_passEleVetoLoose[40];   //[nTaus]
   Bool_t          tau_passEleVetoMedium[40];   //[nTaus]
   Bool_t          tau_passEleVetoTight[40];   //[nTaus]
   Bool_t          tau_passMuVetoLoose[40];   //[nTaus]
   Bool_t          tau_passMuVetoMedium[40];   //[nTaus]
   Bool_t          tau_passMuVetoTight[40];   //[nTaus]
   UInt_t          tau_ID[40];   //[nTaus]
   Float_t         tau_combinedIsoDeltaBetaCorr3Hits[40];   //[nTaus]
   Float_t         tau_chargedIsoPtSum[40];   //[nTaus]
   Float_t         tau_neutralIsoPtSum[40];   //[nTaus]
   Float_t         tau_puCorrPtSum[40];   //[nTaus]
   Float_t         tau_eleVetoMVA[40];   //[nTaus]
   Int_t           tau_eleVetoCategory[40];   //[nTaus]
   Float_t         tau_muonVetoMVA[40];   //[nTaus]
   Float_t         tau_isoMVAnewDMwLT[40];   //[nTaus]
   Float_t         tau_isoMVAnewDMwoLT[40];   //[nTaus]
   Float_t         tau_leadCandPt[40];   //[nTaus]
   Int_t           tau_leadCandID[40];   //[nTaus]
   Float_t         tau_leadChargedHadrCandPt[40];   //[nTaus]
   Int_t           tau_leadChargedHadrCandID[40];   //[nTaus]
   UInt_t          nIsoPFCandidates;
   Float_t         isoPFCandidatePt[40];   //[nIsoPFCandidates]
   Float_t         isoPFCandidateEta[40];   //[nIsoPFCandidates]
   Float_t         isoPFCandidatePhi[40];   //[nIsoPFCandidates]
   Float_t         isoPFCandidateIso04[40];   //[nIsoPFCandidates]
   Float_t         isoPFCandidateD0[40];   //[nIsoPFCandidates]
   Int_t           isoPFCandidatePdgId[40];   //[nIsoPFCandidates]
   Int_t           nPhotons;
   Float_t         phoE[40];   //[nPhotons]
   Float_t         phoPt[40];   //[nPhotons]
   Float_t         phoEta[40];   //[nPhotons]
   Float_t         phoPhi[40];   //[nPhotons]
   Float_t         phoSigmaIetaIeta[40];   //[nPhotons]
   Float_t         phoFull5x5SigmaIetaIeta[40];   //[nPhotons]
   Float_t         phoR9[40];   //[nPhotons]
   Float_t         pho_HoverE[40];   //[nPhotons]
   Float_t         pho_sumChargedHadronPtAllVertices[40][500];   //[nPhotons]
   Float_t         pho_sumChargedHadronPt[40];   //[nPhotons]
   Float_t         pho_sumNeutralHadronEt[40];   //[nPhotons]
   Float_t         pho_sumPhotonEt[40];   //[nPhotons]
   Float_t         pho_sumWorstVertexChargedHadronPt[40];   //[nPhotons]
   Float_t         pho_pfIsoChargedHadronIso[40];   //[nPhotons]
   Float_t         pho_pfIsoChargedHadronIsoWrongVtx[40];   //[nPhotons]
   Float_t         pho_pfIsoNeutralHadronIso[40];   //[nPhotons]
   Float_t         pho_pfIsoPhotonIso[40];   //[nPhotons]
   Float_t         pho_pfIsoModFrixione[40];   //[nPhotons]
   Float_t         pho_pfIsoSumPUPt[40];   //[nPhotons]
   Bool_t          pho_isConversion[40];   //[nPhotons]
   Bool_t          pho_passEleVeto[40];   //[nPhotons]
   Float_t         pho_RegressionE[40];   //[nPhotons]
   Float_t         pho_RegressionEUncertainty[40];   //[nPhotons]
   Float_t         pho_IDMVA[40];   //[nPhotons]
   Float_t         pho_superClusterEnergy[40];   //[nPhotons]
   Float_t         pho_superClusterRawEnergy[40];   //[nPhotons]
   Float_t         pho_superClusterEta[40];   //[nPhotons]
   Float_t         pho_superClusterPhi[40];   //[nPhotons]
   Float_t         pho_superClusterSeedX[40];   //[nPhotons]
   Float_t         pho_superClusterSeedY[40];   //[nPhotons]
   Float_t         pho_superClusterSeedZ[40];   //[nPhotons]
   Float_t         pho_superClusterSeedT[40];   //[nPhotons]
   Float_t         pho_superClusterX[40];   //[nPhotons]
   Float_t         pho_superClusterY[40];   //[nPhotons]
   Float_t         pho_superClusterZ[40];   //[nPhotons]
   Bool_t          pho_hasPixelSeed[40];   //[nPhotons]
   Bool_t          pho_passHLTFilter[40][100];   //[nPhotons]
   Int_t           pho_convType[40];   //[nPhotons]
   Float_t         pho_convTrkZ[40];   //[nPhotons]
   Float_t         pho_convTrkClusZ[40];   //[nPhotons]
   Float_t         pho_vtxSumPx[40][500];   //[nPhotons]
   Float_t         pho_vtxSumPy[40][500];   //[nPhotons]
   Int_t           nJets;
   Float_t         jetE[100];   //[nJets]
   Float_t         jetPt[100];   //[nJets]
   Float_t         jetEta[100];   //[nJets]
   Float_t         jetPhi[100];   //[nJets]
   Float_t         jetCSV[100];   //[nJets]
   Float_t         jetCISV[100];   //[nJets]
   Float_t         jetMass[100];   //[nJets]
   Float_t         jetJetArea[100];   //[nJets]
   Float_t         jetPileupE[100];   //[nJets]
   Float_t         jetPileupId[100];   //[nJets]
   Int_t           jetPileupIdFlag[100];   //[nJets]
   Bool_t          jetPassIDLoose[100];   //[nJets]
   Bool_t          jetPassIDTight[100];   //[nJets]
   Bool_t          jetPassMuFrac[100];   //[nJets]
   Bool_t          jetPassEleFrac[100];   //[nJets]
   Int_t           jetPartonFlavor[100];   //[nJets]
   Int_t           jetHadronFlavor[100];   //[nJets]
   Float_t         jetChargedEMEnergyFraction[100];   //[nJets]
   Float_t         jetNeutralEMEnergyFraction[100];   //[nJets]
   Float_t         jetChargedHadronEnergyFraction[100];   //[nJets]
   Float_t         jetNeutralHadronEnergyFraction[100];   //[nJets]
   Float_t         jetMuonEnergyFraction[100];   //[nJets]
   Float_t         jetHOEnergyFraction[100];   //[nJets]
   Float_t         jetHFHadronEnergyFraction[100];   //[nJets]
   Float_t         jetHFEMEnergyFraction[100];   //[nJets]
   Float_t         jetAllMuonPt[100];   //[nJets]
   Float_t         jetAllMuonEta[100];   //[nJets]
   Float_t         jetAllMuonPhi[100];   //[nJets]
   Float_t         jetAllMuonM[100];   //[nJets]
   UInt_t          nFatJets;
   Float_t         fatJetE[40];   //[nFatJets]
   Float_t         fatJetPt[40];   //[nFatJets]
   Float_t         fatJetEta[40];   //[nFatJets]
   Float_t         fatJetPhi[40];   //[nFatJets]
   Float_t         fatJetPrunedM[40];   //[nFatJets]
   Float_t         fatJetTrimmedM[40];   //[nFatJets]
   Float_t         fatJetFilteredM[40];   //[nFatJets]
   Float_t         fatJetTau1[40];   //[nFatJets]
   Float_t         fatJetTau2[40];   //[nFatJets]
   Float_t         fatJetTau3[40];   //[nFatJets]
   Float_t         metPt;
   Float_t         metPhi;
   Float_t         sumMET;
   Float_t         metType0Pt;
   Float_t         metType0Phi;
   Float_t         metType1Pt;
   Float_t         metType1Phi;
   Float_t         metType0Plus1Pt;
   Float_t         metType0Plus1Phi;
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
   //Bool_t          Flag_badChargedCandidateFilter;
   //Bool_t          Flag_badMuonFilter;
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
   Bool_t          HLTDecision[300];
   Int_t           HLTPrescale[300];
   Int_t           nGenJets;
   Float_t         genJetE[100];   //[nGenJets]
   Float_t         genJetPt[100];   //[nGenJets]
   Float_t         genJetEta[100];   //[nGenJets]
   Float_t         genJetPhi[100];   //[nGenJets]
   Float_t         genMetPt;
   Float_t         genMetPhi;
   Float_t         genVertexZ;
   Float_t         genVertexX;
   Float_t         genVertexY;
   Float_t         genVertexT;
   Float_t         genWeight;
   UInt_t          genSignalProcessID;
   Float_t         genQScale;
   Float_t         genAlphaQCD;
   Float_t         genAlphaQED;
   //std::vector<std::string> *lheComments;
   std::vector<float>   *scaleWeights;
   std::vector<float>   *pdfWeights;
   std::vector<float>   *allTrackdT;
   std::vector<float>   *allTrackPt;
   std::vector<float>   *allTrackdZ;
   std::vector<int>   *allTrackPvIndex;
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

   // List of branches
   TBranch        *b_isData;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_runNum;   //!
   TBranch        *b_lumiNum;   //!
   TBranch        *b_eventNum;   //!
   TBranch        *b_pvX;   //!
   TBranch        *b_pvY;   //!
   TBranch        *b_pvZ;   //!
   TBranch        *b_fixedGridRhoAll;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetAllCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_beamSpotX;   //!
   TBranch        *b_beamSpotY;   //!
   TBranch        *b_beamSpotZ;   //!
   TBranch        *b_nPVAll;   //!
   TBranch        *b_pvAllX;   //!
   TBranch        *b_pvIndex;   //!
   TBranch        *b_pvAllY;   //!
   TBranch        *b_pvAllZ;   //!
   TBranch        *b_pvAllT;   //!
   TBranch        *b_pvAllLogSumPtSq;   //!
   TBranch        *b_pvAllSumPx;   //!
   TBranch        *b_pvNtrack_orig;   //!
   TBranch        *b_pvAllSumPt;   //!
   TBranch        *b_pvAllSumPy;   //!
   TBranch        *b_pvAllLogSumPtSq_dt;   //!
   TBranch        *b_pvAllSumPx_dt;   //!
   TBranch        *b_pvNtrack_orig_dt;   //!
   TBranch        *b_pvAllSumPt_dt;   //!
   TBranch        *b_pvAllSumPy_dt;   //!

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
   TBranch        *b_ele_IDMVATrig;   //!
   TBranch        *b_ele_IDMVANonTrig;   //!
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
   TBranch        *b_pho_superClusterSeedX;   //!
   TBranch        *b_pho_superClusterSeedY;   //!
   TBranch        *b_pho_superClusterSeedZ;   //!
   TBranch        *b_pho_superClusterSeedT;   //!
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
   TBranch        *b_nFatJets;   //!
   TBranch        *b_fatJetE;   //!
   TBranch        *b_fatJetPt;   //!
   TBranch        *b_fatJetEta;   //!
   TBranch        *b_fatJetPhi;   //!
   TBranch        *b_fatJetPrunedM;   //!
   TBranch        *b_fatJetTrimmedM;   //!
   TBranch        *b_fatJetFilteredM;   //!
   TBranch        *b_fatJetTau1;   //!
   TBranch        *b_fatJetTau2;   //!
   TBranch        *b_fatJetTau3;   //!
   TBranch        *b_metPt;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_sumMET;   //!
   TBranch        *b_metType0Pt;   //!
   TBranch        *b_metType0Phi;   //!
   TBranch        *b_metType1Pt;   //!
   TBranch        *b_metType1Phi;   //!
   TBranch        *b_metType0Plus1Pt;   //!
   TBranch        *b_metType0Plus1Phi;   //!
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
   //TBranch        *b_Flag_badChargedCandidateFilter;   //!
   //TBranch        *b_Flag_badMuonFilter;   //!
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
   TBranch        *b_HLTDecision;   //!
   TBranch        *b_HLTPrescale;   //!
   TBranch        *b_nGenJets;   //!
   TBranch        *b_genJetE;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genMetPt;   //!
   TBranch        *b_genMetPhi;   //!
   TBranch        *b_genVertexZ;   //!
   TBranch        *b_genVertexX;   //!
   TBranch        *b_genVertexY;   //!
   TBranch        *b_genVertexT;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genSignalProcessID;   //!
   TBranch        *b_genQScale;   //!
   TBranch        *b_genAlphaQCD;   //!
   TBranch        *b_genAlphaQED;   //!
   //TBranch        *b_lheComments;   //!
   TBranch        *b_scaleWeights;   //!
   TBranch        *b_pdfWeights;   //!
   TBranch        *b_allTrackdT;   //!
   TBranch        *b_allTrackPt;   //!
   TBranch        *b_allTrackdZ;   //!
   TBranch        *b_allTrackPvIndex;   //!
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

   RazorEventsUpgradeTiming(TTree *tree=0);
   virtual ~RazorEventsUpgradeTiming();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RazorEventsUpgradeTiming_cxx
RazorEventsUpgradeTiming::RazorEventsUpgradeTiming(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/run2/Run2RazorNtupleV3.6/MC/RunIISpring16/v1/sixie/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Run2RazorNtuplerV3p6_ToCERN_MC_25ns_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1_v1_v1/160925_172232/0000/razorNtuple_756.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/run2/Run2RazorNtupleV3.6/MC/RunIISpring16/v1/sixie/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Run2RazorNtuplerV3p6_ToCERN_MC_25ns_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1_v1_v1/160925_172232/0000/razorNtuple_756.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/run2/Run2RazorNtupleV3.6/MC/RunIISpring16/v1/sixie/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Run2RazorNtuplerV3p6_ToCERN_MC_25ns_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1_v1_v1/160925_172232/0000/razorNtuple_756.root:/ntuples");
      dir->GetObject("RazorEventsUpgradeTiming",tree);

   }
   Init(tree);
}

RazorEventsUpgradeTiming::~RazorEventsUpgradeTiming()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RazorEventsUpgradeTiming::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RazorEventsUpgradeTiming::LoadTree(Long64_t entry)
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

void RazorEventsUpgradeTiming::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
  //lheComments = 0;
   scaleWeights = 0;
   pdfWeights = 0;
   allTrackdT = 0;
   allTrackPt = 0;
   allTrackdZ = 0;
   allTrackPvIndex = 0;
   alphasWeights = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("runNum", &runNum, &b_runNum);
   fChain->SetBranchAddress("lumiNum", &lumiNum, &b_lumiNum);
   fChain->SetBranchAddress("eventNum", &eventNum, &b_eventNum);
   fChain->SetBranchAddress("pvX", &pvX, &b_pvX);
   fChain->SetBranchAddress("pvY", &pvY, &b_pvY);
   fChain->SetBranchAddress("pvZ", &pvZ, &b_pvZ);
   fChain->SetBranchAddress("fixedGridRhoAll", &fixedGridRhoAll, &b_fixedGridRhoAll);
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("fixedGridRhoFastjetAllCalo", &fixedGridRhoFastjetAllCalo, &b_fixedGridRhoFastjetAllCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, &b_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
//   fChain->SetBranchAddress("beamSpotX", &beamSpotX, &b_beamSpotX);
//   fChain->SetBranchAddress("beamSpotY", &beamSpotY, &b_beamSpotY);
//   fChain->SetBranchAddress("beamSpotZ", &beamSpotZ, &b_beamSpotZ);
   fChain->SetBranchAddress("nPVAll", &nPVAll, &b_nPVAll);
   fChain->SetBranchAddress("pvIndex", pvIndex, &b_pvIndex);
   fChain->SetBranchAddress("pvAllX", pvAllX, &b_pvAllX);
   fChain->SetBranchAddress("pvAllY", pvAllY, &b_pvAllY);
   fChain->SetBranchAddress("pvAllZ", pvAllZ, &b_pvAllZ);
   fChain->SetBranchAddress("pvAllT", pvAllT, &b_pvAllT);
   fChain->SetBranchAddress("pvAllLogSumPtSq_dt", pvAllLogSumPtSq_dt, &b_pvAllLogSumPtSq_dt);
   fChain->SetBranchAddress("pvAllSumPx_dt", pvAllSumPx_dt, &b_pvAllSumPx_dt);
   fChain->SetBranchAddress("pvNtrack_dt", pvNtrack_orig_dt, &b_pvNtrack_orig_dt);
   fChain->SetBranchAddress("pvAllSumPt_dt", pvAllSumPt_dt, &b_pvAllSumPt_dt);
   fChain->SetBranchAddress("pvAllSumPy_dt", pvAllSumPy_dt, &b_pvAllSumPy_dt);
   fChain->SetBranchAddress("pvAllLogSumPtSq", pvAllLogSumPtSq, &b_pvAllLogSumPtSq);
   fChain->SetBranchAddress("pvAllSumPx", pvAllSumPx, &b_pvAllSumPx);
   fChain->SetBranchAddress("pvNtrack", pvNtrack_orig, &b_pvNtrack_orig);
   fChain->SetBranchAddress("pvAllSumPt", pvAllSumPt, &b_pvAllSumPt);
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
   fChain->SetBranchAddress("ele_IDMVATrig", ele_IDMVATrig, &b_ele_IDMVATrig);
   fChain->SetBranchAddress("ele_IDMVANonTrig", ele_IDMVANonTrig, &b_ele_IDMVANonTrig);
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
   fChain->SetBranchAddress("nIsoPFCandidates", &nIsoPFCandidates, &b_nIsoPFCandidates);
   fChain->SetBranchAddress("isoPFCandidatePt", isoPFCandidatePt, &b_isoPFCandidatePt);
   fChain->SetBranchAddress("isoPFCandidateEta", isoPFCandidateEta, &b_isoPFCandidateEta);
   fChain->SetBranchAddress("isoPFCandidatePhi", isoPFCandidatePhi, &b_isoPFCandidatePhi);
   fChain->SetBranchAddress("isoPFCandidateIso04", isoPFCandidateIso04, &b_isoPFCandidateIso04);
   fChain->SetBranchAddress("isoPFCandidateD0", isoPFCandidateD0, &b_isoPFCandidateD0);
   fChain->SetBranchAddress("isoPFCandidatePdgId", isoPFCandidatePdgId, &b_isoPFCandidatePdgId);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
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
   fChain->SetBranchAddress("pho_superClusterSeedX", pho_superClusterSeedX, &b_pho_superClusterSeedX);
   fChain->SetBranchAddress("pho_superClusterSeedY", pho_superClusterSeedY, &b_pho_superClusterSeedY);
   fChain->SetBranchAddress("pho_superClusterSeedZ", pho_superClusterSeedZ, &b_pho_superClusterSeedZ);
   fChain->SetBranchAddress("pho_superClusterSeedT", pho_superClusterSeedT, &b_pho_superClusterSeedT);
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
   fChain->SetBranchAddress("nFatJets", &nFatJets, &b_nFatJets);
   fChain->SetBranchAddress("fatJetE", fatJetE, &b_fatJetE);
   fChain->SetBranchAddress("fatJetPt", fatJetPt, &b_fatJetPt);
   fChain->SetBranchAddress("fatJetEta", fatJetEta, &b_fatJetEta);
   fChain->SetBranchAddress("fatJetPhi", fatJetPhi, &b_fatJetPhi);
   fChain->SetBranchAddress("fatJetPrunedM", fatJetPrunedM, &b_fatJetPrunedM);
   fChain->SetBranchAddress("fatJetTrimmedM", fatJetTrimmedM, &b_fatJetTrimmedM);
   fChain->SetBranchAddress("fatJetFilteredM", fatJetFilteredM, &b_fatJetFilteredM);
   fChain->SetBranchAddress("fatJetTau1", fatJetTau1, &b_fatJetTau1);
   fChain->SetBranchAddress("fatJetTau2", fatJetTau2, &b_fatJetTau2);
   fChain->SetBranchAddress("fatJetTau3", fatJetTau3, &b_fatJetTau3);
   fChain->SetBranchAddress("metPt", &metPt, &b_metPt);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("sumMET", &sumMET, &b_sumMET);
   fChain->SetBranchAddress("metType0Pt", &metType0Pt, &b_metType0Pt);
   fChain->SetBranchAddress("metType0Phi", &metType0Phi, &b_metType0Phi);
   fChain->SetBranchAddress("metType1Pt", &metType1Pt, &b_metType1Pt);
   fChain->SetBranchAddress("metType1Phi", &metType1Phi, &b_metType1Phi);
   fChain->SetBranchAddress("metType0Plus1Pt", &metType0Plus1Pt, &b_metType0Plus1Pt);
   fChain->SetBranchAddress("metType0Plus1Phi", &metType0Plus1Phi, &b_metType0Plus1Phi);
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
   //fChain->SetBranchAddress("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter, &b_Flag_badChargedCandidateFilter);
   //fChain->SetBranchAddress("Flag_badMuonFilter", &Flag_badMuonFilter, &b_Flag_badMuonFilter);
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
//   fChain->SetBranchAddress("HLTDecision", HLTDecision, &b_HLTDecision);
//   fChain->SetBranchAddress("HLTPrescale", HLTPrescale, &b_HLTPrescale);
   fChain->SetBranchAddress("nGenJets", &nGenJets, &b_nGenJets);
   fChain->SetBranchAddress("genJetE", genJetE, &b_genJetE);
   fChain->SetBranchAddress("genJetPt", genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetEta", genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPhi", genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genMetPt", &genMetPt, &b_genMetPt);
   fChain->SetBranchAddress("genMetPhi", &genMetPhi, &b_genMetPhi);
   fChain->SetBranchAddress("genVertexZ", &genVertexZ, &b_genVertexZ);
   fChain->SetBranchAddress("genVertexX", &genVertexX, &b_genVertexX);
   fChain->SetBranchAddress("genVertexY", &genVertexY, &b_genVertexY);
   fChain->SetBranchAddress("genVertexT", &genVertexT, &b_genVertexT);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genSignalProcessID", &genSignalProcessID, &b_genSignalProcessID);
   fChain->SetBranchAddress("genQScale", &genQScale, &b_genQScale);
   fChain->SetBranchAddress("genAlphaQCD", &genAlphaQCD, &b_genAlphaQCD);
   fChain->SetBranchAddress("genAlphaQED", &genAlphaQED, &b_genAlphaQED);
   //fChain->SetBranchAddress("lheComments", &lheComments, &b_lheComments);
   fChain->SetBranchAddress("scaleWeights", &scaleWeights, &b_scaleWeights);
   fChain->SetBranchAddress("pdfWeights", &pdfWeights, &b_pdfWeights);
   fChain->SetBranchAddress("allTrackdT", &allTrackdT, &b_allTrackdT);
   fChain->SetBranchAddress("allTrackPt", &allTrackPt, &b_allTrackPt);
   fChain->SetBranchAddress("allTrackdZ", &allTrackdZ, &b_allTrackdZ);
   fChain->SetBranchAddress("allTrackPvIndex", &allTrackPvIndex, &b_allTrackPvIndex);
   fChain->SetBranchAddress("alphasWeights", &alphasWeights, &b_alphasWeights);
   fChain->SetBranchAddress("nGenParticle", &nGenParticle, &b_nGenParticle);
   fChain->SetBranchAddress("gParticleMotherId", gParticleMotherId, &b_gParticleMotherId);
   fChain->SetBranchAddress("gParticleMotherIndex", gParticleMotherIndex, &b_gParticleMotherIndex);
   fChain->SetBranchAddress("gParticleId", gParticleId, &b_gParticleId);
   fChain->SetBranchAddress("gParticleStatus", gParticleStatus, &b_gParticleStatus);
   fChain->SetBranchAddress("gParticleE", gParticleE, &b_gParticleE);
   fChain->SetBranchAddress("gParticlePt", gParticlePt, &b_gParticlePt);
   fChain->SetBranchAddress("gParticleEta", gParticleEta, &b_gParticleEta);
   fChain->SetBranchAddress("gParticlePhi", gParticlePhi, &b_gParticlePhi);
   Notify();
}

Bool_t RazorEventsUpgradeTiming::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RazorEventsUpgradeTiming::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RazorEventsUpgradeTiming::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RazorEventsUpgradeTiming_cxx
