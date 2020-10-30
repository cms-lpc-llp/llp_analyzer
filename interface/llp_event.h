//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 29 16:56:19 2020 by ROOT version 6.20/07
// from TTree llp/selected AOD information for llp analyses
// found on file: displacedJetMuon_ntupler_1.root
//////////////////////////////////////////////////////////

#ifndef llp_event_h
#define llp_event_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <string>
#include <vector>

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
   Int_t           BunchXing[16];   //[nBunchXing]
   Int_t           nPU[16];   //[nBunchXing]
   Float_t         nPUmean[16];   //[nBunchXing]
   Int_t           nMuons;
   Float_t         muonE[50];   //[nMuons]
   Float_t         muonPt[50];   //[nMuons]
   Float_t         muonEta[50];   //[nMuons]
   Float_t         muonPhi[50];   //[nMuons]
   Int_t           muonCharge[50];   //[nMuons]
   Bool_t          muonIsLoose[50];   //[nMuons]
   Bool_t          muonIsMedium[50];   //[nMuons]
   Bool_t          muonIsTight[50];   //[nMuons]
   Float_t         muon_d0[50];   //[nMuons]
   Float_t         muon_dZ[50];   //[nMuons]
   Float_t         muon_ip3d[50];   //[nMuons]
   Float_t         muon_ip3dSignificance[50];   //[nMuons]
   UInt_t          muonType[50];   //[nMuons]
   UInt_t          muonQuality[50];   //[nMuons]
   Float_t         muon_pileupIso[50];   //[nMuons]
   Float_t         muon_chargedIso[50];   //[nMuons]
   Float_t         muon_photonIso[50];   //[nMuons]
   Float_t         muon_neutralHadIso[50];   //[nMuons]
   Float_t         muon_ptrel[50];   //[nMuons]
   Float_t         muon_chargedMiniIso[50];   //[nMuons]
   Float_t         muon_photonAndNeutralHadronMiniIso[50];   //[nMuons]
   Float_t         muon_chargedPileupMiniIso[50];   //[nMuons]
   Float_t         muon_activityMiniIsoAnnulus[50];   //[nMuons]
   Bool_t          muon_passSingleMuTagFilter[50];   //[nMuons]
   Bool_t          muon_passHLTFilter[50][100];   //[nMuons]
   Float_t         muon_validFractionTrackerHits[50];   //[nMuons]
   Bool_t          muon_isGlobal[50];   //[nMuons]
   Float_t         muon_normChi2[50];   //[nMuons]
   Float_t         muon_chi2LocalPosition[50];   //[nMuons]
   Float_t         muon_kinkFinder[50];   //[nMuons]
   Float_t         muon_segmentCompatability[50];   //[nMuons]
   Bool_t          muonIsICHEPMedium[50];   //[nMuons]
   Int_t           nElectrons;
   Float_t         eleE[50];   //[nElectrons]
   Float_t         elePt[50];   //[nElectrons]
   Float_t         eleEta[50];   //[nElectrons]
   Float_t         elePhi[50];   //[nElectrons]
   Float_t         eleCharge[50];   //[nElectrons]
   Float_t         eleEta_SC[50];   //[nElectrons]
   Float_t         eleSigmaIetaIeta[50];   //[nElectrons]
   Float_t         eleFull5x5SigmaIetaIeta[50];   //[nElectrons]
   Float_t         eleR9[50];   //[nElectrons]
   Float_t         ele_dEta[50];   //[nElectrons]
   Float_t         ele_dPhi[50];   //[nElectrons]
   Float_t         ele_HoverE[50];   //[nElectrons]
   Float_t         ele_d0[50];   //[nElectrons]
   Float_t         ele_dZ[50];   //[nElectrons]
   Float_t         ele_ip3d[50];   //[nElectrons]
   Float_t         ele_ip3dSignificance[50];   //[nElectrons]
   Float_t         ele_pileupIso[50];   //[nElectrons]
   Float_t         ele_chargedIso[50];   //[nElectrons]
   Float_t         ele_photonIso[50];   //[nElectrons]
   Float_t         ele_neutralHadIso[50];   //[nElectrons]
   Int_t           ele_MissHits[50];   //[nElectrons]
   Bool_t          ele_PassConvVeto[50];   //[nElectrons]
   Float_t         ele_OneOverEminusOneOverP[50];   //[nElectrons]
   Float_t         ele_IDMVAGeneralPurpose[50];   //[nElectrons]
   Int_t           ele_IDMVACategoryGeneralPurpose[50];   //[nElectrons]
   Float_t         ele_IDMVAHZZ[50];   //[nElectrons]
   Int_t           ele_IDMVACategoryHZZ[50];   //[nElectrons]
   Float_t         ele_RegressionE[50];   //[nElectrons]
   Float_t         ele_CombineP4[50];   //[nElectrons]
   Float_t         ele_ptrel[50];   //[nElectrons]
   Float_t         ele_chargedMiniIso[50];   //[nElectrons]
   Float_t         ele_photonAndNeutralHadronMiniIso[50];   //[nElectrons]
   Float_t         ele_chargedPileupMiniIso[50];   //[nElectrons]
   Float_t         ele_activityMiniIsoAnnulus[50];   //[nElectrons]
   Bool_t          ele_passSingleEleTagFilter[50];   //[nElectrons]
   Bool_t          ele_passTPOneTagFilter[50];   //[nElectrons]
   Bool_t          ele_passTPTwoTagFilter[50];   //[nElectrons]
   Bool_t          ele_passTPOneProbeFilter[50];   //[nElectrons]
   Bool_t          ele_passTPTwoProbeFilter[50];   //[nElectrons]
   Bool_t          ele_passHLTFilter[50][100];   //[nElectrons]
   Bool_t          ele_passCutBasedIDVeto[50];   //[nElectrons]
   Bool_t          ele_passCutBasedIDLoose[50];   //[nElectrons]
   Bool_t          ele_passCutBasedIDMedium[50];   //[nElectrons]
   Bool_t          ele_passCutBasedIDTight[50];   //[nElectrons]
   Bool_t          ele_passMVAIsoIDWP80[50];   //[nElectrons]
   Bool_t          ele_passMVAIsoIDWP90[50];   //[nElectrons]
   Bool_t          ele_passMVAIsoIDWP9HZZ[50];   //[nElectrons]
   Bool_t          ele_passMVAIsoIDWPLoose[50];   //[nElectrons]
   Bool_t          ele_passMVANoIsoIDWP80[50];   //[nElectrons]
   Bool_t          ele_passMVANoIsoIDWP90[50];   //[nElectrons]
   Bool_t          ele_passMVANoIsoIDWPLoose[50];   //[nElectrons]
   Int_t           nTaus;
   Float_t         tauE[50];   //[nTaus]
   Float_t         tauPt[50];   //[nTaus]
   Float_t         tauEta[50];   //[nTaus]
   Float_t         tauPhi[50];   //[nTaus]
   Bool_t          tau_IsLoose[50];   //[nTaus]
   Bool_t          tau_IsMedium[50];   //[nTaus]
   Bool_t          tau_IsTight[50];   //[nTaus]
   Bool_t          tau_passEleVetoLoose[50];   //[nTaus]
   Bool_t          tau_passEleVetoMedium[50];   //[nTaus]
   Bool_t          tau_passEleVetoTight[50];   //[nTaus]
   Bool_t          tau_passMuVetoLoose[50];   //[nTaus]
   Bool_t          tau_passMuVetoMedium[50];   //[nTaus]
   Bool_t          tau_passMuVetoTight[50];   //[nTaus]
   UInt_t          tau_ID[50];   //[nTaus]
   Float_t         tau_combinedIsoDeltaBetaCorr3Hits[50];   //[nTaus]
   Float_t         tau_chargedIsoPtSum[50];   //[nTaus]
   Float_t         tau_neutralIsoPtSum[50];   //[nTaus]
   Float_t         tau_puCorrPtSum[50];   //[nTaus]
   Float_t         tau_eleVetoMVA[50];   //[nTaus]
   Int_t           tau_eleVetoCategory[50];   //[nTaus]
   Float_t         tau_muonVetoMVA[50];   //[nTaus]
   Float_t         tau_isoMVAnewDMwLT[50];   //[nTaus]
   Float_t         tau_isoMVAnewDMwoLT[50];   //[nTaus]
   Float_t         tau_leadCandPt[50];   //[nTaus]
   Int_t           tau_leadCandID[50];   //[nTaus]
   Float_t         tau_leadChargedHadrCandPt[50];   //[nTaus]
   Int_t           tau_leadChargedHadrCandID[50];   //[nTaus]
   UInt_t          nPFCandidates;
   Int_t           PFCandidatePdgId[2000];   //[nPFCandidates]
   Float_t         PFCandidatePt[2000];   //[nPFCandidates]
   Float_t         PFCandidateEta[2000];   //[nPFCandidates]
   Float_t         PFCandidatePhi[2000];   //[nPFCandidates]
   Int_t           PFCandidateTrackIndex[2000];   //[nPFCandidates]
   Int_t           PFCandidatePVIndex[2000];   //[nPFCandidates]
   Int_t           nPhotons;
   Int_t           nPhotons_overlap;
   Float_t         phoE[50];   //[nPhotons]
   Float_t         phoPt[50];   //[nPhotons]
   Float_t         phoEta[50];   //[nPhotons]
   Float_t         phoPhi[50];   //[nPhotons]
   Float_t         phoSigmaIetaIeta[50];   //[nPhotons]
   Float_t         phoFull5x5SigmaIetaIeta[50];   //[nPhotons]
   Float_t         phoR9[50];   //[nPhotons]
   Float_t         pho_sminor[50];   //[nPhotons]
   Float_t         pho_smajor[50];   //[nPhotons]
   Float_t         pho_HoverE[50];   //[nPhotons]
   Float_t         pho_sumChargedHadronPtAllVertices[50][1000];   //[nPhotons]
   Float_t         pho_sumChargedHadronPt[50];   //[nPhotons]
   Float_t         pho_sumNeutralHadronEt[50];   //[nPhotons]
   Float_t         pho_sumPhotonEt[50];   //[nPhotons]
   Float_t         pho_ecalPFClusterIso[50];   //[nPhotons]
   Float_t         pho_hcalPFClusterIso[50];   //[nPhotons]
   Float_t         pho_trkSumPtHollowConeDR03[50];   //[nPhotons]
   Float_t         pho_sumWorstVertexChargedHadronPt[50];   //[nPhotons]
   Float_t         pho_pfIsoChargedHadronIso[50];   //[nPhotons]
   Float_t         pho_pfIsoChargedHadronIsoWrongVtx[50];   //[nPhotons]
   Float_t         pho_pfIsoNeutralHadronIso[50];   //[nPhotons]
   Float_t         pho_pfIsoPhotonIso[50];   //[nPhotons]
   Float_t         pho_pfIsoModFrixione[50];   //[nPhotons]
   Float_t         pho_pfIsoSumPUPt[50];   //[nPhotons]
   Bool_t          pho_isConversion[50];   //[nPhotons]
   Bool_t          pho_passEleVeto[50];   //[nPhotons]
   Float_t         pho_RegressionE[50];   //[nPhotons]
   Float_t         pho_RegressionEUncertainty[50];   //[nPhotons]
   Float_t         pho_IDMVA[50];   //[nPhotons]
   Float_t         pho_superClusterEnergy[50];   //[nPhotons]
   Float_t         pho_superClusterRawEnergy[50];   //[nPhotons]
   Float_t         pho_superClusterEta[50];   //[nPhotons]
   Float_t         pho_superClusterPhi[50];   //[nPhotons]
   Float_t         pho_superClusterX[50];   //[nPhotons]
   Float_t         pho_superClusterY[50];   //[nPhotons]
   Float_t         pho_superClusterZ[50];   //[nPhotons]
   Bool_t          pho_hasPixelSeed[50];   //[nPhotons]
   Bool_t          pho_passHLTFilter[50][100];   //[nPhotons]
   Bool_t          pho_passCutBasedIDLoose[50];   //[nPhotons]
   Bool_t          pho_passCutBasedIDMedium[50];   //[nPhotons]
   Bool_t          pho_passCutBasedIDTight[50];   //[nPhotons]
   Bool_t          pho_passMVAIDWP80[50];   //[nPhotons]
   Bool_t          pho_passMVAIDWP90[50];   //[nPhotons]
   Int_t           pho_convType[50];   //[nPhotons]
   Float_t         pho_convTrkZ[50];   //[nPhotons]
   Float_t         pho_convTrkClusZ[50];   //[nPhotons]
   Float_t         pho_vtxSumPx[50][1000];   //[nPhotons]
   Float_t         pho_vtxSumPy[50][1000];   //[nPhotons]
   Bool_t          pho_isStandardPhoton[50];   //[nPhotons]
   Float_t         pho_seedRecHitSwitchToGain6[50];   //[nPhotons]
   Float_t         pho_seedRecHitSwitchToGain1[50];   //[nPhotons]
   Float_t         pho_anyRecHitSwitchToGain6[50];   //[nPhotons]
   Float_t         pho_anyRecHitSwitchToGain1[50];   //[nPhotons]
   Int_t           nCscWireDigis;
   Int_t           nCscStripDigis;
   Int_t           nCscSeg;
   Float_t         cscSegPhi[100];   //[nCscSeg]
   Float_t         cscSegEta[100];   //[nCscSeg]
   Float_t         cscSegX[100];   //[nCscSeg]
   Float_t         cscSegY[100];   //[nCscSeg]
   Float_t         cscSegZ[100];   //[nCscSeg]
   Float_t         cscSegT[100];   //[nCscSeg]
   Float_t         cscSegChi2[100];   //[nCscSeg]
   Int_t           cscSegChamber[100];   //[nCscSeg]
   Int_t           cscSegStation[100];   //[nCscSeg]
   Int_t           cscSegNRecHits[100];   //[nCscSeg]
   Int_t           ncscRechits;
   Float_t         cscRechitsPhi[1000];   //[ncscRechits]
   Float_t         cscRechitsEta[1000];   //[ncscRechits]
   Float_t         cscRechitsX[1000];   //[ncscRechits]
   Float_t         cscRechitsY[1000];   //[ncscRechits]
   Float_t         cscRechitsZ[1000];   //[ncscRechits]
   Float_t         cscRechitsE[1000];   //[ncscRechits]
   Float_t         cscRechitsTpeak[1000];   //[ncscRechits]
   Float_t         cscRechitsTwire[1000];   //[ncscRechits]
   Int_t           cscRechitsQuality[1000];   //[ncscRechits]
   Int_t           cscRechitsChamber[1000];   //[ncscRechits]
   Int_t           cscRechitsStation[1000];   //[ncscRechits]
   Int_t           cscRechitsClusterId[1000];   //[ncscRechits]
   Int_t           cscRechitsChannels[1000];   //[ncscRechits]
   UInt_t          cscRechitsNStrips[1000];   //[ncscRechits]
   Int_t           cscRechitsHitWire[1000];   //[ncscRechits]
   Int_t           cscRechitsWGroupsBX[1000];   //[ncscRechits]
   UInt_t          cscRechitsNWireGroups[1000];   //[ncscRechits]
   Int_t           cscRechitsDetId[1000];   //[ncscRechits]
   Int_t           nCscRechitClusters;
   Float_t         cscRechitCluster_match_cscSegCluster_minDeltaR[50];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_cscSegCluster_index[50];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gParticle_minDeltaR[50];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_gParticle_index[50];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_gParticle_id[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterX[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterY[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterZ[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTime[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTimeTotal[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTimeSpread[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterGenMuonDeltaR[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMajorAxis[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMinorAxis[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterEtaPhiSpread[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterPhiSpread[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterEtaSpread[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterXSpread[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterXYSpread[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterYSpread[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterZSpread[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterPhi[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterEta[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoPt[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoE[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoPt[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoE[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterCaloJetVeto[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterSize[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus11[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus12[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus13[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus21[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus22[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus31[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus32[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus41[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus42[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus11[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus12[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus13[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus21[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus22[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus31[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus32[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus41[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus42[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMe11Ratio[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMe12Ratio[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNStation[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterMaxStation[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMaxStationRatio[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNChamber[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterMaxChamber[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMaxChamberRatio[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterVertexR[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterVertexZ[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterVertexDis[50];   //[nCscRechitClusters]
   Float_t         cscRechitClusterVertexChi2[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN1[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN5[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN10[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN15[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN20[50];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN[50];   //[nCscRechitClusters]
   Int_t           nCscSegClusters;
   Float_t         cscSegCluster_match_gParticle_minDeltaR[50];   //[nCscSegClusters]
   Int_t           cscSegCluster_match_gParticle_index[50];   //[nCscSegClusters]
   Int_t           cscSegCluster_match_gParticle_id[50];   //[nCscSegClusters]
   Float_t         cscSegClusterX[50];   //[nCscSegClusters]
   Float_t         cscSegClusterY[50];   //[nCscSegClusters]
   Float_t         cscSegClusterZ[50];   //[nCscSegClusters]
   Float_t         cscSegClusterTime[50];   //[nCscSegClusters]
   Float_t         cscSegClusterTimeSpread[50];   //[nCscSegClusters]
   Float_t         cscSegClusterGenMuonDeltaR[50];   //[nCscSegClusters]
   Float_t         cscSegClusterMajorAxis[50];   //[nCscSegClusters]
   Float_t         cscSegClusterMinorAxis[50];   //[nCscSegClusters]
   Float_t         cscSegClusterEtaPhiSpread[50];   //[nCscSegClusters]
   Float_t         cscSegClusterPhiSpread[50];   //[nCscSegClusters]
   Float_t         cscSegClusterEtaSpread[50];   //[nCscSegClusters]
   Float_t         cscSegClusterXSpread[50];   //[nCscSegClusters]
   Float_t         cscSegClusterYSpread[50];   //[nCscSegClusters]
   Float_t         cscSegClusterZSpread[50];   //[nCscSegClusters]
   Float_t         cscSegClusterPhi[50];   //[nCscSegClusters]
   Float_t         cscSegClusterEta[50];   //[nCscSegClusters]
   Float_t         cscSegClusterJetVetoPt[50];   //[nCscSegClusters]
   Float_t         cscSegClusterJetVetoE[50];   //[nCscSegClusters]
   Float_t         cscSegClusterMuonVetoPt[50];   //[nCscSegClusters]
   Float_t         cscSegClusterMuonVetoE[50];   //[nCscSegClusters]
   Float_t         cscSegClusterCaloJetVeto[50];   //[nCscSegClusters]
   Int_t           cscSegClusterSize[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus11[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus12[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus13[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus21[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus22[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus31[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus32[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus41[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus42[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus11[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus12[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus13[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus21[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus22[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus31[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus32[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus41[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus42[50];   //[nCscSegClusters]
   Float_t         cscSegClusterMe11Ratio[50];   //[nCscSegClusters]
   Float_t         cscSegClusterMe12Ratio[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNStation[50];   //[nCscSegClusters]
   Int_t           cscSegClusterMaxStation[50];   //[nCscSegClusters]
   Float_t         cscSegClusterMaxStationRatio[50];   //[nCscSegClusters]
   Int_t           cscSegClusterNChamber[50];   //[nCscSegClusters]
   Int_t           cscSegClusterMaxChamber[50];   //[nCscSegClusters]
   Float_t         cscSegClusterMaxChamberRatio[50];   //[nCscSegClusters]
   Float_t         cscSegClusterVertexR[50];   //[nCscSegClusters]
   Float_t         cscSegClusterVertexZ[50];   //[nCscSegClusters]
   Float_t         cscSegClusterVertexDis[50];   //[nCscSegClusters]
   Float_t         cscSegClusterVertexChi2[50];   //[nCscSegClusters]
   Int_t           cscSegClusterVertexN1[50];   //[nCscSegClusters]
   Int_t           cscSegClusterVertexN5[50];   //[nCscSegClusters]
   Int_t           cscSegClusterVertexN10[50];   //[nCscSegClusters]
   Int_t           cscSegClusterVertexN15[50];   //[nCscSegClusters]
   Int_t           cscSegClusterVertexN20[50];   //[nCscSegClusters]
   Int_t           cscSegClusterVertexN[50];   //[nCscSegClusters]
   Int_t           nDtRechits;
   Float_t         dtRechitX[2000];   //[nDtRechits]
   Float_t         dtRechitY[2000];   //[nDtRechits]
   Float_t         dtRechitZ[2000];   //[nDtRechits]
   Float_t         dtRechitEta[2000];   //[nDtRechits]
   Float_t         dtRechitPhi[2000];   //[nDtRechits]
   Float_t         dtRechitTime[2000];   //[nDtRechits]
   Int_t           dtRechitStation[2000];   //[nDtRechits]
   Int_t           dtRechitWheel[2000];   //[nDtRechits]
   Int_t           nDtRechitClusters;
   Float_t         dtRechitCluster_match_gParticle_minDeltaR[50];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_gParticle_index[50];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_gParticle_id[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterX[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterY[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterZ[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTime[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTimeSpread[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterGenMuonDeltaR[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMajorAxis[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMinorAxis[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterEtaPhiSpread[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterPhiSpread[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterEtaSpread[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterXSpread[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterYSpread[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterZSpread[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterPhi[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterEta[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoPt[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoE[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoPt[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoE[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterCaloJetVeto[50];   //[nDtRechitClusters]
   Int_t           dtRechitClusterSize[50];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNStation[50];   //[nDtRechitClusters]
   Int_t           dtRechitClusterMaxStation[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMaxStationRatio[50];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNChamber[50];   //[nDtRechitClusters]
   Int_t           dtRechitClusterMaxChamber[50];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMaxChamberRatio[50];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation1[50];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation2[50];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation3[50];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation4[50];   //[nDtRechitClusters]
   Int_t           nRpc;
   Float_t         rpcPhi[1000];   //[nRpc]
   Float_t         rpcEta[1000];   //[nRpc]
   Float_t         rpcX[1000];   //[nRpc]
   Float_t         rpcY[1000];   //[nRpc]
   Float_t         rpcZ[1000];   //[nRpc]
   Float_t         rpcT[1000];   //[nRpc]
   Int_t           rpcBx[1000];   //[nRpc]
   Float_t         rpcTError[1000];   //[nRpc]
   Int_t           nDtSeg;
   Float_t         dtSegPhi[500];   //[nDtSeg]
   Float_t         dtSegEta[500];   //[nDtSeg]
   Float_t         dtSegX[500];   //[nDtSeg]
   Float_t         dtSegY[500];   //[nDtSeg]
   Float_t         dtSegZ[500];   //[nDtSeg]
   Int_t           dtSegStation[500];   //[nDtSeg]
   Int_t           dtSegWheel[500];   //[nDtSeg]
   Float_t         dtSegTime[500];   //[nDtSeg]
   Float_t         dtSegTimeError[500];   //[nDtSeg]
   Int_t           nDtSegClusters;
   Float_t         dtSegCluster_match_gParticle_minDeltaR[50];   //[nDtSegClusters]
   Int_t           dtSegCluster_match_gParticle_index[50];   //[nDtSegClusters]
   Int_t           dtSegCluster_match_gParticle_id[50];   //[nDtSegClusters]
   Float_t         dtSegClusterX[50];   //[nDtSegClusters]
   Float_t         dtSegClusterY[50];   //[nDtSegClusters]
   Float_t         dtSegClusterZ[50];   //[nDtSegClusters]
   Float_t         dtSegClusterTime[50];   //[nDtSegClusters]
   Float_t         dtSegClusterTimeSpread[50];   //[nDtSegClusters]
   Float_t         dtSegClusterGenMuonDeltaR[50];   //[nDtSegClusters]
   Float_t         dtSegClusterMajorAxis[50];   //[nDtSegClusters]
   Float_t         dtSegClusterMinorAxis[50];   //[nDtSegClusters]
   Float_t         dtSegClusterEtaPhiSpread[50];   //[nDtSegClusters]
   Float_t         dtSegClusterPhiSpread[50];   //[nDtSegClusters]
   Float_t         dtSegClusterEtaSpread[50];   //[nDtSegClusters]
   Float_t         dtSegClusterXSpread[50];   //[nDtSegClusters]
   Float_t         dtSegClusterYSpread[50];   //[nDtSegClusters]
   Float_t         dtSegClusterZSpread[50];   //[nDtSegClusters]
   Float_t         dtSegClusterPhi[50];   //[nDtSegClusters]
   Float_t         dtSegClusterEta[50];   //[nDtSegClusters]
   Float_t         dtSegClusterJetVetoPt[50];   //[nDtSegClusters]
   Float_t         dtSegClusterJetVetoE[50];   //[nDtSegClusters]
   Float_t         dtSegClusterMuonVetoPt[50];   //[nDtSegClusters]
   Float_t         dtSegClusterMuonVetoE[50];   //[nDtSegClusters]
   Float_t         dtSegClusterCaloJetVeto[50];   //[nDtSegClusters]
   Int_t           dtSegClusterSize[50];   //[nDtSegClusters]
   Int_t           dtSegClusterNSegmentStation1[50];   //[nDtSegClusters]
   Int_t           dtSegClusterNSegmentStation2[50];   //[nDtSegClusters]
   Int_t           dtSegClusterNSegmentStation3[50];   //[nDtSegClusters]
   Int_t           dtSegClusterNSegmentStation4[50];   //[nDtSegClusters]
   Int_t           dtSegClusterNStation[50];   //[nDtSegClusters]
   Int_t           dtSegClusterMaxStation[50];   //[nDtSegClusters]
   Float_t         dtSegClusterMaxStationRatio[50];   //[nDtSegClusters]
   Int_t           dtSegClusterNChamber[50];   //[nDtSegClusters]
   Int_t           dtSegClusterMaxChamber[50];   //[nDtSegClusters]
   Float_t         dtSegClusterMaxChamberRatio[50];   //[nDtSegClusters]
   Int_t           nHORechits;
   Float_t         hoRechit_Eta[1000];   //[nHORechits]
   Float_t         hoRechit_Phi[1000];   //[nHORechits]
   Float_t         hoRechit_E[1000];   //[nHORechits]
   Float_t         hoRechit_X[1000];   //[nHORechits]
   Float_t         hoRechit_Y[1000];   //[nHORechits]
   Float_t         hoRechit_Z[1000];   //[nHORechits]
   Float_t         hoRechit_T[1000];   //[nHORechits]
   Int_t           nRechits;
   Float_t         ecalRechit_Eta[2000];   //[nRechits]
   Float_t         ecalRechit_Phi[2000];   //[nRechits]
   Float_t         ecalRechit_E[2000];   //[nRechits]
   Float_t         ecalRechit_T[2000];   //[nRechits]
   Float_t         ecalRechit_E_Error[2000];   //[nRechits]
   Float_t         ecalRechit_T_Error[2000];   //[nRechits]
   Bool_t          ecalRechit_kSaturatedflag[2000];   //[nRechits]
   Bool_t          ecalRechit_kLeadingEdgeRecoveredflag[2000];   //[nRechits]
   Bool_t          ecalRechit_kPoorRecoflag[2000];   //[nRechits]
   Bool_t          ecalRechit_kWeirdflag[2000];   //[nRechits]
   Bool_t          ecalRechit_kDiWeirdflag[2000];   //[nRechits]
   Int_t           nHBHERechits;
   Float_t         hbheRechit_Eta[2000];   //[nHBHERechits]
   Float_t         hbheRechit_Phi[2000];   //[nHBHERechits]
   Float_t         hbheRechit_E[2000];   //[nHBHERechits]
   Float_t         hbheRechit_X[2000];   //[nHBHERechits]
   Float_t         hbheRechit_Y[2000];   //[nHBHERechits]
   Float_t         hbheRechit_Z[2000];   //[nHBHERechits]
   Float_t         hbheRechit_T[2000];   //[nHBHERechits]
   Int_t           hbheRechit_iEta[2000];   //[nHBHERechits]
   Int_t           hbheRechit_iPhi[2000];   //[nHBHERechits]
   Int_t           hbheRechit_depth[2000];   //[nHBHERechits]
   Int_t           nTracks;
   Float_t         track_Pt[2000];   //[nTracks]
   Float_t         track_Eta[2000];   //[nTracks]
   Float_t         track_Phi[2000];   //[nTracks]
   Int_t           track_charge[2000];   //[nTracks]
   Int_t           track_bestVertexIndex[2000];   //[nTracks]
   Int_t           track_nMissingInnerHits[2000];   //[nTracks]
   Int_t           track_nMissingOuterHits[2000];   //[nTracks]
   Int_t           track_nPixelHits[2000];   //[nTracks]
   Int_t           track_nHits[2000];   //[nTracks]
   Float_t         track_dxyToBS[2000];   //[nTracks]
   Float_t         track_dxyErr[2000];   //[nTracks]
   Float_t         track_dzToPV[2000];   //[nTracks]
   Float_t         track_dzErr[2000];   //[nTracks]
   Float_t         track_chi2[2000];   //[nTracks]
   Int_t           track_ndof[2000];   //[nTracks]
   Int_t           nJets;
   Float_t         jetE[50];   //[nJets]
   Float_t         jetPt[50];   //[nJets]
   Float_t         jetEta[50];   //[nJets]
   Float_t         jetPhi[50];   //[nJets]
   Float_t         jetCSV[50];   //[nJets]
   Float_t         jetCISV[50];   //[nJets]
   Float_t         jetCMVA[50];   //[nJets]
   Float_t         jetProbb[50];   //[nJets]
   Float_t         jetProbc[50];   //[nJets]
   Float_t         jetProbudsg[50];   //[nJets]
   Float_t         jetProbbb[50];   //[nJets]
   Float_t         jetMass[50];   //[nJets]
   Float_t         jetJetArea[50];   //[nJets]
   Float_t         jetPileupE[50];   //[nJets]
   Float_t         jetPileupId[50];   //[nJets]
   Int_t           jetPileupIdFlag[50];   //[nJets]
   Bool_t          jetPassIDLoose[50];   //[nJets]
   Bool_t          jetPassIDTight[50];   //[nJets]
   Bool_t          jetPassMuFrac[50];   //[nJets]
   Bool_t          jetPassEleFrac[50];   //[nJets]
   Int_t           jetPartonFlavor[50];   //[nJets]
   Int_t           jetHadronFlavor[50];   //[nJets]
   Float_t         jetElectronEnergyFraction[50];   //[nJets]
   Float_t         jetPhotonEnergyFraction[50];   //[nJets]
   Float_t         jetChargedHadronEnergyFraction[50];   //[nJets]
   Float_t         jetNeutralHadronEnergyFraction[50];   //[nJets]
   Float_t         jetMuonEnergyFraction[50];   //[nJets]
   Float_t         jetHOEnergyFraction[50];   //[nJets]
   Float_t         jetHFHadronEnergyFraction[50];   //[nJets]
   Float_t         jetHFEMEnergyFraction[50];   //[nJets]
   Int_t           jetChargedHadronMultiplicity[50];   //[nJets]
   Int_t           jetNeutralHadronMultiplicity[50];   //[nJets]
   Int_t           jetPhotonMultiplicity[50];   //[nJets]
   Int_t           jetElectronMultiplicity[50];   //[nJets]
   Int_t           jetMuonMultiplicity[50];   //[nJets]
   Int_t           jetNSV[50];   //[nJets]
   Int_t           jetNSVCand[50];   //[nJets]
   Int_t           jetNVertexTracks[50];   //[nJets]
   Int_t           jetNSelectedTracks[50];   //[nJets]
   Float_t         jetDRSVJet[50];   //[nJets]
   Float_t         jetFlightDist2D[50];   //[nJets]
   Float_t         jetFlightDist2DError[50];   //[nJets]
   Float_t         jetFlightDist3D[50];   //[nJets]
   Float_t         jetFlightDist3DError[50];   //[nJets]
   Float_t         jetSV_x[50];   //[nJets]
   Float_t         jetSV_y[50];   //[nJets]
   Float_t         jetSV_z[50];   //[nJets]
   Int_t           jetSVNTracks[50];   //[nJets]
   Float_t         jetSVMass[50];   //[nJets]
   Float_t         jetAllMuonPt[50];   //[nJets]
   Float_t         jetAllMuonEta[50];   //[nJets]
   Float_t         jetAllMuonPhi[50];   //[nJets]
   Float_t         jetAllMuonM[50];   //[nJets]
   Float_t         jetPtWeightedDZ[50];   //[nJets]
   Int_t           jetNRechits[50];   //[nJets]
   Float_t         jetRechitE[50];   //[nJets]
   Float_t         jetRechitT[50];   //[nJets]
   Float_t         jetRechitT_rms[50];   //[nJets]
   Float_t         jetRechitE_Error[50];   //[nJets]
   Float_t         jetRechitT_Error[50];   //[nJets]
   Float_t         jetAlphaMax[50];   //[nJets]
   Float_t         jetBetaMax[50];   //[nJets]
   Float_t         jetGammaMax_ET[50];   //[nJets]
   Float_t         jetGammaMax_EM[50];   //[nJets]
   Float_t         jetGammaMax_Hadronic[50];   //[nJets]
   Float_t         jetGammaMax[50];   //[nJets]
   Float_t         jetPtAllTracks[50];   //[nJets]
   Float_t         jetPtAllPVTracks[50];   //[nJets]
   Float_t         jetMedianTheta2D[50];   //[nJets]
   Float_t         jetMedianIP[50];   //[nJets]
   Float_t         jetMinDeltaRAllTracks[50];   //[nJets]
   Float_t         jetMinDeltaRPVTracks[50];   //[nJets]
   Float_t         jet_sig_et1[50];   //[nJets]
   Float_t         jet_sig_et2[50];   //[nJets]
   Float_t         jet_energy_frac[50];   //[nJets]
   Float_t         jetAlphaMax_wp[50];   //[nJets]
   Float_t         jetBetaMax_wp[50];   //[nJets]
   Float_t         jetGammaMax_ET_wp[50];   //[nJets]
   Float_t         jetGammaMax_EM_wp[50];   //[nJets]
   Float_t         jetGammaMax_Hadronic_wp[50];   //[nJets]
   Float_t         jetGammaMax_wp[50];   //[nJets]
   Float_t         jetPtAllTracks_wp[50];   //[nJets]
   Float_t         jetPtAllPVTracks_wp[50];   //[nJets]
   Float_t         jetMedianTheta2D_wp[50];   //[nJets]
   Float_t         jetMedianIP_wp[50];   //[nJets]
   Float_t         jetMinDeltaRAllTracks_wp[50];   //[nJets]
   Float_t         jetMinDeltaRPVTracks_wp[50];   //[nJets]
   Int_t           jetNPFCands[50];   //[nJets]
   Int_t           jetPFCandIndex[50][5000];   //[nJets]
   Int_t          nFatJets;
   Float_t         fatJetE[50];   //[nFatJets]
   Float_t         fatJetPt[50];   //[nFatJets]
   Float_t         fatJetEta[50];   //[nFatJets]
   Float_t         fatJetPhi[50];   //[nFatJets]
   Float_t         fatJetCorrectedPt[50];   //[nFatJets]
   Float_t         fatJetPrunedM[50];   //[nFatJets]
   Float_t         fatJetTrimmedM[50];   //[nFatJets]
   Float_t         fatJetFilteredM[50];   //[nFatJets]
   Float_t         fatJetSoftDropM[50];   //[nFatJets]
   Float_t         fatJetCorrectedSoftDropM[50];   //[nFatJets]
   Float_t         fatJetUncorrectedSoftDropM[50];   //[nFatJets]
   Float_t         fatJetTau1[50];   //[nFatJets]
   Float_t         fatJetTau2[50];   //[nFatJets]
   Float_t         fatJetTau3[50];   //[nFatJets]
   Float_t         fatJetMaxSubjetCSV[50];   //[nFatJets]
   Bool_t          fatJetPassIDLoose[50];   //[nFatJets]
   Bool_t          fatJetPassIDTight[50];   //[nFatJets]
   Int_t           fatJetNPFCands[50];   //[nFatJets]
   Int_t           fatJetPFCandIndex[50][5000];   //[nFatJets]
   Float_t         metPt;
   Float_t         metPhi;
   Float_t         sumMET;
   Float_t         metUncorrectedPt;
   Float_t         metUncorrectedPhi;
   Float_t         metType1Pt;
   Float_t         metType1Phi;
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
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_globalTightHalo2016Filter;
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
   Bool_t          Flag2_globalSuperTightHalo2016Filter;
   Bool_t          Flag2_globalTightHalo2016Filter;
   Bool_t          Flag2_goodVertices;
   Bool_t          Flag2_BadChargedCandidateFilter;
   Bool_t          Flag2_BadPFMuonFilter;
   Bool_t          Flag2_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag2_HBHENoiseFilter;
   Bool_t          Flag2_HBHEIsoNoiseFilter;
   Bool_t          Flag2_ecalBadCalibFilter;
   Bool_t          Flag2_eeBadScFilter;
   Bool_t          HLTDecision[1201];
   Int_t           HLTPrescale[1201];
   Int_t           nGenJets;
   Float_t         genJetE[50];   //[nGenJets]
   Float_t         genJetPt[50];   //[nGenJets]
   Float_t         genJetEta[50];   //[nGenJets]
   Float_t         genJetPhi[50];   //[nGenJets]
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
   std::string     *lheComments;
   std::vector<float>   *scaleWeights;
   std::vector<float>   *pdfWeights;
   std::vector<float>   *alphasWeights;
   Int_t           nGenParticle;
   Int_t           gParticleMotherId[500];   //[nGenParticle]
   Int_t           gParticleMotherIndex[500];   //[nGenParticle]
   Int_t           gParticleId[500];   //[nGenParticle]
   Int_t           gParticleStatus[500];   //[nGenParticle]
   Float_t         gParticleE[500];   //[nGenParticle]
   Float_t         gParticlePt[500];   //[nGenParticle]
   Float_t         gParticleEta[500];   //[nGenParticle]
   Float_t         gParticlePhi[500];   //[nGenParticle]
   Float_t         gParticleProdVertexX[500];   //[nGenParticle]
   Float_t         gParticleProdVertexY[500];   //[nGenParticle]
   Float_t         gParticleProdVertexZ[500];   //[nGenParticle]
   Float_t         gParticleDecayVertexX[500];   //[nGenParticle]
   Float_t         gParticleDecayVertexY[500];   //[nGenParticle]
   Float_t         gParticleDecayVertexZ[500];   //[nGenParticle]
   Float_t         gLLP_decay_vertex_x[2];
   Float_t         gLLP_decay_vertex_y[2];
   Float_t         gLLP_decay_vertex_z[2];
   Float_t         gLLP_beta[2];
   Float_t         gLLP_e[2];
   Float_t         gLLP_pt[2];
   Float_t         gLLP_eta[2];
   Float_t         gLLP_phi[2];
   Bool_t          gLLP_csc[2];
   Bool_t          gLLP_dt[2];
   Float_t         gLLP_travel_time[2];
   Int_t           gLLP_daughter_id[4];
   Float_t         gLLP_daughter_pt[4];
   Float_t         gLLP_daughter_eta[4];
   Float_t         gLLP_daughter_phi[4];
   Float_t         gLLP_daughter_eta_ecalcorr[4];
   Float_t         gLLP_daughter_phi_ecalcorr[4];
   Float_t         gLLP_daughter_e[4];
   Float_t         gLLP_daughter_mass[4];
   Int_t           gLLP_grandaughter_id[4];
   Float_t         gLLP_grandaughter_pt[4];
   Float_t         gLLP_grandaughter_eta[4];
   Float_t         gLLP_grandaughter_phi[4];
   Float_t         gLLP_grandaughter_eta_ecalcorr[4];
   Float_t         gLLP_grandaughter_phi_ecalcorr[4];
   Float_t         gLLP_grandaughter_e[4];
   Float_t         gLLP_grandaughter_mass[4];

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
   TBranch        *b_ele_passCutBasedIDVeto;   //!
   TBranch        *b_ele_passCutBasedIDLoose;   //!
   TBranch        *b_ele_passCutBasedIDMedium;   //!
   TBranch        *b_ele_passCutBasedIDTight;   //!
   TBranch        *b_ele_passMVAIsoIDWP80;   //!
   TBranch        *b_ele_passMVAIsoIDWP90;   //!
   TBranch        *b_ele_passMVAIsoIDWP9HZZ;   //!
   TBranch        *b_ele_passMVAIsoIDWPLoose;   //!
   TBranch        *b_ele_passMVANoIsoIDWP80;   //!
   TBranch        *b_ele_passMVANoIsoIDWP90;   //!
   TBranch        *b_ele_passMVANoIsoIDWPLoose;   //!
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
   TBranch        *b_nPFCandidates;   //!
   TBranch        *b_PFCandidatePdgId;   //!
   TBranch        *b_PFCandidatePt;   //!
   TBranch        *b_PFCandidateEta;   //!
   TBranch        *b_PFCandidatePhi;   //!
   TBranch        *b_PFCandidateTrackIndex;   //!
   TBranch        *b_PFCandidatePVIndex;   //!
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
   TBranch        *b_pho_passCutBasedIDLoose;   //!
   TBranch        *b_pho_passCutBasedIDMedium;   //!
   TBranch        *b_pho_passCutBasedIDTight;   //!
   TBranch        *b_pho_passMVAIDWP80;   //!
   TBranch        *b_pho_passMVAIDWP90;   //!
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
   TBranch        *b_nCscWireDigis;   //!
   TBranch        *b_nCscStripDigis;   //!
   TBranch        *b_nCscSeg;   //!
   TBranch        *b_cscSegPhi;   //!
   TBranch        *b_cscSegEta;   //!
   TBranch        *b_cscSegX;   //!
   TBranch        *b_cscSegY;   //!
   TBranch        *b_cscSegZ;   //!
   TBranch        *b_cscSegT;   //!
   TBranch        *b_cscSegChi2;   //!
   TBranch        *b_cscSegChamber;   //!
   TBranch        *b_cscSegStation;   //!
   TBranch        *b_cscSegNRecHits;   //!
   TBranch        *b_ncscRechits;   //!
   TBranch        *b_cscRechitsPhi;   //!
   TBranch        *b_cscRechitsEta;   //!
   TBranch        *b_cscRechitsX;   //!
   TBranch        *b_cscRechitsY;   //!
   TBranch        *b_cscRechitsZ;   //!
   TBranch        *b_cscRechitsE;   //!
   TBranch        *b_cscRechitsTpeak;   //!
   TBranch        *b_cscRechitsTwire;   //!
   TBranch        *b_cscRechitsQuality;   //!
   TBranch        *b_cscRechitsChamber;   //!
   TBranch        *b_cscRechitsStation;   //!
   TBranch        *b_cscRechitsClusterId;   //!
   TBranch        *b_cscRechitsChannels;   //!
   TBranch        *b_cscRechitsNStrips;   //!
   TBranch        *b_cscRechitsHitWire;   //!
   TBranch        *b_cscRechitsWGroupsBX;   //!
   TBranch        *b_cscRechitsNWireGroups;   //!
   TBranch        *b_cscRechitsDetId;   //!
   TBranch        *b_nCscRechitClusters;   //!
   TBranch        *b_cscRechitCluster_match_cscSegCluster_minDeltaR;   //!
   TBranch        *b_cscRechitCluster_match_cscSegCluster_index;   //!
   TBranch        *b_cscRechitCluster_match_gParticle_minDeltaR;   //!
   TBranch        *b_cscRechitCluster_match_gParticle_index;   //!
   TBranch        *b_cscRechitCluster_match_gParticle_id;   //!
   TBranch        *b_cscRechitClusterX;   //!
   TBranch        *b_cscRechitClusterY;   //!
   TBranch        *b_cscRechitClusterZ;   //!
   TBranch        *b_cscRechitClusterTime;   //!
   TBranch        *b_cscRechitClusterTimeTotal;   //!
   TBranch        *b_cscRechitClusterTimeSpread;   //!
   TBranch        *b_cscRechitClusterGenMuonDeltaR;   //!
   TBranch        *b_cscRechitClusterMajorAxis;   //!
   TBranch        *b_cscRechitClusterMinorAxis;   //!
   TBranch        *b_cscRechitClusterEtaPhiSpread;   //!
   TBranch        *b_cscRechitClusterPhiSpread;   //!
   TBranch        *b_cscRechitClusterEtaSpread;   //!
   TBranch        *b_cscRechitClusterXSpread;   //!
   TBranch        *b_cscRechitClusterXYSpread;   //!
   TBranch        *b_cscRechitClusterYSpread;   //!
   TBranch        *b_cscRechitClusterZSpread;   //!
   TBranch        *b_cscRechitClusterPhi;   //!
   TBranch        *b_cscRechitClusterEta;   //!
   TBranch        *b_cscRechitClusterJetVetoPt;   //!
   TBranch        *b_cscRechitClusterJetVetoE;   //!
   TBranch        *b_cscRechitClusterMuonVetoPt;   //!
   TBranch        *b_cscRechitClusterMuonVetoE;   //!
   TBranch        *b_cscRechitClusterCaloJetVeto;   //!
   TBranch        *b_cscRechitClusterSize;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus11;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus12;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus13;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus21;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus22;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus31;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus32;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus41;   //!
   TBranch        *b_cscRechitClusterNRechitChamberPlus42;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus11;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus12;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus13;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus21;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus22;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus31;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus32;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus41;   //!
   TBranch        *b_cscRechitClusterNRechitChamberMinus42;   //!
   TBranch        *b_cscRechitClusterMe11Ratio;   //!
   TBranch        *b_cscRechitClusterMe12Ratio;   //!
   TBranch        *b_cscRechitClusterNStation;   //!
   TBranch        *b_cscRechitClusterMaxStation;   //!
   TBranch        *b_cscRechitClusterMaxStationRatio;   //!
   TBranch        *b_cscRechitClusterNChamber;   //!
   TBranch        *b_cscRechitClusterMaxChamber;   //!
   TBranch        *b_cscRechitClusterMaxChamberRatio;   //!
   TBranch        *b_cscRechitClusterVertexR;   //!
   TBranch        *b_cscRechitClusterVertexZ;   //!
   TBranch        *b_cscRechitClusterVertexDis;   //!
   TBranch        *b_cscRechitClusterVertexChi2;   //!
   TBranch        *b_cscRechitClusterVertexN1;   //!
   TBranch        *b_cscRechitClusterVertexN5;   //!
   TBranch        *b_cscRechitClusterVertexN10;   //!
   TBranch        *b_cscRechitClusterVertexN15;   //!
   TBranch        *b_cscRechitClusterVertexN20;   //!
   TBranch        *b_cscRechitClusterVertexN;   //!
   TBranch        *b_nCscSegClusters;   //!
   TBranch        *b_cscSegCluster_match_gParticle_minDeltaR;   //!
   TBranch        *b_cscSegCluster_match_gParticle_index;   //!
   TBranch        *b_cscSegCluster_match_gParticle_id;   //!
   TBranch        *b_cscSegClusterX;   //!
   TBranch        *b_cscSegClusterY;   //!
   TBranch        *b_cscSegClusterZ;   //!
   TBranch        *b_cscSegClusterTime;   //!
   TBranch        *b_cscSegClusterTimeSpread;   //!
   TBranch        *b_cscSegClusterGenMuonDeltaR;   //!
   TBranch        *b_cscSegClusterMajorAxis;   //!
   TBranch        *b_cscSegClusterMinorAxis;   //!
   TBranch        *b_cscSegClusterEtaPhiSpread;   //!
   TBranch        *b_cscSegClusterPhiSpread;   //!
   TBranch        *b_cscSegClusterEtaSpread;   //!
   TBranch        *b_cscSegClusterXSpread;   //!
   TBranch        *b_cscSegClusterYSpread;   //!
   TBranch        *b_cscSegClusterZSpread;   //!
   TBranch        *b_cscSegClusterPhi;   //!
   TBranch        *b_cscSegClusterEta;   //!
   TBranch        *b_cscSegClusterJetVetoPt;   //!
   TBranch        *b_cscSegClusterJetVetoE;   //!
   TBranch        *b_cscSegClusterMuonVetoPt;   //!
   TBranch        *b_cscSegClusterMuonVetoE;   //!
   TBranch        *b_cscSegClusterCaloJetVeto;   //!
   TBranch        *b_cscSegClusterSize;   //!
   TBranch        *b_cscSegClusterNSegmentChamberPlus11;   //!
   TBranch        *b_cscSegClusterNSegmentChamberPlus12;   //!
   TBranch        *b_cscSegClusterNSegmentChamberPlus13;   //!
   TBranch        *b_cscSegClusterNSegmentChamberPlus21;   //!
   TBranch        *b_cscSegClusterNSegmentChamberPlus22;   //!
   TBranch        *b_cscSegClusterNSegmentChamberPlus31;   //!
   TBranch        *b_cscSegClusterNSegmentChamberPlus32;   //!
   TBranch        *b_cscSegClusterNSegmentChamberPlus41;   //!
   TBranch        *b_cscSegClusterNSegmentChamberPlus42;   //!
   TBranch        *b_cscSegClusterNSegmentChamberMinus11;   //!
   TBranch        *b_cscSegClusterNSegmentChamberMinus12;   //!
   TBranch        *b_cscSegClusterNSegmentChamberMinus13;   //!
   TBranch        *b_cscSegClusterNSegmentChamberMinus21;   //!
   TBranch        *b_cscSegClusterNSegmentChamberMinus22;   //!
   TBranch        *b_cscSegClusterNSegmentChamberMinus31;   //!
   TBranch        *b_cscSegClusterNSegmentChamberMinus32;   //!
   TBranch        *b_cscSegClusterNSegmentChamberMinus41;   //!
   TBranch        *b_cscSegClusterNSegmentChamberMinus42;   //!
   TBranch        *b_cscSegClusterMe11Ratio;   //!
   TBranch        *b_cscSegClusterMe12Ratio;   //!
   TBranch        *b_cscSegClusterNStation;   //!
   TBranch        *b_cscSegClusterMaxStation;   //!
   TBranch        *b_cscSegClusterMaxStationRatio;   //!
   TBranch        *b_cscSegClusterNChamber;   //!
   TBranch        *b_cscSegClusterMaxChamber;   //!
   TBranch        *b_cscSegClusterMaxChamberRatio;   //!
   TBranch        *b_cscSegClusterVertexR;   //!
   TBranch        *b_cscSegClusterVertexZ;   //!
   TBranch        *b_cscSegClusterVertexDis;   //!
   TBranch        *b_cscSegClusterVertexChi2;   //!
   TBranch        *b_cscSegClusterVertexN1;   //!
   TBranch        *b_cscSegClusterVertexN5;   //!
   TBranch        *b_cscSegClusterVertexN10;   //!
   TBranch        *b_cscSegClusterVertexN15;   //!
   TBranch        *b_cscSegClusterVertexN20;   //!
   TBranch        *b_cscSegClusterVertexN;   //!
   TBranch        *b_nDtRechits;   //!
   TBranch        *b_dtRechitX;   //!
   TBranch        *b_dtRechitY;   //!
   TBranch        *b_dtRechitZ;   //!
   TBranch        *b_dtRechitEta;   //!
   TBranch        *b_dtRechitPhi;   //!
   TBranch        *b_dtRechitTime;   //!
   TBranch        *b_dtRechitStation;   //!
   TBranch        *b_dtRechitWheel;   //!
   TBranch        *b_nDtRechitClusters;   //!
   TBranch        *b_dtRechitCluster_match_gParticle_minDeltaR;   //!
   TBranch        *b_dtRechitCluster_match_gParticle_index;   //!
   TBranch        *b_dtRechitCluster_match_gParticle_id;   //!
   TBranch        *b_dtRechitClusterX;   //!
   TBranch        *b_dtRechitClusterY;   //!
   TBranch        *b_dtRechitClusterZ;   //!
   TBranch        *b_dtRechitClusterTime;   //!
   TBranch        *b_dtRechitClusterTimeSpread;   //!
   TBranch        *b_dtRechitClusterGenMuonDeltaR;   //!
   TBranch        *b_dtRechitClusterMajorAxis;   //!
   TBranch        *b_dtRechitClusterMinorAxis;   //!
   TBranch        *b_dtRechitClusterEtaPhiSpread;   //!
   TBranch        *b_dtRechitClusterPhiSpread;   //!
   TBranch        *b_dtRechitClusterEtaSpread;   //!
   TBranch        *b_dtRechitClusterXSpread;   //!
   TBranch        *b_dtRechitClusterYSpread;   //!
   TBranch        *b_dtRechitClusterZSpread;   //!
   TBranch        *b_dtRechitClusterPhi;   //!
   TBranch        *b_dtRechitClusterEta;   //!
   TBranch        *b_dtRechitClusterJetVetoPt;   //!
   TBranch        *b_dtRechitClusterJetVetoE;   //!
   TBranch        *b_dtRechitClusterMuonVetoPt;   //!
   TBranch        *b_dtRechitClusterMuonVetoE;   //!
   TBranch        *b_dtRechitClusterCaloJetVeto;   //!
   TBranch        *b_dtRechitClusterSize;   //!
   TBranch        *b_dtRechitClusterNStation;   //!
   TBranch        *b_dtRechitClusterMaxStation;   //!
   TBranch        *b_dtRechitClusterMaxStationRatio;   //!
   TBranch        *b_dtRechitClusterNChamber;   //!
   TBranch        *b_dtRechitClusterMaxChamber;   //!
   TBranch        *b_dtRechitClusterMaxChamberRatio;   //!
   TBranch        *b_dtRechitClusterNSegmentStation1;   //!
   TBranch        *b_dtRechitClusterNSegmentStation2;   //!
   TBranch        *b_dtRechitClusterNSegmentStation3;   //!
   TBranch        *b_dtRechitClusterNSegmentStation4;   //!
   TBranch        *b_nRpc;   //!
   TBranch        *b_rpcPhi;   //!
   TBranch        *b_rpcEta;   //!
   TBranch        *b_rpcX;   //!
   TBranch        *b_rpcY;   //!
   TBranch        *b_rpcZ;   //!
   TBranch        *b_rpcT;   //!
   TBranch        *b_rpcBx;   //!
   TBranch        *b_rpcTError;   //!
   TBranch        *b_nDtSeg;   //!
   TBranch        *b_dtSegPhi;   //!
   TBranch        *b_dtSegEta;   //!
   TBranch        *b_dtSegX;   //!
   TBranch        *b_dtSegY;   //!
   TBranch        *b_dtSegZ;   //!
   TBranch        *b_dtSegStation;   //!
   TBranch        *b_dtSegWheel;   //!
   TBranch        *b_dtSegTime;   //!
   TBranch        *b_dtSegTimeError;   //!
   TBranch        *b_nDtSegClusters;   //!
   TBranch        *b_dtSegCluster_match_gParticle_minDeltaR;   //!
   TBranch        *b_dtSegCluster_match_gParticle_index;   //!
   TBranch        *b_dtSegCluster_match_gParticle_id;   //!
   TBranch        *b_dtSegClusterX;   //!
   TBranch        *b_dtSegClusterY;   //!
   TBranch        *b_dtSegClusterZ;   //!
   TBranch        *b_dtSegClusterTime;   //!
   TBranch        *b_dtSegClusterTimeSpread;   //!
   TBranch        *b_dtSegClusterGenMuonDeltaR;   //!
   TBranch        *b_dtSegClusterMajorAxis;   //!
   TBranch        *b_dtSegClusterMinorAxis;   //!
   TBranch        *b_dtSegClusterEtaPhiSpread;   //!
   TBranch        *b_dtSegClusterPhiSpread;   //!
   TBranch        *b_dtSegClusterEtaSpread;   //!
   TBranch        *b_dtSegClusterXSpread;   //!
   TBranch        *b_dtSegClusterYSpread;   //!
   TBranch        *b_dtSegClusterZSpread;   //!
   TBranch        *b_dtSegClusterPhi;   //!
   TBranch        *b_dtSegClusterEta;   //!
   TBranch        *b_dtSegClusterJetVetoPt;   //!
   TBranch        *b_dtSegClusterJetVetoE;   //!
   TBranch        *b_dtSegClusterMuonVetoPt;   //!
   TBranch        *b_dtSegClusterMuonVetoE;   //!
   TBranch        *b_dtSegClusterCaloJetVeto;   //!
   TBranch        *b_dtSegClusterSize;   //!
   TBranch        *b_dtSegClusterNSegmentStation1;   //!
   TBranch        *b_dtSegClusterNSegmentStation2;   //!
   TBranch        *b_dtSegClusterNSegmentStation3;   //!
   TBranch        *b_dtSegClusterNSegmentStation4;   //!
   TBranch        *b_dtSegClusterNStation;   //!
   TBranch        *b_dtSegClusterMaxStation;   //!
   TBranch        *b_dtSegClusterMaxStationRatio;   //!
   TBranch        *b_dtSegClusterNChamber;   //!
   TBranch        *b_dtSegClusterMaxChamber;   //!
   TBranch        *b_dtSegClusterMaxChamberRatio;   //!
   TBranch        *b_nHORechits;   //!
   TBranch        *b_hoRechit_Eta;   //!
   TBranch        *b_hoRechit_Phi;   //!
   TBranch        *b_hoRechit_E;   //!
   TBranch        *b_hoRechit_X;   //!
   TBranch        *b_hoRechit_Y;   //!
   TBranch        *b_hoRechit_Z;   //!
   TBranch        *b_hoRechit_T;   //!
   TBranch        *b_nRechits;   //!
   TBranch        *b_ecalRechit_Eta;   //!
   TBranch        *b_ecalRechit_Phi;   //!
   TBranch        *b_ecalRechit_E;   //!
   TBranch        *b_ecalRechit_T;   //!
   TBranch        *b_ecalRechit_E_Error;   //!
   TBranch        *b_ecalRechit_T_Error;   //!
   TBranch        *b_ecalRechit_kSaturatedflag;   //!
   TBranch        *b_ecalRechit_kLeadingEdgeRecoveredflag;   //!
   TBranch        *b_ecalRechit_kPoorRecoflag;   //!
   TBranch        *b_ecalRechit_kWeirdflag;   //!
   TBranch        *b_ecalRechit_kDiWeirdflag;   //!
   TBranch        *b_nHBHERechits;   //!
   TBranch        *b_hbheRechit_Eta;   //!
   TBranch        *b_hbheRechit_Phi;   //!
   TBranch        *b_hbheRechit_E;   //!
   TBranch        *b_hbheRechit_X;   //!
   TBranch        *b_hbheRechit_Y;   //!
   TBranch        *b_hbheRechit_Z;   //!
   TBranch        *b_hbheRechit_T;   //!
   TBranch        *b_hbheRechit_iEta;   //!
   TBranch        *b_hbheRechit_iPhi;   //!
   TBranch        *b_hbheRechit_depth;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_track_Pt;   //!
   TBranch        *b_track_Eta;   //!
   TBranch        *b_track_Phi;   //!
   TBranch        *b_track_charge;   //!
   TBranch        *b_track_bestVertexIndex;   //!
   TBranch        *b_track_nMissingInnerHits;   //!
   TBranch        *b_track_nMissingOuterHits;   //!
   TBranch        *b_track_nPixelHits;   //!
   TBranch        *b_track_nHits;   //!
   TBranch        *b_track_dxyToBS;   //!
   TBranch        *b_track_dxyErr;   //!
   TBranch        *b_track_dzToPV;   //!
   TBranch        *b_track_dzErr;   //!
   TBranch        *b_track_chi2;   //!
   TBranch        *b_track_ndof;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetCSV;   //!
   TBranch        *b_jetCISV;   //!
   TBranch        *b_jetCMVA;   //!
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
   TBranch        *b_jetElectronEnergyFraction;   //!
   TBranch        *b_jetPhotonEnergyFraction;   //!
   TBranch        *b_jetChargedHadronEnergyFraction;   //!
   TBranch        *b_jetNeutralHadronEnergyFraction;   //!
   TBranch        *b_jetMuonEnergyFraction;   //!
   TBranch        *b_jetHOEnergyFraction;   //!
   TBranch        *b_jetHFHadronEnergyFraction;   //!
   TBranch        *b_jetHFEMEnergyFraction;   //!
   TBranch        *b_jetChargedHadronMultiplicity;   //!
   TBranch        *b_jetNeutralHadronMultiplicity;   //!
   TBranch        *b_jetPhotonMultiplicity;   //!
   TBranch        *b_jetElectronMultiplicity;   //!
   TBranch        *b_jetMuonMultiplicity;   //!
   TBranch        *b_jetNSV;   //!
   TBranch        *b_jetNSVCand;   //!
   TBranch        *b_jetNVertexTracks;   //!
   TBranch        *b_jetNSelectedTracks;   //!
   TBranch        *b_jetDRSVJet;   //!
   TBranch        *b_jetFlightDist2D;   //!
   TBranch        *b_jetFlightDist2DError;   //!
   TBranch        *b_jetFlightDist3D;   //!
   TBranch        *b_jetFlightDist3DError;   //!
   TBranch        *b_jetSV_x;   //!
   TBranch        *b_jetSV_y;   //!
   TBranch        *b_jetSV_z;   //!
   TBranch        *b_jetSVNTracks;   //!
   TBranch        *b_jetSVMass;   //!
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
   TBranch        *b_jetAlphaMax_wp;   //!
   TBranch        *b_jetBetaMax_wp;   //!
   TBranch        *b_jetGammaMax_ET_wp;   //!
   TBranch        *b_jetGammaMax_EM_wp;   //!
   TBranch        *b_jetGammaMax_Hadronic_wp;   //!
   TBranch        *b_jetGammaMax_wp;   //!
   TBranch        *b_jetPtAllTracks_wp;   //!
   TBranch        *b_jetPtAllPVTracks_wp;   //!
   TBranch        *b_jetMedianTheta2D_wp;   //!
   TBranch        *b_jetMedianIP_wp;   //!
   TBranch        *b_jetMinDeltaRAllTracks_wp;   //!
   TBranch        *b_jetMinDeltaRPVTracks_wp;   //!
   TBranch        *b_jetNPFCands;   //!
   TBranch        *b_jetPFCandIndex;   //!
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
   TBranch        *b_fatJetNPFCands;   //!
   TBranch        *b_fatJetPFCandIndex;   //!
   TBranch        *b_metPt;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_sumMET;   //!
   TBranch        *b_metUncorrectedPt;   //!
   TBranch        *b_metUncorrectedPhi;   //!
   TBranch        *b_metType1Pt;   //!
   TBranch        *b_metType1Phi;   //!
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
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
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
   TBranch        *b_Flag2_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag2_globalTightHalo2016Filter;   //!
   TBranch        *b_Flag2_goodVertices;   //!
   TBranch        *b_Flag2_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag2_BadPFMuonFilter;   //!
   TBranch        *b_Flag2_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag2_HBHENoiseFilter;   //!
   TBranch        *b_Flag2_HBHEIsoNoiseFilter;   //!
   TBranch        *b_Flag2_ecalBadCalibFilter;   //!
   TBranch        *b_Flag2_eeBadScFilter;   //!
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
   TBranch        *b_gParticleProdVertexX;   //!
   TBranch        *b_gParticleProdVertexY;   //!
   TBranch        *b_gParticleProdVertexZ;   //!
   TBranch        *b_gParticleDecayVertexX;   //!
   TBranch        *b_gParticleDecayVertexY;   //!
   TBranch        *b_gParticleDecayVertexZ;   //!
   TBranch        *b_gLLP_decay_vertex_x;   //!
   TBranch        *b_gLLP_decay_vertex_y;   //!
   TBranch        *b_gLLP_decay_vertex_z;   //!
   TBranch        *b_gLLP_beta;   //!
   TBranch        *b_gLLP_e;   //!
   TBranch        *b_gLLP_pt;   //!
   TBranch        *b_gLLP_eta;   //!
   TBranch        *b_gLLP_phi;   //!
   TBranch        *b_gLLP_csc;   //!
   TBranch        *b_gLLP_dt;   //!
   TBranch        *b_gLLP_travel_time;   //!
   TBranch        *b_gLLP_daughter_id;   //!
   TBranch        *b_gLLP_daughter_pt;   //!
   TBranch        *b_gLLP_daughter_eta;   //!
   TBranch        *b_gLLP_daughter_phi;   //!
   TBranch        *b_gLLP_daughter_eta_ecalcorr;   //!
   TBranch        *b_gLLP_daughter_phi_ecalcorr;   //!
   TBranch        *b_gLLP_daughter_e;   //!
   TBranch        *b_gLLP_daughter_mass;   //!
   TBranch        *b_gLLP_grandaughter_id;   //!
   TBranch        *b_gLLP_grandaughter_pt;   //!
   TBranch        *b_gLLP_grandaughter_eta;   //!
   TBranch        *b_gLLP_grandaughter_phi;   //!
   TBranch        *b_gLLP_grandaughter_eta_ecalcorr;   //!
   TBranch        *b_gLLP_grandaughter_phi_ecalcorr;   //!
   TBranch        *b_gLLP_grandaughter_e;   //!
   TBranch        *b_gLLP_grandaughter_mass;   //!

   llp_event(TTree *tree=0);
   virtual ~llp_event();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

