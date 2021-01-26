//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 23 13:32:01 2020 by ROOT version 6.10/09
// from TTree llp/selected AOD information for llp analyses
// found on file: /mnt/hadoop/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p16/MC_Fall18/v1/sixie/WminusH_HToSSTobbbb_ms55_pl10000_ev11000/Run2_Run2_displacedJetMuonNtupler_V1p16_MC_Fall18_apresyan-crab_PrivateProduction_Fall18_DR_step2_WminusH_HToSSTobbbb_ms55_pl10000_batch1_v1_v1_v1/200324_213938/0000/displacedJetMuon_ntupler_1.root
//////////////////////////////////////////////////////////

#ifndef llp_event_h
#define llp_event_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
// Header file for the classes stored in the TTree if any.
#include "string"
#include "vector"
using namespace std;
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
   Float_t         pvAllX[85];   //[nPVAll]
   Float_t         pvAllY[85];   //[nPVAll]
   Float_t         pvAllZ[85];   //[nPVAll]
   Float_t         pvAllLogSumPtSq[85];   //[nPVAll]
   Float_t         pvAllSumPx[85];   //[nPVAll]
   Float_t         pvAllSumPy[85];   //[nPVAll]
   Int_t           nBunchXing;
   Int_t           BunchXing[16];   //[nBunchXing]
   Int_t           nPU[16];   //[nBunchXing]
   Float_t         nPUmean[16];   //[nBunchXing]
   Int_t           nMuons;
   Float_t         muonE[200];   //[nMuons]
   Float_t         muonPt[200];   //[nMuons]
   Float_t         muonEta[200];   //[nMuons]
   Float_t         muonPhi[200];   //[nMuons]
   Int_t           muonCharge[200];   //[nMuons]
   Bool_t          muonIsLoose[200];   //[nMuons]
   Bool_t          muonIsMedium[200];   //[nMuons]
   Bool_t          muonIsTight[200];   //[nMuons]
   Float_t         muon_d0[200];   //[nMuons]
   Float_t         muon_dZ[200];   //[nMuons]
   Float_t         muon_ip3d[200];   //[nMuons]
   Float_t         muon_ip3dSignificance[200];   //[nMuons]
   UInt_t          muonType[200];   //[nMuons]
   UInt_t          muonQuality[200];   //[nMuons]
   Float_t         muon_pileupIso[200];   //[nMuons]
   Float_t         muon_chargedIso[200];   //[nMuons]
   Float_t         muon_photonIso[200];   //[nMuons]
   Float_t         muon_neutralHadIso[200];   //[nMuons]
   Float_t         muon_ptrel[200];   //[nMuons]
   Float_t         muon_chargedMiniIso[200];   //[nMuons]
   Float_t         muon_photonAndNeutralHadronMiniIso[200];   //[nMuons]
   Float_t         muon_chargedPileupMiniIso[200];   //[nMuons]
   Float_t         muon_activityMiniIsoAnnulus[200];   //[nMuons]
   Bool_t          muon_passSingleMuTagFilter[200];   //[nMuons]
   Bool_t          muon_passHLTFilter[200][100];   //[nMuons]
   Float_t         muon_validFractionTrackerHits[200];   //[nMuons]
   Bool_t          muon_isGlobal[200];   //[nMuons]
   Float_t         muon_normChi2[200];   //[nMuons]
   Float_t         muon_chi2LocalPosition[200];   //[nMuons]
   Float_t         muon_kinkFinder[200];   //[nMuons]
   Float_t         muon_segmentCompatability[200];   //[nMuons]
   Bool_t          muonIsICHEPMedium[200];   //[nMuons]
   Int_t           nElectrons;
   Float_t         eleE[20];   //[nElectrons]
   Float_t         elePt[20];   //[nElectrons]
   Float_t         eleEta[20];   //[nElectrons]
   Float_t         elePhi[20];   //[nElectrons]
   Float_t         eleCharge[20];   //[nElectrons]
   Float_t         eleEta_SC[20];   //[nElectrons]
   Float_t         eleE_SC[20];   //[nElectrons]
   Float_t         eleSigmaIetaIeta[20];   //[nElectrons]
   Float_t         eleFull5x5SigmaIetaIeta[20];   //[nElectrons]
   Float_t         eleR9[20];   //[nElectrons]
   Float_t         ele_dEta[20];   //[nElectrons]
   Float_t         ele_dPhi[20];   //[nElectrons]
   Float_t         ele_HoverE[20];   //[nElectrons]
   Float_t         ele_d0[20];   //[nElectrons]
   Float_t         ele_dZ[20];   //[nElectrons]
   Bool_t          ele_passCutBasedIDLoose[20];   //[nElectrons]
   Bool_t          ele_passCutBasedIDTight[20];   //[nElectrons]
   Float_t         ele_ip3d[20];   //[nElectrons]
   Float_t         ele_ip3dSignificance[20];   //[nElectrons]
   Float_t         ele_pileupIso[20];   //[nElectrons]
   Float_t         ele_chargedIso[20];   //[nElectrons]
   Float_t         ele_photonIso[20];   //[nElectrons]
   Float_t         ele_neutralHadIso[20];   //[nElectrons]
   Int_t           ele_MissHits[20];   //[nElectrons]
   Bool_t          ele_PassConvVeto[20];   //[nElectrons]
   Float_t         ele_OneOverEminusOneOverP[20];   //[nElectrons]
   Float_t         ele_IDMVAGeneralPurpose[20];   //[nElectrons]
   Int_t           ele_IDMVACategoryGeneralPurpose[20];   //[nElectrons]
   Float_t         ele_IDMVAHZZ[20];   //[nElectrons]
   Int_t           ele_IDMVACategoryHZZ[20];   //[nElectrons]
   Float_t         ele_RegressionE[20];   //[nElectrons]
   Float_t         ele_CombineP4[20];   //[nElectrons]
   Float_t         ele_ptrel[20];   //[nElectrons]
   Float_t         ele_chargedMiniIso[20];   //[nElectrons]
   Float_t         ele_photonAndNeutralHadronMiniIso[20];   //[nElectrons]
   Float_t         ele_chargedPileupMiniIso[20];   //[nElectrons]
   Float_t         ele_activityMiniIsoAnnulus[20];   //[nElectrons]
   Bool_t          ele_passSingleEleTagFilter[20];   //[nElectrons]
   Bool_t          ele_passTPOneTagFilter[20];   //[nElectrons]
   Bool_t          ele_passTPTwoTagFilter[20];   //[nElectrons]
   Bool_t          ele_passTPOneProbeFilter[20];   //[nElectrons]
   Bool_t          ele_passTPTwoProbeFilter[20];   //[nElectrons]
   Bool_t          ele_passHLTFilter[20][100];   //[nElectrons]
   Int_t           nTaus;
   Float_t         tauE[1];   //[nTaus]
   Float_t         tauPt[1];   //[nTaus]
   Float_t         tauEta[1];   //[nTaus]
   Float_t         tauPhi[1];   //[nTaus]
   Bool_t          tau_IsLoose[1];   //[nTaus]
   Bool_t          tau_IsMedium[1];   //[nTaus]
   Bool_t          tau_IsTight[1];   //[nTaus]
   Bool_t          tau_passEleVetoLoose[1];   //[nTaus]
   Bool_t          tau_passEleVetoMedium[1];   //[nTaus]
   Bool_t          tau_passEleVetoTight[1];   //[nTaus]
   Bool_t          tau_passMuVetoLoose[1];   //[nTaus]
   Bool_t          tau_passMuVetoMedium[1];   //[nTaus]
   Bool_t          tau_passMuVetoTight[1];   //[nTaus]
   UInt_t          tau_ID[1];   //[nTaus]
   Float_t         tau_combinedIsoDeltaBetaCorr3Hits[1];   //[nTaus]
   Float_t         tau_chargedIsoPtSum[1];   //[nTaus]
   Float_t         tau_neutralIsoPtSum[1];   //[nTaus]
   Float_t         tau_puCorrPtSum[1];   //[nTaus]
   Float_t         tau_eleVetoMVA[1];   //[nTaus]
   Int_t           tau_eleVetoCategory[1];   //[nTaus]
   Float_t         tau_muonVetoMVA[1];   //[nTaus]
   Float_t         tau_isoMVAnewDMwLT[1];   //[nTaus]
   Float_t         tau_isoMVAnewDMwoLT[1];   //[nTaus]
   Float_t         tau_leadCandPt[1];   //[nTaus]
   Int_t           tau_leadCandID[1];   //[nTaus]
   Float_t         tau_leadChargedHadrCandPt[1];   //[nTaus]
   Int_t           tau_leadChargedHadrCandID[1];   //[nTaus]
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
   Int_t           nCscWireDigis;
   Int_t           nCscStripDigis;
   Int_t           nCscSeg;
   Float_t         cscSegPhi[1000];   //[nCscSeg]
   Float_t         cscSegEta[1000];   //[nCscSeg]
   Float_t         cscSegX[1000];   //[nCscSeg]
   Float_t         cscSegY[1000];   //[nCscSeg]
   Float_t         cscSegZ[1000];   //[nCscSeg]
   Float_t         cscSegT[1000];   //[nCscSeg]
   Float_t         cscSegChi2[1000];   //[nCscSeg]
   Int_t           cscSegChamber[1000];   //[nCscSeg]
   Int_t           cscSegStation[1000];   //[nCscSeg]
   Int_t           cscSegNRecHits[1000];   //[nCscSeg]
   Int_t           ncscRechits;
   Float_t         cscRechitsPhi[50000];   //[ncscRechits]
   Float_t         cscRechitsEta[50000];   //[ncscRechits]
   Float_t         cscRechitsX[50000];   //[ncscRechits]
   Float_t         cscRechitsY[50000];   //[ncscRechits]
   Float_t         cscRechitsZ[50000];   //[ncscRechits]
   Float_t         cscRechitsE[50000];   //[ncscRechits]
   Float_t         cscRechitsTpeak[50000];   //[ncscRechits]
   Float_t         cscRechitsTwire[50000];   //[ncscRechits]
   Int_t           cscRechitsQuality[50000];   //[ncscRechits]
   Int_t           cscRechitsChamber[50000];   //[ncscRechits]
   Int_t           cscRechitsStation[50000];   //[ncscRechits]
   Int_t           cscRechitsChannels[50000];   //[ncscRechits]
   UInt_t          cscRechitsNStrips[50000];   //[ncscRechits]
   Int_t           cscRechitsHitWire[50000];   //[ncscRechits]
   Int_t           cscRechitsWGroupsBX[50000];   //[ncscRechits]
   UInt_t          cscRechitsNWireGroups[50000];   //[ncscRechits]
   Int_t           cscRechitsDetId[50000];   //[ncscRechits]
   Int_t           nCscRechitClusters;
   Float_t         cscRechitCluster_match_cscSegCluster_minDeltaR[20];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_cscSegCluster_index[20];   //[nCscRechitClusters]
   Float_t         cscRechitCluster_match_gParticle_minDeltaR[20];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_gParticle_index[20];   //[nCscRechitClusters]
   Int_t           cscRechitCluster_match_gParticle_id[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterX[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterY[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterZ[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTime[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTimeSpread[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterGenMuonDeltaR[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMajorAxis[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMinorAxis[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterEtaPhiSpread[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterPhiSpread[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterEtaSpread[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterXSpread[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterYSpread[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterZSpread[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterPhi[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterEta[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoPt[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoE[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoPt[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoE[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterCaloJetVeto[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterSize[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus11[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus12[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus13[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus21[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus22[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus31[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus32[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus41[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus42[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus11[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus12[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus13[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus21[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus22[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus31[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus32[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus41[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus42[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMe11Ratio[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMe12Ratio[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNStation[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterMaxStation[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMaxStationRatio[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNChamber[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterMaxChamber[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMaxChamberRatio[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterVertexR[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterVertexZ[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterVertexDis[20];   //[nCscRechitClusters]
   Float_t         cscRechitClusterVertexChi2[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN1[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN5[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN10[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN15[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN20[20];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN[20];   //[nCscRechitClusters]
   Int_t           nCscSegClusters;
   Float_t         cscSegCluster_match_gParticle_minDeltaR[10];   //[nCscSegClusters]
   Int_t           cscSegCluster_match_gParticle_index[10];   //[nCscSegClusters]
   Int_t           cscSegCluster_match_gParticle_id[10];   //[nCscSegClusters]
   Float_t         cscSegClusterX[10];   //[nCscSegClusters]
   Float_t         cscSegClusterY[10];   //[nCscSegClusters]
   Float_t         cscSegClusterZ[10];   //[nCscSegClusters]
   Float_t         cscSegClusterTime[10];   //[nCscSegClusters]
   Float_t         cscSegClusterTimeSpread[10];   //[nCscSegClusters]
   Float_t         cscSegClusterGenMuonDeltaR[10];   //[nCscSegClusters]
   Float_t         cscSegClusterMajorAxis[10];   //[nCscSegClusters]
   Float_t         cscSegClusterMinorAxis[10];   //[nCscSegClusters]
   Float_t         cscSegClusterEtaPhiSpread[10];   //[nCscSegClusters]
   Float_t         cscSegClusterPhiSpread[10];   //[nCscSegClusters]
   Float_t         cscSegClusterEtaSpread[10];   //[nCscSegClusters]
   Float_t         cscSegClusterXSpread[10];   //[nCscSegClusters]
   Float_t         cscSegClusterYSpread[10];   //[nCscSegClusters]
   Float_t         cscSegClusterZSpread[10];   //[nCscSegClusters]
   Float_t         cscSegClusterPhi[10];   //[nCscSegClusters]
   Float_t         cscSegClusterEta[10];   //[nCscSegClusters]
   Float_t         cscSegClusterJetVetoPt[10];   //[nCscSegClusters]
   Float_t         cscSegClusterJetVetoE[10];   //[nCscSegClusters]
   Float_t         cscSegClusterMuonVetoPt[10];   //[nCscSegClusters]
   Float_t         cscSegClusterMuonVetoE[10];   //[nCscSegClusters]
   Float_t         cscSegClusterCaloJetVeto[10];   //[nCscSegClusters]
   Int_t           cscSegClusterSize[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus11[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus12[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus13[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus21[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus22[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus31[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus32[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus41[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus42[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus11[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus12[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus13[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus21[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus22[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus31[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus32[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus41[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus42[10];   //[nCscSegClusters]
   Float_t         cscSegClusterMe11Ratio[10];   //[nCscSegClusters]
   Float_t         cscSegClusterMe12Ratio[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNStation[10];   //[nCscSegClusters]
   Int_t           cscSegClusterMaxStation[10];   //[nCscSegClusters]
   Float_t         cscSegClusterMaxStationRatio[10];   //[nCscSegClusters]
   Int_t           cscSegClusterNChamber[10];   //[nCscSegClusters]
   Int_t           cscSegClusterMaxChamber[10];   //[nCscSegClusters]
   Float_t         cscSegClusterMaxChamberRatio[10];   //[nCscSegClusters]
   Float_t         cscSegClusterVertexR[10];   //[nCscSegClusters]
   Float_t         cscSegClusterVertexZ[10];   //[nCscSegClusters]
   Float_t         cscSegClusterVertexDis[10];   //[nCscSegClusters]
   Float_t         cscSegClusterVertexChi2[10];   //[nCscSegClusters]
   Int_t           cscSegClusterVertexN1[10];   //[nCscSegClusters]
   Int_t           cscSegClusterVertexN5[10];   //[nCscSegClusters]
   Int_t           cscSegClusterVertexN10[10];   //[nCscSegClusters]
   Int_t           cscSegClusterVertexN15[10];   //[nCscSegClusters]
   Int_t           cscSegClusterVertexN20[10];   //[nCscSegClusters]
   Int_t           cscSegClusterVertexN[10];   //[nCscSegClusters]
   Int_t           nDtRechits;
   Float_t         dtRechitX[50000];   //[nDtRechits]
   Float_t         dtRechitY[50000];   //[nDtRechits]
   Float_t         dtRechitZ[50000];   //[nDtRechits]
   Float_t         dtRechitEta[50000];   //[nDtRechits]
   Float_t         dtRechitPhi[50000];   //[nDtRechits]
   Float_t         dtRechitTime[50000];   //[nDtRechits]
   Int_t           dtRechitStation[50000];   //[nDtRechits]
   Int_t           dtRechitWheel[50000];   //[nDtRechits]
   Int_t           nDtRechitClusters;
   Float_t         dtRechitCluster_match_gParticle_minDeltaR[20];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_gParticle_index[20];   //[nDtRechitClusters]
   Int_t           dtRechitCluster_match_gParticle_id[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterX[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterY[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterZ[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTime[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTimeSpread[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterGenMuonDeltaR[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMajorAxis[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMinorAxis[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterEtaPhiSpread[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterPhiSpread[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterEtaSpread[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterXSpread[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterYSpread[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterZSpread[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterPhi[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterEta[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoPt[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoE[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoPt[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoE[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterCaloJetVeto[20];   //[nDtRechitClusters]
   Int_t           dtRechitClusterSize[20];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNStation[20];   //[nDtRechitClusters]
   Int_t           dtRechitClusterMaxStation[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMaxStationRatio[20];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNChamber[20];   //[nDtRechitClusters]
   Int_t           dtRechitClusterMaxChamber[20];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMaxChamberRatio[20];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation1[20];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation2[20];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation3[20];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation4[20];   //[nDtRechitClusters]
   Int_t           nRpc;
   Float_t         rpcPhi[500];   //[nRpc]
   Float_t         rpcEta[500];   //[nRpc]
   Float_t         rpcX[500];   //[nRpc]
   Float_t         rpcY[500];   //[nRpc]
   Float_t         rpcZ[500];   //[nRpc]
   Float_t         rpcT[500];   //[nRpc]
   Int_t           rpcBx[500];   //[nRpc]
   Float_t         rpcTError[500];   //[nRpc]
   Int_t           nDtSeg;
   Float_t         dtSegPhi[1000];   //[nDtSeg]
   Float_t         dtSegEta[1000];   //[nDtSeg]
   Float_t         dtSegX[1000];   //[nDtSeg]
   Float_t         dtSegY[1000];   //[nDtSeg]
   Float_t         dtSegZ[1000];   //[nDtSeg]
   Int_t           dtSegStation[1000];   //[nDtSeg]
   Int_t           dtSegWheel[1000];   //[nDtSeg]
   Float_t         dtSegTime[1000];   //[nDtSeg]
   Float_t         dtSegTimeError[1000];   //[nDtSeg]
   Int_t           nDtSegClusters;
   Float_t         dtSegCluster_match_gParticle_minDeltaR[10];   //[nDtSegClusters]
   Int_t           dtSegCluster_match_gParticle_index[10];   //[nDtSegClusters]
   Int_t           dtSegCluster_match_gParticle_id[10];   //[nDtSegClusters]
   Float_t         dtSegClusterX[10];   //[nDtSegClusters]
   Float_t         dtSegClusterY[10];   //[nDtSegClusters]
   Float_t         dtSegClusterZ[10];   //[nDtSegClusters]
   Float_t         dtSegClusterTime[10];   //[nDtSegClusters]
   Float_t         dtSegClusterTimeSpread[10];   //[nDtSegClusters]
   Float_t         dtSegClusterGenMuonDeltaR[10];   //[nDtSegClusters]
   Float_t         dtSegClusterMajorAxis[10];   //[nDtSegClusters]
   Float_t         dtSegClusterMinorAxis[10];   //[nDtSegClusters]
   Float_t         dtSegClusterEtaPhiSpread[10];   //[nDtSegClusters]
   Float_t         dtSegClusterPhiSpread[10];   //[nDtSegClusters]
   Float_t         dtSegClusterEtaSpread[10];   //[nDtSegClusters]
   Float_t         dtSegClusterXSpread[10];   //[nDtSegClusters]
   Float_t         dtSegClusterYSpread[10];   //[nDtSegClusters]
   Float_t         dtSegClusterZSpread[10];   //[nDtSegClusters]
   Float_t         dtSegClusterPhi[10];   //[nDtSegClusters]
   Float_t         dtSegClusterEta[10];   //[nDtSegClusters]
   Float_t         dtSegClusterJetVetoPt[10];   //[nDtSegClusters]
   Float_t         dtSegClusterJetVetoE[10];   //[nDtSegClusters]
   Float_t         dtSegClusterMuonVetoPt[10];   //[nDtSegClusters]
   Float_t         dtSegClusterMuonVetoE[10];   //[nDtSegClusters]
   Float_t         dtSegClusterCaloJetVeto[10];   //[nDtSegClusters]
   Int_t           dtSegClusterSize[10];   //[nDtSegClusters]
   Int_t           dtSegClusterNSegmentStation1[10];   //[nDtSegClusters]
   Int_t           dtSegClusterNSegmentStation2[10];   //[nDtSegClusters]
   Int_t           dtSegClusterNSegmentStation3[10];   //[nDtSegClusters]
   Int_t           dtSegClusterNSegmentStation4[10];   //[nDtSegClusters]
   Int_t           dtSegClusterNStation[10];   //[nDtSegClusters]
   Int_t           dtSegClusterMaxStation[10];   //[nDtSegClusters]
   Float_t         dtSegClusterMaxStationRatio[10];   //[nDtSegClusters]
   Int_t           dtSegClusterNChamber[10];   //[nDtSegClusters]
   Int_t           dtSegClusterMaxChamber[10];   //[nDtSegClusters]
   Float_t         dtSegClusterMaxChamberRatio[10];   //[nDtSegClusters]
   Int_t           nHORechits;
   Float_t         hoRechit_Eta[500];   //[nHORechits]
   Float_t         hoRechit_Phi[500];   //[nHORechits]
   Float_t         hoRechit_E[500];   //[nHORechits]
   Float_t         hoRechit_X[500];   //[nHORechits]
   Float_t         hoRechit_Y[500];   //[nHORechits]
   Float_t         hoRechit_Z[500];   //[nHORechits]
   Float_t         hoRechit_T[500];   //[nHORechits]
   Int_t           nRechits;
   Float_t         ecalRechit_Eta[687];   //[nRechits]
   Float_t         ecalRechit_Phi[687];   //[nRechits]
   Float_t         ecalRechit_E[687];   //[nRechits]
   Float_t         ecalRechit_T[687];   //[nRechits]
   Float_t         ecalRechit_E_Error[687];   //[nRechits]
   Float_t         ecalRechit_T_Error[687];   //[nRechits]
   Bool_t          ecalRechit_kSaturatedflag[687];   //[nRechits]
   Bool_t          ecalRechit_kLeadingEdgeRecoveredflag[687];   //[nRechits]
   Bool_t          ecalRechit_kPoorRecoflag[687];   //[nRechits]
   Bool_t          ecalRechit_kWeirdflag[687];   //[nRechits]
   Bool_t          ecalRechit_kDiWeirdflag[687];   //[nRechits]
   Int_t           nHBHERechits;
   Float_t         hbheRechit_Eta[417];   //[nHBHERechits]
   Float_t         hbheRechit_Phi[417];   //[nHBHERechits]
   Float_t         hbheRechit_E[417];   //[nHBHERechits]
   Float_t         hbheRechit_X[417];   //[nHBHERechits]
   Float_t         hbheRechit_Y[417];   //[nHBHERechits]
   Float_t         hbheRechit_Z[417];   //[nHBHERechits]
   Float_t         hbheRechit_T[417];   //[nHBHERechits]
   Int_t           hbheRechit_iEta[417];   //[nHBHERechits]
   Int_t           hbheRechit_iPhi[417];   //[nHBHERechits]
   Int_t           hbheRechit_depth[417];   //[nHBHERechits]
   Int_t           nTracks;
   Float_t         track_Pt[2002];   //[nTracks]
   Float_t         track_Eta[2002];   //[nTracks]
   Float_t         track_Phi[2002];   //[nTracks]
   Int_t           track_bestVertexIndex[2002];   //[nTracks]
   Int_t           track_nMissingInnerHits[2002];   //[nTracks]
   Int_t           track_nMissingOuterHits[2002];   //[nTracks]
   Float_t         track_dxyToBS[2002];   //[nTracks]
   Float_t         track_dxyErr[2002];   //[nTracks]
   Float_t         track_dzToPV[2002];   //[nTracks]
   Float_t         track_dzErr[2002];   //[nTracks]
   Int_t           nJets;
   Float_t         jetE[100];   //[nJets]
   Float_t         jetPt[100];   //[nJets]
   Float_t         jetEta[100];   //[nJets]
   Float_t         jetPhi[100];   //[nJets]
   Float_t         jetCSV[100];   //[nJets]
   Float_t         jetCISV[100];   //[nJets]
   Float_t         jetProbb[100];   //[nJets]
   Float_t         jetProbc[100];   //[nJets]
   Float_t         jetProbudsg[100];   //[nJets]
   Float_t         jetProbbb[100];   //[nJets]
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
   Float_t          jetElectronEnergyFraction[100];   //[nJets]
  Float_t           jetPhotonEnergyFraction[100];   //[nJets]
  Int_t           jetChargedHadronMultiplicity[100];   //[nJets]
  Int_t           jetNeutralHadronMultiplicity[100];   //[nJets]
  Int_t           jetPhotonMultiplicity[100];   //[nJets]
  Int_t           jetElectronMultiplicity[100];   //[nJets]
  Int_t           jetMuonMultiplicity[100];   //[nJets]

   Float_t         jetAllMuonPt[100];   //[nJets]
   Float_t         jetAllMuonEta[100];   //[nJets]
   Float_t         jetAllMuonPhi[100];   //[nJets]
   Float_t         jetAllMuonM[100];   //[nJets]
   Float_t         jetPtWeightedDZ[100];   //[nJets]
   Int_t           jetNRechits[100];   //[nJets]
   Float_t         jetRechitE[100];   //[nJets]
   Float_t         jetRechitT[100];   //[nJets]
   Float_t         jetRechitT_rms[100];   //[nJets]
   Float_t         jetRechitE_Error[100];   //[nJets]
   Float_t         jetRechitT_Error[100];   //[nJets]
   Float_t         jetAlphaMax[100];   //[nJets]
   Float_t         jetBetaMax[100];   //[nJets]
   Float_t         jetGammaMax_ET[100];   //[nJets]
   Float_t         jetGammaMax_EM[100];   //[nJets]
   Float_t         jetGammaMax_Hadronic[100];   //[nJets]
   Float_t         jetGammaMax[100];   //[nJets]
   Float_t         jetPtAllTracks[100];   //[nJets]
   Float_t         jetPtAllPVTracks[100];   //[nJets]
   Float_t         jetMedianTheta2D[100];   //[nJets]
   Float_t         jetMedianIP[100];   //[nJets]
   Float_t         jetMinDeltaRAllTracks[100];   //[nJets]
   Float_t         jetMinDeltaRPVTracks[100];   //[nJets]
   Float_t         jet_sig_et1[100];   //[nJets]
   Float_t         jet_sig_et2[100];   //[nJets]
   Float_t         jet_energy_frac[100];   //[nJets]
   Float_t         jetAlphaMax_wp[100];   //[nJets]
   Float_t         jetBetaMax_wp[100];   //[nJets]
   Float_t         jetGammaMax_ET_wp[100];   //[nJets]
   Float_t         jetGammaMax_EM_wp[100];   //[nJets]
   Float_t         jetGammaMax_Hadronic_wp[100];   //[nJets]
   Float_t         jetGammaMax_wp[100];   //[nJets]
   Float_t         jetPtAllTracks_wp[100];   //[nJets]
   Float_t         jetPtAllPVTracks_wp[100];   //[nJets]
   Float_t         jetMedianTheta2D_wp[100];   //[nJets]
   Float_t         jetMedianIP_wp[100];   //[nJets]
   Float_t         jetMinDeltaRAllTracks_wp[100];   //[nJets]
   Float_t         jetMinDeltaRPVTracks_wp[100];   //[nJets]
   UInt_t          nFatJets;
   Float_t         fatJetE[6];   //[nFatJets]
   Float_t         fatJetPt[6];   //[nFatJets]
   Float_t         fatJetEta[6];   //[nFatJets]
   Float_t         fatJetPhi[6];   //[nFatJets]
   Float_t         fatJetCorrectedPt[6];   //[nFatJets]
   Float_t         fatJetPrunedM[6];   //[nFatJets]
   Float_t         fatJetTrimmedM[6];   //[nFatJets]
   Float_t         fatJetFilteredM[6];   //[nFatJets]
   Float_t         fatJetSoftDropM[6];   //[nFatJets]
   Float_t         fatJetCorrectedSoftDropM[6];   //[nFatJets]
   Float_t         fatJetUncorrectedSoftDropM[6];   //[nFatJets]
   Float_t         fatJetTau1[6];   //[nFatJets]
   Float_t         fatJetTau2[6];   //[nFatJets]
   Float_t         fatJetTau3[6];   //[nFatJets]
   Float_t         fatJetMaxSubjetCSV[6];   //[nFatJets]
   Bool_t          fatJetPassIDLoose[6];   //[nFatJets]
   Bool_t          fatJetPassIDTight[6];   //[nFatJets]
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
   Float_t         genJetE[100];   //[nGenJets]
   Float_t         genJetPt[100];   //[nGenJets]
   Float_t         genJetEta[100];   //[nGenJets]
   Float_t         genJetPhi[100];   //[nGenJets]
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
   string          *lheComments;
   vector<float>   *scaleWeights;
   vector<float>   *pdfWeights;
   vector<float>   *alphasWeights;
   Int_t           nGenParticle;
   Int_t           gParticleMotherId[5000];   //[nGenParticle]
   Int_t           gParticleMotherIndex[5000];   //[nGenParticle]
   Int_t           gParticleId[5000];   //[nGenParticle]
   Int_t           gParticleStatus[5000];   //[nGenParticle]
   Float_t         gParticleE[5000];   //[nGenParticle]
   Float_t         gParticlePt[5000];   //[nGenParticle]
   Float_t         gParticleEta[5000];   //[nGenParticle]
   Float_t         gParticlePhi[5000];   //[nGenParticle]
   Float_t         gParticleProdVertexX[5000];   //[nGenParticle]
   Float_t         gParticleProdVertexY[5000];   //[nGenParticle]
   Float_t         gParticleProdVertexZ[5000];   //[nGenParticle]
   Float_t         gParticleDecayVertexX[5000];   //[nGenParticle]
   Float_t         gParticleDecayVertexY[5000];   //[nGenParticle]
   Float_t         gParticleDecayVertexZ[5000];   //[nGenParticle]
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
   TBranch        *b_eleE_SC;   //!
   TBranch        *b_eleSigmaIetaIeta;   //!
   TBranch        *b_eleFull5x5SigmaIetaIeta;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_ele_dEta;   //!
   TBranch        *b_ele_dPhi;   //!
   TBranch        *b_ele_HoverE;   //!
   TBranch        *b_ele_d0;   //!
   TBranch        *b_ele_dZ;   //!
   TBranch        *b_ele_passCutBasedIDTight;   //!
   TBranch        *b_ele_passCutBasedIDLoose;   //!

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
   TBranch        *b_cscRechitClusterTimeSpread;   //!
   TBranch        *b_cscRechitClusterGenMuonDeltaR;   //!
   TBranch        *b_cscRechitClusterMajorAxis;   //!
   TBranch        *b_cscRechitClusterMinorAxis;   //!
   TBranch        *b_cscRechitClusterEtaPhiSpread;   //!
   TBranch        *b_cscRechitClusterPhiSpread;   //!
   TBranch        *b_cscRechitClusterEtaSpread;   //!
   TBranch        *b_cscRechitClusterXSpread;   //!
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
   TBranch        *b_track_bestVertexIndex;   //!
   TBranch        *b_track_nMissingInnerHits;   //!
   TBranch        *b_track_nMissingOuterHits;   //!
   TBranch        *b_track_dxyToBS;   //!
   TBranch        *b_track_dxyErr;   //!
   TBranch        *b_track_dzToPV;   //!
   TBranch        *b_track_dzErr;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
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

   TBranch        *b_jetElectronEnergyFraction;   //[nJets]
  TBranch        *b_jetPhotonEnergyFraction;   //[nJets]
  TBranch        *b_jetChargedHadronMultiplicity;   //[nJets]
  TBranch        *b_jetNeutralHadronMultiplicity;   //[nJets]
  TBranch        *b_jetPhotonMultiplicity;   //[nJets]
  TBranch        *b_jetElectronMultiplicity;   //[nJets]
  TBranch        *b_jetMuonMultiplicity;   //[nJets]
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
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef llp_event_cxx
llp_event::llp_event(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/mnt/hadoop/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p16/MC_Fall18/v1/sixie/WminusH_HToSSTobbbb_ms55_pl10000_ev11000/Run2_Run2_displacedJetMuonNtupler_V1p16_MC_Fall18_apresyan-crab_PrivateProduction_Fall18_DR_step2_WminusH_HToSSTobbbb_ms55_pl10000_batch1_v1_v1_v1/200324_213938/0000/displacedJetMuon_ntupler_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/mnt/hadoop/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p16/MC_Fall18/v1/sixie/WminusH_HToSSTobbbb_ms55_pl10000_ev11000/Run2_Run2_displacedJetMuonNtupler_V1p16_MC_Fall18_apresyan-crab_PrivateProduction_Fall18_DR_step2_WminusH_HToSSTobbbb_ms55_pl10000_batch1_v1_v1_v1/200324_213938/0000/displacedJetMuon_ntupler_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/mnt/hadoop/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p16/MC_Fall18/v1/sixie/WminusH_HToSSTobbbb_ms55_pl10000_ev11000/Run2_Run2_displacedJetMuonNtupler_V1p16_MC_Fall18_apresyan-crab_PrivateProduction_Fall18_DR_step2_WminusH_HToSSTobbbb_ms55_pl10000_batch1_v1_v1_v1/200324_213938/0000/displacedJetMuon_ntupler_1.root:/ntuples");
      dir->GetObject("llp",tree);

   }
   Init(tree);
}

llp_event::~llp_event()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t llp_event::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t llp_event::LoadTree(Long64_t entry)
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

void llp_event::Init(TTree *tree)
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
   fChain->SetBranchAddress("eleE_SC", eleE_SC, &b_eleE_SC);
   fChain->SetBranchAddress("eleSigmaIetaIeta", eleSigmaIetaIeta, &b_eleSigmaIetaIeta);
   fChain->SetBranchAddress("eleFull5x5SigmaIetaIeta", eleFull5x5SigmaIetaIeta, &b_eleFull5x5SigmaIetaIeta);
   fChain->SetBranchAddress("eleR9", eleR9, &b_eleR9);
   fChain->SetBranchAddress("ele_dEta", ele_dEta, &b_ele_dEta);
   fChain->SetBranchAddress("ele_dPhi", ele_dPhi, &b_ele_dPhi);
   fChain->SetBranchAddress("ele_HoverE", ele_HoverE, &b_ele_HoverE);
   fChain->SetBranchAddress("ele_d0", ele_d0, &b_ele_d0);
   fChain->SetBranchAddress("ele_dZ", ele_dZ, &b_ele_dZ);
   fChain->SetBranchAddress("ele_passCutBasedIDLoose", ele_passCutBasedIDLoose, &b_ele_passCutBasedIDLoose);
   fChain->SetBranchAddress("ele_passCutBasedIDTight", ele_passCutBasedIDTight, &b_ele_passCutBasedIDTight);
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
   fChain->SetBranchAddress("nTaus", &nTaus, &b_nTaus);
   fChain->SetBranchAddress("tauE", &tauE, &b_tauE);
   fChain->SetBranchAddress("tauPt", &tauPt, &b_tauPt);
   fChain->SetBranchAddress("tauEta", &tauEta, &b_tauEta);
   fChain->SetBranchAddress("tauPhi", &tauPhi, &b_tauPhi);
   fChain->SetBranchAddress("tau_IsLoose", &tau_IsLoose, &b_tau_IsLoose);
   fChain->SetBranchAddress("tau_IsMedium", &tau_IsMedium, &b_tau_IsMedium);
   fChain->SetBranchAddress("tau_IsTight", &tau_IsTight, &b_tau_IsTight);
   fChain->SetBranchAddress("tau_passEleVetoLoose", &tau_passEleVetoLoose, &b_tau_passEleVetoLoose);
   fChain->SetBranchAddress("tau_passEleVetoMedium", &tau_passEleVetoMedium, &b_tau_passEleVetoMedium);
   fChain->SetBranchAddress("tau_passEleVetoTight", &tau_passEleVetoTight, &b_tau_passEleVetoTight);
   fChain->SetBranchAddress("tau_passMuVetoLoose", &tau_passMuVetoLoose, &b_tau_passMuVetoLoose);
   fChain->SetBranchAddress("tau_passMuVetoMedium", &tau_passMuVetoMedium, &b_tau_passMuVetoMedium);
   fChain->SetBranchAddress("tau_passMuVetoTight", &tau_passMuVetoTight, &b_tau_passMuVetoTight);
   fChain->SetBranchAddress("tau_ID", &tau_ID, &b_tau_ID);
   fChain->SetBranchAddress("tau_combinedIsoDeltaBetaCorr3Hits", &tau_combinedIsoDeltaBetaCorr3Hits, &b_tau_combinedIsoDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_chargedIsoPtSum", &tau_chargedIsoPtSum, &b_tau_chargedIsoPtSum);
   fChain->SetBranchAddress("tau_neutralIsoPtSum", &tau_neutralIsoPtSum, &b_tau_neutralIsoPtSum);
   fChain->SetBranchAddress("tau_puCorrPtSum", &tau_puCorrPtSum, &b_tau_puCorrPtSum);
   fChain->SetBranchAddress("tau_eleVetoMVA", &tau_eleVetoMVA, &b_tau_eleVetoMVA);
   fChain->SetBranchAddress("tau_eleVetoCategory", &tau_eleVetoCategory, &b_tau_eleVetoCategory);
   fChain->SetBranchAddress("tau_muonVetoMVA", &tau_muonVetoMVA, &b_tau_muonVetoMVA);
   fChain->SetBranchAddress("tau_isoMVAnewDMwLT", &tau_isoMVAnewDMwLT, &b_tau_isoMVAnewDMwLT);
   fChain->SetBranchAddress("tau_isoMVAnewDMwoLT", &tau_isoMVAnewDMwoLT, &b_tau_isoMVAnewDMwoLT);
   fChain->SetBranchAddress("tau_leadCandPt", &tau_leadCandPt, &b_tau_leadCandPt);
   fChain->SetBranchAddress("tau_leadCandID", &tau_leadCandID, &b_tau_leadCandID);
   fChain->SetBranchAddress("tau_leadChargedHadrCandPt", &tau_leadChargedHadrCandPt, &b_tau_leadChargedHadrCandPt);
   fChain->SetBranchAddress("tau_leadChargedHadrCandID", &tau_leadChargedHadrCandID, &b_tau_leadChargedHadrCandID);
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
   fChain->SetBranchAddress("nCscWireDigis", &nCscWireDigis, &b_nCscWireDigis);
   fChain->SetBranchAddress("nCscStripDigis", &nCscStripDigis, &b_nCscStripDigis);
   fChain->SetBranchAddress("nCscSeg", &nCscSeg, &b_nCscSeg);
   fChain->SetBranchAddress("cscSegPhi", cscSegPhi, &b_cscSegPhi);
   fChain->SetBranchAddress("cscSegEta", cscSegEta, &b_cscSegEta);
   fChain->SetBranchAddress("cscSegX", cscSegX, &b_cscSegX);
   fChain->SetBranchAddress("cscSegY", cscSegY, &b_cscSegY);
   fChain->SetBranchAddress("cscSegZ", cscSegZ, &b_cscSegZ);
   fChain->SetBranchAddress("cscSegT", cscSegT, &b_cscSegT);
   fChain->SetBranchAddress("cscSegChi2", cscSegChi2, &b_cscSegChi2);
   fChain->SetBranchAddress("cscSegChamber", cscSegChamber, &b_cscSegChamber);
   fChain->SetBranchAddress("cscSegStation", cscSegStation, &b_cscSegStation);
   fChain->SetBranchAddress("cscSegNRecHits", cscSegNRecHits, &b_cscSegNRecHits);
   fChain->SetBranchAddress("ncscRechits", &ncscRechits, &b_ncscRechits);
   fChain->SetBranchAddress("cscRechitsPhi", cscRechitsPhi, &b_cscRechitsPhi);
   fChain->SetBranchAddress("cscRechitsEta", cscRechitsEta, &b_cscRechitsEta);
   fChain->SetBranchAddress("cscRechitsX", cscRechitsX, &b_cscRechitsX);
   fChain->SetBranchAddress("cscRechitsY", cscRechitsY, &b_cscRechitsY);
   fChain->SetBranchAddress("cscRechitsZ", cscRechitsZ, &b_cscRechitsZ);
   fChain->SetBranchAddress("cscRechitsE", cscRechitsE, &b_cscRechitsE);
   fChain->SetBranchAddress("cscRechitsTpeak", cscRechitsTpeak, &b_cscRechitsTpeak);
   fChain->SetBranchAddress("cscRechitsTwire", cscRechitsTwire, &b_cscRechitsTwire);
   fChain->SetBranchAddress("cscRechitsQuality", cscRechitsQuality, &b_cscRechitsQuality);
   fChain->SetBranchAddress("cscRechitsChamber", cscRechitsChamber, &b_cscRechitsChamber);
   fChain->SetBranchAddress("cscRechitsStation", cscRechitsStation, &b_cscRechitsStation);
   fChain->SetBranchAddress("cscRechitsChannels", cscRechitsChannels, &b_cscRechitsChannels);
   fChain->SetBranchAddress("cscRechitsNStrips", cscRechitsNStrips, &b_cscRechitsNStrips);
   fChain->SetBranchAddress("cscRechitsHitWire", cscRechitsHitWire, &b_cscRechitsHitWire);
   fChain->SetBranchAddress("cscRechitsWGroupsBX", cscRechitsWGroupsBX, &b_cscRechitsWGroupsBX);
   fChain->SetBranchAddress("cscRechitsNWireGroups", cscRechitsNWireGroups, &b_cscRechitsNWireGroups);
   fChain->SetBranchAddress("cscRechitsDetId", cscRechitsDetId, &b_cscRechitsDetId);
   fChain->SetBranchAddress("nCscRechitClusters", &nCscRechitClusters, &b_nCscRechitClusters);
   fChain->SetBranchAddress("cscRechitCluster_match_cscSegCluster_minDeltaR", cscRechitCluster_match_cscSegCluster_minDeltaR, &b_cscRechitCluster_match_cscSegCluster_minDeltaR);
   fChain->SetBranchAddress("cscRechitCluster_match_cscSegCluster_index", cscRechitCluster_match_cscSegCluster_index, &b_cscRechitCluster_match_cscSegCluster_index);
   fChain->SetBranchAddress("cscRechitCluster_match_gParticle_minDeltaR", cscRechitCluster_match_gParticle_minDeltaR, &b_cscRechitCluster_match_gParticle_minDeltaR);
   fChain->SetBranchAddress("cscRechitCluster_match_gParticle_index", cscRechitCluster_match_gParticle_index, &b_cscRechitCluster_match_gParticle_index);
   fChain->SetBranchAddress("cscRechitCluster_match_gParticle_id", cscRechitCluster_match_gParticle_id, &b_cscRechitCluster_match_gParticle_id);
   fChain->SetBranchAddress("cscRechitClusterX", cscRechitClusterX, &b_cscRechitClusterX);
   fChain->SetBranchAddress("cscRechitClusterY", cscRechitClusterY, &b_cscRechitClusterY);
   fChain->SetBranchAddress("cscRechitClusterZ", cscRechitClusterZ, &b_cscRechitClusterZ);
   fChain->SetBranchAddress("cscRechitClusterTime", cscRechitClusterTime, &b_cscRechitClusterTime);
   fChain->SetBranchAddress("cscRechitClusterTimeSpread", cscRechitClusterTimeSpread, &b_cscRechitClusterTimeSpread);
   fChain->SetBranchAddress("cscRechitClusterGenMuonDeltaR", cscRechitClusterGenMuonDeltaR, &b_cscRechitClusterGenMuonDeltaR);
   fChain->SetBranchAddress("cscRechitClusterMajorAxis", cscRechitClusterMajorAxis, &b_cscRechitClusterMajorAxis);
   fChain->SetBranchAddress("cscRechitClusterMinorAxis", cscRechitClusterMinorAxis, &b_cscRechitClusterMinorAxis);
   fChain->SetBranchAddress("cscRechitClusterEtaPhiSpread", cscRechitClusterEtaPhiSpread, &b_cscRechitClusterEtaPhiSpread);
   fChain->SetBranchAddress("cscRechitClusterPhiSpread", cscRechitClusterPhiSpread, &b_cscRechitClusterPhiSpread);
   fChain->SetBranchAddress("cscRechitClusterEtaSpread", cscRechitClusterEtaSpread, &b_cscRechitClusterEtaSpread);
   fChain->SetBranchAddress("cscRechitClusterXSpread", cscRechitClusterXSpread, &b_cscRechitClusterXSpread);
   fChain->SetBranchAddress("cscRechitClusterYSpread", cscRechitClusterYSpread, &b_cscRechitClusterYSpread);
   fChain->SetBranchAddress("cscRechitClusterZSpread", cscRechitClusterZSpread, &b_cscRechitClusterZSpread);
   fChain->SetBranchAddress("cscRechitClusterPhi", cscRechitClusterPhi, &b_cscRechitClusterPhi);
   fChain->SetBranchAddress("cscRechitClusterEta", cscRechitClusterEta, &b_cscRechitClusterEta);
   fChain->SetBranchAddress("cscRechitClusterJetVetoPt", cscRechitClusterJetVetoPt, &b_cscRechitClusterJetVetoPt);
   fChain->SetBranchAddress("cscRechitClusterJetVetoE", cscRechitClusterJetVetoE, &b_cscRechitClusterJetVetoE);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoPt", cscRechitClusterMuonVetoPt, &b_cscRechitClusterMuonVetoPt);
   fChain->SetBranchAddress("cscRechitClusterMuonVetoE", cscRechitClusterMuonVetoE, &b_cscRechitClusterMuonVetoE);
   fChain->SetBranchAddress("cscRechitClusterCaloJetVeto", cscRechitClusterCaloJetVeto, &b_cscRechitClusterCaloJetVeto);
   fChain->SetBranchAddress("cscRechitClusterSize", cscRechitClusterSize, &b_cscRechitClusterSize);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus11", cscRechitClusterNRechitChamberPlus11, &b_cscRechitClusterNRechitChamberPlus11);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus12", cscRechitClusterNRechitChamberPlus12, &b_cscRechitClusterNRechitChamberPlus12);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus13", cscRechitClusterNRechitChamberPlus13, &b_cscRechitClusterNRechitChamberPlus13);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus21", cscRechitClusterNRechitChamberPlus21, &b_cscRechitClusterNRechitChamberPlus21);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus22", cscRechitClusterNRechitChamberPlus22, &b_cscRechitClusterNRechitChamberPlus22);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus31", cscRechitClusterNRechitChamberPlus31, &b_cscRechitClusterNRechitChamberPlus31);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus32", cscRechitClusterNRechitChamberPlus32, &b_cscRechitClusterNRechitChamberPlus32);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus41", cscRechitClusterNRechitChamberPlus41, &b_cscRechitClusterNRechitChamberPlus41);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus42", cscRechitClusterNRechitChamberPlus42, &b_cscRechitClusterNRechitChamberPlus42);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus11", cscRechitClusterNRechitChamberMinus11, &b_cscRechitClusterNRechitChamberMinus11);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus12", cscRechitClusterNRechitChamberMinus12, &b_cscRechitClusterNRechitChamberMinus12);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus13", cscRechitClusterNRechitChamberMinus13, &b_cscRechitClusterNRechitChamberMinus13);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus21", cscRechitClusterNRechitChamberMinus21, &b_cscRechitClusterNRechitChamberMinus21);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus22", cscRechitClusterNRechitChamberMinus22, &b_cscRechitClusterNRechitChamberMinus22);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus31", cscRechitClusterNRechitChamberMinus31, &b_cscRechitClusterNRechitChamberMinus31);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus32", cscRechitClusterNRechitChamberMinus32, &b_cscRechitClusterNRechitChamberMinus32);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus41", cscRechitClusterNRechitChamberMinus41, &b_cscRechitClusterNRechitChamberMinus41);
   fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus42", cscRechitClusterNRechitChamberMinus42, &b_cscRechitClusterNRechitChamberMinus42);
   fChain->SetBranchAddress("cscRechitClusterMe11Ratio", cscRechitClusterMe11Ratio, &b_cscRechitClusterMe11Ratio);
   fChain->SetBranchAddress("cscRechitClusterMe12Ratio", cscRechitClusterMe12Ratio, &b_cscRechitClusterMe12Ratio);
   fChain->SetBranchAddress("cscRechitClusterNStation", cscRechitClusterNStation, &b_cscRechitClusterNStation);
   fChain->SetBranchAddress("cscRechitClusterMaxStation", cscRechitClusterMaxStation, &b_cscRechitClusterMaxStation);
   fChain->SetBranchAddress("cscRechitClusterMaxStationRatio", cscRechitClusterMaxStationRatio, &b_cscRechitClusterMaxStationRatio);
   fChain->SetBranchAddress("cscRechitClusterNChamber", cscRechitClusterNChamber, &b_cscRechitClusterNChamber);
   fChain->SetBranchAddress("cscRechitClusterMaxChamber", cscRechitClusterMaxChamber, &b_cscRechitClusterMaxChamber);
   fChain->SetBranchAddress("cscRechitClusterMaxChamberRatio", cscRechitClusterMaxChamberRatio, &b_cscRechitClusterMaxChamberRatio);
   fChain->SetBranchAddress("cscRechitClusterVertexR", cscRechitClusterVertexR, &b_cscRechitClusterVertexR);
   fChain->SetBranchAddress("cscRechitClusterVertexZ", cscRechitClusterVertexZ, &b_cscRechitClusterVertexZ);
   fChain->SetBranchAddress("cscRechitClusterVertexDis", cscRechitClusterVertexDis, &b_cscRechitClusterVertexDis);
   fChain->SetBranchAddress("cscRechitClusterVertexChi2", cscRechitClusterVertexChi2, &b_cscRechitClusterVertexChi2);
   fChain->SetBranchAddress("cscRechitClusterVertexN1", cscRechitClusterVertexN1, &b_cscRechitClusterVertexN1);
   fChain->SetBranchAddress("cscRechitClusterVertexN5", cscRechitClusterVertexN5, &b_cscRechitClusterVertexN5);
   fChain->SetBranchAddress("cscRechitClusterVertexN10", cscRechitClusterVertexN10, &b_cscRechitClusterVertexN10);
   fChain->SetBranchAddress("cscRechitClusterVertexN15", cscRechitClusterVertexN15, &b_cscRechitClusterVertexN15);
   fChain->SetBranchAddress("cscRechitClusterVertexN20", cscRechitClusterVertexN20, &b_cscRechitClusterVertexN20);
   fChain->SetBranchAddress("cscRechitClusterVertexN", cscRechitClusterVertexN, &b_cscRechitClusterVertexN);
   fChain->SetBranchAddress("nCscSegClusters", &nCscSegClusters, &b_nCscSegClusters);
   fChain->SetBranchAddress("cscSegCluster_match_gParticle_minDeltaR", cscSegCluster_match_gParticle_minDeltaR, &b_cscSegCluster_match_gParticle_minDeltaR);
   fChain->SetBranchAddress("cscSegCluster_match_gParticle_index", cscSegCluster_match_gParticle_index, &b_cscSegCluster_match_gParticle_index);
   fChain->SetBranchAddress("cscSegCluster_match_gParticle_id", cscSegCluster_match_gParticle_id, &b_cscSegCluster_match_gParticle_id);
   fChain->SetBranchAddress("cscSegClusterX", cscSegClusterX, &b_cscSegClusterX);
   fChain->SetBranchAddress("cscSegClusterY", cscSegClusterY, &b_cscSegClusterY);
   fChain->SetBranchAddress("cscSegClusterZ", cscSegClusterZ, &b_cscSegClusterZ);
   fChain->SetBranchAddress("cscSegClusterTime", cscSegClusterTime, &b_cscSegClusterTime);
   fChain->SetBranchAddress("cscSegClusterTimeSpread", cscSegClusterTimeSpread, &b_cscSegClusterTimeSpread);
   fChain->SetBranchAddress("cscSegClusterGenMuonDeltaR", cscSegClusterGenMuonDeltaR, &b_cscSegClusterGenMuonDeltaR);
   fChain->SetBranchAddress("cscSegClusterMajorAxis", cscSegClusterMajorAxis, &b_cscSegClusterMajorAxis);
   fChain->SetBranchAddress("cscSegClusterMinorAxis", cscSegClusterMinorAxis, &b_cscSegClusterMinorAxis);
   fChain->SetBranchAddress("cscSegClusterEtaPhiSpread", cscSegClusterEtaPhiSpread, &b_cscSegClusterEtaPhiSpread);
   fChain->SetBranchAddress("cscSegClusterPhiSpread", cscSegClusterPhiSpread, &b_cscSegClusterPhiSpread);
   fChain->SetBranchAddress("cscSegClusterEtaSpread", cscSegClusterEtaSpread, &b_cscSegClusterEtaSpread);
   fChain->SetBranchAddress("cscSegClusterXSpread", cscSegClusterXSpread, &b_cscSegClusterXSpread);
   fChain->SetBranchAddress("cscSegClusterYSpread", cscSegClusterYSpread, &b_cscSegClusterYSpread);
   fChain->SetBranchAddress("cscSegClusterZSpread", cscSegClusterZSpread, &b_cscSegClusterZSpread);
   fChain->SetBranchAddress("cscSegClusterPhi", cscSegClusterPhi, &b_cscSegClusterPhi);
   fChain->SetBranchAddress("cscSegClusterEta", cscSegClusterEta, &b_cscSegClusterEta);
   fChain->SetBranchAddress("cscSegClusterJetVetoPt", cscSegClusterJetVetoPt, &b_cscSegClusterJetVetoPt);
   fChain->SetBranchAddress("cscSegClusterJetVetoE", cscSegClusterJetVetoE, &b_cscSegClusterJetVetoE);
   fChain->SetBranchAddress("cscSegClusterMuonVetoPt", cscSegClusterMuonVetoPt, &b_cscSegClusterMuonVetoPt);
   fChain->SetBranchAddress("cscSegClusterMuonVetoE", cscSegClusterMuonVetoE, &b_cscSegClusterMuonVetoE);
   fChain->SetBranchAddress("cscSegClusterCaloJetVeto", cscSegClusterCaloJetVeto, &b_cscSegClusterCaloJetVeto);
   fChain->SetBranchAddress("cscSegClusterSize", cscSegClusterSize, &b_cscSegClusterSize);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberPlus11", cscSegClusterNSegmentChamberPlus11, &b_cscSegClusterNSegmentChamberPlus11);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberPlus12", cscSegClusterNSegmentChamberPlus12, &b_cscSegClusterNSegmentChamberPlus12);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberPlus13", cscSegClusterNSegmentChamberPlus13, &b_cscSegClusterNSegmentChamberPlus13);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberPlus21", cscSegClusterNSegmentChamberPlus21, &b_cscSegClusterNSegmentChamberPlus21);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberPlus22", cscSegClusterNSegmentChamberPlus22, &b_cscSegClusterNSegmentChamberPlus22);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberPlus31", cscSegClusterNSegmentChamberPlus31, &b_cscSegClusterNSegmentChamberPlus31);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberPlus32", cscSegClusterNSegmentChamberPlus32, &b_cscSegClusterNSegmentChamberPlus32);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberPlus41", cscSegClusterNSegmentChamberPlus41, &b_cscSegClusterNSegmentChamberPlus41);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberPlus42", cscSegClusterNSegmentChamberPlus42, &b_cscSegClusterNSegmentChamberPlus42);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberMinus11", cscSegClusterNSegmentChamberMinus11, &b_cscSegClusterNSegmentChamberMinus11);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberMinus12", cscSegClusterNSegmentChamberMinus12, &b_cscSegClusterNSegmentChamberMinus12);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberMinus13", cscSegClusterNSegmentChamberMinus13, &b_cscSegClusterNSegmentChamberMinus13);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberMinus21", cscSegClusterNSegmentChamberMinus21, &b_cscSegClusterNSegmentChamberMinus21);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberMinus22", cscSegClusterNSegmentChamberMinus22, &b_cscSegClusterNSegmentChamberMinus22);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberMinus31", cscSegClusterNSegmentChamberMinus31, &b_cscSegClusterNSegmentChamberMinus31);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberMinus32", cscSegClusterNSegmentChamberMinus32, &b_cscSegClusterNSegmentChamberMinus32);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberMinus41", cscSegClusterNSegmentChamberMinus41, &b_cscSegClusterNSegmentChamberMinus41);
   fChain->SetBranchAddress("cscSegClusterNSegmentChamberMinus42", cscSegClusterNSegmentChamberMinus42, &b_cscSegClusterNSegmentChamberMinus42);
   fChain->SetBranchAddress("cscSegClusterMe11Ratio", cscSegClusterMe11Ratio, &b_cscSegClusterMe11Ratio);
   fChain->SetBranchAddress("cscSegClusterMe12Ratio", cscSegClusterMe12Ratio, &b_cscSegClusterMe12Ratio);
   fChain->SetBranchAddress("cscSegClusterNStation", cscSegClusterNStation, &b_cscSegClusterNStation);
   fChain->SetBranchAddress("cscSegClusterMaxStation", cscSegClusterMaxStation, &b_cscSegClusterMaxStation);
   fChain->SetBranchAddress("cscSegClusterMaxStationRatio", cscSegClusterMaxStationRatio, &b_cscSegClusterMaxStationRatio);
   fChain->SetBranchAddress("cscSegClusterNChamber", cscSegClusterNChamber, &b_cscSegClusterNChamber);
   fChain->SetBranchAddress("cscSegClusterMaxChamber", cscSegClusterMaxChamber, &b_cscSegClusterMaxChamber);
   fChain->SetBranchAddress("cscSegClusterMaxChamberRatio", cscSegClusterMaxChamberRatio, &b_cscSegClusterMaxChamberRatio);
   fChain->SetBranchAddress("cscSegClusterVertexR", cscSegClusterVertexR, &b_cscSegClusterVertexR);
   fChain->SetBranchAddress("cscSegClusterVertexZ", cscSegClusterVertexZ, &b_cscSegClusterVertexZ);
   fChain->SetBranchAddress("cscSegClusterVertexDis", cscSegClusterVertexDis, &b_cscSegClusterVertexDis);
   fChain->SetBranchAddress("cscSegClusterVertexChi2", cscSegClusterVertexChi2, &b_cscSegClusterVertexChi2);
   fChain->SetBranchAddress("cscSegClusterVertexN1", cscSegClusterVertexN1, &b_cscSegClusterVertexN1);
   fChain->SetBranchAddress("cscSegClusterVertexN5", cscSegClusterVertexN5, &b_cscSegClusterVertexN5);
   fChain->SetBranchAddress("cscSegClusterVertexN10", cscSegClusterVertexN10, &b_cscSegClusterVertexN10);
   fChain->SetBranchAddress("cscSegClusterVertexN15", cscSegClusterVertexN15, &b_cscSegClusterVertexN15);
   fChain->SetBranchAddress("cscSegClusterVertexN20", cscSegClusterVertexN20, &b_cscSegClusterVertexN20);
   fChain->SetBranchAddress("cscSegClusterVertexN", cscSegClusterVertexN, &b_cscSegClusterVertexN);
   fChain->SetBranchAddress("nDtRechits", &nDtRechits, &b_nDtRechits);
   fChain->SetBranchAddress("dtRechitX", dtRechitX, &b_dtRechitX);
   fChain->SetBranchAddress("dtRechitY", dtRechitY, &b_dtRechitY);
   fChain->SetBranchAddress("dtRechitZ", dtRechitZ, &b_dtRechitZ);
   fChain->SetBranchAddress("dtRechitEta", dtRechitEta, &b_dtRechitEta);
   fChain->SetBranchAddress("dtRechitPhi", dtRechitPhi, &b_dtRechitPhi);
   fChain->SetBranchAddress("dtRechitTime", dtRechitTime, &b_dtRechitTime);
   fChain->SetBranchAddress("dtRechitStation", dtRechitStation, &b_dtRechitStation);
   fChain->SetBranchAddress("dtRechitWheel", dtRechitWheel, &b_dtRechitWheel);
   fChain->SetBranchAddress("nDtRechitClusters", &nDtRechitClusters, &b_nDtRechitClusters);
   fChain->SetBranchAddress("dtRechitCluster_match_gParticle_minDeltaR", dtRechitCluster_match_gParticle_minDeltaR, &b_dtRechitCluster_match_gParticle_minDeltaR);
   fChain->SetBranchAddress("dtRechitCluster_match_gParticle_index", dtRechitCluster_match_gParticle_index, &b_dtRechitCluster_match_gParticle_index);
   fChain->SetBranchAddress("dtRechitCluster_match_gParticle_id", dtRechitCluster_match_gParticle_id, &b_dtRechitCluster_match_gParticle_id);
   fChain->SetBranchAddress("dtRechitClusterX", dtRechitClusterX, &b_dtRechitClusterX);
   fChain->SetBranchAddress("dtRechitClusterY", dtRechitClusterY, &b_dtRechitClusterY);
   fChain->SetBranchAddress("dtRechitClusterZ", dtRechitClusterZ, &b_dtRechitClusterZ);
   fChain->SetBranchAddress("dtRechitClusterTime", dtRechitClusterTime, &b_dtRechitClusterTime);
   fChain->SetBranchAddress("dtRechitClusterTimeSpread", dtRechitClusterTimeSpread, &b_dtRechitClusterTimeSpread);
   fChain->SetBranchAddress("dtRechitClusterGenMuonDeltaR", dtRechitClusterGenMuonDeltaR, &b_dtRechitClusterGenMuonDeltaR);
   fChain->SetBranchAddress("dtRechitClusterMajorAxis", dtRechitClusterMajorAxis, &b_dtRechitClusterMajorAxis);
   fChain->SetBranchAddress("dtRechitClusterMinorAxis", dtRechitClusterMinorAxis, &b_dtRechitClusterMinorAxis);
   fChain->SetBranchAddress("dtRechitClusterEtaPhiSpread", dtRechitClusterEtaPhiSpread, &b_dtRechitClusterEtaPhiSpread);
   fChain->SetBranchAddress("dtRechitClusterPhiSpread", dtRechitClusterPhiSpread, &b_dtRechitClusterPhiSpread);
   fChain->SetBranchAddress("dtRechitClusterEtaSpread", dtRechitClusterEtaSpread, &b_dtRechitClusterEtaSpread);
   fChain->SetBranchAddress("dtRechitClusterXSpread", dtRechitClusterXSpread, &b_dtRechitClusterXSpread);
   fChain->SetBranchAddress("dtRechitClusterYSpread", dtRechitClusterYSpread, &b_dtRechitClusterYSpread);
   fChain->SetBranchAddress("dtRechitClusterZSpread", dtRechitClusterZSpread, &b_dtRechitClusterZSpread);
   fChain->SetBranchAddress("dtRechitClusterPhi", dtRechitClusterPhi, &b_dtRechitClusterPhi);
   fChain->SetBranchAddress("dtRechitClusterEta", dtRechitClusterEta, &b_dtRechitClusterEta);
   fChain->SetBranchAddress("dtRechitClusterJetVetoPt", dtRechitClusterJetVetoPt, &b_dtRechitClusterJetVetoPt);
   fChain->SetBranchAddress("dtRechitClusterJetVetoE", dtRechitClusterJetVetoE, &b_dtRechitClusterJetVetoE);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoPt", dtRechitClusterMuonVetoPt, &b_dtRechitClusterMuonVetoPt);
   fChain->SetBranchAddress("dtRechitClusterMuonVetoE", dtRechitClusterMuonVetoE, &b_dtRechitClusterMuonVetoE);
   fChain->SetBranchAddress("dtRechitClusterCaloJetVeto", dtRechitClusterCaloJetVeto, &b_dtRechitClusterCaloJetVeto);
   fChain->SetBranchAddress("dtRechitClusterSize", dtRechitClusterSize, &b_dtRechitClusterSize);
   fChain->SetBranchAddress("dtRechitClusterNStation", dtRechitClusterNStation, &b_dtRechitClusterNStation);
   fChain->SetBranchAddress("dtRechitClusterMaxStation", dtRechitClusterMaxStation, &b_dtRechitClusterMaxStation);
   fChain->SetBranchAddress("dtRechitClusterMaxStationRatio", dtRechitClusterMaxStationRatio, &b_dtRechitClusterMaxStationRatio);
   fChain->SetBranchAddress("dtRechitClusterNChamber", dtRechitClusterNChamber, &b_dtRechitClusterNChamber);
   fChain->SetBranchAddress("dtRechitClusterMaxChamber", dtRechitClusterMaxChamber, &b_dtRechitClusterMaxChamber);
   fChain->SetBranchAddress("dtRechitClusterMaxChamberRatio", dtRechitClusterMaxChamberRatio, &b_dtRechitClusterMaxChamberRatio);
   fChain->SetBranchAddress("dtRechitClusterNSegmentStation1", dtRechitClusterNSegmentStation1, &b_dtRechitClusterNSegmentStation1);
   fChain->SetBranchAddress("dtRechitClusterNSegmentStation2", dtRechitClusterNSegmentStation2, &b_dtRechitClusterNSegmentStation2);
   fChain->SetBranchAddress("dtRechitClusterNSegmentStation3", dtRechitClusterNSegmentStation3, &b_dtRechitClusterNSegmentStation3);
   fChain->SetBranchAddress("dtRechitClusterNSegmentStation4", dtRechitClusterNSegmentStation4, &b_dtRechitClusterNSegmentStation4);
   fChain->SetBranchAddress("nRpc", &nRpc, &b_nRpc);
   fChain->SetBranchAddress("rpcPhi", rpcPhi, &b_rpcPhi);
   fChain->SetBranchAddress("rpcEta", rpcEta, &b_rpcEta);
   fChain->SetBranchAddress("rpcX", rpcX, &b_rpcX);
   fChain->SetBranchAddress("rpcY", rpcY, &b_rpcY);
   fChain->SetBranchAddress("rpcZ", rpcZ, &b_rpcZ);
   fChain->SetBranchAddress("rpcT", rpcT, &b_rpcT);
   fChain->SetBranchAddress("rpcBx", rpcBx, &b_rpcBx);
   fChain->SetBranchAddress("rpcTError", rpcTError, &b_rpcTError);
   fChain->SetBranchAddress("nDtSeg", &nDtSeg, &b_nDtSeg);
   fChain->SetBranchAddress("dtSegPhi", dtSegPhi, &b_dtSegPhi);
   fChain->SetBranchAddress("dtSegEta", dtSegEta, &b_dtSegEta);
   fChain->SetBranchAddress("dtSegX", dtSegX, &b_dtSegX);
   fChain->SetBranchAddress("dtSegY", dtSegY, &b_dtSegY);
   fChain->SetBranchAddress("dtSegZ", dtSegZ, &b_dtSegZ);
   fChain->SetBranchAddress("dtSegStation", dtSegStation, &b_dtSegStation);
   fChain->SetBranchAddress("dtSegWheel", dtSegWheel, &b_dtSegWheel);
   fChain->SetBranchAddress("dtSegTime", dtSegTime, &b_dtSegTime);
   fChain->SetBranchAddress("dtSegTimeError", dtSegTimeError, &b_dtSegTimeError);
   fChain->SetBranchAddress("nDtSegClusters", &nDtSegClusters, &b_nDtSegClusters);
   fChain->SetBranchAddress("dtSegCluster_match_gParticle_minDeltaR", dtSegCluster_match_gParticle_minDeltaR, &b_dtSegCluster_match_gParticle_minDeltaR);
   fChain->SetBranchAddress("dtSegCluster_match_gParticle_index", dtSegCluster_match_gParticle_index, &b_dtSegCluster_match_gParticle_index);
   fChain->SetBranchAddress("dtSegCluster_match_gParticle_id", dtSegCluster_match_gParticle_id, &b_dtSegCluster_match_gParticle_id);
   fChain->SetBranchAddress("dtSegClusterX", dtSegClusterX, &b_dtSegClusterX);
   fChain->SetBranchAddress("dtSegClusterY", dtSegClusterY, &b_dtSegClusterY);
   fChain->SetBranchAddress("dtSegClusterZ", dtSegClusterZ, &b_dtSegClusterZ);
   fChain->SetBranchAddress("dtSegClusterTime", dtSegClusterTime, &b_dtSegClusterTime);
   fChain->SetBranchAddress("dtSegClusterTimeSpread", dtSegClusterTimeSpread, &b_dtSegClusterTimeSpread);
   fChain->SetBranchAddress("dtSegClusterGenMuonDeltaR", dtSegClusterGenMuonDeltaR, &b_dtSegClusterGenMuonDeltaR);
   fChain->SetBranchAddress("dtSegClusterMajorAxis", dtSegClusterMajorAxis, &b_dtSegClusterMajorAxis);
   fChain->SetBranchAddress("dtSegClusterMinorAxis", dtSegClusterMinorAxis, &b_dtSegClusterMinorAxis);
   fChain->SetBranchAddress("dtSegClusterEtaPhiSpread", dtSegClusterEtaPhiSpread, &b_dtSegClusterEtaPhiSpread);
   fChain->SetBranchAddress("dtSegClusterPhiSpread", dtSegClusterPhiSpread, &b_dtSegClusterPhiSpread);
   fChain->SetBranchAddress("dtSegClusterEtaSpread", dtSegClusterEtaSpread, &b_dtSegClusterEtaSpread);
   fChain->SetBranchAddress("dtSegClusterXSpread", dtSegClusterXSpread, &b_dtSegClusterXSpread);
   fChain->SetBranchAddress("dtSegClusterYSpread", dtSegClusterYSpread, &b_dtSegClusterYSpread);
   fChain->SetBranchAddress("dtSegClusterZSpread", dtSegClusterZSpread, &b_dtSegClusterZSpread);
   fChain->SetBranchAddress("dtSegClusterPhi", dtSegClusterPhi, &b_dtSegClusterPhi);
   fChain->SetBranchAddress("dtSegClusterEta", dtSegClusterEta, &b_dtSegClusterEta);
   fChain->SetBranchAddress("dtSegClusterJetVetoPt", dtSegClusterJetVetoPt, &b_dtSegClusterJetVetoPt);
   fChain->SetBranchAddress("dtSegClusterJetVetoE", dtSegClusterJetVetoE, &b_dtSegClusterJetVetoE);
   fChain->SetBranchAddress("dtSegClusterMuonVetoPt", dtSegClusterMuonVetoPt, &b_dtSegClusterMuonVetoPt);
   fChain->SetBranchAddress("dtSegClusterMuonVetoE", dtSegClusterMuonVetoE, &b_dtSegClusterMuonVetoE);
   fChain->SetBranchAddress("dtSegClusterCaloJetVeto", dtSegClusterCaloJetVeto, &b_dtSegClusterCaloJetVeto);
   fChain->SetBranchAddress("dtSegClusterSize", dtSegClusterSize, &b_dtSegClusterSize);
   fChain->SetBranchAddress("dtSegClusterNSegmentStation1", dtSegClusterNSegmentStation1, &b_dtSegClusterNSegmentStation1);
   fChain->SetBranchAddress("dtSegClusterNSegmentStation2", dtSegClusterNSegmentStation2, &b_dtSegClusterNSegmentStation2);
   fChain->SetBranchAddress("dtSegClusterNSegmentStation3", dtSegClusterNSegmentStation3, &b_dtSegClusterNSegmentStation3);
   fChain->SetBranchAddress("dtSegClusterNSegmentStation4", dtSegClusterNSegmentStation4, &b_dtSegClusterNSegmentStation4);
   fChain->SetBranchAddress("dtSegClusterNStation", dtSegClusterNStation, &b_dtSegClusterNStation);
   fChain->SetBranchAddress("dtSegClusterMaxStation", dtSegClusterMaxStation, &b_dtSegClusterMaxStation);
   fChain->SetBranchAddress("dtSegClusterMaxStationRatio", dtSegClusterMaxStationRatio, &b_dtSegClusterMaxStationRatio);
   fChain->SetBranchAddress("dtSegClusterNChamber", dtSegClusterNChamber, &b_dtSegClusterNChamber);
   fChain->SetBranchAddress("dtSegClusterMaxChamber", dtSegClusterMaxChamber, &b_dtSegClusterMaxChamber);
   fChain->SetBranchAddress("dtSegClusterMaxChamberRatio", dtSegClusterMaxChamberRatio, &b_dtSegClusterMaxChamberRatio);
   fChain->SetBranchAddress("nHORechits", &nHORechits, &b_nHORechits);
   fChain->SetBranchAddress("hoRechit_Eta", &hoRechit_Eta, &b_hoRechit_Eta);
   fChain->SetBranchAddress("hoRechit_Phi", &hoRechit_Phi, &b_hoRechit_Phi);
   fChain->SetBranchAddress("hoRechit_E", &hoRechit_E, &b_hoRechit_E);
   fChain->SetBranchAddress("hoRechit_X", &hoRechit_X, &b_hoRechit_X);
   fChain->SetBranchAddress("hoRechit_Y", &hoRechit_Y, &b_hoRechit_Y);
   fChain->SetBranchAddress("hoRechit_Z", &hoRechit_Z, &b_hoRechit_Z);
   fChain->SetBranchAddress("hoRechit_T", &hoRechit_T, &b_hoRechit_T);
   fChain->SetBranchAddress("nRechits", &nRechits, &b_nRechits);
   fChain->SetBranchAddress("ecalRechit_Eta", ecalRechit_Eta, &b_ecalRechit_Eta);
   fChain->SetBranchAddress("ecalRechit_Phi", ecalRechit_Phi, &b_ecalRechit_Phi);
   fChain->SetBranchAddress("ecalRechit_E", ecalRechit_E, &b_ecalRechit_E);
   fChain->SetBranchAddress("ecalRechit_T", ecalRechit_T, &b_ecalRechit_T);
   fChain->SetBranchAddress("ecalRechit_E_Error", ecalRechit_E_Error, &b_ecalRechit_E_Error);
   fChain->SetBranchAddress("ecalRechit_T_Error", ecalRechit_T_Error, &b_ecalRechit_T_Error);
   fChain->SetBranchAddress("ecalRechit_kSaturatedflag", ecalRechit_kSaturatedflag, &b_ecalRechit_kSaturatedflag);
   fChain->SetBranchAddress("ecalRechit_kLeadingEdgeRecoveredflag", ecalRechit_kLeadingEdgeRecoveredflag, &b_ecalRechit_kLeadingEdgeRecoveredflag);
   fChain->SetBranchAddress("ecalRechit_kPoorRecoflag", ecalRechit_kPoorRecoflag, &b_ecalRechit_kPoorRecoflag);
   fChain->SetBranchAddress("ecalRechit_kWeirdflag", ecalRechit_kWeirdflag, &b_ecalRechit_kWeirdflag);
   fChain->SetBranchAddress("ecalRechit_kDiWeirdflag", ecalRechit_kDiWeirdflag, &b_ecalRechit_kDiWeirdflag);
   fChain->SetBranchAddress("nHBHERechits", &nHBHERechits, &b_nHBHERechits);
   fChain->SetBranchAddress("hbheRechit_Eta", hbheRechit_Eta, &b_hbheRechit_Eta);
   fChain->SetBranchAddress("hbheRechit_Phi", hbheRechit_Phi, &b_hbheRechit_Phi);
   fChain->SetBranchAddress("hbheRechit_E", hbheRechit_E, &b_hbheRechit_E);
   fChain->SetBranchAddress("hbheRechit_X", hbheRechit_X, &b_hbheRechit_X);
   fChain->SetBranchAddress("hbheRechit_Y", hbheRechit_Y, &b_hbheRechit_Y);
   fChain->SetBranchAddress("hbheRechit_Z", hbheRechit_Z, &b_hbheRechit_Z);
   fChain->SetBranchAddress("hbheRechit_T", hbheRechit_T, &b_hbheRechit_T);
   fChain->SetBranchAddress("hbheRechit_iEta", hbheRechit_iEta, &b_hbheRechit_iEta);
   fChain->SetBranchAddress("hbheRechit_iPhi", hbheRechit_iPhi, &b_hbheRechit_iPhi);
   fChain->SetBranchAddress("hbheRechit_depth", hbheRechit_depth, &b_hbheRechit_depth);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("track_Pt", track_Pt, &b_track_Pt);
   fChain->SetBranchAddress("track_Eta", track_Eta, &b_track_Eta);
   fChain->SetBranchAddress("track_Phi", track_Phi, &b_track_Phi);
   fChain->SetBranchAddress("track_bestVertexIndex", track_bestVertexIndex, &b_track_bestVertexIndex);
   fChain->SetBranchAddress("track_nMissingInnerHits", track_nMissingInnerHits, &b_track_nMissingInnerHits);
   fChain->SetBranchAddress("track_nMissingOuterHits", track_nMissingOuterHits, &b_track_nMissingOuterHits);
   fChain->SetBranchAddress("track_dxyToBS", track_dxyToBS, &b_track_dxyToBS);
   fChain->SetBranchAddress("track_dxyErr", track_dxyErr, &b_track_dxyErr);
   fChain->SetBranchAddress("track_dzToPV", track_dzToPV, &b_track_dzToPV);
   fChain->SetBranchAddress("track_dzErr", track_dzErr, &b_track_dzErr);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("jetE", jetE, &b_jetE);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
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

   fChain->SetBranchAddress("jetElectronEnergyFraction", jetElectronEnergyFraction, &b_jetElectronEnergyFraction);
   fChain->SetBranchAddress("jetPhotonEnergyFraction", jetPhotonEnergyFraction, &b_jetPhotonEnergyFraction);
   fChain->SetBranchAddress("jetChargedHadronMultiplicity", jetChargedHadronMultiplicity, &b_jetChargedHadronMultiplicity);
   fChain->SetBranchAddress("jetNeutralHadronMultiplicity", jetNeutralHadronMultiplicity, &b_jetNeutralHadronMultiplicity);
   fChain->SetBranchAddress("jetPhotonMultiplicity", jetPhotonMultiplicity, &b_jetPhotonMultiplicity);
   fChain->SetBranchAddress("jetMuonMultiplicity", jetMuonMultiplicity, &b_jetMuonMultiplicity);

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
   fChain->SetBranchAddress("jetAlphaMax_wp", jetAlphaMax_wp, &b_jetAlphaMax_wp);
   fChain->SetBranchAddress("jetBetaMax_wp", jetBetaMax_wp, &b_jetBetaMax_wp);
   fChain->SetBranchAddress("jetGammaMax_ET_wp", jetGammaMax_ET_wp, &b_jetGammaMax_ET_wp);
   fChain->SetBranchAddress("jetGammaMax_EM_wp", jetGammaMax_EM_wp, &b_jetGammaMax_EM_wp);
   fChain->SetBranchAddress("jetGammaMax_Hadronic_wp", jetGammaMax_Hadronic_wp, &b_jetGammaMax_Hadronic_wp);
   fChain->SetBranchAddress("jetGammaMax_wp", jetGammaMax_wp, &b_jetGammaMax_wp);
   fChain->SetBranchAddress("jetPtAllTracks_wp", jetPtAllTracks_wp, &b_jetPtAllTracks_wp);
   fChain->SetBranchAddress("jetPtAllPVTracks_wp", jetPtAllPVTracks_wp, &b_jetPtAllPVTracks_wp);
   fChain->SetBranchAddress("jetMedianTheta2D_wp", jetMedianTheta2D_wp, &b_jetMedianTheta2D_wp);
   fChain->SetBranchAddress("jetMedianIP_wp", jetMedianIP_wp, &b_jetMedianIP_wp);
   fChain->SetBranchAddress("jetMinDeltaRAllTracks_wp", jetMinDeltaRAllTracks_wp, &b_jetMinDeltaRAllTracks_wp);
   fChain->SetBranchAddress("jetMinDeltaRPVTracks_wp", jetMinDeltaRPVTracks_wp, &b_jetMinDeltaRPVTracks_wp);
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
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
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
   fChain->SetBranchAddress("Flag2_globalSuperTightHalo2016Filter", &Flag2_globalSuperTightHalo2016Filter, &b_Flag2_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag2_globalTightHalo2016Filter", &Flag2_globalTightHalo2016Filter, &b_Flag2_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag2_goodVertices", &Flag2_goodVertices, &b_Flag2_goodVertices);
   fChain->SetBranchAddress("Flag2_BadChargedCandidateFilter", &Flag2_BadChargedCandidateFilter, &b_Flag2_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag2_BadPFMuonFilter", &Flag2_BadPFMuonFilter, &b_Flag2_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag2_EcalDeadCellTriggerPrimitiveFilter", &Flag2_EcalDeadCellTriggerPrimitiveFilter, &b_Flag2_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag2_HBHENoiseFilter", &Flag2_HBHENoiseFilter, &b_Flag2_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag2_HBHEIsoNoiseFilter", &Flag2_HBHEIsoNoiseFilter, &b_Flag2_HBHEIsoNoiseFilter);
   fChain->SetBranchAddress("Flag2_ecalBadCalibFilter", &Flag2_ecalBadCalibFilter, &b_Flag2_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag2_eeBadScFilter", &Flag2_eeBadScFilter, &b_Flag2_eeBadScFilter);
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
   fChain->SetBranchAddress("gParticleEta", gParticleEta, &b_gParticleEta);
   fChain->SetBranchAddress("gParticlePhi", gParticlePhi, &b_gParticlePhi);
   fChain->SetBranchAddress("gParticleProdVertexX", gParticleProdVertexX, &b_gParticleProdVertexX);
   fChain->SetBranchAddress("gParticleProdVertexY", gParticleProdVertexY, &b_gParticleProdVertexY);
   fChain->SetBranchAddress("gParticleProdVertexZ", gParticleProdVertexZ, &b_gParticleProdVertexZ);
   fChain->SetBranchAddress("gParticleDecayVertexX", gParticleDecayVertexX, &b_gParticleDecayVertexX);
   fChain->SetBranchAddress("gParticleDecayVertexY", gParticleDecayVertexY, &b_gParticleDecayVertexY);
   fChain->SetBranchAddress("gParticleDecayVertexZ", gParticleDecayVertexZ, &b_gParticleDecayVertexZ);
   fChain->SetBranchAddress("gLLP_decay_vertex_x", gLLP_decay_vertex_x, &b_gLLP_decay_vertex_x);
   fChain->SetBranchAddress("gLLP_decay_vertex_y", gLLP_decay_vertex_y, &b_gLLP_decay_vertex_y);
   fChain->SetBranchAddress("gLLP_decay_vertex_z", gLLP_decay_vertex_z, &b_gLLP_decay_vertex_z);
   fChain->SetBranchAddress("gLLP_beta", gLLP_beta, &b_gLLP_beta);
   fChain->SetBranchAddress("gLLP_e", gLLP_e, &b_gLLP_e);
   fChain->SetBranchAddress("gLLP_pt", gLLP_pt, &b_gLLP_pt);
   fChain->SetBranchAddress("gLLP_eta", gLLP_eta, &b_gLLP_eta);
   fChain->SetBranchAddress("gLLP_phi", gLLP_phi, &b_gLLP_phi);
   fChain->SetBranchAddress("gLLP_csc", gLLP_csc, &b_gLLP_csc);
   fChain->SetBranchAddress("gLLP_dt", gLLP_dt, &b_gLLP_dt);
   fChain->SetBranchAddress("gLLP_travel_time", gLLP_travel_time, &b_gLLP_travel_time);
   fChain->SetBranchAddress("gLLP_daughter_id", gLLP_daughter_id, &b_gLLP_daughter_id);
   fChain->SetBranchAddress("gLLP_daughter_pt", gLLP_daughter_pt, &b_gLLP_daughter_pt);
   fChain->SetBranchAddress("gLLP_daughter_eta", gLLP_daughter_eta, &b_gLLP_daughter_eta);
   fChain->SetBranchAddress("gLLP_daughter_phi", gLLP_daughter_phi, &b_gLLP_daughter_phi);
   fChain->SetBranchAddress("gLLP_daughter_eta_ecalcorr", gLLP_daughter_eta_ecalcorr, &b_gLLP_daughter_eta_ecalcorr);
   fChain->SetBranchAddress("gLLP_daughter_phi_ecalcorr", gLLP_daughter_phi_ecalcorr, &b_gLLP_daughter_phi_ecalcorr);
   fChain->SetBranchAddress("gLLP_daughter_e", gLLP_daughter_e, &b_gLLP_daughter_e);
   fChain->SetBranchAddress("gLLP_daughter_mass", gLLP_daughter_mass, &b_gLLP_daughter_mass);
   fChain->SetBranchAddress("gLLP_grandaughter_id", gLLP_grandaughter_id, &b_gLLP_grandaughter_id);
   fChain->SetBranchAddress("gLLP_grandaughter_pt", gLLP_grandaughter_pt, &b_gLLP_grandaughter_pt);
   fChain->SetBranchAddress("gLLP_grandaughter_eta", gLLP_grandaughter_eta, &b_gLLP_grandaughter_eta);
   fChain->SetBranchAddress("gLLP_grandaughter_phi", gLLP_grandaughter_phi, &b_gLLP_grandaughter_phi);
   fChain->SetBranchAddress("gLLP_grandaughter_eta_ecalcorr", gLLP_grandaughter_eta_ecalcorr, &b_gLLP_grandaughter_eta_ecalcorr);
   fChain->SetBranchAddress("gLLP_grandaughter_phi_ecalcorr", gLLP_grandaughter_phi_ecalcorr, &b_gLLP_grandaughter_phi_ecalcorr);
   fChain->SetBranchAddress("gLLP_grandaughter_e", gLLP_grandaughter_e, &b_gLLP_grandaughter_e);
   fChain->SetBranchAddress("gLLP_grandaughter_mass", gLLP_grandaughter_mass, &b_gLLP_grandaughter_mass);
   Notify();
}

Bool_t llp_event::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void llp_event::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t llp_event::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef llp_event_cxx
