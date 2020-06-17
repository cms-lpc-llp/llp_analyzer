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
   Float_t         muon_validFractionTrackerHits[2000];   //[nMuons]
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
   Float_t         ele_OneOverEminusOneOverP[2000];   //[nElectrons]
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
   Float_t         isoPFCandidateIso04[2000];   //[nIsoPFCandidates]
   Float_t         isoPFCandidateD0[2000];   //[nIsoPFCandidates]
   Int_t           isoPFCandidatePdgId[2000];   //[nIsoPFCandidates]
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
   Int_t           pho_convType[2000];   //[nPhotons]
   Float_t         pho_convTrkZ[2000];   //[nPhotons]
   Float_t         pho_convTrkClusZ[2000];   //[nPhotons]
   Float_t         pho_vtxSumPx[2000][200];   //[nPhotons]
   Float_t         pho_vtxSumPy[2000][200];   //[nPhotons]
   Bool_t          pho_isStandardPhoton[2000];   //[nPhotons]
   Float_t         pho_seedRecHitSwitchToGain6[2000];   //[nPhotons]
   Float_t         pho_seedRecHitSwitchToGain1[2000];   //[nPhotons]
   Float_t         pho_anyRecHitSwitchToGain6[2000];   //[nPhotons]
   Float_t         pho_anyRecHitSwitchToGain1[2000];   //[nPhotons]
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
  //

  Int_t           nCscClusters;
  Float_t         cscClusterX[2000];   //[ncscClusters]
  Float_t         cscClusterY[2000];   //[ncscClusters]
  Float_t         cscClusterZ[2000];   //[ncscClusters]
  Float_t         cscClusterTime[2000];   //[ncscClusters]
  Float_t         cscClusterTimeSpread[2000];   //[ncscClusters]
  Float_t         cscClusterGenMuonDeltaR[2000];   //[ncscClusters]
  Float_t         cscClusterMajorAxis[2000];   //[ncscClusters]
  Float_t         cscClusterMinorAxis[2000];   //[ncscClusters]
  Float_t         cscClusterEtaPhiSpread[2000];   //[ncscClusters]
  Float_t         cscClusterPhiSpread[2000];   //[ncscClusters]
  Float_t         cscClusterEtaSpread[2000];   //[ncscClusters]
  Float_t         cscClusterXSpread[2000];   //[ncscClusters]
  Float_t         cscClusterYSpread[2000];   //[ncscClusters]
  Float_t         cscClusterZSpread[2000];   //[ncscClusters]
  Float_t         cscClusterPhi[2000];   //[ncscClusters]
  Float_t         cscClusterEta[2000];   //[ncscClusters]
  Float_t         cscClusterJetVetoPt[2000];   //[ncscClusters]
  Float_t         cscClusterJetVetoE[2000];   //[ncscClusters]
  Float_t         cscClusterMuonVetoPt[2000];   //[ncscClusters]
  Float_t         cscClusterMuonVetoE[2000];   //[ncscClusters]
  Float_t         cscClusterCaloJetVeto[2000];   //[ncscClusters]
  Int_t           cscClusterSize[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberPlus11[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberPlus12[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberPlus13[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberPlus21[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberPlus22[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberPlus31[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberPlus32[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberPlus41[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberPlus42[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberMinus11[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberMinus12[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberMinus13[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberMinus21[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberMinus22[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberMinus31[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberMinus32[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberMinus41[2000];   //[ncscClusters]
  Int_t           cscClusterNSegmentChamberMinus42[2000];   //[ncscClusters]
  Float_t         cscClusterMe11Ratio[2000];   //[ncscClusters]
  Float_t         cscClusterMe12Ratio[2000];   //[ncscClusters]
  Int_t           cscClusterNStation[2000];   //[ncscClusters]
  Int_t           cscClusterMaxStation[2000];   //[ncscClusters]
  Float_t         cscClusterMaxStationRatio[2000];   //[ncscClusters]
  Int_t           cscClusterNChamber[2000];   //[ncscClusters]
  Int_t           cscClusterMaxChamber[2000];   //[ncscClusters]
  Float_t         cscClusterMaxChamberRatio[2000];   //[ncscClusters]
  Float_t         cscClusterVertexR[2000];   //[ncscClusters]
  Float_t         cscClusterVertexZ[2000];   //[ncscClusters]
  Float_t         cscClusterVertexDis[2000];   //[ncscClusters]
  Float_t         cscClusterVertexChi2[2000];   //[ncscClusters]
  Int_t           cscClusterVertexN1[2000];   //[ncscClusters]
  Int_t           cscClusterVertexN5[2000];   //[ncscClusters]
  Int_t           cscClusterVertexN10[2000];   //[ncscClusters]
  Int_t           cscClusterVertexN15[2000];   //[ncscClusters]
  Int_t           cscClusterVertexN20[2000];   //[ncscClusters]
  Int_t           cscClusterVertexN[2000];   //[ncscClusters]
  Int_t           cscRechitCluster_match_gParticle_id[2000];
  Int_t           dtRechitCluster_match_gParticle_id[2000];
  Int_t           nDtClusters;
  Float_t         dtClusterX[2000];   //[ndtClusters]
  Float_t         dtClusterY[2000];   //[ndtClusters]
  Float_t         dtClusterZ[2000];   //[ndtClusters]
  Float_t         dtClusterTime[2000];   //[ndtClusters]
  Float_t         dtClusterTimeSpread[2000];   //[ndtClusters]
  Float_t         dtClusterGenMuonDeltaR[2000];   //[ndtClusters]
  Float_t         dtClusterMajorAxis[2000];   //[ndtClusters]
  Float_t         dtClusterMinorAxis[2000];   //[ndtClusters]
  Float_t         dtClusterEtaPhiSpread[2000];   //[ndtClusters]
  Float_t         dtClusterPhiSpread[2000];   //[ndtClusters]
  Float_t         dtClusterEtaSpread[2000];   //[ndtClusters]
  Float_t         dtClusterXSpread[2000];   //[ndtClusters]
  Float_t         dtClusterYSpread[2000];   //[ndtClusters]
  Float_t         dtClusterZSpread[2000];   //[ndtClusters]
  Float_t         dtClusterPhi[2000];   //[ndtClusters]
  Float_t         dtClusterEta[2000];   //[ndtClusters]
  Float_t         dtClusterJetVetoPt[2000];   //[ndtClusters]
  Float_t         dtClusterJetVetoE[2000];   //[ndtClusters]
  Float_t         dtClusterMuonVetoPt[2000];   //[ndtClusters]
  Float_t         dtClusterMuonVetoE[2000];   //[ndtClusters]
  Float_t         dtClusterCaloJetVeto[2000];   //[ndtClusters]
  Int_t           dtClusterSize[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberPlus11[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberPlus12[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberPlus13[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberPlus21[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberPlus22[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberPlus31[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberPlus32[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberPlus41[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberPlus42[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberMinus11[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberMinus12[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberMinus13[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberMinus21[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberMinus22[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberMinus31[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberMinus32[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberMinus41[2000];   //[ndtClusters]
  Int_t           dtClusterNSegmentChamberMinus42[2000];   //[ndtClusters]
  Float_t         dtClusterMe11Ratio[2000];   //[ndtClusters]
  Float_t         dtClusterMe12Ratio[2000];   //[ndtClusters]
  Int_t           dtClusterNStation[2000];   //[ndtClusters]
  Int_t           dtClusterMaxStation[2000];   //[ndtClusters]
  Float_t         dtClusterMaxStationRatio[2000];   //[ndtClusters]
  Int_t           dtClusterNChamber[2000];   //[ndtClusters]
  Int_t           dtClusterMaxChamber[2000];   //[ndtClusters]
  Float_t         dtClusterMaxChamberRatio[2000];   //[ndtClusters]
  Float_t         dtClusterVertexR[2000];   //[ndtClusters]
  Float_t         dtClusterVertexZ[2000];   //[ndtClusters]
  Float_t         dtClusterVertexDis[2000];   //[ndtClusters]
  Float_t         dtClusterVertexChi2[2000];   //[ndtClusters]
  Int_t           dtClusterVertexN1[2000];   //[ndtClusters]
  Int_t           dtClusterVertexN5[2000];   //[ndtClusters]
  Int_t           dtClusterVertexN10[2000];   //[ndtClusters]
  Int_t           dtClusterVertexN15[2000];   //[ndtClusters]
  Int_t           dtClusterVertexN20[2000];   //[ndtClusters]
  Int_t           dtClusterVertexN[2000];   //[ndtClusters]
   Int_t           nCscRechitClusters;
   Float_t         cscRechitClusterX[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterY[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterZ[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTime[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterTimeSpread[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterGenMuonDeltaR[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMajorAxis[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMinorAxis[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterEtaPhiSpread[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterPhiSpread[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterEtaSpread[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterXSpread[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterYSpread[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterZSpread[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterPhi[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterEta[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoPt[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterJetVetoE[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoPt[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMuonVetoE[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterCaloJetVeto[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterSize[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus11[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus12[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus13[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus21[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus22[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus31[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus32[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus41[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberPlus42[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus11[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus12[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus13[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus21[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus22[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus31[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus32[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus41[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNRechitChamberMinus42[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMe11Ratio[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMe12Ratio[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNStation[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterMaxStation[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMaxStationRatio[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterNChamber[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterMaxChamber[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterMaxChamberRatio[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterVertexR[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterVertexZ[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterVertexDis[2000];   //[nCscRechitClusters]
   Float_t         cscRechitClusterVertexChi2[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN1[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN5[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN10[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN15[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN20[2000];   //[nCscRechitClusters]
   Int_t           cscRechitClusterVertexN[2000];   //[nCscRechitClusters]

   Int_t           nDtRechitClusters;
   Float_t         dtRechitClusterX[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterY[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterZ[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTime[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterTimeSpread[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterGenMuonDeltaR[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMajorAxis[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMinorAxis[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterEtaPhiSpread[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterPhiSpread[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterEtaSpread[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterXSpread[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterYSpread[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterZSpread[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterPhi[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterEta[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoPt[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterJetVetoE[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoPt[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMuonVetoE[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterCaloJetVeto[2000];   //[nDtRechitClusters]
   Int_t           dtRechitClusterSize[2000];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNStation[2000];   //[nDtRechitClusters]
   Int_t           dtRechitClusterMaxStation[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMaxStationRatio[2000];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNChamber[2000];   //[nDtRechitClusters]
   Int_t           dtRechitClusterMaxChamber[2000];   //[nDtRechitClusters]
   Float_t         dtRechitClusterMaxChamberRatio[2000];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation1[2000];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation2[2000];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation3[2000];   //[nDtRechitClusters]
   Int_t           dtRechitClusterNSegmentStation4[2000];   //[nDtRechitClusters]

   /*Int_t           nCscSegClusters;
   Float_t         cscSegClusterX[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterY[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterZ[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterTime[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterTimeSpread[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterGenMuonDeltaR[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterMajorAxis[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterMinorAxis[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterEtaPhiSpread[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterPhiSpread[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterEtaSpread[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterXSpread[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterYSpread[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterZSpread[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterPhi[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterEta[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterJetVetoPt[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterJetVetoE[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterMuonVetoPt[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterMuonVetoE[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterCaloJetVeto[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterSize[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus11[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus12[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus13[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus21[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus22[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus31[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus32[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus41[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberPlus42[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus11[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus12[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus13[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus21[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus22[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus31[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus32[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus41[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNSegmentChamberMinus42[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterMe11Ratio[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterMe12Ratio[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNStation[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterMaxStation[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterMaxStationRatio[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterNChamber[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterMaxChamber[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterMaxChamberRatio[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterVertexR[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterVertexZ[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterVertexDis[2000];   //[ncscSegClusters]
   Float_t         cscSegClusterVertexChi2[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterVertexN1[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterVertexN5[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterVertexN10[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterVertexN15[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterVertexN20[2000];   //[ncscSegClusters]
   Int_t           cscSegClusterVertexN[2000];   //[ncscSegClusters]
   */

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
   Int_t           rpcBx[2000];   //[nRpc]
   Float_t         rpcT[2000];   //[nRpc]
   Float_t         rpcTError[2000];   //[nRpc]
   Int_t           nHORechits;
   Float_t         hoRechit_Phi[2000];   //[nHORechits]
   Float_t         hoRechit_Eta[2000];   //[nHORechits]
   Float_t         hoRechit_X[2000];   //[nHORechits]
   Float_t         hoRechit_Y[2000];   //[nHORechits]
   Float_t         hoRechit_Z[2000];   //[nHORechits]
   Float_t         hoRechit_T[2000];   //[nHORechits]
   Float_t         hoRechit_E[2000];   //[nHORechits]
   Int_t           nDt;
   Float_t         dtPhi[2000];   //[nDt]
   Float_t         dtEta[2000];   //[nDt]
   Float_t         dtX[2000];   //[nDt]
   Float_t         dtY[2000];   //[nDt]
   Float_t         dtZ[2000];   //[nDt]
   Float_t         dtDirX[2000];   //[nDt]
   Float_t         dtDirY[2000];   //[nDt]
   Float_t         dtDirZ[2000];   //[nDt]
   Float_t         dtNRecHits[2000];   //[nCsc]
   Float_t         dtNRecHits_flag[2000];   //[nCsc]
   Float_t         dtT[2000];   //[nDt]
   Float_t         dtTError[2000];   //[nDt]

   Int_t           nDtRechits;
   Float_t         dtRechitPhi[2000];   //[nDtRechits]
   Float_t         dtRechitEta[2000];   //[nDtRechits]
   Float_t         dtRechitX[2000];   //[nDtRechits]
   Float_t         dtRechitY[2000];   //[nDtRechits]
   Float_t         dtRechitZ[2000];   //[nDtRechits]
   Int_t           dtRechitStation[2000];   //[nDtRechits]
   Int_t           dtRechitWheel[2000];   //[nDtRechits]

   Int_t           nCaloJets;
   Float_t         calojetE[2000];   //[nJets]
   Float_t         calojetPt[2000];   //[nJets]
   Float_t         calojetEta[2000];   //[nJets]
   Float_t         calojetPhi[2000];   //[nJets]
   Int_t           nJets;
   Float_t         jetE[2000];   //[nJets]
   Float_t         jetPt[2000];   //[nJets]
   Float_t         jetEta[2000];   //[nJets]
   Float_t         jetPhi[2000];   //[nJets]
   Float_t         jetCSV[2000];   //[nJets]
   Float_t         jetCISV[2000];   //[nJets]
   Float_t         jetProbb[2000];   //[nJets]
   Float_t         jetProbc[2000];   //[nJets]
   Float_t         jetProbudsg[2000];   //[nJets]
   Float_t         jetProbbb[2000];   //[nJets]
   Float_t         jetMass[2000];   //[nJets]
   Float_t         jetJetArea[2000];   //[nJets]
   Float_t         jetPileupE[2000];   //[nJets]
   Float_t         jetPileupId[2000];   //[nJets]
   Int_t           jetPileupIdFlag[2000];   //[nJets]
   Bool_t          jetPassIDLoose[2000];   //[nJets]
   Bool_t          jetPassIDTight[2000];   //[nJets]
   Bool_t          jetPassMuFrac[2000];   //[nJets]
   Bool_t          jetPassEleFrac[2000];   //[nJets]
   Int_t           jetPartonFlavor[2000];   //[nJets]
   Int_t           jetHadronFlavor[2000];   //[nJets]
   Float_t         jetChargedEMEnergyFraction[2000];   //[nJets]
   Float_t         jetNeutralEMEnergyFraction[2000];   //[nJets]
   Float_t         jetChargedHadronEnergyFraction[2000];   //[nJets]
   Float_t         jetNeutralHadronEnergyFraction[2000];   //[nJets]
   Float_t         jetMuonEnergyFraction[2000];   //[nJets]
   Float_t         jetHOEnergyFraction[2000];   //[nJets]
   Float_t         jetHFHadronEnergyFraction[2000];   //[nJets]
   Float_t         jetHFEMEnergyFraction[2000];   //[nJets]
   Float_t         jetAllMuonPt[2000];   //[nJets]
   Float_t         jetAllMuonEta[2000];   //[nJets]
   Float_t         jetAllMuonPhi[2000];   //[nJets]
   Float_t         jetAllMuonM[2000];   //[nJets]
   Float_t         jetPtWeightedDZ[2000];   //[nJets]
   Int_t           jetNRechits[2000];   //[nJets]
   Float_t         jetRechitE[2000];   //[nJets]
   Float_t         jetRechitT[2000];   //[nJets]
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
   Bool_t          HLTDecision[601];
   Int_t           HLTPrescale[601];
   Int_t           nGenJets;
   Float_t         genJetE[4000];   //[nGenJets]
   Float_t         genJetPt[4000];   //[nGenJets]
   Float_t         genJetEta[4000];   //[nGenJets]
   Float_t         genJetPhi[4000];   //[nGenJets]
   Float_t         genJetMET[4000];   //[nGenJets]

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
   Float_t         gLLP_eta[2];
   Float_t         gLLP_beta[2];
   Float_t         gLLP_phi[2];
   Float_t         gLLP_decay_vertex_x[2];
   Float_t         gLLP_decay_vertex_y[2];
   Float_t         gLLP_decay_vertex_z[2];

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
   TBranch        *b_cscRechitCluster_match_gParticle_id;
   TBranch        *b_dtRechitCluster_match_gParticle_id;
   TBranch        *b_nDtClusters;   //!
   TBranch        *b_dtClusterX;   //!
   TBranch        *b_dtClusterY;   //!
   TBranch        *b_dtClusterZ;   //!
   TBranch        *b_dtClusterTime;   //!
   TBranch        *b_dtClusterTimeSpread;   //!
   TBranch        *b_dtClusterGenMuonDeltaR;   //!
   TBranch        *b_dtClusterMajorAxis;   //!
   TBranch        *b_dtClusterMinorAxis;   //!
   TBranch        *b_dtClusterEtaPhiSpread;   //!
   TBranch        *b_dtClusterPhiSpread;   //!
   TBranch        *b_dtClusterEtaSpread;   //!
   TBranch        *b_dtClusterXSpread;   //!
   TBranch        *b_dtClusterYSpread;   //!
   TBranch        *b_dtClusterZSpread;   //!
   TBranch        *b_dtClusterPhi;   //!
   TBranch        *b_dtClusterEta;   //!
   TBranch        *b_dtClusterJetVetoPt;   //!
   TBranch        *b_dtClusterJetVetoE;   //!
   TBranch        *b_dtClusterMuonVetoPt;   //!
   TBranch        *b_dtClusterMuonVetoE;   //!
   TBranch        *b_dtClusterCaloJetVeto;   //!
   TBranch        *b_dtClusterSize;   //!
   TBranch        *b_dtClusterNSegmentChamberPlus11;   //!
   TBranch        *b_dtClusterNSegmentChamberPlus12;   //!
   TBranch        *b_dtClusterNSegmentChamberPlus13;   //!
   TBranch        *b_dtClusterNSegmentChamberPlus21;   //!
   TBranch        *b_dtClusterNSegmentChamberPlus22;   //!
   TBranch        *b_dtClusterNSegmentChamberPlus31;   //!
   TBranch        *b_dtClusterNSegmentChamberPlus32;   //!
   TBranch        *b_dtClusterNSegmentChamberPlus41;   //!
   TBranch        *b_dtClusterNSegmentChamberPlus42;   //!
   TBranch        *b_dtClusterNSegmentChamberMinus11;   //!
   TBranch        *b_dtClusterNSegmentChamberMinus12;   //!
   TBranch        *b_dtClusterNSegmentChamberMinus13;   //!
   TBranch        *b_dtClusterNSegmentChamberMinus21;   //!
   TBranch        *b_dtClusterNSegmentChamberMinus22;   //!
   TBranch        *b_dtClusterNSegmentChamberMinus31;   //!
   TBranch        *b_dtClusterNSegmentChamberMinus32;   //!
   TBranch        *b_dtClusterNSegmentChamberMinus41;   //!
   TBranch        *b_dtClusterNSegmentChamberMinus42;   //!
   TBranch        *b_dtClusterMe11Ratio;   //!
   TBranch        *b_dtClusterMe12Ratio;   //!
   TBranch        *b_dtClusterNStation;   //!
   TBranch        *b_dtClusterMaxStation;   //!
   TBranch        *b_dtClusterMaxStationRatio;   //!
   TBranch        *b_dtClusterNChamber;   //!
   TBranch        *b_dtClusterMaxChamber;   //!
   TBranch        *b_dtClusterMaxChamberRatio;   //!
   TBranch        *b_dtClusterVertexR;   //!
   TBranch        *b_dtClusterVertexZ;   //!
   TBranch        *b_dtClusterVertexDis;   //!
   TBranch        *b_dtClusterVertexChi2;   //!
   TBranch        *b_dtClusterVertexN1;   //!
   TBranch        *b_dtClusterVertexN5;   //!
   TBranch        *b_dtClusterVertexN10;   //!
   TBranch        *b_dtClusterVertexN15;   //!
   TBranch        *b_dtClusterVertexN20;   //!
   TBranch        *b_dtClusterVertexN;   //!
   
   TBranch        *b_nCscRechitClusters;   //!
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

   TBranch        *b_nDtRechitClusters;   //!
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

   /*
   TBranch        *b_nCscSegClusters;   //!
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
   */   TBranch        *b_nCsc;   //!
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
   TBranch        *b_rpcBx;   //!
   TBranch        *b_rpcTError;   //!
   TBranch        *b_nHORechits;   //!
   TBranch        *b_hoRechit_Phi;   //!
   TBranch        *b_hoRechit_Eta;   //!
   TBranch        *b_hoRechit_X;   //!
   TBranch        *b_hoRechit_Y;   //!
   TBranch        *b_hoRechit_Z;   //!
   TBranch        *b_hoRechit_E;   //!   
   TBranch        *b_hoRechit_T;   //!
   TBranch        *b_nDtRechits;   //!
   TBranch        *b_dtRechitPhi;   //!
   TBranch        *b_dtRechitEta;   //!
   TBranch        *b_dtRechitX;   //!
   TBranch        *b_dtRechitY;   //!
   TBranch        *b_dtRechitZ;   //!
   TBranch        *b_dtRechitStation;   //!   
   TBranch        *b_dtRechitWheel;   //!
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
   TBranch        *b_calojetE;   //!
   TBranch        *b_calojetPt;   //!
   TBranch        *b_calojetEta;   //!
   TBranch        *b_calojetPhi;   //!
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
   TBranch        *b_jetAllMuonPt;   //!
   TBranch        *b_jetAllMuonEta;   //!
   TBranch        *b_jetAllMuonPhi;   //!
   TBranch        *b_jetAllMuonM;   //!
   TBranch        *b_jetPtWeightedDZ;   //!
   TBranch        *b_jetNRechits;   //!
   TBranch        *b_jetRechitE;   //!
   TBranch        *b_jetRechitT;   //!
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
   TBranch        *b_genJetMET;   //!

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
   TBranch        *b_lheComments;
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
   TBranch        *b_gLLP_eta;
   TBranch        *b_gLLP_beta;
   TBranch        *b_gLLP_phi;
   TBranch        *b_gLLP_decay_vertex_x;
   TBranch        *b_gLLP_decay_vertex_y;
   TBranch        *b_gLLP_decay_vertex_z;


   llp_event(TTree *tree=0);
   virtual ~llp_event();
   //virtual Int_t    Cut(Long64_t entry);
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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("llp_ntupler_126.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("llp_ntupler_126.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("llp_ntupler_126.root:/ntuples");
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
  //std::cout << "proper load Entry" << std::endl;
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
   lheComments = 0;
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
   fChain->SetBranchAddress("isoPFCandidatePt", &isoPFCandidatePt, &b_isoPFCandidatePt);
   fChain->SetBranchAddress("isoPFCandidateEta", &isoPFCandidateEta, &b_isoPFCandidateEta);
   fChain->SetBranchAddress("isoPFCandidatePhi", &isoPFCandidatePhi, &b_isoPFCandidatePhi);
   fChain->SetBranchAddress("isoPFCandidateIso04", &isoPFCandidateIso04, &b_isoPFCandidateIso04);
   fChain->SetBranchAddress("isoPFCandidateD0", &isoPFCandidateD0, &b_isoPFCandidateD0);
   fChain->SetBranchAddress("isoPFCandidatePdgId", &isoPFCandidatePdgId, &b_isoPFCandidatePdgId);
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
   /*
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

   fChain->SetBranchAddress("nCscSegClusters", &nCscSegClusters, &b_nCscSegClusters);
   fChain->SetBranchAddress("nDtClusters", &nDtClusters, &b_nDtClusters);
   fChain->SetBranchAddress("dtClusterX", dtClusterX, &b_dtClusterX);
   fChain->SetBranchAddress("dtClusterY", dtClusterY, &b_dtClusterY);
   fChain->SetBranchAddress("dtClusterZ", dtClusterZ, &b_dtClusterZ);
   fChain->SetBranchAddress("dtClusterTime", dtClusterTime, &b_dtClusterTime);
   fChain->SetBranchAddress("dtClusterTimeSpread", dtClusterTimeSpread, &b_dtClusterTimeSpread);
   fChain->SetBranchAddress("dtClusterGenMuonDeltaR", dtClusterGenMuonDeltaR, &b_dtClusterGenMuonDeltaR);
   fChain->SetBranchAddress("dtClusterMajorAxis", dtClusterMajorAxis, &b_dtClusterMajorAxis);
   fChain->SetBranchAddress("dtClusterMinorAxis", dtClusterMinorAxis, &b_dtClusterMinorAxis);
   fChain->SetBranchAddress("dtClusterEtaPhiSpread", dtClusterEtaPhiSpread, &b_dtClusterEtaPhiSpread);
   fChain->SetBranchAddress("dtClusterPhiSpread", dtClusterPhiSpread, &b_dtClusterPhiSpread);
   fChain->SetBranchAddress("dtClusterEtaSpread", dtClusterEtaSpread, &b_dtClusterEtaSpread);
   fChain->SetBranchAddress("dtClusterXSpread", dtClusterXSpread, &b_dtClusterXSpread);
   fChain->SetBranchAddress("dtClusterYSpread", dtClusterYSpread, &b_dtClusterYSpread);
   fChain->SetBranchAddress("dtClusterZSpread", dtClusterZSpread, &b_dtClusterZSpread);
   fChain->SetBranchAddress("dtClusterPhi", dtClusterPhi, &b_dtClusterPhi);
   fChain->SetBranchAddress("dtClusterEta", dtClusterEta, &b_dtClusterEta);
   fChain->SetBranchAddress("dtClusterJetVetoPt", dtClusterJetVetoPt, &b_dtClusterJetVetoPt);
   fChain->SetBranchAddress("dtClusterJetVetoE", dtClusterJetVetoE, &b_dtClusterJetVetoE);
   fChain->SetBranchAddress("dtClusterMuonVetoPt", dtClusterMuonVetoPt, &b_dtClusterMuonVetoPt);
   fChain->SetBranchAddress("dtClusterMuonVetoE", dtClusterMuonVetoE, &b_dtClusterMuonVetoE);
   fChain->SetBranchAddress("dtClusterCaloJetVeto", dtClusterCaloJetVeto, &b_dtClusterCaloJetVeto);
   fChain->SetBranchAddress("dtClusterSize", dtClusterSize, &b_dtClusterSize);
   fChain->SetBranchAddress("dtClusterNSegmentChamberPlus11", dtClusterNSegmentChamberPlus11, &b_dtClusterNSegmentChamberPlus11);
   fChain->SetBranchAddress("dtClusterNSegmentChamberPlus12", dtClusterNSegmentChamberPlus12, &b_dtClusterNSegmentChamberPlus12);
   fChain->SetBranchAddress("dtClusterNSegmentChamberPlus13", dtClusterNSegmentChamberPlus13, &b_dtClusterNSegmentChamberPlus13);
   fChain->SetBranchAddress("dtClusterNSegmentChamberPlus21", dtClusterNSegmentChamberPlus21, &b_dtClusterNSegmentChamberPlus21);
   fChain->SetBranchAddress("dtClusterNSegmentChamberPlus22", dtClusterNSegmentChamberPlus22, &b_dtClusterNSegmentChamberPlus22);
   fChain->SetBranchAddress("dtClusterNSegmentChamberPlus31", dtClusterNSegmentChamberPlus31, &b_dtClusterNSegmentChamberPlus31);
   fChain->SetBranchAddress("dtClusterNSegmentChamberPlus32", dtClusterNSegmentChamberPlus32, &b_dtClusterNSegmentChamberPlus32);
   fChain->SetBranchAddress("dtClusterNSegmentChamberPlus41", dtClusterNSegmentChamberPlus41, &b_dtClusterNSegmentChamberPlus41);
   fChain->SetBranchAddress("dtClusterNSegmentChamberPlus42", dtClusterNSegmentChamberPlus42, &b_dtClusterNSegmentChamberPlus42);
   fChain->SetBranchAddress("dtClusterNSegmentChamberMinus11", dtClusterNSegmentChamberMinus11, &b_dtClusterNSegmentChamberMinus11);
   fChain->SetBranchAddress("dtClusterNSegmentChamberMinus12", dtClusterNSegmentChamberMinus12, &b_dtClusterNSegmentChamberMinus12);
   fChain->SetBranchAddress("dtClusterNSegmentChamberMinus13", dtClusterNSegmentChamberMinus13, &b_dtClusterNSegmentChamberMinus13);
   fChain->SetBranchAddress("dtClusterNSegmentChamberMinus21", dtClusterNSegmentChamberMinus21, &b_dtClusterNSegmentChamberMinus21);
   fChain->SetBranchAddress("dtClusterNSegmentChamberMinus22", dtClusterNSegmentChamberMinus22, &b_dtClusterNSegmentChamberMinus22);
   fChain->SetBranchAddress("dtClusterNSegmentChamberMinus31", dtClusterNSegmentChamberMinus31, &b_dtClusterNSegmentChamberMinus31);
   fChain->SetBranchAddress("dtClusterNSegmentChamberMinus32", dtClusterNSegmentChamberMinus32, &b_dtClusterNSegmentChamberMinus32);
   fChain->SetBranchAddress("dtClusterNSegmentChamberMinus41", dtClusterNSegmentChamberMinus41, &b_dtClusterNSegmentChamberMinus41);
   fChain->SetBranchAddress("dtClusterNSegmentChamberMinus42", dtClusterNSegmentChamberMinus42, &b_dtClusterNSegmentChamberMinus42);
   fChain->SetBranchAddress("dtClusterMe11Ratio", dtClusterMe11Ratio, &b_dtClusterMe11Ratio);
   fChain->SetBranchAddress("dtClusterMe12Ratio", dtClusterMe12Ratio, &b_dtClusterMe12Ratio);
   fChain->SetBranchAddress("dtClusterNStation", dtClusterNStation, &b_dtClusterNStation);
   fChain->SetBranchAddress("dtClusterMaxStation", dtClusterMaxStation, &b_dtClusterMaxStation);
   fChain->SetBranchAddress("dtClusterMaxStationRatio", dtClusterMaxStationRatio, &b_dtClusterMaxStationRatio);
   fChain->SetBranchAddress("dtClusterNChamber", dtClusterNChamber, &b_dtClusterNChamber);
   fChain->SetBranchAddress("dtClusterMaxChamber", dtClusterMaxChamber, &b_dtClusterMaxChamber);
   fChain->SetBranchAddress("dtClusterMaxChamberRatio", dtClusterMaxChamberRatio, &b_dtClusterMaxChamberRatio);
   fChain->SetBranchAddress("dtClusterVertexR", dtClusterVertexR, &b_dtClusterVertexR);
   fChain->SetBranchAddress("dtClusterVertexZ", dtClusterVertexZ, &b_dtClusterVertexZ);
   fChain->SetBranchAddress("dtClusterVertexDis", dtClusterVertexDis, &b_dtClusterVertexDis);
   fChain->SetBranchAddress("dtClusterVertexChi2", dtClusterVertexChi2, &b_dtClusterVertexChi2);
   fChain->SetBranchAddress("dtClusterVertexN1", dtClusterVertexN1, &b_dtClusterVertexN1);
   fChain->SetBranchAddress("dtClusterVertexN5", dtClusterVertexN5, &b_dtClusterVertexN5);
   fChain->SetBranchAddress("dtClusterVertexN10", dtClusterVertexN10, &b_dtClusterVertexN10);
   fChain->SetBranchAddress("dtClusterVertexN15", dtClusterVertexN15, &b_dtClusterVertexN15);
   fChain->SetBranchAddress("dtClusterVertexN20", dtClusterVertexN20, &b_dtClusterVertexN20);
   fChain->SetBranchAddress("dtClusterVertexN", dtClusterVertexN, &b_dtClusterVertexN);
   fChain->SetBranchAddress("nCscSegClusters", &nCscSegClusters, &b_nCscSegClusters);
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
   */
   
   fChain->SetBranchAddress("nCsc", &nCsc, &b_nCsc);
   fChain->SetBranchAddress("cscPhi", cscPhi, &b_cscPhi);
   fChain->SetBranchAddress("cscEta", cscEta, &b_cscEta);
   fChain->SetBranchAddress("cscX", cscX, &b_cscX);
   fChain->SetBranchAddress("cscY", cscY, &b_cscY);
   fChain->SetBranchAddress("cscZ", cscZ, &b_cscZ);
   fChain->SetBranchAddress("cscDirectionX", cscDirectionX, &b_cscDirectionX);
   fChain->SetBranchAddress("cscDirectionY", cscDirectionY, &b_cscDirectionY);
   fChain->SetBranchAddress("cscDirectionZ", cscDirectionZ, &b_cscDirectionZ);

   fChain->SetBranchAddress("nCscRechitClusters", &nCscRechitClusters, &b_nCscRechitClusters);
   fChain->SetBranchAddress("cscRechitCluster_match_gParticle_id", &cscRechitCluster_match_gParticle_id, &b_cscRechitCluster_match_gParticle_id);
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

   fChain->SetBranchAddress("nDtRechitClusters", &nDtRechitClusters, &b_nDtRechitClusters);
   fChain->SetBranchAddress("dtRechitCluster_match_gParticle_id", &dtRechitCluster_match_gParticle_id, &b_dtRechitCluster_match_gParticle_id);
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
   fChain->SetBranchAddress("dtRechitClusterNSegmentStation1", dtRechitClusterNSegmentStation1, &b_dtRechitClusterNSegmentStation1);
   fChain->SetBranchAddress("dtRechitClusterNSegmentStation2", dtRechitClusterNSegmentStation2, &b_dtRechitClusterNSegmentStation2);
   fChain->SetBranchAddress("dtRechitClusterNSegmentStation3", dtRechitClusterNSegmentStation3, &b_dtRechitClusterNSegmentStation3);
   fChain->SetBranchAddress("dtRechitClusterNSegmentStation4", dtRechitClusterNSegmentStation4, &b_dtRechitClusterNSegmentStation4);
   fChain->SetBranchAddress("dtRechitClusterNStation", dtRechitClusterNStation, &b_dtRechitClusterNStation);
   fChain->SetBranchAddress("dtRechitClusterMaxStation", dtRechitClusterMaxStation, &b_dtRechitClusterMaxStation);
   fChain->SetBranchAddress("dtRechitClusterMaxStationRatio", dtRechitClusterMaxStationRatio, &b_dtRechitClusterMaxStationRatio);
   fChain->SetBranchAddress("dtRechitClusterNChamber", dtRechitClusterNChamber, &b_dtRechitClusterNChamber);
   fChain->SetBranchAddress("dtRechitClusterMaxChamber", dtRechitClusterMaxChamber, &b_dtRechitClusterMaxChamber);
   fChain->SetBranchAddress("dtRechitClusterMaxChamberRatio", dtRechitClusterMaxChamberRatio, &b_dtRechitClusterMaxChamberRatio);

   fChain->SetBranchAddress("cscNRecHits", cscNRecHits, &b_cscNRecHits);
   fChain->SetBranchAddress("cscNRecHits_flag", cscNRecHits_flag, &b_cscNRecHits_flag);
   fChain->SetBranchAddress("cscT", cscT, &b_cscT);
   fChain->SetBranchAddress("cscChi2", cscChi2, &b_cscChi2);
   fChain->SetBranchAddress("nRpc", &nRpc, &b_nRpc);
   fChain->SetBranchAddress("rpcPhi", rpcPhi, &b_rpcPhi);
   fChain->SetBranchAddress("rpcEta", rpcEta, &b_rpcEta);
   fChain->SetBranchAddress("rpcX", rpcX, &b_rpcX);
   fChain->SetBranchAddress("rpcY", rpcY, &b_rpcY);
   fChain->SetBranchAddress("rpcZ", rpcZ, &b_rpcZ);
   fChain->SetBranchAddress("rpcT", rpcT, &b_rpcT);
   fChain->SetBranchAddress("rpcBx", rpcBx, &b_rpcBx);
   fChain->SetBranchAddress("rpcTError", rpcTError, &b_rpcTError);
   fChain->SetBranchAddress("nHORechits", &nHORechits, &b_nHORechits);
   fChain->SetBranchAddress("hoRechit_Phi", hoRechit_Phi, &b_hoRechit_Phi);
   fChain->SetBranchAddress("hoRechit_Eta", hoRechit_Eta, &b_hoRechit_Eta);
   fChain->SetBranchAddress("hoRechit_X", hoRechit_X, &b_hoRechit_X);
   fChain->SetBranchAddress("hoRechit_Y", hoRechit_Y, &b_hoRechit_Y);
   fChain->SetBranchAddress("hoRechit_Z", hoRechit_Z, &b_hoRechit_Z);
   fChain->SetBranchAddress("hoRechit_T", hoRechit_T, &b_hoRechit_T);
   fChain->SetBranchAddress("hoRechit_E", hoRechit_E, &b_hoRechit_E);
   fChain->SetBranchAddress("nDtRechits", &nDtRechits, &b_nDtRechits);
   fChain->SetBranchAddress("dtRechitPhi", dtRechitPhi, &b_dtRechitPhi);
   fChain->SetBranchAddress("dtRechitEta", dtRechitEta, &b_dtRechitEta);
   fChain->SetBranchAddress("dtRechitX", dtRechitX, &b_dtRechitX);
   fChain->SetBranchAddress("dtRechitY", dtRechitY, &b_dtRechitY);
   fChain->SetBranchAddress("dtRechitZ", dtRechitZ, &b_dtRechitZ);
   fChain->SetBranchAddress("dtRechitStation", dtRechitStation, &b_dtRechitStation);
   fChain->SetBranchAddress("dtRechitWhee;", dtRechitWheel, &b_dtRechitWheel);
   fChain->SetBranchAddress("nDt", &nDt, &b_nDt);
   fChain->SetBranchAddress("dtPhi", dtPhi, &b_dtPhi);
   fChain->SetBranchAddress("dtEta", dtEta, &b_dtEta);
   fChain->SetBranchAddress("dtX", dtX, &b_dtX);
   fChain->SetBranchAddress("dtY", dtY, &b_dtY);
   fChain->SetBranchAddress("dtZ", dtZ, &b_dtZ);
   fChain->SetBranchAddress("dtDirX", dtDirX, &b_dtDirX);
   fChain->SetBranchAddress("dtDirY", dtDirY, &b_dtDirY);
   fChain->SetBranchAddress("dtDirZ", dtDirZ, &b_dtDirZ);
   fChain->SetBranchAddress("dtT", dtT, &b_dtT);
   fChain->SetBranchAddress("dtTError", dtTError, &b_dtTError);

   fChain->SetBranchAddress("nCaloJets", &nCaloJets, &b_nCaloJets);
   fChain->SetBranchAddress("calojetE", calojetE, &b_calojetE);
   fChain->SetBranchAddress("calojetPt", calojetPt, &b_calojetPt);
   fChain->SetBranchAddress("calojetEta", calojetEta, &b_calojetEta);
   fChain->SetBranchAddress("calojetPhi", calojetPhi, &b_calojetPhi);
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
   fChain->SetBranchAddress("jetAllMuonPt", jetAllMuonPt, &b_jetAllMuonPt);
   fChain->SetBranchAddress("jetAllMuonEta", jetAllMuonEta, &b_jetAllMuonEta);
   fChain->SetBranchAddress("jetAllMuonPhi", jetAllMuonPhi, &b_jetAllMuonPhi);
   fChain->SetBranchAddress("jetAllMuonM", jetAllMuonM, &b_jetAllMuonM);
   fChain->SetBranchAddress("jetPtWeightedDZ", jetPtWeightedDZ, &b_jetPtWeightedDZ);
   fChain->SetBranchAddress("jetNRechits", jetNRechits, &b_jetNRechits);
   fChain->SetBranchAddress("jetRechitE", jetRechitE, &b_jetRechitE);
   fChain->SetBranchAddress("jetRechitT", jetRechitT, &b_jetRechitT);
   fChain->SetBranchAddress("nFatJets", &nFatJets, &b_nFatJets);
   fChain->SetBranchAddress("fatJetE", &fatJetE, &b_fatJetE);
   fChain->SetBranchAddress("fatJetPt", &fatJetPt, &b_fatJetPt);
   fChain->SetBranchAddress("fatJetEta", &fatJetEta, &b_fatJetEta);
   fChain->SetBranchAddress("fatJetPhi", &fatJetPhi, &b_fatJetPhi);
   fChain->SetBranchAddress("fatJetCorrectedPt", &fatJetCorrectedPt, &b_fatJetCorrectedPt);
   fChain->SetBranchAddress("fatJetPrunedM", &fatJetPrunedM, &b_fatJetPrunedM);
   fChain->SetBranchAddress("fatJetTrimmedM", &fatJetTrimmedM, &b_fatJetTrimmedM);
   fChain->SetBranchAddress("fatJetFilteredM", &fatJetFilteredM, &b_fatJetFilteredM);
   fChain->SetBranchAddress("fatJetSoftDropM", &fatJetSoftDropM, &b_fatJetSoftDropM);
   fChain->SetBranchAddress("fatJetCorrectedSoftDropM", &fatJetCorrectedSoftDropM, &b_fatJetCorrectedSoftDropM);
   fChain->SetBranchAddress("fatJetUncorrectedSoftDropM", &fatJetUncorrectedSoftDropM, &b_fatJetUncorrectedSoftDropM);
   fChain->SetBranchAddress("fatJetTau1", &fatJetTau1, &b_fatJetTau1);
   fChain->SetBranchAddress("fatJetTau2", &fatJetTau2, &b_fatJetTau2);
   fChain->SetBranchAddress("fatJetTau3", &fatJetTau3, &b_fatJetTau3);
   fChain->SetBranchAddress("fatJetMaxSubjetCSV", &fatJetMaxSubjetCSV, &b_fatJetMaxSubjetCSV);
   fChain->SetBranchAddress("fatJetPassIDLoose", &fatJetPassIDLoose, &b_fatJetPassIDLoose);
   fChain->SetBranchAddress("fatJetPassIDTight", &fatJetPassIDTight, &b_fatJetPassIDTight);
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
   fChain->SetBranchAddress("genJetMET", genJetMET, &b_genJetMET);

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
   fChain->SetBranchAddress("gLLP_eta", gLLP_eta, &b_gLLP_eta);
   fChain->SetBranchAddress("gLLP_phi", gLLP_phi, &b_gLLP_phi);
   fChain->SetBranchAddress("gLLP_beta", gLLP_beta, &b_gLLP_beta);

   fChain->SetBranchAddress("gLLP_decay_vertex_x", gLLP_decay_vertex_x, &b_gLLP_decay_vertex_x);
   fChain->SetBranchAddress("gLLP_decay_vertex_y", gLLP_decay_vertex_y, &b_gLLP_decay_vertex_y);
   fChain->SetBranchAddress("gLLP_decay_vertex_z", gLLP_decay_vertex_z, &b_gLLP_decay_vertex_z);

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
