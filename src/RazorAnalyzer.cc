#include "RazorAnalyzer.h"
#include "TLorentzVector.h"
#include "Hemisphere.hh"
#include "Davismt2.hh"

using namespace std;

RazorAnalyzer::RazorAnalyzer(TTree *tree) : llp_event(tree)
{
    //turn off all branches
    fChain->SetBranchStatus("*", 0);
}

RazorAnalyzer::~RazorAnalyzer()
{

}

void RazorAnalyzer::Analyze(bool isData, int option, string outputFileName, string label) {
    cout << "Analyze method called on base RazorAnalyzer instance.  Parameters were: " << isData << " " << option << " " << outputFileName << " " << label << endl;
}

//NOTE: the functions below need to be maintained by hand.  If variables are added or removed from the ntuple, these functions need to be updated to reflect the changes.

void RazorAnalyzer::EnableAll(){
    EnableEventInfo();
    // EnablePVAll();
    EnableMuons();
    EnableElectrons();
    // EnableTaus();
    // EnablePhotons();
    EnableCSC();
    EnableJets();
    EnableCaloJets();
    // EnableFatJets();
    EnableMet();
    EnablePileup();
    EnableMC();
    EnableGenParticles();
    EnableLLP();
    // EnableEcalRechits();
}

void RazorAnalyzer::EnableAllWithEcalRechits(){
    EnableEventInfo();
    EnablePVAll();
    EnableMuons();
    EnableElectrons();
    EnableTaus();
    EnablePhotons();
    EnableCSC();
    EnableJets();
    EnableFatJets();
    EnableMet();
    EnablePileup();
    EnableMC();
    EnableGenParticles();
    EnableLLP();
    EnableEcalRechits();
}

void RazorAnalyzer::EnableEventInfo(){
    fChain->SetBranchStatus("nPV", 1);
    fChain->SetBranchStatus("pvX", 1);
    fChain->SetBranchStatus("pvY", 1);
    fChain->SetBranchStatus("pvZ", 1);
    fChain->SetBranchStatus("isData", 1);
    fChain->SetBranchStatus("runNum", 1);
    fChain->SetBranchStatus("lumiNum", 1);
    fChain->SetBranchStatus("eventNum", 1);
    fChain->SetBranchStatus("eventTime", 1);
    fChain->SetBranchStatus("nSlimmedSecondV", 1);
    fChain->SetBranchStatus("fixedGridRhoAll", 1);
    fChain->SetBranchStatus("fixedGridRhoFastjetAll", 1);
    fChain->SetBranchStatus("fixedGridRhoFastjetAllCalo", 1);
    fChain->SetBranchStatus("fixedGridRhoFastjetCentralCalo", 1);
    fChain->SetBranchStatus("fixedGridRhoFastjetCentralChargedPileUp", 1);
    fChain->SetBranchStatus("fixedGridRhoFastjetCentralNeutral", 1);
    fChain->SetBranchStatus("HLTDecision", 1);
    fChain->SetBranchStatus("HLTPrescale", 1);
}

void RazorAnalyzer::EnablePVAll() {
    fChain->SetBranchStatus("nPVAll",1);
    fChain->SetBranchStatus("pvAllX",1);
    fChain->SetBranchStatus("pvAllY",1);
    fChain->SetBranchStatus("pvAllZ",1);
    // fChain->SetBranchStatus("pvAllLogSumPtSq",1);
    // fChain->SetBranchStatus("pvAllSumPx",1);
    // fChain->SetBranchStatus("pvAllSumPy",1);
}

void RazorAnalyzer::EnablePileup(){
    fChain->SetBranchStatus("nBunchXing", 1);
    fChain->SetBranchStatus("BunchXing", 1);
    fChain->SetBranchStatus("nPU", 1);
    fChain->SetBranchStatus("nPUmean", 1);
}

void RazorAnalyzer::EnableMuons(){
    fChain->SetBranchStatus("nMuons", 1);
    fChain->SetBranchStatus("muonE", 1);
    fChain->SetBranchStatus("muonPt", 1);
    fChain->SetBranchStatus("muonEta", 1);
    fChain->SetBranchStatus("muonPhi", 1);
    fChain->SetBranchStatus("muonCharge", 1);
    fChain->SetBranchStatus("muonIsLoose", 1);
    fChain->SetBranchStatus("muonIsMedium", 1);
    fChain->SetBranchStatus("muonIsTight", 1);
    fChain->SetBranchStatus("muon_d0", 1);
    fChain->SetBranchStatus("muon_dZ", 1);
    fChain->SetBranchStatus("muon_ip3d", 1);
    fChain->SetBranchStatus("muon_ip3dSignificance", 1);
    fChain->SetBranchStatus("muonType", 1);
    fChain->SetBranchStatus("muonQuality", 1);
    fChain->SetBranchStatus("muon_pileupIso", 1);
    fChain->SetBranchStatus("muon_chargedIso", 1);
    fChain->SetBranchStatus("muon_photonIso", 1);
    fChain->SetBranchStatus("muon_neutralHadIso", 1);
    fChain->SetBranchStatus("muon_ptrel", 1);
    fChain->SetBranchStatus("muon_chargedMiniIso", 1);
    fChain->SetBranchStatus("muon_photonAndNeutralHadronMiniIso", 1);
    fChain->SetBranchStatus("muon_chargedPileupMiniIso", 1);
    fChain->SetBranchStatus("muon_activityMiniIsoAnnulus", 1);
    fChain->SetBranchStatus("muon_validFractionTrackerHits", 1);
    fChain->SetBranchStatus("muon_isGlobal", 1);
    fChain->SetBranchStatus("muon_normChi2", 1);
    fChain->SetBranchStatus("muon_chi2LocalPosition", 1);
    fChain->SetBranchStatus("muon_kinkFinder", 1);
    fChain->SetBranchStatus("muon_segmentCompatability", 1);
    fChain->SetBranchStatus("muonIsICHEPMedium", 1);
    fChain->SetBranchStatus("muon_passSingleMuTagFilter", 1);
    fChain->SetBranchStatus("muon_passHLTFilter", 1);
}

void RazorAnalyzer::EnableElectrons(){
    fChain->SetBranchStatus("nElectrons", 1);
    fChain->SetBranchStatus("eleE", 1);
    fChain->SetBranchStatus("elePt", 1);
    fChain->SetBranchStatus("eleEta", 1);
    fChain->SetBranchStatus("elePhi", 1);
    fChain->SetBranchStatus("eleCharge", 1);
    fChain->SetBranchStatus("eleEta_SC", 1);
    fChain->SetBranchStatus("eleSigmaIetaIeta", 1);
    fChain->SetBranchStatus("eleFull5x5SigmaIetaIeta", 1);
    fChain->SetBranchStatus("eleR9", 1);
    fChain->SetBranchStatus("ele_dEta", 1);
    fChain->SetBranchStatus("ele_dPhi", 1);
    fChain->SetBranchStatus("ele_HoverE", 1);
    fChain->SetBranchStatus("ele_d0", 1);
    fChain->SetBranchStatus("ele_dZ", 1);
    fChain->SetBranchStatus("ele_ip3d", 1);
    fChain->SetBranchStatus("ele_ip3dSignificance", 1);
    fChain->SetBranchStatus("ele_pileupIso", 1);
    fChain->SetBranchStatus("ele_chargedIso", 1);
    fChain->SetBranchStatus("ele_photonIso", 1);
    fChain->SetBranchStatus("ele_neutralHadIso", 1);
    fChain->SetBranchStatus("ele_MissHits", 1);
    fChain->SetBranchStatus("ele_PassConvVeto", 1);
    fChain->SetBranchStatus("ele_OneOverEminusOneOverP", 1);
    fChain->SetBranchStatus("ele_IDMVAHZZ", 1);
    fChain->SetBranchStatus("ele_IDMVAGeneralPurpose", 1);
    fChain->SetBranchStatus("ele_RegressionE", 1);
    fChain->SetBranchStatus("ele_CombineP4", 1);
    fChain->SetBranchStatus("ele_ptrel", 1);
    fChain->SetBranchStatus("ele_chargedMiniIso", 1);
    fChain->SetBranchStatus("ele_photonAndNeutralHadronMiniIso", 1);
    fChain->SetBranchStatus("ele_chargedPileupMiniIso", 1);
    fChain->SetBranchStatus("ele_activityMiniIsoAnnulus", 1);
    fChain->SetBranchStatus("ele_passSingleEleTagFilter", 1);
    fChain->SetBranchStatus("ele_passTPOneTagFilter", 1);
    fChain->SetBranchStatus("ele_passTPTwoTagFilter", 1);
    fChain->SetBranchStatus("ele_passTPOneProbeFilter", 1);
    fChain->SetBranchStatus("ele_passTPTwoProbeFilter", 1);
    fChain->SetBranchStatus("ele_passHLTFilter", 1);
}

void RazorAnalyzer::EnableTaus(){
    fChain->SetBranchStatus("nTaus", 1);
    fChain->SetBranchStatus("tauE", 1);
    fChain->SetBranchStatus("tauPt", 1);
    fChain->SetBranchStatus("tauEta", 1);
    fChain->SetBranchStatus("tauPhi", 1);
    fChain->SetBranchStatus("tau_IsLoose", 1);
    fChain->SetBranchStatus("tau_IsMedium", 1);
    fChain->SetBranchStatus("tau_IsTight", 1);
    fChain->SetBranchStatus("tau_passEleVetoLoose", 1);
    fChain->SetBranchStatus("tau_passEleVetoMedium", 1);
    fChain->SetBranchStatus("tau_passEleVetoTight", 1);
    fChain->SetBranchStatus("tau_passMuVetoLoose", 1);
    fChain->SetBranchStatus("tau_passMuVetoMedium", 1);
    fChain->SetBranchStatus("tau_passMuVetoTight", 1);
    fChain->SetBranchStatus("tau_ID", 1);
    fChain->SetBranchStatus("tau_combinedIsoDeltaBetaCorr3Hits", 1);
    fChain->SetBranchStatus("tau_chargedIsoPtSum", 1);
    fChain->SetBranchStatus("tau_neutralIsoPtSum", 1);
    fChain->SetBranchStatus("tau_puCorrPtSum", 1);
    fChain->SetBranchStatus("tau_eleVetoMVA", 1);
    fChain->SetBranchStatus("tau_eleVetoCategory", 1);
    fChain->SetBranchStatus("tau_muonVetoMVA", 1);
    fChain->SetBranchStatus("tau_isoMVAnewDMwLT", 1);
    fChain->SetBranchStatus("tau_isoMVAnewDMwoLT", 1);
    fChain->SetBranchStatus("tau_leadCandPt", 1);
    fChain->SetBranchStatus("tau_leadCandID", 1);
    fChain->SetBranchStatus("tau_leadChargedHadrCandPt", 1);
    fChain->SetBranchStatus("tau_leadChargedHadrCandID", 1);
}

void RazorAnalyzer::EnableIsoPFCandidates(){
    // fChain->SetBranchStatus("nIsoPFCandidates", 1);
    // fChain->SetBranchStatus("isoPFCandidatePt", 1);
    // fChain->SetBranchStatus("isoPFCandidateEta", 1);
    // fChain->SetBranchStatus("isoPFCandidatePhi", 1);
    // fChain->SetBranchStatus("isoPFCandidateIso04", 1);
    // fChain->SetBranchStatus("isoPFCandidateD0", 1);
    // fChain->SetBranchStatus("isoPFCandidatePdgId", 1);
}

void RazorAnalyzer::EnablePhotons(){
    fChain->SetBranchStatus("nPhotons", 1);
    fChain->SetBranchStatus("phoE", 1);
    fChain->SetBranchStatus("phoPt", 1);
    fChain->SetBranchStatus("phoEta", 1);
    fChain->SetBranchStatus("phoPhi", 1);
    fChain->SetBranchStatus("phoSigmaIetaIeta", 1);
    fChain->SetBranchStatus("phoFull5x5SigmaIetaIeta", 1);
    fChain->SetBranchStatus("phoR9", 1);
    fChain->SetBranchStatus("pho_HoverE", 1);
    fChain->SetBranchStatus("pho_sumChargedHadronPt", 1);
    fChain->SetBranchStatus("pho_sumNeutralHadronEt", 1);
    fChain->SetBranchStatus("pho_sumPhotonEt", 1);
    fChain->SetBranchStatus("pho_ecalPFClusterIso", 1);
    fChain->SetBranchStatus("pho_hcalPFClusterIso", 1);
    fChain->SetBranchStatus("pho_trkSumPtHollowConeDR03", 1);
    fChain->SetBranchStatus("pho_sumWorstVertexChargedHadronPt", 1);
    fChain->SetBranchStatus("pho_pfIsoChargedHadronIso", 1);
    fChain->SetBranchStatus("pho_pfIsoChargedHadronIsoWrongVtx", 1);
    fChain->SetBranchStatus("pho_pfIsoNeutralHadronIso", 1);
    fChain->SetBranchStatus("pho_pfIsoPhotonIso", 1);
    fChain->SetBranchStatus("pho_pfIsoModFrixione", 1);
    fChain->SetBranchStatus("pho_pfIsoSumPUPt", 1);
    fChain->SetBranchStatus("pho_isConversion", 1);
    fChain->SetBranchStatus("pho_passEleVeto", 1);
    fChain->SetBranchStatus("pho_RegressionE", 1);
    fChain->SetBranchStatus("pho_RegressionEUncertainty", 1);
    fChain->SetBranchStatus("pho_IDMVA", 1);
    fChain->SetBranchStatus("pho_superClusterEnergy", 1);
    fChain->SetBranchStatus("pho_superClusterRawEnergy", 1);
    fChain->SetBranchStatus("pho_superClusterEta", 1);
    fChain->SetBranchStatus("pho_superClusterPhi", 1);
    fChain->SetBranchStatus("pho_superClusterX", 1);
    fChain->SetBranchStatus("pho_superClusterY", 1);
    fChain->SetBranchStatus("pho_superClusterZ", 1);
    fChain->SetBranchStatus("pho_hasPixelSeed", 1);
    fChain->SetBranchStatus("pho_isStandardPhoton", 1);
    fChain->SetBranchStatus("pho_passHLTFilter", 1);
    fChain->SetBranchStatus("pho_convType", 1);
    fChain->SetBranchStatus("pho_convTrkZ", 1);
    fChain->SetBranchStatus("pho_convTrkClusZ", 1);
    fChain->SetBranchStatus("pho_vtxSumPx", 1);
    fChain->SetBranchStatus("pho_vtxSumPy", 1);
    fChain->SetBranchStatus("pho_seedRecHitSwitchToGain6", 1);
    fChain->SetBranchStatus("pho_seedRecHitSwitchToGain1", 1);
    fChain->SetBranchStatus("pho_anyRecHitSwitchToGain6", 1);
    fChain->SetBranchStatus("pho_anyRecHitSwitchToGain1", 1);
};

void RazorAnalyzer::EnableCaloJets(){
  fChain->SetBranchStatus("nCaloJets", 1);
  fChain->SetBranchStatus("calojetE", 1);
  fChain->SetBranchStatus("calojetPt", 1);
  fChain->SetBranchStatus("calojetEta", 1);
  fChain->SetBranchStatus("calojetPhi", 1);
};

void RazorAnalyzer::EnableJets(){
    fChain->SetBranchStatus("nJets", 1);
    fChain->SetBranchStatus("jetE", 1);
    fChain->SetBranchStatus("jetPt", 1);
    fChain->SetBranchStatus("jetEta", 1);
    fChain->SetBranchStatus("jetPhi", 1);
    fChain->SetBranchStatus("jetCSV", 1);
    fChain->SetBranchStatus("jetCISV", 1);
    fChain->SetBranchStatus("jetMass", 1);
    fChain->SetBranchStatus("jetJetArea", 1);
    fChain->SetBranchStatus("jetPileupE", 1);
    fChain->SetBranchStatus("jetPileupId", 1);
    fChain->SetBranchStatus("jetPileupIdFlag", 1);
    fChain->SetBranchStatus("jetPassIDLoose", 1);
    fChain->SetBranchStatus("jetPassIDTight", 1);
    fChain->SetBranchStatus("jetPassMuFrac", 1);
    fChain->SetBranchStatus("jetPassEleFrac", 1);
    fChain->SetBranchStatus("jetPartonFlavor", 1);
    fChain->SetBranchStatus("jetHadronFlavor", 1);
    fChain->SetBranchStatus("jetChargedEMEnergyFraction", 1);
    fChain->SetBranchStatus("jetNeutralEMEnergyFraction", 1);
    fChain->SetBranchStatus("jetChargedHadronEnergyFraction", 1);
    fChain->SetBranchStatus("jetNeutralHadronEnergyFraction", 1);
    fChain->SetBranchStatus("jetMuonEnergyFraction", 1);
    fChain->SetBranchStatus("jetHOEnergyFraction", 1);
    fChain->SetBranchStatus("jetHFHadronEnergyFraction", 1);
    fChain->SetBranchStatus("jetHFEMEnergyFraction", 1);
    fChain->SetBranchStatus("jetAllMuonPt", 1);
    fChain->SetBranchStatus("jetAllMuonEta", 1);
    fChain->SetBranchStatus("jetAllMuonPhi", 1);
    fChain->SetBranchStatus("jetAllMuonM", 1);
    fChain->SetBranchStatus("jetPtWeightedDZ", 1);
    fChain->SetBranchStatus("jetNRechits", 1);
    fChain->SetBranchStatus("jetRechitE", 1);
    fChain->SetBranchStatus("jetRechitT", 1);
};

void RazorAnalyzer::EnableFatJets(){
    fChain->SetBranchStatus("nFatJets", 1);
    fChain->SetBranchStatus("fatJetE", 1);
    fChain->SetBranchStatus("fatJetPt", 1);
    fChain->SetBranchStatus("fatJetEta", 1);
    fChain->SetBranchStatus("fatJetPhi", 1);
    fChain->SetBranchStatus("fatJetCorrectedPt", 1);
    fChain->SetBranchStatus("fatJetPrunedM", 1);
    fChain->SetBranchStatus("fatJetTrimmedM", 1);
    fChain->SetBranchStatus("fatJetFilteredM", 1);
    fChain->SetBranchStatus("fatJetSoftDropM", 1);
    fChain->SetBranchStatus("fatJetCorrectedSoftDropM", 1);
    fChain->SetBranchStatus("fatJetUncorrectedSoftDropM", 1);
    fChain->SetBranchStatus("fatJetTau1", 1);
    fChain->SetBranchStatus("fatJetTau2", 1);
    fChain->SetBranchStatus("fatJetTau3", 1);
    fChain->SetBranchStatus("fatJetMaxSubjetCSV", 1);
    fChain->SetBranchStatus("fatJetPassIDLoose", 1);
    fChain->SetBranchStatus("fatJetPassIDTight", 1);
}

void RazorAnalyzer::EnableMet(){
    fChain->SetBranchStatus("metPt", 1);
    fChain->SetBranchStatus("metPhi", 1);
    fChain->SetBranchStatus("metNoHFPt", 1);
    fChain->SetBranchStatus("metNoHFPhi", 1);
    fChain->SetBranchStatus("metType0Pt", 1);
    fChain->SetBranchStatus("metType0Phi", 1);
    fChain->SetBranchStatus("metType1Pt", 1);
    fChain->SetBranchStatus("metType1Phi", 1);
    //fChain->SetBranchStatus("metEGCleanPt", 1);
    //fChain->SetBranchStatus("metEGCleanPhi", 1);
    //fChain->SetBranchStatus("metMuEGCleanPt", 1);
    //fChain->SetBranchStatus("metMuEGCleanPhi", 1);
    //fChain->SetBranchStatus("metMuEGCleanCorrPt", 1);
    //fChain->SetBranchStatus("metMuEGCleanCorrPhi", 1);
    //fChain->SetBranchStatus("metUncorrectedPt", 1);
    //fChain->SetBranchStatus("metUncorrectedPhi", 1);
    fChain->SetBranchStatus("metType0Plus1Pt", 1);
    fChain->SetBranchStatus("metType0Plus1Phi", 1);
    fChain->SetBranchStatus("metPuppiPt", 1);
    fChain->SetBranchStatus("metPuppiPhi", 1);
    fChain->SetBranchStatus("metCaloPt", 1);
    fChain->SetBranchStatus("metCaloPhi", 1);
    fChain->SetBranchStatus("sumMET", 1);
    fChain->SetBranchStatus("Flag_HBHENoiseFilter", 1);
    fChain->SetBranchStatus("Flag_HBHETightNoiseFilter", 1);
    fChain->SetBranchStatus("Flag_HBHEIsoNoiseFilter", 1);
    fChain->SetBranchStatus("Flag_badChargedCandidateFilter", 1);
    fChain->SetBranchStatus("Flag_BadChargedCandidateFilter", 1);
    fChain->SetBranchStatus("Flag_badMuonFilter", 1);
    fChain->SetBranchStatus("Flag_BadPFMuonFilter", 1);
    fChain->SetBranchStatus("Flag_badGlobalMuonFilter", 1);
    fChain->SetBranchStatus("Flag_duplicateMuonFilter", 1);
    fChain->SetBranchStatus("Flag_CSCTightHaloFilter", 1);
    fChain->SetBranchStatus("Flag_hcalLaserEventFilter", 1);
    fChain->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1);
    fChain->SetBranchStatus("Flag_EcalDeadCellBoundaryEnergyFilter", 1);
    fChain->SetBranchStatus("Flag_goodVertices", 1);
    fChain->SetBranchStatus("Flag_trackingFailureFilter", 1);
    fChain->SetBranchStatus("Flag_eeBadScFilter", 1);
    fChain->SetBranchStatus("Flag_ecalLaserCorrFilter", 1);
    fChain->SetBranchStatus("Flag_trkPOGFilters", 1);
    fChain->SetBranchStatus("Flag_trkPOG_manystripclus53X", 1);
    fChain->SetBranchStatus("Flag_trkPOG_toomanystripclus53X", 1);
    fChain->SetBranchStatus("Flag_trkPOG_logErrorTooManyClusters", 1);
    fChain->SetBranchStatus("Flag_METFilters", 1);
    fChain->SetBranchStatus("Flag_ecalBadCalibFilter", 1);

    // fChain->SetBranchStatus("metType1PtJetResUp", 1);
    // fChain->SetBranchStatus("metType1PtJetResDown", 1);
    // fChain->SetBranchStatus("metType1PtJetEnUp", 1);
    // fChain->SetBranchStatus("metType1PtJetEnDown", 1);
    // fChain->SetBranchStatus("metType1PtMuonEnUp", 1);
    // fChain->SetBranchStatus("metType1PtMuonEnDown", 1);
    // fChain->SetBranchStatus("metType1PtElectronEnUp", 1);
    // fChain->SetBranchStatus("metType1PtElectronEnDown", 1);
    // fChain->SetBranchStatus("metType1PtTauEnUp", 1);
    // fChain->SetBranchStatus("metType1PtTauEnDown", 1);
    // fChain->SetBranchStatus("metType1PtUnclusteredEnUp", 1);
    // fChain->SetBranchStatus("metType1PtUnclusteredEnDown", 1);
    // fChain->SetBranchStatus("metType1PtPhotonEnUp", 1);
    // fChain->SetBranchStatus("metType1PtPhotonEnDown", 1);
    // fChain->SetBranchStatus("metType1PhiJetResUp", 1);
    // fChain->SetBranchStatus("metType1PhiJetResDown", 1);
    // fChain->SetBranchStatus("metType1PhiJetEnUp", 1);
    // fChain->SetBranchStatus("metType1PhiJetEnDown", 1);
    // fChain->SetBranchStatus("metType1PhiMuonEnUp", 1);
    // fChain->SetBranchStatus("metType1PhiMuonEnDown", 1);
    // fChain->SetBranchStatus("metType1PhiElectronEnUp", 1);
    // fChain->SetBranchStatus("metType1PhiElectronEnDown", 1);
    // fChain->SetBranchStatus("metType1PhiTauEnUp", 1);
    // fChain->SetBranchStatus("metType1PhiTauEnDown", 1);
    // fChain->SetBranchStatus("metType1PhiUnclusteredEnUp", 1);
    // fChain->SetBranchStatus("metType1PhiUnclusteredEnDown", 1);
    // fChain->SetBranchStatus("metType1PhiPhotonEnUp", 1);
    // fChain->SetBranchStatus("metType1PhiPhotonEnDown", 1);
}

void RazorAnalyzer::EnableRazor()
{

};

void RazorAnalyzer::EnableCSC()
{
    fChain->SetBranchStatus("nCsc", 1);
    fChain->SetBranchStatus("cscPhi", 1);
    fChain->SetBranchStatus("cscEta", 1);
    fChain->SetBranchStatus("cscX", 1);
    fChain->SetBranchStatus("cscY", 1);
    fChain->SetBranchStatus("cscZ", 1);
    fChain->SetBranchStatus("cscDirectionX", 1);
    fChain->SetBranchStatus("cscDirectionY", 1);
    fChain->SetBranchStatus("cscDirectionZ", 1);
    fChain->SetBranchStatus("cscNRecHits", 1);
    fChain->SetBranchStatus("cscNRecHits_flag", 1);
    fChain->SetBranchStatus("cscT", 1);
    fChain->SetBranchStatus("cscChi2", 1);
};

void RazorAnalyzer::EnableMC(){
    fChain->SetBranchStatus("nGenJets", 1);
    fChain->SetBranchStatus("genJetE", 1);
    fChain->SetBranchStatus("genJetPt", 1);
    fChain->SetBranchStatus("genJetEta", 1);
    fChain->SetBranchStatus("genJetPhi", 1);
    fChain->SetBranchStatus("genMetPtCalo", 1);
    fChain->SetBranchStatus("genMetPhiCalo", 1);
    fChain->SetBranchStatus("genMetPtTrue", 1);
    fChain->SetBranchStatus("genMetPhiTrue", 1);
    fChain->SetBranchStatus("genVertexX", 1);
    fChain->SetBranchStatus("genVertexY", 1);
    fChain->SetBranchStatus("genVertexZ", 1);
    fChain->SetBranchStatus("genVertexT", 1);
    fChain->SetBranchStatus("genWeight", 1);
    fChain->SetBranchStatus("genSignalProcessID", 1);
    fChain->SetBranchStatus("genQScale", 1);
    fChain->SetBranchStatus("genAlphaQCD", 1);
    fChain->SetBranchStatus("genAlphaQED", 1);
    //fChain->SetBranchStatus("lheComments", 1);
    fChain->SetBranchStatus("scaleWeights", 1);
    fChain->SetBranchStatus("pdfWeights", 1);
    fChain->SetBranchStatus("alphasWeights", 1);
}

void RazorAnalyzer::EnableGenParticles(){
    fChain->SetBranchStatus("nGenParticle", 1);
    fChain->SetBranchStatus("gParticleMotherId", 1);
    fChain->SetBranchStatus("gParticleMotherIndex", 1);
    fChain->SetBranchStatus("gParticleId", 1);
    fChain->SetBranchStatus("gParticleStatus", 1);
    fChain->SetBranchStatus("gParticleE", 1);
    fChain->SetBranchStatus("gParticlePt", 1);
    fChain->SetBranchStatus("gParticleEta", 1);
    fChain->SetBranchStatus("gParticlePhi", 1);
    //fChain->SetBranchStatus("gParticleDecayVertexX", 1);
    //fChain->SetBranchStatus("gParticleDecayVertexY", 1);
    //fChain->SetBranchStatus("gParticleDecayVertexZ", 1);
}
void RazorAnalyzer::EnableLLP(){
    fChain->SetBranchStatus("gLLP_eta", 1);
    fChain->SetBranchStatus("gLLP_phi", 1);
    fChain->SetBranchStatus("gLLP_beta", 1);

    fChain->SetBranchStatus("gLLP_decay_vertex_x", 1);
    fChain->SetBranchStatus("gLLP_decay_vertex_y", 1);
    fChain->SetBranchStatus("gLLP_decay_vertex_z", 1);
    //fChain->SetBranchStatus("gParticleDecayVertexX", 1);
    //fChain->SetBranchStatus("gParticleDecayVertexY", 1);
    //fChain->SetBranchStatus("gParticleDecayVertexZ", 1);
}

void RazorAnalyzer::EnableEcalRechits(){
    fChain->SetBranchStatus("ecalRechit_Eta", 1);
    fChain->SetBranchStatus("ecalRechit_Phi", 1);
    fChain->SetBranchStatus("ecalRechit_X", 1);
    fChain->SetBranchStatus("ecalRechit_Y", 1);
    fChain->SetBranchStatus("ecalRechit_Z", 1);
    fChain->SetBranchStatus("ecalRechit_E", 1);
    fChain->SetBranchStatus("ecalRechit_T", 1);
    fChain->SetBranchStatus("ecalRechit_ID", 1);
    fChain->SetBranchStatus("ecalRechit_FlagOOT", 1);
    fChain->SetBranchStatus("ecalRechit_GainSwitch1", 1);
    fChain->SetBranchStatus("ecalRechit_GainSwitch6", 1);
    fChain->SetBranchStatus("ecalRechit_transpCorr", 1);
}


//////////////////////////////
//ELECTRON
//////////////////////////////

float RazorAnalyzer::GetElectronScaleCorrection( double pt, double eta ) {
  double scaleCorr = 1.0;
  // if ( pt > 0 && fabs(eta) < 1.5) {
  //   scaleCorr = 1.015;
  // } else {
  //   scaleCorr = 1.05;
  // }
  //some dummy code so that it compiler doesn't complain
  if (pt > 0 || eta > 0) {
    scaleCorr = 1.0;
  }
  return scaleCorr;
}

float RazorAnalyzer::GetElectronEffectiveAreaMean(int i, bool use25nsCuts ){

    double effArea = 0.0;
    //Effective areas below are for the sum of Neutral Hadrons + Photons
    if (use25nsCuts) {
        // These are the Spring15 25ns effective areas reported in this presentation:
        // https://indico.cern.ch/event/369239/contributions/874575/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
        if (fabs(eleEta_SC[i]) < 1.0) {
            effArea = 0.1752;
        } else if (fabs(eleEta_SC[i]) < 1.479) {
            effArea = 0.1862;
        } else if (fabs(eleEta_SC[i]) < 2.0) {
            effArea = 0.1411;
        } else if (fabs(eleEta_SC[i]) < 2.2) {
            effArea = 0.1534;
        } else if (fabs(eleEta_SC[i]) < 2.3) {
            effArea = 0.1903;
        } else if (fabs(eleEta_SC[i]) < 2.4) {
            effArea = 0.2243;
        } else if (fabs(eleEta_SC[i]) < 2.5) {
            effArea = 0.2687;
        }
        return effArea;
    }
    else {
        // These are the Spring15 50ns effective areas reported in this presentation:
        // https://indico.cern.ch/event/369235/contributions/874560/attachments/734635/1007867/Rami_EffAreas.pdf
        if (fabs(eleEta_SC[i]) < 0.8) {
            effArea = 0.0973;
        } else if (fabs(eleEta_SC[i]) < 1.3) {
            effArea = 0.0954;
        } else if (fabs(eleEta_SC[i]) < 2.0) {
            effArea = 0.0632;
        } else if (fabs(eleEta_SC[i]) < 2.2) {
            effArea = 0.0727;
        } else {
            effArea = 0.1337;
        }
    }
    return effArea;
}

float RazorAnalyzer::GetElectronEffectiveArea90(int i, string EraName){

  double effArea = 0.0;
  //Effective areas derived on Spring15 MC
  // We will keep using these old ones for Moriond 2017 analyses
  //Reference: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF
  if (EraName == "Spring15") {
    if (fabs(eleEta_SC[i]) < 1.0) {
      effArea = 0.1752;
    } else if (fabs(eleEta_SC[i]) < 1.479) {
      effArea = 0.1862;
    } else if (fabs(eleEta_SC[i]) < 2.0) {
      effArea = 0.1411;
    } else if (fabs(eleEta_SC[i]) < 2.2) {
      effArea = 0.1534;
    } else if (fabs(eleEta_SC[i]) < 2.3) {
      effArea = 0.1903;
    } else if (fabs(eleEta_SC[i]) < 2.4) {
      effArea = 0.2243;
    } else if (fabs(eleEta_SC[i]) < 2.5) {
      effArea = 0.2687;
    }
  }
  else if (EraName == "Spring16" || EraName == "Summer16") {
    //New effective areas derived from Spring 16 MC
    //Reference: https://indico.cern.ch/event/482673/contributions/2187022/attachments/1282446/1905912/talk_electron_ID_spring16.pdf
    //Effective areas below are for the sum of Neutral Hadrons + Photons
    if (fabs(eleEta_SC[i]) < 1.0) {
      effArea = 0.1703;
    } else if (fabs(eleEta_SC[i]) < 1.479) {
      effArea = 0.1715;
    } else if (fabs(eleEta_SC[i]) < 2.0) {
      effArea = 0.1213;
    } else if (fabs(eleEta_SC[i]) < 2.2) {
      effArea = 0.1230;
    } else if (fabs(eleEta_SC[i]) < 2.3) {
      effArea = 0.1635;
    } else if (fabs(eleEta_SC[i]) < 2.4) {
      effArea = 0.1937;
    } else if (fabs(eleEta_SC[i]) < 2.5) {
      effArea = 0.2393;
    }
  }
  else if (EraName == "2017_94X") {
    //Effective areas derived from 94X samples
    //Reference numbers: https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
    //Details found here: https://indico.cern.ch/event/697576/contributions/2940576/attachments/1620927/2578913/eleIdTuning.pdf
    if (fabs(eleEta_SC[i]) < 1.0) {
      effArea = 0.1440;
    } else if (fabs(eleEta_SC[i]) < 1.479) {
      effArea = 0.1562;
    } else if (fabs(eleEta_SC[i]) < 2.0) {
      effArea = 0.1032;
    } else if (fabs(eleEta_SC[i]) < 2.2) {
      effArea = 0.0859;
    } else if (fabs(eleEta_SC[i]) < 2.3) {
      effArea = 0.1116;
    } else if (fabs(eleEta_SC[i]) < 2.4) {
      effArea = 0.1321;
    } else if (fabs(eleEta_SC[i]) < 2.5) {
      effArea = 0.1654;
    }
  }

  return effArea;
}

bool RazorAnalyzer::isEGammaPOGVetoElectron(int i, bool applyID, bool applyIso, bool use25nsCuts, string EraName){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGVetoElectronID(i,use25nsCuts, EraName)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGVetoElectronIso(i,use25nsCuts)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isEGammaPOGLooseElectron(int i, bool applyID, bool applyIso, bool use25nsCuts, string EraName){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGLooseElectronID(i,use25nsCuts, EraName)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGLooseElectronIso(i,use25nsCuts)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isEGammaPOGMediumElectron(int i, bool applyID, bool applyIso, bool use25nsCuts, string EraName){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGMediumElectronID(i,use25nsCuts, EraName)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGMediumElectronIso(i,use25nsCuts)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isEGammaPOGTightElectron(int i, bool applyID, bool applyIso, bool use25nsCuts, string EraName){
  bool pass = true;
  if (applyID) {
    if (!passEGammaPOGTightElectronID(i,use25nsCuts, EraName)) pass = false;
  }
  if (applyIso) {
    if (!passEGammaPOGTightElectronIso(i,use25nsCuts)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isVetoElectron(int i, bool applyID, bool applyIso, string EraName){
  return isMVANonTrigVetoElectron(i, applyID, applyIso, EraName);
}

bool RazorAnalyzer::isLooseElectron(int i, bool applyID, bool applyIso, bool use25nsCuts, string EraName){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
  if (applyID) {
    if (!passEGammaPOGLooseElectronID(i,use25nsCuts, EraName)) pass = false;
  }
  if (applyIso) {
    if (!((ele_chargedMiniIso[i] + fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveAreaMean(i)*pow(dr/0.3,2)))/elePt[i] < 0.1)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isMediumElectron(int i, bool applyID, bool applyIso, bool use25nsCuts, string EraName){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
  if (applyID) {
    if (!passEGammaPOGMediumElectronID(i,use25nsCuts, EraName)) pass = false;
  }
  if (applyIso) {
    if (!((ele_chargedMiniIso[i] +  fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveAreaMean(i)*pow(dr/0.3,2)))/elePt[i] < 0.1)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isTightElectron(int i, bool applyID, bool applyIso, bool use25nsCuts, string EraName){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
  if (applyID) {
    if (!passEGammaPOGTightElectronID(i,use25nsCuts, EraName)) pass = false;
  }
  if (applyIso) {
    if (!((ele_chargedMiniIso[i] + fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveAreaMean(i)*pow(dr/0.3,2)))/elePt[i] < 0.1)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isMVANonTrigVetoElectron(int i, bool applyID, bool applyIso, string EraName){

  bool pass = true;
  if (applyID) {
    if (!passMVANonTrigVetoElectronID(i,EraName)) pass = false;
  }
  if (applyIso) {
    if (!passMVANonTrigVetoElectronIso(i)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::passEGammaPOGVetoElectronID(int i, bool use25nsCuts, string EraName){
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;

    if (EraName == "Spring15") {
      //For Moriond 2017, SUSY group decided to keep using Spring15 cut-based ID
      if(fabs(eleEta_SC[i]) < 1.479) {
	if ( fabs(ele_dEta[i]) < 0.0152
	     && fabs(ele_dPhi[i]) < 0.216
	     && eleFull5x5SigmaIetaIeta[i] < 0.0114
	     && ele_HoverE[i] < 0.181
	     && fabs(ele_d0[i]) < 0.0564
	     && fabs(ele_dZ[i]) < 0.472
	     && fabs(ele_OneOverEminusOneOverP[i]) < 0.207
	     && ele_PassConvVeto[i]
	     && ele_MissHits[i] <= 2
	     ) {
	  pass = true;
	}
      } else {
	if (fabs(ele_dEta[i]) < 0.0113
	    && fabs(ele_dPhi[i]) < 0.237
	    && eleFull5x5SigmaIetaIeta[i] < 0.0352
	    && ele_HoverE[i] < 0.116
	    && fabs(ele_d0[i]) < 0.222
	    && fabs(ele_dZ[i]) < 0.921
	    && fabs(ele_OneOverEminusOneOverP[i]) < 0.174
	    && ele_PassConvVeto[i]
	    && ele_MissHits[i] <= 3
	    ) {
	  pass = true;
	}
      }
    }

    else if (EraName == "Summer16") {

      // Veto ID recommended for analyses performed on 2016 data using 8XX releases.
      // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
      // Based on these slides: https://indico.cern.ch/event/732971/contributions/3022843/attachments/1658685/2656462/eleIdTuning.pdf

      if(fabs(eleEta_SC[i]) < 1.479) {
        if ( fabs(ele_dEta[i]) < 0.00749
	     && fabs(ele_dPhi[i]) < 0.228
	     && eleFull5x5SigmaIetaIeta[i] < 0.0115
	     && ele_HoverE[i] < 0.356
	     && fabs(ele_d0[i]) < 0.05
	     && fabs(ele_dZ[i]) < 0.10
	     && fabs(ele_OneOverEminusOneOverP[i]) < 0.299
	     && ele_PassConvVeto[i]
	     && ele_MissHits[i] <= 2
	     ) {
	  pass = true;
        }
      } else {
        if (fabs(ele_dEta[i]) < 0.00895
	    && fabs(ele_dPhi[i]) < 0.213
	    && eleFull5x5SigmaIetaIeta[i] < 0.037
	    && ele_HoverE[i] < 0.211
	    && fabs(ele_d0[i]) < 0.1
	    && fabs(ele_dZ[i]) < 0.2
	    && fabs(ele_OneOverEminusOneOverP[i]) < 0.15
	    && ele_PassConvVeto[i]
	    && ele_MissHits[i] <= 3
	    ) {
	  pass = true;
        }
      }
    }
    else if (EraName == "2017_94X") {

      // Veto ID recommended for analyses performed on 2017 data using 94X releases.
      // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
      // Based on these slides: https://indico.cern.ch/event/732971/contributions/3022843/attachments/1658685/2656462/eleIdTuning.pdf

      if(fabs(eleEta_SC[i]) < 1.479) {
	if ( fabs(ele_dEta[i]) < 0.00463
	     && fabs(ele_dPhi[i]) < 0.148
	     && eleFull5x5SigmaIetaIeta[i] < 0.0126
	     && ele_HoverE[i] < 0.05 + 1.16 / eleE[i] + 0.0324*fixedGridRhoFastjetAll / eleE[i]
	     && fabs(ele_d0[i]) < 0.05
	     && fabs(ele_dZ[i]) < 0.10
	     && fabs(ele_OneOverEminusOneOverP[i]) < 0.209
	     && ele_PassConvVeto[i]
	     && ele_MissHits[i] <= 2
	     ) {
	  pass = true;
	}
      } else {
	if (fabs(ele_dEta[i]) < 0.00814
	    && fabs(ele_dPhi[i]) < 0.19
	    && eleFull5x5SigmaIetaIeta[i] < 0.0457
	    && ele_HoverE[i] < 0.05 + 2.54 / eleE[i] + 0.183*fixedGridRhoFastjetAll / eleE[i]
	    && fabs(ele_d0[i]) < 0.1
	    && fabs(ele_dZ[i]) < 0.2
	    && fabs(ele_OneOverEminusOneOverP[i]) < 0.132
	    && ele_PassConvVeto[i]
	    && ele_MissHits[i] <= 3
	    ) {
	  pass = true;
	}
      }
    }
    return pass;
}

bool RazorAnalyzer::passEGammaPOGLooseElectronID(int i, bool use25nsCuts, string EraName){
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;

    if (EraName == "Spring15") {
      //For Moriond 2017, SUSY group decided to keep using Spring15 cut-based ID
      if(fabs(eleEta_SC[i]) < 1.479) {
	if ( fabs(ele_dEta[i]) < 0.0105
	     && fabs(ele_dPhi[i]) < 0.115
	     && eleFull5x5SigmaIetaIeta[i] < 0.0103
	     && ele_HoverE[i] < 0.104
	     && fabs(ele_d0[i]) < 0.0261
	     && fabs(ele_dZ[i]) < 0.41
	     && fabs(ele_OneOverEminusOneOverP[i]) < 0.102
	     && ele_PassConvVeto[i]
	     && ele_MissHits[i] <= 2
	     ) {
	  pass = true;
	}
      } else {
	if (fabs(ele_dEta[i]) < 0.00814
	    && fabs(ele_dPhi[i]) < 0.182
	    && eleFull5x5SigmaIetaIeta[i] < 0.0301
	    && ele_HoverE[i] < 0.0897
	    && fabs(ele_d0[i]) < 0.118
	    && fabs(ele_dZ[i]) < 0.822
	    && fabs(ele_OneOverEminusOneOverP[i]) < 0.126
	    && ele_PassConvVeto[i]
	    && ele_MissHits[i] <= 1
	    ) {
	  pass = true;
	}
      }
    }
    else if (EraName == "Summer16") {

      // Loose ID recommended for analyses performed on 2016 data using 8XX releases.
      // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
      // Based on these slides: https://indico.cern.ch/event/482677/contributions/2259342/attachments/1316731/1972911/talk_electron_ID_spring16_update.pdf
      if(fabs(eleEta_SC[i]) < 1.479) {
	if ( fabs(ele_dEta[i]) < 0.00477
	     && fabs(ele_dPhi[i]) < 0.222
	     && eleFull5x5SigmaIetaIeta[i] < 0.011
	     && ele_HoverE[i] < 0.298
	     && fabs(ele_d0[i]) < 0.05
	     && fabs(ele_dZ[i]) < 0.10
	     && fabs(ele_OneOverEminusOneOverP[i]) < 0.241
	     && ele_PassConvVeto[i]
	     && ele_MissHits[i] <= 1
	     ) {
	  pass = true;
	}
      } else {
	if (fabs(ele_dEta[i]) < 0.00868
	    && fabs(ele_dPhi[i]) < 0.213
	    && eleFull5x5SigmaIetaIeta[i] < 0.0314
	    && ele_HoverE[i] < 0.101
	    && fabs(ele_d0[i]) < 0.1
	    && fabs(ele_dZ[i]) < 0.2
	    && fabs(ele_OneOverEminusOneOverP[i]) < 0.14
	    && ele_PassConvVeto[i]
	    && ele_MissHits[i] <= 1
	    ) {
	  pass = true;
	}
      }
    }
    else if (EraName == "2017_94X") {

      // Loose ID recommended for analyses performed on 2017 data using 94X releases.
      // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
      // Based on these slides: https://indico.cern.ch/event/732971/contributions/3022843/attachments/1658685/2656462/eleIdTuning.pdf

      if(fabs(eleEta_SC[i]) < 1.479) {
	if ( fabs(ele_dEta[i]) < 0.00377
	     && fabs(ele_dPhi[i]) < 0.0884
	     && eleFull5x5SigmaIetaIeta[i] < 0.0112
	     && ele_HoverE[i] < 0.05 + 1.16 / eleE[i] + 0.0324*fixedGridRhoFastjetAll / eleE[i]
	     && fabs(ele_d0[i]) < 0.05
	     && fabs(ele_dZ[i]) < 0.10
	     && fabs(ele_OneOverEminusOneOverP[i]) < 0.193
	     && ele_PassConvVeto[i]
	     && ele_MissHits[i] <= 1
	     ) {
	  pass = true;
	}
      } else {
	if (fabs(ele_dEta[i]) < 0.00674
	    && fabs(ele_dPhi[i]) < 0.169
	    && eleFull5x5SigmaIetaIeta[i] < 0.0425
	    && ele_HoverE[i] < 0.0441 + 2.54 / eleE[i] + 0.183*fixedGridRhoFastjetAll / eleE[i]
	    && fabs(ele_d0[i]) < 0.1
	    && fabs(ele_dZ[i]) < 0.2
	    && fabs(ele_OneOverEminusOneOverP[i]) < 0.111
	    && ele_PassConvVeto[i]
	    && ele_MissHits[i] <= 1
	    ) {
	  pass = true;
	}
      }
    }

    return pass;
}

bool RazorAnalyzer::passEGammaPOGMediumElectronID(int i, bool use25nsCuts, string EraName){
  if (!use25nsCuts) {
    std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
    return false;
  }
  bool pass = false;

  if (EraName == "Spring15") {
    //For Moriond 2017, SUSY group decided to keep using Spring15 cut-based ID
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.0103
	   && fabs(ele_dPhi[i]) < 0.0336
	   && eleFull5x5SigmaIetaIeta[i] < 0.0101
	   && ele_HoverE[i] < 0.0876
	   && fabs(ele_d0[i]) < 0.0118
	   && fabs(ele_dZ[i]) < 0.373
	   && fabs(ele_OneOverEminusOneOverP[i]) < 0.0174
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.00733
	  && fabs(ele_dPhi[i]) < 0.114
	  && eleFull5x5SigmaIetaIeta[i] < 0.0283
	  && ele_HoverE[i] < 0.0678
	  && fabs(ele_d0[i]) < 0.0739
	  && fabs(ele_dZ[i]) < 0.602
	  && fabs(ele_OneOverEminusOneOverP[i]) < 0.0898
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 1
	  ) {
	pass = true;
      }
    }
  }
  else if (EraName == "Summer16") {
    // Medium ID recommended for analyses performed on 2016 data using 8XX releases.
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
    // Based on these slides: https://indico.cern.ch/event/482677/contributions/2259342/attachments/1316731/1972911/talk_electron_ID_spring16_update.pdf

    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.00311
    	   && fabs(ele_dPhi[i]) < 0.103
    	   && eleFull5x5SigmaIetaIeta[i] < 0.00998
    	   && ele_HoverE[i] < 0.253
    	   && fabs(ele_d0[i]) < 0.05
    	   && fabs(ele_dZ[i]) < 0.10
    	   && fabs(ele_OneOverEminusOneOverP[i]) < 0.134
    	   && ele_PassConvVeto[i]
    	   && ele_MissHits[i] <= 1
           ) {
    	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.00609
    	  && fabs(ele_dPhi[i]) < 0.045
    	  && eleFull5x5SigmaIetaIeta[i] < 0.0298
    	  && ele_HoverE[i] < 0.0878
    	  && fabs(ele_d0[i]) < 0.1
    	  && fabs(ele_dZ[i]) < 0.2
    	  && fabs(ele_OneOverEminusOneOverP[i]) < 0.13
    	  && ele_PassConvVeto[i]
    	  && ele_MissHits[i] <= 1
    	  ) {
    	pass = true;
      }
    }
  }
  else if (EraName == "2017_94X") {

    // Tight ID recommended for analyses performed on 2017 data using 94X releases.
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
    // Based on these slides: https://indico.cern.ch/event/482677/contributions/2259342/attachments/1316731/1972911/talk_electron_ID_spring16_update.pdf

    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.0032
	   && fabs(ele_dPhi[i]) < 0.0547
	   && eleFull5x5SigmaIetaIeta[i] < 0.0106
	   && ele_HoverE[i] < 0.046 + 1.16 / eleE[i] + 0.0324*fixedGridRhoFastjetAll / eleE[i]
	   && fabs(ele_d0[i]) < 0.05
	   && fabs(ele_dZ[i]) < 0.10
	   && fabs(ele_OneOverEminusOneOverP[i]) < 0.184
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 1
	   ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.00632
	  && fabs(ele_dPhi[i]) < 0.0394
	  && eleFull5x5SigmaIetaIeta[i] < 0.0387
	  && ele_HoverE[i] < 0.0275 + 2.52 / eleE[i] + 0.183*fixedGridRhoFastjetAll / eleE[i]
	  && fabs(ele_d0[i]) < 0.1
	  && fabs(ele_dZ[i]) < 0.2
	  && fabs(ele_OneOverEminusOneOverP[i]) < 0.0721
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 1
	  ) {
	pass = true;
      }
    }
  }
  return pass;
}

bool RazorAnalyzer::passEGammaPOGTightElectronID(int i, bool use25nsCuts, string EraName){
  if (!use25nsCuts) {
    std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
    return false;
  }
  bool pass = false;

  if (EraName == "Spring15") {
    //For Moriond 2017, SUSY group decided to keep using Spring15 cut-based ID
    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.00926
	   && fabs(ele_dPhi[i]) < 0.0336
	   && eleFull5x5SigmaIetaIeta[i] < 0.0101
	   && ele_HoverE[i] < 0.0597
	   && fabs(ele_d0[i]) < 0.0111
	   && fabs(ele_dZ[i]) < 0.0466
	   && fabs(ele_OneOverEminusOneOverP[i]) < 0.012
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.00724
	  && fabs(ele_dPhi[i]) < 0.0918
	  && eleFull5x5SigmaIetaIeta[i] < 0.0279
	  && ele_HoverE[i] < 0.0615
	  && fabs(ele_d0[i]) < 0.0351
	  && fabs(ele_dZ[i]) < 0.417
	  && fabs(ele_OneOverEminusOneOverP[i]) < 0.00999
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 1
	  ) {
	pass = true;
      }
    }
  }
  else if (EraName == "Summer16") {
    // Tight ID recommended for analyses performed on 2016 data using 8XX releases.
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
    // Based on these slides: https://indico.cern.ch/event/482677/contributions/2259342/attachments/1316731/1972911/talk_electron_ID_spring16_update.pdf

    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.00308
	   && fabs(ele_dPhi[i]) < 0.0816
	   && eleFull5x5SigmaIetaIeta[i] < 0.00998
	   && ele_HoverE[i] < 0.0414
	   && fabs(ele_d0[i]) < 0.05
	   && fabs(ele_dZ[i]) < 0.10
	   && fabs(ele_OneOverEminusOneOverP[i]) < 0.0129
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 1
           ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.00605
	  && fabs(ele_dPhi[i]) < 0.0394
	  && eleFull5x5SigmaIetaIeta[i] < 0.0292
	  && ele_HoverE[i] < 0.0641
	  && fabs(ele_d0[i]) < 0.1
	  && fabs(ele_dZ[i]) < 0.2
	  && fabs(ele_OneOverEminusOneOverP[i]) < 0.0129
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 1
	  ) {
	pass = true;
      }
    }
  }
  else if (EraName == "2017_94X") {

    // Tight ID recommended for analyses performed on 2017 data using 94X releases.
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
    // Based on these slides: https://indico.cern.ch/event/732971/contributions/3022843/attachments/1658685/2656462/eleIdTuning.pdf

    if(fabs(eleEta_SC[i]) < 1.479) {
      if ( fabs(ele_dEta[i]) < 0.00255
	   && fabs(ele_dPhi[i]) < 0.022
	   && eleFull5x5SigmaIetaIeta[i] < 0.0104
	   && ele_HoverE[i] < 0.026 + 1.15 / eleE[i] + 0.0324*fixedGridRhoFastjetAll / eleE[i]
	   && fabs(ele_d0[i]) < 0.05
	   && fabs(ele_dZ[i]) < 0.10
	   && fabs(ele_OneOverEminusOneOverP[i]) < 0.159
	   && ele_PassConvVeto[i]
	   && ele_MissHits[i] <= 1
	   ) {
	pass = true;
      }
    } else {
      if (fabs(ele_dEta[i]) < 0.00501
	  && fabs(ele_dPhi[i]) < 0.0236
	  && eleFull5x5SigmaIetaIeta[i] < 0.0353
	  && ele_HoverE[i] < 0.0188 + 2.06 / eleE[i] + 0.183*fixedGridRhoFastjetAll / eleE[i]
	  && fabs(ele_d0[i]) < 0.1
	  && fabs(ele_dZ[i]) < 0.2
	  && fabs(ele_OneOverEminusOneOverP[i]) < 0.0197
	  && ele_PassConvVeto[i]
	  && ele_MissHits[i] <= 1
	  ) {
	pass = true;
      }
    }
  }

  return pass;
}

bool RazorAnalyzer::passMVANonTrigVetoElectronID(int i, string EraName){

  bool pass = false;
  if ( passMVAVetoElectronID(i, EraName)
       && fabs(ele_ip3dSignificance[i]) < 4
      ) {
    pass = true;
  }
  //cout << pass << "\n";

  return pass;

}

bool RazorAnalyzer::passMVAVetoElectronID(int i, string EraName){

  //Giovanni Zevi Della Porta
  //tried to match working points of 80X electron MVA with 74X electron MVA
  //pT < 10 bin should use "HZZ" MVA
  //pT > 10 bin should use "GeneralPurpose" MVA
  //https://indico.cern.ch/event/590228/contributions/2380031/attachments/1375541/2088587/EGMSUS_newIDs_17Nov16.pdf
  Int_t subdet = 0;
  if (fabs(eleEta_SC[i]) < 0.8) subdet = 0;
  else if (fabs(eleEta_SC[i]) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0 && elePt[i] <= 15.0) ptBin = 1;
  if (elePt[i] > 15.0 && elePt[i] <= 25.0) ptBin = 2;
  if (elePt[i] > 25.0 ) ptBin = 3;

  double MVACut = -999;

  if (EraName == "17Nov2017Rereco") {
    //Vinay Hegde
    //https://indico.cern.ch/event/719317/contributions/2963816/attachments/1630110/2598062/MVAidSUSY_10Apr18_SUSYMeeting.pdf
    //Working Point is called SUSY VLoose
    //Used for 2017 Data (17Nov2017Rereco)
    //pt 5-10
    if (subdet == 0 && ptBin == 0) MVACut = 0.488;
    if (subdet == 1 && ptBin == 0) MVACut = -0.045;
    if (subdet == 2 && ptBin == 0) MVACut = 0.176;
    //For pT 10-25
    if (subdet == 0 && (ptBin == 1 || ptBin == 2) ) MVACut = -0.788 - 0.00987 * ( elePt[i] - 10);
    if (subdet == 1 && (ptBin == 1 || ptBin == 2) ) MVACut = -0.85 - 0.005 * ( elePt[i] - 10);
    if (subdet == 2 && (ptBin == 1 || ptBin == 2) ) MVACut = -0.81 - 0.00513 * ( elePt[i] - 10);
    //For pT>25
    if (subdet == 0 && ptBin == 3) MVACut = -0.64;
    if (subdet == 1 && ptBin == 3) MVACut = -0.775;
    if (subdet == 2 && ptBin == 3) MVACut = -0.733;

  } else {
     //Giovanni Zevi Della Porta
    //tried to match working points of 80X electron MVA with 74X electron MVA
    //Working Point is called SUSY VLoose
    //Used for 2016 Data
    //pT < 10 bin should use "HZZ" MVA
    //pT > 10 bin should use "GeneralPurpose" MVA
    //https://indico.cern.ch/event/590228/contributions/2380031/attachments/1375541/2088587/EGMSUS_newIDs_17Nov16.pdf

    //pt 5-10
    if (subdet == 0 && ptBin == 0) MVACut = 0.46;
    if (subdet == 1 && ptBin == 0) MVACut = -0.03;
    if (subdet == 2 && ptBin == 0) MVACut = 0.06;
    //For pT 10-15
    if (subdet == 0 && ptBin == 1) MVACut = -0.48;
    if (subdet == 1 && ptBin == 1) MVACut = -0.67;
    if (subdet == 2 && ptBin == 1) MVACut = -0.49;
    //For pT 15-25
    if (subdet == 0 && ptBin == 2) MVACut = -0.48 - 0.037 * ( elePt[i] - 15);
    if (subdet == 1 && ptBin == 2) MVACut = -0.67 - 0.024 * ( elePt[i] - 15);
    if (subdet == 2 && ptBin == 2) MVACut = -0.49 - 0.034 * ( elePt[i] - 15);
    //For pT>25
    if (subdet == 0 && ptBin == 3) MVACut = -0.85;
    if (subdet == 1 && ptBin == 3) MVACut = -0.91;
    if (subdet == 2 && ptBin == 3) MVACut = -0.83;
  }

  double mvaVar = ele_IDMVAGeneralPurpose[i];
  if (ptBin == 0) mvaVar = ele_IDMVAHZZ[i];
  if (! (mvaVar == -99)){
    cout << ptBin << " " << subdet << " : " << ele_IDMVAGeneralPurpose[i] << " " << ele_IDMVAHZZ[i] << " --> " << mvaVar << " : cut = " <<  MVACut << " | pass = ";

  }

  bool pass = false;
  if (mvaVar > MVACut ) {
    pass = true;
  }
  //cout << pass << "\n";

  return pass;

}

bool RazorAnalyzer::passMVALooseElectronID(int i, string EraName){

  //Tuned to match signal efficiency of the EGM Cut-based Veto Working Point
  //pT < 10 bin should use "HZZ" MVA
  //pT > 10 bin should use "GeneralPurpose" MVA
  //
  Int_t subdet = 0;
  if (fabs(eleEta_SC[i]) < 0.8) subdet = 0;
  else if (fabs(eleEta_SC[i]) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0 && elePt[i] <= 15.0) ptBin = 1;
  if (elePt[i] > 15.0 && elePt[i] <= 25.0) ptBin = 2;
  if (elePt[i] > 25.0 ) ptBin = 3;

  double MVACut = -999;

  if (EraName == "17Nov2017Rereco") {
    //Working point tuned for WWZ analysis (it should be retuned for 2017 dataset)
    //For 2017 Data
    //pt 5-10
    if (subdet == 0 && ptBin == 0) MVACut = 0.83;
    if (subdet == 1 && ptBin == 0) MVACut = 0.76;
    if (subdet == 2 && ptBin == 0) MVACut = 0.80;
    //For pT 10-15
    if (subdet == 0 && ptBin == 1) MVACut = 0.66;
    if (subdet == 1 && ptBin == 1) MVACut = 0.44;
    if (subdet == 2 && ptBin == 1) MVACut = 0.49;
    //For pT 15-25
    if (subdet == 0 && ptBin == 2) MVACut = 0.66 - 0.064 * ( elePt[i] - 15);
    if (subdet == 1 && ptBin == 2) MVACut = 0.44 - 0.057 * ( elePt[i] - 15);
    if (subdet == 2 && ptBin == 2) MVACut = 0.49 - 0.030 * ( elePt[i] - 15);
    //For pT>25
    if (subdet == 0 && ptBin == 3) MVACut = 0.02;
    if (subdet == 1 && ptBin == 3) MVACut = -0.13;
    if (subdet == 2 && ptBin == 3) MVACut = 0.19;
  } else {
    //Working point tuned for WWZ analysis
    //For 2016 Data
    //pt 5-10
    if (subdet == 0 && ptBin == 0) MVACut = 0.83;
    if (subdet == 1 && ptBin == 0) MVACut = 0.76;
    if (subdet == 2 && ptBin == 0) MVACut = 0.80;
    //For pT 10-15
    if (subdet == 0 && ptBin == 1) MVACut = 0.66;
    if (subdet == 1 && ptBin == 1) MVACut = 0.44;
    if (subdet == 2 && ptBin == 1) MVACut = 0.49;
    //For pT 15-25
    if (subdet == 0 && ptBin == 2) MVACut = 0.66 - 0.064 * ( elePt[i] - 15);
    if (subdet == 1 && ptBin == 2) MVACut = 0.44 - 0.057 * ( elePt[i] - 15);
    if (subdet == 2 && ptBin == 2) MVACut = 0.49 - 0.030 * ( elePt[i] - 15);
    //For pT>25
    if (subdet == 0 && ptBin == 3) MVACut = 0.02;
    if (subdet == 1 && ptBin == 3) MVACut = -0.13;
    if (subdet == 2 && ptBin == 3) MVACut = 0.19;
  }

  double mvaVar = ele_IDMVAGeneralPurpose[i];
  if (ptBin == 0) mvaVar = ele_IDMVAHZZ[i];

  //cout << ptBin << " " << subdet << " : " << ele_IDMVAGeneralPurpose[i] << " " << ele_IDMVAHZZ[i] << " --> " << mvaVar << " : cut = " <<  MVACut << " | pass = ";

  bool pass = false;
  if (mvaVar > MVACut) {
    pass = true;
  }
  //cout << pass << "\n";

  return pass;

}

bool RazorAnalyzer::passEGammaPOGVetoElectronIso(int i, bool use25nsCuts){
    // Recommended for analyses performed on 2016 data using 8XX releases.
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;

    if(fabs(eleEta_SC[i]) < 1.479) {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.175
           ) {
            pass = true;
        }
    } else {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.159
           ) {
            pass = true;
        }
    }
    return pass;
}

bool RazorAnalyzer::passEGammaPOGLooseElectronIso(int i, bool use25nsCuts){
    // Recommended for analyses performed on 2016 data using 8XX releases.
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;

    if(fabs(eleEta_SC[i]) < 1.479) {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0994
           ) {
            pass = true;
        }
    } else {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.107
           ) {
            pass = true;
        }
    }
    return pass;
}

bool RazorAnalyzer::passEGammaPOGMediumElectronIso(int i, bool use25nsCuts){
    // Recommended for analyses performed on 2016 data using 8XX releases.
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;

    if(fabs(eleEta_SC[i]) < 1.479) {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0695
           ) {
            pass = true;
        }
    } else {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0821
           ) {
            pass = true;
        }
    }
    return pass;
}

bool RazorAnalyzer::passEGammaPOGTightElectronIso(int i, bool use25nsCuts){
    // Recommended for analyses performed on 2016 data using 8XX releases.
    if (!use25nsCuts) {
        std::cerr << "Error: 50ns cuts are not implemented for this electron ID" << std::endl;
        return false;
    }
    bool pass = false;

    if(fabs(eleEta_SC[i]) < 1.479) {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0588
           ) {
            pass = true;
        }
    } else {
        if ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) / elePt[i] < 0.0571
           ) {
            pass = true;
        }
    }
    return pass;
}

bool RazorAnalyzer::passMVANonTrigVetoElectronIso(int i){

  bool pass = false;
  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
  if (  ( (elePt[i] > 20 && (ele_chargedMiniIso[i] + fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveAreaMean(i)*pow(dr/0.3,2)))/elePt[i] < 0.2 )
	  ||
	   (elePt[i] <= 20 && (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - GetElectronEffectiveAreaMean(i)*fixedGridRhoFastjetAll)) < 5)
	  )
	) {
    pass = true;
  }
  return pass;
}

bool RazorAnalyzer::passHZZElectronPreselection(int i){

  Int_t subdet = 0;
  if (fabs(eleEta_SC[i]) < 0.8) subdet = 0;
  else if (fabs(eleEta_SC[i]) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0) ptBin = 1;

  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = 0.47;
  if (subdet == 1 && ptBin == 0) MVACut = 0.004;
  if (subdet == 2 && ptBin == 0) MVACut = 0.295;
  if (subdet == 0 && ptBin == 1) MVACut = 0.5;
  if (subdet == 1 && ptBin == 1) MVACut = 0.12;
  if (subdet == 2 && ptBin == 1) MVACut = 0.6;

  bool pass = false;
  if (ele_IDMVAHZZ[i] > MVACut
      && fabs(ele_d0[i]) < 0.5
      && fabs(ele_dZ[i]) < 1.0
      ) {
    pass = true;
  }

  return pass;
}


bool RazorAnalyzer::isHZZElectron(int i){

  Int_t subdet = 0;
  if (fabs(eleEta_SC[i]) < 0.8) subdet = 0;
  else if (fabs(eleEta_SC[i]) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (elePt[i] > 10.0) ptBin = 1;

  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = 0.47;
  if (subdet == 1 && ptBin == 0) MVACut = 0.004;
  if (subdet == 2 && ptBin == 0) MVACut = 0.295;
  if (subdet == 0 && ptBin == 1) MVACut = 0.5;
  if (subdet == 1 && ptBin == 1) MVACut = 0.12;
  if (subdet == 2 && ptBin == 1) MVACut = 0.6;

  bool pass = false;
  if (ele_IDMVAHZZ[i] > MVACut
      && fabs(ele_d0[i]) < 0.5
      && fabs(ele_dZ[i]) < 1.0
      && fabs(ele_d0[i]) < 0.05
      && ( (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i] < 0.4)
      ) {
    pass = true;
  }

  return pass;
}




bool RazorAnalyzer::matchTagElectronHLTFilters(int i){
  bool match = false;
  if (
      //Data filters
      ele_passHLTFilter[i][1]
      || ele_passHLTFilter[i][5] || ele_passHLTFilter[i][6] || ele_passHLTFilter[i][12] || ele_passHLTFilter[i][13]
      || ele_passHLTFilter[i][49] || ele_passHLTFilter[i][53] || ele_passHLTFilter[i][57] || ele_passHLTFilter[i][60]
      //MC filters
      || ele_passHLTFilter[i][64]
      || ele_passHLTFilter[i][3] || ele_passHLTFilter[i][8] || ele_passHLTFilter[i][10] || ele_passHLTFilter[i][15]
       ) {
    match = true;
  }

  return match;
}

bool RazorAnalyzer::matchProbeElectronHLTFilters(int i){
  bool match = false;
  if (
      ele_passHLTFilter[i][50] || ele_passHLTFilter[i][51]
      || ele_passHLTFilter[i][61]     || ele_passHLTFilter[i][62]
       ) {
    match = true;
  }

  return match;
}

bool RazorAnalyzer::matchProbeSCHLTFilters(int i){
  bool match = false;
  if (
      ele_passHLTFilter[i][54] || ele_passHLTFilter[i][55]
      || ele_passHLTFilter[i][58]    || ele_passHLTFilter[i][59]
       ) {
    match = true;
  }

  return match;
}

bool RazorAnalyzer::matchElectronHLTFilters(int i, string HLTFilter, string analysisTag) {
  if (analysisTag == "2015") return matchElectronHLTFilters2015(i, HLTFilter);
  else if (analysisTag == "2016") return matchElectronHLTFilters2016(i, HLTFilter);
  else {
    cout << "Analysis Tag " << analysisTag << " is not supported. Returning false.\n";
    return false;
  }
}


bool RazorAnalyzer::matchElectronHLTFilters2016(int i, string HLTFilter){
  bool match = false;

  if (HLTFilter == "SingleElectron") {
    if (
	//Data filters
	ele_passHLTFilter[i][1] || ele_passHLTFilter[i][5] || ele_passHLTFilter[i][6] || ele_passHLTFilter[i][12] || ele_passHLTFilter[i][13]
	//MC filters
	|| ele_passHLTFilter[i][3] || ele_passHLTFilter[i][8] || ele_passHLTFilter[i][10] || ele_passHLTFilter[i][15]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "Ele23Loose") {
    if (
	//Data filters
	ele_passHLTFilter[i][1]
	 ) {
      match = true;
    }
  }
  if (HLTFilter == "Ele27Loose") {
    if (
	//Data filters
	ele_passHLTFilter[i][5] || ele_passHLTFilter[i][7]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "Ele27Tight") {
    if (
	//Data filters
	ele_passHLTFilter[i][6]  || ele_passHLTFilter[i][8]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "Ele32Loose") {
    if (
	//Data filters
	ele_passHLTFilter[i][12]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "Ele32Tight") {
    if (
	//Data filters
	ele_passHLTFilter[i][13]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "Ele105") {
    if ( ele_passHLTFilter[i][17] ) match = true;
  }
  if (HLTFilter == "Ele115") {
    if ( ele_passHLTFilter[i][19] ) match = true;
  }

  if (HLTFilter == "DoubleElectronLeg1") {
    if (
	//Data filters
	ele_passHLTFilter[i][24] || ele_passHLTFilter[i][28]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "DoubleElectronLeg2") {
    if (
	//Data filters
	ele_passHLTFilter[i][25] || ele_passHLTFilter[i][29]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "DoubleElectronLeg2DZ") {
    if (
	//Data filters
	ele_passHLTFilter[i][26] || ele_passHLTFilter[i][30]
	 ) {
      match = true;
    }
  }

  return match;
}

bool RazorAnalyzer::matchElectronHLTFilters2015(int i, string HLTFilter){
  bool match = false;

  if (HLTFilter == "SingleElectron") {
    if (
	//Data filters
	ele_passHLTFilter[i][1] || ele_passHLTFilter[i][5] || ele_passHLTFilter[i][6] || ele_passHLTFilter[i][12] || ele_passHLTFilter[i][13]
	//MC filters
	|| ele_passHLTFilter[i][3] || ele_passHLTFilter[i][8] || ele_passHLTFilter[i][10] || ele_passHLTFilter[i][15]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "Ele23Loose") {
    if (
	//Data filters
	ele_passHLTFilter[i][1]
	//MC filters
	|| ele_passHLTFilter[i][3] || ele_passHLTFilter[i][64]
	 ) {
      match = true;
    }
  }
  if (HLTFilter == "Ele27Loose") {
    if (
	//Data filters
	ele_passHLTFilter[i][5]
	//MC filters
	|| ele_passHLTFilter[i][8]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "Ele27Tight") {
    if (
	//Data filters
	ele_passHLTFilter[i][6]
	//MC filters
	|| ele_passHLTFilter[i][10]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "Ele32Tight") {
    if (
	//Data filters
	ele_passHLTFilter[i][13]
	//MC filters
	|| ele_passHLTFilter[i][15]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "Ele105") {
    if ( ele_passHLTFilter[i][17] ) match = true;
  }
  if (HLTFilter == "Ele115") {
    if ( ele_passHLTFilter[i][19] ) match = true;
  }

  if (HLTFilter == "DoubleElectronLeg1") {
    if (
	//Data filters
	ele_passHLTFilter[i][24] || ele_passHLTFilter[i][28]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "DoubleElectronLeg2") {
    if (
	//Data filters
	ele_passHLTFilter[i][25] || ele_passHLTFilter[i][29]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "DoubleElectronLeg2DZ") {
    if (
	//Data filters
	ele_passHLTFilter[i][26] || ele_passHLTFilter[i][30]
	 ) {
      match = true;
    }
  }

  return match;
}

//////////////////////////////
// MUON
//////////////////////////////

float RazorAnalyzer::GetMuonEffectiveAreaMean(int i, string type ){

  double effArea = 0.0;
  //Effective areas below are for the sum of Neutral Hadrons + Photons
  if (type == "neutral") {
    if (fabs(muonEta[i]) < 0.8) {
      effArea = 0.0735;
    } else if (fabs(muonEta[i]) < 1.3) {
      effArea = 0.0619;
    } else if (fabs(muonEta[i]) < 2.0) {
      effArea = 0.0465;
    } else if (fabs(muonEta[i]) < 2.2) {
      effArea = 0.0433;
    } else {
      effArea = 0.0577;
    }
  }
  if (type == "charged") {
    if (fabs(muonEta[i]) < 0.8) {
      effArea = 0.0106;
    } else if (fabs(muonEta[i]) < 1.3) {
      effArea = 0.0096;
    } else if (fabs(muonEta[i]) < 2.0) {
      effArea = 0.0079;
    } else if (fabs(muonEta[i]) < 2.2) {
      effArea = 0.0058;
    } else {
      effArea = 0.0053;
    }
  }
  return effArea;
}


float RazorAnalyzer::GetMuonEffectiveArea90(int i, string EraName ){

  double effArea = 0.0;
  //Effective areas below are for the sum of Neutral Hadrons + Photons
  if (EraName == "Spring15" || EraName == "Spring16" || EraName == "Summer16") {
  //values from https://github.com/cms-data/PhysicsTools-NanoAOD/blob/10e7935ba38c2172ebb75979dcc1d8174b0566cd/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt
    if (fabs(muonEta[i]) < 0.8) {
      effArea = 0.0735;
    } else if (fabs(muonEta[i]) < 1.3) {
      effArea = 0.0619;
    } else if (fabs(muonEta[i]) < 2.0) {
      effArea = 0.0465;
    } else if (fabs(muonEta[i]) < 2.2) {
      effArea = 0.0433;
    } else {
      effArea = 0.0577;
    }
  }
  else if ( EraName == "2017_94X") {
  //values from https://github.com/cms-data/PhysicsTools-NanoAOD/blob/10e7935ba38c2172ebb75979dcc1d8174b0566cd/effAreaMuons_cone03_pfNeuHadronsAndPhotons_94X.txt
    if (fabs(muonEta[i]) < 0.8) {
      effArea = 0.0566;
    } else if (fabs(muonEta[i]) < 1.3) {
      effArea = 0.0562;
    } else if (fabs(muonEta[i]) < 2.0) {
      effArea = 0.0363;
    } else if (fabs(muonEta[i]) < 2.2) {
      effArea = 0.019;
    } else {
      effArea = 0.0064;
    }
  }
  return effArea;
}



bool RazorAnalyzer::isMuonPOGLooseMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!(muonIsLoose[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!((muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.2)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isMuonPOGMediumMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!(muonIsICHEPMedium[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!((muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.12)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isMuonPOGTightMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  if (applyID) {
    if (!(muonIsTight[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!((muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.12)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isVetoMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/muonPt[i]));
  if (applyID) {
    if (!(muonIsLoose[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!(
	  ( muonPt[i] > 20 && (muon_chargedMiniIso[i] + fmax(0.0, muon_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetMuonEffectiveAreaMean(i,"neutral")*pow(dr/0.3,2)) )/muonPt[i] < 0.2 )
	  ||
	  ( muonPt[i] <= 20 && (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - fixedGridRhoFastjetAll*GetMuonEffectiveAreaMean(i,"neutral") )) < 10 )
	  )) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isLooseMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/muonPt[i]));
  if (applyID) {
    if (!(muonIsLoose[i] && fabs(muon_ip3dSignificance[i]) < 4)) pass = false;
  }
  if (applyIso) {
    if (!( (muon_chargedMiniIso[i] + fmax(0.0, muon_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetMuonEffectiveAreaMean(i,"neutral")*pow(dr/0.3,2)))/muonPt[i] < 0.2)) pass = false;
  }
  return pass;
}

bool RazorAnalyzer::isTightMuon(int i, bool applyID, bool applyIso){
  bool pass = true;
  double dr = fmax(0.05,fmin(0.2, 10/muonPt[i]));
  if (applyID) {
    if (!(muonIsMedium[i] && fabs(muon_ip3dSignificance[i]) < 4 && fabs(muon_d0[i]) < 0.2)) pass = false;
  }
  if (applyIso) {
    if (!( (muon_chargedMiniIso[i] + fmax(0.0, muon_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetMuonEffectiveAreaMean(i,"neutral")*pow(dr/0.3,2)) )/muonPt[i] < 0.2)) pass = false;
  }
  return pass;
}




bool RazorAnalyzer::passHZZMuonPreselection(int i){
  bool pass = false;
  if (muonIsLoose[i]
      && fabs(muon_d0[i]) < 0.5
      && fabs(muon_dZ[i]) < 1.0
      ) {
    pass = true;
  }
  return pass;
}


bool RazorAnalyzer::isHZZMuon(int i){
  bool pass = false;
  if (muonIsLoose[i]
      && fabs(muon_d0[i]) < 0.5
      && fabs(muon_dZ[i]) < 1.0
      && fabs(muon_ip3dSignificance[i]) < 4
      && ( ( (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.4 ))
      ) {
    pass = true;
  }
  return pass;
}



bool RazorAnalyzer::matchTagMuonHLTFilters(int i){
  return  matchMuonHLTFilters(i, "SingleMuon");
}


bool RazorAnalyzer::matchMuonHLTFilters(int i, string HLTFilter){
  bool match = false;

  if (HLTFilter == "SingleMuon") {
    if (
	//Data filters
	muon_passHLTFilter[i][1]      // IsoMu17_eta2p1
	|| muon_passHLTFilter[i][3]   // IsoMu20
	|| muon_passHLTFilter[i][4]   // HLT_IsoMu20_eta2p1
	|| muon_passHLTFilter[i][6]   // HLT_IsoMu24_eta2p1
	|| muon_passHLTFilter[i][8]   // IsoMu27
	|| muon_passHLTFilter[i][9]   // IsoTkMu20
	|| muon_passHLTFilter[i][10]  // IsoTkMu20_eta2p1
	|| muon_passHLTFilter[i][11]  // IsoTkMu24_eta2p1
	|| muon_passHLTFilter[i][12]  // IsoTkMu27
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "IsoMu20") {
    if ( muon_passHLTFilter[i][3] ) match = true;
  }
  if (HLTFilter == "IsoMu22") {
    if ( muon_passHLTFilter[i][7] || muon_passHLTFilter[i][9] ) match = true;
  }
  if (HLTFilter == "IsoMu24") {
    if ( muon_passHLTFilter[i][11] || muon_passHLTFilter[i][12] ) match = true;
  }
  if (HLTFilter == "IsoTkMu20") {
    if ( muon_passHLTFilter[i][16] ) match = true;
  }
  if (HLTFilter == "IsoTkMu22") {
    if ( muon_passHLTFilter[i][17] || muon_passHLTFilter[i][18] ) match = true;
  }
  if (HLTFilter == "IsoTkMu24") {
    if ( muon_passHLTFilter[i][19] ) match = true;
  }
  if (HLTFilter == "IsoMu27") {
    if ( muon_passHLTFilter[i][14] ) match = true;
  }
  if (HLTFilter == "IsoTkMu27") {
    if ( muon_passHLTFilter[i][21] ) match = true;
  }
  if (HLTFilter == "Mu50") {
    if ( muon_passHLTFilter[i][22] ) match = true;
  }
  if (HLTFilter == "Mu55") {
    if ( muon_passHLTFilter[i][25] ) match = true;
  }
  if (HLTFilter == "Mu45_eta2p1") {
    if ( muon_passHLTFilter[i][26] ) match = true;
  }
  if (HLTFilter == "TkMu50") {
    if ( muon_passHLTFilter[i][23] ) match = true;
  }


  if (HLTFilter == "DoubleMuonLeg1") {
    if (
	//Data filters : 17_8 Triggers
	muon_passHLTFilter[i][19] || muon_passHLTFilter[i][27]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "DoubleMuonLeg2") {
    if (
	//Data filters : 17_8 triggers
	muon_passHLTFilter[i][25] || muon_passHLTFilter[i][28]
	 ) {
      match = true;
    }
  }

  if (HLTFilter == "DoubleMuonLeg2DZ") {
    if (
	//Data filters : 17_8 triggers
	muon_passHLTFilter[i][26] || muon_passHLTFilter[i][29]
	 ) {
      match = true;
    }
  }

  return match;
}

//////////////////////////////
// TAU
//////////////////////////////

bool RazorAnalyzer::isLooseTau(int i){
  bool pass = false;
  if (tau_IsLoose[i]) {
    pass = true;
  }

  return pass;
}

bool RazorAnalyzer::isMediumTau(int i){
  bool pass = false;
  if (tau_IsMedium[i]) {
    pass = true;
  }

  return pass;
}

bool RazorAnalyzer::isTightTau(int i){
  bool pass = false;
  if (tau_IsTight[i]) {
    pass = true;
  }

  return pass;
}

//////////////////////////////
// JET
//////////////////////////////

//From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80X
bool RazorAnalyzer::isCSVL(int i, string dataset)
{
  if (dataset == "74X")
  {
    return (jetCISV[i] > 0.605);
  }
  else if (dataset == "80X" )
  {
    return (jetCISV[i] > 0.5426);
  }
  else if (dataset == "94X" )
  {
    return (jetCISV[i] > 0.5803);
  }
  else
  {
    std::cerr << "reached end of comparison and found not compatible tag in isCSVL" << std::endl;
    return false;
  }
};

bool RazorAnalyzer::isCSVM(int i, string dataset)
{
  if (dataset == "74X")
  {
    return jetCISV[i] > 0.890;
  }
  else if (dataset == "80X" )
  {
    return jetCISV[i] > 0.8484;
  }
  else if (dataset == "94X" )
  {
    return jetCISV[i] > 0.8838;
  }
  else
  {
    std::cerr << "reached end of comparison and found not compatible tag in isCSVL" << std::endl;
    return false;
  }
};

bool RazorAnalyzer::isCSVT(int i, string dataset)
{
  if (dataset == "74X")
  {
    return jetCISV[i] > 0.970;
  }
  else if (dataset == "80X" )
  {
    return jetCISV[i] > 0.9535;
  }
  else if (dataset == "94X" )
  {
    return jetCISV[i] > 0.9693;
  }
  else
  {
    std::cerr << "reached end of comparison and found not compatible tag in isCSVL" << std::endl;
    return false;
  }

};

//Jet Energy Corrections
double RazorAnalyzer::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
						 double rho, double jetArea,
						 int run,
						 std::vector<std::pair<int,int> > JetCorrectionsIOV,
						 std::vector<FactorizedJetCorrector*> jetcorrector,
						 int jetCorrectionLevel,
						 bool printDebug) {

  int foundIndex = -1;
  for (unsigned int i=0; i<JetCorrectionsIOV.size(); i++) {
    if (run >= JetCorrectionsIOV[i].first && run <= JetCorrectionsIOV[i].second) {
      foundIndex = i;
    }
  }
  if (foundIndex == -1) {
    cout << "Warning: run = " << run << " was not found in any valid IOV range. use default index = 0 for Jet energy corrections. \n";
    foundIndex = 0;
  }

  if (!jetcorrector[foundIndex]) {
    cout << "WWARNING: Jet corrector pointer is null. Returning JEC = 0. \n";
    return 0;
  }

  jetcorrector[foundIndex]->setJetEta(jetEta);
  jetcorrector[foundIndex]->setJetPt(jetRawPt);
  jetcorrector[foundIndex]->setJetPhi(jetPhi);
  jetcorrector[foundIndex]->setJetE(jetE);
  jetcorrector[foundIndex]->setRho(rho);
  jetcorrector[foundIndex]->setJetA(jetArea);

  std::vector<float> corrections;
  corrections = jetcorrector[foundIndex]->getSubCorrections();

  if (printDebug) cout << "Computing Jet Energy Corrections for jet with raw momentum: " << jetRawPt << " " << jetEta << " " << jetPhi << "\n";

  double cumulativeCorrection = 1.0;
  for (UInt_t j=0; j<corrections.size(); ++j) {

    //only correct up to the required level. if -1, then do all correction levels
    if (jetCorrectionLevel >= 0 && int(j) > jetCorrectionLevel) continue;

    double currentCorrection = corrections.at(j)/cumulativeCorrection;
    cumulativeCorrection = corrections.at(j);
    if (printDebug) cout << "Correction Level " << j << " : current correction = " << currentCorrection << " , cumulative correction = " << cumulativeCorrection << "\n";
  }
  if (printDebug) cout << "Final Cumulative Correction: " << cumulativeCorrection << "\n";

  return cumulativeCorrection;

}

//Jet Energy Corrections
double RazorAnalyzer::JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
						 double rho, double jetArea,
						 FactorizedJetCorrector *jetcorrector,
						 int jetCorrectionLevel,
						 bool printDebug) {
  if (!jetcorrector) {
    cout << "WWARNING: Jet corrector pointer is null. Returning JEC = 0. \n";
    return 0;
  }

  jetcorrector->setJetEta(jetEta);
  jetcorrector->setJetPt(jetRawPt);
  jetcorrector->setJetPhi(jetPhi);
  jetcorrector->setJetE(jetE);
  jetcorrector->setRho(rho);
  jetcorrector->setJetA(jetArea);

  std::vector<float> corrections;
  corrections = jetcorrector->getSubCorrections();

  if (printDebug) cout << "Computing Jet Energy Corrections for jet with raw momentum: " << jetRawPt << " " << jetEta << " " << jetPhi << "\n";

  double cumulativeCorrection = 1.0;
  for (UInt_t j=0; j<corrections.size(); ++j) {

    //only correct up to the required level. if -1, then do all correction levels
    if (jetCorrectionLevel >= 0 && int(j) > jetCorrectionLevel) continue;

    double currentCorrection = corrections.at(j)/cumulativeCorrection;
    cumulativeCorrection = corrections.at(j);
    if (printDebug) cout << "Correction Level " << j << " : current correction = " << currentCorrection << " , cumulative correction = " << cumulativeCorrection << "\n";
  }
  if (printDebug) cout << "Final Cumulative Correction: " << cumulativeCorrection << "\n";

  return cumulativeCorrection;

}


//compute the smeared jet pt (if option = "up" or "down", will change the smear factor by +/- 1 sigma )
//NOTE: these are Run 1 recommendations and should be replaced as soon as a Run 2 prescription is available.
double RazorAnalyzer::JetEnergySmearingFactor( double jetPt, double jetEta, double NPU, SimpleJetResolution *JetResolutionCalculator, TRandom3 *random ) {

  std::vector<float> fJetEta, fJetPtNPU;
  fJetEta.push_back(jetEta);
  fJetPtNPU.push_back(jetPt);
  fJetPtNPU.push_back(NPU);
  double MCJetResolution = JetResolutionCalculator->resolution(fJetEta,fJetPtNPU);

  double c = 1;
  if (fabs(jetEta) < 0.5) c = 1.079;
  else if(fabs(jetEta) < 1.1) c = 1.099;
  else if(fabs(jetEta) < 1.7) c = 1.121;
  else if(fabs(jetEta) < 2.3) c = 1.208;
  else if(fabs(jetEta) < 2.8) c = 1.254;
  else if(fabs(jetEta) < 3.2) c = 1.395;
  else if(fabs(jetEta) < 5.0) c = 1.056;

  double sigma = sqrt( c*c - 1) * MCJetResolution;

  return fmax( 1.0 + random->Gaus(0, sigma) , 0);

}

//return smearing factor for up/down shift in JER
double RazorAnalyzer::UpDownJetEnergySmearingFactor(double unsmearedPt, double jetEta, double NPU, SimpleJetResolution *JetResolutionCalculator, double smearedPt, string option){
    //get jet resolution
    std::vector<float> fJetEta, fJetPtNPU;
    fJetEta.push_back(jetEta);
    fJetPtNPU.push_back(unsmearedPt);
    fJetPtNPU.push_back(NPU);
    double MCJetResolution = JetResolutionCalculator->resolution(fJetEta,fJetPtNPU);

    //get sigma used to smear the jet
    double c = 1;
    if (fabs(jetEta) < 0.5) c = 1.079;
    else if(fabs(jetEta) < 1.1) c = 1.099;
    else if(fabs(jetEta) < 1.7) c = 1.121;
    else if(fabs(jetEta) < 2.3) c = 1.208;
    else if(fabs(jetEta) < 2.8) c = 1.254;
    else if(fabs(jetEta) < 3.2) c = 1.395;
    else if(fabs(jetEta) < 5.0) c = 1.056;
    double sigma = sqrt( c*c - 1) * MCJetResolution;
    //get number of sigmas the jet was smeared
    double z = (smearedPt / unsmearedPt - 1) /sigma;

    if(option == "up"){ //get c plus 1 sigma
        double cUp = 1.0;
        if (fabs(jetEta) < 0.5) cUp = 1.105;
        else if(fabs(jetEta) < 1.1) cUp = 1.127;
        else if(fabs(jetEta) < 1.7) cUp = 1.150;
        else if(fabs(jetEta) < 2.3) cUp = 1.254;
        else if(fabs(jetEta) < 2.8) cUp = 1.316;
        else if(fabs(jetEta) < 3.2) cUp = 1.458;
        else if(fabs(jetEta) < 5.0) cUp = 1.247;
        double sigmaUp = sqrt( cUp*cUp - 1) * MCJetResolution;
        return 1.0 + z*sigmaUp;
    }
    else if(option == "down"){ //get c minus 1 sigma
        double cDown = 1.0;
        if (fabs(jetEta) < 0.5) cDown = 1.053;
        else if(fabs(jetEta) < 1.1) cDown = 1.071;
        else if(fabs(jetEta) < 1.7) cDown = 1.092;
        else if(fabs(jetEta) < 2.3) cDown = 1.162;
        else if(fabs(jetEta) < 2.8) cDown = 1.192;
        else if(fabs(jetEta) < 3.2) cDown = 1.332;
        else if(fabs(jetEta) < 5.0) cDown = 0.865;
        double sigmaDown = sqrt( cDown*cDown - 1) * MCJetResolution;
        return 1.0 + z*sigmaDown;
    }
    else{
        std::cout << "Error in UpDownJetEnergySmear: please specify option='up' or 'down'.  Returning 1.0" << std::endl;
    }
    return 1.0;
}


//b-tagging scale factors from https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation53XReReco/SFb-pt_WITHttbar_payload_EPS13.txt
//(if option = "up" or "down", will change the scale factor by +/- 1 sigma)
//NOTE: These are the Run 1 recommended scale factors and should be replaced as soon as Run 2 factors are available.
double RazorAnalyzer::BTagScaleFactor(double jetPt, bool CSVM, string option){
    double tmpBTagCorrFactor = 1.0;
    //nominal correction factor
    double tmpCorrFactor = 0.938887 + 0.00017124 * jetPt + (-2.76366e-07) * jetPt * jetPt ;

    if(option == "up" || option == "down"){
        double uncertainty = 0.0;
        if (jetPt < 30) uncertainty = 0.0415707;
        else if (jetPt < 40) uncertainty = 0.0204209;
        else if (jetPt < 50) uncertainty = 0.0223227;
        else if (jetPt < 60) uncertainty = 0.0206655;
        else if (jetPt < 70) uncertainty = 0.0199325;
        else if (jetPt < 80) uncertainty = 0.0174121;
        else if (jetPt < 100) uncertainty = 0.0202332;
        else if (jetPt < 120) uncertainty = 0.0182446;
        else if (jetPt < 160) uncertainty = 0.0159777;
        else if (jetPt < 210) uncertainty = 0.0218531;
        else if (jetPt < 260) uncertainty = 0.0204688;
        else if (jetPt < 320) uncertainty = 0.0265191;
        else if (jetPt < 400) uncertainty = 0.0313175;
        else if (jetPt < 500) uncertainty = 0.0415417;
        else if (jetPt < 600) uncertainty = 0.0740446;
        else if (jetPt < 800) uncertainty = 0.0596716;
        else uncertainty = 2*0.0596716;

        if (option == "up") tmpCorrFactor += uncertainty;
        else if (option == "down") tmpCorrFactor -= uncertainty;
    }

    double MCEff = 1.0;
    if (jetPt < 50) MCEff = 0.65;
    else if (jetPt < 80) MCEff = 0.70;
    else if (jetPt < 120) MCEff = 0.73;
    else if (jetPt < 210) MCEff = 0.73;
    else MCEff = 0.66;

    //If pass CSV Medium
    if (CSVM) {
        tmpBTagCorrFactor = tmpCorrFactor;
    } else {
        tmpBTagCorrFactor = (1/MCEff - tmpCorrFactor) / (1/MCEff - 1);
    }
    return tmpBTagCorrFactor;
}

//////////////////////////////
// PHOTON
//////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
bool RazorAnalyzer::photonPassesElectronVeto(int i){
    //use presence of a pixel seed as proxy for an electron veto
    return (pho_passEleVeto[i]);
}
/*
double RazorAnalyzer::getPhotonSminorSmajor(int ind_pho, bool isSminor){
	double etaAverage = 0.0;
        double phiAverage = 0.0;
        double mTotalWeight = 0.0;
        double tmpSumE = 0.0;
        double phoSetaeta = 0.0;
	double phoSphiphi = 0.0;
        double phoSetaphi = 0.0;
        double phoSminor = 0.0;
        double phoSmajor = 0.0;


	if(ecalRechit_ID->empty()) return 0.0;


	uint seedhitIndex =  (*pho_SeedRechitIndex)[ind_pho];
	bool isFromEB = bool( (*ecalRechit_ID)[seedhitIndex] < 840000000 );


	for (uint k=0; k<(*pho_EcalRechitIndex)[ind_pho].size(); ++k)
        {
		uint rechitIndex = (*pho_EcalRechitIndex)[ind_pho][k];
		if((*ecalRechit_E)[rechitIndex] < 1.0) continue;
		tmpSumE += (*ecalRechit_E)[rechitIndex];
	}
	for (uint k=0; k<(*pho_EcalRechitIndex)[ind_pho].size(); ++k)
        {
                uint rechitIndex = (*pho_EcalRechitIndex)[ind_pho][k];
                if((*ecalRechit_E)[rechitIndex] < 1.0) continue;
                double thisWeight = max(4.2+log(((*ecalRechit_E)[rechitIndex])/tmpSumE),0.0);
                mTotalWeight += thisWeight;
	        float thisIPhiIY =  1.0 * iPhi_or_iY_from_detID( (*ecalRechit_ID)[rechitIndex] , isFromEB);
                float thisIEtaIX =  1.0 * iEta_or_iX_from_detID( (*ecalRechit_ID)[rechitIndex] , isFromEB);

                etaAverage += thisWeight * (thisIEtaIX) ;
                phiAverage += thisWeight * (thisIPhiIY) ;
        }

	etaAverage = etaAverage / mTotalWeight;
        phiAverage = phiAverage / mTotalWeight;

	for (uint k=0; k<(*pho_EcalRechitIndex)[ind_pho].size(); ++k)
        {
                uint rechitIndex = (*pho_EcalRechitIndex)[ind_pho][k];
                if((*ecalRechit_E)[rechitIndex] < 1.0) continue;
                double thisWeight = max(4.2+log(((*ecalRechit_E)[rechitIndex])/tmpSumE),0.0);
		float thisIPhiIY =  1.0 * iPhi_or_iY_from_detID( (*ecalRechit_ID)[rechitIndex] , isFromEB);
                float thisIEtaIX =  1.0 * iEta_or_iX_from_detID( (*ecalRechit_ID)[rechitIndex] , isFromEB);
		phoSetaeta += thisWeight * (thisIEtaIX - etaAverage) * (thisIEtaIX - etaAverage);
                phoSphiphi += thisWeight * (thisIPhiIY - phiAverage) * (thisIPhiIY - phiAverage);
                phoSetaphi += thisWeight * (thisIEtaIX - etaAverage) * (thisIPhiIY - phiAverage);
        }
	if(tmpSumE>0.0)
        {
		phoSetaeta = phoSetaeta / mTotalWeight;
                phoSphiphi = phoSphiphi / mTotalWeight;
                phoSetaphi = phoSetaphi / mTotalWeight;
                phoSminor = 0.5 * (phoSetaeta + phoSphiphi - pow(pow(phoSetaeta - phoSphiphi,2.0) + 4.0*pow(phoSetaphi,2.0),0.5));
                phoSmajor = 0.5 * (phoSetaeta + phoSphiphi + pow(pow(phoSetaeta - phoSphiphi,2.0) + 4.0*pow(phoSetaphi,2.0),0.5));
        }

	if (isSminor) return phoSminor;
	else return phoSmajor;
}
*/

//Spring 15 values
void RazorAnalyzer::getPhotonEffAreaRun2( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
  if( fabs( eta ) < 1.0 )
    {
      effAreaChHad = 0.0157;
      effAreaNHad  = 0.0143;
      effAreaPho   = 0.0725;
    }
  else if( fabs( eta ) < 1.479 )
    {
      effAreaChHad = 0.0143;
      effAreaNHad  = 0.0210;
      effAreaPho   = 0.0604;
    }
  else if( fabs( eta ) < 2.0 )
    {
      effAreaChHad = 0.0115;
      effAreaNHad  = 0.0148;
      effAreaPho   = 0.0320;
    }
  else if( fabs( eta ) < 2.2 )
    {
      effAreaChHad = 0.0094;
      effAreaNHad  = 0.0082;
      effAreaPho   = 0.0512;
    }
  else if( fabs( eta ) < 2.3 )
    {
      effAreaChHad = 0.0095;
      effAreaNHad  = 0.0124;
      effAreaPho   = 0.0766;
    }
  else if( fabs( eta ) < 2.4 )
    {
      effAreaChHad = 0.0068;
      effAreaNHad  = 0.0186;
      effAreaPho   = 0.0949;
    }
  else
    {
      effAreaChHad = 0.0053;
      effAreaNHad  = 0.0320;
      effAreaPho   = 0.1160;
    }
};

// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
void RazorAnalyzer::getPhotonEffArea90( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
  if( fabs( eta ) < 1.0 )
    {
      effAreaChHad = 0.0360;
      effAreaNHad  = 0.0597;
      effAreaPho   = 0.1210;
    }
  else if( fabs( eta ) < 1.479 )
    {
      effAreaChHad = 0.0377;
      effAreaNHad  = 0.0807;
      effAreaPho   = 0.1107;
    }
  else if( fabs( eta ) < 2.0 )
    {
      effAreaChHad = 0.0306;
      effAreaNHad  = 0.0629;
      effAreaPho   = 0.0699;
    }
  else if( fabs( eta ) < 2.2 )
    {
      effAreaChHad = 0.0283;
      effAreaNHad  = 0.0197;
      effAreaPho   = 0.1056;
    }
  else if( fabs( eta ) < 2.3 )
    {
      effAreaChHad = 0.0254;
      effAreaNHad  = 0.0184;
      effAreaPho   = 0.1457;
    }
  else if( fabs( eta ) < 2.4 )
    {
      effAreaChHad = 0.0217;
      effAreaNHad  = 0.0284;
      effAreaPho   = 0.1719;
    }
  else
    {
      effAreaChHad = 0.0167;
      effAreaNHad  = 0.0591;
      effAreaPho   = 0.1998;
    }
};

//https://indico.cern.ch/event/732974/contributions/3071231/attachments/1685034/2709197/Photon_ID_Finalv2.pdf
void RazorAnalyzer::getPhotonEffAreaPFClusterIso( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
      if(fabs (eta) < 0.8 ) {
		effAreaChHad = 0.037;
		effAreaNHad = 0.089;
		effAreaPho = 0.19;
	}
      else {
		effAreaChHad = 0.031;
                effAreaNHad = 0.15;
                effAreaPho = 0.19;
	}
 };

//photon ID and isolation cuts from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
bool RazorAnalyzer::photonPassesIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut, bool useEffectiveArea90, bool usePrivatePF, bool usePFClusterIso){
    //get effective area for isolation calculations
    double effAreaChargedHadrons = 0.0;
    double effAreaNeutralHadrons = 0.0;
    double effAreaPhotons = 0.0;

    double effAreaChargedHadrons_ClusterIso = 0.0;
    double effAreaNeutralHadrons_ClusterIso = 0.0;
    double effAreaPhotons_ClusterIso = 0.0;

    //get the effective areas. results are passed to variables by reference
    if(usePFClusterIso)
      {
        getPhotonEffAreaPFClusterIso( pho_superClusterEta[i] , effAreaChargedHadrons_ClusterIso, effAreaNeutralHadrons_ClusterIso, effAreaPhotons_ClusterIso);
      }
    if (useEffectiveArea90)
      {
	getPhotonEffArea90( pho_superClusterEta[i] , effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons);
      }
    else
      {
	getPhotonEffAreaRun2( pho_superClusterEta[i] , effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons);
      }

    //Rho corrected PF charged hadron isolation
    double PFIsoCorrected_ChHad = fmax(pho_pfIsoChargedHadronIso[i] - fixedGridRhoFastjetAll*effAreaChargedHadrons, 0.);
    if(usePrivatePF) PFIsoCorrected_ChHad = fmax(pho_sumChargedHadronPt[i] - fixedGridRhoFastjetAll*effAreaChargedHadrons, 0.);
    if(usePFClusterIso) PFIsoCorrected_ChHad = fmax(pho_trkSumPtHollowConeDR03[i] - fixedGridRhoFastjetAll*effAreaChargedHadrons_ClusterIso, 0.);
    if(PFIsoCorrected_ChHad > PFChHadIsoCut) return false;

    //Rho corrected PF neutral hadron isolation
    double PFIsoCorrected_NeuHad = fmax(pho_pfIsoNeutralHadronIso[i] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.);
    if(usePrivatePF) PFIsoCorrected_NeuHad = fmax(pho_sumNeutralHadronEt[i] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.);
    if(usePFClusterIso) PFIsoCorrected_NeuHad = fmax(pho_sumNeutralHadronEt[i] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.);
    if(PFIsoCorrected_NeuHad > PFNeuHadIsoCut) return false;

    //Rho corrected PF photon isolation
    double PFIsoCorrected_Photons = fmax(pho_pfIsoPhotonIso[i] - fixedGridRhoFastjetAll*effAreaPhotons, 0.);
    if(usePrivatePF) PFIsoCorrected_Photons = fmax(pho_sumPhotonEt[i] - fixedGridRhoFastjetAll*effAreaPhotons, 0.);
    if(usePFClusterIso) PFIsoCorrected_Photons = fmax(pho_ecalPFClusterIso[i] - fixedGridRhoFastjetAll*effAreaPhotons_ClusterIso, 0.);
    if(PFIsoCorrected_Photons > PFPhotIsoCut) return false;

    //photon passed all cuts
    return true;
}


// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
bool RazorAnalyzer::photonPassLooseIDWithoutEleVeto(int i, bool use25nsCuts ){

  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(pho_HoverE[i] > 0.0597) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.01031) pass = false;
    } else {
      if(pho_HoverE[i] > 0.0481) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.03013) pass = false;
    }
  }  else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  return pass;
}

// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
bool RazorAnalyzer::photonPassMediumIDWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(pho_HoverE[i] > 0.0396) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.01022) pass = false;
    } else {
      if(pho_HoverE[i] > 0.0219) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.03001) pass = false;
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  return pass;
}

// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
bool RazorAnalyzer::photonPassTightIDWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(pho_HoverE[i] > 0.0269) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.00994) pass = false;
    } else {
      if(pho_HoverE[i] > 0.0213) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.03000) pass = false;
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  return pass;
}


// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
bool RazorAnalyzer::photonPassLooseDelayedIDWithoutEleVeto(int i, bool use25nsCuts ){

  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(phoFull5x5SigmaIetaIeta[i] > 0.01031) pass = false;
    } else {
      if(phoFull5x5SigmaIetaIeta[i] > 0.03013) pass = false;
    }
  }  else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  //double phoSminor = getPhotonSminorSmajor(i, true);
  //if(phoSminor < 0.15 || phoSminor > 0.7) pass = false;

  return pass;
}

// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
bool RazorAnalyzer::photonPassMediumDelayedIDWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(phoFull5x5SigmaIetaIeta[i] > 0.01022) pass = false;
    } else {
      if(phoFull5x5SigmaIetaIeta[i] > 0.03001) pass = false;
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  //double phoSminor = getPhotonSminorSmajor(i, true);
  //if(phoSminor < 0.15 || phoSminor > 0.5) pass = false;

  return pass;
}

// 80X values from EGamma Presentation
// https://indico.cern.ch/event/491517/contributions/2349134/attachments/1359450/2056689/CutBasedPhotonID_24-10-2016.pdf
bool RazorAnalyzer::photonPassTightDelayedIDWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(phoFull5x5SigmaIetaIeta[i] > 0.00994) pass = false;
    } else {
      if(phoFull5x5SigmaIetaIeta[i] > 0.03000) pass = false;
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  //double phoSminor = getPhotonSminorSmajor(i, true);
  //if(phoSminor < 0.15 || phoSminor > 0.3) pass = false;

  return pass;
}


bool RazorAnalyzer::photonPassLooseID(int i, bool use25nsCuts){

  bool pass = true;

  if (!photonPassLooseIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::photonPassMediumID(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassMediumIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::photonPassTightID(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassTightIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

// 80X-v2.2 Cuts from EGamma Presentation
// https://indico.cern.ch/event/491548/contributions/2384977/attachments/1377936/2093213/CutBasedPhotonID_25-11-2016.pdf
bool RazorAnalyzer::photonPassLooseIso(int i, bool use25nsCuts, bool usePrivatePF, bool usePFClusterIso){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(!usePFClusterIso) return photonPassesIsolation(i, 1.295, 10.910 + 0.0148*phoPt[i] + 0.000017*phoPt[i]*phoPt[i], 3.630 + 0.0047*phoPt[i], true, usePrivatePF, usePFClusterIso );
      else return photonPassesIsolation(i, 8.5 + 0.00091*phoPt[i], 10.910 + 0.0148*phoPt[i] + 0.000017*phoPt[i]*phoPt[i], 8.0+0.00092*phoPt[i], true, usePrivatePF, usePFClusterIso );
    } else {
      return photonPassesIsolation(i, 1.011, 5.931 + 0.0163*phoPt[i] + 0.000014*phoPt[i]*phoPt[i], 6.641 + 0.0034*phoPt[i], true, usePrivatePF, usePFClusterIso);
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    return false;
  }

}

// 80X-v2.2 Cuts from EGamma Presentation
// https://indico.cern.ch/event/491548/contributions/2384977/attachments/1377936/2093213/CutBasedPhotonID_25-11-2016.pdf
//bool RazorAnalyzer::photonPassMediumIso(int i, bool use25nsCuts, bool usePrivatePF, bool usePFClusterIso){
bool RazorAnalyzer::photonPassMediumIso(int i, bool use25nsCuts, bool usePrivatePF){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      return photonPassesIsolation(i, 0.441, 2.725 + 0.0148*phoPt[i] + 0.000017*phoPt[i]*phoPt[i], 2.571 + 0.0047*phoPt[i], true, usePrivatePF);
    } else {
      return photonPassesIsolation(i, 0.442, 1.715 + 0.0163*phoPt[i] + 0.000014*phoPt[i]*phoPt[i], 3.863 + 0.0034*phoPt[i], true, usePrivatePF);
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    return false;
  }

}

// 80X-v2.2 Cuts from EGamma Presentation
// https://indico.cern.ch/event/491548/contributions/2384977/attachments/1377936/2093213/CutBasedPhotonID_25-11-2016.pdf
// for OOT photon: https://indico.cern.ch/event/732974/contributions/3071231/attachments/1685034/2709197/Photon_ID_Finalv2.pdf
bool RazorAnalyzer::photonPassTightIso(int i, bool use25nsCuts, bool usePrivatePF, bool usePFClusterIso){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(!usePFClusterIso)  return photonPassesIsolation(i, 0.202, 0.264 + 0.0148*phoPt[i] + 0.000017*phoPt[i]*phoPt[i], 2.362 + 0.0047*phoPt[i], true, usePrivatePF);
      else return photonPassesIsolation(i, 5.5 + 0.00094*phoPt[i], 0.264 + 0.0148*phoPt[i] + 0.000017*phoPt[i]*phoPt[i], 5.0+0.00092*phoPt[i], true, usePrivatePF, true );
    } else {
      return photonPassesIsolation(i, 0.034, 0.586 + 0.0163*phoPt[i] + 0.000014*phoPt[i]*phoPt[i], 2.617 + 0.0034*phoPt[i], true, usePrivatePF);
    }
  }  else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    return false;
  }

}


bool RazorAnalyzer::isLoosePhoton(int i, bool use25nsCuts){

  bool pass = true;
  if(!isLoosePhotonWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::isMediumPhoton(int i, bool use25nsCuts){
  bool pass = true;

  if(!isMediumPhotonWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::isTightPhoton(int i, bool use25nsCuts){
  bool pass = true;
  if (!isTightPhotonWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::isLoosePhotonWithoutEleVeto(int i, bool use25nsCuts){

  bool pass = true;

  if (!photonPassLooseIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassLooseIso(i,use25nsCuts)) pass = false;

  return pass;
}

bool RazorAnalyzer::isMediumPhotonWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassMediumIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassMediumIso(i,use25nsCuts)) pass = false;

  return pass;
}

bool RazorAnalyzer::isTightPhotonWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassTightIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassTightIso(i,use25nsCuts)) pass = false;

  return pass;
}


bool RazorAnalyzer::isLooseDelayedPhotonWithoutEleVeto(int i, bool use25nsCuts){

  bool pass = true;

  if (!photonPassLooseDelayedIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassLooseIso(i,use25nsCuts,false,true)) pass = false;

  return pass;
}

bool RazorAnalyzer::isMediumDelayedPhotonWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassMediumDelayedIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassMediumIso(i,use25nsCuts,false)) pass = false;

  return pass;
}

bool RazorAnalyzer::isTightDelayedPhotonWithoutEleVeto(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassTightDelayedIDWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassTightIso(i,use25nsCuts,false,true)) pass = false;

  return pass;
}


bool RazorAnalyzer::matchPhotonHLTFilters(int i, string HLTFilter){
  bool match = false;

  if (HLTFilter == "DiPhoton30_18_WithPixMatch_Leg1") {
    if (
	//Data filters
	pho_passHLTFilter[i][8]
	//MC filters

	 ) {
      match = true;
    }
  }

  if (HLTFilter == "DiPhoton30_18_WithPixMatch_Leg2") {
    if (
	//Data filters
	pho_passHLTFilter[i][9]
	 ) {
      match = true;
    }
  }


  if (HLTFilter == "Photon42_Photon25_Mass15_Leg1") {
    if (
	//Data filters
	pho_passHLTFilter[i][4]
	//MC filters

	 ) {
      match = true;
    }
  }

  if (HLTFilter == "Photon42_Photon25_Mass15_Leg2") {
    if (
	//Data filters
	pho_passHLTFilter[i][5]
	 ) {
      match = true;
    }
  }


  return match;
}

TLorentzVector RazorAnalyzer::GetCorrectedMomentum( TVector3 vtx, TVector3 phoPos, double phoE )
{
  TVector3 phoDir = phoPos - vtx;
  TVector3 phoP3  = phoDir.Unit()*phoE;
  return TLorentzVector( phoP3, phoE);
};

bool RazorAnalyzer::photonPassLooseIDWithoutEleVetoExo15004( int i )
{
  bool pass = true;
  if ( fabs(pho_superClusterEta[i]) < 1.4442 )
    {
      if( pho_HoverE[i] > 0.05 ) pass = false;
      if( phoFull5x5SigmaIetaIeta[i] > 0.0105 ) pass = false;
    }
  else
    {
      if( pho_HoverE[i] > 0.05 ) pass = false;
      if( phoFull5x5SigmaIetaIeta[i] > 0.0280 ) pass = false;
    }
  return pass;
};

bool RazorAnalyzer::photonPassesIsolationExo15004(int i, double PFChHadIsoCut, double PFPhotIsoCut )
{
  double effAreaPhotons = 0.0;
  double eta = pho_superClusterEta[i];
  getPhotonEffAreaExo15004( eta, effAreaPhotons );

  //Rho corrected PF charged hadron isolation
  //double PFIsoCorrected_ChHad = pho_sumChargedHadronPt[i];//No PU correction (Caltech Original)
  double PFIsoCorrected_ChHad = pho_pfIsoChargedHadronIso[i];//(Exo15004 default pfIso)
  if(PFIsoCorrected_ChHad > PFChHadIsoCut) return false;

  //Rho corrected PF photon isolation
  //double PFIsoCorrected_Photons = pho_sumPhotonEt[i] - fixedGridRhoFastjetAll*effAreaPhotons;//PU corr can go neg!(Caltech Original)
  double PFIsoCorrected_Photons = pho_pfIsoPhotonIso[i] - fixedGridRhoAll*effAreaPhotons;//PU corr can go neg!(Exo15004 default pfIso)
  if(PFIsoCorrected_Photons > PFPhotIsoCut) return false;

  //photon passed all cuts
  return true;
};

void RazorAnalyzer::getPhotonEffAreaExo15004( float eta, double& effAreaPho )
{
  if ( fabs( eta ) < 0.9 )
    {
      effAreaPho = 0.17;
    }
  else if ( fabs( eta ) < 1.4442 )
    {
      effAreaPho = 0.14;
    }
  else if ( fabs( eta ) < 2.0 )
    {
      effAreaPho = 0.11;
    }
  else if ( fabs( eta ) < 2.2 )
    {
      effAreaPho = 0.14;
    }
  else
    {
      effAreaPho = 0.22;
    }
};


bool RazorAnalyzer::photonPassLooseIsoExo15004(int i)
{
  if( fabs(pho_superClusterEta[i]) < 1.4442 )
    {
      return photonPassesIsolationExo15004(i, 5, (2.75 - 2.5) + 0.0045*phoPt[i] );
    }
  else if ( fabs(pho_superClusterEta[i]) < 2.0 )
    {
      return photonPassesIsolationExo15004(i, 5, (2.0 - 2.5) + 0.003*phoPt[i] );
    }
  else
    {
      return photonPassesIsolationExo15004(i, 5, (2.0 - 2.5) + 0.003*phoPt[i] );
    }
};

//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////Photon ID 2017 92X///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


// 92X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_92X_and_later
void RazorAnalyzer::getPhotonEffArea90_2017( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
  if( fabs( eta ) < 1.0 )
    {
      effAreaChHad = 0.0385;
      effAreaNHad  = 0.0636;
      effAreaPho   = 0.1240;
    }
  else if( fabs( eta ) < 1.479 )
    {
      effAreaChHad = 0.0468;
      effAreaNHad  = 0.1103;
      effAreaPho   = 0.1093;
    }
  else if( fabs( eta ) < 2.0 )
    {
      effAreaChHad = 0.0435;
      effAreaNHad  = 0.0759;
      effAreaPho   = 0.0631;
    }
  else if( fabs( eta ) < 2.2 )
    {
      effAreaChHad = 0.0378;
      effAreaNHad  = 0.0236;
      effAreaPho   = 0.0779;
    }
  else if( fabs( eta ) < 2.3 )
    {
      effAreaChHad = 0.0338;
      effAreaNHad  = 0.0151;
      effAreaPho   = 0.0999;
    }
  else if( fabs( eta ) < 2.4 )
    {
      effAreaChHad = 0.0314;
      effAreaNHad  = 0.00007;
      effAreaPho   = 0.1155;
    }
  else
    {
      effAreaChHad = 0.0269;
      effAreaNHad  = 0.0132;
      effAreaPho   = 0.1373;
    }
};

void RazorAnalyzer::getPhotonEffAreaRun2_2017( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
  if( fabs( eta ) < 1.0 )
    {
      effAreaChHad = 0.0385;
      effAreaNHad  = 0.0636;
      effAreaPho   = 0.1240;
    }
  else if( fabs( eta ) < 1.479 )
    {
      effAreaChHad = 0.0468;
      effAreaNHad  = 0.1103;
      effAreaPho   = 0.1093;
    }
  else if( fabs( eta ) < 2.0 )
    {
      effAreaChHad = 0.0435;
      effAreaNHad  = 0.0759;
      effAreaPho   = 0.0631;
    }
  else if( fabs( eta ) < 2.2 )
    {
      effAreaChHad = 0.0378;
      effAreaNHad  = 0.0236;
      effAreaPho   = 0.0779;
    }
  else if( fabs( eta ) < 2.3 )
    {
      effAreaChHad = 0.0338;
      effAreaNHad  = 0.0151;
      effAreaPho   = 0.0999;
    }
  else if( fabs( eta ) < 2.4 )
    {
      effAreaChHad = 0.0314;
      effAreaNHad  = 0.00007;
      effAreaPho   = 0.1155;
    }
  else
    {
      effAreaChHad = 0.0269;
      effAreaNHad  = 0.0132;
      effAreaPho   = 0.1373;
    }
};


//photon ID and isolation cuts from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
bool RazorAnalyzer::photonPassesIsolation_2017(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut, bool useEffectiveArea90, bool usePrivatePF){
    //get effective area for isolation calculations
    double effAreaChargedHadrons = 0.0;
    double effAreaNeutralHadrons = 0.0;
    double effAreaPhotons = 0.0;

    //get the effective areas. results are passed to variables by reference
    if (useEffectiveArea90)
      {
	getPhotonEffArea90_2017( pho_superClusterEta[i] , effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons);
      }
    else
      {
	getPhotonEffAreaRun2_2017( pho_superClusterEta[i] , effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons);
      }

    //Rho corrected PF charged hadron isolation
    double PFIsoCorrected_ChHad = fmax(pho_pfIsoChargedHadronIso[i] - fixedGridRhoFastjetAll*effAreaChargedHadrons, 0.);
    if(usePrivatePF) PFIsoCorrected_ChHad = fmax(pho_sumChargedHadronPt[i] - fixedGridRhoFastjetAll*effAreaChargedHadrons, 0.);
    if(PFIsoCorrected_ChHad > PFChHadIsoCut) return false;

    //Rho corrected PF neutral hadron isolation
    double PFIsoCorrected_NeuHad = fmax(pho_pfIsoNeutralHadronIso[i] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.);
    if(usePrivatePF) PFIsoCorrected_NeuHad = fmax(pho_sumNeutralHadronEt[i] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.);
    if(PFIsoCorrected_NeuHad > PFNeuHadIsoCut) return false;

    //Rho corrected PF photon isolation
    double PFIsoCorrected_Photons = fmax(pho_pfIsoPhotonIso[i] - fixedGridRhoFastjetAll*effAreaPhotons, 0.);
    if(usePrivatePF) PFIsoCorrected_Photons = fmax(pho_sumPhotonEt[i] - fixedGridRhoFastjetAll*effAreaPhotons, 0.);
    if(PFIsoCorrected_Photons > PFPhotIsoCut) return false;

    //photon passed all cuts
    return true;
}


// 92X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_92X_and_later
bool RazorAnalyzer::photonPassLooseIDWithoutEleVeto_2017(int i, bool use25nsCuts ){

  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(pho_HoverE[i] > 0.105) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0103) pass = false;
    } else {
      if(pho_HoverE[i] > 0.029) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0276) pass = false;
    }
  }  else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  return pass;
}

// 92X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_92X_and_later
bool RazorAnalyzer::photonPassMediumIDWithoutEleVeto_2017(int i, bool use25nsCuts){
  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(pho_HoverE[i] > 0.035) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0103) pass = false;
    } else {
      if(pho_HoverE[i] > 0.027) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0271) pass = false;
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  return pass;
}

// 92X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_92X_and_later
bool RazorAnalyzer::photonPassTightIDWithoutEleVeto_2017(int i, bool use25nsCuts){
  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(pho_HoverE[i] > 0.020) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0103) pass = false;
    } else {
      if(pho_HoverE[i] > 0.025) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0271) pass = false;
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  return pass;
}

bool RazorAnalyzer::photonPassLooseID_2017(int i, bool use25nsCuts){

  bool pass = true;

  if (!photonPassLooseIDWithoutEleVeto_2017(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::photonPassMediumID_2017(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassMediumIDWithoutEleVeto_2017(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::photonPassTightID_2017(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassTightIDWithoutEleVeto_2017(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

// 92X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_92X_and_later
bool RazorAnalyzer::photonPassLooseIso_2017(int i, bool use25nsCuts, bool usePrivatePF){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      return photonPassesIsolation_2017(i, 2.839, 9.188 + 0.0126*phoPt[i] + 0.000026*phoPt[i]*phoPt[i], 2.956 + 0.0035*phoPt[i], true, usePrivatePF );
    } else {
      return photonPassesIsolation_2017(i, 2.150, 10.471 + 0.0119*phoPt[i] + 0.000025*phoPt[i]*phoPt[i], 4.895 + 0.0040*phoPt[i], true, usePrivatePF);
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    return false;
  }

}

// 92X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_92X_and_later
bool RazorAnalyzer::photonPassMediumIso_2017(int i, bool use25nsCuts, bool usePrivatePF){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      return photonPassesIsolation_2017(i, 1.416, 2.491 + 0.0126*phoPt[i] + 0.000026*phoPt[i]*phoPt[i], 2.952 + 0.0040*phoPt[i], true, usePrivatePF);
    } else {
      return photonPassesIsolation_2017(i, 1.012, 9.131 + 0.0119*phoPt[i] + 0.000025*phoPt[i]*phoPt[i], 4.095 + 0.0040*phoPt[i], true, usePrivatePF);
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    return false;
  }

}

// 92X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_92X_and_later
bool RazorAnalyzer::photonPassTightIso_2017(int i, bool use25nsCuts, bool usePrivatePF){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation_2017(i, 1.158, 1.267 + 0.0126*phoPt[i] + 0.000026*phoPt[i]*phoPt[i], 2.065 + 0.0035*phoPt[i], true, usePrivatePF);
    } else {
      return photonPassesIsolation_2017(i, 0.575, 8.916 + 0.0119*phoPt[i] + 0.000025*phoPt[i]*phoPt[i], 3.272 + 0.0040*phoPt[i], true, usePrivatePF);
    }
  }  else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    return false;
  }

}


bool RazorAnalyzer::isLoosePhotonWithoutEleVeto_2017(int i, bool use25nsCuts){

  bool pass = true;

  if (!photonPassLooseIDWithoutEleVeto_2017(i,use25nsCuts)) pass = false;
  if(!photonPassLooseIso_2017(i,use25nsCuts)) pass = false;

  return pass;
}

bool RazorAnalyzer::isMediumPhotonWithoutEleVeto_2017(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassMediumIDWithoutEleVeto_2017(i,use25nsCuts)) pass = false;
  if(!photonPassMediumIso_2017(i,use25nsCuts)) pass = false;

  return pass;
}

bool RazorAnalyzer::isTightPhotonWithoutEleVeto_2017(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassTightIDWithoutEleVeto_2017(i,use25nsCuts)) pass = false;
  if(!photonPassTightIso_2017(i,use25nsCuts)) pass = false;

  return pass;
}

bool RazorAnalyzer::isLoosePhoton_2017(int i, bool use25nsCuts){

  bool pass = true;
  if(!isLoosePhotonWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}


bool RazorAnalyzer::isMediumPhoton_2017(int i, bool use25nsCuts){
  bool pass = true;

  if(!isMediumPhotonWithoutEleVeto_2017(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::isTightPhoton_2017(int i, bool use25nsCuts){
  bool pass = true;
  if (!isTightPhotonWithoutEleVeto_2017(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}


//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////Photon ID 2017 94X///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


// 94X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
void RazorAnalyzer::getPhotonEffArea90_94X( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
  if( fabs( eta ) < 1.0 )
    {
      effAreaChHad = 0.0112;
      effAreaNHad  = 0.0668;
      effAreaPho   = 0.1113;
    }
  else if( fabs( eta ) < 1.479 )
    {
      effAreaChHad = 0.0108;
      effAreaNHad  = 0.1054;
      effAreaPho   = 0.0953;
    }
  else if( fabs( eta ) < 2.0 )
    {
      effAreaChHad = 0.0106;
      effAreaNHad  = 0.0786;
      effAreaPho   = 0.0619;
    }
  else if( fabs( eta ) < 2.2 )
    {
      effAreaChHad = 0.01002;
      effAreaNHad  = 0.0233;
      effAreaPho   = 0.0837;
    }
  else if( fabs( eta ) < 2.3 )
    {
      effAreaChHad = 0.0098;
      effAreaNHad  = 0.0078;
      effAreaPho   = 0.1070;
    }
  else if( fabs( eta ) < 2.4 )
    {
      effAreaChHad = 0.0089;
      effAreaNHad  = 0.0028;
      effAreaPho   = 0.1212;
    }
  else
    {
      effAreaChHad = 0.0087;
      effAreaNHad  = 0.0137;
      effAreaPho   = 0.1466;
    }
};

void RazorAnalyzer::getPhotonEffAreaRun2_94X( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho )
{
  if( fabs( eta ) < 1.0 )
    {
      effAreaChHad = 0.0112;
      effAreaNHad  = 0.0668;
      effAreaPho   = 0.1113;
    }
  else if( fabs( eta ) < 1.479 )
    {
      effAreaChHad = 0.0108;
      effAreaNHad  = 0.1054;
      effAreaPho   = 0.0953;
    }
  else if( fabs( eta ) < 2.0 )
    {
      effAreaChHad = 0.0106;
      effAreaNHad  = 0.0786;
      effAreaPho   = 0.0619;
    }
  else if( fabs( eta ) < 2.2 )
    {
      effAreaChHad = 0.01002;
      effAreaNHad  = 0.0233;
      effAreaPho   = 0.0837;
    }
  else if( fabs( eta ) < 2.3 )
    {
      effAreaChHad = 0.0098;
      effAreaNHad  = 0.0078;
      effAreaPho   = 0.1070;
    }
  else if( fabs( eta ) < 2.4 )
    {
      effAreaChHad = 0.0089;
      effAreaNHad  = 0.0028;
      effAreaPho   = 0.1212;
    }
  else
    {
      effAreaChHad = 0.0087;
      effAreaNHad  = 0.0137;
      effAreaPho   = 0.1466;
    }
};


//photon ID and isolation cuts from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
bool RazorAnalyzer::photonPassesIsolation_94X(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut, bool useEffectiveArea90, bool usePrivatePF){
    //get effective area for isolation calculations
    double effAreaChargedHadrons = 0.0;
    double effAreaNeutralHadrons = 0.0;
    double effAreaPhotons = 0.0;

    //get the effective areas. results are passed to variables by reference
    if (useEffectiveArea90)
      {
	getPhotonEffArea90_94X( pho_superClusterEta[i] , effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons);
      }
    else
      {
	getPhotonEffAreaRun2_94X( pho_superClusterEta[i] , effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons);
      }

    //Rho corrected PF charged hadron isolation
    double PFIsoCorrected_ChHad = fmax(pho_pfIsoChargedHadronIso[i] - fixedGridRhoFastjetAll*effAreaChargedHadrons, 0.);
    if(usePrivatePF) PFIsoCorrected_ChHad = fmax(pho_sumChargedHadronPt[i] - fixedGridRhoFastjetAll*effAreaChargedHadrons, 0.);
    if(PFIsoCorrected_ChHad > PFChHadIsoCut) return false;

    //Rho corrected PF neutral hadron isolation
    double PFIsoCorrected_NeuHad = fmax(pho_pfIsoNeutralHadronIso[i] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.);
    if(usePrivatePF) PFIsoCorrected_NeuHad = fmax(pho_sumNeutralHadronEt[i] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.);
    if(PFIsoCorrected_NeuHad > PFNeuHadIsoCut) return false;

    //Rho corrected PF photon isolation
    double PFIsoCorrected_Photons = fmax(pho_pfIsoPhotonIso[i] - fixedGridRhoFastjetAll*effAreaPhotons, 0.);
    if(usePrivatePF) PFIsoCorrected_Photons = fmax(pho_sumPhotonEt[i] - fixedGridRhoFastjetAll*effAreaPhotons, 0.);
    if(PFIsoCorrected_Photons > PFPhotIsoCut) return false;

    //photon passed all cuts
    return true;
}


// 94X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
bool RazorAnalyzer::photonPassLooseIDWithoutEleVeto_94X(int i, bool use25nsCuts ){

  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(pho_HoverE[i] > 0.04596 ) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0106 ) pass = false;
    } else {
      if(pho_HoverE[i] > 0.0590 ) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0272) pass = false;
    }
  }  else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  return pass;
}

// 94X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
bool RazorAnalyzer::photonPassMediumIDWithoutEleVeto_94X(int i, bool use25nsCuts){
  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(pho_HoverE[i] > 0.02197 ) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.01015 ) pass = false;
    } else {
      if(pho_HoverE[i] > 0.0326 ) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0272) pass = false;
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  return pass;
}

// 94X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
bool RazorAnalyzer::photonPassTightIDWithoutEleVeto_94X(int i, bool use25nsCuts){
  bool pass = true;

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      if(pho_HoverE[i] > 0.02148 ) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.00996 ) pass = false;
    } else {
      if(pho_HoverE[i] > 0.0321 ) pass = false;
      if(phoFull5x5SigmaIetaIeta[i] > 0.0271 ) pass = false;
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    pass = false;
  }

  return pass;
}

bool RazorAnalyzer::photonPassLooseID_94X(int i, bool use25nsCuts){

  bool pass = true;

  if (!photonPassLooseIDWithoutEleVeto_94X(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::photonPassMediumID_94X(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassMediumIDWithoutEleVeto_94X(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::photonPassTightID_94X(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassTightIDWithoutEleVeto_94X(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

// 94X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
bool RazorAnalyzer::photonPassLooseIso_94X(int i, bool use25nsCuts, bool usePrivatePF){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      return photonPassesIsolation_94X(i, 1.694, 24.032 + 0.01512*phoPt[i] + 2.259e-05*phoPt[i]*phoPt[i], 2.876 + 0.004017*phoPt[i], true, usePrivatePF );
    } else {
      return photonPassesIsolation_94X(i, 2.089, 19.722 + 0.0117*phoPt[i] + 2.3e-05*phoPt[i]*phoPt[i], 4.162 + 0.0037*phoPt[i], true, usePrivatePF);
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    return false;
  }

}

// 94X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
bool RazorAnalyzer::photonPassMediumIso_94X(int i, bool use25nsCuts, bool usePrivatePF){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
      return photonPassesIsolation_94X(i, 1.141, 1.189 + 0.01512*phoPt[i] + 2.259e-05*phoPt[i]*phoPt[i], 2.08 + 0.004017*phoPt[i], true, usePrivatePF);
    } else {
      return photonPassesIsolation_94X(i, 1.051, 2.718 + 0.0117*phoPt[i] + 2.3e-05*phoPt[i]*phoPt[i], 3.867 + 0.0037*phoPt[i], true, usePrivatePF);
    }
  } else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    return false;
  }

}

// 94X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
bool RazorAnalyzer::photonPassTightIso_94X(int i, bool use25nsCuts, bool usePrivatePF){

  if (use25nsCuts) {
    if(fabs(pho_superClusterEta[i]) < 1.479){
        return photonPassesIsolation_94X(i, 0.65, 0.317 + 0.01512*phoPt[i] + 2.259e-05*phoPt[i]*phoPt[i], 2.044 + 0.004017*phoPt[i], true, usePrivatePF);
    } else {
      return photonPassesIsolation_94X(i, 0.517, 2.716 + 0.0117*phoPt[i] + 2.3e-05*phoPt[i]*phoPt[i], 3.032 + 0.0037*phoPt[i], true, usePrivatePF);
    }
  }  else {
    cout << "Warning: you are not using 25nsCuts. return false.\n";
    return false;
  }

}


bool RazorAnalyzer::isLoosePhotonWithoutEleVeto_94X(int i, bool use25nsCuts){

  bool pass = true;

  if (!photonPassLooseIDWithoutEleVeto_94X(i,use25nsCuts)) pass = false;
  if(!photonPassLooseIso_94X(i,use25nsCuts)) pass = false;

  return pass;
}

bool RazorAnalyzer::isMediumPhotonWithoutEleVeto_94X(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassMediumIDWithoutEleVeto_94X(i,use25nsCuts)) pass = false;
  if(!photonPassMediumIso_94X(i,use25nsCuts)) pass = false;

  return pass;
}

bool RazorAnalyzer::isTightPhotonWithoutEleVeto_94X(int i, bool use25nsCuts){
  bool pass = true;

  if (!photonPassTightIDWithoutEleVeto_94X(i,use25nsCuts)) pass = false;
  if(!photonPassTightIso_94X(i,use25nsCuts)) pass = false;

  return pass;
}

bool RazorAnalyzer::isLoosePhoton_94X(int i, bool use25nsCuts){

  bool pass = true;
  if(!isLoosePhotonWithoutEleVeto(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}


bool RazorAnalyzer::isMediumPhoton_94X(int i, bool use25nsCuts){
  bool pass = true;

  if(!isMediumPhotonWithoutEleVeto_94X(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

bool RazorAnalyzer::isTightPhoton_94X(int i, bool use25nsCuts){
  bool pass = true;
  if (!isTightPhotonWithoutEleVeto_94X(i,use25nsCuts)) pass = false;
  if(!photonPassesElectronVeto(i)) pass = false;

  return pass;
}

/*
// 92X values from EGamma twiki
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_92X_and_later
void RazorAnalyzer::getPhotonEffAreaExo15004( float eta, double& effAreaPho )
{
  if ( fabs( eta ) < 1.0 )
    {
      effAreaPho = 0.1240;
    }
  else if ( fabs( eta ) < 1.479 )
    {
      effAreaPho = 0.1093;
    }
  else if ( fabs( eta ) < 2.0 )
    {
      effAreaPho = 0.0631;
    }
  else if ( fabs( eta ) < 2.2 )
    {
      effAreaPho = 0.0779;
    }
  else if ( fabs( eta ) < 2.3 )
    {
      effAreaPho = 0.0999;
    }
  else if ( fabs( eta ) < 2.4 )
    {
      effAreaPho = 0.1155;
    }
  else
    {
      effAreaPho = 0.1373;
    }
};
*/

//////////////////////////////
// GEN
//////////////////////////////

//Finds closes gen electron and returns index to gParticle
int RazorAnalyzer::findClosestGenElectron(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nGenParticle; j++){
    if (gParticleStatus[j] != 1) continue;
    if (abs(gParticleId[j]) != 11) continue;
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1
	 && deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < minDR
	 ) {
      matchedIndex = j;
      minDR = deltaR( eta, phi, gParticleEta[j], gParticlePhi[j]);
    }
  }

  return matchedIndex;
};


//Finds closes gen muon and returns index to gParticle
int RazorAnalyzer::findClosestGenMuon(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nGenParticle; j++){
    if (gParticleStatus[j] != 1) continue;
    if (abs(gParticleId[j]) != 13) continue;
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1
	 && deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < minDR
	 ) {
      matchedIndex = j;
      minDR = deltaR( eta, phi, gParticleEta[j], gParticlePhi[j]);
    }
  }

  return matchedIndex;
};

//Finds closest gen jet and returns index to gParticle
int RazorAnalyzer::findClosestGenJet(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nGenJets; j++){
    //    if (gParticleStatus[j] != 1) continue;
    //if (abs(gParticleId[j]) != 13) continue;

    if ( deltaR(eta, phi, genJetEta[j], genJetPhi[j]) < 0.3
         && deltaR(eta, phi, genJetEta[j], genJetPhi[j]) < minDR
         ) {
      matchedIndex = j;
      minDR = deltaR( eta, phi, genJetEta[j], genJetPhi[j]);
    }

  }

  return matchedIndex;
};



//Checks if the gParticle is a tau and that comes from a W or a Z
bool RazorAnalyzer::isGenTau(int index){
  return ( abs(gParticleId[index]) == 15 && gParticleStatus[index] == 2 &&
	   (abs(gParticleMotherId[index]) == 23 ||abs(gParticleMotherId[index]) == 24)
	   );
};


//Checks if the gParticle is a tau and that comes from a W or a Z
bool RazorAnalyzer::isGenLeptonicTau(int index){
  if (abs(gParticleId[index]) == 15 && gParticleStatus[index] == 2
      && (abs(gParticleMotherId[index]) == 24 || abs(gParticleMotherId[index]) == 23)
      ) {

    for(int k = 0; k < nGenParticle; k++){
      if ( (abs(gParticleId[k]) == 11 || abs(gParticleId[k]) == 13) && gParticleMotherIndex[k] == index) {
	return true;
      }
    }
  }
  return false;
};

//Finds closes gen tau and returns index to gParticle
int RazorAnalyzer::findClosestGenTau(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nGenParticle; j++){
    if (abs(gParticleId[j]) != 15) continue;
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1
	 && deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < minDR
	 ) {
      matchedIndex = j;
      minDR = deltaR( eta, phi, gParticleEta[j], gParticlePhi[j]);
    }
  }

  return matchedIndex;
};

int RazorAnalyzer::findClosestRecoTau(double eta, double phi) {
  int matchedIndex = -1;
  float minDR = 999.0;
  for(int j = 0; j < nTaus; j++){
    if ( deltaR(eta, phi, tauEta[j], tauPhi[j]) < 0.1
	 && deltaR(eta, phi, tauEta[j], tauPhi[j]) < minDR
	 ) {
      matchedIndex = j;
      minDR = deltaR(eta, phi, tauEta[j], tauPhi[j]);
    }
  }

  return matchedIndex;
};

//Finds closest gen tau and checks its mother, if matched return pdgID
int RazorAnalyzer::GetTauMatchedID(double eta, double phi){
  int matchedID = 0;
  int matchedIndex = findClosestGenTau(eta, phi);

  //find muon if no tau was found
  if (matchedIndex < 0) matchedIndex = findClosestGenMuon(eta,phi);
  if (matchedIndex >= 0) {
    return gParticleId[matchedIndex];
  }

  //find electron if no tau or muon was found
  if (matchedIndex < 0) matchedIndex = findClosestGenElectron(eta,phi);
  if (matchedIndex >= 0) {
    return gParticleId[matchedIndex];
  }

  //if nothing was found
  if(matchedIndex < 0) return 0;//No Match -> ID == 0
  if (gParticleMotherId[matchedIndex] > 50) {
    matchedID = gParticleMotherId[matchedIndex];
  } else if (abs(gParticleMotherId[matchedIndex]) == 23 ||
	     abs(gParticleMotherId[matchedIndex]) == 24) {
    matchedID = gParticleId[matchedIndex];
  }

  return matchedID;
};

//Returns index of the closest parton. If no match is found returns zero.
int RazorAnalyzer::findClosestParton(float eta, float phi){
  float minDRToParton = 9999;
  int partonIndex = -1;
  for(int j = 0; j < nGenParticle; j++){
    //only look for outgoing partons
    if  (!( ( (abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) || abs(gParticleId[j]) == 21)
	    && gParticleStatus[j] == 23)
	 ) continue;
    double tmpDR = deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]);
    if ( tmpDR < minDRToParton ) {
      minDRToParton = tmpDR;
      partonIndex = j;
    }
  }

  return partonIndex;
};


//Checks if a gen muon is in the given eta and phi direction
bool RazorAnalyzer::matchesGenMuon(double eta, double phi){
  bool result = false;
  for(int j = 0; j < nGenParticle; j++){
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1 &&
	 abs(gParticleId[j]) == 13 && gParticleStatus[j] == 1 &&
	 (abs(gParticleMotherId[j]) == 23 ||abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 15)
	 ) {
      result = true;
      break;
    }
  }
  return result;
};

//Checks if a gen electron is in the given eta and phi direction
bool RazorAnalyzer::matchesGenElectron(double eta, double phi){
  bool result = false;
  for(int j = 0; j < nGenParticle; j++){
    if ( deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]) < 0.1 &&
	 abs(gParticleId[j]) == 11 && gParticleStatus[j] == 1 &&
	 (abs(gParticleMotherId[j]) == 23 ||abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 15)
	 ) {
      result = true;
      break;
    }
  }
  return result;
};

// Returns true if the gen particle at the specified index
// decays hadronically into a quark with the specified status code.
bool RazorAnalyzer::isHadronicDecay(int index, int status) {
    for ( int j = 0; j < nGenParticle; j++ ) {
        if ( gParticleMotherIndex[j] == index && gParticleStatus[j] == status ) {
            if ( abs(gParticleId[j]) > 0 && abs(gParticleId[j]) < 5 ) return true;
        }
    }
    return false;
}

// Generic matching to hard process particles by deltaR
// eta, phi: coordinates of reco-level particle
// id: MC ID of the gen-level particle to match
// status: pythia status code of the gen-level particle (default 22)
// r: maximum deltaR needed for match
// returns the index of the matched particle in the gen particles collection,
// or -1 if no match
int RazorAnalyzer::getMatchingHardProcessParticleIndex(double eta, double phi,
        int id, int status, double r) {
    int matchedIndex = -1;
    float minDeltaR = -1;
    for ( int j = 0; j < nGenParticle; j++ ) {
        if (abs(gParticleId[j]) != abs(id)) continue;
        if (gParticleStatus[j] != status) continue;
        float dR = deltaR(eta, phi, gParticleEta[j], gParticlePhi[j]);
        if (dR > r) continue;
        if (minDeltaR < 0 || dR < minDeltaR) {
            minDeltaR = dR;
            matchedIndex = j;
        }
    }
    return matchedIndex;
}

// Gets index of matching gen W, if any
int RazorAnalyzer::getMatchingGenWIndex(double eta, double phi, double r) {
    int index = getMatchingHardProcessParticleIndex(eta, phi, 24, 22, r);
    return index;
}

// Gets index of matching gen top, if any
int RazorAnalyzer::getMatchingGenTopIndex(double eta, double phi, double r) {
    int index = getMatchingHardProcessParticleIndex(eta, phi, 6, 22, r);
    return index;
}


//Compute the genHT variable
double RazorAnalyzer::getGenHT(){
  double genHT = 0;
  for(int j = 0; j < nGenParticle; j++){
    if ( (gParticleStatus[j] == 23 || gParticleStatus[j] == 22)  &&
	 ( gParticleId[j] == 21 || ( abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) < 6 )) &&
	 ( gParticleMotherIndex[j] == -1 || (gParticleMotherIndex[j] >= 0 && gParticleStatus[gParticleMotherIndex[j]] == 21))
	 ) {
      genHT += gParticlePt[j];
    }
  }

  return genHT;
};

//Returns the number of gen-level ISR jets
int RazorAnalyzer::getNISR( std::vector<FactorizedJetCorrector*> &JetCorrector, std::vector<std::pair<int,int> > &JetCorrectorIOV ) {
    int NISRJets = 0;
    for(int i = 0; i < nJets; i++) {

        //Jet Corrections
        double JEC = JetEnergyCorrectionFactor( jetPt[i], jetEta[i], jetPhi[i], jetE[i],
                fixedGridRhoAll, jetJetArea[i], runNum,
                JetCorrectorIOV, JetCorrector );
        TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );

        //these are the cuts Ana/Manuel told me to use
        if ( thisJet.Pt() > 30 && fabs( thisJet.Eta()) < 2.4 && jetPassIDLoose[i]) {

            //try to match to gen partons
            //Follow prescription here: https://github.com/manuelfs/babymaker/blob/0136340602ee28caab14e3f6b064d1db81544a0a/bmaker/plugins/bmaker_full.cc#L1268-L1295
            bool match = false;
            for(int g = 0; g < nGenParticle; g++){

                double dR = deltaR( gParticleEta[g], gParticlePhi[g] , thisJet.Eta() , thisJet.Phi());
                if (dR > 0.3) continue;

                //check match against leptons
                if (abs(gParticleId[g]) == 11 || abs(gParticleId[g]) == 13 || abs(gParticleId[g]) == 15 ) {
                    match = true;
                }

                //check match against prompt photons
                if (abs(gParticleId[g]) == 22 &&
                        ( (gParticleStatus[g] == 1 && gParticleMotherId[g] != 22) || gParticleStatus[g] == 22 || gParticleStatus[g] == 23) &&
                        (abs(gParticleMotherId[g]) == 25 || abs(gParticleMotherId[g]) == 21 || abs(gParticleMotherId[g]) == 2212 ||
                         (abs(gParticleMotherId[g]) >= 1 && abs(gParticleMotherId[g]) <= 6) )
                   ) {
                    match = true;
                }

                //match to quarks
                if (gParticleStatus[g] == 23 && abs(gParticleId[g]) <= 5 &&
                        ( abs(gParticleMotherId[g]) == 6 ||  abs(gParticleMotherId[g]) == 23 ||  abs(gParticleMotherId[g]) == 24
                          ||  abs(gParticleMotherId[g]) == 25 ||  abs(gParticleMotherId[g]) > 1e6)) {
                    match = true;
                }
            }
            if (!match) NISRJets++;
        }
    }
    return NISRJets;
}


//////////////////////////////
// MISC
//////////////////////////////

double RazorAnalyzer::deltaPhi(double phi1, double phi2) {
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
    dphi -= TMath::TwoPi();
  while (dphi <= -TMath::Pi())
    dphi += TMath::TwoPi();
  return dphi;
}

double RazorAnalyzer::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1,phi2);
  double deta = eta1 - eta2;
  return sqrt( dphi*dphi + deta*deta);
}

TLorentzVector RazorAnalyzer::makeTLorentzVector(double pt, double eta, double phi, double energy){
    TLorentzVector vec;
    vec.SetPtEtaPhiE(pt, eta, phi, energy);
    return vec;
}

TLorentzVector RazorAnalyzer::makeTLorentzVectorPtEtaPhiM(double pt, double eta, double phi, double mass){
    TLorentzVector vec;
    vec.SetPtEtaPhiM(pt, eta, phi, mass);
    return vec;
}

// Returns true if a muon or electron passing the veto selection
// is within deltaR < dR of the given eta and phi coordinates
bool RazorAnalyzer::matchesVetoLepton(float eta, float phi, float dR) {
    for ( int i = 0; i < nMuons; i++ ) {
        if ( muonPt[i] < 5 || fabs(muonEta[i]) > 2.4 ) continue;
        if ( isVetoMuon(i) && deltaR(eta, phi, muonEta[i], muonPhi[i]) < dR ) {
            return true;
        }
    }
    for ( int i = 0; i < nElectrons; i++ ) {
        if ( elePt[i] < 5 || fabs(eleEta[i]) > 2.5 ) continue;
        if ( isVetoElectron(i) && deltaR(eta, phi, eleEta[i], elePhi[i]) < dR ) {
            return true;
        }
    }
    return false;
}

double RazorAnalyzer::GetAlphaT(vector<TLorentzVector> jets)
{
    int nJets = jets.size();
    vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations
    vector<TLorentzVector> possibleHem2s;
    double alphaT = 0;

    if(nJets < 2) return alphaT;

    int nComb = pow(2, nJets); // # possible combinations

    // steal from the getHemispheres method

    //step 1: store all possible partitions of the input jets
    int j_count;
    for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
        TLorentzVector j_temp1, j_temp2;
        int itemp = i;
        j_count = nComb/2;
        int count = 0;

        while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2
            if(itemp/j_count == 1){
                j_temp1 += jets[count];
            } else {
                j_temp2 += jets[count];
            }
            itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count
            j_count /= 2;
            count++;
        }
        possibleHem1s.push_back(j_temp1);
        possibleHem2s.push_back(j_temp2);
    }

    //step 2: Select combination that mininize |ET1 - ET2|
    double eMin = -1;
    TLorentzVector myHem1;
    TLorentzVector myHem2;

    for(size_t i=0; i < possibleHem1s.size(); i++)
    {
        double eTemp = fabs(possibleHem1s[i].Et() - possibleHem2s[i].Et());
        if (eMin < 0 || eTemp < eMin)
        {
            eMin = eTemp;
            myHem1 = possibleHem1s[i];
            myHem2 = possibleHem2s[i];
        }
    }

    float MhtX = 0., MhtY = 0.;
    float HT = 0.;
    for (auto& obj : jets) { HT += obj.Pt(); MhtX += obj.Px(); MhtY += obj.Py(); }

      TLorentzVector MyMHT;
      MyMHT.SetPxPyPzE(-MhtX, -MhtY, 0, sqrt(pow(MhtX,2) + pow(MhtY,2)));

    float MHT = MyMHT.Pt();

    // Calculate alphaT
    alphaT = 0.5 * (1-eMin/HT)/sqrt(1-pow(MHT/HT,2));

    return alphaT;
}

double RazorAnalyzer::GetDPhiMin(vector<TLorentzVector> jets)
    // This variable is used in the alphaT analysis
{
    double dPhiMin = -1.;
    float HT = 0.;
    float MhtX = 0.;
    float MhtY = 0.;
    // Search for min dPhi between recomputed missing HT and test jets
    for (auto& obj : jets) { HT += obj.Pt(); MhtX += obj.Px(); MhtY += obj.Py(); }
    TLorentzVector MyMHT;
    MyMHT.SetPxPyPzE(-MhtX, -MhtY, 0, sqrt(pow(MhtX,2) + pow(MhtY,2)));

    for (auto& obj : jets)
    {
    // Recompute MHT by ignoring a test jet
        float recomputedMHTX = MhtX - obj.Px();
        float recomputedMHTY = MhtY - obj.Py();
        TLorentzVector recomputedMHT;
        recomputedMHT.SetPxPyPzE(-recomputedMHTX, -recomputedMHTY, 0, sqrt(pow(recomputedMHTX,2) + pow(recomputedMHTY,2)));
        double phiTemp = fabs(recomputedMHT.Phi() - obj.Phi());
        if (dPhiMin < 0 || phiTemp < dPhiMin)   dPhiMin = phiTemp;
    }

    return dPhiMin;
}

vector<TLorentzVector> RazorAnalyzer::getHemispheres(vector<TLorentzVector> jets){
    int nJets = jets.size();
    vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations
    vector<TLorentzVector> possibleHem2s;

    if(nJets < 2){ //return empty hemispheres if there are fewer than 2 jets provided
        TLorentzVector emptyHem1, emptyHem2;
        vector<TLorentzVector> emptyHemsOut;
        emptyHemsOut.push_back(emptyHem1);
        emptyHemsOut.push_back(emptyHem2);
        return emptyHemsOut;
    }

    //stolen from https://github.com/pierinim/BSMatLHC/blob/master/BSMApp/src/CMS/CMSHemisphere.cc
    int nComb = pow(2, nJets);

    //step 1: store all possible partitions of the input jets
    int j_count;
    for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
        TLorentzVector j_temp1, j_temp2;
        int itemp = i;
        j_count = nComb/2;
        int count = 0;
        while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2
            if(itemp/j_count == 1){
                j_temp1 += jets[count];
            } else {
                j_temp2 += jets[count];
            }
            itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count
            j_count /= 2;
            count++;
        }
        possibleHem1s.push_back(j_temp1);
        possibleHem2s.push_back(j_temp2);
    }

    //step 2: choose the partition that minimizes m1^2 + m2^2
    double mMin = -1;
    TLorentzVector myHem1;
    TLorentzVector myHem2;
    for(size_t i=0; i < possibleHem1s.size(); i++){
        double mTemp = possibleHem1s[i].M2() + possibleHem2s[i].M2();
        if(mMin < 0 || mTemp < mMin){
            mMin = mTemp;
            myHem1 = possibleHem1s[i];
            myHem2 = possibleHem2s[i];
        }
    }

    //return the hemispheres in decreasing order of pt
    vector<TLorentzVector> hemsOut;
    if(myHem1.Pt() > myHem2.Pt()){
        hemsOut.push_back(myHem1);
        hemsOut.push_back(myHem2);
    } else {
        hemsOut.push_back(myHem2);
        hemsOut.push_back(myHem1);
    }

    return hemsOut;
}

std::vector< std::vector<int> > RazorAnalyzer::getHemispheresV2( std::vector<TLorentzVector> jets )
{
  //returns vector with original indices to jets
  int nJets = jets.size();
  vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations
  std::vector< std::vector<int> > index1;
  vector<TLorentzVector> possibleHem2s;
  std::vector< std::vector<int> > index2;

  if(nJets < 2){ //return empty hemispheres if there are fewer than 2 jets provided
    std::vector<int> emptyIndex1, emptyIndex2;
    std::vector< std::vector<int> > void_return;
    void_return.push_back( emptyIndex1 );
    void_return.push_back( emptyIndex2 );
    return void_return;
  }

  //stolen from https://github.com/pierinim/BSMatLHC/blob/master/BSMApp/src/CMS/CMSHemisphere.cc
  int nComb = pow(2, nJets);
  //std::cout << "njets: " << nJets << " ncomb: " << nComb << std::endl;
  //step 1: store all possible partitions of the input jets
  int j_count;
  for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
    //std::cout << "=iter: " << i << std::endl;
    TLorentzVector j_temp1, j_temp2;
    std::vector<int> tmp_index1, tmp_index2;
    int itemp = i;
    j_count = nComb/2;
    int count = 0;
    while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2
      //std::cout << "j_count: " << j_count << " itemp: " << itemp << " count: " << count << std::endl;
      if(itemp/j_count == 1){
	j_temp1 += jets[count];
	tmp_index1.push_back( count );
      } else {
	j_temp2 += jets[count];
	tmp_index2.push_back( count );
      }
      itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count
      j_count /= 2;
      count++;
    }
    possibleHem1s.push_back(j_temp1);
    index1.push_back( tmp_index1 );
    possibleHem2s.push_back(j_temp2);
    index2.push_back( tmp_index2 );
  }

  //step 2: choose the partition that minimizes m1^2 + m2^2
  double mMin = -1;
  TLorentzVector myHem1;
  TLorentzVector myHem2;
  int partition_index = -1;
  for(size_t i=0; i < possibleHem1s.size(); i++){
    double mTemp = possibleHem1s[i].M2() + possibleHem2s[i].M2();
    if(mMin < 0 || mTemp < mMin){
      mMin = mTemp;
      myHem1 = possibleHem1s[i];
      myHem2 = possibleHem2s[i];
      partition_index = i;
    }
  }

  //return the hemispheres in decreasing order of pt
  vector<TLorentzVector> hemsOut;
  std::vector< std::vector<int> > index_out;
  if(myHem1.Pt() > myHem2.Pt()){
    hemsOut.push_back(myHem1);
    hemsOut.push_back(myHem2);
    index_out.push_back( index1[partition_index] );
    index_out.push_back( index2[partition_index] );
  } else {
    hemsOut.push_back(myHem2);
    hemsOut.push_back(myHem1);
    index_out.push_back( index2[partition_index] );
    index_out.push_back( index1[partition_index] );
  }

  return index_out;
};


double RazorAnalyzer::computeMR(TLorentzVector hem1, TLorentzVector hem2){
    return sqrt(pow(hem1.P() + hem2.P(), 2) - pow(hem1.Pz() + hem2.Pz(), 2));
}

double RazorAnalyzer::computeRsq(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector pfMet){
    double mR = computeMR(hem1, hem2);
    double term1 = pfMet.Pt()/2*(hem1.Pt() + hem2.Pt());
    double term2 = pfMet.Px()/2*(hem1.Px() + hem2.Px()) + pfMet.Py()/2*(hem1.Py() + hem2.Py()); //dot product of MET with (p1T + p2T)
    double mTR = sqrt(term1 - term2);
    return (mTR / mR) * (mTR / mR);
}


double RazorAnalyzer::GetMT( TLorentzVector visible, TVector3 met )
{
  TVector3 vis( visible.Px(), visible.Py(), visible.Pz() );
  //return sqrt( visible.M2() + 2.0*( vis.Pt()*met.Pt() - vis.Dot( met ) ) );
  return sqrt( 2.0*( vis.Pt()*met.Pt() - vis.Dot( met ) ) );
};

double RazorAnalyzer::GetMT( TLorentzVector visible, TLorentzVector met )
{
  TVector3 _met( met.Px(), met.Py(), met.Pz() );
  return GetMT( visible, _met );
};


double RazorAnalyzer::GetMTEnergy( TLorentzVector visible, TVector3 met )
{
  TVector3 vis( visible.Px(), visible.Py(), visible.Pz() );
  //return sqrt( visible.M2() + 2.0*( visible.E()*met.Pt() - vis.Dot( met ) ) );
  return sqrt( 2.0*( visible.E()*met.Pt() - vis.Dot( met ) ) );
};

double RazorAnalyzer::GetMTEnergy( TLorentzVector visible, TLorentzVector met )
{
  TVector3 _met( met.Px(), met.Py(), met.Pz() );
  return GetMTEnergy( visible, _met );
};


double RazorAnalyzer::GetDphi( TLorentzVector visible, TVector3 met )
{
  TVector3 vis( visible.Px(), visible.Py(), visible.Pz() );
  return vis.DeltaPhi( met );
};

double RazorAnalyzer::GetDphi( TLorentzVector visible, TLorentzVector met )
{
  TVector3 _met( met.Px(), met.Py(), met.Pz() );
  return GetDphi( visible, _met );
};

//auxiliary functions for RazorInclusive and MatchedRazorInclusive analyses
bool RazorAnalyzer::passesHadronicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    if(MR < 0 || Rsq < 0) passes = false;
    //temporarily disable these
    //if(MR < 400 || Rsq < 0.25) passes = false;
    //if(MR < 450 && Rsq < 0.3) passes = false;
    return passes;
}

bool RazorAnalyzer::passesLeptonicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    if(MR < 0 || Rsq < 0) passes = false;
    //temporarily disable these
    // if(MR < 300 || Rsq < 0.15) passes = false;
    // if(MR < 350 && Rsq < 0.2) passes = false;
    return passes;
}

//Checks if ToSubtract matches any particles in Collection, and subtracts the momentum of ToSubtract from the closest one
int RazorAnalyzer::SubtractParticleFromCollection(TLorentzVector ToSubtract, vector<TLorentzVector>& Collection, float deltaRMatch){
    //if Collection has any elements within R<deltaRMatch of the vector ToSubtract, find the closest such element
    //otherwise, return -1
    double closestDR = -1;
    int closestDRIndex = -1;
    for(uint i = 0; i < Collection.size(); i++){
        double thisDR = Collection[i].DeltaR(ToSubtract);
        if(closestDR < 0 || thisDR < closestDR){
            closestDR = thisDR;
            closestDRIndex = i;
        }
    }
    if(closestDR < 0 || closestDR > deltaRMatch){ //if we didn't look at any objects or we didn't get one within 0.4, return -1
        return -1;
    }

    //subtract the momentum (magnitude) of ToSubtract from that of the closest vector in Collection

    //if we're subtracting away everything, just set the vector to 0 and return the index of the changed vector
    if(ToSubtract.P() >= Collection[closestDRIndex].P()){
        Collection[closestDRIndex].SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0);
        return closestDRIndex;
    }
    //otherwise, scale the 4-momentum by the appropriate factor
    double scalePFactor = (Collection[closestDRIndex].P() - ToSubtract.P())/Collection[closestDRIndex].P(); // = new P / old P
    Collection[closestDRIndex].SetPxPyPzE(Collection[closestDRIndex].Px()*scalePFactor, Collection[closestDRIndex].Py()*scalePFactor, Collection[closestDRIndex].Pz()*scalePFactor, Collection[closestDRIndex].E()*scalePFactor);
    return closestDRIndex;
}

double RazorAnalyzer::calcMT2(float testMass, bool massive, vector<TLorentzVector> jets, TLorentzVector MET, int hemi_seed, int hemi_association)
{
  //computes MT2 using a test mass of testMass, with hemispheres made massless if massive is set to false
  //hemispheres are clustered by finding the grouping of input jets that minimizes the Lund distance

  if(jets.size() < 2) return -9999; //need at least two jets for the calculation
  vector<float> px, py, pz, E;
  for(uint i = 0; i < jets.size(); i++){
    //push 4vector components onto individual lists
    px.push_back(jets[i].Px());
    py.push_back(jets[i].Py());
    pz.push_back(jets[i].Pz());
    E.push_back(jets[i].E());
  }

  //form the hemispheres using the provided Hemisphere class
  Hemisphere* hemis = new Hemisphere(px, py, pz, E, hemi_seed, hemi_association);
  vector<int> grouping = hemis->getGrouping();
  TLorentzVector pseudojet1(0.,0.,0.,0.);
  TLorentzVector pseudojet2(0.,0.,0.,0.);

  //make the hemisphere vectors
  for(uint i=0; i<px.size(); ++i){
    if(grouping[i]==1){
      pseudojet1.SetPx(pseudojet1.Px() + px[i]);
      pseudojet1.SetPy(pseudojet1.Py() + py[i]);
      pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
      pseudojet1.SetE( pseudojet1.E()  + E[i]);
    }else if(grouping[i] == 2){
      pseudojet2.SetPx(pseudojet2.Px() + px[i]);
      pseudojet2.SetPy(pseudojet2.Py() + py[i]);
      pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
      pseudojet2.SetE( pseudojet2.E()  + E[i]);
    }
  }
  delete hemis;

  //now compute MT2 using the Davismt2 class

  //these arrays contain (mass, px, py) for the pseudojets and the MET
  double pa[3];
  double pb[3];
  double pmiss[3];

  pmiss[0] = 0;
  pmiss[1] = static_cast<double> (MET.Px());
  pmiss[2] = static_cast<double> (MET.Py());

  pa[0] = static_cast<double> (massive ? pseudojet1.M() : 0);
  pa[1] = static_cast<double> (pseudojet1.Px());
  pa[2] = static_cast<double> (pseudojet1.Py());

  pb[0] = static_cast<double> (massive ? pseudojet2.M() : 0);
  pb[1] = static_cast<double> (pseudojet2.Px());
  pb[2] = static_cast<double> (pseudojet2.Py());

  Davismt2 *mt2 = new Davismt2();
  mt2->set_momenta(pa, pb, pmiss);
  mt2->set_mn(testMass);
  Float_t MT2=mt2->get_mt2();
  delete mt2;
  return MT2;
};

////conversion between DetId <-> ieta/ix/iphi/iy

int RazorAnalyzer::detID_from_iEtaiPhi(int iEta_or_iX=1, int iPhi_or_iY=1, bool isEB = true, bool isEEMinus = false)
{
	uint32_t detID = 0;
	int Ecal = 3;
	int EcalBarrel=1;
	int EcalEndcap=2;
	int iz = isEEMinus?-1:1;

	if(isEB)
	{
		detID = ((Ecal&0xF)<<28)|((EcalBarrel&0x7)<<25);
		detID |= ((iEta_or_iX>0)?(0x10000|(iEta_or_iX<<9)):((-iEta_or_iX)<<9))|(iPhi_or_iY&0x1FF);
	}
	else
	{
		detID = ((Ecal&0xF)<<28)|((EcalEndcap&0x7)<<25);

		detID |=(iPhi_or_iY&0x7f)|((iEta_or_iX&0x7f)<<7)|((iz>0)?(0x4000):(0));
	}

	return int(detID);
};


int RazorAnalyzer::iEta_or_iX_from_detID(int detID=1, bool isEB = true)
{
	int iEta_or_iX = 0;
	uint32_t id_ = uint32_t(detID);
	if(isEB)
	{
		int zside = (id_&0x10000)?(1):(-1);
		int ietaAbs = (id_>>9)&0x7F;
		iEta_or_iX =  zside*ietaAbs;
	}
	else
	{
		iEta_or_iX = (id_>>7)&0x7F;
	}
	return iEta_or_iX;

};

int RazorAnalyzer::iPhi_or_iY_from_detID(int detID=1, bool isEB = true)
{
	int iPhi_or_iY = 0;
	uint32_t id_ = uint32_t(detID);
	if(isEB)
	{
		iPhi_or_iY =  id_&0x1FF;
	}
	else
	{
		iPhi_or_iY = id_&0x7F;
	}
	return iPhi_or_iY;
};
