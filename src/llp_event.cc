#ifndef llp_event_STANDALONE
#include "LLPAnalysis/llpAnalyzer/interface/llp_event.h"

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
	fChain->SetBranchAddress("ele_passCutBasedIDVeto",    ele_passCutBasedIDVeto,      &b_ele_passCutBasedIDVetog); 
	fChain->SetBranchAddress("ele_passCutBasedIDLoose",   ele_passCutBasedIDLoose,     &b_ele_passCutBasedIDLooseg); 
	fChain->SetBranchAddress("ele_passCutBasedIDMedium",  ele_passCutBasedIDMedium,    &b_ele_passCutBasedIDMediumg); 
	fChain->SetBranchAddress("ele_passCutBasedIDTight",   ele_passCutBasedIDTight,     &b_ele_passCutBasedIDTightg); 
	fChain->SetBranchAddress("ele_passMVAIsoIDWP80",      ele_passMVAIsoIDWP80,        &b_ele_passMVAIsoIDWP80g); 
	fChain->SetBranchAddress("ele_passMVAIsoIDWP90",      ele_passMVAIsoIDWP90,        &b_ele_passMVAIsoIDWP90g); 
	fChain->SetBranchAddress("ele_passMVAIsoIDWPHZZ",     ele_passMVAIsoIDWPHZZ,       &b_ele_passMVAIsoIDWPHZZg); 
	fChain->SetBranchAddress("ele_passMVAIsoIDWPLoose",   ele_passMVAIsoIDWPLoose,     &b_ele_passMVAIsoIDWPLooseg); 
	fChain->SetBranchAddress("ele_passMVANoIsoIDWP80",    ele_passMVANoIsoIDWP80,      &b_ele_passMVANoIsoIDWP80g); 
	fChain->SetBranchAddress("ele_passMVANoIsoIDWP90",    ele_passMVANoIsoIDWP90,      &b_ele_passMVANoIsoIDWP90g); 
	fChain->SetBranchAddress("ele_passMVANoIsoIDWPLoose", ele_passMVANoIsoIDWPLoose,   &b_ele_passMVANoIsoIDWPLooseg); 
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
	fChain->SetBranchAddress("pho_passCutBasedIDLoose",  pho_passCutBasedIDLoose,  &b_pho_passCutBasedIDLoose);
	fChain->SetBranchAddress("pho_passCutBasedIDMedium", pho_passCutBasedIDMedium, &b_pho_passCutBasedIDMedium);
	fChain->SetBranchAddress("pho_passCutBasedIDTight",  pho_passCutBasedIDTight,  &b_pho_passCutBasedIDTight);
	fChain->SetBranchAddress("pho_passMVAIDWP80",  pho_passMVAIDWP80,  &b_pho_passMVAIDWP80);
	fChain->SetBranchAddress("pho_passMVAIDWP90",  pho_passMVAIDWP90,  &b_pho_passMVAIDWP90);
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
	fChain->SetBranchAddress("nCsc", &nCsc, &b_nCsc);
	fChain->SetBranchAddress("cscPhi", cscPhi, &b_cscPhi);
	fChain->SetBranchAddress("cscEta", cscEta, &b_cscEta);
	fChain->SetBranchAddress("cscX", cscX, &b_cscX);
	fChain->SetBranchAddress("cscY", cscY, &b_cscY);
	fChain->SetBranchAddress("cscZ", cscZ, &b_cscZ);
	fChain->SetBranchAddress("cscDirectionX", cscDirectionX, &b_cscDirectionX);
	fChain->SetBranchAddress("cscDirectionY", cscDirectionY, &b_cscDirectionY);
	fChain->SetBranchAddress("cscDirectionZ", cscDirectionZ, &b_cscDirectionZ);
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
	fChain->SetBranchAddress("rpcTError", rpcTError, &b_rpcTError);
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

	fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
	fChain->SetBranchAddress("jetE", &jetE, &b_jetE);
	fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
	fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);

	fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
	fChain->SetBranchAddress("jetCSV", &jetCSV, &b_jetCSV);
	fChain->SetBranchAddress("jetCISV", &jetCISV, &b_jetCISV);
	fChain->SetBranchAddress("jetProbb", &jetProbb, &b_jetProbb);
	fChain->SetBranchAddress("jetProbc", &jetProbc, &b_jetProbc);
	fChain->SetBranchAddress("jetProbudsg", &jetProbudsg, &b_jetProbudsg);
	fChain->SetBranchAddress("jetProbbb", &jetProbbb, &b_jetProbbb);
	fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
	fChain->SetBranchAddress("jetJetArea", &jetJetArea, &b_jetJetArea);
	fChain->SetBranchAddress("jetPileupE", &jetPileupE, &b_jetPileupE);
	fChain->SetBranchAddress("jetPileupId", &jetPileupId, &b_jetPileupId);
	fChain->SetBranchAddress("jetPileupIdFlag", &jetPileupIdFlag, &b_jetPileupIdFlag);
	fChain->SetBranchAddress("jetPassIDLoose", &jetPassIDLoose, &b_jetPassIDLoose);
	fChain->SetBranchAddress("jetPassIDTight", &jetPassIDTight, &b_jetPassIDTight);
	fChain->SetBranchAddress("jetPassMuFrac", &jetPassMuFrac, &b_jetPassMuFrac);
	fChain->SetBranchAddress("jetPassEleFrac", &jetPassEleFrac, &b_jetPassEleFrac);
	fChain->SetBranchAddress("jetPartonFlavor", &jetPartonFlavor, &b_jetPartonFlavor);
	fChain->SetBranchAddress("jetHadronFlavor", &jetHadronFlavor, &b_jetHadronFlavor);
	fChain->SetBranchAddress("jetElectronEnergyFraction", &jetElectronEnergyFraction, &b_jetElectronEnergyFraction);
	fChain->SetBranchAddress("jetPhotonEnergyFraction", &jetPhotonEnergyFraction, &b_jetPhotonEnergyFraction);
	fChain->SetBranchAddress("jetChargedHadronEnergyFraction", &jetChargedHadronEnergyFraction, &b_jetChargedHadronEnergyFraction);
	fChain->SetBranchAddress("jetNeutralHadronEnergyFraction", &jetNeutralHadronEnergyFraction, &b_jetNeutralHadronEnergyFraction);
	fChain->SetBranchAddress("jetMuonEnergyFraction", &jetMuonEnergyFraction, &b_jetMuonEnergyFraction);
	fChain->SetBranchAddress("jetHOEnergyFraction", &jetHOEnergyFraction, &b_jetHOEnergyFraction);
	fChain->SetBranchAddress("jetHFHadronEnergyFraction", &jetHFHadronEnergyFraction, &b_jetHFHadronEnergyFraction);
	fChain->SetBranchAddress("jetHFEMEnergyFraction", &jetHFEMEnergyFraction, &b_jetHFEMEnergyFraction);
	fChain->SetBranchAddress("jetChargedHadronMultiplicity", &jetChargedHadronMultiplicity, &b_jetChargedHadronMultiplicity);
	fChain->SetBranchAddress("jetNeutralHadronMultiplicity", &jetNeutralHadronMultiplicity, &b_jetNeutralHadronMultiplicity);
	fChain->SetBranchAddress("jetPhotonMultiplicity", &jetPhotonMultiplicity, &b_jetPhotonMultiplicity);
	fChain->SetBranchAddress("jetElectronMultiplicity", &jetElectronMultiplicity, &b_jetElectronMultiplicity);
	fChain->SetBranchAddress("jetMuonMultiplicity", &jetMuonMultiplicity, &b_jetMuonMultiplicity);
	fChain->SetBranchAddress("jetAllMuonPt", &jetAllMuonPt, &b_jetAllMuonPt);
	fChain->SetBranchAddress("jetAllMuonEta", &jetAllMuonEta, &b_jetAllMuonEta);
	fChain->SetBranchAddress("jetAllMuonPhi", &jetAllMuonPhi, &b_jetAllMuonPhi);
	fChain->SetBranchAddress("jetAllMuonM", &jetAllMuonM, &b_jetAllMuonM);
	fChain->SetBranchAddress("jetPtWeightedDZ", &jetPtWeightedDZ, &b_jetPtWeightedDZ);
	fChain->SetBranchAddress("jetNRechits", &jetNRechits, &b_jetNRechits);
	fChain->SetBranchAddress("jetRechitE", &jetRechitE, &b_jetRechitE);
	fChain->SetBranchAddress("jetRechitT", &jetRechitT, &b_jetRechitT);
	fChain->SetBranchAddress("jetRechitT_rms", &jetRechitT_rms, &b_jetRechitT_rms);

	fChain->SetBranchAddress("jetRechitE_Error", &jetRechitE_Error, &b_jetRechitE_Error);
	fChain->SetBranchAddress("jetRechitT_Error", &jetRechitT_Error, &b_jetRechitT_Error);
	fChain->SetBranchAddress("jetAlphaMax", &jetAlphaMax, &b_jetAlphaMax);
	fChain->SetBranchAddress("jetBetaMax", &jetBetaMax, &b_jetBetaMax);
	fChain->SetBranchAddress("jetGammaMax_ET", &jetGammaMax_ET, &b_jetGammaMax_ET);
	fChain->SetBranchAddress("jetGammaMax_EM", &jetGammaMax_EM, &b_jetGammaMax_EM);
	fChain->SetBranchAddress("jetGammaMax_Hadronic", &jetGammaMax_Hadronic, &b_jetGammaMax_Hadronic);
	fChain->SetBranchAddress("jetGammaMax", &jetGammaMax, &b_jetGammaMax);
	fChain->SetBranchAddress("jetPtAllTracks", &jetPtAllTracks, &b_jetPtAllTracks);
	fChain->SetBranchAddress("jetPtAllPVTracks", &jetPtAllPVTracks, &b_jetPtAllPVTracks);
	fChain->SetBranchAddress("jetMedianTheta2D", &jetMedianTheta2D, &b_jetMedianTheta2D);
	fChain->SetBranchAddress("jetMedianIP", &jetMedianIP, &b_jetMedianIP);
	fChain->SetBranchAddress("jetMinDeltaRAllTracks", &jetMinDeltaRAllTracks, &b_jetMinDeltaRAllTracks);
	fChain->SetBranchAddress("jetMinDeltaRPVTracks", &jetMinDeltaRPVTracks, &b_jetMinDeltaRPVTracks);
	fChain->SetBranchAddress("jet_sig_et1", &jet_sig_et1, &b_jet_sig_et1);
	fChain->SetBranchAddress("jet_sig_et2", &jet_sig_et2, &b_jet_sig_et2);
	fChain->SetBranchAddress("jet_energy_frac", &jet_energy_frac, &b_jet_energy_frac);

	fChain->SetBranchAddress("jetAlphaMax_wp", &jetAlphaMax_wp, &b_jetAlphaMax_wp);
	fChain->SetBranchAddress("jetBetaMax_wp", &jetBetaMax_wp, &b_jetBetaMax_wp);
	fChain->SetBranchAddress("jetGammaMax_ET_wp", &jetGammaMax_ET_wp, &b_jetGammaMax_ET_wp);
	fChain->SetBranchAddress("jetGammaMax_EM_wp", &jetGammaMax_EM_wp, &b_jetGammaMax_EM_wp);
	fChain->SetBranchAddress("jetGammaMax_Hadronic_wp", &jetGammaMax_Hadronic_wp, &b_jetGammaMax_Hadronic_wp);
	fChain->SetBranchAddress("jetGammaMax_wp", &jetGammaMax_wp, &b_jetGammaMax_wp);
	fChain->SetBranchAddress("jetPtAllTracks_wp", &jetPtAllTracks_wp, &b_jetPtAllTracks_wp);
	fChain->SetBranchAddress("jetPtAllPVTracks_wp", &jetPtAllPVTracks_wp, &b_jetPtAllPVTracks_wp);
	fChain->SetBranchAddress("jetMedianTheta2D_wp", &jetMedianTheta2D_wp, &b_jetMedianTheta2D_wp);
	fChain->SetBranchAddress("jetMedianIP_wp", &jetMedianIP_wp, &b_jetMedianIP_wp);
	fChain->SetBranchAddress("jetMinDeltaRAllTracks_wp", &jetMinDeltaRAllTracks_wp, &b_jetMinDeltaRAllTracks_wp);
	fChain->SetBranchAddress("jetMinDeltaRPVTracks_wp", &jetMinDeltaRPVTracks_wp, &b_jetMinDeltaRPVTracks_wp);

	fChain->SetBranchAddress("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction, &b_jetChargedEMEnergyFraction);
	fChain->SetBranchAddress("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction, &b_jetNeutralEMEnergyFraction);

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
	fChain->SetBranchAddress("Flag2_globalSuperTightHalo2016Filter",              &Flag2_globalSuperTightHalo2016Filter,            &b_Flag2_globalSuperTightHalo2016Filter);
	fChain->SetBranchAddress("Flag2_globalTightHalo2016Filter",                   &Flag2_globalTightHalo2016Filter,                 &b_Flag2_globalTightHalo2016Filter);
	fChain->SetBranchAddress("Flag2_goodVertices",                                &Flag2_goodVertices,                              &b_Flag2_goodVertices);
	fChain->SetBranchAddress("Flag2_BadChargedCandidateFilter",                   &Flag2_BadChargedCandidateFilter,                 &b_Flag2_BadChargedCandidateFilter);
	fChain->SetBranchAddress("Flag2_BadPFMuonFilter",                             &Flag2_BadPFMuonFilter,                           &b_Flag2_BadPFMuonFilter);
	fChain->SetBranchAddress("Flag2_EcalDeadCellTriggerPrimitiveFilter",          &Flag2_EcalDeadCellTriggerPrimitiveFilter,        &b_Flag2_EcalDeadCellTriggerPrimitiveFilter);
	fChain->SetBranchAddress("Flag2_HBHENoiseFilter",                             &Flag2_HBHENoiseFilter,                           &b_Flag2_HBHENoiseFilter);
	fChain->SetBranchAddress("Flag2_HBHEIsoNoiseFilter",                          &Flag2_HBHEIsoNoiseFilter,                        &b_Flag2_HBHEIsoNoiseFilter);
	fChain->SetBranchAddress("Flag2_ecalBadCalibFilter",                          &Flag2_ecalBadCalibFilter,                        &b_Flag2_ecalBadCalibFilter);
	fChain->SetBranchAddress("Flag2_eeBadScFilter",                               &Flag2_eeBadScFilter,                             &b_Flag2_eeBadScFilter);
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

	fChain->SetBranchAddress("gLLP_travel_time", gLLP_travel_time, &b_gLLP_travel_time);
	fChain->SetBranchAddress("gLLP_e", gLLP_e, &b_gLLP_e);
	fChain->SetBranchAddress("gLLP_pt", gLLP_pt, &b_gLLP_pt);
	fChain->SetBranchAddress("gLLP_eta", gLLP_eta, &b_gLLP_eta);
	fChain->SetBranchAddress("gLLP_phi", gLLP_phi, &b_gLLP_phi);
	fChain->SetBranchAddress("gLLP_beta", gLLP_beta, &b_gLLP_beta);

	fChain->SetBranchAddress("gLLP_decay_vertex_x", gLLP_decay_vertex_x, &b_gLLP_decay_vertex_x);
	fChain->SetBranchAddress("gLLP_decay_vertex_y", gLLP_decay_vertex_y, &b_gLLP_decay_vertex_y);
	fChain->SetBranchAddress("gLLP_decay_vertex_z", gLLP_decay_vertex_z, &b_gLLP_decay_vertex_z);

	fChain->SetBranchAddress("gLLP_prod_vertex_x", gLLP_prod_vertex_x, &b_gLLP_prod_vertex_x);
	fChain->SetBranchAddress("gLLP_prod_vertex_y", gLLP_prod_vertex_y, &b_gLLP_prod_vertex_y);
	fChain->SetBranchAddress("gLLP_prod_vertex_z", gLLP_prod_vertex_z, &b_gLLP_prod_vertex_z);
	/*
	   fChain->SetBranchAddress("gen_time", gen_time, &b_gen_time);
	   fChain->SetBranchAddress("photon_travel_time", photon_travel_time, &b_photon_travel_time);
	   fChain->SetBranchAddress("gLLP_daughter_travel_time", gLLP_daughter_travel_time, &b_gLLP_daughter_travel_time);
	   fChain->SetBranchAddress("gLLP_daughter_e", gLLP_daughter_e, &b_gLLP_daughter_e);
	   fChain->SetBranchAddress("gLLP_daughter_pt", gLLP_daughter_pt, &b_gLLP_daughter_pt);
	   fChain->SetBranchAddress("gLLP_daughter_eta", gLLP_daughter_eta, &b_gLLP_daughter_eta);
	   fChain->SetBranchAddress("gLLP_daughter_phi", gLLP_daughter_phi, &b_gLLP_daughter_phi);
	   fChain->SetBranchAddress("gLLP_daughter_eta_ecalcorr", gLLP_daughter_eta_ecalcorr, &b_gLLP_daughter_eta_ecalcorr);
	   fChain->SetBranchAddress("gLLP_daughter_phi_ecalcorr", gLLP_daughter_phi_ecalcorr, &b_gLLP_daughter_phi_ecalcorr);
	   fChain->SetBranchAddress("gLLP_min_delta_r_match_jet", gLLP_min_delta_r_match_jet, &b_gLLP_min_delta_r_match_jet);
	   fChain->SetBranchAddress("gLLP_daughter_match_jet_index", gLLP_daughter_match_jet_index, &b_gLLP_daughter_match_jet_index);
	   */
	//daughters
	fChain->SetBranchAddress("gLLP_daughter_EB", gLLP_daughter_EB, &b_gLLP_daughter_EB );
	fChain->SetBranchAddress("gLLP_daughter_ETL", gLLP_daughter_ETL, &b_gLLP_daughter_ETL);

	fChain->SetBranchAddress("gLLP_daughter_photon_travel_time_EB", gLLP_daughter_photon_travel_time_EB, &b_gLLP_daughter_photon_travel_time_EB);
	fChain->SetBranchAddress("gLLP_daughter_photon_travel_time_ETL", gLLP_daughter_photon_travel_time_ETL, &b_gLLP_daughter_photon_travel_time_ETL);

	fChain->SetBranchAddress("gLLP_daughter_travel_time_EB", gLLP_daughter_travel_time_EB, &b_gLLP_daughter_travel_time_EB);
	fChain->SetBranchAddress("gLLP_daughter_travel_time_ETL", gLLP_daughter_travel_time_ETL, &b_gLLP_daughter_travel_time_ETL);

	fChain->SetBranchAddress("gen_time_daughter_EB", gen_time_daughter_EB, &b_gen_time_daughter_EB);
	fChain->SetBranchAddress("gen_time_daughter_ETL", gen_time_daughter_ETL, &b_gen_time_daughter_ETL);

	fChain->SetBranchAddress("gLLP_daughter_id", gLLP_daughter_id, &b_gLLP_daughter_id);
	fChain->SetBranchAddress("gLLP_daughter_pt", gLLP_daughter_pt, &b_gLLP_daughter_pt);
	fChain->SetBranchAddress("gLLP_daughter_eta", gLLP_daughter_eta, &b_gLLP_daughter_eta);
	fChain->SetBranchAddress("gLLP_daughter_phi", gLLP_daughter_phi, &b_gLLP_daughter_phi);
	fChain->SetBranchAddress("gLLP_daughter_eta_ecalcorr", gLLP_daughter_eta_ecalcorr, &b_gLLP_daughter_eta_ecalcorr);
	fChain->SetBranchAddress("gLLP_daughter_phi_ecalcorr", gLLP_daughter_phi_ecalcorr, &b_gLLP_daughter_phi_ecalcorr);
	fChain->SetBranchAddress("gLLP_daughter_e", gLLP_daughter_e, &b_gLLP_daughter_e);
	fChain->SetBranchAddress("gLLP_daughter_mass", gLLP_daughter_mass, &b_gLLP_daughter_mass);

	fChain->SetBranchAddress("gLLP_daughter_match_jet_index", gLLP_daughter_match_jet_index, &b_gLLP_daughter_match_jet_index);
	fChain->SetBranchAddress("gLLP_daughter_min_delta_r_match_jet", gLLP_daughter_min_delta_r_match_jet, &b_gLLP_daughter_min_delta_r_match_jet);

	//grandaughters
	fChain->SetBranchAddress("gLLP_grandaughter_EB", gLLP_grandaughter_EB, &b_gLLP_grandaughter_EB );
	fChain->SetBranchAddress("gLLP_grandaughter_ETL", gLLP_grandaughter_ETL, &b_gLLP_grandaughter_ETL);

	fChain->SetBranchAddress("gLLP_grandaughter_photon_travel_time_EB", gLLP_grandaughter_photon_travel_time_EB, &b_gLLP_grandaughter_photon_travel_time_EB);
	fChain->SetBranchAddress("gLLP_grandaughter_photon_travel_time_ETL", gLLP_grandaughter_photon_travel_time_ETL, &b_gLLP_grandaughter_photon_travel_time_ETL);

	fChain->SetBranchAddress("gLLP_grandaughter_travel_time_EB", gLLP_grandaughter_travel_time_EB, &b_gLLP_grandaughter_travel_time_EB);
	fChain->SetBranchAddress("gLLP_grandaughter_travel_time_ETL", gLLP_grandaughter_travel_time_ETL, &b_gLLP_grandaughter_travel_time_ETL);

	fChain->SetBranchAddress("gen_time_grandaughter_EB", gen_time_grandaughter_EB, &b_gen_time_grandaughter_EB);
	fChain->SetBranchAddress("gen_time_grandaughter_ETL", gen_time_grandaughter_ETL, &b_gen_time_grandaughter_ETL);

	fChain->SetBranchAddress("gLLP_grandaughter_id", gLLP_grandaughter_id, &b_gLLP_grandaughter_id);
	fChain->SetBranchAddress("gLLP_grandaughter_pt", gLLP_grandaughter_pt, &b_gLLP_grandaughter_pt);
	fChain->SetBranchAddress("gLLP_grandaughter_eta", gLLP_grandaughter_eta, &b_gLLP_grandaughter_eta);
	fChain->SetBranchAddress("gLLP_grandaughter_phi", gLLP_grandaughter_phi, &b_gLLP_grandaughter_phi);
	fChain->SetBranchAddress("gLLP_grandaughter_eta_ecalcorr", gLLP_grandaughter_eta_ecalcorr, &b_gLLP_grandaughter_eta_ecalcorr);
	fChain->SetBranchAddress("gLLP_grandaughter_phi_ecalcorr", gLLP_grandaughter_phi_ecalcorr, &b_gLLP_grandaughter_phi_ecalcorr);
	fChain->SetBranchAddress("gLLP_grandaughter_e", gLLP_grandaughter_e, &b_gLLP_grandaughter_e);
	fChain->SetBranchAddress("gLLP_grandaughter_mass", gLLP_grandaughter_mass, &b_gLLP_grandaughter_mass);

	fChain->SetBranchAddress("gLLP_grandaughter_match_jet_index", gLLP_grandaughter_match_jet_index, &b_gLLP_grandaughter_match_jet_index);
	fChain->SetBranchAddress("gLLP_grandaughter_min_delta_r_match_jet", gLLP_grandaughter_min_delta_r_match_jet, &b_gLLP_grandaughter_min_delta_r_match_jet);

	fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
	fChain->SetBranchAddress("track_Pt", track_Pt, &b_track_Pt);
	fChain->SetBranchAddress("track_Eta", track_Eta, &b_track_Eta);
	fChain->SetBranchAddress("track_Phi", track_Phi, &b_track_Phi);
	fChain->SetBranchAddress("track_charge", track_charge, &b_track_charge);
	fChain->SetBranchAddress("track_bestVertexIndex", track_bestVertexIndex, &b_track_bestVertexIndex);
	fChain->SetBranchAddress("track_nMissingInnerHits", track_nMissingInnerHits, &b_track_nMissingInnerHits);
	fChain->SetBranchAddress("track_nMissingOuterHits", track_nMissingOuterHits, &b_track_nMissingOuterHits);
	fChain->SetBranchAddress("track_nPixelHits", track_nPixelHits, &b_track_nPixelHits);
	fChain->SetBranchAddress("track_nHits", track_nHits, &b_track_nHits); 
	fChain->SetBranchAddress("track_angle", track_angle, &b_track_angle);
	fChain->SetBranchAddress("track_dxyToBS", track_dxyToBS, &b_track_dxyToBS);
	fChain->SetBranchAddress("track_dxyErr", track_dxyErr, &b_track_dxyErr);
	fChain->SetBranchAddress("track_dzToPV", track_dzToPV, &b_track_dzToPV);
	fChain->SetBranchAddress("track_dzErr", track_dzErr, &b_track_dzErr);
	fChain->SetBranchAddress("track_chi2", track_chi2, &b_track_chi2);
	fChain->SetBranchAddress("track_ndof", track_ndof, &b_track_ndof);

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

#endif
 
