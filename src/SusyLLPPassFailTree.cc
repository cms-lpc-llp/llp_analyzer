#include "RazorHelper.h"
#include "SusyLLPPassFailTree.h"
#include "assert.h"
#include "TTree.h"

// Constructor
SusyLLPPassFailTree::SusyLLPPassFailTree()
{
	InitVariables();
};

SusyLLPPassFailTree::~SusyLLPPassFailTree()
{
	if (f_) f_->Close();
};

//Initialize
void SusyLLPPassFailTree::InitVariables()
{ 
	//evnt info
	runNum=0; lumiSec=0; evtNum=0; category=0;
	npv=0; npu=0; rho=-1; weight=-1;
	pileupWeight = 1; pileupWeightUp = 1; pileupWeightDown = 1;
	met=-1; metPhi=-1;
	HT=-1;
	jetMet_dPhi=-999.; jetMet_dPhiMin=999.; jetMet_dPhiMin4=999.;
	jetMet_dPhiStar=-999.; jetMet_dPhiStarMin=999.; 

	jetPho_dPhi=-999.; jetPho_dPhiMin=999.; jetPho_dPhiMin4=999.;
	jetPho_dPhiStar=-999.; jetPho_dPhiStarMin=999.; 

	jet2_dPhi=999.;

	// met filters
	Flag2_globalSuperTightHalo2016Filter = 0;
	Flag2_globalTightHalo2016Filter = 0;
	Flag2_goodVertices = 0;
	Flag2_BadChargedCandidateFilter = 0;
	Flag2_BadPFMuonFilter = 0;
	Flag2_EcalDeadCellTriggerPrimitiveFilter = 0;
	Flag2_HBHENoiseFilter = 0;
	Flag2_HBHEIsoNoiseFilter = 0;
	Flag2_ecalBadCalibFilter = 0;
	Flag2_eeBadScFilter = 0;


	//leptons
	nLeptons = 0;
	//MT = 0;
	for( int i = 0; i < N_MAX_LEPTONS; i++ )
	{
		lepE[i]      = -999.;
		lepPt[i]     = -999.;
		lepEta[i]    = -999.;
		lepPhi[i]    = -999.;
		lepPdgId[i]  = -999;
		lepDZ[i]     = -999.;
		// lepLoosePassId[i] = false;
		// lepMediumPassId[i] = false;
		// lepTightPassId[i] = false;
		lepPassId[i] = false;
	}

	//Z-candidate
	ZMass = -999.; ZPt = -999.; ZEta = -999.; ZPhi = -999.;
	MT = -999.;
	ZleptonIndex1 = -999; ZleptonIndex2 = -999;

	//HH-candidate
	dr_max = 999; dm_hh = 999; mh1 = -1; mh2 = -1; avg_mh = -1;

	//muons
	nMuons = 0;
	for( int i = 0; i < N_MAX_MUONS; i++ )
	{
		muonE[i] = -999;
		muonPt[i] = -999;
		muonEta[i] = -999;
		muonPhi[i] = -999;
		muonCharge[i] = 0;
		muonIsLoose[i] = 0;
		muonIsMedium[i] = 0;
		muonIsTight[i] = 0;
		muon_pileupIso[i] = -999;
		muon_chargedIso[i] = -999;
		muon_photonIso[i] = -999;
		muon_neutralHadIso[i] = -999;
	}

	//electrons
	nElectrons = 0;
	for( int i = 0; i < N_MAX_ELECTRONS; i++ )
	{
		eleE[i] = -999;
		elePt[i] = -999;
		eleEta[i] = -999;
		elePhi[i] = -999;
		eleCharge[i] = 0;
		ele_passCutBasedIDVeto[i] = 0;
		ele_passCutBasedIDLoose[i] = 0;
		ele_passCutBasedIDMedium[i] = 0;
		ele_passCutBasedIDTight[i] = 0;
	}

	//taus
	nTaus = 0; 
	for( int i = 0; i < N_MAX_TAUS; i++ )
	{
		tauE[i] = -999;
		tauPt[i] = -999;
		tauEta[i] = -999;
		tauPhi[i] = -999;
		tau_IsLoose[i] = 0;
		tau_IsMedium[i] = 0;
		tau_IsTight[i] = 0;
		tau_passEleVetoLoose[i] = 0;
		tau_passEleVetoMedium[i] = 0;
		tau_passEleVetoTight[i] = 0;
		tau_passMuVetoLoose[i] = 0;
		tau_passMuVetoMedium[i] = 0;
		tau_passMuVetoTight[i] = 0;
	}

	//photons
	nPhotons = 0;
	for( int i = 0; i < N_MAX_PHOTONS; i++ )
	{
		phoE[i] = -999;
		phoPt[i] = -999;
		phoEta[i] = -999;
		phoPhi[i] = -999;
		pho_passCutBasedIDLoose[i] = 0;
		pho_passCutBasedIDMedium[i] = 0;
		pho_passCutBasedIDTight[i] = 0;
		pho_passMVAIDWP80[i] = 0;
		pho_passMVAIDWP90[i] = 0;
	}

	//jets
	nJets = 0;
	//nBJets = 0;
	//nCaloJets = 0;
	for( int i = 0; i < N_MAX_JETS; i++ )
	{
		jetE[i] = 0.0;
		jetPt[i] = 0.0;
		jetPass[i] = 0;
		
		jetEta[i] = 0.0;
		jetPhi[i] = 0.0;
		jetCSV[i] = 0.0;
		jetCISV[i] = 0.0;
		jetMass[i] =  -99.0;
		jetJetArea[i] = -99.0;
		jetPileupE[i] = -99.0;
		jetPileupId[i] = -99.0;
		jetPileupIdFlag[i] = -1;
		jetPassIDLoose[i] = false;
		jetPassIDTight[i] = false;
		jetPassMuFrac[i] = false;
		jetPassEleFrac[i] = false;
		jetPartonFlavor[i] = 0;
		jetHadronFlavor[i] = 0;
		jetElectronEnergyFraction[i] = -99.0;
		jetPhotonEnergyFraction[i] = -99.0;
		jetChargedHadronEnergyFraction[i] = -99.0;
		jetNeutralHadronEnergyFraction[i] = -99.0;
		jetMuonEnergyFraction[i] = -99.0;
		jetHOEnergyFraction[i] = -99.0;
		jetHFHadronEnergyFraction[i] = -99.0;
		jetHFEMEnergyFraction[i] = -99.0;
		jetChargedHadronMultiplicity[i] = 0;
		jetNeutralHadronMultiplicity[i] = 0;
		jetPhotonMultiplicity[i] = 0;
		jetElectronMultiplicity[i] = 0;
		jetMuonMultiplicity[i] = 0;
		jetQGLikelihood[i] = -999;
		jetNSV[i] = 0;
		jetNSVCand[i] = 0;
		jetNVertexTracks[i] = 0;
		jetNSelectedTracks[i] = 0;
		jetDRSVJet[i] = -999;
		jetFlightDist2D[i] = -999;
		jetFlightDist2DError[i] = -999;
		jetFlightDist3D[i] = -999;
		jetFlightDist3DError[i] = -999;
		jetSV_x[i] = -999;
		jetSV_y[i] = -999;
		jetSV_z[i] = -999;
		jetSVNTracks[i] = 0;
		jetSVMass[i] = -999;
		jetAllMuonPt[i] = 0.0;
		jetAllMuonEta[i] = 0.0;
		jetAllMuonPhi[i] = 0.0;
		jetAllMuonM[i] = 0.0;
		jetPtWeightedDZ[i] = 0.0;
		jetNRechits[i] = 0;
		jetRechitE[i] = 0.0;
		jetRechitT[i] = 0.0;
		jetRechitT_rms[i] = 0.0;

		jetRechitE_Error[i] = 0.0;
		jetRechitT_Error[i] = 0.0;
		jetGammaMax[i] = -99.0;
		jetGammaMax_EM[i] = -99.0;
		jetGammaMax_Hadronic[i] = -99.0;
		jetGammaMax_ET[i] = -99.0;
		jetAlphaMax[i] = -99.0;
		jetBetaMax[i] = -99.0;
		jetPtAllTracks[i] = -99.0;
		jetPtAllPVTracks[i] = -99.0;
		jetMedianTheta2D[i] = -99.0;
		jetMedianIP[i] = -99.0;
		jetMinDeltaRAllTracks[i] =-99.0;
		jetMinDeltaRPVTracks[i] = -99.0;
		jet_sig_et1[i] = -99.0;
		jet_sig_et2[i] = -99.0;
		jet_energy_frac[i] = 0.0;

		jetGammaMax_wp[i] = -99.0;
		jetGammaMax_EM_wp[i] = -99.0;
		jetGammaMax_Hadronic_wp[i] = -99.0;
		jetGammaMax_ET_wp[i] = -99.0;
		jetAlphaMax_wp[i] = -99.0;
		jetBetaMax_wp[i] = -99.0;
		jetPtAllTracks_wp[i] = -99.0;
		jetPtAllPVTracks_wp[i] = -99.0;
		jetMedianTheta2D_wp[i] = -99.0;
		jetMedianIP_wp[i] = -99.0;
		jetMinDeltaRAllTracks_wp[i] =-99.0;
		jetMinDeltaRPVTracks_wp[i] = -99.0;

		jetTime[i]   = -999.;
		ecalNRechits[i] = -999;
		ecalRechitE[i] = -999.;

		jetChargedEMEnergyFraction[i] = -999.;
		jetNeutralEMEnergyFraction[i] = -999.;

		jetEcalE[i] = -999.;
		jetHcalE[i] = -999.;

		jetChargedMultiplicity[i] = 0;
		//jetNHits[i] = 0;
		//jetNPixelHits[i] = 0;
		jetNPixelHitsMedian[i] = 0;
		jetNHitsMedian[i] = 0;
	}

	//triggers
	for(int i = 0; i <NTriggersMAX; i++){
		HLTDecision[i] = false;
	}

	//Tracks
	nTracks=0;
	for(int i = 0; i <RECHITARRAYSIZE; i++){
		track_Pt[i] = -999;
		track_Eta[i] = -999;
		track_Phi[i] = -999;

		track_charge[i] = 0;
		track_bestVertexIndex[i] = 0;
		track_nMissingInnerHits[i] = 0;
		track_nMissingOuterHits[i] = 0;
		track_nPixelHits[i] = 0;
		track_nHits[i] = 0; 

		track_angle[i] = -999;
		track_dxyToBS[i] = -999;
		track_dxyErr[i] = -999;
		track_dzToPV[i] = -999;
		track_dzErr[i] = -999;
		track_chi2[i] = -999;

		track_ndof[i] = 0;
	}

	//MC
	nGenJets = 0;
	for ( int i = 0; i < OBJSIZE; i++ )
	{
		genJetE[i] = -999.;
		genJetPt[i] = -999.;
		genJetEta[i] = -999.;
		genJetPhi[i] = -999.;
		genJetMET[i] = -999.;
		//genJet_match_jet_index[i] = 666;
		//genJet_min_delta_r_match_jet[i] = -666.;
	}

	genMetPtCalo  = -999.;
	genMetPhiCalo = -999.;
	genMetPtTrue  = -999.;
	genMetPhiTrue = -999.;

	genVertexX = -999.;
	genVertexY = -999.;
	genVertexZ = -999.;
	genVertexT = -999.;

	genWeight = -999.;
	/*
	   genSignalProcessID = 0;
	   genQScale = -999.;
	   genAlphaQCD = -999.;
	   genAlphaQED = -999.;

	   scaleWeights->clear();
	   pdfWeights->clear();
	   alphasWeights->clear();
	   */
	//gen info
	for(int i = 0; i < GENPARTICLEARRAYSIZE; i++)
	{
		gParticleMotherId[i] = -99999;
		gParticleMotherIndex[i] = -99999;
		gParticleId[i] = -99999;
		gParticleStatus[i] = -99999;
		gParticleE[i] = -99999.0;
		gParticlePt[i] = -99999.0;
		gParticlePx[i] = -99999.0;
		gParticlePy[i] = -99999.0;
		gParticlePz[i] = -99999.0;
		gParticleEta[i] = -99999.0;
		gParticlePhi[i] = -99999.0;

		gParticleProdVertexX[i] = -99999.0;
		gParticleProdVertexY[i] = -99999.0;
		gParticleProdVertexZ[i] = -99999.0;
		gParticleDecayVertexX[i] = -99999.0;
		gParticleDecayVertexY[i] = -99999.0;
		gParticleDecayVertexZ[i] = -99999.0;

	}
	for ( int i = 0; i < LLP_ARRAY_SIZE; i++ )
	{
		gLLP_prod_vertex_x[i] = -666.;
		gLLP_prod_vertex_y[i] = -666.;
		gLLP_prod_vertex_z[i] = -666.;
		gLLP_decay_vertex_x[i] = -666.;
		gLLP_decay_vertex_y[i] = -666.;
		gLLP_decay_vertex_z[i] = -666.;
		gLLP_beta[i] = -666.;
		gLLP_pt[i] = -666.;
		gLLP_eta[i] = -666.;
		gLLP_e[i] = -666.;
		gLLP_phi[i] = -666.;
		gLLP_csc[i] = false;
		gLLP_dt[i] = false;
		gLLP_travel_time[i] = -666.;
	}

	for ( int i = 0; i < LLP_DAUGHTER_ARRAY_SIZE; i++ )
	{
		gLLP_daughter_id[i] = 0;
		gLLP_daughter_pt[i] = -666.;
		gLLP_daughter_eta[i] = -666.;
		gLLP_daughter_phi[i] = -666.;
		gLLP_daughter_eta_ecalcorr[i] = -666.;
		gLLP_daughter_phi_ecalcorr[i] = -666.;
		gLLP_daughter_e[i] = -666.;
		gLLP_daughter_mass[i] = -666.;

		gLLP_daughter_travel_time[i] = -666.;
		gen_time[i] = -666.;
		gen_time_pv[i] = -666.;
		photon_travel_time[i] = -666.;
		photon_travel_time_pv[i] = -666.;
	}

	//grandaughters
	for ( int i = 0; i < LLP_GRAND_DAUGHTER_ARRAY_SIZE; i++ ) {
		// gLLP_grandaughter_EB[i] = false;
		// gLLP_grandaughter_ETL[i] = false;

		// gLLP_grandaughter_photon_travel_time_EB[i] = -666.;
		// gLLP_grandaughter_photon_travel_time_ETL[i] = -666.;

		// gLLP_grandaughter_travel_time_EB[i] = -666.;
		// gLLP_grandaughter_travel_time_ETL[i] = -666.;

		// gen_time_grandaughter_EB[i] = -666.;
		// gen_time_grandaughter_ETL[i] = -666.;

		gLLP_grandaughter_id[i] = 0;
		gLLP_grandaughter_pt[i] = -666.;
		gLLP_grandaughter_eta[i] = -666.;
		gLLP_grandaughter_phi[i] = -666.;
		gLLP_grandaughter_eta_ecalcorr[i] = -666.;
		gLLP_grandaughter_phi_ecalcorr[i] = -666.;
		gLLP_grandaughter_e[i] = -666.;
		gLLP_grandaughter_mass[i] = -666.;
	}


};

//SetBranchAddress
void SusyLLPPassFailTree::InitTree()
{
	assert(tree_);
	InitVariables();

	//event info
	tree_->SetBranchAddress("runNum",      &runNum);
	tree_->SetBranchAddress("lumiSec",     &lumiSec);
	tree_->SetBranchAddress("evtNum",      &evtNum);
	tree_->SetBranchAddress("category",    &category);
	tree_->SetBranchAddress("npv",         &npv);
	tree_->SetBranchAddress("npu",         &npu);
	tree_->SetBranchAddress("pileupWeight",      &pileupWeight);
	tree_->SetBranchAddress("pileupWeightUp",      &pileupWeightUp);
	tree_->SetBranchAddress("pileupWeightDown",      &pileupWeightDown);
	tree_->SetBranchAddress("weight",      &weight);
	tree_->SetBranchAddress("rho",         &rho);
	tree_->SetBranchAddress("met",         &met);
	tree_->SetBranchAddress("metPhi",      &metPhi);

	// met filters
	tree_->SetBranchAddress("Flag2_globalSuperTightHalo2016Filter", &Flag2_globalSuperTightHalo2016Filter);
	tree_->SetBranchAddress("Flag2_globalTightHalo2016Filter", &Flag2_globalTightHalo2016Filter);
	tree_->SetBranchAddress("Flag2_goodVertices", &Flag2_goodVertices);
	tree_->SetBranchAddress("Flag2_BadChargedCandidateFilter", &Flag2_BadChargedCandidateFilter);
	tree_->SetBranchAddress("Flag2_BadPFMuonFilter", &Flag2_BadPFMuonFilter);
	tree_->SetBranchAddress("Flag2_EcalDeadCellTriggerPrimitiveFilter", &Flag2_EcalDeadCellTriggerPrimitiveFilter);
	tree_->SetBranchAddress("Flag2_HBHENoiseFilter", &Flag2_HBHENoiseFilter);
	tree_->SetBranchAddress("Flag2_HBHEIsoNoiseFilter", &Flag2_HBHEIsoNoiseFilter);
	tree_->SetBranchAddress("Flag2_ecalBadCalibFilter", &Flag2_ecalBadCalibFilter);
	tree_->SetBranchAddress("Flag2_eeBadScFilter", &Flag2_eeBadScFilter);

	//Leptons
	tree_->SetBranchAddress("nLeptons",    &nLeptons);
	tree_->SetBranchAddress("lepE",        lepE);
	tree_->SetBranchAddress("lepPt",       lepPt);
	tree_->SetBranchAddress("lepEta",      lepEta);
	tree_->SetBranchAddress("lepPhi",      lepPhi);
	tree_->SetBranchAddress("lepPdgId",  lepPdgId);
	tree_->SetBranchAddress("lepDZ",     lepDZ);
	// tree_->SetBranchAddress("lepLoosePassId", lepLoosePassId);
	// tree_->SetBranchAddress("lepMediumPassId", lepMediumPassId);
	// tree_->SetBranchAddress("lepTightPassId", lepTightPassId);
	tree_->SetBranchAddress("lepPassId", lepPassId);
	//tree_->SetBranchAddress("MT",    &MT);

	//Z-candidate
	tree_->SetBranchAddress("ZMass",       &ZMass);
	tree_->SetBranchAddress("ZPt",         &ZPt);
	tree_->SetBranchAddress("ZEta",        &ZEta);
	tree_->SetBranchAddress("ZPhi",        &ZPhi);
	tree_->SetBranchAddress("ZleptonIndex1", &ZleptonIndex1);
	tree_->SetBranchAddress("ZleptonIndex2", &ZleptonIndex2);
	tree_->SetBranchAddress("MT", &MT);

	//HH-candidate
	tree_->SetBranchAddress("dr_max", &dr_max);
	tree_->SetBranchAddress("dm_hh", &dm_hh);
	tree_->SetBranchAddress("avg_mh", &avg_mh);
	tree_->SetBranchAddress("mh1", &mh1);
	tree_->SetBranchAddress("mh2", &mh2);

	/*
	   tree_->SetBranchAddress("nElectrons",     &nElectrons);
	   tree_->SetBranchAddress("nPhotons",     &nPhotons);
	   tree_->SetBranchAddress("nTaus",     &nTaus);
	   tree_->SetBranchAddress("nMuons",     &nMuons);
	   */


	//muons
	tree_->SetBranchAddress("nMuons",             &nMuons); 
	tree_->SetBranchAddress("muonE",              muonE); 
	tree_->SetBranchAddress("muonPt",             muonPt); 
	tree_->SetBranchAddress("muonEta",            muonEta); 
	tree_->SetBranchAddress("muonPhi",            muonPhi); 
	tree_->SetBranchAddress("muonCharge",         muonCharge); 
	tree_->SetBranchAddress("muonIsLoose",        muonIsLoose); 
	tree_->SetBranchAddress("muonIsMedium",       muonIsMedium); 
	tree_->SetBranchAddress("muonIsTight",        muonIsTight); 
	tree_->SetBranchAddress("muon_pileupIso",     muon_pileupIso); 
	tree_->SetBranchAddress("muon_chargedIso",    muon_chargedIso); 
	tree_->SetBranchAddress("muon_photonIso",     muon_photonIso); 
	tree_->SetBranchAddress("muon_neutralHadIso", muon_neutralHadIso); 

	//electrons
	tree_->SetBranchAddress("nElectrons",                  &nElectrons); 
	tree_->SetBranchAddress("eleE",                        eleE); 
	tree_->SetBranchAddress("elePt",                       elePt); 
	tree_->SetBranchAddress("eleEta",                      eleEta); 
	tree_->SetBranchAddress("elePhi",                      elePhi); 
	tree_->SetBranchAddress("eleCharge",                   eleCharge); 
	tree_->SetBranchAddress("ele_passCutBasedIDVeto",      ele_passCutBasedIDVeto); 
	tree_->SetBranchAddress("ele_passCutBasedIDLoose",     ele_passCutBasedIDLoose); 
	tree_->SetBranchAddress("ele_passCutBasedIDMedium",    ele_passCutBasedIDMedium); 
	tree_->SetBranchAddress("ele_passCutBasedIDTight  ",   ele_passCutBasedIDTight  ); 

	//taus
	tree_->SetBranchAddress("nTaus",                       &nTaus);
	tree_->SetBranchAddress("tauE",                        tauE);
	tree_->SetBranchAddress("tauPt",                       tauPt); 
	tree_->SetBranchAddress("tauEta",                      tauEta); 
	tree_->SetBranchAddress("tauPhi",                      tauPhi); 
	tree_->SetBranchAddress("tau_IsLoose",                 tau_IsLoose); 
	tree_->SetBranchAddress("tau_IsMedium",                tau_IsMedium); 
	tree_->SetBranchAddress("tau_IsTight",                 tau_IsTight); 
	tree_->SetBranchAddress("tau_passEleVetoLoose",        tau_passEleVetoLoose); 
	tree_->SetBranchAddress("tau_passEleVetoMedium",       tau_passEleVetoMedium);
	tree_->SetBranchAddress("tau_passEleVetoTight",        tau_passEleVetoTight); 
	tree_->SetBranchAddress("tau_passMuVetoLoose",         tau_passMuVetoLoose); 
	tree_->SetBranchAddress("tau_passMuVetoMedium",        tau_passMuVetoMedium); 
	tree_->SetBranchAddress("tau_passMuVetoTight  ",       tau_passMuVetoTight  );

	//photons
	tree_->SetBranchAddress("nPhotons",                   &nPhotons); 
	tree_->SetBranchAddress("phoE",                       phoE); 
	tree_->SetBranchAddress("phoPt",                      phoPt); 
	tree_->SetBranchAddress("phoEta",                     phoEta); 
	tree_->SetBranchAddress("phoPhi",                     phoPhi); 
	tree_->SetBranchAddress("pho_passMVAIDWP80",          pho_passMVAIDWP80); 
	tree_->SetBranchAddress("pho_passMVAIDWP90",          pho_passMVAIDWP90); 
	tree_->SetBranchAddress("pho_passCutBasedIDLoose",    pho_passCutBasedIDLoose); 
	tree_->SetBranchAddress("pho_passCutBasedIDMedium",   pho_passCutBasedIDMedium); 
	tree_->SetBranchAddress("pho_passCutBasedIDTight ",   pho_passCutBasedIDTight ); 

	//jets
	tree_->SetBranchAddress("nJets",     &nJets);
	tree_->SetBranchAddress("jetE", jetE);
	tree_->SetBranchAddress("jetPt", jetPt);
	tree_->SetBranchAddress("jetEta", jetEta);

	tree_->SetBranchAddress("jetEcalE", jetEcalE);
	tree_->SetBranchAddress("jetHcalE", jetHcalE);

	tree_->SetBranchAddress("jetPhi", jetPhi);
	tree_->SetBranchAddress("jetCSV", jetCSV);
	tree_->SetBranchAddress("jetCISV", jetCISV);
	tree_->SetBranchAddress("jetProbb", jetProbb);
	tree_->SetBranchAddress("jetProbc", jetProbc);
	tree_->SetBranchAddress("jetProbudsg", jetProbudsg);
	tree_->SetBranchAddress("jetProbbb", jetProbbb);
	tree_->SetBranchAddress("jetMass", jetMass);
	tree_->SetBranchAddress("jetJetArea", jetJetArea);
	tree_->SetBranchAddress("jetPileupE", jetPileupE);
	tree_->SetBranchAddress("jetPileupId", jetPileupId);
	tree_->SetBranchAddress("jetPileupIdFlag", jetPileupIdFlag);
	tree_->SetBranchAddress("jetPassIDLoose", jetPassIDLoose);
	tree_->SetBranchAddress("jetPassIDTight", jetPassIDTight);
	tree_->SetBranchAddress("jetPassMuFrac", jetPassMuFrac);
	tree_->SetBranchAddress("jetPassEleFrac", jetPassEleFrac);
	tree_->SetBranchAddress("jetPartonFlavor", jetPartonFlavor);
	tree_->SetBranchAddress("jetHadronFlavor", jetHadronFlavor);
	tree_->SetBranchAddress("jetElectronEnergyFraction", jetElectronEnergyFraction);
	tree_->SetBranchAddress("jetPhotonEnergyFraction", jetPhotonEnergyFraction);
	tree_->SetBranchAddress("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction);
	tree_->SetBranchAddress("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction);
	tree_->SetBranchAddress("jetMuonEnergyFraction", jetMuonEnergyFraction);
	tree_->SetBranchAddress("jetHOEnergyFraction", jetHOEnergyFraction);
	tree_->SetBranchAddress("jetHFHadronEnergyFraction", jetHFHadronEnergyFraction);
	tree_->SetBranchAddress("jetHFEMEnergyFraction",jetHFEMEnergyFraction);
	tree_->SetBranchAddress("jetChargedHadronMultiplicity", jetChargedHadronMultiplicity);
	tree_->SetBranchAddress("jetNeutralHadronMultiplicity", jetNeutralHadronMultiplicity);
	tree_->SetBranchAddress("jetPhotonMultiplicity", jetPhotonMultiplicity);
	tree_->SetBranchAddress("jetElectronMultiplicity", jetElectronMultiplicity);
	tree_->SetBranchAddress("jetMuonMultiplicity", jetMuonMultiplicity);
	tree_->SetBranchAddress("jetQGLikelihood",          jetQGLikelihood);
	tree_->SetBranchAddress("jetNSV",                   jetNSV);
	tree_->SetBranchAddress("jetNSVCand",               jetNSVCand);
	tree_->SetBranchAddress("jetNVertexTracks",         jetNVertexTracks);
	tree_->SetBranchAddress("jetNSelectedTracks",       jetNSelectedTracks);
	tree_->SetBranchAddress("jetDRSVJet",               jetDRSVJet);
	tree_->SetBranchAddress("jetFlightDist2D",          jetFlightDist2D);
	tree_->SetBranchAddress("jetFlightDist2DError",     jetFlightDist2DError);
	tree_->SetBranchAddress("jetFlightDist3D",          jetFlightDist3D);
	tree_->SetBranchAddress("jetFlightDist3DError",     jetFlightDist3DError );
	tree_->SetBranchAddress("jetSV_x",                  jetSV_x);
	tree_->SetBranchAddress("jetSV_y",                  jetSV_y);
	tree_->SetBranchAddress("jetSV_z",                  jetSV_z);
	tree_->SetBranchAddress("jetSVNTracks",             jetSVNTracks);
	tree_->SetBranchAddress("jetSVMass",                jetSVMass);
	tree_->SetBranchAddress("jetAllMuonPt", jetAllMuonPt);
	tree_->SetBranchAddress("jetAllMuonEta", jetAllMuonEta);
	tree_->SetBranchAddress("jetAllMuonPhi", jetAllMuonPhi);
	tree_->SetBranchAddress("jetAllMuonM", jetAllMuonM);
	tree_->SetBranchAddress("jetPtWeightedDZ", jetPtWeightedDZ);
	tree_->SetBranchAddress("jetNRechits", jetNRechits);
	tree_->SetBranchAddress("jetRechitE", jetRechitE);
	tree_->SetBranchAddress("jetRechitT", jetRechitT);
	tree_->SetBranchAddress("jetRechitT_rms", jetRechitT_rms);

	tree_->SetBranchAddress("jetRechitE_Error", jetRechitE_Error);
	tree_->SetBranchAddress("jetRechitT_Error", jetRechitT_Error);
	tree_->SetBranchAddress("jetAlphaMax",jetAlphaMax);
	tree_->SetBranchAddress("jetBetaMax",jetBetaMax);
	tree_->SetBranchAddress("jetGammaMax_ET",jetGammaMax_ET);
	tree_->SetBranchAddress("jetGammaMax_EM",jetGammaMax_EM);
	tree_->SetBranchAddress("jetGammaMax_Hadronic",jetGammaMax_Hadronic);
	tree_->SetBranchAddress("jetGammaMax",jetGammaMax);
	tree_->SetBranchAddress("jetPtAllTracks",jetPtAllTracks);
	tree_->SetBranchAddress("jetPtAllPVTracks",jetPtAllPVTracks);
	tree_->SetBranchAddress("jetMedianTheta2D",jetMedianTheta2D);
	tree_->SetBranchAddress("jetMedianIP",jetMedianIP);
	tree_->SetBranchAddress("jetMinDeltaRAllTracks",jetMinDeltaRAllTracks);
	tree_->SetBranchAddress("jetMinDeltaRPVTracks",jetMinDeltaRPVTracks);
	tree_->SetBranchAddress("jet_sig_et1",jet_sig_et1);
	tree_->SetBranchAddress("jet_sig_et2",jet_sig_et2);
	tree_->SetBranchAddress("jet_energy_frac",jet_energy_frac);

	tree_->SetBranchAddress("jetAlphaMax_wp",jetAlphaMax_wp);
	tree_->SetBranchAddress("jetBetaMax_wp",jetBetaMax_wp);
	tree_->SetBranchAddress("jetGammaMax_ET_wp",jetGammaMax_ET_wp);
	tree_->SetBranchAddress("jetGammaMax_EM_wp",jetGammaMax_EM_wp);
	tree_->SetBranchAddress("jetGammaMax_Hadronic_wp",jetGammaMax_Hadronic_wp);
	tree_->SetBranchAddress("jetGammaMax_wp",jetGammaMax_wp);
	tree_->SetBranchAddress("jetPtAllTracks_wp",jetPtAllTracks_wp);
	tree_->SetBranchAddress("jetPtAllPVTracks_wp",jetPtAllPVTracks_wp);
	tree_->SetBranchAddress("jetMedianTheta2D_wp",jetMedianTheta2D_wp);
	tree_->SetBranchAddress("jetMedianIP_wp",jetMedianIP_wp);
	tree_->SetBranchAddress("jetMinDeltaRAllTracks_wp",jetMinDeltaRAllTracks_wp);
	tree_->SetBranchAddress("jetMinDeltaRPVTracks_wp",jetMinDeltaRPVTracks_wp);

	tree_->SetBranchAddress("jetTime",   jetTime);
	tree_->SetBranchAddress("ecalNRechits",   ecalNRechits);
	tree_->SetBranchAddress("ecalRechitE", ecalRechitE);

	tree_->SetBranchAddress("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction);
	tree_->SetBranchAddress("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction);

	tree_->SetBranchAddress("jetChargedMultiplicity",   jetChargedMultiplicity);
	//tree_->SetBranchAddress("jetNHits",   jetNHits);
	//tree_->SetBranchAddress("jetNPixelHits",   jetNPixelHits);
	tree_->SetBranchAddress("jetNPixelHitsMedian",   jetNPixelHitsMedian);
	tree_->SetBranchAddress("jetNHitsMedian",   jetNHitsMedian);
	// tree_->SetBranchAddress("jetLoosePassId", jetLoosePassId);
	// tree_->SetBranchAddress("jetTightPassId", jetTightPassId);

	// triggers
	tree_->SetBranchAddress("HLTDecision",   &HLTDecision);

	//MC
	tree_->SetBranchAddress("genMetPtTrue",         &genMetPtTrue);
	tree_->SetBranchAddress("genMetPhiTrue",      &genMetPhiTrue);
	tree_->SetBranchAddress("genMetPtCalo",      &genMetPtCalo);
	tree_->SetBranchAddress("genMetPhiCalo",      &genMetPhiCalo);

	//Tracks
	tree_->SetBranchAddress("nTracks", &nTracks);
	tree_->SetBranchAddress("track_Pt", track_Pt);
	tree_->SetBranchAddress("track_Eta", track_Eta);
	tree_->SetBranchAddress("track_Phi", track_Phi);
	tree_->SetBranchAddress("track_charge", track_charge);
	tree_->SetBranchAddress("track_bestVertexIndex", track_bestVertexIndex);
	tree_->SetBranchAddress("track_nMissingInnerHits", track_nMissingInnerHits);
	tree_->SetBranchAddress("track_nMissingOuterHits", track_nMissingOuterHits);
	tree_->SetBranchAddress("track_nPixelHits", track_nPixelHits);
	tree_->SetBranchAddress("track_nHits", track_nHits); 
	tree_->SetBranchAddress("track_angle", track_angle);
	tree_->SetBranchAddress("track_dxyToBS", track_dxyToBS);
	tree_->SetBranchAddress("track_dxyErr", track_dxyErr);
	tree_->SetBranchAddress("track_dzToPV", track_dzToPV);
	tree_->SetBranchAddress("track_dzErr", track_dzErr);
	tree_->SetBranchAddress("track_chi2", track_chi2);
	tree_->SetBranchAddress("track_ndof", track_ndof);
	   tree_->SetBranchAddress("nGenParticle",      &nGenParticle);
	   tree_->SetBranchAddress("gParticleId",      &gParticleId);
	   tree_->SetBranchAddress("gParticleStatus",      &gParticleStatus);
	   tree_->SetBranchAddress("gParticleMotherId",      &gParticleMotherId);
	   tree_->SetBranchAddress("gParticleE",      &gParticleE);
	   tree_->SetBranchAddress("gParticlePt",      &gParticlePt);
	   tree_->SetBranchAddress("gParticleEta",      &gParticleEta);
	   tree_->SetBranchAddress("gParticlePhi",      &gParticlePhi);

	/*
	   tree_->SetBranchAddress("nGenJets",      &nGenJets);
	   tree_->SetBranchAddress("genJetE",      &genJetE);
	   tree_->SetBranchAddress("genJetPt",      &genJetPt);
	   tree_->SetBranchAddress("genJetEta",      &genJetEta);
	   tree_->SetBranchAddress("genJetPhi",      &genJetPhi);
	   tree_->SetBranchAddress("genJetMET",      &genJetMET);


	//gen info
	tree_->SetBranchAddress("gLLP_eta",    gLLP_eta);
	tree_->SetBranchAddress("gLLP_phi",    gLLP_phi);
	tree_->SetBranchAddress("gLLP_csc",    gLLP_csc);
	tree_->SetBranchAddress("gLLP_ctau",    gLLP_ctau);
	tree_->SetBranchAddress("gLLP_beta",    gLLP_beta);

	tree_->SetBranchAddress("gLLP_decay_vertex_r",    gLLP_decay_vertex_r);
	tree_->SetBranchAddress("gLLP_decay_vertex_x",    gLLP_decay_vertex_x);
	tree_->SetBranchAddress("gLLP_decay_vertex_y",    gLLP_decay_vertex_y);
	tree_->SetBranchAddress("gLLP_decay_vertex_z",    gLLP_decay_vertex_z);

	// gLLP
	tree_->SetBranchAddress("gLLP_travel_time",    gLLP_travel_time);
	tree_->SetBranchAddress("gLLP_e",    gLLP_e);
	tree_->SetBranchAddress("gLLP_pt",    gLLP_pt);
	tree_->SetBranchAddress("gLLP_eta",    gLLP_eta);
	tree_->SetBranchAddress("gLLP_phi",    gLLP_phi);
	tree_->SetBranchAddress("gLLP_beta",    gLLP_beta);
	tree_->SetBranchAddress("gLLP_csc",    gLLP_csc);
	tree_->SetBranchAddress("gLLP_ctau",    gLLP_ctau);

	tree_->SetBranchAddress("gLLP_decay_vertex_r",    gLLP_decay_vertex_r);
	tree_->SetBranchAddress("gLLP_decay_vertex_x",    gLLP_decay_vertex_x);
	tree_->SetBranchAddress("gLLP_decay_vertex_y",    gLLP_decay_vertex_y);
	tree_->SetBranchAddress("gLLP_decay_vertex_z",    gLLP_decay_vertex_z);

	tree_->SetBranchAddress("gLLP_prod_vertex_x",    gLLP_prod_vertex_x);
	tree_->SetBranchAddress("gLLP_prod_vertex_y",    gLLP_prod_vertex_y);
	tree_->SetBranchAddress("gLLP_prod_vertex_z",    gLLP_prod_vertex_z);
	 *
	// gLLP daughters
	tree_->SetBranchAddress("gen_time",    gen_time);
	tree_->SetBranchAddress("photon_travel_time",    photon_travel_time);
	tree_->SetBranchAddress("gLLP_daughter_travel_time",    gLLP_daughter_travel_time);
	tree_->SetBranchAddress("gLLP_daughter_e",    gLLP_daughter_e);
	tree_->SetBranchAddress("gLLP_daughter_pt",    gLLP_daughter_pt);
	tree_->SetBranchAddress("gLLP_daughter_eta",    gLLP_daughter_eta);
	tree_->SetBranchAddress("gLLP_daughter_phi",    gLLP_daughter_phi);
	tree_->SetBranchAddress("gLLP_daughter_eta_ecalcorr",    gLLP_daughter_eta_ecalcorr);
	tree_->SetBranchAddress("gLLP_daughter_phi_ecalcorr",    gLLP_daughter_phi_ecalcorr);
	tree_->SetBranchAddress("gLLP_min_delta_r_match_jet",    gLLP_min_delta_r_match_jet);
	tree_->SetBranchAddress("gLLP_daughter_match_jet_index",   gLLP_daughter_match_jet_index);
	 *  
	//daughters
	tree_->SetBranchAddress("gLLP_daughter_EB", gLLP_daughter_EB);
	tree_->SetBranchAddress("gLLP_daughter_ETL", gLLP_daughter_ETL);

	tree_->SetBranchAddress("gLLP_daughter_photon_travel_time_EB", gLLP_daughter_photon_travel_time_EB);
	tree_->SetBranchAddress("gLLP_daughter_photon_travel_time_ETL", gLLP_daughter_photon_travel_time_ETL);

	tree_->SetBranchAddress("gLLP_daughter_travel_time_EB", gLLP_daughter_travel_time_EB);
	tree_->SetBranchAddress("gLLP_daughter_travel_time_ETL", gLLP_daughter_travel_time_ETL);

	tree_->SetBranchAddress("gen_time_daughter_EB", gen_time_daughter_EB);
	tree_->SetBranchAddress("gen_time_daughter_ETL", gen_time_daughter_ETL);

	tree_->SetBranchAddress("gLLP_daughter_id", gLLP_daughter_id);
	tree_->SetBranchAddress("gLLP_daughter_pt", gLLP_daughter_pt);
	tree_->SetBranchAddress("gLLP_daughter_eta", gLLP_daughter_eta);
	tree_->SetBranchAddress("gLLP_daughter_phi", gLLP_daughter_phi);
	tree_->SetBranchAddress("gLLP_daughter_eta_ecalcorr", gLLP_daughter_eta_ecalcorr);
	tree_->SetBranchAddress("gLLP_daughter_phi_ecalcorr", gLLP_daughter_phi_ecalcorr);
	tree_->SetBranchAddress("gLLP_daughter_e", gLLP_daughter_e);
	tree_->SetBranchAddress("gLLP_daughter_mass", gLLP_daughter_mass);

	tree_->SetBranchAddress("gLLP_daughter_match_jet_index", gLLP_daughter_match_jet_index);
	tree_->SetBranchAddress("gLLP_daughter_min_delta_r_match_jet", gLLP_daughter_min_delta_r_match_jet);

	//grandaughters
	tree_->SetBranchAddress("gLLP_grandaughter_EB", gLLP_grandaughter_EB);
	tree_->SetBranchAddress("gLLP_grandaughter_ETL", gLLP_grandaughter_ETL);

	tree_->SetBranchAddress("gLLP_grandaughter_photon_travel_time_EB", gLLP_grandaughter_photon_travel_time_EB);
	tree_->SetBranchAddress("gLLP_grandaughter_photon_travel_time_ETL", gLLP_grandaughter_photon_travel_time_ETL);

	tree_->SetBranchAddress("gLLP_grandaughter_travel_time_EB", gLLP_grandaughter_travel_time_EB);
	tree_->SetBranchAddress("gLLP_grandaughter_travel_time_ETL", gLLP_grandaughter_travel_time_ETL);

	tree_->SetBranchAddress("gen_time_grandaughter_EB", gen_time_grandaughter_EB);
	tree_->SetBranchAddress("gen_time_grandaughter_ETL", gen_time_grandaughter_ETL);

	tree_->SetBranchAddress("gLLP_grandaughter_id", gLLP_grandaughter_id);
	tree_->SetBranchAddress("gLLP_grandaughter_pt", gLLP_grandaughter_pt);
	tree_->SetBranchAddress("gLLP_grandaughter_eta", gLLP_grandaughter_eta);
	tree_->SetBranchAddress("gLLP_grandaughter_phi", gLLP_grandaughter_phi);
	tree_->SetBranchAddress("gLLP_grandaughter_eta_ecalcorr", gLLP_grandaughter_eta_ecalcorr);
	tree_->SetBranchAddress("gLLP_grandaughter_phi_ecalcorr", gLLP_grandaughter_phi_ecalcorr);
	tree_->SetBranchAddress("gLLP_grandaughter_e", gLLP_grandaughter_e);
	tree_->SetBranchAddress("gLLP_grandaughter_mass", gLLP_grandaughter_mass);

	tree_->SetBranchAddress("gLLP_grandaughter_match_jet_index", gLLP_grandaughter_match_jet_index);
	tree_->SetBranchAddress("gLLP_grandaughter_min_delta_r_match_jet", gLLP_grandaughter_min_delta_r_match_jet);
	*/
};

//LoadTree
void SusyLLPPassFailTree::LoadTree(const char* file)
{
	f_ = TFile::Open(file);
	assert(f_);
	tree_ = dynamic_cast<TTree*>(f_->Get("SusyLLPPassFailTree"));
	InitTree();
	assert(tree_);
};

//CreateTree
void SusyLLPPassFailTree::CreateTree()
{
	//tree
	tree_ = new TTree("SusyLLPPassFailTree","SusyLLPPassFailTree");
	f_ = 0;

	//event info
	//tree_->Branch("runNum",      &runNum,     "runNum/i");      // event run number
	//tree_->Branch("lumiSec",     &lumiSec,    "lumiSec/i");     // event lumi section
	tree_->Branch("evtNum",      &evtNum,     "evtNum/i");      // event number
	//tree_->Branch("category",    &category,   "category/i");    // dilepton category
	//tree_->Branch("npv",         &npv,        "npv/i");         // number of primary vertices
	//tree_->Branch("npu",         &npu,        "npu/i");         // number of in-time PU events (MC)
	tree_->Branch("pileupWeight",      &pileupWeight,     "pileupWeight/F");
	//tree_->Branch("pileupWeightUp",      &pileupWeightUp,     "pileupWeightUp/F");
	//tree_->Branch("pileupWeightDown",      &pileupWeightDown,     "pileupWeightDown/F");
	tree_->Branch("weight",      &weight,     "weight/F");

	//tree_->Branch("rho",         &rho,        "rho/F");
	//tree_->Branch("met",         &met,        "met/F");         // MET
	//tree_->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)
	//tree_->Branch("HT",    &HT,    "HT/F");

	//tree_->Branch("jetMet_dPhi",    &jetMet_dPhi,    "jetMet_dPhi/F");
	//tree_->Branch("jetMet_dPhiStar",    &jetMet_dPhiStar,    "jetMet_dPhiStar/F");
	//tree_->Branch("jetMet_dPhiMin",    &jetMet_dPhiMin,    "jetMet_dPhiMin/F");
	//tree_->Branch("jetMet_dPhiStarMin",    &jetMet_dPhiStarMin,    "jetMet_dPhiStarMin/F");
	//tree_->Branch("jetMet_dPhiMin4",    &jetMet_dPhiMin4,    "jetMet_dPhiMin4/F");

	//tree_->Branch("jetPho_dPhi",    &jetPho_dPhi,    "jetPho_dPhi/F");
	//tree_->Branch("jetPho_dPhiStar",    &jetPho_dPhiStar,    "jetPho_dPhiStar/F");
	//tree_->Branch("jetPho_dPhiMin",    &jetPho_dPhiMin,    "jetPho_dPhiMin/F");
	//tree_->Branch("jetPho_dPhiStarMin",    &jetPho_dPhiStarMin,    "jetPho_dPhiStarMin/F");
	//tree_->Branch("jetPho_dPhiMin4",    &jetPho_dPhiMin4,    "jetPho_dPhiMin4/F");

	//tree_->Branch("jet2_dPhi",    &jet2_dPhi,    "jet2_dPhi/F");

	// met filters
	//tree_->Branch("Flag2_globalSuperTightHalo2016Filter", &Flag2_globalSuperTightHalo2016Filter, "Flag2_globalSuperTightHalo2016Filter/O");
	//tree_->Branch("Flag2_globalTightHalo2016Filter", &Flag2_globalTightHalo2016Filter, "Flag2_globalTightHalo2016Filter/O");
	//tree_->Branch("Flag2_goodVertices", &Flag2_goodVertices, "Flag2_goodVertices/O");
	//tree_->Branch("Flag2_BadChargedCandidateFilter", &Flag2_BadChargedCandidateFilter, "Flag2_BadChargedCandidateFilter/O");
	//tree_->Branch("Flag2_BadPFMuonFilter", &Flag2_BadPFMuonFilter, "Flag2_BadPFMuonFilter/O");
	//tree_->Branch("Flag2_EcalDeadCellTriggerPrimitiveFilter", &Flag2_EcalDeadCellTriggerPrimitiveFilter, "Flag2_EcalDeadCellTriggerPrimitiveFilter/O");
	//tree_->Branch("Flag2_HBHENoiseFilter", &Flag2_HBHENoiseFilter, "Flag2_HBHENoiseFilter/O");
	//tree_->Branch("Flag2_HBHEIsoNoiseFilter", &Flag2_HBHEIsoNoiseFilter, "Flag2_HBHEIsoNoiseFilter/O");
	//tree_->Branch("Flag2_ecalBadCalibFilter", &Flag2_ecalBadCalibFilter, "Flag2_ecalBadCalibFilter/O");
	//tree_->Branch("Flag2_eeBadScFilter", &Flag2_eeBadScFilter, "Flag2_eeBadScFilter/O");

	//leptons
	tree_->Branch("nLeptons",  &nLeptons, "nLeptons/I");
	//tree_->Branch("lepPt",     lepPt,     "lepPt[nLeptons]/F");

	//Z-candidate
	//tree_->Branch("MT",      &MT,  "MT/F");
	//tree_->Branch("ZMass",      &ZMass,  "ZMass/F");
	//tree_->Branch("ZPt",        &ZPt,    "ZPt/F");

	//muons
	tree_->Branch("nMuons", &nMuons,"nMuons/I");
	//tree_->Branch("muonPt", muonPt,"muonPt[nMuons]/F");

	//electrons
	tree_->Branch("nElectrons", &nElectrons,"nElectrons/I");
	//tree_->Branch("elePt", elePt,"elePt[nElectrons]/F");

	//taus
	//tree_->Branch("nTaus", &nTaus,"nTaus/I");
	//tree_->Branch("tauPt", tauPt,"tauPt[nTaus]/F");

	//photons
	tree_->Branch("nPhotons", &nPhotons,"nPhotons/I");
	//tree_->Branch("phoPt", phoPt,"phoPt[nPhotons]/F");

	//jets
	tree_->Branch("nJets", &nJets,"nJets/I");
	tree_->Branch("jetPt", jetPt,"jetPt[nJets]/F");
	tree_->Branch("jetPass", jetPass,"jetPass[nJets]/I");

	//tree_->Branch("jetChargeHadronEnergyFraction",   jetChargeHadronEnergyFraction,   "jetChargeHadronEnergyFraction[nJets]/F");
	//tree_->Branch("jetTime",   jetTime,   "jetTime[nJets]/F");
	//tree_->Branch("jetGammaMax_ET",jetGammaMax_ET,"jetGammaMax_ET[nJets]/F");
	//tree_->Branch("jetMinDeltaRPVTracks",jetMinDeltaRPVTracks,"jetMinDeltaRPVTracks[nJets]/F");


	//HLT
	//tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[1201]/O"); //hardcoded
	//tree_->Branch("HLTPrescale", HLTPrescale, "HLTPrescale[1201]/I"); //hardcoded


	//MC

};
