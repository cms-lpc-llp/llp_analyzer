#include "RazorHelper.h"
#include "SusyLLPTree.h"
#include "assert.h"
#include "TTree.h"

#define N_MAX_LLPS 2
#define N_MAX_LLP_DAUGHTERS 4
#define N_MAX_LLP_GRAND_DAUGHTERS 4
#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_MUONS 100
#define N_MAX_ELECTRONS 100
#define N_MAX_TAUS 100
#define N_MAX_PHOTONS 100
#define NTriggersMAX 602 //Number of trigger in the .dat file

// Constructor
SusyLLPTree::SusyLLPTree()
{
  InitVariables();
};

SusyLLPTree::~SusyLLPTree()
{
  if (f_) f_->Close();
};

//Initialize
void SusyLLPTree::InitVariables()
{ 
    //evnt info
    runNum=0; lumiSec=0; evtNum=0; category=0;
    npv=0; npu=0; rho=-1; weight=-1;
    pileupWeight = 1; pileupWeightUp = 1; pileupWeightDown = 1;
    met=-1; metPhi=-1;
    HT=-1;
    jetMet_dPhi=-999.; jetMet_dPhiMin=999.; jetMet_dPhiMin4=999.;
    jetMet_dPhiStar=-999.; jetMet_dPhiStarMin=999.; 

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
/*
      jetE[i]      = -999.;
      jetEt[i] = -999.;
      jetPt[i]     = -999.;
      jetEta[i]    = -999.;
      jetPhi[i]    = -999.;
      jetTime[i]   = -999.;
      // jetLoosePassId[i] = false;
      jetPassId[i] = false;
      jetCSVT[i] = false;
      jetCISV[i] = -1.;
      matched[i] = false;
      jet_matched_gLLP0_daughter[i] = false;
      jet_matched_gLLP1_daughter[i] = false;
      jet_matched_gLLP0_grandaughter[i] = false;
      jet_matched_gLLP1_grandaughter[i] = false;
      jet_energy_frac[i] = -999.;
      jet_sig_et1[i] = -999.;
      jet_sig_et2[i] = -999.;
      ecalNRechits[i] = -999;
      ecalRechitE[i] = -999.;
      jetChargedEMEnergyFraction[i] = -999.;
      jetNeutralEMEnergyFraction[i] = -999.;
      jetChargedHadronEnergyFraction[i] = -999.;
      jetNeutralHadronEnergyFraction[i] = -999.;
      jetHoverE[i] = -999.;
      jetGammaMax_ET[i] = -999.;
      jetMinDeltaRPVTracks[i] = -999.;
      jetPtAllPVTracks[i] = -999.;
      jetMinDeltaRAllTracks[i] = -999.;
      jetPtAllTracks[i] = -999.;

      jetGammaMax_ET_wp[i] = -999.;
      jetMinDeltaRPVTracks_wp[i] = -999.;
      jetPtAllPVTracks_wp[i] = -999.;
      jetMinDeltaRAllTracks_wp[i] = -999.;
      jetPtAllTracks_wp[i] = -999.;
*/
/*
      calojetE[i] = -999.;
      calojetEt[i] = -999.;
      calojetPt[i] = -999.;
      calojetEta[i] = -999.;
      calojetPhi[i] = -999.;
      calojetTime[i] = -999.;
      calojetPassId[i] = -999.;
      calojet_EMEnergyFraction[i] = -999.;
      calojet_HadronicEnergyFraction[i] = -999.;
      //calojetChargedEMEnergyFraction[i] = -999.;
      //calojetNeutralEMEnergyFraction[i] = -999.;
      //calojetChargedHadronEnergyFraction[i] = -999.;
      //calojetNeutralHadronEnergyFraction[i] = -999.;
      calojetGammaMax_ET[i] = -999.;
      calojetMinDeltaRPVTracks[i] = -999.;
      calojetPtAllPVTracks[i] = -999.;
      calojetNRechits[i] = -999.;
      calojetRechitE[i] = -999.;
*/
    }
    
    //triggers
    for(int i = 0; i <NTriggersMAX; i++){
      HLTDecision[i] = false;
    }

    //gLLP
    gLLP0_EB = false;
    gLLP1_EB = false;
    gLLP0_EE = false;
    gLLP1_EE = false;

    for(int i = 0; i <N_MAX_LLPS; i++){
      gLLP_travel_time[i] = -999.;
      gLLP_e[i] = -999.;
      gLLP_pt[i] = -999.;
      gLLP_eta[i] = -999.;
      gLLP_beta[i] = -999.;
      gLLP_phi[i] = -999.;
      gLLP_decay_vertex_r[i] = -999.;
      gLLP_decay_vertex_x[i] = -999.;
      gLLP_decay_vertex_y[i] = -999.;
      gLLP_decay_vertex_z[i] = -999.;
      gLLP_prod_vertex_x[i] = -999.;
      gLLP_prod_vertex_y[i] = -999.;
      gLLP_prod_vertex_z[i] = -999.;

    }
      gLLP_dr = -999.;
      gLLP_grandaughter_dr1 = -999.;
      gLLP_grandaughter_dr2 = -999.;
      gLLP_grandaughter_matched_jet_dr1 = -999.; 
      gLLP_grandaughter_matched_jet_dr2 = -999.; 
/*
    //gLLP daughters
    for(int i = 0; i <N_MAX_LLP_DAUGHTERS; i++){
      gen_time[i] = -999.;
      photon_travel_time[i] = -999.;
      gLLP_daughter_travel_time[i] = -999.;
      gLLP_daughter_e[i] = -999.;
      gLLP_daughter_pt[i] = -999.;
      gLLP_daughter_eta[i] = -999.;
      gLLP_daughter_phi[i] = -999.;
      gLLP_daughter_eta_ecalcorr[i] = -999.;
      gLLP_daughter_phi_ecalcorr[i] = -999.;
      gLLP_min_delta_r_match_jet[i] = -999.;
      gLLP_daughter_match_jet_index[i] = 999;

    }
*/
    //gLLP daughters
    for(int i = 0; i <N_MAX_LLP_DAUGHTERS; i++){
    gLLP_daughter_EB[i] = false;
    gLLP_daughter_ETL[i] = false;

    gLLP_daughter_photon_travel_time_EB[i] = -666.;
    gLLP_daughter_photon_travel_time_ETL[i] = -666.;

    gLLP_daughter_travel_time_EB[i] = -666.;
    gLLP_daughter_travel_time_ETL[i] = -666.;

    gen_time_daughter_EB[i] = -666.;
    gen_time_daughter_ETL[i] = -666.;

    gLLP_daughter_id[i] = 0;
    gLLP_daughter_pt[i] = -666.;
    gLLP_daughter_eta[i] = -666.;
    gLLP_daughter_phi[i] = -666.;
    gLLP_daughter_eta_ecalcorr[i] = -666.;
    gLLP_daughter_phi_ecalcorr[i] = -666.;
    gLLP_daughter_e[i] = -666.;
    gLLP_daughter_mass[i] = -666.;

    //gLLP_daughter_match_calojet_index[i] = -666;
    gLLP_daughter_match_jet_index[i] = -666;
    //gLLP_daughter_min_delta_r_match_calojet[i] = 666.;
    gLLP_daughter_min_delta_r_match_jet[i] = 666.;
    }
 //grandaughters
  for(int i = 0; i <N_MAX_LLP_GRAND_DAUGHTERS; i++)
  {
    gLLP_grandaughter_EB[i] = false;
    gLLP_grandaughter_ETL[i] = false;

    gLLP_grandaughter_photon_travel_time_EB[i] = -666.;
    gLLP_grandaughter_photon_travel_time_ETL[i] = -666.;

    gLLP_grandaughter_travel_time_EB[i] = -666.;
    gLLP_grandaughter_travel_time_ETL[i] = -666.;

    gen_time_grandaughter_EB[i] = -666.;
    gen_time_grandaughter_ETL[i] = -666.;

    gLLP_grandaughter_id[i] = 0;
    gLLP_grandaughter_pt[i] = -666.;
    gLLP_grandaughter_eta[i] = -666.;
    gLLP_grandaughter_phi[i] = -666.;
    gLLP_grandaughter_eta_ecalcorr[i] = -666.;
    gLLP_grandaughter_phi_ecalcorr[i] = -666.;
    gLLP_grandaughter_e[i] = -666.;
    gLLP_grandaughter_mass[i] = -666.;

    //gLLP_grandaughter_match_calojet_index[i] = -666;
    gLLP_grandaughter_match_jet_index[i] = -666;
    //gLLP_grandaughter_min_delta_r_match_calojet[i] = 666.;
    gLLP_grandaughter_min_delta_r_match_jet[i] = 666.;

  }

};

//SetBranchAddress
void SusyLLPTree::InitTree()
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
    /*
    tree_->SetBranchAddress("nJets",     &nJets);
    tree_->SetBranchAddress("nBJets",     &nBJets);
    tree_->SetBranchAddress("jetE",      jetE);
    tree_->SetBranchAddress("jetEt",      jetEt);

    tree_->SetBranchAddress("jetPt",     jetPt);
    tree_->SetBranchAddress("jetEta",    jetEta);
    tree_->SetBranchAddress("jetPhi",    jetPhi);
    tree_->SetBranchAddress("jetTime",   jetTime);
    tree_->SetBranchAddress("jetPassId", jetPassId);
    tree_->SetBranchAddress("jetCSVT", &jetCSVT);
    tree_->SetBranchAddress("jetCISV", jetCISV);
    tree_->SetBranchAddress("matched", matched);
    tree_->SetBranchAddress("jet_matched_gLLP0_daughter", jet_matched_gLLP0_daughter);
    tree_->SetBranchAddress("jet_matched_gLLP1_daughter", jet_matched_gLLP1_daughter);
    tree_->SetBranchAddress("jet_matched_gLLP0_grandaughter", jet_matched_gLLP0_grandaughter);
    tree_->SetBranchAddress("jet_matched_gLLP1_grandaughter", jet_matched_gLLP1_grandaughter);
    tree_->SetBranchAddress("jet_energy_frac", jet_energy_frac);
    tree_->SetBranchAddress("jet_sig_et1", jet_sig_et1);
    tree_->SetBranchAddress("jet_sig_et2", jet_sig_et2);

    tree_->SetBranchAddress("ecalNRechits",   ecalNRechits);
    tree_->SetBranchAddress("ecalRechitE", ecalRechitE);
    tree_->SetBranchAddress("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction);
    tree_->SetBranchAddress("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction);
    tree_->SetBranchAddress("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction);
    tree_->SetBranchAddress("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction);
    tree_->SetBranchAddress("jetGammaMax_ET", jetGammaMax_ET);
    tree_->SetBranchAddress("jetMinDeltaRPVTracks", jetMinDeltaRPVTracks);
    tree_->SetBranchAddress("jetPtAllPVTracks", jetPtAllPVTracks);
    tree_->SetBranchAddress("jetMinDeltaRAllTracks", jetMinDeltaRAllTracks);
    tree_->SetBranchAddress("jetPtAllTracks", jetPtAllTracks);

    tree_->SetBranchAddress("jetGammaMax_ET_wp", jetGammaMax_ET_wp);
    tree_->SetBranchAddress("jetMinDeltaRPVTracks_wp", jetMinDeltaRPVTracks_wp);
    tree_->SetBranchAddress("jetPtAllPVTracks_wp", jetPtAllPVTracks_wp);
    tree_->SetBranchAddress("jetMinDeltaRAllTracks_wp", jetMinDeltaRAllTracks_wp);
    tree_->SetBranchAddress("jetPtAllTracks_wp", jetPtAllTracks_wp);
*/
/*
    //calojets
    tree_->SetBranchAddress("nCaloJets",     &nCaloJets);
    tree_->SetBranchAddress("calojetE",      calojetE);
    tree_->SetBranchAddress("calojetEt",      calojetEt);
    tree_->SetBranchAddress("calojetPt",     calojetPt);
    tree_->SetBranchAddress("calojetEta",    calojetEta);
    tree_->SetBranchAddress("calojetPhi",    calojetPhi);
    tree_->SetBranchAddress("calojetTime",   calojetTime);
    tree_->SetBranchAddress("calojetPassId", calojetPassId);
    tree_->SetBranchAddress("calojetRechitE",   calojetRechitE);
    // tree_->SetBranchAddress("calojetRechitT_rms", calojetRechitT_rms);
    tree_->SetBranchAddress("calojet_EMEnergyFraction", calojet_EMEnergyFraction);
    tree_->SetBranchAddress("calojet_HadronicEnergyFraction", calojet_HadronicEnergyFraction);
    //tree_->SetBranchAddress("calojetChargedEMEnergyFraction", calojetChargedEMEnergyFraction);
    //tree_->SetBranchAddress("calojetNeutralEMEnergyFraction", calojetNeutralEMEnergyFraction);
    //tree_->SetBranchAddress("calojetChargedHadronEnergyFraction", calojetChargedHadronEnergyFraction);
    //tree_->SetBranchAddress("calojetNeutralHadronEnergyFraction", calojetNeutralHadronEnergyFraction);

    tree_->SetBranchAddress("calojetGammaMax_ET", calojetGammaMax_ET);
    tree_->SetBranchAddress("calojetMinDeltaRPVTracks", calojetMinDeltaRPVTracks);
    tree_->SetBranchAddress("calojetPtAllPVTracks", calojetPtAllPVTracks);
*/  
    // tree_->SetBranchAddress("jetLoosePassId", jetLoosePassId);
    // tree_->SetBranchAddress("jetTightPassId", jetTightPassId);

    // triggers
    tree_->SetBranchAddress("HLTDecision",   &HLTDecision);

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
/*
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
*/  
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
 //tree_->SetBranchAddress("gLLP_daughter_match_calojet_index", gLLP_daughter_match_calojet_index);
 //tree_->SetBranchAddress("gLLP_daughter_min_delta_r_match_calojet", gLLP_daughter_min_delta_r_match_calojet);

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
 //tree_->SetBranchAddress("gLLP_grandaughter_match_calojet_index", gLLP_grandaughter_match_calojet_index);
 //tree_->SetBranchAddress("gLLP_grandaughter_min_delta_r_match_calojet", gLLP_grandaughter_min_delta_r_match_calojet);

};

//LoadTree
void SusyLLPTree::LoadTree(const char* file)
{
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->Get("SusyLLPTree"));
    InitTree();
    assert(tree_);
};

//CreateTree
void SusyLLPTree::CreateTree()
{
    //tree
    tree_ = new TTree("SusyLLPTree","SusyLLPTree");
    f_ = 0;

    //event info
    tree_->Branch("runNum",      &runNum,     "runNum/i");      // event run number
    tree_->Branch("lumiSec",     &lumiSec,    "lumiSec/i");     // event lumi section
    tree_->Branch("evtNum",      &evtNum,     "evtNum/i");      // event number
    tree_->Branch("category",    &category,   "category/i");    // dilepton category
    tree_->Branch("npv",         &npv,        "npv/i");         // number of primary vertices
    tree_->Branch("npu",         &npu,        "npu/i");         // number of in-time PU events (MC)
    tree_->Branch("pileupWeight",      &pileupWeight,     "pileupWeight/F");
    tree_->Branch("pileupWeightUp",      &pileupWeightUp,     "pileupWeightUp/F");
    tree_->Branch("pileupWeightDown",      &pileupWeightDown,     "pileupWeightDown/F");
    tree_->Branch("weight",      &weight,     "weight/F");
    tree_->Branch("rho",         &rho,        "rho/F");
    tree_->Branch("met",         &met,        "met/F");         // MET
    tree_->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)
    tree_->Branch("HT",    &HT,    "HT/F");
    tree_->Branch("jetMet_dPhi",    &jetMet_dPhi,    "jetMet_dPhi/F");
    tree_->Branch("jetMet_dPhiStar",    &jetMet_dPhiStar,    "jetMet_dPhiStar/F");
    tree_->Branch("jetMet_dPhiMin",    &jetMet_dPhiMin,    "jetMet_dPhiMin/F");
    tree_->Branch("jetMet_dPhiStarMin",    &jetMet_dPhiStarMin,    "jetMet_dPhiStarMin/F");
    tree_->Branch("jetMet_dPhiMin4",    &jetMet_dPhiMin4,    "jetMet_dPhiMin4/F");

    // met filters
    tree_->Branch("Flag2_globalSuperTightHalo2016Filter", &Flag2_globalSuperTightHalo2016Filter, "Flag2_globalSuperTightHalo2016Filter/O");
    tree_->Branch("Flag2_globalTightHalo2016Filter", &Flag2_globalTightHalo2016Filter, "Flag2_globalTightHalo2016Filter/O");
    tree_->Branch("Flag2_goodVertices", &Flag2_goodVertices, "Flag2_goodVertices/O");
    tree_->Branch("Flag2_BadChargedCandidateFilter", &Flag2_BadChargedCandidateFilter, "Flag2_BadChargedCandidateFilter/O");
    tree_->Branch("Flag2_BadPFMuonFilter", &Flag2_BadPFMuonFilter, "Flag2_BadPFMuonFilter/O");
    tree_->Branch("Flag2_EcalDeadCellTriggerPrimitiveFilter", &Flag2_EcalDeadCellTriggerPrimitiveFilter, "Flag2_EcalDeadCellTriggerPrimitiveFilter/O");
    tree_->Branch("Flag2_HBHENoiseFilter", &Flag2_HBHENoiseFilter, "Flag2_HBHENoiseFilter/O");
    tree_->Branch("Flag2_HBHEIsoNoiseFilter", &Flag2_HBHEIsoNoiseFilter, "Flag2_HBHEIsoNoiseFilter/O");
    tree_->Branch("Flag2_ecalBadCalibFilter", &Flag2_ecalBadCalibFilter, "Flag2_ecalBadCalibFilter/O");
    tree_->Branch("Flag2_eeBadScFilter", &Flag2_eeBadScFilter, "Flag2_eeBadScFilter/O");

    //leptons
    tree_->Branch("nLeptons",  &nLeptons, "nLeptons/I");
    tree_->Branch("lepE",      lepE,      "lepE[nLeptons]/F");
    tree_->Branch("lepPt",     lepPt,     "lepPt[nLeptons]/F");
    tree_->Branch("lepEta",    lepEta,    "lepEta[nLeptons]/F");
    tree_->Branch("lepPhi",    lepPhi,    "lepPhi[nLeptons]/F");
    tree_->Branch("lepPdgId",  lepPdgId,  "lepPdgId[nLeptons]/I");
    tree_->Branch("lepDZ",     lepDZ,     "lepDZ[nLeptons]/F");
    tree_->Branch("lepPassId", lepPassId, "lepPassId[nLeptons]/O");
    // tree_->Branch("lepLoosePassId", lepLoosePassId, "lepLoosePassId[nLeptons]/O");
    // tree_->Branch("lepMediumPassId", lepMediumPassId, "lepMediumPassId[nLeptons]/O");
    // tree_->Branch("lepTightPassId", lepTightPassId, "lepTightPassId[nLeptons]/O");

    //Z-candidate
    tree_->Branch("MT",      &MT,  "MT/F");
    tree_->Branch("ZMass",      &ZMass,  "ZMass/F");
    tree_->Branch("ZPt",        &ZPt,    "ZPt/F");
    tree_->Branch("ZEta",       &ZEta,   "ZEta/F");
    tree_->Branch("ZPhi",       &ZPhi,   "ZPhi/F");
    tree_->Branch("ZleptonIndex1", &ZleptonIndex1, "ZleptonIndex1/I");
    tree_->Branch("ZleptonIndex2", &ZleptonIndex2, "ZleptonIndex2/I");

    //hh-candidate
    tree_->Branch("dr_max",      &dr_max,  "dr_max/F");
    tree_->Branch("dm_hh",      &dm_hh,  "dm_hh/F");
    tree_->Branch("avg_mh",      &avg_mh,  "avg_mh/F");
    tree_->Branch("mh1",      &mh1,  "mh1/F");
    tree_->Branch("mh2",      &mh2,  "mh2/F");

/*    
  tree_->Branch("nElectrons", &nElectrons,"nElectrons/I");
  tree_->Branch("nPhotons", &nPhotons,"nPhotons/I");
  tree_->Branch("nTaus", &nTaus,"nTaus/I");
  tree_->Branch("nMuons", &nMuons,"nMuons/I");
*/

  //muons
  tree_->Branch("nMuons", &nMuons,"nMuons/I");
  tree_->Branch("muonE", muonE,"muonE[nMuons]/F");
  tree_->Branch("muonPt", muonPt,"muonPt[nMuons]/F");
  tree_->Branch("muonEta", muonEta,"muonEta[nMuons]/F");
  tree_->Branch("muonPhi", muonPhi,"muonPhi[nMuons]/F");
  tree_->Branch("muonCharge", muonCharge, "muonCharge[nMuons]/I");
  tree_->Branch("muonIsLoose", muonIsLoose,"muonIsLoose[nMuons]/O");
  tree_->Branch("muonIsMedium", muonIsMedium,"muonIsMedium[nMuons]/O");
  tree_->Branch("muonIsTight", muonIsTight,"muonIsTight[nMuons]/O");
  tree_->Branch("muon_pileupIso", muon_pileupIso, "muon_pileupIso[nMuons]/F");
  tree_->Branch("muon_chargedIso", muon_chargedIso, "muon_chargedIso[nMuons]/F");
  tree_->Branch("muon_photonIso", muon_photonIso, "muon_photonIso[nMuons]/F");
  tree_->Branch("muon_neutralHadIso", muon_neutralHadIso, "muon_neutralHadIso[nMuons]/F");

  //electrons
  tree_->Branch("nElectrons", &nElectrons,"nElectrons/I");
  tree_->Branch("eleE", eleE,"eleE[nElectrons]/F");
  tree_->Branch("elePt", elePt,"elePt[nElectrons]/F");
  tree_->Branch("eleEta", eleEta,"eleEta[nElectrons]/F");
  tree_->Branch("elePhi", elePhi,"elePhi[nElectrons]/F");
  tree_->Branch("eleCharge", eleCharge, "eleCharge[nElectrons]/F");
  tree_->Branch("ele_passCutBasedIDVeto", ele_passCutBasedIDVeto, "ele_passCutBasedIDVeto[nElectrons]/O");
  tree_->Branch("ele_passCutBasedIDLoose", ele_passCutBasedIDLoose, "ele_passCutBasedIDLoose[nElectrons]/O");
  tree_->Branch("ele_passCutBasedIDMedium", ele_passCutBasedIDMedium, "ele_passCutBasedIDMedium[nElectrons]/O");
  tree_->Branch("ele_passCutBasedIDTight", ele_passCutBasedIDTight, "ele_passCutBasedIDTight[nElectrons]/O");
  //tree_->Branch("ele_passMVAIsoIDWP80", ele_passMVAIsoIDWP80, "ele_passMVAIsoIDWP80[nElectrons]/O");
  //tree_->Branch("ele_passMVAIsoIDWP90", ele_passMVAIsoIDWP90, "ele_passMVAIsoIDWP90[nElectrons]/O"); 

  //taus
  tree_->Branch("nTaus", &nTaus,"nTaus/I");
  tree_->Branch("tauE", tauE,"tauE[nTaus]/F");
  tree_->Branch("tauPt", tauPt,"tauPt[nTaus]/F");
  tree_->Branch("tauEta", tauEta,"tauEta[nTaus]/F");
  tree_->Branch("tauPhi", tauPhi,"tauPhi[nTaus]/F");
  tree_->Branch("tau_IsLoose", tau_IsLoose, "tau_IsLoose[nTaus]/O");
  tree_->Branch("tau_IsMedium", tau_IsMedium, "tau_IsMedium[nTaus]/O");
  tree_->Branch("tau_IsTight", tau_IsTight, "tau_IsTight[nTaus]/O");
  tree_->Branch("tau_passEleVetoLoose", tau_passEleVetoLoose, "tau_passEleVetoLoose[nTaus]/O");
  tree_->Branch("tau_passEleVetoMedium", tau_passEleVetoMedium, "tau_passEleVetoMedium[nTaus]/O");
  tree_->Branch("tau_passEleVetoTight", tau_passEleVetoTight, "tau_passEleVetoTight[nTaus]/O");
  tree_->Branch("tau_passMuVetoLoose", tau_passMuVetoLoose, "tau_passMuVetoLoose[nTaus]/O");
  tree_->Branch("tau_passMuVetoMedium", tau_passMuVetoMedium, "tau_passMuVetoMedium[nTaus]/O");
  tree_->Branch("tau_passMuVetoTight", tau_passMuVetoTight, "tau_passMuVetoTight[nTaus]/O");

  //photons
  tree_->Branch("nPhotons", &nPhotons,"nPhotons/I");
  //tree_->Branch("nPhotons_overlap", &nPhotons_overlap,"nPhotons_overlap/I");
  tree_->Branch("phoE", phoE,"phoE[nPhotons]/F");
  tree_->Branch("phoPt", phoPt,"phoPt[nPhotons]/F");
  tree_->Branch("phoEta", phoEta,"phoEta[nPhotons]/F");
  tree_->Branch("phoPhi", phoPhi,"phoPhi[nPhotons]/F");
  tree_->Branch("pho_passCutBasedIDLoose", pho_passCutBasedIDLoose, "pho_passCutBasedIDLoose[nPhotons]/O");
  tree_->Branch("pho_passCutBasedIDMedium", pho_passCutBasedIDMedium, "pho_passCutBasedIDMedium[nPhotons]/O");
  tree_->Branch("pho_passCutBasedIDTight", pho_passCutBasedIDTight, "pho_passCutBasedIDTight[nPhotons]/O");
  tree_->Branch("pho_passMVAIDWP80", pho_passMVAIDWP80, "pho_passMVAIDWP80[nPhotons]/O");
  tree_->Branch("pho_passMVAIDWP90", pho_passMVAIDWP90, "pho_passMVAIDWP90[nPhotons]/O");

    //jets
  tree_->Branch("nJets", &nJets,"nJets/I");
  tree_->Branch("jetE", jetE,"jetE[nJets]/F");
  tree_->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  tree_->Branch("jetEta", jetEta,"jetEta[nJets]/F");

  tree_->Branch("jetPhi", jetPhi,"jetPhi[nJets]/F");
  tree_->Branch("jetCSV", jetCSV,"jetCSV[nJets]/F");
  tree_->Branch("jetCISV", jetCISV,"jetCISV[nJets]/F");
  tree_->Branch("jetProbb", jetProbb,"jetProbb[nJets]/F");
  tree_->Branch("jetProbc", jetProbc,"jetProbc[nJets]/F");
  tree_->Branch("jetProbudsg", jetProbudsg,"jetProbudsg[nJets]/F");
  tree_->Branch("jetProbbb", jetProbbb,"jetProbbb[nJets]/F");
  tree_->Branch("jetMass", jetMass, "jetMass[nJets]/F");
  tree_->Branch("jetJetArea", jetJetArea, "jetJetArea[nJets]/F");
  tree_->Branch("jetPileupE", jetPileupE, "jetPileupE[nJets]/F");
  tree_->Branch("jetPileupId", jetPileupId, "jetPileupId[nJets]/F");
  tree_->Branch("jetPileupIdFlag", jetPileupIdFlag, "jetPileupIdFlag[nJets]/I");
  tree_->Branch("jetPassIDLoose", jetPassIDLoose, "jetPassIDLoose[nJets]/O");
  tree_->Branch("jetPassIDTight", jetPassIDTight, "jetPassIDTight[nJets]/O");
  tree_->Branch("jetPassMuFrac", jetPassMuFrac, "jetPassMuFrac[nJets]/O");
  tree_->Branch("jetPassEleFrac", jetPassEleFrac, "jetPassEleFrac[nJets]/O");
  tree_->Branch("jetPartonFlavor", jetPartonFlavor, "jetPartonFlavor[nJets]/I");
  tree_->Branch("jetHadronFlavor", jetHadronFlavor, "jetHadronFlavor[nJets]/I");
  tree_->Branch("jetElectronEnergyFraction", jetElectronEnergyFraction, "jetElectronEnergyFraction[nJets]/F");
  tree_->Branch("jetPhotonEnergyFraction", jetPhotonEnergyFraction, "jetPhotonEnergyFraction[nJets]/F");
  tree_->Branch("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction, "jetChargedHadronEnergyFraction[nJets]/F");
  tree_->Branch("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction, "jetNeutralHadronEnergyFraction[nJets]/F");
  tree_->Branch("jetMuonEnergyFraction", jetMuonEnergyFraction, "jetMuonEnergyFraction[nJets]/F");
  tree_->Branch("jetHOEnergyFraction", jetHOEnergyFraction, "jetHOEnergyFraction[nJets]/F");
  tree_->Branch("jetHFHadronEnergyFraction", jetHFHadronEnergyFraction, "jetHFHadronEnergyFraction[nJets]/F");
  tree_->Branch("jetHFEMEnergyFraction",jetHFEMEnergyFraction, "jetHFEMEnergyFraction[nJets]/F");
  tree_->Branch("jetChargedHadronMultiplicity", jetChargedHadronMultiplicity, "jetChargedHadronMultiplicity[nJets]/I");
  tree_->Branch("jetNeutralHadronMultiplicity", jetNeutralHadronMultiplicity, "jetNeutralHadronMultiplicity[nJets]/I");
  tree_->Branch("jetPhotonMultiplicity", jetPhotonMultiplicity, "jetPhotonMultiplicity[nJets]/I");
  tree_->Branch("jetElectronMultiplicity", jetElectronMultiplicity, "jetElectronMultiplicity[nJets]/I");
  tree_->Branch("jetMuonMultiplicity", jetMuonMultiplicity, "jetMuonMultiplicity[nJets]/I");
  tree_->Branch("jetAllMuonPt", jetAllMuonPt,"jetAllMuonPt[nJets]/F");
  tree_->Branch("jetAllMuonEta", jetAllMuonEta,"jetAllMuonEta[nJets]/F");
  tree_->Branch("jetAllMuonPhi", jetAllMuonPhi,"jetAllMuonPhi[nJets]/F");
  tree_->Branch("jetAllMuonM", jetAllMuonM,"jetAllMuonM[nJets]/F");
  tree_->Branch("jetPtWeightedDZ", jetPtWeightedDZ,"jetPtWeightedDZ[nJets]/F");
  tree_->Branch("jetNRechits", jetNRechits,"jetNRechits[nJets]/I");
  tree_->Branch("jetRechitE", jetRechitE,"jetRechitE[nJets]/F");
  tree_->Branch("jetRechitT", jetRechitT,"jetRechitT[nJets]/F");
  tree_->Branch("jetRechitT_rms", jetRechitT_rms,"jetRechitT_rms[nJets]/F");

  tree_->Branch("jetRechitE_Error", jetRechitE_Error,"jetRechitE_Error[nJets]/F");
  tree_->Branch("jetRechitT_Error", jetRechitT_Error,"jetRechitT_Error[nJets]/F");
  tree_->Branch("jetAlphaMax",jetAlphaMax,"jetAlphaMax[nJets]/F");
  tree_->Branch("jetBetaMax",jetBetaMax,"jetBetaMax[nJets]/F");
  tree_->Branch("jetGammaMax_ET",jetGammaMax_ET,"jetGammaMax_ET[nJets]/F");
  tree_->Branch("jetGammaMax_EM",jetGammaMax_EM,"jetGammaMax_EM[nJets]/F");
  tree_->Branch("jetGammaMax_Hadronic",jetGammaMax_Hadronic,"jetGammaMax_Hadronic[nJets]/F");
  tree_->Branch("jetGammaMax",jetGammaMax,"jetGammaMax[nJets]/F");
  tree_->Branch("jetPtAllTracks",jetPtAllTracks,"jetPtAllTracks[nJets]/F");
  tree_->Branch("jetPtAllPVTracks",jetPtAllPVTracks,"jetPtAllPVTracks[nJets]/F");
  tree_->Branch("jetMedianTheta2D",jetMedianTheta2D,"jetMedianTheta2D[nJets]/F");
  tree_->Branch("jetMedianIP",jetMedianIP,"jetMedianIP[nJets]/F");
  tree_->Branch("jetMinDeltaRAllTracks",jetMinDeltaRAllTracks,"jetMinDeltaRAllTracks[nJets]/F");
  tree_->Branch("jetMinDeltaRPVTracks",jetMinDeltaRPVTracks,"jetMinDeltaRPVTracks[nJets]/F");
  tree_->Branch("jet_sig_et1",jet_sig_et1,"jet_sig_et1[nJets]/F");
  tree_->Branch("jet_sig_et2",jet_sig_et2,"jet_sig_et2[nJets]/F");
  tree_->Branch("jet_energy_frac",jet_energy_frac,"jet_energy_frac[nJets]/F");

  tree_->Branch("jetAlphaMax_wp",jetAlphaMax_wp,"jetAlphaMax_wp[nJets]/F");
  tree_->Branch("jetBetaMax_wp",jetBetaMax_wp,"jetBetaMax_wp[nJets]/F");
  tree_->Branch("jetGammaMax_ET_wp",jetGammaMax_ET_wp,"jetGammaMax_ET_wp[nJets]/F");
  tree_->Branch("jetGammaMax_EM_wp",jetGammaMax_EM_wp,"jetGammaMax_EM_wp[nJets]/F");
  tree_->Branch("jetGammaMax_Hadronic_wp",jetGammaMax_Hadronic_wp,"jetGammaMax_Hadronic_wp[nJets]/F");
  tree_->Branch("jetGammaMax_wp",jetGammaMax_wp,"jetGammaMax_wp[nJets]/F");
  tree_->Branch("jetPtAllTracks_wp",jetPtAllTracks_wp,"jetPtAllTracks_wp[nJets]/F");
  tree_->Branch("jetPtAllPVTracks_wp",jetPtAllPVTracks_wp,"jetPtAllPVTracks_wp[nJets]/F");
  tree_->Branch("jetMedianTheta2D_wp",jetMedianTheta2D_wp,"jetMedianTheta2D_wp[nJets]/F");
  tree_->Branch("jetMedianIP_wp",jetMedianIP_wp,"jetMedianIP_wp[nJets]/F");
  tree_->Branch("jetMinDeltaRAllTracks_wp",jetMinDeltaRAllTracks_wp,"jetMinDeltaRAllTracks_wp[nJets]/F");
  tree_->Branch("jetMinDeltaRPVTracks_wp",jetMinDeltaRPVTracks_wp,"jetMinDeltaRPVTracks_wp[nJets]/F");

    tree_->Branch("jetTime",   jetTime,   "jetTime[nJets]/F");
    tree_->Branch("ecalNRechits",   ecalNRechits,   "ecalNRechits[nJets]/F");
    tree_->Branch("ecalRechitE", ecalRechitE, "ecalRechitE[nJets]/F");

    tree_->Branch("jetChargedEMEnergyFraction",   jetChargedEMEnergyFraction,   "jetChargedEMEnergyFraction[nJets]/F");
    tree_->Branch("jetNeutralEMEnergyFraction",   jetNeutralEMEnergyFraction,   "jetNeutralEMEnergyFraction[nJets]/F");
/*
    tree_->Branch("nJets",     &nJets,    "nJets/I");
    tree_->Branch("nBJets",     &nBJets,    "nBJets/I");
    tree_->Branch("jetE",      jetE,      "jetE[nJets]/F");
    tree_->Branch("jetEt",      jetEt,      "jetEt[nJets]/F");
    tree_->Branch("jetPt",     jetPt,     "jetPt[nJets]/F");
    tree_->Branch("jetEta",    jetEta,    "jetEta[nJets]/F");
    tree_->Branch("jetPhi",    jetPhi,    "jetPhi[nJets]/F");
    tree_->Branch("jetTime",   jetTime,   "jetTime[nJets]/F");
    tree_->Branch("jetPassId", jetPassId, "jetPassId[nJets]/O");
    tree_->Branch("jetCSVT", jetCSVT, "jetCSVT[nJets]/O");
    tree_->Branch("jetCISV", jetCISV, "jetCISV[nJets]/F");
    tree_->Branch("matched", matched, "matched[nJets]/O");
    tree_->Branch("jet_matched_gLLP0_daughter", jet_matched_gLLP0_daughter, "jet_matched_gLLP0_daughter[nJets]/O");
    tree_->Branch("jet_matched_gLLP1_daughter", jet_matched_gLLP1_daughter, "jet_matched_gLLP1_daughter[nJets]/O");
    tree_->Branch("jet_matched_gLLP0_grandaughter", jet_matched_gLLP0_grandaughter, "jet_matched_gLLP0_grandaughter[nJets]/O");
    tree_->Branch("jet_matched_gLLP1_grandaughter", jet_matched_gLLP1_grandaughter, "jet_matched_gLLP1_grandaughter[nJets]/O");
    tree_->Branch("jet_energy_frac", jet_energy_frac, "jet_energy_frac[nJets]/F");
    tree_->Branch("jet_sig_et1", jet_sig_et1, "jet_sig_et1[nJets]/F");
    tree_->Branch("jet_sig_et2", jet_sig_et2, "jet_sig_et2[nJets]/F");


    tree_->Branch("jetGammaMax_ET", jetGammaMax_ET, "jetGammaMax_ET[nJets]/F");
    tree_->Branch("jetMinDeltaRPVTracks", jetMinDeltaRPVTracks, "jetMinDeltaRPVTracks[nJets]/F");
    tree_->Branch("jetPtAllPVTracks", jetPtAllPVTracks, "jetPtAllPVTracks[nJets]/F");
    tree_->Branch("jetMinDeltaRAllTracks", jetMinDeltaRAllTracks, "jetMinDeltaRAllTracks[nJets]/F");
    tree_->Branch("jetPtAllTracks", jetPtAllTracks, "jetPtAllTracks[nJets]/F");

    tree_->Branch("jetGammaMax_ET_wp", jetGammaMax_ET_wp, "jetGammaMax_ET_wp[nJets]/F");
    tree_->Branch("jetMinDeltaRPVTracks_wp", jetMinDeltaRPVTracks_wp, "jetMinDeltaRPVTracks_wp[nJets]/F");
    tree_->Branch("jetPtAllPVTracks_wp", jetPtAllPVTracks_wp, "jetPtAllPVTracks_wp[nJets]/F");
    tree_->Branch("jetMinDeltaRAllTracks_wp", jetMinDeltaRAllTracks_wp, "jetMinDeltaRAllTracks_wp[nJets]/F");
    tree_->Branch("jetPtAllTracks_wp", jetPtAllTracks_wp, "jetPtAllTracks_wp[nJets]/F");

    tree_->Branch("ecalNRechits",   ecalNRechits,   "ecalNRechits[nJets]/F");
    tree_->Branch("ecalRechitE", ecalRechitE, "ecalRechitE[nJets]/F");
    tree_->Branch("jetChargedEMEnergyFraction",   jetChargedEMEnergyFraction,   "jetChargedEMEnergyFraction[nJets]/F");
    tree_->Branch("jetNeutralEMEnergyFraction",   jetNeutralEMEnergyFraction,   "jetNeutralEMEnergyFraction[nJets]/F");
    tree_->Branch("jetChargedHadronEnergyFraction",   jetChargedHadronEnergyFraction,   "jetChargedHadronEnergyFraction[nJets]/F");
    tree_->Branch("jetNeutralHadronEnergyFraction",   jetNeutralHadronEnergyFraction,   "jetNeutralHadronEnergyFraction[nJets]/F");
    tree_->Branch("jetHoverE",   &jetHoverE,   "jetHoverE[nJets]/F");
    // tree_->Branch("jetLoosePassId", jetLoosePassId, "jetLoosePassId[nJets]/O");
    // tree_->Branch("jetTightPassId", jetTightPassId, "jetTightPassId[nJets]/O");
*/

/*
    //calojet
    tree_->Branch("nCaloJets",     &nCaloJets,    "nCaloJets/I");
    tree_->Branch("calojetE",      calojetE,      "calojetE[nCaloJets]/F");
    tree_->Branch("calojetEt",      calojetEt,      "calojetEt[nCaloJets]/F");
    tree_->Branch("calojetPt",     calojetPt,     "calojetPt[nCaloJets]/F");
    tree_->Branch("calojetEta",    calojetEta,    "calojetEta[nCaloJets]/F");
    tree_->Branch("calojetPhi",    calojetPhi,    "calojetPhi[nCaloJets]/F");
    tree_->Branch("calojetTime",   calojetTime,   "calojetTime[nCaloJets]/F");
    tree_->Branch("calojetPassId", calojetPassId, "calojetPassId[nCaloJets]/O");
    //tree_->Branch("calojetChargedEMEnergyFraction",   calojetChargedEMEnergyFraction,   "calojetChargedEMEnergyFraction[nCaloJets]/F");
    //tree_->Branch("calojetNeutralEMEnergyFraction",   calojetNeutralEMEnergyFraction,   "calojetNeutralEMEnergyFraction[nCaloJets]/F");
    //tree_->Branch("calojetChargedHadronEnergyFraction",   calojetChargedHadronEnergyFraction,   "calojetChargedHadronEnergyFraction[nCaloJets]/F");
    //tree_->Branch("calojetNeutralHadronEnergyFraction",   calojetNeutralHadronEnergyFraction,   "calojetNeutralHadronEnergyFraction[nCaloJets]/F");
    tree_->Branch("calojet_EMEnergyFraction", calojet_EMEnergyFraction, "calojet_EMEnergyFraction[nCaloJets]/F");
    tree_->Branch("calojet_HadronicEnergyFraction",   calojet_HadronicEnergyFraction,   "calojet_HadronicEnergyFraction[nCaloJets]/F");
    tree_->Branch("calojetGammaMax_ET", calojetGammaMax_ET, "calojetGammaMax_ET[nCaloJets]/F");
    tree_->Branch("calojetMinDeltaRPVTracks", calojetMinDeltaRPVTracks, "calojetMinDeltaRPVTracks[nCaloJets]/F");
    tree_->Branch("calojetPtAllPVTracks", calojetPtAllPVTracks, "calojetPtAllPVTracks[nCaloJets]/F");
    tree_->Branch("calojetNRechits",   calojetNRechits,   "calojetNRechits[nCaloJets]/F");
    tree_->Branch("calojetRechitE", calojetRechitE, "calojetRechitE[nCaloJets]/F");
*/
    //HLT
    tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[602]/O"); //hardcoded

  //gLLP branches
  tree_->Branch("gLLP0_EB", &gLLP0_EB, "gLLP0_EB/O");
  tree_->Branch("gLLP1_EB", &gLLP1_EB, "gLLP1_EB/O");
  tree_->Branch("gLLP0_EE", &gLLP0_EE, "gLLP0_EE/O");
  tree_->Branch("gLLP1_EE", &gLLP1_EE, "gLLP1_EE/O");

  tree_->Branch("gLLP_travel_time",          gLLP_travel_time,          "gLLP_travel_time[2]/F");
  tree_->Branch("gLLP_e",          gLLP_e,          "gLLP_e[2]/F");
  tree_->Branch("gLLP_pt",          gLLP_pt,          "gLLP_pt[2]/F");
  tree_->Branch("gLLP_eta",          gLLP_eta,          "gLLP_eta[2]/F");
  tree_->Branch("gLLP_beta",          gLLP_beta,          "gLLP_beta[2]/F");
  tree_->Branch("gLLP_phi",          gLLP_phi,          "gLLP_phi[2]/F");
  tree_->Branch("gLLP_ctau",          gLLP_ctau,          "gLLP_ctau[2]/F");
  tree_->Branch("gLLP_dr",          &gLLP_dr,          "gLLP_dr/F");
  tree_->Branch("gLLP_grandaughter_dr1",          &gLLP_grandaughter_dr1,          "gLLP_grandaughter_dr1/F");
  tree_->Branch("gLLP_grandaughter_dr2",          &gLLP_grandaughter_dr2,          "gLLP_grandaughter_dr2/F");
  tree_->Branch("gLLP_grandaughter_matched_jet_dr1",          &gLLP_grandaughter_matched_jet_dr1,          "gLLP_grandaughter_matched_jet_dr1/F");
  tree_->Branch("gLLP_grandaughter_matched_jet_dr2",          &gLLP_grandaughter_matched_jet_dr2,          "gLLP_grandaughter_matched_jet_dr2/F");

  tree_->Branch("gLLP_decay_vertex_r",          gLLP_decay_vertex_r,          "gLLP_decay_vertex_r[2]/F");
  tree_->Branch("gLLP_decay_vertex_x",          gLLP_decay_vertex_x,          "gLLP_decay_vertex_x[2]/F");
  tree_->Branch("gLLP_decay_vertex_y",          gLLP_decay_vertex_y,          "gLLP_decay_vertex_y[2]/F");
  tree_->Branch("gLLP_decay_vertex_z",          gLLP_decay_vertex_z,          "gLLP_decay_vertex_z[2]/F");

  tree_->Branch("gLLP_prod_vertex_x",          gLLP_prod_vertex_x,          "gLLP_prod_vertex_x[2]/F");
  tree_->Branch("gLLP_prod_vertex_y",          gLLP_prod_vertex_y,          "gLLP_prod_vertex_y[2]/F");
  tree_->Branch("gLLP_prod_vertex_z",          gLLP_prod_vertex_z,          "gLLP_prod_vertex_z[2]/F");
/*
  //gLLP daughter branches
  tree_->Branch("gen_time",          gen_time,          "gen_time[4]/F");
  tree_->Branch("photon_travel_time",          photon_travel_time,          "photon_travel_time[4]/F");
  tree_->Branch("gLLP_daughter_travel_time",          gLLP_daughter_travel_time,          "gLLP_daughter_travel_time[4]/F");
  tree_->Branch("gLLP_daughter_e",          gLLP_daughter_e,          "gLLP_daughter_e[4]/F");
  tree_->Branch("gLLP_daughter_pt",          gLLP_daughter_pt,          "gLLP_daughter_pt[4]/F");
  tree_->Branch("gLLP_daughter_eta",          gLLP_daughter_eta,          "gLLP_daughter_eta[4]/F");
  tree_->Branch("gLLP_daughter_phi",          gLLP_daughter_phi,          "gLLP_daughter_phi[4]/F");
  tree_->Branch("gLLP_daughter_eta_ecalcorr",          gLLP_daughter_eta_ecalcorr,          "gLLP_daughter_eta_ecalcorr[4]/F");
  tree_->Branch("gLLP_daughter_phi_ecalcorr",          gLLP_daughter_phi_ecalcorr,          "gLLP_daughter_phi_ecalcorr[4]/F");
  tree_->Branch("gLLP_min_delta_r_match_jet",          gLLP_min_delta_r_match_jet,          "gLLP_min_delta_r_match_jet[4]/F");
  tree_->Branch("gLLP_daughter_match_jet_index",          gLLP_daughter_match_jet_index,          "gLLP_daughter_match_jet_index[4]/i");
*/
   //daughters
   tree_->Branch("gLLP_daughter_EB", gLLP_daughter_EB, "gLLP_daughter_EB[4]/O");
   tree_->Branch("gLLP_daughter_ETL", gLLP_daughter_ETL, "gLLP_daughter_ETL[4]/O");

   tree_->Branch("gLLP_daughter_photon_travel_time_EB", gLLP_daughter_photon_travel_time_EB, "gLLP_daughter_photon_travel_time_EB[4]/F");
   tree_->Branch("gLLP_daughter_photon_travel_time_ETL", gLLP_daughter_photon_travel_time_ETL, "gLLP_daughter_photon_travel_time_ETL[4]/F");

   tree_->Branch("gLLP_daughter_travel_time_EB", gLLP_daughter_travel_time_EB, "gLLP_daughter_travel_time_EB[4]/F");
   tree_->Branch("gLLP_daughter_travel_time_ETL", gLLP_daughter_travel_time_ETL, "gLLP_daughter_travel_time_ETL[4]/F");

   tree_->Branch("gen_time_daughter_EB", gen_time_daughter_EB, "gen_time_daughter_EB[4]/F");
   tree_->Branch("gen_time_daughter_ETL", gen_time_daughter_ETL, "gen_time_daughter_ETL[4]/F");

   tree_->Branch("gLLP_daughter_id", gLLP_daughter_id, "gLLP_daughter_id[4]/I");
   tree_->Branch("gLLP_daughter_pt", gLLP_daughter_pt, "gLLP_daughter_pt[4]/F");
   tree_->Branch("gLLP_daughter_eta", gLLP_daughter_eta, "gLLP_daughter_eta[4]/F");
   tree_->Branch("gLLP_daughter_phi", gLLP_daughter_phi, "gLLP_daughter_phi[4]/F");
   tree_->Branch("gLLP_daughter_eta_ecalcorr", gLLP_daughter_eta_ecalcorr, "gLLP_daughter_eta_ecalcorr[4]/F");
   tree_->Branch("gLLP_daughter_phi_ecalcorr", gLLP_daughter_phi_ecalcorr, "gLLP_daughter_phi_ecalcorr[4]/F");
   tree_->Branch("gLLP_daughter_e", gLLP_daughter_e, "gLLP_daughter_e[4]/F");
   tree_->Branch("gLLP_daughter_mass", gLLP_daughter_mass, "gLLP_daughter_mass[4]/F");

   tree_->Branch("gLLP_daughter_min_delta_r_match_calojet", gLLP_daughter_min_delta_r_match_calojet, "gLLP_daughter_min_delta_r_match_calojet[4]/F");
   tree_->Branch("gLLP_daughter_match_calojet_index", gLLP_daughter_match_calojet_index, "gLLP_daughter_match_calojet_index[4]/i");
   tree_->Branch("gLLP_daughter_match_jet_index", gLLP_daughter_match_jet_index, "gLLP_daughter_match_jet_index[4]/i");
   tree_->Branch("gLLP_daughter_min_delta_r_match_jet", gLLP_daughter_min_delta_r_match_jet, "gLLP_daughter_min_delta_r_match_jet[4]/F");

   //grandaughters
   tree_->Branch("gLLP_grandaughter_EB", gLLP_grandaughter_EB, "gLLP_grandaughter_EB[4]/O");
   tree_->Branch("gLLP_grandaughter_ETL", gLLP_grandaughter_ETL, "gLLP_grandaughter_ETL[4]/O");

   tree_->Branch("gLLP_grandaughter_photon_travel_time_EB", gLLP_grandaughter_photon_travel_time_EB, "gLLP_grandaughter_photon_travel_time_EB[4]/F");
   tree_->Branch("gLLP_grandaughter_photon_travel_time_ETL", gLLP_grandaughter_photon_travel_time_ETL, "gLLP_grandaughter_photon_travel_time_ETL[4]/F");

   tree_->Branch("gLLP_grandaughter_travel_time_EB", gLLP_grandaughter_travel_time_EB, "gLLP_grandaughter_travel_time_EB[4]/F");
   tree_->Branch("gLLP_grandaughter_travel_time_ETL", gLLP_grandaughter_travel_time_ETL, "gLLP_grandaughter_travel_time_ETL[4]/F");

   tree_->Branch("gen_time_grandaughter_EB", gen_time_grandaughter_EB, "gen_time_grandaughter_EB[4]/F");
   tree_->Branch("gen_time_grandaughter_ETL", gen_time_grandaughter_ETL, "gen_time_grandaughter_ETL[4]/F");

   tree_->Branch("gLLP_grandaughter_id", gLLP_grandaughter_id, "gLLP_grandaughter_id[4]/I");
   tree_->Branch("gLLP_grandaughter_pt", gLLP_grandaughter_pt, "gLLP_grandaughter_pt[4]/F");
   tree_->Branch("gLLP_grandaughter_eta", gLLP_grandaughter_eta, "gLLP_grandaughter_eta[4]/F");
   tree_->Branch("gLLP_grandaughter_phi", gLLP_grandaughter_phi, "gLLP_grandaughter_phi[4]/F");
   tree_->Branch("gLLP_grandaughter_eta_ecalcorr", gLLP_grandaughter_eta_ecalcorr, "gLLP_grandaughter_eta_ecalcorr[4]/F");
   tree_->Branch("gLLP_grandaughter_phi_ecalcorr", gLLP_grandaughter_phi_ecalcorr, "gLLP_grandaughter_phi_ecalcorr[4]/F");
   tree_->Branch("gLLP_grandaughter_e", gLLP_grandaughter_e, "gLLP_grandaughter_e[4]/F");
   tree_->Branch("gLLP_grandaughter_mass", gLLP_grandaughter_mass, "gLLP_grandaughter_mass[4]/F");

   //tree_->Branch("gLLP_grandaughter_min_delta_r_match_calojet", gLLP_grandaughter_min_delta_r_match_calojet, "gLLP_grandaughter_min_delta_r_match_calojet[4]/F");
   //tree_->Branch("gLLP_grandaughter_match_calojet_index", gLLP_grandaughter_match_calojet_index, "gLLP_grandaughter_match_calojet_index[4]/i");
   tree_->Branch("gLLP_grandaughter_match_jet_index", gLLP_grandaughter_match_jet_index, "gLLP_grandaughter_match_jet_index[4]/i");
   tree_->Branch("gLLP_grandaughter_min_delta_r_match_jet", gLLP_grandaughter_min_delta_r_match_jet, "gLLP_grandaughter_min_delta_r_match_jet[4]/F");

};
