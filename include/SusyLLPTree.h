// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef SusyLLPTree_H
#define SusyLLPTree_H

#define N_MAX_LLPS 2
#define N_MAX_LLP_DAUGHTERS 4
#define N_MAX_LLP_GRAND_DAUGHTERS 4
#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 2000
#define NTriggersMAX 601 //Number of trigger in the .dat file
#define N_CSC_CUT 20
#define JET_PT_CUT 10
#define MUON_PT_CUT 20

#include <iostream>
#include <string>
#include <sys/stat.h>
#include "assert.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TRandom.h"
#include "FactorizedJetCorrector.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "SimpleJetResolution.h"
#include "BTagCalibrationStandalone.h"
#include "TTree.h"

#include "RazorAnalyzer.h"

#include "RazorHelper.h"

class SusyLLPTree
{

public:
  SusyLLPTree();
  ~SusyLLPTree();

  //tree
  TTree *tree_;
  TFile *f_;

  //event info
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t rho, weight;
  Float_t met, metPhi;

  //leptons
  int nLeptons;
  float lepE[N_MAX_LEPTONS];
  float lepPt[N_MAX_LEPTONS];
  float lepEta[N_MAX_LEPTONS];
  float lepPhi[N_MAX_LEPTONS];
  int  lepPdgId[N_MAX_LEPTONS];
  float lepDZ[N_MAX_LEPTONS];
  // bool lepLoosePassId[N_MAX_LEPTONS];
  // bool lepMediumPassId[N_MAX_LEPTONS];
  // bool lepTightPassId[N_MAX_LEPTONS];
  bool lepPassId[N_MAX_LEPTONS];

  //Z-candidate
  float MT;
  float ZMass;
  float ZPt;
  float ZEta;
  float ZPhi;
  int ZleptonIndex1;
  int ZleptonIndex2;

  //jets
  int nJets;
  float jetE[N_MAX_JETS];
  float jetEt[N_MAX_JETS];
  float jetPt[N_MAX_JETS];
  float jetEta[N_MAX_JETS];
  float jetPhi[N_MAX_JETS];
  float jetTime[N_MAX_JETS];
  float ecalNRechits[N_MAX_JETS];
  float ecalRechitE[N_MAX_JETS];
  float jetChargedEMEnergyFraction[N_MAX_JETS];
  float jetNeutralEMEnergyFraction[N_MAX_JETS];
  float jetChargedHadronEnergyFraction[N_MAX_JETS];
  float jetNeutralHadronEnergyFraction[N_MAX_JETS];
  float jetGammaMax_ET[N_MAX_JETS];
  float jetMinDeltaRPVTracks[N_MAX_JETS];
  float jetPtAllPVTracks[N_MAX_JETS];
  float jetMinDeltaRAllTracks[N_MAX_JETS];
  float jetPtAllTracks[N_MAX_JETS];
  float jetGammaMax_ET_wp[N_MAX_JETS];
  float jetMinDeltaRPVTracks_wp[N_MAX_JETS];
  float jetPtAllPVTracks_wp[N_MAX_JETS];
  float jetMinDeltaRAllTracks_wp[N_MAX_JETS];
  float jetPtAllTracks_wp[N_MAX_JETS];
  // bool jetLoosePassId[N_MAX_JETS];
  bool jetPassId[N_MAX_JETS];
  bool matched[N_MAX_JETS];
  bool jet_matched_gLLP0_daughter[N_MAX_JETS];
  bool jet_matched_gLLP1_daughter[N_MAX_JETS];
  bool jet_matched_gLLP0_grandaughter[N_MAX_JETS];
  bool jet_matched_gLLP1_grandaughter[N_MAX_JETS];
  // bool jetTightPassId[N_MAX_JETS];
  float jet_energy_frac[N_MAX_JETS];
  float jet_sig_et1[N_MAX_JETS];
  float jet_sig_et2[N_MAX_JETS];

  //calojet
  int nCaloJets;
  float calojetE[N_MAX_JETS];
  float calojetEt[N_MAX_JETS];
  float calojetPt[N_MAX_JETS];
  float calojetEta[N_MAX_JETS];
  float calojetPhi[N_MAX_JETS];
  float calojetTime[N_MAX_JETS];
  float calojetNRechits[N_MAX_JETS];
  float calojetRechitE[N_MAX_JETS];
  //float calojetChargedEMEnergyFraction[N_MAX_JETS];
  //float calojetNeutralEMEnergyFraction[N_MAX_JETS];
  //float calojetChargedHadronEnergyFraction[N_MAX_JETS];
  //float calojetNeutralHadronEnergyFraction[N_MAX_JETS];
  float calojet_EMEnergyFraction[N_MAX_JETS];
  float calojet_HadronicEnergyFraction[N_MAX_JETS];
  float calojetGammaMax_ET[N_MAX_JETS];
  float calojetMinDeltaRPVTracks[N_MAX_JETS];
  float calojetPtAllPVTracks[N_MAX_JETS];
  // bool jetLoosePassId[N_MAX_JETS];
  bool calojetPassId[N_MAX_JETS];
  // // bool jetTightPassId[N_MAX_JETS];

  //gLLP
  float gLLP_csc[N_MAX_LLPS];
  float gLLP_ctau[N_MAX_LLPS];
 
   //gLLP
   bool          gLLP0_EB;
   bool          gLLP1_EB;
   bool          gLLP0_EE;
   bool          gLLP1_EE;

   float         gLLP_travel_time[N_MAX_LLPS];
   float         gLLP_e[N_MAX_LLPS];
   float         gLLP_pt[N_MAX_LLPS];
   float         gLLP_eta[N_MAX_LLPS];
   float         gLLP_beta[N_MAX_LLPS];
   float         gLLP_phi[N_MAX_LLPS];
   float         gLLP_decay_vertex_x[N_MAX_LLPS];
   float         gLLP_decay_vertex_y[N_MAX_LLPS];
   float         gLLP_decay_vertex_z[N_MAX_LLPS];
   float         gLLP_prod_vertex_x[N_MAX_LLPS];
   float         gLLP_prod_vertex_y[N_MAX_LLPS];
   float         gLLP_prod_vertex_z[N_MAX_LLPS];
/*
   //gLLP_daughter
   float         gen_time[N_MAX_LLP_DAUGHTERS];
   float         photon_travel_time[N_MAX_LLP_DAUGHTERS];
   float         gLLP_daughter_travel_time[N_MAX_LLP_DAUGHTERS];
   float         gLLP_daughter_e[N_MAX_LLP_DAUGHTERS];
   float         gLLP_daughter_pt[N_MAX_LLP_DAUGHTERS];
   float         gLLP_daughter_eta[N_MAX_LLP_DAUGHTERS];
   float         gLLP_daughter_phi[N_MAX_LLP_DAUGHTERS];
   float         gLLP_daughter_eta_ecalcorr[N_MAX_LLP_DAUGHTERS];
   float         gLLP_daughter_phi_ecalcorr[N_MAX_LLP_DAUGHTERS];
   float         gLLP_min_delta_r_match_jet[N_MAX_LLP_DAUGHTERS];
   unsigned int          gLLP_daughter_match_jet_index[N_MAX_LLP_DAUGHTERS];
*/
 //daughters
 bool gLLP_daughter_EB[N_MAX_LLP_DAUGHTERS]; 
 bool gLLP_daughter_ETL[N_MAX_LLP_DAUGHTERS];

 float gLLP_daughter_photon_travel_time_EB[N_MAX_LLP_DAUGHTERS];
 float gLLP_daughter_photon_travel_time_ETL[N_MAX_LLP_DAUGHTERS];

 float gLLP_daughter_travel_time_EB[N_MAX_LLP_DAUGHTERS];
 float gLLP_daughter_travel_time_ETL[N_MAX_LLP_DAUGHTERS];

 float gen_time_daughter_EB[N_MAX_LLP_DAUGHTERS];
 float gen_time_daughter_ETL[N_MAX_LLP_DAUGHTERS];

 int   gLLP_daughter_id[N_MAX_LLP_DAUGHTERS];
 float gLLP_daughter_pt[N_MAX_LLP_DAUGHTERS];
 float gLLP_daughter_eta[N_MAX_LLP_DAUGHTERS];
 float gLLP_daughter_phi[N_MAX_LLP_DAUGHTERS];
 float gLLP_daughter_eta_ecalcorr[N_MAX_LLP_DAUGHTERS];
 float gLLP_daughter_phi_ecalcorr[N_MAX_LLP_DAUGHTERS];
 float gLLP_daughter_e[N_MAX_LLP_DAUGHTERS];
 float gLLP_daughter_mass[N_MAX_LLP_DAUGHTERS];

 unsigned int gLLP_daughter_match_jet_index[N_MAX_LLP_DAUGHTERS];
 float gLLP_daughter_min_delta_r_match_jet[N_MAX_LLP_DAUGHTERS];
 unsigned int gLLP_daughter_match_calojet_index[N_MAX_LLP_DAUGHTERS];
 float gLLP_daughter_min_delta_r_match_calojet[N_MAX_LLP_DAUGHTERS];

 //grandaughters
 bool gLLP_grandaughter_EB[N_MAX_LLP_GRAND_DAUGHTERS]; 
 bool gLLP_grandaughter_ETL[N_MAX_LLP_GRAND_DAUGHTERS];

 float gLLP_grandaughter_photon_travel_time_EB[N_MAX_LLP_GRAND_DAUGHTERS];
 float gLLP_grandaughter_photon_travel_time_ETL[N_MAX_LLP_GRAND_DAUGHTERS];

 float gLLP_grandaughter_travel_time_EB[N_MAX_LLP_GRAND_DAUGHTERS];
 float gLLP_grandaughter_travel_time_ETL[N_MAX_LLP_GRAND_DAUGHTERS];

 float gen_time_grandaughter_EB[N_MAX_LLP_GRAND_DAUGHTERS];
 float gen_time_grandaughter_ETL[N_MAX_LLP_GRAND_DAUGHTERS];

 int   gLLP_grandaughter_id[N_MAX_LLP_GRAND_DAUGHTERS];
 float gLLP_grandaughter_pt[N_MAX_LLP_GRAND_DAUGHTERS];
 float gLLP_grandaughter_eta[N_MAX_LLP_GRAND_DAUGHTERS];
 float gLLP_grandaughter_phi[N_MAX_LLP_GRAND_DAUGHTERS];
 float gLLP_grandaughter_eta_ecalcorr[N_MAX_LLP_GRAND_DAUGHTERS];
 float gLLP_grandaughter_phi_ecalcorr[N_MAX_LLP_GRAND_DAUGHTERS];
 float gLLP_grandaughter_e[N_MAX_LLP_GRAND_DAUGHTERS];
 float gLLP_grandaughter_mass[N_MAX_LLP_GRAND_DAUGHTERS];

 unsigned int gLLP_grandaughter_match_jet_index[N_MAX_LLP_GRAND_DAUGHTERS];
 float gLLP_grandaughter_min_delta_r_match_jet[N_MAX_LLP_GRAND_DAUGHTERS];
 unsigned int gLLP_grandaughter_match_calojet_index[N_MAX_LLP_GRAND_DAUGHTERS];
 float gLLP_grandaughter_min_delta_r_match_calojet[N_MAX_LLP_GRAND_DAUGHTERS];

  //MET
 float metType1Pt;
 float metType1Px;
 float metType1Py;
 float metType1Eta;
 float metType1Phi; 

  //HLT 
  bool HLTDecision[NTriggersMAX];

  UInt_t wzevtNum,trig, trig_lepId, trig_lepId_dijet; //number of events that pass each criteria

  //functions
  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();

};
#endif
