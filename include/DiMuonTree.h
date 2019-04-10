#ifndef DiMuonTree_H
#define DiMuonTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

  class DiMuonTree {

    public:

      /// bit map
      /// DON'T CHANGE ORDER

      //*******************************************
      //=== DiMuonTriggerBits  ====
      //*******************************************
      enum DiMuonTriggerBits { kDiMuonTrigger_DiMuon                                    = 0x000001
      };

      /// variables
      UInt_t                  f_i_evt;
      Float_t                 f_weight;
      UInt_t                  f_run_number;
      UInt_t                  f_lumi_section_number;
      UInt_t                  f_event_number;
      //DiMuon
      Float_t                 f_mumu_mass;
      Float_t                 f_mumu_pt;
      Float_t                 f_mumu_eta;
      Float_t                 f_mumu_phi;
      //Muons
      Float_t                 f_mu_pt[2];
      Float_t                 f_mu_eta[2];
      Float_t                 f_mu_phi[2];
      Int_t                   f_mu_type[2];
      UInt_t                  f_mu_pass_sel[2];
      Float_t                 f_mu_chargedMiniIso[2];
      Float_t                 f_mu_photonAndNeutralHadronMiniIso[2];
      Float_t                 f_mu_d0[2];
      Float_t                 f_mu_dZ[2];
      Float_t                 f_mu_ip3d[2];
      Float_t                 f_mu_ip3dSignificance[2];
      //
      Float_t                 f_met_pt;
      Float_t                 f_met_phi;

    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;

      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;

      /// default constructor
      DiMuonTree()  {};
      /// default destructor
      ~DiMuonTree(){
        if (f_) f_->Close();
      };

      /// initialize varibles and fill list of available variables
      void InitVariables()
      {
        f_i_evt = 0;
        ResetVariables();
      }
      void ResetVariables() {
        f_weight			                     = 0.0;
        f_run_number		                     = 0.0;
        f_lumi_section_number	                     = 0.0;
        f_event_number		                     = 0.0;
        //
        f_mumu_mass 	                     = 0.0;
        f_mumu_pt 		                     = 0.0;
        f_mumu_eta 		                     = 0.0;
        f_mumu_phi 		                     = 0.0;
        //
        for ( int i = 0; i < 2; i++ )
        {
          f_mu_pt[i]  = 0.0;
          f_mu_eta[i] = 0.0;
          f_mu_phi[i] = 0.0;
          f_mu_type[i] = 0.0;
          f_mu_pass_sel[i] = 0.0;
          f_mu_chargedMiniIso[2] = 0.0;
         f_mu_photonAndNeutralHadronMiniIso[2] = 0.0;
         f_mu_d0[2] = 0.0;
         f_mu_dZ[2] = 0.0;
         f_mu_ip3d[2] = 0.0;
         f_mu_ip3dSignificance[2] = 0.0;
        }
        //
        f_met_pt 		                     = 0.0;
        f_met_phi 			                 = 0.0;
      }

      /// load a loadTree
      void LoadTree(const char* file){
        f_ = TFile::Open(file);
        assert(f_);
        tree_ = dynamic_cast<TTree*>(f_->Get("DiMuon"));
        assert(tree_);
      }

      /// create a JetTree
      void CreateTree(){
        tree_ = new TTree("DiMuon","DiMuon");
        f_ = 0;

        //book the branches
        tree_->Branch("i_evt",&f_i_evt,"i_evt/i");
        tree_->Branch("weight",&f_weight,"weight/F");
        tree_->Branch("run",&f_run_number,"run/i");
        tree_->Branch("lumi",&f_lumi_section_number,"lumi/i");
        tree_->Branch("event",&f_event_number,"event/i");
        //
        tree_->Branch("mumu_mass",&f_mumu_mass,"mumu_mass/F");
        tree_->Branch("mumu_pt",&f_mumu_pt,"mumu_pt/F");
        tree_->Branch("mumu_eta",&f_mumu_eta,"mumu_eta/F");
        tree_->Branch("mumu_phi",&f_mumu_phi,"mumu_phi/F");
        //
        tree_->Branch("mu_pt",f_mu_pt,"mu_pt[2]/F");
        tree_->Branch("mu_eta",f_mu_eta,"mu_eta[2]/F");
        tree_->Branch("mu_phi",f_mu_phi,"mu_phi[2]/F");
        tree_->Branch("mu_type",f_mu_type,"mu_type[2]/I");
        tree_->Branch("mu_pass_sel",f_mu_pass_sel,"mu_pass_sel[2]/i");
        tree_->Branch("mu_chargedMiniIso",f_mu_chargedMiniIso,"mu_chargedMiniIso[2]/F");
        tree_->Branch("mu_photonAndNeutralHadronMiniIso",f_mu_photonAndNeutralHadronMiniIso,"mu_photonAndNeutralHadronMiniIso[2]/F");
        tree_->Branch("mu_d0",f_mu_d0,"mu_d0[2]/F");
        tree_->Branch("mu_dZ",f_mu_dZ,"mu_dZ[2]/F");
        tree_->Branch("mu_ip3d",f_mu_ip3d,"mu_ip3d[2]/F");
        tree_->Branch("mu_ip3dSignificance",f_mu_ip3dSignificance,"mu_ip3dSignificance[2]/F");
        //
        tree_->Branch("met_pt",&f_met_pt,"met_pt/F");
        tree_->Branch("met_phi",&f_met_phi,"met_phi/F");
      }

      // initialze a JetTree
      void InitTree(){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;
        tree_->SetBranchAddress("i_evt",&f_i_evt);
        tree_->SetBranchAddress("weight",&f_weight);
        tree_->SetBranchAddress("run",&f_run_number);
        tree_->SetBranchAddress("lumi",&f_lumi_section_number);
        tree_->SetBranchAddress("event",&f_event_number);
        //
        tree_->SetBranchAddress("mumu_mass",&f_mumu_mass);
        tree_->SetBranchAddress("mumu_pt",&f_mumu_pt);
        tree_->SetBranchAddress("mumu_eta",&f_mumu_eta);
        tree_->SetBranchAddress("mumu_phi",&f_mumu_phi);
        //
        tree_->SetBranchAddress("mu_pt",f_mu_pt);
        tree_->SetBranchAddress("mu_eta",f_mu_eta);
        tree_->SetBranchAddress("mu_phi",f_mu_phi);
        tree_->SetBranchAddress("mu_type",f_mu_type);
        tree_->SetBranchAddress("mu_pass_sel",f_mu_pass_sel);
        tree_->SetBranchAddress("mu_chargedMiniIso",f_mu_chargedMiniIso);
        tree_->SetBranchAddress("mu_photonAndNeutralHadronMiniIso",f_mu_photonAndNeutralHadronMiniIso);
        tree_->SetBranchAddress("mu_d0",f_mu_d0);
        tree_->SetBranchAddress("mu_dZ",f_mu_dZ);
        tree_->SetBranchAddress("mu_ip3d",f_mu_ip3d);
        tree_->SetBranchAddress("mu_ip3dSignificance",f_mu_ip3dSignificance);
        //
        tree_->SetBranchAddress("met_pt",&f_met_pt);
        tree_->SetBranchAddress("met_phi",&f_met_phi);

        gErrorIgnoreLevel = currentState;
      }

  };


#endif
