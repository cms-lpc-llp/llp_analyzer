#ifndef JetTree_H
#define JetTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

  class JetTree {

    public:

      /// bit map
      /// DON'T CHANGE ORDER

      //*******************************************
      //=== JetTriggerBits  ====
      //*******************************************
      enum JetTriggerBits { kJetTrigger_Jet                                    = 0x000001
      };

      /// variables
      UInt_t                  fi_evt;
      Float_t                 fWeight;
      UInt_t                  fRunNumber;
      UInt_t                  fLumiSectionNumber;
      UInt_t                  fEventNumber;
      Bool_t                  fJetEventNumberParity;
      Float_t                 fJetGenPt;
      Float_t                 fJetGenEta;
      Float_t                 fJetGenPhi;
      Float_t                 fJetRawPt;
      Float_t                 fJetPt;
      Float_t                 fJetEta;
      Float_t                 fJetPhi;
      UInt_t                  fJetTriggerBit;
      Float_t                 fRho;
      Float_t                 fNVertices;
      Float_t                 fJetCSV;
      Float_t                 fJetCISV;
      Float_t                 fJetArea;
      Float_t                 fJetPileupId;
      Int_t                   fJetPartonFlavor;
      Float_t                 fJetJEC;
      //All Photons Match To the Jet (Take Seed RecHit as a reference)
      Int_t                   fJetNPhotons;
      Float_t                 fJetPhotonE[1000];
      Float_t                 fJetPhotonEta[1000];
      Float_t                 fJetPhotonPhi[1000];
      Float_t                 fJetPhotonSeedRecHitE[1000];
      Float_t                 fJetPhotonSeedRecHitEta[1000];
      Float_t                 fJetPhotonSeedRecHitPhi[1000];
      Float_t                 fJetPhotonSeedRecHitTime[1000];
      //All RecHits matched to the Jet through the Photons
      Int_t                   fJetNRecHits;
      Float_t                 fJetRecHitE[1000];
      Float_t                 fJetRecHitEta[1000];
      Float_t                 fJetRecHitPhi[1000];
      Float_t                 fJetRecHitTime[1000];

    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;

      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;

      /// default constructor
      JetTree()  {};
      /// default destructor
      ~JetTree(){
        if (f_) f_->Close();
      };

      /// initialize varibles and fill list of available variables
      void InitVariables()
      {
        fi_evt = 0;
        ResetVariables();
      }
      void ResetVariables() {
        fWeight			                     = 0.0;
        fRunNumber		                     = 0.0;
        fLumiSectionNumber	                     = 0.0;
        fEventNumber		                     = 0.0;
        fJetEventNumberParity 	                     = 0.0;
        fJetGenPt 		                     = 0.0;
        fJetGenEta 		                     = 0.0;
        fJetGenPhi 		                     = 0.0;
        fJetRawPt 		                     = 0.0;
        fJetPt 			                     = 0.0;
        fJetEta 		                     = 0.0;
        fJetPhi 		                     = 0.0;
        fJetTriggerBit		                     = 0.0;
        fRho  			                     = 0.0;
        fNVertices 		                     = 0.0;
        fJetCSV                                      = 0.0;
        fJetCISV                                     = 0.0;
        fJetArea                                     = 0.0;
        fJetPileupId                                 = 0.0;
        fJetPartonFlavor                             = 0;
        fJetJEC                                      = 0;
        fJetNPhotons                                 = 0;
        fJetNRecHits                                 = 0;
	for ( int i = 0; i < 1000; i++ ) {
    //All Photons
	  fJetPhotonE[i] = 0.0;
	  fJetPhotonEta[i] = 0.0;
	  fJetPhotonPhi[i] = 0.0;
    fJetPhotonSeedRecHitE[i]  = 0.0;
    fJetPhotonSeedRecHitEta[i]  = 0.0;
    fJetPhotonSeedRecHitPhi[i]  = 0.0;
	  fJetPhotonSeedRecHitTime[i] = 0.0;
    //All RecHits
    fJetRecHitE[i]    = 0.0;
    fJetRecHitEta[i]  = 0.0;
    fJetRecHitPhi[i]  = 0.0;
	  fJetRecHitTime[i] = 0.0;
	}
      }

      /// load a JetTree
      void LoadTree(const char* file){
        f_ = TFile::Open(file);
        assert(f_);
        tree_ = dynamic_cast<TTree*>(f_->Get("Jets"));
        assert(tree_);
      }

      /// create a JetTree
      void CreateTree(){
        tree_ = new TTree("Jets","Jets");
        f_ = 0;

        //book the branches
        tree_->Branch("i_evt",&fi_evt,"i_evt/i");
        tree_->Branch("weight",&fWeight,"weight/F");
        tree_->Branch("run",&fRunNumber,"run/i");
        tree_->Branch("lumi",&fLumiSectionNumber,"lumi/i");
        tree_->Branch("event",&fEventNumber,"event/i");
        tree_->Branch("EventNumberParity",&fJetEventNumberParity,"EventNumberParity/O");
        tree_->Branch("genpt",&fJetGenPt,"genpt/F");
        tree_->Branch("geneta",&fJetGenEta,"geneta/F");
        tree_->Branch("genphi",&fJetGenPhi,"genphi/F");
        tree_->Branch("rawpt",&fJetRawPt,"rawpt/F");
        tree_->Branch("pt",&fJetPt,"pt/F");
        tree_->Branch("eta",&fJetEta,"eta/F");
        tree_->Branch("phi",&fJetPhi,"phi/F");
        tree_->Branch("triggerBit",&fJetTriggerBit,"triggerBit/i");
        tree_->Branch("rho",&fRho,"rho/F");
        tree_->Branch("vertices",&fNVertices,"vertices/F");
        tree_->Branch("CSV",&fJetCSV,"CSV/F");
        tree_->Branch("CISV",&fJetCISV,"CISV/F");
        tree_->Branch("Area",&fJetArea,"Area/F");
        tree_->Branch("PileupId",&fJetPileupId,"PileupId/F");
        tree_->Branch("PartonFlavor",&fJetPartonFlavor,"PartonFlavor/I");
        tree_->Branch("JEC",&fJetJEC,"JEC/F");
        //All Photons Matched to Jets
        tree_->Branch("JetNPhotons", &fJetNPhotons, "JetNPhotons/I");
        tree_->Branch("JetPhotonE", fJetPhotonE, "JetPhotonE[JetNPhotons]/F");
        tree_->Branch("JetPhotonEta", fJetPhotonEta, "JetPhotonEta[JetNPhotons]/F");
        tree_->Branch("JetPhotonPhi", fJetPhotonPhi, "JetPhotonPhi[JetNPhotons]/F");
        tree_->Branch("JetPhotonSeedRecHitE", fJetPhotonSeedRecHitE, "JetPhotonSeedRecHitE[JetNPhotons]/F");
        tree_->Branch("JetPhotonSeedRecHitEta", fJetPhotonSeedRecHitEta, "JetPhotonSeedRecHitEta[JetNPhotons]/F");
        tree_->Branch("JetPhotonSeedRecHitPhi", fJetPhotonSeedRecHitPhi, "JetPhotonSeedRecHitPhi[JetNPhotons]/F");
        tree_->Branch("JetPhotonSeedRecHitTime", fJetPhotonSeedRecHitTime, "JetPhotonSeedRecHitTime[JetNPhotons]/F");
        //All RecHits Matched to Jets
        tree_->Branch("JetNRecHits", &fJetNRecHits, "JetNRecHits/I");
        tree_->Branch("JetRecHitE", fJetRecHitE, "JetRecHitE[JetNRecHits]/F");
        tree_->Branch("JetRecHitEta", fJetRecHitEta, "JetRecHitEta[JetNRecHits]/F");
        tree_->Branch("JetRecHitPhi", fJetRecHitPhi, "JetRecHitPhi[JetNRecHits]/F");
        tree_->Branch("JetRecHitTime", fJetRecHitTime, "JetRecHitTime[JetNRecHits]/F");
      }

      // initialze a JetTree
      void InitTree(){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;

        tree_->SetBranchAddress("i_evt",&fi_evt);
        tree_->SetBranchAddress("weight",&fWeight);
        tree_->SetBranchAddress("run",&fRunNumber);
        tree_->SetBranchAddress("lumi",&fLumiSectionNumber);
        tree_->SetBranchAddress("event",&fEventNumber);
        tree_->SetBranchAddress("EventNumberParity",&fJetEventNumberParity);
        tree_->SetBranchAddress("genpt",&fJetGenPt);
        tree_->SetBranchAddress("geneta",&fJetGenEta);
        tree_->SetBranchAddress("genphi",&fJetGenPhi);
        tree_->SetBranchAddress("rawpt",&fJetRawPt);
        tree_->SetBranchAddress("pt",&fJetPt);
        tree_->SetBranchAddress("eta",&fJetEta);
        tree_->SetBranchAddress("phi",&fJetPhi);
        tree_->SetBranchAddress("triggerBit",&fJetTriggerBit);
        tree_->SetBranchAddress("rho",&fRho);
        tree_->SetBranchAddress("vertices",&fNVertices);
        tree_->SetBranchAddress("CSV",&fJetCSV);
        tree_->SetBranchAddress("CISV",&fJetCISV);
        tree_->SetBranchAddress("Area",&fJetArea);
        tree_->SetBranchAddress("PileupId",&fJetPileupId);
        tree_->SetBranchAddress("PartonFlavor",&fJetPartonFlavor);
        tree_->SetBranchAddress("JEC",&fJetJEC);
        tree_->SetBranchAddress("JetNPhotons", &fJetNPhotons);
        tree_->SetBranchAddress("JetPhotonE", fJetPhotonE);
        tree_->SetBranchAddress("JetPhotonEta", fJetPhotonEta);
        tree_->SetBranchAddress("JetPhotonPhi", fJetPhotonPhi);
        tree_->SetBranchAddress("JetPhotonSeedRecHitE", fJetPhotonSeedRecHitE);
        tree_->SetBranchAddress("JetPhotonSeedRecHitEta", fJetPhotonSeedRecHitEta);
        tree_->SetBranchAddress("JetPhotonSeedRecHitPhi", fJetPhotonSeedRecHitPhi);
        tree_->SetBranchAddress("JetPhotonSeedRecHitTime", fJetPhotonSeedRecHitTime);
        //All RecHits Matched to Jet
        tree_->SetBranchAddress("JetNRecHits", &fJetNRecHits);
        tree_->SetBranchAddress("JetRecHitE", fJetRecHitE);
        tree_->SetBranchAddress("JetRecHitEta", fJetRecHitEta);
        tree_->SetBranchAddress("JetRecHitPhi", fJetRecHitPhi);
        tree_->SetBranchAddress("JetRecHitTime", fJetRecHitTime);
        gErrorIgnoreLevel = currentState;
      }

  };


#endif
