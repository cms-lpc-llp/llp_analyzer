#ifndef PhotonTree_H
#define PhotonTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

  class PhotonTree {

    public:

      //******************************************* 
      //=== tree versions ===
      //*******************************************
      enum PhotonTreeVersion { kPhotonTreeStd, kPhotonTreeLight };

      /// variables
      Float_t                 fWeight;
      Float_t		      fpvZ_New;
      Float_t		      fpvdT_New;
      Float_t		      fpvdZ_New;
      Float_t		      fpvZ_Gen;
      Float_t		      fPV_ndof;
      UInt_t                  fRunNumber;
      UInt_t                  fLumiSectionNumber;
      UInt_t                  fEventNumber;
      Bool_t                  fPhoEventNumberParity;
      Int_t                   fNearstHighPtGenPdgId;
      Int_t                   fNearstHighPtGenMotherId;
      Int_t                   fMotherPdgId;
      Float_t                 fDRToClosestParton;
      Float_t                 fPhoGenE;  
      Float_t                 fPhoGenPt; 
      Float_t                 fPhoGenEta; 
      Float_t                 fPhoGenPhi; 
      Float_t                 fRho; 
      Float_t                 fNVertices;

      Float_t                 fPhoPt; 
      Float_t                 fPhoEta; 
      Float_t                 fPhoPhi; 
      Float_t                 fPhoSigmaIetaIeta;
      Float_t                 fPhoFull5x5SigmaIetaIeta;
      Float_t                 fPhoR9;
      Float_t                 fPhoHOverE;
      Float_t                 fPhoChargedHadronIso;
      Float_t                 fPhoChargedHadronIso_NewPV_NoTiming;
      Float_t                 fPhoChargedHadronIso_NewPV_Timing50_TrkVtx;
      Float_t                 fPhoChargedHadronIso_NewPV_Timing80_TrkVtx;
      Float_t                 fPhoChargedHadronIso_NewPV_Timing100_TrkVtx;
      Float_t                 fPhoChargedHadronIso_NewPV_Timing120_TrkVtx;
      Float_t                 fPhoChargedHadronIso_NewPV_Timing50_TrkPho;
      Float_t                 fPhoChargedHadronIso_NewPV_Timing80_TrkPho;
      Float_t                 fPhoChargedHadronIso_NewPV_Timing100_TrkPho;
      Float_t                 fPhoChargedHadronIso_NewPV_Timing120_TrkPho;
 
      Float_t                 fPhoNeutralHadronIso;
      Float_t                 fPhoPhotonIso;
      Bool_t                  fPhoIsConversion;
      Bool_t                  fPhoPassEleVeto;
      Float_t                 fPhoRegressionE;
      Float_t                 fPhoRegressionSigmaE;
      Float_t                 fPhoIDMVA;
      Float_t                 fPhoSCEta;
      Bool_t                  fPhoHasPixelSeed;

      //ID bits
      Bool_t                  fPhoIsLoose;
      Bool_t                  fPhoIsMedium;
      Bool_t                  fPhoIsTight;
      Bool_t                  fPhoPassLooseID;
      Bool_t                  fPhoPassMediumID;
      Bool_t                  fPhoPassTightID;
      Bool_t                  fPhoPassLooseIso;
      Bool_t                  fPhoPassMediumIso;
      Bool_t                  fPhoPassTightIso;

    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;
  
      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;

      /// default constructor  
      PhotonTree()  {};
      /// default destructor
      ~PhotonTree(){ 
        if (f_) f_->Close();  
      };
    
      /// initialize varibles and fill list of available variables
      void InitVariables() {
	fWeight                  = 0.0;
	fpvZ_New                  = 0.0;
	fpvdZ_New                  = 0.0;
	fpvdT_New                  = 0.0;
	fpvZ_Gen                  = 0.0;
	fPV_ndof                  = 0.0;
	fRunNumber               = 0.0;
	fLumiSectionNumber       = 0.0;
	fEventNumber             = 0.0;
	fPhoEventNumberParity     = false;
	fMotherPdgId                   = 0.0;
	fNearstHighPtGenPdgId                   = 0.0;
	fNearstHighPtGenMotherId                   = 0.0;
	fDRToClosestParton       = 0.0;
	fPhoGenE                 = 0.0;
	fPhoGenPt                = 0.0;
	fPhoGenEta               = 0.0;
	fPhoGenPhi               = 0.0;
	fRho                     = 0.0;
	fNVertices               = 0.0;
	fPhoPt                   = 0.0;
	fPhoEta                  = 0.0;
	fPhoPhi                  = 0.0;
	fPhoSigmaIetaIeta        = 0.0;
	fPhoFull5x5SigmaIetaIeta = 0.0;
	fPhoR9                   = 0.0;
	fPhoHOverE               = 0.0;
	fPhoChargedHadronIso     = 0.0;
	fPhoChargedHadronIso_NewPV_NoTiming     = 0.0;
	fPhoChargedHadronIso_NewPV_Timing50_TrkVtx     = 0.0;
	fPhoChargedHadronIso_NewPV_Timing80_TrkVtx     = 0.0;
	fPhoChargedHadronIso_NewPV_Timing100_TrkVtx     = 0.0;
	fPhoChargedHadronIso_NewPV_Timing120_TrkVtx     = 0.0;
	fPhoChargedHadronIso_NewPV_Timing50_TrkPho     = 0.0;
	fPhoChargedHadronIso_NewPV_Timing80_TrkPho     = 0.0;
	fPhoChargedHadronIso_NewPV_Timing100_TrkPho     = 0.0;
	fPhoChargedHadronIso_NewPV_Timing120_TrkPho     = 0.0;
	
	fPhoNeutralHadronIso     = 0.0;
	fPhoPhotonIso            = 0.0;
	fPhoIsConversion         = false;
	fPhoPassEleVeto          = false;
	fPhoRegressionE          = 0.0;
	fPhoRegressionSigmaE     = 0.0;
	fPhoIDMVA                = 0.0;
	fPhoSCEta                = 0.0;
	fPhoHasPixelSeed         = false;
	fPhoIsLoose              = false;
	fPhoIsMedium             = false;
	fPhoIsTight              = false;
	fPhoPassLooseID         = false;
	fPhoPassMediumID        = false;
	fPhoPassTightID         = false;
	fPhoPassLooseIso         = false;
	fPhoPassMediumIso        = false;
	fPhoPassTightIso         = false;
      }
    
      /// load a PhotonTree
      void LoadTree(const char* file){
        f_ = TFile::Open(file);
        assert(f_);
        tree_ = dynamic_cast<TTree*>(f_->Get("Photons"));
        assert(tree_);
      }
    
      /// create a PhotonTree
      void CreateTree(int version=kPhotonTreeStd){
        tree_ = new TTree("Photons","Photons");
        f_ = 0;

        //book the branches
	tree_->Branch("weight",&fWeight,"weight/F");
	tree_->Branch("pvZ_New",&fpvZ_New,"pvZ_New/F");
	tree_->Branch("pvdZ_New",&fpvdZ_New,"pvdZ_New/F");
	tree_->Branch("pvdT_New",&fpvdT_New,"pvdT_New/F");
	tree_->Branch("pvZ_Gen",&fpvZ_Gen,"pvZ_Gen/F");
	tree_->Branch("PV_ndof",&fPV_ndof,"PV_ndof/F");
	tree_->Branch("run",&fRunNumber,"run/i");
	tree_->Branch("lumi",&fLumiSectionNumber,"lumi/i");
	tree_->Branch("event",&fEventNumber,"event/i");
	tree_->Branch("PhoEventNumberParity",&fPhoEventNumberParity,"PhoEventNumberParity/O");
	tree_->Branch("MotherPdgId",&fMotherPdgId,"MotherPdgId/I");
	tree_->Branch("NearstHighPtGenPdgId",&fNearstHighPtGenPdgId,"NearstHighPtGenPdgId/I");
	tree_->Branch("NearstHighPtGenMotherId",&fNearstHighPtGenMotherId,"NearstHighPtGenMotherId/I");
	tree_->Branch("DRToClosestParton",&fDRToClosestParton,"DRToClosestParton/F");
	tree_->Branch("PhoGenE",&fPhoGenE,"PhoGenE/F");
	tree_->Branch("PhoGenPt",&fPhoGenPt,"PhoGenPt/F");
	tree_->Branch("PhoGenEta",&fPhoGenEta,"PhoGenEta/F");
	tree_->Branch("PhoGenPhi",&fPhoGenPhi,"PhoGenPhi/F");
	tree_->Branch("Rho",&fRho,"Rho/F");
	tree_->Branch("NVertices",&fNVertices,"NVertices/F");
	tree_->Branch("PhoPt",&fPhoPt,"PhoPt/F");
	tree_->Branch("PhoEta",&fPhoEta,"PhoEta/F");
	tree_->Branch("PhoPhi",&fPhoPhi,"PhoPhi/F");
	tree_->Branch("PhoSigmaIetaIeta",&fPhoSigmaIetaIeta,"PhoSigmaIetaIeta/F");
	tree_->Branch("PhoFull5x5SigmaIetaIeta",&fPhoFull5x5SigmaIetaIeta,"PhoFull5x5SigmaIetaIeta/F");
	tree_->Branch("PhoR9",&fPhoR9,"PhoR9/F");
	tree_->Branch("PhoHOverE",&fPhoHOverE,"PhoHOverE/F");
	tree_->Branch("PhoChargedHadronIso",&fPhoChargedHadronIso,"PhoChargedHadronIso/F");
	tree_->Branch("PhoChargedHadronIso_NewPV_NoTiming",&fPhoChargedHadronIso_NewPV_NoTiming,"PhoChargedHadronIso_NewPV_NoTiming/F");
	tree_->Branch("PhoChargedHadronIso_NewPV_Timing50_TrkVtx",&fPhoChargedHadronIso_NewPV_Timing50_TrkVtx,"PhoChargedHadronIso_NewPV_Timing50_TrkVtx/F");
	tree_->Branch("PhoChargedHadronIso_NewPV_Timing80_TrkVtx",&fPhoChargedHadronIso_NewPV_Timing80_TrkVtx,"PhoChargedHadronIso_NewPV_Timing80_TrkVtx/F");
	tree_->Branch("PhoChargedHadronIso_NewPV_Timing100_TrkVtx",&fPhoChargedHadronIso_NewPV_Timing100_TrkVtx,"PhoChargedHadronIso_NewPV_Timing100_TrkVtx/F");
	tree_->Branch("PhoChargedHadronIso_NewPV_Timing120_TrkVtx",&fPhoChargedHadronIso_NewPV_Timing120_TrkVtx,"PhoChargedHadronIso_NewPV_Timing120_TrkVtx/F");
	tree_->Branch("PhoChargedHadronIso_NewPV_Timing50_TrkPho",&fPhoChargedHadronIso_NewPV_Timing50_TrkPho,"PhoChargedHadronIso_NewPV_Timing50_TrkPho/F");
	tree_->Branch("PhoChargedHadronIso_NewPV_Timing80_TrkPho",&fPhoChargedHadronIso_NewPV_Timing80_TrkPho,"PhoChargedHadronIso_NewPV_Timing80_TrkPho/F");
	tree_->Branch("PhoChargedHadronIso_NewPV_Timing100_TrkPho",&fPhoChargedHadronIso_NewPV_Timing100_TrkPho,"PhoChargedHadronIso_NewPV_Timing100_TrkPho/F");
	tree_->Branch("PhoChargedHadronIso_NewPV_Timing120_TrkPho",&fPhoChargedHadronIso_NewPV_Timing120_TrkPho,"PhoChargedHadronIso_NewPV_Timing120_TrkPho/F");

	tree_->Branch("PhoNeutralHadronIso",&fPhoNeutralHadronIso,"PhoNeutralHadronIso/F");
	tree_->Branch("PhoPhotonIso",&fPhoPhotonIso,"PhoPhotonIso/F");
	tree_->Branch("PhoIsConversion",&fPhoIsConversion,"PhoIsConversion/O");
	tree_->Branch("PhoPassEleVeto",&fPhoPassEleVeto,"PhoPassEleVeto/O");
	tree_->Branch("PhoRegressionE",&fPhoRegressionE,"PhoRegressionE/F");
	tree_->Branch("PhoRegressionSigmaE",&fPhoRegressionSigmaE,"PhoRegressionSigmaE/F");
	tree_->Branch("PhoIDMVA",&fPhoIDMVA,"PhoIDMVA/F");
	tree_->Branch("PhoSCEta",&fPhoSCEta,"PhoSCEta/F");
	tree_->Branch("PhoHasPixelSeed",&fPhoHasPixelSeed,"PhoHasPixelSeed/O");
	tree_->Branch("PhoIsLoose",&fPhoIsLoose,"PhoIsLoose/O");
	tree_->Branch("PhoIsMedium",&fPhoIsMedium,"PhoIsMedium/O");
	tree_->Branch("PhoIsTight",&fPhoIsTight,"PhoIsTight/O");
	tree_->Branch("PhoPassLooseID",&fPhoPassLooseID,"PhoPassLooseID/O");
	tree_->Branch("PhoPassMediumID",&fPhoPassMediumID,"PhoPassMediumID/O");
	tree_->Branch("PhoPassTightID",&fPhoPassTightID,"PhoPassTightID/O");
	tree_->Branch("PhoPassLooseIso",&fPhoPassLooseIso,"PhoPassLooseIso/O");
	tree_->Branch("PhoPassMediumIso",&fPhoPassMediumIso,"PhoPassMediumIso/O");
	tree_->Branch("PhoPassTightIso",&fPhoPassTightIso,"PhoPassTightIso/O");

	if (version == kPhotonTreeStd ) {

        }
      } 

      // initialze a PhotonTree
      void InitTree(int version=kPhotonTreeStd){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;
        // gErrorIgnoreLevel = kError;
        gErrorIgnoreLevel = kBreak;

	tree_->SetBranchAddress("weight",&fWeight);
	tree_->SetBranchAddress("pvZ_New",&fpvZ_New);
	tree_->SetBranchAddress("pvdZ_New",&fpvdZ_New);
	tree_->SetBranchAddress("pvdT_New",&fpvdT_New);
	tree_->SetBranchAddress("pvZ_Gen",&fpvZ_Gen);
	tree_->SetBranchAddress("PV_ndof",&fPV_ndof);
	tree_->SetBranchAddress("run",&fRunNumber);
	tree_->SetBranchAddress("lumi",&fLumiSectionNumber);
	tree_->SetBranchAddress("event",&fEventNumber);
	tree_->SetBranchAddress("PhoEventNumberParity",&fPhoEventNumberParity);
	tree_->SetBranchAddress("MotherPdgId",&fMotherPdgId);
	tree_->SetBranchAddress("NearstHighPtGenPdgId",&fNearstHighPtGenPdgId);
	tree_->SetBranchAddress("NearstHighPtGenMotherId",&fNearstHighPtGenMotherId);
	tree_->SetBranchAddress("DRToClosestParton",&fDRToClosestParton);
	tree_->SetBranchAddress("PhoGenE",&fPhoGenE);
	tree_->SetBranchAddress("PhoGenPt",&fPhoGenPt);
	tree_->SetBranchAddress("PhoGenEta",&fPhoGenEta);
	tree_->SetBranchAddress("PhoGenPhi",&fPhoGenPhi);
	tree_->SetBranchAddress("Rho",&fRho);
	tree_->SetBranchAddress("NVertices",&fNVertices);
	tree_->SetBranchAddress("PhoPt",&fPhoPt);
	tree_->SetBranchAddress("PhoEta",&fPhoEta);
	tree_->SetBranchAddress("PhoPhi",&fPhoPhi);
	tree_->SetBranchAddress("PhoSigmaIetaIeta",&fPhoSigmaIetaIeta);
	tree_->SetBranchAddress("PhoFull5x5SigmaIetaIeta",&fPhoFull5x5SigmaIetaIeta);
	tree_->SetBranchAddress("PhoR9",&fPhoR9);
	tree_->SetBranchAddress("PhoHOverE",&fPhoHOverE);
	tree_->SetBranchAddress("PhoChargedHadronIso",&fPhoChargedHadronIso);
	tree_->SetBranchAddress("PhoChargedHadronIso_NewPV_NoTiming",&fPhoChargedHadronIso_NewPV_NoTiming);
	tree_->SetBranchAddress("PhoChargedHadronIso_NewPV_Timing50_TrkVtx",&fPhoChargedHadronIso_NewPV_Timing50_TrkVtx);
	tree_->SetBranchAddress("PhoChargedHadronIso_NewPV_Timing80_TrkVtx",&fPhoChargedHadronIso_NewPV_Timing80_TrkVtx);
	tree_->SetBranchAddress("PhoChargedHadronIso_NewPV_Timing100_TrkVtx",&fPhoChargedHadronIso_NewPV_Timing100_TrkVtx);
	tree_->SetBranchAddress("PhoChargedHadronIso_NewPV_Timing120_TrkVtx",&fPhoChargedHadronIso_NewPV_Timing120_TrkVtx);
	tree_->SetBranchAddress("PhoChargedHadronIso_NewPV_Timing50_TrkPho",&fPhoChargedHadronIso_NewPV_Timing50_TrkPho);
	tree_->SetBranchAddress("PhoChargedHadronIso_NewPV_Timing80_TrkPho",&fPhoChargedHadronIso_NewPV_Timing80_TrkPho);
	tree_->SetBranchAddress("PhoChargedHadronIso_NewPV_Timing100_TrkPho",&fPhoChargedHadronIso_NewPV_Timing100_TrkPho);
	tree_->SetBranchAddress("PhoChargedHadronIso_NewPV_Timing120_TrkPho",&fPhoChargedHadronIso_NewPV_Timing120_TrkPho);
	
	tree_->SetBranchAddress("PhoNeutralHadronIso",&fPhoNeutralHadronIso);
	tree_->SetBranchAddress("PhoPhotonIso",&fPhoPhotonIso);
	tree_->SetBranchAddress("PhoIsConversion",&fPhoIsConversion);
	tree_->SetBranchAddress("PhoPassEleVeto",&fPhoPassEleVeto);
	tree_->SetBranchAddress("PhoRegressionE",&fPhoRegressionE);
	tree_->SetBranchAddress("PhoRegressionSigmaE",&fPhoRegressionSigmaE);
	tree_->SetBranchAddress("PhoIDMVA",&fPhoIDMVA);
	tree_->SetBranchAddress("PhoSCEta",&fPhoSCEta);
	tree_->SetBranchAddress("PhoHasPixelSeed",&fPhoHasPixelSeed);
	tree_->SetBranchAddress("PhoIsLoose",&fPhoIsLoose);
	tree_->SetBranchAddress("PhoIsMedium",&fPhoIsMedium);
	tree_->SetBranchAddress("PhoIsTight",&fPhoIsTight);
	tree_->SetBranchAddress("PhoPassLooseID",&fPhoPassLooseID);
	tree_->SetBranchAddress("PhoPassMediumID",&fPhoPassMediumID);
	tree_->SetBranchAddress("PhoPassTightID",&fPhoPassTightID);
	tree_->SetBranchAddress("PhoPassLooseIso",&fPhoPassLooseIso);
	tree_->SetBranchAddress("PhoPassMediumIso",&fPhoPassMediumIso);
	tree_->SetBranchAddress("PhoPassTightIso",&fPhoPassTightIso);

        if (version == kPhotonTreeStd ) {
        }
        gErrorIgnoreLevel = currentState;
      }

  }; 


#endif
