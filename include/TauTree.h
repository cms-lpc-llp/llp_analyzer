#ifndef TauTree_H
#define TauTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

class TauTree {
  
 public:
  
  //******************************************* 
  //=== tree versions ===
  //*******************************************
  enum TauTreeVersion {kTauTreeStd, kTauTreeLight};
  
  /// variables
  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Bool_t                  fTauEventNumberParity;
  Float_t                 fRho; 
  Float_t                 fNVertices; 
  Float_t                 fTauGenPt; 
  Float_t                 fTauGenEta; 
  Float_t                 fTauGenPhi; 
  Float_t                 fTauPt; 
  Float_t                 fTauEta; 
  Float_t                 fTauPhi; 
  Int_t                   fPdgId;
  Float_t                 fDRToClosestParton;
  
  //ID bits
  Bool_t                  fTauIsTight;
  Bool_t                  fTauIsMedium; 
  Bool_t                  fTauIsLoose;
  
  // Typical Selection Working Points
  Bool_t                  fPassVetoSelection;
  Bool_t                  fPassLooseSelection;
  Bool_t                  fPassTightSelection;
  
  //additional variables      
  UInt_t  	          fTauTypeBits;
  Bool_t  		  fIsAllArbitrated;
  Float_t                 fTauTkNchi2; 
  Float_t                 fTauGlobalNchi2; 
  Float_t                 fTauNValidHits; 
  Float_t                 fTauNTrackerHits; 
  Float_t                 fTauNPixelHits; 
  Float_t                 fTauNMatches; 
  Float_t                 fTauD0;    
  Float_t                 fTauIP3d; 
  Float_t                 fTauIP3dSig; 
  Float_t                 fTauTrkKink; 
  Float_t                 fTauGlobalKink; 
  Float_t                 fTauSegmentCompatibility; 
  Float_t                 fTauCaloCompatibility; 
  Float_t                 fTauHadEnergy; 
  Float_t                 fTauHoEnergy; 
  Float_t                 fTauEmEnergy; 
  Float_t                 fTauHadS9Energy; 
  Float_t                 fTauHoS9Energy; 
  Float_t                 fTauEmS9Energy; 
  Float_t                 fTauTrkIso03; 
  Float_t                 fTauEMIso03; 
  Float_t                 fTauHadIso03; 
  Float_t                 fTauTrkIso05; 
  Float_t                 fTauEMIso05; 
  Float_t                 fTauHadIso05; 
  Float_t                 fTauPFIso04; 
  Float_t                 fChargedIso_DR0p0To0p1;
  Float_t                 fChargedIso_DR0p1To0p2;
  Float_t                 fChargedIso_DR0p2To0p3;
  Float_t                 fChargedIso_DR0p3To0p4;
  Float_t                 fChargedIso_DR0p4To0p5;
  Float_t                 fGammaIso_DR0p0To0p1;
  Float_t                 fGammaIso_DR0p1To0p2;
  Float_t                 fGammaIso_DR0p2To0p3;
  Float_t                 fGammaIso_DR0p3To0p4;
  Float_t                 fGammaIso_DR0p4To0p5;
  Float_t                 fNeutralHadronIso_DR0p0To0p1;
  Float_t                 fNeutralHadronIso_DR0p1To0p2;
  Float_t                 fNeutralHadronIso_DR0p2To0p3;
  Float_t                 fNeutralHadronIso_DR0p3To0p4;
  Float_t                 fNeutralHadronIso_DR0p4To0p5;
  Bool_t                  fTauPassTriggerDenominator;

 private:
  //Variables
  float _minPt;
  float _maxEta;
  
 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  
  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;
  
  /// default constructor  
  TauTree()  {};
  /// default destructor
  ~TauTree(){ 
    if (f_) f_->Close();  
  };
   
  float GetMinPt(){return _minPt;};
  float GetMaxEta(){return _maxEta;};

  void SetMinPt(float minPt){_minPt = minPt;};
  void SetMaxEta(float maxEta){_maxEta = maxEta;};

  //initialize varibles and fill list of available variables
  void InitVariables() {
    fWeight			   = 0.0;
    fRunNumber		           = 0.0;
    fLumiSectionNumber	           = 0.0;
    fEventNumber		   = 0.0;
    fTauEventNumberParity	   = 0.0;
    fRho 			   = 0.0;
    fNVertices 		           = 0.0;
    fTauGenPt 		           = 0.0;
    fTauGenEta 		           = 0.0;
    fTauGenPhi 		           = 0.0;
    fTauPt  		           = 0.0;
    fTauEta 		           = 0.0;
    fTauPhi 		           = 0.0;
    fPdgId		           = 13;
    fDRToClosestParton             = 9999;
    fTauIsTight                    = false;
    fTauIsLoose                    = false;
    fPassVetoSelection             = false;
    fPassLooseSelection            = false;
    fPassTightSelection            = false;
    fTauTypeBits	           = 0.0;
    fIsAllArbitrated	           = 0.0;
    fTauTkNchi2 	           = 0.0;
    fTauGlobalNchi2 	           = 0.0;
    fTauNValidHits  	           = 0.0;
    fTauNTrackerHits 	           = 0.0;
    fTauNPixelHits  	           = 0.0;
    fTauNMatches 	           = 0.0;
    fTauD0  		           = 0.0;
    fTauIP3d 		           = 0.0;
    fTauIP3dSig 	           = 0.0;
    fTauTrkKink 	           = 0.0;
    fTauGlobalKink  	           = 0.0;
    fTauSegmentCompatibility       = 0.0;
    fTauCaloCompatibility 	   = 0.0;
    fTauHadEnergy 		   = 0.0;
    fTauHoEnergy 		   = 0.0;
    fTauEmEnergy 		   = 0.0;
    fTauHadS9Energy 		   = 0.0;
    fTauHoS9Energy  		   = 0.0;
    fTauEmS9Energy  		   = 0.0;
    fTauTrkIso03 		   = 0.0;
    fTauEMIso03 		   = 0.0;
    fTauHadIso03 		   = 0.0;
    fTauTrkIso05 		   = 0.0;
    fTauEMIso05 		   = 0.0;
    fTauHadIso05 		   = 0.0;
    fTauPFIso04 		   = 0.0;
    fChargedIso_DR0p0To0p1 	   = 0.0;
    fChargedIso_DR0p1To0p2 	   = 0.0;
    fChargedIso_DR0p2To0p3 	   = 0.0;
    fChargedIso_DR0p3To0p4 	   = 0.0;
    fChargedIso_DR0p4To0p5 	   = 0.0;
    fGammaIso_DR0p0To0p1	   = 0.0;
    fGammaIso_DR0p1To0p2	   = 0.0;
    fGammaIso_DR0p2To0p3	   = 0.0;
    fGammaIso_DR0p3To0p4	   = 0.0;
    fGammaIso_DR0p4To0p5	   = 0.0;
    fNeutralHadronIso_DR0p0To0p1   = 0.0;
    fNeutralHadronIso_DR0p1To0p2   = 0.0;
    fNeutralHadronIso_DR0p2To0p3   = 0.0;
    fNeutralHadronIso_DR0p3To0p4   = 0.0;
    fNeutralHadronIso_DR0p4To0p5   = 0.0;
    fTauPassTriggerDenominator     = kFALSE;
  }
  
  /// load a TauTree
  void LoadTree(const char* file){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->Get("Taus"));
    assert(tree_);
  }
  
  /// create a TauTree
  void CreateTree(TauTreeVersion version = kTauTreeStd){
    tree_ = new TTree("Taus","Taus");
    f_ = 0;
    
    //book the branches
    tree_->Branch("weight",&fWeight,"weight/F");
    tree_->Branch("run",&fRunNumber,"run/i");
    tree_->Branch("lumi",&fLumiSectionNumber,"lumi/i");
    tree_->Branch("event",&fEventNumber,"event/i");
    tree_->Branch("EventNumberParity",&fTauEventNumberParity,"EventNumberParity/O"); 
    tree_->Branch("Rho",&fRho,"Rho/F"); 
    tree_->Branch("NVertices",&fNVertices,"NVertices/F"); 
    tree_->Branch("pt",&fTauPt,"pt/F"); 
    tree_->Branch("eta",&fTauEta,"eta/F"); 
    tree_->Branch("phi",&fTauPhi,"phi/F"); 
    tree_->Branch("genpt",&fTauGenPt,"genpt/F"); 
    tree_->Branch("geneta",&fTauGenEta,"geneta/F"); 
    tree_->Branch("genphi",&fTauGenPhi,"genphi/F"); 
    tree_->Branch("pdgid",&fPdgId,"pdgid/I"); 
    tree_->Branch("DRToClosestParton",&fDRToClosestParton,"DRToClosestParton/F"); 
    tree_->Branch("IsTight",&fTauIsTight,"IsTight/O"); 
    tree_->Branch("IsLoose",&fTauIsLoose,"IsLoose/O"); 
    tree_->Branch("PassVetoSelection",&fPassVetoSelection,"PassVetoSelection/O"); 
    tree_->Branch("PassLooseSelection",&fPassLooseSelection,"PassLooseSelection/O"); 	
    tree_->Branch("PassTightSelection",&fPassTightSelection,"PassTightSelection/O"); 
    tree_->Branch("PFIso04",&fTauPFIso04,"PFIso04/F"); 
    tree_->Branch("D0",&fTauD0,"D0/F"); 
    tree_->Branch("IP3d",&fTauIP3d,"IP3d/F"); 
    tree_->Branch("IP3dSig",&fTauIP3dSig,"IP3dSig/F"); 
    if (version == kTauTreeStd ) {
      tree_->Branch("typeBits",&fTauTypeBits,"typeBits/i"); 
      tree_->Branch("isAllArbitrated",&fIsAllArbitrated,"isAllArbitrated/i"); 
      tree_->Branch("TkNchi2",&fTauTkNchi2,"TkNchi2/F"); 
      tree_->Branch("GlobalNchi2",&fTauGlobalNchi2,"GlobalNchi2/F"); 
      tree_->Branch("NValidHits",&fTauNValidHits,"NValidHits/F"); 
      tree_->Branch("NTrackerHits",&fTauNTrackerHits,"NTrackerHits/F"); 
      tree_->Branch("NPixelHits",&fTauNPixelHits,"NPixelHits/F"); 
      tree_->Branch("NMatches",&fTauNMatches,"NMatches/F"); 
      tree_->Branch("TrkKink",&fTauTrkKink,"TrkKink/F"); 
      tree_->Branch("GlobalKink",&fTauGlobalKink,"GlobalKink/F"); 
      tree_->Branch("SegmentCompatibility",&fTauSegmentCompatibility,"SegmentCompatibility/F"); 
      tree_->Branch("CaloCompatibility",&fTauCaloCompatibility,"CaloCompatibility/F"); 
      tree_->Branch("HadEnergy",&fTauHadEnergy,"HadEnergy/F"); 
      tree_->Branch("HoEnergy",&fTauHoEnergy,"HoEnergy/F"); 
      tree_->Branch("EmEnergy",&fTauEmEnergy,"EmEnergy/F"); 
      tree_->Branch("HadS9Energy",&fTauHadS9Energy,"HadS9Energy/F"); 
      tree_->Branch("HoS9Energy",&fTauHoS9Energy,"HoS9Energy/F"); 
      tree_->Branch("EmS9Energy",&fTauEmS9Energy,"EmS9Energy/F"); 
      tree_->Branch("TrkIso03",&fTauTrkIso03,"TrkIso03/F"); 
      tree_->Branch("EMIso03",&fTauEMIso03,"EMIso03/F"); 
      tree_->Branch("HadIso03",&fTauHadIso03,"HadIso03/F"); 
      tree_->Branch("TrkIso05",&fTauTrkIso05,"TrkIso05/F"); 
      tree_->Branch("EMIso05",&fTauEMIso05,"EMIso05/F"); 
      tree_->Branch("HadIso05",&fTauHadIso05,"HadIso05/F");        
      tree_->Branch("ChargedIso_DR0p0To0p1",&fChargedIso_DR0p0To0p1,"ChargedIso_DR0p0To0p1/F");
      tree_->Branch("ChargedIso_DR0p1To0p2",&fChargedIso_DR0p1To0p2,"ChargedIso_DR0p1To0p2/F");
      tree_->Branch("ChargedIso_DR0p2To0p3",&fChargedIso_DR0p2To0p3,"ChargedIso_DR0p2To0p3/F");
      tree_->Branch("ChargedIso_DR0p3To0p4",&fChargedIso_DR0p3To0p4,"ChargedIso_DR0p3To0p4/F");
      tree_->Branch("ChargedIso_DR0p4To0p5",&fChargedIso_DR0p4To0p5,"ChargedIso_DR0p4To0p5/F");
      tree_->Branch("GammaIso_DR0p0To0p1",&fGammaIso_DR0p0To0p1,"GammaIso_DR0p0To0p1/F");
      tree_->Branch("GammaIso_DR0p1To0p2",&fGammaIso_DR0p1To0p2,"GammaIso_DR0p1To0p2/F");
      tree_->Branch("GammaIso_DR0p2To0p3",&fGammaIso_DR0p2To0p3,"GammaIso_DR0p2To0p3/F");
      tree_->Branch("GammaIso_DR0p3To0p4",&fGammaIso_DR0p3To0p4,"GammaIso_DR0p3To0p4/F");
      tree_->Branch("GammaIso_DR0p4To0p5",&fGammaIso_DR0p4To0p5,"GammaIso_DR0p4To0p5/F");
      tree_->Branch("NeutralHadronIso_DR0p0To0p1",&fNeutralHadronIso_DR0p0To0p1,"NeutralHadronIso_DR0p0To0p1/F");
      tree_->Branch("NeutralHadronIso_DR0p1To0p2",&fNeutralHadronIso_DR0p1To0p2,"NeutralHadronIso_DR0p1To0p2/F");
      tree_->Branch("NeutralHadronIso_DR0p2To0p3",&fNeutralHadronIso_DR0p2To0p3,"NeutralHadronIso_DR0p2To0p3/F");
      tree_->Branch("NeutralHadronIso_DR0p3To0p4",&fNeutralHadronIso_DR0p3To0p4,"NeutralHadronIso_DR0p3To0p4/F");
      tree_->Branch("NeutralHadronIso_DR0p4To0p5",&fNeutralHadronIso_DR0p4To0p5,"NeutralHadronIso_DR0p4To0p5/F");
      tree_->Branch("PassTriggerDenominator",&fTauPassTriggerDenominator,"PassTriggerDenominator/O");
    }
  } 
  
  // initialze a TauTree
  void InitTree(TauTreeVersion version = kTauTreeStd){
    assert(tree_);
    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;
    // gErrorIgnoreLevel = kError;
    gErrorIgnoreLevel = kBreak;
    
    tree_->SetBranchAddress("weight",&fWeight);
    tree_->SetBranchAddress("run",&fRunNumber);
    tree_->SetBranchAddress("lumi",&fLumiSectionNumber);
    tree_->SetBranchAddress("event",&fEventNumber);
    tree_->SetBranchAddress("EventNumberParity",&fTauEventNumberParity);
    tree_->SetBranchAddress("Rho",&fRho);
    tree_->SetBranchAddress("NVertices",&fNVertices);
    tree_->SetBranchAddress("pt",&fTauPt);
    tree_->SetBranchAddress("eta",&fTauEta);
    tree_->SetBranchAddress("phi",&fTauPhi);
    tree_->SetBranchAddress("genpt",&fTauGenPt);
    tree_->SetBranchAddress("geneta",&fTauGenEta);
    tree_->SetBranchAddress("genphi",&fTauGenPhi);
    tree_->SetBranchAddress("pdgid",&fPdgId);
    tree_->SetBranchAddress("DRToClosestParton",&fDRToClosestParton);
    tree_->SetBranchAddress("IsTight",&fTauIsTight); 
    tree_->SetBranchAddress("IsLoose",&fTauIsLoose); 
    tree_->SetBranchAddress("PassVetoSelection",&fPassVetoSelection); 
    tree_->SetBranchAddress("PassLooseSelection",&fPassLooseSelection); 
    tree_->SetBranchAddress("PassTightSelection",&fPassTightSelection); 
    tree_->SetBranchAddress("PFIso04",&fTauPFIso04);
    tree_->SetBranchAddress("D0",&fTauD0);
    tree_->SetBranchAddress("IP3d",&fTauIP3d);
    tree_->SetBranchAddress("IP3dSig",&fTauIP3dSig);
    
    if (version == kTauTreeStd ) {
      tree_->SetBranchAddress("typeBits",&fTauTypeBits);
      tree_->SetBranchAddress("isAllArbitrated",&fIsAllArbitrated);
      tree_->SetBranchAddress("TkNchi2",&fTauTkNchi2);
      tree_->SetBranchAddress("GlobalNchi2",&fTauGlobalNchi2);
      tree_->SetBranchAddress("NValidHits",&fTauNValidHits);
      tree_->SetBranchAddress("NTrackerHits",&fTauNTrackerHits);
      tree_->SetBranchAddress("NPixelHits",&fTauNPixelHits);
      tree_->SetBranchAddress("NMatches",&fTauNMatches);
      tree_->SetBranchAddress("TrkKink",&fTauTrkKink);
      tree_->SetBranchAddress("GlobalKink",&fTauGlobalKink);
      tree_->SetBranchAddress("SegmentCompatibility",&fTauSegmentCompatibility);
      tree_->SetBranchAddress("CaloCompatibility",&fTauCaloCompatibility);
      tree_->SetBranchAddress("HadEnergy",&fTauHadEnergy);
      tree_->SetBranchAddress("HoEnergy",&fTauHoEnergy);
      tree_->SetBranchAddress("EmEnergy",&fTauEmEnergy);
      tree_->SetBranchAddress("HadS9Energy",&fTauHadS9Energy);
      tree_->SetBranchAddress("HoS9Energy",&fTauHoS9Energy);
      tree_->SetBranchAddress("EmS9Energy",&fTauEmS9Energy);
      tree_->SetBranchAddress("TrkIso03",&fTauTrkIso03);
      tree_->SetBranchAddress("EMIso03",&fTauEMIso03);
      tree_->SetBranchAddress("HadIso03",&fTauHadIso03);
      tree_->SetBranchAddress("TrkIso05",&fTauTrkIso05);
      tree_->SetBranchAddress("EMIso05",&fTauEMIso05);
      tree_->SetBranchAddress("HadIso05",&fTauHadIso05);       
      tree_->SetBranchAddress("ChargedIso_DR0p0To0p1",&fChargedIso_DR0p0To0p1);
      tree_->SetBranchAddress("ChargedIso_DR0p1To0p2",&fChargedIso_DR0p1To0p2);
      tree_->SetBranchAddress("ChargedIso_DR0p2To0p3",&fChargedIso_DR0p2To0p3);
      tree_->SetBranchAddress("ChargedIso_DR0p3To0p4",&fChargedIso_DR0p3To0p4);
      tree_->SetBranchAddress("ChargedIso_DR0p4To0p5",&fChargedIso_DR0p4To0p5);
      tree_->SetBranchAddress("GammaIso_DR0p0To0p1",&fGammaIso_DR0p0To0p1);
      tree_->SetBranchAddress("GammaIso_DR0p1To0p2",&fGammaIso_DR0p1To0p2);
      tree_->SetBranchAddress("GammaIso_DR0p2To0p3",&fGammaIso_DR0p2To0p3);
      tree_->SetBranchAddress("GammaIso_DR0p3To0p4",&fGammaIso_DR0p3To0p4);
      tree_->SetBranchAddress("GammaIso_DR0p4To0p5",&fGammaIso_DR0p4To0p5);
      tree_->SetBranchAddress("NeutralHadronIso_DR0p0To0p1",&fNeutralHadronIso_DR0p0To0p1);
      tree_->SetBranchAddress("NeutralHadronIso_DR0p1To0p2",&fNeutralHadronIso_DR0p1To0p2);
      tree_->SetBranchAddress("NeutralHadronIso_DR0p2To0p3",&fNeutralHadronIso_DR0p2To0p3);
      tree_->SetBranchAddress("NeutralHadronIso_DR0p3To0p4",&fNeutralHadronIso_DR0p3To0p4);
      tree_->SetBranchAddress("NeutralHadronIso_DR0p4To0p5",&fNeutralHadronIso_DR0p4To0p5);
      tree_->SetBranchAddress("PassTriggerDenominator",&fTauPassTriggerDenominator);
    }
    gErrorIgnoreLevel = currentState;
  }
  
}; 


#endif
