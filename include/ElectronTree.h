#ifndef ElectronTree_H
#define ElectronTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

class ElectronTree {

 public:

  /// bit map
  /// DON'T CHANGE ORDER
  //*******************************************
  //=== ElectronTriggerBits  ====
  //*******************************************
  enum ElectronTriggerBits { kEleTrigger_Ele27Loose                             = 0x000001, 
			     kEleTrigger_Ele27Tight                             = 0x000002, 
			     kEleTrigger_Ele32Tight                             = 0x000004,
			     kEleTrigger_Ele105                                 = 0x000008,
			     kEleTrigger_Ele115                                 = 0x000010,
  };

  //******************************************* 
  //=== tree versions ===
  //*******************************************
  enum ElectronTreeVersion { kEleTreeStd, kEleTreeLight };

  /// variables
  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Bool_t                  fEleEventNumberParity;
  Float_t                 fEleGenPt; 
  Float_t                 fEleGenEta; 
  Float_t                 fEleGenPhi;
  Float_t                 fElePt; 
  Float_t                 fEleEta; 
  Float_t                 fElePhi; 
  Float_t                 fEleSCEt; 
  Float_t                 fEleSCEta; 
  Float_t                 fEleSCPhi; 
  Float_t                 fEleEcalEnergy; 
  Bool_t                  fEleIsEcalDriven;
  UInt_t                  fEleTriggerBit;
  UInt_t                  fNPU;
  Float_t                 fRho; 
  Float_t                 fRhoNeutralCentral; 
  Float_t                 fNVertices; 
  Int_t                   fPdgId;
  Float_t                 fDRToClosestParton;
  Float_t                 fActivity;

  // Typical Selection Working Points
  Bool_t                  fPassVetoSelection;
  Bool_t                  fPassLooseSelection;
  Bool_t                  fPassTightSelection;
  Bool_t                  fPassMVAVetoSelection;

  // Conversion and IP
  Float_t                 fEleD0; 
  Float_t                 fEleDZ; 
  Float_t                 fEleIP3d; 
  Float_t                 fEleIP3dSig; 
  Bool_t                  fElePassConversion;
  Float_t                 fEleConvDCot;
  Float_t                 fEleConvDist;
  Float_t                 fEleNMissHits;

  // E/P variables
  Float_t                 fEleNBrem; 
  Float_t                 fEleFBrem; 
  Float_t                 fEleEOverP; 
  Float_t                 fEleESeedClusterOverPIn; 
  Float_t                 fEleESeedClusterOverPout; 
  Float_t                 fEleEEleClusterOverPout; 
  Float_t                 fEleOneOverEMinusOneOverP; 

  // track cluster matching
  Float_t                 fEleDEtaIn; 
  Float_t                 fEleDPhiIn; 
  Float_t                 fEledEtaCalo;
  Float_t                 fEledPhiCalo;

  //shower shape
  Float_t                 fEleSigmaIEtaIEta; 
  Float_t                 fEleSigmaIPhiIPhi; 
  Float_t                 fEleSigmaIEtaIPhi;
  Float_t                 fEleSCEtaWidth;
  Float_t                 fEleSCPhiWidth;
  Float_t                 fEleR9;
  Float_t                 fElePreShowerOverRaw;
  Float_t                 fEleHoverE; 

  //track quality
  Float_t                 fEleGsfTrackChi2OverNdof;    
  Float_t                 fEleKFTrackChi2OverNDoF;
  Float_t                 fEleKFTrackNHits;
  Float_t                 fEleKFTrackNLayersWithMeasurement;
  Float_t                 fEleOneMinusSeedE1x5OverE5x5;

  //Isolation Variables
  Float_t                 fElePFMVA;
  Float_t                 fEleTrkIso03; 
  Float_t                 fEleEMIso03; 
  Float_t                 fEleHadIso03; 
  Float_t                 fEleTrkIso04; 
  Float_t                 fEleEMIso04; 
  Float_t                 fEleHadIso04;     
  Float_t                 fElePFIso04; 
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
  Float_t                 fPtRel;
  Float_t                 fMiniIsoCharged;
  Float_t                 fMiniIsoNeutral;
  Float_t                 fMiniIso;
  Float_t                 fMiniIsoDBCorr;

  Bool_t                  fElePassTriggerDenominator;

  //Regression Variables
  Bool_t                  fIsEB;           
  Bool_t                  fIsEE;           
  Float_t                 fSCRawEnergy;
  Float_t                 fNClusters;
  Float_t                 fEtaSeed;
  Float_t                 fPhiSeed;
  Float_t                 fESeed;
  Float_t                 fE3x3Seed;
  Float_t                 fE5x5Seed;
  Float_t                 fEMaxSeed;
  Float_t                 fE2ndSeed;
  Float_t                 fETopSeed;
  Float_t                 fEBottomSeed;
  Float_t                 fELeftSeed;
  Float_t                 fERightSeed;
  Float_t                 fE2x5MaxSeed;
  Float_t                 fE2x5TopSeed;
  Float_t                 fE2x5BottomSeed;
  Float_t                 fE2x5LeftSeed;
  Float_t                 fE2x5RightSeed;
  Float_t                 fIEtaSeed;
  Float_t                 fIPhiSeed;
  Float_t                 fEtaCrySeed;
  Float_t                 fPhiCrySeed;
  Float_t                 fEcalEnergyError;
  Float_t                 fGsfTrackPIn;
  Float_t                 fTrackMomentumError;
  Float_t                 fpmeangsf;
  Float_t                 fpmeankf;
  Float_t                 fCharge;
  Float_t                 fGeneratedEnergy;
  Float_t                 fGeneratedEnergyStatus1;
  Float_t                 fGeneratedEnergyStatus3;
  Float_t                 fEleMomentum_Regression_V0;
  Float_t                 fEleMomentum_Regression_V1;
  Float_t                 fEleMomentum_Regression_V2;
  Float_t                 fEleMomentumError_Regression_V0;
  Float_t                 fEleMomentumError_Regression_V1;
  Float_t                 fEleMomentumError_Regression_V2;
  Float_t                 fEleClassification;

  //MVA Variables
  Float_t                 fIDMVAHZZ;
  Float_t                 fIDMVAGeneralPurpose;

 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  
  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  /// default constructor  
  ElectronTree()  {};
  /// default destructor
  ~ElectronTree(){ 
    if (f_) f_->Close();  
  };
    
  /// initialize varibles and fill list of available variables
  void InitVariables() {
    fWeight			       = 0.0;
    fRunNumber		       = 0.0;
    fLumiSectionNumber	       = 0.0;
    fEventNumber		       = 0.0;
    fEleEventNumberParity 	       = 0.0;
    fEleGenPt        	       = 0.0;
    fEleGenEta 		       = 0.0;
    fEleGenPhi 		       = 0.0;
    fElePt 			       = 0.0;
    fEleEta 		       = 0.0;
    fElePhi 		       = 0.0;
    fEleSCEt 		       = 0.0;
    fEleSCEta 		       = 0.0;
    fEleSCPhi 		       = 0.0;
    fEleEcalEnergy 		       = 0.0;
    fEleIsEcalDriven	       = 0.0;
    fEleTriggerBit	       = 0.0;
    fNPU                       = 0;
    fRho  		       = 0.0;
    fRhoNeutralCentral         = 0.0;
    fNVertices 		       = 0.0;
    fPdgId    		       = 11;
    fDRToClosestParton         = 9999;
    fActivity                  = 9999;
    fPassVetoSelection                 = false;
    fPassLooseSelection                = false;
    fPassTightSelection                = false;
    fPassMVAVetoSelection       = false;
    fEleD0 			       = 0.0;
    fEleDZ 			       = 0.0;
    fEleIP3d 		       = 0.0;
    fEleIP3dSig 		       = 0.0;
    fElePassConversion 	       = 0.0;
    fEleConvDCot		       = 0.0;
    fEleConvDist		       = 0.0;
    fEleNMissHits 		       = 0.0;
    fEleNBrem 		       = 0.0;
    fEleFBrem 		       = 0.0;
    fEleEOverP 		       = 0.0;
    fEleESeedClusterOverPIn        = 0.0;
    fEleESeedClusterOverPout       = 0.0;
    fEleEEleClusterOverPout        = 0.0;
    fEleOneOverEMinusOneOverP      = 0.0;
    fEleDEtaIn 		       = 0.0;
    fEleDPhiIn 		       = 0.0;
    fEledEtaCalo		       = 0.0;
    fEledPhiCalo		       = 0.0;
    fEleSigmaIEtaIEta 	       = 0.0;
    fEleSigmaIPhiIPhi 	       = 0.0;
    fEleSigmaIEtaIPhi	       = 0.0;
    fEleSCEtaWidth		       = 0.0;
    fEleSCPhiWidth		       = 0.0;
    fEleR9			       = 0.0;
    fElePreShowerOverRaw	       = 0.0;
    fEleHoverE 		       = 0.0;
    fEleGsfTrackChi2OverNdof       = 0.0;
    fEleKFTrackChi2OverNDoF	       = 0.0;
    fEleKFTrackNHits	       = 0.0;
    fEleKFTrackNLayersWithMeasurement = 0.0;
    fEleOneMinusSeedE1x5OverE5x5   = 0.0;
    fElePFMVA		       = 0.0;
    fEleTrkIso03  		       = 0.0;
    fEleEMIso03 		       = 0.0;
    fEleHadIso03  		       = 0.0;
    fEleTrkIso04  		       = 0.0;
    fEleEMIso04 		       = 0.0;
    fEleHadIso04     	       = 0.0;
    fElePFIso04 		       = 0.0;
    fChargedIso_DR0p0To0p1	       = 0.0;
    fChargedIso_DR0p1To0p2	       = 0.0;
    fChargedIso_DR0p2To0p3	       = 0.0;
    fChargedIso_DR0p3To0p4	       = 0.0;
    fChargedIso_DR0p4To0p5	       = 0.0;
    fGammaIso_DR0p0To0p1	       = 0.0;
    fGammaIso_DR0p1To0p2	       = 0.0;
    fGammaIso_DR0p2To0p3	       = 0.0;
    fGammaIso_DR0p3To0p4	       = 0.0;
    fGammaIso_DR0p4To0p5	       = 0.0;
    fNeutralHadronIso_DR0p0To0p1   = 0.0;
    fNeutralHadronIso_DR0p1To0p2   = 0.0;
    fNeutralHadronIso_DR0p2To0p3   = 0.0;
    fNeutralHadronIso_DR0p3To0p4   = 0.0;
    fNeutralHadronIso_DR0p4To0p5   = 0.0;
    fPtRel                         = 0.0;
    fMiniIsoCharged                = 0.0;
    fMiniIsoNeutral                = 0.0;
    fMiniIso                       = 0.0;
    fMiniIsoDBCorr                 = 0.0;
    fElePassTriggerDenominator     = 0.0;
    fIsEB                          = 0.0;           
    fIsEE                          = 0.0;           
    fSCRawEnergy                   = 0.0;
    fNClusters                     = 0.0;
    fEtaSeed                       = 0.0;
    fPhiSeed                       = 0.0;
    fESeed                         = 0.0;
    fE3x3Seed                      = 0.0;
    fE5x5Seed                      = 0.0;
    fEMaxSeed                      = 0.0;
    fE2ndSeed                      = 0.0;
    fETopSeed                      = 0.0;
    fEBottomSeed                   = 0.0;
    fELeftSeed                     = 0.0;
    fERightSeed                    = 0.0;
    fE2x5MaxSeed                   = 0.0;
    fE2x5TopSeed                   = 0.0;
    fE2x5BottomSeed                = 0.0;
    fE2x5LeftSeed                  = 0.0;
    fE2x5RightSeed                 = 0.0;
    fIEtaSeed                      = 0.0;
    fIPhiSeed                      = 0.0;
    fEtaCrySeed                    = 0.0;
    fPhiCrySeed                    = 0.0;
    fEcalEnergyError               = 0.0;
    fGsfTrackPIn                   = 0.0;
    fTrackMomentumError            = 0.0;
    fpmeangsf                      = 0.0;
    fpmeankf                       = 0.0;
    fCharge                        = 0.0;
    fGeneratedEnergy               = 0.0;
    fGeneratedEnergyStatus1        = 0.0;
    fGeneratedEnergyStatus3        = 0.0;
    fEleMomentum_Regression_V0     = 0.0;
    fEleMomentum_Regression_V1     = 0.0;
    fEleMomentum_Regression_V2     = 0.0;
    fEleMomentumError_Regression_V0= 0.0;
    fEleMomentumError_Regression_V1= 0.0;
    fEleMomentumError_Regression_V2= 0.0;
    fEleClassification             = 0.0;
    fIDMVAHZZ                     = -9999.0;
    fIDMVAGeneralPurpose                  = -9999.0;

  }
    
  /// load a ElectronTree
  void LoadTree(const char* file){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->Get("Electrons"));
    assert(tree_);
  }
    
  /// create a ElectronTree
  void CreateTree(int version=kEleTreeStd){
    tree_ = new TTree("Electrons","Electrons");
    f_ = 0;

    //book the branches
    tree_->Branch("weight",&fWeight,"weight/F");
    tree_->Branch("run",&fRunNumber,"run/i");
    tree_->Branch("lumi",&fLumiSectionNumber,"lumi/i");
    tree_->Branch("event",&fEventNumber,"event/i");
    tree_->Branch("EventNumberParity",&fEleEventNumberParity,"EventNumberParity/O"); 
    tree_->Branch("rho",&fRho,"rho/F"); 
    tree_->Branch("NPU",&fNPU,"NPU/i"); 
    tree_->Branch("rhoNeutralCentral",&fRhoNeutralCentral,"rhoNeutralCentral/F"); 
    tree_->Branch("pt",&fElePt,"pt/F"); 
    tree_->Branch("eta",&fEleEta,"eta/F"); 
    tree_->Branch("phi",&fElePhi,"phi/F");
    tree_->Branch("Charge",&fCharge,"Charge/F");
    tree_->Branch("genpt",&fEleGenPt,"genpt/F"); 
    tree_->Branch("geneta",&fEleGenEta,"geneta/F"); 
    tree_->Branch("genphi",&fEleGenPhi,"genphi/F");
    tree_->Branch("pdgid",&fPdgId,"pdgid/I");
    tree_->Branch("DRToClosestParton",&fDRToClosestParton,"DRToClosestParton/F");
    tree_->Branch("Activity",&fActivity,"Activity/F");
    tree_->Branch("scEta",&fEleSCEta,"scEta/F"); 
    tree_->Branch("d0",&fEleD0,"d0/F"); 
    tree_->Branch("dz",&fEleDZ,"dz/F"); 
    tree_->Branch("missHits",&fEleNMissHits,"NMissHits/F"); 
    tree_->Branch("passConv",&fElePassConversion,"passConv/O"); 
    tree_->Branch("deta",&fEleDEtaIn,"deta/F"); 
    tree_->Branch("dphi",&fEleDPhiIn,"dphi/F"); 
    tree_->Branch("see",&fEleSigmaIEtaIEta,"see/F"); 
    tree_->Branch("IoEmIoP",&fEleOneOverEMinusOneOverP,"IoEmIoP/F");  
    tree_->Branch("R9",&fEleR9,"R9/F"); 
    tree_->Branch("HoE",&fEleHoverE,"HoE/F"); 
    tree_->Branch("pfIso04",&fElePFIso04,"pfIso04/F"); 
    tree_->Branch("IDMVAHZZ",&fIDMVAHZZ,"IDMVAHZZ/F");
    tree_->Branch("IDMVAGeneralPurpose",&fIDMVAGeneralPurpose,"IDMVAGeneralPurpose/F");
    tree_->Branch("vertices",&fNVertices,"vertices/F"); 
    tree_->Branch("PassVetoSelection",&fPassVetoSelection,"PassVetoSelection/O"); 
    tree_->Branch("PassLooseSelection",&fPassLooseSelection,"PassLooseSelection/O"); 
    tree_->Branch("PassTightSelection",&fPassTightSelection,"PassTightSelection/O"); 
    tree_->Branch("PassMVAVetoSelection",&fPassMVAVetoSelection,"PassMVAVetoSelection/O"); 
    tree_->Branch("PtRel",&fPtRel,"PtRel/F"); 
    tree_->Branch("MiniIsoCharged",&fMiniIsoCharged,"MiniIsoCharged/F"); 
    tree_->Branch("MiniIsoNeutral",&fMiniIsoNeutral,"MiniIsoNeutral/F"); 
    tree_->Branch("MiniIso",&fMiniIso,"MiniIso/F"); 
    tree_->Branch("MiniIsoDBCorr",&fMiniIsoDBCorr,"MiniIsoDBCorr/F"); 
    tree_->Branch("triggerBit",&fEleTriggerBit,"triggerBit/i"); 
    tree_->Branch("ip3ds",&fEleIP3dSig,"ip3ds/F"); 


    if (version == kEleTreeStd ) {
      tree_->Branch("scEt",&fEleSCEt,"scEt/F"); 
      tree_->Branch("scPhi",&fEleSCPhi,"scPhi/F"); 
      tree_->Branch("ecalenergy",&fEleEcalEnergy,"ecalenergy/F"); 
      tree_->Branch("ecaldriven",&fEleIsEcalDriven,"ecaldriven/O"); 
      tree_->Branch("ip3d",&fEleIP3d,"ip3d/F"); 
      tree_->Branch("dcot",&fEleConvDCot,"dcot/F"); 
      tree_->Branch("dist",&fEleConvDist,"dist/F"); 
      tree_->Branch("nbrems",&fEleNBrem,"nbrems/F"); 
      tree_->Branch("fbrem",&fEleFBrem,"fbrem/F"); 
      tree_->Branch("EoP",&fEleEOverP,"EoP/F"); 
      tree_->Branch("EoPin",&fEleESeedClusterOverPIn,"EoPin/F"); 
      tree_->Branch("ESeedoPout",&fEleESeedClusterOverPout,"ESeedoPout/F"); 
      tree_->Branch("EEleoPout",&fEleEEleClusterOverPout,"EEleoPout/F"); 
      tree_->Branch("detacalo",&fEledEtaCalo,"detacalo/F"); 
      tree_->Branch("dphicalo",&fEledPhiCalo,"dphicalo/F"); 
      tree_->Branch("spp",&fEleSigmaIPhiIPhi,"spp/F"); 
      tree_->Branch("sep",&fEleSigmaIEtaIPhi,"sep/F"); 
      tree_->Branch("etawidth",&fEleSCEtaWidth,"etawidth/F"); 
      tree_->Branch("phiwidth",&fEleSCPhiWidth,"phiwidth/F"); 
      tree_->Branch("PreShowerOverRaw",&fElePreShowerOverRaw,"PreShowerOverRaw/F"); 
      tree_->Branch("gsfchi2",&fEleGsfTrackChi2OverNdof,"gsfchi2/F"); 
      tree_->Branch("kfchi2",&fEleKFTrackChi2OverNDoF,"kfchi2/F"); 
      tree_->Branch("kfhits",&fEleKFTrackNHits,"kfhits/F"); 
      tree_->Branch("kflayers",&fEleKFTrackNLayersWithMeasurement,"kflayers/F"); 
      tree_->Branch("OneMinusSeedE1x5OverE5x5",&fEleOneMinusSeedE1x5OverE5x5,"OneMinusSeedE1x5OverE5x5/F"); 
      tree_->Branch("PFMVA",&fElePFMVA,"PFMVA/F"); 
      tree_->Branch("trkIso03",&fEleTrkIso03,"trkIso03/F"); 
      tree_->Branch("ecalIso03",&fEleEMIso03,"ecalIso03/F"); 
      tree_->Branch("hcalIso03",&fEleHadIso03,"hcalIso03/F"); 
      tree_->Branch("trkIso04",&fEleTrkIso04,"trkIso04/F"); 
      tree_->Branch("ecalIso04",&fEleEMIso04,"ecalIso04/F"); 
      tree_->Branch("hcalIso04",&fEleHadIso04,"hcalIso04/F"); 
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
      tree_->Branch("PassTriggerDenominator",&fElePassTriggerDenominator,"PassTriggerDenominator/O"); 
      tree_->Branch("IsEB",&fIsEB,"IsEB/O");
      tree_->Branch("IsEE",&fIsEE,"IsEE/O");
      tree_->Branch("SCRawEnergy",&fSCRawEnergy,"SCRawEnergy/F");
      tree_->Branch("NClusters",&fNClusters,"NClusters/F");
      tree_->Branch("EtaSeed",&fEtaSeed,"EtaSeed/F");
      tree_->Branch("PhiSeed",&fPhiSeed,"PhiSeed/F");
      tree_->Branch("ESeed",&fESeed,"ESeed/F");
      tree_->Branch("E3x3Seed",&fE3x3Seed,"E3x3Seed/F");
      tree_->Branch("E5x5Seed",&fE5x5Seed,"E5x5Seed/F");
      tree_->Branch("EMaxSeed",&fEMaxSeed,"EMaxSeed/F");
      tree_->Branch("E2ndSeed",&fE2ndSeed,"E2ndSeed/F");
      tree_->Branch("ETopSeed",&fETopSeed,"ETopSeed/F");
      tree_->Branch("EBottomSeed",&fEBottomSeed,"EBottomSeed/F");
      tree_->Branch("ELeftSeed",&fELeftSeed,"ELeftSeed/F");
      tree_->Branch("ERightSeed",&fERightSeed,"ERightSeed/F");
      tree_->Branch("E2x5MaxSeed",&fE2x5MaxSeed,"E2x5MaxSeed/F");
      tree_->Branch("E2x5TopSeed",&fE2x5TopSeed,"E2x5TopSeed/F");
      tree_->Branch("E2x5BottomSeed",&fE2x5BottomSeed,"E2x5BottomSeed/F");
      tree_->Branch("E2x5LeftSeed",&fE2x5LeftSeed,"E2x5LeftSeed/F");
      tree_->Branch("E2x5RightSeed",&fE2x5RightSeed,"E2x5RightSeed/F");
      tree_->Branch("IEtaSeed",&fIEtaSeed,"IEtaSeed/F");
      tree_->Branch("IPhiSeed",&fIPhiSeed,"IPhiSeed/F");
      tree_->Branch("EtaCrySeed",&fEtaCrySeed,"EtaCrySeed/F");
      tree_->Branch("PhiCrySeed",&fPhiCrySeed,"PhiCrySeed/F");
      tree_->Branch("ecalenergyerror",&fEcalEnergyError,"ecalenergyerror/F");
      tree_->Branch("pmodegsf",&fGsfTrackPIn,"pmodegsf/F");
      tree_->Branch("perror",&fTrackMomentumError,"perror/F");
      tree_->Branch("pmeangsf",&fpmeangsf,"pmodegsf/F");
      tree_->Branch("pmeankf",&fpmeankf,"pmeankf/F");
      tree_->Branch("GeneratedEnergy",&fGeneratedEnergy,"GeneratedEnergy/F");
      tree_->Branch("GeneratedEnergyStatus1",&fGeneratedEnergyStatus1,"GeneratedEnergyStatus1/F");
      tree_->Branch("GeneratedEnergyStatus3",&fGeneratedEnergyStatus3,"GeneratedEnergyStatus3/F");
      tree_->Branch("EleMomentum_Regression_V0",&fEleMomentum_Regression_V0,"EleMomentum_Regression_V0/F");
      tree_->Branch("EleMomentum_Regression_V1",&fEleMomentum_Regression_V1,"EleMomentum_Regression_V1/F");
      tree_->Branch("EleMomentum_Regression_V2",&fEleMomentum_Regression_V2,"EleMomentum_Regression_V2/F");
      tree_->Branch("EleMomentumError_Regression_V0",&fEleMomentumError_Regression_V0,"EleMomentumError_Regression_V0/F");
      tree_->Branch("EleMomentumError_Regression_V1",&fEleMomentumError_Regression_V1,"EleMomentumError_Regression_V1/F");
      tree_->Branch("EleMomentumError_Regression_V2",&fEleMomentumError_Regression_V2,"EleMomentumError_Regression_V2/F");
      tree_->Branch("EleClassification",&fEleClassification,"EleClassification/F");
    }

  } 

  // initialze a ElectronTree
  void InitTree(int version=kEleTreeStd){
    assert(tree_);
    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;

    tree_->SetBranchAddress("weight",&fWeight);
    tree_->SetBranchAddress("run",&fRunNumber);
    tree_->SetBranchAddress("lumi",&fLumiSectionNumber);
    tree_->SetBranchAddress("event",&fEventNumber);
    tree_->SetBranchAddress("EventNumberParity",&fEleEventNumberParity);
    tree_->SetBranchAddress("rho",&fRho);
    tree_->SetBranchAddress("NPU",&fNPU);
    tree_->SetBranchAddress("rhoNeutralCentral",&fRhoNeutralCentral);
    tree_->SetBranchAddress("pt",&fElePt);
    tree_->SetBranchAddress("eta",&fEleEta);
    tree_->SetBranchAddress("phi",&fElePhi);
    tree_->SetBranchAddress("genpt",&fEleGenPt);
    tree_->SetBranchAddress("geneta",&fEleGenEta);
    tree_->SetBranchAddress("genphi",&fEleGenPhi);
    tree_->SetBranchAddress("pdgid",&fPdgId);
    tree_->SetBranchAddress("DRToClosestParton",&fDRToClosestParton);
    tree_->SetBranchAddress("Activity",&fActivity);
    tree_->SetBranchAddress("Charge",&fCharge);
    tree_->SetBranchAddress("scEta",&fEleSCEta);
    tree_->SetBranchAddress("d0",&fEleD0);
    tree_->SetBranchAddress("dz",&fEleDZ);
    tree_->SetBranchAddress("vertices",&fNVertices);
    tree_->SetBranchAddress("passConv",&fElePassConversion);
    tree_->SetBranchAddress("missHits",&fEleNMissHits);
    tree_->SetBranchAddress("IoEmIoP",&fEleOneOverEMinusOneOverP);
    tree_->SetBranchAddress("deta",&fEleDEtaIn);
    tree_->SetBranchAddress("dphi",&fEleDPhiIn);
    tree_->SetBranchAddress("see",&fEleSigmaIEtaIEta);
    tree_->SetBranchAddress("R9",&fEleR9);
    tree_->SetBranchAddress("HoE",&fEleHoverE);
    tree_->SetBranchAddress("pfIso04",&fElePFIso04);
    tree_->SetBranchAddress("IDMVAHZZ",&fIDMVAHZZ);
    tree_->SetBranchAddress("IDMVAGeneralPurpose",&fIDMVAGeneralPurpose);
    tree_->SetBranchAddress("PassVetoSelection",&fPassVetoSelection);
    tree_->SetBranchAddress("PassLooseSelection",&fPassLooseSelection);
    tree_->SetBranchAddress("PassTightSelection",&fPassTightSelection);
    tree_->SetBranchAddress("PassMVAVetoSelection",&fPassMVAVetoSelection);
    //tree_->SetBranchAddress("PassMVANonTrigVetoSelection",&fPassMVAVetoSelection);
    tree_->SetBranchAddress("PtRel",&fPtRel);
    tree_->SetBranchAddress("MiniIsoCharged",&fMiniIsoCharged);
    tree_->SetBranchAddress("MiniIsoNeutral",&fMiniIsoNeutral);
    tree_->SetBranchAddress("MiniIso",&fMiniIso);
    tree_->SetBranchAddress("MiniIsoDBCorr",&fMiniIsoDBCorr);
    tree_->SetBranchAddress("triggerBit",&fEleTriggerBit);
    tree_->SetBranchAddress("ip3ds",&fEleIP3dSig);
    
    if (version == kEleTreeStd ) {
      tree_->SetBranchAddress("scEt",&fEleSCEt);
      tree_->SetBranchAddress("scPhi",&fEleSCPhi);
      tree_->SetBranchAddress("ecalenergy",&fEleEcalEnergy);
      tree_->SetBranchAddress("ecaldriven",&fEleIsEcalDriven);
      tree_->SetBranchAddress("ip3d",&fEleIP3d);
      tree_->SetBranchAddress("dcot",&fEleConvDCot);
      tree_->SetBranchAddress("dist",&fEleConvDist);
      tree_->SetBranchAddress("nbrems",&fEleNBrem);
      tree_->SetBranchAddress("fbrem",&fEleFBrem);
      tree_->SetBranchAddress("EoP",&fEleEOverP);
      tree_->SetBranchAddress("EoPin",&fEleESeedClusterOverPIn);
      tree_->SetBranchAddress("ESeedoPout",&fEleESeedClusterOverPout);
      tree_->SetBranchAddress("EEleoPout",&fEleEEleClusterOverPout);
      tree_->SetBranchAddress("detacalo",&fEledEtaCalo);
      tree_->SetBranchAddress("dphicalo",&fEledPhiCalo);
      tree_->SetBranchAddress("spp",&fEleSigmaIPhiIPhi);
      tree_->SetBranchAddress("sep",&fEleSigmaIEtaIPhi);
      tree_->SetBranchAddress("etawidth",&fEleSCEtaWidth);
      tree_->SetBranchAddress("phiwidth",&fEleSCPhiWidth);
      tree_->SetBranchAddress("PreShowerOverRaw",&fElePreShowerOverRaw);
      tree_->SetBranchAddress("gsfchi2",&fEleGsfTrackChi2OverNdof);
      tree_->SetBranchAddress("kfchi2",&fEleKFTrackChi2OverNDoF);
      tree_->SetBranchAddress("kfhits",&fEleKFTrackNHits);
      tree_->SetBranchAddress("kflayers",&fEleKFTrackNLayersWithMeasurement);
      tree_->SetBranchAddress("OneMinusSeedE1x5OverE5x5",&fEleOneMinusSeedE1x5OverE5x5);
      tree_->SetBranchAddress("PFMVA",&fElePFMVA);
      tree_->SetBranchAddress("trkIso03",&fEleTrkIso03);
      tree_->SetBranchAddress("ecalIso03",&fEleEMIso03);
      tree_->SetBranchAddress("hcalIso03",&fEleHadIso03);
      tree_->SetBranchAddress("trkIso04",&fEleTrkIso04);
      tree_->SetBranchAddress("ecalIso04",&fEleEMIso04);
      tree_->SetBranchAddress("hcalIso04",&fEleHadIso04);
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
      tree_->SetBranchAddress("PassTriggerDenominator",&fElePassTriggerDenominator);
      tree_->SetBranchAddress("IsEB",&fIsEB);
      tree_->SetBranchAddress("IsEE",&fIsEE);
      tree_->SetBranchAddress("SCRawEnergy",&fSCRawEnergy);
      tree_->SetBranchAddress("NClusters",&fNClusters);
      tree_->SetBranchAddress("EtaSeed",&fEtaSeed);
      tree_->SetBranchAddress("PhiSeed",&fPhiSeed);
      tree_->SetBranchAddress("ESeed",&fESeed);
      tree_->SetBranchAddress("E3x3Seed",&fE3x3Seed);
      tree_->SetBranchAddress("E5x5Seed",&fE5x5Seed);
      tree_->SetBranchAddress("EMaxSeed",&fEMaxSeed);
      tree_->SetBranchAddress("E2ndSeed",&fE2ndSeed);
      tree_->SetBranchAddress("ETopSeed",&fETopSeed);
      tree_->SetBranchAddress("EBottomSeed",&fEBottomSeed);
      tree_->SetBranchAddress("ELeftSeed",&fELeftSeed);
      tree_->SetBranchAddress("ERightSeed",&fERightSeed);
      tree_->SetBranchAddress("E2x5MaxSeed",&fE2x5MaxSeed);
      tree_->SetBranchAddress("E2x5TopSeed",&fE2x5TopSeed);
      tree_->SetBranchAddress("E2x5BottomSeed",&fE2x5BottomSeed);
      tree_->SetBranchAddress("E2x5LeftSeed",&fE2x5LeftSeed);
      tree_->SetBranchAddress("E2x5RightSeed",&fE2x5RightSeed);
      tree_->SetBranchAddress("IEtaSeed",&fIEtaSeed);
      tree_->SetBranchAddress("IPhiSeed",&fIPhiSeed);
      tree_->SetBranchAddress("EtaCrySeed",&fEtaCrySeed);
      tree_->SetBranchAddress("PhiCrySeed",&fPhiCrySeed);
      tree_->SetBranchAddress("ecalenergyerror",&fEcalEnergyError);
      tree_->SetBranchAddress("pmodegsf",&fGsfTrackPIn);
      tree_->SetBranchAddress("perror",&fTrackMomentumError);
      tree_->SetBranchAddress("pmeangsf",&fpmeankf);
      tree_->SetBranchAddress("pmeankf",&fpmeankf);
      tree_->SetBranchAddress("GeneratedEnergy",&fGeneratedEnergy);
      tree_->SetBranchAddress("GeneratedEnergyStatus1",&fGeneratedEnergyStatus1);
      tree_->SetBranchAddress("GeneratedEnergyStatus3",&fGeneratedEnergyStatus3);
      tree_->SetBranchAddress("EleMomentum_Regression_V0",&fEleMomentum_Regression_V0);
      tree_->SetBranchAddress("EleMomentum_Regression_V1",&fEleMomentum_Regression_V1);
      tree_->SetBranchAddress("EleMomentum_Regression_V2",&fEleMomentum_Regression_V2);
      tree_->SetBranchAddress("EleMomentumError_Regression_V0",&fEleMomentumError_Regression_V0);
      tree_->SetBranchAddress("EleMomentumError_Regression_V1",&fEleMomentumError_Regression_V1);
      tree_->SetBranchAddress("EleMomentumError_Regression_V2",&fEleMomentumError_Regression_V2);
      tree_->SetBranchAddress("EleClassification",&fEleClassification);
    }

    gErrorIgnoreLevel = currentState;
  }

}; 


#endif
