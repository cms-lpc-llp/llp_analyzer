#ifndef ControlSampleEvents_H
#define ControlSampleEvents_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"
#include <Rtypes.h>
#include "TLorentzVector.h"
#include "TH1F.h"

class ControlSampleEvents {
  
 public:
  
  /// bit map
  /// DON'T CHANGE ORDER
  enum TreeType { kTreeType_Default = 0,	              // dummy	      
		  kTreeType_OneLepton_Full = 1,               // lepton NOT added to MET
		  kTreeType_OneLepton_Reduced = 2,            // lepton NOT added to MET
		  kTreeType_OneLeptonAdd2MET_Full = 3,        // lepton added to MET
		  kTreeType_OneLeptonAdd2MET_Reduced = 4,     // lepton added to MET
		  kTreeType_Dilepton_Full = 5,                // lepton NOT added to MET
		  kTreeType_Dilepton_Reduced = 6,             // lepton NOT added to MET
		  kTreeType_DileptonAdd2MET_Full = 7,         // lepton added to MET
		  kTreeType_DileptonAdd2MET_Reduced = 8,      // lepton added to MET
		  kTreeType_Photon_Full = 9,                  // photon added to MET
		  kTreeType_Photon_Reduced = 10,              // photon added to MET
		  kTreeType_ZeroLepton_Full = 11,             // No Leptons, No Photons
		  kTreeType_ZeroLepton_Reduced = 12,          // No Leptons, No Photons
		  kTreeType_OneVetoLepton_Full = 13,          // lepton NOT added to MET
		  kTreeType_OneVetoLepton_Reduced = 14,       // lepton NOT added to MET
		  kTreeType_TightPlusVetoLepton_Full = 15,    // lepton NOT added to MET
		  kTreeType_TightPlusVetoLepton_Reduced = 16, // lepton NOT added to MET
		  kTreeType_OneTauLepton_Full = 17,           // lepton NOT added to MET
		  kTreeType_OneTauLepton_Reduced = 18,        // lepton NOT added to MET
  };
  
  //*******************************************
  //=== Process IDs  ====
  //*******************************************
  enum BkgProcessId { kData = 0,
		      kQCD = 1,
		      kWJets = 2,
		      kZJets = 3,
		      kTTJets = 4,
		      kSingleT = 5,
		      kVV = 6,
		      kHiggs = 7,
		      kSUSY = 99,
		      kUnknown = 999
  };

  /// variables
  Int_t                   option;
  Float_t                 genWeight;
  Float_t                 weight;
  Float_t                 pileupWeight;
  Float_t                 topPtWeight;
  Float_t                 ptTop;
  Float_t                 ptAntitop;
  Float_t                 trigWeight1L;
  Float_t                 muonEffWeight;
  Float_t                 eleEffWeight;
  Float_t                 muonRecoEffWeight;
  Float_t                 eleRecoEffWeight;
  Float_t                 btagW;
  Float_t                 btagW_up;
  Float_t                 btagW_down;
  UInt_t                  run;
  UInt_t                  lumi;
  UInt_t                  event;
  UInt_t                  processID;
  UInt_t                  NPU_0;
  UInt_t                  NPU_Minus1;
  UInt_t                  NPU_Plus1;
  UInt_t                  NPV;
  Float_t                 Rho;
  Float_t                 genHT;
  Int_t                   NISRJets;
  TLorentzVector          genlep1;
  TLorentzVector          genlep2;
  Int_t                   genlep1Type;
  Int_t                   genlep2Type;
  Bool_t                  HLTDecision[300];
  Int_t                   HLTPrescale[300];
  Float_t                 lep1Pt;
  Float_t                 lep1Eta;
  TLorentzVector          lep1;
  TLorentzVector          lep2;
  Int_t                   lep1Type;
  Int_t                   lep2Type;
  Int_t                   lep1MatchedGenLepIndex;
  Int_t                   lep2MatchedGenLepIndex;      
  Bool_t                  lep1PassVeto;
  Bool_t                  lep1PassLoose;
  Bool_t                  lep1PassMedium;
  Bool_t                  lep1PassTight;
  Bool_t                  lep1PassVetoID;
  Bool_t                  lep1PassLooseID;
  Bool_t                  lep1PassMediumID;
  Bool_t                  lep1PassTightID;
  Bool_t                  lep1PassVetoIso;
  Bool_t                  lep1PassLooseIso;
  Bool_t                  lep1PassMediumIso;
  Bool_t                  lep1PassTightIso;
  Float_t                 lep1MinDRToBJet;
  Float_t                 lep1Activity;
  Bool_t                  lep2PassVeto;
  Bool_t                  lep2PassLoose;
  Bool_t                  lep2PassMedium;
  Bool_t                  lep2PassTight;
  Bool_t                  lep2PassVetoID;
  Bool_t                  lep2PassLooseID;
  Bool_t                  lep2PassMediumID;
  Bool_t                  lep2PassTightID;
  Bool_t                  lep2PassVetoIso;
  Bool_t                  lep2PassLooseIso;
  Bool_t                  lep2PassMediumIso;
  Bool_t                  lep2PassTightIso;
  Float_t                 lep2MinDRToBJet;     
  Float_t                 lep2Activity;
  Int_t                   NGenBJets;
  TLorentzVector          bjet1;
  TLorentzVector          bjet2;
  Bool_t                  bjet1PassLoose;
  Bool_t                  bjet1PassMedium;
  Bool_t                  bjet1PassTight;
  Bool_t                  bjet2PassLoose;
  Bool_t                  bjet2PassMedium;
  Bool_t                  bjet2PassTight;
  TLorentzVector          jet1;
  TLorentzVector          jet2;      
  Bool_t                  jet1PassCSVLoose;
  Bool_t                  jet1PassCSVMedium;
  Bool_t                  jet1PassCSVTight;
  Bool_t                  jet2PassCSVLoose;
  Bool_t                  jet2PassCSVMedium;
  Bool_t                  jet2PassCSVTight;
  Int_t                   nWTags;
  Int_t                   nTopTags;
  Float_t                 wTagScaleFactor;
  Float_t                 topTagScaleFactor;
  Float_t                 MR;
  Float_t                 Rsq;
  Float_t                 RsqnoHF;
  Float_t                 MR_NoDilepton;
  Float_t                 Rsq_NoDilepton;
  Float_t                 MR_NoLeadJet;
  Float_t                 Rsq_NoLeadJet;
  Float_t                 MET;
  Float_t                 METRaw;
  Float_t                 METnoHF;
  Float_t                 MET_NoDilepton;
  Float_t                 MET_NoLeadJet;
  Float_t                 minDPhi;
  Float_t                 minDPhiN;
  Float_t                 dPhiRazor;
  UInt_t                  NJets40;
  UInt_t                  NJets80;
  UInt_t                  NBJetsLoose;
  UInt_t                  NBJetsMedium;
  UInt_t                  NBJetsTight;
  Float_t                 HT;
  Float_t                 MHT;
  Float_t                 MHTnoHF;
  Float_t                 lep1MT;
  Float_t                 lep1MTnoHF;
  Float_t                 mll;
  Float_t                 u1, u2; // Parallel and perpendicular component of the recoil vector to Z
  Bool_t                  Flag_HBHENoiseFilter;//
  Bool_t                  Flag_HBHEIsoNoiseFilter;//
  Bool_t                  Flag_badChargedCandidateFilter;//
  Bool_t                  Flag_badMuonFilter;//
  Bool_t                  Flag_badGlobalMuonFilter;//
  Bool_t                  Flag_duplicateMuonFilter;//
  Bool_t                  Flag_CSCTightHaloFilter;
  Bool_t                  Flag_hcalLaserEventFilter; //
  Bool_t                  Flag_EcalDeadCellTriggerPrimitiveFilter; //
  Bool_t                  Flag_goodVertices;
  Bool_t                  Flag_trackingFailureFilter;//
  Bool_t                  Flag_eeBadScFilter; //
  Bool_t                  Flag_ecalLaserCorrFilter;
  Bool_t                  Flag_trkPOGFilters; //
  Bool_t                  Flag_trkPOG_manystripclus53X;
  Bool_t                  Flag_trkPOG_toomanystripclus53X;
  Bool_t                  Flag_trkPOG_logErrorTooManyClusters;
  Bool_t                  Flag_METFilters;
  UInt_t                  nVetoMuons;
  UInt_t                  nMediumMuons;
  UInt_t                  nLooseMuons;
  UInt_t                  nTightMuons;
  UInt_t                  nSelectedPhotons;
  Float_t                 recoZmass;
  Float_t                 recoZpt;
  Float_t                 recoWpt;
  Float_t                 recoWphi;
  Float_t                 genZpt;
  Float_t                 genZphi;
  Float_t                 genWpt;
  Float_t                 genWphi;
  Float_t                 MT2;
  Float_t                 genJetMR;
  Float_t                 genJetHT;

  TLorentzVector          pho1;
  TLorentzVector          pho2;
  Bool_t                  pho1PassLoose;
  Bool_t                  pho1PassTight;
  Bool_t                  pho1PassMedium;
  Bool_t                  pho1HLTFilter[50];
  Bool_t                  HLT_Dimuon;
  Bool_t                  HLT_SingleMu;
  Bool_t                  HLT_Photon;
  Bool_t                  HLT_Razor;
  Bool_t                  HLT_Photon36;
  Bool_t                  HLT_Photon50;
  Bool_t                  HLT_Photon75;
  Bool_t                  HLT_Photon90;
  Bool_t                  HLT_Photon120;
  Bool_t                  HLT_Photon165;
  Float_t                 MR_NoZ;
  Float_t                 Rsq_NoZ;
  Float_t                 MR_NoW;
  Float_t                 Rsq_NoW;
  Float_t                 MR_NoPho;
  Float_t                 Rsq_NoPho;
  Float_t                 HT_NoZ;
  Float_t                 HT_NoW;
  Float_t                 HT_NoPho;
  Float_t                 dPhiRazor_NoZ;
  Float_t                 dPhiRazor_NoW;
  Float_t                 dPhiRazor_NoPho;
  Float_t                 MET_NoZ;
  Float_t                 MET_NoW;
  Float_t                 MET_NoPho;
  Float_t                 METPhi_NoPho;
  Float_t                 METPhi_NoW;
  Float_t                 METPhi_NoZ;
  Float_t                 METPhi;
  Float_t                 METnoHFPhi;
  Float_t                 METRawPhi;
  UInt_t                  NJets_NoZ;
  UInt_t                  NJets_NoW;
  UInt_t                  NJets_NoPho;
  UInt_t                  NJets80_NoZ;
  UInt_t                  NJets80_NoW;
  UInt_t                  NJets80_NoPho;
  Int_t                   pho1_motherID;
  Float_t                 pho1_sigmaietaieta; 
  Float_t                 pho1_chargediso;
  Float_t                 pho1_photoniso;
  Float_t                 pho1_neutralhadroniso;
  Float_t                 pho1_hOverE;
  Float_t                 pho1_pfiso;
  Float_t                 minDRGenPhotonToParton;

  float metType1PtJetResUp;
  float metType1PtJetResDown;
  float metType1PtJetEnUp;
  float metType1PtJetEnDown;
  float metType1PtMuonEnUp;
  float metType1PtMuonEnDown;
  float metType1PtElectronEnUp;
  float metType1PtElectronEnDown;
  float metType1PtTauEnUp;
  float metType1PtTauEnDown;
  float metType1PtUnclusteredEnUp;
  float metType1PtUnclusteredEnDown;
  float metType1PtPhotonEnUp;
  float metType1PtPhotonEnDown;       
  float metType1PtMETUncertaintySize;
  float metType1PtJetResUpSmear;
  float metType1PtJetResDownSmear;
  float metType1PtMETFullUncertaintySize;
  
  float metType1PhiJetResUp;
  float metType1PhiJetResDown;
  float metType1PhiJetEnUp;
  float metType1PhiJetEnDown;
  float metType1PhiMuonEnUp;
  float metType1PhiMuonEnDown;
  float metType1PhiElectronEnUp;
  float metType1PhiElectronEnDown;
  float metType1PhiTauEnUp;
  float metType1PhiTauEnDown;
  float metType1PhiUnclusteredEnUp;
  float metType1PhiUnclusteredEnDown;
  float metType1PhiPhotonEnUp;
  float metType1PhiPhotonEnDown;
  float metType1PhiMETUncertaintySize;
  float metType1PhiJetResUpSmear;
  float metType1PhiJetResDownSmear;
  float metType1PhiMETFullUncertaintySize;

 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  
  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  /// default constructor  
  ControlSampleEvents()  {
    genlep1Ptr  = &genlep1;
    genlep2Ptr  = &genlep2;
    lep1Ptr     = &lep1;
    lep2Ptr     = &lep2;
    bjet1Ptr    = &bjet1;
    bjet2Ptr    = &bjet2;       
    jet1Ptr     = &jet1;
    jet2Ptr     = &jet2;       
    pho1Ptr     = &pho1;
    pho2Ptr     = &pho2;
  };

  /// default destructor
  ~ControlSampleEvents(){ 
    if (f_) f_->Close();  
  };
    
  /// initialize varibles and fill list of available variables
  void InitVariables() {
    option               = -1;
    genWeight            = 0.0;
    weight               = 0.0;
    pileupWeight         = 0.0;
    topPtWeight          = 1.0;
    ptTop                = -1.0;
    ptAntitop            = -1.0;
    trigWeight1L         = 1.0;
    muonEffWeight        = 1.0;
    eleEffWeight         = 1.0;
    muonRecoEffWeight    = 1.0;
    eleRecoEffWeight     = 1.0;
    btagW                = 1.0;
    btagW_up             = 1.0;
    btagW_down           = 1.0;
    run                  = 0;
    lumi                 = 0;
    event                = 0;
    processID            = ControlSampleEvents::kUnknown;
    NPU_0                = 0;
    NPU_Minus1           = 0;
    NPU_Plus1            = 0;
    NPV                  = 0;
    Rho                  = 0.0;
    genHT                = 0.0;
    NISRJets             = 0;
    genlep1              = TLorentzVector();
    genlep2              = TLorentzVector();
    genlep1Type          = 0.0;
    genlep2Type          = 0.0;
    for(int i=0;i<300;++i) HLTDecision[i] = false;
    for(int i=0;i<300;++i) HLTPrescale[i] = 0;
    lep1Pt               = 0.0;
    lep1Eta              = 0.0;
    lep1                 = TLorentzVector();
    lep2                 = TLorentzVector();
    lep1Type             = 0.0;
    lep2Type             = 0.0;
    lep1MatchedGenLepIndex = -1;
    lep2MatchedGenLepIndex = -1;
    lep1PassVeto         = 0.0;
    lep1PassLoose        = 0.0;
    lep1PassMedium       = 0.0;
    lep1PassTight        = 0.0;
    lep1PassVetoID       = 0.0;
    lep1PassLooseID      = 0.0;
    lep1PassMediumID     = 0.0;
    lep1PassTightID      = 0.0;
    lep1PassVetoIso      = 0.0;
    lep1PassLooseIso     = 0.0;
    lep1PassMediumIso    = 0.0;
    lep1PassTightIso     = 0.0;
    lep1MinDRToBJet      = 0.0;
    lep1Activity         = 0.0;
    lep2PassVeto         = 0.0;
    lep2PassLoose        = 0.0;
    lep2PassMedium       = 0.0;
    lep2PassTight        = 0.0;
    lep2PassVetoID       = 0.0;
    lep2PassLooseID      = 0.0;
    lep2PassMediumID     = 0.0;
    lep2PassTightID      = 0.0;
    lep2PassVetoIso      = 0.0;
    lep2PassLooseIso     = 0.0;
    lep2PassMediumIso    = 0.0;
    lep2PassTightIso     = 0.0;
    lep2MinDRToBJet      = 0.0;
    lep2Activity         = 0.0;
    NGenBJets            = -1;
    bjet1                = TLorentzVector();
    bjet2                = TLorentzVector();
    bjet1PassLoose       = 0.0;
    bjet1PassMedium      = 0.0;
    bjet1PassTight       = 0.0;
    bjet2PassLoose       = 0.0;
    bjet2PassMedium      = 0.0;
    bjet2PassTight       = 0.0;
    jet1                 = TLorentzVector();
    jet2                 = TLorentzVector();
    jet1PassCSVLoose     = 0.0;
    jet1PassCSVMedium    = 0.0;
    jet1PassCSVTight     = 0.0;
    jet2PassCSVLoose     = 0.0;
    jet2PassCSVMedium    = 0.0;
    jet2PassCSVTight     = 0.0;
    nWTags               = 0;
    nTopTags             = 0;
    wTagScaleFactor      = 1.0; 
    topTagScaleFactor    = 1.0;
    MR                   = 0.0;
    Rsq                  = 0.0;
    RsqnoHF              = 0.0;
    MR_NoDilepton        = 0.0;
    Rsq_NoDilepton       = 0.0;
    MR_NoLeadJet         = 0.0;
    Rsq_NoLeadJet        = 0.0;
    MET                  = 0.0;
    METnoHF              = 0.0;
    METRaw               = 0.0;
    MET_NoDilepton       = 0.0;
    MET_NoLeadJet        = 0.0;
    minDPhi              = 0.0;
    minDPhiN             = 0.0;
    dPhiRazor            = 0.0;
    NJets40              = 0;
    NJets80              = 0;
    NBJetsLoose          = 0;
    NBJetsMedium         = 0;
    NBJetsTight          = 0;
    HT                   = 0.0;      
    MHT                  = 0.;
    MHTnoHF              = 0.;
    lep1MT               = 0.0;  
    lep1MTnoHF           = 0.0;  
    mll                  = 0.0;
    Flag_HBHENoiseFilter = 0.0;//
    Flag_HBHEIsoNoiseFilter = 0.0;//
    Flag_badChargedCandidateFilter = 0.0;//
    Flag_badMuonFilter = 0.0;//
    Flag_badGlobalMuonFilter = 0.0;//
    Flag_duplicateMuonFilter = 0.0;//
    Flag_CSCTightHaloFilter = 0.0;
    Flag_hcalLaserEventFilter = 0.0; //
    Flag_EcalDeadCellTriggerPrimitiveFilter = 0.0; //
    Flag_goodVertices = 0.0;
    Flag_trackingFailureFilter = 0.0;//
    Flag_eeBadScFilter = 0.0; //
    Flag_ecalLaserCorrFilter = 0.0;
    Flag_trkPOGFilters = 0.0; //
    Flag_trkPOG_manystripclus53X = 0.0;
    Flag_trkPOG_toomanystripclus53X = 0.0;
    Flag_trkPOG_logErrorTooManyClusters = 0.0;
    Flag_METFilters = 0.;
    nVetoMuons      = 0;
    nLooseMuons     = 0;
    nMediumMuons    = 0;
    nTightMuons     = 0;
    nSelectedPhotons = 0;
    recoZmass       = -99.;
    recoZpt         = -99.;
    recoWpt         = -99.;
    recoWphi        = -99.;
    genWpt          = -99.;
    genWphi         = -99.;
    genZpt          = -99.;
    genZphi         = -99.;
    MT2             = -99.;
    genJetHT        = -99.;
    genJetMR        = -99.;

    pho1 = TLorentzVector();
    pho2 = TLorentzVector();

    pho1PassLoose = false;
    pho1PassTight = false;
    pho1PassMedium = false;
    HLT_Dimuon = false;
    HLT_SingleMu = false;
    HLT_Photon = false;
    HLT_Razor = false;
    HLT_Photon36 = false;
    HLT_Photon50 = false;
    HLT_Photon75 = false;
    HLT_Photon90 = false;
    HLT_Photon120 = false;
    HLT_Photon165 = false;
    MR_NoZ = 0.0 ; 
    Rsq_NoZ = 0.0 ; 
    MR_NoW = 0.0 ; 
    Rsq_NoW = 0.0 ; 
    MR_NoPho = 0.0 ; 
    Rsq_NoPho = 0.0 ; 
    HT_NoZ = 0.0 ; 
    HT_NoW = 0.0 ; 
    HT_NoPho = 0.0 ; 
    dPhiRazor_NoZ = 0.0 ; 
    dPhiRazor_NoW = 0.0 ; 
    dPhiRazor_NoPho = 0.0 ; 
    MET_NoZ = 0.0 ; 
    MET_NoW = 0.0 ; 
    MET_NoPho = 0.0 ;
    METPhi_NoPho = -99.;
    METPhi_NoZ = -99.;
    METPhi_NoW = -99.;
    METPhi = -99.;
    METRawPhi = -99.;
    METnoHFPhi = -99.;
    u1 = -99.;
    u2 = -99.;
    NJets_NoZ = 0; 
    NJets_NoW = 0; 
    NJets_NoPho = 0 ; 
    NJets80_NoZ = 0 ; 
    NJets80_NoW = 0 ; 
    NJets80_NoPho = 0 ; 
    pho1_motherID = 0;
    pho1_sigmaietaieta = -999.;
    pho1_chargediso = -999.;
    pho1_photoniso = -999.;
    pho1_neutralhadroniso = -999.;
    pho1_hOverE = -999.;
    pho1_pfiso = -999.;
    minDRGenPhotonToParton = -999;

    metType1PtJetResUp=-999.;
    metType1PtJetResDown=-999.;
    metType1PtJetEnUp=-999.;
    metType1PtJetEnDown=-999.;
    metType1PtMuonEnUp=-999.;
    metType1PtMuonEnDown=-999.;
    metType1PtElectronEnUp=-999.;
    metType1PtElectronEnDown=-999.;
    metType1PtTauEnUp=-999.;
    metType1PtTauEnDown=-999.;
    metType1PtUnclusteredEnUp=-999.;
    metType1PtUnclusteredEnDown=-999.;
    metType1PtPhotonEnUp=-999.;
    metType1PtPhotonEnDown=-999.;
    metType1PtMETUncertaintySize=-999.;
    metType1PtJetResUpSmear=-999.;
    metType1PtJetResDownSmear=-999.;
    metType1PtMETFullUncertaintySize=-999.;
  
    metType1PhiJetResUp=-999.;
    metType1PhiJetResDown=-999.;
    metType1PhiJetEnUp=-999.;
    metType1PhiJetEnDown=-999.;
    metType1PhiMuonEnUp=-999.;
    metType1PhiMuonEnDown=-999.;
    metType1PhiElectronEnUp=-999.;
    metType1PhiElectronEnDown=-999.;
    metType1PhiTauEnUp=-999.;
    metType1PhiTauEnDown=-999.;
    metType1PhiUnclusteredEnUp=-999.;
    metType1PhiUnclusteredEnDown=-999.;
    metType1PhiPhotonEnUp=-999.;
    metType1PhiPhotonEnDown=-999.;
    metType1PhiMETUncertaintySize=-999.;
    metType1PhiJetResUpSmear=-999.;
    metType1PhiJetResDownSmear=-999.;
    metType1PhiMETFullUncertaintySize=-999.;


  }
    
  /// load a ControlSampleEvents
  void LoadTree(const char* file, int treeType = kTreeType_Default){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->Get("ControlSampleEvent"));
    InitTree(treeType);
    assert(tree_);
  }
    
  /// create a ControlSampleEvents
  void CreateTree(int treeType = kTreeType_Default){
    tree_ = new TTree("ControlSampleEvent","ControlSampleEvent");
    f_ = 0;

    //book the branches that go in all types of trees
    tree_->Branch("option",&option,"option/I");
    tree_->Branch("genWeight",&genWeight,"genWeight/F");
    tree_->Branch("weight",&weight,"weight/F");
    tree_->Branch("pileupWeight",&pileupWeight,"pileupWeight/F");
    tree_->Branch("btagW",&btagW,"btagW/F");
    tree_->Branch("btagW_up",&btagW_up,"btagW_up/F");
    tree_->Branch("btagW_down",&btagW_down,"btagW_down/F");
    tree_->Branch("run",&run,"run/i");
    tree_->Branch("lumi",&lumi,"lumi/i");
    tree_->Branch("event",&event,"event/i");
    tree_->Branch("NPU_0",&NPU_0,"NPU_0/i");
    tree_->Branch("NPV",&NPV,"NPV/i");
    tree_->Branch("MR",&MR,"MR/F");
    tree_->Branch("Rsq",&Rsq,"Rsq/F");
    tree_->Branch("MET",&MET,"MET/F");
    tree_->Branch("NJets40",&NJets40,"NJets40/i");
    tree_->Branch("NJets80",&NJets80,"NJets80/i");
    tree_->Branch("NGenBJets",&NGenBJets,"NGenBJets/I");
    tree_->Branch("NBJetsLoose",&NBJetsLoose,"NBJetsLoose/i");
    tree_->Branch("NBJetsMedium",&NBJetsMedium,"NBJetsMedium/i");
    tree_->Branch("NBJetsTight",&NBJetsTight,"NBJetsTight/i");
    tree_->Branch("genHT",&genHT,"genHT/F");
    tree_->Branch("NISRJets",&NISRJets,"NISRJets/I");
    tree_->Branch("nWTags", &nWTags, "nWTags/I");
    tree_->Branch("nTopTags", &nTopTags, "nTopTags/I");
    tree_->Branch("wTagScaleFactor", &wTagScaleFactor, "wTagScaleFactor/F");
    tree_->Branch("topTagScaleFactor", &topTagScaleFactor, "topTagScaleFactor/F");
    // noise filters
    tree_->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter,"Flag_HBHENoiseFilter/O");
    tree_->Branch("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter,"Flag_HBHEIsoNoiseFilter/O");
    tree_->Branch("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter,"Flag_badChargedCandidateFilter/O");
    tree_->Branch("Flag_badMuonFilter", &Flag_badMuonFilter,"Flag_badMuonFilter/O");
    tree_->Branch("Flag_badGlobalMuonFilter", &Flag_badGlobalMuonFilter,"Flag_badGlobalMuonFilter/O");
    tree_->Branch("Flag_duplicateMuonFilter", &Flag_duplicateMuonFilter,"Flag_duplicateMuonFilter/O");
    tree_->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter,"Flag_CSCTightHaloFilter/O");
    tree_->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter,"Flag_hcalLaserEventFilter/O");
    tree_->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter,"Flag_EcalDeadCellTriggerPrimitiveFilter/O");
    tree_->Branch("Flag_goodVertices", &Flag_goodVertices,"Flag_goodVertices/O");
    tree_->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter,"Flag_trackingFailureFilter/O");
    tree_->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter,"Flag_eeBadScFilter/O");
    tree_->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter,"Flag_ecalLaserCorrFilter/O");
    tree_->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters,"Flag_trkPOGFilters/O");
    tree_->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X,"Flag_trkPOG_manystripclus53X/O");
    tree_->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X,"Flag_trkPOG_toomanystripclus53X/O");
    tree_->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters,"Flag_trkPOG_logErrorTooManyClusters/O");
    tree_->Branch("Flag_METFilters", &Flag_METFilters,"Flag_METFilters/O");	

    // book the branches that go into only One Lepton trees
    if (treeType == kTreeType_OneLepton_Reduced 
	|| treeType == kTreeType_ZeroLepton_Reduced 
	|| treeType == kTreeType_OneVetoLepton_Reduced 
	|| treeType == kTreeType_OneTauLepton_Reduced 
	) {
      tree_->Branch("genlep1Type",&genlep1Type,"genlep1Type/I");
      tree_->Branch("lep1Type",&lep1Type,"lep1Type/I");
      tree_->Branch("lep1PassVeto",&lep1PassVeto,"lep1PassVeto/O");
      tree_->Branch("lep1PassLoose",&lep1PassLoose,"lep1PassLoose/O");
      tree_->Branch("lep1PassMedium",&lep1PassMedium,"lep1PassMedium/O");
      tree_->Branch("lep1PassTight",&lep1PassTight,"lep1PassTight/O");
      tree_->Branch("lep1MatchedGenLepIndex",&lep1MatchedGenLepIndex,"lep1MatchedGenLepIndex/I");
      tree_->Branch("lep1MT",&lep1MT,"lep1MT/F");
      tree_->Branch("lep1MTnoHF",&lep1MTnoHF,"lep1MTnoHF/F");
      tree_->Branch("topPtWeight",&topPtWeight,"topPtWeight/F");
      tree_->Branch("ptTop",&ptTop,"ptTop/F");
      tree_->Branch("ptAntitop",&ptAntitop,"ptAntitop/F");
      tree_->Branch("trigWeight1L",&trigWeight1L,"trigWeight1L/F");
      tree_->Branch("muonEffWeight",&muonEffWeight,"muonEffWeight/F");
      tree_->Branch("eleEffWeight",&eleEffWeight,"eleEffWeight/F");
      tree_->Branch("muonRecoEffWeight",&muonRecoEffWeight,"muonRecoEffWeight/F");
      tree_->Branch("eleRecoEffWeight",&eleRecoEffWeight,"eleRecoEffWeight/F");
      tree_->Branch("MET",&MET,"MET/F");
      tree_->Branch("METPhi",&METPhi,"METPhi/F");
      tree_->Branch("HLTDecision",&HLTDecision,"HLTDecision[300]/O");
      tree_->Branch("HT",&HT,"HT/F");
      tree_->Branch("lep1Pt",&lep1Pt,"lep1Pt/F");
      tree_->Branch("lep1Eta",&lep1Eta,"lep1Eta/F");
      tree_->Branch("METnoHF", &METnoHF, "METnoHF/F");
      tree_->Branch("RsqnoHF", &RsqnoHF, "RsqnoHF/F");
      tree_->Branch("MHT", &MHT, "MHT/F");
      tree_->Branch("MHTnoHF", &MHTnoHF, "MHTnoHF/F");
      tree_->Branch("lep1" ,&lep1Ptr);

      tree_->Branch("metType1PtJetResUp", &metType1PtJetResUp, "metType1PtJetResUp/F");
      tree_->Branch("metType1PtJetResDown", &metType1PtJetResDown, "metType1PtJetResDown/F");
      tree_->Branch("metType1PtJetEnUp", &metType1PtJetEnUp, "metType1PtJetEnUp/F");
      tree_->Branch("metType1PtJetEnDown", &metType1PtJetEnDown, "metType1PtJetEnDown/F");
      tree_->Branch("metType1PtMuonEnUp", &metType1PtMuonEnUp, "metType1PtMuonEnUp/F");
      tree_->Branch("metType1PtMuonEnDown", &metType1PtMuonEnDown, "metType1PtMuonEnDown/F");
      tree_->Branch("metType1PtElectronEnUp", &metType1PtElectronEnUp, "metType1PtElectronEnUp/F");
      tree_->Branch("metType1PtElectronEnDown", &metType1PtElectronEnDown, "metType1PtElectronEnDown/F");
      tree_->Branch("metType1PtTauEnUp", &metType1PtTauEnUp, "metType1PtTauEnUp/F");
      tree_->Branch("metType1PtTauEnDown", &metType1PtTauEnDown, "metType1PtTauEnDown/F");
      tree_->Branch("metType1PtUnclusteredEnUp", &metType1PtUnclusteredEnUp, "metType1PtUnclusteredEnUp/F");
      tree_->Branch("metType1PtUnclusteredEnDown", &metType1PtUnclusteredEnDown, "metType1PtUnclusteredEnDown/F");
      tree_->Branch("metType1PtPhotonEnUp", &metType1PtPhotonEnUp, "metType1PtPhotonEnUp/F");
      tree_->Branch("metType1PtPhotonEnDown", &metType1PtPhotonEnDown, "metType1PtPhotonEnDown/F");
      
      tree_->Branch("metType1PhiJetResUp", &metType1PhiJetResUp, "metType1PhiJetResUp/F");
      tree_->Branch("metType1PhiJetResDown", &metType1PhiJetResDown, "metType1PhiJetResDown/F");
      tree_->Branch("metType1PhiJetEnUp", &metType1PhiJetEnUp, "metType1PhiJetEnUp/F");
      tree_->Branch("metType1PhiJetEnDown", &metType1PhiJetEnDown, "metType1PhiJetEnDown/F");
      tree_->Branch("metType1PhiMuonEnUp", &metType1PhiMuonEnUp, "metType1PhiMuonEnUp/F");
      tree_->Branch("metType1PhiMuonEnDown", &metType1PhiMuonEnDown, "metType1PhiMuonEnDown/F");
      tree_->Branch("metType1PhiElectronEnUp", &metType1PhiElectronEnUp, "metType1PhiElectronEnUp/F");
      tree_->Branch("metType1PhiElectronEnDown", &metType1PhiElectronEnDown, "metType1PhiElectronEnDown/F");
      tree_->Branch("metType1PhiTauEnUp", &metType1PhiTauEnUp, "metType1PhiTauEnUp/F");
      tree_->Branch("metType1PhiTauEnDown", &metType1PhiTauEnDown, "metType1PhiTauEnDown/F");
      tree_->Branch("metType1PhiUnclusteredEnUp", &metType1PhiUnclusteredEnUp, "metType1PhiUnclusteredEnUp/F");
      tree_->Branch("metType1PhiUnclusteredEnDown", &metType1PhiUnclusteredEnDown, "metType1PhiUnclusteredEnDown/F");
      tree_->Branch("metType1PhiPhotonEnUp", &metType1PhiPhotonEnUp, "metType1PhiPhotonEnUp/F");
      tree_->Branch("metType1PhiPhotonEnDown", &metType1PhiPhotonEnDown, "metType1PhiPhotonEnDown/F");
   }
  
    if (treeType == kTreeType_OneLepton_Full
	|| treeType == kTreeType_ZeroLepton_Full 
	|| treeType == kTreeType_OneVetoLepton_Full 
	|| treeType == kTreeType_OneTauLepton_Full 
	) {
      tree_->Branch("NPU_Minus1",&NPU_Minus1,"NPU_Minus1/i");
      tree_->Branch("NPU_Plus1",&NPU_Plus1,"NPU_Plus1/i");
      tree_->Branch("event",&event,"event/i");
      tree_->Branch("processID",&processID,"processID/i");
      tree_->Branch("HLTDecision",&HLTDecision,"HLTDecision[300]/O");
      tree_->Branch("genlep1Type",&genlep1Type,"genlep1Type/I");
      tree_->Branch("lep1Type",&lep1Type,"lep1Type/I");
      tree_->Branch("lep1PassVeto",&lep1PassVeto,"lep1PassVeto/O");
      tree_->Branch("lep1PassLoose",&lep1PassLoose,"lep1PassLoose/O");
      tree_->Branch("lep1PassMedium",&lep1PassMedium,"lep1PassMedium/O");
      tree_->Branch("lep1PassTight",&lep1PassTight,"lep1PassTight/O");
      tree_->Branch("lep1MatchedGenLepIndex",&lep1MatchedGenLepIndex,"lep1MatchedGenLepIndex/I");
      tree_->Branch("lep1MT",&lep1MT,"lep1MT/F");
      tree_->Branch("lep1Activity",&lep1Activity,"lep1Activity/F");
      tree_->Branch("topPtWeight",&topPtWeight,"topPtWeight/F");
      tree_->Branch("ptTop",&ptTop,"ptTop/F");
      tree_->Branch("ptAntitop",&ptAntitop,"ptAntitop/F");
      tree_->Branch("trigWeight1L",&trigWeight1L,"trigWeight1L/F");
      tree_->Branch("eleEffWeight",&eleEffWeight,"eleEffWeight/F");
      tree_->Branch("muonEffWeight",&muonEffWeight,"muonEffWeight/F");
      tree_->Branch("eleRecoEffWeight",&eleRecoEffWeight,"eleRecoEffWeight/F");
      tree_->Branch("muonRecoEffWeight",&muonRecoEffWeight,"muonRecoEffWeight/F");
      tree_->Branch("MET",&MET,"MET/F");
      tree_->Branch("METPhi",&METPhi,"METPhi/F");
      tree_->Branch("METnoHF", &METnoHF, "METnoHF/F");
      tree_->Branch("METnoHFPhi",&METnoHFPhi,"METnoHFPhi/F");
      tree_->Branch("METRaw",&METRaw,"METRaw/F");
      tree_->Branch("METRawPhi",&METRawPhi,"METRawPhi/F");

      tree_->Branch("bjet1PassLoose",&bjet1PassLoose,"bjet1PassLoose/O");
      tree_->Branch("bjet1PassMedium",&bjet1PassMedium,"bjet1PassMedium/O");
      tree_->Branch("bjet1PassTight",&bjet1PassTight,"bjet1PassTight/O");
      tree_->Branch("bjet2PassLoose",&bjet2PassLoose,"bjet2PassLoose/O");
      tree_->Branch("bjet2PassMedium",&bjet2PassMedium,"bjet2PassMedium/O");
      tree_->Branch("bjet2PassTight",&bjet2PassTight,"bjet2PassTight/O");
      tree_->Branch("jet1PassCSVLoose",&jet1PassCSVLoose,"jet1PassCSVLoose/O");
      tree_->Branch("jet1PassCSVMedium",&jet1PassCSVMedium,"jet1PassCSVMedium/O");
      tree_->Branch("jet1PassCSVTight",&jet1PassCSVTight,"jet1PassCSVTight/O");
      tree_->Branch("jet2PassCSVLoose",&jet2PassCSVLoose,"jet2PassCSVLoose/O");
      tree_->Branch("jet2PassCSVMedium",&jet2PassCSVMedium,"jet2PassCSVMedium/O");
      tree_->Branch("jet2PassCSVTight",&jet2PassCSVTight,"jet2PassCSVTight/O");
      tree_->Branch("dPhiRazor",&dPhiRazor,"dPhiRazor/F");
      tree_->Branch("HT",&HT,"HT/F");	  
      tree_->Branch("genWpt",&genWpt,"genWpt/F");
      tree_->Branch("genWphi",&genWphi,"genWphi/F");
      tree_->Branch("genlep1", "TLorentzVector", &genlep1Ptr);
      tree_->Branch("lep1",    "TLorentzVector", &lep1Ptr);
      tree_->Branch("bjet1",   "TLorentzVector", &bjet1Ptr);
      tree_->Branch("bjet2",   "TLorentzVector", &bjet2Ptr);
      tree_->Branch("jet1",    "TLorentzVector", &jet1Ptr);
      tree_->Branch("jet2",    "TLorentzVector", &jet2Ptr);
    }
  
    if (treeType == kTreeType_Dilepton_Full
	|| treeType == kTreeType_TightPlusVetoLepton_Full
	) {
      tree_->Branch("NPU_Minus1",&NPU_Minus1,"NPU_Minus1/i");
      tree_->Branch("NPU_Plus1",&NPU_Plus1,"NPU_Plus1/i");
      tree_->Branch("event",&event,"event/i");
      tree_->Branch("processID",&processID,"processID/i");
      tree_->Branch("HLTDecision",&HLTDecision,"HLTDecision[300]/O");
      tree_->Branch("genlep1Type",&genlep1Type,"genlep1Type/I");
      tree_->Branch("lep1Type",&lep1Type,"lep1Type/I");
      tree_->Branch("lep1MatchedGenLepIndex",&lep1MatchedGenLepIndex,"lep1MatchedGenLepIndex/I");
      tree_->Branch("genlep2Type",&genlep2Type,"genlep2Type/I");
      tree_->Branch("lep2Type",&lep2Type,"lep2Type/I");
      tree_->Branch("lep2MatchedGenLepIndex",&lep2MatchedGenLepIndex,"lep2MatchedGenLepIndex/I");
      tree_->Branch("lep1PassVeto",&lep1PassVeto,"lep1PassVeto/O");
      tree_->Branch("lep1PassLoose",&lep1PassLoose,"lep1PassLoose/O");
      tree_->Branch("lep1PassMedium",&lep1PassMedium,"lep1PassMedium/O");
      tree_->Branch("lep1PassTight",&lep1PassTight,"lep1PassTight/O");
      tree_->Branch("lep1PassVetoID",&lep1PassVetoID,"lep1PassVetoID/O");
      tree_->Branch("lep1PassLooseID",&lep1PassLooseID,"lep1PassLooseID/O");
      tree_->Branch("lep1PassMediumID",&lep1PassMediumID,"lep1PassMediumID/O");
      tree_->Branch("lep1PassTightID",&lep1PassTightID,"lep1PassTightID/O");
      tree_->Branch("lep1PassVetoIso",&lep1PassVetoIso,"lep1PassVetoIso/O");
      tree_->Branch("lep1PassLooseIso",&lep1PassLooseIso,"lep1PassLooseIso/O");
      tree_->Branch("lep1PassMediumIso",&lep1PassMediumIso,"lep1PassMediumIso/O");
      tree_->Branch("lep1PassTightIso",&lep1PassTightIso,"lep1PassTightIso/O");
      tree_->Branch("lep1MinDRToBJet",&lep1MinDRToBJet,"lep1MinDRToBJet/F");
      tree_->Branch("lep1Activity",&lep1Activity,"lep1Activity/F");
      tree_->Branch("lep2PassVeto",&lep2PassVeto,"lep2PassVeto/O");
      tree_->Branch("lep2PassLoose",&lep2PassLoose,"lep2PassLoose/O");
      tree_->Branch("lep2PassMedium",&lep2PassMedium,"lep2PassMedium/O");
      tree_->Branch("lep2PassTight",&lep2PassTight,"lep2PassTight/O");
      tree_->Branch("lep2PassVetoID",&lep2PassVetoID,"lep2PassVetoID/O");
      tree_->Branch("lep2PassLooseID",&lep2PassLooseID,"lep2PassLooseID/O");
      tree_->Branch("lep2PassMediumID",&lep2PassMediumID,"lep2PassMediumID/O");
      tree_->Branch("lep2PassTightID",&lep2PassTightID,"lep2PassTightID/O");
      tree_->Branch("lep2PassVetoIso",&lep2PassVetoIso,"lep2PassVetoIso/O");
      tree_->Branch("lep2PassLooseIso",&lep2PassLooseIso,"lep2PassLooseIso/O");
      tree_->Branch("lep2PassMediumIso",&lep2PassMediumIso,"lep2PassMediumIso/O");
      tree_->Branch("lep2PassTightIso",&lep2PassTightIso,"lep2PassTightIso/O");
      tree_->Branch("lep2MinDRToBJet",&lep2MinDRToBJet,"lep2MinDRToBJet/F");
      tree_->Branch("lep2Activity",&lep2Activity,"lep2Activity/F");
      tree_->Branch("topPtWeight",&topPtWeight,"topPtWeight/F");
      tree_->Branch("ptTop",&ptTop,"ptTop/F");
      tree_->Branch("ptAntitop",&ptAntitop,"ptAntitop/F");
      tree_->Branch("trigWeight1L",&trigWeight1L,"trigWeight1L/F");
      tree_->Branch("eleEffWeight",&eleEffWeight,"eleEffWeight/F");
      tree_->Branch("muonEffWeight",&muonEffWeight,"muonEffWeight/F");
      tree_->Branch("eleRecoEffWeight",&eleRecoEffWeight,"eleRecoEffWeight/F");
      tree_->Branch("muonRecoEffWeight",&muonRecoEffWeight,"muonRecoEffWeight/F");
      tree_->Branch("bjet1PassLoose",&bjet1PassLoose,"bjet1PassLoose/O");
      tree_->Branch("bjet1PassMedium",&bjet1PassMedium,"bjet1PassMedium/O");
      tree_->Branch("bjet1PassTight",&bjet1PassTight,"bjet1PassTight/O");
      tree_->Branch("bjet2PassLoose",&bjet2PassLoose,"bjet2PassLoose/O");
      tree_->Branch("bjet2PassMedium",&bjet2PassMedium,"bjet2PassMedium/O");
      tree_->Branch("bjet2PassTight",&bjet2PassTight,"bjet2PassTight/O");
      tree_->Branch("jet1PassCSVLoose",&jet1PassCSVLoose,"jet1PassCSVLoose/O");
      tree_->Branch("jet1PassCSVMedium",&jet1PassCSVMedium,"jet1PassCSVMedium/O");
      tree_->Branch("jet1PassCSVTight",&jet1PassCSVTight,"jet1PassCSVTight/O");
      tree_->Branch("jet2PassCSVLoose",&jet2PassCSVLoose,"jet2PassCSVLoose/O");
      tree_->Branch("jet2PassCSVMedium",&jet2PassCSVMedium,"jet2PassCSVMedium/O");
      tree_->Branch("jet2PassCSVTight",&jet2PassCSVTight,"jet2PassCSVTight/O");
      tree_->Branch("lep1MT",&lep1MT,"lep1MT/F");
      tree_->Branch("mll",&mll,"mll/F");      
      tree_->Branch("MET",&MET,"MET/F");
      tree_->Branch("METPhi",&METPhi,"METPhi/F");
      tree_->Branch("METnoHF", &METnoHF, "METnoHF/F");
      tree_->Branch("METnoHFPhi", &METnoHFPhi, "METnoHFPhi/F");
      tree_->Branch("METRaw", &METRaw, "METRaw/F");
      tree_->Branch("METRawPhi", &METRawPhi, "METRawPhi/F");
      tree_->Branch("dPhiRazor",&dPhiRazor,"dPhiRazor/F");
      tree_->Branch("HT",&HT,"HT/F");	  
      tree_->Branch("genZpt",&genZpt,"genZpt/F");
      tree_->Branch("genZphi",&genZphi,"genZphi/F");
      tree_->Branch("genlep1", "TLorentzVector", &genlep1Ptr);
      tree_->Branch("genlep2", "TLorentzVector", &genlep2Ptr);
      tree_->Branch("lep1",    "TLorentzVector", &lep1Ptr);
      tree_->Branch("lep2",    "TLorentzVector", &lep2Ptr);
      tree_->Branch("bjet1",   "TLorentzVector", &bjet1Ptr);
      tree_->Branch("bjet2",   "TLorentzVector", &bjet2Ptr);
      tree_->Branch("jet1",    "TLorentzVector", &jet1Ptr);
      tree_->Branch("jet2",    "TLorentzVector", &jet2Ptr);
    }

    if (treeType == kTreeType_OneLeptonAdd2MET_Full ) {
      tree_->Branch("lep1",    "TLorentzVector", &lep1Ptr);
      tree_->Branch("lep1Type", &lep1Type, "lep1Type/I");
      tree_->Branch("lep1PassVeto",&lep1PassVeto,"lep1PassVeto/O");
      tree_->Branch("lep1PassLoose",&lep1PassLoose,"lep1PassLoose/O");
      tree_->Branch("lep1PassMedium",&lep1PassMedium,"lep1PassMedium/O");
      tree_->Branch("lep1PassTight",&lep1PassTight,"lep1PassTight/O");
      tree_->Branch("topPtWeight",&topPtWeight,"topPtWeight/F");
      tree_->Branch("ptTop",&ptTop,"ptTop/F");
      tree_->Branch("ptAntitop",&ptAntitop,"ptAntitop/F");
      tree_->Branch("trigWeight1L",&trigWeight1L,"trigWeight1L/F");
      tree_->Branch("eleEffWeight",&eleEffWeight,"eleEffWeight/F");
      tree_->Branch("muonEffWeight",&muonEffWeight,"muonEffWeight/F");
      tree_->Branch("eleRecoEffWeight",&eleRecoEffWeight,"eleRecoEffWeight/F");
      tree_->Branch("muonRecoEffWeight",&muonRecoEffWeight,"muonRecoEffWeight/F");
      tree_->Branch("MET",&MET,"MET/F");
      tree_->Branch("METPhi",&METPhi,"METPhi/F");
      tree_->Branch("MR_NoW",&MR_NoW,"MR_NoW/F");
      tree_->Branch("Rsq_NoW",&Rsq_NoW,"Rsq_NoW/F");
      tree_->Branch("MET_NoW",&MET_NoW,"MET_NoW/F");
      tree_->Branch("METPhi_NoW",&METPhi_NoW,"METPhi_NoW/F");
      tree_->Branch("HT_NoW",&HT_NoW,"HT_NoW/F");
      tree_->Branch("dPhiRazor_NoW",&dPhiRazor_NoW,"dPhiRazor_NoW/F");
      tree_->Branch("NJets_NoW",&NJets_NoW,"NJets_NoW/i");
      tree_->Branch("NJets80_NoW",&NJets80_NoW,"NJets80_NoW/i");
      tree_->Branch("nVetoMuons",&nVetoMuons,"nVetoMuons/i");
      tree_->Branch("nLooseMuons",&nLooseMuons,"nLooseMuons/i");
      tree_->Branch("nMediumMuons",&nMediumMuons,"nMediumMuons/i");
      tree_->Branch("nTightMuons",&nTightMuons,"nTightMuons/i");
      tree_->Branch("lep1MT",&lep1MT,"lep1MT/F");
      tree_->Branch("HLT_SingleMu",&HLT_SingleMu,"HLT_SingleMu/O");
      tree_->Branch("genWpt",&genWpt,"genWpt/F");
      tree_->Branch("genWphi",&genWphi,"genWphi/F");
      tree_->Branch("recoWpt",&recoWpt,"recoWpt/F");
      tree_->Branch("recoWphi",&recoWphi,"recoWphi/F");
      tree_->Branch("HLTDecision",&HLTDecision,"HLTDecision[300]/O");
    }
  
    if (treeType == kTreeType_DileptonAdd2MET_Full ) {
      tree_->Branch("lep1",    "TLorentzVector", &lep1Ptr);
      tree_->Branch("lep2",    "TLorentzVector", &lep2Ptr);
      tree_->Branch("lep1Type",&lep1Type,"lep1Type/I");
      tree_->Branch("lep2Type",&lep2Type,"lep2Type/I");
      tree_->Branch("lep1PassVeto",&lep1PassVeto,"lep1PassVeto/O");
      tree_->Branch("lep1PassLoose",&lep1PassLoose,"lep1PassLoose/O");
      tree_->Branch("lep1PassMedium",&lep1PassMedium,"lep1PassMedium/O");
      tree_->Branch("lep1PassTight",&lep1PassTight,"lep1PassTight/O");
      tree_->Branch("lep2PassVeto",&lep2PassVeto,"lep2PassVeto/O");
      tree_->Branch("lep2PassLoose",&lep2PassLoose,"lep2PassLoose/O");
      tree_->Branch("lep2PassMedium",&lep2PassMedium,"lep2PassMedium/O");
      tree_->Branch("lep2PassTight",&lep2PassTight,"lep2PassTight/O");
      tree_->Branch("topPtWeight",&topPtWeight,"topPtWeight/F");
      tree_->Branch("ptTop",&ptTop,"ptTop/F");
      tree_->Branch("ptAntitop",&ptAntitop,"ptAntitop/F");
      tree_->Branch("trigWeight1L",&trigWeight1L,"trigWeight1L/F");
      tree_->Branch("eleEffWeight",&eleEffWeight,"eleEffWeight/F");
      tree_->Branch("muonEffWeight",&muonEffWeight,"muonEffWeight/F");
      tree_->Branch("eleRecoEffWeight",&eleRecoEffWeight,"eleRecoEffWeight/F");
      tree_->Branch("muonRecoEffWeight",&muonRecoEffWeight,"muonRecoEffWeight/F");
      tree_->Branch("MET",&MET,"MET/F");
      tree_->Branch("u1",&u1,"u1/F");
      tree_->Branch("u2",&u2,"u2/F");
      tree_->Branch("METPhi",&METPhi,"METPhi/F");
      tree_->Branch("MR_NoZ",&MR_NoZ,"MR_NoZ/F");
      tree_->Branch("Rsq_NoZ",&Rsq_NoZ,"Rsq_NoZ/F");
      tree_->Branch("MET_NoZ",&MET_NoZ,"MET_NoZ/F");
      tree_->Branch("METPhi_NoZ",&METPhi_NoZ,"METPhi_NoZ/F");
      tree_->Branch("HT_NoZ",&HT_NoZ,"HT_NoZ/F");
      tree_->Branch("dPhiRazor_NoZ",&dPhiRazor_NoZ,"dPhiRazor_NoZ/F");
      tree_->Branch("recoZmass",&recoZmass,"recoZmass/F");
      tree_->Branch("recoZpt",&recoZpt,"recoZpt/F");
      tree_->Branch("genZpt",&genZpt,"genZpt/F");
      tree_->Branch("genZphi",&genZphi,"genZphi/F");
      tree_->Branch("NJets_NoZ",&NJets_NoZ,"NJets_NoZ/i");
      tree_->Branch("NJets80_NoZ",&NJets80_NoZ,"NJets80_NoZ/i");
      tree_->Branch("nVetoMuons",&nVetoMuons,"nVetoMuons/i");
      tree_->Branch("nLooseMuons",&nLooseMuons,"nLooseMuons/i");
      tree_->Branch("nMediumMuons",&nMediumMuons,"nMediumMuons/i");
      tree_->Branch("nTightMuons",&nTightMuons,"nTightMuons/i");
      tree_->Branch("HLT_Dimuon",&HLT_Dimuon,"HLT_Dimuon/O");
      tree_->Branch("HLTDecision",&HLTDecision,"HLTDecision[300]/O");
    }
    
    // fill the photon tree
    if (treeType == kTreeType_Photon_Full) {
      tree_->Branch("HLTDecision",&HLTDecision,"HLTDecision[300]/O");
      tree_->Branch("HLTPrescale",&HLTPrescale,"HLTPrescale[300]/I");
      tree_->Branch("pho1HLTFilter",&pho1HLTFilter,"HLTDecision[50]/O");
      tree_->Branch("pho1","TLorentzVector", &pho1Ptr);
      tree_->Branch("pho2","TLorentzVector", &pho2Ptr);
      tree_->Branch("jet1",    "TLorentzVector", &jet1Ptr);
      tree_->Branch("jet2",    "TLorentzVector", &jet2Ptr);
	  
      tree_->Branch("MR_NoPho",&MR_NoPho,"MR_NoPho/F");
      tree_->Branch("Rsq_NoPho",&Rsq_NoPho,"Rsq_NoPho/F");
      tree_->Branch("MET_NoPho",&MET_NoPho,"MET_NoPho/F");
      tree_->Branch("METPhi_NoPho",&METPhi_NoPho,"METPhi_NoPho/F");
      tree_->Branch("HT_NoPho",&HT_NoPho,"HT_NoPho/F");
      tree_->Branch("dPhiRazor_NoPho",&dPhiRazor_NoPho,"dPhiRazor_NoPho/F");
      tree_->Branch("NJets_NoPho",&NJets_NoPho,"NJets_NoPho/i");
      tree_->Branch("NJets80_NoPho",&NJets80_NoPho,"NJets80_NoPho/i");

      tree_->Branch("nSelectedPhotons",&nSelectedPhotons,"nSelectedPhotons/i");
      tree_->Branch("NJets_NoPho",&NJets_NoPho,"NJets_NoPho/i");
      tree_->Branch("NJets80_NoPho",&NJets80_NoPho,"NJets80_NoPho/i");
      tree_->Branch("pho1_motherID",&pho1_motherID,"pho1_motherID/I");
      tree_->Branch("pho1_sigmaietaieta",&pho1_sigmaietaieta,"pho1_sigmaietaieta/F");
      tree_->Branch("pho1_chargediso",&pho1_chargediso,"pho1_chargediso/F");
      tree_->Branch("pho1_photoniso",&pho1_photoniso,"pho1_photoniso/F");
      tree_->Branch("pho1_neutralhadroniso",&pho1_neutralhadroniso,"pho1_neutralhadroniso/F");
      tree_->Branch("pho1_hOverE",&pho1_hOverE,"pho1_hOverE/F");
      tree_->Branch("pho1_pfiso",&pho1_pfiso,"pho1_pfiso/F");
      tree_->Branch("pho1PassTight",&pho1PassTight,"pho1PassTight/O");
      tree_->Branch("minDRGenPhotonToParton",&minDRGenPhotonToParton,"minDRGenPhotonToParton/F");

      tree_->Branch("genJetHT",&genJetHT,"genJetHT/F");
      tree_->Branch("genJetMR",&genJetMR,"genJetMR/F");
    }  
  }
  

  // initialze a ControlSampleEvents
  void InitTree(int treeType = kTreeType_Default){
    assert(tree_);
    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;

    tree_->SetBranchAddress("option",&option);
    tree_->SetBranchAddress("genWeight",&genWeight);
    tree_->SetBranchAddress("weight",&weight);
    tree_->SetBranchAddress("pileupWeight",&pileupWeight);
    tree_->SetBranchAddress("btagW",&btagW);
    tree_->SetBranchAddress("btagW_up",&btagW_up);
    tree_->SetBranchAddress("btagW_down",&btagW_down);
    tree_->SetBranchAddress("run",&run);
    tree_->SetBranchAddress("lumi",&lumi);
    tree_->SetBranchAddress("NPU_0",&NPU_0);
    tree_->SetBranchAddress("NPV",&NPV);
    tree_->SetBranchAddress("MR",&MR);
    tree_->SetBranchAddress("Rsq",&Rsq);
    tree_->SetBranchAddress("NJets40",&NJets40);
    tree_->SetBranchAddress("NJets80",&NJets80);
    tree_->SetBranchAddress("NGenBJets",&NGenBJets);
    tree_->SetBranchAddress("NBJetsLoose",&NBJetsLoose);
    tree_->SetBranchAddress("NBJetsMedium",&NBJetsMedium);
    tree_->SetBranchAddress("NBJetsTight",&NBJetsTight);
    tree_->SetBranchAddress("genHT",&genHT);
    tree_->SetBranchAddress("NISRJets",&NISRJets);

    tree_->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter);
    tree_->SetBranchAddress("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter);
    tree_->SetBranchAddress("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter);
    tree_->SetBranchAddress("Flag_badMuonFilter", &Flag_badMuonFilter);
    tree_->SetBranchAddress("Flag_badGlobalMuonFilter", &Flag_badGlobalMuonFilter);
    tree_->SetBranchAddress("Flag_duplicateMuonFilter", &Flag_duplicateMuonFilter);
    tree_->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter);
    tree_->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter);
    tree_->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter);
    tree_->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices);
    tree_->SetBranchAddress("Flag_trackingFailureFilter", &Flag_trackingFailureFilter);
    tree_->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter);
    tree_->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter);
    tree_->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters);
    tree_->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X);
    tree_->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X);
    tree_->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters);
    tree_->SetBranchAddress("Flag_METFilters", &Flag_METFilters);

    // book the branches that go into only One Lepton trees
    if (treeType == kTreeType_OneLepton_Reduced
	|| treeType == kTreeType_ZeroLepton_Reduced
	|| treeType == kTreeType_OneVetoLepton_Reduced
	|| treeType == kTreeType_OneTauLepton_Reduced
	) {
      tree_->SetBranchStatus("*", 0);
      tree_->SetBranchStatus("lep1Type", 1);
      tree_->SetBranchStatus("lep1PassTight", 1);
      tree_->SetBranchStatus("lep1MT", 1);
      tree_->SetBranchStatus("lep1Pt", 1);
      tree_->SetBranchStatus("lep1Eta", 1);
      tree_->SetBranchStatus("MET", 1);
      tree_->SetBranchStatus("weight", 1);
      tree_->SetBranchStatus("pileupWeight", 1);
      tree_->SetBranchStatus("topPtWeight", 1);
      tree_->SetBranchStatus("ptTop", 1);
      tree_->SetBranchStatus("ptAntitop", 1);
      tree_->SetBranchStatus("trigWeight1L", 1);
      tree_->SetBranchStatus("eleEffWeight", 1);
      tree_->SetBranchStatus("muonEffWeight", 1);
      tree_->SetBranchStatus("eleRecoEffWeight", 1);
      tree_->SetBranchStatus("muonRecoEffWeight", 1);
      tree_->SetBranchStatus("MR", 1);
      tree_->SetBranchStatus("Rsq", 1);
      tree_->SetBranchStatus("NJets40", 1);
      tree_->SetBranchStatus("NJets80", 1);
      tree_->SetBranchStatus("NBJetsLoose", 1);
      tree_->SetBranchStatus("NBJetsMedium", 1);
      tree_->SetBranchStatus("HLTDecision", 1);
      
      tree_->SetBranchAddress("genlep1Type",&genlep1Type);
      tree_->SetBranchAddress("lep1Type",&lep1Type);
      tree_->SetBranchAddress("lep1PassVeto",&lep1PassVeto);
      tree_->SetBranchAddress("lep1PassLoose",&lep1PassLoose);
      tree_->SetBranchAddress("lep1PassMedium",&lep1PassMedium);
      tree_->SetBranchAddress("lep1PassTight",&lep1PassTight);
      tree_->SetBranchAddress("lep1MatchedGenLepIndex",&lep1MatchedGenLepIndex);
      tree_->SetBranchAddress("lep1MT",&lep1MT);	
      tree_->SetBranchAddress("topPtWeight", &topPtWeight);
      tree_->SetBranchAddress("ptTop", &ptTop);
      tree_->SetBranchAddress("ptAntitop", &ptAntitop);
      tree_->SetBranchAddress("trigWeight1L", &trigWeight1L);
      tree_->SetBranchAddress("eleEffWeight", &eleEffWeight);
      tree_->SetBranchAddress("muonEffWeight", &muonEffWeight);
      tree_->SetBranchAddress("eleRecoEffWeight", &eleRecoEffWeight);
      tree_->SetBranchAddress("muonRecoEffWeight", &muonRecoEffWeight);
      tree_->SetBranchAddress("MET",&MET);
      tree_->SetBranchAddress("METPhi",&METPhi);
      tree_->SetBranchAddress("METnoHF",&METnoHF);
      tree_->SetBranchAddress("HT",&HT);
      tree_->SetBranchAddress("HLTDecision",&HLTDecision);
      tree_->SetBranchAddress("lep1Pt",&lep1Pt);
      tree_->SetBranchAddress("lep1Eta",&lep1Eta);
      tree_->SetBranchAddress("MHT",&MHT);
      tree_->SetBranchAddress("MHTnoHF",&MHTnoHF);
      tree_->SetBranchAddress("lep1" ,&lep1Ptr);
   }

    if (treeType == kTreeType_OneLepton_Full
	|| treeType == kTreeType_ZeroLepton_Full	
	|| treeType == kTreeType_OneVetoLepton_Full
	|| treeType == kTreeType_OneTauLepton_Full
	) {
      tree_->SetBranchAddress("NPU_Minus1",&NPU_Minus1);
      tree_->SetBranchAddress("NPU_Plus1",&NPU_Plus1);
      tree_->SetBranchAddress("event",&event);
      tree_->SetBranchAddress("processID",&processID);
      tree_->SetBranchAddress("HLTDecision",&HLTDecision);
      tree_->SetBranchAddress("genlep1Type",&genlep1Type);
      tree_->SetBranchAddress("lep1Type",&lep1Type);
      tree_->SetBranchAddress("lep1PassVeto",&lep1PassVeto);
      tree_->SetBranchAddress("lep1PassLoose",&lep1PassLoose);
      tree_->SetBranchAddress("lep1PassMedium",&lep1PassMedium);
      tree_->SetBranchAddress("lep1PassTight",&lep1PassTight);
      tree_->SetBranchAddress("lep1MatchedGenLepIndex",&lep1MatchedGenLepIndex);
      tree_->SetBranchAddress("lep1MT",&lep1MT);
      tree_->SetBranchAddress("lep1Activity",&lep1Activity);
      tree_->SetBranchAddress("topPtWeight", &topPtWeight);
      tree_->SetBranchAddress("ptTop", &ptTop);
      tree_->SetBranchAddress("ptAntitop", &ptAntitop);
      tree_->SetBranchAddress("trigWeight1L", &trigWeight1L);
      tree_->SetBranchAddress("eleEffWeight", &eleEffWeight);
      tree_->SetBranchAddress("muonEffWeight", &muonEffWeight);
      tree_->SetBranchAddress("eleRecoEffWeight", &eleRecoEffWeight);
      tree_->SetBranchAddress("muonRecoEffWeight", &muonRecoEffWeight);
      tree_->SetBranchAddress("MET",&MET);
      tree_->SetBranchAddress("METPhi",&METPhi);
      tree_->SetBranchAddress("METnoHF",&METnoHF);
      tree_->SetBranchAddress("METnoHFPhi",&METnoHFPhi);
      tree_->SetBranchAddress("METRaw",&METRaw);
      tree_->SetBranchAddress("METRawPhi",&METRawPhi);
      tree_->SetBranchAddress("bjet1PassLoose",&bjet1PassLoose);
      tree_->SetBranchAddress("bjet1PassMedium",&bjet1PassMedium);
      tree_->SetBranchAddress("bjet1PassTight",&bjet1PassTight);
      tree_->SetBranchAddress("bjet2PassLoose",&bjet2PassLoose);
      tree_->SetBranchAddress("bjet2PassMedium",&bjet2PassMedium);
      tree_->SetBranchAddress("bjet2PassTight",&bjet2PassTight);
      tree_->SetBranchAddress("jet1PassCSVLoose",&jet1PassCSVLoose);
      tree_->SetBranchAddress("jet1PassCSVMedium",&jet1PassCSVMedium);
      tree_->SetBranchAddress("jet1PassCSVTight",&jet1PassCSVTight);
      tree_->SetBranchAddress("jet2PassCSVLoose",&jet2PassCSVLoose);
      tree_->SetBranchAddress("jet2PassCSVMedium",&jet2PassCSVMedium);
      tree_->SetBranchAddress("jet2PassCSVTight",&jet2PassCSVTight);
      tree_->SetBranchAddress("dPhiRazor",&dPhiRazor);
      tree_->SetBranchAddress("HT",&HT);
      tree_->SetBranchAddress("genWpt",&genWpt);
      tree_->SetBranchAddress("genWphi",&genWphi);
      tree_->SetBranchAddress("genWpt",  &genWpt);
      tree_->SetBranchAddress("genWphi",  &genWphi);
      tree_->SetBranchAddress("genlep1",  &genlep1Ptr);
      tree_->SetBranchAddress("lep1",     &lep1Ptr);
      tree_->SetBranchAddress("bjet1",    &bjet1Ptr);
      tree_->SetBranchAddress("bjet2",    &bjet2Ptr);
      tree_->SetBranchAddress("jet1",     &jet1Ptr);
      tree_->SetBranchAddress("jet2",     &jet2Ptr);
    }
  
    if (treeType == kTreeType_Dilepton_Full
	|| treeType == kTreeType_TightPlusVetoLepton_Full
	) {
      tree_->SetBranchAddress("NPU_Minus1",&NPU_Minus1);
      tree_->SetBranchAddress("NPU_Plus1",&NPU_Plus1);
      tree_->SetBranchAddress("event",&event);
      tree_->SetBranchAddress("processID",&processID);
      tree_->SetBranchAddress("HLTDecision",&HLTDecision);
      tree_->SetBranchAddress("genlep1Type",&genlep1Type);
      tree_->SetBranchAddress("lep1Type",&lep1Type);
      tree_->SetBranchAddress("lep1MatchedGenLepIndex",&lep1MatchedGenLepIndex);
      tree_->SetBranchAddress("genlep2Type",&genlep2Type);
      tree_->SetBranchAddress("lep2Type",&lep2Type);
      tree_->SetBranchAddress("lep2MatchedGenLepIndex",&lep2MatchedGenLepIndex);
      tree_->SetBranchAddress("lep1PassVeto",&lep1PassVeto);
      tree_->SetBranchAddress("lep1PassLoose",&lep1PassLoose);
      tree_->SetBranchAddress("lep1PassMedium",&lep1PassMedium);
      tree_->SetBranchAddress("lep1PassTight",&lep1PassTight);
      tree_->SetBranchAddress("lep1PassVetoID",&lep1PassVetoID);
      tree_->SetBranchAddress("lep1PassLooseID",&lep1PassLooseID);
      tree_->SetBranchAddress("lep1PassMediumID",&lep1PassMediumID);
      tree_->SetBranchAddress("lep1PassTightID",&lep1PassTightID);
      tree_->SetBranchAddress("lep1PassVetoIso",&lep1PassVetoIso);
      tree_->SetBranchAddress("lep1PassLooseIso",&lep1PassLooseIso);
      tree_->SetBranchAddress("lep1PassMediumIso",&lep1PassMediumIso);
      tree_->SetBranchAddress("lep1PassTightIso",&lep1PassTightIso);
      tree_->SetBranchAddress("lep1MinDRToBJet",&lep1MinDRToBJet);
      tree_->SetBranchAddress("lep1Activity",&lep1Activity);
      tree_->SetBranchAddress("lep2PassVeto",&lep2PassVeto);
      tree_->SetBranchAddress("lep2PassLoose",&lep2PassLoose);
      tree_->SetBranchAddress("lep2PassMedium",&lep2PassMedium);
      tree_->SetBranchAddress("lep2PassTight",&lep2PassTight);
      tree_->SetBranchAddress("lep2PassVetoID",&lep2PassVetoID);
      tree_->SetBranchAddress("lep2PassLooseID",&lep2PassLooseID);
      tree_->SetBranchAddress("lep2PassMediumID",&lep2PassMediumID);
      tree_->SetBranchAddress("lep2PassTightID",&lep2PassTightID);
      tree_->SetBranchAddress("lep2PassVetoIso",&lep2PassVetoIso);
      tree_->SetBranchAddress("lep2PassLooseIso",&lep2PassLooseIso);
      tree_->SetBranchAddress("lep2PassMediumIso",&lep2PassMediumIso);
      tree_->SetBranchAddress("lep2PassTightIso",&lep2PassTightIso);
      tree_->SetBranchAddress("lep2MinDRToBJet",&lep2MinDRToBJet);
      tree_->SetBranchAddress("lep2Activity",&lep2Activity);
      tree_->SetBranchAddress("topPtWeight", &topPtWeight);
      tree_->SetBranchAddress("ptTop", &ptTop);
      tree_->SetBranchAddress("ptAntitop", &ptAntitop);
      tree_->SetBranchAddress("trigWeight1L", &trigWeight1L);
      tree_->SetBranchAddress("eleEffWeight", &eleEffWeight);
      tree_->SetBranchAddress("muonEffWeight", &muonEffWeight);
      tree_->SetBranchAddress("eleRecoEffWeight", &eleRecoEffWeight);
      tree_->SetBranchAddress("muonRecoEffWeight", &muonRecoEffWeight);
      tree_->SetBranchAddress("bjet1PassLoose",&bjet1PassLoose);
      tree_->SetBranchAddress("bjet1PassMedium",&bjet1PassMedium);
      tree_->SetBranchAddress("bjet1PassTight",&bjet1PassTight);
      tree_->SetBranchAddress("bjet2PassLoose",&bjet2PassLoose);
      tree_->SetBranchAddress("bjet2PassMedium",&bjet2PassMedium);
      tree_->SetBranchAddress("bjet2PassTight",&bjet2PassTight);
      tree_->SetBranchAddress("jet1PassCSVLoose",&jet1PassCSVLoose);
      tree_->SetBranchAddress("jet1PassCSVMedium",&jet1PassCSVMedium);
      tree_->SetBranchAddress("jet1PassCSVTight",&jet1PassCSVTight);
      tree_->SetBranchAddress("jet2PassCSVLoose",&jet2PassCSVLoose);
      tree_->SetBranchAddress("jet2PassCSVMedium",&jet2PassCSVMedium);
      tree_->SetBranchAddress("jet2PassCSVTight",&jet2PassCSVTight);
      tree_->SetBranchAddress("lep1MT",&lep1MT);
      tree_->SetBranchAddress("mll",&mll);
      tree_->SetBranchAddress("MET",&MET);
      tree_->SetBranchAddress("METPhi",&METPhi);
      tree_->SetBranchAddress("METnoHF",&METnoHF);
      tree_->SetBranchAddress("METnoHFPhi",&METnoHFPhi);
      tree_->SetBranchAddress("METRaw",&METRaw);
      tree_->SetBranchAddress("METRawPhi",&METRawPhi);
      tree_->SetBranchAddress("dPhiRazor",&dPhiRazor);
      tree_->SetBranchAddress("HT",&HT);
      tree_->SetBranchAddress("genZpt",&genZpt);
      tree_->SetBranchAddress("genZphi",&genZphi);
      tree_->SetBranchAddress("genlep1",  &genlep1Ptr);
      tree_->SetBranchAddress("genlep2",  &genlep2Ptr);
      tree_->SetBranchAddress("lep1",     &lep1Ptr);
      tree_->SetBranchAddress("lep2",     &lep2Ptr);
      tree_->SetBranchAddress("bjet1",    &bjet1Ptr);
      tree_->SetBranchAddress("bjet2",    &bjet2Ptr);
      tree_->SetBranchAddress("jet1",     &jet1Ptr);
      tree_->SetBranchAddress("jet2",     &jet2Ptr);
    }



    if (treeType == kTreeType_OneLeptonAdd2MET_Full ) {
      tree_->SetBranchAddress("lep1",&lep1Ptr);
      tree_->SetBranchAddress("lep1Type",&lep1Type);
      tree_->SetBranchAddress("topPtWeight", &topPtWeight);
      tree_->SetBranchAddress("ptTop", &ptTop);
      tree_->SetBranchAddress("ptAntitop", &ptAntitop);
      tree_->SetBranchAddress("trigWeight1L", &trigWeight1L);
      tree_->SetBranchAddress("eleEffWeight", &eleEffWeight);
      tree_->SetBranchAddress("muonEffWeight", &muonEffWeight);
      tree_->SetBranchAddress("eleRecoEffWeight", &eleRecoEffWeight);
      tree_->SetBranchAddress("muonRecoEffWeight", &muonRecoEffWeight);
      tree_->SetBranchAddress("MR_NoW", &MR_NoW);
      tree_->SetBranchAddress("Rsq_NoW",&Rsq_NoW);
      tree_->SetBranchAddress("MET",&MET);
      tree_->SetBranchAddress("METPhi",&METPhi);
      tree_->SetBranchAddress("MET_NoW",&MET_NoW);
      tree_->SetBranchAddress("METPhi_NoW",&METPhi_NoW);
      tree_->SetBranchAddress("HT_NoW",&HT_NoW);
      tree_->SetBranchAddress("dPhiRazor_NoW",&dPhiRazor_NoW);
      tree_->SetBranchAddress("NJets_NoW",&NJets_NoW);
      tree_->SetBranchAddress("NJets80_NoW",&NJets80_NoW);
      tree_->SetBranchAddress("nVetoMuons",&nVetoMuons);
      tree_->SetBranchAddress("nLooseMuons",&nLooseMuons);
      tree_->SetBranchAddress("nMediumMuons",&nMediumMuons);
      tree_->SetBranchAddress("nTightMuons",&nTightMuons);
      tree_->SetBranchAddress("lep1MT",&lep1MT);
      tree_->SetBranchAddress("HLT_SingleMu",&HLT_SingleMu);
      tree_->SetBranchAddress("genWpt",&genWpt);
      tree_->SetBranchAddress("genWphi",&genWphi);
      tree_->SetBranchAddress("recoWpt",&recoWpt);
      tree_->SetBranchAddress("recoWphi",&recoWphi);
      tree_->SetBranchAddress("HLTDecision",&HLTDecision);
    }

    if (treeType == kTreeType_DileptonAdd2MET_Full ) {
      tree_->SetBranchAddress("lep1",          &lep1Ptr);
      tree_->SetBranchAddress("lep2",          &lep2Ptr);
      tree_->SetBranchAddress("lep1Type",&lep1Type);
      tree_->SetBranchAddress("lep2Type",&lep2Type);
      tree_->SetBranchAddress("topPtWeight", &topPtWeight);
      tree_->SetBranchAddress("ptTop", &ptTop);
      tree_->SetBranchAddress("ptAntitop", &ptAntitop);
      tree_->SetBranchAddress("trigWeight1L", &trigWeight1L);
      tree_->SetBranchAddress("eleEffWeight", &eleEffWeight);
      tree_->SetBranchAddress("muonEffWeight", &muonEffWeight);
      tree_->SetBranchAddress("eleRecoEffWeight", &eleRecoEffWeight);
      tree_->SetBranchAddress("muonRecoEffWeight", &muonRecoEffWeight);
      tree_->SetBranchAddress("MR_NoZ",        &MR_NoZ);
      tree_->SetBranchAddress("Rsq_NoZ",       &Rsq_NoZ);
      tree_->SetBranchAddress("MET",	       &MET);
      tree_->SetBranchAddress("MET_NoZ",       &MET_NoZ);
      tree_->SetBranchAddress("METPhi_NoZ",    &METPhi_NoZ);
      tree_->SetBranchAddress("u1",         &u1);
      tree_->SetBranchAddress("u2",         &u2);
      tree_->SetBranchAddress("METPhi",         &METPhi);
      tree_->SetBranchAddress("HT_NoZ",        &HT_NoZ);
      tree_->SetBranchAddress("dPhiRazor_NoZ", &dPhiRazor_NoZ);
      tree_->SetBranchAddress("recoZmass",     &recoZmass);
      tree_->SetBranchAddress("recoZpt",       &recoZpt);
      tree_->SetBranchAddress("genZpt",        &genZpt);
      tree_->SetBranchAddress("genZphi",       &genZphi);
      tree_->SetBranchAddress("NJets_NoZ",     &NJets_NoZ);
      tree_->SetBranchAddress("NJets80_NoZ",   &NJets80_NoZ);
      tree_->SetBranchAddress("nVetoMuons",    &nVetoMuons);
      tree_->SetBranchAddress("nLooseMuons",   &nLooseMuons);
      tree_->SetBranchAddress("nMediumMuons",  &nMediumMuons);
      tree_->SetBranchAddress("nTightMuons",   &nTightMuons);
      tree_->SetBranchAddress("HLT_Dimuon",    &HLT_Dimuon);
      tree_->SetBranchAddress("HLTDecision",&HLTDecision);
    }
    
    // fill the photon tree
    if (treeType == kTreeType_Photon_Full) {
      tree_->SetBranchAddress("HLTDecision",&HLTDecision);
      tree_->SetBranchAddress("HLTPrescale",&HLTPrescale);
      tree_->SetBranchAddress("pho1HLTFilter",&pho1HLTFilter);
      tree_->SetBranchAddress("pho1", &pho1Ptr);
      tree_->SetBranchAddress("pho2", &pho2Ptr);
      tree_->SetBranchAddress("jet1" ,&jet1Ptr);
      tree_->SetBranchAddress("jet2" ,&jet2Ptr);
	  
      tree_->SetBranchAddress("MR_NoPho",&MR_NoPho);
      tree_->SetBranchAddress("Rsq_NoPho",&Rsq_NoPho);
      tree_->SetBranchAddress("MET_NoPho",&MET_NoPho);
      tree_->SetBranchAddress("METPhi_NoPho",&METPhi_NoPho);
      tree_->SetBranchAddress("HT_NoPho",&HT_NoPho);
      tree_->SetBranchAddress("dPhiRazor_NoPho",&dPhiRazor_NoPho);
      tree_->SetBranchAddress("NJets_NoPho",&NJets_NoPho);
      tree_->SetBranchAddress("NJets80_NoPho",&NJets80_NoPho);

      tree_->SetBranchAddress("nSelectedPhotons",&nSelectedPhotons);
      tree_->SetBranchAddress("NJets_NoPho",&NJets_NoPho);
      tree_->SetBranchAddress("NJets80_NoPho",&NJets80_NoPho);
      tree_->SetBranchAddress("pho1_motherID",&pho1_motherID);
      tree_->SetBranchAddress("pho1_sigmaietaieta",&pho1_sigmaietaieta);
      tree_->SetBranchAddress("pho1_chargediso",&pho1_chargediso);
      tree_->SetBranchAddress("pho1_photoniso",&pho1_photoniso);
      tree_->SetBranchAddress("pho1_neutralhadroniso",&pho1_neutralhadroniso);
      tree_->SetBranchAddress("pho1_hOverE",&pho1_hOverE);
      tree_->SetBranchAddress("pho1_pfiso",&pho1_pfiso);
      tree_->SetBranchAddress("pho1PassTight",&pho1PassTight);
      tree_->SetBranchAddress("minDRGenPhotonToParton",&minDRGenPhotonToParton);
      tree_->SetBranchAddress("genJetMR",&genJetMR);
      tree_->SetBranchAddress("genJetHT",&genJetHT);     
    }
    
    gErrorIgnoreLevel = currentState;
  }

  //check if the current event satisfies the selection criteria for the given control sample
  //Supports the following control regions:
  //TTBarSingleLepton
  //TTBarDilepton
  //WSingleLepton
  //ZLLDilepton
  //ZNuNuDilepton
  //ZNuNuSingleLepton
  //ZNuNuPhoton
  bool inControlSample(string sampleName, string option = "", bool isRunOne = true){
      using namespace std;

      if(sampleName == "TTBarSingleLepton"){
          //HLT
          if(isRunOne){
              bool passedTrigger = HLTDecision[0] || HLTDecision[1] || HLTDecision[8] || HLTDecision[9];
              if(!passedTrigger) return false;
          }

          //lepton = tight ele or mu with pt > 30
          if(abs(lep1Type) != 11 && abs(lep1Type) != 13) return false;
          if(!lep1PassTight) return false;
          if(lep1.Pt() < 30) return false;

          //MET and MT cuts
          if(MET < 30) return false;
          if(lep1MT < 30 || lep1MT > 100) return false;

          //b-tag requirement
          if(option == "TwoLooseBTag"){
              if(NBJetsLoose < 2) return false;
          }
          else{
              if(NBJetsMedium < 1) return false;
          }

          //razor baseline cut
          if(MR < 300 || Rsq < 0.15) return false;

          //passes selection
          return true;
      }
      else if(sampleName == "TTBarDilepton"){
          //HLT
          if(isRunOne){
              bool passedTrigger = HLTDecision[3] || HLTDecision[4] || HLTDecision[12] || HLTDecision[6] || HLTDecision[7];
              if(!passedTrigger) return false;
          }

          //leptons = two loose ele or mu with pt > 25, mass outside Z window
          if(abs(lep1Type) != 11 && abs(lep1Type) != 13) return false;
          if(abs(lep2Type) != 11 && abs(lep2Type) != 13) return false;
          if(!lep1PassLoose) return false;
          if(!lep2PassLoose) return false;
          if(lep1.Pt() < 25) return false;
          if(lep2.Pt() < 25) return false;
          float mLL = (lep1+lep2).M();
          if(mLL < 20) return false;
          if(abs(lep1Type) == abs(lep2Type) && mLL > 76 && mLL < 106) return false;

          //MET cut
          if(MET < 40) return false;
          
          //b-tag requirement
          if(option == "TwoLooseBTag"){
              if(NBJetsLoose < 2) return false;
          }
          else{
              if(NBJetsMedium < 1) return false;
          }

          //razor baseline cut
          if(MR < 300 || Rsq < 0.15) return false;

          //passes selection
          return true;
      }
      else if(sampleName == "WSingleLepton"){
          //HLT
          if(isRunOne){
              bool passedTrigger = HLTDecision[0] || HLTDecision[1] || HLTDecision[8] || HLTDecision[9]; 
              if(!passedTrigger) return false;
          }

          //lepton = tight ele or mu with pt > 30
          if(abs(lep1Type) != 11 && abs(lep1Type) != 13) return false;
          if(!lep1PassTight) return false;
          if(lep1.Pt() < 30) return false;

          //MET and MT cuts
          if(MET < 30) return false;
          if(lep1MT < 30 || lep1MT > 100) return false;

          //b-tag requirement
          if(NBJetsMedium > 0) return false;

          //razor baseline cut
          if(MR < 300 || Rsq < 0.15) return false;

          //passes selection
          return true;
      }
      else if(sampleName == "ZLLDilepton"){
          //HLT
          if(isRunOne){
              bool passedTrigger = HLTDecision[3] || HLTDecision[4] || HLTDecision[12];
              if(!passedTrigger) return false;
          }

          //leptons = two loose ele or mu with pt > 25, mass inside Z window
          if(abs(lep1Type) != 11 && abs(lep1Type) != 13) return false;
          if(abs(lep2Type) != 11 && abs(lep2Type) != 13) return false;
          if(!lep1PassLoose) return false;
          if(!lep2PassLoose) return false;
          if(lep1.Pt() < 25) return false;
          if(lep2.Pt() < 25) return false;
          float mLL = (lep1+lep2).M();
          if(mLL < 80 || mLL > 110) return false;

          //b-tag requirement
          if(NBJetsMedium > 0) return false;

          //razor baseline cut
          if(MR < 300 || Rsq < 0.15) return false;

          //passes selection
          return true;
      }
      else if(sampleName == "ZNuNuDilepton"){
          //HLT
          if(isRunOne){
              bool passedTrigger = HLTDecision[3] || HLTDecision[4] || HLTDecision[12];
              if(!passedTrigger) return false;
          }
          //leptons = two loose ele or mu with pt > 25, mass inside Z window
          //if(abs(lep1Type) != 11 && abs(lep1Type) != 13) return false;
          //if(abs(lep2Type) != 11 && abs(lep2Type) != 13) return false;
          //if(!lep1PassLoose) return false;
          //if(!lep2PassLoose) return false;
          if(lep1.Pt() < 25) return false;
          if(lep2.Pt() < 25) return false;
          float mLL = (lep1+lep2).M();
          if(mLL < 80 || mLL > 110) return false;

          //b-tag requirement
          if(NBJetsMedium > 0) return false; 

          //razor baseline cut
          if(MR_NoZ < 300 || Rsq_NoZ < 0.15) return false;

          //passes selection
          return true;
      }
      else if(sampleName == "ZNuNuSingleLepton"){
          //HLT
          if(isRunOne){
              bool passedTrigger = HLTDecision[0] || HLTDecision[1] || HLTDecision[8] || HLTDecision[9]; 
              if(!passedTrigger) return false;
          }

          //lepton = tight ele or mu with pt > 30
          if(abs(lep1Type) != 11 && abs(lep1Type) != 13) return false;
          if(!lep1PassTight) return false;
          if(lep1.Pt() < 30) return false;

          //MET and MT cuts
          if(MET < 30) return false;
          if(lep1MT < 30 || lep1MT > 100) return false;

          //b-tag requirement
          if(NBJetsMedium > 0) return false;

          //razor baseline cut
          if(MR_NoW < 300 || Rsq_NoW < 0.15) return false;

          //passes selection
          return true;
      }
      else if(sampleName == "ZNuNuPhoton"){
          if(isRunOne){
              //HLT
              bool passedTrigger = HLTDecision[29] || HLTDecision[30] || HLTDecision[31] || HLTDecision[32] || HLTDecision[33] || HLTDecision[34];
              if(!passedTrigger) return false;
          }

          //photon ID
          if(!pho1PassMedium) return false;
          if(pho1.Pt() < 80) return false;

          //razor baseline cut
          if(MR_NoPho < 300 || Rsq_NoPho < 0.15) return false;

          //passes selection
          return true;
      }
      else{
          std::cout << "Warning: control sample " << sampleName << " is not recognized." << std::endl;
      }

      return false;
  }

  //compute b-tagging scale factor
  double getRunOneBTagMediumScaleFactor(){
      double btagScaleFactor = 1.0;
      double bjet1EventScaleFactor = 1.0;
      double bjet2EventScaleFactor = 1.0;
      if (bjet1.Pt() > 20) {
          double bjet1ScaleFactor = 0.938887 + 0.00017124 * bjet1.Pt() + (-2.76366e-07) * bjet1.Pt() * bjet1.Pt() ;
          double MCEff = 1.0;
          if (bjet1.Pt() < 50) MCEff = 0.65;
          else if (bjet1.Pt() < 80) MCEff = 0.70;
          else if (bjet1.Pt() < 120) MCEff = 0.73;
          else if (bjet1.Pt() < 210) MCEff = 0.73;
          else MCEff = 0.66;				 
          if (bjet1PassMedium) bjet1EventScaleFactor = bjet1ScaleFactor;
          else bjet1EventScaleFactor = ( 1/MCEff - bjet1ScaleFactor) / ( 1/MCEff - 1);
      }
      if (bjet2.Pt() > 20) {
          double bjet2ScaleFactor = 0.938887 + 0.00017124 * bjet2.Pt() + (-2.76366e-07) * bjet2.Pt() * bjet2.Pt() ;
          double MCEff = 1.0;
          if (bjet2.Pt() < 50) MCEff = 0.65;
          else if (bjet2.Pt() < 80) MCEff = 0.70;
          else if (bjet2.Pt() < 120) MCEff = 0.73;
          else if (bjet2.Pt() < 210) MCEff = 0.73;
          else MCEff = 0.66;		 
          if (bjet2PassMedium) bjet2EventScaleFactor = bjet2ScaleFactor;
          else bjet2EventScaleFactor = ( 1/MCEff - bjet2ScaleFactor) / ( 1/MCEff - 1);
      }
      btagScaleFactor = bjet1EventScaleFactor * bjet2EventScaleFactor;
      return btagScaleFactor;
  }

  //apply all data/MC corrections
  float getMCCorrection(TH1F *pileupWeightHist, string sampleName, bool isRunOne = true){
      float singleMuTriggerSF = 0.97;
      float singleEleTriggerSF = 0.97;
      float doubleMuTriggerSF = 0.97;
      float doubleMuNormalizationSF = 0.97;

      float corrFactor = 1.0;

      //PU reweighting
      corrFactor *= pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(NPU_0));

      if(isRunOne){
          //btagging scale factor
          if(sampleName == "TTBarSingleLepton" || sampleName == "TTBarDilepton" || sampleName == "WSingleLepton" || sampleName == "ZLLDilepton"){
              corrFactor *= getRunOneBTagMediumScaleFactor();
          }

          //trigger scale factors: single muon
          if((HLTDecision[0] || HLTDecision[1]) && (sampleName == "TTBarSingleLepton" || sampleName == "WSingleLepton" || sampleName == "ZNuNuFromW")){
              corrFactor *= singleMuTriggerSF;
          }
          //trigger scale factors: single ele 
          if((HLTDecision[8] || HLTDecision[9]) && (sampleName == "TTBarSingleLepton" || sampleName == "WSingleLepton" || sampleName == "ZNuNuFromW")){
              corrFactor *= singleEleTriggerSF;
          }
          //trigger scale factors: double muon
          else if((HLTDecision[3] || HLTDecision[4]) && (sampleName == "TTBarDilepton" || sampleName == "ZLLDilepton" || sampleName == "ZNuNuDilepton" || sampleName == "ZLLDilepton")){
              corrFactor *= doubleMuTriggerSF;
              corrFactor *= doubleMuNormalizationSF;
          }

          //TODO: muon and electron ID scale factors
      }

      return corrFactor;
  }

 private:
  TLorentzVector* genlep1Ptr;
  TLorentzVector* genlep2Ptr;
  TLorentzVector* lep1Ptr;
  TLorentzVector* lep2Ptr;
  TLorentzVector* bjet1Ptr;
  TLorentzVector* bjet2Ptr;
  TLorentzVector* jet1Ptr;
  TLorentzVector* jet2Ptr;
  TLorentzVector* pho1Ptr;
  TLorentzVector* pho2Ptr;

}; 


#endif

