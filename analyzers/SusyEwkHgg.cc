//LOCAL INCLUDES
#include "SusyEwkHgg.h"
#include "RazorHelper.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
#include "EnergyScaleCorrection_class_2017.hh"
//C++ INCLUDES
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h>
//ROOT INCLUDES
#include <TH1F.h>
#include <TH2D.h>
#include "TRandom3.h"

using namespace std;

enum SusyEwkHggBox {
  Zmm = 0,
  Zee = 1,
  Emu = 2,
  OneMu = 3,
  OneEle = 4,
  HighPt = 5,
  Hbb = 6,
  Zbb = 7,
  HighRes = 8,
  LowRes = 9,
  None = 10
};

struct PhotonCandidate
{
  int   Index;
  TLorentzVector photon;
  TLorentzVector photonSC;
  float scEta;
  float scPhi;
  float SigmaIetaIeta;
  float R9;
  float HoverE;
  float sumChargedHadronPt;
  float sumNeutralHadronEt;
  float sumPhotonEt;
  float sigmaEOverE;
  bool  _passEleVeto;
  bool  _passIso;
  float scale;
};

struct MuonCandidate
{
  int   Index;
  TLorentzVector muon;
  int muonCharge;
  int isTightMuon;
};

struct ElectronCandidate
{
  int   Index;
  TLorentzVector electron;
  int eleCharge;
  int isTightElectron;
};

struct JetCandidate
{
  int   Index;
  TLorentzVector jet;
  bool isCSVL;
  bool isCSVM;
  bool isCSVT;
};

struct BjetCandidate
{
  int   Index;
  TLorentzVector bjet;
  bool isCSVL;
  bool isCSVM;
  bool isCSVT;
};

struct evt
{
  std::string run;
  std::string event;
};

#define _phodebug 0
#define _debug    0
#define _metdebug 1
#define _ecaldebug 0
#define _jetdebug 0
#define _info     1

const double EB_R = 129.0;
const double EE_Z = 317.0;

const double JET_CUT = 30.;
const double BJET_CUT = 20.;
const int NUM_PDF_WEIGHTS = 60;

//Testing branching and merging
void SusyEwkHgg::Analyze(bool isData, int option, string outFileName, string label)
{
  //*****************************************************************************
  //Settings
  //*****************************************************************************
  TRandom3 random(3003);
  bool doPhotonScaleCorrection = true;
  string analysisTag = "Razor2016_MoriondRereco";
  //string analysisTag = "Razor2016_80X";
  if ( label != "") analysisTag = label;

  string dataset = "80X";
  if ( analysisTag == "Razor2016_MoriondRereco" || analysisTag == "Razor2016_03Feb2017Rereco" ) dataset = "80X";
  else if ( analysisTag == "Razor2017_17Nov2017Rereco" || analysisTag == "Razor2017_31Mar2018Rereco" ) dataset = "94X";

  //***************************************************
  //What the options means:
  //option 0 / 10: Regular Selection, No MR Skim
  //option 1 / 11: Regular Selection, MR>150 Skim
  //option 2 / 12: NoEleVeto, No MR Skim
  //option 3 / 13: NoPhotonID, No MR Skim
  //option 4 / 14: NoPhotonIso, No MR Skim
  //option 5 / 15: NoPhotonID, NoPhotonIso, No MR Skim
  //option 6 / 16: NoPhotonID, NoPhotonIso, NoEleVeto, No MR Skim
  //option 7 / 17: TightPhotonID, LooseIso
  //option 7 / 17: LooseID, TightIso
  //option 7 / 17: TightID , TightIso

  bool doMRSkim = false;
  bool doRequireIso = true;
  bool doRequireID = true;
  bool doEleVeto = true;
  bool doRequireTightID = false;
  bool doRequireTightIso = false;

  if (option == 1 || option == 11) doMRSkim = true;
  if (option == 2 || option == 12) doEleVeto = false;
  if (option == 3 || option == 13) doRequireID = false;
  if (option == 4 || option == 14) doRequireIso = false;
  if (option == 5 || option == 15) { doRequireID = false; doRequireIso = false; }
  if (option == 6 || option == 16) { doRequireID = false;doRequireIso = false; doEleVeto = false; }
  if (option == 7 || option == 17) doRequireTightID = true;
  if (option == 8 || option == 18) doRequireTightIso = true;
  if (option == 9 || option == 19) { doRequireTightID = true; doRequireTightIso = true; }

  bool isFastsimSMS = (option >= 10);
  bool is2DMassScan = (option >= 20);
  std::cout << "[INFO]: option = " << option << std::endl;
  std::cout << "[INFO]: analysisTag --> " << analysisTag << std::endl;
  std::cout << "[INFO]: doRequireID --> " << doRequireID << std::endl;
  std::cout << "[INFO]: doRequireIso --> " << doRequireIso << std::endl;
  std::cout << "[INFO]: doEleVeto --> " << doEleVeto << std::endl;
  std::cout << "[INFO]: doMRSkim --> " << doMRSkim << std::endl;
  std::cout << "[INFO]: doRequireTightID --> " << doRequireTightID << std::endl;
  std::cout << "[INFO]: doRequireTightIso --> " << doRequireTightIso << std::endl;
  std::cout << "[INFO]: doPhotonScaleCorrection --> " << doPhotonScaleCorrection << std::endl;
  std::cout << "[INFO]: isFastsimSMS --> " << isFastsimSMS << std::endl;
  std::cout << "[INFO]: is2DMassScan --> " << is2DMassScan << std::endl;

  //initialization: create one TTree for each analysis box
  if ( _info ) std::cout << "Initializing..." << std::endl;

  if ( outFileName.empty() )
    {
      if ( _info ) std::cout << "SusyEwkHgg: Output filename not specified!" << endl << "Using default output name SusyEwkHgg.root" << std::endl;
      outFileName = "SusyEwkHgg.root";
    }
  TFile* outFile = new TFile( outFileName.c_str(), "RECREATE" );
  //---------------------------
  //one tree to hold all events
  //---------------------------
  TTree *razorTree = new TTree("HggRazorLeptons", "Info on selected razor inclusive events");

  //For signal samples, create one output file and tree per signal mass point
  map<int, TFile*> smsFiles;
  map<int, TTree*> smsTrees;
  map<int, TH1F*> smsNEvents;
  map<int, TH1F*> smsSumWeights;
  map<int, TH1F*> smsSumScaleWeights;
  map<int, TH1F*> smsSumPdfWeights;
  map<int, TH1F*> smsNISRJets;
  map<int, TH1F*> smsPtISR;
  map<int, TH1F*> smsNPV;

  map<pair<int,int>, TFile*> smsFiles2D;
  map<pair<int,int>, TTree*> smsTrees2D;
  map<pair<int,int>, TH1F*> smsNEvents2D;
  map<pair<int,int>, TH1F*> smsSumWeights2D;
  map<pair<int,int>, TH1F*> smsSumScaleWeights2D;
  map<pair<int,int>, TH1F*> smsSumPdfWeights2D;
  map<pair<int,int>, TH1F*> smsNISRJets2D;
  map<pair<int,int>, TH1F*> smsPtISR2D;
  map<pair<int,int>, TH1F*> smsNPV2D;


  //Get CMSSW Directory
  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");

  //--------------------------------
  //Initialize helper
  //--------------------------------
  /* Available tag relevant for this analysis are:
    Razor2015_76X
    Razor2016_MoriondRereco
    Razor2016_03Feb2017Rereco
    Razor2017_92X
    Razor2017_17Nov2017Rereco
    Razor2017_31Mar2018Rereco
  */
  RazorHelper *helper =  new RazorHelper(analysisTag, isData, isFastsimSMS);;
  cout << "Done Loading Helper\n";

  //--------------------------------
  //Photon Energy Scale and Resolution Corrections
  //--------------------------------
  std::string photonCorrectionPath = "./";
   if ( cmsswPath != NULL ) photonCorrectionPath = string(cmsswPath) + "/src/RazorAnalyzer/data/PhotonCorrections/";

  EnergyScaleCorrection_class *photonCorrector = 0;
  EnergyScaleCorrection_class_2017 *photonCorrector_2017 = 0;

  if (analysisTag == "Razor2015_76X") {
    photonCorrector = new EnergyScaleCorrection_class(Form("%s/76X_16DecRereco_2015", photonCorrectionPath.c_str()));
  } else if (analysisTag == "Razor2016_80X") {
    photonCorrector = new EnergyScaleCorrection_class(Form("%s/80X_2016", photonCorrectionPath.c_str()));
  } else if (analysisTag == "Razor2016_MoriondRereco") {
    photonCorrector = new EnergyScaleCorrection_class(Form("%s/Winter_2016_reReco_v1_ele", photonCorrectionPath.c_str()));
  } else if (analysisTag == "Razor2016_03Feb2017Rereco") {
    photonCorrector_2017 = new EnergyScaleCorrection_class_2017(Form("%s/Legacy2016_07Aug2017_FineEtaR9_v3_ele_unc", photonCorrectionPath.c_str()));
  } else if (analysisTag == "Razor2017_92X") {
    photonCorrector_2017 = new EnergyScaleCorrection_class_2017(Form("%s/Run2017_17Nov2017_v1_ele_unc", photonCorrectionPath.c_str()));
  } else if (analysisTag == "Razor2017_17Nov2017Rereco") {
    photonCorrector_2017 = new EnergyScaleCorrection_class_2017(Form("%s/Run2017_17Nov2017_v1_ele_unc", photonCorrectionPath.c_str()));
  } else if (analysisTag == "Razor2017_31Mar2018Rereco") {
    photonCorrector_2017 = new EnergyScaleCorrection_class_2017(Form("%s/Run2017_17Nov2017_v1_ele_unc", photonCorrectionPath.c_str()));
  }


  if ( analysisTag != "Razor2017_92X" && analysisTag != "Razor2017_17Nov2017Rereco" && analysisTag != "Razor2017_31Mar2018Rereco" &&  analysisTag != "Razor2016_03Feb2017Rereco"  ) {
    if(!isData) {
      photonCorrector->doScale = false;
      photonCorrector->doSmearings = true;
    } else {
      photonCorrector->doScale = true;
      photonCorrector->doSmearings = false;
    }
  } else {
    if(!isData) {
      photonCorrector_2017->doScale = false;
      photonCorrector_2017->doSmearings = true;
    } else {
      photonCorrector_2017->doScale = true;
      photonCorrector_2017->doSmearings = false;
    }
  }
  
  //--------------------------------
  //Including Jet Energy Corrections
  //--------------------------------
  std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
  std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();

  //---------------
  //btag efficiency
  //---------------
  //Medium
  TH2D *btagMediumEfficiencyHist = 0;
  TH2D *btagMediumCharmEfficiencyHist = 0;
  TH2D *btagMediumLightJetsEfficiencyHist = 0;
  //Loose
  TH2D *btagLooseEfficiencyHist = 0;
  TH2D *btagLooseCharmEfficiencyHist = 0;
  TH2D *btagLooseLightJetsEfficiencyHist = 0;
  if ( !isData )
  {
    //Medium
    TFile *btagEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_BJets_25ns_CSVM_Fullsim_80X.root");
    btagMediumEfficiencyHist = (TH2D*)btagEfficiencyFile->Get("Efficiency_PtEta");
    assert(btagMediumEfficiencyHist);
    TFile *btagCharmEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_CJets_25ns_CSVM_Fullsim_80X.root");
    btagMediumCharmEfficiencyHist = (TH2D*)btagCharmEfficiencyFile->Get("Efficiency_PtEta");
    assert(btagMediumCharmEfficiencyHist);
    TFile *btagLightJetsEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_LightJets_25ns_CSVM_Fullsim_80X.root");
    btagMediumLightJetsEfficiencyHist = (TH2D*)btagLightJetsEfficiencyFile->Get("Efficiency_PtEta");
    assert(btagMediumLightJetsEfficiencyHist);
    //Loose
    TFile *btagLooseEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_BJets_25ns_CSVL_Fullsim_80X.root");
    btagLooseEfficiencyHist = (TH2D*)btagLooseEfficiencyFile->Get("Efficiency_PtEta");
    assert(btagLooseEfficiencyHist);
    TFile *btagLooseCharmEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_CJets_25ns_CSVL_Fullsim_80X.root");
    btagLooseCharmEfficiencyHist = (TH2D*)btagLooseCharmEfficiencyFile->Get("Efficiency_PtEta");
    assert(btagLooseCharmEfficiencyHist);
    TFile *btagLooseLightJetsEfficiencyFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_LightJets_25ns_CSVL_Fullsim_80X.root");
    btagLooseLightJetsEfficiencyHist = (TH2D*)btagLooseLightJetsEfficiencyFile->Get("Efficiency_PtEta");
    assert(btagLooseLightJetsEfficiencyHist);
  }

  //-----------------------
  //B-tagging scale factors
  //-----------------------

  string bTagPathname = "./";
  // if ( cmsswPath != NULL ) bTagPathname = string(cmsswPath) + "/src/RazorAnalyzer/data/ScaleFactors/";
  // else bTagPathname = "data/ScaleFactors/";
  //Fullsim

  TString effMeasType, misMeasType;
  BTagCalibration *btagcalib =0;
  if (analysisTag == "Razor2015_76X") {
    btagcalib = new BTagCalibration("csvv2", Form("%s/CSVv2_76X.csv",bTagPathname.c_str()));
    effMeasType="mujets";
    misMeasType="comb";
  } else if (analysisTag == "Razor2016_80X" || analysisTag == "Razor2016_MoriondRereco" || analysisTag == "Razor2016_03Feb2017Rereco" || analysisTag == "Razor2017_92X" ) {
    if(isFastsimSMS) {
      btagcalib = new BTagCalibration("csvv2", Form("%s/fastsim_csvv2_ttbar_26_1_2017.csv",bTagPathname.c_str()));
      effMeasType="fastsim";
      misMeasType="fastsim";
    } else {
      btagcalib = new BTagCalibration("csvv2", Form("%s/CSVv2_Moriond17_B_H.csv",bTagPathname.c_str()));
      effMeasType="comb";
      misMeasType="incl";
    }
  } else if ( analysisTag == "Razor2017_17Nov2017Rereco" || analysisTag == "Razor2017_31Mar2018Rereco") {
    btagcalib = new BTagCalibration("csvv2", Form("%s/CSVv2_94XSF_V2_B_F.csv",bTagPathname.c_str()));
    effMeasType="comb";
    misMeasType="incl";
  }

  //Medium WP (b-tagging workpoint)
  BTagCalibrationReader btagreaderM(btagcalib,           //calibration instance
				    BTagEntry::OP_MEDIUM, //operating point
				    effMeasType.Data(),             //measurement type
				    "central");           //systematics type
  BTagCalibrationReader btagreaderM_up(btagcalib, BTagEntry::OP_MEDIUM, effMeasType.Data(), "up");   //sys up
  BTagCalibrationReader btagreaderM_do(btagcalib, BTagEntry::OP_MEDIUM, effMeasType.Data(), "down"); //sys down
  BTagCalibrationReader btagreaderMistagM(btagcalib,            //calibration instance
					  BTagEntry::OP_MEDIUM,  //operating point
					  misMeasType.Data(),                //measurement type
					  "central");            //systematics type
  //BTagCalibrationReader btagreaderMistagM_up(&btagcalib, BTagEntry::OP_MEDIUM, "comb", "up");    //sys up
  BTagCalibrationReader btagreaderMistagM_up(btagcalib, BTagEntry::OP_MEDIUM, misMeasType.Data(), "up");    //sys up
  //BTagCalibrationReader btagreaderMistagM_do(&btagcalib, BTagEntry::OP_MEDIUM, "comb", "down");  //sys down
  BTagCalibrationReader btagreaderMistagM_do(btagcalib, BTagEntry::OP_MEDIUM, misMeasType.Data(), "down");  //sys down
  //Loose WP
  BTagCalibrationReader btagreaderL(btagcalib,           //calibration instance
				    BTagEntry::OP_LOOSE,  //operating point
				    effMeasType.Data(),             //measurement type
				    "central");           //systematics type
  BTagCalibrationReader btagreaderL_up(btagcalib, BTagEntry::OP_LOOSE, effMeasType.Data(), "up");  //sys up
  BTagCalibrationReader btagreaderL_do(btagcalib, BTagEntry::OP_LOOSE, effMeasType.Data(), "down");  //sys down
  BTagCalibrationReader btagreaderMistagL(btagcalib,           //calibration instance
					  BTagEntry::OP_LOOSE,  //operating point
					  //"comb",               //measurement type
					  misMeasType.Data(),               //measurement type
					  "central");           //systematics type
  //BTagCalibrationReader btagreaderMistagL_up(&btagcalib, BTagEntry::OP_LOOSE, "comb", "up");    //sys up
  BTagCalibrationReader btagreaderMistagL_up(btagcalib, BTagEntry::OP_LOOSE, misMeasType.Data(), "up");    //sys up
  //BTagCalibrationReader btagreaderMistagL_do(&btagcalib, BTagEntry::OP_LOOSE, "comb", "down");  //sys down
  BTagCalibrationReader btagreaderMistagL_do(btagcalib, BTagEntry::OP_LOOSE, misMeasType.Data(), "down");  //sys down

  //----------
  //pu histo
  //----------
  TH1D* puhisto = new TH1D("pileup", "", 50, 0, 50);

  //histogram containing total number of processed events (for normalization)
  TH1F *histNPV = new TH1F("NPV", "NPV", 2, -0.5, 1.5);
  TH1F *histNISRJets = new TH1F("NISRJets", "NISRJets", 7, -0.5, 6.5);
  float PtISRBins[9] = {0,50, 100,150,200,300,400,600,7000};
  TH1F *histPtISR = new TH1F("PtISR", "PtISR", 8, PtISRBins );
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
  TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
  TH1F *SumPdfWeights = new TH1F("SumPdfWeights", "SumPdfWeights", NUM_PDF_WEIGHTS, -0.5, NUM_PDF_WEIGHTS-0.5);

  //--------------
  //tree variables
  //--------------
  float weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  float triggerEffWeight;
  float triggerEffSFWeight;
  float photonEffSF;
  float photonEffSFUp, photonEffSFDown;
  float photonEffSFSys;
  float pho1EffSF, pho1EffSFUnc, pho1EffSFUp, pho1EffSFDown;
  float pho2EffSF, pho2EffSFUnc, pho2EffSFUp, pho2EffSFDown;
  float leptonEffSF;
  float leptonEffSFUp, leptonEffSFDown;
  float leptonEffSFSys;
  float lep1EffSF, lep1EffSFUnc, lep1EffSFUp, lep1EffSFDown;
  float lep2EffSF, lep2EffSFUnc, lep2EffSFUp, lep2EffSFDown;
  float ISRSystWeightUp, ISRSystWeightDown;
  int   NISRJets;
  float ptISR;
  //For btag scale factor uncertainty
  float btagCorrFactor;
  float sf_btagUp, sf_btagDown;
  float sf_bmistagUp, sf_bmistagDown;
  //For scale variation uncertainties
  float sf_facScaleUp, sf_facScaleDown;
  float sf_renScaleUp, sf_renScaleDown;
  float sf_facRenScaleUp, sf_facRenScaleDown;
  //PDF SF
  std::vector<float> sf_pdf;

  int NPU;
  int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
  float theMR, theMR_JESUp, theMR_JESDown;
  float theRsq, theRsq_JESUp, theRsq_JESDown, t1Rsq, t1Rsq_JESUp, t1Rsq_JESDown;
  float genMetRsq;
  float MET, MET_JESUp, MET_JESDown, t1MET, t1MET_JESUp, t1MET_JESDown;
  float MET_RecV2, t1MET_RecV2;
  float Myt1MET;
  float HT;


  int nSelectedPhotons;
  float mGammaGamma, pTGammaGamma, mGammaGammaSC, pTGammaGammaSC, sigmaMoverM;
  float mbbZ, mbbZ_L, mbbH, mbbH_L;
  float pTbbZ, pTbbZ_L, pTbbH, pTbbH_L;
  bool passedDiphotonTrigger;
  SusyEwkHggBox razorbox = None;

  unsigned int run, lumi, event;

  //selected lepton variables
  int lep1Type = 0;
  int lep1PassSelection = 0;
  float lep1Pt = -999;
  float lep1Eta = -999;
  float lep1Phi = -999;
  int lep2Type = 0;
  int lep2PassSelection = 0;
  float lep2Pt = -999;
  float lep2Eta = -999;
  float lep2Phi = -999;
  float dileptonMass = -999;
  float lep1MT = -999;
  float lep1GenMetMT = -999;

  //selected photon variables
  float Pho_E[2], Pho_scale[2], Pho_Pt[2], Pho_Eta[2], Pho_Phi[2], Pho_SigmaIetaIeta[2], Pho_R9[2], Pho_HoverE[2];
  float PhoSC_E[2], PhoSC_Pt[2], PhoSC_Eta[2], PhoSC_Phi[2];
  float Pho_sumChargedHadronPt[2], Pho_sumNeutralHadronEt[2], Pho_sumPhotonEt[2], Pho_sigmaEOverE[2];
  bool  Pho_passEleVeto[2], Pho_passIso[2];
  int   Pho_motherID[2];

  //jet information
  int n_Jets, n_BJets, nLooseBTaggedJets, nMediumBTaggedJets;
  int n_Jets_JESUp, n_Jets_JESDown;
  float jet_E[50], jet_Pt[50], jet_Eta[50], jet_Phi[50];
  bool jetIsCSVL[50], jetIsCSVM[50], jetIsCSVT[50];
  float bjet_E[50], bjet_Pt[50], bjet_Eta[50], bjet_Phi[50];
  bool bjetIsCSVL[50], bjetIsCSVM[50], bjetIsCSVT[50];

  //ECALGainSwitchFlag
  bool Flag_hasEcalGainSwitch = false;

  //SMS info
  int mChi = 0;
  int mLSP = 0;

  float N2Pt = 0;
  float N2Eta = 0;
  float N2Phi = 0;
  float N2Mass = 0;
  float N3Pt = 0;
  float N3Eta = 0;
  float N3Phi = 0;
  float N3Mass = 0;
  float HPt = 0;
  float HEta = 0;
  float HPhi = 0;
  float HMass = 0;
  float ZPt = 0;
  float ZEta = 0;
  float ZPhi = 0;
  float ZMass = 0;
  float N1FromN2Pt = 0;
  float N1FromN2Eta = 0;
  float N1FromN2Phi = 0;
  float N1FromN2Mass = 0;
  float N1FromN3Pt = 0;
  float N1FromN3Eta = 0;
  float N1FromN3Phi = 0;
  float N1FromN3Mass = 0;


  //------------------------
  //set branches on big tree
  //------------------------
  razorTree->Branch("genMetPt", &genMetPt, "genMetPt/F");
  razorTree->Branch("N2Pt", &N2Pt, "N2Pt/F");
  razorTree->Branch("N2Eta", &N2Eta, "N2Eta/F");
  razorTree->Branch("N2Phi", &N2Phi, "N2Phi/F");
  razorTree->Branch("N2Mass", &N2Mass, "N2Mass/F");
  razorTree->Branch("N3Pt", &N3Pt, "N3Pt/F");
  razorTree->Branch("N3Eta", &N3Eta, "N3Eta/F");
  razorTree->Branch("N3Phi", &N3Phi, "N3Phi/F");
  razorTree->Branch("N3Mass", &N3Mass, "N3Mass/F");
  razorTree->Branch("HPt", &HPt, "HPt/F");
  razorTree->Branch("HEta", &HEta, "HEta/F");
  razorTree->Branch("HPhi", &HPhi, "HPhi/F");
  razorTree->Branch("HMass", &HMass, "HMass/F");
  razorTree->Branch("ZPt", &ZPt, "ZPt/F");
  razorTree->Branch("ZEta", &ZEta, "ZEta/F");
  razorTree->Branch("ZPhi", &ZPhi, "ZPhi/F");
  razorTree->Branch("ZMass", &ZMass, "ZMass/F");
  razorTree->Branch("N1FromN2Pt", &N1FromN2Pt, "N1FromN2Pt/F");
  razorTree->Branch("N1FromN2Eta", &N1FromN2Eta, "N1FromN2Eta/F");
  razorTree->Branch("N1FromN2Phi", &N1FromN2Phi, "N1FromN2Phi/F");
  razorTree->Branch("N1FromN2Mass", &N1FromN2Mass, "N1FromN2Mass/F");
  razorTree->Branch("N1FromN3Pt", &N1FromN3Pt, "N1FromN3Pt/F");
  razorTree->Branch("N1FromN3Eta", &N1FromN3Eta, "N1FromN3Eta/F");
  razorTree->Branch("N1FromN3Phi", &N1FromN3Phi, "N1FromN3Phi/F");
  razorTree->Branch("N1FromN3Mass", &N1FromN3Mass, "N1FromN3Mass/F");

  razorTree->Branch("weight", &weight, "weight/F");
  razorTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
  razorTree->Branch("pileupWeightUp", &pileupWeightUp, "pileupWeightUp/F");
  razorTree->Branch("pileupWeightDown", &pileupWeightDown, "pileupWeightDown/F");
  razorTree->Branch("triggerEffWeight", &triggerEffWeight, "triggerEffWeight/F");
  razorTree->Branch("triggerEffSFWeight", &triggerEffSFWeight, "triggerEffSFWeight/F");
  razorTree->Branch("photonEffSF", &photonEffSF, "photonEffSF/F");
  razorTree->Branch("photonEffSFUp", &photonEffSFUp, "photonEffSFUp/F");
  razorTree->Branch("photonEffSFDown", &photonEffSFDown, "photonEffSFDown/F");
  razorTree->Branch("photonEffSFSys", &photonEffSFSys, "photonEffSFSys/F");
  razorTree->Branch("pho1EffSF", &pho1EffSF, "pho1EffSF/F");
  razorTree->Branch("pho1EffSFUnc", &pho1EffSFUnc, "pho1EffSFUnc/F");
  razorTree->Branch("pho1EffSFUp", &pho1EffSFUp, "pho1EffSFUp/F");
  razorTree->Branch("pho1EffSFDown", &pho1EffSFDown, "pho1EffSFDown/F");
  razorTree->Branch("pho2EffSF", &pho2EffSF, "pho2EffSF/F");
  razorTree->Branch("pho2EffSFUnc", &pho2EffSFUnc, "pho2EffSFUnc/F");
  razorTree->Branch("pho2EffSFUp", &pho2EffSFUp, "pho2EffSFUp/F");
  razorTree->Branch("pho2EffSFDown", &pho2EffSFDown, "pho2EffSFDown/F");
  razorTree->Branch("leptonEffSF", &leptonEffSF, "leptonEffSF/F");
  razorTree->Branch("leptonEffSFUp", &leptonEffSFUp, "leptonEffSFUp/F");
  razorTree->Branch("leptonEffSFDown", &leptonEffSFDown, "leptonEffSFDown/F");
  razorTree->Branch("leptonEffSFSys", &leptonEffSFSys, "leptonEffSFSys/F");
  razorTree->Branch("lep1EffSF", &lep1EffSF, "lep1EffSF/F");
  razorTree->Branch("lep1EffSFUnc", &lep1EffSFUnc, "lep1EffSFUnc/F");
  razorTree->Branch("lep1EffSFUp", &lep1EffSFUp, "lep1EffSFUp/F");
  razorTree->Branch("lep1EffSFDown", &lep1EffSFDown, "lep1EffSFDown/F");
  razorTree->Branch("lep2EffSF", &lep2EffSF, "lep2EffSF/F");
  razorTree->Branch("lep2EffSFUnc", &lep2EffSFUnc, "lep2EffSFUnc/F");
  razorTree->Branch("lep2EffSFUp", &lep2EffSFUp, "lep2EffSFUp/F");
  razorTree->Branch("lep2EffSFDown", &lep2EffSFDown, "lep2EffSFDown/F");
  razorTree->Branch("ISRSystWeightUp", &ISRSystWeightUp, "ISRSystWeightUp/F");
  razorTree->Branch("ISRSystWeightDown", &ISRSystWeightDown, "ISRSystWeightDown/F");
  razorTree->Branch("NISRJets", &NISRJets, "NISRJets/I");
  razorTree->Branch("ptISR", &ptISR, "ptISR/F");

  razorTree->Branch("btagCorrFactor", &btagCorrFactor, "btagCorrFactor/F");
  razorTree->Branch("sf_btagUp", &sf_btagUp, "sf_btagUp/F");
  razorTree->Branch("sf_btagDown", &sf_btagDown, "sf_btagDown/F");
  razorTree->Branch("sf_bmistagUp", &sf_bmistagUp, "sf_bmistagUp/F");
  razorTree->Branch("sf_bmistagDown", &sf_bmistagDown, "sf_bmistagDown/F");

  razorTree->Branch("sf_facScaleUp", &sf_facScaleUp, "sf_facScaleUp/F");
  razorTree->Branch("sf_facScaleDown", &sf_facScaleDown, "sf_facScaleDown/F");
  razorTree->Branch("sf_renScaleUp", &sf_renScaleUp, "sf_renScaleUp/F");
  razorTree->Branch("sf_renScaleDown", &sf_renScaleDown, "sf_renScaleDown/F");
  razorTree->Branch("sf_facRenScaleUp", &sf_facRenScaleUp, "sf_facRenScaleUp/F");
  razorTree->Branch("sf_facRenScaleDown", &sf_facRenScaleDown, "sf_facRenScaleDown/F");
  razorTree->Branch("pdfWeights", "std::vector<float>",&pdfWeights); //get PDF weights directly from RazorEvents
  razorTree->Branch("sf_pdf", "std::vector<float>",&sf_pdf); //sf PDF

  //MET filters
  razorTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
  razorTree->Branch("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, "Flag_HBHEIsoNoiseFilter/O");
  razorTree->Branch("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter, "Flag_badChargedCandidateFilter/O");
  razorTree->Branch("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, "Flag_BadChargedCandidateFilter/O");
  razorTree->Branch("Flag_badMuonFilter", &Flag_badMuonFilter, "Flag_badMuonFilter/O");
  razorTree->Branch("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, "Flag_BadPFMuonFilter/O");
  razorTree->Branch("Flag_badGlobalMuonFilter", &Flag_badGlobalMuonFilter, "Flag_badGlobalMuonFilter/O");
  razorTree->Branch("Flag_duplicateMuonFilter", &Flag_duplicateMuonFilter, "Flag_duplicateMuonFilter/O");
  razorTree->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, "Flag_CSCTightHaloFilter/O");
  razorTree->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, "Flag_hcalLaserEventFilter/O");
  razorTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/O");
  razorTree->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/O");
  razorTree->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, "Flag_trackingFailureFilter/O");
  razorTree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/O");
  razorTree->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, "Flag_ecalLaserCorrFilter/O");
  razorTree->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters, "Flag_trkPOGFilters/O");
  razorTree->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, "Flag_trkPOG_manystripclus53X/O");
  razorTree->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, "Flag_trkPOG_toomanystripclus53X/O");
  razorTree->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, "Flag_trkPOG_logErrorTooManyClusters/O");
  razorTree->Branch("Flag_METFilters", &Flag_METFilters, "Flag_METFilters/O");
  razorTree->Branch("Flag_hasEcalGainSwitch", &Flag_hasEcalGainSwitch, "Flag_hasEcalGainSwitch/O");
  razorTree->Branch("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, "Flag_ecalBadCalibFilter/O");

  razorTree->Branch("run", &run, "run/i");
  razorTree->Branch("lumi", &lumi, "lumi/i");
  razorTree->Branch("event", &event, "event/i");
  razorTree->Branch("passedDiphotonTrigger", &passedDiphotonTrigger, "passedDiphotonTrigger/O");
  razorTree->Branch("NPU", &NPU, "npu/i");
  razorTree->Branch("nPV", &nPV, "nPV/i");
  razorTree->Branch("nLooseBTaggedJets", &nLooseBTaggedJets, "nLooseBTaggedJets/I");
  razorTree->Branch("nMediumBTaggedJets", &nMediumBTaggedJets, "nMediumBTaggedJets/I");
  razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
  razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
  razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
  razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
  razorTree->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
  razorTree->Branch("MR", &theMR, "MR/F");
  razorTree->Branch("MR_JESUp", &theMR_JESUp, "MR_JESUp/F");
  razorTree->Branch("MR_JESDown", &theMR_JESDown, "MR_JESDown/F");
  razorTree->Branch("Rsq", &theRsq, "Rsq/F");
  razorTree->Branch("Rsq_JESUp", &theRsq_JESUp, "Rsq_JESUp/F");
  razorTree->Branch("Rsq_JESDown", &theRsq_JESDown, "Rsq_JESDown/F");
  razorTree->Branch("t1Rsq", &t1Rsq, "t1Rsq/F");
  razorTree->Branch("t1Rsq_JESUp", &t1Rsq_JESUp, "t1Rsq_JESUp/F");
  razorTree->Branch("t1Rsq_JESDown", &t1Rsq_JESDown, "t1Rsq_JESDown/F");
  razorTree->Branch("genMetRsq", &genMetRsq, "genMetRsq/F");
  razorTree->Branch("MET", &MET, "MET/F");
  razorTree->Branch("MET_JESUp", &MET_JESUp, "MET_JESUp/F");
  razorTree->Branch("MET_JESDown", &MET_JESDown, "MET_JESDown/F");
  razorTree->Branch("t1MET", &t1MET, "t1MET/F");
  razorTree->Branch("t1MET_JESUp", &t1MET_JESUp, "t1MET_JESUp/F");
  razorTree->Branch("t1MET_JESDown", &t1MET_JESDown, "t1MET_JESDown/F");
  razorTree->Branch("HT", &HT, "HT/F");

  //2017F MET RECIPE
  razorTree->Branch("MET_RecV2", &MET_RecV2, "MET_RecV2/F");
  razorTree->Branch("t1MET_RecV2", &t1MET_RecV2, "t1MET_RecV2/F");
  razorTree->Branch("Myt1MET", &Myt1MET, "Myt1MET/F");
  

  razorTree->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
  razorTree->Branch("mGammaGamma", &mGammaGamma, "mGammaGamma/F");
  razorTree->Branch("pTGammaGamma", &pTGammaGamma, "pTGammaGamma/F");
  razorTree->Branch("mGammaGammaSC", &mGammaGammaSC, "mGammaGammaSC/F");
  razorTree->Branch("pTGammaGammaSC", &pTGammaGammaSC, "pTGammaGammaSC/F");
  razorTree->Branch("sigmaMoverM", &sigmaMoverM, "sigmaMoverM/F");
  razorTree->Branch("box", &razorbox, "box/I");

  razorTree->Branch("lep1Type", &lep1Type, "lep1Type/I");
  razorTree->Branch("lep1PassSelection", &lep1PassSelection, "lep1PassSelection/I");
  razorTree->Branch("lep1Pt", &lep1Pt, "lep1Pt/F");
  razorTree->Branch("lep1Eta", &lep1Eta, "lep1Eta/F");
  razorTree->Branch("lep1Phi", &lep1Phi, "lep1Phi/F");
  razorTree->Branch("lep2Type", &lep2Type, "lep2Type/I");
  razorTree->Branch("lep2PassSelection", &lep2PassSelection, "lep2PassSelection/I");
  razorTree->Branch("lep2Pt", &lep2Pt, "lep2Pt/F");
  razorTree->Branch("lep2Eta", &lep2Eta, "lep2Eta/F");
  razorTree->Branch("lep2Phi", &lep2Phi, "lep2Phi/F");
  razorTree->Branch("dileptonMass", &dileptonMass, "dileptonMass/F");
  razorTree->Branch("lep1MT", &lep1MT, "lep1MT/F");
  razorTree->Branch("lep1GenMetMT", &lep1GenMetMT, "lep1GenMetMT/F");

  razorTree->Branch("pho1E", &Pho_E[0], "pho1E/F");
  //razorTree->Branch("pho1scale", &Pho_scale[0], "pho1scale/F");//only used for debugging
  razorTree->Branch("pho1Pt", &Pho_Pt[0], "pho1Pt/F");
  razorTree->Branch("pho1Eta", &Pho_Eta[0], "pho1Eta/F");
  razorTree->Branch("pho1Phi", &Pho_Phi[0], "pho1Phi/F");
  razorTree->Branch("pho1SC_E", &PhoSC_E[0], "pho1SC_E/F");
  razorTree->Branch("pho1SC_Pt", &PhoSC_Pt[0], "pho1SC_Pt/F");
  razorTree->Branch("pho1SC_Eta", &PhoSC_Eta[0], "pho1SC_Eta/F");
  razorTree->Branch("pho1SC_Phi", &PhoSC_Phi[0], "pho1SC_Phi/F");
  razorTree->Branch("pho1SigmaIetaIeta", &Pho_SigmaIetaIeta[0], "pho1SigmaIetaIeta/F");
  razorTree->Branch("pho1R9", &Pho_R9[0], "pho1R9/F");
  razorTree->Branch("pho1HoverE", &Pho_HoverE[0], "pho1HoverE/F");
  razorTree->Branch("pho1sumChargedHadronPt", &Pho_sumChargedHadronPt[0], "pho1sumChargedHadronPt/F");
  razorTree->Branch("pho1sumNeutralHadronEt", &Pho_sumNeutralHadronEt[0], "pho1sumNeutralHadronEt/F");
  razorTree->Branch("pho1sumPhotonEt", &Pho_sumPhotonEt[0], "pho1sumPhotonEt/F");
  razorTree->Branch("pho1sigmaEOverE", &Pho_sigmaEOverE[0], "pho1sigmaEOverE/F");
  razorTree->Branch("pho1passEleVeto", &Pho_passEleVeto[0], "pho1passEleVeto/O");
  razorTree->Branch("pho1passIso", &Pho_passIso[0], "pho1passIso/O");
  razorTree->Branch("pho1MotherID", &Pho_motherID[0], "pho1MotherID/I");

  razorTree->Branch("pho2E", &Pho_E[1], "pho2E/F");
  //razorTree->Branch("pho2scale", &Pho_scale[1], "pho2scale/F");//only used for debugging
  razorTree->Branch("pho2Pt", &Pho_Pt[1], "pho2Pt/F");
  razorTree->Branch("pho2Eta", &Pho_Eta[1], "pho2Eta/F");
  razorTree->Branch("pho2Phi", &Pho_Phi[1], "pho2Phi/F");
  razorTree->Branch("pho2SC_E", &PhoSC_E[1], "pho2SC_E/F");
  razorTree->Branch("pho2SC_Pt", &PhoSC_Pt[1], "pho2SC_Pt/F");
  razorTree->Branch("pho2SC_Eta", &PhoSC_Eta[1], "pho2SC_Eta/F");
  razorTree->Branch("pho2SC_Phi", &PhoSC_Phi[1], "pho2SC_Phi/F");
  razorTree->Branch("pho2SigmaIetaIeta", &Pho_SigmaIetaIeta[1], "pho2SigmaIetaIeta/F");
  razorTree->Branch("pho2R9", &Pho_R9[1], "pho2R9/F");
  razorTree->Branch("pho2HoverE", &Pho_HoverE[1], "pho2HoverE/F");
  razorTree->Branch("pho2sumChargedHadronPt", &Pho_sumChargedHadronPt[1], "pho2sumChargedHadronPt/F");
  razorTree->Branch("pho2sumNeutralHadronEt", &Pho_sumNeutralHadronEt[1], "pho2sumNeutralHadronEt/F");
  razorTree->Branch("pho2sumPhotonEt", &Pho_sumPhotonEt[1], "pho2sumPhotonEt/F");
  razorTree->Branch("pho2sigmaEOverE", &Pho_sigmaEOverE[1], "pho2sigmaEOverE/F");
  razorTree->Branch("pho2passEleVeto", &Pho_passEleVeto[1], "pho2passEleVeto/O");
  razorTree->Branch("pho2passIso", &Pho_passIso[1], "pho2passIso/O");
  razorTree->Branch("pho2MotherID", &Pho_motherID[1], "pho2MotherID/I");

  razorTree->Branch("mbbZ", &mbbZ, "mbbZ/F");
  razorTree->Branch("mbbH", &mbbH, "mbbH/F");
  razorTree->Branch("mbbZ_L", &mbbZ_L, "mbbZ_L/F");
  razorTree->Branch("mbbH_L", &mbbH_L, "mbbH_L/F");
  razorTree->Branch("pTbbZ", &pTbbZ, "pTbbZ/F");
  razorTree->Branch("pTbbH", &pTbbH, "pTbbH/F");
  razorTree->Branch("pTbbZ_L", &pTbbZ_L, "pTbbZ_L/F");
  razorTree->Branch("pTbbH_L", &pTbbH_L, "pTbbH_L/F");

  razorTree->Branch("n_Jets", &n_Jets, "n_Jets/I");
  razorTree->Branch("jet_E", jet_E, "jet_E[n_Jets]/F");
  razorTree->Branch("jet_Pt", jet_Pt, "jet_Pt[n_Jets]/F");
  razorTree->Branch("jet_Eta", jet_Eta, "jet_Eta[n_Jets]/F");
  razorTree->Branch("jet_Phi", jet_Phi, "jet_Phi[n_Jets]/F");
  razorTree->Branch("jetIsCSVL", jetIsCSVL, "jetIsCSVL[n_Jets]/O");
  razorTree->Branch("jetIsCSVM", jetIsCSVM, "jetIsCSVM[n_Jets]/O");
  razorTree->Branch("jetIsCSVT", jetIsCSVT, "jetIsCSVT[n_Jets]/O");
  razorTree->Branch("n_BJets", &n_BJets, "n_BJets/I");
  //razorTree->Branch("bjet_E", bjet_E, "bjet_E[n_BJets]/F");
  //razorTree->Branch("bjet_Pt", bjet_Pt, "bjet_Pt[n_BJets]/F");
  //razorTree->Branch("bjet_Eta", bjet_Eta, "bjet_Eta[n_BJets]/F");
  //razorTree->Branch("bjet_Phi", bjet_Phi, "bjet_Phi[n_BJets]/F");
  //razorTree->Branch("bjetIsCSVL", bjetIsCSVL, "bjetIsCSVL[n_BJets]/O");
  //razorTree->Branch("bjetIsCSVM", bjetIsCSVM, "bjetIsCSVM[n_BJets]/O");
  //razorTree->Branch("bjetIsCSVT", bjetIsCSVT, "bjetIsCSVT[n_BJets]/O");
  razorTree->Branch("n_Jets_JESUp", &n_Jets_JESUp, "n_Jets_JESUp/I");
  razorTree->Branch("n_Jets_JESDown", &n_Jets_JESDown, "n_Jets_JESDown/I");
  razorTree->Branch("HLTDecision", HLTDecision, "HLTDecision[300]/O");

  //GenParticles
  razorTree->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
  razorTree->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
  razorTree->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
  razorTree->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
  razorTree->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
  razorTree->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
  razorTree->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
  razorTree->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
  razorTree->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");

  razorTree->Branch("mChi", &mChi, "mChi/I");
  razorTree->Branch("mLSP", &mLSP, "mLSP/I");

  //begin loop
  if ( fChain == 0 ) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  std::cout << "[INFO]: Total Entries = " << fChain->GetEntries() << "\n";
  for ( Long64_t jentry=0; jentry < nentries; jentry++ )
    {
      if(_metdebug||_ecaldebug) std::cout << "[INFO]: Processing entry " << jentry << std::endl;
      //begin event
      if( _info && (jentry % 10000 == 0) ) std::cout << "[INFO]: Processing entry " << jentry << std::endl;
      Long64_t ientry = LoadTree( jentry );
      if ( ientry < 0 ) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;



      //***********************************************
      //Filter out the Pathological Fastsim Events
      //***********************************************
      if (isFastsimSMS) {
	bool isPathologicalFastsimEvent = false;

	for (int i = 0; i < nJets; i++){
	  double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
						 fixedGridRhoFastjetAll, jetJetArea[i], runNum, JetCorrectorIOV,JetCorrector);
	  double jetCorrPt = jetPt[i]*JEC;
	  if (jetCorrPt < 20) continue;
	  if (fabs(jetEta[i]) > 2.5) continue;

	  //Match to Gen Jet
	  bool isMatch = false;
	  for(int j = 0; j < nGenJets; j++){
	    double tmpDR = deltaR( genJetEta[j],genJetPhi[j], jetEta[i],jetPhi[i]);
	    if ( tmpDR < 0.4
		 ) {
	      isMatch = true;
	    }
	  }

	  // these are the pathological fastsim jets
	  if (!isMatch && jetChargedHadronEnergyFraction[i] < 0.1 ) {
	    isPathologicalFastsimEvent = true;
	  }
	}
	//reject event if it's pathological
	if (isPathologicalFastsimEvent) continue;
      }


      //fill normalization histogram
      NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
      weight = genWeight;
      SumWeights->Fill(1.0, weight);

      //reset tree variables
      ISRSystWeightUp    = 1.0;
      ISRSystWeightDown  = 1.0;
      NISRJets           = 0;
      ptISR              = -1;
      pileupWeight       = 1.0;
      pileupWeightUp     = 1.0;
      pileupWeightDown   = 1.0;
      triggerEffWeight   = 1.0;
      triggerEffSFWeight = 1.0;
      leptonEffSF        = 1.0;
      leptonEffSFUp      = 1.0;
      leptonEffSFDown    = 1.0;
      leptonEffSFSys     = 0.0;
      lep1EffSF          = 1.0;
      lep1EffSFUnc       = 0.0;
      lep1EffSFUp        = 1.0;
      lep1EffSFDown      = 1.0;
      lep2EffSF          = 1.0;
      lep2EffSFUnc       = 0.0;
      lep2EffSFUp        = 1.0;
      lep2EffSFDown      = 1.0;
      photonEffSF        = 1.0;
      photonEffSFUp      = 1.0;
      photonEffSFDown    = 1.0;
      photonEffSFSys     = 0.0;
      pho1EffSF          = 1.0;
      pho1EffSFUnc       = 0.0;
      pho1EffSFUp        = 1.0;
      pho1EffSFDown      = 1.0;
      pho2EffSF          = 1.0;
      pho2EffSFUnc       = 0.0;
      pho2EffSFUp        = 1.0;
      pho2EffSFDown      = 1.0;

      btagCorrFactor    = 1.0;
      sf_btagUp         = 1.0;
      sf_btagDown       = 1.0;
      sf_bmistagUp      = 1.0;
      sf_bmistagDown    = 1.0;

      sf_facScaleUp = 1.0;
      sf_facScaleDown = 1.0;
      sf_renScaleUp = 1.0;
      sf_renScaleDown = 1.0;
      sf_facRenScaleUp = 1.0;
      sf_facRenScaleDown = 1.0;

      n_Jets = 0;
      n_BJets = 0;
      n_Jets_JESUp = 0;
      n_Jets_JESDown = 0;
      nLooseBTaggedJets = 0;
      nMediumBTaggedJets = 0;
      nLooseMuons = 0;
      nTightMuons = 0;
      nLooseElectrons = 0;
      nTightElectrons = 0;
      nTightTaus = 0;
      theMR = -666;
      theMR_JESUp   = -666;
      theMR_JESDown = -666;
      theRsq = -666;
      theRsq_JESUp   = -666;
      theRsq_JESDown = -666;
      t1Rsq  = -666;
      t1Rsq_JESUp   = -666;
      t1Rsq_JESDown = -666;
      genMetRsq     = -666;

      nSelectedPhotons = 0;
      mGammaGamma    = -1;
      pTGammaGamma   = -1;
      mGammaGammaSC  = -1;
      pTGammaGammaSC = -1;
      mbbZ   = 0;
      mbbH   = 0;
      mbbZ_L = 0;
      mbbH_L = 0;
      pTbbZ   = 0;
      pTbbH   = 0;
      pTbbZ_L = 0;
      pTbbH_L = 0;
      run = runNum;
      lumi = lumiNum;
      event = eventNum;
      passedDiphotonTrigger = false;
      Flag_hasEcalGainSwitch = false;

      //lepton variables
      lep1Type = 0;
      lep1PassSelection = 0;
      lep1Pt = -999;
      lep1Eta = -999;
      lep1Phi = -999;
      lep2Type = 0;
      lep2PassSelection = 0;
      lep2Pt = -999;
      lep2Eta = -999;
      lep2Phi = -999;
      dileptonMass = -999;
      lep1MT = -999;
      lep1GenMetMT = -999;

      //selected photons  and two leptons variables
      for ( int i = 0; i < 2; i++ )
	{
	  Pho_E[i]                  = -99.;
	  Pho_Pt[i]                 = -99.;
	  Pho_Eta[i]                = -99.;
	  Pho_Phi[i]                = -99.;
	  PhoSC_E[i]                = -99.;
	  PhoSC_Pt[i]               = -99.;
	  PhoSC_Eta[i]              = -99.;
	  PhoSC_Phi[i]              = -99.;
	  Pho_SigmaIetaIeta[i]      = -99.;
	  Pho_R9[i]                 = -99.;
	  Pho_HoverE[i]             = -99.;
	  Pho_sumChargedHadronPt[i] = -99.;
	  Pho_sumNeutralHadronEt[i] = -99.;
	  Pho_sumPhotonEt[i]        = -99.;
	  Pho_sigmaEOverE[i]        = -99.;
	  Pho_passEleVeto[i]        = false;
	  Pho_passIso[i]            = false;
	  Pho_motherID[i]           = 0;
	}

      //jets
      for ( int i = 0; i < 50; i++ )
	{
	  jet_E[i]      = -99.;
	  jet_Pt[i]     = -99.;
	  jet_Eta[i]    = -99.;
	  jet_Phi[i]    = -99.;
	  bjet_E[i]     = -99.;
	  bjet_Pt[i]    = -99.;
	  bjet_Eta[i]   = -99.;
	  bjet_Phi[i]   = -99.;
	  bjetIsCSVL[i] = 0;
	  bjetIsCSVM[i] = 0;
	  bjetIsCSVT[i] = 0;
	}

      mChi = 0;
      mLSP = 0;


      //--------------------------------------------------------------
      //Extract SUSY model parameters from lheComment variable
      //--------------------------------------------------------------
       bool parsedLHE = false;
      if(isFastsimSMS && lheComments){
	//cout << lheComments << " " << *lheComments << "\n";

	//Save some information on signal particles
	bool foundV1 = false;
	bool foundV2 = false;
	TLorentzVector v1;
	TLorentzVector v2;
	for(int g = 0; g < nGenParticle; g++){
	  //cout << gParticleId[g] << " " << gParticleStatus[g] << " " << gParticlePt[g] << " " << gParticleEta[g] << " | " << gParticleMotherId[g] << "\n";
	  if (gParticleStatus[g]  == 62) {
            //62: outgoing subprocess particle with primordial kT included
	    if (!foundV1) {
	      v1.SetPtEtaPhiE( gParticlePt[g], gParticleEta[g], gParticlePhi[g], gParticleE[g]);
	      foundV1 = true;
	    } else if(!foundV2) {
	      v2.SetPtEtaPhiE( gParticlePt[g], gParticleEta[g], gParticlePhi[g], gParticleE[g]);
	      foundV2 = true;
	    } else {
	      cout << "Warning: found more than two status=62 particles\n";
	    }
	  }
	}
	if ( foundV1 && foundV2) {
	  ptISR = (v1+v2).Pt();
	}


	if (!is2DMassScan) {


	  //Save some information on signal particles
	  for(int g = 0; g < nGenParticle; g++){

	    // original N2 produced
	    if (gParticleId[g] == 1000023 && gParticleStatus[g]  == 62) {
	      N2Pt = gParticlePt[g];
	      N2Eta = gParticleEta[g];
	      N2Phi = gParticlePhi[g];
	      N2Mass = sqrt( gParticleE[g]*gParticleE[g] - pow(gParticlePt[g]*cosh(gParticleEta[g]),2));


	    }
	    // original N3 produced
	    if (gParticleId[g] == 1000025 && gParticleStatus[g]  == 62) {
	      N3Pt = gParticlePt[g];
	      N3Eta = gParticleEta[g];
	      N3Phi = gParticlePhi[g];
	      N3Mass = sqrt( gParticleE[g]*gParticleE[g] - pow(gParticlePt[g]*cosh(gParticleEta[g]),2));
	    }


	    // H from N2 -> H N1 decay
	    if (gParticleId[g] == 25 && gParticleMotherId[g]  == 1000023) {
	      HPt = gParticlePt[g];
	      HEta = gParticleEta[g];
	      HPhi = gParticlePhi[g];
	      HMass = sqrt( gParticleE[g]*gParticleE[g] - pow(gParticlePt[g]*cosh(gParticleEta[g]),2));
	    }
	    //N1 from N2 -> H N1 decay
	    if (gParticleId[g] == 1000022 && gParticleMotherId[g]  == 1000023) {
	      N1FromN2Pt = gParticlePt[g];
	      N1FromN2Eta = gParticleEta[g];
	      N1FromN2Phi = gParticlePhi[g];
	      N1FromN2Mass = sqrt( gParticleE[g]*gParticleE[g] - pow(gParticlePt[g]*cosh(gParticleEta[g]),2));
	    }
	    // Z from N3 -> Z N1 decay
	    if (gParticleId[g] == 23 && gParticleMotherId[g]  == 1000025) {
	      ZPt = gParticlePt[g];
	      ZEta = gParticleEta[g];
	      ZPhi = gParticlePhi[g];
	      ZMass = sqrt( gParticleE[g]*gParticleE[g] - pow(gParticlePt[g]*cosh(gParticleEta[g]),2));
	      //cout << "foudn Z : " << ZPt << " " << ZMass << "\n";
	    }
	    //N1 from N3 -> Z N1 decay
	    if (gParticleId[g] == 1000022 && gParticleMotherId[g]  == 1000025) {
	      N1FromN3Pt = gParticlePt[g];
	      N1FromN3Eta = gParticleEta[g];
	      N1FromN3Phi = gParticlePhi[g];
	      N1FromN3Mass = sqrt( gParticleE[g]*gParticleE[g] - pow(gParticlePt[g]*cosh(gParticleEta[g]),2));
	    }

	  }





	  //parse lhe comment string to get Chargino/Neutralino2 masses, LHE:  Les Houches Event file format
	  stringstream parser(*lheComments);
	  string item;
	  getline(parser, item, '_'); //prefix
	  if(getline(parser, item, '_')){ //gluino mass
	    mChi = atoi(item.c_str());
	    if(mChi == 0) { //fix for the case where the model name contains an underscore
	      if(getline(parser, item, '_')){
		mChi = atoi(item.c_str());
	      }

	      parsedLHE = true;


	      if (fabs( mChi - N2Mass ) > 1) {
		cout << "Weird: " << mChi << " | " << N2Mass << " " << N3Mass << " " << N1FromN2Mass << " " << N1FromN3Mass << " " << HMass << " " << ZMass << " : " << *lheComments << "\n";
		cout << "Throwing the event out. \n";
		continue;
	      }

	      if (smsFiles.count(mChi) == 0){ //create file and tree
		//format file name
		string thisFileName = outFileName;
		thisFileName.erase(thisFileName.end()-5, thisFileName.end());
		thisFileName += "_" + to_string(mChi) + ".root";

		smsFiles[mChi] = new TFile(thisFileName.c_str(), "recreate");
		smsTrees[mChi] = razorTree->CloneTree(0);
		smsNEvents[mChi] = new TH1F(Form("NEvents%d", mChi), "NEvents", 1,0.5,1.5);
		smsSumWeights[mChi] = new TH1F(Form("SumWeights%d", mChi), "SumWeights", 1,0.5,1.5);
		smsSumScaleWeights[mChi] = new TH1F(Form("SumScaleWeights%d", mChi), "SumScaleWeights", 6,-0.5,5.5);
		smsSumPdfWeights[mChi] = new TH1F(Form("SumPdfWeights%d", mChi), "SumPdfWeights", NUM_PDF_WEIGHTS,-0.5,NUM_PDF_WEIGHTS-0.5);
		smsNISRJets[mChi] = new TH1F(Form("NISRJets%d", mChi), "NISRJets", 7,-0.5,6.5);
		smsPtISR[mChi] = new TH1F(Form("PtISR%d", mChi), "PtISR", 8, PtISRBins);
		smsNPV[mChi] = new TH1F(Form("NPV%d", mChi), "NPV", 2,-0.5,1.5);
		cout << "Created new output file " << thisFileName << endl;
	      }
	      //Fill NEvents hist
	      smsNEvents[mChi]->Fill(1.0, genWeight);
	      smsSumWeights[mChi]->Fill(1.0, weight);
	      smsNISRJets[mChi]->Fill(min(NISRJets,6), genWeight);
	      smsPtISR[mChi]->Fill(fmin( ptISR , 6999.0), genWeight);
	      smsNPV[mChi]->Fill( (nPV >= 20)?1:0 , genWeight);
	      smsSumScaleWeights[mChi]->Fill(0.0, sf_facScaleUp);
	      smsSumScaleWeights[mChi]->Fill(1.0, sf_facScaleDown);
	      smsSumScaleWeights[mChi]->Fill(2.0, sf_renScaleUp);
	      smsSumScaleWeights[mChi]->Fill(3.0, sf_renScaleDown);
	      smsSumScaleWeights[mChi]->Fill(4.0, sf_facRenScaleUp);
	      smsSumScaleWeights[mChi]->Fill(5.0, sf_facRenScaleDown);

	      for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) {
		smsSumPdfWeights[mChi]->Fill(double(iwgt),(*pdfWeights)[iwgt]);
	      }

	    }
	  }
	} else {

                cout << "======================================================================\n";
                cout << " HPt = " << HPt << " \n";
                cout << "======================================================================\n";

	  //Save some information on signal particles
	  for(int g = 0; g < nGenParticle; g++){

	    // original N2 produced
	    if (gParticleId[g] == 1000023 && gParticleStatus[g]  == 62) {
	      N2Pt = gParticlePt[g];
	      N2Eta = gParticleEta[g];
	      N2Phi = gParticlePhi[g];
	      N2Mass = sqrt( gParticleE[g]*gParticleE[g] - pow(gParticlePt[g]*cosh(gParticleEta[g]),2));


	    }
	    // original N3 produced
	    if (gParticleId[g] == 1000025 && gParticleStatus[g]  == 62) {
	      N3Pt = gParticlePt[g];
	      N3Eta = gParticleEta[g];
	      N3Phi = gParticlePhi[g];
	      N3Mass = sqrt( gParticleE[g]*gParticleE[g] - pow(gParticlePt[g]*cosh(gParticleEta[g]),2));
	    }


	    // H from N2 -> H N1 decay
	    if (gParticleId[g] == 25 && gParticleMotherId[g]  == 1000023) {
	      HPt = gParticlePt[g];
	      HEta = gParticleEta[g];
	      HPhi = gParticlePhi[g];
	      HMass = sqrt( gParticleE[g]*gParticleE[g] - pow(gParticlePt[g]*cosh(gParticleEta[g]),2));
	    }
	    //N1 from N2 -> H N1 decay
	    if (gParticleId[g] == 1000022 && gParticleMotherId[g]  == 1000023) {
	      N1FromN2Pt = gParticlePt[g];
	      N1FromN2Eta = gParticleEta[g];
	      N1FromN2Phi = gParticlePhi[g];
	      N1FromN2Mass = sqrt( gParticleE[g]*gParticleE[g] - pow(gParticlePt[g]*cosh(gParticleEta[g]),2));
	    }
	    // Z from N3 -> Z N1 decay
	    if (gParticleId[g] == 23 && gParticleMotherId[g]  == 1000025) {
	      ZPt = gParticlePt[g];
	      ZEta = gParticleEta[g];
	      ZPhi = gParticlePhi[g];
	      ZMass = sqrt( gParticleE[g]*gParticleE[g] - pow(gParticlePt[g]*cosh(gParticleEta[g]),2));
	      //cout << "foudn Z : " << ZPt << " " << ZMass << "\n";
	    }
	    //N1 from N3 -> Z N1 decay
	    if (gParticleId[g] == 1000022 && gParticleMotherId[g]  == 1000025) {
	      N1FromN3Pt = gParticlePt[g];
	      N1FromN3Eta = gParticleEta[g];
	      N1FromN3Phi = gParticlePhi[g];
	      N1FromN3Mass = sqrt( gParticleE[g]*gParticleE[g] - pow(gParticlePt[g]*cosh(gParticleEta[g]),2));
	    }

	  }


                cout << "======================================================================\n";
                cout << " HPt = " << HPt << " \n";
                cout << "======================================================================\n";

	  //parse lhe comment string to get gluino and LSP masses
	  stringstream parser(*lheComments);
	  string item;
	  getline(parser, item, '_'); //prefix
	  if(getline(parser, item, '_')){ //gluino mass
	    mChi = atoi(item.c_str());
	    if(mChi == 0) { //fix for the case where the model name contains an underscore
	      if(getline(parser, item, '_')){
		mChi = atoi(item.c_str());
		if(mChi == 0) { //fix for the case where the model name contains an underscore
		  if(getline(parser, item, '_')){
		    mChi = atoi(item.c_str());
		  }
		}
	      }
	    }
	    if(getline(parser, item, '_')){ //LSP mass
	      mLSP = atoi(item.c_str());
	      pair<int,int> smsPair = make_pair(mChi, mLSP);

	      parsedLHE = true;
	      if (smsFiles2D.count(smsPair) == 0){ //create file and tree
		//format file name
		string thisFileName = outFileName;
		thisFileName.erase(thisFileName.end()-5, thisFileName.end());
		thisFileName += "_" + to_string(mChi) + "_" + to_string(mLSP) + ".root";

		smsFiles2D[smsPair] = new TFile(thisFileName.c_str(), "recreate");
		smsTrees2D[smsPair] = razorTree->CloneTree(0);
		smsNEvents2D[smsPair] = new TH1F(Form("NEvents%d%d", mChi, mLSP), "NEvents", 1,0.5,1.5);
		smsSumWeights2D[smsPair] = new TH1F(Form("SumWeights%d%d", mChi, mLSP), "SumWeights", 1,0.5,1.5);
		smsSumScaleWeights2D[smsPair] = new TH1F(Form("SumScaleWeights%d%d", mChi, mLSP), "SumScaleWeights", 6,-0.5,5.5);
		smsSumPdfWeights2D[smsPair] = new TH1F(Form("SumPdfWeights%d%d", mChi, mLSP), "SumPdfWeights", NUM_PDF_WEIGHTS,-0.5,NUM_PDF_WEIGHTS-0.5);
		smsNISRJets2D[smsPair] = new TH1F(Form("NISRJets%d%d", mChi, mLSP), "NISRJets", 7,-0.5,6.5);
		smsPtISR2D[smsPair] = new TH1F(Form("PtISR%d%d", mChi, mLSP), "PtISR", 8,PtISRBins);
		smsNPV2D[smsPair] = new TH1F(Form("NPV%d%d", mChi, mLSP), "NPV", 2,-0.5,1.5);

		cout << "Created new output file " << thisFileName << endl;
	      }
	      //Fill NEvents hist
	      smsNEvents2D[smsPair]->Fill(1.0, genWeight);
	      smsSumWeights2D[smsPair]->Fill(1.0, weight);
	      smsNISRJets2D[smsPair]->Fill(min(NISRJets,6), genWeight);
	      smsPtISR2D[smsPair]->Fill(fmin( ptISR , 6999.0), genWeight);
	      smsNPV2D[smsPair]->Fill( (nPV >= 20)?1:0 , genWeight);

	      smsSumScaleWeights2D[smsPair]->Fill(0.0, sf_facScaleUp);
	      smsSumScaleWeights2D[smsPair]->Fill(1.0, sf_facScaleDown);
	      smsSumScaleWeights2D[smsPair]->Fill(2.0, sf_renScaleUp);
	      smsSumScaleWeights2D[smsPair]->Fill(3.0, sf_renScaleDown);
	      smsSumScaleWeights2D[smsPair]->Fill(4.0, sf_facRenScaleUp);
	      smsSumScaleWeights2D[smsPair]->Fill(5.0, sf_facRenScaleDown);

	      for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) {
		smsSumPdfWeights2D[smsPair]->Fill(double(iwgt),(*pdfWeights)[iwgt]);
	      }
	    }
	  }
	}
      } // end if fastsim


      //------------------
      //Pileup reweighting
      //------------------
      pileupWeight = 1.0;
      if( !isData ) {
	//Get number of PU interactions
	for (int i = 0; i < nBunchXing; i++) {
	  if (BunchXing[i] == 0) {
	    NPU = nPUmean[i];
	  }
	}
	puhisto->Fill(NPU);
	pileupWeight = helper->getPileupWeight(NPU);
	pileupWeightUp = helper->getPileupWeightUp(NPU) / pileupWeight;
	pileupWeightDown = helper->getPileupWeightDown(NPU) / pileupWeight;
      }

      /////////////////////////////////
      //Scale and PDF variations
      /////////////////////////////////
      if( !isData ) {
	if ( (*scaleWeights).size() >= 9 )
	  {
	    sf_facScaleUp      = (*scaleWeights)[1]/genWeight;
	    sf_facScaleDown    = (*scaleWeights)[2]/genWeight;
	    sf_renScaleUp      = (*scaleWeights)[3]/genWeight;
	    sf_renScaleDown    = (*scaleWeights)[6]/genWeight;
	    sf_facRenScaleUp   = (*scaleWeights)[4]/genWeight;
	    sf_facRenScaleDown = (*scaleWeights)[8]/genWeight;


	    SumScaleWeights->Fill(0.0, (*scaleWeights)[1]);
	    SumScaleWeights->Fill(1.0, (*scaleWeights)[2]);
	    SumScaleWeights->Fill(2.0, (*scaleWeights)[3]);
	    SumScaleWeights->Fill(3.0, (*scaleWeights)[6]);
	    SumScaleWeights->Fill(4.0, (*scaleWeights)[4]);
	    SumScaleWeights->Fill(5.0, (*scaleWeights)[8]);
	  }

	sf_pdf.erase( sf_pdf.begin(), sf_pdf.end() );
	for ( unsigned int iwgt = 0; iwgt < pdfWeights->size(); ++iwgt )
	  {
	    sf_pdf.push_back( pdfWeights->at(iwgt)/genWeight );
	    SumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
	  }
      }

      if ( _debug ) std::cout << "============" << std::endl;
      if ( _debug ) std::cout << "run == " << run << " && evt == " << event << std::endl;


      //*************************************************************************
      //count ISR Jets
      //*************************************************************************
      if (!isData) {
	for(int i = 0; i < nJets; i++) {

	  //Jet Corrections
	  double JEC = JetEnergyCorrectionFactor( jetPt[i], jetEta[i], jetPhi[i], jetE[i],
						  fixedGridRhoAll, jetJetArea[i], runNum,
						  JetCorrectorIOV, JetCorrector );
	  TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );

	  //these are the cuts Ana/Manuel told me to use
	  if ( thisJet.Pt() > 30 && fabs( thisJet.Eta()) < 2.4 && jetPassIDLoose[i]) {

	    //try to match to gen partons
	    //Follow prescription here: https://github.com/manuelfs/babymaker/blob/0136340602ee28caab14e3f6b064d1db81544a0a/bmaker/plugins/bmaker_full.cc#L1268-L1295
	    bool match = false;
	    for(int g = 0; g < nGenParticle; g++){

	      double dR = deltaR( gParticleEta[g], gParticlePhi[g] , thisJet.Eta() , thisJet.Phi());
	      if (dR > 0.3) continue;

	      //check match against leptons
	      if (abs(gParticleId[g]) == 11 || abs(gParticleId[g]) == 13 || abs(gParticleId[g]) == 15 ) {
		match = true;
	      }

	      //check match against prompt photons
	      if (abs(gParticleId[g]) == 22 &&
		  ( (gParticleStatus[g] == 1 && gParticleMotherId[g] != 22) || gParticleStatus[g] == 22 || gParticleStatus[g] == 23) &&
		  (abs(gParticleMotherId[g]) == 25 || abs(gParticleMotherId[g]) == 21 || abs(gParticleMotherId[g]) == 2212 ||
		   (abs(gParticleMotherId[g]) >= 1 && abs(gParticleMotherId[g]) <= 6) )
		  ) {
		match = true;
	      }

	      //match to quarks
	      if (gParticleStatus[g] == 23 && abs(gParticleId[g]) <= 5 &&
		  ( abs(gParticleMotherId[g]) == 6 ||  abs(gParticleMotherId[g]) == 23 ||  abs(gParticleMotherId[g]) == 24
		    ||  abs(gParticleMotherId[g]) == 25 ||  abs(gParticleMotherId[g]) > 1e6)) {
		match = true;
	      }
	    }
	    if (!match) NISRJets++;
	  }
	}
      }
      //Fill N ISR Jet
      histNISRJets->SetBinContent( min(NISRJets,6)+1, histNISRJets->GetBinContent(min(NISRJets,6)+1) + genWeight);
      histPtISR->Fill( fmin( ptISR , 6999.0), genWeight);
      histNPV->Fill( (nPV >= 20)?1:0 , genWeight);

      //************************************************************************
      //For Debugging
      //************************************************************************
      // if (NISRJets >=1) {
      // 	cout << "\n\n" << "NISRJets = " << NISRJets << "\n";

      // 	for(int i = 0; i < nJets; i++) {

      // 	  //Jet Corrections
      // 	  double JEC = JetEnergyCorrectionFactor( jetPt[i], jetEta[i], jetPhi[i], jetE[i],
      // 						  fixedGridRhoAll, jetJetArea[i], runNum,
      // 						  JetCorrectorIOV, JetCorrector );
      // 	  TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );

      // 	  //these are the cuts Ana/Manuel told me to use
      // 	  if ( thisJet.Pt() > 30 && fabs( thisJet.Eta()) < 2.4 && jetPassIDLoose[i]) {

      // 	    cout << "Jet : " << thisJet.Pt() << " " << thisJet.Eta() << " " << thisJet.Phi() << "\n" ;

      // 	    bool match = false;
      // 	    for(int g = 0; g < nGenParticle; g++){

      // 	      double dR = deltaR( gParticleEta[g], gParticlePhi[g] , thisJet.Eta() , thisJet.Phi());
      // 	      if (dR > 0.3) continue;

      // 	      cout << "Nearby GenParticle " << g << " : " << gParticleId[g] << " " << gParticleStatus[g] << " | " << gParticlePt[g] << " " << gParticleEta[g] << " " << gParticlePhi[g] << " : " << gParticleMotherIndex[g] << " " << gParticleMotherId[g] << "\n";

      // 	      if (abs(gParticleId[g]) == 11 || abs(gParticleId[g]) == 13 || abs(gParticleId[g]) == 15 ) {
      // 		match = true;
      // 		cout << "match lepton\n";
      // 	      }

      // 	      if (abs(gParticleId[g]) == 22 &&
      // 		  ( (gParticleStatus[g] == 1 && gParticleMotherId[g] != 22) || gParticleStatus[g] == 22 || gParticleStatus[g] == 23) &&
      // 		  (abs(gParticleMotherId[g]) == 25 || abs(gParticleMotherId[g]) == 21 || abs(gParticleMotherId[g]) == 2212 ||
      // 		   (abs(gParticleMotherId[g]) >= 1 && abs(gParticleMotherId[g]) <= 6) )
      // 		  ) {
      // 		match = true;
      // 		cout << "match prompt photon\n";
      // 	      }

      // 	      if (gParticleStatus[g] == 23 && abs(gParticleId[g]) <= 5 &&
      // 		  ( abs(gParticleMotherId[g]) == 6 ||  abs(gParticleMotherId[g]) == 23 ||  abs(gParticleMotherId[g]) == 24
      // 		    ||  abs(gParticleMotherId[g]) == 25 ||  abs(gParticleMotherId[g]) > 1e6)) {
      // 		match = true;
      // 		cout << "match parton\n";
      // 	      }
      // 	    }
      // 	  }
      // 	}

      // 	for(int g = 0; g < nGenParticle; g++){
      // 	  cout << "GenParticle " << g << " : " << gParticleId[g] << " " << gParticleStatus[g] << " | " << gParticlePt[g] << " " << gParticleEta[g] << " " << gParticlePhi[g] << " : " << gParticleMotherIndex[g] << " " << gParticleMotherId[g] << "\n";
      // 	}

      // }

      //*************************************************************************



      //*************************************************************************
      //Start Object Selection
      //*************************************************************************
      razorbox = None;

      //-------------
      //tau selection
      //-------------
      //Don't mess with Taus
      /*for( int i = 0; i < nTaus; i++ )
	{
	  if( !isTightTau(i) ) continue;
	  nTightTaus++;
	}
      */
	
	bool ecal_prefiring_affected = false;
      //------------------
      //good photon selection
      //------------------
      vector<TLorentzVector> GoodPhotons;
      vector<double> GoodPhotonSigmaE; // energy uncertainties of selected photons
      vector<bool> GoodPhotonPassesIso; //store whether each photon is isolated
      std::vector< PhotonCandidate > phoCand;//PhotonCandidate defined in RazorAuxPhoton.hh
      int nPhotonsAbove40GeV = 0;
      for(int i = 0; i < nPhotons; i++)
      {
        //ECAL Prefiring Recipe
        if( phoPt[i] > 50. && fabs( pho_superClusterEta[i] ) > 2.25 && fabs ( pho_superClusterEta[i] ) < 3.0 ) ecal_prefiring_affected = true;
        //if( phoPt[i] > 50. && fabs( pho_superClusterEta[i] ) > 2.25 && fabs ( pho_superClusterEta[i] ) < 3.0 ) continue;
        if(_ecaldebug) std::cout << "[ECAL Prefiring: not affected \n"<< std::endl;
        //if( phoPt[i] > 50. && fabs( phoEta[i] ) > 2.25 && fabs ( phoEta[i] ) < 3.0 ) continue;

        if ( (pho_seedRecHitSwitchToGain6[i] ||
              pho_seedRecHitSwitchToGain1[i] ||
              pho_anyRecHitSwitchToGain6[i] ||
              pho_anyRecHitSwitchToGain1[i]
            ))
            {
              Flag_hasEcalGainSwitch = true;
            }

        //R9 requirement, only R9<0.5 photons are allow to be candidates
        //if( phoR9[i] < 0.5 ) continue;
        if ( _phodebug ) std::cout << " pho_P9: " << phoR9[i] << std::endl;

        //ID cuts -- apply isolation after candidate pair selection
        if ( _phodebug ) std::cout << "pho# " << i << " phopt1: " << phoPt[i] << " pho_eta: " << phoEta[i] << std::endl;
        if ( _phodebug ) std::cout << "pho# " << i << " phopt1: " << phoPt[i] << " pho_sc_eta: " << pho_superClusterEta[i] << std::endl;
        //--loose ID (default)--
        if (doRequireID)
        {
          if (analysisTag == "Razor2017_31Mar2018Rereco" )
          {
            if ( !photonPassLooseIDWithoutEleVeto_94X(i) )
            {
              if ( _phodebug ) std::cout << "[DEBUG]: failed run2 ID 2017 94X" << std::endl;
              continue;
            }
          }
          else if (analysisTag == "Razor2017_92X" || analysisTag == "Razor2017_17Nov2017Rereco" )
          {
            if ( !photonPassLooseIDWithoutEleVeto_2017(i) )
            {
              if ( _phodebug ) std::cout << "[DEBUG]: failed run2 ID 2017 92X" << std::endl;
              continue;
            }
          }
          else
          {
            if ( !photonPassLooseIDWithoutEleVeto(i) )
            {
              if ( _phodebug ) std::cout << "[DEBUG]: failed run2 ID 2016" << std::endl;
              continue;
            }
          }
        if ( _phodebug ) std::cout << "pho# " << i << " passes default loose ID" << std::endl;
        }

        //---------------
	      //----tight ID---
        //---------------
        if (doRequireTightID)
        {
          if (analysisTag == "Razor2017_31Mar2018Rereco" )
          {
            if ( !photonPassTightIDWithoutEleVeto_94X(i) )
            {
              if ( _phodebug ) std::cout << "[DEBUG]: failed run2 Tight ID 2017 94X" << std::endl;
              continue;
            }
          }
          else if (analysisTag == "Razor2017_92X" || analysisTag == "Razor2017_17Nov2017Rereco" )
          {
            if ( !photonPassTightIDWithoutEleVeto_2017(i) )
            {
              if ( _phodebug ) std::cout << "[DEBUG]: failed run2 Tight ID 2017 92X" << std::endl;
              continue;
            }
          }
          else
          {
            if ( !photonPassTightIDWithoutEleVeto(i) )
            {
              if ( _phodebug ) std::cout << "[DEBUG]: failed run2 Tight ID 2016" << std::endl;
              continue;
            }
          }
        if ( _phodebug ) std::cout << "pho# " << i << " passes requied tight ID" << std::endl;
        }

        //**********************************************************
        //Isolation, electron veto, and Barrel requirements are introduced here
        //if we want to use the "regular" selection sequence
        //**********************************************************

        //---------------
        //ONLY EB photons
        //---------------
        if (!(fabs(pho_superClusterEta[i]) < 1.4442 )) continue;
        if ( _phodebug ) std::cout << "pho# " << i << " phopt1: " << phoPt[i] << " pho_sc_eta: " << pho_superClusterEta[i] << std::endl;
        if ( _phodebug ) std::cout << "pho# " << i << " passes EB Eta" << std::endl;
        //-------------
        // --ele Veto--
        //-------------
	      if (doEleVeto) if (!(pho_passEleVeto[i])) continue;
        if ( _phodebug ) std::cout << "pho# " << i << " passes Ele Veto" << std::endl;
        //-------------------
        //--loose isolation--
        //-------------------
        if (doRequireIso)
        {
          if (analysisTag == "Razor2017_31Mar2018Rereco" )
          {
            if (!(photonPassLooseIso_94X(i))) continue;
          }
          else if (analysisTag == "Razor2017_92X" || analysisTag == "Razor2017_17Nov2017Rereco" )
          {
            if (!(photonPassLooseIso_2017(i))) continue;
          }
          else
          {
            if (!(photonPassLooseIso(i))) continue;
          }
        if ( _phodebug ) std::cout << "pho# " << i << " passes Iso" << std::endl;
        }
        //--------------------
        //--tight isolation---
        //--------------------
        if (doRequireTightIso)
        {
          if (analysisTag == "Razor2017_31Mar2018Rereco" )
          {
            if (!(photonPassTightIso_94X(i))) continue;
          }
          else if (analysisTag == "Razor2017_92X" || analysisTag == "Razor2017_17Nov2017Rereco" )
          {
            if (!(photonPassTightIso_2017(i))) continue;
          }
          else
          {
            if (!(photonPassTightIso(i))) continue;
          }
        if ( _phodebug ) std::cout << "pho# " << i << " passes tight Iso" << std::endl;
        }
/*
        //*****************************************************************************
        //Photons must be separated from any selected leptons (muons and electrons)
        //*****************************************************************************
        //Remove muon overlaps
        bool overlapm = false;
        for(int j = 0; j < int(GoodMuons.size()); j++)
        {
          TLorentzVector mu = GoodMuons.at(j);
          if (RazorAnalyzer::deltaR(phoEta[i],phoPhi[i],mu.Eta(),mu.Phi()) < 0.5)  overlapm = true;
        }
        if (overlapm) continue;//removing muon overlaps
        if ( _phodebug ) std::cout << "pho# " << i << " passes muon overlap" << std::endl;
        //remove electron overlaps
        bool overlape = false;
        for(int k = 0; k < int(GoodElectrons.size()); k++)
        {
          TLorentzVector ele = GoodElectrons.at(k);
          if (RazorAnalyzer::deltaR(phoEta[i],phoPhi[i],ele.Eta(),ele.Phi()) < 1.0) overlape = true;
        }
        if( doEleVeto && overlape ) continue;//removing electron overlaps only when we don't want Zee candidates
        if ( _phodebug ) std::cout << "pho# " << i << " passes electron overlap" << std::endl;
*/

        //----------------------------------------------
        //Get Scale and Define Corrected Photon momentum
        //----------------------------------------------
        float pho_pt_corr = phoPt[i];
        double scale = 1;
        double smear = 0;
        if ( doPhotonScaleCorrection )
        {
          if ( analysisTag != "Razor2016_03Feb2017Rereco" && analysisTag != "Razor2017_92X" &&  analysisTag != "Razor2017_17Nov2017Rereco" && analysisTag != "Razor2017_31Mar2018Rereco" )
          //if ( analysisTag != "Razor2017_92X" &&  analysisTag != "Razor2017_17Nov2017Rereco"  &&  analysisTag != "Razor2016_MoriondRereco"  )
          {
            scale = photonCorrector->ScaleCorrection(run, (fabs(pho_superClusterEta[i]) < 1.5),
                                                     phoR9[i], pho_superClusterEta[i], phoE[i]/cosh(pho_superClusterEta[i]));
            smear = photonCorrector->getSmearingSigma(run, (fabs(pho_superClusterEta[i]) < 1.5),
                                                      phoR9[i], pho_superClusterEta[i], phoE[i]/cosh(pho_superClusterEta[i]), 0., 0.);
          }
          else
          {
            //old version (92X)
            //scale = photonCorrector_2017->ScaleCorrection_class_2017(run, (fabs(pho_superClusterEta[i]) < 1.5),
            //                                              phoR9[i], pho_superClusterEta[i], phoE[i]/cosh(pho_superClusterEta[i]));
            //smear = photonCorrector_2017->smearingSigma(run, (fabs(pho_superClusterEta[i]) < 1.5),
             //                                              phoR9[i], pho_superClusterEta[i], phoE[i]/cosh(pho_superClusterEta[i]), 0, 0., 0.);
            
             const EnergyScaleCorrection_class_2017::ScaleCorrection_class_2017* scaleCorr = photonCorrector_2017->EnergyScaleCorrection_class_2017::getScaleCorr(run, phoE[i]/cosh(pho_superClusterEta[i]), pho_superClusterEta[i], phoR9[i], 12);
            const EnergyScaleCorrection_class_2017::SmearCorrection_class_2017* smearCorr = photonCorrector_2017->EnergyScaleCorrection_class_2017::getSmearCorr(run, phoE[i]/cosh(pho_superClusterEta[i]), pho_superClusterEta[i], phoR9[i], 12);
            if(scaleCorr!=NULL) scale  = scaleCorr->scale();
            if(smearCorr!=NULL) smear  = smearCorr->sigma(phoE[i]/cosh(pho_superClusterEta[i]));
           
          }
          //apply scale to data and smearing to MC
	        if (isData)
          {
            pho_pt_corr = phoPt[i]*scale;
            /*std::cout << run << " " << (fabs(pho_superClusterEta[i]) < 1.5) << " " <<  phoR9[i] << " "
			    << pho_superClusterEta[i] << " " <<  phoE[i]/cosh(pho_superClusterEta[i]) << " "
			    << scale << std::endl;
           */
           if (_phodebug && scale != 1.0 ) std::cout << "[DEBUG] : Photon Energy Scale Corrections: "
							    << phoPt[i] << " * " << scale << " --> " << pho_pt_corr << "\n";
          }
          else
          {
            pho_pt_corr = phoPt[i]*(1+smear*random.Gaus());
          }
        if ( _phodebug ) std::cout << "pho# " << i << " finishes photon energy corr" << std::endl;
        if ( _phodebug ) std::cout << "pho# " << i << " phopt1: " << phoPt[i] << " pho_sc_eta: " << pho_superClusterEta[i] << std::endl;
        }

        //---------------------------------------------------------
        //Set Corrected 4-Momentum for photon candidate with mass=0
        //---------------------------------------------------------
        TVector3 vec;
        vec.SetPtEtaPhi( pho_pt_corr, phoEta[i], phoPhi[i] );
    	  TLorentzVector thisPhoton;
    	  thisPhoton.SetVectM( vec, .0 );

        if ( phoPt[i] < 20.0 )
        {
          if ( _phodebug ) std::cout << "[DEBUG]: failed pt" << std::endl;
          continue;
        }
        //removing gap photons
        if ( fabs(pho_superClusterEta[i]) > 1.4442 && fabs(pho_superClusterEta[i]) < 1.566 )
        {
          if ( _phodebug ) std::cout << "[INFO]: failed gap" << std::endl;
          continue;
        }

        //count number of photons above 40 GeV
        if( phoPt[i] > 40.0 ) nPhotonsAbove40GeV++;
        //-----------------------------
        //uncorrected photon 4-momentum
        //-----------------------------
        TVector3 vtx( pvX, pvY, pvZ );
        TVector3 phoPos;
        if ( fabs( pho_superClusterEta[i] ) < 1.479 )
        {
          phoPos.SetXYZ( EB_R*cos( pho_superClusterPhi[i]), EB_R*sin( pho_superClusterPhi[i] ), EB_R*sinh( pho_superClusterEta[i] ) );
        }
        else
        {
          double R = fabs( EE_Z/sinh( pho_superClusterEta[i] ) );
          if ( pho_superClusterEta[i] > .0 )
          {
            phoPos.SetXYZ( R*cos( pho_superClusterPhi[i] ), R*sin( pho_superClusterPhi[i] ), EE_Z);
          }
          else
          {
            phoPos.SetXYZ( R*cos( pho_superClusterPhi[i] ), R*sin( pho_superClusterPhi[i] ), -EE_Z);
          }
        }
        //set super-cluster 4-momentum
        TLorentzVector phoSC = GetCorrectedMomentum( vtx, phoPos, pho_RegressionE[i] );

        //------------------------
        //Filling Photon Candidate
        //------------------------
        PhotonCandidate tmp_phoCand;
        tmp_phoCand.Index = i;
        tmp_phoCand.photon = thisPhoton;
        tmp_phoCand.photonSC = phoSC;
        tmp_phoCand.scEta = pho_superClusterEta[i];
        if ( _phodebug ) std::cout << "[INFO] : pho# " << i << " pho_sc_eta: " << pho_superClusterEta[i] << ", tmp_phoCand.scEta = " << tmp_phoCand.scEta << std::endl;
        tmp_phoCand.scPhi = pho_superClusterPhi[i];
        tmp_phoCand.SigmaIetaIeta = phoFull5x5SigmaIetaIeta[i];
        tmp_phoCand.R9 = phoR9[i];
        tmp_phoCand.HoverE = pho_HoverE[i];
        tmp_phoCand.sumChargedHadronPt = pho_pfIsoChargedHadronIso[i];
        tmp_phoCand.sumNeutralHadronEt = pho_pfIsoNeutralHadronIso[i];
        tmp_phoCand.sumPhotonEt = pho_pfIsoPhotonIso[i];
        tmp_phoCand.sigmaEOverE = pho_RegressionEUncertainty[i]/pho_RegressionE[i];
        tmp_phoCand._passEleVeto = pho_passEleVeto[i];
        tmp_phoCand.scale = scale;
	if (analysisTag == "Razor2017_31Mar2018Rereco" )
	{
		tmp_phoCand._passIso = photonPassLooseIso_94X(i);
	}
	else if (analysisTag == "Razor2017_92X" || analysisTag == "Razor2017_17Nov2017Rereco" )
	{
		tmp_phoCand._passIso = photonPassLooseIso_2017(i);
	}
	else
	{
		tmp_phoCand._passIso = photonPassLooseIso(i);
	}
        phoCand.push_back( tmp_phoCand );
        nSelectedPhotons++;
        GoodPhotons.push_back(thisPhoton);
      }//end of loop over photons

      //--------------------------------------
      //Require at least two photon candidates
      //--------------------------------------
      if ( phoCand.size() < 2 ) continue;

      if ( _debug ) std::cout << "[DEBUG]: nphotons--> " << phoCand.size()
			      << " " << nSelectedPhotons << std::endl;

      //----------------------------------------
      //find the "best" photon pair, highest sum Pt!
      //----------------------------------------
      TLorentzVector HiggsCandidate(0,0,0,0);
      TLorentzVector HiggsCandidateSC(0,0,0,0);
      int HiggsPhoIndex1 = -1;
      int HiggsPhoIndex2 = -1;
      double bestSumPt = -99.;
      std::vector< PhotonCandidate > phoSelectedCand;
      PhotonCandidate bestCand[2];

      for ( size_t i = 0; i < phoCand.size(); i++ )
      {
        for ( size_t j = i+1; j < phoCand.size(); j++ )
        {
          PhotonCandidate pho1 = phoCand[i];
          PhotonCandidate pho2 = phoCand[j];
          double diphotonMass = (pho1.photon + pho2.photon).M();
          if ( _debug ) std::cout << "[DEBUG] Diphoton Sum pT: " << pho1.photon.Pt() + pho2.photon.Pt() << std::endl;
          //discard candidate if mgg < 50 GeV
          if( diphotonMass < 50 )
          {
            if ( _debug ) std::cout << "[DEBUG]: Diphoton mass < 50 GeV: mgg-> " << diphotonMass << std::endl;
            if ( _debug ) std::cout << "... pho1Pt: " << pho1.photon.Pt()  << " pho2Pt: " << pho2.photon.Pt()  << std::endl;
            continue;
          }
          //---------------------------------------------
          //if the sum of the photon pT's is larger than
	        //that of the current Higgs candidate,
	        //make this the Higgs candidate
	        //---------------------------------------------
          if( pho1.photon.Pt() + pho2.photon.Pt() > bestSumPt )
          {
            bestSumPt = pho1.photon.Pt() + pho2.photon.Pt();
            pTGammaGamma = pho1.photon.Pt() + pho2.photon.Pt();
            HiggsCandidate = pho1.photon + pho2.photon;
            HiggsCandidateSC = pho1.photonSC + pho2.photonSC;
            if ( pho1.photon.Pt() >= pho2.photon.Pt() )
            {
              bestCand[0] = pho1;
              bestCand[1] = pho2;
              HiggsPhoIndex1 = pho1.Index;
              HiggsPhoIndex2 = pho2.Index;
            }
            else
            {
              bestCand[0] = pho2;
              bestCand[1] = pho1;
              HiggsPhoIndex1 = pho2.Index;
              HiggsPhoIndex2 = pho1.Index;
            }
          }//best pt if
        }//loop j-th photon
      }//loop i-th photon

      //---------------------------------------
      //just use this container for convenience
      //to parse the data into TTree
      //---------------------------------------
      //bestCand[0] is the leading photon
      //bestCand[1] is the subleading photon
      //if( bestCand[0].photon.Pt()/HiggsCandidate.M() < 1./3. || bestCand[1].photon.Pt()/HiggsCandidate.M() < 1./3. ) continue;
      //if( bestCand[0].photon.Pt()/HiggsCandidate.M() < 1./4. && bestCand[1].photon.Pt()/HiggsCandidate.M() < 1./4. ) continue;
      //if( bestCand[0].photon.Pt()/HiggsCandidate.M() < 1./3. ) continue;
      //if( bestCand[1].photon.Pt()/HiggsCandidate.M() < 1./4. ) continue;
      if ( _phodebug ) std::cout << "pho PT :  " << bestCand[0].photon.Pt() << " and  " << bestCand[1].photon.Pt()
                                 << " Mgg : " << HiggsCandidate.M() << " pt1/Mgg = "   << bestCand[0].photon.Pt()/HiggsCandidate.M()
                                 << " pt2/Mgg = " << bestCand[1].photon.Pt()/HiggsCandidate.M() << std::endl;

      phoSelectedCand.push_back(bestCand[0]);
      phoSelectedCand.push_back(bestCand[1]);
      if ( _phodebug ) std::cout << " " << std::endl;
      //-----------------------------------
      //Filling Selected Photon Information
      //-----------------------------------
      TLorentzVector pho_cand_vec[2];
      int _pho_index = 0;
      for ( auto& tmpPho : phoSelectedCand )
      {
        if ( !( tmpPho.Index == HiggsPhoIndex1 || tmpPho.Index == HiggsPhoIndex2 ) ) continue;
        if( _pho_index > 1 ) std::cerr << "[ERROR]: Photon index larger than 1!" << std::endl;
        pho_cand_vec[_pho_index]           = tmpPho.photon;
        Pho_E[_pho_index]                  = tmpPho.photon.E();
        Pho_Pt[_pho_index]                 = tmpPho.photon.Pt();
        Pho_Eta[_pho_index]                = tmpPho.photon.Eta();
        Pho_Phi[_pho_index]                = tmpPho.photon.Phi();
        PhoSC_E[_pho_index]                = tmpPho.photonSC.E();
        PhoSC_Pt[_pho_index]               = tmpPho.photonSC.Pt();
        //PhoSC_Eta[_pho_index]              = tmpPho.photonSC.Eta();
        //if ( _phodebug ) std::cout << " tmpPho.photonSC.Eta(): " << tmpPho.photonSC.Eta() << ", PhoSC_Eta[_pho_index] : " << PhoSC_Eta[_pho_index] << std::endl;
        PhoSC_Eta[_pho_index]              = tmpPho.scEta;
        if ( _phodebug ) std::cout << " tmpPho.scEta: " << tmpPho.scEta << ", PhoSC_Eta[_pho_index] : " << PhoSC_Eta[_pho_index] << std::endl;
        //PhoSC_Phi[_pho_index]              = tmpPho.photonSC.Phi();
        PhoSC_Phi[_pho_index]              = tmpPho.scPhi;
        Pho_SigmaIetaIeta[_pho_index]      = tmpPho.SigmaIetaIeta;
        Pho_R9[_pho_index]                 = tmpPho.R9;
        Pho_HoverE[_pho_index]             = tmpPho.HoverE;
        Pho_sumChargedHadronPt[_pho_index] = tmpPho.sumChargedHadronPt;
        Pho_sumNeutralHadronEt[_pho_index] = tmpPho.sumNeutralHadronEt;
        Pho_sumPhotonEt[_pho_index]        = tmpPho.sumPhotonEt;
        Pho_sigmaEOverE[_pho_index]        = tmpPho.sigmaEOverE;
        Pho_passEleVeto[_pho_index]        = tmpPho._passEleVeto;
        Pho_passIso[_pho_index]            = tmpPho._passIso;
        Pho_scale[_pho_index]              = tmpPho.scale;
        _pho_index++;
      }

      //removing events with less than two good photon candidates
      if ( _pho_index < 2 ) continue;

      if ( _debug )
      {
        std::cout << "[DEBUG]: best photon pair: "
		    << "\n-> pho1Pt: " << Pho_Pt[0]
		    << "\n-> pho2Pt: " << Pho_Pt[1]
		    << std::endl;
      }


      //record higgs candidate info
      mGammaGamma    = HiggsCandidate.M();
      pTGammaGamma   = HiggsCandidate.Pt();
      mGammaGammaSC  = HiggsCandidateSC.M();
      pTGammaGammaSC = HiggsCandidateSC.Pt();
      if ( _debug ) std::cout << "[DEBUG]: mgg-> " << mGammaGamma << " pTgg->" << pTGammaGamma << std::endl;

      //******************************************************
      //compute trigger efficiency weights for MC
      //******************************************************
      triggerEffWeight = 1.0;
      triggerEffSFWeight = 1.0;
      double leadPhoPt=0;
      double leadPhoEta=0;
      double trailingPhoPt=0;
      double trailingPhoEta=0;
      if (Pho_Pt[0] > Pho_Pt[1])
      {
        leadPhoPt = Pho_Pt[0];
        leadPhoEta = Pho_Eta[0];
        trailingPhoPt = Pho_Pt[1];
        trailingPhoEta= Pho_Eta[1];
      }
      else
      {
        leadPhoPt = Pho_Pt[1];
        leadPhoEta = Pho_Eta[1];
        trailingPhoPt = Pho_Pt[0];
        trailingPhoEta= Pho_Eta[0];
      }
      //------------------------------
      //trigger efficiencies
      //------------------------------
      double triggerEffLeadingLeg = helper->getDiphotonTrigLeadingLegEff( leadPhoPt, leadPhoEta );
      double triggerEffTrailingLeg = helper->getDiphotonTrigTrailingLegEff( trailingPhoPt, trailingPhoEta );
      triggerEffWeight = triggerEffLeadingLeg*triggerEffTrailingLeg;
      double triggerEffSFLeadingLeg = helper->getDiphotonTrigLeadingLegEffSF( leadPhoPt, leadPhoEta );
      double triggerEffSFTrailingLeg = helper->getDiphotonTrigTrailingLegEffSF( trailingPhoPt, trailingPhoEta );
      triggerEffSFWeight = triggerEffSFLeadingLeg*triggerEffSFTrailingLeg;
      //******************************************************
      //compute photon efficiency scale factor
      //******************************************************
      if ( analysisTag == "Razor2017_92X" || analysisTag == "Razor2017_17Nov2017Rereco" || analysisTag == "Razor2017_31Mar2018Rereco" )
      {

        pho1EffSF *= helper->getPhotonScaleFactor(leadPhoPt, leadPhoEta, true);
        pho1EffSFUnc += helper->getPhotonScaleFactorError(leadPhoPt, leadPhoEta, true);
        pho1EffSFUp *= pho1EffSF + pho1EffSFUnc; 
        pho1EffSFDown *= pho1EffSF - pho1EffSFUnc; 

        pho2EffSF *= helper->getPhotonScaleFactor(trailingPhoPt, trailingPhoEta, true);
        pho2EffSFUnc += helper->getPhotonScaleFactorError(trailingPhoPt, trailingPhoEta, true);
        pho2EffSFUp *= pho2EffSF + pho2EffSFUnc; 
        pho2EffSFDown *= pho2EffSF - pho2EffSFUnc; 

        photonEffSF *= pho1EffSF*pho2EffSF;
	photonEffSFUp *= fmax(pho1EffSFUp, pho1EffSFDown)*fmax(pho2EffSFUp, pho2EffSFDown);
	photonEffSFDown *= fmin(pho1EffSFUp, pho1EffSFDown)*fmin(pho2EffSFUp, pho2EffSFDown);
	photonEffSFSys += photonEffSF*sqrt( pow(pho1EffSFUnc/pho1EffSF, 2.) + pow(pho2EffSFUnc/pho2EffSF, 2.) );

        /*photonEffSF = helper->getPhotonScaleFactor(leadPhoPt, leadPhoEta, true) *
                      helper->getPhotonScaleFactor(trailingPhoPt, trailingPhoEta, true);*/
      }
      else
      {

        pho1EffSF *= helper->getPhotonScaleFactor(leadPhoPt, leadPhoEta);
        pho1EffSFUnc += helper->getPhotonScaleFactorError(leadPhoPt, leadPhoEta);
        pho1EffSFUp *= pho1EffSF + pho1EffSFUnc; 
        pho1EffSFDown *= pho1EffSF - pho1EffSFUnc; 

        pho2EffSF *= helper->getPhotonScaleFactor(trailingPhoPt, trailingPhoEta);
        pho2EffSFUnc += helper->getPhotonScaleFactorError(trailingPhoPt, trailingPhoEta);
        pho2EffSFUp *= pho2EffSF + pho2EffSFUnc; 
        pho2EffSFDown *= pho2EffSF - pho2EffSFUnc; 

        photonEffSF *= pho1EffSF*pho2EffSF;
	photonEffSFUp *= fmax(pho1EffSFUp, pho1EffSFDown)*fmax(pho2EffSFUp, pho2EffSFDown);
	photonEffSFDown *= fmin(pho1EffSFUp, pho1EffSFDown)*fmin(pho2EffSFUp, pho2EffSFDown);
	photonEffSFSys += photonEffSF*sqrt( pow(pho1EffSFUnc/pho1EffSF, 2.) + pow(pho2EffSFUnc/pho2EffSF, 2.) );

        /*photonEffSF = helper->getPhotonScaleFactor(leadPhoPt, leadPhoEta) *
                      helper->getPhotonScaleFactor(trailingPhoPt, trailingPhoEta);*/
      }

      if (isFastsimSMS)
      {

        pho1EffSF *= helper->getPhotonFastsimToFullsimScaleFactor(leadPhoPt, leadPhoEta);
        pho1EffSFUnc += helper->getPhotonFastsimToFullsimScaleFactorError(leadPhoPt, leadPhoEta);
        pho1EffSFUp *= pho1EffSF + pho1EffSFUnc; 
        pho1EffSFDown *= pho1EffSF - pho1EffSFUnc; 

        pho2EffSF *= helper->getPhotonFastsimToFullsimScaleFactor(trailingPhoPt, trailingPhoEta);
        pho2EffSFUnc += helper->getPhotonFastsimToFullsimScaleFactorError(trailingPhoPt, trailingPhoEta);
        pho2EffSFUp *= pho2EffSF + pho2EffSFUnc; 
        pho2EffSFDown *= pho2EffSF - pho2EffSFUnc; 

        photonEffSF *= pho1EffSF*pho2EffSF;
	photonEffSFUp *= fmax(pho1EffSFUp, pho1EffSFDown)*fmax(pho2EffSFUp, pho2EffSFDown);
	photonEffSFDown *= fmin(pho1EffSFUp, pho1EffSFDown)*fmin(pho2EffSFUp, pho2EffSFDown);
	photonEffSFSys += 0.; //negliegible sys from fastsim to fullsim

        /*photonEffSF *= helper->getPhotonFastsimToFullsimScaleFactor(leadPhoPt, leadPhoEta) *
                       helper->getPhotonFastsimToFullsimScaleFactor(trailingPhoPt, trailingPhoEta);*/
      }

      //***********************************************************
      //get mother ID of photons
      //***********************************************************
      //cout << "Photon1 : " << Pho_Pt[0] << " " << Pho_Eta[0] << " " << Pho_Phi[0] << "\n";
      for(int g = 0; g < nGenParticle; g++)
      {
        if (!(deltaR(gParticleEta[g] , gParticlePhi[g], Pho_Eta[0],Pho_Phi[0]) < 0.5) ) continue;
        // status = 22 for Higgs bosons in MadGraph/Pythia8
        //if(gParticleStatus[g] != 1) continue;
        if(gParticleId[g] != 22) continue;
        if(!( (gParticleStatus[g] == 1 && gParticleMotherId[g] != 22) || gParticleStatus[g] == 22 || gParticleStatus[g] == 23)) continue;
        Pho_motherID[0] = gParticleMotherId[g];
	    }
      for(int g = 0; g < nGenParticle; g++)
      {
        if (!(deltaR(gParticleEta[g] , gParticlePhi[g], Pho_Eta[1],Pho_Phi[1]) < 0.5) ) continue;
        // status = 22 for Higgs bosons in MadGraph/Pythia8
        if(gParticleId[g] != 22) continue;
        if(!( (gParticleStatus[g] == 1 && gParticleMotherId[g] != 22) || gParticleStatus[g] == 22 || gParticleStatus[g] == 23)) continue;
        Pho_motherID[1] = gParticleMotherId[g];
      }

      //--------------
      //good muon selection
      //--------------
      string muonEraName = "";
      if (dataset == "94X") muonEraName = "2017_94X";
      else if (dataset == "80X") muonEraName = "Spring15";
      vector<TLorentzVector> GoodMuons;
      std::vector< MuonCandidate > muCand;
      for( int i = 0; i < nMuons; i++ )
	{
          //TLorentzVector for this muon
          TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]);

	  double dr = 10.0/fmin(fmax(muonPt[i], 50.0),200.0);
	  //use Loose Muon POG ID and |d0|<0.2 and |dZ|<0.5 && miniiso / pt < 0.2
	  if(!(
	       muonIsLoose[i] && fabs(muon_d0[i]) < 0.2 && fabs(muon_dZ[i]) < 0.5
	       && 
	       ((muon_chargedMiniIso[i] + fmax(0.0, muon_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetMuonEffectiveArea90(i,muonEraName)*pow(dr/0.3,2)) )/muonPt[i] < 0.2)
	       )) continue;
	     
	  if(muonPt[i] < 20) continue;
	  if(abs(muonEta[i]) > 2.4) continue;
        //*****************************************************************************
        //Photons must be separated from any selected leptons (muons and electrons)
        //*****************************************************************************
        //Remove muon overlaps
        bool overlapm = false;
        for(int j = 0; j < 2; j++)
        {
          TLorentzVector pho = pho_cand_vec[j];
          if (RazorAnalyzer::deltaR(pho.Eta(),pho.Phi(),thisMuon.Eta(),thisMuon.Phi()) <= 0.5)  overlapm = true;
        }
        if (overlapm) continue;//removing muon overlaps
        if ( _phodebug ) std::cout << "mu# " << i << " passes photon overlap" << std::endl;
	  nLooseMuons++;
          GoodMuons.push_back(thisMuon);
	  if( isTightMuon(i) ) nTightMuons++;
	  //Filling Muon Candidate
	  MuonCandidate tmp_muCand;
	  tmp_muCand.Index = i;
	  tmp_muCand.muon = thisMuon;
	  tmp_muCand.muonCharge = muonCharge[i];
	  tmp_muCand.isTightMuon = isTightMuon(i);
          muCand.push_back( tmp_muCand );
	}

      //------------------
      //good electron selection
      //------------------
      string electronEraName = "";
      if (dataset == "94X") electronEraName = "2017_94X";
      else if (dataset == "80X") electronEraName = "Spring15";
      vector<TLorentzVector> GoodElectrons;
      std::vector< ElectronCandidate > eleCand;
      for( int i = 0; i < nElectrons; i++ )
	{
        //*****************************************************************************
        //Photons must be separated from any selected leptons (muons and electrons)
        //*****************************************************************************
        //remove electron overlaps
        bool overlape = false;
        for(int k = 0; k < 2; k++)
        {
          TLorentzVector pho = pho_cand_vec[k];
          if (RazorAnalyzer::deltaR(pho.Eta(),pho.Phi(),eleEta[i],elePhi[i]) <= 1.0)  overlape = true;
        }
        if( doEleVeto && overlape ) continue;//removing electron overlaps only when we don't want Zee candidates
        if ( _phodebug ) std::cout << "ele# " << i << " passes photon overlap" << std::endl;
          //Remove overlaps
          bool overlap = false;
          for(int j = 0; j < int(GoodMuons.size()); j++){
                  TLorentzVector mu = GoodMuons.at(j);
                  if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],mu.Eta(),mu.Phi()) <= 0.4)  overlap = true;
          }
          if (overlap) continue;
        if ( _phodebug ) std::cout << "ele# " << i << " passes muon overlap" << std::endl;
          //TLorentzVector for this electron
          TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);

	  //DR definition for mini-isolation
	  double dr = fmax(0.05,fmin(0.2, 10/elePt[i]));
        if ( _phodebug ) std::cout << "ele# " << i << " dr = " << dr << std::endl;

        if ( _phodebug ) std::cout << "ele# " << i << " id = " << passEGammaPOGLooseElectronID(i,true,electronEraName) << std::endl;
        if ( _phodebug ) std::cout << "ele# " << i << " iso = " << (ele_chargedMiniIso[i] +  fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveArea90(i,electronEraName)*pow(dr/0.3,2)))/elePt[i]  << std::endl;
	  if(!(passEGammaPOGLooseElectronID(i,true,electronEraName) && 
	       ((ele_chargedMiniIso[i] +  fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveArea90(i,electronEraName)*pow(dr/0.3,2)))/elePt[i] < 0.1) 
		)) continue;
        if ( _phodebug ) std::cout << "ele# " << i << " passes id " << std::endl;
	  //if ( !( passMVALooseElectronID(i,dataset) && passMVANonTrigVetoElectronIso(i) && fabs(ele_ip3dSignificance[i]) < 4. ) ) continue;//Only for electron WP test
	  if( elePt[i] < 20 ) continue;
	  if( abs(eleEta[i]) > 2.4 ) continue;
        if ( _phodebug ) std::cout << "ele# " << i << " passes pt and eta " << std::endl;
	  nLooseElectrons++;
        if ( _phodebug ) std::cout << "ele# " << i << " added to nLooseEle " << std::endl;
          GoodElectrons.push_back(thisElectron);      	  
	  if(passEGammaPOGTightElectronID(i,true,electronEraName) && 
	     ((ele_chargedMiniIso[i] +  fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveArea90(i,electronEraName)*pow(dr/0.3,2)))/elePt[i] < 0.1) 
	     ) nTightElectrons++;
	  //Filling Electron Candidate
	  ElectronCandidate tmp_eleCand;
	  tmp_eleCand.Index = i;
	  tmp_eleCand.electron = thisElectron;
	  tmp_eleCand.eleCharge = eleCharge[i];
	  tmp_eleCand.isTightElectron = (passEGammaPOGTightElectronID(i,true,electronEraName) && 
					 ((ele_chargedMiniIso[i] +  fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - fixedGridRhoFastjetAll*GetElectronEffectiveArea90(i,electronEraName)*pow(dr/0.3,2)))/elePt[i] < 0.1));
          eleCand.push_back( tmp_eleCand );
	}

      //------------
      //BTagged Jets
      //------------
      vector<TLorentzVector> GoodBJets;
      std::vector< BjetCandidate > bjetCand;
      vector<bool> GoodBJetsIsCVSL;
      vector<bool> GoodBJetsIsCVSM;
      vector<bool> GoodBJetsIsCVST;
      vector< pair<TLorentzVector, bool> > GoodCSVLBJets; //contains CSVL jets passing selection.  The bool is true if the jet passes CSVM, false if not
      std::cout << "[INFO]: begin b jet loop \n"<< std::endl;
      for(int i = 0; i < nJets; i++)
      {
        //Jet energy Corrections
        double JEC = JetEnergyCorrectionFactor( jetPt[i], jetEta[i], jetPhi[i], jetE[i],
                                                fixedGridRhoAll, jetJetArea[i], runNum,
                                                JetCorrectorIOV, JetCorrector );
        TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );
        if (_jetdebug) std::cout << "jet no. " << i << " raw info (pt, eta, phi, E):" << jetPt[i] << " , " << jetEta[i] << " , " << jetPhi[i] << " , " << jetE[i] << " ;\n"<< std::endl;
        if (_jetdebug) std::cout << "jet no. " << i << " jec info (pt, eta, phi, E):" << jetPt[i]*JEC << " , " << jetEta[i] << " , " << jetPhi[i] << " , " << jetE[i]*JEC << " ;\n"<< std::endl;
        //ECAL Prefiring Recipe
        if( thisJet.Pt() > 100. && fabs( thisJet.Eta() ) > 2.25 && fabs ( thisJet.Eta() ) < 3.0 ) ecal_prefiring_affected = true;
        //if( thisJet.Pt() > 100. && fabs( thisJet.Eta() ) > 2.25 && fabs ( thisJet.Eta() ) < 3.0 ) continue;
        if(_ecaldebug) std::cout << "[ECAL Prefiring: not affected \n"<< std::endl;
        
	if( thisJet.Pt() < BJET_CUT ) continue;//According to the April 1st 2015 A
        if (_jetdebug) std::cout << "jet no. " << i << " passed b-jet pt 20 cut ;\n"<< std::endl;
        if( fabs( thisJet.Eta() ) >= 2.4 ) continue;
        if (_jetdebug) std::cout << "jet no. " << i << " passed b-jet eta 2.4 cut ;\n"<< std::endl;
        if (!isFastsimSMS)
        {
          if ( !jetPassIDLoose[i] ) continue;
        }
        //Remove muon from jet collection
        bool overlapjm = false;
        for(int j = 0; j < int(GoodMuons.size()); j++)
        {
          TLorentzVector mu = GoodMuons.at(j);
          if (RazorAnalyzer::deltaR( thisJet.Eta(), thisJet.Phi(), mu.Eta(), mu.Phi()) < 0.4 ) overlapjm = true;
        }
        if(overlapjm) continue;
        if (_jetdebug) std::cout << "jet no. " << i << " passed b-jet muon cleaning ;\n"<< std::endl;

        //remove electrons from jet collection
        bool overlapje = false;
        for(int k = 0; k < int(GoodElectrons.size()); k++)
        {
          TLorentzVector ele = GoodElectrons.at(k);
          if (RazorAnalyzer::deltaR( thisJet.Eta(), thisJet.Phi(), ele.Eta(), ele.Phi()) < 0.4 ) overlapje = true;
        }
        if(overlapje) continue;
        if (_jetdebug) std::cout << "jet no. " << i << " passed b-jet ele cleaning ;\n"<< std::endl;

        //exclude higgs-candidate photons from the jet collection
        double deltaRJetPhoton = min( thisJet.DeltaR( pho_cand_vec[0] ), thisJet.DeltaR( pho_cand_vec[1] ) );
        if ( deltaRJetPhoton <= 0.4 ) continue;//According to the April 1st 2015 AN
        if (_jetdebug) std::cout << "jet no. " << i << " passed b-jet photon cleaning ;\n"<< std::endl;
        //-------------------------
        //Get B-tagging scale fact
        //---------------------------
        double jetCorrPt = thisJet.Pt();
        double jetCorrE  = thisJet.E();
        if ( !isData )
        {
          //****************************************************************************
          //Apply b-tagging correction factor
          //****************************************************************************
          if ( !isData && abs(jetEta[i]) < 2.4 && jetCorrPt > BJET_CUT )
          {
            double effMedium = 0;
            double effLoose  = 0;
            BTagEntry::JetFlavor jetType = BTagEntry::FLAV_B;
            if ( abs(jetPartonFlavor[i]) == 5)
            {
              effMedium = btagMediumEfficiencyHist->GetBinContent( btagMediumEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0)),
									                                                 btagMediumEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
              effLoose  = btagLooseEfficiencyHist->GetBinContent( btagLooseEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0)),
									                                                btagLooseEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
              jetType   = BTagEntry::FLAV_B;
            }
            else if ( abs(jetPartonFlavor[i]) == 4)
            {
              effMedium = btagMediumCharmEfficiencyHist->GetBinContent( btagMediumCharmEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0)),
										                                                    btagMediumCharmEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
              effLoose  = btagLooseCharmEfficiencyHist->GetBinContent( btagLooseCharmEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0)),
									                                                     btagLooseCharmEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
              jetType   = BTagEntry::FLAV_C;
            }
            else
            {
              effMedium = btagMediumLightJetsEfficiencyHist->GetBinContent( btagMediumLightJetsEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0)),
										                                                        btagMediumLightJetsEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
		          effLoose  = btagLooseLightJetsEfficiencyHist->GetBinContent( btagLooseLightJetsEfficiencyHist->GetXaxis()->FindFixBin(fmax(fmin(jetCorrPt,199.9),10.0)),
										                                                       btagLooseLightJetsEfficiencyHist->GetYaxis()->FindFixBin(fabs(jetEta[i])));
              jetType  = BTagEntry::FLAV_UDSG;
            }
            //----------------
            //get scale factor
            //----------------
            double jetSF_Loo = -1;
            double jetSF_LooUp = -1;
            double jetSF_LooDown = -1;
            double jetSF_Med = -1;
            double jetSF_MedUp = -1;
            double jetSF_MedDown = -1;
            if ( abs(jetPartonFlavor[i]) == 5 || abs(jetPartonFlavor[i]) == 4 )//c,b quarks
            {
              if (jetCorrPt < 670.) //670 is the largest pt range listed in the CSV text file
              {
                //M
                jetSF_Med     = btagreaderM.eval(jetType, jetEta[i], jetCorrPt);
                jetSF_MedUp   = btagreaderM_up.eval(jetType, jetEta[i], jetCorrPt);
                jetSF_MedDown = btagreaderM_do.eval(jetType, jetEta[i], jetCorrPt);
                //L
                jetSF_Loo     = btagreaderL.eval(jetType, jetEta[i], jetCorrPt);
                jetSF_LooUp   = btagreaderL_up.eval(jetType, jetEta[i], jetCorrPt);
                jetSF_LooDown = btagreaderL_do.eval(jetType, jetEta[i], jetCorrPt);
              }
              else
              {
                //M
                jetSF_Med     = btagreaderM.eval(jetType, jetEta[i], 669);
                jetSF_MedUp   = btagreaderM_up.eval(jetType, jetEta[i], 669);
                jetSF_MedDown = btagreaderM_do.eval(jetType, jetEta[i], 669);
                //L
                jetSF_Loo     = btagreaderL.eval(jetType, jetEta[i], 669);
                jetSF_LooUp   = btagreaderL_up.eval(jetType, jetEta[i], 669);
                jetSF_LooDown = btagreaderL_do.eval(jetType, jetEta[i], 669);
              }
            }
            else//rest of the quarks
            {
              //M
              jetSF_Med     = btagreaderMistagM.eval(jetType, jetEta[i], 100);//fix in eta and pt
              jetSF_MedUp   = btagreaderMistagM_up.eval(jetType, jetEta[i], 100);
              jetSF_MedDown = btagreaderMistagM_do.eval(jetType, jetEta[i], 100);
              //L (to be checked)
              if ( jetCorrPt < 1000 )
              {
                jetSF_Loo     = btagreaderMistagL.eval(jetType, jetEta[i], jetCorrPt);
                jetSF_LooUp   = btagreaderMistagL_up.eval(jetType, jetEta[i], jetCorrPt);
                jetSF_LooDown = btagreaderMistagL_do.eval(jetType, jetEta[i], jetCorrPt);
              }
              else
              {
                jetSF_Loo     = btagreaderMistagL.eval(jetType, jetEta[i], 999);
                jetSF_LooUp   = btagreaderMistagL_up.eval(jetType, jetEta[i], 999);
                jetSF_LooDown = btagreaderMistagL_do.eval(jetType, jetEta[i], 999);
              }
            }
            //-----------------------
            //Loose Working point only
            //-----------------------
            //Apply Scale Factor
            if ( jetSF_Med <= 0 || jetSF_MedUp <= 0 || jetSF_MedDown <= 0  || jetSF_Loo <= 0 || jetSF_LooUp <= 0 || jetSF_LooDown <= 0 )
            {
              std::cout << "Warning: b-tag scale factor is <= 0!" << std::endl;
              std::cout << jetSF_Med << " " << jetSF_MedUp << " " << jetSF_MedDown << " " << jetSF_Loo << " " << jetSF_LooUp << " " << jetSF_LooDown << std::endl;
            }
            else if ( isCSVL(i,dataset) )
            {
              btagCorrFactor *= jetSF_Loo;
              if ( abs( jetPartonFlavor[i] ) == 5 || abs( jetPartonFlavor[i] ) == 4 )
              {
                sf_btagUp   *= jetSF_LooUp/jetSF_Loo;
                sf_btagDown *= jetSF_LooDown/jetSF_Loo;
              }
              else
              {
                sf_bmistagUp   *= jetSF_LooUp/jetSF_Loo;
                sf_bmistagDown *= jetSF_LooDown/jetSF_Loo;
              }
            }
            else
            {
              //only apply the scale factor on the inefficiency, if the corrected efficiency doesn't go above 100%
              //only record up/down systematics if the nominal and up and down corrected systematics do not go above 100%
              double sf = 1.0;
              if (effLoose*jetSF_Loo < 1.0) sf = (1/effLoose - jetSF_Loo) / (1/effLoose - 1);
              btagCorrFactor *= sf;
              if (abs(jetPartonFlavor[i]) == 5 || abs(jetPartonFlavor[i]) == 4)
              {
                if (effLoose*jetSF_Loo < 1.0 && effLoose*jetSF_LooUp < 1.0)
                {
                  sf_btagUp *= (1/effLoose - jetSF_LooUp) / (1/effLoose - 1) / sf;
                }
                if (effLoose*jetSF_Loo < 1.0 && effLoose*jetSF_LooDown < 1.0)
                {
                  sf_btagDown *= (1/effLoose - jetSF_LooDown) / (1/effLoose - 1) / sf;
                }
              }
              else
              {
                if ( effLoose*jetSF_Loo < 1.0 && effLoose*jetSF_LooUp < 1.0 )
                {
                  sf_bmistagUp *= (1/effLoose - jetSF_LooUp) / (1/effLoose - 1) / sf;
                }
                if ( effLoose*jetSF_Loo < 1.0 && effLoose*jetSF_LooDown < 1.0)
                {
                  sf_bmistagDown *= (1/effLoose - jetSF_LooDown) / (1/effLoose - 1) / sf;
                }
              }
            }
          }//Jetcut
        }//isData (Done with b-tagging scale factor)
        if( !isCSVL(i,dataset) ) continue;
        if (_jetdebug) std::cout << "jet no. " << i << " passed b-jet Loose WP ;\n"<< std::endl;
        GoodBJets.push_back(thisJet);
        GoodBJetsIsCVSL.push_back(isCSVL(i,dataset));
        GoodBJetsIsCVSM.push_back(isCSVM(i,dataset));
        GoodBJetsIsCVST.push_back(isCSVT(i,dataset));
        n_BJets++;
        //Filling bjet Candidate
        BjetCandidate tmp_bjetCand;
        tmp_bjetCand.Index  = i;
        tmp_bjetCand.bjet   = thisJet;
        tmp_bjetCand.isCSVL = isCSVL(i,dataset);
        tmp_bjetCand.isCSVM = isCSVM(i,dataset);
        tmp_bjetCand.isCSVT = isCSVT(i,dataset);
        bjetCand.push_back( tmp_bjetCand );
        nLooseBTaggedJets++;
        if( isCSVM(i,dataset) ) nMediumBTaggedJets++;
      }//end of loop over jet for b-jet
      std::cout << "[INFO]: end b jet loop \n"<< std::endl;

      //***************************************************
      //***************************************************
      //Create Event Categorization
      //***************************************************
      //***************************************************
      //-------------------------------
      //1) Look for Zmm Candidate
      //-------------------------------

      //Find two muons with the highest pt
      TLorentzVector ZCandidate(0,0,0,0);
      int ZMuIndex1 = -1;
      int ZMuIndex2 = -1;
      double bestDimuonPt = -1;
      std::vector< MuonCandidate > muSelectedCand;
      MuonCandidate bestMuCand[2];

      if( muCand.size() > 1 )
      {
        for ( size_t i = 0; i < muCand.size(); i++ )
        {
          for ( size_t j = i+1; j < muCand.size(); j++ )
          {
            MuonCandidate mu1 = muCand[i];
            MuonCandidate mu2 = muCand[j];
            //need dimuon mass between [76, 106] GeV
            double dimuonMass = (mu1.muon + mu2.muon).M();
            /*
            if( dimuonMass < 76 || dimuonMass > 106 )
            {
              if ( _debug ) std::cout << "[DEBUG]: Dimuon mass is out of range [76, 106]  GeV: dimuon masss-> " << dimuonMass << std::endl;
              if ( _debug ) std::cout << "... mu1Pt: " << mu1.muon.Pt()  << " mu2Pt: " << mu2.muon.Pt()  << std::endl;
              continue;
            }
            */
            //---------------------------------------------
            //if the sum of the muon pT's is larger than
            //that of the current Z candidate,
            //make this the Z candidate
            //---------------------------------------------
            if( mu1.muon.Pt() + mu2.muon.Pt() > bestDimuonPt )
            {
              bestDimuonPt = mu1.muon.Pt() + mu2.muon.Pt();
              ZCandidate = mu1.muon + mu2.muon;
              if ( mu1.muon.Pt() >= mu2.muon.Pt() )
              {
                if ( _debug ) std::cout << "assign muon candidate, mu1Pt > mu2Pt" << std::endl;
                bestMuCand[0] = mu1;
                bestMuCand[1] = mu2;
                ZMuIndex1 = mu1.Index;
                ZMuIndex2 = mu2.Index;
              }
              else
              {
                if ( _debug ) std::cout << "assign muon candidate, mu2Pt > mu1Pt" << std::endl;
                bestMuCand[0] = mu2;
                bestMuCand[1] = mu1;
                ZMuIndex1 = mu2.Index;
                ZMuIndex2 = mu1.Index;
              }
            }//best pt if
          }
        }
        //---------------------------------------
        //just use this container for convenience
        //to parse the data into TTree
        //---------------------------------------
        muSelectedCand.push_back(bestMuCand[0]);
        muSelectedCand.push_back(bestMuCand[1]);
        //if ( ZCandidate.M() >= 76. && ZCandidate.M() < 106 )
        //{
          //Fill in selected muon info
          razorbox = Zmm;
          lep1Type = 13 * -1 * bestMuCand[0].muonCharge;
          lep1Pt = bestMuCand[0].muon.Pt();
          lep1Eta = bestMuCand[0].muon.Eta();
          lep1Phi = bestMuCand[0].muon.Phi();
          lep1PassSelection = 1 + 2 * bestMuCand[0].isTightMuon;
          lep2Type = 13 * -1 * bestMuCand[1].muonCharge;
          lep2Pt = bestMuCand[1].muon.Pt();
          lep2Eta = bestMuCand[1].muon.Eta();
          lep2Phi = bestMuCand[1].muon.Phi();
          lep2PassSelection = 1 + 2 * bestMuCand[1].isTightMuon;
          //for MC apply lepton eff scale factor
          if (!isData )
          { 
            if ( matchesGenMuon(lep1Eta,lep1Phi)) {
                lep1EffSF *= helper->getLooseMuonScaleFactor( lep1Pt, lep1Eta, true);
                lep1EffSFUnc += helper->getLooseMuonScaleFactorError( lep1Pt, lep1Eta, true);
		leptonEffSF *=  lep1EffSF;
                lep1EffSFUp *= lep1EffSF + lep1EffSFUnc; 
                lep1EffSFDown *= lep1EffSF - lep1EffSFUnc; 
	    }
            if ( matchesGenMuon(lep2Eta,lep2Phi)) {
                lep2EffSF *= helper->getLooseMuonScaleFactor( lep2Pt, lep2Eta, true);
                lep2EffSFUnc += helper->getLooseMuonScaleFactorError( lep2Pt, lep2Eta, true);
		leptonEffSF *=  lep2EffSF;
                lep2EffSFUp *= lep2EffSF + lep2EffSFUnc; 
                lep2EffSFDown *= lep2EffSF - lep2EffSFUnc; 
            }
	    //max product of lep1/2sfup/down
	    leptonEffSFUp *= fmax(lep1EffSFUp, lep1EffSFDown)*fmax(lep2EffSFUp, lep2EffSFDown);
	    leptonEffSFDown *= fmin(lep1EffSFUp, lep1EffSFDown)*fmin(lep2EffSFUp, lep2EffSFDown);
	    leptonEffSFSys += leptonEffSF*sqrt( pow(lep1EffSFUnc/lep1EffSF, 2.) + pow(lep2EffSFUnc/lep2EffSF, 2.) );
          }
          //record Z candidate inf
          dileptonMass   = ZCandidate.M();
          bestDimuonPt   = ZCandidate.Pt();
          if ( _debug ) std::cout << "[DEBUG]: dimuon mass-> " << dileptonMass << " dimuon pT->" << bestDimuonPt << std::endl;
        //}
      }//end if muCand.size() > 1

      //-------------------------------
      //2) Look for Zee Candidate
      //-------------------------------

      if (razorbox == None)
      {
        //Find two electrons with the highest pt
        ZCandidate.SetPxPyPzE(0,0,0,0);
        int ZEleIndex1 = -1;
        int ZEleIndex2 = -1;
        double bestDielectronPt = -1;
        std::vector< ElectronCandidate > eleSelectedCand;
        ElectronCandidate bestEleCand[2];

        if( eleCand.size() > 1 )
        {
          for ( size_t i = 0; i < eleCand.size(); i++ )
          {
            for ( size_t j = i+1; j < eleCand.size(); j++ )
            {
              ElectronCandidate ele1 = eleCand[i];
              ElectronCandidate ele2 = eleCand[j];
              double dielectronMass = (ele1.electron + ele2.electron).M();
              /*
               * if( dielectronMass < 76 || dielectronMass > 106 )
              {
                if ( _debug ) std::cout << "[DEBUG]: Dielectron mass is out of range [76, 106]  GeV: dielectron masss-> " << dielectronMass << std::endl;
                if ( _debug ) std::cout << "... ele1Pt: " << ele1.electron.Pt()  << " ele2Pt: " << ele2.electron.Pt()  << std::endl;
                continue;
              }
              */
              //---------------------------------------------
              //if the sum of the electron pT's is larger than
              //that of the current Z candidate,
              //make this the Z candidate
              //---------------------------------------------
              if( ele1.electron.Pt() + ele2.electron.Pt() > bestDielectronPt )
              {
                bestDielectronPt = ele1.electron.Pt() + ele2.electron.Pt();
                ZCandidate = ele1.electron + ele2.electron;
                if ( ele1.electron.Pt() >= ele2.electron.Pt() )
                {
                  bestEleCand[0] = ele1;
                  bestEleCand[1] = ele2;
                  ZEleIndex1 = ele1.Index;
                  ZEleIndex2 = ele2.Index;
                }
                else
                {
                  if ( _debug ) std::cout << "assign electron candidate, ele2Pt > ele1Pt" << std::endl;
                  bestEleCand[0] = ele2;
                  bestEleCand[1] = ele1;
                  ZEleIndex1 = ele2.Index;
                  ZEleIndex2 = ele1.Index;
                }
              }//best pt if
            }
          }//end electron candidate loop
          //--------------------------------------
          //just use this container for convenienc
          //to parse the data into TTre
          //--------------------------------------
          eleSelectedCand.push_back(bestEleCand[0]);
          eleSelectedCand.push_back(bestEleCand[1]);

          //if ( ZCandidate.M() >= 76. && ZCandidate.M() < 106 )
          //{
            //Fill in selected electron info
            razorbox = Zee;
            lep1Type = 11 * -1 * bestEleCand[0].eleCharge;
            lep1Pt = bestEleCand[0].electron.Pt();
            lep1Eta = bestEleCand[0].electron.Eta();
            lep1Phi = bestEleCand[0].electron.Phi();
            lep1PassSelection = 1 + 2 * bestEleCand[0].isTightElectron;
            lep2Type = 11 * -1 * bestEleCand[1].eleCharge;
            lep2Pt = bestEleCand[1].electron.Pt();
            lep2Eta = bestEleCand[1].electron.Eta();
            lep2Phi = bestEleCand[1].electron.Phi();
            lep2PassSelection = 1 + 2 * bestEleCand[1].isTightElectron;
            //for MC apply lepton eff scale factor
	    if (!isData )
	    { 
            if ( matchesGenElectron(lep1Eta,lep1Phi)) {
                lep1EffSF *= helper->getLooseElectronScaleFactor( lep1Pt, lep1Eta, true);
                lep1EffSFUnc += helper->getLooseElectronScaleFactorError( lep1Pt, lep1Eta, true);
		leptonEffSF *=  lep1EffSF;
                lep1EffSFUp *= lep1EffSF + lep1EffSFUnc; 
                lep1EffSFDown *= lep1EffSF - lep1EffSFUnc; 
	    }
            if ( matchesGenElectron(lep2Eta,lep2Phi)) {
                lep2EffSF *= helper->getLooseElectronScaleFactor( lep2Pt, lep2Eta, true);
                lep2EffSFUnc += helper->getLooseElectronScaleFactorError( lep2Pt, lep2Eta, true);
		leptonEffSF *=  lep2EffSF;
                lep2EffSFUp *= lep2EffSF + lep2EffSFUnc; 
                lep2EffSFDown *= lep2EffSF - lep2EffSFUnc; 
            }
	    //max product of lep1/2sfup/down
	    leptonEffSFUp *= fmax(lep1EffSFUp, lep1EffSFDown)*fmax(lep2EffSFUp, lep2EffSFDown);
	    leptonEffSFDown *= fmin(lep1EffSFUp, lep1EffSFDown)*fmin(lep2EffSFUp, lep2EffSFDown);
	    leptonEffSFSys += leptonEffSF*sqrt( pow(lep1EffSFUnc/lep1EffSF, 2.) + pow(lep2EffSFUnc/lep2EffSF, 2.) );
	    }
            //record Z candidate info
            dileptonMass   = ZCandidate.M();
            bestDielectronPt   = ZCandidate.Pt();
            if ( _debug ) std::cout << "[DEBUG]: dielectron mass-> " << dileptonMass << " dielectron pT->" << bestDielectronPt << std::endl;
          //}//end check dielectron candidate
        }//end if eleCand.size()>1
      }//end checking if there was already a razor box requirement satisfied

      //-------------------------------
      //3) Look for Emu candidate
      //-------------------------------
      if (razorbox == None)
      {
        //Find two electrons with the highest pt
        ZCandidate.SetPxPyPzE(0,0,0,0);
        int EmuMuonIndex = -1;
        int EmuEleIndex = -1;
        double bestDileptonPt = -1;
        std::vector< ElectronCandidate > eleSelectedCandEmu;
        std::vector< MuonCandidate > muSelectedCandEmu;
        ElectronCandidate bestEmuEleCand;
        MuonCandidate bestEmuMuCand;

        if ( eleCand.size() > 0 && muCand.size() > 0)
        {
          for ( size_t i = 0; i < eleCand.size(); i++ )
          {
            for ( size_t j = 0; j < muCand.size(); j++ )
            {
              MuonCandidate mu      = muCand[j];
              ElectronCandidate ele = eleCand[i];
              double dileptonMass = (mu.muon + ele.electron).M();
              //---------------------------------------------
              //if the sum of the leptons pT's is larger than
              //that of the current Z candidate,
              //make this the Z candidate
              //---------------------------------------------
              if ( mu.muon.Pt() + ele.electron.Pt() > bestDileptonPt )
              {
                bestDileptonPt = mu.muon.Pt() + ele.electron.Pt();
                ZCandidate = mu.muon + ele.electron;
                if ( _debug ) std::cout << "assign electron and muon candidates" << std::endl;
                bestEmuMuCand = mu;
                bestEmuEleCand = ele;
                EmuMuonIndex = mu.Index;
                EmuEleIndex = ele.Index;
              } //best pt if
            } // finish muon loop
          } // finish electron loop
          //--------------------------------------
          //just use this container for convenienc
          //to parse the data into TTre
          //--------------------------------------
          muSelectedCandEmu.push_back(bestEmuMuCand);
          eleSelectedCandEmu.push_back(bestEmuEleCand);
          //Fill in selected lepton info
          if ( ZCandidate.M() > 0 )
          {
            razorbox = Emu;
            lep1Type = 13 * -1 * bestEmuMuCand.muonCharge;
            lep1Pt = bestEmuMuCand.muon.Pt();
            lep1Eta = bestEmuMuCand.muon.Eta();
            lep1Phi = bestEmuMuCand.muon.Phi();
            lep1PassSelection = 1 + 2 * bestEmuMuCand.isTightMuon;
            lep2Type = 11 * -1 * bestEmuEleCand.eleCharge;
            lep2Pt = bestEmuEleCand.electron.Pt();
            lep2Eta = bestEmuEleCand.electron.Eta();
            lep2Phi = bestEmuEleCand.electron.Phi();
            lep2PassSelection = 1 + 2 * bestEmuEleCand.isTightElectron;
            //for MC apply lepton eff scale factor
            /*if (!isData )
            {
              if ( matchesGenMuon(lep1Eta,lep1Phi)) leptonEffSF *=  helper->getVetoElectronScaleFactor( lep1Pt, lep1Eta, true);
              if ( matchesGenElectron(lep2Eta,lep2Phi)) leptonEffSF *=  helper->getLooseElectronScaleFactor( lep2Pt, lep2Eta, true);
            }*/
	    if (!isData )
	    { 
            if ( matchesGenMuon(lep1Eta,lep1Phi)) {
                lep1EffSF *= helper->getLooseMuonScaleFactor( lep1Pt, lep1Eta, true);
                lep1EffSFUnc += helper->getLooseMuonScaleFactorError( lep1Pt, lep1Eta, true);
		leptonEffSF *=  lep1EffSF;
                lep1EffSFUp *= lep1EffSF + lep1EffSFUnc; 
                lep1EffSFDown *= lep1EffSF - lep1EffSFUnc; 
	    }
            if ( matchesGenElectron(lep2Eta,lep2Phi)) {
                lep2EffSF *= helper->getLooseElectronScaleFactor( lep2Pt, lep2Eta, true);
                lep2EffSFUnc += helper->getLooseElectronScaleFactorError( lep2Pt, lep2Eta, true);
		leptonEffSF *=  lep2EffSF;
                lep2EffSFUp *= lep2EffSF + lep2EffSFUnc; 
                lep2EffSFDown *= lep2EffSF - lep2EffSFUnc; 
            }
	    //max product of lep1/2sfup/down
	    leptonEffSFUp *= fmax(lep1EffSFUp, lep1EffSFDown)*fmax(lep2EffSFUp, lep2EffSFDown);
	    leptonEffSFDown *= fmin(lep1EffSFUp, lep1EffSFDown)*fmin(lep2EffSFUp, lep2EffSFDown);
	    leptonEffSFSys += leptonEffSF*sqrt( pow(lep1EffSFUnc/lep1EffSF, 2.) + pow(lep2EffSFUnc/lep2EffSF, 2.) );
	    }
            //record Z candidate inf
            dileptonMass   = ZCandidate.M();
            bestDileptonPt = ZCandidate.Pt();
            if ( _debug ) std::cout << "[DEBUG]: dilepton mass-> " << dileptonMass << "dilepton pT->" << bestDileptonPt << std::endl;
          }//end checking that no other razorBox was satisfied.
        }// end if eleCand > 0 and muCand > 0
      }// end of Emu loop (check if razor box is empty)


      //------------------
      // One Muon Category
      //------------------
      TLorentzVector LeptonCandidate;
      if (razorbox == None)
      {
        double bestLeptonPt = -1;
        std::vector< MuonCandidate > muSelectedCandOneMu;
        MuonCandidate bestCandOneMu;
        if ( muCand.size() > 0 )
        {
          for( int i = 0; i < muCand.size(); i++ )
          {
            MuonCandidate mu = muCand[i];
            //-----------------------------------------
            //if the muon's pT is larger than that of
            //the current best Lepton pT, make it the
            //best lepton candidate
            //-----------------------------------------
            if (mu.muon.Pt() > bestLeptonPt)
            {
              bestLeptonPt = mu.muon.Pt();
              bestCandOneMu = mu;
            }
          } // end of muCand loop
          razorbox = OneMu;
          lep1Type = 13 * -1 * bestCandOneMu.muonCharge;
          lep1Pt = bestCandOneMu.muon.Pt();
          lep1Eta = bestCandOneMu.muon.Eta();
          lep1Phi = bestCandOneMu.muon.Phi();
          lep1PassSelection = 1 + 2 * bestCandOneMu.isTightMuon;
          LeptonCandidate.SetPtEtaPhiM( lep1Pt, lep1Eta, lep1Phi, 0.1057 );
          if (!isData )
          {
            if ( matchesGenMuon(lep1Eta,lep1Phi)) {
                lep1EffSF *= helper->getLooseMuonScaleFactor( lep1Pt, lep1Eta, true);
                lep1EffSFUnc += helper->getLooseMuonScaleFactorError( lep1Pt, lep1Eta, true);
		leptonEffSF *=  lep1EffSF;
                lep1EffSFUp *= lep1EffSF + lep1EffSFUnc; 
                lep1EffSFDown *= lep1EffSF - lep1EffSFUnc; 
		
		leptonEffSFUp *=  lep1EffSF + lep1EffSFUnc; 
		leptonEffSFDown *=  lep1EffSF - lep1EffSFUnc; 
                leptonEffSFSys += leptonEffSF*lep1EffSFUnc;
	    }
          }
        } // end if muCand.size() > 0 loop
      }//end checking if another razorBox was already there

      //---------------------
      //One Electron Category
      //---------------------
      if (razorbox == None)
      {
        double bestLeptonPt = -1;
        std::vector< ElectronCandidate > eleSelectedCandOneEle;
        ElectronCandidate bestCandOneEle;

        if ( eleCand.size() > 0 )
        {
          for( int i = 0; i < eleCand.size(); i++ )
          {
            ElectronCandidate ele = eleCand[i];
		  		  if (ele.electron.Pt() > bestLeptonPt)
            {
              bestLeptonPt = ele.electron.Pt();
              bestCandOneEle = ele;
            }
          } //end of eleCand loop
          razorbox = OneEle;
          lep1Type = 11 * -1 * bestCandOneEle.eleCharge;
          lep1Pt = bestCandOneEle.electron.Pt();
          lep1Eta = bestCandOneEle.electron.Eta();
          lep1Phi = bestCandOneEle.electron.Phi();
          lep1PassSelection = 1 + 2 * bestCandOneEle.isTightElectron;
          LeptonCandidate.SetPtEtaPhiM( lep1Pt, lep1Eta, lep1Phi, 0.000511 );

          if (!isData )
          {
            if ( matchesGenElectron(lep1Eta,lep1Phi)) {
                lep1EffSF *= helper->getLooseElectronScaleFactor( lep1Pt, lep1Eta, true);
                lep1EffSFUnc += helper->getLooseElectronScaleFactorError( lep1Pt, lep1Eta, true);
		leptonEffSF *=  lep1EffSF;
                lep1EffSFUp *= lep1EffSF + lep1EffSFUnc; 
                lep1EffSFDown *= lep1EffSF - lep1EffSFUnc; 
		
		leptonEffSFUp *=  lep1EffSF + lep1EffSFUnc; 
		leptonEffSFDown *=  lep1EffSF - lep1EffSFUnc; 
                leptonEffSFSys += leptonEffSF*lep1EffSFUnc;
            }
          }
        } // end of if eleCand.size() > 0 loop
      }//end of one lepton category (check that razor box was not yet assigned)

      //----------------
      //High-pt category
      //----------------
      //if( razorbox == None && HiggsCandidate.Pt() > 110. ) razorbox = HighPt;

      //------------
      //Hbb category
      //------------
      TLorentzVector HbbZbbCandidate(0,0,0,0);
      std::vector< BjetCandidate > bjetSelectedCand;
      BjetCandidate bestBjetCandHbb[2];
      if( razorbox == None )
      {
        //if there are two or more loose b-tags and one medium b-tag, look for b-bbar resonances
        //if( nLooseBTaggedJets > 1 && nMediumBTaggedJets > 0 )
        if( nLooseBTaggedJets > 1 )
        {
          for(int i = 0; i < int(bjetCand.size()); i++)
          {
            for(int j = i+1; j < int(bjetCand.size()); j++)
            {
              BjetCandidate bjet1 = bjetCand[i];
              BjetCandidate bjet2 = bjetCand[j];
              //if neither of the b-jets passes CSVL, continue
              //if( !bjet1.isCSVL && !bjet2.isCSVL ) continue;
              double mbb = (bjet1.bjet + bjet2.bjet).M();
              double pTbb = (bjet1.bjet + bjet2.bjet).Pt();
              //if mbb is closer to the higgs mass than mbbH, make mbbH = mbb
              if( fabs(mbbH - 125.0) > fabs(mbb - 125.0) )
              {
                mbbH = mbb;
                pTbbH = pTbb;
                bestBjetCandHbb[0] = bjet1;
                bestBjetCandHbb[1] = bjet2;
                HbbZbbCandidate = bjet1.bjet + bjet2.bjet;
              }
            }
          }

          if ( HbbZbbCandidate.M() > 95. && HbbZbbCandidate.M() < 140. )
          {
            razorbox = Hbb;
            mbbH     = HbbZbbCandidate.M();
            pTbbH    = HbbZbbCandidate.Pt();
            bjetSelectedCand.push_back(bestBjetCandHbb[0]);
            bjetSelectedCand.push_back(bestBjetCandHbb[1]);
          }
        }
      }//end Hbb category

      //------------
      //Zbb category
      //------------
      HbbZbbCandidate.SetPxPyPzE(0,0,0,0);
      BjetCandidate bestBjetCandZbb[2];
      if( razorbox == None )
      {
        //make sure container is empty before start.
        bjetSelectedCand.clear();
        //if there are two or more loose b-tags and one medium b-tag, look for b-bbar resonances
        //if( nLooseBTaggedJets > 1 && nMediumBTaggedJets > 0 )
        if( nLooseBTaggedJets > 1 )
        {
          for(int i = 0; i < int(bjetCand.size()); i++)
          {
            for(int j = i+1; j < int(bjetCand.size()); j++)
            {
              BjetCandidate bjet1 = bjetCand[i];
              BjetCandidate bjet2 = bjetCand[j];
              //if neither of the b-jets passes CSVL, continue
              //if( !bjet1.isCSVL && !bjet2.isCSVL ) continue;
              double mbb = (bjet1.bjet + bjet2.bjet).M();
              double pTbb = (bjet1.bjet + bjet2.bjet).Pt();
              //if mbb is closer to the higgs mass than mbbH, make mbbH = mbb
              if( fabs(mbbZ - 91.2) > fabs(mbb - 91.2) )
              {
                mbbZ = mbb;
                pTbbZ = pTbb;
                bestBjetCandZbb[0] = bjet1;
                bestBjetCandZbb[1] = bjet2;
                HbbZbbCandidate = bjet1.bjet + bjet2.bjet;
              }
            }
          }
          if ( HbbZbbCandidate.M() > 60. && HbbZbbCandidate.M() < 95. )
          {
            razorbox = Zbb;
            mbbZ     = HbbZbbCandidate.M();
            pTbbZ    = HbbZbbCandidate.Pt();
            bjetSelectedCand.push_back(bestBjetCandZbb[0]);
            bjetSelectedCand.push_back(bestBjetCandZbb[1]);
          }
        }
      }//end Zbb category

      //----------------
      //High-pt category
      //----------------
      if( razorbox == None && HiggsCandidate.Pt() > 110. ) razorbox = HighPt;
      
      //------------------------------------------------
      //I n v a ri a n t   m a s s   r e s o l u t i o n
      //------------------------------------------------
      sigmaMoverM = 0.5*sqrt( Pho_sigmaEOverE[0]*Pho_sigmaEOverE[0] + Pho_sigmaEOverE[1]*Pho_sigmaEOverE[1] );
      //inclusive HggRazor
      if(razorbox == None)
      {
        //HighRes Box
        if( sigmaMoverM < 0.0085 ) razorbox = HighRes;
        //LowRes Box
        else razorbox = LowRes;
      }
      if (_debug) cout << "razorbox = : " << razorbox << "\n";

      //------------------------------------------------------------------------------------------
      //Jets (need to go last since we need to know the bjets forming the Hbb or Zbb candidates)
      //------------------------------------------------------------------------------------------
      //Propagate jet uncertainties to MET
      float MetXCorr_JESUp = 0;
      float MetYCorr_JESUp = 0;
      float MetXCorr_JESDown = 0;
      float MetYCorr_JESDown = 0;
      //variables for 2017F MET RECIPE
      float metX_RecV2 =0;
      float metY_RecV2 =0;

      vector<TLorentzVector> GoodJets;
      std::vector< JetCandidate > jetCand;
      //vector<bool> GoodJetsIsCVSL;
      //vector<bool> GoodJetsIsCVSM;
      //vector<bool> GoodJetsIsCVST;
      vector<TLorentzVector> GoodJetsJESUp;
      vector<TLorentzVector> GoodJetsJESDown;

      std::cout << "[INFO]: begin jet loop \n"<< std::endl;
      for(int i = 0; i < nJets; i++)
      {
        //Jet Energy Corrections
        double JEC = JetEnergyCorrectionFactor( jetPt[i], jetEta[i], jetPhi[i], jetE[i],
						                                    fixedGridRhoAll, jetJetArea[i], runNum,
						                                    JetCorrectorIOV, JetCorrector );
        TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );
        if (_jetdebug) std::cout << "jet no. " << i << " raw info (pt, eta, phi, E):" << jetPt[i] << " , " << jetEta[i] << " , " << jetPhi[i] << " , " << jetE[i] << " ;\n"<< std::endl;
        if (_jetdebug) std::cout << "jet no. " << i << " jec info (pt, eta, phi, E):" << jetPt[i]*JEC << " , " << jetEta[i] << " , " << jetPhi[i] << " , " << jetE[i]*JEC << " ;\n"<< std::endl;
        //2017F MET RECIPE
        TLorentzVector thisJet_Raw = makeTLorentzVector( jetPt[i], jetEta[i], jetPhi[i], jetE[i] );
        //if(thisJet.Pt() < 75. &&  (fabs(thisJet.Eta()) > 2.65 && fabs(thisJet.Eta()) < 3.139) )
        if(thisJet.Pt() < 50. &&  (fabs(thisJet.Eta()) > 2.65 && fabs(thisJet.Eta()) < 3.139) )
        {
                metX_RecV2 += thisJet.Px();
                metY_RecV2 += thisJet.Py();
                if(_metdebug) std::cout << "[satisfy MET_V2]: no. " << i << "\n"<< std::endl;
                if(_metdebug) std::cout << "[satisfy MET_V2]: jetPt " << jetPt[i] << ", jetPt*JEC " << jetPt[i]*JEC << "\n"<< std::endl;
                if(_metdebug) std::cout << "[satisfy MET_V2]: jetJECPx " << thisJet.Px() << ", jetJECPy " << thisJet.Py() << "\n"<< std::endl;
        }
        if(_metdebug) std::cout << "[MET_V2]: metX_RecV2 " << metX_RecV2 << "\n"<< std::endl;
        if(_metdebug) std::cout << "[MET_V2]: metY_RecV2 " << metY_RecV2 << "\n"<< std::endl;
        //ECAL Prefiring Recipe
        if( thisJet.Pt() > 100. && fabs( thisJet.Eta() ) > 2.25 && fabs ( thisJet.Eta() ) < 3.0 ) ecal_prefiring_affected = true;
        //if( thisJet.Pt() > 100. && fabs( thisJet.Eta() ) > 2.25 && fabs ( thisJet.Eta() ) < 3.0 ) continue;
        if(_ecaldebug) std::cout << "[ECAL Prefiring: not affected \n"<< std::endl;

        if( thisJet.Pt() < JET_CUT ) continue;//According to the April 1st 2015 AN
        if (_jetdebug) std::cout << "jet no. " << i << " passed jet pt 30 cut ;\n"<< std::endl;
        if( fabs( thisJet.Eta() ) >= 2.4 ) continue;
        if (_jetdebug) std::cout << "jet no. " << i << " passed jet eta 2.4 cut ;\n"<< std::endl;
        if (!isFastsimSMS)
        {
          if ( !jetPassIDLoose[i] ) continue;
        }

        bool overlapjm = false;
        for(int j = 0; j < int(GoodMuons.size()); j++)
        {
          TLorentzVector mu = GoodMuons.at(j);
          if (RazorAnalyzer::deltaR( thisJet.Eta(), thisJet.Phi(), mu.Eta(), mu.Phi()) < 0.4 ) overlapjm = true;
        }
        if(overlapjm) continue;
        if (_jetdebug) std::cout << "jet no. " << i << " passed jet muon cleaning ;\n"<< std::endl;

        bool overlapje = false;
        for(int k = 0; k < int(GoodElectrons.size()); k++)
        {
          TLorentzVector ele = GoodElectrons.at(k);
          if (RazorAnalyzer::deltaR( thisJet.Eta(), thisJet.Phi(), ele.Eta(), ele.Phi()) < 0.4 ) overlapje = true;
        }
        if(overlapje) continue;
        if (_jetdebug) std::cout << "jet no. " << i << " passed jet ele cleaning ;\n"<< std::endl;
        //exclude selected photons from the jet collection
        double deltaRJetPhoton = min( thisJet.DeltaR( pho_cand_vec[0] ), thisJet.DeltaR( pho_cand_vec[1] ) );
        if ( deltaRJetPhoton <= 0.4 ) continue;
        if (_jetdebug) std::cout << "jet no. " << i << " passed jet photon cleaning ;\n"<< std::endl;
        //Exclude selected b-jets from the jet collection
        bool overlapbjj = false;
        for(int b = 0; b < int(bjetSelectedCand.size()); b++)
        {
          BjetCandidate Bjet = bjetSelectedCand[b];
          if (RazorAnalyzer::deltaR( thisJet.Eta(), thisJet.Phi(), Bjet.bjet.Eta(), Bjet.bjet.Phi()) < 0.4 ) overlapbjj = true;
        }
        if(overlapbjj) continue;
        if (_jetdebug) std::cout << "jet no. " << i << " passed jet b-jet cand cleaning ;\n"<< std::endl;
        JetCandidate  tmp_jetCand;
        tmp_jetCand.jet    = thisJet;
        tmp_jetCand.isCSVL = isCSVL(i,dataset);
        tmp_jetCand.isCSVM = isCSVM(i,dataset);
        tmp_jetCand.isCSVT = isCSVT(i,dataset);
        jetCand.push_back( tmp_jetCand );
        GoodJets.push_back(thisJet);
        n_Jets++;

        double jetCorrPt = thisJet.Pt();
        double jetCorrE  = thisJet.E();
        if ( !isData )
        {
          double unc = helper->getJecUnc( jetCorrPt, jetEta[i], 999 ); //use run=999 by default
          double jetPtJESUp = jetCorrPt*(1+unc);
          double jetPtJESDown = jetCorrPt/(1+unc);
          double jetEJESUp = jetCorrE*(1+unc);
          double jetEJESDown = jetCorrE/(1+unc);
          TLorentzVector thisJetJESUp = makeTLorentzVector(jetPtJESUp, jetEta[i], jetPhi[i], jetEJESUp);
          TLorentzVector thisJetJESDown = makeTLorentzVector(jetPtJESDown, jetEta[i], jetPhi[i], jetEJESDown);
          //Propagate uncertainties to the MET
          if (jetPtJESUp > 20)
          { 
                          MetXCorr_JESUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
                          MetYCorr_JESUp += -1 * (thisJetJESUp.Py() - thisJet.Py());
          }
          if (jetPtJESDown > 20)
          {
                          MetXCorr_JESDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
                          MetYCorr_JESDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
          }
          if ( jetPtJESUp > JET_CUT )
          {
            GoodJetsJESUp.push_back(thisJetJESUp);
            n_Jets_JESUp++;
          }
          if ( jetPtJESDown > JET_CUT )
          {
            GoodJetsJESDown.push_back(thisJetJESDown);
            n_Jets_JESDown++;
          }
        }
      } //loop over jets
      std::cout << "[INFO]: end jet loop \n"<< std::endl;
      //--------------------
      //Store jet info
      //---------------------
      for ( int iJet = 0; iJet < int(jetCand.size()) ; iJet++ )
      {
        jet_E[iJet]     = jetCand[iJet].jet.E();
        jet_Pt[iJet]    = jetCand[iJet].jet.Pt();
        jet_Eta[iJet]   = jetCand[iJet].jet.Eta();
        jet_Phi[iJet]   = jetCand[iJet].jet.Phi();
        jetIsCSVL[iJet] = jetCand[iJet].isCSVL;
        jetIsCSVM[iJet] = jetCand[iJet].isCSVM;
        jetIsCSVT[iJet] = jetCand[iJet].isCSVT;
      }

      //----------------------------------------------------------
      //Add Hbb/Zbb candiates to jet collection stored in the Tree
      //----------------------------------------------------------
      for ( auto bjet : bjetSelectedCand )
      {
        jet_E[n_Jets]     = bjet.bjet.E();
        jet_Pt[n_Jets]    = bjet.bjet.Pt();
        jet_Eta[n_Jets]   = bjet.bjet.Eta();
        jet_Phi[n_Jets]   = bjet.bjet.Phi();
        jetIsCSVL[n_Jets] = bjet.isCSVL;
        jetIsCSVM[n_Jets] = bjet.isCSVM;
        jetIsCSVT[n_Jets] = bjet.isCSVT;
        n_Jets++;
      }

      //Compute the razor variables using the selected jets and the diphoton system
      HT = Pho_Pt[0] + Pho_Pt[1]; //HT = sum of photon pT  + jet pT
      vector<TLorentzVector> ObjectCandidates;

      //Add all jet but Hbb or Zbb candidate jets
      for( auto& jet : GoodJets )
      {
        ObjectCandidates.push_back(jet);
        HT += jet.Pt();
      }

      //add H->gg candidate
      ObjectCandidates.push_back(HiggsCandidate);
      //Add Z candidate
      if ( razorbox == Zmm || razorbox == Zee || razorbox == Emu ) ObjectCandidates.push_back(ZCandidate);
      //Add leptons
      if ( razorbox == OneMu || razorbox == OneEle ) ObjectCandidates.push_back(LeptonCandidate);
      //Add Hbb or Zbb Candidate
      if ( razorbox == Hbb || razorbox == Zbb ) ObjectCandidates.push_back(HbbZbbCandidate);

      TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
      TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM( metType1Pt, 0, metType1Phi, 0 );
      TLorentzVector genMET = makeTLorentzVectorPtEtaPhiM( genMetPt, 0, genMetPhi, 0 );

      MET = metPt;
      t1MET = metType1Pt;
      if(_metdebug) std::cout << "[t1MET]: metType1Pt " << metType1Pt << "\n"<< std::endl;
      if(_metdebug) std::cout << "[t1MET]: t1MET " << t1MET << "\n"<< std::endl;
      if(_metdebug) std::cout << "[t1MET]: t1PFMET.Px() " << t1PFMET.Px() << "\n"<< std::endl;
      if(_metdebug) std::cout << "[t1MET]: t1PFMET.Py() " << t1PFMET.Py() << "\n"<< std::endl;

        //2017F MET RECIPE
        float MetX_RecV2 = PFMET.Px() + metX_RecV2;
        float MetY_RecV2 = PFMET.Py() + metY_RecV2;
        float t1MetX_RecV2 = t1PFMET.Px() + metX_RecV2;
        float t1MetY_RecV2 = t1PFMET.Py() + metY_RecV2;
        //if(_metdebug) std::cout << "[MET_V2]: MetX_RecV2 " << MetX_RecV2 << "\n"<< std::endl;
        //if(_metdebug) std::cout << "[MET_V2]: MetY_RecV2 " << MetY_RecV2 << "\n"<< std::endl;
        if(_metdebug) std::cout << "[t1MET_V2]: t1MetX_RecV2 " << t1MetX_RecV2 << "\n"<< std::endl;
        if(_metdebug) std::cout << "[t1MET_V2]: t1MetY_RecV2 " << t1MetY_RecV2 << "\n"<< std::endl;

        TLorentzVector PFMET_RecV2(MetX_RecV2, MetY_RecV2, 0, sqrt(pow(MetX_RecV2,2) + pow(MetY_RecV2,2) ));
        TLorentzVector t1PFMET_RecV2(t1MetX_RecV2, t1MetY_RecV2, 0, sqrt(pow(t1MetX_RecV2,2) + pow(t1MetY_RecV2,2) ));

        MET_RecV2  = PFMET_RecV2.Pt();
        t1MET_RecV2  = t1PFMET_RecV2.Pt();
        Myt1MET = t1MET_RecV2;
        if(_metdebug) std::cout << "[Myt1MET]: Myt1MET " << Myt1MET << "\n"<< std::endl;
        if(_metdebug) std::cout << "[t1MET_V2]: t1MET_RecV2 " << t1MET_RecV2 << "\n"<< std::endl;
        if(_metdebug) std::cout << "[t1MET_V2]: t1PFMET_RecV2.Px() " << t1PFMET_RecV2.Px() << "\n"<< std::endl;
        if(_metdebug) std::cout << "[t1MET_V2]: t1PFMET_RecV2.Py() " << t1PFMET_RecV2.Py() << "\n"<< std::endl;

      //need to implement the custom type1 MET corrections
      //Note: need to compute lep1MT at this point
      if (razorbox == OneMu)
      {
        lep1MT = sqrt(0.1057*0.1057 + 2*t1PFMET.Pt()*lep1Pt*(1 - cos(deltaPhi(t1PFMET.Phi(),lep1Phi))));
        lep1GenMetMT = sqrt(0.1057*0.1057 + 2*genMET.Pt()*lep1Pt*(1 - cos(deltaPhi(genMET.Phi(),lep1Phi))));
      }
      else if (razorbox == OneEle)
      {
        lep1MT = sqrt(0.000511*0.000511 + 2*t1PFMET.Pt()*lep1Pt*(1 - cos(deltaPhi(t1PFMET.Phi(),lep1Phi))));
        lep1GenMetMT = sqrt(0.000511*0.000511 + 2*genMET.Pt()*lep1Pt*(1 - cos(deltaPhi(genMET.Phi(),lep1Phi))));
      }

      vector<TLorentzVector> hemispheres = getHemispheres(ObjectCandidates);
      theMR  = computeMR(hemispheres[0], hemispheres[1]);
      if ( theMR > 0 )
      {
        theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
        //t1Rsq  = computeRsq(hemispheres[0], hemispheres[1], t1PFMET);
        t1Rsq  = computeRsq(hemispheres[0], hemispheres[1], t1PFMET_RecV2);
        genMetRsq = computeRsq(hemispheres[0], hemispheres[1], genMET);
      }

      //***********************
      //MR skim
      //***********************
      if (doMRSkim)
      {
        if (!(theMR > 150 && GoodJets.size() >= 1)) continue;
      }

      //Calculations for JES systematics
      if( !isData )
      {
        //JES up
        vector<TLorentzVector> ObjectCandidates_JESUp;
        for( auto& jet : GoodJetsJESUp ) ObjectCandidates_JESUp.push_back(jet);
        ObjectCandidates_JESUp.push_back(HiggsCandidate);

        float PFMetXJESUp   = PFMET.Px() + MetXCorr_JESUp;
        float PFMetYJESUp   = PFMET.Py() + MetYCorr_JESUp;
        float t1PFMetXJESUp = t1PFMET.Px() + MetXCorr_JESUp;
        float t1PFMetYJESUp = t1PFMET.Py() + MetYCorr_JESUp;

        TLorentzVector PFMET_JESUp(PFMetXJESUp, PFMetYJESUp, 0, sqrt( pow(PFMetXJESUp,2) + pow(PFMetYJESUp,2) ));
        TLorentzVector t1PFMET_JESUp(t1PFMetXJESUp, t1PFMetYJESUp, 0, sqrt( pow(t1PFMetXJESUp,2) + pow(t1PFMetYJESUp,2) ));
        vector<TLorentzVector> hemispheres_JESUp = getHemispheres(ObjectCandidates_JESUp);
        theMR_JESUp  = computeMR(hemispheres_JESUp[0], hemispheres_JESUp[1]);
        theRsq_JESUp = computeRsq(hemispheres_JESUp[0], hemispheres_JESUp[1], PFMET_JESUp);
        t1Rsq_JESUp  = computeRsq(hemispheres_JESUp[0], hemispheres_JESUp[1], t1PFMET_JESUp);
        MET_JESUp    = PFMET_JESUp.Pt();
        t1MET_JESUp  = t1PFMET_JESUp.Pt();

	      //JES down
	      vector<TLorentzVector> ObjectCandidates_JESDown;
	      for( auto& jet : GoodJetsJESDown ) ObjectCandidates_JESDown.push_back(jet);
	      ObjectCandidates_JESDown.push_back(HiggsCandidate);

        float PFMetXJESDown   = PFMET.Px() + MetXCorr_JESDown;
        float PFMetYJESDown   = PFMET.Py() + MetYCorr_JESDown;
        float t1PFMetXJESDown = t1PFMET.Px() + MetXCorr_JESDown;
        float t1PFMetYJESDown = t1PFMET.Py() + MetYCorr_JESDown;

        TLorentzVector PFMET_JESDown(PFMetXJESDown, PFMetYJESDown, 0, sqrt( pow(PFMetXJESDown,2) + pow(PFMetYJESDown,2) ));
        TLorentzVector t1PFMET_JESDown(t1PFMetXJESDown, t1PFMetYJESDown, 0, sqrt( pow(t1PFMetXJESDown,2) + pow(t1PFMetYJESDown,2) ));
        vector<TLorentzVector> hemispheres_JESDown = getHemispheres(ObjectCandidates_JESDown);
        theMR_JESDown  = computeMR(hemispheres_JESDown[0], hemispheres_JESDown[1]);
        theRsq_JESDown = computeRsq(hemispheres_JESDown[0], hemispheres_JESDown[1], PFMET_JESDown);
        t1Rsq_JESDown  = computeRsq(hemispheres_JESDown[0], hemispheres_JESDown[1], t1PFMET_JESDown);
        MET_JESDown    = PFMET_JESDown.Pt();
        t1MET_JESDown  = t1PFMET_JESDown.Pt();
        
      }

      if ( theMR < 0.0 )
      {
        if ( _debug ) std::cout << "[INFO]: MR < 150 GeV, MR: " << theMR << std::endl;
        for ( auto& jet : ObjectCandidates )
        {
          if ( _debug ) std::cout << "phoPT: " << pTGammaGamma
				      << " jet pt : " << jet.Pt() << " eta: " << jet.Eta() << " phi: " << jet.Phi()
				      << " h1 pt: " << hemispheres[0].Pt() << " h1 eta: " << hemispheres[0].Eta()
				      << " h2 pt: " << hemispheres[1].Pt() << " h2 eta: " << hemispheres[1].Eta() << std::endl;
        }
      }

      //Fill Event
        //if ( _debug ) std::cout << "run == " << run << " && evt == " << event << " && ecal_prefiring_affected == " << ecal_prefiring_affected << std::endl;
	//if(ecal_prefiring_affected) continue;
      if (_debug) cout << "Fill event: " << mChi << " " << theMR << " " << t1Rsq << " " << sigmaMoverM << "\n";

      if (!isFastsimSMS)
      {
        razorTree->Fill();
      }
      else if (parsedLHE)
      {
        if (!is2DMassScan)
        {
          smsTrees[mChi]->Fill();
        }
        else
        {
          pair<int,int> smsPair = make_pair(mChi, mLSP);
          smsTrees2D[smsPair]->Fill();
        }
      }

    }//end of loop

  if ( _info ) std::cout << "[INFO]: Number of events processed: " << NEvents->Integral() << std::endl;

  if(!isFastsimSMS)
  {
    if ( _info ) std::cout << "[INFO]: Writing output trees..." << std::endl;
    outFile->cd();
    razorTree->Write();
    NEvents->Write();
    SumWeights->Write();
    SumScaleWeights->Write();
    SumPdfWeights->Write();
    histNISRJets->Write();
    histPtISR->Write();
    histNPV->Write();
    puhisto->Write();
  }
  else
  {
    if (!is2DMassScan)
    {
      for(auto &filePtr : smsFiles)
      {
        cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
        filePtr.second->cd();
        smsTrees[filePtr.first]->Write();
        smsNEvents[filePtr.first]->Write("NEvents");
        smsSumWeights[filePtr.first]->Write("SumWeights");
        smsSumScaleWeights[filePtr.first]->Write("SumScaleWeights");
        smsSumPdfWeights[filePtr.first]->Write("SumPdfWeights");
        smsNISRJets[filePtr.first]->Write("NISRJets");
        smsPtISR[filePtr.first]->Write("PtISR");
        smsNPV[filePtr.first]->Write("NPV");
      }
    }
    else
    {
      for(auto &filePtr : smsFiles2D)
      {
        cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
        filePtr.second->cd();
        smsTrees2D[filePtr.first]->Write();
        smsNEvents2D[filePtr.first]->Write("NEvents");
        smsSumWeights2D[filePtr.first]->Write("SumWeights");
        smsSumScaleWeights2D[filePtr.first]->Write("SumScaleWeights");
        smsSumPdfWeights2D[filePtr.first]->Write("SumPdfWeights");
        smsNISRJets2D[filePtr.first]->Write("NISRJets");
        smsPtISR2D[filePtr.first]->Write("PtISR");
        smsNPV2D[filePtr.first]->Write("NPV");
      }
    }
  }

  outFile->Close();
  delete photonCorrector;
  delete photonCorrector_2017;
  delete helper;

}
