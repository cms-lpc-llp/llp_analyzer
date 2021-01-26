#include "RazorHelper.h"

// Constructor
RazorHelper::RazorHelper(std::string tag_, bool isData_, bool isFastsim_):
        tag(tag_), isData(isData_), isFastsim(isFastsim_) {
    std::cout << "RazorHelper initializing with tag " << tag << std::endl;

    eleVetoEffSFMinPt = -1;
    eleLooseEffSFMinPt = -1;
    muVetoEffSFMinPt = -1;
    muLooseEffSFMinPt = -1;

    // check that CMSSW is set up
    loadCMSSWPath();
    if (cmsswPath == "") {
        loadTag_Null();
        return;
    }

    // tag for 2015 data
    else if (tag == "Razor2015") {
        loadTag_Razor2015();
    }

    // tag for 2015 76X ReReco data
    else if (tag == "Razor2015_76X") {
        loadTag_Razor2015_76X();
    }

    // tag for 2016 17Aug2017 Rereco
    else if (tag == "Razor2016_07Aug2017Rereco"){
        loadTag_Razor2016_07Aug2017Rereco();
    }

    // tag for 2016 80X Moriond Rereco
    else if (tag == "Razor2016_MoriondRereco") {
        loadTag_Razor2016_MoriondRereco();
    }

    // tag for 2016 80X 03Feb2017 Rereco
    else if (tag == "Razor2016_03Feb2017Rereco") {
        loadTag_Razor2016_03Feb2017Rereco();
    }

    // tag for 2016G 80X data
    else if (tag == "Razor2016G_80X") {
        loadTag_Razor2016G_80X();
    }

    // tag for 2016G unblinded 80X data
    else if (tag == "Razor2016G_SUSYUnblind_80X") {
      //eleVetoEffSFMinPt = 15.01;
      //muVetoEffSFMinPt = 15.01;
      loadTag_Razor2016G_SUSYUnblind_80X();
    }

    // tag for 2016 ICHEP 80X data
    else if (tag == "Razor2016_ICHEP_80X") {
        loadTag_Razor2016_ICHEP_80X();
    }

    // tag for 2017 Prompt Reco
    else if (tag == "Razor2017_92X") {
        loadTag_Razor2017_92X();
    }

    // tag for 2017 17Nov2017 Rereco
    else if (tag == "Razor2017_17Nov2017Rereco") {
        loadTag_Razor2017_17Nov2017Rereco();
    }

    // tag for 2017 31Mar2018 Rereco
    else if (tag == "Razor2017_31Mar2018Rereco") {
        loadTag_Razor2017_31Mar2018Rereco();
    }
    else if (tag == "Razor2018_17SeptEarlyReReco"){
      loadTag_Razor2018_17SeptEarlyReReco();
    }
    else if(tag == "Razor2016_Source2018")
    {
      loadTag_Razor2016_Source2018();
    }
    else if(tag == "Razor2017_Source2018")
    {
      loadTag_Razor2017_Source2018();
    }
   // tag not found
    else {
        std::cout << "Error in RazorHelper::RazorHelper : specified tag " << tag << " is not supported!" << std::endl;
        loadTag_Null();
        return;
    }

}

// Destructor -- close all of the open TFile objects
RazorHelper::~RazorHelper() {
    if (pileupWeightFile) {
        pileupWeightFile->Close();
        delete pileupWeightFile;
    }
    if (eleTightEfficiencyFile) {
        eleTightEfficiencyFile->Close();
        delete eleTightEfficiencyFile;
    }
    if (eleVetoEfficiencyFile) {
        eleVetoEfficiencyFile->Close();
        delete eleVetoEfficiencyFile;
    }
    if (eleLooseEfficiencyFile) {
        eleLooseEfficiencyFile->Close();
        delete eleLooseEfficiencyFile;
    }
    if (eleEffSFFile) {
        eleEffSFFile->Close();
        delete eleEffSFFile;
    }
    if (vetoEleEffSFFile) {
        vetoEleEffSFFile->Close();
        delete vetoEleEffSFFile;
    }
    if (muTightEfficiencyFile) {
        muTightEfficiencyFile->Close();
        delete muTightEfficiencyFile;
    }
    if (muVetoEfficiencyFile) {
        muVetoEfficiencyFile->Close();
        delete muVetoEfficiencyFile;
    }
    if (muLooseEfficiencyFile) {
        muLooseEfficiencyFile->Close();
        delete muLooseEfficiencyFile;
    }
    if (muEffSFFile) {
        muEffSFFile->Close();
        delete muEffSFFile;
    }
    if (vetoMuEffSFFile) {
        vetoMuEffSFFile->Close();
        delete vetoMuEffSFFile;
    }
    if (tauEfficiencyFile) {
        tauEfficiencyFile->Close();
        delete tauEfficiencyFile;
    }
    if (btagEfficiencyFile) {
        btagEfficiencyFile->Close();
        delete btagEfficiencyFile;
    }
    if (btagCharmEfficiencyFile) {
        btagCharmEfficiencyFile->Close();
        delete btagCharmEfficiencyFile;
    }
    if (btagLightJetsEfficiencyFile) {
        btagLightJetsEfficiencyFile->Close();
        delete btagLightJetsEfficiencyFile;
    }
    if (eleTrigSFFile) {
        eleTrigSFFile->Close();
        delete eleTrigSFFile;
    }
    if (muTrigSFFile) {
        muTrigSFFile->Close();
        delete muTrigSFFile;
    }
    if (eleTrigEffFile) {
        eleTrigEffFile->Close();
        delete eleTrigEffFile;
    }
    if (muTrigEffFile) {
        muTrigEffFile->Close();
        delete muTrigEffFile;
    }
    for (uint i=0;i < JetCorrector.size() ; i++) {
      if (JetCorrector[i]) delete JetCorrector[i];
      if (JetResolutionParameters[i]) delete JetResolutionParameters[i];
      if (JetResolutionCalculator[i]) delete JetResolutionCalculator[i];
    }
    if (btagcalib) delete btagcalib;
    if (btagreader) delete btagreader;
    if (btagreader_up) delete btagreader_up;
    if (btagreader_do) delete btagreader_do;
    if (btagreaderMistag) delete btagreaderMistag;
    if (btagreaderMistag_up) delete btagreaderMistag_up;
    if (btagreaderMistag_do) delete btagreaderMistag_do;
    if (btagreaderfastsim) delete btagreaderfastsim;
    if (btagreaderfastsim_up) delete btagreaderfastsim_up;
    if (btagreaderfastsim_do) delete btagreaderfastsim_do;
    if (puppiSoftDropCorrFile) {
        puppiSoftDropCorrFile->Close();
        delete puppiSoftDropCorrFile;
    }
    if (wTopTagEffFile) {
        wTopTagEffFile->Close();
        delete wTopTagEffFile;
    }
}

// Retrieves CMSSW_BASE and stores in variable cmsswPath
void RazorHelper::loadCMSSWPath() {
    char* cmsswPathChar = getenv("CMSSW_BASE");
    if (cmsswPathChar == NULL) {
        std::cout << "Warning in RazorHelper::loadCMSSWPath : CMSSW_BASE not detected." << std::endl;
        cmsswPath = "";
    }
    cmsswPath = std::string(cmsswPathChar);
}

void RazorHelper::loadTag_Null() {
    std::cout << "Warning: initializing all RazorHelper files and histograms to 0" << std::endl;

    // pileup weights
    pileupWeightFile = 0;
    pileupWeightHist = 0;
    pileupWeightSysUpHist = 0;
    pileupWeightSysDownHist = 0;
    // electron efficiencies and scale factors
    eleTightEfficiencyFile = 0;
    eleVetoEfficiencyFile = 0;
    eleLooseEfficiencyFile = 0;
    eleEffSFFile = 0;
    vetoEleEffSFFile = 0;
    eleTightEfficiencyHist = 0;
    eleVetoEfficiencyHist = 0;
    eleLooseEfficiencyHist = 0;
    eleTightEffFastsimSFHist = 0;
    eleVetoEffFastsimSFHist = 0;
    eleLooseEffFastsimSFHist = 0;
    eleTightEffSFHist = 0;
    eleVetoEffSFHist = 0;
    eleLooseEffSFHist = 0;

    // muon efficiencies and scale factors
    muTightEfficiencyFile = 0;
    muVetoEfficiencyFile = 0;
    muLooseEfficiencyFile = 0;
    muEffSFFile = 0;
    vetoMuEffSFFile = 0;
    muTightEfficiencyHist = 0;
    muVetoEfficiencyHist = 0;
    muLooseEfficiencyHist = 0;
    muTightEffFastsimSFHist = 0;
    muVetoEffFastsimSFHist = 0;
    muLooseEffFastsimSFHist = 0;
    muTightEffSFHist = 0;
    muVetoEffSFHist = 0;
    muLooseEffSFHist = 0;

    // tau efficiencies and scale factors
    tauEfficiencyFile = 0;
    tauLooseEfficiencyHist = 0;

    // b-tag efficiencies and scale factors
    btagEfficiencyFile = 0;
    btagCharmEfficiencyFile = 0;
    btagLightJetsEfficiencyFile = 0;
    btagMediumEfficiencyHist = 0;
    btagMediumCharmEfficiencyHist = 0;
    btagMediumLightJetsEfficiencyHist = 0;

    // single lepton trigger scale factors
    eleTrigSFFile = 0;
    eleTrigSFHist = 0;
    muTrigSFFile = 0;
    muTrigSFHist = 0;

    // single lepton trigger efficiencies
    eleTrigEffFile = 0;
    eleTrigEffHist = 0;
    muTrigEffFile = 0;
    muTrigEffHist = 0;

}

void RazorHelper::loadHiggsPt() { //same for all years
    // pileup weights
    std::cout << "RazorHelper: loading higgs pt weight histograms" << std::endl;
    higgsPtWeightFile = TFile::Open("ggH_HiggsPtReweight_NNLOPS.root");
    higgsPtWeightHist = (TH1F*)higgsPtWeightFile->Get("higgsPtReweight");
    for (int i = 0; i < 9; i++)
    {
      higgsPtWeightSysHist[i] = (TH1F*)higgsPtWeightFile->Get(Form("higgsPtReweightSys%d", i));

    }
}

////////////////////////////////////////////////
//  2015
////////////////////////////////////////////////

void RazorHelper::loadTag_Razor2015() {
    loadPileup_Razor2015();
    loadLepton_Razor2015();
    loadJECs_Razor2015();
    loadBTag_Razor2015();
    loadTrigger_Razor2015();
}

void RazorHelper::loadPileup_Razor2015() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PileupWeights/PileupReweight_Spring15MCTo2015Data.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
}

void RazorHelper::loadLepton_Razor2015(){

    // electron efficiencies and scale factors
    std::cout << "RazorHelper: loading electron efficiency histograms" << std::endl;
    // For 2015 the lepton efficiency histograms for tight and veto IDs in fullsim are both
    // stored in the fastsim-to-fullsim scale factor files.
    eleTightEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/ElectronEffFastsimToFullsimCorrectionFactors.2015.root");
    eleVetoEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/ElectronEffFastsimToFullsimCorrectionFactors.2015.root");
    eleEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_TightElectronSelectionEffDenominatorReco_2015Final_Golden.root");
    vetoEleEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_VetoElectronSelectionEffDenominatorReco_2015Final_Golden.root");
    eleTightEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("ElectronEff_Tight_Fullsim");
    eleVetoEfficiencyHist = (TH2D*)eleVetoEfficiencyFile->Get("ElectronEff_Veto_Fullsim");
    eleTightEffFastsimSFHist = (TH2D*)eleTightEfficiencyFile->Get("ElectronTight_FastsimScaleFactor");
    eleVetoEffFastsimSFHist = (TH2D*)eleVetoEfficiencyFile->Get("ElectronVeto_FastsimScaleFactor");
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorReco");
    eleVetoEffSFHist = (TH2D*)vetoEleEffSFFile->Get("ScaleFactor_VetoElectronSelectionEffDenominatorReco");
    // In 2015 it was not necessary to apply any reco tracking scale factors
    // (this was introduced in 2016 to adjust for the impact of the 'HIP' effect)
    eleGSFTrackEffSFHist = 0;
    eleGSFTrackEffHist = 0;

    // muon efficiencies and scale factors
    std::cout << "RazorHelper: loading muon efficiency histograms" << std::endl;
    // For 2015 the lepton efficiency histograms for tight and veto IDs in fullsim are both
    // stored in the fastsim-to-fullsim scale factor files.
    muTightEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/MuonEffFastsimToFullsimCorrectionFactors.2015.root");
    muVetoEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/MuonEffFastsimToFullsimCorrectionFactors.2015.root");
    muEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_TightMuonSelectionEffDenominatorReco_2015Final_Golden.root");
    vetoMuEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_VetoMuonSelectionEffDenominatorReco_2015Final_Golden.root");
    muTightEfficiencyHist = (TH2D*)muTightEfficiencyFile->Get("MuonEff_Tight_Fullsim");
    muVetoEfficiencyHist = (TH2D*)muVetoEfficiencyFile->Get("MuonEff_Veto_Fullsim");
    muTightEffFastsimSFHist = (TH2D*)muTightEfficiencyFile->Get("MuonTight_FastsimScaleFactor");
    muVetoEffFastsimSFHist = (TH2D*)muVetoEfficiencyFile->Get("MuonVeto_FastsimScaleFactor");
    muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorReco");
    muVetoEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorReco");
    // In 2015 it was not necessary to apply any reco tracking scale factors
    // (this was introduced in 2016 to adjust for the impact of the 'HIP' effect)
    muTrackEffSFHist = 0;
    muTrackEffHist = 0;

    // tau efficiencies and scale factors
    std::cout << "RazorHelper: loading tau efficiency histograms" << std::endl;
    tauEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/TauEffFastsimToFullsimCorrectionFactors.2015.root");
    tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");

}

void RazorHelper::loadJECs_Razor2015() {
    std::cout << "RazorHelper: loading jet energy correction constants" << std::endl;
    // initialize
    std::string jecPathname = "./";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();

    if (isData) {

      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Fall15_25nsV2_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }
    else if (isFastsim) {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fastsim_MCRUN2_74_V9_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fastsim_MCRUN2_74_V9_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fastsim_MCRUN2_74_V9_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Fastsim_MCRUN2_74_V9_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));

    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }

}


void RazorHelper::loadBTag_Razor2015() {
    // b-tag efficiencies and scale factors
    std::cout << "RazorHelper: loading btag efficiency histograms" << std::endl;
    btagEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/BTagEffFastsimToFullsimCorrectionFactors.2015.root");
    btagCharmEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/CharmJetBTagEffFastsimToFullsimCorrectionFactors.2015.root");
    btagLightJetsEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/LightJetBTagEffFastsimToFullsimCorrectionFactors.2015.root");
    btagMediumEfficiencyHist = (TH2D*)btagEfficiencyFile->Get("BTagEff_Medium_Fullsim");
    btagMediumCharmEfficiencyHist = (TH2D*)btagCharmEfficiencyFile->Get("BTagEff_Medium_Fullsim");
    btagMediumLightJetsEfficiencyHist = (TH2D*)btagLightJetsEfficiencyFile->Get("BTagEff_Medium_Fullsim");

    // Fullsim
    btagcalib = new BTagCalibration("csvv2","./CSVv2.csv");

    btagreader = new BTagCalibrationReader(btagcalib,               // calibration instance
                                           BTagEntry::OP_MEDIUM,     // operating point
				           "mujets",                 // measurement type
				           "central");               // systematics type
    btagreader_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets", "up");  // sys up
    btagreader_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets", "down");  // sys down
    btagreaderMistag = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "central");
    btagreaderMistag_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "up");    // sys up
    btagreaderMistag_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "down");  // sys down

    // Fastsim
    // btagcalibfastsim = new BTagCalibration("csvv2", Form("%s/CSV_13TEV_Combined_20_11_2015.csv",bTagPathname.c_str()));
    btagcalibfastsim = new BTagCalibration("csvv2", "./CSV_13TEV_Combined_20_11_2015.csv");
    btagreaderfastsim = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "central");
    btagreaderfastsim_up = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "up");
    btagreaderfastsim_do = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "down");
}

void RazorHelper::loadTrigger_Razor2015() {
    // single lepton trigger scale factors
    std::cout << "RazorHelper: loading trigger efficiency histograms" << std::endl;
    eleTrigSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2015Final_Golden.root");
    eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");
    muTrigSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2015Final_Golden.root");
    muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");

    eleTrigEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/Spring15MC/SingleElectronTriggerEfficiencyFromFullsim.root");
    eleTrigEffHist = (TH2D*)eleTrigEffFile->Get("hEffEtaPt");
    muTrigEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/Spring15MC/SingleMuonTriggerEfficiencyFromFullsim.root");
    muTrigEffHist = (TH2D*)muTrigEffFile->Get("hEffEtaPt");

    //get trigger numbers
    dileptonTriggerNums = std::vector<int>(8);
    hadronicTriggerNums = std::vector<int>(11);
    if( isData ) {
        dileptonTriggerNums = { 41,43,30,31,47,48,49,50 };
        singleLeptonTriggerNums = std::vector<int>(13);
        singleLeptonTriggerNums = { 2,7,12,11,15,22,23,24,25,26,27,28,29 };
        hadronicTriggerNums = { 134,135,136,137,138,139,140,141,142,143,144 };
    }
    else {
        dileptonTriggerNums = { 41,43,30,31,47,48,49,50 };
        singleLeptonTriggerNums = std::vector<int>(11);
        singleLeptonTriggerNums = { 2,7,12,11,15,18,19,20,21,28,29 };
        hadronicTriggerNums = { 134,135,136,137,138,139,140,141,142,143,144 };
    }
}

////////////////////////////////////////////////
//  2015 76X ReReco
////////////////////////////////////////////////

void RazorHelper::loadTag_Razor2015_76X() {
    loadPileup_Razor2015_76X();
    loadLepton_Razor2015(); //we'll use the same here for now
    loadPhoton_Razor2015_76X();
    loadJECs_Razor2015_76X();   //
    loadBTag_Razor2015();   //we'll use the same here for now, but this needs to be updated
    loadTrigger_Razor2015();//use the same for now
}

void RazorHelper::loadPileup_Razor2015_76X() {
    // pileup weights
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PileupWeights/PileupReweight2015_7_6.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
}

void RazorHelper::loadPhoton_Razor2015_76X(){
    // photon efficiency scale factors
    std::cout << "RazorHelper: loading photon efficiency scale factor histograms" << std::endl;
    phoEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2015/efficiency_results_PhoLooseEffDenominatorReco_2015.root");
    phoLooseEffSFHist = (TH2D*)phoEffSFFile->Get("ScaleFactor_PhoLooseEffDenominatorReco");
}

void RazorHelper::loadJECs_Razor2015_76X() {
    std::cout << "RazorHelper: loading jet energy correction constants" << std::endl;
    // initialize
    // load JEC parameters
    std::string jecPathname = "./";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();

    if (isData) {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Fall15_25nsV2_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }
    else if (isFastsim) {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fastsim_MCRUN2_74_V9_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fastsim_MCRUN2_74_V9_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fastsim_MCRUN2_74_V9_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Fastsim_MCRUN2_74_V9_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersTemp = std::vector<JetCorrectorParameters> ();
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersTemp.push_back(JetCorrectorParameters(
                  Form("%s/Fall15_25nsV2_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersTemp = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorTemp = new FactorizedJetCorrector(correctionParametersTemp);
      std::string jecUncPath = jecPathname+"/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncTemp = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorTemp = new SimpleJetResolution(*JetResolutionParametersTemp);

      correctionParameters.push_back(correctionParametersTemp);
      JetResolutionParameters.push_back(JetResolutionParametersTemp);
      JetCorrector.push_back( JetCorrectorTemp );
      JetResolutionCalculator.push_back(JetResolutionCalculatorTemp);
      jecUnc.push_back(jecUncTemp);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 0, 99999999 ));
    }

}

////////////////////////////////////////////////
////  2016 17Aug2017 Rereco
//////////////////////////////////////////////////
void RazorHelper::loadTag_Razor2016_Source2018() {
    loadPileup_Razor2016_07Aug2017Rereco();
    // loadLepton_Razor2016_MoriondRereco();
    // loadPhoton_Razor2016_07Aug2017Rereco_DelayedPhoton();
    // loadBTag_Razor2016_MoriondRereco();
    loadTrigger_Razor2016_07Aug2017Rereco();
    loadJECs_Razor2016_07Aug2017Rereco();
    loadHiggsPt();
    // loadAK8JetTag_Razor2016_MoriondRereco();
}
void RazorHelper::loadPileup_Razor2016_Source2018() {
    // pileup weights
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;


    pileupWeightFile = TFile::Open("PileupReweight_source18_target16.root.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
    std::cout << "PileupReweight_source18_target16.root.root\n";

}
void RazorHelper::loadTag_Razor2016_07Aug2017Rereco() {
    loadPileup_Razor2016_07Aug2017Rereco();
    // loadLepton_Razor2016_MoriondRereco();
    // loadPhoton_Razor2016_07Aug2017Rereco_DelayedPhoton();
    // loadBTag_Razor2016_MoriondRereco();
    loadTrigger_Razor2016_07Aug2017Rereco();
    loadJECs_Razor2016_07Aug2017Rereco();
    loadHiggsPt();
    // loadAK8JetTag_Razor2016_MoriondRereco();
}


void RazorHelper::loadPileup_Razor2016_07Aug2017Rereco() {
    // pileup weights
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;


    pileupWeightFile = TFile::Open("PileupReweight_MC_Summer16_ggH_HToSSTobbbb_MH-125_TuneCUETP8M1_13TeV-powheg-pythia8.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
    std::cout << "PileupReweight_MC_Summer16_ggH_HToSSTobbbb_MH-125_TuneCUETP8M1_13TeV-powheg-pythia8.root\n";




}
void RazorHelper::loadTrigger_Razor2016_07Aug2017Rereco() {
    // single lepton trigger scale factors
    std::cout << "RazorHelper: loading trigger efficiency histograms" << std::endl;
    metTriggerSFFile = TFile::Open("METTriggers_SF.root");
    metTriggerSFHist = (TH1F*)metTriggerSFFile->Get("trigger_efficiency_Fall17");

}
void RazorHelper::loadJECs_Razor2016_07Aug2017Rereco() {
    std::cout << "RazorHelper: loading jet energy correction constants, using Summer16_07Aug2017_V11." << std::endl;
    // initialize
    std::string jecPathname = "JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    std::cout << "here1\n";
    if (isData) {
      //IOV: 2016BCD
      std::vector<JetCorrectorParameters> correctionParametersBCD = std::vector<JetCorrectorParameters> ();
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersBCD = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorBCD = new FactorizedJetCorrector(correctionParametersBCD);
      std::string jecUncPathBCD = jecPathname+"/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017BCD_V11_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncBCD = new JetCorrectionUncertainty(jecUncPathBCD);
      SimpleJetResolution* JetResolutionCalculatorBCD = new SimpleJetResolution(*JetResolutionParametersBCD);

      correctionParameters.push_back(correctionParametersBCD);
      JetResolutionParameters.push_back(JetResolutionParametersBCD);
      JetCorrector.push_back( JetCorrectorBCD );
      JetResolutionCalculator.push_back(JetResolutionCalculatorBCD);
      jecUnc.push_back(jecUncBCD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 272007, 276811 ));

      //IOV: 2016EF
      std::vector<JetCorrectorParameters> correctionParametersEF = std::vector<JetCorrectorParameters> ();
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017EF_V11_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017EF_V11_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017EF_V11_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017EF_V11_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersEF = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorEF = new FactorizedJetCorrector(correctionParametersEF);
      std::string jecUncPathEF = jecPathname+"/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017EF_V11_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncEF = new JetCorrectionUncertainty(jecUncPathEF);
      SimpleJetResolution* JetResolutionCalculatorEF = new SimpleJetResolution(*JetResolutionParametersEF);

      correctionParameters.push_back(correctionParametersEF);
      JetResolutionParameters.push_back(JetResolutionParametersEF);
      JetCorrector.push_back( JetCorrectorEF );
      JetResolutionCalculator.push_back(JetResolutionCalculatorEF);
      jecUnc.push_back(jecUncEF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 276831, 278808 ));

      //IOV: 2016GH
      std::vector<JetCorrectorParameters> correctionParametersGH = std::vector<JetCorrectorParameters> ();
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017GH_V11_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017GH_V11_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017GH_V11_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017GH_V11_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersGH = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorGH = new FactorizedJetCorrector(correctionParametersGH);
      std::string jecUncPathGH = jecPathname+"/Summer16_07Aug2017_V11_DATA/Summer16_07Aug2017GH_V11_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncGH = new JetCorrectionUncertainty(jecUncPathGH);
      SimpleJetResolution* JetResolutionCalculatorGH = new SimpleJetResolution(*JetResolutionParametersGH);

      correctionParameters.push_back(correctionParametersGH);
      JetResolutionParameters.push_back(JetResolutionParametersGH);
      JetCorrector.push_back( JetCorrectorGH );
      JetResolutionCalculator.push_back(JetResolutionCalculatorGH);
      jecUnc.push_back(jecUncGH);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 278820, 284044 ));

    }
    else if (isFastsim) {
      std::cout << "Fastsim JEC\n";

      std::vector<JetCorrectorParameters> correctionParametersFastsim = std::vector<JetCorrectorParameters> ();
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersFastsim = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorFastsim = new FactorizedJetCorrector(correctionParametersFastsim);
      std::string jecUncPath = jecPathname+"/Spring16_FastSimV1_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncFastsim = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorFastsim = new SimpleJetResolution(*JetResolutionParametersFastsim);

      correctionParameters.push_back(correctionParametersFastsim);
      JetResolutionParameters.push_back(JetResolutionParametersFastsim);
      JetCorrector.push_back( JetCorrectorFastsim );
      JetResolutionCalculator.push_back(JetResolutionCalculatorFastsim);
      jecUnc.push_back(jecUncFastsim);
      JetCorrectionsIOV.push_back( std::pair<int,int>( -1, 99999999 ));
    }
    else {
      std::cout << "Loading Jet Energy Corrections: Summer16_23Sep2016V6_MC \n";
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));

      JetCorrectorParameters *JetResolutionParametersMC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorMC = new SimpleJetResolution(*JetResolutionParametersMC);

      correctionParameters.push_back(correctionParametersMC);
      JetResolutionParameters.push_back(JetResolutionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorMC);
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( -1, 99999999 ));
    }
}



////////////////////////////////////////////////
//  2016 Moriond Rereco
////////////////////////////////////////////////
void RazorHelper::loadTag_Razor2016_MoriondRereco() {
    loadPileup_Razor2016_MoriondRereco();
    loadLepton_Razor2016_MoriondRereco();
    loadPhoton_Razor2016_MoriondRereco();
    loadBTag_Razor2016_MoriondRereco();
    loadTrigger_Razor2016_MoriondRereco();
    //loadJECs_Razor2016_MoriondRereco();
    loadJECs_Razor2016_03Feb2017Rereco();
    loadAK8JetTag_Razor2016_MoriondRereco();
}


void RazorHelper::loadPileup_Razor2016_MoriondRereco() {
    // pileup weights
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;


    if (!isFastsim) {
      pileupWeightFile = TFile::Open("PileupReweight_Summer16_2016_36p2ifb.root");
      pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
      pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
      pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
      std::cout << "PileupReweight_Summer16_2016_36p2ifb.root\n";
    } else {
      pileupWeightFile = TFile::Open("root://eoscms:////store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PileupWeights/PileupReweight_2016_36p2ifb.root");
      pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
      pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
      pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
      std::cout << "PileupReweight_2016_36p2ifb.root\n";
    }

}

void RazorHelper::loadLepton_Razor2016_MoriondRereco(){

    // electron efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 electron efficiency histograms" << std::endl;
    //eleTightEfficiencyFile = TFile::Open("root://eoscms:////store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    eleTightEfficiencyFile = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    //eleVetoEfficiencyFile = TFile::Open("root://eoscms:////store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    eleVetoEfficiencyFile = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    //eleEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightElectronSelectionEffDenominatorGen_2016_Rereco_Golden.root");
    eleEffSFFile = TFile::Open("efficiency_results_TightElectronSelectionEffDenominatorGen_2016_Rereco_Golden.root");
    //vetoEleEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoElectronSelectionEffDenominatorGen_2016_Rereco_Golden.root");
    vetoEleEffSFFile = TFile::Open("efficiency_results_VetoElectronSelectionEffDenominatorGen_2016_Rereco_Golden.root");
    //eleGSFTrackEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiencySF_muEleTracking_2016_average.root");
    eleGSFTrackEffSFFile = TFile::Open("efficiencySF_muEleTracking_2016_average.root");
    //eleGSFTrackEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Reco_Fullsim.root");
    eleGSFTrackEffFile = TFile::Open("Efficiency_PromptElectron_TTJets_25ns_Reco_Fullsim.root");
    //eleTightEffFastsimSFFile = TFile::Open("root://eoscms:////store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    eleTightEffFastsimSFFile = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    //eleVetoEffFastsimSFFile = TFile::Open("root://eoscms:////store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    eleVetoEffFastsimSFFile = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.2016.root");

    eleTightEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("ElectronEff_Tight_Fullsim");
    eleVetoEfficiencyHist = (TH2D*)eleVetoEfficiencyFile->Get("ElectronEff_Veto_Fullsim");

    // We don't have ID scale factors for Fastsim yet.
    eleTightEffFastsimSFHist =  (TH2D*)eleTightEffFastsimSFFile->Get("ElectronTight_FastsimScaleFactor");
    eleVetoEffFastsimSFHist = (TH2D*)eleVetoEffFastsimSFFile->Get("ElectronEff_Veto_Fullsim");
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorGen");
    //add for now 20190124, no fastsim
    eleLooseEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorGen");
    eleLooseEffFastsimSFHist =  (TH2D*)eleTightEffFastsimSFFile->Get("ElectronTight_FastsimScaleFactor");
    eleVetoEffSFHist = (TH2D*)vetoEleEffSFFile->Get("ScaleFactor_VetoElectronSelectionEffDenominatorGen");
    eleGSFTrackEffHist = (TH2D*)eleGSFTrackEffFile->Get("Efficiency_PtEta");
    // These scale factors are weighted according to the fraction of the 2016 run affected
    // by the 'HIP' issue, under the assumption that tracking scale factors are 1 for runs
    // not affected by the 'HIP'.
    eleGSFTrackEffSFHist = (TH2D*)eleGSFTrackEffSFFile->Get("h2_scaleFactorsEGamma");

    // muon efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 muon efficiency histograms" << std::endl;
    muTightEfficiencyFile = TFile::Open("MuonEffFastsimToFullsimCorrectionFactors.2016.root");
    muVetoEfficiencyFile = TFile::Open("MuonEffFastsimToFullsimCorrectionFactors.2016.root");
    muEffSFFile = TFile::Open("efficiency_results_TightMuonSelectionEffDenominatorGen_2016_Rereco_Golden.root");
    vetoMuEffSFFile = TFile::Open("efficiency_results_VetoMuonSelectionEffDenominatorGen_2016_Rereco_Golden.root");
    muTrackEffSFFile = TFile::Open("efficiencySF_muEleTracking_2016_average.root");
    muTrackEffFile = TFile::Open("Efficiency_PromptMuon_TTJets_25ns_Reco_Fullsim.root");
    muTightEffFastsimSFFile = TFile::Open("MuonEffFastsimToFullsimCorrectionFactors.2016.root");
    muVetoEffFastsimSFFile = TFile::Open("MuonEffFastsimToFullsimCorrectionFactors.2016.root");

    muTightEfficiencyHist = (TH2D*)muTightEfficiencyFile->Get("MuonEff_Tight_Fullsim");
    muVetoEfficiencyHist = (TH2D*)muVetoEfficiencyFile->Get("MuonEff_Veto_Fullsim");
    // We don't have ID scale factors for Fastsim yet.
    muTightEffFastsimSFHist = (TH2D*)muTightEffFastsimSFFile->Get("MuonTight_FastsimScaleFactor");
    muVetoEffFastsimSFHist = (TH2D*)muVetoEffFastsimSFFile->Get("MuonVeto_FastsimScaleFactor");
    muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorGen");
    muVetoEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorGen");
    //add for now 20190124, no fastsim
    muLooseEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorGen");
    muLooseEffFastsimSFHist = (TH2D*)muVetoEffFastsimSFFile->Get("MuonVeto_FastsimScaleFactor");
    muTrackEffHist = (TH2D*)muTrackEffFile->Get("Efficiency_PtEta");
    // These scale factors are weighted according to the fraction of the 2016 run affected
    // by the 'HIP' issue, under the assumption that tracking scale factors are 1 for runs
    // not affected by the 'HIP'.
    muTrackEffSFHist = (TH2D*)muTrackEffSFFile->Get("muon");

    // tau efficiencies and scale factors
    std::cout << "RazorHelper: loading tau efficiency histograms" << std::endl;
    tauEfficiencyFile = TFile::Open("TauEffFastsimToFullsimCorrectionFactors.2016.root");
    tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");

}


void RazorHelper::loadPhoton_Razor2016_MoriondRereco(){
    // photon efficiency scale factors
    std::cout << "RazorHelper: loading photon efficiency scale factor histograms" << std::endl;
    phoEffSFFile = TFile::Open("efficiency_results_PhoLooseEffDenominatorReco_2016_Rereco.root");
    phoLooseEffSFHist = (TH2D*)phoEffSFFile->Get("ScaleFactor_PhoLooseEffDenominatorReco");

    phoEffFastsimSFFile = TFile::Open("PhotonEffFastsimToFullsimCorrectionFactors.2016.root");
    phoLooseEffFastsimSFHist = (TH2D*)phoEffFastsimSFFile->Get("ElectronLoose_FastsimScaleFactor");
}

void RazorHelper::loadPhoton_Razor2016_07Aug2017Rereco_DelayedPhoton(){
    std::cout << "RazorHelper: loading photon efficiency scale factor histograms" << std::endl;
    phoEffSFFile = TFile::Open("efficiency_results_PhoLooseEffDenominatorReco_2016_Rereco_DelayedPhoton.root");
    phoLooseEffSFHist = (TH2D*)phoEffSFFile->Get("ScaleFactor_PhoLooseEffDenominatorReco");
    phoTightEffSFFile = TFile::Open("efficiency_results_PhoTightEffDenominatorReco_2016_Rereco_DelayedPhoton.root");
    phoTightEffSFHist = (TH2D*)phoTightEffSFFile->Get("ScaleFactor_PhoTightEffDenominatorReco");

    phoEffFastsimSFFile = TFile::Open("PhotonEffFastsimToFullsimCorrectionFactors.2016.root");
    phoLooseEffFastsimSFHist = (TH2D*)phoEffFastsimSFFile->Get("ElectronLoose_FastsimScaleFactor");
}

void RazorHelper::loadBTag_Razor2016_MoriondRereco() {
    // b-tag efficiencies and scale factors
    std::cout << "RazorHelper: loading btag efficiency histograms" << std::endl;
    btagEfficiencyFile = TFile::Open("Efficiency_BJets_25ns_CSVM_Fullsim_80X.root");
    btagCharmEfficiencyFile = TFile::Open("Efficiency_CJets_25ns_CSVM_Fullsim_80X.root");
    btagLightJetsEfficiencyFile = TFile::Open("Efficiency_LightJets_25ns_CSVM_Fullsim_80X.root");
    btagMediumEfficiencyHist = (TH2D*)btagEfficiencyFile->Get("Efficiency_PtEta");
    btagMediumCharmEfficiencyHist = (TH2D*)btagCharmEfficiencyFile->Get("Efficiency_PtEta");
    btagMediumLightJetsEfficiencyHist = (TH2D*)btagLightJetsEfficiencyFile->Get("Efficiency_PtEta");

    // Fullsim
    btagcalib = new BTagCalibration("csvv2", "./CSVv2_Moriond17_B_H.csv");
    btagreader = new BTagCalibrationReader(btagcalib,               // calibration instance
                                           BTagEntry::OP_MEDIUM,     // operating point
				           "comb",                 // measurement type
				           "central");               // systematics type
    btagreader_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "up");  // sys up
    btagreader_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "down");  // sys down
    btagreaderMistag = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "incl", "central");
    btagreaderMistag_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "incl", "up");    // sys up
    btagreaderMistag_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "incl", "down");  // sys down

    // Fastsim
    btagcalibfastsim = new BTagCalibration("csvv2", "./fastsim_csvv2_ttbar_26_1_2017.csv");
    btagreaderfastsim = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "central");
    btagreaderfastsim_up = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "up");
    btagreaderfastsim_do = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "down");
}

void RazorHelper::loadTrigger_Razor2016_MoriondRereco() {
    // single lepton trigger scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 trigger efficiency histograms" << std::endl;
    eleTrigSFFile = TFile::Open("efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2016_Rereco_Golden.root");
    eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");

    muTrigSFFile = TFile::Open("efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2016_Rereco_Golden.root");
    muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");

    eleTrigEffFile = TFile::Open("SingleElectronTriggerEfficiency_2016_Rereco_Golden.root");
    eleTrigEffHist = (TH2D*)eleTrigEffFile->Get("hEffEtaPt");

    muTrigEffFile = TFile::Open("SingleMuonTriggerEfficiency_2016_Rereco_Golden.root");
    muTrigEffHist = (TH2D*)muTrigEffFile->Get("hEffEtaPt");

    //diphoton trigger scale factors
    diphotonTrigLeadingLegEffFile = TFile::Open("PhoHLTLeadingLegEffDenominatorLoose_2016_Rereco.root");
    diphotonTrigLeadingLegEffHist = (TH2D*)diphotonTrigLeadingLegEffFile->Get("hEffEtaPt");
    diphotonTrigTrailingLegEffFile = TFile::Open("PhoHLTTrailingLegEffDenominatorLoose_2016_Rereco.root");
    diphotonTrigTrailingLegEffHist = (TH2D*)diphotonTrigTrailingLegEffFile->Get("hEffEtaPt");

    diphotonTrigLeadingLegEffSFFile = TFile::Open("efficiency_results_PhoHLTLeadingLegEffDenominatorLoose_2016_Rereco.root");
    diphotonTrigLeadingLegEffSFHist = (TH2D*)diphotonTrigLeadingLegEffSFFile->Get("ScaleFactor_PhoHLTLeadingLegEffDenominatorLoose");
    diphotonTrigTrailingLegEffSFFile = TFile::Open("efficiency_results_PhoHLTTrailingLegEffDenominatorLoose_2016_Rereco.root");
    diphotonTrigTrailingLegEffSFHist = (TH2D*)diphotonTrigTrailingLegEffSFFile->Get("ScaleFactor_PhoHLTTrailingLegEffDenominatorLoose");

    //get trigger numbers
    // (we are using the same list of trigger numbers for data and MC for 2016)
    dileptonTriggerNums = { 44,45,57,59,64,65,66,67,68 };
    singleLeptonTriggerNums = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43 };
    hadronicTriggerNums = { 163,164,165,166,167,168,169,170,171,172,173,174,175,176 };
}

void RazorHelper::loadTrigger_Razor2016_07Aug2017Rereco_DelayedPhoton(){
std::cout << "RazorHelper: loading 2016 trigger efficiency histograms" << std::endl;
    eleTrigSFFile = TFile::Open("efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2016_Rereco_Golden.root");
    eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");

    muTrigSFFile = TFile::Open("efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2016_Rereco_Golden.root");
    muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");

    eleTrigEffFile = TFile::Open("SingleElectronTriggerEfficiency_2016_Rereco_Golden.root");
    eleTrigEffHist = (TH2D*)eleTrigEffFile->Get("hEffEtaPt");

    muTrigEffFile = TFile::Open("SingleMuonTriggerEfficiency_2016_Rereco_Golden.root");
    muTrigEffHist = (TH2D*)muTrigEffFile->Get("hEffEtaPt");

    diphotonTrigLeadingLegEffFile = TFile::Open("PhoHLTLeadingLegEffDenominatorLoose_2016_Rereco_DelayedPhoton.root");
    diphotonTrigLeadingLegEffHist = (TH2D*)diphotonTrigLeadingLegEffFile->Get("hEffEtaPt");
    diphotonTrigTrailingLegEffFile = TFile::Open("PhoHLTTrailingLegEffDenominatorLoose_2016_Rereco_DelayedPhoton.root");
    diphotonTrigTrailingLegEffHist = (TH2D*)diphotonTrigTrailingLegEffFile->Get("hEffEtaPt");

    diphotonTrigLeadingLegEffSFFile = TFile::Open("efficiency_results_PhoHLTLeadingLegEffDenominatorLoose_2016_Rereco_DelayedPhoton.root");
    diphotonTrigLeadingLegEffSFHist = (TH2D*)diphotonTrigLeadingLegEffSFFile->Get("ScaleFactor_PhoHLTLeadingLegEffDenominatorLoose");
    diphotonTrigTrailingLegEffSFFile = TFile::Open("efficiency_results_PhoHLTTrailingLegEffDenominatorLoose_2016_Rereco_DelayedPhoton.root");
    diphotonTrigTrailingLegEffSFHist = (TH2D*)diphotonTrigTrailingLegEffSFFile->Get("ScaleFactor_PhoHLTTrailingLegEffDenominatorLoose");

    dileptonTriggerNums = { 44,45,57,59,64,65,66,67,68 };
    singleLeptonTriggerNums = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43 };
    hadronicTriggerNums = { 163,164,165,166,167,168,169,170,171,172,173,174,175,176 };
}


void RazorHelper::loadJECs_Razor2016_MoriondRereco() {
    std::cout << "RazorHelper: loading jet energy correction constants, using Summer16_23Sep2016_V3." << std::endl;
    // initialize
    std::string jecPathname = "./";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    std::cout << "here1\n";
    if (isData) {
      //IOV: 2016BCD
      std::vector<JetCorrectorParameters> correctionParametersBCD = std::vector<JetCorrectorParameters> ();
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersBCD = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorBCD = new FactorizedJetCorrector(correctionParametersBCD);
      std::string jecUncPathBCD = jecPathname+"/Summer16_23Sep2016BCDV3_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncBCD = new JetCorrectionUncertainty(jecUncPathBCD);
      SimpleJetResolution* JetResolutionCalculatorBCD = new SimpleJetResolution(*JetResolutionParametersBCD);

      correctionParameters.push_back(correctionParametersBCD);
      JetResolutionParameters.push_back(JetResolutionParametersBCD);
      JetCorrector.push_back( JetCorrectorBCD );
      JetResolutionCalculator.push_back(JetResolutionCalculatorBCD);
      jecUnc.push_back(jecUncBCD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 276811 ));

      //IOV: 2016E
      std::vector<JetCorrectorParameters> correctionParametersEF = std::vector<JetCorrectorParameters> ();
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016EFV3_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016EFV3_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016EFV3_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016EFV3_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersEF = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorEF = new FactorizedJetCorrector(correctionParametersEF);
      std::string jecUncPathEF = jecPathname+"/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016EFV3_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncEF = new JetCorrectionUncertainty(jecUncPathEF);
      SimpleJetResolution* JetResolutionCalculatorEF = new SimpleJetResolution(*JetResolutionParametersEF);

      correctionParameters.push_back(correctionParametersEF);
      JetResolutionParameters.push_back(JetResolutionParametersEF);
      JetCorrector.push_back( JetCorrectorEF );
      JetResolutionCalculator.push_back(JetResolutionCalculatorEF);
      jecUnc.push_back(jecUncEF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 276831, 278801 ));

      //IOV: 2016G
      std::vector<JetCorrectorParameters> correctionParametersG = std::vector<JetCorrectorParameters> ();
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersG = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorG = new FactorizedJetCorrector(correctionParametersG);
      std::string jecUncPathG = jecPathname+"/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncG = new JetCorrectionUncertainty(jecUncPathG);
      SimpleJetResolution* JetResolutionCalculatorG = new SimpleJetResolution(*JetResolutionParametersG);

      correctionParameters.push_back(correctionParametersG);
      JetResolutionParameters.push_back(JetResolutionParametersG);
      JetCorrector.push_back( JetCorrectorG );
      JetResolutionCalculator.push_back(JetResolutionCalculatorG);
      jecUnc.push_back(jecUncG);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 278802, 280385 ));

      //IOV: 2016H
      std::vector<JetCorrectorParameters> correctionParametersH = std::vector<JetCorrectorParameters> ();
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016HV3_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016HV3_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016HV3_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersH = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorH = new FactorizedJetCorrector(correctionParametersH);
      std::string jecUncPathH = jecPathname+"/Summer16_23Sep2016V3_DATA/Summer16_23Sep2016HV3_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncH = new JetCorrectionUncertainty(jecUncPathH);
      SimpleJetResolution* JetResolutionCalculatorH = new SimpleJetResolution(*JetResolutionParametersH);

      correctionParameters.push_back(correctionParametersH);
      JetResolutionParameters.push_back(JetResolutionParametersH);
      JetCorrector.push_back( JetCorrectorH );
      JetResolutionCalculator.push_back(JetResolutionCalculatorH);
      jecUnc.push_back(jecUncH);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 280919, 99999999 ));

    }
    else if (isFastsim) {
      std::cout << "Fastsim JEC\n";

      std::vector<JetCorrectorParameters> correctionParametersFastsim = std::vector<JetCorrectorParameters> ();
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersFastsim = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorFastsim = new FactorizedJetCorrector(correctionParametersFastsim);
      std::string jecUncPath = jecPathname+"/Spring16_FastSimV1_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncFastsim = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorFastsim = new SimpleJetResolution(*JetResolutionParametersFastsim);

      correctionParameters.push_back(correctionParametersFastsim);
      JetResolutionParameters.push_back(JetResolutionParametersFastsim);
      JetCorrector.push_back( JetCorrectorFastsim );
      JetResolutionCalculator.push_back(JetResolutionCalculatorFastsim);
      jecUnc.push_back(jecUncFastsim);
      JetCorrectionsIOV.push_back( std::pair<int,int>( -1, 99999999 ));
    }
    else {
      std::cout << "Loading Jet Energy Corrections: Summer16_23Sep2016V3_MC \n";
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));

      JetCorrectorParameters *JetResolutionParametersMC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorMC = new SimpleJetResolution(*JetResolutionParametersMC);

      correctionParameters.push_back(correctionParametersMC);
      JetResolutionParameters.push_back(JetResolutionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorMC);
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( -1, 99999999 ));
    }
}



////////////////////////////////////////////////
//  2016 03Feb2017 Rereco
////////////////////////////////////////////////
void RazorHelper::loadTag_Razor2016_03Feb2017Rereco() {
    loadPileup_Razor2016_MoriondRereco();
    loadLepton_Razor2016_03Feb2017Rereco();
    loadPhoton_Razor2016_MoriondRereco();
    loadBTag_Razor2016_MoriondRereco();
    loadTrigger_Razor2016_MoriondRereco();
    loadJECs_Razor2016_03Feb2017Rereco();
    loadAK8JetTag_Razor2016_MoriondRereco();
}


void RazorHelper::loadLepton_Razor2016_03Feb2017Rereco(){

    // electron efficiencies and scale factors
    // LAST UPDATED: 24 Feb 2019
    std::cout << "RazorHelper: loading 2016 electron efficiency histograms" << std::endl;
    eleTightEfficiencyFile = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    eleVetoEfficiencyFile = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    eleEffSFFile = TFile::Open("efficiency_results_TightElectronSelectionEffDenominatorGen_2016_Rereco_Golden.root");
    looseEleEffSFFile = TFile::Open("efficiency_results_LooseElectronSelectionEffDenominatorGen_2016_23SepRereco.root");
    vetoEleEffSFFile = TFile::Open("efficiency_results_VetoElectronSelectionEffDenominatorGen_2016_Rereco_Golden.root");
    eleGSFTrackEffSFFile = TFile::Open("efficiencySF_muEleTracking_2016_average.root");
    eleGSFTrackEffFile = TFile::Open("Efficiency_PromptElectron_TTJets_25ns_Reco_Fullsim.root");
    eleTightEffFastsimSFFile = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    eleVetoEffFastsimSFFile = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.2016.root");

    eleTightEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("ElectronEff_Tight_Fullsim");
    eleVetoEfficiencyHist = (TH2D*)eleVetoEfficiencyFile->Get("ElectronEff_Veto_Fullsim");

    // We don't have ID scale factors for Fastsim yet.
    eleTightEffFastsimSFHist =  (TH2D*)eleTightEffFastsimSFFile->Get("ElectronTight_FastsimScaleFactor");
    eleVetoEffFastsimSFHist = (TH2D*)eleVetoEffFastsimSFFile->Get("ElectronEff_Veto_Fullsim");
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorGen");
    eleLooseEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_LooseElectronSelectionEffDenominatorGen");
    eleLooseEffFastsimSFHist =  (TH2D*)eleTightEffFastsimSFFile->Get("ElectronTight_FastsimScaleFactor");
    eleVetoEffSFHist = (TH2D*)vetoEleEffSFFile->Get("ScaleFactor_VetoElectronSelectionEffDenominatorGen");
    eleGSFTrackEffHist = (TH2D*)eleGSFTrackEffFile->Get("Efficiency_PtEta");
    // These scale factors are weighted according to the fraction of the 2016 run affected
    // by the 'HIP' issue, under the assumption that tracking scale factors are 1 for runs
    // not affected by the 'HIP'.
    eleGSFTrackEffSFHist = (TH2D*)eleGSFTrackEffSFFile->Get("h2_scaleFactorsEGamma");

    // muon efficiencies and scale factors
    // LAST UPDATED: 24 Feb 2019
    std::cout << "RazorHelper: loading 2016 muon efficiency histograms" << std::endl;
    muTightEfficiencyFile = TFile::Open("MuonEffFastsimToFullsimCorrectionFactors.2016.root");
    muVetoEfficiencyFile = TFile::Open("MuonEffFastsimToFullsimCorrectionFactors.2016.root");
    muEffSFFile = TFile::Open("efficiency_results_TightMuonSelectionEffDenominatorGen_2016_Rereco_Golden.root");
    vetoMuEffSFFile = TFile::Open("efficiency_results_VetoMuonSelectionEffDenominatorGen_2016_23SepRereco.root");
    muTrackEffSFFile = TFile::Open("efficiencySF_muEleTracking_2016_average.root");
    muTrackEffFile = TFile::Open("Efficiency_PromptMuon_TTJets_25ns_Reco_Fullsim.root");
    muTightEffFastsimSFFile = TFile::Open("MuonEffFastsimToFullsimCorrectionFactors.2016.root");
    muVetoEffFastsimSFFile = TFile::Open("MuonEffFastsimToFullsimCorrectionFactors.2016.root");

    muTightEfficiencyHist = (TH2D*)muTightEfficiencyFile->Get("MuonEff_Tight_Fullsim");
    muVetoEfficiencyHist = (TH2D*)muVetoEfficiencyFile->Get("MuonEff_Veto_Fullsim");
    // We don't have ID scale factors for Fastsim yet.
    muTightEffFastsimSFHist = (TH2D*)muTightEffFastsimSFFile->Get("MuonTight_FastsimScaleFactor");
    muVetoEffFastsimSFHist = (TH2D*)muVetoEffFastsimSFFile->Get("MuonVeto_FastsimScaleFactor");
    muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorGen");
    muVetoEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorGen");

    muLooseEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorGen");
    muLooseEffFastsimSFHist = (TH2D*)muVetoEffFastsimSFFile->Get("MuonVeto_FastsimScaleFactor");
    muTrackEffHist = (TH2D*)muTrackEffFile->Get("Efficiency_PtEta");
    // These scale factors are weighted according to the fraction of the 2016 run affected
    // by the 'HIP' issue, under the assumption that tracking scale factors are 1 for runs
    // not affected by the 'HIP'.
    muTrackEffSFHist = (TH2D*)muTrackEffSFFile->Get("muon");

    // tau efficiencies and scale factors
    std::cout << "RazorHelper: loading tau efficiency histograms" << std::endl;
    tauEfficiencyFile = TFile::Open("TauEffFastsimToFullsimCorrectionFactors.2016.root");
    tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");

}


void RazorHelper::loadJECs_Razor2016_03Feb2017Rereco() {
    std::cout << "RazorHelper: loading jet energy correction constants, using Summer16_23Sep2016_V4." << std::endl;
    // initialize
    std::string jecPathname = "./";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    std::cout << "here1\n";
    if (isData) {
      //IOV: 2016BCD
      std::vector<JetCorrectorParameters> correctionParametersBCD = std::vector<JetCorrectorParameters> ();
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersBCD = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorBCD = new FactorizedJetCorrector(correctionParametersBCD);
      std::string jecUncPathBCD = jecPathname+"/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncBCD = new JetCorrectionUncertainty(jecUncPathBCD);
      SimpleJetResolution* JetResolutionCalculatorBCD = new SimpleJetResolution(*JetResolutionParametersBCD);

      correctionParameters.push_back(correctionParametersBCD);
      JetResolutionParameters.push_back(JetResolutionParametersBCD);
      JetCorrector.push_back( JetCorrectorBCD );
      JetResolutionCalculator.push_back(JetResolutionCalculatorBCD);
      jecUnc.push_back(jecUncBCD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 276811 ));

      //IOV: 2016E
      std::vector<JetCorrectorParameters> correctionParametersEF = std::vector<JetCorrectorParameters> ();
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersEF.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersEF = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorEF = new FactorizedJetCorrector(correctionParametersEF);
      std::string jecUncPathEF = jecPathname+"/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncEF = new JetCorrectionUncertainty(jecUncPathEF);
      SimpleJetResolution* JetResolutionCalculatorEF = new SimpleJetResolution(*JetResolutionParametersEF);

      correctionParameters.push_back(correctionParametersEF);
      JetResolutionParameters.push_back(JetResolutionParametersEF);
      JetCorrector.push_back( JetCorrectorEF );
      JetResolutionCalculator.push_back(JetResolutionCalculatorEF);
      jecUnc.push_back(jecUncEF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 276831, 278801 ));

      //IOV: 2016G
      std::vector<JetCorrectorParameters> correctionParametersG = std::vector<JetCorrectorParameters> ();
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersG.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersG = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorG = new FactorizedJetCorrector(correctionParametersG);
      std::string jecUncPathG = jecPathname+"/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncG = new JetCorrectionUncertainty(jecUncPathG);
      SimpleJetResolution* JetResolutionCalculatorG = new SimpleJetResolution(*JetResolutionParametersG);

      correctionParameters.push_back(correctionParametersG);
      JetResolutionParameters.push_back(JetResolutionParametersG);
      JetCorrector.push_back( JetCorrectorG );
      JetResolutionCalculator.push_back(JetResolutionCalculatorG);
      jecUnc.push_back(jecUncG);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 278802, 280385 ));

      //IOV: 2016H
      std::vector<JetCorrectorParameters> correctionParametersH = std::vector<JetCorrectorParameters> ();
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersH.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersH = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorH = new FactorizedJetCorrector(correctionParametersH);
      std::string jecUncPathH = jecPathname+"/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncH = new JetCorrectionUncertainty(jecUncPathH);
      SimpleJetResolution* JetResolutionCalculatorH = new SimpleJetResolution(*JetResolutionParametersH);

      correctionParameters.push_back(correctionParametersH);
      JetResolutionParameters.push_back(JetResolutionParametersH);
      JetCorrector.push_back( JetCorrectorH );
      JetResolutionCalculator.push_back(JetResolutionCalculatorH);
      jecUnc.push_back(jecUncH);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 280919, 99999999 ));

    }
    else if (isFastsim) {
      std::cout << "Fastsim JEC\n";

      std::vector<JetCorrectorParameters> correctionParametersFastsim = std::vector<JetCorrectorParameters> ();
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersFastsim = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorFastsim = new FactorizedJetCorrector(correctionParametersFastsim);
      std::string jecUncPath = jecPathname+"/Spring16_FastSimV1_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncFastsim = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorFastsim = new SimpleJetResolution(*JetResolutionParametersFastsim);

      correctionParameters.push_back(correctionParametersFastsim);
      JetResolutionParameters.push_back(JetResolutionParametersFastsim);
      JetCorrector.push_back( JetCorrectorFastsim );
      JetResolutionCalculator.push_back(JetResolutionCalculatorFastsim);
      jecUnc.push_back(jecUncFastsim);
      JetCorrectionsIOV.push_back( std::pair<int,int>( -1, 99999999 ));
    }
    else {
      std::cout << "Loading Jet Energy Corrections: Summer16_23Sep2016V4_MC \n";
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));

      JetCorrectorParameters *JetResolutionParametersMC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorMC = new SimpleJetResolution(*JetResolutionParametersMC);

      correctionParameters.push_back(correctionParametersMC);
      JetResolutionParameters.push_back(JetResolutionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorMC);
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( -1, 99999999 ));
    }
}
void RazorHelper::loadAK8JetTag_Razor2016_MoriondRereco() {
    cout << "RazorHelper: loading 2016 top/W-tagging TF1s and histograms" << endl;
    puppiSoftDropCorrFile = TFile::Open("puppiCorr.root");
    puppiSoftDropCorr_Gen = (TF1*)puppiSoftDropCorrFile->Get("puppiJECcorr_gen");
    puppiSoftDropCorr_RecoCentral = (TF1*)puppiSoftDropCorrFile->Get("puppiJECcorr_reco_0eta1v3");
    puppiSoftDropCorr_RecoForward = (TF1*)puppiSoftDropCorrFile->Get("puppiJECcorr_reco_1v3eta2v5");

    wTopTagEffFile = TFile::Open("AK8WTopTagEff.root");
    wTagEffFullsim = (TH1F*)wTopTagEffFile->Get("WTagEffFullsim");
    wTagEffFastsim = (TH1F*)wTopTagEffFile->Get("WTagEffFastsim");
    wTagEffFastsimSF = (TH1F*)wTopTagEffFile->Get("WTagEffFastsimSF");

    topTagEffFullsim = (TH1F*)wTopTagEffFile->Get("TopTagEffFullsim");
    topTagEffFastsim = (TH1F*)wTopTagEffFile->Get("TopTagEffFastsim");
    topTagEffFastsimSF = (TH1F*)wTopTagEffFile->Get("TopTagEffFastsimSF");
}


////////////////////////////////////////////////
//  2016 Prompt Reco
////////////////////////////////////////////////

void RazorHelper::loadPhoton_Razor2016(){
    // photon efficiency scale factors
    std::cout << "RazorHelper: loading photon efficiency scale factor histograms" << std::endl;
    phoEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/efficiency_results_PhoLooseEffDenominatorReco_2016_ICHEP.root");
    phoLooseEffSFHist = (TH2D*)phoEffSFFile->Get("ScaleFactor_PhoLooseEffDenominatorReco");

    phoEffFastsimSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/PhotonEffFastsimToFullsimCorrectionFactors.2016.root");
    phoLooseEffFastsimSFHist = (TH2D*)phoEffFastsimSFFile->Get("ElectronLoose_FastsimScaleFactor");
}


void RazorHelper::loadBTag_Razor2016() {
    // b-tag efficiencies and scale factors
    std::cout << "RazorHelper: loading btag efficiency histograms" << std::endl;
    btagEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_BJets_25ns_CSVM_Fullsim_80X.root");
    btagCharmEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_CJets_25ns_CSVM_Fullsim_80X.root");
    btagLightJetsEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_LightJets_25ns_CSVM_Fullsim_80X.root");
    btagMediumEfficiencyHist = (TH2D*)btagEfficiencyFile->Get("Efficiency_PtEta");
    btagMediumCharmEfficiencyHist = (TH2D*)btagCharmEfficiencyFile->Get("Efficiency_PtEta");
    btagMediumLightJetsEfficiencyHist = (TH2D*)btagLightJetsEfficiencyFile->Get("Efficiency_PtEta");

    // Fullsim
    btagcalib = new BTagCalibration("csvv2", "./CSVv2_Moriond17_B_H.csv");
    btagreader = new BTagCalibrationReader(btagcalib,               // calibration instance
                                           BTagEntry::OP_MEDIUM,     // operating point
				           "mujets",                 // measurement type
				           "central");               // systematics type
    btagreader_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets", "up");  // sys up
    btagreader_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets", "down");  // sys down
    btagreaderMistag = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "central");
    btagreaderMistag_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "up");    // sys up
    btagreaderMistag_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "down");  // sys down

    // Fastsim
    btagcalibfastsim = new BTagCalibration("csvv2", "./fastsim_csvv2_ttbar_26_1_2017.csv");
    btagreaderfastsim = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "central");
    btagreaderfastsim_up = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "up");
    btagreaderfastsim_do = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "down");
}

void RazorHelper::loadJECs_Razor2016() {
    std::cout << "RazorHelper: loading jet energy correction constants, using Spring16_25nsV6." << std::endl;
    // initialize
    std::string jecPathname = "./";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();

    if (isData) {
      //IOV: 2016BCD
      std::vector<JetCorrectorParameters> correctionParametersBCD = std::vector<JetCorrectorParameters> ();
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10BCD_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10BCD_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10BCD_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersBCD.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10BCD_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));

      JetCorrectorParameters *JetResolutionParametersBCD = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorBCD = new FactorizedJetCorrector(correctionParametersBCD);
      std::string jecUncPathBCD = jecPathname+"/Spring16_PromptReco_V10_DATA/Spring16_25nsV10BCD_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncBCD = new JetCorrectionUncertainty(jecUncPathBCD);
      SimpleJetResolution* JetResolutionCalculatorBCD = new SimpleJetResolution(*JetResolutionParametersBCD);

      correctionParameters.push_back(correctionParametersBCD);
      JetResolutionParameters.push_back(JetResolutionParametersBCD);
      JetCorrector.push_back( JetCorrectorBCD );
      JetResolutionCalculator.push_back(JetResolutionCalculatorBCD);
      jecUnc.push_back(jecUncBCD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 276811 ));

      //IOV: 2016E
      std::vector<JetCorrectorParameters> correctionParametersE = std::vector<JetCorrectorParameters> ();
      correctionParametersE.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10E_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersE.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10E_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersE.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10E_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersE.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10E_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersE = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorE = new FactorizedJetCorrector(correctionParametersE);
      std::string jecUncPathE = jecPathname+"/Spring16_PromptReco_V10_DATA/Spring16_25nsV10E_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncE = new JetCorrectionUncertainty(jecUncPathE);
      SimpleJetResolution* JetResolutionCalculatorE = new SimpleJetResolution(*JetResolutionParametersE);

      correctionParameters.push_back(correctionParametersE);
      JetResolutionParameters.push_back(JetResolutionParametersE);
      JetCorrector.push_back( JetCorrectorE );
      JetResolutionCalculator.push_back(JetResolutionCalculatorE);
      jecUnc.push_back(jecUncE);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 276831, 277420 ));

      //IOV: 2016F
      std::vector<JetCorrectorParameters> correctionParametersF = std::vector<JetCorrectorParameters> ();
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10F_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10F_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10F_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10F_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersF = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorF = new FactorizedJetCorrector(correctionParametersF);
      std::string jecUncPathF = jecPathname+"/Spring16_PromptReco_V10_DATA/Spring16_25nsV10F_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncF = new JetCorrectionUncertainty(jecUncPathF);
      SimpleJetResolution* JetResolutionCalculatorF = new SimpleJetResolution(*JetResolutionParametersF);

      correctionParameters.push_back(correctionParametersF);
      JetResolutionParameters.push_back(JetResolutionParametersF);
      JetCorrector.push_back( JetCorrectorF );
      JetResolutionCalculator.push_back(JetResolutionCalculatorF);
      jecUnc.push_back(jecUncF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 277772, 276801 ));

      //IOV: 2016GH
      std::vector<JetCorrectorParameters> correctionParametersGH = std::vector<JetCorrectorParameters> ();
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10p2_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10p2_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10p2_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersGH.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_PromptReco_V10_DATA/Spring16_25nsV10p2_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersGH = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorGH = new FactorizedJetCorrector(correctionParametersGH);
      std::string jecUncPathGH = jecPathname+"/Spring16_PromptReco_V10_DATA/Spring16_25nsV10p2_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncGH = new JetCorrectionUncertainty(jecUncPathGH);
      SimpleJetResolution* JetResolutionCalculatorGH = new SimpleJetResolution(*JetResolutionParametersGH);

      correctionParameters.push_back(correctionParametersGH);
      JetResolutionParameters.push_back(JetResolutionParametersGH);
      JetCorrector.push_back( JetCorrectorGH );
      JetResolutionCalculator.push_back(JetResolutionCalculatorGH);
      jecUnc.push_back(jecUncGH);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 278802, 99999999 ));

    }
    else if (isFastsim) {
       std::vector<JetCorrectorParameters> correctionParametersFastsim = std::vector<JetCorrectorParameters> ();
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersFastsim = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorFastsim = new FactorizedJetCorrector(correctionParametersFastsim);
      std::string jecUncPath = jecPathname+"/Spring16_FastSimV1_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncFastsim = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorFastsim = new SimpleJetResolution(*JetResolutionParametersFastsim);

      correctionParameters.push_back(correctionParametersFastsim);
      JetResolutionParameters.push_back(JetResolutionParametersFastsim);
      JetCorrector.push_back( JetCorrectorFastsim );
      JetResolutionCalculator.push_back(JetResolutionCalculatorFastsim);
      jecUnc.push_back(jecUncFastsim);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 99999999 ));
    }
    else {
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersMC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorMC = new SimpleJetResolution(*JetResolutionParametersMC);

      correctionParameters.push_back(correctionParametersMC);
      JetResolutionParameters.push_back(JetResolutionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorMC);
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 99999999 ));
    }

}

////////////////////////////////////////////////
//  2016 ICHEP
////////////////////////////////////////////////

void RazorHelper::loadTag_Razor2016_ICHEP_80X() {
    loadPileup_Razor2016_ICHEP();
    loadLepton_Razor2016_ICHEP();
    loadPhoton_Razor2016(); // same as 2016 inclusive
    loadBTag_Razor2016(); // same as 2016 inclusive
    loadTrigger_Razor2016_ICHEP();
    loadJECs_Razor2016(); // same as 2016 inclusive
}

void RazorHelper::loadPileup_Razor2016_ICHEP() {
    // pileup weights
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 ICHEP pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PileupWeights/PileupReweight2016_ICHEP.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
}

void RazorHelper::loadLepton_Razor2016_ICHEP(){

    // electron efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 ICHEP electron efficiency histograms" << std::endl;
    eleTightEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Tight_Fullsim.root");
    eleEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightElectronSelectionEffDenominatorReco_2016_ICHEP_Golden.root");
    eleGSFTrackEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden_12p9_ICHEP/efficiencySF_GsfTracking_2016ICHEP_Golden.root");
    eleGSFTrackEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Reco_Fullsim.root");

    eleTightEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("Efficiency_PtEta");
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorReco");
    // We did not produce veto electron efficiencies or scale factors for this run period
    eleVetoEfficiencyFile = 0;
    eleVetoEfficiencyHist = 0;
    vetoEleEffSFFile = 0;
    eleVetoEffSFHist = 0;
    // We don't have any fastsim scale factors for this run period
    eleTightEffFastsimSFHist = 0;
    eleVetoEffFastsimSFHist = 0;
    // We don't apply any tracking scale factors for this (post-HIP-effect) run period
    eleGSFTrackEffSFHist = (TH2D*)eleGSFTrackEffSFFile->Get("EGamma_SF2D");
    eleGSFTrackEffHist = (TH2D*)eleGSFTrackEffFile->Get("Efficiency_PtEta");

    // muon efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 ICHEP muon efficiency histograms" << std::endl;
    muTightEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Tight_Fullsim.root");
    muEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightMuonSelectionEffDenominatorReco_2016_ICHEP_Golden.root");
    muTrackEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden_12p9_ICHEP/efficiencySF_TrackReconstruction_2016ICHEP.root");
    muTrackEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Reco_Fullsim.root");

    muTightEfficiencyHist = (TH2D*)muTightEfficiencyFile->Get("Efficiency_PtEta");
    muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorReco");
    // We did not produce veto electron efficiencies or scale factors for this run period
    muVetoEfficiencyFile = 0;
    muVetoEfficiencyHist = 0;
    vetoMuEffSFFile = 0;
    muVetoEffSFHist = 0;
    // We don't have any fastsim scale factors for this run period
    muTightEffFastsimSFHist = 0;
    muVetoEffFastsimSFHist = 0;
    // We don't apply any tracking scale factors for this (post-HIP-effect) run period
    muTrackEffSFHist = (TH2D*)muTrackEffSFFile->Get("muon");
    muTrackEffHist = (TH2D*)muTrackEffFile->Get("Efficiency_PtEta");

    // tau efficiencies and scale factors
    std::cout << "RazorHelper: loading tau efficiency histograms" << std::endl;
    tauEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/TauEffFastsimToFullsimCorrectionFactors.root");
    tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");

}

void RazorHelper::loadTrigger_Razor2016_ICHEP() {
    // single lepton trigger scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016 ICHEP trigger efficiency histograms" << std::endl;
    eleTrigSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2016_ICHEP_Golden.root");
    eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");
    muTrigSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2016_ICHEP_Golden.root");
    muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");

    eleTrigEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleElectronTriggerEfficiency_2016_ICHEP_Golden.root");
    eleTrigEffHist = (TH2D*)eleTrigEffFile->Get("hEffEtaPt");

    muTrigEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleMuonTriggerEfficiency_2016_ICHEP_Golden.root");
    muTrigEffHist = (TH2D*)muTrigEffFile->Get("hEffEtaPt");

    //diphoton trigger scale factors
    diphotonTrigLeadingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTLeadingLegEffDenominatorLoose_2016_ICHEP.root");
    diphotonTrigLeadingLegEffHist = (TH2D*)diphotonTrigLeadingLegEffFile->Get("hEffEtaPt");
    diphotonTrigTrailingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTTrailingLegEffDenominatorLoose_2016_ICHEP.root");
    diphotonTrigTrailingLegEffHist = (TH2D*)diphotonTrigTrailingLegEffFile->Get("hEffEtaPt");

    //get trigger numbers
    dileptonTriggerNums = { 44,45,57,59,64,65,66,67,68 };
    singleLeptonTriggerNums = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43 };
    hadronicTriggerNums = { 164,165,166,167,168,169,170,171,172,173,174,175,176 };
}



////////////////////////////////////////////////
//  2016 G
////////////////////////////////////////////////

void RazorHelper::loadTag_Razor2016G_80X() {
    loadPileup_Razor2016G();
    loadLepton_Razor2016G();
    loadPhoton_Razor2016(); // same as 2016 inclusive
    loadBTag_Razor2016(); // same as 2016 inclusive
    loadTrigger_Razor2016G();
    loadJECs_Razor2016(); // same as 2016 inclusive
}

void RazorHelper::loadPileup_Razor2016G() {
    // pileup weights
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016G pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PileupWeights/PileupReweight2016G.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
}

void RazorHelper::loadLepton_Razor2016G(){

    // electron efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016G electron efficiency histograms" << std::endl;
    eleTightEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Tight_Fullsim.root");
    eleVetoEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Veto_Fullsim.root");
    eleEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightElectronSelectionEffDenominatorReco_2016G_Golden.root");
    vetoEleEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoElectronSelectionEffDenominatorReco_2016G_Golden.root");

    eleTightEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("Efficiency_PtEta");
    eleVetoEfficiencyHist = (TH2D*)eleVetoEfficiencyFile->Get("Efficiency_PtEta");
    // We don't have any fastsim scale factors for this run period
    eleTightEffFastsimSFHist = 0;
    eleVetoEffFastsimSFHist = 0;
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorReco");
    eleVetoEffSFHist = (TH2D*)vetoEleEffSFFile->Get("ScaleFactor_VetoElectronSelectionEffDenominatorReco");
    // We don't apply any tracking scale factors for this (post-HIP-effect) run period
    eleGSFTrackEffSFHist = 0;
    eleGSFTrackEffHist = 0;

    // muon efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016G muon efficiency histograms" << std::endl;
    muTightEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Tight_Fullsim.root");
    muVetoEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Veto_Fullsim.root");
    muEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightMuonSelectionEffDenominatorReco_2016G_Golden.root");
    vetoMuEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoMuonSelectionEffDenominatorReco_2016G_Golden.root");
    muTightEfficiencyHist = (TH2D*)muTightEfficiencyFile->Get("Efficiency_PtEta");
    muVetoEfficiencyHist = (TH2D*)muVetoEfficiencyFile->Get("Efficiency_PtEta");
    // We don't have any fastsim scale factors for this run period
    muTightEffFastsimSFHist = 0;
    muVetoEffFastsimSFHist = 0;
    muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorReco");
    muVetoEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorReco");
    // We don't apply any tracking scale factors for this (post-HIP-effect) run period
    muTrackEffSFHist = 0;
    muTrackEffHist = 0;

    // tau efficiencies and scale factors
    std::cout << "RazorHelper: loading tau efficiency histograms" << std::endl;
    tauEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/TauEffFastsimToFullsimCorrectionFactors.root");
    tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");

}

void RazorHelper::loadTrigger_Razor2016G() {
    // single lepton trigger scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016G trigger efficiency histograms" << std::endl;
    eleTrigSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2016G_Golden.root");
    eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");
    muTrigSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2016G_Golden.root");
    muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");

    eleTrigEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleElectronTriggerEfficiency_2016G_Golden.root");
    eleTrigEffHist = (TH2D*)eleTrigEffFile->Get("hEffEtaPt");

    muTrigEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleMuonTriggerEfficiency_2016G_Golden.root");
    muTrigEffHist = (TH2D*)muTrigEffFile->Get("hEffEtaPt");

    //diphoton trigger scale factors
    diphotonTrigLeadingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTLeadingLegEffDenominatorLoose_2016_ICHEP.root");
    diphotonTrigLeadingLegEffHist = (TH2D*)diphotonTrigLeadingLegEffFile->Get("hEffEtaPt");
    diphotonTrigTrailingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTTrailingLegEffDenominatorLoose_2016_ICHEP.root");
    diphotonTrigTrailingLegEffHist = (TH2D*)diphotonTrigTrailingLegEffFile->Get("hEffEtaPt");

    //get trigger numbers
    dileptonTriggerNums = { 44,45,57,59,64,65,66,67,68 };
    singleLeptonTriggerNums = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43 };
    hadronicTriggerNums = { 164,165,166,167,168,169,170,171,172,173,174,175,176 };
}

////////////////////////////////////////////////
//  2016 G Unblinded Data
////////////////////////////////////////////////

void RazorHelper::loadTag_Razor2016G_SUSYUnblind_80X() {
    loadPileup_Razor2016G_SUSYUnblind();
    loadLepton_Razor2016G_SUSYUnblind();
    loadPhoton_Razor2016_MoriondRereco(); // same as 2016 inclusive
    loadBTag_Razor2016G_SUSYUnblind();
    loadTrigger_Razor2016G_SUSYUnblind();
    loadJECs_Razor2016_MoriondRereco(); // same as 2016 inclusive
}

void RazorHelper::loadPileup_Razor2016G_SUSYUnblind() {
    // pileup weights
    // LAST UPDATED: 3 November 2016
    std::cout << "RazorHelper: loading 2016G_SUSYUnblind pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PileupWeights/PileupReweight2016G_SUSYUnblind.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
}

void RazorHelper::loadLepton_Razor2016G_SUSYUnblind(){

    // electron efficiencies and scale factors
    // LAST UPDATED: 3 November 2016
    std::cout << "RazorHelper: loading 2016G_SUSYUnblind electron efficiency histograms" << std::endl;
    eleTightEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Tight_Fullsim.root");
    eleVetoEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Veto_Fullsim.root");
    eleEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightElectronSelectionEffDenominatorGen_2016G_Rereco_SUSYUnblind_Golden.root");
    vetoEleEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoElectronSelectionEffDenominatorGen_2016G_Rereco_SUSYUnblind_Golden.root");

    eleTightEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("Efficiency_PtEta");
    eleVetoEfficiencyHist = (TH2D*)eleVetoEfficiencyFile->Get("Efficiency_PtEta");
    // We don't have any fastsim scale factors for this run period
    eleTightEffFastsimSFHist = 0;
    eleVetoEffFastsimSFHist = 0;
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorGen");
    eleVetoEffSFHist = (TH2D*)vetoEleEffSFFile->Get("ScaleFactor_VetoElectronSelectionEffDenominatorGen");
    // We don't apply any tracking scale factors for this (post-HIP-effect) run period
    eleGSFTrackEffSFHist = 0;
    eleGSFTrackEffHist = 0;

    // muon efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016G_SUSYUnblind muon efficiency histograms" << std::endl;
    muTightEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Tight_Fullsim.root");
    muVetoEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptMuon_TTJets_25ns_Veto_Fullsim.root");
    muEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_TightMuonSelectionEffDenominatorGen_2016G_Rereco_SUSYUnblind_Golden.root");
    vetoMuEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_VetoMuonSelectionEffDenominatorGen_2016G_SUSYUnblind_Golden.root");
    muTightEfficiencyHist = (TH2D*)muTightEfficiencyFile->Get("Efficiency_PtEta");
    muVetoEfficiencyHist = (TH2D*)muVetoEfficiencyFile->Get("Efficiency_PtEta");
    // We don't have any fastsim scale factors for this run period
    muTightEffFastsimSFHist = 0;
    muVetoEffFastsimSFHist = 0;
    muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorGen");
    muVetoEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorGen");
    // We don't apply any tracking scale factors for this (post-HIP-effect) run period
    muTrackEffSFHist = 0;
    muTrackEffHist = 0;

    // tau efficiencies and scale factors
    std::cout << "RazorHelper: loading tau efficiency histograms" << std::endl;
    tauEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/FastsimToFullsim/TauEffFastsimToFullsimCorrectionFactors.root");
    tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");

}

void RazorHelper::loadBTag_Razor2016G_SUSYUnblind() {
    // b-tag efficiencies and scale factors
    std::cout << "RazorHelper: loading btag efficiency histograms" << std::endl;
    btagEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_BJets_25ns_CSVM_Fullsim_80X.root");
    btagCharmEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_CJets_25ns_CSVM_Fullsim_80X.root");
    btagLightJetsEfficiencyFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/Efficiency_LightJets_25ns_CSVM_Fullsim_80X.root");
    btagMediumEfficiencyHist = (TH2D*)btagEfficiencyFile->Get("Efficiency_PtEta");
    btagMediumCharmEfficiencyHist = (TH2D*)btagCharmEfficiencyFile->Get("Efficiency_PtEta");
    btagMediumLightJetsEfficiencyHist = (TH2D*)btagLightJetsEfficiencyFile->Get("Efficiency_PtEta");

    // Fullsim
    btagcalib = new BTagCalibration("csvv2", "./CSVv2_Moriond17_G_H.csv");
    btagreader = new BTagCalibrationReader(btagcalib,               // calibration instance
                                           BTagEntry::OP_MEDIUM,     // operating point
				           "comb",                 // measurement type
				           "central");               // systematics type
    btagreader_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets", "up");  // sys up
    btagreader_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets", "down");  // sys down
    btagreaderMistag = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "incl", "central");
    btagreaderMistag_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "incl", "up");    // sys up
    btagreaderMistag_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "incl", "down");  // sys down

    // Fastsim
    btagcalibfastsim = new BTagCalibration("csvv2", "./fastsim_csvv2_ttbar_26_1_2017.csv");
    btagreaderfastsim = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "central");
    btagreaderfastsim_up = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "up");
    btagreaderfastsim_do = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "down");
}

void RazorHelper::loadTrigger_Razor2016G_SUSYUnblind() {
    // single lepton trigger scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2016G_SUSYUnblind trigger efficiency histograms" << std::endl;
    eleTrigSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2016G_Rereco_SUSYUnblind_Golden.root");
    eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");
    muTrigSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2016G_Rereco_SUSYUnblind_Golden.root");
    muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");

    eleTrigEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleElectronTriggerEfficiency_2016G_Rereco_SUSYUnblind_Golden.root");
    eleTrigEffHist = (TH2D*)eleTrigEffFile->Get("hEffEtaPt");

    muTrigEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleMuonTriggerEfficiency_2016G_Rereco_SUSYUnblind_Golden.root");
    muTrigEffHist = (TH2D*)muTrigEffFile->Get("hEffEtaPt");

    //diphoton trigger scale factors
    diphotonTrigLeadingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTLeadingLegEffDenominatorLoose_2016_ICHEP.root");
    diphotonTrigLeadingLegEffHist = (TH2D*)diphotonTrigLeadingLegEffFile->Get("hEffEtaPt");
    diphotonTrigTrailingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTTrailingLegEffDenominatorLoose_2016_ICHEP.root");
    diphotonTrigTrailingLegEffHist = (TH2D*)diphotonTrigTrailingLegEffFile->Get("hEffEtaPt");

    //get trigger numbers
    dileptonTriggerNums = { 44,45,57,59,64,65,66,67,68 };
    singleLeptonTriggerNums = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43 };
    hadronicTriggerNums = { 164,165,166,167,168,169,170,171,172,173,174,175,176 };
}


////////////////////////////////////////////////
//  2017 PromptReco
////////////////////////////////////////////////
void RazorHelper::loadTag_Razor2017_92X() {
  loadPileup_Razor2017_92X();
  loadLepton_Razor2016_MoriondRereco();
  loadPhoton_Razor2017_92X();
  loadBTag_Razor2016_MoriondRereco();
  loadTrigger_Razor2017_92X();
  loadJECs_Razor2017_17Nov2017Rereco();
}

void RazorHelper::loadPileup_Razor2017_92X() {
    // pileup weights
    // LAST UPDATED: 28 October 2017
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;

    if (!isFastsim) {
      pileupWeightFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/PileupWeights/PileupReweight_2017_41p2ifb.root");
      pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
      pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
      pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
      std::cout << "PileupReweight_2017_41p2ifb.root\n";
    } else {
      // Will do something for Fastsim in the future
      pileupWeightFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/PileupWeights/PileupReweight_2017_41p2ifb.root");
      pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
      pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
      pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
      std::cout << "PileupReweight_2017_41p2ifb.root\n";
    }


}

void RazorHelper::loadTrigger_Razor2017_92X() {
    // single lepton trigger scale factors
    // LAST UPDATED: 30 July 2017
    std::cout << "RazorHelper: loading 2016 trigger efficiency histograms" << std::endl;
    eleTrigSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2016_Rereco_Golden.root");
    eleTrigSFHist = (TH2D*)eleTrigSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");

    muTrigSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2016_Rereco_Golden.root");
    muTrigSFHist = (TH2D*)muTrigSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");

    eleTrigEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleElectronTriggerEfficiency_2016_Rereco_Golden.root");
    eleTrigEffHist = (TH2D*)eleTrigEffFile->Get("hEffEtaPt");

    muTrigEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/SingleMuonTriggerEfficiency_2016_Rereco_Golden.root");
    muTrigEffHist = (TH2D*)muTrigEffFile->Get("hEffEtaPt");

    //diphoton trigger scale factors
    diphotonTrigLeadingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTLeadingLegEffDenominatorLoose_2016_Rereco.root");
    diphotonTrigLeadingLegEffHist = (TH2D*)diphotonTrigLeadingLegEffFile->Get("hEffEtaPt");
    diphotonTrigTrailingLegEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/PhoHLTTrailingLegEffDenominatorLoose_2016_Rereco.root");
    diphotonTrigTrailingLegEffHist = (TH2D*)diphotonTrigTrailingLegEffFile->Get("hEffEtaPt");

    diphotonTrigLeadingLegEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/efficiency_results_PhoHLTLeadingLegEffDenominatorLoose_2016_Rereco.root");
    diphotonTrigLeadingLegEffSFHist = (TH2D*)diphotonTrigLeadingLegEffSFFile->Get("ScaleFactor_PhoHLTLeadingLegEffDenominatorLoose");
    diphotonTrigTrailingLegEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/efficiency_results_PhoHLTTrailingLegEffDenominatorLoose_2016_Rereco.root");
    diphotonTrigTrailingLegEffSFHist = (TH2D*)diphotonTrigTrailingLegEffSFFile->Get("ScaleFactor_PhoHLTTrailingLegEffDenominatorLoose");

    //get trigger index numbers
    std::cout << "RazorHelper: loading 2017 trigger indices" << std::endl;
    dileptonTriggerNums = { 16, 17,18, 19, 20, 21, 22, 23,24 };
    singleLeptonTriggerNums = { 1,2,3,4,5,6,7,12,13,14,15 };
    hadronicTriggerNums = { 106, 107, 108, 109, 110, 111 };
}


void RazorHelper::loadPhoton_Razor2017_92X(){
    // photon efficiency scale factors
    // use avaerage results for run 2017BCDEF for now
    std::cout << "RazorHelper: loading photon efficiency scale factor histograms" << std::endl;
    phoEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2017/efficiency_results_PhoLooseEffDenominatorReco_2017BCDEF_94X.root");
    phoLooseEffSFHist = (TH2D*)phoEffSFFile->Get("EGamma_SF2D");

    // results for 2017MC is not available yet, use 2016 version for now
    phoEffFastsimSFFile = TFile::Open("PhotonEffFastsimToFullsimCorrectionFactors.2016.root");
    phoLooseEffFastsimSFHist = (TH2D*)phoEffFastsimSFFile->Get("ElectronLoose_FastsimScaleFactor");
}

////////////////////////////////////////////////
//  2017 17Nov2017 Rereco
////////////////////////////////////////////////
void RazorHelper::loadTag_Razor2017_17Nov2017Rereco() {
  loadPileup_Razor2017_17Nov2017Rereco();
  loadLepton_Razor2017_17Nov2017Rereco();
  // loadPhoton_Razor2017_92X();
  // loadBTag_Razor2017_17Nov2017Rereco();
  // loadTrigger_Razor2017_92X();
  loadTrigger_Razor2017_17Nov2017Rereco();
  loadJECs_Razor2017_17Nov2017Rereco();
  loadHiggsPt();
}
void RazorHelper::loadTag_Razor2017_Source2018() {
  loadPileup_Razor2017_Source2018();
  loadLepton_Razor2017_17Nov2017Rereco();
  // loadPhoton_Razor2017_92X();
  // loadBTag_Razor2017_17Nov2017Rereco();
  // loadTrigger_Razor2017_92X();
  loadTrigger_Razor2017_17Nov2017Rereco();
  loadJECs_Razor2017_17Nov2017Rereco();
  loadHiggsPt();
}

void RazorHelper::loadPileup_Razor2017_Source2018() {
    // pileup weights
    // LAST UPDATED: 06 July 2018
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;


      // Will do something for Fastsim in the future
    pileupWeightFile = TFile::Open("PileupReweight_source18_target17.root.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
    std::cout << "PileupReweight_source18_target17.root.root\n";

}
void RazorHelper::loadPileup_Razor2017_17Nov2017Rereco() {
    // pileup weights
    // LAST UPDATED: 06 July 2018
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;


      // Will do something for Fastsim in the future
    pileupWeightFile = TFile::Open("PileupReweight_MC_Fall17_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
    std::cout << "PileupReweight_MC_Fall17_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root\n";

}
void RazorHelper::loadTrigger_Razor2017_17Nov2017Rereco() {
    // single lepton trigger scale factors
    std::cout << "RazorHelper: loading trigger efficiency histograms" << std::endl;
    metTriggerSFFile = TFile::Open("METTriggers_SF.root");
    metTriggerSFHist = (TH1F*)metTriggerSFFile->Get("trigger_efficiency_Fall17");

}
void RazorHelper::loadLepton_Razor2017_17Nov2017Rereco(){

    // LAST UPDATED: 18 August 2020

    std::cout << "RazorHelper: loading 2017 muon efficiency histograms" << std::endl;


    muTriggerEfficiencyFile = TFile::Open("MuonEfficienciesAndSF_RunBtoF_Nov17Nov2017.root");
    muIdEfficiencyFile = TFile::Open("MuonRunBCDEF_mc_ID.2017.root");
    muIsoEfficiencyFile = TFile::Open("MuonRunBCDEF_mc_ISO.2017.root");
    muIdSFFile = TFile::Open("MuonRunBCDEF_SF_ID.2017.root");
    muIsoSFFile = TFile::Open("MuonRunBCDEF_SF_ISO.2017.root");


    muTriggerEfficiencyDir = muTriggerEfficiencyFile ->GetDirectory("IsoMu27_PtEtaBins");
    muTriggerSFHist =  (TH2D*)muTriggerEfficiencyDir->Get("pt_abseta_ratio");//pt on xaxis and eta on yaxis
    muTriggerEfficiencyHist =  (TH2D*)muTriggerEfficiencyDir->GetDirectory("efficienciesMC")->Get("pt_abseta_MC");//pt on xaxis and eta on yaxis


    muLooseIdEfficiencyHist = (TH2D*)muIdEfficiencyFile->Get("NUM_LooseID_DEN_genTracks_pt_abseta");
    muTightIdEfficiencyHist = (TH2D*)muIdEfficiencyFile->Get("NUM_TightID_DEN_genTracks_pt_abseta");
    muLooseIsoEfficiencyHist = (TH2D*)muIsoEfficiencyFile->Get("NUM_LooseRelIso_DEN_LooseID_pt_abseta");
    muTightIsoEfficiencyHist = (TH2D*)muIsoEfficiencyFile->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

    muLooseIdSFHist = (TH2D*)muIdSFFile->Get("NUM_LooseID_DEN_genTracks_pt_abseta");
    muTightIdSFHist = (TH2D*)muIdSFFile->Get("NUM_TightID_DEN_genTracks_pt_abseta");
    muLooseIsoSFHist = (TH2D*)muIsoSFFile->Get("NUM_LooseRelIso_DEN_LooseID_pt_abseta");
    muTightIsoSFHist = (TH2D*)muIsoSFFile->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");




}

void RazorHelper::loadBTag_Razor2017_17Nov2017Rereco() {
    // b-tag efficiencies and scale factors
    std::cout << "RazorHelper: loading btag efficiency histograms for tag 17Nov2017Rereco" << std::endl;
    btagEfficiencyFile = TFile::Open("Efficiency_BJets_25ns_CSVM_Fullsim_80X.root");
    btagCharmEfficiencyFile = TFile::Open("Efficiency_CJets_25ns_CSVM_Fullsim_80X.root");
    btagLightJetsEfficiencyFile = TFile::Open("Efficiency_LightJets_25ns_CSVM_Fullsim_80X.root");
    btagMediumEfficiencyHist = (TH2D*)btagEfficiencyFile->Get("Efficiency_PtEta");
    btagMediumCharmEfficiencyHist = (TH2D*)btagCharmEfficiencyFile->Get("Efficiency_PtEta");
    btagMediumLightJetsEfficiencyHist = (TH2D*)btagLightJetsEfficiencyFile->Get("Efficiency_PtEta");

    // Fullsim
   btagcalib = new BTagCalibration("csvv2", "./CSVv2_94XSF_V2_B_F.csv");
   btagreader = new BTagCalibrationReader( btagcalib,               // calibration instance
                                           BTagEntry::OP_MEDIUM,     // operating point
				           "comb",                 // measurement type
				           "central");               // systematics type
    btagreader_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "up");  // sys up
    btagreader_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "comb", "down");  // sys down
    btagreaderMistag = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "incl", "central");
    btagreaderMistag_up = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "incl", "up");    // sys up
    btagreaderMistag_do = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "incl", "down");  // sys down

    // Fastsim
    btagcalibfastsim = new BTagCalibration("csvv2", "./csvv2_13TEV_17_6_3_2019.csv");
    btagreaderfastsim = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "central");
    btagreaderfastsim_up = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "up");
    btagreaderfastsim_do = new BTagCalibrationReader(btagcalibfastsim, BTagEntry::OP_MEDIUM, "fastsim", "down");

}

void RazorHelper::loadJECs_Razor2017_17Nov2017Rereco() {
    std::cout << "RazorHelper: loading jet energy correction constants, using Fall17_17Nov2017_V32." << std::endl;
    // initialize
    std::string jecPathname = "JEC/";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    std::cout << "here1\n";
    if (isData) {
      //IOV: 2017B
      std::vector<JetCorrectorParameters> correctionParametersB = std::vector<JetCorrectorParameters> ();
      correctionParametersB.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017B_V32_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersB.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017B_V32_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersB.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017B_V32_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersB.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersB = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorB = new FactorizedJetCorrector(correctionParametersB);
      std::string jecUncPathB = jecPathname+"/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017B_V32_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncB = new JetCorrectionUncertainty(jecUncPathB);
      SimpleJetResolution* JetResolutionCalculatorB = new SimpleJetResolution(*JetResolutionParametersB);

      correctionParameters.push_back(correctionParametersB);
      JetResolutionParameters.push_back(JetResolutionParametersB);
      JetCorrector.push_back( JetCorrectorB );
      JetResolutionCalculator.push_back(JetResolutionCalculatorB);
      jecUnc.push_back(jecUncB);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 297046, 299329 ));

      //IOV: 2017C
      std::vector<JetCorrectorParameters> correctionParametersC = std::vector<JetCorrectorParameters> ();
      correctionParametersC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017C_V32_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017C_V32_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017C_V32_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017C_V32_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorC = new FactorizedJetCorrector(correctionParametersC);
      std::string jecUncPathC = jecPathname+"/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017C_V32_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncC = new JetCorrectionUncertainty(jecUncPathC);
      SimpleJetResolution* JetResolutionCalculatorC = new SimpleJetResolution(*JetResolutionParametersC);

      correctionParameters.push_back(correctionParametersC);
      JetResolutionParameters.push_back(JetResolutionParametersC);
      JetCorrector.push_back( JetCorrectorC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorC);
      jecUnc.push_back(jecUncC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 299368, 302029 ));

      //IOV: 2017DE
      std::vector<JetCorrectorParameters> correctionParametersD = std::vector<JetCorrectorParameters> ();
      correctionParametersD.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017DE_V32_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersD.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017DE_V32_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersD.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017DE_V32_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersD.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017DE_V32_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersD = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorD = new FactorizedJetCorrector(correctionParametersD);
      std::string jecUncPathD = jecPathname+"/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017DE_V32_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncD = new JetCorrectionUncertainty(jecUncPathD);
      SimpleJetResolution* JetResolutionCalculatorD = new SimpleJetResolution(*JetResolutionParametersD);

      correctionParameters.push_back(correctionParametersD);
      JetResolutionParameters.push_back(JetResolutionParametersD);
      JetCorrector.push_back( JetCorrectorD );
      JetResolutionCalculator.push_back(JetResolutionCalculatorD);
      jecUnc.push_back(jecUncD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 302030, 304797 ));



      //IOV: 2017F
      std::vector<JetCorrectorParameters> correctionParametersF = std::vector<JetCorrectorParameters> ();
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017F_V32_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017F_V32_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017F_V32_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017F_V32_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersF = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorF = new FactorizedJetCorrector(correctionParametersF);
      std::string jecUncPathF = jecPathname+"/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017F_V32_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncF = new JetCorrectionUncertainty(jecUncPathF);
      SimpleJetResolution* JetResolutionCalculatorF = new SimpleJetResolution(*JetResolutionParametersF);

      correctionParameters.push_back(correctionParametersF);
      JetResolutionParameters.push_back(JetResolutionParametersF);
      JetCorrector.push_back( JetCorrectorF );
      JetResolutionCalculator.push_back(JetResolutionCalculatorF);
      jecUnc.push_back(jecUncF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 305040, 306462 ));


    }
    else if (isFastsim) {
      std::cout << "Fastsim JEC\n";

      std::vector<JetCorrectorParameters> correctionParametersFastsim = std::vector<JetCorrectorParameters> ();
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Spring16_FastSimV1_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersFastsim = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorFastsim = new FactorizedJetCorrector(correctionParametersFastsim);
      std::string jecUncPath = jecPathname+"/Spring16_FastSimV1_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncFastsim = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorFastsim = new SimpleJetResolution(*JetResolutionParametersFastsim);

      correctionParameters.push_back(correctionParametersFastsim);
      JetResolutionParameters.push_back(JetResolutionParametersFastsim);
      JetCorrector.push_back( JetCorrectorFastsim );
      JetResolutionCalculator.push_back(JetResolutionCalculatorFastsim);
      jecUnc.push_back(jecUncFastsim);
      JetCorrectionsIOV.push_back( std::pair<int,int>( -1, 99999999 ));
    }
    else {
      std::cout << "Loading Jet Energy Corrections: Fall17_17Nov2017V32_MC \n";
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));

      JetCorrectorParameters *JetResolutionParametersMC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorMC = new SimpleJetResolution(*JetResolutionParametersMC);

      correctionParameters.push_back(correctionParametersMC);
      JetResolutionParameters.push_back(JetResolutionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorMC);
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( -1, 99999999 ));
    }
}


////////////////////////////////////////////////
//  2017 31Mar2018 Rereco
////////////////////////////////////////////////
void RazorHelper::loadTag_Razor2017_31Mar2018Rereco() {
  loadPileup_Razor2017_17Nov2017Rereco();
  loadLepton_Razor2017_31Mar2018Rereco();
  loadPhoton_Razor2017_31Mar2018Rereco();
  loadBTag_Razor2017_17Nov2017Rereco();
  loadTrigger_Razor2017_92X();
  loadJECs_Razor2017_31Mar2018Rereco();
}

void RazorHelper::loadPhoton_Razor2017_31Mar2018Rereco(){
//identical to loadPhoton_Razor2017_92X, would check if there's new version released
    // photon efficiency scale factors
    // use avaerage results for run 2017BCDEF for now
    std::cout << "RazorHelper: loading photon efficiency scale factor histograms" << std::endl;
    phoEffSFFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2017/efficiency_results_PhoLooseEffDenominatorReco_2017BCDEF_94X.root");
    phoLooseEffSFHist = (TH2D*)phoEffSFFile->Get("EGamma_SF2D");

    // results for 2017MC is not available yet, use 2016 version for now
    phoEffFastsimSFFile = TFile::Open("PhotonEffFastsimToFullsimCorrectionFactors.2016.root");
    phoLooseEffFastsimSFHist = (TH2D*)phoEffFastsimSFFile->Get("ElectronLoose_FastsimScaleFactor");
}

void RazorHelper::loadLepton_Razor2017_31Mar2018Rereco(){

    // electron efficiencies and scale factors
    // LAST UPDATED: 31 August 2018
    std::cout << "RazorHelper: loading 2017 electron efficiency histograms" << std::endl;
    //eleTightEfficiencyFile = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    //eleLooseEfficiencyFile = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    //eleVetoEfficiencyFile = TFile::Open("ElectronMVAIDScaleFactor_SUSYVLoose_2017_17Nov2017Rereco.root");
    eleGSFTrackEffFile = TFile::Open("root://eoscms:///store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2016_Golden/Efficiency_PromptElectron_TTJets_25ns_Reco_Fullsim.root");
    eleEffSFFile = TFile::Open("ElectronScaleFactors_Run2017_31Mar2018.root");
    looseEleEffSFFile = TFile::Open("efficiency_results_LooseElectronSelectionEffDenominatorGen_2017_31Mar2018_Golden.root");
    //vetoEleEffSFFile = TFile::Open("efficiency_results_VetoElectronSelectionEffDenominatorGen_2017_31Mar2018_Golden.root");
    eleGSFTrackEffSFFile = TFile::Open("ElectronRecoEffScaleFactors_Run2017.root");
    eleTightEffFastsimSFFile = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.2016.root");
    eleLooseEffFastsimSFFile = TFile::Open("detailed_ele_full_fast_sf_17.root");
    eleVetoEffFastsimSFFile = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.2016.root");

    // eleTightEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("ElectronEff_Tight_Fullsim");
    // eleLooseEfficiencyHist = (TH2D*)eleTightEfficiencyFile->Get("ElectronEff_Loose_Fullsim");
    // eleVetoEfficiencyHist = (TH2D*)eleVetoEfficiencyFile->Get("ElectronEff_Veto_Fullsim");
    // eleGSFTrackEffHist = (TH2D*)eleGSFTrackEffFile->Get("Efficiency_PtEta");
    // We don't have ID scale factors for Fastsim yet.
    eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorGen");
    //for now 20190124
    eleLooseEffSFHist = (TH2D*)looseEleEffSFFile->Get("ScaleFactor_LooseElectronSelectionEffDenominatorGen");
    //eleLooseEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_LooseElectronSelectionEffDenominatorGen");
    //eleVetoEffSFHist = (TH2D*)vetoEleEffSFFile->Get("ScaleFactor_VetoElectronSelectionEffDenominatorGen");
    eleGSFTrackEffSFHist = (TH2D*)eleGSFTrackEffSFFile->Get("h2_scaleFactorsEGamma");
    eleTightEffFastsimSFHist =  (TH2D*)eleTightEffFastsimSFFile->Get("ElectronTight_FastsimScaleFactor");
    //for now 20190124
    eleLooseEffFastsimSFHist =  (TH2D*)eleLooseEffFastsimSFFile->Get("CutBasedLooseNoIso94XV2_sf");
    //eleLooseEffFastsimSFHist =  (TH2D*)eleLooseEffFastsimSFFile->Get("ElectronLoose_FastsimScaleFactor");
    eleVetoEffFastsimSFHist = (TH2D*)eleVetoEffFastsimSFFile->Get("ElectronEff_Veto_Fullsim");

    // muon efficiencies and scale factors
    // LAST UPDATED: 18 October 2016
    std::cout << "RazorHelper: loading 2017 muon efficiency histograms" << std::endl;
    muTightEfficiencyFile = TFile::Open("MuonIsoScaleFactor_2017_17Nov2017Rereco.root");
    muVetoEfficiencyFile = TFile::Open("MuonIsoScaleFactor_2017_17Nov2017Rereco.root");
    //muEffSFFile = TFile::Open("efficiency_results_TightMuonSelectionEffDenominatorGen_2017_31Mar2018_Golden.root");
    vetoMuEffSFFile = TFile::Open("efficiency_results_VetoMuonSelectionEffDenominatorGen_2017_31Mar2018_Golden.root");
    muTrackEffSFFile = TFile::Open("efficiencySF_muEleTracking_2016_average.root");
    muTrackEffFile = TFile::Open("Efficiency_PromptMuon_TTJets_25ns_Reco_Fullsim.root");
    muTightEffFastsimSFFile = TFile::Open("MuonEffFastsimToFullsimCorrectionFactors.2016.root");
    muVetoEffFastsimSFFile = TFile::Open("MuonEffFastsimToFullsimCorrectionFactors.2016.root");

    muTightEfficiencyHist = (TH2D*)muTightEfficiencyFile->Get("MuonEff_Tight_Fullsim");
    muVetoEfficiencyHist = (TH2D*)muVetoEfficiencyFile->Get("MuonEff_Veto_Fullsim");
    // We don't have ID scale factors for Fastsim yet.
    muTightEffFastsimSFHist = (TH2D*)muTightEffFastsimSFFile->Get("MuonTight_FastsimScaleFactor");
    muVetoEffFastsimSFHist = (TH2D*)muVetoEffFastsimSFFile->Get("MuonVeto_FastsimScaleFactor");
    //muTightEffSFHist = (TH2D*)muEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorGen");
    //for now 20190124
    muLooseEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorGen");
    muLooseEffFastsimSFHist = (TH2D*)muVetoEffFastsimSFFile->Get("MuonVeto_FastsimScaleFactor");
    muVetoEffSFHist = (TH2D*)vetoMuEffSFFile->Get("ScaleFactor_VetoMuonSelectionEffDenominatorGen");
    muTrackEffHist = (TH2D*)muTrackEffFile->Get("Efficiency_PtEta");
    // These scale factors are weighted according to the fraction of the 2016 run affected
    // by the 'HIP' issue, under the assumption that tracking scale factors are 1 for runs
    // not affected by the 'HIP'.
    muTrackEffSFHist = (TH2D*)muTrackEffSFFile->Get("muon");

    // tau efficiencies and scale factors
    std::cout << "RazorHelper: loading tau efficiency histograms" << std::endl;
    tauEfficiencyFile = TFile::Open("TauEffFastsimToFullsimCorrectionFactors.2016.root");
    tauLooseEfficiencyHist = (TH2D*)tauEfficiencyFile->Get("TauEff_Loose_Fullsim");

}


void RazorHelper::loadJECs_Razor2017_31Mar2018Rereco() {
    std::cout << "RazorHelper: loading jet energy correction constants, using Fall17_17Nov2017_V32." << std::endl;
    // initialize
    std::string jecPathname = "./";
    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    std::cout << "here1\n";
    if (isData) {
      //IOV: 2017B
      std::vector<JetCorrectorParameters> correctionParametersB = std::vector<JetCorrectorParameters> ();
      correctionParametersB.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017B_V32_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersB.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017B_V32_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersB.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017B_V32_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersB.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersB = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorB = new FactorizedJetCorrector(correctionParametersB);
      std::string jecUncPathB = jecPathname+"/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017B_V32_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncB = new JetCorrectionUncertainty(jecUncPathB);
      SimpleJetResolution* JetResolutionCalculatorB = new SimpleJetResolution(*JetResolutionParametersB);

      correctionParameters.push_back(correctionParametersB);
      JetResolutionParameters.push_back(JetResolutionParametersB);
      JetCorrector.push_back( JetCorrectorB );
      JetResolutionCalculator.push_back(JetResolutionCalculatorB);
      jecUnc.push_back(jecUncB);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 1, 299329 ));

      //IOV: 2017C
      std::vector<JetCorrectorParameters> correctionParametersC = std::vector<JetCorrectorParameters> ();
      correctionParametersC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017C_V32_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017C_V32_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017C_V32_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017C_V32_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorC = new FactorizedJetCorrector(correctionParametersC);
      std::string jecUncPathC = jecPathname+"/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017C_V32_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncC = new JetCorrectionUncertainty(jecUncPathC);
      SimpleJetResolution* JetResolutionCalculatorC = new SimpleJetResolution(*JetResolutionParametersC);

      correctionParameters.push_back(correctionParametersC);
      JetResolutionParameters.push_back(JetResolutionParametersC);
      JetCorrector.push_back( JetCorrectorC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorC);
      jecUnc.push_back(jecUncC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 299368, 302029 ));

      //IOV: 2017DE
      std::vector<JetCorrectorParameters> correctionParametersDE = std::vector<JetCorrectorParameters> ();
      correctionParametersDE.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017DE_V32_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersDE.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017DE_V32_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersDE.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017DE_V32_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersDE.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017DE_V32_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersDE = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorDE = new FactorizedJetCorrector(correctionParametersDE);
      std::string jecUncPathDE = jecPathname+"/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017DE_V32_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncDE = new JetCorrectionUncertainty(jecUncPathDE);
      SimpleJetResolution* JetResolutionCalculatorDE = new SimpleJetResolution(*JetResolutionParametersDE);

      correctionParameters.push_back(correctionParametersDE);
      JetResolutionParameters.push_back(JetResolutionParametersDE);
      JetCorrector.push_back( JetCorrectorDE );
      JetResolutionCalculator.push_back(JetResolutionCalculatorDE);
      jecUnc.push_back(jecUncDE);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 302030, 304797 ));

      //IOV: 2017F
      std::vector<JetCorrectorParameters> correctionParametersF = std::vector<JetCorrectorParameters> ();
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017F_V32_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017F_V32_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017F_V32_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersF.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017F_V32_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersF = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorF = new FactorizedJetCorrector(correctionParametersF);
      std::string jecUncPathF = jecPathname+"/Fall17_17Nov2017_V32_DATA/Fall17_17Nov2017F_V32_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncF = new JetCorrectionUncertainty(jecUncPathF);
      SimpleJetResolution* JetResolutionCalculatorF = new SimpleJetResolution(*JetResolutionParametersF);

      correctionParameters.push_back(correctionParametersF);
      JetResolutionParameters.push_back(JetResolutionParametersF);
      JetCorrector.push_back( JetCorrectorF );
      JetResolutionCalculator.push_back(JetResolutionCalculatorF);
      jecUnc.push_back(jecUncF);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 305040, 99999999 ));


    }
    else if (isFastsim) {
      std::cout << "Fastsim JEC\n";

      std::vector<JetCorrectorParameters> correctionParametersFastsim = std::vector<JetCorrectorParameters> ();
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_FastsimV1_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_FastsimV1_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersFastsim.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_FastsimV1_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersFastsim = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorFastsim = new FactorizedJetCorrector(correctionParametersFastsim);
      std::string jecUncPath = jecPathname+"/Fall17_FastsimV1_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncFastsim = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorFastsim = new SimpleJetResolution(*JetResolutionParametersFastsim);

      correctionParameters.push_back(correctionParametersFastsim);
      JetResolutionParameters.push_back(JetResolutionParametersFastsim);
      JetCorrector.push_back( JetCorrectorFastsim );
      JetResolutionCalculator.push_back(JetResolutionCalculatorFastsim);
      jecUnc.push_back(jecUncFastsim);
      JetCorrectionsIOV.push_back( std::pair<int,int>( -1, 99999999 ));
    }
    else {
      std::cout << "Loading Jet Energy Corrections: Fall17_17Nov2017V32_MC \n";
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));

      JetCorrectorParameters *JetResolutionParametersMC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorMC = new SimpleJetResolution(*JetResolutionParametersMC);

      correctionParameters.push_back(correctionParametersMC);
      JetResolutionParameters.push_back(JetResolutionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorMC);
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( -1, 99999999 ));
    }
}


////////////////////////////////////////////////
//  2018
////////////////////////////////////////////////




void RazorHelper::loadTag_Razor2018_17SeptEarlyReReco() {
  loadPileup_Razor2018_17SeptEarlyReReco();
  // loadLepton_Razor2016_MoriondRereco();
  // loadPhoton_Razor2017_92X();
  // loadBTag_Razor2016_MoriondRereco();
  loadTrigger_Razor2018_17SeptEarlyReReco();
  loadJECs_Razor2018_17SeptEarlyReReco();
  loadHiggsPt();
}

void RazorHelper::loadPileup_Razor2018_17SeptEarlyReReco() {
    // pileup weights
    // LAST UPDATED: 9 Dec 2019
    std::cout << "RazorHelper: loading pileup weight histograms" << std::endl;
    pileupWeightFile = TFile::Open("PileupReweight_MC_Fall18_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root");
    // pileupWeightFile = TFile::Open("PileupReweight_source16_target18.root");

    // pileupWeightFile = TFile::Open("/storage/user/christiw/login-1/christiw/LLP/CMSSW_9_4_4/src/llp_analyzer/data/PileupWeights/PileupReweight_MC_ggH_HToSSTobbbb_ms55_pl1000_RunIIFall18.root");
    pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
    pileupWeightSysUpHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysUp");
    pileupWeightSysDownHist = (TH1F*)pileupWeightFile->Get("PileupReweightSysDown");
    std::cout << "PileupReweight_MC_ggH_HToSSTobbbb_ms55_pl1000_RunIIFall18.root\n";
}
void RazorHelper::loadTrigger_Razor2018_17SeptEarlyReReco() {
    // single lepton trigger scale factors
    std::cout << "RazorHelper: loading trigger efficiency histograms" << std::endl;
    metTriggerSFFile = TFile::Open("METTriggers_SF.root");
    metTriggerSFHist = (TH1F*)metTriggerSFFile->Get("trigger_efficiency_Fall18");

}


void RazorHelper::loadJECs_Razor2018_17SeptEarlyReReco() {
    std::cout << "RazorHelper: loading jet energy correction constants, using 2018_17SeptEarlyReReco" << std::endl;
    // initialize
    std::string jecPathname = "JEC/";
    //
    // std::string jecPathname = getenv("CMSSW_BASE");
    // if(jecPathname != NULL) jecPathname = jecPathname + "/src/llp_analyzer/data/JEC/";
    // else cout<<"CMSSW ENV NOT SET"<<endl;
    // if(jecPathname != NULL and option == 1) jecPathname = "JEC/"; //run on condor if option == 1


    correctionParameters = std::vector<std::vector<JetCorrectorParameters> >();
    JetResolutionParameters = std::vector<JetCorrectorParameters*>();
    JetCorrector = std::vector<FactorizedJetCorrector*>();
    jecUnc = std::vector<JetCorrectionUncertainty*>();
    JetResolutionCalculator = std::vector<SimpleJetResolution*>();
    JetCorrectionsIOV = std::vector<std::pair<int,int> >();
    if (isData) {
      //IOV: 2018A
      std::vector<JetCorrectorParameters> correctionParametersA = std::vector<JetCorrectorParameters> ();
      correctionParametersA.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunA_V19_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersA.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunA_V19_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersA.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunA_V19_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersA.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunA_V19_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersA = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorA = new FactorizedJetCorrector(correctionParametersA);
      std::string jecUncPathA = jecPathname+"/Autumn18_RunsABCD_V19_DATA/Autumn18_RunA_V19_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncA = new JetCorrectionUncertainty(jecUncPathA);
      SimpleJetResolution* JetResolutionCalculatorA = new SimpleJetResolution(*JetResolutionParametersA);

      correctionParameters.push_back(correctionParametersA);
      JetResolutionParameters.push_back(JetResolutionParametersA);
      JetCorrector.push_back( JetCorrectorA );
      JetResolutionCalculator.push_back(JetResolutionCalculatorA);
      jecUnc.push_back(jecUncA);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 315252, 316995 ));

      //IOV: 2018B
      std::vector<JetCorrectorParameters> correctionParametersB = std::vector<JetCorrectorParameters> ();
      correctionParametersB.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunB_V19_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersB.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunB_V19_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersB.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunB_V19_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersB.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunB_V19_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersB = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorB = new FactorizedJetCorrector(correctionParametersB);
      std::string jecUncPathB = jecPathname+"/Autumn18_RunsABCD_V19_DATA/Autumn18_RunB_V19_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncB = new JetCorrectionUncertainty(jecUncPathB);
      SimpleJetResolution* JetResolutionCalculatorB = new SimpleJetResolution(*JetResolutionParametersB);

      correctionParameters.push_back(correctionParametersB);
      JetResolutionParameters.push_back(JetResolutionParametersB);
      JetCorrector.push_back( JetCorrectorB );
      JetResolutionCalculator.push_back(JetResolutionCalculatorB);
      jecUnc.push_back(jecUncB);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 317080, 319310 ));

      //IOV: 2018C
      std::vector<JetCorrectorParameters> correctionParametersC = std::vector<JetCorrectorParameters> ();
      correctionParametersC.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunC_V19_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersC.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunC_V19_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersC.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunC_V19_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersC.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunC_V19_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorC = new FactorizedJetCorrector(correctionParametersC);
      std::string jecUncPathC = jecPathname+"/Autumn18_RunsABCD_V19_DATA/Autumn18_RunC_V19_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncC = new JetCorrectionUncertainty(jecUncPathC);
      SimpleJetResolution* JetResolutionCalculatorC = new SimpleJetResolution(*JetResolutionParametersC);

      correctionParameters.push_back(correctionParametersC);
      JetResolutionParameters.push_back(JetResolutionParametersC);
      JetCorrector.push_back( JetCorrectorC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorC);
      jecUnc.push_back(jecUncC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 319337, 320065 ));

      //IOV: 2018D
      std::vector<JetCorrectorParameters> correctionParametersD = std::vector<JetCorrectorParameters> ();
      correctionParametersD.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunD_V19_DATA_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersD.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunD_V19_DATA_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersD.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunD_V19_DATA_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersD.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_RunsABCD_V19_DATA/Autumn18_RunD_V19_DATA_L2L3Residual_AK4PFchs.txt", jecPathname.c_str())));
      JetCorrectorParameters *JetResolutionParametersD = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorD = new FactorizedJetCorrector(correctionParametersD);
      std::string jecUncPathD = jecPathname+"/Autumn18_RunsABCD_V19_DATA/Autumn18_RunD_V19_DATA_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncD = new JetCorrectionUncertainty(jecUncPathD);
      SimpleJetResolution* JetResolutionCalculatorD = new SimpleJetResolution(*JetResolutionParametersD);

      correctionParameters.push_back(correctionParametersD);
      JetResolutionParameters.push_back(JetResolutionParametersD);
      JetCorrector.push_back( JetCorrectorD );
      JetResolutionCalculator.push_back(JetResolutionCalculatorD);
      jecUnc.push_back(jecUncD);
      JetCorrectionsIOV.push_back( std::pair<int,int>( 320673, 325175 ));

    }

    else {
      std::cout << "Loading Jet Energy Corrections: Autumn18_V19_MC \n";
      std::vector<JetCorrectorParameters> correctionParametersMC = std::vector<JetCorrectorParameters> ();
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_V19_MC/Autumn18_V19_MC_L1FastJet_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_V19_MC/Autumn18_V19_MC_L2Relative_AK4PFchs.txt", jecPathname.c_str())));
      correctionParametersMC.push_back(JetCorrectorParameters(
                  Form("%s/Autumn18_V19_MC/Autumn18_V19_MC_L3Absolute_AK4PFchs.txt", jecPathname.c_str())));

      JetCorrectorParameters *JetResolutionParametersMC = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",jecPathname.c_str()));
      FactorizedJetCorrector *JetCorrectorMC = new FactorizedJetCorrector(correctionParametersMC);
      std::string jecUncPath = jecPathname+"/Autumn18_V19_MC/Autumn18_V19_MC_Uncertainty_AK4PFchs.txt";
      JetCorrectionUncertainty *jecUncMC = new JetCorrectionUncertainty(jecUncPath);
      SimpleJetResolution* JetResolutionCalculatorMC = new SimpleJetResolution(*JetResolutionParametersMC);

      correctionParameters.push_back(correctionParametersMC);
      JetResolutionParameters.push_back(JetResolutionParametersMC);
      JetCorrector.push_back( JetCorrectorMC );
      JetResolutionCalculator.push_back(JetResolutionCalculatorMC);
      jecUnc.push_back(jecUncMC);
      JetCorrectionsIOV.push_back( std::pair<int,int>( -1, 99999999 ));
    }
}



////////////////////////////////////////////////
//  Utilities
////////////////////////////////////////////////




enum TheRunEra{y2016B,y2016C,y2016D,y2016E,y2016F,y2016G,y2016H,y2017B,y2017C,y2017D,y2017E,y2017F,y2018A,y2018B,y2018C,y2018D,y2016MC,y2017MC,y2018MC};

std::pair<double,double> RazorHelper::METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int runnb, int year, bool isMC, int npv){

  std::pair<double,double>  TheXYCorr_Met_MetPhi(uncormet,uncormet_phi);

  if(npv>100) npv=100;
  int runera =-1;
  bool usemetv2 =false;
  if(isMC && year == 2016) runera = y2016MC;
  else if(isMC && year == 2017) {runera = y2017MC; usemetv2 =true;}
  else if(isMC && year == 2018) runera = y2018MC;

  else if(!isMC && runnb >=272007 &&runnb<=275376  ) runera = y2016B;
  else if(!isMC && runnb >=275657 &&runnb<=276283  ) runera = y2016C;
  else if(!isMC && runnb >=276315 &&runnb<=276811  ) runera = y2016D;
  else if(!isMC && runnb >=276831 &&runnb<=277420  ) runera = y2016E;
  else if(!isMC && runnb >=277772 &&runnb<=278808  ) runera = y2016F;
  else if(!isMC && runnb >=278820 &&runnb<=280385  ) runera = y2016G;
  else if(!isMC && runnb >=280919 &&runnb<=284044  ) runera = y2016H;

  else if(!isMC && runnb >=297020 &&runnb<=299329 ){ runera = y2017B; usemetv2 =true;}
  else if(!isMC && runnb >=299337 &&runnb<=302029 ){ runera = y2017C; usemetv2 =true;}
  else if(!isMC && runnb >=302030 &&runnb<=303434 ){ runera = y2017D; usemetv2 =true;}
  else if(!isMC && runnb >=303435 &&runnb<=304826 ){ runera = y2017E; usemetv2 =true;}
  else if(!isMC && runnb >=304911 &&runnb<=306462 ){ runera = y2017F; usemetv2 =true;}

  else if(!isMC && runnb >=315252 &&runnb<=316995 ) runera = y2018A;
  else if(!isMC && runnb >=316998 &&runnb<=319312 ) runera = y2018B;
  else if(!isMC && runnb >=319313 &&runnb<=320393 ) runera = y2018C;
  else if(!isMC && runnb >=320394 &&runnb<=325273 ) runera = y2018D;

  else {
    //Couldn't find data/MC era => no correction applied
    return TheXYCorr_Met_MetPhi;
  }

  double METxcorr(0.),METycorr(0.);

  if(!usemetv2){//Current recommendation for 2016 and 2018
    if(runera==y2016B) METxcorr = -(-0.0478335*npv -0.108032);
    if(runera==y2016B) METycorr = -(0.125148*npv +0.355672);
    if(runera==y2016C) METxcorr = -(-0.0916985*npv +0.393247);
    if(runera==y2016C) METycorr = -(0.151445*npv +0.114491);
    if(runera==y2016D) METxcorr = -(-0.0581169*npv +0.567316);
    if(runera==y2016D) METycorr = -(0.147549*npv +0.403088);
    if(runera==y2016E) METxcorr = -(-0.065622*npv +0.536856);
    if(runera==y2016E) METycorr = -(0.188532*npv +0.495346);
    if(runera==y2016F) METxcorr = -(-0.0313322*npv +0.39866);
    if(runera==y2016F) METycorr = -(0.16081*npv +0.960177);
    if(runera==y2016G) METxcorr = -(0.040803*npv -0.290384);
    if(runera==y2016G) METycorr = -(0.0961935*npv +0.666096);
    if(runera==y2016H) METxcorr = -(0.0330868*npv -0.209534);
    if(runera==y2016H) METycorr = -(0.141513*npv +0.816732);
    if(runera==y2017B) METxcorr = -(-0.259456*npv +1.95372);
    if(runera==y2017B) METycorr = -(0.353928*npv -2.46685);
    if(runera==y2017C) METxcorr = -(-0.232763*npv +1.08318);
    if(runera==y2017C) METycorr = -(0.257719*npv -1.1745);
    if(runera==y2017D) METxcorr = -(-0.238067*npv +1.80541);
    if(runera==y2017D) METycorr = -(0.235989*npv -1.44354);
    if(runera==y2017E) METxcorr = -(-0.212352*npv +1.851);
    if(runera==y2017E) METycorr = -(0.157759*npv -0.478139);
    if(runera==y2017F) METxcorr = -(-0.232733*npv +2.24134);
    if(runera==y2017F) METycorr = -(0.213341*npv +0.684588);
    if(runera==y2018A) METxcorr = -(0.362865*npv -1.94505);
    if(runera==y2018A) METycorr = -(0.0709085*npv -0.307365);
    if(runera==y2018B) METxcorr = -(0.492083*npv -2.93552);
    if(runera==y2018B) METycorr = -(0.17874*npv -0.786844);
    if(runera==y2018C) METxcorr = -(0.521349*npv -1.44544);
    if(runera==y2018C) METycorr = -(0.118956*npv -1.96434);
    if(runera==y2018D) METxcorr = -(0.531151*npv -1.37568);
    if(runera==y2018D) METycorr = -(0.0884639*npv -1.57089);
    if(runera==y2016MC) METxcorr = -(-0.195191*npv -0.170948);
    if(runera==y2016MC) METycorr = -(-0.0311891*npv +0.787627);
    if(runera==y2017MC) METxcorr = -(-0.217714*npv +0.493361);
    if(runera==y2017MC) METycorr = -(0.177058*npv -0.336648);
    if(runera==y2018MC) METxcorr = -(0.296713*npv -0.141506);
    if(runera==y2018MC) METycorr = -(0.115685*npv +0.0128193);
  }
  else {//these are the corrections for v2 MET recipe (currently recommended for 2017)
    if(runera==y2016B) METxcorr = -(-0.0374977*npv +0.00488262);
    if(runera==y2016B) METycorr = -(0.107373*npv +-0.00732239);
    if(runera==y2016C) METxcorr = -(-0.0832562*npv +0.550742);
    if(runera==y2016C) METycorr = -(0.142469*npv +-0.153718);
    if(runera==y2016D) METxcorr = -(-0.0400931*npv +0.753734);
    if(runera==y2016D) METycorr = -(0.127154*npv +0.0175228);
    if(runera==y2016E) METxcorr = -(-0.0409231*npv +0.755128);
    if(runera==y2016E) METycorr = -(0.168407*npv +0.126755);
    if(runera==y2016F) METxcorr = -(-0.0161259*npv +0.516919);
    if(runera==y2016F) METycorr = -(0.141176*npv +0.544062);
    if(runera==y2016G) METxcorr = -(0.0583851*npv +-0.0987447);
    if(runera==y2016G) METycorr = -(0.0641427*npv +0.319112);
    if(runera==y2016H) METxcorr = -(0.0706267*npv +-0.13118);
    if(runera==y2016H) METycorr = -(0.127481*npv +0.370786);
    if(runera==y2017B) METxcorr = -(-0.19563*npv +1.51859);
    if(runera==y2017B) METycorr = -(0.306987*npv +-1.84713);
    if(runera==y2017C) METxcorr = -(-0.161661*npv +0.589933);
    if(runera==y2017C) METycorr = -(0.233569*npv +-0.995546);
    if(runera==y2017D) METxcorr = -(-0.180911*npv +1.23553);
    if(runera==y2017D) METycorr = -(0.240155*npv +-1.27449);
    if(runera==y2017E) METxcorr = -(-0.149494*npv +0.901305);
    if(runera==y2017E) METycorr = -(0.178212*npv +-0.535537);
    if(runera==y2017F) METxcorr = -(-0.165154*npv +1.02018);
    if(runera==y2017F) METycorr = -(0.253794*npv +0.75776);
    if(runera==y2018A) METxcorr = -(0.362642*npv +-1.55094);
    if(runera==y2018A) METycorr = -(0.0737842*npv +-0.677209);
    if(runera==y2018B) METxcorr = -(0.485614*npv +-2.45706);
    if(runera==y2018B) METycorr = -(0.181619*npv +-1.00636);
    if(runera==y2018C) METxcorr = -(0.503638*npv +-1.01281);
    if(runera==y2018C) METycorr = -(0.147811*npv +-1.48941);
    if(runera==y2018D) METxcorr = -(0.520265*npv +-1.20322);
    if(runera==y2018D) METycorr = -(0.143919*npv +-0.979328);
    if(runera==y2016MC) METxcorr = -(-0.159469*npv +-0.407022);
    if(runera==y2016MC) METycorr = -(-0.0405812*npv +0.570415);
    if(runera==y2017MC) METxcorr = -(-0.182569*npv +0.276542);
    if(runera==y2017MC) METycorr = -(0.155652*npv +-0.417633);
    if(runera==y2018MC) METxcorr = -(0.299448*npv +-0.13866);
    if(runera==y2018MC) METycorr = -(0.118785*npv +0.0889588);
  }

  double CorrectedMET_x = uncormet *cos( uncormet_phi)+METxcorr;
  double CorrectedMET_y = uncormet *sin( uncormet_phi)+METycorr;

  double CorrectedMET = sqrt(CorrectedMET_x*CorrectedMET_x+CorrectedMET_y*CorrectedMET_y);
  double CorrectedMETPhi;
  if(CorrectedMET_x==0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::Pi();
  else if(CorrectedMET_x==0 && CorrectedMET_y<0 )CorrectedMETPhi = -TMath::Pi();
  else if(CorrectedMET_x >0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x);
  else if(CorrectedMET_x <0&& CorrectedMET_y>0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) + TMath::Pi();
  else if(CorrectedMET_x <0&& CorrectedMET_y<0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) - TMath::Pi();
  else CorrectedMETPhi =0;

  TheXYCorr_Met_MetPhi.first= CorrectedMET;
  TheXYCorr_Met_MetPhi.second= CorrectedMETPhi;
  return TheXYCorr_Met_MetPhi;

}

double RazorHelper::getHiggsPtWeight(float higgsPt) {
    if (higgsPtWeightHist) {
        return higgsPtWeightHist->GetBinContent(higgsPtWeightHist->GetXaxis()->FindFixBin(higgsPt));
    }
    else {
        std::cout << "RazorHelper error: higgsPt weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}
double RazorHelper::getHiggsPtWeightSys(float higgsPt, int sys_i) {
    if (higgsPtWeightSysHist[sys_i]) {
        return higgsPtWeightSysHist[sys_i]->GetBinContent(higgsPtWeightSysHist[sys_i]->GetXaxis()->FindFixBin(higgsPt));
    }
    else {
        std::cout << "RazorHelper error: higgsPt weight systematics requested, but no histogram available!" << std::endl;
        return 0;
    }
}


double RazorHelper::getPileupWeight(int NPU) {
    if (pileupWeightHist) {
        return pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "RazorHelper error: pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

double RazorHelper::getPileupWeightUp(int NPU) {
    if (pileupWeightSysUpHist) {
        return pileupWeightSysUpHist->GetBinContent(pileupWeightSysUpHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "RazorHelper error: 'up' pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

double RazorHelper::getPileupWeightDown(int NPU) {
    if (pileupWeightSysDownHist) {
        return pileupWeightSysDownHist->GetBinContent(pileupWeightSysDownHist->GetXaxis()->FindFixBin(NPU));
    }
    else {
        std::cout << "RazorHelper error: 'down' pileup weight requested, but no histogram available!" << std::endl;
        return 0;
    }
}

double RazorHelper::getMetTriggerSF(float met) {
    if (metTriggerSFHist) {
        return metTriggerSFHist->GetBinContent(metTriggerSFHist->GetXaxis()->FindFixBin(met));
    }
    else {
        std::cout << "RazorHelper error: MET trigger SF requested, but no histogram available!" << std::endl;
        return 0;
    }
}


// Get scale factor from a histogram with pt on the x-axis and eta on the y-axis
double RazorHelper::lookupPtEtaScaleFactor(TH2D *hist, double pt, double eta, double ptmin, double ptmax, bool useAbsEta) {
    if (hist) {
      // constrain to histogram bounds
      ptmax = hist->GetXaxis()->GetXmax() * 0.999;
      double etaToUse = eta;
      if( useAbsEta ) {
	       etaToUse = fabs(eta);
      }

      return hist->GetBinContent(
				 hist->GetXaxis()->FindFixBin( fmax( fmin(pt,ptmax), ptmin ) ),
				 hist->GetYaxis()->FindFixBin( etaToUse )
				 );
    }
    else {
        std::cout << "Error: expected a histogram, got a null pointer" << std::endl;
        return 0;
    }
}

// Get scale factor uncertainty for a histogram with pt on the x-axis and eta on the y-axis
double RazorHelper::lookupPtEtaScaleFactorError(TH2D *hist, double pt, double eta, double ptmin, double ptmax, bool useAbsEta) {
    if (hist) {
        // constrain to histogram bounds
        if( ptmax > hist->GetXaxis()->GetXmax() * 0.999 ) {
            ptmax = hist->GetXaxis()->GetXmax() * 0.999;
        }
        double etaToUse = eta;
        if( useAbsEta ) {
            etaToUse = fabs(eta);
        }
        return hist->GetBinError(
                hist->GetXaxis()->FindFixBin( fmax( fmin(pt,ptmax), ptmin ) ),
                hist->GetYaxis()->FindFixBin( etaToUse )
                );
    }
    else {
        std::cout << "Error: expected a histogram, got a null pointer" << std::endl;
        return 0;
    }
}

// Get scale factor from histogram with eta on the x-axis and pt on the y-axis
double RazorHelper::lookupEtaPtScaleFactor(TH2D *hist, double pt, double eta, double ptmin, double ptmax, bool useAbsEta) {
    if (hist) {
        // constrain to histogram bounds
        if( ptmax > hist->GetYaxis()->GetXmax() * 0.999 ) {
            ptmax = hist->GetYaxis()->GetXmax() * 0.999;
        }
        double etaToUse = eta;
        if( useAbsEta ) {
            etaToUse = fabs(eta);
        }
        return hist->GetBinContent(
                hist->GetXaxis()->FindFixBin( etaToUse ),
                hist->GetYaxis()->FindFixBin( fmax( fmin(pt,ptmax), ptmin ) )
                );
    }
    else {
        std::cout << "Error: expected a histogram, got a null pointer" << std::endl;
        return 0;
    }
}


// Get scale factor error from histogram with eta on the x-axis and pt on the y-axis
double RazorHelper::lookupEtaPtScaleFactorError(TH2D *hist, double pt, double eta, double ptmin, double ptmax, bool useAbsEta) {
    if (hist) {
        // constrain to histogram bounds
        if( ptmax > hist->GetYaxis()->GetXmax() * 0.999 ) {
            ptmax = hist->GetYaxis()->GetXmax() * 0.999;
        }
        double etaToUse = eta;
        if( useAbsEta ) {
            etaToUse = fabs(eta);
        }
        return hist->GetBinError(
                hist->GetXaxis()->FindFixBin( etaToUse ),
                hist->GetYaxis()->FindFixBin( fmax( fmin(pt,ptmax), ptmin ) )
                );
    }
    else {
        std::cout << "Error: expected a histogram, got a null pointer" << std::endl;
        return 0;
    }
}


// Gets the correct event-level scale factor depending on whether the object passes or fails selection
double RazorHelper::getPassOrFailScaleFactor(double eff, double sf, bool passes) {
  if (passes) {
    return sf;
  }
  else if (eff*sf >= 1.0 || eff >= 1.0) { //provide safety against infinite or negative weights
    return 1.0;
  }
  else if (eff == 0 || sf == 0) { //provide safety against infinite or negative weights
    return 0.0;
  }
  return (1/eff - sf) / (1/eff - 1);
}

// Helper function for computing scale factors from histograms
// The smear parameter is used to add a fractional uncertainty in quadrature with the one from the histogram
// Returns a vector of five scale factors:
// 1) central value of scale factor
// 2,3) scale factor plus/minus uncertainty from fullsim scale factor
// 4,5) scale factor plus/minus uncertainty from fastsim scale factor
std::vector<double> RazorHelper::getLeptonScaleFactors(TH2D *effHist, TH2D *sfHist,
        TH2D *fastsimHist, double pt, double eta, bool passes, double smear) {

    double eff = lookupPtEtaScaleFactor( effHist, pt, eta );
    double eff_fastsimUp = eff;
    double eff_fastsimDown = eff;
    if (isFastsim) { //correct efficiency for Fastsim
        double sf = lookupPtEtaScaleFactor( fastsimHist, pt, eta );
        double sfErr = lookupPtEtaScaleFactorError( fastsimHist, pt, eta );
        eff *= sf;
        eff_fastsimUp *= (sf + sfErr);
        eff_fastsimDown *= (sf - sfErr);
    }
    double effSF = lookupPtEtaScaleFactor( sfHist, pt, eta );
    double effSFErr = lookupPtEtaScaleFactorError( sfHist, pt, eta );
    // add smearing factor to error
    effSFErr = sqrt(effSFErr*effSFErr + (smear*effSF)*(smear*effSF));
    double effSFUp = effSF + effSFErr;
    double effSFDown = effSF - effSFErr;

    double tmpSF = getPassOrFailScaleFactor( eff, effSF, passes );
    double tmpSFUp = getPassOrFailScaleFactor( eff, effSFUp, passes );
    double tmpSFDown = getPassOrFailScaleFactor( eff, effSFDown, passes );
    double tmpSF_fastsimUp = getPassOrFailScaleFactor( eff_fastsimUp, effSF, passes );
    double tmpSF_fastsimDown = getPassOrFailScaleFactor( eff_fastsimDown, effSF, passes );

    std::vector<double> out { tmpSF, tmpSFUp, tmpSFDown, tmpSF_fastsimUp, tmpSF_fastsimDown };
    return out;
}

// This function updates the provided scale factor values in-place by multiplying them by the
// SFs corresponding to the provided values of pt and eta.
void RazorHelper::updateScaleFactors(TH2D *effHist, TH2D *sfHist, TH2D *fastsimHist,
        float pt, float eta, bool passes, float &sf, float &sfUp, float &sfDown,
        float &sfFastsimUp, float &sfFastsimDown, float smear) {

    //retrieve scale factors and deal correctly with passing/failing selection
    std::vector<double> scaleFactors = getLeptonScaleFactors( effHist, sfHist, fastsimHist, pt, eta, passes, smear );

    //propagate to input values
    sf *= scaleFactors[0];
    sfUp *= scaleFactors[1]/scaleFactors[0];
    sfDown *= scaleFactors[2]/scaleFactors[0];
    sfFastsimUp *= scaleFactors[3]/scaleFactors[0];
    sfFastsimDown *= scaleFactors[4]/scaleFactors[0];
}

void RazorHelper::updateTightMuonScaleFactors(float pt, float eta, bool isTight,
        float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown) {

    updateScaleFactors( muTightEfficiencyHist, muTightEffSFHist, muTightEffFastsimSFHist,
            pt, eta, isTight, sf, sfUp, sfDown, sfFastsimUp, sfFastsimDown );

}

void RazorHelper::updateVetoMuonScaleFactors(float pt, float eta, bool isVeto,
        float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown) {
    updateScaleFactors( muVetoEfficiencyHist, muVetoEffSFHist, muVetoEffFastsimSFHist,
            fmax( pt, muVetoEffSFMinPt), eta, isVeto, sf, sfUp, sfDown, sfFastsimUp, sfFastsimDown );

}

void RazorHelper::updateTightElectronScaleFactors(float pt, float eta, bool isTight,
        float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown) {

    updateScaleFactors( eleTightEfficiencyHist, eleTightEffSFHist, eleTightEffFastsimSFHist,
            pt, eta, isTight, sf, sfUp, sfDown, sfFastsimUp, sfFastsimDown, 0.02 );

}

void RazorHelper::updateLooseElectronScaleFactors(float pt, float eta, bool isLoose,
        float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown) {

    updateScaleFactors( eleLooseEfficiencyHist, eleLooseEffSFHist, eleLooseEffFastsimSFHist,
            pt, eta, isLoose, sf, sfUp, sfDown, sfFastsimUp, sfFastsimDown, 0.02 );

}

void RazorHelper::updateVetoElectronScaleFactors(float pt, float eta, bool isVeto,
        float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown) {

    updateScaleFactors( eleVetoEfficiencyHist, eleVetoEffSFHist, eleVetoEffFastsimSFHist,
            fmax( pt, eleVetoEffSFMinPt), eta, isVeto, sf, sfUp, sfDown, sfFastsimUp, sfFastsimDown, 0.02 );

}

// Computes a single lepton scale factor
double RazorHelper::getLeptonScaleFactor(TH2D *effHist, TH2D *sfHist, TH2D *fastsimHist,
        double pt, double eta, bool passes) {
    double eff = lookupPtEtaScaleFactor( effHist, pt, eta );
    if (isFastsim) { //correct efficiency for Fastsim
        double sf = lookupPtEtaScaleFactor( fastsimHist, pt, eta );
        eff *= sf;
    }
    double effSF = lookupPtEtaScaleFactor( sfHist, pt, eta );

    return getPassOrFailScaleFactor( eff, effSF, passes );
}

double RazorHelper::getLeptonScaleFactor(TH2D *sfHist, double pt, double eta) {
    return lookupPtEtaScaleFactor( sfHist, pt, eta );
}

double RazorHelper::getLeptonScaleFactorError(TH2D *effHist, TH2D *sfHist, TH2D *fastsimHist,
       double pt, double eta, bool passes) {
    if( effHist == NULL ) return -999;
    if( !passes ) return -999;
    double unc = lookupPtEtaScaleFactorError( sfHist, pt, eta );
    if (isFastsim) { //correct efficiency for Fastsim
        unc = lookupPtEtaScaleFactorError( fastsimHist, pt, eta );
    }

    return unc;
}

double RazorHelper::getMuonScaleFactor(float pt, float eta, bool isTrigger, bool id, bool isTight) {

    if (isTrigger) return getLeptonScaleFactor( muTriggerSFHist, pt, eta );
    else if (id && !isTight) return getLeptonScaleFactor( muLooseIdSFHist, pt, eta );
    else if (id && isTight) return getLeptonScaleFactor( muTightIdSFHist, pt, eta );
    else if (!id && !isTight) return getLeptonScaleFactor( muLooseIsoSFHist, pt, eta );
    else return getLeptonScaleFactor( muTightIsoSFHist, pt, eta );
}


double RazorHelper::getMuonMCEfficiency(float pt, float eta, bool isTrigger, bool id, bool isTight) {
    if (isTrigger) return getLeptonScaleFactor( muTriggerEfficiencyHist, pt, eta );
    else if (id && !isTight) return getLeptonScaleFactor( muLooseIdEfficiencyHist, pt, eta );
    else if (id && isTight) return getLeptonScaleFactor( muTightIdEfficiencyHist, pt, eta );
    else if (!id && !isTight) return getLeptonScaleFactor( muLooseIsoEfficiencyHist, pt, eta );
    else return getLeptonScaleFactor( muTightIsoEfficiencyHist, pt, eta );
}


double RazorHelper::getTightMuonScaleFactor(float pt, float eta, bool isTight) {
    return getLeptonScaleFactor( muTightEfficiencyHist, muTightEffSFHist, muTightEffFastsimSFHist,
            pt, eta, isTight );
}

double RazorHelper::getVetoMuonScaleFactor(float pt, float eta, bool isVeto) {
    return getLeptonScaleFactor( muVetoEfficiencyHist, muVetoEffSFHist, muVetoEffFastsimSFHist,
            fmax( pt, muVetoEffSFMinPt), eta, isVeto );
}

double RazorHelper::getVetoMuonScaleFactorError(float pt, float eta, bool isVeto) {
    return getLeptonScaleFactorError( muVetoEfficiencyHist, muVetoEffSFHist, muVetoEffFastsimSFHist,
            fmax( pt, muVetoEffSFMinPt), eta, isVeto );
}

double RazorHelper::getLooseMuonScaleFactor(float pt, float eta, bool isLoose) {
    return getLeptonScaleFactor( muLooseEfficiencyHist, muLooseEffSFHist, muLooseEffFastsimSFHist,
            fmax( pt, muLooseEffSFMinPt), eta, isLoose );
}

double RazorHelper::getLooseMuonScaleFactorError(float pt, float eta, bool isLoose) {
    return getLeptonScaleFactorError( muLooseEfficiencyHist, muLooseEffSFHist, muLooseEffFastsimSFHist,
            fmax( pt, muLooseEffSFMinPt), eta, isLoose );
}

double RazorHelper::getMuonTrackScaleFactor(float pt, float eta, bool isReconstructed) {
    double eff = lookupPtEtaScaleFactor( muTrackEffHist, pt, eta, 10.1, 199.9, false ); //use signed eta
    double sf = lookupPtEtaScaleFactor( muTrackEffSFHist, pt, eta, 10.1, 199.9, false ); //use signed eta
    return getPassOrFailScaleFactor( eff, sf, isReconstructed );
}

double RazorHelper::getTightElectronScaleFactor(float pt, float eta, bool isTight) {
    return getLeptonScaleFactor( eleTightEfficiencyHist, eleTightEffSFHist, eleTightEffFastsimSFHist,
            pt, eta, isTight );
}

double RazorHelper::getLooseElectronScaleFactor(float pt, float eta, bool isLoose) {
    return getLeptonScaleFactor( eleLooseEfficiencyHist, eleLooseEffSFHist, eleLooseEffFastsimSFHist,
            pt, eta, isLoose );
}

double RazorHelper::getLooseElectronScaleFactorError(float pt, float eta, bool isLoose) {
    return getLeptonScaleFactorError( eleLooseEfficiencyHist, eleLooseEffSFHist, eleLooseEffFastsimSFHist,
            pt, eta, isLoose );
}

double RazorHelper::getVetoElectronScaleFactor(float pt, float eta, bool isVeto) {
    return getLeptonScaleFactor( eleVetoEfficiencyHist, eleVetoEffSFHist, eleVetoEffFastsimSFHist,
            fmax( pt, eleVetoEffSFMinPt), eta, isVeto );
}

double RazorHelper::getEleGSFTrackScaleFactor(float pt, float eta, bool isReconstructed) {
    double eff = lookupPtEtaScaleFactor( eleGSFTrackEffHist, pt, eta, 10.1, 199.9, false ); //use signed eta
    // note that the electron GSF tracking efficiency histogram has eta on the x-axis
    double sf = lookupEtaPtScaleFactor( eleGSFTrackEffSFHist, pt, eta, 10.1, 199.9, false ); //use signed eta
    return getPassOrFailScaleFactor( eff, sf, isReconstructed );
}


double RazorHelper::getPhotonScaleFactor(float pt, float eta, bool invert) {
  double sf = 1.0;
  if (phoLooseEffSFHist)
    {
      if( invert )
	{
	  sf = lookupEtaPtScaleFactor( phoLooseEffSFHist, pt, eta, 20.01, 99.9, false );
	}
      else
	{
	  sf = lookupPtEtaScaleFactor( phoLooseEffSFHist, pt, eta, 20.01, 99.9 );
	}
    }
  else { std::cout << "[WARNING] Could not load phoLooseEffSFHist.\n"; }
  return sf;
}

double RazorHelper::getPhotonScaleFactor_Tight(float pt, float eta, bool invert) {
  double sf = 1.0;
  if (phoTightEffSFHist)
    {
      if( invert )
        {
          sf = lookupEtaPtScaleFactor( phoTightEffSFHist, pt, eta, 20.01, 99.9, false );
        }
      else
        {
          sf = lookupPtEtaScaleFactor( phoTightEffSFHist, pt, eta, 20.01, 99.9 );
        }
    }
  else { std::cout << "[WARNING] Could not load phoTightEffSFHist.\n"; }
  return sf;
}

double RazorHelper::getPhotonScaleFactorError(float pt, float eta, bool invert) {
  double unc = 1.0;
  if (phoLooseEffSFHist)
    {
      if( invert )
	{
	  unc = lookupEtaPtScaleFactorError( phoLooseEffSFHist, pt, eta, 20.01, 99.9, false );
	}
      else
	{
	  unc = lookupPtEtaScaleFactorError( phoLooseEffSFHist, pt, eta, 20.01, 99.9 );
	}
    }
  else { std::cout << "[WARNING] Could not load phoLooseEffSFHist.\n"; }
  return unc;
}

double RazorHelper::getPhotonFastsimToFullsimScaleFactor(float pt, float eta) {
  double sf = 1.0;
  if (phoLooseEffFastsimSFHist) sf = lookupPtEtaScaleFactor( phoLooseEffFastsimSFHist, pt, eta, 20.01, 299.9 );
  else { std::cout << "[WARNING] Could not load phoLooseEffFastsimSFHist.\n"; }
  return sf;
}

double RazorHelper::getPhotonFastsimToFullsimScaleFactorError(float pt, float eta) {
  double unc = 1.0;
  if (phoLooseEffFastsimSFHist) unc = lookupPtEtaScaleFactorError( phoLooseEffFastsimSFHist, pt, eta, 20.01, 299.9 );
  else { std::cout << "[WARNING] Could not load phoLooseEffFastsimSFHist.\n"; }
  return unc;
}

// Returns the trigger scale factor corresponding to the given values of pt and eta.
// For fastsim events, the scale factor will be multiplied by the trigger efficiency
// given in fastsimHist.
double RazorHelper::getTriggerScaleFactor(TH2D *sfHist, TH2D *fastsimHist, float pt, float eta,
        bool isTight, bool passedTrigger, float fastsimPtCut, float ptCut) {
    double trigSF = lookupPtEtaScaleFactor( sfHist, pt, eta, 199.9, ptCut );
    if (isFastsim) {
        // note that the fastsim trigger histograms have pt on the y-axis
        trigSF *= lookupEtaPtScaleFactor( fastsimHist, pt, eta, 999.9, fastsimPtCut );
    }
    if (passedTrigger && isTight){
        return trigSF;
    }
    return 1.0;
}


double RazorHelper::getTriggerEfficiency(TH2D *effHist, float pt, float eta,
        bool isTight, bool passedTrigger, float ptCut) {
    double trigEff = lookupEtaPtScaleFactor( effHist, pt, eta, 199.9, ptCut );
    if (passedTrigger && isTight){
        return trigEff;
    }
    return 1.0;
}

double RazorHelper::getSingleMuTriggerScaleFactor(float pt, float eta, bool isTight, bool passedTrigger) {
    return getTriggerScaleFactor( muTrigSFHist, muTrigEffHist, pt, eta,
            isTight, passedTrigger, 15.01, 20.01 );
}

double RazorHelper::getSingleMuTriggerEfficiency(float pt, float eta, bool isTight, bool passedTrigger) {
    if (muTrigEffHist) {
        return getTriggerEfficiency(muTrigEffHist, pt, eta, isTight, passedTrigger, 20.01);
    }
    std::cout << "[WARNING] Single muon trigger efficiency requested, but no histogram is available.  Returning 0" << std::endl;
    return 0;
}

double RazorHelper::getSingleEleTriggerEfficiency(float pt, float eta, bool isTight, bool passedTrigger) {
    if (eleTrigEffHist) {
        return getTriggerEfficiency(eleTrigEffHist, pt, eta, isTight, passedTrigger, 25.01);
    }
    std::cout << "[WARNING] Single electron trigger efficiency requested, but no histogram is available.  Returning 0" << std::endl;
    return 0;
}

double RazorHelper::getSingleEleTriggerScaleFactor(float pt, float eta, bool isTight, bool passedTrigger) {
    return getTriggerScaleFactor( eleTrigSFHist, eleTrigEffHist, pt, eta,
            isTight, passedTrigger, 25.01, 25.01 );
}

// Update the provided trigger scale factors in-place by multiplying them by the SFs
// corresponding to the given values of pt and eta
void RazorHelper::updateTriggerScaleFactors(TH2D *sfHist, TH2D *fastsimHist,
        float pt, float eta, bool isTight, bool passedTrigger, float &sf, float &sfUp, float &sfDown,
        float fastsimPtCut, float extraSyst) {
    double trigSF = lookupPtEtaScaleFactor( sfHist, pt, eta );
    double trigSFErr = lookupPtEtaScaleFactorError( sfHist, pt, eta );
    double trigSFUp = trigSF + trigSFErr + extraSyst;
    double trigSFDown = trigSF - trigSFErr - extraSyst;
    if (passedTrigger && isTight){
        sf *= trigSF;
        sfUp *= trigSFUp/trigSF;
        sfDown *= trigSFDown/trigSF;
    }
    if (isFastsim) {
        if (passedTrigger && isTight) {
            // note that the fastsim trigger histograms have pt on the y-axis
            double trigSFFastsim = lookupEtaPtScaleFactor( fastsimHist, pt, eta, 999.9, fastsimPtCut );
            sf *= trigSFFastsim;
        }
    }
}

void RazorHelper::updateSingleMuTriggerScaleFactors(float pt, float eta, bool isTight,
        bool passedTrigger, float &sf, float &sfUp, float &sfDown) {
    updateTriggerScaleFactors( muTrigSFHist, muTrigEffHist, pt, eta, isTight,
            passedTrigger, sf, sfUp, sfDown, 15.01 );
}

void RazorHelper::updateSingleEleTriggerScaleFactors(float pt, float eta, bool isTight,
        bool passedTrigger, float &sf, float &sfUp, float &sfDown) {
    updateTriggerScaleFactors( eleTrigSFHist, eleTrigEffHist, pt, eta, isTight,
            passedTrigger, sf, sfUp, sfDown, 25.01, 0.02 );
}

double RazorHelper::getDiphotonTrigLeadingLegEff(float pt, float eta) {
  double sf = 1.0;
  if (diphotonTrigLeadingLegEffHist) sf = lookupEtaPtScaleFactor( diphotonTrigLeadingLegEffHist, pt, eta, 20.01, 99.9 );
  return sf;
}

double RazorHelper::getDiphotonTrigLeadingLegEffSF(float pt, float eta) {
  double sf = 1.0;
  if (diphotonTrigLeadingLegEffSFHist) sf = lookupPtEtaScaleFactor( diphotonTrigLeadingLegEffSFHist, pt, eta, 20.01, 99.9 );
  return sf;
}

double RazorHelper::getDiphotonTrigTrailingLegEff(float pt, float eta) {
  double sf = 1.0;
  if (diphotonTrigTrailingLegEffHist) sf = lookupEtaPtScaleFactor( diphotonTrigTrailingLegEffHist, pt, eta, 20.01, 99.9 );
  return sf;
}

double RazorHelper::getDiphotonTrigTrailingLegEffSF(float pt, float eta) {
  double sf = 1.0;
  if (diphotonTrigTrailingLegEffSFHist) sf = lookupPtEtaScaleFactor( diphotonTrigTrailingLegEffSFHist, pt, eta, 20.01, 99.9 );
  return sf;
}

// Retrieve jet energy uncertainty as a function of pt and eta
double RazorHelper::getJecUnc( float pt, float eta , int run) {

  int foundIndex = -1;
  for (uint i=0; i<JetCorrectionsIOV.size(); i++) {
    if (run >= JetCorrectionsIOV[i].first && run <= JetCorrectionsIOV[i].second) {
      foundIndex = i;
    }
  }
  if (foundIndex == -1) {
    std::cout << "Warning: run = " << run << " was not found in any valid IOV range. use default index = 0 for Jet energy corrections. \n";
    foundIndex = 0;
  }

  jecUnc[foundIndex]->setJetPt(pt);
  jecUnc[foundIndex]->setJetEta(eta);
  return jecUnc[foundIndex]->getUncertainty(true);
}

// Get one b-tag scale factor
double RazorHelper::getBTagScaleFactor(float pt, float eta, int flavor, bool isCSVM) {
    // Get efficiency and jet flavor
    double effMedium = 0;
    BTagEntry::JetFlavor jetType = BTagEntry::FLAV_B;
    if ( abs(flavor) == 5) {
        effMedium = lookupPtEtaScaleFactor( btagMediumEfficiencyHist, pt, eta );
        jetType = BTagEntry::FLAV_B;
    }
    else if ( abs(flavor) == 4) {
        effMedium = lookupPtEtaScaleFactor( btagMediumCharmEfficiencyHist, pt, eta );
        jetType = BTagEntry::FLAV_C;
    }
    else {
        effMedium = lookupPtEtaScaleFactor( btagMediumLightJetsEfficiencyHist, pt, eta );
        jetType = BTagEntry::FLAV_UDSG;
    }

    // Get scale factor
    if (pt >= 670) pt = 669; //670 is the largest pt range listed in the CSV text file
    double jet_scalefactor = -1;
    if ( abs(flavor) == 5 || abs(flavor) == 4 ) {
        jet_scalefactor = btagreader->eval(jetType, eta, pt);
    }
    else {
        jet_scalefactor = 0.907317;
    }

    //correct efficiency for Fastsim
    //Do this only for b-jets for now
    if (isFastsim && abs(flavor) == 5) {
        jet_scalefactor *= btagreaderfastsim->eval(jetType, eta, pt);
    }

    return getPassOrFailScaleFactor( effMedium, jet_scalefactor, isCSVM );
}

// Update the provided b-tag scale factors in-place by multiplying by the SFs corresponding
// to the given values of pt, eta, and flavor.
void RazorHelper::updateBTagScaleFactors(float pt, float eta, int flavor, bool isCSVM,
        float &sf, float &sfUp, float &sfDown, float &sfFastsimUp, float &sfFastsimDown,
        float &sfMistagUp, float &sfMistagDown) {
    // Get efficiency and jet flavor
    double effMedium = 0;
    BTagEntry::JetFlavor jetType = BTagEntry::FLAV_B;
    if ( abs(flavor) == 5) {
        effMedium = lookupPtEtaScaleFactor( btagMediumEfficiencyHist, pt, eta );
        jetType = BTagEntry::FLAV_B;
    }
    else if ( abs(flavor) == 4) {
        effMedium = lookupPtEtaScaleFactor( btagMediumCharmEfficiencyHist, pt, eta );
        jetType = BTagEntry::FLAV_C;
    }
    else {
        effMedium = lookupPtEtaScaleFactor( btagMediumLightJetsEfficiencyHist, pt, eta );
        jetType = BTagEntry::FLAV_UDSG;
    }

    // Get scale factor
    if (pt >= 670) pt = 669; //670 is the largest pt range listed in the CSV text file
    double jet_scalefactor = -1;
    double jet_scalefactorUp = -1;
    double jet_scalefactorDown = -1;
    if ( abs(flavor) == 5 || abs(flavor) == 4 ) {
        jet_scalefactor = btagreader->eval(jetType, eta, pt);
        jet_scalefactorUp = btagreader_up->eval(jetType, eta, pt);
        jet_scalefactorDown = btagreader_do->eval(jetType, eta, pt);
    }
    else {
        jet_scalefactor = 0.907317;
        jet_scalefactorUp = 1.257317;
        jet_scalefactorDown = 0.557317;
    }
    double jet_scalefactorFastsimUp = jet_scalefactor;
    double jet_scalefactorFastsimDown = jet_scalefactor;

    //correct efficiency for Fastsim
    //Do this only for b-jets for now
    if (isFastsim && abs(flavor) == 5) {
        double jet_scalefactorFastsim = btagreaderfastsim->eval(jetType, eta, pt);
        jet_scalefactor *= jet_scalefactorFastsim;
        jet_scalefactorUp *= jet_scalefactorFastsim;
        jet_scalefactorDown *= jet_scalefactorFastsim;
        jet_scalefactorFastsimUp *= btagreaderfastsim_up->eval(jetType, eta, pt);
        jet_scalefactorFastsimDown *= btagreaderfastsim_do->eval(jetType, eta, pt);
    }

    //apply and propagate scale factor
    double tmpSF = getPassOrFailScaleFactor( effMedium, jet_scalefactor, isCSVM );
    sf *= tmpSF;
    if (abs(flavor) == 5 || abs(flavor) == 4) {
        sfUp *= getPassOrFailScaleFactor( effMedium, jet_scalefactorUp, isCSVM ) / tmpSF;
        sfDown *= getPassOrFailScaleFactor( effMedium, jet_scalefactorDown, isCSVM ) / tmpSF;
        sfFastsimUp *= getPassOrFailScaleFactor( effMedium, jet_scalefactorFastsimUp, isCSVM ) / tmpSF;
        sfFastsimDown *= getPassOrFailScaleFactor( effMedium, jet_scalefactorFastsimDown, isCSVM ) / tmpSF;
    }
    else {
        sfMistagUp *= getPassOrFailScaleFactor( effMedium, jet_scalefactorUp, isCSVM ) / tmpSF;
        sfMistagDown *= getPassOrFailScaleFactor( effMedium, jet_scalefactorDown, isCSVM ) / tmpSF;
    }
}

// top pt reweighting from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
float RazorHelper::getTopPtWeight( float ptT, float ptTbar ) {
    // do not extrapolate correction beyond pt = 800 GeV
    if( ptT > 800 ) ptT = 800;
    if( ptTbar > 800 ) ptTbar = 800;

    // fit parameters recommended for 13 TeV data
    float a = 0.0615;
    float b = -0.0005;

    // weight from top
    float weightT = exp( a + b*ptT );
    // weight from antitop
    float weightTbar = exp( a + b*ptTbar );

    // combine into total weight
    return sqrt( weightT * weightTbar );
}

// electron scale corrections, derived privately on 2016 data
float RazorHelper::getElectronScaleCorrection( float eta ) {
    if ( fabs(eta) < 0.4 ) {
        return 1/0.993387;
    }
    else if ( fabs(eta) < 0.8 ) {
        return 1/0.993516;
    }
    else if ( fabs(eta) < 1.4442 ) {
        return 1/0.990877;
    }
    else {
        return 1/0.998137;
    }
}

float RazorHelper::getElectronResCorrection( float eta ) {
    if ( fabs(eta) < 0.4 ) {
        return 0.0137622;
    }
    else if ( fabs(eta) < 0.8 ) {
        return 0.961012;
    }
    else if ( fabs(eta) < 1.4442 ) {
        return 1.01832;
    }
    else {
        return 1.29274;
    }
}

float RazorHelper::getCorrectedElectronPt( float pt, float eta ) {
    return gRandom->Gaus( pt*getElectronScaleCorrection(eta), getElectronResCorrection(eta) );
}

float RazorHelper::getSoftDropMassCorrectionForWTag(float pt, float eta) {
    // The mass correction is explained here:
    // https://github.com/cms-jet/PuppiSoftdropMassCorr

    float genCorr = puppiSoftDropCorr_Gen->Eval(pt);
    float recoCorr = 1.0;
    if( fabs(eta) < 1.3 ) {
        recoCorr = puppiSoftDropCorr_RecoCentral->Eval(pt);
    }
    else {
        recoCorr = puppiSoftDropCorr_RecoForward->Eval(pt);
    }
    return genCorr * recoCorr;
}

bool RazorHelper::isWTaggedAK8Jet(RazorAnalyzer *ra, uint iJet, bool isData, int updown) {
    // updown: int indicating upward/downward variation to apply.
    //  if positive, will vary soft drop mass upward according to the uncertainty.
    //  if negative, will vary it downward.
    //  if zero, will use the nominal mass

    // See comments at CalcAK8JetInfo()
    float softDropMass = ra->fatJetUncorrectedSoftDropM[iJet];
    if (!isData) {
        softDropMass *= getSoftDropMassCorrectionForWTag(
                ra->fatJetCorrectedPt[iJet], ra->fatJetEta[iJet]);
        if (updown > 0) {
            softDropMass *= 1.0094;
        }
        else if (updown < 0) {
            softDropMass /= 1.0094;
        }
    }
    if (softDropMass < 65 || softDropMass > 105) return false;
    if (ra->fatJetTau2[iJet] / ra->fatJetTau1[iJet] > 0.4) return false;
    return true;
}

bool RazorHelper::isTopTaggedAK8Jet(RazorAnalyzer *ra, uint iJet) {
    // See comments at CalcAK8JetInfo()
    float softDropMass = ra->fatJetCorrectedSoftDropM[iJet];
    if (softDropMass < 105 || softDropMass > 210) return false;
    if (ra->fatJetTau3[iJet] / ra->fatJetTau2[iJet] > 0.46) return false;
    if (ra->fatJetMaxSubjetCSV[iJet] < 0.5426) return false;
    return true;
}

// Retrieve efficiency for W/top tag.  updown controls systematic variations for fastsim:
// updown > 0: vary efficiency up
// updown < 0: vary efficiency down
// updown = 0: use nominal efficiency
float RazorHelper::getTagEfficiency(TH1F *effHist, float genPt, int updown) {
    if ( genPt > effHist->GetXaxis()->GetXmax() ) {
        genPt = effHist->GetXaxis()->GetXmax() * 0.999;
    }
    float eff = effHist->GetBinContent(effHist->FindFixBin(genPt));
    float effErr = effHist->GetBinError(effHist->FindFixBin(genPt));
    if (updown > 0) eff += effErr;
    else if (updown < 0) eff -= effErr;
    return eff;
}

float RazorHelper::getWTagEfficiency(float genWPt, int updown) {
    if ( isFastsim ) {
        return getTagEfficiency(wTagEffFastsim, genWPt, updown);
    }
    return getTagEfficiency(wTagEffFullsim, genWPt, updown);
}

float RazorHelper::getWTagFastsimSF(float genWPt, int updown) {
    return getTagEfficiency(wTagEffFastsimSF, genWPt, updown);
}

float RazorHelper::getTopTagEfficiency(float genTopPt, int updown) {
    if ( isFastsim ) {
        return getTagEfficiency(topTagEffFastsim, genTopPt, updown);
    }
    return getTagEfficiency(topTagEffFullsim, genTopPt, updown);
}

float RazorHelper::getTopTagFastsimSF(float genTopPt, int updown) {
    return getTagEfficiency(topTagEffFastsimSF, genTopPt, updown);
}

RazorHelper::AK8JetInfo RazorHelper::CalcAK8JetInfo(RazorAnalyzer *ra, bool isData) {
    // For the 2016 inclusive razor analysis,
    // we use the top/W tagging guidelines listed here:
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetTopTagging
    // I'm hard-coding the cuts and scale factor/uncertainty
    // values here for now; if they change often in the future
    // it might be worth breaking them out into a separate struct
    // that is passed into this function.

    const int AK8_PT_CUT = 200;
    const float AK8_ETA_CUT = 2.4;

    const float W_TAG_SF = 1.0;
    const float W_TAG_SF_UP = 1.06;
    const float W_TAG_SF_DOWN = 0.94;

    const int TOP_TAG_PT_CUT = 400;
    const float TOP_TAG_SF = 1.05;
    const float TOP_TAG_SF_UP = 1.12;
    const float TOP_TAG_SF_DOWN = 1.01;

    AK8JetInfo jetInfo;

    // Loop over AK8 jets to count tags and update scale factors
    for (unsigned int iJet = 0; iJet < ra->nFatJets; iJet++) {
        // Baseline cuts on pt and eta
        if ( ra->fatJetCorrectedPt[iJet] < AK8_PT_CUT ) continue;
        if ( fabs(ra->fatJetEta[iJet]) > AK8_ETA_CUT ) continue;
        if ( !ra->fatJetPassIDLoose[iJet] ) continue;
        if ( ra->matchesVetoLepton(ra->fatJetEta[iJet], ra->fatJetPhi[iJet], 0.8) ) continue;

        // W tagging
        bool isWTagged = isWTaggedAK8Jet(ra, iJet, isData);
        if ( isWTagged ) {
            jetInfo.nWTags++;
        }
        // Apply scale factor if it matches a gen W
        int genWIndex = ra->getMatchingGenWIndex(
                ra->fatJetEta[iJet], ra->fatJetPhi[iJet]);
        bool matchesGenW = (genWIndex >= 0);
        // Note: only consider hadronic Ws
        if (matchesGenW && ra->isHadronicDecay(genWIndex)) {
            float genWPt = ra->gParticlePt[genWIndex];
            float eff = getWTagEfficiency(genWPt);
            jetInfo.wTagScaleFactor *= getPassOrFailScaleFactor(
                    eff, W_TAG_SF, isWTagged);
            jetInfo.wTagScaleFactor_Tau21Up *= getPassOrFailScaleFactor(
                    eff, W_TAG_SF_UP, isWTagged);
            jetInfo.wTagScaleFactor_Tau21Down *= getPassOrFailScaleFactor(
                    eff, W_TAG_SF_DOWN, isWTagged);
            if (isFastsim) {
                jetInfo.wTagScaleFactor_FastsimEffUp *= getPassOrFailScaleFactor(
                        eff, W_TAG_SF, isWTagged);
                jetInfo.wTagScaleFactor_FastsimEffDown *= getPassOrFailScaleFactor(
                        eff, W_TAG_SF, isWTagged);

                float fastsimSF = getWTagFastsimSF(genWPt);
                jetInfo.wTagScaleFactor *= getPassOrFailScaleFactor(
                        eff*fastsimSF, fastsimSF, isWTagged);

                float fastsimSFUp = getWTagFastsimSF(genWPt, 1);
                float fastsimSFDown = getWTagFastsimSF(genWPt, -1);
                jetInfo.wTagScaleFactor_FastsimEffUp *= getPassOrFailScaleFactor(
                        eff*fastsimSFUp, fastsimSFUp, isWTagged);
                jetInfo.wTagScaleFactor_FastsimEffDown *= getPassOrFailScaleFactor(
                        eff*fastsimSFDown, fastsimSFDown, isWTagged);
            }
        }
        // Up/down variations of PUPPI soft drop mass
        if ( isWTaggedAK8Jet(ra, iJet, isData, 1) ) {
            jetInfo.nWTags_SDMassUp++;
        }
        if ( isWTaggedAK8Jet(ra, iJet, isData, -1) ) {
            jetInfo.nWTags_SDMassDown++;
        }

        // Top tagging
        if ( ra->fatJetCorrectedPt[iJet] < TOP_TAG_PT_CUT ) continue;
        bool isTopTagged = isTopTaggedAK8Jet(ra, iJet);
        if ( isTopTagged ) {
            jetInfo.nTopTags++;
        }
        // Apply scale factor if it matches a gen top
        int genTopIndex = ra->getMatchingGenTopIndex(
                ra->fatJetEta[iJet], ra->fatJetPhi[iJet]);
        bool matchesGenTop = (genTopIndex >= 0);
        // Note: consider tops regardless of decay mode
        if (matchesGenTop) {
            float genTopPt = ra->gParticlePt[genTopIndex];
            float eff = getTopTagEfficiency(genTopPt);
            jetInfo.topTagScaleFactor *= getPassOrFailScaleFactor(
                    eff, TOP_TAG_SF, isTopTagged);
            jetInfo.topTagScaleFactor_Tau32Up *= getPassOrFailScaleFactor(
                    eff, TOP_TAG_SF_UP, isTopTagged);
            jetInfo.topTagScaleFactor_Tau32Down *= getPassOrFailScaleFactor(
                    eff, TOP_TAG_SF_DOWN, isTopTagged);
            if (isFastsim) {
                jetInfo.topTagScaleFactor_FastsimEffUp *= getPassOrFailScaleFactor(
                        eff, TOP_TAG_SF, isTopTagged);
                jetInfo.topTagScaleFactor_FastsimEffDown *= getPassOrFailScaleFactor(
                        eff, TOP_TAG_SF, isTopTagged);

                float fastsimSF = getTopTagFastsimSF(genTopPt);
                jetInfo.topTagScaleFactor *= getPassOrFailScaleFactor(
                        eff*fastsimSF, fastsimSF, isTopTagged);

                float fastsimSFUp = getTopTagFastsimSF(genTopPt, 1);
                float fastsimSFDown = getTopTagFastsimSF(genTopPt, -1);
                jetInfo.topTagScaleFactor_FastsimEffUp *= getPassOrFailScaleFactor(
                        eff*fastsimSFUp, fastsimSFUp, isTopTagged);
                jetInfo.topTagScaleFactor_FastsimEffDown *= getPassOrFailScaleFactor(
                        eff*fastsimSFDown, fastsimSFDown, isTopTagged);
            }
        }
    }

    return jetInfo;
}
