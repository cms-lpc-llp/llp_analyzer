#include "LLPAnalysis/llpAnalyzer/interface/SusyLLP.h"
#include "LLPAnalysis/llpAnalyzer/interface/RazorHelper.h"
#include "LLPAnalysis/llpAnalyzer/interface/SusyLLPTree.h"
#include "LLPAnalysis/llpAnalyzer/interface/JetCorrectorParameters.h"
#include "LLPAnalysis/llpAnalyzer/interface/JetCorrectionUncertainty.h"
#include "LLPAnalysis/llpAnalyzer/interface/BTagCalibrationStandalone.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"

#define xclean 1
#define presel 0
#define presel_jet 1

#define _debug 0
#define _debug_pdf 0
#define _debug_dnn 0
#define _debug_pf 0
#define _debug_lab 0
#define _debug_npu 0
#define _debug_sync 0
#define _debug_met 0
#define _debug_jet 0
#define _debug_trk 0
#define _debug_lep 0
#define _debug_match 0
#define _debug_avgH 0
#define _debug_trg 0
#define _debug_pre 0
#define _debug_csc 0

#define N_MAX_LLP_DAUGHTERS 4
#define N_MAX_LLP_GRAND_DAUGHTERS 4
#define N_MAX_LLPS 2
#define N_MAX_LEPTONS 20
#define N_MAX_JETS 20
#define NTriggersMAX 1201 //Number of trigger in the .dat file
using namespace std;

struct leptons
{
	TLorentzVector lepton;
	int pdgId;
	// float dZ;
	// bool passLooseId;
	// bool passMediumId;
	// bool passTightId;
	// bool passId;
};

struct taus
{
	TLorentzVector tau;
};

struct photons
{
	TLorentzVector photon;
};

struct fatjets
{
	TLorentzVector fatjet;
	float CorrectedPt;

	//ecal rechits
	int fatjetNRecHitsEcal;
	float fatjetEnergyRecHitsEcal;
	float fatjetTimeRecHitsEcal;

	//hcal hbhe rechits
	int fatjetNRecHitsHcal;
	float fatjetEnergyRecHitsHcal;
	float fatjetTimeRecHitsHcal;

};


struct jets
{
	TLorentzVector jet;
	//pat::Jet jet;
	float time;
	int jetNeutralHadronMultiplicity;
	int jetChargedHadronMultiplicity;
	int jetMuonMultiplicity;
	int jetElectronMultiplicity;
	int jetPhotonMultiplicity;
	float jetNeutralHadronEnergyFraction;
	float jetChargedHadronEnergyFraction;
	float jetMuonEnergyFraction;
	float jetElectronEnergyFraction;
	float jetPhotonEnergyFraction;
	float jetCSV;
	int ecalNRechits;
	float ecalRechitE;

	float jetAlphaMax;
	float jetBetaMax;
	float jetGammaMax;
	float jetGammaMax_Hadronic;
	float jetGammaMax_EM;
	float jetGammaMax_ET;
	float jetPtAllTracks;
	float jetPtAllPVTracks;
	float jetMinDeltaRAllTracks;
	float jetMinDeltaRPVTracks;

	float jetChargedEMEnergyFraction;
	float jetNeutralEMEnergyFraction;

	int jetChargedMultiplicity;
	//int jetNHits;
	//int jetNPixelHits;
	float jetNPixelHitsMedian;
	float jetNHitsMedian;

	//int   jetNSV;
	//int   jetNSVCand;
	int   jetNVertexTracks;
	int   jetNSelectedTracks;
	float jetDRSVJet;
	//float jetFlightDist2D;
	//float jetFlightDist2DError;
	//float jetFlightDist3D;
	//float jetFlightDist3DError;
	//float jetSV_x;
	//float jetSV_y;
	//float jetSV_z;
	//int   jetSVNTracks;
	float jetSVMass;

	//ecal rechits
	int jetNRecHitsEcal;
	float jetEnergyRecHitsEcal;
	float jetTimeRecHitsEcal;
	float jetTimeRmsEcal;
	float jetTimeRmsEcal_fix;

	//hcal hbhe rechits
	int jetNRecHitsHcal;
	float jetEnergyRecHitsHcal;
	float jetTimeRecHitsHcal;

	double jetsig1EB;
	double jetsig2EB;
	double jetptDEB;
	double jetsig1PF;
	double jetsig2PF;
	double jetptDPF;

	//float dnn_score_v1;
	//float dnn_score;
	float dnn_score_v3;
	float dnn_score_v3_miniAOD;

	//muon sys var
	float jetCscEF;
	float jetDtEF;

	//max trk pt
	int maxpvid;
	int maxsvid;
	float maxpvpt;
	float maxsvpt;
};

//pt comparison
//not used so far
struct greater_than_pt
{
	inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
};

//lepton highest pt comparator
struct largest_pt_lep
{
	inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt_lep;

//jet highest pt comparator
struct largest_pt_jet
{
	inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;

//fatjet highest pt comparator
struct largest_pt_fatjet
{
	inline bool operator() (const fatjets& p1, const fatjets& p2){return p1.fatjet.Pt() > p2.fatjet.Pt();}
} my_largest_pt_fatjet;

//Analyze
void SusyLLP::Analyze(bool isData, int options, string outputfilename, string analysisTag, string process)
{
	//initialization: create one TTree for each analysis box
	cout << "Initializing..." << endl;
	cout << "IsData = " << isData << "\n";
	cout << "options = " << options << "\n";

	//---------------------------
	//-----------NN Setup----------
	//---------------------------
	///	
	//	//-----------v1-----------
	//	std::string basePathV1 = std::string(std::getenv("CMSSW_BASE")) + "/src/LLPAnalysis/llpAnalyzer/nn_inference";
	//	
	//	std::string graphPathV1 = basePathV1 + "/graph_NoHCAL_NoSi.pb";
	//	std::string inputTensorNameV1 = "input_input";
	//	std::string outputTensorNameV1 = "FCN/output/Softmax";//"FCN/dense_4/Softmax";//or Softmax?
	//	
	//	// threading setup
	//	// to enable tensorflow-native multi-threading, change to "tensorflow" and increase nThreads
	//	std::string threadPoolV1 = "no_threads";
	//	int nThreadsV1 = 1;
	//	
	//	std::vector<std::string> inputFeaturesV1 = { "Jet_0_nTrackConstituents","Jet_0_nSelectedTracks", "Jet_0_timeRecHitsEB", "Jet_0_energyRecHitsEB", "Jet_0_nRecHitsEB", "Jet_0_cHadEFrac", "Jet_0_nHadEFrac", "Jet_0_eleEFrac", "Jet_0_photonEFrac", "Jet_0_ptAllTracks", "Jet_0_ptAllPVTracks", "Jet_0_alphaMax", "Jet_0_betaMax", "Jet_0_gammaMax", "Jet_0_gammaMaxEM", "Jet_0_gammaMaxHadronic", "Jet_0_gammaMaxET","Jet_0_minDeltaRAllTracks","Jet_0_minDeltaRPVTracks",};
	//
	//	int nInputsV1 = inputFeaturesV1.size();
	//	std::vector<float> inputValuesV1(nInputsV1);
	//	
	//	// setup TensorFlow objects
	//	tensorflow::setLogging();
	//	tensorflow::GraphDef* graphDefV1 = tensorflow::loadGraphDef(graphPathV1);
	//	tensorflow::Session* sessionV1 = tensorflow::createSession(graphDefV1, nThreadsV1);
	//		
	//	// register an input tensor (1 x nInputs) that is filled during the event loop
	//	tensorflow::Tensor inputTensorV1(tensorflow::DT_FLOAT, {1, nInputsV1});
	//	
	//
	//	//-----------v2-----------
	//	std::string basePath = std::string(std::getenv("CMSSW_BASE")) + "/src/LLPAnalysis/llpAnalyzer/nn_inference/tagger_AK4_v2";
	//	
	//	std::string graphPath = basePath + "/graph.pb";
	//	std::string inputTensorName = "input_input";
	//	std::string outputTensorName = "FCN/output/Softmax";//"FCN/dense_4/Softmax";//or Softmax?
	//	
	//	// threading setup
	//	// to enable tensorflow-native multi-threading, change to "tensorflow" and increase nThreads
	//	std::string threadPool = "no_threads";
	//	int nThreads = 1;
	//	
	//	std::vector<std::string> inputFeatures = { "Jet_nTrackConstituents", "Jet_nSelectedTracks", "Jet_timeRecHitsEB", "Jet_eFracRecHitsEB", "Jet_nRecHitsEB", "Jet_sig1EB", "Jet_sig2EB", "Jet_ptDEB", "Jet_sig1PF", "Jet_sig2PF", "Jet_ptDPF", "Jet_cHadEFrac", "Jet_nHadEFrac", "Jet_eleEFrac", "Jet_photonEFrac", "Jet_ptAllTracks", "Jet_ptAllPVTracks", "Jet_alphaMax", "Jet_betaMax", "Jet_gammaMax", "Jet_gammaMaxEM", "Jet_gammaMaxHadronic", "Jet_gammaMaxET", "Jet_minDeltaRAllTracks", "Jet_minDeltaRPVTracks",};
	//
	//	int nInputs = inputFeatures.size();
	//	std::vector<float> inputValues(nInputs);
	//	
	//	// setup TensorFlow objects
	//	tensorflow::setLogging();
	//	tensorflow::GraphDef* graphDef = tensorflow::loadGraphDef(graphPath);
	//	tensorflow::Session* session = tensorflow::createSession(graphDef, nThreads);
	//		
	//	// register an input tensor (1 x nInputs) that is filled during the event loop
	//	tensorflow::Tensor inputTensor(tensorflow::DT_FLOAT, {1, nInputs});
	//	
	//-----------v3-----------
	std::string basePathV3 = std::string(std::getenv("CMSSW_BASE")) + "/src/LLPAnalysis/llpAnalyzer/nn_inference/tagger_AK4_v3";

	std::string graphPathV3 = basePathV3 + "/graph.pb";
	std::string inputTensorNameV3 = "input_input";
	std::string outputTensorNameV3 = "FCN/output/Softmax";//"FCN/dense_4/Softmax";//or Softmax?

	// threading setup
	// to enable tensorflow-native multi-threading, change to "tensorflow" and increase nThreads
	std::string threadPoolV3 = "no_threads";
	int nThreadsV3 = 1;

	//std::vector<std::string> inputFeaturesV3 = { "Jet_nTrackConstituents", "Jet_nSelectedTracks", "Jet_timeRecHitsEB", "Jet_eFracRecHitsEB", "Jet_nRecHitsEB", "Jet_sig1EB", "Jet_sig2EB", "Jet_ptDEB", "Jet_cHadEFrac", "Jet_nHadEFrac", "Jet_eleEFrac", "Jet_photonEFrac", "Jet_ptAllTracks", "Jet_ptAllPVTracks", "Jet_alphaMax", "Jet_betaMax", "Jet_gammaMax", "Jet_gammaMaxEM", "Jet_gammaMaxHadronic", "Jet_gammaMaxET", "Jet_minDeltaRAllTracks", "Jet_minDeltaRPVTracks",};
	std::vector<std::string> inputFeaturesV3 = { "Jet_nTrackConstituents", "Jet_nSelectedTracks", "Jet_timeRecHitsEB", "Jet_eFracRecHitsEB", "Jet_nRecHitsEB", "Jet_sig1EB", "Jet_sig2EB", "Jet_ptDEB", "Jet_cHadEFrac", "Jet_nHadEFrac", "Jet_eleEFrac", "Jet_photonEFrac", "Jet_ptAllTracks", "Jet_ptAllPVTracks", "Jet_alphaMax", "Jet_betaMax", "Jet_gammaMax", "Jet_gammaMaxEM", "Jet_gammaMaxHadronic", "Jet_gammaMaxET", "Jet_minDeltaRAllTracks", "Jet_minDeltaRPVTracks",};

	int nInputsV3 = inputFeaturesV3.size();
	std::vector<float> inputValuesV3(nInputsV3);

	// setup TensorFlow objects
	tensorflow::setLogging();
	tensorflow::GraphDef* graphDefV3 = tensorflow::loadGraphDef(graphPathV3);
	tensorflow::Session* sessionV3 = tensorflow::createSession(graphDefV3, nThreadsV3);

	// register an input tensor (1 x nInputs) that is filled during the event loop
	tensorflow::Tensor inputTensorV3(tensorflow::DT_FLOAT, {1, nInputsV3});


	//-----------v3 miniAOD-----------
	std::string basePathV3miniAOD = std::string(std::getenv("CMSSW_BASE")) + "/src/LLPAnalysis/llpAnalyzer/nn_inference/tagger_AK4_miniAOD_v3";

	std::string graphPathV3miniAOD = basePathV3miniAOD + "/graph.pb";
	std::string inputTensorNameV3miniAOD = "input_input";
	std::string outputTensorNameV3miniAOD = "FCN/output/Softmax";//"FCN/dense_4/Softmax";//or Softmax?

	// threading setup
	// to enable tensorflow-native multi-threading, change to "tensorflow" and increase nThreads
	std::string threadPoolV3miniAOD = "no_threads";
	int nThreadsV3miniAOD = 1;

	//'Jet_nTrackConstituents', 'Jet_nSelectedTracks', 'Jet_timeRecHitsEB', 'Jet_eFracRecHitsEB', 'Jet_nRecHitsEB', 'Jet_sig1EB', 'Jet_sig2EB', 'Jet_ptDEB', 'Jet_cHadEFrac', 'Jet_nHadEFrac', 'Jet_eleEFrac', 'Jet_photonEFrac'
	std::vector<std::string> inputFeaturesV3miniAOD = { "Jet_nTrackConstituents", "Jet_nSelectedTracks", "Jet_timeRecHitsEB", "Jet_eFracRecHitsEB", "Jet_nRecHitsEB", "Jet_sig1EB", "Jet_sig2EB", "Jet_ptDEB", "Jet_cHadEFrac", "Jet_nHadEFrac", "Jet_eleEFrac", "Jet_photonEFrac",};

	int nInputsV3miniAOD = inputFeaturesV3miniAOD.size();
	std::vector<float> inputValuesV3miniAOD(nInputsV3miniAOD);

	// setup TensorFlow objects
	tensorflow::setLogging();
	tensorflow::GraphDef* graphDefV3miniAOD = tensorflow::loadGraphDef(graphPathV3miniAOD);
	tensorflow::Session* sessionV3miniAOD = tensorflow::createSession(graphDefV3miniAOD, nThreadsV3miniAOD);

	// register an input tensor (1 x nInputs) that is filled during the event loop
	tensorflow::Tensor inputTensorV3miniAOD(tensorflow::DT_FLOAT, {1, nInputsV3miniAOD});

	//---------------------------
	//-----------option----------
	//---------------------------
	int option;
	std::string label;
	bool signalScan;
	bool signalHZScan;
	bool signalHDecayScan;

	//HUNDRED'S DIGIT
	//option of run condor or locally
	if (options < 200){
		option = 1; // used when running condor
	}
	else{
		option = 0;// used when running locally
		//options need to be larger than 200, if do local test run
		cout << "option = 0, running locally, load aux locally \n";
	}

	//TEN'S DIGIT
	// label of signal / bkg
	if ((options/10)%10 == 1){
		label = "wH";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 2){
		//label = "zH";
		label = "MR_JetHT";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 3){
		//label = "bkg_wH";
		label = "MR_SingleMuon";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 4){
		//label = "bkg_zH";
		label = "MR_SingleElectron";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 5){
		label = "HH";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 6){
		label = "bkg_HH";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 7){
		label = "MR_EMU";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 8){
		label = "MR_PHO";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 9){
		label = "MR_ZLL";
		cout << "process label: " << label << "\n";
	}
	else{
		cout << "What process it is? Label not defined. \n";
	}

	//UNIT'S DIGIT
	// signalScan option
	if(options%10==3){
		signalScan = false;
		signalHZScan = false;
		signalHDecayScan = true;
	}
	else if(options%10==2){
		signalScan = false;
		signalHZScan = true;
		signalHDecayScan = false;
	}
	else if(options%10==1){
		signalScan = true;
		signalHZScan = false;
		signalHDecayScan = false;
	}
	else{
		signalScan = false;
		signalHZScan = false;
		signalHDecayScan = false;
	}

	// DATA or MC
	if( isData )
	{
		std::cout << "[INFO]: running on data with label: " << label << " and option: " << option << " and signalScan is " << signalScan << " and signalHZScan is " << signalHZScan << " and signalHDecayScan is " << signalHDecayScan<< std::endl;
	}
	else
	{
		std::cout << "[INFO]: running on MC with label: " << label << " and option: " << option << " and signalScan is " << signalScan << " and signalHZScan is " << signalHZScan << " and signalHDecayScan is " << signalHDecayScan<< std::endl;
	}

	const float ELE_MASS = 0.000511;
	const float MU_MASS  = 0.105658;
	const float TAU_MASS  = 1.77686;
	const float Z_MASS   = 91.2;
	//const float H_MASS   = 125.0;

	string dataset = "94X";
	//Analysis Tag
	//reference in RazorHelper
	if (analysisTag == ""){
		analysisTag = "Razor2016_80X";
	}

	int wzId;
	int NTrigger;//Number of trigger in trigger paths
	//int elePt_cut = 0;
	//int muonPt_cut = 0;
	//uint nLepton_cut = 0;



	if (label == "HH" || label == "bkg_HH" ){
		NTrigger = 1;
		//muonPt_cut = 15;
		//elePt_cut = 15;
		//nLepton_cut = 0;
	}
	if (label == "MR_EMU" ){
		NTrigger = 10;
	}
	if (label == "MR_PHO" ){
		NTrigger = 16;
	}
	if (label == "MR_SingleMuon" ){
		NTrigger = 8;
	}
	if (label == "MR_SingleElectron" ){
		NTrigger = 2;
	}
	if (label == "MR_ZLL" ){
		NTrigger = 10;
	}
	if (label == "MR_JetHT" ){
		NTrigger = 11;
	}
	if (label == "zH" || label == "bkg_zH" ){
		NTrigger = 4;
		//muonPt_cut = 15;
		//elePt_cut = 15;
		//nLepton_cut = 2;
	}
	//else{}
	if (label == "wH" || label == "bkg_wH" ){
		NTrigger = 2;
		//muonPt_cut = 27;
		//elePt_cut = 32;
		//nLepton_cut = 1;
	}

	int trigger_paths[NTrigger];
	if (label == "wH" || label == "bkg_wH"){
		wzId = 24;
		trigger_paths[0] = 87;
		trigger_paths[1] = 135;
		// trigger_paths[2] = 310;
	}
	else if (label == "zH" || label == "bkg_zH"){
		wzId = 23;
		trigger_paths[0] = 177;
		trigger_paths[1] = 362;
		// trigger_paths[2] = 310;
		trigger_paths[2] = 87;
		trigger_paths[3] = 135;
	}
	else if (label == "HH" || label == "bkg_HH"){
		wzId = 25;
		trigger_paths[0] = 310;
	}
	else if (label == "MR_EMU"){
		wzId = 25;
		trigger_paths[0] = 372; //HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL
		trigger_paths[1] = 373; //HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ
		trigger_paths[2] = 374; //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
		trigger_paths[3] = 375; //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
		trigger_paths[4] = 376; //HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL
		trigger_paths[5] = 377; //HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL
		trigger_paths[6] = 378; //HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL
		trigger_paths[7] = 379; //HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL
		trigger_paths[8] = 615; //HLT_Mu27_Ele37_CaloIdL_MW
		trigger_paths[9] = 616; //HLT_Mu37_Ele27_CaloIdL_MW
	}
	else if (label == "MR_PHO"){
		wzId = 25;
		trigger_paths[15] = 402;  //HLT_Photon22
		trigger_paths[14] = 403;  //HLT_Photon30
		trigger_paths[13] = 774;  //HLT_Photon33
		trigger_paths[12] = 404;  //HLT_Photon36
		trigger_paths[11] = 405;  //HLT_Photon50
		trigger_paths[10] = 406;  //HLT_Photon75
		trigger_paths[9] = 407;  //HLT_Photon90
		trigger_paths[8] = 408;  //HLT_Photon120
		trigger_paths[7] = 536;  //HLT_Photon125
		trigger_paths[6] =  58;  //HLT_Photon150
		trigger_paths[5] = 775;  //HLT_Photon200
		trigger_paths[4] = 409;  //HLT_Photon175
		trigger_paths[3] = 335;  //HLT_Photon250_NoHE
		trigger_paths[2] = 336;  //HLT_Photon300_NoHE
		trigger_paths[1] = 583;  //HLT_Photon500
		trigger_paths[0] = 584;  //HLT_Photon600
	}
	else if (label == "MR_SingleMuon"){
		wzId = 25;
		trigger_paths[7] = 115; //HLT_IsoMu17_eta2p1
		trigger_paths[6] = 120; //HLT_IsoMu18
		trigger_paths[5] = 129; //HLT_IsoMu20 
		trigger_paths[4] = 133; //HLT_IsoMu22
		trigger_paths[3] = 135; //HLT_IsoMu24
		trigger_paths[2] = 136; //HLT_IsoMu27
		trigger_paths[1] = 644; //HLT_IsoMu24_eta2p1
		trigger_paths[0] = 645; //HLT_IsoMu30
	}
	else if (label == "MR_SingleElectron"){
		wzId = 25;
		trigger_paths[1] = 87; //HLT_Ele32_WPTight_Gsf
		trigger_paths[0] = 88; //HLT_Ele32_eta2p1_WPLoose_Gsf
	}
	else if (label == "MR_ZLL"){
		wzId = 25;
		trigger_paths[9] = 87; //HLT_Ele32_WPTight_Gsf
		trigger_paths[8] = 88; //HLT_Ele32_eta2p1_WPLoose_Gsf
		trigger_paths[7] = 115; //HLT_IsoMu17_eta2p1
		trigger_paths[6] = 120; //HLT_IsoMu18
		trigger_paths[5] = 129; //HLT_IsoMu20 
		trigger_paths[4] = 133; //HLT_IsoMu22
		trigger_paths[3] = 135; //HLT_IsoMu24
		trigger_paths[2] = 136; //HLT_IsoMu27
		trigger_paths[1] = 644; //HLT_IsoMu24_eta2p1
		trigger_paths[0] = 645; //HLT_IsoMu30
	}
	else if (label == "MR_JetHT"){
		wzId = 25;
		trigger_paths[10] = 239; //HLT_PFJet40
		trigger_paths[9] = 240; //HLT_PFJet60
		trigger_paths[8] = 241; //HLT_PFJet80
		trigger_paths[7] = 242; //HLT_PFJet140
		trigger_paths[6] = 243; //HLT_PFJet200
		trigger_paths[5] = 244; //HLT_PFJet260
		trigger_paths[4] = 245; //HLT_PFJet320
		trigger_paths[3] = 246; //HLT_PFJet400
		trigger_paths[2] = 247; //HLT_PFJet450
		trigger_paths[1] = 248; //HLT_PFJet500
		trigger_paths[0] = 666; //HLT_PFJet550
	}

	//--------------------------------
	//Initialize helper
	//--------------------------------
	RazorHelper *helper = 0;
	if (process == ""){
		process = "ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Summer16_2016";
	}

	helper = new RazorHelper(analysisTag, isData, false, process.c_str());

	//-----------------------------------------------
	//Set up Output File
	//-----------------------------------------------
	string outfilename = outputfilename;
	if (outfilename == "") outfilename = "SusyLLPTree.root";
	TFile *outFile;
	if(!signalScan && !signalHZScan && !signalHDecayScan) outFile = new TFile(outfilename.c_str(), "RECREATE");

	SusyLLPTree *llp_tree = new SusyLLPTree;
	llp_tree->CreateTree();
	llp_tree->tree_->SetAutoFlush(0);
	llp_tree->InitTree();

	//histogram containing total number of processed events (for normalization)
	TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
	//cut flow 
	TH1F *NEventsMet200 = new TH1F("NEventsMet200", "NEventsMet200", 1, 1, 2);
	//TH1F *NEventsTrg = new TH1F("NEventsTrg", "NEventsTrg", 1, 1, 2);
	//TH1F *NEventsFlag = new TH1F("NEventsFlag", "NEventsFlag", 1, 1, 2);
	//TH1F *NEventsLepVeto = new TH1F("NEventsLepVeto", "NEventsLepVeto", 1, 1, 2);
	//TH1F *NEventsPhoVeto = new TH1F("NEventsPhoVeto", "NEventsPhoVeto", 1, 1, 2);
	//TH1F *NEventsTauVeto = new TH1F("NEventsTauVeto", "NEventsTauVeto", 1, 1, 2);
	//TH1F *NEventsMDPhi = new TH1F("NEventsMDPhi", "NEventsMDPhi", 1, 1, 2);

	//for signals, need one output file for each signal point
	map<pair<int,int>, TFile*> Files2D;
	map<pair<int,int>, TTree*> Trees2D;
	map<pair<int,int>, TH1F*> NEvents2D;
	map<pair<int,int>, TH1F*> NEventsMet2002D;
	//H Decay NEvents
	map<pair<int,int>, TH1F*> NEventsH2D;
	map<pair<int,int>, TH1F*> NEventsHbb2D;
	map<pair<int,int>, TH1F*> NEventsHgg2D;
	map<pair<int,int>, TH1F*> NEventsHcc2D;
	map<pair<int,int>, TH1F*> NEventsHZZ2D;
	map<pair<int,int>, TH1F*> NEventsHWW2D;
	map<pair<int,int>, TH1F*> NEventsHmm2D;
	map<pair<int,int>, TH1F*> NEventsHtt2D;
	
	//HH decay NEvents
	map<pair<int,int>, TH1F*> NEventsHH2D;
	map<pair<int,int>, TH1F*> NEventsHHbb2D;
	map<pair<int,int>, TH1F*> NEventsHHgg2D;
	map<pair<int,int>, TH1F*> NEventsHHcc2D;
	map<pair<int,int>, TH1F*> NEventsHHZZ2D;
	map<pair<int,int>, TH1F*> NEventsHHWW2D;
	map<pair<int,int>, TH1F*> NEventsHHmm2D;
	map<pair<int,int>, TH1F*> NEventsHHtt2D;
	
	//HZ decay Events
	map<pair<int,int>, TH1F*> NEventsZH2D;
	map<pair<int,int>, TH1F*> NEventsZHbb2D;
	map<pair<int,int>, TH1F*> NEventsZHgg2D;
	map<pair<int,int>, TH1F*> NEventsZHcc2D;
	map<pair<int,int>, TH1F*> NEventsZHZZ2D;
	map<pair<int,int>, TH1F*> NEventsZHWW2D;
	map<pair<int,int>, TH1F*> NEventsZHmm2D;
	map<pair<int,int>, TH1F*> NEventsZHtt2D;
	

	//Scale and PDF variations
	map<int, TH1D*> smsSumWeights;
	map<int, TH1D*> smsSumScaleWeights;
	map<int, TH1D*> smsSumPdfWeights;
	map<int, TH1D*> smsSumAlphasWeights;
	map<pair<int,int>, TH1D*> smsSumWeights2D;
	map<pair<int,int>, TH1D*> smsSumScaleWeights2D;
	map<pair<int,int>, TH1D*> smsSumPdfWeights2D;
	map<pair<int,int>, TH1D*> smsSumAlphasWeights2D;

	TH1D *SumWeights = new TH1D("SumWeights", "SumWeights", 1, -0.5, 0.5);
	TH1D *SumScaleWeights = new TH1D("SumScaleWeights", "SumScaleWeights", 9, -0.5, 8.5);
	TH1D *SumPdfWeights = new TH1D("SumPdfWeights", "SumPdfWeights", 100, -0.5, 99.5);
	TH1D *SumAlphasWeights = new TH1D("SumAlphasWeights", "SumAlphasWeights", 2, -0.5, 1.5);


	// ************************************************************************
	//Look over Input File Events
	// ************************************************************************
	if (fChain == 0) return;
	cout << "Total Events: " << fChain->GetEntries() << "\n";
	Long64_t nbytes = 0, nb = 0;

	for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {

		//begin event
		if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		if(_debug) std::cout << "passed break " << jentry << std::endl;
		//GetEntry(ientry);
		nb = fChain->GetEntry(jentry); 
		if(_debug) std::cout << "passed GetEntry " << jentry << std::endl;
		nbytes += nb;
		if(_debug) std::cout << "passed nbytes " << jentry << std::endl;

		//fill normalization histogram
		if(_debug) std::cout << "deb0 " << jentry << std::endl;
		llp_tree->InitVariables();
		if(_debug) std::cout << "deb1 " << jentry << std::endl;

		//std::cout << *lheComments<<endl;
		//if(getline(parser, item, '_')) std::cout << item.c_str()<<endl;
		float weight;
		//float genWeight=1;
		////PDF SF
		//std::vector<float> sf_pdf;
		//scaleWeights = new std::vector<float>; scaleWeights->clear();
  		//pdfWeights = new std::vector<float>; pdfWeights->clear();
  		//alphasWeights = new std::vector<float>; alphasWeights->clear();
  		////For scale variation uncertainties
  		//float sf_facScaleUp, sf_facScaleDown;
  		//float sf_renScaleUp, sf_renScaleDown;	
  		//float sf_facRenScaleUp, sf_facRenScaleDown;

		if (!isData && signalScan)
		{

			string mh_substring = lheComments->substr(lheComments->find("MH-")+3);
			int mh = stoi(mh_substring.substr(0,mh_substring.find('_')));
			string mx_substring = lheComments->substr(lheComments->find("MS-")+3);
			int mx = stoi(mx_substring.substr(0,mx_substring.find('_')));
			string ctau_substring = lheComments->substr(lheComments->find("ctauS-")+6);
			int ctau = stoi(ctau_substring.substr(0,ctau_substring.find('_')));
			llp_tree->mH = mh;
			llp_tree->mX = mx;
			llp_tree->ctau = ctau;


			pair<int,int> signalPair = make_pair(mx, ctau);

			if (Files2D.count(signalPair) == 0){ //create file and tree
				//format file name
				string thisFileName = outfilename;
				thisFileName.erase(thisFileName.end()-5, thisFileName.end());
				thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";

				Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
				Trees2D[signalPair] =  llp_tree->tree_->CloneTree(0);
				NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1,0.5,1.5);
				NEventsMet2002D[signalPair] = new TH1F(Form("NEventsMet200%d%d", mx, ctau), "NEventsMet200", 1,0.5,1.5);

				smsSumWeights2D[signalPair] = new TH1D(Form("SumWeights%d%d", mh, mx), "SumWeights", 1 ,-0.5,0.5);
				smsSumScaleWeights2D[signalPair] = new TH1D(Form("SumScaleWeights%d%d", mh, mx), "SumScaleWeights", 9 ,-0.5,8.5);
				smsSumPdfWeights2D[signalPair] = new TH1D(Form("SumPdfWeights%d%d", mh, mx), "SumPdfWeights", 100 ,-0.5,99.5);
				smsSumAlphasWeights2D[signalPair] = new TH1D(Form("SumAlphasWeights%d%d", mh, mx), "SumAlphasWeights", 2 ,-0.5,1.5);


				cout << "Created new output file " << thisFileName << endl;
			}
			//Fill NEvents hist
			NEvents2D[signalPair]->Fill(1.0, genWeight);
			if(metType1Pt>200.) NEventsMet2002D[signalPair]->Fill(1.0, genWeight);
			smsSumWeights2D[signalPair]->Fill(1.0, weight);

			smsSumScaleWeights2D[signalPair]->Fill(0.0, llp_tree->sf_facScaleUp);
			smsSumScaleWeights2D[signalPair]->Fill(1.0, llp_tree->sf_facScaleDown);
			smsSumScaleWeights2D[signalPair]->Fill(2.0, llp_tree->sf_renScaleUp);
			smsSumScaleWeights2D[signalPair]->Fill(3.0, llp_tree->sf_renScaleDown);
			smsSumScaleWeights2D[signalPair]->Fill(4.0, llp_tree->sf_facRenScaleUp);
			smsSumScaleWeights2D[signalPair]->Fill(5.0, llp_tree->sf_facRenScaleDown);
			//it should goes to 8? or we only care about first 6

			for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) {
				smsSumPdfWeights2D[signalPair]->Fill(double(iwgt),(*pdfWeights)[iwgt]);
				if(_debug_pdf) cout << "pdf weight component: " << (*pdfWeights)[iwgt] << endl;
			}

			if(_debug_pdf) cout <<"alphas size: "<< alphasWeights->size() <<endl;
			for (unsigned int iwgt=0; iwgt<alphasWeights->size(); ++iwgt) {
				smsSumAlphasWeights2D[signalPair]->Fill(double(iwgt),(*alphasWeights)[iwgt]);
				if(_debug_pdf) cout<<"alphas: "<<(*alphasWeights)[iwgt] <<endl;
			}


		}
		else if (!isData && (signalHZScan || signalHDecayScan) )
		{

			//TChiHZ_HToBB_LLN2N3_150_3000
			//std::cout << *lheComments<<endl;
			int mchi = 0;
			int ctau = 0;	
			stringstream parser(*lheComments);
			string item;
			getline(parser, item, '_'); //prefix
			if(getline(parser, item, '_')) 
			{
				//std::cout << item.c_str()<<endl; //HToBB
				if(getline(parser, item, '_')) 
				{
					//std::cout << item.c_str()<<endl; //LLN2N3
					if(getline(parser, item, '_')) 
					{
						//std::cout << item.c_str()<<endl; //150
						mchi = atoi(item.c_str());
						if(getline(parser, item, '_')) 
						{
							//std::cout << item.c_str()<<endl; //3000
							ctau = atoi(item.c_str());
						}
					}
				}
			}
			//std::cout << mchi<<endl;
			//std::cout << ctau<<endl;
			llp_tree->mH = mchi;
			llp_tree->ctau = ctau;


			pair<int,int> signalPair = make_pair(mchi, ctau);

			if (Files2D.count(signalPair) == 0){ //create file and tree
				//format file name
				string thisFileName = outfilename;
				thisFileName.erase(thisFileName.end()-5, thisFileName.end());
				thisFileName += "_" + to_string(mchi) + "_" + to_string(ctau) + ".root";

				Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
				Trees2D[signalPair] =  llp_tree->tree_->CloneTree(0);
				NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mchi, ctau), "NEvents", 1,0.5,1.5);
				NEventsMet2002D[signalPair] = new TH1F(Form("NEventsMet200%d%d", mchi, ctau), "NEventsMet200", 1,0.5,1.5);

				// H decay
				NEventsH2D[signalPair] = new TH1F(Form("NEventsH%d%d", mchi, ctau), "NEventsH", 1,0.5,1.5);
				NEventsHbb2D[signalPair] = new TH1F(Form("NEventsHbb%d%d", mchi, ctau), "NEventsHbb", 1,0.5,1.5);
				NEventsHcc2D[signalPair] = new TH1F(Form("NEventsHcc%d%d", mchi, ctau), "NEventsHcc", 1,0.5,1.5);
				NEventsHgg2D[signalPair] = new TH1F(Form("NEventsHgg%d%d", mchi, ctau), "NEventsHgg", 1,0.5,1.5);
				NEventsHWW2D[signalPair] = new TH1F(Form("NEventsHWW%d%d", mchi, ctau), "NEventsHWW", 1,0.5,1.5);
				NEventsHZZ2D[signalPair] = new TH1F(Form("NEventsHZZ%d%d", mchi, ctau), "NEventsHZZ", 1,0.5,1.5);
				NEventsHmm2D[signalPair] = new TH1F(Form("NEventsHmm%d%d", mchi, ctau), "NEventsHmm", 1,0.5,1.5);
				NEventsHtt2D[signalPair] = new TH1F(Form("NEventsHtt%d%d", mchi, ctau), "NEventsHtt", 1,0.5,1.5);

				// HH
				NEventsHH2D[signalPair] = new TH1F(Form("NEventsHH%d%d", mchi, ctau), "NEventsHH", 1,0.5,1.5);
				NEventsHHbb2D[signalPair] = new TH1F(Form("NEventsHHbb%d%d", mchi, ctau), "NEventsHHbb", 1,0.5,1.5);
				NEventsHHcc2D[signalPair] = new TH1F(Form("NEventsHHcc%d%d", mchi, ctau), "NEventsHHcc", 1,0.5,1.5);
				NEventsHHgg2D[signalPair] = new TH1F(Form("NEventsHHgg%d%d", mchi, ctau), "NEventsHHgg", 1,0.5,1.5);
				NEventsHHWW2D[signalPair] = new TH1F(Form("NEventsHHWW%d%d", mchi, ctau), "NEventsHHWW", 1,0.5,1.5);
				NEventsHHZZ2D[signalPair] = new TH1F(Form("NEventsHHZZ%d%d", mchi, ctau), "NEventsHHZZ", 1,0.5,1.5);
				NEventsHHmm2D[signalPair] = new TH1F(Form("NEventsHHmm%d%d", mchi, ctau), "NEventsHHmm", 1,0.5,1.5);
				NEventsHHtt2D[signalPair] = new TH1F(Form("NEventsHHtt%d%d", mchi, ctau), "NEventsHHtt", 1,0.5,1.5);

				//HZ
				NEventsZH2D[signalPair] = new TH1F(Form("NEventsZH%d%d", mchi, ctau), "NEventsZH", 1,0.5,1.5);
				NEventsZHbb2D[signalPair] = new TH1F(Form("NEventsZHbb%d%d", mchi, ctau), "NEventsZHbb", 1,0.5,1.5);
				NEventsZHcc2D[signalPair] = new TH1F(Form("NEventsZHcc%d%d", mchi, ctau), "NEventsZHcc", 1,0.5,1.5);
				NEventsZHgg2D[signalPair] = new TH1F(Form("NEventsZHgg%d%d", mchi, ctau), "NEventsZHgg", 1,0.5,1.5);
				NEventsZHWW2D[signalPair] = new TH1F(Form("NEventsZHWW%d%d", mchi, ctau), "NEventsZHWW", 1,0.5,1.5);
				NEventsZHZZ2D[signalPair] = new TH1F(Form("NEventsZHZZ%d%d", mchi, ctau), "NEventsZHZZ", 1,0.5,1.5);
				NEventsZHmm2D[signalPair] = new TH1F(Form("NEventsZHmm%d%d", mchi, ctau), "NEventsZHmm", 1,0.5,1.5);
				NEventsZHtt2D[signalPair] = new TH1F(Form("NEventsZHtt%d%d", mchi, ctau), "NEventsZHtt", 1,0.5,1.5);

				smsSumWeights2D[signalPair] = new TH1D(Form("SumWeights%d%d", mchi, ctau), "SumWeights", 1 ,-0.5,0.5);
				smsSumScaleWeights2D[signalPair] = new TH1D(Form("SumScaleWeights%d%d", mchi, ctau), "SumScaleWeights", 9 ,-0.5,8.5);
				smsSumPdfWeights2D[signalPair] = new TH1D(Form("SumPdfWeights%d%d", mchi, ctau), "SumPdfWeights", 100 ,-0.5,99.5);
				smsSumAlphasWeights2D[signalPair] = new TH1D(Form("SumAlphasWeights%d%d", mchi, ctau), "SumAlphasWeights", 2 ,-0.5,1.5);

				cout << "Created new output file " << thisFileName << endl;
			}
			//Fill NEvents hist
			NEvents2D[signalPair]->Fill(1.0, genWeight);
			if(metType1Pt>200.) NEventsMet2002D[signalPair]->Fill(1.0, genWeight);

			// H decay
			if((gLLP_daughter_id[1]==25)||(gLLP_daughter_id[3]==25)) NEventsH2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==5 && abs(gLLP_grandaughter_id[1])==5)||(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==5 && abs(gLLP_grandaughter_id[3])==5)) NEventsHbb2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==4 && abs(gLLP_grandaughter_id[1])==4)||(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==4 && abs(gLLP_grandaughter_id[3])==4)) NEventsHcc2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==21 && abs(gLLP_grandaughter_id[1])==21)||(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==21 && abs(gLLP_grandaughter_id[3])==21)) NEventsHgg2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==23 && abs(gLLP_grandaughter_id[1])==23)||(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==23 && abs(gLLP_grandaughter_id[3])==23)) NEventsHZZ2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==24 && abs(gLLP_grandaughter_id[1])==24)||(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==24 && abs(gLLP_grandaughter_id[3])==24)) NEventsHWW2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==13 && abs(gLLP_grandaughter_id[1])==13)||(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==13 && abs(gLLP_grandaughter_id[3])==13)) NEventsHmm2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==15 && abs(gLLP_grandaughter_id[1])==15)||(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==15 && abs(gLLP_grandaughter_id[3])==15)) NEventsHtt2D[signalPair]->Fill(1.0, genWeight);

			//HH
			if((gLLP_daughter_id[1]==25)&&(gLLP_daughter_id[3]==25)) NEventsHH2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==5 && abs(gLLP_grandaughter_id[1])==5)&&(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==5 && abs(gLLP_grandaughter_id[3])==5)) NEventsHHbb2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==4 && abs(gLLP_grandaughter_id[1])==4)&&(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==4 && abs(gLLP_grandaughter_id[3])==4)) NEventsHHcc2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==21 && abs(gLLP_grandaughter_id[1])==21)&&(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==21 && abs(gLLP_grandaughter_id[3])==21)) NEventsHHgg2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==23 && abs(gLLP_grandaughter_id[1])==23)&&(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==23 && abs(gLLP_grandaughter_id[3])==23)) NEventsHHZZ2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==24 && abs(gLLP_grandaughter_id[1])==24)&&(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==24 && abs(gLLP_grandaughter_id[3])==24)) NEventsHHWW2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==13 && abs(gLLP_grandaughter_id[1])==13)&&(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==13 && abs(gLLP_grandaughter_id[3])==13)) NEventsHHmm2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==15 && abs(gLLP_grandaughter_id[1])==15)&&(gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==15 && abs(gLLP_grandaughter_id[3])==15)) NEventsHHtt2D[signalPair]->Fill(1.0, genWeight);

			//HZ
			if((gLLP_daughter_id[1]==25 && gLLP_daughter_id[3]==23)||(gLLP_daughter_id[1]==23&&gLLP_daughter_id[3]==25)) NEventsZH2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==5 && abs(gLLP_grandaughter_id[1])==5 && gLLP_daughter_id[3]==23)||(gLLP_daughter_id[1]==23&&gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==5 && abs(gLLP_grandaughter_id[3])==5)) NEventsZHbb2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==4 && abs(gLLP_grandaughter_id[1])==4 && gLLP_daughter_id[3]==23)||(gLLP_daughter_id[1]==23&&gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==4 && abs(gLLP_grandaughter_id[3])==4)) NEventsZHcc2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==21 && abs(gLLP_grandaughter_id[1])==21 && gLLP_daughter_id[3]==23)||(gLLP_daughter_id[1]==23&&gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==21 && abs(gLLP_grandaughter_id[3])==21)) NEventsZHgg2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==24 && abs(gLLP_grandaughter_id[1])==24 && gLLP_daughter_id[3]==23)||(gLLP_daughter_id[1]==23&&gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==24 && abs(gLLP_grandaughter_id[3])==24)) NEventsZHWW2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==23 && abs(gLLP_grandaughter_id[1])==23 && gLLP_daughter_id[3]==23)||(gLLP_daughter_id[1]==23&&gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==23 && abs(gLLP_grandaughter_id[3])==23)) NEventsZHZZ2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==13 && abs(gLLP_grandaughter_id[1])==13 && gLLP_daughter_id[3]==23)||(gLLP_daughter_id[1]==23&&gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==13 && abs(gLLP_grandaughter_id[3])==13)) NEventsZHmm2D[signalPair]->Fill(1.0, genWeight);
			if((gLLP_daughter_id[1]==25 && abs(gLLP_grandaughter_id[0])==15 && abs(gLLP_grandaughter_id[1])==15 && gLLP_daughter_id[3]==23)||(gLLP_daughter_id[1]==23&&gLLP_daughter_id[3]==25 && abs(gLLP_grandaughter_id[2])==15 && abs(gLLP_grandaughter_id[3])==15)) NEventsZHtt2D[signalPair]->Fill(1.0, genWeight);

			smsSumWeights2D[signalPair]->Fill(1.0, weight);

			smsSumScaleWeights2D[signalPair]->Fill(0.0, llp_tree->sf_facScaleUp);
			smsSumScaleWeights2D[signalPair]->Fill(1.0, llp_tree->sf_facScaleDown);
			smsSumScaleWeights2D[signalPair]->Fill(2.0, llp_tree->sf_renScaleUp);
			smsSumScaleWeights2D[signalPair]->Fill(3.0, llp_tree->sf_renScaleDown);
			smsSumScaleWeights2D[signalPair]->Fill(4.0, llp_tree->sf_facRenScaleUp);
			smsSumScaleWeights2D[signalPair]->Fill(5.0, llp_tree->sf_facRenScaleDown);
			//it should goes to 8? or we only care about first 6

			for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) {
				smsSumPdfWeights2D[signalPair]->Fill(double(iwgt),(*pdfWeights)[iwgt]);
			}

			if(_debug_pdf) cout<<"alphas size: "<< alphasWeights->size() <<endl;
			for (unsigned int iwgt=0; iwgt<alphasWeights->size(); ++iwgt) {
				smsSumAlphasWeights2D[signalPair]->Fill(double(iwgt),(*alphasWeights)[iwgt]);
			}




		}

		if (label =="bkg_wH"|| label == "bkg_zH" || label == "bkg_HH" ){
		//if (label =="bkg_wH"|| label == "bkg_zH" || label == "bkg_HH" || label == "HH"){
			if (isData)
			{
				NEvents->Fill(1);
				if(metType1Pt>200.) NEventsMet200->Fill(1);
				llp_tree->weight = 1;
			}
			else
			{
				//NEvents->Fill(genWeight);
				//llp_tree->weight = genWeight;
				NEvents->Fill(1);
				if(metType1Pt>200.) NEventsMet200->Fill(1);
				llp_tree->weight = 1;

				//signal not mass scan
				//NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
      				weight = genWeight;
      				SumWeights->Fill(1.0, weight);
				if ( (*scaleWeights).size() >= 9 ) 
				{
					llp_tree->sf_facScaleUp      = (*scaleWeights)[1]/genWeight;
					llp_tree->sf_facScaleDown    = (*scaleWeights)[2]/genWeight;
					llp_tree->sf_renScaleUp      = (*scaleWeights)[3]/genWeight;
					llp_tree->sf_renScaleDown    = (*scaleWeights)[6]/genWeight;
					llp_tree->sf_facRenScaleUp   = (*scaleWeights)[4]/genWeight;
					llp_tree->sf_facRenScaleDown = (*scaleWeights)[8]/genWeight;


					SumScaleWeights->Fill(0.0, (*scaleWeights)[1]);
					SumScaleWeights->Fill(1.0, (*scaleWeights)[2]);
					SumScaleWeights->Fill(2.0, (*scaleWeights)[3]);
					SumScaleWeights->Fill(3.0, (*scaleWeights)[6]);
					SumScaleWeights->Fill(4.0, (*scaleWeights)[4]);
					SumScaleWeights->Fill(5.0, (*scaleWeights)[8]);
				}

				llp_tree->sf_pdf.erase( llp_tree->sf_pdf.begin(), llp_tree->sf_pdf.end() );
				if(_debug_pdf) cout << "pdf weight size: " << pdfWeights->size() << endl;
				if(_debug_pdf) cout << "scale weight size: " << scaleWeights->size() << endl;
				if(_debug_pdf) cout << "alphas weight size: " << alphasWeights->size() << endl;
				for ( unsigned int iwgt = 0; iwgt < pdfWeights->size(); ++iwgt ) 
				{
					llp_tree->sf_pdf.push_back( pdfWeights->at(iwgt)/genWeight );
					SumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
					if(_debug_pdf) cout << "pdf weight component: " << (*pdfWeights)[iwgt] << endl;
				}

				for ( unsigned int iwgt = 0; iwgt < alphasWeights->size(); ++iwgt ) 
				{
					SumAlphasWeights->Fill(double(iwgt),(*alphasWeights)[iwgt]);
				}

			}
		}
		//else if(label =="MR_EMU"|| label == "MR_PHO" || label == "MR_ZLL")
		//
		else if( (label.find("MR") != std::string::npos) ){
			NEvents->Fill(1);
			if(metType1Pt>200.) NEventsMet200->Fill(1);
		}
		else{
			//generatedEvents->Fill(1);
			llp_tree->weight = 1;
		}
		if(_debug) std::cout << "deb2 " << jentry << std::endl;
		if(!isData && analysisTag=="CT2018_17SeptEarlyReReco"){
			//if( eventNum% 2 ==0 ) continue;
		}
		//event info
		llp_tree->runNum = runNum;
		llp_tree->lumiSec = lumiNum;
		llp_tree->evtNum = eventNum;
		llp_tree->pvX = pvX;
		llp_tree->pvY = pvY;
		llp_tree->pvZ = pvZ;
		//PVAll
		llp_tree->npvall = nPVAll;
		//std::cout << nPVAll <<endl;
		//std::cout << llp_tree->npvall <<endl;
		int pvid0 = -1;
		float pvD2 = pvX*pvX+pvY*pvY+pvZ*pvZ;
		float pvAllD2 = 0;
		float delta_pvD2 = 999.;
		float min_delta_pvD2 = 999.;
		for(int i=0;i<nPVAll;i++)
		{
			llp_tree->pvAllX[i] = pvAllX[i];
			llp_tree->pvAllY[i] = pvAllY[i];
			llp_tree->pvAllZ[i] = pvAllZ[i];
			pvAllD2 = pvAllX[i]*pvAllX[i]+pvAllY[i]*pvAllY[i]+pvAllZ[i]*pvAllZ[i];
			delta_pvD2 = abs(pvD2-pvAllD2);
			if(delta_pvD2<min_delta_pvD2)
			{
				min_delta_pvD2 = delta_pvD2;
				pvid0 = i;
			}
		}
		//std::cout << "pvid0" <<pvid0 <<endl;
		//std::cout << "min delta pv d2" <<min_delta_pvD2 <<endl;
		llp_tree->pv_index = pvid0;
		if(_debug) std::cout << "deb3 " << jentry << std::endl;
		if(_debug) std::cout << "nBunchXing " << nBunchXing << std::endl;
		if (label == "zH" || label == "wH" ){
			NEvents->Fill(1);
			if(metType1Pt>200.) NEventsMet200->Fill(1);
			bool wzFlag = false;
			for (int i=0; i < nGenParticle; ++i)
			{
				// if (abs(gParticleId[i]) == wzId && gParticleStatus[i] == 22)
				//if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == wzId)
				if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] && abs(gParticleMotherId[i]) == wzId)
				{
					wzFlag = true;
				}

			}
			if ( wzFlag == false ) continue;

		}
		else if (label == "HH"){
			NEvents->Fill(1);
			if(metType1Pt>200.) NEventsMet200->Fill(1);
				llp_tree->weight = 1;

				//signal not mass scan
				//NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
      				weight = genWeight;
      				SumWeights->Fill(1.0, weight);
				if ( (*scaleWeights).size() >= 9 ) 
				{
					llp_tree->sf_facScaleUp      = (*scaleWeights)[1]/genWeight;
					llp_tree->sf_facScaleDown    = (*scaleWeights)[2]/genWeight;
					llp_tree->sf_renScaleUp      = (*scaleWeights)[3]/genWeight;
					llp_tree->sf_renScaleDown    = (*scaleWeights)[6]/genWeight;
					llp_tree->sf_facRenScaleUp   = (*scaleWeights)[4]/genWeight;
					llp_tree->sf_facRenScaleDown = (*scaleWeights)[8]/genWeight;


					SumScaleWeights->Fill(0.0, (*scaleWeights)[1]);
					SumScaleWeights->Fill(1.0, (*scaleWeights)[2]);
					SumScaleWeights->Fill(2.0, (*scaleWeights)[3]);
					SumScaleWeights->Fill(3.0, (*scaleWeights)[6]);
					SumScaleWeights->Fill(4.0, (*scaleWeights)[4]);
					SumScaleWeights->Fill(5.0, (*scaleWeights)[8]);
				}

				llp_tree->sf_pdf.erase( llp_tree->sf_pdf.begin(), llp_tree->sf_pdf.end() );
				if(_debug_pdf) cout << "pdf weight size: " << pdfWeights->size() << endl;
				if(_debug_pdf) cout << "scale weight size: " << scaleWeights->size() << endl;
				if(_debug_pdf) cout << "alphas weight size: " << alphasWeights->size() << endl;
				for ( unsigned int iwgt = 0; iwgt < pdfWeights->size(); ++iwgt ) 
				{
					llp_tree->sf_pdf.push_back( pdfWeights->at(iwgt)/genWeight );
					SumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
					if(_debug_pdf) cout << "pdf weight component: " << (*pdfWeights)[iwgt] << endl;
				}

				for ( unsigned int iwgt = 0; iwgt < alphasWeights->size(); ++iwgt ) 
				{
					SumAlphasWeights->Fill(double(iwgt),(*alphasWeights)[iwgt]);
				}

			bool wzFlag = false;
			for (int i=0; i < nGenParticle; ++i)
			{
				if ((abs(gParticleId[i]) == 5) && gParticleStatus[i] && abs(gParticleMotherId[i]) == wzId)
				{
					wzFlag = true;
				}

			}
			if(_debug) std::cout << "wzFlag " << wzFlag << std::endl;
			//if ( wzFlag == false ) continue;

		}

		if(!isData)
		{

			//   for(int i=0; i < nGenJets; i++)
			//   {
			//   llp_tree->genJetE[llp_tree->nGenJets] = genJetE[i];
			//   llp_tree->genJetPt[llp_tree->nGenJets] = genJetPt[i];
			//   llp_tree->genJetEta[llp_tree->nGenJets] = genJetEta[i];
			//   llp_tree->genJetPhi[llp_tree->nGenJets] = genJetPhi[i];
			//// llp_tree->nGenJets++;
			//}

			llp_tree->genMetPtTrue = genMetPtTrue;
			llp_tree->genMetPhiTrue = genMetPhiTrue;
			llp_tree->genMetPtCalo = genMetPtCalo;
			llp_tree->genMetPhiCalo = genMetPhiCalo;


			for(int i = 0; i < 2;i++)
			{
				//std::cout << "i " << i << std::endl;
				//std::cout << "gLLP_beta[i] " << gLLP_beta[i] << std::endl;
				//std::cout << "gLLP_eta[i] " << gLLP_eta[i] << std::endl;
				//std::cout << "gLLP_phi[i] " << gLLP_phi[i] << std::endl;
				//std::cout << "gLLP_decay_vertex_x[i] " << gLLP_decay_vertex_x[i] << std::endl;
				//std::cout << "gLLP_decay_vertex_y[i] " << gLLP_decay_vertex_y[i] << std::endl;
				//std::cout << "gLLP_decay_vertex_z[i] " << gLLP_decay_vertex_z[i] << std::endl;
				llp_tree->gLLP_e[i] = gLLP_e[i];
				llp_tree->gLLP_pt[i] = gLLP_pt[i];
				llp_tree->gLLP_eta[i] = gLLP_eta[i];
				llp_tree->gLLP_phi[i] = gLLP_phi[i];
				llp_tree->gLLP_travel_time[i] = gLLP_travel_time[i];
				llp_tree->gLLP_decay_vertex_r[i] = sqrt(gLLP_decay_vertex_x[i]*gLLP_decay_vertex_x[i]+gLLP_decay_vertex_y[i]*gLLP_decay_vertex_y[i]);
				//std::cout << "llp_tree->gLLP_decay_vertex_r[i] " << llp_tree->gLLP_decay_vertex_r[i] << std::endl;
				llp_tree->gLLP_decay_vertex_x[i] = gLLP_decay_vertex_x[i];
				llp_tree->gLLP_decay_vertex_y[i] = gLLP_decay_vertex_y[i];
				llp_tree->gLLP_decay_vertex_z[i] = gLLP_decay_vertex_z[i];
				//acceptance cut
				if(llp_tree->gLLP_decay_vertex_r[i] > 30. && llp_tree->gLLP_decay_vertex_r[i] < 184. && abs(llp_tree->gLLP_decay_vertex_z[i]) < 376. && abs(llp_tree->gLLP_eta[i])!=666) {
					if(i<1) llp_tree->gLLP0_EB = true;
					else llp_tree->gLLP1_EB = true;
				}

				float beta = gLLP_beta[i];
				if(beta<0){
					TLorentzVector gLLP = makeTLorentzVector( gLLP_pt[i], gLLP_eta[i], gLLP_phi[i], gLLP_e[i] );

					//std::cout << "gLLP.Px " << gLLP.Px() << std::endl;
					//std::cout << "gLLP.Py " << gLLP.Py() << std::endl;
					//std::cout << "gLLP.Pz " << gLLP.Pz() << std::endl;
					float Px = gLLP.Px();
					float Py = gLLP.Py();
					float Pz = gLLP.Pz();

					//float P2 = sqrt(Px*Px+Py*Py+Pz*Pz);
					//std::cout << "P2 " << P2 << std::endl;
					//std::cout << "gLLP_e[i] " << gLLP_e[i] << std::endl;
					beta = sqrt(Px*Px+Py*Py+Pz*Pz)/gLLP_e[i];
				}
				float gLLP_decay_vertex = sqrt(pow(llp_tree->gLLP_decay_vertex_r[i], 2) + pow(llp_tree->gLLP_decay_vertex_z[i],2));
				//std::cout << "gLLP_decay_vertex " << gLLP_decay_vertex << std::endl;
				float gamma = 1.0/sqrt(1-beta*beta);
				//std::cout << "beta " << beta << std::endl;
				//std::cout << "gamma " << gamma << std::endl;
				llp_tree->gLLP_ctau[i] = gLLP_decay_vertex/(beta * gamma);
				llp_tree->gLLP_beta[i] = gLLP_beta[i];

				//if (abs(llp_tree->gLLP_eta[i]) < 2.4 && abs(llp_tree->gLLP_eta[i]) > 0.9
				//  && abs(llp_tree->gLLP_decay_vertex_z[i])<1100 && abs(llp_tree->gLLP_decay_vertex_z[i])>568
				//  && llp_tree->gLLP_decay_vertex_r[i] < 695.5) llp_tree->gLLP_csc[i] = true;

				if( abs(llp_tree->gLLP_decay_vertex_z[i])<268.3 && llp_tree->gLLP_decay_vertex_r[i] < 129.0 && abs(llp_tree->gLLP_eta[i]) < 1.48 ) llp_tree->gLLP_eb[i] = true;

			}//gLLP 




			if(_debug) std::cout << "nBunchXing " << nBunchXing << std::endl;
			for (int i=0; i < nBunchXing; i++)
			{
				if (BunchXing[i] == 0)
				{
					llp_tree->npu = int(nPUmean[i]);
					if(_debug_npu) 
					{
						std::cout << "npu " << llp_tree->npu << std::endl;
						std::cout << "nPUmean[i] " << nPUmean[i] << std::endl;
					}
				}
			}

			llp_tree->pileupWeight = 1;
			llp_tree->pileupWeight = helper->getPileupWeight(llp_tree->npu);
			llp_tree->pileupWeightUp = helper->getPileupWeightUp(llp_tree->npu) / llp_tree->pileupWeight;
			llp_tree->pileupWeightDown = helper->getPileupWeightDown(llp_tree->npu) / llp_tree->pileupWeight;
			if(_debug_npu) 
			{
				std::cout << "pileupWeightUp " << llp_tree->pileupWeightUp << std::endl;
				std::cout << "pileupWeightDown " << llp_tree->pileupWeightDown << std::endl;
			}
		}

		if(_debug_met) std::cout << "npu " << llp_tree->npu << std::endl;
		if(_debug && llp_tree->npu != 0 ) std::cout << "npu " << llp_tree->npu << std::endl;
		//get NPU
		llp_tree->npv = nPV;
		llp_tree->rho = fixedGridRhoFastjetAll;
		llp_tree->met = metType1Pt;
		if(_debug_met) std::cout << "met " << llp_tree->met << std::endl;
		//if( llp_tree->met < 120. ) continue;
		//if( llp_tree->met < 150. ) continue;
		if(_debug_lab) std::cout << "label " << label.c_str() << std::endl;
		if(_debug_lab) std::cout << "met " << llp_tree->met << std::endl;
		if( (label.find("MR") == std::string::npos) && llp_tree->met < 200. ) continue;
		//if( presel && (label.find("MR") == std::string::npos) && llp_tree->met < 200. ) continue;
		if( (label=="MR_EMU") && llp_tree->met < 30. ) continue;
		if( (label.find("MR_Single") != std::string::npos) && llp_tree->met < 40. ) continue;
		if( (label.find("MR_ZLL") != std::string::npos) && isData && llp_tree->met >= 30. ) continue;
		if( (label.find("MR_JetHT") != std::string::npos) && llp_tree->met >= 30. ) continue;
		if( (label=="MR_PHO") && llp_tree->met >= 30. ) continue;
		if(_debug_lab) std::cout << "label " << label.c_str() << "passed "<< std::endl;
		if(_debug_lab) std::cout << "met " << llp_tree->met << "passed "<< std::endl;
		if(_debug_met) std::cout << "metType1Pt passed" << metType1Pt << std::endl;
		llp_tree->metPhi = metType1Phi;
		if(_debug) std::cout << "npv " << llp_tree->npv << std::endl;
		TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM( metType1Pt, 0, metType1Phi, 0 );

		if(_debug) std::cout << "isData " << isData << std::endl;
		//Znunu
		if(!isData && (label.find("MR_ZLL") != std::string::npos) )
		{
			if(_debug) std::cout << "before nu cuts " << std::endl;
			int count_nu=0;
			int count_nu25=0;
			if(_debug) std::cout << "nGenParticle " << nGenParticle << std::endl;
			for(int i = 0; i < nGenParticle; i++)
			{
				//(abs(gParticleId)==12||abs(gParticleId)==14||abs(gParticleId)==16)&&gParticleMotherId==23
				if ( !( abs(gParticleId[i])==12 || abs(gParticleId[i])==14 || abs(gParticleId[i])==16 ) ) continue;
				if (abs(gParticleMotherId[i])!=23 ) continue;
				if (abs(gParticleStatus[i])!=1) continue;
				count_nu++;
				if(gParticlePt[i]<25) continue;
				count_nu25++;
			}
			if(_debug) std::cout << "count_nu " << count_nu << std::endl;
			if(_debug) std::cout << "count_nu25 " << count_nu25 << std::endl;
			if(count_nu!=2) continue;
			if(count_nu25!=2) continue;
			if(_debug) std::cout << "passed nu cuts " << std::endl;
		}

		//met filters
		llp_tree->Flag2_globalSuperTightHalo2016Filter          = Flag2_globalSuperTightHalo2016Filter;
		llp_tree->Flag2_globalTightHalo2016Filter               = Flag2_globalTightHalo2016Filter;
		llp_tree->Flag2_goodVertices                            = Flag2_goodVertices;
		llp_tree->Flag2_BadChargedCandidateFilter               = Flag2_BadChargedCandidateFilter;
		llp_tree->Flag2_BadPFMuonFilter                         = Flag2_BadPFMuonFilter;
		llp_tree->Flag2_EcalDeadCellTriggerPrimitiveFilter      = Flag2_EcalDeadCellTriggerPrimitiveFilter;
		llp_tree->Flag2_HBHENoiseFilter                         = Flag2_HBHENoiseFilter;
		llp_tree->Flag2_HBHEIsoNoiseFilter                      = Flag2_HBHEIsoNoiseFilter;
		llp_tree->Flag2_ecalBadCalibFilter                      = Flag2_ecalBadCalibFilter;
		llp_tree->Flag2_eeBadScFilter                           = Flag2_eeBadScFilter;

		if(presel && !(Flag2_globalSuperTightHalo2016Filter && Flag2_BadPFMuonFilter && Flag2_EcalDeadCellTriggerPrimitiveFilter && Flag2_HBHENoiseFilter && Flag2_HBHEIsoNoiseFilter )) continue;
		if(presel &&  (analysisTag.find("2016") == std::string::npos) && !Flag2_ecalBadCalibFilter) continue;
		if(presel && isData && !Flag2_eeBadScFilter) continue;
		//Triggers
		for(int i = 0; i < NTriggersMAX; i++){
			llp_tree->HLTDecision[i] = HLTDecision[i];
			llp_tree->HLTPrescale[i] = HLTPrescale[i];
		}
		if(_debug_trg) std::cout << "begin: 310 " << HLTDecision[310] << std::endl;
		if(_debug_trg) std::cout << "begin: 310 " << llp_tree->HLTDecision[310] << std::endl;
		//if( (label.find("MR") == std::string::npos) && !HLTDecision[310]) continue; 
		//if( presel && (label.find("MR") == std::string::npos) && !HLTDecision[467]) continue; 
		//analysisTag: l=CT2016_07Aug2017Rereco  l=CT2017_17Nov2017Rereco l=CT2018_17SeptEarlyReReco 
		//trigger: 2016:  467; 2017: 467 || 717 || 710; 2018: 467 || 717 || 710
		if( presel && (label.find("MR") == std::string::npos) && analysisTag=="CT2016_07Aug2017Rereco" && !HLTDecision[467]) continue; 
		if( presel && (label.find("MR") == std::string::npos) && analysisTag=="CT2017_17Nov2017Rereco" && !(HLTDecision[467]||HLTDecision[717]||HLTDecision[710]) ) continue; 
		if( presel && (label.find("MR") == std::string::npos) && analysisTag=="CT2018_17SeptEarlyReReco" && !(HLTDecision[467]||HLTDecision[717]||HLTDecision[710]) ) continue; 

		if(_debug_lab) std::cout << "label " << label.c_str() << std::endl;
		bool triggered = false;
		if( (label.find("MR_PHO") != std::string::npos) || (label.find("MR_JetHT") != std::string::npos) ){
			if(_debug_lab) std::cout << "label " << label.c_str() << std::endl;
			triggered = false;
			for(int i = 0; i < NTrigger; i++){
				int trigger_temp = trigger_paths[i];
				if(_debug_pre) std::cout << "trig " << trigger_temp << std::endl;
				if(_debug_pre) std::cout << "trig dec " << HLTDecision[trigger_temp] << std::endl;
				if(_debug_pre) std::cout << "trig pre " << HLTPrescale[trigger_temp] << std::endl;
				triggered = triggered || HLTDecision[trigger_temp];
				if(_debug_pre) std::cout << "triggered " << triggered << std::endl;
				if(triggered){
					//assign weight based on prescale of the highest threshold passed trigger
					llp_tree->weight = 1*HLTPrescale[trigger_temp];
					if(_debug_pre) std::cout << "weight " << llp_tree->weight << std::endl;
					break;
				}
				else{
					llp_tree->weight = 1;
					if(_debug_pre) std::cout << "weight " << llp_tree->weight << std::endl;
				}
			}//end for loop
		}//MR_PHO/JetHT re-weight based on Prescale
		else if( (label.find("MR_EMU") != std::string::npos) || (label.find("MR_Single") != std::string::npos) || (label=="MR_ZLL")){
			if(_debug_lab) std::cout << "label " << label.c_str() << std::endl;
			triggered = false;
			int tempW = 999999999;
			for(int i = 0; i < NTrigger; i++){
				int trigger_temp = trigger_paths[i];
				if(_debug_pre) std::cout << "trig " << trigger_temp << std::endl;
				if(_debug_pre) std::cout << "trig dec " << HLTDecision[trigger_temp] << std::endl;
				if(_debug_pre) std::cout << "trig pre " << HLTPrescale[trigger_temp] << std::endl;
				triggered = triggered || HLTDecision[trigger_temp];
				if(_debug_pre) std::cout << "triggered " << triggered << std::endl;
				if(triggered && tempW!=1){
					//if weight is already 1, no need to compare again
					//assign weight based on smallest prescale of passed trigger
					if(_debug_pre) std::cout << "tempW " << tempW << std::endl;
					if(tempW>=HLTPrescale[trigger_temp] && HLTPrescale[trigger_temp]>0){
						//make sure trigger is not disabled
						tempW = HLTPrescale[trigger_temp];
						if(_debug_pre) std::cout << "tempW " << tempW << std::endl;
					}
				}
			}//end for loop
			if(triggered){
				llp_tree->weight = tempW;
				if(_debug_pre) std::cout << "weight " << llp_tree->weight << std::endl;
			}
			else{
				llp_tree->weight = 1;
				if(_debug_pre) std::cout << "weight " << llp_tree->weight << std::endl;
			}
		}//MR_EMU/SingleMuon/SingleElectron/ZLL re-weight based on Prescale	

		if( (label.find("MR_ZLL") != std::string::npos) && !isData) triggered=true; 
		//if( (label.find("MR") != std::string::npos) && !triggered) continue; 

		if(_debug) std::cout << "passed trigger " << std::endl;

		//bool triggered = false;
		//for(int i = 0; i < NTrigger; i++)
		//{
		//	if(_debug_trg) std::cout << "i " << i << ", NTrigger "<< NTrigger << std::endl;

		// int trigger_temp = trigger_paths[i];
		//	if(_debug_trg) std::cout << "temp  " << trigger_paths[i] << ", triggered "<< triggered << std::endl;

		//	triggered = triggered || HLTDecision[trigger_temp];
		//	if(_debug_trg) std::cout << "is triggered ?"<< triggered << std::endl;

		//}
		//if (triggered) trig->Fill(1);
		//if(_debug) std::cout << "triggered " << triggered << std::endl;


		// ************************************************************************
		//Start Object Selection
		// ************************************************************************
		//sync 
		if(_debug_sync)
		{
			std::cout << "nMuons " << nMuons << std::endl;
			std::cout << "nElectron " << nElectrons << std::endl;
			std::cout << "nPhoton " << nPhotons << std::endl;
			std::cout << "nTau " << nTaus << std::endl;
		}


		std::vector<leptons> Leptons;
		//-------------------------------
		//Muons
		//-------------------------------
		//twiki (https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation)
		for( int i = 0; i < nMuons; i++ )
		{
			if(_debug_sync)
			{
				std::cout << "nMuons " << nMuons << std::endl;
				std::cout << "iMuon " << i << ", Pt " << muonPt[i] << ", Eta " << muonEta[i]<< ", Phi " << muonPhi[i]<< ", E " << muonE[i] << std::endl;
				std::cout << "iMuon " << i << ", muon_chargedIso[i] " << muon_chargedIso[i] << ", muon_neutralHadIso[i] " << muon_neutralHadIso[i]<< ", muon_photonIso[i] " << muon_photonIso[i]<< ", muon_pileupIso[i] " << muon_pileupIso[i] << std::endl;
				std::cout << "iMuon " << i << ", muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i] " << muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i] << ", std::max( stuff, 0)  " << std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) << ", muon_chargedIso[i] + std::max( stuff, 0) " << (muon_chargedIso[i] + std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) ) << ", (muon_chargedIso[i] + std::max( stuff, 0))/muonPt " << (muon_chargedIso[i] + std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) )/muonPt[i] << std::endl;
			}

			if(muonPt[i] <= 10 ) continue;
			//if( (label=="MR_EMU") && muonPt[i] < 25. ) continue;
			if( (label.find("MR") != std::string::npos) && muonPt[i] < 25.) continue; 
			if(fabs(muonEta[i]) > 2.4) continue;

			if(!muonIsLoose[i]) continue;
			if( (muon_chargedIso[i] + std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) )/muonPt[i] >= 0.25) continue;

			if(xclean){
				//remove overlaps
				bool overlap = false;
				for(auto& lep : Leptons)
				{
					if (RazorAnalyzerLLP::deltaR(muonEta[i],muonPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
				}
				if(overlap) continue;
			}

			leptons tmpMuon;
			tmpMuon.lepton.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i], MU_MASS);
			tmpMuon.pdgId = 13 * -1 * muonCharge[i];

			Leptons.push_back(tmpMuon);

			llp_tree->muon_pileupIso[llp_tree->nMuons] = muon_pileupIso[i];
			llp_tree->muon_chargedIso[llp_tree->nMuons] = muon_chargedIso[i];
			llp_tree->muon_photonIso[llp_tree->nMuons] = muon_photonIso[i];
			llp_tree->muon_neutralHadIso[llp_tree->nMuons] = muon_neutralHadIso[i];

			llp_tree->muonIsLoose[llp_tree->nMuons] = muonIsLoose[i];

			llp_tree->muonPt[llp_tree->nMuons] = muonPt[i];
			llp_tree->muonEta[llp_tree->nMuons] = muonEta[i];
			llp_tree->muonE[llp_tree->nMuons] = muonE[i];
			llp_tree->muonPhi[llp_tree->nMuons] = muonPhi[i];

			llp_tree->nMuons++;
		}
		if( (label=="MR_SingleMuon") && llp_tree->nMuons != 1 ) continue;
		if(_debug) std::cout << "nMuons " << llp_tree->nMuons << std::endl;

		//-------------------------------
		//Electrons
		//-------------------------------
		for( int i = 0; i < nElectrons; i++ )
		{
			if(_debug_sync)
			{
				std::cout << "nElectrons " << nElectrons << std::endl;
				std::cout << "iElectron " << i << ", Pt " << elePt[i] << ", Eta " << eleEta[i]<< ", Phi " << elePhi[i]<< ", E " << eleE[i] << std::endl;
				std::cout << "iElectron " << i << ", ele_passCutBasedIDVeto[i] " << ele_passCutBasedIDVeto[i]<< std::endl;
			}
			///
			//			if(eventNum==38906||eventNum==71362||eventNum==75125||eventNum==3877||eventNum==20903)
			//			{
			//				std::cout << "eventNum " << eventNum << std::endl;	
			//				std::cout << "elePt[i] " << elePt[i] << std::endl;	
			//				std::cout << "eleEta[i] " << eleEta[i] << std::endl;	
			//				std::cout << "ele_passCutBasedIDVeto[i] " << ele_passCutBasedIDVeto[i] << std::endl;	
			//			}
			/////
			if(elePt[i] <= 10 ) continue;
			//if( (label=="MR_EMU") && elePt[i] < 25. ) continue;
			if( (label.find("MR") != std::string::npos) && elePt[i] < 25.) continue; 
			if(fabs(eleEta[i]) > 2.5) continue;

			if(!ele_passCutBasedIDVeto[i]) continue;

			if(xclean){
				//remove overlaps
				bool overlap = false;
				for(auto& lep : Leptons)
				{
					if (RazorAnalyzerLLP::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
				}
				if(overlap) continue;
			}

			leptons tmpElectron;
			tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);

			Leptons.push_back(tmpElectron);

			llp_tree->ele_passCutBasedIDVeto[llp_tree->nElectrons] = ele_passCutBasedIDVeto[i];

			llp_tree->elePt[llp_tree->nElectrons] = elePt[i];
			llp_tree->eleEta[llp_tree->nElectrons] = eleEta[i];
			llp_tree->eleE[llp_tree->nElectrons] = eleE[i];
			llp_tree->elePhi[llp_tree->nElectrons] = elePhi[i];

			llp_tree->nElectrons++;
		}
		if( (label=="MR_SingleElectron") && llp_tree->nElectrons != 1 ) continue;
		if(_debug) std::cout << "nElectrons " << llp_tree->nElectrons << std::endl;



		std::vector<taus> Taus;
		//-------------------------------
		//Taus
		//-------------------------------
		for( int i = 0; i < nTaus; i++ )
		{
			if(_debug_sync)
			{
				std::cout << "nTaus " << nTaus << std::endl;
				std::cout << "iTau " << i << ", Pt " << tauPt[i] << ", Eta " << tauEta[i]<< ", Phi " << tauPhi[i]<< ", E " << tauE[i] << std::endl;
				std::cout << "iTau " << i << ", tau_IsLoose[i] " << tau_IsLoose[i]<< std::endl;
			}

			//if(tauPt[i] <= 20 ) continue;
			if(tauPt[i] <= 18 ) continue;
			if(fabs(tauEta[i]) > 2.3) continue;

			if(!tau_IsLoose[i]) continue;
			if(tau_ID[i]%2==0) continue;
			//std::cout << "iTau " << i << ", tau_ID[i] " << tau_ID[i]<< std::endl;

			if(xclean){
				//remove overlaps
				bool overlap = false;
				for(auto& lep : Leptons)
				{
					if (RazorAnalyzerLLP::deltaR(tauEta[i],tauPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
				}
				if(overlap) continue;

				bool overlap_tau = false;
				for(auto& tau : Taus)
				{
					if (RazorAnalyzerLLP::deltaR(tauEta[i],tauPhi[i],tau.tau.Eta(),tau.tau.Phi()) < 0.3) overlap_tau = true;
				}
				if(overlap_tau) continue;
			}

			taus tmpTau;
			tmpTau.tau.SetPtEtaPhiM(tauPt[i],tauEta[i], tauPhi[i], TAU_MASS);

			Taus.push_back(tmpTau);

			llp_tree->tau_IsLoose[llp_tree->nTaus] = tau_IsLoose[i];

			llp_tree->tauPt[llp_tree->nTaus] = tauPt[i];
			llp_tree->tauEta[llp_tree->nTaus] = tauEta[i];
			llp_tree->tauE[llp_tree->nTaus] = tauE[i];
			llp_tree->tauPhi[llp_tree->nTaus] = tauPhi[i];
			llp_tree->tau_ID[llp_tree->nTaus] = tau_ID[i];

			llp_tree->nTaus++;
		}
		if(_debug) std::cout << "nTaus " << llp_tree->nTaus << std::endl;

		std::vector<photons> Photons;
		//-------------------------------
		//Photons
		//-------------------------------
		for( int i = 0; i < nPhotons; i++ )
		{
			if(_debug_sync)
			{
				std::cout << "nPhotons " << nPhotons << std::endl;
				std::cout << "iPhoton " << i << ", Pt " << phoPt[i] << ", Eta " << phoEta[i]<< ", Phi " << phoPhi[i]<< ", E " << phoE[i] << std::endl;
				std::cout << "iPhoton " << i << ", pho_passCutBasedIDLoose[i] " << pho_passCutBasedIDLoose[i]<< std::endl;
			}

			if(phoPt[i] <= 15 ) continue;
			if(fabs(phoEta[i]) > 2.5) continue;

			if(!pho_passCutBasedIDLoose[i]) continue;

			if(xclean){
				//remove overlaps
				bool overlap = false;
				for(auto& lep : Leptons)
				{
					if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
				}
				if(overlap) continue;

				bool overlap_tau = false;
				for(auto& tau : Taus)
				{
					if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],tau.tau.Eta(),tau.tau.Phi()) < 0.3) overlap_tau = true;
				}
				if(overlap_tau) continue;

				bool overlap_pho = false;
				for(auto& pho : Photons)
				{
					if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],pho.photon.Eta(),pho.photon.Phi()) < 0.3) overlap_pho = true;
				}
				if(overlap_pho) continue;
			}

			photons tmpPhoton;
			tmpPhoton.photon.SetPtEtaPhiM(phoPt[i],phoEta[i], phoPhi[i], 0.);

			Photons.push_back(tmpPhoton);

			llp_tree->pho_passCutBasedIDLoose[llp_tree->nPhotons] = pho_passCutBasedIDLoose[i];

			llp_tree->phoPt[llp_tree->nPhotons] = phoPt[i];
			llp_tree->phoEta[llp_tree->nPhotons] = phoEta[i];
			llp_tree->phoE[llp_tree->nPhotons] = phoE[i];
			llp_tree->phoPhi[llp_tree->nPhotons] = phoPhi[i];

			if(phoPt[i]>50. && fabs(phoEta[i])>2.25 && fabs(phoEta[i])<3.0) llp_tree->IsPrefiringAffected = true;

			llp_tree->nPhotons++;
		}
		if(_debug) std::cout << "nPhotons " << llp_tree->nPhotons << std::endl;
		if( (label=="MR_PHO") && llp_tree->nPhotons != 1 ) continue;

		//-------------------------------
		//Leptons
		//-------------------------------
		TLorentzVector lepp4;
		sort(Leptons.begin(), Leptons.end(), my_largest_pt_lep);
		for ( auto &tmp : Leptons )
		{
			llp_tree->lepE[llp_tree->nLeptons]      = tmp.lepton.E();
			llp_tree->lepPt[llp_tree->nLeptons]     = tmp.lepton.Pt();
			llp_tree->lepEta[llp_tree->nLeptons]    = tmp.lepton.Eta();
			llp_tree->lepPhi[llp_tree->nLeptons]    = tmp.lepton.Phi();
			llp_tree->lepPdgId[llp_tree->nLeptons]  = tmp.pdgId;
			//llp_tree->lepDZ[llp_tree->nLeptons]     = tmp.dZ;
			//llp_tree->lepPassId[llp_tree->nLeptons] = tmp.passId;
			if(_debug) std::cout << "lepE " << tmp.lepton.E() << std::endl;

			// std::cout << "lepton pdg " << llp_tree->lepPdgId[llp_tree->nLeptons] << std::endl;
			llp_tree->nLeptons++;

			lepp4 += tmp.lepton;
		}
		float dPhi = RazorAnalyzerLLP::deltaPhi(llp_tree->metPhi, lepp4.Phi());
		llp_tree->MT = sqrt(2*(llp_tree->met)*lepp4.Pt()*(1-cos(dPhi)));
		if(_debug_lep) std::cout<<"MT "<<llp_tree->MT<<std::endl;
		//if (triggered) trig_lepId->Fill(1);
		if(_debug) std::cout << "nLeptons " << llp_tree->nLeptons << std::endl;

		//-------------------------------
		// reconstruct Z
		//-------------------------------
		double ZMass = -999;
		double ZPt = -999;
		double tmpDistToZPole = 9999;
		pair<uint,uint> ZCandidateLeptonIndex;
		bool foundZ = false;
		TLorentzVector ZCandidate;
		if(_debug) std::cout << "Leptons.size() " << Leptons.size() << std::endl;
		for( uint i = 0; i < Leptons.size(); i++ )
		{
			for( uint j = 0; j < Leptons.size(); j++ )
			{
				if (!( Leptons[i].pdgId == -1*Leptons[j].pdgId )) continue;// same flavor opposite charge
				double tmpMass = (Leptons[i].lepton+Leptons[j].lepton).M();
				//select the pair closest to Z pole mass
				if ( fabs( tmpMass - Z_MASS) < tmpDistToZPole)
				{
					tmpDistToZPole = tmpMass;
					if (Leptons[i].pdgId > 0)
					{
						ZCandidateLeptonIndex = pair<int,int>(i,j);
					}
					else
					{
						ZCandidateLeptonIndex = pair<int,int>(j,i);
					}
					ZMass = tmpMass;
					ZPt = (Leptons[i].lepton+Leptons[j].lepton).Pt();
					ZCandidate = Leptons[i].lepton+Leptons[j].lepton;
					foundZ = true;
				}
			}
		}
		if (foundZ && fabs(ZMass-Z_MASS) < 30.0)
		{
			llp_tree->ZMass = ZMass;
			llp_tree->ZPt   = ZPt;
			llp_tree->ZEta  = ZCandidate.Eta();
			llp_tree->ZPhi  = ZCandidate.Phi();
			llp_tree->ZleptonIndex1 = ZCandidateLeptonIndex.first;
			llp_tree->ZleptonIndex2 = ZCandidateLeptonIndex.second;
		} // endif foundZ
		if( (label=="MR_ZLL") && isData && !(foundZ && fabs(ZMass-Z_MASS) < 30.0) ) continue;

		//if( presel && (label.find("MR") == std::string::npos) && !(llp_tree->nLeptons==0 && llp_tree->nMuons==0 && llp_tree->nElectrons==0 && llp_tree->nTaus==0 && llp_tree->nPhotons==0) ) continue; 
		//if( presel && (label.find("MR") == std::string::npos) && !(llp_tree->nLeptons==0 && llp_tree->nMuons==0 && llp_tree->nElectrons==0 && llp_tree->nTaus==0 && llp_tree->nPhotons==0) ) continue; 


		//-----------------------------------------------
		//Select Jets
		//-----------------------------------------------
		//std::vector<double> jetPtVector;
		//std::vector<double> jetCISVVector;
		std::vector<jets> Jets;
		//auto highest = [](auto a, auto b) { return a > b; };
		//cout <<"nJets :" << nJets << std::endl;

		if(_debug) std::cout << "nJets " << nJets << std::endl;
		if(_debug_pf) std::cout << "nJets " << nJets << std::endl;
		if(_debug_jet) std::cout << "nJets " << nJets << std::endl;
		if(_debug_jet) std::cout << "jetE 0 " << jetE[0] << std::endl;
		//if(_debug_jet) std::cout << "jetChargedEMEnergyFraction 0 " << jetChargedEMEnergyFraction[0] << std::endl;
		//if(_debug_jet) std::cout << "jetNeutralEMEnergyFraction 0 " << jetNeutralEMEnergyFraction[0] << std::endl;
		if(_debug_jet) std::cout << "jetGammaMax_ET 0 " << jetGammaMax_ET[0] << std::endl;

		//pt 30
		float jetMet_dPhiMin_eta_2p4_temp = 999;
		float jetMet_dPhiMin_eta_3_temp = 999;
		float jetMet_dPhiMin_eta_all_temp = 999;
		//pt 20
		float jetMet_dPhiMin_pt_20_eta_2p4_temp = 999;
		float jetMet_dPhiMin_pt_20_eta_3_temp = 999;
		float jetMet_dPhiMin_pt_20_eta_all_temp = 999;

		float ht = 0.;

		//all pt
		int nAK4Jets_in_HEM = 0;
		int nAK4Jets_in_HEM_eta_2p4 = 0;
		int nAK4Jets_in_HEM_eta_2p5 = 0;
		int nAK4Jets_in_HEM_eta_3 = 0;
		//pt 20
		int nAK4Jets_in_HEM_pt_20 = 0;
		int nAK4Jets_in_HEM_pt_20_eta_2p4 = 0;
		int nAK4Jets_in_HEM_pt_20_eta_2p5 = 0;
		int nAK4Jets_in_HEM_pt_20_eta_3 = 0;
		//pt 30
		int nAK4Jets_in_HEM_pt_30 = 0;
		int nAK4Jets_in_HEM_pt_30_eta_2p4 = 0;
		int nAK4Jets_in_HEM_pt_30_eta_2p5 = 0;
		int nAK4Jets_in_HEM_pt_30_eta_3 = 0;

		//pick max sum track pt inside jet wrt. pv
		float tmp_sumJetTrkPt[nJets][nPVAll];
		for(int iJet=0; iJet<nJets;iJet++)
		{
			for(int iPV=0; iPV<nPVAll;iPV++)
			{
				tmp_sumJetTrkPt[iJet][iPV] = 0;
			}
		}
		for(int i = 0; i < nJets; i++)
		{
			//HEM: reject events with jets in problematic region
			//Affected runs: 2018, after  >=319077 
			//if(Jets->at(j).eta>-3. and Jets->at(j).eta<-1.3 and Jets->at(j).phi>-1.57 and Jets->at(j).phi<-0.87)
			//add versions of pt>20/30, |eta|<2.4/3/5.2(all)
			if(jetEta[i]>-3. && jetEta[i]<-1.3 && jetPhi[i]>-1.57 && jetPhi[i]<-0.87)
			{
				nAK4Jets_in_HEM++;

			}
			if(jetEta[i]>-3. && jetEta[i]<-1.3 && jetPhi[i]>-1.57 && jetPhi[i]<-0.87 && fabs(jetEta[i])<2.4)
			{
				nAK4Jets_in_HEM_eta_2p4++;
			}
			if(jetEta[i]>-3. && jetEta[i]<-1.3 && jetPhi[i]>-1.57 && jetPhi[i]<-0.87 && fabs(jetEta[i])<2.5)
			{
				nAK4Jets_in_HEM_eta_2p5++;
			}
			if(jetEta[i]>-3. && jetEta[i]<-1.3 && jetPhi[i]>-1.57 && jetPhi[i]<-0.87 && fabs(jetEta[i])<3)
			{
				nAK4Jets_in_HEM_eta_3++;
			}
			//pt 20
			if(jetEta[i]>-3. && jetEta[i]<-1.3 && jetPhi[i]>-1.57 && jetPhi[i]<-0.87 && jetPt[i]>20.)
			{
				nAK4Jets_in_HEM_pt_20++;

			}
			if(jetEta[i]>-3. && jetEta[i]<-1.3 && jetPhi[i]>-1.57 && jetPhi[i]<-0.87 && fabs(jetEta[i])<2.4 && jetPt[i]>20.)
			{
				nAK4Jets_in_HEM_pt_20_eta_2p4++;
			}
			if(jetEta[i]>-3. && jetEta[i]<-1.3 && jetPhi[i]>-1.57 && jetPhi[i]<-0.87 && fabs(jetEta[i])<2.5 && jetPt[i]>20.)
			{
				nAK4Jets_in_HEM_pt_20_eta_2p5++;
			}
			if(jetEta[i]>-3. && jetEta[i]<-1.3 && jetPhi[i]>-1.57 && jetPhi[i]<-0.87 && fabs(jetEta[i])<3 && jetPt[i]>20.)
			{
				nAK4Jets_in_HEM_pt_20_eta_3++;
			}
			//pt 30
			if(jetEta[i]>-3. && jetEta[i]<-1.3 && jetPhi[i]>-1.57 && jetPhi[i]<-0.87 && jetPt[i]>30.)
			{
				nAK4Jets_in_HEM_pt_30++;

			}
			if(jetEta[i]>-3. && jetEta[i]<-1.3 && jetPhi[i]>-1.57 && jetPhi[i]<-0.87 && fabs(jetEta[i])<2.4 && jetPt[i]>30.)
			{
				nAK4Jets_in_HEM_pt_30_eta_2p4++;
			}
			if(jetEta[i]>-3. && jetEta[i]<-1.3 && jetPhi[i]>-1.57 && jetPhi[i]<-0.87 && fabs(jetEta[i])<2.5 && jetPt[i]>30.)
			{
				nAK4Jets_in_HEM_pt_30_eta_2p5++;
			}
			if(jetEta[i]>-3. && jetEta[i]<-1.3 && jetPhi[i]>-1.57 && jetPhi[i]<-0.87 && fabs(jetEta[i])<3 && jetPt[i]>30.)
			{
				nAK4Jets_in_HEM_pt_30_eta_3++;
			}


			ht += jetPt[i];

			// pt>30, |eta|<2.4 , Min Delta Phi (jet, met)
			if(jetMet_dPhiMin_eta_2p4_temp > abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi)) && fabs(jetEta[i])<2.4 && jetPt[i]>30.)
			{
				jetMet_dPhiMin_eta_2p4_temp = abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi));
			}
			if(jetMet_dPhiMin_eta_3_temp > abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi)) && fabs(jetEta[i])<3 && jetPt[i]>30.)
			{
				jetMet_dPhiMin_eta_3_temp = abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi));
			}
			if(jetMet_dPhiMin_eta_all_temp > abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi)) && jetPt[i]>30.)
			{
				jetMet_dPhiMin_eta_all_temp = abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi));
			}

			// pt>20, |eta|<3 , Min Delta Phi (jet, met)
			// pt>20, all eta , Min Delta Phi (jet, met)
			if(jetMet_dPhiMin_pt_20_eta_2p4_temp > abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi)) && fabs(jetEta[i])<2.4 && jetPt[i]>20.)
			{
				jetMet_dPhiMin_pt_20_eta_2p4_temp = abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi));
			}
			if(jetMet_dPhiMin_pt_20_eta_3_temp > abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi)) && fabs(jetEta[i])<3 && jetPt[i]>20.)
			{
				jetMet_dPhiMin_pt_20_eta_3_temp = abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi));
			}
			if(jetMet_dPhiMin_pt_20_eta_all_temp > abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi)) && jetPt[i]>20.)
			{
				jetMet_dPhiMin_pt_20_eta_all_temp = abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi));
			}
			//if(jetMet_dPhiMin_eta_all_temp > abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi)) )
			//{
			//	jetMet_dPhiMin_eta_all_temp = abs(RazorAnalyzerLLP::deltaPhi(jetPhi[i],metType1Phi));
			//}

			if(xclean){
				//------------------------------------------------------------
				//exclude selected muons and electrons from the jet collection
				//------------------------------------------------------------
				double deltaR = -1;
				for(auto& lep : Leptons){
					double thisDR = RazorAnalyzerLLP::deltaR(jetEta[i],jetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
					if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
				}
				if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

				//------------------------------------------------------------
				//exclude selected taus from the jet collection
				//------------------------------------------------------------
				double deltaR_tau = -1;
				for(auto& tau : Taus){
					double thisDR_tau = RazorAnalyzerLLP::deltaR(jetEta[i],jetPhi[i],tau.tau.Eta(),tau.tau.Phi());
					if(deltaR_tau < 0 || thisDR_tau < deltaR_tau) deltaR_tau = thisDR_tau;
				}
				if(deltaR_tau > 0 && deltaR_tau < 0.4) continue; //jet matches a selected tau

				//------------------------------------------------------------
				//exclude selected photons from the jet collection
				//------------------------------------------------------------
				double deltaR_pho = -1;
				for(auto& pho : Photons){
					double thisDR_pho = RazorAnalyzerLLP::deltaR(jetEta[i],jetPhi[i],pho.photon.Eta(),pho.photon.Phi());
					if(deltaR_pho < 0 || thisDR_pho < deltaR_pho) deltaR_pho = thisDR_pho;
				}	
				if(deltaR_pho > 0 && deltaR_pho < 0.4) continue; //jet matches a selected photon
			}

			//------------------------------------------------------------
			//Apply Jet Energy and Resolution Corrections
			//------------------------------------------------------------

			TLorentzVector thisJet = makeTLorentzVector( jetPt[i], jetEta[i], jetPhi[i], jetE[i] );
			//double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
			//		fixedGridRhoFastjetAll, jetJetArea[i] , JetCorrector);

			//double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
			//		fixedGridRhoFastjetAll, jetJetArea[i],
			//		runNum,
			//		JetCorrectorIOV,JetCorrector);
			////cout <<"JEC :" << JEC << std::endl;

			//TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );

			if( thisJet.Pt() < 30 ) continue;//According to the April 1st 2015 AN
			//if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
			if( fabs( thisJet.Eta() ) >= 1. ) continue;
			//if( fabs( thisJet.Eta() ) >= 1.48 ) continue;
			//if( fabs( thisJet.Eta() ) >= 2.4 ) continue;
			//if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
			// if ( !jetPassIDLoose[i] ) continue;
			// if (!(jetRechitE[i] > 0.0)) continue;
			// if(jetNRechits[i]<10) continue;
			// if(jetRechitT[i] < 0.0) continue;
			// if ((jetChargedHadronEnergyFraction[i]+jetChargedEMEnergyFraction[i]) > 0.4) continue;
			// if ((jetChargedHadronEnergyFraction[i]+jetNeutralHadronEnergyFraction[i])/(jetChargedEMEnergyFraction[i]+jetNeutralEMEnergyFraction[i]) < 0.2) continue;

			// std::cout <<jetRechitT[i] << "," << jetRechitE[i] <<  "," << jetNRechits[i] << std::endl;
			//if( jetRechitT[i]<=-1 ) continue;
			//if(_debug_dnn) std::cout << "Jet Time " << jetRechitT[i] << std::endl;
			//std::cout << "Jet Time " << jetRechitT[i] << std::endl;
			//std::cout << "Jet Time rms" << jetRechitT_rms[i] << std::endl;
			double jetTimeRecHitsECAL = -100;
			double jetTimeRmsECAL = 0;
			double jetTimeRmsECAL_fix = 0;
			if (isnan(jetRechitT[i]) || jetRechitE[i] == 0) {
				jetTimeRecHitsECAL = -100;
			} else {
				jetTimeRecHitsECAL = jetRechitT[i];
				jetTimeRmsECAL = jetRechitT_rms[i];
			}
			if(_debug_dnn) std::cout << "Jet Time " << jetTimeRecHitsECAL << std::endl;
			if( presel_jet && jetTimeRecHitsECAL<=-1 ) continue;
			if( presel_jet && jetMuonEnergyFraction[i]>=0.6 ) continue;
			if( presel_jet && jetElectronEnergyFraction[i]>=0.6 ) continue;
			if( presel_jet && jetPhotonEnergyFraction[i]>=0.8 ) continue;

			// ***********************************
			//Compute Rechit Quantities
			// ***********************************
			double jetEnergyRecHitsECAL = 0;
			double jetEnergyRecHitsHCAL = 0;
			double jetCscEnergyRecHitsECAL = 0;
			double jetDtEnergyRecHitsECAL = 0;
			double jetCscEnergyRecHitsEF = 0;
			double jetDtEnergyRecHitsEF = 0;
			//double jetTimeRecHitsECAL = -100;
			double jetTimeRecHitsHCAL = -100;
			double tmpJetTimeEnergyRecHitsHCAL = 0;
			int jetNRecHitsECAL = 0;
			int jetNRecHitsHCAL = 0;
			std::vector<double> ebrechitphi;
			std::vector<double> ebrechiteta;
			std::vector<double> ebrechitet;
			std::vector<double> ebrechitetsq;
			bool jetEcalRechit2Csc[nRechits];
			bool jetEcalRechit2Dt[nRechits];

			if(_debug_pf) std::cout << "this jet pt " << thisJet.Pt() << std::endl;

			//Loop over ECAL rechits
			for (int q=0; q < nRechits; q++) {

				jetEcalRechit2Csc[q] = false;
				if(_debug_csc) cout << "flag " << jetEcalRechit2Csc[q] << "\n";
				jetEcalRechit2Dt[q] = false;

				if (ecalRechit_E[q] <= 0.5) continue;
				if (abs(ecalRechit_Eta[q]) >= 1.48) continue;
				double tmpDR = RazorAnalyzerLLP::deltaR(thisJet.Eta(), thisJet.Phi(), ecalRechit_Eta[q], ecalRechit_Phi[q]);
				if (tmpDR > 0.4) continue;			  
				if (ecalRechit_kSaturatedflag[q] || 
						ecalRechit_kLeadingEdgeRecoveredflag[q] || 
						ecalRechit_kPoorRecoflag[q] ||
						ecalRechit_kWeirdflag[q] || 
						ecalRechit_kDiWeirdflag[q]) continue;
				if (ecalRechit_T_Error[q] < 0 || ecalRechit_T_Error[q] > 100) continue;
				if (abs(ecalRechit_T[q]) > 12.5) continue;
				if (abs(ecalRechit_Eta[q]) > 1.5) continue;

				//cout << "Rechit " << q << " : " << ecalRechit_E[q] << "\n";

				ebrechitphi.push_back(ecalRechit_Phi[q]);
				ebrechiteta.push_back(ecalRechit_Eta[q]);
				ebrechitet.push_back(ecalRechit_E[q]/cosh(ecalRechit_Eta[q]));
				ebrechitetsq.push_back( pow(ecalRechit_E[q]/cosh(ecalRechit_Eta[q]),2) );

				jetEnergyRecHitsECAL += ecalRechit_E[q];
				jetNRecHitsECAL++;

				//CSC EF
				for (int kc=0; kc<nCscSeg; kc++) {
					double tmpcscDPhi = abs(RazorAnalyzerLLP::deltaPhi(ecalRechit_Phi[q], cscSegPhi[kc]));
					if(tmpcscDPhi>0.04) continue;
					jetEcalRechit2Csc[q] = true;
					if(_debug_csc) cout << "Rechit " << q << "flag " << jetEcalRechit2Csc[q] << "\n";
				}//CSC EF

				//DT EF
				for (int kd=0; kd<nDtSeg; kd++) {
					double tmpdtDPhi = abs(RazorAnalyzerLLP::deltaPhi(ecalRechit_Phi[q], cscSegPhi[kd]));
					if(tmpdtDPhi>0.04) continue;
					jetEcalRechit2Dt[q] = true;
				}//DT EF

			}  

			//ecal rechit match to csc
			for (int q=0; q < nRechits; q++) {
				if(jetEcalRechit2Csc[q]) jetCscEnergyRecHitsECAL += ecalRechit_E[q]; 
				if(jetEcalRechit2Dt[q]) jetDtEnergyRecHitsECAL += ecalRechit_E[q]; 
			}  
			if(_debug_csc) cout << "CscSeg E " << jetCscEnergyRecHitsECAL  << "\n";

			double jetsig1EB(-1.),jetsig2EB(-1.);
			RazorAnalyzerLLP::jet_second_moments(ebrechitet,ebrechiteta,ebrechitphi,jetsig1EB,jetsig2EB);
			double jetptDEB = -1;
			if ( accumulate(ebrechitet.begin(),ebrechitet.end(),0) > 0) {
				jetptDEB = sqrt( accumulate(ebrechitetsq.begin(),ebrechitetsq.end(),0)) / accumulate(ebrechitet.begin(),ebrechitet.end(),0) ;
			}	
			if(_debug_pf) std::cout << "this jet ptDEB " << jetptDEB << std::endl;
			if (jetNRecHitsECAL == 0) {
				jetEnergyRecHitsECAL = -1;
				jetNRecHitsECAL = -1;
				jetCscEnergyRecHitsECAL = -1;
				jetDtEnergyRecHitsECAL = -1;
				jetCscEnergyRecHitsEF = -1;
				jetDtEnergyRecHitsEF = -1;
			}
			else {
				jetCscEnergyRecHitsEF = jetCscEnergyRecHitsECAL/jetEnergyRecHitsECAL;
				jetDtEnergyRecHitsEF = jetDtEnergyRecHitsECAL/jetEnergyRecHitsECAL;
			}
			jetTimeRmsECAL_fix = (jetNRecHitsECAL>0) ? jetRechitT_rms[i]/sqrt(jetNRecHitsECAL) : -1.;
			//if (isnan(jetRechitT[i]) || jetRechitE[i] == 0) {
			//  jetTimeRecHitsECAL = -100;
			//} else {
			//  jetTimeRecHitsECAL = jetRechitT[i];
			//}
			///
			//			if(eventNum==39746)
			//			{
			//				cout << "ECAL energy, N, time " << jetEnergyRecHitsECAL << " , " << jetNRecHitsECAL << " : " << jetTimeRecHitsECAL << "\n";			
			//			}
			///
			//Loop over HCAL rechits
			for (int q=0; q < nHBHERechits; q++) {
				if (hbheRechit_E[q] <= 0.1) continue;
				double tmpDR = RazorAnalyzerLLP::deltaR(thisJet.Eta(), thisJet.Phi(), hbheRechit_Eta[q], hbheRechit_Phi[q]);
				if (tmpDR > 0.4) continue;

				jetEnergyRecHitsHCAL += hbheRechit_E[q];
				jetNRecHitsHCAL++;
				tmpJetTimeEnergyRecHitsHCAL += hbheRechit_E[q]*hbheRechit_T[q];
			}

			if (jetEnergyRecHitsHCAL > 0) {
				jetTimeRecHitsHCAL = tmpJetTimeEnergyRecHitsHCAL / jetEnergyRecHitsHCAL;
			} else {
				jetEnergyRecHitsHCAL = -1;
				jetTimeRecHitsHCAL = -100;
			}
			if(_debug_pf) std::cout << "this jet time hcal " << jetTimeRecHitsHCAL << std::endl;

			std::vector<double> pfcandphi;
			std::vector<double> pfcandeta;
			std::vector<double> pfcandpt;
			std::vector<double> pfcandptsq;

			//Loop over PF candidates
			for (int q=0; q < jetNPFCands[i]; q++) {
				int thisIndex = jetPFCandIndex[i][q];
				if (abs(PFCandidateEta[thisIndex]) >= 1.48) continue;
				double tmpDR = RazorAnalyzerLLP::deltaR(thisJet.Eta(), thisJet.Phi(), PFCandidateEta[thisIndex], PFCandidatePhi[thisIndex]);
				if (tmpDR > 0.4) continue;			  

				pfcandphi.push_back(PFCandidatePhi[thisIndex]);
				pfcandeta.push_back(PFCandidateEta[thisIndex]);
				pfcandpt.push_back(PFCandidatePt[thisIndex]);
				pfcandptsq.push_back( pow(PFCandidatePt[thisIndex],2) );

			}

			double jetsig1PF(-1.),jetsig2PF(-1.);
			RazorAnalyzerLLP::jet_second_moments(pfcandpt,pfcandeta,pfcandphi,jetsig1PF,jetsig2PF);
			double jetptDPF = -1;
			if ( accumulate(pfcandpt.begin(),pfcandpt.end(),0) > 0) {
				jetptDPF = sqrt( accumulate(pfcandptsq.begin(),pfcandptsq.end(),0)) / accumulate(pfcandpt.begin(),pfcandpt.end(),0) ;
			}	
			if(_debug_pf) std::cout << "this jet sig1PF " << jetsig1PF << std::endl;
			if(_debug_pf) std::cout << "this jet sig2PF " << jetsig2PF << std::endl;
			if(_debug_pf) std::cout << "this jet ptDPF " << jetptDPF << std::endl;


			// cout << thisJet.Pt() << " " << thisJet.Eta() << " " << thisJet.Phi() << " | "
			//      << jetEnergyRecHitsECAL << " " << jetEnergyRecHitsHCAL << " : " << jetNRecHitsECAL << " " 
			//      << jetNRecHitsHCAL << " | " 
			//      << jetTimeRecHitsECAL << " " << jetTimeRecHitsHCAL << "\n";

			// ***********************************
			//Evaluate NN tagger
			// ***********************************
			//std::vector<std::string> inputFeaturesV1 = { "Jet_0_nTrackConstituents","Jet_0_nSelectedTracks", "Jet_0_timeRecHitsEB", "Jet_0_energyRecHitsEB", "Jet_0_nRecHitsEB", "Jet_0_cHadEFrac", "Jet_0_nHadEFrac", "Jet_0_eleEFrac", "Jet_0_photonEFrac", "Jet_0_ptAllTracks", "Jet_0_ptAllPVTracks", "Jet_0_alphaMax", "Jet_0_betaMax", "Jet_0_gammaMax", "Jet_0_gammaMaxEM", "Jet_0_gammaMaxHadronic", "Jet_0_gammaMaxET","Jet_0_minDeltaRAllTracks","Jet_0_minDeltaRPVTracks",};
			///
			//			inputValuesV1[0] = jetChargedHadronMultiplicity[i]+jetElectronMultiplicity[i]+jetMuonMultiplicity[i];
			//			inputValuesV1[1] = jetNSelectedTracks[i];
			//			  //std::cout<< " input value 1: " << jetNSelectedTracks[i] <<std::endl;
			//			inputValuesV1[2] = jetTimeRecHitsECAL;
			//			inputValuesV1[3] = (jetEnergyRecHitsECAL == 0) ? -1 : sqrt(jetEnergyRecHitsECAL);
			//			inputValuesV1[4] = jetNRecHitsECAL;
			//			inputValuesV1[5] = jetChargedHadronEnergyFraction[i];
			//			inputValuesV1[6] = jetNeutralHadronEnergyFraction[i];
			//			inputValuesV1[7] = jetElectronEnergyFraction[i];
			//			inputValuesV1[8] = jetPhotonEnergyFraction[i];
			//			inputValuesV1[9] = (jetPtAllTracks[i] == -99) ? -1 : jetPtAllTracks[i];
			//			inputValuesV1[10] = (jetPtAllPVTracks[i] == -99 || jetPtAllPVTracks[i] == 0) ? -1 : jetPtAllPVTracks[i];
			//			inputValuesV1[11] = (jetAlphaMax[i] == -99) ? -100 : jetAlphaMax[i];
			//			inputValuesV1[12] = (jetBetaMax[i] == -99) ? -100 : jetBetaMax[i];
			//			inputValuesV1[13] = (jetGammaMax[i] == -99) ? -100 : jetGammaMax[i];
			//			inputValuesV1[14] = (jetGammaMax_EM[i] == -99) ? -100 : jetGammaMax_EM[i];
			//			inputValuesV1[15] = (jetGammaMax_Hadronic[i] == -99) ? -100 : jetGammaMax_Hadronic[i];
			//			inputValuesV1[16] = (jetGammaMax_ET[i] == -99) ? -100 : jetGammaMax_ET[i];
			//			inputValuesV1[17] = (jetMinDeltaRAllTracks[i] == -99 || jetMinDeltaRAllTracks[i] == 15) ? 999 : jetMinDeltaRAllTracks[i];
			//			inputValuesV1[18] = (jetMinDeltaRPVTracks[i] == -99 || jetMinDeltaRPVTracks[i] == 15) ? 999 : jetMinDeltaRPVTracks[i];
			//
			//			// fill the input tensor using a data pointer that is shifted consecutively
			//			float* dV1 = inputTensorV1.flat<float>().data();
			//			for (float vV1 : inputValuesV1) {
			//			  //std::cout<< " input value: " << v <<std::endl;
			//			  *dV1 = vV1;
			//			  dV1++;
			//			}
			//
			//			// run the inference
			//			std::vector<tensorflow::Tensor> outputsV1;		
			//			tensorflow::run(sessionV1, {{inputTensorNameV1, inputTensorV1}}, {outputTensorNameV1}, &outputsV1, threadPoolV1);
			//			
			//			// the result
			//			double outputValueV1 = outputsV1[0].matrix<float>()(0, 1);
			//			//std::cout << "output value: " << outputValue << std::endl;
			//			//std::cout << "\n" << std::endl;
			//			
			//	//std::vector<std::string> inputFeatures = { "Jet_nTrackConstituents", "Jet_nSelectedTracks", "Jet_timeRecHitsEB", "Jet_eFracRecHitsEB", "Jet_nRecHitsEB", "Jet_sig1EB", "Jet_sig2EB", "Jet_ptDEB", "Jet_sig1PF", "Jet_sig2PF", "Jet_ptDPF", "Jet_cHadEFrac", "Jet_nHadEFrac", "Jet_eleEFrac", "Jet_photonEFrac", "Jet_ptAllTracks", "Jet_ptAllPVTracks", "Jet_alphaMax", "Jet_betaMax", "Jet_gammaMax", "Jet_gammaMaxEM", "Jet_gammaMaxHadronic", "Jet_gammaMaxET", "Jet_minDeltaRAllTracks", "Jet_minDeltaRPVTracks",};
			//			inputValues[0] = jetChargedHadronMultiplicity[i]+jetElectronMultiplicity[i]+jetMuonMultiplicity[i];
			//			inputValues[1] = jetNSelectedTracks[i];
			//			  //std::cout<< " input value 1: " << jetNSelectedTracks[i] <<std::endl;
			//			inputValues[2] = jetTimeRecHitsECAL;
			//			inputValues[3] = (jetEnergyRecHitsECAL == 0) ? -1 : (jetEnergyRecHitsECAL/jetE[i]);
			//			inputValues[4] = jetNRecHitsECAL;
			//			inputValues[5] = jetsig1EB;
			//			inputValues[6] = jetsig2EB;
			//			inputValues[7] = jetptDEB;
			//			inputValues[8] = jetsig1PF;
			//			inputValues[9] = jetsig2PF;
			//			inputValues[10] = jetptDPF;
			//			inputValues[11] = jetChargedHadronEnergyFraction[i];
			//			inputValues[12] = jetNeutralHadronEnergyFraction[i];
			//			inputValues[13] = jetElectronEnergyFraction[i];
			//			inputValues[14] = jetPhotonEnergyFraction[i];
			//			inputValues[15] = (jetPtAllTracks[i] == -99) ? -1 : jetPtAllTracks[i];
			//			inputValues[16] = (jetPtAllPVTracks[i] == -99 || jetPtAllPVTracks[i] == 0) ? -1 : jetPtAllPVTracks[i];
			//			inputValues[17] = (jetAlphaMax[i] == -99) ? -100 : jetAlphaMax[i];
			//			inputValues[18] = (jetBetaMax[i] == -99) ? -100 : jetBetaMax[i];
			//			inputValues[19] = (jetGammaMax[i] == -99) ? -100 : jetGammaMax[i];
			//			inputValues[20] = (jetGammaMax_EM[i] == -99) ? -100 : jetGammaMax_EM[i];
			//			inputValues[21] = (jetGammaMax_Hadronic[i] == -99) ? -100 : jetGammaMax_Hadronic[i];
			//			inputValues[22] = (jetGammaMax_ET[i] == -99) ? -100 : jetGammaMax_ET[i];
			//			inputValues[23] = (jetMinDeltaRAllTracks[i] == -99 || jetMinDeltaRAllTracks[i] == 15) ? 999 : jetMinDeltaRAllTracks[i];
			//			inputValues[24] = (jetMinDeltaRPVTracks[i] == -99 || jetMinDeltaRPVTracks[i] == 15) ? 999 : jetMinDeltaRPVTracks[i];
			//
			//			// fill the input tensor using a data pointer that is shifted consecutively
			//			float* d = inputTensor.flat<float>().data();
			//			for (float v : inputValues) {
			//			  //std::cout<< " input value: " << v <<std::endl;
			//			  *d = v;
			//			  d++;
			//			}
			//
			//			// run the inference
			//			std::vector<tensorflow::Tensor> outputs;		
			//			tensorflow::run(session, {{inputTensorName, inputTensor}}, {outputTensorName}, &outputs, threadPool);
			//			
			//			// the result
			//			double outputValue = outputs[0].matrix<float>()(0, 1);
			//			//std::cout << "output value: " << outputValue << std::endl;
			//			//std::cout << "\n" << std::endl;
			///			

			//std::vector<std::string> inputFeaturesV3 = { "Jet_nTrackConstituents", "Jet_nSelectedTracks", "Jet_timeRecHitsEB", "Jet_eFracRecHitsEB", "Jet_nRecHitsEB", "Jet_sig1EB", "Jet_sig2EB", "Jet_ptDEB", "Jet_cHadEFrac", "Jet_nHadEFrac", "Jet_eleEFrac", "Jet_photonEFrac", "Jet_ptAllTracks", "Jet_ptAllPVTracks", "Jet_alphaMax", "Jet_betaMax", "Jet_gammaMax", "Jet_gammaMaxEM", "Jet_gammaMaxHadronic", "Jet_gammaMaxET", "Jet_minDeltaRAllTracks", "Jet_minDeltaRPVTracks",};
			inputValuesV3[0] = jetChargedHadronMultiplicity[i]+jetElectronMultiplicity[i]+jetMuonMultiplicity[i];
			inputValuesV3[1] = jetNSelectedTracks[i];
			//std::cout<< " input value 1: " << jetNSelectedTracks[i] <<std::endl;
			inputValuesV3[2] = jetTimeRecHitsECAL;
			inputValuesV3[3] = (jetEnergyRecHitsECAL == 0) ? -1 : (jetEnergyRecHitsECAL/jetE[i]);
			inputValuesV3[4] = jetNRecHitsECAL;
			inputValuesV3[5] = (jetsig1EB<=0) ? -1: jetsig1EB;
			inputValuesV3[6] = (jetsig2EB<=0) ? -1: jetsig2EB;
			inputValuesV3[7] = jetptDEB;
			inputValuesV3[8] = jetChargedHadronEnergyFraction[i];
			inputValuesV3[9] = jetNeutralHadronEnergyFraction[i];
			inputValuesV3[10] = jetElectronEnergyFraction[i];
			inputValuesV3[11] = jetPhotonEnergyFraction[i];
			inputValuesV3[12] = (jetPtAllTracks[i] == -99 || jetPtAllTracks[i] == 0) ? -1 : jetPtAllTracks[i];
			inputValuesV3[13] = (jetPtAllPVTracks[i] == -99 || jetPtAllPVTracks[i] == 0) ? -1 : jetPtAllPVTracks[i];
			inputValuesV3[14] = (jetAlphaMax[i] == -99 || jetPtAllTracks[i] == 0) ? -100 : jetAlphaMax[i];
			inputValuesV3[15] = (jetBetaMax[i] == -99 || jetPtAllTracks[i] == 0) ? -100 : jetBetaMax[i];
			inputValuesV3[16] = (jetGammaMax[i] == -99 || jetPtAllTracks[i] == 0) ? -100 : jetGammaMax[i];
			//inputValuesV3[15] = (jetBetaMax[i] == -99 || jetPt[i] == 0) ? -100 : jetBetaMax[i];
			//inputValuesV3[16] = (jetGammaMax[i] == -99 || jetE[i] == 0) ? -100 : jetGammaMax[i];
			inputValuesV3[17] = (jetGammaMax_EM[i] == -99 || jetPtAllTracks[i] == 0 ||jetE[i]*(jetPhotonEnergyFraction[i]+jetElectronEnergyFraction[i]) == 0 || isnan(jetGammaMax_EM[i])) ? -100 : jetGammaMax_EM[i];
			inputValuesV3[18] = (jetGammaMax_Hadronic[i] == -99 || jetPtAllTracks[i] == 0 || jetE[i]*(jetNeutralHadronEnergyFraction[i]+jetChargedHadronEnergyFraction[i]) == 0) ? -100 : jetGammaMax_Hadronic[i];
			inputValuesV3[19] = (jetGammaMax_ET[i] == -99 || jetPtAllTracks[i] == 0 || thisJet.Et() == 0) ? -100 : jetGammaMax_ET[i];
			///
			//		inputValuesV3[14] = (jetAlphaMax[i] == -99 || jetAlphaMax[i] == 0) ? -100 : jetAlphaMax[i];
			//		inputValuesV3[15] = (jetBetaMax[i] == -99 || jetBetaMax[i] == 0) ? -100 : jetBetaMax[i];
			//		inputValuesV3[16] = (jetGammaMax[i] == -99 || jetGammaMax[i] == 0) ? -100 : jetGammaMax[i];
			//		inputValuesV3[17] = (jetGammaMax_EM[i] == -99 || jetGammaMax_EM[i] == 0 || isnan(jetGammaMax_EM[i])) ? -100 : jetGammaMax_EM[i];
			//		inputValuesV3[18] = (jetGammaMax_Hadronic[i] == -99 || jetGammaMax_Hadronic[i] == 0) ? -100 : jetGammaMax_Hadronic[i];
			//		inputValuesV3[19] = (jetGammaMax_ET[i] == -99 || jetGammaMax_ET[i] == 0) ? -100 : jetGammaMax_ET[i];
			///	
			inputValuesV3[20] = (jetMinDeltaRAllTracks[i] == -99 || jetMinDeltaRAllTracks[i] == 15) ? 999 : jetMinDeltaRAllTracks[i];
			inputValuesV3[21] = (jetMinDeltaRPVTracks[i] == -99 || jetMinDeltaRPVTracks[i] == 15) ? 999 : jetMinDeltaRPVTracks[i];
			///
			//			if(eventNum==5337||eventNum==39906 || eventNum==24897)
			//			{
			//				
			//			  std::cout<< " evt "<<eventNum <<std::endl;
			//			 for(int rr=0;rr<22;rr++)
			//			{
			//			  std::cout<< " input value "<<rr<<" : " << inputValuesV3[rr] <<std::endl;
			//			}
			//			std::cout<< " input value 0 : " << jetChargedHadronMultiplicity[i]+jetElectronMultiplicity[i]+jetMuonMultiplicity[i]<<std::endl;
			//			std::cout<< " input value 1 : " << jetNSelectedTracks[i]<<std::endl;
			//			std::cout<< " input value 2 : " << jetTimeRecHitsECAL<<std::endl;
			//			std::cout<< " input value 3 : " <<  (jetEnergyRecHitsECAL/jetE[i])<<std::endl;
			//			std::cout<< " input value 4 : " << jetNRecHitsECAL<<std::endl;
			//			std::cout<< " input value 5 : " << jetsig1EB<<std::endl;
			//			std::cout<< " input value 6 : " << jetsig2EB<<std::endl;
			//			std::cout<< " input value 7 : " << jetptDEB<<std::endl;
			//			std::cout<< " input value 8 : " << jetChargedHadronEnergyFraction[i]<<std::endl;
			//			std::cout<< " input value 9 : " << jetNeutralHadronEnergyFraction[i]<<std::endl;
			//			std::cout<< " input value 10 : " << jetElectronEnergyFraction[i]<<std::endl;
			//			std::cout<< " input value 11 : " << jetPhotonEnergyFraction[i]<<std::endl;
			//			std::cout<< " input value 12 : " <<  jetPtAllTracks[i]<<std::endl;
			//			std::cout<< " input value 13 : " <<  jetPtAllPVTracks[i]<<std::endl;
			//			std::cout<< " input value 14 : " <<  jetAlphaMax[i]<<std::endl;
			//			std::cout<< " input value 15 : " <<  jetBetaMax[i]<<std::endl;
			//			std::cout<< " input value 16 : " <<  jetGammaMax[i]<<std::endl;
			//			std::cout<< " input value 17 : " <<  jetGammaMax_EM[i]<<std::endl;
			//			std::cout<< " input value 18 : " <<  jetGammaMax_Hadronic[i]<<std::endl;
			//			std::cout<< " input value 19 : " <<  jetGammaMax_ET[i]<<std::endl;
			//			std::cout<< " input value 20 : " <<  jetMinDeltaRAllTracks[i]<<std::endl;
			//			std::cout<< " input value 21 : " << jetMinDeltaRPVTracks[i] <<std::endl;
			//
			//			}

			// fill the input tensor using a data pointer that is shifted consecutively
			float* dV3 = inputTensorV3.flat<float>().data();
			for (float vV3 : inputValuesV3) {
				//std::cout<< " input value: " << v <<std::endl;
				*dV3 = vV3;
				dV3++;
			}

			// run the inference
			std::vector<tensorflow::Tensor> outputsV3;		
			tensorflow::run(sessionV3, {{inputTensorNameV3, inputTensorV3}}, {outputTensorNameV3}, &outputsV3, threadPoolV3);

			// the result
			double outputValueV3 = outputsV3[0].matrix<float>()(0, 1);
			//std::cout << "output value: " << outputValue << std::endl;
			//std::cout << "\n" << std::endl;

			//---------v3 miniAOD--------
			inputValuesV3miniAOD[0] = jetChargedHadronMultiplicity[i]+jetElectronMultiplicity[i]+jetMuonMultiplicity[i];
			inputValuesV3miniAOD[1] = jetNSelectedTracks[i];
			//std::cout<< " input value 1: " << jetNSelectedTracks[i] <<std::endl;
			inputValuesV3miniAOD[2] = jetTimeRecHitsECAL;
			inputValuesV3miniAOD[3] = (jetEnergyRecHitsECAL == 0) ? -1 : (jetEnergyRecHitsECAL/jetE[i]);
			inputValuesV3miniAOD[4] = jetNRecHitsECAL;
			inputValuesV3miniAOD[5] = jetsig1EB;
			inputValuesV3miniAOD[6] = jetsig2EB;
			inputValuesV3miniAOD[7] = jetptDEB;
			inputValuesV3miniAOD[8] = jetChargedHadronEnergyFraction[i];
			inputValuesV3miniAOD[9] = jetNeutralHadronEnergyFraction[i];
			inputValuesV3miniAOD[10] = jetElectronEnergyFraction[i];
			inputValuesV3miniAOD[11] = jetPhotonEnergyFraction[i];

			// fill the input tensor using a data pointer that is shifted consecutively
			float* dV3miniAOD = inputTensorV3miniAOD.flat<float>().data();
			for (float vV3miniAOD : inputValuesV3miniAOD) {
				//std::cout<< " input value: " << v <<std::endl;
				*dV3miniAOD = vV3miniAOD;
				dV3miniAOD++;
			}

			// run the inference
			std::vector<tensorflow::Tensor> outputsV3miniAOD;		
			tensorflow::run(sessionV3miniAOD, {{inputTensorNameV3miniAOD, inputTensorV3miniAOD}}, {outputTensorNameV3miniAOD}, &outputsV3miniAOD, threadPoolV3miniAOD);

			// the result
			double outputValueV3miniAOD = outputsV3miniAOD[0].matrix<float>()(0, 1);
			//std::cout << "output value: " << outputValue << std::endl;
			//std::cout << "\n" << std::endl;

			jets tmpJet;
			tmpJet.jet    = thisJet;
			tmpJet.time   = jetRechitT[i];
			tmpJet.ecalNRechits = jetNRechits[i];
			tmpJet.ecalRechitE = jetRechitE[i];
			tmpJet.jetNeutralHadronMultiplicity = jetNeutralHadronMultiplicity[i];
			tmpJet.jetChargedHadronMultiplicity = jetChargedHadronMultiplicity[i];
			tmpJet.jetMuonMultiplicity = jetMuonMultiplicity[i];
			tmpJet.jetElectronMultiplicity = jetElectronMultiplicity[i];
			tmpJet.jetPhotonMultiplicity = jetPhotonMultiplicity[i];
			tmpJet.jetNeutralHadronEnergyFraction = jetNeutralHadronEnergyFraction[i];
			tmpJet.jetChargedHadronEnergyFraction = jetChargedHadronEnergyFraction[i];
			tmpJet.jetMuonEnergyFraction = jetMuonEnergyFraction[i];
			tmpJet.jetElectronEnergyFraction = jetElectronEnergyFraction[i];
			tmpJet.jetPhotonEnergyFraction = jetPhotonEnergyFraction[i];
			tmpJet.jetCSV = jetCISV[i];


			tmpJet.jetAlphaMax = jetAlphaMax[i];
			tmpJet.jetBetaMax = jetBetaMax[i];
			tmpJet.jetGammaMax = jetGammaMax[i];
			tmpJet.jetGammaMax_Hadronic = jetGammaMax_Hadronic[i];
			tmpJet.jetGammaMax_EM = jetGammaMax_EM[i];
			tmpJet.jetGammaMax_ET = jetGammaMax_ET[i];
			tmpJet.jetPtAllPVTracks = jetPtAllPVTracks[i];
			tmpJet.jetPtAllTracks = jetPtAllTracks[i];
			tmpJet.jetMinDeltaRPVTracks = jetMinDeltaRPVTracks[i];
			tmpJet.jetMinDeltaRAllTracks = jetMinDeltaRAllTracks[i];

			tmpJet.jetChargedEMEnergyFraction = jetElectronEnergyFraction[i];
			tmpJet.jetNeutralEMEnergyFraction = jetPhotonEnergyFraction[i];

			tmpJet.jetNVertexTracks = jetNVertexTracks[i];
			tmpJet.jetNSelectedTracks = jetNSelectedTracks[i];
			tmpJet.jetDRSVJet = jetDRSVJet[i];
			tmpJet.jetSVMass = jetSVMass[i];

			tmpJet.jetChargedMultiplicity = jetChargedHadronMultiplicity[i]+jetElectronMultiplicity[i]+jetMuonMultiplicity[i];
			//ecal rechits
			tmpJet.jetNRecHitsEcal = jetNRecHitsECAL;
			tmpJet.jetEnergyRecHitsEcal = jetEnergyRecHitsECAL;
			tmpJet.jetTimeRecHitsEcal = jetTimeRecHitsECAL;
			tmpJet.jetTimeRmsEcal = jetTimeRmsECAL;
			tmpJet.jetTimeRmsEcal_fix = jetTimeRmsECAL_fix;

			//hcal hbhe rechits
			tmpJet.jetNRecHitsHcal = jetNRecHitsHCAL;
			tmpJet.jetEnergyRecHitsHcal = jetEnergyRecHitsHCAL;
			tmpJet.jetTimeRecHitsHcal = jetTimeRecHitsHCAL;

			tmpJet.jetsig1EB = jetsig1EB;
			tmpJet.jetsig2EB = jetsig2EB;
			tmpJet.jetptDEB = jetptDEB;
			tmpJet.jetsig1PF = jetsig1PF;
			tmpJet.jetsig2PF = jetsig2PF;
			tmpJet.jetptDPF = jetptDPF;

			//tmpJet.dnn_score_v1 = outputValueV1;
			//tmpJet.dnn_score = outputValue;
			tmpJet.dnn_score_v3 = outputValueV3;
			tmpJet.dnn_score_v3_miniAOD = outputValueV3miniAOD;

			tmpJet.jetCscEF = jetCscEnergyRecHitsEF;
			tmpJet.jetDtEF = jetDtEnergyRecHitsEF;

			if(_debug_trk) std::cout << "nTracks" << nTracks << std::endl;
			std::vector<float> nPixelHits;
			std::vector<float> nHits;
			//int nChargedMultiplicity=0;
			int pvid=-1;
			for(int iTrack=0; iTrack<nTracks; iTrack++)
			{
				double thisDR_trk = RazorAnalyzerLLP::deltaR(jetEta[i],jetPhi[i],track_Eta[iTrack], track_Phi[iTrack]);
				if(thisDR_trk<0.4)
				{
					//nChargedMultiplicity++;
					nPixelHits.push_back(track_nPixelHits[iTrack]);
					nHits.push_back(track_nHits[iTrack]);
					//track check other PV
					pvid = track_bestVertexIndex[iTrack];
					tmp_sumJetTrkPt[i][pvid] += track_Pt[iTrack];
					

				}

			}
			//pick max sum track pt inside jet wrt. pv
			int maxpvid = -1;
			int maxsvid = -1;
			float maxpvpt = 0;
			float maxsvpt = 0;
			for(int ipv=0; ipv<nPVAll; ipv++)
			{
				if(tmp_sumJetTrkPt[i][ipv] > maxpvpt)
				{
					maxpvid = ipv;
					maxpvpt = tmp_sumJetTrkPt[i][ipv];
				}

				if(ipv!=pvid0)
				{
					if(tmp_sumJetTrkPt[i][ipv] > maxsvpt)
					{
						maxsvid = ipv;
						maxsvpt = tmp_sumJetTrkPt[i][ipv];		
					}
				}
			}
			//std::cout << "maxpvid" <<maxpvid <<endl;
			//std::cout << "maxpvpt" <<maxpvpt <<endl;
			tmpJet.maxpvid = maxpvid;
			tmpJet.maxsvid = maxsvid;
			tmpJet.maxpvpt = maxpvpt;
			tmpJet.maxsvpt = maxsvpt;

			std::sort(nPixelHits.begin(), nPixelHits.end());
			float nPixelHitsMedian = 0;
			if (nPixelHits.size() > 0) {
				if (nPixelHits.size() % 2 ==0) nPixelHitsMedian = ((nPixelHits[nPixelHits.size()/2 -1] + nPixelHits[nPixelHits.size()/2]) /2);
				else nPixelHitsMedian = nPixelHits[nPixelHits.size()/2];
			}    
			std::sort(nHits.begin(), nHits.end());
			float nHitsMedian = 0;
			if (nHits.size() > 0) {
				if (nHits.size() % 2 ==0) nHitsMedian = ((nHits[nHits.size()/2 -1] + nHits[nHits.size()/2]) /2);
				else nHitsMedian = nHits[nHits.size()/2];
			}

			//tmpJet.jetChargedMultiplicity = nChargedMultiplicity;
			tmpJet.jetNPixelHitsMedian = nPixelHitsMedian;
			tmpJet.jetNHitsMedian = nHitsMedian;


			Jets.push_back(tmpJet);

		}
		llp_tree->jetMet_dPhiMin_eta_2p4 = jetMet_dPhiMin_eta_2p4_temp;
		llp_tree->jetMet_dPhiMin_eta_3 = jetMet_dPhiMin_eta_3_temp;
		llp_tree->jetMet_dPhiMin_eta_all = jetMet_dPhiMin_eta_all_temp;
		if( presel && (label.find("MR") == std::string::npos) && llp_tree->jetMet_dPhiMin_eta_all < 0.5 ) continue;
		llp_tree->jetMet_dPhiMin_pt_20_eta_2p4 = jetMet_dPhiMin_pt_20_eta_2p4_temp;
		llp_tree->jetMet_dPhiMin_pt_20_eta_3 = jetMet_dPhiMin_pt_20_eta_3_temp;
		llp_tree->jetMet_dPhiMin_pt_20_eta_all = jetMet_dPhiMin_pt_20_eta_all_temp;
		llp_tree->HT = ht;
		llp_tree->nCHSJets_in_HEM = nAK4Jets_in_HEM;
		llp_tree->nCHSJets_in_HEM_eta_2p4 = nAK4Jets_in_HEM_eta_2p4;
		llp_tree->nCHSJets_in_HEM_eta_2p5 = nAK4Jets_in_HEM_eta_2p5;
		llp_tree->nCHSJets_in_HEM_eta_3 = nAK4Jets_in_HEM_eta_3;
		llp_tree->nCHSJets_in_HEM_pt_20 = nAK4Jets_in_HEM_pt_20;
		llp_tree->nCHSJets_in_HEM_pt_20_eta_2p4 = nAK4Jets_in_HEM_pt_20_eta_2p4;
		llp_tree->nCHSJets_in_HEM_pt_20_eta_2p5 = nAK4Jets_in_HEM_pt_20_eta_2p5;
		llp_tree->nCHSJets_in_HEM_pt_20_eta_3 = nAK4Jets_in_HEM_pt_20_eta_3;
		llp_tree->nCHSJets_in_HEM_pt_30 = nAK4Jets_in_HEM_pt_30;
		llp_tree->nCHSJets_in_HEM_pt_30_eta_2p4 = nAK4Jets_in_HEM_pt_30_eta_2p4;
		llp_tree->nCHSJets_in_HEM_pt_30_eta_2p5 = nAK4Jets_in_HEM_pt_30_eta_2p5;
		llp_tree->nCHSJets_in_HEM_pt_30_eta_3 = nAK4Jets_in_HEM_pt_30_eta_3;




		//sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

		//if (Jets.size()>0)
		//{
		//	llp_tree->jetMet_dPhi = RazorAnalyzerLLP::deltaPhi(jetPhi[0],metType1Phi);
		//	//TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM( metType1Pt, 0, metType1Phi, 0 );
		//	TLorentzVector jet0 = makeTLorentzVectorPtEtaPhiM( jetPt[0], 0, jetPhi[0], 0 );
		//	llp_tree->jetMet_dPhiStar = RazorAnalyzerLLP::deltaPhi(jetPhi[0],  (t1PFMET+jet0).Phi() );
		//}
		//else{
		//	llp_tree->jetMet_dPhi = -999.;
		//	llp_tree->jetMet_dPhiStar = -999.;
		//}

		////float jetMet_dPhiMin_temp = 999 ; 
		////float jetMet_dPhiStarMin_temp = 999 ; 
		////float jetMet_dPhiMin4_temp = 999 ; 

		//for ( auto &tmp : Jets )
		//{
		//	llp_tree->jetNeutralHadronMultiplicity[llp_tree->nJets] = tmp.jetNeutralHadronMultiplicity;
		//	llp_tree->jetChargedHadronMultiplicity[llp_tree->nJets] = tmp.jetChargedHadronMultiplicity;
		//	llp_tree->jetMuonMultiplicity[llp_tree->nJets] = tmp.jetMuonMultiplicity;
		//	llp_tree->jetElectronMultiplicity[llp_tree->nJets] = tmp.jetElectronMultiplicity;
		//	llp_tree->jetPhotonMultiplicity[llp_tree->nJets] = tmp.jetPhotonMultiplicity;
		//	llp_tree->jetNeutralHadronEnergyFraction[llp_tree->nJets] = tmp.jetNeutralHadronEnergyFraction;
		//	llp_tree->jetChargedHadronEnergyFraction[llp_tree->nJets] = tmp.jetChargedHadronEnergyFraction;
		//	llp_tree->jetMuonEnergyFraction[llp_tree->nJets] = tmp.jetMuonEnergyFraction;
		//	llp_tree->jetElectronEnergyFraction[llp_tree->nJets] = tmp.jetElectronEnergyFraction;
		//	llp_tree->jetPhotonEnergyFraction[llp_tree->nJets] = tmp.jetPhotonEnergyFraction;
		//	llp_tree->jetCSV[llp_tree->nJets] = tmp.jetCSV;

		//	llp_tree->jetPt[llp_tree->nJets] = tmp.jet.Pt();
		//	llp_tree->jetEta[llp_tree->nJets] = tmp.jet.Eta();
		//	llp_tree->jetE[llp_tree->nJets] = tmp.jet.E();
		//	llp_tree->jetPhi[llp_tree->nJets] = tmp.jet.Phi();
		//	llp_tree->jetTime[llp_tree->nJets] = tmp.time;
		//	llp_tree->ecalNRechits[llp_tree->nJets] = tmp.ecalNRechits;
		//	llp_tree->ecalRechitE[llp_tree->nJets] = tmp.ecalRechitE;

		//	llp_tree->jetAlphaMax[llp_tree->nJets] = tmp.jetAlphaMax;
		//	llp_tree->jetBetaMax[llp_tree->nJets] = tmp.jetBetaMax;
		//	llp_tree->jetGammaMax[llp_tree->nJets] = tmp.jetGammaMax;
		//	llp_tree->jetGammaMax_Hadronic[llp_tree->nJets] = tmp.jetGammaMax_Hadronic;
		//	llp_tree->jetGammaMax_EM[llp_tree->nJets] = tmp.jetGammaMax_EM;
		//	llp_tree->jetGammaMax_ET[llp_tree->nJets] = tmp.jetGammaMax_ET;
		//	llp_tree->jetPtAllPVTracks[llp_tree->nJets] = tmp.jetPtAllPVTracks;
		//	llp_tree->jetPtAllTracks[llp_tree->nJets] = tmp.jetPtAllTracks;
		//	llp_tree->jetMinDeltaRAllTracks[llp_tree->nJets] = tmp.jetMinDeltaRAllTracks;
		//	llp_tree->jetMinDeltaRPVTracks[llp_tree->nJets] = tmp.jetMinDeltaRPVTracks;

		//	llp_tree->jetChargedEMEnergyFraction[llp_tree->nJets] = tmp.jetChargedEMEnergyFraction;
		//	llp_tree->jetNeutralEMEnergyFraction[llp_tree->nJets] = tmp.jetNeutralEMEnergyFraction;

		//	llp_tree->jetEcalE[llp_tree->nJets] = tmp.jet.E()*(tmp.jetElectronEnergyFraction+tmp.jetPhotonEnergyFraction);
		//	llp_tree->jetHcalE[llp_tree->nJets] = tmp.jet.E()*(tmp.jetNeutralHadronEnergyFraction+tmp.jetChargedHadronEnergyFraction);

		//	llp_tree->jetNVertexTracks[llp_tree->nJets] = tmp.jetNVertexTracks;
		//	llp_tree->jetNSelectedTracks[llp_tree->nJets] = tmp.jetNSelectedTracks;
		//	llp_tree->jetDRSVJet[llp_tree->nJets] = tmp.jetDRSVJet;
		//	llp_tree->jetSVMass[llp_tree->nJets] = tmp.jetSVMass;

		//	llp_tree->jetChargedMultiplicity[llp_tree->nJets] = tmp.jetChargedMultiplicity;
		//	llp_tree->jetNPixelHitsMedian[llp_tree->nJets] = tmp.jetNPixelHitsMedian;
		//	llp_tree->jetNHitsMedian[llp_tree->nJets] = tmp.jetNHitsMedian;

		//	//ecal rechits
		//	llp_tree->jetNRecHitsEcal[llp_tree->nJets] = tmp.jetNRecHitsEcal;
		//	llp_tree->jetEnergyRecHitsEcal[llp_tree->nJets] = tmp.jetEnergyRecHitsEcal;
		//	llp_tree->jetTimeRecHitsEcal[llp_tree->nJets] = tmp.jetTimeRecHitsEcal;

		//	//hcal hbhe rechits
		//	llp_tree->jetNRecHitsHcal[llp_tree->nJets] = tmp.jetNRecHitsHcal;
		//	llp_tree->jetEnergyRecHitsHcal[llp_tree->nJets] = tmp.jetEnergyRecHitsHcal;
		//	llp_tree->jetTimeRecHitsHcal[llp_tree->nJets] = tmp.jetTimeRecHitsHcal;

		//	llp_tree->jet_sig_et1[llp_tree->nJets] = tmp.jetsig1EB;
		//	llp_tree->jet_sig_et2[llp_tree->nJets] = tmp.jetsig2EB;
		//	llp_tree->jet_pt_deb[llp_tree->nJets] = tmp.jetptDEB;
		//	llp_tree->jet_sig_pt1[llp_tree->nJets] = tmp.jetsig1PF;
		//	llp_tree->jet_sig_pt2[llp_tree->nJets] = tmp.jetsig2PF;
		//	llp_tree->jet_pt_dpf[llp_tree->nJets] = tmp.jetptDPF;
		//}
		//llp_tree->HT = ht;
		//llp_tree->nCHSJets_in_HEM = nAK4Jets_in_HEM;
		//llp_tree->nCHSJets_in_HEM_eta_2p4 = nAK4Jets_in_HEM_eta_2p4;
		//llp_tree->nCHSJets_in_HEM_eta_2p5 = nAK4Jets_in_HEM_eta_2p5;
		//llp_tree->nCHSJets_in_HEM_eta_3 = nAK4Jets_in_HEM_eta_3;
		//llp_tree->nCHSJets_in_HEM_pt_20 = nAK4Jets_in_HEM_pt_20;
		//llp_tree->nCHSJets_in_HEM_pt_20_eta_2p4 = nAK4Jets_in_HEM_pt_20_eta_2p4;
		//llp_tree->nCHSJets_in_HEM_pt_20_eta_2p5 = nAK4Jets_in_HEM_pt_20_eta_2p5;
		//llp_tree->nCHSJets_in_HEM_pt_20_eta_3 = nAK4Jets_in_HEM_pt_20_eta_3;
		//llp_tree->nCHSJets_in_HEM_pt_30 = nAK4Jets_in_HEM_pt_30;
		//llp_tree->nCHSJets_in_HEM_pt_30_eta_2p4 = nAK4Jets_in_HEM_pt_30_eta_2p4;
		//llp_tree->nCHSJets_in_HEM_pt_30_eta_2p5 = nAK4Jets_in_HEM_pt_30_eta_2p5;
		//llp_tree->nCHSJets_in_HEM_pt_30_eta_3 = nAK4Jets_in_HEM_pt_30_eta_3;




		sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

		if (Jets.size()>0)
		{
			llp_tree->jetMet_dPhi = RazorAnalyzerLLP::deltaPhi(jetPhi[0],metType1Phi);
			//TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM( metType1Pt, 0, metType1Phi, 0 );
			TLorentzVector jet0 = makeTLorentzVectorPtEtaPhiM( jetPt[0], 0, jetPhi[0], 0 );
			llp_tree->jetMet_dPhiStar = RazorAnalyzerLLP::deltaPhi(jetPhi[0],  (t1PFMET+jet0).Phi() );
		}
		else{
			llp_tree->jetMet_dPhi = -999.;
			llp_tree->jetMet_dPhiStar = -999.;
		}

		float jetMet_dPhiMin_temp = 999 ; 
		float jetMet_dPhiStarMin_temp = 999 ; 
		float jetMet_dPhiMin4_temp = 999 ; 

		for ( auto &tmp : Jets )
		{
			llp_tree->jetNeutralHadronMultiplicity[llp_tree->nJets] = tmp.jetNeutralHadronMultiplicity;
			llp_tree->jetChargedHadronMultiplicity[llp_tree->nJets] = tmp.jetChargedHadronMultiplicity;
			llp_tree->jetMuonMultiplicity[llp_tree->nJets] = tmp.jetMuonMultiplicity;
			llp_tree->jetElectronMultiplicity[llp_tree->nJets] = tmp.jetElectronMultiplicity;
			llp_tree->jetPhotonMultiplicity[llp_tree->nJets] = tmp.jetPhotonMultiplicity;
			llp_tree->jetNeutralHadronEnergyFraction[llp_tree->nJets] = tmp.jetNeutralHadronEnergyFraction;
			llp_tree->jetChargedHadronEnergyFraction[llp_tree->nJets] = tmp.jetChargedHadronEnergyFraction;
			llp_tree->jetMuonEnergyFraction[llp_tree->nJets] = tmp.jetMuonEnergyFraction;
			llp_tree->jetElectronEnergyFraction[llp_tree->nJets] = tmp.jetElectronEnergyFraction;
			llp_tree->jetPhotonEnergyFraction[llp_tree->nJets] = tmp.jetPhotonEnergyFraction;
			llp_tree->jetCSV[llp_tree->nJets] = tmp.jetCSV;

			llp_tree->jetPt[llp_tree->nJets] = tmp.jet.Pt();
			llp_tree->jetEta[llp_tree->nJets] = tmp.jet.Eta();
			llp_tree->jetE[llp_tree->nJets] = tmp.jet.E();
			llp_tree->jetPhi[llp_tree->nJets] = tmp.jet.Phi();
			llp_tree->jetTime[llp_tree->nJets] = tmp.time;
			llp_tree->ecalNRechits[llp_tree->nJets] = tmp.ecalNRechits;
			llp_tree->ecalRechitE[llp_tree->nJets] = tmp.ecalRechitE;

			llp_tree->jetAlphaMax[llp_tree->nJets] = tmp.jetAlphaMax;
			llp_tree->jetBetaMax[llp_tree->nJets] = tmp.jetBetaMax;
			llp_tree->jetGammaMax[llp_tree->nJets] = tmp.jetGammaMax;
			llp_tree->jetGammaMax_Hadronic[llp_tree->nJets] = tmp.jetGammaMax_Hadronic;
			llp_tree->jetGammaMax_EM[llp_tree->nJets] = tmp.jetGammaMax_EM;
			llp_tree->jetGammaMax_ET[llp_tree->nJets] = tmp.jetGammaMax_ET;
			llp_tree->jetPtAllPVTracks[llp_tree->nJets] = tmp.jetPtAllPVTracks;
			llp_tree->jetPtAllTracks[llp_tree->nJets] = tmp.jetPtAllTracks;
			llp_tree->jetMinDeltaRAllTracks[llp_tree->nJets] = tmp.jetMinDeltaRAllTracks;
			llp_tree->jetMinDeltaRPVTracks[llp_tree->nJets] = tmp.jetMinDeltaRPVTracks;

			llp_tree->jetChargedEMEnergyFraction[llp_tree->nJets] = tmp.jetChargedEMEnergyFraction;
			llp_tree->jetNeutralEMEnergyFraction[llp_tree->nJets] = tmp.jetNeutralEMEnergyFraction;

			llp_tree->jetEcalE[llp_tree->nJets] = tmp.jet.E()*(tmp.jetElectronEnergyFraction+tmp.jetPhotonEnergyFraction);
			llp_tree->jetHcalE[llp_tree->nJets] = tmp.jet.E()*(tmp.jetNeutralHadronEnergyFraction+tmp.jetChargedHadronEnergyFraction);

			llp_tree->jetNVertexTracks[llp_tree->nJets] = tmp.jetNVertexTracks;
			llp_tree->jetNSelectedTracks[llp_tree->nJets] = tmp.jetNSelectedTracks;
			llp_tree->jetDRSVJet[llp_tree->nJets] = tmp.jetDRSVJet;
			llp_tree->jetSVMass[llp_tree->nJets] = tmp.jetSVMass;

			llp_tree->jetChargedMultiplicity[llp_tree->nJets] = tmp.jetChargedMultiplicity;
			llp_tree->jetNPixelHitsMedian[llp_tree->nJets] = tmp.jetNPixelHitsMedian;
			llp_tree->jetNHitsMedian[llp_tree->nJets] = tmp.jetNHitsMedian;

			//ecal rechits
			llp_tree->jetNRecHitsEcal[llp_tree->nJets] = tmp.jetNRecHitsEcal;
			llp_tree->jetEnergyRecHitsEcal[llp_tree->nJets] = tmp.jetEnergyRecHitsEcal;
			llp_tree->jetTimeRecHitsEcal[llp_tree->nJets] = tmp.jetTimeRecHitsEcal;
			llp_tree->jetTimeRmsEcal[llp_tree->nJets] = tmp.jetTimeRmsEcal;
			llp_tree->jetTimeRmsEcal_fix[llp_tree->nJets] = tmp.jetTimeRmsEcal_fix;

			//hcal hbhe rechits
			llp_tree->jetNRecHitsHcal[llp_tree->nJets] = tmp.jetNRecHitsHcal;
			llp_tree->jetEnergyRecHitsHcal[llp_tree->nJets] = tmp.jetEnergyRecHitsHcal;
			llp_tree->jetTimeRecHitsHcal[llp_tree->nJets] = tmp.jetTimeRecHitsHcal;

			llp_tree->jet_sig_et1[llp_tree->nJets] = tmp.jetsig1EB;
			llp_tree->jet_sig_et2[llp_tree->nJets] = tmp.jetsig2EB;
			llp_tree->jet_pt_deb[llp_tree->nJets] = tmp.jetptDEB;
			llp_tree->jet_sig_pt1[llp_tree->nJets] = tmp.jetsig1PF;
			llp_tree->jet_sig_pt2[llp_tree->nJets] = tmp.jetsig2PF;
			llp_tree->jet_pt_dpf[llp_tree->nJets] = tmp.jetptDPF;

			//llp_tree->jetDNNScoreV1[llp_tree->nJets] = tmp.dnn_score_v1;
			//llp_tree->jetDNNScore[llp_tree->nJets] = tmp.dnn_score;
			llp_tree->jetDNNScoreV3[llp_tree->nJets] = tmp.dnn_score_v3;
			llp_tree->jetDNNScoreV3miniAOD[llp_tree->nJets] = tmp.dnn_score_v3_miniAOD;

			llp_tree->jetCscEF[llp_tree->nJets] = tmp.jetCscEF;
			llp_tree->jetDtEF[llp_tree->nJets] = tmp.jetDtEF;

			//std::cout << "jetEta " << tmp.jet.Eta() << std::endl;
			//std::cout << "jetEta " << llp_tree->jetEta[llp_tree->nJets] << std::endl;

			if(jetMet_dPhiMin4_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi)) && llp_tree->nJets < 4)
			{
				jetMet_dPhiMin4_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi));
			}
			if (jetMet_dPhiMin_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi)))
			{
				jetMet_dPhiMin_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi));
			}     
			TLorentzVector jet_temp = makeTLorentzVectorPtEtaPhiM( tmp.jet.Pt(), 0, tmp.jet.Phi(), 0 );
			if (jetMet_dPhiStarMin_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(), (t1PFMET+jet_temp).Phi() )))
			{
				jetMet_dPhiStarMin_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(), (t1PFMET+jet_temp).Phi() ));
			}     

			//if(tmp.jet.Pt()>100. && fabs(tmp.jet.Eta())>2.25 && fabs(tmp.jet.Eta())<3.0) llp_tree->IsPrefiringAffected = true;

			//Int_t maxSumJetTrkPt_pvid[N_MAX_JETS];
			//Float_t maxSumJetTrkPt_pv[N_MAX_JETS];
			llp_tree->maxSumJetTrkPt_pvid[llp_tree->nJets] = tmp.maxpvid;
			llp_tree->maxSumJetTrkPt_svid[llp_tree->nJets] = tmp.maxsvid;
			llp_tree->maxSumJetTrkPt_pv[llp_tree->nJets] = tmp.maxsvpt;
			llp_tree->maxSumJetTrkPt_sv[llp_tree->nJets] = tmp.maxsvpt;


			llp_tree->nJets++;
		}

		//fatjets
		std::vector<fatjets> FatJets;
		for(int i = 0; i < nFatJets; i++)
		{

			//------------------------------------------------------------
			//exclude selected muons and electrons from the jet collection
			//------------------------------------------------------------
			double deltaR = -1;
			for(auto& lep : Leptons){
				double thisDR = RazorAnalyzerLLP::deltaR(fatJetEta[i],fatJetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
				if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
			}
			if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

			//------------------------------------------------------------
			//exclude selected taus from the jet collection
			//------------------------------------------------------------
			double deltaR_tau = -1;
			for(auto& tau : Taus){
				double thisDR_tau = RazorAnalyzerLLP::deltaR(fatJetEta[i],fatJetPhi[i],tau.tau.Eta(),tau.tau.Phi());
				if(deltaR_tau < 0 || thisDR_tau < deltaR_tau) deltaR_tau = thisDR_tau;
			}
			if(deltaR_tau > 0 && deltaR_tau < 0.4) continue; //jet matches a selected tau

			//------------------------------------------------------------
			//exclude selected photons from the jet collection
			//------------------------------------------------------------
			double deltaR_pho = -1;
			for(auto& pho : Photons){
				double thisDR_pho = RazorAnalyzerLLP::deltaR(fatJetEta[i],fatJetPhi[i],pho.photon.Eta(),pho.photon.Phi());
				if(deltaR_pho < 0 || thisDR_pho < deltaR_pho) deltaR_pho = thisDR_pho;
			}
			if(deltaR_pho > 0 && deltaR_pho < 0.4) continue; //jet matches a selected photon

			//------------------------------------------------------------
			//Apply Jet Energy and Resolution Corrections
			//------------------------------------------------------------

			TLorentzVector thisFatJet = makeTLorentzVector( fatJetPt[i], fatJetEta[i], fatJetPhi[i], fatJetE[i] );

			if( thisFatJet.Pt() < 30 ) continue;//According to the April 1st 2015 AN
			if( fabs( thisFatJet.Eta() ) >= 1.48 ) continue;

			//std::cout <<fatJetPt[i] << "," << fatJetE[i] <<  "," << fatJetEta[i] << std::endl;

			// ***********************************
			//Compute Rechit Quantities
			// ***********************************
			double fatjetEnergyRecHitsECAL = 0;
			double fatjetEnergyRecHitsHCAL = 0;
			double fatjetTimeRecHitsECAL = -100;
			double fatjetTimeRecHitsHCAL = -100;
			double tmpFatJetTimeEnergyRecHitsECAL = 0;
			double tmpFatJetTimeEnergyRecHitsHCAL = 0;
			int fatjetNRecHitsECAL = 0;
			int fatjetNRecHitsHCAL = 0;

			//Loop over ECAL rechits
			for (int q=0; q < nRechits; q++) {
				if (ecalRechit_E[q] <= 0.5) continue;
				double tmpDR = RazorAnalyzerLLP::deltaR(thisFatJet.Eta(), thisFatJet.Phi(), ecalRechit_Eta[q], ecalRechit_Phi[q]);
				if (tmpDR > 0.8) continue;			  
				if (ecalRechit_kSaturatedflag[q] || 
						ecalRechit_kLeadingEdgeRecoveredflag[q] || 
						ecalRechit_kPoorRecoflag[q] ||
						ecalRechit_kWeirdflag[q] || 
						ecalRechit_kDiWeirdflag[q]) continue;
				if (ecalRechit_T_Error[q] < 0 || ecalRechit_T_Error[q] > 100) continue;
				if (abs(ecalRechit_T[q]) > 12.5) continue;
				if (abs(ecalRechit_Eta[q]) > 1.5) continue;

				//cout << "Rechit " << q << " : " << ecalRechit_E[q] << "\n";

				fatjetEnergyRecHitsECAL += ecalRechit_E[q];
				fatjetNRecHitsECAL++;
				tmpFatJetTimeEnergyRecHitsECAL += ecalRechit_E[q]*ecalRechit_T[q];
			}  

			if (fatjetEnergyRecHitsECAL > 0) {
				fatjetTimeRecHitsECAL = tmpFatJetTimeEnergyRecHitsECAL / fatjetEnergyRecHitsECAL;
			} else {
				fatjetEnergyRecHitsECAL = -1;
				fatjetTimeRecHitsECAL = -100;
			}

			//if (fatjetNRecHitsECAL == 0) {
			//  fatjetEnergyRecHitsECAL = -1;
			//}
			//if (isnan(fatjetRechitT[i]) || jetRechitE[i] == 0) {
			//  fatjetTimeRecHitsECAL = -100;
			//} else {
			//  fatjetTimeRecHitsECAL = jetRechitT[i];
			//}

			//Loop over HCAL rechits
			for (int q=0; q < nHBHERechits; q++) {
				if (hbheRechit_E[q] <= 0.1) continue;
				double tmpDR = RazorAnalyzerLLP::deltaR(thisFatJet.Eta(), thisFatJet.Phi(), hbheRechit_Eta[q], hbheRechit_Phi[q]);
				if (tmpDR > 0.8) continue;

				fatjetEnergyRecHitsHCAL += hbheRechit_E[q];
				fatjetNRecHitsHCAL++;
				tmpFatJetTimeEnergyRecHitsHCAL += hbheRechit_E[q]*hbheRechit_T[q];
			}

			if (fatjetEnergyRecHitsHCAL > 0) {
				fatjetTimeRecHitsHCAL = tmpFatJetTimeEnergyRecHitsHCAL / fatjetEnergyRecHitsHCAL;
			} else {
				fatjetEnergyRecHitsHCAL = -1;
				fatjetTimeRecHitsHCAL = -100;
			}


			fatjets tmpFatJet;
			tmpFatJet.fatjet    = thisFatJet;
			tmpFatJet.CorrectedPt   = fatJetCorrectedPt[i];

			//ecal rechits
			tmpFatJet.fatjetNRecHitsEcal = fatjetNRecHitsECAL;
			tmpFatJet.fatjetEnergyRecHitsEcal = fatjetEnergyRecHitsECAL;
			tmpFatJet.fatjetTimeRecHitsEcal = fatjetTimeRecHitsECAL;

			//hcal hbhe rechits
			tmpFatJet.fatjetNRecHitsHcal = fatjetNRecHitsHCAL;
			tmpFatJet.fatjetEnergyRecHitsHcal = fatjetEnergyRecHitsHCAL;
			tmpFatJet.fatjetTimeRecHitsHcal = fatjetTimeRecHitsHCAL;


			FatJets.push_back(tmpFatJet);
		}

		sort(FatJets.begin(), FatJets.end(), my_largest_pt_fatjet);
		for ( auto &tmp : FatJets )
		{

			llp_tree->fatJetPt[llp_tree->nFatJets] = tmp.fatjet.Pt();
			llp_tree->fatJetEta[llp_tree->nFatJets] = tmp.fatjet.Eta();
			llp_tree->fatJetE[llp_tree->nFatJets] = tmp.fatjet.E();
			llp_tree->fatJetPhi[llp_tree->nFatJets] = tmp.fatjet.Phi();
			llp_tree->fatJetCorrectedPt[llp_tree->nFatJets] = tmp.CorrectedPt;

			//ecal rechits
			llp_tree->fatjetNRecHitsEcal[llp_tree->nFatJets] = tmp.fatjetNRecHitsEcal;
			llp_tree->fatjetEnergyRecHitsEcal[llp_tree->nFatJets] = tmp.fatjetEnergyRecHitsEcal;
			llp_tree->fatjetTimeRecHitsEcal[llp_tree->nFatJets] = tmp.fatjetTimeRecHitsEcal;

			//hcal hbhe rechits
			llp_tree->fatjetNRecHitsHcal[llp_tree->nFatJets] = tmp.fatjetNRecHitsHcal;
			llp_tree->fatjetEnergyRecHitsHcal[llp_tree->nFatJets] = tmp.fatjetEnergyRecHitsHcal;
			llp_tree->fatjetTimeRecHitsHcal[llp_tree->nFatJets] = tmp.fatjetTimeRecHitsHcal;

			llp_tree->nFatJets++;
		}

		//match fatjet (above 250) to ak4 jets with dR 0.4
		//jetIn250AK8[i] = -1;
		for(int i = 0; i < llp_tree->nFatJets; i++){
			//cout << "fatjet i pt " << i << " : " << llp_tree->fatJetPt[i] << "\n";
			if(llp_tree->fatJetPt[i]<250) continue;
			for(int j = 0; j < llp_tree->nJets; j++){
				//cout << "ak4jet j pt " << j << " : " << llp_tree->jetPt[j] << "\n";
				double tmpdR = RazorAnalyzerLLP::deltaR(llp_tree->fatJetEta[i], llp_tree->fatJetPhi[i], llp_tree->jetEta[j], llp_tree->jetPhi[j]);
				//cout << "ak4jet j dR " << j << " : " << tmpdR << "\n";
				if(tmpdR<0.4){
					llp_tree->jetIn250AK8[j] = 1;
				}
				//cout << "ak4jet j bool " << j << " : " << llp_tree->jetIn250AK8[j] << "\n";
			}
		}


		//gLLP grandaughters
		double ecal_radius = 129.0;
		//double hcal_radius = 179.0;
		double EB_z = 268.36447217; // 129*sinh(1.479)
		double EE_z = 298.5; //where Ecal Endcap starts in z direction
		if(!isData){
			llp_tree->genVertexX = genVertexX;
			llp_tree->genVertexY = genVertexY;
			llp_tree->genVertexZ = genVertexZ;
			llp_tree->genVertexT = genVertexT;

			int countZ =0;
			int countH =0;
			//gLLP daughters
			for(int i = 0; i <4; i++){
				llp_tree->gLLP_daughter_id[i] = gLLP_daughter_id[i];
				//label HH/HZ/ZZ
				if(gLLP_daughter_id[i]==23) countZ++;
				if(gLLP_daughter_id[i]==25) countH++;
				llp_tree->gLLP_daughter_mass[i] = gLLP_daughter_mass[i];
				llp_tree->gLLP_daughter_e[i] = gLLP_daughter_e[i];
				llp_tree->gLLP_daughter_pt[i] = gLLP_daughter_pt[i];
				llp_tree->gLLP_daughter_eta[i] = gLLP_daughter_eta[i];
				llp_tree->gLLP_daughter_phi[i] = gLLP_daughter_phi[i];

				TLorentzVector tmpdau = makeTLorentzVectorPtEtaPhiM(gLLP_daughter_pt[i], gLLP_daughter_eta[i], gLLP_daughter_phi[i], gLLP_daughter_mass[i]);
				int llp_index=0;
				if(i>1) llp_index=1;
				double radius = sqrt( pow(gLLP_decay_vertex_x[llp_index],2) + pow(gLLP_decay_vertex_y[llp_index],2) );
				llp_tree->gLLP_daughter_travel_time_EB[i] = (1./30.)*fabs(ecal_radius-radius)/(tmpdau.Pt()/tmpdau.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
				//Calculate dt from generation point to ECAL face
				double x_ecal = gLLP_decay_vertex_x[llp_index] + 30. * (tmpdau.Px()/tmpdau.E())*(llp_tree->gLLP_daughter_travel_time_EB[i]);
				double y_ecal = gLLP_decay_vertex_y[llp_index] + 30. * (tmpdau.Py()/tmpdau.E())*(llp_tree->gLLP_daughter_travel_time_EB[i]);
				double z_ecal = gLLP_decay_vertex_z[llp_index] + 30. * (tmpdau.Pz()/tmpdau.E())*(llp_tree->gLLP_daughter_travel_time_EB[i]);


				if( fabs(z_ecal) < EB_z && radius <= ecal_radius &&  fabs(gLLP_decay_vertex_z[llp_index]) < EE_z) {
					llp_tree->gLLP_daughter_photon_travel_time_EB[i] = (1./30) * sqrt(pow(ecal_radius,2)+pow(z_ecal,2));
					if(genVertexT==-999) genVertexT=0;
					llp_tree->gen_time_daughter_EB[i] = gLLP_travel_time[llp_index] + (llp_tree->gLLP_daughter_travel_time_EB[i]) - llp_tree->gLLP_daughter_photon_travel_time_EB[i] + genVertexT;

				} else {
					llp_tree->gLLP_daughter_travel_time_EB[i] = -666;
					llp_tree->gen_time_daughter_EB[i] = -666.;
					llp_tree->gLLP_daughter_photon_travel_time_EB[i] = -666.;
				}

				// Correction of eta and phi based on ecal points
				double phi = atan((y_ecal-genVertexY)/(x_ecal-genVertexX));
				if  (x_ecal < 0.0) {
					phi = TMath::Pi() + phi;
				}
				phi = deltaPhi(phi,0.0);

				double theta = atan(sqrt(pow(x_ecal-genVertexX,2)+pow(y_ecal-genVertexY,2))/abs(z_ecal-genVertexZ));
				double eta = -1.0*TMath::Sign(1.0, z_ecal-genVertexZ)*log(tan(theta/2));
				llp_tree->gLLP_daughter_eta_ecalcorr[i] = eta;
				llp_tree->gLLP_daughter_phi_ecalcorr[i] = phi;

				//match to AK4
				double min_delta_r = 666.;
				int match_jet_index = -666;
				for ( int i_jet = 0; i_jet < llp_tree->nJets; i_jet++ )
				{
					double current_delta_r = deltaR(llp_tree->gLLP_daughter_eta_ecalcorr[i], llp_tree->gLLP_daughter_phi_ecalcorr[i], llp_tree->jetEta[i_jet], llp_tree->jetPhi[i_jet]);
					if(_debug_match) std::cout << " i_jet " << i_jet << ", match_jet_index " << match_jet_index << ", min_delta_r " << min_delta_r <<", current_delta_r "<< current_delta_r << std::endl;
					if ( current_delta_r < min_delta_r )
					{
						min_delta_r = current_delta_r;
						match_jet_index = i_jet;
					}
					if(_debug_match) std::cout << " i_jet " << i_jet << ", match_jet_index " << match_jet_index << ", min_delta_r " << min_delta_r << std::endl;
				}//end matching to jets 

				if ( min_delta_r < 0.4 )
					//if ( min_delta_r < 20 )
				{
					llp_tree->gLLP_daughter_match_jet_index[i] = match_jet_index;
					llp_tree->gLLP_daughter_min_delta_r_match_jet[i] = min_delta_r;
					//llp_tree->matched[match_jet_index] = true;
					if(i<2) llp_tree->jet_matched_gLLP0_daughter[match_jet_index] = true;
					else llp_tree->jet_matched_gLLP1_daughter[match_jet_index] = true;
					if(_debug_match) std::cout << " i " << i << ", match_jet_index " << match_jet_index << ", min_delta_r " << min_delta_r << std::endl;
				}

				//match to AK8
				double fatjet_min_delta_r = 666.;
				int match_fatjet_index = -666;
				for ( int i_fatjet = 0; i_fatjet < llp_tree->nFatJets; i_fatjet++ )
				{
					double current_delta_r = deltaR(llp_tree->gLLP_daughter_eta_ecalcorr[i], llp_tree->gLLP_daughter_phi_ecalcorr[i], llp_tree->fatJetEta[i_fatjet], llp_tree->fatJetPhi[i_fatjet]);
					if(_debug_match) std::cout << " i_fatjet " << i_fatjet << ", match_fatjet_index " << match_fatjet_index << ", fatjet_min_delta_r " << fatjet_min_delta_r <<", current_delta_r "<< current_delta_r << std::endl;
					if ( current_delta_r < fatjet_min_delta_r )
					{
						fatjet_min_delta_r = current_delta_r;
						match_fatjet_index = i_fatjet;
					}
					if(_debug_match) std::cout << " i_fatjet " << i_fatjet << ", match_fatjet_index " << match_fatjet_index << ", fatjet_min_delta_r " << fatjet_min_delta_r << std::endl;
				}//end matching to jets 

				if ( fatjet_min_delta_r < 0.8 )
				{
					llp_tree->gLLP_daughter_match_fatjet_index[i] = match_fatjet_index;
					llp_tree->gLLP_daughter_min_delta_r_match_fatjet[i] = fatjet_min_delta_r;
					//llp_tree->matched[match_fatjet_index] = true;
					if(i<2) llp_tree->fatjet_matched_gLLP0_daughter[match_fatjet_index] = true;
					else llp_tree->fatjet_matched_gLLP1_daughter[match_fatjet_index] = true;
					if(_debug_match) std::cout << " i " << i << ", match_fatjet_index " << match_fatjet_index << ", fatjet_min_delta_r " << fatjet_min_delta_r << std::endl;
				}


			}//end of loop over daughters
			//label HH/HZ/ZZ
			if(countH==2) llp_tree->sig_label = 0; //HH
			if(countH==1 && countZ==1) llp_tree->sig_label = 1; //HZ
			if(countZ==2) llp_tree->sig_label = 2; //ZZ



			//gLLP grandaughters
			for(int i = 0; i <4; i++){
				llp_tree->gLLP_grandaughter_id[i] = gLLP_grandaughter_id[i];
				llp_tree->gLLP_grandaughter_mass[i] = gLLP_grandaughter_mass[i];
				llp_tree->gLLP_grandaughter_e[i] = gLLP_grandaughter_e[i];
				llp_tree->gLLP_grandaughter_pt[i] = gLLP_grandaughter_pt[i];
				llp_tree->gLLP_grandaughter_eta[i] = gLLP_grandaughter_eta[i];
				llp_tree->gLLP_grandaughter_phi[i] = gLLP_grandaughter_phi[i];
				//llp_tree->gLLP_grandaughter_eta_ecalcorr[i] = gLLP_grandaughter_eta_ecalcorr[i];
				//llp_tree->gLLP_grandaughter_phi_ecalcorr[i] = gLLP_grandaughter_phi_ecalcorr[i];

				TLorentzVector tmpdau = makeTLorentzVectorPtEtaPhiM(gLLP_grandaughter_pt[i], gLLP_grandaughter_eta[i], gLLP_grandaughter_phi[i], gLLP_grandaughter_mass[i]);
				int llp_index=0;
				if(i>1) llp_index=1;
				double radius = sqrt( pow(gLLP_decay_vertex_x[llp_index],2) + pow(gLLP_decay_vertex_y[llp_index],2) );
				llp_tree->gLLP_grandaughter_travel_time_EB[i] = (1./30.)*fabs(ecal_radius-radius)/(tmpdau.Pt()/tmpdau.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
				//Calculate dt from generation point to ECAL face
				double x_ecal = gLLP_decay_vertex_x[llp_index] + 30. * (tmpdau.Px()/tmpdau.E())*(llp_tree->gLLP_grandaughter_travel_time_EB[i]);
				double y_ecal = gLLP_decay_vertex_y[llp_index] + 30. * (tmpdau.Py()/tmpdau.E())*(llp_tree->gLLP_grandaughter_travel_time_EB[i]);
				double z_ecal = gLLP_decay_vertex_z[llp_index] + 30. * (tmpdau.Pz()/tmpdau.E())*(llp_tree->gLLP_grandaughter_travel_time_EB[i]);


				if( fabs(z_ecal) < EB_z && radius <= ecal_radius &&  fabs(gLLP_decay_vertex_z[llp_index]) < EE_z) {
					llp_tree->gLLP_grandaughter_photon_travel_time_EB[i] = (1./30) * sqrt(pow(ecal_radius,2)+pow(z_ecal,2));
					//llp_tree->photon_travel_time_dau_pv[i] = (1./30) * sqrt(pow(x_ecal-genVertexX,2) + pow(y_ecal-genVertexY,2) + pow(z_ecal-genVertexZ,2));
					//llp_tree->gen_time_dau_pv[i] =  gLLP_travel_time[llp_index] + (llp_tree->gLLP_grandaughter_travel_time_EB[i]) - photon_travel_time_dau_pv[i] + genVertexT;
					if(genVertexT==-999) genVertexT=0;
					llp_tree->gen_time_grandaughter_EB[i] = gLLP_travel_time[llp_index] + (llp_tree->gLLP_grandaughter_travel_time_EB[i]) - llp_tree->gLLP_grandaughter_photon_travel_time_EB[i] + genVertexT;

				} else {
					llp_tree->gLLP_grandaughter_travel_time_EB[i] = -666;
					//llp_tree->gen_time_dau_pv[i] = -666.;
					llp_tree->gen_time_grandaughter_EB[i] = -666.;
					llp_tree->gLLP_grandaughter_photon_travel_time_EB[i] = -666.;
					//llp_tree->photon_travel_time_dau_pv[i] = -666.;
				}

				// Correction of eta and phi based on ecal points
				double phi = atan((y_ecal-genVertexY)/(x_ecal-genVertexX));
				if  (x_ecal < 0.0) {
					phi = TMath::Pi() + phi;
				}
				phi = deltaPhi(phi,0.0);

				//std::cout << "phi " << phi 
				//	<< ", x_ecal " << x_ecal
				//	<< ", y_ecal " << y_ecal
				//	<< ", z_ecal " << z_ecal
				//	<< ", genVertexX " << genVertexX
				//	<< ", genVertexY " << genVertexY
				//	<< ", genVertexZ " << genVertexZ
				//	<< ", genVertexT " << genVertexT
				//	<< std::endl;
				double theta = atan(sqrt(pow(x_ecal-genVertexX,2)+pow(y_ecal-genVertexY,2))/abs(z_ecal-genVertexZ));
				double eta = -1.0*TMath::Sign(1.0, z_ecal-genVertexZ)*log(tan(theta/2));
				llp_tree->gLLP_grandaughter_eta_ecalcorr[i] = eta;
				llp_tree->gLLP_grandaughter_phi_ecalcorr[i] = phi;

				//llp_tree->gLLP_grandaughter_EB[i] = gLLP_grandaughter_EB[i];
				//llp_tree->gLLP_grandaughter_ETL[i] = gLLP_grandaughter_ETL[i];
				//llp_tree->gLLP_grandaughter_photon_travel_time_EB[i] = gLLP_grandaughter_photon_travel_time_EB[i];
				//llp_tree->gLLP_grandaughter_photon_travel_time_ETL[i] = gLLP_grandaughter_photon_travel_time_ETL[i];
				//llp_tree->gLLP_grandaughter_travel_time_EB[i] = gLLP_grandaughter_travel_time_EB[i];
				//llp_tree->gLLP_grandaughter_travel_time_ETL[i] = gLLP_grandaughter_travel_time_ETL[i];
				//llp_tree->gen_time_grandaughter_EB[i] = gen_time_grandaughter_EB[i];
				//llp_tree->gen_time_grandaughter_ETL[i] = gen_time_grandaughter_ETL[i];
				//llp_tree->gLLP_grandaughter_match_jet_index[i] = gLLP_grandaughter_match_jet_index[i];

				if(_debug_match) std::cout << "evt: "<<llp_tree->evtNum <<" n_jet " << llp_tree->nJets << std::endl; 
				double min_delta_r = 666.;
				int match_jet_index = -666;
				for ( int i_jet = 0; i_jet < llp_tree->nJets; i_jet++ )
				{
					double current_delta_r = deltaR(llp_tree->gLLP_grandaughter_eta_ecalcorr[i], llp_tree->gLLP_grandaughter_phi_ecalcorr[i], llp_tree->jetEta[i_jet], llp_tree->jetPhi[i_jet]);
					if(_debug_match) std::cout << " i_jet " << i_jet << ", match_jet_index " << match_jet_index << ", min_delta_r " << min_delta_r <<", current_delta_r "<< current_delta_r << std::endl;
					if ( current_delta_r < min_delta_r )
					{
						min_delta_r = current_delta_r;
						match_jet_index = i_jet;
					}
					if(_debug_match) std::cout << " i_jet " << i_jet << ", match_jet_index " << match_jet_index << ", min_delta_r " << min_delta_r << std::endl;
				}//end matching to jets 

				if ( min_delta_r < 0.4 )
					//if ( min_delta_r < 20 )
				{
					llp_tree->gLLP_grandaughter_match_jet_index[i] = match_jet_index;
					llp_tree->gLLP_grandaughter_min_delta_r_match_jet[i] = min_delta_r;
					//llp_tree->matched[match_jet_index] = true;
					if(i<2) llp_tree->jet_matched_gLLP0_grandaughter[match_jet_index] = true;
					else llp_tree->jet_matched_gLLP1_grandaughter[match_jet_index] = true;
					if(_debug_match) std::cout << " i " << i << ", match_jet_index " << match_jet_index << ", min_delta_r " << min_delta_r << std::endl;
				}
				//llp_tree->gLLP_grandaughter_min_delta_r_match_jet[i] = gLLP_grandaughter_min_delta_r_match_jet[i];

			}//loop over grandaughters

			//tagged?
			for(int i=0;i<4;i++)
			{
				int llp_index=0;
				if(i>1) llp_index=1;
				if(llp_tree->gLLP_grandaughter_match_jet_index[i]>-666)
				{
					int jet_index = llp_tree->gLLP_grandaughter_match_jet_index[i];
					if(llp_tree->jetTime[jet_index] > 0.09 
							&& llp_tree->jetMinDeltaRPVTracks[jet_index] > 0.06
							&& llp_tree->jetGammaMax_ET[jet_index] < 0.16
							&& llp_tree->jetChargedHadronEnergyFraction[jet_index] < 0.06)
					{
						llp_tree->gLLP_tagged[llp_index] = true;
					}
				}
			}

		}//finished gen part

		if( (label=="MR_JetHT") && llp_tree->nJets != 2 ) continue;
		if(label=="MR_JetHT") 
		{
			llp_tree->jet2_dPhi = abs(RazorAnalyzerLLP::deltaPhi(llp_tree->jetPhi[0], llp_tree->jetPhi[1]));
		}
		if( (label=="MR_PHO") && llp_tree->nJets != 1 ) continue;
		if(_debug) std::cout << "nJets in tree " << llp_tree->nJets << std::endl;
		llp_tree->jetMet_dPhiMin = jetMet_dPhiMin_temp;
		llp_tree->jetMet_dPhiStarMin = jetMet_dPhiStarMin_temp;
		llp_tree->jetMet_dPhiMin4 = jetMet_dPhiMin4_temp;


		if(label=="MR_PHO"){
			TLorentzVector photon0 = makeTLorentzVectorPtEtaPhiM(llp_tree->phoPt[0],llp_tree->phoEta[0], llp_tree->phoPhi[0], 0.);
			if (Jets.size()>0)
			{
				llp_tree->jetPho_dPhi = RazorAnalyzerLLP::deltaPhi(jetPhi[0], llp_tree->phoPhi[0]);
				TLorentzVector jet0 = makeTLorentzVectorPtEtaPhiM( jetPt[0], 0, jetPhi[0], 0 );
				llp_tree->jetPho_dPhiStar = RazorAnalyzerLLP::deltaPhi(jetPhi[0],  (photon0+jet0).Phi() );
			}
			else{
				llp_tree->jetPho_dPhi = -999.;
				llp_tree->jetPho_dPhiStar = -999.;
			}

			float jetPho_dPhiMin_temp = 999 ; 
			float jetPho_dPhiStarMin_temp = 999 ; 
			float jetPho_dPhiMin4_temp = 999 ; 


			int MR_PHO_nJets = 0;
			for ( auto &tmp : Jets )
			{

				if(jetPho_dPhiMin4_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),llp_tree->phoPhi[0])) && MR_PHO_nJets < 4)
				{
					jetPho_dPhiMin4_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),llp_tree->phoPhi[0]));
				}
				if (jetPho_dPhiMin_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),llp_tree->phoPhi[0])))
				{
					jetPho_dPhiMin_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),llp_tree->phoPhi[0]));
				}     
				TLorentzVector jet_temp = makeTLorentzVectorPtEtaPhiM( tmp.jet.Pt(), 0, tmp.jet.Phi(), 0 );
				if (jetPho_dPhiStarMin_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(), (photon0+jet_temp).Phi() )))
				{
					jetPho_dPhiStarMin_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(), (photon0+jet_temp).Phi() ));
				}     

				MR_PHO_nJets++;

			}
			if(_debug) std::cout << "nJets in tree " << MR_PHO_nJets << std::endl;
			llp_tree->jetPho_dPhiMin = jetPho_dPhiMin_temp;
			llp_tree->jetPho_dPhiStarMin = jetPho_dPhiStarMin_temp;
			llp_tree->jetPho_dPhiMin4 = jetPho_dPhiMin4_temp;


		}//MR_PHO dPhi


		//llp_tree->tree_->Fill();
		if(_debug) std::cout << "nJets in tree " << llp_tree->nJets << std::endl;

		if(!isData && signalScan)
		{
			pair<int,int> smsPair = make_pair(llp_tree->mX, llp_tree->ctau);
			Trees2D[smsPair]->Fill();
		}
		else if(!isData && (signalHZScan || signalHDecayScan) )
		{
			pair<int,int> smsPair = make_pair(llp_tree->mH, llp_tree->ctau);
			Trees2D[smsPair]->Fill();
		}
		else
		{
			llp_tree->tree_->Fill();
		}

	}//end fill


	if(!isData && signalScan)
	{
		for(auto &filePtr : Files2D)
		{
			cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
			filePtr.second->cd();
			Trees2D[filePtr.first]->Write();
			NEvents2D[filePtr.first]->Write("NEvents");
			NEventsMet2002D[filePtr.first]->Write("NEventsMet200");
			smsSumWeights2D[filePtr.first]->Write("SumWeights");
			smsSumScaleWeights2D[filePtr.first]->Write("SumScaleWeights");
			smsSumPdfWeights2D[filePtr.first]->Write("SumPdfWeights");
			filePtr.second->Close();

		}
	}
	else if(!isData && signalHZScan )
	{
		for(auto &filePtr : Files2D)
		{
			cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
			filePtr.second->cd();
			Trees2D[filePtr.first]->Write();
			NEvents2D[filePtr.first]->Write("NEvents");
			NEventsMet2002D[filePtr.first]->Write("NEventsMet200");
			smsSumWeights2D[filePtr.first]->Write("SumWeights");
			smsSumScaleWeights2D[filePtr.first]->Write("SumScaleWeights");
			smsSumPdfWeights2D[filePtr.first]->Write("SumPdfWeights");
			filePtr.second->Close();

		}
	}
	else if(!isData && signalHDecayScan)
	{
		for(auto &filePtr : Files2D)
		{
			cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
			filePtr.second->cd();
			Trees2D[filePtr.first]->Write();
			NEvents2D[filePtr.first]->Write("NEvents");
			NEventsMet2002D[filePtr.first]->Write("NEventsMet200");

			//H decay
			NEventsH2D[filePtr.first]->Write("NEventsH");	
			NEventsHbb2D[filePtr.first]->Write("NEventsHbb");	
			NEventsHcc2D[filePtr.first]->Write("NEventsHcc");	
			NEventsHgg2D[filePtr.first]->Write("NEventsHgg");	
			NEventsHWW2D[filePtr.first]->Write("NEventsHWW");	
			NEventsHZZ2D[filePtr.first]->Write("NEventsHZZ");	
			NEventsHmm2D[filePtr.first]->Write("NEventsHmm");	
			NEventsHtt2D[filePtr.first]->Write("NEventsHtt");	
			//HH
			NEventsHH2D[filePtr.first]->Write("NEventsHH");	
			NEventsHHbb2D[filePtr.first]->Write("NEventsHHbb");	
			NEventsHHcc2D[filePtr.first]->Write("NEventsHHcc");	
			NEventsHHgg2D[filePtr.first]->Write("NEventsHHgg");	
			NEventsHHWW2D[filePtr.first]->Write("NEventsHHWW");	
			NEventsHHZZ2D[filePtr.first]->Write("NEventsHHZZ");	
			NEventsHHmm2D[filePtr.first]->Write("NEventsHHmm");	
			NEventsHHtt2D[filePtr.first]->Write("NEventsHHtt");	
			//HZ
			NEventsZH2D[filePtr.first]->Write("NEventsZH");	
			NEventsZHbb2D[filePtr.first]->Write("NEventsZHbb");	
			NEventsZHcc2D[filePtr.first]->Write("NEventsZHcc");	
			NEventsZHgg2D[filePtr.first]->Write("NEventsZHgg");	
			NEventsZHWW2D[filePtr.first]->Write("NEventsZHWW");	
			NEventsZHZZ2D[filePtr.first]->Write("NEventsZHZZ");	
			NEventsZHmm2D[filePtr.first]->Write("NEventsZHmm");	
			NEventsZHtt2D[filePtr.first]->Write("NEventsZHtt");	

			smsSumWeights2D[filePtr.first]->Write("SumWeights");
			smsSumScaleWeights2D[filePtr.first]->Write("SumScaleWeights");
			smsSumPdfWeights2D[filePtr.first]->Write("SumPdfWeights");
			filePtr.second->Close();
		}
	} 
	else
	{
		cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
		cout << "Writing output trees..." << endl;		
		outFile->cd();
		llp_tree->tree_->Write();
		cout << "Writing NEvents..." << endl;		
		NEvents->Write();
		NEventsMet200->Write();
		cout << "Writing SumWeights..." << endl;		
		SumWeights->Write();
		cout << "SumWeights " << SumWeights->GetBinContent(1) << "\n";
		cout << "Writing SumScaleWeights..." << endl;		
    		SumScaleWeights->Write();
		cout << "SumScaleWeights " << SumScaleWeights->GetBinContent(1) << "\n";
		cout << "Writing SumPdfWeights..." << endl;		
    		SumPdfWeights->Write();
		cout << "SumPdfWeights " << SumPdfWeights->GetBinContent(1) << "\n";
		//cout << "Writing SumAlphasWeights..." << endl;		
    		//SumAlphasWeights->Write();
		//cout << "SumAlphasWeights " << SumAlphasWeights->GetBinContent(1) << "\n";
		//outFile->Write();
		outFile->Close();
	}

	//if (helper) delete helper; //for some reason this is causing crashes. something is 
	//                             not done corrector in the RazorHelper destructor
}
