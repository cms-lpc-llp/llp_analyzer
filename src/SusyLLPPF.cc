#include "LLPAnalysis/llpAnalyzer/interface/SusyLLPPF.h"
#include "LLPAnalysis/llpAnalyzer/interface/RazorHelper.h"
#include "LLPAnalysis/llpAnalyzer/interface/SusyLLPPFTree.h"
#include "LLPAnalysis/llpAnalyzer/interface/JetCorrectorParameters.h"
#include "LLPAnalysis/llpAnalyzer/interface/JetCorrectionUncertainty.h"
#include "LLPAnalysis/llpAnalyzer/interface/BTagCalibrationStandalone.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"

#define dnn_thre 0.996
#define dr_obj_jet 0.4
#define dr_obj_tau 0.5
#define dr_obj_emp 0.4 //e mu pho
//MetAna setting dr_obj_tau/emp 0.3

#define _debug 0
#define _debug_tau 0
#define _debug_pho 0
#define _debug_ee 0
#define _debug_nj 0
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
#define _debug_pre 1

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



struct ak4jets
{
	TLorentzVector jet;
	bool PassFail;
	float dnn_score_v3;
};

//pt comparison
//not used so far
struct greater_than_pt
{
	inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
};

//lepton highest pt comparator
struct largest_pt_pf_lep
{
	inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt_pf_lep;

//jet highest pt comparator
struct largest_pt_pf_jet
{
	inline bool operator() (const ak4jets& p1, const ak4jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_pf_jet;


//Analyze
void SusyLLPPF::Analyze(bool isData, int options, string outputfilename, string analysisTag, string process)
{
	//initialization: create one TTree for each analysis box
	cout << "Initializing..." << endl;
	cout << "IsData = " << isData << "\n";
	cout << "options = " << options << "\n";

	//---------------------------
	//-----------NN Setup----------
	//---------------------------
	//-----------v3-----------
	std::string basePathV3 = std::string(std::getenv("CMSSW_BASE")) + "/src/LLPAnalysis/llpAnalyzer/nn_inference/tagger_AK4_v3";
	
	std::string graphPathV3 = basePathV3 + "/graph.pb";
	std::string inputTensorNameV3 = "input_input";
	std::string outputTensorNameV3 = "FCN/output/Softmax";//"FCN/dense_4/Softmax";//or Softmax?
	
	// threading setup
	// to enable tensorflow-native multi-threading, change to "tensorflow" and increase nThreads
	std::string threadPoolV3 = "no_threads";
	int nThreadsV3 = 1;
	
	std::vector<std::string> inputFeaturesV3 = { "Jet_nTrackConstituents", "Jet_nSelectedTracks", "Jet_timeRecHitsEB", "Jet_eFracRecHitsEB", "Jet_nRecHitsEB", "Jet_sig1EB", "Jet_sig2EB", "Jet_ptDEB", "Jet_cHadEFrac", "Jet_nHadEFrac", "Jet_eleEFrac", "Jet_photonEFrac", "Jet_ptAllTracks", "Jet_ptAllPVTracks", "Jet_alphaMax", "Jet_betaMax", "Jet_gammaMax", "Jet_gammaMaxEM", "Jet_gammaMaxHadronic", "Jet_gammaMaxET", "Jet_minDeltaRAllTracks", "Jet_minDeltaRPVTracks",};

	int nInputsV3 = inputFeaturesV3.size();
	std::vector<float> inputValuesV3(nInputsV3);
	
	// setup TensorFlow objects
	tensorflow::setLogging();
	tensorflow::GraphDef* graphDefV3 = tensorflow::loadGraphDef(graphPathV3);
	tensorflow::Session* sessionV3 = tensorflow::createSession(graphDefV3, nThreadsV3);
		
	// register an input tensor (1 x nInputs) that is filled during the event loop
	tensorflow::Tensor inputTensorV3(tensorflow::DT_FLOAT, {1, nInputsV3});


	//---------------------------
	//-----------option----------
	//---------------------------
	int option;
	std::string label;
	bool signalScan;

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
		label = "MR_ZMM";
		//label = "wH";
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
		label = "MR_ZEE";
		cout << "process label: " << label << "\n";
	}
	else{
		cout << "What process it is? Label not defined. \n";
	}

	//UNIT'S DIGIT
	// signalScan option
	if(options%10==1){
		signalScan = true;
	}
	else{
		signalScan = false;
	}

	// DATA or MC
	if( isData )
	{
		std::cout << "[INFO]: running on data with label: " << label << " and option: " << option << " and signalScan is " << signalScan << std::endl;
	}
	else
	{
		std::cout << "[INFO]: running on MC with label: " << label << " and option: " << option << " and signalScan is " << signalScan << std::endl;
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
	if (label == "MR_SingleMuon" || label == "MR_ZMM"){
		NTrigger = 2;
		//NTrigger = 8;
	}
	if (label == "MR_SingleElectron" || label == "MR_ZEE"){
		NTrigger = 3;
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
	else if (label == "MR_SingleMuon" || label == "MR_ZMM"){
		wzId = 25;
		//trigger_paths[7] = 115; //HLT_IsoMu17_eta2p1
		//trigger_paths[6] = 120; //HLT_IsoMu18
		//trigger_paths[5] = 129; //HLT_IsoMu20 
		//trigger_paths[4] = 133; //HLT_IsoMu22
		//trigger_paths[3] = 135; //HLT_IsoMu24
		//trigger_paths[2] = 136; //HLT_IsoMu27
		//trigger_paths[1] = 644; //HLT_IsoMu24_eta2p1
		//trigger_paths[0] = 645; //HLT_IsoMu30
		trigger_paths[1] = 135; //HLT_IsoMu24
		trigger_paths[0] = 136; //HLT_IsoMu27
	}
	else if (label == "MR_SingleElectron" || label == "MR_ZEE"){
		wzId = 25;
		trigger_paths[2] = 87; //HLT_Ele32_WPTight_Gsf
		trigger_paths[1] = 88; //HLT_Ele32_eta2p1_WPLoose_Gsf
		//87  HLT_Ele32_WPTight_Gsf
		trigger_paths[0] = 627; 
		//627  HLT_Ele35_WPTight_Gsf
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
	if (outfilename == "") outfilename = "SusyLLPPFTree.root";
	TFile *outFile;
	if(!signalScan) outFile = new TFile(outfilename.c_str(), "RECREATE");

	SusyLLPPFTree *llp_tree = new SusyLLPPFTree;
	llp_tree->CreateTree();
	llp_tree->tree_->SetAutoFlush(0);
	llp_tree->InitTree();

	//histogram containing total number of processed events (for normalization)
	TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

	//for signals, need one output file for each signal point
	map<pair<int,int>, TFile*> Files2D;
	map<pair<int,int>, TTree*> Trees2D;
	map<pair<int,int>, TH1F*> NEvents2D;

	//*************************************************************************
	//Look over Input File Events
	//*************************************************************************
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

		if (label =="bkg_wH"|| label == "bkg_zH" || label == "bkg_HH"){
			if (isData)
			{
				NEvents->Fill(1);
				llp_tree->weight = 1;
			}
			else
			{
				//NEvents->Fill(genWeight);
				//llp_tree->weight = genWeight;
				NEvents->Fill(1);
				llp_tree->weight = 1;
			}
		}
		//else if(label =="MR_EMU"|| label == "MR_PHO" || label == "MR_ZLL"){
		else if( (label.find("MR") != std::string::npos) ){
			NEvents->Fill(1);
		}
		else{
			//generatedEvents->Fill(1);
			llp_tree->weight = 1;
		}
		if(_debug) std::cout << "deb2 " << jentry << std::endl;
		//if(!isData){
			//if( eventNum% 2 ==0 ) continue;
		//}
		//event info
		llp_tree->runNum = runNum;
		//llp_tree->lumiSec = lumiNum;
		llp_tree->evtNum = eventNum;
		if(_debug) std::cout << "deb3 " << jentry << std::endl;
		if(_debug) std::cout << "nBunchXing " << nBunchXing << std::endl;
		if (label == "zH" || label == "wH" ){
			NEvents->Fill(1);
			bool wzFlag = false;
			for (int i=0; i < nGenParticle; ++i)
			{
				if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] && abs(gParticleMotherId[i]) == wzId)
				{
					wzFlag = true;
				}

			}
			if ( wzFlag == false ) continue;

		}
		else if (label == "HH"){
			NEvents->Fill(1);
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



			if(_debug) std::cout << "nBunchXing " << nBunchXing << std::endl;
			int NPU=0;
			for (int i=0; i < nBunchXing; i++)
			{
				if (BunchXing[i] == 0)
				{
					//llp_tree->npu = int(nPUmean[i]);
					NPU = int(nPUmean[i]);
					if(_debug_npu) 
					{
						std::cout << "npu " << llp_tree->npu << std::endl;
						std::cout << "nPUmean[i] " << nPUmean[i] << std::endl;
					}
				}
			}

			llp_tree->pileupWeight = 1;
			llp_tree->pileupWeight = helper->getPileupWeight(NPU);
			//llp_tree->pileupWeight = helper->getPileupWeight(llp_tree->npu);
			//llp_tree->pileupWeightUp = helper->getPileupWeightUp(llp_tree->npu) / llp_tree->pileupWeight;
			//llp_tree->pileupWeightDown = helper->getPileupWeightDown(llp_tree->npu) / llp_tree->pileupWeight;
			if(_debug_npu) 
			{
				std::cout << "pileupWeightUp " << llp_tree->pileupWeightUp << std::endl;
				std::cout << "pileupWeightDown " << llp_tree->pileupWeightDown << std::endl;
			}
		}

		if(_debug_met) std::cout << "npu " << llp_tree->npu << std::endl;
		if(_debug && llp_tree->npu != 0 ) std::cout << "npu " << llp_tree->npu << std::endl;
		//get NPU
		//llp_tree->npv = nPV;
		//llp_tree->rho = fixedGridRhoFastjetAll;
		llp_tree->met = metType1Pt;
		if(_debug_met) std::cout << "met " << llp_tree->met << std::endl;
		//if( llp_tree->met < 120. ) continue;
		//if( llp_tree->met < 150. ) continue;
		if(_debug_lab) std::cout << "label " << label.c_str() << std::endl;
		if(_debug_lab) std::cout << "met " << llp_tree->met << std::endl;
		//if( (label.find("MR") == std::string::npos) && llp_tree->met < 200. ) continue;
		//if( (label=="MR_EMU") && llp_tree->met < 30. ) continue;
		//if( (label.find("MR_Single") != std::string::npos) && llp_tree->met < 40. ) continue;
		//if( (label.find("MR_ZLL") != std::string::npos) && isData && llp_tree->met >= 30. ) continue;
		//if( (label.find("MR_JetHT") != std::string::npos) && llp_tree->met >= 30. ) continue;
		//if( (label=="MR_PHO") && llp_tree->met >= 30. ) continue;
		if( (label.find("MR") == std::string::npos) && metType1Pt < 200. ) continue;
		if( (label=="MR_EMU") && metType1Pt < 30. ) continue;
		if( (label.find("MR_Single") != std::string::npos) && metType1Pt < 40. ) continue;
		if( (label=="MR_ZLL" || label=="MR_ZEE" || label=="MR_ZMM") && isData && metType1Pt >= 30. ) continue;
		if( (label.find("MR_JetHT") != std::string::npos) && metType1Pt >= 30. ) continue;
		if( (label=="MR_PHO") && metType1Pt >= 30. ) continue;
		if(_debug_lab) std::cout << "label " << label.c_str() << "passed "<< std::endl;
		if(_debug_lab) std::cout << "met " << llp_tree->met << "passed "<< std::endl;
		if(_debug_met) std::cout << "metType1Pt passed" << metType1Pt << std::endl;
		//llp_tree->metPhi = metType1Phi;
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

		if(!(Flag2_globalSuperTightHalo2016Filter && Flag2_BadPFMuonFilter && Flag2_EcalDeadCellTriggerPrimitiveFilter && Flag2_HBHENoiseFilter && Flag2_HBHEIsoNoiseFilter )) continue;
		if( (analysisTag.find("2016") == std::string::npos) && !Flag2_ecalBadCalibFilter) continue;
		if(isData && !Flag2_eeBadScFilter) continue;
		//Triggers
		if(_debug_trg) std::cout << "begin: 310 " << HLTDecision[310] << std::endl;
		if(_debug_trg) std::cout << "begin: 310 " << llp_tree->HLTDecision[310] << std::endl;
		//if( (label.find("MR") == std::string::npos) && !HLTDecision[310]) continue; 
		if( (label.find("MR") == std::string::npos) && !HLTDecision[467]) continue; 
		//if(!HLTDecision[467]) continue; 

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
		else if( (label.find("MR_EMU") != std::string::npos) || (label.find("MR_Single") != std::string::npos) || (label=="MR_ZLL") || (label=="MR_ZMM") || (label=="MR_ZEE")){
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
		if( (label.find("MR") != std::string::npos) && !triggered) continue; 

		if(_debug) std::cout << "passed trigger " << std::endl;
		if(_debug_pho||_debug_ee) std::cout << "passed trigger " << std::endl;

		

		//*************************************************************************
		//Start Object Selection
		//*************************************************************************
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
			if(fabs(muonEta[i]) > 2.4) continue;

			if(!muonIsLoose[i]) continue;
			if( (muon_chargedIso[i] + std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) )/muonPt[i] >= 0.25) continue;


			if( (label=="MR_SingleMuon" || label=="MR_ZMM") && (muonPt[i] < 30. || !muonIsTight[i] || (muon_chargedIso[i] + std::max(muon_neutralHadIso[i]      + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) )/muonPt[i] >= 0.15 )) continue; 
			if( label=="MR_EMU" && muonPt[i] < 30.) continue; 

			//remove overlaps
			bool overlap = false;
			for(auto& lep : Leptons)
			{
				if (RazorAnalyzerLLP::deltaR(muonEta[i],muonPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < dr_obj_emp) overlap = true;
				//if (RazorAnalyzerLLP::deltaR(muonEta[i],muonPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
			}
			if(overlap) continue;

			leptons tmpMuon;
			tmpMuon.lepton.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i], MU_MASS);
			tmpMuon.pdgId = 13 * -1 * muonCharge[i];

			Leptons.push_back(tmpMuon);
			llp_tree->nMuons++;
		}
		if( (label=="MR_SingleMuon") && llp_tree->nMuons != 1 ) continue;
		if( (label=="MR_ZMM") && llp_tree->nMuons != 2 ) continue;
		//if( (label=="MR_SingleElectron")  && llp_tree->nMuons != 0 ) continue;
		//if( (label=="MR_EMU")  && llp_tree->nMuons != 1 ) continue;
		if(_debug) std::cout << "nMuons " << llp_tree->nMuons << std::endl;
		//if( (label=="MR_PHO") && llp_tree->nMuons != 0 ) continue;
		//if( (label.find("MR_SingleElectron") != std::string::npos) && llp_tree->nMuons != 0 ) continue;
		//if( (label=="MR_SingleElectron" || label=="MR_JetHT")  && llp_tree->nMuons != 0 ) continue;
		//if( (label=="MR_JetHT")  && llp_tree->nMuons != 0 ) continue;

		if(_debug_ee) std::cout << "nMuons " << llp_tree->nMuons << std::endl;
		if(_debug_ee) std::cout << "nElectrons " << nElectrons << std::endl;
		//-------------------------------
		//Electrons
		//-------------------------------
		for( int i = 0; i < nElectrons; i++ )
		{
			if(_debug_sync||_debug_ee)
			{
				//std::cout << "nElectrons " << nElectrons << std::endl;
				std::cout << "iElectron " << i << ", Pt " << elePt[i] << ", Eta " << eleEta[i]<< ", Phi " << elePhi[i]<< ", E " << eleE[i] << std::endl;
				std::cout << "iElectron " << i << ", ele_passCutBasedIDVeto[i] " << ele_passCutBasedIDVeto[i]<< std::endl;
			}

			if(elePt[i] <= 10 ) continue;
			if(fabs(eleEta[i]) > 2.5) continue;
			if(!ele_passCutBasedIDVeto[i]) continue;


			if( (label=="MR_SingleElectron" || label=="MR_ZLL" || label=="MR_ZEE") && (elePt[i] < 37. || !ele_passCutBasedIDTight[i]) ) continue; 
			if( label=="MR_EMU" && elePt[i] < 30.) continue; 
			//remove overlaps
			bool overlap = false;
			for(auto& lep : Leptons)
			{
				if (RazorAnalyzerLLP::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < dr_obj_emp) overlap = true;
				//if (RazorAnalyzerLLP::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
			}
			if(overlap) continue;

			leptons tmpElectron;
			tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);

			Leptons.push_back(tmpElectron);
			llp_tree->nElectrons++;
		}
		if( (label=="MR_SingleElectron") && llp_tree->nElectrons != 1 ) continue;
		if( (label=="MR_ZEE") && llp_tree->nElectrons != 2 ) continue;
		if(_debug||_debug_ee) std::cout << "nElectrons " << llp_tree->nElectrons << std::endl;

		//if( (label=="MR_PHO") && llp_tree->nElectrons != 0 ) continue;
		//if( (label=="MR_SingleMuon")  && llp_tree->nElectrons != 0 ) continue;
		//if( (label=="MR_EMU") && llp_tree->nElectrons != 1 ) continue;
		//if( (label=="MR_JetHT") && llp_tree->nElectrons != 0 ) continue;
		//if( (label.find("MR_SingleMuon") != std::string::npos) && llp_tree->nElectrons != 0 ) continue;
		//if( (label.find("MR_EMU") != std::string::npos) && llp_tree->nElectrons != 1 ) continue;
		//if( (label.find("MR_JetHT") != std::string::npos) && llp_tree->nElectrons != 0 ) continue;


		if( (label=="MR_EMU") && !(llp_tree->nMuons==1 && llp_tree->nElectrons==1) ) continue;


		std::vector<taus> Taus;
		//-------------------------------
		//Taus
		//-------------------------------
		for( int i = 0; i < nTaus; i++ )
		{
			if(_debug_sync||_debug_tau)
			{
				std::cout << "nTaus " << nTaus << std::endl;
				std::cout << "iTau " << i << ", Pt " << tauPt[i] << ", Eta " << tauEta[i]<< ", Phi " << tauPhi[i]<< ", E " << tauE[i] << std::endl;
				std::cout << "iTau " << i << ", tau_IsLoose[i] " << tau_IsLoose[i]<< std::endl;
				std::cout << "iTau " << i << ", tau_ID[i] " << tau_ID[i]<< std::endl;
			}

			if(tauPt[i] <= 18 ) continue;
			if(fabs(tauEta[i]) > 2.3) continue;

			if(!tau_IsLoose[i]) continue;
			if(tau_ID[i]%2==0) continue;

			//remove overlaps
			bool overlap = false;
			for(auto& lep : Leptons)
			{
				if (RazorAnalyzerLLP::deltaR(tauEta[i],tauPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < dr_obj_tau) overlap = true;
				//if (RazorAnalyzerLLP::deltaR(tauEta[i],tauPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
			}
			if(overlap) continue;

			bool overlap_tau = false;
			for(auto& tau : Taus)
			{
				if (RazorAnalyzerLLP::deltaR(tauEta[i],tauPhi[i],tau.tau.Eta(),tau.tau.Phi()) < dr_obj_tau) overlap_tau = true;
				//if (RazorAnalyzerLLP::deltaR(tauEta[i],tauPhi[i],tau.tau.Eta(),tau.tau.Phi()) < 0.3) overlap_tau = true;
			}
			if(overlap_tau) continue;

			taus tmpTau;
			tmpTau.tau.SetPtEtaPhiM(tauPt[i],tauEta[i], tauPhi[i], TAU_MASS);

			Taus.push_back(tmpTau);
			llp_tree->nTaus++;
		}
		if(_debug) std::cout << "nTaus " << llp_tree->nTaus << std::endl;
		//if( (label.find("MR") == std::string::npos) && (llp_tree->nTaus != 0)) continue; 

		if(_debug_pho) std::cout << "nMuons " << llp_tree->nMuons << std::endl;
		if(_debug_pho) std::cout << "nElectrons " << llp_tree->nElectrons << std::endl;
		if(_debug_pho) std::cout << "nTaus " << llp_tree->nTaus << std::endl;
		if(_debug_pho) std::cout << "nPhotons " << llp_tree->nPhotons << std::endl;
		std::vector<photons> Photons;
		//-------------------------------
		//Photons
		//-------------------------------
		for( int i = 0; i < nPhotons; i++ )
		{
			if(_debug_sync||_debug_pho)
			{
				std::cout << "nPhotons " << nPhotons << std::endl;
				std::cout << "iPhoton " << i << ", Pt " << phoPt[i] << ", Eta " << phoEta[i]<< ", Phi " << phoPhi[i]<< ", E " << phoE[i] << std::endl;
				std::cout << "iPhoton " << i << ", pho_passCutBasedIDLoose[i] " << pho_passCutBasedIDLoose[i]<< std::endl;
			}

			if(phoPt[i] <= 15 ) continue;
			if( label=="MR_SinglePhoton" && phoPt[i] < 25.) continue; 
			//if( (label.find("MR") != std::string::npos) && phoPt[i] < 25.) continue; 
			if(fabs(phoEta[i]) > 2.5) continue;

			if(!pho_passCutBasedIDLoose[i]) continue;

			//remove overlaps
			bool overlap = false;
			for(auto& lep : Leptons)
			{
				if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < dr_obj_emp) overlap = true;
				//if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
			}
			if(overlap) continue;

			bool overlap_tau = false;
			for(auto& tau : Taus)
			{
				//if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],tau.tau.Eta(),tau.tau.Phi()) < dr_obj_emp) overlap_tau = true;
				if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],tau.tau.Eta(),tau.tau.Phi()) < dr_obj_tau) overlap_tau = true;
				//if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],tau.tau.Eta(),tau.tau.Phi()) < 0.3) overlap_tau = true;
			}
			if(overlap_tau) continue;

			bool overlap_pho = false;
			for(auto& pho : Photons)
			{
				if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],pho.photon.Eta(),pho.photon.Phi()) < dr_obj_emp) overlap_pho = true;
				//if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],pho.photon.Eta(),pho.photon.Phi()) < 0.3) overlap_pho = true;
			}
			if(overlap_pho) continue;

			photons tmpPhoton;
			tmpPhoton.photon.SetPtEtaPhiM(phoPt[i],phoEta[i], phoPhi[i], 0.);

			Photons.push_back(tmpPhoton);
			llp_tree->nPhotons++;
		}
		if(_debug) std::cout << "nPhotons " << llp_tree->nPhotons << std::endl;
		if( (label=="MR_PHO") && llp_tree->nPhotons != 1 ) continue;

		//if( (label.find("MR_Single") != std::string::npos) && llp_tree->nPhotons != 0 ) continue;
		//if( (label.find("MR_EMU") != std::string::npos) && llp_tree->nPhotons != 0 ) continue;
		//if( (label.find("MR_ZLL") != std::string::npos) && llp_tree->nPhotons != 0 ) continue;
		//if( (label.find("MR_JetHT") != std::string::npos) && llp_tree->nPhotons != 0 ) continue;

		//-------------------------------
		//Leptons
		//-------------------------------
		TLorentzVector lepp4;
		sort(Leptons.begin(), Leptons.end(), my_largest_pt_pf_lep);
		for ( auto &tmp : Leptons )
		{
			// std::cout << "lepton pdg " << llp_tree->lepPdgId[llp_tree->nLeptons] << std::endl;
			llp_tree->nLeptons++;

			lepp4 += tmp.lepton;
		}
		float dPhi = RazorAnalyzerLLP::deltaPhi(llp_tree->metPhi, lepp4.Phi());
		//float MT = sqrt(2*(llp_tree->met)*lepp4.Pt()*(1-cos(dPhi)));
		llp_tree->MT = sqrt(2*(llp_tree->met)*lepp4.Pt()*(1-cos(dPhi)));
		if(_debug_lep) std::cout<<"MT "<<llp_tree->MT<<std::endl;
		if( (label=="MR_SingleMuon" || label=="MR_SingleElectron") && llp_tree->MT>=100 ) continue;
		//if (triggered) trig_lepId->Fill(1);
		if(_debug) std::cout << "nLeptons " << llp_tree->nLeptons << std::endl;
		//if( (label.find("MR_Single") != std::string::npos) && MT>=100 ) continue;
		//if( (label=="MR_PHO") && llp_tree->nLeptons != 0 ) continue;
		//if( (label=="MR_JetHT") && llp_tree->nLeptons != 0 ) continue;


		//-------------------------------
		// reconstruct Z
		//-------------------------------
		double ZMass = -999;
		//double ZPt = -999;
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
					//ZPt = (Leptons[i].lepton+Leptons[j].lepton).Pt();
					ZCandidate = Leptons[i].lepton+Leptons[j].lepton;
					foundZ = true;
				}
			}
		}
		if(_debug_ee) std::cout << "ZMass " << ZMass << std::endl;
		if(foundZ) llp_tree->ZMass = ZMass;
		if( (label=="MR_ZLL" || label=="MR_ZMM" || label=="MR_ZEE") && isData && !(foundZ && fabs(ZMass-Z_MASS) < 30.0) ) continue;
		//if( (label=="MR_SingleMuon") && !(MT<100 && llp_tree->nMuons==1 && llp_tree->nElectrons==0 && llp_tree->nTaus==0 && llp_tree->nPhotons==0) ) continue;
		//if( (label=="MR_SingleElectron") && !(MT<100 && llp_tree->nMuons==0 && llp_tree->nElectrons==1 && llp_tree->nTaus==0 && llp_tree->nPhotons==0) ) continue;
		//if( (label=="MR_EMU") && !( llp_tree->nMuons==1 && llp_tree->nElectrons==1 && llp_tree->nTaus==0 && llp_tree->nPhotons==0) ) continue;
		//if( (label=="MR_PHO") && !( llp_tree->nMuons==0 && llp_tree->nElectrons==0 && llp_tree->nTaus==0 && llp_tree->nPhotons==1) ) continue;

		//-----------------------------------------------
		//Select Jets
		//-----------------------------------------------
		//std::vector<double> jetPtVector;
		//std::vector<double> jetCISVVector;
		std::vector<ak4jets> AK4Jets;

		if(_debug) std::cout << "nJets " << nJets << std::endl;
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
		//std::vector<jets> Jets;
		if(_debug_nj) std::cout << "len Jets " << AK4Jets.size() << std::endl;
		if(_debug_nj) std::cout << "nJets " << nJets << std::endl;
		//Jets.clear();
		//Jets.empty();
		if(_debug_nj) std::cout << "len Jets " << AK4Jets.size() << std::endl;

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
			if(_debug_nj) std::cout << "len Jets " << AK4Jets.size() << std::endl;


			//------------------------------------------------------------
			//exclude selected muons and electrons from the jet collection
			//------------------------------------------------------------
			double deltaR = -1;
			for(auto& lep : Leptons){
				double thisDR = RazorAnalyzerLLP::deltaR(jetEta[i],jetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
				if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
			}
			if(deltaR > 0 && deltaR < dr_obj_jet) continue; //jet matches a selected lepton

			//------------------------------------------------------------
			//exclude selected taus from the jet collection
			//------------------------------------------------------------
			double deltaR_tau = -1;
			for(auto& tau : Taus){
				double thisDR_tau = RazorAnalyzerLLP::deltaR(jetEta[i],jetPhi[i],tau.tau.Eta(),tau.tau.Phi());
				if(deltaR_tau < 0 || thisDR_tau < deltaR_tau) deltaR_tau = thisDR_tau;
			}
			if(deltaR_tau > 0 && deltaR_tau < dr_obj_jet) continue; //jet matches a selected tau

			//------------------------------------------------------------
			//exclude selected photons from the jet collection
			//------------------------------------------------------------
			double deltaR_pho = -1;
			for(auto& pho : Photons){
				double thisDR_pho = RazorAnalyzerLLP::deltaR(jetEta[i],jetPhi[i],pho.photon.Eta(),pho.photon.Phi());
				if(deltaR_pho < 0 || thisDR_pho < deltaR_pho) deltaR_pho = thisDR_pho;
			}	
			if(deltaR_pho > 0 && deltaR_pho < dr_obj_jet) continue; //jet matches a selected photon

			//------------------------------------------------------------
			//Apply Jet Energy and Resolution Corrections
			//------------------------------------------------------------

			TLorentzVector thisJet = makeTLorentzVector( jetPt[i], jetEta[i], jetPhi[i], jetE[i] );
			if(_debug_nj) std::cout << "len Jets " << AK4Jets.size() << std::endl;

			if( thisJet.Pt() < 30 ) continue;//According to the April 1st 2015 AN
			if( fabs( thisJet.Eta() ) >= 1.48 ) continue;
			if(_debug_nj) std::cout << "len Jets " << AK4Jets.size() << std::endl;
			//if( jetRechitT[i]<=-1 ) continue;
			//if(_debug_dnn) std::cout << "Jet Time " << jetRechitT[i] << std::endl;
			double jetTimeRecHitsECAL = -100;
			if (isnan(jetRechitT[i]) || jetRechitE[i] == 0) {
			  jetTimeRecHitsECAL = -100;
			} else {
			  jetTimeRecHitsECAL = jetRechitT[i];
			}
			if(_debug_dnn) std::cout << "Jet Time " << jetTimeRecHitsECAL << std::endl;
			if( jetTimeRecHitsECAL<=-1 ) continue;
			if(_debug_dnn) std::cout << "Jet Time " << jetTimeRecHitsECAL << std::endl;
			if( jetMuonEnergyFraction[i]>=0.6 ) continue;
			if( jetElectronEnergyFraction[i]>=0.6 ) continue;
			if( jetPhotonEnergyFraction[i]>=0.8 ) continue;

			//************************************
			//Compute Rechit Quantities
			//************************************
			double jetEnergyRecHitsECAL = 0;
			double jetEnergyRecHitsHCAL = 0;
			//double jetTimeRecHitsECAL = -100;
			double jetTimeRecHitsHCAL = -100;
			double tmpJetTimeEnergyRecHitsHCAL = 0;
 			int jetNRecHitsECAL = 0;
 			int jetNRecHitsHCAL = 0;
			std::vector<double> ebrechitphi;
			std::vector<double> ebrechiteta;
			std::vector<double> ebrechitet;
			std::vector<double> ebrechitetsq;

			if(_debug_pf) std::cout << "this jet pt " << thisJet.Pt() << std::endl;

			//Loop over ECAL rechits
			for (int q=0; q < nRechits; q++) {
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
			}  
			double jetsig1EB(-1.),jetsig2EB(-1.);
			RazorAnalyzerLLP::jet_second_moments(ebrechitet,ebrechiteta,ebrechitphi,jetsig1EB,jetsig2EB);
			double jetptDEB = -1;
			if ( accumulate(ebrechitet.begin(),ebrechitet.end(),0) > 0) {
				jetptDEB = sqrt( accumulate(ebrechitetsq.begin(),ebrechitetsq.end(),0)) / accumulate(ebrechitet.begin(),ebrechitet.end(),0) ;
			}	
			if(_debug_dnn) std::cout << "this jet sig1EB " << jetsig1EB << std::endl;
			if(_debug_dnn) std::cout << "this jet sig2EB " << jetsig2EB << std::endl;
			if(_debug_dnn) std::cout << "this jet ptDEB " << jetptDEB << std::endl;
			if(_debug_pf) std::cout << "this jet ptDEB " << jetptDEB << std::endl;
			if (jetNRecHitsECAL == 0) {
			  jetEnergyRecHitsECAL = -1;
			  jetNRecHitsECAL = -1;
			}
			//if (isnan(jetRechitT[i]) || jetRechitE[i] == 0) {
			//  jetTimeRecHitsECAL = -100;
			//} else {
			//  jetTimeRecHitsECAL = jetRechitT[i];
			//}

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

			
			if(_debug_nj) std::cout << "len Jets " << AK4Jets.size() << std::endl;
			//************************************
			//Evaluate NN tagger
			//************************************
			inputValuesV3[0] = jetChargedHadronMultiplicity[i]+jetElectronMultiplicity[i]+jetMuonMultiplicity[i];
			inputValuesV3[1] = jetNSelectedTracks[i];
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
			inputValuesV3[17] = (jetGammaMax_EM[i] == -99 || jetPtAllTracks[i] == 0 ||jetE[i]*(jetPhotonEnergyFraction[i]+jetElectronEnergyFraction[i]) == 0 || isnan(jetGammaMax_EM[i])) ? -100 : jetGammaMax_EM[i];
			inputValuesV3[18] = (jetGammaMax_Hadronic[i] == -99 || jetPtAllTracks[i] == 0 || jetE[i]*(jetNeutralHadronEnergyFraction[i]+jetChargedHadronEnergyFraction[i]) == 0) ? -100 : jetGammaMax_Hadronic[i];
			inputValuesV3[19] = (jetGammaMax_ET[i] == -99 || jetPtAllTracks[i] == 0 || thisJet.Et() == 0) ? -100 : jetGammaMax_ET[i];
			inputValuesV3[20] = (jetMinDeltaRAllTracks[i] == -99 || jetMinDeltaRAllTracks[i] == 15) ? 999 : jetMinDeltaRAllTracks[i];
			inputValuesV3[21] = (jetMinDeltaRPVTracks[i] == -99 || jetMinDeltaRPVTracks[i] == 15) ? 999 : jetMinDeltaRPVTracks[i];
			//std::cout<< " input value 1: " << jetNSelectedTracks[i] <<std::endl;
			
			if(_debug_dnn) {
			  std::cout<< " =========E V E N T : "<<eventNum <<std::endl;
			 for(int rr=0;rr<22;rr++)
			{
			  std::cout<< " input value "<<rr<<" : " << inputValuesV3[rr] <<std::endl;
			}}
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
			if(_debug_nj) std::cout << "output value: " << outputValueV3 << std::endl;
			if(_debug_nj) std::cout << "len Jets " << AK4Jets.size() << std::endl;

			ak4jets tmpJet;
			tmpJet.jet    = thisJet;
			//tmpJet.dnn_score_v3 = 0;
			tmpJet.dnn_score_v3 = outputValueV3;
			bool dnn_wp_pass = 0;

			if(outputValueV3>dnn_thre)
			{
				dnn_wp_pass=1;
			}
			else
			{
				dnn_wp_pass=0;
			}

			tmpJet.PassFail = dnn_wp_pass;
			AK4Jets.push_back(tmpJet);

			if(_debug_nj) std::cout << "pushed into AK4Jets " << std::endl;
			if(_debug_nj) std::cout << "len AK4Jets " << AK4Jets.size() << std::endl;

		}
		llp_tree->jetMet_dPhiMin_eta_2p4 = jetMet_dPhiMin_eta_2p4_temp;
		llp_tree->jetMet_dPhiMin_eta_3 = jetMet_dPhiMin_eta_3_temp;
		llp_tree->jetMet_dPhiMin_eta_all = jetMet_dPhiMin_eta_all_temp;
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


		if(_debug_nj) std::cout << "len AK4Jets " << AK4Jets.size() << std::endl;
		sort(AK4Jets.begin(), AK4Jets.end(), my_largest_pt_pf_jet);
		if(_debug_nj)  std::cout << "sorted AK4Jets " << std::endl;

  		float jetMet_dPhiMin_temp = 999 ; 
		llp_tree->nJets = 0;
		for ( auto &tmp : AK4Jets )
		{
			llp_tree->jetPt[llp_tree->nJets] = tmp.jet.Pt();
			llp_tree->jetPhi[llp_tree->nJets] = tmp.jet.Phi();
			llp_tree->jetEta[llp_tree->nJets] = tmp.jet.Eta();
			llp_tree->jetPass[llp_tree->nJets] = tmp.PassFail;
			llp_tree->jetDNNScoreV3[llp_tree->nJets] = tmp.dnn_score_v3;

			if (jetMet_dPhiMin_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi)))
			{
				jetMet_dPhiMin_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi));
			}     
			llp_tree->nJets++;
		}
			if(_debug_nj) std::cout << "filled Jets " << std::endl;
		//if(llp_tree->nJets ==0 ) continue;
		if( (label=="MR_JetHT") && llp_tree->nJets != 2 ) continue;
		float jet2_dPhi_temp = 0;
		if(label=="MR_JetHT") 
		{
			//jet2_dPhi_temp = abs(RazorAnalyzerLLP::deltaPhi(llp_tree->jetPhi[0], llp_tree->jetPhi[1]));
			llp_tree->jet2_dPhi = abs(RazorAnalyzerLLP::deltaPhi(llp_tree->jetPhi[0], llp_tree->jetPhi[1]));
		}
		if( (label=="MR_PHO") && llp_tree->nJets != 1 ) continue;
		//if(label=="MR_JetHT" && llp_tree->jet2_dPhi<=2.8) continue;
		//if(label=="MR_JetHT" && jet2_dPhi_temp<=2.8) continue;
		if(_debug) std::cout << "nJets in tree " << llp_tree->nJets << std::endl;
		llp_tree->jetMet_dPhiMin = jetMet_dPhiMin_temp;
		//llp_tree->jetMet_dPhiStarMin = jetMet_dPhiStarMin_temp;
		//llp_tree->jetMet_dPhiMin4 = jetMet_dPhiMin4_temp;


		if(label=="MR_PHO"){
			TLorentzVector photon0 = makeTLorentzVectorPtEtaPhiM(llp_tree->phoPt[0],llp_tree->phoEta[0], llp_tree->phoPhi[0], 0.);
			float jetPho_dPhiMin_temp = 999 ; 

			int MR_PHO_nJets = 0;
			for ( auto &tmp : AK4Jets )
			{

				if (jetPho_dPhiMin_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),llp_tree->phoPhi[0])))
				{
					jetPho_dPhiMin_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),llp_tree->phoPhi[0]));
				}     
				MR_PHO_nJets++;

			}
			if(_debug) std::cout << "nJets in tree " << MR_PHO_nJets << std::endl;
			llp_tree->jetPho_dPhiMin = jetPho_dPhiMin_temp;
			//llp_tree->jetPho_dPhiStarMin = jetPho_dPhiStarMin_temp;
			//llp_tree->jetPho_dPhiMin4 = jetPho_dPhiMin4_temp;

		//if( (label=="MR_PHO") && jetPho_dPhiMin_temp ==999 ) continue;
		//if( (label=="MR_PHO") && jetPho_dPhiMin_temp <=2.8 ) continue;

		}//MR_PHO dPhi


		//llp_tree->tree_->Fill();
		if(_debug) std::cout << "nJets in tree " << llp_tree->nJets << std::endl;

		if(!isData && signalScan)
		{
			pair<int,int> smsPair = make_pair(llp_tree->mX, llp_tree->ctau);
			Trees2D[smsPair]->Fill();
		}
		else
		{
			llp_tree->tree_->Fill();
		}

	}//end fill


	if(!isData && signalScan)
	{
	}
	else
	{
		cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
		cout << "Writing output trees..." << endl;		
		outFile->cd();
		llp_tree->tree_->Write();
		NEvents->Write();
		// outFile->Write();
		outFile->Close();
	}

}
