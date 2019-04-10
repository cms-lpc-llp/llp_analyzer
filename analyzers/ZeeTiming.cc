#include "ZeeTiming.h"
#include "JetCorrectorParameters.h"
#include "RazorHelper.h"

//C++ includes
#include <sys/stat.h>
#include <random>


//ROOT includes
#include "TH1F.h"





using namespace std;
const double SPEED_OF_LIGHT = 29.9792458; // speed of light in cm / ns


const int N_E_divide = 19;
double E_divide[N_E_divide] = {43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0, 67.0, 70.0, 73.0, 78.0, 84.0, 91.0, 100.0, 115.0, 140.0, 190.0, 1000.0};
double timecorr_shift[N_E_divide] = {277.32228254, 275.67500044, 272.10244004, 287.15428262, 294.29875716, 287.44070277, 282.95400642, 279.80934988, 282.70083441, 277.48972614, 275.41563559, 274.64132692, 268.13864834, 266.62392304, 270.24297452, 260.13781686, 263.16768474, 242.21320465, 181.26510246};

double timecorr_smear_aa = 6591.9*6591.9 - 6536.8*6536.8;
double timecorr_smear_bb = 2.0*211.1*211.1 - 2.0*96.2*96.2;

float ZeeTiming::getTimeCalibConstant(TTree *tree, vector <uint> & start_run, vector <uint> & end_run, uint run, uint detID) {
  float timeCalib = 0.0;
  
  int N_entries = tree->GetEntries();
  int i_entry=0;
  for(uint i=0;i<start_run.size();i++)
    {
      if(run>= start_run[i] && run<= end_run[i])
	{
	  i_entry = i;
	  break;
	}
    }
  
  if(i_entry> N_entries) return timeCalib;
  tree->GetEntry(i_entry);
  std::vector<int>::iterator p_id;
  p_id = std::find(detID_all->begin(), detID_all->end(), detID);
  if (p_id == detID_all->end()) return timeCalib;
  uint idx = std::distance(detID_all->begin(), p_id);
  
  if(idx<=IC_time_all->size()) timeCalib = IC_time_all->at(idx);	
  
  return timeCalib;
};

float ZeeTiming::getPedestalNoise(TTree *tree, vector <uint> & start_time, vector <uint> & end_time, uint time, uint detID) {
  float pedestalNoise = 1.0;
  
  int N_entries = tree->GetEntries();
  int i_entry=0;
  for(uint i=0;i<start_time.size();i++)
    {
      if(time>= start_time[i] && time<= end_time[i])
	{
	  i_entry = i;
	  break;
	}
    }
  
  if(i_entry> N_entries) return pedestalNoise;
  tree->GetEntry(i_entry);
  std::vector<int>::iterator p_id;
  p_id = std::find(detID_all->begin(), detID_all->end(), detID);
  if (p_id == detID_all->end()) return pedestalNoise;
  uint idx = std::distance(detID_all->begin(), p_id);
  
  if(idx<=rms_G12_all->size()) pedestalNoise = rms_G12_all->at(idx);	
  
  return pedestalNoise;
};


float ZeeTiming::getADCToGeV( uint run, int isFromEB) {
  double ADCToGeV = 1.0;
  //EB
  if (isFromEB == 1) {
    if (run >= 1 && run <= 271950) ADCToGeV = 0.039680;
    else if (run >= 271951 && run <= 277366) ADCToGeV = 0.039798;
    else if (run >= 277367 && run <= 281825) ADCToGeV = 0.039436;
    else if (run >= 281826 && run <= 999999) ADCToGeV = 0.039298;
  }   
  //EE
  else if (isFromEB == 0) {
    if (run >= 1 && run <= 271950) ADCToGeV = 0.067230;
    else if (run >= 271951 && run <= 277366) ADCToGeV = 0.067370;
    else if (run >= 277367 && run <= 281825) ADCToGeV = 0.066764;
    else if (run >= 281826 && run <= 999999) ADCToGeV = 0.065957;
  }
  return ADCToGeV;
}


void ZeeTiming::Analyze(bool isData, int option, string outFileName, string label)
{

  //*****************************************************************************
  //Settings
  //*****************************************************************************
  TRandom3 random(3003);
  //bool doPhotonScaleCorrection = true;

  string analysisTag = "Razor2016_07Aug2017Rereco";
  if ( label != "") analysisTag = label;

  RazorHelper *helper = 0;
  helper = new RazorHelper(analysisTag, isData, false);


  //*****************************************************************************
  //Load Pedestals
  //*****************************************************************************
  vector <uint> start_time;//start run of all IOV 
  vector <uint> end_time;//end run of all IOV
  start_time_tmp=0; 
  end_time_tmp=0;
  rms_G12_all=0;
  detID_all=0 ;

  TFile *f_pedestal = 0;
  TTree *tree_pedestal = 0;

  if(isData)
  { 
	f_pedestal = TFile::Open("/mnt/hadoop/store/group/phys_susy/razor/Run2Analysis/EcalTiming/EcalPedestals_Legacy2016_time_v1/tree_EcalPedestals_Legacy2016_time_v1_G12rmsonly.root","READ");
	//f_pedestal = TFile::Open("tree_EcalPedestals_Legacy2016_time_v1_G12rmsonly.root","READ");
	tree_pedestal = (TTree*)f_pedestal->Get("pedestal");
	tree_pedestal->SetBranchAddress("start_time_second", &start_time_tmp);
	tree_pedestal->SetBranchAddress("end_time_second", &end_time_tmp);
	tree_pedestal->SetBranchAddress("rms_G12", &rms_G12_all);
	tree_pedestal->SetBranchAddress("detID", &detID_all);
	int N_entries_pedestal = tree_pedestal->GetEntries();

	cout << "Total Pedestal IOVs: " << N_entries_pedestal << "\n";
	for(int i=0;i<N_entries_pedestal;i++) {
		cout << "Loading Pedestal IOV " << i << "\n";
		tree_pedestal->GetEntry(i);
		start_time.push_back(start_time_tmp);
		end_time.push_back(end_time_tmp);
	}
  }
	
  //load timing intercalibration
  //TFile * file_IC_timing_map = new TFile("IC_map_SingleElectron_2016.root","READ");
  ////TFile * file_IC_timing_map = new TFile("IC_average_timing_2016.root","READ");
 // TFile * file_IC_timing_map = new TFile("IC_map_SingleElectron_2016_test_iter19.root","READ");
  //TH2F * h2_IC_timing_map = (TH2F*)file_IC_timing_map->Get("h2_IC_map_iter19");
  ////TH2F * h2_IC_timing_map = (TH2F*)file_IC_timing_map->Get("h2_timeAll_iEta_vs_iPhi_data");

  // //test 
  // uint test_time = 1464000000;
  // //cout<<"EB test..."<<endl;
  // for(int ieta=-85;ieta<=85 && ieta!=0;ieta++) {
  //   for(int iphi=1;iphi<=360;iphi++) {
  //     int detID = detID_from_iEtaiPhi(ieta, iphi, true, false);
  //     cout<<test_time<<"  "<<ieta<<"  "<<iphi<<"  "<<detID;			
  //     float pedestalRMS = getPedestalNoise(tree_pedestal, start_time,end_time,test_time, detID);
  //     cout << "   " << pedestalRMS << endl;			
  //   }
  // }
  
  
  //*****************************************************************************
  //Open Output File
  //*****************************************************************************
  if ( outFileName.empty() ) {
    std::cout << "ZeeTiming: Output filename not specified!" << endl << "Using default output name ZeeTiming.root" << std::endl;
    outFileName = "ZeeTiming.root";
  }
  TFile* outFile = new TFile( outFileName.c_str(), "RECREATE" );




  //---------------------------
  //one tree to hold all events
  //---------------------------
  TTree *outputTree = new TTree("ZeeTiming", "Info on selected razor inclusive events");

  float weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  float mass;
  float t1_TOF2, t2_TOF2;
  float t1, t2;
  float t1_SmearToData, t2_SmearToData;
  float t1_seed_SmearToData, t2_seed_SmearToData;
  float t1_subseed_SmearToData, t2_subseed_SmearToData;
  float t1_ShiftToData, t2_ShiftToData;
  float t1_seed_ShiftToData, t2_seed_ShiftToData;
  float t1_subseed_ShiftToData, t2_subseed_ShiftToData;
  float t1_seed, t2_seed;
  float t1_subseed, t2_subseed;
  float seed1_pedestal, seed2_pedestal;
  float subseed1_pedestal, subseed2_pedestal;
  float seed1_transpCorr, seed2_transpCorr;
  float subseed1_transpCorr, subseed2_transpCorr;
  float t1raw_seed, t2raw_seed;
  float t1calib_seed, t2calib_seed;
  float t1raw_subseed, t2raw_subseed;
  float t1calib_subseed, t2calib_subseed;
  float ele1E, ele1Pt, ele1Eta, ele1Phi, ele1seedE, ele1subseedE;
  float ele2E, ele2Pt, ele2Eta, ele2Phi, ele2seedE, ele2subseedE;
  int  ele1seedIEta, ele1seedIPhi, ele1seedIX, ele1seedIY;
  int  ele1subseedIEta, ele1subseedIPhi, ele1subseedIX, ele1subseedIY;
  bool ele1IsEB;
  int  ele2seedIEta, ele2seedIPhi, ele2seedIX, ele2seedIY; 
  int  ele2subseedIEta, ele2subseedIPhi, ele2subseedIX, ele2subseedIY; 
  bool ele2IsEB;
  int NPU;
  vector<float> *ele1Rechit_E;
  vector<float> *ele1Rechit_rawT;
  vector<float> *ele1Rechit_calibT;
  vector<float> *ele1Rechit_t;
  vector<float> *ele1Rechit_pedestal;
  vector<float> *ele1Rechit_Eta;
  vector<float> *ele1Rechit_IEtaIX;
  vector<float> *ele1Rechit_Phi;
  vector<float> *ele1Rechit_IPhiIY;
  vector<float> *ele1Rechit_transpCorr;

  vector<float> *ele2Rechit_E;
  vector<float> *ele2Rechit_rawT;
  vector<float> *ele2Rechit_calibT;
  vector<float> *ele2Rechit_t;
  vector<float> *ele2Rechit_pedestal;
  vector<float> *ele2Rechit_Eta;
  vector<float> *ele2Rechit_IEtaIX;
  vector<float> *ele2Rechit_Phi;
  vector<float> *ele2Rechit_IPhiIY;
  vector<float> *ele2Rechit_transpCorr;
	
	
  ele1Rechit_E = new std::vector<float>; ele1Rechit_E->clear();
  ele1Rechit_rawT = new std::vector<float>; ele1Rechit_rawT->clear();
  ele1Rechit_calibT = new std::vector<float>; ele1Rechit_calibT->clear();
  ele1Rechit_t = new std::vector<float>; ele1Rechit_t->clear();
  ele1Rechit_pedestal = new std::vector<float>; ele1Rechit_pedestal->clear();
  ele1Rechit_Eta = new std::vector<float>; ele1Rechit_Eta->clear();
  ele1Rechit_IEtaIX = new std::vector<float>; ele1Rechit_IEtaIX->clear();
  ele1Rechit_Phi = new std::vector<float>; ele1Rechit_Phi->clear();
  ele1Rechit_IPhiIY = new std::vector<float>; ele1Rechit_IPhiIY->clear();
  ele1Rechit_transpCorr = new std::vector<float>; ele1Rechit_transpCorr->clear();

  ele2Rechit_E = new std::vector<float>; ele2Rechit_E->clear();
  ele2Rechit_rawT = new std::vector<float>; ele2Rechit_rawT->clear();
  ele2Rechit_calibT = new std::vector<float>; ele2Rechit_calibT->clear();
  ele2Rechit_t = new std::vector<float>; ele2Rechit_t->clear();
  ele2Rechit_pedestal = new std::vector<float>; ele2Rechit_pedestal->clear();
  ele2Rechit_Eta = new std::vector<float>; ele2Rechit_Eta->clear();
  ele2Rechit_IEtaIX = new std::vector<float>; ele2Rechit_IEtaIX->clear();
  ele2Rechit_Phi = new std::vector<float>; ele2Rechit_Phi->clear();
  ele2Rechit_IPhiIY = new std::vector<float>; ele2Rechit_IPhiIY->clear();
  ele2Rechit_transpCorr = new std::vector<float>; ele2Rechit_transpCorr->clear();

  //int nPV;
  unsigned int run, lumi, event;


  outputTree->Branch("weight", &weight, "weight/F");
  outputTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
  outputTree->Branch("pileupWeightUp", &pileupWeightUp, "pileupWeightUp/F");
  outputTree->Branch("pileupWeightDown", &pileupWeightDown, "pileupWeightDown/F");
  outputTree->Branch("run", &run, "run/i");
  outputTree->Branch("lumi", &lumi, "lumi/i");
  outputTree->Branch("event", &event, "event/i");
  outputTree->Branch("eventTime", &eventTime, "eventTime/i");
  outputTree->Branch("NPU", &NPU, "npu/i");
  outputTree->Branch("nPV", &nPV, "nPV/i");
  outputTree->Branch("mass", &mass, "mass/F");
  outputTree->Branch("t1", &t1, "t1/F");
  outputTree->Branch("t2", &t2, "t2/F");
  outputTree->Branch("t1_SmearToData", &t1_SmearToData, "t1_SmearToData/F");
  outputTree->Branch("t1_seed_SmearToData", &t1_seed_SmearToData, "t1_seed_SmearToData/F");
  outputTree->Branch("t1_subseed_SmearToData", &t1_subseed_SmearToData, "t1_subseed_SmearToData/F");
  outputTree->Branch("t1_ShiftToData", &t1_ShiftToData, "t1_ShiftToData/F");
  outputTree->Branch("t1_seed_ShiftToData", &t1_seed_ShiftToData, "t1_seed_ShiftToData/F");
  outputTree->Branch("t1_subseed_ShiftToData", &t1_subseed_ShiftToData, "t1_subseed_ShiftToData/F");
  outputTree->Branch("t2_SmearToData", &t2_SmearToData, "t2_SmearToData/F");
  outputTree->Branch("t2_seed_SmearToData", &t2_seed_SmearToData, "t2_seed_SmearToData/F");
  outputTree->Branch("t2_subseed_SmearToData", &t2_subseed_SmearToData, "t2_subseed_SmearToData/F");
  outputTree->Branch("t2_ShiftToData", &t2_ShiftToData, "t2_ShiftToData/F");
  outputTree->Branch("t2_seed_ShiftToData", &t2_seed_ShiftToData, "t2_seed_ShiftToData/F");
  outputTree->Branch("t2_subseed_ShiftToData", &t2_subseed_ShiftToData, "t2_subseed_ShiftToData/F");
  outputTree->Branch("t1_TOF2", &t1_TOF2, "t1_TOF2/F");
  outputTree->Branch("t2_TOF2", &t2_TOF2, "t2_TOF2/F");
  outputTree->Branch("t1_seed", &t1_seed, "t1_seed/F");
  outputTree->Branch("t1_subseed", &t1_subseed, "t1_subseed/F");
  outputTree->Branch("t2_seed", &t2_seed, "t2_seed/F");
  outputTree->Branch("t2_subseed", &t2_subseed, "t2_subseed/F");
  outputTree->Branch("seed1_pedestal", &seed1_pedestal, "seed1_pedestal/F");
  outputTree->Branch("subseed1_pedestal", &subseed1_pedestal, "subseed1_pedestal/F");
  outputTree->Branch("seed2_pedestal", &seed2_pedestal, "seed2_pedestal/F");
  outputTree->Branch("subseed2_pedestal", &subseed2_pedestal, "subseed2_pedestal/F");
  outputTree->Branch("seed1_transpCorr", &seed1_transpCorr, "seed1_transpCorr/F");
  outputTree->Branch("subseed1_transpCorr", &subseed1_transpCorr, "subseed1_transpCorr/F");
  outputTree->Branch("seed2_transpCorr", &seed2_transpCorr, "seed2_transpCorr/F");
  outputTree->Branch("subseed2_transpCorr", &subseed2_transpCorr, "subseed2_transpCorr/F");
  outputTree->Branch("t1raw_seed", &t1raw_seed, "t1raw_seed/F");
  outputTree->Branch("t1calib_seed", &t1calib_seed, "t1calib_seed/F");
  outputTree->Branch("t1raw_subseed", &t1raw_subseed, "t1raw_subseed/F");
  outputTree->Branch("t1calib_subseed", &t1calib_subseed, "t1calib_subseed/F");
  outputTree->Branch("t2raw_seed", &t2raw_seed, "t2raw_seed/F");
  outputTree->Branch("t2calib_seed", &t2calib_seed, "t2calib_seed/F");
  outputTree->Branch("t2raw_subseed", &t2raw_subseed, "t2raw_subseed/F");
  outputTree->Branch("t2calib_subseed", &t2calib_subseed, "t2calib_subseed/F");
  outputTree->Branch("ele1E", &ele1E, "ele1E/F");
  outputTree->Branch("ele1seedE", &ele1seedE, "ele1seedE/F");
  outputTree->Branch("ele1subseedE", &ele1subseedE, "ele1subseedE/F");
  outputTree->Branch("ele1Pt", &ele1Pt, "ele1Pt/F");
  outputTree->Branch("ele1Eta", &ele1Eta, "ele1Eta/F");
  outputTree->Branch("ele1Phi", &ele1Phi, "ele1Phi/F");
  outputTree->Branch("ele1IsEB", &ele1IsEB, "ele1IsEB/O");
  outputTree->Branch("ele1seedIEta", &ele1seedIEta, "ele1seedIEta/I");
  outputTree->Branch("ele1subseedIEta", &ele1subseedIEta, "ele1subseedIEta/I");
  outputTree->Branch("ele1seedIPhi", &ele1seedIPhi, "ele1seedIPhi/I");
  outputTree->Branch("ele1subseedIPhi", &ele1subseedIPhi, "ele1subseedIPhi/I");
  outputTree->Branch("ele1seedIX", &ele1seedIX, "ele1seedIX/I");
  outputTree->Branch("ele1subseedIX", &ele1subseedIX, "ele1subseedIX/I");
  outputTree->Branch("ele1seedIY", &ele1seedIY, "ele1seedIY/I");
  outputTree->Branch("ele1subseedIY", &ele1subseedIY, "ele1subseedIY/I");
  outputTree->Branch("ele2E", &ele2E, "ele2E/F");
  outputTree->Branch("ele2seedE", &ele2seedE, "ele2seedE/F");
  outputTree->Branch("ele2subseedE", &ele2subseedE, "ele2subseedE/F");
  outputTree->Branch("ele2Pt", &ele2Pt, "ele2Pt/F");
  outputTree->Branch("ele2Eta", &ele2Eta, "ele2Eta/F");
  outputTree->Branch("ele2Phi", &ele2Phi, "ele2Phi/F");
  outputTree->Branch("ele2IsEB", &ele2IsEB, "ele2IsEB/O");
  outputTree->Branch("ele2seedIEta", &ele2seedIEta, "ele2seedIEta/I");
  outputTree->Branch("ele2subseedIEta", &ele2subseedIEta, "ele2subseedIEta/I");
  outputTree->Branch("ele2seedIPhi", &ele2seedIPhi, "ele2seedIPhi/I");
  outputTree->Branch("ele2subseedIPhi", &ele2subseedIPhi, "ele2subseedIPhi/I");
  outputTree->Branch("ele2seedIX", &ele2seedIX, "ele2seedIX/I");
  outputTree->Branch("ele2subseedIX", &ele2subseedIX, "ele2subseedIX/I");
  outputTree->Branch("ele2seedIY", &ele2seedIY, "ele2seedIY/I");
  outputTree->Branch("ele2subseedIY", &ele2subseedIY, "ele2subseedIY/I");
  outputTree->Branch("ele1Rechit_E", "std::vector<float>",&ele1Rechit_E);
  outputTree->Branch("ele1Rechit_rawT", "std::vector<float>",&ele1Rechit_rawT);
  outputTree->Branch("ele1Rechit_calibT", "std::vector<float>",&ele1Rechit_calibT);
  outputTree->Branch("ele1Rechit_t", "std::vector<float>",&ele1Rechit_t);
  outputTree->Branch("ele1Rechit_pedestal", "std::vector<float>",&ele1Rechit_pedestal);
  outputTree->Branch("ele1Rechit_Eta", "std::vector<float>",&ele1Rechit_Eta);
  outputTree->Branch("ele1Rechit_IEtaIX", "std::vector<float>",&ele1Rechit_IEtaIX);
  outputTree->Branch("ele1Rechit_Phi", "std::vector<float>",&ele1Rechit_Phi);
  outputTree->Branch("ele1Rechit_IPhiIY", "std::vector<float>",&ele1Rechit_IPhiIY);
  outputTree->Branch("ele1Rechit_transpCorr", "std::vector<float>",&ele1Rechit_transpCorr);

  outputTree->Branch("ele2Rechit_E", "std::vector<float>",&ele2Rechit_E);
  outputTree->Branch("ele2Rechit_rawT", "std::vector<float>",&ele2Rechit_rawT);
  outputTree->Branch("ele2Rechit_calibT", "std::vector<float>",&ele2Rechit_calibT);
  outputTree->Branch("ele2Rechit_t", "std::vector<float>",&ele2Rechit_t);
  outputTree->Branch("ele2Rechit_pedestal", "std::vector<float>",&ele2Rechit_pedestal);
  outputTree->Branch("ele2Rechit_Eta", "std::vector<float>",&ele2Rechit_Eta);
  outputTree->Branch("ele2Rechit_IEtaIX", "std::vector<float>",&ele2Rechit_IEtaIX);
  outputTree->Branch("ele2Rechit_Phi", "std::vector<float>",&ele2Rechit_Phi);
  outputTree->Branch("ele2Rechit_IPhiIY", "std::vector<float>",&ele2Rechit_IPhiIY);
  outputTree->Branch("ele2Rechit_transpCorr", "std::vector<float>",&ele2Rechit_transpCorr);


  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);


  //begin loop
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //begin event
    if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    //initialize branches
    weight = 0;
    pileupWeight = 0;
    pileupWeightUp = 0;
    pileupWeightDown = 0;
    mass = 0;
    t1 = -999;
    t2 = -999;
    t1_SmearToData = -999;
    t1_seed_SmearToData = -999;
    t1_subseed_SmearToData = -999;
    t1_ShiftToData = -999;
    t1_seed_ShiftToData = -999;
    t1_subseed_ShiftToData = -999;
    t2_SmearToData = -999;
    t2_seed_SmearToData = -999;
    t2_subseed_SmearToData = -999;
    t2_ShiftToData = -999;
    t2_seed_ShiftToData = -999;
    t2_subseed_ShiftToData = -999;
    t1_TOF2 = -999;
    t2_TOF2 = -999;
    t1_seed = -999;
    t1_subseed = -999;
    t2_seed = -999;
    t2_subseed = -999;
    seed1_pedestal = -999;
    subseed1_pedestal = -999;
    seed2_pedestal = -999;
    subseed2_pedestal = -999;
    seed1_transpCorr = -999;
    subseed1_transpCorr = -999;
    seed2_transpCorr = -999;
    subseed2_transpCorr = -999;
    t1raw_seed = -999;
    t1calib_seed = -999;
    t1raw_subseed = -999;
    t1calib_subseed = -999;
    t2raw_seed = -999;
    t2calib_seed = -999;
    t2raw_subseed = -999;
    t2calib_subseed = -999;
    ele1E = -999; 
    ele1seedE = -999; 
    ele1subseedE = -999; 
    ele1Pt = -999;
    ele1Eta = -999;
    ele1Phi = -999;
    ele2E = -999;
    ele2seedE = -999; 
    ele2subseedE = -999; 
    ele2Pt = -999;
    ele2Eta = -999;
    ele2Phi = -999;
    ele1seedIEta = -999;
    ele1subseedIEta = -999;
    ele1seedIPhi = -999;
    ele1subseedIPhi = -999;
    ele1seedIX = -999;
    ele1subseedIX = -999;
    ele1seedIY = -999;
    ele1subseedIY = -999;
    ele1IsEB = 0;
    ele2seedIEta = -999;
    ele2subseedIEta = -999;
    ele2seedIPhi = -999;
    ele2subseedIPhi = -999;
    ele2seedIX = -999;
    ele2subseedIX = -999;
    ele2seedIY = -999; 
    ele2subseedIY = -999; 
    ele2IsEB = 0;
    NPU = 0;   
    ele1Rechit_E->clear();
    ele1Rechit_rawT->clear();
    ele1Rechit_calibT->clear();
    ele1Rechit_t->clear();
    ele1Rechit_pedestal->clear();
    ele1Rechit_Eta->clear();
    ele1Rechit_IEtaIX->clear();
    ele1Rechit_Phi->clear();
    ele1Rechit_IPhiIY->clear();
    ele1Rechit_transpCorr->clear();

    ele2Rechit_E->clear();
    ele2Rechit_rawT->clear();
    ele2Rechit_calibT->clear();
    ele2Rechit_t->clear();
    ele2Rechit_pedestal->clear();
    ele2Rechit_Eta->clear();
    ele2Rechit_IEtaIX->clear();
    ele2Rechit_Phi->clear();
    ele2Rechit_IPhiIY->clear();
    ele2Rechit_transpCorr->clear();


    //fill normalization histogram
    NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
    weight = genWeight;


    //get NPU
     if( !isData )
    {
    for (int i=0; i < nBunchXing; ++i) {
      if (BunchXing[i] == 0) {
        NPU = nPUmean[i];
      }
    }
    pileupWeight = helper->getPileupWeight(NPU);
    pileupWeightUp = helper->getPileupWeightUp(NPU) / pileupWeight;
    pileupWeightDown = helper->getPileupWeightDown(NPU) / pileupWeight;


    }
   
    run = runNum;
    lumi = lumiNum;
    event = eventNum;
    
    double pvX = 0;


    int nEle = 0;
    TLorentzVector ele1 = makeTLorentzVector(0,0,0,0);
    TLorentzVector ele2 = makeTLorentzVector(0,0,0,0);
    double ele1_time = 0;
    double ele2_time= 0;
    double ele1_time_TOF2 = 0;
    double ele2_time_TOF2 = 0;
    double ele1_seedtime = 0;
    double ele1_subseedtime = 0;
    double ele2_seedtime = 0;
    double ele2_subseedtime = 0;
    double ele1_seedtimeraw = 0;
    double ele1_seedtimecalib = 0;
    double ele1_subseedtimeraw = 0;
    double ele1_subseedtimecalib = 0;
    double ele2_seedtimeraw = 0;
    double ele2_seedtimecalib = 0;
    double ele2_subseedtimeraw = 0;
    double ele2_subseedtimecalib = 0;

    for(int i = 0; i < nElectrons; i++){
      //if(elePt[i] < 35) continue;
      if(fabs(eleEta[i]) > 2.5) continue;
      if(fabs(eleEta[i]) > 1.4442 && fabs(eleEta[i]) < 1.566) continue;
      if(!(isEGammaPOGTightElectron(i))) continue;
      if(ecalRechit_ID->empty() ) continue;
      
      TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
          
      //rough definition
      uint seedhitIndex =  (*ele_SeedRechitIndex)[i];
      bool isFromEB = bool( (*ecalRechit_ID)[seedhitIndex] < 840000000 );
      double rawseedHitTime =  (*ecalRechit_T)[seedhitIndex];
      double eleseedE =  (*ecalRechit_E)[seedhitIndex];
      if(eleseedE < 10.0) continue;

      //find the subseed: the highest energy crystal among the neighboring crystals of the seed 
      double IC_timing_seed_zhicai = 0.0;

      int seed_iEtaiX = 0; 
      int seed_iPhiiY = 0; 
      if(isFromEB)
      {
		seed_iEtaiX = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
		seed_iPhiiY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
		//if(isData) IC_timing_seed_zhicai = -1.0* h2_IC_timing_map->GetBinContent(seed_iEtaiX+85+1,seed_iPhiiY);
      }
      else
      {
		seed_iEtaiX = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
                seed_iPhiiY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
      }

      double elesubseedE = -999.0;
      uint subseedhitIndex = 0;
      double IC_timing_subseed_zhicai = 0.0;

      for (uint k=0; k<(*ele_EcalRechitIndex)[i].size(); ++k) {
        uint rechitIndex = (*ele_EcalRechitIndex)[i][k];
	if (rechitIndex == seedhitIndex) continue;

	double thisRechitE =(*ecalRechit_E)[rechitIndex];
	if (thisRechitE < elesubseedE) continue;

	bool isThisfromEB = bool( (*ecalRechit_ID)[rechitIndex] < 840000000 );
	int this_iEtaiX = 0;
	int this_iPhiiY = 0;
        if(isThisfromEB)
	{
		this_iEtaiX = iEta_or_iX_from_detID( (*ecalRechit_ID)[rechitIndex] , true);
                this_iPhiiY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[rechitIndex] , true);
	}
	else
	{
		this_iEtaiX = iEta_or_iX_from_detID( (*ecalRechit_ID)[rechitIndex] , false);
                this_iPhiiY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[rechitIndex] , false);	
	}
	int distance_from_seed = (this_iEtaiX-seed_iEtaiX)*(this_iEtaiX-seed_iEtaiX) + (this_iPhiiY-seed_iPhiiY)*(this_iPhiiY-seed_iPhiiY);
	if(distance_from_seed >  1) continue;//only accept neighboring crystals
	
	if(thisRechitE > elesubseedE)
	{
		elesubseedE = thisRechitE;
		subseedhitIndex = rechitIndex;
	}

      }


      if(isFromEB && isData)
      {
		int subseed_iEtaiX = iEta_or_iX_from_detID( (*ecalRechit_ID)[subseedhitIndex] , true);
		int subseed_iPhiiY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[subseedhitIndex] , true);
		//IC_timing_subseed_zhicai = -1.0* h2_IC_timing_map->GetBinContent(subseed_iEtaiX+85+1,subseed_iPhiiY);
      }
	
      double rawsubseedHitTime =  (*ecalRechit_T)[subseedhitIndex];
	
      double calibseedHitTime =  rawseedHitTime + IC_timing_seed_zhicai;
      double calibsubseedHitTime =  rawsubseedHitTime + IC_timing_subseed_zhicai;
      //apply TOF correction
      double TOFCorrectedseedHitTime = rawseedHitTime + IC_timing_seed_zhicai + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[seedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[seedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT;
      double TOFCorrectedsubseedHitTime = rawsubseedHitTime + IC_timing_subseed_zhicai + (std::sqrt(pow((*ecalRechit_X)[subseedhitIndex],2)+pow((*ecalRechit_Y)[subseedhitIndex],2)+pow((*ecalRechit_Z)[subseedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[subseedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[subseedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[subseedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT;


      // cout << "Ele: " << i << " : " << elePt[i] << " " << eleEta[i] << " : \n" 
      //	<< "  runNum: " << runNum << "  detID: " << (*ecalRechit_ID)[seedhitIndex] << "  IC_time_SeptRereco: " << IC_time_SeptRereco << "  IC_time_LagacyRereco: " << IC_time_LagacyRereco << " \n"
      //   	<< " " << rawseedHitTime << " -> " << calibratedseedHitTime << " -> " << TOFCorrectedseedHitTime << " "
      //	<< "\n";
 


      double seedPedNoise = isData ? getPedestalNoise(tree_pedestal, start_time,end_time, eventTime, (*ecalRechit_ID)[seedhitIndex]) : 1.0;
      double subseedPedNoise = isData ? getPedestalNoise(tree_pedestal, start_time,end_time, eventTime, (*ecalRechit_ID)[subseedhitIndex]) : 1.0;
   
      double tmpSumWeightedTime = 0;
      double tmpSumWeightedTime_TOF2 = 0;
      double tmpSumWeight = 0;
      double tmpSumWeight_TOF2 = 0;

      vector<float> *eleRechit_E;
      vector<float> *eleRechit_rawT;
      vector<float> *eleRechit_calibT;
      vector<float> *eleRechit_t;
      vector<float> *eleRechit_pedestal;
      vector<float> *eleRechit_Eta;
      vector<float> *eleRechit_IEtaIX;
      vector<float> *eleRechit_Phi;
      vector<float> *eleRechit_IPhiIY;
      vector<float> *eleRechit_transpCorr;

      eleRechit_E = new std::vector<float>;     eleRechit_E->clear();
      eleRechit_rawT = new std::vector<float>;     eleRechit_rawT->clear();
      eleRechit_calibT = new std::vector<float>;     eleRechit_calibT->clear();
      eleRechit_t = new std::vector<float>;     eleRechit_t->clear();
      eleRechit_pedestal = new std::vector<float>;     eleRechit_pedestal->clear();
      eleRechit_Eta = new std::vector<float>;     eleRechit_Eta->clear();
      eleRechit_IEtaIX = new std::vector<float>;     eleRechit_IEtaIX->clear();
      eleRechit_Phi = new std::vector<float>;     eleRechit_Phi->clear();
      eleRechit_IPhiIY = new std::vector<float>;     eleRechit_IPhiIY->clear();
      eleRechit_transpCorr = new std::vector<float>;     eleRechit_transpCorr->clear();

      for (uint k=0; k<(*ele_EcalRechitIndex)[i].size(); ++k) {	
	uint rechitIndex = (*ele_EcalRechitIndex)[i][k];

   	if((*ecalRechit_E)[rechitIndex] < 1.0) continue;
 
	bool isEBrechit = bool( abs((*ecalRechit_Eta)[rechitIndex]) < 1.5 );

	double rawT_this = (*ecalRechit_T)[rechitIndex];
	double IC_timing_this_zhicai = 0.0;	

	if(isEBrechit && isData){
		int rec_iEtaiX = iEta_or_iX_from_detID( (*ecalRechit_ID)[rechitIndex] , true);
		int rec_iPhiiY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[rechitIndex] , true);
		//IC_timing_this_zhicai = -1.0* h2_IC_timing_map->GetBinContent(rec_iEtaiX+85+1,rec_iPhiiY);
	}

	double calibT_this = rawT_this + IC_timing_seed_zhicai;	
	//apply TOF correction
	double corrT = rawT_this + IC_timing_seed_zhicai + (std::sqrt(pow((*ecalRechit_X)[rechitIndex],2)+pow((*ecalRechit_Y)[rechitIndex],2)+pow((*ecalRechit_Z)[rechitIndex],2))-std::sqrt(pow((*ecalRechit_X)[rechitIndex]-pvX,2)+pow((*ecalRechit_Y)[rechitIndex]-pvY,2)+pow((*ecalRechit_Z)[rechitIndex]-pvZ,2)))/SPEED_OF_LIGHT;

  	eleRechit_E->push_back((*ecalRechit_E)[rechitIndex]); 
        eleRechit_Eta->push_back((*ecalRechit_Eta)[rechitIndex]); 
        eleRechit_Phi->push_back((*ecalRechit_Phi)[rechitIndex]); 
        if(isData) eleRechit_transpCorr->push_back((*ecalRechit_transpCorr)[rechitIndex]); 
        else eleRechit_transpCorr->push_back(1.0); 

        eleRechit_rawT->push_back(rawT_this); 
        eleRechit_calibT->push_back(calibT_this); 
        eleRechit_t->push_back(corrT); 

        if (isEBrechit) {
          eleRechit_IEtaIX->push_back(iEta_or_iX_from_detID( (*ecalRechit_ID)[rechitIndex] , true));
          eleRechit_IPhiIY->push_back(iPhi_or_iY_from_detID( (*ecalRechit_ID)[rechitIndex] , true));
        } else {
          eleRechit_IEtaIX->push_back(iEta_or_iX_from_detID( (*ecalRechit_ID)[rechitIndex] , false));
          eleRechit_IPhiIY->push_back(iPhi_or_iY_from_detID( (*ecalRechit_ID)[rechitIndex] , false));
        }	

	double pedNoise = isData ? getPedestalNoise(tree_pedestal, start_time,end_time, eventTime, (*ecalRechit_ID)[rechitIndex]) : 0.042; // 42 MeV for MC
	//double pedNoise = 1;
	double ADCToGeV = isData ? getADCToGeV(runNum, isFromEB) : 1.0;
	double sigmaE = pedNoise * ADCToGeV;
        eleRechit_pedestal->push_back(sigmaE); 

	double sigmaT2 = 0.5*N_EB*N_EB / ((*ecalRechit_E)[rechitIndex] * (*ecalRechit_E)[rechitIndex] / (sigmaE*sigmaE)) + C_EB * C_EB;
	
	if(!isData) sigmaT2 = 0.5*N_EB_MC*N_EB_MC / ((*ecalRechit_E)[rechitIndex] * (*ecalRechit_E)[rechitIndex] / (sigmaE*sigmaE)) + C_EB_MC * C_EB_MC;
 	//if(!isData) sigmaT2 = 1.0 / ((*ecalRechit_E)[rechitIndex] * (*ecalRechit_E)[rechitIndex] / (sigmaE*sigmaE));

	tmpSumWeightedTime += corrT * ( 1.0 / sigmaT2 );
	tmpSumWeightedTime_TOF2 += rawT_this * ( 1.0 / sigmaT2 );
	tmpSumWeight += ( 1.0 / sigmaT2 );
	tmpSumWeight_TOF2 += ( 1.0 / sigmaT2 );
	// cout << "\n";
      }
      double weightedTime = tmpSumWeightedTime / tmpSumWeight;
      double weightedTime_TOF2 = tmpSumWeightedTime_TOF2 / tmpSumWeight_TOF2 + (std::sqrt(pow((*ecalRechit_X)[seedhitIndex],2)+pow((*ecalRechit_Y)[seedhitIndex],2)+pow((*ecalRechit_Z)[seedhitIndex],2))-std::sqrt(pow((*ecalRechit_X)[seedhitIndex]-pvX,2)+pow((*ecalRechit_Y)[seedhitIndex]-pvY,2)+pow((*ecalRechit_Z)[seedhitIndex]-pvZ,2)))/SPEED_OF_LIGHT;

            
      if (thisElectron.Pt() > ele1.Pt()) {
	///give ele1 to ele2
	ele2Rechit_E->clear();
        ele2Rechit_rawT->clear();
        ele2Rechit_calibT->clear();
        ele2Rechit_t->clear();
        ele2Rechit_pedestal->clear();
        ele2Rechit_Eta->clear();
        ele2Rechit_IEtaIX->clear();
        ele2Rechit_Phi->clear();
        ele2Rechit_IPhiIY->clear();
        ele2Rechit_transpCorr->clear();

	for(int i=0;i<ele1Rechit_E->size();i++)
	{
		ele2Rechit_E->push_back(ele1Rechit_E->at(i));
		ele2Rechit_rawT->push_back(ele1Rechit_rawT->at(i));
		ele2Rechit_calibT->push_back(ele1Rechit_calibT->at(i));
		ele2Rechit_t->push_back(ele1Rechit_t->at(i));
		ele2Rechit_pedestal->push_back(ele1Rechit_pedestal->at(i));
		ele2Rechit_Eta->push_back(ele1Rechit_Eta->at(i));
		ele2Rechit_IEtaIX->push_back(ele1Rechit_IEtaIX->at(i));
		ele2Rechit_Phi->push_back(ele1Rechit_Phi->at(i));
		ele2Rechit_IPhiIY->push_back(ele1Rechit_IPhiIY->at(i));
		ele2Rechit_transpCorr->push_back(ele1Rechit_transpCorr->at(i));
	}

	ele2=ele1;
        ele2E=ele1E, ele2Pt=ele1Pt, ele2Eta=ele1Eta, ele2Phi=ele1Phi, ele2seedE=ele1seedE, ele2subseedE=ele1subseedE;
        ele2seedIEta=ele1seedIEta, ele2seedIPhi=ele1seedIPhi, ele2seedIX=ele1seedIX, ele2seedIY=ele1seedIY;
        ele2subseedIEta=ele1subseedIEta, ele2subseedIPhi=ele1subseedIPhi, ele2subseedIX=ele1subseedIX, ele2subseedIY=ele1subseedIY;
        ele2IsEB=ele1IsEB;
	seed2_pedestal=seed1_pedestal;
	subseed2_pedestal=subseed1_pedestal;
	ele2_time = ele1_time;
	ele2_time_TOF2 = ele1_time_TOF2;
	ele2_seedtime = ele1_seedtime;
	ele2_subseedtime = ele1_subseedtime;
	ele2_seedtimeraw = ele1_seedtimeraw;
	ele2_seedtimecalib = ele1_seedtimecalib;
	ele2_subseedtimeraw = ele1_subseedtimeraw;
	ele2_subseedtimecalib = ele1_subseedtimecalib;

	///give ele to ele1
	//rechits..
        ele1Rechit_E->clear();
        ele1Rechit_rawT->clear();
        ele1Rechit_calibT->clear();
        ele1Rechit_t->clear();
        ele1Rechit_pedestal->clear();
        ele1Rechit_Eta->clear();
        ele1Rechit_IEtaIX->clear();
        ele1Rechit_Phi->clear();
        ele1Rechit_IPhiIY->clear();
        ele1Rechit_transpCorr->clear();

	for(int i=0;i<eleRechit_E->size();i++)
	{
		ele1Rechit_E->push_back(eleRechit_E->at(i));
		ele1Rechit_rawT->push_back(eleRechit_rawT->at(i));
		ele1Rechit_calibT->push_back(eleRechit_calibT->at(i));
		ele1Rechit_t->push_back(eleRechit_t->at(i));
		ele1Rechit_pedestal->push_back(eleRechit_pedestal->at(i));
		ele1Rechit_Eta->push_back(eleRechit_Eta->at(i));
		ele1Rechit_IEtaIX->push_back(eleRechit_IEtaIX->at(i));
		ele1Rechit_Phi->push_back(eleRechit_Phi->at(i));
		ele1Rechit_IPhiIY->push_back(eleRechit_IPhiIY->at(i));
		ele1Rechit_transpCorr->push_back(eleRechit_transpCorr->at(i));
	}
	
	/////
	ele1 = thisElectron;
	seed1_pedestal = seedPedNoise; 
	subseed1_pedestal = subseedPedNoise; 
	if(isData) seed1_transpCorr = (*ecalRechit_transpCorr)[seedhitIndex]; 
	if(isData) subseed1_transpCorr = (*ecalRechit_transpCorr)[subseedhitIndex]; 
	ele1_time = weightedTime;
	ele1_time_TOF2 = weightedTime_TOF2;
	ele1_seedtime = TOFCorrectedseedHitTime;
	ele1_subseedtime = TOFCorrectedsubseedHitTime;
	ele1_seedtimeraw = rawseedHitTime;
	ele1_seedtimecalib = calibseedHitTime;
	ele1_subseedtimeraw = rawsubseedHitTime;
	ele1_subseedtimecalib = calibsubseedHitTime;
	ele1seedE = eleseedE;
	ele1subseedE = elesubseedE;
	ele1IsEB = bool( abs(eleEta_SC[i]) < 1.5 );
	if (ele1IsEB) {
	  ele1seedIEta = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
	  ele1subseedIEta = iEta_or_iX_from_detID( (*ecalRechit_ID)[subseedhitIndex] , true);
	  ele1seedIPhi = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
	  ele1subseedIPhi = iPhi_or_iY_from_detID( (*ecalRechit_ID)[subseedhitIndex] , true);
	  ele1seedIX = -999;
	  ele1subseedIX = -999;
	  ele1seedIY = -999;
	  ele1subseedIY = -999;
	} else {
	  ele1seedIEta = -999;
	  ele1subseedIEta = -999;
	  ele1seedIPhi = -999;
	  ele1subseedIPhi = -999;
	  ele1seedIX = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
	  ele1subseedIX = iEta_or_iX_from_detID( (*ecalRechit_ID)[subseedhitIndex] , false);
	  ele1seedIY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
	  ele1subseedIY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[subseedhitIndex] , false);
	}
	
      } else if (thisElectron.Pt() > ele2.Pt()) {
	///give ele to ele2
	//rechits..
        ele2Rechit_E->clear();
        ele2Rechit_rawT->clear();
        ele2Rechit_calibT->clear();
        ele2Rechit_t->clear();
        ele2Rechit_pedestal->clear();
        ele2Rechit_Eta->clear();
        ele2Rechit_IEtaIX->clear();
        ele2Rechit_Phi->clear();
        ele2Rechit_IPhiIY->clear();
        ele2Rechit_transpCorr->clear();

	for(int i=0;i<eleRechit_E->size();i++)
	{
		ele2Rechit_E->push_back(eleRechit_E->at(i));
		ele2Rechit_rawT->push_back(eleRechit_rawT->at(i));
		ele2Rechit_calibT->push_back(eleRechit_calibT->at(i));
		ele2Rechit_t->push_back(eleRechit_t->at(i));
		ele2Rechit_pedestal->push_back(eleRechit_pedestal->at(i));
		ele2Rechit_Eta->push_back(eleRechit_Eta->at(i));
		ele2Rechit_IEtaIX->push_back(eleRechit_IEtaIX->at(i));
		ele2Rechit_Phi->push_back(eleRechit_Phi->at(i));
		ele2Rechit_IPhiIY->push_back(eleRechit_IPhiIY->at(i));
		ele2Rechit_transpCorr->push_back(eleRechit_transpCorr->at(i));
	}
	
	ele2 = thisElectron;
	seed2_pedestal = seedPedNoise; 
	subseed2_pedestal = subseedPedNoise; 
	if(isData) seed2_transpCorr = (*ecalRechit_transpCorr)[seedhitIndex]; 
	if(isData) subseed2_transpCorr = (*ecalRechit_transpCorr)[subseedhitIndex]; 
	ele2_time = weightedTime;
	ele2_time_TOF2 = weightedTime_TOF2;
	ele2_seedtime = TOFCorrectedseedHitTime; 
	ele2_subseedtime = TOFCorrectedsubseedHitTime; 
 	ele2_seedtimeraw = rawseedHitTime;
 	ele2_seedtimecalib = calibseedHitTime;
 	ele2_subseedtimeraw = rawsubseedHitTime;
 	ele2_subseedtimecalib = calibsubseedHitTime;
	ele2seedE = eleseedE;
	ele2subseedE = elesubseedE;
 	ele2IsEB = bool( abs(eleEta_SC[i]) < 1.5 );
	if (ele2IsEB) {
	  ele2seedIEta = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
	  ele2subseedIEta = iEta_or_iX_from_detID( (*ecalRechit_ID)[subseedhitIndex] , true);
	  ele2seedIPhi = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , true);
	  ele2subseedIPhi = iPhi_or_iY_from_detID( (*ecalRechit_ID)[subseedhitIndex] , true);
	  ele2seedIX = -999;
	  ele2subseedIX = -999;
	  ele2seedIY = -999;
	  ele2subseedIY = -999;
	} else {
	  ele2seedIEta = -999;
	  ele2subseedIEta = -999;
	  ele2seedIPhi = -999;
	  ele2subseedIPhi = -999;
	  ele2seedIX = iEta_or_iX_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
	  ele2subseedIX = iEta_or_iX_from_detID( (*ecalRechit_ID)[subseedhitIndex] , false);
	  ele2seedIY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[seedhitIndex] , false);
	  ele2subseedIY = iPhi_or_iY_from_detID( (*ecalRechit_ID)[subseedhitIndex] , false);
	}

      }	
   
	delete eleRechit_E;
        delete eleRechit_rawT;
        delete eleRechit_calibT;
        delete eleRechit_t;
        delete eleRechit_pedestal;
        delete eleRechit_Eta;
        delete eleRechit_IEtaIX;
        delete eleRechit_Phi;
        delete eleRechit_IPhiIY;
        delete eleRechit_transpCorr;
      	
	nEle++;

    }

 
    if (nEle >= 2) {
      ele1E = ele1.E();
      ele1Pt = ele1.Pt();
      ele1Eta = ele1.Eta();
      ele1Phi = ele1.Phi();
      ele2E = ele2.E();
      ele2Pt = ele2.Pt();
      ele2Eta = ele2.Eta();
      ele2Phi = ele2.Phi();

      mass = (ele1+ele2).M();
      t1 = ele1_time;
      t1_TOF2 = ele1_time_TOF2;
      t2 = ele2_time;
      t2_TOF2 = ele2_time_TOF2;
      t1_seed = ele1_seedtime;
      t1_subseed = ele1_subseedtime;
      t2_seed = ele2_seedtime;
      t2_subseed = ele2_subseedtime;
      t1raw_seed = ele1_seedtimeraw;
      t1calib_seed = ele1_seedtimecalib;
      t1raw_subseed = ele1_subseedtimeraw;
      t1calib_subseed = ele1_subseedtimecalib;
      t2raw_seed = ele2_seedtimeraw;
      t2calib_seed = ele2_seedtimecalib;
      t2raw_subseed = ele2_subseedtimeraw;
      t2calib_subseed = ele2_subseedtimecalib;
       //cout << "ele2: " << ele2.Pt() << " " << ele2_seedtime << "\n";
       //cout << "ele2: " << ele2.Pt() << " " << ele2_subseedtime << "\n";

if(!isData)
{
        double TR_SMEAR1 = 0.0;
        double TR_SMEAR2 = 0.0;
        double TR_SHIFT1 = 0.0;
        double TR_SHIFT2 = 0.0;
        int E_bin1 = 0;
        int E_bin2 = 0;
        for(int ipt = 0; ipt <N_E_divide; ipt++)
        {
                if(ele1E>E_divide[ipt]) E_bin1 ++;
                if(ele2E>E_divide[ipt]) E_bin2 ++;
        }

        if(E_bin1 >= N_E_divide) E_bin1 = N_E_divide-1;
        if(E_bin2 >= N_E_divide) E_bin2 = N_E_divide-1;

        TR_SHIFT1 = 0.001*timecorr_shift[E_bin1];
        TR_SHIFT2 = 0.001*timecorr_shift[E_bin2];

        if(ele1Pt>0.0) TR_SMEAR1 = 0.001*sqrt((timecorr_smear_aa/(ele1Pt*ele1Pt) + timecorr_smear_bb)/2.0);
        if(ele2Pt>0.0) TR_SMEAR2 = 0.001*sqrt((timecorr_smear_aa/(ele2Pt*ele2Pt) + timecorr_smear_bb)/2.0);

        std::random_device rd;
        std::mt19937 e2(rd());
        std::normal_distribution<> dist1(t1, TR_SMEAR1);
        std::normal_distribution<> dist1_seed(t1_seed, TR_SMEAR1);
        std::normal_distribution<> dist1_subseed(t1_subseed, TR_SMEAR1);
        std::normal_distribution<> dist2(t2, TR_SMEAR2);
        std::normal_distribution<> dist2_seed(t2_seed, TR_SMEAR2);
        std::normal_distribution<> dist2_subseed(t2_subseed, TR_SMEAR2);
        t1_SmearToData = dist1(e2) + TR_SHIFT1;
        t1_seed_SmearToData = dist1_seed(e2) + TR_SHIFT1;
        t1_subseed_SmearToData = dist1_subseed(e2) + TR_SHIFT1;
        t1_ShiftToData = t1 + TR_SHIFT1;
        t1_seed_ShiftToData = t1_seed + TR_SHIFT1;
        t1_subseed_ShiftToData = t1_subseed + TR_SHIFT1;
        t2_SmearToData = dist2(e2) + TR_SHIFT2;
        t2_seed_SmearToData = dist2_seed(e2) + TR_SHIFT2;
        t2_subseed_SmearToData = dist2_subseed(e2) + TR_SHIFT2;
        t2_ShiftToData = t2 + TR_SHIFT2;
        t2_seed_ShiftToData = t2_seed + TR_SHIFT2;
        t2_subseed_ShiftToData = t2_subseed + TR_SHIFT2;
}
else
{
        t1_SmearToData = t1;
        t1_ShiftToData = t1;
        t2_SmearToData = t2;
        t2_ShiftToData = t2;

  	t1_seed_SmearToData = t1_seed;
  	t1_subseed_SmearToData = t1_subseed;
        t1_seed_ShiftToData = t1_seed;
        t1_subseed_ShiftToData = t1_subseed;
        t2_seed_SmearToData = t2_seed;
        t2_subseed_SmearToData = t2_subseed;
        t2_seed_ShiftToData = t2_seed;
        t2_subseed_ShiftToData = t2_subseed;

}
	//Fill Event
    if (mass > 60 && mass < 150) {
      outputTree->Fill();
    }

    }

    

  }//end of event loop
  
  cout << "Writing output trees..." << endl;
  outputTree->Write();
  NEvents->Write();

  outFile->Close();
}
