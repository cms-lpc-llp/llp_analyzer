#include "RazorInclusive.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>
#include <assert.h>
#include <sstream>

//ROOT includes
#include "TH1F.h"
#include "TH2D.h"

using namespace std;

void RazorInclusive::Analyze(bool isData, int option, string outFileName, string label)
{
  //initialization: create one TTree for each analysis box 
  bool combineTrees = true;
  bool isFastsimSMS = (option == 1);
  cout << "Initializing..." << endl;
  if(isFastsimSMS) cout << "Will split up the output according to SUSY particle masses" << endl;
  bool printdebug = false;

  //Pileup Weights
  // TFile *pileupWeightFile = 0;
  // TH1D *pileupWeightHist = 0;
  // No pileup reweighting procedure for Run2 yet
  // pileupWeightFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run1Analysis/ScaleFactors/Run1/Run1PileupWeights.root");
  // pileupWeightHist = (TH1D*)pileupWeightFile->Get("PUWeight_Run1");
  // assert(pileupWeightHist);

  //Lepton Efficiency Correction Factors
  // TH2D *eleLooseEffSFHist = 0;
  // TH2D *eleTightEffSFHist = 0;
  // No Scale factors for Run2 yet
  // TFile *eleEffSFFile = TFile::Open("root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/Run2/ElectronSelection.root");
  // eleLooseEffSFHist = (TH2D*)eleEffSFFile->Get("sfLOOSE");
  // assert(eleLooseEffSFHist);
  // eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("sfTIGHT");
  // assert(eleTightEffSFHist);
  
  if (outFileName.empty()){
    cout << "RazorInclusive: Output filename not specified!" << endl << "Using default output name RazorInclusive.root" << endl;
    outFileName = "RazorInclusive.root";
  }
  TFile outFile(outFileName.c_str(), "RECREATE");
    
  //one tree to hold all events
  TTree *razorTree = new TTree("RazorInclusive", "Info on selected razor inclusive events");

  //For signal samples, create one output file and tree per signal mass point
  map<pair<int,int>, TFile*> smsFiles;
  map<pair<int,int>, TTree*> smsTrees;
  map<pair<int,int>, TH1F*> smsNEvents;
    
  //initialize jet energy corrections
  TRandom3 *random = new TRandom3(33333); //Artur wants this number 33333
  std::vector<JetCorrectorParameters> correctionParameters;
  //get correct directory for JEC files (different for lxplus and t3-higgs)
  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/RazorAnalyzer/data/JEC/";
  cout << "Getting JEC parameters from " << pathname << endl;
  
  if (isData) {
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt", pathname.c_str())));
  } else if (isFastsimSMS) {
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Fastsim_MCRUN2_74_V9_L1FastJet_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Fastsim_MCRUN2_74_V9_L2Relative_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Fastsim_MCRUN2_74_V9_L3Absolute_AK4PFchs.txt", pathname.c_str())));
  } else {
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));
  }

  
  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);
  JetCorrectorParameters *JetResolutionParameters = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",pathname.c_str()));
  SimpleJetResolution *JetResolutionCalculator = new SimpleJetResolution(*JetResolutionParameters);

  //separate trees for individual boxes
  map<string, TTree*> razorBoxes;
  vector<string> boxNames;
  boxNames.push_back("MuEle");
  boxNames.push_back("MuMu");
  boxNames.push_back("EleEle");
  boxNames.push_back("MuMultiJet");
  boxNames.push_back("MuJet");
  boxNames.push_back("EleMultiJet");
  boxNames.push_back("EleJet");
  boxNames.push_back("LooseLeptonMultiJet");
  boxNames.push_back("MultiJet");
  boxNames.push_back("DiJet");

  for(size_t i = 0; i < boxNames.size(); i++){
    razorBoxes[boxNames[i]] = new TTree(boxNames[i].c_str(), boxNames[i].c_str());
  }

  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

  //tree variables
  int nVtx;
  int nSelectedJets, nBTaggedJets;
  int nJets80;
  int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
  int nVetoMuons, nVetoElectrons, nLooseTaus;
  float dPhiRazor;
  float theMR;
  float theRsq;
  float MRForSF;
  float RsqForSF;
  float met;
  float HT;
  float mT, mTLoose, leadingJetPt, subleadingJetPt, leadingTightMuPt, leadingTightElePt;
  float weight = 1.0;
  int mGluino, mLSP;
  float leadingGenJetPt, subleadingGenJetPt;
  //float pileupWeight = 1.0;
  //float lepEffCorrFactor = 1.0;
  //float lepTrigCorrFactor = 1.0;
  //float btagCorrFactor = 1.0;
  //bool  hltDecision[100];
  int nGenMuons, nGenElectrons, nGenTauMuons, nGenTauElectrons, nGenHadTaus, nGenTaus;
  float leadingGenMuonPt, leadingGenElectronPt, leadingGenTauPt;
  float leadingGenMuonEta, leadingGenElectronEta, leadingGenTauEta;
  float minDRGenLeptonToGenParton;
  float genmet;
  float genHt;
  float minDPhiMetToJets;
  float leadingMuonPt;
  float allMuonPt;

  RazorBox box;

  //set branches on big tree
  if(combineTrees){
    razorTree->Branch("nVtx", &nVtx, "nVtx/I");
    razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
    razorTree->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
    razorTree->Branch("nJets80", &nJets80, "nJets80/I");
    razorTree->Branch("nVetoMuons", &nVetoMuons, "nVetoMuons/I");
    razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
    razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
    razorTree->Branch("nVetoElectrons", &nVetoElectrons, "nVetoElectrons/I");
    razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
    razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
    razorTree->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
    razorTree->Branch("MR", &theMR, "MR/F");
    razorTree->Branch("dPhiRazor", &dPhiRazor, "dPhiRazor/F");
    razorTree->Branch("Rsq", &theRsq, "Rsq/F");
    razorTree->Branch("MRForSF", &MRForSF, "MRForSF/F");
    razorTree->Branch("RsqForSF", &RsqForSF, "RsqForSF/F");
    razorTree->Branch("met", &met, "met/F");
    razorTree->Branch("HT", &HT, "HT/F");
    razorTree->Branch("mT", &mT, "mT/F");
    razorTree->Branch("mTLoose", &mTLoose, "mTLoose/F");//for LooseLepton boxes
    razorTree->Branch("leadingTightMuPt", &leadingTightMuPt, "leadingTightMuPt/F");
    razorTree->Branch("leadingTightElePt", &leadingTightElePt, "leadingTightElePt/F");
    razorTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/F");
    razorTree->Branch("subleadingJetPt", &subleadingJetPt, "subleadingJetPt/F");
    razorTree->Branch("box", &box, "box/I");
    razorTree->Branch("HLTDecision", &HLTDecision, "HLTDecision[300]/O");
    razorTree->Branch("minDPhiMetToJets", &minDPhiMetToJets, "minDPhiMetToJets/F");
    razorTree->Branch("leadingMuonPt", &leadingMuonPt, "leadingMuonPt/F");
    razorTree->Branch("allMuonPt", &allMuonPt, "allMuonPt/F");


    //MET Filters
    razorTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
    razorTree->Branch("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, "Flag_HBHEIsoNoiseFilter/O");
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


    if (!isData) {    
      razorTree->Branch("weight", &weight, "weight/F");
      razorTree->Branch("leadingGenJetPt", &leadingGenJetPt, "leadingGenJetPt/F");
      razorTree->Branch("subleadingGenJetPt", &subleadingGenJetPt, "subleadingGenJetPt/F");
      //   razorTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
      //   razorTree->Branch("lepEffCorrFactor", &lepEffCorrFactor, "lepEffCorrFactor/F");
      //   razorTree->Branch("lepTrigCorrFactor", &lepTrigCorrFactor, "lepTrigCorrFactor/F");
      //   razorTree->Branch("btagCorrFactor", &btagCorrFactor, "btagCorrFactor/F");
      if (isFastsimSMS) {
	razorTree->Branch("mGluino", &mGluino, "mGluino/I");
	razorTree->Branch("mLSP", &mLSP, "mLSP/I");
      }

      //For extra MC information
      razorTree->Branch("nGenMuons", &nGenMuons, "nGenMuons/I");
      razorTree->Branch("nGenElectrons", &nGenElectrons, "nGenElectrons/I");
      razorTree->Branch("nGenTauMuons", &nGenTauMuons, "nGenTauMuons/I");
      razorTree->Branch("nGenTauElectrons", &nGenTauElectrons, "nGenTauElectrons/I");
      razorTree->Branch("nGenHadTaus", &nGenHadTaus, "nGenHadTaus/I");
      razorTree->Branch("nGenTaus", &nGenTaus, "nGenTaus/I");
      razorTree->Branch("leadingGenMuonPt", &leadingGenMuonPt, "leadingGenMuonPt/F");
      razorTree->Branch("leadingGenElectronPt", &leadingGenElectronPt, "leadingGenElectronPt/F");
      razorTree->Branch("leadingGenTauPt", &leadingGenTauPt, "leadingGenTauPt/F");
      razorTree->Branch("leadingGenMuonEta", &leadingGenMuonEta, "leadingGenMuonEta/F");
      razorTree->Branch("leadingGenElectronEta", &leadingGenElectronEta, "leadingGenElectronEta/F");
      razorTree->Branch("leadingGenTauEta", &leadingGenTauEta, "leadingGenTauEta/F");
      razorTree->Branch("minDRGenLeptonToGenParton", &minDRGenLeptonToGenParton, "minDRGenLeptonToGenParton/F");
      razorTree->Branch("genmet", &genmet, "genmet/F");
      razorTree->Branch("genHt", &genHt, "genHt/F");
     
    } else {
      razorTree->Branch("run", &runNum, "run/i");
      razorTree->Branch("lumi", &lumiNum, "lumi/i");
      razorTree->Branch("event", &eventNum, "event/i");
    }
  }
  //set branches on all trees
  else{ 
    for(auto& box : razorBoxes){
      box.second->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
      box.second->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
      box.second->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
      box.second->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
      box.second->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
      box.second->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
      box.second->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
      box.second->Branch("dPhiRazor", &dPhiRazor, "dPhiRazor/F");
      box.second->Branch("MR", &theMR, "MR/F");
      box.second->Branch("met", &met, "met/F");
      box.second->Branch("HT", &HT, "HT/F");
      box.second->Branch("Rsq", &theRsq, "Rsq/F");
      box.second->Branch("mT", &mT, "mT/F");
      box.second->Branch("mTLoose", &mTLoose, "mTLoose/F");
      box.second->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/F");
      box.second->Branch("subleadingJetPt", &subleadingJetPt, "subleadingJetPt/F");
      box.second->Branch("leadingTightMuPt", &leadingTightMuPt, "leadingTightMuPt/F");
      box.second->Branch("leadingTightElePt", &leadingTightElePt, "leadingTightElePt/F");
      if (!isData) {    
	box.second->Branch("weight", &weight, "weight/F");
	box.second->Branch("leadingGenJetPt", &leadingGenJetPt, "leadingGenJetPt/F");
	box.second->Branch("subleadingGenJetPt", &subleadingGenJetPt, "subleadingGenJetPt/F");
	//   box.second->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
	//   box.second->Branch("lepEffCorrFactor", &lepEffCorrFactor, "lepEffCorrFactor/F");
	//   box.second->Branch("lepTrigCorrFactor", &lepTrigCorrFactor, "lepTrigCorrFactor/F");
	//   box.second->Branch("btagCorrFactor", &btagCorrFactor, "btagCorrFactor/F");
	if(isFastsimSMS){
	  box.second->Branch("mGluino", &mGluino, "mGluino/I");
	  box.second->Branch("mLSP", &mLSP, "mLSP/I");
	}
      } else {
	box.second->Branch("run", &runNum, "run/i");
	box.second->Branch("lumi", &lumiNum, "lumi/i");
	box.second->Branch("event", &eventNum, "event/i");
      }
    }
  }

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
    printdebug = false;

    //fill normalization histogram
    NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);

    //reset tree variables
    nVtx = nPV;
    nSelectedJets = 0;
    nBTaggedJets = 0;
    nVetoMuons = 0;
    nLooseMuons = 0;
    nTightMuons = 0;
    nVetoElectrons = 0;
    nLooseElectrons = 0;
    nTightElectrons = 0;
    nLooseTaus = 0;
    nTightTaus = 0;
    theMR = -1;
    theRsq = -1;
    MRForSF = -1;
    RsqForSF = -1;
    mT = -1;
    mTLoose = -1;
    leadingJetPt = -1;
    subleadingJetPt = -1;
    leadingTightMuPt = -1;
    leadingTightElePt = -1;
    if(combineTrees) box = NONE;
    weight = genWeight;
    if(isFastsimSMS){
      mGluino = -1;
      mLSP = -1;
    }
    nGenMuons = 0;
    nGenElectrons = 0;
    nGenTauMuons = 0;
    nGenTauElectrons = 0;
    nGenHadTaus = 0;
    nGenTaus = 0;
    leadingGenMuonPt = 0;
    leadingGenElectronPt = 0;
    leadingGenTauPt = 0;
    leadingGenMuonEta = -999;
    leadingGenElectronEta = -999;
    leadingGenTauEta = -999;
    minDRGenLeptonToGenParton = 9999;
    genmet = -999;
    genHt  = 0.;

    /////////////////////////////////
    //SMS information
    /////////////////////////////////

    bool parsedLHE = false;
    if(isFastsimSMS && lheComments){
      //parse lhe comment string to get gluino and LSP masses
      stringstream parser((*lheComments));
      string item;
      getline(parser, item, '_'); //prefix
      if(getline(parser, item, '_')){ //gluino mass 
	mGluino = atoi(item.c_str());
	if(getline(parser, item, '_')){ //LSP mass 
	  mLSP = atoi(item.c_str());
	  pair<int,int> smsPair = make_pair(mGluino, mLSP);
	  parsedLHE = true;
	  if (smsFiles.count(smsPair) == 0){ //create file and tree
	    //format file name
	    string thisFileName = outFileName;
	    thisFileName.erase(thisFileName.end()-5, thisFileName.end());
	    thisFileName += "_" + to_string(mGluino) + "_" + to_string(mLSP) + ".root";

	    smsFiles[smsPair] = new TFile(thisFileName.c_str(), "recreate");
	    smsTrees[smsPair] = razorTree->CloneTree(0);
	    smsNEvents[smsPair] = new TH1F(Form("NEvents%d%d", mGluino, mLSP), "NEvents", 1,1,2);
	    cout << "Created new output file " << thisFileName << endl;
	  }
	  //Fill NEvents hist 
	  smsNEvents[smsPair]->Fill(1.0);
	}
      }
    }

    //*****************************************
    //MC Information
    //*****************************************
    genmet = genMetPt;
    // calculate gen Ht
    for(int j = 0; j < nGenParticle; j++){
      
      if ( (abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) ||  (abs(gParticleId[j]) == 21) )
	if (gParticleStatus[j] == 23 )
          genHt += gParticlePt[j];
    }

    for(int j = 0; j < nGenParticle; j++){
      if (abs(gParticleId[j]) == 11 && gParticleStatus[j] == 1 	      
	  //&& abs(gParticleEta[j]) < 2.5 && gParticlePt[j] > 5
	  ) {
	if (  (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23) ) {
	  nGenElectrons++;
	  if (gParticlePt[j] > leadingGenElectronPt) {
	    leadingGenElectronPt = gParticlePt[j];
	    leadingGenElectronEta = gParticleEta[j];
	  }
	}
	if ( abs(gParticleMotherId[j]) == 15) {
	  nGenTauElectrons++;
	}
      }
      if (abs(gParticleId[j]) == 13 && gParticleStatus[j] == 1  
	  //&& abs(gParticleEta[j]) < 2.4 && gParticlePt[j] > 5
	  ) {
	if ( (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23)) {
	  nGenMuons++;
	  if (gParticlePt[j] > leadingGenMuonPt) {
	    leadingGenMuonPt = gParticlePt[j];
	    leadingGenMuonEta = gParticleEta[j];
	  }
	}
	if ( abs(gParticleMotherId[j]) == 15) {
	  nGenTauMuons++;
	}
      }
      if (abs(gParticleId[j]) == 15 && gParticleStatus[j] == 2 
	  && (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23)
	  //&& abs(gParticleEta[j]) < 2.4 && gParticlePt[j] > 20
	  ) {
	nGenTaus++;
	bool isLeptonicTau = false;
	for(int k = 0; k < nGenParticle; k++){
	  if ( (abs(gParticleId[k]) == 11 || abs(gParticleId[k]) == 13) && gParticleMotherIndex[k] == j) {
	    isLeptonicTau = true;
	    break;
	  }
	}
	if (!isLeptonicTau) nGenHadTaus++;
	if (!isLeptonicTau && gParticlePt[j] > leadingGenTauPt) {
	  leadingGenTauPt = gParticlePt[j];
	  leadingGenTauEta = gParticleEta[j];
	}
	   
      }
    }
	
    for(int j = 0; j < nGenParticle; j++){
      float minDRToGenLepton = 9999;
      //int closestLeptonIndex = -1;

      //only look for outgoing partons
      if  (!( ((abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) || abs(gParticleId[j]) == 21) 
	      && gParticleStatus[j] == 23)
	   ) continue;
	       		      
      //look for closest lepton
      for(int k = 0; k < nGenParticle; k++){
	if ( 
	    (abs(gParticleId[k]) == 11 && gParticleStatus[k] == 1 	      
	     && (abs(gParticleMotherId[k]) == 24 || abs(gParticleMotherId[k]) == 23)
	     && abs(gParticleEta[k]) < 2.5 && gParticlePt[k] > 5)
	    ||
	    (abs(gParticleId[k]) == 13 && gParticleStatus[k] == 1  
	     && (abs(gParticleMotherId[k]) == 24 || abs(gParticleMotherId[k]) == 23)
	     && abs(gParticleEta[k]) < 2.4 && gParticlePt[k] > 5
	     )
	    ||
	    (abs(gParticleId[k]) == 15 && gParticleStatus[k] == 2 
	     && (abs(gParticleMotherId[k]) == 24 || abs(gParticleMotherId[k]) == 23)
	     && abs(gParticleEta[k]) < 2.4 && gParticlePt[k] > 20
	     )
	     ) {	 
	  double tmpDR = deltaR( gParticleEta[j], gParticlePhi[j], gParticleEta[k], gParticlePhi[k]);
	  if ( tmpDR < minDRToGenLepton ) {
	    minDRToGenLepton = tmpDR;
	    //closestLeptonIndex = k;
	  }
	}
      }

      // cout << "Parton " << j << " : " << gParticleId[j] << " | " << gParticlePt[j] << " " << gParticleEta[j] << " " << gParticlePhi[j] << " | " << gParticleMotherId[j] 
      //      << " : " << minDRToGenLepton 
      //      << " | " ;
      // if (closestLeptonIndex >= 0 && minDRToGenLepton < 0.4) {
      //   cout << gParticleId[closestLeptonIndex] << " " << gParticlePt[closestLeptonIndex] << " " << gParticleEta[closestLeptonIndex] << " " << gParticlePhi[closestLeptonIndex] << " " << gParticleMotherId[closestLeptonIndex];
      // }
      // cout << "\n";

      if ( minDRToGenLepton < minDRGenLeptonToGenParton) minDRGenLeptonToGenParton = minDRToGenLepton;
    }



    //*****************************************
    //Triggers
    //*****************************************
    bool passedDileptonTrigger = false;
    bool passedSingleLeptonTrigger = false;
    bool passedHadronicTrigger= false;

    if (isData) {
      /* 2015 data dilepton triggers:
	 30    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
	 31    HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
	 41    HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
	 43    HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ
	 47    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
	 48    HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
	 49    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL
	 50    HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL
      passedDileptonTrigger = bool( HLTDecision[41] || HLTDecision[43] 
				    || HLTDecision[30] || HLTDecision[31] 
				    || HLTDecision[47] || HLTDecision[48] || HLTDecision[49] || HLTDecision[50] );
      */
      /* 2016 data dilepton triggers:

	 44    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
	 45    HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
	 57    HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
	 59    HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ
	 64    HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL
	 65    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
	 66    HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
	 67    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL
	 68    HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL
      */
      passedDileptonTrigger = bool( HLTDecision[44] || HLTDecision[45] 
				    || HLTDecision[57] || HLTDecision[59] 
				    || HLTDecision[64] || HLTDecision[65] || HLTDecision[66] || HLTDecision[67] || HLTDecision[68]);
      /* 2015 data single lepton triggers:
	 2    HLT_Mu50
	 7    HLT_IsoMu20
	 11    HLT_IsoMu27
	 12    HLT_IsoTkMu20
	 15    HLT_IsoTkMu27
	 22    HLT_Ele23_WPLoose_Gsf
	 23    HLT_Ele27_WPLoose_Gsf
	 24    HLT_Ele27_eta2p1_WPLoose_Gsf
	 25    HLT_Ele27_eta2p1_WPTight_Gsf
	 26    HLT_Ele32_eta2p1_WPLoose_Gsf
	 27    HLT_Ele32_eta2p1_WPTight_Gsf
	 28    HLT_Ele105_CaloIdVT_GsfTrkIdT
	 29    HLT_Ele115_CaloIdVT_GsfTrkIdT      
      passedSingleLeptonTrigger = bool(HLTDecision[2] || HLTDecision[7] || HLTDecision[12] || HLTDecision[11] || HLTDecision[15]
				       || HLTDecision[22] || HLTDecision[23] || HLTDecision[24] || HLTDecision[25] || 
				       HLTDecision[26] || HLTDecision[27] ||
				       HLTDecision[28] || HLTDecision[29]);
      */
      /* 2016 data single lepton triggers:
	 4    HLT_Mu50
	 6    HLT_Mu45_eta2p1
	 12    HLT_IsoMu18
	 13    HLT_IsoMu20
	 18    HLT_IsoMu27
	 19    HLT_IsoTkMu18
	 20    HLT_IsoTkMu20
	 24    HLT_IsoTkMu27
	 27    HLT_Ele22_eta2p1_WPLoose_Gsf
	 28    HLT_Ele22_eta2p1_WPTight_Gsf
	 29    HLT_Ele23_WPLoose_Gsf
	 30    HLT_Ele24_eta2p1_WPLoose_Gsf
	 31    HLT_Ele25_WPTight_Gsf
	 32    HLT_Ele25_eta2p1_WPLoose_Gsf
	 33    HLT_Ele25_eta2p1_WPTight_Gsf
	 34    HLT_Ele27_WPLoose_Gsf
	 35    HLT_Ele27_WPTight_Gsf
	 36    HLT_Ele27_eta2p1_WPLoose_Gsf
	 37    HLT_Ele27_eta2p1_WPTight_Gsf
	 38    HLT_Ele32_eta2p1_WPLoose_Gsf
	 39    HLT_Ele32_eta2p1_WPTight_Gsf
	 42    HLT_Ele105_CaloIdVT_GsfTrkIdT
	 43    HLT_Ele115_CaloIdVT_GsfTrkIdT
      */
      passedSingleLeptonTrigger = bool(HLTDecision[4] || HLTDecision[6] || HLTDecision[12] || HLTDecision[13] || HLTDecision[18] || HLTDecision[19] ||
				       HLTDecision[20] || HLTDecision[24] ||
				       HLTDecision[27] || HLTDecision[28] || HLTDecision[29] || HLTDecision[30] || HLTDecision[31] || HLTDecision[32] || HLTDecision[33] ||
				       HLTDecision[34] || HLTDecision[35] || HLTDecision[36] || HLTDecision[37] || 
				       HLTDecision[38] || HLTDecision[39] ||
				       HLTDecision[42] || HLTDecision[43]);
      /* 2015 data hadronic triggers:
	 134    HLT_RsqMR240_Rsq0p09_MR200
	 135    HLT_RsqMR240_Rsq0p09_MR200_4jet
	 136    HLT_RsqMR270_Rsq0p09_MR200
	 137    HLT_RsqMR270_Rsq0p09_MR200_4jet
	 138    HLT_RsqMR260_Rsq0p09_MR200
	 139    HLT_RsqMR260_Rsq0p09_MR200_4jet
	 140    HLT_RsqMR300_Rsq0p09_MR200
	 141    HLT_RsqMR300_Rsq0p09_MR200_4jet
	 142    HLT_Rsq0p25
	 143    HLT_Rsq0p30
	 144    HLT_Rsq0p36
      passedHadronicTrigger = bool(HLTDecision[134] || HLTDecision[135] || HLTDecision[136] 
				   || HLTDecision[137] || HLTDecision[138] || HLTDecision[139] 
				   || HLTDecision[140] || HLTDecision[141] || HLTDecision[142] 
				   || HLTDecision[143] || HLTDecision[144]);
      */
      /* 2016 data hadronic trigger:
	 164    HLT_RsqMR240_Rsq0p09_MR200
	 165    HLT_RsqMR240_Rsq0p09_MR200_4jet
	 166    HLT_RsqMR270_Rsq0p09_MR200
	 167    HLT_RsqMR270_Rsq0p09_MR200_4jet
	 168    HLT_RsqMR260_Rsq0p09_MR200
	 169    HLT_RsqMR260_Rsq0p09_MR200_4jet
	 170    HLT_RsqMR300_Rsq0p09_MR200
	 171    HLT_RsqMR300_Rsq0p09_MR200_4jet
	 172    HLT_Rsq0p25
	 173    HLT_Rsq0p30
	 174    HLT_Rsq0p36
	 175    HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200
	 176    HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200
      */
      passedHadronicTrigger = bool(HLTDecision[164] || HLTDecision[165] || HLTDecision[166] 
				   || HLTDecision[167] || HLTDecision[168] || HLTDecision[160] 
				   || HLTDecision[170] || HLTDecision[171] || HLTDecision[172] 
				   || HLTDecision[173] || HLTDecision[174] || HLTDecision[175]
				   || HLTDecision[176] );
    } else {
      passedDileptonTrigger = bool(HLTDecision[41] || HLTDecision[43] 
				   || HLTDecision[30] || HLTDecision[31] 
				   || HLTDecision[47] || HLTDecision[48] || HLTDecision[49] || HLTDecision[50] );
      passedSingleLeptonTrigger = bool( HLTDecision[2] || HLTDecision[7] || HLTDecision[12] 
					|| HLTDecision[11] || HLTDecision[15] 
					|| HLTDecision[18] || HLTDecision[19] || HLTDecision[20] 
					|| HLTDecision[21] || HLTDecision[28] || HLTDecision[29] || HLTDecision[160]);
      passedHadronicTrigger = bool(HLTDecision[134] || HLTDecision[135] || HLTDecision[136] 
				   || HLTDecision[137] || HLTDecision[138] || HLTDecision[139] 
				   || HLTDecision[140] || HLTDecision[141] || HLTDecision[142] 
				   || HLTDecision[143] || HLTDecision[144]);    
    }

    //ignore trigger for Fastsim
    if(isFastsimSMS){
      passedDileptonTrigger = true;
      passedSingleLeptonTrigger = true;
      passedHadronicTrigger = true;
    }

    //*****************************************
    //Get Pileup Information
    //*****************************************
    double NPU = 0;
    for (int i=0; i < nBunchXing; ++i) {
      if (BunchXing[i] == 0) {
	NPU = nPUmean[i];
      }
    }

    // if (pileupWeightHist) {
    //   pileupWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(NPU));
    // }

    //*****************************************
    //Select Leptons
    //*****************************************
    //lepEffCorrFactor  = 1.0;
    leadingMuonPt = 0;
    double allMuonPx = 0;
    double allMuonPy = 0;
    //TLorentzVector allMuonSystem; allMuonSystem.SetPxPyPzE(0,0,0,0);
    vector<TLorentzVector> GoodLeptons; //leptons used to compute hemispheres
    TLorentzVector leadingTightMu, leadingTightEle; //used for mT calculation    
    for(int i = 0; i < nMuons; i++){

      if(muonPt[i] < 5) continue;
      if(abs(muonEta[i]) > 2.4) continue;
     
      //Calculate MC->Data Scale Factors
      if (RazorAnalyzer::matchesGenMuon(muonEta[i],muonPhi[i])) {	
	//apply muon efficiency scale factors
      }
      
      TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 

      if (muonPt[i] > leadingMuonPt) leadingMuonPt = muonPt[i];
      allMuonPx += muonPt[i]*cos(muonPhi[i]);
      allMuonPy += muonPt[i]*sin(muonPhi[i]);

      if(isVetoMuon(i)) nVetoMuons++;
      if(isLooseMuon(i) && muonPt[i] >= 20 ) nLooseMuons++;
      if(isTightMuon(i) && muonPt[i] >= 20){
	nTightMuons++;
	if(muonPt[i] > leadingTightMuPt){
	  leadingTightMu = thisMuon;
	  leadingTightMuPt = muonPt[i];
	}
      }

      if(!isVetoMuon(i)) continue;  
      GoodLeptons.push_back(thisMuon); 
    }
    allMuonPt = sqrt( pow(allMuonPx,2) + pow(allMuonPy,2));

    for(int i = 0; i < nElectrons; i++){
      if(elePt[i] < 5) continue;
      if(fabs(eleEta[i]) > 2.5) continue;

      if (RazorAnalyzer::matchesGenElectron(eleEta[i],elePhi[i])) {
	//No Efficiency Scale Factors for Run2 Yet

	//   double effLoose = getElectronEfficiency("loose",elePt[i],eleEta[i]);	  
	//   double effLooseSF = eleLooseEffSFHist->GetBinContent( eleLooseEffSFHist->GetXaxis()->FindFixBin(fabs(eleEta[i])) , 
	// 							eleLooseEffSFHist->GetYaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)));
	//   double tmpLooseSF = 1.0;
	//   if( isLooseElectron(i) ) {
	//     tmpLooseSF = effLooseSF;
	//   } else {
	//     tmpLooseSF = ( 1/effLoose - effLooseSF) / ( 1/effLoose - 1);
	//   }

	//   if (tmpLooseSF != tmpLooseSF) cout << tmpLooseSF << " " << effLoose << " " << effLooseSF << "\n";

	//   double effTight = getElectronEfficiency("tight",elePt[i],eleEta[i]);	  
	//   double effTightSF = eleTightEffSFHist->GetBinContent( eleTightEffSFHist->GetXaxis()->FindFixBin(fabs(eleEta[i])) , 
	// 							eleTightEffSFHist->GetYaxis()->FindFixBin(fmax(fmin(elePt[i],199.9),10.01)));
	//   double tmpTightSF = 1.0;
	//   if( isTightElectron(i) ) {
	//     tmpTightSF = effTightSF;
	//   } else {
	//     tmpTightSF = ( 1/effTight - effTightSF) / ( 1/effTight - 1);
	//   }
	//   if (tmpTightSF != tmpTightSF) cout << tmpTightSF << " " << effTight << " " << effTightSF << "\n";

	//   lepEffCorrFactor *= tmpLooseSF;
	//   lepEffCorrFactor *= tmpTightSF;
      }

      //remove overlaps
      bool overlap = false;
      for(auto& lep : GoodLeptons){
	if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
      }
      if(overlap) continue;

      //TLorentzVector for this electron
      TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);

      if(isVetoElectron(i)) nVetoElectrons++;
      if( isLooseElectron(i) && elePt[i] > 25 ) nLooseElectrons++;
      if( isTightElectron(i) && elePt[i] > 25 ){
	nTightElectrons++;
	if(elePt[i] > leadingTightElePt){
	  leadingTightEle = thisElectron;
	  leadingTightElePt = elePt[i];
	}
      }

      if(!isVetoElectron(i)) continue; 
      GoodLeptons.push_back(thisElectron);            
    }

    //******************************
    //Only Do Taus for Run2
    //******************************
    for(int i = 0; i < nTaus; i++){	 
      if (tauPt[i] < 20) continue;
      if (fabs(tauEta[i]) > 2.4) continue;

      //remove overlaps
      bool overlap = false;
      for(auto& lep : GoodLeptons){
	if (RazorAnalyzer::deltaR(tauEta[i],tauPhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
      }
      if(overlap) continue;

      if(isLooseTau(i)) nLooseTaus++;
      if(isTightTau(i)) nTightTaus++;

      if (!isLooseTau(i)) continue;
      TLorentzVector thisTau; thisTau.SetPtEtaPhiM(tauPt[i], tauEta[i], tauPhi[i], 1.777);
      GoodLeptons.push_back(thisTau);  
    }

    vector<TLorentzVector> GoodJets;
    int numJetsAbove80GeV = 0;

    // //***********************************************
    // //use genjets instead , for debugging
    // //***********************************************
    leadingGenJetPt = 0;
    subleadingGenJetPt = 0;
    for(int j = 0; j < nGenJets; j++){

      if(genJetPt[j] < 40) continue;
      if(fabs(genJetEta[j]) > 3.0) continue;

      //exclude selected muons and electrons from the jet collection
      double deltaR = -1;
      TLorentzVector thisJet = makeTLorentzVector(genJetPt[j], genJetEta[j], genJetPhi[j], genJetE[j]);
      for(auto& lep : GoodLeptons){
     	double thisDR = thisJet.DeltaR(lep);
     	if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      }
      if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

      if(genJetPt[j] > leadingGenJetPt){
	subleadingGenJetPt = leadingGenJetPt;
	leadingGenJetPt = genJetPt[j];
      }
      else if(genJetPt[j] > subleadingGenJetPt){
	subleadingGenJetPt = genJetPt[j];
      }

      //   if(genJetPt[j] > 80) numJetsAbove80GeV++;
      //   GoodJets.push_back(thisJet);
      //   nSelectedJets++;

      //   bool isBJet = false;
      //   for(int i = 0; i < nJets; i++){
      // 	double tmpDR = RazorAnalyzer::deltaR( genJetEta[j], genJetPhi[j], jetEta[i], jetPhi[i] );
      // 	if ( tmpDR < 0.4 && abs(jetPartonFlavor[i]) == 5) isBJet = true;
      //   }
      //   if(isBJet) nBTaggedJets++;

    }


    //initialize B-Tagging Correction Factor
    //btagCorrFactor = 1.0;

    //***********************************************
    //Variables for Type1 Met Correction
    //***********************************************
    double MetX_Type1Corr = 0;
    double MetY_Type1Corr = 0;

    //cout << "Njets: " << nJets << "\n";

    //***********************************************
    //Select Jets
    //***********************************************
    for(int i = 0; i < nJets; i++){

      //*****************************************************************
      //exclude selected muons and electrons from the jet collection
      //*****************************************************************
      double deltaR = -1;
      for(auto& lep : GoodLeptons){
	double thisDR = RazorAnalyzer::deltaR(jetEta[i],jetPhi[i],lep.Eta(),lep.Phi());  
	if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      }
      if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton


      //*****************************************************************
      //apply Jet ID
      //*****************************************************************
      if (!jetPassIDTight[i]) continue;

      //*****************************************************************
      //Apply Jet Energy and Resolution Corrections
      //*****************************************************************
      double tmpRho = fixedGridRhoFastjetAll;
      double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
					     tmpRho, jetJetArea[i], 
					     JetCorrector);   
      double JECLevel1 = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
						   tmpRho, jetJetArea[i], 
						   JetCorrector, 0);   

      double jetEnergySmearFactor = 1.0;
      if (!isData) {
	jetEnergySmearFactor = JetEnergySmearingFactor( jetPt[i]*JEC, jetEta[i], NPU, JetResolutionCalculator, random);
	jetEnergySmearFactor = 1.0; //turn this off for now
      }

      TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);
      TLorentzVector L1CorrJet = makeTLorentzVector(jetPt[i]*JECLevel1, jetEta[i], jetPhi[i], jetE[i]*JECLevel1);
      TLorentzVector UnCorrJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);      
      double jetCorrPt = jetPt[i]*JEC*jetEnergySmearFactor;
      //double jetCorrE = jetE[i]*JEC*jetEnergySmearFactor;

      //*******************************
      //B-Tagging Correction Factor : no Run2 corrections yet
      //*******************************
      // if (abs(jetPartonFlavor[i]) == 5 &&jetCorrPt > 20) {
      // 	double tmpBTagCorrFactor = 1.0;

      // 	double tmpCorrFactor = 0.938887 + 0.00017124 * jetCorrPt + (-2.76366e-07) * jetCorrPt * jetCorrPt ;
      // 	double MCEff = 1.0;
      // 	if (jetCorrPt < 50) MCEff = 0.65;
      // 	else if (jetCorrPt < 80) MCEff = 0.70;
      // 	else if (jetCorrPt < 120) MCEff = 0.73;
      // 	else if (jetCorrPt < 210) MCEff = 0.73;
      // 	else MCEff = 0.66;				 

      // 	//if pass CSV Medium
      // 	if( isCSVM(i) ) {
      // 	  tmpBTagCorrFactor = tmpCorrFactor;
      // 	} else {
      // 	  tmpBTagCorrFactor = ( 1/MCEff - tmpCorrFactor) / ( 1/MCEff - 1);
      // 	}

      // 	btagCorrFactor *= tmpBTagCorrFactor;
      // }

      //*******************************
      //Custom Calculated Type1 Met Correction
      //*******************************
      if (jetCorrPt > 15 && 
	  jetChargedEMEnergyFraction[i] + jetNeutralEMEnergyFraction[i] <= 0.9
	  ) {
	MetX_Type1Corr += -1 * ( thisJet.Px() - L1CorrJet.Px()  );
	MetY_Type1Corr += -1 * ( thisJet.Py() - L1CorrJet.Py()  );
      }

      //*******************************************************
      //apply  Pileup Jet ID : No working point yet for Run2
      //*******************************************************
      // int level = 2; //loose jet ID
      // if (!((jetPileupIdFlag[i] & (1 << level)) != 0)) continue;

      //*******************************************************
      //apply Jet cuts
      //*******************************************************
      if(jetCorrPt < 40) continue;
      if(fabs(jetEta[i]) > 3.0) continue;

      if(jetCorrPt > 80) numJetsAbove80GeV++;
      GoodJets.push_back(thisJet);
      nSelectedJets++;

      if( isCSVM(i) ){ 
	nBTaggedJets++;
      }
    }

    nJets80 = numJetsAbove80GeV;
    //if(numJetsAbove80GeV < 2) continue; //event fails to have two 80 GeV jets    

    //get leading and subleading jet pt
    for(auto &jet : GoodJets){
      if(jet.Pt() > leadingJetPt){
	subleadingJetPt = leadingJetPt;
	leadingJetPt = jet.Pt();
      }
      else if(jet.Pt() > subleadingJetPt){
	subleadingJetPt = jet.Pt();
      }
    }

    //*************************************************************
    //Compute the razor variables using the selected jets and possibly leptons
    //*************************************************************
    vector<TLorentzVector> GoodPFObjects;
    for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
    for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);

    //*************************************************************
    //Make Alternative GoodPFObject List for Scale Factor Lookup
    //*************************************************************
    vector<TLorentzVector> GoodPFObjectsForSF;
    for(auto& jet : GoodJets) GoodPFObjectsForSF.push_back(jet);
    for(auto& lep : GoodLeptons) GoodPFObjectsForSF.push_back(lep);
    
    //Find Gen-Leptons within acceptance (consider electrons or muons only, including from tau decay)
    for(int j = 0; j < nGenParticle; j++){

      //Find Gen-Leptons within acceptance
      if ( 
	  //gen electron
	  (abs(gParticleId[j]) == 11 && gParticleStatus[j] == 1 	      
	   && abs(gParticleEta[j]) < 2.5 && gParticlePt[j] > 5
	   && ( abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23 || abs(gParticleMotherId[j]) == 15)
	   )
	  ||
	  //gen muon
	  ( abs(gParticleId[j]) == 13 && gParticleStatus[j] == 1  
	    && abs(gParticleEta[j]) < 2.4 && gParticlePt[j] > 5
	    && ( (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23) || abs(gParticleMotherId[j]) == 15 )
	    )
	   ) {

	bool isSelected = false;
	
	//Find if the gen lepton has been selected as a lepton
	for(auto& lep : GoodLeptons) {
	  if ( deltaR( gParticleEta[j], gParticlePhi[j], lep.Eta(), lep.Phi()) < 0.1 ) isSelected = true;
	}

	//Find if the gen lepton was within a selected jet
	for(auto& jet : GoodJets) {
	  if ( deltaR( gParticleEta[j], gParticlePhi[j], jet.Eta(), jet.Phi()) < 0.4 ) isSelected = true;
	}
	
	//if the gen lepton was not within a selected lepton nor a selected jet, then include it in the GoodPFObjectsForSF list
	if (!isSelected) {
	  TLorentzVector thisLepton; 
	  thisLepton.SetPtEtaPhiE(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
	  GoodPFObjectsForSF.push_back(thisLepton);

	  // ******************************************************
	  //For debug info
	  // ******************************************************
	  // cout << "\nFound gen lepton that was not selected: " << gParticleId[j] << " " << gParticlePt[j] << " " << gParticleEta[j] << " " << gParticlePhi[j] << " : " << gParticleMotherId[j] << "\n";
	  // for(int i = 0; i < nMuons; i++){
	  //   if(muonPt[i] < 5) continue;
	  //   if(abs(muonEta[i]) > 2.4) continue;
     
	  //   cout << "Muon Candidate : " << muonPt[i] << " " << muonEta[i] << " " << muonPhi[i] << " " 
	  // 	 << deltaR( gParticleEta[j], gParticlePhi[j], muonEta[i], muonPhi[i]) << " : "
	  // 	 << muonType[i] << " : "
	  // 	 << isTightMuon(i) << " " << isLooseMuon(i) << " " << isVetoMuon(i) << " : "
	  // 	 << isTightMuon(i,true,false) << " " << isLooseMuon(i,true,false) << " " << isVetoMuon(i,true,false) << " : "
	  // 	 << isTightMuon(i,false,true) << " " << isLooseMuon(i,false,true) << " " << isVetoMuon(i,false,true) << " : "
	  // 	 << "\n";
	  // }

	  // for(int i = 0; i < nElectrons; i++){
	  //   if(elePt[i] < 5) continue;
	  //   if(fabs(eleEta[i]) > 2.5) continue;
	  //   cout << "Ele Candidate : " << elePt[i] << " " << eleEta[i] << " " << elePhi[i] << " " 
	  // 	 << deltaR( gParticleEta[j], gParticlePhi[j], eleEta[i], elePhi[i] ) << " : "
	  // 	 << isTightElectron(i) << " " << isLooseElectron(i) << " " << isVetoElectron(i) << " "
	  // 	 << isTightElectron(i,true,false) << " " << isLooseElectron(i,true,false) << " " << isVetoElectron(i,true,false) << " "
	  // 	 << isTightElectron(i,false,true) << " " << isLooseElectron(i,false,true) << " " << isVetoElectron(i,false,true) << " "
	  // 	 << "\n";
	  // }
	  
	  // for(int i = 0; i < nTaus; i++){	 
	  //   if (tauPt[i] < 20) continue;
	  //   if (fabs(tauEta[i]) > 2.4) continue;

	  //   cout << "Tau Candidate : " << tauPt[i] << " " << tauEta[i] << " " << tauPhi[i] << " " 
	  // 	 << deltaR( gParticleEta[j], gParticlePhi[j], tauEta[i], tauPhi[i] ) << " : "
	  // 	 << isTightTau(i) << " " << isLooseTau(i) << " "
	  // 	 << "\n";
	  // }
	  
	  // for(int i = 0; i < nJets; i++){
	  //   cout << "Jet Candidate : " << jetPt[i] << " " << jetEta[i] << " " << jetPhi[i] << " " 
	  // 	 << deltaR( gParticleEta[j], gParticlePhi[j], jetEta[i], jetPhi[i] ) << " : "
	  // 	 << jetPassIDTight[i] << " : " << isCSVM(i) << " "
	  // 	 << "\n"; 
	  // }

	  // for(auto& lep : GoodLeptons) {
	  //   cout << "GoodLepton : " << lep.M() << " " << lep.Pt() << " " << lep.Eta() << " " << lep.Phi() << " : " << deltaR( gParticleEta[j], gParticlePhi[j], lep.Eta(), lep.Phi()) << "\n";
	  // }	  

	  // for(auto& jet : GoodJets) {
	  //   cout << "GoodJet : " << jet.M() << " " << jet.Pt() << " " << jet.Eta() << " " << jet.Phi() << " : " << deltaR( gParticleEta[j], gParticlePhi[j], jet.Eta(), jet.Phi()) << "\n";
	  // }
	  
	}

      }
    }


    //*************************************************************
    //Apply Type1 Met Correction
    //*************************************************************
    double PFMetCustomType1CorrectedX = metPt*cos(metPhi) + MetX_Type1Corr;
    double PFMetCustomType1CorrectedY = metPt*sin(metPhi) + MetY_Type1Corr;
    TLorentzVector PFMETCustomType1Corrected; 
    PFMETCustomType1Corrected.SetPxPyPzE(PFMetCustomType1CorrectedX, PFMetCustomType1CorrectedY, 0, 
					 sqrt( pow(PFMetCustomType1CorrectedX,2) + pow(PFMetCustomType1CorrectedY,2)));      
    TLorentzVector PFMETUnCorr = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
    TLorentzVector PFMETType1 = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);
    TLorentzVector PFMETType0Plus1 = makeTLorentzVectorPtEtaPhiM(metType0Plus1Pt, 0, metType0Plus1Phi, 0);
    TLorentzVector PFMETNoHF = makeTLorentzVectorPtEtaPhiM(metNoHFPt, 0, metNoHFPhi, 0);
    TLorentzVector MyMET = PFMETCustomType1Corrected; //This is the MET that will be used below.

    HT = 0;
    for(auto& obj : GoodPFObjects) HT += obj.Pt();
    
    //compute minimum dPhi Met to Jets
    minDPhiMetToJets = 9999;
    for (int k=0; k < int(GoodJets.size()); ++k) {
      if (k>=4) break;
      if ( fabs(deltaPhi(MyMET.Phi(), GoodJets[k].Phi())) < minDPhiMetToJets ) {
	minDPhiMetToJets = fabs(deltaPhi(MyMET.Phi(), GoodJets[k].Phi()));
      }
    }

    if ( GoodPFObjects.size() < 20) {
      vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
      theMR = computeMR(hemispheres[0], hemispheres[1]); 
      theRsq = computeRsq(hemispheres[0], hemispheres[1], MyMET);
      dPhiRazor = deltaPhi(hemispheres[0].Phi(),hemispheres[1].Phi());

      //Compute alternative MR and Rsq to be used for Scale Factor lookup
      vector<TLorentzVector> hemispheresForSF = getHemispheres(GoodPFObjectsForSF);
      MRForSF = computeMR(hemispheresForSF[0], hemispheresForSF[1]); 
      RsqForSF = computeRsq(hemispheresForSF[0], hemispheresForSF[1], MyMET);
    } else {
      theMR = -999;
      theRsq = -999;
      dPhiRazor = -999;
      cout << "WARNING: Event has  more than 20 objects\n";
      cout << "NElectrons =  " << nVetoElectrons << "\n";
      cout << "NMuons =  " << nVetoMuons << "\n";
      cout << "NTaus = " << nLooseTaus << "\n";
      for(auto& lep : GoodLeptons) {
	cout << "lepton : " << lep.Pt() << " " << lep.Eta() << " " << lep.Phi() << " " << lep.M() << "\n";
      }
      for(auto& jet : GoodJets) {
	cout << "jet: " << jet.Pt() << " " << jet.Eta() << " " << jet.Phi() << " " << jet.M() << "\n";
      }      
    }
    met = MyMET.Pt();

    //save transverse mass 
    if(nTightMuons + nTightElectrons > 0){
      TLorentzVector leadingLepton;
      if(leadingTightMuPt > leadingTightElePt) leadingLepton = leadingTightMu;
      else leadingLepton = leadingTightEle;

      float deltaPhiLepMet = leadingLepton.DeltaPhi(MyMET);
      mT = sqrt(2*leadingLepton.Pt()*MyMET.Pt()*( 1.0 - cos( deltaPhiLepMet ) ) ); 
    }
    //mT with leading lepton, regardless of quality
    if(GoodLeptons.size() > 0){
      //get the highest pt lepton
      float maxLepPt = -1;
      int maxLepIndex = -1;
      for(uint i = 0; i < GoodLeptons.size(); i++){
	if(GoodLeptons[i].Pt() > maxLepPt){
	  maxLepPt = GoodLeptons[i].Pt();
	  maxLepIndex = i;
	}
      }
      if(maxLepIndex >= 0){
	float deltaPhiLepMet = GoodLeptons[maxLepIndex].DeltaPhi(MyMET);
	mTLoose = sqrt(2*GoodLeptons[maxLepIndex].Pt()*MyMET.Pt()*( 1.0 - cos( deltaPhiLepMet ) ) );
      }
    }

    //**********************************************************************
    //Apply ECAL Dead Cells Filter : Not fixed in miniAOD yet
    //**********************************************************************
    //if (Flag_EcalDeadCellTriggerPrimitiveFilter == false) continue;

    //**********************************************************************
    //Compute correction factor weight : no corr factors for Run2 yet
    //**********************************************************************
    // weight *= pileupWeight;
    // weight *= lepEffCorrFactor;
    // weight *= btagCorrFactor;    

    //**********************************************************************
    //Categorize Events into Razor Boxes 
    //**********************************************************************

    //MuEle Box
    if(passedDileptonTrigger && nTightElectrons > 0 && nTightMuons > 0 ){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = MuEle;
	  razorTree->Fill();
	}
	else razorBoxes["MuEle"]->Fill();	
      }
    }
    //MuMu Box
    else if(passedDileptonTrigger && nTightMuons > 1 ){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = MuMu;
	  razorTree->Fill();
	}
	else razorBoxes["MuMu"]->Fill();	
      }
    }
    //EleEle Box
    else if(passedDileptonTrigger && nTightElectrons > 1 ){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = EleEle;
	  razorTree->Fill();
	}
	else razorBoxes["EleEle"]->Fill();
      }
    }
    //MuSixJet Box
    else if(passedSingleLeptonTrigger && nTightMuons > 0 && nSelectedJets > 5){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = MuSixJet;
	  razorTree->Fill();
	}
	else razorBoxes["MuSixJet"]->Fill();	
      }     
    }
    //MuFourJet Box
    else if(passedSingleLeptonTrigger && nTightMuons > 0 && nSelectedJets > 3){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = MuFourJet;
	  razorTree->Fill();
	}
	else razorBoxes["MuFourJet"]->Fill();	
      }     
    }
    //MuJet Box
    else if(passedSingleLeptonTrigger && nTightMuons > 0 ){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = MuJet;
	  razorTree->Fill();
	}
	else razorBoxes["MuJet"]->Fill();
      }     
    }
    //EleSixJet Box
    else if(passedSingleLeptonTrigger && nTightElectrons > 0 && nSelectedJets > 5 ) {
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = EleSixJet;
	  razorTree->Fill();
	}
	else razorBoxes["EleSixJet"]->Fill();	
      }     
    }
    //EleMultiJet Box
    else if(passedSingleLeptonTrigger && nTightElectrons > 0 && nSelectedJets > 3 ){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = EleFourJet;
	  razorTree->Fill();
	}
	else razorBoxes["EleFourJet"]->Fill();	
      }     
    }
    //EleJet Box
    else if(passedSingleLeptonTrigger && nTightElectrons > 0 ){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = EleJet;
	  razorTree->Fill();
	}
	else razorBoxes["EleJet"]->Fill();
      }     
    }
    //Soft Lepton + SixJet Box
    else if(passedHadronicTrigger && nJets80 >= 2 && nLooseTaus + nVetoElectrons + nVetoMuons > 0 && nSelectedJets > 5){
      if(passesHadronicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = LooseLeptonSixJet;
	  razorTree->Fill();
	}
	else razorBoxes["LooseLeptonSixJet"]->Fill();
      }     
    }
    //Soft Lepton + FourJet Box
    else if(passedHadronicTrigger && nJets80 >= 2 && nLooseTaus + nVetoElectrons + nVetoMuons > 0 && nSelectedJets > 3){
      if(passesHadronicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = LooseLeptonFourJet;
	  razorTree->Fill();
	}
	else razorBoxes["LooseLeptonFourJet"]->Fill();
      }     
    }
    //SixJet Box
    else if(passedHadronicTrigger && nJets80 >= 2 && nSelectedJets > 5 && nLooseTaus+nVetoElectrons+nVetoMuons+nTightElectrons+nTightMuons == 0){
      if(passesHadronicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = SixJet;
	  razorTree->Fill();
	}
	else razorBoxes["SixJet"]->Fill();
      }
    }
    //MultiJet Box
    else if(passedHadronicTrigger && nJets80 >= 2 && nSelectedJets > 3 && nLooseTaus+nVetoElectrons+nVetoMuons+nTightElectrons+nTightMuons == 0){
      if(passesHadronicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = FourJet;
	  razorTree->Fill();
	}
	else razorBoxes["FourJet"]->Fill();
      }
    }
    //Loose Lepton + DiJet Box
    else if(passedHadronicTrigger && nJets80 >= 2 && nLooseTaus + nVetoElectrons + nVetoMuons > 0){
      if(passesHadronicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = LooseLeptonDiJet;
	  razorTree->Fill();
	}
	else razorBoxes["LooseLeptonDiJet"]->Fill();
      }     
    } else if (passedHadronicTrigger && nJets80 >= 2  && nLooseTaus+nVetoElectrons+nVetoMuons+nTightElectrons+nTightMuons == 0) {
      if(passesHadronicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = DiJet;
	  razorTree->Fill();
	}
	else razorBoxes["DiJet"]->Fill();
      }     
    }

    //Fill Fastsim tree
    if(parsedLHE){
      pair<int,int> smsPair = make_pair(mGluino, mLSP);
      smsTrees[smsPair]->Fill();
    }

    //******************************
    //Print Debug
    //******************************
    //printdebug = true;
    if (printdebug) {
      cout << "\nEvent : " << runNum << " " << lumiNum << " " << eventNum << "\n";
      cout << "MR Rsq = " << theMR << " " <<  theRsq << " " << MyMET.Pt()  << "\n";
      cout << "Triggers : " << HLTDecision[134] << " " << HLTDecision[135] << " "<< HLTDecision[136] << " "<< HLTDecision[137] << " : "
	   << HLTDecision[142] << " " << HLTDecision[143] << " "<< HLTDecision[144] << " : " 
	   << HLTDecision[119] << " " << HLTDecision[120] << " "<< HLTDecision[130] << " " 
	   << HLTDecision[131] << " " << HLTDecision[132] << " "<< HLTDecision[133] << " " 
	   << "\n";
      cout << "NElectrons =  " << nTightElectrons << " " << nVetoElectrons << "\n";
      cout << "NMuons =  " << nTightMuons << " " << nVetoMuons << "\n";
      cout << "NTaus = " << nLooseTaus << "\n";
      
      for(auto& lep : GoodLeptons) {
	cout << "good lepton: " << lep.Pt() << " " << lep.Eta() << " " << lep.Phi() << " " << lep.M() << "\n";
      }
      for(auto& jet : GoodJets) {
	cout << "good jet: " << jet.Pt() << " " << jet.Eta() << " " << jet.Phi() << " " << jet.M() << "\n";
      }
    
      for(int i = 0; i < nMuons; i++){
	if(muonPt[i] < 5) continue;
        if(abs(muonEta[i]) > 2.4) continue;
	cout << "Muon " << i << " : " << muonPt[i] << " " << muonEta[i] << " " << muonPhi[i] << " : " << isVetoMuon(i) << " " << isLooseMuon(i)  << " " << isTightMuon(i) << "\n";
      }

      for(int i = 0; i < nElectrons; i++){
        if(elePt[i] < 5) continue;
        if(fabs(eleEta[i]) > 2.5) continue;
	cout << "Electron " << i << " : " << elePt[i] << " " << eleEta[i] << " " << elePhi[i] << " : " << isVetoElectron(i) << " " << isLooseElectron(i)  << " " << isTightElectron(i) << "\n";	
      }
      
      for(int i = 0; i < nTaus; i++){	 
        if (tauPt[i] < 20) continue;
        if (fabs(tauEta[i]) > 2.4) continue;
	cout << "Tau " << i << " : " << tauPt[i] << " " << tauEta[i] << " " << tauPhi[i] << " : " << isLooseTau(i) << " " << isTightTau(i) << "\n";
      }
      
      for(int i = 0; i < nJets; i++){
	if(fabs(jetEta[i]) > 3.0) continue;
	if (!jetPassIDTight[i]) continue;
	double tmpRho = fixedGridRhoFastjetAll;
	double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
					       tmpRho, jetJetArea[i], 
					       JetCorrector);   
        double jetEnergySmearFactor = 1.0;
	TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);
	cout << "jet " << i << " : " << thisJet.Pt() << " " << thisJet.Eta() << " " << thisJet.Phi() << " : " 
	     << isCSVL(i) << " " << isCSVM(i) << " " << isCSVT(i) << " : "
	     << jetMuonEnergyFraction[i] << " " << jetAllMuonPt[i] << " " << jetAllMuonEta[i] << " " << jetAllMuonPhi[i] << " : "
	     << "\n";       
      }

      for(int j = 0; j < nGenParticle; j++){
	cout << "GenParticle " << j << " : " << gParticleId[j] << " " << gParticleStatus[j] << " " << gParticlePt[j] << " " << gParticleEta[j] << " " << gParticlePhi[j] << " : " << gParticleMotherId[j] << "\n";     
      }
    } //endif print debug
    
    
  }//end of event loop

  if(!isFastsimSMS){
    cout << "Writing output trees..." << endl;
    outFile.cd();
    if(combineTrees) razorTree->Write();
    else for(auto& box : razorBoxes) box.second->Write();
    NEvents->Write();
  }
  else{
    for(auto &filePtr : smsFiles){
      cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
      filePtr.second->cd();
      smsTrees[filePtr.first]->Write();
      smsNEvents[filePtr.first]->Write("NEvents");
    }
  }

  outFile.Close();
  if(isFastsimSMS){
    for(auto &f : smsFiles){
      f.second->Close();
    }
  }
}
