#include "HZZRazor.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>
#include <assert.h>

//ROOT includes
#include "TH1F.h"
#include "TH2D.h"

using namespace std;

void HZZRazor::Analyze(bool IsData, int option, string outFileName, string label)
{
  //initialization: create one TTree for each analysis box 
  cout << "Initializing..." << endl;
  bool printdebug = false;

  // //Pileup Weights
  // TFile *pileupWeightFile = 0;
  // TH1D *pileupWeightHist = 0;
  //   pileupWeightFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Run1PileupWeights.root", "READ");
  //   pileupWeightHist = (TH1D*)pileupWeightFile->Get("PUWeight_Run1");
  //   assert(pileupWeightHist);

  // //Lepton Efficiency Correction Factors
  // TH2D *eleLooseEffSFHist = 0;
  // TH2D *eleTightEffSFHist = 0;
  //   TFile *eleEffSFFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/ScaleFactors/Run1/ElectronSelection_Run2012ReReco_53X.root","READ");
  //   eleLooseEffSFHist = (TH2D*)eleEffSFFile->Get("sfLOOSE");
  //   assert(eleLooseEffSFHist);
  //   eleTightEffSFHist = (TH2D*)eleEffSFFile->Get("sfTIGHT");
  //   assert(eleTightEffSFHist);

  if (outFileName.empty()){
    cout << "HZZRazor: Output filename not specified!" << endl << "Using default output name HZZRazor.root" << endl;
    outFileName = "HZZRazor.root";
  }
  TFile outFile(outFileName.c_str(), "RECREATE");
    
  //one tree to hold all events
  TTree *razorTree = new TTree("HZZRazor", "Info on selected razor inclusive events");
    
  //initialize jet energy corrections
  TRandom3 *random = new TRandom3(33333); //Artur wants this number 33333
  std::vector<JetCorrectorParameters> correctionParameters;
  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/RazorAnalyzer/data/JEC/";
  cout << "Getting JEC parameters from " << pathname << endl;  
  if (IsData) {
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_DATA_L1FastJet_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_DATA_L2Relative_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_DATA_L3Absolute_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_DATA_L2L3Residual_AK4PFchs.txt", pathname.c_str())));
  } else {
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_50nsV4_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));
  }
  
  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);
  JetCorrectorParameters *JetResolutionParameters = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",pathname.c_str()));
  SimpleJetResolution *JetResolutionCalculator = new SimpleJetResolution(*JetResolutionParameters);


  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

  //tree variables
  float weight = 1.0;
  unsigned int lumi, run, event;
  uint npu;
  int nLooseBTaggedJets, nMediumBTaggedJets;
  int nSelectedMuons, nSelectedElectrons;

  float MR;
  float Rsq;  
  float met;
  float metphi;
  float t1Rsq;  
  float t1met;
  float t1metphi;
  float HT;

  float m4l, pt4l, eta4l, phi4l;
  float lep1E, lep1Pt, lep1Eta, lep1Phi;
  float lep2E, lep2Pt, lep2Eta, lep2Phi;
  float lep3E, lep3Pt, lep3Eta, lep3Phi;
  float lep4E, lep4Pt, lep4Eta, lep4Phi;
  float mZ1, mZ2;

  //jet information
  int n_Jets = 0;
  float jet_E[10], jet_Pt[10], jet_Eta[10], jet_Phi[10], jet_CSV[10];
  bool jet_LooseBTag[10], jet_MediumBTag[10];
  float mHem1, ptHem1, etaHem1, phiHem1, mHem2, ptHem2, etaHem2, phiHem2;

  //float pileupWeight = 1.0;
  //float lepEffCorrFactor = 1.0;
  float btagCorrFactor = 1.0;

  //set branches on big tree
  razorTree->Branch("lumi", &lumi, "lumi/i");
  razorTree->Branch("run", &run, "run/i");
  razorTree->Branch("event", &event, "event/i");
  razorTree->Branch("npu", &npu, "npu/i");
  razorTree->Branch("nLooseBTaggedJets", &nLooseBTaggedJets, "nLooseBTaggedJets/I");
  razorTree->Branch("nMediumBTaggedJets", &nMediumBTaggedJets, "nMediumBTaggedJets/I");
  razorTree->Branch("nMuons", &nSelectedMuons, "nMuons/I");
  razorTree->Branch("nElectrons", &nSelectedElectrons, "nElectrons/I");
  razorTree->Branch("MR", &MR, "MR/F");
  razorTree->Branch("Rsq", &Rsq, "Rsq/F");
  razorTree->Branch("t1Rsq", &t1Rsq, "t1Rsq/F");
  razorTree->Branch("MET", &met, "MET/F");
  razorTree->Branch("METPhi", &metphi, "METPhi/F");
  razorTree->Branch("t1MET", &t1met, "t1MET/F");
  razorTree->Branch("t1METPhi", &t1metphi, "t1METPhi/F");
  razorTree->Branch("HT", &HT, "HT/F");

  razorTree->Branch("m4l", &m4l, "m4l/F");
  razorTree->Branch("pt4l", &pt4l, "pt4l/F");
  razorTree->Branch("eta4l", &eta4l, "eta4l/F");
  razorTree->Branch("phi4l", &phi4l, "phi4l/F");    
  razorTree->Branch("lep1E", &lep1E, "lep1E/F");
  razorTree->Branch("lep1Pt", &lep1Pt, "lep1Pt/F");
  razorTree->Branch("lep1Eta", &lep1Eta, "lep1Eta/F");
  razorTree->Branch("lep1Phi", &lep1Phi, "lep1Phi/F");
  razorTree->Branch("lep2E", &lep2E, "lep2E/F");
  razorTree->Branch("lep2Pt", &lep2Pt, "lep2Pt/F");
  razorTree->Branch("lep2Eta", &lep2Eta, "lep2Eta/F");
  razorTree->Branch("lep2Phi", &lep2Phi, "lep2Phi/F");
  razorTree->Branch("lep3E", &lep3E, "lep3E/F");
  razorTree->Branch("lep3Pt", &lep3Pt, "lep3Pt/F");
  razorTree->Branch("lep3Eta", &lep3Eta, "lep3Eta/F");
  razorTree->Branch("lep3Phi", &lep3Phi, "lep3Phi/F");
  razorTree->Branch("lep4E", &lep4E, "lep4E/F");
  razorTree->Branch("lep4Pt", &lep4Pt, "lep4Pt/F");
  razorTree->Branch("lep4Eta", &lep4Eta, "lep4Eta/F");
  razorTree->Branch("lep4Phi", &lep4Phi, "lep4Phi/F");

  razorTree->Branch("mZ1", &mZ1, "mZ1/F");
  razorTree->Branch("mZ2", &mZ2, "mZ2/F");
    
  razorTree->Branch("n_Jets", &n_Jets, "n_Jets/I");
  razorTree->Branch("jet_E", jet_E, "jet_E[n_Jets]/F");
  razorTree->Branch("jet_Pt", jet_Pt, "jet_Pt[n_Jets]/F");
  razorTree->Branch("jet_Eta", jet_Eta, "jet_Eta[n_Jets]/F");
  razorTree->Branch("jet_Phi", jet_Phi, "jet_Phi[n_Jets]/F");
  razorTree->Branch("jet_CSV", jet_CSV, "jet_CSVi[n_Jets]/F");
  razorTree->Branch("jet_LooseBTag", jet_LooseBTag, "jet_LooseBTag[n_Jets]/F");
  razorTree->Branch("jet_MediumBTag", jet_MediumBTag, "jet_MediumBTag[n_Jets]/F");
  razorTree->Branch("mHem1", &mHem1, "mHem1/F");
  razorTree->Branch("ptHem1", &ptHem1, "ptHem1/F");
  razorTree->Branch("etaHem1", &etaHem1, "etaHem1/F");
  razorTree->Branch("phiHem1", &phiHem1, "phiHem1/F");
  razorTree->Branch("mHem2", &mHem2, "mHem2/F");
  razorTree->Branch("ptHem2", &ptHem2, "ptHem2/F");
  razorTree->Branch("etaHem2", &etaHem2, "etaHem2/F");
  razorTree->Branch("phiHem2", &phiHem2, "phiHem2/F");

  //MET Filters
  // razorTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
  // razorTree->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, "Flag_CSCTightHaloFilter/O");
  // razorTree->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, "Flag_hcalLaserEventFilter/O");
  // razorTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/O");
  // razorTree->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/O");
  // razorTree->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, "Flag_trackingFailureFilter/O");
  // razorTree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/O");
  // razorTree->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, "Flag_ecalLaserCorrFilter/O");
  // razorTree->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters, "Flag_trkPOGFilters/O");
  // razorTree->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, "Flag_trkPOG_manystripclus53X/O");
  // razorTree->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, "Flag_trkPOG_toomanystripclus53X/O");
  // razorTree->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, "Flag_trkPOG_logErrorTooManyClusters/O");
  // razorTree->Branch("Flag_METFilters", &Flag_METFilters, "Flag_METFilters/O");
  razorTree->Branch("HLTDecision", &HLTDecision, "HLTDecision[100]/O");

  if (!IsData) {    
    razorTree->Branch("weight", &weight, "weight/F");
    //   razorTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
    //   razorTree->Branch("lepEffCorrFactor", &lepEffCorrFactor, "lepEffCorrFactor/F");
    //   razorTree->Branch("btagCorrFactor", &btagCorrFactor, "btagCorrFactor/F");
  } 
  //set branches on all trees
 

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
    NEvents->Fill(1.0);

    //reset tree variables
    run = runNum;
    lumi = lumiNum;
    event = eventNum;
    nLooseBTaggedJets = 0;
    nMediumBTaggedJets = 0;
    MR = -1;
    met = -1;
    Rsq = -1;
    t1met = -1;
    t1Rsq = -1;
    weight = 1.0;
    m4l = -999;
    pt4l = -999;
    eta4l = -999;
    phi4l = -999;
    lep1E = -999;
    lep1Pt = -999;
    lep1Eta = -999;
    lep1Phi = -999;
    lep2E = -999;
    lep2Pt = -999;
    lep2Eta = -999;
    lep2Phi = -999;
    lep3E = -999;
    lep3Pt = -999;
    lep3Eta = -999;
    lep3Phi = -999;
    lep4E =-999; 
    lep4Pt = -999;
    lep4Eta = -999;
    lep4Phi = -999;
    mZ1 = -999;
    mZ2 = -999;

    // //*****************************************
    // //TODO: triggers!
    // //*****************************************
    // bool passedDileptonTrigger = false;
    // bool passedSingleLeptonTrigger = false;
    // bool passedLeptonicTrigger = false;
    // bool passedHadronicTrigger= false;
    // if (isRunOne) {
    //   if (HLTDecision[46] || HLTDecision[47] ||HLTDecision[48] ||HLTDecision[49] ||HLTDecision[50]) passedHadronicTrigger = true;
    //   if (HLTDecision[3] || HLTDecision[4] || HLTDecision[6] || HLTDecision[7]|| HLTDecision[12]) passedDileptonTrigger = true;
    //   if (isData) {
    // 	if (HLTDecision[0] || HLTDecision[1] ||HLTDecision[8] ||HLTDecision[9]) passedSingleLeptonTrigger = true;
    //   } else {
    // 	if (HLTDecision[0] || HLTDecision[1] || HLTDecision[9]) passedSingleLeptonTrigger = true;
    //   }
    //   passedLeptonicTrigger = passedSingleLeptonTrigger || passedDileptonTrigger;
    //   if(!(passedLeptonicTrigger || passedHadronicTrigger)) continue; //ensure event passed a trigger
    // } else {
    //   passedDileptonTrigger = true;
    //   passedSingleLeptonTrigger = true;
    //   passedLeptonicTrigger = passedSingleLeptonTrigger || passedDileptonTrigger;
    //   passedHadronicTrigger = true;
    // }
        
    //*****************************************
    //Get Pileup Information
    //*****************************************
    npu = 0;
    for (int i=0; i < nBunchXing; ++i) {
      if (BunchXing[i] == 0) {
	npu = nPUmean[i];
      }
      if (BunchXing[i] == -1) {
	npu = nPUmean[i];
      }
      if (BunchXing[i] == 1) {
	npu = nPUmean[i];
      }	  
    }


    //*****************************************
    //Select Leptons
    //*****************************************
    // cout << "\n\n***********************************************\n";
    // cout << "Event: " << runNum << " " << eventNum << "\n";
    vector<TLorentzVector> preselLeptons;
    vector<TLorentzVector> GoodLeptons;
    vector<int> preselType;
    vector<int> preselIndex;
    vector<int> preselCharge;

    for(int i = 0; i < nMuons; i++){

      if(muonPt[i] < 5) continue;
      if(abs(muonEta[i]) > 2.4) continue;

      //Calculate MC->Data Scale Factors
      if (RazorAnalyzer::matchesGenMuon(muonEta[i],muonPhi[i])) {	
	//apply muon efficiency scale factors
      }

      //cout << "Muon " << i << " " << muonPt[i] << " " << muonEta[i] << " " << muonPhi[i] << " : " << passHZZMuonPreselection(i) << " " << isHZZMuon(i) << " | " << "\n";

      if(!passHZZMuonPreselection(i)) continue;  
      if (isHZZMuon(i)) nSelectedMuons++;

      preselType.push_back(13);
      preselIndex.push_back(i);
      preselCharge.push_back(muonCharge[i]);
      TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]);
      preselLeptons.push_back(thisMuon);
      if(isHZZMuon(i)) GoodLeptons.push_back(thisMuon);
    }
    for(int i = 0; i < nElectrons; i++){
      if(elePt[i] < 5) continue;
      if(fabs(eleEta[i]) > 2.5) continue;

      //remove overlaps
      bool overlap = false;
      for(auto& lep : preselLeptons){
	if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.Eta(),lep.Phi()) < 0.05) overlap = true;
      }
      //cout << "Electron " << i << " " << elePt[i] << " " << eleEta[i] << " " << elePhi[i] << " : overlap=" << overlap << " " << passHZZElectronPreselection(i) << " " << isHZZElectron(i) << "\n";

      if(overlap) continue;

      if(!passHZZElectronPreselection(i)) continue; 
      if (isHZZElectron(i)) nSelectedElectrons++;

      preselType.push_back(11);
      preselIndex.push_back(i);
      preselCharge.push_back(muonCharge[i]);
      TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
      preselLeptons.push_back(thisElectron); 
      if(isHZZElectron(i)) GoodLeptons.push_back(thisElectron);
    }

    
    //***********************************************
    //select Z1 candidate
    //***********************************************
    bool foundZ1Candidate = false;
    pair<int,int> Z1LeptonPairIndices;  Z1LeptonPairIndices.first = -1; Z1LeptonPairIndices.second = -1;
    double bestZ1Mass = 0;
    for (int i = 0 ; i < int(preselType.size()); ++i) {
      for (int j = i+1 ; j < int(preselType.size()); ++j) {

	//leptons must pass full selection criteria
	if (!( (preselType[i] == 13 && isHZZMuon(preselIndex[i]) )
	       || (preselType[i] == 11 && isHZZElectron(preselIndex[i]) )
	       )) continue;
	if (!( (preselType[j] == 13 && isHZZMuon(preselIndex[j]) )
	       || (preselType[j] == 11 && isHZZElectron(preselIndex[j]) )
	       )) continue;
	
	//look for opp sign, same flavor pairs
	if (preselCharge[i]*preselCharge[j] < 0 && preselType[i] == preselType[j]
	    && ( (preselLeptons[i]+preselLeptons[j]).M() > 40 && (preselLeptons[i]+preselLeptons[j]).M() < 120 ) 
	    && ( bestZ1Mass == 0 || fabs( (preselLeptons[i]+preselLeptons[j]).M() - 91.1876 ) < fabs( bestZ1Mass - 91.1876))
	    ) {
	  foundZ1Candidate = true;
	  bestZ1Mass =  (preselLeptons[i]+preselLeptons[j]).M();
	  Z1LeptonPairIndices.first = i;
	  Z1LeptonPairIndices.second = j;
	}       
      }
    }
    TLorentzVector Z1Candidate; Z1Candidate.SetPxPyPzE(0,0,0,0); 
    if (foundZ1Candidate) {
      Z1Candidate = preselLeptons[Z1LeptonPairIndices.first]+preselLeptons[Z1LeptonPairIndices.second];
      //cout << "Z1Candidate: " << Z1Candidate.M() << " " << Z1Candidate.Pt() << " " << Z1Candidate.Eta() << " " << Z1Candidate.Phi() << "\n";
    } else {
      continue;
    }
    
    //***********************************************
    //select Z2 candidate
    //***********************************************
    bool foundZ2Candidate = false;
    pair<int,int> Z2LeptonPairIndices;  Z2LeptonPairIndices.first = -1; Z2LeptonPairIndices.second = -1;
    double leptonPtSum = 0;

    for (int i = 0 ; i < int(preselType.size()); ++i) {
      for (int j = i+1 ; j < int(preselType.size()); ++j) {
	
	//leptons must pass full selection criteria
	if (!( (preselType[i] == 13 && isHZZMuon(preselIndex[i]) )
	       || (preselType[i] == 11 && isHZZElectron(preselIndex[i]) )
	       )) continue;
	if (!( (preselType[j] == 13 && isHZZMuon(preselIndex[j]) )
	       || (preselType[j] == 11 && isHZZElectron(preselIndex[j]) )
	       )) continue;
	
	//cannot be one of the Z1 leptons
	if (i == Z1LeptonPairIndices.first || i == Z1LeptonPairIndices.second) continue;
	if (j == Z1LeptonPairIndices.first || j == Z1LeptonPairIndices.second) continue;

	if (preselCharge[i]*preselCharge[j] < 0 && preselType[i] == preselType[j]
	    && ( (preselLeptons[i]+preselLeptons[j]).M() > 4 && (preselLeptons[i]+preselLeptons[j]).M() < 120 ) 
	    && ( leptonPtSum == 0 || (preselLeptons[i].Pt() + preselLeptons[j].Pt()) > leptonPtSum)
	    ) {
	  foundZ2Candidate = true;
	  leptonPtSum = preselLeptons[i].Pt() + preselLeptons[j].Pt();
	  Z2LeptonPairIndices.first = i;
	  Z2LeptonPairIndices.second = j;
	}
      }
    }
    TLorentzVector Z2Candidate; Z2Candidate.SetPxPyPzE(0,0,0,0); 
    if (foundZ2Candidate) {
      Z2Candidate = preselLeptons[Z2LeptonPairIndices.first]+preselLeptons[Z2LeptonPairIndices.second];
      //cout << "Z2Candidate: " << Z2Candidate.M() << " " << Z2Candidate.Pt() << " " << Z2Candidate.Eta() << " " << Z2Candidate.Phi() << "\n";
    } else {
      continue;
    }

    TLorentzVector ZZCandidate; ZZCandidate.SetPxPyPzE(0,0,0,0); 
    if (foundZ1Candidate && foundZ2Candidate) {
      ZZCandidate = Z1Candidate + Z2Candidate;      
      //cout << "ZZCandidate: " << ZZCandidate.M() << " " << ZZCandidate.Pt() << " " << ZZCandidate.Eta() << " " << ZZCandidate.Phi() << "\n";
    } else {
      continue;
    }
 
    //***********************************************
    //Get leading and subleading leptons
    //***********************************************
    double leadingLeptonPt = preselLeptons[Z1LeptonPairIndices.first].Pt();
    double subLeadingLeptonPt = preselLeptons[Z1LeptonPairIndices.second].Pt();
    if (leadingLeptonPt < subLeadingLeptonPt) {
      double tmp = leadingLeptonPt;
      leadingLeptonPt = subLeadingLeptonPt;
      subLeadingLeptonPt = tmp;
    }
    if ( preselLeptons[Z2LeptonPairIndices.first].Pt() > leadingLeptonPt ) {
      subLeadingLeptonPt = leadingLeptonPt;
      leadingLeptonPt = preselLeptons[Z2LeptonPairIndices.first].Pt();
    } else if ( preselLeptons[Z2LeptonPairIndices.first].Pt() > subLeadingLeptonPt) {
      subLeadingLeptonPt = preselLeptons[Z2LeptonPairIndices.first].Pt();
    }
    if ( preselLeptons[Z2LeptonPairIndices.second].Pt() > leadingLeptonPt ) {
      subLeadingLeptonPt = leadingLeptonPt;
      leadingLeptonPt = preselLeptons[Z2LeptonPairIndices.second].Pt();
    } else if ( preselLeptons[Z2LeptonPairIndices.second].Pt() > subLeadingLeptonPt) {
      subLeadingLeptonPt = preselLeptons[Z2LeptonPairIndices.second].Pt();
    }
     

    //***********************************************
    //Selection Cuts
    //***********************************************
    //Two leptons must have pt>20, 10
    if (!(leadingLeptonPt > 20 && subLeadingLeptonPt > 10)) continue;

    //Same Flavor, Opp charge pairs must not have mass < 4
    if ( preselType[Z2LeptonPairIndices.first] == preselType[Z1LeptonPairIndices.first] && 
	 preselCharge[Z2LeptonPairIndices.first]*preselCharge[Z1LeptonPairIndices.first] < 0 ) {
      if ( (preselLeptons[Z2LeptonPairIndices.first] + preselLeptons[Z1LeptonPairIndices.first]).M() < 4) continue;
    }
    if ( preselType[Z2LeptonPairIndices.second] == preselType[Z1LeptonPairIndices.first] && 
	 preselCharge[Z2LeptonPairIndices.second]*preselCharge[Z1LeptonPairIndices.first] < 0 ) {
      if ( (preselLeptons[Z2LeptonPairIndices.second] + preselLeptons[Z1LeptonPairIndices.first]).M() < 4) continue;
    }
    if ( preselType[Z2LeptonPairIndices.first] == preselType[Z1LeptonPairIndices.second] && 
	 preselCharge[Z2LeptonPairIndices.first]*preselCharge[Z1LeptonPairIndices.second] < 0 ) {
      if ( (preselLeptons[Z2LeptonPairIndices.first] + preselLeptons[Z1LeptonPairIndices.second]).M() < 4) continue;
    }
    if ( preselType[Z2LeptonPairIndices.second] == preselType[Z1LeptonPairIndices.second] && 
	 preselCharge[Z2LeptonPairIndices.second]*preselCharge[Z1LeptonPairIndices.second] < 0 ) {
      if ( (preselLeptons[Z2LeptonPairIndices.second] + preselLeptons[Z1LeptonPairIndices.second]).M() < 4) continue;
    }
    
    //Z2Mass > 12, Let's make this cut later in the ntuple
    //if (Z2Candidate.M() < 12) continue;

    //M_4l cuts, let's make a loose cut on this and cut tighter in the ntuple
    if ( ZZCandidate.M() < 50 ) continue;
    

    //***********************************************
    //Fill ZZ candidate information
    //***********************************************
    m4l = ZZCandidate.M();
    pt4l = ZZCandidate.Pt();
    eta4l = ZZCandidate.Eta();
    phi4l = ZZCandidate.Phi();
    lep1E = preselLeptons[Z1LeptonPairIndices.first].E();
    lep1Pt = preselLeptons[Z1LeptonPairIndices.first].Pt();
    lep1Eta = preselLeptons[Z1LeptonPairIndices.first].Eta();
    lep1Phi = preselLeptons[Z1LeptonPairIndices.first].Phi();
    lep2E = preselLeptons[Z1LeptonPairIndices.second].E();
    lep2Pt = preselLeptons[Z1LeptonPairIndices.second].Pt();
    lep2Eta = preselLeptons[Z1LeptonPairIndices.second].Eta();
    lep2Phi = preselLeptons[Z1LeptonPairIndices.second].Phi();
    lep3E = preselLeptons[Z2LeptonPairIndices.first].E();
    lep3Pt = preselLeptons[Z2LeptonPairIndices.first].Pt();
    lep3Eta = preselLeptons[Z2LeptonPairIndices.first].Eta();
    lep3Phi = preselLeptons[Z2LeptonPairIndices.first].Phi();
    lep4E = preselLeptons[Z2LeptonPairIndices.second].E();
    lep4Pt = preselLeptons[Z2LeptonPairIndices.second].Pt();
    lep4Eta = preselLeptons[Z2LeptonPairIndices.second].Eta();
    lep4Phi = preselLeptons[Z2LeptonPairIndices.second].Phi();
    mZ1 = Z1Candidate.M();
    mZ2 = Z2Candidate.M();


    vector<TLorentzVector> GoodJets;
    //initialize B-Tagging Correction Factor
    btagCorrFactor = 1.0;

    //***********************************************
    //Variables for Type1 Met Correction
    //***********************************************
    double MetX_Type1Corr = 0;
    double MetY_Type1Corr = 0;



    n_Jets = 0;
    for(int i=0; i<10; ++i) {
      jet_E[n_Jets] = -999;
      jet_Pt[n_Jets] = -999;
      jet_Eta[n_Jets] = -999;
      jet_Phi[n_Jets] = -999;
    }
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
      if(deltaR > 0 && deltaR < 0.2) continue; //jet matches a selected lepton
      

      //*****************************************************************
      //apply Jet ID
      //*****************************************************************
      if (!jetPassIDLoose[i]) continue;

      //*****************************************************************
      //Apply Jet Energy and Resolution Corrections
      //*****************************************************************
      double tmpRho = fixedGridRhoFastjetAll;
      double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
    					     tmpRho, jetJetArea[i], 
    					     JetCorrector);   

      double jetEnergySmearFactor = 1.0;
      if (!isData) {
	jetEnergySmearFactor = JetEnergySmearingFactor( jetPt[i]*JEC, jetEta[i], npu, JetResolutionCalculator, random);
      }
      
      TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);
      TLorentzVector UnCorrJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);      
      double jetCorrPt = jetPt[i]*JEC*jetEnergySmearFactor;

      //*******************************
      //B-Tagging Correction Factor
      //*******************************
      if (abs(jetPartonFlavor[i]) == 5 &&jetCorrPt > 20) {
	double tmpBTagCorrFactor = 1.0;
	
	double tmpCorrFactor = 0.938887 + 0.00017124 * jetCorrPt + (-2.76366e-07) * jetCorrPt * jetCorrPt ;
	double MCEff = 1.0;
	if (jetCorrPt < 50) MCEff = 0.65;
	else if (jetCorrPt < 80) MCEff = 0.70;
	else if (jetCorrPt < 120) MCEff = 0.73;
	else if (jetCorrPt < 210) MCEff = 0.73;
	else MCEff = 0.66;				 
	
	//if pass CSV Medium
	if( isCSVM(i)) {
	  tmpBTagCorrFactor = tmpCorrFactor;
	} else {
	  tmpBTagCorrFactor = ( 1/MCEff - tmpCorrFactor) / ( 1/MCEff - 1);
	}

	btagCorrFactor *= tmpBTagCorrFactor;
      }

      //*******************************
      //Add to Type1 Met Correction
      //*******************************
      if (jetPt[i]*JEC*jetEnergySmearFactor > 20) {
	MetX_Type1Corr += -1 * ( thisJet.Px() - UnCorrJet.Px()  );
	MetY_Type1Corr += -1 * ( thisJet.Py() - UnCorrJet.Py()  );
      }

      //*******************************************************
      //apply  Pileup Jet ID
      //*******************************************************
      bool passPUJetID = true;
      // int level = 2; //loose jet ID
      // if (!((jetPileupIdFlag[i] & (1 << level)) != 0)) passPUJetID = false;
      
      //*******************************************************
      //apply Jet cuts
      //*******************************************************
      if(jetCorrPt < 30) continue;
      if(fabs(jetEta[i]) > 3.0) continue;
            
      if (passPUJetID || !passPUJetID) {
	//cout << "Jet " << i << " : " << thisJet.Pt() << " " << thisJet.Eta() << " " << thisJet.Phi() << " | " << isOldCSVL(i) << " " << isOldCSVM(i) << " | " << passPUJetID << "\n";
      }

      GoodJets.push_back(thisJet);

      jet_E[n_Jets] = thisJet.E();
      jet_Pt[n_Jets] = thisJet.Pt();
      jet_Eta[n_Jets] = thisJet.Eta();
      jet_Phi[n_Jets] = thisJet.Phi();
      jet_CSV[n_Jets] = jetCISV[i];
      jet_LooseBTag[n_Jets] = isCSVL(i);
      jet_MediumBTag[n_Jets] = isCSVM(i);
      n_Jets++;

      if(isCSVL(i)){ 
    	nLooseBTaggedJets++;
      }
      if(isCSVM(i)){ 
    	nMediumBTaggedJets++;
      }
    }


    //Compute the razor variables using the selected jets and possibly leptons
    vector<TLorentzVector> GoodPFObjects;
    GoodPFObjects.push_back(ZZCandidate);
    for(auto& jet : GoodJets) {
      GoodPFObjects.push_back(jet);
    }


    //*************************************************************
    //Apply Type1 Met Correction
    //*************************************************************
    double PFMetX = metPt*cos(metPhi) + MetX_Type1Corr;
    double PFMetY = metPt*sin(metPhi) + MetY_Type1Corr;
    TLorentzVector PFMETCorr; 
    PFMETCorr.SetPxPyPzE(PFMetX, PFMetY, 0, sqrt(PFMetX*PFMetX + PFMetY*PFMetY));      
    TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
    TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM(metType0Plus1Pt, 0, metType0Plus1Phi, 0);

    HT = 0;
    for(auto& obj : GoodPFObjects) HT += obj.Pt();


    vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
    std::vector< std::vector<int> > hemisphereObjectIndices = getHemispheresV2( GoodPFObjects );
    MR = computeMR(hemispheres[0], hemispheres[1]); 
    met = metPt;
    metphi = metPhi;
    Rsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
    t1met = metType0Plus1Pt;
    t1metphi = metType0Plus1Phi;
    t1Rsq = computeRsq(hemispheres[0], hemispheres[1], t1PFMET);

    int HiggsHemIndex= -1;
    for( auto& tmp : hemisphereObjectIndices[0] ) {
      if ( tmp == 0 ) HiggsHemIndex = 0;
    }
    for( auto& tmp : hemisphereObjectIndices[1] ) {
      if ( tmp == 0 ) HiggsHemIndex = 1;
    }
    
    if( HiggsHemIndex == 0 ) {
      mHem1   = hemispheres[0].M();
      ptHem1  = hemispheres[0].Pt();
      etaHem1 = hemispheres[0].Eta();
      phiHem1 = hemispheres[0].Phi();
      mHem2   = hemispheres[1].M();
      ptHem2  = hemispheres[1].Pt();
      etaHem2 = hemispheres[1].Eta();
      phiHem2 = hemispheres[1].Phi();
    }
    else if( HiggsHemIndex == 1 ) {
      mHem1   = hemispheres[1].M();
      ptHem1  = hemispheres[1].Pt();
      etaHem1 = hemispheres[1].Eta();
      phiHem1 = hemispheres[1].Phi();
      mHem2   = hemispheres[0].M();
      ptHem2  = hemispheres[0].Pt();
      etaHem2 = hemispheres[0].Eta();
      phiHem2 = hemispheres[0].Phi();
    }	

    //**********************************************************************
    //Apply ECAL Dead Cells Filter
    //**********************************************************************
    // if (isRunOne) {
    //   if (isData) {
    // 	if (!Flag_HBHENoiseFilter || !Flag_CSCTightHaloFilter || !Flag_eeBadScFilter ) {
    // 	  cout << "Fail noise filter\n";
    // 	  continue;
    // 	}
    //   }
    // } else {
    //   if (Flag_EcalDeadCellTriggerPrimitiveFilter == false) continue;
    // }

    //**********************************************************************
    //Compute correction factor weight
    //**********************************************************************
    // weight *= pileupWeight;
    // weight *= lepEffCorrFactor;
    // weight *= btagCorrFactor;    


    //**********************************************************************
    //Fill Tree
    //**********************************************************************
    // cout << "MR = " << MR << " , R^2 = " << t1Rsq << " , HT = " << HT << " t1met = " << t1met << "\n";
    // cout << "Fill\n";
    razorTree->Fill();

    //******************************
    //Print Debug
    //******************************
    if (printdebug) {
      cout << "\nNew Event\n";
      for(int j = 0; j < nGenParticle; j++){
	cout << "GenParticle " << j << " : " << gParticleId[j] << " " << gParticleStatus[j] << " " << gParticlePt[j] << " " << gParticleEta[j] << " " << gParticlePhi[j] << " : " << gParticleMotherId[j] << "\n";
      }
    }



  }//end of event loop
  
  cout << "Writing output trees..." << endl;
  razorTree->Write();
  NEvents->Write();

  outFile.Close();
}
