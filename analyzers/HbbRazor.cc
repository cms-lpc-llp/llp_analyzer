#include "HbbRazor.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>
#include <assert.h>
#include <iostream>
#include <iomanip>


//ROOT includes
#include "TH1F.h"
#include "TH2D.h"

using namespace std;

struct greater_than_pt{
    inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){
        return p1.Pt() > p2.Pt();
    }
};

void HbbRazor::Analyze(bool isData, int option, string outFileName, string label)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  bool combineTrees = true;
  bool isRunOne = false;
  bool printdebug = false;

  // //Pileup Weights
  // TFile *pileupWeightFile = 0;
  // TH1D *pileupWeightHist = 0;
  // if (isRunOne) {
  //   pileupWeightFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Run1PileupWeights.root", "READ");
  //   pileupWeightHist = (TH1D*)pileupWeightFile->Get("PUWeight_Run1");
  //   assert(pileupWeightHist);
  // }
  //debug
  if (outFileName.empty()){
    cout << "HbbRazor: Output filename not specified!" << endl << "Using default output name HbbRazor.root" << endl;
    outFileName = "HbbRazor.root";
  }
  TFile outFile(outFileName.c_str(), "RECREATE");
    
  //one tree to hold all events
  TTree *razorTree = new TTree("HbbRazor", "Info on selected razor inclusive events");
    
  //initialize jet energy corrections
  TRandom3 *random = new TRandom3(33333); //Artur wants this number 33333
  std::vector<JetCorrectorParameters> correctionParameters;
  //get correct directory for JEC files (different for lxplus and t3-higgs)
  struct stat sb;
  string dir;
  if(stat("/afs/cern.ch/work/s/sixie/public", &sb) == 0 && S_ISDIR(sb.st_mode)){ //check if Si's directory exists
    if (isRunOne) 
      {
	//dir = "/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data";
	dir = "data/JEC";
      } 
    else
      {
	//dir = "/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/data";
	dir = "data/JEC";
      }
    cout << "Getting JEC parameters from " << dir << endl;
  }
  else{ //we are on t3-higgs (for running locally on your laptop we need a separate solution)
    dir = Form("%s/src/RazorAnalyzer/data/", getenv("CMSSW_BASE"));
    cout << "Getting JEC parameters from " << dir << endl;
  }

  if (isRunOne) {
    if (isData) {
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Winter14_V8_DATA_L1FastJet_AK5PF.txt", dir.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Winter14_V8_DATA_L2Relative_AK5PF.txt", dir.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Winter14_V8_DATA_L3Absolute_AK5PF.txt", dir.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Winter14_V8_DATA_L2L3Residual_AK5PF.txt", dir.c_str())));
    } else {
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer13_MC_L1FastJet_AK5PF.txt", dir.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer13_V4_MC_L2Relative_AK5PF.txt", dir.c_str())));
      correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer13_V4_MC_L3Absolute_AK5PF.txt", dir.c_str())));
    }
  } else {
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt", dir.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt", dir.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt", dir.c_str())));
  }
  
  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);
  JetCorrectorParameters *JetResolutionParameters = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",dir.c_str()));
  SimpleJetResolution *JetResolutionCalculator = new SimpleJetResolution(*JetResolutionParameters);


  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

  //tree variables
  int nBTaggedJets;
  int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons;
  int nLooseTaus;
  float dPhiRazor;
  float theMR;
  float theRsq;  
  float met;
  float HT;
  float weight = 1.0;
  // float pileupWeight = 1.0;
  float lepEffCorrFactor = 1.0;
  //float lepTrigCorrFactor = 1.0;
  float btagCorrFactor = 1.0;
  float b1pt = 0;
  float b1eta = 0;
  float b1phi = 0;
  float b2pt = 0;
  float b2eta = 0;
  float b2phi = 0;
  float mbb = 0;
  float ptbb = 0;
  float etabb = 0;
  //jet information
  int n_Jets = 0;
  float jet_E[10], jet_Pt[10], jet_Eta[10], jet_Phi[10];

  bool  passed_DiPFJet80_DiPFJet30_BTagCSVd07d05, passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05, passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03, passed_DiJet80Eta2p6_BTagIP3DFastPVLoose, passed_QuadJet45, passed_QuadJet50, passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200, passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;
  float theRsq_t0, met_t0, theRsq_t1, met_t1, theRsq_t01, met_t01;
  int nVtx, nPU_mean;


  //set branches on big tree
  if(combineTrees){
    razorTree->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
    razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
    razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
    razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
    razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
    razorTree->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
    razorTree->Branch("MR", &theMR, "MR/F");
    razorTree->Branch("dPhiRazor", &dPhiRazor, "dPhiRazor/F");
    razorTree->Branch("Rsq", &theRsq, "Rsq/F");
    razorTree->Branch("met", &met, "met/F");
    razorTree->Branch("HT", &HT, "HT/F");
    razorTree->Branch("mbb", &mbb, "mbb/F");
    razorTree->Branch("ptbb", &ptbb, "ptbb/F");
    razorTree->Branch("etabb", &etabb, "etabb/F");
    razorTree->Branch("b1pt", &b1pt, "b1pt/F");
    razorTree->Branch("b1eta", &b1eta, "b1eta/F");
    razorTree->Branch("b1phi", &b1phi, "b1phi/F");
    razorTree->Branch("b2pt", &b2pt, "b2pt/F");
    razorTree->Branch("b2eta", &b2eta, "b2eta/F");
    razorTree->Branch("b2phi", &b2phi, "b2phi/F");
    razorTree->Branch("b2phi", &b2phi, "b2phi/F");
    razorTree->Branch("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05", &passed_DiPFJet80_DiPFJet30_BTagCSVd07d05, "passed_DiPFJet80_DiPFJet30_BTagCSVd07d05/O");
    razorTree->Branch("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05", &passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05, "passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05/O");
    razorTree->Branch("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03", &passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03, "passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03/O");
    razorTree->Branch("passed_DiJet80Eta2p6_BTagIP3DFastPVLoose", &passed_DiJet80Eta2p6_BTagIP3DFastPVLoose, "passed_DiJet80Eta2p6_BTagIP3DFastPVLoose/O");
    razorTree->Branch("passed_QuadJet45", &passed_QuadJet45, "passed_QuadJet45/O");
    razorTree->Branch("passed_QuadJet50", &passed_QuadJet50, "passed_QuadJet50/O");
    razorTree->Branch("passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200", &passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200, "passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200/O");
    razorTree->Branch("passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200, "passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200/O");
    razorTree->Branch("n_Jets", &n_Jets, "n_Jets/I");
    razorTree->Branch("jet_E", jet_E, "jet_E[n_Jets]/F");
    razorTree->Branch("jet_Pt", jet_Pt, "jet_Pt[n_Jets]/F");
    razorTree->Branch("jet_Eta", jet_Eta, "jet_Eta[n_Jets]/F");
    razorTree->Branch("jet_Phi", jet_Phi, "jet_Phi[n_Jets]/F");
    
    razorTree->Branch("nVtx", &nVtx, "nVtx/I");
    razorTree->Branch("nPU_mean", &nPU_mean, "nPU_mean/I");
    razorTree->Branch("theRsq_t0",  &theRsq_t0, "theRsq_t0/F");
    razorTree->Branch("theRsq_t1",  &theRsq_t1, "theRsq_t1/F");
    razorTree->Branch("theRsq_t01", &theRsq_t01,"theRsq_t01/F");
    razorTree->Branch("met_t0",  &met_t0, "met_t0/F");
    razorTree->Branch("met_t1",  &met_t1, "met_t1/F");
    razorTree->Branch("met_t01", &met_t01,"met_t01/F");
    
    if (!isData) {    
      razorTree->Branch("weight", &weight, "weight/F");
    } else {
      razorTree->Branch("run", &runNum, "run/i");
      razorTree->Branch("lumi", &lumiNum, "lumi/i");
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
    NEvents->Fill(1.0);

    //reset tree variables
    nBTaggedJets = 0;
    nLooseMuons = 0;
    nTightMuons = 0;
    nLooseElectrons = 0;
    nTightElectrons = 0;
    nLooseTaus = 0;
    theMR = -1;
    theRsq = -1;
    weight = 1.0;
    passed_DiPFJet80_DiPFJet30_BTagCSVd07d05 = false;
    passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05 = false;
    passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03 = false;
    passed_DiJet80Eta2p6_BTagIP3DFastPVLoose = false;
    passed_QuadJet45 = false;
    passed_QuadJet50 = false;
    passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200 = false;
    passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200 = false;
    theRsq_t0  = 0;
    theRsq_t1  = 0;
    theRsq_t01 = 0;
    met_t0  = 0;
    met_t1  = 0;
    met_t01 = 0;
    nVtx = 0;
    nPU_mean = 0;
    
    //*****************************************
    //TODO: triggers!
    //*****************************************

    // apply dijet triggers
    if(HLTDecision[91] == 1 ) passed_DiPFJet80_DiPFJet30_BTagCSVd07d05 = true;
    if(HLTDecision[92] == 1 ) passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05 = true;
    if(HLTDecision[93] == 1 ) passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03 = true;
    if(HLTDecision[94] == 1 ) passed_DiJet80Eta2p6_BTagIP3DFastPVLoose = true;
    if(HLTDecision[95] == 1 ) passed_QuadJet45 = true;
    if(HLTDecision[96] == 1 ) passed_QuadJet50 = true;
    if(HLTDecision[175] == 1 ) passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200 = true;
    if(HLTDecision[176] == 1 ) passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200 = true;

    // PU information
    nVtx = nPV;
    if(!isData)
      for(int i=0; i<nBunchXing; i++)
	if(BunchXing[i]==0) nPU_mean = nPUmean[i];
    

    vector<TLorentzVector> GoodJets;
    int numJetsAbove80GeV = 0;

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

    vector<TLorentzVector> BJetCandidates; //candidates for H->bb

    //***********************************************
    //Select Jets
    //***********************************************
    for(int i = 0; i < nJets; i++){      

      //*****************************************************************
      //apply Jet ID
      //*****************************************************************
      if (isRunOne) {
	if (!jetPassIDTight[i]) continue;
      }

      //*****************************************************************
      //Apply Jet Energy and Resolution Corrections
      //*****************************************************************
      double tmpRho = fixedGridRhoFastjetAll;
      if (isRunOne) tmpRho = fixedGridRhoAll;
      double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
    					     tmpRho, jetJetArea[i], 
    					     JetCorrector);   

      double jetEnergySmearFactor = 1.0;
      if (!isData) {
	if (isRunOne) {
	  jetEnergySmearFactor = JetEnergySmearingFactor( jetPt[i]*JEC, jetEta[i], nPU_mean, JetResolutionCalculator, random);
	}
      }
      
      TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);
      TLorentzVector UnCorrJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);      
      double jetCorrPt = jetPt[i]*JEC*jetEnergySmearFactor;
      //double jetCorrE = jetE[i]*JEC*jetEnergySmearFactor;

      // cout << "Jet " << i << " : " << jetCorrPt << " " << jetEta[i] << " " << jetPhi[i] << " : " <<jetPartonFlavor[i] 
      // 	   << "\n";
      
      //*******************************
      //B-Tagging Correction Factor
      //*******************************
      if (isRunOne) {
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
	  if(isCSVM(i)) {
	    tmpBTagCorrFactor = tmpCorrFactor;
	  } else {
	    tmpBTagCorrFactor = ( 1/MCEff - tmpCorrFactor) / ( 1/MCEff - 1);
	  }

	  btagCorrFactor *= tmpBTagCorrFactor;
	}
      }

      //*******************************
      //Add to Type1 Met Correction
      //*******************************
      if (isRunOne) {
	if (jetPt[i]*JEC*jetEnergySmearFactor > 20) {
	  MetX_Type1Corr += -1 * ( thisJet.Px() - UnCorrJet.Px()  );
	  MetY_Type1Corr += -1 * ( thisJet.Py() - UnCorrJet.Py()  );
	}
      }

      //*******************************************************
      //apply  Pileup Jet ID
      //*******************************************************
      if (isRunOne) {
	int level = 2; //loose jet ID
	if (!((jetPileupIdFlag[i] & (1 << level)) != 0)) continue;
      }
      
      //*******************************************************
      //apply Jet cuts
      //*******************************************************
      if(jetCorrPt < 40) continue;
      if(fabs(jetEta[i]) > 2.5) continue;
            
      if(jetCorrPt > 80) numJetsAbove80GeV++;
      GoodJets.push_back(thisJet);

      if(isCSVT(i)){ 
    	nBTaggedJets++;
	TLorentzVector thisBJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);       
	BJetCandidates.push_back(thisBJet);
      }
    }

    //*************************************************************
    //find Higgs->bb candidates
    //*************************************************************
    // for (int i=0; i<BJetCandidates.size(); ++i) {
    //   cout << "BJet " << i << " : " << BJetCandidates[i].Pt() << " " << BJetCandidates[i].Eta()  << " " << BJetCandidates[i].Phi() 
    // 	   << "\n";
    // }

    b1pt = 0;
    b1eta = 0;
    b1phi = 0;
    b2pt = 0;
    b2eta = 0;
    b2phi = 0;
    mbb = 0;
    ptbb = 0;
    etabb = 0;
    bool foundHiggsCandidate = false;
    pair<TLorentzVector,TLorentzVector> HiggsToBBCandidate;
    double highestHiggsPt = 0;
    for (int i=0; i<int(BJetCandidates.size()); ++i) {
      for (int j=i+1; j< int (BJetCandidates.size()); ++j) {
	//cout << "BJet Pair " << i << " " << j << " : " << (BJetCandidates[i]+BJetCandidates[j]).Pt() << "\n";
	if ( (BJetCandidates[i]+BJetCandidates[j]).Pt() > highestHiggsPt) {
	  //cout << "highest pt\n";
	  foundHiggsCandidate = true;
	  highestHiggsPt = (BJetCandidates[i]+BJetCandidates[j]).Pt();
	  HiggsToBBCandidate.first = BJetCandidates[i];
	  HiggsToBBCandidate.second = BJetCandidates[j];
	}
      }
    }

    TLorentzVector HiggsCandidate;
    if (foundHiggsCandidate) {
      HiggsCandidate = HiggsToBBCandidate.first + HiggsToBBCandidate.second;
      mbb = HiggsCandidate.M();
      ptbb = HiggsCandidate.Pt();
      etabb = HiggsCandidate.Eta();
      b1pt = HiggsToBBCandidate.first.Pt();
      b1eta = HiggsToBBCandidate.first.Eta();
      b1phi = HiggsToBBCandidate.first.Phi();
      b2pt = HiggsToBBCandidate.second.Pt();
      b2eta = HiggsToBBCandidate.second.Eta();
      b2phi = HiggsToBBCandidate.second.Phi();
    }

    //if(numJetsAbove80GeV < 2) continue; //event fails to have two 80 GeV jets

    //Compute the razor variables using the selected jets and possibly leptons
    vector<TLorentzVector> GoodPFObjects;
    if (foundHiggsCandidate) {
      GoodPFObjects.push_back(HiggsCandidate);
    }
    for(auto& jet : GoodJets) {
      if (foundHiggsCandidate) {
       	if ( (jet.Pt() == HiggsToBBCandidate.first.Pt() && jet.Eta() == HiggsToBBCandidate.first.Eta() && jet.Phi() == HiggsToBBCandidate.first.Phi())
	     || 
	     (jet.Pt() == HiggsToBBCandidate.second.Pt() && jet.Eta() == HiggsToBBCandidate.second.Eta() && jet.Phi() == HiggsToBBCandidate.second.Phi())
	     ) {
	  continue;
	}
      }
      GoodPFObjects.push_back(jet);
    }

    //*************************************************************
    //Apply Type1 Met Correction
    //*************************************************************
    double PFMetX = metPt*cos(metPhi) + MetX_Type1Corr;
    double PFMetY = metPt*sin(metPhi) + MetY_Type1Corr;
    TLorentzVector PFMET; 
    if (isRunOne) {
      PFMET.SetPxPyPzE(PFMetX, PFMetY, 0, sqrt(PFMetX*PFMetX + PFMetY*PFMetY));      
    } else {
      PFMET.SetPxPyPzE(PFMetX, PFMetY, 0, sqrt(PFMetX*PFMetX + PFMetY*PFMetY));      
    }
    TLorentzVector PFMETUnCorr = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);

    HT = 0;
    for(auto& obj : GoodPFObjects) HT += obj.Pt();


    vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
    theMR = computeMR(hemispheres[0], hemispheres[1]); 
    theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
    dPhiRazor = deltaPhi(hemispheres[0].Phi(),hemispheres[1].Phi());
    met = metPt;

    TLorentzVector PFMET_t0; 
    TLorentzVector PFMET_t1; 
    TLorentzVector PFMET_t01; 
    
    met_t0 = metType0Pt;
    met_t1 = metType1Pt;
    met_t01 = metType0Plus1Pt;

    PFMET_t0.SetPxPyPzE(metType0Pt*cos(metType0Phi), metType0Pt*sin(metType0Phi), 0, metType0Pt);
    PFMET_t1.SetPxPyPzE(metType1Pt*cos(metType1Phi), metType1Pt*sin(metType1Phi), 0, metType1Pt);
    PFMET_t01.SetPxPyPzE(metType0Plus1Pt*cos(metType0Plus1Phi), metType0Plus1Pt*sin(metType0Plus1Phi), 0, metType0Plus1Pt);
    
    theRsq_t0  = computeRsq(hemispheres[0], hemispheres[1], PFMET_t0);
    theRsq_t1  = computeRsq(hemispheres[0], hemispheres[1], PFMET_t0);
    theRsq_t01 = computeRsq(hemispheres[0], hemispheres[1], PFMET_t01);

    sort(GoodJets.begin(), GoodJets.end(), greater_than_pt());

    for(int i=0; i<int(GoodJets.size()); i++)
      {
	jet_E[n_Jets] = GoodJets[i].E();
	jet_Pt[n_Jets] = GoodJets[i].Pt();
	jet_Eta[n_Jets] = GoodJets[i].Eta();
	jet_Phi[n_Jets] = GoodJets[i].Phi();
	n_Jets++;	
      }
      
    //**********************************************************************
    //Apply ECAL Dead Cells Filter
    //**********************************************************************
    if (isRunOne) {
      if (isData) {
	if (!Flag_HBHENoiseFilter || !Flag_CSCTightHaloFilter || !Flag_eeBadScFilter ) {
	  cout << "Fail noise filter\n";
	  continue;
	}
      }
    } else {
      if (Flag_EcalDeadCellTriggerPrimitiveFilter == false) continue;
    }

    //**********************************************************************
    //Compute correction factor weight
    //**********************************************************************
    if (isRunOne) {
      // weight *= pileupWeight;
      weight *= lepEffCorrFactor;
      weight *= btagCorrFactor;    
    }


    //**********************************************************************
    //Categorize Events into Razor Boxes 
    //**********************************************************************
    if (foundHiggsCandidate 
	&& mbb > 0 && theMR > 0 && theRsq > 0.0
	) {
      razorTree->Fill();
    }


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
  if(combineTrees) razorTree->Write();
  NEvents->Write();

  outFile.Close();
}
