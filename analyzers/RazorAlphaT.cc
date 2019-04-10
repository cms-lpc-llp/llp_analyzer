#include "RazorAlphaT.h"
#include "JetCorrectorParameters.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

void RazorAlphaT::Analyze(bool isData, int option, string outFileName, string label)
{
  //initialization: create one TTree for each analysis box 
  cout << "Initializing..." << endl;
  bool combineTrees = true;
  if ( outFileName.empty() )
    {
      cout << "RazorAlphaT: Output filename not specified!" << endl << "Using default output name RazorAlphaT.root" << endl;
      outFileName = "RazorAlphaT.root";
    }
  if (combineTrees) cout << "Using combineTrees" << endl;
  else cout << "Using razorBoxes" << endl;
  
  TFile outFile(outFileName.c_str(), "RECREATE");
  
  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  if (cmsswPath == NULL) {
    cout << "Warning: CMSSW_BASE not detected. Exiting..." << endl;
    return;
  }
  
  //*******************************************
  //Jet Energy Corrections
  //*******************************************
  string pathname;
  if (cmsswPath != NULL) pathname = string(cmsswPath) + "/src/RazorAnalyzer/data/JEC/";
  else {
        cout << "Warning: CMSSW_BASE not detected.  Looking for JEC parameters in data/JEC";
        pathname = "data/JEC/";
    }

  cout << "Getting JEC parameters from " << pathname << endl;
  std::vector<JetCorrectorParameters> correctionParameters;
  if (isData) {
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt", pathname.c_str())));
  }
  else {
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
    correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));
  }
  
  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector( correctionParameters );
  
  //one tree to hold all events
  TTree *razorTree = new TTree("RazorAlphaT", "Info on selected razor DM events");
  
  // tree to compute the cut efficiency
  //TTree *effTree = new TTree("CutEfficiency","Efficiencies of cuts on DM events"); 

  //separate trees for individual boxes
  map<string, TTree*> razorBoxes;
  vector<string> boxNames;
  boxNames.push_back("MuMu");
  boxNames.push_back("EleEle");
  boxNames.push_back("Mu");
  boxNames.push_back("Ele");
  boxNames.push_back("MultiJet");
  boxNames.push_back("TwoBJet");
  boxNames.push_back("OneBJet");
  for(size_t i = 0; i < boxNames.size(); i++){
    razorBoxes[boxNames[i]] = new TTree(boxNames[i].c_str(), boxNames[i].c_str());
  }
  
  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  
  //tree variables
  int nSelectedJets, nBTaggedJetsL, nBTaggedJetsM, nBTaggedJetsT;
  int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
  UInt_t run, lumi, event;
  float MuonPt[5], MuonEta[5], MuonPhi[5], MuonE[5];
  float ElePt[5], EleEta[5], ElePhi[5], EleE[5];
  float MR;
  float HT, MHT;
  float alphaT, dPhiMin;
  float Rsq, t1Rsq, RsqCorr, t1RsqCorr;
  float t1metPt, t1metPhi, metPtCorr, metPhiCorr, t1metPtCorr, t1metPhiCorr;
  RazorBox box;
  float JetPt_uncorr[30], JetEta_uncorr[30], JetPhi_uncorr[30], JetE_uncorr[30];
  float JetPt[30], JetEta[30], JetPhi[30], JetE[30];
  bool hasMatchingGenJet[30];
  int matchingGenJetIndex[30];
  float leadingJetPt, leadingJetEta, subLeadingJetPt;
  
  
    //set branches on big tree
  if(combineTrees)
    {
    razorTree->Branch("run",&run,"run/i");
    razorTree->Branch("lumi",&lumi,"lumi/i");
    razorTree->Branch("event",&event,"event/i");
   
      razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
      razorTree->Branch("nBTaggedJetsL", &nBTaggedJetsL, "nBTaggedJetsL/I");
      razorTree->Branch("nBTaggedJetsM", &nBTaggedJetsM, "nBTaggedJetsM/I");
      razorTree->Branch("nBTaggedJetsT", &nBTaggedJetsT, "nBTaggedJetsT/I");
      razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
      razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
      razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
      razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
      razorTree->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
      
      razorTree->Branch("JetE_uncorr", JetE_uncorr, "JetE_uncorr[nSelectedJets]/F");
      razorTree->Branch("JetPt_uncorr", JetPt_uncorr, "JetPt_uncorr[nSelectedJets]/F");
      razorTree->Branch("JetPhi_uncorr", JetPhi_uncorr, "JetPhi_uncorr[nSelectedJets]/F");
      razorTree->Branch("JetEta_uncorr", JetEta_uncorr, "JetEta_uncorr[nSelectedJets]/F");

      razorTree->Branch("JetE", JetE, "JetE[nSelectedJets]/F");
      razorTree->Branch("JetPt", JetPt, "JetPt[nSelectedJets]/F");
      razorTree->Branch("JetEta", JetEta, "JetEta[nSelectedJets]/F");
      razorTree->Branch("JetPhi", JetPhi, "JetPhi[nSelectedJets]/F");
      
      razorTree->Branch("alphaT", &alphaT, "alphaT/F");
      razorTree->Branch("dPhiMin", &dPhiMin, "dPhiMin/F");
      
      razorTree->Branch("MR", &MR, "MR/F");
      razorTree->Branch("Rsq", &Rsq, "Rsq/F");
      razorTree->Branch("t1Rsq", &t1Rsq, "t1Rsq/F");
      razorTree->Branch("metPt", &metPt, "metPt/F");
      razorTree->Branch("metPhi", &metPhi, "metPhi/F");
      razorTree->Branch("t1metPt", &t1metPt, "t1metPt/F");
      razorTree->Branch("t1metPhi", &t1metPhi, "t1metPhi/F");
      razorTree->Branch("HT", &HT, "HT/F");
      razorTree->Branch("MHT", &MHT, "MHT/F");
      
      razorTree->Branch("leadingJetPt",&leadingJetPt,"leadingJetPt/F");
      razorTree->Branch("leadingJetEta",&leadingJetEta,"leadingJetEta/F");
      razorTree->Branch("subLeadingJetPt",&subLeadingJetPt,"subLeadingJetPt/F");
      
      razorTree->Branch("RsqCorr", &RsqCorr, "RsqCorr/F");
      razorTree->Branch("t1RsqCorr", &t1RsqCorr, "t1RsqCorr/F");
      razorTree->Branch("metPtCorr", &metPtCorr, "metPtCorr/F");
      razorTree->Branch("metPhiCorr", &metPhiCorr, "metPhiCorr/F");
      razorTree->Branch("t1metPtCorr", &t1metPtCorr, "t1metPtCorr/F");
      razorTree->Branch("t1metPhiCorr", &t1metPhiCorr, "t1metPhiCorr/F");
	
      razorTree->Branch("box", &box, "box/I");
      
      razorTree->Branch("HLTDecision", HLTDecision, "HLTDecision[160]/O");
    }
    //set branches on all trees

  else{ 
    for(auto& box : razorBoxes){
      box.second->Branch("nBTaggedJetsL", &nBTaggedJetsL, "nBTaggedJetsL/I");
      box.second->Branch("nBTaggedJetsM", &nBTaggedJetsM, "nBTaggedJetsM/I");
      box.second->Branch("nBTaggedJetsT", &nBTaggedJetsT, "nBTaggedJetsT/I");
      box.second->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
      box.second->Branch("MuonPt", MuonPt,"MuonPt[nLooseMuons]/F");
      box.second->Branch("MuonEta", MuonEta,"MuonEta[nLooseMuons]/F");
      box.second->Branch("MuonPhi", MuonPhi,"MuonPhi[nLooseMuons]/F");
      box.second->Branch("MuonE", MuonE,"MuonE[nLooseMuons]/F");
      box.second->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
      box.second->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
      box.second->Branch("ElePt", ElePt,"ElePt[nLooseElectrons]/F");
      box.second->Branch("EleEta", EleEta,"EleEta[nLooseElectrons]/F");
      box.second->Branch("ElePhi", ElePhi,"ElePhi[nLooseElectrons]/F");
      box.second->Branch("EleE", EleE,"EleE[nLooseElectrons]/F");
      box.second->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
      box.second->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
      
      //ADDED LINES
      
      box.second->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
      box.second->Branch("JetE_uncorr", JetE_uncorr, "JetE_uncorr[nSelectedJets]/F");
      box.second->Branch("JetPt_uncorr", JetPt_uncorr, "JetPt_uncorr[nSelectedJets]/F");
      box.second->Branch("JetPhi_uncorr", JetPhi_uncorr, "JetPhi_uncorr[nSelectedJets]/F");
      box.second->Branch("JetEta_uncorr", JetEta_uncorr, "JetEta_uncorr[nSelectedJets]/F");
      
      box.second->Branch("alphaT", &alphaT, "alphaT/F");
      box.second->Branch("dPhiMin", &dPhiMin, "dPhiMin/F");
      box.second->Branch("MR", &MR, "MR/F");
      box.second->Branch("Rsq", &Rsq, "Rsq/F");
      box.second->Branch("t1Rsq", &t1Rsq, "t1Rsq/F");
      box.second->Branch("metPt", &metPt, "metPt/F");
      box.second->Branch("metPhi", &metPhi, "metPhi/F");
      box.second->Branch("t1metPt", &t1metPt, "t1metPt/F");
      box.second->Branch("t1metPhi", &t1metPhi, "t1metPhi/F");
      
      box.second->Branch("RsqCorr", &RsqCorr, "RsqCorr/F");
      box.second->Branch("t1RsqCorr", &t1RsqCorr, "t1RsqCorr/F");
      box.second->Branch("metPtCorr", &metPtCorr, "metPtCorr/F");
      box.second->Branch("metPhiCorr", &metPhiCorr, "metPhiCorr/F");
      box.second->Branch("t1metPtCorr", &t1metPtCorr, "t1metPtCorr/F");
      box.second->Branch("t1metPhiCorr", &t1metPhiCorr, "t1metPhiCorr/F");
      
      box.second->Branch("genMetPt", &genMetPt, "genMetPt/F");
      box.second->Branch("genMetPhi", &genMetPhi, "genMetPhi/F");
      box.second->Branch("metPt", &metPt, "metPt/F");
      box.second->Branch("metPhi", &metPhi, "metPhi/F");
      
      box.second->Branch("nGenJets",&nGenJets, "nGenJets/I");
      box.second->Branch("genJetE", genJetE, "genJetE[nGenJets]/F");
      box.second->Branch("genJetPt", genJetPt, "genJetPt[nGenJets]/F");
      box.second->Branch("genJetEta", genJetEta, "genJetEta[nGenJets]/F");
      box.second->Branch("genJetPhi", genJetPhi, "genJetPhi[nGenJets]/F");
      

      box.second->Branch("JetE", JetE, "JetE[nSelectedJets]/F");
      box.second->Branch("JetPt", JetPt, "JetPt[nSelectedJets]/F");
      box.second->Branch("JetEta", JetEta, "JetEta[nSelectedJets]/F");
      box.second->Branch("JetPhi", JetPhi, "JetPhi[nSelectedJets]/F");
      
      box.second->Branch("hasMatchingGenJet", hasMatchingGenJet, "hasMatchingGenJet[nSelectedJets]/O");
      box.second->Branch("matchingGenJetIndex", matchingGenJetIndex, "matchingGenJetIndex[nSelectedJets]/I");
      
      box.second->Branch("HLTDecision", HLTDecision, "HLTDecision[120]/O");
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

        //fill normalization histogram
        NEvents->Fill(1.0);
      
        // event info
        run = runNum;
        lumi = lumiNum;
        event = eventNum;

        //reset tree variables
        nSelectedJets = 0;
        nBTaggedJetsL = 0;
	    nBTaggedJetsM = 0;
	    nBTaggedJetsT =0;
        nLooseMuons = 0;
        nTightMuons = 0;
        nLooseElectrons = 0;
        nTightElectrons = 0;
        nTightTaus = 0;
        alphaT    = -1.0;
        dPhiMin   = -1.0;
        MR        = -1.0;
        Rsq       = -1.0;
	    t1Rsq     = -1.0;
	    RsqCorr   = -1.0;
	    t1RsqCorr = -1.0;
      leadingJetPt = -1;
      leadingJetEta = 0.0;
      subLeadingJetPt = -1;
	  
	for ( int j = 0; j < 30; j++ )
          {
	    JetE_uncorr[j]   = 0.0;
	    JetPt_uncorr[j]  = 0.0;
	    JetPhi_uncorr[j] = 0.0;
	    JetEta_uncorr[j] = 0.0;
            JetE[j ]  = 0.0;
	    JetPt[j]  = 0.0;
	    JetPhi[j] = 0.0;
	    JetEta[j] = 0.0;
          }
        if(combineTrees) box = NONE;

	//DEBUGGING GENJETS OUTPUT
	/*	if (jentry < 5){
	  cout << "nGenJets: " << nGenJets << endl;
	  for (int j = 0; j < nGenJets; j++){
	    cout << "genJetE: " << genJetE[j] << endl;
	    cout << "genJetPt: " << genJetPt[j] << endl;
	    cout << "genJetEta: " << genJetEta[j] << endl;
	    cout << "genJetPhi: " << genJetPhi[j] << endl;
	  }
	}
	*/
        //TODO: triggers!
        bool passedLeptonicTrigger = true;
        bool passedHadronicTrigger= true;
        if ( !(passedLeptonicTrigger || passedHadronicTrigger) ) continue; //ensure event passed a trigger
        
        vector<TLorentzVector> LooseMu; //Muons use to compute MET
        for(int i = 0; i < nMuons; i++){
	  if ( muonPt[i] < 10.0 || fabs(muonEta[i]) > 2.4 ) continue;  
	  if ( isLooseMuon(i) )
	    {
	      nLooseMuons++;
	      TLorentzVector thisMuon = 
		makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
	      LooseMu.push_back(thisMuon);
	    }
	  if( isTightMuon(i) ) nTightMuons++;
        }
	

	//------------------------------
	// Reco+ID Muons
	//------------------------------
	vector<TLorentzVector> LooseEle; //Electrons use to compute MET
        for(int i = 0; i < nElectrons; i++){
	  if ( elePt[i] < 10.0 || fabs(eleEta[i]) > 2.4 ) continue;
	  if( isLooseElectron(i) )
	    {
	      nLooseElectrons++;
	      TLorentzVector thisElectron = 
		makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
	      LooseEle.push_back(thisElectron);
	    }
	  if( isTightElectron(i) ) nTightElectrons++;
        }
	
	//----------------------
	//ID Tau and Veto
	//----------------------
	for(int i = 0; i < nTaus; i++){
	  if(!isTightTau(i)) continue; 
	  nTightTaus++;
	}
	if(nTightTaus > 0)continue;
	
        
	//------------------------
	//Jet ID + Correction
	//------------------------
	vector<TLorentzVector> GoodJets;
	vector<TLorentzVector> GoodJets_uncorr;
        int numJetsAbove80GeV = 0;
	for( int i = 0; i < nJets; i++ )
	  {
	    if(jetPt[i] < 40.0 || fabs(jetEta[i]) > 3.0) continue; 
	    
	    //ADDED LINES
	    //int level = 2; //3rd bit of jetPileupIdFlag
	    /*
	      In the case of v1p5 ntuples the jetPassIDLoose flag is not implemented,
	      please comment out the continue statement below
	    */
	    
	    if ( !jetPassIDLoose[i] ) continue;
	    // if ( !((jetPileupIdFlag[i] & (1 << level)) != 0) ) continue;	  
	    
	    double deltaR = -1.0;
	    TLorentzVector thisJet_uncorr = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);
	    //Obtain JEC
	    double JEC = JetEnergyCorrectionFactor( jetPt[i], JetEta[i], jetPhi[i], jetE[i],
						    fixedGridRhoAll, jetJetArea[i],
						    JetCorrector );
	    if (jentry == 1) cout << "FixedGridRhoAll: " << fixedGridRhoAll << endl << "jetJetArea[" << i << "]: " << jetJetArea[i] << endl; 
	    //Apply JEC
	    TLorentzVector thisJet = makeTLorentzVector(JEC*jetPt[i], jetEta[i], jetPhi[i], JEC*jetE[i]);
	    
	    //Removing Muons from the jet collection 
	    for ( auto& Mu : LooseMu)
	      {
		double thisDR = thisJet_uncorr.DeltaR(Mu);
		if ( deltaR < 0.0 || thisDR < deltaR ) deltaR = thisDR;
	      }
	    if ( deltaR > 0.0 && deltaR < 0.3 ) continue; //jet matches a Loose Muon
	    //Removing Electrons form the jet collection
	    deltaR = -1.0;
	    for ( auto& Ele : LooseEle )
	      {
		double thisDR = thisJet_uncorr.DeltaR(Ele);
		if(deltaR < 0.0 || thisDR < deltaR) deltaR = thisDR;
	      }
	    if ( deltaR > 0.0 && deltaR < 0.3 ) continue; //jet matches a Loose Electron
	    
	    if ( thisJet.Pt() > 80.0 ) numJetsAbove80GeV++;
	    GoodJets_uncorr.push_back(thisJet_uncorr);
	    
	    GoodJets.push_back(thisJet);
	    nSelectedJets++;
	    //b-tagging 
	    if(isCSVL(i)) nBTaggedJetsL++;
	    if(isCSVM(i)) nBTaggedJetsM++;
	    if(isCSVT(i)) nBTaggedJetsT++;
	  }
	
        if ( numJetsAbove80GeV < 2 ) continue; //event fails to have two 80 GeV jets
	
	int jIndex = 0;
	for (auto& tmpJet: GoodJets_uncorr)
	  {
	    JetPt_uncorr[jIndex] = tmpJet.Pt();
	    JetE_uncorr[jIndex] = tmpJet.E();
	    JetPhi_uncorr[jIndex] = tmpJet.Phi();
	    JetEta_uncorr[jIndex] = tmpJet.Eta();
	    jIndex=jIndex+1;
	  }  

	//Find whether there is a matching genJet for the reconstructed jet
	jIndex = 0;
	for ( auto& tmpJet : GoodJets )
	  { 
	    matchingGenJetIndex[jIndex] = findClosestGenJet(JetEta_uncorr[jIndex], JetPhi_uncorr[jIndex]);
	    if ( matchingGenJetIndex[jIndex] == -1 )
	      {
		hasMatchingGenJet[jIndex] = false;
	      }
	    else
	      {
		hasMatchingGenJet[jIndex] = true;
	      }
	    jIndex++;
	  }

	//output information for corrected jets
	jIndex = 0;
	for (auto& tmpJet: GoodJets)
	  {	    
	    JetPt[jIndex]  = tmpJet.Pt();
	    JetE[jIndex]   = tmpJet.E();
	    JetEta[jIndex] = tmpJet.Eta();
	    JetPhi[jIndex] = tmpJet.Phi();
	
        if (JetPt[jIndex]> leadingJetPt && JetPt[jIndex]<14000) 
        {
	    subLeadingJetPt = leadingJetPt;
	    leadingJetPt = JetPt[jIndex];
        leadingJetEta = JetEta[jIndex];
	    }
	    else if (JetPt[jIndex]>subLeadingJetPt && JetPt[jIndex]<14000) 
        {
	     subLeadingJetPt = JetPt[jIndex];
	    }
	    jIndex++;
	  }
    
    // Compute the variables alpha T and dPhiMin using the selected jets
    alphaT = GetAlphaT(GoodJets);
    dPhiMin = GetDPhiMin(GoodJets);  

	//Compute the razor variables using the selected jets and possibly leptons
	TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
	TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM( metType1Pt, 0, metType1Phi, 0 );
	t1metPt  = metType1Pt;
	t1metPhi = metType1Phi;
	
    vector<TLorentzVector> hemispheres = getHemispheres( GoodJets );
    
    MR    = computeMR(hemispheres[0], hemispheres[1]); 
	Rsq   = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	t1Rsq = computeRsq(hemispheres[0], hemispheres[1], t1PFMET);
   
   
    // Compute HT and MHT
    float MhtX = 0., MhtY = 0.;
    HT = 0.; 
    for (auto& obj : GoodJets) { HT += obj.Pt(); MhtX += obj.Px(); MhtY += obj.Py(); }

      TLorentzVector MyMHT;
      MyMHT.SetPxPyPzE(-MhtX, -MhtY, 0, sqrt(pow(MhtX,2) + pow(MhtY,2)));

      MHT = MyMHT.Pt();
    
	//MuMu Box
        if ( passedLeptonicTrigger && nLooseMuons > 1 && nLooseElectrons == 0 && nBTaggedJetsL == 0 )
	  {
	    TLorentzVector pfmet = makeTLorentzVector(metPt, 0, metPhi, 0);
	    TLorentzVector t1pfmet = makeTLorentzVector(metType1Pt, 0, metType1Phi, 0);
	    int it = 0;
	    for ( auto& Mu : LooseMu )
	      {
		pfmet.SetPx( pfmet.Px() + Mu.Px() );
		pfmet.SetPy( pfmet.Py() + Mu.Py() );
		t1pfmet.SetPx( t1pfmet.Px() + Mu.Px() );
                t1pfmet.SetPy( t1pfmet.Py() + Mu.Py() );
		MuonPt[it]  = Mu.Pt();
		MuonEta[it] = Mu.Eta();
		MuonPhi[it] = Mu.Phi();
		MuonE[it]   = Mu.E();
		it++;
	      }
	    
	    metPtCorr  = pfmet.Pt();
	    metPhiCorr = pfmet.Phi();
	    t1metPtCorr  = t1pfmet.Pt();
	    t1metPhiCorr = t1pfmet.Phi();
	    RsqCorr   = computeRsq(hemispheres[0], hemispheres[1], pfmet);
	    t1RsqCorr = computeRsq(hemispheres[0], hemispheres[1], t1pfmet);
	    if(combineTrees)
	      {
		box = MuMu;
		razorTree->Fill();
	      }
	    else razorBoxes["MuMu"]->Fill();
	  }
        //EleEle Box
        else if ( passedLeptonicTrigger && nLooseElectrons > 1 && nLooseMuons == 0 && nBTaggedJetsL == 0 )
	  {
	    TLorentzVector pfmet = makeTLorentzVector(metPt, 0, metPhi, 0);
	    TLorentzVector t1pfmet = makeTLorentzVector(metType1Pt, 0, metType1Phi, 0);
	    int it = 0;
	    for ( auto& Ele : LooseEle )
	      {
		pfmet.SetPx(pfmet.Px() + Ele.Px());
		pfmet.SetPy(pfmet.Py() + Ele.Py());
		t1pfmet.SetPx(t1pfmet.Px() + Ele.Px());
                t1pfmet.SetPy(t1pfmet.Py() + Ele.Py());
		ElePt[it] = Ele.Pt();
		EleEta[it] = Ele.Eta();
		ElePhi[it] = Ele.Phi();
		EleE[it] = Ele.E();
	      }
	    
	    metPtCorr  = pfmet.Pt();
            metPhiCorr = pfmet.Phi();
            t1metPtCorr  = t1pfmet.Pt();
            t1metPhiCorr = t1pfmet.Phi();
            RsqCorr   = computeRsq(hemispheres[0], hemispheres[1], pfmet);
            t1RsqCorr =computeRsq(hemispheres[0], hemispheres[1], t1pfmet);
	    if( combineTrees )
	      {
		box = EleEle;
		razorTree->Fill();
	      }
	    else razorBoxes["EleEle"]->Fill();
	  }
	//Mu Box
        else if ( passedLeptonicTrigger && nLooseMuons == 1 && nLooseElectrons == 0 && nBTaggedJetsL == 0 )
	  {
	    TLorentzVector pfmet = makeTLorentzVector(metPt, 0, metPhi, 0);
	    TLorentzVector t1pfmet = makeTLorentzVector(metType1Pt, 0, metType1Phi, 0);
	    int it = 0;
	    for ( auto& Mu : LooseMu )
	      {
		pfmet.SetPx( pfmet.Px() + Mu.Px() );
		pfmet.SetPy( pfmet.Py() + Mu.Py() );
		t1pfmet.SetPx( t1pfmet.Px() + Mu.Px() );
                t1pfmet.SetPy( t1pfmet.Py() + Mu.Py() );
		MuonPt[it] = Mu.Pt();
		MuonEta[it] = Mu.Eta();
		MuonPhi[it] = Mu.Phi();
		MuonE[it] = Mu.E();
		it++;
	      }
	    
	    metPtCorr  = pfmet.Pt();
            metPhiCorr = pfmet.Phi();
            t1metPtCorr  = t1pfmet.Pt();
            t1metPhiCorr = t1pfmet.Phi();
            RsqCorr   = computeRsq(hemispheres[0], hemispheres[1], pfmet);
            t1RsqCorr =computeRsq(hemispheres[0], hemispheres[1], t1pfmet);
	    
	    if(combineTrees){
	      box = MuJet;
	      razorTree->Fill();
	    }
	    else razorBoxes["Mu"]->Fill();
	  }
	//Ele Box
        else if( passedLeptonicTrigger && nLooseElectrons == 1 && nLooseMuons == 0 && nBTaggedJetsL == 0 )
	  {
	    TLorentzVector pfmet   = makeTLorentzVector(metPt, 0, metPhi, 0);
	    TLorentzVector t1pfmet = makeTLorentzVector(metType1Pt, 0, metType1Phi, 0);
	    int it = 0;
	    for(auto& Ele : LooseEle){
	      pfmet.SetPx( pfmet.Px() + Ele.Px() );
	      pfmet.SetPy( pfmet.Py() + Ele.Py() );
	      t1pfmet.SetPx( t1pfmet.Px() + Ele.Px() );
              t1pfmet.SetPy( t1pfmet.Py() + Ele.Py() );
	      ElePt[it] = Ele.Pt();
	      EleEta[it] = Ele.Eta();
	      ElePhi[it] = Ele.Phi();
	      EleE[it] = Ele.E();
	    }
	    
	    metPtCorr  = pfmet.Pt();
            metPhiCorr = pfmet.Phi();
            t1metPtCorr  = t1pfmet.Pt();
            t1metPhiCorr = t1pfmet.Phi();
            RsqCorr   = computeRsq(hemispheres[0], hemispheres[1], pfmet);
            t1RsqCorr =computeRsq(hemispheres[0], hemispheres[1], t1pfmet);
	    
	    if(combineTrees){
	      box = EleJet;
	      razorTree->Fill();
	    }
	    else razorBoxes["Ele"]->Fill();
	  }
        //MultiJet Box
        else if ( passedHadronicTrigger && nBTaggedJetsL == 0 && nLooseMuons == 0 && nLooseElectrons == 0 )
	  {
	    if(combineTrees)
	      {
		box = MultiJet;
		razorTree->Fill();
	      }
	    else { razorBoxes["MultiJet"]->Fill(); }
	  }
        //TwoBJet Box
        else if ( passedHadronicTrigger && nBTaggedJetsT > 1 && nLooseMuons == 0 && nLooseElectrons == 0 )
	  {
	    if ( combineTrees )
	      {
		box = TwoBJet;
		razorTree->Fill();
	      }
	    else razorBoxes["TwoBJet"]->Fill();
	  }
        //OneBJet Box
        else if ( passedHadronicTrigger && nBTaggedJetsT == 1 && nLooseMuons == 0 && nLooseElectrons == 0 )
	  {
	    if(combineTrees){
	      box = OneBJet;
	      razorTree->Fill();
	    }
	    else razorBoxes["OneBJet"]->Fill();
	  }       
    }//end of event loop
    
    cout << "Writing output trees..." << endl;
    if(combineTrees) razorTree->Write();
    else for(auto& box : razorBoxes) box.second->Write();
    NEvents->Write();
    
    outFile.Close();
}
