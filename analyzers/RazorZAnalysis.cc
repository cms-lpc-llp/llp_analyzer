#include "RazorZAnalysis.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>

//ROOT includes
#include "TH1F.h"

using namespace std;

void RazorZAnalysis::Analyze(bool isData, int option, string outFileName, string label)
{
  //initialization: create one TTree for each analysis box 
  cout << "Initializing..." << endl;
  bool combineTrees = true;
  if (outFileName.empty()){
    cout << "RazorInclusive: Output filename not specified!" << endl << "Using default output name RazorInclusive.root" << endl;
    outFileName = "RazorInclusive.root";
  }
  TFile outFile(outFileName.c_str(), "RECREATE");
    
  //one tree to hold all events
  TTree *razorTree = new TTree("RazorInclusive", "Info on selected razor inclusive events");
    
  //initialize jet energy corrections
  std::vector<JetCorrectorParameters> correctionParameters;
  //get correct directory for JEC files (different for lxplus and t3-higgs)
  struct stat sb;
  string dir;
  if(stat("/afs/cern.ch/work/s/sixie/public", &sb) == 0 && S_ISDIR(sb.st_mode)){ //check if Si's directory exists
    dir = "/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data";
    cout << "Getting JEC parameters from " << dir << endl;
  }
  else{ //we are on t3-higgs (for running locally on your laptop we need a separate solution)
    dir = Form("%s/src/RazorAnalyzer/data/", getenv("CMSSW_BASE"));
    cout << "Getting JEC parameters from " << dir << endl;
  }
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/PHYS14_V2_MC_L1FastJet_AK4PFchs.txt", dir.c_str())));
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/PHYS14_V2_MC_L2Relative_AK4PFchs.txt", dir.c_str())));
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/PHYS14_V2_MC_L3Absolute_AK4PFchs.txt", dir.c_str())));    

  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);

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
  int nSelectedJets, nBTaggedJets;
  int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
  int nVetoMuons, nVetoElectrons, nLooseTaus;
  int nIsoChargedPFCandidate;
  int nJets40, nJets80;
  float dPhiRazor;
  float theMR;
  float theRsq;  
  float dileptonMass;
  float dileptonPt;
  float met;
  float HT;
  bool passTrigger;
  float npu;
  UInt_t run;
  UInt_t lumi;
  UInt_t event;
  RazorBox box;

  //set branches on big tree
  if(combineTrees){
    razorTree->Branch("npu", &npu, "npu/F");
    razorTree->Branch("run", &run, "run/i");
    razorTree->Branch("lumi", &lumi, "lumi/i");
    razorTree->Branch("event", &event, "event/i");
    razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
    razorTree->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
    razorTree->Branch("nVetoMuons", &nVetoMuons, "nVetoMuons/I");
    razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
    razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
    razorTree->Branch("nVetoElectrons", &nVetoElectrons, "nVetoElectrons/I");
    razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
    razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
    razorTree->Branch("nIsoChargedPFCandidate", &nIsoChargedPFCandidate, "nIsoChargedPFCandidate/I");
    razorTree->Branch("nJets40", &nJets40, "nJets40/I");
    razorTree->Branch("nJets80", &nJets80, "nJets80/I");
    razorTree->Branch("MR", &theMR, "MR/F");
    razorTree->Branch("dPhiRazor", &dPhiRazor, "dPhiRazor/F");
    razorTree->Branch("Rsq", &theRsq, "Rsq/F");
    razorTree->Branch("dileptonMass", &dileptonMass, "dileptonMass/F");
    razorTree->Branch("dileptonPt", &dileptonPt, "dileptonPt/F");
    razorTree->Branch("met", &met, "met/F");
    razorTree->Branch("HT", &HT, "HT/F");
    razorTree->Branch("passTrigger", &passTrigger, "passTrigger/O");
    razorTree->Branch("box", &box, "box/I");
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
      box.second->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
      box.second->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
      box.second->Branch("dPhiRazor", &dPhiRazor, "dPhiRazor/F");
      box.second->Branch("MR", &theMR, "MR/F");
      box.second->Branch("met", &met, "met/F");
      box.second->Branch("HT", &HT, "HT/F");
      box.second->Branch("Rsq", &theRsq, "Rsq/F");
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

    //reset tree variables
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
    nIsoChargedPFCandidate = 0;
    nJets40 = 0;
    nJets80 = 0;
    dPhiRazor = 0;
    theMR = -1;
    theRsq = -1;
    dileptonMass = 0;
    dileptonPt = 0;
    met = 0;
    HT = 0;
    passTrigger = false;
    if(combineTrees) box = NONE;

    
    //get NPU
    for (int i=0; i < nBunchXing; ++i) {
      if (BunchXing[i] == 0) {
	npu = nPUmean[i];
      }   
    }
    run = runNum;
    lumi = lumiNum;
    event = eventNum;
    

    //TODO: triggers!
    bool passedLeptonicTrigger = true;
    bool passedHadronicTrigger= true;
    if(!(passedLeptonicTrigger || passedHadronicTrigger)) continue; //ensure event passed a trigger
        
    vector<TLorentzVector> GoodLeptons; //leptons used to compute hemispheres
    vector<TLorentzVector> LooseLeptons; //leptons used to define the Z
    for(int i = 0; i < nMuons; i++){

      if(muonPt[i] < 5) continue;
      if(abs(muonEta[i]) > 2.4) continue;
 
     //remove overlaps
      bool overlap = false;
      for(auto& lep : GoodLeptons){
	if (RazorAnalyzer::deltaR(muonEta[i],muonPhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
      }
      if(overlap) continue;

      if(isVetoMuon(i)) nVetoMuons++;
      if(isLooseMuon(i) && muonPt[i] >= 10) nLooseMuons++;
      if(isTightMuon(i) && muonPt[i] >= 10) nTightMuons++;

      if(!isVetoMuon(i)) continue;  
      TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
      GoodLeptons.push_back(thisMuon);
      if (isLooseMuon(i) && muonPt[i] >= 10) LooseLeptons.push_back(thisMuon);
    }
    for(int i = 0; i < nElectrons; i++){
      if(elePt[i] < 5) continue;
      if(fabs(eleEta[i]) > 2.5) continue;
      if(isMVANonTrigVetoElectron(i)) nVetoElectrons++;
      if(isLooseElectron(i) && elePt[i] > 10 ) nLooseElectrons++;
      if(isTightElectron(i) && elePt[i] > 10 ) nTightElectrons++;

      //remove overlaps
      bool overlap = false;
      for(auto& lep : GoodLeptons){
	if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
      }
      if(overlap) continue;

      if(!isMVANonTrigVetoElectron(i)) continue; 
      TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
      GoodLeptons.push_back(thisElectron);  
      if (isLooseElectron(i) && elePt[i] > 10 ) LooseLeptons.push_back(thisElectron);
    }

    // for(uint i = 0; i < nIsoPFCandidates; i++){

    //   //remove overlaps
    //   bool overlap = false;
    //   for(auto& lep : GoodLeptons){
    // 	if (RazorAnalyzer::deltaR(isoPFCandidateEta[i],isoPFCandidatePhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
    //   }
    //   if(overlap) continue;
      
    //   if (isoPFCandidatePt[i] > 20 && fabs(isoPFCandidateEta[i]) < 2.4 
    // 	  && isoPFCandidateIso04[i] / isoPFCandidatePt[i] < 0.4) {
    // 	nIsoChargedPFCandidate++;
    //   }
    // }
    

      
    vector<TLorentzVector> GoodJets;
    int numJetsAbove80GeV = 0;



    //***********************************************
    //Select Jets
    //***********************************************
    for(int i = 0; i < nJets; i++){

      double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
    					     fixedGridRhoFastjetAll, jetJetArea[i], 
    					     JetCorrector);   
      double jetCorrPt = jetPt[i]*JEC;
      double jetCorrE = jetE[i]*JEC;

      if(jetCorrPt < 40) continue;
      if(fabs(jetEta[i]) > 3.0) continue;

      //apply jet iD
      int level = 2; //loose jet ID
      if (!jetPassIDTight[i]) continue;
      if (!((jetPileupIdFlag[i] & (1 << level)) != 0)) continue;      

      //exclude selected muons and electrons from the jet collection
      double deltaR = -1;
      TLorentzVector thisJet = makeTLorentzVector(jetCorrPt, jetEta[i], jetPhi[i], jetCorrE);
      for(auto& lep : GoodLeptons){
    	double thisDR = thisJet.DeltaR(lep);
    	if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      }
      if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton
            
      if(jetCorrPt > 80) numJetsAbove80GeV++;
      GoodJets.push_back(thisJet);
      nSelectedJets++;

      if(isCSVM(i)){ 
    	nBTaggedJets++;
      }
    }

    nJets80 = numJetsAbove80GeV;
    nJets40 = nSelectedJets;


    //**********************************************************************
    //Event Selection
    //**********************************************************************
    if(nJets40 < 1) continue; //event fails to have two 80 GeV jets
    if(!(nLooseElectrons >= 2 || nLooseMuons >= 2)) continue;


    HT = 0;
    theMR = 0;
    theRsq = 0;
    dPhiRazor = 0;
    met = metPt;

    
    //**********************************************************************
    //Compute dilepton kinematics
    //**********************************************************************
    int lep1Index = -1;
    int lep2Index = -1;
    for (int k=0; k < int(LooseLeptons.size()); ++k) {
      if (lep1Index < 0) {
	lep1Index = k;
      } else if ( LooseLeptons[k].Pt() > LooseLeptons[lep1Index].Pt()) {
	lep2Index = lep1Index;
	lep1Index = k;
      } else if ( lep2Index < 0 ) {
	lep2Index = k;
      } else if ( LooseLeptons[k].Pt() > LooseLeptons[lep2Index].Pt() ) {
	lep2Index = k;
      }
    }    
    if (lep1Index >= 0 && lep2Index >= 0) {
      dileptonMass = (LooseLeptons[lep1Index] + LooseLeptons[lep2Index]).M();
      dileptonPt = (LooseLeptons[lep1Index] + LooseLeptons[lep2Index]).Pt();
    } 
      
    // cuts on Z leptons
    if (!(lep1Index >= 0 && lep2Index >= 0)) continue;
    if (!(LooseLeptons[lep1Index].Pt() > 25 &&  LooseLeptons[lep2Index].Pt() > 25)) continue;
    if (!(dileptonMass > 20)) continue;

    //**********************************************************************
    //Compute razor stuff
    //**********************************************************************
    if (nJets40 >= 1) {
      //Compute the razor variables using the selected jets and possibly leptons
      vector<TLorentzVector> GoodPFObjects;
      for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
      if(passedLeptonicTrigger) for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);
      TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
      
      HT = 0;
      for(auto& obj : GoodPFObjects) HT += obj.Pt();
      
      
      vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
      theMR = computeMR(hemispheres[0], hemispheres[1]); 
      theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
      dPhiRazor = deltaPhi(hemispheres[0].Phi(),hemispheres[1].Phi());
      met = metPt;
    }


    //**********************************************************************
    //Apply ECAL Dead Cells Filter
    //**********************************************************************
    if (Flag_HBHENoiseFilter == false) continue;
    if (Flag_CSCTightHaloFilter == false) continue;
    if (Flag_hcalLaserEventFilter == false) continue;
    if (Flag_EcalDeadCellTriggerPrimitiveFilter == false) continue;
    if (Flag_trackingFailureFilter == false) continue;
    if (Flag_eeBadScFilter == false) continue;
    //if (Flag_trkPOGFilters == false) continue;


 
    // for(int i = 0; i < nMuons; i++){

    //   if(muonPt[i] < 5) continue;
    //   if(abs(muonEta[i]) > 2.4) continue;  
    //   if(!isVetoMuon(i)) continue;  
    //   //remove overlaps

    //   cout << "mu: " << muonPt[i] << " " << muonEta[i] << " " << muonPhi[i] << " : " << isVetoMuon(i) << " " << isLooseMuon(i) << " " << isTightMuon(i) <<  "\n";
    // }
    // for(int i = 0; i < nElectrons; i++){
    //   if(elePt[i] < 5) continue;
    //   if(fabs(eleEta[i]) > 2.5) continue;
    
    //   //remove overlaps

    //   if(!isMVANonTrigVetoElectron(i)) continue; 
    //   cout << "ele: " << elePt[i] << " " << eleEta[i] << " " << elePhi[i]  << " : " << isMVANonTrigVetoElectron(i) << " " << isLooseElectron(i) << " " << isTightElectron(i) << "\n";
    // }



    //Save Trigger Info
    passTrigger = ( HLTDecision[3] || HLTDecision[4] || HLTDecision[12] );


    //MuEle Box
    if(passedLeptonicTrigger && nTightElectrons > 0 && nLooseMuons > 0 ){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = MuEle;
	  razorTree->Fill();
	}
	else razorBoxes["MuEle"]->Fill();	
      }
    }
    //MuMu Box
    else if(passedLeptonicTrigger && nTightMuons > 0 && nLooseMuons > 1){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = MuMu;
	  razorTree->Fill();
	}
	else razorBoxes["MuMu"]->Fill();	
      }
    }
    //EleEle Box
    else if(passedLeptonicTrigger && nTightElectrons > 0 && nLooseElectrons > 1 ){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = EleEle;
	  razorTree->Fill();
	}
	else razorBoxes["EleEle"]->Fill();
      }
    }
    //MuMultiJet Box
    else if(passedLeptonicTrigger && nTightMuons > 0 && nSelectedJets > 3){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = MuMultiJet;
	  razorTree->Fill();
	}
	else razorBoxes["MuMultiJet"]->Fill();	
      }     
    }
    //MuJet Box
    else if(passedLeptonicTrigger && nTightMuons > 0 ){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = MuJet;
	  razorTree->Fill();
	}
	else razorBoxes["MuJet"]->Fill();
      }     
    }
    //EleMultiJet Box
    else if(passedLeptonicTrigger && nTightElectrons > 0 && nSelectedJets > 3 ){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = EleMultiJet;
	  razorTree->Fill();
	}
	else razorBoxes["EleMultiJet"]->Fill();	
      }     
    }
    //EleJet Box
    else if(passedLeptonicTrigger && nTightElectrons > 0 ){
      if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = EleJet;
	  razorTree->Fill();
	}
	else razorBoxes["EleJet"]->Fill();
      }     
    }
    //Soft Lepton + MultiJet Box
    else if(passedHadronicTrigger && nLooseTaus + nVetoElectrons + nVetoMuons > 0 && nSelectedJets > 3){
      if(passesHadronicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = LooseLeptonMultiJet;
	  razorTree->Fill();
	}
	else razorBoxes["LooseLeptonMultiJet"]->Fill();
      }     
    }
    //MultiJet Box
    else if(passedHadronicTrigger && nSelectedJets > 3){
      if(passesHadronicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = MultiJet;
	  razorTree->Fill();
	}
	else razorBoxes["MultiJet"]->Fill();
      }   
    //Loose Lepton + DiJet Box
    } else if(passedHadronicTrigger && nLooseTaus + nVetoElectrons + nVetoMuons > 0){
      if(passesHadronicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = LooseLeptonDiJet;
	  razorTree->Fill();
	}
	else razorBoxes["LooseLeptonDiJet"]->Fill();
      }     
    } else if (passedHadronicTrigger) {
      if(passesHadronicRazorBaseline(theMR, theRsq)){ 
	if(combineTrees){
	  box = DiJet;
	  razorTree->Fill();
	}
	else razorBoxes["DiJet"]->Fill();
      }     
    }
  }//end of event loop
  
  cout << "Writing output trees..." << endl;
  if(combineTrees) razorTree->Write();
  else for(auto& box : razorBoxes) box.second->Write();
  NEvents->Write();

  outFile.Close();
}
