
#include "VetoLeptonEfficiencyControlRegion.h"
#include "JetCorrectorParameters.h"
#include "ControlSampleEvents.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

//**********************************************************************************
//Do Event Selection type selection : Dilepton selection
//**********************************************************************************
void VetoLeptonEfficiencyControlRegion::Analyze(bool isData, int option, string outputfilename, string label)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    std::vector<JetCorrectorParameters> correctionParameters;
    correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data/PHYS14_V2_MC_L1FastJet_AK4PFchs.txt"));
    correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data/PHYS14_V2_MC_L2Relative_AK4PFchs.txt"));
    correctionParameters.push_back(JetCorrectorParameters("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/data/PHYS14_V2_MC_L3Absolute_AK4PFchs.txt"));    
    FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);

    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "RazorVetoLeptonStudy.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    ControlSampleEvents *events = new ControlSampleEvents;
    events->CreateTree(ControlSampleEvents::kTreeType_Default);
    events->tree_->SetAutoFlush(0);

    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  
    //begin loop
    if (fChain == 0) return;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {

      //begin event
      if(jentry % 1000 == 0) cout << "Processing entry " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

        //fill normalization histogram
        NEvents->Fill(1.0);
     
	//event info
	events->weight = 1.0;
	events->run = runNum;
	events->lumi = lumiNum;
	events->event = eventNum;
	events->processID = 0;

 
	//get NPU
	for (int i=0; i < nBunchXing; ++i) {
	  if (BunchXing[i] == 0) {
	    events->NPU_0 = nPU[i];
	  }
	  if (BunchXing[i] == -1) {
	    events->NPU_Minus1 = nPU[i];
	  }
	  if (BunchXing[i] == 1) {
	    events->NPU_Plus1 = nPU[i];
	  }	  
	}

 
        //TODO: triggers!
        bool passedLeptonicTrigger = true;
        bool passedHadronicTrigger= true;
        if(!(passedLeptonicTrigger || passedHadronicTrigger)) continue; //ensure event passed a trigger
        

	//******************************************
	//find generated leptons
	//******************************************
	vector<int> genLeptonIndex;

	for(int j = 0; j < nGenParticle; j++){
	  //look for electrons
	  if (abs(gParticleId[j]) == 11 && gParticleStatus[j] == 1 	      
	      && abs(gParticleEta[j]) < 2.5 && gParticlePt[j] > 5
	      ) {
	    if ( abs(gParticleMotherId[j]) == 24 
		 || abs(gParticleMotherId[j]) == 23 
		 || ( (abs(gParticleMotherId[j]) == 15 || abs(gParticleMotherId[j]) == 13 || abs(gParticleMotherId[j]) == 11 )
		       && gParticleMotherIndex[j] >= 0 
		      && (abs(gParticleMotherId[gParticleMotherIndex[j]]) == 24 || 
			  abs(gParticleMotherId[gParticleMotherIndex[j]]) == 23)
		     )
		 )  {	      
	      genLeptonIndex.push_back(j);	      
	    }	    
	  }

	  //look for muons
	  if (abs(gParticleId[j]) == 13 && gParticleStatus[j] == 1  
	      && abs(gParticleEta[j]) < 2.4 && gParticlePt[j] > 5
	      ) {
	    if ( abs(gParticleMotherId[j]) == 24 
		 || abs(gParticleMotherId[j]) == 23 
		 || ( (abs(gParticleMotherId[j]) == 15 || abs(gParticleMotherId[j]) == 13 || abs(gParticleMotherId[j]) == 11 )
		      && gParticleMotherIndex[j] >= 0 
		      && (abs(gParticleMotherId[gParticleMotherIndex[j]]) == 24 || 
			  abs(gParticleMotherId[gParticleMotherIndex[j]]) == 23)
		      )
		 ) {	     
	      genLeptonIndex.push_back(j);
	    }	    
	  }

	  //look for hadronic taus
	  if (abs(gParticleId[j]) == 15 && gParticleStatus[j] == 2 
	      && (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23)
	      && abs(gParticleEta[j]) < 2.4 && gParticlePt[j] > 20
	      ) {
	    bool isLeptonicTau = false;
	    for(int k = 0; k < nGenParticle; k++){
	      if ( (abs(gParticleId[k]) == 11 || abs(gParticleId[k]) == 13) && gParticleMotherIndex[k] == j) {
		isLeptonicTau = true;
		break;
	      }

	      //handle weird case
	      //the status 2 tau has a status 23 tau as mother. the e or mu have the status 23 tau as its mother
	      if ( gParticleMotherIndex[j] >= 0 && gParticleId[gParticleMotherIndex[j]] == gParticleId[j] 
		   && gParticleStatus[gParticleMotherIndex[j]] == 23
		   && (abs(gParticleId[k]) == 11 || abs(gParticleId[k]) == 13)
		   && gParticleMotherIndex[k] == gParticleMotherIndex[j]) {
		isLeptonicTau = true;
		break;
	      }

	    }
	    if (!isLeptonicTau) {
	      genLeptonIndex.push_back(j);
	    }	   
	  }

	} //loop over gen particles
	


	//******************************************
	//sort gen leptons by pt
	//******************************************
	int tempIndex = -1;
	for(uint i = 0; i < genLeptonIndex.size() ; i++) {
	  for (uint j=0; j < genLeptonIndex.size()-1; j++) {
	    if (gParticlePt[genLeptonIndex[j+1]] > gParticlePt[genLeptonIndex[j]]) { 
	      tempIndex = genLeptonIndex[j];             // swap elements
	      genLeptonIndex[j] = genLeptonIndex[j+1];
	      genLeptonIndex[j+1] = tempIndex;
	    }
	  }
	}



	// for (int i=0;i<genLeptonIndex.size();i++) cout << "Lepton " << i << " : " << gParticleId[genLeptonIndex[i]] << " | " 
	// 					       << gParticlePt[genLeptonIndex[i]] << " "
	// 					       << gParticleEta[genLeptonIndex[i]] << " "
	// 					       << gParticlePhi[genLeptonIndex[i]] << " " 
	// 					       << " \n";
	// cout << "\n";


	// if (genLeptonIndex.size() >= 3) {
	//   cout << "\n";
	//   cout << "\n";
	//   for (int i=0;i<genLeptonIndex.size();i++) cout << "Lepton " << i << " : " << gParticleId[genLeptonIndex[i]] << " | " 
	// 					     << gParticlePt[genLeptonIndex[i]] << " "
	// 					     << gParticleEta[genLeptonIndex[i]] << " "
	// 					     << gParticlePhi[genLeptonIndex[i]] << " " 
	// 					     << " \n";
						  
	//   cout << "\n";

	//   for(int j = 0; j < nGenParticle; j++){
	//     cout << "Particle " << j << " : " << gParticleId[j] << " " << gParticleStatus[j] << " | "
	// 	 << gParticlePt[j] << " "
	// 	 << gParticleEta[j] << " "
	// 	 << gParticlePhi[j] << " "
	// 	 << " | " << gParticleMotherId[j] << " , " << gParticleMotherIndex[j] 
	// 	 << "\n";
	//   }
	// } 

	events->genlep1.SetPtEtaPhiM(0,0,0,0);
	events->genlep2.SetPtEtaPhiM(0,0,0,0);
	events->genlep1Type=0;
	events->genlep2Type=0;
	for (uint i=0;i<genLeptonIndex.size();i++) {
	  if(i==0) {
	    double mass = 0.000511;
	    if (abs(gParticleId[genLeptonIndex[i]]) == 13) mass = 0.1057;
	    if (abs(gParticleId[genLeptonIndex[i]]) == 15) mass = 1.777;	    
	    events->genlep1.SetPtEtaPhiM(gParticlePt[genLeptonIndex[i]],gParticleEta[genLeptonIndex[i]],gParticlePhi[genLeptonIndex[i]],mass);
	    events->genlep1Type = gParticleId[genLeptonIndex[i]];
	  }
	  if(i==1) {
	    double mass = 0.000511;
	    if (abs(gParticleId[genLeptonIndex[i]]) == 13) mass = 0.1057;
	    if (abs(gParticleId[genLeptonIndex[i]]) == 15) mass = 1.777;	    
	    events->genlep2.SetPtEtaPhiM(gParticlePt[genLeptonIndex[i]],gParticleEta[genLeptonIndex[i]],gParticlePhi[genLeptonIndex[i]],mass);
	    events->genlep2Type = gParticleId[genLeptonIndex[i]];
	  }				      
	}


	//***************************************************************************
	//Find the Good Leptons
	//***************************************************************************
	vector<TLorentzVector> GoodElectronsAndMuons;//leptons used to compute hemispheres
	vector<TLorentzVector> GoodLeptons;//leptons used to compute hemispheres

        for(int i = 0; i < nMuons; i++){
            if(muonPt[i] < 5) continue;
            if(abs(muonEta[i]) > 2.4) continue;                                  	   
	    if(!isVetoMuon(i)) continue;  
	    TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
            GoodLeptons.push_back(thisMuon);
	    GoodElectronsAndMuons.push_back(thisMuon);
        }

        for(int i = 0; i < nElectrons; i++){
            if(elePt[i] < 5) continue;
            if(fabs(eleEta[i]) > 2.5) continue;

	    //don't count electrons that were already selected as muons
	    bool alreadySelected = false;
	    for (uint j=0; j<GoodLeptons.size(); j++) {
	      if ( deltaR(GoodLeptons[j].Eta(), GoodLeptons[j].Phi(), eleEta[i],elePhi[i]) < 0.1) alreadySelected = true;
	    }
	    if (alreadySelected) continue;
	    
	    if(!isMVANonTrigVetoElectron(i)) continue; 
	    TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
            GoodLeptons.push_back(thisElectron);        
	    GoodElectronsAndMuons.push_back(thisElectron);
       }     
	for(int i = 0; i < nTaus; i++){	 
	  if (tauPt[i] < 20) continue;
	  if (fabs(tauEta[i]) > 2.4) continue;	  
	  if (!isLooseTau(i)) continue;
	  TLorentzVector thisTau; thisTau.SetPtEtaPhiM(tauPt[i], tauEta[i], tauPhi[i], 1.777);
	  GoodLeptons.push_back(thisTau);  
	}
	


	//********************************
	//Find Lepton 1
	//********************************
	bool foundLep1 = false;
	events->lep1.SetPtEtaPhiM(0,0,0,0);
	events->lep1Type = 0;
	events->lep1MatchedGenLepIndex = -1;
	events->lep1PassVeto = false;
	events->lep1PassLoose = false;
	events->lep1PassTight = false;
	events->lep1PassVetoID = false;
	events->lep1PassLooseID = false;
	events->lep1PassTightID = false;
	events->lep1PassVetoIso = false;
	events->lep1PassLooseIso = false;
	events->lep1PassTightIso = false;
	events->lep1MinDRToBJet = 0;

        for(int i = 0; i < nMuons; i++){

	  if(abs(muonEta[i]) > 2.4) continue;
	  if(!(muonPt[i] >= 10)) continue;
	  if (!isVetoMuon(i)) continue;

	  //for dilepton selection, require first lepton to be tight
	  if (option == 0) {
	    if(!isTightMuon(i)) continue;
	  }

	  if (!foundLep1
	      || muonPt[i] > events->lep1.Pt()
	      ) {
	    foundLep1 = true;
	    
	    events->lep1.SetPtEtaPhiM(muonPt[i],muonEta[i],muonPhi[i], 0.1057);
	    events->lep1Type = 13 * -1 * muonCharge[i];
	    events->lep1PassTight = isTightMuon(i);
	    events->lep1PassLoose = isLooseMuon(i);
	    events->lep1PassVeto = isVetoMuon(i);
	    events->lep1PassTightID = isTightMuon(i, true, false);
	    events->lep1PassLooseID = isLooseMuon(i, true, false);
	    events->lep1PassVetoID = isVetoMuon(i, true, false);
	    events->lep1PassTightIso = isTightMuon(i, false, true);
	    events->lep1PassLooseIso = isLooseMuon(i, false, true);
	    events->lep1PassVetoIso = isVetoMuon(i, false, true);
	  }		  
        }

        for(int i = 0; i < nElectrons; i++){

	  if(fabs(eleEta[i]) > 2.5) continue;
	  if(!(elePt[i] > 10)) continue;
	  if(!isMVANonTrigVetoElectron(i)) continue;

	  if(!(isTightElectron(i) && elePt[i] > 10)) continue;
	        
	   if (!foundLep1
	       || elePt[i] > events->lep1.Pt()
	       ) {
	     foundLep1 = true;
	     
	     events->lep1.SetPtEtaPhiM(elePt[i],eleEta[i],elePhi[i], 0.000511);
	     events->lep1Type = 11 * -1 * eleCharge[i];
	     events->lep1PassTight = isTightElectron(i);
	     events->lep1PassLoose = isLooseElectron(i);
	     events->lep1PassVeto = isMVANonTrigVetoElectron(i);
	     events->lep1PassTightID = isTightElectron(i, true, false);
	     events->lep1PassLooseID = isLooseElectron(i, true, false);
	     events->lep1PassVetoID = passMVANonTrigVetoElectronID(i);
	     events->lep1PassTightIso = isTightElectron(i, false, true);
	     events->lep1PassLooseIso = isLooseElectron(i, false, true);
	     events->lep1PassVetoIso = passMVANonTrigVetoElectronIso(i);
	     
	   }
        }   
	
        for(int i = 0; i < nTaus; i++){

	  //skip taus for lepton 1 if using dilepton selection
	  if (option == 0) continue;

	  if (tauPt[i] < 20) continue;
	  if (abs(tauEta[i]) > 2.4) continue;
	  if (!isLooseTau(i)) continue;

	  //veto electrons and muons
	  bool overlap = false;
	  for(auto& lep : GoodElectronsAndMuons){
	    if( deltaR(tauEta[i],tauPhi[i],lep.Eta(),lep.Phi()) < 0.1) overlap = true;
	  }
	  if (overlap) continue;

	  if (!foundLep1
	      || tauPt[i] > events->lep1.Pt()
	      ) {
	    foundLep1 = true;	   

	    events->lep1.SetPtEtaPhiM(tauPt[i],tauEta[i],tauPhi[i], 1.777);
	    events->lep1Type = 15;
	    events->lep1PassTight = isTightTau(i);
	    events->lep1PassLoose = isLooseTau(i);
	    events->lep1PassVeto = isLooseTau(i);
	    events->lep1PassTightID = false;
	    events->lep1PassLooseID = false;
	    events->lep1PassVetoID =  false;
	    events->lep1PassTightIso = false;
	    events->lep1PassLooseIso = false;
	    events->lep1PassVetoIso =  false;	    	  
	  }
        }



	if (events->lep1.Pt() > 0) {
	  events->lep1MatchedGenLepIndex = -1;
	  if (deltaR(events->lep1.Eta(),events->lep1.Phi(),events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
	      //&& abs(events->lep1Type) == abs(events->genlep1Type)
	      ) {
	    events->lep1MatchedGenLepIndex = 1;
	  }
	  if (deltaR(events->lep1.Eta(),events->lep1.Phi(),events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
	      //&& abs(events->lep1Type) == abs(events->genlep2Type)
	      ) {
	    events->lep1MatchedGenLepIndex = 2;
	  }
	}
	
	//********************************
	//Find Lepton 2
	//********************************
	bool foundLep2 = false;
	double lep2Iso = 0;
	events->lep2.SetPtEtaPhiM(0,0,0,0);
	events->lep2Type = 0;
	events->lep2MatchedGenLepIndex = -1;
	events->lep2PassVeto = false;
	events->lep2PassLoose = false;
	events->lep2PassTight = false;
	events->lep2PassVetoID = false;
	events->lep2PassLooseID = false;
	events->lep2PassTightID = false;
	events->lep2PassVetoIso = false;
	events->lep2PassLooseIso = false;
	events->lep2PassTightIso = false;
	events->lep2MinDRToBJet = 0;


        for(int i = 0; i < nMuons; i++){

	  if(abs(muonEta[i]) > 2.4) continue;
	  if(!(isVetoMuon(i, true, false) && muonPt[i] >= 5)) continue;
	  
	  //don't overlap with lepton 1
	  if (events->lep1.Pt() > 0) {
	    if (deltaR(events->lep1.Eta(),events->lep1.Phi(),muonEta[i],muonPhi[i]) < 0.1) continue;
	  }
	  
	  if (!foundLep2	      
	      || muonPt[i] > events->lep2.Pt()
	      ) {
	    foundLep2 = true;
	    
	    events->lep2.SetPtEtaPhiM(muonPt[i],muonEta[i],muonPhi[i], 0.1057);
	    events->lep2Type = 13 * -1 * muonCharge[i];
	    events->lep2PassTight = isTightMuon(i);
	    events->lep2PassLoose = isLooseMuon(i);
	    events->lep2PassVeto = isVetoMuon(i);
	    events->lep2PassTightID = isTightMuon(i, true, false);
	    events->lep2PassLooseID = isLooseMuon(i, true, false);
	    events->lep2PassVetoID = isVetoMuon(i, true, false);
	    events->lep2PassTightIso = isTightMuon(i, false, true);
	    events->lep2PassLooseIso = isLooseMuon(i, false, true);
	    events->lep2PassVetoIso = isVetoMuon(i, false, true);
	    lep2Iso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
	  }		  
        }

        for(int i = 0; i < nElectrons; i++){

	  if(fabs(eleEta[i]) > 2.5) continue;
	  if(!(passMVANonTrigVetoElectronID(i) && elePt[i] > 5)) continue;

	  //don't overlap with lepton 1
	  if (events->lep1.Pt() > 0) {
	    if (deltaR(events->lep1.Eta(),events->lep1.Phi(),eleEta[i],elePhi[i]) < 0.1) continue;
	  }

	  if (!foundLep2
	      || elePt[i] > events->lep2.Pt()
	      ) {
	    foundLep2 = true;
	    
	    events->lep2.SetPtEtaPhiM(elePt[i],eleEta[i],elePhi[i], 0.000511);
	    events->lep2Type = 11 * -1 * eleCharge[i];
	    events->lep2PassTight = isTightElectron(i);
	    events->lep2PassLoose = isLooseElectron(i);
	    events->lep2PassVeto = isMVANonTrigVetoElectron(i);
	    events->lep2PassTightID = isTightElectron(i, true, false);
	    events->lep2PassLooseID = isLooseElectron(i, true, false);
	    events->lep2PassVetoID = passMVANonTrigVetoElectronID(i);
	    events->lep2PassTightIso = isTightElectron(i, false, true);
	    events->lep2PassLooseIso = isLooseElectron(i, false, true);
	    events->lep2PassVetoIso = passMVANonTrigVetoElectronIso(i);
	    
	  }
        }   

        for(int i = 0; i < nTaus; i++){

	  if (tauPt[i] < 20) continue;
	  if (abs(tauEta[i]) > 2.4) continue;
	  if (!isLooseTau(i)) continue;

	  //don't overlap with lepton 1
	  if (events->lep1.Pt() > 0) {
	    if (deltaR(events->lep1.Eta(),events->lep1.Phi(),tauEta[i],tauPhi[i]) < 0.1) continue;
	  }

	  //veto electrons and muons
	  bool overlap = false;
	  for(auto& lep : GoodLeptons){
	    if( deltaR(tauEta[i],tauPhi[i],lep.Eta(),lep.Phi()) < 0.1) overlap = true;
	  }
	  if (overlap) continue;


	  if (!foundLep2
	      || tauPt[i] > events->lep2.Pt()
	      ) {
	    foundLep2 = true;
	    
	    events->lep2.SetPtEtaPhiM(tauPt[i],tauEta[i],tauPhi[i], 1.777);
	    events->lep2Type = 15;
	    events->lep2PassTight = isTightTau(i);
	    events->lep2PassLoose = isLooseTau(i);
	    events->lep2PassVeto = isLooseTau(i);
	    events->lep2PassTightID = false;
	    events->lep2PassLooseID = false;
	    events->lep2PassVetoID =  false;
	    events->lep2PassTightIso = false;
	    events->lep2PassLooseIso = false;
	    events->lep2PassVetoIso =  false;	    	  
	  }
        }

	if (events->lep2.Pt() > 0) {
	  events->lep2MatchedGenLepIndex = -1;
	  if (deltaR(events->lep2.Eta(),events->lep2.Phi(),events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
	      //&& abs(events->lep2Type) == abs(events->genlep1Type)
	      ) {
	    events->lep2MatchedGenLepIndex = 1;
	  }
	  if (deltaR(events->lep2.Eta(),events->lep2.Phi(),events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
	      //&& abs(events->lep2Type) == abs(events->genlep2Type)
	      ) {
	    events->lep2MatchedGenLepIndex = 2;
	  }
	}
	

	//***************************************************************************
	//Find the Jets
	//***************************************************************************
	bool bjet1Found = false;
	bool bjet2Found = false;
        vector<TLorentzVector> GoodJets;
        int numJetsAbove80GeV = 0;
	int numJetsAbove40GeV = 0;
	int nBJetsLoose20GeV = 0;
	int nBJetsMedium20GeV = 0;
	int nBJetsTight20GeV = 0;
	events->bjet1.SetPtEtaPhiM(0,0,0,0);
	events->bjet2.SetPtEtaPhiM(0,0,0,0);
	events->bjet1PassLoose = false;
	events->bjet1PassMedium = false;
	events->bjet1PassTight = false;
	events->bjet2PassLoose = false;
	events->bjet2PassMedium = false;
	events->bjet2PassTight = false;       
	events->minDPhi = 9999;
	events->minDPhiN = 9999;

        for(int i = 0; i < nJets; i++){

	  //exclude selected muons and electrons from the jet collection
	  double dR = -1;
	  for(auto& lep : GoodLeptons){
	    double thisDR = deltaR(jetEta[i],jetPhi[i],lep.Eta(),lep.Phi());
	    if(dR < 0 || thisDR < dR) dR = thisDR;
	  }
	  if(dR > 0 && dR < 0.4) continue; //jet matches a selected lepton
	  

	  double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
						 fixedGridRhoFastjetAll, jetJetArea[i], 
						 JetCorrector);   
	  TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC);
	  
	  //find the best 2 btagged jets
	  if (isCSVL(i)) {
	    if (!bjet1Found) {
	      bjet1Found = true;
	      events->bjet1.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), thisJet.M());
	      events->bjet1PassLoose = false;
	      events->bjet1PassMedium = false;
	      events->bjet1PassTight = false;
	      if(isCSVL(i)) events->bjet1PassLoose = true;	      
	      if(isCSVM(i)) events->bjet1PassMedium = true;
	      if(isCSVT(i)) events->bjet1PassTight = true;	      
	    } else if ( (isCSVT(i) && !events->bjet1PassTight)
			|| (isCSVM(i) && !events->bjet1PassMedium)			
			|| ( ( (isCSVT(i) && events->bjet1PassTight)
			       ||  (isCSVM(i) && events->bjet1PassMedium)
			       ||  (isCSVL(i) && events->bjet1PassLoose) )
			     && thisJet.Pt() > events->bjet1.Pt() )
			) {
	      events->bjet2.SetPtEtaPhiM(events->bjet1.Pt(), events->bjet1.Eta(), events->bjet1.Phi(), events->bjet1.M());
	      events->bjet2PassLoose = events->bjet1PassLoose;
	      events->bjet2PassMedium = events->bjet1PassMedium;
	      events->bjet2PassTight = events->bjet1PassTight;
	    
	      events->bjet1.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), thisJet.M());
	      events->bjet1PassLoose = false;
	      events->bjet1PassMedium = false;
	      events->bjet1PassTight = false;
	      if(isCSVL(i)) events->bjet1PassLoose = true;	      
	      if(isCSVM(i)) events->bjet1PassMedium = true;
	      if(isCSVT(i)) events->bjet1PassTight = true;
	    } else {
	      if (!bjet2Found 
		  || (isCSVT(i) && !events->bjet2PassTight)
		  || (isCSVM(i) && !events->bjet2PassMedium)			
		  || ( ( (isCSVT(i) && events->bjet2PassTight)
			 ||  (isCSVM(i) && events->bjet2PassMedium)
			 ||  (isCSVL(i) && events->bjet2PassLoose) )
		       && thisJet.Pt() > events->bjet2.Pt() )
		  ) {
		bjet2Found = true;
		events->bjet2.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), thisJet.M());
		events->bjet2PassLoose = false;
		events->bjet2PassMedium = false;
		events->bjet2PassTight = false;
		if(isCSVL(i)) events->bjet2PassLoose = true;	      
		if(isCSVM(i)) events->bjet2PassMedium = true;
		if(isCSVT(i)) events->bjet2PassTight = true;	      
	      }
	    }
	  }

	  if (jetPt[i]*JEC > 20 && isCSVL(i)) nBJetsLoose20GeV++;
	  if (jetPt[i]*JEC > 20 && isCSVM(i)) nBJetsMedium20GeV++;
	  if (jetPt[i]*JEC > 20 && isCSVT(i)) nBJetsTight20GeV++;

	  if(jetPt[i]*JEC < 30) continue;
	  if(fabs(jetEta[i]) > 3.0) continue;
	  
	  if ( (events->bjet1.Pt() <= 0 || deltaR(events->bjet1.Eta(),events->bjet1.Phi(),thisJet.Eta(),thisJet.Phi()) > 0.01)
	       &&
	       (events->bjet2.Pt() <= 0 || deltaR(events->bjet2.Eta(),events->bjet2.Phi(),thisJet.Eta(),thisJet.Phi()) > 0.01) 
	       ) {
	    numJetsAbove40GeV++;
	  }

	  if(jetPt[i]*JEC > 80) numJetsAbove80GeV++;
	  
	  GoodJets.push_back(thisJet);	  
	  
        } //loop over jets

	//*****************************************************
	//Fill lepton minDR to closest bjets
	//*****************************************************
	events->lep1MinDRToBJet = 9999;
	if (events->lep1.Pt() > 0) {
	  if (events->bjet1.Pt() > 0 && events->bjet1PassLoose) 
	    events->lep1MinDRToBJet = deltaR(events->bjet1.Eta(),events->bjet1.Phi(),events->lep1.Eta(),events->lep1.Phi());
	  if (events->bjet2.Pt() > 0 && events->bjet2PassLoose && 
	      deltaR(events->bjet2.Eta(),events->bjet2.Phi(),events->lep1.Eta(),events->lep1.Phi()) < events->lep1MinDRToBJet)
	    events->lep1MinDRToBJet = deltaR(events->bjet2.Eta(),events->bjet2.Phi(),events->lep1.Eta(),events->lep1.Phi());
	}
	events->lep2MinDRToBJet = 9999;
	if (events->lep2.Pt() > 0) {
	  if (events->bjet1.Pt() > 0 && events->bjet1PassLoose) 
	    events->lep2MinDRToBJet = deltaR(events->bjet1.Eta(),events->bjet1.Phi(),events->lep2.Eta(),events->lep2.Phi());
	  if (events->bjet2.Pt() > 0 && events->bjet2PassLoose && 
	      deltaR(events->bjet2.Eta(),events->bjet2.Phi(),events->lep2.Eta(),events->lep2.Phi()) < events->lep2MinDRToBJet)
	    events->lep2MinDRToBJet = deltaR(events->bjet2.Eta(),events->bjet2.Phi(),events->lep2.Eta(),events->lep2.Phi());
	}
	

	//*****************************************************
	//sort good jets
	//*****************************************************
	TLorentzVector tmpjet;
	for (int i=0;i<int(GoodJets.size());i++) {
	  for (int j=0;j<int(GoodJets.size()-1);j++) {
	    if (GoodJets[j+1].Pt() > GoodJets[j].Pt()) {
	      tmpjet = GoodJets[j]; //swap elements
	      GoodJets[j] = GoodJets[j+1];
	      GoodJets[j+1] = tmpjet;
	    }
	  }
	}

	//Fill two leading jets
	events->jet1.SetPtEtaPhiM(0,0,0,0);
	events->jet2.SetPtEtaPhiM(0,0,0,0);
	events->jet1PassCSVLoose = false;
	events->jet1PassCSVMedium = false;
	events->jet1PassCSVTight = false;
	events->jet2PassCSVLoose = false;
	events->jet2PassCSVMedium = false;
	events->jet2PassCSVTight = false; 
      

	//compute minDPhi variables
	double JER = 0.1; //average jet energy resolution
	for (int i=0;i<int(GoodJets.size());i++) {
	  if (i>2) continue; //only consider first 3 jets for minDPhi 	  
	  double dPhi = fmin( fabs(deltaPhi(metPhi,GoodJets[i].Phi())) , fabs( 3.1415926 - fabs(deltaPhi(metPhi,GoodJets[i].Phi()))));
	  if (dPhi < events->minDPhi) events->minDPhi = dPhi;
	  double deltaT = 0;
	  for(auto& jet2 : GoodJets) {
	    deltaT += pow(JER*jet2.Pt()*sin(fabs(deltaPhi(GoodJets[i].Phi(),jet2.Phi()))),2);
	  }
	  double dPhiN = dPhi / atan2( sqrt(deltaT) , metPt);
	  if (dPhiN < events->minDPhiN) events->minDPhiN = dPhiN;
	}

	//Make Good Jet Collection, excluding the leading jet
	vector<TLorentzVector> GoodJets_NoLeadJet;
	int leadJetIndex = -1;
	double leadJetPt = 0;
	for (int i=0;i<int(GoodJets.size());i++) {
	  if (GoodJets[i].Pt() > leadJetPt) {
	    leadJetIndex = i;
	    leadJetPt = GoodJets[i].Pt();
	  }
	}
	int numJetsAbove80GeV_NoLeadJet = 0;
	for (int i=0;i<int(GoodJets.size());i++) {
  	  if (i != leadJetIndex) {
	    GoodJets_NoLeadJet.push_back(GoodJets[i]);
	    if (GoodJets[i].Pt() > 80) numJetsAbove80GeV_NoLeadJet++;
	  }
	}

        //Compute the razor variables using the selected jets and possibly leptons
        vector<TLorentzVector> GoodPFObjects;
        vector<TLorentzVector> GoodPFObjects_NoLeadJet;
        for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
        for(auto& jet : GoodJets_NoLeadJet) GoodPFObjects_NoLeadJet.push_back(jet);
        if(passedLeptonicTrigger) for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);
        if(passedLeptonicTrigger) for(auto& lep : GoodLeptons) GoodPFObjects_NoLeadJet.push_back(lep);
        TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
        TLorentzVector PFMET_NoLeadJet = PFMET; if (leadJetIndex >= 0) PFMET_NoLeadJet = PFMET + GoodJets[leadJetIndex];

	events->MR = 0;
	events->Rsq = 0;

	//only compute razor variables if we have 2 jets above 80 GeV
	if (numJetsAbove80GeV >= 2) {
	  vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
	  events->MR = computeMR(hemispheres[0], hemispheres[1]); 
	  events->Rsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
	  events->dPhiRazor = deltaPhi(hemispheres[0].Phi(),hemispheres[1].Phi());
	}

	if (numJetsAbove80GeV_NoLeadJet >= 2) {
	  vector<TLorentzVector> hemispheres_NoLeadJet = getHemispheres(GoodPFObjects_NoLeadJet);
	  events->MR_NoLeadJet = computeMR(hemispheres_NoLeadJet[0], hemispheres_NoLeadJet[1]); 
	  events->Rsq_NoLeadJet = computeRsq(hemispheres_NoLeadJet[0], hemispheres_NoLeadJet[1], PFMET_NoLeadJet);
	}
	events->MET = metPt;
	events->MET_NoLeadJet = PFMET_NoLeadJet.Pt();
	events->NJets40 = numJetsAbove40GeV;
	events->NBJetsLoose = nBJetsLoose20GeV;
	events->NBJetsMedium = nBJetsMedium20GeV;
	events->NBJetsTight = nBJetsTight20GeV;

	events->HT = 0;
	for(auto& pfobj : GoodPFObjects) events->HT += pfobj.Pt();
	
	//compute M_T for lep1 and MET
	events->lep1MT = sqrt(events->lep1.M2() + 2*metPt*events->lep1.Pt()*(1 - cos(deltaPhi(metPhi,events->lep1.Phi()))));
	

	//fill event 
	events->tree_->Fill();


	//*****************************************************************************************************************
	//DEBUGGING
	//*****************************************************************************************************************
	bool printdebug = false;
	// if (events->bjet1.Pt() > 30 && events->bjet2.Pt() > 30 && events->bjet1PassMedium && events->bjet2PassMedium &&
	//     events->lep1.Pt() > 25 && events->lep1PassTight && events->lep2.Pt() > 25 
	//     //&& events->lep2Type*events->lep1Type < 0
	//     //&& abs(events->lep1Type) != abs(events->lep2Type) 
	//     && events->lep2PassVeto
	//     //&& abs(events->lep2Type) == 13
	//     && events->lep2MatchedGenLepIndex < 0
	//     //&& !events->lep2PassVeto
	//     ) prindebug = true;

	if (printdebug) { 

	  cout << "Tight dilepton event: " << eventNum << "\n";	    
	  cout << "lep1: " << events->lep1Type << " | " << events->lep1.Pt() << " " << events->lep1.Eta() << " " << events->lep1.Phi() 
	       << " | " << events->lep1PassTight << " " << events->lep1PassLoose << " " << events->lep1PassVeto 
	       << " | " << events->lep1PassTightID << " " << events->lep1PassLooseID << " " << events->lep1PassVetoID 
	       << " | " << events->lep1PassTightIso << " " << events->lep1PassLooseIso << " " << events->lep1PassVetoIso
	       << " | " << events->lep1MatchedGenLepIndex << "\n";
	  cout << "lep2: " << events->lep2Type << " | " << events->lep2.Pt() << " " << events->lep2.Eta() << " " << events->lep2.Phi() 
	       << " | " << events->lep2PassTight << " " << events->lep2PassLoose << " " << events->lep2PassVeto 
	       << " | " << events->lep2PassTightID << " " << events->lep2PassLooseID << " " << events->lep2PassVetoID 
	       << " | " << events->lep2PassTightIso << " " << events->lep2PassLooseIso  << " " << events->lep2PassVetoIso  
	       << " | " << events->lep2MatchedGenLepIndex << " | " << lep2Iso << "\n";	  
	  cout << "btag1: " << events->bjet1.Pt() << " " << events->bjet1.Eta()  << " " << events->bjet1.Phi() << " | " 
	       << events->bjet1PassLoose << " "  << events->bjet1PassMedium << " "  << events->bjet1PassTight << " \n";
	  cout << "btag2: " << events->bjet2.Pt() << " " << events->bjet2.Eta()  << " " << events->bjet2.Phi() << " | " 
	       << events->bjet2PassLoose << " "  << events->bjet2PassMedium << " "  << events->bjet2PassTight << " \n";
	  
	  cout << "genlep1: " << events->genlep1Type << " | " << events->genlep1.Pt() << " " << events->genlep1.Eta() << " " << events->genlep1.Phi() << "\n";
	  cout << "genlep2: " << events->genlep2Type << " | " << events->genlep2.Pt() << " " << events->genlep2.Eta() << " " << events->genlep2.Phi() << "\n";
	      
	  cout << "lep2 DR: " << deltaR(events->lep2.Eta(),events->lep2.Phi(),events->genlep1.Eta(),events->genlep1.Phi()) << " "
	       << deltaR(events->lep2.Eta(),events->lep2.Phi(),events->genlep2.Eta(),events->genlep2.Phi()) << " : " 
	       << events->lep2Type << " " << events->genlep1Type << " " << events->genlep2Type << "\n";
	  cout << "lep2 DR to BJet: " 
	       << deltaR(events->bjet1.Eta(),events->bjet1.Phi(),events->lep2.Eta(),events->lep2.Phi()) << " " 
	       << deltaR(events->bjet2.Eta(),events->bjet2.Phi(),events->lep2.Eta(),events->lep2.Phi()) << " : " 
	       << events->lep2MinDRToBJet  << "\n";
	

	  for(int i = 0; i < nMuons; i++){
            if(muonPt[i] < 5) continue;
            if(abs(muonEta[i]) > 2.4) continue;                                  	   
	    cout << "Muon : " << i << " : " << muonPt[i] << " " << muonEta[i] << " " << muonPhi[i] << " : " << isVetoMuon(i, true, false) << " " << isVetoMuon(i, false, true) << " " << isVetoMuon(i) << "\n";
	  }
	  for(int i = 0; i < nElectrons; i++){
            if(elePt[i] < 5) continue;
            if(fabs(eleEta[i]) > 2.5) continue;
	    cout << "Muon : " << i << " : " << muonPt[i] << " " << muonEta[i] << " " << muonPhi[i] << " : " << passMVANonTrigVetoElectronID(i) << " " << passMVANonTrigVetoElectronIso(i) << " " << isMVANonTrigVetoElectron(i) << "\n";	   	
	  }     
	  


	  for(int j = 0; j < nGenParticle; j++){
	    cout << "Particle " << j << " : " << gParticleId[j] << " " << gParticleStatus[j] << " | "
		 << gParticlePt[j] << " "
		 << gParticleEta[j] << " "
		 << gParticlePhi[j] << " "
		 << " | " << gParticleMotherId[j] << " , " << gParticleMotherIndex[j] 
		 << "\n";	 
	  }
	  cout << "\n\n\n";	    	     
	}



    }//end of event loop


    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}






//ControlSampleEvent->Draw("probeMatchedGenLepIndex","bjet1.Pt()>30 && bjet2.Pt() > 30 && bjet1PassLoose && bjet2PassLoose && tag.Pt() > 25 && tagPassTight && probe.Pt()>25 && probeType*tagType < 0 && abs(tagType) != abs(probeType) && !probePassVeto")
//ControlSampleEvent->Draw("lep2MinDRToBJet","bjet1.Pt()>30 && bjet2.Pt() > 30 && bjet1PassLoose && bjet2PassLoose && lep1.Pt() > 25 && lep1PassTight && lep2.Pt()>25 && lep2Type*lep1Type < 0 && abs(lep1Type) != abs(lep2Type) && lep2PassVetoID && lep2MatchedGenLepIndex > 0 ")
//TagAndProbeEvent->Draw("probeMinDRToBJet","bjet1.Pt()>30 && bjet2.Pt() > 30 && bjet1PassLoose && bjet2PassLoose && tag.Pt() > 25 && tagPassTight && probe.Pt()>25 && probeType*tagType < 0 && abs(tagType) != abs(probeType) && probePassVetoID && probeMatchedGenLepIndex > 0 ")

//ControlSampleEvent->Draw("lep2MatchedGenLepIndex > 0","bjet1.Pt()>30 && bjet2.Pt() > 30 && bjet1PassLoose && bjet2PassLoose && lep1.Pt() > 25 && lep1PassTight && lep2.Pt()>25 && lep2Type*lep1Type < 0 && abs(lep1Type) != abs(lep2Type) && lep2PassVetoID  ")
//TagAndProbeEvent->Draw("probeMatchedGenLepIndex > 0","bjet1.Pt()>30 && bjet2.Pt() > 30 && bjet1PassLoose && bjet2PassLoose && tag.Pt() > 25 && tagPassTight && probe.Pt()>25 && probeType*tagType < 0 && abs(tagType) != abs(probeType) && probePassVetoID ")
//ControlSampleEvent->Draw("lep1MinDRToBJet","bjet1.Pt()>30 && bjet2.Pt() > 30 && bjet1PassLoose && bjet2PassLoose && lep1.Pt() > 25 && lep1PassVetoID && lep1MatchedGenLepIndex > 0")
