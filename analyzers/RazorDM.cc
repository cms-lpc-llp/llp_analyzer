#include "RazorDM.h"
#include "JetCorrectorParameters.h"
#include "RazorHelper.h"
#include "TMath.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

void RazorDM::Analyze(bool isData, int option, string outFileName, string label)
{
  //initialization: create one TTree for each analysis box 
  cout << "Initializing..." << endl;
  bool isFastsimSMS = (option == 1 || option == 11);
  bool combineTrees = true;
  if ( outFileName.empty() )
  {
      cout << "RazorDM: Output filename not specified!" << endl << "Using default output name RazorDM.root" << endl;
      outFileName = "RazorDM.root";
  }
  
  if (combineTrees) cout << "Using combineTrees" << endl;
  else cout << "Using razorBoxes" << endl;
 
  int fillCounter = 0; 
  TFile outFile(outFileName.c_str(), "RECREATE");

  //one tree to hold all events
  TTree *razorTree = new TTree("RazorDM", "Info on selected razor DM events");

  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  if (cmsswPath == NULL) 
  {
      cout << "Warning: CMSSW_BASE not detected. Exiting..." << endl;
      return;
  }
  
  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
    
  //Initialize helper
  RazorHelper *helper = 0;
  string analysisTag = "Razor2016_MoriondRereco";
  if ( label != "") analysisTag = label; 
  helper = new RazorHelper(analysisTag, isData, isFastsimSMS);

  // Get jet corrector
  std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
  std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();
 
  //tree variables
  int nSelectedJets, nBTaggedJetsL, nBTaggedJetsM, nBTaggedJetsT, numJetsAbove80GeV;
  int nSecondaryVertices;
  int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus, nLooseTaus, nLoosePhotons, nTightPhotons;
  UInt_t run, lumi, event;
  float MR, deltaPhi;
  float HT, MHT;
  float alphaT, dPhiMin;
  float Rsq;
  float t1metPt, t1metPhi, caloMet;
  float JetPt_uncorr[30], JetEta_uncorr[30], JetPhi_uncorr[30], JetE_uncorr[30];
  float LeadJetNeutralHadronFraction, LeadJetChargedHadronFraction;
  float JetPt[30], JetEta[30], JetPhi[30], JetE[30];
  float leadingJetPt, leadingJetEta, subLeadingJetPt;
  
  
    //set branches on big tree
  if(combineTrees)
  {
      razorTree->Branch("run",&run,"run/i");
      razorTree->Branch("lumi",&lumi,"lumi/i");
      razorTree->Branch("event",&event,"event/i");
      razorTree->Branch("nSecondaryVertices", &nSecondaryVertices, "nSecondaryVertices/i");
   
      razorTree->Branch("numJetsAbove80GeV", &numJetsAbove80GeV, "numJetsAbove80GeV/I");
      razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
      razorTree->Branch("nBTaggedJetsL", &nBTaggedJetsL, "nBTaggedJetsL/I");
      razorTree->Branch("nBTaggedJetsM", &nBTaggedJetsM, "nBTaggedJetsM/I");
      razorTree->Branch("nBTaggedJetsT", &nBTaggedJetsT, "nBTaggedJetsT/I");
      razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
      razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
      razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
      razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
      razorTree->Branch("nLoosePhotons", &nLoosePhotons, "nLoosePhotons/I");
      razorTree->Branch("nTightPhotons", &nTightPhotons, "nTightPhotons/I");
      razorTree->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
      razorTree->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
      
      razorTree->Branch("JetE_uncorr", JetE_uncorr, "JetE_uncorr[nSelectedJets]/F");
      razorTree->Branch("JetPt_uncorr", JetPt_uncorr, "JetPt_uncorr[nSelectedJets]/F");
      razorTree->Branch("JetPhi_uncorr", JetPhi_uncorr, "JetPhi_uncorr[nSelectedJets]/F");
      razorTree->Branch("JetEta_uncorr", JetEta_uncorr, "JetEta_uncorr[nSelectedJets]/F");

      razorTree->Branch("JetE", JetE, "JetE[nSelectedJets]/F");
      razorTree->Branch("JetPt", JetPt, "JetPt[nSelectedJets]/F");
      razorTree->Branch("JetEta", JetEta, "JetEta[nSelectedJets]/F");
      razorTree->Branch("JetPhi", JetPhi, "JetPhi[nSelectedJets]/F");
      razorTree->Branch("LeadJetNeutralHadronFraction", &LeadJetNeutralHadronFraction, "LeadJetNeutralHadronFraction/F"); 
      razorTree->Branch("LeadJetChargedHadronFraction", &LeadJetChargedHadronFraction, "LeadJetChargedHadronFraction/F"); 
      
      razorTree->Branch("alphaT", &alphaT, "alphaT/F");
      razorTree->Branch("dPhiMin", &dPhiMin, "dPhiMin/F");
      
      razorTree->Branch("MR", &MR, "MR/F");
      razorTree->Branch("Rsq", &Rsq, "Rsq/F");
      razorTree->Branch("deltaPhi", &deltaPhi, "deltaPhi/F");
      razorTree->Branch("metPt", &metPt, "metPt/F");
      razorTree->Branch("metPhi", &metPhi, "metPhi/F");
      razorTree->Branch("caloMet", &caloMet, "caloMet/F");
      razorTree->Branch("t1metPt", &t1metPt, "t1metPt/F");
      razorTree->Branch("t1metPhi", &t1metPhi, "t1metPhi/F");
      razorTree->Branch("HT", &HT, "HT/F");
      razorTree->Branch("MHT", &MHT, "MHT/F");
      
      razorTree->Branch("leadingJetPt",&leadingJetPt,"leadingJetPt/F");
      razorTree->Branch("leadingJetEta",&leadingJetEta,"leadingJetEta/F");
      razorTree->Branch("subLeadingJetPt",&subLeadingJetPt,"subLeadingJetPt/F");
      
      razorTree->Branch("HLTDecision", HLTDecision, "HLTDecision[160]/O");
    }
    
    //begin loop
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    cout << "NEntries = " << (Long64_t) nentries << endl;
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      //begin event
      if(jentry % 1000 == 0) cout << "Processing entry " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //fill normalization histogram
      NEvents->Fill(1.0);
      
      // event info
      run = runNum;
      lumi = lumiNum;
      event = eventNum;
      nSecondaryVertices = nSlimmedSecondV;
      
      //reset tree variables
      nSelectedJets = 0;
      numJetsAbove80GeV = 0;
      nBTaggedJetsL = 0;
	  nBTaggedJetsM = 0;
	  nBTaggedJetsT =0;
      nLooseMuons = 0;
      nTightMuons = 0;
      nLooseElectrons = 0;
      nTightElectrons = 0;
      nLoosePhotons = 0;
      nTightPhotons = 0;
      nLooseTaus = 0;
      nTightTaus = 0;
      alphaT    = -1.0;
      dPhiMin   = -1.0;
      MR        = -1.0;
      Rsq       = -1.0;
	  deltaPhi     = -1.0;
      leadingJetPt = -1;
      leadingJetEta = 0.0;
      subLeadingJetPt = -1;
	  LeadJetChargedHadronFraction = -1.; 
	  LeadJetNeutralHadronFraction = -1.; 
	  
	  for ( int j = 0; j < 30; j++ )
      {
          JetE_uncorr[j]   = 0.0;
	      JetPt_uncorr[j]  = 0.0;
	      JetPhi_uncorr[j] = 0.0;
	      JetEta_uncorr[j] = 0.0;
          JetE[j]  = 0.0;
	      JetPt[j]  = 0.0;
	      JetPhi[j] = 0.0;
	      JetEta[j] = 0.0;
      }
      //TODO: triggers!
      //------------------------------
	  // Reco+ID Muons
	  //------------------------------
      bool passedLeptonicTrigger = true;
      bool passedHadronicTrigger= true;
      if ( !(passedLeptonicTrigger || passedHadronicTrigger) ) continue; //ensure event passed a trigger
      vector<TLorentzVector> LooseMu; //Muons use to compute MET
      vector<TLorentzVector> GoodLeptons; //leptons used to compute hemispheres
      for (int i = 0; i < nMuons; i++)
      {
          if ( muonPt[i] <= 10.0 || fabs(muonEta[i]) >= 2.4 ) continue;  
	      TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
	      if ( isLooseMuon(i) )
          {
              nLooseMuons++;
	          LooseMu.push_back(thisMuon);
              GoodLeptons.push_back(thisMuon); 
	      }
          if (isTightMuon(i)) nTightMuons++;
      }
      
      //------------------------------
	  // Reco+ID Electron
	  //------------------------------
	  vector<TLorentzVector> LooseEle; //Electrons use to compute MET
      for(int i = 0; i < nElectrons; i++)
      {
        if (elePt[i] <= 10.0 || fabs(eleEta[i]) >= 2.5) continue;
        if (fabs(eleEta[i]) > 1.442 && fabs(eleEta[i]) < 1.566) continue;
        bool overlap = false;
        for(auto& lep : GoodLeptons)
        {
	        if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
        }
        if (overlap) continue;
	    TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
	    if (isLooseElectron(i))
        {
            nLooseElectrons++;
	        LooseEle.push_back(thisElectron);
            GoodLeptons.push_back(thisElectron);            
	    }
	    if (isTightElectron(i)) nTightElectrons++;
      }
	
	  //----------------------
	  //ID Tau and Veto
	  //----------------------
	  for(int i = 0; i < nTaus; i++)
      {
          if (tauPt[i] <= 18 || fabs(tauEta[i]) >= 2.3) continue;
          
          //remove overlaps
          bool overlap = false;
          for (auto& lep : GoodLeptons)
          {
              if (RazorAnalyzer::deltaR(tauEta[i],tauPhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
          }
          if (overlap) continue;
          if (!(tau_ID[i] & 1)) continue;
          if (tau_combinedIsoDeltaBetaCorr3Hits[i] > 5.) continue;
    
          TLorentzVector thisTau; thisTau.SetPtEtaPhiM(tauPt[i], tauEta[i], tauPhi[i], 1.777);
          GoodLeptons.push_back(thisTau);  
          nLooseTaus++;
          nTightTaus++;
	  }
	
      //------------------------
      // Photon ID and Veto
      // ----------------------
      for (int i = 0; i < nPhotons; i++)
      {
          if (phoPt[i] <= 20) continue;
          if (fabs(phoEta[i]) >= 2.5) continue;
          
          //remove overlaps
          bool overlap = false;
          for (auto& lep : GoodLeptons)
          {
              if (RazorAnalyzer::deltaR(phoEta[i],phoPhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
          }
          if (overlap) continue;
          if (fabs(pho_superClusterEta[i]) <= 1.479)
          {
              if (pho_HoverE[i] > 0.05) continue;
              if (phoFull5x5SigmaIetaIeta[i]    > 0.0102) continue;
              if (!photonPassesIsolation(i, 3.32, 1.92 + 0.014*phoPt[i] + 1.9e-5*TMath::Power(phoPt[i],2), 0.81 + 0.0053*phoPt[i], true)) continue;
          } 
          else 
          {
              if (pho_HoverE[i] > 0.05) continue;
              if (phoFull5x5SigmaIetaIeta[i]    > 0.0274) continue;
              if (!photonPassesIsolation(i, 1.97, 11.86 + 0.0139*phoPt[i] + 2.5e-5*TMath::Power(phoPt[i],2), 0.83 + 0.0034*phoPt[i], true)) continue;
          }
	      TLorentzVector thisPhoton = makeTLorentzVector(phoPt[i], phoEta[i], phoPhi[i], phoE[i]);
          nLoosePhotons++;
          GoodLeptons.push_back(thisPhoton);  
      } 


	  //------------------------
	  //Jet ID + Correction
	  //------------------------
	  vector<TLorentzVector> GoodJets;
	  vector<TLorentzVector> GoodJets_uncorr;
      vector<float> JetChargedHadronFraction;
      vector<float> JetNeutralHadronFraction;
	  for (int i = 0; i < nJets; i++)
	  {
      
        //*****************************************************************
        //exclude selected muons and electrons from the jet collection
        //*****************************************************************
	    
	    double deltaR = -1.0;
        for (auto& lep : GoodLeptons)
        {
	        double thisDR = RazorAnalyzer::deltaR(jetEta[i],jetPhi[i],lep.Eta(),lep.Phi());  
	        if (deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
        }
        
        if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton
	    if ( !jetPassIDLoose[i] ) continue;

	    //Obtain JEC
        double tmpRho = fixedGridRhoFastjetAll;
        double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
					       tmpRho, jetJetArea[i], runNum,
					       JetCorrectorIOV,JetCorrector);   
        double JECLevel1 = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
						     tmpRho, jetJetArea[i], runNum,
						     JetCorrectorIOV,JetCorrector, 0);   
        double jetEnergySmearFactor = 1.0;

//        if (jentry == 1) cout << "FixedGridRhoAll: " << fixedGridRhoAll << endl << "jetJetArea[" << i << "]: " << jetJetArea[i] << endl; 
	    
        //Apply JEC
        TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);
        TLorentzVector L1CorrJet = makeTLorentzVector(jetPt[i]*JECLevel1, jetEta[i], jetPhi[i], jetE[i]*JECLevel1);
        TLorentzVector UnCorrJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);      
        double jetCorrPt = jetPt[i]*JEC*jetEnergySmearFactor;

	    //b-tagging
        if (thisJet.Pt() > 15 && fabs(jetEta[i]) < 2.5)
        {
            if (jetCISV[i] > 0.46)  nBTaggedJetsL++;
            if (jetCISV[i] > 0.80)  nBTaggedJetsM++;
            if (jetCISV[i] > 0.935)  nBTaggedJetsT++;
        } 
        
        if(jetCorrPt <= 30) continue;
        if(fabs(jetEta[i]) >= 3.0) continue;
	    
	    if (thisJet.Pt() > 80.0) numJetsAbove80GeV++;
	    GoodJets_uncorr.push_back(UnCorrJet);
	    GoodJets.push_back(thisJet);
        JetNeutralHadronFraction.push_back(jetNeutralHadronEnergyFraction[i]);
        JetChargedHadronFraction.push_back(jetChargedHadronEnergyFraction[i]);
	    
        nSelectedJets++;
      }
	
	
	  int jIndex = 0;
	  for (auto& tmpJet: GoodJets_uncorr)
	  {
	    JetPt_uncorr[jIndex] = tmpJet.Pt();
	    JetE_uncorr[jIndex] = tmpJet.E();
	    JetPhi_uncorr[jIndex] = tmpJet.Phi();
	    JetEta_uncorr[jIndex] = tmpJet.Eta();
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
            LeadJetNeutralHadronFraction = JetNeutralHadronFraction.at(jIndex);
            LeadJetChargedHadronFraction = JetChargedHadronFraction.at(jIndex);
	    }
	    else if (JetPt[jIndex]>subLeadingJetPt && JetPt[jIndex]<14000) 
        {
	        subLeadingJetPt = JetPt[jIndex];
	    }
        auto sortJets = []( TLorentzVector a, TLorentzVector b )
        { 
            return (a.Pt() > b.Pt()); 
        };
        
        std::sort(GoodJets.begin() , GoodJets.end(), sortJets);
	    
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
      caloMet = metCaloPt;	
      vector<TLorentzVector> hemispheres = getHemispheres( GoodJets );
    
      MR    = computeMR(hemispheres[0], hemispheres[1]); 
	  Rsq   = computeRsq(hemispheres[0], hemispheres[1], t1PFMET);
      deltaPhi = fabs(hemispheres[0].DeltaPhi(hemispheres[1])); 

      // Compute HT and MHT
      float MhtX = 0., MhtY = 0.;
      HT = 0.; 
      for (auto& obj : GoodJets) 
      { 
        HT += obj.Pt(); MhtX += obj.Px(); MhtY += obj.Py(); 
      }

      TLorentzVector MyMHT;
      MyMHT.SetPxPyPzE(-MhtX, -MhtY, 0, sqrt(pow(MhtX,2) + pow(MhtY,2)));

      MHT = MyMHT.Pt();


      // EVENT SELECTION
//      if (MR < 200) continue;
//      if (Rsq < 0.35) continue;
//      if (deltaPhi > 2.5) continue;
      if (nLooseTaus > 0 || nLooseElectrons > 0 || nLooseMuons > 0 || nLoosePhotons > 0) continue;
      if (numJetsAbove80GeV < 2) continue; //event fails to have two 80 GeV jets

    // Monojet cuts:
      if (t1metPt < 200) continue;
      if (leadingJetPt < 100) continue;
      if (dPhiMin < 0.5) continue;
      if (fabs(leadingJetEta) > 2.5) continue;
      if (nBTaggedJetsM > 0) continue;
      if (LeadJetChargedHadronFraction < 0.1) continue;
      if (LeadJetNeutralHadronFraction > 0.8) continue;
      cout << "Run " << run << "\tEvt " << event << endl;       
      razorTree->Fill();
      fillCounter++; 
    } //end of event loop
    cout << "Writing output trees..." << endl;
    outFile.cd();
    if(combineTrees) razorTree->Write();
    NEvents->Write();
    outFile.Close();
    cout << "Filled: " << fillCounter << endl;
}
