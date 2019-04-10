#include "RazorMetAna.h"

//C++ includes

//ROOT includes
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

bool passesHadronicRazorBaselineCopy(double MR, double Rsq); //makes sense to put these into RazorAnalyzer
bool passesLeptonicRazorBaselineCopy(double MR, double Rsq);

void RazorMetAna::Analyze(bool isData, int option, string outFileName, string label)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    string outfilename = outFileName;
    if (outFileName == "") outfilename = "RazorMET.root";
    TFile outFile(outFileName.c_str(), "RECREATE");
    
    //one tree to hold all events
    TTree *razorTree = new TTree("RazorMet", "Info on selected razor inclusive events");
    
    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

    // histograms
    //TH1F *h_CaloMET     = new TH1F ("h_CaloMET","h_CaloMET",100,0,200);

    TH1F *h_MR_MultiJet = new TH1F ("h_MR_MultiJet","h_MR_MultiJet",100,0,4000);
    TH1F *h_Rsq_MultiJet = new TH1F ("h_Rsq_MultiJet","h_Rsq_MultiJet",50,0,1.5);

    TH1F *h_MR_MuEle = new TH1F ("h_MR_MuEle","h_MR_MuEle",100,0,4000);
    TH1F *h_Rsq_MuEle = new TH1F ("h_Rsq_MuEle","h_Rsq_MuEle",50,0,1.5);

    TH1F *h_MR_MuMu = new TH1F ("h_MR_MuMu","h_MR_MuMu",100,0,4000);
    TH1F *h_Rsq_MuMu = new TH1F ("h_Rsq_MuMu","h_Rsq_MuMu",50,0,1.5);

    TH1F *h_MR_EleEle = new TH1F ("h_MR_EleEle","h_MR_EleEle",100,0,4000);
    TH1F *h_Rsq_EleEle = new TH1F ("h_Rsq_EleEle","h_Rsq_EleEle",50,0,1.5);

    TH1F *h_MR_MuMultiJet = new TH1F ("h_MR_MuMultiJet","h_MR_MuMultiJet",100,0,4000);
    TH1F *h_Rsq_MuMultiJet = new TH1F ("h_Rsq_MuMultiJet","h_Rsq_MuMultiJet",50,0,1.5);

    TH1F *h_MR_MuJet = new TH1F ("h_MR_MuJet","h_MR_MuJet",100,0,4000);
    TH1F *h_Rsq_MuJet = new TH1F ("h_Rsq_MuJet","h_Rsq_MuJet",50,0,1.5);

    TH1F *h_MR_EleMultiJet = new TH1F ("h_MR_EleMultiJet","h_MR_EleMultiJet",100,0,4000);
    TH1F *h_Rsq_EleMultiJet = new TH1F ("h_Rsq_EleMultiJet","h_Rsq_EleMultiJet",50,0,1.5);

    TH1F *h_MR_EleJet = new TH1F ("h_MR_EleJet","h_MR_EleJet",100,0,4000);
    TH1F *h_Rsq_EleJet = new TH1F ("h_Rsq_EleJet","h_Rsq_EleJet",50,0,1.5);

    TH1F *h_ZMass = new TH1F ("h_ZMass","h_ZMass",50,70,120);
    TH1F *h_upar = new TH1F ("h_upar","h_upar",50,-100,100);
    TH1F *h_uperp = new TH1F ("h_uperp","h_uperp",50,-100,100);
    TH2F *h_upar_NVtx = new TH2F ("h_upar_NVtx","h_upar_NVtx",25,0,50, 100, -50,50);
    TH2F *h_uperp_NVtx = new TH2F ("h_uperp_NVtx","h_uperp_NVtx",25,0,50, 100, -50,50);

    //tree variables
    int nSelectedJets, nBTaggedJets;
    int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nSelectedTaus;
    float theMR;
    float theRsq;

    razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");

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
        nLooseMuons = 0;
        nTightMuons = 0;
        nLooseElectrons = 0;
        nTightElectrons = 0;
        nSelectedTaus = 0;
        theMR = -1;
        theRsq = -1;

        //TODO: triggers!
        bool passedLeptonicTrigger = true;
        bool passedHadronicTrigger= true;
        if(!(passedLeptonicTrigger || passedHadronicTrigger)) continue; //ensure event passed a trigger

        vector<TLorentzVector> GoodLeptons; //leptons used to compute hemispheres
        for(int i = 0; i < nMuons; i++){
            if(!isLooseMuon(i)) continue;  
            if(muonPt[i] < 10) continue;
            if(abs(muonEta[i]) > 2.4) continue;

            nLooseMuons++;
            TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
            GoodLeptons.push_back(thisMuon);

            if(isTightMuon(i)){ 
                nTightMuons++;
            }
        }
        for(int i = 0; i < nElectrons; i++){
            if(!isLooseElectron(i)) continue; 
            if(elePt[i] < 10) continue;

            nLooseElectrons++;
            TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
            GoodLeptons.push_back(thisElectron);

            if(isTightElectron(i)){ 
                nTightElectrons++;
            }
        }
        for(int i = 0; i < nTaus; i++){
            if(!isTightTau(i)) continue; 

            nSelectedTaus++;
        }
	
        vector<TLorentzVector> GoodJets;
        int numJetsAbove80GeV = 0;
        for(int i = 0; i < nJets; i++){
            if(jetPt[i] < 40) continue;
            if(fabs(jetEta[i]) > 3.0) continue;

            //exclude selected muons and electrons from the jet collection
            double deltaR = -1;
            TLorentzVector thisJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);
            for(auto& lep : GoodLeptons){
                double thisDR = thisJet.DeltaR(lep);
                if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
            }
            if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton
            
            if(jetPt[i] > 80) numJetsAbove80GeV++;
            GoodJets.push_back(thisJet);
            nSelectedJets++;

            if(isCSVM(i)){ 
                nBTaggedJets++;
            }
        }

        if(numJetsAbove80GeV < 2) continue; //event fails to have two 80 GeV jets

        //Compute the razor variables using the selected jets and possibly leptons
        vector<TLorentzVector> GoodPFObjects;
        for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
        if(passedLeptonicTrigger) for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);
        TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);

        vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
        theMR = computeMR(hemispheres[0], hemispheres[1]); 
        theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);

        //MuEle Box
        if(passedLeptonicTrigger && nTightElectrons > 0 && nLooseMuons > 0 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaselineCopy(theMR, theRsq)){ 
		if(theRsq > 0.25 && theMR > 1000) 
		  h_MR_MuEle->Fill(theMR);
		if(theRsq > 0.10 && theMR > 1500) 
		    h_Rsq_MuEle->Fill(theRsq);
            }
        }
        //MuMu Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nLooseMuons > 1 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaselineCopy(theMR, theRsq)){ 
		if(theRsq > 0.25 && theMR > 1000) 
		    h_MR_MuMu->Fill(theMR);
		if(theRsq > 0.10 && theMR > 1500) 
		    h_Rsq_MuMu->Fill(theRsq);
            }
        }
        //EleEle Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nLooseElectrons > 1 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaselineCopy(theMR, theRsq)){ 
		if(theRsq > 0.25 && theMR > 1000) 
		    h_MR_EleEle->Fill(theMR);
		if(theRsq > 0.10 && theMR > 1500) 
		    h_Rsq_EleEle->Fill(theRsq);
            }
        }
        //MuMultiJet Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(passesLeptonicRazorBaselineCopy(theMR, theRsq)){ 
		if(theRsq > 0.25 && theMR > 1000) 
		    h_MR_MuMultiJet->Fill(theMR);
		if(theRsq > 0.10 && theMR > 1500) 
		    h_Rsq_MuMultiJet->Fill(theRsq);
            }
        }
        //MuJet Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaselineCopy(theMR, theRsq)){ 
		if(theRsq > 0.25 && theMR > 1000) 
		    h_MR_MuJet->Fill(theMR);
		if(theRsq > 0.10 && theMR > 1500) 
		    h_Rsq_MuJet->Fill(theRsq);
            }
        }
        //EleMultiJet Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(passesLeptonicRazorBaselineCopy(theMR, theRsq)){ 
		if(theRsq > 0.25 && theMR > 1000) 
		    h_MR_EleMultiJet->Fill(theMR);
		if(theRsq > 0.10 && theMR > 1500) 
		    h_Rsq_EleMultiJet->Fill(theRsq);
            }
        }
        //EleJet Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaselineCopy(theMR, theRsq)){ 
		if(theRsq > 0.25 && theMR > 1000) 
		    h_MR_EleJet->Fill(theMR);
		if(theRsq > 0.10 && theMR > 1500) 
		    h_Rsq_EleJet->Fill(theRsq);
            }
        }
        //MultiJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 0 && nSelectedJets > 3){
	    if(passesHadronicRazorBaselineCopy(theMR, theRsq)){  
		if(theRsq > 0.25 && theMR > 1000) 
		    h_MR_MultiJet->Fill(theMR);
		if(theRsq > 0.10 && theMR > 1500) 
		    h_Rsq_MultiJet->Fill(theRsq);
	      }
	  }

	// MET resolution stuff
        vector<TLorentzVector> ZMuons; //leptons used to compute hemispheres
        for(int i = 0; i < nMuons; i++){
	  if(!isLooseMuon(i)) continue;  
	  if(muonPt[i] < 10) continue;
	  if(abs(muonEta[i]) > 2.4) continue;
	    
	  TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
	  ZMuons.push_back(thisMuon);	  
        }
	
	// if(ZMuons.size()==2)
	//   {
	//     TLorentzVector m1 = makeTLorentzVector(ZMuons[0].Pt(), ZMuons[0].Eta(), ZMuons[0].Phi(), ZMuons[0].E()); 
	//     TLorentzVector m2 = makeTLorentzVector(ZMuons[1].Pt(), ZMuons[1].Eta(), ZMuons[1].Phi(), ZMuons[1].E()); 

	//     TLorentzVector theZ = m1+m2;
	//     h_ZMass->Fill(theZ.M());

	//     TLorentzVector qT = theZ.Perp();
	//     TLorentzVector met = makeTLorentzVector(metPt, 0, metPhi, 0);
	//     TLorentzVector uT = -(qT + met);
	    
	//     Float_t u_par   = - (uT.Vect()).Dot(qT.Vect().Unit());
	//     TVector3 u_Vec_perp = - (uT.Vect()).Cross(qT.Vect().Unit());

	//     Float_t sign = TMath::Sin(qT.Vect().DeltaPhi(uT.Vect()));
	//     Float_t u_perp = (sign>0) ? u_Vec_perp.Mag() : - u_Vec_perp.Mag();

	//     h_upar->Fill(u_par+qT.Mag());
	//     h_uperp->Fill(u_perp);

	//     h_upar_NVtx->Fill(nPV, u_par+qT.Mag());
	//     h_uperp_NVtx->Fill(nPV, u_perp);
	//   }

    }//end of event loop

    cout << "Writing output trees..." << endl;
    razorTree->Write();   

    // Write histograms
    NEvents->Write();

    h_MR_MultiJet->Write();
    h_Rsq_MultiJet->Write();
    
    h_MR_MuEle->Write();
    h_Rsq_MuEle->Write();

    h_MR_MuMu->Write();
    h_Rsq_MuMu->Write();

    h_MR_EleEle->Write();
    h_Rsq_EleEle->Write();

    h_MR_MuMultiJet->Write();
    h_Rsq_MuMultiJet->Write();

    h_MR_MuJet->Write();
    h_Rsq_MuJet->Write();

    h_MR_EleMultiJet->Write();
    h_Rsq_EleMultiJet->Write();

    h_MR_EleJet->Write();
    h_Rsq_EleJet->Write();

    h_ZMass->Write();
    h_upar->Write();
    h_upar_NVtx->Write();
    h_uperp->Write();
    h_uperp_NVtx->Write();

    outFile.Close();
}

bool passesHadronicRazorBaselineCopy(double MR, double Rsq){
    bool passes = true;
    if(MR < 0 || Rsq < 0) passes = false;
    //temporarily disable these
    //if(MR < 400 || Rsq < 0.25) passes = false;
    //if(MR < 450 && Rsq < 0.3) passes = false;
    return passes;
}

bool passesLeptonicRazorBaselineCopy(double MR, double Rsq){
    bool passes = true;
    if(MR < 300 || Rsq < 0.15) passes = false;
    if(MR < 350 && Rsq < 0.2) passes = false;
    return passes;
}
