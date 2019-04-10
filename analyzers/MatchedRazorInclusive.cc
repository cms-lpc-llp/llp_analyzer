//MatchedRazorInclusive implements the inclusive razor selection, computes MR from all jets matching in DeltaR with gen-level jets, and partitions the selected events into boxes according to lepton and jet content.

#include "MatchedRazorInclusive.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

void MatchedRazorInclusive::Analyze(bool isData, int option, string outFileName, string label)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    bool combineTrees = false;
    if (outFileName.empty()){
        cout << "MatchedRazorInclusive: Output filename not specified!" << endl << "Using default output name MatchedRazorInclusive.root" << endl;
        outFileName = "MatchedRazorInclusive.root";
    }
    TFile outFile(outFileName.c_str(), "RECREATE");
    
    //one tree to hold all events
    TTree *razorTree = new TTree("MatchedRazorInclusive", "Info on selected razor inclusive events");
    
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
    boxNames.push_back("MultiJet");
    boxNames.push_back("TwoBJet");
    boxNames.push_back("OneBJet");
    boxNames.push_back("ZeroBJet");
    for(size_t i = 0; i < boxNames.size(); i++){
        razorBoxes[boxNames[i]] = new TTree(boxNames[i].c_str(), boxNames[i].c_str());
    }

    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

    //tree variables
    int nSelectedJets, nBTaggedJets, nSelectedGenJets;
    int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
    float theMR;
    float theRsq;
    RazorBox box;

    //set branches on big tree
    if(combineTrees){
        razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
        razorTree->Branch("nSelectedGenJets", &nSelectedGenJets, "nSelectedGenJets/I");
        razorTree->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
        razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
        razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
        razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
        razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
        razorTree->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
        razorTree->Branch("MR", &theMR, "MR/F");
        razorTree->Branch("Rsq", &theRsq, "Rsq/F");
        razorTree->Branch("box", &box, "box/I");
    }
    //set branches on all trees
    else{ 
        for(auto& box : razorBoxes){
            box.second->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
            box.second->Branch("nSelectedGenJets", &nSelectedGenJets, "nSelectedGenJets/I");
            box.second->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
            box.second->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
            box.second->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
            box.second->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
            box.second->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
            box.second->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
            box.second->Branch("MR", &theMR, "MR/F");
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
        nSelectedGenJets = 0;
        nBTaggedJets = 0;
        nLooseMuons = 0;
        nTightMuons = 0;
        nLooseElectrons = 0;
        nTightElectrons = 0;
        nTightTaus = 0;
        theMR = -1;
        theRsq = -1;
        if(combineTrees) box = NONE;

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

            nTightTaus++;
        }

        //select genjets
        vector<TLorentzVector> GoodGenJets;
        for(int i = 0; i < nGenJets; i++){
            if(genJetPt[i] < 40) continue;
            if(fabs(genJetEta[i]) > 3.0) continue;
            TLorentzVector thisGenJet = makeTLorentzVector(genJetPt[i], genJetEta[i], genJetPhi[i], genJetE[i]);
            GoodGenJets.push_back(thisGenJet);
            nSelectedGenJets++;
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

            //exclude jet if it doesn't match a selected genjet
            bool matchedGenJet = false;
            for(auto& genjet : GoodGenJets){
                double thisDR = thisJet.DeltaR(genjet);
                if(thisDR < 0.4){
                    matchedGenJet = true;
                    break;
                }
            }
            if(!matchedGenJet) continue;
            
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
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = MuEle;
                    razorTree->Fill();
                }
                else razorBoxes["MuEle"]->Fill();
            }
        }
        //MuMu Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nLooseMuons > 1 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = MuMu;
                    razorTree->Fill();
                }
                else razorBoxes["MuMu"]->Fill();
            }
        }
        //EleEle Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nLooseElectrons > 1 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = EleEle;
                    razorTree->Fill();
                }
                else razorBoxes["EleEle"]->Fill();
            }
        }
        //MuMultiJet Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = MuMultiJet;
                    razorTree->Fill();
                }
                else razorBoxes["MuMultiJet"]->Fill();
            }
        }
        //MuJet Box
        else if(passedLeptonicTrigger && nTightMuons > 0 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = MuJet;
                    razorTree->Fill();
                }
                else razorBoxes["MuJet"]->Fill();
            }
        }
        //EleMultiJet Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nBTaggedJets > 0 && nSelectedJets > 3){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = EleMultiJet;
                    razorTree->Fill();
                }
                else razorBoxes["EleMultiJet"]->Fill();
            }
        }
        //EleJet Box
        else if(passedLeptonicTrigger && nTightElectrons > 0 && nBTaggedJets > 0){
            if(passesLeptonicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = EleJet;
                    razorTree->Fill();
                }
                else razorBoxes["EleJet"]->Fill();
            }
        }
        //MultiJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 0 && nSelectedJets > 3){
            if(passesHadronicRazorBaseline(theMR, theRsq)){  
                if(combineTrees){
                    box = MultiJet;
                    razorTree->Fill();
                }
                else { razorBoxes["MultiJet"]->Fill(); }
            }
        }
        //TwoBJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 1){
            if(passesHadronicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = TwoBJet;
                    razorTree->Fill();
                }
                else razorBoxes["TwoBJet"]->Fill();
            }
        }
        //OneBJet Box
        else if(passedHadronicTrigger && nBTaggedJets > 0){
            if(passesHadronicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = OneBJet;
                    razorTree->Fill();
                }
                else razorBoxes["OneBJet"]->Fill();
            }
        }
        //ZeroBJetBox
        else if(passedHadronicTrigger){
            if(passesHadronicRazorBaseline(theMR, theRsq)){ 
                if(combineTrees){
                    box = ZeroBJet;
                    razorTree->Fill();
                }
                else razorBoxes["ZeroBJet"]->Fill();
            }
        }
    }//end of event loop

    cout << "Writing output trees..." << endl;
    if(combineTrees) razorTree->Write();
    else for(auto& box : razorBoxes) box.second->Write();
    NEvents->Write();

    outFile.Close();
}


