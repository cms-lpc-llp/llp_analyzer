#include "TauNtupler.h"
#include "TauTree.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

void TauNtupler::Analyze(bool isData, int Option, string outputfilename, string label)
{
  //initialization: create one TTree for each analysis box 
  std::cout << "Initializing..." << std::endl;
  std::string outfilename = outputfilename;
  if (outfilename == "") outfilename = "TauNtuple.root";
  TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
  TauTree *tauTree = new TauTree;
  tauTree->CreateTree(TauTree::TauTreeVersion::kTauTreeLight);
  tauTree->tree_->SetAutoFlush(0);
  tauTree->SetMinPt(20.0);
  tauTree->SetMaxEta(2.4);
    
    std::cout << "Run With Option = " << Option << "\n";
    
    UInt_t NTausFilled = 0;
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
      
      
      //****************************************
      //Tree entries based on reco objects
      //****************************************
      if (Option < 10 ) {
	
	for(int i = 0; i < nTaus; i++){
	  if(tauPt[i] < tauTree->GetMinPt()) continue;
	  
	  //***********************
	  //Fill Tau Variables
	  //***********************
   
	  tauTree->fWeight = 1;
	  tauTree->fRunNumber = runNum;
	  tauTree->fLumiSectionNumber = lumiNum;
	  tauTree->fEventNumber = eventNum;
	  tauTree->fTauEventNumberParity = (eventNum % 2 == 0);
	  tauTree->fTauPt = tauPt[i]; 
	  tauTree->fTauEta = tauEta[i]; 
	  tauTree->fTauPhi = tauPhi[i]; 
	  tauTree->fRho = 0; 
	  tauTree->fNVertices = nPV; 
	  
	  tauTree->fTauIsLoose = isLooseTau(i);
	  tauTree->fTauIsTight = isTightTau(i);
	  tauTree->fPassLooseSelection = isLooseTau(i);
	  tauTree->fPassTightSelection = isTightTau(i);

	  //Match to Gen particle
	  int matchedIndex = findClosestGenTau(tauEta[i], tauPhi[i]);
	  if (matchedIndex >= 0) {
            tauTree->fTauGenPt = gParticlePt[matchedIndex];
            tauTree->fTauGenEta = gParticleEta[matchedIndex];
            tauTree->fTauGenPhi = gParticlePhi[matchedIndex];
	  }
	  //Find matchedID
	  int matchedID = GetTauMatchedID(tauEta[i], tauPhi[i]);
	  tauTree->fPdgId = matchedID;

	  //select only fakes
	  if (Option == 0) {
	    //if (!(matchedID == 0 || abs(matchedID) > 50)) continue;
	    if (abs(matchedID) == 15) continue;
	  }
	  //select only real prompt
	  if (Option == 1) {
	    if (!(abs(matchedID) == 15)) continue;
	  }	   

	  //Find Closest Parton
	  int partonIndex = findClosestParton(tauEta[i], tauPhi[i]);
	  tauTree->fDRToClosestParton = deltaR(tauEta[i], tauPhi[i], 
					       gParticleEta[partonIndex], gParticlePhi[partonIndex]);
	  
	  //***********************
	  //Fill Tau
	  //***********************
	  NTausFilled++;
	  tauTree->tree_->Fill();
	}
      }
      
      
      //********************************************
      //Tree entries based on gen-level objects
      //********************************************
      else if (Option >= 10 && Option < 20) {
	
	for(int i = 0; i < nGenParticle; i++){
	  //select prompt hadronic taus
	  if (Option == 11) {
	    if (!( isGenTau(i) && !isGenLeptonicTau(i))) continue;
	  }
	  
	  if(gParticlePt[i] < tauTree->GetMinPt()) continue;
	  //***********************
	  //Fill Tau Variables
	  //***********************
	  
	  tauTree->fWeight = 1;
	  tauTree->fRunNumber = runNum;
	  tauTree->fLumiSectionNumber = lumiNum;
	  tauTree->fEventNumber = eventNum;
	  tauTree->fTauEventNumberParity = (eventNum % 2 == 0);
	  tauTree->fTauGenPt = gParticlePt[i];
	  tauTree->fTauGenEta = gParticleEta[i];
	  tauTree->fTauGenPhi = gParticlePhi[i];
	  tauTree->fRho = 0; 
	  tauTree->fNVertices = nPV; 
	  tauTree->fPdgId = gParticleId[i];

	  //Find Closest Parton
	  int partonIndex = findClosestParton(gParticleEta[i], gParticlePhi[i]);
          tauTree->fDRToClosestParton = deltaR(gParticleEta[i], gParticlePhi[i],
                                               gParticleEta[partonIndex], gParticlePhi[partonIndex]);

	  //Find Associated Reco Tau
	  int matchedIndex = findClosestRecoTau(gParticleEta[i], gParticlePhi[i]);
	  
	  if (matchedIndex >= 0) {
	    tauTree->fTauPt = tauPt[matchedIndex]; 
	    tauTree->fTauEta = tauEta[matchedIndex]; 
	    tauTree->fTauPhi = tauPhi[matchedIndex]; 
	    tauTree->fTauIsLoose = isLooseTau(matchedIndex);
	    tauTree->fTauIsTight = isTightTau(matchedIndex);
	    tauTree->fPassLooseSelection = isLooseTau(matchedIndex);
	    tauTree->fPassTightSelection = isTightTau(matchedIndex);
	  } else {
	    tauTree->fTauPt = 0;
	    tauTree->fTauEta = 0;
	    tauTree->fTauPhi = 0;
	    tauTree->fTauIsLoose = false;
	    tauTree->fTauIsTight = false;
	    tauTree->fPassLooseSelection = false;
	    tauTree->fPassTightSelection = false;
	  }
	  
	  //***********************
	  //Fill Tau
	  //***********************
	  NTausFilled++;
	  tauTree->tree_->Fill();


	}
	  
      }
	
    }//end of event loop

    cout << "Filled Total of " << NTausFilled << " Taus\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}



