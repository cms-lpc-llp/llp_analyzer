#include "JetNtupler.h"
#include "JetTree.h"
#include "JetCorrectorParameters.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

void JetNtupler::Analyze(bool isData, int Option, string outputfilename, string label)
{
    //initialization: create one TTree for each analysis box 
  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/RazorAnalyzer/data/JEC/";
  cout << "Getting JEC parameters from " << pathname << endl;
  
  std::vector<JetCorrectorParameters> correctionParameters;
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));
   
  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);
  
    cout << "Initializing..." << endl;
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "JetNtuple.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    JetTree *jetTree = new JetTree;
    jetTree->CreateTree();
    jetTree->tree_->SetAutoFlush(0);
    
    cout << "Run With Option = " << Option << "\n";
    
    UInt_t NJetsFilled = 0;
 
    //begin loop
    if (fChain == 0) return;
    UInt_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;

    cout << "nentries = " << nentries << "\n";
    for (UInt_t jentry=0; jentry<nentries;jentry++) {
      //begin event
      if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //****************************************
      //Tree entries based on reco objects
      //****************************************
      if (Option < 10 ) {

	for(int i = 0; i < nJets; i++){

          double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
						 fixedGridRhoFastjetAll, jetJetArea[i], 
						 JetCorrector);   
	  
	  if(jetPt[i]*JEC < 20) continue;


	  //***********************
	  //Fill Jet Variables
	  //***********************
   
	  jetTree->fWeight = 1;
	  jetTree->fRunNumber = runNum;
	  jetTree->fLumiSectionNumber = lumiNum;
	  jetTree->fEventNumber = eventNum;
	  jetTree->fJetEventNumberParity = (eventNum % 2 == 0);
	  jetTree->fJetRawPt = jetPt[i];
	  jetTree->fJetPt = jetPt[i]*JEC; 
	  jetTree->fJetEta = jetEta[i]; 
	  jetTree->fJetPhi = jetPhi[i]; 
	  jetTree->fRho = fixedGridRhoFastjetAll; 
	  jetTree->fNVertices = nPV; 
	  jetTree->fJetCSV = jetCSV[i];
	  jetTree->fJetCISV = jetCISV[i];
	  jetTree->fJetArea = jetJetArea[i];
	  jetTree->fJetPartonFlavor = jetPartonFlavor[i];
	  jetTree->fJetPileupId = jetPileupId[i];
	  jetTree->fJetPartonFlavor = jetPartonFlavor[i];
	  jetTree->fJetJEC = JEC;

	  //Match to Gen Jet
	  int matchedIndex = -1;
	  float minDR = 9999;

	  for(int j = 0; j < nGenJets; j++){
	    double tmpDR = deltaR( genJetEta[j],genJetPhi[j], jetEta[i],jetPhi[i]);
	    if ( tmpDR < 0.4
		 && tmpDR < minDR
		 ) {		
	      matchedIndex = j;
	      minDR = tmpDR;
	    }
	  }

	  if (matchedIndex >= 0) {
	    jetTree->fJetGenPt = genJetPt[matchedIndex];
	    jetTree->fJetGenEta = genJetEta[matchedIndex];
	    jetTree->fJetGenPhi = genJetPhi[matchedIndex];	    
	  } else {
	    jetTree->fJetGenPt = -1;
	    jetTree->fJetGenEta = 0;
	    jetTree->fJetGenPhi = 0;
	  }

	 
	  //select only fakes
	  if (Option == 0) {
	    if (matchedIndex != 0) continue;
	  }
	  if (Option == 1) {
	    if (!(abs(jetPartonFlavor[i]) >= 1 && abs(jetPartonFlavor[i]) <= 3)) continue;
	  }
	  if (Option == 5) {
	    if (abs(jetPartonFlavor[i]) != 5) continue;
	  }
	  if (Option == 4) {
	    if (abs(jetPartonFlavor[i]) != 4) continue;
	  }
	  if (Option == 21) {
	    if (abs(jetPartonFlavor[i]) != 21) continue;
	  }	   


	  // //Find Closest Parton
	  // float minDRToParton = 9999;
	  // for(int j = 0; j < nGenParticle; j++){
	      
	  //   //only look for outgoing partons
	  //   if  (!( ((abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) || abs(gParticleId[j]) == 21) 
	  // 	    && gParticleStatus[j] == 23)
	  // 	 ) continue;
	      
	  //   double tmpDR = deltaR( gParticleEta[j], gParticlePhi[j], jetEta[i], jetPhi[i]);
	  //   if ( tmpDR < minDRToParton ) minDRToParton = tmpDR;
	  // }
	  


	  //***********************
	  //Fill Jet
	  //***********************
	  NJetsFilled++;
	  jetTree->tree_->Fill();
	}
      }
    }//end of event loop

    cout << "Filled Total of " << NJetsFilled << " Jets\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}



