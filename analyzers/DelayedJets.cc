#include "DelayedJets.h"
#include "JetTree.h"
#include "JetCorrectorParameters.h"

//C++ includes

//ROOT includes
#include "TH1F.h"


struct MatchedPhoton
{
  TLorentzVector pho;
  float SeedRecHitE;
  float SeedRecHitEta;
  float SeedRecHitPhi;
  float SeedRecHitTime;
};

struct MatchedRecHit
{
  float E;
  float Eta;
  float Phi;
  float Time;
};

using namespace std;

void DelayedJets::Analyze(bool isData, int Option, string outputfilename, string label)
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

  std::cout << "Initializing..." << std::endl;
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "DelayedJets.root";
  TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
  JetTree *jetTree = new JetTree;
  jetTree->CreateTree();
  jetTree->tree_->SetAutoFlush(0);
  jetTree->InitVariables();

  std::cout << "Run With Option = " << Option << std::endl;

  UInt_t NJetsFilled = 0;

  //begin loop
  if (fChain == 0) return;
  UInt_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  std::cout << "nentries = " << nentries << std::endl;
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

  jetTree->ResetVariables();
	jetTree->fJetNPhotons = 0;
  jetTree->fJetNRecHits = 0;
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


	//Match to Photons
	int matchedIndex = -1;
	float minDR = 9999;
	std::vector<MatchedPhoton> matchedPhotons;
  std::vector<MatchedRecHit> matchedRecHits;

	for( int j = 0; j < nPhotons; j++ )
	  {
	    double tmpDR = deltaR( phoEta[j], phoPhi[j], jetEta[i],jetPhi[i]);
	    if ( tmpDR < 0.4 )
      {
        MatchedPhoton pho;
        TVector3 vec;
        vec.SetPtEtaPhi( phoPt[j], phoEta[j], phoPhi[j] );
        TLorentzVector thisPhoton;
        thisPhoton.SetVectM( vec, .0 );
        pho.pho  = thisPhoton;
				//std::cout << jentry << "nPhotons: "  << nPhotons <<  " j: "  << j << " index: "<< pho_SeedRechitIndex->at(j)  << std::endl;
				if ( ecalRechit_T->size() != 0 )
        {
          pho.SeedRecHitE = ecalRechit_E->at(pho_SeedRechitIndex->at(j));
          pho.SeedRecHitEta = ecalRechit_Eta->at(pho_SeedRechitIndex->at(j));
          pho.SeedRecHitPhi = ecalRechit_Phi->at(pho_SeedRechitIndex->at(j));
          pho.SeedRecHitTime = ecalRechit_T->at(pho_SeedRechitIndex->at(j));
          //RecHits
          MatchedRecHit recHit;
          for ( uint it_rechit = 0; it_rechit < pho_EcalRechitIndex->at(j).size(); it_rechit++ )
          {
            recHit.E    = ecalRechit_E->at(pho_EcalRechitIndex->at(j).at(it_rechit));
            recHit.Eta  = ecalRechit_Eta->at(pho_EcalRechitIndex->at(j).at(it_rechit));
            recHit.Phi  = ecalRechit_Phi->at(pho_EcalRechitIndex->at(j).at(it_rechit));
            recHit.Time = ecalRechit_T->at(pho_EcalRechitIndex->at(j).at(it_rechit));
            matchedRecHits.push_back(recHit);
            jetTree->fJetNRecHits++;
          }
        }
        std::cout << pho_EcalRechitIndex->at(j).at(0) << std::endl;
        matchedPhotons.push_back(pho);
        jetTree->fJetNPhotons++;
        //matchedIndex = j;
        //minDR = tmpDR;
      }
    }

	int ctr = 0;
	for ( auto tmp : matchedPhotons)
	  {
	    jetTree->fJetPhotonE[ctr] = tmp.pho.E();
	    jetTree->fJetPhotonEta[ctr] = tmp.pho.Eta();
	    jetTree->fJetPhotonPhi[ctr] = tmp.pho.Phi();
      jetTree->fJetPhotonSeedRecHitE[ctr] = tmp.SeedRecHitE;
      jetTree->fJetPhotonSeedRecHitEta[ctr] = tmp.SeedRecHitEta;
      jetTree->fJetPhotonSeedRecHitPhi[ctr] = tmp.SeedRecHitPhi;
	    jetTree->fJetPhotonSeedRecHitTime[ctr] = tmp.SeedRecHitTime;
	    ctr++;
	  }

  int ctr_rechit = 0;
  for ( auto tmp : matchedRecHits)
  {
    jetTree->fJetRecHitE[ctr_rechit]    = tmp.E;
    jetTree->fJetRecHitEta[ctr_rechit]  = tmp.Eta;
    jetTree->fJetRecHitPhi[ctr_rechit]  = tmp.Phi;
    jetTree->fJetRecHitTime[ctr_rechit] = tmp.Time;
    ctr_rechit++;
  }

  if ( matchedIndex >= 0 )
	  {
	    std::cout << "found photon candidate, candidate time: " << pho_SeedRechitIndex->at(matchedIndex) << std::endl;
	  }
	else
	  {

	  }


	/*
	//Match to Gen Jet
	matchedIndex = -1;
	minDR = 9999;

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
	*/

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
  jetTree->fi_evt++;
}
}
  }//end of event loop

  cout << "Filled Total of " << NJetsFilled << " Jets\n";
  cout << "Writing output trees..." << endl;
  outFile->Write();
  outFile->Close();

};
