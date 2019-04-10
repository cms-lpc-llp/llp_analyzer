#include "ElectronNtupler.h"
#include "ElectronTree.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

void ElectronNtupler::Analyze(bool isData, int Option, string outputfilename, string label)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "ElectronNtuple.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    ElectronTree *eleTree = new ElectronTree;
    eleTree->CreateTree(ElectronTree::kEleTreeLight);
    eleTree->tree_->SetAutoFlush(0);
    
    cout << "Run With Option = " << Option << "\n";
    
    UInt_t NElectronsFilled = 0;
 
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


      //Get NPU
      int npu = 0;
      //Get number of PU interactions
      for (int i = 0; i < nBunchXing; i++) {
	if (BunchXing[i] == 0) {
	  npu = nPUmean[i];
	}
      }

      //****************************************
      //Tree entries based on reco objects
      //****************************************
      //initialize
      eleTree->fEleTriggerBit = 0;

      if (Option < 10 ) {

        for(int i = 0; i < nElectrons; i++){
	  if(elePt[i] < 5) continue;

	  //***********************
	  //Fill Electron Variables
	  //***********************
   
	  eleTree->fWeight = 1;
	  eleTree->fRunNumber = runNum;
	  eleTree->fLumiSectionNumber = lumiNum;
	  eleTree->fEventNumber = eventNum;
	  eleTree->fEleEventNumberParity = (eventNum % 2 == 0);
	  eleTree->fCharge = eleCharge[i] ;
	  eleTree->fElePt = elePt[i]; 
	  eleTree->fEleEta = eleEta[i]; 
	  eleTree->fElePhi = elePhi[i]; 
	  eleTree->fEleSCEta = eleEta_SC[i]; 
	  eleTree->fEleTriggerBit = 0;
	  eleTree->fNPU = npu; 
	  eleTree->fRho = fixedGridRhoFastjetAll; 
	  eleTree->fRhoNeutralCentral = fixedGridRhoFastjetCentralNeutral;
	  eleTree->fNVertices = nPV; 
	  eleTree->fActivity = ele_activityMiniIsoAnnulus[i];
	  eleTree->fEleD0 = ele_d0[i]; 
	  eleTree->fEleDZ = ele_dZ[i]; 
	  eleTree->fEleIP3d = ele_ip3d[i];
	  eleTree->fEleIP3dSig = ele_ip3dSignificance[i];
	  eleTree->fElePassConversion = ele_PassConvVeto[i];
	  eleTree->fEleNMissHits =ele_MissHits[i];
	  eleTree->fEleOneOverEMinusOneOverP = ele_OneOverEminusOneOverP[i];
	  eleTree->fEleDEtaIn = ele_dEta[i];					   
	  eleTree->fEleDPhiIn = ele_dPhi[i];
	  eleTree->fEleSigmaIEtaIEta = eleFull5x5SigmaIetaIeta[i];
	  eleTree->fEleR9 = eleR9[i];
	  eleTree->fEleHoverE = ele_HoverE[i];
	  eleTree->fElePFIso04 = (ele_chargedIso[i] + fmax(0.0,  ele_photonIso[i] + ele_neutralHadIso[i] - 0.5*ele_pileupIso[i])) / elePt[i];
	  eleTree->fIDMVAGeneralPurpose = ele_IDMVAGeneralPurpose[i];
	  eleTree->fIDMVAHZZ = ele_IDMVAHZZ[i];
	  eleTree->fPassVetoSelection = isEGammaPOGVetoElectron(i);
	  eleTree->fPassLooseSelection = isLooseElectron(i);
	  eleTree->fPassTightSelection = isTightElectron(i);
	  eleTree->fPassMVAVetoSelection = isVetoElectron(i);
	  eleTree->fPtRel = ele_ptrel[i];
	  eleTree->fMiniIsoCharged = ele_chargedMiniIso[i];
	  eleTree->fMiniIsoNeutral = ele_photonAndNeutralHadronMiniIso[i];
	  eleTree->fMiniIso = ele_chargedMiniIso[i] + ele_photonAndNeutralHadronMiniIso[i];
	  eleTree->fMiniIsoDBCorr = ele_chargedMiniIso[i] + fmax(0.0, ele_photonAndNeutralHadronMiniIso[i] - 0.5*ele_chargedPileupMiniIso[i]);

	  UInt_t tmpTriggerBit = 0;
	  if (matchElectronHLTFilters(i, "Ele27Loose", "2016")) tmpTriggerBit |= ElectronTree::kEleTrigger_Ele27Loose;
	  if (matchElectronHLTFilters(i, "Ele27Tight", "2016")) tmpTriggerBit |= ElectronTree::kEleTrigger_Ele27Tight;
	  if (matchElectronHLTFilters(i, "Ele32Tight", "2016")) tmpTriggerBit |= ElectronTree::kEleTrigger_Ele32Tight;
	  if (matchElectronHLTFilters(i, "Ele105", "2016")) tmpTriggerBit |= ElectronTree::kEleTrigger_Ele105;
	  if (matchElectronHLTFilters(i, "Ele115", "2016")) tmpTriggerBit |= ElectronTree::kEleTrigger_Ele115;
	  eleTree->fEleTriggerBit = tmpTriggerBit;

	  //Match to Gen particle
	  int matchedIndex = -1;
	  float minDR = 9999;

	  for(int j = 0; j < nGenParticle; j++){
	    if (abs(gParticleId[j]) != 11) continue;	      
	    if ( deltaR( eleEta[i], elePhi[i], gParticleEta[j], gParticlePhi[j]) < 0.1
		 && deltaR( eleEta[i], elePhi[i], gParticleEta[j], gParticlePhi[j]) < minDR
		 ) {		
	      matchedIndex = j;
	      minDR = deltaR( eleEta[i], elePhi[i], gParticleEta[j], gParticlePhi[j]);
	    }
	  }

	  int matchedID = 0;
	  if (matchedIndex >= 0) {
	    eleTree->fEleGenPt = gParticlePt[matchedIndex];
	    eleTree->fEleGenEta = gParticleEta[matchedIndex];
	    eleTree->fEleGenPhi = gParticlePhi[matchedIndex];
	    if (abs(gParticleMotherId[matchedIndex]) > 50 || 
		abs(gParticleMotherId[matchedIndex]) == 15 || 
		abs(gParticleMotherId[matchedIndex]) == 13 ) {
	      matchedID = gParticleMotherId[matchedIndex];
	    } else if (abs(gParticleMotherId[matchedIndex]) == 23 || abs(gParticleMotherId[matchedIndex]) == 24 ||
		       (abs(gParticleMotherId[matchedIndex]) >= 1 && abs(gParticleMotherId[matchedIndex]) <= 5) ||
		       abs(gParticleMotherId[matchedIndex]) == 21 || abs(gParticleMotherId[matchedIndex]) == 2212		       
		       ) {
	      matchedID = gParticleId[matchedIndex];
	    }
	  }
	  eleTree->fPdgId = matchedID;

	  //select only fakes
	  if (Option == 0) {
	    if (!(matchedID == 0 || abs(matchedID) > 50)) continue;
	    //if (abs(matchedID) == 11) continue;
	  }
	  //select only real prompt
	  if (Option == 1) {
	    if (!(abs(matchedID) == 11)) continue;
	  }	   


	  //Find Closest Parton
	  float minDRToParton = 9999;
	  for(int j = 0; j < nGenParticle; j++){
	      
	    //only look for outgoing partons
	    if  (!( ((abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) || abs(gParticleId[j]) == 21) 
		    && gParticleStatus[j] == 23)
		 ) continue;
	      
	    double tmpDR = deltaR( gParticleEta[j], gParticlePhi[j], eleEta[i], elePhi[i]);
	    if ( tmpDR < minDRToParton ) minDRToParton = tmpDR;
	  }
	  eleTree->fDRToClosestParton = minDRToParton;



	  if (
	      //matchedID == 0 && 
	      
	      ele_IDMVAGeneralPurpose[i] > 0.95) {
	    std::cout << "DEBUG: \n";
	    std::cout << elePt[i] << " " << eleEta[i] << " " << elePhi[i] << " : " <<  ele_IDMVAGeneralPurpose[i] << " " <<  ele_IDMVAHZZ[i] << "\n";
	    cout << "match : " << matchedIndex << " ";
	    if (matchedIndex >= 0) cout << gParticleMotherId[matchedIndex];
	    cout << "\n";

	    if (abs(gParticleMotherId[matchedIndex]) > 50 || 
	  	  abs(gParticleMotherId[matchedIndex]) == 15 || 
	  	  abs(gParticleMotherId[matchedIndex]) == 13) {
	  	cout << "matched HF\n";
	    } else if (abs(gParticleMotherId[matchedIndex]) == 23 || abs(gParticleMotherId[matchedIndex]) == 24 ||
		       (abs(gParticleMotherId[matchedIndex]) >= 1 && abs(gParticleMotherId[matchedIndex]) <= 5) ||
		       abs(gParticleMotherId[matchedIndex]) == 21 || abs(gParticleMotherId[matchedIndex]) == 2212		       
		       ) {
	  	cout << "match prompt\n";
	    } else {
	  	cout << "nothing\n";
	    }

	    for(int j = 0; j < nGenParticle; j++){
	  	std::cout << "particle " << j << " : " << gParticleId[j] << " " << gParticleStatus[j]  << " " 
			  << gParticleMotherId[j] << " "
	  		  << gParticlePt[j] << " " <<  gParticleEta[j] << " " << gParticlePhi[j] << " " 
	  	     <<  deltaR( eleEta[i], elePhi[i], gParticleEta[j], gParticlePhi[j]) << "\n";
	    }    
	    std::cout << "\n";
	  }
	

	  //***********************
	  //Fill Electron
	  //***********************
	  NElectronsFilled++;
	  eleTree->tree_->Fill();	   
        }
      }

      //********************************************
      //Tree entries based on gen-level objects
      //********************************************
      else if (Option >= 10 && Option < 20) {

	for(int i = 0; i < nGenParticle; i++){
	    
	  //select prompt electrons
	  if (Option == 11) {
	    if (abs(gParticleId[i]) != 11) continue;	      
	    if (!(abs(gParticleMotherId[i]) == 23 || 
		  abs(gParticleMotherId[i]) == 24)
		) continue; 
	  }

	  if(gParticlePt[i] < 5) continue;

	  //***********************
	  //Fill Electron Variables
	  //***********************
   	  eleTree->fWeight = 1;
	  eleTree->fRunNumber = runNum;
	  eleTree->fLumiSectionNumber = lumiNum;
	  eleTree->fEventNumber = eventNum;
	  eleTree->fEleEventNumberParity = (eventNum % 2 == 0);
	  eleTree->fEleGenPt = gParticlePt[i];
	  eleTree->fEleGenEta = gParticleEta[i];
	  eleTree->fEleGenPhi = gParticlePhi[i];
	  eleTree->fNPU = npu; 
	  eleTree->fRho = fixedGridRhoFastjetAll; 
	  eleTree->fRhoNeutralCentral = fixedGridRhoFastjetCentralNeutral;
	  eleTree->fNVertices = nPV; 
	  eleTree->fPdgId = gParticleId[i];


	  //Find Closest Parton
	  float minDRToParton = 9999;
	  for(int j = 0; j < nGenParticle; j++){
	      
	    //only look for outgoing partons
	    if  (!( ((abs(gParticleId[j]) >= 1 && abs(gParticleId[j]) <= 5) || abs(gParticleId[j]) == 21) 
		    && gParticleStatus[j] == 23)
		 ) continue;
	      
	    double tmpDR = deltaR( gParticleEta[j], gParticlePhi[j], gParticleEta[i], gParticlePhi[i]);
	    if ( tmpDR < minDRToParton ) minDRToParton = tmpDR;
	  }
	  eleTree->fDRToClosestParton = minDRToParton;


	  //Find Associated Reco Electron
	  int matchedIndex = -1;
	  float minDR = 9999;
	  for(int j = 0; j < nElectrons; j++){
	    if ( deltaR( eleEta[j], elePhi[j], gParticleEta[i], gParticlePhi[i]) < 0.1
		 && deltaR( eleEta[j], elePhi[j], gParticleEta[i], gParticlePhi[i]) < minDR
		 ) {		
	      matchedIndex = j;
	      minDR = deltaR( eleEta[j], elePhi[j], gParticleEta[i], gParticlePhi[i]);
	    }	    
	  }	 

	  if (matchedIndex >= 0) {
	    eleTree->fActivity = ele_activityMiniIsoAnnulus[matchedIndex];
	    eleTree->fCharge = eleCharge[matchedIndex] ;
	    eleTree->fElePt = elePt[matchedIndex]; 
	    eleTree->fEleEta = eleEta[matchedIndex]; 
	    eleTree->fElePhi = elePhi[matchedIndex]; 
	    eleTree->fEleSCEta = eleEta_SC[matchedIndex]; 
	    eleTree->fEleTriggerBit = 0;
	    eleTree->fEleD0 = ele_d0[matchedIndex]; 
	    eleTree->fEleDZ = ele_dZ[matchedIndex]; 
	    eleTree->fEleIP3d = ele_ip3d[matchedIndex];
	    eleTree->fEleIP3dSig = ele_ip3dSignificance[matchedIndex];
	    eleTree->fElePassConversion = ele_PassConvVeto[matchedIndex];
	    eleTree->fEleNMissHits =ele_MissHits[matchedIndex];
	    eleTree->fEleOneOverEMinusOneOverP = ele_OneOverEminusOneOverP[matchedIndex];
	    eleTree->fEleDEtaIn = ele_dEta[matchedIndex];					   
	    eleTree->fEleDPhiIn = ele_dPhi[matchedIndex];
	    eleTree->fEleSigmaIEtaIEta = eleFull5x5SigmaIetaIeta[matchedIndex];
	    eleTree->fEleR9 = eleR9[matchedIndex];
	    eleTree->fEleHoverE = ele_HoverE[matchedIndex];
	    eleTree->fElePFIso04 = (ele_chargedIso[matchedIndex] + fmax(0.0,  ele_photonIso[matchedIndex] + ele_neutralHadIso[matchedIndex] - 0.5*ele_pileupIso[matchedIndex])) / elePt[matchedIndex];
	    eleTree->fIDMVAGeneralPurpose = ele_IDMVAGeneralPurpose[matchedIndex];
	    eleTree->fIDMVAHZZ = ele_IDMVAHZZ[matchedIndex];
	    eleTree->fPassVetoSelection = isEGammaPOGVetoElectron(matchedIndex);
	    eleTree->fPassLooseSelection = isLooseElectron(matchedIndex);
	    eleTree->fPassTightSelection = isTightElectron(matchedIndex);
	    eleTree->fPassMVAVetoSelection = isVetoElectron(matchedIndex);
	    eleTree->fPtRel = ele_ptrel[matchedIndex];
	    eleTree->fMiniIsoCharged = ele_chargedMiniIso[matchedIndex];
	    eleTree->fMiniIsoNeutral = ele_photonAndNeutralHadronMiniIso[matchedIndex];
	    eleTree->fMiniIso = ele_chargedMiniIso[matchedIndex] + ele_photonAndNeutralHadronMiniIso[matchedIndex];
	    eleTree->fMiniIsoDBCorr = ele_chargedMiniIso[matchedIndex] + fmax(0.0, ele_photonAndNeutralHadronMiniIso[matchedIndex] - 0.5*ele_chargedPileupMiniIso[matchedIndex]);
	    eleTree->fEleTriggerBit = 0;
	    if (matchElectronHLTFilters(matchedIndex, "Ele27Loose", "2016")) eleTree->fEleTriggerBit |= ElectronTree::kEleTrigger_Ele27Loose;
	    if (matchElectronHLTFilters(matchedIndex, "Ele27Tight", "2016")) eleTree->fEleTriggerBit |= ElectronTree::kEleTrigger_Ele27Tight;
	    if (matchElectronHLTFilters(matchedIndex, "Ele32Tight", "2016")) eleTree->fEleTriggerBit |= ElectronTree::kEleTrigger_Ele32Tight;
	    if (matchElectronHLTFilters(matchedIndex, "Ele105", "2016")) eleTree->fEleTriggerBit |= ElectronTree::kEleTrigger_Ele105;
	    if (matchElectronHLTFilters(matchedIndex, "Ele115", "2016")) eleTree->fEleTriggerBit |= ElectronTree::kEleTrigger_Ele115;


	  } else {
	    eleTree->fActivity = 9999;
	    eleTree->fCharge = 0;
	    eleTree->fElePt = 0;
	    eleTree->fEleEta = 0;
	    eleTree->fElePhi = 0;
	    eleTree->fEleSCEta = 0;
	    eleTree->fEleTriggerBit = 0;
	    eleTree->fEleD0 = 0;
	    eleTree->fEleDZ = 0;
	    eleTree->fElePassConversion = false;
	    eleTree->fEleNMissHits = 0;
	    eleTree->fEleOneOverEMinusOneOverP = 0;
	    eleTree->fEleDEtaIn = 0;
	    eleTree->fEleDPhiIn = 0;
	    eleTree->fEleSigmaIEtaIEta = 0;
	    eleTree->fEleR9 = 0;
	    eleTree->fEleHoverE = 0;
	    eleTree->fElePFIso04 = 0;
	    eleTree->fIDMVAGeneralPurpose = 0;
	    eleTree->fIDMVAHZZ = 0;
	    eleTree->fPassVetoSelection = false;
	    eleTree->fPassLooseSelection = false;
	    eleTree->fPassTightSelection = false;
	    eleTree->fPassMVAVetoSelection = false;
	    eleTree->fPtRel = 0;
	    eleTree->fMiniIso = 0;
	    eleTree->fEleTriggerBit = 0;
	  }
	  
	  //***********************
	  //Fill Electron
	  //***********************
	  NElectronsFilled++;
	  eleTree->tree_->Fill();
	}

      }


    }//end of event loop

    cout << "Filled Total of " << NElectronsFilled << " Electrons\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}



