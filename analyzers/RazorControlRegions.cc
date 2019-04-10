#include "RazorControlRegions.h"
#include "JetCorrectorParameters.h"
#include "ControlSampleEvents.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "RazorHelper.h"

//C++ includes
#include <math.h>
#include <assert.h>

//ROOT includes
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>

using namespace std;

const int NUM_PDF_WEIGHTS = 60;

const int JET_PT_CUT = 30;
const int BJET_PT_CUT = 30;
const float JET_ETA_CUT = 2.4;

struct greater_than_pt{
    inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){
        return p1.Pt() > p2.Pt();
    }
};
 
// Usage of input parameters:
// -option: used to select output tree type.  See table below for available options.
// -label: used to specify the analysis tag (Razor2015, Razor2016, etc)
void RazorControlRegions::Analyze(bool isData, int option, string outputfilename, string label)
{
    cout << "Initializing..." << endl;
    cout << "IsData = " << isData << "\n";

    string tag = "Razor2016_MoriondRereco";
    if ( label != "" ) {
        tag = label;
    }

    // initialize helper
    RazorHelper helper(tag, isData, false);

    // retrieve JEC tools
    std::vector<FactorizedJetCorrector*> JetCorrector = helper.getJetCorrector();
    std::vector<std::pair<int,int> > JetCorrectorIOV = helper.getJetCorrectionsIOV();
    std::vector<SimpleJetResolution*> JetResolutionCalculator = helper.getJetResolutionCalculator();

    TRandom3 *random = new TRandom3(33333); 
    bool printSyncDebug = false;

    //*************************************************************************
    //Set up Output File
    //*************************************************************************
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "RazorControlRegions.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    ControlSampleEvents *events = new ControlSampleEvents;
    
    //**********************************************
    //Options
    //tens and ones digits refer to tree type.  
    //1: Single-Lepton
    //2: Single-Lepton Add To MET
    //3: Dilepton
    //4: Dilepton Add To MET
    //5: Photon Add To MET
    //6: Zero Lepton
    //7: Single Veto-Lepton
    //8: Tight-Lepton + Veto-Lepton
    //9: Single Tau
    //11: Single-Lepton Reduced
    //12: Single-Lepton Add To MET Reduced
    //13: Dilepton Reduced
    //14: Dilepton Add To MET Reduced
    //15: Photon Reduced
    //16: Photon Reduced
    //17: Single Veto-Lepton Reduced
    //18: Tight Lepton + Veto-Lepton Reduced
    //19: Single Tau Reduced
    //hundreds digit refers to lepton skim option
    //0: no skim
    //1: 1 lepton skim (pT > 30 GeV)
    //2: 2 lepton skim
    //3: 1 lepton skim (pT > 0 GeV)
    //5: photon skim
    //thousand digit refers to razor skim option
    //0: no skim
    //1: razor skim MR > 300 Rsq > 0.15
    //2: razor skim MR_NoW > 300 && Rsq_NoW > 0.15
    //3: razor skim MR_NoZ > 300 && Rsq_NoZ > 0.15
    //5: razor skim MR_NoPho > 300 && Rsq_NoPho > 0.15
    //6: razorDM skim MR > 150 && Rsq > 0.4
    //7: razorDM skim MR_NoW > 150 && Rsq_NoW > 0.4
    //8: razorDM skim MR_NoZ > 150 && Rsq_NoZ > 0.4
    //9: razorDM skim MR_NoPho > 150 && Rsq_NoPho > 0.4
    //*********************************************
    int razorSkimOption = floor(float(option) / 1000);
    int leptonSkimOption = floor( float(option - razorSkimOption*1000) / 100);
    int treeTypeOption = option - razorSkimOption*1000 - leptonSkimOption*100;

    cout<<"Info: razorSkimOption: "<<razorSkimOption<<", "<<"leptonSkimOption: " <<leptonSkimOption<<" , "<<"treeTypeOption: "<<treeTypeOption<<endl;
    
    if (treeTypeOption == 1)
      events->CreateTree(ControlSampleEvents::kTreeType_OneLepton_Full);
    else if (treeTypeOption == 11)
      events->CreateTree(ControlSampleEvents::kTreeType_OneLepton_Reduced);
    else if (treeTypeOption == 2)
      events->CreateTree(ControlSampleEvents::kTreeType_OneLeptonAdd2MET_Full);
    else if (treeTypeOption == 3)
      events->CreateTree(ControlSampleEvents::kTreeType_Dilepton_Full);
    else if (treeTypeOption == 4)
      events->CreateTree(ControlSampleEvents::kTreeType_DileptonAdd2MET_Full);
    else if (treeTypeOption == 5)
      events->CreateTree(ControlSampleEvents::kTreeType_Photon_Full);
    else if (treeTypeOption == 6)
      events->CreateTree(ControlSampleEvents::kTreeType_ZeroLepton_Full);
    else if (treeTypeOption == 7)
      events->CreateTree(ControlSampleEvents::kTreeType_OneVetoLepton_Full);
    else if (treeTypeOption == 8)
      events->CreateTree(ControlSampleEvents::kTreeType_TightPlusVetoLepton_Full);
    else if (treeTypeOption == 9)
      events->CreateTree(ControlSampleEvents::kTreeType_OneTauLepton_Full);
    else {
      events->CreateTree(ControlSampleEvents::kTreeType_Default);
    }
    events->tree_->SetAutoFlush(0);

    //histogram containing total number of processed events (for normalization)
    int NEventProcessed = 0;
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
    TH1F *SumTopPtWeights = new TH1F("SumTopPtWeights", "SumTopPtWeights", 1, 1, 2);
 
    //*************************************************************************
    //Look over Input File Events
    //*************************************************************************
    if (fChain == 0) return;
    cout << "Total Events: " << fChain->GetEntries() << "\n";
    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {
      //initialize all varabiles
      events->InitVariables();

      //begin event
      if(jentry % 1000 == 0) cout << "Processing entry " << jentry << endl;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


      printSyncDebug = false;
      if (printSyncDebug) {
	cout << "\n****************************************************************\n";
	cout << "Debug Event : " << runNum << " " << lumiNum << " " << eventNum << "\n";
      }


      //fill normalization histogram
      NEventProcessed++;
      NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
      
      //event info
      events->option = option;
      events->genWeight = genWeight;
      events->run = runNum;
      events->lumi = lumiNum;
      events->event = eventNum;
      events->processID = 0;
 
      //get NPU
      for (int i=0; i < nBunchXing; ++i) {
	if (BunchXing[i] == 0) {
	  events->NPU_0 = nPUmean[i];
	}
	if (BunchXing[i] == -1) {
	  events->NPU_Minus1 = nPUmean[i];
	}
	if (BunchXing[i] == 1) {
	  events->NPU_Plus1 = nPUmean[i];
	}	  
      }
      events->NPV = nPV;
      events->genHT = getGenHT();
      events->NISRJets = getNISR( JetCorrector, JetCorrectorIOV );

        
      //************************************************************************************
      //Reconstructed vector boson momentum after pythia parton showering
      //************************************************************************************
      int genWBosonIndex = -1;
      int genZBosonIndex = -1;
      int genWLeptonIndex = -1;
      int genWNeutrinoIndex = -1;
      int genZLepton1Index = -1;
      int genZLepton2Index = -1;
      TLorentzVector genWVector;  genWVector.SetPtEtaPhiE(0,0,0,0);
      TLorentzVector genZVector;  genZVector.SetPtEtaPhiE(0,0,0,0);
      TLorentzVector genLepton;   genLepton.SetPtEtaPhiE(0,0,0,0);
      TLorentzVector genNeutrino; genNeutrino.SetPtEtaPhiE(0,0,0,0);
      TLorentzVector genZLepton1; genZLepton1.SetPtEtaPhiE(0,0,0,0);
      TLorentzVector genZLepton2; genZLepton2.SetPtEtaPhiE(0,0,0,0);

      //First find W or Z boson
      for(int j = 0; j < nGenParticle; j++){
	if ( gParticleStatus[j] == 22 && abs(gParticleId[j]) == 24) genWBosonIndex = j;
	if ( gParticleStatus[j] == 22 && abs(gParticleId[j]) == 23) genZBosonIndex = j;
      }

      //Find ttbar pair
      events->ptTop = -1;
      events->ptAntitop = -1;
      for(int j = 0; j < nGenParticle; j++){
          //top
          if ( gParticleStatus[j] == 22 && gParticleId[j] == 6  && events->ptTop < 0 ) {
              events->ptTop = gParticlePt[j];
          }
          //antitop
          if ( gParticleStatus[j] == 22 && gParticleId[j] == -6 && events->ptAntitop < 0 ) {
              events->ptAntitop = gParticlePt[j];
          }
      }
      // get top pt weight
      if ( events->ptTop > 0 && events->ptAntitop > 0 ) {
          events->topPtWeight = helper.getTopPtWeight( events->ptTop, events->ptAntitop );
          // fill sum of top pt weights
          SumTopPtWeights->SetBinContent( 1, SumTopPtWeights->GetBinContent(1) + events->topPtWeight);
      }

      //Next find the status 23 lepton and neutrinos from W or Z decay
      //If W or Z boson was found in the first step, require that they are daughters of the W or Z boson, 
      //if no W or Z boson was found, then don't require that.
      for(int j = 0; j < nGenParticle; j++){
	if ( gParticleStatus[j] == 23 && (abs(gParticleId[j]) == 11 || abs(gParticleId[j]) == 13 || abs(gParticleId[j]) == 15 )
	     && ( (genWBosonIndex >= 0 && gParticleMotherIndex[j] == genWBosonIndex) || genWBosonIndex == -1)		  
	     ) { 
	  genWLeptonIndex = j;
	}
      	if ( gParticleStatus[j] == 23 && (abs(gParticleId[j]) == 12 || abs(gParticleId[j]) == 14 || abs(gParticleId[j]) == 16 )
	     && ( (genWBosonIndex >= 0 && gParticleMotherIndex[j] == genWBosonIndex) || genWBosonIndex == -1)		  
	     ) { 
	  genWNeutrinoIndex = j;
	}
	if ( gParticleStatus[j] == 23 && ( gParticleId[j] == 11 || gParticleId[j] == 13 || gParticleId[j] == 15 )
	     && ((genZBosonIndex >= 0 && gParticleMotherIndex[j] == genZBosonIndex) || genZBosonIndex == -1)
	     ) { 
	  genZLepton1Index = j;
	}
	if ( gParticleStatus[j] == 23 && ( gParticleId[j] == -11 || gParticleId[j] == -13 || gParticleId[j] == -15 )
	     && ((genZBosonIndex >= 0 && gParticleMotherIndex[j] == genZBosonIndex) || genZBosonIndex == -1)
	     ) { 
	  genZLepton2Index = j;
	}
      }
    
      //Next collect all final state leptons and neutrinos from W or Z decay
      for(int j = 0; j < nGenParticle; j++){

	if (!(gParticleStatus[j] == 1 || (gParticleStatus[j] == 2 && abs(gParticleId[j]) == 15))) continue;
	    
	//if W was found, and lepton/neutrino daughter is stable
	if ( genWBosonIndex >= 0 && (abs(gParticleId[j]) >= 11 && abs(gParticleId[j]) <= 16 )
	     && gParticleMotherIndex[j] == genWBosonIndex ) {
	  TLorentzVector tmpVector; tmpVector.SetPtEtaPhiE(gParticlePt[j],gParticleEta[j],gParticlePhi[j],gParticleE[j]);
	  genWVector = genWVector + tmpVector;	  
	  if ( abs(gParticleId[j]) == 11 || abs(gParticleId[j]) == 13 || abs(gParticleId[j]) == 15) {
	    genLepton = genLepton + tmpVector;
	  } else {
	    genNeutrino = genNeutrino + tmpVector;
	  }
	}
	    
	//if status 22 lepton was found, find any daughters of it
	if ( genWLeptonIndex >= 0 && ( gParticleId[j] == gParticleId[genWLeptonIndex] || gParticleId[j] == 22)
	     && gParticleMotherIndex[j] == genWLeptonIndex ) {
	  TLorentzVector tmpVector; tmpVector.SetPtEtaPhiE(gParticlePt[j],gParticleEta[j],gParticlePhi[j],gParticleE[j]);
	  genWVector = genWVector + tmpVector;
	  genLepton = genLepton + tmpVector;	
	}
	//if status 22 neutrino was found, find any daughters of it
	if ( genWNeutrinoIndex >= 0 && ( gParticleId[j] == gParticleId[genWNeutrinoIndex] )
	     && gParticleMotherIndex[j] == genWNeutrinoIndex ) {
	  TLorentzVector tmpVector; tmpVector.SetPtEtaPhiE(gParticlePt[j],gParticleEta[j],gParticlePhi[j],gParticleE[j]);
	  genWVector = genWVector + tmpVector;	    
	  genNeutrino = genNeutrino + tmpVector;
	}
	 
	//if Z was found, and lepton daughter is stable
	if ( genZBosonIndex >= 0 && (abs(gParticleId[j]) >= 11 && abs(gParticleId[j]) <= 16 )
	     && gParticleMotherIndex[j] == genZBosonIndex ) {
	  TLorentzVector tmpVector; tmpVector.SetPtEtaPhiE(gParticlePt[j],gParticleEta[j],gParticlePhi[j],gParticleE[j]);
	  genZVector = genZVector + tmpVector;	
	  if ( gParticleId[j] > 0) {
	    genZLepton1 = genZLepton1 +  tmpVector;	
	  } else {
	    genZLepton2 = genZLepton2 +  tmpVector;
	  }
	}
	//if status 22 leptons were found, find any daughters of it
	if ( (gParticleId[j] == gParticleId[genZLepton1Index] || gParticleId[j] == 22)
	     && (genZLepton1Index >= 0 && gParticleMotherIndex[j] == genZLepton1Index)
	     ) {
	  TLorentzVector tmpVector; tmpVector.SetPtEtaPhiE(gParticlePt[j],gParticleEta[j],gParticlePhi[j],gParticleE[j]);
	  genZVector = genZVector + tmpVector;
	  genZLepton1 = genZLepton1 +  tmpVector;	
	}	   
	if ( (gParticleId[j] == gParticleId[genZLepton2Index] || gParticleId[j] == 22)
	     && ( genZLepton2Index >= 0 && gParticleMotherIndex[j] == genZLepton2Index)
	     ) {
	  TLorentzVector tmpVector; tmpVector.SetPtEtaPhiE(gParticlePt[j],gParticleEta[j],gParticlePhi[j],gParticleE[j]);
	  genZVector = genZVector + tmpVector;
	  genZLepton2 = genZLepton2 +  tmpVector;
	}
      }
  
      events->genWpt = genWVector.Pt();
      events->genWphi = genWVector.Phi();
      events->genZpt = genZVector.Pt();
      events->genZphi =  genZVector.Phi();


      //******************************************
      //Find Generated leptons
      //******************************************
      vector<int> genLeptonIndex;
      for(int j = 0; j < nGenParticle; j++){

	if ( gParticleStatus[j] == 22 && (abs(gParticleId[j]) == 11 || abs(gParticleId[j]) == 13 || abs(gParticleId[j]) == 15 )) {
	  genWLeptonIndex  = j;
	  //genLepton = makeTLorentzVector(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);      
	}
	if ( gParticleStatus[j] == 22 && (abs(gParticleId[j]) == 12 || abs(gParticleId[j]) == 14 || abs(gParticleId[j]) == 16 )) {
	  genWNeutrinoIndex = j;
	  //genNeutrino = makeTLorentzVector(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);      
	}
	if ( gParticleStatus[j] == 23 && (gParticleId[j] == 11 || gParticleId[j] == 13 || gParticleId[j] == 15 )) {
	  genZLepton1Index = j;
	  //genZLepton1 = makeTLorentzVector(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);      
	}
	if ( gParticleStatus[j] == 23 && (gParticleId[j] == -11 || gParticleId[j] == -13 || gParticleId[j] == -15 )) {
	  genZLepton2Index = j;
	  //genZLepton2 = makeTLorentzVector(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);      
	}
	
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
	

      //sort gen leptons by pt
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

      //*************************************************************************
      //Pileup Weights
      //*************************************************************************
      double NPU = 0;
      double pileupWeight = 1.0;
      if(!isData){
	//Get number of PU interactions
	for (int i = 0; i < nBunchXing; i++) {
	  if (BunchXing[i] == 0) {
	    NPU = nPUmean[i];
	  }
	}
        pileupWeight = helper.getPileupWeight(NPU);
        events->pileupWeight = pileupWeight;
      }


      //*************************************************************************
      //Find Reconstructed Leptons
      //*************************************************************************
      float probabilityToFail1LTrig = 1.0;

      vector<int> VetoLeptonIndex; 
      vector<int> VetoLeptonType;
      vector<int> VetoLeptonPt;
      vector<int> LooseLeptonIndex; 
      vector<int> LooseLeptonType;
      vector<int> LooseLeptonPt;
      vector<int> TightLeptonIndex; 
      vector<int> TightLeptonType;
      vector<int> TightLeptonPt;
      vector<TLorentzVector> GoodLeptons;//leptons used to compute hemispheres
      vector<int> GoodLeptonType;//leptons used to compute hemispheres
      vector<bool> GoodLeptonIsTight;//leptons used to compute hemispheres
      vector<bool> GoodLeptonIsMedium;//leptons used to compute hemispheres
      vector<bool> GoodLeptonIsLoose;//leptons used to compute hemispheres
      vector<bool> GoodLeptonIsVeto;//leptons used to compute hemispheres
      vector<double> GoodLeptonActivity;//leptons used to compute hemispheres

      //*******************************************************
      //Loop Over Muons
      //*******************************************************
      for(int i = 0; i < nMuons; i++){

	if(muonPt[i] < 5) continue;
	if(fabs(muonEta[i]) > 2.4) continue;

	//don't count duplicate muons 
	bool alreadySelected = false;
	for (uint j=0; j<GoodLeptons.size(); j++) {
	  if ( deltaR(GoodLeptons[j].Eta(), GoodLeptons[j].Phi(), muonEta[i],muonPhi[i]) < 0.1) alreadySelected = true;
	}
	if (alreadySelected) continue;

	if(isTightMuon(i) && muonPt[i] >= 20) {
	  TightLeptonType.push_back(13 * -1 * muonCharge[i]);
	  TightLeptonIndex.push_back(i);
	  TightLeptonPt.push_back(muonPt[i]);
	}
	else if(isVetoMuon(i)) {
	  VetoLeptonType.push_back(13 * -1 * muonCharge[i]);
	  VetoLeptonIndex.push_back(i);
	  VetoLeptonPt.push_back(muonPt[i]);
	}

	if (printSyncDebug) cout << "muon " << i << " " << muonPt[i] << " " << muonEta[i] << " " << muonPhi[i] << " : Tight = " << isTightMuon(i) << " Veto = " << isVetoMuon(i) << " \n";
                        	   
	TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
	thisMuon.SetPtEtaPhiM( muonPt[i], muonEta[i], muonPhi[i], 0.1057);

	//*******************************************************
	//For Single and Dilepton Options, use only tight leptons
	//For Veto Options, use only veto leptons
	//*******************************************************
	bool isGoodLepton = false;
	if (treeTypeOption == 1 || treeTypeOption == 2 || treeTypeOption == 11 || treeTypeOption == 12
	    || treeTypeOption == 3 || treeTypeOption == 4 || treeTypeOption == 13 || treeTypeOption == 14
	    ) {
	  if (isTightMuon(i)) isGoodLepton = true;
	}

	if (treeTypeOption == 6 || treeTypeOption == 7 || treeTypeOption == 9 
	    || treeTypeOption == 16 || treeTypeOption == 17 || treeTypeOption == 19 	    
	    ) {
	  if (isVetoMuon(i)) isGoodLepton = true;
	}
	
	if (treeTypeOption == 8 || treeTypeOption == 18) {
	  if (isVetoMuon(i) || isTightMuon(i)) isGoodLepton = true;
	}

	if (isGoodLepton) {
	  GoodLeptons.push_back(thisMuon);
	  GoodLeptonType.push_back(13 * -1 * muonCharge[i]);
	  GoodLeptonIsTight.push_back( isTightMuon(i) );
	  GoodLeptonIsMedium.push_back( isTightMuon(i) ); //tight and medium muon selections are currently the same
	  GoodLeptonIsLoose.push_back( isLooseMuon(i) );
	  GoodLeptonIsVeto.push_back( isVetoMuon(i) );
	  GoodLeptonActivity.push_back( muon_activityMiniIsoAnnulus[i] );
	}

	//*******************************************************
	//Compute Muon Trigger Efficiency Correction Factors
	//*******************************************************
	if (!isData && RazorAnalyzer::matchesGenMuon(muonEta[i], muonPhi[i])
	    && (treeTypeOption == 1 || treeTypeOption == 2 || treeTypeOption == 11 || treeTypeOption == 12
		|| treeTypeOption == 3 || treeTypeOption == 4 || treeTypeOption == 13 || treeTypeOption == 14
		)	       
	    && muonPt[i] > 20
	    ) {

	    //Trigger efficiency correction.  Note that this computes raw trigger efficiency, 
            //not MC-to-data scale factors.  MC-to-data trigger scale factors are not currently supported. 
            probabilityToFail1LTrig *= ( 1 - helper.getSingleMuTriggerEfficiency( muonPt[i], muonEta[i], true, true ) ); //update probability that no lepton fires a 1L trigger	  
	}

      } // loop over muons


      for(int i = 0; i < nElectrons; i++){

        //correct the electron pt in MC
        float eleCorrPt = elePt[i];
        if ( !isData ) {
            //this correction is disabled for now until we have implemented the POG recommended correction
            //eleCorrPt = helper.getCorrectedElectronPt( elePt[i], eleEta[i] );
        }

	if(eleCorrPt < 5) continue;
	if(fabs(eleEta[i]) > 2.5) continue;

	//don't count electrons that were already selected as muons
	bool alreadySelected = false;
	for (uint j=0; j<GoodLeptons.size(); j++) {
	  if ( deltaR(GoodLeptons[j].Eta(), GoodLeptons[j].Phi(), eleEta[i],elePhi[i]) < 0.1) alreadySelected = true;
	}
	if (alreadySelected) continue;

	if( isTightElectron(i) && eleCorrPt > 25 ) {
	  TightLeptonType.push_back(11 * -1 * eleCharge[i]);
	  TightLeptonIndex.push_back(i);
	  TightLeptonPt.push_back(eleCorrPt);
	}
	else if(isVetoElectron(i)) {
	  VetoLeptonType.push_back(11 * -1 * eleCharge[i]);
	  VetoLeptonIndex.push_back(i);
	  VetoLeptonPt.push_back(eleCorrPt);
	}
            
	if (printSyncDebug) cout << "ele " << i << " " << eleCorrPt << " " << eleEta[i] << " " << elePhi[i] << " : Tight = " << isTightElectron(i) << " Veto = " << isVetoElectron(i) << " \n";

	if(!isVetoElectron(i)) continue; 

	TLorentzVector thisElectron;
	if (isData) {
	  thisElectron.SetPtEtaPhiM( eleCorrPt, eleEta[i], elePhi[i], 0.000511);
	} else {
	  thisElectron.SetPtEtaPhiM( eleCorrPt, eleEta[i], elePhi[i], 0.000511);
	}

	//*******************************************************
	//For Single Lepton Options, use only tight leptons
	//For Dilepton Options, use only loose leptons
	//*******************************************************
	bool isGoodLepton = false;
	if (treeTypeOption == 1 || treeTypeOption == 2 || treeTypeOption == 11 || treeTypeOption == 12
	    || treeTypeOption == 3 || treeTypeOption == 4 || treeTypeOption == 13 || treeTypeOption == 14
	    ) {
	  if (isTightElectron(i)) isGoodLepton = true;
	}
	if ( treeTypeOption == 6 || treeTypeOption == 7 || treeTypeOption == 9 
	     || treeTypeOption == 16 || treeTypeOption == 17 || treeTypeOption == 19 
	     ) {
	  if (isVetoElectron(i) ) isGoodLepton = true;
	}
	if (treeTypeOption == 8 || treeTypeOption == 18) {
	  if (isVetoElectron(i) || isTightElectron(i)) isGoodLepton = true;
	}

	if (isGoodLepton) {
	  GoodLeptons.push_back(thisElectron);        
	  GoodLeptonType.push_back(11 * -1 * eleCharge[i]);
	  GoodLeptonIsTight.push_back( isTightElectron(i) );
	  GoodLeptonIsMedium.push_back( isMediumElectron(i) );
	  GoodLeptonIsLoose.push_back( isLooseElectron(i) );
	  GoodLeptonIsVeto.push_back( isVetoElectron(i) );
	  GoodLeptonActivity.push_back( ele_activityMiniIsoAnnulus[i] );
	}

	//*******************************************************
	//Compute Electron Trigger Efficiency Correction Factors
	//*******************************************************
	if (!isData && RazorAnalyzer::matchesGenElectron(eleEta[i],elePhi[i])
	    && (treeTypeOption == 1 || treeTypeOption == 2 || treeTypeOption == 11 || treeTypeOption == 12
		|| treeTypeOption == 3 || treeTypeOption == 4 || treeTypeOption == 13 || treeTypeOption == 14
		)
	    && eleCorrPt > 25	    
	    ) {
	  
	  //Trigger efficiency correction.  Note that this computes raw trigger efficiency, 
            //not MC-to-data scale factors.  MC-to-data trigger scale factors are not currently supported. 
	  probabilityToFail1LTrig *= ( 1 - helper.getSingleEleTriggerEfficiency( eleCorrPt, eleEta[i], true, true ) ); //update probability that no lepton fires a 1L trigger
	}

      } //loop over electrons

      for(int i = 0; i < nTaus; i++){
	if (tauPt[i] < 20) continue;
	if (fabs(tauEta[i]) > 2.4) continue;

	//don't count taus that were already selected as muons or electrons
	bool alreadySelected = false;
	for (uint j=0; j<GoodLeptons.size(); j++) {
	  if ( deltaR(GoodLeptons[j].Eta(), GoodLeptons[j].Phi(), tauEta[i],tauPhi[i]) < 0.1 ) alreadySelected = true;
	}
	if (alreadySelected) continue;

	if(isLooseTau(i)){
	  LooseLeptonType.push_back(15);
	  LooseLeptonIndex.push_back(i);
	  LooseLeptonPt.push_back(tauPt[i]);
	}

	TLorentzVector thisTau;
	thisTau.SetPtEtaPhiM( tauPt[i], tauEta[i], tauPhi[i], 1.777);

	//*******************************************************
	//For Single Lepton Options, use only tight leptons
	//For Dilepton Options, use only loose leptons
	//*******************************************************
	bool isGoodLepton = false;
	if (treeTypeOption == 9 || treeTypeOption == 19 ) {
	  if (isLooseTau(i)) isGoodLepton = true;
	}

	if (isGoodLepton) {
	  GoodLeptons.push_back(thisTau);        
	  GoodLeptonType.push_back(15);
	  GoodLeptonIsTight.push_back( isTightTau(i) );
	  GoodLeptonIsMedium.push_back( isMediumTau(i) );
	  GoodLeptonIsLoose.push_back( isLooseTau(i) );
	  GoodLeptonIsVeto.push_back( isLooseTau(i) );
	  GoodLeptonActivity.push_back( 9999 );
	}	  
      }

      //************************************************************************
      //Sort the collections by pT
      //************************************************************************
      for(uint i = 0; i < GoodLeptons.size() ; i++) {
	for (uint j=0; j < GoodLeptons.size()-1; j++) {
	  if (GoodLeptons[j+1].Pt() > GoodLeptons[j].Pt()) { 

	    // swap elements
	    TLorentzVector tmpV = GoodLeptons[j]; 
	    int tmpType = GoodLeptonType[j];
	    bool tmpIsTight = GoodLeptonIsTight[j];
	    bool tmpIsMedium = GoodLeptonIsMedium[j];
	    bool tmpIsLoose = GoodLeptonIsLoose[j];
	    bool tmpIsVeto = GoodLeptonIsVeto[j];
	    double tmpActivity = GoodLeptonActivity[j];

	    GoodLeptons[j] = GoodLeptons[j+1];
	    GoodLeptonType[j] = GoodLeptonType[j+1];
	    GoodLeptonIsTight[j] = GoodLeptonIsTight[j+1];
	    GoodLeptonIsMedium[j] = GoodLeptonIsMedium[j+1];
	    GoodLeptonIsLoose[j] = GoodLeptonIsLoose[j+1];
	    GoodLeptonIsVeto[j] = GoodLeptonIsVeto[j+1];
	    GoodLeptonActivity[j] = GoodLeptonActivity[j+1];

	    GoodLeptons[j+1] = tmpV;
	    GoodLeptonType[j+1] = tmpType;
	    GoodLeptonIsTight[j+1] = tmpIsTight;
	    GoodLeptonIsMedium[j+1] = tmpIsMedium;
	    GoodLeptonIsLoose[j+1] = tmpIsLoose;
	    GoodLeptonIsVeto[j+1] = tmpIsVeto;	  
	    GoodLeptonActivity[j+1] = tmpActivity;
	  }
	}
      }

      
      //************************************************************************
      //Fill Lepton Information using Good Leptons collection
      //************************************************************************
      events->lep1.SetPtEtaPhiM(0,0,0,0);
      events->lep2.SetPtEtaPhiM(0,0,0,0);
      events->lep1Type = 0;
      events->lep2Type = 0;
      events->lep1MatchedGenLepIndex = -1;
      events->lep2MatchedGenLepIndex = -1;
      events->lep1PassVeto = false;
      events->lep1PassLoose = false;
      events->lep1PassMedium = false;
      events->lep1PassTight = false;
      events->lep2PassVeto = false;
      events->lep2PassLoose = false;
      events->lep2PassTight = false;

      bool lep1Found = false;
      bool lep2Found = false;
      
      //************************************************************************
      //For single lepton and dilepton Control Regions
      //************************************************************************
      if ( treeTypeOption == 1 || treeTypeOption == 2 || treeTypeOption == 11 || treeTypeOption == 12 
	   || treeTypeOption == 3 || treeTypeOption == 4 || treeTypeOption == 13 || treeTypeOption == 14
	   || treeTypeOption == 7 || treeTypeOption == 17
	   || treeTypeOption == 9 || treeTypeOption == 19
	   ) {
	for (uint i=0; i<GoodLeptons.size(); i++) {

	  //For tau control region, consider only tau leptons
	  if (treeTypeOption == 9 || treeTypeOption == 19) {
	    if ( abs(GoodLeptonType[i]) != 15) continue;
	  }

	  if (!lep1Found) {
	    events->lep1Type = GoodLeptonType[i];
	    events->lep1.SetPtEtaPhiM(GoodLeptons[i].Pt(), GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),GoodLeptons[i].M());	  
	    events->lep1MatchedGenLepIndex = -1;
	    if ( abs(GoodLeptonType[i]) == abs(events->genlep1Type) &&
		 deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		 ) {
	      events->lep1MatchedGenLepIndex = 1;
	    }
	    if (abs(GoodLeptonType[i]) == abs(events->genlep2Type)
		&& deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		) {
	      events->lep1MatchedGenLepIndex = 2;
	    }
	    events->lep1PassTight = GoodLeptonIsTight[i];
	    events->lep1PassMedium = GoodLeptonIsMedium[i];
	    events->lep1PassLoose = GoodLeptonIsLoose[i];
	    events->lep1PassVeto = GoodLeptonIsVeto[i];
	    events->lep1Activity = GoodLeptonActivity[i];
	    lep1Found = true;
	  } else {
	    if (!lep2Found) {
	      events->lep2Type = GoodLeptonType[i];
	      events->lep2.SetPtEtaPhiM(GoodLeptons[i].Pt(), GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),GoodLeptons[i].M());	  	  	   
	      events->lep2MatchedGenLepIndex = -1;
	      if (abs(GoodLeptonType[i]) == abs(events->genlep1Type)
		  && deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 1;
	      }
	      if (abs(GoodLeptonType[i]) == abs(events->genlep2Type)
		  && deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 2;
	      }
	      events->lep2PassTight = GoodLeptonIsTight[i];
	      events->lep2PassLoose = GoodLeptonIsLoose[i];
	      events->lep2PassVeto = GoodLeptonIsVeto[i];
	      events->lep2Activity = GoodLeptonActivity[i];
	      lep2Found = true;
	    }
	  }
	} //loop over good leptons	


	//******************************************************
	//compute lepton efficiency scale factors 	
	//******************************************************
	double genElePtCut = 0;
	double genMuonPtCut = 0;
	bool useTightSelection = false;
	if ( treeTypeOption == 1 || treeTypeOption == 2 || treeTypeOption == 11 || treeTypeOption == 12 
	     || treeTypeOption == 3 || treeTypeOption == 4 || treeTypeOption == 13 || treeTypeOption == 14
	     ) {
	  genElePtCut = 25;
	  genMuonPtCut = 20;
	  useTightSelection = true;
	}
	if ( treeTypeOption == 7 || treeTypeOption == 17 ) {
	  genElePtCut = 5;
	  genMuonPtCut = 5;
	    useTightSelection = false;
	}
	  
	for (uint i=0;i<genLeptonIndex.size();i++) {

          // For each gen electron or muon passing the pt requirement, determine if it matches one of the selection electrons or muons
          // Then apply the appropriate scale factor
	  if (abs(gParticleId[genLeptonIndex[i]]) == 11 && gParticlePt[genLeptonIndex[i]] > genElePtCut && abs(gParticleEta[genLeptonIndex[i]]) < 2.5) {
	    bool isSelected = false;
	    for( uint j = 0; j < GoodLeptons.size(); j++ ) {
	      if (!(abs(GoodLeptonType[j]) == 11 && GoodLeptons[j].Pt() > genElePtCut)) continue;
	      if (useTightSelection) {
		if (!GoodLeptonIsTight[j]) continue;
	      } else {
		if (!GoodLeptonIsVeto[j]) continue;
	      }
	      double dR = deltaR( GoodLeptons[j].Eta() , GoodLeptons[j].Phi(), 
				  gParticleEta[genLeptonIndex[i]], gParticlePhi[genLeptonIndex[i]] );
	      if( dR < 0.1 ) {
		isSelected= true;
		break;
	      }
	    }
	    if (useTightSelection) {
	      events->eleEffWeight *= helper.getTightElectronScaleFactor( gParticlePt[genLeptonIndex[i]], 
									  gParticleEta[genLeptonIndex[i]], 
									  isSelected ); 
	    } else {
	      events->eleEffWeight *= helper.getVetoElectronScaleFactor( gParticlePt[genLeptonIndex[i]], 
									 gParticleEta[genLeptonIndex[i]], 
									  isSelected ); 
	    }
	  }

	  if (abs(gParticleId[genLeptonIndex[i]]) == 13 && gParticlePt[genLeptonIndex[i]] > genMuonPtCut && abs(gParticleEta[genLeptonIndex[i]]) < 2.4) {

	    bool isSelected = false;
	    for( uint j = 0; j < GoodLeptons.size(); j++ ) {
	      if (!(abs(GoodLeptonType[j]) == 13 && GoodLeptons[j].Pt() > genMuonPtCut)) continue;
	      if (useTightSelection) {
	      	if (!GoodLeptonIsTight[j]) continue;
	      } else {
	      	if (!GoodLeptonIsVeto[j]) continue;
	      }
	      double dR = deltaR( GoodLeptons[j].Eta() , GoodLeptons[j].Phi(),
	      			  gParticleEta[genLeptonIndex[i]], gParticlePhi[genLeptonIndex[i]] );
	      if( dR < 0.1 ) {
	      	isSelected = true;
	      	break;
	      }
	    }
	    if (useTightSelection) {
	      events->muonEffWeight *= helper.getTightMuonScaleFactor( gParticlePt[genLeptonIndex[i]], 
	    							       gParticleEta[genLeptonIndex[i]], 
	    							       isSelected ); 
	    } else {
	      events->muonEffWeight *= helper.getVetoMuonScaleFactor( gParticlePt[genLeptonIndex[i]], 
	    							      gParticleEta[genLeptonIndex[i]], 
	    							      isSelected ); 
	      
	    }
	  }
	}
	  	
      } // endif using 1-lepton or 2-lepton control region option
      

      //************************************************************************
      //For Loose Lepton + Veto lepton control region option
      //First select a Loose Lepton, assign it to lep1
      //Next select a Veto Lepton
       //************************************************************************
      if ( treeTypeOption == 8 || treeTypeOption == 18 ) {
	
	//look for a loose lepton first
	for (uint i=0; i<GoodLeptons.size(); i++) {
	  if (!GoodLeptonIsLoose[i]) continue;
	  if (!lep1Found) {
	    events->lep1Type = GoodLeptonType[i];
	    events->lep1.SetPtEtaPhiM(GoodLeptons[i].Pt(), GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),GoodLeptons[i].M());	  
	    events->lep1MatchedGenLepIndex = -1;
	    if ( abs(GoodLeptonType[i]) == abs(events->genlep1Type) &&
		 deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		 ) {
	      events->lep1MatchedGenLepIndex = 1;
	    }
	    if (abs(GoodLeptonType[i]) == abs(events->genlep2Type)
		&& deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		) {
	      events->lep1MatchedGenLepIndex = 2;
	    }
	    events->lep1PassTight = GoodLeptonIsTight[i];
	    events->lep1PassMedium = GoodLeptonIsMedium[i];
	    events->lep1PassLoose = GoodLeptonIsLoose[i];
	    events->lep1PassVeto = GoodLeptonIsVeto[i];
	    events->lep1Activity = GoodLeptonActivity[i];
	    lep1Found = true;
	    break;
	  }
	}

	//next look for a veto lepton
	if (lep1Found) {
	  for (uint i=0; i<GoodLeptons.size(); i++) {
	    if (!GoodLeptonIsVeto[i]) continue; 

	    //skip the object already saved as lep1
	    if ( deltaR( GoodLeptons[i].Eta(),GoodLeptons[i].Phi(), events->lep1.Eta(), events->lep1.Phi()) < 0.01) continue;

	    if (!lep2Found) {
	      events->lep2Type = GoodLeptonType[i];
	      events->lep2.SetPtEtaPhiM(GoodLeptons[i].Pt(), GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),GoodLeptons[i].M());	  	  	   
	      events->lep2MatchedGenLepIndex = -1;
	      if (abs(GoodLeptonType[i]) == abs(events->genlep1Type)
		  && deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep1.Eta(),events->genlep1.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 1;
	      }
	      if (abs(GoodLeptonType[i]) == abs(events->genlep2Type)
		  && deltaR(GoodLeptons[i].Eta(),GoodLeptons[i].Phi(),events->genlep2.Eta(),events->genlep2.Phi()) < 0.1
		  ) {
		events->lep2MatchedGenLepIndex = 2;
	      }
	      events->lep2PassTight = GoodLeptonIsTight[i];
	      events->lep2PassMedium = GoodLeptonIsMedium[i];
	      events->lep2PassLoose = GoodLeptonIsLoose[i];
	      events->lep2PassVeto = GoodLeptonIsVeto[i];
	      events->lep2Activity = GoodLeptonActivity[i];
	      lep2Found = true;
	    }
	  }
	}
		
      } //if using loose lepton + veto lepton control region option

      //************************************************************************
      //save this for the 1-lepton mini ntuples
      //************************************************************************
      if (lep1Found) {
	events->lep1Pt = events->lep1.Pt();
	events->lep1Eta = events->lep1.Eta();
      } else {
	events->lep1Pt = -99;
	events->lep1Eta = -99;
      }
      
      //************************************************************************
      //save for dilepton ntuples
      //************************************************************************
      if (lep1Found && lep2Found) {
	events->mll = (events->lep1 + events->lep2).M();
      }
      
      if (printSyncDebug) {
	cout << "\n\n";
	cout << "lep1: " << events->lep1Type << " | " << events->lep1.Pt() << " " << events->lep1.Eta() << " " << events->lep1.Phi() << " | Tight = " << events->lep1PassTight << " Loose = " << events->lep1PassLoose << " Veto = " << events->lep1PassVeto << "\n";
	cout << "lep2: " << events->lep2Type << " | " << events->lep2.Pt() << " " << events->lep2.Eta() << " " << events->lep2.Phi() << " | Tight = " << events->lep2PassTight << " Loose = " << events->lep2PassLoose << " Veto = " << events->lep2PassVeto << "\n";	 
      }
    

      //****************************************************//
      //             Select photons                         //
      //****************************************************//
      vector<TLorentzVector> GoodPhotons;
      vector<int> GoodPhotonIndex;
      int nPhotonsAbove40GeV = 0;
      
      if (treeTypeOption == 5 || treeTypeOption == 15) {
	for(int i = 0; i < nPhotons; i++){

	  float pho_pt_corr = pho_RegressionE[i]/cosh(phoEta[i]);//regression corrected pt
	  TVector3 vec;
	  vec.SetPtEtaPhi( pho_pt_corr, phoEta[i], phoPhi[i] );

	  if(pho_pt_corr < 30) continue;
	  if(fabs(phoEta[i]) > 2.5) continue;
	  
          //NOTE not running these for MC because trigger is broken in 80X samples
          if(isData) {
              if(!(pho_passHLTFilter[i][22] || pho_passHLTFilter[i][23]  || pho_passHLTFilter[i][24]  || pho_passHLTFilter[i][25]  || pho_passHLTFilter[i][26] 
                          || pho_passHLTFilter[i][27]  || pho_passHLTFilter[i][28] || pho_passHLTFilter[i][29] || pho_passHLTFilter[i][30]  
                          || pho_passHLTFilter[i][31]  || pho_passHLTFilter[i][32] || pho_passHLTFilter[i][33] || pho_passHLTFilter[i][34] 
                          || pho_passHLTFilter[i][35]  || pho_passHLTFilter[i][36] || pho_passHLTFilter[i][37] || pho_passHLTFilter[i][38]
                  )) continue;
          }
	  
	  double tmpEffAreaChargedHadrons = 0.0;
	  double tmpEffAreaNeutralHadrons = 0.0;
	  double tmpEffAreaPhotons = 0.0;
	  getPhotonEffArea90( pho_superClusterEta[i], tmpEffAreaChargedHadrons, tmpEffAreaNeutralHadrons, tmpEffAreaPhotons);
	  
	  //Use Loose Photon ID and charged Iso < 2.5 GeV
	  if (!( 		
		photonPassLooseID(i,true) //use this for deriving scale factor
		&& max(pho_sumChargedHadronPt[i] - fixedGridRhoFastjetAll*tmpEffAreaChargedHadrons, 0.) < 2.5
		&&
		
		photonPassesElectronVeto(i) //use this for performing fits
		&& pho_HoverE[i] < 0.05  
		&& ( 
		    (fabs(phoEta[i]) < 1.5  
		     && fmax(pho_pfIsoPhotonIso[i] - fixedGridRhoFastjetAll*tmpEffAreaPhotons, 0.0) < 3.630 + 0.0047*pho_pt_corr
		     && fmax(pho_pfIsoNeutralHadronIso[i] - fixedGridRhoFastjetAll*tmpEffAreaNeutralHadrons, 0.0) < 10.910 + 0.0148*pho_pt_corr + 0.000017*pho_pt_corr*pho_pt_corr
		     ) ||
		    (fabs(phoEta[i]) > 1.5  
		     && fmax(pho_pfIsoPhotonIso[i] - fixedGridRhoFastjetAll*tmpEffAreaPhotons, 0.) < 6.641 + 0.0034*pho_pt_corr
		     && fmax(pho_pfIsoNeutralHadronIso[i] - fixedGridRhoFastjetAll*tmpEffAreaNeutralHadrons, 0.) < 5.931 + 0.0163*pho_pt_corr + 0.000014*pho_pt_corr*pho_pt_corr
		     )
		     )
	       )	      
	      ) continue;
	  
	  if(pho_pt_corr > 40) nPhotonsAbove40GeV++;

	  TLorentzVector thisPhoton;
	  thisPhoton.SetVectM( vec, .0 );

	  GoodPhotons.push_back(thisPhoton);
	  GoodPhotonIndex.push_back(i);
	}
      }
      events->nSelectedPhotons = nPhotonsAbove40GeV;
      
      //Sort Photon Collection
      for(uint i = 0; i < GoodPhotons.size() ; i++) {
	for (uint j=0; j < GoodPhotons.size()-1; j++) {
	  if (GoodPhotons[j+1].Pt() > GoodPhotons[j].Pt()) { 
	    // swap elements
	    TLorentzVector tmpV = GoodPhotons[j]; 
	    int tmpIndex = GoodPhotonIndex[j];

	    GoodPhotons[j] = GoodPhotons[j+1];
	    GoodPhotonIndex[j] = GoodPhotonIndex[j+1];

	    GoodPhotons[j+1] = tmpV;
	    GoodPhotonIndex[j] = tmpIndex;
	  }
	}
      }

      //************************************************************************
      //Select Jets
      //************************************************************************
      bool bjet1Found = false;
      bool bjet2Found = false;
      vector<TLorentzVector> GoodJets;
      int numJetsAbove80GeV = 0;
      int numJetsAbove40GeV = 0;
      double MetX_Type1Corr = 0;
      double MetY_Type1Corr = 0;
      int nBJetsLoose40GeV = 0;
      int nBJetsMedium40GeV = 0;
      int nBJetsTight40GeV = 0;
      events->NGenBJets = 0;
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

      float mhx = 0., mhy = 0., mhx_nohf = 0., mhy_nohf = 0.;

      for(int i = 0; i < nJets; i++){

	//exclude Good leptons from the jet collection
	//NB: Currently the code excludes ALL identified "Good" leptons from the jet collection
	//even if the event is interpreted as a single lepton or two lepton event. This may cause
	//inconsistencies in the future.
	double dR = -1;
	for(auto& lep : GoodLeptons){
	  double thisDR = deltaR(jetEta[i],jetPhi[i],lep.Eta(),lep.Phi());
	  if(dR < 0 || thisDR < dR) dR = thisDR;
	}
	if(dR > 0 && dR < 0.4) continue; //jet matches a selected lepton
	
	
	if (printSyncDebug)  {
	  cout << "jet " << i << " : " << jetPt[i] << " " << jetEta[i] << " " << jetPhi[i] 
	       << " : rho = " << fixedGridRhoAll << " area = " << jetJetArea[i] << " "
	       << " | " 
	       << "correctedPt = " << jetPt[i]*JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
									 fixedGridRhoAll, jetJetArea[i], runNum,
									 JetCorrectorIOV, JetCorrector) << " "
	       << " | passID = " << jetPassIDTight[i] << " passPUJetID = " << bool((jetPileupIdFlag[i] & (1 << 2)) != 0) 
	       << " | csv = " << jetCSV[i] << " passCSVL = " << isCSVL(i) << " passCSVM = " << isCSVM(i) << " " << "\n";
	}

	//*******************************************************
	//apply jet iD
	//*******************************************************
	if (!jetPassIDTight[i]) continue;

	//*******************************************************
	//Correct Jet Energy Scale and Resolution
	//*******************************************************
	double tmpRho = fixedGridRhoFastjetAll;

	double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
					       tmpRho, jetJetArea[i],  runNum, 
					       JetCorrectorIOV,JetCorrector);   
	double JECLevel1 = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
						     tmpRho, jetJetArea[i], runNum,
						     JetCorrectorIOV,JetCorrector, 0);   
	double jetEnergySmearFactor = 1.0;
	if (!isData) {
	  std::vector<float> fJetEta, fJetPtNPU;
	  fJetEta.push_back(jetEta[i]);  
	  fJetPtNPU.push_back(jetPt[i]*JEC); 
	  fJetPtNPU.push_back(events->NPU_0); 
	  if (printSyncDebug) {
	    //cout << "Jet Resolution : " << jetPt[i]*JEC << " " << jetEta[i] << " " << jetPhi[i] << " : " 
	    // << JetResolutionCalculator->resolution(fJetEta,fJetPtNPU) << "\n";
	  }
	  //jetEnergySmearFactor = JetEnergySmearingFactor( jetPt[i]*JEC, jetEta[i], events->NPU_0, JetResolutionCalculator, random);
	}
	if (printSyncDebug) {
	  cout << "Jet Smearing Factor " << jetEnergySmearFactor << "\n";
	}

	TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);
	TLorentzVector L1CorrJet = makeTLorentzVector(jetPt[i]*JECLevel1, jetEta[i], jetPhi[i], jetE[i]*JECLevel1);
	TLorentzVector UnCorrJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);

	//**********************************************************************************************************
	//Add to Type1 Met Correction
	//Note: pT cut should be 10 not 20, but we're saving only 20 GeV jets in the razor ntuple for now
	//**********************************************************************************************************
	if (jetPt[i]*JEC*jetEnergySmearFactor > 20 && 
	    jetChargedEMEnergyFraction[i] + jetNeutralEMEnergyFraction[i] <= 0.9	    
	    ) {
	  MetX_Type1Corr += -1 * ( thisJet.Px() - L1CorrJet.Px()  );
	  MetY_Type1Corr += -1 * ( thisJet.Py() - L1CorrJet.Py()  );
	  if (printSyncDebug) cout << "Met Type1 Corr: " << thisJet.Px() - L1CorrJet.Px() << " " << thisJet.Py() - L1CorrJet.Py() << "\n";
	}


	//**********************************************************************************************
	//apply  Pileup Jet ID, don't do this yet for Run2 because we don't know the cut to make yet
	//**********************************************************************************************
	// int level = 2; //loose jet ID
	//if (!((jetPileupIdFlag[i] & (1 << level)) != 0)) continue;


	if (abs(jetPartonFlavor[i]) == 5) {
	  //count number of bjets in acceptance
	  if (thisJet.Pt() > BJET_PT_CUT && fabs(thisJet.Eta()) < 2.4) events->NGenBJets++;

	  if (!bjet1Found) {
	    bjet1Found = true;
	    events->bjet1.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), thisJet.M());
	    events->bjet1PassLoose = false;
	    events->bjet1PassMedium = false;
	    events->bjet1PassTight = false;
	    if(isCSVL(i)) events->bjet1PassLoose = true;	      
	    if(isCSVM(i)) events->bjet1PassMedium = true;
	    if(isCSVT(i)) events->bjet1PassTight = true;	      
	  } else if ( thisJet.Pt() > events->bjet1.Pt() ) {
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
	    if (!bjet2Found || thisJet.Pt() > events->bjet2.Pt() ) {
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
	} //if it's a bjet

	//---------------
	//b-tag corrector
	//---------------
	double thisJetPt = jetPt[i]*JEC*jetEnergySmearFactor;
	if ( !isData && abs( jetPartonFlavor[i] ) == 5 && fabs( jetEta[i] ) < 2.4 && thisJetPt > BJET_PT_CUT )
	  {
              float tmpSF = 1.0;
              helper.updateBTagScaleFactors( thisJetPt, jetEta[i], 5, isCSVM(i), 
                      events->btagW, events->btagW_up, events->btagW_down, 
                      tmpSF, tmpSF, tmpSF, tmpSF ); //scale factors for fastsim and mistag not used
	  }//btag correction

	/*std::cout << "jetPt: " << thisJetPt << " eta: " << jetEta[i] << " flavor: " << jetPartonFlavor[i]
	  << " events->btagW: " << events->btagW << std::endl;*/
	
	if (jetPt[i]*JEC > BJET_PT_CUT && isCSVL(i) ) nBJetsLoose40GeV++;
	if (jetPt[i]*JEC > BJET_PT_CUT && isCSVM(i) ) nBJetsMedium40GeV++;
	if (jetPt[i]*JEC > BJET_PT_CUT && isCSVT(i) ) nBJetsTight40GeV++;

	if(jetPt[i]*JEC*jetEnergySmearFactor < JET_PT_CUT) continue;

	// calculate MHT
	if(jetPt[i]*JEC*jetEnergySmearFactor > 20)
	  {
	    mhx -= jetPt[i]*JEC*jetEnergySmearFactor*cos(jetPhi[i]);
	    mhy -= jetPt[i]*JEC*jetEnergySmearFactor*sin(jetPhi[i]);
	    
	    if (fabs(jetEta[i]) < 3.0)
	      {
		mhx_nohf -= jetPt[i]*JEC*jetEnergySmearFactor*cos(jetPhi[i]);
		mhy_nohf -= jetPt[i]*JEC*jetEnergySmearFactor*sin(jetPhi[i]);
	      }
	  }
	
	if(fabs(jetEta[i]) > JET_ETA_CUT) continue;
	 
	numJetsAbove40GeV++;
	if(jetPt[i]*JEC*jetEnergySmearFactor > 80) numJetsAbove80GeV++;	  
	GoodJets.push_back(thisJet);	  
	  
      } //loop over jets
    
      /*
      std::cout << "btagW:  " << events->btagW << std::endl;
      std::cout << "btagW_up:  " << events->btagW_up << std::endl;
      std::cout << "btag_down:  " << events->btagW_down << std::endl;
      */

      //***************************************************************
      //W/top tagging
      //***************************************************************
      RazorHelper::AK8JetInfo jetInfo = helper.CalcAK8JetInfo(this, isData);
      events->nWTags = jetInfo.nWTags;
      events->nTopTags = jetInfo.nTopTags;
      events->wTagScaleFactor = jetInfo.wTagScaleFactor;
      events->topTagScaleFactor = jetInfo.topTagScaleFactor;
      
      events->MHT = sqrt(mhx*mhx + mhy*mhy);
      events->MHTnoHF = sqrt(mhx_nohf*mhx_nohf + mhy_nohf*mhy_nohf);

      //sort good jets
      sort(GoodJets.begin(), GoodJets.end(), greater_than_pt());
   
      //Fill two leading jets
      events->jet1.SetPtEtaPhiM(0,0,0,0);
      events->jet2.SetPtEtaPhiM(0,0,0,0);
      events->jet1PassCSVLoose = false;
      events->jet1PassCSVMedium = false;
      events->jet1PassCSVTight = false;
      events->jet2PassCSVLoose = false;
      events->jet2PassCSVMedium = false;
      events->jet2PassCSVTight = false;       
      for (int i=0;i<int(GoodJets.size());i++) {
	if (i==0) {
	  events->jet1.SetPtEtaPhiM(GoodJets[i].Pt(), GoodJets[i].Eta(),GoodJets[i].Phi(),GoodJets[i].M());
	  if(isCSVL(i)) events->jet1PassCSVLoose = true;	      
	  if(isCSVM(i)) events->jet1PassCSVMedium = true;
	  if(isCSVT(i)) events->jet1PassCSVTight = true;	
	}
	if (i==1) {
	  events->jet2.SetPtEtaPhiM(GoodJets[i].Pt(), GoodJets[i].Eta(),GoodJets[i].Phi(),GoodJets[i].M());
	  if(isCSVL(i)) events->jet2PassCSVLoose = true;	      
	  if(isCSVM(i)) events->jet2PassCSVMedium = true;
	  if(isCSVT(i)) events->jet2PassCSVTight = true; 	   	
	}
      }
    
      //***************************************************************
      //compute minDPhi variables
      //***************************************************************
    
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

    
      //***************************************************************************
      //Compute the razor variables using the selected jets and good leptons
      //***************************************************************************
      vector<TLorentzVector> GoodPFObjects;
      for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
      for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);

      double PFMetCustomType1X = metPt*cos(metPhi) + MetX_Type1Corr;
      double PFMetCustomType1Y = metPt*sin(metPhi) + MetY_Type1Corr;

      double PFMetnoHFX = metNoHFPt*cos(metNoHFPhi) + MetX_Type1Corr;
      double PFMetnoHFY = metNoHFPt*sin(metNoHFPhi) + MetY_Type1Corr;
    
      TLorentzVector PFMETUnCorr = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
      TLorentzVector PFMETCustomType1; 
      PFMETCustomType1.SetPxPyPzE(PFMetCustomType1X, PFMetCustomType1Y, 0, 
				  sqrt(PFMetCustomType1X*PFMetCustomType1X + PFMetCustomType1Y*PFMetCustomType1Y));
      TLorentzVector PFMETType1 = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);
      TLorentzVector PFMETType0Plus1 = makeTLorentzVectorPtEtaPhiM(metType0Plus1Pt, 0, metType0Plus1Phi, 0);
      TLorentzVector PFMETMuEGCleanCorr = makeTLorentzVectorPtEtaPhiM(metMuEGCleanCorrPt, 0, metMuEGCleanCorrPhi, 0);

      TLorentzVector PFMETnoHFType1;
      PFMETnoHFType1.SetPxPyPzE(PFMetnoHFX, PFMetnoHFY, 0, sqrt(PFMetnoHFX*PFMetnoHFX + PFMetnoHFY*PFMetnoHFY));
      
      TLorentzVector MyMET;
      if(isData) {
          MyMET = PFMETMuEGCleanCorr; //Using this for V3.13 ntuples run on 2017 Feb 03 reMiniAOD
      }
      else {
          MyMET = PFMETCustomType1;
      }
	
      if (printSyncDebug) {
	cout << "UnCorrectedMET: " << PFMETUnCorr.Pt() << " " << PFMETUnCorr.Phi() << "\n";
	cout << "Corrected PFMET: " << PFMETCustomType1.Pt() << " " << PFMETCustomType1.Phi() << " | X,Y Correction :  " << MetX_Type1Corr << " " << MetY_Type1Corr << "\n";
      }

      events->MR = 0;
      events->Rsq = 0;
	
      //only compute razor variables if we have 2 jets
      //Use Type 1 PFMet for MET
      if (GoodPFObjects.size() >= 2 
	  && GoodJets.size() < 20 //this is a protection again crazy events
	  ) {
	vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
	events->MR = computeMR(hemispheres[0], hemispheres[1]); 
	events->Rsq = computeRsq(hemispheres[0], hemispheres[1], MyMET);
	events->RsqnoHF = computeRsq(hemispheres[0], hemispheres[1], PFMETnoHFType1);
        events->dPhiRazor = fabs(hemispheres[0].DeltaPhi(hemispheres[1]));
      }
	

      if (printSyncDebug)  {
	cout << "MR = " << events->MR << " Rsq = " << events->Rsq << " | "
	     << " Mll = " << (events->lep1 + events->lep2).M() << " | " 
	     << " NJets80 = " << numJetsAbove80GeV << " NJets40 = " << numJetsAbove40GeV << " GoodPFObjects.size() = " << GoodPFObjects.size() << " "
	     << " MET = " << MyMET.Pt() << " MetPhi = " << MyMET.Phi() << " nBTagsMedium = " << nBJetsMedium40GeV << "\n";
      }

      events->MET = MyMET.Pt();
      events->METPhi = MyMET.Phi();
      events->METnoHF = PFMETnoHFType1.Pt();
      events->METnoHFPhi = PFMETnoHFType1.Phi();
      events->METRaw = PFMETUnCorr.Pt();
      events->METRawPhi = PFMETUnCorr.Phi();
      events->NJets40 = numJetsAbove40GeV;
      events->NJets80 = numJetsAbove80GeV;
      events->NBJetsLoose = nBJetsLoose40GeV;
      events->NBJetsMedium = nBJetsMedium40GeV;
      events->NBJetsTight = nBJetsTight40GeV;

      ///
      events-> metType1PtJetResUp          = metType1PtJetResUp;	      
      events-> metType1PtJetResDown        = metType1PtJetResDown;
      events-> metType1PtJetEnUp           = metType1PtJetEnUp;	 
      events-> metType1PtJetEnDown         = metType1PtJetEnDown;
      events-> metType1PtMuonEnUp          = metType1PtMuonEnUp;
      events-> metType1PtMuonEnDown        = metType1PtMuonEnDown;
      events-> metType1PtElectronEnUp      = metType1PtElectronEnUp;
      events-> metType1PtElectronEnDown    = metType1PtElectronEnDown;
      events-> metType1PtTauEnUp           = metType1PtTauEnUp;
      events-> metType1PtTauEnDown         = metType1PtTauEnDown;
      events-> metType1PtUnclusteredEnUp   = metType1PtUnclusteredEnUp;
      events-> metType1PtUnclusteredEnDown = metType1PtUnclusteredEnDown;
      events-> metType1PtPhotonEnUp        = metType1PtPhotonEnUp;
      events-> metType1PtPhotonEnDown      = metType1PtPhotonEnDown;

      events-> metType1PhiJetResUp          = metType1PhiJetResUp;	      
      events-> metType1PhiJetResDown        = metType1PhiJetResDown;
      events-> metType1PhiJetEnUp           = metType1PhiJetEnUp;	 
      events-> metType1PhiJetEnDown         = metType1PhiJetEnDown;
      events-> metType1PhiMuonEnUp          = metType1PhiMuonEnUp;
      events-> metType1PhiMuonEnDown        = metType1PhiMuonEnDown;
      events-> metType1PhiElectronEnUp      = metType1PhiElectronEnUp;
      events-> metType1PhiElectronEnDown    = metType1PhiElectronEnDown;
      events-> metType1PhiTauEnUp           = metType1PhiTauEnUp;
      events-> metType1PhiTauEnDown         = metType1PhiTauEnDown;
      events-> metType1PhiUnclusteredEnUp   = metType1PhiUnclusteredEnUp;
      events-> metType1PhiUnclusteredEnDown = metType1PhiUnclusteredEnDown;
      events-> metType1PhiPhotonEnUp        = metType1PhiPhotonEnUp;
      events-> metType1PhiPhotonEnDown      = metType1PhiPhotonEnDown;
      
      ///

      events->HT = 0;
      for(auto& jet : GoodJets) events->HT += jet.Pt();
	
      //compute M_T for lep1 and MET
      events->lep1MT = sqrt(events->lep1.M2() + 2*MyMET.Pt()*events->lep1.Pt()*(1 - cos(deltaPhi(MyMET.Phi(),events->lep1.Phi()))));
      events->lep1MTnoHF = sqrt(events->lep1.M2() + 2*PFMETnoHFType1.Pt()*events->lep1.Pt()*(1 - cos(deltaPhi(PFMETnoHFType1.Phi(),events->lep1.Phi()))));
	
      //save HLT Decisions
      for(int k=0; k<200; ++k) {
	events->HLTDecision[k] = HLTDecision[k];	
	events->HLTPrescale[k] = HLTPrescale[k];	
      }
    	
      //MET Filter
      events->Flag_HBHENoiseFilter = Flag_HBHENoiseFilter;
      events->Flag_HBHEIsoNoiseFilter = Flag_HBHEIsoNoiseFilter;
      events->Flag_badChargedCandidateFilter = Flag_badChargedCandidateFilter;
      events->Flag_badMuonFilter = Flag_badMuonFilter;
      events->Flag_badGlobalMuonFilter = Flag_badGlobalMuonFilter;
      events->Flag_duplicateMuonFilter = Flag_duplicateMuonFilter;
      events->Flag_CSCTightHaloFilter = Flag_CSCTightHaloFilter;
      events->Flag_hcalLaserEventFilter = Flag_hcalLaserEventFilter;
      events->Flag_EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilter;
      events->Flag_goodVertices = Flag_goodVertices;
      events->Flag_trackingFailureFilter = Flag_trackingFailureFilter;
      events->Flag_eeBadScFilter = Flag_eeBadScFilter;
      events->Flag_ecalLaserCorrFilter = Flag_ecalLaserCorrFilter;
      events->Flag_trkPOGFilters = Flag_trkPOGFilters;
      events->Flag_trkPOG_manystripclus53X = Flag_trkPOG_manystripclus53X;
      events->Flag_trkPOG_toomanystripclus53X = true;
      events->Flag_trkPOG_logErrorTooManyClusters = Flag_trkPOG_logErrorTooManyClusters;
      events->Flag_METFilters = Flag_METFilters;	


      //***********************************************************************//
      //    Compute razor vars adding photon to the MET
      //***********************************************************************//
	
      if(GoodPhotons.size()>0){

	events->pho1 = GoodPhotons[0];
	if(GoodPhotons.size()>1)
	  events->pho2 = GoodPhotons[1];
	  

	// match photons to gen particles to remove double counting between QCD and GJet samples
	double phoMinDR = 9999;
	for(int g = 0; g < nGenParticle; g++){
	  if (!(deltaR(gParticleEta[g] , gParticlePhi[g], GoodPhotons[0].Eta(),GoodPhotons[0].Phi()) < 0.1) ) continue;
	  if(!(gParticleStatus[g] == 1 || gParticleStatus[g] == 23)) continue;
	  if(!(gParticleId[g] == 22 || abs(gParticleId[g]) == 11)) continue;
	  
	  double tmpMinDR = deltaR(GoodPhotons[0].Eta(), GoodPhotons[0].Phi(), gParticleEta[g] , gParticlePhi[g]);
	  if ( tmpMinDR < phoMinDR) {
	    phoMinDR = tmpMinDR;
	    if (gParticleStatus[g] == 23 || abs(gParticleId[g]) == 11 ) {
	      events->pho1_motherID = 22;
		} else {
	      events->pho1_motherID = gParticleMotherId[g];
	    }

	    //cout << "Found gen photon: " << gParticlePt[g] << " " << gParticleEta[g] << " " << gParticlePhi[g] << "\n";

	    //find dR to closest parton
	    double tmpMinDRPhotonToParton = 9999;
	    for(int l = 0; l < nGenParticle; l++){
	      if (l == g) continue; //extra protection, don't compare against the gen photon
	      if (!( abs(gParticleId[l]) <= 5 || gParticleId[l] == 21 )) continue;
	      if (gParticleStatus[l] != 23) continue;

	      //cout << "parton " << l << " : " << gParticleId[l] << " " << gParticleStatus[l] << " : " << gParticlePt[l] << " " << gParticleEta[l] << " " << gParticlePhi[l] << "\n";

	      if ( deltaR(gParticleEta[g] , gParticlePhi[g], gParticleEta[l] , gParticlePhi[l]) < tmpMinDRPhotonToParton) {
		tmpMinDRPhotonToParton = deltaR(gParticleEta[g] , gParticlePhi[g], gParticleEta[l] , gParticlePhi[l]);
	      }
	    }
	    events->minDRGenPhotonToParton = tmpMinDRPhotonToParton;
	  }
	}

	if (events->minDRGenPhotonToParton < 0) {
	  for(int g = 0; g < nGenParticle; g++){
	    if ( deltaR(GoodPhotons[0].Eta(), GoodPhotons[0].Phi(), gParticleEta[g] , gParticlePhi[g]) < 0.5) {
	      cout << "genparticle " << g << " : " << gParticleId[g] << " : " << gParticleStatus[g] << " : " << gParticlePt[g] << " " << gParticleEta[g] << " " << gParticlePhi[g] << " : " <<  deltaR(GoodPhotons[0].Eta(), GoodPhotons[0].Phi(), gParticleEta[g] , gParticlePhi[g]) << "\n";
	    }
	  }
	}

	for (int q=0; q<50; q++) events->pho1HLTFilter[q] = pho_passHLTFilter[GoodPhotonIndex[0]][q];
	events->pho1_sigmaietaieta = phoFull5x5SigmaIetaIeta[GoodPhotonIndex[0]];
	double effAreaChargedHadrons = 0.0;
	double effAreaNeutralHadrons = 0.0;
	double effAreaPhotons = 0.0;
	getPhotonEffArea90( pho_superClusterEta[GoodPhotonIndex[0]] , effAreaChargedHadrons, effAreaNeutralHadrons, effAreaPhotons);
	events->pho1_chargediso = max(pho_sumChargedHadronPt[GoodPhotonIndex[0]] - fixedGridRhoFastjetAll*effAreaChargedHadrons, 0.);
	events->pho1_photoniso = max(pho_sumPhotonEt[GoodPhotonIndex[0]] - fixedGridRhoFastjetAll*effAreaPhotons, 0.);
	events->pho1_neutralhadroniso = max(pho_sumNeutralHadronEt[GoodPhotonIndex[0]] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.);
	events->pho1_hOverE = pho_HoverE[GoodPhotonIndex[0]];
	events->pho1_pfiso = max(pho_sumChargedHadronPt[GoodPhotonIndex[0]] - fixedGridRhoFastjetAll*effAreaChargedHadrons, 0.) +
	  max(pho_sumNeutralHadronEt[GoodPhotonIndex[0]] - fixedGridRhoFastjetAll*effAreaNeutralHadrons, 0.) +
	  max(pho_sumPhotonEt[GoodPhotonIndex[0]] - fixedGridRhoFastjetAll*effAreaPhotons, 0.);
	events->pho1PassTight = photonPassTightID(GoodPhotonIndex[0], true);

	//compute MET with leading photon added
	TLorentzVector m1 = GoodPhotons[0];
	TLorentzVector m2 = MyMET;
	TLorentzVector photonPlusMet_perp = makeTLorentzVectorPtEtaPhiM((m1 + m2).Pt(), 0., (m1 + m2).Phi(), 0.0);

	events->MET_NoPho = photonPlusMet_perp.Pt();
	events->METPhi_NoPho = photonPlusMet_perp.Phi();

	//remove leading photon from collection of selected jets
	vector<TLorentzVector> GoodJetsNoLeadPhoton = GoodJets;
	int subtractedIndex = SubtractParticleFromCollection(GoodPhotons[0], GoodJetsNoLeadPhoton);
	if(subtractedIndex >= 0){
	  if(GoodJetsNoLeadPhoton[subtractedIndex].Pt() < JET_PT_CUT){ //erase this jet
	    GoodJetsNoLeadPhoton.erase(GoodJetsNoLeadPhoton.begin()+subtractedIndex);
	  }
	}	
	    
	//count the number of jets above 80 GeV now
	int numJets80_noPho = 0.;
	for(auto& jet : GoodJetsNoLeadPhoton){
	  if(jet.Pt() > 80) numJets80_noPho++;
	}
	events->NJets80_NoPho = numJets80_noPho;
	    
	//count jets and compute HT
	events->NJets_NoPho = GoodJetsNoLeadPhoton.size();
	float ht_noPho = 0.;
	for(auto& pf : GoodJetsNoLeadPhoton) ht_noPho += pf.Pt();
	events->HT_NoPho = ht_noPho;
	    
	if(GoodJetsNoLeadPhoton.size() >= 2 && GoodJetsNoLeadPhoton.size() <20){
	  //remake the hemispheres using the new jet collection
	  vector<TLorentzVector> hemispheresNoLeadPhoton = getHemispheres(GoodJetsNoLeadPhoton);
	  TLorentzVector PFMET_NOPHO = makeTLorentzVectorPtEtaPhiM(events->MET_NoPho, 0, events->METPhi_NoPho, 0);
	  events->MR_NoPho = computeMR(hemispheresNoLeadPhoton[0], hemispheresNoLeadPhoton[1]); 
	  events->Rsq_NoPho = computeRsq(hemispheresNoLeadPhoton[0], hemispheresNoLeadPhoton[1], PFMET_NOPHO);
	  events->dPhiRazor_NoPho = fabs(hemispheresNoLeadPhoton[0].DeltaPhi(hemispheresNoLeadPhoton[1]));
	}

      }

      // GenJET MR and HT
      vector<TLorentzVector> GenJetObjects;
      for(int j = 0; j < nGenJets; j++){
	if (genJetPt[j] > JET_PT_CUT && fabs(genJetEta[j]) < 3) {
	  TLorentzVector thisGenJet = makeTLorentzVector(genJetPt[j], genJetEta[j], genJetPhi[j], genJetE[j]);
	  GenJetObjects.push_back(thisGenJet);
	  events->genJetHT += genJetPt[j];
	}
      }
      if (GenJetObjects.size() >= 2 ) {
	vector<TLorentzVector> tmpHemispheres = getHemispheres(GenJetObjects);
	events->genJetMR = computeMR(tmpHemispheres[0], tmpHemispheres[1]); 
      }
    
      //***********************************************************************//
      //    Compute razor vars adding good leptons to the MET
      //***********************************************************************//
      //remove selected muons from collection of selected jets and add them to the MET
      vector<TLorentzVector> GoodJetsNoLeptons = GoodJets;
      TLorentzVector TotalLepVec;
      int NLeptonsAdded = 0;
      for(auto& lep : GoodLeptons){
	NLeptonsAdded++;
	TotalLepVec = TotalLepVec + lep; //add this muon's momentum to the sum
	int subtractedIndex = SubtractParticleFromCollection(lep, GoodJetsNoLeptons);
	if(subtractedIndex >= 0){
	  if(GoodJetsNoLeptons[subtractedIndex].Pt() < JET_PT_CUT){ //erase this jet
	    GoodJetsNoLeptons.erase(GoodJetsNoLeptons.begin()+subtractedIndex);
	  }
	}
	
	//only count 1 lepton for 1-lepton option. 
	if (treeTypeOption == 2) {
	  if (NLeptonsAdded >= 1) break;
	}

	//only count 2 lepton for 2-lepton option. 
	if (treeTypeOption == 4) {
	  if (NLeptonsAdded >= 2) break;
	}

     }

      //make the MET vector with the muons (or gen muons) added
      TLorentzVector LepPlusMet_perp = makeTLorentzVector((TotalLepVec + MyMET).Pt(), 0., (TotalLepVec + MyMET).Phi(), 0.);

      //count jets and compute HT
      int njets80NoLep = 0;
      float ht_NoLep   = 0.;
      for(auto& jet : GoodJetsNoLeptons){
	ht_NoLep += jet.Pt();
	if(jet.Pt() > 80) njets80NoLep++;
      }


      //**********************************************************
      //1-lepton Add to MET Variables
      //**********************************************************
      if (treeTypeOption == 2) 
      {
	events->MET_NoW = LepPlusMet_perp.Pt();
	events->METPhi_NoW = LepPlusMet_perp.Phi();
	events->NJets_NoW = GoodJetsNoLeptons.size();
	events->HT_NoW = ht_NoLep;
	events->NJets80_NoW = njets80NoLep;
    events->METPhi = MyMET.Phi();

	//compute razor variables 
	if(events->NJets_NoW > 1 && GoodJets.size()<20){
	  vector<TLorentzVector> hemispheresNoW = getHemispheres(GoodJetsNoLeptons);
	  events->Rsq_NoW = computeRsq(hemispheresNoW[0], hemispheresNoW[1], LepPlusMet_perp);
	  events->MR_NoW = computeMR(hemispheresNoW[0], hemispheresNoW[1]); 
	  events->dPhiRazor_NoW = fabs(hemispheresNoW[0].DeltaPhi(hemispheresNoW[1])); 
	}

	//more W kinematics
	if (GoodLeptons.size() > 0) {
	  events->recoWpt = (GoodLeptons[0] + MyMET).Pt();
	  events->recoWphi = (GoodLeptons[0] + MyMET).Phi();
	}
      }


      //**********************************************************
      //2-lepton Add to MET Variables
      //**********************************************************
      if (treeTypeOption == 4) {
	events->MET_NoZ = LepPlusMet_perp.Pt();
	events->METPhi_NoZ = LepPlusMet_perp.Phi();
	events->NJets_NoZ = GoodJetsNoLeptons.size();
	events->HT_NoZ = ht_NoLep;
	events->NJets80_NoZ = njets80NoLep;
    events->METPhi = MyMET.Phi();

    // Compute the recoil for Z
    TVector2 vZPt((TotalLepVec.Pt())*cos(TotalLepVec.Phi()),(TotalLepVec.Pt())*sin(TotalLepVec.Phi())); // Z bosont pT
    TVector2 vMet((MyMET.Pt())*cos(MyMET.Phi()), (MyMET.Pt())*sin(MyMET.Phi())); // MET vector
    TVector2 vU = -1.0*(vZPt+vMet); // recoil vector
    events->u1 = ((TotalLepVec.Px())*(vU.Px()) + (TotalLepVec.Py())*(vU.Py()))/(TotalLepVec.Pt());  // u1 = (pT . u)/|pT|
    events->u2 = ((TotalLepVec.Px())*(vU.Py()) - (TotalLepVec.Py())*(vU.Px()))/(TotalLepVec.Pt());  // u2 = (pT x u)/|pT|

	//Z kinematics
	if(NLeptonsAdded >= 2){
	  events->recoZpt = TotalLepVec.Pt();
	  events->recoZmass = TotalLepVec.M();
	} else {
	  events->recoZpt = -1;
	  events->recoZmass = -1;
	}

	//compute razor variables 
	if(events->NJets_NoZ > 1 && GoodJets.size()<20)
	  {
	    vector<TLorentzVector> hemispheresNoZ = getHemispheres(GoodJetsNoLeptons);
	    events->Rsq_NoZ = computeRsq(hemispheresNoZ[0], hemispheresNoZ[1], LepPlusMet_perp);
	    events->MR_NoZ = computeMR(hemispheresNoZ[0], hemispheresNoZ[1]); 
	    events->dPhiRazor_NoZ = fabs(hemispheresNoZ[0].DeltaPhi(hemispheresNoZ[1])); 
	  }
      }
    

      //****************************************************************************
      //Compute All Event Weights
      //****************************************************************************

      events->weight = events->genWeight
	* pileupWeight
	* events->muonEffWeight * events->eleEffWeight 
	* events->btagW;

      // lepton trigger weight
      if ( (treeTypeOption == 1 || treeTypeOption == 2 || treeTypeOption == 11 || treeTypeOption == 12
                  || treeTypeOption == 3 || treeTypeOption == 4 || treeTypeOption == 13 || treeTypeOption == 14 ) ) {
          events->trigWeight1L = 1 - probabilityToFail1LTrig; //efficiency = 1 - P( all leptons fail )

	  //if there are no leptons in the event, then probabilityToFail1LTrig remains at 1.0 and trigWeight1L is 0
	  //we will give such events a weight of 1.0, to preserve information of other weights
          if ( events->trigWeight1L == 0 ) events->trigWeight1L = 1.0; 

          events->weight *= events->trigWeight1L;
      }
      

      //****************************************************************************
      //Event Skimming
      //****************************************************************************
      bool passSkim = true;

      // Dilepton skim
      if (leptonSkimOption == 2) {
	if (!(
	      (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	      && (abs(events->lep2Type) == 11  || abs(events->lep2Type) == 13 )
	      && events->lep1.Pt() > 0 && events->lep2.Pt() > 0
	      )
	    ) passSkim = false;
      }

      //single lepton skim
      if (leptonSkimOption == 1) {
	if (!( 
	      (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	      && events->lep1.Pt() > 15
	       )
	    ) passSkim = false;
      }
      //single lepton skim
      if (leptonSkimOption == 3) {
	if (!( events->lep1.Pt() > 0
	       )
	    ) passSkim = false;
      }

      //Photon Skim
      if(leptonSkimOption == 5)
	{
	  if (!(
		GoodPhotons.size() > 0
		&& GoodPhotons[0].Pt() > 25	
		)
	      ) passSkim = false;	
	}
      
      //razor skim
      if (razorSkimOption == 1) {
	if ( !( 
	       (events->MR > 300 && events->Rsq > 0.1) || GoodJets.size() >= 20
		)
	     ) passSkim = false;
      }

      //1-lepton add to MET razor skim
      if(razorSkimOption == 2)
	{
	  if (!(
		events->MR_NoW > 300 && events->Rsq_NoW > 0.15
		)
	      ) passSkim = false;		
	}
      

      //Dilepton add to MET razor skim
      if(razorSkimOption == 3) 
	{
	  if (!(
		events->MR_NoZ > 300 && events->Rsq_NoZ > 0.15
		)
	      ) passSkim = false;
	}
      
      //Photon Skim
      if(razorSkimOption == 5)
	{
	  if (!(
		events->MR_NoPho > 300
		&& events->Rsq_NoPho > 0.15
		)
	      ) passSkim = false;	
	}
      
      //RazorDM Skim
      if(razorSkimOption == 6)
	{
	  if (!(
		events->MR > 150 && events->Rsq > 0.4
		)
	      ) passSkim = false;	
	}

      //1-lepton add to MET razorDM skim
      if(razorSkimOption == 7)
	{
	  if (!(
		events->MR_NoW > 150 && events->Rsq_NoW > 0.4
		)
	      ) passSkim = false;		
	}
      
      //Dilepton add to MET razorDM skim
      if(razorSkimOption == 8) 
	{
	  if (!(
		events->MR_NoZ > 150 && events->Rsq_NoZ > 0.4
		)
	      ) passSkim = false;
	}
      
      //Photon razorDM Skim
      if(razorSkimOption == 9)
	{
	  if (!(
		events->MR_NoPho > 150 && events->Rsq_NoPho > 0.4
		)
	      ) passSkim = false;	
	}
      

      //fill event 
      if (passSkim) {
	if (printSyncDebug) cout<<"Filling the tree... " <<endl;
	events->tree_->Fill();
      }


    }//end of event loop


    cout << "Filled Total of " << NEventProcessed << " Events\n";
    cout << "Sum of event weights = " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}

