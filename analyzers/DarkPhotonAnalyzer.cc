//LOCAL INCLUDES
#include "DarkPhotonAnalyzer.h"
#include "RazorHelper.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
//C++ INCLUDES
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h> 
//ROOT INCLUDES
#include <TH1F.h>
#include <TMath.h>
#include <TH2D.h>
#include "TRandom3.h"

using namespace std;

enum HggRazorLeptonsBox {
  Zmm = 0,
  Zee = 1,
  OneMu = 2,
  OneEle = 3,
  None = 10
};

const int NUM_PDF_WEIGHTS = 60;

void DarkPhotonAnalyzer::Analyze(bool isData, int Option, string outputFilename, string label)
{
  cout << "Initializing..." << endl;

  //--------------------------------
  //Initialize helper
  //--------------------------------
  string analysisTag = "Razor2016_80X";
  if ( label != "") analysisTag = label;
  bool isFastsimSMS = false;
  RazorHelper *helper = 0;
  if (analysisTag == "Razor2015_76X") helper = new RazorHelper("Razor2015_76X", isData, isFastsimSMS);
  else if (analysisTag == "Razor2016_80X") helper = new RazorHelper("Razor2016_80X", isData, isFastsimSMS);
  else helper = new RazorHelper(analysisTag, isData, isFastsimSMS);


  //--------------------------------
  //Initialize Output
  //--------------------------------
  string outfilename = outputFilename;
  if (outfilename == "") outfilename = "/afs/cern.ch/user/j/jbamber/scratch1/CMSSW_7_4_15/src/ROOT_outputs/DarkPhoton.root";
  TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");

  cout << "Run With Option = " << Option << "\n";

  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
  TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
  TH1F *SumPdfWeights = new TH1F("SumPdfWeights", "SumPdfWeights", NUM_PDF_WEIGHTS, -0.5, NUM_PDF_WEIGHTS-0.5);
    
  //--------------
  //tree variables
  //--------------
  unsigned int run, lumi, event;
  float weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  float triggerEffWeight;
  float triggerEffSFWeight;
  float photonEffSF;

  int NPU;
  float MET, METPhi;
  float PhotonPt, PhotonEta, PhotonPhi;
  float MT;
  float j1_Eta;
  float j2_Eta;
  //selected lepton variables
  int lep1Type = 0;
  int lep1PassSelection = 0;
  float lep1Pt = -999;
  float lep1Eta = -999;
  float lep1Phi = -999;  
  int lep2Type = 0;
  int lep2PassSelection = 0;
  float lep2Pt = -999;
  float lep2Eta = -999;
  float lep2Phi = -999;
  float dileptonMass = -999;
  float lep1MT = -999;
  float lep1GenMetMT = -999;
  //
  int Iso_lepton1;
  int Iso_lepton2;
  // working variables:
  float Copy_Jet_PT[900];
  
  HggRazorLeptonsBox razorbox = None;
  
  TTree *outputTree = new TTree("BkgTree", "Info on selected razor inclusive events");

  outputTree->Branch("weight", &weight, "weight/F");
  outputTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
  outputTree->Branch("pileupWeightUp", &pileupWeightUp, "pileupWeightUp/F");
  outputTree->Branch("pileupWeightDown", &pileupWeightDown, "pileupWeightDown/F");
  outputTree->Branch("triggerEffWeight", &triggerEffWeight, "triggerEffWeight/F");
  outputTree->Branch("triggerEffSFWeight", &triggerEffSFWeight, "triggerEffSFWeight/F");
  outputTree->Branch("photonEffSF", &photonEffSF, "photonEffSF/F");
  outputTree->Branch("run", &run, "run/i");
  outputTree->Branch("lumi", &lumi, "lumi/i");
  outputTree->Branch("event", &event, "event/i");
  outputTree->Branch("NPU", &NPU, "npu/i");
  outputTree->Branch("MET", &MET, "MET/F");
  outputTree->Branch("METPhi", &METPhi, "METPhi/F");
  outputTree->Branch("PhotonPt", &PhotonPt, "PhotonPt/F");
  outputTree->Branch("PhotonEta", &PhotonEta, "PhotonEta/F");
  outputTree->Branch("PhotonPhi", &PhotonPhi, "PhotonPhi/F");
  outputTree->Branch("MT", &MT, "MT/F");
  outputTree->Branch("j1_Eta", &j1_Eta, "j1_Eta/F");
  outputTree->Branch("j2_Eta", &j2_Eta, "j2_Eta/F");
  // lepton stuff
  outputTree->Branch("lep1Type", &lep1Type, "lep1Type/I");
  outputTree->Branch("lep1PassSelection", &lep1PassSelection, "lep1PassSelection/I");
  outputTree->Branch("lep1Pt", &lep1Pt, "lep1Pt/F");
  outputTree->Branch("lep1Eta", &lep1Eta, "lep1Eta/F");
  outputTree->Branch("lep1Phi", &lep1Phi, "lep1Phi/F");
  outputTree->Branch("lep2Type", &lep2Type, "lep2Type/I");
  outputTree->Branch("lep2PassSelection", &lep2PassSelection, "lep2PassSelection/I");
  outputTree->Branch("lep2Pt", &lep2Pt, "lep2Pt/F");
  outputTree->Branch("lep2Eta", &lep2Eta, "lep2Eta/F");
  outputTree->Branch("lep2Phi", &lep2Phi, "lep2Phi/F");
  outputTree->Branch("dileptonMass", &dileptonMass, "dileptonMass/F");
  outputTree->Branch("lep1MT", &lep1MT, "lep1MT/F");
  outputTree->Branch("lep1GenMetMT", &lep1GenMetMT, "lep1GenMetMT/F");
  //
  outputTree->Branch("box", &razorbox, "box/I");
  //
  outputTree->Branch("Iso_lepton1", &Iso_lepton1, "Iso_lepton1/i");	// isolated lepton? (no add. requirements on Eta/PT) 1=YES, 0=NO
  outputTree->Branch("Iso_lepton2", &Iso_lepton2, "Iso_lepton2/i");	// isolated lepton? (WITH add. requirements on Eta/PT) 1=YES, 0=NO
  //begin loop
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    double leadingPhotonPt = -999;
    double leadingPhotonEta = -999;
    double leadingPhotonPhi = -999;

    // define variables for finding the highest PT photon
    int i_best=0;
    
    //Select highest pT photon
    for(int i = 0; i < nPhotons; i++){

      //require medium photon selection
      if (!isMediumPhoton(i)) continue;

      //find the highest pT photon here
      //*******
      if (phoPt[i] > phoPt[i_best]) {
      	  i_best = i;
      }
      //*******
    }

    //Fill MET, MetPhi, MT here
    MET = metPt;
    METPhi = metPhi;
    PhotonPt = phoPt[i_best];
    PhotonEta = phoEta[i_best];
    PhotonPhi = phoPhi[i_best];
    
    MT = sqrt(2*MET*PhotonPt*(1 - cos(PhotonPhi - METPhi)));
    
    //Fill j1, j2 eta values here
    
    // ## Find two jets with highest pT: j1, j2 with j1_pT > j2_pT
	// make array copy of Jet
	for(int i=0; i<900; i++) {
		Copy_Jet_PT[i] = jetPt[i];	
	}
	// find indices
	int j1_index = TMath::LocMax(900,jetPt);	// find j1 index
	Copy_Jet_PT[j1_index] = 0;				// set j1 value in Copy to zero
	int j2_index = TMath::LocMax(900,Copy_Jet_PT);	// find j2 index
	// assign variables
	j1_Eta = jetEta[j1_index];
	j2_Eta = jetEta[j2_index];
	
	//Fill & select lepton PT, Eta and Phi for the best Z->2lepton leptons
	
	float Mz = 91.188;
	razorbox = None;
	//-------------------------------
	//1) Look for Zmm Candidate
	//-------------------------------
	double bestDimuonPt = -1;
	TLorentzVector ZCandidate;
	for( int i = 0; i < nMuons; i++ )	{
		//if(!isVetoMuon(i)) continue;  
		if(muonPt[i] < 20) continue;
		if(abs(muonEta[i]) > 2.4) continue;
		for( int j = i+1; j < nMuons; j++ )	{
			//if(!isVetoMuon(j)) continue;  
			if(muonPt[j] < 20) continue;
			if(abs(muonEta[j]) > 2.4) continue;
		
			TLorentzVector tmpMuon1;
			tmpMuon1.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i],0.1057);
			TLorentzVector tmpMuon2;
			tmpMuon2.SetPtEtaPhiM(muonPt[j],muonEta[j], muonPhi[j],0.1057);
			double tmpMass = (tmpMuon1+tmpMuon2).M();	    
			double tmpDileptonPt = (tmpMuon1+tmpMuon2).Pt();
			
			//if ( _debug ) cout << "Zmm candidate: " << tmpMass << " " << tmpDileptonPt << "\n";
			
			if ( tmpMass > (Mz-15) && tmpMass < (Mz+15) && tmpDileptonPt > bestDimuonPt)  {
				bestDimuonPt = tmpDileptonPt;
				razorbox = Zmm;
				lep1Type = 13 * -1 * muonCharge[i];
				lep1Pt = muonPt[i];
				lep1Eta = muonEta[i];
				lep1Phi = muonPhi[i];
				//lep1PassSelection = 1 + 2 * isTightMuon(i);
				lep2Type = 13 * -1 * muonCharge[j];
				lep2Pt = muonPt[j];
				lep2Eta = muonEta[j];
				lep2Phi = muonPhi[j];
				//lep2PassSelection = 1 + 2 * isTightMuon(j);
				dileptonMass = tmpMass;
				ZCandidate = tmpMuon1 + tmpMuon2;
			
				//for MC apply lepton eff scale factor
				/*if (!isData ) {
					if ( matchesGenMuon(lep1Eta,lep1Phi)) leptonEffSF *=  helper->getVetoMuonScaleFactor( lep1Pt, lep1Eta, true);		
					if ( matchesGenMuon(lep2Eta,lep2Phi)) leptonEffSF *=  helper->getVetoMuonScaleFactor( lep2Pt, lep2Eta, true);			
				}*/
			}
		}
	}


	//-------------------------------
	//2) Look for Zee Candidate
	//-------------------------------
	if (razorbox == None) {
		double bestDielectronPt = -1;
		for( int i = 0; i < nElectrons; i++ )	{
			//if(!isVetoElectron(i)) continue;  
			if(elePt[i] < 20) continue;
			if(abs(eleEta[i]) > 2.4) continue;
			for( int j = i+1; j < nElectrons; j++ )	{
				//if(!isVetoElectron(j)) continue;  
				if(elePt[j] < 20) continue;
				if(abs(eleEta[j]) > 2.4) continue;
			
				TLorentzVector tmpElectron1;
				tmpElectron1.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i],0.000511);
				TLorentzVector tmpElectron2;
				tmpElectron2.SetPtEtaPhiM(elePt[j],eleEta[j], elePhi[j],0.000511);
				double tmpMass = (tmpElectron1+tmpElectron2).M();	    
				double tmpDileptonPt = (tmpElectron1+tmpElectron2).Pt();
		
				//if ( _debug ) cout << "Zee candidate: " << tmpMass << " " << tmpDileptonPt << "\n";
		
				if ( tmpMass > (Mz-15) && tmpMass < (Mz+15) && tmpDileptonPt > bestDielectronPt)  {
					bestDielectronPt = tmpDileptonPt;
					razorbox = Zee;
					lep1Type = 11 * -1 * eleCharge[i];
					lep1Pt = elePt[i];
					lep1Eta = eleEta[i];
					lep1Phi = elePhi[i];
					//lep1PassSelection = 1 + 2 * isTightElectron(i);
					lep2Type = 11 * -1 * eleCharge[j];
					lep2Pt = elePt[j];
					lep2Eta = eleEta[j];
					lep2Phi = elePhi[j];
					//lep2PassSelection = 1 + 2 * isTightElectron(j);
					dileptonMass = tmpMass;
					ZCandidate = tmpElectron1 + tmpElectron2;
		
					//for MC apply lepton eff scale factor
					/*if (!isData ) {
						if ( matchesGenElectron(lep1Eta,lep1Phi)) leptonEffSF *=  helper->getVetoElectronScaleFactor( lep1Pt, lep1Eta, true);		
						if ( matchesGenElectron(lep2Eta,lep2Phi)) leptonEffSF *=  helper->getVetoElectronScaleFactor( lep2Pt, lep2Eta, true);			
					}*/
				}
			}
		}
	}
    
	// Isolated lepton criterion
	Iso_lepton1 = 0;
	Iso_lepton2 = 0;
	Float_t Phi_l;
	Float_t Eta_l;
	Float_t DR;
	Float_t e_Eta;
	Float_t Pi = TMath::Pi();
	for (Int_t j=0; j<700; j++) {	
		Phi_l = TMath::Abs(elePhi[j] - PhotonPhi);
		if (Phi_l > Pi) {
			Phi_l = 2*Pi - Phi_l;
		}
		Eta_l = TMath::Abs(eleEta[j] - PhotonEta);
		DR = sqrt(pow(Phi_l,2) + pow(Eta_l,2));
		if (DR > 0.3) {									// Isolated electron criterion
			Iso_lepton1 = 1;
			e_Eta = TMath::Abs(eleEta[j]);
			if ((elePt[j]>10) && ( ((e_Eta < 1.4442) || (e_Eta > 1.566)) && (e_Eta<2.5))) {
				Iso_lepton2 = 1;
			}
		}
	}
	for (Int_t j=0; j<700; j++) {	
		Phi_l = TMath::Abs(muonPhi[j] - PhotonPhi);
		if (Phi_l > Pi) {
			Phi_l = 2*Pi - Phi_l;
		}
		Eta_l = TMath::Abs(muonEta[j] - PhotonEta);
		DR = sqrt(pow(Phi_l,2) + pow(Eta_l,2));
		if (DR > 0.3) {									// Isolated muon criterion
			Iso_lepton1 = 1;
			if ((muonPt[j]>10) && (TMath::Abs(muonEta[j])<2.1)) {
				Iso_lepton2 = 1;
			}
		}
	}
	
    //fill normalization histogram    
    NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
    weight = genWeight;
    SumWeights->Fill(1.0, weight);    
    run = runNum;
    lumi = lumiNum; 
    event = eventNum; 

    //------------------
    //Pileup reweighting
    //------------------
    pileupWeight = 1.0;
    if( !isData ) {
      //Get number of PU interactions
      for (int i = 0; i < nBunchXing; i++) {
	if (BunchXing[i] == 0) {
	  NPU = nPUmean[i];
	}
      }
      pileupWeight = helper->getPileupWeight(NPU);
      pileupWeightUp = helper->getPileupWeightUp(NPU) / pileupWeight;
      pileupWeightDown = helper->getPileupWeightDown(NPU) / pileupWeight;	
    }
    
    /////////////////////////////////
    //Scale and PDF variations
    /////////////////////////////////
    if( !isData ) {
      if ( (*scaleWeights).size() >= 9 ) 
	{
	  // sf_facScaleUp      = (*scaleWeights)[1]/genWeight;
	  // sf_facScaleDown    = (*scaleWeights)[2]/genWeight;
	  // sf_renScaleUp      = (*scaleWeights)[3]/genWeight;
	  // sf_renScaleDown    = (*scaleWeights)[6]/genWeight;
	  // sf_facRenScaleUp   = (*scaleWeights)[4]/genWeight;
	  // sf_facRenScaleDown = (*scaleWeights)[8]/genWeight;
	    
	    
	  SumScaleWeights->Fill(0.0, (*scaleWeights)[1]);
	  SumScaleWeights->Fill(1.0, (*scaleWeights)[2]);
	  SumScaleWeights->Fill(2.0, (*scaleWeights)[3]);
	  SumScaleWeights->Fill(3.0, (*scaleWeights)[6]);
	  SumScaleWeights->Fill(4.0, (*scaleWeights)[4]);
	  SumScaleWeights->Fill(5.0, (*scaleWeights)[8]);
	}
	
      // sf_pdf.erase( sf_pdf.begin(), sf_pdf.end() );
      for ( unsigned int iwgt = 0; iwgt < pdfWeights->size(); ++iwgt ) {
      // 	  sf_pdf.push_back( pdfWeights->at(iwgt)/genWeight );
	SumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
      }
    }
    
    // Put Trigger stuff here ?
      
    outputTree->Fill();

  } // loop over events 
    
  cout << "Writing output trees..." << endl;
  outFile->Write();
  outFile->Close();
    
}
