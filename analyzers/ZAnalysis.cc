//LOCAL INCLUDES
#include "ZAnalysis.h"
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
#include <TH2D.h>
#include "TRandom3.h"
#include "AngleConversion.h"

using namespace std;


const Double_t MASS_LOW  = 40;
const Double_t MASS_HIGH = 200;
const Double_t PT_CUT    = 22;
const Double_t ETA_CUT   = 2.4;
const Double_t MUON_MASS = 0.105658369;

const Int_t BOSON_ID  = 23;
const Int_t LEPTON_ID = 13;
const int NUM_PDF_WEIGHTS = 60;

enum { eNone=0, eMuMu=1, eEleEle=2 };  // event category enum

//Testing branching and merging
void ZAnalysis::Analyze(bool isData, int option, string outFileName, string label)
{
  //*****************************************************************************
  //Settings
  //*****************************************************************************
  TRandom3 random(3003);

  string analysisTag = "Razor2016_80X";
  if ( label != "") analysisTag = label;

  std::cout << "[INFO]: option = " << option << std::endl;
  std::cout << "[INFO]: analysisTag --> " << analysisTag << std::endl;

  
  if ( outFileName.empty() )
    {
      std::cout << "ZAnalysis: Output filename not specified!" << endl << "Using default output name ZAnalysis.root" << std::endl;
      outFileName = "ZAnalysis.root";
    }
  TFile* outFile = new TFile( outFileName.c_str(), "RECREATE" );
  //---------------------------
  //one tree to hold all events
  //---------------------------
  TTree *outTree = new TTree("ZAnalysis", "Info on selected razor inclusive events");
  
  //Get CMSSW Directory
  // char* cmsswPath;
  // cmsswPath = getenv("CMSSW_BASE");

  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  if (analysisTag == "Razor2015_76X") helper = new RazorHelper("Razor2015_76X", isData, false);
  else if (analysisTag == "Razor2016_80X") helper = new RazorHelper("Razor2016_80X", isData, false);
  else helper = new RazorHelper(analysisTag, isData, false);
  

  //--------------------------------
  //Including Jet Energy Corrections
  //--------------------------------
  std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
  std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();


  //----------
  //pu histo
  //----------
  TH1D* puhisto = new TH1D("pileup", "", 50, 0, 50);
  
  //histogram containing total number of processed events (for normalization)
  TH1F *histNPV = new TH1F("NPV", "NPV", 2, -0.5, 1.5);
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
  TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
  TH1F *SumPdfWeights = new TH1F("SumPdfWeights", "SumPdfWeights", NUM_PDF_WEIGHTS, -0.5, NUM_PDF_WEIGHTS-0.5);

  //--------------
  //tree variables
  //--------------
  int NPU;
  bool matchGen;
  int category;

  // UInt_t  id_1, id_2;
  // Double_t x_1, x_2, xPDF_1, xPDF_2;
  // Double_t scalePDF, weightPDF;

  TLorentzVector *genV=0;
  float genVPt, genVPhi, genVy, genVMass;
  float genWeight;
  float PUWeight, PUWeightUp, PUWeightDown;
  float weight;

  int genLep1Id, genLep2Id;
  TLorentzVector *genLep1;
  TLorentzVector *genLep2;
  TLorentzVector *genZ;

  int q1, q2;
  TLorentzVector *dilep;
  TLorentzVector *lep1;
  TLorentzVector *lep2;

  //PDF SF
  std::vector<float> sf_pdf;
  

  //------------------------
  //set branches on big tree
  //------------------------
 
    outTree->Branch("runNum",      &runNum,     "runNum/i");      // event run number
    outTree->Branch("lumiSec",     &lumiNum,    "lumiSec/i");     // event lumi section
    outTree->Branch("evtNum",      &eventNum,     "evtNum/i");      // event number
    outTree->Branch("matchGen",    &matchGen,   "matchGen/i");    // event has both leptons matched to MC Z->ll
    outTree->Branch("category",    &category,   "category/i");    // dilepton category
    // outTree->Branch("id_1",        &id_1,       "id_1/i");        // PDF info -- parton ID for parton 1
    // outTree->Branch("id_2",        &id_2,       "id_2/i");        // PDF info -- parton ID for parton 2
    // outTree->Branch("x_1",         &x_1,        "x_1/d");         // PDF info -- x for parton 1
    // outTree->Branch("x_2",         &x_2,        "x_2/d");         // PDF info -- x for parton 2
    // outTree->Branch("xPDF_1",      &xPDF_1,     "xPDF_1/d");      // PDF info -- x*F for parton 1
    // outTree->Branch("xPDF_2",      &xPDF_2,     "xPDF_2/d");      // PDF info -- x*F for parton 2
    // outTree->Branch("scalePDF",    &scalePDF,   "scalePDF/d");    // PDF info -- energy scale of parton interaction
    // outTree->Branch("weightPDF",   &weightPDF,  "weightPDF/d");   // PDF info -- PDF weight
    outTree->Branch("npv",         &nPV,        "npv/i");         // number of primary vertices
    outTree->Branch("npu",         &NPU,        "npu/i");         // number of in-time PU events (MC)
    outTree->Branch("genWeight",   &genWeight,  "genWeight/F");
    outTree->Branch("PUWeight",    &PUWeight,   "PUWeight/F");
    outTree->Branch("PUWeightUp",  &PUWeightUp, "PUWeightUp/F");
    outTree->Branch("PUWeightDown",  &PUWeightDown, "PUWeightDown/F");
    outTree->Branch("weight",      &weight,     "weight/F");      // event weight per 1/fb (MC)

    outTree->Branch("genLep1Id",   &genLep1Id,  "genLep1Id/I");  
    outTree->Branch("genLep2Id",   &genLep2Id,  "genLep2Id/I");  
    outTree->Branch("genZ",        "TLorentzVector", &genZ);      // di-lepton 4-vector
    outTree->Branch("genLep1",     "TLorentzVector", &genLep1);   // tag lepton 4-vector
    outTree->Branch("genLep2",     "TLorentzVector", &genLep2);   // tag lepton 4-vector

    outTree->Branch("q1",          &q1,         "q1/I");          // charge of tag lepton
    outTree->Branch("q2",          &q2,         "q2/I");          // charge of probe lepton
    outTree->Branch("dilep",       "TLorentzVector", &dilep);     // tag lepton 4-vector
    outTree->Branch("lep1",        "TLorentzVector", &lep1);      // tag lepton 4-vector
    outTree->Branch("lep2",        "TLorentzVector", &lep2);      // tag lepton 4-vector

  //begin loop
  if ( fChain == 0 ) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  std::cout << "[INFO]: Total Entries = " << fChain->GetEntries() << "\n";
  for ( Long64_t jentry=0; jentry < nentries; jentry++ ) {
      //begin event
      if( jentry % 10000 == 0 ) std::cout << "[INFO]: Processing entry " << jentry << std::endl;
      Long64_t ientry = LoadTree( jentry );
      if ( ientry < 0 ) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
    
      //fill normalization histogram    
      NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
      weight = genWeight;
      SumWeights->Fill(1.0, weight);
      
      //reset tree variables
      runNum = -999;
      lumiNum = -999;
      eventNum = -999;
      matchGen = -999;
      category = -999;
      // id_1 = -999;
      // id_2 = -999;
      // x_1 = -999;
      // x_2 = -999;
      // xPDF_1 = -999;
      // xPDF_2 = -999;
      // scalePDF = -999;
      // weightPDF = -999;
      nPV = -999;
      NPU = -999;
      genWeight = -999;
      PUWeight = -999;
      weight = 0;
      genLep1Id = -999;
      genLep2Id = -999;
      genZ = 0;
      genLep1 = 0;
      genLep2 = 0;
      q1 = -999;
      q2 = -999;
      dilep = 0;
      lep1 = 0;
      lep2 = 0;



      //------------------
      //Pileup reweighting
      //------------------
      PUWeight = 1.0;
      if( !isData ) {
	//Get number of PU interactions
	for (int i = 0; i < nBunchXing; i++) {
	  if (BunchXing[i] == 0) {
	    NPU = nPUmean[i];
	  }
	}
	puhisto->Fill(NPU);
	PUWeight = helper->getPileupWeight(NPU);
	PUWeightUp = helper->getPileupWeightUp(NPU) / PUWeight;
	PUWeightDown = helper->getPileupWeightDown(NPU) / PUWeight;	
      }
      
      /////////////////////////////////
      //Scale and PDF variations
      /////////////////////////////////

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
      
      sf_pdf.erase( sf_pdf.begin(), sf_pdf.end() );
      for ( unsigned int iwgt = 0; iwgt < pdfWeights->size(); ++iwgt ) {
	sf_pdf.push_back( pdfWeights->at(iwgt)/genWeight );
	SumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
      }
      
      //*************************************************************************
      //Gen-Level information
      //*************************************************************************
      double genZPt, genZEta, genZPhi, genZEnergy;
      double genLep1Pt, genLep1Eta, genLep1Phi, genLep1Energy;
      double genLep2Pt, genLep2Eta, genLep2Phi, genLep2Energy;
      for(int j = 0; j < nGenParticle; j++){
	// cout << "particle: " << j << " : " << gParticleId[j] << " " << gParticleStatus[j] << " : " 
	//      << gParticlePt[j] << " " 
	//      << gParticleEta[j] << " " 
	//      << gParticlePhi[j] << " " 
	//      << " | " << gParticleMotherId[j] << "\n";

	//Z Boson
	if (gParticleStatus[j] == 22 && gParticleId[j] == 23) {
	  genZPt = gParticlePt[j];
	  genZEta = gParticleEta[j];	 
	  genZPhi = gParticlePhi[j];
	  genZEnergy = gParticleE[j];
	}

	//if leptons
	if ( (abs(gParticleId[j]) == 11 || abs(gParticleId[j]) == 13 || abs(gParticleId[j]) == 15)
	     && gParticleStatus[j] == 1
	     ) {	   	    

	  //ZLepton1
	  if (gParticleMotherId[j] == 23 && gParticleId[j] > 0) {
	    genLep1Id = gParticleId[j];
	    genLep1Pt = gParticlePt[j];
	    genLep1Eta = gParticleEta[j];	   
	    genLep1Phi = gParticlePhi[j];	   
	    genLep1Energy = gParticleE[j];	   
	  }
	  //ZLepton2
	  if (gParticleMotherId[j] == 23 && gParticleId[j] < 0) {
	    genLep2Id = gParticleId[j];
	    genLep2Pt = gParticlePt[j];
	    genLep2Eta = gParticleEta[j];	   
	    genLep2Phi = gParticlePhi[j];	   
	    genLep2Energy = gParticleE[j];   
	  }
	 

	} // endif leptons

      }//loop over gen particles




      //*************************************************************************
      //Start Object Selection
      //*************************************************************************

      //-------------------------------
      //1) Look for Zmm Candidate
      //-------------------------------
      int lep1Type = 0;
      double lep1Pt = 0;
      double lep1Eta = -999;
      double lep1Phi = -999;
      int lep2Type = 0;
      double lep2Pt = 0;
      double lep2Eta = -999;
      double lep2Phi = -999;

      double bestDiMuonMass = 0;
      TLorentzVector ZCandidate;
      int NumberOfVetoLeptons = 0;
      for( int i = 0; i < nMuons; i++ )	{

	// Count number of Veto Leptons so that events with 
	// more than 2 veto leptons can be rejected
	if( isMuonPOGLooseMuon(i, true, false) && 
	    ((muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.15)
	    ) NumberOfVetoLeptons++;
	
        //Muon Selection Requirements
	if(!( isMuonPOGMediumMuon(i, true, false) && 
	      ((muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i] < 0.15)	  
	      )) continue;  
	if(muonPt[i] < 20) continue;
	if(abs(muonEta[i]) > 2.4) continue;
	
	for( int j = i+1; j < nMuons; j++ )	{	  
	  if(!( isMuonPOGMediumMuon(j, true, false) && 
		((muon_chargedIso[j] + fmax(0.0,  muon_photonIso[j] + muon_neutralHadIso[j] - 0.5*muon_pileupIso[j])) / muonPt[j] < 0.15)
		)) continue;  
	  if(muonPt[j] < 20) continue;
	  if(abs(muonEta[j]) > 2.4) continue;
	
	  TLorentzVector tmpMuon1;
	  tmpMuon1.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i],0.1057);
	  TLorentzVector tmpMuon2;
	  tmpMuon2.SetPtEtaPhiM(muonPt[j],muonEta[j], muonPhi[j],0.1057);
	  double tmpMass = (tmpMuon1+tmpMuon2).M();	    	  
	  
	  if ( muonPt[j] > lep2Pt ) {
	    bestDiMuonMass = tmpMass;
	    category = eMuMu;
	    lep1Type = 13 * -1 * muonCharge[i];
	    lep1Pt = muonPt[i];
	    lep1Eta = muonEta[i];
	    lep1Phi = muonPhi[i];
	    lep2Type = 13 * -1 * muonCharge[j];
	    lep2Pt = muonPt[j];
	    lep2Eta = muonEta[j];
	    lep2Phi = muonPhi[j];
	    ZCandidate = tmpMuon1 + tmpMuon2;
	    matchGen = matchesGenMuon(lep1Eta,lep1Phi) && matchesGenMuon(lep2Eta,lep2Phi);	      
	  } //if better Z candidate
	}// loop 2nd muon
      } //loop 1st muon
      

      //-------------------------------
      //2) Look for Zee Candidate
      //-------------------------------
      if (category == eNone) {
	double bestDielectronMass = 0;
	for( int i = 0; i < nElectrons; i++ )	{
	  if (isEGammaPOGVetoElectron(i,true,true,true,"Summer16")) NumberOfVetoLeptons++;

	  if(!isEGammaPOGMediumElectron(i,true,true,true,"Summer16")) continue;  
	  if(elePt[i] < 20) continue;
	  if(abs(eleEta[i]) > 2.5) continue;
	  for( int j = i+1; j < nElectrons; j++ )	{
	    if(!isEGammaPOGMediumElectron(i,true,true,true,"Summer16")) continue; 
	    if(elePt[j] < 20) continue;
	    if(abs(eleEta[j]) > 2.5) continue;
	    
	    TLorentzVector tmpElectron1;
	    tmpElectron1.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i],0.000511);
	    TLorentzVector tmpElectron2;
	    tmpElectron2.SetPtEtaPhiM(elePt[j],eleEta[j], elePhi[j],0.000511);
	    double tmpMass = (tmpElectron1+tmpElectron2).M();	    

	    if ( elePt[j] > lep2Pt ) {
	      bestDielectronMass = tmpMass;
	      category = eEleEle;
	      lep1Type = 11 * -1 * eleCharge[i];
	      lep1Pt = elePt[i];
	      lep1Eta = eleEta[i];
	      lep1Phi = elePhi[i];
	      lep2Type = 11 * -1 * eleCharge[j];
	      lep2Pt = elePt[j];
	      lep2Eta = eleEta[j];
	      lep2Phi = elePhi[j];
	      ZCandidate = tmpElectron1 + tmpElectron2;
	      matchGen = matchesGenElectron(lep1Eta,lep1Phi) && matchesGenElectron(lep2Eta,lep2Phi);	      
	    } //if better Z mass match
	  } //loop 2nd electron
	} //loop 1st electron
      }

      //-------------------------------
      // Save Reco Lepton Information
      //-------------------------------
      lep1 = new TLorentzVector;
      lep2 = new TLorentzVector;
      dilep = new TLorentzVector;
      if ( abs(lep1Type) == 11) {
	lep1->SetPtEtaPhiM(lep1Pt, lep1Eta, lep1Phi,0.000511);
	lep2->SetPtEtaPhiM(lep2Pt, lep2Eta, lep2Phi,0.000511);
      } else {
	lep1->SetPtEtaPhiM(lep1Pt, lep1Eta, lep1Phi,0.1057);
	lep2->SetPtEtaPhiM(lep2Pt, lep2Eta, lep2Phi,0.1057);
      }
      dilep->SetPtEtaPhiM( (*lep1+*lep2).Pt(), (*lep1+*lep2).Eta(), (*lep1+*lep2).Phi(), (*lep1+*lep2).M());
      q1 = -1 * lep1Type / abs(lep1Type);
      q2 = -1 * lep2Type / abs(lep2Type);


    
      //Fill Event
      outTree->Fill();

      //end of event loop
  }
  
  std::cout << "[INFO]: Number of events processed: " << NEvents->Integral() << std::endl;
  
  std::cout << "[INFO]: Writing output trees..." << std::endl;    
  outFile->cd();
  outTree->Write();
  NEvents->Write();
  SumWeights->Write();
  SumScaleWeights->Write();
  SumPdfWeights->Write();
  histNPV->Write();
  puhisto->Write();
  
  outFile->Close();
  delete helper;

}
