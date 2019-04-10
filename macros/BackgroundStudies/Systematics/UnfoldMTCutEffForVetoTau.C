//================================================================================================
//
// Simple Example
//
//root -l RazorAnalyzer/macros/ObjectStudies/MakeElectronEfficiencyPlots.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/ElectronNtuple/ElectronNtuple_PromptGenLevel_TTJets_25ns.root",-1,"Electron")'
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TH1D.h>                
#include <TH2F.h>                
#include <TCanvas.h>                
#include <TLegend.h> 
#include <THStack.h> 

#include "RazorAnalyzer/include/ControlSampleEvents.h"
#include "RazorAnalyzer/macros/tdrstyle.C"
#include "RazorAnalyzer/macros/CMS_lumi.C"
#include "RazorAnalyzer/include/RecoilCorrector.hh"

#endif



//=== MAIN MACRO ================================================================================================= 


void RunUnfoldMTCutEffForVetoTaus( vector<string> datafiles, vector<vector<string> > bkgfiles, 
				      vector<string> bkgLabels, vector<int> bkgColors, double lumi) {
  
  // string Label = "";
  // if (label != "") Label = "_" + label;

  //*******************************************************************
  // Settings 
  //*******************************************************************
  float MR = 0;
  float Rsq = 0;
  float mll = 0;

  bool printdebug = false;

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  const int NMRBins = 10;
  const int NRsqBins = 9;
  const int NLepPtBins = 8;
  const int NLepEtaBins = 6;
  double MRBins[NMRBins] = {300, 350, 400, 450, 500, 550, 700, 900, 1200, 4000};
  double RsqBins[NRsqBins] = {0.15,0.175,0.20,0.225, 0.25,0.30,0.41,0.52,1.5};  
  double LepPtBins[NLepPtBins] = {5,10,15,20,30,40,100,1000};  
  double LepEtaBins[NLepEtaBins] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

  assert ( bkgfiles.size() == bkgLabels.size() );
  assert ( bkgfiles.size() == bkgColors.size() );

  vector<vector<string> > inputfiles;
  vector<string> processLabels;
  vector<int> color;

  inputfiles.push_back(datafiles);
  processLabels.push_back("Data");
  color.push_back(kBlack);

  assert(bkgfiles.size() == bkgLabels.size());
  assert(bkgfiles.size() == bkgColors.size());
  for (int i=0; i < int(bkgfiles.size()); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processLabels.push_back(bkgLabels[i]);
     color.push_back(bkgColors[i]);
  }

 

  TH1D *histLep1Pt = new TH1D("histLep1Pt", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1Pt_FScaleUp = new TH1D("histLep1Pt_FScaleUp", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1Pt_FScaleDown = new TH1D("histLep1Pt_FScaleDown", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1Pt_RScaleUp = new TH1D("histLep1Pt_RScaleUp", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1Pt_RScaleDown = new TH1D("histLep1Pt_RScaleDown", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);

  TH1D *histLep1PtPassMTCut = new TH1D("histLep1PtPassMTCut", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassMTCut_JESUp = new TH1D("histLep1PtPassMTCut_JESUp", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassMTCut_JESDown = new TH1D("histLep1PtPassMTCut_JESDown", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassMTCut_LESUp = new TH1D("histLep1PtPassMTCut_LESUp", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassMTCut_LESDown = new TH1D("histLep1PtPassMTCut_LESDown", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassMTCut_FScaleUp = new TH1D("histLep1PtPassMTCut_FScaleUp", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassMTCut_FScaleDown = new TH1D("histLep1PtPassMTCut_FScaleDown", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassMTCut_RScaleUp = new TH1D("histLep1PtPassMTCut_RScaleUp", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassMTCut_RScaleDown = new TH1D("histLep1PtPassMTCut_RScaleDown", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);


  TH1D *histLep1Eta = new TH1D("histLep1Eta", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1Eta_FScaleUp = new TH1D("histLep1Eta_FScaleUp", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1Eta_FScaleDown = new TH1D("histLep1Eta_FScaleDown", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1Eta_RScaleUp = new TH1D("histLep1Eta_RScaleUp", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1Eta_RScaleDown = new TH1D("histLep1Eta_RScaleDown", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);

  TH1D *histLep1EtaPassMTCut = new TH1D("histLep1EtaPassMTCut", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassMTCut_JESUp = new TH1D("histLep1EtaPassMTCut_JESUp", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassMTCut_JESDown = new TH1D("histLep1EtaPassMTCut_JESDown", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassMTCut_LESUp = new TH1D("histLep1EtaPassMTCut_LESUp", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassMTCut_LESDown = new TH1D("histLep1EtaPassMTCut_LESDown", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassMTCut_FScaleUp = new TH1D("histLep1EtaPassMTCut_FScaleUp", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassMTCut_FScaleDown = new TH1D("histLep1EtaPassMTCut_FScaleDown", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassMTCut_RScaleUp = new TH1D("histLep1EtaPassMTCut_RScaleUp", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassMTCut_RScaleDown = new TH1D("histLep1EtaPassMTCut_RScaleDown", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
 
  double dataYield = 0;
  double MCYield = 0;
  double MCTTBarYield = 0;

  vector<pair<UInt_t,UInt_t> > RunAndEvent;

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  for (uint i=0; i < inputfiles.size(); ++i) {
    for (uint j=0; j < inputfiles[i].size(); ++j) {

      TFile *f = TFile::Open(inputfiles[i][j].c_str(), "READ");
      assert(f);
      TTree *tree = (TTree*)f->Get("RazorInclusive");
      assert(tree);

      float weight = 0;
      int nBTaggedJets = 0;
      int nSelectedJets = 0;
      int nJets80 = 0;
      int nVetoMuons = -1;
      int nVetoElectrons = -1;
      int nLooseTaus = -1;
      float leadingGenLeptonPt = 0;
      float leadingGenLeptonEta = 0;
      int leadingGenLeptonType = 0;
      float MR = 0;
      float Rsq = 0;
      float mTLoose = 0;
      float mTLoose_JESUp = 0;
      float mTLoose_JESDown = 0;
      float mTLoose_EESUp = 0;
      float mTLoose_EESDown = 0;
      float mTLoose_MESUp = 0;
      float mTLoose_MESDown = 0;
      bool  HLTDecision[150];
      float sf_facScaleUp = 0;
      float sf_facScaleDown = 0;
      float sf_renScaleUp = 0;
      float sf_renScaleDown = 0;

      tree->SetBranchStatus("*", 0);
      tree->SetBranchStatus("weight", 1);
      tree->SetBranchStatus("nBTaggedJets",1);
      tree->SetBranchStatus("nSelectedJets", 1);
      tree->SetBranchStatus("nJets80", 1);
      tree->SetBranchStatus("nVetoMuons", 1);
      tree->SetBranchStatus("nVetoElectrons", 1);
      tree->SetBranchStatus("nLooseTaus", 1);
      tree->SetBranchStatus("leadingGenLeptonPt", 1);
      tree->SetBranchStatus("leadingGenLeptonEta", 1);
      tree->SetBranchStatus("leadingGenLeptonType", 1);
      tree->SetBranchStatus("MR", 1);
      tree->SetBranchStatus("Rsq", 1);
      tree->SetBranchStatus("mTLoose", 1);
      tree->SetBranchStatus("mTLoose_JESUp", 1);
      tree->SetBranchStatus("mTLoose_JESDown", 1);
      tree->SetBranchStatus("mTLoose_EESUp", 1);
      tree->SetBranchStatus("mTLoose_EESDown", 1);
      tree->SetBranchStatus("mTLoose_MESUp", 1);
      tree->SetBranchStatus("mTLoose_MESDown", 1);
      tree->SetBranchStatus("HLTDecision", 1);
      tree->SetBranchStatus("sf_facScaleUp", 1);
      tree->SetBranchStatus("sf_facScaleDown", 1);
      tree->SetBranchStatus("sf_renScaleUp", 1);
      tree->SetBranchStatus("sf_renScaleDown", 1);
      
      tree->SetBranchAddress("weight",&weight);
      tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
      tree->SetBranchAddress("nSelectedJets",&nSelectedJets);
      tree->SetBranchAddress("nJets80",&nJets80);
      tree->SetBranchAddress("nVetoMuons",&nVetoMuons);
      tree->SetBranchAddress("nVetoElectrons",&nVetoElectrons);
      tree->SetBranchAddress("nLooseTaus",&nLooseTaus);
      tree->SetBranchAddress("leadingGenLeptonPt",&leadingGenLeptonPt);
      tree->SetBranchAddress("leadingGenLeptonEta",&leadingGenLeptonEta);
      tree->SetBranchAddress("leadingGenLeptonType",&leadingGenLeptonType);
      tree->SetBranchAddress("MR",&MR);
      tree->SetBranchAddress("Rsq",&Rsq);
      tree->SetBranchAddress("mTLoose",&mTLoose);
      tree->SetBranchAddress("mTLoose_JESUp",&mTLoose_JESUp);
      tree->SetBranchAddress("mTLoose_JESDown",&mTLoose_JESDown);
      tree->SetBranchAddress("mTLoose_EESUp",&mTLoose_EESUp);
      tree->SetBranchAddress("mTLoose_EESDown",&mTLoose_EESDown);
      tree->SetBranchAddress("mTLoose_MESUp",&mTLoose_MESUp);
      tree->SetBranchAddress("mTLoose_MESDown",&mTLoose_MESDown);
      tree->SetBranchAddress("HLTDecision",&HLTDecision);
      tree->SetBranchAddress("sf_facScaleUp",&sf_facScaleUp);
      tree->SetBranchAddress("sf_facScaleDown",&sf_facScaleDown);
      tree->SetBranchAddress("sf_renScaleUp",&sf_renScaleUp);
      tree->SetBranchAddress("sf_renScaleDown",&sf_renScaleDown);

      bool isData = false;
      if ( processLabels[i] == "Data") isData = true;

      cout << "process: " << processLabels[i] << " | file " << inputfiles[i][j] << " | Total Entries: " << tree->GetEntries() << "\n";
      for(UInt_t ientry=0; ientry < tree->GetEntries(); ientry++) {       	
	tree->GetEntry(ientry);
      
	if (ientry % 100000 == 0) cout << "Event " << ientry << endl;      

	//******************************
	//Trigger Selection
	//******************************
	bool passTrigger = false;

	//Razor Hadronic Triggers
	if ( HLTDecision[134] || HLTDecision[135] || HLTDecision[136] 
	     || HLTDecision[137] || HLTDecision[138]
	     || HLTDecision[139] || HLTDecision[140]
	     || HLTDecision[141] || HLTDecision[142]
	     || HLTDecision[143] || HLTDecision[144]
	     ) passTrigger = true;
	if (!passTrigger) continue;
	
	//******************************
	//Selection Cuts 
	//******************************
	if (!( nLooseTaus >= 1 )) continue;
	if (!(nJets80 >= 2)) continue;
	if (!(nSelectedJets >= 4)) continue;
	if (!(MR > 400 && Rsq > 0.25)) continue;
	if (!TMath::Finite(weight)) continue;
 
	//******************************
	//Fill histograms
	//******************************
	if (isData) {
	} else {
	  //cout << leadingGenLeptonPt << " " << weight << "\n";
	  histLep1Pt->Fill(leadingGenLeptonPt, weight*lumi);
	  histLep1Pt_FScaleUp->Fill(leadingGenLeptonPt, sf_facScaleUp*weight*lumi);
	  histLep1Pt_FScaleDown->Fill(leadingGenLeptonPt, sf_facScaleDown*weight*lumi);
	  histLep1Pt_RScaleUp->Fill(leadingGenLeptonPt, sf_renScaleUp*weight*lumi);
	  histLep1Pt_RScaleDown->Fill(leadingGenLeptonPt, sf_renScaleDown*weight*lumi);
	  histLep1Eta->Fill(fabs(leadingGenLeptonEta), weight*lumi);
	  histLep1Eta_FScaleUp->Fill(fabs(leadingGenLeptonEta), sf_facScaleUp*weight*lumi);
	  histLep1Eta_FScaleDown->Fill(fabs(leadingGenLeptonEta), sf_facScaleDown*weight*lumi);
	  histLep1Eta_RScaleUp->Fill(fabs(leadingGenLeptonEta), sf_renScaleUp*weight*lumi);
	  histLep1Eta_RScaleDown->Fill(fabs(leadingGenLeptonEta), sf_renScaleDown*weight*lumi);
	  if (mTLoose > 30 && mTLoose < 100) {
	    histLep1PtPassMTCut->Fill(leadingGenLeptonPt, weight*lumi);
	    histLep1PtPassMTCut_FScaleUp->Fill(leadingGenLeptonPt, sf_facScaleUp*weight*lumi);
	    histLep1PtPassMTCut_FScaleDown->Fill(leadingGenLeptonPt, sf_facScaleDown*weight*lumi);
	    histLep1PtPassMTCut_RScaleUp->Fill(leadingGenLeptonPt, sf_renScaleUp*weight*lumi);
	    histLep1PtPassMTCut_RScaleDown->Fill(leadingGenLeptonPt, sf_renScaleDown*weight*lumi);
	    histLep1EtaPassMTCut->Fill(fabs(leadingGenLeptonEta), weight*lumi);
	    histLep1EtaPassMTCut_FScaleUp->Fill(fabs(leadingGenLeptonEta), sf_facScaleUp*weight*lumi);
	    histLep1EtaPassMTCut_FScaleDown->Fill(fabs(leadingGenLeptonEta), sf_facScaleDown*weight*lumi);
	    histLep1EtaPassMTCut_RScaleUp->Fill(fabs(leadingGenLeptonEta), sf_renScaleUp*weight*lumi);
	    histLep1EtaPassMTCut_RScaleDown->Fill(fabs(leadingGenLeptonEta), sf_renScaleDown*weight*lumi);
	  }
	  if (mTLoose_JESUp > 30 && mTLoose_JESUp < 100) {
	    histLep1PtPassMTCut_JESUp->Fill(leadingGenLeptonPt, weight*lumi);
	    histLep1EtaPassMTCut_JESUp->Fill(fabs(leadingGenLeptonEta), weight*lumi);
	  }
	  if (mTLoose_JESDown > 30 && mTLoose_JESDown < 100) {
	    histLep1PtPassMTCut_JESDown->Fill(leadingGenLeptonPt, weight*lumi);
	    histLep1EtaPassMTCut_JESDown->Fill(fabs(leadingGenLeptonEta), weight*lumi);
	  }	  
	  
	}
      }

      f->Close();

    }
  }
 



  //--------------------------------------------------------------------------------------------------------------
  // Make Plots
  //==============================================================================================================


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open("VetoTauMTCutEfficiency.root", "RECREATE");
  file->cd();
 
  file->WriteTObject(histLep1Pt, "histLep1Pt", "WriteDelete");
  file->WriteTObject(histLep1Pt_FScaleUp, "histLep1Pt_FScaleUp","WriteDelete");
  file->WriteTObject(histLep1Pt_FScaleDown, "histLep1Pt_FScaleDown","WriteDelete");
  file->WriteTObject(histLep1Pt_RScaleUp, "histLep1Pt_RScaleUp","WriteDelete");
  file->WriteTObject(histLep1Pt_RScaleDown, "histLep1Pt_RScaleDown","WriteDelete");

  file->WriteTObject(histLep1PtPassMTCut, "histLep1PtPassMTCut", "WriteDelete");
  file->WriteTObject(histLep1PtPassMTCut_JESUp, "histLep1PtPassMTCut_JESUp", "WriteDelete");
  file->WriteTObject(histLep1PtPassMTCut_JESDown, "histLep1PtPassMTCut_JESDOwn", "WriteDelete");
  file->WriteTObject(histLep1PtPassMTCut_FScaleUp, "histLep1PtPassMTCut_FScaleUp", "WriteDelete");
  file->WriteTObject(histLep1PtPassMTCut_FScaleDown, "histLep1PtPassMTCut_FScaleDown", "WriteDelete");
  file->WriteTObject(histLep1PtPassMTCut_RScaleUp, "histLep1PtPassMTCut_RScaleUp", "WriteDelete");
  file->WriteTObject(histLep1PtPassMTCut_RScaleDown, "histLep1PtPassMTCut_RScaleDown", "WriteDelete");

  file->WriteTObject(histLep1Eta, "histLep1Eta", "WriteDelete");
  file->WriteTObject(histLep1Eta_FScaleUp, "histLep1Eta_FScaleUp","WriteDelete");
  file->WriteTObject(histLep1Eta_FScaleDown, "histLep1Eta_FScaleDown","WriteDelete");
  file->WriteTObject(histLep1Eta_RScaleUp, "histLep1Eta_RScaleUp","WriteDelete");
  file->WriteTObject(histLep1Eta_RScaleDown, "histLep1Eta_RScaleDown","WriteDelete");

  file->WriteTObject(histLep1EtaPassMTCut, "histLep1EtaPassMTCut", "WriteDelete");
  file->WriteTObject(histLep1EtaPassMTCut_JESUp, "histLep1EtaPassMTCut_JESUp", "WriteDelete");
  file->WriteTObject(histLep1EtaPassMTCut_JESDown, "histLep1EtaPassMTCut_JESDOwn", "WriteDelete");
  file->WriteTObject(histLep1EtaPassMTCut_FScaleUp, "histLep1EtaPassMTCut_FScaleUp", "WriteDelete");
  file->WriteTObject(histLep1EtaPassMTCut_FScaleDown, "histLep1EtaPassMTCut_FScaleDown", "WriteDelete");
  file->WriteTObject(histLep1EtaPassMTCut_RScaleUp, "histLep1EtaPassMTCut_RScaleUp", "WriteDelete");
  file->WriteTObject(histLep1EtaPassMTCut_RScaleDown, "histLep1EtaPassMTCut_RScaleDown", "WriteDelete");

  file->Close();
  delete file;       
    

}



void RunUnfoldDPhiCutEffForVetoTaus( vector<string> datafiles, vector<vector<string> > bkgfiles, 
			  vector<string> bkgLabels, vector<int> bkgColors, double lumi) {
  
  // string Label = "";
  // if (label != "") Label = "_" + label;

  //*******************************************************************
  // Settings 
  //*******************************************************************
  float MR = 0;
  float Rsq = 0;
  float mll = 0;

  bool printdebug = false;

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  const int NMRBins = 10;
  const int NRsqBins = 9;
  const int NLepPtBins = 8;
  const int NLepEtaBins = 6;
  double MRBins[NMRBins] = {300, 350, 400, 450, 500, 550, 700, 900, 1200, 4000};
  double RsqBins[NRsqBins] = {0.15,0.175,0.20,0.225, 0.25,0.30,0.41,0.52,1.5};  
  double LepPtBins[NLepPtBins] = {5,10,15,20,30,40,100,1000};  
  double LepEtaBins[NLepEtaBins] = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 };

  assert ( bkgfiles.size() == bkgLabels.size() );
  assert ( bkgfiles.size() == bkgColors.size() );

  vector<vector<string> > inputfiles;
  vector<string> processLabels;
  vector<int> color;

  inputfiles.push_back(datafiles);
  processLabels.push_back("Data");
  color.push_back(kBlack);

  assert(bkgfiles.size() == bkgLabels.size());
  assert(bkgfiles.size() == bkgColors.size());
  for (int i=0; i < int(bkgfiles.size()); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processLabels.push_back(bkgLabels[i]);
     color.push_back(bkgColors[i]);
  }

  TH1D *histLep1Pt = new TH1D("histLep1Pt", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1Pt_FScaleUp = new TH1D("histLep1Pt_FScaleUp", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1Pt_FScaleDown = new TH1D("histLep1Pt_FScaleDown", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1Pt_RScaleUp = new TH1D("histLep1Pt_RScaleUp", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1Pt_RScaleDown = new TH1D("histLep1Pt_RScaleDown", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);

  TH1D *histLep1PtPassDPhiCut = new TH1D("histLep1PtPassDPhiCut", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassDPhiCut_JESUp = new TH1D("histLep1PtPassDPhiCut_JESUp", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassDPhiCut_JESDown = new TH1D("histLep1PtPassDPhiCut_JESDown", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassDPhiCut_LESUp = new TH1D("histLep1PtPassDPhiCut_LESUp", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassDPhiCut_LESDown = new TH1D("histLep1PtPassDPhiCut_LESDown", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassDPhiCut_FScaleUp = new TH1D("histLep1PtPassDPhiCut_FScaleUp", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassDPhiCut_FScaleDown = new TH1D("histLep1PtPassDPhiCut_FScaleDown", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassDPhiCut_RScaleUp = new TH1D("histLep1PtPassDPhiCut_RScaleUp", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
  TH1D *histLep1PtPassDPhiCut_RScaleDown = new TH1D("histLep1PtPassDPhiCut_RScaleDown", "; Lepton p_{T} [GeV/c] ; Number of Events", NLepPtBins-1, LepPtBins);
 
  TH1D *histLep1Eta = new TH1D("histLep1Eta", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1Eta_FScaleUp = new TH1D("histLep1Eta_FScaleUp", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1Eta_FScaleDown = new TH1D("histLep1Eta_FScaleDown", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1Eta_RScaleUp = new TH1D("histLep1Eta_RScaleUp", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1Eta_RScaleDown = new TH1D("histLep1Eta_RScaleDown", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);

  TH1D *histLep1EtaPassDPhiCut = new TH1D("histLep1EtaPassDPhiCut", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassDPhiCut_JESUp = new TH1D("histLep1EtaPassDPhiCut_JESUp", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassDPhiCut_JESDown = new TH1D("histLep1EtaPassDPhiCut_JESDown", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassDPhiCut_LESUp = new TH1D("histLep1EtaPassDPhiCut_LESUp", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassDPhiCut_LESDown = new TH1D("histLep1EtaPassDPhiCut_LESDown", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassDPhiCut_FScaleUp = new TH1D("histLep1EtaPassDPhiCut_FScaleUp", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassDPhiCut_FScaleDown = new TH1D("histLep1EtaPassDPhiCut_FScaleDown", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassDPhiCut_RScaleUp = new TH1D("histLep1EtaPassDPhiCut_RScaleUp", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
  TH1D *histLep1EtaPassDPhiCut_RScaleDown = new TH1D("histLep1EtaPassDPhiCut_RScaleDown", "; Lepton #eta; Number of Events", NLepEtaBins-1, LepEtaBins);
 
  double dataYield = 0;
  double MCYield = 0;
  double MCTTBarYield = 0;

  vector<pair<UInt_t,UInt_t> > RunAndEvent;

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  for (uint i=0; i < inputfiles.size(); ++i) {
    for (uint j=0; j < inputfiles[i].size(); ++j) {

      TFile *f = TFile::Open(inputfiles[i][j].c_str(), "READ");
      assert(f);
      TTree *tree = (TTree*)f->Get("RazorInclusive");
      assert(tree);

      float weight = 0;
      int box = 0;
      int nBTaggedJets = 0;
      int nSelectedJets = 0;
      int nJets80 = 0;
      int nVetoMuons = -1;
      int nVetoElectrons = -1;
      int nLooseTaus = -1;
      float leadingGenLeptonPt = 0;
      float leadingGenLeptonEta = 0;
      int leadingGenLeptonType = 0;
      float MR = 0;
      float Rsq = 0;
      float dPhiRazor = 0;
      float dPhiRazor_JESUp = 0;
      float dPhiRazor_JESDown = 0;
      float dPhiRazor_EESUp = 0;
      float dPhiRazor_EESDown = 0;
      float dPhiRazor_MESUp = 0;
      float dPhiRazor_MESDown = 0;
      bool  HLTDecision[150];
      float sf_facScaleUp = 0;
      float sf_facScaleDown = 0;
      float sf_renScaleUp = 0;
      float sf_renScaleDown = 0;

      tree->SetBranchStatus("*", 0);
      tree->SetBranchStatus("weight", 1);
      tree->SetBranchStatus("box", 1);
      tree->SetBranchStatus("nBTaggedJets",1);
      tree->SetBranchStatus("nSelectedJets", 1);
      tree->SetBranchStatus("nJets80", 1);
      tree->SetBranchStatus("nVetoMuons", 1);
      tree->SetBranchStatus("nVetoElectrons", 1);
      tree->SetBranchStatus("nLooseTaus", 1);
      tree->SetBranchStatus("leadingGenLeptonPt", 1);
      tree->SetBranchStatus("leadingGenLeptonEta", 1);
      tree->SetBranchStatus("leadingGenLeptonType", 1);
      tree->SetBranchStatus("MR", 1);
      tree->SetBranchStatus("Rsq", 1);
      tree->SetBranchStatus("dPhiRazor", 1);
      tree->SetBranchStatus("dPhiRazor_JESUp", 1);
      tree->SetBranchStatus("dPhiRazor_JESDown", 1);
      tree->SetBranchStatus("dPhiRazor_EESUp", 1);
      tree->SetBranchStatus("dPhiRazor_EESDown", 1);
      tree->SetBranchStatus("dPhiRazor_MESUp", 1);
      tree->SetBranchStatus("dPhiRazor_MESDown", 1);
      tree->SetBranchStatus("HLTDecision", 1);
      tree->SetBranchStatus("sf_facScaleUp", 1);
      tree->SetBranchStatus("sf_facScaleDown", 1);
      tree->SetBranchStatus("sf_renScaleUp", 1);
      tree->SetBranchStatus("sf_renScaleDown", 1);
      
      tree->SetBranchAddress("weight",&weight);
      tree->SetBranchAddress("box",&box);
      tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
      tree->SetBranchAddress("nSelectedJets",&nSelectedJets);
      tree->SetBranchAddress("nJets80",&nJets80);
      tree->SetBranchAddress("nVetoMuons",&nVetoMuons);
      tree->SetBranchAddress("nVetoElectrons",&nVetoElectrons);
      tree->SetBranchAddress("nLooseTaus",&nLooseTaus);
      tree->SetBranchAddress("leadingGenLeptonPt",&leadingGenLeptonPt);
      tree->SetBranchAddress("leadingGenLeptonEta",&leadingGenLeptonEta);
      tree->SetBranchAddress("leadingGenLeptonType",&leadingGenLeptonType);
      tree->SetBranchAddress("MR",&MR);
      tree->SetBranchAddress("Rsq",&Rsq);
      tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
      tree->SetBranchAddress("dPhiRazor_JESUp",&dPhiRazor_JESUp);
      tree->SetBranchAddress("dPhiRazor_JESDown",&dPhiRazor_JESDown);
      tree->SetBranchAddress("dPhiRazor_EESUp",&dPhiRazor_EESUp);
      tree->SetBranchAddress("dPhiRazor_EESDown",&dPhiRazor_EESDown);
      tree->SetBranchAddress("dPhiRazor_MESUp",&dPhiRazor_MESUp);
      tree->SetBranchAddress("dPhiRazor_MESDown",&dPhiRazor_MESDown);
      tree->SetBranchAddress("HLTDecision",&HLTDecision);
      tree->SetBranchAddress("sf_facScaleUp",&sf_facScaleUp);
      tree->SetBranchAddress("sf_facScaleDown",&sf_facScaleDown);
      tree->SetBranchAddress("sf_renScaleUp",&sf_renScaleUp);
      tree->SetBranchAddress("sf_renScaleDown",&sf_renScaleDown);

      bool isData = false;
      if ( processLabels[i] == "Data") isData = true;

      cout << "process: " << processLabels[i] << " | file " << inputfiles[i][j] << " | Total Entries: " << tree->GetEntries() << "\n";
      for(UInt_t ientry=0; ientry < tree->GetEntries(); ientry++) {
	tree->GetEntry(ientry);
      
	if (ientry % 100000 == 0) cout << "Event " << ientry << endl;      

	//******************************
	//Trigger Selection
	//******************************
	bool passTrigger = false;

	//Razor Hadronic Triggers
	if ( HLTDecision[134] || HLTDecision[135] || HLTDecision[136] 
	     || HLTDecision[137] || HLTDecision[138]
	     || HLTDecision[139] || HLTDecision[140]
	     || HLTDecision[141] || HLTDecision[142]
	     || HLTDecision[143] || HLTDecision[144]
	     ) passTrigger = true;
	if (!passTrigger) continue;
	
	//******************************
	//Selection Cuts 
	//******************************
	if (!( box == 11 || box == 12 )) continue;
	if (!(nJets80 >= 2)) continue;
	if (!(MR > 400 && Rsq > 0.25)) continue;
	if (!TMath::Finite(weight)) continue; 
	if (!( abs(leadingGenLeptonType) == 15)) continue;

	//******************************
	//Fill histograms
	//******************************
	if (isData) {
	} else {
	  //cout << leadingGenLeptonPt << " " << weight << "\n";
	  histLep1Pt->Fill(leadingGenLeptonPt, weight*lumi);
	  histLep1Pt_FScaleUp->Fill(leadingGenLeptonPt, sf_facScaleUp*weight*lumi);
	  histLep1Pt_FScaleDown->Fill(leadingGenLeptonPt, sf_facScaleDown*weight*lumi);
	  histLep1Pt_RScaleUp->Fill(leadingGenLeptonPt, sf_renScaleUp*weight*lumi);
	  histLep1Pt_RScaleDown->Fill(leadingGenLeptonPt, sf_renScaleDown*weight*lumi);
	  histLep1Eta->Fill(fabs(leadingGenLeptonEta), weight*lumi);
	  histLep1Eta_FScaleUp->Fill(fabs(leadingGenLeptonEta), sf_facScaleUp*weight*lumi);
	  histLep1Eta_FScaleDown->Fill(fabs(leadingGenLeptonEta), sf_facScaleDown*weight*lumi);
	  histLep1Eta_RScaleUp->Fill(fabs(leadingGenLeptonEta), sf_renScaleUp*weight*lumi);
	  histLep1Eta_RScaleDown->Fill(fabs(leadingGenLeptonEta), sf_renScaleDown*weight*lumi);
	  if (fabs(dPhiRazor) < 2.8) {
	    histLep1PtPassDPhiCut->Fill(leadingGenLeptonPt, weight*lumi);
	    histLep1PtPassDPhiCut_FScaleUp->Fill(leadingGenLeptonPt, sf_facScaleUp*weight*lumi);
	    histLep1PtPassDPhiCut_FScaleDown->Fill(leadingGenLeptonPt, sf_facScaleDown*weight*lumi);
	    histLep1PtPassDPhiCut_RScaleUp->Fill(leadingGenLeptonPt, sf_renScaleUp*weight*lumi);
	    histLep1PtPassDPhiCut_RScaleDown->Fill(leadingGenLeptonPt, sf_renScaleDown*weight*lumi);
	    histLep1EtaPassDPhiCut->Fill(fabs(leadingGenLeptonEta), weight*lumi);
	    histLep1EtaPassDPhiCut_FScaleUp->Fill(fabs(leadingGenLeptonEta), sf_facScaleUp*weight*lumi);
	    histLep1EtaPassDPhiCut_FScaleDown->Fill(fabs(leadingGenLeptonEta), sf_facScaleDown*weight*lumi);
	    histLep1EtaPassDPhiCut_RScaleUp->Fill(fabs(leadingGenLeptonEta), sf_renScaleUp*weight*lumi);
	    histLep1EtaPassDPhiCut_RScaleDown->Fill(fabs(leadingGenLeptonEta), sf_renScaleDown*weight*lumi);
	  }
	  if (fabs(dPhiRazor_JESUp) < 2.8) {
	    histLep1PtPassDPhiCut_JESUp->Fill(leadingGenLeptonPt, weight*lumi);
	    histLep1EtaPassDPhiCut_JESUp->Fill(fabs(leadingGenLeptonEta), weight*lumi);
	  }
	  if (fabs(dPhiRazor_JESDown) < 2.8) {
	    histLep1PtPassDPhiCut_JESDown->Fill(leadingGenLeptonPt, weight*lumi);
	    histLep1EtaPassDPhiCut_JESDown->Fill(fabs(leadingGenLeptonEta), weight*lumi);
	  }
	  if (abs(leadingGenLeptonType) == 11) {
	    if (fabs(dPhiRazor_EESUp) < 2.8) {
	      histLep1PtPassDPhiCut_LESUp->Fill(leadingGenLeptonPt, weight*lumi);
	      histLep1EtaPassDPhiCut_LESUp->Fill(fabs(leadingGenLeptonEta), weight*lumi);
	    }
	    if (fabs(dPhiRazor_EESDown) < 2.8) {
	      histLep1PtPassDPhiCut_LESDown->Fill(leadingGenLeptonPt, weight*lumi);
	      histLep1EtaPassDPhiCut_LESDown->Fill(fabs(leadingGenLeptonEta), weight*lumi);
	    }
	  }
	  if (abs(leadingGenLeptonType) == 13) {
	    if ( fabs(dPhiRazor_MESUp) < 2.8) {
	      histLep1PtPassDPhiCut_LESUp->Fill(leadingGenLeptonPt, weight*lumi);
	      histLep1EtaPassDPhiCut_LESUp->Fill(fabs(leadingGenLeptonEta), weight*lumi);
	    }
	    if (fabs(dPhiRazor_MESDown) < 2.8) {
	      histLep1PtPassDPhiCut_LESDown->Fill(leadingGenLeptonPt, weight*lumi);
	      histLep1EtaPassDPhiCut_LESDown->Fill(fabs(leadingGenLeptonEta), weight*lumi);
	    }
	  }
	  
	}
      }

      f->Close();

    }
  }
 



  //--------------------------------------------------------------------------------------------------------------
  // Make Plots
  //==============================================================================================================


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open("DPhiCutEfficiencyForLostTau.root", "RECREATE");
  file->cd();
 
  file->WriteTObject(histLep1Pt, "histLep1Pt", "WriteDelete");
  file->WriteTObject(histLep1Pt_FScaleUp, "histLep1Pt_FScaleUp","WriteDelete");
  file->WriteTObject(histLep1Pt_FScaleDown, "histLep1Pt_FScaleDown","WriteDelete");
  file->WriteTObject(histLep1Pt_RScaleUp, "histLep1Pt_RScaleUp","WriteDelete");
  file->WriteTObject(histLep1Pt_RScaleDown, "histLep1Pt_RScaleDown","WriteDelete");

  file->WriteTObject(histLep1PtPassDPhiCut, "histLep1PtPassDPhiCut", "WriteDelete");
  file->WriteTObject(histLep1PtPassDPhiCut_JESUp, "histLep1PtPassDPhiCut_JESUp", "WriteDelete");
  file->WriteTObject(histLep1PtPassDPhiCut_JESDown, "histLep1PtPassDPhiCut_JESDOwn", "WriteDelete");
  file->WriteTObject(histLep1PtPassDPhiCut_LESUp, "histLep1PtPassDPhiCut_LESUp", "WriteDelete");
  file->WriteTObject(histLep1PtPassDPhiCut_LESDown, "histLep1PtPassDPhiCut_LESDOwn", "WriteDelete");
  file->WriteTObject(histLep1PtPassDPhiCut_FScaleUp, "histLep1PtPassDPhiCut_FScaleUp", "WriteDelete");
  file->WriteTObject(histLep1PtPassDPhiCut_FScaleDown, "histLep1PtPassDPhiCut_FScaleDown", "WriteDelete");
  file->WriteTObject(histLep1PtPassDPhiCut_RScaleUp, "histLep1PtPassDPhiCut_RScaleUp", "WriteDelete");
  file->WriteTObject(histLep1PtPassDPhiCut_RScaleDown, "histLep1PtPassDPhiCut_RScaleDown", "WriteDelete");

 file->WriteTObject(histLep1Eta, "histLep1Eta", "WriteDelete");
  file->WriteTObject(histLep1Eta_FScaleUp, "histLep1Eta_FScaleUp","WriteDelete");
  file->WriteTObject(histLep1Eta_FScaleDown, "histLep1Eta_FScaleDown","WriteDelete");
  file->WriteTObject(histLep1Eta_RScaleUp, "histLep1Eta_RScaleUp","WriteDelete");
  file->WriteTObject(histLep1Eta_RScaleDown, "histLep1Eta_RScaleDown","WriteDelete");

  file->WriteTObject(histLep1EtaPassDPhiCut, "histLep1EtaPassDPhiCut", "WriteDelete");
  file->WriteTObject(histLep1EtaPassDPhiCut_JESUp, "histLep1EtaPassDPhiCut_JESUp", "WriteDelete");
  file->WriteTObject(histLep1EtaPassDPhiCut_JESDown, "histLep1EtaPassDPhiCut_JESDOwn", "WriteDelete");
  file->WriteTObject(histLep1EtaPassDPhiCut_LESUp, "histLep1EtaPassDPhiCut_LESUp", "WriteDelete");
  file->WriteTObject(histLep1EtaPassDPhiCut_LESDown, "histLep1EtaPassDPhiCut_LESDOwn", "WriteDelete");
  file->WriteTObject(histLep1EtaPassDPhiCut_FScaleUp, "histLep1EtaPassDPhiCut_FScaleUp", "WriteDelete");
  file->WriteTObject(histLep1EtaPassDPhiCut_FScaleDown, "histLep1EtaPassDPhiCut_FScaleDown", "WriteDelete");
  file->WriteTObject(histLep1EtaPassDPhiCut_RScaleUp, "histLep1EtaPassDPhiCut_RScaleUp", "WriteDelete");
  file->WriteTObject(histLep1EtaPassDPhiCut_RScaleDown, "histLep1EtaPassDPhiCut_RScaleDown", "WriteDelete");

  file->Close();
  delete file;       
    

}




void ComputeMTCutEfficiencyVsPt() {

  TFile *f = new TFile("VetoTauMTCutEfficiency.root","UPDATE");
  TH1D *num = (TH1D*)f->Get("histLep1PtPassMTCut");
  TH1D *den = (TH1D*)f->Get("histLep1Pt");
  TH1D *num_JESUp = (TH1D*)f->Get("histLep1PtPassMTCut_JESUp");
  TH1D *num_JESDown = (TH1D*)f->Get("histLep1PtPassMTCut_JESDOwn");
  TH1D *num_FScaleUp = (TH1D*)f->Get("histLep1PtPassMTCut_FScaleUp");
  TH1D *num_FScaleDown = (TH1D*)f->Get("histLep1PtPassMTCut_FScaleDown");
  TH1D *num_RScaleUp = (TH1D*)f->Get("histLep1PtPassMTCut_RScaleUp");
  TH1D *num_RScaleDown = (TH1D*)f->Get("histLep1PtPassMTCut_RScaleDown");
  TH1D *den_FScaleUp = (TH1D*)f->Get("histLep1Pt_FScaleUp");
  TH1D *den_FScaleDown = (TH1D*)f->Get("histLep1Pt_FScaleDown");
  TH1D *den_RScaleUp = (TH1D*)f->Get("histLep1Pt_RScaleUp");
  TH1D *den_RScaleDown = (TH1D*)f->Get("histLep1Pt_RScaleDown");

  TH1D *eff = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsPt");
  TH1D *eff_JESUp = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsPt_JESUp");
  TH1D *eff_JESDown = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsPt_JESDown");
  TH1D *eff_FScaleUp = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsPt_FScaleUp");
  TH1D *eff_FScaleDown = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsPt_FScaleDown");
  TH1D *eff_RScaleUp = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsPt_RScaleUp");
  TH1D *eff_RScaleDown = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsPt_RScaleDown");

  for (int i=1; i<eff->GetXaxis()->GetNbins()+1; i++) {
    double e = num->GetBinContent(i) / den->GetBinContent(i);
    double e_JESUp = num_JESUp->GetBinContent(i) / den->GetBinContent(i);
    double e_JESDown = num_JESDown->GetBinContent(i) / den->GetBinContent(i);
     double e_FScaleUp = num_FScaleUp->GetBinContent(i) / den_FScaleUp->GetBinContent(i);
    double e_FScaleDown = num_FScaleDown->GetBinContent(i) / den_FScaleDown->GetBinContent(i);
    double e_RScaleUp = num_RScaleUp->GetBinContent(i) / den_RScaleUp->GetBinContent(i);
    double e_RScaleDown = num_RScaleDown->GetBinContent(i) / den_RScaleDown->GetBinContent(i);

    double unc_JES = (e_JESUp - e_JESDown) / 0.5*(e_JESUp + e_JESDown);
    double unc_FScale = (e_FScaleUp - e_FScaleDown) / 0.5*(e_FScaleUp + e_FScaleDown);
    double unc_RScale = (e_RScaleUp - e_RScaleDown) / 0.5*(e_RScaleUp + e_RScaleDown);
    double unc_total = sqrt( pow(unc_JES,2) + pow(unc_FScale,2) + pow(unc_RScale,2) );

    cout << "Bin " << i << " : " << e << " : " << unc_JES << " "  << unc_FScale << " " << unc_RScale << " : " << unc_total << "\n";    

    eff->SetBinContent(i, e);
    eff->SetBinError(i, unc_total * e);
    
    eff_JESUp->SetBinContent(i, num_JESUp->GetBinContent(i) / den->GetBinContent(i));
    eff_JESDown->SetBinContent(i, num_JESDown->GetBinContent(i) / den->GetBinContent(i));
    eff_FScaleUp->SetBinContent(i, num_FScaleUp->GetBinContent(i) / den_FScaleUp->GetBinContent(i));
    eff_FScaleDown->SetBinContent(i, num_FScaleDown->GetBinContent(i) / den_FScaleDown->GetBinContent(i));
    eff_RScaleUp->SetBinContent(i, num_RScaleUp->GetBinContent(i) / den_RScaleUp->GetBinContent(i));
    eff_RScaleDown->SetBinContent(i, num_RScaleDown->GetBinContent(i) / den_RScaleDown->GetBinContent(i));
  }

  f->WriteTObject(eff,"VetoTauMTCutEfficiencyVsPt","WriteDelete");
  f->WriteTObject(eff_JESUp,"VetoTauMTCutEfficiencyVsPt_JESUp","WriteDelete");
  f->WriteTObject(eff_JESDown,"VetoTauMTCutEfficiencyVsPt_JESDown","WriteDelete");
  f->WriteTObject(eff_FScaleUp,"VetoTauMTCutEfficiencyVsPt_FScaleUp","WriteDelete");
  f->WriteTObject(eff_FScaleDown,"VetoTauMTCutEfficiencyVsPt_FScaleDown","WriteDelete");
  f->WriteTObject(eff_RScaleUp,"VetoTauMTCutEfficiencyVsPt_RScaleUp","WriteDelete");
  f->WriteTObject(eff_RScaleDown,"VetoTauMTCutEfficiencyVsPt_RScaleDown","WriteDelete");
  f->Close();
  
}

void ComputeMTCutEfficiencyVsEta() {

  TFile *f = new TFile("VetoTauMTCutEfficiency.root","UPDATE");
  TH1D *num = (TH1D*)f->Get("histLep1EtaPassMTCut");
  TH1D *den = (TH1D*)f->Get("histLep1Eta");
  TH1D *num_JESUp = (TH1D*)f->Get("histLep1EtaPassMTCut_JESUp");
  TH1D *num_JESDown = (TH1D*)f->Get("histLep1EtaPassMTCut_JESDOwn");
  TH1D *num_FScaleUp = (TH1D*)f->Get("histLep1EtaPassMTCut_FScaleUp");
  TH1D *num_FScaleDown = (TH1D*)f->Get("histLep1EtaPassMTCut_FScaleDown");
  TH1D *num_RScaleUp = (TH1D*)f->Get("histLep1EtaPassMTCut_RScaleUp");
  TH1D *num_RScaleDown = (TH1D*)f->Get("histLep1EtaPassMTCut_RScaleDown");
  TH1D *den_FScaleUp = (TH1D*)f->Get("histLep1Eta_FScaleUp");
  TH1D *den_FScaleDown = (TH1D*)f->Get("histLep1Eta_FScaleDown");
  TH1D *den_RScaleUp = (TH1D*)f->Get("histLep1Eta_RScaleUp");
  TH1D *den_RScaleDown = (TH1D*)f->Get("histLep1Eta_RScaleDown");

  TH1D *eff = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsEta");
  TH1D *eff_JESUp = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsEta_JESUp");
  TH1D *eff_JESDown = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsEta_JESDown");
  TH1D *eff_FScaleUp = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsEta_FScaleUp");
  TH1D *eff_FScaleDown = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsEta_FScaleDown");
  TH1D *eff_RScaleUp = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsEta_RScaleUp");
  TH1D *eff_RScaleDown = (TH1D*)num->Clone("VetoTauMTCutEfficiencyVsEta_RScaleDown");

  for (int i=1; i<eff->GetXaxis()->GetNbins()+1; i++) {
    double e = num->GetBinContent(i) / den->GetBinContent(i);
    double e_JESUp = num_JESUp->GetBinContent(i) / den->GetBinContent(i);
    double e_JESDown = num_JESDown->GetBinContent(i) / den->GetBinContent(i);
     double e_FScaleUp = num_FScaleUp->GetBinContent(i) / den_FScaleUp->GetBinContent(i);
    double e_FScaleDown = num_FScaleDown->GetBinContent(i) / den_FScaleDown->GetBinContent(i);
    double e_RScaleUp = num_RScaleUp->GetBinContent(i) / den_RScaleUp->GetBinContent(i);
    double e_RScaleDown = num_RScaleDown->GetBinContent(i) / den_RScaleDown->GetBinContent(i);

    double unc_JES = (e_JESUp - e_JESDown) / 0.5*(e_JESUp + e_JESDown);
    double unc_FScale = (e_FScaleUp - e_FScaleDown) / 0.5*(e_FScaleUp + e_FScaleDown);
    double unc_RScale = (e_RScaleUp - e_RScaleDown) / 0.5*(e_RScaleUp + e_RScaleDown);
    double unc_total = sqrt( pow(unc_JES,2) + pow(unc_FScale,2) + pow(unc_RScale,2) );

    cout << "Bin " << i << " : " << e << " : " << unc_JES << " "  << unc_FScale << " " << unc_RScale << " : " << unc_total << "\n";    

    eff->SetBinContent(i, e);
    eff->SetBinError(i, unc_total * e);
    
    eff_JESUp->SetBinContent(i, num_JESUp->GetBinContent(i) / den->GetBinContent(i));
    eff_JESDown->SetBinContent(i, num_JESDown->GetBinContent(i) / den->GetBinContent(i));
    eff_FScaleUp->SetBinContent(i, num_FScaleUp->GetBinContent(i) / den_FScaleUp->GetBinContent(i));
    eff_FScaleDown->SetBinContent(i, num_FScaleDown->GetBinContent(i) / den_FScaleDown->GetBinContent(i));
    eff_RScaleUp->SetBinContent(i, num_RScaleUp->GetBinContent(i) / den_RScaleUp->GetBinContent(i));
    eff_RScaleDown->SetBinContent(i, num_RScaleDown->GetBinContent(i) / den_RScaleDown->GetBinContent(i));
  }

  f->WriteTObject(eff,"VetoTauMTCutEfficiencyVsEta","WriteDelete");
  f->WriteTObject(eff_JESUp,"VetoTauMTCutEfficiencyVsEta_JESUp","WriteDelete");
  f->WriteTObject(eff_JESDown,"VetoTauMTCutEfficiencyVsEta_JESDown","WriteDelete");
  f->WriteTObject(eff_FScaleUp,"VetoTauMTCutEfficiencyVsEta_FScaleUp","WriteDelete");
  f->WriteTObject(eff_FScaleDown,"VetoTauMTCutEfficiencyVsEta_FScaleDown","WriteDelete");
  f->WriteTObject(eff_RScaleUp,"VetoTauMTCutEfficiencyVsEta_RScaleUp","WriteDelete");
  f->WriteTObject(eff_RScaleDown,"VetoTauMTCutEfficiencyVsEta_RScaleDown","WriteDelete");
  f->Close();
  
}




void ComputeDPhiCutEfficiencyVsPt() {

  TFile *f = new TFile("DPhiCutEfficiencyForLostTau.root","UPDATE");
  TH1D *num = (TH1D*)f->Get("histLep1PtPassDPhiCut");
  TH1D *den = (TH1D*)f->Get("histLep1Pt");
  TH1D *num_JESUp = (TH1D*)f->Get("histLep1PtPassDPhiCut_JESUp");
  TH1D *num_JESDown = (TH1D*)f->Get("histLep1PtPassDPhiCut_JESDOwn");
  TH1D *num_LESUp = (TH1D*)f->Get("histLep1PtPassDPhiCut_LESUp");
  TH1D *num_LESDown = (TH1D*)f->Get("histLep1PtPassDPhiCut_LESDOwn");
  TH1D *num_FScaleUp = (TH1D*)f->Get("histLep1PtPassDPhiCut_FScaleUp");
  TH1D *num_FScaleDown = (TH1D*)f->Get("histLep1PtPassDPhiCut_FScaleDown");
  TH1D *num_RScaleUp = (TH1D*)f->Get("histLep1PtPassDPhiCut_RScaleUp");
  TH1D *num_RScaleDown = (TH1D*)f->Get("histLep1PtPassDPhiCut_RScaleDown");
  TH1D *den_FScaleUp = (TH1D*)f->Get("histLep1Pt_FScaleUp");
  TH1D *den_FScaleDown = (TH1D*)f->Get("histLep1Pt_FScaleDown");
  TH1D *den_RScaleUp = (TH1D*)f->Get("histLep1Pt_RScaleUp");
  TH1D *den_RScaleDown = (TH1D*)f->Get("histLep1Pt_RScaleDown");

  TH1D *eff = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsPt");
  TH1D *eff_JESUp = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsPt_JESUp");
  TH1D *eff_JESDown = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsPt_JESDown");
  TH1D *eff_LESUp = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsPt_LESUp");
  TH1D *eff_LESDown = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsPt_LESDown");
  TH1D *eff_FScaleUp = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsPt_FScaleUp");
  TH1D *eff_FScaleDown = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsPt_FScaleDown");
  TH1D *eff_RScaleUp = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsPt_RScaleUp");
  TH1D *eff_RScaleDown = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsPt_RScaleDown");

  for (int i=1; i<eff->GetXaxis()->GetNbins()+1; i++) {
    double e = num->GetBinContent(i) / den->GetBinContent(i);
    double e_JESUp = num_JESUp->GetBinContent(i) / den->GetBinContent(i);
    double e_JESDown = num_JESDown->GetBinContent(i) / den->GetBinContent(i);
    double e_LESUp = num_LESUp->GetBinContent(i) / den->GetBinContent(i);
    double e_LESDown = num_LESDown->GetBinContent(i) / den->GetBinContent(i);
    double e_FScaleUp = num_FScaleUp->GetBinContent(i) / den_FScaleUp->GetBinContent(i);
    double e_FScaleDown = num_FScaleDown->GetBinContent(i) / den_FScaleDown->GetBinContent(i);
    double e_RScaleUp = num_RScaleUp->GetBinContent(i) / den_RScaleUp->GetBinContent(i);
    double e_RScaleDown = num_RScaleDown->GetBinContent(i) / den_RScaleDown->GetBinContent(i);

    double unc_JES = (e_JESUp - e_JESDown) / 0.5*(e_JESUp + e_JESDown);
    double unc_LES = (e_LESUp - e_LESDown) / 0.5*(e_LESUp + e_LESDown);
    double unc_FScale = (e_FScaleUp - e_FScaleDown) / 0.5*(e_FScaleUp + e_FScaleDown);
    double unc_RScale = (e_RScaleUp - e_RScaleDown) / 0.5*(e_RScaleUp + e_RScaleDown);
    double unc_total = sqrt( pow(unc_JES,2) + pow(unc_LES,2) + pow(unc_FScale,2) + pow(unc_RScale,2) );

    cout << "Bin " << i << " : " << e << " : " << unc_JES << " " << unc_LES << " " << unc_FScale << " " << unc_RScale << " : " << unc_total << "\n";    

    eff->SetBinContent(i, e);
    eff->SetBinError(i, unc_total * e);
    
    eff_JESUp->SetBinContent(i, num_JESUp->GetBinContent(i) / den->GetBinContent(i));
    eff_JESDown->SetBinContent(i, num_JESDown->GetBinContent(i) / den->GetBinContent(i));
    eff_LESUp->SetBinContent(i, num_LESUp->GetBinContent(i) / den->GetBinContent(i));
    eff_LESDown->SetBinContent(i, num_LESDown->GetBinContent(i) / den->GetBinContent(i));
    eff_FScaleUp->SetBinContent(i, num_FScaleUp->GetBinContent(i) / den_FScaleUp->GetBinContent(i));
    eff_FScaleDown->SetBinContent(i, num_FScaleDown->GetBinContent(i) / den_FScaleDown->GetBinContent(i));
    eff_RScaleUp->SetBinContent(i, num_RScaleUp->GetBinContent(i) / den_RScaleUp->GetBinContent(i));
    eff_RScaleDown->SetBinContent(i, num_RScaleDown->GetBinContent(i) / den_RScaleDown->GetBinContent(i));
  }

  f->WriteTObject(eff,"VetoTauDPhiCutEfficiencyVsPt","WriteDelete");
  f->WriteTObject(eff_JESUp,"VetoTauDPhiCutEfficiencyVsPt_JESUp","WriteDelete");
  f->WriteTObject(eff_JESDown,"VetoTauDPhiCutEfficiencyVsPt_JESDown","WriteDelete");
  f->WriteTObject(eff_LESUp,"VetoTauDPhiCutEfficiencyVsPt_LESUp","WriteDelete");
  f->WriteTObject(eff_LESDown,"VetoTauDPhiCutEfficiencyVsPt_LESDown","WriteDelete");
  f->WriteTObject(eff_FScaleUp,"VetoTauDPhiCutEfficiencyVsPt_FScaleUp","WriteDelete");
  f->WriteTObject(eff_FScaleDown,"VetoTauDPhiCutEfficiencyVsPt_FScaleDown","WriteDelete");
  f->WriteTObject(eff_RScaleUp,"VetoTauDPhiCutEfficiencyVsPt_RScaleUp","WriteDelete");
  f->WriteTObject(eff_RScaleDown,"VetoTauDPhiCutEfficiencyVsPt_RScaleDown","WriteDelete");
  f->Close();
  
}



void ComputeDPhiCutEfficiencyVsEta() {

  TFile *f = new TFile("DPhiCutEfficiencyForLostTau.root","UPDATE");
  TH1D *num = (TH1D*)f->Get("histLep1EtaPassDPhiCut");
  TH1D *den = (TH1D*)f->Get("histLep1Eta");
  TH1D *num_JESUp = (TH1D*)f->Get("histLep1EtaPassDPhiCut_JESUp");
  TH1D *num_JESDown = (TH1D*)f->Get("histLep1EtaPassDPhiCut_JESDOwn");
  TH1D *num_LESUp = (TH1D*)f->Get("histLep1EtaPassDPhiCut_LESUp");
  TH1D *num_LESDown = (TH1D*)f->Get("histLep1EtaPassDPhiCut_LESDOwn");
  TH1D *num_FScaleUp = (TH1D*)f->Get("histLep1EtaPassDPhiCut_FScaleUp");
  TH1D *num_FScaleDown = (TH1D*)f->Get("histLep1EtaPassDPhiCut_FScaleDown");
  TH1D *num_RScaleUp = (TH1D*)f->Get("histLep1EtaPassDPhiCut_RScaleUp");
  TH1D *num_RScaleDown = (TH1D*)f->Get("histLep1EtaPassDPhiCut_RScaleDown");
  TH1D *den_FScaleUp = (TH1D*)f->Get("histLep1Eta_FScaleUp");
  TH1D *den_FScaleDown = (TH1D*)f->Get("histLep1Eta_FScaleDown");
  TH1D *den_RScaleUp = (TH1D*)f->Get("histLep1Eta_RScaleUp");
  TH1D *den_RScaleDown = (TH1D*)f->Get("histLep1Eta_RScaleDown");

  TH1D *eff = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsEta");
  TH1D *eff_JESUp = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsEta_JESUp");
  TH1D *eff_JESDown = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsEta_JESDown");
  TH1D *eff_LESUp = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsEta_LESUp");
  TH1D *eff_LESDown = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsEta_LESDown");
  TH1D *eff_FScaleUp = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsEta_FScaleUp");
  TH1D *eff_FScaleDown = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsEta_FScaleDown");
  TH1D *eff_RScaleUp = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsEta_RScaleUp");
  TH1D *eff_RScaleDown = (TH1D*)num->Clone("VetoTauDPhiCutEfficiencyVsEta_RScaleDown");

  for (int i=1; i<eff->GetXaxis()->GetNbins()+1; i++) {
    double e = num->GetBinContent(i) / den->GetBinContent(i);
    double e_JESUp = num_JESUp->GetBinContent(i) / den->GetBinContent(i);
    double e_JESDown = num_JESDown->GetBinContent(i) / den->GetBinContent(i);
    double e_LESUp = num_LESUp->GetBinContent(i) / den->GetBinContent(i);
    double e_LESDown = num_LESDown->GetBinContent(i) / den->GetBinContent(i);
    double e_FScaleUp = num_FScaleUp->GetBinContent(i) / den_FScaleUp->GetBinContent(i);
    double e_FScaleDown = num_FScaleDown->GetBinContent(i) / den_FScaleDown->GetBinContent(i);
    double e_RScaleUp = num_RScaleUp->GetBinContent(i) / den_RScaleUp->GetBinContent(i);
    double e_RScaleDown = num_RScaleDown->GetBinContent(i) / den_RScaleDown->GetBinContent(i);

    double unc_JES = (e_JESUp - e_JESDown) / 0.5*(e_JESUp + e_JESDown);
    double unc_LES = (e_LESUp - e_LESDown) / 0.5*(e_LESUp + e_LESDown);
    double unc_FScale = (e_FScaleUp - e_FScaleDown) / 0.5*(e_FScaleUp + e_FScaleDown);
    double unc_RScale = (e_RScaleUp - e_RScaleDown) / 0.5*(e_RScaleUp + e_RScaleDown);
    double unc_total = sqrt( pow(unc_JES,2) + pow(unc_LES,2) + pow(unc_FScale,2) + pow(unc_RScale,2) );

    cout << "Bin " << i << " : " << e << " : " << unc_JES << " " << unc_LES << " " << unc_FScale << " " << unc_RScale << " : " << unc_total << "\n";    

    eff->SetBinContent(i, e);
    eff->SetBinError(i, unc_total * e);
    
    eff_JESUp->SetBinContent(i, num_JESUp->GetBinContent(i) / den->GetBinContent(i));
    eff_JESDown->SetBinContent(i, num_JESDown->GetBinContent(i) / den->GetBinContent(i));
    eff_LESUp->SetBinContent(i, num_LESUp->GetBinContent(i) / den->GetBinContent(i));
    eff_LESDown->SetBinContent(i, num_LESDown->GetBinContent(i) / den->GetBinContent(i));
    eff_FScaleUp->SetBinContent(i, num_FScaleUp->GetBinContent(i) / den_FScaleUp->GetBinContent(i));
    eff_FScaleDown->SetBinContent(i, num_FScaleDown->GetBinContent(i) / den_FScaleDown->GetBinContent(i));
    eff_RScaleUp->SetBinContent(i, num_RScaleUp->GetBinContent(i) / den_RScaleUp->GetBinContent(i));
    eff_RScaleDown->SetBinContent(i, num_RScaleDown->GetBinContent(i) / den_RScaleDown->GetBinContent(i));
  }

  f->WriteTObject(eff,"VetoTauDPhiCutEfficiencyVsEta","WriteDelete");
  f->WriteTObject(eff_JESUp,"VetoTauDPhiCutEfficiencyVsEta_JESUp","WriteDelete");
  f->WriteTObject(eff_JESDown,"VetoTauDPhiCutEfficiencyVsEta_JESDown","WriteDelete");
  f->WriteTObject(eff_LESUp,"VetoTauDPhiCutEfficiencyVsEta_LESUp","WriteDelete");
  f->WriteTObject(eff_LESDown,"VetoTauDPhiCutEfficiencyVsEta_LESDown","WriteDelete");
  f->WriteTObject(eff_FScaleUp,"VetoTauDPhiCutEfficiencyVsEta_FScaleUp","WriteDelete");
  f->WriteTObject(eff_FScaleDown,"VetoTauDPhiCutEfficiencyVsEta_FScaleDown","WriteDelete");
  f->WriteTObject(eff_RScaleUp,"VetoTauDPhiCutEfficiencyVsEta_RScaleUp","WriteDelete");
  f->WriteTObject(eff_RScaleDown,"VetoTauDPhiCutEfficiencyVsEta_RScaleDown","WriteDelete");
  f->Close();
  
}





void UnfoldMTCutEffForVetoTau( int option = -1) {

  vector<string> datafiles;
  vector<vector<string> > bkgfiles;
  vector<string> processLabels;
  vector<int> colors;

  vector<string> bkgfiles_ttbar;
  vector<string> bkgfiles_wjets;
  vector<string> bkgfiles_singletop;
  vector<string> bkgfiles_dy;  
  vector<string> bkgfiles_vv; 
  vector<string> bkgfiles_qcd;
  vector<string> bkgfiles_znunu;

  // bkgfiles_wjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108/FullRazorInclusive_TTJets_1pb_weighted.root");
  // bkgfiles_ttbar.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108/FullRazorInclusive_WJetsToLNu_HTBinned_1pb_weighted.root");

  bkgfiles_wjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108_NEW/FullRazorInclusive_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
  bkgfiles_wjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108_NEW/FullRazorInclusive_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
  bkgfiles_wjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108_NEW/FullRazorInclusive_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
  bkgfiles_wjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108_NEW/FullRazorInclusive_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
  bkgfiles_wjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108_NEW/FullRazorInclusive_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
  bkgfiles_wjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108_NEW/FullRazorInclusive_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
  bkgfiles_wjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108_NEW/FullRazorInclusive_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
  bkgfiles_ttbar.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108_NEW/FullRazorInclusive_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
  bkgfiles_ttbar.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108_NEW/FullRazorInclusive_TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
  bkgfiles_ttbar.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108_NEW/FullRazorInclusive_TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");


  bkgfiles.push_back(bkgfiles_wjets);
  bkgfiles.push_back(bkgfiles_ttbar);

  processLabels.push_back("WJets");  
  processLabels.push_back("TTJets");  
  
  colors.push_back(kRed);
  colors.push_back(kGreen+2);  
  
  double lumi = 2185;


  //*********************************************************************
  //Run
  //*********************************************************************
  //RunUnfoldMTCutEffForVetoTaus(datafiles, bkgfiles,processLabels,  colors, lumi);
   ComputeMTCutEfficiencyVsPt();
   ComputeMTCutEfficiencyVsEta();

  //RunUnfoldDPhiCutEffForVetoTaus(datafiles, bkgfiles,processLabels,  colors, lumi);
   ComputeDPhiCutEfficiencyVsPt();
   ComputeDPhiCutEfficiencyVsEta();

}

