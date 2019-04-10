

//================================================================================================
//
// Simple Example
//
//root -l RazorAnalyzer/macros/ObjectStudies/MakeElectronEfficiencyPlots.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/ElectronNtuple/ElectronNtuple_PromptGenLevel_TTJets_25ns.root",-1,"Electron")'
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"
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
#include <TH1F.h>                
#include <TH2F.h>                
#include <TCanvas.h>                
#include <TLegend.h> 
#include <TPaveText.h>
#include <THStack.h> 
#include <TFractionFitter.h> 

#include "RazorAnalyzer/include/ControlSampleEvents.h"
#include "RazorAnalyzer/macros/tdrstyle.C"
#include "RazorAnalyzer/macros/CMS_lumi.C"

#endif

using namespace RooFit ;

void RunPhotonControlSample_FitSigmaIEtaIEta(  vector<string> datafiles,  vector<string> fakeTemplateFiles,  vector<string> promptTemplateFiles,  double lumi, string option, int channelOption = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  float MR = 0;
  float Rsq = 0;
  float mll = 0;

  bool printdebug = false;

  TFile *PromptFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/PhotonTemplates_NoR9Triggers_Pt185MR300Rsq0p15.root", "READ");
  TH1F *promptHist = (TH1F*)PromptFile->Get("PhotonSigmaIEtaIEtaTemplate_Prompt_Barrel"); 
  TH1F *promptHist_EE = (TH1F*)PromptFile->Get("PhotonSigmaIEtaIEtaTemplate_Prompt_Endcap");
  TH1F *fakeHist_EB = (TH1F*)PromptFile->Get("PhotonSigmaIEtaIEtaTemplate_Fake_Barrel"); 
  TH1F *fakeHist_EE = (TH1F*)PromptFile->Get("PhotonSigmaIEtaIEtaTemplate_Fake_Endcap");

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  double MRBins[] = {300, 400, 500, 700, 900, 1200, 4000};
  double RsqBins[] = {0.15, 0.2, 0.25, 0.30, 0.41, 0.52, 1.5};

  const int NMRBins = sizeof(MRBins)/sizeof(double)-1;
  const int NRsqBins = sizeof(RsqBins)/sizeof(double)-1;

  double sigmaietaieta_bins_EB[100] = {0.};
  double sigmaietaieta_bins_EE[100] = {0.};

  sigmaietaieta_bins_EB[21] = 0.0103;
  sigmaietaieta_bins_EE[55] = 0.0271;

  const int NBinsSigmaietaieta = sizeof(sigmaietaieta_bins_EB)/sizeof(double)-1;

  for(int a = 1; a < 100; a++)
    {
      if(a!=21)
	sigmaietaieta_bins_EB[a] = sigmaietaieta_bins_EB[a-1] + 0.0005;
      if(a!=55)
	sigmaietaieta_bins_EE[a] = sigmaietaieta_bins_EE[a-1] + 0.0005;
    }

  vector<vector<string> > inputfiles;
  vector<string> processLabels;
  vector<int> color;

  inputfiles.push_back(datafiles);
  processLabels.push_back("Data");
  color.push_back(kBlack);
  inputfiles.push_back(fakeTemplateFiles);
  processLabels.push_back("FakeTemplate");
  color.push_back(kBlack);
  inputfiles.push_back(promptTemplateFiles);
  processLabels.push_back("PromptTemplate");
  color.push_back(kBlack);
  

  vector<TH1D*> histSigmaIetaIeta_EB;
  vector<TH1D*> histSigmaIetaIeta_EE;
  vector<TH1D*> histSigmaIetaIetaTemplate_EB;
  vector<TH1D*> histSigmaIetaIetaTemplate_EE;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    histSigmaIetaIeta_EB.push_back(new TH1D(Form("histSigmaIetaIeta_EB_%s",processLabels[i].c_str()), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB));
    histSigmaIetaIeta_EE.push_back(new TH1D(Form("histSigmaIetaIeta_EE_%s",processLabels[i].c_str()), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EE));
    histSigmaIetaIeta_EB[i]->Sumw2();
    histSigmaIetaIeta_EE[i]->Sumw2();

    histSigmaIetaIetaTemplate_EB.push_back(new TH1D(Form("histSigmaIetaIetaTemplate_EB_%s",processLabels[i].c_str()), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB));
    histSigmaIetaIetaTemplate_EE.push_back(new TH1D(Form("histSigmaIetaIetaTemplate_EE_%s",processLabels[i].c_str()), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EE));
    histSigmaIetaIetaTemplate_EB[i]->Sumw2();
    histSigmaIetaIetaTemplate_EE[i]->Sumw2();
  }

  TH1D *histSigmaIetaIetaFakeTemplate_EB_binned[NMRBins][NRsqBins];
  TH1D *histSigmaIetaIetaFakeTemplate_EE_binned[NMRBins][NRsqBins];

  TH1D *histSigmaIetaIetaFakeTemplate_EB_MRbinned[NMRBins];
  TH1D *histSigmaIetaIetaFakeTemplate_EE_MRbinned[NMRBins];

  TH1D *histSigmaIetaIetaFakeTemplate_EB_Rsqbinned[NRsqBins];
  TH1D *histSigmaIetaIetaFakeTemplate_EE_Rsqbinned[NRsqBins];

  for (int a=0;a<NMRBins;a++) {
    histSigmaIetaIetaFakeTemplate_EB_MRbinned[a]  = new TH1D (Form("histSigmaIetaIetaFakeTemplate_EB_MRbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB);
    histSigmaIetaIetaFakeTemplate_EE_MRbinned[a]  = new TH1D (Form("histSigmaIetaIetaFakeTemplate_EE_MRbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EE);
  }

  for (int a=0;a<NRsqBins;a++) {
    histSigmaIetaIetaFakeTemplate_EB_Rsqbinned[a] = new TH1D (Form("histSigmaIetaIetaFakeTemplate_EB_Rsqbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB);
    histSigmaIetaIetaFakeTemplate_EE_Rsqbinned[a] = new TH1D (Form("histSigmaIetaIetaFakeTemplate_EE_Rsqbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EE);
  }

  for (int a=0;a<NMRBins;a++) {
    for (int b=0;b<NRsqBins;b++){
      histSigmaIetaIetaFakeTemplate_EB_binned[a][b] = new TH1D (Form("histSigmaIetaIetaFakeTemplate_EB_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB);
      histSigmaIetaIetaFakeTemplate_EB_binned[a][b] -> Sumw2();

      histSigmaIetaIetaFakeTemplate_EE_binned[a][b] = new TH1D (Form("histSigmaIetaIetaFakeTemplate_EE_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB);
      histSigmaIetaIetaFakeTemplate_EE_binned[a][b] -> Sumw2();
    }
  }

  TH1D *histSigmaIetaIetaPromptTemplate_EB_binned[NMRBins][NRsqBins];
  TH1D *histSigmaIetaIetaPromptTemplate_EE_binned[NMRBins][NRsqBins];

  TH1D *histSigmaIetaIetaPromptTemplate_EB_MRbinned[NMRBins];
  TH1D *histSigmaIetaIetaPromptTemplate_EE_MRbinned[NMRBins];

  TH1D *histSigmaIetaIetaPromptTemplate_EB_Rsqbinned[NRsqBins];
  TH1D *histSigmaIetaIetaPromptTemplate_EE_Rsqbinned[NRsqBins];

  for (int a=0;a<NMRBins;a++) {
    histSigmaIetaIetaPromptTemplate_EB_MRbinned[a]  = new TH1D (Form("histSigmaIetaIetaPromptTemplate_EB_MRbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB);
    histSigmaIetaIetaPromptTemplate_EE_MRbinned[a]  = new TH1D (Form("histSigmaIetaIetaPromptTemplate_EE_MRbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EE);
  }

  for (int a=0;a<NRsqBins;a++) {
    histSigmaIetaIetaPromptTemplate_EB_Rsqbinned[a] = new TH1D (Form("histSigmaIetaIetaPromptTemplate_EB_Rsqbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB);
    histSigmaIetaIetaPromptTemplate_EE_Rsqbinned[a] = new TH1D (Form("histSigmaIetaIetaPromptTemplate_EE_Rsqbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EE);
  }

  for (int a=0;a<NMRBins;a++) {
    for (int b=0;b<NRsqBins;b++){
      histSigmaIetaIetaPromptTemplate_EB_binned[a][b] = new TH1D (Form("histSigmaIetaIetaPromptTemplate_EB_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB);
      histSigmaIetaIetaPromptTemplate_EB_binned[a][b] -> Sumw2();

      histSigmaIetaIetaPromptTemplate_EE_binned[a][b] = new TH1D (Form("histSigmaIetaIetaPromptTemplate_EE_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB);
      histSigmaIetaIetaPromptTemplate_EE_binned[a][b] -> Sumw2();
    }
  }



  TH1D *histSigmaIetaIeta_EB_binned[NMRBins][NRsqBins];
  TH1D *histSigmaIetaIeta_EE_binned[NMRBins][NRsqBins];

  TH1D *histSigmaIetaIeta_EB_MRbinned[NMRBins];
  TH1D *histSigmaIetaIeta_EE_MRbinned[NMRBins];

  TH1D *histSigmaIetaIeta_EB_Rsqbinned[NRsqBins];
  TH1D *histSigmaIetaIeta_EE_Rsqbinned[NRsqBins];

  for (int a=0;a<NMRBins;a++) {
    histSigmaIetaIeta_EB_MRbinned[a]  = new TH1D (Form("histSigmaIetaIeta_EB_MRbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB);
    histSigmaIetaIeta_EE_MRbinned[a]  = new TH1D (Form("histSigmaIetaIeta_EE_MRbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EE);
  }

  for (int a=0;a<NRsqBins;a++) {
    histSigmaIetaIeta_EB_Rsqbinned[a] = new TH1D (Form("histSigmaIetaIeta_EB_Rsqbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB);
    histSigmaIetaIeta_EE_Rsqbinned[a] = new TH1D (Form("histSigmaIetaIeta_EE_Rsqbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EE);
  }

  for (int a=0;a<NMRBins;a++) {
    for (int b=0;b<NRsqBins;b++){
      histSigmaIetaIeta_EB_binned[a][b] = new TH1D (Form("histSigmaIetaIeta_EB_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB);
      histSigmaIetaIeta_EB_binned[a][b] -> Sumw2();

      histSigmaIetaIeta_EE_binned[a][b] = new TH1D (Form("histSigmaIetaIeta_EE_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", NBinsSigmaietaieta, sigmaietaieta_bins_EB);
      histSigmaIetaIeta_EE_binned[a][b] -> Sumw2();
    }
  }

  TH2D *histSigmaIetaIeta_EB_MRRsq= new TH2D("histSigmaIetaIeta_EB_MRRsq", "; MR; Rsq", NMRBins, MRBins, NRsqBins, RsqBins);
  TH2D *histSigmaIetaIeta_EE_MRRsq= new TH2D("histSigmaIetaIeta_EE_MRRsq", "; MR; Rsq", NMRBins, MRBins, NRsqBins, RsqBins);
  TH1D *histSigmaIetaIeta_EB_MR= new TH1D("histSigmaIetaIeta_EB_MR", "; MR", NMRBins, MRBins);
  TH1D *histSigmaIetaIeta_EB_Rsq= new TH1D("histSigmaIetaIeta_EB_Rsq", "; Rsq", NRsqBins, RsqBins);
  TH1D *histSigmaIetaIeta_EE_MR= new TH1D("histSigmaIetaIeta_EE_MR", "; MR", NMRBins, MRBins);
  TH1D *histSigmaIetaIeta_EE_Rsq= new TH1D("histSigmaIetaIeta_EE_Rsq", "; Rsq", NRsqBins, RsqBins);

  histSigmaIetaIeta_EB_MR->Sumw2();
  histSigmaIetaIeta_EB_Rsq->Sumw2();
  histSigmaIetaIeta_EE_MR->Sumw2();
  histSigmaIetaIeta_EE_Rsq->Sumw2();


  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  for (uint i=0; i < inputfiles.size(); ++i) {

    //for duplicate event checking
    map<pair<uint,uint>, bool > processedRunEvents;

    for (uint j=0; j < inputfiles[i].size(); ++j) {
      ControlSampleEvents *events = new ControlSampleEvents;
      events->LoadTree(inputfiles[i][j].c_str(), ControlSampleEvents::kTreeType_Photon_Full);

      bool isData = false;
      if ( processLabels[i] == "Data") isData = true;
    
      cout << "process: " << processLabels[i] << " | file " << inputfiles[i][j] << " | Total Entries: " << events->tree_->GetEntries() << "\n";
      // for(UInt_t ientry=0; ientry < 10000; ientry++) {       	
      for(UInt_t ientry=0; ientry < events->tree_->GetEntries(); ientry++) {       	
	events->tree_->GetEntry(ientry);
      

	if (ientry % 100000 == 0) cout << "Event " << ientry << endl;      

	double puWeight = 1;      
	double weight = 1;
	if ( processLabels[i] == "PromptTemplate") {
	  weight = lumi * events->weight;
	}

	if (isnan(events->weight) || isinf(events->weight)) {
	  continue;
	  cout << "...bad event: " << weight << " " << "\n";
	}

	//******************************
	//Trigger Selection
	//******************************
	bool passTrigger = false;

	if ( processLabels[i] == "Data" || processLabels[i] == "FakeTemplate" ) {	  
	  double dataWeight = 1;
	  if (events->pho1.Pt() > 185) {
	    dataWeight = 1;
	    if (events->HLTDecision[102] && events->pho1HLTFilter[36]) passTrigger = true;
	  } 
          else if (events->pho1.Pt() > 138) {
            dataWeight = events->HLTPrescale[100];
            if (events->HLTDecision[100] && events->pho1HLTFilter[34]) passTrigger = true;
          } else if (events->pho1.Pt() > 108) {
            dataWeight = events->HLTPrescale[99];
            if (events->HLTDecision[99] && events->pho1HLTFilter[33]) passTrigger = true;
          } else if (events->pho1.Pt() > 88) {
            dataWeight = events->HLTPrescale[98];
            if (events->HLTDecision[98] && events->pho1HLTFilter[32]) passTrigger = true;
          } else if (events->pho1.Pt() > 60){
            dataWeight = events->HLTPrescale[97];
            if (events->HLTDecision[97] && events->pho1HLTFilter[31]) passTrigger = true;
          } else {
            dataWeight = events->HLTPrescale[96];
            if (events->HLTDecision[96] && events->pho1HLTFilter[30]) passTrigger = true;
          }

	  weight = dataWeight;

	} else {
	  //For Prompt Template from MC, don't require triggers because 
	  //the MC samples don't have proper trigger infomation
	  passTrigger = true;
	}

	if (!passTrigger) continue;

	//******************************
	//Selection Cuts 
	//******************************
	//Photon selection
	if (! (events->pho1.Pt() > 185)) continue;

	if ( processLabels[i] == "Data") {

	  // select photons in EB
	  if(fabs(events->pho1.Eta()) < 1.479 && events->pho1_sigmaietaieta < 0.015 && events->pho1_chargediso < 2.5) {
	    histSigmaIetaIeta_EB[i]->Fill(events->pho1_sigmaietaieta);
	      
	    for(int ii = 0; ii<NMRBins; ii++) {
	      for(int jj = 0; jj<NRsqBins; jj++) {
		if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] ) {
		  if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
		    histSigmaIetaIeta_EB_binned[ii][jj]->Fill(events->pho1_sigmaietaieta, weight);
		  }
		}
	      }
	    }
	    
	    for(int ii = 0; ii<NMRBins; ii++) {
	      if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] && events->Rsq_NoPho > 0.15 ) {
		histSigmaIetaIeta_EB_MRbinned[ii]->Fill(events->pho1_sigmaietaieta, weight);
	      }
	    }
	    for(int jj = 0; jj<NRsqBins; jj++) {
	      if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] && events->MR_NoPho > 400. ){
		histSigmaIetaIeta_EB_Rsqbinned[jj]->Fill(events->pho1_sigmaietaieta, weight);
	      }					      
	    }
	  }
	  // select photons in EE
	  if(fabs(events->pho1.Eta()) > 1.479 && events->pho1_sigmaietaieta > 0.015 && events->pho1_chargediso < 2.5) {
	    histSigmaIetaIeta_EE[i]->Fill(events->pho1_sigmaietaieta);
	      
	    for(int ii = 0; ii<NMRBins; ii++) {
	      for(int jj = 0; jj<NRsqBins; jj++) {
		if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] ) {
		  if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
		    histSigmaIetaIeta_EE_binned[ii][jj]->Fill(events->pho1_sigmaietaieta, weight);
		  }
		}
	      }
	    }
	    
	    for(int ii = 0; ii<NMRBins; ii++) {
	      if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] && events->Rsq_NoPho > 0.15 ) {
		histSigmaIetaIeta_EE_MRbinned[ii]->Fill(events->pho1_sigmaietaieta, weight);
	      }
	    }
	    for(int jj = 0; jj<NRsqBins; jj++) {
	      if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] && events->MR_NoPho > 400. ){
		histSigmaIetaIeta_EE_Rsqbinned[jj]->Fill(events->pho1_sigmaietaieta, weight);
	      }					      
	    }
	  }	  
	} // end if Data

	if ( processLabels[i] == "FakeTemplate") {
	  if( fabs(events->pho1.Eta()) < 1.479 ) { // make the fake template for barrel photons
	    if( events->pho1_chargediso > 5.0 && events->pho1_sigmaietaieta < 0.015 ) {
	      histSigmaIetaIetaTemplate_EB[i]->Fill(events->pho1_sigmaietaieta);
	      for(int ii = 0; ii<NMRBins; ii++) {
		for(int jj = 0; jj<NRsqBins; jj++) {
		  if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] ) {
		    if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
		      histSigmaIetaIetaFakeTemplate_EB_binned[ii][jj]->Fill(events->pho1_sigmaietaieta, weight);
		    }
		  }
		}
	      }
	      
	      for(int ii = 0; ii<NMRBins; ii++) {
		if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] && events->Rsq_NoPho > 0.15 ) {
		  histSigmaIetaIetaFakeTemplate_EB_MRbinned[ii]->Fill(events->pho1_sigmaietaieta, weight);
		}
	      }
	      for(int jj = 0; jj<NRsqBins; jj++) {
		if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] && events->MR_NoPho > 400. ){
		  histSigmaIetaIetaFakeTemplate_EB_Rsqbinned[jj]->Fill(events->pho1_sigmaietaieta, weight);
		}					      
	      }
	    }
	  }
	    
	  if( fabs(events->pho1.Eta()) > 1.479 ) { // make the fake template for endcap photons
	    if( events->pho1_chargediso > 5.0 && events->pho1_sigmaietaieta > 0.015 ) {
	      histSigmaIetaIetaTemplate_EE[i]->Fill(events->pho1_sigmaietaieta);
	      for(int ii = 0; ii<NMRBins; ii++) {
		for(int jj = 0; jj<NRsqBins; jj++) {
		  if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] ) {
		    if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
		      histSigmaIetaIetaFakeTemplate_EE_binned[ii][jj]->Fill(events->pho1_sigmaietaieta, weight);
		    }
		  }
		}
	      }
	      
	      for(int ii = 0; ii<NMRBins; ii++) {
		if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] && events->Rsq_NoPho > 0.15 ) {
		  histSigmaIetaIetaFakeTemplate_EE_MRbinned[ii]->Fill(events->pho1_sigmaietaieta, weight);
		}
	      }
	      for(int jj = 0; jj<NRsqBins; jj++) {
		if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] && events->MR_NoPho > 400. ){
		  histSigmaIetaIetaFakeTemplate_EE_Rsqbinned[jj]->Fill(events->pho1_sigmaietaieta, weight);
		}
	      }
	    }
	  }
	} // end if Fake Template

	if ( processLabels[i] == "PromptTemplate") {
	  if( fabs(events->pho1.Eta()) < 1.479 ) { // make the prompt template for barrel photons
	    if( events->pho1_chargediso < 2.5 && events->pho1_sigmaietaieta < 0.015 ) {
	      histSigmaIetaIetaTemplate_EB[i]->Fill(events->pho1_sigmaietaieta);
	      for(int ii = 0; ii<NMRBins; ii++) {
		for(int jj = 0; jj<NRsqBins; jj++) {
		  if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] ) {
		    if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
		      histSigmaIetaIetaPromptTemplate_EB_binned[ii][jj]->Fill(events->pho1_sigmaietaieta, weight);
		    }
		  }
		}
	      }
	      
	      for(int ii = 0; ii<NMRBins; ii++) {
		if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] && events->Rsq_NoPho > 0.15 ) {
		  histSigmaIetaIetaPromptTemplate_EB_MRbinned[ii]->Fill(events->pho1_sigmaietaieta, weight);
		}
	      }
	      for(int jj = 0; jj<NRsqBins; jj++) {
		if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] && events->MR_NoPho > 400. ){
		  histSigmaIetaIetaPromptTemplate_EB_Rsqbinned[jj]->Fill(events->pho1_sigmaietaieta, weight);
		}					      
	      }
	    }
	  }
	    
	  if( fabs(events->pho1.Eta()) > 1.479 ) { // make the prompt template for endcap photons
	    if( events->pho1_chargediso < 2.5 && events->pho1_sigmaietaieta > 0.015 ) {
	      histSigmaIetaIetaTemplate_EE[i]->Fill(events->pho1_sigmaietaieta);
	      for(int ii = 0; ii<NMRBins; ii++) {
		for(int jj = 0; jj<NRsqBins; jj++) {
		  if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] ) {
		    if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
		      histSigmaIetaIetaPromptTemplate_EE_binned[ii][jj]->Fill(events->pho1_sigmaietaieta, weight);
		    }
		  }
		}
	      }
	      
	      for(int ii = 0; ii<NMRBins; ii++) {
		if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] && events->Rsq_NoPho > 0.15 ) {
		  histSigmaIetaIetaPromptTemplate_EE_MRbinned[ii]->Fill(events->pho1_sigmaietaieta, weight);
		}
	      }
	      for(int jj = 0; jj<NRsqBins; jj++) {
		if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] && events->MR_NoPho > 400. ){
		  histSigmaIetaIetaPromptTemplate_EE_Rsqbinned[jj]->Fill(events->pho1_sigmaietaieta, weight);
		}
	      }
	    }
	  }
	}
	  	          
      } //loop over events
    } //loop over input files
  } //loop over input file groups

  //****************************
  //Add CMS and Lumi Labels
  //****************************
  lumi_13TeV = "26.4 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.15;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  lumiTextSize = 0.4;

  double rangeMax=0.02;
  double rangeMin=0.005;
  /*
    RooRealVar SIeta("SIeta","#sigma_{i#etai#eta}",rangeMin,rangeMax);
    RooDataHist dataSR("dataSR","dataSR",SIeta,Import(*histSigmaIetaIeta_EB[0]));
   
    RooDataHist MCprompt("MCprompt","MCprompt",SIeta,Import(*promptHist));
    RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaTemplate_EB[0]));
   
    RooHistPdf PromptPDF("PromptPDF","PromptPDF",SIeta,MCprompt,0);
    RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",SIeta,MCnonPrompt,0);
   
    RooRealVar Fitfrac("Fitfrac","Fitfrac",0.7,0.,1.5);
   
    RooAddPdf PDF("PDF","PSR-NPSB-PDF", RooArgList(PromptPDF,nonPromptPDF),Fitfrac);
  
    TPaveText *tpav_txt = new TPaveText(0.43043478,0.6548342,0.7652174,0.8510471,"brNDC");
   
    tpav_txt->SetBorderSize(0);
    tpav_txt->SetFillStyle(0);
    tpav_txt->SetTextAlign(11);
    tpav_txt->SetTextFont(42);
    tpav_txt->SetTextSize(0.03);
    // double pureMCPurity=hpromptSR_Scut->Integral()/hpromptPlusnonPromptSR_Scut->Integral();
   
   
    // char mcPurity[100];
    // sprintf(mcPurity,"Purity(MC, EB): %4.3f ",pureMCPurity); 
    // tpav_txt->AddText(mcPurity);
   
   
    RooPlot* frame4 = SIeta.frame(Title("CMS #it{Preliminary}                        2.3 fb^{-1}, 13 TeV"));

    dataSR.plotOn(frame4,Name("dataSR1"));  
    PDF.fitTo(dataSR);
    PDF.plotOn(frame4);
   
   
  
    PDF.plotOn(frame4,Name("PromptPDF"),Components(RooArgList(PromptPDF)), LineColor(kRed),LineStyle(2), LineWidth(2));
    PDF.plotOn(frame4,Name("nonPromptPDF"),Components(RooArgList(nonPromptPDF)), LineColor(kGreen),LineStyle(2), LineWidth(3));
    PDF.plotOn(frame4,Name("PDF"),Components(RooArgList(PDF)), LineColor(kBlue),LineStyle(1), LineWidth(2));
   


    TLegend *leg=new TLegend(00.5402358,0.4731495,0.7441038,0.6734398,NULL,"brNDC");
    leg->SetTextFont(62);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(3);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->SetShadowColor(0);
    leg->SetDrawOption(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.02);

  
    leg->AddEntry(frame4->findObject("PromptPDF"),"Prompt PDF(from MC),EB","l");
    leg->AddEntry(frame4->findObject("nonPromptPDF"),"NonPrompt PDF(from MC),EB","l");
    leg->AddEntry(frame4->findObject("PDF"),"P+NP PDF(from MC),EB","l");
    leg->AddEntry(frame4->findObject("dataSR1"),"Data,EB","P");

  
  
  


    SIeta.setRange("SR",rangeMin,0.0107);
    RooAbsReal* fracP = PromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
    RooAbsReal* Ip=PromptPDF.createIntegral(SIeta);
    RooAbsReal* fracNP = nonPromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
    RooAbsReal* Inp=nonPromptPDF.createIntegral(SIeta);

    double f=Fitfrac.getVal();
    double fError=Fitfrac.getError();
    double dataAreaFull=histSigmaIetaIeta_EB[0]->Integral();
    double promptIntFull=f*dataAreaFull;
    double nonpromptIntFull=(1-f)*dataAreaFull;

    double fracP_val=fracP->getVal();
  
    double fracNP_val=fracNP->getVal();
  
    double promptIntS=fracP_val*f;
    double nonpromptIntS=fracNP_val*(1-f);
    double PurityFromFit=promptIntS/(promptIntS+nonpromptIntS);
  
    float SigmaErr=(fError/f)*PurityFromFit;

    cout<<"Beta error: "<<SigmaErr<<endl;

    cout<<"fracP:  "<<fracP_val<<endl;
    cout<<"fracNP:  "<<fracNP_val<<endl;
    cout<<"Purity From Fit: "<<PurityFromFit<<endl;

  


    char RfracFit[100];
    sprintf(RfracFit,"P(Fit to EB Data): %4.3f +- %4.3f ",PurityFromFit,SigmaErr); 
    tpav_txt->AddText(RfracFit);

  

  
    
    float chosenEff = 0.0107;
    TLine *line =new TLine(chosenEff,0,chosenEff,100);
    line->SetLineColor(kOrange+8);
    line->SetLineWidth(2);
    


    TCanvas *cFit=new TCanvas("cFit","cFit",850,850);
   
    cFit->cd();
    gPad->Update();
    gPad->SetLogy();
    cFit->Range(0.01,0.01,0.02,1000);
    frame4->Draw();
    tpav_txt->Draw();
    line->Draw();
    leg->Draw();
    cFit->SaveAs("PurityFitEB.png");
    cFit->SaveAs("PurityFitEB.pdf");
    cFit->SaveAs("PurityFitEB.gif");
  */


  // now do the fits in the MR/Rsq bins   
  for(int ii=0; ii<NMRBins; ii++)
    for(int jj=0; jj<NRsqBins; jj++) {
  // for(int ii=0; ii<1; ii++)
  //   for(int jj=0; jj<1; jj++) {
      RooRealVar SIeta("SIeta","#sigma_{i#etai#eta}",rangeMin,rangeMax);
      RooDataHist dataSR("dataSR","dataSR",SIeta,Import(*histSigmaIetaIeta_EB_binned[ii][jj]));
   
      RooDataHist MCprompt("MCprompt","MCprompt",SIeta,Import(*histSigmaIetaIetaPromptTemplate_EB_binned[ii][jj]));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaTemplate_EB[0]));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*fakeHist_EB));
      RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaFakeTemplate_EB_binned[ii][jj]));

      RooHistPdf PromptPDF("PromptPDF","PromptPDF",SIeta,MCprompt,0);
      RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",SIeta,MCnonPrompt,0);
   
      RooRealVar Fitfrac("Fitfrac","Fitfrac",0.95,0.50,1.0);
   
      RooAddPdf PDF("PDF","PSR-NPSB-PDF", RooArgList(PromptPDF,nonPromptPDF),Fitfrac);
  
      TPaveText *tpav_txt = new TPaveText(0.5,0.7,0.85,0.9,"brNDC");
   
      tpav_txt->SetBorderSize(0);
      tpav_txt->SetFillStyle(0);
      tpav_txt->SetTextAlign(11);
      tpav_txt->SetTextFont(42);
      tpav_txt->SetTextSize(0.03);
      // double pureMCPurity=hpromptSR_Scut->Integral()/hpromptPlusnonPromptSR_Scut->Integral();
   
   
      // char mcPurity[100];
      // sprintf(mcPurity,"Purity(MC, EB): %4.3f ",pureMCPurity); 
      // tpav_txt->AddText(mcPurity);
   
   
      RooPlot* frame4 = SIeta.frame(Title(""));
      frame4->SetTitle("");

      dataSR.plotOn(frame4,Name("dataSR1"));  
      PDF.fitTo(dataSR);
      PDF.plotOn(frame4);      
  
      PDF.plotOn(frame4,Name("PromptPDF"),Components(RooArgList(PromptPDF)), LineColor(kRed),LineStyle(2), LineWidth(2));
      PDF.plotOn(frame4,Name("nonPromptPDF"),Components(RooArgList(nonPromptPDF)), LineColor(kGreen),LineStyle(2), LineWidth(3));
      PDF.plotOn(frame4,Name("PDF"),Components(RooArgList(PDF)), LineColor(kBlue),LineStyle(1), LineWidth(2));   

      TLegend *leg=new TLegend(0.6,0.55,0.85,0.75,NULL,"brNDC");
      leg->SetTextFont(62);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(3);
      leg->SetFillColor(0);
      leg->SetFillStyle(1001);
      leg->SetShadowColor(0);
      leg->SetDrawOption(0);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.02);
  
      leg->AddEntry(frame4->findObject("PromptPDF"),"Prompt Template, EB","l");
      leg->AddEntry(frame4->findObject("nonPromptPDF"),"Fak Template, EB","l");
      leg->AddEntry(frame4->findObject("PDF"),"Total Fit, EB","l");
      leg->AddEntry(frame4->findObject("dataSR1"),"Data, EB","P");      

      SIeta.setRange("SR",rangeMin,0.0103);
      RooAbsReal* fracP = PromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
      RooAbsReal* Ip=PromptPDF.createIntegral(SIeta);
      RooAbsReal* fracNP = nonPromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
      RooAbsReal* Inp=nonPromptPDF.createIntegral(SIeta);

      double f=Fitfrac.getVal();
      double fError=Fitfrac.getError();
      double dataAreaFull=histSigmaIetaIeta_EB_binned[ii][jj]->Integral();
      double promptIntFull=f*dataAreaFull;
      double nonpromptIntFull=(1-f)*dataAreaFull;

      double fracP_val=fracP->getVal();
  
      double fracNP_val=fracNP->getVal();
  
      double promptIntS=fracP_val*f;
      double nonpromptIntS=fracNP_val*(1-f);
      double PurityFromFit=promptIntS/(promptIntS+nonpromptIntS);
  
      float SigmaErr=(fError/f)*PurityFromFit;

      cout<<"Beta error: "<<SigmaErr<<endl;

      cout<<"fracP:  "<<fracP_val<<endl;
      cout<<"fracNP:  "<<fracNP_val<<endl;
      cout<<"Purity From Fit: "<<PurityFromFit<<endl;
  
      histSigmaIetaIeta_EB_MRRsq->SetBinContent(ii+1, jj+1, PurityFromFit);

      char RfracFit[100];
      sprintf(RfracFit,"Purity in EB: %4.3f +- %4.3f ",PurityFromFit,SigmaErr); 
      tpav_txt->AddText(RfracFit);

      float chosenEff = 0.0103;
      TLine *line =new TLine(chosenEff,0,chosenEff,100);
      line->SetLineColor(kOrange+8);
      line->SetLineWidth(2);    

      TCanvas *cFit=new TCanvas("cFit","cFit",850,850);
   
      cFit->cd();
      gPad->Update();
      gPad->SetLogy();
      cFit->Range(0.01,0.01,0.02,1000);
      frame4->GetYaxis()->SetTitle("Number of events");
      frame4->GetYaxis()->SetTitleOffset(1.35);
      frame4->Draw();
      tpav_txt->Draw();
      line->Draw();
      leg->Draw();
      CMS_lumi(cFit,4,0);
      cFit->SaveAs(Form("PurityFitEB_%d_%d.png", ii, jj));
      cFit->SaveAs(Form("PurityFitEB_%d_%d.pdf", ii, jj));
    }
  
  //////////////////////////
  // now do the fits in the MR bins   
  for(int ii=0; ii<NMRBins; ii++)
  // for(int ii=0; ii<1; ii++)
    {     
      cout<<"MR BINNED FIT: "<<ii<<endl;

      RooRealVar SIeta("SIeta","#sigma_{i#etai#eta}",rangeMin,rangeMax);
      RooDataHist dataSR("dataSR","dataSR",SIeta,Import(*histSigmaIetaIeta_EB_MRbinned[ii]));
      
      RooDataHist MCprompt("MCprompt","MCprompt",SIeta,Import(*histSigmaIetaIetaPromptTemplate_EB_MRbinned[ii]));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaTemplate_EB[0]));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*fakeHist_EB));
      RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaFakeTemplate_EB_MRbinned[ii]));
   
      RooHistPdf PromptPDF("PromptPDF","PromptPDF",SIeta,MCprompt,0);
      RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",SIeta,MCnonPrompt,0);
      
      RooRealVar Fitfrac("Fitfrac","Fitfrac",0.7,0.,1.0);
      
      RooAddPdf PDF("PDF","PSR-NPSB-PDF", RooArgList(PromptPDF,nonPromptPDF),Fitfrac);
      
      TPaveText *tpav_txt = new TPaveText(0.5,0.7,0.85,0.9,"brNDC");
      
      tpav_txt->SetBorderSize(0);
      tpav_txt->SetFillStyle(0);
      tpav_txt->SetTextAlign(11);
      tpav_txt->SetTextFont(42);
      tpav_txt->SetTextSize(0.03);
      // double pureMCPurity=hpromptSR_Scut->Integral()/hpromptPlusnonPromptSR_Scut->Integral();
      
   
      // char mcPurity[100];
      // sprintf(mcPurity,"Purity(MC, EB): %4.3f ",pureMCPurity); 
      // tpav_txt->AddText(mcPurity);
   
   
      RooPlot* frame4 = SIeta.frame(Title(""));
      frame4->SetTitle("");
      
      dataSR.plotOn(frame4,Name("dataSR1"));  
      PDF.fitTo(dataSR);
      PDF.plotOn(frame4);
      
      PDF.plotOn(frame4,Name("PromptPDF"),Components(RooArgList(PromptPDF)), LineColor(kRed),LineStyle(2), LineWidth(2));
      PDF.plotOn(frame4,Name("nonPromptPDF"),Components(RooArgList(nonPromptPDF)), LineColor(kGreen),LineStyle(2), LineWidth(3));
      PDF.plotOn(frame4,Name("PDF"),Components(RooArgList(PDF)), LineColor(kBlue),LineStyle(1), LineWidth(2));
      
      TLegend *leg=new TLegend(0.6,0.55,0.85,0.75,NULL,"brNDC");
      leg->SetTextFont(62);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(3);
      leg->SetFillColor(0);
      leg->SetFillStyle(1001);
      leg->SetShadowColor(0);
      leg->SetDrawOption(0);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.02);
      
      
      leg->AddEntry(frame4->findObject("PromptPDF"),"Prompt Template, EB","l");
      leg->AddEntry(frame4->findObject("nonPromptPDF"),"Fake Template, EB","l");
      leg->AddEntry(frame4->findObject("PDF"),"Total Fit, EB","l");
      leg->AddEntry(frame4->findObject("dataSR1"),"Data,EB","P");      

      SIeta.setRange("SR",rangeMin,0.0103);
      RooAbsReal* fracP = PromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
      RooAbsReal* Ip=PromptPDF.createIntegral(SIeta);
      RooAbsReal* fracNP = nonPromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
      RooAbsReal* Inp=nonPromptPDF.createIntegral(SIeta);
      
      double f=Fitfrac.getVal();
      double fError=Fitfrac.getError();
      double dataAreaFull=histSigmaIetaIeta_EB_MRbinned[ii]->Integral();
      double promptIntFull=f*dataAreaFull;
      double nonpromptIntFull=(1-f)*dataAreaFull;
      
      double fracP_val=fracP->getVal();
      
      double fracNP_val=fracNP->getVal();
      
      double promptIntS=fracP_val*f;
      double nonpromptIntS=fracNP_val*(1-f);
      double PurityFromFit=promptIntS/(promptIntS+nonpromptIntS);
      
      float SigmaErr=(fError/f)*PurityFromFit;
      
      cout<<"Beta error: "<<SigmaErr<<endl;
      
      cout<<"fracP:  "<<fracP_val<<endl;
      cout<<"fracNP:  "<<fracNP_val<<endl;
      cout<<"Purity From Fit: "<<PurityFromFit<<endl;
      
      
      histSigmaIetaIeta_EB_MR->SetBinContent(ii+1, PurityFromFit);
      histSigmaIetaIeta_EB_MR->SetBinError(ii+1, SigmaErr);
      
      TCanvas *c1=new TCanvas("c1","c1",850,850);      

      TH1F *h1 = new TH1F("h1","h1",1000,0, 4000);
      for(int a = 0; a < 1000; a++) { h1->SetBinContent(a, 0.95); h1->SetBinError(a, 0.05); }
      h1->SetFillStyle(3356);
      h1->SetFillColor(4);
      histSigmaIetaIeta_EB_MR->SetBinContent(ii+1, PurityFromFit);
      histSigmaIetaIeta_EB_MR->SetBinError(ii+1, SigmaErr);
      histSigmaIetaIeta_EB_MR->SetMarkerStyle(8);
      histSigmaIetaIeta_EB_MR->SetStats(0);
      histSigmaIetaIeta_EB_MR->SetLineColor(2);
      histSigmaIetaIeta_EB_MR->SetLineWidth(2);
      histSigmaIetaIeta_EB_MR->SetTitle("");
      histSigmaIetaIeta_EB_MR->GetXaxis()->SetTitle("MR [GeV]");
      histSigmaIetaIeta_EB_MR->GetYaxis()->SetTitle("Photon Purity");      
      histSigmaIetaIeta_EB_MR->GetYaxis()->SetTitleOffset(1.35);
      histSigmaIetaIeta_EB_MR->GetYaxis()->SetRangeUser(0.75, 1.05);
      histSigmaIetaIeta_EB_MR->GetXaxis()->SetRangeUser(400, 1200.);
      histSigmaIetaIeta_EB_MR->Draw("e");
      h1->Draw("e3 same");       
      histSigmaIetaIeta_EB_MR->Draw("e same");

      CMS_lumi(c1,4,0);

      c1->cd();
      c1->Update();
      gPad->Update();

      c1->SaveAs("Purity_EB_MR.png");
      c1->SaveAs("Purity_EB_MR.pdf");
      c1->SaveAs("Purity_EB_MR.root");
      
      char RfracFit[100];
      sprintf(RfracFit,"#splitline{%4.0f<M_{R}<%4.0f}{Purity in EB : %4.3f +- %4.3f}", MRBins[ii], MRBins[ii+1], PurityFromFit,SigmaErr); 
      tpav_txt->AddText(RfracFit);              
    
      float chosenEff = 0.0103;
      TLine *line =new TLine(chosenEff,0,chosenEff,100);
      line->SetLineColor(kOrange+8);
      line->SetLineWidth(2);            
      
      TCanvas *cFit=new TCanvas("cFit","cFit",850,850);
      
      cFit->cd();
      gPad->Update();
      gPad->SetLogy();
      cFit->Range(0.01,0.01,0.02,1000);
      frame4->GetYaxis()->SetTitle("Number of events");
      frame4->GetYaxis()->SetTitleOffset(1.35);
      frame4->Draw();
      tpav_txt->Draw();
      line->Draw();
      leg->Draw();
      CMS_lumi(cFit,4,0);
      cFit->SaveAs(Form("PurityFitEB_MRbinned_%d.png", ii));
      cFit->SaveAs(Form("PurityFitEB_MRbinned_%d.pdf", ii));
    }   
  

  //   // now do the fits in the Rsq bins   
  for(int ii=0; ii<NRsqBins; ii++)
    {     
      cout<<"Rsq BINNED FIT: "<<ii<<endl;
      
      RooRealVar SIeta("SIeta","#sigma_{i#etai#eta}",rangeMin,rangeMax);
      RooDataHist dataSR("dataSR","dataSR",SIeta,Import(*histSigmaIetaIeta_EB_Rsqbinned[ii]));
      
      RooDataHist MCprompt("MCprompt","MCprompt",SIeta,Import(*histSigmaIetaIetaPromptTemplate_EB_Rsqbinned[ii]));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaTemplate_EB[0]));
       // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*fakeHist_EB));
      RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaFakeTemplate_EB_Rsqbinned[ii]));
     
      RooHistPdf PromptPDF("PromptPDF","PromptPDF",SIeta,MCprompt,0);
      RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",SIeta,MCnonPrompt,0);
      
      RooRealVar Fitfrac("Fitfrac","Fitfrac",0.7,0.,1.0);
      
      RooAddPdf PDF("PDF","PSR-NPSB-PDF", RooArgList(PromptPDF,nonPromptPDF),Fitfrac);
      
      TPaveText *tpav_txt = new TPaveText(0.5,0.7,0.85,0.9,"brNDC");
      
      tpav_txt->SetBorderSize(0);
      tpav_txt->SetFillStyle(0);
      tpav_txt->SetTextAlign(11);
      tpav_txt->SetTextFont(42);
      tpav_txt->SetTextSize(0.03);
      // double pureMCPurity=hpromptSR_Scut->Integral()/hpromptPlusnonPromptSR_Scut->Integral();
      
      
      // char mcPurity[100];
      // sprintf(mcPurity,"Purity(MC, EB): %4.3f ",pureMCPurity); 
      // tpav_txt->AddText(mcPurity);
      
      
      RooPlot* frame4 = SIeta.frame(Title(""));
      frame4->SetTitle("");
  
      dataSR.plotOn(frame4,Name("dataSR1"));  
      PDF.fitTo(dataSR);
      PDF.plotOn(frame4);            
      
      PDF.plotOn(frame4,Name("PromptPDF"),Components(RooArgList(PromptPDF)), LineColor(kRed),LineStyle(2), LineWidth(2));
      PDF.plotOn(frame4,Name("nonPromptPDF"),Components(RooArgList(nonPromptPDF)), LineColor(kGreen),LineStyle(2), LineWidth(3));
      PDF.plotOn(frame4,Name("PDF"),Components(RooArgList(PDF)), LineColor(kBlue),LineStyle(1), LineWidth(2));            
      
      TLegend *leg=new TLegend(0.6,0.55,0.85,0.75,NULL,"brNDC");
      leg->SetTextFont(62);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(3);
      leg->SetFillColor(0);
      leg->SetFillStyle(1001);
      leg->SetShadowColor(0);
      leg->SetDrawOption(0);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.02);
            
      leg->AddEntry(frame4->findObject("PromptPDF"),"Prompt Template, EB","l");
      leg->AddEntry(frame4->findObject("nonPromptPDF"),"Fake Template, EB","l");
      leg->AddEntry(frame4->findObject("PDF"),"Total Fit, EB","l");
      leg->AddEntry(frame4->findObject("dataSR1"),"Data, EB","P");      
      
      SIeta.setRange("SR",rangeMin,0.0103);
      RooAbsReal* fracP = PromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
      RooAbsReal* Ip=PromptPDF.createIntegral(SIeta);
      RooAbsReal* fracNP = nonPromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
      RooAbsReal* Inp=nonPromptPDF.createIntegral(SIeta);
      
      double f=Fitfrac.getVal();
      double fError=Fitfrac.getError();
      double dataAreaFull=histSigmaIetaIeta_EB_Rsqbinned[ii]->Integral();
      double promptIntFull=f*dataAreaFull;
      double nonpromptIntFull=(1-f)*dataAreaFull;
      
      double fracP_val=fracP->getVal();
      
      double fracNP_val=fracNP->getVal();
      
      double promptIntS=fracP_val*f;
      double nonpromptIntS=fracNP_val*(1-f);
      double PurityFromFit=promptIntS/(promptIntS+nonpromptIntS);
      
      float SigmaErr=(fError/f)*PurityFromFit;
      
      cout<<"Beta error: "<<SigmaErr<<endl;
      
      cout<<"fracP:  "<<fracP_val<<endl;
      cout<<"fracNP:  "<<fracNP_val<<endl;
      cout<<"Purity From Fit: "<<PurityFromFit<<endl;
            
      histSigmaIetaIeta_EB_Rsq->SetBinContent(ii+1, PurityFromFit);
      histSigmaIetaIeta_EB_Rsq->SetBinError(ii+1, SigmaErr);

      TCanvas *c1=new TCanvas("c1","c1",850,850);      
      c1->cd();
      gPad->Update();

      TH1F *h1 = new TH1F("h1","h1",1000,0, 1.5);
      for(int a = 0; a < 1000; a++) { h1->SetBinContent(a, 0.95); h1->SetBinError(a, 0.05); }
      h1->SetFillStyle(3356);
      h1->SetFillColor(4);
      histSigmaIetaIeta_EB_Rsq->SetBinContent(ii+1, PurityFromFit);
      histSigmaIetaIeta_EB_Rsq->SetBinError(ii+1, SigmaErr);
      histSigmaIetaIeta_EB_Rsq->SetMarkerStyle(8);
      histSigmaIetaIeta_EB_Rsq->SetStats(0);
      histSigmaIetaIeta_EB_Rsq->SetLineColor(2);
      histSigmaIetaIeta_EB_Rsq->SetLineWidth(2);
      histSigmaIetaIeta_EB_Rsq->SetTitle("");
      histSigmaIetaIeta_EB_Rsq->GetXaxis()->SetTitle("Rsq");
      histSigmaIetaIeta_EB_Rsq->GetYaxis()->SetTitle("Photon Purity");      
      histSigmaIetaIeta_EB_Rsq->GetYaxis()->SetRangeUser(0.75, 1.05);
      histSigmaIetaIeta_EB_Rsq->GetYaxis()->SetTitleOffset(1.35);
      histSigmaIetaIeta_EB_Rsq->GetXaxis()->SetRangeUser(0.15, 0.52);
      histSigmaIetaIeta_EB_Rsq->Draw("e");
      h1->Draw("e3 same");       
      histSigmaIetaIeta_EB_Rsq->Draw("e same");
      CMS_lumi(c1,4,0);
      c1->SaveAs("Purity_EB_Rsq.png");
      c1->SaveAs("Purity_EB_Rsq.pdf");
      c1->SaveAs("Purity_EB_Rsq.root");
            
      char RfracFit[100];
      sprintf(RfracFit,"#splitline{%4.2f<R^{2}<%4.2f}{Purity in EB : %4.3f +- %4.3f}", RsqBins[ii], RsqBins[ii+1], PurityFromFit,SigmaErr); 
      tpav_txt->AddText(RfracFit);
    
      float chosenEff = 0.0103;
      TLine *line =new TLine(chosenEff,0,chosenEff,100);
      line->SetLineColor(kOrange+8);
      line->SetLineWidth(2);
      
      TCanvas *cFit=new TCanvas("cFit","cFit",850,850);
      
      cFit->cd();
      gPad->Update();
      gPad->SetLogy();
      cFit->Range(0.01,0.01,0.02,1000);
      frame4->GetYaxis()->SetTitle("Number of events");
      frame4->GetYaxis()->SetTitleOffset(1.35);
      frame4->Draw();
      tpav_txt->Draw();
      line->Draw();
      leg->Draw();
      CMS_lumi(cFit,4,0);
      cFit->SaveAs(Form("PurityFitEB_Rsqbinned_%d.png", ii));
      cFit->SaveAs(Form("PurityFitEB_Rsqbinned_%d.pdf", ii));
    }   

  // //////////////////////////

  // ///////////////////////
  // ///// ENDCAP //////////
  // ///////////////////////

  // TObjArray *mc_ee = new TObjArray(2);        // template histograms are put in this array
  // mc_ee->Add(histSigmaIetaIetaTemplate_EE[0]); //fakes
  // mc_ee->Add(promptHist_EE); //prompt
  
  // TFractionFitter* fit_ee = new TFractionFitter(histSigmaIetaIeta_EE[0], mc_ee); // initialise
  
  // lowbound   = histSigmaIetaIeta_EE[0]->GetXaxis()->FindBin(0.0195)+1;
  // upperbound = histSigmaIetaIeta_EE[0]->GetXaxis()->FindBin(0.0271)-1;

  // // fit_ee->SetRangeX(lowbound, NBinsSigmaietaieta);       

  // status = fit_ee->Fit();               // perform the fit
  // std::cout << "fit status EE: " << status << std::endl;
  // if (status == 0) {                       // check on fit status
  //   TCanvas * c1 = new TCanvas("c", "c", 600, 600) ;
  //   c1->cd();
  //   TPad pad1("pad1","pad1",0,0.4,1,1);
  //   pad1.SetBottomMargin(0);
  //   pad1.SetLogy();
  //   pad1.Draw();
  //   pad1.cd();
    
  //   TH1F* result = (TH1F*) fit_ee->GetPlot();
  //   result->SetLineColor(2);
    
  //   fit_ee->GetResult( 0, p0, errP0);
  //   cout<<"Fakes Fraction: "<< p0 <<"+/-"<<errP0<<endl;
  //   fit_ee->GetResult( 1, p1, errP1);
  //   cout<<"Prompt Fraction: "<< p1 <<"+/-"<<errP1<<endl;
    
  //   histSigmaIetaIeta_EE[0]->Draw("Ep");
  //   result->Draw("samehist");
    
  //   TH1F* resultfake = (TH1F*) histSigmaIetaIetaTemplate_EE[0]->Clone();
  //   resultfake->Scale((p0*histSigmaIetaIeta_EE[0]->Integral())/resultfake->Integral()); // normalize the Template of the fakes
  //   resultfake->SetLineColor(3);
    
  //   TH1F* resultprompt = (TH1F*) promptHist_EE->Clone();
  //   resultprompt->Scale((p1*histSigmaIetaIeta_EE[0]->Integral())/resultprompt->Integral()); // normalize the Template of the prompts
  //   resultprompt->SetLineColor(4);
    
  //   Double_t TotalFakes  = resultfake->Integral(0, upperbound);
  //   Double_t TotalPrompt = resultprompt->Integral(0, upperbound);
  //   Double_t Total       = result->Integral(0, upperbound);
    
  //   cout<<"Fake amount passing sigma_ietaieta cut, INCLUSIVE: "<<TotalFakes<<endl;
  //   cout<<"Prompt amount passing sigma_ietaieta cut: "<<TotalPrompt<<endl;
  //   cout<<"Total amount passing sigma_ietaieta cut: "<<Total<<endl;
    
  //   resultfake->Draw("samehist");
  //   resultprompt->Draw("samehist");
    
  //   pad1.Modified();
  //   gPad->Update();
    
  //   c1->cd();
    
  //   TH1F *resultcopy = (TH1F*)result->Clone();
  //   resultcopy -> Divide(histSigmaIetaIeta_EE[0]);
  //   resultcopy->GetYaxis()->SetTitle("Total Fit/Data");
  //   resultcopy->SetMinimum(0.2);
  //   resultcopy->SetMaximum(2.0);
    
  //   TPad pad2("pad2","pad2",0,0.0,1,0.4);
  //   pad2.SetTopMargin(0);
  //   pad2.SetTopMargin(0.008);
  //   pad2.SetBottomMargin(0.25);
  //   pad2.SetGridy();
  //   pad2.Draw();
  //   pad2.cd();
  //   resultcopy->Draw("pe");
  //   pad2.Modified();
  //   gPad->Update();
    
  //   c1->SetLogy();
  //   c1->SaveAs("Fit_ee.png");
  // }

  rangeMax = 0.05;
  rangeMin = 0.01;

  //  now do the fits in the MR/Rsq bins in EE
  for(int ii=0; ii<NMRBins; ii++)
    for(int jj=0; jj<NRsqBins; jj++) {
       
      cout<<"BINNED FIT in EE: "<<ii<<" "<<jj<<endl;
              
      RooRealVar SIeta("SIeta","#sigma_{i#etai#eta}",rangeMin,rangeMax);
      RooDataHist dataSR("dataSR","dataSR",SIeta,Import(*histSigmaIetaIeta_EE_binned[ii][jj]));
	 
      RooDataHist MCprompt("MCprompt","MCprompt",SIeta,Import(*histSigmaIetaIetaPromptTemplate_EE_binned[ii][jj]));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaTemplate_EE[0]));
       // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*fakeHist_EE));
      RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaFakeTemplate_EE_binned[ii][jj]));
  
      RooHistPdf PromptPDF("PromptPDF","PromptPDF",SIeta,MCprompt,0);
      RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",SIeta,MCnonPrompt,0);
   
      RooRealVar Fitfrac("Fitfrac","Fitfrac",0.7,0.,1.0);
   
      RooAddPdf PDF("PDF","PSR-NPSB-PDF", RooArgList(PromptPDF,nonPromptPDF),Fitfrac);

      TPaveText *tpav_txt = new TPaveText(0.5,0.7,0.85,0.9,"brNDC");
   
      tpav_txt->SetBorderSize(0);
      tpav_txt->SetFillStyle(0);
      tpav_txt->SetTextAlign(11);
      tpav_txt->SetTextFont(42);
      tpav_txt->SetTextSize(0.03);
      // double pureMCPurity=hpromptSR_Scut->Integral()/hpromptPlusnonPromptSR_Scut->Integral();
   
   
      // char mcPurity[100];
      // sprintf(mcPurity,"Purity(MC, EB): %4.3f ",pureMCPurity); 
      // tpav_txt->AddText(mcPurity);
   
      RooPlot* frame4 = SIeta.frame(Title(""));
      frame4->SetTitle("");

      dataSR.plotOn(frame4,Name("dataSR1"));  

      PDF.fitTo(dataSR);
      PDF.plotOn(frame4);      

      PDF.plotOn(frame4,Name("PromptPDF"),Components(RooArgList(PromptPDF)), LineColor(kRed),LineStyle(2), LineWidth(2));
      PDF.plotOn(frame4,Name("nonPromptPDF"),Components(RooArgList(nonPromptPDF)), LineColor(kGreen),LineStyle(2), LineWidth(3));
      PDF.plotOn(frame4,Name("PDF"),Components(RooArgList(PDF)), LineColor(kBlue),LineStyle(1), LineWidth(2));   

      TLegend *leg=new TLegend(0.6,0.55,0.85,0.75,NULL,"brNDC");
      leg->SetTextFont(62);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(3);
      leg->SetFillColor(0);
      leg->SetFillStyle(1001);
      leg->SetShadowColor(0);
      leg->SetDrawOption(0);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.02);

      leg->AddEntry(frame4->findObject("PromptPDF"),"Prompt Template, EE","l");
      leg->AddEntry(frame4->findObject("nonPromptPDF"),"Fak Template, EE","l");
      leg->AddEntry(frame4->findObject("PDF"),"Total Fit, EE","l");
      leg->AddEntry(frame4->findObject("dataSR1"),"Data, EE","P");      

      SIeta.setRange("SR",rangeMin,0.0271);
      RooAbsReal* fracP = PromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
      RooAbsReal* Ip=PromptPDF.createIntegral(SIeta);
      RooAbsReal* fracNP = nonPromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
      RooAbsReal* Inp=nonPromptPDF.createIntegral(SIeta);

      double f=Fitfrac.getVal();
      double fError=Fitfrac.getError();
      double dataAreaFull=histSigmaIetaIeta_EE_binned[ii][jj]->Integral();
      double promptIntFull=f*dataAreaFull;
      double nonpromptIntFull=(1-f)*dataAreaFull;

      double fracP_val=fracP->getVal();
  
      double fracNP_val=fracNP->getVal();
  
      double promptIntS=fracP_val*f;
      double nonpromptIntS=fracNP_val*(1-f);
      double PurityFromFit=promptIntS/(promptIntS+nonpromptIntS);
  
      float SigmaErr=(fError/f)*PurityFromFit;

      cout<<"Beta error: "<<SigmaErr<<endl;

      cout<<"fracP:  "<<fracP_val<<endl;
      cout<<"fracNP:  "<<fracNP_val<<endl;
      cout<<"Purity From Fit: "<<PurityFromFit<<endl;
  
      histSigmaIetaIeta_EB_MRRsq->SetBinContent(ii+1, jj+1, PurityFromFit);

      char RfracFit[100];
      sprintf(RfracFit,"Purity in EE: %4.3f +- %4.3f ",PurityFromFit,SigmaErr); 
      tpav_txt->AddText(RfracFit);

      float chosenEff = 0.0271;
      TLine *line =new TLine(chosenEff,0,chosenEff,100);
      line->SetLineColor(kOrange+8);
      line->SetLineWidth(2);    

      TCanvas *cFit=new TCanvas("cFit","cFit",850,850);
   
      cFit->cd();
      gPad->Update();
      gPad->SetLogy();
      cFit->Range(0.01,0.01,0.02,1000);
      frame4->GetYaxis()->SetTitle("Number of events");
      frame4->GetYaxis()->SetTitleOffset(1.35);
      frame4->Draw();
      tpav_txt->Draw();
      line->Draw();
      leg->Draw();
      CMS_lumi(cFit,4,0);
      cFit->SaveAs(Form("PurityFitEE_%d_%d.png", ii, jj));
      cFit->SaveAs(Form("PurityFitEE_%d_%d.pdf", ii, jj));


    }  

    // now do the fits in the MR bins   
  for(int ii=0; ii<NMRBins; ii++)
    {     
      cout<<"MR BINNED FIT: "<<ii<<endl;
      
      RooRealVar SIeta("SIeta","#sigma_{i#etai#eta}",rangeMin,rangeMax);
      RooDataHist dataSR("dataSR","dataSR",SIeta,Import(*histSigmaIetaIeta_EE_MRbinned[ii]));
      
      RooDataHist MCprompt("MCprompt","MCprompt",SIeta,Import(*histSigmaIetaIetaPromptTemplate_EE_MRbinned[ii]));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaTemplate_EE[0]));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*fakeHist_EE));
      RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaFakeTemplate_EE_MRbinned[ii]));
      
      RooHistPdf PromptPDF("PromptPDF","PromptPDF",SIeta,MCprompt,0);
      RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",SIeta,MCnonPrompt,0);
      
      RooRealVar Fitfrac("Fitfrac","Fitfrac",0.7,0.,1.0);
      
      RooAddPdf PDF("PDF","PSR-NPSB-PDF", RooArgList(PromptPDF,nonPromptPDF),Fitfrac);
      
      TPaveText *tpav_txt = new TPaveText(0.5,0.7,0.85,0.9,"brNDC");
      
      tpav_txt->SetBorderSize(0);
      tpav_txt->SetFillStyle(0);
      tpav_txt->SetTextAlign(11);
      tpav_txt->SetTextFont(42);
      tpav_txt->SetTextSize(0.03);
      // double pureMCPurity=hpromptSR_Scut->Integral()/hpromptPlusnonPromptSR_Scut->Integral();
      
   
      // char mcPurity[100];
      // sprintf(mcPurity,"Purity(MC, EB): %4.3f ",pureMCPurity); 
      // tpav_txt->AddText(mcPurity);
   
   
      RooPlot* frame4 = SIeta.frame(Title(""));
      frame4->SetTitle("");
      
      dataSR.plotOn(frame4,Name("dataSR1"));  
      PDF.fitTo(dataSR);
      PDF.plotOn(frame4);
      
      PDF.plotOn(frame4,Name("PromptPDF"),Components(RooArgList(PromptPDF)), LineColor(kRed),LineStyle(2), LineWidth(2));
      PDF.plotOn(frame4,Name("nonPromptPDF"),Components(RooArgList(nonPromptPDF)), LineColor(kGreen),LineStyle(2), LineWidth(3));
      PDF.plotOn(frame4,Name("PDF"),Components(RooArgList(PDF)), LineColor(kBlue),LineStyle(1), LineWidth(2));
      
      TLegend *leg=new TLegend(0.6,0.55,0.85,0.75,NULL,"brNDC");
      leg->SetTextFont(62);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(3);
      leg->SetFillColor(0);
      leg->SetFillStyle(1001);
      leg->SetShadowColor(0);
      leg->SetDrawOption(0);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.02);
      
      
      leg->AddEntry(frame4->findObject("PromptPDF"),"Prompt Template, EE","l");
      leg->AddEntry(frame4->findObject("nonPromptPDF"),"Fake Template, EE","l");
      leg->AddEntry(frame4->findObject("PDF"),"Total Fit, EE","l");
      leg->AddEntry(frame4->findObject("dataSR1"),"Data, EE","P");      

      SIeta.setRange("SR",rangeMin,0.0271);
      RooAbsReal* fracP = PromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
      RooAbsReal* Ip=PromptPDF.createIntegral(SIeta);
      RooAbsReal* fracNP = nonPromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
      RooAbsReal* Inp=nonPromptPDF.createIntegral(SIeta);
      
      double f=Fitfrac.getVal();
      double fError=Fitfrac.getError();
      double dataAreaFull=histSigmaIetaIeta_EE_MRbinned[ii]->Integral();
      double promptIntFull=f*dataAreaFull;
      double nonpromptIntFull=(1-f)*dataAreaFull;
      
      double fracP_val=fracP->getVal();
      
      double fracNP_val=fracNP->getVal();
      
      double promptIntS=fracP_val*f;
      double nonpromptIntS=fracNP_val*(1-f);
      double PurityFromFit=promptIntS/(promptIntS+nonpromptIntS);
      
      float SigmaErr=(fError/f)*PurityFromFit;
      
      cout<<"Beta error: "<<SigmaErr<<endl;
      
      cout<<"fracP:  "<<fracP_val<<endl;
      cout<<"fracNP:  "<<fracNP_val<<endl;
      cout<<"Purity From Fit: "<<PurityFromFit<<endl;
      

      TCanvas *c1=new TCanvas("c1","c1",850,850);      
      c1->cd();
      gPad->Update();

      TH1F *h1 = new TH1F("h1","h1",1000,0, 4000.0);
      for(int a = 0; a < 1000; a++) { h1->SetBinContent(a, 0.95); h1->SetBinError(a, 0.05); }
      h1->SetFillStyle(3356);
      h1->SetFillColor(4);
      histSigmaIetaIeta_EE_MR->SetBinContent(ii+1, PurityFromFit);
      histSigmaIetaIeta_EE_MR->SetBinError(ii+1, SigmaErr);
      histSigmaIetaIeta_EE_MR->SetMarkerStyle(8);
      histSigmaIetaIeta_EE_MR->SetStats(0);
      histSigmaIetaIeta_EE_MR->SetLineColor(2);
      histSigmaIetaIeta_EE_MR->SetLineWidth(2);
      histSigmaIetaIeta_EE_MR->SetTitle("");
      histSigmaIetaIeta_EE_MR->GetXaxis()->SetTitle("MR [GeV]");
      histSigmaIetaIeta_EE_MR->GetYaxis()->SetTitle("Photon Purity");      
      histSigmaIetaIeta_EE_MR->GetYaxis()->SetRangeUser(0.75, 1.05);
      histSigmaIetaIeta_EE_MR->GetYaxis()->SetTitleOffset(1.35);
      histSigmaIetaIeta_EE_MR->GetXaxis()->SetRangeUser(400, 1200);
      histSigmaIetaIeta_EE_MR->Draw("e");
      h1->Draw("e3 same");       
      histSigmaIetaIeta_EE_MR->Draw("e same");
      CMS_lumi(c1,4,0);
      c1->SaveAs("Purity_EE_MR.png");
      c1->SaveAs("Purity_EE_MR.pdf");
      c1->SaveAs("Purity_EE_MR.root");
      
      char RfracFit[100];
      sprintf(RfracFit,"#splitline{%4.0f<M_{R}<%4.0f}{Purity in EE : %4.3f +- %4.3f}", MRBins[ii], MRBins[ii+1], PurityFromFit,SigmaErr); 
      tpav_txt->AddText(RfracFit);              
    
      float chosenEff = 0.0271;
      TLine *line =new TLine(chosenEff,0,chosenEff,100);
      line->SetLineColor(kOrange+8);
      line->SetLineWidth(2);            
      
      TCanvas *cFit=new TCanvas("cFit","cFit",850,850);
      
      cFit->cd();
      gPad->Update();
      gPad->SetLogy();
      cFit->Range(0.01,0.01,0.02,1000);
      frame4->GetYaxis()->SetTitle("Number of events");
      frame4->GetYaxis()->SetTitleOffset(1.35);
      frame4->Draw();
      tpav_txt->Draw();
      line->Draw();
      leg->Draw();
      CMS_lumi(cFit,4,0);
      cFit->SaveAs(Form("PurityFitEE_MRbinned_%d.png", ii));
      cFit->SaveAs(Form("PurityFitEE_MRbinned_%d.pdf", ii));
    }   

    // now do the fits in the Rsq bins   
  for(int ii=0; ii<NRsqBins; ii++)
    {     
      cout<<"Rsq BINNED FIT: "<<ii<<endl;
      
      RooRealVar SIeta("SIeta","#sigma_{i#etai#eta}",rangeMin,rangeMax);
      RooDataHist dataSR("dataSR","dataSR",SIeta,Import(*histSigmaIetaIeta_EE_Rsqbinned[ii]));
      
      RooDataHist MCprompt("MCprompt","MCprompt",SIeta,Import(*histSigmaIetaIetaPromptTemplate_EE_Rsqbinned[ii]));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaTemplate_EE[0]));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*fakeHist_EE));
       RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",SIeta,Import(*histSigmaIetaIetaFakeTemplate_EE_Rsqbinned[ii]));
     
      RooHistPdf PromptPDF("PromptPDF","PromptPDF",SIeta,MCprompt,0);
      RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",SIeta,MCnonPrompt,0);
      
      RooRealVar Fitfrac("Fitfrac","Fitfrac",0.7,0.,1.0);
      
      RooAddPdf PDF("PDF","PSR-NPSB-PDF", RooArgList(PromptPDF,nonPromptPDF),Fitfrac);
      
      TPaveText *tpav_txt = new TPaveText(0.5,0.7,0.85,0.9,"brNDC");
      
      tpav_txt->SetBorderSize(0);
      tpav_txt->SetFillStyle(0);
      tpav_txt->SetTextAlign(11);
      tpav_txt->SetTextFont(42);
      tpav_txt->SetTextSize(0.03);
      // double pureMCPurity=hpromptSR_Scut->Integral()/hpromptPlusnonPromptSR_Scut->Integral();
      
   
      // char mcPurity[100];
      // sprintf(mcPurity,"Purity(MC, EB): %4.3f ",pureMCPurity); 
      // tpav_txt->AddText(mcPurity);
   
   
      RooPlot* frame4 = SIeta.frame(Title(""));
      frame4->SetTitle("");
      
      dataSR.plotOn(frame4,Name("dataSR1"));  
      PDF.fitTo(dataSR);
      PDF.plotOn(frame4);
      
      PDF.plotOn(frame4,Name("PromptPDF"),Components(RooArgList(PromptPDF)), LineColor(kRed),LineStyle(2), LineWidth(2));
      PDF.plotOn(frame4,Name("nonPromptPDF"),Components(RooArgList(nonPromptPDF)), LineColor(kGreen),LineStyle(2), LineWidth(3));
      PDF.plotOn(frame4,Name("PDF"),Components(RooArgList(PDF)), LineColor(kBlue),LineStyle(1), LineWidth(2));
      
      TLegend *leg=new TLegend(0.6,0.55,0.85,0.75,NULL,"brNDC");
      leg->SetTextFont(62);
      leg->SetLineColor(1);
      leg->SetLineStyle(1);
      leg->SetLineWidth(3);
      leg->SetFillColor(0);
      leg->SetFillStyle(1001);
      leg->SetShadowColor(0);
      leg->SetDrawOption(0);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.02);
      
      
      leg->AddEntry(frame4->findObject("PromptPDF"),"Prompt Template, EE","l");
      leg->AddEntry(frame4->findObject("nonPromptPDF"),"Fake Template, EE","l");
      leg->AddEntry(frame4->findObject("PDF"),"Total Fit, EE","l");
      leg->AddEntry(frame4->findObject("dataSR1"),"Data, EE","P");      

      SIeta.setRange("SR",rangeMin,0.0271);
      RooAbsReal* fracP = PromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
      RooAbsReal* Ip=PromptPDF.createIntegral(SIeta);
      RooAbsReal* fracNP = nonPromptPDF.createIntegral(SIeta,NormSet(SIeta),Range("SR"));
      RooAbsReal* Inp=nonPromptPDF.createIntegral(SIeta);
      
      double f=Fitfrac.getVal();
      double fError=Fitfrac.getError();
      double dataAreaFull=histSigmaIetaIeta_EE_Rsqbinned[ii]->Integral();
      double promptIntFull=f*dataAreaFull;
      double nonpromptIntFull=(1-f)*dataAreaFull;
      
      double fracP_val=fracP->getVal();
      
      double fracNP_val=fracNP->getVal();
      
      double promptIntS=fracP_val*f;
      double nonpromptIntS=fracNP_val*(1-f);
      double PurityFromFit=promptIntS/(promptIntS+nonpromptIntS);
      
      float SigmaErr=(fError/f)*PurityFromFit;
      
      cout<<"Beta error: "<<SigmaErr<<endl;
      
      cout<<"fracP:  "<<fracP_val<<endl;
      cout<<"fracNP:  "<<fracNP_val<<endl;
      cout<<"Purity From Fit: "<<PurityFromFit<<endl;
      
      histSigmaIetaIeta_EE_Rsq->SetBinContent(ii+1, PurityFromFit);
      histSigmaIetaIeta_EE_Rsq->SetBinError(ii+1, SigmaErr);

      TCanvas *c1=new TCanvas("c1","c1",850,850);      
      c1->cd();
      gPad->Update();

      TH1F *h1 = new TH1F("h1","h1",1000,0, 1.5);
      for(int a = 0; a < 1000; a++) { h1->SetBinContent(a, 0.95); h1->SetBinError(a, 0.05); }
      h1->SetFillStyle(3356);
      h1->SetFillColor(4);
      histSigmaIetaIeta_EE_Rsq->SetBinContent(ii+1, PurityFromFit);
      histSigmaIetaIeta_EE_Rsq->SetBinError(ii+1, SigmaErr);
      histSigmaIetaIeta_EE_Rsq->SetMarkerStyle(8);
      histSigmaIetaIeta_EE_Rsq->SetStats(0);
      histSigmaIetaIeta_EE_Rsq->SetLineColor(2);
      histSigmaIetaIeta_EE_Rsq->SetLineWidth(2);
      histSigmaIetaIeta_EE_Rsq->SetTitle("");
      histSigmaIetaIeta_EE_Rsq->GetXaxis()->SetTitle("Rsq");
      histSigmaIetaIeta_EE_Rsq->GetYaxis()->SetTitle("Photon Purity");      
      histSigmaIetaIeta_EE_Rsq->GetYaxis()->SetRangeUser(0.75, 1.05);
      histSigmaIetaIeta_EE_Rsq->GetYaxis()->SetTitleOffset(1.35);
      histSigmaIetaIeta_EE_Rsq->GetXaxis()->SetRangeUser(0.15, 0.52);
      histSigmaIetaIeta_EE_Rsq->Draw("e");
      h1->Draw("e3 same");       
      histSigmaIetaIeta_EE_Rsq->Draw("e same");
      CMS_lumi(c1,4,0);
      c1->SaveAs("Purity_EE_Rsq.png");
      c1->SaveAs("Purity_EE_Rsq.pdf");
      c1->SaveAs("Purity_EE_Rsq.root");
           
      
      char RfracFit[100];
      sprintf(RfracFit,"#splitline{%4.2f<R^{2}<%4.2f}{Purity in EE : %4.3f +- %4.3f}", RsqBins[ii], RsqBins[ii+1], PurityFromFit,SigmaErr); 
      tpav_txt->AddText(RfracFit);              
    
      float chosenEff = 0.0271;
      TLine *line =new TLine(chosenEff,0,chosenEff,100);
      line->SetLineColor(kOrange+8);
      line->SetLineWidth(2);            
      
      TCanvas *cFit=new TCanvas("cFit","cFit",850,850);
      
      cFit->cd();
      gPad->Update();
      gPad->SetLogy();
      cFit->Range(0.01,0.01,0.02,1000);
      frame4->GetYaxis()->SetTitle("Number of events");
      frame4->GetYaxis()->SetTitleOffset(1.35);
      frame4->Draw();
      tpav_txt->Draw();
      line->Draw();
      leg->Draw();
      CMS_lumi(cFit,4,0);
      cFit->SaveAs(Form("PurityFitEE_Rsqbinned_%d.png", ii));
      cFit->SaveAs(Form("PurityFitEE_Rsqbinned_%d.pdf", ii));
    }   



  ///////////////////////
  ///////////////////////
  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("PhotonControlRegionPlots"+Label+".root").c_str(), "RECREATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histSigmaIetaIeta_EB[i], Form("histSigmaIetaIeta_EB_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histSigmaIetaIeta_EE[i], Form("histSigmaIetaIeta_EE_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histSigmaIetaIetaTemplate_EB[i], Form("histSigmaIetaIetaTemplate_EB_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histSigmaIetaIetaTemplate_EE[i], Form("histSigmaIetaIetaTemplate_EE_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histSigmaIetaIeta_EB_MRRsq, "histSigmaIetaIeta_EB_MRRsq", "WriteDelete");
    file->WriteTObject(histSigmaIetaIeta_EE_MRRsq, "histSigmaIetaIeta_EE_MRRsq", "WriteDelete");
    file->WriteTObject(histSigmaIetaIeta_EB_MR, "histSigmaIetaIeta_EB_MR", "WriteDelete");
    file->WriteTObject(histSigmaIetaIeta_EB_Rsq, "histSigmaIetaIeta_EB_Rsq", "WriteDelete");
    file->WriteTObject(histSigmaIetaIeta_EE_MR, "histSigmaIetaIeta_EE_MR", "WriteDelete");
    file->WriteTObject(histSigmaIetaIeta_EE_Rsq, "histSigmaIetaIeta_EE_Rsq", "WriteDelete");

    for(int ii=0; ii<NMRBins; ii++)
      for(int jj=0; jj<NRsqBins; jj++) {
	file->WriteTObject(histSigmaIetaIeta_EB_binned[ii][jj], Form("histSigmaIetaIeta_EB_binned_%s_%d_%d",processLabels[i].c_str(), ii, jj), "WriteDelete");
	file->WriteTObject(histSigmaIetaIeta_EE_binned[ii][jj], Form("histSigmaIetaIeta_EE_binned_%s_%d_%d",processLabels[i].c_str(), ii, jj), "WriteDelete");
      }
    for(int ii=0; ii<NMRBins; ii++) {
      file->WriteTObject(histSigmaIetaIeta_EB_MRbinned[ii], Form("histSigmaIetaIeta_EB_MRbinned_%s_%d",processLabels[i].c_str(), ii), "WriteDelete");
      file->WriteTObject(histSigmaIetaIeta_EE_MRbinned[ii], Form("histSigmaIetaIeta_EE_MRbinned_%s_%d",processLabels[i].c_str(), ii), "WriteDelete");
    }
    for(int ii=0; ii<NRsqBins; ii++) {
      file->WriteTObject(histSigmaIetaIeta_EB_Rsqbinned[ii], Form("histSigmaIetaIeta_EB_Rsqbinned_%s_%d",processLabels[i].c_str(), ii), "WriteDelete");
      file->WriteTObject(histSigmaIetaIeta_EE_Rsqbinned[ii], Form("histSigmaIetaIeta_EE_Rsqbinned_%s_%d",processLabels[i].c_str(), ii), "WriteDelete");
    }
  }

  file->Close();
  delete file;
}

void PhotonControlSample_FitSigmaIEtaIEta( int option = 10) {

  vector<string> datafiles;
  vector<string> fakeTemplateFiles;
  vector<string> promptTemplateFiles;
  vector<vector<string> > bkgfiles;
  vector<string> processLabels;
  vector<int> colors;

  string datafile = "";


  //No Skims  
  if (option >= 10) {
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDCuts/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_2016_GoodLumiGolden_26p4ifb.root");
    //  
  } else {
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDCuts/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_2016_GoodLumiGolden_26p4ifb_MR300Skim.root");  
    fakeTemplateFiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDIsoCuts/Skim/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_2016_GoodLumiGolden_MR300Rsq0p15Skim.root");
    promptTemplateFiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDCuts/Skim/RunTwoRazorControlRegions_PhotonFull_GJets_HTBinned_1pb_weighted_MR300Rsq0p15Skim.root");
  }


  double lumi = 26400;

  //*********************************************************************
  //GJets Control Region
  //*********************************************************************
  if (option == 0) {
    RunPhotonControlSample_FitSigmaIEtaIEta(datafiles, fakeTemplateFiles, promptTemplateFiles, lumi,"MR300Rsq0p15",0,"MR300Rsq0p15");
  }
  if (option == 10) {
    RunPhotonControlSample_FitSigmaIEtaIEta(datafiles, fakeTemplateFiles, promptTemplateFiles, lumi,"Inclusive",0,"Inclusive");
  }


}



//**********************
//With Photon ID + Iso Cuts
//**********************
// YieldPho36_58To70 : 1.30574e+07
// YieldPho50_58To70 : 1.30185e+07
// Ratio : 0.997014

// YieldPho50_85To95 : 2.09166e+06
// YieldPho75_85To95 : 1.9931e+06
// Ratio : 0.95288

// YieldPho75_105To115 : 767340
// YieldPho90_105To115 : 768620
// Ratio : 1.00167

// YieldPho90_135To145 : 239080
// YieldPho120_135To145 : 230440
// Ratio : 0.963861

// YieldPho120_185To200 : 71070
// YieldPho165_185To200 : 70911
// Ratio : 0.997763

//**********************
//After Razor Cuts
//**********************
// YieldPho36_58To70 : 18000
// YieldPho50_58To70 : 10620
// Ratio : 0.59

// YieldPho50_85To95 : 10520
// YieldPho75_85To95 : 9120
// Ratio : 0.86692

// YieldPho75_105To115 : 11520
// YieldPho90_105To115 : 12480
// Ratio : 1.08333

// YieldPho90_135To145 : 7990
// YieldPho120_135To145 : 7545
// Ratio : 0.944305

// YieldPho120_185To200 : 7230
// YieldPho165_185To200 : 6711
// Ratio : 0.928216


/* The template for prompts in barrel was made with this snippet on command line
   TChain * chain = new TChain("ControlSampleEvent");
   chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
   chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
   chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
   chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
   chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
   double sigmaietaieta_bins_EB[100] = {0.};
   double sigmaietaieta_bins_EE[100] = {0.};

   sigmaietaieta_bins_EB[21] = 0.0103;
   sigmaietaieta_bins_EE[55] = 0.0271;
   const int NBinsSigmaietaieta = sizeof(sigmaietaieta_bins_EB)/sizeof(double)-1;
 
   for(int a = 1; a < 100; a++) { if(a!=21) sigmaietaieta_bins_EB[a] = sigmaietaieta_bins_EB[a-1] + 0.0005; }
   for(int a = 1; a < 100; a++) { if(a!=55) sigmaietaieta_bins_EE[a] = sigmaietaieta_bins_EE[a-1] + 0.0005; }


   h1 = new TH1F("h1","h1",NBinsSigmaietaieta,sigmaietaieta_bins_EB);
   h1->Sumw2();

   chain->Draw("pho1_sigmaietaieta>>h1","weight*((HLTDecision[88]==1 || HLTDecision[89]==1 || HLTDecision[90]==1 || HLTDecision[91] ==1 ||HLTDecision[92]==1 || HLTDecision[93]==1)&&fabs(pho1.Eta())<1.479&&pho1_chargediso<2.5&&pho1_sigmaietaieta<0.015)");
   h1->SaveAs("a.root")


   The template for prompts in endcap was made with this snippet on command line
   TChain * chain = new TChain("ControlSampleEvent");
   chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
   chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
   chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
   chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
   chain->Add("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/PhotonFull_1p23_2015Final_HEEleVetoCut/RunTwoRazorControlRegions_PhotonFull_GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
   h1 = new TH1F("h1","h1",NBinsSigmaietaieta,sigmaietaieta_bins_EE);
   h1->Sumw2();

   chain->Draw("pho1_sigmaietaieta>>h1","weight*((HLTDecision[88]==1 || HLTDecision[89]==1 || HLTDecision[90]==1 || HLTDecision[91] ==1 ||HLTDecision[92]==1 || HLTDecision[93]==1)&&fabs(pho1.Eta())>1.479&&pho1_chargediso<2.5&&pho1_sigmaietaieta>0.015)");
   h1->SaveAs("a.root")

*/
