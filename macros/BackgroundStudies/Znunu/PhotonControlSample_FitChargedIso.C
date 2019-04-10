

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


//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1D *hist) {
  Double_t norm = 0;
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

}



void RunPhotonControlSample_FitChargedIso(  vector<string> datafiles,  vector<string> fakeTemplateFiles,  vector<string> promptTemplateFiles,  double lumi, string option, int channelOption = -1, string label = "") {
  
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
  TH1F *promptHist = (TH1F*)PromptFile->Get("PhotonChargedIsoTemplate_Prompt_Barrel"); 
  TH1F *promptHist_EE = (TH1F*)PromptFile->Get("PhotonChargedIsoTemplate_Prompt_Endcap");
  TH1F *fakeHist_EB = (TH1F*)PromptFile->Get("PhotonChargedIsoTemplate_Fake_Barrel"); 
  TH1F *fakeHist_EE = (TH1F*)PromptFile->Get("PhotonChargedIsoTemplate_Fake_Endcap");

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  double MRBins[] = {300, 400, 500, 700, 900, 1200, 4000};
  double RsqBins[] = {0.15, 0.2, 0.25, 0.30, 0.41, 0.52, 1.5};

  const int NMRBins = sizeof(MRBins)/sizeof(double)-1;
  const int NRsqBins = sizeof(RsqBins)/sizeof(double)-1;

  double chargediso_bins_EB[100] = {0.};
  double chargediso_bins_EE[100] = {0.};

  const int NBinsChargediso = sizeof(chargediso_bins_EB)/sizeof(double)-1;

  for(int a = 1; a < 100; a++) {
    chargediso_bins_EB[a] = chargediso_bins_EB[a-1] + 0.1;
    chargediso_bins_EE[a] = chargediso_bins_EE[a-1] + 0.1;
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
  

  vector<TH1D*> histChargedIso_EB;
  vector<TH1D*> histChargedIso_EE;
  vector<TH1D*> histChargedIsoTemplate_EB;
  vector<TH1D*> histChargedIsoTemplate_EE;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    histChargedIso_EB.push_back(new TH1D(Form("histChargedIso_EB_%s",processLabels[i].c_str()), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB));
    histChargedIso_EE.push_back(new TH1D(Form("histChargedIso_EE_%s",processLabels[i].c_str()), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EE));
    histChargedIso_EB[i]->Sumw2();
    histChargedIso_EE[i]->Sumw2();

    histChargedIsoTemplate_EB.push_back(new TH1D(Form("histChargedIsoTemplate_EB_%s",processLabels[i].c_str()), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB));
    histChargedIsoTemplate_EE.push_back(new TH1D(Form("histChargedIsoTemplate_EE_%s",processLabels[i].c_str()), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EE));
    histChargedIsoTemplate_EB[i]->Sumw2();
    histChargedIsoTemplate_EE[i]->Sumw2();
  }

  TH1D *histChargedIsoFakeTemplate_EB = new TH1D ( "histChargedIsoFakeTemplate_EB", "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
  TH1D *histChargedIsoFakeTemplate_EE = new TH1D ( "histChargedIsoFakeTemplate_EE", "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EE);

  TH1D *histChargedIsoFakeTemplate_EB_binned[NMRBins][NRsqBins];
  TH1D *histChargedIsoFakeTemplate_EE_binned[NMRBins][NRsqBins];

  TH1D *histChargedIsoFakeTemplate_EB_MRbinned[NMRBins];
  TH1D *histChargedIsoFakeTemplate_EE_MRbinned[NMRBins];

  TH1D *histChargedIsoFakeTemplate_EB_Rsqbinned[NRsqBins];
  TH1D *histChargedIsoFakeTemplate_EE_Rsqbinned[NRsqBins];

  for (int a=0;a<NMRBins;a++) {
    histChargedIsoFakeTemplate_EB_MRbinned[a]  = new TH1D (Form("histChargedIsoFakeTemplate_EB_MRbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
    histChargedIsoFakeTemplate_EE_MRbinned[a]  = new TH1D (Form("histChargedIsoFakeTemplate_EE_MRbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EE);
  }

  for (int a=0;a<NRsqBins;a++) {
    histChargedIsoFakeTemplate_EB_Rsqbinned[a] = new TH1D (Form("histChargedIsoFakeTemplate_EB_Rsqbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
    histChargedIsoFakeTemplate_EE_Rsqbinned[a] = new TH1D (Form("histChargedIsoFakeTemplate_EE_Rsqbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EE);
  }

  for (int a=0;a<NMRBins;a++) {
    for (int b=0;b<NRsqBins;b++){
      histChargedIsoFakeTemplate_EB_binned[a][b] = new TH1D (Form("histChargedIsoFakeTemplate_EB_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
      histChargedIsoFakeTemplate_EB_binned[a][b] -> Sumw2();

      histChargedIsoFakeTemplate_EE_binned[a][b] = new TH1D (Form("histChargedIsoFakeTemplate_EE_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
      histChargedIsoFakeTemplate_EE_binned[a][b] -> Sumw2();
    }
  }

  TH1D *histChargedIsoPromptTemplate_EB = new TH1D ( "histChargedIsoPromptTemplate_EB", "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
  TH1D *histChargedIsoPromptTemplate_EE = new TH1D ( "histChargedIsoPromptTemplate_EE", "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EE);

  TH1D *histChargedIsoPromptTemplate_EB_binned[NMRBins][NRsqBins];
  TH1D *histChargedIsoPromptTemplate_EE_binned[NMRBins][NRsqBins];

  TH1D *histChargedIsoPromptTemplate_EB_MRbinned[NMRBins];
  TH1D *histChargedIsoPromptTemplate_EE_MRbinned[NMRBins];

  TH1D *histChargedIsoPromptTemplate_EB_Rsqbinned[NRsqBins];
  TH1D *histChargedIsoPromptTemplate_EE_Rsqbinned[NRsqBins];

  for (int a=0;a<NMRBins;a++) {
    histChargedIsoPromptTemplate_EB_MRbinned[a]  = new TH1D (Form("histChargedIsoPromptTemplate_EB_MRbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
    histChargedIsoPromptTemplate_EE_MRbinned[a]  = new TH1D (Form("histChargedIsoPromptTemplate_EE_MRbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EE);
  }

  for (int a=0;a<NRsqBins;a++) {
    histChargedIsoPromptTemplate_EB_Rsqbinned[a] = new TH1D (Form("histChargedIsoPromptTemplate_EB_Rsqbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
    histChargedIsoPromptTemplate_EE_Rsqbinned[a] = new TH1D (Form("histChargedIsoPromptTemplate_EE_Rsqbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EE);
  }

  for (int a=0;a<NMRBins;a++) {
    for (int b=0;b<NRsqBins;b++){
      histChargedIsoPromptTemplate_EB_binned[a][b] = new TH1D (Form("histChargedIsoPromptTemplate_EB_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
      histChargedIsoPromptTemplate_EB_binned[a][b] -> Sumw2();

      histChargedIsoPromptTemplate_EE_binned[a][b] = new TH1D (Form("histChargedIsoPromptTemplate_EE_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
      histChargedIsoPromptTemplate_EE_binned[a][b] -> Sumw2();
    }
  }



  TH1D *histChargedIso_EB_binned[NMRBins][NRsqBins];
  TH1D *histChargedIso_EE_binned[NMRBins][NRsqBins];

  TH1D *histChargedIso_EB_MRbinned[NMRBins];
  TH1D *histChargedIso_EE_MRbinned[NMRBins];

  TH1D *histChargedIso_EB_Rsqbinned[NRsqBins];
  TH1D *histChargedIso_EE_Rsqbinned[NRsqBins];

  for (int a=0;a<NMRBins;a++) {
    histChargedIso_EB_MRbinned[a]  = new TH1D (Form("histChargedIso_EB_MRbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
    histChargedIso_EE_MRbinned[a]  = new TH1D (Form("histChargedIso_EE_MRbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EE);
  }

  for (int a=0;a<NRsqBins;a++) {
    histChargedIso_EB_Rsqbinned[a] = new TH1D (Form("histChargedIso_EB_Rsqbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
    histChargedIso_EE_Rsqbinned[a] = new TH1D (Form("histChargedIso_EE_Rsqbinned_%d", a), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EE);
  }

  for (int a=0;a<NMRBins;a++) {
    for (int b=0;b<NRsqBins;b++){
      histChargedIso_EB_binned[a][b] = new TH1D (Form("histChargedIso_EB_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
      histChargedIso_EB_binned[a][b] -> Sumw2();

      histChargedIso_EE_binned[a][b] = new TH1D (Form("histChargedIso_EE_binned_%d_%d", a, b), "; Sigma_ietaieta; Number of Events", NBinsChargediso, chargediso_bins_EB);
      histChargedIso_EE_binned[a][b] -> Sumw2();
    }
  }

  TH2D *histChargedIso_EB_MRRsq= new TH2D("histChargedIso_EB_MRRsq", "; MR; Rsq", NMRBins, MRBins, NRsqBins, RsqBins);
  TH2D *histChargedIso_EE_MRRsq= new TH2D("histChargedIso_EE_MRRsq", "; MR; Rsq", NMRBins, MRBins, NRsqBins, RsqBins);
  TH1D *histChargedIso_EB_MR= new TH1D("histChargedIso_EB_MR", "; MR", NMRBins, MRBins);
  TH1D *histChargedIso_EB_Rsq= new TH1D("histChargedIso_EB_Rsq", "; Rsq", NRsqBins, RsqBins);
  TH1D *histChargedIso_EE_MR= new TH1D("histChargedIso_EE_MR", "; MR", NMRBins, MRBins);
  TH1D *histChargedIso_EE_Rsq= new TH1D("histChargedIso_EE_Rsq", "; Rsq", NRsqBins, RsqBins);

  histChargedIso_EB_MR->Sumw2();
  histChargedIso_EB_Rsq->Sumw2();
  histChargedIso_EE_MR->Sumw2();
  histChargedIso_EE_Rsq->Sumw2();


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

	// if ( processLabels[i] == "Data" || processLabels[i] == "FakeTemplate" ) {	  
	if ( processLabels[i] == "Data" ) {	  
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
	if (! (events->pho1_hOverE < 0.05)) continue;

	//apply cuts on ecal iso & neutral hadron iso
	if (fabs(events->pho1.Eta()) < 1.5) {
	  if (!( events->pho1_photoniso < 2.554 + 0.0047*events->pho1.Pt()
		 && events->pho1_neutralhadroniso < 4.50 + 0.0148*events->pho1.Pt() + 0.000017*events->pho1.Pt()*events->pho1.Pt()
		 )) continue;
	} else {
	  if (!(events->pho1_photoniso < 2.75 + 0.0034*events->pho1.Pt()
		&& events->pho1_neutralhadroniso < 0.432 + 0.0163*events->pho1.Pt() + 0.000014*events->pho1.Pt()*events->pho1.Pt()
		)) continue;
	}


	if ( processLabels[i] == "Data") {

	  // select photons in EB
	  if(fabs(events->pho1.Eta()) < 1.479 && events->pho1_sigmaietaieta < 0.01042 && events->pho1_chargediso < 10 ) {
	    histChargedIso_EB[i]->Fill(events->pho1_chargediso);
	      
	    for(int ii = 0; ii<NMRBins; ii++) {
	      for(int jj = 0; jj<NRsqBins; jj++) {
		if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] ) {
		  if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
		    histChargedIso_EB_binned[ii][jj]->Fill(events->pho1_chargediso, weight);
		  }
		}
	      }
	    }
	    
	    for(int ii = 0; ii<NMRBins; ii++) {
	      if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] && events->Rsq_NoPho > 0.15 ) {
		histChargedIso_EB_MRbinned[ii]->Fill(events->pho1_chargediso, weight);
	      }
	    }
	    for(int jj = 0; jj<NRsqBins; jj++) {
	      if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] && events->MR_NoPho > 400. ){
		histChargedIso_EB_Rsqbinned[jj]->Fill(events->pho1_chargediso, weight);
	      }					      
	    }
	  }
	  // select photons in EE
	  if(fabs(events->pho1.Eta()) > 1.479 && events->pho1_sigmaietaieta < 0.02649 && events->pho1_chargediso < 10) {
	    histChargedIso_EE[i]->Fill(events->pho1_chargediso);
	    
	    for(int ii = 0; ii<NMRBins; ii++) {
	      for(int jj = 0; jj<NRsqBins; jj++) {
		if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] ) {
		  if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
		    histChargedIso_EE_binned[ii][jj]->Fill(events->pho1_chargediso, weight);
		  }
		}
	      }
	    }
	    
	    for(int ii = 0; ii<NMRBins; ii++) {
	      if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] && events->Rsq_NoPho > 0.15 ) {
		histChargedIso_EE_MRbinned[ii]->Fill(events->pho1_chargediso, weight);
	      }
	    }
	    for(int jj = 0; jj<NRsqBins; jj++) {
	      if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] && events->MR_NoPho > 400. ){
		histChargedIso_EE_Rsqbinned[jj]->Fill(events->pho1_chargediso, weight);
	      }					      
	    }
	  }	  
	} // end if Data

	if ( processLabels[i] == "FakeTemplate") {
	  if( fabs(events->pho1.Eta()) < 1.479 ) { // make the fake template for barrel photons
	    if( events->pho1_sigmaietaieta > 0.011 && events->pho1_sigmaietaieta < 0.015 && events->pho1_chargediso < 10 ) {
	      histChargedIsoTemplate_EB[i]->Fill(events->pho1_chargediso);
	      histChargedIsoFakeTemplate_EB->Fill(events->pho1_chargediso, weight);

	      for(int ii = 0; ii<NMRBins; ii++) {
		for(int jj = 0; jj<NRsqBins; jj++) {
		  if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] ) {
		    if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
		      histChargedIsoFakeTemplate_EB_binned[ii][jj]->Fill(events->pho1_chargediso, weight);
		    }
		  }
		}
	      }
	      
	      for(int ii = 0; ii<NMRBins; ii++) {
		if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] && events->Rsq_NoPho > 0.15 ) {
		  histChargedIsoFakeTemplate_EB_MRbinned[ii]->Fill(events->pho1_chargediso, weight);
		}
	      }
	      for(int jj = 0; jj<NRsqBins; jj++) {
		if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] && events->MR_NoPho > 400. ){
		  histChargedIsoFakeTemplate_EB_Rsqbinned[jj]->Fill(events->pho1_chargediso, weight);
		}					      
	      }
	    }
	  }
	    
	  if( fabs(events->pho1.Eta()) > 1.479 ) { // make the fake template for endcap photons
	    if( events->pho1_sigmaietaieta > 0.03 && events->pho1_sigmaietaieta < 0.035 && events->pho1_chargediso < 10 ) {
	      histChargedIsoTemplate_EE[i]->Fill(events->pho1_chargediso);
	      histChargedIsoFakeTemplate_EE->Fill(events->pho1_chargediso, weight);

	      for(int ii = 0; ii<NMRBins; ii++) {
		for(int jj = 0; jj<NRsqBins; jj++) {
		  if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] ) {
		    if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
		      histChargedIsoFakeTemplate_EE_binned[ii][jj]->Fill(events->pho1_chargediso, weight);
		    }
		  }
		}
	      }
	      
	      for(int ii = 0; ii<NMRBins; ii++) {
		if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] && events->Rsq_NoPho > 0.15 ) {
		  histChargedIsoFakeTemplate_EE_MRbinned[ii]->Fill(events->pho1_chargediso, weight);
		}
	      }
	      for(int jj = 0; jj<NRsqBins; jj++) {
		if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] && events->MR_NoPho > 400. ){
		  histChargedIsoFakeTemplate_EE_Rsqbinned[jj]->Fill(events->pho1_chargediso, weight);
		}
	      }
	    }
	  }
	} // end if Fake Template

	if ( processLabels[i] == "PromptTemplate") {
	  if( fabs(events->pho1.Eta()) < 1.479 ) { // make the prompt template for barrel photons
	    if( events->pho1_sigmaietaieta < 0.0103 && events->pho1_chargediso < 10 ) {
	      histChargedIsoTemplate_EB[i]->Fill(events->pho1_chargediso);
	      histChargedIsoPromptTemplate_EB->Fill(events->pho1_chargediso, weight);

	      for(int ii = 0; ii<NMRBins; ii++) {
		for(int jj = 0; jj<NRsqBins; jj++) {
		  if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] ) {
		    if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
		      histChargedIsoPromptTemplate_EB_binned[ii][jj]->Fill(events->pho1_chargediso, weight);
		    }
		  }
		}
	      }
	      
	      for(int ii = 0; ii<NMRBins; ii++) {
		if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] && events->Rsq_NoPho > 0.15 ) {
		  histChargedIsoPromptTemplate_EB_MRbinned[ii]->Fill(events->pho1_chargediso, weight);
		}
	      }
	      for(int jj = 0; jj<NRsqBins; jj++) {
		if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] && events->MR_NoPho > 400. ){
		  histChargedIsoPromptTemplate_EB_Rsqbinned[jj]->Fill(events->pho1_chargediso, weight);
		}					      
	      }
	    }
	  }
	    
	  if( fabs(events->pho1.Eta()) > 1.479 ) { // make the prompt template for endcap photons
	    if( events->pho1_sigmaietaieta < 0.0271 && events->pho1_chargediso < 10 ) {
	      histChargedIsoTemplate_EE[i]->Fill(events->pho1_chargediso);
	      histChargedIsoPromptTemplate_EE->Fill(events->pho1_chargediso, weight);

	      for(int ii = 0; ii<NMRBins; ii++) {
		for(int jj = 0; jj<NRsqBins; jj++) {
		  if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] ) {
		    if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] ){
		      histChargedIsoPromptTemplate_EE_binned[ii][jj]->Fill(events->pho1_chargediso, weight);
		    }
		  }
		}
	      }
	      
	      for(int ii = 0; ii<NMRBins; ii++) {
		if(events->MR_NoPho > MRBins[ii] && events->MR_NoPho < MRBins[ii+1] && events->Rsq_NoPho > 0.15 ) {
		  histChargedIsoPromptTemplate_EE_MRbinned[ii]->Fill(events->pho1_chargediso, weight);
		}
	      }
	      for(int jj = 0; jj<NRsqBins; jj++) {
		if(events->Rsq_NoPho > RsqBins[jj] && events->Rsq_NoPho < RsqBins[jj+1] && events->MR_NoPho > 400. ){
		  histChargedIsoPromptTemplate_EE_Rsqbinned[jj]->Fill(events->pho1_chargediso, weight);
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

  double rangeMin=0.0;
  double rangeMax=10.0;




  // // now do the fits in the MR/Rsq bins   
  // for(int ii=0; ii<NMRBins; ii++)
  //   for(int jj=0; jj<NRsqBins; jj++) {
  //     RooRealVar varChargedIso("varChargedIso","#sigma_{i#etai#eta}",rangeMin,rangeMax);
  //     RooDataHist dataSR("dataSR","dataSR",varChargedIso,Import(*histChargedIso_EB_binned[ii][jj]));
   
  //     RooDataHist MCprompt("MCprompt","MCprompt",varChargedIso,Import(*histChargedIsoPromptTemplate_EB));
  //     RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoFakeTemplate_EB));
  //     // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoTemplate_EB[0]));
  //     // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoFakeTemplate_EB_binned[ii][jj]));

  //     RooHistPdf PromptPDF("PromptPDF","PromptPDF",varChargedIso,MCprompt,0);
  //     RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",varChargedIso,MCnonPrompt,0);
   
  //     RooRealVar Fitfrac("Fitfrac","Fitfrac",0.95,0.50,1.0);
   
  //     RooAddPdf PDF("PDF","PSR-NPSB-PDF", RooArgList(PromptPDF,nonPromptPDF),Fitfrac);
  
  //     TPaveText *tpav_txt = new TPaveText(0.5,0.7,0.85,0.9,"brNDC");
   
  //     tpav_txt->SetBorderSize(0);
  //     tpav_txt->SetFillStyle(0);
  //     tpav_txt->SetTextAlign(11);
  //     tpav_txt->SetTextFont(42);
  //     tpav_txt->SetTextSize(0.03);
  //     // double pureMCPurity=hpromptSR_Scut->Integral()/hpromptPlusnonPromptSR_Scut->Integral();
   
   
  //     // char mcPurity[100];
  //     // sprintf(mcPurity,"Purity(MC, EB): %4.3f ",pureMCPurity); 
  //     // tpav_txt->AddText(mcPurity);
   
   
  //     RooPlot* frame4 = varChargedIso.frame(Title(""));
  //     frame4->SetTitle("");

  //     dataSR.plotOn(frame4,Name("dataSR1"));  
  //     PDF.fitTo(dataSR);
  //     PDF.plotOn(frame4);      
  
  //     PDF.plotOn(frame4,Name("PromptPDF"),Components(RooArgList(PromptPDF)), LineColor(kRed),LineStyle(2), LineWidth(2));
  //     PDF.plotOn(frame4,Name("nonPromptPDF"),Components(RooArgList(nonPromptPDF)), LineColor(kGreen),LineStyle(2), LineWidth(3));
  //     PDF.plotOn(frame4,Name("PDF"),Components(RooArgList(PDF)), LineColor(kBlue),LineStyle(1), LineWidth(2));   

  //     TLegend *leg=new TLegend(0.6,0.55,0.85,0.75,NULL,"brNDC");
  //     leg->SetTextFont(62);
  //     leg->SetLineColor(1);
  //     leg->SetLineStyle(1);
  //     leg->SetLineWidth(3);
  //     leg->SetFillColor(0);
  //     leg->SetFillStyle(1001);
  //     leg->SetShadowColor(0);
  //     leg->SetDrawOption(0);
  //     leg->SetBorderSize(0);
  //     leg->SetTextSize(0.02);
  
  //     leg->AddEntry(frame4->findObject("PromptPDF"),"Prompt Template, EB","l");
  //     leg->AddEntry(frame4->findObject("nonPromptPDF"),"Fak Template, EB","l");
  //     leg->AddEntry(frame4->findObject("PDF"),"Total Fit, EB","l");
  //     leg->AddEntry(frame4->findObject("dataSR1"),"Data, EB","P");      

  //     varChargedIso.setRange("SR",rangeMin,2.5);
  //     RooAbsReal* fracP = PromptPDF.createIntegral(varChargedIso,NormSet(varChargedIso),Range("SR"));
  //     RooAbsReal* Ip=PromptPDF.createIntegral(varChargedIso);
  //     RooAbsReal* fracNP = nonPromptPDF.createIntegral(varChargedIso,NormSet(varChargedIso),Range("SR"));
  //     RooAbsReal* Inp=nonPromptPDF.createIntegral(varChargedIso);

  //     double f=Fitfrac.getVal();
  //     double fError=Fitfrac.getError();
  //     double dataAreaFull=histChargedIso_EB_binned[ii][jj]->Integral();
  //     double promptIntFull=f*dataAreaFull;
  //     double nonpromptIntFull=(1-f)*dataAreaFull;

  //     double fracP_val=fracP->getVal();
  
  //     double fracNP_val=fracNP->getVal();
  
  //     double promptIntS=fracP_val*f;
  //     double nonpromptIntS=fracNP_val*(1-f);
  //     double PurityFromFit=promptIntS/(promptIntS+nonpromptIntS);
  
  //     float SigmaErr=(fError/f)*PurityFromFit;

  //     cout<<"Beta error: "<<SigmaErr<<endl;

  //     cout<<"fracP:  "<<fracP_val<<endl;
  //     cout<<"fracNP:  "<<fracNP_val<<endl;
  //     cout<<"Purity From Fit: "<<PurityFromFit<<endl;
  
  //     histChargedIso_EB_MRRsq->SetBinContent(ii+1, jj+1, PurityFromFit);

  //     char RfracFit[100];
  //     sprintf(RfracFit,"Purity in EB: %4.3f +- %4.3f ",PurityFromFit,SigmaErr); 
  //     tpav_txt->AddText(RfracFit);

  //     float chosenEff = 0.0103;
  //     TLine *line =new TLine(chosenEff,0,chosenEff,100);
  //     line->SetLineColor(kOrange+8);
  //     line->SetLineWidth(2);    

  //     TCanvas *cFit=new TCanvas("cFit","cFit",850,850);
   
  //     cFit->cd();
  //     gPad->Update();
  //     gPad->SetLogy();
  //     cFit->Range(0.01,0.01,0.02,1000);
  //     frame4->GetYaxis()->SetTitle("Number of events");
  //     frame4->GetYaxis()->SetTitleOffset(1.35);
  //     frame4->Draw();
  //     tpav_txt->Draw();
  //     line->Draw();
  //     leg->Draw();
  //     CMS_lumi(cFit,4,0);
  //     cFit->SaveAs(Form("PurityFitEB_%d_%d.png", ii, jj));
  //     cFit->SaveAs(Form("PurityFitEB_%d_%d.pdf", ii, jj));
  //   }
  
  //////////////////////////
  // now do the fits in the MR bins   
  for(int ii=0; ii<NMRBins; ii++)
  // for(int ii=0; ii<1; ii++)
    {     
      cout<<"MR BINNED FIT: "<<ii<<endl;

      RooRealVar varChargedIso("varChargedIso","#sigma_{i#etai#eta}",rangeMin,rangeMax);
      RooDataHist dataSR("dataSR","dataSR",varChargedIso,Import(*histChargedIso_EB_MRbinned[ii]));
      
      RooDataHist MCprompt("MCprompt","MCprompt",varChargedIso,Import(*histChargedIsoPromptTemplate_EB));
      RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoFakeTemplate_EB));
     // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoTemplate_EB[0]));
       // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoFakeTemplate_EB_MRbinned[ii]));
   
      RooHistPdf PromptPDF("PromptPDF","PromptPDF",varChargedIso,MCprompt,0);
      RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",varChargedIso,MCnonPrompt,0);
      
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
   
   
      RooPlot* frame4 = varChargedIso.frame(Title(""));
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

      varChargedIso.setRange("SR",rangeMin,1.3);
      RooAbsReal* fracP = PromptPDF.createIntegral(varChargedIso,NormSet(varChargedIso),Range("SR"));
      RooAbsReal* Ip=PromptPDF.createIntegral(varChargedIso);
      RooAbsReal* fracNP = nonPromptPDF.createIntegral(varChargedIso,NormSet(varChargedIso),Range("SR"));
      RooAbsReal* Inp=nonPromptPDF.createIntegral(varChargedIso);
      
      double f=Fitfrac.getVal();
      double fError=Fitfrac.getError();
      double dataAreaFull=histChargedIso_EB_MRbinned[ii]->Integral();
      double promptIntFull=f*dataAreaFull;
      double nonpromptIntFull=(1-f)*dataAreaFull;
      
      double fracP_val=fracP->getVal();
      
      double fracNP_val=fracNP->getVal();
      
      double promptIntS=fracP_val*f;
      double nonpromptIntS=fracNP_val*(1-f);
      double PurityFromFit=promptIntS/(promptIntS+nonpromptIntS);
      
      float SigmaErr=(fError/f)*PurityFromFit;
      
      cout << "fitFrac = " << f << "\n";
      cout<<"Beta error: "<<SigmaErr<<endl;            
      cout<<"fracP:  "<<fracP_val<<endl;
      cout<<"fracNP:  "<<fracNP_val<<endl;
      cout<<"Purity From Fit: " << fracP_val << " * " << f << " / ( " << fracNP_val << " * " << (1-f) << " + " << fracP_val << " * " << f << " ) = "
	  <<PurityFromFit<<endl;
      
      
      histChargedIso_EB_MR->SetBinContent(ii+1, PurityFromFit);
      histChargedIso_EB_MR->SetBinError(ii+1, SigmaErr);
      
      TCanvas *c1=new TCanvas("c1","c1",850,850);      

      TH1F *h1 = new TH1F("h1","h1",1000,0, 4000);
      for(int a = 0; a < 1000; a++) { h1->SetBinContent(a, 0.95); h1->SetBinError(a, 0.05); }
      h1->SetFillStyle(3356);
      h1->SetFillColor(4);
      histChargedIso_EB_MR->SetBinContent(ii+1, PurityFromFit);
      histChargedIso_EB_MR->SetBinError(ii+1, SigmaErr);
      histChargedIso_EB_MR->SetMarkerStyle(8);
      histChargedIso_EB_MR->SetStats(0);
      histChargedIso_EB_MR->SetLineColor(2);
      histChargedIso_EB_MR->SetLineWidth(2);
      histChargedIso_EB_MR->SetTitle("");
      histChargedIso_EB_MR->GetXaxis()->SetTitle("MR [GeV]");
      histChargedIso_EB_MR->GetYaxis()->SetTitle("Photon Purity");      
      histChargedIso_EB_MR->GetYaxis()->SetTitleOffset(1.35);
      histChargedIso_EB_MR->GetYaxis()->SetRangeUser(0.75, 1.05);
      histChargedIso_EB_MR->GetXaxis()->SetRangeUser(400, 1200.);
      histChargedIso_EB_MR->Draw("e");
      h1->Draw("e3 same");       
      histChargedIso_EB_MR->Draw("e same");

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
      
      RooRealVar varChargedIso("varChargedIso","#sigma_{i#etai#eta}",rangeMin,rangeMax);
      RooDataHist dataSR("dataSR","dataSR",varChargedIso,Import(*histChargedIso_EB_Rsqbinned[ii]));
      
      RooDataHist MCprompt("MCprompt","MCprompt",varChargedIso,Import(*histChargedIsoPromptTemplate_EB));
      RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoFakeTemplate_EB));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoTemplate_EB[0]));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoFakeTemplate_EB_Rsqbinned[ii]));
     
      RooHistPdf PromptPDF("PromptPDF","PromptPDF",varChargedIso,MCprompt,0);
      RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",varChargedIso,MCnonPrompt,0);
      
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
      
      
      RooPlot* frame4 = varChargedIso.frame(Title(""));
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
      
      varChargedIso.setRange("SR",rangeMin,1.3);
      RooAbsReal* fracP = PromptPDF.createIntegral(varChargedIso,NormSet(varChargedIso),Range("SR"));
      RooAbsReal* Ip=PromptPDF.createIntegral(varChargedIso);
      RooAbsReal* fracNP = nonPromptPDF.createIntegral(varChargedIso,NormSet(varChargedIso),Range("SR"));
      RooAbsReal* Inp=nonPromptPDF.createIntegral(varChargedIso);
      
      double f=Fitfrac.getVal();
      double fError=Fitfrac.getError();
      double dataAreaFull=histChargedIso_EB_Rsqbinned[ii]->Integral();
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
            
      histChargedIso_EB_Rsq->SetBinContent(ii+1, PurityFromFit);
      histChargedIso_EB_Rsq->SetBinError(ii+1, SigmaErr);

      TCanvas *c1=new TCanvas("c1","c1",850,850);      
      c1->cd();
      gPad->Update();

      TH1F *h1 = new TH1F("h1","h1",1000,0, 1.5);
      for(int a = 0; a < 1000; a++) { h1->SetBinContent(a, 0.95); h1->SetBinError(a, 0.05); }
      h1->SetFillStyle(3356);
      h1->SetFillColor(4);
      histChargedIso_EB_Rsq->SetBinContent(ii+1, PurityFromFit);
      histChargedIso_EB_Rsq->SetBinError(ii+1, SigmaErr);
      histChargedIso_EB_Rsq->SetMarkerStyle(8);
      histChargedIso_EB_Rsq->SetStats(0);
      histChargedIso_EB_Rsq->SetLineColor(2);
      histChargedIso_EB_Rsq->SetLineWidth(2);
      histChargedIso_EB_Rsq->SetTitle("");
      histChargedIso_EB_Rsq->GetXaxis()->SetTitle("Rsq");
      histChargedIso_EB_Rsq->GetYaxis()->SetTitle("Photon Purity");      
      histChargedIso_EB_Rsq->GetYaxis()->SetRangeUser(0.75, 1.05);
      histChargedIso_EB_Rsq->GetYaxis()->SetTitleOffset(1.35);
      histChargedIso_EB_Rsq->GetXaxis()->SetRangeUser(0.15, 0.52);
      histChargedIso_EB_Rsq->Draw("e");
      h1->Draw("e3 same");       
      histChargedIso_EB_Rsq->Draw("e same");
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



  rangeMin = 0.0;
  rangeMax = 10.0;

  // //  now do the fits in the MR/Rsq bins in EE
  // for(int ii=0; ii<NMRBins; ii++)
  //   for(int jj=0; jj<NRsqBins; jj++) {
       
  //     cout<<"BINNED FIT in EE: "<<ii<<" "<<jj<<endl;
              
  //     RooRealVar varChargedIso("varChargedIso","#sigma_{i#etai#eta}",rangeMin,rangeMax);
  //     RooDataHist dataSR("dataSR","dataSR",varChargedIso,Import(*histChargedIso_EE_binned[ii][jj]));
	 
  //     RooDataHist MCprompt("MCprompt","MCprompt",varChargedIso,Import(*histChargedIsoPromptTemplate_EE));
  //     RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoFakeTemplate_EE));
  //     // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoTemplate_EE[0]));
  //     // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoFakeTemplate_EE_binned[ii][jj]));
  
  //     RooHistPdf PromptPDF("PromptPDF","PromptPDF",varChargedIso,MCprompt,0);
  //     RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",varChargedIso,MCnonPrompt,0);
   
  //     RooRealVar Fitfrac("Fitfrac","Fitfrac",0.7,0.,1.0);
   
  //     RooAddPdf PDF("PDF","PSR-NPSB-PDF", RooArgList(PromptPDF,nonPromptPDF),Fitfrac);

  //     TPaveText *tpav_txt = new TPaveText(0.5,0.7,0.85,0.9,"brNDC");
   
  //     tpav_txt->SetBorderSize(0);
  //     tpav_txt->SetFillStyle(0);
  //     tpav_txt->SetTextAlign(11);
  //     tpav_txt->SetTextFont(42);
  //     tpav_txt->SetTextSize(0.03);
  //     // double pureMCPurity=hpromptSR_Scut->Integral()/hpromptPlusnonPromptSR_Scut->Integral();
   
   
  //     // char mcPurity[100];
  //     // sprintf(mcPurity,"Purity(MC, EB): %4.3f ",pureMCPurity); 
  //     // tpav_txt->AddText(mcPurity);
   
  //     RooPlot* frame4 = varChargedIso.frame(Title(""));
  //     frame4->SetTitle("");

  //     dataSR.plotOn(frame4,Name("dataSR1"));  

  //     PDF.fitTo(dataSR);
  //     PDF.plotOn(frame4);      

  //     PDF.plotOn(frame4,Name("PromptPDF"),Components(RooArgList(PromptPDF)), LineColor(kRed),LineStyle(2), LineWidth(2));
  //     PDF.plotOn(frame4,Name("nonPromptPDF"),Components(RooArgList(nonPromptPDF)), LineColor(kGreen),LineStyle(2), LineWidth(3));
  //     PDF.plotOn(frame4,Name("PDF"),Components(RooArgList(PDF)), LineColor(kBlue),LineStyle(1), LineWidth(2));   

  //     TLegend *leg=new TLegend(0.6,0.55,0.85,0.75,NULL,"brNDC");
  //     leg->SetTextFont(62);
  //     leg->SetLineColor(1);
  //     leg->SetLineStyle(1);
  //     leg->SetLineWidth(3);
  //     leg->SetFillColor(0);
  //     leg->SetFillStyle(1001);
  //     leg->SetShadowColor(0);
  //     leg->SetDrawOption(0);
  //     leg->SetBorderSize(0);
  //     leg->SetTextSize(0.02);

  //     leg->AddEntry(frame4->findObject("PromptPDF"),"Prompt Template, EE","l");
  //     leg->AddEntry(frame4->findObject("nonPromptPDF"),"Fak Template, EE","l");
  //     leg->AddEntry(frame4->findObject("PDF"),"Total Fit, EE","l");
  //     leg->AddEntry(frame4->findObject("dataSR1"),"Data, EE","P");      

  //     varChargedIso.setRange("SR",rangeMin,2.5); 
  //     RooAbsReal* fracP = PromptPDF.createIntegral(varChargedIso,NormSet(varChargedIso),Range("SR"));
  //     RooAbsReal* Ip=PromptPDF.createIntegral(varChargedIso);
  //     RooAbsReal* fracNP = nonPromptPDF.createIntegral(varChargedIso,NormSet(varChargedIso),Range("SR"));
  //     RooAbsReal* Inp=nonPromptPDF.createIntegral(varChargedIso);

  //     double f=Fitfrac.getVal();
  //     double fError=Fitfrac.getError();
  //     double dataAreaFull=histChargedIso_EE_binned[ii][jj]->Integral();
  //     double promptIntFull=f*dataAreaFull;
  //     double nonpromptIntFull=(1-f)*dataAreaFull;

  //     double fracP_val=fracP->getVal();
  
  //     double fracNP_val=fracNP->getVal();
  
  //     double promptIntS=fracP_val*f;
  //     double nonpromptIntS=fracNP_val*(1-f);
  //     double PurityFromFit=promptIntS/(promptIntS+nonpromptIntS);
  
  //     float SigmaErr=(fError/f)*PurityFromFit;

  //     cout<<"Beta error: "<<SigmaErr<<endl;

  //     cout<<"fracP:  "<<fracP_val<<endl;
  //     cout<<"fracNP:  "<<fracNP_val<<endl;
  //     cout<<"Purity From Fit: "<<PurityFromFit<<endl;
  
  //     histChargedIso_EB_MRRsq->SetBinContent(ii+1, jj+1, PurityFromFit);

  //     char RfracFit[100];
  //     sprintf(RfracFit,"Purity in EE: %4.3f +- %4.3f ",PurityFromFit,SigmaErr); 
  //     tpav_txt->AddText(RfracFit);

  //     float chosenEff = 0.0271;
  //     TLine *line =new TLine(chosenEff,0,chosenEff,100);
  //     line->SetLineColor(kOrange+8);
  //     line->SetLineWidth(2);    

  //     TCanvas *cFit=new TCanvas("cFit","cFit",850,850);
   
  //     cFit->cd();
  //     gPad->Update();
  //     gPad->SetLogy();
  //     cFit->Range(0.01,0.01,0.02,1000);
  //     frame4->GetYaxis()->SetTitle("Number of events");
  //     frame4->GetYaxis()->SetTitleOffset(1.35);
  //     frame4->Draw();
  //     tpav_txt->Draw();
  //     line->Draw();
  //     leg->Draw();
  //     CMS_lumi(cFit,4,0);
  //     cFit->SaveAs(Form("PurityFitEE_%d_%d.png", ii, jj));
  //     cFit->SaveAs(Form("PurityFitEE_%d_%d.pdf", ii, jj));


  //   }  

    // now do the fits in the MR bins   
  for(int ii=0; ii<NMRBins; ii++)
    {     
      cout<<"MR BINNED FIT: "<<ii<<endl;
      
      RooRealVar varChargedIso("varChargedIso","#sigma_{i#etai#eta}",rangeMin,rangeMax);
      RooDataHist dataSR("dataSR","dataSR",varChargedIso,Import(*histChargedIso_EE_MRbinned[ii]));
      
      RooDataHist MCprompt("MCprompt","MCprompt",varChargedIso,Import(*histChargedIsoPromptTemplate_EE));
      RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoFakeTemplate_EE));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoTemplate_EE[0]));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoFakeTemplate_EE_MRbinned[ii]));
      
      RooHistPdf PromptPDF("PromptPDF","PromptPDF",varChargedIso,MCprompt,0);
      RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",varChargedIso,MCnonPrompt,0);
      
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
   
   
      RooPlot* frame4 = varChargedIso.frame(Title(""));
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

      varChargedIso.setRange("SR",rangeMin,0.2);
      RooAbsReal* fracP = PromptPDF.createIntegral(varChargedIso,NormSet(varChargedIso),Range("SR"));
      RooAbsReal* Ip=PromptPDF.createIntegral(varChargedIso);
      RooAbsReal* fracNP = nonPromptPDF.createIntegral(varChargedIso,NormSet(varChargedIso),Range("SR"));
      RooAbsReal* Inp=nonPromptPDF.createIntegral(varChargedIso);
      
      double f=Fitfrac.getVal();
      double fError=Fitfrac.getError();
      double dataAreaFull=histChargedIso_EE_MRbinned[ii]->Integral();
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
      histChargedIso_EE_MR->SetBinContent(ii+1, PurityFromFit);
      histChargedIso_EE_MR->SetBinError(ii+1, SigmaErr);
      histChargedIso_EE_MR->SetMarkerStyle(8);
      histChargedIso_EE_MR->SetStats(0);
      histChargedIso_EE_MR->SetLineColor(2);
      histChargedIso_EE_MR->SetLineWidth(2);
      histChargedIso_EE_MR->SetTitle("");
      histChargedIso_EE_MR->GetXaxis()->SetTitle("MR [GeV]");
      histChargedIso_EE_MR->GetYaxis()->SetTitle("Photon Purity");      
      histChargedIso_EE_MR->GetYaxis()->SetRangeUser(0.75, 1.05);
      histChargedIso_EE_MR->GetYaxis()->SetTitleOffset(1.35);
      histChargedIso_EE_MR->GetXaxis()->SetRangeUser(400, 1200);
      histChargedIso_EE_MR->Draw("e");
      h1->Draw("e3 same");       
      histChargedIso_EE_MR->Draw("e same");
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
      
      RooRealVar varChargedIso("varChargedIso","#sigma_{i#etai#eta}",rangeMin,rangeMax);
      RooDataHist dataSR("dataSR","dataSR",varChargedIso,Import(*histChargedIso_EE_Rsqbinned[ii]));
      
      RooDataHist MCprompt("MCprompt","MCprompt",varChargedIso,Import(*histChargedIsoPromptTemplate_EE));
      RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoFakeTemplate_EE));
      // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoTemplate_EE[0]));
       // RooDataHist MCnonPrompt("MCnonPrompt","MCnonPrompt",varChargedIso,Import(*histChargedIsoFakeTemplate_EE_Rsqbinned[ii]));
     
      RooHistPdf PromptPDF("PromptPDF","PromptPDF",varChargedIso,MCprompt,0);
      RooHistPdf nonPromptPDF("nonPromptPDF","nonPromptPDF",varChargedIso,MCnonPrompt,0);
      
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
   
   
      RooPlot* frame4 = varChargedIso.frame(Title(""));
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

      varChargedIso.setRange("SR",rangeMin,0.2);
      RooAbsReal* fracP = PromptPDF.createIntegral(varChargedIso,NormSet(varChargedIso),Range("SR"));
      RooAbsReal* Ip=PromptPDF.createIntegral(varChargedIso);
      RooAbsReal* fracNP = nonPromptPDF.createIntegral(varChargedIso,NormSet(varChargedIso),Range("SR"));
      RooAbsReal* Inp=nonPromptPDF.createIntegral(varChargedIso);
      
      double f=Fitfrac.getVal();
      double fError=Fitfrac.getError();
      double dataAreaFull=histChargedIso_EE_Rsqbinned[ii]->Integral();
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
      
      histChargedIso_EE_Rsq->SetBinContent(ii+1, PurityFromFit);
      histChargedIso_EE_Rsq->SetBinError(ii+1, SigmaErr);

      TCanvas *c1=new TCanvas("c1","c1",850,850);      
      c1->cd();
      gPad->Update();

      TH1F *h1 = new TH1F("h1","h1",1000,0, 0.4);
      for(int a = 0; a < 1000; a++) { h1->SetBinContent(a, 0.95); h1->SetBinError(a, 0.05); }
      h1->SetFillStyle(3356);
      h1->SetFillColor(4);
      histChargedIso_EE_Rsq->SetBinContent(ii+1, PurityFromFit);
      histChargedIso_EE_Rsq->SetBinError(ii+1, SigmaErr);
      histChargedIso_EE_Rsq->SetMarkerStyle(8);
      histChargedIso_EE_Rsq->SetStats(0);
      histChargedIso_EE_Rsq->SetLineColor(2);
      histChargedIso_EE_Rsq->SetLineWidth(2);
      histChargedIso_EE_Rsq->SetTitle("");
      histChargedIso_EE_Rsq->GetXaxis()->SetTitle("Rsq");
      histChargedIso_EE_Rsq->GetYaxis()->SetTitle("Photon Purity");      
      histChargedIso_EE_Rsq->GetYaxis()->SetRangeUser(0.75, 1.05);
      histChargedIso_EE_Rsq->GetYaxis()->SetTitleOffset(1.35);
      histChargedIso_EE_Rsq->GetXaxis()->SetRangeUser(0.15, 0.52);
      histChargedIso_EE_Rsq->Draw("e");
      h1->Draw("e3 same");       
      histChargedIso_EE_Rsq->Draw("e same");
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

  file->WriteTObject(histChargedIso_EB[0], Form("histChargedIso_EB_%s",processLabels[0].c_str()), "WriteDelete");
  file->WriteTObject(histChargedIso_EE[0], Form("histChargedIso_EE_%s",processLabels[0].c_str()), "WriteDelete");
  file->WriteTObject(histChargedIso_EB_MRRsq, "histChargedIso_EB_MRRsq", "WriteDelete");
  file->WriteTObject(histChargedIso_EE_MRRsq, "histChargedIso_EE_MRRsq", "WriteDelete");
  file->WriteTObject(histChargedIso_EB_MR, "histChargedIso_EB_MR", "WriteDelete");
  file->WriteTObject(histChargedIso_EB_Rsq, "histChargedIso_EB_Rsq", "WriteDelete");
  file->WriteTObject(histChargedIso_EE_MR, "histChargedIso_EE_MR", "WriteDelete");
  file->WriteTObject(histChargedIso_EE_Rsq, "histChargedIso_EE_Rsq", "WriteDelete");
  
  NormalizeHist(histChargedIsoFakeTemplate_EB);
  NormalizeHist(histChargedIsoFakeTemplate_EE);
  NormalizeHist(histChargedIsoPromptTemplate_EB);
  NormalizeHist(histChargedIsoPromptTemplate_EE);
  file->WriteTObject(histChargedIsoFakeTemplate_EB, "histChargedIsoFakeTemplate_EB", "WriteDelete");
  file->WriteTObject(histChargedIsoFakeTemplate_EE, "histChargedIsoFakeTemplate_EE", "WriteDelete");
  file->WriteTObject(histChargedIsoPromptTemplate_EB, "histChargedIsoPromptTemplate_EB", "WriteDelete");
  file->WriteTObject(histChargedIsoPromptTemplate_EE, "histChargedIsoPromptTemplate_EE", "WriteDelete");
        
  for(int ii=0; ii<NMRBins; ii++)
    for(int jj=0; jj<NRsqBins; jj++) {
      file->WriteTObject(histChargedIso_EB_binned[ii][jj], Form("histChargedIso_EB_binned_%s_%d_%d",processLabels[0].c_str(), ii, jj), "WriteDelete");
      file->WriteTObject(histChargedIso_EE_binned[ii][jj], Form("histChargedIso_EE_binned_%s_%d_%d",processLabels[0].c_str(), ii, jj), "WriteDelete");
      NormalizeHist(histChargedIsoFakeTemplate_EB_binned[ii][jj]);
      NormalizeHist(histChargedIsoFakeTemplate_EE_binned[ii][jj]);
      NormalizeHist(histChargedIsoPromptTemplate_EB_binned[ii][jj]);
      NormalizeHist(histChargedIsoPromptTemplate_EE_binned[ii][jj]);
      file->WriteTObject(histChargedIsoFakeTemplate_EB_binned[ii][jj], Form("histChargedIsoFakeTemplate_EB_binned_%d_%d", ii, jj), "WriteDelete");
      file->WriteTObject(histChargedIsoFakeTemplate_EE_binned[ii][jj], Form("histChargedIsoFakeTemplate_EE_binned_%d_%d", ii, jj), "WriteDelete");
      file->WriteTObject(histChargedIsoPromptTemplate_EB_binned[ii][jj], Form("histChargedIsoPromptTemplate_EB_binned_%d_%d", ii, jj), "WriteDelete");
      file->WriteTObject(histChargedIsoPromptTemplate_EE_binned[ii][jj], Form("histChargedIsoPromptTemplate_EE_binned_%d_%d", ii, jj), "WriteDelete");
    }
  for(int ii=0; ii<NMRBins; ii++) {
    file->WriteTObject(histChargedIso_EB_MRbinned[ii], Form("histChargedIso_EB_MRbinned_%s_%d",processLabels[0].c_str(), ii), "WriteDelete");
    file->WriteTObject(histChargedIso_EE_MRbinned[ii], Form("histChargedIso_EE_MRbinned_%s_%d",processLabels[0].c_str(), ii), "WriteDelete");
    NormalizeHist(histChargedIsoFakeTemplate_EB_MRbinned[ii]);
    NormalizeHist(histChargedIsoFakeTemplate_EE_MRbinned[ii]);
    NormalizeHist(histChargedIsoPromptTemplate_EB_MRbinned[ii]);
    NormalizeHist(histChargedIsoPromptTemplate_EE_MRbinned[ii]);
    file->WriteTObject(histChargedIsoFakeTemplate_EB_MRbinned[ii], Form("histChargedIsoFakeTemplate_EB_MRbinned_%d", ii), "WriteDelete");
    file->WriteTObject(histChargedIsoFakeTemplate_EE_MRbinned[ii], Form("histChargedIsoFakeTemplate_EE_MRbinned_%d", ii), "WriteDelete");
    file->WriteTObject(histChargedIsoPromptTemplate_EB_MRbinned[ii], Form("histChargedIsoPromptTemplate_EB_MRbinned_%d", ii), "WriteDelete");
    file->WriteTObject(histChargedIsoPromptTemplate_EE_MRbinned[ii], Form("histChargedIsoPromptTemplate_EE_MRbinned_%d", ii), "WriteDelete");
  }
  for(int ii=0; ii<NRsqBins; ii++) {
    file->WriteTObject(histChargedIso_EB_Rsqbinned[ii], Form("histChargedIso_EB_Rsqbinned_%s_%d",processLabels[0].c_str(), ii), "WriteDelete");
    file->WriteTObject(histChargedIso_EE_Rsqbinned[ii], Form("histChargedIso_EE_Rsqbinned_%s_%d",processLabels[0].c_str(), ii), "WriteDelete");
    NormalizeHist(histChargedIsoFakeTemplate_EB_Rsqbinned[ii]);
    NormalizeHist(histChargedIsoFakeTemplate_EE_Rsqbinned[ii]);
    NormalizeHist(histChargedIsoPromptTemplate_EB_Rsqbinned[ii]);
    NormalizeHist(histChargedIsoPromptTemplate_EE_Rsqbinned[ii]);
    file->WriteTObject(histChargedIsoFakeTemplate_EB_Rsqbinned[ii], Form("histChargedIsoFakeTemplate_EB_Rsqbinned_%d", ii), "WriteDelete");
    file->WriteTObject(histChargedIsoFakeTemplate_EE_Rsqbinned[ii], Form("histChargedIsoFakeTemplate_EE_Rsqbinned_%d", ii), "WriteDelete");
    file->WriteTObject(histChargedIsoPromptTemplate_EB_Rsqbinned[ii], Form("histChargedIsoPromptTemplate_EB_Rsqbinned_%d", ii), "WriteDelete");
    file->WriteTObject(histChargedIsoPromptTemplate_EE_Rsqbinned[ii], Form("histChargedIsoPromptTemplate_EE_Rsqbinned_%d", ii), "WriteDelete");
  }

  file->WriteTObject(histChargedIsoTemplate_EB[0], Form("histChargedIsoTemplate_EB_%s",processLabels[0].c_str()), "WriteDelete");
  file->WriteTObject(histChargedIsoTemplate_EE[0], Form("histChargedIsoTemplate_EE_%s",processLabels[0].c_str()), "WriteDelete");


  file->Close();
  delete file;
}

void PhotonControlSample_FitChargedIso( int option = 10) {

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
    // datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDCuts/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_2016_GoodLumiGolden_26p4ifb_MR300Skim.root");  
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDIsoCuts/Skim/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_2016_GoodLumiGolden_26p4ifb_MR300Rsq0p15Skim.root");  
    fakeTemplateFiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDIsoCuts/Skim/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_2016_GoodLumiGolden_26p4ifb_MR300Rsq0p15Skim.root");
     // fakeTemplateFiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDIsoCuts/Skim/RunTwoRazorControlRegions_PhotonFull_QCD_HTBinned_1pb_weighted_MR300Rsq0p15Skim.root");
    promptTemplateFiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDIsoCuts/Skim/RunTwoRazorControlRegions_PhotonFull_GJets_HTBinned_1pb_weighted_MR300Rsq0p15Skim.root");
  }


  double lumi = 26400;

  //*********************************************************************
  //GJets Control Region
  //*********************************************************************
  if (option == 0) {
    RunPhotonControlSample_FitChargedIso(datafiles, fakeTemplateFiles, promptTemplateFiles, lumi,"MR300Rsq0p15",0,"MR300Rsq0p15");
  }
  if (option == 10) {
    RunPhotonControlSample_FitChargedIso(datafiles, fakeTemplateFiles, promptTemplateFiles, lumi,"Inclusive",0,"Inclusive");
  }


}

