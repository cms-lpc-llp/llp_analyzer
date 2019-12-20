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

#endif




void PlotDataAndStackedBkg( vector<TH1D*> hist , vector<string> processLabels, vector<int> color,  bool hasData, string varName, string label ) {

  TCanvas *cv =0;
  TLegend *legend = 0;

  cv = new TCanvas("cv","cv", 800,600);
  cv->SetHighLightColor(2);
  cv->SetFillColor(0);
  cv->SetBorderMode(0);
  cv->SetBorderSize(2);
  cv->SetLeftMargin(0.16);
  cv->SetRightMargin(0.3);
  cv->SetTopMargin(0.07);
  cv->SetBottomMargin(0.12);
  cv->SetFrameBorderMode(0);  

  TPad *pad1 = new TPad("pad1","pad1", 0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(0.04);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.60,0.50,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stack = new THStack();
  TH1D *histDataOverMC = (TH1D*)hist[0]->Clone("histDataOverMC");

  if (hasData) {
    for (int i = hist.size()-1; i >= 1; --i) {
      hist[i]->SetFillColor(color[i]);
      hist[i]->SetFillStyle(1001);
      
      if ( hist[i]->Integral() > 0) {
  	stack->Add(hist[i]);
      }
    }
  } else {
    for (int i = hist.size()-1; i >= 0; --i) {
      hist[i]->SetFillColor(color[i]);
      hist[i]->SetFillStyle(1001);
      
      if ( hist[i]->Integral() > 0) {
  	stack->Add(hist[i]);
      }
    }
  }

  for (uint i = 0 ; i < hist.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(hist[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(hist[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stack->GetHists()->GetEntries() > 0) {
    stack->Draw("hist");
    stack->GetHistogram()->GetXaxis()->SetTitle(((TH1D*)(stack->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stack->GetHistogram()->GetYaxis()->SetTitle(((TH1D*)(stack->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stack->GetHistogram()->GetYaxis()->SetTitleOffset(0.8);
    stack->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    stack->GetHistogram()->GetXaxis()->SetTitleSize(0.15);
    stack->SetMaximum( 1.2* fmax( stack->GetMaximum(), hist[0]->GetMaximum()) );
    stack->SetMinimum( 0.1 );

    if (hasData) {
      hist[0]->SetLineWidth(2);
      hist[0]->SetLineColor(color[0]);
      hist[0]->Draw("e1same");
    }
    legend->Draw();
  }
  cv->cd();
  cv->Update();


  TPad *pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
  pad2->Draw();
  pad2->cd();
    
  for (int b=0; b<histDataOverMC->GetXaxis()->GetNbins()+2; ++b) {
    double data = 0;
    if (hasData) {
      data = hist[0]->GetBinContent(b);
    }
    double MC = 0;
    double MCErrSqr = 0;
    if (hasData) {
      for (uint i = 1 ; i < hist.size(); ++i) {
	MC += hist[i]->GetBinContent(b);
	MCErrSqr += pow(hist[i]->GetBinError(b),2);
      }
    } else {
      MC = 1;
    }
      
    if (MC > 0) {
      histDataOverMC->SetBinContent(b, data / MC);
      histDataOverMC->SetBinError(b, (data / MC)*sqrt(1/data + MCErrSqr/pow(MC,2) ));
    } else {
      histDataOverMC->SetBinContent(b, 0);
      histDataOverMC->SetBinError(b, 0);
    }
    //cout << "bin " << b << " : " << histDataOverMC->GetBinContent(b) << " " << histDataOverMC->GetBinError(b) << "\n";
  }

  histDataOverMC->GetYaxis()->SetTitle("Data/MC");
  histDataOverMC->GetYaxis()->SetNdivisions(306);
  histDataOverMC->GetYaxis()->SetTitleSize(0.10);
  histDataOverMC->GetYaxis()->SetTitleOffset(0.3);
  histDataOverMC->GetYaxis()->SetRangeUser(0.0,1.5);
  histDataOverMC->GetYaxis()->SetLabelSize(0.10);
  histDataOverMC->GetXaxis()->SetLabelSize(0.125);
  histDataOverMC->GetXaxis()->SetTitleSize(0.15);
  histDataOverMC->GetXaxis()->SetTitleOffset(1.0);
  histDataOverMC->SetStats(false);
  histDataOverMC->Draw("e1");
  
  pad1->SetLogy(false);
  cv->SaveAs(Form("Razor_VetoLeptonCR_%s%s.gif",varName.c_str(), label.c_str()));
  
  pad1->SetLogy(true);
  cv->SaveAs(Form("Razor_VetoLeptonCR_%s%s_Logy.gif",varName.c_str(),label.c_str()));

}

//=== MAIN MACRO ================================================================================================= 


void RunSelectVetoLeptonControlSample( vector<string> datafiles, vector<vector<string> > bkgfiles, vector<string> bkgLabels, 
				       vector<int> bkgColors, 
				       double lumi, string option, int channelOption = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  float MR = 0;
  float Rsq = 0;
  float mll = 0;


  bool printdebug = false;

  TFile *NVtxWeightFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/data/NVtxReweight_ZToMuMu_2015C_25ns.root", "READ");
  TH1D *NVtxWeightHist = (TH1D*)NVtxWeightFile->Get("NVtxReweight");
  assert(NVtxWeightHist);
  
  TFile *pileupWeightFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Run1PileupWeights.root", "READ");
  TH1D *pileupWeightHist = (TH1D*)pileupWeightFile->Get("PUWeight_Run1");
  assert(pileupWeightHist);

  TFile *eleEffSFFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/ScaleFactors/Run1/ElectronSelection_Run2012ReReco_53X.root","READ");
  TH2D *eleEffSFHist = (TH2D*)eleEffSFFile->Get("sfLOOSE");
  assert(eleEffSFHist);

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  const int NMRBins = 10;
  const int NRsqBins = 9;
  double MRBins[NMRBins] = {300, 350, 400, 450, 500, 550, 700, 900, 1200, 4000};
  double RsqBins[NRsqBins] = {0.15,0.175,0.20,0.225, 0.25,0.30,0.41,0.52,1.5};  

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

  vector<TH1D*> histMR;
  vector<TH1D*> histRsq;
  vector<TH1D*> histLep1Pt;
  vector<TH1D*> histLep1Eta;
  vector<TH1D*> histLep1Activity;
  vector<TH1D*> histLep1MT;
  vector<TH1D*> histMET;
  vector<TH1D*> histNJets40;
  vector<TH1D*> histNJets80;
  vector<TH1D*> histNBtags;
  vector<TH2F*> histMRVsRsq;

  assert (inputfiles.size() == processLabels.size());
  bool addedData = false;
  int DataHistIndex = -1;
  for (uint i=0; i < inputfiles.size(); ++i) {   
    histMR.push_back(new TH1D(Form("histMR_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", 40, 400, 2400));
    histRsq.push_back(new TH1D(Form("histRsq_%s",processLabels[i].c_str()), "; R^{2} ; Number of Events", 25, 0.15, 1.65));
    histLep1Pt.push_back(new TH1D(Form("histLep1Pt_%s",processLabels[i].c_str()), "; Lepton p_{T} [GeV/c] ; Number of Events", 20, 0, 200));
    histLep1Eta.push_back(new TH1D(Form("histLep1Eta_%s",processLabels[i].c_str()), "; Lepton #eta ; Number of Events", 10, -2.5, 2.5));
    histLep1Activity.push_back(new TH1D(Form("histLep1Activity_%s",processLabels[i].c_str()), "; Lepton Activity ; Number of Events", 100, 0, 200));
    histLep1MT.push_back(new TH1D(Form("histLep1MT_%s",processLabels[i].c_str()), "; Lep1MT [GeV/c] ; Number of Events", 40, 0, 200));
    histMET.push_back(new TH1D(Form("histMET_%s",processLabels[i].c_str()), "; MET [GeV/c] ; Number of Events", 40, 0, 1000));
    histNJets40.push_back(new TH1D(Form("histNJets40_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 40); Number of Events", 15, -0.5, 14.5));
    histNJets80.push_back(new TH1D(Form("histNJets80_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 80); Number of Events", 10, -0.5, 9.5));
    histNBtags.push_back(new TH1D(Form("histNBtags_%s",processLabels[i].c_str()), "; Number of B tags; Number of Events", 10, -0.5,9.5));
    histMRVsRsq.push_back(new TH2F(Form("histMRVsRsq_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; R^{2}; Number of Events", NMRBins-1, MRBins, NRsqBins-1, RsqBins));
    histMR[i]->Sumw2();
    histRsq[i]->Sumw2();
    histLep1Pt[i]->Sumw2();
    histLep1Eta[i]->Sumw2();
    histLep1Activity[i]->Sumw2();
    histLep1MT[i]->Sumw2();
    histNJets40[i]->Sumw2();
    histNJets80[i]->Sumw2();
    histNBtags[i]->Sumw2();  
  }
 
  double dataYield = 0;
  double MCYield = 0;
  double MCWJetsYield = 0;

  vector<pair<UInt_t,UInt_t> > RunAndEvent;

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  for (uint i=0; i < inputfiles.size(); ++i) {
    cout << "p: " << i << " : " << inputfiles[i].size() << "\n";
    for (uint j=0; j < inputfiles[i].size(); ++j) {
      ControlSampleEvents *events = new ControlSampleEvents;
      events->LoadTree(inputfiles[i][j].c_str(),ControlSampleEvents::kTreeType_OneLepton_Full);

      bool isData = false;
      if ( processLabels[i] == "Data") isData = true;

      cout << "process: " << processLabels[i] << " | file " << inputfiles[i][j] << " | Total Entries: " << events->tree_->GetEntries() << "\n";
      for(UInt_t ientry=0; ientry < events->tree_->GetEntries(); ientry++) {       	
	events->tree_->GetEntry(ientry);
      
	if (ientry % 1000000 == 0) cout << "Event " << ientry << endl;      

	double puWeight = 1;      
	double weight = 1;

	if (!isData) {
	  //puWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(events->NPU_0));
	  puWeight = NVtxWeightHist->GetBinContent(NVtxWeightHist->GetXaxis()->FindFixBin(events->NPV)); //temporarily use nvtx weighting
	  puWeight = 1.0;
	  weight = lumi * events->weight * puWeight;
	}


	//******************************
	//Trigger Selection : Razor Hadronic Triggers
	//******************************
	bool passTrigger = false;
	if (isData) {
	  if (events->HLTDecision[132] ==true || events->HLTDecision[133] ==true) passTrigger = true;
	} else {
	  if (events->HLTDecision[136] ==true || events->HLTDecision[137] ==true) passTrigger = true;
	}
	if (!passTrigger) continue;
 

	//******************************
	//Selection Cuts 
	//******************************
	if (!( abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13 ) ) continue;

	//lepton selection
	if (! (events->lep1PassVeto) ) continue;

	//MET cut
	if (!(events->NJets80 >= 2 && events->MR > 400 && events->Rsq > 0.25)) continue;

                
	//******************************
	//ChannelOptions
	//******************************
	// Electron Channel
	if (channelOption == 11 &&
	    !(abs(events->lep1Type) == 11)
	    ) continue;
      
	// Muon Channel
	if (channelOption == 13 &&
	    !(abs(events->lep1Type) == 13)
	    ) continue;
      
 
	// //******************************
	// //Apply Scale Factors
	// //******************************
	// if (!(hasData && i==0)) {
	// 	double triggerEffScaleFactor = 1.0;

	// 	double leptonEffScaleFactor = 1.0;
	// 	if (abs(events->lep1Type) == 11) {
	// 	  leptonEffScaleFactor = 0.96; //approximate guess for the average looking at all the scale factors in pt-eta bins
	// 	} else if (abs(events->lep1Type) == 13) {
	// 	  leptonEffScaleFactor = 1.00; //approximate guess for the average looking at all the scale factors in pt-eta bins
	// 	}
		
	// 	if (abs(events->lep1Type) == 11) {
	// 	  triggerEffScaleFactor = 0.97;
	// 	} else if (abs(events->lep1Type) == 13) {
	// 	  triggerEffScaleFactor= 0.97; //approximate guess for the average looking at all the scale factors in pt-eta bins
	// 	}

	// 	//b-tagging scale factors
	// 	double btagScaleFactor = 1.0;
	// 	double bjet1EventScaleFactor = 1.0;
	// 	double bjet2EventScaleFactor = 1.0;
	// 	if (events->bjet1.Pt() > 20) {
	// 	  double bjet1ScaleFactor = 0.938887 + 0.00017124 * events->bjet1.Pt() + (-2.76366e-07) * events->bjet1.Pt() * events->bjet1.Pt() ;
	// 	  double MCEff = 1.0;
	// 	  if (events->bjet1.Pt() < 50) MCEff = 0.65;
	// 	  else if (events->bjet1.Pt() < 80) MCEff = 0.70;
	// 	  else if (events->bjet1.Pt() < 120) MCEff = 0.73;
	// 	  else if (events->bjet1.Pt() < 210) MCEff = 0.73;
	// 	  else MCEff = 0.66;				 
	// 	  if (events->bjet1PassMedium) bjet1EventScaleFactor = bjet1ScaleFactor;
	// 	  else bjet1EventScaleFactor = ( 1/MCEff - bjet1ScaleFactor) / ( 1/MCEff - 1);
	// 	}
	// 	if (events->bjet2.Pt() > 20) {
	// 	  double bjet2ScaleFactor = 0.938887 + 0.00017124 * events->bjet2.Pt() + (-2.76366e-07) * events->bjet2.Pt() * events->bjet1.Pt() ;
	// 	  double MCEff = 1.0;
	// 	  if (events->bjet2.Pt() < 50) MCEff = 0.65;
	// 	  else if (events->bjet2.Pt() < 80) MCEff = 0.70;
	// 	  else if (events->bjet2.Pt() < 120) MCEff = 0.73;
	// 	  else if (events->bjet2.Pt() < 210) MCEff = 0.73;
	// 	  else MCEff = 0.66;				 
	// 	  if (events->bjet2PassMedium) bjet2EventScaleFactor = bjet2ScaleFactor;
	// 	  else bjet2EventScaleFactor = ( 1/MCEff - bjet2ScaleFactor) / ( 1/MCEff - 1);
	// 	}
	// 	btagScaleFactor = bjet1EventScaleFactor * bjet2EventScaleFactor;


	// 	weight *= leptonEffScaleFactor;
	// 	weight *= btagScaleFactor;


	// 	if (processLabels[i] == "TTJets") {
	// 	  weight *= TTBarScaleFactorsHist->GetBinContent( TTBarScaleFactorsHist->GetXaxis()->FindFixBin(fmin(fmax(events->MR,300.1),699.9)) ,  
	// 							  TTBarScaleFactorsHist->GetYaxis()->FindFixBin(fmin(fmax(events->Rsq,0.1501),1.499)) );
	// 	}

	// }



	//******************************
	//Fill histograms
	//******************************
	if (isData) {

	  dataYield += 1.0;
	  histLep1Pt[i]->Fill(events->lep1.Pt());
	  histLep1Eta[i]->Fill(events->lep1.Eta());
	  histLep1Activity[i]->Fill(events->lep1Activity);
	  histLep1MT[i]->Fill(events->lep1MT);
	  histMET[i]->Fill(events->MET);
	  histNJets40[i]->Fill( events->NJets40 );
	  histNJets80[i]->Fill( events->NJets80 );
	  histNBtags[i]->Fill( events->NBJetsMedium );

	  if (events->NJets80 >= 2 && events->MR > 0) {
	    histMR[i]->Fill(events->MR);
	    histRsq[i]->Fill(events->Rsq);
	  }

	  histMRVsRsq[i]->Fill(events->MR,events->Rsq);


	} else {
	  MCYield += weight;
	  histLep1Pt[i]->Fill(events->lep1.Pt(), weight);
	  histLep1Eta[i]->Fill(events->lep1.Eta(), weight);
	  histLep1Activity[i]->Fill(events->lep1Activity, weight);
	  histLep1MT[i]->Fill(events->lep1MT, weight);
	  histMET[i]->Fill(events->MET, weight);
	  histNJets40[i]->Fill( events->NJets40 , weight);
	  histNJets80[i]->Fill( events->NJets80 , weight);
	  histNBtags[i]->Fill( events->NBJetsMedium , weight);
	
	  if (events->NJets80 >= 2 && events->MR > 0) {
	    histMR[i]->Fill(events->MR, weight );
	    histRsq[i]->Fill(events->Rsq, weight );
	  }

	  histMRVsRsq[i]->Fill(events->MR,events->Rsq, weight);

	}
      }
      delete events;
    }
  }



  //--------------------------------------------------------------------------------------------------------------
  // Compute Expected Statistical Uncertainty
  //==============================================================================================================
  // TH2F *statUnc_MRVsRsq = (TH2F*)(histMRVsRsq[0]->Clone("statUnc_MRVsRsq"));
  // for (int i=0; i<statUnc_MRVsRsq->GetXaxis()->GetNbins()+1;i++) {
  //   for (int j=0; j<statUnc_MRVsRsq->GetYaxis()->GetNbins()+1;j++) {
  //     if (statUnc_MRVsRsq->GetBinContent(i,j) > 1) {
  // 	statUnc_MRVsRsq->SetBinContent(i,j,1.0 / sqrt(statUnc_MRVsRsq->GetBinContent(i,j)));  	
  //     } else {
  // 	statUnc_MRVsRsq->SetBinContent(i,j,0);
  //     }
  //   }
  // }

 
  //--------------------------------------------------------------------------------------------------------------
  // Make Plots
  //==============================================================================================================


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;
  TLegend *legend = 0;

  //*******************************************************************************************
  //MR
  //*******************************************************************************************
  PlotDataAndStackedBkg( histMR, processLabels, color, true, "MR", Label);
  PlotDataAndStackedBkg( histRsq, processLabels, color, true, "Rsq", Label);
  PlotDataAndStackedBkg( histNJets40, processLabels, color, true, "NJets40", Label);
  PlotDataAndStackedBkg( histNJets80, processLabels, color, true, "NJets80", Label);
  PlotDataAndStackedBkg( histNBtags, processLabels, color, true, "NBtags", Label);
  PlotDataAndStackedBkg( histLep1Pt, processLabels, color, true, "Lep1Pt", Label);
  PlotDataAndStackedBkg( histLep1Eta, processLabels, color, true, "Lep1Eta", Label);
  PlotDataAndStackedBkg( histLep1Activity, processLabels, color, true, "Lep1Activity", Label);
  PlotDataAndStackedBkg( histLep1MT, processLabels, color, true, "Lep1MT", Label);
  PlotDataAndStackedBkg( histMET, processLabels, color, true, "MET", Label);



  //--------------------------------------------------------------------------------------------------------------
  // Tables
  //==============================================================================================================
  cout << "For Luminosity = " << lumi << " pb^-1\n";
  //cout << "TTJets: " << 

  cout << "Yield \n";
  cout << "Data: " << dataYield << "\n";
  cout << "MC: " << MCYield << "\n";

  // //--------------------------------------------------------------------------------------------------------------
  // // Output
  // //==============================================================================================================
  TFile *file = TFile::Open(("VetoLeptonControlRegionPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNJets40[i], Form("histNJets40_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNBtags[i], Form("histNBtags_%s",processLabels[i].c_str()), "WriteDelete");
  }


  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNJets40[i], Form("histNJets40_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNBtags[i], Form("histNBtags_%s",processLabels[i].c_str()), "WriteDelete");
  }
 
   
  file->Close();
  delete file;       
    

}







void SelectVetoLeptonControlSample( int option = -1) {

  vector<string> datafiles;
  vector<vector<string> > bkgfiles;
  vector<string> processLabels;
  vector<int> colors;


  datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/VetoLepton_1p17/RazorSkim/RunTwoRazorControlRegions_VetoLeptonFull_HTMHT_Run2015B_GoodLumi.root_RazorSkim.root");
  datafiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/VetoLepton_1p17/RazorSkim/RunTwoRazorControlRegions_VetoLeptonFull_HTMHT_Run2015C_GoodLumi.root_RazorSkim.root");
 
  vector<string> bkgfiles_ttbar;
  bkgfiles_ttbar.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/VetoLepton_1p17/RazorSkim/RunTwoRazorControlRegions_VetoLeptonFull_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root_RazorSkim.root");
  vector<string> bkgfiles_wjets;
  bkgfiles_wjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/VetoLepton_1p17/RazorSkim/RunTwoRazorControlRegions_VetoLeptonFull_WJetsToLNuHTBinned_1pb_weighted.root_RazorSkim.root");
  vector<string> bkgfiles_singletop;
  bkgfiles_singletop.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/VetoLepton_1p17/RazorSkim/RunTwoRazorControlRegions_VetoLeptonFull_ST_1pb_weighted.root_RazorSkim.root");
  vector<string> bkgfiles_qcd;
  bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/VetoLepton_1p17/RazorSkim/RunTwoRazorControlRegions_VetoLeptonFull_QCDHTBinned_1pb_weighted.root_RazorSkim.root");
  vector<string> bkgfiles_znunu;
  bkgfiles_znunu.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/VetoLepton_1p17/RazorSkim/RunTwoRazorControlRegions_VetoLeptonFull_ZJetsToNuNuHTBinned_1pb_weighted.root_RazorSkim.root");
  
  bkgfiles.push_back(bkgfiles_ttbar);
  bkgfiles.push_back(bkgfiles_wjets);
  bkgfiles.push_back(bkgfiles_singletop);
  bkgfiles.push_back(bkgfiles_qcd);
  bkgfiles.push_back(bkgfiles_znunu);

  processLabels.push_back("TTJets");  
  processLabels.push_back("WJets");  
  processLabels.push_back("SingleTop");
  processLabels.push_back("QCD");  
  processLabels.push_back("ZNuNu");  
  //processLabels.push_back("DY");
  //processLabels.push_back("VV");
  //processLabels.push_back("TT+V");
  
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kOrange+1);
  colors.push_back(kViolet);
  colors.push_back(kCyan+2);

  


  //*********************************************************************
  //Veto Electron Selection
  //*********************************************************************
  if (option == 0) {
    // RunSelectVetoLeptonControlSample(datafiles, bkgfiles, processLabels,  colors, 16.1, "",11,"VetoElectron");
    // RunSelectVetoLeptonControlSample(datafiles, bkgfiles, processLabels,  colors, 39.9, "",11,"VetoElectron");
    RunSelectVetoLeptonControlSample(datafiles, bkgfiles, processLabels,  colors, 36.0, "",11,"VetoElectron");
  }
 
  //*********************************************************************
  //Veto Muon Control Region
  //*********************************************************************
  if (option == 1) {
    // RunSelectVetoLeptonControlSample(datafiles, bkgfiles,processLabels,  colors, 16.1, "",13,"VetoMuon");
    // RunSelectVetoLeptonControlSample(datafiles, bkgfiles,processLabels,  colors, 39.9, "",13,"VetoMuon");
    RunSelectVetoLeptonControlSample(datafiles, bkgfiles,processLabels,  colors, 36.0, "",13,"VetoMuon");
  }
 
  //*********************************************************************
  //Veto Electron and Veto Muon Combined Control Region
  //*********************************************************************
  if (option == -1) {
    // RunSelectVetoLeptonControlSample(datafiles, bkgfiles, processLabels,  colors, 16.1, "",0,"VetoLepton");
    // RunSelectVetoLeptonControlSample(datafiles, bkgfiles, processLabels,  colors, 39.9, "",0,"VetoLepton");
    RunSelectVetoLeptonControlSample(datafiles, bkgfiles, processLabels,  colors, 36.0, "",0,"VetoLepton");
  }

}

