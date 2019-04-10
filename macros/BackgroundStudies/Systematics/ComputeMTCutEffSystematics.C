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


//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
TH1D* NormalizeHist(TH1D *originalHist) {
  TH1D* hist = (TH1D*)originalHist->Clone((string(originalHist->GetName())+"_normalized").c_str());
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return hist;
}


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

  //****************************
  //Add CMS and Lumi Labels
  //****************************
  // lumi_13TeV = "42 pb^{-1}";
  lumi_13TeV = "157 pb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  CMS_lumi(pad1,4,0);

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
  histDataOverMC->GetYaxis()->SetRangeUser(0.5,1.5);
  histDataOverMC->GetYaxis()->SetLabelSize(0.10);
  histDataOverMC->GetXaxis()->SetLabelSize(0.125);
  histDataOverMC->GetXaxis()->SetTitleSize(0.15);
  histDataOverMC->GetXaxis()->SetTitleOffset(1.0);
  histDataOverMC->SetStats(false);
  histDataOverMC->Draw("e1");
  
  pad1->SetLogy(false);
  cv->SaveAs(Form("Razor_WJetsCR_%s%s.gif",varName.c_str(), label.c_str()));
  
  pad1->SetLogy(true);
  cv->SaveAs(Form("Razor_WJetsCR_%s%s_Logy.gif",varName.c_str(),label.c_str()));

}

//=== MAIN MACRO ================================================================================================= 


void doComputeMTCutEffSystematics( string file, string label = "") {
  
  TH1D *histLep1MT = new TH1D("histLep1MT", "; Lep1MT [GeV/c] ; Number of Events", 100, 0, 200);
  TH1D *histLep1MTScaleVariation1 = new TH1D("histLep1MTScaleVariation1", "; Lep1MT [GeV/c] ; Number of Events", 100, 0, 200);
  TH1D *histLep1MTScaleVariation2 = new TH1D("histLep1MTScaleVariation2", "; Lep1MT [GeV/c] ; Number of Events", 100, 0, 200);
  TH1D *histLep1MTScaleVariation3 = new TH1D("histLep1MTScaleVariation3", "; Lep1MT [GeV/c] ; Number of Events", 100, 0, 200);
  TH1D *histLep1MTScaleVariation4 = new TH1D("histLep1MTScaleVariation4", "; Lep1MT [GeV/c] ; Number of Events", 100, 0, 200);
  TH1D *histLep1MTScaleVariation6 = new TH1D("histLep1MTScaleVariation6", "; Lep1MT [GeV/c] ; Number of Events", 100, 0, 200);
  TH1D *histLep1MTScaleVariation8 = new TH1D("histLep1MTScaleVariation8", "; Lep1MT [GeV/c] ; Number of Events", 100, 0, 200);

  histLep1MT->Sumw2();
  histLep1MTScaleVariation1->Sumw2();
  histLep1MTScaleVariation2->Sumw2();
  histLep1MTScaleVariation3->Sumw2();
  histLep1MTScaleVariation4->Sumw2();
  histLep1MTScaleVariation6->Sumw2();
  histLep1MTScaleVariation8->Sumw2();
  
  vector<double> TotalYield;
  vector<double> PassMTCutYield;
  for (int i=0; i<7;++i) {
    TotalYield.push_back(0.0);
    PassMTCutYield.push_back(0.0);
  }


  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  ControlSampleEvents *events = new ControlSampleEvents;
  events->LoadTree(file.c_str(),ControlSampleEvents::kTreeType_OneLepton_Reduced);
  
  cout << " | file " << file << " | Total Entries: " << events->tree_->GetEntries() << "\n";
  for(UInt_t ientry=0; ientry < events->tree_->GetEntries(); ientry++) {       	
    events->tree_->GetEntry(ientry);
    
    if (ientry % 1000000 == 0) cout << "Event " << ientry << endl;      
    
    
   	
    //******************************
    //Trigger Selection
    //******************************
    bool passTrigger = false;
    
    //Single Lepton Triggers:
    //Use Single Lepton Triggers
    if ( events->HLTDecision[3] || events->HLTDecision[8] || events->HLTDecision[12] 
	 || events->HLTDecision[11] || events->HLTDecision[14]
	 )  
      passTrigger = true;
    
    //use MC triggers
    if ( events->HLTDecision[17] || events->HLTDecision[18] || events->HLTDecision[19] || 
	 events->HLTDecision[20] ||
	 events->HLTDecision[26] || events->HLTDecision[27]	  
	 ) passTrigger = true;
    
    
    if (!passTrigger) continue;
	

    //******************************
    //Selection Cuts 
    //******************************
    if (!( abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13 ) ) continue;
    
    
    //lepton selection
    if (! (events->lep1Pt > 30) ) continue;
    if (! (events->lep1PassTight) ) continue;
    
    //MR Bins
    //if (!(events->MR > 400 && events->MR <= 600)) continue;
    //if (!(events->MR > 600 && events->MR <= 800)) continue;
    //if (!(events->MR > 800 && events->MR <= 1000)) continue;
    //if (!(events->MR > 1000)) continue;
 
	
    //******************************
    //Fill histograms
    //******************************
    histLep1MT->Fill(events->lep1MT);
    histLep1MTScaleVariation1->Fill(events->lep1MT,events->ScaleVariationWeight1);
    histLep1MTScaleVariation2->Fill(events->lep1MT,events->ScaleVariationWeight2);
    histLep1MTScaleVariation3->Fill(events->lep1MT,events->ScaleVariationWeight3);
    histLep1MTScaleVariation4->Fill(events->lep1MT,events->ScaleVariationWeight4);
    histLep1MTScaleVariation6->Fill(events->lep1MT,events->ScaleVariationWeight6);
    histLep1MTScaleVariation8->Fill(events->lep1MT,events->ScaleVariationWeight8);

    //cout << events->lep1MT << " " << events->ScaleVariationWeight1 << "\n";

    TotalYield[0] += 1.0;
    TotalYield[1] += events->ScaleVariationWeight1;
    TotalYield[2] += events->ScaleVariationWeight2;
    TotalYield[3] += events->ScaleVariationWeight3;
    TotalYield[4] += events->ScaleVariationWeight4;
    TotalYield[5] += events->ScaleVariationWeight6;
    TotalYield[6] += events->ScaleVariationWeight8;
    if (events->lep1MT > 100) {
      PassMTCutYield[0] += 1.0;
      PassMTCutYield[1] += events->ScaleVariationWeight1;
      PassMTCutYield[2] += events->ScaleVariationWeight2;
      PassMTCutYield[3] += events->ScaleVariationWeight3;
      PassMTCutYield[4] += events->ScaleVariationWeight4;
      PassMTCutYield[5] += events->ScaleVariationWeight6;
      PassMTCutYield[6] += events->ScaleVariationWeight8;      
    }
    
  }

  delete events;

  //--------------------------------------------------------------------------------------------------------------
  // Make Plots
  //==============================================================================================================
  histLep1MT = NormalizeHist(histLep1MT);
  histLep1MTScaleVariation1 = NormalizeHist(histLep1MTScaleVariation1);
  histLep1MTScaleVariation2 = NormalizeHist(histLep1MTScaleVariation2);
  histLep1MTScaleVariation3 = NormalizeHist(histLep1MTScaleVariation3);
  histLep1MTScaleVariation4 = NormalizeHist(histLep1MTScaleVariation4);
  histLep1MTScaleVariation6 = NormalizeHist(histLep1MTScaleVariation6);
  histLep1MTScaleVariation8 = NormalizeHist(histLep1MTScaleVariation8);


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;
  TLegend *legend = 0;

  cv = new TCanvas("cv","cv",800,600);
  legend = new TLegend(0.60,0.50,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  histLep1MT->Draw("hist");
  histLep1MTScaleVariation1->SetLineColor(kRed);
  histLep1MTScaleVariation2->SetLineColor(kBlue);
  histLep1MTScaleVariation1->Draw("histsame");
  histLep1MTScaleVariation2->Draw("histsame");

  cv->SaveAs("MTScaleVariation.gif");
  

  // //--------------------------------------------------------------------------------------------------------------
  // // // Output
  // // //==============================================================================================================
  // TFile *file = TFile::Open(("WJetsControlRegionPlots"+Label+".root").c_str(), "UPDATE");
  // file->cd();
  // for(int i=0; i<int(inputfiles.size()); i++) {
  //   file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
  //   file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
  //   file->WriteTObject(histNJets40[i], Form("histNJets40_%s",processLabels[i].c_str()), "WriteDelete");
  //   file->WriteTObject(histNBtags[i], Form("histNBtags_%s",processLabels[i].c_str()), "WriteDelete");
  // }
  // for(int i=0; i<int(inputfiles.size()); i++) {
  //   file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
  //   file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
  //   file->WriteTObject(histNJets40[i], Form("histNJets40_%s",processLabels[i].c_str()), "WriteDelete");
  //   file->WriteTObject(histNBtags[i], Form("histNBtags_%s",processLabels[i].c_str()), "WriteDelete");
  // }   
  // file->Close();
  // delete file;       
    
  double SumErrSqr = 0;
  for (int i=0; i<7;++i) {
    cout << i << " : " << PassMTCutYield[i] / TotalYield[i] << "\n";
  }
  SumErrSqr += pow( (PassMTCutYield[1] / TotalYield[1] - PassMTCutYield[0] / TotalYield[0]) / (PassMTCutYield[0] / TotalYield[0]), 2);
  SumErrSqr += pow( (PassMTCutYield[2] / TotalYield[2] - PassMTCutYield[0] / TotalYield[0]) / (PassMTCutYield[0] / TotalYield[0]), 2);
  SumErrSqr += pow( (PassMTCutYield[3] / TotalYield[3] - PassMTCutYield[0] / TotalYield[0]) / (PassMTCutYield[0] / TotalYield[0]), 2);
  SumErrSqr += pow( (PassMTCutYield[4] / TotalYield[4] - PassMTCutYield[0] / TotalYield[0]) / (PassMTCutYield[0] / TotalYield[0]), 2);
  SumErrSqr += pow( (PassMTCutYield[5] / TotalYield[5] - PassMTCutYield[0] / TotalYield[0]) / (PassMTCutYield[0] / TotalYield[0]), 2);
  SumErrSqr += pow( (PassMTCutYield[6] / TotalYield[6] - PassMTCutYield[0] / TotalYield[0]) / (PassMTCutYield[0] / TotalYield[0]), 2);

  cout << "Total Systematic: " << sqrt(SumErrSqr) << "\n";


}







void ComputeMTCutEffSystematics( int option = -1) {

  doComputeMTCutEffSystematics( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/OneLeptonReducedSystematics_1p20/RunTwoRazorControlRegions_OneLeptonReduced_WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root" );
  doComputeMTCutEffSystematics( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/OneLeptonReducedSystematics_1p20/RunTwoRazorControlRegions_OneLeptonReduced_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root" );

}

//Systematics:

//W+Jets: 
//MR > 600-800: 1.2%
//MR > 800-1000: 0.7%
//MR > 1000: 0.3%

//TTBar: 
//MR no cut: 0.3%
//MR > 600-800: 0.2%
//MR > 800-1000: 0.2%
//MR > 1000: 0.2%
