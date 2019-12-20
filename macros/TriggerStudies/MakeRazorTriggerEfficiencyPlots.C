//================================================================================================
//
// Simple Example
//
//root -l 'RazorAnalyzer/macros/ObjectStudies/MakeElectronEfficiencyPlots.C(2)'. Currently 2 is the only working option for 2016 data.
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TStyle.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TH1F.h>                
#include <TH2F.h>                
#include <TCanvas.h>    
#include <TChain.h>
#include <TGraphAsymmErrors.h>                
#include <TLegend.h>                
#include <TLatex.h>                
#include <TLine.h>                

#include "../ObjectStudies/EfficiencyUtils.hh"

#endif


//Axis
const float axisTitleSize = 0.06;
const float axisTitleOffset = .8;

const float axisTitleSizeRatioX   = 0.18;
const float axisLabelSizeRatioX   = 0.12;
const float axisTitleOffsetRatioX = 0.94;

const float axisTitleSizeRatioY   = 0.15;
const float axisLabelSizeRatioY   = 0.108;
const float axisTitleOffsetRatioY = 0.32;

//Margins
const float leftMargin   = 0.12;
const float rightMargin  = 0.05;
const float topMargin    = 0.07;
const float bottomMargin = 0.12;

//CMS STANDARD
TString CMSText = "CMS";
TString extraText   = "Preliminary";
//TString lumiText = "10.3 fb^{-1} (13 TeV)"; //2018Av1v2
TString lumiText = "13.5 fb^{-1} (13 TeV)"; //2018Av1v2v3
//TString lumiText = "6.8 fb^{-1} (13 TeV)"; //2018B
//TString lumiText = "3.8 fb^{-1} (13 TeV)"; //2018C
//TString lumiText = "24.3 fb^{-1} (13 TeV)"; //2018ABC
//TString lumiText = "35.9 fb^{-1} (13 TeV)";
//TString lumiText = "Simulation (13 TeV)";

bool AddCMS( TCanvas* C );

const int JET_CUT = 120;


bool PassSelection( bool *HLTDecision, int year, int wp, double HLTMR, double HLTRSQ ) {

    bool pass = false;

    //***********************************
    //2016 Trigger Bits
    //***********************************
    if (year == 2016) {
      //trigger
      if (wp == 1) {
        pass = HLTDecision[163]; // PFHT900
      }
      if (wp == 2) {
        pass = HLTDecision[151]; // PFHT200
      }
      if (wp == 3) {
        pass = HLTDecision[172]; // Rsq0p25
      }
      if (wp == 4) {
        pass = HLTDecision[166];  // RsqMR270
      }
      if (wp == 5) {
        pass = HLTDecision[168];
      }
      if (wp == 7) {
        pass = HLTDecision[170];
      }
      if (wp == 13) {
        pass = HLTDecision[164] || HLTDecision[166];
      }
      if (wp == 16) {
        pass = HLTDecision[166] || HLTDecision[164] || HLTDecision[165] || HLTDecision[167];
      }
      if (wp == 20) {
        pass = HLTDecision[166] || HLTDecision[164] || HLTDecision[165] || HLTDecision[167] || HLTDecision[162] || HLTDecision[163]	 
	  || HLTDecision[125] || HLTDecision[126] ||  HLTDecision[134] || HLTDecision[135]
	  ;
      }
    }

    //***********************************
    //2017 Trigger Bits
    //***********************************
    if (year == 2017) {
      if (wp == 6) {
	pass = HLTDecision[108] || HLTDecision[110] ; // dijet || quad-jet
	//pass = HLTDecision[108] || HLTDecision[110] || HLTDecision[106] || HLTDecision[78] || HLTDecision[79] || HLTDecision[80] || HLTDecision[81] || HLTDecision[82] || HLTDecision[83] || HLTDecision[84] || HLTDecision[85] || HLTDecision[86] || HLTDecision[87] || HLTDecision[88];
      }
      if (wp == 7) {
	pass = HLTDecision[109] || HLTDecision[111] ;
      }
      if (wp == 8) {
	pass = bool(HLTMR>200 && HLTRSQ>0.09 && (HLTMR+300)*(HLTRSQ+0.25)>300);
      }    
      if (wp == 9) {
	pass = HLTDecision[108] || HLTDecision[110] || HLTDecision[106] ; // rsq-only || dijet || quad-jet 
      }
    }
    
    return pass;  
}




void plotRazorTriggerEfficiency() {
    TCanvas *cv =0;
    TLegend *legend =0;
    TLatex *boxLabel = 0;

    TFile *fileTrigger_RsqMR260_Rsq0p09_MR200_TTJets = new TFile("Efficiency_RazorTrigger_RsqMR260_Rsq0p09_MR200_All_TTJets_25ns.root","READ");
    TFile *fileTrigger_RsqMR300_Rsq0p09_MR200_TTJets = new TFile("Efficiency_RazorTrigger_RsqMR300_Rsq0p09_MR200_All_TTJets_25ns.root","READ");
    TFile *fileTrigger_RsqMR260_Rsq0p09_MR200_WJets = new TFile("Efficiency_RazorTrigger_RsqMR260_Rsq0p09_MR200_All_WJets_25ns.root","READ");
    TFile *fileTrigger_RsqMR300_Rsq0p09_MR200_WJets = new TFile("Efficiency_RazorTrigger_RsqMR300_Rsq0p09_MR200_All_WJets_25ns.root","READ");
    TFile *fileTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData = new TFile("Efficiency_RazorTrigger_RsqMR240_Rsq0p09_MR200_All_AllLeptonData_2016B.root","READ");
    TFile *fileTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData = new TFile("Efficiency_RazorTrigger_RsqMR270_Rsq0p09_MR200_All_AllLeptonData_2016B.root","READ");

    TH2F* effPtTrigger_RsqMR260_Rsq0p09_MR200_TTJets = (TH2F*)fileTrigger_RsqMR260_Rsq0p09_MR200_TTJets->Get("Efficiency_MRRsq");
    TH2F* effPtTrigger_RsqMR300_Rsq0p09_MR200_TTJets = (TH2F*)fileTrigger_RsqMR300_Rsq0p09_MR200_TTJets->Get("Efficiency_MRRsq");
    TH2F* effPtTrigger_RsqMR260_Rsq0p09_MR200_WJets = (TH2F*)fileTrigger_RsqMR260_Rsq0p09_MR200_WJets->Get("Efficiency_MRRsq");
    TH2F* effPtTrigger_RsqMR300_Rsq0p09_MR200_WJets = (TH2F*)fileTrigger_RsqMR300_Rsq0p09_MR200_WJets->Get("Efficiency_MRRsq");
    TH2F* effPtTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData = (TH2F*)fileTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->Get("Efficiency_MRRsq");
    TH2F* effPtTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData = (TH2F*)fileTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData->Get("Efficiency_MRRsq");

    TGraphAsymmErrors* effMRTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData = (TGraphAsymmErrors*)fileTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->Get("Efficiency_MR");
    TGraphAsymmErrors* effRsqTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData = (TGraphAsymmErrors*)fileTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->Get("Efficiency_Rsq");

    TLine *vline = 0;
    TLine *hline = 0;

    cv = new TCanvas("cv","cv", 800,600);
    effMRTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->SetTitle("; M_{R} [GeV/c^{2}]; Efficiency");
    effMRTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->Draw("ap");
    effMRTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->SetMarkerStyle(20);
    effMRTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->SetMarkerSize(1);
    effMRTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetYaxis()->SetRangeUser(0.8,1.05);
    effMRTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetXaxis()->SetRangeUser(300,1100);
    effMRTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetXaxis()->SetTitleSize(0.05);
    effMRTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetYaxis()->SetTitleSize(0.05);
    effMRTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetXaxis()->SetTitleOffset(0.75);
    effMRTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetYaxis()->SetTitleOffset(0.90);
    vline = new TLine(500,0.8,500,1.05);
    vline->SetLineStyle(2);
    vline->SetLineWidth(3);
    vline->SetLineColor(kBlue);
    hline = new TLine(300,1.0,1100,1.0);
    hline->SetLineStyle(3);
    hline->SetLineWidth(3);
    hline->SetLineColor(kBlack);
    vline->Draw();
    hline->Draw();
    cv->SaveAs("RazorTriggerEfficiencyVsMR_RsqMR240_Rsq0p09_MR200_AllLeptonData_2016B.gif");
    cv->SaveAs("RazorTriggerEfficiencyVsMR_RsqMR240_Rsq0p09_MR200_AllLeptonData_2016B.pdf");


    cv = new TCanvas("cv","cv", 800,600);
    effRsqTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->SetTitle("; R^{2};  Efficiency");
    effRsqTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->Draw("ap");
    effRsqTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->SetMarkerStyle(20);
    effRsqTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->SetMarkerSize(1);
    effRsqTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetYaxis()->SetRangeUser(0.8,1.05);
    effRsqTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetXaxis()->SetRangeUser(0.15,1.5);
    effRsqTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetXaxis()->SetTitleSize(0.05);
    effRsqTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetYaxis()->SetTitleSize(0.05);
    effRsqTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetXaxis()->SetTitleOffset(0.75);
    effRsqTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetYaxis()->SetTitleOffset(0.90);
    vline = new TLine(0.25,0.8,0.25,1.05);
    vline->SetLineStyle(2);
    vline->SetLineWidth(3);
    vline->SetLineColor(kBlue);
    hline = new TLine(0.15,1.0,1.5,1.0);
    hline->SetLineStyle(3);
    hline->SetLineWidth(3);
    hline->SetLineColor(kBlack);
    vline->Draw();
    hline->Draw();
    cv->SaveAs("RazorTriggerEfficiencyVsRsq_RsqMR240_Rsq0p09_MR200_AllLeptonData_2016B.gif");
    cv->SaveAs("RazorTriggerEfficiencyVsRsq_RsqMR240_Rsq0p09_MR200_AllLeptonData_2016B.pdf");



    // cv = new TCanvas("cv","cv", 800,600);
    // cv->SetRightMargin(0.15);
    // effPtTrigger_RsqMR260_Rsq0p09_MR200_TTJets->SetStats(0);
    // effPtTrigger_RsqMR260_Rsq0p09_MR200_TTJets->SetTitle("; M_{R} [GeV/c^{2}; R^{2}; Efficiency");
    // effPtTrigger_RsqMR260_Rsq0p09_MR200_TTJets->GetYaxis()->SetRangeUser(0,0.8);

    // effPtTrigger_RsqMR260_Rsq0p09_MR200_TTJets->Draw("colz");

    // boxLabel = new TLatex();
    // boxLabel->SetNDC();
    // boxLabel->SetTextSize(0.040);
    // boxLabel->SetTextFont(42);
    // boxLabel->SetTextColor(kBlack);
    // boxLabel->DrawLatex(0.2,0.92,"Hadronic & Loose Lepton Categories ( t#bar{t} Monte Carlo )");
    // boxLabel->DrawLatex(0.4,0.85,"RsqMR260_Rsq0p09_MR200 Trigger");
    // boxLabel->Draw();

    // cv->SaveAs("RazorTriggerEfficiency_RsqMR260_Rsq0p09_MR200_TTJetsHadronicBoxes_MC.gif");
    // cv->SaveAs("RazorTriggerEfficiency_RsqMR260_Rsq0p09_MR200_TTJetsHadronicBoxes_MC.pdf");


    // cv = new TCanvas("cv","cv", 800,600);
    // cv->SetRightMargin(0.15);
    // effPtTrigger_RsqMR300_Rsq0p09_MR200_TTJets->SetStats(0);
    // effPtTrigger_RsqMR300_Rsq0p09_MR200_TTJets->SetTitle("; M_{R} [GeV/c^{2}; R^{2}; Efficiency");
    // effPtTrigger_RsqMR300_Rsq0p09_MR200_TTJets->Draw("colz");

    // boxLabel = new TLatex();
    // boxLabel->SetNDC();
    // boxLabel->SetTextSize(0.040);
    // boxLabel->SetTextFont(42);
    // boxLabel->SetTextColor(kBlack);
    // boxLabel->DrawLatex(0.2,0.92,"Hadronic & Loose Lepton Categories ( t#bar{t} Monte Carlo )");
    // boxLabel->DrawLatex(0.4,0.85,"RsqMR300_Rsq0p09_MR200 Trigger");
    // boxLabel->Draw();

    // cv->SaveAs("RazorTriggerEfficiency_RsqMR300_Rsq0p09_MR200_TTJetsHadronicBoxes_MC.gif");
    // cv->SaveAs("RazorTriggerEfficiency_RsqMR300_Rsq0p09_MR200_TTJetsHadronicBoxes_MC.pdf");


    // cv = new TCanvas("cv","cv", 800,600);
    // cv->SetRightMargin(0.15);
    // effPtTrigger_RsqMR260_Rsq0p09_MR200_WJets->SetStats(0);
    // effPtTrigger_RsqMR260_Rsq0p09_MR200_WJets->SetTitle("; M_{R} [GeV/c^{2}; R^{2}; Efficiency");
    // effPtTrigger_RsqMR260_Rsq0p09_MR200_WJets->Draw("colz");

    // boxLabel = new TLatex();
    // boxLabel->SetNDC();
    // boxLabel->SetTextSize(0.040);
    // boxLabel->SetTextFont(42);
    // boxLabel->SetTextColor(kBlack);
    // boxLabel->DrawLatex(0.2,0.92,"1-Lepton Categories ( W+Jets Monte Carlo )");
    // boxLabel->DrawLatex(0.4,0.85,"RsqMR260_Rsq0p09_MR200 Trigger");
    // boxLabel->Draw();

    // cv->SaveAs("RazorTriggerEfficiency_RsqMR260_Rsq0p09_MR200_WJetsOneLeptonBoxes_MC.gif");
    // cv->SaveAs("RazorTriggerEfficiency_RsqMR260_Rsq0p09_MR200_WJetsOneLeptonBoxes_MC.pdf");


    // cv = new TCanvas("cv","cv", 800,600);
    // cv->SetRightMargin(0.15);
    // effPtTrigger_RsqMR300_Rsq0p09_MR200_WJets->SetStats(0);
    // effPtTrigger_RsqMR300_Rsq0p09_MR200_WJets->SetTitle("; M_{R} [GeV/c^{2}; R^{2}; Efficiency");
    // effPtTrigger_RsqMR300_Rsq0p09_MR200_WJets->Draw("colz");

    // boxLabel = new TLatex();
    // boxLabel->SetNDC();
    // boxLabel->SetTextSize(0.040);
    // boxLabel->SetTextFont(42);
    // boxLabel->SetTextColor(kBlack);
    // boxLabel->DrawLatex(0.2,0.92,"1-Lepton Lepton Categories ( W+Jets Monte Carlo )");
    // boxLabel->DrawLatex(0.4,0.85,"RsqMR300_Rsq0p09_MR200 Trigger");
    // boxLabel->Draw();

    // cv->SaveAs("RazorTriggerEfficiency_RsqMR300_Rsq0p09_MR200_WJetsOneLeptonBoxes_MC.gif");
    // cv->SaveAs("RazorTriggerEfficiency_RsqMR300_Rsq0p09_MR200_WJetsOneLeptonBoxes_MC.pdf");


    //*******************************************************************************  
    // Data Efficiency Plots
    //*******************************************************************************


    // cv = new TCanvas("cv","cv", 800,600);
    // cv->SetRightMargin(0.15);
    // effPtTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->SetStats(0);
    // effPtTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->SetTitle("; M_{R} [GeV/c^{2}]; R^{2}; Efficiency");
    // //effPtTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->GetYaxis()->SetRangeUser(0,0.8);
    // effPtTrigger_RsqMR240_Rsq0p09_MR200_AllLeptonData->Draw("colz");

    // boxLabel = new TLatex();
    // boxLabel->SetNDC();
    // boxLabel->SetTextSize(0.040);
    // boxLabel->SetTextFont(42);
    // boxLabel->SetTextColor(kBlack);
    // boxLabel->DrawLatex(0.2,0.92,"1-Lepton Categories ( Lepton Trigger Data )");
    // boxLabel->DrawLatex(0.4,0.85,"RsqMR240_Rsq0p09_MR200 Trigger");
    // boxLabel->Draw();

    // cv->SaveAs("RazorTriggerEfficiency_RsqMR240_Rsq0p09_MR200_AllLeptonData_2015B.gif");
    // cv->SaveAs("RazorTriggerEfficiency_RsqMR240_Rsq0p09_MR200_AllLeptonData_2015B.pdf");


    // cv = new TCanvas("cv","cv", 800,600);
    // cv->SetRightMargin(0.15);
    // effPtTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData->SetStats(0);
    // effPtTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData->SetTitle("; M_{R} [GeV/c^{2}; R^{2}; Efficiency");
    // effPtTrigger_RsqMR270_Rsq0p09_MR200_AllLeptonData->Draw("colz");

    // boxLabel = new TLatex();
    // boxLabel->SetNDC();
    // boxLabel->SetTextSize(0.040);
    // boxLabel->SetTextFont(42);
    // boxLabel->SetTextColor(kBlack);
    // boxLabel->DrawLatex(0.2,0.92,"1-Lepton Lepton Categories ( Lepton Trigger Data )");
    // boxLabel->DrawLatex(0.4,0.85,"RsqMR270_Rsq0p09_MR200 Trigger");
    // boxLabel->Draw();

    // cv->SaveAs("RazorTriggerEfficiency_RsqMR270_Rsq0p09_MR200_AllLeptonData_2015B.gif");
    // cv->SaveAs("RazorTriggerEfficiency_RsqMR270_Rsq0p09_MR200_AllLeptonData_2015B.pdf");


}



//=== MAIN MACRO ================================================================================================= 

void ProduceRazorTriggerEfficiencyPlots(const string inputfile, int year, int wp, int option = -1, string label = "") {

  cout << "TEST\n\n";

    string Label = "";
    if (label != "") Label = "_" + label;

    //--------------------------------------------------------------------------------------------------------------
    // Settings 
    //============================================================================================================== 
    bool printdebug = false;


    //*****************************************************************************************
    //Make some histograms
    //*****************************************************************************************
    TH1F *histDenominatorMR = new TH1F ("histDenominatorMR",";M_{R} [GeV/c^{2}]; Number of Events", 1000, 0, 2000);
    TH1F *histNumeratorMR = new TH1F ("histNumeratorMR",";M_{R} [GeV/c^{2}]; Number of Events", 1000, 0 , 2000);
    TH1F *histDenominatorRsq = new TH1F ("histDenominatorRsq",";R^{2}; Number of Events", 400, 0 , 2.0);
    TH1F *histNumeratorRsq = new TH1F ("histNumeratorRsq",";R^{2}; Number of Events", 400, 0 , 2.0);

    // TH1F *histDenominatorMR = new TH1F ("histDenominatorMR",";M_{R} [GeV/c^{2}]; Number of Events", 40, 0, 1000);
    // TH1F *histNumeratorMR = new TH1F ("histNumeratorMR",";M_{R} [GeV/c^{2}]; Number of Events", 40, 0 , 1000);
    // TH1F *histDenominatorRsq = new TH1F ("histDenominatorRsq",";R^{2}; Number of Events", 30, 0 , 1.5);
    // TH1F *histNumeratorRsq = new TH1F ("histNumeratorRsq",";R^{2}; Number of Events", 30, 0 , 1.5);

    TH2F *histDenominatorMRRsq = 0;
    TH2F *histNumeratorMRRsq = 0;
    // histDenominatorMRRsq = new TH2F ("histDenominatorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 1000, 0 , 1000, 400, 0, 10);
    // histNumeratorMRRsq = new TH2F ("histNumeratorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 1000, 0 , 1000, 400, 0, 10);

    // // if (option >= 10) {
    histDenominatorMRRsq = new TH2F ("histDenominatorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 100, 150 , 1900, 100, 0.0, 2.5);
    histNumeratorMRRsq = new TH2F ("histNumeratorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 100, 150 , 1900, 100, 0.0, 2.5);
    // // } else {
    // histDenominatorMRRsq = new TH2F ("histDenominatorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 1000, 0 , 1000, 400, 0, 10);
    // histNumeratorMRRsq = new TH2F ("histNumeratorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 1000, 0 , 1000, 400, 0, 10);
    //histDenominatorMRRsq = new TH2F ("histDenominatorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 40, 0 , 1000, 30, 0, 1.5);
    //histNumeratorMRRsq = new TH2F ("histNumeratorMRRsq",";M_{R} [GeV/c^{2}] ; R^{2} ; Number of Events", 40, 0 , 1000, 30, 0, 1.5);

    // vector<double> MyMRBins;
    // vector<double> MyRsqBins;
    // // MyMRBins.push_back(200);
    // // MyMRBins.push_back(225); 
    // // MyMRBins.push_back(250); 
    // // MyMRBins.push_back(275); 
    // // MyMRBins.push_back(300); 
    // // MyMRBins.push_back(325); 
    // // MyMRBins.push_back(350); 
    // // MyMRBins.push_back(375); 
    // MyMRBins.push_back(400); 
    // MyMRBins.push_back(425); 
    // MyMRBins.push_back(450); 
    // MyMRBins.push_back(500); 
    // MyMRBins.push_back(550); 
    // MyMRBins.push_back(600); 
    // MyMRBins.push_back(800); 
    // MyMRBins.push_back(1000); 
    // // MyRsqBins.push_back(0.10);
    // // MyRsqBins.push_back(0.15);
    // // MyRsqBins.push_back(0.20);
    // MyRsqBins.push_back(0.25);
    // MyRsqBins.push_back(0.30);
    // MyRsqBins.push_back(0.40);
    // MyRsqBins.push_back(0.60);
    // MyRsqBins.push_back(10);

    // histNumeratorMRRsq = rebin(histNumeratorMRRsq, MyMRBins, MyRsqBins);
    // histDenominatorMRRsq = rebin(histDenominatorMRRsq, MyMRBins, MyRsqBins);


    // }

    //*******************************************************************************************
    //Read file
    //*******************************************************************************************                
    TFile* inputFile = TFile::Open(inputfile.c_str(),"READ");
    assert(inputFile);
    inputFile->ls();
    TTree* tree = (TTree*)inputFile->Get("RazorInclusive");
    //TTree* tree = (TTree*)inputFile->Get("ControlSampleEvent");
    assert(tree);

    float weight = 0;
    int nvtx = 0;
    int box = -1;
    int nBTaggedJets = 0;
    int nSelectedJets = 0;
    int NJets80 = 0;
    float dPhiRazor = 0;
    float MR = 0;
    float Rsq = 0;
    float MET = 0;
    bool  HLTDecision[300];
    UInt_t run = 0;
    UInt_t lumi = 0;
    UInt_t event = 0;
    bool Flag_HBHENoiseFilter = false;
    bool Flag_goodVertices = false;
    bool Flag_eeBadScFilter = false;
    float leadingMuonPt = 0;
    float allMuonPt = 0;
    float HLTMR = 0;
    float HLTRSQ = 0;
    float leadingJetPt = 0;
    float subleadingJetPt = 0;

    tree->SetBranchAddress("weight",&weight);
    tree->SetBranchAddress("box",&box);
    tree->SetBranchAddress("nVtx",&nvtx);
    tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
    tree->SetBranchAddress("nSelectedJets",&nSelectedJets);
    tree->SetBranchAddress("nJets80",&NJets80);
    tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
    tree->SetBranchAddress("MR",&MR);
    tree->SetBranchAddress("Rsq",&Rsq);
    tree->SetBranchAddress("met",&MET);
    tree->SetBranchAddress("HLTDecision",&HLTDecision);
    tree->SetBranchAddress("run",&run);
    tree->SetBranchAddress("lumi",&lumi);
    tree->SetBranchAddress("event",&event);
    tree->SetBranchAddress("Flag_HBHENoiseFilter",&Flag_HBHENoiseFilter);
    tree->SetBranchAddress("Flag_goodVertices",&Flag_goodVertices);
    tree->SetBranchAddress("Flag_eeBadScFilter",&Flag_eeBadScFilter);
    tree->SetBranchAddress("HLTMR",&HLTMR);
    tree->SetBranchAddress("HLTRSQ",&HLTRSQ);
    tree->SetBranchAddress("leadingJetPt",&leadingJetPt);
    tree->SetBranchAddress("subleadingJetPt",&subleadingJetPt);
    //  tree->SetBranchAddress("leadingMuonPt",&leadingMuonPt);
    //  tree->SetBranchAddress("allMuonPt",&allMuonPt);

    cout << "Total Entries: " << tree->GetEntries() << "\n";

    cout << "TEST\n\n\n\nTEST\n";

    //for duplicate event checking
    //  map<pair<uint,uint>, bool > processedRunEvents;

    for(UInt_t ientry=0; ientry < tree->GetEntries(); ientry++) {
        tree->GetEntry(ientry);
        if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
        //Cuts
        if (!(NJets80 >= 2)) continue;
        if (!(leadingJetPt>JET_CUT && subleadingJetPt>JET_CUT)) continue;

       
        //Select Electron Triggered Events
        if (option == 1 || option == 11) {
            
	  bool passedSingleLeptonTrigger = false;

	  if (year == 2016) {
	    passedSingleLeptonTrigger = bool(HLTDecision[4] || HLTDecision[13] || HLTDecision[18] || HLTDecision[20] 
		 || HLTDecision[24] || HLTDecision[29] || HLTDecision[34] || HLTDecision[36] || HLTDecision[37]
		 || HLTDecision[38] || HLTDecision[39] || HLTDecision[42] || HLTDecision[42]);
	  } else if (year == 2017) {
	    passedSingleLeptonTrigger = bool(HLTDecision[12] || HLTDecision[13] || HLTDecision[14] || HLTDecision[15]);
	  }
	  if (!(passedSingleLeptonTrigger)) continue;
        
	  //Only Select Electron Boxes
	    if (!(box == 6 || box == 7 || box ==  8)) continue;

        }

       
        //Cleaning Cuts
        if (!(Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter)) continue;
        //    if (fabs(dPhiRazor) > 2.8) continue;


        //**** MR - Rsq ****
        //if (!(MR > 150 && Rsq > 0.15)) continue;
        //if (!(MR > 600)) continue;
        histDenominatorMRRsq->Fill(MR,Rsq);
        if(PassSelection(HLTDecision,year,wp,HLTMR, HLTRSQ)) {
            // cout << MR << " " << Rsq << " " << weight << "\n";
            histNumeratorMRRsq->Fill(MR,Rsq);
        }
      

        //**** MR ****
        if (Rsq>0.4) 
        { 
            histDenominatorMR->Fill(MR);
            //Numerator
            if(PassSelection(HLTDecision,year,wp,HLTMR, HLTRSQ)) {
                histNumeratorMR->Fill(MR);
            } else {
                if (MR > 700) {
                    //cout << "Fail Event: " << run << " " << lumi << " " << event << " : " << MR << " " << Rsq << " " << bool(Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter) << "\n";
                }
            }
        }

        //**** Rsq ****
      if (MR>600) 
        {
            histDenominatorRsq->Fill(Rsq);      
            //Numerator
            if(PassSelection(HLTDecision,year,wp,HLTMR, HLTRSQ)) {
                histNumeratorRsq->Fill(Rsq);
            } else {
                if (Rsq > 0.6) { 
                    //  cout << "Fail Event: " << run << " " << lumi << " " << event << " : " << MR << " " << Rsq << " " << bool(Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter) << "\n";
                }
            }
        }


    }

    //--------------------------------------------------------------------------------------------------------------
    // Make Efficiency Plots
    //==============================================================================================================
    vector<double> MRBins; 
    MRBins.push_back(150);
    MRBins.push_back(175);
    MRBins.push_back(200);
    MRBins.push_back(225);
    MRBins.push_back(250);
    MRBins.push_back(275);
    MRBins.push_back(300);
    MRBins.push_back(325);
    MRBins.push_back(350);
    MRBins.push_back(375);
    MRBins.push_back(400);
    MRBins.push_back(450);
    MRBins.push_back(500);
    MRBins.push_back(600);
    MRBins.push_back(650);
    MRBins.push_back(700);
    MRBins.push_back(750);
    MRBins.push_back(800);
    MRBins.push_back(850);
    MRBins.push_back(900);
    MRBins.push_back(950);
    MRBins.push_back(1000);
    // MRBins.push_back(1050);
    // MRBins.push_back(1100);
    // MRBins.push_back(1200);
    // MRBins.push_back(1300);
    // MRBins.push_back(1400);
    // MRBins.push_back(1600);
    // MRBins.push_back(1900);
    vector<double> RsqBins; 
    RsqBins.push_back(0.0);
    RsqBins.push_back(0.025);
    RsqBins.push_back(0.05);
    RsqBins.push_back(0.075);
    RsqBins.push_back(0.10);
    RsqBins.push_back(0.125);
    RsqBins.push_back(0.15);
    RsqBins.push_back(0.20);
    RsqBins.push_back(0.30);
    RsqBins.push_back(0.40);
    RsqBins.push_back(0.50);
    RsqBins.push_back(0.60);
    RsqBins.push_back(0.70);
    RsqBins.push_back(0.80);
    RsqBins.push_back(0.90);
    RsqBins.push_back(1.0);
    RsqBins.push_back(1.2);
    RsqBins.push_back(2.0);

    vector<double> MR2DBins; 
    //  MR2DBins.push_back(600);
    //  MR2DBins.push_back(650);
    //  MR2DBins.push_back(700);
    //  MR2DBins.push_back(750);
    //  MR2DBins.push_back(800);
    //  MR2DBins.push_back(850);
    //  MR2DBins.push_back(900);
    //  MR2DBins.push_back(950);
    //  MR2DBins.push_back(1000);
    //  MR2DBins.push_back(1050);
    //  MR2DBins.push_back(1100);
    //  MR2DBins.push_back(1200);
    //  MR2DBins.push_back(1300);
    //  MR2DBins.push_back(1400);
    //  MR2DBins.push_back(1600);
    //  MR2DBins.push_back(1900);
    MR2DBins.push_back(150);
    MR2DBins.push_back(175);
    MR2DBins.push_back(200);
    MR2DBins.push_back(225);
    MR2DBins.push_back(250);
    MR2DBins.push_back(275);
    MR2DBins.push_back(300);
    MR2DBins.push_back(325);
    MR2DBins.push_back(350);
    MR2DBins.push_back(375);
    MR2DBins.push_back(400);
    MR2DBins.push_back(450);
    MR2DBins.push_back(500);
    MR2DBins.push_back(600);
    MR2DBins.push_back(700);
    MR2DBins.push_back(900);
    MR2DBins.push_back(1200);
    // vector<double> MRBins; 
    // MRBins.push_back(500);
    // vector<double> RsqBins; 
    // RsqBins.push_back(0.25);

    // TGraphAsymmErrors *efficiency_MR = createEfficiencyGraph(histNumeratorMR, histDenominatorMR, "Efficiency_MR" , vector<double>() ,  -99, -99, 0.9, 1);
    // TGraphAsymmErrors5*efficiency_Rsq = createEfficiencyGraph(histNumeratorRsq, histDenominatorRsq, "Efficiency_Rsq" , vector<double>() ,  -99, -99, 0.9, 1);
    TH2F *efficiency_MRRsq = createEfficiencyHist2D(histNumeratorMRRsq, histDenominatorMRRsq, "Efficiency_MRRsq" , MR2DBins ,RsqBins);  
    TGraphAsymmErrors *efficiency_MR = createEfficiencyGraph(histNumeratorMR, histDenominatorMR, "Efficiency_MR" , MRBins ,  -99, -99, 0.9, 1);
    TGraphAsymmErrors *efficiency_Rsq = createEfficiencyGraph(histNumeratorRsq, histDenominatorRsq, "Efficiency_Rsq" , RsqBins ,  -99, -99, 0.9, 1);
    //TH2F *efficiency_MRRsq = createEfficiencyHist2D(histNumeratorMRRsq, histDenominatorMRRsq, "Efficiency_MRRsq" , MRBins , RsqBins);  


    //--------------------------------------------------------------------------------------------------------------
    // Draw
    //==============================================================================================================
    TCanvas *cv =0;

    cv = new TCanvas("cv","cv",800,600);
    gPad->SetGridx();
    efficiency_MR->Draw("AP");
    efficiency_MR->SetTitle("");
    efficiency_MR->GetYaxis()->SetRangeUser(0.0,1.1);
    efficiency_MR->GetXaxis()->SetTitle("M_{R} [GeV/c^{2}]");
    efficiency_MR->GetYaxis()->SetTitle("Efficiency");
    efficiency_MR->GetYaxis()->SetTitleOffset(1.2);
    efficiency_MR->SetLineWidth(3);  
    efficiency_MR->SetMarkerSize(2);
    AddCMS(cv);
    cv->SaveAs(("Efficiency"+Label+"_MR.png").c_str());
    cv->SaveAs(("Efficiency"+Label+"_MR.C").c_str());

    cv = new TCanvas("cv","cv",800,600);
    gPad->SetGridx();
    efficiency_Rsq->Draw("AP");
    efficiency_Rsq->SetTitle("");
    efficiency_Rsq->GetXaxis()->SetRangeUser(0.0,1.1);
    efficiency_Rsq->GetYaxis()->SetRangeUser(0.0,1.1);
    efficiency_Rsq->GetXaxis()->SetTitle("R^{2}");
    efficiency_Rsq->GetYaxis()->SetTitle("Efficiency");
    efficiency_Rsq->GetYaxis()->SetTitleOffset(1.2);
    efficiency_Rsq->SetLineWidth(3);  
    efficiency_Rsq->SetMarkerSize(2);
    AddCMS(cv);
    cv->SaveAs(("Efficiency"+Label+"_Rsq.png").c_str());
    cv->SaveAs(("Efficiency"+Label+"_Rsq.C").c_str());

    cv = new TCanvas("cv","cv",800,600);
    //gPad->SetLogx();
    //  gPad->SetLogz();
    gPad->SetRightMargin(0.12);
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("0.1e");
    efficiency_MRRsq->SetMarkerSize(1.2);
    efficiency_MRRsq->SetMaximum(1);
    efficiency_MRRsq->SetMinimum(0);
    efficiency_MRRsq->Draw("COLZ");
    efficiency_MRRsq->SetTitle("");
    efficiency_MRRsq->GetXaxis()->SetTitle("M_{R}");
    efficiency_MRRsq->GetYaxis()->SetTitle("R^{2}");
    efficiency_MRRsq->GetZaxis()->SetTitle("Efficiency");
    efficiency_MRRsq->GetYaxis()->SetTitleOffset(1.2);
    efficiency_MRRsq->GetYaxis()->SetRangeUser(0.0,1.2);
    efficiency_MRRsq->SetLineWidth(3); 
    AddCMS(cv);
    cv->SaveAs(("Efficiency"+Label+"_MRRsq.png").c_str());
    cv->SaveAs(("Efficiency"+Label+"_MRRsq.C").c_str());


    //--------------------------------------------------------------------------------------------------------------
    // Output
    //==============================================================================================================
    TFile *file = TFile::Open(("Efficiency"+Label+".root").c_str(), "UPDATE");
    file->cd();
    file->WriteTObject(efficiency_MR, "Efficiency_MR", "WriteDelete");
    file->WriteTObject(efficiency_Rsq, "Efficiency_Rsq", "WriteDelete");
    file->WriteTObject(efficiency_MRRsq, "Efficiency_MRRsq", "WriteDelete");

    file->Close();
    //delete file;       

}

void MakeRazorTriggerEfficiencyPlots( int option = 0) {

    if (option == 1) {  
        //***************************************
        // TTbar MC : Use Hadronic Boxes
        //***************************************
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 1, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_TTJets_25ns");   
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_4jet_TTJets_25ns");   
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 3, 0, "RazorTrigger_RsqMR300_Rsq0p09_MR200_TTJets_25ns");   
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 4, 0, "RazorTrigger_RsqMR300_Rsq0p09_MR200_4jet_TTJets_25ns");   
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/MC/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root", 10, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_All_TTJets_25ns");   
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 11, 0, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_TTJets_25ns");   

        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/MC/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root", 10, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_All_TTJets_25ns");   

        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_25ns.root", 10, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_All_TTJets_25ns"); 
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 10, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_All_TTJets_25ns"); 
        //ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 10, 0, "RazorTrigger_RsqMR260_Rsq0p09_MR200_All_TTJets_25ns"); 

        //***************************************
        // W+Jets MC: Use Single Lepton Boxes
        //***************************************
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_WJetsToLNu_HTBinned_25ns_1pb_weighted.root", 1, 1, "RazorTrigger_RsqMR260_Rsq0p09_MR200_WJets_25ns");   
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_WJetsToLNu_HTBinned_25ns_1pb_weighted.root", 2, 1, "RazorTrigger_RsqMR260_Rsq0p09_MR200_4jet_WJets_25ns");   
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_WJetsToLNu_HTBinned_25ns_1pb_weighted.root", 3, 1, "RazorTrigger_RsqMR300_Rsq0p09_MR200_WJets_25ns");   
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_WJetsToLNu_HTBinned_25ns_1pb_weighted.root", 4, 1, "RazorTrigger_RsqMR300_Rsq0p09_MR200_4jet_WJets_25ns");   
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_WJetsToLNu_HTBinned_25ns_1pb_weighted.root", 10, 1, "RazorTrigger_RsqMR260_Rsq0p09_MR200_All_WJets_25ns");   
        // ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/Spring15_25ns/RazorInclusive_WJetsToLNu_HTBinned_25ns_1pb_weighted.root", 11, 1, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_WJets_25ns");   

    }

    if (option == 2) {

      //**************************************************************************
      //2015 Data    
      //**************************************************************************
      //ProduceRazorTriggerEfficiencyPlots("/afs/cern.ch/user/q/qnguyen/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForARCReview20151129/RazorSkim/RazorInclusive_SingleLepton_Run2015D_2093pb_GoodLumiGolden_RazorSkim.root", 15, 11, "RazorTrigger_RsqMR240_Rsq0p09_MR200_All_AllLeptonData_2015D");

      //**************************************************************************
      //2016 Data (with Razor Hadronic Triggers)
      //**************************************************************************
      //ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_29Aug2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_Data_RazorSkim.root", 2016, 16, 11, "RazorTrigger_RsqMR270_All_SingleElectronData");

      //**************************************************************************
      //2016 Data (with Razor Hadronic Triggers OR PFHT triggers OR PFJet triggers)
      //**************************************************************************     
      //ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_03Oct2017_NoCuts/Signal/FullRazorInclusive_Razor2016_MoriondRereco_Data_NoDuplicates_GoodLumiGolden.root", 2016, 20, 11, "RazorTrigger_RsqMR270ORPFHT900ORPFJet450_All_SingleElectronData");

      //**************************************************************************
      //2017 Data (For Razor Hadronic Triggers
      //**************************************************************************
      // //wp=6 main trigger
/*       
       ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2017/V4p0_20171031/FullRazorInclusive_SingleElectron_Run2017C-PromptReco-v1_GoldenLumi.root", 2017, 6, 11, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_SingleElectronData_2017Cv1");
       ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2017/V4p0_20171031/FullRazorInclusive_SingleElectron_Run2017C-PromptReco-v2_GoldenLumi.root", 2017, 6, 11, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_SingleElectronData_2017Cv2");
       ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2017/V4p0_20171031/FullRazorInclusive_SingleElectron_Run2017C-PromptReco-v3_GoldenLumi.root", 2017, 6, 11, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_SingleElectronData_2017Cv3");
       ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2017/V4p0_20171031/FullRazorInclusive_SingleElectron_Run2017D-PromptReco-v1_GoldenLumi.root", 2017, 6, 11, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_SingleElectronData_2017Dv1");
       ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2017/V4p0_20171031/FullRazorInclusive_SingleElectron_Run2017E-PromptReco-v1_GoldenLumi.root", 2017, 6, 11, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_SingleElectronData_2017Ev1");
       ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2017/V4p0_20171031/FullRazorInclusive_SingleElectron_Run2017F-PromptReco-v1_GoldenLumi.root", 2017, 6, 11, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_SingleElectronData_2017Fv1");
       ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2017/V4p0_20171031/FullRazorInclusive_SingleElectron_Run2017A-F-PromptReco_GoldenLumi.root", 2017, 6, 11, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_SingleElectronData_2017A-F");
       
       ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2017/V4p0_20171031/FullRazorInclusive_SingleElectron_Run2017B-F-PromptReco_GoldenLumi.root", 2017, 6, 11, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_SingleElectronData_2017B-F");

*/

       //ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2017/V4p0_20171031/FullRazorInclusive_SingleElectron_Run2017C-F-PromptReco_GoldenLumi.root", 2017, 6, 11, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_SingleElectronData_2017C-F");
       //ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2017/V4p0_20171127/FullRazorInclusive_SingleElectron_Run2017C-F-PromptReco_GoldenLumi.root", 2017, 6, 11, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_SingleElectronData_2017C-F");
       //ProduceRazorTriggerEfficiencyPlots("/eos/cms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2018/V5p0_20180524//FullRazorInclusive_EGamma_2018A_PromptReco_GoodLumi.root", 2017, 6, 11, "RazorTrigger_RsqMR300_Rsq0p09_MR200_All_EGammaData_2018A");

       //wp=6 dijet||quadjet 
       //10.3 fb-1 
       //ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2018/V5p0_20180601/FullRazorInclusive_EGamma_2018A_v1v2_PromptReco_GoodLumi.root", 2017, 6, 11, "RazorTrigger_EGammaData_2018Av1v2_dijet_quadjet_JetPt120");
       ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2018/V5p0_20180810/FullRazorInclusive_EGamma_2018A_PromptReco_v6_GoodLumi.root", 2017, 6, 11, "RazorTrigger_EGammaData_2018A_dijet_quadjet_JetPt120");
       //24.3 fb-1 
       //ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2018/V5p0_20180810/FullRazorInclusive_EGamma_2018ABCD_PromptReco_v7_GoodLumi.root", 2017, 6, 11, "RazorTrigger_EGammaData_2018ABCD_dijet_quadjet_JetPt120");
       //ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2018/V5p0_20180810/FullRazorInclusive_EGamma_2018BCD_PromptReco_v7_GoodLumi.root", 2017, 6, 11, "RazorTrigger_EGammaData_2018BCD_dijet_quadjet_JetPt120");
       //6.8 fb-1 
       //ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2018/V5p0_20180810/FullRazorInclusive_EGamma_2018B_PromptReco_v7_GoodLumi.root", 2017, 6, 11, "RazorTrigger_EGammaData_2018B_dijet_quadjet_JetPt120");
       //3.8 fb-1 
       //ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2018/V5p0_20180810/FullRazorInclusive_EGamma_2018C_PromptReco_v7_GoodLumi.root", 2017, 6, 11, "RazorTrigger_EGammaData_2018C_dijet_quadjet_JetPt120");
       //ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2018/V5p0_20180810/FullRazorInclusive_EGamma_2018D_PromptReco_v7_GoodLumi.root", 2017, 6, 11, "RazorTrigger_EGammaData_2018D_dijet_quadjet_JetPt120");
       //wp=9 dijet||quadjet||rsq-only 
       //ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2018/V5p0_20180601/FullRazorInclusive_EGamma_2018A_v1v2_PromptReco_GoodLumi.root", 2017, 9, 11, "RazorTrigger_EGammaData_2018Av1v2_dijet_qudjet_rsq-onlyJetPt120");
       
      // //wp=7 backup trigger
       //ProduceRazorTriggerEfficiencyPlots("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2017/V4p0_20171031/FullRazorInclusive_SingleElectron_Run2017A-F-PromptReco_GoldenLumi.root", 2017, 7, 11, "RazorTrigger_RsqMR320_Rsq0p09_MR200_All_SingleElectronData_2017A-F");


    }
   

    if (option == 0) {
        plotRazorTriggerEfficiency();
    }

}


bool AddCMS( TCanvas* C )
{
  C->cd();
  float lumix = 0.925;
  float lumiy = 0.945;
  float lumifont = 42;
  
  float cmsx = 0.225;
  float cmsy = 0.940;
  float cmsTextFont   = 61;  // default is helvetic-bold
  float extrax = cmsx + 0.198;
  float extray = cmsy;
  float extraTextFont = 52;  // default is helvetica-italics
  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;
  float cmsSize = 0.06;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    
  float extraTextSize = extraOverCmsTextSize*cmsSize;
  latex.SetTextFont(lumifont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(cmsSize);    
  latex.DrawLatex(lumix, lumiy,lumiText);

  latex.SetTextFont(cmsTextFont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(cmsSize);
  latex.DrawLatex(cmsx, cmsy, CMSText);
   
  latex.SetTextFont(extraTextFont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(extraTextSize);
  latex.DrawLatex(extrax, extray, extraText);
  return true;
};
