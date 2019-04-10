//================================================================================================
//
// Simple Example
//
//root -l RazorAnalyzer/macros/ObjectStudies/MakeBTaggingEfficiencyPlots.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/JetNtuple/JetNtuple_Prompt_TTJets_25ns.root",-1,"BTagging")'
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
#include <TH1F.h>                
#include <TCanvas.h>                
#include <TLegend.h>                
#include <TGraphAsymmErrors.h>                

#include "RazorAnalyzer/macros/ObjectStudies/EfficiencyUtils.hh"
#include "RazorAnalyzer/include/JetTree.h"

#endif


bool PassSelection( JetTree* JetTree ) {

  bool pass = false;

  //Medium WP
  if (JetTree->fJetCISV > 0.890) {
    pass = true;
  }

  return pass;
}


void plotBTaggingEfficiency() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileCSVMediumBJets = new TFile("Efficiency_BJets_25ns.root","READ");
  TFile *fileCSVMediumCJets = new TFile("Efficiency_CharmJets_25ns.root","READ");
  TFile *fileCSVMediumLightJets = new TFile("Efficiency_LightJets_25ns.root","READ");
 
  TGraphAsymmErrors* effPtCSVMediumBJets = (TGraphAsymmErrors*)fileCSVMediumBJets->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtCSVMediumCJets = (TGraphAsymmErrors*)fileCSVMediumCJets->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtCSVMediumLightJets = (TGraphAsymmErrors*)fileCSVMediumLightJets->Get("Efficiency_Pt");
  TGraphAsymmErrors* effEtaCSVMediumBJets = (TGraphAsymmErrors*)fileCSVMediumBJets->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaCSVMediumCJets = (TGraphAsymmErrors*)fileCSVMediumCJets->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaCSVMediumLightJets = (TGraphAsymmErrors*)fileCSVMediumLightJets->Get("Efficiency_Eta");
  TGraphAsymmErrors* effNPVCSVMediumBJets = (TGraphAsymmErrors*)fileCSVMediumBJets->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNPVCSVMediumCJets = (TGraphAsymmErrors*)fileCSVMediumCJets->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNPVCSVMediumLightJets = (TGraphAsymmErrors*)fileCSVMediumLightJets->Get("Efficiency_NPV");


  //*********************************************************************
  //B-tag efficiency Vs Pt
  //*********************************************************************
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtCSVMediumBJets, "B-Jets", "LP");

  effPtCSVMediumBJets->SetLineWidth(3);
  effPtCSVMediumBJets->SetLineColor(kBlack);
  effPtCSVMediumBJets->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
  effPtCSVMediumBJets->GetYaxis()->SetTitle("Selection Efficiency");
  effPtCSVMediumBJets->GetYaxis()->SetTitleOffset(1.2);

  effPtCSVMediumBJets->Draw("AP");
  
  legend->Draw();

  cv->SaveAs("BTaggingEfficiencyVsPt.gif");
  cv->SaveAs("BTaggingEfficiencyVsPt.pdf");


  //*********************************************************************
  //B-tag efficiency Vs Eta
  //*********************************************************************
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effEtaCSVMediumBJets, "B-Jets", "LP");

  effEtaCSVMediumBJets->SetLineWidth(3);
  effEtaCSVMediumBJets->SetLineColor(kBlack);
  effEtaCSVMediumBJets->GetXaxis()->SetTitle("Jet #eta");
  effEtaCSVMediumBJets->GetYaxis()->SetTitle("Selection Efficiency");
  effEtaCSVMediumBJets->GetYaxis()->SetTitleOffset(1.2);

  effEtaCSVMediumBJets->Draw("AP");
  
  legend->Draw();

  cv->SaveAs("BTaggingEfficiencyVsEta.gif");
  cv->SaveAs("BTaggingEfficiencyVsEta.pdf");

  //*********************************************************************
  //B-tag efficiency Vs NPV
  //*********************************************************************
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effNPVCSVMediumBJets, "B-Jets", "LP");

  effNPVCSVMediumBJets->SetLineWidth(3);
  effNPVCSVMediumBJets->SetLineColor(kBlack);
  effNPVCSVMediumBJets->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
  effNPVCSVMediumBJets->GetYaxis()->SetTitle("Selection Efficiency");
  effNPVCSVMediumBJets->GetYaxis()->SetTitleOffset(1.2);
  effNPVCSVMediumBJets->GetXaxis()->SetRangeUser(5,35);

  effNPVCSVMediumBJets->Draw("AP");
  
  legend->Draw();

  cv->SaveAs("BTaggingEfficiencyVsNPV.gif");
  cv->SaveAs("BTaggingEfficiencyVsNPV.pdf");



  //*********************************************************************
  //Mis Tag efficiency Vs Pt
  //*********************************************************************


  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtCSVMediumCJets, "Charm Jets", "LP");
  legend->AddEntry(effPtCSVMediumLightJets, "Light Jets", "LP");

  effPtCSVMediumCJets->SetLineWidth(3);
  effPtCSVMediumCJets->SetLineColor(kBlue);
  effPtCSVMediumCJets->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
  effPtCSVMediumCJets->GetYaxis()->SetTitle("Selection Efficiency");
  effPtCSVMediumCJets->GetYaxis()->SetTitleOffset(1.2);
  effPtCSVMediumCJets->GetYaxis()->SetRangeUser(0,0.3);

  effPtCSVMediumLightJets->SetLineWidth(3);
  effPtCSVMediumLightJets->SetLineColor(kRed);

  effPtCSVMediumCJets->Draw("AP");
  effPtCSVMediumLightJets->Draw("Psame");
  
  legend->Draw();
  cv->SaveAs("BTaggingMistagEfficiencyVsPt.gif");
  cv->SaveAs("BTaggingMistagEfficiencyVsPt.pdf");

  //*********************************************************************
  //Mis Tag efficiency Vs Eta
  //*********************************************************************


  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effEtaCSVMediumCJets, "Charm Jets", "LP");
  legend->AddEntry(effEtaCSVMediumLightJets, "Light Jets", "LP");

  effEtaCSVMediumCJets->SetLineWidth(3);
  effEtaCSVMediumCJets->SetLineColor(kBlue);
  effEtaCSVMediumCJets->GetXaxis()->SetTitle("Jet #eta");
  effEtaCSVMediumCJets->GetYaxis()->SetTitle("Selection Efficiency");
  effEtaCSVMediumCJets->GetYaxis()->SetTitleOffset(1.2);
  effEtaCSVMediumCJets->GetYaxis()->SetRangeUser(0,0.3);

  effEtaCSVMediumLightJets->SetLineWidth(3);
  effEtaCSVMediumLightJets->SetLineColor(kRed);

  effEtaCSVMediumCJets->Draw("AP");
  effEtaCSVMediumLightJets->Draw("Psame");
  
  legend->Draw();
  cv->SaveAs("BTaggingMistagEfficiencyVsEta.gif");
  cv->SaveAs("BTaggingMistagEfficiencyVsEta.pdf");

  //*********************************************************************
  //Mis Tag efficiency Vs NPV
  //*********************************************************************


  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effNPVCSVMediumCJets, "Charm Jets", "LP");
  legend->AddEntry(effNPVCSVMediumLightJets, "Light Jets", "LP");

  effNPVCSVMediumCJets->SetLineWidth(3);
  effNPVCSVMediumCJets->SetLineColor(kBlue);
  effNPVCSVMediumCJets->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
  effNPVCSVMediumCJets->GetYaxis()->SetTitle("Selection Efficiency");
  effNPVCSVMediumCJets->GetYaxis()->SetTitleOffset(1.2);
  effNPVCSVMediumCJets->GetYaxis()->SetRangeUser(0,0.3);
  effNPVCSVMediumCJets->GetXaxis()->SetRangeUser(5,35);

  effNPVCSVMediumLightJets->SetLineWidth(3);
  effNPVCSVMediumLightJets->SetLineColor(kRed);

  effNPVCSVMediumCJets->Draw("AP");
  effNPVCSVMediumLightJets->Draw("Psame");
  
  legend->Draw();
  cv->SaveAs("BTaggingMistagEfficiencyVsNPV.gif");
  cv->SaveAs("BTaggingMistagEfficiencyVsNPV.pdf");



  return;

 //  cv = new TCanvas("cv","cv", 800,600);

 //  legend = new TLegend(0.50,0.34,0.90,0.54);
 //  legend->SetTextSize(0.03);
 //  legend->SetBorderSize(0);
 //  legend->SetFillStyle(0);
 //  legend->AddEntry(effEtaVeto, "Veto", "LP");
 //  legend->AddEntry(effEtaLoose, "Loose", "LP");
 //  legend->AddEntry(effEtaTight, "Tight", "LP");

 //  effEtaVeto->SetLineWidth(3);
 //  effEtaVeto->SetLineColor(kBlack);
 //  effEtaVeto->GetXaxis()->SetTitle("BTagging #eta");
 //  effEtaVeto->GetYaxis()->SetTitle("Selection Efficiency");
 //  effEtaVeto->GetYaxis()->SetTitleOffset(1.2);

 //  effEtaLoose->SetLineWidth(3);
 //  effEtaLoose->SetLineColor(kBlue);
 //  effEtaTight->SetLineWidth(3);
 //  effEtaTight->SetLineColor(kRed);

 //  effEtaVeto->Draw("AP");
 //  effEtaLoose->Draw("Psame");
 //  effEtaTight->Draw("Psame");
  
 //  legend->Draw();  
 //  cv->SaveAs("BTaggingSelectionEfficiencyVsEta.gif");
 //  cv->SaveAs("BTaggingSelectionEfficiencyVsEta.pdf");



 // cv = new TCanvas("cv","cv", 800,600);

 //  legend = new TLegend(0.50,0.75,0.90,0.90);
 //  legend->SetTextSize(0.03);
 //  legend->SetBorderSize(0);
 //  legend->SetFillStyle(0);
 //  legend->AddEntry(effNpvVeto, "Veto", "LP");
 //  legend->AddEntry(effNpvLoose, "Loose", "LP");
 //  legend->AddEntry(effNpvTight, "Tight", "LP");

 //  effNpvVeto->SetLineWidth(3);
 //  effNpvVeto->SetLineColor(kBlack);
 //  effNpvVeto->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
 //  effNpvVeto->GetYaxis()->SetTitle("Selection Efficiency");
 //  effNpvVeto->GetYaxis()->SetTitleOffset(1.2);
 //  effNpvVeto->GetXaxis()->SetRangeUser(5,35);
 //  effNpvVeto->GetYaxis()->SetRangeUser(0.5,1.0);

 //  effNpvLoose->SetLineWidth(3);
 //  effNpvLoose->SetLineColor(kBlue);
 //  effNpvTight->SetLineWidth(3);
 //  effNpvTight->SetLineColor(kRed);

 //  effNpvVeto->Draw("AP");
 //  effNpvLoose->Draw("Psame");
 //  effNpvTight->Draw("Psame");
  
 //  legend->Draw();  
 //  cv->SaveAs("BTaggingSelectionEfficiencyVsNpv.gif");
 //  cv->SaveAs("BTaggingSelectionEfficiencyVsNpv.pdf");






 //  cv = new TCanvas("cv","cv", 800,600);

 //  legend = new TLegend(0.50,0.70,0.90,0.90);
 //  legend->SetTextSize(0.03);
 //  legend->SetBorderSize(0);
 //  legend->SetFillStyle(0);
 //  legend->AddEntry(effFakePtVeto, "Veto", "LP");
 //  legend->AddEntry(effFakePtLoose, "Loose", "LP");
 //  legend->AddEntry(effFakePtTight, "Tight", "LP");

 //  effFakePtVeto->SetLineWidth(3);
 //  effFakePtVeto->SetLineColor(kBlack);
 //  effFakePtVeto->GetXaxis()->SetTitle("BTagging p_{T} [GeV/c]");
 //  effFakePtVeto->GetYaxis()->SetTitle("Selection Efficiency");
 //  effFakePtVeto->GetYaxis()->SetTitleOffset(1.35);
 //  effFakePtVeto->GetYaxis()->SetRangeUser(0,0.20);

 //  effFakePtLoose->SetLineWidth(3);
 //  effFakePtLoose->SetLineColor(kBlue);
 //  effFakePtTight->SetLineWidth(3);
 //  effFakePtTight->SetLineColor(kRed);

 //  effFakePtVeto->Draw("AP");
 //  effFakePtLoose->Draw("Psame");
 //  effFakePtTight->Draw("Psame");
  
 //  legend->Draw();  
 //  cv->SaveAs("BTaggingSelectionFakeRateVsPt.gif");
 //  cv->SaveAs("BTaggingSelectionFakeRateVsPt.pdf");

 
 //  cv = new TCanvas("cv","cv", 800,600);

 //  legend = new TLegend(0.50,0.70,0.90,0.90);
 //  legend->SetTextSize(0.03);
 //  legend->SetBorderSize(0);
 //  legend->SetFillStyle(0);
 //  legend->AddEntry(effFakeEtaVeto, "Veto", "LP");
 //  legend->AddEntry(effFakeEtaLoose, "Loose", "LP");
 //  legend->AddEntry(effFakeEtaTight, "Tight", "LP");

 //  effFakeEtaVeto->SetLineWidth(3);
 //  effFakeEtaVeto->SetLineColor(kBlack);
 //  effFakeEtaVeto->GetXaxis()->SetTitle("BTagging #eta");
 //  effFakeEtaVeto->GetYaxis()->SetTitle("Selection Efficiency");
 //  effFakeEtaVeto->GetYaxis()->SetTitleOffset(1.35);
 //  effFakeEtaVeto->GetYaxis()->SetRangeUser(0,0.10);

 //  effFakeEtaLoose->SetLineWidth(3);
 //  effFakeEtaLoose->SetLineColor(kBlue);
 //  effFakeEtaTight->SetLineWidth(3);
 //  effFakeEtaTight->SetLineColor(kRed);

 //  effFakeEtaVeto->Draw("AP");
 //  effFakeEtaLoose->Draw("Psame");
 //  effFakeEtaTight->Draw("Psame");
  
 //  legend->Draw();  
 //  cv->SaveAs("BTaggingSelectionFakeRateVsEta.gif");
 //  cv->SaveAs("BTaggingSelectionFakeRateVsEta.pdf");

 


 // cv = new TCanvas("cv","cv", 800,600);

 //  legend = new TLegend(0.50,0.75,0.90,0.90);
 //  legend->SetTextSize(0.03);
 //  legend->SetBorderSize(0);
 //  legend->SetFillStyle(0);
 //  legend->AddEntry(effFakeNpvVeto, "Veto", "LP");
 //  legend->AddEntry(effFakeNpvLoose, "Loose", "LP");
 //  legend->AddEntry(effFakeNpvTight, "Tight", "LP");

 //  effFakeNpvVeto->SetLineWidth(3);
 //  effFakeNpvVeto->SetLineColor(kBlack);
 //  effFakeNpvVeto->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
 //  effFakeNpvVeto->GetYaxis()->SetTitle("Selection Efficiency");
 //  effFakeNpvVeto->GetYaxis()->SetTitleOffset(1.35);
 //  effFakeNpvVeto->GetXaxis()->SetRangeUser(5,35);
 //  effFakeNpvVeto->GetYaxis()->SetRangeUser(0,0.07);

 //  effFakeNpvLoose->SetLineWidth(3);
 //  effFakeNpvLoose->SetLineColor(kBlue);
 //  effFakeNpvTight->SetLineWidth(3);
 //  effFakeNpvTight->SetLineColor(kRed);

 //  effFakeNpvVeto->Draw("AP");
 //  effFakeNpvLoose->Draw("Psame");
 //  effFakeNpvTight->Draw("Psame");
  
 //  legend->Draw();  
 //  cv->SaveAs("BTaggingSelectionFakeRateVsNpv.gif");
 //  cv->SaveAs("BTaggingSelectionFakeRateVsNpv.pdf");




}




//=== MAIN MACRO ================================================================================================= 

void ProduceBTaggingEfficiencyPlots(const string inputfile, int option = -1, string label = "") {


  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 300);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 300);
  TH1F *histDenominatorEta = new TH1F ("histDenominatorEta",";Electron Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histNumeratorEta = new TH1F ("histNumeratorEta",";Electron Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histDenominatorPhi = new TH1F ("histDenominatorPhi",";Electron Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histNumeratorPhi = new TH1F ("histNumeratorPhi",";Electron Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histDenominatorRho = new TH1F ("histDenominatorRho",";Electron Rho; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorRho = new TH1F ("histNumeratorRho",";Electron Rho; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpv = new TH1F ("histDenominatorNpv",";Electron Npv; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpv = new TH1F ("histNumeratorNpv",";Electron Npv; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpu = new TH1F ("histDenominatorNpu",";Electron Npu; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpu = new TH1F ("histNumeratorNpu",";Electron Npu; Number of Events", 50, 0 , 100);

  TH2F *histDenominatorPtEta = new TH2F ("histDenominatorPtEta",";Jet p_{T} [GeV/c] ; Jet #eta; Number of Events", 40, 0, 200, 60, -3.0, 3.0);
  TH2F *histNumeratorPtEta = new TH2F ("histNumeratorPtEta",";Jet p_{T} [GeV/c] ; Jet #eta; Number of Events", 40, 0, 200, 60, -3.0, 3.0);

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  JetTree *jetTree = new JetTree;
  jetTree->LoadTree(inputfile.c_str());
  jetTree->InitTree();
  int NCounts = 0;

  cout << "Total Entries: " << jetTree->tree_->GetEntries() << "\n";
  int nentries = jetTree->tree_->GetEntries();
  for(UInt_t ientry=0; ientry < jetTree->tree_->GetEntries(); ientry++) {       	

    jetTree->tree_->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //Selection options
    if (option == 5) {
      if (abs(jetTree->fJetPartonFlavor) != 5) continue;
    }
    if (option == 4) {
      if (abs(jetTree->fJetPartonFlavor) != 4) continue;
    }
    if (option == 0) {
      if (!(abs(jetTree->fJetPartonFlavor) == 21 || (abs(jetTree->fJetPartonFlavor) >= 1 && abs(jetTree->fJetPartonFlavor) <= 3))) continue;
    }

    //Cuts
    if (jetTree->fJetGenPt < 30) continue;
    if (abs(jetTree->fJetGenEta) > 2.4) continue;

    NCounts++;

    // if (NCounts > 1000000 || ientry > 10000000) break;

    //**** PT - ETA ****
    histDenominatorPtEta->Fill(jetTree->fJetGenPt,jetTree->fJetGenEta);
    if(PassSelection(jetTree)) {
      histNumeratorPtEta->Fill(jetTree->fJetGenPt,jetTree->fJetGenEta);
    }


    //**** PT ****
      histDenominatorPt->Fill(jetTree->fJetGenPt);

      //Numerator
      if(PassSelection(jetTree)) {
        histNumeratorPt->Fill(jetTree->fJetGenPt);        
      }


    //**** Eta ****
    if (fabs(jetTree->fJetGenPt) > 30) {
      histDenominatorEta->Fill(jetTree->fJetGenEta);

      //Numerator
      if(PassSelection(jetTree)) {
        histNumeratorEta->Fill(jetTree->fJetGenEta);        
      }

    }

    //**** Phi ****
    if (fabs(jetTree->fJetGenEta) < 2.4) {
      histDenominatorPhi->Fill(jetTree->fJetGenPhi);

      //Numerator
      if(PassSelection(jetTree)) {
        histNumeratorPhi->Fill(jetTree->fJetGenPhi);        
      }

    }

    //**** Rho ****
    if (fabs(jetTree->fJetGenEta) < 2.4) {
      histDenominatorRho->Fill(jetTree->fRho);

      //Numerator
      if(PassSelection(jetTree)) {
        histNumeratorRho->Fill(jetTree->fRho);        
      }

    }
    //**** Npv ****
    if (fabs(jetTree->fJetGenEta) < 2.4) {
      histDenominatorNpv->Fill(jetTree->fNVertices);

      //Numerator
      if(PassSelection(jetTree)) {
        histNumeratorNpv->Fill(jetTree->fNVertices);        
      }

    }

    // //**** Npu ****
    // if (fabs(jetTree->fJetGenEta) < 2.4) {
    //   histDenominatorNpu->Fill(jetTree->);

    //   //Numerator
    //   if(PassSelection(jetTree)) {
    //     histNumeratorNpu->Fill(jetTree->);        
    //   }

    // }


  }


  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency Plots
  //==============================================================================================================

  TGraphAsymmErrors *efficiency_pt = createEfficiencyGraph(histNumeratorPt, histDenominatorPt, "Efficiency_Pt" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_eta = createEfficiencyGraph(histNumeratorEta, histDenominatorEta, "Efficiency_Eta" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_phi = createEfficiencyGraph(histNumeratorPhi, histDenominatorPhi, "Efficiency_Phi" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_rho = createEfficiencyGraph(histNumeratorRho, histDenominatorRho, "Efficiency_Rho" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_npv = createEfficiencyGraph(histNumeratorNpv, histDenominatorNpv, "Efficiency_Npv" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_npu = createEfficiencyGraph(histNumeratorNpu, histDenominatorNpu, "Efficiency_Npu" , vector<double>() ,  -99, -99, 0, 1);  
  TH2F *efficiency_pteta = createEfficiencyHist2D(histNumeratorPtEta, histDenominatorPtEta, "Efficiency_PtEta" , vector<double>() ,vector<double>());  


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;

  cv = new TCanvas("cv","cv",800,600);
  efficiency_pt->Draw("AP");
  efficiency_pt->SetTitle("");
  efficiency_pt->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Pt.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_eta->Draw("AP");
  efficiency_eta->SetTitle("");
  efficiency_eta->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Eta.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_phi->Draw("AP");
  efficiency_phi->SetTitle("");
  efficiency_phi->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Phi.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_rho->Draw("AP");
  efficiency_rho->SetTitle("");
  efficiency_rho->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Rho.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npv->Draw("AP");
  efficiency_npv->SetTitle("");
  efficiency_npv->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Npv.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npu->Draw("AP");
  efficiency_npu->SetTitle("");
  efficiency_npu->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Npu.gif").c_str());


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("Efficiency"+Label+".root").c_str(), "UPDATE");
  file->cd();
  file->WriteTObject(efficiency_pt, "Efficiency_Pt", "WriteDelete");
  file->WriteTObject(efficiency_eta, "Efficiency_Eta", "WriteDelete");
  file->WriteTObject(efficiency_phi, "Efficiency_Phi", "WriteDelete");
  file->WriteTObject(efficiency_rho, "Efficiency_Rho", "WriteDelete");
  file->WriteTObject(efficiency_npv, "Efficiency_NPV", "WriteDelete");
  file->WriteTObject(efficiency_npu, "Efficiency_NPU", "WriteDelete");
  file->WriteTObject(efficiency_pteta, "Efficiency_PtEta", "WriteDelete");

  file->Close();
  delete file;       

}


void MakeFastsimToFullSimCorrectionFactors() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileFullsimMedium = new TFile("Efficiency_BJets_25ns_Fullsim.root","READ");
  TFile *fileFastsimMedium = new TFile("Efficiency_BJets_25ns_Fastsim.root","READ");
  TH2F* histFullsimMedium = (TH2F*)fileFullsimMedium->Get("Efficiency_PtEta");
  TH2F* histFastsimMedium = (TH2F*)fileFastsimMedium->Get("Efficiency_PtEta");

  TH2F* histSFMedium = (TH2F*)histFullsimMedium->Clone("BTagMedium_FastsimScaleFactor");
  histSFMedium->GetXaxis()->SetTitle("Jet p_{T} [GeV/c]");
  histSFMedium->GetYaxis()->SetTitle("Jet #eta");

  //Medium WP
  for (int b=1; b<histSFMedium->GetXaxis()->GetNbins()+1 ; ++b) {
    for (int c=1; c<histSFMedium->GetYaxis()->GetNbins()+1 ; ++c) {
      double sf = histFullsimMedium->GetBinContent(b,c) / histFastsimMedium->GetBinContent(b,c);
      double sferr = 0;
      if ( histFullsimMedium->GetBinContent(b,c) > 0 && histFastsimMedium->GetBinContent(b,c) > 0) {
	sferr = sf*sqrt( pow(histFullsimMedium->GetBinError(b,c)/histFullsimMedium->GetBinContent(b,c),2) +
				pow(histFastsimMedium->GetBinError(b,c)/histFastsimMedium->GetBinContent(b,c),2) );
      }
      histSFMedium->SetBinContent(b,c,sf);
      histSFMedium->SetBinError(b,c,sferr);
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open("BTagEffFastsimToFullsimCorrectionFactors.root", "UPDATE");
  file->cd();
  file->WriteTObject(histSFMedium, "BTagMedium_FastsimScaleFactor", "WriteDelete");  
  file->WriteTObject(histFullsimMedium, "BTagEff_Medium_Fullsim", "WriteDelete");
  file->WriteTObject(histFastsimMedium, "BTagEff_Medium_Fastsim", "WriteDelete");
  file->Close();
  delete file;      

}





void MakeBTaggingEfficiencyPlots( int Option = 0) {
 
  if (Option==1) {
    ProduceBTaggingEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/JetNtuple/V1p20/JetNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fullsim.root", 5 , "BJets_25ns_Fullsim");
    // ProduceBTaggingEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/JetNtuple/JetNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Spring15_25ns.root", 4 , "CharmJets_25ns");
    // ProduceBTaggingEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/JetNtuple/JetNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Spring15_25ns.root", 0 , "LightJets_25ns");    
  }
  
  if (Option==2) {
    ProduceBTaggingEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/JetNtuple/V1p20/JetNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fastsim.root", 5 , "BJets_25ns_Fastsim");
  }
  
  
  //plotBTaggingEfficiency();
  if (Option == 100) {
    MakeFastsimToFullSimCorrectionFactors();
  }

}
