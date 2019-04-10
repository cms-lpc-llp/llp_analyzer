//================================================================================================
//
// Simple Example
//
//root -l RazorAnalyzer/macros/ObjectStudies/MakeMuonEfficiencyPlots.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/MuonNtuple/MuonNtuple_PromptGenLevel_TTJets_25ns.root",-1,"Muon")'
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
#include "RazorAnalyzer/include/MuonTree.h"

#endif

double GetEffArea( double eta ) {

  double effArea = 0.0; 
  //Effective areas using mean region

  if (fabs(eta) < 0.8) {
      effArea = 0.0735;
    } else if (fabs(eta) < 1.3) {
      effArea = 0.0619;	
    } else if (fabs(eta) < 2.0) {
      effArea = 0.0465;	
    } else if (fabs(eta) < 2.2) {
      effArea = 0.0433;	
    } else {
      effArea = 0.0577;	
    }  

  return effArea;

}

bool PassSelection( MuonTree* MuTree, int wp ) {

  bool pass = false;
  double dr = fmax(0.05,fmin(0.2, 10/MuTree->fMuPt));

  //improved isolation

  //**********************************
  //Tight Selection
  //**********************************
  if (wp == 3) {
    if (MuTree->fMuPt > 0 
	&& MuTree->fPassTightSelection
	//&& MuTree->fMuIsTight && fabs(MuTree->fMuIP3dSig)<4 && fabs(MuTree->fMuD0) < 0.2 && MuTree->fMuPFIso04 < 0.12
	) {
      pass = true;
    }
  } 

  //**********************************
  //Loose Selection
  //**********************************
  if (wp == 2) {
    if (MuTree->fMuPt > 0 
	//&& MuTree->fPassLooseSelection
	&& MuTree->fMuIsLoose && fabs(MuTree->fMuIP3dSig)<4
	//&& (MuTree->fMiniIsoCharged + fmax(0.0, MuTree->fMiniIso - MuTree->fRho*GetEffArea(MuTree->fMuEta)*pow(dr/0.3,2))) / MuTree->fMuPt < 0.2
	&& (MuTree->fMiniIsoCharged + fmax(0.0, MuTree->fMiniIso - MuTree->fRhoNeutralCentral*GetEffArea(MuTree->fMuEta)*pow(dr/0.3,2))) / MuTree->fMuPt < 0.2
	) {
      pass = true;
    }
  }
  
  //**********************************
  //Veto Selection
  //**********************************
  if (wp == 1 || wp == 10) {
    if (MuTree->fMuPt > 20) {
      if (MuTree->fMuPt > 0 
	  && MuTree->fPassVetoSelection
  	  //&& MuTree->fMuIsLoose && fabs(MuTree->fMuIP3dSig)<4 
  	  //&& MuTree->fMiniIso < 0.2
  	  ) {
  	pass = true;
      }
    } else {
      if (MuTree->fMuPt > 0 
	  && MuTree->fPassVetoSelection
	  //&& MuTree->fMuIsLoose && fabs(MuTree->fMuIP3dSig)<4 && MuTree->fMuPFIso04*MuTree->fMuPt < 10
	  ) {
  	pass = true;
      }
    }
  }

  //**********************************
  //Relative Isolation Veto Selection
  //**********************************
  if (wp == 11) {
    if (MuTree->fMuPt > 20) {
      if (MuTree->fMuPt > 0 
	  && MuTree->fMuIsLoose && fabs(MuTree->fMuIP3dSig)<4 
	  && MuTree->fMuPFIso04 < 0.4
	  ) {
	pass = true;
      }
    } else {
      if (MuTree->fMuPt > 0 && MuTree->fMuIsLoose && fabs(MuTree->fMuIP3dSig)<4 && MuTree->fMuPFIso04 < 0.4) {
	pass = true;
      }
    }
  }

  if (wp == 100) {
    // cout << "wp = 100\n";
    // cout << MuTree->fMuTriggerBit << "\n";
    // cout << bool((MuTree->fMuTriggerBit & MuonTree::kMuTrigger_IsoMu20) == MuonTree::kMuTrigger_IsoMu20) << " " 
    // 	 << bool ((MuTree->fMuTriggerBit & MuonTree::kMuTrigger_IsoTkMu20) == MuonTree::kMuTrigger_IsoTkMu20) << " "
    // 	 << bool( (MuTree->fMuTriggerBit & MuonTree::kMuTrigger_IsoMu20) == MuonTree::kMuTrigger_IsoMu20
    // 		  || (MuTree->fMuTriggerBit & MuonTree::kMuTrigger_IsoTkMu20) == MuonTree::kMuTrigger_IsoTkMu20) << "\n";

    if ( (MuTree->fMuTriggerBit & MuonTree::kMuTrigger_IsoMu20) == MuonTree::kMuTrigger_IsoMu20
	 || (MuTree->fMuTriggerBit & MuonTree::kMuTrigger_IsoTkMu20) == MuonTree::kMuTrigger_IsoTkMu20
	 ) pass = true;
  }

  if (wp == 101) {
    if ( (MuTree->fMuTriggerBit & MuonTree::kMuTrigger_IsoMu27) == MuonTree::kMuTrigger_IsoMu27
	 || (MuTree->fMuTriggerBit & MuonTree::kMuTrigger_IsoTkMu27) == MuonTree::kMuTrigger_IsoTkMu27
	 ) pass = true;
  }

  if (wp == 102) {
    if ( (MuTree->fMuTriggerBit & MuonTree::kMuTrigger_Mu50) == MuonTree::kMuTrigger_Mu50
	 ) pass = true;
  }

  if (wp == 103) {
    if ( (MuTree->fMuTriggerBit & MuonTree::kMuTrigger_Mu50_eta2p1) == MuonTree::kMuTrigger_Mu50_eta2p1
	 ) pass = true;
  }

  if (wp == 110) {
    if ( (MuTree->fMuTriggerBit & MuonTree::kMuTrigger_IsoMu27) == MuonTree::kMuTrigger_IsoMu27
	 || (MuTree->fMuTriggerBit & MuonTree::kMuTrigger_IsoTkMu27) == MuonTree::kMuTrigger_IsoTkMu27
	 || (MuTree->fMuTriggerBit & MuonTree::kMuTrigger_Mu50) == MuonTree::kMuTrigger_Mu50
	 ) pass = true;
  }


  // reco only
  if (wp ==0) {
     if (MuTree->fMuPt > 0) pass = true;
  }

  return pass;

}


void plotMuonEfficiency() {

  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileRelIso = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/MyNotes/notes/AN-14-276/trunk/data/Efficiency_Muon_NumeratorRelIso0p4_DenominatorLooseIDAndIPCut.root","READ");
  TFile *fileImprovedIso = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/MyNotes/notes/AN-14-276/trunk/data/Efficiency_Muon_NumeratorImprovedIso_DenominatorLooseIDAndIPCut.root","READ");
  TFile *fileVeto = new TFile("Efficiency_PromptMuon_TTJets_50ns_Veto.root","READ");
  TFile *fileLoose = new TFile("Efficiency_PromptMuon_TTJets_50ns_Loose.root","READ");
  TFile *fileTight = new TFile("Efficiency_PromptMuon_TTJets_50ns_Tight.root","READ");
  TFile *fileFakesVeto = new TFile("Efficiency_FakeMuon_TTJets_50ns_Veto.root","READ");
  TFile *fileFakesLoose = new TFile("Efficiency_FakeMuon_TTJets_50ns_Loose.root","READ");
  TFile *fileFakesTight = new TFile("Efficiency_FakeMuon_TTJets_50ns_Tight.root","READ");

  TGraphAsymmErrors* effPtRelIso = (TGraphAsymmErrors*)fileRelIso->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtImprovedIso = (TGraphAsymmErrors*)fileImprovedIso->Get("Efficiency_Pt");

  TGraphAsymmErrors* effPtVeto = (TGraphAsymmErrors*)fileVeto->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtLoose = (TGraphAsymmErrors*)fileLoose->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtTight = (TGraphAsymmErrors*)fileTight->Get("Efficiency_Pt");
  TGraphAsymmErrors* effEtaVeto = (TGraphAsymmErrors*)fileVeto->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaLoose = (TGraphAsymmErrors*)fileLoose->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaTight = (TGraphAsymmErrors*)fileTight->Get("Efficiency_Eta");
  TGraphAsymmErrors* effNpvVeto = (TGraphAsymmErrors*)fileVeto->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNpvLoose = (TGraphAsymmErrors*)fileLoose->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNpvTight = (TGraphAsymmErrors*)fileTight->Get("Efficiency_NPV");

  TGraphAsymmErrors* effFakePtVeto = (TGraphAsymmErrors*)fileFakesVeto->Get("Efficiency_Pt");
  TGraphAsymmErrors* effFakePtLoose = (TGraphAsymmErrors*)fileFakesLoose->Get("Efficiency_Pt");
  TGraphAsymmErrors* effFakePtTight = (TGraphAsymmErrors*)fileFakesTight->Get("Efficiency_Pt");
  TGraphAsymmErrors* effFakeEtaVeto = (TGraphAsymmErrors*)fileFakesVeto->Get("Efficiency_Eta");
  TGraphAsymmErrors* effFakeEtaLoose = (TGraphAsymmErrors*)fileFakesLoose->Get("Efficiency_Eta");
  TGraphAsymmErrors* effFakeEtaTight = (TGraphAsymmErrors*)fileFakesTight->Get("Efficiency_Eta");
  TGraphAsymmErrors* effFakeNpvVeto = (TGraphAsymmErrors*)fileFakesVeto->Get("Efficiency_NPV");
  TGraphAsymmErrors* effFakeNpvLoose = (TGraphAsymmErrors*)fileFakesLoose->Get("Efficiency_NPV");
  TGraphAsymmErrors* effFakeNpvTight = (TGraphAsymmErrors*)fileFakesTight->Get("Efficiency_NPV");



  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtRelIso, "Relative Iso < 0.4", "LP");
  legend->AddEntry(effPtImprovedIso, "Improved Veto Isolation", "LP");

  effPtRelIso->SetLineWidth(3);
  effPtRelIso->SetLineColor(kBlack);
  effPtRelIso->GetXaxis()->SetTitle("Muon p_{T} [GeV/c]");
  effPtRelIso->GetYaxis()->SetTitle("Isolation Efficiency");
  effPtRelIso->GetYaxis()->SetTitleOffset(1.2);

  effPtRelIso->SetLineWidth(3);
  effPtRelIso->SetLineColor(kRed);
  effPtImprovedIso->SetLineWidth(3);
  effPtImprovedIso->SetLineColor(kBlue);

  effPtRelIso->Draw("AP");
  effPtImprovedIso->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("MuonIsolationEfficiencyVsPt.gif");
  cv->SaveAs("MuonIsolationEfficiencyVsPt.pdf");



  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtVeto, "Veto", "LP");
  legend->AddEntry(effPtLoose, "Loose", "LP");
  legend->AddEntry(effPtTight, "Tight", "LP");

  effPtTight->SetLineWidth(3);
  effPtTight->SetLineColor(kRed);
  effPtTight->GetXaxis()->SetTitle("Muon p_{T} [GeV/c]");
  effPtTight->GetYaxis()->SetTitle("Selection Efficiency");
  effPtTight->GetYaxis()->SetTitleOffset(1.2);

  effPtLoose->SetLineWidth(3);
  effPtLoose->SetLineColor(kBlue);
  effPtVeto->SetLineWidth(3);
  effPtVeto->SetLineColor(kBlack);

  effPtTight->Draw("AP");
  effPtVeto->Draw("Psame");
  effPtLoose->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("MuonSelectionEfficiencyVsPt.gif");
  cv->SaveAs("MuonSelectionEfficiencyVsPt.pdf");



  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effEtaVeto, "Veto", "LP");
  legend->AddEntry(effEtaLoose, "Loose", "LP");
  legend->AddEntry(effEtaTight, "Tight", "LP");

  effEtaTight->SetLineWidth(3);
  effEtaTight->SetLineColor(kRed);
  effEtaTight->GetXaxis()->SetTitle("Muon #eta");
  effEtaTight->GetYaxis()->SetTitle("Selection Efficiency");
  effEtaTight->GetYaxis()->SetTitleOffset(1.2);

  effEtaLoose->SetLineWidth(3);
  effEtaLoose->SetLineColor(kBlue);
  effEtaVeto->SetLineWidth(3);
  effEtaVeto->SetLineColor(kBlack);

  effEtaTight->Draw("AP");
  effEtaVeto->Draw("Psame");
  effEtaLoose->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("MuonSelectionEfficiencyVsEta.gif");
  cv->SaveAs("MuonSelectionEfficiencyVsEta.pdf");




  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effNpvVeto, "Veto", "LP");
  legend->AddEntry(effNpvLoose, "Loose", "LP");
  legend->AddEntry(effNpvTight, "Tight", "LP");

  effNpvTight->SetLineWidth(3);
  effNpvTight->SetLineColor(kRed);
  effNpvTight->GetXaxis()->SetTitle("Number of Reconstructed Vertices");
  effNpvTight->GetYaxis()->SetTitle("Selection Efficiency");
  effNpvTight->GetYaxis()->SetTitleOffset(1.2);
  effNpvTight->GetXaxis()->SetRangeUser(5,35);

  effNpvLoose->SetLineWidth(3);
  effNpvLoose->SetLineColor(kBlue);
  effNpvVeto->SetLineWidth(3);
  effNpvVeto->SetLineColor(kBlack);

  effNpvTight->Draw("AP");
  effNpvVeto->Draw("Psame");
  effNpvLoose->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("MuonSelectionEfficiencyVsNpv.gif");
  cv->SaveAs("MuonSelectionEfficiencyVsNpv.pdf");




  //***************************************************************
  //Fake Muons : Efficiency Vs Pt
  //***************************************************************
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.70,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effFakePtVeto, "Veto", "LP");
  legend->AddEntry(effFakePtLoose, "Loose", "LP");
  legend->AddEntry(effFakePtTight, "Tight", "LP");

  effFakePtTight->SetLineWidth(3);
  effFakePtTight->SetLineColor(kRed);
  effFakePtTight->GetXaxis()->SetTitle("Muon p_{T} [GeV/c]");
  effFakePtTight->GetYaxis()->SetTitle("Selection Efficiency");
  effFakePtTight->GetYaxis()->SetTitleOffset(1.35);

  effFakePtLoose->SetLineWidth(3);
  effFakePtLoose->SetLineColor(kBlue);
  effFakePtVeto->SetLineWidth(3);
  effFakePtVeto->SetLineColor(kBlack);

  effFakePtTight->Draw("AP");
  effFakePtLoose->Draw("Psame");
  effFakePtVeto->Draw("Psame");
  
  effFakePtTight->GetYaxis()->SetRangeUser(0,0.05);

  legend->Draw();  
  //cv->SetLogy();
  cv->SaveAs("MuonSelectionFakeRateVsPt.gif");
  cv->SaveAs("MuonSelectionFakeRateVsPt.pdf");


  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.70,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effFakeEtaVeto, "Veto", "LP");
  legend->AddEntry(effFakeEtaLoose, "Loose", "LP");
  legend->AddEntry(effFakeEtaTight, "Tight", "LP");

  effFakeEtaTight->SetLineWidth(3);
  effFakeEtaTight->SetLineColor(kRed);
  effFakeEtaTight->GetXaxis()->SetTitle("Muon #eta");
  effFakeEtaTight->GetYaxis()->SetTitle("Selection Efficiency");
  effFakeEtaTight->GetYaxis()->SetTitleOffset(1.35);

  effFakeEtaLoose->SetLineWidth(3);
  effFakeEtaLoose->SetLineColor(kBlue);
  effFakeEtaVeto->SetLineWidth(3);
  effFakeEtaVeto->SetLineColor(kBlack);

  effFakeEtaTight->Draw("AP");
  effFakeEtaLoose->Draw("Psame");
  effFakeEtaVeto->Draw("Psame");
  
  effFakeEtaTight->GetYaxis()->SetRangeUser(0,0.04);

  legend->Draw();  
  //cv->SetLogy();
  cv->SaveAs("MuonSelectionFakeRateVsEta.gif");
  cv->SaveAs("MuonSelectionFakeRateVsEta.pdf");



  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.70,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effFakeNpvVeto, "Veto", "LP");
  legend->AddEntry(effFakeNpvLoose, "Loose", "LP");
  legend->AddEntry(effFakeNpvTight, "Tight", "LP");

  effFakeNpvTight->SetLineWidth(3);
  effFakeNpvTight->SetLineColor(kRed);
  effFakeNpvTight->GetXaxis()->SetTitle("Number of Reconstructed Vertices");
  effFakeNpvTight->GetYaxis()->SetTitle("Selection Efficiency");
  effFakeNpvTight->GetYaxis()->SetTitleOffset(1.35);
  effFakeNpvTight->GetXaxis()->SetRangeUser(5,35);
  effFakeNpvTight->GetYaxis()->SetRangeUser(0,0.04);

  effFakeNpvLoose->SetLineWidth(3);
  effFakeNpvLoose->SetLineColor(kBlue);
  effFakeNpvVeto->SetLineWidth(3);
  effFakeNpvVeto->SetLineColor(kBlack);

  effFakeNpvTight->Draw("AP");
  effFakeNpvLoose->Draw("Psame");
  effFakeNpvVeto->Draw("Psame"); 

  legend->Draw();  
  //cv->SetLogy();
  cv->SaveAs("MuonSelectionFakeRateVsNpv.gif");
  cv->SaveAs("MuonSelectionFakeRateVsNpv.pdf");


}


void plotMuonTriggerEfficiency() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileTriggerIsoMu20 = new TFile("Efficiency_PromptMuon_TTJets_50ns_MuTriggerIsoMu20.root","READ");
  TFile *fileTriggerIsoMu27 = new TFile("Efficiency_PromptMuon_TTJets_50ns_MuTriggerIsoMu27.root","READ");
  TFile *fileTriggerMu50 = new TFile("Efficiency_PromptMuon_TTJets_50ns_MuTriggerMu50.root","READ");
  TFile *fileTriggerMuCombined = new TFile("Efficiency_PromptMuon_TTJets_50ns_MuTriggerCombined.root","READ");
  
  TGraphAsymmErrors* effPtTriggerIsoMu20 = (TGraphAsymmErrors*)fileTriggerIsoMu20->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtTriggerIsoMu27 = (TGraphAsymmErrors*)fileTriggerIsoMu27->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtTriggerMu50 = (TGraphAsymmErrors*)fileTriggerMu50->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtTriggerMuCombined = (TGraphAsymmErrors*)fileTriggerMuCombined->Get("Efficiency_Pt");

  TGraphAsymmErrors* effEtaTriggerIsoMu20 = (TGraphAsymmErrors*)fileTriggerIsoMu20->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaTriggerIsoMu27 = (TGraphAsymmErrors*)fileTriggerIsoMu27->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaTriggerMu50 = (TGraphAsymmErrors*)fileTriggerMu50->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaTriggerMuCombined = (TGraphAsymmErrors*)fileTriggerMuCombined->Get("Efficiency_Eta");

  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtTriggerIsoMu20, "IsoMu20", "LP");
  legend->AddEntry(effPtTriggerIsoMu27, "IsoMu27", "LP");
  legend->AddEntry(effPtTriggerMu50, "Mu50", "LP");
  legend->AddEntry(effPtTriggerMuCombined, "IsoMu27 OR Mu50", "LP");

  effPtTriggerIsoMu20->SetLineWidth(3);
  effPtTriggerIsoMu20->SetLineColor(kBlack);
  effPtTriggerIsoMu20->GetXaxis()->SetTitle("Muon p_{T} [GeV/c]");
  effPtTriggerIsoMu20->GetYaxis()->SetTitle("Trigger Efficiency");
  effPtTriggerIsoMu20->GetYaxis()->SetTitleOffset(1.2);

  effPtTriggerIsoMu27->SetLineWidth(3);
  effPtTriggerIsoMu27->SetLineColor(kBlue);
  effPtTriggerMu50->SetLineWidth(3);
  effPtTriggerMu50->SetLineColor(kRed);
  effPtTriggerMuCombined->SetLineWidth(3);
  effPtTriggerMuCombined->SetLineColor(kGreen+2);

  effPtTriggerIsoMu20->Draw("AP");
  effPtTriggerIsoMu27->Draw("Psame");
  effPtTriggerMu50->Draw("Psame");
  effPtTriggerMuCombined->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("SingleMuTriggerEfficiencyVsPt.gif");
  cv->SaveAs("SingleMuTriggerEfficiencyVsPt.pdf");



  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effEtaTriggerIsoMu20, "IsoMu20", "LP");
  legend->AddEntry(effEtaTriggerIsoMu27, "IsoMu27", "LP");
  legend->AddEntry(effEtaTriggerMu50, "Mu50", "LP");
  legend->AddEntry(effEtaTriggerMuCombined, "IsoMu27 OR Mu50", "LP");

  effEtaTriggerIsoMu20->SetLineWidth(3);
  effEtaTriggerIsoMu20->SetLineColor(kBlack);
  effEtaTriggerIsoMu20->GetXaxis()->SetTitle("Muon #eta");
  effEtaTriggerIsoMu20->GetYaxis()->SetTitle("Trigger Efficiency");
  effEtaTriggerIsoMu20->GetYaxis()->SetTitleOffset(1.2);

  effEtaTriggerIsoMu27->SetLineWidth(3);
  effEtaTriggerIsoMu27->SetLineColor(kBlue);
  effEtaTriggerMu50->SetLineWidth(3);
  effEtaTriggerMu50->SetLineColor(kRed);
  effEtaTriggerMuCombined->SetLineWidth(3);
  effEtaTriggerMuCombined->SetLineColor(kGreen+2);

  effEtaTriggerIsoMu20->Draw("AP");
  effEtaTriggerIsoMu27->Draw("Psame");
  effEtaTriggerMu50->Draw("Psame");
  effEtaTriggerMuCombined->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("SingleMuTriggerEfficiencyVsEta.gif");
  cv->SaveAs("SingleMuTriggerEfficiencyVsEta.pdf");

}


void plotMuonMiniIsoEfficiency() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *filePromptEACorr = new TFile("Efficiency_PromptMuon_TTJets_25_MiniIsolationEACorrRho.root","READ");
  TFile *filePromptEACorrNeutralCentral = new TFile("Efficiency_PromptMuon_TTJets_25_MiniIsolationEACorrRhoNeutralCentral.root","READ");
  TFile *fileFakeEACorr = new TFile("Efficiency_FakeMuon_TTJets_25_MiniIsolationEACorrRho.root","READ");
  TFile *fileFakeEACorrNeutralCentral = new TFile("Efficiency_FakeMuon_TTJets_25_MiniIsolationEACorrRhoNeutralCentral.root","READ");

  TGraphAsymmErrors* effNPVPromptEACorr = (TGraphAsymmErrors*)filePromptEACorr->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNPVPromptEACorrNeutralCentral = (TGraphAsymmErrors*)filePromptEACorrNeutralCentral->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNPVFakeEACorr = (TGraphAsymmErrors*)fileFakeEACorr->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNPVFakeEACorrNeutralCentral = (TGraphAsymmErrors*)fileFakeEACorrNeutralCentral->Get("Efficiency_NPV");

 
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.60,0.90,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effNPVPromptEACorr, "EA Corr (Std Rho)", "LP");
  legend->AddEntry(effNPVPromptEACorrNeutralCentral, "EA Corr (Rho CentralNeutral)", "LP");

  effNPVPromptEACorr->SetLineWidth(3);
  effNPVPromptEACorr->SetLineColor(kBlack);
  effNPVPromptEACorr->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
  effNPVPromptEACorr->GetYaxis()->SetTitle("MiniIsolation Cut Efficiency");
  effNPVPromptEACorr->GetYaxis()->SetTitleOffset(1.4);

  effNPVPromptEACorr->SetLineWidth(3);
  effNPVPromptEACorr->SetLineColor(kBlue);
  effNPVPromptEACorrNeutralCentral->SetLineWidth(3);
  effNPVPromptEACorrNeutralCentral->SetLineColor(kRed);

  effNPVPromptEACorr->Draw("AP");
  effNPVPromptEACorrNeutralCentral->Draw("Psame");

  effNPVPromptEACorr->GetYaxis()->SetRangeUser(0.90,1.0);
  effNPVPromptEACorr->GetXaxis()->SetRangeUser(2,30);
  
  legend->Draw();  
  cv->SaveAs("MiniIsoEfficiencyVsNPV_PromptMuons.gif");


  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.65,0.90,0.85);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effNPVFakeEACorr, "EA Corr (Std Rho)", "LP");
  legend->AddEntry(effNPVFakeEACorrNeutralCentral, "EA Corr (Rho CentralNeutral)", "LP");

  effNPVFakeEACorr->SetLineWidth(3);
  effNPVFakeEACorr->SetLineColor(kBlack);
  effNPVFakeEACorr->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
  effNPVFakeEACorr->GetYaxis()->SetTitle("MiniIsolation Cut Efficiency");
  effNPVFakeEACorr->GetYaxis()->SetTitleOffset(1.4);

  effNPVFakeEACorr->SetLineWidth(3);
  effNPVFakeEACorr->SetLineColor(kBlue);
  effNPVFakeEACorrNeutralCentral->SetLineWidth(3);
  effNPVFakeEACorrNeutralCentral->SetLineColor(kRed);

  effNPVFakeEACorr->Draw("AP");
  effNPVFakeEACorrNeutralCentral->Draw("Psame");

  effNPVFakeEACorr->GetYaxis()->SetRangeUser(0.0,0.20);
  effNPVFakeEACorr->GetXaxis()->SetRangeUser(2,30);
  
  legend->Draw();  
  cv->SaveAs("MiniIsoEfficiencyVsNPV_FakeMuons.gif");




}



void MakeFastsimToFullSimCorrectionFactors() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileFullsimTight = new TFile("Efficiency_PromptMuon_TTJets_25ns_Tight_Fullsim.root","READ");
  TFile *fileFullsimVeto = new TFile("Efficiency_PromptMuon_TTJets_25ns_Veto_Fullsim.root","READ");
  TFile *fileFastsimTight = new TFile("Efficiency_PromptMuon_TTJets_25ns_Tight_Fastsim.root","READ");
  TFile *fileFastsimVeto = new TFile("Efficiency_PromptMuon_TTJets_25ns_Veto_Fastsim.root","READ");

  TH2F* histFullsimTight = (TH2F*)fileFullsimTight->Get("Efficiency_PtEta");
  TH2F* histFullsimVeto = (TH2F*)fileFullsimVeto->Get("Efficiency_PtEta");
  TH2F* histFastsimTight = (TH2F*)fileFastsimTight->Get("Efficiency_PtEta");
  TH2F* histFastsimVeto = (TH2F*)fileFastsimVeto->Get("Efficiency_PtEta");

  TH2F* histSFTight = (TH2F*)histFullsimTight->Clone("MuonTight_FastsimScaleFactor");
  TH2F* histSFVeto = (TH2F*)histFullsimVeto->Clone("MuonVeto_FastsimScaleFactor");
  histSFTight->GetXaxis()->SetTitle("Muon p_{T} [GeV/c]");
  histSFTight->GetYaxis()->SetTitle("Muon #eta");
  histSFVeto->GetXaxis()->SetTitle("Muon p_{T} [GeV/c]");
  histSFVeto->GetYaxis()->SetTitle("Muon #eta");

  //Tight WP
  for (int b=1; b<histSFTight->GetXaxis()->GetNbins()+1 ; ++b) {
    for (int c=1; c<histSFTight->GetYaxis()->GetNbins()+1 ; ++c) {
      double sf = histFullsimTight->GetBinContent(b,c) / histFastsimTight->GetBinContent(b,c);
      double sferr = 0;
      if ( histFullsimTight->GetBinContent(b,c) > 0 && histFastsimTight->GetBinContent(b,c) > 0) {
	sferr = sf*sqrt( pow(histFullsimTight->GetBinError(b,c)/histFullsimTight->GetBinContent(b,c),2) +
				pow(histFastsimTight->GetBinError(b,c)/histFastsimTight->GetBinContent(b,c),2) );
      }
      histSFTight->SetBinContent(b,c,sf);
      histSFTight->SetBinError(b,c,sferr);
    }
  }

  //Veto WP
  for (int b=1; b<histSFVeto->GetXaxis()->GetNbins()+1 ; ++b) {
    for (int c=1; c<histSFVeto->GetYaxis()->GetNbins()+1 ; ++c) {
      double sf = histFullsimVeto->GetBinContent(b,c) / histFastsimVeto->GetBinContent(b,c);
      double sferr = 0;
      if ( histFullsimVeto->GetBinContent(b,c) > 0 && histFastsimVeto->GetBinContent(b,c) > 0) {
	sferr = sf*sqrt( pow(histFullsimVeto->GetBinError(b,c)/histFullsimVeto->GetBinContent(b,c),2) +
			 pow(histFastsimVeto->GetBinError(b,c)/histFastsimVeto->GetBinContent(b,c),2) );
      }
      histSFVeto->SetBinContent(b,c,sf);
      histSFVeto->SetBinError(b,c,sferr);
    }
  }

   //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open("MuonEffFastsimToFullsimCorrectionFactors.root", "UPDATE");
  file->cd();
  file->WriteTObject(histSFTight, "MuonTight_FastsimScaleFactor", "WriteDelete");  
  file->WriteTObject(histSFVeto, "MuonVeto_FastsimScaleFactor", "WriteDelete");
  file->WriteTObject(histFullsimTight, "MuonEff_Tight_Fullsim", "WriteDelete");
  file->WriteTObject(histFullsimVeto, "MuonEff_Veto_Fullsim", "WriteDelete");
  file->WriteTObject(histFastsimTight, "MuonEff_Tight_Fastsim", "WriteDelete");
  file->WriteTObject(histFastsimVeto, "MuonEff_Veto_Fastsim", "WriteDelete");
  file->Close();
  delete file;      

}




//=== MAIN MACRO ================================================================================================= 

void ProduceMuonEfficiencyPlots(const string inputfile, int wp = 0,  int option = -1, string label = "") {

  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt",";Muon p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt",";Muon p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorEta = new TH1F ("histDenominatorEta",";Muon Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histNumeratorEta = new TH1F ("histNumeratorEta",";Muon Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histDenominatorPhi = new TH1F ("histDenominatorPhi",";Muon Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histNumeratorPhi = new TH1F ("histNumeratorPhi",";Muon Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histDenominatorRho = new TH1F ("histDenominatorRho",";Muon Rho; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorRho = new TH1F ("histNumeratorRho",";Muon Rho; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpv = new TH1F ("histDenominatorNpv",";Muon Npv; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpv = new TH1F ("histNumeratorNpv",";Muon Npv; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpu = new TH1F ("histDenominatorNpu",";Muon Npu; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpu = new TH1F ("histNumeratorNpu",";Muon Npu; Number of Events", 50, 0 , 100);

  TH2F *histDenominatorPtEta = new TH2F ("histDenominatorPtEta",";Muon p_{T} [GeV/c] ; Muon #eta; Number of Events", 34, 0 , 170, 60, -3.0, 3.0);
  TH2F *histNumeratorPtEta = new TH2F ("histNumeratorPtEta",";Muon p_{T} [GeV/c] ; Muon #eta; Number of Events", 34, 0 , 170, 60, -3.0, 3.0);

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  MuonTree *MuTree = new MuonTree;
  MuTree->LoadTree(inputfile.c_str());
  MuTree->InitTree(MuonTree::kMuTreeLight);

  cout << "Total Entries: " << MuTree->tree_->GetEntries() << "\n";
  int nentries = MuTree->tree_->GetEntries();
  for(UInt_t ientry=0; ientry < MuTree->tree_->GetEntries(); ientry++) {       	
    MuTree->tree_->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //Cuts
    if (MuTree->fMuGenPt < 5) continue;
    if (abs(MuTree->fMuGenEta) > 2.4) continue;

    if (!(MuTree->fMuPt > 5)) continue;

    //For Iso only efficiency require ID cuts   
    if (wp == 10 || wp == 11) {
      if (!(MuTree->fMuPt > 0 && MuTree->fMuIsLoose  && fabs(MuTree->fMuIP3dSig)<4)) continue;
    } 

    //For Trigger efficiency require tight selection
    if (wp >= 100 ){
      if (!(MuTree->fPassTightSelection)) continue;
    }

    // //Require Loose muons
    // if (!(MuTree->fMuIsLoose && fabs(MuTree->fMuIP3dSig)<4)) continue;


    if (option==0) {
      //**** PT - ETA ****
      histDenominatorPtEta->Fill(MuTree->fMuGenPt,MuTree->fMuGenEta);
      if(PassSelection(MuTree,wp)) {
	histNumeratorPtEta->Fill(MuTree->fMuGenPt,MuTree->fMuGenEta);
      }


      //**** PT ****
      histDenominatorPt->Fill(MuTree->fMuGenPt);

      //Numerator
      if(PassSelection(MuTree,wp)) {
        histNumeratorPt->Fill(MuTree->fMuGenPt);        
      }


      //**** Eta ****
      if (fabs(MuTree->fMuGenPt) > 30) {
	histDenominatorEta->Fill(MuTree->fMuGenEta);

	//Numerator
	if(PassSelection(MuTree,wp)) {
	  histNumeratorEta->Fill(MuTree->fMuGenEta);        
	}

      }

      //**** Phi ****
      if (fabs(MuTree->fMuGenEta) < 2.4) {
	histDenominatorPhi->Fill(MuTree->fMuGenPhi);

	//Numerator
	if(PassSelection(MuTree,wp)) {
	  histNumeratorPhi->Fill(MuTree->fMuGenPhi);        
	}

      }

      //**** Rho ****
      if (fabs(MuTree->fMuGenEta) < 2.4) {
	histDenominatorRho->Fill(MuTree->fRho);

	//Numerator
	if(PassSelection(MuTree,wp)) {
	  histNumeratorRho->Fill(MuTree->fRho);        
	}

      }
      //**** Npv ****
      if (fabs(MuTree->fMuGenEta) < 2.4) {
	histDenominatorNpv->Fill(MuTree->fNVertices);

	//Numerator
	if(PassSelection(MuTree,wp)) {
	  histNumeratorNpv->Fill(MuTree->fNVertices);        
	}

      }

      // //**** Npu ****
      // if (fabs(MuTree->fMuGenEta) < 2.4) {
      //   histDenominatorNpu->Fill(MuTree->);

      //   //Numerator
      //   if(PassSelection(MuTree,wp)) {
      //     histNumeratorNpu->Fill(MuTree->);        
      //   }

      // }
    }
    if (option==1) {
      //**** PT - ETA ****
      histDenominatorPtEta->Fill(MuTree->fMuPt,MuTree->fMuEta);
      if(PassSelection(MuTree,wp)) {
	histNumeratorPtEta->Fill(MuTree->fMuPt,MuTree->fMuEta);
      }


      //**** PT ****
      histDenominatorPt->Fill(MuTree->fMuPt);

      //Numerator
      if(PassSelection(MuTree,wp)) {
        histNumeratorPt->Fill(MuTree->fMuPt);        
      }


      //**** Eta ****
      if (fabs(MuTree->fMuPt) > 55) {
	histDenominatorEta->Fill(MuTree->fMuEta);

	//Numerator
	if(PassSelection(MuTree,wp)) {
	  histNumeratorEta->Fill(MuTree->fMuEta);        
	}

      }

      //**** Phi ****
      if (fabs(MuTree->fMuEta) < 2.4) {
	histDenominatorPhi->Fill(MuTree->fMuPhi);

	//Numerator
	if(PassSelection(MuTree,wp)) {
	  histNumeratorPhi->Fill(MuTree->fMuPhi);        
	}

      }

      //**** Rho ****
      if (fabs(MuTree->fMuEta) < 2.4) {
	histDenominatorRho->Fill(MuTree->fRho);

	//Numerator
	if(PassSelection(MuTree,wp)) {
	  histNumeratorRho->Fill(MuTree->fRho);        
	}

      }
      //**** Npv ****
      if (fabs(MuTree->fMuEta) < 2.4) {
	histDenominatorNpv->Fill(MuTree->fNVertices);

	//Numerator
	if(PassSelection(MuTree,wp)) {
	  histNumeratorNpv->Fill(MuTree->fNVertices);        
	}

      }

    
    }



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


void MakeMuonEfficiencyPlots(int option = 0) {
  
  if (option == 1) {
    ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/V1p20/MuonNtuple_PromptGenLevel_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fullsim.root", 1, 0, "PromptMuon_TTJets_25ns_Veto_Fullsim");
    // ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/V1p20/MuonNtuple_PromptGenLevel_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fullsim.root", 2, 0, "PromptMuon_TTJets_25ns_Loose_Fullsim");
    ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/V1p20/MuonNtuple_PromptGenLevel_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fullsim.root", 3, 0, "PromptMuon_TTJets_25ns_Tight_Fullsim");
    
    // ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/MuonNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fullsim.root", 1, 1, "FakeMuon_TTJets_25ns_Fullsim_Veto");
    // ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/MuonNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fullsim.root", 2, 1, "FakeMuon_TTJets_25ns_Fullsim_Loose");
    // ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/MuonNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fullsim.root", 3, 1, "FakeMuon_TTJets_25ns_Fullsim_Tight");
  }
  if (option == 2) {
    ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/V1p20/MuonNtuple_PromptGenLevel_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fastsim.root", 1, 0, "PromptMuon_TTJets_25ns_Veto_Fastsim");
    ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/V1p20/MuonNtuple_PromptGenLevel_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fastsim.root", 3, 0, "PromptMuon_TTJets_25ns_Tight_Fastsim");     
  }

  if (option == 3) {
    ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/MuonNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_50ns.root", 100, 1, "PromptMuon_TTJets_50ns_MuTriggerIsoMu20");
    ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/MuonNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_50ns.root", 101, 1, "PromptMuon_TTJets_50ns_MuTriggerIsoMu27");
    ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/MuonNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_50ns.root", 102, 1, "PromptMuon_TTJets_50ns_MuTriggerMu50");
    ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/MuonNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_50ns.root", 110, 1, "PromptMuon_TTJets_50ns_MuTriggerCombined");
    plotMuonTriggerEfficiency();
    return;
  }

  if (option == 4) {
    // ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/V1p18/MuonNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "PromptMuon_TTJets_25_MiniIsolationEACorrRho");
    // ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/V1p18/MuonNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "FakeMuon_TTJets_25_MiniIsolationEACorrRho");   
    // ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/V1p18/MuonNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "PromptMuon_TTJets_25_MiniIsolationEACorrRhoNeutralCentral");
    //  ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/V1p18/MuonNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "FakeMuon_TTJets_25_MiniIsolationEACorrRhoNeutralCentral");   
     plotMuonMiniIsoEfficiency();
  }


  if (option == 10) {   
    ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/work/public/Run2SUSY/MuonNtuple/MuonNtuple_PromptGenLevel_TTJets_20bx25.root", 10, 0, "Muon_NumeratorImprovedIso_DenominatorLooseIDAndIPCut");
    ProduceMuonEfficiencyPlots("/afs/cern.ch/user/s/sixie/work/public/Run2SUSY/MuonNtuple/MuonNtuple_PromptGenLevel_TTJets_20bx25.root", 11, 0, "Muon_NumeratorRelIso0p4_DenominatorLooseIDAndIPCut");
  }

  plotMuonEfficiency();

  if (option == 100) {
    MakeFastsimToFullSimCorrectionFactors();
  }
}
