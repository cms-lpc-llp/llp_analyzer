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
#include <TH1F.h>                
#include <TCanvas.h>                
#include <TGraphAsymmErrors.h>                
#include <TLegend.h>                

#include "RazorAnalyzer/macros/ObjectStudies/EfficiencyUtils.hh"
#include "RazorAnalyzer/include/ElectronTree.h"

#endif

double GetEffArea( double scEta ) {

  double effArea = 0.0;
  //Effective areas below are for the sum of Neutral Hadrons + Photons
  if (fabs(scEta) < 0.8) {
    effArea = 0.0973;
  } else if (fabs(scEta) < 1.3) {
    effArea = 0.0954;
  } else if (fabs(scEta) < 2.0) {
    effArea = 0.0632;	
  } else if (fabs(scEta) < 2.2) {
    effArea = 0.0727;	
  } else {
    effArea = 0.1337;	
  } 

  // //Effective areas using 90% region
  // if (fabs(scEta) < 1.0) {
  //   effArea = 0.1752;
  // } else if (fabs(scEta) < 1.479) {
  //   effArea = 0.1862;
  // } else if (fabs(scEta) < 2.0) {
  //   effArea = 0.1411;	
  // } else if (fabs(scEta) < 2.2) {
  //   effArea = 0.1534;	
  // } else if (fabs(scEta) < 2.3) {
  //   effArea = 0.1903;	
  // } else if (fabs(scEta) < 2.4) {
  //   effArea = 0.2243;	
  // } else {
  //   effArea = 0.2687;	
  // } 
  return effArea;

}



Bool_t passIDMVANonTrigVetoID( ElectronTree *eleTree) {

  Int_t subdet = 0;  
  if (fabs(eleTree->fEleSCEta) < 0.8) subdet = 0;
  else if (fabs(eleTree->fEleSCEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (eleTree->fElePt > 10.0) ptBin = 1;

  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = -0.11;
  if (subdet == 1 && ptBin == 0) MVACut = -0.55;
  if (subdet == 2 && ptBin == 0) MVACut = -0.60;
  if (subdet == 0 && ptBin == 1) MVACut = -0.16;
  if (subdet == 1 && ptBin == 1) MVACut = -0.65;
  if (subdet == 2 && ptBin == 1) MVACut = -0.74;

  bool pass = false;
  if (eleTree->fIDMVANonTrig > MVACut 
      && fabs(eleTree->fEleIP3dSig) < 4
      ) {
    pass = true;
  }   
  return pass;

}


bool PassVetoMiniIso( ElectronTree* eleTree ) {

  bool pass = false;
  double dr = fmax(0.05,fmin(0.2, 10/eleTree->fElePt));
  pass = bool( (eleTree->fMiniIso - eleTree->fRho*GetEffArea(eleTree->fEleSCEta)*pow(dr/0.3,2)) / eleTree->fElePt < 0.1 );
  return pass;
}



TH2F* compareEfficiency( TH2F *h1, TH2F *h2, bool useSigmas = false) {

  TH2F *result = (TH2F*)h1->Clone();
  for (int i=0; i<h1->GetXaxis()->GetNbins()+2; i++) {
    for (int j=0; j<h1->GetYaxis()->GetNbins()+2; j++) {
      double ratio = 0;
      if (h2->GetBinContent(i,j) > 0) ratio = h1->GetBinContent(i,j) / h2->GetBinContent(i,j);
      double err = 0;
      if (h1->GetBinContent(i,j) > 0  && h2->GetBinContent(i,j)) 
	err = ratio*sqrt( pow( h1->GetBinError(i,j)/h1->GetBinContent(i,j) , 2) + pow( h2->GetBinError(i,j)/h2->GetBinContent(i,j) , 2));

      if (useSigmas) {
	result->SetBinContent( i,j, (ratio - 1)/err );
      } else {
	result->SetBinContent( i,j,ratio -1 );
      }
    }
  }

  return result;
}



void compareEfficiencies() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileTTJetsAMCAtNLO = new TFile("Efficiency_PromptElectron_TTJets_aMCAtNLO_25ns_MiniIsoCut.root","READ");
  TFile *fileTTJetsMadgraph = new TFile("Efficiency_PromptElectron_TTJets_25ns_MiniIsoCut.root","READ");
  TFile *fileTTJetsHT2500ToInf = new TFile("Efficiency_PromptElectron_TTJetsHT2500ToInf_25ns_MiniIsoCut.root","READ");
  TFile *fileTTJetsHT600To800 = new TFile("Efficiency_PromptElectron_TTJetsHT600To800_25ns_MiniIsoCut.root","READ");
  TFile *fileWJetsHT2500ToInf = new TFile("Efficiency_PromptElectron_WJetsHT2500ToInf_25ns_MiniIsoCut.root","READ");
  TFile *fileWJetsHT400To600 = new TFile("Efficiency_PromptElectron_WJetsHT400To600_25ns_MiniIsoCut.root","READ");
  TFile *fileDYJetsHT600ToInf = new TFile("Efficiency_PromptElectron_DYJetsHT600ToInf_25ns_MiniIsoCut.root","READ");
  TFile *fileDYJetsHT200To400 = new TFile("Efficiency_PromptElectron_DYJetsHT200To400_25ns_MiniIsoCut.root","READ");

  TFile *fileTTJetsMadgraphAll = new TFile("Efficiency_PromptElectron_TTJetsMadgraphAll_25ns_MiniIsoCut.root","READ");
  TFile *fileWJetsMadgraphAll = new TFile("Efficiency_PromptElectron_WJetsMadgraphAll_25ns_MiniIsoCut.root","READ");
  TFile *fileDYJetsMadgraphAll = new TFile("Efficiency_PromptElectron_DYJetsMadgraphAll_25ns_MiniIsoCut.root","READ");


  TH2F* effTTJetsAMCAtNLO = (TH2F*)fileTTJetsAMCAtNLO->Get("Efficiency_PtActivity");
  TH2F* effTTJetsMadgraph = (TH2F*)fileTTJetsMadgraph->Get("Efficiency_PtActivity");
  TH2F* effTTJetsHT2500ToInf = (TH2F*)fileTTJetsHT2500ToInf->Get("Efficiency_PtActivity");
  TH2F* effTTJetsHT600To800 = (TH2F*)fileTTJetsHT600To800->Get("Efficiency_PtActivity");
  TH2F* effWJetsHT2500ToInf = (TH2F*)fileWJetsHT2500ToInf->Get("Efficiency_PtActivity");
  TH2F* effWJetsHT400To600 = (TH2F*)fileWJetsHT400To600->Get("Efficiency_PtActivity");
  TH2F* effDYJetsHT600ToInf = (TH2F*)fileDYJetsHT600ToInf->Get("Efficiency_PtActivity");
  TH2F* effDYJetsHT200To400 = (TH2F*)fileDYJetsHT200To400->Get("Efficiency_PtActivity");
  TH2F* effTTJetsMadgraphAll = (TH2F*)fileTTJetsMadgraphAll->Get("Efficiency_PtActivity");
  TH2F* effWJetsMadgraphAll = (TH2F*)fileWJetsMadgraphAll->Get("Efficiency_PtActivity");
  TH2F* effDYJetsMadgraphAll = (TH2F*)fileDYJetsMadgraphAll->Get("Efficiency_PtActivity");

  TH2F* effDifferenceRatio_TTJets_DifferentHTBins = compareEfficiency( effTTJetsHT600To800 , effTTJetsHT2500ToInf , false);
  TH2F* effDifferenceSigmas_TTJets_DifferentHTBins = compareEfficiency( effTTJetsHT600To800 , effTTJetsHT2500ToInf, true );
  TH2F* effDifferenceRatio_WJets_DifferentHTBins = compareEfficiency( effWJetsHT400To600 , effWJetsHT2500ToInf , false);
  TH2F* effDifferenceSigmas_WJets_DifferentHTBins = compareEfficiency( effWJetsHT400To600 , effWJetsHT2500ToInf, true );
  TH2F* effDifferenceRatio_DYJets_DifferentHTBins = compareEfficiency( effDYJetsHT200To400 , effDYJetsHT600ToInf , false);
  TH2F* effDifferenceSigmas_DYJets_DifferentHTBins = compareEfficiency( effDYJetsHT200To400 , effDYJetsHT600ToInf, true );

  TH2F* effDifferenceRatio_TTJetsVsDYJets = compareEfficiency( effTTJetsMadgraphAll , effDYJetsMadgraphAll , false);
  TH2F* effDifferenceSigmas_TTJetsVsDYJets = compareEfficiency( effTTJetsMadgraphAll, effDYJetsMadgraphAll, true );
  TH2F* effDifferenceRatio_WJetsVsDYJets = compareEfficiency( effWJetsMadgraphAll , effDYJetsMadgraphAll , false);
  TH2F* effDifferenceSigmas_WJetsVsDYJets = compareEfficiency( effWJetsMadgraphAll, effDYJetsMadgraphAll, true );

  TH2F* effDifferenceRatio_NLO = compareEfficiency( effTTJetsAMCAtNLO , effTTJetsMadgraph , false);
  TH2F* effDifferenceSigmas_NLO = compareEfficiency( effTTJetsAMCAtNLO , effTTJetsMadgraph, true );

  //color palette
  const Int_t Number = 5;
  Double_t Red[Number]    =  { 0.00, 0.85, 1.00, 1.00, 1.00};
  Double_t Green[Number]  =  { 0.00, 0.85, 1.00, 0.85, 0.00};
  Double_t Blue[Number]   =  { 1.00, 1.00, 1.00, 0.85, 0.00};
  Double_t Length[Number] =  { 0.00, 0.2 , 0.5,  0.8 , 1.00};
  Int_t nb=200;
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
  gStyle->SetPaintTextFormat("4.2f");

  //**************************************************
  //TTJets HT Bins
  //**************************************************

  cv = new TCanvas("cv","cv", 800,600);

  effDifferenceRatio_TTJets_DifferentHTBins->Draw("colz,text");
  effDifferenceRatio_TTJets_DifferentHTBins->SetStats(0);
  effDifferenceRatio_TTJets_DifferentHTBins->SetMaximum(0.2);
  effDifferenceRatio_TTJets_DifferentHTBins->SetMinimum(-0.2);
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("MiniIsoEfficiency_TTJets_CompareHTBins_Ratio.gif");


  cv = new TCanvas("cv","cv", 800,600);

  effDifferenceSigmas_TTJets_DifferentHTBins->Draw("colz,text");
  effDifferenceSigmas_TTJets_DifferentHTBins->SetStats(0);
  effDifferenceSigmas_TTJets_DifferentHTBins->SetMaximum(5);
  effDifferenceSigmas_TTJets_DifferentHTBins->SetMinimum(-5);
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("MiniIsoEfficiency_TTJets_CompareHTBins_Sigmas.gif");



  //**************************************************
  //TTJets NLO
  //**************************************************

  cv = new TCanvas("cv","cv", 800,600);

  effDifferenceRatio_NLO->Draw("colz,text");
  effDifferenceRatio_NLO->SetStats(0);
  effDifferenceRatio_NLO->SetMaximum(0.2);
  effDifferenceRatio_NLO->SetMinimum(-0.2);
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("MiniIsoEfficiency_NLO_Ratio.gif");


  cv = new TCanvas("cv","cv", 800,600);

  effDifferenceSigmas_NLO->Draw("colz,text");
  effDifferenceSigmas_NLO->SetStats(0);
  effDifferenceSigmas_NLO->SetMaximum(5);
  effDifferenceSigmas_NLO->SetMinimum(-5);
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("MiniIsoEfficiency_NLO_Sigmas.gif");

  //**************************************************
  //WJets HT Bins
  //**************************************************

  cv = new TCanvas("cv","cv", 800,600);

  effDifferenceRatio_WJets_DifferentHTBins->Draw("colz,text");
  effDifferenceRatio_WJets_DifferentHTBins->SetStats(0);
  effDifferenceRatio_WJets_DifferentHTBins->SetMaximum(0.2);
  effDifferenceRatio_WJets_DifferentHTBins->SetMinimum(-0.2);
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("MiniIsoEfficiency_WJets_CompareHTBins_Ratio.gif");


  cv = new TCanvas("cv","cv", 800,600);

  effDifferenceSigmas_WJets_DifferentHTBins->Draw("colz,text");
  effDifferenceSigmas_WJets_DifferentHTBins->SetStats(0);
  effDifferenceSigmas_WJets_DifferentHTBins->SetMaximum(5);
  effDifferenceSigmas_WJets_DifferentHTBins->SetMinimum(-5);
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("MiniIsoEfficiency_WJets_CompareHTBins_Sigmas.gif");


  //**************************************************
  //DYJets HT Bins
  //**************************************************

  cv = new TCanvas("cv","cv", 800,600);

  effDifferenceRatio_DYJets_DifferentHTBins->Draw("colz,text");
  effDifferenceRatio_DYJets_DifferentHTBins->SetStats(0);
  effDifferenceRatio_DYJets_DifferentHTBins->SetMaximum(0.2);
  effDifferenceRatio_DYJets_DifferentHTBins->SetMinimum(-0.2);
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("MiniIsoEfficiency_DYJets_CompareHTBins_Ratio.gif");


  cv = new TCanvas("cv","cv", 800,600);

  effDifferenceSigmas_DYJets_DifferentHTBins->Draw("colz,text");
  effDifferenceSigmas_DYJets_DifferentHTBins->SetStats(0);
  effDifferenceSigmas_DYJets_DifferentHTBins->SetMaximum(5);
  effDifferenceSigmas_DYJets_DifferentHTBins->SetMinimum(-5);
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("MiniIsoEfficiency_DYJets_CompareHTBins_Sigmas.gif");


  //**************************************************
  //Compare TTJets vs DY
  //**************************************************

  cv = new TCanvas("cv","cv", 800,600);

  effDifferenceRatio_TTJetsVsDYJets->Draw("colz,text");
  effDifferenceRatio_TTJetsVsDYJets->SetStats(0);
  effDifferenceRatio_TTJetsVsDYJets->SetMaximum(0.2);
  effDifferenceRatio_TTJetsVsDYJets->SetMinimum(-0.2);
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("MiniIsoEfficiency_TTJetsVsDYJets_Ratio.gif");


  cv = new TCanvas("cv","cv", 800,600);

  effDifferenceSigmas_TTJetsVsDYJets->Draw("colz,text");
  effDifferenceSigmas_TTJetsVsDYJets->SetStats(0);
  effDifferenceSigmas_TTJetsVsDYJets->SetMaximum(5);
  effDifferenceSigmas_TTJetsVsDYJets->SetMinimum(-5);
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("MiniIsoEfficiency_TTJetsVsDYJets_Sigmas.gif");

  //**************************************************
  //Compare WJets vs DY
  //**************************************************

  cv = new TCanvas("cv","cv", 800,600);

  effDifferenceRatio_WJetsVsDYJets->Draw("colz,text");
  effDifferenceRatio_WJetsVsDYJets->SetStats(0);
  effDifferenceRatio_WJetsVsDYJets->SetMaximum(0.2);
  effDifferenceRatio_WJetsVsDYJets->SetMinimum(-0.2);
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("MiniIsoEfficiency_WJetsVsDYJets_Ratio.gif");


  cv = new TCanvas("cv","cv", 800,600);

  effDifferenceSigmas_WJetsVsDYJets->Draw("colz,text");
  effDifferenceSigmas_WJetsVsDYJets->SetStats(0);
  effDifferenceSigmas_WJetsVsDYJets->SetMaximum(5);
  effDifferenceSigmas_WJetsVsDYJets->SetMinimum(-5);
  cv->SetLogy();
  cv->SetLogx();
  cv->SaveAs("MiniIsoEfficiency_WJetsVsDYJets_Sigmas.gif");



}


//=== MAIN MACRO ================================================================================================= 

void MakeElectronEfficiencyVsActivityPlots(const string inputfile, int wp, int option = -1, string label = "") {
 
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************

  TH2F *histDenominatorPtActivity = new TH2F ("histDenominatorPtActivity",";Electron p_{T} [GeV/c] ; Electron Activity / pT ; Number of Events", 300, 0 , 300, 500, 0, 50);
  TH2F *histNumeratorPtActivity = new TH2F ("histNumeratorPtActivity",";Electronp p_{T} [GeV/c] ; Electron Activity / pT ; Number of Events", 300, 0 , 300, 500, 0, 50);

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  ElectronTree *EleTree = new ElectronTree;
  EleTree->LoadTree(inputfile.c_str());
  EleTree->InitTree(ElectronTree::kEleTreeLight);

  cout << "Total Entries: " << EleTree->tree_->GetEntries() << "\n";
  int nentries = EleTree->tree_->GetEntries();
  for(UInt_t ientry=0; ientry < EleTree->tree_->GetEntries(); ientry++) {       	
    EleTree->tree_->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //Cuts
    if (EleTree->fEleGenPt < 5) continue;
    if (abs(EleTree->fEleGenEta) > 2.5) continue;
    if (!(EleTree->fElePt > 0)) continue;
    if (!passIDMVANonTrigVetoID(EleTree)) continue;


    //**** PT - ETA ****
    histDenominatorPtActivity->Fill(EleTree->fElePt,EleTree->fActivity);
    if( PassVetoMiniIso( EleTree ) ) {
      histNumeratorPtActivity->Fill(EleTree->fElePt,EleTree->fActivity);
    }
    
    
  }


  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency Plots
  //==============================================================================================================
  vector<double> ptBins;
  ptBins.push_back(5);
  ptBins.push_back(10);
  ptBins.push_back(15);
  ptBins.push_back(20);
  ptBins.push_back(30);
  ptBins.push_back(50);
  ptBins.push_back(100);
  ptBins.push_back(200);
  ptBins.push_back(300);
  vector<double> ActivityBins;
  ActivityBins.push_back(0);
  ActivityBins.push_back(0.1);
  ActivityBins.push_back(0.2);
  ActivityBins.push_back(0.3);
  ActivityBins.push_back(0.4);  
  ActivityBins.push_back(0.5);  
  ActivityBins.push_back(1.0); 
  ActivityBins.push_back(2.0);
  ActivityBins.push_back(5.0);
  ActivityBins.push_back(10.0);
  ActivityBins.push_back(20.0);
  ActivityBins.push_back(50.0);

  TH2F *efficiency_PtActivity = createEfficiencyHist2D(histNumeratorPtActivity, histDenominatorPtActivity, "Efficiency_PtActivity" , ptBins , ActivityBins);  


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;



  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("Efficiency"+Label+".root").c_str(), "UPDATE");
  file->cd();
  file->WriteTObject(efficiency_PtActivity, "Efficiency_PtActivity", "WriteDelete");

  file->Close();
  delete file;       

}

void ElectronActivityStudy( int option = 1) {

  if (option == 1) {
   
    MakeElectronEfficiencyVsActivityPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p17/ElectronNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", 1, 0, "PromptElectron_TTJets_25ns_MiniIsoCut");
    // MakeElectronEfficiencyVsActivityPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p17/ElectronNtuple_Prompt_TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", 1, 0, "PromptElectron_TTJetsHT2500ToInf_25ns_MiniIsoCut");
    // MakeElectronEfficiencyVsActivityPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p17/ElectronNtuple_Prompt_TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", 1, 0, "PromptElectron_TTJetsHT600To800_25ns_MiniIsoCut");
    // MakeElectronEfficiencyVsActivityPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p17/ElectronNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root", 1, 0, "PromptElectron_TTJets_aMCAtNLO_25ns_MiniIsoCut");

    // MakeElectronEfficiencyVsActivityPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p17/ElectronNtuple_Prompt_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", 1, 0, "PromptElectron_WJetsHT2500ToInf_25ns_MiniIsoCut");
    // MakeElectronEfficiencyVsActivityPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p17/ElectronNtuple_Prompt_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", 1, 0, "PromptElectron_WJetsHT400To600_25ns_MiniIsoCut");

    // MakeElectronEfficiencyVsActivityPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p17/ElectronNtuple_Prompt_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", 1, 0, "PromptElectron_DYJetsHT600ToInf_25ns_MiniIsoCut");
    // MakeElectronEfficiencyVsActivityPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p17/ElectronNtuple_Prompt_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", 1, 0, "PromptElectron_DYJetsHT200To400_25ns_MiniIsoCut");

    // MakeElectronEfficiencyVsActivityPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p17/ElectronNtuple_Prompt_TTJetsMadgraphAll.root", 1, 0, "PromptElectron_TTJetsMadgraphAll_25ns_MiniIsoCut");
    // MakeElectronEfficiencyVsActivityPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p17/ElectronNtuple_Prompt_WJetsMadgraphAll.root", 1, 0, "PromptElectron_WJetsMadgraphAll_25ns_MiniIsoCut");
    // MakeElectronEfficiencyVsActivityPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p17/ElectronNtuple_Prompt_DYJetsMadgraphAll.root", 1, 0, "PromptElectron_DYJetsMadgraphAll_25ns_MiniIsoCut");

  }

  if (option ==0) {
    compareEfficiencies();
  }



}
