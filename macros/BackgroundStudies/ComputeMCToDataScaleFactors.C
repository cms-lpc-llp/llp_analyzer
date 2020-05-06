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





//=== MAIN MACRO ================================================================================================= 


void PlotMCToDataScaleFactors_TTBarSingleLepton() {
  
  TFile *TTBarDileptonFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/ScaleFactors/Run1/TTBarDileptonScaleFactors.root", "READ");
  TFile *TTBarSingleMuFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/TTBarSingleLeptonControlRegionPlots_MR300Rsq0p15_OneMediumBTag_SingleMu.root", "READ");
  TFile *TTBarSingleEleFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/TTBarSingleLeptonControlRegionPlots_MR300Rsq0p15_OneMediumBTag_SingleEle.root", "READ");

  TH2D* TTBarDileptonHist = (TH2D*)TTBarDileptonFile->Get("TTBarDileptonScaleFactor");
  TH2D* TTBarSingleMuHist = (TH2D*)TTBarSingleMuFile->Get("MCToDataScaleFactor_MRVsRsq");
  TH2D* TTBarSingleEleHist = (TH2D*)TTBarSingleEleFile->Get("MCToDataScaleFactor_MRVsRsq");

  
  TCanvas *cv = 0;
  cv = new TCanvas("cv","cv",800,600);

  TTBarSingleMuHist->SetMaximum(1.5);
  TTBarSingleMuHist->SetMinimum(0.0);
  TTBarSingleMuHist->Draw("colz");
  TTBarSingleMuHist->SetStats(false);
  TTBarSingleMuHist->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetRightMargin(0.15);
  cv->SaveAs("TTBarSingleMuScaleFactor.gif");



  cv = new TCanvas("cv","cv",800,600);

  TTBarSingleEleHist->SetMaximum(1.5);
  TTBarSingleEleHist->SetMinimum(0.0);
  TTBarSingleEleHist->Draw("colz");
  TTBarSingleEleHist->SetStats(false);
  TTBarSingleEleHist->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetRightMargin(0.15);
  cv->SaveAs("TTBarSingleEleScaleFactor.gif");



  TH2D* CompareTTBarSingleEleSingleMu = (TH2D*)TTBarSingleMuHist->Clone("CompareTTBarSingleEleSingleMu");
  for (int i=0; i<CompareTTBarSingleEleSingleMu->GetXaxis()->GetNbins()+1;i++) {
    for (int j=0; j<CompareTTBarSingleEleSingleMu->GetYaxis()->GetNbins()+1;j++) {
      CompareTTBarSingleEleSingleMu->SetBinContent(i,j,  (TTBarSingleMuHist->GetBinContent(i,j) - TTBarSingleEleHist->GetBinContent(i,j)) / sqrt( pow(TTBarSingleMuHist->GetBinError(i,j), 2) + pow(TTBarSingleEleHist->GetBinError(i,j), 2)));
    }
  }
  cv = new TCanvas("cv","cv",800,600);

  CompareTTBarSingleEleSingleMu->SetMaximum(5);
  CompareTTBarSingleEleSingleMu->SetMinimum(-5);
  CompareTTBarSingleEleSingleMu->Draw("colz");
  CompareTTBarSingleEleSingleMu->SetStats(false);
  CompareTTBarSingleEleSingleMu->GetZaxis()->SetTitle("Difference in #sigma");
  cv->SetRightMargin(0.15);
  cv->SaveAs("CompareTTBarSingleEleSingleMu.gif");


  //********************************************************************************************
  //Combine SingleMu and SingleEle
  //********************************************************************************************
  TH2D* TTBarSingleLep = (TH2D*)TTBarSingleMuHist->Clone("TTBarSingleLep");
  TH2D* TTBarSingleMuDataMinusBkgHist = (TH2D*)TTBarSingleMuFile->Get("DataMinusBkg_MRVsRsq");
  TH2D* TTBarSingleEleDataMinusBkgHist = (TH2D*)TTBarSingleEleFile->Get("DataMinusBkg_MRVsRsq");
  TH2D* TTBarSingleMuMCHist = (TH2D*)TTBarSingleMuFile->Get("histMRVsRsq_TTJets");
  TH2D* TTBarSingleEleMCHist = (TH2D*)TTBarSingleEleFile->Get("histMRVsRsq_TTJets");

  assert(TTBarSingleLep);
  assert(TTBarSingleMuDataMinusBkgHist);
  assert(TTBarSingleEleDataMinusBkgHist);
  assert(TTBarSingleMuMCHist );
  assert(TTBarSingleEleMCHist );

  for (int i=0; i<CompareTTBarSingleEleSingleMu->GetXaxis()->GetNbins()+1;i++) {
    for (int j=0; j<CompareTTBarSingleEleSingleMu->GetYaxis()->GetNbins()+1;j++) {

      double mc = TTBarSingleMuMCHist->GetBinContent(i,j) + TTBarSingleEleMCHist->GetBinContent(i,j);
      double mcErr = sqrt( pow(TTBarSingleMuMCHist->GetBinError(i,j),2) + pow(TTBarSingleEleMCHist->GetBinError(i,j),2));

      double data = TTBarSingleMuDataMinusBkgHist->GetBinContent(i,j) + TTBarSingleEleDataMinusBkgHist->GetBinContent(i,j);
      double dataErr = sqrt( pow(TTBarSingleMuDataMinusBkgHist->GetBinError(i,j),2) + pow(TTBarSingleEleDataMinusBkgHist->GetBinError(i,j),2));

      // cout << "bin " << i << j << " : " << TTBarSingleMuDataMinusBkgHist->GetBinContent(i,j) << " " << TTBarSingleEleDataMinusBkgHist->GetBinContent(i,j)  << " | "
      // 	   << TTBarSingleMuMCHist->GetBinContent(i,j) << " " <<  TTBarSingleEleMCHist->GetBinContent(i,j) << " "
      // 	   << "\n";
      if (i != 0 && j!=0) {
	cout << "bin " << i << j << " : " << data << " +/- " << dataErr << " , " << mc << " +/- " << mcErr 
	     << " | " << data/mc << " +/- " << (data/mc)*sqrt( pow(mcErr/mc,2)+pow(dataErr/data,2) ) << "\n";

	TTBarSingleLep->SetBinContent(i,j, data/mc);
	TTBarSingleLep->SetBinError(i,j, (data/mc)*sqrt( pow(mcErr/mc,2)+pow(dataErr/data,2) ) );
      } else {
	TTBarSingleLep->SetBinContent(i,j, 1.0);
	TTBarSingleLep->SetBinError(i,j, 0);
      }
    }
  }
  cv = new TCanvas("cv","cv",800,600);
  TTBarSingleLep->SetMaximum(1.5);
  TTBarSingleLep->SetMinimum(0.0);
  TTBarSingleLep->Draw("colz");
  TTBarSingleLep->SetStats(false);
  TTBarSingleLep->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetRightMargin(0.15);
  cv->SetLogx();
  cv->SetLogy();
  cv->SaveAs("TTBarSingleLep.gif");
  


  //********************************************************************************************
  //Compare Dilepton and Single Lepton Scale Factor For TTBar
  //********************************************************************************************

  TH2D* CompareTTBarSingleLeptonWithDileptonSigmas = (TH2D*)TTBarDileptonHist->Clone("CompareTTBarSingleLeptonWithDileptonSigmas");
  TH2D* CompareTTBarSingleLeptonWithDileptonFraction = (TH2D*)TTBarDileptonHist->Clone("CompareTTBarSingleLeptonWithDileptonFraction");
  for (int i=0; i<CompareTTBarSingleLeptonWithDileptonSigmas->GetXaxis()->GetNbins()+1;i++) {
    for (int j=0; j<CompareTTBarSingleLeptonWithDileptonSigmas->GetYaxis()->GetNbins()+1;j++) {
      if (i!=0 && j!=0) {
	CompareTTBarSingleLeptonWithDileptonSigmas->SetBinContent(i,j, (TTBarDileptonHist->GetBinContent(i,j) - TTBarSingleLep->GetBinContent(i,j)) / sqrt( pow(TTBarDileptonHist->GetBinError(i,j), 2) + pow(TTBarSingleLep->GetBinError(i,j), 2)));
	CompareTTBarSingleLeptonWithDileptonFraction->SetBinContent(i,j, (TTBarDileptonHist->GetBinContent(i,j) - TTBarSingleLep->GetBinContent(i,j)) / TTBarDileptonHist->GetBinContent(i,j) );
	
	cout << "Bin " << i << " " << j << " : " << TTBarDileptonHist->GetBinContent(i,j) << " +/- " << TTBarDileptonHist->GetBinError(i,j) << " , " << TTBarSingleLep->GetBinContent(i,j) << " +/- " << TTBarSingleLep->GetBinError(i,j) << " : " 
	     << (TTBarDileptonHist->GetBinContent(i,j) - TTBarSingleLep->GetBinContent(i,j)) / sqrt( pow(TTBarDileptonHist->GetBinError(i,j), 2) + pow(TTBarSingleLep->GetBinError(i,j), 2)) << " " << (TTBarDileptonHist->GetBinContent(i,j) - TTBarSingleLep->GetBinContent(i,j)) / TTBarDileptonHist->GetBinContent(i,j) 
	     << " \n";
      } else {
	CompareTTBarSingleLeptonWithDileptonSigmas->SetBinContent(i,j, 0);
	CompareTTBarSingleLeptonWithDileptonFraction->SetBinContent(i,j, 0);	
      }

    }
  }
  cv = new TCanvas("cv","cv",800,600);
  CompareTTBarSingleLeptonWithDileptonSigmas->SetMaximum(5);
  CompareTTBarSingleLeptonWithDileptonSigmas->SetMinimum(-5);
  CompareTTBarSingleLeptonWithDileptonSigmas->SetStats(false);
  CompareTTBarSingleLeptonWithDileptonSigmas->Draw("colz,text");
  CompareTTBarSingleLeptonWithDileptonSigmas->GetZaxis()->SetTitle("Difference in #sigma");
  CompareTTBarSingleLeptonWithDileptonSigmas->Draw("colz,text");
  cv->SetRightMargin(0.15);
  cv->SetLogy();
  cv->SaveAs("CompareTTBarSingleLeptonWithDileptonSigmas.gif");

 
  cv = new TCanvas("cv","cv",800,600);
  CompareTTBarSingleLeptonWithDileptonFraction->SetMaximum(1.0);
  CompareTTBarSingleLeptonWithDileptonFraction->SetMinimum(-1.0);
  CompareTTBarSingleLeptonWithDileptonFraction->SetStats(false);
  CompareTTBarSingleLeptonWithDileptonFraction->Draw("colz,text");
  CompareTTBarSingleLeptonWithDileptonFraction->GetZaxis()->SetTitle("Fractional Difference");
  CompareTTBarSingleLeptonWithDileptonFraction->Draw("colz,text");
  gPad->Update();
  cv->SetRightMargin(0.15);
  cv->SetLogy();
  cv->SaveAs("CompareTTBarSingleLeptonWithDileptonFraction.gif");

  return;

  TFile *file = TFile::Open("TTBarSingleLeptonScaleFactors.root", "RECREATE");
  file->cd();  
  file->WriteTObject(TTBarSingleLep, "TTBarSingleLeptonScaleFactor", "WriteDelete");
  file->WriteTObject(TTBarSingleMuHist, "TTBarSingleMuScaleFactor", "WriteDelete");
  file->WriteTObject(TTBarSingleEleHist, "TTBarSingleEleScaleFactor", "WriteDelete");
  file->Close();
  delete file;


}


void PlotMCToDataScaleFactors_WJets( ) {
  

  TFile *WJetsSingleMuFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/WJetsControlRegionPlots_MR300Rsq0p15_ZeroMediumBTag_SingleMu.root", "READ");
  TFile *WJetsSingleEleFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/WJetsControlRegionPlots_MR300Rsq0p15_ZeroMediumBTag_SingleEle.root", "READ");

  TH2D* WJetsSingleMuHist = (TH2D*)WJetsSingleMuFile->Get("MCToDataScaleFactor_MRVsRsq");
  TH2D* WJetsSingleEleHist = (TH2D*)WJetsSingleEleFile->Get("MCToDataScaleFactor_MRVsRsq");

  
  TCanvas *cv = 0;
  cv = new TCanvas("cv","cv",800,600);
  WJetsSingleMuHist->SetMaximum(1.5);
  WJetsSingleMuHist->SetMinimum(0.0);
  WJetsSingleMuHist->Draw("colz");
  WJetsSingleMuHist->SetStats(false);
  WJetsSingleMuHist->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetRightMargin(0.15);
  cv->SaveAs("WJetsSingleMuScaleFactor.gif");


  cv = new TCanvas("cv","cv",800,600);

  WJetsSingleEleHist->SetMaximum(1.5);
  WJetsSingleEleHist->SetMinimum(0.0);
  WJetsSingleEleHist->Draw("colz");
  WJetsSingleEleHist->SetStats(false);
  WJetsSingleEleHist->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetRightMargin(0.15);
  cv->SaveAs("WJetsSingleEleScaleFactor.gif");


  
  TH2D* CompareWJetsSingleEleSingleMu = (TH2D*)WJetsSingleMuHist->Clone("CompareWJetsSingleEleSingleMu");
  for (int i=0; i<CompareWJetsSingleEleSingleMu->GetXaxis()->GetNbins()+1;i++) {
    for (int j=0; j<CompareWJetsSingleEleSingleMu->GetYaxis()->GetNbins()+1;j++) {
      CompareWJetsSingleEleSingleMu->SetBinContent(i,j, (WJetsSingleMuHist->GetBinContent(i,j) - WJetsSingleEleHist->GetBinContent(i,j)) / sqrt( pow(WJetsSingleMuHist->GetBinError(i,j), 2) + pow(WJetsSingleEleHist->GetBinError(i,j), 2)));

      cout << "bin " << i << " " << j << " : " << WJetsSingleMuHist->GetBinContent(i,j) << " +/- " << WJetsSingleMuHist->GetBinError(i,j) << " , " << WJetsSingleEleHist->GetBinContent(i,j) << " +/- " << WJetsSingleEleHist->GetBinError(i,j) << " | " << (WJetsSingleMuHist->GetBinContent(i,j) - WJetsSingleEleHist->GetBinContent(i,j)) / sqrt( pow(WJetsSingleMuHist->GetBinError(i,j), 2) + pow(WJetsSingleEleHist->GetBinError(i,j), 2)) << "\n";

    }
  }
  cv = new TCanvas("cv","cv",800,600);

  CompareWJetsSingleEleSingleMu->SetMaximum(5);
  CompareWJetsSingleEleSingleMu->SetMinimum(-5);
  CompareWJetsSingleEleSingleMu->Draw("colz");
  CompareWJetsSingleEleSingleMu->SetStats(false);
  CompareWJetsSingleEleSingleMu->GetZaxis()->SetTitle("Difference in #sigma");
  cv->SetRightMargin(0.15);
  cv->SaveAs("CompareWJetsSingleEleSingleMu.gif");


  //********************************************************************************************
  //Combine SingleMu and SingleEle
  //********************************************************************************************
  TH2D* WJetsSingleLep = (TH2D*)WJetsSingleMuHist->Clone("WJetsSingleLep");
  TH2D* WJetsSingleMuDataMinusBkgHist = (TH2D*)WJetsSingleMuFile->Get("DataMinusBkg_MRVsRsq");
  TH2D* WJetsSingleEleDataMinusBkgHist = (TH2D*)WJetsSingleEleFile->Get("DataMinusBkg_MRVsRsq");
  TH2D* WJetsSingleMuMCHist = (TH2D*)WJetsSingleMuFile->Get("histMRVsRsq_WJets");
  TH2D* WJetsSingleEleMCHist = (TH2D*)WJetsSingleEleFile->Get("histMRVsRsq_WJets");

  for (int i=0; i<CompareWJetsSingleEleSingleMu->GetXaxis()->GetNbins()+1;i++) {
    for (int j=0; j<CompareWJetsSingleEleSingleMu->GetYaxis()->GetNbins()+1;j++) {

      double mc = WJetsSingleMuMCHist->GetBinContent(i,j) + WJetsSingleEleMCHist->GetBinContent(i,j);
      double mcErr = sqrt( pow(WJetsSingleMuMCHist->GetBinError(i,j),2) + pow(WJetsSingleEleMCHist->GetBinError(i,j),2));

      double data = WJetsSingleMuDataMinusBkgHist->GetBinContent(i,j) + WJetsSingleEleDataMinusBkgHist->GetBinContent(i,j);
      double dataErr = sqrt( pow(WJetsSingleMuDataMinusBkgHist->GetBinError(i,j),2) + pow(WJetsSingleEleDataMinusBkgHist->GetBinError(i,j),2));

      cout << "bin " << i << j << " : " << WJetsSingleMuDataMinusBkgHist->GetBinContent(i,j) << " " << WJetsSingleEleDataMinusBkgHist->GetBinContent(i,j)  << " | "
	   << WJetsSingleMuMCHist->GetBinContent(i,j) << " " <<  WJetsSingleEleMCHist->GetBinContent(i,j) << " "
	   << "\n";
     // cout << "bin " << i << j << " : " << data << " +/- " << dataErr << " , " << mc << " +/- " << mcErr << "\n";

      WJetsSingleLep->SetBinContent(i,j, data/mc);
      WJetsSingleLep->SetBinError(i,j, (data/mc)*sqrt( pow(mcErr/mc,2)+pow(dataErr/data,2) ) );
    }
  }
  cv = new TCanvas("cv","cv",800,600);
  WJetsSingleLep->SetMaximum(1.5);
  WJetsSingleLep->SetMinimum(0.0);
  WJetsSingleLep->Draw("colz");
  WJetsSingleLep->SetStats(false);
  WJetsSingleLep->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetRightMargin(0.15);
  cv->SetLogx();
  cv->SetLogy();
  cv->SaveAs("WJetsSingleLep.gif");
  

  TFile *file = TFile::Open("WJetsSingleLeptonScaleFactors.root", "RECREATE");
  file->cd();  
  file->WriteTObject(WJetsSingleLep, "WJetsSingleLeptonScaleFactor", "WriteDelete");
  file->WriteTObject(WJetsSingleMuHist, "WJetsSingleMuScaleFactor", "WriteDelete");
  file->WriteTObject(WJetsSingleEleHist, "WJetsSingleEleScaleFactor", "WriteDelete");
  file->Close();
  delete file;

}





void PlotMCToDataScaleFactors_ZToLL( ) {
  
  // string Label = "";
  // if (label != "") Label = "_" + label;


  TFile *ZToEEFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/ZToLLControlRegionPlots_MR300Rsq0p05_ee.root", "READ");
  TFile *ZToMMFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/ZToLLControlRegionPlots_MR300Rsq0p05_mumu.root", "READ");

  TH2D* ZToEEHist = (TH2D*)ZToEEFile->Get("MCToDataScaleFactor_MRVsRsq");
  TH2D* ZToMMHist = (TH2D*)ZToMMFile->Get("MCToDataScaleFactor_MRVsRsq");

  
  TCanvas *cv = 0;
  cv = new TCanvas("cv","cv",800,600);

  ZToEEHist->SetMaximum(1.5);
  ZToEEHist->SetMinimum(0.0);
  ZToEEHist->Draw("colz");
  ZToEEHist->SetStats(false);
  ZToEEHist->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetRightMargin(0.15);
  cv->SaveAs("ZToEEScaleFactor.gif");


  cv = new TCanvas("cv","cv",800,600);
  ZToMMHist->SetMaximum(1.5);
  ZToMMHist->SetMinimum(0.0);
  ZToMMHist->Draw("colz");
  ZToMMHist->SetStats(false);
  ZToMMHist->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetRightMargin(0.15);
  cv->SaveAs("ZToMMScaleFactor.gif");


  TH2D* CompareZToLL_EEVsMM = (TH2D*)ZToMMHist->Clone("CompareZToLL_EEVsMM");
  for (int i=0; i<CompareZToLL_EEVsMM->GetXaxis()->GetNbins()+1;i++) {
    for (int j=0; j<CompareZToLL_EEVsMM->GetYaxis()->GetNbins()+1;j++) {
      CompareZToLL_EEVsMM->SetBinContent(i,j, (ZToMMHist->GetBinContent(i,j) - ZToEEHist->GetBinContent(i,j)) / sqrt( pow(ZToMMHist->GetBinError(i,j), 2) + pow(ZToEEHist->GetBinError(i,j), 2)));
    }
  }
  cv = new TCanvas("cv","cv",800,600);
  CompareZToLL_EEVsMM->SetMaximum(5);
  CompareZToLL_EEVsMM->SetMinimum(-5);
  CompareZToLL_EEVsMM->Draw("");
  CompareZToLL_EEVsMM->SetStats(false);
  CompareZToLL_EEVsMM->GetZaxis()->SetTitle("Difference (#sigma)");
  CompareZToLL_EEVsMM->Draw("colz");
  cv->SetRightMargin(0.15);
  cv->SetLogy();
  cv->SaveAs("CompareZToLL_EEVsMM.gif");

  //********************************************************************************************
  //Combine SingleMu and SingleEle
  //********************************************************************************************
  TH2D* ZToLLDilepton = (TH2D*)ZToMMHist->Clone("ZToLLDilepton");
  TH2D* ZToMMDataMinusBkgHist = (TH2D*)ZToMMFile->Get("DataMinusBkg_MRVsRsq");
  TH2D* ZToEEDataMinusBkgHist = (TH2D*)ZToEEFile->Get("DataMinusBkg_MRVsRsq");
  TH2D* ZToMMMCHist = (TH2D*)ZToMMFile->Get("histMRVsRsq_DY");
  TH2D* ZToEEMCHist = (TH2D*)ZToEEFile->Get("histMRVsRsq_DY");

  for (int i=0; i<CompareZToLL_EEVsMM->GetXaxis()->GetNbins()+1;i++) {
    for (int j=0; j<CompareZToLL_EEVsMM->GetYaxis()->GetNbins()+1;j++) {

      double mc = ZToMMMCHist->GetBinContent(i,j) + ZToEEMCHist->GetBinContent(i,j);
      double mcErr = sqrt( pow(ZToMMMCHist->GetBinError(i,j),2) + pow(ZToEEMCHist->GetBinError(i,j),2));

      double data = ZToMMDataMinusBkgHist->GetBinContent(i,j) + ZToEEDataMinusBkgHist->GetBinContent(i,j);
      double dataErr = sqrt( pow(ZToMMDataMinusBkgHist->GetBinError(i,j),2) + pow(ZToEEDataMinusBkgHist->GetBinError(i,j),2));

      cout << "bin " << i << j << " : " << ZToMMDataMinusBkgHist->GetBinContent(i,j) << " " << ZToEEDataMinusBkgHist->GetBinContent(i,j)  << " | "
  	   << ZToMMMCHist->GetBinContent(i,j) << " " <<  ZToEEMCHist->GetBinContent(i,j) << " "
  	   << "\n";
     // cout << "bin " << i << j << " : " << data << " +/- " << dataErr << " , " << mc << " +/- " << mcErr << "\n";

      ZToLLDilepton->SetBinContent(i,j, data/mc);
      ZToLLDilepton->SetBinError(i,j, (data/mc)*sqrt( pow(mcErr/mc,2)+pow(dataErr/data,2) ) );
    }
  }
  cv = new TCanvas("cv","cv",800,600);
  ZToLLDilepton->SetMaximum(1.5);
  ZToLLDilepton->SetMinimum(0.0);
  ZToLLDilepton->Draw("colz,texte1");
  ZToLLDilepton->SetStats(false);
  ZToLLDilepton->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetLogy();
  cv->SetRightMargin(0.15);
  cv->SaveAs("ZToLLDileptonScaleFactor.gif");
  


  TFile *file = TFile::Open("ZToLLScaleFactors.root", "RECREATE");
  file->cd();
  
  file->WriteTObject(ZToLLDilepton, "ZToLLDileptonScaleFactor", "WriteDelete");
  file->WriteTObject(ZToEEHist, "ZToEEScaleFactor", "WriteDelete");
  file->WriteTObject(ZToMMHist, "ZToMMScaleFactor", "WriteDelete");
  file->Close();
  delete file;


}


void PlotMCToDataScaleFactors_TTBarDilepton( ) {
  
  // string Label = "";
  // if (label != "") Label = "_" + label;


  TFile *TTBarDileptonToEMFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/TTBarDileptonControlRegionPlots_MR300Rsq0p15_emu.root", "READ");
  TFile *TTBarDileptonToEEFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/TTBarDileptonControlRegionPlots_MR300Rsq0p15_ee.root", "READ");
  TFile *TTBarDileptonToMMFile = TFile::Open("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/TTBarDileptonControlRegionPlots_MR300Rsq0p15_mumu.root", "READ");

  TH2D* TTBarDileptonToEMHist = (TH2D*)TTBarDileptonToEMFile->Get("MCToDataScaleFactor_MRVsRsq");
  TH2D* TTBarDileptonToEEHist = (TH2D*)TTBarDileptonToEEFile->Get("MCToDataScaleFactor_MRVsRsq");
  TH2D* TTBarDileptonToMMHist = (TH2D*)TTBarDileptonToMMFile->Get("MCToDataScaleFactor_MRVsRsq");
  
  TCanvas *cv = 0;

  cv = new TCanvas("cv","cv",800,600);
  TTBarDileptonToEMHist->SetMaximum(1.5);
  TTBarDileptonToEMHist->SetMinimum(0.0);
  TTBarDileptonToEMHist->Draw("colz");
  TTBarDileptonToEMHist->SetStats(false);
  TTBarDileptonToEMHist->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetRightMargin(0.15);
  cv->SaveAs("TTBarDileptonToEMHist.gif");


  cv = new TCanvas("cv","cv",800,600);
  TTBarDileptonToEEHist->SetMaximum(1.5);
  TTBarDileptonToEEHist->SetMinimum(0.0);
  TTBarDileptonToEEHist->Draw("colz");
  TTBarDileptonToEEHist->SetStats(false);
  TTBarDileptonToEEHist->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetRightMargin(0.15);
  cv->SaveAs("TTBarDileptonToEEHist.gif");

  cv = new TCanvas("cv","cv",800,600);
  TTBarDileptonToMMHist->SetMaximum(1.5);
  TTBarDileptonToMMHist->SetMinimum(0.0);
  TTBarDileptonToMMHist->Draw("colz");
  TTBarDileptonToMMHist->SetStats(false);
  TTBarDileptonToMMHist->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetRightMargin(0.15);
  cv->SaveAs("TTBarDileptonToMMHist.gif");



  TH2D* CompareTTBarDilepton_EMVsMM = (TH2D*)TTBarDileptonToEMHist->Clone("CompareTTBarDilepton_EMVsMM");
  TH2D* CompareTTBarDilepton_EMVsEE = (TH2D*)TTBarDileptonToEMHist->Clone("CompareTTBarDilepton_EMVsEE");
  TH2D* CompareTTBarDilepton_EEVsMM = (TH2D*)TTBarDileptonToEMHist->Clone("CompareTTBarDilepton_EEVsMM");
  for (int i=0; i<CompareTTBarDilepton_EMVsMM->GetXaxis()->GetNbins()+1;i++) {
    for (int j=0; j<CompareTTBarDilepton_EMVsMM->GetYaxis()->GetNbins()+1;j++) {
      CompareTTBarDilepton_EEVsMM->SetBinContent(i,j, (TTBarDileptonToMMHist->GetBinContent(i,j) - TTBarDileptonToEEHist->GetBinContent(i,j)) / sqrt( pow(TTBarDileptonToMMHist->GetBinError(i,j), 2) + pow(TTBarDileptonToEEHist->GetBinError(i,j), 2)));
      CompareTTBarDilepton_EMVsMM->SetBinContent(i,j, (TTBarDileptonToEMHist->GetBinContent(i,j) - TTBarDileptonToMMHist->GetBinContent(i,j)) / sqrt( pow(TTBarDileptonToEMHist->GetBinError(i,j), 2) + pow(TTBarDileptonToMMHist->GetBinError(i,j), 2)));
      CompareTTBarDilepton_EMVsEE->SetBinContent(i,j, (TTBarDileptonToEMHist->GetBinContent(i,j) - TTBarDileptonToEEHist->GetBinContent(i,j)) / sqrt( pow(TTBarDileptonToEMHist->GetBinError(i,j), 2) + pow(TTBarDileptonToEEHist->GetBinError(i,j), 2)));
      cout << "Bin " << i << " " << j << " : " << TTBarDileptonToMMHist->GetBinContent(i,j) << " +/- " << TTBarDileptonToEMHist->GetBinError(i,j) 
	   << " , " << TTBarDileptonToEEHist->GetBinContent(i,j) << " +/- " << TTBarDileptonToEEHist->GetBinError(i,j) << " | " 
	   << (TTBarDileptonToMMHist->GetBinContent(i,j) - TTBarDileptonToEEHist->GetBinContent(i,j)) / sqrt( pow(TTBarDileptonToMMHist->GetBinError(i,j), 2) + pow(TTBarDileptonToEEHist->GetBinError(i,j), 2)) 
	   << "\n";

    }
  }
  cv = new TCanvas("cv","cv",800,600);
  CompareTTBarDilepton_EMVsMM->SetMaximum(5);
  CompareTTBarDilepton_EMVsMM->SetMinimum(-5);
  CompareTTBarDilepton_EMVsMM->Draw("colz");
  CompareTTBarDilepton_EMVsMM->SetStats(false);
  CompareTTBarDilepton_EMVsMM->GetZaxis()->SetTitle("Difference in #sigma");
  cv->SetRightMargin(0.15);
  cv->SaveAs("CompareTTBarDilepton_EMVsMM.gif");

  cv = new TCanvas("cv","cv",800,600);
  CompareTTBarDilepton_EMVsEE->SetMaximum(5);
  CompareTTBarDilepton_EMVsEE->SetMinimum(-5);
  CompareTTBarDilepton_EMVsEE->Draw("colz");
  CompareTTBarDilepton_EMVsEE->SetStats(false);
  CompareTTBarDilepton_EMVsEE->GetZaxis()->SetTitle("Difference in #sigma");
  cv->SetRightMargin(0.15);
  cv->SaveAs("CompareTTBarDilepton_EMVsEE.gif");

  cv = new TCanvas("cv","cv",800,600);
  CompareTTBarDilepton_EEVsMM->SetMaximum(5);
  CompareTTBarDilepton_EEVsMM->SetMinimum(-5);
  CompareTTBarDilepton_EEVsMM->Draw("colz");
  CompareTTBarDilepton_EEVsMM->SetStats(false);
  CompareTTBarDilepton_EEVsMM->GetZaxis()->SetTitle("Difference in #sigma");
  cv->SetRightMargin(0.15);
  cv->SaveAs("CompareTTBarDilepton_EEVsMM.gif");



  //********************************************************************************************
  //Combine EE,MM,EM
  //********************************************************************************************
  TH2D* TTBarDilepton = (TH2D*)TTBarDileptonToEMHist->Clone("TTBarDilepton");
  TH2D* TTBarDileptonToEEDataMinusBkgHist = (TH2D*)TTBarDileptonToEEFile->Get("DataMinusBkg_MRVsRsq");
  TH2D* TTBarDileptonToMMDataMinusBkgHist = (TH2D*)TTBarDileptonToMMFile->Get("DataMinusBkg_MRVsRsq");
  TH2D* TTBarDileptonToEMDataMinusBkgHist = (TH2D*)TTBarDileptonToEMFile->Get("DataMinusBkg_MRVsRsq");
  TH2D* TTBarDileptonToEEMCHist = (TH2D*)TTBarDileptonToEEFile->Get("histMRVsRsq_TTJets");
  TH2D* TTBarDileptonToMMMCHist = (TH2D*)TTBarDileptonToMMFile->Get("histMRVsRsq_TTJets");
  TH2D* TTBarDileptonToEMMCHist = (TH2D*)TTBarDileptonToEMFile->Get("histMRVsRsq_TTJets");


  for (int i=0; i<CompareTTBarDilepton_EMVsMM->GetXaxis()->GetNbins()+1;i++) {
    for (int j=0; j<CompareTTBarDilepton_EMVsMM->GetYaxis()->GetNbins()+1;j++) {

      double mc = TTBarDileptonToEEMCHist->GetBinContent(i,j) + TTBarDileptonToMMMCHist->GetBinContent(i,j) + TTBarDileptonToEMMCHist->GetBinContent(i,j);
      double mcErr = sqrt( pow(TTBarDileptonToEEMCHist->GetBinError(i,j),2) + pow(TTBarDileptonToMMMCHist->GetBinError(i,j),2)+ pow(TTBarDileptonToEMMCHist->GetBinError(i,j),2) );

      double data = TTBarDileptonToEEDataMinusBkgHist->GetBinContent(i,j) + TTBarDileptonToMMDataMinusBkgHist->GetBinContent(i,j) + TTBarDileptonToEMDataMinusBkgHist->GetBinContent(i,j);
      double dataErr = sqrt( pow(TTBarDileptonToEEDataMinusBkgHist->GetBinError(i,j),2) + pow(TTBarDileptonToMMDataMinusBkgHist->GetBinError(i,j),2) + pow(TTBarDileptonToEMDataMinusBkgHist->GetBinError(i,j),2));

      cout << "bin " << i << j << " : " 
	   << TTBarDileptonToEEDataMinusBkgHist->GetBinContent(i,j) << " +/-" << TTBarDileptonToEEDataMinusBkgHist->GetBinError(i,j) << " " 
	   << TTBarDileptonToMMDataMinusBkgHist->GetBinContent(i,j) << " +/-" << TTBarDileptonToMMDataMinusBkgHist->GetBinError(i,j) << " "
	   << TTBarDileptonToEMDataMinusBkgHist->GetBinContent(i,j)  << " +/- " << TTBarDileptonToEMDataMinusBkgHist->GetBinError(i,j) << " | "
	//<< TTBarDileptonToEEMCHist->GetBinContent(i,j) << " " <<  TTBarDileptonToMMMCHist->GetBinContent(i,j) << " " << TTBarDileptonToEMMCHist->GetBinContent(i,j) << " " 
	   << data << " +/- " << dataErr << " : " << mc << " +/- " << mcErr 
	   << " | " << data/mc << " +/- " << (data/mc)*sqrt( pow(mcErr/mc,2)+pow(dataErr/data,2) )
  	   << "\n";

      TTBarDilepton->SetBinContent(i,j, data/mc);
      TTBarDilepton->SetBinError(i,j, (data/mc)*sqrt( pow(mcErr/mc,2)+pow(dataErr/data,2) ) );
    }
  }
  cv = new TCanvas("cv","cv",800,600);
  TTBarDilepton->SetMaximum(1.5);
  TTBarDilepton->SetMinimum(0.0);
  TTBarDilepton->Draw("colz,texte1");
  TTBarDilepton->SetStats(false);
  TTBarDilepton->GetZaxis()->SetTitle("MC to Data Scale Factor");
  cv->SetRightMargin(0.15);
  cv->SetLogy();
  cv->SaveAs("TTBarDileptonScaleFactor.gif");
  


  TFile *file = TFile::Open("TTBarDileptonScaleFactors.root", "RECREATE");
  file->cd();  
  file->WriteTObject(TTBarDilepton, "TTBarDileptonScaleFactor", "WriteDelete");
  file->WriteTObject(TTBarDileptonToEEHist, "TTBarDileptonToEEScaleFactor", "WriteDelete");
  file->WriteTObject(TTBarDileptonToMMHist, "TTBarDileptonToMMScaleFactor", "WriteDelete");
  file->WriteTObject(TTBarDileptonToEMHist, "TTBarDileptonToEMScaleFactor", "WriteDelete");
  file->Close();
  delete file;


}



void ComputeMCToDataScaleFactors( int option = 0 ) {
    
  if (option == 0) PlotMCToDataScaleFactors_ZToLL();
  if (option == 1) PlotMCToDataScaleFactors_TTBarDilepton();
  if (option == 2) PlotMCToDataScaleFactors_WJets();
  if (option == 3) PlotMCToDataScaleFactors_TTBarSingleLepton();

}
