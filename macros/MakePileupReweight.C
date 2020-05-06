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
#include <TRandom3.h> 
#include <TLatex.h> 

#include "cms_lpc_llp/llp_analyzer/include/ControlSampleEvents.h"
#include "cms_lpc_llp/llp_analyzer/macros/tdrstyle.C"
#include "cms_lpc_llp/llp_analyzer/macros/CMS_lumi.C"

#endif

//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
TH1F* NormalizeHist(TH1F *originalHist) {
  TH1F* hist = (TH1F*)originalHist->Clone((string(originalHist->GetName())+"_normalized").c_str());
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



void MakePileupReweight(int option = 0, int year = 2016, string process = "ZJetsToNuNu_HT-100ToInf_13TeV-madgraph") {


  TFile *pileupTargetFile = 0;
  TFile *pileupSourceFile = 0;
  TFile *file = 0 ;

/*
  //For 2016 Data
   pileupSourceFile = new TFile("RazorAnalyzer/data/PileupWeights/PileupSource_MC80X_Summer16.root", "READ");
   if (option == 0) pileupTargetFile = new TFile("RazorAnalyzer/data/PileupWeights/PileupTarget_2016_36p2ifb.root", "READ");
   else if (option == 1) pileupTargetFile = new TFile("RazorAnalyzer/data/PileupWeights/PileupTarget_2016_36p2ifb_SysUp.root", "READ");
   else if (option == 2) pileupTargetFile = new TFile("RazorAnalyzer/data/PileupWeights/PileupTarget_2016_36p2ifb_SysDown.root", "READ");
   else {
     return;
   }
   file = TFile::Open("PileupReweight_Summer16_2016_36p2ifb.root", "UPDATE");

  //For 2017 Data
  pileupSourceFile = new TFile("RazorAnalyzer/data/PileupWeights/PileupSource_MC_Fall2017.root", "READ");
  // /afs/cern.ch/user/j/jmao/work/public/releases/CMSSW_9_2_1/src/RazorAnalyzer/data/PileupWeights/PileupSource_MC_Fall2017.root
  if (option == 0) pileupTargetFile = new TFile("RazorAnalyzer/data/PileupWeights/PileupTarget_2017Rereco_41p2ifb.root", "READ");
  else if (option == 1) pileupTargetFile = new TFile("RazorAnalyzer/data/PileupWeights/PileupTarget_2017Rereco_41p2ifb_SysUp.root", "READ");
  else if (option == 2) pileupTargetFile = new TFile("RazorAnalyzer/data/PileupWeights/PileupTarget_2017Rereco_41p2ifb_SysDown.root", "READ");
  else {
    return;
  }
  file = TFile::Open("PileupReweight_2017Rereco_41p2ifb.root", "UPDATE");
*/

  //For calo timing analysis
  string campaign = "Summer16";
  if (year==2016) campaign = "Summer16";
  else if (year==2017) campaign = "Fall17";
  else if (year==2018) campaign = "Autumn18";

  string pathname = getenv("CMSSW_BASE");

   pileupSourceFile = new TFile(Form("%s/src/cms_lpc_llp/llp_analyzer/data/PileupWeights/PileupSource_%s_%s_%d_calo.root", pathname.c_str(), process.c_str(), campaign.c_str(), year), "READ");
   if (option == 0) pileupTargetFile = new TFile(Form("%s/src/cms_lpc_llp/llp_analyzer/data/PileupWeights/PileupTarget_%d_calo.root", pathname.c_str(), year), "READ");
   else if (option == 1) pileupTargetFile = new TFile(Form("%s/src/cms_lpc_llp/llp_analyzer/data/PileupWeights/PileupTarget_%d_calo_sysUp.root", pathname.c_str(), year), "READ");
   else if (option == 2) pileupTargetFile = new TFile(Form("%s/src/cms_lpc_llp/llp_analyzer/data/PileupWeights/PileupTarget_%d_calo_sysDown.root", pathname.c_str(), year), "READ");
   else {
     return;
   }
   file = TFile::Open(Form("%s/src/cms_lpc_llp/llp_analyzer/data/PileupWeights/PileupReweight_%s_%s_%d_calo.root", pathname.c_str(), process.c_str(), campaign.c_str(), year), "UPDATE");


  TH1F *pileupTargetHist = (TH1F*)pileupTargetFile->Get("pileup");
  assert(pileupTargetHist);
  std::cout << "pileupTargetHist " << pileupTargetHist->Integral() << std::endl;

  //TH1F *pileupSourceHist = (TH1F*)pileupSourceFile->Get("PileupSourceHist");
  TH1F *pileupSourceHist = (TH1F*)pileupSourceFile->Get("PUMean"); 
  //TH1F *pileupSourceHist = (TH1F*)pileupSourceFile->Get("PUMean_Razor2017_92X"); 
  assert(pileupSourceHist);
  std::cout << "pileupSourceFile " << pileupSourceHist->Integral() << std::endl;
  std::cout << "FILES RETRIEVED" << std::endl;
  //*******************************************************************************************
  //Make NVtx Reweighting Function
  //*******************************************************************************************
  TH1F *PileupTargetNormalized = NormalizeHist( pileupTargetHist );
  TH1F *PileupSourceNormalized = NormalizeHist( pileupSourceHist );
  std::cout << "HISTOS NORMALIZED" << std::endl;
  TH1F *PileupReweight = new TH1F ("PileupReweight",";NPU;Weight",100,0,100);

  for (int i=0; i<PileupReweight->GetXaxis()->GetNbins()+2; i++) {

    double data = 0;
    double bkg = 0;
    if (PileupSourceNormalized->GetBinContent(i) > 0) {
      PileupReweight->SetBinContent(i,PileupTargetNormalized->GetBinContent(i)/PileupSourceNormalized->GetBinContent(i));
    } else if (PileupTargetNormalized->GetBinContent(i) == 0){
      PileupReweight->SetBinContent(i,0.0);
    } else {
      if (i == 1) {
	PileupReweight->SetBinContent(i,1);
      } else {
	PileupReweight->SetBinContent(i,PileupReweight->GetBinContent(i-1));
      }
    }

    cout << "Bin " << i << " : " << PileupReweight->GetXaxis()->GetBinCenter(i) << " : " << PileupReweight->GetBinCenter(i) << " : " << PileupTargetNormalized->GetBinContent(i) << " / " << PileupSourceNormalized->GetBinContent(i) << " = " << PileupReweight->GetBinContent(i) << "\n";
  }

  
  double k = 0;
  for (int i=0; i<PileupSourceNormalized->GetXaxis()->GetNbins()+2; i++) {
    double weight = PileupReweight->GetBinContent(PileupReweight->GetXaxis()->FindFixBin(i-1));
    cout << i << " : " << i-1 << " : " << weight << " : " << PileupSourceNormalized->GetBinContent(i) << " -> " << PileupSourceNormalized->GetBinContent(i)*weight << "\n";
    k += PileupSourceNormalized->GetBinContent(i)*weight;    
  }
  cout << "int = " << k << "\n";

  file->cd();
  if (option == 0) file->WriteTObject(PileupReweight, "PileupReweight", "WriteDelete");
  else if (option == 1) file->WriteTObject(PileupReweight, "PileupReweightSysUp", "WriteDelete");
  else if (option == 2) file->WriteTObject(PileupReweight, "PileupReweightSysDown", "WriteDelete");
  file->Close();
  delete file;

}





