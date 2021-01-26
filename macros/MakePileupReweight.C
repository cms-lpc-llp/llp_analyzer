//================================================================================================
//
// Simple Example
//
//root -l macros/MakePileupReweight.C(0,2016)
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

#include "llp_analyzer/include/ControlSampleEvents.h"
#include "llp_analyzer/macros/tdrstyle.C"
#include "llp_analyzer/macros/CMS_lumi.C"

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



void MakePileupReweight(int option = 0, int year = 2018) {


  TFile *pileupTargetFile = 0;
  TFile *pileupSourceFile = 0;
  TFile *file = 0 ;
  char* cmsswPath;
  string pathname = getenv("CMSSW_BASE");


  if (year == 2018)
  {
    pileupSourceFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupSource_MC_Fall18_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root", pathname.c_str()), "READ");
    if (option == 0) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2018.root", pathname.c_str()), "READ");
    else if (option == 1) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2018_SysUp.root", pathname.c_str()), "READ");
    else if (option == 2) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2018_SysDown.root", pathname.c_str()), "READ");
    else {
      return;
    }
    file = TFile::Open(Form("%s/src/llp_analyzer/data/PileupWeights/PileupReweight_MC_Fall18_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root", pathname.c_str()), "UPDATE");
  }
  else if (year == 2017)
  {
    // pileupSourceFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupSource_MC_Fall17_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root", pathname.c_str()), "READ");
    pileupSourceFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupSource_MC_Fall17_ZToMuMu_NNPDF31_13TeV-powheg_M_50_120.root", pathname.c_str()), "READ");

    if (option == 0) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2017EOYReReco_41p5ifb.root", pathname.c_str()), "READ");
    else if (option == 1) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2017EOYReReco_41p5ifb_SysUp.root", pathname.c_str()), "READ");
    else if (option == 2) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2017EOYReReco_41p5ifb_SysDown.root", pathname.c_str()), "READ");
    else {
      return;
    }
    // file = TFile::Open(Form("%s/src/llp_analyzer/data/PileupWeights/PileupReweight_MC_Fall17_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root", pathname.c_str()), "UPDATE");
    file = TFile::Open(Form("%s/src/llp_analyzer/data/PileupWeights/PileupReweight_MC_Fall17_ZToMuMu_NNPDF31_13TeV-powheg_M_50_120.root", pathname.c_str()), "UPDATE");
  }
  else if (year == 2016)
  {
    pileupSourceFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupSource_MC_Summer16_ggH_HToSSTobbbb_MH-125_TuneCUETP8M1_13TeV-powheg-pythia8.root", pathname.c_str()), "READ");
    if (option == 0) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2016Rereco_07Aug17.root", pathname.c_str()), "READ");
    else if (option == 1) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2016Rereco_07Aug17_SysUp.root", pathname.c_str()), "READ");
    else if (option == 2) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2016Rereco_07Aug17_SysDown.root", pathname.c_str()), "READ");
    else {
      return;
    }
    file = TFile::Open(Form("%s/src/llp_analyzer/data/PileupWeights/PileupReweight_MC_Summer16_ggH_HToSSTobbbb_MH-125_TuneCUETP8M1_13TeV-powheg-pythia8.root", pathname.c_str()), "UPDATE");
  }
  else if (year == 2019)//source is 2018, target is 2016
  {
    pileupSourceFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupSource_MC_Fall18_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root", pathname.c_str()), "READ");
    if (option == 0) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2016Rereco_07Aug17.root", pathname.c_str()), "READ");
    else if (option == 1) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2016Rereco_07Aug17_SysUp.root", pathname.c_str()), "READ");
    else if (option == 2) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2016Rereco_07Aug17_SysDown.root", pathname.c_str()), "READ");
    else {
      return;
    }
    file = TFile::Open(Form("%s/src/llp_analyzer/data/PileupWeights/PileupReweight_source18_target16.root", pathname.c_str()), "UPDATE");
  }
  else if (year == 2020)//source is 2018, target is 2017
  {
    pileupSourceFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupSource_MC_Fall18_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root", pathname.c_str()), "READ");
    if (option == 0) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2017EOYReReco_41p5ifb.root", pathname.c_str()), "READ");
    else if (option == 1) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2017EOYReReco_41p5ifb_SysUp.root", pathname.c_str()), "READ");
    else if (option == 2) pileupTargetFile = new TFile(Form("%s/src/llp_analyzer/data/PileupWeights/PileupTarget_2017EOYReReco_41p5ifb_SysDown.root", pathname.c_str()), "READ");
    else {
      return;
    }
    file = TFile::Open(Form("%s/src/llp_analyzer/data/PileupWeights/PileupReweight_source18_target17.root", pathname.c_str()), "UPDATE");
  }
  else{
    cout<<"[ERROR]: YEAR IS NOT CORRECT"<<endl;
  }

  TH1F *pileupTargetHist = (TH1F*)pileupTargetFile->Get("pileup");
  assert(pileupTargetHist);
  std::cout << "pileupTargetHist " << pileupTargetHist->Integral() << std::endl;

  TH1F *pileupSourceHist = (TH1F*)pileupSourceFile->Get("PUMean");
  assert(pileupSourceHist);
  std::cout << "pileupSourceFile " << pileupSourceHist->Integral() << std::endl;
  std::cout << "FILES RETRIEVED" << std::endl;
  //*******************************************************************************************
  //Make NVtx Reweighting Function
  //*******************************************************************************************
  TH1F *PileupTargetNormalized = NormalizeHist( pileupTargetHist );
  TH1F *PileupSourceNormalized = NormalizeHist( pileupSourceHist );
  std::cout << "HISTOS NORMALIZED" << std::endl;
  TH1F *PileupReweight = new TH1F ("PileupReweight",";NPU;Weight",200,0,200);

  for (int i=0; i<PileupReweight->GetXaxis()->GetNbins()+2; i++) {

    double data = 0;
    double bkg = 0;
    if (PileupSourceNormalized->GetBinContent(i) > 0) {
      PileupReweight->SetBinContent(i,PileupTargetNormalized->GetBinContent(i)/PileupSourceNormalized->GetBinContent(i));
    }
    else if (PileupTargetNormalized->GetBinContent(i) == 0){ //if both 0
      PileupReweight->SetBinContent(i,0.0);
    }
    else {//source is 0, target nonzero
      if (i == 1) {//Sys up
	       PileupReweight->SetBinContent(i,1);
      }
      else {
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
