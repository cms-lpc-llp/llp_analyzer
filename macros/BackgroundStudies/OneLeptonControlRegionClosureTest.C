
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <vector>
#include <map>
#include <iostream>
#include "RazorAnalyzer/macros/tdrstyle.C"
#include "RazorAnalyzer/macros/CMS_lumi.C"

const Int_t NComponents = 10;
int color[NComponents] = {kRed, kGreen+2, kBlue, kViolet, kAzure+10, kBlack, kOrange+1, kGray, kBlack, kBlack};


//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
TH1D* NormalizeHist(TH1D *originalHist) {
  TH1D* hist = (TH1D*)originalHist->Clone((string(originalHist->GetName())+"_normalized").c_str());
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=1; int(b)<hist->GetXaxis()->GetNbins()+1; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=1; int(b)<hist->GetXaxis()->GetNbins()+1; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return hist;
}




//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void OneLeptonControlRegionClosureTest( ) {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
   
  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  TH1D* histPtW = 0;
  histPtW =  new TH1D( "histPtW",";Recoil [GeV/c];Number of Events", 60, 100, 700);
  histPtW->SetStats(false);
  histPtW->Sumw2();
  TH1D* histPtW_Corrected = 0;
  histPtW_Corrected =  new TH1D( "histPtW_Corrected","; Recoil [GeV/c];Number of Events", 60, 100, 700);
  histPtW_Corrected->SetStats(false);
  histPtW_Corrected->Sumw2();


  TH1D* histMR_lowMT_Uncorrected = 0;
  histMR_lowMT_Uncorrected =  new TH1D( "histMR_lowMT_Uncorrected",";MR [GeV/c^{2}];Number of Events", 40, 400, 2400);
  histMR_lowMT_Uncorrected->SetStats(false);
  histMR_lowMT_Uncorrected->Sumw2();

  TH1D* histMR_lowMT_Corrected = 0;
  histMR_lowMT_Corrected =  new TH1D( "histMR_lowMT_Corrected",";MR [GeV/c^{2}];Number of Events", 40, 400, 2400);
  histMR_lowMT_Corrected->SetStats(false);
  histMR_lowMT_Corrected->Sumw2();

  TH1D* histMR_highMT_Uncorrected = 0;
  histMR_highMT_Uncorrected =  new TH1D( "histMR_highMT_Uncorrected",";MR [GeV/c^{2}];Number of Events", 40, 400, 2400);
  histMR_highMT_Uncorrected->SetStats(false);
  histMR_highMT_Uncorrected->Sumw2();

  TH1D* histMR_highMT_Corrected = 0;
  histMR_highMT_Corrected =  new TH1D( "histMR_highMT_Corrected",";MR [GeV/c^{2}];Number of Events", 40, 400, 2400);
  histMR_highMT_Corrected->SetStats(false);
  histMR_highMT_Corrected->Sumw2();

  // TH1D* histMR_highMT_Prediction = 0;
  // histMR_highMT_Prediction =  new TH1D( "histMR_highMT_Prediction",";MR [GeV/c^{2}];Number of Events", 60, 0, 3000);
  // histMR_highMT_Prediction->SetStats(false);
  // histMR_highMT_Prediction->Sumw2();


  //*******************************************************************************************
  //Load Input File
  //*******************************************************************************************  
  // TFile* inputFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForMoriond20160119/RazorInclusive_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_25ns.root","READ");
  TFile* inputFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForMoriond20160119/RazorInclusive_WJetsToLNu_HTBinned_1pb_weighted.root","READ");
  assert(inputFile);
  TTree* tree = 0;
  tree = (TTree*)inputFile->Get("RazorInclusive");
  assert(tree);    
  
  float weight = 0;
  int box = -1;
  int nBTaggedJets = 0;
  float dPhiRazor = 0;
  float MR = 0;
  float Rsq = 0;
  float mT = 0;
  float mTLoose = 0;
  float ptW = 0;

  tree->SetBranchAddress("weight",&weight);    
  tree->SetBranchAddress("mT",&mT);    
  tree->SetBranchAddress("MR",&MR);    
  tree->SetBranchAddress("ptW",&ptW);    
  
  cout << "Process : Total Events: " << tree->GetEntries() << "\n";
  for (int n=0;n<tree->GetEntries();n++) { 
    
    tree->GetEntry(n);
    if (n % 1000000 == 0) cout << "Processing Event " << n << "\n";       
    if (ptW <= 100) continue;

    double w = weight * 2200;
    double correction = fmax(1.0 - 0.3*(ptW / 600),0);
    histPtW->Fill(ptW,w);
    if (mT < 100) {
      histMR_lowMT_Uncorrected->Fill(MR,w);
      histMR_lowMT_Corrected->Fill(MR,w*correction);
    }
    if (mT > 120) {
      histMR_highMT_Uncorrected->Fill(MR,w);
      histMR_highMT_Corrected->Fill(MR,w*correction);
    }
    histPtW_Corrected->Fill(ptW, weight*correction);    
  }

	
  //Normalize Histograms
  histPtW = NormalizeHist(histPtW);
  histPtW_Corrected = NormalizeHist(histPtW_Corrected);
  // histMR_lowMT_Uncorrected = NormalizeHist(histMR_lowMT_Uncorrected);
  // histMR_lowMT_Corrected = NormalizeHist(histMR_lowMT_Corrected);
  // histMR_highMT_Uncorrected = NormalizeHist(histMR_highMT_Uncorrected);
  // histMR_highMT_Corrected = NormalizeHist(histMR_highMT_Corrected);
  
  //Derive Correction from low MT region
  TH1D* recoilCorrection  = (TH1D*)histMR_lowMT_Corrected->Clone("recoilCorrection");
  for (UInt_t b=1; int(b)<recoilCorrection->GetXaxis()->GetNbins()+1; ++b) {
    cout << "bin " << b << " : " << histMR_lowMT_Corrected->GetBinContent(b) / histMR_lowMT_Uncorrected->GetBinContent(b) << " | " 
	 << histMR_highMT_Corrected->GetBinContent(b) / histMR_highMT_Uncorrected->GetBinContent(b) << " | "
	 << histMR_highMT_Uncorrected->GetBinContent(b) << " "
	 << "\n";
    recoilCorrection->SetBinContent(b,histMR_lowMT_Corrected->GetBinContent(b) / histMR_lowMT_Uncorrected->GetBinContent(b) );
  }
  

  //Apply low MT correction to high MT region
  TH1D* histMR_highMT_Prediction  = (TH1D*)histMR_highMT_Corrected->Clone("histMR_highMT_Prediction");
  TH1D* histMR_highMT_PredictionOverCorrected  = (TH1D*)histMR_highMT_Corrected->Clone("histMR_highMT_PredictionOverCorrected");
  for (UInt_t b=1; int(b)<histMR_highMT_Prediction->GetXaxis()->GetNbins()+1; ++b) {
    cout << b << " : " << histMR_highMT_Uncorrected->GetBinContent(b) << " " << recoilCorrection->GetBinContent(b) << " --> "
	 << histMR_highMT_Uncorrected->GetBinContent(b) * recoilCorrection->GetBinContent(b) << " | "
	 << histMR_highMT_Corrected->GetBinContent(b) << " " 
	 << "\n";
    histMR_highMT_Prediction->SetBinContent(b,histMR_highMT_Uncorrected->GetBinContent(b) * recoilCorrection->GetBinContent(b) );
    histMR_highMT_Prediction->SetBinError(b,histMR_highMT_Uncorrected->GetBinError(b) * recoilCorrection->GetBinContent(b) );
  }

  // inputFile->Close();
  // delete inputFile;
  
  // histMR_highMT_Prediction = NormalizeHist(histMR_highMT_Prediction);
  // histMR_lowMT_Uncorrected = NormalizeHist(histMR_lowMT_Uncorrected);
  // histMR_lowMT_Corrected = NormalizeHist(histMR_lowMT_Corrected);
  // histMR_highMT_Uncorrected = NormalizeHist(histMR_highMT_Uncorrected);
  // histMR_highMT_Corrected = NormalizeHist(histMR_highMT_Corrected);
  

  for (UInt_t b=1; int(b)<histMR_highMT_PredictionOverCorrected->GetXaxis()->GetNbins()+1; ++b) {
    histMR_highMT_PredictionOverCorrected->SetBinContent(b,histMR_highMT_Prediction->GetBinContent(b) / histMR_highMT_Corrected->GetBinContent(b));
    histMR_highMT_PredictionOverCorrected->SetBinError(b, sqrt( pow(histMR_highMT_Prediction->GetBinError(b)/histMR_highMT_Prediction->GetBinContent(b),2) + pow(histMR_highMT_Corrected->GetBinError(b) / histMR_highMT_Corrected->GetBinContent(b),2)) );    
  }

   //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open("OneLeptonCRClosureTest.root", "UPDATE");
  file->cd();

  file->WriteTObject(histPtW, "histPtW", "WriteDelete");
  file->WriteTObject(histPtW_Corrected, "histPtW_Corrected", "WriteDelete");
  file->WriteTObject(histMR_lowMT_Uncorrected, "histMR_lowMT_Uncorrected", "WriteDelete");
  file->WriteTObject(histMR_lowMT_Corrected, "histMR_lowMT_Corrected", "WriteDelete");
  file->WriteTObject(histMR_highMT_Uncorrected, "histMR_highMT_Uncorrected", "WriteDelete");
  file->WriteTObject(histMR_highMT_Corrected, "histMR_highMT_Corrected", "WriteDelete");
  file->WriteTObject(histMR_highMT_Prediction, "histMR_highMT_Prediction", "WriteDelete");
  file->WriteTObject(histMR_highMT_PredictionOverCorrected, "histMR_highMT_PredictionOverCorrected", "WriteDelete");

 }


