//--------------------------------------------------------------
//
// make ratio plots for Razor sideband and dijet control region
//
//--------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TROOT.h"                        // access to gROOT, entry point to ROOT system
#include "TSystem.h"                      // interface to OS
#include "TStyle.h"                       // class to handle ROOT plotting styles
#include "TFile.h"                        // file handle class
#include "TTree.h"                        // class to access ntuples
#include "TEventList.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2Poly.h"
#include "TBenchmark.h"                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#endif

#include "CalStyleRemix.hh"

void processFile(TString inputFileName, TString outputFileName, int physProc);
void calculateMCScaleFactors() {
  cout << "W+jets" << endl;
  processFile("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_27Nov2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_WJets_1pb_weighted.root", "WJets.root", 0);
  cout << "Z+jets" << endl;
  processFile("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_27Nov2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_ZInv_1pb_weighted.root", "ZInv.root", 1);
  cout << "tt+jets" << endl;
  processFile("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_27Nov2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_TTJets_1pb_weighted.root", "TTJets.root", 2);
}

void processFile(TString inputFileName, TString outputFileName, int physProc) {

  TFile *inputFile = TFile::Open(inputFileName,"read");
  TTree *inputTree = (TTree*) inputFile->Get("RazorInclusive");

  float MR,Rsq,dPhiRazor,leadingJetPt;
  int nSelectedJets,nBTaggedJets, NISRJets;
  float met, metOverCaloMet;
  int box;
  float topPtWeight;
  float weight;

  inputTree->SetBranchAddress("MR",             &MR);
  inputTree->SetBranchAddress("Rsq",            &Rsq);
  inputTree->SetBranchAddress("dPhiRazor",      &dPhiRazor);
  inputTree->SetBranchAddress("met",            &met);
  inputTree->SetBranchAddress("weight",         &weight);
  inputTree->SetBranchAddress("box",            &box);
  inputTree->SetBranchAddress("leadingJetPt",   &leadingJetPt);
  inputTree->SetBranchAddress("nSelectedJets",  &nSelectedJets);
  inputTree->SetBranchAddress("nBTaggedJets",   &nBTaggedJets);
  inputTree->SetBranchAddress("NISRJets",       &NISRJets);
  inputTree->SetBranchAddress("metOverCaloMet", &metOverCaloMet);
  inputTree->SetBranchAddress("topPtWeight",    &topPtWeight);

  inputTree->GetEntry(0);

  TString cut="(box==11||box==12||box==14)*(MR>400 && Rsq>0.1)*(Flag_HBHENoiseFilter && Flag_HBHEIsoNoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter && Flag_badChargedCandidateFilter && Flag_badMuonFilter)";

  TFile *fScaleFactors = TFile::Open("../../../data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Razor2016_MoriondRereco.root");
  TFile *fNJetScaleFactors = TFile::Open("../../../data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_Razor2016_MoriondRereco.root");
  TFile *fBTagScaleFactors = TFile::Open("../../../data/ScaleFactors/RazorMADD2015/RazorBTagScaleFactors_Razor2016_MoriondRereco.root");

  TH2Poly *hScaleFactors=0; 
  TH1F *hNJetScaleFactors=0; 
  vector<TH1F*> hBTagScaleFactorsDiJet;
  vector<TH1F*> hBTagScaleFactorsMultiJet;
  vector<TH1F*> hBTagScaleFactorsSevenJet;

  if (physProc==0) {
    hScaleFactors= (TH2Poly*) fScaleFactors->Get("WJetsScaleFactors");
    hNJetScaleFactors= (TH1F*) fNJetScaleFactors->Get("WJetsScaleFactors");
  }
  else if (physProc==1) {
    hScaleFactors= (TH2Poly*) fScaleFactors->Get("GJetsInvScaleFactors");
    hNJetScaleFactors= (TH1F*) fNJetScaleFactors->Get("GJetsInvScaleFactors");
    for (int nb = 0; nb <= 3; nb++) {
        hBTagScaleFactorsDiJet.push_back(
                (TH1F*)fBTagScaleFactors->Get(Form(
                        "MRInvDiJet%dBScaleFactors", min(nb, 2))));
        hBTagScaleFactorsMultiJet.push_back(
                (TH1F*)fBTagScaleFactors->Get(Form(
                        "MRInvMultiJet%dBScaleFactors", min(nb, 2))));
        hBTagScaleFactorsSevenJet.push_back(
                (TH1F*)fBTagScaleFactors->Get(Form(
                        "MRInvSevenJet%dBScaleFactors", min(nb, 2))));
    }
  }
  else if (physProc==2) {
    hScaleFactors= (TH2Poly*) fScaleFactors->Get("TTJetsScaleFactors");
    hNJetScaleFactors= (TH1F*) fNJetScaleFactors->Get("TTJetsScaleFactors");
  }

  if (physProc != 1) {
    for (int nb = 0; nb <= 3; nb++) {
        hBTagScaleFactorsDiJet.push_back(
                (TH1F*)fBTagScaleFactors->Get(Form(
                        "MRDiJet%dBScaleFactors", min(nb, 2))));
        hBTagScaleFactorsMultiJet.push_back(
                (TH1F*)fBTagScaleFactors->Get(Form(
                        "MRMultiJet%dBScaleFactors", nb)));
        hBTagScaleFactorsSevenJet.push_back(
                (TH1F*)fBTagScaleFactors->Get(Form(
                        "MRSevenJet%dBScaleFactors", nb)));
    }
  }

  TFile *outputFile = new TFile(outputFileName,"recreate");
  TTree *outputTree = (TTree*) inputTree->GetTree()->CloneTree(0);

  float mcScaleFactor;
  outputTree->Branch("mcScaleFactor",&mcScaleFactor, "mcScaleFactor/F");

  inputTree->Draw(">>elist1",cut);
  
  TEventList *list = (TEventList*)gDirectory->Get("elist1");

  for (Int_t i=0; i<list->GetN(); i++) {
    inputTree->GetEntry(list->GetEntry(i));

    double tNJets=min((double)nSelectedJets, hNJetScaleFactors->GetXaxis()->GetXmax()*0.999);
    tNJets=max(tNJets, hNJetScaleFactors->GetXaxis()->GetXmin()*1.001);

    double tMR=min((double)MR, hScaleFactors->GetXaxis()->GetXmax()*0.999);
    tMR=max(tMR, hScaleFactors->GetXaxis()->GetXmin()*1.001);

    double tRsq=min((double)Rsq, hScaleFactors->GetYaxis()->GetXmax()*0.999);
    tRsq=max(tRsq, hScaleFactors->GetYaxis()->GetXmin()*1.001);

    int tNBtags=min(nBTaggedJets, 3);

    double scaleFactor = hScaleFactors->GetBinContent(hScaleFactors->FindBin(tMR, tRsq));
    double njetScaleFactor = hNJetScaleFactors->GetBinContent(hNJetScaleFactors->FindFixBin(tNJets));

    TH1F *btagHist = 0;
    if (tNJets >= 7) {
        btagHist = hBTagScaleFactorsSevenJet[tNBtags];
    }
    else if (tNJets >= 4) {
        btagHist = hBTagScaleFactorsMultiJet[tNBtags];
    }
    else {
        btagHist = hBTagScaleFactorsDiJet[tNBtags];
    }
    double btagScaleFactor = btagHist->GetBinContent(
            btagHist->FindFixBin(tMR));

    if (physProc==2) {
        scaleFactor *= topPtWeight;
    }

    if (scaleFactor==0 || scaleFactor>2) {
      cout << scaleFactor << endl;
      scaleFactor=1;
    }
    if (njetScaleFactor==0 || njetScaleFactor>2) {
      cout << njetScaleFactor << endl;
      njetScaleFactor=1;
    }

    mcScaleFactor=scaleFactor*njetScaleFactor*btagScaleFactor;

    outputTree->Fill();
  }

  outputFile->Write();
  outputFile->Close();

}
