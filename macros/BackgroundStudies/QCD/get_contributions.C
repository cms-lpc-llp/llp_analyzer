//--------------------------------------------------------------
//
// make plots for Razor-sidebands
//
//--------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TBenchmark.h>                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#endif

#include <CalStyleRemix.hh>

void get_contributions() {

  //--------------------------------------------------------------
  //
  // setup
  //
  //--------------------------------------------------------------

  TCanvas *c = MakeCanvas("c","c",800,600);

  TFile *fD = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p13_05Mar2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root","read");
  TFile *fQ = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_03Feb2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_QCD_1pb_weighted.root","read");
  //TFile *fT = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_03Feb2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_TTJetsInclusive_1pb_weighted.root","read");
  TFile *fT = TFile::Open("TTJets.root","read");
  //TFile *fW = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_03Feb2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_WJets_1pb_weighted.root","read");
  TFile *fW = TFile::Open("WJets.root","read");
  //TFile *fZ = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_03Feb2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_ZInv_1pb_weighted.root","read")
  TFile *fZ = TFile::Open("ZInv.root","read");

  TString cut_str="(box==11||box==12)*(MR>400 && Rsq>0.15)*(Flag_HBHENoiseFilter && Flag_HBHEIsoNoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter && Flag_badChargedCandidateFilter && Flag_badMuonFilter)*(nBTaggedJets==1)";
  TString cut_str_rsq="(box==11||box==12)*(MR>400 && Rsq>0.15)*(Flag_HBHENoiseFilter && Flag_HBHEIsoNoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter && Flag_badChargedCandidateFilter && Flag_badMuonFilter)*(Rsq<0.25)*(nBTaggedJets==1)";

  //TString cut_str_mc="*weight*35800*mcScaleFactor";
  TString cut_str_mc="*weight*35800";
  TString cut_str_mc_qcd="*weight*35800";
 
  TString cut_str_mr1="*(MR>500 && MR<600)";
  TString cut_str_mr2="*(MR>600 && MR<700)";
  TString cut_str_mr3="*(MR>700 && MR<800)";
  TString cut_str_mr4="*(MR>800 && MR<900)";
  TString cut_str_mr5="*(MR>900 && MR<1000)";

  TString cut_str_ljpt1="*(leadingJetPt>80  && leadingJetPt<200)";
  TString cut_str_ljpt2="*(leadingJetPt>200 && leadingJetPt<250)";
  TString cut_str_ljpt3="*(leadingJetPt>250 && leadingJetPt<300)";
  TString cut_str_ljpt4="*(leadingJetPt>300 && leadingJetPt<400)";
  TString cut_str_ljpt5="*(leadingJetPt>400 && leadingJetPt<500)";

  TString cut_str_rsq1="*(Rsq>0.15 && Rsq<0.20)";
  TString cut_str_rsq2="*(Rsq>0.20 && Rsq<0.25)";
  TString cut_str_rsq3="*(Rsq>0.25 && Rsq<0.30)";

  // draw MR
  Float_t nbin=5, xmin=500, xmax=1000;

  TH1F *hMR_D = new TH1F("hMR_D", "hMR_D", nbin, xmin, xmax); hMR_D->Sumw2();
  TH1F *hMR_Q = new TH1F("hMR_Q", "hMR_Q", nbin, xmin, xmax); hMR_Q->Sumw2();
  TH1F *hMR_T = new TH1F("hMR_T", "hMR_T", nbin, xmin, xmax); hMR_T->Sumw2();
  TH1F *hMR_W = new TH1F("hMR_W", "hMR_W", nbin, xmin, xmax); hMR_W->Sumw2();
  TH1F *hMR_Z = new TH1F("hMR_Z", "hMR_Z", nbin, xmin, xmax); hMR_Z->Sumw2();

  TH1F *hMR_pass_D = new TH1F("hMR_pass_D", "hMR_pass_D", nbin, xmin, xmax); hMR_pass_D->Sumw2();
  TH1F *hMR_pass_Q = new TH1F("hMR_pass_Q", "hMR_pass_Q", nbin, xmin, xmax); hMR_pass_Q->Sumw2();
  TH1F *hMR_pass_T = new TH1F("hMR_pass_T", "hMR_pass_T", nbin, xmin, xmax); hMR_pass_T->Sumw2();
  TH1F *hMR_pass_W = new TH1F("hMR_pass_W", "hMR_pass_W", nbin, xmin, xmax); hMR_pass_W->Sumw2();
  TH1F *hMR_pass_Z = new TH1F("hMR_pass_Z", "hMR_pass_Z", nbin, xmin, xmax); hMR_pass_Z->Sumw2();

  TH1F *hMR_fail_D = new TH1F("hMR_fail_D", "hMR_fail_D", nbin, xmin, xmax); hMR_fail_D->Sumw2();
  TH1F *hMR_fail_Q = new TH1F("hMR_fail_Q", "hMR_fail_Q", nbin, xmin, xmax); hMR_fail_Q->Sumw2();
  TH1F *hMR_fail_T = new TH1F("hMR_fail_T", "hMR_fail_T", nbin, xmin, xmax); hMR_fail_T->Sumw2();
  TH1F *hMR_fail_W = new TH1F("hMR_fail_W", "hMR_fail_W", nbin, xmin, xmax); hMR_fail_W->Sumw2();
  TH1F *hMR_fail_Z = new TH1F("hMR_fail_Z", "hMR_fail_Z", nbin, xmin, xmax); hMR_fail_Z->Sumw2();

  // draw Rsq
  nbin=8; xmin=0.15; xmax=0.35;

  TH1F *hRsq_D = new TH1F("hRsq_D", "hRsq_D", nbin, xmin, xmax); hRsq_D->Sumw2();
  TH1F *hRsq_Q = new TH1F("hRsq_Q", "hRsq_Q", nbin, xmin, xmax); hRsq_Q->Sumw2();
  TH1F *hRsq_T = new TH1F("hRsq_T", "hRsq_T", nbin, xmin, xmax); hRsq_T->Sumw2();
  TH1F *hRsq_W = new TH1F("hRsq_W", "hRsq_W", nbin, xmin, xmax); hRsq_W->Sumw2();
  TH1F *hRsq_Z = new TH1F("hRsq_Z", "hRsq_Z", nbin, xmin, xmax); hRsq_Z->Sumw2();

  TH1F *hRsq_pass_D = new TH1F("hRsq_pass_D", "hRsq_pass_D", nbin, xmin, xmax); hRsq_pass_D->Sumw2();
  TH1F *hRsq_pass_Q = new TH1F("hRsq_pass_Q", "hRsq_pass_Q", nbin, xmin, xmax); hRsq_pass_Q->Sumw2();
  TH1F *hRsq_pass_T = new TH1F("hRsq_pass_T", "hRsq_pass_T", nbin, xmin, xmax); hRsq_pass_T->Sumw2();
  TH1F *hRsq_pass_W = new TH1F("hRsq_pass_W", "hRsq_pass_W", nbin, xmin, xmax); hRsq_pass_W->Sumw2();
  TH1F *hRsq_pass_Z = new TH1F("hRsq_pass_Z", "hRsq_pass_Z", nbin, xmin, xmax); hRsq_pass_Z->Sumw2();

  TH1F *hRsq_fail_D = new TH1F("hRsq_fail_D", "hRsq_fail_D", nbin, xmin, xmax); hRsq_fail_D->Sumw2();
  TH1F *hRsq_fail_Q = new TH1F("hRsq_fail_Q", "hRsq_fail_Q", nbin, xmin, xmax); hRsq_fail_Q->Sumw2();
  TH1F *hRsq_fail_T = new TH1F("hRsq_fail_T", "hRsq_fail_T", nbin, xmin, xmax); hRsq_fail_T->Sumw2();
  TH1F *hRsq_fail_W = new TH1F("hRsq_fail_W", "hRsq_fail_W", nbin, xmin, xmax); hRsq_fail_W->Sumw2();
  TH1F *hRsq_fail_Z = new TH1F("hRsq_fail_Z", "hRsq_fail_Z", nbin, xmin, xmax); hRsq_fail_Z->Sumw2();

  // draw leading jet pt                                                                                                
  Float_t xbins_ljpt[10] = {80, 100, 150, 200, 250, 300, 350, 400, 450, 500};
  nbin=9;

  TH1F *hLjpt_D = new TH1F("hLjpt_D", "hLjpt_D", nbin, &xbins_ljpt[0]); hLjpt_D->Sumw2();
  TH1F *hLjpt_Q = new TH1F("hLjpt_Q", "hLjpt_Q", nbin, &xbins_ljpt[0]); hLjpt_Q->Sumw2();
  TH1F *hLjpt_T = new TH1F("hLjpt_T", "hLjpt_T", nbin, &xbins_ljpt[0]); hLjpt_T->Sumw2();
  TH1F *hLjpt_W = new TH1F("hLjpt_W", "hLjpt_W", nbin, &xbins_ljpt[0]); hLjpt_W->Sumw2();
  TH1F *hLjpt_Z = new TH1F("hLjpt_Z", "hLjpt_Z", nbin, &xbins_ljpt[0]); hLjpt_Z->Sumw2();

  TH1F *hLjpt_pass_D = new TH1F("hLjpt_pass_D", "hLjpt_pass_D", nbin, &xbins_ljpt[0]); hLjpt_pass_D->Sumw2();
  TH1F *hLjpt_pass_Q = new TH1F("hLjpt_pass_Q", "hLjpt_pass_Q", nbin, &xbins_ljpt[0]); hLjpt_pass_Q->Sumw2();
  TH1F *hLjpt_pass_T = new TH1F("hLjpt_pass_T", "hLjpt_pass_T", nbin, &xbins_ljpt[0]); hLjpt_pass_T->Sumw2();
  TH1F *hLjpt_pass_W = new TH1F("hLjpt_pass_W", "hLjpt_pass_W", nbin, &xbins_ljpt[0]); hLjpt_pass_W->Sumw2();
  TH1F *hLjpt_pass_Z = new TH1F("hLjpt_pass_Z", "hLjpt_pass_Z", nbin, &xbins_ljpt[0]); hLjpt_pass_Z->Sumw2();

  TH1F *hLjpt_fail_D = new TH1F("hLjpt_fail_D", "hLjpt_fail_D", nbin, &xbins_ljpt[0]); hLjpt_fail_D->Sumw2();
  TH1F *hLjpt_fail_Q = new TH1F("hLjpt_fail_Q", "hLjpt_fail_Q", nbin, &xbins_ljpt[0]); hLjpt_fail_Q->Sumw2();
  TH1F *hLjpt_fail_T = new TH1F("hLjpt_fail_T", "hLjpt_fail_T", nbin, &xbins_ljpt[0]); hLjpt_fail_T->Sumw2();
  TH1F *hLjpt_fail_W = new TH1F("hLjpt_fail_W", "hLjpt_fail_W", nbin, &xbins_ljpt[0]); hLjpt_fail_W->Sumw2();
  TH1F *hLjpt_fail_Z = new TH1F("hLjpt_fail_Z", "hLjpt_fail_Z", nbin, &xbins_ljpt[0]); hLjpt_fail_Z->Sumw2();

  // draw dPhiRazor
  nbin=8; 
  Float_t xbins[9] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 2.8, 3.0, 3.2};

  TH1F *hDPhiR_D = new TH1F("hDPhiR_D", "hDPhiR_D", nbin, &xbins[0]); hDPhiR_D->Sumw2();
  TH1F *hDPhiR_Q = new TH1F("hDPhiR_Q", "hDPhiR_Q", nbin, &xbins[0]); hDPhiR_Q->Sumw2();
  TH1F *hDPhiR_T = new TH1F("hDPhiR_T", "hDPhiR_T", nbin, &xbins[0]); hDPhiR_T->Sumw2();
  TH1F *hDPhiR_W = new TH1F("hDPhiR_W", "hDPhiR_W", nbin, &xbins[0]); hDPhiR_W->Sumw2();
  TH1F *hDPhiR_Z = new TH1F("hDPhiR_Z", "hDPhiR_Z", nbin, &xbins[0]); hDPhiR_Z->Sumw2();

  TH1F *hDPhiR_rsq1_D = new TH1F("hDPhiR_rsq1_D", "hDPhiR_rsq1_D", nbin, &xbins[0]); hDPhiR_rsq1_D->Sumw2();
  TH1F *hDPhiR_rsq1_Q = new TH1F("hDPhiR_rsq1_Q", "hDPhiR_rsq1_Q", nbin, &xbins[0]); hDPhiR_rsq1_Q->Sumw2();
  TH1F *hDPhiR_rsq1_T = new TH1F("hDPhiR_rsq1_T", "hDPhiR_rsq1_T", nbin, &xbins[0]); hDPhiR_rsq1_T->Sumw2();
  TH1F *hDPhiR_rsq1_W = new TH1F("hDPhiR_rsq1_W", "hDPhiR_rsq1_W", nbin, &xbins[0]); hDPhiR_rsq1_W->Sumw2();
  TH1F *hDPhiR_rsq1_Z = new TH1F("hDPhiR_rsq1_Z", "hDPhiR_rsq1_Z", nbin, &xbins[0]); hDPhiR_rsq1_Z->Sumw2();

  TH1F *hDPhiR_rsq2_D = new TH1F("hDPhiR_rsq2_D", "hDPhiR_rsq2_D", nbin, &xbins[0]); hDPhiR_rsq2_D->Sumw2();
  TH1F *hDPhiR_rsq2_Q = new TH1F("hDPhiR_rsq2_Q", "hDPhiR_rsq2_Q", nbin, &xbins[0]); hDPhiR_rsq2_Q->Sumw2();
  TH1F *hDPhiR_rsq2_T = new TH1F("hDPhiR_rsq2_T", "hDPhiR_rsq2_T", nbin, &xbins[0]); hDPhiR_rsq2_T->Sumw2();
  TH1F *hDPhiR_rsq2_W = new TH1F("hDPhiR_rsq2_W", "hDPhiR_rsq2_W", nbin, &xbins[0]); hDPhiR_rsq2_W->Sumw2();
  TH1F *hDPhiR_rsq2_Z = new TH1F("hDPhiR_rsq2_Z", "hDPhiR_rsq2_Z", nbin, &xbins[0]); hDPhiR_rsq2_Z->Sumw2();
  
  TH1F *hDPhiR_rsq3_D = new TH1F("hDPhiR_rsq3_D", "hDPhiR_rsq3_D", nbin, &xbins[0]); hDPhiR_rsq3_D->Sumw2();
  TH1F *hDPhiR_rsq3_Q = new TH1F("hDPhiR_rsq3_Q", "hDPhiR_rsq3_Q", nbin, &xbins[0]); hDPhiR_rsq3_Q->Sumw2();
  TH1F *hDPhiR_rsq3_T = new TH1F("hDPhiR_rsq3_T", "hDPhiR_rsq3_T", nbin, &xbins[0]); hDPhiR_rsq3_T->Sumw2();
  TH1F *hDPhiR_rsq3_W = new TH1F("hDPhiR_rsq3_W", "hDPhiR_rsq3_W", nbin, &xbins[0]); hDPhiR_rsq3_W->Sumw2();
  TH1F *hDPhiR_rsq3_Z = new TH1F("hDPhiR_rsq3_Z", "hDPhiR_rsq3_Z", nbin, &xbins[0]); hDPhiR_rsq3_Z->Sumw2();

  TH1F *hDPhiR_ljpt1_D = new TH1F("hDPhiR_ljpt1_D", "hDPhiR_ljpt1_D", nbin, &xbins[0]); hDPhiR_ljpt1_D->Sumw2();
  TH1F *hDPhiR_ljpt1_Q = new TH1F("hDPhiR_ljpt1_Q", "hDPhiR_ljpt1_Q", nbin, &xbins[0]); hDPhiR_ljpt1_Q->Sumw2();
  TH1F *hDPhiR_ljpt1_T = new TH1F("hDPhiR_ljpt1_T", "hDPhiR_ljpt1_T", nbin, &xbins[0]); hDPhiR_ljpt1_T->Sumw2();
  TH1F *hDPhiR_ljpt1_W = new TH1F("hDPhiR_ljpt1_W", "hDPhiR_ljpt1_W", nbin, &xbins[0]); hDPhiR_ljpt1_W->Sumw2();
  TH1F *hDPhiR_ljpt1_Z = new TH1F("hDPhiR_ljpt1_Z", "hDPhiR_ljpt1_Z", nbin, &xbins[0]); hDPhiR_ljpt1_Z->Sumw2();

  TH1F *hDPhiR_ljpt2_D = new TH1F("hDPhiR_ljpt2_D", "hDPhiR_ljpt2_D", nbin, &xbins[0]); hDPhiR_ljpt2_D->Sumw2();
  TH1F *hDPhiR_ljpt2_Q = new TH1F("hDPhiR_ljpt2_Q", "hDPhiR_ljpt2_Q", nbin, &xbins[0]); hDPhiR_ljpt2_Q->Sumw2();
  TH1F *hDPhiR_ljpt2_T = new TH1F("hDPhiR_ljpt2_T", "hDPhiR_ljpt2_T", nbin, &xbins[0]); hDPhiR_ljpt2_T->Sumw2();
  TH1F *hDPhiR_ljpt2_W = new TH1F("hDPhiR_ljpt2_W", "hDPhiR_ljpt2_W", nbin, &xbins[0]); hDPhiR_ljpt2_W->Sumw2();
  TH1F *hDPhiR_ljpt2_Z = new TH1F("hDPhiR_ljpt2_Z", "hDPhiR_ljpt2_Z", nbin, &xbins[0]); hDPhiR_ljpt2_Z->Sumw2();

  TH1F *hDPhiR_ljpt3_D = new TH1F("hDPhiR_ljpt3_D", "hDPhiR_ljpt3_D", nbin, &xbins[0]); hDPhiR_ljpt3_D->Sumw2();
  TH1F *hDPhiR_ljpt3_Q = new TH1F("hDPhiR_ljpt3_Q", "hDPhiR_ljpt3_Q", nbin, &xbins[0]); hDPhiR_ljpt3_Q->Sumw2();
  TH1F *hDPhiR_ljpt3_T = new TH1F("hDPhiR_ljpt3_T", "hDPhiR_ljpt3_T", nbin, &xbins[0]); hDPhiR_ljpt3_T->Sumw2();
  TH1F *hDPhiR_ljpt3_W = new TH1F("hDPhiR_ljpt3_W", "hDPhiR_ljpt3_W", nbin, &xbins[0]); hDPhiR_ljpt3_W->Sumw2();
  TH1F *hDPhiR_ljpt3_Z = new TH1F("hDPhiR_ljpt3_Z", "hDPhiR_ljpt3_Z", nbin, &xbins[0]); hDPhiR_ljpt3_Z->Sumw2();

  TH1F *hDPhiR_ljpt4_D = new TH1F("hDPhiR_ljpt4_D", "hDPhiR_ljpt4_D", nbin, &xbins[0]); hDPhiR_ljpt4_D->Sumw2();
  TH1F *hDPhiR_ljpt4_Q = new TH1F("hDPhiR_ljpt4_Q", "hDPhiR_ljpt4_Q", nbin, &xbins[0]); hDPhiR_ljpt4_Q->Sumw2();
  TH1F *hDPhiR_ljpt4_T = new TH1F("hDPhiR_ljpt4_T", "hDPhiR_ljpt4_T", nbin, &xbins[0]); hDPhiR_ljpt4_T->Sumw2();
  TH1F *hDPhiR_ljpt4_W = new TH1F("hDPhiR_ljpt4_W", "hDPhiR_ljpt4_W", nbin, &xbins[0]); hDPhiR_ljpt4_W->Sumw2();
  TH1F *hDPhiR_ljpt4_Z = new TH1F("hDPhiR_ljpt4_Z", "hDPhiR_ljpt4_Z", nbin, &xbins[0]); hDPhiR_ljpt4_Z->Sumw2();

  TH1F *hDPhiR_ljpt5_D = new TH1F("hDPhiR_ljpt5_D", "hDPhiR_ljpt5_D", nbin, &xbins[0]); hDPhiR_ljpt5_D->Sumw2();
  TH1F *hDPhiR_ljpt5_Q = new TH1F("hDPhiR_ljpt5_Q", "hDPhiR_ljpt5_Q", nbin, &xbins[0]); hDPhiR_ljpt5_Q->Sumw2();
  TH1F *hDPhiR_ljpt5_T = new TH1F("hDPhiR_ljpt5_T", "hDPhiR_ljpt5_T", nbin, &xbins[0]); hDPhiR_ljpt5_T->Sumw2();
  TH1F *hDPhiR_ljpt5_W = new TH1F("hDPhiR_ljpt5_W", "hDPhiR_ljpt5_W", nbin, &xbins[0]); hDPhiR_ljpt5_W->Sumw2();
  TH1F *hDPhiR_ljpt5_Z = new TH1F("hDPhiR_ljpt5_Z", "hDPhiR_ljpt5_Z", nbin, &xbins[0]); hDPhiR_ljpt5_Z->Sumw2();

  TH1F *hDPhiR_mr1_D = new TH1F("hDPhiR_mr1_D", "hDPhiR_mr1_D", nbin, &xbins[0]); hDPhiR_mr1_D->Sumw2();
  TH1F *hDPhiR_mr1_Q = new TH1F("hDPhiR_mr1_Q", "hDPhiR_mr1_Q", nbin, &xbins[0]); hDPhiR_mr1_Q->Sumw2();
  TH1F *hDPhiR_mr1_T = new TH1F("hDPhiR_mr1_T", "hDPhiR_mr1_T", nbin, &xbins[0]); hDPhiR_mr1_T->Sumw2();
  TH1F *hDPhiR_mr1_W = new TH1F("hDPhiR_mr1_W", "hDPhiR_mr1_W", nbin, &xbins[0]); hDPhiR_mr1_W->Sumw2();
  TH1F *hDPhiR_mr1_Z = new TH1F("hDPhiR_mr1_Z", "hDPhiR_mr1_Z", nbin, &xbins[0]); hDPhiR_mr1_Z->Sumw2();

  TH1F *hDPhiR_mr2_D = new TH1F("hDPhiR_mr2_D", "hDPhiR_mr2_D", nbin, &xbins[0]); hDPhiR_mr2_D->Sumw2();
  TH1F *hDPhiR_mr2_Q = new TH1F("hDPhiR_mr2_Q", "hDPhiR_mr2_Q", nbin, &xbins[0]); hDPhiR_mr2_Q->Sumw2();
  TH1F *hDPhiR_mr2_T = new TH1F("hDPhiR_mr2_T", "hDPhiR_mr2_T", nbin, &xbins[0]); hDPhiR_mr2_T->Sumw2();
  TH1F *hDPhiR_mr2_W = new TH1F("hDPhiR_mr2_W", "hDPhiR_mr2_W", nbin, &xbins[0]); hDPhiR_mr2_W->Sumw2();
  TH1F *hDPhiR_mr2_Z = new TH1F("hDPhiR_mr2_Z", "hDPhiR_mr2_Z", nbin, &xbins[0]); hDPhiR_mr2_Z->Sumw2();

  TH1F *hDPhiR_mr3_D = new TH1F("hDPhiR_mr3_D", "hDPhiR_mr3_D", nbin, &xbins[0]); hDPhiR_mr3_D->Sumw2();
  TH1F *hDPhiR_mr3_Q = new TH1F("hDPhiR_mr3_Q", "hDPhiR_mr3_Q", nbin, &xbins[0]); hDPhiR_mr3_Q->Sumw2();
  TH1F *hDPhiR_mr3_T = new TH1F("hDPhiR_mr3_T", "hDPhiR_mr3_T", nbin, &xbins[0]); hDPhiR_mr3_T->Sumw2();
  TH1F *hDPhiR_mr3_W = new TH1F("hDPhiR_mr3_W", "hDPhiR_mr3_W", nbin, &xbins[0]); hDPhiR_mr3_W->Sumw2();
  TH1F *hDPhiR_mr3_Z = new TH1F("hDPhiR_mr3_Z", "hDPhiR_mr3_Z", nbin, &xbins[0]); hDPhiR_mr3_Z->Sumw2();

  TH1F *hDPhiR_mr4_D = new TH1F("hDPhiR_mr4_D", "hDPhiR_mr4_D", nbin, &xbins[0]); hDPhiR_mr4_D->Sumw2();
  TH1F *hDPhiR_mr4_Q = new TH1F("hDPhiR_mr4_Q", "hDPhiR_mr4_Q", nbin, &xbins[0]); hDPhiR_mr4_Q->Sumw2();
  TH1F *hDPhiR_mr4_T = new TH1F("hDPhiR_mr4_T", "hDPhiR_mr4_T", nbin, &xbins[0]); hDPhiR_mr4_T->Sumw2();
  TH1F *hDPhiR_mr4_W = new TH1F("hDPhiR_mr4_W", "hDPhiR_mr4_W", nbin, &xbins[0]); hDPhiR_mr4_W->Sumw2();
  TH1F *hDPhiR_mr4_Z = new TH1F("hDPhiR_mr4_Z", "hDPhiR_mr4_Z", nbin, &xbins[0]); hDPhiR_mr4_Z->Sumw2();

  TH1F *hDPhiR_mr5_D = new TH1F("hDPhiR_mr5_D", "hDPhiR_mr5_D", nbin, &xbins[0]); hDPhiR_mr5_D->Sumw2();
  TH1F *hDPhiR_mr5_Q = new TH1F("hDPhiR_mr5_Q", "hDPhiR_mr5_Q", nbin, &xbins[0]); hDPhiR_mr5_Q->Sumw2();
  TH1F *hDPhiR_mr5_T = new TH1F("hDPhiR_mr5_T", "hDPhiR_mr5_T", nbin, &xbins[0]); hDPhiR_mr5_T->Sumw2();
  TH1F *hDPhiR_mr5_W = new TH1F("hDPhiR_mr5_W", "hDPhiR_mr5_W", nbin, &xbins[0]); hDPhiR_mr5_W->Sumw2();
  TH1F *hDPhiR_mr5_Z = new TH1F("hDPhiR_mr5_Z", "hDPhiR_mr5_Z", nbin, &xbins[0]); hDPhiR_mr5_Z->Sumw2();

  // open trees
  //TTree *tD = (TTree*) fD->Get("QCDTree");
  //TTree *tQ = (TTree*) fQ->Get("QCDTree");
  //TTree *tT = (TTree*) fT->Get("QCDTree");
  //TTree *tW = (TTree*) fW->Get("QCDTree");
  //TTree *tZ = (TTree*) fZ->Get("QCDTree");

  TTree *tD = (TTree*) fD->Get("RazorInclusive");
  TTree *tQ = (TTree*) fQ->Get("RazorInclusive");
  TTree *tT = (TTree*) fT->Get("RazorInclusive");
  TTree *tW = (TTree*) fW->Get("RazorInclusive");
  TTree *tZ = (TTree*) fZ->Get("RazorInclusive");

  //--------------------------------------------------------------
  //
  // hopefully no configuration below
  //
  //--------------------------------------------------------------

  // draw MR (again)
  tD->Draw("MR>>hMR_D", cut_str_rsq);
  tQ->Draw("MR>>hMR_Q", cut_str_rsq+cut_str_mc_qcd);
  tT->Draw("MR>>hMR_T", cut_str_rsq+cut_str_mc);
  tW->Draw("MR>>hMR_W", cut_str_rsq+cut_str_mc);
  tZ->Draw("MR>>hMR_Z", cut_str_rsq+cut_str_mc);
  cout << hMR_D->GetEffectiveEntries() << ", " << hMR_Q->GetEffectiveEntries() << ", " << hMR_T->GetEffectiveEntries() << ", " << hMR_W->GetEffectiveEntries() << ", " << hMR_Z->GetEffectiveEntries() << endl;
}
/*
void foobar() {
  tD->Draw("MR>>hMR_pass_D", cut_str_rsq+"*(abs(dPhiRazor)>2.8)");
  tQ->Draw("MR>>hMR_pass_Q", cut_str_rsq+"*(abs(dPhiRazor)>2.8)"+cut_str_mc_qcd);
  tT->Draw("MR>>hMR_pass_T", cut_str_rsq+"*(abs(dPhiRazor)>2.8)"+cut_str_mc);
  tW->Draw("MR>>hMR_pass_W", cut_str_rsq+"*(abs(dPhiRazor)>2.8)"+cut_str_mc);
  tZ->Draw("MR>>hMR_pass_Z", cut_str_rsq+"*(abs(dPhiRazor)>2.8)"+cut_str_mc);
  
  tD->Draw("MR>>hMR_fail_D", cut_str_rsq+"*(abs(dPhiRazor)<2.8)");
  tQ->Draw("MR>>hMR_fail_Q", cut_str_rsq+"*(abs(dPhiRazor)<2.8)"+cut_str_mc_qcd);
  tT->Draw("MR>>hMR_fail_T", cut_str_rsq+"*(abs(dPhiRazor)<2.8)"+cut_str_mc);
  tW->Draw("MR>>hMR_fail_W", cut_str_rsq+"*(abs(dPhiRazor)<2.8)"+cut_str_mc);
  tZ->Draw("MR>>hMR_fail_Z", cut_str_rsq+"*(abs(dPhiRazor)<2.8)"+cut_str_mc);

  // draw Rsq (again)
  tD->Draw("Rsq>>hRsq_D", cut_str_rsq);
  tQ->Draw("Rsq>>hRsq_Q", cut_str+cut_str_mc_qcd);
  tT->Draw("Rsq>>hRsq_T", cut_str+cut_str_mc);
  tW->Draw("Rsq>>hRsq_W", cut_str+cut_str_mc);
  tZ->Draw("Rsq>>hRsq_Z", cut_str+cut_str_mc);
  
  tD->Draw("Rsq>>hRsq_pass_D", cut_str_rsq+"*(abs(dPhiRazor)>2.8)");
  tQ->Draw("Rsq>>hRsq_pass_Q", cut_str+"*(abs(dPhiRazor)>2.8)"+cut_str_mc_qcd);
  tT->Draw("Rsq>>hRsq_pass_T", cut_str+"*(abs(dPhiRazor)>2.8)"+cut_str_mc);
  tW->Draw("Rsq>>hRsq_pass_W", cut_str+"*(abs(dPhiRazor)>2.8)"+cut_str_mc);
  tZ->Draw("Rsq>>hRsq_pass_Z", cut_str+"*(abs(dPhiRazor)>2.8)"+cut_str_mc);
  
  tD->Draw("Rsq>>hRsq_fail_D", cut_str_rsq+"*(abs(dPhiRazor)<2.8)");
  tQ->Draw("Rsq>>hRsq_fail_Q", cut_str+"*(abs(dPhiRazor)<2.8)"+cut_str_mc_qcd);
  tT->Draw("Rsq>>hRsq_fail_T", cut_str+"*(abs(dPhiRazor)<2.8)"+cut_str_mc);
  tW->Draw("Rsq>>hRsq_fail_W", cut_str+"*(abs(dPhiRazor)<2.8)"+cut_str_mc);
  tZ->Draw("Rsq>>hRsq_fail_Z", cut_str+"*(abs(dPhiRazor)<2.8)"+cut_str_mc);

  // draw lead jet pT (again)
  tD->Draw("leadingJetPt>>hLjpt_D", cut_str_rsq);
  tQ->Draw("leadingJetPt>>hLjpt_Q", cut_str_rsq+cut_str_mc_qcd);
  tT->Draw("leadingJetPt>>hLjpt_T", cut_str_rsq+cut_str_mc);
  tW->Draw("leadingJetPt>>hLjpt_W", cut_str_rsq+cut_str_mc);
  tZ->Draw("leadingJetPt>>hLjpt_Z", cut_str_rsq+cut_str_mc);

  tD->Draw("leadingJetPt>>hLjpt_pass_D", cut_str_rsq+"*(abs(dPhiRazor)>2.8)");
  tQ->Draw("leadingJetPt>>hLjpt_pass_Q", cut_str_rsq+"*(abs(dPhiRazor)>2.8)"+cut_str_mc_qcd);
  tT->Draw("leadingJetPt>>hLjpt_pass_T", cut_str_rsq+"*(abs(dPhiRazor)>2.8)"+cut_str_mc);
  tW->Draw("leadingJetPt>>hLjpt_pass_W", cut_str_rsq+"*(abs(dPhiRazor)>2.8)"+cut_str_mc);
  tZ->Draw("leadingJetPt>>hLjpt_pass_Z", cut_str_rsq+"*(abs(dPhiRazor)>2.8)"+cut_str_mc);

  tD->Draw("leadingJetPt>>hLjpt_fail_D", cut_str_rsq+"*(abs(dPhiRazor)<2.8)");
  tQ->Draw("leadingJetPt>>hLjpt_fail_Q", cut_str_rsq+"*(abs(dPhiRazor)<2.8)"+cut_str_mc_qcd);
  tT->Draw("leadingJetPt>>hLjpt_fail_T", cut_str_rsq+"*(abs(dPhiRazor)<2.8)"+cut_str_mc);
  tW->Draw("leadingJetPt>>hLjpt_fail_W", cut_str_rsq+"*(abs(dPhiRazor)<2.8)"+cut_str_mc);
  tZ->Draw("leadingJetPt>>hLjpt_fail_Z", cut_str_rsq+"*(abs(dPhiRazor)<2.8)"+cut_str_mc);

  // draw dPhiR (again)
  tD->Draw("abs(dPhiRazor)>>hDPhiR_D", cut_str_rsq);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_Q", cut_str_rsq+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_T", cut_str_rsq+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_W", cut_str_rsq+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_Z", cut_str_rsq+cut_str_mc);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_rsq1_D", cut_str+cut_str_rsq1);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_rsq1_Q", cut_str+cut_str_rsq1+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_rsq1_T", cut_str+cut_str_rsq1+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_rsq1_W", cut_str+cut_str_rsq1+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_rsq1_Z", cut_str+cut_str_rsq1+cut_str_mc);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_rsq2_D", cut_str+cut_str_rsq2);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_rsq2_Q", cut_str+cut_str_rsq2+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_rsq2_T", cut_str+cut_str_rsq2+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_rsq2_W", cut_str+cut_str_rsq2+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_rsq2_Z", cut_str+cut_str_rsq2+cut_str_mc);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_rsq3_D", cut_str+cut_str_rsq3);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_rsq3_Q", cut_str+cut_str_rsq3+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_rsq3_T", cut_str+cut_str_rsq3+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_rsq3_W", cut_str+cut_str_rsq3+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_rsq3_Z", cut_str+cut_str_rsq3+cut_str_mc);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_ljpt1_D", cut_str_rsq+cut_str_ljpt1);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_ljpt1_Q", cut_str_rsq+cut_str_ljpt1+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_ljpt1_T", cut_str_rsq+cut_str_ljpt1+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_ljpt1_W", cut_str_rsq+cut_str_ljpt1+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_ljpt1_Z", cut_str_rsq+cut_str_ljpt1+cut_str_mc);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_ljpt2_D", cut_str_rsq+cut_str_ljpt2);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_ljpt2_Q", cut_str_rsq+cut_str_ljpt2+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_ljpt2_T", cut_str_rsq+cut_str_ljpt2+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_ljpt2_W", cut_str_rsq+cut_str_ljpt2+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_ljpt2_Z", cut_str_rsq+cut_str_ljpt2+cut_str_mc);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_ljpt3_D", cut_str_rsq+cut_str_ljpt3);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_ljpt3_Q", cut_str_rsq+cut_str_ljpt3+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_ljpt3_T", cut_str_rsq+cut_str_ljpt3+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_ljpt3_W", cut_str_rsq+cut_str_ljpt3+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_ljpt3_Z", cut_str_rsq+cut_str_ljpt3+cut_str_mc);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_ljpt4_D", cut_str_rsq+cut_str_ljpt4);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_ljpt4_Q", cut_str_rsq+cut_str_ljpt4+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_ljpt4_T", cut_str_rsq+cut_str_ljpt4+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_ljpt4_W", cut_str_rsq+cut_str_ljpt4+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_ljpt4_Z", cut_str_rsq+cut_str_ljpt4+cut_str_mc);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_ljpt5_D", cut_str_rsq+cut_str_ljpt5);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_ljpt5_Q", cut_str_rsq+cut_str_ljpt5+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_ljpt5_T", cut_str_rsq+cut_str_ljpt5+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_ljpt5_W", cut_str_rsq+cut_str_ljpt5+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_ljpt5_Z", cut_str_rsq+cut_str_ljpt5+cut_str_mc);

  tD->Draw("abs(dPhiRazor)>>hDPhiR_mr1_D", cut_str_rsq+cut_str_mr1);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_mr1_Q", cut_str_rsq+cut_str_mr1+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_mr1_T", cut_str_rsq+cut_str_mr1+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_mr1_W", cut_str_rsq+cut_str_mr1+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_mr1_Z", cut_str_rsq+cut_str_mr1+cut_str_mc);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_mr2_D", cut_str_rsq+cut_str_mr2);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_mr2_Q", cut_str_rsq+cut_str_mr2+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_mr2_T", cut_str_rsq+cut_str_mr2+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_mr2_W", cut_str_rsq+cut_str_mr2+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_mr2_Z", cut_str_rsq+cut_str_mr2+cut_str_mc);
  
  tD->Draw("abs(dPhiRazor)>>hDPhiR_mr3_D", cut_str_rsq+cut_str_mr3);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_mr3_Q", cut_str_rsq+cut_str_mr3+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_mr3_T", cut_str_rsq+cut_str_mr3+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_mr3_W", cut_str_rsq+cut_str_mr3+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_mr3_Z", cut_str_rsq+cut_str_mr3+cut_str_mc);

  tD->Draw("abs(dPhiRazor)>>hDPhiR_mr4_D", cut_str_rsq+cut_str_mr4);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_mr4_Q", cut_str_rsq+cut_str_mr4+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_mr4_T", cut_str_rsq+cut_str_mr4+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_mr4_W", cut_str_rsq+cut_str_mr4+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_mr4_Z", cut_str_rsq+cut_str_mr4+cut_str_mc);

  tD->Draw("abs(dPhiRazor)>>hDPhiR_mr5_D", cut_str_rsq+cut_str_mr5);
  tQ->Draw("abs(dPhiRazor)>>hDPhiR_mr5_Q", cut_str_rsq+cut_str_mr5+cut_str_mc_qcd);
  tT->Draw("abs(dPhiRazor)>>hDPhiR_mr5_T", cut_str_rsq+cut_str_mr5+cut_str_mc);
  tW->Draw("abs(dPhiRazor)>>hDPhiR_mr5_W", cut_str_rsq+cut_str_mr5+cut_str_mc);
  tZ->Draw("abs(dPhiRazor)>>hDPhiR_mr5_Z", cut_str_rsq+cut_str_mr5+cut_str_mc);

  //--------------------------------------------------------------
  //
  // configure draw options
  //
  //--------------------------------------------------------------

  //data
  InitData(hMR_D, "", "", kBlack);
  InitData(hMR_pass_D, "", "", kBlack);
  InitData(hMR_fail_D, "", "", kBlack);
  InitData(hRsq_D, "", "", kBlack);
  InitData(hRsq_pass_D, "", "", kBlack);
  InitData(hRsq_fail_D, "", "", kBlack);
  InitData(hLjpt_D, "", "", kBlack);
  InitData(hLjpt_pass_D, "", "", kBlack);
  InitData(hLjpt_fail_D, "", "", kBlack);
  InitData(hDPhiR_D, "", "", kBlack);

  InitData(hDPhiR_rsq1_D, "", "", kBlack);
  InitData(hDPhiR_rsq2_D, "", "", kBlack);
  InitData(hDPhiR_rsq3_D, "", "", kBlack);

  InitData(hDPhiR_ljpt1_D, "", "", kBlack);
  InitData(hDPhiR_ljpt2_D, "", "", kBlack);
  InitData(hDPhiR_ljpt3_D, "", "", kBlack);
  InitData(hDPhiR_ljpt4_D, "", "", kBlack);
  InitData(hDPhiR_ljpt5_D, "", "", kBlack);

  InitData(hDPhiR_ljpt1_D, "", "", kBlack);
  InitData(hDPhiR_ljpt2_D, "", "", kBlack);
  InitData(hDPhiR_ljpt3_D, "", "", kBlack);
  InitData(hDPhiR_ljpt4_D, "", "", kBlack);
  InitData(hDPhiR_ljpt5_D, "", "", kBlack);

  // qcd
  InitHist(hMR_Q, "", "", kMagenta);
  InitHist(hMR_pass_Q, "", "", kMagenta);
  InitHist(hMR_fail_Q, "", "", kMagenta);
  InitHist(hRsq_Q, "", "", kMagenta);
  InitHist(hRsq_pass_Q, "", "", kMagenta);
  InitHist(hRsq_fail_Q, "", "", kMagenta);
  InitHist(hLjpt_Q, "", "", kMagenta);
  InitHist(hLjpt_pass_Q, "", "", kMagenta);
  InitHist(hLjpt_fail_Q, "", "", kMagenta);
  InitHist(hDPhiR_Q, "", "", kMagenta);  

  InitHist(hDPhiR_rsq1_Q, "", "", kMagenta);  
  InitHist(hDPhiR_rsq2_Q, "", "", kMagenta);  
  InitHist(hDPhiR_rsq3_Q, "", "", kMagenta);  

  InitHist(hDPhiR_ljpt1_Q, "", "", kMagenta);  
  InitHist(hDPhiR_ljpt2_Q, "", "", kMagenta);  
  InitHist(hDPhiR_ljpt3_Q, "", "", kMagenta);  
  InitHist(hDPhiR_ljpt4_Q, "", "", kMagenta);  
  InitHist(hDPhiR_ljpt5_Q, "", "", kMagenta);  

  InitHist(hDPhiR_mr1_Q, "", "", kMagenta);  
  InitHist(hDPhiR_mr2_Q, "", "", kMagenta);  
  InitHist(hDPhiR_mr3_Q, "", "", kMagenta);  
  InitHist(hDPhiR_mr4_Q, "", "", kMagenta);  
  InitHist(hDPhiR_mr5_Q, "", "", kMagenta);

  // tt
  InitHist(hMR_T, "", "", kGreen+2);
  InitHist(hMR_pass_T, "", "", kGreen+2);
  InitHist(hMR_fail_T, "", "", kGreen+2);
  InitHist(hRsq_T, "", "", kGreen+2);
  InitHist(hRsq_pass_T, "", "", kGreen+2);
  InitHist(hRsq_fail_T, "", "", kGreen+2);
  InitHist(hLjpt_T, "", "", kGreen+2);
  InitHist(hLjpt_pass_T, "", "", kGreen+2);
  InitHist(hLjpt_fail_T, "", "", kGreen+2);
  InitHist(hDPhiR_T, "", "", kGreen+2);  

  InitHist(hDPhiR_rsq1_T, "", "", kGreen+2);  
  InitHist(hDPhiR_rsq2_T, "", "", kGreen+2);  
  InitHist(hDPhiR_rsq3_T, "", "", kGreen+2);  

  InitHist(hDPhiR_ljpt1_T, "", "", kGreen+2);  
  InitHist(hDPhiR_ljpt2_T, "", "", kGreen+2);  
  InitHist(hDPhiR_ljpt3_T, "", "", kGreen+2);  
  InitHist(hDPhiR_ljpt4_T, "", "", kGreen+2);  
  InitHist(hDPhiR_ljpt5_T, "", "", kGreen+2);  

  InitHist(hDPhiR_mr1_T, "", "", kGreen+2);  
  InitHist(hDPhiR_mr2_T, "", "", kGreen+2);  
  InitHist(hDPhiR_mr3_T, "", "", kGreen+2);  
  InitHist(hDPhiR_mr4_T, "", "", kGreen+2);  
  InitHist(hDPhiR_mr5_T, "", "", kGreen+2);

  // z
  InitHist(hMR_Z, "", "", kBlue+1);
  InitHist(hMR_pass_Z, "", "", kBlue+1);
  InitHist(hMR_fail_Z, "", "", kBlue+1);
  InitHist(hRsq_Z, "", "", kBlue+1);
  InitHist(hRsq_pass_Z, "", "", kBlue+1);
  InitHist(hRsq_fail_Z, "", "", kBlue+1);
  InitHist(hLjpt_Z, "", "", kBlue+1);
  InitHist(hLjpt_pass_Z, "", "", kBlue+1);
  InitHist(hLjpt_fail_Z, "", "", kBlue+1);
  InitHist(hDPhiR_Z, "", "", kBlue+1);  

  InitHist(hDPhiR_rsq1_Z, "", "", kBlue+1);  
  InitHist(hDPhiR_rsq2_Z, "", "", kBlue+1);  
  InitHist(hDPhiR_rsq3_Z, "", "", kBlue+1);  

  InitHist(hDPhiR_ljpt1_Z, "", "", kBlue+1);  
  InitHist(hDPhiR_ljpt2_Z, "", "", kBlue+1);  
  InitHist(hDPhiR_ljpt3_Z, "", "", kBlue+1);  
  InitHist(hDPhiR_ljpt4_Z, "", "", kBlue+1);  
  InitHist(hDPhiR_ljpt5_Z, "", "", kBlue+1);  

  InitHist(hDPhiR_mr1_Z, "", "", kBlue+1);  
  InitHist(hDPhiR_mr2_Z, "", "", kBlue+1);  
  InitHist(hDPhiR_mr3_Z, "", "", kBlue+1);  
  InitHist(hDPhiR_mr4_Z, "", "", kBlue+1);  
  InitHist(hDPhiR_mr5_Z, "", "", kBlue+1);

  // w
  InitHist(hMR_W, "", "", kRed+1);
  InitHist(hMR_pass_W, "", "", kRed+1);
  InitHist(hMR_fail_W, "", "", kRed+1);
  InitHist(hRsq_W, "", "", kRed+1);
  InitHist(hRsq_pass_W, "", "", kRed+1);
  InitHist(hRsq_fail_W, "", "", kRed+1);
  InitHist(hLjpt_W, "", "", kRed+1);
  InitHist(hLjpt_pass_W, "", "", kRed+1);
  InitHist(hLjpt_fail_W, "", "", kRed+1);
  InitHist(hDPhiR_W, "", "", kRed+1);  

  InitHist(hDPhiR_rsq1_W, "", "", kRed+1);  
  InitHist(hDPhiR_rsq2_W, "", "", kRed+1);  
  InitHist(hDPhiR_rsq3_W, "", "", kRed+1);  

  InitHist(hDPhiR_ljpt1_W, "", "", kRed+1);  
  InitHist(hDPhiR_ljpt2_W, "", "", kRed+1);  
  InitHist(hDPhiR_ljpt3_W, "", "", kRed+1);  
  InitHist(hDPhiR_ljpt4_W, "", "", kRed+1);  
  InitHist(hDPhiR_ljpt5_W, "", "", kRed+1);  

  InitHist(hDPhiR_mr1_W, "", "", kRed+1);  
  InitHist(hDPhiR_mr2_W, "", "", kRed+1);  
  InitHist(hDPhiR_mr3_W, "", "", kRed+1);  
  InitHist(hDPhiR_mr4_W, "", "", kRed+1);  
  InitHist(hDPhiR_mr5_W, "", "", kRed+1); 

  hMR_Q->Add(hMR_T);
  hMR_Q->Add(hMR_W);
  hMR_Q->Add(hMR_Z);
  hMR_T->Add(hMR_W);
  hMR_T->Add(hMR_Z);
  hMR_W->Add(hMR_Z);

  Double_t scale=hMR_D->Integral()/hMR_Q->Integral();
  //hMR_Q->Scale(scale);
  //hMR_T->Scale(scale);
  //hMR_W->Scale(scale);
  //hMR_Z->Scale(scale);

  hMR_pass_Q->Add(hMR_pass_T);
  hMR_pass_Q->Add(hMR_pass_W);
  hMR_pass_Q->Add(hMR_pass_Z);
  hMR_pass_T->Add(hMR_pass_W);
  hMR_pass_T->Add(hMR_pass_Z);
  hMR_pass_W->Add(hMR_pass_Z);

  scale=hMR_pass_D->Integral()/hMR_pass_Q->Integral();
  //hMR_pass_Q->Scale(scale);
  //hMR_pass_T->Scale(scale);
  //hMR_pass_W->Scale(scale);
  //hMR_pass_Z->Scale(scale);

  hMR_fail_Q->Add(hMR_fail_T);
  hMR_fail_Q->Add(hMR_fail_W);
  hMR_fail_Q->Add(hMR_fail_Z);
  hMR_fail_T->Add(hMR_fail_W);
  hMR_fail_T->Add(hMR_fail_Z);
  hMR_fail_W->Add(hMR_fail_Z);

  scale=hMR_fail_D->Integral()/hMR_fail_Q->Integral();
  //hMR_fail_Q->Scale(scale);
  //hMR_fail_T->Scale(scale);
  //hMR_fail_W->Scale(scale);
  //hMR_fail_Z->Scale(scale);

  hRsq_Q->Add(hRsq_T);
  hRsq_Q->Add(hRsq_W);
  hRsq_Q->Add(hRsq_Z);
  hRsq_T->Add(hRsq_W);
  hRsq_T->Add(hRsq_Z);
  hRsq_W->Add(hRsq_Z);

  scale=hRsq_D->Integral()/hRsq_Q->Integral();
  //hRsq_Q->Scale(scale);
  //hRsq_T->Scale(scale);
  //hRsq_W->Scale(scale);
  //hRsq_Z->Scale(scale);

  hRsq_pass_Q->Add(hRsq_pass_T);
  hRsq_pass_Q->Add(hRsq_pass_W);
  hRsq_pass_Q->Add(hRsq_pass_Z);
  hRsq_pass_T->Add(hRsq_pass_W);
  hRsq_pass_T->Add(hRsq_pass_Z);
  hRsq_pass_W->Add(hRsq_pass_Z);

  scale=hRsq_pass_D->Integral()/hRsq_pass_Q->Integral();
  //hRsq_pass_Q->Scale(scale);
  //hRsq_pass_T->Scale(scale);
  //hRsq_pass_W->Scale(scale);
  //hRsq_pass_Z->Scale(scale);

  hRsq_fail_Q->Add(hRsq_fail_T);
  hRsq_fail_Q->Add(hRsq_fail_W);
  hRsq_fail_Q->Add(hRsq_fail_Z);
  hRsq_fail_T->Add(hRsq_fail_W);
  hRsq_fail_T->Add(hRsq_fail_Z);
  hRsq_fail_W->Add(hRsq_fail_Z);

  scale=hRsq_fail_D->Integral()/hRsq_fail_Q->Integral();
  //hRsq_fail_Q->Scale(scale);
  //hRsq_fail_T->Scale(scale);
  //hRsq_fail_W->Scale(scale);
  //hRsq_fail_Z->Scale(scale);

  hLjpt_Q->Add(hLjpt_T);
  hLjpt_Q->Add(hLjpt_W);
  hLjpt_Q->Add(hLjpt_Z);
  hLjpt_T->Add(hLjpt_W);
  hLjpt_T->Add(hLjpt_Z);
  hLjpt_W->Add(hLjpt_Z);

  scale=hLjpt_D->Integral()/hLjpt_Q->Integral();
  //hLjpt_Q->Scale(scale);
  //hLjpt_T->Scale(scale);
  //hLjpt_W->Scale(scale);
  //hLjpt_Z->Scale(scale);

  hLjpt_pass_Q->Add(hLjpt_pass_T);
  hLjpt_pass_Q->Add(hLjpt_pass_W);
  hLjpt_pass_Q->Add(hLjpt_pass_Z);
  hLjpt_pass_T->Add(hLjpt_pass_W);
  hLjpt_pass_T->Add(hLjpt_pass_Z);
  hLjpt_pass_W->Add(hLjpt_pass_Z);

  scale=hLjpt_pass_D->Integral()/hLjpt_pass_Q->Integral();
  //hLjpt_pass_Q->Scale(scale);
  //hLjpt_pass_T->Scale(scale);
  //hLjpt_pass_W->Scale(scale);
  //hLjpt_pass_Z->Scale(scale);

  hLjpt_fail_Q->Add(hLjpt_fail_T);
  hLjpt_fail_Q->Add(hLjpt_fail_W);
  hLjpt_fail_Q->Add(hLjpt_fail_Z);
  hLjpt_fail_T->Add(hLjpt_fail_W);
  hLjpt_fail_T->Add(hLjpt_fail_Z);
  hLjpt_fail_W->Add(hLjpt_fail_Z);

  scale=hLjpt_fail_D->Integral()/hLjpt_fail_Q->Integral();
  //hLjpt_fail_Q->Scale(scale);
  //hLjpt_fail_T->Scale(scale);
  //hLjpt_fail_W->Scale(scale);
  //hLjpt_fail_Z->Scale(scale);

  hDPhiR_Q->Add(hDPhiR_T);
  hDPhiR_Q->Add(hDPhiR_W);
  hDPhiR_Q->Add(hDPhiR_Z);
  hDPhiR_T->Add(hDPhiR_W);
  hDPhiR_T->Add(hDPhiR_Z);
  hDPhiR_W->Add(hDPhiR_Z);

  scale=hDPhiR_D->Integral()/hDPhiR_Q->Integral();
  //hDPhiR_Q->Scale(scale);
  //hDPhiR_T->Scale(scale);
  //hDPhiR_W->Scale(scale);
  //hDPhiR_Z->Scale(scale);

  hDPhiR_rsq1_Q->Add(hDPhiR_rsq1_T);
  hDPhiR_rsq1_Q->Add(hDPhiR_rsq1_W);
  hDPhiR_rsq1_Q->Add(hDPhiR_rsq1_Z);
  hDPhiR_rsq1_T->Add(hDPhiR_rsq1_W);
  hDPhiR_rsq1_T->Add(hDPhiR_rsq1_Z);
  hDPhiR_rsq1_W->Add(hDPhiR_rsq1_Z);

  scale=hDPhiR_rsq1_D->Integral()/hDPhiR_rsq1_Q->Integral();
  //hDPhiR_rsq1_Q->Scale(scale);
  //hDPhiR_rsq1_T->Scale(scale);
  //hDPhiR_rsq1_W->Scale(scale);
  //hDPhiR_rsq1_Z->Scale(scale);

  hDPhiR_rsq2_Q->Add(hDPhiR_rsq2_T);
  hDPhiR_rsq2_Q->Add(hDPhiR_rsq2_W);
  hDPhiR_rsq2_Q->Add(hDPhiR_rsq2_Z);
  hDPhiR_rsq2_T->Add(hDPhiR_rsq2_W);
  hDPhiR_rsq2_T->Add(hDPhiR_rsq2_Z);
  hDPhiR_rsq2_W->Add(hDPhiR_rsq2_Z);

  scale=hDPhiR_rsq2_D->Integral()/hDPhiR_rsq2_Q->Integral();
  //hDPhiR_rsq2_Q->Scale(scale);
  //hDPhiR_rsq2_T->Scale(scale);
  //hDPhiR_rsq2_W->Scale(scale);
  //hDPhiR_rsq2_Z->Scale(scale);

  hDPhiR_rsq3_Q->Add(hDPhiR_rsq3_T);
  hDPhiR_rsq3_Q->Add(hDPhiR_rsq3_W);
  hDPhiR_rsq3_Q->Add(hDPhiR_rsq3_Z);
  hDPhiR_rsq3_T->Add(hDPhiR_rsq3_W);
  hDPhiR_rsq3_T->Add(hDPhiR_rsq3_Z);
  hDPhiR_rsq3_W->Add(hDPhiR_rsq3_Z);

  scale=hDPhiR_rsq3_D->Integral()/hDPhiR_rsq3_Q->Integral();
  //hDPhiR_rsq3_Q->Scale(scale);
  //hDPhiR_rsq3_T->Scale(scale);
  //hDPhiR_rsq3_W->Scale(scale);
  //hDPhiR_rsq3_Z->Scale(scale);

  hDPhiR_ljpt1_Q->Add(hDPhiR_ljpt1_T);
  hDPhiR_ljpt1_Q->Add(hDPhiR_ljpt1_W);
  hDPhiR_ljpt1_Q->Add(hDPhiR_ljpt1_Z);
  hDPhiR_ljpt1_T->Add(hDPhiR_ljpt1_W);
  hDPhiR_ljpt1_T->Add(hDPhiR_ljpt1_Z);
  hDPhiR_ljpt1_W->Add(hDPhiR_ljpt1_Z);

  scale=hDPhiR_ljpt1_D->Integral()/hDPhiR_ljpt1_Q->Integral();
  //hDPhiR_ljpt1_Q->Scale(scale);
  //hDPhiR_ljpt1_T->Scale(scale);
  //hDPhiR_ljpt1_W->Scale(scale);
  //hDPhiR_ljpt1_Z->Scale(scale);

  hDPhiR_ljpt2_Q->Add(hDPhiR_ljpt2_T);
  hDPhiR_ljpt2_Q->Add(hDPhiR_ljpt2_W);
  hDPhiR_ljpt2_Q->Add(hDPhiR_ljpt2_Z);
  hDPhiR_ljpt2_T->Add(hDPhiR_ljpt2_W);
  hDPhiR_ljpt2_T->Add(hDPhiR_ljpt2_Z);
  hDPhiR_ljpt2_W->Add(hDPhiR_ljpt2_Z);

  scale=hDPhiR_ljpt2_D->Integral()/hDPhiR_ljpt2_Q->Integral();
  //hDPhiR_ljpt2_Q->Scale(scale);
  //hDPhiR_ljpt2_T->Scale(scale);
  //hDPhiR_ljpt2_W->Scale(scale);
  //hDPhiR_ljpt2_Z->Scale(scale);

  hDPhiR_ljpt3_Q->Add(hDPhiR_ljpt3_T);
  hDPhiR_ljpt3_Q->Add(hDPhiR_ljpt3_W);
  hDPhiR_ljpt3_Q->Add(hDPhiR_ljpt3_Z);
  hDPhiR_ljpt3_T->Add(hDPhiR_ljpt3_W);
  hDPhiR_ljpt3_T->Add(hDPhiR_ljpt3_Z);
  hDPhiR_ljpt3_W->Add(hDPhiR_ljpt3_Z);

  scale=hDPhiR_ljpt3_D->Integral()/hDPhiR_ljpt3_Q->Integral();
  //hDPhiR_ljpt3_Q->Scale(scale);
  //hDPhiR_ljpt3_T->Scale(scale);
  //hDPhiR_ljpt3_W->Scale(scale);
  //hDPhiR_ljpt3_Z->Scale(scale);

  hDPhiR_ljpt4_Q->Add(hDPhiR_ljpt4_T);
  hDPhiR_ljpt4_Q->Add(hDPhiR_ljpt4_W);
  hDPhiR_ljpt4_Q->Add(hDPhiR_ljpt4_Z);
  hDPhiR_ljpt4_T->Add(hDPhiR_ljpt4_W);
  hDPhiR_ljpt4_T->Add(hDPhiR_ljpt4_Z);
  hDPhiR_ljpt4_W->Add(hDPhiR_ljpt4_Z);

  scale=hDPhiR_ljpt4_D->Integral()/hDPhiR_ljpt4_Q->Integral();
  //hDPhiR_ljpt4_Q->Scale(scale);
  //hDPhiR_ljpt4_T->Scale(scale);
  //hDPhiR_ljpt4_W->Scale(scale);
  //hDPhiR_ljpt4_Z->Scale(scale);

  hDPhiR_ljpt5_Q->Add(hDPhiR_ljpt5_T);
  hDPhiR_ljpt5_Q->Add(hDPhiR_ljpt5_W);
  hDPhiR_ljpt5_Q->Add(hDPhiR_ljpt5_Z);
  hDPhiR_ljpt5_T->Add(hDPhiR_ljpt5_W);
  hDPhiR_ljpt5_T->Add(hDPhiR_ljpt5_Z);
  hDPhiR_ljpt5_W->Add(hDPhiR_ljpt5_Z);

  scale=hDPhiR_ljpt5_D->Integral()/hDPhiR_ljpt5_Q->Integral();
  //hDPhiR_ljpt5_Q->Scale(scale);
  //hDPhiR_ljpt5_T->Scale(scale);
  //hDPhiR_ljpt5_W->Scale(scale);
  //hDPhiR_ljpt5_Z->Scale(scale);

  hDPhiR_mr1_Q->Add(hDPhiR_mr1_T);
  hDPhiR_mr1_Q->Add(hDPhiR_mr1_W);
  hDPhiR_mr1_Q->Add(hDPhiR_mr1_Z);
  hDPhiR_mr1_T->Add(hDPhiR_mr1_W);
  hDPhiR_mr1_T->Add(hDPhiR_mr1_Z);
  hDPhiR_mr1_W->Add(hDPhiR_mr1_Z);

  scale=hDPhiR_mr1_D->Integral()/hDPhiR_mr1_Q->Integral();
  //hDPhiR_mr1_Q->Scale(scale);
  //hDPhiR_mr1_T->Scale(scale);
  //hDPhiR_mr1_W->Scale(scale);
  //hDPhiR_mr1_Z->Scale(scale);

  hDPhiR_mr2_Q->Add(hDPhiR_mr2_T);
  hDPhiR_mr2_Q->Add(hDPhiR_mr2_W);
  hDPhiR_mr2_Q->Add(hDPhiR_mr2_Z);
  hDPhiR_mr2_T->Add(hDPhiR_mr2_W);
  hDPhiR_mr2_T->Add(hDPhiR_mr2_Z);
  hDPhiR_mr2_W->Add(hDPhiR_mr2_Z);

  scale=hDPhiR_mr2_D->Integral()/hDPhiR_mr2_Q->Integral();
  //hDPhiR_mr2_Q->Scale(scale);
  //hDPhiR_mr2_T->Scale(scale);
  //hDPhiR_mr2_W->Scale(scale);
  //hDPhiR_mr2_Z->Scale(scale);

  hDPhiR_mr3_Q->Add(hDPhiR_mr3_T);
  hDPhiR_mr3_Q->Add(hDPhiR_mr3_W);
  hDPhiR_mr3_Q->Add(hDPhiR_mr3_Z);
  hDPhiR_mr3_T->Add(hDPhiR_mr3_W);
  hDPhiR_mr3_T->Add(hDPhiR_mr3_Z);
  hDPhiR_mr3_W->Add(hDPhiR_mr3_Z);

  scale=hDPhiR_mr3_D->Integral()/hDPhiR_mr3_Q->Integral();
  //hDPhiR_mr3_Q->Scale(scale);
  //hDPhiR_mr3_T->Scale(scale);
  //hDPhiR_mr3_W->Scale(scale);
  //hDPhiR_mr3_Z->Scale(scale);

  hDPhiR_mr4_Q->Add(hDPhiR_mr4_T);
  hDPhiR_mr4_Q->Add(hDPhiR_mr4_W);
  hDPhiR_mr4_Q->Add(hDPhiR_mr4_Z);
  hDPhiR_mr4_T->Add(hDPhiR_mr4_W);
  hDPhiR_mr4_T->Add(hDPhiR_mr4_Z);
  hDPhiR_mr4_W->Add(hDPhiR_mr4_Z);

  scale=hDPhiR_mr4_D->Integral()/hDPhiR_mr4_Q->Integral();
  //hDPhiR_mr4_Q->Scale(scale);
  //hDPhiR_mr4_T->Scale(scale);
  //hDPhiR_mr4_W->Scale(scale);
  //hDPhiR_mr4_Z->Scale(scale);

  hDPhiR_mr5_Q->Add(hDPhiR_mr5_T);
  hDPhiR_mr5_Q->Add(hDPhiR_mr5_W);
  hDPhiR_mr5_Q->Add(hDPhiR_mr5_Z);
  hDPhiR_mr5_T->Add(hDPhiR_mr5_W);
  hDPhiR_mr5_T->Add(hDPhiR_mr5_Z);
  hDPhiR_mr5_W->Add(hDPhiR_mr5_Z);

  scale=hDPhiR_mr5_D->Integral()/hDPhiR_mr5_Q->Integral();
  //hDPhiR_mr5_Q->Scale(scale);
  //hDPhiR_mr5_T->Scale(scale);
  //hDPhiR_mr5_W->Scale(scale);
  //hDPhiR_mr5_Z->Scale(scale);

  TLegend *leg = new TLegend(0.67,0.57,0.89,0.86);
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetLineColor(0);
  leg->AddEntry(hMR_D, "Data", "pel");
  leg->AddEntry(hMR_Q, "QCD", "f");
  leg->AddEntry(hMR_T, "TT+jets", "f");
  leg->AddEntry(hMR_W, "W+jets", "f");
  leg->AddEntry(hMR_Z, "Zvv+jets", "f");  

  //--------------------------------------------------------------
  //
  // draw
  //
  //--------------------------------------------------------------

  hMR_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hMR_Q->GetMaximum(), hMR_D->GetMaximum()));
  hMR_Q->GetXaxis()->SetTitle("M_{R}");
  hMR_Q->GetYaxis()->SetTitle("Events");
  hMR_Q->SetTitle("");
  hMR_Q->Draw("hist");
  hMR_T->Draw("histsame");
  hMR_W->Draw("histsame");
  hMR_Z->Draw("histsame");
  hMR_D->Draw("same e");
  leg->Draw();
  c->SaveAs("MR_multijet.png");

  hMR_pass_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hMR_pass_Q->GetMaximum(), hMR_pass_D->GetMaximum()));
  hMR_pass_Q->GetXaxis()->SetTitle("M_{R}");
  hMR_pass_Q->GetYaxis()->SetTitle("Events");
  hMR_pass_Q->SetTitle("");
  hMR_pass_Q->Draw("hist");
  hMR_pass_T->Draw("histsame");
  hMR_pass_W->Draw("histsame");
  hMR_pass_Z->Draw("histsame");
  hMR_pass_D->Draw("same e");
  leg->Draw();
  c->SaveAs("MR_pass_multijet.png");

  hMR_fail_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hMR_fail_Q->GetMaximum(), hMR_fail_D->GetMaximum()));
  hMR_fail_Q->GetXaxis()->SetTitle("M_{R}");
  hMR_fail_Q->GetYaxis()->SetTitle("Events");
  hMR_fail_Q->SetTitle("");
  hMR_fail_Q->Draw("hist");
  hMR_fail_T->Draw("histsame");
  hMR_fail_W->Draw("histsame");
  hMR_fail_Z->Draw("histsame");
  hMR_fail_D->Draw("same e");
  leg->Draw();
  c->SaveAs("MR_fail_multijet.png");
  
  hRsq_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hRsq_Q->GetMaximum(), hRsq_D->GetMaximum()));
  hRsq_Q->GetXaxis()->SetTitle("R^{2}");
  hRsq_Q->GetYaxis()->SetTitle("Events");
  hRsq_Q->SetTitle("");
  hRsq_Q->Draw("hist");
  hRsq_T->Draw("histsame");
  hRsq_W->Draw("histsame");
  hRsq_Z->Draw("histsame");
  hRsq_D->Draw("same e");
  leg->Draw();
  c->SaveAs("Rsq_multijet.png");

  hRsq_pass_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hRsq_pass_Q->GetMaximum(), hRsq_pass_D->GetMaximum()));
  hRsq_pass_Q->GetXaxis()->SetTitle("R^{2}");
  hRsq_pass_Q->GetYaxis()->SetTitle("Events");
  hRsq_pass_Q->SetTitle("");
  hRsq_pass_Q->Draw("hist");
  hRsq_pass_T->Draw("histsame");
  hRsq_pass_W->Draw("histsame");
  hRsq_pass_Z->Draw("histsame");
  hRsq_pass_D->Draw("same e");
  leg->Draw();
  c->SaveAs("Rsq_pass_multijet.png");

  hRsq_fail_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hRsq_fail_Q->GetMaximum(), hRsq_fail_D->GetMaximum()));
  hRsq_fail_Q->GetXaxis()->SetTitle("R^{2}");
  hRsq_fail_Q->GetYaxis()->SetTitle("Events");
  hRsq_fail_Q->SetTitle("");
  hRsq_fail_Q->Draw("hist");
  hRsq_fail_T->Draw("histsame");
  hRsq_fail_W->Draw("histsame");
  hRsq_fail_Z->Draw("histsame");
  hRsq_fail_D->Draw("same e");
  leg->Draw();
  c->SaveAs("Rsq_fail_multijet.png");

  hLjpt_Q->GetYaxis()->SetRangeUser(0.0,1.3*TMath::Max(hLjpt_Q->GetMaximum(), hLjpt_D->GetMaximum()));
  hLjpt_Q->GetXaxis()->SetTitle("lead. jet p_{T}");
  hLjpt_Q->GetYaxis()->SetTitle("Events");
  hLjpt_Q->SetTitle("");
  hLjpt_Q->Draw("hist");
  hLjpt_T->Draw("histsame");
  hLjpt_W->Draw("histsame");
  hLjpt_Z->Draw("histsame");
  hLjpt_D->Draw("same e");
  leg->Draw();
  c->SaveAs("Ljpt_multijet.png");

  hLjpt_pass_Q->GetYaxis()->SetRangeUser(0.0,1.3*TMath::Max(hLjpt_pass_Q->GetMaximum(), hLjpt_pass_D->GetMaximum()));
  hLjpt_pass_Q->GetXaxis()->SetTitle("lead. jet p_{T}");
  hLjpt_pass_Q->GetYaxis()->SetTitle("Events");
  hLjpt_pass_Q->SetTitle("");
  hLjpt_pass_Q->Draw("hist");
  hLjpt_pass_T->Draw("histsame");
  hLjpt_pass_W->Draw("histsame");
  hLjpt_pass_Z->Draw("histsame");
  hLjpt_pass_D->Draw("same e");
  leg->Draw();
  c->SaveAs("Ljpt_pass_multijet.png");

  hLjpt_fail_Q->GetYaxis()->SetRangeUser(0.0,1.3*TMath::Max(hLjpt_fail_Q->GetMaximum(), hLjpt_fail_D->GetMaximum()));
  hLjpt_fail_Q->GetXaxis()->SetTitle("lead. jet p_{T}");
  hLjpt_fail_Q->GetYaxis()->SetTitle("Events");
  hLjpt_fail_Q->SetTitle("");
  hLjpt_fail_Q->Draw("hist");
  hLjpt_fail_T->Draw("histsame");
  hLjpt_fail_W->Draw("histsame");
  hLjpt_fail_Z->Draw("histsame");
  hLjpt_fail_D->Draw("same e");
  leg->Draw();
  c->SaveAs("Ljpt_fail_multijet.png");

  leg->SetX1NDC(0.39); leg->SetX2NDC(0.61);

  hDPhiR_Q->GetYaxis()->SetRangeUser(0.0,2.0*TMath::Max(hDPhiR_Q->GetMaximum(), hDPhiR_D->GetMaximum()));
  hDPhiR_Q->GetXaxis()->SetTitle("#Delta#phi_{R}");
  hDPhiR_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_Q->SetTitle("");

  hDPhiR_Q->Draw("hist");
  hDPhiR_T->Draw("histsame");
  hDPhiR_W->Draw("histsame");
  hDPhiR_Z->Draw("histsame");
  hDPhiR_D->Draw("same e");

  leg->Draw();

  c->SaveAs("DPhiR_multijet.png");

  hDPhiR_rsq1_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_rsq1_Q->GetMaximum(), hDPhiR_rsq1_D->GetMaximum()));
  hDPhiR_rsq1_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_rsq1_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_rsq1_Q->SetTitle("0.15 < R^2 < 0.20");
  
  hDPhiR_rsq1_Q->Draw("hist");
  hDPhiR_rsq1_T->Draw("histsame");
  hDPhiR_rsq1_W->Draw("histsame");
  hDPhiR_rsq1_Z->Draw("histsame");
  hDPhiR_rsq1_D->Draw("same e");
  
  leg->Draw();

  c->SaveAs("DPhiR_rsq1_multijet.png");

  hDPhiR_rsq2_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_rsq2_Q->GetMaximum(), hDPhiR_rsq2_D->GetMaximum()));
  hDPhiR_rsq2_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_rsq2_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_rsq2_Q->SetTitle("0.20 < R^2 < 0.25");
  
  hDPhiR_rsq2_Q->Draw("hist");
  hDPhiR_rsq2_T->Draw("histsame");
  hDPhiR_rsq2_W->Draw("histsame");
  hDPhiR_rsq2_Z->Draw("histsame");
  hDPhiR_rsq2_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_rsq2_multijet.png");

  hDPhiR_rsq3_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_rsq3_Q->GetMaximum(), hDPhiR_rsq3_D->GetMaximum()));
  hDPhiR_rsq3_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_rsq3_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_rsq3_Q->SetTitle("0.25 < R^2 < 0.30");
  
  hDPhiR_rsq3_Q->Draw("hist");
  hDPhiR_rsq3_T->Draw("histsame");
  hDPhiR_rsq3_W->Draw("histsame");
  hDPhiR_rsq3_Z->Draw("histsame");
  hDPhiR_rsq3_D->Draw("same e");

  leg->Draw();

  c->SaveAs("DPhiR_rsq3_multijet.png");
  
  hDPhiR_ljpt1_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_ljpt1_Q->GetMaximum(), hDPhiR_ljpt1_D->GetMaximum()));
  hDPhiR_ljpt1_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_ljpt1_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_ljpt1_Q->SetTitle("80 < lead. jet p_{T} < 200");
  
  hDPhiR_ljpt1_Q->Draw("hist");
  hDPhiR_ljpt1_T->Draw("histsame");
  hDPhiR_ljpt1_W->Draw("histsame");
  hDPhiR_ljpt1_Z->Draw("histsame");
  hDPhiR_ljpt1_D->Draw("same e");
  
  leg->Draw();

  c->SaveAs("DPhiR_ljpt1_multijet.png");

  hDPhiR_ljpt2_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_ljpt2_Q->GetMaximum(), hDPhiR_ljpt2_D->GetMaximum()));
  hDPhiR_ljpt2_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_ljpt2_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_ljpt2_Q->SetTitle("200 < lead. jet p_{T} < 250");
  
  hDPhiR_ljpt2_Q->Draw("hist");
  hDPhiR_ljpt2_T->Draw("histsame");
  hDPhiR_ljpt2_W->Draw("histsame");
  hDPhiR_ljpt2_Z->Draw("histsame");
  hDPhiR_ljpt2_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_ljpt2_multijet.png");
  
  hDPhiR_ljpt3_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_ljpt3_Q->GetMaximum(), hDPhiR_ljpt3_D->GetMaximum()));
  hDPhiR_ljpt3_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_ljpt3_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_ljpt3_Q->SetTitle("250 < lead. jet p_{T} < 300");
  
  hDPhiR_ljpt3_Q->Draw("hist");
  hDPhiR_ljpt3_T->Draw("histsame");
  hDPhiR_ljpt3_W->Draw("histsame");
  hDPhiR_ljpt3_Z->Draw("histsame");
  hDPhiR_ljpt3_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_ljpt3_multijet.png");
  
  hDPhiR_ljpt4_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_ljpt4_Q->GetMaximum(), hDPhiR_ljpt4_D->GetMaximum()));
  hDPhiR_ljpt4_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_ljpt4_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_ljpt4_Q->SetTitle("300 < lead. jet p_{T} < 400");
  
  hDPhiR_ljpt4_Q->Draw("hist");
  hDPhiR_ljpt4_T->Draw("histsame");
  hDPhiR_ljpt4_W->Draw("histsame");
  hDPhiR_ljpt4_Z->Draw("histsame");
  hDPhiR_ljpt4_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_ljpt4_multijet.png");
  
  hDPhiR_ljpt5_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_ljpt5_Q->GetMaximum(), hDPhiR_ljpt5_D->GetMaximum()));
  hDPhiR_ljpt5_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_ljpt5_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_ljpt5_Q->SetTitle("400 < lead. jet p_{T} < 500");
  
  hDPhiR_ljpt5_Q->Draw("hist");
  hDPhiR_ljpt5_T->Draw("histsame");
  hDPhiR_ljpt5_W->Draw("histsame");
  hDPhiR_ljpt5_Z->Draw("histsame");
  hDPhiR_ljpt5_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_ljpt5_multijet.png");

  hDPhiR_mr1_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_mr1_Q->GetMaximum(), hDPhiR_mr1_D->GetMaximum()));
  hDPhiR_mr1_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_mr1_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_mr1_Q->SetTitle("500 < MR < 600");
  
  hDPhiR_mr1_Q->Draw("hist");
  hDPhiR_mr1_T->Draw("histsame");
  hDPhiR_mr1_W->Draw("histsame");
  hDPhiR_mr1_Z->Draw("histsame");
  hDPhiR_mr1_D->Draw("same e");
  
  leg->Draw();

  c->SaveAs("DPhiR_mr1_multijet.png");

  hDPhiR_mr2_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_mr2_Q->GetMaximum(), hDPhiR_mr2_D->GetMaximum()));
  hDPhiR_mr2_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_mr2_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_mr2_Q->SetTitle("600 < MR < 700");
  
  hDPhiR_mr2_Q->Draw("hist");
  hDPhiR_mr2_T->Draw("histsame");
  hDPhiR_mr2_W->Draw("histsame");
  hDPhiR_mr2_Z->Draw("histsame");
  hDPhiR_mr2_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_mr2_multijet.png");
  
  hDPhiR_mr3_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_mr3_Q->GetMaximum(), hDPhiR_mr3_D->GetMaximum()));
  hDPhiR_mr3_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_mr3_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_mr3_Q->SetTitle("700 < MR < 800");
  
  hDPhiR_mr3_Q->Draw("hist");
  hDPhiR_mr3_T->Draw("histsame");
  hDPhiR_mr3_W->Draw("histsame");
  hDPhiR_mr3_Z->Draw("histsame");
  hDPhiR_mr3_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_mr3_multijet.png");
  
  hDPhiR_mr4_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_mr4_Q->GetMaximum(), hDPhiR_mr4_D->GetMaximum()));
  hDPhiR_mr4_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_mr4_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_mr4_Q->SetTitle("800 < MR < 900");
  
  hDPhiR_mr4_Q->Draw("hist");
  hDPhiR_mr4_T->Draw("histsame");
  hDPhiR_mr4_W->Draw("histsame");
  hDPhiR_mr4_Z->Draw("histsame");
  hDPhiR_mr4_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_mr4_multijet.png");
  
  hDPhiR_mr5_Q->GetYaxis()->SetRangeUser(0.0,1.2*TMath::Max(hDPhiR_mr5_Q->GetMaximum(), hDPhiR_mr5_D->GetMaximum()));
  hDPhiR_mr5_Q->GetXaxis()->SetTitle("#Delta#phi_{razor}");
  hDPhiR_mr5_Q->GetYaxis()->SetTitle("Events");
  hDPhiR_mr5_Q->SetTitle("900 < MR < 1000");
  
  hDPhiR_mr5_Q->Draw("hist");
  hDPhiR_mr5_T->Draw("histsame");
  hDPhiR_mr5_W->Draw("histsame");
  hDPhiR_mr5_Z->Draw("histsame");
  hDPhiR_mr5_D->Draw("same e");
  
  leg->Draw();
  
  c->SaveAs("DPhiR_mr5_multijet.png");

}
*/
