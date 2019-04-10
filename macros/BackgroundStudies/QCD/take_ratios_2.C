#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TAxis.h>
#include <TMath.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
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

void take_ratios_2() {

  //--------------------------------------------------------------
  //
  // setup
  //
  //--------------------------------------------------------------

  TFile *fD = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Razor/HTMHT_Run2015D_Golden.root","read");
  TFile *fQ = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Razor/QCD_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2137pb_skim.root","read");
  TFile *fT = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Razor/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_2137pb_skim.root","read");
  TFile *fW = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Razor/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2137pb_skim.root","read");
  TFile *fZ = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Razor/ZJetsToNuNu_13TeV-madgraph_2137pb_skim.root","read");

  TString cut_str="weight*(box==11||box==12)*puWeight*(Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter)";
  TString cut_str_dat="weight*(box==11||box==12)*puWeight*(Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter)";

  const Int_t nbinx=13, nbiny=2;
  Float_t xmin=80, ymin=0;
  Float_t xmax=700, ymax=2;
  Float_t xbins[nbinx+1] = { xmin, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 600, xmax };
  Float_t ybins[nbiny+1] = { ymin, 1, ymax };

  TH2F *dPhiPF_Q = new TH2F("dPhiPF_Q", "dPhiPF_Q", nbinx, &xbins[0], nbiny, &ybins[0]); dPhiPF_Q->Sumw2();
  TH2F *dPhiPF_D = new TH2F("dPhiPF_D", "dPhiPF_D", nbinx, &xbins[0], nbiny, &ybins[0]); dPhiPF_D->Sumw2();
  TH2F *dPhiPF_T = new TH2F("dPhiPF_T", "dPhiPF_T", nbinx, &xbins[0], nbiny, &ybins[0]); dPhiPF_T->Sumw2();
  TH2F *dPhiPF_W = new TH2F("dPhiPF_W", "dPhiPF_W", nbinx, &xbins[0], nbiny, &ybins[0]); dPhiPF_W->Sumw2();
  TH2F *dPhiPF_Z = new TH2F("dPhiPF_Z", "dPhiPF_Z", nbinx, &xbins[0], nbiny, &ybins[0]); dPhiPF_Z->Sumw2();

  // open trees
  TTree *tQ = (TTree*) fQ->Get("QCDTree");
  TTree *tD = (TTree*) fD->Get("QCDTree");
  TTree *tT = (TTree*) fT->Get("QCDTree");
  TTree *tW = (TTree*) fW->Get("QCDTree");
  TTree *tZ = (TTree*) fZ->Get("QCDTree");

  //--------------------------------------------------------------
  //
  // hopefully no configuration below
  //
  //--------------------------------------------------------------

  TCanvas *c = new TCanvas("c","c",800,600);
  gStyle->SetOptStat(0);

  // draw Rsq (again)
  tQ->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt>>dPhiPF_Q", cut_str);
  tD->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt>>dPhiPF_D", cut_str);
  tT->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt>>dPhiPF_T", cut_str);
  tW->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt>>dPhiPF_W", cut_str);
  tZ->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt>>dPhiPF_Z", cut_str);

  dPhiPF_D->Add(dPhiPF_T, -1);
  dPhiPF_D->Add(dPhiPF_W, -1);
  dPhiPF_D->Add(dPhiPF_Z, -1);

  TGraphAsymmErrors *qcd_with_mr = new TGraphAsymmErrors();
  qcd_with_mr->SetMarkerStyle(20);
  qcd_with_mr->SetMarkerColor(kViolet);
  qcd_with_mr->SetLineColor(kViolet-1);
  
  for (Int_t j=0; j<dPhiPF_Q->GetNbinsX(); j++) {
    Float_t rsq=dPhiPF_Q->GetXaxis()->GetBinCenter(j+1);
    Float_t dnP=dPhiPF_Q->GetBinError(dPhiPF_Q->GetBin(j+1,1))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(j+1,1));
    Float_t dnF=dPhiPF_Q->GetBinError(dPhiPF_Q->GetBin(j+1,2))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(j+1,2));
    Float_t nPF=0;
    if (dnP>0 && dnF>0) nPF=dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(j+1,1))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(j+1,2));
    Float_t drsq=0.005;
    Float_t dnPF_u=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
    Float_t dnPF_l=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
    dnPF_l = (dnPF_l<nPF ? dnPF_l : nPF);
    qcd_with_mr->SetPoint(j, rsq, nPF);
    qcd_with_mr->SetPointError(j, drsq, drsq, dnPF_l, dnPF_u);
  }
  
  TGraphAsymmErrors *dat_with_mr = new TGraphAsymmErrors();
  dat_with_mr->SetMarkerStyle(20);
  
  for (Int_t j=0; j<dPhiPF_D->GetNbinsX(); j++) { 
    Float_t rsq=dPhiPF_D->GetXaxis()->GetBinCenter(j+1);
    Float_t dnP=dPhiPF_D->GetBinError(dPhiPF_D->GetBin(j+1,1))/dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(j+1,1));
    Float_t dnF=dPhiPF_D->GetBinError(dPhiPF_D->GetBin(j+1,2))/dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(j+1,2));
    Float_t nPF=0;
    if (dnP>0 && dnF>0) nPF=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(j+1,1))/dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(j+1,2));
    Float_t drsq=0.005;
    Float_t dnPF_u=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
    Float_t dnPF_l=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
    dnPF_l = (dnPF_l<nPF ? dnPF_l : nPF);
    dat_with_mr->SetPoint(j, rsq, nPF);
    dat_with_mr->SetPointError(j, drsq, drsq, dnPF_l, dnPF_u);
  }

  TF1 *qcd_fxn = new TF1("qcd_fxn","[0]*x^[1]+[2]",80,3000);
  qcd_fxn->SetLineColor(kViolet);
  TF1 *dat_fxn = new TF1("dat_fxn","[0]*x^[1]+[2]",80,3000);
  dat_fxn->SetParameter(0,100);
  dat_fxn->SetParameter(1,-1);
  dat_fxn->SetParameter(2,0.1);
  dat_fxn->SetLineColor(kBlack);
  
  qcd_with_mr->GetXaxis()->SetTitle("lead. jet p_{T}");
  qcd_with_mr->GetYaxis()->SetTitle("N_{pass}/N_{fail}");
  qcd_with_mr->SetTitle("");
  qcd_with_mr->Draw("ap");
  qcd_with_mr->Fit(qcd_fxn,"R");
  dat_with_mr->Fit(dat_fxn,"R");
  dat_with_mr->Draw("psame");
  c->SaveAs("npf_vs_ljpt.png");

}
