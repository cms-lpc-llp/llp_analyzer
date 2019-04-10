//--------------------------------------------------------------
//
// make ratio plots for Razor sideband and dijet control region
//
//--------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TAxis.h>
#include <TMath.h>
#include <TH3F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TBenchmark.h>                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#endif

#include <CalStyleRemix.hh>

void makepaperplots_dijet() {

  enum {mr=0, ljpt, rsq};
  enum {razor=0, dijet};
  //--------------------------------------------------------------
  //
  // RAZOR-TRIGGERED region setup
  //
  //--------------------------------------------------------------

  // for mr

  //Int_t binIn=mr;
  //const Int_t nbinx=1, nbiny=4, nbinz=2;
  //Float_t xmin=0.15, xmax=0.25; 
  //Float_t ymin=400,  ymax=3000;
  //Float_t zmin=0,    zmax=2;
  //Float_t xbins[nbinx+1] = {xmin, xmax};
  //Float_t ybins[nbiny+1] = {ymin, 600, 800, 1000, ymax};
  //Float_t zbins[nbinz+1] = {zmin, 1, zmax};
  //
  //TString pname = "npf_vs_mr_dijet_2b_fit.pdf";
  //TString pname2 = "npf_vs_mr_dijet_2b_fit.C";

  // for rsq
  Int_t binIn=rsq;
  //const Int_t nbinx=1, nbiny=10, nbinz=2;
  const Int_t nbinx=1, nbiny=8, nbinz=2;
  Float_t xmin=400, xmax=3000; 
  Float_t ymin=0.15,  ymax=0.35; 
  Float_t zmin=0,    zmax=2;
  Float_t xbins[nbinx+1] = {xmin, xmax};
  //Float_t ybins[nbiny+1] = { ymin, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, ymax};
  Float_t ybins[nbiny+1] = { ymin, 0.16, 0.17, 0.18, 0.19, 0.20, 0.225, 0.25, ymax};
  Float_t zbins[nbinz+1] = {zmin, 1, zmax};
  TString pname = "npf_vs_rsq_dijet_0b.pdf";
  TString pname2 = "npf_vs_rsq_dijet_0b.C";

  TFile *fD = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p13_05Mar2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root","read");
  TFile *fQ = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_03Feb2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_QCD_1pb_weighted.root","read");
  TFile *fT = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_03Feb2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_TTJetsInclusive_1pb_weighted.root","read");
  TFile *fW = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_03Feb2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_WJets_1pb_weighted.root","read");
  TFile *fZ = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_03Feb2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_ZInv_1pb_weighted.root","read");

  //TFile *fD = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p3_2015JECs/FullRazorInclusive_Data_NoDuplicates_GoodLumiGolden.root","read");
  //TFile *fQ = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p2_JEC2015V6/FullRazorInclusive_QCD_1pb_weighted.root","read");
  //TFile *fT = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p2_JEC2015V6/FullRazorInclusive_TTJets_1pb_weighted.root","read");
  //TFile *fW = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p2_JEC2015V6/FullRazorInclusive_WJets_1pb_weighted.root","read");
  //TFile *fZ = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p2_JEC2015V6/FullRazorInclusive_ZInv_1pb_weighted.root","read");

    TString cut_str="(box==14)*(MR>400 && Rsq>0.15)*(Flag_HBHENoiseFilter && Flag_HBHEIsoNoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter && Flag_badChargedCandidateFilter && Flag_badMuonFilter)*(nBTaggedJets==0)";
    TString cut_str_dat="(box==14)*(MR>400 && Rsq>0.15)*(Flag_HBHENoiseFilter && Flag_HBHEIsoNoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter && Flag_badChargedCandidateFilter && Flag_badMuonFilter)*(Rsq<0.25)*(nBTaggedJets==0)";
    TString cut_str_mc="*weight*35800";



  //--------------------------------------------------------------
  //
  // hopefully no configuration below
  //
  //--------------------------------------------------------------

  TH3F *dPhiPF_Q = new TH3F("dPhiPF_Q", "dPhiPF_Q", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_Q->Sumw2();
  TH3F *dPhiPF_D = new TH3F("dPhiPF_D", "dPhiPF_D", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_D->Sumw2();
  TH3F *dPhiPF_T = new TH3F("dPhiPF_T", "dPhiPF_T", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_T->Sumw2();
  TH3F *dPhiPF_W = new TH3F("dPhiPF_W", "dPhiPF_W", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_W->Sumw2();
  TH3F *dPhiPF_Z = new TH3F("dPhiPF_Z", "dPhiPF_Z", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_Z->Sumw2();

  //Float_t ybins[nbiny+1] = {ymin, 500, 600, 700, 800, 1000, ymax};
  Float_t ybins2[nbiny+3] = {float(ymin - 0.1), ymin, 600, 800, 1000, ymax, float(ymax + 0.1)};
  TH1F *fxn_plus_err = new TH1F("fxn_plus_err","fxn_plus_err", nbiny+2, &ybins2[0]); fxn_plus_err->Sumw2();  

  // open trees
  TTree *tQ = (TTree*) fQ->Get("RazorInclusive");
  TTree *tD = (TTree*) fD->Get("RazorInclusive");
  TTree *tT = (TTree*) fT->Get("RazorInclusive");
  TTree *tW = (TTree*) fW->Get("RazorInclusive");
  TTree *tZ = (TTree*) fZ->Get("RazorInclusive");

  TCanvas *c = MakeCanvas("c","c",800,600);

  if (binIn==mr) {
    tQ->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_Q", cut_str+cut_str_mc);
    tD->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_D", cut_str_dat);
    tT->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_T", cut_str+cut_str_mc);
    tW->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_W", cut_str+cut_str_mc);
    tZ->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_Z", cut_str+cut_str_mc);
  } else if (binIn==rsq) {
    tQ->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_Q", cut_str+cut_str_mc);
    tD->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_D", cut_str_dat);
    tT->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_T", cut_str+cut_str_mc);
    tW->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_W", cut_str+cut_str_mc);
    tZ->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_Z", cut_str+cut_str_mc);
  }

  TF1 *qcd_fxn = new TF1("qcd_fxn","[0]", 400, 2500);
  //qcd_fxn->SetParameter(0,0.043);
  //qcd_fxn->SetParameter(0,0.048);
  qcd_fxn->SetParameter(0,0.068);

  //qcd_fxn->SetParameter(0,3.1e7);
  //qcd_fxn->SetParameter(1,-3.1);
  //qcd_fxn->SetParameter(2,0.062);

  fxn_plus_err->SetBinContent(1,qcd_fxn->Eval(ymin-25));
  fxn_plus_err->SetBinError(1,qcd_fxn->Eval(ymin-25)*0.8);

  fxn_plus_err->SetBinContent(nbiny+2,qcd_fxn->Eval(ymax+25));
  fxn_plus_err->SetBinError(nbiny+2,qcd_fxn->Eval(ymax+25)*0.8);

  Double_t wtf=0, lesswtf=0;

  vector<TGraphAsymmErrors *> qcd_with_mr; 
  for (Int_t i=0; i<dPhiPF_Q->GetNbinsX(); i++) {
    qcd_with_mr.push_back(new TGraphAsymmErrors());

    Int_t k=0;
  
    for (Int_t j=0; j<dPhiPF_Q->GetNbinsY(); j++) {
  
      Float_t rsq=dPhiPF_Q->GetYaxis()->GetBinCenter(j+1);
  
      Float_t dnP=dPhiPF_Q->GetBinError(dPhiPF_Q->GetBin(i+1,j+1,1))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1));
      Float_t dnF=dPhiPF_Q->GetBinError(dPhiPF_Q->GetBin(i+1,j+1,2))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2));
  
      Float_t nPF=0;
      if (dnP>0 && dnF>0) {
	nPF=dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2));
	wtf+=qcd_fxn->Eval(rsq)*(dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1))+dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2)));
	lesswtf+=(dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1))+dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2)));

	Float_t drsq=0.5*dPhiPF_Q->GetYaxis()->GetBinWidth(j+1);
	Float_t dnPF_u=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
	Float_t dnPF_l=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
	dnPF_l = (dnPF_l<nPF ? dnPF_l : nPF);

	fxn_plus_err->SetBinContent(k+2,qcd_fxn->Eval(rsq));
	fxn_plus_err->SetBinError(k+2,qcd_fxn->Eval(rsq)*0.85);
  
	qcd_with_mr[i]->SetPoint(k, rsq, nPF);
	qcd_with_mr[i]->SetPointError(k, drsq, drsq, dnPF_l, dnPF_u);
	k++;

      }
    }
  }

  vector<TGraphAsymmErrors *> dat_with_mr; 
  Float_t xrsq=0, drsq=0, dnP=0, dnF=0, dnP_T=0, dnF_T=0, dnP_W=0, dnF_W=0, dnP_Z=0, dnF_Z=0, nPass=0, nFail=0, nPF=0;
  for (Int_t i=0; i<dPhiPF_D->GetNbinsX(); i++) {
    dat_with_mr.push_back(new TGraphAsymmErrors());

    Int_t k=0;

    for (Int_t j=0; j<dPhiPF_D->GetNbinsY(); j++) {
      xrsq=0; drsq=0;
      dnP=0; dnF=0; 
      dnP_T=0; dnF_T=0; 
      dnP_W=0; dnF_W=0; 
      dnP_Z=0; dnF_Z=0; 
      nPass=0; nFail=0; nPF=0;

      xrsq=dPhiPF_D->GetYaxis()->GetBinCenter(j+1);

      dnP=dPhiPF_D->GetBinError(dPhiPF_D->GetBin(i+1,j+1,1))/dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,1));
      dnF=dPhiPF_D->GetBinError(dPhiPF_D->GetBin(i+1,j+1,2))/dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,2));

      dnP=dnP*dnP;
      dnF=dnF*dnF;

      dnP_T=dPhiPF_T->GetBinError(dPhiPF_T->GetBin(i+1,j+1,1))/dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,1));
      dnF_T=dPhiPF_T->GetBinError(dPhiPF_T->GetBin(i+1,j+1,2))/dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,2));

      dnP_W=dPhiPF_W->GetBinError(dPhiPF_W->GetBin(i+1,j+1,1))/dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,1));
      dnF_W=dPhiPF_W->GetBinError(dPhiPF_W->GetBin(i+1,j+1,2))/dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,2));

      dnP_Z=dPhiPF_Z->GetBinError(dPhiPF_Z->GetBin(i+1,j+1,1))/dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,1));
      dnF_Z=dPhiPF_Z->GetBinError(dPhiPF_Z->GetBin(i+1,j+1,2))/dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,2));

      nPass=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,1)) - dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,1)) 
	- dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,1)) - dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,1));
      
      nFail=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,2)) - dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,2)) 
	- dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,2)) - dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,2));
      
      dnP+=dnP_T*dnP_T+dnP_W*dnP_W+dnP_Z*dnP_Z;
      dnF+=dnF_T*dnF_T+dnF_W*dnF_W+dnP_Z*dnP_Z;

      cout << nPass << ", " << nFail << endl;

      if (nPass>0 && nFail>0) {
	nPF=nPass/nFail;

	drsq=0.5*dPhiPF_D->GetYaxis()->GetBinWidth(j+1);
	Float_t dnPF_u=nPF*TMath::Sqrt( dnP + dnF );
	Float_t dnPF_l=nPF*TMath::Sqrt( dnP + dnF );
	dnPF_l = (dnPF_l<nPF ? dnPF_l : nPF);

	// if (k==0) {
	//   dat_with_mr[i]->SetPoint(k, xrsq-100, nPF);
	//   dat_with_mr[i]->SetPointError(k, drsq, drsq, dnPF_l, dnPF_u);
	//   k++;
	//   dat_with_mr[i]->SetPoint(k, xrsq, nPF);
	//   dat_with_mr[i]->SetPointError(k, drsq, drsq, dnPF_l, dnPF_u);
	//   k++;
	// } else {
	  dat_with_mr[i]->SetPoint(k, xrsq, nPF);
	  dat_with_mr[i]->SetPointError(k, drsq, drsq, dnPF_l, dnPF_u);
	  k++;
	// }


      }
    }
  }


  Int_t i=0;

  qcd_fxn->SetLineColor(kBlue+1);

  //dat_with_mr[i]->Fit("qcd_fxn");

  qcd_with_mr[i]->SetMarkerStyle(24);
  //qcd_with_mr[i]->SetMarkerStyle(20);
  //qcd_with_mr[i]->SetLineColor(kRed);
  //qcd_with_mr[i]->SetMarkerColor(kRed);

  if (binIn==mr) qcd_with_mr[i]->GetXaxis()->SetTitle("M_{R} [GeV]");  
  else if (binIn==rsq) qcd_with_mr[i]->GetXaxis()->SetTitle("R^{2}");
  qcd_with_mr[i]->GetYaxis()->SetTitle("Translation Factor #zeta");
  qcd_with_mr[i]->GetYaxis()->SetTitleOffset(1.2);
  qcd_with_mr[i]->GetYaxis()->SetRangeUser(0, 1.0);
  if (binIn==mr) qcd_with_mr[i]->GetYaxis()->SetRangeUser(0, 0.75);
  qcd_with_mr[i]->GetXaxis()->SetNdivisions(508);
  qcd_with_mr[i]->GetYaxis()->SetNdivisions(508);
  qcd_with_mr[i]->SetTitle("");
  qcd_with_mr[i]->Draw("ap e1");
  dat_with_mr[i]->Draw("p e1 same");

  qcd_with_mr[i]->GetXaxis()->SetRangeUser(399.9,3000.1);

  fxn_plus_err->SetFillColor(kAzure+7);
  fxn_plus_err->SetFillStyle(3254);
  fxn_plus_err->SetMarkerStyle(0);
  fxn_plus_err->SetLineColor(kBlue+1);
  fxn_plus_err->SetLineWidth(2);

  fxn_plus_err->Draw("same f e3");
  qcd_with_mr[i]->Draw("same p e1");
  dat_with_mr[i]->Draw("same p e1");
  qcd_fxn->Draw("same l");


  TLegend *leg = new TLegend(0.40,0.70,0.75,0.86);
  leg->SetTextSize(0.035); leg->SetTextFont(42);
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetLineColor(0);
  
  leg->AddEntry(dat_with_mr[0], "Data Control Region", "pel");
  leg->AddEntry(qcd_with_mr[0], "QCD MC Simulation", "pel");
  leg->AddEntry(fxn_plus_err, "Functional Form Model", "lf");

  //qcd_fxn->SetParameter(0,9.1e12);
  //qcd_fxn->SetParameter(1,-5.0);
  //qcd_fxn->SetParameter(2,0.048);

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.035);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.485, 0.66, "#zeta = 0.043");
  //tex->DrawLatex(0.485, 0.66, "#zeta = 0.048");
  //tex->DrawLatex(0.485, 0.66, "#zeta = 0.068");

  leg->Draw();

  tex = new TLatex(0.18,0.93,"CMS");
  tex->SetNDC();
  tex->SetTextFont(62);   tex->SetTextSize(0.065);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.71,0.93,"35.8 fb^{-1} (13 TeV)");
  tex->SetNDC();
  tex->SetTextFont(42);   tex->SetTextSize(0.05);
  tex->SetLineWidth(2);
  tex->Draw();
  c->SaveAs(pname);
  c->SaveAs(pname2);

}
