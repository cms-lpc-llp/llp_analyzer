//--------------------------------------------------------------
//
// make ratio plots for Razor-sidebands and Dijets
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

#include <CalStyleRemix.hh>

void take_ratios_Z() {

  //--------------------------------------------------------------
  //
  // setup
  //
  //--------------------------------------------------------------

  TFile *fD = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Zjets/DoubleMuon_Run2015D_Golden.root","read");
  TFile *fZ = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Zjets/DYJetsToLL_amcatnlo_2137pb_weighted.root","read");
  TFile *fW = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Zjets/VV.root","read");
  TFile *fT = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Zjets/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_2137pb_leptonic.root","read");

  TString cut_str="weight*(box==1)*puWeight*(zMass>60 && zMass<120)*passedDileptonTrigger*(NJets40>1)*(MR>300)";
  TString cut_str_dat=cut_str+"*(1)";

  enum {mr=0, ljpt, rsq};

  // for mr
  //Int_t binIn=mr;
  //const Int_t nbinx=1, nbiny=6, nbinz=2;
  //Float_t xmin=0, xmax=1;
  //Float_t ymin=300, ymax=3000;
  //Float_t zmin=0, zmax=2;
  //Float_t xbins[nbinx+1] = { xmin, xmax };
  //Float_t ybins[nbiny+1] = { ymin, 500, 1000, 1500, 2000, 2500, ymax };
  //Float_t zbins[nbinz+1] = { zmin, 1, zmax };
  //TString pname = "npf_vs_mr_zjets.png";

  // for leading jet pt
  //Int_t binIn=ljpt;
  //const Int_t nbinx=1, nbiny=11, nbinz=2;
  //Float_t xmin=40, xmax=500;
  //Float_t ymin=40, ymax=500;
  //Float_t zmin=0, zmax=2;
  //Float_t xbins[nbinx+1] = { xmin, xmax };
  //Float_t ybins[nbiny+1] = { ymin, 60, 80, 100, 150, 200, 250, 300, 350, 400, 450, ymax };
  //Float_t zbins[nbinz+1] = { zmin, 1, zmax };
  //TString pname = "npf_vs_ljpt_zjets.png";

  // for rsq
  Int_t binIn=rsq;
  const Int_t nbinx=1, nbiny=4, nbinz=2;
  Float_t xmin=0,   xmax=3000;
  Float_t ymin=0.0, ymax=0.2;
  Float_t zmin=0,   zmax=2;
  Float_t xbins[nbinx+1] = { xmin, xmax };
  Float_t ybins[nbiny+1] = { ymin, 0.05, 0.1, 0.15, ymax };
  Float_t zbins[nbinz+1] = { zmin, 1, zmax };
  TString pname = "npf_vs_rsq_zjets.png";
  
  TH3F *dPhiPF_D = new TH3F("dPhiPF_D", "dPhiPF_D", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_D->Sumw2();
  TH3F *dPhiPF_T = new TH3F("dPhiPF_T", "dPhiPF_T", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_T->Sumw2();
  TH3F *dPhiPF_W = new TH3F("dPhiPF_W", "dPhiPF_W", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_W->Sumw2();
  TH3F *dPhiPF_Z = new TH3F("dPhiPF_Z", "dPhiPF_Z", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_Z->Sumw2();

  // open trees
  TTree *tD = (TTree*) fD->Get("QCDTree");
  TTree *tT = (TTree*) fT->Get("QCDTree");
  TTree *tW = (TTree*) fW->Get("QCDTree");
  TTree *tZ = (TTree*) fZ->Get("QCDTree");

  //--------------------------------------------------------------
  //
  // hopefully no configuration below
  //
  //--------------------------------------------------------------

  TCanvas *c = MakeCanvas("c","c",800,600);

  if (binIn==mr) {
    tD->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_D", cut_str_dat);
    tT->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_T", cut_str);
    tW->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_W", cut_str);
    tZ->Draw("(abs(dPhiRazor)>2.8)+0.5:MR:Rsq>>dPhiPF_Z", cut_str);
  }
  else if (binIn==ljpt) {
    tD->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt:leadingJetPt>>dPhiPF_D", cut_str_dat);
    tT->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt:leadingJetPt>>dPhiPF_T", cut_str);
    tW->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt:leadingJetPt>>dPhiPF_W", cut_str);
    tZ->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt:leadingJetPt>>dPhiPF_Z", cut_str);
  }
  else if (binIn==rsq) {
    tD->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_D", cut_str_dat);
    tT->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_T", cut_str);
    tW->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_W", cut_str);
    tZ->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_Z", cut_str);
  }

  //dPhiPF_D->Add(dPhiPF_T, -1);
  //dPhiPF_D->Add(dPhiPF_W, -1);

  vector<TGraphAsymmErrors *> qcd_with_mr; 
  for (Int_t i=0; i<dPhiPF_Z->GetNbinsX(); i++) {
    qcd_with_mr.push_back(new TGraphAsymmErrors());
    qcd_with_mr[i]->SetMarkerStyle(20);
    qcd_with_mr[i]->SetMarkerColor(kViolet);
    qcd_with_mr[i]->SetLineColor(kViolet-1);

    Int_t k=0;
  
    for (Int_t j=0; j<dPhiPF_Z->GetNbinsY(); j++) {
  
      Float_t rsq=dPhiPF_Z->GetYaxis()->GetBinCenter(j+1);

      Float_t dnP=dPhiPF_Z->GetBinError(dPhiPF_Z->GetBin(i+1,j+1,1))/dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,1));
      Float_t dnF=dPhiPF_Z->GetBinError(dPhiPF_Z->GetBin(i+1,j+1,2))/dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,2));

      Float_t nPF=0;
      if (dnP>0 && dnF>0) {
	nPF=dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,1))/dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,2));

	Float_t drsq=0.5*dPhiPF_Z->GetYaxis()->GetBinWidth(j+1);
	Float_t dnPF_u=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
	Float_t dnPF_l=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
	dnPF_l = (dnPF_l<nPF ? dnPF_l : nPF);
  
	qcd_with_mr[i]->SetPoint(k, rsq, nPF);
	qcd_with_mr[i]->SetPointError(k, drsq, drsq, dnPF_l, dnPF_u);
	k++;
      }
    }
  }

  vector<TGraphAsymmErrors *> dat_with_mr; 
  Float_t xrsq=0, drsq=0, dnP=0, dnF=0, dnP_T=0, dnF_T=0, dnP_W=0, dnF_W=0, nPass=0, nFail=0, nPF=0;
  for (Int_t i=0; i<dPhiPF_D->GetNbinsX(); i++) {
    dat_with_mr.push_back(new TGraphAsymmErrors());
    dat_with_mr[i]->SetMarkerStyle(20);

    Int_t k=0;

    for (Int_t j=0; j<dPhiPF_D->GetNbinsY(); j++) {
      xrsq=0; drsq=0;
      dnP=0; dnF=0; 
      dnP_T=0; dnF_T=0; 
      dnP_W=0; dnF_W=0; 
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

      nPass=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,1)) - dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,1)) 
	- dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,1));

      nFail=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,2)) - dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,2)) 
	- dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,2));

      dnP+=dnP_T*dnP_T+dnP_W*dnP_W;
      dnF+=dnF_T*dnF_T+dnF_W*dnF_W;

      if (nPass>0 && nFail>0) {
	nPF=nPass/nFail;
  
	drsq=0.5*dPhiPF_D->GetYaxis()->GetBinWidth(j+1);;

	Float_t dnPF_u=nPF*TMath::Sqrt( dnP + dnF );
	Float_t dnPF_l=nPF*TMath::Sqrt( dnP + dnF );
	dnPF_l = (dnPF_l<nPF ? dnPF_l : nPF);
	
	dat_with_mr[i]->SetPoint(k, xrsq, nPF);
	dat_with_mr[i]->SetPointError(k, drsq, drsq, dnPF_l, dnPF_u);
	k++;
      }
    }
  }


  TLegend *leg = new TLegend(0.39,0.70,0.70,0.86);
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetLineColor(0);
  leg->AddEntry(dat_with_mr[0], "Data-(tt+ewk)", "pel");
  leg->AddEntry(qcd_with_mr[0], "Z+jets", "pel");

  Int_t i=0;

  qcd_with_mr[i]->GetYaxis()->SetTitle("N_{pass}/N_{fail}");
  if (binIn==mr) qcd_with_mr[i]->GetXaxis()->SetTitle("M_{R}");
  else if (binIn==ljpt) {
    qcd_with_mr[i]->GetXaxis()->SetTitle("lead. jet p_{T}");
  }
  else if (binIn==rsq) {
    qcd_with_mr[i]->GetXaxis()->SetTitle("R^{2}");
    qcd_with_mr[i]->GetYaxis()->SetRangeUser(0, 1.5);
  }
  qcd_with_mr[i]->GetXaxis()->SetNdivisions(508);
  qcd_with_mr[i]->GetYaxis()->SetNdivisions(508);
  qcd_with_mr[i]->SetTitle("");
  qcd_with_mr[i]->Draw("ap");
  dat_with_mr[i]->Draw("psame");
  
  leg->Draw();
  c->SaveAs(pname);

}

  //cout << dPhiPF_Q->GetNbinsX() << endl;
  //cout << dPhiPF_Q->GetNbinsY() << endl;
  //cout << dPhiPF_Q->GetNbinsZ() << endl;

  //TH1F* qcd_mr  = (TH1F*) dPhiPF_Q->ProjectionX("qcd_mr");
  //TH1F* qcd_rsq = (TH1F*) dPhiPF_Q->ProjectionY("qcd_rsq");
  //TH1F* qcd_pf  = (TH1F*) dPhiPF_Q->ProjectionZ("qcd_pf");
  //
  //TH1F* dat_mr  = (TH1F*) dPhiPF_D->ProjectionX("dat_mr");
  //TH1F* dat_rsq = (TH1F*) dPhiPF_D->ProjectionY("dat_rsq");
  //TH1F* dat_pf  = (TH1F*) dPhiPF_D->ProjectionZ("dat_pf");
