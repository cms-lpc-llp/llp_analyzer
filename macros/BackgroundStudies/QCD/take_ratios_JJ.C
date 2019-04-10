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
#include <TBenchmark.h>                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#endif

#include <CalStyleRemix.hh>

void take_ratios() {

  enum {mr=0, ljpt, rsq};
  enum {razor=0, dijet};
  //--------------------------------------------------------------
  //
  // RAZOR-TRIGGERED region setup
  //
  //--------------------------------------------------------------
  Int_t sample=razor;
  // for mr
  Int_t binIn=mr;
  //const Int_t nbinx=1, nbiny=20, nbinz=2;
  const Int_t nbinx=1, nbiny=7, nbinz=2;
  //const Int_t nbinx=1, nbiny=1, nbinz=2;
  //const Int_t nbinx=1, nbiny=7, nbinz=2;
  Float_t xmin=0.15, xmax=0.25; 
  Float_t ymin=500,  ymax=4000; 
  //Float_t ymin=400,  ymax=1000; 
  Float_t zmin=0,    zmax=2;
  Float_t xbins[nbinx+1] = {xmin, xmax};
  Float_t ybins[nbiny+1] = {ymin, 575, 650, 750, 900, 1000, 1500, ymax};
  //Float_t ybins[nbiny+1] = {ymin, ymax};
  //Float_t ybins[nbiny+1] = {ymin, 450, 500, 550, 600, 700, 800, ymax};
  Float_t zbins[nbinz+1] = {zmin, 1, zmax};
  //TString pname = "npf_vs_mr_multijet.png";
  TString pname = "npf_vs_mr_dijet.png";

  // for leading jet pt
  //Int_t binIn=ljpt;
  //const Int_t nbinx=1, nbiny=15, nbinz=2;
  //Float_t xmin=40, xmax=500;
  //Float_t ymin=40, ymax=500;
  //Float_t zmin=0, zmax=2;
  //Float_t xbins[nbinx+1] = { xmin, xmax };
  //Float_t ybins[nbiny+1] = { ymin, 60, 80, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, ymax };
  //Float_t zbins[nbinz+1] = { zmin, 1, zmax };
  //TString pname = "npf_vs_ljpt_razor.png";

  // for rsq
  //Int_t binIn=rsq;
  ////const Int_t nbinx=1, nbiny=10, nbinz=2;
  //const Int_t nbinx=1, nbiny=7, nbinz=2;
  //Float_t xmin=400, xmax=3000; 
  //Float_t ymin=0.15,  ymax=0.35; 
  //Float_t zmin=0,    zmax=2;
  //Float_t xbins[nbinx+1] = {xmin, xmax};
  ////Float_t ybins[nbiny+1] = { ymin, 0.17, 0.19, 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, ymax};
  //Float_t ybins[nbiny+1] = { ymin, 0.16, 0.17, 0.18, 0.19, 0.20, 0.25, ymax};
  //Float_t zbins[nbinz+1] = {zmin, 1, zmax};
  //TString pname = "npf_vs_rsq_razor.png";

  TFile *fD = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_02Jan2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_HTMHT_2016G_23Sep2016_SUSYUnblind.root","read");
  TFile *fQ = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_02Jan2017/Signal/FullRazorInclusive_Razor2016G_SUSYUnblind_80X_QCD_1pb_weighted.root","read");
  TFile *fT = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_02Jan2017/Signal/FullRazorInclusive_Razor2016G_SUSYUnblind_80X_TTJetsInclusive_1pb_weighted.root","read");
  TFile *fW = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_02Jan2017/Signal/FullRazorInclusive_Razor2016G_SUSYUnblind_80X_WJets_1pb_weighted.root","read");
  TFile *fZ = TFile::Open("root://eoscms//store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p8_02Jan2017/Signal/FullRazorInclusive_Razor2016G_SUSYUnblind_80X_ZInv_1pb_weighted.root","read");

  //TString cut_str="(box==11||box==12)*(MR>400 && Rsq>0.15)*(Flag_HBHENoiseFilter && Flag_HBHEIsoNoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter && Flag_badChargedCandidateFilter && Flag_badMuonFilter)*(Rsq<0.25)";

  TString cut_str="(box==14)*(MR>400 && Rsq>0.15)*(Flag_HBHENoiseFilter && Flag_HBHEIsoNoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter && Flag_badChargedCandidateFilter && Flag_badMuonFilter)*(Rsq<0.25)";

  TString cut_str_dat=cut_str;
  TString cut_str_mc="*weight*4400";

  //--------------------------------------------------------------
  //
  // DIJET region setup
  //
  //--------------------------------------------------------------
  //Int_t sample=dijet;
  // for mr
  //Int_t binIn=mr;
  //const Int_t nbinx=1, nbiny=8, nbinz=2;
  //Float_t xmin=0.0, xmax=10; 
  //Float_t ymin=0,   ymax=3000; 
  //Float_t zmin=0,   zmax=2;
  //Float_t xbins[nbinx+1] = {xmin, xmax};
  //Float_t ybins[nbiny+1] = {ymin, 400, 500, 600, 700, 800, 1000, 1500, ymax};
  //Float_t zbins[nbinz+1] = {zmin, 1, zmax};
  //TString pname = "npf_vs_mr_dijet.png";

  // for leading jet pt
  //Int_t binIn=ljpt;
  //const Int_t nbinx=1, nbiny=13, nbinz=2;
  //Float_t xmin=80, xmax=700;
  //Float_t ymin=80, ymax=700;
  //Float_t zmin=0, zmax=2;
  //Float_t xbins[nbinx+1] = { xmin, xmax };
  //Float_t ybins[nbiny+1] = { ymin, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 600, ymax };
  //Float_t zbins[nbinz+1] = { zmin, 1, zmax };
  //TString pname = "npf_vs_ljpt_dijet.png";

  // for rsq
  //Int_t binIn=rsq;
  //const Int_t nbinx=1, nbiny=4, nbinz=2;
  //Float_t xmin=400, xmax=3000; 
  //Float_t ymin=0,  ymax=0.2; 
  //Float_t zmin=0,    zmax=2;
  //Float_t xbins[nbinx+1] = {xmin, xmax};
  //Float_t ybins[nbiny+1] = { ymin, 0.05, 0.1, 0.15, ymax};
  //Float_t zbins[nbinz+1] = {zmin, 1, zmax};
  //TString pname = "npf_vs_rsq_dijet.png";
  //
  //TFile *fD = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Dijet/JetHT_Run2015D_PRv4_Golden_skim.root","read");
  //TFile *fQ = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Dijet/QCD_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2137pb_skim.root","read");
  //TFile *fT = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Dijet/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_2137pb_skim.root","read");
  //TFile *fW = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Dijet/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2137pb_skim.root","read");
  //TFile *fZ = TFile::Open("root://eoscms//store/user/jlawhorn/RazorQCD_Dijet/ZJetsToNuNu_13TeV-madgraph_2137pb_skim.root","read");
  //
  //TString cut_str="weight*(box==100)*puWeight*(Flag_HBHENoiseFilter && Flag_goodVertices && Flag_eeBadScFilter)*(MR>300)";
  //TString cut_str_dat=cut_str+"*(HLTDecision[105]*HLTPrescale[105])";

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

  TH1F *fxn_plus_err = new TH1F("fxn_plus_err","fxn_plus_err", nbiny, &ybins[0]); fxn_plus_err->Sumw2();

  // open trees
  //TTree *tQ = (TTree*) fQ->Get("QCDTree");
  //TTree *tD = (TTree*) fD->Get("QCDTree");
  //TTree *tT = (TTree*) fT->Get("QCDTree");
  //TTree *tW = (TTree*) fW->Get("QCDTree");
  //TTree *tZ = (TTree*) fZ->Get("QCDTree");

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
  } else if (binIn==ljpt) {
    tQ->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt:leadingJetPt>>dPhiPF_Q", cut_str+cut_str_mc);
    tD->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt:leadingJetPt>>dPhiPF_D", cut_str_dat);
    tT->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt:leadingJetPt>>dPhiPF_T", cut_str+cut_str_mc);
    tW->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt:leadingJetPt>>dPhiPF_W", cut_str+cut_str_mc);
    tZ->Draw("(abs(dPhiRazor)>2.8)+0.5:leadingJetPt:leadingJetPt>>dPhiPF_Z", cut_str+cut_str_mc);
  } else if (binIn==rsq) {
    tQ->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_Q", cut_str+cut_str_mc);
    tD->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_D", cut_str_dat);
    tT->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_T", cut_str+cut_str_mc);
    tW->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_W", cut_str+cut_str_mc);
    tZ->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_Z", cut_str+cut_str_mc);
  }


  //9.10981e+12   4.90629e+13   4.34389e+06  -6.75286e-17
  // 2  p1          -5.00116e+00   8.35076e-01   6.67793e-05  -3.23865e-03
  // 3  p2           4.79514e-02 

  //TF1 *qcd_fxn = new TF1("qcd_fxn","[0]*x^[1]+[2]", 400, 2500);
  TF1 *qcd_fxn = new TF1("qcd_fxn","[0]", 500, 4000);
  //TF1 *qcd_fxn = new TF1("qcd_fxn","[0]", 0, 2500);
  //qcd_fxn->SetParameter(0,0.242121);

  //qcd_fxn->SetParameter(0,0.103);
  qcd_fxn->SetParameter(0,0.0500);

  //qcd_fxn->SetParameter(0,9.1e12);
  //qcd_fxn->SetParameter(1,-5.0);
  //qcd_fxn->SetParameter(2,0.048);

  Double_t wtf=0, lesswtf=0;

  //cout << qcd_fxn->Integral(400,2500) << endl;
  //cout << "----" << endl;
  //TF1 *qcd_fxn_2 = new TF1("qcd_fxn_2","[0]*x^[1]+[2]", 400, 2500);
  //qcd_fxn_2->SetParameter(0,3.1e7*1.67);
  //qcd_fxn_2->SetParameter(1,-3.1);
  //qcd_fxn_2->SetParameter(2,0.062*1.67);

  //TF1 *qcd_fxn = new TF1("qcd_fxn","[0]",0.15,0.35);
  //qcd_fxn->SetParameter(0,0.19);

  vector<TGraphAsymmErrors *> qcd_with_mr; 
  for (Int_t i=0; i<dPhiPF_Q->GetNbinsX(); i++) {
    qcd_with_mr.push_back(new TGraphAsymmErrors());
    qcd_with_mr[i]->SetMarkerStyle(20);
    qcd_with_mr[i]->SetMarkerColor(kViolet);
    qcd_with_mr[i]->SetLineColor(kViolet-1);

    Int_t k=0;
  
    for (Int_t j=0; j<dPhiPF_Q->GetNbinsY(); j++) {
  
      Float_t rsq=dPhiPF_Q->GetYaxis()->GetBinCenter(j+1);
  
      Float_t dnP=dPhiPF_Q->GetBinError(dPhiPF_Q->GetBin(i+1,j+1,1))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1));
      Float_t dnF=dPhiPF_Q->GetBinError(dPhiPF_Q->GetBin(i+1,j+1,2))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2));
  
      Float_t nPF=0;
      if (dnP>0 && dnF>0) {
	nPF=dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2));
	if (nPF<0.4) {
	  wtf+=nPF;
	  lesswtf+=1;
	}
	//wtf+=qcd_fxn->Eval(rsq)*(dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1))+dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2)));
	//lesswtf+=(dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1))+dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2)));
	cout << (nPF - qcd_fxn->Eval(rsq))/qcd_fxn->Eval(rsq) << endl;
	
	Float_t drsq=0.5*dPhiPF_Q->GetYaxis()->GetBinWidth(j+1);
	Float_t dnPF_u=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
	Float_t dnPF_l=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
	dnPF_l = (dnPF_l<nPF ? dnPF_l : nPF);

	fxn_plus_err->SetBinContent(k+1,qcd_fxn->Eval(rsq));
	//fxn_plus_err->SetBinError(k+1,qcd_fxn->Eval(rsq)*0.8);
	fxn_plus_err->SetBinError(k+1,qcd_fxn->Eval(rsq));
	//fxn_plus_err->SetBinContent(k+1,0.192443);
	//fxn_plus_err->SetBinError(k+1,0.192443*0.87);
  
	qcd_with_mr[i]->SetPoint(k, rsq, nPF);
	qcd_with_mr[i]->SetPointError(k, drsq, drsq, dnPF_l, dnPF_u);
	k++;
      }
    }
    //cout << qcd_with_mr[i]->GetMean() << endl;
  }

  cout << "----" << endl;
  cout << wtf  << " / " << lesswtf << " = " << wtf/lesswtf << endl;
  cout << "----" << endl;
  wtf=0; lesswtf=0;
  vector<TGraphAsymmErrors *> dat_with_mr; 
  Float_t xrsq=0, drsq=0, dnP=0, dnF=0, dnP_T=0, dnF_T=0, dnP_W=0, dnF_W=0, dnP_Z=0, dnF_Z=0, nPass=0, nFail=0, nPF=0;
  for (Int_t i=0; i<dPhiPF_D->GetNbinsX(); i++) {
    dat_with_mr.push_back(new TGraphAsymmErrors());
    dat_with_mr[i]->SetMarkerStyle(20);

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

      if (sample==razor) {
	nPass=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,1)) - dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,1)) 
	  - dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,1)) - dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,1));
	
	nFail=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,2)) - dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,2)) 
	  - dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,2)) - dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,2));
	
	dnP+=dnP_T*dnP_T+dnP_W*dnP_W+dnP_Z*dnP_Z;
	dnF+=dnF_T*dnF_T+dnF_W*dnF_W+dnP_Z*dnP_Z;
      }
      else if (sample==dijet) {
	nPass=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,1));
	nFail=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,2));
      }
      //cout << nPass << " / " << nFail << " = " << nPass/nFail << endl;
      if (nPass>0 && nFail>0) {
	nPF=nPass/nFail;
	wtf+=nPF;
	lesswtf+=1;
	//cout << (nPF - qcd_fxn->Eval(xrsq))/qcd_fxn->Eval(xrsq) << endl;

	Float_t drsq=0.5*dPhiPF_D->GetYaxis()->GetBinWidth(j+1);
	Float_t dnPF_u=nPF*TMath::Sqrt( dnP + dnF );
	Float_t dnPF_l=nPF*TMath::Sqrt( dnP + dnF );
	dnPF_l = (dnPF_l<nPF ? dnPF_l : nPF);
	
	dat_with_mr[i]->SetPoint(k, xrsq, nPF);
	dat_with_mr[i]->SetPointError(k, drsq, drsq, dnPF_l, dnPF_u);
	k++;
      }
    }
  }

  //cout << "----" << endl;
  cout << wtf  << " / " << lesswtf << " = " << wtf/lesswtf << endl;
  cout << "----" << endl;

  TLegend *leg = new TLegend(0.50,0.70,0.75,0.86);
  //TLegend *leg = new TLegend(0.25,0.70,0.50,0.86);
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetLineColor(0);
  if (sample==razor) leg->AddEntry(dat_with_mr[0], "Data-(tt+ewk)", "pel");
  else leg->AddEntry(dat_with_mr[0], "Data", "pel");
  leg->AddEntry(qcd_with_mr[0], "QCD", "pel");
  leg->AddEntry(qcd_fxn, "Fit result", "l");

  Int_t i=0;
  
  if (binIn==mr) qcd_with_mr[i]->GetXaxis()->SetTitle("M_{R}");  
  else if (binIn==ljpt) qcd_with_mr[i]->GetXaxis()->SetTitle("lead. jet p_{T}");
  else if (binIn==rsq) qcd_with_mr[i]->GetXaxis()->SetTitle("R^{2}");
  qcd_with_mr[i]->GetYaxis()->SetTitle("N_{pass}/N_{fail}");
  qcd_with_mr[i]->GetYaxis()->SetRangeUser(0, 1.0);
  //if (sample==razor && binIn==mr) qcd_with_mr[i]->GetYaxis()->SetRangeUser(0, 1.5);
  if (sample==razor && binIn==mr) qcd_with_mr[i]->GetYaxis()->SetRangeUser(0, 0.5);
  else if (sample==dijet && binIn==ljpt) qcd_with_mr[i]->GetYaxis()->SetRangeUser(0, 0.5);
  else if (sample==dijet && binIn==rsq) qcd_with_mr[i]->GetYaxis()->SetRangeUser(0, 1.0);
  else if (sample==razor && binIn==ljpt) qcd_with_mr[i]->GetYaxis()->SetRangeUser(0, 1.5);
  qcd_with_mr[i]->GetXaxis()->SetNdivisions(508);
  qcd_with_mr[i]->GetYaxis()->SetNdivisions(508);
  qcd_with_mr[i]->SetTitle("");
  //qcd_with_mr[i]->SetTitle("80 < lead. jet p_T < 200");
  //qcd_with_mr[i]->SetTitle("200 < lead. jet p_T < 300");
  //qcd_with_mr[i]->SetTitle("300 < lead. jet p_T");
  qcd_with_mr[i]->Draw("ap e1");
  dat_with_mr[i]->Draw("p e1 same");
  //if (binIn==mr) qcd_with_mr[i]->Fit(qcd_fxn,"R");
  //qcd_with_mr[i]->Fit(qcd_fxn,"R");
  fxn_plus_err->SetFillColor(kYellow);
  fxn_plus_err->SetFillStyle(1001);
  fxn_plus_err->SetMarkerStyle(0);
  //fxn_plus_err->SetLineColor(kOrange);
  //fxn_plus_err->SetMarkerColor(kOrange);
  fxn_plus_err->Draw("same f e3");
  qcd_with_mr[i]->Draw("same p e1");
  dat_with_mr[i]->Draw("same p e1");
  qcd_fxn->Draw("same");
  //qcd_fxn_2->Draw("same");

  leg->Draw();
  c->SaveAs(pname);

}
