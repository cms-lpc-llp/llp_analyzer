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
#include <TH2D.h>
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

void doTheThing(int box, int nBtags);

void getQCDTranslationFactors2D() {

  doTheThing(0,0);
  doTheThing(0,1);
  doTheThing(0,2);
  doTheThing(0,3);

  doTheThing(1,0);
  doTheThing(1,1);
  doTheThing(1,2);

}
void doTheThing(int box, int nBtags) {

  enum {multijet=0, dijet};

  //int box=dijet;
  //int nBtags=2;

  const Int_t nbinx=5, nbiny=7, nbinz=2;
  //const Int_t nbinx=5, nbiny=5, nbinz=2;
  Float_t xmin=400, xmax=3000; 
  Float_t ymin=0.15,  ymax=0.3; 
  //Float_t ymin=0.15,  ymax=0.25; 
  Float_t zmin=0,    zmax=2;
  Float_t xbins[nbinx+1] = {xmin, 500, 600, 800, 1000, xmax};
  Float_t ybins[nbiny+1] = { ymin, 0.16, 0.18, 0.20, 0.225, 0.25, 0.275, ymax};
  //Float_t ybins[nbiny+1] = { ymin, 0.16, 0.18, 0.20, 0.225, ymax};
  Float_t zbins[nbinz+1] = {zmin, 1, zmax};
  char pname[100];
  char gname[100];
  char dname[100];
  char fname[100];

  TFile *fD = TFile::Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_29Aug2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root","read");
  TFile *fQ = TFile::Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_12Sep2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_QCD_1pb_weighted.root","read");
  TFile *fT = TFile::Open("TTJets.root","read");
  TFile *fW = TFile::Open("WJets.root","read");
  TFile *fZ = TFile::Open("ZInv.root","read");

  //TFile *fT = TFile::Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_12Sep2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root","read");
  //TFile *fW = TFile::Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_12Sep2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_WJets_1pb_weighted.root","read");
  //TFile *fZ = TFile::Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_12Sep2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_ZInv_1pb_weighted.root","read");

  TString cut_str="*(MR>400 && Rsq>0.15)*(Flag_HBHENoiseFilter && Flag_HBHEIsoNoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter && Flag_badChargedCandidateFilter && Flag_badMuonFilter)";

  TString box_cut="(box==11||box==12)";
  if (box==dijet) box_cut="(box==14)";

  TString b_cut;
  if (nBtags==0) b_cut="*(nBTaggedJets==0)";
  else if (nBtags==1) b_cut="*(nBTaggedJets==1)";
  else if (nBtags==2 && box==dijet) b_cut="*(nBTaggedJets>1)";
  else if (nBtags==2) b_cut="*(nBTaggedJets==2)";
  else if (nBtags==3 && box==dijet) cout << "wtf b tags" << endl;
  else if (nBtags==3) b_cut="*(nBTaggedJets>2)";
  else b_cut="";

  TString cut_sf="*mcScaleFactor";
  TString cut_weight="*weight*35800";

  TH3F *dPhiPF_Q = new TH3F("dPhiPF_Q", "dPhiPF_Q", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_Q->Sumw2();
  TH3F *dPhiPF_D = new TH3F("dPhiPF_D", "dPhiPF_D", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_D->Sumw2();
  TH3F *dPhiPF_T = new TH3F("dPhiPF_T", "dPhiPF_T", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_T->Sumw2();
  TH3F *dPhiPF_W = new TH3F("dPhiPF_W", "dPhiPF_W", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_W->Sumw2();
  TH3F *dPhiPF_Z = new TH3F("dPhiPF_Z", "dPhiPF_Z", nbinx, &xbins[0], nbiny, &ybins[0], nbinz, &zbins[0]); dPhiPF_Z->Sumw2();

  // open trees
  TTree *tQ = (TTree*) fQ->Get("RazorInclusive");
  TTree *tD = (TTree*) fD->Get("RazorInclusive");
  TTree *tT = (TTree*) fT->Get("RazorInclusive");
  TTree *tW = (TTree*) fW->Get("RazorInclusive");
  TTree *tZ = (TTree*) fZ->Get("RazorInclusive");

  TCanvas *c = MakeCanvas("c","c",800,600);
  TLegend *leg = new TLegend(0.65, 0.7, 0.9, 0.9);
  leg->SetShadowColor(0);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  tD->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_D", box_cut+cut_str+b_cut);
  tQ->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_Q", box_cut+cut_str+b_cut+cut_weight);
  tT->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_T", box_cut+cut_str+b_cut+cut_weight+cut_sf);
  tW->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_W", box_cut+cut_str+b_cut+cut_weight+cut_sf);
  tZ->Draw("(abs(dPhiRazor)>2.8)+0.5:Rsq:MR>>dPhiPF_Z", box_cut+cut_str+b_cut+cut_weight+cut_sf);

  TFile *outfile = new TFile("qcdTranslationFactors_new.root","update");

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

	Float_t drsq=0.5*dPhiPF_Q->GetYaxis()->GetBinWidth(j+1);
	Float_t dnPF_u=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
	Float_t dnPF_l=nPF*TMath::Sqrt( dnP*dnP + dnF*dnF );
	dnPF_l = (dnPF_l<nPF ? dnPF_l : nPF);

	qcd_with_mr[i]->SetPoint(k, rsq, nPF);
	qcd_with_mr[i]->SetPointError(k, drsq, drsq, dnPF_l, dnPF_u);
	k++;

      }
    }
  }

  if (box==multijet)
    sprintf(pname, "npf_2d_multijet_%ib", nBtags);
  else if (box==dijet)
    sprintf(pname, "npf_2d_dijet_%ib", nBtags);

  vector<TGraphAsymmErrors *> dat_with_mr; 
  TH2D *hTranslationFactors2D = new TH2D(pname,pname, nbinx, &xbins[0], nbiny, &ybins[0]);

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

      //cout << dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,1)) << ", " << dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,2)) << endl;

      nPass=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,1)) - dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,1)) 
	- dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,1)) - dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,1));
      
      nFail=dPhiPF_D->GetBinContent(dPhiPF_D->GetBin(i+1,j+1,2)) - dPhiPF_T->GetBinContent(dPhiPF_T->GetBin(i+1,j+1,2)) 
	- dPhiPF_W->GetBinContent(dPhiPF_W->GetBin(i+1,j+1,2)) - dPhiPF_Z->GetBinContent(dPhiPF_Z->GetBin(i+1,j+1,2));
      
      dnP+=dnP_T*dnP_T+dnP_W*dnP_W+dnP_Z*dnP_Z;
      dnF+=dnF_T*dnF_T+dnF_W*dnF_W+dnF_Z*dnF_Z;

      if (nPass>0 && nFail>0) {
	nPF=nPass/nFail;
	drsq=0.5*dPhiPF_D->GetYaxis()->GetBinWidth(j+1);
	Float_t dnPF_u=nPF*TMath::Sqrt( dnP + dnF );
	Float_t dnPF_l=nPF*TMath::Sqrt( dnP + dnF );
	dnPF_l = (dnPF_l<nPF ? dnPF_l : nPF);
	
	dat_with_mr[i]->SetPoint(k, xrsq, nPF);
	dat_with_mr[i]->SetPointError(k, drsq, drsq, dnPF_l, dnPF_u);

	//cout << i << ", " << j << ", " << xrsq << ", "<< nPF << endl;
	
	hTranslationFactors2D->SetBinContent(i+1,j+1,nPF);

	k++;

      }
    }
  }

  gStyle->SetPalette(56);

  hTranslationFactors2D->Draw("col text");

  if (box==multijet) 
    sprintf(pname, "npf_2d_multijet_%ib.pdf", nBtags);
  else if (box==dijet) 
    sprintf(pname, "npf_2d_dijet_%ib.pdf", nBtags);

  c->SaveAs(pname);

  vector<double> slopes;
  vector<double> slopes_err;
  cout << pname << endl;

  vector<TGraphAsymmErrors *> fit_with_mr; 

  for (Int_t i=0; i<dPhiPF_D->GetNbinsX(); i++) {

    fit_with_mr.push_back(new TGraphAsymmErrors());

    if (box==multijet) {
      sprintf(pname, "npf_multijet_%ib_%i.pdf", nBtags, i);
      sprintf(gname, "mc_npf_multijet_%ib_%i", nBtags, i);
      sprintf(dname, "dat_npf_multijet_%ib_%i", nBtags, i);
      sprintf(fname, "fit_npf_multijet_%ib_%i", nBtags, i);
    }
    else if (box==dijet) {
      sprintf(pname, "npf_dijet_%ib_%i.pdf", nBtags, i);
      sprintf(gname, "mc_npf_dijet_%ib_%i", nBtags, i);
      sprintf(dname, "dat_npf_dijet_%ib_%i", nBtags, i);
      sprintf(fname, "fit_npf_dijet_%ib_%i", nBtags, i);
    }

    TF1 *qcd_fxn = new TF1(fname,"[0]+[1]*x",80,3000);

    dat_with_mr[i]->Fit(fname,"0");

    double p0 = qcd_fxn->GetParameter(0);
    double e0 = qcd_fxn->GetParError(0);

    double p1 = qcd_fxn->GetParameter(1);
    double e1 = qcd_fxn->GetParError(1);

    for (Int_t j=0; j<dat_with_mr[i]->GetN(); j++) {
      double x, y;
      dat_with_mr[i]->GetPoint(j,x,y);

      double vn=p0+p1*x;
      double vu=(p0+e0)+(p1+e1)*x;
      double vd=(p0-e0)+(p1-e1)*x;

      fit_with_mr[i]->SetPoint(j, x+0.002, vn);
      fit_with_mr[i]->SetPointError(j, 0, 0, vn-vd, vu-vn);

      //cout << vn << ", " << vu << ", " << vd << ", " << (p0+e0)+(p1-e1)*x << ", " << (p0-e0)+(p1+e1)*x << endl;
      
    }

    slopes.push_back(qcd_fxn->GetParameter(1));
    slopes_err.push_back(qcd_fxn->GetParError(1));

    qcd_with_mr[i]->SetMarkerStyle(24);
    fit_with_mr[i]->SetMarkerStyle(22);

    qcd_with_mr[i]->GetXaxis()->SetTitle("R^{2}");
    qcd_with_mr[i]->GetYaxis()->SetTitle("Translation Factor #zeta");
    qcd_with_mr[i]->GetYaxis()->SetTitleOffset(1.2);
    qcd_with_mr[i]->GetYaxis()->SetRangeUser(0, 2.0);
    qcd_with_mr[i]->GetXaxis()->SetNdivisions(508);
    qcd_with_mr[i]->GetYaxis()->SetNdivisions(508);
    qcd_with_mr[i]->SetTitle("");
    qcd_with_mr[i]->Draw("ap e1");
    dat_with_mr[i]->Draw("p e1 same");
    
    qcd_with_mr[i]->GetXaxis()->SetRangeUser(399.9,3000.1);

    fit_with_mr[i]->SetMarkerColor(kRed);
    fit_with_mr[i]->SetLineColor(kRed);
    fit_with_mr[i]->Draw("same p e1");    
    qcd_with_mr[i]->Draw("same p e1");
    dat_with_mr[i]->Draw("same p e1");

    leg->Clear();
    leg->AddEntry(dat_with_mr[i], "Data", "lp");
    leg->AddEntry(qcd_with_mr[i], "MC", "lp");
    leg->AddEntry(fit_with_mr[i], "Fit", "lp");

    leg->Draw();

    qcd_fxn->Write(fname);
    qcd_with_mr[i]->Write(gname);
    dat_with_mr[i]->Write(dname);

    c->SaveAs(pname);

  }

  if (box==multijet) sprintf(pname, "npf_multijet_%ib",nBtags);
  else if (box==dijet) sprintf(pname, "npf_dijet_%ib",nBtags);

  outfile->Write();
  outfile->Close();

  //cout << pname << endl;
  //for (uint i=0; i<slopes.size(); i++) {
  //cout << xbins[i] << ", " << xbins[i+1] << ": " << slopes.at(i) << " pm " << slopes_err.at(i) << endl;
  //}

}
