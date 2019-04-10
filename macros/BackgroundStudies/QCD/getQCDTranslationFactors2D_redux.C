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
#include "TAxis.h"
#include "TMath.h"
#include "TH3F.h"
#include "TH2D.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TBenchmark.h"                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#endif

#include "CalStyleRemix.hh"

enum category {dijet=0, multijet, uberjet};

void doTheThing(int box, int nBtags, 
		TString extracuts, TString tag,
		TString var1, vector<float> bins1, int nbins1,
		TString var2, vector<float> bins2, int nbins2,
		TFile *outfile
		);

void getQCDTranslationFactors2D_redux() {

  TFile *outf = new TFile("test.root","recreate");

  TString getLoRegion = "*(MR>400 && MR<1600 && Rsq>0.2 && Rsq<0.3)";

  const int nbins_rsq_lo = 5;
  vector<float> rsq_lo_bins;
  rsq_lo_bins.push_back(0.20);
  rsq_lo_bins.push_back(0.22);
  rsq_lo_bins.push_back(0.24);
  rsq_lo_bins.push_back(0.26);
  rsq_lo_bins.push_back(0.28);
  rsq_lo_bins.push_back(0.30);

  const int nbins_mr_lo = 6;
  vector<float> mr_lo_bins;
  mr_lo_bins.push_back(400);
  mr_lo_bins.push_back(600);
  mr_lo_bins.push_back(800);
  mr_lo_bins.push_back(1000);
  mr_lo_bins.push_back(1200);
  mr_lo_bins.push_back(1400);
  mr_lo_bins.push_back(1600);

  const int nbins_rsq_lo2 = 4;
  vector<float> rsq_lo_bins2;
  rsq_lo_bins2.push_back(0.20);
  rsq_lo_bins2.push_back(0.225);
  rsq_lo_bins2.push_back(0.25);
  rsq_lo_bins2.push_back(0.275);
  rsq_lo_bins2.push_back(0.30);

  const int nbins_mr_lo2 = 4;
  vector<float> mr_lo_bins2;
  mr_lo_bins2.push_back(400);
  mr_lo_bins2.push_back(600);
  mr_lo_bins2.push_back(800);
  mr_lo_bins2.push_back(1000);
  mr_lo_bins2.push_back(1600);

  TString getHiRegion = "*(MR>1600 && MR<4000 && Rsq>0.1 && Rsq<0.2)";

  const int nbins_rsq_hi = 5;
  vector<float> rsq_hi_bins;
  rsq_hi_bins.push_back(0.10);
  rsq_hi_bins.push_back(0.12);
  rsq_hi_bins.push_back(0.14);
  rsq_hi_bins.push_back(0.16);
  rsq_hi_bins.push_back(0.18);
  rsq_hi_bins.push_back(0.20);

  const int nbins_mr_hi = 5;
  vector<float> mr_hi_bins;
  mr_hi_bins.push_back(1600);
  mr_hi_bins.push_back(1800);
  mr_hi_bins.push_back(2000);
  mr_hi_bins.push_back(2500);
  mr_hi_bins.push_back(3000);
  mr_hi_bins.push_back(4000);

  const int nbins_rsq_hi2 = 3;
  vector<float> rsq_hi_bins2;
  rsq_hi_bins2.push_back(0.10);
  rsq_hi_bins2.push_back(0.13);
  rsq_hi_bins2.push_back(0.20);

  const int nbins_mr_hi2 = 2;
  vector<float> mr_hi_bins2;
  mr_hi_bins2.push_back(1600);
  mr_hi_bins2.push_back(2000);
  mr_hi_bins2.push_back(4000);

  //enum category {dijet=0, multijet, uberjet};
//  void doTheThing(int box, int nBtags,
//		  TString extracuts,
//		  TString var1, vector<float> bins1, int nbins1,
//		  TString var2, vector<float> bins2, int nbins2,
//                TFile *outfile
//		  );


  //dijet low region
  doTheThing(0, 0, 
	     getLoRegion, "lo_",
	     "MR",  mr_lo_bins,  nbins_mr_lo,
	     "Rsq", rsq_lo_bins, nbins_rsq_lo,
	     outf);

  doTheThing(0, 1, 
	     getLoRegion, "lo_",
	     "MR",  mr_lo_bins,  nbins_mr_lo,
	     "Rsq", rsq_lo_bins, nbins_rsq_lo,
	     outf);

  doTheThing(0, 2, 
	     getLoRegion, "lo_",
	     "MR",  mr_lo_bins,  nbins_mr_lo,
	     "Rsq", rsq_lo_bins, nbins_rsq_lo,
	     outf);

  //dijet high region
  doTheThing(0, 0, 
	     getHiRegion, "hi_",
	     "MR",  mr_hi_bins,  nbins_mr_hi,
	     "Rsq", rsq_hi_bins, nbins_rsq_hi,
	     outf);

  doTheThing(0, 1, 
	     getHiRegion, "hi_",
	     "MR",  mr_hi_bins,  nbins_mr_hi,
	     "Rsq", rsq_hi_bins, nbins_rsq_hi,
	     outf);

  doTheThing(0, 2, 
	     getHiRegion, "hi_",
	     "MR",  mr_hi_bins,  nbins_mr_hi,
	     "Rsq", rsq_hi_bins, nbins_rsq_hi,
	     outf);

  //multijet low region
  doTheThing(1, 0, 
	     getLoRegion, "lo_",
	     "MR",  mr_lo_bins,  nbins_mr_lo,
	     "Rsq", rsq_lo_bins, nbins_rsq_lo,
	     outf);

  doTheThing(1, 1, 
	     getLoRegion, "lo_",
	     "MR",  mr_lo_bins,  nbins_mr_lo,
	     "Rsq", rsq_lo_bins, nbins_rsq_lo,
	     outf);

  doTheThing(1, 2, 
	     getLoRegion, "lo_",
	     "MR",  mr_lo_bins,  nbins_mr_lo,
	     "Rsq", rsq_lo_bins, nbins_rsq_lo,
	     outf);

  doTheThing(1, 3, 
	     getLoRegion, "lo_",
	     "MR",  mr_lo_bins,  nbins_mr_lo,
	     "Rsq", rsq_lo_bins, nbins_rsq_lo,
	     outf);


  //multijet hi region
  doTheThing(1, 0, 
	     getHiRegion, "hi_",
	     "MR",  mr_hi_bins,  nbins_mr_hi,
	     "Rsq", rsq_hi_bins, nbins_rsq_hi,
	     outf);

  doTheThing(1, 1, 
	     getHiRegion, "hi_",
	     "MR",  mr_hi_bins,  nbins_mr_hi,
	     "Rsq", rsq_hi_bins, nbins_rsq_hi,
	     outf);

  doTheThing(1, 2, 
	     getHiRegion, "hi_",
	     "MR",  mr_hi_bins,  nbins_mr_hi,
	     "Rsq", rsq_hi_bins, nbins_rsq_hi,
	     outf);

  doTheThing(1, 3, 
	     getHiRegion, "hi_",
	     "MR",  mr_hi_bins,  nbins_mr_hi,
	     "Rsq", rsq_hi_bins, nbins_rsq_hi,
	     outf);

  //uberjet low region
  doTheThing(2, 0, 
	     getLoRegion, "lo_",
	     "MR",  mr_lo_bins2,  nbins_mr_lo2,
	     "Rsq", rsq_lo_bins2, nbins_rsq_lo2,
	     outf);

  doTheThing(2, 1, 
	     getLoRegion, "lo_",
	     "MR",  mr_lo_bins2,  nbins_mr_lo2,
	     "Rsq", rsq_lo_bins2, nbins_rsq_lo2,
	     outf);

  doTheThing(2, 2, 
	     getLoRegion, "lo_",
	     "MR",  mr_lo_bins2,  nbins_mr_lo2,
	     "Rsq", rsq_lo_bins2, nbins_rsq_lo2,
	     outf);

  doTheThing(2, 3, 
	     getLoRegion, "lo_",
	     "MR",  mr_lo_bins2,  nbins_mr_lo2,
	     "Rsq", rsq_lo_bins2, nbins_rsq_lo2,
	     outf);


  //uberjet hi region
  doTheThing(2, 0, 
	     getHiRegion, "hi_",
	     "MR",  mr_hi_bins2,  nbins_mr_hi2,
	     "Rsq", rsq_hi_bins2, nbins_rsq_hi2,
	     outf);

  doTheThing(2, 1, 
	     getHiRegion, "hi_",
	     "MR",  mr_hi_bins2,  nbins_mr_hi2,
	     "Rsq", rsq_hi_bins2, nbins_rsq_hi2,
	     outf);

  doTheThing(2, 2, 
	     getHiRegion, "hi_",
	     "MR",  mr_hi_bins2,  nbins_mr_hi2,
	     "Rsq", rsq_hi_bins2, nbins_rsq_hi2,
	     outf);

  doTheThing(2, 3, 
	     getHiRegion, "hi_",
	     "MR",  mr_hi_bins2,  nbins_mr_hi2,
	     "Rsq", rsq_hi_bins2, nbins_rsq_hi2,
	     outf);

  outf->Write();
  outf->Close();
  
}

void doTheThing(int box, int nBtags,		
		TString extracuts, TString tag,
		TString var1, vector<float> bins1, int nbins1,
		TString var2, vector<float> bins2, int nbins2,
		TFile *outfile) {

  const int nbinx=nbins1, nbiny=nbins2, nbinz=2;
  Float_t zbins[nbinz+1] = {0, 1, 2};

  char pname[100], gname[100], dname[100], fname[100];

  TFile *fD = TFile::Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_05Oct2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_Data_NoDuplicates_GoodLumiGolden.root","read");
  TFile *fQ = TFile::Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_05Oct2017/Signal/FullRazorInclusive_Razor2016_MoriondRereco_QCD_1pb_weighted.root","read");
  TFile *fT = TFile::Open("TTJets.root","read");
  TFile *fW = TFile::Open("WJets.root","read");
  TFile *fZ = TFile::Open("ZInv.root","read");

  TString base_cut="*(Flag_HBHENoiseFilter && Flag_HBHEIsoNoiseFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_CSCTightHaloFilter && Flag_badChargedCandidateFilter && Flag_badMuonFilter)";

  TString box_cut="(box==11||box==12||box==14)";
  if (box==dijet) box_cut="(box==14)";
  else if (box==multijet) box_cut="(box==11||box==12)*(nSelectedJets<7)";
  else if (box==uberjet) box_cut="(box==11||box==12)*(nSelectedJets>=7)";

  TString b_cut;
  if (nBtags==0) b_cut="*(nBTaggedJets==0)";
  else if (nBtags==1) b_cut="*(nBTaggedJets==1)";
  else if (nBtags==2 && box==dijet) b_cut="*(nBTaggedJets>1)";
  else if (nBtags==2) b_cut="*(nBTaggedJets==2)";
  else if (nBtags==3 && box==dijet) cout << "wtf b tags" << endl;
  else if (nBtags==3) b_cut="*(nBTaggedJets>2)";
  else b_cut="";

  TString cut_sf="*mcScaleFactor";
  TString cut_norm="*weight*35800";

  TH3F *dPhiPF_Q = new TH3F("dPhiPF_Q"+tag, "dPhiPF_Q"+tag, nbinx, &bins1[0], nbiny, &bins2[0], nbinz, &zbins[0]); dPhiPF_Q->Sumw2();
  TH3F *dPhiPF_D = new TH3F("dPhiPF_D"+tag, "dPhiPF_D"+tag, nbinx, &bins1[0], nbiny, &bins2[0], nbinz, &zbins[0]); dPhiPF_D->Sumw2();
  TH3F *dPhiPF_T = new TH3F("dPhiPF_T"+tag, "dPhiPF_T"+tag, nbinx, &bins1[0], nbiny, &bins2[0], nbinz, &zbins[0]); dPhiPF_T->Sumw2();
  TH3F *dPhiPF_W = new TH3F("dPhiPF_W"+tag, "dPhiPF_W"+tag, nbinx, &bins1[0], nbiny, &bins2[0], nbinz, &zbins[0]); dPhiPF_W->Sumw2();
  TH3F *dPhiPF_Z = new TH3F("dPhiPF_Z"+tag, "dPhiPF_Z"+tag, nbinx, &bins1[0], nbiny, &bins2[0], nbinz, &zbins[0]); dPhiPF_Z->Sumw2();
  
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

  //std::cout << "(abs(dPhiRazor)>2.8)+0.5:"+var2+":"+var1+">>dPhiPF_Q" << std::endl;
  //std::cout << box_cut+extracuts+base_cut+b_cut+cut_norm << std::endl;

  tD->Draw("(abs(dPhiRazor)>2.8)+0.5:"+var2+":"+var1+">>dPhiPF_D"+tag, box_cut+extracuts+base_cut+b_cut);
  tQ->Draw("(abs(dPhiRazor)>2.8)+0.5:"+var2+":"+var1+">>dPhiPF_Q"+tag, box_cut+extracuts+base_cut+b_cut+cut_norm);
  tT->Draw("(abs(dPhiRazor)>2.8)+0.5:"+var2+":"+var1+">>dPhiPF_T"+tag, box_cut+extracuts+base_cut+b_cut+cut_norm+cut_sf);
  tW->Draw("(abs(dPhiRazor)>2.8)+0.5:"+var2+":"+var1+">>dPhiPF_W"+tag, box_cut+extracuts+base_cut+b_cut+cut_norm+cut_sf);
  tZ->Draw("(abs(dPhiRazor)>2.8)+0.5:"+var2+":"+var1+">>dPhiPF_Z"+tag, box_cut+extracuts+base_cut+b_cut+cut_norm+cut_sf);
  
  //std::cout << dPhiPF_Q->GetEntries() << std::endl;

  //TFile *outfile = new TFile("qcdTranslationFactors_4to6jets.root","update");

  outfile->cd();

  vector<TGraphAsymmErrors *> qcd_with_mr; 

  for (Int_t i=0; i<dPhiPF_Q->GetNbinsX(); i++) {
    qcd_with_mr.push_back(new TGraphAsymmErrors());
    //std::cout << "-------" << std::endl;
    //std::cout << dPhiPF_Q->GetXaxis()->GetBinCenter(i+1) << std::endl;
    Int_t k=0;  
    for (Int_t j=0; j<dPhiPF_Q->GetNbinsY(); j++) {
      Float_t rsq=dPhiPF_Q->GetYaxis()->GetBinCenter(j+1);
      //std::cout << "   " << rsq << std::endl;
      Float_t dnP=dPhiPF_Q->GetBinError(dPhiPF_Q->GetBin(i+1,j+1,1))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,1));
      Float_t dnF=dPhiPF_Q->GetBinError(dPhiPF_Q->GetBin(i+1,j+1,2))/dPhiPF_Q->GetBinContent(dPhiPF_Q->GetBin(i+1,j+1,2));

      std::cout << dnP << ", " << dnF << std::endl;
  
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
  else if (box==uberjet)
    sprintf(pname, "npf_2d_uberjet_%ib", nBtags);

  vector<TGraphAsymmErrors *> dat_with_mr; 

  //TH3F *dPhiPF_Q = new TH3F("dPhiPF_Q", "dPhiPF_Q", nbinx, &bins1[0], nbiny, &bins2[0], nbinz, &zbins[0]); dPhiPF_Q->Sumw2();

  TH2D *hTranslationFactors2D = new TH2D(pname+tag,pname+tag, nbinx, &bins1[0], nbiny, &bins2[0]);

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
  else if (box==uberjet) 
    sprintf(pname, "npf_2d_uberjet_%ib.pdf", nBtags);

  c->SaveAs(tag+pname);

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
    else if (box==uberjet) {
      sprintf(pname, "npf_uberjet_%ib_%i.pdf", nBtags, i);
      sprintf(gname, "mc_npf_uberjet_%ib_%i", nBtags, i);
      sprintf(dname, "dat_npf_uberjet_%ib_%i", nBtags, i);
      sprintf(fname, "fit_npf_uberjet_%ib_%i", nBtags, i);
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

    //outfile->cd();

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

    qcd_fxn->Write(tag+fname);
    qcd_with_mr[i]->Write(tag+gname);
    dat_with_mr[i]->Write(tag+dname);

    c->SaveAs(tag+pname);

  }

  //if (box==multijet) sprintf(pname, "npf_multijet_%ib",nBtags);
  //else if (box==dijet) sprintf(pname, "npf_dijet_%ib",nBtags);
  //else if (box==uberjet) sprintf(pname, "npf_uberjet_%ib",nBtags);

  //outfile->Write();
  //outfile->Close();

  //cout << pname << endl;
  //for (uint i=0; i<slopes.size(); i++) {
  //cout << xbins[i] << ", " << xbins[i+1] << ": " << slopes.at(i) << " pm " << slopes_err.at(i) << endl;
  //}

}
