#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>
#include <string.h>
//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>

//#include <tree.hh>

const bool _debug = true;

//Axis
const float axisTitleSize = 0.06;
const float axisTitleOffset = .8;

const float axisTitleSizeRatioX   = 0.18;
const float axisLabelSizeRatioX   = 0.12;
const float axisTitleOffsetRatioX = 0.94;

const float axisTitleSizeRatioY   = 0.15;
const float axisLabelSizeRatioY   = 0.108;
const float axisTitleOffsetRatioY = 0.32;

//Margins
const float leftMargin   = 0.12;
const float rightMargin  = 0.05;
const float topMargin    = 0.07;
const float bottomMargin = 0.12;

//CMS STANDARD
TString CMSText = "CMS";
TString extraText   = "Preliminary";
//TString lumiText = "2.32 fb^{-1} (13 TeV)";
//TString lumiText = "35.9 fb^{-1} (13 TeV)";
TString lumiText = "Simulation (13 TeV)";

bool AddCMS( TCanvas* C );

using namespace std;

int main( int argc, char** argv )
{
  ifstream file;// read input list directly
  std::string line;
  int n = 0;
  std::cout << "[Usage]: ./PlotOnOytputs inputList outputFile\n" << std::endl;
  srand(time(NULL));
  gROOT->Reset();
  gStyle->SetOptStat(0);
  if ( argc == 3 )
  { 
        file.open (argv[1], ios::in | ios::binary);
        std::cout << "[INFO]: Opening file " << argv[1] << " ......" << std::endl;
        std::cout<< std::endl;
        if ( !file.is_open () )
        {
                std::cerr << "!! File open error:" << argv[1] << "; make sure the file is in the correct location" << std::endl;
                return 1;
        } else {
                while(getline(file,line)) ++n;
                std::cout << "n = " << n << "\n" << std::endl;
        }
        //std::string outputFile = argv[2];
  }

  std::ifstream ifs( argv[1], std::ifstream::in );
  assert(ifs);
  //while(std::getline(ifs,line)) n++;

  TFile* fin[n]; 
  TTree* tree[n]; 

  int i = 0;
  std::string process, rootFileName;
  std::string labelProcess[n];
  while ( ifs.good() ){
          ifs >> process >> rootFileName;
          if ( ifs.eof() ) continue;
          if ( process.find("#") != std::string::npos ) continue;
          if ( _debug ) std::cout << process << " " << rootFileName << std::endl;
          labelProcess[i] = process;
          fin[i] = new TFile( rootFileName.c_str(), "READ");
          //assert( fin[i] );
          if ( _debug ) std::cout << "[INFO]: file: " << rootFileName << " passed check\n\n"<< std::endl;

          //------------------------
          //Getting TTree and Histos
          //------------------------
          tree[i] = (TTree*)fin[i]->Get("HggRazorLeptons");

          i++;
          //std::cout << "i = " << i << "\n" << std::endl;
  }

  
  //********************************************************
  //Print output
  //********************************************************
  std::string outputFile = argv[2];
  TFile* fout = new TFile( outputFile.c_str(), "RECREATE");

  int m = 2*n;
  std::cout << "m = " << m << "\n" << std::endl;

  TH1F* hist_MR[m];
  TH1F* hist_Rsq[m]; 
  TH1F* hist_HPt[m]; //Higgs Pt
  TH1F* hist_pTGammaGamma[m]; //reconstructed Higgs Pt
  TH1F* hist_HT[m];  //HT = sum of photon pT  + jet pT
  TH1F* hist_MET[m];

  //MR leptonic
  TCanvas* c_MR_leptonic = new TCanvas( "c_MR_leptonic", "c_MR_leptonic", 800, 700 );
  c_MR_leptonic->SetHighLightColor(2);
  c_MR_leptonic->SetFillColor(0);
  c_MR_leptonic->SetBorderMode(0);
  c_MR_leptonic->SetBorderSize(2);
  c_MR_leptonic->SetLeftMargin( leftMargin );
  c_MR_leptonic->SetRightMargin( 1.6*rightMargin );
  c_MR_leptonic->SetTopMargin( topMargin );
  c_MR_leptonic->SetBottomMargin( bottomMargin );
  c_MR_leptonic->SetFrameBorderMode(0);
  c_MR_leptonic->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_MR_str = "hist_MR_"+Convert.str();
          hist_MR[i] = new TH1F(hist_MR_str.c_str(),"",30,0,2500);
          std::string draw_MR_str = "MR>>"+hist_MR_str;
          tree[i]->Draw(draw_MR_str.c_str(),"box > -1 && box <5 && jet_Pt != 0","HISTSAME");
          if(i<3) hist_MR[i]->SetLineColor(i+2);
          else hist_MR[i]->SetLineColor(6);
          hist_MR[i]->SetLineWidth(3);
          double norm = hist_MR[i]->Integral();
          hist_MR[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_MR[0]->SetTitle(";M_{R} [GeV];Fraction of events");
  hist_MR[0]->GetYaxis()->SetRangeUser(0,0.18);
  hist_MR[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_MR[0]->GetXaxis()->SetTitleOffset(1.4);
   TLatex *   tex = new TLatex(0.27,0.88,"pp #rightarrow #tilde{#chi}^{0,#pm}_{i} #tilde{#chi}^{0,#pm}_{j} #rightarrow  #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} + X_{soft}; #tilde{#chi}^{0}_{1} #rightarrow Z #tilde{G} (100%)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.81,"m_{#tilde{#chi}^{0}_{2}} #approx m_{#tilde{#chi}^{#pm}_{1}} #approx m_{#tilde{#chi}^{0}_{1}};  m_{#tilde{G}} = 1 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.74,"m_{#chi} = 127 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(2);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.67,"m_{#chi} = 300 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(3);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.60,"m_{#chi} = 600 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(4);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.53,"m_{#chi} = 900 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(6);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.66,0.64,"Leptonic bins");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_MR_leptonic->Write();
  c_MR_leptonic->Write();
  c_MR_leptonic->Update();
  AddCMS(c_MR_leptonic);
  c_MR_leptonic->SaveAs("c_MR_leptonic.pdf");
  c_MR_leptonic->SaveAs("c_MR_leptonic.png");
  c_MR_leptonic->SaveAs("c_MR_leptonic.C");


  //MR hadronic
  TCanvas* c_MR_hadronic = new TCanvas( "c_MR_hadronic", "c_MR_hadronic", 800, 700 );
  c_MR_hadronic->SetHighLightColor(2);
  c_MR_hadronic->SetFillColor(0);
  c_MR_hadronic->SetBorderMode(0);
  c_MR_hadronic->SetBorderSize(2);
  c_MR_hadronic->SetLeftMargin( leftMargin );
  c_MR_hadronic->SetRightMargin( 1.6*rightMargin );
  c_MR_hadronic->SetTopMargin( topMargin );
  c_MR_hadronic->SetBottomMargin( bottomMargin );
  c_MR_hadronic->SetFrameBorderMode(0);
  c_MR_hadronic->SetFrameBorderMode(0);
  for(int i = n; i < m; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_MR_str = "hist_MR_"+Convert.str();
          hist_MR[i] = new TH1F(hist_MR_str.c_str(),"",30,0,2500);
          std::string draw_MR_str = "MR>>"+hist_MR_str;
          tree[i-n]->Draw(draw_MR_str.c_str(),"box > 4 && box <10 && jet_Pt != 0","HISTSAME");
          if(i<3+n) hist_MR[i]->SetLineColor(i-n+2);
          else hist_MR[i]->SetLineColor(6);
          hist_MR[i]->SetLineWidth(3);
          double norm = hist_MR[i]->Integral();
          hist_MR[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_MR[n]->SetTitle(";M_{R} [GeV];Fraction of events");
  hist_MR[n]->GetYaxis()->SetRangeUser(0,0.18);
  hist_MR[n]->GetYaxis()->SetTitleOffset(1.4);
  hist_MR[n]->GetXaxis()->SetTitleOffset(1.4);
      tex = new TLatex(0.27,0.88,"pp #rightarrow #tilde{#chi}^{0,#pm}_{i} #tilde{#chi}^{0,#pm}_{j} #rightarrow  #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} + X_{soft}; #tilde{#chi}^{0}_{1} #rightarrow Z #tilde{G} (100%)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.81,"m_{#tilde{#chi}^{0}_{2}} #approx m_{#tilde{#chi}^{#pm}_{1}} #approx m_{#tilde{#chi}^{0}_{1}};  m_{#tilde{G}} = 1 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.74,"m_{#chi} = 127 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(2);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.67,"m_{#chi} = 300 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(3);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.60,"m_{#chi} = 600 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(4);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.53,"m_{#chi} = 900 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(6);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.66,0.64,"Hadronic bins");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_MR_hadronic->Write();
  c_MR_hadronic->Write();
  c_MR_hadronic->Update();
  AddCMS(c_MR_hadronic);
  c_MR_hadronic->SaveAs("c_MR_hadronic.pdf");
  c_MR_hadronic->SaveAs("c_MR_hadronic.png");
  c_MR_hadronic->SaveAs("c_MR_hadronic.C");


  //MR comparison
  TCanvas* c_MR_comparison = new TCanvas( "c_MR_comparison", "c_MR_comparison", 800, 700 );
  c_MR_comparison->SetHighLightColor(2);
  c_MR_comparison->SetFillColor(0);
  c_MR_comparison->SetBorderMode(0);
  c_MR_comparison->SetBorderSize(2);
  c_MR_comparison->SetLeftMargin( leftMargin );
  c_MR_comparison->SetRightMargin( 1.6*rightMargin );
  c_MR_comparison->SetTopMargin( topMargin );
  c_MR_comparison->SetBottomMargin( bottomMargin );
  c_MR_comparison->SetFrameBorderMode(0);
  c_MR_comparison->SetFrameBorderMode(0);
  for(int i = 2; i < 3; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_MR_str = "hist_MR_"+Convert.str();
          hist_MR[i] = new TH1F(hist_MR_str.c_str(),"",30,0,2500);
          std::string draw_MR_str = "MR>>"+hist_MR_str;
          tree[i]->Draw(draw_MR_str.c_str(),"box > -1 && box <5 && jet_Pt != 0","HISTSAME");
          hist_MR[i]->SetLineColor(4);
          hist_MR[i]->SetLineWidth(3);
          double norm = hist_MR[i]->Integral();
          hist_MR[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  for(int i = 2+n; i < 3+n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_MR_str = "hist_MR_"+Convert.str();
          hist_MR[i] = new TH1F(hist_MR_str.c_str(),"",30,0,2500);
          std::string draw_MR_str = "MR>>"+hist_MR_str;
          tree[i-n]->Draw(draw_MR_str.c_str(),"box > 4 && box <10 && jet_Pt != 0","HISTSAME");
          hist_MR[i]->SetLineColor(6);
          hist_MR[i]->SetLineWidth(3);
          double norm = hist_MR[i]->Integral();
          hist_MR[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_MR[2]->SetTitle(";M_{R} [GeV];Fraction of events");
  hist_MR[2]->GetYaxis()->SetRangeUser(0,0.12);
  hist_MR[2]->GetYaxis()->SetTitleOffset(1.4);
  hist_MR[2]->GetXaxis()->SetTitleOffset(1.4);
      tex = new TLatex(0.27,0.88,"pp #rightarrow #tilde{#chi}^{0,#pm}_{i} #tilde{#chi}^{0,#pm}_{j} #rightarrow  #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} + X_{soft}; #tilde{#chi}^{0}_{1} #rightarrow Z #tilde{G} (100%)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.81,"m_{#tilde{#chi}^{0}_{2}} #approx m_{#tilde{#chi}^{#pm}_{1}} #approx m_{#tilde{#chi}^{0}_{1}};  m_{#tilde{G}} = 1 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.74,"Leptonic bins m_{#chi} = 600 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(4);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.67,"Hadronic bins m_{#chi} = 600 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(6);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_MR_comparison->Write();
  c_MR_comparison->Write();
  c_MR_comparison->Update();
  AddCMS(c_MR_comparison);
  c_MR_comparison->SaveAs("c_MR_comparison.pdf");
  c_MR_comparison->SaveAs("c_MR_comparison.png");
  c_MR_comparison->SaveAs("c_MR_comparison.C");


  //Rsq leptonic
  TCanvas* c_Rsq_leptonic = new TCanvas( "c_Rsq_leptonic", "c_Rsq_leptonic", 800, 700 );
  c_Rsq_leptonic->SetHighLightColor(2);
  c_Rsq_leptonic->SetFillColor(0);
  c_Rsq_leptonic->SetBorderMode(0);
  c_Rsq_leptonic->SetBorderSize(2);
  c_Rsq_leptonic->SetLeftMargin( leftMargin );
  c_Rsq_leptonic->SetRightMargin( 1.6*rightMargin );
  c_Rsq_leptonic->SetTopMargin( topMargin );
  c_Rsq_leptonic->SetBottomMargin( bottomMargin );
  c_Rsq_leptonic->SetFrameBorderMode(0);
  c_Rsq_leptonic->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_Rsq_str = "hist_Rsq_"+Convert.str();
          hist_Rsq[i] = new TH1F(hist_Rsq_str.c_str(),"",30,0,1.);
          std::string draw_Rsq_str = "Rsq>>"+hist_Rsq_str;
          tree[i]->Draw(draw_Rsq_str.c_str(),"box > -1 && box <5 && jet_Pt != 0","HISTSAME");
          if(i<3) hist_Rsq[i]->SetLineColor(i+2);
          else hist_Rsq[i]->SetLineColor(6);
          hist_Rsq[i]->SetLineWidth(3);
          double norm = hist_Rsq[i]->Integral();
          hist_Rsq[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_Rsq[0]->SetTitle(";R^{2};Fraction of events");
  hist_Rsq[0]->GetYaxis()->SetRangeUser(0,0.40);
  hist_Rsq[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_Rsq[0]->GetXaxis()->SetTitleOffset(1.4);
      tex = new TLatex(0.27,0.88,"pp #rightarrow #tilde{#chi}^{0,#pm}_{i} #tilde{#chi}^{0,#pm}_{j} #rightarrow  #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} + X_{soft}; #tilde{#chi}^{0}_{1} #rightarrow Z #tilde{G} (100%)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.81,"m_{#tilde{#chi}^{0}_{2}} #approx m_{#tilde{#chi}^{#pm}_{1}} #approx m_{#tilde{#chi}^{0}_{1}};  m_{#tilde{G}} = 1 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.74,"m_{#chi} = 127 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(2);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.67,"m_{#chi} = 300 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(3);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.60,"m_{#chi} = 600 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(4);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.53,"m_{#chi} = 900 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(6);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.66,0.64,"Leptonic bins");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_Rsq_leptonic->Write();
  c_Rsq_leptonic->Write();
  c_Rsq_leptonic->Update();
  AddCMS(c_Rsq_leptonic);
  c_Rsq_leptonic->SaveAs("c_Rsq_leptonic.pdf");
  c_Rsq_leptonic->SaveAs("c_Rsq_leptonic.png");
  c_Rsq_leptonic->SaveAs("c_Rsq_leptonic.C");


  //Rsq hadronic
  TCanvas* c_Rsq_hadronic = new TCanvas( "c_Rsq_hadronic", "c_Rsq_hadronic", 800, 700 );
  c_Rsq_hadronic->SetHighLightColor(2);
  c_Rsq_hadronic->SetFillColor(0);
  c_Rsq_hadronic->SetBorderMode(0);
  c_Rsq_hadronic->SetBorderSize(2);
  c_Rsq_hadronic->SetLeftMargin( leftMargin );
  c_Rsq_hadronic->SetRightMargin( 1.6*rightMargin );
  c_Rsq_hadronic->SetTopMargin( topMargin );
  c_Rsq_hadronic->SetBottomMargin( bottomMargin );
  c_Rsq_hadronic->SetFrameBorderMode(0);
  c_Rsq_hadronic->SetFrameBorderMode(0);
  for(int i = n; i < m; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_Rsq_str = "hist_Rsq_"+Convert.str();
          hist_Rsq[i] = new TH1F(hist_Rsq_str.c_str(),"",30,0,1.);
          std::string draw_Rsq_str = "Rsq>>"+hist_Rsq_str;
          tree[i-n]->Draw(draw_Rsq_str.c_str(),"box > 4 && box <10 && jet_Pt != 0","HISTSAME");
          if(i<3+n) hist_Rsq[i]->SetLineColor(i-n+2);
          else hist_Rsq[i]->SetLineColor(6);
          hist_Rsq[i]->SetLineWidth(3);
          double norm = hist_Rsq[i]->Integral();
          hist_Rsq[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_Rsq[n]->SetTitle(";R^{2};Fraction of events");
  hist_Rsq[n]->GetYaxis()->SetRangeUser(0,0.25);
  hist_Rsq[n]->GetYaxis()->SetTitleOffset(1.4);
  hist_Rsq[n]->GetXaxis()->SetTitleOffset(1.4);
      tex = new TLatex(0.27,0.88,"pp #rightarrow #tilde{#chi}^{0,#pm}_{i} #tilde{#chi}^{0,#pm}_{j} #rightarrow  #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} + X_{soft}; #tilde{#chi}^{0}_{1} #rightarrow Z #tilde{G} (100%)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.81,"m_{#tilde{#chi}^{0}_{2}} #approx m_{#tilde{#chi}^{#pm}_{1}} #approx m_{#tilde{#chi}^{0}_{1}};  m_{#tilde{G}} = 1 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.74,"m_{#chi} = 127 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(2);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.67,"m_{#chi} = 300 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(3);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.60,"m_{#chi} = 600 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(4);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.53,"m_{#chi} = 900 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(6);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.66,0.64,"Hadronic bins");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_Rsq_hadronic->Write();
  c_Rsq_hadronic->Write();
  c_Rsq_hadronic->Update();
  AddCMS(c_Rsq_hadronic);
  c_Rsq_hadronic->SaveAs("c_Rsq_hadronic.pdf");
  c_Rsq_hadronic->SaveAs("c_Rsq_hadronic.png");
  c_Rsq_hadronic->SaveAs("c_Rsq_hadronic.C");


  //Rsq comparison
  TCanvas* c_Rsq_comparison = new TCanvas( "c_Rsq_comparison", "c_Rsq_comparison", 800, 700 );
  c_Rsq_comparison->SetHighLightColor(2);
  c_Rsq_comparison->SetFillColor(0);
  c_Rsq_comparison->SetBorderMode(0);
  c_Rsq_comparison->SetBorderSize(2);
  c_Rsq_comparison->SetLeftMargin( leftMargin );
  c_Rsq_comparison->SetRightMargin( 1.6*rightMargin );
  c_Rsq_comparison->SetTopMargin( topMargin );
  c_Rsq_comparison->SetBottomMargin( bottomMargin );
  c_Rsq_comparison->SetFrameBorderMode(0);
  c_Rsq_comparison->SetFrameBorderMode(0);
  for(int i = 2; i < 3; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_Rsq_str = "hist_Rsq_"+Convert.str();
          hist_Rsq[i] = new TH1F(hist_Rsq_str.c_str(),"",30,0,1.);
          std::string draw_Rsq_str = "Rsq>>"+hist_Rsq_str;
          tree[i]->Draw(draw_Rsq_str.c_str(),"box > -1 && box <5 && jet_Pt != 0","HISTSAME");
          hist_Rsq[i]->SetLineColor(4);
          hist_Rsq[i]->SetLineWidth(3);
          double norm = hist_Rsq[i]->Integral();
          hist_Rsq[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  for(int i = 2+n; i < 3+n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_Rsq_str = "hist_Rsq_"+Convert.str();
          hist_Rsq[i] = new TH1F(hist_Rsq_str.c_str(),"",30,0,1.);
          std::string draw_Rsq_str = "Rsq>>"+hist_Rsq_str;
          tree[i-n]->Draw(draw_Rsq_str.c_str(),"box > 4 && box <10 && jet_Pt != 0","HISTSAME");
          hist_Rsq[i]->SetLineColor(6);
          hist_Rsq[i]->SetLineWidth(3);
          double norm = hist_Rsq[i]->Integral();
          hist_Rsq[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_Rsq[2]->SetTitle(";R^{2};Fraction of events");
  hist_Rsq[2]->GetYaxis()->SetRangeUser(0,0.12);
  hist_Rsq[2]->GetYaxis()->SetTitleOffset(1.4);
  hist_Rsq[2]->GetXaxis()->SetTitleOffset(1.4);
      tex = new TLatex(0.27,0.88,"pp #rightarrow #tilde{#chi}^{0,#pm}_{i} #tilde{#chi}^{0,#pm}_{j} #rightarrow  #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} + X_{soft}; #tilde{#chi}^{0}_{1} #rightarrow Z #tilde{G} (100%)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.81,"m_{#tilde{#chi}^{0}_{2}} #approx m_{#tilde{#chi}^{#pm}_{1}} #approx m_{#tilde{#chi}^{0}_{1}};  m_{#tilde{G}} = 1 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.74,"Leptonic bins m_{#chi} = 600 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(4);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.67,"Hadronic bins m_{#chi} = 600 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(6);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_Rsq_comparison->Write();
  c_Rsq_comparison->Write();
  c_Rsq_comparison->Update();
  AddCMS(c_Rsq_comparison);
  c_Rsq_comparison->SaveAs("c_Rsq_comparison.pdf");
  c_Rsq_comparison->SaveAs("c_Rsq_comparison.png");
  c_Rsq_comparison->SaveAs("c_Rsq_comparison.C");


  //pTGammaGamma leptonic
  TCanvas* c_pTGammaGamma_leptonic = new TCanvas( "c_pTGammaGamma_leptonic", "c_pTGammaGamma_leptonic", 800, 700 );
  c_pTGammaGamma_leptonic->SetHighLightColor(2);
  c_pTGammaGamma_leptonic->SetFillColor(0);
  c_pTGammaGamma_leptonic->SetBorderMode(0);
  c_pTGammaGamma_leptonic->SetBorderSize(2);
  c_pTGammaGamma_leptonic->SetLeftMargin( leftMargin );
  c_pTGammaGamma_leptonic->SetRightMargin( 1.6*rightMargin );
  c_pTGammaGamma_leptonic->SetTopMargin( topMargin );
  c_pTGammaGamma_leptonic->SetBottomMargin( bottomMargin );
  c_pTGammaGamma_leptonic->SetFrameBorderMode(0);
  c_pTGammaGamma_leptonic->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pTGammaGamma_str = "hist_pTGammaGamma_"+Convert.str();
          hist_pTGammaGamma[i] = new TH1F(hist_pTGammaGamma_str.c_str(),"",30,0,2000);
          std::string draw_pTGammaGamma_str = "pTGammaGamma>>"+hist_pTGammaGamma_str;
          tree[i]->Draw(draw_pTGammaGamma_str.c_str(),"box > -1 && box <5 && jet_Pt != 0","HISTSAME");
          if(i<3) hist_pTGammaGamma[i]->SetLineColor(i+2);
          else hist_pTGammaGamma[i]->SetLineColor(6);
          hist_pTGammaGamma[i]->SetLineWidth(3);
          double norm = hist_pTGammaGamma[i]->Integral();
          hist_pTGammaGamma[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pTGammaGamma[0]->SetTitle(";p_{T}^{#gamma #gamma} [GeV];Fraction of events");
  hist_pTGammaGamma[0]->GetYaxis()->SetRangeUser(0,0.36);
  hist_pTGammaGamma[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pTGammaGamma[0]->GetXaxis()->SetTitleOffset(1.4);
      tex = new TLatex(0.27,0.88,"pp #rightarrow #tilde{#chi}^{0,#pm}_{i} #tilde{#chi}^{0,#pm}_{j} #rightarrow  #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} + X_{soft}; #tilde{#chi}^{0}_{1} #rightarrow Z #tilde{G} (100%)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.81,"m_{#tilde{#chi}^{0}_{2}} #approx m_{#tilde{#chi}^{#pm}_{1}} #approx m_{#tilde{#chi}^{0}_{1}};  m_{#tilde{G}} = 1 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.74,"m_{#chi} = 127 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(2);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.67,"m_{#chi} = 300 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(3);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.60,"m_{#chi} = 600 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(4);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.53,"m_{#chi} = 900 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(6);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.66,0.64,"Leptonic bins");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pTGammaGamma_leptonic->Write();
  c_pTGammaGamma_leptonic->Write();
  c_pTGammaGamma_leptonic->Update();
  AddCMS(c_pTGammaGamma_leptonic);
  c_pTGammaGamma_leptonic->SaveAs("c_pTGammaGamma_leptonic.pdf");
  c_pTGammaGamma_leptonic->SaveAs("c_pTGammaGamma_leptonic.png");
  c_pTGammaGamma_leptonic->SaveAs("c_pTGammaGamma_leptonic.C");


  //pTGammaGamma hadronic
  TCanvas* c_pTGammaGamma_hadronic = new TCanvas( "c_pTGammaGamma_hadronic", "c_pTGammaGamma_hadronic", 800, 700 );
  c_pTGammaGamma_hadronic->SetHighLightColor(2);
  c_pTGammaGamma_hadronic->SetFillColor(0);
  c_pTGammaGamma_hadronic->SetBorderMode(0);
  c_pTGammaGamma_hadronic->SetBorderSize(2);
  c_pTGammaGamma_hadronic->SetLeftMargin( leftMargin );
  c_pTGammaGamma_hadronic->SetRightMargin( 1.6*rightMargin );
  c_pTGammaGamma_hadronic->SetTopMargin( topMargin );
  c_pTGammaGamma_hadronic->SetBottomMargin( bottomMargin );
  c_pTGammaGamma_hadronic->SetFrameBorderMode(0);
  c_pTGammaGamma_hadronic->SetFrameBorderMode(0);
  for(int i = n; i < m; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pTGammaGamma_str = "hist_pTGammaGamma_"+Convert.str();
          hist_pTGammaGamma[i] = new TH1F(hist_pTGammaGamma_str.c_str(),"",30,0,2000);
          std::string draw_pTGammaGamma_str = "pTGammaGamma>>"+hist_pTGammaGamma_str;
          tree[i-n]->Draw(draw_pTGammaGamma_str.c_str(),"box > 4 && box <10 && jet_Pt != 0","HISTSAME");
          if(i<3+n) hist_pTGammaGamma[i]->SetLineColor(i-n+2);
          else hist_pTGammaGamma[i]->SetLineColor(6);
          hist_pTGammaGamma[i]->SetLineWidth(3);
          double norm = hist_pTGammaGamma[i]->Integral();
          hist_pTGammaGamma[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pTGammaGamma[n]->SetTitle(";p_{T}^{#gamma #gamma} [GeV];Fraction of events");
  hist_pTGammaGamma[n]->GetYaxis()->SetRangeUser(0,0.36);
  hist_pTGammaGamma[n]->GetYaxis()->SetTitleOffset(1.4);
  hist_pTGammaGamma[n]->GetXaxis()->SetTitleOffset(1.4);
      tex = new TLatex(0.27,0.88,"pp #rightarrow #tilde{#chi}^{0,#pm}_{i} #tilde{#chi}^{0,#pm}_{j} #rightarrow  #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} + X_{soft}; #tilde{#chi}^{0}_{1} #rightarrow Z #tilde{G} (100%)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.81,"m_{#tilde{#chi}^{0}_{2}} #approx m_{#tilde{#chi}^{#pm}_{1}} #approx m_{#tilde{#chi}^{0}_{1}};  m_{#tilde{G}} = 1 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.74,"m_{#chi} = 127 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(2);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.67,"m_{#chi} = 300 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(3);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.60,"m_{#chi} = 600 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(4);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.71,0.53,"m_{#chi} = 900 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(6);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.66,0.64,"Hadronic bins");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pTGammaGamma_hadronic->Write();
  c_pTGammaGamma_hadronic->Write();
  c_pTGammaGamma_hadronic->Update();
  AddCMS(c_pTGammaGamma_hadronic);
  c_pTGammaGamma_hadronic->SaveAs("c_pTGammaGamma_hadronic.pdf");
  c_pTGammaGamma_hadronic->SaveAs("c_pTGammaGamma_hadronic.png");
  c_pTGammaGamma_hadronic->SaveAs("c_pTGammaGamma_hadronic.C");


  //pTGammaGamma comparison
  TCanvas* c_pTGammaGamma_comparison = new TCanvas( "c_pTGammaGamma_comparison", "c_pTGammaGamma_comparison", 800, 700 );
  c_pTGammaGamma_comparison->SetHighLightColor(2);
  c_pTGammaGamma_comparison->SetFillColor(0);
  c_pTGammaGamma_comparison->SetBorderMode(0);
  c_pTGammaGamma_comparison->SetBorderSize(2);
  c_pTGammaGamma_comparison->SetLeftMargin( leftMargin );
  c_pTGammaGamma_comparison->SetRightMargin( 1.6*rightMargin );
  c_pTGammaGamma_comparison->SetTopMargin( topMargin );
  c_pTGammaGamma_comparison->SetBottomMargin( bottomMargin );
  c_pTGammaGamma_comparison->SetFrameBorderMode(0);
  c_pTGammaGamma_comparison->SetFrameBorderMode(0);
  for(int i = 2; i < 3; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pTGammaGamma_str = "hist_pTGammaGamma_"+Convert.str();
          hist_pTGammaGamma[i] = new TH1F(hist_pTGammaGamma_str.c_str(),"",30,0,2000);
          std::string draw_pTGammaGamma_str = "pTGammaGamma>>"+hist_pTGammaGamma_str;
          tree[i]->Draw(draw_pTGammaGamma_str.c_str(),"box > -1 && box <5 && jet_Pt != 0","HISTSAME");
          hist_pTGammaGamma[i]->SetLineColor(4);
          hist_pTGammaGamma[i]->SetLineWidth(3);
          double norm = hist_pTGammaGamma[i]->Integral();
          hist_pTGammaGamma[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  for(int i = 2+n; i < 3+n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pTGammaGamma_str = "hist_pTGammaGamma_"+Convert.str();
          hist_pTGammaGamma[i] = new TH1F(hist_pTGammaGamma_str.c_str(),"",30,0,2000);
          std::string draw_pTGammaGamma_str = "pTGammaGamma>>"+hist_pTGammaGamma_str;
          tree[i-n]->Draw(draw_pTGammaGamma_str.c_str(),"box > 4 && box <10 && jet_Pt != 0","HISTSAME");
          hist_pTGammaGamma[i]->SetLineColor(6);
          hist_pTGammaGamma[i]->SetLineWidth(3);
          double norm = hist_pTGammaGamma[i]->Integral();
          hist_pTGammaGamma[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pTGammaGamma[2]->SetTitle(";p_{T}^{#gamma #gamma} [GeV];Fraction of events");
  hist_pTGammaGamma[2]->GetYaxis()->SetRangeUser(0,0.16);
  hist_pTGammaGamma[2]->GetYaxis()->SetTitleOffset(1.4);
  hist_pTGammaGamma[2]->GetXaxis()->SetTitleOffset(1.4);
      tex = new TLatex(0.27,0.88,"pp #rightarrow #tilde{#chi}^{0,#pm}_{i} #tilde{#chi}^{0,#pm}_{j} #rightarrow  #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} + X_{soft}; #tilde{#chi}^{0}_{1} #rightarrow Z #tilde{G} (100%)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.81,"m_{#tilde{#chi}^{0}_{2}} #approx m_{#tilde{#chi}^{#pm}_{1}} #approx m_{#tilde{#chi}^{0}_{1}};  m_{#tilde{G}} = 1 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.74,"Leptonic bins m_{#chi} = 600 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(4);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
      tex = new TLatex(0.52,0.67,"Hadronic bins m_{#chi} = 600 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextColor(6);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pTGammaGamma_comparison->Write();
  c_pTGammaGamma_comparison->Write();
  c_pTGammaGamma_comparison->Update();
  AddCMS(c_pTGammaGamma_comparison);
  c_pTGammaGamma_comparison->SaveAs("c_pTGammaGamma_comparison.pdf");
  c_pTGammaGamma_comparison->SaveAs("c_pTGammaGamma_comparison.png");
  c_pTGammaGamma_comparison->SaveAs("c_pTGammaGamma_comparison.C");


  fout->Close();


}

bool AddCMS( TCanvas* C )
{
  C->cd();
  float lumix = 0.925;
  float lumiy = 0.945;
  float lumifont = 42;
  
  float cmsx = 0.225;
  float cmsy = 0.940;
  float cmsTextFont   = 61;  // default is helvetic-bold
  float extrax = cmsx + 0.198;
  float extray = cmsy;
  float extraTextFont = 52;  // default is helvetica-italics
  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;
  float cmsSize = 0.06;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    
  float extraTextSize = extraOverCmsTextSize*cmsSize;
  latex.SetTextFont(lumifont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(cmsSize);    
  latex.DrawLatex(lumix, lumiy,lumiText);

  latex.SetTextFont(cmsTextFont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(cmsSize);
  latex.DrawLatex(cmsx, cmsy, CMSText);
   
  latex.SetTextFont(extraTextFont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(extraTextSize);
  latex.DrawLatex(extrax, extray, extraText);
  return true;
};
