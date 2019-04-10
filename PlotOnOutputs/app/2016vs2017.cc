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
TString lumiText = "Simulation (13 TeV)";

bool AddCMS( TCanvas* C );

using namespace std;

int main( int argc, char** argv )
{
  ifstream file;// read input list directly
  std::string line;
  int n = 0;
  std::cout << "[Usage]: ./2016vs2017 inputList outputFile\n" << std::endl;
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

  int m = 1*n;
  std::cout << "m = " << m << "\n" << std::endl;

  //mgg, with the blinded region. ptgg, pho1Pt, pho1Eta, pho1Phi, pho2Pt, pho2Eta, pho2Phi

  TH1F* hist_pTGammaGamma[m]; //reconstructed ptgg
  TH1F* hist_mGammaGamma[m]; //reconstructed mgg
  TH1F* hist_ptggOvermgg[m]; //reconstructed mgg
  TH1F* hist_pho1Pt[m]; 
  TH1F* hist_pho1Eta[m]; 
  TH1F* hist_pho1Phi[m]; 
  TH1F* hist_pho2Pt[m]; 
  TH1F* hist_pho2Eta[m]; 
  TH1F* hist_pho2Phi[m]; 


  //pTGammaGamma 
  TCanvas* c_pTGammaGamma = new TCanvas( "c_pTGammaGamma", "c_pTGammaGamma", 800, 700 );
  c_pTGammaGamma->SetHighLightColor(2);
  c_pTGammaGamma->SetFillColor(0);
  c_pTGammaGamma->SetBorderMode(0);
  c_pTGammaGamma->SetBorderSize(2);
  c_pTGammaGamma->SetLeftMargin( leftMargin );
  c_pTGammaGamma->SetRightMargin( 1.6*rightMargin );
  c_pTGammaGamma->SetTopMargin( topMargin );
  c_pTGammaGamma->SetBottomMargin( bottomMargin );
  c_pTGammaGamma->SetFrameBorderMode(0);
  c_pTGammaGamma->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pTGammaGamma_str = "hist_pTGammaGamma_"+Convert.str();
          hist_pTGammaGamma[i] = new TH1F(hist_pTGammaGamma_str.c_str(),"",30,0,200);
          std::string draw_pTGammaGamma_str = "pTGammaGamma>>"+hist_pTGammaGamma_str;
          if(i<1) tree[i]->Draw(draw_pTGammaGamma_str.c_str(),"( ( (mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.) ) && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 )","HISTSAME");
          else tree[i]->Draw(draw_pTGammaGamma_str.c_str(),"( ( (mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.) ) && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 )","*ESAME"); 
          hist_pTGammaGamma[i]->SetLineColor(i+2);
          hist_pTGammaGamma[i]->SetLineWidth(3);
          double norm = hist_pTGammaGamma[i]->Integral();
          hist_pTGammaGamma[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pTGammaGamma[0]->SetTitle(";p_{T}^{#gamma #gamma} [GeV];Fraction of events");
  hist_pTGammaGamma[0]->GetYaxis()->SetRangeUser(0,0.2);
  hist_pTGammaGamma[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pTGammaGamma[0]->GetXaxis()->SetTitleOffset(1.4);
   TLatex *   tex = new TLatex(0.27,0.88,"data 2016");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.78,"data 2017");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pTGammaGamma->Write();
  c_pTGammaGamma->Write();
  c_pTGammaGamma->Update();
  AddCMS(c_pTGammaGamma);
  c_pTGammaGamma->SaveAs("c_pTGammaGamma.pdf");
  c_pTGammaGamma->SaveAs("c_pTGammaGamma.png");
  c_pTGammaGamma->SaveAs("c_pTGammaGamma.C");

  //mGammaGamma 
  TCanvas* c_mGammaGamma = new TCanvas( "c_mGammaGamma", "c_mGammaGamma", 800, 700 );
  c_mGammaGamma->SetHighLightColor(2);
  c_mGammaGamma->SetFillColor(0);
  c_mGammaGamma->SetBorderMode(0);
  c_mGammaGamma->SetBorderSize(2);
  c_mGammaGamma->SetLeftMargin( leftMargin );
  c_mGammaGamma->SetRightMargin( 1.6*rightMargin );
  c_mGammaGamma->SetTopMargin( topMargin );
  c_mGammaGamma->SetBottomMargin( bottomMargin );
  c_mGammaGamma->SetFrameBorderMode(0);
  c_mGammaGamma->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_mGammaGamma_str = "hist_mGammaGamma_"+Convert.str();
          hist_mGammaGamma[i] = new TH1F(hist_mGammaGamma_str.c_str(),"",30,105,180);
          std::string draw_mGammaGamma_str = "mGammaGamma>>"+hist_mGammaGamma_str;
          if(i<1) tree[i]->Draw(draw_mGammaGamma_str.c_str(),"( ( (mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.) ) && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 )","HISTSAME");
          else tree[i]->Draw(draw_mGammaGamma_str.c_str(),"( ( (mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.) ) && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 )","*ESAME"); 
          hist_mGammaGamma[i]->SetLineColor(i+2);
          hist_mGammaGamma[i]->SetLineWidth(3);
          double norm = hist_mGammaGamma[i]->Integral();
          hist_mGammaGamma[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_mGammaGamma[0]->SetTitle(";m_{#gamma #gamma} [GeV];Fraction of events");
  hist_mGammaGamma[0]->GetYaxis()->SetRangeUser(0,0.13);
  hist_mGammaGamma[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_mGammaGamma[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.88,"data 2016");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.78,"data 2017");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_mGammaGamma->Write();
  c_mGammaGamma->Write();
  c_mGammaGamma->Update();
  AddCMS(c_mGammaGamma);
  c_mGammaGamma->SaveAs("c_mGammaGamma.pdf");
  c_mGammaGamma->SaveAs("c_mGammaGamma.png");
  c_mGammaGamma->SaveAs("c_mGammaGamma.C");

  //ptggOvermgg 
  TCanvas* c_ptggOvermgg = new TCanvas( "c_ptggOvermgg", "c_ptggOvermgg", 800, 700 );
  c_ptggOvermgg->SetHighLightColor(2);
  c_ptggOvermgg->SetFillColor(0);
  c_ptggOvermgg->SetBorderMode(0);
  c_ptggOvermgg->SetBorderSize(2);
  c_ptggOvermgg->SetLeftMargin( leftMargin );
  c_ptggOvermgg->SetRightMargin( 1.6*rightMargin );
  c_ptggOvermgg->SetTopMargin( topMargin );
  c_ptggOvermgg->SetBottomMargin( bottomMargin );
  c_ptggOvermgg->SetFrameBorderMode(0);
  c_ptggOvermgg->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_ptggOvermgg_str = "hist_ptggOvermgg_"+Convert.str();
          hist_ptggOvermgg[i] = new TH1F(hist_ptggOvermgg_str.c_str(),"",30,0,2);
          std::string draw_ptggOvermgg_str = "pTGammaGamma/mGammaGamma>>"+hist_ptggOvermgg_str;
          if(i<1) tree[i]->Draw(draw_ptggOvermgg_str.c_str(),"( ( (mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.) ) && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 )","HISTSAME");
          else tree[i]->Draw(draw_ptggOvermgg_str.c_str(),"( ( (mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.) ) && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 )","*ESAME"); 
          hist_ptggOvermgg[i]->SetLineColor(i+2);
          hist_ptggOvermgg[i]->SetLineWidth(3);
          double norm = hist_ptggOvermgg[i]->Integral();
          hist_ptggOvermgg[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_ptggOvermgg[0]->SetTitle(";p_{T}^{#gamma #gamma}/m_{#gamma #gamma};Fraction of events");
  hist_ptggOvermgg[0]->GetYaxis()->SetRangeUser(0,0.25);
  hist_ptggOvermgg[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_ptggOvermgg[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.88,"data 2016");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.78,"data 2017");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_ptggOvermgg->Write();
  c_ptggOvermgg->Write();
  c_ptggOvermgg->Update();
  AddCMS(c_ptggOvermgg);
  c_ptggOvermgg->SaveAs("c_ptggOvermgg.pdf");
  c_ptggOvermgg->SaveAs("c_ptggOvermgg.png");
  c_ptggOvermgg->SaveAs("c_ptggOvermgg.C");
  
  //pho1Pt 
  TCanvas* c_pho1Pt = new TCanvas( "c_pho1Pt", "c_pho1Pt", 800, 700 );
  c_pho1Pt->SetHighLightColor(2);
  c_pho1Pt->SetFillColor(0);
  c_pho1Pt->SetBorderMode(0);
  c_pho1Pt->SetBorderSize(2);
  c_pho1Pt->SetLeftMargin( leftMargin );
  c_pho1Pt->SetRightMargin( 1.6*rightMargin );
  c_pho1Pt->SetTopMargin( topMargin );
  c_pho1Pt->SetBottomMargin( bottomMargin );
  c_pho1Pt->SetFrameBorderMode(0);
  c_pho1Pt->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho1Pt_str = "hist_pho1Pt_"+Convert.str();
          hist_pho1Pt[i] = new TH1F(hist_pho1Pt_str.c_str(),"",30,0,200);
          std::string draw_pho1Pt_str = "pho1Pt>>"+hist_pho1Pt_str;
          if(i<1) tree[i]->Draw(draw_pho1Pt_str.c_str(),"( ( (mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.) ) && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 )","HISTSAME");
          else tree[i]->Draw(draw_pho1Pt_str.c_str(),"( ( (mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.) ) && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 )","*ESAME"); 
          hist_pho1Pt[i]->SetLineColor(i+2);
          hist_pho1Pt[i]->SetLineWidth(3);
          double norm = hist_pho1Pt[i]->Integral();
          hist_pho1Pt[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho1Pt[0]->SetTitle(";p_{T}^{1} [GeV];Fraction of events");
  hist_pho1Pt[0]->GetYaxis()->SetRangeUser(0,0.2);
  hist_pho1Pt[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho1Pt[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.88,"data 2016");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.78,"data 2017");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho1Pt->Write();
  c_pho1Pt->Write();
  c_pho1Pt->Update();
  AddCMS(c_pho1Pt);
  c_pho1Pt->SaveAs("c_pho1Pt.pdf");
  c_pho1Pt->SaveAs("c_pho1Pt.png");
  c_pho1Pt->SaveAs("c_pho1Pt.C");

  //pho1Phi 
  TCanvas* c_pho1Phi = new TCanvas( "c_pho1Phi", "c_pho1Phi", 800, 700 );
  c_pho1Phi->SetHighLightColor(2);
  c_pho1Phi->SetFillColor(0);
  c_pho1Phi->SetBorderMode(0);
  c_pho1Phi->SetBorderSize(2);
  c_pho1Phi->SetLeftMargin( leftMargin );
  c_pho1Phi->SetRightMargin( 1.6*rightMargin );
  c_pho1Phi->SetTopMargin( topMargin );
  c_pho1Phi->SetBottomMargin( bottomMargin );
  c_pho1Phi->SetFrameBorderMode(0);
  c_pho1Phi->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho1Phi_str = "hist_pho1Phi_"+Convert.str();
          hist_pho1Phi[i] = new TH1F(hist_pho1Phi_str.c_str(),"",100,-3.,3.);
          std::string draw_pho1Phi_str = "pho1Phi>>"+hist_pho1Phi_str;
          if(i<1) tree[i]->Draw(draw_pho1Phi_str.c_str(),"(mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.)", "HISTSAME");
          else tree[i]->Draw(draw_pho1Phi_str.c_str(),"(mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.)", "*ESAME");
          hist_pho1Phi[i]->SetLineColor(i+2);
          hist_pho1Phi[i]->SetLineWidth(3);
          double norm = hist_pho1Phi[i]->Integral();
          hist_pho1Phi[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho1Phi[0]->SetTitle(";phi_{#gamma}^{2};Fraction of events");
  hist_pho1Phi[0]->GetYaxis()->SetRangeUser(0,0.015);
  hist_pho1Phi[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho1Phi[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.88,"data 2016");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.78,"data 2017");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho1Phi->Write();
  c_pho1Phi->Write();
  c_pho1Phi->Update();
  AddCMS(c_pho1Phi);
  c_pho1Phi->SaveAs("c_pho1Phi.pdf");
  c_pho1Phi->SaveAs("c_pho1Phi.png");
  c_pho1Phi->SaveAs("c_pho1Phi.C");
  
  //pho1Eta 
  TCanvas* c_pho1Eta = new TCanvas( "c_pho1Eta", "c_pho1Eta", 800, 700 );
  c_pho1Eta->SetHighLightColor(2);
  c_pho1Eta->SetFillColor(0);
  c_pho1Eta->SetBorderMode(0);
  c_pho1Eta->SetBorderSize(2);
  c_pho1Eta->SetLeftMargin( leftMargin );
  c_pho1Eta->SetRightMargin( 1.6*rightMargin );
  c_pho1Eta->SetTopMargin( topMargin );
  c_pho1Eta->SetBottomMargin( bottomMargin );
  c_pho1Eta->SetFrameBorderMode(0);
  c_pho1Eta->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho1Eta_str = "hist_pho1Eta_"+Convert.str();
          hist_pho1Eta[i] = new TH1F(hist_pho1Eta_str.c_str(),"",100,-3.,3.);
          std::string draw_pho1Eta_str = "pho1Eta>>"+hist_pho1Eta_str;
          if(i<1) tree[i]->Draw(draw_pho1Eta_str.c_str(),"(mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.)", "HISTSAME");
          else tree[i]->Draw(draw_pho1Eta_str.c_str(),"(mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.)", "*ESAME");
          hist_pho1Eta[i]->SetLineColor(i+2);
          hist_pho1Eta[i]->SetLineWidth(3);
          double norm = hist_pho1Eta[i]->Integral();
          hist_pho1Eta[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho1Eta[0]->SetTitle(";eta_{#gamma}^{2};Fraction of events");
  hist_pho1Eta[0]->GetYaxis()->SetRangeUser(0,0.03);
  hist_pho1Eta[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho1Eta[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.88,"data 2016");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.78,"data 2017");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho1Eta->Write();
  c_pho1Eta->Write();
  c_pho1Eta->Update();
  AddCMS(c_pho1Eta);
  c_pho1Eta->SaveAs("c_pho1Eta.pdf");
  c_pho1Eta->SaveAs("c_pho1Eta.png");
  c_pho1Eta->SaveAs("c_pho1Eta.C");
  
  //pho2Pt 
  TCanvas* c_pho2Pt = new TCanvas( "c_pho2Pt", "c_pho2Pt", 800, 700 );
  c_pho2Pt->SetHighLightColor(2);
  c_pho2Pt->SetFillColor(0);
  c_pho2Pt->SetBorderMode(0);
  c_pho2Pt->SetBorderSize(2);
  c_pho2Pt->SetLeftMargin( leftMargin );
  c_pho2Pt->SetRightMargin( 1.6*rightMargin );
  c_pho2Pt->SetTopMargin( topMargin );
  c_pho2Pt->SetBottomMargin( bottomMargin );
  c_pho2Pt->SetFrameBorderMode(0);
  c_pho2Pt->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho2Pt_str = "hist_pho2Pt_"+Convert.str();
          hist_pho2Pt[i] = new TH1F(hist_pho2Pt_str.c_str(),"",30,0,200);
          std::string draw_pho2Pt_str = "pho2Pt>>"+hist_pho2Pt_str;
          if(i<1) tree[i]->Draw(draw_pho2Pt_str.c_str(),"( ( (mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.) ) && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 )","HISTSAME");
          else tree[i]->Draw(draw_pho2Pt_str.c_str(),"( ( (mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.) ) && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 )","*ESAME"); 
          hist_pho2Pt[i]->SetLineColor(i+2);
          hist_pho2Pt[i]->SetLineWidth(3);
          double norm = hist_pho2Pt[i]->Integral();
          hist_pho2Pt[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho2Pt[0]->SetTitle(";p_{T}^{2} [GeV];Fraction of events");
  hist_pho2Pt[0]->GetYaxis()->SetRangeUser(0,0.3);
  hist_pho2Pt[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho2Pt[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.88,"data 2016");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.78,"data 2017");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho2Pt->Write();
  c_pho2Pt->Write();
  c_pho2Pt->Update();
  AddCMS(c_pho2Pt);
  c_pho2Pt->SaveAs("c_pho2Pt.pdf");
  c_pho2Pt->SaveAs("c_pho2Pt.png");
  c_pho2Pt->SaveAs("c_pho2Pt.C");


  //pho2Phi 
  TCanvas* c_pho2Phi = new TCanvas( "c_pho2Phi", "c_pho2Phi", 800, 700 );
  c_pho2Phi->SetHighLightColor(2);
  c_pho2Phi->SetFillColor(0);
  c_pho2Phi->SetBorderMode(0);
  c_pho2Phi->SetBorderSize(2);
  c_pho2Phi->SetLeftMargin( leftMargin );
  c_pho2Phi->SetRightMargin( 1.6*rightMargin );
  c_pho2Phi->SetTopMargin( topMargin );
  c_pho2Phi->SetBottomMargin( bottomMargin );
  c_pho2Phi->SetFrameBorderMode(0);
  c_pho2Phi->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho2Phi_str = "hist_pho2Phi_"+Convert.str();
          hist_pho2Phi[i] = new TH1F(hist_pho2Phi_str.c_str(),"",100,-3.,3.);
          std::string draw_pho2Phi_str = "pho2Phi>>"+hist_pho2Phi_str;
          if(i<1) tree[i]->Draw(draw_pho2Phi_str.c_str(),"(mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.)", "HISTSAME");
          else tree[i]->Draw(draw_pho2Phi_str.c_str(),"(mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.)", "*ESAME");
          hist_pho2Phi[i]->SetLineColor(i+2);
          hist_pho2Phi[i]->SetLineWidth(3);
          double norm = hist_pho2Phi[i]->Integral();
          hist_pho2Phi[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho2Phi[0]->SetTitle(";phi_{#gamma}^{2};Fraction of events");
  hist_pho2Phi[0]->GetYaxis()->SetRangeUser(0,0.015);
  hist_pho2Phi[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho2Phi[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.88,"data 2016");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.78,"data 2017");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho2Phi->Write();
  c_pho2Phi->Write();
  c_pho2Phi->Update();
  AddCMS(c_pho2Phi);
  c_pho2Phi->SaveAs("c_pho2Phi.pdf");
  c_pho2Phi->SaveAs("c_pho2Phi.png");
  c_pho2Phi->SaveAs("c_pho2Phi.C");
  
  //pho2Eta 
  TCanvas* c_pho2Eta = new TCanvas( "c_pho2Eta", "c_pho2Eta", 800, 700 );
  c_pho2Eta->SetHighLightColor(2);
  c_pho2Eta->SetFillColor(0);
  c_pho2Eta->SetBorderMode(0);
  c_pho2Eta->SetBorderSize(2);
  c_pho2Eta->SetLeftMargin( leftMargin );
  c_pho2Eta->SetRightMargin( 1.6*rightMargin );
  c_pho2Eta->SetTopMargin( topMargin );
  c_pho2Eta->SetBottomMargin( bottomMargin );
  c_pho2Eta->SetFrameBorderMode(0);
  c_pho2Eta->SetFrameBorderMode(0);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho2Eta_str = "hist_pho2Eta_"+Convert.str();
          hist_pho2Eta[i] = new TH1F(hist_pho2Eta_str.c_str(),"",100,-3.,3.);
          std::string draw_pho2Eta_str = "pho2Eta>>"+hist_pho2Eta_str;
          if(i<1) tree[i]->Draw(draw_pho2Eta_str.c_str(),"(mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.)", "HISTSAME");
          else tree[i]->Draw(draw_pho2Eta_str.c_str(),"(mGammaGamma>105. && mGammaGamma<120.) || (mGammaGamma>130. && mGammaGamma<180.)", "*ESAME");
          hist_pho2Eta[i]->SetLineColor(i+2);
          hist_pho2Eta[i]->SetLineWidth(3);
          double norm = hist_pho2Eta[i]->Integral();
          hist_pho2Eta[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho2Eta[0]->SetTitle(";eta_{#gamma}^{2};Fraction of events");
  hist_pho2Eta[0]->GetYaxis()->SetRangeUser(0,0.03);
  hist_pho2Eta[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho2Eta[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.88,"data 2016");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.78,"data 2017");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho2Eta->Write();
  c_pho2Eta->Write();
  c_pho2Eta->Update();
  AddCMS(c_pho2Eta);
  c_pho2Eta->SaveAs("c_pho2Eta.pdf");
  c_pho2Eta->SaveAs("c_pho2Eta.png");
  c_pho2Eta->SaveAs("c_pho2Eta.C");
  
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
