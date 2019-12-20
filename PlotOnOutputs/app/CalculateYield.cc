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
//TString lumiText = "Simulation (13 TeV)";
TString lumiText = "Weighted (13 TeV)";

bool AddCMS( TCanvas* C );

using namespace std;

int main( int argc, char** argv )
{
  ifstream file;// read input list directly
  std::string line;
  int n = 0;
  std::cout << "[Usage]: ./CalculateYield inputList output.root\n" << std::endl;
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
  
  TH1F* hist_mGammaGamma[n]; //reconstructed Higgs mass

  //float lumi = 1000.; // 1/fb 
  
  //mGammaGamma leptonic
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
  TLegend *leg = new TLegend(0.71,0.74,0.91,0.94);
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_mGammaGamma_str = "hist_mGammaGamma_"+Convert.str();
          hist_mGammaGamma[i] = new TH1F(hist_mGammaGamma_str.c_str(),"",50,100,180);
          std::string draw_mGammaGamma_str = "mGammaGamma>>"+hist_mGammaGamma_str;
          //tree[i]->Draw(draw_mGammaGamma_str.c_str(),"( (mGammaGamma>100. && mGammaGamma<180.) && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && pho1R9>0.5 && pho2R9>0.5 && n_Jets>=1 ) ","HISTSAME");
          tree[i]->Draw(draw_mGammaGamma_str.c_str(),"weight*( (mGammaGamma>100. && mGammaGamma<180.) && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && pho1R9>0.5 && pho2R9>0.5 && n_Jets>=1 ) ","HISTSAME");
          if(i<3) hist_mGammaGamma[i]->SetLineColor(i+2);
          else hist_mGammaGamma[i]->SetLineColor(i+3);
          hist_mGammaGamma[i]->SetLineWidth(3);
          //double norm = hist_mGammaGamma[i]->Integral();
          //hist_mGammaGamma[i]->Scale(1/norm);
          leg->AddEntry(hist_mGammaGamma[i],labelProcess[i].c_str(),"lp");
          std::cout << "i = " << i << "\n" << std::endl;
          std::cout << "Yield of "<< labelProcess[i].c_str() << " : " << (hist_mGammaGamma[i]->GetSum())*1000. << "\n" << std::endl;
  }
  leg->Draw("SAME");
  hist_mGammaGamma[0]->SetTitle(";m_{#gamma #gamma} [GeV];# events");
  //hist_mGammaGamma[0]->GetYaxis()->SetRangeUser(0,0.36);
  hist_mGammaGamma[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_mGammaGamma[0]->GetXaxis()->SetTitleOffset(1.4);  
  c_mGammaGamma->Write();
  c_mGammaGamma->Write();
  c_mGammaGamma->Update();
  AddCMS(c_mGammaGamma);
  c_mGammaGamma->SaveAs("c_mGammaGamma.pdf");
  c_mGammaGamma->SaveAs("c_mGammaGamma.png");
  c_mGammaGamma->SaveAs("c_mGammaGamma.C");

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
