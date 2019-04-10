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
//TString lumiText = "36 /fb (13 TeV)";
TString lumiText = "Simulation (13 TeV)";

bool AddCMS( TCanvas* C );

using namespace std;

int main( int argc, char** argv )
{
  ifstream file;// read input list directly
  std::string line;
  int n = 0;
  std::cout << "[Usage]: ./MetRecipe inputList outputFile\n" << std::endl;
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
          if(process=="original") tree[i] = (TTree*)fin[i]->Get("HggRazorLeptons");
          if(process=="corrected") tree[i] = (TTree*)fin[i]->Get("HggRazorLeptons");

          /*
          {
                  TDirectory *dir = (TDirectory*)f->Get("/afs/cern.ch/user/m/mschoene/public/data_2016_sync_sep06.root:/diPhoton_data/HT0toInf_j0toInf_b0toInf");
                  tree[i] = dir->GetObject("tree_diPhoton_data_HT0toInf_j0toInf_b0toInf",tree);
          }
          */

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


  //numbers 
  TH1F* hist_MR[m];
  TH1F* hist_R2[m];
  TH1F* hist_t1MET[m];

  TH1F* hist_ratio;
  
  TString cut = "mGammaGamma > 76. && mGammaGamma < 106. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 0 && pho2passEleVeto == 0 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && (  HLTDecision[22] || HLTDecision[23] || HLTDecision[24]  )";
  //TString cut = "mGammaGamma > 76. && mGammaGamma < 106. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 0 && pho2passEleVeto == 0 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && (Flag_HBHENoiseFilter == 1 && Flag_goodVertices == 1 && Flag_eeBadScFilter == 1 && Flag_HBHEIsoNoiseFilter == 1 && Flag_CSCTightHaloFilter == 1 ) && (  HLTDecision[22] || HLTDecision[23] || HLTDecision[24]  )";

  //Mr numbers 
  TCanvas* c_MR = new TCanvas( "c_MR", "c_MR", 800, 700 );
  c_MR->SetHighLightColor(2);
  c_MR->SetFillColor(0);
  c_MR->SetBorderMode(0);
  c_MR->SetBorderSize(2);
  c_MR->SetLeftMargin( leftMargin );
  c_MR->SetRightMargin( 1.6*rightMargin );
  c_MR->SetTopMargin( topMargin );
  c_MR->SetBottomMargin( bottomMargin );
  c_MR->SetFrameBorderMode(0);
  c_MR->SetFrameBorderMode(0);
  TPad* pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_MR_str = "hist_MR_"+Convert.str();
          hist_MR[i] = new TH1F(hist_MR_str.c_str(),"",100,0,1500);
                  std::string draw_MR_str = "MR>>"+hist_MR_str;
                  tree[i]->Draw(draw_MR_str.c_str(),cut,"HISTSAME");
          hist_MR[i]->SetLineColor(i+2);
          hist_MR[i]->SetLineWidth(3);
          double norm = hist_MR[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_MR[0]->SetTitle(";MR [GeV];# of events");
  hist_MR[0]->GetYaxis()->SetRangeUser(0,2e5);
  hist_MR[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_MR[0]->GetXaxis()->SetTitleOffset(1.4);
 TLatex* tex = new TLatex(0.27,0.78,"original");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"recipe applied");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_MR->cd();
  TPad* pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_MR[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_MR[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("original/corr ");
  hist_ratio->GetYaxis()->SetNdivisions(505);
  hist_ratio->GetYaxis()->SetTitleSize(20);
  hist_ratio->GetYaxis()->SetTitleFont(43);
  hist_ratio->GetYaxis()->SetTitleOffset(1.55);
  hist_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hist_ratio->GetYaxis()->SetLabelSize(15);
  hist_ratio->GetXaxis()->SetTitleSize(20);
  hist_ratio->GetXaxis()->SetTitleFont(43);
  hist_ratio->GetXaxis()->SetTitleOffset(4.);
  hist_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hist_ratio->GetXaxis()->SetLabelSize(15);
  c_MR->cd();
  c_MR->Write();
  c_MR->Write();
  c_MR->Update();
  AddCMS(c_MR);
  c_MR->SaveAs("c_MR_metrecipe.pdf");
  c_MR->SaveAs("c_MR_metrecipe.png");
  c_MR->SaveAs("c_MR_metrecipe.C");

  //R2 numbers 
  TCanvas* c_R2 = new TCanvas( "c_R2", "c_R2", 800, 700 );
  c_R2->SetHighLightColor(2);
  c_R2->SetFillColor(0);
  c_R2->SetBorderMode(0);
  c_R2->SetBorderSize(2);
  c_R2->SetLeftMargin( leftMargin );
  c_R2->SetRightMargin( 1.6*rightMargin );
  c_R2->SetTopMargin( topMargin );
  c_R2->SetBottomMargin( bottomMargin );
  c_R2->SetFrameBorderMode(0);
  c_R2->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_R2_str = "hist_R2_"+Convert.str();
          hist_R2[i] = new TH1F(hist_R2_str.c_str(),"",100,0,1.5);
                  std::string draw_R2_str = "t1Rsq>>"+hist_R2_str;
                  tree[i]->Draw(draw_R2_str.c_str(),cut,"HISTSAME");
          hist_R2[i]->SetLineColor(i+2);
          hist_R2[i]->SetLineWidth(3);
          double norm = hist_R2[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_R2[0]->SetTitle(";Rsq;# of events");
  hist_R2[0]->GetYaxis()->SetRangeUser(0,6e4);
  hist_R2[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_R2[0]->GetXaxis()->SetTitleOffset(1.4);
 tex = new TLatex(0.27,0.78,"original");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"recipe applied");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_R2->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_R2[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_R2[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("original/corr ");
  hist_ratio->GetYaxis()->SetNdivisions(505);
  hist_ratio->GetYaxis()->SetTitleSize(20);
  hist_ratio->GetYaxis()->SetTitleFont(43);
  hist_ratio->GetYaxis()->SetTitleOffset(1.55);
  hist_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hist_ratio->GetYaxis()->SetLabelSize(15);
  hist_ratio->GetXaxis()->SetTitleSize(20);
  hist_ratio->GetXaxis()->SetTitleFont(43);
  hist_ratio->GetXaxis()->SetTitleOffset(4.);
  hist_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hist_ratio->GetXaxis()->SetLabelSize(15);
  c_R2->cd();
  c_R2->Write();
  c_R2->Write();
  c_R2->Update();
  AddCMS(c_R2);
  c_R2->SaveAs("c_R2_metrecipe.pdf");
  c_R2->SaveAs("c_R2_metrecipe.png");
  c_R2->SaveAs("c_R2_metrecipe.C");

  //t1MET numbers 
  TCanvas* c_t1MET = new TCanvas( "c_t1MET", "c_t1MET", 800, 700 );
  c_t1MET->SetHighLightColor(2);
  c_t1MET->SetFillColor(0);
  c_t1MET->SetBorderMode(0);
  c_t1MET->SetBorderSize(2);
  c_t1MET->SetLeftMargin( leftMargin );
  c_t1MET->SetRightMargin( 1.6*rightMargin );
  c_t1MET->SetTopMargin( topMargin );
  c_t1MET->SetBottomMargin( bottomMargin );
  c_t1MET->SetFrameBorderMode(0);
  c_t1MET->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_t1MET_str = "hist_t1MET_"+Convert.str();
          hist_t1MET[i] = new TH1F(hist_t1MET_str.c_str(),"",60,0,300);
                  std::string draw_t1MET_str = "t1MET>>"+hist_t1MET_str;
                  tree[i]->Draw(draw_t1MET_str.c_str(),cut,"HISTSAME");
          hist_t1MET[i]->SetLineColor(i+2);
          hist_t1MET[i]->SetLineWidth(3);
          double norm = hist_t1MET[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_t1MET[0]->SetTitle(";t1MET (GeV);# of events");
  hist_t1MET[0]->GetYaxis()->SetRangeUser(0,1e5);
  hist_t1MET[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_t1MET[0]->GetXaxis()->SetTitleOffset(1.4);
 tex = new TLatex(0.27,0.78,"original");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"recipe applied");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_t1MET->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_t1MET[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0);
  hist_ratio->SetMaximum(5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_t1MET[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("original/corr ");
  hist_ratio->GetYaxis()->SetNdivisions(505);
  hist_ratio->GetYaxis()->SetTitleSize(20);
  hist_ratio->GetYaxis()->SetTitleFont(43);
  hist_ratio->GetYaxis()->SetTitleOffset(1.55);
  hist_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hist_ratio->GetYaxis()->SetLabelSize(15);
  hist_ratio->GetXaxis()->SetTitleSize(20);
  hist_ratio->GetXaxis()->SetTitleFont(43);
  hist_ratio->GetXaxis()->SetTitleOffset(4.);
  hist_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hist_ratio->GetXaxis()->SetLabelSize(15);
  c_t1MET->cd();
  c_t1MET->Write();
  c_t1MET->Write();
  c_t1MET->Update();
  AddCMS(c_t1MET);
  c_t1MET->SaveAs("c_t1MET_metrecipe.pdf");
  c_t1MET->SaveAs("c_t1MET_metrecipe.png");
  c_t1MET->SaveAs("c_t1MET_metrecipe.C");


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
