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
TString lumiText = "36 /fb (13 TeV)";
//TString lumiText = "Simulation (13 TeV)";

bool AddCMS( TCanvas* C );

using namespace std;

int main( int argc, char** argv )
{
  ifstream file;// read input list directly
  std::string line;
  int n = 0;
  std::cout << "[Usage]: ./RazorVsMt2 inputList outputFile\n" << std::endl;
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
          if(process=="Razor") tree[i] = (TTree*)fin[i]->Get("HggRazorLeptons");
          if(process=="Mt2") 
                  tree[i] = (TTree*)fin[i]->Get("diPhoton_data/HT0toInf_j0toInf_b0toInf/tree_diPhoton_data_HT0toInf_j0toInf_b0toInf");
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

  //fractions
  TH1F* hist_pTGammaGamma[m]; //reconstructed ptgg
  TH1F* hist_mGammaGamma[m]; //reconstructed mgg
  TH1F* hist_ptggOvermgg[m]; //reconstructed mgg
  TH1F* hist_pho1Pt[m]; 
  TH1F* hist_pho1Eta[m]; 
  TH1F* hist_pho1Phi[m]; 
  TH1F* hist_pho1R9[m]; 
  TH1F* hist_pho1SigmaIetaIeta[m]; 
  TH1F* hist_pho2Pt[m]; 
  TH1F* hist_pho2Eta[m]; 
  TH1F* hist_pho2Phi[m]; 
  TH1F* hist_pho2R9[m]; 
  TH1F* hist_pho2SigmaIetaIeta[m];

  //numbers 
  TH1F* hist_pTGG[m]; //reconstructed ptgg
  TH1F* hist_mGG[m]; //reconstructed mgg
  TH1F* hist_G1Pt[m]; 
  TH1F* hist_G1Eta[m]; 
  TH1F* hist_G1Phi[m]; 
  TH1F* hist_G1R9[m]; 
  TH1F* hist_G1SigmaIetaIeta[m]; 
  TH1F* hist_G2Pt[m]; 
  TH1F* hist_G2Eta[m]; 
  TH1F* hist_G2Phi[m]; 
  TH1F* hist_G2R9[m]; 
  TH1F* hist_G2SigmaIetaIeta[m]; 

  TH1F* hist_ratio;


  //mGammaGamma fraction 
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
  c_mGammaGamma->cd();
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_mGammaGamma_str = "hist_mGammaGamma_"+Convert.str();
          hist_mGammaGamma[i] = new TH1F(hist_mGammaGamma_str.c_str(),"",30,103,160);
          if(i<1)
          { 
                  std::string draw_mGammaGamma_str = "mGammaGamma>>"+hist_mGammaGamma_str;
                  tree[i]->Draw(draw_mGammaGamma_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_mGammaGamma_str = "h_mass>>"+hist_mGammaGamma_str;
                  tree[i]->Draw(draw_mGammaGamma_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
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
   TLatex * tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_mGammaGamma->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->SetBottomMargin(0);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_mGammaGamma[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_mGammaGamma[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_mGammaGamma->cd();
  c_mGammaGamma->Write();
  c_mGammaGamma->Write();
  c_mGammaGamma->Update();
  AddCMS(c_mGammaGamma);
  c_mGammaGamma->SaveAs("c_mGammaGamma_com.pdf");
  c_mGammaGamma->SaveAs("c_mGammaGamma_com.png");
  c_mGammaGamma->SaveAs("c_mGammaGamma_com.C");

  //mGammaGamma numbers 
  TCanvas* c_mGG = new TCanvas( "c_mGG", "c_mGG", 800, 700 );
  c_mGG->SetHighLightColor(2);
  c_mGG->SetFillColor(0);
  c_mGG->SetBorderMode(0);
  c_mGG->SetBorderSize(2);
  c_mGG->SetLeftMargin( leftMargin );
  c_mGG->SetRightMargin( 1.6*rightMargin );
  c_mGG->SetTopMargin( topMargin );
  c_mGG->SetBottomMargin( bottomMargin );
  c_mGG->SetFrameBorderMode(0);
  c_mGG->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_mGG_str = "hist_mGG_"+Convert.str();
          hist_mGG[i] = new TH1F(hist_mGG_str.c_str(),"",30,103,160);
          if(i<1)
          { 
                  std::string draw_mGG_str = "mGammaGamma>>"+hist_mGG_str;
                  tree[i]->Draw(draw_mGG_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_mGG_str = "h_mass>>"+hist_mGG_str;
                  tree[i]->Draw(draw_mGG_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_mGG[i]->SetLineColor(i+2);
          hist_mGG[i]->SetLineWidth(3);
          double norm = hist_mGG[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_mGG[0]->SetTitle(";m_{#gamma #gamma} [GeV];# of events");
  hist_mGG[0]->GetYaxis()->SetRangeUser(0,27000);
  hist_mGG[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_mGG[0]->GetXaxis()->SetTitleOffset(1.4);
 tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_mGG->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_mGG[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_mGG[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_mGG->cd();
  c_mGG->Write();
  c_mGG->Write();
  c_mGG->Update();
  AddCMS(c_mGG);
  c_mGG->SaveAs("c_mGG_com.pdf");
  c_mGG->SaveAs("c_mGG_com.png");
  c_mGG->SaveAs("c_mGG_com.C");

  //pTGammaGamma fraction 
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
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pTGammaGamma_str = "hist_pTGammaGamma_"+Convert.str();
          hist_pTGammaGamma[i] = new TH1F(hist_pTGammaGamma_str.c_str(),"",30,0,200);
          if(i<1)
          { 
                  std::string draw_pTGammaGamma_str = "pTGammaGamma>>"+hist_pTGammaGamma_str;
                  tree[i]->Draw(draw_pTGammaGamma_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_pTGammaGamma_str = "h_pt>>"+hist_pTGammaGamma_str;
                  tree[i]->Draw(draw_pTGammaGamma_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_pTGammaGamma[i]->SetLineColor(i+2);
          hist_pTGammaGamma[i]->SetLineWidth(3);
          double norm = hist_pTGammaGamma[i]->Integral();
          hist_pTGammaGamma[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pTGammaGamma[0]->SetTitle(";pT_{#gamma #gamma} [GeV];Fraction of events");
  hist_pTGammaGamma[0]->GetYaxis()->SetRangeUser(0,0.2);
  hist_pTGammaGamma[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pTGammaGamma[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pTGammaGamma->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_pTGammaGamma[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_pTGammaGamma[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_pTGammaGamma->cd();
  c_pTGammaGamma->Write();
  c_pTGammaGamma->Write();
  c_pTGammaGamma->Update();
  AddCMS(c_pTGammaGamma);
  c_pTGammaGamma->SaveAs("c_pTGammaGamma_com.pdf");
  c_pTGammaGamma->SaveAs("c_pTGammaGamma_com.png");
  c_pTGammaGamma->SaveAs("c_pTGammaGamma_com.C");

  //pTGammaGamma numbers 
  TCanvas* c_pTGG = new TCanvas( "c_pTGG", "c_pTGG", 800, 700 );
  c_pTGG->SetHighLightColor(2);
  c_pTGG->SetFillColor(0);
  c_pTGG->SetBorderMode(0);
  c_pTGG->SetBorderSize(2);
  c_pTGG->SetLeftMargin( leftMargin );
  c_pTGG->SetRightMargin( 1.6*rightMargin );
  c_pTGG->SetTopMargin( topMargin );
  c_pTGG->SetBottomMargin( bottomMargin );
  c_pTGG->SetFrameBorderMode(0);
  c_pTGG->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pTGG_str = "hist_pTGG_"+Convert.str();
          hist_pTGG[i] = new TH1F(hist_pTGG_str.c_str(),"",30,0,200);
          if(i<1)
          { 
                  std::string draw_pTGG_str = "pTGammaGamma>>"+hist_pTGG_str;
                  tree[i]->Draw(draw_pTGG_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_pTGG_str = "h_pt>>"+hist_pTGG_str;
                  tree[i]->Draw(draw_pTGG_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_pTGG[i]->SetLineColor(i+2);
          hist_pTGG[i]->SetLineWidth(3);
          double norm = hist_pTGG[i]->Integral();
          //hist_pTGG[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pTGG[0]->SetTitle(";pT_{#gamma #gamma} [GeV];# of events");
  hist_pTGG[0]->GetYaxis()->SetRangeUser(0,60000);
  hist_pTGG[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pTGG[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pTGG->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_pTGG[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_pTGG[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_pTGG->cd();
  c_pTGG->Write();
  c_pTGG->Write();
  c_pTGG->Update();
  AddCMS(c_pTGG);
  c_pTGG->SaveAs("c_pTGG_com.pdf");
  c_pTGG->SaveAs("c_pTGG_com.png");
  c_pTGG->SaveAs("c_pTGG_com.C");

  //pho1Pt fraction 
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
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho1Pt_str = "hist_pho1Pt_"+Convert.str();
          hist_pho1Pt[i] = new TH1F(hist_pho1Pt_str.c_str(),"",30,0,200);
          if(i<1)
          { 
                  std::string draw_pho1Pt_str = "pho1Pt>>"+hist_pho1Pt_str;
                  tree[i]->Draw(draw_pho1Pt_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_pho1Pt_str = "ptGamma0>>"+hist_pho1Pt_str;
                  tree[i]->Draw(draw_pho1Pt_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_pho1Pt[i]->SetLineColor(i+2);
          hist_pho1Pt[i]->SetLineWidth(3);
          double norm = hist_pho1Pt[i]->Integral();
          hist_pho1Pt[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho1Pt[0]->SetTitle(";Pt_{#gamma 1} [GeV];Fraction of events");
  hist_pho1Pt[0]->GetYaxis()->SetRangeUser(0,0.2);
  hist_pho1Pt[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho1Pt[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho1Pt->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_pho1Pt[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_pho1Pt[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_pho1Pt->cd();
  c_pho1Pt->Write();
  c_pho1Pt->Write();
  c_pho1Pt->Update();
  AddCMS(c_pho1Pt);
  c_pho1Pt->SaveAs("c_pho1Pt_com.pdf");
  c_pho1Pt->SaveAs("c_pho1Pt_com.png");
  c_pho1Pt->SaveAs("c_pho1Pt_com.C");

  //pho1Pt numbers 
  TCanvas* c_G1Pt = new TCanvas( "c_G1Pt", "c_G1Pt", 800, 700 );
  c_G1Pt->SetHighLightColor(2);
  c_G1Pt->SetFillColor(0);
  c_G1Pt->SetBorderMode(0);
  c_G1Pt->SetBorderSize(2);
  c_G1Pt->SetLeftMargin( leftMargin );
  c_G1Pt->SetRightMargin( 1.6*rightMargin );
  c_G1Pt->SetTopMargin( topMargin );
  c_G1Pt->SetBottomMargin( bottomMargin );
  c_G1Pt->SetFrameBorderMode(0);
  c_G1Pt->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_G1Pt_str = "hist_G1Pt_"+Convert.str();
          hist_G1Pt[i] = new TH1F(hist_G1Pt_str.c_str(),"",30,0,200);
          if(i<1)
          { 
                  std::string draw_G1Pt_str = "pho1Pt>>"+hist_G1Pt_str;
                  tree[i]->Draw(draw_G1Pt_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_G1Pt_str = "ptGamma0>>"+hist_G1Pt_str;
                  tree[i]->Draw(draw_G1Pt_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_G1Pt[i]->SetLineColor(i+2);
          hist_G1Pt[i]->SetLineWidth(3);
          double norm = hist_G1Pt[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_G1Pt[0]->SetTitle(";Pt_{#gamma 1} [GeV];# of events");
  hist_G1Pt[0]->GetYaxis()->SetRangeUser(0,60000);
  hist_G1Pt[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_G1Pt[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_G1Pt->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_G1Pt[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_G1Pt[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_G1Pt->cd();
  c_G1Pt->Write();
  c_G1Pt->Write();
  c_G1Pt->Update();
  AddCMS(c_G1Pt);
  c_G1Pt->SaveAs("c_G1Pt_com.pdf");
  c_G1Pt->SaveAs("c_G1Pt_com.png");
  c_G1Pt->SaveAs("c_G1Pt_com.C");

  //pho2Pt fraction 
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
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho2Pt_str = "hist_pho2Pt_"+Convert.str();
          hist_pho2Pt[i] = new TH1F(hist_pho2Pt_str.c_str(),"",30,0,200);
          if(i<1)
          { 
                  std::string draw_pho2Pt_str = "pho2Pt>>"+hist_pho2Pt_str;
                  tree[i]->Draw(draw_pho2Pt_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_pho2Pt_str = "ptGamma1>>"+hist_pho2Pt_str;
                  tree[i]->Draw(draw_pho2Pt_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_pho2Pt[i]->SetLineColor(i+2);
          hist_pho2Pt[i]->SetLineWidth(3);
          double norm = hist_pho2Pt[i]->Integral();
          hist_pho2Pt[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho2Pt[0]->SetTitle(";Pt_{#gamma 2} [GeV];Fraction of events");
  hist_pho2Pt[0]->GetYaxis()->SetRangeUser(0,0.3);
  hist_pho2Pt[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho2Pt[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho2Pt->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_pho2Pt[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_pho2Pt[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_pho2Pt->cd();
  c_pho2Pt->Write();
  c_pho2Pt->Write();
  c_pho2Pt->Update();
  AddCMS(c_pho2Pt);
  c_pho2Pt->SaveAs("c_pho2Pt_com.pdf");
  c_pho2Pt->SaveAs("c_pho2Pt_com.png");
  c_pho2Pt->SaveAs("c_pho2Pt_com.C");

  //pho2Pt numbers 
  TCanvas* c_G2Pt = new TCanvas( "c_G2Pt", "c_G2Pt", 800, 700 );
  c_G2Pt->SetHighLightColor(2);
  c_G2Pt->SetFillColor(0);
  c_G2Pt->SetBorderMode(0);
  c_G2Pt->SetBorderSize(2);
  c_G2Pt->SetLeftMargin( leftMargin );
  c_G2Pt->SetRightMargin( 1.6*rightMargin );
  c_G2Pt->SetTopMargin( topMargin );
  c_G2Pt->SetBottomMargin( bottomMargin );
  c_G2Pt->SetFrameBorderMode(0);
  c_G2Pt->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_G2Pt_str = "hist_G2Pt_"+Convert.str();
          hist_G2Pt[i] = new TH1F(hist_G2Pt_str.c_str(),"",30,0,200);
          if(i<1)
          { 
                  std::string draw_G2Pt_str = "pho2Pt>>"+hist_G2Pt_str;
                  tree[i]->Draw(draw_G2Pt_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_G2Pt_str = "ptGamma1>>"+hist_G2Pt_str;
                  tree[i]->Draw(draw_G2Pt_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_G2Pt[i]->SetLineColor(i+2);
          hist_G2Pt[i]->SetLineWidth(3);
          double norm = hist_G2Pt[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_G2Pt[0]->SetTitle(";Pt_{#gamma 2} [GeV];# of events");
  hist_G2Pt[0]->GetYaxis()->SetRangeUser(0,80000);
  hist_G2Pt[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_G2Pt[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_G2Pt->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_G2Pt[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_G2Pt[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_G2Pt->cd();
  c_G2Pt->Write();
  c_G2Pt->Write();
  c_G2Pt->Update();
  AddCMS(c_G2Pt);
  c_G2Pt->SaveAs("c_G2Pt_com.pdf");
  c_G2Pt->SaveAs("c_G2Pt_com.png");
  c_G2Pt->SaveAs("c_G2Pt_com.C");

  //pho1Phi fraction 
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
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho1Phi_str = "hist_pho1Phi_"+Convert.str();
          hist_pho1Phi[i] = new TH1F(hist_pho1Phi_str.c_str(),"",30,-3,3);
          if(i<1)
          { 
                  std::string draw_pho1Phi_str = "pho1Phi>>"+hist_pho1Phi_str;
                  tree[i]->Draw(draw_pho1Phi_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_pho1Phi_str = "phiGamma0>>"+hist_pho1Phi_str;
                  tree[i]->Draw(draw_pho1Phi_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_pho1Phi[i]->SetLineColor(i+2);
          hist_pho1Phi[i]->SetLineWidth(3);
          double norm = hist_pho1Phi[i]->Integral();
          hist_pho1Phi[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho1Phi[0]->SetTitle(";Phi_{#gamma 1} [GeV];Fraction of events");
  hist_pho1Phi[0]->GetYaxis()->SetRangeUser(0,0.13);
  hist_pho1Phi[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho1Phi[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho1Phi->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_pho1Phi[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_pho1Phi[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_pho1Phi->cd();
  c_pho1Phi->Write();
  c_pho1Phi->Write();
  c_pho1Phi->Update();
  AddCMS(c_pho1Phi);
  c_pho1Phi->SaveAs("c_pho1Phi_com.pdf");
  c_pho1Phi->SaveAs("c_pho1Phi_com.png");
  c_pho1Phi->SaveAs("c_pho1Phi_com.C");

  //pho1Phi numbers 
  TCanvas* c_G1Phi = new TCanvas( "c_G1Phi", "c_G1Phi", 800, 700 );
  c_G1Phi->SetHighLightColor(2);
  c_G1Phi->SetFillColor(0);
  c_G1Phi->SetBorderMode(0);
  c_G1Phi->SetBorderSize(2);
  c_G1Phi->SetLeftMargin( leftMargin );
  c_G1Phi->SetRightMargin( 1.6*rightMargin );
  c_G1Phi->SetTopMargin( topMargin );
  c_G1Phi->SetBottomMargin( bottomMargin );
  c_G1Phi->SetFrameBorderMode(0);
  c_G1Phi->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_G1Phi_str = "hist_G1Phi_"+Convert.str();
          hist_G1Phi[i] = new TH1F(hist_G1Phi_str.c_str(),"",30,-3,3);
          if(i<1)
          { 
                  std::string draw_G1Phi_str = "pho1Phi>>"+hist_G1Phi_str;
                  tree[i]->Draw(draw_G1Phi_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_G1Phi_str = "phiGamma0>>"+hist_G1Phi_str;
                  tree[i]->Draw(draw_G1Phi_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_G1Phi[i]->SetLineColor(i+2);
          hist_G1Phi[i]->SetLineWidth(3);
          double norm = hist_G1Phi[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_G1Phi[0]->SetTitle(";Phi_{#gamma 1} [GeV];# of events");
  hist_G1Phi[0]->GetYaxis()->SetRangeUser(0,10000);
  hist_G1Phi[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_G1Phi[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_G1Phi->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_G1Phi[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_G1Phi[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_G1Phi->cd();
  c_G1Phi->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_G1Phi[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_G1Phi[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_G1Phi->cd();
  c_G1Phi->Write();
  c_G1Phi->Write();
  c_G1Phi->Update();
  AddCMS(c_G1Phi);
  c_G1Phi->SaveAs("c_G1Phi_com.pdf");
  c_G1Phi->SaveAs("c_G1Phi_com.png");
  c_G1Phi->SaveAs("c_G1Phi_com.C");

  //pho2Phi fraction 
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
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho2Phi_str = "hist_pho2Phi_"+Convert.str();
          hist_pho2Phi[i] = new TH1F(hist_pho2Phi_str.c_str(),"",30,-3,3);
          if(i<1)
          { 
                  std::string draw_pho2Phi_str = "pho2Phi>>"+hist_pho2Phi_str;
                  tree[i]->Draw(draw_pho2Phi_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_pho2Phi_str = "phiGamma1>>"+hist_pho2Phi_str;
                  tree[i]->Draw(draw_pho2Phi_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_pho2Phi[i]->SetLineColor(i+2);
          hist_pho2Phi[i]->SetLineWidth(3);
          double norm = hist_pho2Phi[i]->Integral();
          hist_pho2Phi[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho2Phi[0]->SetTitle(";Phi_{#gamma 2} [GeV];Fraction of events");
  hist_pho2Phi[0]->GetYaxis()->SetRangeUser(0,0.13);
  hist_pho2Phi[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho2Phi[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho2Phi->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_pho2Phi[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_pho2Phi[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_pho2Phi->cd();
  c_pho2Phi->Write();
  c_pho2Phi->Write();
  c_pho2Phi->Update();
  AddCMS(c_pho2Phi);
  c_pho2Phi->SaveAs("c_pho2Phi_com.pdf");
  c_pho2Phi->SaveAs("c_pho2Phi_com.png");
  c_pho2Phi->SaveAs("c_pho2Phi_com.C");

  //pho2Phi numbers 
  TCanvas* c_G2Phi = new TCanvas( "c_G2Phi", "c_G2Phi", 800, 700 );
  c_G2Phi->SetHighLightColor(2);
  c_G2Phi->SetFillColor(0);
  c_G2Phi->SetBorderMode(0);
  c_G2Phi->SetBorderSize(2);
  c_G2Phi->SetLeftMargin( leftMargin );
  c_G2Phi->SetRightMargin( 1.6*rightMargin );
  c_G2Phi->SetTopMargin( topMargin );
  c_G2Phi->SetBottomMargin( bottomMargin );
  c_G2Phi->SetFrameBorderMode(0);
  c_G2Phi->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_G2Phi_str = "hist_G2Phi_"+Convert.str();
          hist_G2Phi[i] = new TH1F(hist_G2Phi_str.c_str(),"",30,-3,3);
          if(i<1)
          { 
                  std::string draw_G2Phi_str = "pho2Phi>>"+hist_G2Phi_str;
                  tree[i]->Draw(draw_G2Phi_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_G2Phi_str = "phiGamma1>>"+hist_G2Phi_str;
                  tree[i]->Draw(draw_G2Phi_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_G2Phi[i]->SetLineColor(i+2);
          hist_G2Phi[i]->SetLineWidth(3);
          double norm = hist_G2Phi[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_G2Phi[0]->SetTitle(";Phi_{#gamma 2} [GeV];# of events");
  hist_G2Phi[0]->GetYaxis()->SetRangeUser(0,10000);
  hist_G2Phi[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_G2Phi[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_G2Phi->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_G2Phi[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_G2Phi[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_G2Phi->cd();
  c_G2Phi->Write();
  c_G2Phi->Write();
  c_G2Phi->Update();
  AddCMS(c_G2Phi);
  c_G2Phi->SaveAs("c_G2Phi_com.pdf");
  c_G2Phi->SaveAs("c_G2Phi_com.png");
  c_G2Phi->SaveAs("c_G2Phi_com.C");

  //pho1Eta fraction 
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
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho1Eta_str = "hist_pho1Eta_"+Convert.str();
          hist_pho1Eta[i] = new TH1F(hist_pho1Eta_str.c_str(),"",30,-3,3);
          if(i<1)
          { 
                  std::string draw_pho1Eta_str = "pho1Eta>>"+hist_pho1Eta_str;
                  tree[i]->Draw(draw_pho1Eta_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_pho1Eta_str = "etaGamma0>>"+hist_pho1Eta_str;
                  tree[i]->Draw(draw_pho1Eta_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_pho1Eta[i]->SetLineColor(i+2);
          hist_pho1Eta[i]->SetLineWidth(3);
          double norm = hist_pho1Eta[i]->Integral();
          hist_pho1Eta[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho1Eta[0]->SetTitle(";Eta_{#gamma 1} [GeV];Fraction of events");
  hist_pho1Eta[0]->GetYaxis()->SetRangeUser(0,0.13);
  hist_pho1Eta[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho1Eta[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho1Eta->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_pho1Eta[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_pho1Eta[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_pho1Eta->cd();
  c_pho1Eta->Write();
  c_pho1Eta->Write();
  c_pho1Eta->Update();
  AddCMS(c_pho1Eta);
  c_pho1Eta->SaveAs("c_pho1Eta_com.pdf");
  c_pho1Eta->SaveAs("c_pho1Eta_com.png");
  c_pho1Eta->SaveAs("c_pho1Eta_com.C");

  //pho1Eta numbers 
  TCanvas* c_G1Eta = new TCanvas( "c_G1Eta", "c_G1Eta", 800, 700 );
  c_G1Eta->SetHighLightColor(2);
  c_G1Eta->SetFillColor(0);
  c_G1Eta->SetBorderMode(0);
  c_G1Eta->SetBorderSize(2);
  c_G1Eta->SetLeftMargin( leftMargin );
  c_G1Eta->SetRightMargin( 1.6*rightMargin );
  c_G1Eta->SetTopMargin( topMargin );
  c_G1Eta->SetBottomMargin( bottomMargin );
  c_G1Eta->SetFrameBorderMode(0);
  c_G1Eta->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_G1Eta_str = "hist_G1Eta_"+Convert.str();
          hist_G1Eta[i] = new TH1F(hist_G1Eta_str.c_str(),"",30,-3,3);
          if(i<1)
          { 
                  std::string draw_G1Eta_str = "pho1Eta>>"+hist_G1Eta_str;
                  tree[i]->Draw(draw_G1Eta_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_G1Eta_str = "etaGamma0>>"+hist_G1Eta_str;
                  tree[i]->Draw(draw_G1Eta_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_G1Eta[i]->SetLineColor(i+2);
          hist_G1Eta[i]->SetLineWidth(3);
          double norm = hist_G1Eta[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_G1Eta[0]->SetTitle(";Eta_{#gamma 1} [GeV];# of events");
  hist_G1Eta[0]->GetYaxis()->SetRangeUser(0,22000);
  hist_G1Eta[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_G1Eta[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_G1Eta->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_G1Eta[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_G1Eta[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_G1Eta->cd();
  c_G1Eta->Write();
  c_G1Eta->Write();
  c_G1Eta->Update();
  AddCMS(c_G1Eta);
  c_G1Eta->SaveAs("c_G1Eta_com.pdf");
  c_G1Eta->SaveAs("c_G1Eta_com.png");
  c_G1Eta->SaveAs("c_G1Eta_com.C");

  //pho2Eta fraction 
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
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho2Eta_str = "hist_pho2Eta_"+Convert.str();
          hist_pho2Eta[i] = new TH1F(hist_pho2Eta_str.c_str(),"",30,-3,3);
          if(i<1)
          { 
                  std::string draw_pho2Eta_str = "pho2Eta>>"+hist_pho2Eta_str;
                  tree[i]->Draw(draw_pho2Eta_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_pho2Eta_str = "etaGamma1>>"+hist_pho2Eta_str;
                  tree[i]->Draw(draw_pho2Eta_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_pho2Eta[i]->SetLineColor(i+2);
          hist_pho2Eta[i]->SetLineWidth(3);
          double norm = hist_pho2Eta[i]->Integral();
          hist_pho2Eta[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho2Eta[0]->SetTitle(";Eta_{#gamma 2} [GeV];Fraction of events");
  hist_pho2Eta[0]->GetYaxis()->SetRangeUser(0,0.13);
  hist_pho2Eta[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho2Eta[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho2Eta->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_pho2Eta[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_pho2Eta[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_pho2Eta->cd();
  c_pho2Eta->Write();
  c_pho2Eta->Write();
  c_pho2Eta->Update();
  AddCMS(c_pho2Eta);
  c_pho2Eta->SaveAs("c_pho2Eta_com.pdf");
  c_pho2Eta->SaveAs("c_pho2Eta_com.png");
  c_pho2Eta->SaveAs("c_pho2Eta_com.C");

  //pho2Eta numbers 
  TCanvas* c_G2Eta = new TCanvas( "c_G2Eta", "c_G2Eta", 800, 700 );
  c_G2Eta->SetHighLightColor(2);
  c_G2Eta->SetFillColor(0);
  c_G2Eta->SetBorderMode(0);
  c_G2Eta->SetBorderSize(2);
  c_G2Eta->SetLeftMargin( leftMargin );
  c_G2Eta->SetRightMargin( 1.6*rightMargin );
  c_G2Eta->SetTopMargin( topMargin );
  c_G2Eta->SetBottomMargin( bottomMargin );
  c_G2Eta->SetFrameBorderMode(0);
  c_G2Eta->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_G2Eta_str = "hist_G2Eta_"+Convert.str();
          hist_G2Eta[i] = new TH1F(hist_G2Eta_str.c_str(),"",30,-3,3);
          if(i<1)
          { 
                  std::string draw_G2Eta_str = "pho2Eta>>"+hist_G2Eta_str;
                  tree[i]->Draw(draw_G2Eta_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho2Pt/mGammaGamma>1./3. || pho1Pt/mGammaGamma>1./3.) && pho2Pt/mGammaGamma>1./4. && pho1Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_G2Eta_str = "etaGamma1>>"+hist_G2Eta_str;
                  tree[i]->Draw(draw_G2Eta_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_G2Eta[i]->SetLineColor(i+2);
          hist_G2Eta[i]->SetLineWidth(3);
          double norm = hist_G2Eta[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_G2Eta[0]->SetTitle(";Eta_{#gamma 2} [GeV];# of events");
  hist_G2Eta[0]->GetYaxis()->SetRangeUser(0,22000);
  hist_G2Eta[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_G2Eta[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_G2Eta->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_G2Eta[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_G2Eta[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_G2Eta->cd();
  c_G2Eta->Write();
  c_G2Eta->Write();
  c_G2Eta->Update();
  AddCMS(c_G2Eta);
  c_G2Eta->SaveAs("c_G2Eta_com.pdf");
  c_G2Eta->SaveAs("c_G2Eta_com.png");
  c_G2Eta->SaveAs("c_G2Eta_com.C");

  //pho1R9 fraction 
  TCanvas* c_pho1R9 = new TCanvas( "c_pho1R9", "c_pho1R9", 800, 700 );
  c_pho1R9->SetHighLightColor(2);
  c_pho1R9->SetFillColor(0);
  c_pho1R9->SetBorderMode(0);
  c_pho1R9->SetBorderSize(2);
  c_pho1R9->SetLeftMargin( leftMargin );
  c_pho1R9->SetRightMargin( 1.6*rightMargin );
  c_pho1R9->SetTopMargin( topMargin );
  c_pho1R9->SetBottomMargin( bottomMargin );
  c_pho1R9->SetFrameBorderMode(0);
  c_pho1R9->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho1R9_str = "hist_pho1R9_"+Convert.str();
          hist_pho1R9[i] = new TH1F(hist_pho1R9_str.c_str(),"",30,0,1.5);
          if(i<1)
          { 
                  std::string draw_pho1R9_str = "pho1R9>>"+hist_pho1R9_str;
                  tree[i]->Draw(draw_pho1R9_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_pho1R9_str = "r9Gamma0>>"+hist_pho1R9_str;
                  tree[i]->Draw(draw_pho1R9_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_pho1R9[i]->SetLineColor(i+2);
          hist_pho1R9[i]->SetLineWidth(3);
          double norm = hist_pho1R9[i]->Integral();
          hist_pho1R9[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho1R9[0]->SetTitle(";R9_{#gamma 1} [GeV];Fraction of events");
  hist_pho1R9[0]->GetYaxis()->SetRangeUser(0,0.8);
  hist_pho1R9[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho1R9[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho1R9->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_pho1R9[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_pho1R9[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_pho1R9->cd();
  c_pho1R9->Write();
  c_pho1R9->Write();
  c_pho1R9->Update();
  AddCMS(c_pho1R9);
  c_pho1R9->SaveAs("c_pho1R9_com.pdf");
  c_pho1R9->SaveAs("c_pho1R9_com.png");
  c_pho1R9->SaveAs("c_pho1R9_com.C");

  //pho1R9 numbers 
  TCanvas* c_G1R9 = new TCanvas( "c_G1R9", "c_G1R9", 800, 700 );
  c_G1R9->SetHighLightColor(2);
  c_G1R9->SetFillColor(0);
  c_G1R9->SetBorderMode(0);
  c_G1R9->SetBorderSize(2);
  c_G1R9->SetLeftMargin( leftMargin );
  c_G1R9->SetRightMargin( 1.6*rightMargin );
  c_G1R9->SetTopMargin( topMargin );
  c_G1R9->SetBottomMargin( bottomMargin );
  c_G1R9->SetFrameBorderMode(0);
  c_G1R9->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_G1R9_str = "hist_G1R9_"+Convert.str();
          hist_G1R9[i] = new TH1F(hist_G1R9_str.c_str(),"",30,0,1.5);
          if(i<1)
          { 
                  std::string draw_G1R9_str = "pho1R9>>"+hist_G1R9_str;
                  tree[i]->Draw(draw_G1R9_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_G1R9_str = "r9Gamma0>>"+hist_G1R9_str;
                  tree[i]->Draw(draw_G1R9_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_G1R9[i]->SetLineColor(i+2);
          hist_G1R9[i]->SetLineWidth(3);
          double norm = hist_G1R9[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_G1R9[0]->SetTitle(";R9_{#gamma 1} [GeV];# of events");
  hist_G1R9[0]->GetYaxis()->SetRangeUser(0,200000);
  hist_G1R9[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_G1R9[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_G1R9->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_G1R9[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_G1R9[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_G1R9->cd();
  c_G1R9->Write();
  c_G1R9->Write();
  c_G1R9->Update();
  AddCMS(c_G1R9);
  c_G1R9->SaveAs("c_G1R9_com.pdf");
  c_G1R9->SaveAs("c_G1R9_com.png");
  c_G1R9->SaveAs("c_G1R9_com.C");

  //pho2R9 fraction 
  TCanvas* c_pho2R9 = new TCanvas( "c_pho2R9", "c_pho2R9", 800, 700 );
  c_pho2R9->SetHighLightColor(2);
  c_pho2R9->SetFillColor(0);
  c_pho2R9->SetBorderMode(0);
  c_pho2R9->SetBorderSize(2);
  c_pho2R9->SetLeftMargin( leftMargin );
  c_pho2R9->SetRightMargin( 1.6*rightMargin );
  c_pho2R9->SetTopMargin( topMargin );
  c_pho2R9->SetBottomMargin( bottomMargin );
  c_pho2R9->SetFrameBorderMode(0);
  c_pho2R9->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho2R9_str = "hist_pho2R9_"+Convert.str();
          hist_pho2R9[i] = new TH1F(hist_pho2R9_str.c_str(),"",30,0,1.5);
          if(i<1)
          { 
                  std::string draw_pho2R9_str = "pho2R9>>"+hist_pho2R9_str;
                  tree[i]->Draw(draw_pho2R9_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_pho2R9_str = "r9Gamma1>>"+hist_pho2R9_str;
                  tree[i]->Draw(draw_pho2R9_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_pho2R9[i]->SetLineColor(i+2);
          hist_pho2R9[i]->SetLineWidth(3);
          double norm = hist_pho2R9[i]->Integral();
          hist_pho2R9[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho2R9[0]->SetTitle(";R9_{#gamma 2} [GeV];Fraction of events");
  hist_pho2R9[0]->GetYaxis()->SetRangeUser(0,0.8);
  hist_pho2R9[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho2R9[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho2R9->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_pho2R9[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_pho2R9[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_pho2R9->cd();
  c_pho2R9->Write();
  c_pho2R9->Write();
  c_pho2R9->Update();
  AddCMS(c_pho2R9);
  c_pho2R9->SaveAs("c_pho2R9_com.pdf");
  c_pho2R9->SaveAs("c_pho2R9_com.png");
  c_pho2R9->SaveAs("c_pho2R9_com.C");

  //pho2R9 numbers 
  TCanvas* c_G2R9 = new TCanvas( "c_G2R9", "c_G2R9", 800, 700 );
  c_G2R9->SetHighLightColor(2);
  c_G2R9->SetFillColor(0);
  c_G2R9->SetBorderMode(0);
  c_G2R9->SetBorderSize(2);
  c_G2R9->SetLeftMargin( leftMargin );
  c_G2R9->SetRightMargin( 1.6*rightMargin );
  c_G2R9->SetTopMargin( topMargin );
  c_G2R9->SetBottomMargin( bottomMargin );
  c_G2R9->SetFrameBorderMode(0);
  c_G2R9->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_G2R9_str = "hist_G2R9_"+Convert.str();
          hist_G2R9[i] = new TH1F(hist_G2R9_str.c_str(),"",30,0,1.5);
          if(i<1)
          { 
                  std::string draw_G2R9_str = "pho2R9>>"+hist_G2R9_str;
                  tree[i]->Draw(draw_G2R9_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_G2R9_str = "r9Gamma1>>"+hist_G2R9_str;
                  tree[i]->Draw(draw_G2R9_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_G2R9[i]->SetLineColor(i+2);
          hist_G2R9[i]->SetLineWidth(3);
          double norm = hist_G2R9[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_G2R9[0]->SetTitle(";R9_{#gamma 2} [GeV];# of events");
  hist_G2R9[0]->GetYaxis()->SetRangeUser(0,100000);
  hist_G2R9[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_G2R9[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_G2R9->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_G2R9[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_G2R9[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_G2R9->cd();
  c_G2R9->Write();
  c_G2R9->Write();
  c_G2R9->Update();
  AddCMS(c_G2R9);
  c_G2R9->SaveAs("c_G2R9_com.pdf");
  c_G2R9->SaveAs("c_G2R9_com.png");
  c_G2R9->SaveAs("c_G2R9_com.C");

  //pho1SigmaIetaIeta fraction 
  TCanvas* c_pho1SigmaIetaIeta = new TCanvas( "c_pho1SigmaIetaIeta", "c_pho1SigmaIetaIeta", 800, 700 );
  c_pho1SigmaIetaIeta->SetHighLightColor(2);
  c_pho1SigmaIetaIeta->SetFillColor(0);
  c_pho1SigmaIetaIeta->SetBorderMode(0);
  c_pho1SigmaIetaIeta->SetBorderSize(2);
  c_pho1SigmaIetaIeta->SetLeftMargin( leftMargin );
  c_pho1SigmaIetaIeta->SetRightMargin( 1.6*rightMargin );
  c_pho1SigmaIetaIeta->SetTopMargin( topMargin );
  c_pho1SigmaIetaIeta->SetBottomMargin( bottomMargin );
  c_pho1SigmaIetaIeta->SetFrameBorderMode(0);
  c_pho1SigmaIetaIeta->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho1SigmaIetaIeta_str = "hist_pho1SigmaIetaIeta_"+Convert.str();
          hist_pho1SigmaIetaIeta[i] = new TH1F(hist_pho1SigmaIetaIeta_str.c_str(),"",30,0,0.015);
          if(i<1)
          { 
                  std::string draw_pho1SigmaIetaIeta_str = "pho1SigmaIetaIeta>>"+hist_pho1SigmaIetaIeta_str;
                  tree[i]->Draw(draw_pho1SigmaIetaIeta_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_pho1SigmaIetaIeta_str = "sigmaIetaIetaGamma0>>"+hist_pho1SigmaIetaIeta_str;
                  tree[i]->Draw(draw_pho1SigmaIetaIeta_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_pho1SigmaIetaIeta[i]->SetLineColor(i+2);
          hist_pho1SigmaIetaIeta[i]->SetLineWidth(3);
          double norm = hist_pho1SigmaIetaIeta[i]->Integral();
          hist_pho1SigmaIetaIeta[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho1SigmaIetaIeta[0]->SetTitle(";SigmaIetaIeta_{#gamma 1} [GeV];Fraction of events");
  hist_pho1SigmaIetaIeta[0]->GetYaxis()->SetRangeUser(0,0.8);
  hist_pho1SigmaIetaIeta[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho1SigmaIetaIeta[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho1SigmaIetaIeta->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_pho1SigmaIetaIeta[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_pho1SigmaIetaIeta[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_pho1SigmaIetaIeta->cd();
  c_pho1SigmaIetaIeta->Write();
  c_pho1SigmaIetaIeta->Write();
  c_pho1SigmaIetaIeta->Update();
  AddCMS(c_pho1SigmaIetaIeta);
  c_pho1SigmaIetaIeta->SaveAs("c_pho1SigmaIetaIeta_com.pdf");
  c_pho1SigmaIetaIeta->SaveAs("c_pho1SigmaIetaIeta_com.png");
  c_pho1SigmaIetaIeta->SaveAs("c_pho1SigmaIetaIeta_com.C");

  //pho1SigmaIetaIeta numbers 
  TCanvas* c_G1SigmaIetaIeta = new TCanvas( "c_G1SigmaIetaIeta", "c_G1SigmaIetaIeta", 800, 700 );
  c_G1SigmaIetaIeta->SetHighLightColor(2);
  c_G1SigmaIetaIeta->SetFillColor(0);
  c_G1SigmaIetaIeta->SetBorderMode(0);
  c_G1SigmaIetaIeta->SetBorderSize(2);
  c_G1SigmaIetaIeta->SetLeftMargin( leftMargin );
  c_G1SigmaIetaIeta->SetRightMargin( 1.6*rightMargin );
  c_G1SigmaIetaIeta->SetTopMargin( topMargin );
  c_G1SigmaIetaIeta->SetBottomMargin( bottomMargin );
  c_G1SigmaIetaIeta->SetFrameBorderMode(0);
  c_G1SigmaIetaIeta->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_G1SigmaIetaIeta_str = "hist_G1SigmaIetaIeta_"+Convert.str();
          hist_G1SigmaIetaIeta[i] = new TH1F(hist_G1SigmaIetaIeta_str.c_str(),"",30,0,0.015);
          if(i<1)
          { 
                  std::string draw_G1SigmaIetaIeta_str = "pho1SigmaIetaIeta>>"+hist_G1SigmaIetaIeta_str;
                  tree[i]->Draw(draw_G1SigmaIetaIeta_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_G1SigmaIetaIeta_str = "sigmaIetaIetaGamma0>>"+hist_G1SigmaIetaIeta_str;
                  tree[i]->Draw(draw_G1SigmaIetaIeta_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_G1SigmaIetaIeta[i]->SetLineColor(i+2);
          hist_G1SigmaIetaIeta[i]->SetLineWidth(3);
          double norm = hist_G1SigmaIetaIeta[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_G1SigmaIetaIeta[0]->SetTitle(";SigmaIetaIeta_{#gamma 1} [GeV];# of events");
  hist_G1SigmaIetaIeta[0]->GetYaxis()->SetRangeUser(0,150000);
  hist_G1SigmaIetaIeta[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_G1SigmaIetaIeta[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_G1SigmaIetaIeta->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_G1SigmaIetaIeta[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_G1SigmaIetaIeta[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_G1SigmaIetaIeta->cd();
  c_G1SigmaIetaIeta->Write();
  c_G1SigmaIetaIeta->Write();
  c_G1SigmaIetaIeta->Update();
  AddCMS(c_G1SigmaIetaIeta);
  c_G1SigmaIetaIeta->SaveAs("c_G1SigmaIetaIeta_com.pdf");
  c_G1SigmaIetaIeta->SaveAs("c_G1SigmaIetaIeta_com.png");
  c_G1SigmaIetaIeta->SaveAs("c_G1SigmaIetaIeta_com.C");

  //pho2SigmaIetaIeta fraction 
  TCanvas* c_pho2SigmaIetaIeta = new TCanvas( "c_pho2SigmaIetaIeta", "c_pho2SigmaIetaIeta", 800, 700 );
  c_pho2SigmaIetaIeta->SetHighLightColor(2);
  c_pho2SigmaIetaIeta->SetFillColor(0);
  c_pho2SigmaIetaIeta->SetBorderMode(0);
  c_pho2SigmaIetaIeta->SetBorderSize(2);
  c_pho2SigmaIetaIeta->SetLeftMargin( leftMargin );
  c_pho2SigmaIetaIeta->SetRightMargin( 1.6*rightMargin );
  c_pho2SigmaIetaIeta->SetTopMargin( topMargin );
  c_pho2SigmaIetaIeta->SetBottomMargin( bottomMargin );
  c_pho2SigmaIetaIeta->SetFrameBorderMode(0);
  c_pho2SigmaIetaIeta->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_pho2SigmaIetaIeta_str = "hist_pho2SigmaIetaIeta_"+Convert.str();
          hist_pho2SigmaIetaIeta[i] = new TH1F(hist_pho2SigmaIetaIeta_str.c_str(),"",30,0,0.015);
          if(i<1)
          { 
                  std::string draw_pho2SigmaIetaIeta_str = "pho2SigmaIetaIeta>>"+hist_pho2SigmaIetaIeta_str;
                  tree[i]->Draw(draw_pho2SigmaIetaIeta_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_pho2SigmaIetaIeta_str = "sigmaIetaIetaGamma1>>"+hist_pho2SigmaIetaIeta_str;
                  tree[i]->Draw(draw_pho2SigmaIetaIeta_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_pho2SigmaIetaIeta[i]->SetLineColor(i+2);
          hist_pho2SigmaIetaIeta[i]->SetLineWidth(3);
          double norm = hist_pho2SigmaIetaIeta[i]->Integral();
          hist_pho2SigmaIetaIeta[i]->Scale(1/norm);
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_pho2SigmaIetaIeta[0]->SetTitle(";SigmaIetaIeta_{#gamma 2} [GeV];Fraction of events");
  hist_pho2SigmaIetaIeta[0]->GetYaxis()->SetRangeUser(0,0.8);
  hist_pho2SigmaIetaIeta[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_pho2SigmaIetaIeta[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_pho2SigmaIetaIeta->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_pho2SigmaIetaIeta[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_pho2SigmaIetaIeta[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_pho2SigmaIetaIeta->cd();
  c_pho2SigmaIetaIeta->Write();
  c_pho2SigmaIetaIeta->Write();
  c_pho2SigmaIetaIeta->Update();
  AddCMS(c_pho2SigmaIetaIeta);
  c_pho2SigmaIetaIeta->SaveAs("c_pho2SigmaIetaIeta_com.pdf");
  c_pho2SigmaIetaIeta->SaveAs("c_pho2SigmaIetaIeta_com.png");
  c_pho2SigmaIetaIeta->SaveAs("c_pho2SigmaIetaIeta_com.C");

  //pho2SigmaIetaIeta numbers 
  TCanvas* c_G2SigmaIetaIeta = new TCanvas( "c_G2SigmaIetaIeta", "c_G2SigmaIetaIeta", 800, 700 );
  c_G2SigmaIetaIeta->SetHighLightColor(2);
  c_G2SigmaIetaIeta->SetFillColor(0);
  c_G2SigmaIetaIeta->SetBorderMode(0);
  c_G2SigmaIetaIeta->SetBorderSize(2);
  c_G2SigmaIetaIeta->SetLeftMargin( leftMargin );
  c_G2SigmaIetaIeta->SetRightMargin( 1.6*rightMargin );
  c_G2SigmaIetaIeta->SetTopMargin( topMargin );
  c_G2SigmaIetaIeta->SetBottomMargin( bottomMargin );
  c_G2SigmaIetaIeta->SetFrameBorderMode(0);
  c_G2SigmaIetaIeta->SetFrameBorderMode(0);
  pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->Draw();
  pad1->cd();
  for(int i = 0; i < n; i++){
          ostringstream Convert;
          Convert<<i;
          std::string hist_G2SigmaIetaIeta_str = "hist_G2SigmaIetaIeta_"+Convert.str();
          hist_G2SigmaIetaIeta[i] = new TH1F(hist_G2SigmaIetaIeta_str.c_str(),"",30,0,0.015);
          if(i<1)
          { 
                  std::string draw_G2SigmaIetaIeta_str = "pho2SigmaIetaIeta>>"+hist_G2SigmaIetaIeta_str;
                  tree[i]->Draw(draw_G2SigmaIetaIeta_str.c_str(),"( mGammaGamma>103. && mGammaGamma<160. && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 ) && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && pho1R9>0.5 && pho2R9>0.5 && HLTDecision[82]","HISTSAME");
          }
          else 
          {
                  std::string draw_G2SigmaIetaIeta_str = "sigmaIetaIetaGamma1>>"+hist_G2SigmaIetaIeta_str;
                  tree[i]->Draw(draw_G2SigmaIetaIeta_str.c_str(),"((( ptGamma0/h_mass) > 1./3.) && ((ptGamma1/h_mass )> 1./4.) && fabs(etaGamma0)<=1.4442  && fabs(etaGamma1)<=1.4442 && r9Gamma0 >= 0.5 && r9Gamma1>=0.5 && h_mass>=103 && h_mass<=160 )","*ESAME");
          } 
          hist_G2SigmaIetaIeta[i]->SetLineColor(i+2);
          hist_G2SigmaIetaIeta[i]->SetLineWidth(3);
          double norm = hist_G2SigmaIetaIeta[i]->Integral();
          std::cout << "i = " << i << "\n" << std::endl;
  }
  hist_G2SigmaIetaIeta[0]->SetTitle(";SigmaIetaIeta_{#gamma 2} [GeV];# of events");
  hist_G2SigmaIetaIeta[0]->GetYaxis()->SetRangeUser(0,150000);
  hist_G2SigmaIetaIeta[0]->GetYaxis()->SetTitleOffset(1.4);
  hist_G2SigmaIetaIeta[0]->GetXaxis()->SetTitleOffset(1.4);
   tex = new TLatex(0.27,0.78,"data 2016 Razor");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
   tex = new TLatex(0.27,0.68,"data 2016 Mt2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetTextColor(3);
   tex->SetLineWidth(2);
   tex->Draw("SAME");
  c_G2SigmaIetaIeta->cd();
  pad2 = new TPad("pad2","pad2",0,0.05,1,0.27);
  pad2->Draw();
  pad2->cd();
  hist_ratio = (TH1F*)hist_G2SigmaIetaIeta[0]->Clone("ratio");
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetMinimum(0.5);
  hist_ratio->SetMaximum(1.5);
  hist_ratio->Sumw2();
  hist_ratio->SetStats(0);
  hist_ratio->Divide(hist_G2SigmaIetaIeta[1]);
  hist_ratio->SetMarkerStyle(21);
  hist_ratio->Draw("ep");
  hist_ratio->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("Razor/Mt2 ");
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
  c_G2SigmaIetaIeta->cd();
  c_G2SigmaIetaIeta->Write();
  c_G2SigmaIetaIeta->Write();
  c_G2SigmaIetaIeta->Update();
  AddCMS(c_G2SigmaIetaIeta);
  c_G2SigmaIetaIeta->SaveAs("c_G2SigmaIetaIeta_com.pdf");
  c_G2SigmaIetaIeta->SaveAs("c_G2SigmaIetaIeta_com.png");
  c_G2SigmaIetaIeta->SaveAs("c_G2SigmaIetaIeta_com.C");

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
