
#include "TPaveText.h"
#include <string>
#include <iostream>
#include <math.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>

void TMVA_var_plots_PU0_vs_PU140(string pu_ = "PU0")
{

  TFile *file_PU0 =  new TFile("HggRazorUpgradeTiming_PU0_Timing_vtx.root");
  TFile *file_PU140 =  new TFile("HggRazorUpgradeTiming_PU140_Timing_vtx.root");

  TTree *tree_PU0 = (TTree*)file_PU0->Get("TreeS");
  TTree *tree_PU140 = (TTree*)file_PU140->Get("TreeS");

  int N_Entries_PU0 = tree_PU0->GetEntries();
  int N_Entries_PU140 = tree_PU140->GetEntries();

   double MaxY=0;
   TCanvas * myC = new TCanvas("c1","c1",100,100,800,700);
   gStyle->SetOptStat(0);

   TH1F * h_ptasym_PU0 = new TH1F("h_ptasym_PU0","h_ptasym_PU0",100,-1.1,1.1); 
   TH1F * h_ptasym_PU140 = new TH1F("h_ptasym_PU140","h_ptasym_PU140",100,-1.1,1.1); 
   tree_PU0->Draw("ptasym>>h_ptasym_PU0");
   tree_PU140->Draw("ptasym>>h_ptasym_PU140");
   h_ptasym_PU0->Scale(1.0/N_Entries_PU0);  
   h_ptasym_PU140->Scale(1.0/N_Entries_PU140);  
   MaxY=0;
   h_ptasym_PU140->SetLineColor(2); 
   h_ptasym_PU140->SetTitle("");
   h_ptasym_PU140->GetXaxis()->SetTitle("ptasym"); 
   h_ptasym_PU140->GetYaxis()->SetTitle("Events"); 
   h_ptasym_PU140->GetYaxis()->SetTitleOffset(1.3);
   h_ptasym_PU140->GetYaxis()->SetRangeUser(1.0e-4,1.0);
   h_ptasym_PU0->SetLineColor(4); 
   h_ptasym_PU140->Draw();
   h_ptasym_PU0->Draw("same");
   TLegend *leg_ptasym_sig = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_ptasym_sig->SetBorderSize(0);
        leg_ptasym_sig->SetTextSize(0.04);
        leg_ptasym_sig->SetLineColor(1);
        leg_ptasym_sig->SetLineStyle(1);
        leg_ptasym_sig->SetLineWidth(1);
        leg_ptasym_sig->SetFillColor(0);
        leg_ptasym_sig->SetFillStyle(1001);
        leg_ptasym_sig->AddEntry(h_ptasym_PU140, "PU140" ,"l");
        leg_ptasym_sig->AddEntry(h_ptasym_PU0, "PU0" ,"l");
        leg_ptasym_sig->Draw();
   myC->SetLogy(1);
   myC->SaveAs("~/www/sharebox/tomyself/tmp/ptasym_PU0_vs_PU140.pdf");	
   myC->SaveAs("~/www/sharebox/tomyself/tmp/ptasym_PU0_vs_PU140.png");	


   TH1F * h_ptbal_PU0 = new TH1F("h_ptbal_PU0","h_ptbal_PU0",200,-50,150); 
   TH1F * h_ptbal_PU140 = new TH1F("h_ptbal_PU140","h_ptbal_PU140",200,-50,150); 
   tree_PU0->Draw("ptbal>>h_ptbal_PU0");
   tree_PU140->Draw("ptbal>>h_ptbal_PU140");
   h_ptbal_PU0->Scale(1.0/N_Entries_PU0);  
   h_ptbal_PU140->Scale(1.0/N_Entries_PU140);  
   MaxY=0;
   h_ptbal_PU140->SetLineColor(2); 
   h_ptbal_PU140->SetTitle("");
   h_ptbal_PU140->GetXaxis()->SetTitle("ptbal"); 
   h_ptbal_PU140->GetYaxis()->SetTitle("Events"); 
   h_ptbal_PU140->GetYaxis()->SetTitleOffset(1.3);
   h_ptbal_PU140->GetYaxis()->SetRangeUser(1.0e-4,1.0);
   h_ptbal_PU0->SetLineColor(4); 
   h_ptbal_PU140->Draw();
   h_ptbal_PU0->Draw("same");
   TLegend *leg_ptbal_sig = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_ptbal_sig->SetBorderSize(0);
        leg_ptbal_sig->SetTextSize(0.04);
        leg_ptbal_sig->SetLineColor(1);
        leg_ptbal_sig->SetLineStyle(1);
        leg_ptbal_sig->SetLineWidth(1);
        leg_ptbal_sig->SetFillColor(0);
        leg_ptbal_sig->SetFillStyle(1001);
        leg_ptbal_sig->AddEntry(h_ptbal_PU140, "PU140" ,"l");
        leg_ptbal_sig->AddEntry(h_ptbal_PU0, "PU0" ,"l");
        leg_ptbal_sig->Draw();
   myC->SetLogy(1);
   myC->SaveAs("~/www/sharebox/tomyself/tmp/ptbal_PU0_vs_PU140.pdf");	
   myC->SaveAs("~/www/sharebox/tomyself/tmp/ptbal_PU0_vs_PU140.png");	


   TH1F * h_logsumpt2_PU0 = new TH1F("h_logsumpt2_PU0","h_logsumpt2_PU0",100,-3,10); 
   TH1F * h_logsumpt2_PU140 = new TH1F("h_logsumpt2_PU140","h_logsumpt2_PU140",100,-3,10); 
   tree_PU0->Draw("logsumpt2>>h_logsumpt2_PU0");
   tree_PU140->Draw("logsumpt2>>h_logsumpt2_PU140");
   h_logsumpt2_PU0->Scale(1.0/N_Entries_PU0);  
   h_logsumpt2_PU140->Scale(1.0/N_Entries_PU140);  
   MaxY=0;
   h_logsumpt2_PU140->SetLineColor(2); 
   h_logsumpt2_PU140->SetTitle("");
   h_logsumpt2_PU140->GetXaxis()->SetTitle("logsumpt2"); 
   h_logsumpt2_PU140->GetYaxis()->SetTitle("Events"); 
   h_logsumpt2_PU140->GetYaxis()->SetTitleOffset(1.3);
   h_logsumpt2_PU140->GetYaxis()->SetRangeUser(1.0e-4,1.0);
   h_logsumpt2_PU0->SetLineColor(4); 
   h_logsumpt2_PU140->Draw();
   h_logsumpt2_PU0->Draw("same");
   TLegend *leg_logsumpt2_sig = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_logsumpt2_sig->SetBorderSize(0);
        leg_logsumpt2_sig->SetTextSize(0.04);
        leg_logsumpt2_sig->SetLineColor(1);
        leg_logsumpt2_sig->SetLineStyle(1);
        leg_logsumpt2_sig->SetLineWidth(1);
        leg_logsumpt2_sig->SetFillColor(0);
        leg_logsumpt2_sig->SetFillStyle(1001);
        leg_logsumpt2_sig->AddEntry(h_logsumpt2_PU140, "PU140" ,"l");
        leg_logsumpt2_sig->AddEntry(h_logsumpt2_PU0, "PU0" ,"l");
        leg_logsumpt2_sig->Draw();
   myC->SetLogy(1);
   myC->SaveAs("~/www/sharebox/tomyself/tmp/logsumpt2_PU0_vs_PU140.pdf");	
   myC->SaveAs("~/www/sharebox/tomyself/tmp/logsumpt2_PU0_vs_PU140.png");	

   TH1F * h_chi2_pho_vtx_PU0 = new TH1F("h_chi2_pho_vtx_PU0","h_chi2_pho_vtx_PU0",100,0,50); 
   TH1F * h_chi2_pho_vtx_PU140 = new TH1F("h_chi2_pho_vtx_PU140","h_chi2_pho_vtx_PU140",100,0,50); 
   tree_PU0->Draw("chi2_pho_vtx>>h_chi2_pho_vtx_PU0");
   tree_PU140->Draw("chi2_pho_vtx>>h_chi2_pho_vtx_PU140");
   h_chi2_pho_vtx_PU0->Scale(1.0/N_Entries_PU0);  
   h_chi2_pho_vtx_PU140->Scale(1.0/N_Entries_PU140);  
   MaxY=0;
   h_chi2_pho_vtx_PU140->SetLineColor(2); 
   h_chi2_pho_vtx_PU140->SetTitle("");
   h_chi2_pho_vtx_PU140->GetXaxis()->SetTitle("chi2_pho_vtx"); 
   h_chi2_pho_vtx_PU140->GetYaxis()->SetTitle("Events"); 
   h_chi2_pho_vtx_PU140->GetYaxis()->SetTitleOffset(1.3);
   h_chi2_pho_vtx_PU140->GetYaxis()->SetRangeUser(1.0e-4,1.0);
   h_chi2_pho_vtx_PU0->SetLineColor(4); 
   h_chi2_pho_vtx_PU140->Draw();
   h_chi2_pho_vtx_PU0->Draw("same");
   TLegend *leg_chi2_pho_vtx_sig = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_chi2_pho_vtx_sig->SetBorderSize(0);
        leg_chi2_pho_vtx_sig->SetTextSize(0.04);
        leg_chi2_pho_vtx_sig->SetLineColor(1);
        leg_chi2_pho_vtx_sig->SetLineStyle(1);
        leg_chi2_pho_vtx_sig->SetLineWidth(1);
        leg_chi2_pho_vtx_sig->SetFillColor(0);
        leg_chi2_pho_vtx_sig->SetFillStyle(1001);
        leg_chi2_pho_vtx_sig->AddEntry(h_chi2_pho_vtx_PU140, "PU140" ,"l");
        leg_chi2_pho_vtx_sig->AddEntry(h_chi2_pho_vtx_PU0, "PU0" ,"l");
        leg_chi2_pho_vtx_sig->Draw();
   myC->SetLogy(1);
   myC->SaveAs("~/www/sharebox/tomyself/tmp/chi2_pho_vtx_PU0_vs_PU140.pdf");	
   myC->SaveAs("~/www/sharebox/tomyself/tmp/chi2_pho_vtx_PU0_vs_PU140.png");	

   TH1F * h_logchi2_pho_vtx_PU0 = new TH1F("h_logchi2_pho_vtx_PU0","h_logchi2_pho_vtx_PU0",100,-10,10); 
   TH1F * h_logchi2_pho_vtx_PU140 = new TH1F("h_logchi2_pho_vtx_PU140","h_logchi2_pho_vtx_PU140",100,-10,10); 
   tree_PU0->Draw("logchi2_pho_vtx>>h_logchi2_pho_vtx_PU0");
   tree_PU140->Draw("logchi2_pho_vtx>>h_logchi2_pho_vtx_PU140");
   h_logchi2_pho_vtx_PU0->Scale(1.0/N_Entries_PU0);  
   h_logchi2_pho_vtx_PU140->Scale(1.0/N_Entries_PU140);  
   MaxY=0;
   h_logchi2_pho_vtx_PU140->SetLineColor(2); 
   h_logchi2_pho_vtx_PU140->SetTitle("");
   h_logchi2_pho_vtx_PU140->GetXaxis()->SetTitle("logchi2_pho_vtx"); 
   h_logchi2_pho_vtx_PU140->GetYaxis()->SetTitle("Events"); 
   h_logchi2_pho_vtx_PU140->GetYaxis()->SetTitleOffset(1.3);
   h_logchi2_pho_vtx_PU140->GetYaxis()->SetRangeUser(1.0e-4,1.0);
   h_logchi2_pho_vtx_PU0->SetLineColor(4); 
   h_logchi2_pho_vtx_PU140->Draw();
   h_logchi2_pho_vtx_PU0->Draw("same");
   TLegend *leg_logchi2_pho_vtx_sig = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_logchi2_pho_vtx_sig->SetBorderSize(0);
        leg_logchi2_pho_vtx_sig->SetTextSize(0.04);
        leg_logchi2_pho_vtx_sig->SetLineColor(1);
        leg_logchi2_pho_vtx_sig->SetLineStyle(1);
        leg_logchi2_pho_vtx_sig->SetLineWidth(1);
        leg_logchi2_pho_vtx_sig->SetFillColor(0);
        leg_logchi2_pho_vtx_sig->SetFillStyle(1001);
        leg_logchi2_pho_vtx_sig->AddEntry(h_logchi2_pho_vtx_PU140, "PU140" ,"l");
        leg_logchi2_pho_vtx_sig->AddEntry(h_logchi2_pho_vtx_PU0, "PU0" ,"l");
        leg_logchi2_pho_vtx_sig->Draw();
   myC->SetLogy(1);
   myC->SaveAs("~/www/sharebox/tomyself/tmp/logchi2_pho_vtx_PU0_vs_PU140.pdf");	
   myC->SaveAs("~/www/sharebox/tomyself/tmp/logchi2_pho_vtx_PU0_vs_PU140.png");	



}
