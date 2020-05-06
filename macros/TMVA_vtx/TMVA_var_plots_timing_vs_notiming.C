
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

void TMVA_var_plots_timing_vs_notiming(string pu_ = "PU0")
{

   TFile *file_NoTiming = new TFile(("HggRazorUpgradeTiming_"+pu_+"_NoTiming_vtx.root").c_str(),"READ");
   TFile *file_Timing = new TFile(("HggRazorUpgradeTiming_"+pu_+"_Timing_vtx.root").c_str(),"READ");
   TTree *tree_sig_NoTiming = (TTree*)file_NoTiming->Get("TreeS");
   TTree *tree_bkg_NoTiming = (TTree*)file_NoTiming->Get("TreeB");

   TTree *tree_sig_Timing = (TTree*)file_Timing->Get("TreeS");
   TTree *tree_bkg_Timing = (TTree*)file_Timing->Get("TreeB");

   double MaxY=0;

   TCanvas * myC = new TCanvas("c1","c1",100,100,800,700);
   
   gStyle->SetOptStat(0);
   

   TH1F * h_ptasym_sig_Timing = new TH1F("h_ptasym_sig_Timing","h_ptasym_sig_Timing",100,-1.1,1.1); 
   TH1F * h_ptasym_sig_NoTiming = new TH1F("h_ptasym_sig_NoTiming","h_ptasym_sig_NoTiming",100,-1.1,1.1); 
   tree_sig_Timing->Draw("ptasym>>h_ptasym_sig_Timing");
   tree_sig_NoTiming->Draw("ptasym>>h_ptasym_sig_NoTiming");
   MaxY=0;
   if(h_ptasym_sig_Timing->GetMaximum()>MaxY) MaxY = h_ptasym_sig_Timing->GetMaximum();
   if(h_ptasym_sig_NoTiming->GetMaximum()>MaxY) MaxY = h_ptasym_sig_NoTiming->GetMaximum();
   h_ptasym_sig_NoTiming->SetLineColor(2); 
   h_ptasym_sig_NoTiming->SetTitle("");
   h_ptasym_sig_NoTiming->GetXaxis()->SetTitle("ptasym"); 
   h_ptasym_sig_NoTiming->GetYaxis()->SetTitle("Events"); 
   h_ptasym_sig_NoTiming->GetYaxis()->SetTitleOffset(1.3);
   h_ptasym_sig_NoTiming->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_ptasym_sig_Timing->SetLineColor(4); 
   h_ptasym_sig_NoTiming->Draw();
   h_ptasym_sig_Timing->Draw("same");
   TLegend *leg_ptasym_sig = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_ptasym_sig->SetBorderSize(0);
        leg_ptasym_sig->SetTextSize(0.04);
        leg_ptasym_sig->SetLineColor(1);
        leg_ptasym_sig->SetLineStyle(1);
        leg_ptasym_sig->SetLineWidth(1);
        leg_ptasym_sig->SetFillColor(0);
        leg_ptasym_sig->SetFillStyle(1001);
        leg_ptasym_sig->AddEntry(h_ptasym_sig_NoTiming, "without timing" ,"l");
        leg_ptasym_sig->AddEntry(h_ptasym_sig_Timing, "with timing" ,"l");
        leg_ptasym_sig->Draw();
   myC->SetLogy(1);
   myC->SaveAs(("plots/ptasym_"+pu_+"_sig.pdf").c_str());	
   myC->SaveAs(("plots/ptasym_"+pu_+"_sig.png").c_str());	

   TH1F * h_ptasym_bkg_Timing = new TH1F("h_ptasym_bkg_Timing","h_ptasym_bkg_Timing",100,-1.1,1.1); 
   TH1F * h_ptasym_bkg_NoTiming = new TH1F("h_ptasym_bkg_NoTiming","h_ptasym_bkg_NoTiming",100,-1.1,1.1); 
   tree_bkg_Timing->Draw("ptasym>>h_ptasym_bkg_Timing");
   tree_bkg_NoTiming->Draw("ptasym>>h_ptasym_bkg_NoTiming");
   MaxY=0;
   if(h_ptasym_bkg_Timing->GetMaximum()>MaxY) MaxY = h_ptasym_bkg_Timing->GetMaximum();
   if(h_ptasym_bkg_NoTiming->GetMaximum()>MaxY) MaxY = h_ptasym_bkg_NoTiming->GetMaximum();
   h_ptasym_bkg_NoTiming->SetLineColor(2); 
   h_ptasym_bkg_NoTiming->SetTitle("");
   h_ptasym_bkg_NoTiming->GetXaxis()->SetTitle("ptasym"); 
   h_ptasym_bkg_NoTiming->GetYaxis()->SetTitle("Events"); 
   h_ptasym_bkg_NoTiming->GetYaxis()->SetTitleOffset(1.3);
   h_ptasym_bkg_NoTiming->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_ptasym_bkg_Timing->SetLineColor(4); 
   h_ptasym_bkg_NoTiming->Draw();
   h_ptasym_bkg_Timing->Draw("same");
   TLegend *leg_ptasym_bkg = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_ptasym_bkg->SetBorderSize(0);
        leg_ptasym_bkg->SetTextSize(0.04);
        leg_ptasym_bkg->SetLineColor(1);
        leg_ptasym_bkg->SetLineStyle(1);
        leg_ptasym_bkg->SetLineWidth(1);
        leg_ptasym_bkg->SetFillColor(0);
        leg_ptasym_bkg->SetFillStyle(1001);
        leg_ptasym_bkg->AddEntry(h_ptasym_bkg_NoTiming, "without timing" ,"l");
        leg_ptasym_bkg->AddEntry(h_ptasym_bkg_Timing, "with timing" ,"l");
        leg_ptasym_bkg->Draw();
   myC->SetLogy(1);
   myC->SaveAs(("plots/ptasym_"+pu_+"_bkg.pdf").c_str());	
   myC->SaveAs(("plots/ptasym_"+pu_+"_bkg.png").c_str());	



   TH1F * h_ptbal_sig_Timing = new TH1F("h_ptbal_sig_Timing","h_ptbal_sig_Timing",200,-50,150); 
   TH1F * h_ptbal_sig_NoTiming = new TH1F("h_ptbal_sig_NoTiming","h_ptbal_sig_NoTiming",200,-50,150); 
   tree_sig_Timing->Draw("ptbal>>h_ptbal_sig_Timing");
   tree_sig_NoTiming->Draw("ptbal>>h_ptbal_sig_NoTiming");
   MaxY=0;
   if(h_ptbal_sig_Timing->GetMaximum()>MaxY) MaxY = h_ptbal_sig_Timing->GetMaximum();
   if(h_ptbal_sig_NoTiming->GetMaximum()>MaxY) MaxY = h_ptbal_sig_NoTiming->GetMaximum();
   h_ptbal_sig_NoTiming->SetLineColor(2); 
   h_ptbal_sig_NoTiming->SetTitle("");
   h_ptbal_sig_NoTiming->GetXaxis()->SetTitle("ptbal"); 
   h_ptbal_sig_NoTiming->GetYaxis()->SetTitle("Events"); 
   h_ptbal_sig_NoTiming->GetYaxis()->SetTitleOffset(1.3);
   h_ptbal_sig_NoTiming->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_ptbal_sig_Timing->SetLineColor(4); 
   h_ptbal_sig_NoTiming->Draw();
   h_ptbal_sig_Timing->Draw("same");
   TLegend *leg_ptbal_sig = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_ptbal_sig->SetBorderSize(0);
        leg_ptbal_sig->SetTextSize(0.04);
        leg_ptbal_sig->SetLineColor(1);
        leg_ptbal_sig->SetLineStyle(1);
        leg_ptbal_sig->SetLineWidth(1);
        leg_ptbal_sig->SetFillColor(0);
        leg_ptbal_sig->SetFillStyle(1001);
        leg_ptbal_sig->AddEntry(h_ptbal_sig_NoTiming, "without timing" ,"l");
        leg_ptbal_sig->AddEntry(h_ptbal_sig_Timing, "with timing" ,"l");
        leg_ptbal_sig->Draw();
   myC->SetLogy(1);
   myC->SaveAs(("plots/ptbal_"+pu_+"_sig.pdf").c_str());	
   myC->SaveAs(("plots/ptbal_"+pu_+"_sig.png").c_str());	

   TH1F * h_ptbal_bkg_Timing = new TH1F("h_ptbal_bkg_Timing","h_ptbal_bkg_Timing",200,-50,150); 
   TH1F * h_ptbal_bkg_NoTiming = new TH1F("h_ptbal_bkg_NoTiming","h_ptbal_bkg_NoTiming",200,-50,150); 
   tree_bkg_Timing->Draw("ptbal>>h_ptbal_bkg_Timing");
   tree_bkg_NoTiming->Draw("ptbal>>h_ptbal_bkg_NoTiming");
   MaxY=0;
   if(h_ptbal_bkg_Timing->GetMaximum()>MaxY) MaxY = h_ptbal_bkg_Timing->GetMaximum();
   if(h_ptbal_bkg_NoTiming->GetMaximum()>MaxY) MaxY = h_ptbal_bkg_NoTiming->GetMaximum();
   h_ptbal_bkg_NoTiming->SetLineColor(2); 
   h_ptbal_bkg_NoTiming->SetTitle("");
   h_ptbal_bkg_NoTiming->GetXaxis()->SetTitle("ptbal"); 
   h_ptbal_bkg_NoTiming->GetYaxis()->SetTitle("Events"); 
   h_ptbal_bkg_NoTiming->GetYaxis()->SetTitleOffset(1.3);
   h_ptbal_bkg_NoTiming->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_ptbal_bkg_Timing->SetLineColor(4); 
   h_ptbal_bkg_NoTiming->Draw();
   h_ptbal_bkg_Timing->Draw("same");
   TLegend *leg_ptbal_bkg = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_ptbal_bkg->SetBorderSize(0);
        leg_ptbal_bkg->SetTextSize(0.04);
        leg_ptbal_bkg->SetLineColor(1);
        leg_ptbal_bkg->SetLineStyle(1);
        leg_ptbal_bkg->SetLineWidth(1);
        leg_ptbal_bkg->SetFillColor(0);
        leg_ptbal_bkg->SetFillStyle(1001);
        leg_ptbal_bkg->AddEntry(h_ptbal_bkg_NoTiming, "without timing" ,"l");
        leg_ptbal_bkg->AddEntry(h_ptbal_bkg_Timing, "with timing" ,"l");
        leg_ptbal_bkg->Draw();
   myC->SetLogy(1);
   myC->SaveAs(("plots/ptbal_"+pu_+"_bkg.pdf").c_str());	
   myC->SaveAs(("plots/ptbal_"+pu_+"_bkg.png").c_str());	


   TH1F * h_logsumpt2_sig_Timing = new TH1F("h_logsumpt2_sig_Timing","h_logsumpt2_sig_Timing",100,-3,10); 
   TH1F * h_logsumpt2_sig_NoTiming = new TH1F("h_logsumpt2_sig_NoTiming","h_logsumpt2_sig_NoTiming",100,-3,10); 
   tree_sig_Timing->Draw("logsumpt2>>h_logsumpt2_sig_Timing");
   tree_sig_NoTiming->Draw("logsumpt2>>h_logsumpt2_sig_NoTiming");
   MaxY=0;
   if(h_logsumpt2_sig_Timing->GetMaximum()>MaxY) MaxY = h_logsumpt2_sig_Timing->GetMaximum();
   if(h_logsumpt2_sig_NoTiming->GetMaximum()>MaxY) MaxY = h_logsumpt2_sig_NoTiming->GetMaximum();
   h_logsumpt2_sig_NoTiming->SetLineColor(2); 
   h_logsumpt2_sig_NoTiming->SetTitle("");
   h_logsumpt2_sig_NoTiming->GetXaxis()->SetTitle("logsumpt2"); 
   h_logsumpt2_sig_NoTiming->GetYaxis()->SetTitle("Events"); 
   h_logsumpt2_sig_NoTiming->GetYaxis()->SetTitleOffset(1.3);
   h_logsumpt2_sig_NoTiming->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_logsumpt2_sig_Timing->SetLineColor(4); 
   h_logsumpt2_sig_NoTiming->Draw();
   h_logsumpt2_sig_Timing->Draw("same");
   TLegend *leg_logsumpt2_sig = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_logsumpt2_sig->SetBorderSize(0);
        leg_logsumpt2_sig->SetTextSize(0.04);
        leg_logsumpt2_sig->SetLineColor(1);
        leg_logsumpt2_sig->SetLineStyle(1);
        leg_logsumpt2_sig->SetLineWidth(1);
        leg_logsumpt2_sig->SetFillColor(0);
        leg_logsumpt2_sig->SetFillStyle(1001);
        leg_logsumpt2_sig->AddEntry(h_logsumpt2_sig_NoTiming, "without timing" ,"l");
        leg_logsumpt2_sig->AddEntry(h_logsumpt2_sig_Timing, "with timing" ,"l");
        leg_logsumpt2_sig->Draw();
   myC->SetLogy(1);
   myC->SaveAs(("plots/logsumpt2_"+pu_+"_sig.pdf").c_str());	
   myC->SaveAs(("plots/logsumpt2_"+pu_+"_sig.png").c_str());	

   TH1F * h_logsumpt2_bkg_Timing = new TH1F("h_logsumpt2_bkg_Timing","h_logsumpt2_bkg_Timing",100,-3,10); 
   TH1F * h_logsumpt2_bkg_NoTiming = new TH1F("h_logsumpt2_bkg_NoTiming","h_logsumpt2_bkg_NoTiming",100,-3,10); 
   tree_bkg_Timing->Draw("logsumpt2>>h_logsumpt2_bkg_Timing");
   tree_bkg_NoTiming->Draw("logsumpt2>>h_logsumpt2_bkg_NoTiming");
   MaxY=0;
   if(h_logsumpt2_bkg_Timing->GetMaximum()>MaxY) MaxY = h_logsumpt2_bkg_Timing->GetMaximum();
   if(h_logsumpt2_bkg_NoTiming->GetMaximum()>MaxY) MaxY = h_logsumpt2_bkg_NoTiming->GetMaximum();
   h_logsumpt2_bkg_NoTiming->SetLineColor(2); 
   h_logsumpt2_bkg_NoTiming->SetTitle("");
   h_logsumpt2_bkg_NoTiming->GetXaxis()->SetTitle("logsumpt2"); 
   h_logsumpt2_bkg_NoTiming->GetYaxis()->SetTitle("Events"); 
   h_logsumpt2_bkg_NoTiming->GetYaxis()->SetTitleOffset(1.3);
   h_logsumpt2_bkg_NoTiming->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_logsumpt2_bkg_Timing->SetLineColor(4); 
   h_logsumpt2_bkg_NoTiming->Draw();
   h_logsumpt2_bkg_Timing->Draw("same");
   TLegend *leg_logsumpt2_bkg = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_logsumpt2_bkg->SetBorderSize(0);
        leg_logsumpt2_bkg->SetTextSize(0.04);
        leg_logsumpt2_bkg->SetLineColor(1);
        leg_logsumpt2_bkg->SetLineStyle(1);
        leg_logsumpt2_bkg->SetLineWidth(1);
        leg_logsumpt2_bkg->SetFillColor(0);
        leg_logsumpt2_bkg->SetFillStyle(1001);
        leg_logsumpt2_bkg->AddEntry(h_logsumpt2_bkg_NoTiming, "without timing" ,"l");
        leg_logsumpt2_bkg->AddEntry(h_logsumpt2_bkg_Timing, "with timing" ,"l");
        leg_logsumpt2_bkg->Draw();
   myC->SetLogy(1);
   myC->SaveAs(("plots/logsumpt2_"+pu_+"_bkg.pdf").c_str());	
   myC->SaveAs(("plots/logsumpt2_"+pu_+"_bkg.png").c_str());	

}
