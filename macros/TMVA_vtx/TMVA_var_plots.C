
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

void TMVA_var_plots(string inputFilename="HggRazorUpgradeTiming_PU0_NoTiming_vtx.root", string plot_suffix="PU0_NoTiming")
{

   TFile *file_in = new TFile(inputFilename.c_str(),"READ");
   TTree *treeS = (TTree*)file_in->Get("TreeS");
   TTree *treeB = (TTree*)file_in->Get("TreeB");

   double MaxY=0;
   

   TCanvas * myC = new TCanvas("c1","c1",100,100,800,700);


   gStyle->SetOptStat(0);
   

   TH1F * h_ptasym_sig = new TH1F("h_ptasym_sig","h_ptasym_sig",100,-1.1,1.1); 
   TH1F * h_ptasym_bkg = new TH1F("h_ptasym_bkg","h_ptasym_bkg",100,-1.1,1.1); 
   treeS->Draw("ptasym>>h_ptasym_sig");
   treeB->Draw("ptasym>>h_ptasym_bkg");
   MaxY=0;
   if(h_ptasym_sig->GetMaximum()>MaxY) MaxY = h_ptasym_sig->GetMaximum();
   if(h_ptasym_bkg->GetMaximum()>MaxY) MaxY = h_ptasym_bkg->GetMaximum();
   h_ptasym_bkg->SetLineColor(2); 
   h_ptasym_bkg->SetTitle(("ptasym_"+plot_suffix).c_str()); 
   h_ptasym_bkg->GetXaxis()->SetTitle("ptasym"); 
   h_ptasym_bkg->GetYaxis()->SetTitle("Events"); 
   h_ptasym_bkg->GetYaxis()->SetTitleOffset(1.3);
   h_ptasym_bkg->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_ptasym_sig->SetLineColor(4); 
   h_ptasym_bkg->Draw();
   h_ptasym_sig->Draw("same");
        TLegend *leg_ptasym = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_ptasym->SetBorderSize(0);
        leg_ptasym->SetTextSize(0.04);
        leg_ptasym->SetLineColor(1);
        leg_ptasym->SetLineStyle(1);
        leg_ptasym->SetLineWidth(1);
        leg_ptasym->SetFillColor(0);
        leg_ptasym->SetFillStyle(1001);
        leg_ptasym->AddEntry(h_ptasym_bkg, "non-gen-matched vtx" ,"l");
        leg_ptasym->AddEntry(h_ptasym_sig, "gen-matched vtx" ,"l");
        leg_ptasym->Draw();

   myC->SetLogy(1);
   myC->SaveAs(("plots/ptasym_"+plot_suffix+".pdf").c_str());	
   myC->SaveAs(("plots/ptasym_"+plot_suffix+".png").c_str());	

   TH1F * h_ptbal_sig = new TH1F("h_ptbal_sig","h_ptbal_sig",200,-50,150); 
   TH1F * h_ptbal_bkg = new TH1F("h_ptbal_bkg","h_ptbal_bkg",200,-50,150); 
   treeS->Draw("ptbal>>h_ptbal_sig");
   treeB->Draw("ptbal>>h_ptbal_bkg");
   MaxY=0;
   if(h_ptbal_sig->GetMaximum()>MaxY) MaxY = h_ptbal_sig->GetMaximum();
   if(h_ptbal_bkg->GetMaximum()>MaxY) MaxY = h_ptbal_bkg->GetMaximum();
   h_ptbal_bkg->SetLineColor(2); 
   h_ptbal_bkg->SetTitle(("ptbal_"+plot_suffix).c_str()); 
   h_ptbal_bkg->GetXaxis()->SetTitle("ptbal"); 
   h_ptbal_bkg->GetYaxis()->SetTitle("Events"); 
   h_ptbal_bkg->GetYaxis()->SetTitleOffset(1.3);
   h_ptbal_bkg->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_ptbal_sig->SetLineColor(4); 
   h_ptbal_bkg->Draw();
   h_ptbal_sig->Draw("same");
        TLegend *leg_ptbal = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_ptbal->SetBorderSize(0);
        leg_ptbal->SetTextSize(0.04);
        leg_ptbal->SetLineColor(1);
        leg_ptbal->SetLineStyle(1);
        leg_ptbal->SetLineWidth(1);
        leg_ptbal->SetFillColor(0);
        leg_ptbal->SetFillStyle(1001);
        leg_ptbal->AddEntry(h_ptbal_bkg, "non-gen-matched vtx" ,"l");
        leg_ptbal->AddEntry(h_ptbal_sig, "gen-matched vtx" ,"l");
        leg_ptbal->Draw();
   myC->SaveAs(("plots/ptbal_"+plot_suffix+".pdf").c_str());	
   myC->SaveAs(("plots/ptbal_"+plot_suffix+".png").c_str());	

   TH1F * h_logsumpt2_sig = new TH1F("h_logsumpt2_sig","h_logsumpt2_sig",100,-3,10); 
   TH1F * h_logsumpt2_bkg = new TH1F("h_logsumpt2_bkg","h_logsumpt2_bkg",100,-3,10); 
   treeS->Draw("logsumpt2>>h_logsumpt2_sig");
   treeB->Draw("logsumpt2>>h_logsumpt2_bkg");
   MaxY=0;
   if(h_logsumpt2_sig->GetMaximum()>MaxY) MaxY = h_logsumpt2_sig->GetMaximum();
   if(h_logsumpt2_bkg->GetMaximum()>MaxY) MaxY = h_logsumpt2_bkg->GetMaximum();
   h_logsumpt2_bkg->SetLineColor(2); 
   h_logsumpt2_bkg->SetTitle(("logsumpt2_"+plot_suffix).c_str()); 
   h_logsumpt2_bkg->GetXaxis()->SetTitle("logsumpt2"); 
   h_logsumpt2_bkg->GetYaxis()->SetTitle("Events"); 
   h_logsumpt2_bkg->GetYaxis()->SetTitleOffset(1.3);
   h_logsumpt2_bkg->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_logsumpt2_sig->SetLineColor(4); 
   h_logsumpt2_bkg->Draw();
   h_logsumpt2_sig->Draw("same");
        TLegend *leg_logsumpt2 = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_logsumpt2->SetBorderSize(0);
        leg_logsumpt2->SetTextSize(0.04);
        leg_logsumpt2->SetLineColor(1);
        leg_logsumpt2->SetLineStyle(1);
        leg_logsumpt2->SetLineWidth(1);
        leg_logsumpt2->SetFillColor(0);
        leg_logsumpt2->SetFillStyle(1001);
        leg_logsumpt2->AddEntry(h_logsumpt2_bkg, "non-gen-matched vtx" ,"l");
        leg_logsumpt2->AddEntry(h_logsumpt2_sig, "gen-matched vtx" ,"l");
        leg_logsumpt2->Draw();
   myC->SaveAs(("plots/logsumpt2_"+plot_suffix+".pdf").c_str());	
   myC->SaveAs(("plots/logsumpt2_"+plot_suffix+".png").c_str());	

   TH1F * h_logAbsSumPt_sig = new TH1F("h_logAbsSumPt_sig","h_logAbsSumPt_sig",100,-3,10); 
   TH1F * h_logAbsSumPt_bkg = new TH1F("h_logAbsSumPt_bkg","h_logAbsSumPt_bkg",100,-3,10); 
   treeS->Draw("log(vtxPt)>>h_logAbsSumPt_sig");
   treeB->Draw("log(vtxPt)>>h_logAbsSumPt_bkg");
   MaxY=0;
   if(h_logAbsSumPt_sig->GetMaximum()>MaxY) MaxY = h_logAbsSumPt_sig->GetMaximum();
   if(h_logAbsSumPt_bkg->GetMaximum()>MaxY) MaxY = h_logAbsSumPt_bkg->GetMaximum();
   h_logAbsSumPt_bkg->SetLineColor(2); 
   h_logAbsSumPt_bkg->SetTitle(("logAbsSumPt_"+plot_suffix).c_str()); 
   h_logAbsSumPt_bkg->GetXaxis()->SetTitle("logAbsSumPt"); 
   h_logAbsSumPt_bkg->GetYaxis()->SetTitle("Events"); 
   h_logAbsSumPt_bkg->GetYaxis()->SetTitleOffset(1.3);
   h_logAbsSumPt_bkg->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_logAbsSumPt_sig->SetLineColor(4); 
   h_logAbsSumPt_bkg->Draw();
   h_logAbsSumPt_sig->Draw("same");
        TLegend *leg_logAbsSumPt = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_logAbsSumPt->SetBorderSize(0);
        leg_logAbsSumPt->SetTextSize(0.04);
        leg_logAbsSumPt->SetLineColor(1);
        leg_logAbsSumPt->SetLineStyle(1);
        leg_logAbsSumPt->SetLineWidth(1);
        leg_logAbsSumPt->SetFillColor(0);
        leg_logAbsSumPt->SetFillStyle(1001);
        leg_logAbsSumPt->AddEntry(h_logAbsSumPt_bkg, "non-gen-matched vtx" ,"l");
        leg_logAbsSumPt->AddEntry(h_logAbsSumPt_sig, "gen-matched vtx" ,"l");
        leg_logAbsSumPt->Draw();
   myC->SaveAs(("plots/logAbsSumPt_"+plot_suffix+".pdf").c_str());	
   myC->SaveAs(("plots/logAbsSumPt_"+plot_suffix+".png").c_str());	


   TH1F * h_logsumAbsPt_sig = new TH1F("h_logsumAbsPt_sig","h_logsumAbsPt_sig",100,-3,10); 
   TH1F * h_logsumAbsPt_bkg = new TH1F("h_logsumAbsPt_bkg","h_logsumAbsPt_bkg",100,-3,10); 
   treeS->Draw("log(vtxSumPt)>>h_logsumAbsPt_sig");
   treeB->Draw("log(vtxSumPt)>>h_logsumAbsPt_bkg");
   MaxY=0;
   if(h_logsumAbsPt_sig->GetMaximum()>MaxY) MaxY = h_logsumAbsPt_sig->GetMaximum();
   if(h_logsumAbsPt_bkg->GetMaximum()>MaxY) MaxY = h_logsumAbsPt_bkg->GetMaximum();
   h_logsumAbsPt_bkg->SetLineColor(2); 
   h_logsumAbsPt_bkg->SetTitle(("logsumAbsPt_"+plot_suffix).c_str()); 
   h_logsumAbsPt_bkg->GetXaxis()->SetTitle("logsumAbsPt"); 
   h_logsumAbsPt_bkg->GetYaxis()->SetTitle("Events"); 
   h_logsumAbsPt_bkg->GetYaxis()->SetTitleOffset(1.3);
   h_logsumAbsPt_bkg->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_logsumAbsPt_sig->SetLineColor(4); 
   h_logsumAbsPt_bkg->Draw();
   h_logsumAbsPt_sig->Draw("same");
        TLegend *leg_logsumAbsPt = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_logsumAbsPt->SetBorderSize(0);
        leg_logsumAbsPt->SetTextSize(0.04);
        leg_logsumAbsPt->SetLineColor(1);
        leg_logsumAbsPt->SetLineStyle(1);
        leg_logsumAbsPt->SetLineWidth(1);
        leg_logsumAbsPt->SetFillColor(0);
        leg_logsumAbsPt->SetFillStyle(1001);
        leg_logsumAbsPt->AddEntry(h_logsumAbsPt_bkg, "non-gen-matched vtx" ,"l");
        leg_logsumAbsPt->AddEntry(h_logsumAbsPt_sig, "gen-matched vtx" ,"l");
        leg_logsumAbsPt->Draw();
   myC->SaveAs(("plots/logsumAbsPt_"+plot_suffix+".pdf").c_str());	
   myC->SaveAs(("plots/logsumAbsPt_"+plot_suffix+".png").c_str());	



   TH1F * h_limPullToConv_sig = new TH1F("h_limPullToConv_sig","h_limPullToConv_sig",100,-1,11); 
   TH1F * h_limPullToConv_bkg = new TH1F("h_limPullToConv_bkg","h_limPullToConv_bkg",100,-1,11); 
   treeS->Draw("limPullToConv>>h_limPullToConv_sig");
   treeB->Draw("limPullToConv>>h_limPullToConv_bkg");
   MaxY=0;
   if(h_limPullToConv_sig->GetMaximum()>MaxY) MaxY = h_limPullToConv_sig->GetMaximum();
   if(h_limPullToConv_bkg->GetMaximum()>MaxY) MaxY = h_limPullToConv_bkg->GetMaximum();
   h_limPullToConv_bkg->SetLineColor(2); 
   h_limPullToConv_bkg->SetTitle(("limPullToConv_"+plot_suffix).c_str()); 
   h_limPullToConv_bkg->GetXaxis()->SetTitle("limPullToConv"); 
   h_limPullToConv_bkg->GetYaxis()->SetTitle("Events"); 
   h_limPullToConv_bkg->GetYaxis()->SetTitleOffset(1.3);
   h_limPullToConv_bkg->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_limPullToConv_sig->SetLineColor(4); 
   h_limPullToConv_bkg->Draw();
   h_limPullToConv_sig->Draw("same");
        TLegend *leg_limPullToConv = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_limPullToConv->SetBorderSize(0);
        leg_limPullToConv->SetTextSize(0.04);
        leg_limPullToConv->SetLineColor(1);
        leg_limPullToConv->SetLineStyle(1);
        leg_limPullToConv->SetLineWidth(1);
        leg_limPullToConv->SetFillColor(0);
        leg_limPullToConv->SetFillStyle(1001);
        leg_limPullToConv->AddEntry(h_limPullToConv_bkg, "non-gen-matched vtx" ,"l");
        leg_limPullToConv->AddEntry(h_limPullToConv_sig, "gen-matched vtx" ,"l");
        leg_limPullToConv->Draw();
   myC->SetLogy(1);
   myC->SaveAs(("plots/limPullToConv_"+plot_suffix+".pdf").c_str());	
   myC->SaveAs(("plots/limPullToConv_"+plot_suffix+".png").c_str());	

   TH1F * h_nConv_sig = new TH1F("h_nConv_sig","h_nConv_sig",100,-1,3); 
   TH1F * h_nConv_bkg = new TH1F("h_nConv_bkg","h_nConv_bkg",100,-1,3); 
   treeS->Draw("nConv>>h_nConv_sig");
   treeB->Draw("nConv>>h_nConv_bkg");
   MaxY=0;
   if(h_nConv_sig->GetMaximum()>MaxY) MaxY = h_nConv_sig->GetMaximum();
   if(h_nConv_bkg->GetMaximum()>MaxY) MaxY = h_nConv_bkg->GetMaximum();
   h_nConv_bkg->SetLineColor(2); 
   h_nConv_bkg->SetTitle(("nConv_"+plot_suffix).c_str()); 
   h_nConv_bkg->GetXaxis()->SetTitle("nConv"); 
   h_nConv_bkg->GetYaxis()->SetTitle("Events"); 
   h_nConv_bkg->GetYaxis()->SetTitleOffset(1.3);
   h_nConv_bkg->GetYaxis()->SetRangeUser(0.1,12*MaxY);
   h_nConv_sig->SetLineColor(4); 
   h_nConv_bkg->Draw();
   h_nConv_sig->Draw("same");
        TLegend *leg_nConv = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_nConv->SetBorderSize(0);
        leg_nConv->SetTextSize(0.04);
        leg_nConv->SetLineColor(1);
        leg_nConv->SetLineStyle(1);
        leg_nConv->SetLineWidth(1);
        leg_nConv->SetFillColor(0);
        leg_nConv->SetFillStyle(1001);
        leg_nConv->AddEntry(h_nConv_bkg, "non-gen-matched vtx" ,"l");
        leg_nConv->AddEntry(h_nConv_sig, "gen-matched vtx" ,"l");
        leg_nConv->Draw();
   myC->SetLogy(1);
   myC->SaveAs(("plots/nConv_"+plot_suffix+".pdf").c_str());	
   myC->SaveAs(("plots/nConv_"+plot_suffix+".png").c_str());	


   TH2F * h_vtx_gg_Pt_sig = new TH2F("h_vtx_gg_Pt_sig","h_vtx_gg_Pt_sig",100,0,100,30,0,30); 
   TH2F * h_vtx_gg_Pt_bkg = new TH2F("h_vtx_gg_Pt_bkg","h_vtx_gg_Pt_bkg",100,0,100,30,0,30); 
   treeS->Draw("vtxPt:diphoPt>>h_vtx_gg_Pt_sig");
   treeB->Draw("vtxPt:diphoPt>>h_vtx_gg_Pt_bkg");
   h_vtx_gg_Pt_bkg->SetTitle(""); 
   h_vtx_gg_Pt_bkg->GetXaxis()->SetTitle("Pt-diphoton"); 
   h_vtx_gg_Pt_bkg->GetYaxis()->SetTitle("Pt-vertex (associated tracks)"); 
   h_vtx_gg_Pt_bkg->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/vtx_gg_Pt_"+plot_suffix+"_bkg.pdf").c_str());	
   myC->SaveAs(("plots/vtx_gg_Pt_"+plot_suffix+"_bkg.png").c_str());	
   h_vtx_gg_Pt_sig->SetTitle(""); 
   h_vtx_gg_Pt_sig->GetXaxis()->SetTitle("Pt-diphoton"); 
   h_vtx_gg_Pt_sig->GetYaxis()->SetTitle("Pt-vertex (associated tracks)"); 
   h_vtx_gg_Pt_sig->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/vtx_gg_Pt_"+plot_suffix+"_sig.pdf").c_str());	
   myC->SaveAs(("plots/vtx_gg_Pt_"+plot_suffix+"_sig.png").c_str());	

   TH2F * h_dz_logsumpt2_Pt_sig = new TH2F("h_dz_logsumpt2_Pt_sig","h_dz_logsumpt2_Pt_sig",100,-0.06,0.06,100,-3,10); 
   TH2F * h_dz_logsumpt2_Pt_bkg = new TH2F("h_dz_logsumpt2_Pt_bkg","h_dz_logsumpt2_Pt_bkg",100,-20,20,100,-3,10); 
   treeS->Draw("logsumpt2:vtxdZ>>h_dz_logsumpt2_Pt_sig");
   treeB->Draw("logsumpt2:vtxdZ>>h_dz_logsumpt2_Pt_bkg");
   h_dz_logsumpt2_Pt_bkg->SetTitle(""); 
   h_dz_logsumpt2_Pt_bkg->GetXaxis()->SetTitle("vertexZ - genVertexZ"); 
   h_dz_logsumpt2_Pt_bkg->GetYaxis()->SetTitle("logsumpt2"); 
   h_dz_logsumpt2_Pt_bkg->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/dz_logsumpt2_Pt_"+plot_suffix+"_bkg.pdf").c_str());	
   myC->SaveAs(("plots/dz_logsumpt2_Pt_"+plot_suffix+"_bkg.png").c_str());	
   h_dz_logsumpt2_Pt_sig->SetTitle(""); 
   h_dz_logsumpt2_Pt_sig->GetXaxis()->SetTitle("vertexZ - genVertexZ"); 
   h_dz_logsumpt2_Pt_sig->GetYaxis()->SetTitle("logsumpt2"); 
   h_dz_logsumpt2_Pt_sig->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/dz_logsumpt2_Pt_"+plot_suffix+"_sig.pdf").c_str());	
   myC->SaveAs(("plots/dz_logsumpt2_Pt_"+plot_suffix+"_sig.png").c_str());	

   TH2F * h_dz_logAbsSumPt_Pt_sig = new TH2F("h_dz_logAbsSumPt_Pt_sig","h_dz_logAbsSumPt_Pt_sig",100,-0.06,0.06,100,-3,10); 
   TH2F * h_dz_logAbsSumPt_Pt_bkg = new TH2F("h_dz_logAbsSumPt_Pt_bkg","h_dz_logAbsSumPt_Pt_bkg",100,-20,20,100,-3,10); 
   treeS->Draw("log(vtxPt):vtxdZ>>h_dz_logAbsSumPt_Pt_sig");
   treeB->Draw("log(vtxPt):vtxdZ>>h_dz_logAbsSumPt_Pt_bkg");
   h_dz_logAbsSumPt_Pt_bkg->SetTitle(""); 
   h_dz_logAbsSumPt_Pt_bkg->GetXaxis()->SetTitle("vertexZ - genVertexZ"); 
   h_dz_logAbsSumPt_Pt_bkg->GetYaxis()->SetTitle("logAbsSumPt"); 
   h_dz_logAbsSumPt_Pt_bkg->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/dz_logAbsSumPt_Pt_"+plot_suffix+"_bkg.pdf").c_str());	
   myC->SaveAs(("plots/dz_logAbsSumPt_Pt_"+plot_suffix+"_bkg.png").c_str());	
   h_dz_logAbsSumPt_Pt_sig->SetTitle(""); 
   h_dz_logAbsSumPt_Pt_sig->GetXaxis()->SetTitle("vertexZ - genVertexZ"); 
   h_dz_logAbsSumPt_Pt_sig->GetYaxis()->SetTitle("logAbsSumPt"); 
   h_dz_logAbsSumPt_Pt_sig->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/dz_logAbsSumPt_Pt_"+plot_suffix+"_sig.pdf").c_str());	
   myC->SaveAs(("plots/dz_logAbsSumPt_Pt_"+plot_suffix+"_sig.png").c_str());	

   TH2F * h_dz_logsumAbsPt_Pt_sig = new TH2F("h_dz_logsumAbsPt_Pt_sig","h_dz_logsumAbsPt_Pt_sig",100,-0.06,0.06,100,-3,10); 
   TH2F * h_dz_logsumAbsPt_Pt_bkg = new TH2F("h_dz_logsumAbsPt_Pt_bkg","h_dz_logsumAbsPt_Pt_bkg",100,-20,20,100,-3,10); 
   treeS->Draw("log(vtxSumPt):vtxdZ>>h_dz_logsumAbsPt_Pt_sig");
   treeB->Draw("log(vtxSumPt):vtxdZ>>h_dz_logsumAbsPt_Pt_bkg");
   h_dz_logsumAbsPt_Pt_bkg->SetTitle(""); 
   h_dz_logsumAbsPt_Pt_bkg->GetXaxis()->SetTitle("vertexZ - genVertexZ"); 
   h_dz_logsumAbsPt_Pt_bkg->GetYaxis()->SetTitle("logsumAbsPt"); 
   h_dz_logsumAbsPt_Pt_bkg->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/dz_logsumAbsPt_Pt_"+plot_suffix+"_bkg.pdf").c_str());	
   myC->SaveAs(("plots/dz_logsumAbsPt_Pt_"+plot_suffix+"_bkg.png").c_str());	
   h_dz_logsumAbsPt_Pt_sig->SetTitle(""); 
   h_dz_logsumAbsPt_Pt_sig->GetXaxis()->SetTitle("vertexZ - genVertexZ"); 
   h_dz_logsumAbsPt_Pt_sig->GetYaxis()->SetTitle("logsumAbsPt"); 
   h_dz_logsumAbsPt_Pt_sig->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/dz_logsumAbsPt_Pt_"+plot_suffix+"_sig.pdf").c_str());	
   myC->SaveAs(("plots/dz_logsumAbsPt_Pt_"+plot_suffix+"_sig.png").c_str());	


   double y_max_h_dzdtchi2_sig = 3.0e6;
   if(plot_suffix.find("PU0")!=std::string::npos) y_max_h_dzdtchi2_sig = 3000.0;
   TH2F * h_dzdtchi2_logsumpt2_Pt_sig = new TH2F("h_dzdtchi2_logsumpt2_Pt_sig","h_dzdtchi2_logsumpt2_Pt_sig",100,0,20,100,-3,10); 
   TH2F * h_dzdtchi2_logsumpt2_Pt_bkg = new TH2F("h_dzdtchi2_logsumpt2_Pt_bkg","h_dzdtchi2_logsumpt2_Pt_bkg",1000,0, y_max_h_dzdtchi2_sig, 100,-3,10); 
   treeS->Draw("logsumpt2:(vtxdT*vtxdT/(0.06*0.06)+vtxdZ*vtxdZ/(0.01*0.01))>>h_dzdtchi2_logsumpt2_Pt_sig");
   treeB->Draw("logsumpt2:(vtxdT*vtxdT/(0.06*0.06)+vtxdZ*vtxdZ/(0.01*0.01))>>h_dzdtchi2_logsumpt2_Pt_bkg");
   h_dzdtchi2_logsumpt2_Pt_bkg->SetTitle(""); 
   h_dzdtchi2_logsumpt2_Pt_bkg->GetXaxis()->SetTitle("#chi^{2}(dz, dt)"); 
   h_dzdtchi2_logsumpt2_Pt_bkg->GetYaxis()->SetTitle("logsumpt2"); 
   h_dzdtchi2_logsumpt2_Pt_bkg->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/dzdtchi2_logsumpt2_Pt_"+plot_suffix+"_bkg.pdf").c_str());	
   myC->SaveAs(("plots/dzdtchi2_logsumpt2_Pt_"+plot_suffix+"_bkg.png").c_str());	
   h_dzdtchi2_logsumpt2_Pt_sig->SetTitle(""); 
   h_dzdtchi2_logsumpt2_Pt_sig->GetXaxis()->SetTitle("#chi^{2}(dz, dt)"); 
   h_dzdtchi2_logsumpt2_Pt_sig->GetYaxis()->SetTitle("logsumpt2"); 
   h_dzdtchi2_logsumpt2_Pt_sig->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/dzdtchi2_logsumpt2_Pt_"+plot_suffix+"_sig.pdf").c_str());	
   myC->SaveAs(("plots/dzdtchi2_logsumpt2_Pt_"+plot_suffix+"_sig.png").c_str());	


   TH2F * h_pvNtrack_logsumpt2_Pt_sig = new TH2F("h_pvNtrack_logsumpt2_Pt_sig","h_pvNtrack_logsumpt2_Pt_sig",150,0,150,100,-3,10); 
   TH2F * h_pvNtrack_logsumpt2_Pt_bkg = new TH2F("h_pvNtrack_logsumpt2_Pt_bkg","h_pvNtrack_logsumpt2_Pt_bkg",150,0,150,100,-3,10); 
   treeS->Draw("logsumpt2:pvNtrack>>h_pvNtrack_logsumpt2_Pt_sig");
   treeB->Draw("logsumpt2:pvNtrack>>h_pvNtrack_logsumpt2_Pt_bkg");
   h_pvNtrack_logsumpt2_Pt_bkg->SetTitle(""); 
   h_pvNtrack_logsumpt2_Pt_bkg->GetXaxis()->SetTitle("# of associated tracks"); 
   h_pvNtrack_logsumpt2_Pt_bkg->GetYaxis()->SetTitle("logsumpt2"); 
   h_pvNtrack_logsumpt2_Pt_bkg->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/pvNtrack_logsumpt2_Pt_"+plot_suffix+"_bkg.pdf").c_str());	
   myC->SaveAs(("plots/pvNtrack_logsumpt2_Pt_"+plot_suffix+"_bkg.png").c_str());	
   h_pvNtrack_logsumpt2_Pt_sig->SetTitle(""); 
   h_pvNtrack_logsumpt2_Pt_sig->GetXaxis()->SetTitle("# of associated tracks"); 
   h_pvNtrack_logsumpt2_Pt_sig->GetYaxis()->SetTitle("logsumpt2"); 
   h_pvNtrack_logsumpt2_Pt_sig->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/pvNtrack_logsumpt2_Pt_"+plot_suffix+"_sig.pdf").c_str());	
   myC->SaveAs(("plots/pvNtrack_logsumpt2_Pt_"+plot_suffix+"_sig.png").c_str());	



   TH1F * h_pvNtrack_sig = new TH1F("h_pvNtrack_sig","h_pvNtrack_sig",140,0,140); 
   TH1F * h_pvNtrack_bkg = new TH1F("h_pvNtrack_bkg","h_pvNtrack_bkg",140,0,140); 
   treeS->Draw("pvNtrack>>h_pvNtrack_sig");
   treeB->Draw("pvNtrack>>h_pvNtrack_bkg");
   MaxY=0;
   if(h_pvNtrack_sig->GetMaximum()>MaxY) MaxY = h_pvNtrack_sig->GetMaximum();
   if(h_pvNtrack_bkg->GetMaximum()>MaxY) MaxY = h_pvNtrack_bkg->GetMaximum();
   h_pvNtrack_bkg->SetLineColor(2); 
   h_pvNtrack_bkg->SetTitle(("pvNtrack_"+plot_suffix).c_str()); 
   h_pvNtrack_bkg->GetXaxis()->SetTitle("# of associated tracks"); 
   h_pvNtrack_bkg->GetYaxis()->SetTitle("Events"); 
   h_pvNtrack_bkg->GetYaxis()->SetTitleOffset(1.3);
   h_pvNtrack_bkg->GetYaxis()->SetRangeUser(1,12*MaxY);
   h_pvNtrack_sig->SetLineColor(4); 
   h_pvNtrack_bkg->Draw();
   h_pvNtrack_sig->Draw("same");
        TLegend *leg_pvNtrack = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_pvNtrack->SetBorderSize(0);
        leg_pvNtrack->SetTextSize(0.04);
        leg_pvNtrack->SetLineColor(1);
        leg_pvNtrack->SetLineStyle(1);
        leg_pvNtrack->SetLineWidth(1);
        leg_pvNtrack->SetFillColor(0);
        leg_pvNtrack->SetFillStyle(1001);
        leg_pvNtrack->AddEntry(h_pvNtrack_bkg, "non-gen-matched vtx" ,"l");
        leg_pvNtrack->AddEntry(h_pvNtrack_sig, "gen-matched vtx" ,"l");
        leg_pvNtrack->Draw();

   myC->SetLogy(1);
   myC->SaveAs(("plots/pvNtrack_"+plot_suffix+".pdf").c_str());	
   myC->SaveAs(("plots/pvNtrack_"+plot_suffix+".png").c_str());	

}


