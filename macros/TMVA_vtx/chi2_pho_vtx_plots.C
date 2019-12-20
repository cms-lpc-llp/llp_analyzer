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
#include <TPad.h>


void chi2_pho_vtx_plots(string inputFilename="../../HggRazorUpgradeTiming_PU0_Timing.root", string plot_suffix="PU0_Timing")
{

	
	TFile f_in(inputFilename.c_str());
	
	TTree *t_in = (TTree*)f_in.Get("HggRazor");


	double MaxY=0;
   	TCanvas * myC = new TCanvas("c1","c1",100,100,800,700);
   
        myC->cd();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);

	gStyle->SetOptStat(0);
        gStyle->SetPadTickY(1);

	
   TH1F * h_virtualVtxdZ_smalldEta = new TH1F("h_virtualVtxdZ_smalldEta","h_virtualVtxdZ_smalldEta",40,-100,100);
   TH1F * h_virtualVtxdZ_largedEta = new TH1F("h_virtualVtxdZ_largedEta","h_virtualVtxdZ_largedEta",40,-100,100);
   t_in->Draw("(10.0*virtualVtxdZ)>>h_virtualVtxdZ_smalldEta","abs(pho1Eta-pho2Eta)<0.8");
   t_in->Draw("(10.0*virtualVtxdZ)>>h_virtualVtxdZ_largedEta","abs(pho1Eta-pho2Eta)>0.8");
   MaxY=0;
   if(h_virtualVtxdZ_smalldEta->GetMaximum()>MaxY) MaxY = h_virtualVtxdZ_smalldEta->GetMaximum();
   if(h_virtualVtxdZ_largedEta->GetMaximum()>MaxY) MaxY = h_virtualVtxdZ_largedEta->GetMaximum();
   h_virtualVtxdZ_largedEta->SetLineColor(4);
   h_virtualVtxdZ_largedEta->SetTitle("");
   h_virtualVtxdZ_largedEta->GetXaxis()->SetTitle("z - z_{gen} / mm");
   h_virtualVtxdZ_largedEta->GetYaxis()->SetTitle("Events/5mm");
   h_virtualVtxdZ_largedEta->GetYaxis()->SetTitleOffset(1.5);
   h_virtualVtxdZ_largedEta->GetYaxis()->SetRangeUser(0.0,1.2*MaxY);
   h_virtualVtxdZ_smalldEta->SetLineColor(2);
   h_virtualVtxdZ_largedEta->Draw();
   h_virtualVtxdZ_smalldEta->Draw("same");
        TLegend *leg_virtualVtxdZ = new TLegend(0.65, 0.78, 0.85, 0.89, NULL,"brNDC");
        leg_virtualVtxdZ->SetBorderSize(0);
        leg_virtualVtxdZ->SetTextSize(0.04);
        leg_virtualVtxdZ->SetLineColor(1);
        leg_virtualVtxdZ->SetLineStyle(1);
        leg_virtualVtxdZ->SetLineWidth(1);
        leg_virtualVtxdZ->SetFillColor(0);
        leg_virtualVtxdZ->SetFillStyle(1001);
        leg_virtualVtxdZ->AddEntry(h_virtualVtxdZ_largedEta, "|#Delta#eta_{#gamma#gamma}|>0.8" ,"l");
        leg_virtualVtxdZ->AddEntry(h_virtualVtxdZ_smalldEta, "|#Delta#eta_{#gamma#gamma}|<0.8" ,"l");
        leg_virtualVtxdZ->Draw();
   myC->SetLogy(0);
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/virtualVtxdZ_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/virtualVtxdZ_"+plot_suffix+".C").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/virtualVtxdZ_"+plot_suffix+".png").c_str());

   double int_smallEta_middle = h_virtualVtxdZ_smalldEta->Integral(17,21);
   double int_smallEta_all = 1.0*h_virtualVtxdZ_smalldEta->GetEntries();
   cout<<"smalldEta: "<<int_smallEta_middle<<" out of "<<int_smallEta_all<<"  =  "<<100.0*int_smallEta_middle/int_smallEta_all<<"%"<<endl;      

   double int_largeEta_middle = h_virtualVtxdZ_largedEta->Integral(17,21);
   double int_largeEta_all = 1.0*h_virtualVtxdZ_largedEta->GetEntries();
   cout<<"largedEta: "<<int_largeEta_middle<<" out of "<<int_largeEta_all<<"  =  "<<100.0*int_largeEta_middle/int_largeEta_all<<"%"<<endl;      


   TH1F * h_chi2_pho_vtx_smalldEta_Sig = new TH1F("h_chi2_pho_vtx_smalldEta_Sig","h_chi2_pho_vtx_smalldEta_Sig",100,0,50);
   TH1F * h_chi2_pho_vtx_smalldEta_Bkg = new TH1F("h_chi2_pho_vtx_smalldEta_Bkg","h_chi2_pho_vtx_smalldEta_Bkg",100,0,50);
   t_in->Draw("chi2_pho_vtx_gen>>h_chi2_pho_vtx_smalldEta_Sig","abs(pho1Eta-pho2Eta)<0.8 && chi2_pho_vtx_gen < 50");
   t_in->Draw("chi2_pho_vtx>>h_chi2_pho_vtx_smalldEta_Bkg","abs(pho1Eta-pho2Eta)<0.8 && chi2_pho_vtx < 50 && isMatchPv==0 ");
   MaxY=0;

   h_chi2_pho_vtx_smalldEta_Sig->Scale(1.0/h_chi2_pho_vtx_smalldEta_Sig->GetEntries());
   h_chi2_pho_vtx_smalldEta_Bkg->Scale(1.0/h_chi2_pho_vtx_smalldEta_Bkg->GetEntries());

   if(h_chi2_pho_vtx_smalldEta_Sig->GetMaximum()>MaxY) MaxY = h_chi2_pho_vtx_smalldEta_Sig->GetMaximum();
   if(h_chi2_pho_vtx_smalldEta_Bkg->GetMaximum()>MaxY) MaxY = h_chi2_pho_vtx_smalldEta_Bkg->GetMaximum();
   h_chi2_pho_vtx_smalldEta_Bkg->SetLineColor(2);
   h_chi2_pho_vtx_smalldEta_Bkg->SetTitle("");
   h_chi2_pho_vtx_smalldEta_Bkg->GetXaxis()->SetTitle("Vertex #chi^{2}");
   h_chi2_pho_vtx_smalldEta_Bkg->GetYaxis()->SetTitle("pdf(#chi^{2})");
   h_chi2_pho_vtx_smalldEta_Bkg->GetYaxis()->SetTitleOffset(1.5);
   h_chi2_pho_vtx_smalldEta_Bkg->GetYaxis()->SetRangeUser(1.0e-4,1.0);
   h_chi2_pho_vtx_smalldEta_Sig->SetLineColor(4);
   h_chi2_pho_vtx_smalldEta_Bkg->Draw("");
   h_chi2_pho_vtx_smalldEta_Sig->Draw("same");
        TLegend *leg_chi2_pho_vtx_smalldEta = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_chi2_pho_vtx_smalldEta->SetBorderSize(0);
        leg_chi2_pho_vtx_smalldEta->SetTextSize(0.04);
        leg_chi2_pho_vtx_smalldEta->SetLineColor(1);
        leg_chi2_pho_vtx_smalldEta->SetLineStyle(1);
        leg_chi2_pho_vtx_smalldEta->SetLineWidth(1);
        leg_chi2_pho_vtx_smalldEta->SetFillColor(0);
        leg_chi2_pho_vtx_smalldEta->SetFillStyle(1001);
        leg_chi2_pho_vtx_smalldEta->AddEntry(h_chi2_pho_vtx_smalldEta_Sig, "H->#gamma#gamma vertex" ,"l");
        leg_chi2_pho_vtx_smalldEta->AddEntry(h_chi2_pho_vtx_smalldEta_Bkg, "PU vertex" ,"l");
        leg_chi2_pho_vtx_smalldEta->Draw();
   myC->SetLogy(1);
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/chi2_pho_vtx_smalldEta_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/chi2_pho_vtx_smalldEta_"+plot_suffix+".C").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/chi2_pho_vtx_smalldEta_"+plot_suffix+".png").c_str());


   TH1F * h_chi2_pho_vtx_largedEta_Sig = new TH1F("h_chi2_pho_vtx_largedEta_Sig","h_chi2_pho_vtx_largedEta_Sig",100,0,50);
   TH1F * h_chi2_pho_vtx_largedEta_Bkg = new TH1F("h_chi2_pho_vtx_largedEta_Bkg","h_chi2_pho_vtx_largedEta_Bkg",100,0,50);
   t_in->Draw("chi2_pho_vtx_gen>>h_chi2_pho_vtx_largedEta_Sig","abs(pho1Eta-pho2Eta)>0.8 && chi2_pho_vtx_gen < 50");
   t_in->Draw("chi2_pho_vtx>>h_chi2_pho_vtx_largedEta_Bkg","abs(pho1Eta-pho2Eta)>0.8 && chi2_pho_vtx < 50 && isMatchPv==0 ");
   MaxY=0;

   h_chi2_pho_vtx_largedEta_Sig->Scale(1.0/h_chi2_pho_vtx_largedEta_Sig->GetEntries());
   h_chi2_pho_vtx_largedEta_Bkg->Scale(1.0/h_chi2_pho_vtx_largedEta_Bkg->GetEntries());

   if(h_chi2_pho_vtx_largedEta_Sig->GetMaximum()>MaxY) MaxY = h_chi2_pho_vtx_largedEta_Sig->GetMaximum();
   if(h_chi2_pho_vtx_largedEta_Bkg->GetMaximum()>MaxY) MaxY = h_chi2_pho_vtx_largedEta_Bkg->GetMaximum();
   h_chi2_pho_vtx_largedEta_Bkg->SetLineColor(2);
   h_chi2_pho_vtx_largedEta_Bkg->SetTitle("");
   h_chi2_pho_vtx_largedEta_Bkg->GetXaxis()->SetTitle("Vertex #chi^{2}");
   h_chi2_pho_vtx_largedEta_Bkg->GetYaxis()->SetTitle("pdf(#chi^{2})");
   h_chi2_pho_vtx_largedEta_Bkg->GetYaxis()->SetTitleOffset(1.5);
   h_chi2_pho_vtx_largedEta_Bkg->GetYaxis()->SetRangeUser(1.0e-4,1.0);
   h_chi2_pho_vtx_largedEta_Sig->SetLineColor(4);
   h_chi2_pho_vtx_largedEta_Bkg->Draw("");
   h_chi2_pho_vtx_largedEta_Sig->Draw("same");
        TLegend *leg_chi2_pho_vtx_largedEta = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_chi2_pho_vtx_largedEta->SetBorderSize(0);
        leg_chi2_pho_vtx_largedEta->SetTextSize(0.04);
        leg_chi2_pho_vtx_largedEta->SetLineColor(1);
        leg_chi2_pho_vtx_largedEta->SetLineStyle(1);
        leg_chi2_pho_vtx_largedEta->SetLineWidth(1);
        leg_chi2_pho_vtx_largedEta->SetFillColor(0);
        leg_chi2_pho_vtx_largedEta->SetFillStyle(1001);
        leg_chi2_pho_vtx_largedEta->AddEntry(h_chi2_pho_vtx_largedEta_Sig, "H->#gamma#gamma vertex" ,"l");
        leg_chi2_pho_vtx_largedEta->AddEntry(h_chi2_pho_vtx_largedEta_Bkg, "PU vertex" ,"l");
        leg_chi2_pho_vtx_largedEta->Draw();
   myC->SetLogy(1);
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/chi2_pho_vtx_largedEta_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/chi2_pho_vtx_largedEta_"+plot_suffix+".C").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/chi2_pho_vtx_largedEta_"+plot_suffix+".png").c_str());

   TH1F * h_chi2_min_pho_vtx_smalldEta_Sig = new TH1F("h_chi2_min_pho_vtx_smalldEta_Sig","h_chi2_min_pho_vtx_smalldEta_Sig",100,0,50);
   TH1F * h_chi2_min_pho_vtx_smalldEta_Bkg = new TH1F("h_chi2_min_pho_vtx_smalldEta_Bkg","h_chi2_min_pho_vtx_smalldEta_Bkg",100,0,50);
   t_in->Draw("chi2_min_pho_vtx_gen>>h_chi2_min_pho_vtx_smalldEta_Sig","abs(pho1Eta-pho2Eta)<0.8 && chi2_min_pho_vtx_gen < 50");
   t_in->Draw("chi2_min_pho_vtx>>h_chi2_min_pho_vtx_smalldEta_Bkg","abs(pho1Eta-pho2Eta)<0.8 && chi2_min_pho_vtx < 50 && isMatchPv==0 ");
   MaxY=0;

   h_chi2_min_pho_vtx_smalldEta_Sig->Scale(1.0/h_chi2_min_pho_vtx_smalldEta_Sig->GetEntries());
   h_chi2_min_pho_vtx_smalldEta_Bkg->Scale(1.0/h_chi2_min_pho_vtx_smalldEta_Bkg->GetEntries());

   if(h_chi2_min_pho_vtx_smalldEta_Sig->GetMaximum()>MaxY) MaxY = h_chi2_min_pho_vtx_smalldEta_Sig->GetMaximum();
   if(h_chi2_min_pho_vtx_smalldEta_Bkg->GetMaximum()>MaxY) MaxY = h_chi2_min_pho_vtx_smalldEta_Bkg->GetMaximum();
   h_chi2_min_pho_vtx_smalldEta_Bkg->SetLineColor(2);
   h_chi2_min_pho_vtx_smalldEta_Bkg->SetTitle("");
   h_chi2_min_pho_vtx_smalldEta_Bkg->GetXaxis()->SetTitle("Vertex #chi^{2}");
   h_chi2_min_pho_vtx_smalldEta_Bkg->GetYaxis()->SetTitle("pdf(#chi^{2})");
   h_chi2_min_pho_vtx_smalldEta_Bkg->GetYaxis()->SetTitleOffset(1.5);
   h_chi2_min_pho_vtx_smalldEta_Bkg->GetYaxis()->SetRangeUser(1.0e-4,1.0);
   h_chi2_min_pho_vtx_smalldEta_Sig->SetLineColor(4);
   h_chi2_min_pho_vtx_smalldEta_Bkg->Draw("");
   h_chi2_min_pho_vtx_smalldEta_Sig->Draw("same");
        TLegend *leg_chi2_min_pho_vtx_smalldEta = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_chi2_min_pho_vtx_smalldEta->SetBorderSize(0);
        leg_chi2_min_pho_vtx_smalldEta->SetTextSize(0.04);
        leg_chi2_min_pho_vtx_smalldEta->SetLineColor(1);
        leg_chi2_min_pho_vtx_smalldEta->SetLineStyle(1);
        leg_chi2_min_pho_vtx_smalldEta->SetLineWidth(1);
        leg_chi2_min_pho_vtx_smalldEta->SetFillColor(0);
        leg_chi2_min_pho_vtx_smalldEta->SetFillStyle(1001);
        leg_chi2_min_pho_vtx_smalldEta->AddEntry(h_chi2_min_pho_vtx_smalldEta_Sig, "H->#gamma#gamma vertex" ,"l");
        leg_chi2_min_pho_vtx_smalldEta->AddEntry(h_chi2_min_pho_vtx_smalldEta_Bkg, "PU vertex" ,"l");
        leg_chi2_min_pho_vtx_smalldEta->Draw();
   myC->SetLogy(1);
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/chi2_min_pho_vtx_smalldEta_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/chi2_min_pho_vtx_smalldEta_"+plot_suffix+".C").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/chi2_min_pho_vtx_smalldEta_"+plot_suffix+".png").c_str());


   TH1F * h_chi2_min_pho_vtx_largedEta_Sig = new TH1F("h_chi2_min_pho_vtx_largedEta_Sig","h_chi2_min_pho_vtx_largedEta_Sig",100,0,50);
   TH1F * h_chi2_min_pho_vtx_largedEta_Bkg = new TH1F("h_chi2_min_pho_vtx_largedEta_Bkg","h_chi2_min_pho_vtx_largedEta_Bkg",100,0,50);
   t_in->Draw("chi2_min_pho_vtx_gen>>h_chi2_min_pho_vtx_largedEta_Sig","abs(pho1Eta-pho2Eta)>0.8 && chi2_min_pho_vtx_gen < 50");
   t_in->Draw("chi2_min_pho_vtx>>h_chi2_min_pho_vtx_largedEta_Bkg","abs(pho1Eta-pho2Eta)>0.8 && chi2_min_pho_vtx < 50 && isMatchPv==0 ");
   MaxY=0;

   h_chi2_min_pho_vtx_largedEta_Sig->Scale(1.0/h_chi2_min_pho_vtx_largedEta_Sig->GetEntries());
   h_chi2_min_pho_vtx_largedEta_Bkg->Scale(1.0/h_chi2_min_pho_vtx_largedEta_Bkg->GetEntries());

   if(h_chi2_min_pho_vtx_largedEta_Sig->GetMaximum()>MaxY) MaxY = h_chi2_min_pho_vtx_largedEta_Sig->GetMaximum();
   if(h_chi2_min_pho_vtx_largedEta_Bkg->GetMaximum()>MaxY) MaxY = h_chi2_min_pho_vtx_largedEta_Bkg->GetMaximum();
   h_chi2_min_pho_vtx_largedEta_Bkg->SetLineColor(2);
   h_chi2_min_pho_vtx_largedEta_Bkg->SetTitle("");
   h_chi2_min_pho_vtx_largedEta_Bkg->GetXaxis()->SetTitle("Vertex #chi^{2}");
   h_chi2_min_pho_vtx_largedEta_Bkg->GetYaxis()->SetTitle("pdf(#chi^{2})");
   h_chi2_min_pho_vtx_largedEta_Bkg->GetYaxis()->SetTitleOffset(1.5);
   h_chi2_min_pho_vtx_largedEta_Bkg->GetYaxis()->SetRangeUser(1.0e-4,1.0);
   h_chi2_min_pho_vtx_largedEta_Sig->SetLineColor(4);
   h_chi2_min_pho_vtx_largedEta_Bkg->Draw("");
   h_chi2_min_pho_vtx_largedEta_Sig->Draw("same");
        TLegend *leg_chi2_min_pho_vtx_largedEta = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_chi2_min_pho_vtx_largedEta->SetBorderSize(0);
        leg_chi2_min_pho_vtx_largedEta->SetTextSize(0.04);
        leg_chi2_min_pho_vtx_largedEta->SetLineColor(1);
        leg_chi2_min_pho_vtx_largedEta->SetLineStyle(1);
        leg_chi2_min_pho_vtx_largedEta->SetLineWidth(1);
        leg_chi2_min_pho_vtx_largedEta->SetFillColor(0);
        leg_chi2_min_pho_vtx_largedEta->SetFillStyle(1001);
        leg_chi2_min_pho_vtx_largedEta->AddEntry(h_chi2_min_pho_vtx_largedEta_Sig, "H->#gamma#gamma vertex" ,"l");
        leg_chi2_min_pho_vtx_largedEta->AddEntry(h_chi2_min_pho_vtx_largedEta_Bkg, "PU vertex" ,"l");
        leg_chi2_min_pho_vtx_largedEta->Draw();
   myC->SetLogy(1);
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/chi2_min_pho_vtx_largedEta_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/chi2_min_pho_vtx_largedEta_"+plot_suffix+".C").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/chi2_min_pho_vtx_largedEta_"+plot_suffix+".png").c_str());

  //chi2 rank
  TH1F * h_vertex_rank_fraction_smalldEta = new TH1F("h_vertex_rank_fraction_smalldEta","h_vertex_rank_fraction_smalldEta",60,0,60);
  TH1F * h_vertex_min_rank_fraction_smalldEta = new TH1F("h_vertex_min_rank_fraction_smalldEta","h_vertex_min_rank_fraction_smalldEta",60,0,60);
  TH1F * h_vertex_rank_fraction_largedEta = new TH1F("h_vertex_rank_fraction_largedEta","h_vertex_rank_fraction_largedEta",60,0,60);
  TH1F * h_vertex_min_rank_fraction_largedEta = new TH1F("h_vertex_min_rank_fraction_largedEta","h_vertex_min_rank_fraction_largedEta",60,0,60);

  int rank_chi2_genVtx = 0;
  int rank_chi2_min_genVtx = 0;
  int N_Total_smalldEta = 0;
  int N_Total_largedEta = 0;
  int N_genVtx_BelowRank_smalldEta[70] ;	
  int N_genVtx_BelowRank_largedEta[70] ;	
  int N_genVtx_BelowRank_smalldEta_min[70] ;	
  int N_genVtx_BelowRank_largedEta_min[70] ;
  for(int i=0;i<70;i++)
  {
	N_genVtx_BelowRank_smalldEta[i]=0;
	N_genVtx_BelowRank_largedEta[i]=0;
	N_genVtx_BelowRank_smalldEta_min[i]=0;
	N_genVtx_BelowRank_largedEta_min[i]=0;
  }	
  float pho1Eta = 0.0;
  float pho2Eta = 0.0;
  t_in->SetBranchAddress( "pho1Eta",  &pho1Eta);
  t_in->SetBranchAddress( "pho2Eta",  &pho2Eta);
  t_in->SetBranchAddress( "rank_chi2_genVtx",  &rank_chi2_genVtx);
  t_in->SetBranchAddress( "rank_chi2_min_genVtx",  &rank_chi2_min_genVtx);
  int N_Entries = t_in->GetEntries(); 
  for(int i=0;i<N_Entries;i++)
  {
	t_in->GetEntry(i);
	if(std::abs(pho1Eta-pho2Eta)<0.8)
	{
	N_Total_smalldEta ++;
	for(int iR=rank_chi2_genVtx;iR<70;iR++)
	{
	N_genVtx_BelowRank_smalldEta[iR] ++;
	}
	for(int iR=rank_chi2_min_genVtx;iR<70;iR++)
        {
        N_genVtx_BelowRank_smalldEta_min[iR] ++;
        }
	}
	
	else
	{
	N_Total_largedEta ++;
	for(int iR=rank_chi2_genVtx;iR<70;iR++)
	{
	N_genVtx_BelowRank_largedEta[iR] ++;
	}
	for(int iR=rank_chi2_min_genVtx;iR<70;iR++)
        {
        N_genVtx_BelowRank_largedEta_min[iR] ++;
        }
	}
  }
 for(int i=0;i<60;i++)
  {
   h_vertex_rank_fraction_smalldEta->SetBinContent(i+1, (1.0*N_genVtx_BelowRank_smalldEta[i]/(N_Total_smalldEta*1.0)));
   h_vertex_min_rank_fraction_smalldEta->SetBinContent(i+1, (1.0*N_genVtx_BelowRank_smalldEta_min[i]/(N_Total_smalldEta*1.0)));
   h_vertex_rank_fraction_largedEta->SetBinContent(i+1, (1.0*N_genVtx_BelowRank_largedEta[i]/(N_Total_largedEta*1.0)));
   h_vertex_min_rank_fraction_largedEta->SetBinContent(i+1, (1.0*N_genVtx_BelowRank_largedEta_min[i]/(N_Total_largedEta*1.0)));
  }

   h_vertex_rank_fraction_largedEta->SetLineColor(2);
   h_vertex_rank_fraction_largedEta->SetTitle("");
   h_vertex_rank_fraction_largedEta->GetXaxis()->SetTitle("R (vertex rank by #chi^{2})");
   h_vertex_rank_fraction_largedEta->GetYaxis()->SetTitle("Event Fraction with R(H->#gamma#gamma)<R");
   h_vertex_rank_fraction_largedEta->GetYaxis()->SetTitleOffset(1.5);
   h_vertex_rank_fraction_largedEta->GetYaxis()->SetRangeUser(0.0,1.0);
   h_vertex_rank_fraction_smalldEta->SetLineColor(4);
   h_vertex_rank_fraction_largedEta->Draw("");
   h_vertex_rank_fraction_smalldEta->Draw("same");
        TLegend *leg_vertex_rank_fraction = new TLegend(0.58, 0.2, 0.88, 0.4, NULL,"brNDC");
        leg_vertex_rank_fraction->SetBorderSize(0);
        leg_vertex_rank_fraction->SetTextSize(0.04);
        leg_vertex_rank_fraction->SetLineColor(1);
        leg_vertex_rank_fraction->SetLineStyle(1);
        leg_vertex_rank_fraction->SetLineWidth(1);
        leg_vertex_rank_fraction->SetFillColor(0);
        leg_vertex_rank_fraction->SetFillStyle(1001);
        leg_vertex_rank_fraction->AddEntry(h_vertex_rank_fraction_smalldEta, "|#Delta#eta_{#gamma#gamma}|<0.8" ,"l");
        leg_vertex_rank_fraction->AddEntry(h_vertex_rank_fraction_largedEta, "|#Delta#eta_{#gamma#gamma}|>0.8" ,"l");
        leg_vertex_rank_fraction->Draw();
   myC->SetLogy(0);
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/vertex_rank_fraction_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/vertex_rank_fraction_"+plot_suffix+".C").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/vertex_rank_fraction_"+plot_suffix+".png").c_str());

   h_vertex_min_rank_fraction_largedEta->SetLineColor(2);
   h_vertex_min_rank_fraction_largedEta->SetTitle("");
   h_vertex_min_rank_fraction_largedEta->GetXaxis()->SetTitle("R (vertex rank by #chi^{2})");
   h_vertex_min_rank_fraction_largedEta->GetYaxis()->SetTitle("Event Fraction with R(H->#gamma#gamma)<R");
   h_vertex_min_rank_fraction_largedEta->GetYaxis()->SetTitleOffset(1.5);
   h_vertex_min_rank_fraction_largedEta->GetYaxis()->SetRangeUser(0.0,1.0);
   h_vertex_min_rank_fraction_smalldEta->SetLineColor(4);
   h_vertex_min_rank_fraction_largedEta->Draw("");
   h_vertex_min_rank_fraction_smalldEta->Draw("same");
        TLegend *leg_vertex_min_rank_fraction = new TLegend(0.58, 0.2, 0.88, 0.4, NULL,"brNDC");
        leg_vertex_min_rank_fraction->SetBorderSize(0);
        leg_vertex_min_rank_fraction->SetTextSize(0.04);
        leg_vertex_min_rank_fraction->SetLineColor(1);
        leg_vertex_min_rank_fraction->SetLineStyle(1);
        leg_vertex_min_rank_fraction->SetLineWidth(1);
        leg_vertex_min_rank_fraction->SetFillColor(0);
        leg_vertex_min_rank_fraction->SetFillStyle(1001);
        leg_vertex_min_rank_fraction->AddEntry(h_vertex_min_rank_fraction_smalldEta, "|#Delta#eta_{#gamma#gamma}|<0.8" ,"l");
        leg_vertex_min_rank_fraction->AddEntry(h_vertex_min_rank_fraction_largedEta, "|#Delta#eta_{#gamma#gamma}|>0.8" ,"l");
        leg_vertex_min_rank_fraction->Draw();
   myC->SetLogy(0);
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/vertex_min_rank_fraction_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/vertex_min_rank_fraction_"+plot_suffix+".C").c_str());
   myC->SaveAs(("~/www/sharebox/tomyself/tmp/vertex_min_rank_fraction_"+plot_suffix+".png").c_str());



}
