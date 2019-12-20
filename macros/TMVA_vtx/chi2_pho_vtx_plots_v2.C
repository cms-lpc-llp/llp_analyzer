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


void chi2_pho_vtx_plots_v2()//(string inputFilename="../../HggRazorUpgradeTiming_PU0_Timing.root", string plot_suffix="PU0_Timing")
{
	string plot_suffix = "PU0_vs_PU140";
	
	TFile f_in_PU0("../../HggRazorUpgradeTiming_PU0_Timing.root");
	TFile f_in_PU140("../../HggRazorUpgradeTiming_PU140_Timing.root");
	TTree *t_in_PU0 = (TTree*)f_in_PU0.Get("HggRazor");
	TTree *t_in_PU140 = (TTree*)f_in_PU140.Get("HggRazor");


	double MaxY=0;
   	TCanvas * myC = new TCanvas("c1","c1",100,100,800,700);
   
        myC->cd();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);

	gStyle->SetOptStat(0);


	
   TH1F * h_virtualVtxdZ_smalldEta = new TH1F("h_virtualVtxdZ_smalldEta","h_virtualVtxdZ_smalldEta",40,-100,100);
   TH1F * h_virtualVtxdZ_largedEta = new TH1F("h_virtualVtxdZ_largedEta","h_virtualVtxdZ_largedEta",40,-100,100);
   t_in_PU0->Draw("(10.0*virtualVtxdZ)>>h_virtualVtxdZ_smalldEta","abs(pho1Eta-pho2Eta)<0.8");
   t_in_PU0->Draw("(10.0*virtualVtxdZ)>>h_virtualVtxdZ_largedEta","abs(pho1Eta-pho2Eta)>0.8");
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
   myC->SaveAs(("plots/virtualVtxdZ_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("plots/virtualVtxdZ_"+plot_suffix+".png").c_str());

   double int_smallEta_middle = h_virtualVtxdZ_smalldEta->Integral(17,21);
   double int_smallEta_all = 1.0*h_virtualVtxdZ_smalldEta->GetEntries();
   cout<<"smalldEta: "<<int_smallEta_middle<<" out of "<<int_smallEta_all<<"  =  "<<100.0*int_smallEta_middle/int_smallEta_all<<"%"<<endl;      

   double int_largeEta_middle = h_virtualVtxdZ_largedEta->Integral(17,21);
   double int_largeEta_all = 1.0*h_virtualVtxdZ_largedEta->GetEntries();
   cout<<"largedEta: "<<int_largeEta_middle<<" out of "<<int_largeEta_all<<"  =  "<<100.0*int_largeEta_middle/int_largeEta_all<<"%"<<endl;      


   TH1F * h_chi2_pho_vtx_smalldEta_Sig = new TH1F("h_chi2_pho_vtx_smalldEta_Sig","h_chi2_pho_vtx_smalldEta_Sig",100,0,50);
   TH1F * h_chi2_pho_vtx_smalldEta_Bkg = new TH1F("h_chi2_pho_vtx_smalldEta_Bkg","h_chi2_pho_vtx_smalldEta_Bkg",100,0,50);
   t_in_PU0->Draw("chi2_pho_vtx_gen>>h_chi2_pho_vtx_smalldEta_Sig","abs(pho1Eta-pho2Eta)<0.8 && chi2_pho_vtx_gen < 50");
   t_in_PU140->Draw("chi2_pho_vtx>>h_chi2_pho_vtx_smalldEta_Bkg","abs(pho1Eta-pho2Eta)<0.8 && chi2_pho_vtx < 50 && isMatchPv==0 ");
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
   h_chi2_pho_vtx_smalldEta_Bkg->GetYaxis()->SetRangeUser(1.0e-3,1.0);
   h_chi2_pho_vtx_smalldEta_Sig->SetLineColor(4);
   h_chi2_pho_vtx_smalldEta_Bkg->Draw();
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
   myC->SaveAs(("plots/chi2_pho_vtx_smalldEta_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("plots/chi2_pho_vtx_smalldEta_"+plot_suffix+".png").c_str());

}
