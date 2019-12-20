
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include <string>
#include <iostream>
#include <math.h>



void vtxEff_Plot()
{

   TFile *f_NoTiming = new TFile("../../retrain_NoTiming/HggRazorUpgradeTiming_PU140_NoTiming.root","READ");
   TFile *f_TrackVertexTimingOnly = new TFile("../../retrain_TrackVertexTimingOnly/HggRazorUpgradeTiming_PU140_Timing.root","READ");
   TFile *f_PhotonTimingOnly = new TFile("../../retrain_PhotonTimingOnly/HggRazorUpgradeTiming_PU140_NoTiming.root","READ");
   TFile *f_TrackVertexPhotonTiming = new TFile("../../retrain_TrackVertexPhotonTiming/HggRazorUpgradeTiming_PU140_Timing.root","READ");

   TTree *tree_NoTiming = (TTree*)f_NoTiming->Get("HggRazor");
   TTree *tree_TrackVertexTimingOnly = (TTree*)f_TrackVertexTimingOnly->Get("HggRazor");
   TTree *tree_PhotonTimingOnly = (TTree*)f_PhotonTimingOnly->Get("HggRazor");
   TTree *tree_TrackVertexPhotonTiming = (TTree*)f_TrackVertexPhotonTiming->Get("HggRazor");

   int NvtxTotal[4];
   int NvtxSelected[4];
   int NvtxSelected2[4];
   for(int i=0;i<4;i++)
	{
	NvtxTotal[i] = 0;
	NvtxSelected[i] = 0;
	NvtxSelected2[i] = 0;
	}

   TH1F * h_NoTiming = new TH1F("h_NoTiming","h_NoTiming",20,0,20); 
   tree_NoTiming->Draw("(10.0*vtxdZ)>>h_NoTiming");  
   NvtxTotal[0] = tree_NoTiming->GetEntries();
   NvtxSelected[0] = tree_NoTiming->GetEntries("vtxdZ<1.0");
   NvtxSelected2[0] = tree_NoTiming->GetEntries("vtxdZ<0.2");

   TH1F * h_TrackVertexTimingOnly = new TH1F("h_TrackVertexTimingOnly","h_TrackVertexTimingOnly",20,0,20); 
   tree_TrackVertexTimingOnly->Draw("(10.0*vtxdZ)>>h_TrackVertexTimingOnly");  
   NvtxTotal[1] = tree_TrackVertexTimingOnly->GetEntries();
   NvtxSelected[1] = tree_TrackVertexTimingOnly->GetEntries("vtxdZ<1.0");
   NvtxSelected2[1] = tree_TrackVertexTimingOnly->GetEntries("vtxdZ<0.2");

   TH1F * h_PhotonTimingOnly = new TH1F("h_PhotonTimingOnly","h_PhotonTimingOnly",20,0,20); 
   tree_PhotonTimingOnly->Draw("(10.0*vtxdZ)>>h_PhotonTimingOnly");  
   NvtxTotal[2] = tree_PhotonTimingOnly->GetEntries();
   NvtxSelected[2] = tree_PhotonTimingOnly->GetEntries("vtxdZ<1.0");
   NvtxSelected2[2] = tree_PhotonTimingOnly->GetEntries("vtxdZ<0.2");

   TH1F * h_TrackVertexPhotonTiming = new TH1F("h_TrackVertexPhotonTiming","h_TrackVertexPhotonTiming",20,0,20); 
   tree_TrackVertexPhotonTiming->Draw("(10.0*vtxdZ)>>h_TrackVertexPhotonTiming");  
   NvtxTotal[3] = tree_TrackVertexPhotonTiming->GetEntries();
   NvtxSelected[3] = tree_TrackVertexPhotonTiming->GetEntries("vtxdZ<1.0");
   NvtxSelected2[3] = tree_TrackVertexPhotonTiming->GetEntries("vtxdZ<0.2");


TCanvas * myC = new TCanvas("c1","c1",100,100,800,700);

        myC->cd();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);

        gStyle->SetOptStat(0);

   double MaxY=0;
   if(h_NoTiming->GetMaximum()>MaxY) MaxY = h_NoTiming->GetMaximum();
   if(h_TrackVertexTimingOnly->GetMaximum()>MaxY) MaxY = h_TrackVertexTimingOnly->GetMaximum();
   if(h_PhotonTimingOnly->GetMaximum()>MaxY) MaxY = h_PhotonTimingOnly->GetMaximum();
   if(h_TrackVertexPhotonTiming->GetMaximum()>MaxY) MaxY = h_TrackVertexPhotonTiming->GetMaximum();

   h_NoTiming->Scale(1.0/NvtxTotal[0]);
   h_TrackVertexTimingOnly->Scale(1.0/NvtxTotal[1]);
   h_PhotonTimingOnly->Scale(1.0/NvtxTotal[2]);
   h_TrackVertexPhotonTiming->Scale(1.0/NvtxTotal[3]);

   h_NoTiming->SetLineColor(1);
   h_NoTiming->SetTitle("");
   h_NoTiming->GetXaxis()->SetTitle("|z - z_{gen}| / mm");
   h_NoTiming->GetYaxis()->SetTitle("Event pdf/1mm");
   h_NoTiming->GetYaxis()->SetTitleOffset(1.5);
   h_NoTiming->GetYaxis()->SetRangeUser(0,1.0);
   h_NoTiming->Draw();
   h_TrackVertexTimingOnly->SetLineColor(2);
   h_PhotonTimingOnly->SetLineColor(4);
   h_TrackVertexPhotonTiming->SetLineColor(3);

   h_TrackVertexTimingOnly->Draw("same");
   h_PhotonTimingOnly->Draw("same");
   h_TrackVertexPhotonTiming->Draw("same");

        TLegend *leg_vtxdZ = new TLegend(0.25, 0.63, 0.7, 0.89, NULL,"brNDC");
        leg_vtxdZ->SetBorderSize(0);
        leg_vtxdZ->SetTextSize(0.04);
        leg_vtxdZ->SetLineColor(1);
        leg_vtxdZ->SetLineStyle(1);
        leg_vtxdZ->SetLineWidth(1);
        leg_vtxdZ->SetFillColor(0);
        leg_vtxdZ->SetFillStyle(1001);
        leg_vtxdZ->AddEntry(h_NoTiming, "No Timing" ,"l");
        leg_vtxdZ->AddEntry(h_TrackVertexTimingOnly, "Track/Vertex Timing Only" ,"l");
        leg_vtxdZ->AddEntry(h_PhotonTimingOnly, "Photon Timing Only" ,"l");
        leg_vtxdZ->AddEntry(h_TrackVertexPhotonTiming, "Track/Vertex and Photon Timing" ,"l");
        leg_vtxdZ->Draw();
   myC->SetLogy(0);
   myC->SetLogx(0);
   myC->SaveAs("~/www/sharebox/tomyself/tmp/vtxdZ_PU140.pdf");
   myC->SaveAs("~/www/sharebox/tomyself/tmp/vtxdZ_PU140.png");
 
   FILE* m_outfile = fopen("vtx_eff_table.tex", "a");

   fprintf(m_outfile, "retrain (10mm)...");

   for(int i=0;i<4;i++)
	{
	fprintf(m_outfile, "& %4.1f\\%% ", (100.0*NvtxSelected[i])/(NvtxTotal[i]*1.0));
	}
   fprintf(m_outfile, " \\\\ \\hline \n");
   fprintf(m_outfile, "retrain (2mm)...");
   for(int i=0;i<4;i++)
	{
	fprintf(m_outfile, "& %4.1f\\%% ", (100.0*NvtxSelected2[i])/(NvtxTotal[i]*1.0));
	}

   fprintf(m_outfile, " \\\\ \\hline \n");
}
