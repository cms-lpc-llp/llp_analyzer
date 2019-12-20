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


void trackdT(string inputFilename="../../HggRazorUpgradeTiming_PU140_NoTiming.root", string plot_suffix="PU140_NoTiming")
{


	TFile f_in(inputFilename.c_str());
	
	TTree *t_in = (TTree*)f_in.Get("HggRazor");


	double MaxY=0;
   	TCanvas * myC = new TCanvas("c1","c1",100,100,800,700);
   
        myC->cd();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);


	gStyle->SetOptStat(0);
	
   TH1F * h_allTrackdT_Sig = new TH1F("h_allTrackdT_Sig","h_allTrackdT_Sig",600,-0.3,0.3);
   TH1F * h_allTrackdT_Bkg = new TH1F("h_allTrackdT_Bkg","h_allTrackdT_Bkg",600,-0.3,0.3);
   t_in->Draw("allTrackdT_Sig>>h_allTrackdT_Sig");
   t_in->Draw("allTrackdT_Bkg>>h_allTrackdT_Bkg");
   MaxY=0;
   h_allTrackdT_Sig->Scale(1.0/h_allTrackdT_Sig->GetEntries());
   h_allTrackdT_Bkg->Scale(1.0/h_allTrackdT_Bkg->GetEntries());

   if(h_allTrackdT_Sig->GetMaximum()>MaxY) MaxY = h_allTrackdT_Sig->GetMaximum();
   if(h_allTrackdT_Bkg->GetMaximum()>MaxY) MaxY = h_allTrackdT_Bkg->GetMaximum();
   h_allTrackdT_Bkg->SetLineColor(2);
   h_allTrackdT_Bkg->SetTitle("");
   h_allTrackdT_Bkg->GetXaxis()->SetTitle("t_{track} - t_{vertex} / ns");
   h_allTrackdT_Bkg->GetYaxis()->SetTitle("Events/0.001ns");
   h_allTrackdT_Bkg->GetYaxis()->SetTitleOffset(1.8);
   h_allTrackdT_Bkg->GetYaxis()->SetRangeUser(0.0,1.2*MaxY);
   h_allTrackdT_Sig->SetLineColor(4);
   h_allTrackdT_Bkg->Draw();
   h_allTrackdT_Sig->Draw("same");
        TLegend *leg_allTrackdT = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_allTrackdT->SetBorderSize(0);
        leg_allTrackdT->SetTextSize(0.04);
        leg_allTrackdT->SetLineColor(1);
        leg_allTrackdT->SetLineStyle(1);
        leg_allTrackdT->SetLineWidth(1);
        leg_allTrackdT->SetFillColor(0);
        leg_allTrackdT->SetFillStyle(1001);
        leg_allTrackdT->AddEntry(h_allTrackdT_Bkg, "non-gen-matched vtx" ,"l");
        leg_allTrackdT->AddEntry(h_allTrackdT_Sig, "gen-matched vtx" ,"l");
        leg_allTrackdT->Draw();

   myC->SetLogy(0);
   myC->SaveAs(("plots/allTrackdT_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("plots/allTrackdT_"+plot_suffix+".png").c_str());

   double int_Sig_middle = h_allTrackdT_Sig->Integral(240,360);
   double int_Sig_middleL = h_allTrackdT_Sig->Integral(240,300);
   double int_Bkg_middle = h_allTrackdT_Bkg->Integral(240,360);
   double int_Bkg_middleL = h_allTrackdT_Bkg->Integral(240,300);

   double int_Sig_all = h_allTrackdT_Sig->Integral(1,600);
   double int_Sig_allL = h_allTrackdT_Sig->Integral(1,300);
   double int_Bkg_all = h_allTrackdT_Bkg->Integral(1,600);
   double int_Bkg_allL = h_allTrackdT_Bkg->Integral(1,300);

   cout<<"all Pt:"<<endl;
   cout<<"gen-matched: "<<int_Sig_middle<<" out of "<<int_Sig_all<<"  =  "<<100.0*int_Sig_middle/int_Sig_all<<"%"<<endl;	
   cout<<"non-gen-matched: "<<int_Bkg_middle<<" out of "<<int_Bkg_all<<"  =  "<<100.0*int_Bkg_middle/int_Bkg_all<<"%"<<endl;	

   cout<<"gen-matched: "<<int_Sig_middleL<<" out of "<<int_Sig_allL<<"  =  "<<100.0*int_Sig_middleL/int_Sig_allL<<"%"<<endl;	

   cout<<"non-gen-matched: "<<int_Bkg_middleL<<" out of "<<int_Bkg_allL<<"  =  "<<100.0*int_Bkg_middleL/int_Bkg_allL<<"%"<<endl;	


   TH1F * h_allTrackdT_NoScale_Sig = new TH1F("h_allTrackdT_NoScale_Sig","h_allTrackdT_NoScale_Sig",600,-0.3,0.3);
   TH1F * h_allTrackdT_NoScale_Bkg = new TH1F("h_allTrackdT_NoScale_Bkg","h_allTrackdT_NoScale_Bkg",600,-0.3,0.3);
   t_in->Draw("allTrackdT_Sig>>h_allTrackdT_NoScale_Sig");
   t_in->Draw("allTrackdT_Bkg>>h_allTrackdT_NoScale_Bkg");
   MaxY=0;
//   h_allTrackdT_NoScale_Sig->Scale(1.0/h_allTrackdT_NoScale_Sig->GetEntries());
//   h_allTrackdT_NoScale_Bkg->Scale(1.0/h_allTrackdT_NoScale_Bkg->GetEntries());

   if(h_allTrackdT_NoScale_Sig->GetMaximum()>MaxY) MaxY = h_allTrackdT_NoScale_Sig->GetMaximum();
   if(h_allTrackdT_NoScale_Bkg->GetMaximum()>MaxY) MaxY = h_allTrackdT_NoScale_Bkg->GetMaximum();
   h_allTrackdT_NoScale_Bkg->SetLineColor(2);
   h_allTrackdT_NoScale_Bkg->SetTitle("");
   h_allTrackdT_NoScale_Bkg->GetXaxis()->SetTitle("t_{track} - t_{vertex} / ns");
   h_allTrackdT_NoScale_Bkg->GetYaxis()->SetTitle("Events/0.001ns");
   h_allTrackdT_NoScale_Bkg->GetYaxis()->SetTitleOffset(1.8);
   h_allTrackdT_NoScale_Bkg->GetYaxis()->SetRangeUser(0.0,1.2*MaxY);
   h_allTrackdT_NoScale_Sig->SetLineColor(4);
   h_allTrackdT_NoScale_Bkg->Draw();
   h_allTrackdT_NoScale_Sig->Draw("same");
        TLegend *leg_allTrackdT_NoScale = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_allTrackdT_NoScale->SetBorderSize(0);
        leg_allTrackdT_NoScale->SetTextSize(0.04);
        leg_allTrackdT_NoScale->SetLineColor(1);
        leg_allTrackdT_NoScale->SetLineStyle(1);
        leg_allTrackdT_NoScale->SetLineWidth(1);
        leg_allTrackdT_NoScale->SetFillColor(0);
        leg_allTrackdT_NoScale->SetFillStyle(1001);
        leg_allTrackdT_NoScale->AddEntry(h_allTrackdT_NoScale_Bkg, "non-gen-matched vtx" ,"l");
        leg_allTrackdT_NoScale->AddEntry(h_allTrackdT_NoScale_Sig, "gen-matched vtx" ,"l");
        leg_allTrackdT_NoScale->Draw();

   myC->SetLogy(0);
   myC->SaveAs(("plots/allTrackdT_NoScale_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("plots/allTrackdT_NoScale_"+plot_suffix+".png").c_str());

    int_Sig_middle = h_allTrackdT_NoScale_Sig->Integral(240,360);
    int_Sig_middleL = h_allTrackdT_NoScale_Sig->Integral(240,300);
    int_Bkg_middle = h_allTrackdT_NoScale_Bkg->Integral(240,360);
    int_Bkg_middleL = h_allTrackdT_NoScale_Bkg->Integral(240,300);

    int_Sig_all = h_allTrackdT_NoScale_Sig->Integral(1,600);
    int_Sig_allL = h_allTrackdT_NoScale_Sig->Integral(1,300);
    int_Bkg_all = h_allTrackdT_NoScale_Bkg->Integral(1,600);
    int_Bkg_allL = h_allTrackdT_NoScale_Bkg->Integral(1,300);

    
   cout<<"all Pt no scale:"<<endl;
   cout<<"gen-matched: "<<int_Sig_middle<<" out of "<<int_Sig_all<<"  =  "<<100.0*int_Sig_middle/int_Sig_all<<"%"<<endl;	
   cout<<"non-gen-matched: "<<int_Bkg_middle<<" out of "<<int_Bkg_all<<"  =  "<<100.0*int_Bkg_middle/int_Bkg_all<<"%"<<endl;	

   cout<<"gen-matched: "<<int_Sig_middleL<<" out of "<<int_Sig_allL<<"  =  "<<100.0*int_Sig_middleL/int_Sig_allL<<"%"<<endl;	

   cout<<"non-gen-matched: "<<int_Bkg_middleL<<" out of "<<int_Bkg_allL<<"  =  "<<100.0*int_Bkg_middleL/int_Bkg_allL<<"%"<<endl;	
 



   TH1F * h_allTrackdT_highPt_Sig = new TH1F("h_allTrackdT_highPt_Sig","h_allTrackdT_highPt_Sig",600,-0.3,0.3);
   TH1F * h_allTrackdT_highPt_Bkg = new TH1F("h_allTrackdT_highPt_Bkg","h_allTrackdT_highPt_Bkg",600,-0.3,0.3);
   t_in->Draw("allTrackdT_Sig>>h_allTrackdT_highPt_Sig","allTrackPt_Sig>10.0");
   t_in->Draw("allTrackdT_Bkg>>h_allTrackdT_highPt_Bkg","allTrackPt_Bkg>10.0");
   MaxY=0;
   h_allTrackdT_highPt_Sig->Scale(1.0/h_allTrackdT_highPt_Sig->GetEntries());
   h_allTrackdT_highPt_Bkg->Scale(1.0/h_allTrackdT_highPt_Bkg->GetEntries());

   if(h_allTrackdT_highPt_Sig->GetMaximum()>MaxY) MaxY = h_allTrackdT_highPt_Sig->GetMaximum();
   if(h_allTrackdT_highPt_Bkg->GetMaximum()>MaxY) MaxY = h_allTrackdT_highPt_Bkg->GetMaximum();
   h_allTrackdT_highPt_Bkg->SetLineColor(2);
   h_allTrackdT_highPt_Bkg->SetTitle("");
   h_allTrackdT_highPt_Bkg->GetXaxis()->SetTitle("t_{track} - t_{vertex} / ns");
   h_allTrackdT_highPt_Bkg->GetYaxis()->SetTitle("Events/0.001ns");
   h_allTrackdT_highPt_Bkg->GetYaxis()->SetTitleOffset(1.8);
   h_allTrackdT_highPt_Bkg->GetYaxis()->SetRangeUser(0.0,1.2*MaxY);
   h_allTrackdT_highPt_Sig->SetLineColor(4);
   h_allTrackdT_highPt_Bkg->Draw();
   h_allTrackdT_highPt_Sig->Draw("same");
        TLegend *leg_allTrackdT_highPt = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_allTrackdT_highPt->SetBorderSize(0);
        leg_allTrackdT_highPt->SetTextSize(0.04);
        leg_allTrackdT_highPt->SetLineColor(1);
        leg_allTrackdT_highPt->SetLineStyle(1);
        leg_allTrackdT_highPt->SetLineWidth(1);
        leg_allTrackdT_highPt->SetFillColor(0);
        leg_allTrackdT_highPt->SetFillStyle(1001);
        leg_allTrackdT_highPt->AddEntry(h_allTrackdT_highPt_Bkg, "non-gen-matched vtx" ,"l");
        leg_allTrackdT_highPt->AddEntry(h_allTrackdT_highPt_Sig, "gen-matched vtx" ,"l");
        leg_allTrackdT_highPt->Draw();

   myC->SetLogy(0);
   myC->SaveAs(("plots/allTrackdT_highPt_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("plots/allTrackdT_highPt_"+plot_suffix+".png").c_str());

   int_Sig_middle = h_allTrackdT_highPt_Sig->Integral(240,360);
   int_Bkg_middle = h_allTrackdT_highPt_Bkg->Integral(240,360);

   int_Sig_all = h_allTrackdT_highPt_Sig->Integral(1,600);
   int_Bkg_all = h_allTrackdT_highPt_Bkg->Integral(1,600);

   cout<<"Pt>10 GeV: "<<endl;
   cout<<"gen-matched: "<<int_Sig_middle<<" out of "<<int_Sig_all<<"  =  "<<100.0*int_Sig_middle/int_Sig_all<<"%"<<endl;	
   cout<<"non-gen-matched: "<<int_Bkg_middle<<" out of "<<int_Bkg_all<<"  =  "<<100.0*int_Bkg_middle/int_Bkg_all<<"%"<<endl;	


   TH1F * h_allTrackdT_NoScale_highPt_Sig = new TH1F("h_allTrackdT_NoScale_highPt_Sig","h_allTrackdT_NoScale_highPt_Sig",600,-0.3,0.3);
   TH1F * h_allTrackdT_NoScale_highPt_Bkg = new TH1F("h_allTrackdT_NoScale_highPt_Bkg","h_allTrackdT_NoScale_highPt_Bkg",600,-0.3,0.3);
   t_in->Draw("allTrackdT_Sig>>h_allTrackdT_NoScale_highPt_Sig","allTrackPt_Sig>10.0");
   t_in->Draw("allTrackdT_Bkg>>h_allTrackdT_NoScale_highPt_Bkg","allTrackPt_Bkg>10.0");
   MaxY=0;
   //h_allTrackdT_NoScale_highPt_Sig->Scale(1.0/h_allTrackdT_NoScale_highPt_Sig->GetEntries());
   //h_allTrackdT_NoScale_highPt_Bkg->Scale(1.0/h_allTrackdT_NoScale_highPt_Bkg->GetEntries());

   if(h_allTrackdT_NoScale_highPt_Sig->GetMaximum()>MaxY) MaxY = h_allTrackdT_NoScale_highPt_Sig->GetMaximum();
   if(h_allTrackdT_NoScale_highPt_Bkg->GetMaximum()>MaxY) MaxY = h_allTrackdT_NoScale_highPt_Bkg->GetMaximum();
   h_allTrackdT_NoScale_highPt_Bkg->SetLineColor(2);
   h_allTrackdT_NoScale_highPt_Bkg->SetTitle("");
   h_allTrackdT_NoScale_highPt_Bkg->GetXaxis()->SetTitle("t_{track} - t_{vertex} / ns");
   h_allTrackdT_NoScale_highPt_Bkg->GetYaxis()->SetTitle("Events/0.001ns");
   h_allTrackdT_NoScale_highPt_Bkg->GetYaxis()->SetTitleOffset(1.8);
   h_allTrackdT_NoScale_highPt_Bkg->GetYaxis()->SetRangeUser(0.0,1.2*MaxY);
   h_allTrackdT_NoScale_highPt_Sig->SetLineColor(4);
   h_allTrackdT_NoScale_highPt_Bkg->Draw();
   h_allTrackdT_NoScale_highPt_Sig->Draw("same");
        TLegend *leg_allTrackdT_NoScale_highPt = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_allTrackdT_NoScale_highPt->SetBorderSize(0);
        leg_allTrackdT_NoScale_highPt->SetTextSize(0.04);
        leg_allTrackdT_NoScale_highPt->SetLineColor(1);
        leg_allTrackdT_NoScale_highPt->SetLineStyle(1);
        leg_allTrackdT_NoScale_highPt->SetLineWidth(1);
        leg_allTrackdT_NoScale_highPt->SetFillColor(0);
        leg_allTrackdT_NoScale_highPt->SetFillStyle(1001);
        leg_allTrackdT_NoScale_highPt->AddEntry(h_allTrackdT_NoScale_highPt_Bkg, "non-gen-matched vtx" ,"l");
        leg_allTrackdT_NoScale_highPt->AddEntry(h_allTrackdT_NoScale_highPt_Sig, "gen-matched vtx" ,"l");
        leg_allTrackdT_NoScale_highPt->Draw();

   myC->SetLogy(0);
   myC->SaveAs(("plots/allTrackdT_NoScale_highPt_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("plots/allTrackdT_NoScale_highPt_"+plot_suffix+".png").c_str());

   int_Sig_middle = h_allTrackdT_NoScale_highPt_Sig->Integral(240,360);
   int_Bkg_middle = h_allTrackdT_NoScale_highPt_Bkg->Integral(240,360);

   int_Sig_all = h_allTrackdT_NoScale_highPt_Sig->Integral(1,600);
   int_Bkg_all = h_allTrackdT_NoScale_highPt_Bkg->Integral(1,600);

   cout<<"Pt>10 GeV - no scale: "<<endl;
   cout<<"gen-matched: "<<int_Sig_middle<<" out of "<<int_Sig_all<<"  =  "<<100.0*int_Sig_middle/int_Sig_all<<"%"<<endl;	
   cout<<"non-gen-matched: "<<int_Bkg_middle<<" out of "<<int_Bkg_all<<"  =  "<<100.0*int_Bkg_middle/int_Bkg_all<<"%"<<endl;	



   TH1F * h_allTrackdT_NoScale_PtWeight_Sig = new TH1F("h_allTrackdT_NoScale_PtWeight_Sig","h_allTrackdT_NoScale_PtWeight_Sig",600,-0.3,0.3);
   TH1F * h_allTrackdT_NoScale_PtWeight_Bkg = new TH1F("h_allTrackdT_NoScale_PtWeight_Bkg","h_allTrackdT_NoScale_PtWeight_Bkg",600,-0.3,0.3);
   t_in->Draw("allTrackdT_Sig>>h_allTrackdT_NoScale_PtWeight_Sig","allTrackPt_Sig*allTrackPt_Sig*1.0");
   t_in->Draw("allTrackdT_Bkg>>h_allTrackdT_NoScale_PtWeight_Bkg","allTrackPt_Bkg*allTrackPt_Bkg*1.0");
   MaxY=0;
   //h_allTrackdT_NoScale_PtWeight_Sig->Scale(1.0/h_allTrackdT_NoScale_PtWeight_Sig->GetEntries());
   //h_allTrackdT_NoScale_PtWeight_Bkg->Scale(1.0/h_allTrackdT_NoScale_PtWeight_Bkg->GetEntries());

   if(h_allTrackdT_NoScale_PtWeight_Sig->GetMaximum()>MaxY) MaxY = h_allTrackdT_NoScale_PtWeight_Sig->GetMaximum();
   if(h_allTrackdT_NoScale_PtWeight_Bkg->GetMaximum()>MaxY) MaxY = h_allTrackdT_NoScale_PtWeight_Bkg->GetMaximum();
   h_allTrackdT_NoScale_PtWeight_Bkg->SetLineColor(2);
   h_allTrackdT_NoScale_PtWeight_Bkg->SetTitle("");
   h_allTrackdT_NoScale_PtWeight_Bkg->GetXaxis()->SetTitle("t_{track} - t_{vertex} / ns");
   h_allTrackdT_NoScale_PtWeight_Bkg->GetYaxis()->SetTitle("#Sigma P^{2}_{t-track}/GeV^{2}/0.001ns");
   h_allTrackdT_NoScale_PtWeight_Bkg->GetYaxis()->SetTitleOffset(1.8);
   //h_allTrackdT_NoScale_PtWeight_Bkg->GetYaxis()->SetRangeUser(1.0,1.2*MaxY);
   h_allTrackdT_NoScale_PtWeight_Bkg->GetYaxis()->SetRangeUser(1.0,1.3e5);
   h_allTrackdT_NoScale_PtWeight_Sig->SetLineColor(4);
   h_allTrackdT_NoScale_PtWeight_Bkg->Draw("h");
   h_allTrackdT_NoScale_PtWeight_Sig->Draw("hsame");
        TLegend *leg_allTrackdT_NoScale_PtWeight = new TLegend(0.5, 0.78, 0.8, 0.89, NULL,"brNDC");
        leg_allTrackdT_NoScale_PtWeight->SetBorderSize(0);
        leg_allTrackdT_NoScale_PtWeight->SetTextSize(0.04);
        leg_allTrackdT_NoScale_PtWeight->SetLineColor(1);
        leg_allTrackdT_NoScale_PtWeight->SetLineStyle(1);
        leg_allTrackdT_NoScale_PtWeight->SetLineWidth(1);
        leg_allTrackdT_NoScale_PtWeight->SetFillColor(0);
        leg_allTrackdT_NoScale_PtWeight->SetFillStyle(1001);
        leg_allTrackdT_NoScale_PtWeight->AddEntry(h_allTrackdT_NoScale_PtWeight_Bkg, "non-gen-matched vtx" ,"l");
        leg_allTrackdT_NoScale_PtWeight->AddEntry(h_allTrackdT_NoScale_PtWeight_Sig, "gen-matched vtx" ,"l");
        leg_allTrackdT_NoScale_PtWeight->Draw();

   myC->SetLogy(0);
   myC->SaveAs(("plots/allTrackdT_NoScale_PtWeight_"+plot_suffix+".pdf").c_str());
   myC->SaveAs(("plots/allTrackdT_NoScale_PtWeight_"+plot_suffix+".png").c_str());

   int_Sig_middle = h_allTrackdT_NoScale_PtWeight_Sig->Integral(240,360);
   int_Bkg_middle = h_allTrackdT_NoScale_PtWeight_Bkg->Integral(240,360);

   int_Sig_all = h_allTrackdT_NoScale_PtWeight_Sig->Integral(1,600);
   int_Bkg_all = h_allTrackdT_NoScale_PtWeight_Bkg->Integral(1,600);

   cout<<"Pt>10 GeV - no scale: "<<endl;
   cout<<"gen-matched: "<<int_Sig_middle<<" out of "<<int_Sig_all<<"  =  "<<100.0*int_Sig_middle/int_Sig_all<<"%"<<endl;	
   cout<<"non-gen-matched: "<<int_Bkg_middle<<" out of "<<int_Bkg_all<<"  =  "<<100.0*int_Bkg_middle/int_Bkg_all<<"%"<<endl;	
 

   TH2F * h_track_Pt_deltaT_Sig = new TH2F("h_track_Pt_deltaT_Sig","h_track_Pt_deltaT_Sig",600,-0.3,0.3,600,0,30); 
   TH2F * h_track_Pt_deltaT_Bkg = new TH2F("h_track_Pt_deltaT_Bkg","h_track_Pt_deltaT_Bkg",600,-0.3,0.3,600,0,30); 
   t_in->Draw("allTrackPt_Sig:allTrackdT_Sig>>h_track_Pt_deltaT_Sig");
   t_in->Draw("allTrackPt_Bkg:allTrackdT_Bkg>>h_track_Pt_deltaT_Bkg");
   h_track_Pt_deltaT_Bkg->SetTitle(""); 
   h_track_Pt_deltaT_Bkg->GetXaxis()->SetTitle("t_{track} - t_{vertex} / ns"); 
   h_track_Pt_deltaT_Bkg->GetYaxis()->SetTitle("P_{t-track} / GeV"); 
   h_track_Pt_deltaT_Bkg->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/track_Pt_deltaT_"+plot_suffix+"_Bkg.pdf").c_str());    
   myC->SaveAs(("plots/track_Pt_deltaT_"+plot_suffix+"_Bkg.png").c_str());    
   h_track_Pt_deltaT_Sig->SetTitle(""); 
   h_track_Pt_deltaT_Sig->GetXaxis()->SetTitle("t_{track} - t_{vertex} / ns"); 
   h_track_Pt_deltaT_Sig->GetYaxis()->SetTitle("P_{t-track} / GeV"); 
   h_track_Pt_deltaT_Sig->Draw("colz");
   myC->SetLogy(0);
   myC->SaveAs(("plots/track_Pt_deltaT_"+plot_suffix+"_Sig.pdf").c_str());    
   myC->SaveAs(("plots/track_Pt_deltaT_"+plot_suffix+"_Sig.png").c_str()); 


}
