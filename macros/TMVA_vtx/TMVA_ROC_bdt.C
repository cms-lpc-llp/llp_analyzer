//C++ INCLUDES
#include <sys/stat.h>
#include <map>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <assert.h>
//ROOT INCLUDES
#include <TH1F.h>
#include <TH2D.h>
#include "TRandom3.h"
#include "TMVA/Reader.h"
#include <TTree.h>
#include <TPad.h>
#include <TGraphErrors.h>




const int MaxEvent = 1000000;
const int nEffBin = 200;


double SigBDT[MaxEvent], BkgBDT[MaxEvent];

double BkgEff[nEffBin], SigEff[nEffBin], errBkgEff[nEffBin], errSigEff[nEffBin];


void getEff(string inputFilename, string xmlFileName)
{
	for(int i=0;i<MaxEvent;i++)
	{
	SigBDT[i]=0.0;
	BkgBDT[i]=0.0;
	}
	for(int i=0;i<nEffBin;i++)
	{
	BkgEff[i]=0.0;
	SigEff[i]=0.0;
	errBkgEff[i]=0.0;
	errSigEff[i]=0.0;
	}
	
	float ptasym = 0.;
        float ptbal = 0.;
        float logsumpt2 = 0.;
        float pull_conv = 0.;
        float nConv = 0.;

        TMVA::Reader *vtxmvareader = 0;
        vtxmvareader = new TMVA::Reader( "!Color:Silent" );
        vtxmvareader->AddVariable("ptasym", &ptasym );
        vtxmvareader->AddVariable("ptbal", &ptbal );
        vtxmvareader->AddVariable("logsumpt2", &logsumpt2 );
        vtxmvareader->AddVariable("limPullToConv", &pull_conv );
        vtxmvareader->AddVariable("nConv", &nConv );
        vtxmvareader->BookMVA("BDT",xmlFileName.c_str());

        double MinBDT=9999.9;
	double MaxBDT=-9999.9;

	TFile *file_in = new TFile(inputFilename.c_str(),"READ");
   	
	TTree *tree_S_in = (TTree*)file_in->Get("TreeS");
	TTree *tree_S = tree_S_in->CopyTree("event%2==1");
   	tree_S->SetBranchAddress( "ptasym",         &ptasym);
   	tree_S->SetBranchAddress( "ptbal",          &ptbal);
   	tree_S->SetBranchAddress( "logsumpt2",      &logsumpt2);
   	tree_S->SetBranchAddress( "limPullToConv",  &pull_conv);
   	tree_S->SetBranchAddress( "nConv",          &nConv);
	int NEntries_S = tree_S->GetEntries();
	for(int i=0;i<NEntries_S;i++)
	{
	tree_S->GetEvent(i);
        SigBDT[i] = vtxmvareader->EvaluateMVA("BDT");
	if(SigBDT[i]<MinBDT) MinBDT = SigBDT[i];
	if(SigBDT[i]>MaxBDT) MaxBDT = SigBDT[i];
	}
	
	TTree *tree_B = (TTree*)file_in->Get("TreeB");
   	tree_B->SetBranchAddress( "ptasym",         &ptasym);
   	tree_B->SetBranchAddress( "ptbal",          &ptbal);
   	tree_B->SetBranchAddress( "logsumpt2",      &logsumpt2);
   	tree_B->SetBranchAddress( "limPullToConv",  &pull_conv);
   	tree_B->SetBranchAddress( "nConv",          &nConv);
	int NEntries_B = tree_B->GetEntries();
	for(int i=0;i<NEntries_B;i++)
	{
	tree_B->GetEvent(i);
        BkgBDT[i] = vtxmvareader->EvaluateMVA("BDT");
	if(BkgBDT[i]<MinBDT) MinBDT = BkgBDT[i];
	if(BkgBDT[i]>MaxBDT) MaxBDT = BkgBDT[i];
	}
	
	for(int i=0;i<nEffBin;i++)
	{
		double bdtCut = MinBDT + (MaxBDT-MinBDT)*i/(nEffBin*1.0);
		int SigPass = 0;
		int BkgPass = 0;
		for(int ip=0;ip<NEntries_S;ip++)
		{
		  if(SigBDT[ip]>bdtCut) SigPass++;
		}
		for(int ip=0;ip<NEntries_B;ip++)
		{
		  if(BkgBDT[ip]>bdtCut) BkgPass++;
		}
		BkgEff[i] = (BkgPass*1.0)/(NEntries_B*1.0);
		SigEff[i] = (SigPass*1.0)/(NEntries_S*1.0);
		//if(BkgPass>0 && NEntries_B) errBkgEff[i] = BkgEff[i]*sqrt(1.0/BkgPass + 1.0/NEntries_B);
	        //if(SigPass>0 && NEntries_S) errSigEff[i] = SigEff[i]*sqrt(1.0/SigPass + 1.0/NEntries_S);
		//cout<<"i: "<<i<<"  bdtCut: "<<bdtCut<<"   BkgEff: "<<BkgEff[i]<<" +- "<<errBkgEff[i]<<"   SigEff: "<<SigEff[i]<<" +- "<<errSigEff[i]<<endl;
	}
};

void TMVA_ROC_bdt()
{
	TCanvas * myC = new TCanvas("c1","c1",100,100,700,700);	
	myC->cd();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);

	getEff("HggRazorUpgradeTiming_PU140_NoTiming_vtx.root", "weights_PU140_NoTiming/TMVAClassification_BDT.weights.xml");
	TGraph *gr_1 = new TGraph(nEffBin,BkgEff,SigEff);
	gr_1->SetLineWidth(2);
	gr_1->GetXaxis()->SetTitle("non-gen-matched vtx Eff");
	gr_1->GetYaxis()->SetTitle("gen-matched vtx Eff");
	gr_1->GetXaxis()->SetRangeUser(0.0,1.05);
	gr_1->GetYaxis()->SetRangeUser(0.0,1.05);
	gr_1->SetTitle("");
	gr_1->GetYaxis()->SetTitleOffset(1.0);
	gr_1->GetYaxis()->SetTitleSize(0.06);
	gr_1->GetYaxis()->SetLabelSize(0.042);
	gr_1->GetYaxis()->SetLabelOffset(0.0);
	gr_1->GetXaxis()->SetTitleOffset(1.0);
	gr_1->GetXaxis()->SetTitleSize(0.06);
	gr_1->GetXaxis()->SetLabelSize(0.042);
	gr_1->GetXaxis()->SetLabelOffset(0.0);
	gr_1->Draw("AL");

	getEff("HggRazorUpgradeTiming_PU140_Timing_vtx.root", "weights_PU140_Timing_vtx/TMVAClassification_BDT.weights.xml");
	TGraph *gr_2 = new TGraph(nEffBin,BkgEff,SigEff);
	gr_2->SetLineColor( 3 );
  	gr_2->SetLineWidth( 2 );
	gr_2->Draw("same");

	getEff("HggRazorUpgradeTiming_PU140_NoTiming_vtx.root", "../../data/TMVAClassification_BDTVtxId_SL_2016.xml");
	TGraph *gr_5 = new TGraph(nEffBin,BkgEff,SigEff);
	gr_5->SetLineColor( 1 );
  	gr_5->SetLineWidth( 2 );
	gr_5->Draw("same");

	for(int i=0;i<nEffBin;i++)
	{
	BkgEff[i]=i*1.0/nEffBin;
	SigEff[i]=i*1.0/nEffBin;
	}
   	TGraph *gr_6 = new TGraph(nEffBin,BkgEff,SigEff);
	gr_6->SetLineColor( 7 );
  	gr_6->SetLineWidth( 2 );
  	gr_6->SetLineStyle( 7 );
	gr_6->Draw("same");

 	TLegend *leg = new TLegend(0.5, 0.25, 0.8, 0.4, NULL,"brNDC");
  	leg->SetBorderSize(0);
  	leg->SetTextSize(0.04);
  	leg->SetLineColor(1);
  	leg->SetLineStyle(1);
  	leg->SetLineWidth(1);
  	leg->SetFillColor(0);
  	leg->SetFillStyle(1001);
  	leg->AddEntry(gr_5, "flashgg default" ,"l");
  	leg->AddEntry(gr_1, "no timing" ,"l");
  	leg->AddEntry(gr_2, "track/vtx timing only" ,"l");
  	leg->Draw();



	myC->SaveAs("plots/ROC_vtx_PU140.pdf");
	myC->SaveAs("plots/ROC_vtx_PU140.png");
}
