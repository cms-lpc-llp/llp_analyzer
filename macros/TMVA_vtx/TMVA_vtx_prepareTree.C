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

void TMVA_vtx_prepareTree(string inputFilename = "../../HggRazorUpgradeTiming_PU0_NoTiming.root", string outputFilename = "HggRazorUpgradeTiming_PU0_NoTiming_vtx.root")
{
	unsigned int event;
	unsigned int event_out;
	
	float chi2_pho_vtx = 0.0;
	float logchi2_pho_vtx = 0.0;
	Int_t rank_chi2_pho_vtx = 0;
	float ptasym = 0.;
  	float ptbal = 0.;
  	float logsumpt2 = 0.;
  	float pull_conv = 0.;
  	float nConv = 0.;
  	float weight = 1.0;

	float vtxdX = 0.0;
	float vtxdY = 0.0;
	float vtxdZ = 0.0;
	float vtxdT = 0.0;

	float vtxPt = 0.0;
	float vtxSumPt = 0.0;
	float diphoPt = 0.0;	

	Int_t pvNtrack = 0;

	float vtxdX_all[500];
	float vtxdY_all[500];
	float vtxdZ_all[500];
	float vtxdT_all[500];
	
	float vtxPt_all[500];
	float vtxSumPt_all[500];
	Int_t pvNtrack_all[500];
	
	Int_t nPVAll = 0;
	Int_t isMatchPv[500];
  	Int_t isMaxbdtPv[500];
   	float chi2_pho_vtx_all[500];// = {0.};
   	Int_t rank_chi2_pho_vtx_all[500];// = {0.};
   	float ptasym_all[500];// = {0.};
  	float ptbal_all[500];// = {0.};
  	float logsumpt2_all[500];// = {0.};
  	float pull_conv_all[500];// = {0.};
  	float nConv_all[500];// = {0.};

	for(int i=0;i<500;i++)
  	{
        isMatchPv[i]=0;
        isMaxbdtPv[i]=0;
        ptasym_all[i]=0.0;
        ptbal_all[i]=0.0;
        logsumpt2_all[i]=0.0;
        pull_conv_all[i]=0.0;
        nConv_all[i]=0.0;

	vtxdX_all[i] = 0.0;
	vtxdY_all[i] = 0.0;
	vtxdZ_all[i] = 0.0;
	vtxdT_all[i] = 0.0;

	vtxPt_all[i] = 0.0;
	vtxSumPt_all[i] = 0.0;
	pvNtrack_all[i] = 0;
	
  	}


   TFile *file_in = new TFile(inputFilename.c_str(),"READ");
   TTree *tree_in = (TTree*)file_in->Get("HggRazor");
    
   tree_in->SetBranchAddress( "event", 		&event);
   tree_in->SetBranchAddress( "nPVAll", 	&nPVAll);
   tree_in->SetBranchAddress( "isMatchPv",   	isMatchPv);
   if(inputFilename.find("NoTiming")!=std::string::npos) tree_in->SetBranchAddress( "chi2_min_pho_vtx", chi2_pho_vtx_all);
   else tree_in->SetBranchAddress( "chi2_pho_vtx", chi2_pho_vtx_all);
   tree_in->SetBranchAddress( "ptasym",   	ptasym_all);
   tree_in->SetBranchAddress( "ptbal",   	ptbal_all);
   tree_in->SetBranchAddress( "logsumpt2",   	logsumpt2_all);
   tree_in->SetBranchAddress( "limPullToConv",  pull_conv_all);
   tree_in->SetBranchAddress( "nConv",   	nConv_all);
   
   tree_in->SetBranchAddress( "vtxdX_all",   	vtxdX_all);
   tree_in->SetBranchAddress( "vtxdY_all",   	vtxdY_all);
   tree_in->SetBranchAddress( "vtxdZ_all",   	vtxdZ_all);
   tree_in->SetBranchAddress( "vtxdT_all",   	vtxdT_all);
   tree_in->SetBranchAddress( "vertexpt",   	vtxPt_all);
   tree_in->SetBranchAddress( "vertexsumpt",   	vtxSumPt_all);
   tree_in->SetBranchAddress( "diphotonpt",   	&diphoPt);
   tree_in->SetBranchAddress( "pvNtrack",       pvNtrack_all);


   TFile *file_out = new TFile(outputFilename.c_str(),"RECREATE");
   TTree *tree_out_S = new TTree("TreeS", "Tree_Signal");
   TTree *tree_out_B = new TTree("TreeB", "Tree_Background");

   tree_out_S->Branch( "event",		&event_out, "event/i");
   tree_out_S->Branch( "chi2_pho_vtx",	&chi2_pho_vtx, "chi2_pho_vtx/F");
   tree_out_S->Branch( "rank_chi2_pho_vtx",	&rank_chi2_pho_vtx, "rank_chi2_pho_vtx/I");
   tree_out_S->Branch( "logchi2_pho_vtx",	&logchi2_pho_vtx, "logchi2_pho_vtx/F");
   tree_out_S->Branch( "ptasym",	&ptasym, "ptasym/F");
   tree_out_S->Branch( "ptbal", 	&ptbal, "ptbal/F");
   tree_out_S->Branch( "logsumpt2", 	&logsumpt2, "logsumpt2/F");
   tree_out_S->Branch( "limPullToConv", &pull_conv, "limPullToConv/F");
   tree_out_S->Branch( "nConv", 	&nConv, "nConv/F");
   tree_out_S->Branch( "vtxdX", 	&vtxdX, "vtxdX/F");
   tree_out_S->Branch( "vtxdY", 	&vtxdY, "vtxdY/F");
   tree_out_S->Branch( "vtxdZ", 	&vtxdZ, "vtxdZ/F");
   tree_out_S->Branch( "vtxdT", 	&vtxdT, "vtxdT/F");
   tree_out_S->Branch( "vtxPt", 	&vtxPt, "vtxPt/F");
   tree_out_S->Branch( "vtxSumPt", 	&vtxSumPt, "vtxSumPt/F");
   tree_out_S->Branch( "diphoPt", 	&diphoPt, "diphoPt/F");
   tree_out_S->Branch( "pvNtrack", 	&pvNtrack, "pvNtrack/I");
   //tree_out_S->Branch( "weight", 	&weight, "weight/F");

   tree_out_B->Branch( "event",		&event_out, "event/i");
   tree_out_B->Branch( "chi2_pho_vtx",	&chi2_pho_vtx, "chi2_pho_vtx/F");
   tree_out_B->Branch( "rank_chi2_pho_vtx",	&rank_chi2_pho_vtx, "rank_chi2_pho_vtx/I");
   tree_out_B->Branch( "logchi2_pho_vtx",	&logchi2_pho_vtx, "logchi2_pho_vtx/F");
   tree_out_B->Branch( "ptasym",	&ptasym, "ptasym/F");
   tree_out_B->Branch( "ptbal", 	&ptbal, "ptbal/F");
   tree_out_B->Branch( "logsumpt2", 	&logsumpt2, "logsumpt2/F");
   tree_out_B->Branch( "limPullToConv", &pull_conv, "limPullToConv/F");
   tree_out_B->Branch( "nConv", 	&nConv, "nConv/F");
   tree_out_B->Branch( "weight", 	&weight, "weight/F");
   tree_out_B->Branch( "vtxdX", 	&vtxdX, "vtxdX/F");
   tree_out_B->Branch( "vtxdY", 	&vtxdY, "vtxdY/F");
   tree_out_B->Branch( "vtxdZ", 	&vtxdZ, "vtxdZ/F");
   tree_out_B->Branch( "vtxdT", 	&vtxdT, "vtxdT/F");
   tree_out_B->Branch( "vtxPt", 	&vtxPt, "vtxPt/F");
   tree_out_B->Branch( "vtxSumPt", 	&vtxSumPt, "vtxSumPt/F");
   tree_out_B->Branch( "diphoPt", 	&diphoPt, "diphoPt/F");
   tree_out_B->Branch( "pvNtrack", 	&pvNtrack, "pvNtrack/I");
 

   int NEntries =  tree_in->GetEntries();
   for(int i=0;i<NEntries;i++)
   {
	tree_in->GetEntry(i);
	event_out = event;
	for(int jj=0;jj<nPVAll;jj++)
	{
		rank_chi2_pho_vtx_all[jj] = 1;
		for(int jjj=0;(jjj<nPVAll && jjj!=jj);jjj++)
		{
			if(chi2_pho_vtx_all[jjj]<chi2_pho_vtx_all[jj])
			{
				rank_chi2_pho_vtx_all[jj] ++;
			}
		}
	}
	
	for(int j=0;j<nPVAll;j++)
	{	
		chi2_pho_vtx = chi2_pho_vtx_all[j];
		rank_chi2_pho_vtx = rank_chi2_pho_vtx_all[j];
		if(chi2_pho_vtx>0.0) logchi2_pho_vtx = log(chi2_pho_vtx);
		ptasym = ptasym_all[j];
		ptbal = ptbal_all[j];
		logsumpt2 = logsumpt2_all[j];
		pull_conv = pull_conv_all[j];
		nConv = nConv_all[j];
		vtxdX = vtxdX_all[j];
		vtxdY = vtxdY_all[j];
		vtxdZ = vtxdZ_all[j];
		vtxdT = vtxdT_all[j];

		vtxPt = vtxPt_all[j];
		vtxSumPt = vtxSumPt_all[j];

		pvNtrack = pvNtrack_all[j];

		if(!isinf(logsumpt2))
		{
		if(isMatchPv[j]==0)
		{
			tree_out_B->Fill();
		}		
		else
		{
			tree_out_S->Fill();
		}
		}
	}
   }

   tree_out_B->Write();
   tree_out_S->Write();	

}
