
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



void vtxEff(string inputDir = "../../maxBDTpv/", string method = "max BDT (default)")
{

   int NvtxTotal[6];
   int NvtxSelected[6];
   for(int i=0;i<6;i++)
	{
	NvtxTotal[i] = 0;
	NvtxSelected[i] = 0;
	}
   TFile *file_PU0_NoTiming = new TFile((inputDir+"HggRazorUpgradeTiming_PU0_NoTiming.root").c_str(),"READ");
   TTree *tree_PU0_NoTiming = (TTree*)file_PU0_NoTiming->Get("HggRazor");
   NvtxTotal[0] = tree_PU0_NoTiming->GetEntries();
   NvtxSelected[0] = tree_PU0_NoTiming->GetEntries("vtxdZ<1.0");

   TFile *file_PU0_Timing = new TFile((inputDir+"HggRazorUpgradeTiming_PU0_Timing.root").c_str(),"READ");
   TTree *tree_PU0_Timing = (TTree*)file_PU0_Timing->Get("HggRazor");
   NvtxTotal[1] = tree_PU0_Timing->GetEntries();
   NvtxSelected[1] = tree_PU0_Timing->GetEntries("vtxdZ<1.0");

   TFile *file_PU140_NoTiming = new TFile((inputDir+"HggRazorUpgradeTiming_PU140_NoTiming.root").c_str(),"READ");
   TTree *tree_PU140_NoTiming = (TTree*)file_PU140_NoTiming->Get("HggRazor");
   NvtxTotal[2] = tree_PU140_NoTiming->GetEntries();
   NvtxSelected[2] = tree_PU140_NoTiming->GetEntries("vtxdZ<1.0");

   TFile *file_PU140_Timing = new TFile((inputDir+"HggRazorUpgradeTiming_PU140_Timing.root").c_str(),"READ");
   TTree *tree_PU140_Timing = (TTree*)file_PU140_Timing->Get("HggRazor");
   NvtxTotal[3] = tree_PU140_Timing->GetEntries();
   NvtxSelected[3] = tree_PU140_Timing->GetEntries("vtxdZ<1.0");

   TFile *file_PU200_NoTiming = new TFile((inputDir+"HggRazorUpgradeTiming_PU200_NoTiming.root").c_str(),"READ");
   TTree *tree_PU200_NoTiming = (TTree*)file_PU200_NoTiming->Get("HggRazor");
   NvtxTotal[4] = tree_PU200_NoTiming->GetEntries();
   NvtxSelected[4] = tree_PU200_NoTiming->GetEntries("vtxdZ<1.0");

   TFile *file_PU200_Timing = new TFile((inputDir+"HggRazorUpgradeTiming_PU200_Timing.root").c_str(),"READ");
   TTree *tree_PU200_Timing = (TTree*)file_PU200_Timing->Get("HggRazor");
   NvtxTotal[5] = tree_PU200_Timing->GetEntries();
   NvtxSelected[5] = tree_PU200_Timing->GetEntries("vtxdZ<1.0");
 
   FILE* m_outfile = fopen("vtx_eff_table.tex", "a");

   fprintf(m_outfile, "%s", method.c_str());

   for(int i=0;i<6;i++)
	{
	fprintf(m_outfile, "& %4.1f\\%% ", (100.0*NvtxSelected[i])/(NvtxTotal[i]*1.0));
	}

   fprintf(m_outfile, " \\\\ \\hline \n");

}
