#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <TPad.h>

void     CalStyleRemix();
TCanvas* MakeCanvas   (const char* name, const char *title, int dX = 500, int dY = 500);
void     InitSubPad   (TPad* pad, int i);
void     InitHist     (TH1 *hist, const char *xtit, const char *ytit  = "Number of Entries",
		       Int_t color = kBlack);
void     InitData     (TH1 *hist, const char *xtit, const char *ytit  = "Number of Entries",
		       Int_t color = kBlack);
void     InitGraph    (TGraph *graph, const char *xtit, const char *ytit  = "Number of Entries",
		       Int_t color = kBlack);
void     SetStyle     ();
