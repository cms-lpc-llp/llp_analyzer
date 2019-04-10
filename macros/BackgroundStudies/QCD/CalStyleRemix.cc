#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TStyle.h>
#include "CalStyleRemix.hh"

void CalStyleRemix() {
  const char* author   = "$Author: jlawhorn $$";
  printf(" Caltech root style REMIX(%s).\n",author);
  printf("\n");
  printf(" Use: MakeCanvas(name,title)\n");
  printf("      InitSubPad(pad,nPad)\n");
  printf("      InitGraph(graph,xTitle,yTitle,color)\n");
  printf("      InitHist(hist,xTitle,yTitle,color)\n");
  printf("\n");
  SetStyle();
}

TCanvas* MakeCanvas(const char* name, const char *title, int dX, int dY)
{
  // Start with a canvas
  TCanvas *canvas = new TCanvas(name,title,0,0,dX,dY);
  // General overall stuff
  canvas->SetFillColor      (0);
  canvas->SetBorderMode     (0);
  canvas->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canvas->SetLeftMargin     (0.18);
  canvas->SetRightMargin    (0.05);
  canvas->SetTopMargin      (0.08);
  canvas->SetBottomMargin   (0.15);
  // Setup a frame which makes sense
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);

  return canvas;
}

void InitSubPad(TPad* pad, int i)
{
  //printf("Pad: %p, index: %d\n",pad,i);

  pad->cd(i);
  TPad *tmpPad = (TPad*) pad->GetPad(i);
  tmpPad->SetLeftMargin  (0.18);
  tmpPad->SetTopMargin   (0.05);
  tmpPad->SetRightMargin (0.07);
  tmpPad->SetBottomMargin(0.15);
  return;
}

void InitHist(TH1 *hist, const char *xtit, const char *ytit, Int_t color)
{
  hist->SetXTitle(xtit);
  hist->SetYTitle(ytit);
  hist->SetLineColor(color);
  hist->SetFillColor(color);
  hist->SetFillStyle(1001);
  hist->SetTitleSize  (0.055,"Y");
  hist->SetTitleOffset(1.600,"Y");
  hist->SetLabelOffset(0.014,"Y");
  hist->SetLabelSize  (0.050,"Y");
  hist->SetLabelFont  (42   ,"Y");
  hist->SetTitleSize  (0.055,"X");
  hist->SetTitleOffset(1.300,"X");
  hist->SetLabelOffset(0.014,"X");
  hist->SetLabelSize  (0.050,"X");
  hist->SetLabelFont  (42   ,"X");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize (0.6);
  // Strangely enough this cannot be set anywhere else??
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->SetTitle("");  
  return;
}

void InitData(TH1 *hist, const char *xtit, const char *ytit, Int_t color)
{
  hist->GetXaxis()->SetTitle(xtit);
  hist->GetYaxis()->SetTitle(ytit);
  hist->SetLineColor(color);
  hist->SetTitleSize  (0.055,"Y");
  hist->SetTitleOffset(1.600,"Y");
  hist->SetLabelOffset(0.014,"Y");
  hist->SetLabelSize  (0.050,"Y");
  hist->SetLabelFont  (42   ,"Y");
  hist->SetTitleSize  (0.055,"X");
  hist->SetTitleOffset(1.300,"X");
  hist->SetLabelOffset(0.014,"X");
  hist->SetLabelSize  (0.050,"X");
  hist->SetLabelFont  (42   ,"X");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(color);
  //hist->SetMarkerSize (0.6);
  // Strangely enough this cannot be set anywhere else??
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->SetTitle("");  
  return;
}

void InitGraph(TGraph *graph, const char *xtit, const char *ytit, Int_t color)
{
  graph->GetXaxis()->SetTitle(xtit);
  graph->GetYaxis()->SetTitle(ytit);
  graph->SetLineColor(color);
  //graph->SetTitleSize  (0.055,"Y");
  //graph->SetTitleOffset(1.600,"Y");
  //graph->SetLabelOffset(0.014,"Y");
  //graph->SetLabelSize  (0.050,"Y");
  //graph->SetLabelFont  (42   ,"Y");
  //graph->SetTitleSize  (0.055,"X");
  //graph->SetTitleOffset(1.300,"X");
  //graph->SetLabelOffset(0.014,"X");
  //graph->SetLabelSize  (0.050,"X");
  //graph->SetLabelFont  (42   ,"X");
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(color);
  graph->SetMarkerSize (0.6);
  // Strangely enough this cannot be set anywhere else??
  graph->GetYaxis()->SetTitleFont(42);
  graph->GetXaxis()->SetTitleFont(42);
  graph->SetTitle("");  
  return;
}

void SetStyle()
{
  TStyle *CalStyle = new TStyle("Caltech-Style","The Even Perfecter Style for Plots ;-)");
  gStyle = CalStyle;

  CalStyle->SetPalette(56);

  // Canvas
  CalStyle->SetCanvasColor     (0);
  CalStyle->SetCanvasBorderSize(10);
  CalStyle->SetCanvasBorderMode(0);
  CalStyle->SetCanvasDefH      (700);
  CalStyle->SetCanvasDefW      (700);
  CalStyle->SetCanvasDefX      (100);
  CalStyle->SetCanvasDefY      (100);

  // Pads
  CalStyle->SetPadColor       (0);
  CalStyle->SetPadBorderSize  (10);
  CalStyle->SetPadBorderMode  (0);
  CalStyle->SetPadBottomMargin(0.13);
  CalStyle->SetPadTopMargin   (0.08);
  CalStyle->SetPadLeftMargin  (0.15);
  CalStyle->SetPadRightMargin (0.05);
  CalStyle->SetPadGridX       (0);
  CalStyle->SetPadGridY       (0);
  CalStyle->SetPadTickX       (0);
  CalStyle->SetPadTickY       (0);

  // Frames
  CalStyle->SetFrameFillStyle ( 0);
  CalStyle->SetFrameFillColor ( 0);
  CalStyle->SetFrameLineColor ( 1);
  CalStyle->SetFrameLineStyle ( 0);
  CalStyle->SetFrameLineWidth ( 1);
  CalStyle->SetFrameBorderSize(10);
  CalStyle->SetFrameBorderMode( 0);

  // Histograms
  CalStyle->SetHistFillColor(2);
  CalStyle->SetHistFillStyle(0);
  CalStyle->SetHistLineColor(1);
  CalStyle->SetHistLineStyle(0);
  CalStyle->SetHistLineWidth(2);
  CalStyle->SetNdivisions(505);

  // Functions
  CalStyle->SetFuncColor(1);
  CalStyle->SetFuncStyle(0);
  CalStyle->SetFuncWidth(2);

  // Various
  CalStyle->SetMarkerStyle(20);
  CalStyle->SetMarkerColor(kBlack);
  CalStyle->SetMarkerSize (1.2);

  CalStyle->SetTitleBorderSize(0);
  CalStyle->SetTitleFillColor (0);
  CalStyle->SetTitleX         (0.2);

  CalStyle->SetTitleSize  (0.055,"X");
  CalStyle->SetTitleOffset(1.200,"X");
  CalStyle->SetLabelOffset(0.005,"X");
  CalStyle->SetLabelSize  (0.050,"X");
  CalStyle->SetLabelFont  (42   ,"X");

  CalStyle->SetStripDecimals(kFALSE);

  CalStyle->SetTitleSize  (0.055,"Y");
  CalStyle->SetTitleOffset(1.600,"Y");
  CalStyle->SetLabelOffset(0.010,"Y");
  CalStyle->SetLabelSize  (0.050,"Y");
  CalStyle->SetLabelFont  (42   ,"Y");

  CalStyle->SetTextSize   (0.055);
  CalStyle->SetTextFont   (42);
  CalStyle->SetPaintTextFormat("4.2f");

  CalStyle->SetStatFont   (42);
  CalStyle->SetTitleFont  (42);
  CalStyle->SetTitleFont  (42,"X");
  CalStyle->SetTitleFont  (42,"Y");

  CalStyle->SetOptStat    (0);
  return;
}
