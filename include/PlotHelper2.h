#ifndef PlotHelper6612
#define PlotHelper6612

// Plot helper, 6612 version
// Intended to help organizing plots and to avoid cluttering
// The main idea is to write a plot both to a pdf and to a TFile at the same time,
//    so that we can check quickly using pdf and come back to TFile later if needed
//    the helper functions here will be mostly on *.ps files
// Added feature compared to 6152: being able to insert a table of contents (and it's clickable!)
//    automatic home button
//    automatic page number
// Developer: Yi Chen, 6612

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
using namespace std;

#include "TCanvas.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TVirtualPad.h"
#include "TVirtualPS.h"
#include "TStyle.h"

class PsFileHelper
{
private:
   string FileName;
   string Option;
   bool Status;   // false - file not opened; true - file good and ongoing
private:
   bool AutomaticHomeButton;
   double HomeButtonX;
   double HomeButtonY;
   double HomeButtonSize;
   string HomeButtonDestination;
   bool PrintPageNumber;
   int NextPageNumber;
public:
   PsFileHelper();
   PsFileHelper(string filename);
   PsFileHelper(string filename, string option);
   ~PsFileHelper();
   void Open(string filename);
   void Close(bool Convert = true);
   string GetFileName();
   void SetOption(string option);
   string GetOption();
   template <class PlotType> void AddPlot(PlotType *Histogram, string PlotOption = "",
      bool LogY = false, bool LogZ = false, bool Grid = false, bool LogX = false);
   template <class PlotType> void Add4PanelPlot(PlotType *Histogram1, PlotType *Histogram2,
      PlotType *Histogram3, PlotType *Histogram4, string PlotOption = "",
      bool LogY = false, bool LogZ = false, bool Grid = false, bool LogX = false);
   template <class PlotType> void AddPlot(PlotType &Histogram, string PlotOption = "",
      bool LogY = false, bool LogZ = false, bool Grid = false, bool LogX = false);
   template <class PlotType> void Add4PanelPlot(PlotType &Histogram1, PlotType &Histogram2,
      PlotType &Histogram3, PlotType &Histogram4, string PlotOption = "",
      bool LogY = false, bool LogZ = false, bool Grid = false, bool LogX = false);
   template <class PlotType> void AddNormalizedPlot(PlotType *Histogram, string PlotOption = "",
      bool LogY = false, bool LogZ = false, bool Grid = false);
   template <class PlotType> void AddNormalizedPlot(PlotType &Histogram, string PlotOption = "",
      bool LogY = false, bool LogZ = false, bool Grid = false);
   template <class PlotType> void AddPlotWithText(PlotType *Histogram, string Text,
      string PlotOption = "", double X = 0.1, double Y = 0.9, double TextSize = 0.03);
   template <class PlotType> void AddPlotWithText(PlotType &Histogram, string Text,
      string PlotOption = "", double X = 0.1, double Y = 0.9, double TextSize = 0.03);
   void AddHistogramFromFile(TFile &File, string HistogramName,
      string PlotOption = "", bool LogY = false, bool LogZ = false, bool Grid = false, bool LogX = false);
   void Add4PanelHistogramFromFile(TFile &File, string HistogramName1, string HistogramName2,
      string HistogramName3, string HistogramName4,
      string PlotOption = "", bool LogY = false, bool LogZ = false, bool Grid = false, bool LogX = false);
   void AddGraphFromFile(TFile &File, string GraphName,
      string PlotOption = "", bool LogY = false, bool LogZ = false, bool Grid = false, bool LogX = false);
   void AddCanvas(TCanvas *Canvas);
   void AddCanvas(TCanvas &Canvas);
   void AddCanvasWithText(TCanvas *Canvas, string Text, double X = 0.1, double Y = 0.9, double TextSize = 0.03);
   void AddCanvasWithText(TCanvas &Canvas, string Text, double X = 0.1, double Y = 0.9, double TextSize = 0.03);
   void AddTextPage(string Text, double X = 0.15, double Y = 0.5, double TextSize = 0.05);
   void AddTextPage(vector<string> Text, double X = 0.1, double Y = 0.9, double TextSize = 0.04);
   void AddTimeStampPage();
   void AddTableOfContentPage(vector<string> Items, vector<string> Destinations,
      double X = 0.11, double Y = 0.8, double TextSize = 0.03);
public:
   void PrintStringToPSFile(string String);
   void PrintLineToPSFile(string Line);
   void InsertNamedDestination(string Name);
   void InsertBallLinkAbsolute(int X, int Y, int Radius, string Destination);
   void InsertHomeButtonAbsolute(double X, double Y, double Size, string Destination);
   void SetAutomaticHomeButton(bool newvalue = true, string Destination = "HomePage",
      double X = 50, double Y = 50, double Size = 75);
   void SetPageNumber(bool printnumber);
   void InsertPageNumber(TCanvas &Canvas, int PageNumber = -1);
};

PsFileHelper::PsFileHelper()
{
   Status = false;
   FileName = "";
   Option = "Landscape";

   AutomaticHomeButton = false;
   HomeButtonX = 50;
   HomeButtonY = 50;
   HomeButtonSize = 75;
   HomeButtonDestination = "";
   PrintPageNumber = true;
   NextPageNumber = 1;
}

PsFileHelper::PsFileHelper(string filename)
{
   Option = "Landscape";

   AutomaticHomeButton = false;
   HomeButtonX = 50;
   HomeButtonY = 50;
   HomeButtonSize = 75;
   HomeButtonDestination = "";
   PrintPageNumber = true;
   NextPageNumber = 1;
   
   Open(filename);
}

PsFileHelper::PsFileHelper(string filename, string option)
{
   Option = option;
   
   AutomaticHomeButton = false;
   HomeButtonX = 50;
   HomeButtonY = 50;
   HomeButtonSize = 75;
   HomeButtonDestination = "";
   PrintPageNumber = true;
   NextPageNumber = 1;

   Open(filename);
}

PsFileHelper::~PsFileHelper()
{
   if(Status == true)
      Close(true);
}

void PsFileHelper::Open(string filename)
{
   FileName = filename;

   Status = true;

   TCanvas canvas;
   canvas.Print((FileName + "[").c_str(), Option.c_str());

   TVirtualPS *PSHandle = (TVirtualPS *)gROOT->GetListOfSpecials()->FindObject(FileName.c_str());
   if(PSHandle != NULL)
   {
      PSHandle->PrintStr("\n");
      PSHandle->PrintStr("%% pdfmark prolog: tell ps reader to ignore pdfmark if it's not supported\n");
      PSHandle->PrintStr("/pdfmark where {pop} {userdict /pdfmark /cleartomark load put} ifelse\n");
      PSHandle->PrintStr("\n");
   }

   NextPageNumber = 1;
}

void PsFileHelper::Close(bool Convert)
{
   TCanvas canvas;
   canvas.Print((FileName + "]").c_str(), Option.c_str());

   if(Convert == true)
   {
      gROOT->ProcessLine((".! ps2pdf14 " + FileName).c_str());
      gROOT->ProcessLine((".! rm " + FileName).c_str());
   }

   Status = false;
   FileName = "";
   NextPageNumber = -1;
}

string PsFileHelper::GetFileName()
{
   return FileName;
}

void PsFileHelper::SetOption(string option)
{
   Option = option;
}

string PsFileHelper::GetOption()
{
   return Option;
}

template <class PlotType> void PsFileHelper::AddPlot(PlotType *Histogram, string PlotOption, bool LogY, bool LogZ, bool Grid, bool LogX)
{
   if(Histogram == NULL)
      return;

   TCanvas canvas;

   Histogram->Draw(PlotOption.c_str());

   if(LogX == true)
      canvas.SetLogx();
   if(LogY == true)
      canvas.SetLogy();
   if(LogZ == true)
      canvas.SetLogz();

   if(Grid == true)
   {
      canvas.SetGridx();
      canvas.SetGridy();
   }

   AddCanvas(canvas);
}

template <class PlotType> void PsFileHelper::Add4PanelPlot(PlotType *Histogram1, PlotType *Histogram2,
   PlotType *Histogram3, PlotType *Histogram4, string PlotOption, bool LogY, bool LogZ, bool Grid, bool LogX)
{
   TCanvas canvas;

   canvas.Divide(2, 2);

   PlotType *Histograms[4] = {Histogram1, Histogram2, Histogram3, Histogram4};

   for(int i = 0; i < 4; i++)
   {
      if(Histograms[i] == NULL)
         continue;

      TVirtualPad *Pad = canvas.cd(i + 1);

      Histograms[i]->Draw(PlotOption.c_str());

      if(LogX == true)
         Pad->SetLogx();
      if(LogY == true)
         Pad->SetLogy();
      if(LogZ == true)
         Pad->SetLogz();

      if(Grid == true)
      {
         Pad->SetGridx();
         Pad->SetGridy();
      }
   }

   AddCanvas(canvas);
}

template <class PlotType> void PsFileHelper::AddPlot(PlotType &Histogram, string PlotOption, bool LogY, bool LogZ, bool Grid, bool LogX)
{
   AddPlot(&Histogram, PlotOption, LogY, LogZ, Grid, LogX);
}

template <class PlotType> void Add4PanelPlot(PlotType &Histogram1, PlotType &Histogram2,
   PlotType &Histogram3, PlotType &Histogram4, string PlotOption, bool LogY, bool LogZ, bool Grid, bool LogX)
{
   Add4PanelPlot(&Histogram1, &Histogram2, &Histogram3, &Histogram4, PlotOption, LogY, LogZ, Grid, LogX);
}

template <class PlotType> void PsFileHelper::AddNormalizedPlot(PlotType *Histogram, string PlotOption,
   bool LogY, bool LogZ, bool Grid)
{
   if(Histogram == NULL)
      return;

   TCanvas canvas;

   Histogram->DrawNormalized(PlotOption.c_str());

   if(Histogram->GetMaximum() > 1e-10)
   {
      if(LogY == true)
         canvas.SetLogy();
      if(LogZ == true)
         canvas.SetLogz();
   }
   if(Grid == true)
   {
      canvas.SetGridx();
      canvas.SetGridy();
   }

   AddCanvas(canvas);
}

template <class PlotType> void PsFileHelper::AddNormalizedPlot(PlotType &Histogram, string PlotOption,
   bool LogY, bool LogZ, bool Grid)
{
   AddNormalizedPlot(&Histogram, PlotOption, LogY, LogZ, Grid);
}

template <class PlotType> void PsFileHelper::AddPlotWithText(PlotType *Histogram, string Text,
   string PlotOption, double X, double Y, double TextSize)
{
   TCanvas canvas;

   Histogram->Draw(PlotOption.c_str());

   AddCanvasWithText(canvas, Text, X, Y, TextSize);
}

template <class PlotType> void PsFileHelper::AddPlotWithText(PlotType &Histogram, string Text,
   string PlotOption, double X, double Y, double TextSize)
{
   AddPlotWithText(&Histogram, Text, PlotOption, X, Y, TextSize);
}

void PsFileHelper::AddHistogramFromFile(TFile &File, string HistogramName, string PlotOption, bool LogY,
   bool LogZ, bool Grid, bool LogX)
{
   TH1 *Histogram = (TH1 *)File.Get(HistogramName.c_str());
   if(Histogram == NULL)
      return;

   AddPlot(Histogram, PlotOption, LogY, LogZ, Grid, LogX);
}

void PsFileHelper::Add4PanelHistogramFromFile(TFile &File, string HistogramName1, string HistogramName2,
   string HistogramName3, string HistogramName4, string PlotOption, bool LogY, bool LogZ, bool Grid, bool LogX)
{
   TH1 *Histograms[4];

   Histograms[0] = (TH1 *)File.Get(HistogramName1.c_str());
   Histograms[1] = (TH1 *)File.Get(HistogramName2.c_str());
   Histograms[2] = (TH1 *)File.Get(HistogramName3.c_str());
   Histograms[3] = (TH1 *)File.Get(HistogramName4.c_str());

   Add4PanelPlot(Histograms[0], Histograms[1], Histograms[2], Histograms[3], PlotOption, LogY, LogZ, Grid, LogX);
}

void PsFileHelper::AddGraphFromFile(TFile &File, string GraphName, string PlotOption, bool LogY, bool LogZ,
   bool Grid, bool LogX)
{
   TGraph *Graph = (TGraph *)File.Get(GraphName.c_str());
   if(Graph == NULL)
      return;

   AddPlot(Graph, PlotOption, LogY, LogZ, Grid, LogX);
}

void PsFileHelper::AddCanvas(TCanvas *Canvas)
{
   if(Canvas == NULL)
      return;

   if(PrintPageNumber == true)
      InsertPageNumber(*Canvas, NextPageNumber);
   Canvas->Print(FileName.c_str(), Option.c_str());
   NextPageNumber = NextPageNumber + 1;

   if(AutomaticHomeButton == true)
      InsertHomeButtonAbsolute(HomeButtonX, HomeButtonY, HomeButtonSize, HomeButtonDestination);
}

void PsFileHelper::AddCanvas(TCanvas &Canvas)
{
   AddCanvas(&Canvas);
}

void PsFileHelper::AddCanvasWithText(TCanvas *Canvas, string Text, double X, double Y, double TextSize)
{
   Canvas->cd();

   TLatex text(X, Y, Text.c_str());
   text.SetTextSize(TextSize);
   text.Draw();

   if(PrintPageNumber == true)
      InsertPageNumber(*Canvas, NextPageNumber);
   Canvas->Print(FileName.c_str(), Option.c_str());
   NextPageNumber = NextPageNumber + 1;

   if(AutomaticHomeButton == true)
      InsertHomeButtonAbsolute(HomeButtonX, HomeButtonY, HomeButtonSize, HomeButtonDestination);
}

void PsFileHelper::AddCanvasWithText(TCanvas &Canvas, string Text, double X, double Y, double TextSize)
{
   AddCanvasWithText(&Canvas, Text, X, Y, TextSize);
}

void PsFileHelper::AddTextPage(string Text, double X, double Y, double TextSize)
{
   TCanvas canvas;

   TLatex text(X, Y, Text.c_str());
   text.SetTextFont(42);
   text.SetTextSize(TextSize);
   text.Draw();

   if(PrintPageNumber == true)
      InsertPageNumber(canvas, NextPageNumber);
   canvas.Print(FileName.c_str(), Option.c_str());
   NextPageNumber = NextPageNumber + 1;

   if(AutomaticHomeButton == true)
      InsertHomeButtonAbsolute(HomeButtonX, HomeButtonY, HomeButtonSize, HomeButtonDestination);
}

void PsFileHelper::AddTextPage(vector<string> Text, double X, double Y, double TextSize)
{
   TCanvas canvas;

   vector<TLatex *> texts;

   for(int i = 0; i < (int)Text.size(); i++)
   {
      texts.push_back(new TLatex(X, Y - i * TextSize * 1.75, Text[i].c_str()));
      texts[i]->SetName(Form("TextLine%d", i));
      texts[i]->SetTextFont(42);
      texts[i]->SetTextSize(TextSize);
   }

   for(int i = 0; i < (int)Text.size(); i++)
      texts[i]->Draw();

   if(PrintPageNumber == true)
      InsertPageNumber(canvas, NextPageNumber);
   canvas.Print(FileName.c_str(), Option.c_str());
   NextPageNumber = NextPageNumber + 1;

   if(AutomaticHomeButton == true)
      InsertHomeButtonAbsolute(HomeButtonX, HomeButtonY, HomeButtonSize, HomeButtonDestination);

   for(int i = 0; i < (int)Text.size(); i++)
      delete texts[i];
   texts.clear();
}

void PsFileHelper::AddTimeStampPage()
{
   time_t CurrentTime = time(NULL);

   string str = "Generated at ";
   str = str + ctime(&CurrentTime);

   AddTextPage(str, 0.03);
}

void PsFileHelper::AddTableOfContentPage(vector<string> Items, vector<string> Destinations, double X, double Y,
   double TextSize)
{
   TCanvas canvas;

   vector<TLatex *> texts;

   texts.push_back(new TLatex(0.37, Y + 0.1, "Table of Contents"));
   texts[0]->SetName("TextLine00000");
   texts[0]->SetTextSize(0.05);

   for(int i = 0; i < (int)Items.size(); i++)
   {
      texts.push_back(new TLatex(X, Y - i * TextSize * 1.5, Items[i].c_str()));
      texts[i+1]->SetName(Form("TextLine%d", i));
      texts[i+1]->SetTextSize(TextSize);
   }

   for(int i = 0; i < (int)texts.size(); i++)
      texts[i]->Draw();

   if(PrintPageNumber == true)
      InsertPageNumber(canvas, NextPageNumber);
   canvas.Print(FileName.c_str(), Option.c_str());
   NextPageNumber = NextPageNumber + 1;

   if(AutomaticHomeButton == true)
      InsertHomeButtonAbsolute(HomeButtonX, HomeButtonY, HomeButtonSize, HomeButtonDestination);

   for(int i = 0; i < (int)Items.size() && i < (int)Destinations.size(); i++)
   {
      if(Destinations[i] == "")
         continue;

      double TextWidth = 0.8;

      double MagicNumber = 4.235;   // will need to change to something better
      double X1 = canvas.XtoPixel(X - TextSize * 0.2) * MagicNumber;
      double X2 = canvas.XtoPixel(X + TextWidth + TextSize * 0.2) * MagicNumber;
      double Y1 = canvas.YtoPixel(1 - (Y - i * TextSize * 1.5) + TextSize * 0.2) * MagicNumber;
      double Y2 = canvas.YtoPixel(1 - (Y - i * TextSize * 1.5) - TextSize * 1) * MagicNumber;

      stringstream box;
      box << "[ /Rect [" << X1 << " " << Y1 << " " << X2 << " " << Y2 << "]";

      PrintLineToPSFile("");
      PrintLineToPSFile(box.str());
      PrintLineToPSFile("/Color [0 0 1]");
      PrintLineToPSFile("/Dest /" + Destinations[i]);
      PrintLineToPSFile("/Subtype /Link");
      PrintLineToPSFile("/ANN pdfmark");
      PrintLineToPSFile("");
   }

   for(int i = 0; i < (int)texts.size(); i++)
      delete texts[i];
   texts.clear();
}

void PsFileHelper::PrintStringToPSFile(string String)
{
   TVirtualPS *PSHandle = (TVirtualPS *)gROOT->GetListOfSpecials()->FindObject(FileName.c_str());
   if(PSHandle == NULL)   // something's wrong.  file not opened yet?   anyways do nothing.
      return;

   PSHandle->PrintStr(String.c_str());
}

void PsFileHelper::PrintLineToPSFile(string Line)
{
   PrintStringToPSFile(Line + "\n");
}

void PsFileHelper::InsertNamedDestination(string Name)
{
   PrintLineToPSFile("");
   // PrintLineToPSFile("%% Defining named destination: " + Name);
   PrintLineToPSFile("[ /Dest /" + Name);
   PrintLineToPSFile("   /DEST pdfmark");
   PrintLineToPSFile("");
}

void PsFileHelper::InsertBallLinkAbsolute(int X, int Y, int Radius, string Destination)
{
   stringstream str1;
   str1 << "newpath " << X << " " << Y << " " << Radius << " 0 360 arc";
   stringstream str2;
   str2 << "[ /Rect [" << X - Radius << " " << Y - Radius
      << " " << X + Radius << " " << Y + Radius << "]";

   PrintLineToPSFile("");
   PrintLineToPSFile(str1.str());
   PrintLineToPSFile("fill");
   PrintLineToPSFile("");
   PrintLineToPSFile(str2.str());
   PrintLineToPSFile("/Color [0 0 1]");
   PrintLineToPSFile("/Dest /" + Destination);
   PrintLineToPSFile("/Subtype /Link");
   PrintLineToPSFile("/ANN pdfmark");
   PrintLineToPSFile("");
}

void PsFileHelper::InsertHomeButtonAbsolute(double X, double Y, double Size, string Destination)
{
   double H = 0.4;   // height and width of the small box in the bottom, in units of "Size"
   double W = 0.75;

   double X2 = 0.7;   // chimney left side x (lower side y is always 0.5), and width/height in units of "Size"
   double W2 = 0.125;
   double H2 = 0.4;

   double LowerLeftX = X - Size / 2;
   double LowerLeftY = Y - Size / 2;

   double PointsX[11] = {0};
   double PointsY[11] = {0};

   PointsX[0] = 0;              PointsY[0] = H;
   PointsX[1] = 0.5;            PointsY[1] = 1;

   PointsX[2] = X2;             PointsY[2] = (1 - X2) / (1 - 0.5) * (1 - H) + H;
   PointsX[3] = X2;             PointsY[3] = H2 + H;
   PointsX[4] = X2 + W2;        PointsY[4] = H2 + H;
   PointsX[5] = X2 + W2;        PointsY[5] = (1 - X2 - W2) / (1 - 0.5) * (1 - H) + H;

   PointsX[6] = 1;              PointsY[6] = H;
   PointsX[7] = 0.5 + W / 2;    PointsY[7] = H;
   PointsX[8] = 0.5 + W / 2;    PointsY[8] = 0;
   PointsX[9] = 0.5 - W / 2;    PointsY[9] = 0;
   PointsX[10] = 0.5 - W / 2;   PointsY[10] = H;

   PrintLineToPSFile("");
   PrintLineToPSFile("newpath");
   for(int i = 0; i < 11; i++)
   {
      char TempString[1000];
      sprintf(TempString, "%.1f %.1f ", LowerLeftX + PointsX[i] * Size, LowerLeftY + PointsY[i] * Size);
      
      if(i == 0)
         PrintLineToPSFile(string(TempString) + " moveto");
      else
         PrintLineToPSFile(string(TempString) + " lineto");
      
   }
   PrintLineToPSFile("closepath fill");

   stringstream box;
   box << "[ /Rect [" << LowerLeftX << " " << LowerLeftY << " "
      << LowerLeftX + Size << " " << LowerLeftY + Size << "]";

   PrintLineToPSFile("");
   PrintLineToPSFile(box.str());
   PrintLineToPSFile("/Color [0 0 1]");
   PrintLineToPSFile("/Dest /" + Destination);
   PrintLineToPSFile("/Subtype /Link");
   PrintLineToPSFile("/ANN pdfmark");
   PrintLineToPSFile("");
}

void PsFileHelper::SetAutomaticHomeButton(bool newvalue, string Destination, double X, double Y, double Size)
{
   AutomaticHomeButton = newvalue;
   HomeButtonDestination = Destination;
   HomeButtonX = X;
   HomeButtonY = Y;
   HomeButtonSize = Size;
}

void PsFileHelper::SetPageNumber(bool printnumber)
{
   PrintPageNumber = printnumber;
}

void PsFileHelper::InsertPageNumber(TCanvas &Canvas, int PageNumber)
{
   if(PageNumber < 0)
      PageNumber = NextPageNumber;

   Canvas.cd();

   TLatex Latex;
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.02);
   Latex.SetNDC(true);
   Latex.DrawLatex(0.98, 0.01, Form("%d", PageNumber));
}

#endif


