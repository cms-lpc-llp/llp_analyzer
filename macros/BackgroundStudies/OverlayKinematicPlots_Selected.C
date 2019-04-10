#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdexcept>
using namespace std;

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TFrame.h"

#include "SimpleTable.h"
#include "TKey.h"

#include "PlotHelper2.h"
#include "SetStyle.h"

#define LUMI 5000

int main();
void OverlayHistogram(PsFileHelper &PsFile, vector<TFile *> Files, vector<double> Xsection, string HistogramName);

int main()
{
   SetStyle();

   PsFileHelper PsFile("DataMCComparison.ps");
   PsFile.AddTextPage("Data-MC comparison!");

   //Read the input files
   vector<TFile *> Files;
   string line;
   ifstream inputFile("inputFilesListQCD25vs50.txt");
    if (!inputFile.is_open())
        throw runtime_error("File Not Found!");    
    Files.push_back(new TFile("test/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola.root") ); // dummy placeholder, data will be index 0
    while (getline(inputFile, line))
      Files.push_back(new TFile(line.c_str()));
       
   //Get CrossSection
   vector<double> Xsection;
   SimpleTable xstab("data/xSections.dat");
   Xsection.push_back(1.0); //data weight is 1.0
   for(int i=1; i<Files.size(); i++)
     {
       string datasetName = Files[i]->GetName();
       datasetName.erase( datasetName.begin(),  datasetName.begin() + datasetName.find_first_of("/") + 1 ); 
       datasetName.erase( datasetName.find_last_of("."), datasetName.find_last_of(".") + 3 );

       Xsection.push_back(xstab.Get(datasetName.c_str()));
       cout << "Cross Section for " << datasetName<<" is = "<<Xsection[i] << "\n";
     }

   // Make the plots
   OverlayHistogram(PsFile, Files, Xsection, "h_MR_MultiJet");
   OverlayHistogram(PsFile, Files, Xsection, "h_Rsq_MultiJet");

   OverlayHistogram(PsFile, Files, Xsection, "h_MR_EleJet");
   OverlayHistogram(PsFile, Files, Xsection, "h_Rsq_EleJet");

   for(int i = 0; i < (int)Files.size(); i++)
   {
      Files[i]->Close();
      delete Files[i];
   }
   Files.clear();

   PsFile.AddTimeStampPage();
   PsFile.Close();
}

void OverlayHistogram(PsFileHelper &PsFile, vector<TFile *> Files, vector<double > Xsection, string HistogramName)
{
  // Get the histograms from ROOT files
   vector<TH1F *> Histograms;
   for(int i = 0; i < (int)Files.size(); i++)
     {
       if(!((TH1F *)Files[i]->Get(HistogramName.c_str()))) {cout<<"the histogram: " << HistogramName.c_str() << " is not present in file number "<< i<<endl; return;}
       Histograms.push_back((TH1F *)Files[i]->Get(HistogramName.c_str()));
     }

  // Get the normalization factor for each process
   vector<Double_t> Factor;
   for(int i = 0; i < (int)Files.size(); i++)
     Factor.push_back( Xsection[i] * LUMI / ((TH1F *)Files[i]->Get("NEvents"))->GetEntries() );

   // Scale MC to Data luminosity
   Factor[0] = 1.0; // data
   for(int i = 1; i < (int)Files.size(); i++)
     if(Histograms[i]->Integral())
       Histograms[i]->Scale( Factor[i] );
   
   // Set up the colors/styles
   for(int i = 1; i < (int)Files.size(); i++)
     {
       string datasetName = Files[i]->GetName(); 

       if ( datasetName.find("TTJets") != std::string::npos )
	 { Histograms[i]->SetLineColor(2); Histograms[i]->SetFillColor(2); } 
       if ( datasetName.find("QCD") != std::string::npos )
	 { Histograms[i]->SetLineColor(3); Histograms[i]->SetFillColor(3); } 
       if ( datasetName.find("DYJetsToLL") != std::string::npos )
	 { Histograms[i]->SetLineColor(4); Histograms[i]->SetFillColor(4); } 
       if ( datasetName.find("WJetsToLNu") != std::string::npos )
	 { Histograms[i]->SetLineColor(5); Histograms[i]->SetFillColor(5); } 
       if ( datasetName.find("ZJetsToNuNu") != std::string::npos )
	 { Histograms[i]->SetLineColor(6); Histograms[i]->SetFillColor(6); } 
     }
  
   // Set up the Canvas and split into two pads
   TCanvas *C = new TCanvas("C", "",0,0,800,600);
   
   // Group backgrounds together
   int indexQCD = 0; int indexDYJetsToLL = 0; int indexZJetsToNuNu = 0; int indexWJetsToLNu = 0;
   int counterQCD = 0; int counterDYJetsToLL = 0; int counterZJetsToNuNu = 0; int counterWJetsToLNu = 0;
   for(int i = 1; i < (int)Files.size(); i++)
     {
       string datasetName = Files[i]->GetName(); 

       if ( datasetName.find("QCD") != std::string::npos && counterQCD > 0 )
	 { Histograms[indexQCD]->Add(Histograms[i], 1.0); counterQCD++; }
       if ( datasetName.find("QCD") != std::string::npos && counterQCD == 0 )
	 { counterQCD++; indexQCD = i; }

       if ( datasetName.find("DYJetsToLL") != std::string::npos && counterDYJetsToLL > 0 )
	 { Histograms[indexDYJetsToLL]->Add(Histograms[i], 1.0); counterDYJetsToLL++; }
       if ( datasetName.find("DYJetsToLL") != std::string::npos && counterDYJetsToLL == 0 )
	 { counterDYJetsToLL++; indexDYJetsToLL = i; }

       if ( datasetName.find("WJetsToLNu") != std::string::npos && counterWJetsToLNu > 0 )
	 { Histograms[indexWJetsToLNu]->Add(Histograms[i], 1.0); counterWJetsToLNu++; }
       if ( datasetName.find("WJetsToLNu") != std::string::npos && counterWJetsToLNu == 0 )
	 { counterWJetsToLNu++; indexWJetsToLNu = i; }

       if ( datasetName.find("ZJetsToNuNu") != std::string::npos && counterZJetsToNuNu > 0 )
	 { Histograms[indexZJetsToNuNu]->Add(Histograms[i], 1.0); counterZJetsToNuNu++; }
       if ( datasetName.find("ZJetsToNuNu") != std::string::npos && counterZJetsToNuNu == 0 )
	 { counterZJetsToNuNu++; indexZJetsToNuNu = i; }
     }

   // Set up the legend
   TLegend legend(0.70, 0.70, 0.90, 0.90);
   legend.AddEntry(Histograms[1], "TTbar", "f");
   legend.AddEntry(Histograms[indexQCD], "QCD", "f");
   legend.AddEntry(Histograms[indexDYJetsToLL], "DYJetsToLL", "f");
   legend.AddEntry(Histograms[indexWJetsToLNu], "WJetsToLNu", "f");
   legend.AddEntry(Histograms[indexZJetsToNuNu], "ZJetsToNuNu", "f");
   legend.SetFillColor(0);
   
   cout<<" Top "<<Histograms[1]->Integral()<<"\n"
       <<" Znn "<<Histograms[indexZJetsToNuNu]->Integral()<<"\n"
       <<" DY "<<Histograms[indexDYJetsToLL]->Integral()<<"\n"
       <<" WJ "<<Histograms[indexWJetsToLNu]->Integral()<<"\n"
       <<" QCD "<<Histograms[indexQCD]->Integral()<<"\n"
       <<endl;

   // make the stack of backgrounds 
   THStack hs("hs",HistogramName.c_str());
   hs.Add(Histograms[1]);
   hs.Add(Histograms[indexDYJetsToLL]);
   hs.Add(Histograms[indexWJetsToLNu]);
   hs.Add(Histograms[indexZJetsToNuNu]);
   hs.Add(Histograms[indexQCD]);

   hs.Draw("bar");
   hs.SetMinimum(0.000001);
   hs.GetHistogram()->SetMinimum(0.000001);
   legend.Draw("same");
   C->Update();
   
   PsFile.AddCanvas(C);      
}
