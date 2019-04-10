//================================================================================================
//
// Skim
//
//________________________________________________________________________________________________
//

//#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TH1F.h>                
#include <TCanvas.h>                
#include <TLegend.h> 
#include <THStack.h> 
#include <TKey.h> 
#include <TApplication.h>
//#endif


//=== MAIN MACRO ================================================================================================= 

int DoMixSamples( string inputfile1, string inputfile2, double weightFactor1, double weightFactor2, string outputfile) {

  //create output file
  TFile *outputFile = new TFile(outputfile.c_str(), "RECREATE");

  //loop over all TTrees in the file
  TFile *inputFile1 = TFile::Open(inputfile1.c_str(), "READ");
  assert(inputFile1);
  inputFile1->cd();
  inputFile1->Purge(); //purge unwanted TTree cycles in file
  TFile *inputFile2 = TFile::Open(inputfile2.c_str(), "READ");
  assert(inputFile2);
  inputFile2->cd();
  inputFile2->Purge(); //purge unwanted TTree cycles in file

  //  TIter nextkey(inputFile1->GetListOfKeys());
  //TKey *key;


  TTree *inputTree1 = (TTree*)inputFile1->Get("HggRazorLeptons");
  TTree *inputTree2 = (TTree*)inputFile2->Get("HggRazorLeptons");
  cout << "Processing tree " << inputTree1->GetName() << endl;
  
  //create new tree
  outputFile->cd();
  TTree *outputTree1 = inputTree1->CloneTree(0);  
  TTree *outputTree2 = inputTree2->CloneTree(0);  
  cout << "Events in the ntuple: " << inputTree1->GetEntries() << endl;
  
  uint run;
  uint event;
  float weight;
 
  //First File
  inputTree1->SetBranchAddress("weight", &weight);
    
  for (int n=0;n<inputTree1->GetEntries();n++) { 
    if (n%10000==0) cout << "Processed Event " << n << "\n";
    inputTree1->GetEntry(n);
    weight = weight * weightFactor1;
    outputTree1->Fill(); 
  }
    
  //Second File
  inputTree2->SetBranchAddress("weight", &weight);	    
  for (int n=0;n<inputTree2->GetEntries();n++) { 
    if (n%1000000==0) cout << "Processed Event " << n << "\n";
    inputTree2->GetEntry(n);
    weight = weight * weightFactor2;
    outputTree2->Fill(); 
  }
    
  TList *list = new TList;
  list->Add(outputTree1);
  list->Add(outputTree2);
  TTree *outputTree = TTree::MergeTrees(list);
  outputTree->SetName("HggRazorLeptons");
  outputTree->Write();

  cout << "Number of Input Events From File 1: " << inputTree1->GetEntries() << "\n";
  cout << "Number of Input Events From File 2: " << inputTree2->GetEntries() << "\n";
  cout << "Number of Output Events In File: " << outputTree->GetEntries() << "\n";



  //Merge the Histograms
  TH1F *NEvents1 = (TH1F*)inputFile1->Get("NEvents");
  TH1F *SumWeights1 = (TH1F*)inputFile1->Get("SumWeights");
  TH1F *SumScaleWeights1 = (TH1F*)inputFile1->Get("SumScaleWeights");
  TH1F *SumPdfWeights1 = (TH1F*)inputFile1->Get("SumPdfWeights");
  TH1F *NISRJets1 = (TH1F*)inputFile1->Get("NISRJets");
  TH1F *PtISR1 = (TH1F*)inputFile1->Get("PtISR");
  TH1F *NPV1 = (TH1F*)inputFile1->Get("NPV");

  TH1F *NEvents2 = (TH1F*)inputFile1->Get("NEvents");
  TH1F *SumWeights2 = (TH1F*)inputFile1->Get("SumWeights");
  TH1F *SumScaleWeights2 = (TH1F*)inputFile1->Get("SumScaleWeights");
  TH1F *SumPdfWeights2 = (TH1F*)inputFile1->Get("SumPdfWeights");
  TH1F *NISRJets2 = (TH1F*)inputFile1->Get("NISRJets");
  TH1F *PtISR2 = (TH1F*)inputFile1->Get("PtISR");
  TH1F *NPV2 = (TH1F*)inputFile1->Get("NPV");

  TH1F *NEvents = (TH1F*)NEvents1->Clone("NEvents");
  TH1F *SumWeights = (TH1F*)SumWeights1->Clone("SumWeights");
  TH1F *SumScaleWeights = (TH1F*)SumScaleWeights1->Clone("SumScaleWeights");
  TH1F *SumPdfWeights = (TH1F*)SumPdfWeights1->Clone("SumPdfWeights");
  TH1F *NISRJets = (TH1F*)NISRJets1->Clone("NISRJets");
  TH1F *PtISR = (TH1F*)PtISR1->Clone("PtISR");
  TH1F *NPV = (TH1F*)NPV1->Clone("NPV");

  for (int i=0; i< NEvents->GetXaxis()->GetNbins()+2; i++) {
    NEvents->SetBinContent( i, weightFactor1 * NEvents1->GetBinContent(i) + 
			    weightFactor2 * NEvents2->GetBinContent(i));
  }

  for (int i=0; i< SumWeights->GetXaxis()->GetNbins()+2; i++) {
    SumWeights->SetBinContent( i, weightFactor1 * SumWeights1->GetBinContent(i) + 
			    weightFactor2 * SumWeights2->GetBinContent(i));
  }

  for (int i=0; i< SumScaleWeights->GetXaxis()->GetNbins()+2; i++) {
    SumScaleWeights->SetBinContent( i, weightFactor1 * SumScaleWeights1->GetBinContent(i) + 
			    weightFactor2 * SumScaleWeights2->GetBinContent(i));
  }

  for (int i=0; i< SumPdfWeights->GetXaxis()->GetNbins()+2; i++) {
    SumPdfWeights->SetBinContent( i, weightFactor1 * SumPdfWeights1->GetBinContent(i) + 
			    weightFactor2 * SumPdfWeights2->GetBinContent(i));
  }

  for (int i=0; i< NISRJets->GetXaxis()->GetNbins()+2; i++) {
    NISRJets->SetBinContent( i, weightFactor1 * NISRJets1->GetBinContent(i) + 
			    weightFactor2 * NISRJets2->GetBinContent(i));
  }

  for (int i=0; i< PtISR->GetXaxis()->GetNbins()+2; i++) {
    PtISR->SetBinContent( i, weightFactor1 * PtISR1->GetBinContent(i) + 
			    weightFactor2 * PtISR2->GetBinContent(i));
  }

  for (int i=0; i< NPV->GetXaxis()->GetNbins()+2; i++) {
    NPV->SetBinContent( i, weightFactor1 * NPV1->GetBinContent(i) + 
			    weightFactor2 * NPV2->GetBinContent(i));
  }

  NEvents->Write(); 
  SumWeights->Write(); 
  SumScaleWeights->Write();
  SumPdfWeights->Write();
  NISRJets->Write();
  PtISR->Write();
  NPV->Write();
	
  inputFile1->Close();
  inputFile2->Close();
  cout << "Closing output file." << endl;
  outputFile->Purge(); 
  outputFile->Close();
  delete outputFile;
  delete list;
  //gApplication->Terminate();
  return 0;
}


void MixSignalSamples() {

 cout << "test\n";

  vector<double> ZBRScan;
  vector<string> ZBRScanLabel;
  // ZBRScan.push_back( 0.1 ); ZBRScanLabel.push_back("BR1090");
  // ZBRScan.push_back( 0.2 ); ZBRScanLabel.push_back("BR2080");
  // ZBRScan.push_back( 0.3 ); ZBRScanLabel.push_back("BR3070");
  // ZBRScan.push_back( 0.4 ); ZBRScanLabel.push_back("BR4060");
  ZBRScan.push_back( 0.5 ); ZBRScanLabel.push_back("BR5050");
  // ZBRScan.push_back( 0.6 ); ZBRScanLabel.push_back("BR6040");
  // ZBRScan.push_back( 0.7 ); ZBRScanLabel.push_back("BR7030");
  // ZBRScan.push_back( 0.8 ); ZBRScanLabel.push_back("BR8020");
  // ZBRScan.push_back( 0.9 ); ZBRScanLabel.push_back("BR9010");

 
  int samples[36] = { 127,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,625,650,675,700,725,750,775,800,825,850,875,900,925,950,975,1000};
  for (uint s = 0; s < ZBRScan.size(); s++) {
    for (int i=1; i<36; i++) {
      cout << "Sample : " << i << " : " << samples[i] << "\n";
      DoMixSamples(Form("/eos/cms/store/group/phys_susy/razor/Run2Analysis/SusyEwkHgg/signal/jobs_0908/TChiHHjobs/combined/TChiHH_%i_1pb_weighted.root",samples[i]),Form("/eos/cms/store/group/phys_susy/razor/Run2Analysis/SusyEwkHgg/signal/jobs_0908/TChiHZjobs/combined/TChiHZ_%i_1pb_weighted.root",samples[i]), ZBRScan[s]*ZBRScan[s], 2*ZBRScan[s]*(1-ZBRScan[s]), Form("/eos/cms/store/group/phys_susy/razor/Run2Analysis/SusyEwkHgg/signal/jobs_0908/TChiHHHZMixedSamples/TChiHHHZ_%s_%i_1pb_weighted.root",ZBRScanLabel[s].c_str(),samples[i]));   
    }
  
  }
  
}


