//================================================================================================
//
// Skim
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
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

#include "RazorAnalyzer/include/ControlSampleEvents.h"

#endif


//=== MAIN MACRO ================================================================================================= 


void SkimRazorControlSample_TightPlusVetoIDLepton( string inputfile, string outputfile) {
  
 
  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  TFile *outputFile = new TFile(outputfile.c_str(), "RECREATE");

  ControlSampleEvents *events = new ControlSampleEvents;
  events->LoadTree(inputfile.c_str());

  //create new normalized tree
  outputFile->cd();
  TTree *outputTree = events->tree_->CloneTree(0);  
  cout << "Events in the original ntuple: " << events->tree_->GetEntries() << "\n";

  for(UInt_t ientry=0; ientry < events->tree_->GetEntries(); ientry++) {       	
    events->tree_->GetEntry(ientry);      
    if (ientry % 1000000 == 0) cout << "Event " << ientry << endl;      
    
    if ( events->lep1.Pt() > 20 && events->lep2.Pt() > 20 && events->lep1PassTight && events->lep2PassVetoID
	 ) {
      outputTree->Fill();
    }
  }

  cout << "Events Passing Skimming: " << outputTree->GetEntries() << "\n";
  cout << "Skim Efficiency : " <<  double(outputTree->GetEntries()) / events->tree_->GetEntries() << "\n";
  outputTree->Write();
  outputFile->Close();
  delete outputFile;

}

