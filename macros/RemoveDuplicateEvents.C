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

int RemoveDuplicateEvents( string inputfile, string outputfile) {

  //create output file
  TFile *outputFile = new TFile(outputfile.c_str(), "RECREATE");

  //loop over all TTrees in the file
  TFile *inputFile = TFile::Open(inputfile.c_str(), "READ");
  assert(inputFile);
  inputFile->cd();
  inputFile->Purge(); //purge unwanted TTree cycles in file
  TIter nextkey(inputFile->GetListOfKeys());
  TKey *key;
  while((key = (TKey*)nextkey())){
    string className = key->GetClassName();
    cout << "Getting key from file.  Class type: " << className << endl;
    if(className.compare("TTree") != 0){
      cout << "Skipping key (not a TTree)" << endl;
      continue;
    }

    TTree *inputTree = (TTree*)key->ReadObj();
    cout << "Processing tree " << inputTree->GetName() << endl;

    //create new tree
    outputFile->cd();
    TTree *outputTree = inputTree->CloneTree(0);  
    cout << "Events in the ntuple: " << inputTree->GetEntries() << endl;

    uint run;
    uint event;
    inputTree->SetBranchAddress("run", &run);
    inputTree->SetBranchAddress("event", &event);
	    
    map<pair<uint,uint>, bool > processedRunEvents;
    for (int n=0;n<inputTree->GetEntries();n++) { 
      if (n%1000000==0) cout << "Processed Event " << n << "\n";
      inputTree->GetEntry(n);

      if(processedRunEvents.find(make_pair(run, event)) == processedRunEvents.end()){ 
          //not duplicate
          processedRunEvents[make_pair(run, event)] = true;
          outputTree->Fill(); 
      }
      else{
          cout << "Duplicate event " << run << " " << event << endl;
      }
    }

    cout << "Number of Input Events: " << inputTree->GetEntries() << "\n";
    cout << "Number of Output Events: " << outputTree->GetEntries() << "\n";

    //save
    outputTree->Write();
    inputFile->cd();
  }
	
  inputFile->Close();
  cout << "Closing output file." << endl;
  outputFile->Close();
  delete outputFile;
  gApplication->Terminate();
  return 0;
}

// void RemoveDuplicateEvents( string inputfile, string outputfile) {
  
 
//   //create output file
//   TFile *outputFile = new TFile(outputfile.c_str(), "RECREATE");

//   //loop over all TTrees in the file and add the weight branch to each of them
//   TFile inputFile(inputfile.c_str(), "READ");
//   inputFile.cd();
//   inputFile.Purge(); //purge unwanted TTree cycles in file
//   TIter nextkey(inputFile.GetListOfKeys());
//   TKey *key;
//   cout << "here1\n";

//   TTree *inputTree = (TTree*)inputFile.Get("zz4lTree/probe_tree");
//   cout << "Processing tree " << inputTree->GetName() << endl;
  
//   //create new normalized tree
//   outputFile->cd();
//   TTree *outputTree = inputTree->CloneTree(0);  
//   cout << "Events in the ntuple: " << inputTree->GetEntries() << endl;
  
//   uint run;
//   uint event;
//   inputTree->SetBranchAddress("run", &run);
//   inputTree->SetBranchAddress("event", &event);
  
//   vector<pair<uint,uint> > processedRunEventList;
  
//   for (int n=0;n<inputTree->GetEntries();n++) { 
//     if (n%1000000==0) cout << "Processed Event " << n << "\n";
//     inputTree->GetEntry(n);

//     bool isDuplicate = false;
//     for( int i=0; i < int(processedRunEventList.size()) ; ++i) {
//       if (processedRunEventList[i].first == run && processedRunEventList[i].second == event) isDuplicate = true;
//     }
//     if (!isDuplicate) {
//       pair<uint,uint> newRunEvent = pair<uint,uint>(run, event);
//       processedRunEventList.push_back(newRunEvent);
//       outputTree->Fill(); 
//     }
//   }

//   cout << "Number of Input Events: " << inputTree->GetEntries() << "\n";
//   cout << "Number of Output Events: " << outputTree->GetEntries() << "\n";

//   //save
//   outputTree->Write();
//   inputFile.cd();
	
//   inputFile.Close();
//   cout << "Closing output file." << endl;
//   outputFile->Close();
//   delete outputFile;

// }

// void RemoveDuplicateEvents( string inputfile) {
  
 
//   //loop over all TTrees in the file and add the weight branch to each of them
//   TFile inputFile(inputfile.c_str(), "READ");
//   inputFile.cd();
//   inputFile.Purge(); //purge unwanted TTree cycles in file
//   TIter nextkey(inputFile.GetListOfKeys());
//   TKey *key;
//   cout << "here1\n";

//   TTree *inputTree = (TTree*)inputFile.Get("probe_tree");
//   cout << "Processing tree " << inputTree->GetName() << endl;
  
//   //create new normalized tree
//   uint run;
//   uint lumi;
//   uint event;
//   float mass;
//   inputTree->SetBranchAddress("run", &run);
//   inputTree->SetBranchAddress("lumi", &lumi);
//   inputTree->SetBranchAddress("event", &event);
//   inputTree->SetBranchAddress("mass", &mass);
  
//   vector<pair<uint,uint> > processedRunEventList;
  
//   for (int n=0;n<inputTree->GetEntries();n++) { 
//     if (n%1000000==0) cout << "Processed Event " << n << "\n";
//     inputTree->GetEntry(n);

//     if (mass > 115 && mass < 135) {
//       cout << run << " " << lumi << " " << event << " " << mass << "\n";
//     }

//   }

//   cout << "Number of Input Events: " << inputTree->GetEntries() << "\n";

//   //save
//   inputFile.cd();
	
//   inputFile.Close();
  
// }

