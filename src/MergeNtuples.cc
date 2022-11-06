#include <fstream>
#include <sstream>
#include <iterator>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TKey.h"
#include <assert.h>
#include <TRandom3.h>
#include "TTreeFormula.h"
#include <iostream>
using namespace std;


std::string ParseCommandLine( int argc, char* argv[], std::string opt )
{
  for (int i = 1; i < argc; i++ )
    {
      std::string tmp( argv[i] );
      if ( tmp.find( opt ) != std::string::npos )
        {
          if ( tmp.find( "=" )  != std::string::npos ) return tmp.substr( tmp.find_last_of("=") + 1 );
	  if ( tmp.find( "--" ) != std::string::npos ) return "yes";
	}
    }
  
  return "";
};




//get list of files to open, add normalization branch to the tree in each file
int main(int argc, char* argv[]) {

    //parse input list to get names of ROOT files
    if(argc < 3){
        cerr << "usage NormalizeNtuple [inputfile1] [inputfile2] [outputfile]" << endl;
        return -1;
    }
    string inputfilename1(argv[1]);
    string inputfilename2(argv[2]);
    string outputfilename(argv[3]);

  
  
 
    //create output file
    TFile *outputFile = new TFile(outputfilename.c_str(), "RECREATE");

    //open files
    TFile *inputFile1 = TFile::Open(inputfilename1.c_str(), "READ");
    TTree *inputTree1 = (TTree*)inputFile1->Get("MuonSystem");
    if (!inputTree1) cout << "Input Tree not found in file " << inputfilename1 << "\n";
    assert(inputTree1);
    TFile *inputFile2 = TFile::Open(inputfilename2.c_str(), "READ");
    TTree *inputTree2 = (TTree*)inputFile2->Get("tree");
    if (!inputTree2) cout << "Input Tree not found in file " << inputfilename2 << "\n";
    assert(inputTree2);

    uint NEventsTree1 = inputTree1->GetEntries();
    uint NEventsTree2 = inputTree2->GetEntries();
    //NEventsTree1 = 20000;

    //*****************************************************************************************
    //Make map of event number in tree 1 to event number in tree 2
    //*****************************************************************************************
    int numberOfNonMatchedEvents = 0;
    std::map<uint,uint> EventIndexToEventIndexMap;
    //unsigned long long _event;
    uint _event;
    uint _run;
    uint run;
    unsigned long long event;
    inputTree1->SetBranchAddress("runNum", &_run);
    inputTree1->SetBranchAddress("evtNum", &_event);
    inputTree2->SetBranchAddress("run", &run);
    inputTree2->SetBranchAddress("event", &event);

    

    //loop over tree2
    std::vector<std::pair<uint,unsigned long long> > eventList2;
    for (uint m=0; m < NEventsTree2;m++) { 
      inputTree2->GetEntry(m);    
      std::pair<uint,unsigned long long> p (run, event);
      eventList2.push_back(p);
    }
    
    cout << "Total Entries: " << NEventsTree1 << "\n";
    //loop over tree1
    for (uint n=0;n<NEventsTree1;n++) { 
      if ( n % 100 == 0) cout << "Event " << n << "\n";
      inputTree1->GetEntry(n);    

      //cout << "Event " << n << " : " << _run << " : " << _event << "\n";

      //loop over tree2
      bool matchFound = false;
      for (uint m=0; m < eventList2.size();m++) { 

	if (eventList2[m].first != _run) continue;
	if (eventList2[m].second != _event) continue;

	//if matched, then add entry to the map
	EventIndexToEventIndexMap[n] = m;
	//cout << "Match : " << n << " --> " << m << "\n";
	matchFound = true;
	break;
      }            
      if (!matchFound) {
	numberOfNonMatchedEvents++;
	//cout << "Event " << _run << ":" << _event << " Not Matched\n";
      }
      
    }
 
 
    // for (uint n=0;n<NEventsTree1;n++) { 
    //   cout << "Map [" << n << "] --> " << EventIndexToEventIndexMap[n] << "\n";
    // }

    cout << numberOfNonMatchedEvents << " Events out of " << inputTree1->GetEntries()  << " were not matched\n";


    //*****************************************************************************************
    //Produce Output Tree
    //*****************************************************************************************       
    cout << "Processing tree " << inputTree1->GetName() << endl;

    //create output tree
    outputFile->cd();
    TTree *outputTree = inputTree1->CloneTree(0);  
    cout << "Events in the ntuple: " << inputTree1->GetEntries() << endl;

    //add branches to output
    const uint ArraySizeLimit = 100;
    uint nTau = 0;
    float tauPt[ArraySizeLimit];
    for (uint i = 0 ; i < ArraySizeLimit; i++) {
      tauPt[i] = -999;
    }
    outputTree->Branch("nTau", &nTau, "nTau/i");
    outputTree->Branch("tauPt", &tauPt, "tauPt[nTau]/F");
    inputTree2->SetBranchAddress("nTau", &nTau);
    inputTree2->SetBranchAddress("tauPt", &tauPt);


    //loop over Tree1 and add all the branches from tree2
    for (uint n=0;n<NEventsTree1;n++) { 
      if (n%1000==0) cout << "Processed Event " << n << "\n";
      inputTree1->GetEntry(n);
      
      inputTree2->GetEntry(EventIndexToEventIndexMap[n]);
      //cout << "Event : " << run << " " << event << " : " << nTau << "\n";

      outputTree->Fill();		
    }
    
    //save
    outputTree->Write();
    inputFile1->cd();
    inputFile1->Close();
    inputFile2->cd();
    inputFile2->Close();
    cout << "Closing output file." << endl;
    outputFile->Close();
    delete outputFile;
}
