#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"

using namespace std;

int main(int argc, char* argv[]);
void skiming(string filename);

int main(int argc, char* argv[])
{
  string inputFileName;
  inputFileName = argv[1];

  skiming(inputFileName);
}

void skiming(string filename)
  {
    //Get old file, old tree and set top branch address
    TFile *oldfile = new TFile(filename.c_str());
    TTree *oldtree = (TTree*)oldfile->Get("RazorInclusive");
    Long64_t nentries = oldtree->GetEntries();

    TH1F NEvents = *(TH1F *)oldfile->Get("NEvents");
    
    int nSelectedPhotons = 0.;
    
    oldtree->SetBranchAddress("nSelectedPhotons",&nSelectedPhotons);
    
    //Create a new file + a clone of old tree in new file
    TFile *newfile = new TFile(Form("%s_skim.root", filename.c_str()),"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    
    for (Long64_t i=0;i<nentries; i++) {
      oldtree->GetEntry(i);
      
      if(i % 10000000 == 0) cout << "Processing entry " << i << " out of "<<nentries<< endl;
      
      if (nSelectedPhotons > 0) newtree->Fill();
      
    }

    NEvents.Write();
    newtree->AutoSave();
    delete oldfile;
    delete newfile;
  }
