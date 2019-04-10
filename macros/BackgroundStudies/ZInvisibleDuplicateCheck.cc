#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TTreeFormula.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPad.h"
#include "TChain.h"

using namespace std;

void ZInvisibleDuplicateCheck(){
    gROOT->SetBatch();

    //make a map to hold the event numbers
    map<string, bool> eventID;

    //ifstream ifs("lists/razorNtuplerV1p6-Run1/Data_SinglePhoton_Run2012B.cern.txt");
    ifstream ifs("lists/razorNtuplerV1p6-Run1/Data_SinglePhoton_Run2012C.cern.txt");
    //ifstream ifs("lists/razorNtuplerV1p6-Run1/Data_SinglePhotonParked_Run2012D.cern.txt");
    //ifstream ifs("lists/razorNtuplerV1p6-Run1/Photon.cern.txt");
    string curFilename;
    vector<string> inputLines;
    TChain *theChain = new TChain("ntuples/RazorEvents");
    while(getline(ifs, curFilename)){
        theChain->Add(curFilename.c_str());
        cout << curFilename << endl;
    }
    int run, lumi, event;
    theChain->SetBranchAddress("runNum", &run);
    theChain->SetBranchAddress("lumiNum", &lumi);
    theChain->SetBranchAddress("eventNum", &event);

    string space = " ";
    for(int iEvent = 0; iEvent < theChain->GetEntries(); iEvent++){
        if(iEvent % 1000000 == 0) cout << "Processing event " << iEvent << endl;
        theChain->GetEntry(iEvent);
        string id = to_string(run) + space + to_string(lumi) + space + to_string(event);
        if(eventID.find(id) == eventID.end()){
            eventID[id] = true; 
        }
        else{
            cout << "Duplicate found! " << id << endl;
        }
    }
    delete theChain;
}

int main(){
    ZInvisibleDuplicateCheck();
    return 0;
}


