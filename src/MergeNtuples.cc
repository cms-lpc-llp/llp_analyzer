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
    if(argc < 4){
        cerr << "usage MergeNtuple [inputfile1] [inputfile2] [outputfile] [isData]" << endl;
        return -1;
    }
    string inputfilename1(argv[1]);
    string inputfilename2(argv[2]);
    string outputfilename(argv[3]);
    string isdata(argv[4]); //1 (yes) or 0 (no)

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
    //NEventsTree1 = 10000;

    //*****************************************************************************************
    //Make map of event number in tree 1 to event number in tree 2
    //*****************************************************************************************
    int numberOfNonMatchedEvents = 0;
    std::map<uint,uint> EventIndexToEventIndexMap;
    unsigned long long _event;
    uint _run;
    uint run;
    unsigned long long event;
    float met;
    bool Flag2_all;
    int nCscRings;
    int nDtRings;
    int nCscRechitClusters;
    int nDtRechitClusters;
    bool METNoMuTrigger;
    int   gLLP_csc[2];
    int   gLLP_dt[2];
    inputTree1->SetBranchAddress("runNum", &_run);
    inputTree1->SetBranchAddress("evtNum", &_event);
    inputTree1->SetBranchAddress("met", &met);
    inputTree1->SetBranchAddress("Flag2_all", &Flag2_all);
    inputTree1->SetBranchAddress("nCscRings", &nCscRings);
    inputTree1->SetBranchAddress("nDtRings", &nDtRings);
    inputTree1->SetBranchAddress("nCscRechitClusters", &nCscRechitClusters);
    inputTree1->SetBranchAddress("nDtRechitClusters" , &nDtRechitClusters);
    inputTree1->SetBranchAddress("METNoMuTrigger", &METNoMuTrigger);
    inputTree1->SetBranchAddress("gLLP_csc",    gLLP_csc);
    inputTree1->SetBranchAddress("gLLP_dt",    gLLP_dt);
    inputTree2->SetBranchAddress("run", &run);
    inputTree2->SetBranchAddress("event", &event);

    //loop over tree2
    std::vector<std::pair<uint,unsigned long long> > eventList2;
    std::vector<bool> matchedevent;
    for (uint m=0; m < NEventsTree2;m++) { 
      inputTree2->GetEntry(m);    
      std::pair<uint,unsigned long long> p (run, event);
      eventList2.push_back(p);
    }
    
    cout << "Total Entries: " << NEventsTree1 << "\n";
    //loop over tree1
    for (uint n=0;n<NEventsTree1;n++) { 
      if ( n % 1000 == 0) cout << "Event " << n << "\n";
      inputTree1->GetEntry(n);    

      //loop over tree2
      bool matchFound = false;
      for (uint m=0; m < eventList2.size();m++) 
      { 
          	if (eventList2[m].first != _run) continue;
          	if (eventList2[m].second != _event) continue;
          	//if matched, then add entry to the map
          	EventIndexToEventIndexMap[n] = m;
          	//cout << "Match : " << n << " --> " << m << "\n";
          	matchFound = true;
          	break;
      }            

      if (!matchFound) 
      { 
          numberOfNonMatchedEvents++;
          matchedevent.push_back(false);
        	cout << "Event " << _run << ":" << _event <<" Not Matched\n";
      }
      else{
          matchedevent.push_back(true);
      }
      
    }
 
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
    float tauM[ArraySizeLimit];
    float tauPt[ArraySizeLimit];
    float tauEta[ArraySizeLimit];
    float tauPhi[ArraySizeLimit];
    float tauNeutralIso[ArraySizeLimit];
    float tauChargedIso[ArraySizeLimit];
    float tauDeltaR[ArraySizeLimit];
    int   tauDecayMode[ArraySizeLimit];
    bool  tau_IsVVVLoose[ArraySizeLimit];
    bool  tau_IsVVLoose[ArraySizeLimit];
    bool  tau_IsVLoose[ArraySizeLimit];
    bool  tau_IsLoose[ArraySizeLimit];
    bool  tau_IsMedium[ArraySizeLimit];
    bool  tau_IsTight[ArraySizeLimit];
    bool  tau_IsVTight[ArraySizeLimit];
    bool  tau_IsVVTight[ArraySizeLimit];

    for (uint i = 0 ; i < ArraySizeLimit; i++) {
      tauM[i]          = -999.;
      tauPt[i]         = -999.;
      tauEta[i]        = -999.;
      tauPhi[i]        = -999.;
      tauDeltaR[i]     =  999.;
      tauChargedIso[i] = 0.0;
      tauNeutralIso[i] = 0.0;
      tauDecayMode[i]  = 0.0;
      tau_IsVVVLoose[i]= false;
      tau_IsVVLoose[i] = false;
      tau_IsVLoose[i]  = false;
      tau_IsLoose[i]   = false;
      tau_IsMedium[i]  = false;
      tau_IsTight[i]   = false;
      tau_IsVTight[i]  = false;
      tau_IsVVTight[i] = false;
    }
    outputTree->Branch("nTau",   &nTau, "nTau/i");
    outputTree->Branch("tauM",   &tauM, "tauM[nTau]/F");
    outputTree->Branch("tauPt",  &tauPt, "tauPt[nTau]/F");
    outputTree->Branch("tauEta", &tauEta,                "tauEta[nTau]/F");
    outputTree->Branch("tauPhi", &tauPhi,                "tauPhi[nTau]/F");
    outputTree->Branch("tauDeltaR", &tauDeltaR,           "tauDeltaR[nTau]/F");
    outputTree->Branch("tauChargedIso", &tauChargedIso,   "tauChargedIso[nTau]/F");
    outputTree->Branch("tauNeutralIso", &tauNeutralIso,   "tauNeutralIso[nTau]/F");
    outputTree->Branch("tauDecayMode", &tauDecayMode,     "tauDecayMode[nTau]/I");
    outputTree->Branch("tau_IsVVVLoose", &tau_IsVVVLoose, "tau_IsVVVLoose[nTau]/O");
    outputTree->Branch("tau_IsVVLoose", &tau_IsVVLoose,   "tau_IsVVLoose[nTau]/O");
    outputTree->Branch("tau_IsVLoose", &tau_IsVLoose,     "tau_IsVLoose[nTau]/O");
    outputTree->Branch("tau_IsLoose", &tau_IsLoose,       "tau_IsLoose[nTau]/O");
    outputTree->Branch("tau_IsMedium", &tau_IsMedium,     "tau_IsMedium[nTau]/O");
    outputTree->Branch("tau_IsTight", &tau_IsTight,       "tau_IsTight[nTau]/O");
    outputTree->Branch("tau_IsVTight", &tau_IsVTight,     "tau_IsVTight[nTau]/O");
    outputTree->Branch("tau_IsVVTight", &tau_IsVVTight,   "tau_IsVVTight[nTau]/O");

    inputTree2->SetBranchAddress("nTau", &nTau);
    inputTree2->SetBranchAddress("tauM", &tauM);
    inputTree2->SetBranchAddress("tauPt", &tauPt);
    inputTree2->SetBranchAddress("tauEta", &tauEta);
    inputTree2->SetBranchAddress("tauPhi", &tauPhi);
    inputTree2->SetBranchAddress("tauDeltaR", &tauDeltaR);
    inputTree2->SetBranchAddress("tauChargedIso", &tauChargedIso);
    inputTree2->SetBranchAddress("tauNeutralIso", &tauNeutralIso);
    inputTree2->SetBranchAddress("tauDecayMode", &tauDecayMode);
    inputTree2->SetBranchAddress("tau_IsVVVLoose", &tau_IsVVVLoose);
    inputTree2->SetBranchAddress("tau_IsVVLoose", &tau_IsVVLoose);
    inputTree2->SetBranchAddress("tau_IsVLoose", &tau_IsVLoose);
    inputTree2->SetBranchAddress("tau_IsLoose", &tau_IsLoose);
    inputTree2->SetBranchAddress("tau_IsMedium", &tau_IsMedium);
    inputTree2->SetBranchAddress("tau_IsTight", &tau_IsTight);
    inputTree2->SetBranchAddress("tau_IsVTight", &tau_IsVTight);
    inputTree2->SetBranchAddress("tau_IsVVTight", &tau_IsVVTight);

    //*****************************************************************************************
    //Produce Output Histogram
    //***************************************************************************************** 
    //histogram containing total number of processed events (for normalizations)
    TH1F *CutFlow       = new TH1F("CutFlow"   , "CutFlow"         , 10, 0, 10);
    TH1F *Acceptance    = new TH1F("Acceptance", "Acceptance"      , 10, 0, 10);
    TH1F *AcceptanceCsc = new TH1F("AcceptanceCsc", "AcceptanceCsc", 10, 0, 10);
    TH1F *AcceptanceDt  = new TH1F("AcceptanceDt", "AcceptanceDt"  , 10, 0, 10);

    //loop over Tree1 and add all the branches from tree2
    for (uint n=0;n<NEventsTree1;n++) { 
      if (n%1000==0) cout << "Processed Event " << n << "\n";

      //Check if found a match
      if(matchedevent[n] == false) continue;

      //Get entries
      inputTree1->GetEntry(n);
      inputTree2->GetEntry(EventIndexToEventIndexMap[n]);

      //Cuts for skim
      //Bin 1: Total number of events or acceptance
      CutFlow->Fill(0); 
      bool acceptance=false;
      bool acceptancecsc=false;    
      bool acceptancedt=false;      
      if (isdata != "yes")
      {
        if (gLLP_csc[0]==true || gLLP_csc[1]==true) acceptancecsc=true;
        if (gLLP_dt[0]==true  || gLLP_dt[1]==true)  acceptancedt=true;
        if (acceptancecsc==true || acceptancedt==true) acceptance=true;
        if (acceptance)    Acceptance->Fill(0);
        if (acceptancecsc) AcceptanceCsc->Fill(0);
        if (acceptancedt)  AcceptanceDt->Fill(0);
      } 
      //Bin 2: MET Trigger + MET cut
      if ( METNoMuTrigger==false ) continue;
      if ( met < 200 ) continue;
      CutFlow->Fill(1); 
      if(acceptance)    Acceptance->Fill(1);
      if(acceptancecsc) AcceptanceCsc->Fill(1);
      if(acceptancedt)  AcceptanceDt->Fill(1);

      //Bin 3: MET filters 
      if (Flag2_all==false) continue;
      CutFlow->Fill(2); 
      if(acceptance) Acceptance->Fill(2);
      if(acceptancecsc) AcceptanceCsc->Fill(2);
      if(acceptancedt) AcceptanceDt->Fill(2);

      //Bin 4: Number CSC+DT rings < 10
      if ( (nCscRings + nDtRings) > 10) continue;
      CutFlow->Fill(3);
      if(acceptance) Acceptance->Fill(3);
      if(acceptancecsc) AcceptanceCsc->Fill(3);
      if(acceptancedt) AcceptanceDt->Fill(3);

      //Bin 5: At least 1 cluster (CSC or DT)
      if ( (nCscRechitClusters + nDtRechitClusters) == 0 ) continue;
      CutFlow->Fill(4);
      if(acceptance) Acceptance->Fill(4);
      if(acceptancecsc) AcceptanceCsc->Fill(4);
      if(acceptancedt) AcceptanceDt->Fill(4);      

      //Fill out tree if conditions are met
      outputTree->Fill();		
    }
    //save information
    outputTree->Write();
    CutFlow->Write();
    if (isdata != "yes") Acceptance->Write();
    if (isdata != "yes") AcceptanceCsc->Write();
    if (isdata != "yes") AcceptanceDt->Write();
    inputFile1->cd();
    inputFile1->Close();
    inputFile2->cd();
    inputFile2->Close();
    cout << "Closing output file." << endl;
    outputFile->Close();
    delete outputFile;
}
