#include "DummyAnalyzer.h"

void DummyAnalyzer::Analyze(bool isData, int option, string outputFileName, string label)
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //Dummy example: print out the MET and the number of jets
        cout << "MET = " << metPt << "; Number of jets = " << nJets << endl;
    }
}
