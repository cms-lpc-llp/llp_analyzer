#define RazorEvents_cxx
#include "RazorEvents.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
//C++ INCLUDES
#include <iostream>

void RazorEvents::Loop()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    std::cout << "====== event " << jentry << " ======" << std::endl;
    for ( int i = 0; i < nJets; i++ )
      {
	std::cout << "jet " << i << "uE: " << jetE[i] 
		  << " eta: " << jetEta[i] << " phi: " << jetPhi[i] 
		  << " isLoose: " << jetPassIDLoose[i] 
		  << " isTight: " << jetPassIDTight[i] << std::endl;
      }
  }
}
