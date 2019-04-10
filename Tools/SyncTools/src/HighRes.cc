#define HighRes_cxx
#include "HighRes.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <fstream>

void HighRes::Loop()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  std::ofstream ofs ( "HighResCP.txt", std::fstream::out );
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    // A p p l y   s e l e c t i o n   c u t s 
    //----------------------------------------
    if ( MR>250. && t1Rsq>0.05 && pho1Pt>25. && pho2Pt>25. && (pho1Pt>40. || pho2Pt>40.)
	 && pTGammaGamma>20.
	 && fabs(Pho1Eta)<1.44 && fabs(Pho2Eta)<1.44 && mGammaGamma>103. && mGammaGamma<160. )
      {
	ofs << run << " " << event << "\n";
      }
  }
  ofs.close();
};
