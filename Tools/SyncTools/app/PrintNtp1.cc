//C++ INCLUDE
#include <iostream>
//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
//LOCAL INCLUDES
#include "ntp1.hh"

int main( int argc, char* argv[] )
{
  
  TFile* f = new TFile( argv[1] );
  TTree* t1 = (TTree*)f->Get("ntp1");
  ntp1* tree = new ntp1( t1 );
  tree->Loop();

  return 0;
};
