//C++ INCLUDE
#include <iostream>
//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
//LOCAL INCLUDES
#include "RazorEvents.hh"

int main( int argc, char* argv[] )
{
  
  TFile* f = new TFile( argv[1] );
  TDirectory* dir = (TDirectory*)f->Get("ntuples");
  TTree* t1 = (TTree*)dir->Get("RazorEvents");
  RazorEvents* tree = new RazorEvents( t1 );
  tree->Loop();

  return 0;
};
