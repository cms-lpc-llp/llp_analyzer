//C++ INCLUDES
#include <iostream>

//ROOT  INCLUDES
#include <TFile.h>
#include <TTree.h>

//LOCAL INCLUDES
#include "SusyHggTree.hh"
#include "HighRes.hh"
#include "CommandLineInput.hh"

#define _debug  1

int main( int argc, char* argv[] ) 
{
  if ( argc < 3 )
    {
      std::cerr << "[ERROR]: Not Enough Arguments" << std::endl;
      return -1;
    }
  
  // G e t t i n g   F i l e   N a m e s   f r o m   C om m a n d   L i n e 
  //---------------------------------------------------------------------
  std::string fname1 = ParseCommandLine( argc, argv, "-file1=" );
  if ( _debug ) std::cout << "[DEBUG]: file name: " << fname1 << std::endl;
  std::string fname2 = ParseCommandLine( argc, argv, "-file2=" );
  if ( _debug ) std::cout << "[DEBUG]: file name: " << fname2 << std::endl;

  if ( fname1 == "" )
    {
      std::cerr << "[ERROR]: file name not provided" << std::endl;
      return -1;
    }
  
  if ( fname1.find( ".root" ) == std::string::npos )
    {
      std::cerr << "[ERROR]: please provide a valid root file" << std::endl;
      return -1;
    }
  
  // G e t t i n g  type f r o m   C o m m a n d   L i n e
  //------------------------------------------------------
  std::string type = ParseCommandLine( argc, argv, "-type=" );
  if ( _debug ) std::cout << "[INFO]: type: " << type << std::endl;
  
  if ( type != "HighRes" && type != "SusyHggTree" && type != "both" && type != "Both" )
    {
      std::cerr << "[ERROR]: type provided not recognized" << std::endl;
      return -1;
    }
  
  if ( type == "both" || type == "Both" )
    {
      if ( fname2 == "" )
	{
	  std::cerr << "[ERROR]: please provide HighRes ROOT file" << std::endl;
	  return -1;
	}
    }
  // C r e a t i n g   l i s t   o f   s e l e c t e d   e v e n t s
  //----------------------------------------------------------------
  TFile* f;
  TTree* T;
  if ( type == "SusyHggTree" || type == "both" || type == "Both" )
    {
      // R e t r i e v i n g    TTree
      //------------------------------- 
      f = new TFile( fname1.c_str() );
      T = (TTree*)f->Get("SusyHggTree");
      
      SusyHggTree* sht = new SusyHggTree( T );
      sht->Loop();
      sht->Loop2();
      delete f;
    }
  
  if ( type == "HighRes" || type == "both" || type == "Both" )
    {
      // R e t r i e v i n g    TTree                                                                          
      //-------------------------------                                                              
      if ( type != "HighRes" )
	{
	  f = new TFile( fname2.c_str() );
	}
      else
	{
	  f = new TFile( fname1.c_str() );
	}
      T = (TTree*)f->Get("HighRes");
      
      HighRes* hr = new HighRes( T );
      hr->Loop();
      delete f;
    }
  
  return 0;	
}
