//C++ INCLUDES
#include <iostream>
//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
//LOCAL INCLUDES
#include "ntp1.hh"
#include "aux.hh"

#define _debug 1

int main( int argc, char* argv[] )
{

  if( ParseCommandLine( argc, argv, "-help" ) != "" )
    {
      std::cout << "======================[HELP]=======================" << std::endl;
      Usage();
      return 0;
    }
  
  if ( argc < 3 )
    {
      std::cerr << "[ERROR]: insufficient command line arguments provided" << std::endl;
      Usage();
      return 0;
    }
  
  std::string fname1 = ParseCommandLine( argc, argv, "-file1=" );
  std::string fname2 = ParseCommandLine( argc, argv, "-file2=" );

  if ( fname1 == "" || fname2 == "" )
    {
      std::cerr <<  "[ERROR]: please provide two files with run and events" << std::endl;
      return -1;
    }
  
  // D e f i n e   o u t p u t   f i l e   n a m e s
  //-----------------------------------------------
  std::string outname1;
  std::string outname2;
  //First output file
  if ( fname1.find( "." ) != std::string::npos )
    {
      outname1 = fname1.substr( 0, fname1.find_last_of( '.' ) ) + "_Missing.txt";
      if ( _debug ) std::cout << "[DEBUG]: output file name: " << outname1 << std::endl;
    }
  else
    {
      outname1 = fname1 + "_Missing.txt";
      if ( _debug ) std::cout << "[DEBUG]: output file name: " << outname1 << std::endl;
    }
  
  //First output file                                                                             
  if ( fname2.find( "." ) != std::string::npos )
    {
      outname2 = fname2.substr( 0, fname2.find_last_of( '.' ) ) + "_Missing.txt";
      if ( _debug ) std::cout << "[DEBUG]: output file name: " << outname2 << std::endl;
    }
  else
    {
      outname2 = fname2 + "_Missing.txt";
      if ( _debug ) std::cout << "[DEBUG]: output file name: " << outname2 << std::endl;
    }
  
  
  std::map < std::string, evt > map1 = MakeMap( fname1 );
  std::map < std::string, evt > map2 = MakeMap( fname2 );
  CreateDiff( map1, map2, outname2 );//writes in outname2, missing run,event map2 
  CreateDiff( map2, map1, outname1 );//writes in outname1, missing run,event map1 
  
  return 0;
};
