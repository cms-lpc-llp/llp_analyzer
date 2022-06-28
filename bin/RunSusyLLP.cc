/*
 * Minimal example showing how to evaluate data contained in a flat tree structure using TensorFlow.
 * By default, the inference code uses single threading and no batching. The thread model is
 * configurable, however, note that this is done differently depending on the version of TensorFlow,
 * which changed significantly as of version 2.
 *
 * Author: Marcel Rieger
 */

// Adapted to LLP displaced jets in calorimeter by Lisa Benato


#include "assert.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include <TString.h>
#include <TChain.h>
#include "TH1F.h" 
#include "TH1D.h" 
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "LLPAnalysis/llpAnalyzer/interface/SusyLLP.h"

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

void usage()
{
  std::cerr << "Usage: RunSusyLLP  <input list>  [options]\n[options]:\n"
	    << "-d  --isData\n"
	    << "-f=  --outputFile=<output filename> (optional)\n"
	    << "-n=  --optionNumber=<option number> (optional)\n"
	    << "-l=  --optionLabel=<option Label> (optional)\n" 
	    << "-p=  --optionProcess=<option Process> (optional)\n" 
	    << "-s=  --optionSmearTag=<option smear Tag> (optional)\n" 
	    << "-h  --help"
	    << std::endl;
};




int main(int argc, char* argv[]){

  //get input files and analysis type from command line
  if ( ParseCommandLine( argc, argv, "--help" ) != ""  || ParseCommandLine( argc, argv, "-h" ) != ""  || argc < 2 )
    {
      usage();
      return -1;
    }
  
  //----------------------------------------
  //Getting <input list> and <analysis type>
  //----------------------------------------
  string inputFileName(argv[1]);
  
  //--------------------------------
  //G e t t i n g   d a t a  f l a g 
  //--------------------------------
  std::string _isData = ParseCommandLine( argc, argv, "--isData" );
  std::string _d = ParseCommandLine( argc, argv, "-d" );
  bool isData = false;
  if ( _isData == "yes" || _d == "yes" ) isData = true;

  //---------------------------------------------
  //G e t t i n g   o u t p u t F i l e   N a m e  
  //---------------------------------------------
  std::string _outFile = ParseCommandLine( argc, argv, "--outputFile=" );
  std::string _f = ParseCommandLine( argc, argv, "-f=" );
  string outputFileName = "";
  if ( _outFile != "" )
    {
      outputFileName = _outFile;
    }
  else if ( _f != "" )
    {
      outputFileName = _f;
    }
  else
    {
      std::cerr << "[WARNING]: output ROOT file not provided, using default output" << std::endl;
    }
  
  //-----------------------------------------
  //G e t t i n g   o p t i o n   n u m b e r
  //-----------------------------------------
  int option = -1;
  std::string _optionNumber = ParseCommandLine( argc, argv, "--optionNumber=" );
  std::string _n = ParseCommandLine( argc, argv, "-n=" );
  if ( _optionNumber != "" )
    {
      option = atoi( _optionNumber.c_str() );
    }
  else if ( _n != "" )
    {
      option = atoi( _n.c_str() );
    } 
  else
    {
      std::cerr << "[WARNING]: option number not provided, using default option number" << std::endl;
    }
  
  string label = "";
  std::string _optionLabel = ParseCommandLine( argc, argv, "--optionLabel=" );
  std::string _l = ParseCommandLine( argc, argv, "-l=" );
  if ( _optionLabel != "" ) 
    {
      label = _optionLabel;
    }
  else if ( _l != "" )
    {
      label = _l;
    }
  else
    {
      std::cerr << "[WARNING]: optional label not provided, using default optional label" << std::endl;
    }
  
  string process = "";
  std::string _optionProcess = ParseCommandLine( argc, argv, "--optionProcess=" );
  std::string _p = ParseCommandLine( argc, argv, "-p=" );
  if ( _optionProcess != "" ) 
    {
      process = _optionProcess;
    }
  else if ( _p != "" )
    {
      process = _p;
    }
  else
    {
      std::cerr << "[WARNING]: optional process not provided, using default optional process" << std::endl;
    }
  
  string smearTag = "";
  std::string _optionSmearTag = ParseCommandLine( argc, argv, "--optionSmearTag=" );
  std::string _s = ParseCommandLine( argc, argv, "-s=" );
  if ( _optionSmearTag != "" ) 
    {
      smearTag = _optionSmearTag;
    }
  else if ( _s != "" )
    {
      smearTag = _s;
    }
  else
    {
      std::cerr << "[WARNING]: optional smearTag not provided, using default optional smearTag" << std::endl;
    }
  
  
  std::cout << "[INFO]: <input list> --> " << inputFileName << std::endl;
  std::cout << "[INFO]: isData --> " << isData << std::endl;
  std::cout << "[INFO]: outputFileName --> " << outputFileName << std::endl;
  std::cout << "[INFO]: option --> " << option << std::endl;
  std::cout << "[INFO]: optionalLabel --> " << label << std::endl;
  std::cout << "[INFO]: optionalProcess --> " << process << std::endl;
  std::cout << "[INFO]: optionalSmearTag --> " << smearTag << std::endl;
    
    //build the TChain
    //tree name is set give the structure in the first root file, see while loop below
    TChain* theChain = new TChain();
    string curFileName;
    ifstream inputFile(inputFileName.c_str());
    int NFilesLoaded = 0;
    if ( !inputFile ){
      cerr << "Error: input file not found!" << endl;
      return -1;
    }
    
    while ( getline(inputFile, curFileName) )
      {
	if ( NFilesLoaded == 0 )
	  {
	    /*
	      checks root file structure and add first file
	    */
	    TFile* f_0 = TFile::Open( curFileName.c_str() );
	    if( f_0->GetDirectory("ntuples") )
	      {
		theChain->SetName("ntuples/llp");
		std::cout << "[INFO]: default configuration for tchain" << std::endl;
	      }
	    else
	      {
		theChain->SetName("llp");
		std::cout << "[INFO]: alternative configuration for tchain"<< std::endl;
	      }
	    theChain->Add( curFileName.c_str() );
	    delete f_0;
	  }
	else
	  {
	    //Addind remaining files after file structure is decided
	    theChain->Add( curFileName.c_str() );
	  }
        NFilesLoaded++;
      }     
    std::cout << "Loaded Total of " << NFilesLoaded << " files\n";
    if ( theChain == NULL ) return -1;
    

    SusyLLP analyzer(theChain);
    
    //------ EXECUTE ------//
    cout << "Executing SusyLLP..." << endl;
    analyzer.EnableAll();
    //analyzer.Analyze(isData, option, outputFileName, label, process);
    analyzer.Analyze(isData, option, outputFileName, label, process, smearTag);
    cout << "Process completed!" << endl;
    cerr << "------------------------------" << endl; //present so that an empty .err file corresponds to a failed job
   
    return 0;
}




 
