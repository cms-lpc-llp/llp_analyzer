//C++ INCLUDES
#include <iostream>
#include <map>
#include <string>
//ROOT INCLUDES
//LOCAL INCLUDES
#include "ReadEventlist.hh"

std::string ParseCommandLine( int argc, char* argv[], std::string opt );
void Usage();

int main( int argc, char* argv[] )
{
  //Asking for help
  if( ParseCommandLine( argc, argv, "-help" ) != "" )
    {
      std::cout << "======================[HELP]=======================" << std::endl;
      Usage();
      return 0;
    }
  
  if ( argc <=3 )
    {
      std::cerr << "[ERROR]: insufficient command line arguments provided" << std::endl;
      Usage();
      return 0;
    }
  std::map< std::string , RunAndEvent > mymap;
  //std::string fname( argv[1] );
  //std::string list_name( argv[2] );
  //std::string root_output( argv[3] );
  
  //Required input from command line
  std::string pick_events = "";
  std::string list_name   = "";
  std::string root_output = "";
  
  //Vecbos Default
  std::string tree_name    = "ntp1";
  std::string run_branch   = "runNumber";
  std::string event_branch = "eventNumber";
  
  pick_events = ParseCommandLine( argc, argv, "-event_list=" );
  list_name   = ParseCommandLine( argc, argv, "-list_of_ntuples=" );
  root_output = ParseCommandLine( argc, argv, "-output_name=" );
  if ( ParseCommandLine( argc, argv, "-tree_name=" ) != "" ) tree_name = ParseCommandLine( argc, argv, "-tree_name=" );
  if ( ParseCommandLine( argc, argv, "-run_branch=" ) != "" ) run_branch = ParseCommandLine( argc, argv, "-run_branch=" );
  if ( ParseCommandLine( argc, argv, "-event_branch=" ) != "" ) event_branch = ParseCommandLine( argc, argv, "-event_branch=" );
  
  std::cout << "[INFO]: pick events file = " << pick_events << std::endl;
  std::cout << "[INFO]: list of ntuple file = " << list_name << std::endl;
  std::cout << "[INFO]: output file name = " << root_output << std::endl;
  std::cout << "[INFO]: TTree name = " << tree_name << std::endl;
  std::cout << "[INFO]: run TBranch name = " << run_branch << std::endl;
  std::cout << "[INFO]: event TBranch name = " << event_branch << std::endl;
  
  if ( pick_events == "" || list_name == "" || root_output == "" )
    {
      std::cerr << "[ERROR]: insufficient/incorrect command line arguments provided" << std::endl;
      Usage();
      return 0;
    }
  
  if( tree_name == "" )
    {
      std::cerr << "[ERROR]: overwritting of <tree_name> failed" << std::endl;
      Usage();
      return 0;
    }
  
  if( run_branch == "" )
    {
      std::cerr << "[ERROR]: overwritting of <run_branch> failed" << std::endl;
      Usage();
      return 0;
    }
  
  if( event_branch == "" )
    {
      std::cerr << "[ERROR]: overwritting of <event_branch> failed" << std::endl;
      Usage();
      return 0;
    }
  
  //Initializing Map
  FillMap( mymap, pick_events );
  
  /*
    for ( auto& tmp : mymap )
    {
    std::cout << tmp.second.run << " "  << tmp.second.event << std::endl;
    }
  */
  //Skimming Tree
  bool _isRazorEvent = false;
  if ( tree_name == "RazorEvents" ) _isRazorEvent = true;
  SkimTree( mymap, list_name, tree_name, run_branch, event_branch, root_output, _isRazorEvent );
  return 0;
};

std::string ParseCommandLine( int argc, char* argv[], std::string opt )
{
  for (int i = 1; i < argc; i++ )
    {
      std::string tmp( argv[i] );
      if ( tmp.find( opt ) != std::string::npos )
	{
	  return tmp.substr( tmp.find_last_of("=") + 1 );
	}
    }

  return "";
};


void Usage()
{
  std::cout << "===================================================" << std::endl;
  std::cout << "[USAGE]: for general information, please see README" << std::endl;
  std::cout << "===================================================" << std::endl;
  
  std::cout << "[USAGE]: ./CreateSkimmed <options>" << std::endl;
  std::cout << "[USAGE]: available <options> are the following:" << std::endl;
  std::cout << "[USAGE]: --help: print this message" << std::endl;
  std::cout << "[USAGE]: --event_list=<run_event_file_name>;\nplease provide in <run_event_file_name> with file name containing the run and event you want to select" << std::endl;
  std::cout << "[USAGE]: --list_of_ntuples==<ntuple_list_file_name>;\nplease provide in <ntuple_list_file_name> with file name containing the original ntuple list" << std::endl;
  std::cout << "[USAGE]: --output_name==<output_root_file_name>;\nplease provide in <output_root_file_name> with the root file output name for the selected events" << std::endl;
  std::cout << "[USAGE]: --tree_name==<tree_name>;\nplease provide in <tree_name> the correct TTree name in the ntuple " << std::endl;
  std::cout << "[USAGE]: --run_branch==<run_branch_name>;\nplease provide in <run_branch_name> the correct TBranch name used for run" << std::endl;
  std::cout << "[USAGE]: --event_branch==<event_branch_name>;\nplease provide in <event_branch_name> the correct TBranch name used for event" << std::endl;
  
};
