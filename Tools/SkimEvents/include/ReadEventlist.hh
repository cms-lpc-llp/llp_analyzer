#ifndef READ_EVENT_LIST_HH
#define READ_EVENT_LIST_HH 1

//C++ INCLUDES
#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
//ROOT INCLUDES
//LOCAL INCLUDES

struct RunAndEvent
{
  int run;
  unsigned long event;  
  
  /*
  bool operator==( const RunAndEvent& tmp ) 
  {
    if ( this->event == tmp.event && this->run == tmp.run ) return true;
    return false;
  };
  */
};


bool FillMap( std::map< std::string, RunAndEvent >& mymap, std::string fname );
void SkimTree( std::map< std::string, RunAndEvent > mymap, std::string list_name, 
	       std::string tree_name, std::string run_branch, std::string event_branch, std::string root_output,
	       bool Tdir = false );

#endif
