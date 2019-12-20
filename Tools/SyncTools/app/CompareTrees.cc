//C++ INCLUDES
#include <iostream>
//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
//LOCAL INCLUDES
#include "ntp1.hh"
#include "aux.hh"

int main( int argc, char* argv[] )
{
  TFile* f1 = new TFile("data/pick_events_small_Alex_WH_ZH_HtoG_M-125_Total.root");
  TTree* t1 = (TTree*)f1->Get("ntp1");
  ntp1* tree_1 = new ntp1( t1 );

  TFile* f2 = new TFile("data/pick_events_small_Alex_WH_ZH_HtoG_M-125_Total.root");
  TTree* t2 = (TTree*)f2->Get("ntp1");
  ntp1* tree_2 = new ntp1( t2 );
  
  //std::string root_file1 = "data/pick_events_small_Alex_WH_ZH_HtoG_M-125_Total.root";
  std::string root_file1 = "data/test.root";
  std::string root_file2 = "data/test_cp.root";
  //std::string root_file2 = "data/pick_events_small_Alex_WH_ZH_HtoG_M-125_Total.root";
  //std::string root_file1 = "data/pick_events_small_Alex_Photon_RunA_Total.root";                      
  //std::string root_file2 = "data/default_data_PhotonRunA.root";
  std::string tree_name1 = "ControlSampleEvent";
  std::string tree_name2 = "ControlSampleEvent";
  //std::string tree_name2 = "ntp1"; 
  std::string var_list   = "input/input_ControlSample_diffTree.txt";
  //std::string var_list   = "input/input_variables.txt";  
  //CreateDiff( root_file1, root_file2, tree_name1, tree_name2 );
  CreateDiff( root_file1, root_file2, tree_name1, tree_name2, var_list, false );
  //CreateDiff( root_file1, root_file2, tree_name1, tree_name2, var_list, true ); 
  return 0;
};
