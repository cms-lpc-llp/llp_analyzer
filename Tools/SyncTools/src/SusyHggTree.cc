#define SusyHggTree_cxx
#include "SusyHggTree.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include <string>
#include <sstream>
#include <map>



void SusyHggTree::Loop()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  std::ofstream ofs ( "SusyHgg_Inverted.txt", std::ifstream::out );
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if (
	MR>250. && t1Rsq>0.05 && pho1_pt>25. && pho2_pt>25. && (pho1_pt>40. || pho2_pt>40.)
	&& pho1_pass_iso==1 && pho2_pass_iso==1 && ptgg>20.
	&& fabs(pho1_eta)<1.44 && fabs(pho2_eta)<1.44 && mgg>103. && mgg<160.
	)
      {
	ofs << run << " " << evt << std::endl;
      }
  }
  ofs.close();
};


void SusyHggTree::Loop2()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  
  std::ifstream ifs ( "/Users/cmorgoth/Work/git/SyncTools/SusyHgg_Inverted.txt", std::fstream::in );
  std::string r_run, e_evt;
  std::map< std::string, EVT > mymap;
  if ( ifs.is_open() )
    {
      while ( ifs.good() )
        {
          ifs >> r_run >> e_evt;
          std::string tmp = r_run + e_evt;
          EVT tmp_evt;
          tmp_evt.run = r_run;
          tmp_evt.event = e_evt;
          if ( mymap.find( tmp ) == mymap.end() )
            {
              mymap[tmp] = tmp_evt;
            }
        }
    }
  else
    {
      std::cout << "[ERROR]: unable to open file" << std::endl;
    }
  
  std::cout << "[INFO]: map size: " << mymap.size() << std::endl;
  int iso_ctr = 0;
  int id_ctr = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) break;
   nb = fChain->GetEntry(jentry);   nbytes += nb;
   // if (Cut(ientry) < 0) continue;
   
   std::stringstream ss;
   ss << run << evt;
   if ( mymap.find( ss.str() ) == mymap.end() )continue;                                                                                       
   //if ( !( run == 206859 && event == 24345 ) ) continue;                                                                                       
    
   std::cout << "============" << std::endl;
   std::cout << "run == " << run << " && evt == " << evt << std::endl;

   if ( pho1_pass_iso == 0 || pho2_pass_iso == 0 )
     {
       std::cout << "[INFO]: failed ISO: pho1, pho2 -> " <<  pho1_pass_iso << ", " << pho2_pass_iso << std::endl;
       iso_ctr++;
     }
   
   if ( pho1_pass_id == 0 || pho2_pass_id == 0 )
     {
       std::cout << "[INFO]: failed ID: pho1, pho2 -> " <<  pho1_pass_id << ", " << pho2_pass_id << std::endl;
       std::cout << "pho1pt: " << pho1_pt << " pho1eta: " << pho1_eta << std::endl;
       std::cout << "pho2pt: " << pho2_pt << " pho2eta: " << pho2_eta << std::endl;
       passID( pho1_sc_eta, false, pho1_sieie, pho1_HE );
       passID( pho2_sc_eta, false, pho2_sieie, pho2_HE );
       id_ctr++;
     }

   if ( MR < 150. )
     {
       std::cout << "[INFO]: failed MR: MR -> " <<  MR << std::endl;
     }

   if ( mgg > 180. )
     {
       std::cout << "[INFO]: failed mgg: mgg -> " <<  mgg << std::endl;
     }

   if ( mgg > 180. )
     {
       std::cout << "[INFO]: failed mgg: mgg -> " <<  mgg << std::endl;
     }
   
   /*
     if ( pho1_pass_iso == 1 && pho2_pass_iso == 1 && ( pho1_pt > 40. || pho2_pt > 40. ) && pho1_pass_id == 1 && pho2_pass_id ==1
      && mgg > 100. && mgg < 180. && MR > 150. && pho1_eleveto == 0 && pho2_eleveto == 0  && pho1_pt > 24.0  && pho2_pt > 24.0 && ptgg > 20.0 )
      {
      ofs << run << " " << evt << std::endl;
      }
   */
  }

  std::cout << "======Summary=====" << std::endl;
  std::cout << "Failed Iso: " << iso_ctr << std::endl;
  std::cout << "Failed ID: " << id_ctr << std::endl;
};

bool  SusyHggTree::passID(float eta, bool csev, float sieie, float hoe)
{
  bool _isEB = false;
  if ( fabs( eta ) < 1.44 ) _isEB = true;

  if ( _isEB )
    {
      if ( hoe > 0.05 )
	{
	  std::cout << "EB, failed HOE: " << hoe << std::endl; 
	  return false;
	}
      
      if ( sieie > 0.012 )
	{
	  std::cout << "EB, failed sIeIe: " << sieie << std::endl; 
	  return false;
	}
    }
  else
    {
      if ( hoe > 0.05 )
	{
	  std::cout << "EE, failed HOE: " << hoe << std::endl; 
	  return false;
	}
      
      if ( sieie > 0.034 )
	{
	  std::cout << "EE, failed sIeIe: " << sieie << std::endl; 
	  return false;
	}
    }
			      
  return true;
};
