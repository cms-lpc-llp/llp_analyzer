#define ntp1_cxx
#include "ntp1.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
//C++ INCLUDES
#include <fstream>
#include <cmath>
#include <math.h>

void ntp1::Loop()
{
  std::ofstream ofs ("test.txt", std::ofstream::out);
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //ofs << runNumber << " " << eventNumber << "\n";
    std::cout << "====== event " << jentry << " ======" << std::endl;
    int jet_ctr = 0;
    for( int i = 0; i < nAK5PFPUcorrJet; i++ )
      {
	double pt = sqrt( uncorrpxAK5PFPUcorrJet[i]*uncorrpxAK5PFPUcorrJet[i] + uncorrpyAK5PFPUcorrJet[i]*uncorrpyAK5PFPUcorrJet[i] );
	//PrintJetInfo( i );
	if ( pt > 20.0 )
	  {
	    std::cout << "jet " << jet_ctr;
	    PrintJetIdInfo( i );
	    PrintFractions( i );
	    jet_ctr++;
	  }
      }
  }
  ofs.close();
  
};

void ntp1::PrintJetInfo( int i )
{
  double pt = sqrt( uncorrpxAK5PFPUcorrJet[i]*uncorrpxAK5PFPUcorrJet[i] + uncorrpyAK5PFPUcorrJet[i]*uncorrpyAK5PFPUcorrJet[i] );
  float csv = combinedSecondaryVertexBJetTagsAK5PFPUcorrJet[i];
  if ( pt > 20.0 ) std::cout << "uE: " << uncorrenergyAK5PFPUcorrJet[i] << " pt: " << pt 
	    << " eta: " << etaAK5PFPUcorrJet[i] << " phi: " << phiAK5PFPUcorrJet[i]
	    << " jetArea: " << areaAK5PFPUcorrJet[i]  << " CSV: " << csv 
	    << " isLoose: " << isLoosePFPUcorrJet( i ) << " isMedium: " << isMediumPFPUcorrJet( i )
	    << " isTight: " << isTightPFPUcorrJet( i ) << std::endl;
};

void ntp1::PrintJetIdInfo( int i )
{
  double pt = sqrt( uncorrpxAK5PFPUcorrJet[i]*uncorrpxAK5PFPUcorrJet[i] + uncorrpyAK5PFPUcorrJet[i]*uncorrpyAK5PFPUcorrJet[i] );
  if ( pt > 20.0 ) std::cout << "uE: " << uncorrenergyAK5PFPUcorrJet[i] 
			     << " eta: " << etaAK5PFPUcorrJet[i] << " phi: " << phiAK5PFPUcorrJet[i] 
			     << " isLoose: " << isLoosePFPUcorrJet( i )
			     << " isTight: " << isTightPFPUcorrJet( i ) << std::endl;
};

bool ntp1::isLoosePFPUcorrJet(int i)
{
  double px = pxAK5PFPUcorrJet[i];
  double py = pyAK5PFPUcorrJet[i];
  double pt = sqrt( px*px + py*py );
  float UE = uncorrenergyAK5PFPUcorrJet[i];
  float eta = etaAK5PFPUcorrJet[i];
  
  float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[i]/UE;
  float neutralEMFrac = (photonEnergyAK5PFPUcorrJet[i]+HFEMEnergyAK5PFPUcorrJet[i])/UE;
  int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[i] + neutralHadronMultiplicityAK5PFPUcorrJet[i]+
    photonMultiplicityAK5PFPUcorrJet[i] + electronMultiplicityAK5PFPUcorrJet[i] +
    muonMultiplicityAK5PFPUcorrJet[i] + HFHadronMultiplicityAK5PFPUcorrJet[i] +
    HFEMMultiplicityAK5PFPUcorrJet[i];
  float muonFrac = muonEnergyAK5PFPUcorrJet[i]/UE;
  float chargedEMFrac = electronEnergyAK5PFPUcorrJet[i]/UE;
  float chargedHadFrac = chargedHadronEnergyAK5PFPUcorrJet[i]/UE;
  int chargedMult = chargedHadronMultiplicityAK5PFPUcorrJet[i]+
    electronMultiplicityAK5PFPUcorrJet[i]+
    muonMultiplicityAK5PFPUcorrJet[i]; 
  if ( neutralHadFrac < 0.99 && neutralEMFrac < 0.99 && nConstituents > 1 && muonFrac < 0.8  && chargedEMFrac < 0.9 )
    {
      if ( fabs(eta) > 2.4 )
        {
          return true;
        }
      else if ( fabs(eta) < 2.4 && chargedHadFrac > .0 && chargedMult > 0 && chargedEMFrac < 0.99)
        {
          return true;
        }
      else
        {
          return false;
        }
    }
  /*
    std::cout << "[DEBUG]: " << UE << " " << pt << " " << eta << " " << neutralHadFrac << " " << neutralEMFrac << " " 
  << nConstituents << " " << muonFrac << " " << chargedEMFrac
          << " " << chargedHadFrac << " " << chargedMult << std::endl;
  */
  return false;
};

bool ntp1::isMediumPFPUcorrJet(int i)
{
  float UE = uncorrenergyAK5PFPUcorrJet[i];
  float eta = etaAK5PFPUcorrJet[i];
  
  float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[i]/UE;
  float neutralEMFrac = (photonEnergyAK5PFPUcorrJet[i]+HFEMEnergyAK5PFPUcorrJet[i])/UE;
  int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[i] + neutralHadronMultiplicityAK5PFPUcorrJet[i]+
    photonMultiplicityAK5PFPUcorrJet[i] + electronMultiplicityAK5PFPUcorrJet[i] +
    muonMultiplicityAK5PFPUcorrJet[i] + HFHadronMultiplicityAK5PFPUcorrJet[i] +
    HFEMMultiplicityAK5PFPUcorrJet[i];
  float muonFrac = muonEnergyAK5PFPUcorrJet[i]/UE;
  float chargedEMFrac = electronEnergyAK5PFPUcorrJet[i]/UE;
  float chargedHadFrac = chargedHadronEnergyAK5PFPUcorrJet[i]/UE;
  int chargedMult = chargedHadronMultiplicityAK5PFPUcorrJet[i]+
    electronMultiplicityAK5PFPUcorrJet[i]+
    muonMultiplicityAK5PFPUcorrJet[i];
  if ( neutralHadFrac < 0.95 && neutralEMFrac < 0.95 && nConstituents > 1 && muonFrac < 0.8  && chargedEMFrac < 0.9 )
    {
      if ( fabs(eta) > 2.4 )
        {
          return true;
        }
      else if ( fabs(eta) < 2.4 && chargedHadFrac > .0 && chargedMult > 0 && chargedEMFrac < 0.99)
        {
          return true;
        }
      else
        {
          return false;
	}
    }
  
  return false;
};

bool ntp1::isTightPFPUcorrJet(int i)
{
  float UE = uncorrenergyAK5PFPUcorrJet[i];
  float eta = etaAK5PFPUcorrJet[i];
  
  float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[i]/UE;
  float neutralEMFrac = (photonEnergyAK5PFPUcorrJet[i]+HFEMEnergyAK5PFPUcorrJet[i])/UE;
  int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[i] + neutralHadronMultiplicityAK5PFPUcorrJet[i]+
    photonMultiplicityAK5PFPUcorrJet[i] + electronMultiplicityAK5PFPUcorrJet[i] +
    muonMultiplicityAK5PFPUcorrJet[i] + HFHadronMultiplicityAK5PFPUcorrJet[i] +
    HFEMMultiplicityAK5PFPUcorrJet[i];
  float muonFrac = muonEnergyAK5PFPUcorrJet[i]/UE;
  float chargedEMFrac = electronEnergyAK5PFPUcorrJet[i]/UE;
  float chargedHadFrac = chargedHadronEnergyAK5PFPUcorrJet[i]/UE;
  int chargedMult = chargedHadronMultiplicityAK5PFPUcorrJet[i]+
    electronMultiplicityAK5PFPUcorrJet[i]+
    muonMultiplicityAK5PFPUcorrJet[i];
  if ( neutralHadFrac < 0.90 && neutralEMFrac < 0.90 && nConstituents > 1 && muonFrac < 0.8  && chargedEMFrac < 0.9 )
    {
      if ( fabs(eta) > 2.4 )
        {
          return true;
        }
      else if ( fabs(eta) < 2.4 && chargedHadFrac > .0 && chargedMult > 0 && chargedEMFrac < 0.99)
        {
          return true;
        }
      else
        {
          return false;
	}
    }

  return false;
};

void ntp1::PrintFractions( int i )
{
  float UE = uncorrenergyAK5PFPUcorrJet[i];
  float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[i]/UE;
  float neutralEMFrac = (photonEnergyAK5PFPUcorrJet[i]+HFEMEnergyAK5PFPUcorrJet[i])/UE;
  int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[i] + neutralHadronMultiplicityAK5PFPUcorrJet[i]+
    photonMultiplicityAK5PFPUcorrJet[i] + electronMultiplicityAK5PFPUcorrJet[i] +
    muonMultiplicityAK5PFPUcorrJet[i] + HFHadronMultiplicityAK5PFPUcorrJet[i] +
    HFEMMultiplicityAK5PFPUcorrJet[i];
  float muonFrac = muonEnergyAK5PFPUcorrJet[i]/UE;
  float chargedEMFrac = electronEnergyAK5PFPUcorrJet[i]/UE;
  float chargedHadFrac = chargedHadronEnergyAK5PFPUcorrJet[i]/UE;
  int chargedMult = chargedHadronMultiplicityAK5PFPUcorrJet[i]+
    electronMultiplicityAK5PFPUcorrJet[i]+
    muonMultiplicityAK5PFPUcorrJet[i];
  
  std::cout << "neutralHadFrac: " << neutralHadFrac << " neutralEMFrac: "<< neutralEMFrac
	    << " nConstituents: " << nConstituents << " muonFrac: " << muonFrac 
	    << " chargedEMFrac: " << chargedEMFrac << " chargedHadFrac: " << chargedHadFrac
	    << " chargedMult: " << chargedMult << std::endl;
};
