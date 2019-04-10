#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>
#include <string.h>
//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>

//#include <tree.hh>

const bool _debug = true;

//Axis
const float axisTitleSize = 0.06;
const float axisTitleOffset = .8;

const float axisTitleSizeRatioX   = 0.18;
const float axisLabelSizeRatioX   = 0.12;
const float axisTitleOffsetRatioX = 0.94;

const float axisTitleSizeRatioY   = 0.15;
const float axisLabelSizeRatioY   = 0.108;
const float axisTitleOffsetRatioY = 0.32;

//Margins
const float leftMargin   = 0.12;
const float rightMargin  = 0.05;
const float topMargin    = 0.07;
const float bottomMargin = 0.12;

//CMS STANDARD
TString CMSText = "CMS";
TString extraText   = "Preliminary";
//TString lumiText = "2.32 fb^{-1} (13 TeV)";
//TString lumiText = "35.9 fb^{-1} (13 TeV)";
TString lumiText = "Simulation (13 TeV)";

bool AddCMS( TCanvas* C );

using namespace std;

int main( int argc, char** argv )
{
  ifstream file;// read input list directly
  std::string line;
  int n = 0;
  std::cout << "[Usage]: ./PlotOnOytputs inputList \n" << std::endl;
  srand(time(NULL));
  gROOT->Reset();
  gStyle->SetOptStat(0);
  if ( argc == 2 )
  { 
        file.open (argv[1], ios::in | ios::binary);
        std::cout << "[INFO]: Opening file " << argv[1] << " ......" << std::endl;
        std::cout<< std::endl;
        if ( !file.is_open () )
        {
                std::cerr << "!! File open error:" << argv[1] << "; make sure the file is in the correct location" << std::endl;
                return 1;
        } else {
                while(getline(file,line)) ++n;
                std::cout << "n = " << n << "\n" << std::endl;
        }
        //std::string outputFile = argv[2];
  }

  std::ifstream ifs( argv[1], std::ifstream::in );
  assert(ifs);
  //while(std::getline(ifs,line)) n++;

  std::string rootFileName;
  TFile* fin; 
  TTree* tree; 

  while ( ifs.good() ){
          ifs >> rootFileName;
          if ( ifs.eof() ) continue;
          if ( rootFileName.find("#") != std::string::npos ) continue;
          if ( _debug ) std::cout << rootFileName << std::endl;
          fin = new TFile( rootFileName.c_str(), "UPDATE");
          //assert( fin );
          if ( _debug ) std::cout << "[INFO]: file: " << rootFileName << " passed check\n\n"<< std::endl;

          //------------------------
          //Getting TTree and Histos
          //------------------------
          tree = (TTree*)fin->Get("HggRazorLeptons");

  }

  long nentries = tree->GetEntries();

  int box = 10;
  float mGammaGamma = -99.;
  float pTGammaGamma = -99.;
  float MR = -99.;
  float t1Rsq = -99.;

  bool pho1passIso = false;
  bool pho2passIso = false;
  bool pho1passEleVeto = false;
  bool pho2passEleVeto = false;
  float pho1SC_Eta = -99.;
  float pho2SC_Eta = -99.;
  float pho1Pt = -99.;
  float pho2Pt = -99.;
  float pho1R9 = -99.;
  float pho2R9 = -99.;


  tree->SetBranchAddress("box", &box);
  tree->SetBranchAddress("mGammaGamma", &mGammaGamma);
  tree->SetBranchAddress("pTGammaGamma", &pTGammaGamma);
  tree->SetBranchAddress("MR", &MR);
  tree->SetBranchAddress("t1Rsq", &t1Rsq);
  
  tree->SetBranchAddress("pho1passIso", &pho1passIso);
  tree->SetBranchAddress("pho2passIso", &pho2passIso);
  tree->SetBranchAddress("pho1passEleVeto", &pho1passEleVeto);
  tree->SetBranchAddress("pho2passEleVeto", &pho2passEleVeto);
  tree->SetBranchAddress("pho1SC_Eta", &pho1SC_Eta);
  tree->SetBranchAddress("pho2SC_Eta", &pho2SC_Eta);
  tree->SetBranchAddress("pho1Pt", &pho1Pt);
  tree->SetBranchAddress("pho2Pt", &pho2Pt);
  tree->SetBranchAddress("pho1R9", &pho1R9);
  tree->SetBranchAddress("pho2R9", &pho2R9);
  
  //********************************************************
  //Print output
  //********************************************************
  //std::string outputFile = argv[2];
  //TFile* fout = new TFile( outputFile.c_str(), "RECREATE");
  TTree* razortree;
  //razortree = tree->CloneTree();

  int binN = -99;
  TString category = "None";

  TString cut = ""; 
 
  TString baseline_cut = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 ";
  
  //razortree = tree->CopyTree(cut);
  //razortree = tree->CopyTree(baseline_cut);
 
  //razortree->Branch("binN", &binN, "binN/I");          
  //razortree->Branch("category", &category); 
  
  //tree->SetBranchAddress("binN", &binN);
  tree->Branch("binN", &binN, "binN/I");          
  
  TString bin_cut[30];

  bin_cut[0]  =  "&& box==5 && MR>150. && MR<10000. && t1Rsq>0.260 "; 
  bin_cut[1]  =  "&& box==5 && MR>150. && MR<250.   && t1Rsq>0.170  && t1Rsq<0.260 "; 
  bin_cut[2]  =  "&& box==5 && MR>250. && MR<10000. && t1Rsq>0.170  && t1Rsq<0.260 "; 
  bin_cut[3]  =  "&& box==5 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.110 "; 
  bin_cut[4]  =  "&& box==5 && MR>150. && MR<350.   && t1Rsq>0.110  && t1Rsq<0.170 "; 
  bin_cut[5]  =  "&& box==5 && MR>350. && MR<10000. && t1Rsq>0.110  && t1Rsq<0.170 "; 

  bin_cut[6]  =  "&& box==6 && MR>150. && MR<10000. && t1Rsq>0.080 "; 
  bin_cut[7]  =  "&& box==6 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.080 "; 
 
  bin_cut[8]  =  "&& box==7 && MR>150. && MR<10000. && t1Rsq>0.090 "; 
  bin_cut[9]  =  "&& box==7 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.035 "; 
  bin_cut[10] =  "&& box==7 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.035 "; 
 
  bin_cut[11] =  "&& box==8 && MR>150. && MR<10000. && t1Rsq>0.325 "; 
  bin_cut[12] =  "&& box==8 && MR>150. && MR<10000. && t1Rsq>0.225  && t1Rsq<0.285 "; 
  bin_cut[13] =  "&& box==8 && MR>150. && MR<10000. && t1Rsq>0.285  && t1Rsq<0.325 "; 
  bin_cut[14] =  "&& box==8 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.185 "; 
  bin_cut[15] =  "&& box==8 && MR>150. && MR<200.   && t1Rsq>0.185  && t1Rsq<0.225 "; 
  bin_cut[16] =  "&& box==8 && MR>200. && MR<10000. && t1Rsq>0.185  && t1Rsq<0.225 "; 
 
  bin_cut[17] =  "&& box==9 && MR>150. && MR<10000. && t1Rsq>0.325 "; 
  bin_cut[18] =  "&& box==9 && MR>150. && MR<10000. && t1Rsq>0.225  && t1Rsq<0.285 "; 
  bin_cut[19] =  "&& box==9 && MR>150. && MR<10000. && t1Rsq>0.285  && t1Rsq<0.325 "; 
  bin_cut[20] =  "&& box==9 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.185 "; 
  bin_cut[21] =  "&& box==9 && MR>150. && MR<200.   && t1Rsq>0.185  && t1Rsq<0.225 "; 
  bin_cut[22] =  "&& box==9 && MR>200. && MR<10000. && t1Rsq>0.185  && t1Rsq<0.225 "; 
  
  bin_cut[23] =  "&& box==3 && MR>150. && MR<10000. && t1Rsq>0.                    && pTGammaGamma>=110. "; 
  bin_cut[24] =  "&& box==3 && MR>150. && MR<10000. && t1Rsq>0.                    && pTGammaGamma<110. "; 
  
  bin_cut[25] =  "&& box==3 && MR>150. && MR<10000. && t1Rsq>0.                    && pTGammaGamma>=110. "; 
  bin_cut[25] =  "&& box==3 && MR>150. && MR<10000. && t1Rsq>0.125                 && pTGammaGamma<110. "; 
  bin_cut[26] =  "&& box==3 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.055 && pTGammaGamma<110. "; 
  bin_cut[27] =  "&& box==3 && MR>150. && MR<10000. && t1Rsq>0.055  && t1Rsq<0.125 && pTGammaGamma<110. "; 
  
  bin_cut[29] =  "&& (box==0 || box==1 || box==2) && MR>0. && MR<10000. && t1Rsq>0. "; 
 
          for(int i=0; i<30; i++)
          {
                  cut = baseline_cut + bin_cut[i];

  //std::cout << "[INFO]: CUT--> = " << cut<< "\n" << std::endl;
  //std::cout << "if(" << cut << ") binN = " << i << ";"<< "\n" << std::endl;
          }


  for(int iEntry=0;iEntry<nentries;iEntry++){

          tree->GetEntry(iEntry);

          binN = -99;

/*
          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==5 && MR>150. && MR<10000. && t1Rsq>0.260 ) binN = 0;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==5 && MR>150. && MR<250.   && t1Rsq>0.170  && t1Rsq<0.260 ) binN = 1;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==5 && MR>250. && MR<10000. && t1Rsq>0.170  && t1Rsq<0.260 ) binN = 2;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==5 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.110 ) binN = 3;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==5 && MR>150. && MR<350.   && t1Rsq>0.110  && t1Rsq<0.170 ) binN = 4;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==5 && MR>350. && MR<10000. && t1Rsq>0.110  && t1Rsq<0.170 ) binN = 5;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==6 && MR>150. && MR<10000. && t1Rsq>0.080 ) binN = 6;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==6 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.080 ) binN = 7;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==7 && MR>150. && MR<10000. && t1Rsq>0.090 ) binN = 8;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==7 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.035 ) binN = 9;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==7 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.035 ) binN = 10;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==8 && MR>150. && MR<10000. && t1Rsq>0.325 ) binN = 11;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==8 && MR>150. && MR<10000. && t1Rsq>0.225  && t1Rsq<0.285 ) binN = 12;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==8 && MR>150. && MR<10000. && t1Rsq>0.285  && t1Rsq<0.325 ) binN = 13;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==8 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.185 ) binN = 14;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==8 && MR>150. && MR<200.   && t1Rsq>0.185  && t1Rsq<0.225 ) binN = 15;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==8 && MR>200. && MR<10000. && t1Rsq>0.185  && t1Rsq<0.225 ) binN = 16;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==9 && MR>150. && MR<10000. && t1Rsq>0.325 ) binN = 17;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==9 && MR>150. && MR<10000. && t1Rsq>0.225  && t1Rsq<0.285 ) binN = 18;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==9 && MR>150. && MR<10000. && t1Rsq>0.285  && t1Rsq<0.325 ) binN = 19;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==9 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.185 ) binN = 20;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==9 && MR>150. && MR<200.   && t1Rsq>0.185  && t1Rsq<0.225 ) binN = 21;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==9 && MR>200. && MR<10000. && t1Rsq>0.185  && t1Rsq<0.225 ) binN = 22;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==3 && MR>150. && MR<10000. && t1Rsq>0.                    && pTGammaGamma>=110. ) binN = 23;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==3 && MR>150. && MR<10000. && t1Rsq>0.                    && pTGammaGamma<110. ) binN = 24;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==3 && MR>150. && MR<10000. && t1Rsq>0.125                 && pTGammaGamma<110. ) binN = 25;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==3 && MR>150. && MR<10000. && t1Rsq>0.     && t1Rsq<0.055 && pTGammaGamma<110. ) binN = 26;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && box==3 && MR>150. && MR<10000. && t1Rsq>0.055  && t1Rsq<0.125 && pTGammaGamma<110. ) binN = 27;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 ) binN = 28;

          if(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 && (box==0 || box==1 || box==2) && MR>0. && MR<10000. && t1Rsq>0. ) binN = 29;
*/
          
          
          if(!(mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5 ))
          {
                  category = "None";
                  binN = -99;
          }
          else if(box==5)
          {
                  category = "highpt";
                  if(MR>150. && MR<10000. && t1Rsq>0.260 ) binN = 0;
                  else if(MR>150. && MR<250. && t1Rsq>0.170 && t1Rsq<0.260) binN = 1;
                  else if(MR>250. && MR<10000. && t1Rsq>0.170 && t1Rsq<0.260) binN = 2;
                  else if(MR>150. && MR<10000. && t1Rsq>0. && t1Rsq<0.110) binN = 3;
                  else if(MR>150. && MR<350. && t1Rsq>0.110 && t1Rsq<0.170) binN = 4;
                  else if(MR>350. && MR<10000. && t1Rsq>0.110 && t1Rsq<0.170) binN = 5;
          }
          else if(box==6)
          {
                  category = "hbb";
                  if(MR>150. && MR<10000. && t1Rsq>0.080 ) binN = 6;
                  else if(MR>150. && MR<10000. && t1Rsq>0. && t1Rsq<0.080) binN = 7;
          }
          else if(box==7)
          {
                  category = "zbb";
                  if(MR>150. && MR<10000. && t1Rsq>0.090 ) binN = 8;
                  else if(MR>150. && MR<10000. && t1Rsq>0. && t1Rsq<0.035) binN = 9;
                  else if(MR>150. && MR<10000. && t1Rsq>0.035 && t1Rsq<0.090) binN = 10;
          }
          else if(box==8)
          {
                  category = "highres";
                  if(MR>150. && MR<10000. && t1Rsq>0.325 ) binN = 11;
                  else if(MR>150. && MR<10000. && t1Rsq>0.225 && t1Rsq<0.285) binN = 12;
                  else if(MR>150. && MR<10000. && t1Rsq>0.285 && t1Rsq<0.325) binN = 13;
                  else if(MR>150. && MR<10000. && t1Rsq>0. && t1Rsq<0.185) binN = 14;
                  else if(MR>150. && MR<200. && t1Rsq>0.185 && t1Rsq<0.225) binN = 15;
                  else if(MR>200. && MR<10000. && t1Rsq>0.185 && t1Rsq<0.225) binN = 16;
          }
          else if(box==9)
          {
                  category = "lowres";
                  if(MR>150. && MR<10000. && t1Rsq>0.325 ) binN = 17;
                  else if(MR>150. && MR<10000. && t1Rsq>0.225 && t1Rsq<0.285) binN = 18;
                  else if(MR>150. && MR<10000. && t1Rsq>0.285 && t1Rsq<0.325) binN = 19;
                  else if(MR>150. && MR<10000. && t1Rsq>0. && t1Rsq<0.185) binN = 20;
                  else if(MR>150. && MR<200. && t1Rsq>0.185 && t1Rsq<0.225) binN = 21;
                  else if(MR>200. && MR<10000. && t1Rsq>0.185 && t1Rsq<0.225) binN = 22;
          }
          else if(box==3 && pTGammaGamma>=110.) 
          {
                  category = "muhighpt";
                  if(MR>150. && MR<10000. && t1Rsq>0.) binN  = 23; 
          }
          else if(box==3 && pTGammaGamma<110.)
          { 
                  category = "mulowpt";
                  if(MR>150. && MR<10000. && t1Rsq>0.) binN  = 24; 
          }
          else if(box==4 && pTGammaGamma>=110.) 
          {
                  category = "elehighpt";
                  if(MR>150. && MR<10000. && t1Rsq>0.) binN  = 25; 
          }
          else if(box==4 && pTGammaGamma<110.)
          { 
                  category = "elelowpt";
                  if(MR>150. && MR<10000. && t1Rsq>0.125) binN  = 26; 
                  else if(MR>150. && MR<10000. && t1Rsq>0. && t1Rsq<0.055) binN = 27;
                  else if(MR>150. && MR<10000. && t1Rsq>0.055 && t1Rsq<0.125) binN = 28;
          }
          else if(box==0||box==1||box==2)
          {
                  category = "twoleptons";
                  binN = 29;
          }

  //std::cout << "i = " << iEntry << " , category = " << category << " , biN = " << binN << "\n" << std::endl;

  //razortree->Fill();
  tree->Fill();

  }
  
          //std::cout << "[INFO]: CUT--> = " << cut<< "\n" << std::endl;

  //razortree->Write();
  tree->Write();

  int m = 2*n;
  std::cout << "m = " << m << "\n" << std::endl;

  


  //fout->Close();
  fin->Close();

}

bool AddCMS( TCanvas* C )
{
  C->cd();
  float lumix = 0.925;
  float lumiy = 0.945;
  float lumifont = 42;
  
  float cmsx = 0.225;
  float cmsy = 0.940;
  float cmsTextFont   = 61;  // default is helvetic-bold
  float extrax = cmsx + 0.198;
  float extray = cmsy;
  float extraTextFont = 52;  // default is helvetica-italics
  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;
  float cmsSize = 0.06;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    
  float extraTextSize = extraOverCmsTextSize*cmsSize;
  latex.SetTextFont(lumifont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(cmsSize);    
  latex.DrawLatex(lumix, lumiy,lumiText);

  latex.SetTextFont(cmsTextFont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(cmsSize);
  latex.DrawLatex(cmsx, cmsy, CMSText);
   
  latex.SetTextFont(extraTextFont);
  latex.SetTextAlign(31); 
  latex.SetTextSize(extraTextSize);
  latex.DrawLatex(extrax, extray, extraText);
  return true;
};
