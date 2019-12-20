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
  std::cout << "[Usage]: ./PlotOnOytputs inputList output.root\n" << std::endl;
  srand(time(NULL));
  gROOT->Reset();
  gStyle->SetOptStat(0);
  std::string outputFile;
  if ( argc == 3 )
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
        outputFile = argv[2];
        std::cout << "[INFO]: Output file " << argv[2] << " ......" << std::endl;
  }

  std::ifstream ifs( argv[1], std::ifstream::in );
  assert(ifs);

  std::string rootFileName;
  TFile* fin; 
  TTree* tree;

  while ( ifs.good() ){
          ifs >> rootFileName;
          if ( ifs.eof() ) continue;
          if ( rootFileName.find("#") != std::string::npos ) continue;
          if ( _debug ) std::cout << rootFileName << std::endl;
          fin = new TFile( rootFileName.c_str());
          if ( _debug ) std::cout << "[INFO]: file: " << rootFileName << " passed check\n\n"<< std::endl;

          //------------------------
          //Getting TTree and Histos
          //------------------------
          tree = (TTree*)fin->Get("HggRazorLeptons");

  }

  long nentries = tree->GetEntries();
  bool HLTDecision[300];

  tree->Branch("HLTDecision", HLTDecision, "HLTDecision[300]/O");
  //tree->SetBranchAddress("HLTDecision", &HLTDecision);
  
  //********************************************************
  //Print output
  //********************************************************
  TFile* fout = new TFile( outputFile.c_str(), "RECREATE");
  TTree* razortree;
  //2016 data
  //razortree = (TTree*)tree->CopyTree("HLTDecision[82]");
  //razortree = (TTree*)tree->CopyTree("HLTDecision[82] && Flag_goodVertices==1 && Flag_CSCTightHaloFilter==1 && Flag_HBHENoiseFilter==1 && Flag_HBHEIsoNoiseFilter==1 && Flag_EcalDeadCellTriggerPrimitiveFilter==1 && Flag_badMuonFilter==1  && Flag_badChargedCandidateFilter==1 && Flag_eeBadScFilter==1 ");
  //2016 Mc
  //razortree = (TTree*)tree->CopyTree("HLTDecision[82] && Flag_goodVertices==1 && Flag_CSCTightHaloFilter==1 && Flag_HBHENoiseFilter==1 && Flag_HBHEIsoNoiseFilter==1 && Flag_EcalDeadCellTriggerPrimitiveFilter==1 && Flag_badMuonFilter==1  && Flag_badChargedCandidateFilter==1 ");
  //2017 data
  //razortree = (TTree*)tree->CopyTree("HLTDecision[54]");
  //razortree = (TTree*)tree->CopyTree("HLTDecision[54] && Flag_goodVertices==1 && Flag_CSCTightHaloFilter==1 && Flag_HBHENoiseFilter==1 && Flag_HBHEIsoNoiseFilter==1 && Flag_EcalDeadCellTriggerPrimitiveFilter==1 && Flag_BadPFMuonFilter==1  && Flag_BadChargedCandidateFilter==1 && Flag_eeBadScFilter==1 && Flag_ecalBadCalibFilter==1");
  //2017 Mc
  //razortree = (TTree*)tree->CopyTree("HLTDecision[54] && Flag_goodVertices==1 && Flag_CSCTightHaloFilter==1 && Flag_HBHENoiseFilter==1 && Flag_HBHEIsoNoiseFilter==1 && Flag_EcalDeadCellTriggerPrimitiveFilter==1 && Flag_badMuonFilter==1  && Flag_badChargedCandidateFilter==1 ");
  //V4p4
  razortree = (TTree*)tree->CopyTree("HLTDecision[54] && Flag_HBHENoiseFilter == 1 && Flag_CSCTightHaloFilter == 1 && Flag_goodVertices == 1 && Flag_HBHEIsoNoiseFilter == 1 && Flag_EcalDeadCellTriggerPrimitiveFilter==1 && Flag_BadPFMuonFilter==1  && Flag_BadChargedCandidateFilter==1 && Flag_ecalBadCalibFilter==1 ");
  razortree->Write();

  fout->Close();
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
