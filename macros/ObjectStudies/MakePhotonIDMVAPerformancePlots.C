//================================================================================================
//
// HWW selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TGraphAsymmErrors.h>     
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
//#include <MyStyle.h>
#include "TLegend.h"
#include "TEfficiency.h"

#include "RazorAnalyzer/include/PhotonTree.h"

// lumi section selection with JSON files
// #include "MitCommon/DataFormats/interface/Types.h"
// #include "MitAna/DataCont/interface/RunLumiRangeMap.h"
// #include "MitCommon/MathTools/interface/MathUtils.h"
// #include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
// #include "MitHiggs/Utils/interface/EfficiencyUtils.h"
// #include "MitHiggs/Utils/interface/PlotUtils.h"
// #include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"
// #include "HiggsAna/Ntupler/interface/HiggsAnaDefs.hh"

#endif


//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("eleIDdir");
    if (!dir) {
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}



//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, string name ) {
  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double SigEff[nPoints];
  double BkgEff[nPoints];
  double SigEffErrLow[nPoints];
  double SigEffErrHigh[nPoints];
  double BkgEffErrLow[nPoints];
  double BkgEffErrHigh[nPoints];
  double NSigTotal = 0;
  double NBkgTotal = 0;
  
  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
    NBkgTotal += bkgHist->GetBinContent(q);
  }

  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    Double_t nbkg = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
      nbkg += bkgHist->GetBinContent(q);
    }

    Double_t ratio;
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    ratio = n1/n2;


    SigEff[b] = ratio;
    SigEffErrLow[b] = 0;
    SigEffErrHigh[b] = 0;

    n1 = TMath::Nint(nbkg);
    n2 = TMath::Nint(NBkgTotal);
    ratio = n1/n2;
    BkgEff[b] = ratio;
    BkgEffErrLow[b] = 0;
    BkgEffErrHigh[b] = 0;
  }

  TGraphAsymmErrors *tmpSigEffVsBkgEff = new TGraphAsymmErrors (nPoints, BkgEff, SigEff, BkgEffErrLow, BkgEffErrHigh, SigEffErrLow, SigEffErrHigh );
  tmpSigEffVsBkgEff->SetName(name.c_str());
  tmpSigEffVsBkgEff->SetTitle("");
  tmpSigEffVsBkgEff->GetXaxis()->SetTitle("Bkg Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitle("Signal Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsBkgEff->GetXaxis()->SetTitleOffset(1.05);
  tmpSigEffVsBkgEff->SetMarkerSize(0.5);
  tmpSigEffVsBkgEff->SetMarkerStyle(20);

  return tmpSigEffVsBkgEff;
}

//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeCurrentWPSigEffVsBkgEffGraph(Double_t signalEff, Double_t bkgEff, string name ) {
  //Make Met Plots
  double SigEff[1];
  double BkgEff[1];
  double SigEffErrLow[1];
  double SigEffErrHigh[1];
  double BkgEffErrLow[1];
  double BkgEffErrHigh[1];
  double NSigTotal = 0;
  double NBkgTotal = 0;
  double cutValue;

  SigEff[0] = signalEff;
  SigEffErrLow[0] = 0;
  SigEffErrHigh[0] = 0;
  BkgEff[0] = bkgEff;
  BkgEffErrLow[0] = 0;
  BkgEffErrHigh[0] = 0;

  TGraphAsymmErrors *tmpSigEffVsBkgEff = new TGraphAsymmErrors (1, BkgEff, SigEff, BkgEffErrLow, BkgEffErrHigh , SigEffErrLow, SigEffErrHigh );
  tmpSigEffVsBkgEff->SetName(name.c_str());
  tmpSigEffVsBkgEff->SetTitle("");
  tmpSigEffVsBkgEff->GetXaxis()->SetTitle("Bkg Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitle("Signal Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsBkgEff->GetXaxis()->SetTitleOffset(1.05);
  tmpSigEffVsBkgEff->SetMarkerColor(kBlack);
  tmpSigEffVsBkgEff->SetLineColor(kBlack);
  tmpSigEffVsBkgEff->SetMarkerSize(1.5);

  return tmpSigEffVsBkgEff;
}


//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeSigEffVsCutValueGraph(TH1F* signalHist, string name ) {

  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double cutValue[nPoints];
  double cutValueErr[nPoints];
  double SigEff[nPoints];
  double SigEffErrLow[nPoints];
  double SigEffErrHigh[nPoints];
  double NSigTotal = 0;
  
  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }

  for(UInt_t b=0; b < nPoints; ++b) {
    cutValue[b] = signalHist->GetXaxis()->GetBinCenter(b);
    cutValueErr[b] = 0;
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t ratio;
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    ratio = n1/n2;
    SigEff[b] = ratio;
    SigEffErrLow[b] = 0;
    SigEffErrHigh[b] = 0;

  }

  TGraphAsymmErrors *tmpSigEffVsCut = new TGraphAsymmErrors (nPoints, cutValue, SigEff, cutValueErr, cutValueErr, SigEffErrLow, SigEffErrHigh  );
  tmpSigEffVsCut->SetName(name.c_str());
  tmpSigEffVsCut->SetTitle("");
  tmpSigEffVsCut->GetXaxis()->SetTitle("Cut Value");
  tmpSigEffVsCut->GetYaxis()->SetTitle("Efficiency");
  tmpSigEffVsCut->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsCut->GetXaxis()->SetTitleOffset(1.05);

  return tmpSigEffVsCut;
}


//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeCurrentWPSigEffVsCutValueGraph(TH1F* signalHist, string name, Double_t myCutValue ) {
  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double cutValue[1] = {0};
  double cutValueErr[1] = {0};
  double SigEff[1] = {0};
  double SigEffErrLow[1] = {0};
  double SigEffErrHigh[1] = {0};
  double NSigTotal = 0;
  
  Double_t effDiff = 9999;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }

  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t ratio;
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    ratio = n1/n2;
    
      cout << myCutValue << " : " << signalHist->GetXaxis()->GetBinCenter(b) << " , " << cutValue[0] << endl;
    if (fabs(myCutValue - signalHist->GetXaxis()->GetBinCenter(b)) < fabs(myCutValue - cutValue[0])) {
      SigEff[0] = ratio;
      SigEffErrLow[0] = 0;
      SigEffErrHigh[0] = 0;
      cutValue[0] = signalHist->GetXaxis()->GetBinCenter(b);
      cutValueErr[0] = 0;
    }
  }

//   cout << "Final: " << cutValue[0] << " , " << SigEff[0] << endl;

  TGraphAsymmErrors *tmpSigEffVsCut = new TGraphAsymmErrors (1, cutValue, SigEff, cutValueErr, cutValueErr, SigEffErrLow, SigEffErrHigh  );
  tmpSigEffVsCut->SetName(name.c_str());
  tmpSigEffVsCut->SetTitle("");
  tmpSigEffVsCut->GetXaxis()->SetTitle("Cut Value");
  tmpSigEffVsCut->GetYaxis()->SetTitle("Efficiency");
  tmpSigEffVsCut->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsCut->GetXaxis()->SetTitleOffset(1.05);
  tmpSigEffVsCut->SetMarkerColor(kBlack);
  tmpSigEffVsCut->SetLineColor(kBlack);
  tmpSigEffVsCut->SetMarkerSize(1.5);

  return tmpSigEffVsCut;
}



//*************************************************************************************************
//
//*************************************************************************************************
Double_t FindCutValueAtFixedEfficiency(TH1F* signalHist, Double_t targetSignalEff ) {
  //Make Met Plots


  Double_t targetCutValue = -9999;
  Double_t bestCurrentSignalEff = 0;
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double NSigTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t ratio = nsig / NSigTotal;
//     cout << targetSignalEff << " : " << ratio << " , " << signalHist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetSignalEff - ratio) < fabs(targetSignalEff - bestCurrentSignalEff)) {
      targetCutValue = signalHist->GetXaxis()->GetBinCenter(b);
      bestCurrentSignalEff = ratio;
    }
  }

  return targetCutValue;
}


//*************************************************************************************************
//
//*************************************************************************************************
Double_t FindBkgEffAtFixedSignalEfficiency(TH1F* signalHist, TH1F* bkgHist, Double_t targetSignalEff ) {
  //Make Met Plots


  Double_t targetBkgEff = 0;
  Double_t bestCurrentSignalEff = 0;
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double NSigTotal = 0;
  double NBkgTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
    NBkgTotal += bkgHist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t nbkg = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nbkg += bkgHist->GetBinContent(q);
    }

    Double_t ratio = nsig / NSigTotal;
    Double_t bkgEff = nbkg / NBkgTotal;
//     cout << targetSignalEff << " : " << ratio << " , " << bkgEff << " : " << signalHist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetSignalEff - ratio) < fabs(targetSignalEff - bestCurrentSignalEff)) {
      bestCurrentSignalEff = ratio;
      targetBkgEff = bkgEff;
    }
  }

  return targetBkgEff;
}


//*************************************************************************************************
//
//*************************************************************************************************
Double_t FindSigEffAtFixedBkgEfficiency(TH1F* signalHist, TH1F* bkgHist, Double_t targetBkgEff ) {
  //Make Met Plots

  Double_t targetSignalEff = 0;
  Double_t bestCurrentBkgEff = 0;
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double NSigTotal = 0;
  double NBkgTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
    NBkgTotal += bkgHist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t nbkg = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nbkg += bkgHist->GetBinContent(q);
    }

    Double_t sigEff = nsig / NSigTotal;
    Double_t bkgEff = nbkg / NBkgTotal;
//     cout << targetSignalEff << " : " << ratio << " , " << bkgEff << " : " << signalHist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetBkgEff - bkgEff) < fabs(targetBkgEff - bestCurrentBkgEff)) {
      bestCurrentBkgEff = bkgEff;
      targetSignalEff = sigEff;
    }
  }

  return targetSignalEff;
}


Bool_t passPHYS14Tight( PhotonTree *phoTree) {

    return phoTree->fPhoIsTight;
}

Bool_t passPHYS14Loose( PhotonTree *phoTree) {

    return phoTree->fPhoIsLoose;
}



Bool_t passIDMVA( PhotonTree *phoTree) {

  Int_t subdet = 0;  
  if (fabs(phoTree->fPhoSCEta) < 1.479) subdet = 0;
  else subdet = 1;

  Int_t MVABin = -1;
  if (subdet == 0 ) MVABin = 0;
  if (subdet == 1 ) MVABin = 1;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = 0.0;
  if (MVABin == 1) MVACut = 0.0; 

  bool pass = false;
  if (phoTree->fPhoIDMVA > MVACut) {
    pass = true;
  }   
  return pass;
}



//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void PlotPhotonIDMVAPerformancePlots(string InputFile, string Label, Int_t Option)
{  

  string label = "";
  if (Label != "") label = "_" + Label;


  //*****************************************************************************************
  //Plotting Setup
  //*****************************************************************************************
//   vector<Int_t> markers;
//   vector<Int_t> colors;
//   colors.push_back(kRed);     markers.push_back(20);
//   colors.push_back(kCyan);    markers.push_back(21);
// //   colors.push_back(kBlue);    markers.push_back(21);
//   colors.push_back(kMagenta); markers.push_back(22);
//   colors.push_back(kCyan);    markers.push_back(34);
//   colors.push_back(kBlack);   markers.push_back(29);
//   colors.push_back(kGreen);   markers.push_back(33);
//   colors.push_back(kRed-2);   markers.push_back(33);
//   colors.push_back(kOrange);   markers.push_back(33);
//   colors.push_back(kBlue-2);   markers.push_back(33);
//   colors.push_back(kGreen-2);   markers.push_back(33);
//   colors.push_back(kMagenta-2);   markers.push_back(33);


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *PhoIDMVA_Real = new TH1F(("PhoIDMVA_Real"+label).c_str(), "; IDMVA ; Number of Events ",  10000, -2 , 2);
  TH1F *PhoIDMVA_Fake = new TH1F(("PhoIDMVA_Fake"+label).c_str(), "; IDMVA ; Number of Events ",  10000, -2 , 2);


  Double_t RealPhotons = 0;
  Double_t FakePhotons = 0;
  Double_t RealPhotonPassPHYS14Tight = 0;
  Double_t FakePhotonPassPHYS14Tight = 0;
  Double_t RealPhotonPassPHYS14Loose = 0;
  Double_t FakePhotonPassPHYS14Loose = 0;

  Double_t RealPhotonPassIDMVA = 0;
  Double_t FakePhotonPassIDMVA = 0;

  //*****************************************************************************************
  //PhoTree
  //*****************************************************************************************
  PhotonTree *PhoTree = new PhotonTree;
  PhoTree->LoadTree(InputFile.c_str());
  PhoTree->InitTree(PhotonTree::kPhotonTreeLight);

  cout << "Total Entries: " << PhoTree->tree_->GetEntries() << "\n";
  int nentries = PhoTree->tree_->GetEntries();
  nentries = 5000000;
  for(UInt_t ientry=0; ientry < PhoTree->tree_->GetEntries(); ientry++) {       	
    PhoTree->tree_->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
        
    bool isPrompt = false;
    bool isFake = false;
    if ( PhoTree->fMotherPdgId == 25 || (abs(PhoTree->fMotherPdgId) >= 1 && abs(PhoTree->fMotherPdgId) <= 6 ) 
	 || (abs(PhoTree->fMotherPdgId) >= 11 && abs(PhoTree->fMotherPdgId) <= 16 ) 
	 ) {
      isPrompt = true;
      if (RealPhotons >= nentries) continue;
    }
    if ( !( PhoTree->fMotherPdgId == 25 || (abs(PhoTree->fMotherPdgId) >= 1 && abs(PhoTree->fMotherPdgId) <= 6 ) 
	    || (abs(PhoTree->fMotherPdgId) >= 11 && abs(PhoTree->fMotherPdgId) <= 16 ) )
	 ) {
      isFake = true;
      if (FakePhotons >= nentries) continue;
    }
    

    if (PhoTree->fPhoPt < 20) continue;


    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(PhoTree->fPhoEta) < 1.485) subdet = 0;
    else subdet = 1;

    Int_t ptBin = 0;
    if (PhoTree->fPhoPt > 50.0) ptBin = 1;
    if (PhoTree->fPhoPt > 100.0) ptBin = 2;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 3) passCuts = (subdet == 1 && ptBin == 1);    
    if (Option == 4) passCuts = (subdet == 0 && ptBin == 2);
    if (Option == 5) passCuts = (subdet == 1 && ptBin == 2);    
    if (Option == 10) passCuts = (ptBin == 0 );

    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    


    //Some Preselection cuts
    //if (PhoTree->fDRToClosestParton < 1.0) continue;

    

    //Real Photon
    if (isPrompt) {
      RealPhotons += PhoTree->fWeight;
      if (passPHYS14Tight(PhoTree)) RealPhotonPassPHYS14Tight += PhoTree->fWeight;
      if (passPHYS14Loose(PhoTree)) RealPhotonPassPHYS14Loose += PhoTree->fWeight;
      if (passIDMVA(PhoTree)) RealPhotonPassIDMVA += PhoTree->fWeight;      

      PhoIDMVA_Real->Fill(PhoTree->fPhoIDMVA, PhoTree->fWeight);    
      //cout << "prompt : " << PhoTree->fPhoIDMVA << "\n";
    } 
    //FakePhoton
    else if ( isFake ) {
      FakePhotons += PhoTree->fWeight;
      if (passPHYS14Tight(PhoTree)) FakePhotonPassPHYS14Tight += PhoTree->fWeight;
      if (passPHYS14Loose(PhoTree)) FakePhotonPassPHYS14Loose += PhoTree->fWeight;
      if (passIDMVA(PhoTree)) FakePhotonPassIDMVA += PhoTree->fWeight;
      
      PhoIDMVA_Fake->Fill(PhoTree->fPhoIDMVA, PhoTree->fWeight);  
      //cout << "fake : " << PhoTree->fPhoIDMVA << "\n";
    }

  } 
  
  
  //*****************************************************************************************
  //Current Working Points
  //*****************************************************************************************
  cout << "PHYS14 Tight Real Photon Efficiency : " << RealPhotonPassPHYS14Tight << " / " << RealPhotons << " = " << RealPhotonPassPHYS14Tight/RealPhotons << endl;
  cout << "PHYS14 Tight Fake Photon Efficiency : " << FakePhotonPassPHYS14Tight << " / " << FakePhotons << " = " << FakePhotonPassPHYS14Tight/FakePhotons << endl;
  TGraphAsymmErrors* ROC_PHYS14TightWP = MakeCurrentWPSigEffVsBkgEffGraph(RealPhotonPassPHYS14Tight/RealPhotons , FakePhotonPassPHYS14Tight/FakePhotons, "ROC_PHYS14TightWP"+label);

  cout << "PHYS14 Loose Real Photon Efficiency : " << RealPhotonPassPHYS14Loose << " / " << RealPhotons << " = " << RealPhotonPassPHYS14Loose/RealPhotons << endl;
  cout << "PHYS14 Loose Fake Photon Efficiency : " << FakePhotonPassPHYS14Loose << " / " << FakePhotons << " = " << FakePhotonPassPHYS14Loose/FakePhotons << endl;
  TGraphAsymmErrors* ROC_PHYS14LooseWP = MakeCurrentWPSigEffVsBkgEffGraph(RealPhotonPassPHYS14Loose/RealPhotons , FakePhotonPassPHYS14Loose/FakePhotons, "ROC_PHYS14LooseWP"+label);



  Double_t BkgEffPHYS14Tight = FakePhotonPassPHYS14Tight/FakePhotons;
  Double_t SigEffPHYS14Tight = RealPhotonPassPHYS14Tight/RealPhotons;
  Double_t BkgEffPHYS14Loose = FakePhotonPassPHYS14Loose/FakePhotons;
  Double_t SigEffPHYS14Loose = RealPhotonPassPHYS14Loose/RealPhotons;

  cout << "**********************\n";
  Double_t SigEffIDMVA_AtTightBkgEff = FindSigEffAtFixedBkgEfficiency(PhoIDMVA_Real, PhoIDMVA_Fake, BkgEffPHYS14Tight);
  Double_t BkgEffIDMVA_AtTightSigEff = FindBkgEffAtFixedSignalEfficiency(PhoIDMVA_Real, PhoIDMVA_Fake, SigEffPHYS14Tight);
  cout << "Signal Efficiency (wrt PHYS14Tight Cut-based) for : same bkg \n";
  cout << "IDMVA : " << SigEffIDMVA_AtTightBkgEff/SigEffPHYS14Tight <<  endl;
  cout << "Bkg Efficiency (wrt PHYS14Tight Cut-based) for same sig eff \n";
  cout << "IDMVA : " << BkgEffIDMVA_AtTightSigEff/BkgEffPHYS14Tight << endl;
  cout << "**********************\n";

  cout << "**********************\n";
  Double_t BkgEffIDMVA_AtLooseSigEff = FindBkgEffAtFixedSignalEfficiency(PhoIDMVA_Real, PhoIDMVA_Fake, SigEffPHYS14Loose);
  Double_t SigEffIDMVA_AtLooseBkgEff = FindSigEffAtFixedBkgEfficiency(PhoIDMVA_Real, PhoIDMVA_Fake, BkgEffPHYS14Loose);
  cout << "Sig Efficiency (wrt PHYS14Loose Cut-based) for same bkg eff \n";
  cout << "IDMVA : " << SigEffIDMVA_AtLooseBkgEff/SigEffPHYS14Loose << endl;
  cout << "Bkg Efficiency (wrt PHYS14Loose Cut-based) for same sig eff \n";
  cout << "IDMVA : " << BkgEffIDMVA_AtLooseSigEff/BkgEffPHYS14Loose << endl;
  cout << "**********************\n";


  cout << "IDMVA  Real Photon Efficiency : " << RealPhotonPassIDMVA << " / " << RealPhotons << " = " << RealPhotonPassIDMVA/RealPhotons << endl;
  cout << "IDMVA  Fake Photon Efficiency : " << FakePhotonPassIDMVA << " / " << FakePhotons << " = " << FakePhotonPassIDMVA/FakePhotons << endl;
  TGraphAsymmErrors* ROC_IDMVAWP = MakeCurrentWPSigEffVsBkgEffGraph(RealPhotonPassIDMVA/RealPhotons , FakePhotonPassIDMVA/FakePhotons, "ROC_IDMVAWP"+label);



  //*****************************************************************************************
  //Make ROC curves
  //*****************************************************************************************
  TGraphAsymmErrors* ROC_IDMVA = MakeSigEffVsBkgEffGraph(PhoIDMVA_Real, PhoIDMVA_Fake, "ROC_IDMVA"+label );


  //*****************************************************************************************
  //Find Cut with same signal efficiency Make ROC curves
  //*****************************************************************************************
  Double_t CutValue_IDMVA_SameSig = FindCutValueAtFixedEfficiency(PhoIDMVA_Real, SigEffPHYS14Tight );
  Double_t CutValue_IDMVA_SameBkg = FindCutValueAtFixedEfficiency(PhoIDMVA_Fake, BkgEffPHYS14Tight );
  cout << "IDMVA Cut Value @ Same Cut-Based Tight Sig Eff: " << CutValue_IDMVA_SameSig << endl;
  cout << "IDMVA Cut Value @ Same Cut-Based Tight Bkg Eff: " << CutValue_IDMVA_SameBkg << endl;


  TLegend* legend;
  TCanvas* cv;
  string plotname;

  //*****************************************************************************************
  //Plot ROC Curves
  //*****************************************************************************************
  vector<TGraphAsymmErrors*> ROCGraphs;
  vector<string> GraphLabels;
  vector<Int_t> colors;

  //*****************************************************************************************
  //*****************************************************************************************
  ROCGraphs.clear();
  GraphLabels.clear();
  plotname = "PhotonIDMVA"+label;

  ROCGraphs.push_back(ROC_IDMVA);
  GraphLabels.push_back("IDMVA");
  colors.push_back(kGreen+2);
  

  //*****************************************************************************************
  Double_t xmin = 0.0;
  Double_t xmax = 1.0;
  Double_t ymin = 0.0;
  Double_t ymax = 1.0;


  cv = new TCanvas("cv", "cv", 800, 600);

//    legend = new TLegend(0.45,0.20,0.75,0.50);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(ROCGraphs[i],GraphLabels[i].c_str(), "LP");

    ROCGraphs[i]->SetMarkerColor(colors[i]);
    ROCGraphs[i]->SetLineColor(colors[i]);
    ROCGraphs[i]->SetMarkerSize(0.5);
   
    ROCGraphs[i]->GetXaxis()->SetRangeUser(xmin,xmax);    
    ROCGraphs[i]->GetYaxis()->SetRangeUser(ymin,ymax);    
    if (i==0) {
      ROCGraphs[i]->Draw("AP");
    } else {
      ROCGraphs[i]->Draw("Psame");
    }
  }

  legend->AddEntry(ROC_PHYS14TightWP, "PHYS14Tight WP", "P");
  ROC_PHYS14TightWP->SetFillColor(kRed);
  ROC_PHYS14TightWP->SetMarkerColor(kRed);
  ROC_PHYS14TightWP->SetMarkerStyle(34);
  ROC_PHYS14TightWP->SetMarkerSize(2.5);
  ROC_PHYS14TightWP->Draw("Psame");

  legend->AddEntry(ROC_PHYS14LooseWP, "PHYS14Loose WP", "P");
  ROC_PHYS14LooseWP->SetFillColor(kBlue);
  ROC_PHYS14LooseWP->SetMarkerColor(kBlue);
  ROC_PHYS14LooseWP->SetMarkerStyle(34);
  ROC_PHYS14LooseWP->SetMarkerSize(2.5);
  ROC_PHYS14LooseWP->Draw("Psame");

  // legend->AddEntry(ROC_IDMVAWP, "IDMVA WP", "P");
  // ROC_IDMVAWP->SetFillColor(kViolet);
  // ROC_IDMVAWP->SetMarkerColor(kViolet);
  // ROC_IDMVAWP->SetMarkerStyle(34);
  // ROC_IDMVAWP->SetMarkerSize(2.5);
  // ROC_IDMVAWP->Draw("Psame");
  

  legend->Draw();
  
  cv->SaveAs(("ROCGraphs_" + plotname + ".gif").c_str());

  gBenchmark->Show("WWTemplate");       
} 





void MakePhotonIDMVAPerformancePlots() {
 
  // PlotPhotonIDMVAPerformancePlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGJet_FakeGJet_Spring15_50ns.root","GJet_BarrelPt20To50",0);
  //PlotPhotonIDMVAPerformancePlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGJet_FakeGJet_Spring15_50ns.root","GJet_EndcapPt20To50",1);
 
  //PlotPhotonIDMVAPerformancePlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGJet_FakeGJet_Spring15_50ns.root","GJet_BarrelPt50To100",2);
  //PlotPhotonIDMVAPerformancePlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGJet_FakeGJet_Spring15_50ns.root","GJet_EndcapPt50To100",3);
 
  //PlotPhotonIDMVAPerformancePlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGJet_FakeGJet_Spring15_50ns.root","GJet_BarrelPt100",4);
  PlotPhotonIDMVAPerformancePlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/PhotonNtuple_PromptGJet_FakeGJet_Spring15_50ns.root","GJet_EndcapPt100",5);
 


}
