//root -l /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/macros/ObjectStudies/MakeMuonIsoPerformancePlots.C+'("/afs/cern.ch/user/s/sixie/work/public/Run2SUSY/MuonNtuple/MuonNtuple_TTJetsPrompt_T1bbbbFake.root","All",-1)'

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

#include "RazorAnalyzer/include/MuonTree.h"

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
TGraphAsymmErrors* MakeSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, string name, bool cutAbove = true ) {
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

    if (cutAbove) {
      for (UInt_t q=0; q < b ; ++q) {
	nsig += signalHist->GetBinContent(q);
	nbkg += bkgHist->GetBinContent(q);
      }
    } else {
      for (UInt_t q=b; q < nPoints+2; ++q) {
	nsig += signalHist->GetBinContent(q);
	nbkg += bkgHist->GetBinContent(q);
      }
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
    for (UInt_t q=0; q < b; ++q) {
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


Bool_t passCSA14Preselection( MuonTree *muTree) {

    bool pass = false;
 
    pass = true;

    return pass;
}



Bool_t passLooseID( MuonTree *muTree) {
    bool pass = false;
    if (muTree->fMuIsLoose) {
      pass = true;
    }    
    return pass;
}

Bool_t passTightID( MuonTree *muTree) {
    bool pass = false;
    if (muTree->fMuIsTight) {
      pass = true;
    }    
    return pass;
}


Bool_t passVetoIso( MuonTree *muTree) {

    bool pass = false;
    if ( muTree->fMuPFIso04 < 0.4
	 ) {
      pass = true;
    }
  
    return pass;
}



//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void MakeMuonIsoPerformancePlots(string InputFile, string Label, Int_t Option)
{  

  string label = "";
  if (Label != "") label = "_" + Label;


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *MuPFIso_Real = new TH1F(("MuPFIso_Real"+label).c_str(), "; PFIso ; Number of Events ",  10000, 0 , 200);
  TH1F *MuPFRelIso_Real = new TH1F(("MuPFRelIso_Real"+label).c_str(), "; PFRelIso ; Number of Events ",  10000, 0 , 10);
  TH1F *MuPtRel_Real = new TH1F(("MuPtRel_Real"+label).c_str(), "; PtRel ; Number of Events ",  10000, 0 , 10);
  TH1F *MuMiniIso_Real = new TH1F(("MuMiniIso_Real"+label).c_str(), "; MiniIso ; Number of Events ",  10000, 0 , 10);
  TH1F *MuPFIso_Fake = new TH1F(("MuPFIso_Fake"+label).c_str(), "; PFIso ; Number of Events ",  10000, 0 , 200);
  TH1F *MuPFRelIso_Fake = new TH1F(("MuPFRelIso_Fake"+label).c_str(), "; PFRelIso ; Number of Events ",  10000, 0 , 10);
  TH1F *MuPtRel_Fake = new TH1F(("MuPtRel_Fake"+label).c_str(), "; PtRel ; Number of Events ",  10000, 0 , 10);
  TH1F *MuMiniIso_Fake = new TH1F(("MuMiniIso_Fake"+label).c_str(), "; MiniIso ; Number of Events ",  10000, 0 , 10);

  Double_t RealMuons = 0;
  Double_t FakeMuons = 0;
  Double_t RealMuonPassVetoIso = 0;
  Double_t FakeMuonPassVetoIso = 0;
  Double_t RealMuonPassMiniIso = 0;
  Double_t FakeMuonPassMiniIso = 0;

  //*****************************************************************************************
  //MuTree
  //*****************************************************************************************
  MuonTree *MuTree = new MuonTree;
  MuTree->LoadTree(InputFile.c_str());
  MuTree->InitTree(MuonTree::kMuTreeLight);

  cout << "Total Entries: " << MuTree->tree_->GetEntries() << "\n";
  int nentries = MuTree->tree_->GetEntries();
  nentries = 5000000;
  for(UInt_t ientry=0; ientry < MuTree->tree_->GetEntries(); ientry++) {       	
    MuTree->tree_->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
            
    if (abs(MuTree->fPdgId) == 13 && RealMuons >= nentries) continue;
    if ( (abs(MuTree->fPdgId) > 50 || MuTree->fPdgId == 0) && FakeMuons >= nentries) continue;
    

    //Cuts
    if (MuTree->fMuPt < 5) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(MuTree->fMuEta) < 1.485) subdet = 0;
    else subdet = 1;

    Int_t ptBin = 0;
    if (MuTree->fMuPt > 10.0) ptBin = 1;
    if (MuTree->fMuPt > 20.0) ptBin = 2;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 3) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 0 && ptBin == 2);
    if (Option == 5) passCuts = (subdet == 1 && ptBin == 2);
    if (Option == 10) passCuts = (ptBin == 0 );
    if (Option == 11) passCuts = (ptBin == 1 );
    if (Option == 12) passCuts = (ptBin == 2 );

    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    


    //Some Preselection cuts
    if (!passLooseID(MuTree)) continue;


    //Real Muon
    if (abs(MuTree->fPdgId) == 13) {
      //don't consider muons close to partons
      if(!(MuTree->fDRToClosestParton > 1.0)) continue;
      RealMuons += MuTree->fWeight;

      if (passVetoIso(MuTree)) RealMuonPassVetoIso += MuTree->fWeight;
      if (MuTree->fMiniIso < 0.2) RealMuonPassMiniIso += MuTree->fWeight;
      MuPFRelIso_Real->Fill(MuTree->fMuPFIso04, MuTree->fWeight);
      MuPFIso_Real->Fill(MuTree->fMuPFIso04*MuTree->fMuPt, MuTree->fWeight);
      MuPtRel_Real->Fill(MuTree->fPtRel, MuTree->fWeight);
      MuMiniIso_Real->Fill(MuTree->fMiniIso, MuTree->fWeight);

    } 

    //FakeMuon
    else if ( abs(MuTree->fPdgId) > 50 || MuTree->fPdgId == 0) {
      FakeMuons += MuTree->fWeight;

      if (passVetoIso(MuTree)) FakeMuonPassVetoIso += MuTree->fWeight;      
      if (MuTree->fMiniIso < 0.2) FakeMuonPassMiniIso += MuTree->fWeight;
      MuPFRelIso_Fake->Fill(MuTree->fMuPFIso04, MuTree->fWeight);
      MuPFIso_Fake->Fill(MuTree->fMuPFIso04*MuTree->fMuPt, MuTree->fWeight);
      MuPtRel_Fake->Fill(MuTree->fPtRel, MuTree->fWeight);
      MuMiniIso_Fake->Fill(MuTree->fMiniIso, MuTree->fWeight);
 
    }

  } 
  
  
  //*****************************************************************************************
  //Current Working Points
  //*****************************************************************************************
  cout << " VetoIso Real Muon Efficiency : " << RealMuonPassVetoIso << " / " << RealMuons << " = " << RealMuonPassVetoIso/RealMuons << endl;
  cout << " VetoIso Fake Muon Efficiency : " << FakeMuonPassVetoIso << " / " << FakeMuons << " = " << FakeMuonPassVetoIso/FakeMuons << endl;
  TGraphAsymmErrors* ROC_VetoIsoWP = MakeCurrentWPSigEffVsBkgEffGraph(RealMuonPassVetoIso/RealMuons , FakeMuonPassVetoIso/FakeMuons, "ROC_VetoIsoWP"+label);

  cout << "CSA14 VetoIso Real Muon Efficiency : " << RealMuonPassMiniIso << " / " << RealMuons << " = " << RealMuonPassMiniIso/RealMuons << endl;
  cout << "CSA14 VetoIso Fake Muon Efficiency : " << FakeMuonPassMiniIso << " / " << FakeMuons << " = " << FakeMuonPassMiniIso/FakeMuons << endl;
  TGraphAsymmErrors* ROC_MiniIsoWP = MakeCurrentWPSigEffVsBkgEffGraph(RealMuonPassMiniIso/RealMuons , FakeMuonPassMiniIso/FakeMuons, "ROC_MiniIsoWP"+label);

  Double_t BkgEffVetoIso = FakeMuonPassVetoIso/FakeMuons;
  Double_t SigEffVetoIso = RealMuonPassVetoIso/RealMuons;


  //*****************************************************************************************
  //Make ROC curves
  //*****************************************************************************************
  TGraphAsymmErrors* ROC_PFIso = MakeSigEffVsBkgEffGraph(MuPFIso_Real, MuPFIso_Fake, "ROC_PFIso"+label );
  TGraphAsymmErrors* ROC_PFRelIso = MakeSigEffVsBkgEffGraph(MuPFRelIso_Real, MuPFRelIso_Fake, "ROC_PFRelIso"+label );
  TGraphAsymmErrors* ROC_PtRel = MakeSigEffVsBkgEffGraph(MuPtRel_Real, MuPtRel_Fake, "ROC_PtRel"+label,false );
  TGraphAsymmErrors* ROC_MiniIso = MakeSigEffVsBkgEffGraph(MuMiniIso_Real, MuMiniIso_Fake, "ROC_MiniIso"+label );


  //*****************************************************************************************
  //Find Cut with same signal efficiency Make ROC curves
  //*****************************************************************************************
  Double_t CutValue_PFIso_98 = FindCutValueAtFixedEfficiency(MuPFIso_Real, 0.98 );
  Double_t CutValue_PFIso_95 = FindCutValueAtFixedEfficiency(MuPFIso_Real, 0.95 );
  Double_t CutValue_PFIso_90 = FindCutValueAtFixedEfficiency(MuPFIso_Real, 0.90 );
  Double_t CutValue_PFIso_85 = FindCutValueAtFixedEfficiency(MuPFIso_Real, 0.85 );
  Double_t CutValue_PFIso_80 = FindCutValueAtFixedEfficiency(MuPFIso_Real, 0.80 );
  cout << "PFIso Cut Value @ 98% Sig Eff: " << CutValue_PFIso_98 << endl;
  cout << "PFIso Cut Value @ 95% Sig Eff: " << CutValue_PFIso_95 << endl;
  cout << "PFIso Cut Value @ 90% Sig Eff: " << CutValue_PFIso_90 << endl;
  cout << "PFIso Cut Value @ 85% Sig Eff: " << CutValue_PFIso_85 << endl;
  cout << "PFIso Cut Value @ 80% Sig Eff: " << CutValue_PFIso_80 << endl;
  Double_t CutValue_PFRelIso_98 = FindCutValueAtFixedEfficiency(MuPFRelIso_Real, 0.98 );
  Double_t CutValue_PFRelIso_95 = FindCutValueAtFixedEfficiency(MuPFRelIso_Real, 0.95 );
  Double_t CutValue_PFRelIso_90 = FindCutValueAtFixedEfficiency(MuPFRelIso_Real, 0.90 );
  Double_t CutValue_PFRelIso_85 = FindCutValueAtFixedEfficiency(MuPFRelIso_Real, 0.85 );
  Double_t CutValue_PFRelIso_80 = FindCutValueAtFixedEfficiency(MuPFRelIso_Real, 0.80 );
  cout << "PFRelIso Cut Value @ 98% Sig Eff: " << CutValue_PFRelIso_98 << endl;
  cout << "PFRelIso Cut Value @ 95% Sig Eff: " << CutValue_PFRelIso_95 << endl;
  cout << "PFRelIso Cut Value @ 90% Sig Eff: " << CutValue_PFRelIso_90 << endl;
  cout << "PFRelIso Cut Value @ 85% Sig Eff: " << CutValue_PFRelIso_85 << endl;
  cout << "PFRelIso Cut Value @ 80% Sig Eff: " << CutValue_PFRelIso_80 << endl;


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
  plotname = "MuonIso"+label;

  // ROCGraphs.push_back(ROC_PFIso);
  // GraphLabels.push_back("PFIso");
  // colors.push_back(kGreen+2);
  
  ROCGraphs.push_back(ROC_PFRelIso);
  GraphLabels.push_back("PFRelIso");
  colors.push_back(kRed);
  
  ROCGraphs.push_back(ROC_MiniIso);
  GraphLabels.push_back("MiniIso");
  colors.push_back(kBlue);
  


  //*****************************************************************************************
  Double_t xmin = 0.0;
  Double_t xmax = 1.0;
  Double_t ymin = 0.0;
  Double_t ymax = 1.0;



  cv = new TCanvas("cv", "cv", 800, 600);

//    legend = new TLegend(0.45,0.20,0.75,0.50);
  legend = new TLegend(0.50,0.14,0.85,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
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

  legend->AddEntry(ROC_VetoIsoWP, "Razor PHYS14 Veto Iso WP", "P");
  ROC_VetoIsoWP->SetFillColor(kGreen+3);
  ROC_VetoIsoWP->SetMarkerColor(kGreen+3);
  ROC_VetoIsoWP->SetMarkerStyle(34);
  ROC_VetoIsoWP->SetMarkerSize(2.5);
  ROC_VetoIsoWP->Draw("Psame");
 
  legend->AddEntry(ROC_MiniIsoWP, "MiniIso < 0.2", "P");
  ROC_MiniIsoWP->SetFillColor(kBlack);
  ROC_MiniIsoWP->SetMarkerColor(kBlack);
  ROC_MiniIsoWP->SetMarkerStyle(34);
  ROC_MiniIsoWP->SetMarkerSize(2.5);
  ROC_MiniIsoWP->Draw("Psame");


  legend->Draw();
  
  cv->SaveAs(("ROCGraphs_" + plotname + ".gif").c_str());




  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("MuonIsoROCCurves"+Label+".root").c_str(), "RECREATE");
  file->cd();
  file->WriteTObject(ROC_PFRelIso, ("ROC_PFRelIso" + Label).c_str(), "WriteDelete");
  file->WriteTObject(ROC_PFIso, ("ROC_PFIso" + Label).c_str(), "WriteDelete");
  
  file->Close();
  delete file;       

  gBenchmark->Show("WWTemplate");       
} 

