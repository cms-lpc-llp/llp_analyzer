//root -l /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/macros/ObjectStudies/MakeElectronIDMVAPerformancePlots.C+'("/afs/cern.ch/user/s/sixie/work/public/Run2SUSY/ElectronNtuple/ElectronNtuple.TTJets.25ns.root","All",-1)'
//root -l /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/macros/ObjectStudies/MakeElectronIDMVAPerformancePlots.C+'("/afs/cern.ch/user/s/sixie/work/public/Run2SUSY/ElectronNtuple/ElectronNtuple_TTJetsPrompt_T1bbbbFake.root","All",-1)'
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

#include "RazorAnalyzer/include/ElectronTree.h"

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


Bool_t passCSA14Preselection( ElectronTree *eleTree) {

    bool pass = false;
    if(fabs(eleTree->fEleSCEta) < 1.479) {
      if ( (0==0)  
	   && fabs(eleTree->fEleDEtaIn) < 0.0218
	   && fabs(eleTree->fEleDPhiIn) < 0.2781
	   && eleTree->fEleSigmaIEtaIEta < 0.017
	   && eleTree->fEleHoverE < 0.3514
	   && fabs(eleTree->fEleD0) < 0.1069
	   && fabs(eleTree->fEleDZ) < 0.8137
	   && eleTree->fEleOneOverEMinusOneOverP < 0.3309
	   && eleTree->fElePFIso04 < 1.8402
	   && eleTree->fElePassConversion
	   && eleTree->fEleNMissHits < 2
	   ) {
	pass = true;
      }
    } else {
      if ( (0==0)
	   && fabs(eleTree->fEleDEtaIn) < 0.028
	   && fabs(eleTree->fEleDPhiIn) < 0.2592
	   && eleTree->fEleSigmaIEtaIEta < 0.0493
	   && eleTree->fEleHoverE < 0.3204
	   && fabs(eleTree->fEleD0) < 0.3377
	   && fabs(eleTree->fEleDZ) <  0.9532
	   && eleTree->fEleOneOverEMinusOneOverP <  0.1666
	   && eleTree->fElePFIso04 < 1.3637
	   && eleTree->fElePassConversion
	   && eleTree->fEleNMissHits < 3
	   ) {
	pass = true;
      }
    } 
    return pass;
}

Bool_t passCSA14Tight( ElectronTree *eleTree) {

    bool pass = false;
    if(fabs(eleTree->fEleSCEta) < 1.479) {
      if ( fabs(eleTree->fEleDEtaIn) < 0.0091
	   && fabs(eleTree->fEleDPhiIn) < 0.031
	   && eleTree->fEleSigmaIEtaIEta < 0.0106
	   && eleTree->fEleHoverE < 0.0532
	   && fabs(eleTree->fEleD0) < 0.0126
	   && fabs(eleTree->fEleDZ) < 0.0116
	   && eleTree->fEleOneOverEMinusOneOverP < 0.0609
	   //&& eleTree->fElePFIso04 < 0.1649
	   && eleTree->fElePassConversion
	   && eleTree->fEleNMissHits < 1
	   ) {
	pass = true;
      }
    } else {
      if (fabs(eleTree->fEleDEtaIn) < 0.0106
	   && fabs(eleTree->fEleDPhiIn) < 0.0359
	   && eleTree->fEleSigmaIEtaIEta < 0.0305 
	   && eleTree->fEleHoverE < 0.0835 
	  && fabs(eleTree->fEleD0) < 0.0163 
	   && fabs(eleTree->fEleDZ) <  0.5999
	   && eleTree->fEleOneOverEMinusOneOverP <  0.1126
	  // && eleTree->fElePFIso04 < 0.2075
	   && eleTree->fElePassConversion
	   && eleTree->fEleNMissHits < 1
	  ) {
	pass = true;
      }
    } 
    return pass;
}

Bool_t passCSA14Loose( ElectronTree *eleTree) {

    bool pass = false;
    if(fabs(eleTree->fEleSCEta) < 1.479) {
      if ( fabs(eleTree->fEleDEtaIn) < 0.0181
	   && fabs(eleTree->fEleDPhiIn) < 0.0936
	   && eleTree->fEleSigmaIEtaIEta < 0.0123
	   && eleTree->fEleHoverE < 0.141
	   && fabs(eleTree->fEleD0) < 0.0166
	   && fabs(eleTree->fEleDZ) < 0.54342
	   && eleTree->fEleOneOverEMinusOneOverP < 0.1353
	   // && eleTree->fElePFIso04 < 0.24
	   && eleTree->fElePassConversion
	   && eleTree->fEleNMissHits < 1
	   ) {
	pass = true;
      }
    } else {
      if (fabs(eleTree->fEleDEtaIn) < 0.0124
	   && fabs(eleTree->fEleDPhiIn) < 0.0642
	   && eleTree->fEleSigmaIEtaIEta < 0.035 
	   && eleTree->fEleHoverE < 0.1115 
	  && fabs(eleTree->fEleD0) < 0.098 
	   && fabs(eleTree->fEleDZ) <  0.9187
	   && eleTree->fEleOneOverEMinusOneOverP <  0.1443
	  //&& eleTree->fElePFIso04 < 0.3529
	   && eleTree->fElePassConversion
	   && eleTree->fEleNMissHits < 1
	  ) {
	pass = true;
      }
    } 
    return pass;
}

Bool_t passCSA14Veto( ElectronTree *eleTree) {

    bool pass = false;
    if(fabs(eleTree->fEleSCEta) < 1.479) {
      if ( fabs(eleTree->fEleDEtaIn) < 0.02
	   && fabs(eleTree->fEleDPhiIn) < 0.2579
	   && eleTree->fEleSigmaIEtaIEta < 0.0125
	   && eleTree->fEleHoverE < 0.2564
	   && fabs(eleTree->fEleD0) < 0.025
	   && fabs(eleTree->fEleDZ) < 0.5863
	   && eleTree->fEleOneOverEMinusOneOverP < 0.1508
	   // && eleTree->fElePFIso04 < 0.3313 
	   && eleTree->fElePassConversion
	   && eleTree->fEleNMissHits < 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(eleTree->fEleDEtaIn) < 0.0141
	   && fabs(eleTree->fEleDPhiIn) < 0.2591
	   && eleTree->fEleSigmaIEtaIEta < 0.0371 
	   && eleTree->fEleHoverE < 0.1335 
	  && fabs(eleTree->fEleD0) < 0.2232
	   && fabs(eleTree->fEleDZ) <  0.9513
	   && eleTree->fEleOneOverEMinusOneOverP <  0.1542
	  //&& eleTree->fElePFIso04 < 0.3816
	  && eleTree->fElePassConversion
	  && eleTree->fEleNMissHits < 3
	  ) {
	pass = true;
      }
    } 
    return pass;
}



Bool_t passIDMVANonTrigVeto( ElectronTree *eleTree) {

  Int_t subdet = 0;  
  if (fabs(eleTree->fEleSCEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (eleTree->fElePt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = 0.0;
  if (MVABin == 1) MVACut = 0.6; 
  if (MVABin == 2) MVACut = -0.3;
  if (MVABin == 3) MVACut = 0.5;  

    bool pass = false;
    if (eleTree->fIDMVANonTrig > MVACut 

	&& ( (fabs(eleTree->fEleSCEta) < 1.479 && fabs(eleTree->fEleD0) < 0.0166)
	     ||
	     (fabs(eleTree->fEleSCEta) >= 1.479 && fabs(eleTree->fEleD0) < 0.098)
	     )
	//&& eleTree->fElePFIso04 < 0.4
	) {
      pass = true;
    }   
    return pass;
}

Bool_t passIDMVATrigTight( ElectronTree *eleTree) {

    bool pass = false;
    if(fabs(eleTree->fEleSCEta) < 1.479) {
      if ( eleTree->fIDMVATrig > 0.99
	   ) {
	pass = true;
      }
    } else {
      if ( eleTree->fIDMVATrig > 0.99	
	   ) {
	pass = true;
      }
    } 
    return pass;
}


Bool_t passCutBasedIsoOnly(Double_t fElePt, 
                    Double_t fEleEta, 
                    Double_t fElePFIso ) {

  //apply full isolation cut
  Bool_t passIsoCuts = kFALSE;
  if (fElePt >= 10 && fElePt < 20) passIsoCuts = ( fElePFIso  < 0.09 ); 
  if (fElePt >= 20) passIsoCuts = ( fElePFIso  < 0.13 ); 

  return passIsoCuts;
}


Bool_t passImprovedIso( ElectronTree *eleTree) {

    bool pass = false;
    if(fabs(eleTree->fEleSCEta) < 1.479) {
      if ( eleTree->fElePt > 20) {
	if ( 
	    eleTree->fElePFIso04 < 0.3 
	     ) {
	pass = true;
	}
      } else {
	if ( eleTree->fElePFIso04*eleTree->fElePt < 5 ) {
	  pass = true;
	}
      }
    } else {
      if ( eleTree->fElePt > 20) {
        if (
            eleTree->fElePFIso04 < 0.3
	    ) {
	  pass = true;
        }
      } else {
	if ( eleTree->fElePFIso04*eleTree->fElePt < 5 ) {
          pass = true;
	}
      }
    } 
    return pass;
}



//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void MakeElectronIDMVAPerformancePlots(string InputFile, string Label, Int_t Option)
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
  TH1F *EleIDMVATrig_Real = new TH1F(("EleIDMVATrig_Real"+label).c_str(), "; IDMVATrig ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDMVANonTrig_Real = new TH1F(("EleIDMVANonTrig_Real"+label).c_str(), "; IDMVANonTrig ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDMVATrig_Fake = new TH1F(("EleIDMVATrig_Fake"+label).c_str(), "; IDMVATrig ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDMVANonTrig_Fake = new TH1F(("EleIDMVANonTrig_Fake"+label).c_str(), "; IDMVANonTrig ; Number of Events ",  10000, -2 , 2);

  Double_t RealElectrons = 0;
  Double_t FakeElectrons = 0;
  Double_t RealElectronPassCSA14Tight = 0;
  Double_t FakeElectronPassCSA14Tight = 0;
  Double_t RealElectronPassCSA14Loose = 0;
  Double_t FakeElectronPassCSA14Loose = 0;
  Double_t RealElectronPassCSA14Veto = 0;
  Double_t FakeElectronPassCSA14Veto = 0;

  Double_t RealElectronPassIDMVANonTrigVeto = 0;
  Double_t FakeElectronPassIDMVANonTrigVeto = 0;
  Double_t RealElectronPassIDMVATrigTight = 0;
  Double_t FakeElectronPassIDMVATrigTight = 0;

  //*****************************************************************************************
  //EleTree
  //*****************************************************************************************
  ElectronTree *EleTree = new ElectronTree;
  EleTree->LoadTree(InputFile.c_str());
  EleTree->InitTree(ElectronTree::kEleTreeLight);

  cout << "Total Entries: " << EleTree->tree_->GetEntries() << "\n";
  int nentries = EleTree->tree_->GetEntries();
  nentries = 5000000;
  for(UInt_t ientry=0; ientry < EleTree->tree_->GetEntries(); ientry++) {       	
    EleTree->tree_->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
        
    
    if (abs(EleTree->fPdgId) == 11 && RealElectrons >= nentries) continue;
    if ( (abs(EleTree->fPdgId) > 50 || EleTree->fPdgId == 0) && FakeElectrons >= nentries) continue;

    //don't evaluate performance using training events
    if (EleTree->fElePt < 5) continue;


    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(EleTree->fEleEta) < 1.485) subdet = 0;
    else subdet = 1;

    Int_t ptBin = 0;
    if (EleTree->fElePt > 10.0) ptBin = 1;
    if (EleTree->fElePt > 20.0) ptBin = 2;

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
    //if (!passCSA14Preselection(EleTree)) continue;
    //if (!EleTree->fElePassConversion) continue;
    if (!(fabs(EleTree->fEleDZ) < 1)) continue;
    //if (!(EleTree->fElePFIso04 < 0.2)) continue;
    if (!passImprovedIso(EleTree)) continue;
    //if (!(fabs(EleTree->fEleD0) < 0.0166)) continue;

    //Real Electron
    if (abs(EleTree->fPdgId) == 11) {
      RealElectrons += EleTree->fWeight;
      if (passCSA14Tight(EleTree)) RealElectronPassCSA14Tight += EleTree->fWeight;
      if (passCSA14Loose(EleTree)) RealElectronPassCSA14Loose += EleTree->fWeight;
      if (passCSA14Veto(EleTree)) RealElectronPassCSA14Veto += EleTree->fWeight;
      if (passIDMVANonTrigVeto(EleTree)) RealElectronPassIDMVANonTrigVeto += EleTree->fWeight;
      if (passIDMVATrigTight(EleTree)) RealElectronPassIDMVATrigTight += EleTree->fWeight;      

      if ( EleTree->fElePFIso04 < 0.2 ) {
	EleIDMVATrig_Real->Fill(EleTree->fIDMVATrig, EleTree->fWeight);
      } else {
	EleIDMVATrig_Real->Fill(-1, EleTree->fWeight);
      }

      if ( EleTree->fElePFIso04 < 99 &&  
	   ( (fabs(EleTree->fEleSCEta) < 1.479 && fabs(EleTree->fEleD0) < 0.0166)
	     ||
	     (fabs(EleTree->fEleSCEta) >= 1.479 && fabs(EleTree->fEleD0) < 0.098)
	     )
	   ) {
	EleIDMVANonTrig_Real->Fill(EleTree->fIDMVANonTrig, EleTree->fWeight);
      } else {
	EleIDMVANonTrig_Real->Fill(-1, EleTree->fWeight);
      }

    } 
    //FakeElectron
    else if ( abs(EleTree->fPdgId) > 50 || EleTree->fPdgId == 0) {
      FakeElectrons += EleTree->fWeight;
      if (passCSA14Tight(EleTree)) FakeElectronPassCSA14Tight += EleTree->fWeight;
      if (passCSA14Loose(EleTree)) FakeElectronPassCSA14Loose += EleTree->fWeight;
      if (passCSA14Veto(EleTree)) FakeElectronPassCSA14Veto += EleTree->fWeight;
      if (passIDMVANonTrigVeto(EleTree)) FakeElectronPassIDMVANonTrigVeto += EleTree->fWeight;
      if (passIDMVATrigTight(EleTree)) FakeElectronPassIDMVATrigTight += EleTree->fWeight;
      
      if ( EleTree->fElePFIso04 < 0.2 ) {
	EleIDMVATrig_Fake->Fill(EleTree->fIDMVATrig, EleTree->fWeight);
      } else {
	EleIDMVATrig_Real->Fill(-1, EleTree->fWeight);
      }

      if ( EleTree->fElePFIso04 < 99 &&  
	   ( (fabs(EleTree->fEleSCEta) < 1.479 && fabs(EleTree->fEleD0) < 0.0166)
	     ||
	     (fabs(EleTree->fEleSCEta) >= 1.479 && fabs(EleTree->fEleD0) < 0.098)
	     )
	   ) {
	EleIDMVANonTrig_Fake->Fill(EleTree->fIDMVANonTrig, EleTree->fWeight);
      } else {
	EleIDMVANonTrig_Fake->Fill(-1, EleTree->fWeight);
      }
    }

  } 
  
  
  //*****************************************************************************************
  //Current Working Points
  //*****************************************************************************************
  cout << "CSA14 Tight Real Electron Efficiency : " << RealElectronPassCSA14Tight << " / " << RealElectrons << " = " << RealElectronPassCSA14Tight/RealElectrons << endl;
  cout << "CSA14 Tight Fake Electron Efficiency : " << FakeElectronPassCSA14Tight << " / " << FakeElectrons << " = " << FakeElectronPassCSA14Tight/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_CSA14TightWP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassCSA14Tight/RealElectrons , FakeElectronPassCSA14Tight/FakeElectrons, "ROC_CSA14TightWP"+label);

  cout << "CSA14 Loose Real Electron Efficiency : " << RealElectronPassCSA14Loose << " / " << RealElectrons << " = " << RealElectronPassCSA14Loose/RealElectrons << endl;
  cout << "CSA14 Loose Fake Electron Efficiency : " << FakeElectronPassCSA14Loose << " / " << FakeElectrons << " = " << FakeElectronPassCSA14Loose/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_CSA14LooseWP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassCSA14Loose/RealElectrons , FakeElectronPassCSA14Loose/FakeElectrons, "ROC_CSA14LooseWP"+label);

  cout << "CSA14 Veto Real Electron Efficiency : " << RealElectronPassCSA14Veto << " / " << RealElectrons << " = " << RealElectronPassCSA14Veto/RealElectrons << endl;
  cout << "CSA14 Veto Fake Electron Efficiency : " << FakeElectronPassCSA14Veto << " / " << FakeElectrons << " = " << FakeElectronPassCSA14Veto/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_CSA14VetoWP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassCSA14Veto/RealElectrons , FakeElectronPassCSA14Veto/FakeElectrons, "ROC_CSA14VetoWP"+label);



  Double_t BkgEffCSA14Tight = FakeElectronPassCSA14Tight/FakeElectrons;
  Double_t SigEffCSA14Tight = RealElectronPassCSA14Tight/RealElectrons;
  Double_t BkgEffCSA14Veto = FakeElectronPassCSA14Veto/FakeElectrons;
  Double_t SigEffCSA14Veto = RealElectronPassCSA14Veto/RealElectrons;

  cout << "**********************\n";
  Double_t SigEffIDMVATrig_AtTightBkgEff = FindSigEffAtFixedBkgEfficiency(EleIDMVATrig_Real, EleIDMVATrig_Fake, BkgEffCSA14Tight);
  Double_t SigEffIDMVANonTrig_AtTightBkgEff = FindSigEffAtFixedBkgEfficiency(EleIDMVANonTrig_Real, EleIDMVANonTrig_Fake, BkgEffCSA14Tight);
  Double_t BkgEffIDMVATrig_AtTightSigEff = FindBkgEffAtFixedSignalEfficiency(EleIDMVATrig_Real, EleIDMVATrig_Fake, SigEffCSA14Tight);
  Double_t BkgEffIDMVANonTrig_AtTightSigEff = FindBkgEffAtFixedSignalEfficiency(EleIDMVANonTrig_Real, EleIDMVANonTrig_Fake, SigEffCSA14Tight);
  cout << "Signal Efficiency (wrt CSA14Tight Cut-based) for : same bkg \n";
  cout << "IDMVATrig : " << SigEffIDMVATrig_AtTightBkgEff/SigEffCSA14Tight <<  endl;
  cout << "IDMVANonTrig : " << SigEffIDMVANonTrig_AtTightBkgEff/SigEffCSA14Tight <<  endl;
  cout << "Bkg Efficiency (wrt CSA14Veto Cut-based) for same sig eff \n";
  cout << "IDMVATrig : " << BkgEffIDMVATrig_AtTightSigEff/BkgEffCSA14Tight << endl;
  cout << "IDMVANonTrig : " << BkgEffIDMVANonTrig_AtTightSigEff/BkgEffCSA14Tight << endl;
  cout << "**********************\n";

  cout << "**********************\n";
  Double_t BkgEffIDMVATrig_AtVetoSigEff = FindBkgEffAtFixedSignalEfficiency(EleIDMVATrig_Real, EleIDMVATrig_Fake, SigEffCSA14Veto);
  Double_t BkgEffIDMVANonTrig_AtVetoSigEff = FindBkgEffAtFixedSignalEfficiency(EleIDMVANonTrig_Real, EleIDMVANonTrig_Fake, SigEffCSA14Veto);
  Double_t SigEffIDMVATrig_AtVetoBkgEff = FindSigEffAtFixedBkgEfficiency(EleIDMVATrig_Real, EleIDMVATrig_Fake, BkgEffCSA14Veto);
  Double_t SigEffIDMVANonTrig_AtVetoBkgEff = FindSigEffAtFixedBkgEfficiency(EleIDMVANonTrig_Real, EleIDMVANonTrig_Fake, BkgEffCSA14Veto);
  cout << "Sig Efficiency (wrt CSA14Veto Cut-based) for same bkg eff \n";
  cout << "IDMVATrig : " << SigEffIDMVATrig_AtVetoBkgEff/SigEffCSA14Veto << endl;
  cout << "IDMVANonTrig : " << SigEffIDMVANonTrig_AtVetoBkgEff/SigEffCSA14Veto << endl;
  cout << "Bkg Efficiency (wrt CSA14Veto Cut-based) for same sig eff \n";
  cout << "IDMVATrig : " << BkgEffIDMVATrig_AtVetoSigEff/BkgEffCSA14Veto << endl;
  cout << "IDMVANonTrig : " << BkgEffIDMVANonTrig_AtVetoSigEff/BkgEffCSA14Veto << endl;
  cout << "**********************\n";


  cout << "IDMVANonTrig Veto Real Electron Efficiency : " << RealElectronPassIDMVANonTrigVeto << " / " << RealElectrons << " = " << RealElectronPassIDMVANonTrigVeto/RealElectrons << endl;
  cout << "IDMVANonTrig Veto Fake Electron Efficiency : " << FakeElectronPassIDMVANonTrigVeto << " / " << FakeElectrons << " = " << FakeElectronPassIDMVANonTrigVeto/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_IDMVANonTrigVetoWP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassIDMVANonTrigVeto/RealElectrons , FakeElectronPassIDMVANonTrigVeto/FakeElectrons, "ROC_IDMVANonTrigVetoWP"+label);


  cout << "IDMVATrig Tight Real Electron Efficiency : " << RealElectronPassIDMVATrigTight << " / " << RealElectrons << " = " << RealElectronPassIDMVATrigTight/RealElectrons << endl;
  cout << "IDMVATrig Tight Veto Fake Electron Efficiency : " << FakeElectronPassIDMVATrigTight << " / " << FakeElectrons << " = " << FakeElectronPassIDMVATrigTight/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_IDMVATrigTightWP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassIDMVATrigTight/RealElectrons , FakeElectronPassIDMVATrigTight/FakeElectrons, "ROC_IDMVATrigTightWP"+label);



  //*****************************************************************************************
  //Make ROC curves
  //*****************************************************************************************
  TGraphAsymmErrors* ROC_IDMVATrig = MakeSigEffVsBkgEffGraph(EleIDMVATrig_Real, EleIDMVATrig_Fake, "ROC_IDMVATrig"+label );
  TGraphAsymmErrors* ROC_IDMVANonTrig = MakeSigEffVsBkgEffGraph(EleIDMVANonTrig_Real, EleIDMVANonTrig_Fake, "ROC_IDMVANonTrig"+label );


  //*****************************************************************************************
  //Find Cut with same signal efficiency Make ROC curves
  //*****************************************************************************************
  Double_t CutValue_IDMVATrig_SameSig = FindCutValueAtFixedEfficiency(EleIDMVATrig_Real, SigEffCSA14Tight );
  Double_t CutValue_IDMVATrig_SameBkg = FindCutValueAtFixedEfficiency(EleIDMVATrig_Fake, BkgEffCSA14Tight );
  cout << "IDMVATrig Cut Value @ Same Cut-Based Tight Sig Eff: " << CutValue_IDMVATrig_SameSig << endl;
  cout << "IDMVATrig Cut Value @ Same Cut-Based Tight Bkg Eff: " << CutValue_IDMVATrig_SameBkg << endl;

  Double_t CutValue_IDMVANonTrig_SameSig = FindCutValueAtFixedEfficiency(EleIDMVANonTrig_Real, SigEffCSA14Veto );
  Double_t CutValue_IDMVANonTrig_SameBkg = FindCutValueAtFixedEfficiency(EleIDMVANonTrig_Fake, BkgEffCSA14Veto );
  cout << "IDMVANonTrig Cut Value @ Same Cut-Based Veto Sig Eff: " << CutValue_IDMVANonTrig_SameSig << endl;
  cout << "IDMVANonTrig Cut Value @ Same Cut-Based Veto Bkg Eff: " << CutValue_IDMVANonTrig_SameBkg << endl;

  Double_t CutValue_IDMVANonTrig_HalfBkg = FindCutValueAtFixedEfficiency(EleIDMVANonTrig_Fake, 0.5*FakeElectronPassIDMVANonTrigVeto/FakeElectrons );
  cout << "IDMVANonTrig Cut Value @ 50% IDMVA NonTrig Veto Bkg Eff: " << CutValue_IDMVANonTrig_HalfBkg << endl;

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
  plotname = "ElectronIDMVA"+label;

  // ROCGraphs.push_back(ROC_IDMVATrig);
  // GraphLabels.push_back("IDMVATrig");
  // colors.push_back(kBlue);
  
  ROCGraphs.push_back(ROC_IDMVANonTrig);
  GraphLabels.push_back("IDMVANonTrig");
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

  // legend->AddEntry(ROC_CSA14TightWP, "CSA14Tight WP", "P");
  // ROC_CSA14TightWP->SetFillColor(kRed);
  // ROC_CSA14TightWP->SetMarkerColor(kRed);
  // ROC_CSA14TightWP->SetMarkerStyle(34);
  // ROC_CSA14TightWP->SetMarkerSize(2.5);
  // ROC_CSA14TightWP->Draw("Psame");

  legend->AddEntry(ROC_CSA14LooseWP, "CSA14Loose WP", "P");
  ROC_CSA14LooseWP->SetFillColor(kBlue);
  ROC_CSA14LooseWP->SetMarkerColor(kBlue);
  ROC_CSA14LooseWP->SetMarkerStyle(34);
  ROC_CSA14LooseWP->SetMarkerSize(2.5);
  ROC_CSA14LooseWP->Draw("Psame");

  legend->AddEntry(ROC_CSA14VetoWP, "CSA14Veto WP", "P");
  ROC_CSA14VetoWP->SetFillColor(kGreen+3);
  ROC_CSA14VetoWP->SetMarkerColor(kGreen+3);
  ROC_CSA14VetoWP->SetMarkerStyle(34);
  ROC_CSA14VetoWP->SetMarkerSize(2.5);
  ROC_CSA14VetoWP->Draw("Psame");

   legend->AddEntry(ROC_IDMVANonTrigVetoWP, "IDMVANonTrigVeto WP", "P");
   ROC_IDMVANonTrigVetoWP->SetFillColor(kBlack);
   ROC_IDMVANonTrigVetoWP->SetMarkerColor(kBlack);
   ROC_IDMVANonTrigVetoWP->SetMarkerStyle(34);
   ROC_IDMVANonTrigVetoWP->SetMarkerSize(2.5);
   ROC_IDMVANonTrigVetoWP->Draw("Psame");

  // legend->AddEntry(ROC_IDMVATrigTightWP, "IDMVATrigTight WP", "P");
  // ROC_IDMVATrigTightWP->SetFillColor(kRed);
  // ROC_IDMVATrigTightWP->SetMarkerColor(kRed);
  // ROC_IDMVATrigTightWP->SetMarkerStyle(34);
  // ROC_IDMVATrigTightWP->SetMarkerSize(2.5);
  // ROC_IDMVATrigTightWP->Draw("Psame");


  legend->Draw();
  
  cv->SaveAs(("ROCGraphs_" + plotname + ".gif").c_str());

  gBenchmark->Show("WWTemplate");       
} 

