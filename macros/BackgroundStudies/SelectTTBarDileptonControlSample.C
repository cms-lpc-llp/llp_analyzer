//================================================================================================
//
// Simple Example
//
//root -l RazorAnalyzer/macros/ObjectStudies/MakeElectronEfficiencyPlots.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/ElectronNtuple/ElectronNtuple_PromptGenLevel_TTJets_25ns.root",-1,"Electron")'
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TH1D.h>                
#include <TH1F.h>                
#include <TH2F.h>                
#include <TCanvas.h>                
#include <TLegend.h> 
#include <THStack.h> 

#include "RazorAnalyzer/include/ControlSampleEvents.h"
#include "RazorAnalyzer/macros/tdrstyle.C"
#include "RazorAnalyzer/macros/CMS_lumi.C"

#endif

void PlotDataAndStackedBkg( vector<TH1D*> hist , vector<string> processLabels, vector<int> color,  bool hasData, string varName, string label ) {

  TCanvas *cv =0;
  TLegend *legend = 0;

  cv = new TCanvas("cv","cv", 800,700);
  cv->SetHighLightColor(2);
  cv->SetFillColor(0);
  cv->SetBorderMode(0);
  cv->SetBorderSize(2);
  cv->SetLeftMargin(0.16);
  cv->SetRightMargin(0.3);
  cv->SetTopMargin(0.07);
  cv->SetBottomMargin(0.12);
  cv->SetFrameBorderMode(0);  

  TPad *pad1 = new TPad("pad1","pad1", 0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(0.04);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.60,0.50,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stack = new THStack();
  TH1D *histDataOverMC = (TH1D*)hist[0]->Clone("histDataOverMC");

  if (hasData) {
    for (int i = hist.size()-1; i >= 1; --i) {
      hist[i]->SetFillColor(color[i]);
      hist[i]->SetFillStyle(1001);
      
      if ( hist[i]->Integral() > 0) {
  	stack->Add(hist[i]);
      }
    }
  } else {
    for (int i = hist.size()-1; i >= 0; --i) {
      hist[i]->SetFillColor(color[i]);
      hist[i]->SetFillStyle(1001);
      
      if ( hist[i]->Integral() > 0) {
  	stack->Add(hist[i]);
      }
    }
  }

  for (uint i = 0 ; i < hist.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(hist[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(hist[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stack->GetHists()->GetEntries() > 0) {
    stack->Draw("hist");
    stack->GetHistogram()->GetXaxis()->SetTitle(((TH1D*)(stack->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stack->GetHistogram()->GetYaxis()->SetTitle(((TH1D*)(stack->GetHists()->At(0)))->GetYaxis()->GetTitle());    
    stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.0);
    stack->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    stack->GetHistogram()->GetXaxis()->SetTitleSize(0.15);
    stack->SetMaximum( 1.2* fmax( stack->GetMaximum(), hist[0]->GetMaximum()) );
    stack->SetMinimum( 0.1 );

    if (hasData) {
      hist[0]->SetMarkerStyle(20);      
      hist[0]->SetMarkerSize(1);
      hist[0]->SetLineWidth(1);
      hist[0]->SetLineColor(color[0]);
      hist[0]->Draw("pesame");
    }
    legend->Draw();
  }

  //****************************
  //Add CMS and Lumi Labels
  //****************************
  // lumi_13TeV = "42 pb^{-1}";
  lumi_13TeV = "2.1 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  CMS_lumi(pad1,4,0);

  cv->cd();
  cv->Update();


  TPad *pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
    
  for (int b=0; b<histDataOverMC->GetXaxis()->GetNbins()+2; ++b) {
    double data = 0;
    if (hasData) {
      data = hist[0]->GetBinContent(b);
    }
    double MC = 0;
    double MCErrSqr = 0;
    if (hasData) {
      for (uint i = 1 ; i < hist.size(); ++i) {
	MC += hist[i]->GetBinContent(b);
	MCErrSqr += pow(hist[i]->GetBinError(b),2);
      }
    } else {
      MC = 1;
    }
      
    if (MC > 0) {
      histDataOverMC->SetBinContent(b, data / MC);
      histDataOverMC->SetBinError(b, (data / MC)*sqrt(1/data + MCErrSqr/pow(MC,2) ));
    } else {
      histDataOverMC->SetBinContent(b, 0);
      histDataOverMC->SetBinError(b, 0);
    }
    //cout << "bin " << b << " : " << histDataOverMC->GetBinContent(b) << " " << histDataOverMC->GetBinError(b) << "\n";
  }

  histDataOverMC->GetYaxis()->SetTitle("Data/MC");
  histDataOverMC->GetYaxis()->SetNdivisions(306);
  histDataOverMC->GetYaxis()->SetTitleSize(0.10);
  histDataOverMC->GetYaxis()->SetTitleOffset(0.3);
  histDataOverMC->GetYaxis()->SetRangeUser(0.5,1.5);
  histDataOverMC->GetYaxis()->SetLabelSize(0.10);
  histDataOverMC->GetXaxis()->SetLabelSize(0.125);
  histDataOverMC->GetXaxis()->SetTitleSize(0.15);
  histDataOverMC->GetXaxis()->SetTitleOffset(1.0);
  histDataOverMC->SetLineColor(kBlack);
  histDataOverMC->SetMarkerStyle(20);      
  histDataOverMC->SetMarkerSize(1);
  histDataOverMC->SetStats(false);
  histDataOverMC->Draw("pe");
  
  pad1->SetLogy(false);
  cv->SaveAs(Form("Razor_TTBarDileptonCR_%s%s.png",varName.c_str(), label.c_str()));
  cv->SaveAs(Form("Razor_TTBarDileptonCR_%s%s.pdf",varName.c_str(), label.c_str()));
  
  pad1->SetLogy(true);
  cv->SaveAs(Form("Razor_TTBarDileptonCR_%s%s_Logy.png",varName.c_str(),label.c_str()));
  cv->SaveAs(Form("Razor_TTBarDileptonCR_%s%s_Logy.pdf",varName.c_str(),label.c_str()));


 

}

// void PlotDataAndStackedBkg( vector<TH1F*> hist , vector<string> processLabels, vector<int> color,  bool hasData, string varName, string label ) {

//   TCanvas *cv =0;
//   TLegend *legend = 0;

//   cv = new TCanvas("cv","cv", 800,700);
//   cv->SetHighLightColor(2);
//   cv->SetFillColor(0);
//   cv->SetBorderMode(0);
//   cv->SetBorderSize(2);
//   cv->SetLeftMargin(0.16);
//   cv->SetRightMargin(0.3);
//   cv->SetTopMargin(0.07);
//   cv->SetBottomMargin(0.12);
//   cv->SetFrameBorderMode(0);  

//   TPad *pad1 = new TPad("pad1","pad1", 0,0.25,1,1);
//   pad1->SetBottomMargin(0.0);
//   pad1->SetRightMargin(0.04);
//   pad1->Draw();
//   pad1->cd();

//   legend = new TLegend(0.60,0.50,0.90,0.84);
//   legend->SetTextSize(0.04);
//   legend->SetBorderSize(0);
//   legend->SetFillStyle(0);

//   THStack *stack = new THStack();
//   TH1F *histDataOverMC = (TH1F*)hist[0]->Clone("histDataOverMC");

//   if (hasData) {
//     for (int i = hist.size()-1; i >= 1; --i) {
//       hist[i]->SetFillColor(color[i]);
//       hist[i]->SetFillStyle(1001);
      
//       if ( hist[i]->Integral() > 0) {
//   	stack->Add(hist[i]);
//       }
//     }
//   } else {
//     for (int i = hist.size()-1; i >= 0; --i) {
//       hist[i]->SetFillColor(color[i]);
//       hist[i]->SetFillStyle(1001);
      
//       if ( hist[i]->Integral() > 0) {
//   	stack->Add(hist[i]);
//       }
//     }
//   }

//   for (uint i = 0 ; i < hist.size(); ++i) {
//     if (hasData && i==0) {
//       legend->AddEntry(hist[i],(processLabels[i]).c_str(), "LP");
//     } else {
//       legend->AddEntry(hist[i],(processLabels[i]).c_str(), "F");
//     }
//   }

//   if (stack->GetHists()->GetEntries() > 0) {
//     stack->Draw("hist");
//     stack->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stack->GetHists()->At(0)))->GetXaxis()->GetTitle());
//     stack->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stack->GetHists()->At(0)))->GetYaxis()->GetTitle());
//     stack->GetHistogram()->GetYaxis()->SetTitleOffset(0.8);
//     stack->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
//     stack->GetHistogram()->GetXaxis()->SetTitleSize(0.15);
//     stack->SetMaximum( 1.2* fmax( stack->GetMaximum(), hist[0]->GetMaximum()) );
//     stack->SetMinimum( 0.1 );

//     if (hasData) {
//       hist[0]->SetLineWidth(2);
//       hist[0]->SetLineColor(color[0]);
//       hist[0]->Draw("e1same");
//     }
//     legend->Draw();
//   }
//   cv->cd();
//   cv->Update();


//   TPad *pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
//   pad2->SetTopMargin(0.01);
//   pad2->SetBottomMargin(0.37);
//   pad2->SetRightMargin(0.04);
//   pad2->Draw();
//   pad2->cd();
    
//   for (int b=0; b<histDataOverMC->GetXaxis()->GetNbins()+2; ++b) {
//     double data = 0;
//     if (hasData) {
//       data = hist[0]->GetBinContent(b);
//     }
//     double MC = 0;
//     double MCErrSqr = 0;
//     if (hasData) {
//       for (uint i = 1 ; i < hist.size(); ++i) {
// 	MC += hist[i]->GetBinContent(b);
// 	MCErrSqr += pow(hist[i]->GetBinError(b),2);
//       }
//     } else {
//       MC = 1;
//     }
      
//     if (MC > 0) {
//       histDataOverMC->SetBinContent(b, data / MC);
//       histDataOverMC->SetBinError(b, (data / MC)*sqrt(1/data + MCErrSqr/pow(MC,2) ));
//     } else {
//       histDataOverMC->SetBinContent(b, 0);
//       histDataOverMC->SetBinError(b, 0);
//     }
//     //cout << "bin " << b << " : " << histDataOverMC->GetBinContent(b) << " " << histDataOverMC->GetBinError(b) << "\n";
//   }

//   histDataOverMC->GetYaxis()->SetTitle("Data/MC");
//   histDataOverMC->GetYaxis()->SetNdivisions(306);
//   histDataOverMC->GetYaxis()->SetTitleSize(0.10);
//   histDataOverMC->GetYaxis()->SetTitleOffset(0.3);
//   histDataOverMC->GetYaxis()->SetRangeUser(0.5,1.5);
//   histDataOverMC->GetYaxis()->SetLabelSize(0.10);
//   histDataOverMC->GetXaxis()->SetLabelSize(0.125);
//   histDataOverMC->GetXaxis()->SetTitleSize(0.15);
//   histDataOverMC->GetXaxis()->SetTitleOffset(1.0);
//   histDataOverMC->SetStats(false);
//   histDataOverMC->Draw("e1");
  
//   pad1->SetLogy(false);
//   cv->SaveAs(Form("Razor_TTBarDileptonCR_%s%s.gif",varName.c_str(), label.c_str()));
  
//   pad1->SetLogy(true);
//   cv->SaveAs(Form("Razor_TTBarDileptonCR_%s%s_Logy.gif",varName.c_str(),label.c_str()));

// }



//=== MAIN MACRO ================================================================================================= 


void RunSelectTTBarDileptonControlSample(  vector<string> datafiles, vector<vector<string> > bkgfiles, vector<string> bkgLabels, vector<int> bkgColors, double lumi, string option, int channelOption = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  float MR = 0;
  float Rsq = 0;
  float mll = 0;

  // TFile *file = TFile::Open("data.root", "UPDATE");
  // file->cd();
  
  // TTree *tree = new TTree("tree", "tree");
  // tree->Branch("MR",&MR,"MR/F");
  // tree->Branch("Rsq",&Rsq,"Rsq/F");
  // tree->Branch("mll",&mll,"mll/F");



  bool printdebug = false;

  TFile *pileupWeightFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/data/PileupReweight_Spring15MCTo2015Data.root", "READ");
  TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("PileupReweight");
  assert(pileupWeightHist);

  TFile *eleEffSFFile = TFile::Open("root://eoscms//eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_TightElectronSelectionEffDenominatorReco_2015Final_Golden.root","READ");
  TH2D *eleEffSFHist = (TH2D*)eleEffSFFile->Get("ScaleFactor_TightElectronSelectionEffDenominatorReco");
  assert(eleEffSFHist);

  TFile *muonEffSFFile = TFile::Open("root://eoscms//eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_TightMuonSelectionEffDenominatorReco_2015Final_Golden.root","READ");
  TH2D *muonEffSFHist = (TH2D*)muonEffSFFile->Get("ScaleFactor_TightMuonSelectionEffDenominatorReco");
  assert(muonEffSFHist);

  TFile *eleTriggerEffSFFile = TFile::Open("root://eoscms//eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2015Final_Golden.root","READ");
  TH2D *eleTriggerEffSFHist = (TH2D*)eleTriggerEffSFFile->Get("ScaleFactor_EleTriggerEleCombinedEffDenominatorTight");
  assert(eleTriggerEffSFHist);

  TFile *muonTriggerEffSFFile = TFile::Open("root://eoscms//eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/2015_Final_Golden_2093/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2015Final_Golden.root","READ");
  TH2D *muonTriggerEffSFHist = (TH2D*)muonTriggerEffSFFile->Get("ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight");
  assert(muonTriggerEffSFHist);


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  const int NMRBins = 7;
  const int NRsqBins = 7;
  double MRBins[NMRBins] = {300, 350, 400, 450, 500, 550, 700};
  double RsqBins[NRsqBins] = {0.15,0.175,0.20,0.225,0.25,0.30,1.5};

  vector<vector<string> > inputfiles;
  vector<string> processLabels;
  vector<int> color;

  inputfiles.push_back(datafiles);
  processLabels.push_back("Data");
  color.push_back(kBlack);
  
  assert(bkgfiles.size() == bkgLabels.size());
  assert(bkgfiles.size() == bkgColors.size());
  for (int i=0; i < int(bkgfiles.size()); ++i) {
    inputfiles.push_back(bkgfiles[i]);
    processLabels.push_back(bkgLabels[i]);
    color.push_back(bkgColors[i]);
  }

  vector<TH1D*> histMR;
  vector<TH1D*> histRsq;
  vector<TH1D*> histMuonPt;
  vector<TH1D*> histElectronPt;
  vector<TH1D*> histMuonEta;
  vector<TH1D*> histElectronEta;
  vector<TH1D*> histDileptonDeltaPhi;
  vector<TH1D*> histDileptonPt;
  vector<TH1D*> histDileptonEta;
  vector<TH1D*> histMET;
  vector<TH1D*> histDileptonMass;
  vector<TH1D*> histDileptonCharge;
  vector<TH1D*> histNJets40;
  vector<TH1D*> histNJets80;
  vector<TH1D*> histNBtags;
  vector<TH2F*> histMRVsRsq;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    histMR.push_back(new TH1D(Form("histMR_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", 40, 300, 2300));
    histRsq.push_back(new TH1D(Form("histRsq_%s",processLabels[i].c_str()), "; R^{2} ; Number of Events", 50, 0.0, 1.50));
    histMuonPt.push_back(new TH1D(Form("histMuonPt_%s",processLabels[i].c_str()), "; MuonPt [GeV/c] ; Number of Events", 80, 0, 400));
    histElectronPt.push_back(new TH1D(Form("histElectronPt_%s",processLabels[i].c_str()), "; ElectronPt [GeV/c] ; Number of Events", 80, 0, 400));
    histMuonEta.push_back(new TH1D(Form("histMuonEta_%s",processLabels[i].c_str()), "; Muon #eta ; Number of Events", 50, -2.4, 2.4));
    histElectronEta.push_back(new TH1D(Form("histElectronEta_%s",processLabels[i].c_str()), "; Electron #eta ; Number of Events", 50, -2.5, 2.5));
    histDileptonDeltaPhi.push_back(new TH1D(Form("histDileptonDeltaPhi_%s",processLabels[i].c_str()), "; DileptonDeltaPhi [GeV/c] ; Number of Events", 25, 0, 3.1416));
    histDileptonPt.push_back(new TH1D(Form("histDileptonPt_%s",processLabels[i].c_str()), "; DileptonPt [GeV/c] ; Number of Events", 20, 0, 500));
    histDileptonEta.push_back(new TH1D(Form("histDileptonEta_%s",processLabels[i].c_str()), "; DileptonEta [GeV/c] ; Number of Events", 25, -10, 10));
    histMET.push_back(new TH1D(Form("histMET_%s",processLabels[i].c_str()), "; MET [GeV] ; Number of Events", 25, 0, 500));
    histDileptonMass.push_back(new TH1D(Form("histDileptonMass_%s",processLabels[i].c_str()), "; M_{ll} [GeV/c^{2}]; Number of Events", 40, 0, 500));
    histDileptonCharge.push_back(new TH1D(Form("histDileptonCharge_%s",processLabels[i].c_str()), "; Charge; Number of Events", 5, -2.5, 2.5));
    histNJets40.push_back(new TH1D(Form("histNJets40_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 40); Number of Events", 10, -0.5, 9.5));
    histNJets80.push_back(new TH1D(Form("histNJets80_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 80); Number of Events", 10, -0.5, 9.5));
    histNBtags.push_back(new TH1D(Form("histNBtags_%s",processLabels[i].c_str()), "; Number of B tags; Number of Events", 5, -0.5, 4.5));
    histMR[i]->Sumw2();
    histRsq[i]->Sumw2();
    histMuonPt[i]->Sumw2();
    histElectronPt[i]->Sumw2();
    histMuonEta[i]->Sumw2();
    histElectronEta[i]->Sumw2();
    histDileptonDeltaPhi[i]->Sumw2();
    histDileptonPt[i]->Sumw2();
    histDileptonEta[i]->Sumw2();
    histDileptonMass[i]->Sumw2();
    histDileptonCharge[i]->Sumw2();
    histMET[i]->Sumw2();
    histNJets40[i]->Sumw2();
    histNBtags[i]->Sumw2();  

    histMRVsRsq.push_back(new TH2F(Form("histMRVsRsq_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; R^{2}; Number of Events", NMRBins-1, MRBins, NRsqBins-1, RsqBins));
    histMRVsRsq[i]->Sumw2();
  }
 
  double dataYield = 0;
  double MCYield = 0;


  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  for (uint i=0; i < inputfiles.size(); ++i) {

    //for duplicate event checking
    map<pair<uint,uint>, bool > processedRunEvents;

    for (uint j=0; j < inputfiles[i].size(); ++j) {
      ControlSampleEvents *events = new ControlSampleEvents;
      events->LoadTree(inputfiles[i][j].c_str(), ControlSampleEvents::kTreeType_Dilepton_Full);

      bool isData = false;
      if ( processLabels[i] == "Data") isData = true;
    
      cout << "process: " << processLabels[i] << " | file " << inputfiles[i][j] << " | Total Entries: " << events->tree_->GetEntries() << "\n";
      for(UInt_t ientry=0; ientry < events->tree_->GetEntries(); ientry++) {       	
	events->tree_->GetEntry(ientry);
      

	if (ientry % 100000 == 0) cout << "Event " << ientry << endl;      

	double puWeight = 1;      
	double weight = 1;
	if (!isData) {
	  //puWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(events->NPU_0));
	  //weight = lumi * events->weight * puWeight;
	  weight = lumi * events->weight;
	}

	//apply k-factor to ttjets
	//if (processLabels[i] == "TTJets") weight = weight * 1.656;
	//if (processLabels[i] == "WJets") weight = weight * 1.447;


	if (isnan(events->weight) || isinf(events->weight)) {
	  continue;
	  //cout << "...bad event: " << weight << " " << (l1+l2).M() << "\n";
	}


	//******************************
	//Trigger Selection
	//******************************
	bool passTrigger = false;

	//Use Single Lepton Triggers
	//Use Single Lepton Triggers
	if ( events->HLTDecision[2] || events->HLTDecision[7] || events->HLTDecision[12] || events->HLTDecision[11] || events->HLTDecision[15])  
	  passTrigger = true;

	if (isData) {
	  if ( events->HLTDecision[22] || events->HLTDecision[23] || events->HLTDecision[24] || 
	       events->HLTDecision[25] || events->HLTDecision[26] ||
	       events->HLTDecision[27] || events->HLTDecision[28] || events->HLTDecision[29]	  
	       ) passTrigger = true;
	} else {
	  if ( events->HLTDecision[18] || events->HLTDecision[19] || events->HLTDecision[20] || 
	       events->HLTDecision[21] ||
	       events->HLTDecision[28] || events->HLTDecision[29]	  
	       ) passTrigger = true;
	}


	// //MuEG Triggers: Mu17Ele8 , Mu8Ele17
	// if (events->HLTDecision[45] || events->HLTDecision[46]
	//     || events->HLTDecision[47] || events->HLTDecision[48]
	//     ) passTrigger = true;

	// DiMuon Triggers: Mu17Mu8 , Mu17TkMu8
	// if (events->HLTDecision[39] || events->HLTDecision[41] ) passTrigger = true;

	// //DiElectron Triggers:
	// if ( events->HLTDecision[28] 
	// 	  || events->HLTDecision[29]	
	// 	  ) passTrigger = true;

	// //Razor Hadronic Triggers
	// if (isData) {
	// 	if ( events->HLTDecision[132] || events->HLTDecision[133] ) passTrigger = true;
	// } else {
	// 	if ( events->HLTDecision[136] || events->HLTDecision[137] ) passTrigger = true;
	// }

	// //Razor Hadronic Triggers
	// if (isData) {
	// 	if ( events->HLTDecision[132] || events->HLTDecision[133] ) passTrigger = true;
	// } else {
	// 	if ( events->HLTDecision[136] || events->HLTDecision[137] ) passTrigger = true;
	// }

	if (!passTrigger) continue;

	//******************************
	//apply scale corrections
	//******************************
	TLorentzVector l1;
	TLorentzVector l2;
	l1.SetPtEtaPhiM( events->lep1.Pt() , events->lep1.Eta(), events->lep1.Phi(), events->lep1.M());
	l2.SetPtEtaPhiM( events->lep2.Pt() , events->lep2.Eta(), events->lep2.Phi(), events->lep2.M());

	//******************************
	//Selection Cuts 
	//******************************
	if (!( (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	       &&
	       (abs(events->lep2Type) == 11 || abs(events->lep2Type) == 13)
	       )
	    ) continue;


	//lepton selection
	if (! (l1.Pt() > 30 || l2.Pt() > 30)) continue;
	if (! (l1.Pt() > 20 && l2.Pt() > 20
	       && events->lep1PassTight && events->lep2PassTight)
	    ) continue;

  
	//dilepton mass cut
	//if ( (events->lep1+events->lep2).M() < 20) continue;
	if ( (l1+l2).M() < 50) continue;

	//Z-mass window cut
	if ( abs(events->lep1Type) == abs(events->lep2Type) 
	     && 
	     (l1+l2).M() > 76 && (l1+l2).M() < 106
	     ) continue;
 

	//Check for duplicate data events
	if (isData) {
	  //cout << "event: " << events->run << " " << events->event << "\n";
	  if(!(processedRunEvents.find(make_pair(events->run, events->event)) == processedRunEvents.end())) {
	    //cout << "Duplicate event: " << events->run << " " << events->event << "\n";
	    continue;
	  } else {
	    processedRunEvents[make_pair(events->run, events->event)] = true;
	    //cout << processedRunEvents.size() << "\n";
	  }
	}

	//Razor signal region cuts
	if (option == "topEnhanced" || option == "MR300Rsq0p15") {
	  //if (!(events->NJets40 >= 2 )) continue;
	  if (!(events->NBJetsMedium >= 1 )) continue;
	}
      
	if (option == "MR300Rsq0p15" ) {
	  if (!(events->NJets80 >= 2 && events->MR > 300 && events->Rsq > 0.15 )) continue;
	}
      
 
	//******************************
	//ChannelOptions
	//******************************
	//e-mu only
	if (channelOption == 0 &&
	    !((abs(events->lep1Type) == 11 && abs(events->lep2Type) == 13) ||
	      (abs(events->lep1Type) == 13 && abs(events->lep2Type) == 11))
	    ) continue;
      
	//ee or mumu only
	if (channelOption == 1 && 
	    !(abs(events->lep1Type) == abs(events->lep2Type))
	    ) continue;
      
	//ee only
	if (channelOption == 2 &&
	    !(abs(events->lep1Type) == abs(events->lep2Type) && abs(events->lep1Type) == 11)
	    ) continue;

	//mumu only
	if (channelOption == 3 &&
	    !(abs(events->lep1Type) == abs(events->lep2Type) && abs(events->lep1Type) == 13)
	    ) continue;     

	//ee+mm
	if (channelOption == 4 &&
	    !(abs(events->lep1Type) == abs(events->lep2Type))
	    ) continue;     

	//******************************
	//Apply Scale Factors
	//******************************
	if (!isData) {
	  double triggerEffScaleFactor = 1.0;
	
	  double leptonEffScaleFactor = 1.0;
	  // if (abs(events->lep1Type) == 11  ) {
	  //   leptonEffScaleFactor *= eleEffSFHist->GetBinContent( eleEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(l1.Pt(),199.9),20.01)),
	  // 							 eleEffSFHist->GetYaxis()->FindFixBin(fabs(l1.Eta()))
	  // 							 );	 
	  // }
	  // if (abs(events->lep2Type) == 11 ) {
	  //   leptonEffScaleFactor *= eleEffSFHist->GetBinContent(  eleEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(l2.Pt(),199.9),20.01)),
	  // 							  eleEffSFHist->GetYaxis()->FindFixBin(fabs(l2.Eta()))
	  // 							  );
	  // }
	  // if (abs(events->lep1Type) == 13) {
	  //   leptonEffScaleFactor *= muonEffSFHist->GetBinContent( muonEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(l1.Pt(),199.9),20.01)),
	  // 							  muonEffSFHist->GetYaxis()->FindFixBin(fabs(l1.Eta()))
	  // 							  );	 
	  // }
	  // if (abs(events->lep2Type) == 13) {
	  //   leptonEffScaleFactor *= muonEffSFHist->GetBinContent(  muonEffSFHist->GetXaxis()->FindFixBin(fmax(fmin(l2.Pt(),199.9),20.01)),
	  // 							   muonEffSFHist->GetYaxis()->FindFixBin(fabs(l2.Eta()))
	  // 							   );
	  // }
	

	// 	//b-tagging scale factors
	// 	double btagScaleFactor = 1.0;
	// 	double bjet1EventScaleFactor = 1.0;
	// 	double bjet2EventScaleFactor = 1.0;
	// 	if (events->bjet1.Pt() > 20) {
	// 	  double bjet1ScaleFactor = 0.938887 + 0.00017124 * events->bjet1.Pt() + (-2.76366e-07) * events->bjet1.Pt() * events->bjet1.Pt() ;
	// 	  double MCEff = 1.0;
	// 	  if (events->bjet1.Pt() < 50) MCEff = 0.65;
	// 	  else if (events->bjet1.Pt() < 80) MCEff = 0.70;
	// 	  else if (events->bjet1.Pt() < 120) MCEff = 0.73;
	// 	  else if (events->bjet1.Pt() < 210) MCEff = 0.73;
	// 	  else MCEff = 0.66;				 
	// 	  if (events->bjet1PassMedium) bjet1EventScaleFactor = bjet1ScaleFactor;
	// 	  else bjet1EventScaleFactor = ( 1/MCEff - bjet1ScaleFactor) / ( 1/MCEff - 1);
	// 	}
	// 	if (events->bjet2.Pt() > 20) {
	// 	  double bjet2ScaleFactor = 0.938887 + 0.00017124 * events->bjet2.Pt() + (-2.76366e-07) * events->bjet2.Pt() * events->bjet1.Pt() ;
	// 	  double MCEff = 1.0;
	// 	  if (events->bjet2.Pt() < 50) MCEff = 0.65;
	// 	  else if (events->bjet2.Pt() < 80) MCEff = 0.70;
	// 	  else if (events->bjet2.Pt() < 120) MCEff = 0.73;
	// 	  else if (events->bjet2.Pt() < 210) MCEff = 0.73;
	// 	  else MCEff = 0.66;				 
	// 	  if (events->bjet2PassMedium) bjet2EventScaleFactor = bjet2ScaleFactor;
	// 	  else bjet2EventScaleFactor = ( 1/MCEff - bjet2ScaleFactor) / ( 1/MCEff - 1);
	// 	}
	// 	btagScaleFactor = bjet1EventScaleFactor * bjet2EventScaleFactor;

	// 	// cout << events->NBJetsMedium << " : " << events->bjet1.Pt() << " " << events->bjet1PassMedium << " " << events->bjet2.Pt() << " " << events->bjet2PassMedium  
	// 	//      << " : " << bjet1EventScaleFactor << " " << bjet2EventScaleFactor << " " << btagScaleFactor << " "
	// 	//      <<  " \n";


	  // weight = weight * leptonEffScaleFactor;
	  // weight = weight * triggerEffScaleFactor;
	  // 	weight = weight * btagScaleFactor;

	  // 	if (processLabels[i] == "DY") {
	  // 	  weight *= DYScaleFactorsHist->GetBinContent( DYScaleFactorsHist->GetXaxis()->FindFixBin(fmin(fmax(events->MR,300.1),549.9)) ,  
	  // 						       DYScaleFactorsHist->GetYaxis()->FindFixBin(fmin(fmax(events->Rsq,0.0501),1.499)) );	 
	  // 	}
	  
	}


	//******************************
	//Fill histograms
	//******************************
	if (isData) {
	  dataYield += 1.0;
	  histDileptonMass[i]->Fill((l1+l2).M());      
	  histMET[i]->Fill(fmin(events->METnoHF,399.9));


	  if (abs(events->lep1Type) == 13) {
	    histMuonPt[i]->Fill(l1.Pt());
	    histElectronPt[i]->Fill(l2.Pt());
	    histMuonEta[i]->Fill(l1.Eta());
	    histElectronEta[i]->Fill(l2.Eta());
	  } else {
	    histElectronPt[i]->Fill(l1.Pt());
	    histMuonPt[i]->Fill(l2.Pt());
	    histElectronEta[i]->Fill(l1.Eta());
	    histMuonEta[i]->Fill(l2.Eta());
	  }
	  histDileptonDeltaPhi[i]->Fill( acos(cos(l1.Phi() - l2.Phi())) );
	  histDileptonPt[i]->Fill( (l1+l2).Pt() );
	  histDileptonEta[i]->Fill( (l1+l2).Eta() );
	  histDileptonCharge[i]->Fill( events->lep1Type/abs(events->lep1Type) + events->lep2Type/abs(events->lep2Type));
	  histNJets40[i]->Fill( events->NJets40 );
	  histNJets80[i]->Fill( events->NJets80 );
	  histNBtags[i]->Fill( events->NBJetsLoose );

	  if (events->MR > 0) {
	    histMR[i]->Fill(events->MR);
	    histRsq[i]->Fill(events->Rsq);
	  }

	  histMRVsRsq[i]->Fill(events->MR,events->Rsq);

	} else {
	  MCYield += weight;
	  histDileptonMass[i]->Fill((l1+l2).M(), weight );      
	  histMET[i]->Fill(fmin(events->METnoHF,399.9), weight );
	  if (abs(events->lep1Type) == 13) {
	    histMuonPt[i]->Fill(l1.Pt(), weight);
	    histElectronPt[i]->Fill(l2.Pt(), weight);
	    histMuonEta[i]->Fill(l1.Eta(), weight);
	    histElectronEta[i]->Fill(l2.Eta(), weight);
	  } else {
	    histElectronPt[i]->Fill(l1.Pt(), weight);
	    histMuonPt[i]->Fill(l2.Pt(), weight);
	    histElectronEta[i]->Fill(l1.Eta(), weight);
	    histMuonEta[i]->Fill(l2.Eta(), weight);
	  }
	  histDileptonDeltaPhi[i]->Fill( acos(cos(l1.Phi() - l2.Phi())) , weight);
	  histDileptonPt[i]->Fill( (l1+l2).Pt() , weight);
	  histDileptonEta[i]->Fill( (l1+l2).Eta() , weight);

	  histDileptonCharge[i]->Fill( events->lep1Type/abs(events->lep1Type) + events->lep2Type/abs(events->lep2Type) , weight);
	  histNJets40[i]->Fill( events->NJets40 , weight);
	  histNJets80[i]->Fill( events->NJets80 , weight);
	  histNBtags[i]->Fill( events->NBJetsLoose , weight);
	
	  if (events->MR > 0) {
	    histMR[i]->Fill(events->MR, weight );
	    histRsq[i]->Fill(events->Rsq, weight );
	  }

	  histMRVsRsq[i]->Fill(events->MR,events->Rsq, weight);

	}
      }
    }
  }


  //--------------------------------------------------------------------------------------------------------------
  // Compute Expected Statistical Uncertainty
  //==============================================================================================================
  // TH2F *statUnc_MRVsRsq = (TH2F*)(histMRVsRsq[0]->Clone("statUnc_MRVsRsq"));
  // for (int i=0; i<statUnc_MRVsRsq->GetXaxis()->GetNbins()+1;i++) {
  //   for (int j=0; j<statUnc_MRVsRsq->GetYaxis()->GetNbins()+1;j++) {
  //     if (statUnc_MRVsRsq->GetBinContent(i,j) > 1) {
  // 	statUnc_MRVsRsq->SetBinContent(i,j,1.0 / sqrt(statUnc_MRVsRsq->GetBinContent(i,j)));  	
  //     } else {
  // 	statUnc_MRVsRsq->SetBinContent(i,j,0);
  //     }
  //   }
  // }

  //--------------------------------------------------------------------------------------------------------------
  // Subtract Non TTBar Bkg
  //==============================================================================================================
  TH2F *DataMinusBkg_MRVsRsq = (TH2F*)(histMRVsRsq[0]->Clone("DataMinusBkg_MRVsRsq"));
  TH2F *MCToDataScaleFactor_MRVsRsq = (TH2F*)(histMRVsRsq[0]->Clone("MCToDataScaleFactor_MRVsRsq"));

  for (int i=0; i<DataMinusBkg_MRVsRsq->GetXaxis()->GetNbins()+1;i++) {
    for (int j=0; j<DataMinusBkg_MRVsRsq->GetYaxis()->GetNbins()+1;j++) {
      
      double data = histMRVsRsq[0]->GetBinContent(i,j);
      double mc = 0; 
      double mc_StatErr = 0; 
      double bkg = 0;
      double bkg_StatErrSqr = 0;
      double bkg_SysErrSqr = 0;

      for (uint k=1; k < inputfiles.size(); ++k) {

	if (processLabels[k] == "TTJets") {
	  mc = histMRVsRsq[k]->GetBinContent(i,j);
	  mc_StatErr = sqrt(histMRVsRsq[k]->GetBinError(i,j));
	  continue;
	}

	double systematicUncertainty = 0;
	if (processLabels[k] == "VV") systematicUncertainty = 0.2;
	if (processLabels[k] == "SingleTop") systematicUncertainty = 0.2;
	if (processLabels[k] == "TT+V") systematicUncertainty = 0.2;
	if (processLabels[k] == "DY") systematicUncertainty = 0.2;
 
	bkg += histMRVsRsq[k]->GetBinContent(i,j);
	bkg_StatErrSqr += pow(histMRVsRsq[k]->GetBinError(i,j),2);
	bkg_SysErrSqr += pow( histMRVsRsq[k]->GetBinContent(i,j) * systematicUncertainty, 2);
      }

      DataMinusBkg_MRVsRsq->SetBinContent(i,j, data - bkg );
      double dataMinusBkgTotalErr = sqrt(data + bkg_StatErrSqr + bkg_SysErrSqr);
      DataMinusBkg_MRVsRsq->SetBinError(i,j, dataMinusBkgTotalErr );


      cout << "Bin " << DataMinusBkg_MRVsRsq->GetXaxis()->GetBinCenter(i) << " " << DataMinusBkg_MRVsRsq->GetYaxis()->GetBinCenter(j) << " : "
	   << data << " " << mc << " " << bkg << " " << mc_StatErr << " " << bkg_StatErrSqr << " " << bkg_SysErrSqr << "\n";

      MCToDataScaleFactor_MRVsRsq->SetBinContent(i,j, (data - bkg)/mc );
      MCToDataScaleFactor_MRVsRsq->SetBinError(i,j, ((data - bkg)/mc)*sqrt( pow(mc_StatErr/mc,2) + pow(dataMinusBkgTotalErr/(data-bkg),2)) );

    }
  }

  for (int i=0; i<DataMinusBkg_MRVsRsq->GetXaxis()->GetNbins()+1;i++) {
    for (int j=0; j<DataMinusBkg_MRVsRsq->GetYaxis()->GetNbins()+1;j++) {
      cout << "Bin " << DataMinusBkg_MRVsRsq->GetXaxis()->GetBinCenter(i) << " " << DataMinusBkg_MRVsRsq->GetYaxis()->GetBinCenter(j) << " : " << MCToDataScaleFactor_MRVsRsq->GetBinContent(i,j) << " +/- " << MCToDataScaleFactor_MRVsRsq->GetBinError(i,j) << "\n";	
    }
  }



 
  //--------------------------------------------------------------------------------------------------------------
  // Make Plots
  //==============================================================================================================


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;
  TLegend *legend = 0;

  //*******************************************************************************************
  //MR
  //*******************************************************************************************
  PlotDataAndStackedBkg( histMR, processLabels, color, true, "MR", Label);
  PlotDataAndStackedBkg( histRsq, processLabels, color, true, "Rsq", Label);
  PlotDataAndStackedBkg( histDileptonMass, processLabels, color, true, "DileptonMass", Label);
  PlotDataAndStackedBkg( histMuonPt, processLabels, color, true, "MuonPt", Label);
  PlotDataAndStackedBkg( histElectronPt, processLabels, color, true, "ElectronPt", Label);
  PlotDataAndStackedBkg( histMuonEta, processLabels, color, true, "MuonEta", Label);
  PlotDataAndStackedBkg( histElectronEta, processLabels, color, true, "ElectronEta", Label);
  PlotDataAndStackedBkg( histDileptonDeltaPhi, processLabels, color, true, "DileptonDeltaPhi", Label);
  PlotDataAndStackedBkg( histDileptonPt, processLabels, color, true, "DileptonPt", Label);
  PlotDataAndStackedBkg( histDileptonEta, processLabels, color, true, "DileptonEta", Label);
  PlotDataAndStackedBkg( histMET, processLabels, color, true, "MET", Label);
  PlotDataAndStackedBkg( histDileptonCharge, processLabels, color, true, "DileptonCharge", Label);
  PlotDataAndStackedBkg( histNJets40, processLabels, color, true, "NJets40", Label);
  PlotDataAndStackedBkg( histNJets80, processLabels, color, true, "NJets80", Label);
  PlotDataAndStackedBkg( histNBtags, processLabels, color, true, "NBtags", Label);
 

  //--------------------------------------------------------------------------------------------------------------
  // Tables
  //==============================================================================================================
  cout << "For Luminosity = " << lumi << " pb^-1\n";
  // cout << "Yields : MR > 300 && Rsq > 0.1\n";
  //cout << "TTJets: " << 

  cout << "Selected Event Yield \n";
  cout << "Data: " << dataYield << "\n";
  cout << "MC: " << MCYield << "\n";

  // //--------------------------------------------------------------------------------------------------------------
  // // Output
  // //==============================================================================================================
  TFile *file = TFile::Open(("TTBarDileptonControlRegionPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDileptonMass[i], Form("histDileptonMass_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDileptonCharge[i], Form("histDileptonCharge_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNJets40[i], Form("histNJets40_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNBtags[i], Form("histNBtags_%s",processLabels[i].c_str()), "WriteDelete");
  }


  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDileptonMass[i], Form("histDileptonMass_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDileptonCharge[i], Form("histDileptonCharge_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNJets40[i], Form("histNJets40_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histNBtags[i], Form("histNBtags_%s",processLabels[i].c_str()), "WriteDelete");
  }
 
  
  for(int i=0; i<int(histMRVsRsq.size()); i++) {
    file->WriteTObject(histMRVsRsq[i], Form("histMRVsRsq_%s",processLabels[i].c_str()), "WriteDelete");
  }
  //file->WriteTObject(statUnc_MRVsRsq,"statUnc_MRVsRsq_TTJets","WriteDelete");

  file->WriteTObject(DataMinusBkg_MRVsRsq, "DataMinusBkg_MRVsRsq", "WriteDelete");
  file->WriteTObject(MCToDataScaleFactor_MRVsRsq, "MCToDataScaleFactor_MRVsRsq", "WriteDelete");

  file->Close();
  delete file;       

}







void SelectTTBarDileptonControlSample( int option = 0) {

  vector<string> datafiles;
  vector<vector<string> > bkgfiles;
  vector<string> processLabels;
  vector<int> colors;

  string datafile = "";


  //No Skims  
  if (option == 0 || option == 10 || option == 20 || option == 3 || option == 4 || option == 13 || option == 14 || option == 23 || option == 24) {
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleElectron_Run2015D_GoodLumiGolden.root");
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleMuon_Run2015D_GoodLumiGolden.root");
  }
  if (option == 1|| option == 11 || option == 21) {
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleElectron_Run2015D_GoodLumiGolden.root");    
  }
  if (option == 2|| option == 12 || option == 22) {
    datafiles.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleMuon_Run2015D_GoodLumiGolden.root");
  }

  vector<string> bkgfiles_ttbar;
  vector<string> bkgfiles_wjets;
  vector<string> bkgfiles_singletop;
  vector<string> bkgfiles_dy;  
  vector<string> bkgfiles_vv; 
  vector<string> bkgfiles_ttv;

  //bkgfiles_ttbar.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
  //bkgfiles_ttbar.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root");
  bkgfiles_ttbar.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");

  bkgfiles_wjets.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root"); 
  bkgfiles_singletop.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleTop_1pb_weighted.root");  
  bkgfiles_dy.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root");
  //bkgfiles_dy.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root");
  bkgfiles_vv.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_VV_1pb_weighted.root");
  bkgfiles_ttv.push_back("/afs/cern.ch/user/s/sixie/eos2/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_TTV_1pb_weighted.root");    



  bkgfiles.push_back(bkgfiles_ttbar);
  bkgfiles.push_back(bkgfiles_dy);
  bkgfiles.push_back(bkgfiles_vv);
  bkgfiles.push_back(bkgfiles_singletop);
  bkgfiles.push_back(bkgfiles_wjets);
  // bkgfiles.push_back(bkgfiles_ttv);

  processLabels.push_back("TTJets");  
  processLabels.push_back("DY");
  processLabels.push_back("WW");
  processLabels.push_back("SingleTop");
  processLabels.push_back("WJets");
  //processLabels.push_back("TT+V");

  colors.push_back(kAzure+10);
  colors.push_back(kGreen+2);
  colors.push_back(kGray);
  colors.push_back(kBlue);
   colors.push_back(kRed);
  //colors.push_back(kRed);
 
   double lumi = 2185;

  //*********************************************************************
  //E-Mu Control Region
  //*********************************************************************
  // RunSelectTTBarDileptonControlSample(datafiles, bkgfiles, processLabels, colors, 19780,"Inclusive",0,"Inclusive_emu");
  //RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, 19780,"Met>40",0,"MetGreaterThan40_emu");
  //RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, 19780,"topEnhanced",0,"TopEnhanced_emu");
  //RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, 19780,"TwoJet80",0,"TwoJet80_emu");
  
  if (option == 0) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"MR300Rsq0p15",0,"MR300Rsq0p15_emu");
  if (option == 10) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"topEnhanced",0,"TopEnhanced_emu");
  if (option == 20) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"Inclusive",0,"Inclusive_emu");

  //*********************************************************************
  //E-E Control Region
  //*********************************************************************
  //RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, 19789,"Inclusive",2,"Inclusive_ee");
  //RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, 19789,"topEnhanced",2,"TopEnhanced_ee");
  //RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, 19789,"TwoJet80",2,"TwoJet80_ee");
  if (option == 1) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"MR300Rsq0p15",2,"MR300Rsq0p15_ee");
  if (option == 11) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"topEnhanced",2,"TopEnhanced_ee");
  if (option == 21) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"Inclusive",2,"Inclusive_ee");

  //*********************************************************************
  //MuMu Control Region
  //*********************************************************************
  //RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, 19751,"Inclusive",3,"Inclusive_mumu");
  //RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, 19751,"topEnhanced",3,"TopEnhanced_mumu");
  //RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, 19751,"TwoJet80",3,"TwoJet80_mumu");
  if (option == 2) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"MR300Rsq0p15",3,"MR300Rsq0p15_mumu");
  if (option == 12) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"topEnhanced",3,"TopEnhanced_mumu");
  if (option == 22) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"Inclusive",3,"Inclusive_mumu");

  //*********************************************************************
  //EE + MM 
  //*********************************************************************
  if (option == 24) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"Inclusive",4,"Inclusive_eemumu");
  if (option == 14) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"topEnhanced",4,"TopEnhanced_eemumu");

  //*********************************************************************
  //All final states Control Region
  //*********************************************************************
  if (option == 23) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"Inclusive",-1,"Inclusive_all");
  if (option == 13) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"topEnhanced",-1,"TopEnhanced_all");





}


//**********************************
//E-Mu Yields :
//**********************************
//Inclusive Single Lep Triggers ( 30 - 20 )
// Data: 164
// MC: 164.556 (WW only)
// MC: 166.709 (WW+WZ)
// MC: 167.135 (WW+WZ+ZZ)

//MuEG Triggers ( 30 - 20 )
// Data: 148
// MC: 156.16

//Met>40
//Data: 44398
// MC: 44765.7

//TwoJet80
// Data: 4579
// MC: 5488.65

//MR300Rsq0p15
// Data: 579
// MC: 678.504

//**********************
//Mu-Mu Yields
//**********************
//Inclusive
//Data: 7.67962e+06
//MC: 7.87451e+06

//TwoJet80
//Data: 5731
//MC: 6126.44

//**********************
//E-E Yields
//**********************
// Inclusive
//Data: 5.42906e+06
//MC: 5.56361e+06

//TwoJet80
//Data: 263
//MC: 318



// Z->MM
// Data: 7.47318e+06
// MC: 7.62973e+06 Powheg
// MC: 7.62596e+06 Madgraph

// Z->EE
//Data: 5.27829e+06
//MC: 5.37693e+06 Madgraph


