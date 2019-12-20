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
#include <TH2Poly.h>                
#include <TCanvas.h>                
#include <TLegend.h> 
#include <THStack.h> 
#include "TApplication.h"

#include "include/ControlSampleEvents.h"
#include "macros/tdrstyle.C"
#include "macros/CMS_lumi.C"

#endif

void PlotDataAndStackedBkg( vector<TH1D*> hist , vector<string> processLabels, vector<int> color,  bool hasData, string varName, string label ) {
  string plotsDir = "Plots/Razor2016/TTBarDileptonControlRegionMacro";

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

  legend = new TLegend(0.75,0.60,0.90,0.85);
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
    //stack->SetMaximum( 1.2* fmax( stack->GetMaximum(), hist[0]->GetMaximum()) );
    stack->SetMaximum( 10* fmax( stack->GetMaximum(), hist[0]->GetMaximum()) );
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
  lumi_13TeV = "26.4 fb^{-1}";
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
  histDataOverMC->GetYaxis()->SetRangeUser(0.0,3.0);
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
  cv->SaveAs(Form("%s/Razor_TTBarDileptonCrossCheckRegion_%s%s.png", 
              plotsDir.c_str(), varName.c_str(), label.c_str()));
  cv->SaveAs(Form("%s/Razor_TTBarDileptonCrossCheckRegion_%s%s.pdf", 
              plotsDir.c_str(), varName.c_str(), label.c_str()));
  
  pad1->SetLogy(true);
  cv->SaveAs(Form("%s/Razor_TTBarDileptonCrossCheckRegion_%s%s_Logy.png",
              plotsDir.c_str(), varName.c_str(),label.c_str()));
  cv->SaveAs(Form("%s/Razor_TTBarDileptonCrossCheckRegion_%s%s_Logy.pdf",
              plotsDir.c_str(), varName.c_str(),label.c_str()));

}



//=== MAIN MACRO ================================================================================================= 


void RunSelectTTBarDileptonControlSample(  vector<string> datafiles, vector<vector<string> > bkgfiles, vector<string> bkgLabels, vector<int> bkgColors, double lumi, string option, int channelOption = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;

  char* cmsswPath;  
  cmsswPath = getenv("CMSSW_BASE");
  if (cmsswPath == NULL) {
    cout << "Warning: CMSSW_BASE not detected. Exiting..." << endl;
    return;
  }

  TFile *RazorScaleFactorFile = 0;
  TH1F *NJetsSFHist = 0;
  RazorScaleFactorFile = TFile::Open(Form("%s/src/RazorAnalyzer/data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Razor2016.root",cmsswPath), "READ");
  TH2Poly *TTBarSFHist = (TH2Poly*)RazorScaleFactorFile->Get("TTJetsScaleFactors");
  assert(TTBarSFHist);
  TH2Poly *WJetsSFHist = (TH2Poly*)RazorScaleFactorFile->Get("WJetsScaleFactors");
  assert(WJetsSFHist);

  TFile *RazorNJetsScaleFactorFile = TFile::Open(Form("%s/src/RazorAnalyzer/data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_Razor2016.root",cmsswPath), "READ");
  NJetsSFHist = (TH1F*)RazorNJetsScaleFactorFile->Get("NJetsCorrectionScaleFactors");
  assert(NJetsSFHist);

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  const int NMRBinsMultiJet = 6;
  const int NRsqBinsMultiJet = 4;
  double MRBinsMultiJet[NMRBinsMultiJet] = {300, 400, 500, 700, 900, 4000};
  double RsqBinsMultiJet[NRsqBinsMultiJet] = {0.15, 0.20, 0.30, 1.5};
  const int NMRBinsInclusive = 6;
  const int NRsqBinsInclusive = 6;
  double MRBinsInclusive[NMRBinsInclusive] = {300, 400, 500, 700, 900, 4000};
  double RsqBinsInclusive[NRsqBinsInclusive] = {0.15, 0.20, 0.25, 0.30, 0.41, 1.5};
  
  const int NMRBins1D = 8;
  const int NRsqBins1D = 7;
  double MRBins1D[NMRBins1D] = {300, 400, 500, 600, 700, 900, 1200, 4000};
  double RsqBins1D[NRsqBins1D] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 1.5};

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
  vector<TH1D*> histLeptonPt;
  vector<TH1D*> histLeptonEta;
  vector<TH1D*> histLep1MT;
  vector<TH1D*> histMET;
  vector<TH1D*> histNJets40;
  vector<TH1D*> histNJets80;
  vector<TH1D*> histNBtags;
  vector<TH1D*> histMRRsqUnrolled;
  vector<TH2F*> histMRVsRsq;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    histMR.push_back(new TH1D(Form("histMR_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", NMRBins1D-1, MRBins1D));
    histRsq.push_back(new TH1D(Form("histRsq_%s",processLabels[i].c_str()), "; R^{2} ; Number of Events", NRsqBins1D-1, RsqBins1D));

    if (option == "4JetMR300Rsq0p15") {
      histMRRsqUnrolled.push_back(new TH1D(Form("histMRRsqUnrolled_%s",processLabels[i].c_str()), "; Bin Number ; Number of Events", (NMRBinsMultiJet-1)*(NRsqBinsMultiJet-1), 0, (NMRBinsMultiJet-1)*(NRsqBinsMultiJet-1)));
    } else {
      histMRRsqUnrolled.push_back(new TH1D(Form("histMRRsqUnrolled_%s",processLabels[i].c_str()), "; Bin Number ; Number of Events", (NMRBinsInclusive-1)*(NRsqBinsInclusive-1), 0, (NMRBinsInclusive-1)*(NRsqBinsInclusive-1)));
    }
    histLeptonPt.push_back(new TH1D(Form("histLeptonPt_%s",processLabels[i].c_str()), "; LeptonPt [GeV/c] ; Number of Events", 80, 0, 400));    
    histLeptonEta.push_back(new TH1D(Form("histLeptonEta_%s",processLabels[i].c_str()), "; Lepton #eta ; Number of Events", 50, -2.4, 2.4));
    histLep1MT.push_back(new TH1D(Form("histLep1MT_%s",processLabels[i].c_str()), "; Lep1MT [GeV/c] ; Number of Events", 30, 0, 600));
    histMET.push_back(new TH1D(Form("histMET_%s",processLabels[i].c_str()), "; MET [GeV] ; Number of Events", 25, 0, 500));
    histNJets40.push_back(new TH1D(Form("histNJets40_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 40); Number of Events", 10, -0.5, 9.5));
    histNJets80.push_back(new TH1D(Form("histNJets80_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 80); Number of Events", 10, -0.5, 9.5));
    histNBtags.push_back(new TH1D(Form("histNBtags_%s",processLabels[i].c_str()), "; Number of B tags; Number of Events", 5, -0.5, 4.5));
    histMR[i]->Sumw2();
    histRsq[i]->Sumw2();
    histMRRsqUnrolled[i]->Sumw2();
    histLeptonPt[i]->Sumw2();
    histLeptonEta[i]->Sumw2();
    histLep1MT[i]->Sumw2();
    histMET[i]->Sumw2();
    histNJets40[i]->Sumw2();
    histNJets80[i]->Sumw2();
    histNBtags[i]->Sumw2();  

    if (option == "4JetMR300Rsq0p15") { 
      histMRVsRsq.push_back(new TH2F(Form("histMRVsRsq_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; R^{2}; Number of Events", NMRBinsMultiJet-1, MRBinsMultiJet, NRsqBinsMultiJet-1, RsqBinsMultiJet));
    } else {
      histMRVsRsq.push_back(new TH2F(Form("histMRVsRsq_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; R^{2}; Number of Events", NMRBinsInclusive-1, MRBinsInclusive, NRsqBinsInclusive-1, RsqBinsInclusive));
    }
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
	  weight = lumi * events->weight;
	}

	if (isnan(events->weight) || isinf(events->weight)) {
	  cout << "...bad event: " << weight << "\n";
	  continue;
	}


	//******************************
	//Trigger Selection
	//******************************
	bool passTrigger = false;

	//Use Single Lepton Triggers
        if (events->HLTDecision[4] || events->HLTDecision[6] || events->HLTDecision[15] || 
            events->HLTDecision[22] || events->HLTDecision[33] || events->HLTDecision[35] || 
            events->HLTDecision[36] || events->HLTDecision[41] || events->HLTDecision[42]) passTrigger = true;

	if (!isData) {
            passTrigger = true;
	}

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
	if (!(events->lep1PassTight && events->lep2PassTight && l1.Pt() > 30 && l2.Pt() > 30)) continue;
  
	//dilepton mass cut
        if ( events->mll < 20 ) continue;
	if ( events->MET < 40 ) continue;

	//Z-mass window cut
	if ( abs(events->lep1Type) == abs(events->lep2Type) 
	     && 
	     (l1+l2).M() > 76 && (l1+l2).M() < 106
	     ) continue;
 

	//Check for duplicate data events
	if (isData) {
	  if(!(processedRunEvents.find(make_pair(events->run, events->event)) == processedRunEvents.end())) {
	    continue;
	  } else {
	    processedRunEvents[make_pair(events->run, events->event)] = true;
	  }
	}

	//Razor signal region cuts
	if (option == "topEnhanced" || option == "MR300Rsq0p15") {
	  if (!(events->NBJetsMedium >= 1 )) continue;
	}
      
	if (option == "MR300Rsq0p15" ) {
	  if (!(events->MR > 300 && events->Rsq > 0.15 )) continue;
	}

	if (option == "4JetMR300Rsq0p15" ) {
	  if (!(events->NJets40 >= 4 && events->MR > 300 && 
                      events->Rsq > 0.15 && events->NBJetsMedium >= 1)) continue;	  
	}

        if ( option == "2JetMR300Rsq0p15" ) {
          if (!(events->NJets40 >= 2 && events->NJets40 < 4 && events->MR > 300 && 
                      events->Rsq > 0.15 && events->NBJetsMedium >= 1 )) continue;
        }
      
	//MET Filters
	if (!(events->Flag_HBHENoiseFilter && events->Flag_goodVertices && events->Flag_eeBadScFilter)) continue;

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
	  double razorSF = 1.0;
	  if (processLabels[i] == "TTJets") {

	    razorSF = TTBarSFHist->GetBinContent( TTBarSFHist->FindBin(fmin(fmax(events->MR,300.1),3999), fmin(fmax(events->Rsq,0.151),1.49)));
            razorSF *= NJetsSFHist->GetBinContent( NJetsSFHist->FindBin(events->NJets40) );
	  }
	  if (processLabels[i] == "WJets") {
	    razorSF = WJetsSFHist->GetBinContent( WJetsSFHist->FindBin(fmin(fmax(events->MR,300.1),3999), fmin(fmax(events->Rsq,0.151),1.49)));
            razorSF *= NJetsSFHist->GetBinContent( NJetsSFHist->FindBin(events->NJets40) );
	  }
	  weight *= razorSF;

          // top pt reweighting
          if ( events->topPtWeight > 0 ) {
              //weight *= events->topPtWeight;
          }
	}


	//******************************
	//Make lepton 2 disappear
	//******************************		
	double lep1MT = sqrt(events->lep1.M2() + 2*events->MET*events->lep1.Pt()*(1 - cos( acos(cos(events->METPhi - events->lep1.Phi())))));	
	int MRBin = histMRVsRsq[i]->GetXaxis()->FindFixBin( events->MR );
	int RsqBin = histMRVsRsq[i]->GetYaxis()->FindFixBin( events->Rsq );
	int MRRsqBin = 0;
	if (option == "4JetMR300Rsq0p15") {
	  MRRsqBin = (NRsqBinsMultiJet-1)*(MRBin-1) + RsqBin-1;
	} else {
	  MRRsqBin = (NRsqBinsInclusive-1)*(MRBin-1) + RsqBin-1;
	}

	if (isData) {
	  histLeptonPt[i]->Fill(events->lep1.Pt(), 0.5);
	  histLeptonEta[i]->Fill(events->lep1.Eta(), 0.5);
	  histMET[i]->Fill(events->MET, 0.5);
	  histLep1MT[i]->Fill(lep1MT, 0.5);

	  if (lep1MT > 120) {
            dataYield += 0.5;
	    histNJets40[i]->Fill( events->NJets40, 0.5 );
	    histNJets80[i]->Fill( events->NJets80, 0.5 );
	    histNBtags[i]->Fill( events->NBJetsMedium, 0.5 );	    	   
	    histMR[i]->Fill(events->MR, 0.5);
	    histRsq[i]->Fill(events->Rsq, 0.5);	   
	    histMRVsRsq[i]->Fill(events->MR,events->Rsq, 0.5);
	    histMRRsqUnrolled[i]->Fill( MRRsqBin+0.5, 0.5);
	  }

	} else {
	  histLeptonPt[i]->Fill(events->lep1.Pt(), 0.5*weight);
	  histLeptonEta[i]->Fill(events->lep1.Eta(), 0.5*weight);
	  histMET[i]->Fill(events->MET, 0.5*weight);
	  histLep1MT[i]->Fill(lep1MT, 0.5*weight);

	  if (lep1MT > 120) {
            MCYield += weight * 0.5;
	    histNJets40[i]->Fill( events->NJets40, 0.5*weight );
	    histNJets80[i]->Fill( events->NJets80, 0.5*weight );
	    histNBtags[i]->Fill( events->NBJetsMedium, 0.5*weight );	    	   
	    histMR[i]->Fill(events->MR, 0.5*weight);
	    histRsq[i]->Fill(events->Rsq, 0.5*weight);	   
	    histMRVsRsq[i]->Fill(events->MR,events->Rsq, 0.5*weight);
	    histMRRsqUnrolled[i]->Fill( MRRsqBin+0.5, 0.5*weight);
	  }	
	}

	//******************************
	//Make lepton 1 disappear
	//******************************		
	double lep2MT = sqrt(events->lep2.M2() + 2*events->MET*events->lep2.Pt()*(1 - cos( acos(cos(events->METPhi - events->lep2.Phi())))));

	if (isData) {
	  histLeptonPt[i]->Fill(events->lep2.Pt(), 0.5);
	  histLeptonEta[i]->Fill(events->lep2.Eta(), 0.5);
	  histMET[i]->Fill(events->MET, 0.5);
	  histLep1MT[i]->Fill(lep2MT, 0.5);

	  if (lep2MT > 120) {
	    dataYield += 0.5;
	    histNJets40[i]->Fill( events->NJets40, 0.5 );
	    histNJets80[i]->Fill( events->NJets80, 0.5 );
	    histNBtags[i]->Fill( events->NBJetsMedium, 0.5 );	    	   
	    histMR[i]->Fill(events->MR, 0.5);
	    histRsq[i]->Fill(events->Rsq, 0.5);	   
	    histMRVsRsq[i]->Fill(events->MR,events->Rsq, 0.5);
	    histMRRsqUnrolled[i]->Fill( MRRsqBin+0.5, 0.5);
	  }

	} else {
	  histLeptonPt[i]->Fill(events->lep2.Pt(), 0.5*weight);
	  histLeptonEta[i]->Fill(events->lep2.Eta(), 0.5*weight);
	  histMET[i]->Fill(events->MET, 0.5*weight);
	  histLep1MT[i]->Fill(lep2MT, 0.5*weight);

	  if (lep2MT > 120) {
	    MCYield += weight * 0.5;
	    histNJets40[i]->Fill( events->NJets40, 0.5*weight );
	    histNJets80[i]->Fill( events->NJets80, 0.5*weight );
	    histNBtags[i]->Fill( events->NBJetsMedium, 0.5*weight );	    	   
	    histMR[i]->Fill(events->MR, 0.5*weight);
	    histRsq[i]->Fill(events->Rsq, 0.5*weight);	   
	    histMRVsRsq[i]->Fill(events->MR,events->Rsq, 0.5*weight);
	    histMRRsqUnrolled[i]->Fill( MRRsqBin+0.5, 0.5*weight);
	  }	
	}



      } //loop over events
    } //loop over input files
  } //loop over input file groups


  //--------------------------------------------------------------------------------------------------------------
  // Compute Systematic Uncertainties TH2
  //==============================================================================================================
  TH2F * histSystematics = (TH2F*)histMRVsRsq[0]->Clone("TTBar2LSystematics_MRRsq");
  for (int i=1; i<histSystematics->GetXaxis()->GetNbins()+1; ++i) {
    for (int j=1; j<histSystematics->GetYaxis()->GetNbins()+1; ++j) {

      double data = histMRVsRsq[0]->GetBinContent(i,j);
      double MC = 0;
      double MCErrSqr = 0;
      for (uint k = 1 ; k < histMRVsRsq.size(); ++k) {
	MC += histMRVsRsq[k]->GetBinContent(i,j);
	MCErrSqr += pow(histMRVsRsq[k]->GetBinError(i,j),2);
      }
      double dataOverMC = data/MC;
      double dataOverMCErr = 0;
      if (data > 0) {
	dataOverMCErr = (data / MC)*sqrt(1/data + MCErrSqr/pow(MC,2)); 
      } else {
	dataOverMCErr = sqrt( pow(1/MC,2) + MCErrSqr/pow(MC,2));
      }
     
      cout << "Bin " << i << " " << j << " : " << dataOverMC << " +/- " << dataOverMCErr << "\n";
      if ( fabs(dataOverMC - 1) > dataOverMCErr ) {
	cout << "significant\n";
	histSystematics->SetBinContent(i,j,fabs( dataOverMC - 1));
	histSystematics->SetBinError(i,j,dataOverMCErr);
      } else {
	histSystematics->SetBinContent(i,j,0.0);
	histSystematics->SetBinError(i,j,0);
      }
    }
  }


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;
  TLegend *legend = 0;

  //*******************************************************************************************
  //MR
  //*******************************************************************************************
  PlotDataAndStackedBkg( histLeptonPt, processLabels, color, true, "LeptonPt", Label);
  PlotDataAndStackedBkg( histLeptonEta, processLabels, color, true, "LeptonEta", Label);
  PlotDataAndStackedBkg( histMET, processLabels, color, true, "MET", Label);
  PlotDataAndStackedBkg( histLep1MT, processLabels, color, true, "Lep1MT", Label);
  PlotDataAndStackedBkg( histNJets40, processLabels, color, true, "NJets40", Label);
  PlotDataAndStackedBkg( histNJets80, processLabels, color, true, "NJets80", Label);
  PlotDataAndStackedBkg( histNBtags, processLabels, color, true, "NBtags", Label);
  PlotDataAndStackedBkg( histMR, processLabels, color, true, "MR", Label);
  PlotDataAndStackedBkg( histRsq, processLabels, color, true, "Rsq", Label);
  PlotDataAndStackedBkg( histMRRsqUnrolled, processLabels, color, true, "MRRsqUnrolled", Label);


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
  TFile *file = TFile::Open(Form("data/ScaleFactors/RazorMADD2015/TTBarDileptonSystematic_%s_Razor2016.root",label.c_str()), 
          "UPDATE");
  file->cd();
  file->WriteTObject(histSystematics, "TTBarDileptonSystematic", "WriteDelete");

  file->Close();
  delete file;       

  // output file for histograms
  TFile *histFile = TFile::Open("TTJetsDileptonControlRegionHists.root", "RECREATE");
  histFile->cd();
  for (uint i = 0; i < inputfiles.size(); i++) {
      histMRVsRsq[i]->Write();
  }
  histFile->Close();

}







void TTBarDileptonCrossCheck( int option = 0) {
  gROOT->SetBatch();

  vector<string> datafiles;
  vector<vector<string> > bkgfiles;
  vector<string> processLabels;
  vector<int> colors;

  string datafile = "";

  //No Skims  
  datafiles.push_back("Backgrounds/2L/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root");

  vector<string> bkgfiles_ttbar;
  vector<string> bkgfiles_wjets;
  vector<string> bkgfiles_singletop;
  vector<string> bkgfiles_dy;  
  vector<string> bkgfiles_vv; 
  vector<string> bkgfiles_ttv;

  bkgfiles_ttbar.push_back("Backgrounds/2L/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_TTJets_1pb_weighted.root");    
  bkgfiles_wjets.push_back("Backgrounds/2L/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_WJets_1pb_weighted.root"); 
  bkgfiles_singletop.push_back("Backgrounds/2L/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleTop_1pb_weighted.root");  
  bkgfiles_dy.push_back("Backgrounds/2L/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_DYJets_1pb_weighted.root");    
  bkgfiles_vv.push_back("Backgrounds/2L/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_Other_1pb_weighted.root");   

  bkgfiles.push_back(bkgfiles_ttbar);
  bkgfiles.push_back(bkgfiles_dy);
  bkgfiles.push_back(bkgfiles_vv);
  bkgfiles.push_back(bkgfiles_singletop);
  bkgfiles.push_back(bkgfiles_wjets);

  processLabels.push_back("TTJets");  
  processLabels.push_back("DY");
  processLabels.push_back("WW");
  processLabels.push_back("SingleTop");
  processLabels.push_back("WJets");

  colors.push_back(kGreen+2);
  colors.push_back(kBlue+1);
  colors.push_back(kViolet+3);
  colors.push_back(kOrange-3);
  colors.push_back(kRed+1);
 
  double lumi = 26400;

  //*********************************************************************
  //All final states Control Region
  //*********************************************************************
  if (option == 0) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"MR300Rsq0p15",-1,"MR300Rsq0p15_all");

  else if (option == 1) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles,processLabels, colors, lumi,"4JetMR300Rsq0p15",-1,"MultiJet");

  else if (option == 2) RunSelectTTBarDileptonControlSample(datafiles, bkgfiles, processLabels, colors, lumi, "2JetMR300Rsq0p15",-1,"DiJet");

  gApplication->Terminate();

}
