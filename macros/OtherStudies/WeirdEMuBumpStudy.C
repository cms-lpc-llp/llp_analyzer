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
#include <TH1F.h>                
#include <TH2F.h>                
#include <TCanvas.h>                
#include <TLegend.h> 
#include <THStack.h> 

#include "RazorAnalyzer/include/ControlSampleEvents.h"

#endif


//=== MAIN MACRO ================================================================================================= 


void RunSelectTTBarDileptonControlSample( string datafile, vector<string> bkgfiles, vector<string> bkgLabels, vector<int> bkgColors, double lumi, string option, int channelOption = -1, string label = "") {
  
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

  TFile *pileupWeightFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Run1PileupWeights.root", "READ");
  TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("PUWeight_Run1");
  assert(pileupWeightHist);

  TFile *eleEffSFFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/ScaleFactors/ElectronSelection_Run2012ReReco_53X.root","READ");
  TH2D *eleEffSFHist = (TH2D*)eleEffSFFile->Get("sfLOOSE");
  assert(eleEffSFHist);

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  double MRBins[8] = {300, 400, 500, 750, 1000, 1500, 2000, 4000};
  double RsqBins[10] = {0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5};

  vector<string> inputfiles;
  vector<string> processLabels;
  vector<int> color;

  bool hasData = false;
  if (datafile != "") {
    hasData = true;
    inputfiles.push_back(datafile);
    processLabels.push_back("Data");
    color.push_back(kBlack);
  }
  assert(bkgfiles.size() == bkgLabels.size());
  for (int i=0; i < int(bkgfiles.size()); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processLabels.push_back(bkgLabels[i]);
     color.push_back(bkgColors[i]);
  }

  vector<TH1F*> histMR;
  vector<TH1F*> histRsq;
  vector<TH1F*> histMuonPt;
  vector<TH1F*> histElectronPt;
  vector<TH1F*> histDileptonDeltaPhi;
  vector<TH1F*> histDileptonPt;
  vector<TH1F*> histDileptonEta;
  vector<TH1F*> histMET;
  vector<TH1F*> histDileptonMass;
  vector<TH1F*> histDileptonCharge;
  vector<TH1F*> histNJets40;
  vector<TH1F*> histNBtags;
  vector<TH2F*> histMRVsRsq;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    histMR.push_back(new TH1F(Form("histMR_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", 50, 0, 1500));
    histRsq.push_back(new TH1F(Form("histRsq_%s",processLabels[i].c_str()), "; R^{2} ; Number of Events", 25, 0, 0.5));
    histMuonPt.push_back(new TH1F(Form("histMuonPt_%s",processLabels[i].c_str()), "; MuonPt [GeV/c] ; Number of Events", 20, 0, 400));
    histElectronPt.push_back(new TH1F(Form("histElectronPt_%s",processLabels[i].c_str()), "; ElectronPt [GeV/c] ; Number of Events", 20, 0, 400));
    histDileptonDeltaPhi.push_back(new TH1F(Form("histDileptonDeltaPhi_%s",processLabels[i].c_str()), "; DileptonDeltaPhi [GeV/c] ; Number of Events", 25, 0, 3.1416));
    histDileptonPt.push_back(new TH1F(Form("histDileptonPt_%s",processLabels[i].c_str()), "; DileptonPt [GeV/c] ; Number of Events", 50, 0, 500));
    histDileptonEta.push_back(new TH1F(Form("histDileptonEta_%s",processLabels[i].c_str()), "; DileptonEta [GeV/c] ; Number of Events", 25, -10, 10));
    histMET.push_back(new TH1F(Form("histMET_%s",processLabels[i].c_str()), "; MET [GeV] ; Number of Events", 25, 0, 500));
    histDileptonMass.push_back(new TH1F(Form("histDileptonMass_%s",processLabels[i].c_str()), "; M_{ll} [GeV/c^{2}]; Number of Events", 100, 0, 1000));
    histDileptonCharge.push_back(new TH1F(Form("histDileptonCharge_%s",processLabels[i].c_str()), "; Charge; Number of Events", 5, -2.5, 2.5));
    histNJets40.push_back(new TH1F(Form("histNJets40_%s",processLabels[i].c_str()), "; Number of Jets (p_{T} > 40); Number of Events", 10, -0.5, 9.5));
    histNBtags.push_back(new TH1F(Form("histNBtags_%s",processLabels[i].c_str()), "; Number of B tags; Number of Events", 5, -0.5, 4.5));
    histMRVsRsq.push_back(new TH2F(Form("histMRVsRsq_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; R^{2}; Number of Events", 7, MRBins, 9, RsqBins));
  }
 
  double dataYield = 0;
  double MCYield = 0;

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  for (uint i=0; i < inputfiles.size(); ++i) {
    ControlSampleEvents *events = new ControlSampleEvents;
    events->LoadTree(inputfiles[i].c_str());

    cout << "process: " << processLabels[i] << " | Total Entries: " << events->tree_->GetEntries() << "\n";
    for(UInt_t ientry=0; ientry < events->tree_->GetEntries(); ientry++) {       	
      events->tree_->GetEntry(ientry);
      
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;      
      //if (ientry > 1000000) break;

      double puWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(events->NPU_0));
      double weight = lumi * events->weight * puWeight;

      //******************************
      //Trigger Selection
      //******************************
      bool passTrigger = false;

      //DiMuon Triggers: Mu17Mu8 , Mu17TkMu8
      if (events->HLTDecision[3] ==true || events->HLTDecision[4] ==true) passTrigger = true;
      //DiElectron Triggers:
      if (events->HLTDecision[12] ==true) passTrigger = true;
      //MuEG Triggers: Mu17Ele8 , Mu8Ele17
      if (events->HLTDecision[6] ==true || events->HLTDecision[7] ==true) passTrigger = true;

      if (!passTrigger) continue;

      //******************************
      //Selection Cuts 
      //******************************
      if (!( (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	     &&
	     (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	     )
	  ) continue;

      //lepton selection
      if (! (events->lep1.Pt() > 25 && events->lep2.Pt() > 25
	     && events->lep1PassLoose && events->lep2PassLoose)
	  ) continue;

      //dilepton mass cut
      if ( (events->lep1+events->lep2).M() < 20) continue;

      //MET cut
      if (option == "Met>40" ||option == "topEnhanced" || option == "TwoJet80" || option == "MR400Rsq0p10"  || option == "BumpRegion" || option == "LooserBumpRegion"
	  ) {
	if ( events->MET < 40 ) continue;
      }

      //b-tagging      
      if (option == "topEnhanced" || option == "TwoJet80" || option == "MR400Rsq0p10"  || option == "BumpRegion" || option == "LooserBumpRegion"
	  ) {
	if ( !( (events->jet1PassCSVLoose || events->jet2PassCSVLoose)
		&& events->jet1.Pt() > 30 && events->jet2.Pt() > 30)
	     ) continue;
      }
      
      //Z-mass window cut
      if ( abs(events->lep1Type) == abs(events->lep2Type) 
      	   && 
      	   (events->lep1+events->lep2).M() > 76 && (events->lep1+events->lep2).M() < 106
      	   ) continue;


      // //Razor signal region cuts
      if (option == "TwoJet80") {
	if (!(events->MR > 0 )) continue;
      }
      
      if (option == "MR400Rsq0p10" || option == "BumpRegion") {
	if (!(events->MR > 400 && events->Rsq > 0.1 )) continue;
      }
      
      if (option == "LooserBumpRegion") {
	if (!(events->MR > 400 && events->Rsq > 0.075 )) continue;
      }
      
     //isolate the e-mu bump
      if (option == "BumpRegion" || option == "LooserBumpRegion" ) {
	if (!( (events->lep1+events->lep2).M() > 260 && (events->lep1+events->lep2).M() < 290)) continue;
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
      

      //******************************
      //Apply Scale Factors
      //******************************
      if (!(hasData && i==0)) {
	double triggerEffScaleFactor = 1.0;
	if ( (abs(events->lep1Type) == 13 && abs(events->lep2Type) == 13)) triggerEffScaleFactor = 0.9700;	


	double leptonEffScaleFactor = 1.0;
	if (abs(events->lep1Type) == 11) {
	  leptonEffScaleFactor *= eleEffSFHist->GetBinContent( eleEffSFHist->GetXaxis()->FindFixBin(fabs(events->lep1.Eta())) , 
							       eleEffSFHist->GetYaxis()->FindFixBin(events->lep1.Pt()));
	}
	if (abs(events->lep2Type) == 11) {
	  leptonEffScaleFactor *= eleEffSFHist->GetBinContent( eleEffSFHist->GetXaxis()->FindFixBin(fabs(events->lep2.Eta())) , 
							       eleEffSFHist->GetYaxis()->FindFixBin(events->lep2.Pt()));
	}
	
	weight = weight * leptonEffScaleFactor;
	weight = weight * triggerEffScaleFactor;
	//cout << weight << " " << leptonEffScaleFactor << " " << triggerEffScaleFactor << " " <<  events->weight << " " << puWeight << "\n";

	if ((abs(events->lep1Type) == 11 && abs(events->lep2Type) == 13) ||
	    (abs(events->lep1Type) == 13 && abs(events->lep2Type) == 11)) weight = weight * (733.0 / 902.0);
	if (abs(events->lep1Type) == abs(events->lep2Type) && abs(events->lep1Type) == 13) weight = weight * (302.0 / 409.0);
	if (abs(events->lep1Type) == abs(events->lep2Type) && abs(events->lep1Type) == 11) weight = weight * (263.0 / 318.0);


      }


      //******************************
      //Fill histograms
      //******************************
      if (hasData && i==0) {
	dataYield += 1.0;
	histDileptonMass[i]->Fill((events->lep1+events->lep2).M());      
	histMET[i]->Fill(events->MET);

	if (abs(events->lep1Type) == 13) {
	  histMuonPt[i]->Fill(events->lep1.Pt());
	  histElectronPt[i]->Fill(events->lep2.Pt());
	} else {
	  histElectronPt[i]->Fill(events->lep1.Pt());
	  histMuonPt[i]->Fill(events->lep2.Pt());
	}
	histDileptonDeltaPhi[i]->Fill( acos(cos(events->lep1.Phi() - events->lep2.Phi())) );
	histDileptonPt[i]->Fill( (events->lep1+events->lep2).Pt() );
	histDileptonEta[i]->Fill( (events->lep1+events->lep2).Eta() );
	histDileptonCharge[i]->Fill( events->lep1Type/abs(events->lep1Type) + events->lep2Type/abs(events->lep2Type));
	histNJets40[i]->Fill( events->NJets40 );
	histNBtags[i]->Fill( events->NBJetsLoose );

	if (events->MR > 0) {
	  histMR[i]->Fill(events->MR);
	  histRsq[i]->Fill(events->Rsq);
	}

	histMRVsRsq[i]->Fill(events->MR,events->Rsq);

	// //fill a ttree
	// MR = events->MR;
	// Rsq = events->Rsq;
	// mll = (events->lep1+events->lep2).M();
	// tree->Fill();

      } else {
	MCYield += weight;
	histDileptonMass[i]->Fill((events->lep1+events->lep2).M(), weight );      
	histMET[i]->Fill(events->MET, weight );
	if (abs(events->lep1Type) == 13) {
	  histMuonPt[i]->Fill(events->lep1.Pt(), weight);
	  histElectronPt[i]->Fill(events->lep2.Pt(), weight);
	} else {
	  histElectronPt[i]->Fill(events->lep1.Pt(), weight);
	  histMuonPt[i]->Fill(events->lep2.Pt(), weight);
	}
	histDileptonDeltaPhi[i]->Fill( acos(cos(events->lep1.Phi() - events->lep2.Phi())) , weight);
	histDileptonPt[i]->Fill( (events->lep1+events->lep2).Pt() , weight);
	histDileptonEta[i]->Fill( (events->lep1+events->lep2).Eta() , weight);

	histDileptonCharge[i]->Fill( events->lep1Type/abs(events->lep1Type) + events->lep2Type/abs(events->lep2Type) , weight);
	histNJets40[i]->Fill( events->NJets40 , weight);
	histNBtags[i]->Fill( events->NBJetsLoose , weight);
	
	if (events->MR > 0) {
	  histMR[i]->Fill(events->MR, weight );
	  histRsq[i]->Fill(events->Rsq, weight );
	}

	histMRVsRsq[i]->Fill(events->MR,events->Rsq, weight);

      }
    }
  }

  cout << "here1\n";

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
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMR = new THStack();

  if (hasData) {
    for (int i = histMR.size()-1; i >= 1; --i) {
      histMR[i]->SetFillColor(color[i]);
      histMR[i]->SetFillStyle(1001);
      
      if ( histMR[i]->Integral() > 0) {
  	stackMR->Add(histMR[i]);
      }
    }
  } else {
    for (int i = histMR.size()-1; i >= 0; --i) {
      histMR[i]->SetFillColor(color[i]);
      histMR[i]->SetFillStyle(1001);
      
      if ( histMR[i]->Integral() > 0) {
  	stackMR->Add(histMR[i]);
      }
    }
  }

  for (uint i = 0 ; i < histMR.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histMR[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(histMR[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stackMR->GetHists()->GetEntries() > 0) {
    stackMR->Draw();
    stackMR->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackMR->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stackMR->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    stackMR->SetMaximum( 1.2* fmax( stackMR->GetMaximum(), histMR[0]->GetMaximum()) );

    if (hasData) {
      histMR[0]->SetLineWidth(2);
      histMR[0]->SetLineColor(color[0]);
      histMR[0]->Draw("e1same");
    }
    legend->Draw();
    cv->SetLogy(false);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_MR%s.gif",Label.c_str()));
    cv->SetLogy(true);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_MR%s_Logy.gif",Label.c_str()));
  }


  //*******************************************************************************************
  //Rsq
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackRsq = new THStack();

  if (hasData) {
    for (int i = histRsq.size()-1; i >= 1; --i) {
      histRsq[i]->SetFillColor(color[i]);
      histRsq[i]->SetFillStyle(1001);
      
      if ( histRsq[i]->Integral() > 0) {
  	stackRsq->Add(histRsq[i]);
      }
    }
  } else {
    for (int i = histRsq.size()-1; i >= 0; --i) {
      histRsq[i]->SetFillColor(color[i]);
      histRsq[i]->SetFillStyle(1001);
      
      if ( histRsq[i]->Integral() > 0) {
  	stackRsq->Add(histRsq[i]);
      }
    }
  }

  for (uint i = 0 ; i < histRsq.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histRsq[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(histRsq[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stackRsq->GetHists()->GetEntries() > 0) {
    stackRsq->Draw();
    stackRsq->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackRsq->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackRsq->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackRsq->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stackRsq->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    stackRsq->SetMaximum( 1.2* fmax( stackRsq->GetMaximum(), histRsq[0]->GetMaximum()) );

    if (hasData) {
      histRsq[0]->SetLineWidth(2);
      histRsq[0]->SetLineColor(color[0]);
      histRsq[0]->Draw("e1same");
    }
    legend->Draw();
    cv->SetLogy(false);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_Rsq%s.gif",Label.c_str()));
    cv->SetLogy(true);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_Rsq%s_Logy.gif",Label.c_str()));
  }




  //*******************************************************************************************
  //DileptonMass
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDileptonMass = new THStack();

  if (hasData) {
    for (int i = histDileptonMass.size()-1; i >= 1; --i) {
      histDileptonMass[i]->SetFillColor(color[i]);
      histDileptonMass[i]->SetFillStyle(1001);
      
      if ( histDileptonMass[i]->Integral() > 0) {
  	stackDileptonMass->Add(histDileptonMass[i]);
      }
    }
  } else {
    for (int i = histDileptonMass.size()-1; i >= 0; --i) {
      histDileptonMass[i]->SetFillColor(color[i]);
      histDileptonMass[i]->SetFillStyle(1001);
      
      if ( histDileptonMass[i]->Integral() > 0) {
  	stackDileptonMass->Add(histDileptonMass[i]);
      }
    }
  }

  for (uint i = 0 ; i < histDileptonMass.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histDileptonMass[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(histDileptonMass[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stackDileptonMass->GetHists()->GetEntries() > 0) {
    stackDileptonMass->Draw();
    stackDileptonMass->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDileptonMass->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackDileptonMass->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDileptonMass->GetHists()->At(0)))->GetYaxis()->GetTitle());    
    stackDileptonMass->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    stackDileptonMass->SetMaximum( 1.2* fmax( stackDileptonMass->GetMaximum(), histDileptonMass[0]->GetMaximum()) );

    if (hasData) {
      histDileptonMass[0]->SetLineWidth(2);
      histDileptonMass[0]->SetLineColor(color[0]);
      histDileptonMass[0]->Draw("e1same");
    }
    legend->Draw();
    cv->SetLogy(false);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_DileptonMass%s.gif",Label.c_str()));
    cv->SetLogy(true);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_DileptonMass%s_Logy.gif",Label.c_str()));
  }





  //*******************************************************************************************
  //MuonPt
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMuonPt = new THStack();

  if (hasData) {
    for (int i = histMuonPt.size()-1; i >= 1; --i) {
      histMuonPt[i]->SetFillColor(color[i]);
      histMuonPt[i]->SetFillStyle(1001);
      
      if ( histMuonPt[i]->Integral() > 0) {
  	stackMuonPt->Add(histMuonPt[i]);
      }
    }
  } else {
    for (int i = histMuonPt.size()-1; i >= 0; --i) {
      histMuonPt[i]->SetFillColor(color[i]);
      histMuonPt[i]->SetFillStyle(1001);
      
      if ( histMuonPt[i]->Integral() > 0) {
  	stackMuonPt->Add(histMuonPt[i]);
      }
    }
  }

  for (uint i = 0 ; i < histMuonPt.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histMuonPt[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(histMuonPt[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stackMuonPt->GetHists()->GetEntries() > 0) {
    stackMuonPt->Draw();
    stackMuonPt->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMuonPt->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackMuonPt->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMuonPt->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stackMuonPt->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    stackMuonPt->SetMaximum( 1.2* fmax( stackMuonPt->GetMaximum(), histMuonPt[0]->GetMaximum()) );

    if (hasData) {
      histMuonPt[0]->SetLineWidth(2);
      histMuonPt[0]->SetLineColor(color[0]);
      histMuonPt[0]->Draw("e1same");
    }
    legend->Draw();
    cv->SetLogy(false);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_MuonPt%s.gif",Label.c_str()));
    cv->SetLogy(true);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_MuonPt%s_Logy.gif",Label.c_str()));
  }






  //*******************************************************************************************
  //ElectronPt
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackElectronPt = new THStack();

  if (hasData) {
    for (int i = histElectronPt.size()-1; i >= 1; --i) {
      histElectronPt[i]->SetFillColor(color[i]);
      histElectronPt[i]->SetFillStyle(1001);
      
      if ( histElectronPt[i]->Integral() > 0) {
  	stackElectronPt->Add(histElectronPt[i]);
      }
    }
  } else {
    for (int i = histElectronPt.size()-1; i >= 0; --i) {
      histElectronPt[i]->SetFillColor(color[i]);
      histElectronPt[i]->SetFillStyle(1001);
      
      if ( histElectronPt[i]->Integral() > 0) {
  	stackElectronPt->Add(histElectronPt[i]);
      }
    }
  }

  for (uint i = 0 ; i < histElectronPt.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histElectronPt[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(histElectronPt[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stackElectronPt->GetHists()->GetEntries() > 0) {
    stackElectronPt->Draw();
    stackElectronPt->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackElectronPt->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackElectronPt->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackElectronPt->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stackElectronPt->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    stackElectronPt->SetMaximum( 1.2* fmax( stackElectronPt->GetMaximum(), histElectronPt[0]->GetMaximum()) );

    if (hasData) {
      histElectronPt[0]->SetLineWidth(2);
      histElectronPt[0]->SetLineColor(color[0]);
      histElectronPt[0]->Draw("e1same");
    }
    legend->Draw();
    cv->SetLogy(false);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_ElectronPt%s.gif",Label.c_str()));
    cv->SetLogy(true);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_ElectronPt%s_Logy.gif",Label.c_str()));
  }






  //*******************************************************************************************
  //DileptonDeltaPhi
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDileptonDeltaPhi = new THStack();

  if (hasData) {
    for (int i = histDileptonDeltaPhi.size()-1; i >= 1; --i) {
      histDileptonDeltaPhi[i]->SetFillColor(color[i]);
      histDileptonDeltaPhi[i]->SetFillStyle(1001);
      
      if ( histDileptonDeltaPhi[i]->Integral() > 0) {
  	stackDileptonDeltaPhi->Add(histDileptonDeltaPhi[i]);
      }
    }
  } else {
    for (int i = histDileptonDeltaPhi.size()-1; i >= 0; --i) {
      histDileptonDeltaPhi[i]->SetFillColor(color[i]);
      histDileptonDeltaPhi[i]->SetFillStyle(1001);
      
      if ( histDileptonDeltaPhi[i]->Integral() > 0) {
  	stackDileptonDeltaPhi->Add(histDileptonDeltaPhi[i]);
      }
    }
  }

  for (uint i = 0 ; i < histDileptonDeltaPhi.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histDileptonDeltaPhi[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(histDileptonDeltaPhi[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stackDileptonDeltaPhi->GetHists()->GetEntries() > 0) {
    stackDileptonDeltaPhi->Draw();
    stackDileptonDeltaPhi->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDileptonDeltaPhi->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackDileptonDeltaPhi->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDileptonDeltaPhi->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stackDileptonDeltaPhi->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    stackDileptonDeltaPhi->SetMaximum( 1.2* fmax( stackDileptonDeltaPhi->GetMaximum(), histDileptonDeltaPhi[0]->GetMaximum()) );

    if (hasData) {
      histDileptonDeltaPhi[0]->SetLineWidth(2);
      histDileptonDeltaPhi[0]->SetLineColor(color[0]);
      histDileptonDeltaPhi[0]->Draw("e1same");
    }
    legend->Draw();
    cv->SetLogy(false);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_DileptonDeltaPhi%s.gif",Label.c_str()));
    cv->SetLogy(true);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_DileptonDeltaPhi%s_Logy.gif",Label.c_str()));
  }


  //*******************************************************************************************
  //DileptonPt
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDileptonPt = new THStack();

  if (hasData) {
    for (int i = histDileptonPt.size()-1; i >= 1; --i) {
      histDileptonPt[i]->SetFillColor(color[i]);
      histDileptonPt[i]->SetFillStyle(1001);
      
      if ( histDileptonPt[i]->Integral() > 0) {
  	stackDileptonPt->Add(histDileptonPt[i]);
      }
    }
  } else {
    for (int i = histDileptonPt.size()-1; i >= 0; --i) {
      histDileptonPt[i]->SetFillColor(color[i]);
      histDileptonPt[i]->SetFillStyle(1001);
      
      if ( histDileptonPt[i]->Integral() > 0) {
  	stackDileptonPt->Add(histDileptonPt[i]);
      }
    }
  }

  for (uint i = 0 ; i < histDileptonPt.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histDileptonPt[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(histDileptonPt[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stackDileptonPt->GetHists()->GetEntries() > 0) {
    stackDileptonPt->Draw();
    stackDileptonPt->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDileptonPt->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackDileptonPt->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDileptonPt->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stackDileptonPt->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    stackDileptonPt->SetMaximum( 1.2* fmax( stackDileptonPt->GetMaximum(), histDileptonPt[0]->GetMaximum()) );

    if (hasData) {
      histDileptonPt[0]->SetLineWidth(2);
      histDileptonPt[0]->SetLineColor(color[0]);
      histDileptonPt[0]->Draw("e1same");
    }
    legend->Draw();
    cv->SetLogy(false);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_DileptonPt%s.gif",Label.c_str()));
    cv->SetLogy(true);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_DileptonPt%s_Logy.gif",Label.c_str()));
  }


  //*******************************************************************************************
  //DileptonEta
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDileptonEta = new THStack();

  if (hasData) {
    for (int i = histDileptonEta.size()-1; i >= 1; --i) {
      histDileptonEta[i]->SetFillColor(color[i]);
      histDileptonEta[i]->SetFillStyle(1001);
      
      if ( histDileptonEta[i]->Integral() > 0) {
  	stackDileptonEta->Add(histDileptonEta[i]);
      }
    }
  } else {
    for (int i = histDileptonEta.size()-1; i >= 0; --i) {
      histDileptonEta[i]->SetFillColor(color[i]);
      histDileptonEta[i]->SetFillStyle(1001);
      
      if ( histDileptonEta[i]->Integral() > 0) {
  	stackDileptonEta->Add(histDileptonEta[i]);
      }
    }
  }

  for (uint i = 0 ; i < histDileptonEta.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histDileptonEta[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(histDileptonEta[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stackDileptonEta->GetHists()->GetEntries() > 0) {
    stackDileptonEta->Draw();
    stackDileptonEta->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDileptonEta->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackDileptonEta->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDileptonEta->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stackDileptonEta->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    stackDileptonEta->SetMaximum( 1.2* fmax( stackDileptonEta->GetMaximum(), histDileptonEta[0]->GetMaximum()) );

    if (hasData) {
      histDileptonEta[0]->SetLineWidth(2);
      histDileptonEta[0]->SetLineColor(color[0]);
      histDileptonEta[0]->Draw("e1same");
    }
    legend->Draw();
    cv->SetLogy(false);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_DileptonEta%s.gif",Label.c_str()));
    cv->SetLogy(true);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_DileptonEta%s_Logy.gif",Label.c_str()));
  }



  //*******************************************************************************************
  //MET
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMET = new THStack();

  if (hasData) {
    for (int i = histMET.size()-1; i >= 1; --i) {
      histMET[i]->SetFillColor(color[i]);
      histMET[i]->SetFillStyle(1001);
      
      if ( histMET[i]->Integral() > 0) {
  	stackMET->Add(histMET[i]);
      }
    }
  } else {
    for (int i = histMET.size()-1; i >= 0; --i) {
      histMET[i]->SetFillColor(color[i]);
      histMET[i]->SetFillStyle(1001);
      
      if ( histMET[i]->Integral() > 0) {
  	stackMET->Add(histMET[i]);
      }
    }
  }

  for (uint i = 0 ; i < histMET.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histMET[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(histMET[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stackMET->GetHists()->GetEntries() > 0) {
    stackMET->Draw();
    stackMET->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMET->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackMET->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMET->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stackMET->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    stackMET->SetMaximum( 1.2* fmax( stackMET->GetMaximum(), histMET[0]->GetMaximum()) );

    if (hasData) {
      histMET[0]->SetLineWidth(2);
      histMET[0]->SetLineColor(color[0]);
      histMET[0]->Draw("e1same");
    }
    legend->Draw();
    cv->SetLogy(false);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_MET%s.gif",Label.c_str()));
    cv->SetLogy(true);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_MET%s_Logy.gif",Label.c_str()));
  }



  //*******************************************************************************************
  //DileptonCharge
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackDileptonCharge = new THStack();

  if (hasData) {
    for (int i = histDileptonCharge.size()-1; i >= 1; --i) {
      histDileptonCharge[i]->SetFillColor(color[i]);
      histDileptonCharge[i]->SetFillStyle(1001);
      
      if ( histDileptonCharge[i]->Integral() > 0) {
  	stackDileptonCharge->Add(histDileptonCharge[i]);
      }
    }
  } else {
    for (int i = histDileptonCharge.size()-1; i >= 0; --i) {
      histDileptonCharge[i]->SetFillColor(color[i]);
      histDileptonCharge[i]->SetFillStyle(1001);
      
      if ( histDileptonCharge[i]->Integral() > 0) {
  	stackDileptonCharge->Add(histDileptonCharge[i]);
      }
    }
  }

  for (uint i = 0 ; i < histDileptonCharge.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histDileptonCharge[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(histDileptonCharge[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stackDileptonCharge->GetHists()->GetEntries() > 0) {
    stackDileptonCharge->Draw();
    stackDileptonCharge->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDileptonCharge->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackDileptonCharge->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDileptonCharge->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stackDileptonCharge->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    stackDileptonCharge->SetMaximum( 1.2* fmax( stackDileptonCharge->GetMaximum(), histDileptonCharge[0]->GetMaximum()) );

    if (hasData) {
      histDileptonCharge[0]->SetLineWidth(2);
      histDileptonCharge[0]->SetLineColor(color[0]);
      histDileptonCharge[0]->Draw("e1same");
    }
    legend->Draw();
    cv->SetLogy(false);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_DileptonCharge%s.gif",Label.c_str()));
    cv->SetLogy(true);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_DileptonCharge%s_Logy.gif",Label.c_str()));
  }



  //*******************************************************************************************
  //NJets40
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackNJets40 = new THStack();

  if (hasData) {
    for (int i = histNJets40.size()-1; i >= 1; --i) {
      histNJets40[i]->SetFillColor(color[i]);
      histNJets40[i]->SetFillStyle(1001);
      
      if ( histNJets40[i]->Integral() > 0) {
  	stackNJets40->Add(histNJets40[i]);
      }
    }
  } else {
    for (int i = histNJets40.size()-1; i >= 0; --i) {
      histNJets40[i]->SetFillColor(color[i]);
      histNJets40[i]->SetFillStyle(1001);
      
      if ( histNJets40[i]->Integral() > 0) {
  	stackNJets40->Add(histNJets40[i]);
      }
    }
  }

  for (uint i = 0 ; i < histNJets40.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histNJets40[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(histNJets40[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stackNJets40->GetHists()->GetEntries() > 0) {
    stackNJets40->Draw();
    stackNJets40->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackNJets40->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackNJets40->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackNJets40->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stackNJets40->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    stackNJets40->SetMaximum( 1.2* fmax( stackNJets40->GetMaximum(), histNJets40[0]->GetMaximum()) );

    if (hasData) {
      histNJets40[0]->SetLineWidth(2);
      histNJets40[0]->SetLineColor(color[0]);
      histNJets40[0]->Draw("e1same");
    }
    legend->Draw();
    cv->SetLogy(false);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_NJets40%s.gif",Label.c_str()));
    cv->SetLogy(true);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_NJets40%s_Logy.gif",Label.c_str()));
  }


  //*******************************************************************************************
  //NBtags
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackNBtags = new THStack();

  if (hasData) {
    for (int i = histNBtags.size()-1; i >= 1; --i) {
      histNBtags[i]->SetFillColor(color[i]);
      histNBtags[i]->SetFillStyle(1001);
      
      if ( histNBtags[i]->Integral() > 0) {
  	stackNBtags->Add(histNBtags[i]);
      }
    }
  } else {
    for (int i = histNBtags.size()-1; i >= 0; --i) {
      histNBtags[i]->SetFillColor(color[i]);
      histNBtags[i]->SetFillStyle(1001);
      
      if ( histNBtags[i]->Integral() > 0) {
  	stackNBtags->Add(histNBtags[i]);
      }
    }
  }

  for (uint i = 0 ; i < histNBtags.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(histNBtags[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(histNBtags[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stackNBtags->GetHists()->GetEntries() > 0) {
    stackNBtags->Draw();
    stackNBtags->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackNBtags->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackNBtags->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackNBtags->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stackNBtags->GetYaxis()->SetTitleOffset(1.5);
    stackNBtags->SetMaximum( 1.2* fmax( stackNBtags->GetMaximum(), histNBtags[0]->GetMaximum()) );

    if (hasData) {
      histNBtags[0]->SetLineWidth(2);
      histNBtags[0]->SetLineColor(color[0]);
      histNBtags[0]->Draw("e1same");
    }
    legend->Draw();
    cv->SetLogy(false);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_NBtags%s.gif",Label.c_str()));
    cv->SetLogy(true);
    cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_NBtags%s_Logy.gif",Label.c_str()));
  }


  // //*******************************************************************************************
  // //MR Vs Rsq
  // //*******************************************************************************************
  // cv = new TCanvas("cv","cv", 800,600);
  // histMRVsRsq[0]->SetStats(false);
  // histMRVsRsq[0]->Draw("colz");
  // cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_EventCountsTTJets_MRVsRsq%s.gif",Label.c_str()));

  // cv = new TCanvas("cv","cv", 800,600);
  // statUnc_MRVsRsq->SetStats(false);
  // statUnc_MRVsRsq->Draw("colz");
  // cv->SaveAs(Form("Razor_TTBarDileptonControlRegion_StatUncTTJets_MRVsRsq%s.gif",Label.c_str()));

  //--------------------------------------------------------------------------------------------------------------
  // Tables
  //==============================================================================================================
  cout << "For Luminosity = " << lumi << " pb^-1\n";
  cout << "Yields : MR > 300 && Rsq > 0.1\n";
  //cout << "TTJets: " << 

  cout << "Yield inside Z Mass window 60-120\n";
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
 
  
  // for(int i=0; i<int(histMRVsRsq.size()); i++) {
  //   file->WriteTObject(histMRVsRsq[i], Form("histMRVsRsq_%s",processLabels[i].c_str()), "WriteDelete");
  // }
  // file->WriteTObject(statUnc_MRVsRsq,"statUnc_MRVsRsq_TTJets","WriteDelete");

  // file->Close();
  // delete file;       

}







void WeirdEMuBumpStudy( int option = -1, string label = "") {

  string datafile = "";
  //datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RunOneRazorControlRegions_DileptonSkim_DoubleMuParked_GoodLumi.root";
  //string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkim/RunOneRazorControlRegions_DileptonSkimRazorSkim_DoubleElectron_GoodLumi.root";
  //string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkim/RunOneRazorControlRegions_DileptonSkimRazorSkim_MuEG_GoodLumi.root";

  //string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkim/RunOneRazorControlRegions_DileptonSkimRazorSkim_MuEG_GoodLumi.root";
  //string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkim/RunOneRazorControlRegions_DileptonSkimRazorSkim_HTAndHTMHTParked_GoodLumi.root";
  //string datafile = "";
  vector<string> inputfiles;
  vector<string> processLabels;
  vector<int> colors;

  //No Skims
  //datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RunOneRazorControlRegions_DileptonSkim_MuEG_GoodLumi.root";
  //datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/old/RunOneRazorControlRegions_DileptonSkim_DoubleMuParked_GoodLumi.root";
  // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RunOneRazorControlRegions_DileptonSkim_VV_1pb_weighted.root");
  // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RunOneRazorControlRegions_DileptonSkim_TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_1pb_weighted.root");
  // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RunOneRazorControlRegions_DileptonSkim_SingleTop_1pb_weighted.root");
  // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RunOneRazorControlRegions_DileptonSkim_TTV_1pb_weighted.root");
  // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RunOneRazorControlRegions_DileptonSkim_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_1pb_weighted.root");

  //MR300 Skim
  //datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkimMR300/RunOneRazorControlRegions_DileptonSkimRazorSkim_MuEG_GoodLumi.root";
  //datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkimMR300/RunOneRazorControlRegions_DileptonSkimRazorSkim_DoubleMuParked_GoodLumi.root";
  datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkimMR300/RunOneRazorControlRegions_DileptonSkimRazorSkim_DoubleElectron_GoodLumi.root";
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkimMR300/RunOneRazorControlRegions_DileptonSkimRazorSkim_VV_1pb_weighted.root");
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkimMR300/RunOneRazorControlRegions_DileptonSkimRazorSkim_TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_1pb_weighted.root");
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkimMR300/RunOneRazorControlRegions_DileptonSkimRazorSkim_SingleTop_1pb_weighted.root");
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkimMR300/RunOneRazorControlRegions_DileptonSkimRazorSkim_TTV_1pb_weighted.root");
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RunOneRazorControlRegions/RazorSkimMR300/RunOneRazorControlRegions_DileptonSkimRazorSkim_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_1pb_weighted.root");
  



  processLabels.push_back("VV");
  processLabels.push_back("TTJets");  
  processLabels.push_back("SingleTop");
  processLabels.push_back("TT+V");
  processLabels.push_back("DY");

  colors.push_back(kOrange+1);
  colors.push_back(kRed);
  colors.push_back(kBlack);
  colors.push_back(kGray);
  colors.push_back(kGreen+2);
 

  //*********************************************************************
  //E-Mu Control Region
  //*********************************************************************
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"Inclusive",0,"Inclusive_emu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"Met>40",0,"MetGreaterThan40_emu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"topEnhanced",0,"TopEnhanced_emu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"TwoJet80",0,"TwoJet80_emu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"MR400Rsq0p10",0,"MR400Rsq0p10_emu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"BumpRegion",0,"BumpRegion_emu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"LooserBumpRegion",0,"LooserBumpRegion_emu");


  //*********************************************************************
  //MuMu Control Region
  //*********************************************************************
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19751,"Inclusive",3,"Inclusive_mumu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19751,"TwoJet80",3,"TwoJet80_mumu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19751,"MR400Rsq0p10",3,"MR400Rsq0p10_mumu");

  //*********************************************************************
  //E-E Control Region
  //*********************************************************************
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19789,"Inclusive",2,"Inclusive_ee");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19789,"TwoJet80",2,"TwoJet80_ee");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19789,"MR400Rsq0p10",2,"MR400Rsq0p10_ee");


  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19780,"",0,"emu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels,colors, 19789, "",1,"eemumu");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels,colors, 19789,"", 2,"ee");
  //RunSelectTTBarDileptonControlSample(datafile, inputfiles,processLabels, colors, 19751, "TwoJet80", 3,"mumu");

  


}


//**********************
//E-Mu Yields
//**********************
//Inclusive
// Data: 80236
// MC: 77872.1

//Met>40
//Data: 44398
//MC: 44294.6

//Top Enhanced
//Data: 18486
//MC: 20057.3

//TwoJet80
//Data: 6465
//MC: 7527.79


//**********************
//Mu-Mu Yields
//**********************

//TwoJet80
//Data: 5731
//MC: 6126.44

//**********************
//E-E Yields
//**********************


// Z->MM
// Data: 7.47318e+06
// MC: 7.62973e+06 Powheg
// MC: 7.62596e+06 Madgraph

// Z->EE
//Data: 5.27829e+06
//MC: 5.37693e+06 Madgraph


