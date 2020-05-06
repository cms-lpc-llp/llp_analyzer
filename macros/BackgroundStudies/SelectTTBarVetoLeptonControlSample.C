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

const Int_t NComponents = 10;
int color[NComponents] = {kRed, kGreen+2, kBlue, kMagenta, kCyan, kBlack, kOrange, kGray, kBlack, kBlack};
//=== MAIN MACRO ================================================================================================= 


void RunSelectTTBarDileptonVetoLeptonControlSample( vector<string> inputfiles, vector<string> processLabels, int option = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;


  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;
  double lumi = 5000;

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  vector<string> legendLabels;
  legendLabels.push_back("Prompt Lepton");
  legendLabels.push_back("Fake Lepton");


  double MRBins[8] = {300, 400, 500, 750, 1000, 1500, 2000, 4000};
  double RsqBins[10] = {0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5};

  vector<TH1F*> histMR;
  vector<TH1F*> histRsq;
  vector<TH1F*> histDileptonMass;
  vector<TH2F*> histMRVsRsq;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < legendLabels.size(); ++i) {
    histMR.push_back(new TH1F(Form("histMR_%s",legendLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", 100, 0, 3000));
    histRsq.push_back(new TH1F(Form("histRsq_%s",legendLabels[i].c_str()), "; R^{2} ; Number of Events", 100, 0, 2.0));
    histDileptonMass.push_back(new TH1F(Form("histDileptonMass_%s",legendLabels[i].c_str()), "; M_{ll} [GeV/c^{2}]; Number of Events", 100, 0, 1000));
    histMRVsRsq.push_back(new TH2F(Form("histMRVsRsq_%s",legendLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; R^{2}; Number of Events", 7, MRBins, 9, RsqBins));
  }
 
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

      //******************************
      //Selection Cuts 
      //******************************
      if (!( (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	     &&
	     (abs(events->lep1Type) == 11 || abs(events->lep1Type) == 13)
	     )
	  ) continue;

      //lepton 1 selection
      if (! (events->lep1.Pt() > 25 && events->lep1PassTight)
	  ) continue;

      //probe selection
      if (! (events->lep2.Pt() > 25 && events->lep1Type*events->lep2Type < 0
	     && events->lep2PassVetoID)
	  ) continue;

      //b-tagging
      if ( !( (events->bjet1PassMedium && events->bjet2PassMedium)
	      && events->bjet1.Pt() > 30 && events->bjet2.Pt() > 30)
	   ) continue;
      
      //dilepton mass cut
      if ( (events->lep1+events->lep2).M() < 20) continue;

      //Z-mass window cut
      if ( abs(events->lep1Type) == abs(events->lep2Type) 
	   && 
	   (events->lep1+events->lep2).M() > 76 && (events->lep1+events->lep2).M() < 106
	   ) continue;
  
      //Razor signal region cuts
      if (!(events->MR > 300 && events->Rsq > 0.1)) continue;
    
      if (!(events->lep2PassVeto)) continue;

      //option to enhance dilepton ttbar by cutting on NJets
      //if (!(events->NJets40 < 4)) continue;

      //option to enhance dilepton ttbar by requiring lep2 to have higher pT
      //if (!(events->lep2.Pt() > 50)) continue;

      //******************************
      //Options
      //******************************
      //e only
      if (option == 11 &&
	  !(abs(events->lep2Type) == 11)
	  ) continue;
      
     //mu only
      if (option == 13 &&
	  !(abs(events->lep2Type) == 13)
	  ) continue;
      
     //mu only
      if (option == 15 &&
	  !(abs(events->lep2Type) == 15)
	  ) continue;             

      //******************************
      //Fill histograms
      //******************************
      int histIndex = -1;
      if (events->lep2MatchedGenLepIndex>0) histIndex = 0;
      else histIndex = 1;

      histMR[histIndex]->Fill(events->MR,events->weight*lumi);
      histRsq[histIndex]->Fill(events->Rsq,events->weight*lumi);
      histDileptonMass[histIndex]->Fill((events->lep1+events->lep2).M(),events->weight*lumi);      
      histMRVsRsq[histIndex]->Fill(events->MR,events->Rsq, events->weight*lumi);
    }
  }

  cout << "here1\n";

 
 
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
  for (int i = histMR.size()-1; i >= 0; --i) {
   histMR[i]->SetFillColor(color[i]);
    histMR[i]->SetFillStyle(1001);

    if ( histMR[i]->Integral() > 0) {
      stackMR->Add(histMR[i]);
    }
  }
  for (uint i = 0 ; i < histMR.size(); ++i) {
    legend->AddEntry(histMR[i],(legendLabels[i]).c_str(), "F");
  }

  if (stackMR->GetHists()->GetEntries() > 0) {
    stackMR->Draw();
    stackMR->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackMR->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_TTBarVetoLeptonControlRegion_MR%s.gif",Label.c_str()));
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
  for (int i = histRsq.size()-1; i >= 0; --i) {
    histRsq[i]->SetFillColor(color[i]);
    histRsq[i]->SetFillStyle(1001);

    if ( histRsq[i]->Integral() > 0) {
      stackRsq->Add(histRsq[i]);
    }
  }
  for (uint i = 0 ; i < histRsq.size(); ++i) {
    legend->AddEntry(histRsq[i],(legendLabels[i]).c_str(), "F");
  }

  if (stackRsq->GetHists()->GetEntries() > 0) {
    stackRsq->Draw();
    stackRsq->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackRsq->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackRsq->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackRsq->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_TTBarVetoLeptonControlRegion_Rsq%s.gif",Label.c_str()));
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
  for (int i = histDileptonMass.size()-1; i >= 0; --i) {
    histDileptonMass[i]->SetFillColor(color[i]);
    histDileptonMass[i]->SetFillStyle(1001);

    if ( histDileptonMass[i]->Integral() > 0) {
      stackDileptonMass->Add(histDileptonMass[i]);
    }
  }
  for (uint i = 0 ; i < histDileptonMass.size(); ++i) {
    legend->AddEntry(histDileptonMass[i],(legendLabels[i]).c_str(), "F");
  }

  if (stackDileptonMass->GetHists()->GetEntries() > 0) {
    stackDileptonMass->Draw();
    stackDileptonMass->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackDileptonMass->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackDileptonMass->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackDileptonMass->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_TTBarVetoLeptonControlRegion_DileptonMass%s.gif",Label.c_str()));
  }


  //--------------------------------------------------------------------------------------------------------------
  // Tables
  //==============================================================================================================
  cout << "For Luminosity = " << lumi << " pb^-1\n";
  cout << "Yields : MR > 300 && Rsq > 0.1\n";
  //cout << "TTJets: " << 

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("TTBarVetoLeptonControlRegionPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  for(int i=0; i<int(legendLabels.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",legendLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",legendLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDileptonMass[i], Form("histDileptonMass_%s",legendLabels[i].c_str()), "WriteDelete");
  }

 
  file->Close();
  delete file;       

}



void RunSelectTTBarSingleLeptonVetoLeptonControlSample( vector<string> inputfiles, vector<string> processLabels, int option = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;


  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;
  double lumi = 4000;

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  vector<string> legendLabels;
  legendLabels.push_back("Prompt Lepton");
  legendLabels.push_back("Fake Lepton");


  double MRBins[8] = {300, 400, 500, 750, 1000, 1500, 2000, 4000};
  double RsqBins[10] = {0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5};

  vector<TH1F*> histMR;
  vector<TH1F*> histRsq;
  vector<TH2F*> histMRVsRsq;
  vector<TH1F*> histMR_processes;
  vector<TH1F*> histRsq_processes;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < legendLabels.size(); ++i) {
    histMR.push_back(new TH1F(Form("histMR_%s",legendLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", 100, 0, 3000));
    histRsq.push_back(new TH1F(Form("histRsq_%s",legendLabels[i].c_str()), "; R^{2} ; Number of Events", 100, 0, 2.0));
    histMRVsRsq.push_back(new TH2F(Form("histMRVsRsq_%s",legendLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; R^{2}; Number of Events", 7, MRBins, 9, RsqBins));
  }
   for (uint i=0; i < inputfiles.size(); ++i) {
    histMR_processes.push_back(new TH1F(Form("histMR_%s",processLabels[i].c_str()), "; M_{R} [GeV/c^{2}]; Number of Events", 100, 0, 3000));
    histRsq_processes.push_back(new TH1F(Form("histRsq_%s",processLabels[i].c_str()), "; R^{2} ; Number of Events", 100, 0, 2.0));
  }

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

      //skip broken weight events
      if (events->weight == 1) continue;

      //******************************
      //Selection Cuts 
      //******************************
      //lepton selection
      if (! (events->lep1.Pt() > 25 && events->lep1PassVeto)
	  ) continue;

      //No 2nd lepton
      if (events->lep2.Pt() > 0) continue;

      //b-tagging
      if ( !( (events->bjet1PassMedium && events->bjet2PassMedium)
	      && events->bjet1.Pt() > 30 && events->bjet2.Pt() > 30)
	   ) continue;
      

      //Razor signal region cuts
      if (!(events->MR > 300 && events->Rsq > 0.1)) continue;

      //MT Cut
      if (!(events->lep1MT > 30 && events->lep1MT < 100)) continue;
    

      //******************************
      //Options
      //******************************
      //e only
      if (option == 11 &&
	  !(abs(events->lep1Type) == 11)
	  ) continue;
      
     //mu only
      if (option == 13 &&
	  !(abs(events->lep1Type) == 13)
	  ) continue;
      
     //mu only
      if (option == 15 &&
	  !(abs(events->lep1Type) == 15)
	  ) continue;             

      //******************************
      //Fill histograms
      //******************************
      int histIndex = -1;
      if (events->lep1MatchedGenLepIndex>0) histIndex = 0;
      else histIndex = 1;

      histMR[histIndex]->Fill(events->MR,events->weight*lumi);
      histRsq[histIndex]->Fill(events->Rsq,events->weight*lumi);
      histMRVsRsq[histIndex]->Fill(events->MR,events->Rsq, events->weight*lumi);

      histMR_processes[i]->Fill(events->MR,events->weight*lumi);
      histRsq_processes[i]->Fill(events->Rsq,events->weight*lumi);

    }
  }

  cout << "here1\n";

 
 
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
  for (int i = histMR.size()-1; i >= 0; --i) {
   histMR[i]->SetFillColor(color[i]);
    histMR[i]->SetFillStyle(1001);

    if ( histMR[i]->Integral() > 0) {
      stackMR->Add(histMR[i]);
    }
  }
  for (uint i = 0 ; i < histMR.size(); ++i) {
    legend->AddEntry(histMR[i],(legendLabels[i]).c_str(), "F");
  }

  if (stackMR->GetHists()->GetEntries() > 0) {
    stackMR->Draw();
    stackMR->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackMR->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_TTBarSingleLeponVetoLeptonControlRegion_MR%s.gif",Label.c_str()));
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
  for (int i = histRsq.size()-1; i >= 0; --i) {
    histRsq[i]->SetFillColor(color[i]);
    histRsq[i]->SetFillStyle(1001);

    if ( histRsq[i]->Integral() > 0) {
      stackRsq->Add(histRsq[i]);
    }
  }
  for (uint i = 0 ; i < histRsq.size(); ++i) {
    legend->AddEntry(histRsq[i],(legendLabels[i]).c_str(), "F");
  }

  if (stackRsq->GetHists()->GetEntries() > 0) {
    stackRsq->Draw();
    stackRsq->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackRsq->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackRsq->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackRsq->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_TTBarSingleLeponVetoLeptonControlRegion_Rsq%s.gif",Label.c_str()));
  }



  //*******************************************************************************************
  //MR by processes
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMR_processes = new THStack();
  for (int i = histMR_processes.size()-1; i >= 0; --i) {
   histMR_processes[i]->SetFillColor(color[i]);
    histMR_processes[i]->SetFillStyle(1001);

    if ( histMR_processes[i]->Integral() > 0) {
      stackMR_processes->Add(histMR_processes[i]);
    }
  }
  for (uint i = 0 ; i < histMR_processes.size(); ++i) {
    legend->AddEntry(histMR_processes[i],(processLabels[i]).c_str(), "F");
  }

  if (stackMR_processes->GetHists()->GetEntries() > 0) {
    stackMR_processes->Draw();
    stackMR_processes->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR_processes->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackMR_processes->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR_processes->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_TTBarSingleLeponVetoLeptonControlRegion_MR_processes%s.gif",Label.c_str()));
  }

  //*******************************************************************************************
  //Rsq by processes
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackRsq_processes = new THStack();
  for (int i = histRsq_processes.size()-1; i >= 0; --i) {
    histRsq_processes[i]->SetFillColor(color[i]);
    histRsq_processes[i]->SetFillStyle(1001);

    if ( histRsq_processes[i]->Integral() > 0) {
      stackRsq_processes->Add(histRsq_processes[i]);
    }
  }
  for (uint i = 0 ; i < histRsq_processes.size(); ++i) {
    legend->AddEntry(histRsq_processes[i],(processLabels[i]).c_str(), "F");
  }

  if (stackRsq_processes->GetHists()->GetEntries() > 0) {
    stackRsq_processes->Draw();
    stackRsq_processes->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackRsq_processes->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stackRsq_processes->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackRsq_processes->GetHists()->At(0)))->GetYaxis()->GetTitle());
    legend->Draw();
    cv->SaveAs(Form("Razor_TTBarSingleLeponVetoLeptonControlRegion_Rsq_processes%s.gif",Label.c_str()));
  }



  //--------------------------------------------------------------------------------------------------------------
  // Tables
  //==============================================================================================================
  cout << "For Luminosity = " << lumi << " pb^-1\n";
  cout << "Yields : MR > 300 && Rsq > 0.1\n";
  //cout << "TTJets: " << 

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("TTBarSingleLeponVetoLeptonControlRegionPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  for(int i=0; i<int(legendLabels.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",legendLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",legendLabels[i].c_str()), "WriteDelete");
  }
  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR_processes[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq_processes[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
  }

 
  file->Close();
  delete file;       

}





void SelectTTBarVetoLeptonControlSample( int option = -1, string label = "") {

  vector<string> inputfiles;
  vector<string> processLabels;

  // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/TTBarTagAndProbeRegion/TTBarTagAndProbeRegion_TightPlusVetoIDLeptonSkim_TTJets_25ns_weighted.root");  
  // //inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/TTBarTagAndProbeRegion/TTBarTagAndProbeRegion_TightPlusVetoIDLeptonSkim_DYJetsToLL_HT100ToInf_25ns_weighted.root");
  // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/TTBarTagAndProbeRegion/TTBarTagAndProbeRegion_TightPlusVetoIDLeptonSkim_WJetsToLNu_HT100ToInf_25ns_weighted.root");
  // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/TTBarTagAndProbeRegion/TTBarTagAndProbeRegion_TightPlusVetoIDLeptonSkim_SingleTop_25ns_weighted.root"); 
  // inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/TTBarTagAndProbeRegion/TTBarTagAndProbeRegion_Multiboson_25ns_weighted.root");


  // processLabels.push_back("TTJets");
  // //processLabels.push_back("DYJetsToLL");
  // processLabels.push_back("WJetsToLNu");
  // processLabels.push_back("SingleTop");
  // processLabels.push_back("Multiboson");



  // //*********************************************************************
  // //Dilepton Control Region
  // //*********************************************************************
  // RunSelectTTBarDileptonVetoLeptonControlSample(inputfiles,processLabels,11,"e");
  // RunSelectTTBarDileptonVetoLeptonControlSample(inputfiles,processLabels,13,"mu");
  // // RunSelectTTBarDileptonVetoLeptonControlSample(inputfiles,processLabels,15,"tau");
 





  inputfiles.clear();
  processLabels.clear();

  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonEfficiencySingleLeptonControlRegion/skim/VetoLeptonEfficiencySingleLeptonControlRegion_RazorSkim_TTJets_25ns_weighted.root");  
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonEfficiencySingleLeptonControlRegion/skim/VetoLeptonEfficiencySingleLeptonControlRegion_RazorSkim_DYJetsToLL_HT100ToInf_25ns_weighted.root");
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonEfficiencySingleLeptonControlRegion/skim/VetoLeptonEfficiencySingleLeptonControlRegion_RazorSkim_WJetsToLNu_HT100ToInf_25ns_weighted.root");
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonEfficiencySingleLeptonControlRegion/skim/VetoLeptonEfficiencySingleLeptonControlRegion_RazorSkim_SingleTop_25ns_weighted.root"); 
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonEfficiencySingleLeptonControlRegion/skim/VetoLeptonEfficiencySingleLeptonControlRegion_RazorSkim_QCDHT100ToInf_25ns_weighted.root"); 
  inputfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonEfficiencySingleLeptonControlRegion/skim/VetoLeptonEfficiencySingleLeptonControlRegion_Multiboson_25ns_weighted.root");


  processLabels.push_back("TTJets");
  processLabels.push_back("DYJetsToLL");
  processLabels.push_back("WJetsToLNu");
  processLabels.push_back("SingleTop");
  processLabels.push_back("QCD");
  processLabels.push_back("Multiboson");



  //*********************************************************************
  //Single Lepton Control Region
  //*********************************************************************
  RunSelectTTBarSingleLeptonVetoLeptonControlSample(inputfiles,processLabels,11,"e");
  RunSelectTTBarSingleLeptonVetoLeptonControlSample(inputfiles,processLabels,13,"mu");
  //RunSelectTTBarSingleLeptonVetoLeptonControlSample(inputfiles,processLabels,15,"tau");
 


}
