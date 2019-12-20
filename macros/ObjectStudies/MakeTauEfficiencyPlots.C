//================================================================================================
//
// Simple Example
//
//root -l RazorAnalyzer/macros/ObjectStudies/MakeMuonEfficiencyPlots.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/MuonNtuple/MuonNtuple_PromptGenLevel_TTJets_25ns.root",-1,"Muon")'
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
#include <TCanvas.h>                
#include <TLegend.h>                
#include <TGraphAsymmErrors.h>                

#include "RazorAnalyzer/macros/ObjectStudies/EfficiencyUtils.hh"
#include "RazorAnalyzer/include/TauTree.h"

#endif


bool PassSelection( TauTree* tauTree , int wp) {

  bool pass = false;

  if (wp == 1) {
    if (tauTree->fTauPt > 20 && tauTree->fPassLooseSelection) pass = true;  
  }
  
  return pass;
}


void plotTauEfficiency() {

  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileLoose = new TFile("Efficiency_PromptTau_TTJets_25ns_Loose.root","READ");
  TFile *fileFakesLoose = new TFile("Efficiency_FakeTau_TTJets_25ns_Loose.root","READ");


  TGraphAsymmErrors* effPtLoose = (TGraphAsymmErrors*)fileLoose->Get("Efficiency_Pt");
  TGraphAsymmErrors* effEtaLoose = (TGraphAsymmErrors*)fileLoose->Get("Efficiency_Eta");
  TGraphAsymmErrors* effNpvLoose = (TGraphAsymmErrors*)fileLoose->Get("Efficiency_NPV");

  TGraphAsymmErrors* effFakePtLoose = (TGraphAsymmErrors*)fileFakesLoose->Get("Efficiency_Pt");
  TGraphAsymmErrors* effFakeEtaLoose = (TGraphAsymmErrors*)fileFakesLoose->Get("Efficiency_Eta");
  TGraphAsymmErrors* effFakeNpvLoose = (TGraphAsymmErrors*)fileFakesLoose->Get("Efficiency_NPV");




  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtLoose, "Loose", "LP");

  effPtLoose->SetLineWidth(3);
  effPtLoose->SetLineColor(kBlack);
  effPtLoose->GetXaxis()->SetTitle("Tau p_{T} [GeV/c]");
  effPtLoose->GetYaxis()->SetTitle("Selection Efficiency");
  effPtLoose->GetYaxis()->SetTitleOffset(1.2);
  effPtLoose->GetYaxis()->SetRangeUser(0.0,0.6);

  effPtLoose->SetLineWidth(3);
  effPtLoose->SetLineColor(kBlue);

  effPtLoose->Draw("AP");
  
  legend->Draw();  
  cv->SaveAs("TauSelectionEfficiencyVsPt.gif");
  cv->SaveAs("TauSelectionEfficiencyVsPt.pdf");


 
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effEtaLoose, "Loose", "LP");

  effEtaLoose->SetLineWidth(3);
  effEtaLoose->SetLineColor(kBlack);
  effEtaLoose->GetXaxis()->SetTitle("Tau #eta");
  effEtaLoose->GetYaxis()->SetTitle("Selection Efficiency");
  effEtaLoose->GetYaxis()->SetTitleOffset(1.2);
  effEtaLoose->GetYaxis()->SetRangeUser(0.0,0.5);

  effEtaLoose->SetLineWidth(3);
  effEtaLoose->SetLineColor(kBlue);

  effEtaLoose->Draw("AP");
  
  legend->Draw();  
  cv->SaveAs("TauSelectionEfficiencyVsEta.gif");
  cv->SaveAs("TauSelectionEfficiencyVsEta.pdf");




  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.75,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effNpvLoose, "Loose", "LP");

  effNpvLoose->SetLineWidth(3);
  effNpvLoose->SetLineColor(kBlack);
  effNpvLoose->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
  effNpvLoose->GetYaxis()->SetTitle("Selection Efficiency");
  effNpvLoose->GetYaxis()->SetTitleOffset(1.2);
  effNpvLoose->GetXaxis()->SetRangeUser(5,35);
  effNpvLoose->GetYaxis()->SetRangeUser(0.0,0.5);

  effNpvLoose->SetLineWidth(3);
  effNpvLoose->SetLineColor(kBlue);
  effNpvLoose->Draw("AP");
  
  legend->Draw();  
  cv->SaveAs("TauSelectionEfficiencyVsNpv.gif");
  cv->SaveAs("TauSelectionEfficiencyVsNpv.pdf");






  //***************************************************************
  //Fake Taus : Efficiency Vs Pt
  //***************************************************************
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effFakePtLoose, "Loose", "LP");

  effFakePtLoose->SetLineWidth(3);
  effFakePtLoose->SetLineColor(kBlack);
  effFakePtLoose->GetXaxis()->SetTitle("Tau p_{T} [Gev/C]");
  effFakePtLoose->GetYaxis()->SetTitle("Selection Efficiency");
  effFakePtLoose->GetYaxis()->SetTitleOffset(1.35);

  effFakePtLoose->SetLineWidth(3);
  effFakePtLoose->SetLineColor(kBlue);

  effFakePtLoose->Draw("AP");
  
  effFakePtLoose->GetYaxis()->SetRangeUser(0,0.2);

  legend->Draw();  
  //cv->SetLogy();
  cv->SaveAs("TauSelectionFakeRateVsPt.gif");
  cv->SaveAs("TauSelectionFakeRateVsPt.pdf");




  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.7,0.90,0.85);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effFakeEtaLoose, "Loose", "LP");

  effFakeEtaLoose->SetLineWidth(3);
  effFakeEtaLoose->SetLineColor(kBlack);
  effFakeEtaLoose->GetXaxis()->SetTitle("Tau #eta");
  effFakeEtaLoose->GetYaxis()->SetTitle("Selection Efficiency");
  effFakeEtaLoose->GetYaxis()->SetTitleOffset(1.35);

  effFakeEtaLoose->SetLineWidth(3);
  effFakeEtaLoose->SetLineColor(kBlue);

  effFakeEtaLoose->GetYaxis()->SetRangeUser(0,0.10);

  effFakeEtaLoose->Draw("AP");
  
  legend->Draw();  
  cv->SaveAs("TauSelectionFakeRateVsEta.gif");
  cv->SaveAs("TauSelectionFakeRateVsEta.pdf");




 cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.65,0.90,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effFakeNpvLoose, "Loose", "LP");

  effFakeNpvLoose->SetLineWidth(3);
  effFakeNpvLoose->SetLineColor(kBlack);
  effFakeNpvLoose->GetXaxis()->SetTitle("Number of Reconstructed Vertices");
  effFakeNpvLoose->GetYaxis()->SetTitle("Selection Efficiency");
  effFakeNpvLoose->GetYaxis()->SetTitleOffset(1.35);
  effFakeNpvLoose->GetXaxis()->SetRangeUser(5,35);
  effFakeNpvLoose->GetYaxis()->SetRangeUser(0.0,0.10);

  effFakeNpvLoose->SetLineWidth(3);
  effFakeNpvLoose->SetLineColor(kBlue);
  effFakeNpvLoose->Draw("AP");
  
  legend->Draw();  
  cv->SaveAs("TauSelectionFakeRateVsNpv.gif");
  cv->SaveAs("TauSelectionFakeRateVsNpv.pdf");

}



void MakeFastsimToFullSimCorrectionFactors() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileFullsimLoose = new TFile("Efficiency_PromptTau_TTJets_25ns_Loose_Fullsim.root","READ");
  TFile *fileFastsimLoose = new TFile("Efficiency_PromptTau_TTJets_25ns_Loose_Fastsim.root","READ");

  TH2F* histFullsimLoose = (TH2F*)fileFullsimLoose->Get("Efficiency_PtEta");
  TH2F* histFastsimLoose = (TH2F*)fileFastsimLoose->Get("Efficiency_PtEta");

  TH2F* histSFLoose = (TH2F*)histFullsimLoose->Clone("TauLoose_FastsimScaleFactor");
  histSFLoose->GetXaxis()->SetTitle("Tau p_{T} [GeV/c]");
  histSFLoose->GetYaxis()->SetTitle("Tau #eta");

  //Loose WP
  for (int b=1; b<histSFLoose->GetXaxis()->GetNbins()+1 ; ++b) {
    for (int c=1; c<histSFLoose->GetYaxis()->GetNbins()+1 ; ++c) {
      double sf = histFullsimLoose->GetBinContent(b,c) / histFastsimLoose->GetBinContent(b,c);
      double sferr = 0;
      if ( histFullsimLoose->GetBinContent(b,c) > 0 && histFastsimLoose->GetBinContent(b,c) > 0) {
	sferr = sf*sqrt( pow(histFullsimLoose->GetBinError(b,c)/histFullsimLoose->GetBinContent(b,c),2) +
				pow(histFastsimLoose->GetBinError(b,c)/histFastsimLoose->GetBinContent(b,c),2) );
      }
      histSFLoose->SetBinContent(b,c,sf);
      histSFLoose->SetBinError(b,c,sferr);
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open("TauEffFastsimToFullsimCorrectionFactors.root", "UPDATE");
  file->cd();
  file->WriteTObject(histSFLoose, "TauLoose_FastsimScaleFactor", "WriteDelete");  
  file->WriteTObject(histFullsimLoose, "TauEff_Loose_Fullsim", "WriteDelete");
  file->WriteTObject(histFastsimLoose, "TauEff_Loose_Fastsim", "WriteDelete");
  file->Close();
  delete file;      

}





//=== MAIN MACRO ================================================================================================= 

void ProduceTauEfficiencyPlots(const string inputfile, int wp,  int option = -1, string label = "") {


  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt",";Tau p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt",";Tau p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorEta = new TH1F ("histDenominatorEta",";Tau Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histNumeratorEta = new TH1F ("histNumeratorEta",";Tau Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histDenominatorPhi = new TH1F ("histDenominatorPhi",";Tau Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histNumeratorPhi = new TH1F ("histNumeratorPhi",";Tau Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histDenominatorRho = new TH1F ("histDenominatorRho",";Tau Rho; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorRho = new TH1F ("histNumeratorRho",";Tau Rho; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpv = new TH1F ("histDenominatorNpv",";Tau Npv; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpv = new TH1F ("histNumeratorNpv",";Tau Npv; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpu = new TH1F ("histDenominatorNpu",";Tau Npu; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpu = new TH1F ("histNumeratorNpu",";Tau Npu; Number of Events", 50, 0 , 100);

  TH2F *histDenominatorPtEta = new TH2F ("histDenominatorPtEta",";Tau p_{T} [GeV/c] ; Tau #eta; Number of Events", 34, 0 , 170, 60, -3.0, 3.0);
  TH2F *histNumeratorPtEta = new TH2F ("histNumeratorPtEta",";Tau p_{T} [GeV/c] ; Tau #eta; Number of Events", 34, 0 , 170, 60, -3.0, 3.0);

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  TauTree *tauTree = new TauTree;
  tauTree->LoadTree(inputfile.c_str());
  tauTree->InitTree(TauTree::kTauTreeLight);

  cout << "Total Entries: " << tauTree->tree_->GetEntries() << "\n";
  int nentries = tauTree->tree_->GetEntries();
  for(UInt_t ientry=0; ientry < tauTree->tree_->GetEntries(); ientry++) {       	
    tauTree->tree_->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //Cuts
    if (option == 0) {
      if (tauTree->fTauGenPt < 5) continue;
      if (abs(tauTree->fTauGenEta) > 2.4) continue;
    }
   
    if (option >= 10) {
      if (!(tauTree->fTauPt > 20)) continue;
      if (option == 10) {
	if (!(abs(tauTree->fPdgId) == 0)) continue;
      }
      if (option == 11) {
       if (!(abs(tauTree->fPdgId) == 11)) continue;
      }
      if (option == 13) {
	if (!(abs(tauTree->fPdgId) == 13)) continue;
      }
    }

    if (option == 0) {
      //**** PT - ETA ****
      histDenominatorPtEta->Fill(tauTree->fTauGenPt,tauTree->fTauGenEta);
      if(PassSelection(tauTree,wp)) {
	histNumeratorPtEta->Fill(tauTree->fTauGenPt,tauTree->fTauGenEta);
      }


      //**** PT ****
      histDenominatorPt->Fill(tauTree->fTauGenPt);

      //Numerator
      if(PassSelection(tauTree,wp)) {
        histNumeratorPt->Fill(tauTree->fTauGenPt);        
      }


      //**** Eta ****
      if (fabs(tauTree->fTauGenPt) > 30) {
	histDenominatorEta->Fill(tauTree->fTauGenEta);

	//Numerator
	if(PassSelection(tauTree,wp)) {
	  histNumeratorEta->Fill(tauTree->fTauGenEta);        
	}

      }

      //**** Phi ****
      if (fabs(tauTree->fTauGenEta) < 2.4) {
	histDenominatorPhi->Fill(tauTree->fTauGenPhi);

	//Numerator
	if(PassSelection(tauTree,wp)) {
	  histNumeratorPhi->Fill(tauTree->fTauGenPhi);        
	}

      }

      //**** Rho ****
      if (fabs(tauTree->fTauGenEta) < 2.4) {
	histDenominatorRho->Fill(tauTree->fRho);

	//Numerator
	if(PassSelection(tauTree,wp)) {
	  histNumeratorRho->Fill(tauTree->fRho);        
	}

      }
      //**** Npv ****
      if (fabs(tauTree->fTauGenEta) < 2.4) {
	histDenominatorNpv->Fill(tauTree->fNVertices);

	//Numerator
	if(PassSelection(tauTree,wp)) {
	  histNumeratorNpv->Fill(tauTree->fNVertices);        
	}

      }

      // //**** Npu ****
      // if (fabs(tauTree->fTauGenEta) < 2.4) {
      //   histDenominatorNpu->Fill(tauTree->fNPU);

      //   //Numerator
      //   if(PassSelection(tauTree,wp)) {
      //     histNumeratorNpu->Fill(tauTree->fNPU);        
      //   }

      // }
    }

    if (option == 1 || option >= 10) {
    //**** PT - ETA ****
    histDenominatorPtEta->Fill(tauTree->fTauPt,tauTree->fTauEta);
    if(PassSelection(tauTree,wp)) {
      histNumeratorPtEta->Fill(tauTree->fTauPt,tauTree->fTauEta);
    }


    //**** PT ****
      histDenominatorPt->Fill(tauTree->fTauPt);

      //Numerator
      if(PassSelection(tauTree,wp)) {
        histNumeratorPt->Fill(tauTree->fTauPt);        
      }


    //**** Eta ****
    if (fabs(tauTree->fTauPt) > 30) {
      histDenominatorEta->Fill(tauTree->fTauEta);

      //Numerator
      if(PassSelection(tauTree,wp)) {
        histNumeratorEta->Fill(tauTree->fTauEta);        
      }

    }

    //**** Phi ****
    if (fabs(tauTree->fTauEta) < 2.4) {
      histDenominatorPhi->Fill(tauTree->fTauPhi);

      //Numerator
      if(PassSelection(tauTree,wp)) {
        histNumeratorPhi->Fill(tauTree->fTauPhi);        
      }

    }

    //**** Rho ****
    if (fabs(tauTree->fTauEta) < 2.4) {
      histDenominatorRho->Fill(tauTree->fRho);

      //Numerator
      if(PassSelection(tauTree,wp)) {
        histNumeratorRho->Fill(tauTree->fRho);        
      }

    }
    //**** Npv ****
    if (fabs(tauTree->fTauEta) < 2.4) {
      histDenominatorNpv->Fill(tauTree->fNVertices);

      //Numerator
      if(PassSelection(tauTree,wp)) {
        histNumeratorNpv->Fill(tauTree->fNVertices);        
      }

    }

    // //**** Npu ****
    // if (fabs(tauTree->fTauEta) < 2.4) {
    //   histDenominatorNpu->Fill(tauTree->fNPU);

    //   //Numerator
    //   if(PassSelection(tauTree,wp)) {
    //     histNumeratorNpu->Fill(tauTree->fNPU);        
    //   }

    // }
    }


  }


  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency Plots
  //==============================================================================================================

  TGraphAsymmErrors *efficiency_pt = createEfficiencyGraph(histNumeratorPt, histDenominatorPt, "Efficiency_Pt" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_eta = createEfficiencyGraph(histNumeratorEta, histDenominatorEta, "Efficiency_Eta" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_phi = createEfficiencyGraph(histNumeratorPhi, histDenominatorPhi, "Efficiency_Phi" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_rho = createEfficiencyGraph(histNumeratorRho, histDenominatorRho, "Efficiency_Rho" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_npv = createEfficiencyGraph(histNumeratorNpv, histDenominatorNpv, "Efficiency_Npv" , vector<double>() ,  -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_npu = createEfficiencyGraph(histNumeratorNpu, histDenominatorNpu, "Efficiency_Npu" , vector<double>() ,  -99, -99, 0, 1);  
  TH2F *efficiency_pteta = createEfficiencyHist2D(histNumeratorPtEta, histDenominatorPtEta, "Efficiency_PtEta" , vector<double>() ,vector<double>());  


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;

  cv = new TCanvas("cv","cv",800,600);
  efficiency_pt->Draw("AP");
  efficiency_pt->SetTitle("");
  efficiency_pt->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Pt.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_eta->Draw("AP");
  efficiency_eta->SetTitle("");
  efficiency_eta->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Eta.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_phi->Draw("AP");
  efficiency_phi->SetTitle("");
  efficiency_phi->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Phi.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_rho->Draw("AP");
  efficiency_rho->SetTitle("");
  efficiency_rho->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Rho.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npv->Draw("AP");
  efficiency_npv->SetTitle("");
  efficiency_npv->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Npv.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npu->Draw("AP");
  efficiency_npu->SetTitle("");
  efficiency_npu->GetYaxis()->SetRangeUser(0.0,1.0);
  cv->SaveAs(("Efficiency"+Label+"_Npu.gif").c_str());


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("Efficiency"+Label+".root").c_str(), "UPDATE");
  file->cd();
  file->WriteTObject(efficiency_pt, "Efficiency_Pt", "WriteDelete");
  file->WriteTObject(efficiency_eta, "Efficiency_Eta", "WriteDelete");
  file->WriteTObject(efficiency_phi, "Efficiency_Phi", "WriteDelete");
  file->WriteTObject(efficiency_rho, "Efficiency_Rho", "WriteDelete");
  file->WriteTObject(efficiency_npv, "Efficiency_NPV", "WriteDelete");
  file->WriteTObject(efficiency_npu, "Efficiency_NPU", "WriteDelete");
  file->WriteTObject(efficiency_pteta, "Efficiency_PtEta", "WriteDelete");

  file->Close();
  delete file;       

 }

void MakeTauEfficiencyPlots(int option = 0) {

  if (option == 1) {
    ProduceTauEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/TauNtuple/V1p20/TauNtuple_PromptGenLevel_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fullsim.root", 1, 0, "PromptTau_TTJets_25ns_Loose_Fullsim") ;
    // ProduceTauEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/TauNtuple/TauNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Spring15_25ns.root", 1, 10, "FakeTau_TTJets_25ns_Loose") ;
  } 

  if (option == 2) {
    ProduceTauEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/TauNtuple/V1p20/TauNtuple_PromptGenLevel_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fastsim.root", 1, 0, "PromptTau_TTJets_25ns_Loose_Fastsim") ;
  } 

  // plotTauEfficiency();

  if (option == 100) {
    MakeFastsimToFullSimCorrectionFactors();
  }

}
