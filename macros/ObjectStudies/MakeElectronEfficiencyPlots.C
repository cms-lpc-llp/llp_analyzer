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
#include <TCanvas.h>                
#include <TGraphAsymmErrors.h>                
#include <TLegend.h>                

#include "RazorAnalyzer/macros/ObjectStudies/EfficiencyUtils.hh"
#include "RazorAnalyzer/include/ElectronTree.h"

#endif

double GetEffArea( double scEta ) {

  double effArea = 0.0;
  // //Effective areas below are for the sum of Neutral Hadrons + Photons
  // if (fabs(scEta) < 1.0) {
  //   effArea = 0.0960;
  // } else if (fabs(scEta) < 1.479) {
  //   effArea = 0.0947;	
  // } else if (fabs(scEta) < 2.0) {
  //   effArea = 0.0580;	
  // } else if (fabs(scEta) < 2.2) {
  //   effArea = 0.0688;	
  // } else if (fabs(scEta) < 2.3) {
  //   effArea = 0.0967;	
  // } else if (fabs(scEta) < 2.4) {
  //   effArea = 0.1195;	
  // } else if (fabs(scEta) < 2.5) {
  //   effArea = 0.1475;	
  // }  
  
  //Effective areas using 90% region
  if (fabs(scEta) < 1.0) {
    effArea = 0.1752;
  } else if (fabs(scEta) < 1.479) {
    effArea = 0.1862;
  } else if (fabs(scEta) < 2.0) {
    effArea = 0.1411;	
  } else if (fabs(scEta) < 2.2) {
    effArea = 0.1534;	
  } else if (fabs(scEta) < 2.3) {
    effArea = 0.1903;	
  } else if (fabs(scEta) < 2.4) {
    effArea = 0.2243;	
  } else {
    effArea = 0.2687;	
  } 

  return effArea;

}

Bool_t passPreselection( ElectronTree *eleTree) {

    bool pass = false;
    if(fabs(eleTree->fEleSCEta) < 1.479) {
      if ( fabs(eleTree->fEleDEtaIn) < 0.02
	   && fabs(eleTree->fEleDPhiIn) < 0.15
	   && eleTree->fEleHoverE < 0.15
	   ) {
	pass = true;
      }
    } else {
      if (fabs(eleTree->fEleDEtaIn) < 0.02
	  && fabs(eleTree->fEleDPhiIn) < 0.15
	  && eleTree->fEleHoverE < 0.15
	  ) {
	pass = true;
      }
    } 
    return pass;
}



Bool_t passEGammaPOGTight( ElectronTree *eleTree) {

  bool pass = false;
    if(fabs(eleTree->fEleSCEta) < 1.479) {
      if ( fabs(eleTree->fEleDEtaIn) < 0.00926
	   && fabs(eleTree->fEleDPhiIn) < 0.0336
	   && eleTree->fEleSigmaIEtaIEta < 0.0101
	   && eleTree->fEleHoverE < 0.0597
	   && fabs(eleTree->fEleD0) < 0.0111
	   && fabs(eleTree->fEleDZ) < 0.466
	   && eleTree->fEleOneOverEMinusOneOverP < 0.012
	   //&& eleTree->fElePFIso04 < 0.1649
	   && eleTree->fElePassConversion
	   && eleTree->fEleNMissHits <= 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(eleTree->fEleDEtaIn) < 0.00724
	  && fabs(eleTree->fEleDPhiIn) < 0.0918
	  && eleTree->fEleSigmaIEtaIEta < 0.0279 
	   && eleTree->fEleHoverE < 0.0615 
	  && fabs(eleTree->fEleD0) < 0.0351 
	   && fabs(eleTree->fEleDZ) <  0.417
	   && eleTree->fEleOneOverEMinusOneOverP <  0.00999
	  // && eleTree->fElePFIso04 < 0.2075
	   && eleTree->fElePassConversion
	   && eleTree->fEleNMissHits <= 1
	  ) {
	pass = true;
      }
    } 
  return pass;
}

Bool_t passEGammaPOGLoose( ElectronTree *eleTree) {

    bool pass = false;
    if(fabs(eleTree->fEleSCEta) < 1.479) {
      if ( fabs(eleTree->fEleDEtaIn) < 0.0105
	   && fabs(eleTree->fEleDPhiIn) < 0.115
	   && eleTree->fEleSigmaIEtaIEta < 0.0103
	   && eleTree->fEleHoverE < 0.104
	   //&& fabs(eleTree->fEleD0) < 0.0261
	   && fabs(eleTree->fEleDZ) < 0.41
	   && eleTree->fEleOneOverEMinusOneOverP < 0.102
	   // && eleTree->fElePFIso04 < 0.24
	   && eleTree->fElePassConversion
	   && eleTree->fEleNMissHits <= 2
	   ) {
	pass = true;
      }
    } else {
      if (fabs(eleTree->fEleDEtaIn) < 0.00814
	  && fabs(eleTree->fEleDPhiIn) < 0.182
	  && eleTree->fEleSigmaIEtaIEta < 0.0301 
	  && eleTree->fEleHoverE < 0.0897 
	  //&& fabs(eleTree->fEleD0) < 0.118 
	  && fabs(eleTree->fEleDZ) <  0.822
	  && eleTree->fEleOneOverEMinusOneOverP <  0.126
	  //&& eleTree->fElePFIso04 < 0.3529
	  && eleTree->fElePassConversion
	  && eleTree->fEleNMissHits <= 1
	  ) {
	pass = true;
      }
    }
    return pass;
}


Bool_t passEGammaPOGVetoID( ElectronTree *eleTree) {

  bool pass = false;
  if(fabs(eleTree->fEleSCEta) < 1.479) {
    if ( fabs(eleTree->fEleDEtaIn) < 0.0152
	 && fabs(eleTree->fEleDPhiIn) < 0.216
	 && eleTree->fEleSigmaIEtaIEta < 0.0114
	 && eleTree->fEleHoverE < 0.181
	 //&& fabs(eleTree->fEleD0) < 0.0564
	 && fabs(eleTree->fEleDZ) < 0.472
	 && eleTree->fEleOneOverEMinusOneOverP < 0.207
	 // && eleTree->fElePFIso04 < 0.3313 
	 && eleTree->fElePassConversion
	 && eleTree->fEleNMissHits <= 2
	 ) {
      pass = true;
    }
  } else {
    if (fabs(eleTree->fEleDEtaIn) < 0.0113
	&& fabs(eleTree->fEleDPhiIn) < 0.237
	&& eleTree->fEleSigmaIEtaIEta < 0.0352 
	&& eleTree->fEleHoverE < 0.116 
	//&& fabs(eleTree->fEleD0) < 0.222
	&& fabs(eleTree->fEleDZ) <  0.921
	&& eleTree->fEleOneOverEMinusOneOverP <  0.174
	//&& eleTree->fElePFIso04 < 0.3816
	&& eleTree->fElePassConversion
	&& eleTree->fEleNMissHits <= 3
	) {
      pass = true;
    }
  } 
  return pass;
}

Bool_t passEGammaPOGVetoIso( ElectronTree *eleTree) {

    bool pass = false;
    if(fabs(eleTree->fEleSCEta) < 1.479) {
      if ( eleTree->fElePFIso04 < 0.3313   
	   ) {
	pass = true;
      }
    } else {
      if (  eleTree->fElePFIso04 < 0.3816
	  ) {
	pass = true;
      }
    } 
    return pass;
}

Bool_t passImprovedIso( ElectronTree *eleTree) {

    bool pass = false;
    if(eleTree->fElePt > 20) {
      pass = bool(eleTree->fMiniIsoDBCorr / eleTree->fElePt  < 0.2);
    } else {
      pass = bool(eleTree->fElePFIso04*eleTree->fElePt < 5);
    }
    return pass;
}

Bool_t passIDMVANonTrigVeto( ElectronTree *eleTree) {

  Int_t subdet = 0;  
  if (fabs(eleTree->fEleSCEta) < 0.8) subdet = 0;
  else if (fabs(eleTree->fEleSCEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (eleTree->fElePt > 10.0) ptBin = 1;

  Double_t MVACut = -999;
  if (subdet == 0 && ptBin == 0) MVACut = -0.11;
  if (subdet == 1 && ptBin == 0) MVACut = -0.55;
  if (subdet == 2 && ptBin == 0) MVACut = -0.60;
  if (subdet == 0 && ptBin == 1) MVACut = -0.16;
  if (subdet == 1 && ptBin == 1) MVACut = -0.65;
  if (subdet == 2 && ptBin == 1) MVACut = -0.74;

  bool pass = false;
  if (eleTree->fIDMVANonTrig > MVACut 
      && fabs(eleTree->fEleIP3dSig) < 4
      ) {
    pass = true;
  }   
  return pass;

  // Int_t subdet = 0;  
  // if (fabs(eleTree->fEleSCEta) < 1.479) subdet = 0;
  // else subdet = 1;
  // Int_t ptBin = 0;
  // if (eleTree->fElePt > 10.0) ptBin = 1;

  // Int_t MVABin = -1;
  // if (subdet == 0 && ptBin == 0) MVABin = 0;
  // if (subdet == 1 && ptBin == 0) MVABin = 1;
  // if (subdet == 0 && ptBin == 1) MVABin = 2;
  // if (subdet == 1 && ptBin == 1) MVABin = 3;

  // Double_t MVACut = -999;
  // if (MVABin == 0) MVACut = 0.0;
  // if (MVABin == 1) MVACut = 0.6; 
  // if (MVABin == 2) MVACut = -0.3;
  // if (MVABin == 3) MVACut = 0.5;  

  //   bool pass = false;
  //   if (eleTree->fIDMVANonTrig > MVACut 
  // 	&& ( (fabs(eleTree->fEleSCEta) < 1.479 && fabs(eleTree->fEleD0) < 0.0166)
  // 	     ||
  // 	     (fabs(eleTree->fEleSCEta) >= 1.479 && fabs(eleTree->fEleD0) < 0.098)
  // 	     )
  // 	) {
  //     pass = true;
  //   }   
  //   return pass;

}




bool PassSelection( ElectronTree* eleTree, int wp ) {

  bool pass = false;
  double dr = fmax(0.05,fmin(0.2, 10/eleTree->fElePt));

  //Veto
  if (wp == 1) {
    if ( eleTree->fElePt > 0 
	 &&  eleTree->fPassMVANonTrigVetoSelection
	 // && passIDMVANonTrigVeto(eleTree)
	 // && passImprovedIso(eleTree) 
	 ) {
      pass = true;
    }
  }

  if (wp == 2) {
    
 
    if ( eleTree->fElePt > 0 
	 //&& eleTree->fPassLooseSelection
	  && passEGammaPOGLoose(eleTree) 
	 && eleTree->fMiniIso / eleTree->fElePt < 0.1
	 //&& eleTree->fMiniIsoDBCorr / eleTree->fElePt < 0.1
	 //&& (eleTree->fMiniIsoCharged + fmax(0.0, eleTree->fMiniIso - eleTree->fRho*GetEffArea(eleTree->fEleSCEta)*pow(dr/0.3,2))) / eleTree->fElePt < 0.1
	 //&& (eleTree->fMiniIsoCharged + fmax(0.0, eleTree->fMiniIso - eleTree->fRhoNeutralCentral*GetEffArea(eleTree->fEleSCEta)*pow(dr/0.3,2))) / eleTree->fElePt < 0.1
	 
	 ) {
      pass = true;
    }
  }

  if (wp == 3) {
    if ( eleTree->fElePt > 0 
	 && eleTree->fPassTightSelection
	 // && passEGammaPOGTight(eleTree)
	 // && (eleTree->fMiniIso - eleTree->fRho*GetEffArea(eleTree->fEleSCEta)*pow(dr/0.3,2)) / eleTree->fElePt < 0.1
	 ) {
      pass = true;
    }

  }

  //isolation only
  if (wp == 11) {
    if ( eleTree->fElePt > 0	 
	 && passImprovedIso(eleTree) 
	 ) {
      pass = true;
    }
  }

  //trigger
  if (wp == 100) {
    if ( (eleTree->fEleTriggerBit & ElectronTree::kEleTrigger_Ele27Loose) == ElectronTree::kEleTrigger_Ele27Loose) pass = true;
  }
  if (wp == 101) {
    if ( (eleTree->fEleTriggerBit & ElectronTree::kEleTrigger_Ele27Tight) == ElectronTree::kEleTrigger_Ele27Tight) pass = true;
  }
  if (wp == 102) {
    if ( (eleTree->fEleTriggerBit & ElectronTree::kEleTrigger_Ele32Tight) == ElectronTree::kEleTrigger_Ele32Tight) pass = true;
  }
  if (wp == 110) {
    if ( (eleTree->fEleTriggerBit & ElectronTree::kEleTrigger_Ele27Loose) == ElectronTree::kEleTrigger_Ele27Loose
	 || (eleTree->fEleTriggerBit & ElectronTree::kEleTrigger_Ele27Tight) == ElectronTree::kEleTrigger_Ele27Tight
	 || (eleTree->fEleTriggerBit & ElectronTree::kEleTrigger_Ele32Tight) == ElectronTree::kEleTrigger_Ele32Tight
	  ) pass = true;
  }


  return pass;  
}


void plotElectronEfficiency() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileVeto = new TFile("Efficiency_PromptElectron_TTJets_25ns_Veto.root","READ");
  TFile *fileLoose = new TFile("Efficiency_PromptElectron_TTJets_25ns_Loose.root","READ");
  TFile *fileTight = new TFile("Efficiency_PromptElectron_TTJets_25ns_Tight.root","READ");
  TFile *fileFakesVeto = new TFile("Efficiency_FakeElectron_TTJets_25ns_Veto.root","READ");
  TFile *fileFakesLoose = new TFile("Efficiency_FakeElectron_TTJets_25ns_Loose.root","READ");
  TFile *fileFakesTight = new TFile("Efficiency_FakeElectron_TTJets_25ns_Tight.root","READ");

  TFile *fileVetoIsoGivenID = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/MyNotes/notes/AN-14-276/trunk/data//Efficiency_Electron_NumeratorImprovedVetoIso_DenominatorNonTrigMVAID.root","READ");
  TFile *fileCSA14VetoIsoGivenID = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/MyNotes/notes/AN-14-276/trunk/data/Efficiency_Electron_NumeratorCSA14VetoIso_DenominatorNonTrigMVAID.root","READ");



  TGraphAsymmErrors* effPtVeto = (TGraphAsymmErrors*)fileVeto->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtLoose = (TGraphAsymmErrors*)fileLoose->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtTight = (TGraphAsymmErrors*)fileTight->Get("Efficiency_Pt");
  TGraphAsymmErrors* effEtaVeto = (TGraphAsymmErrors*)fileVeto->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaLoose = (TGraphAsymmErrors*)fileLoose->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaTight = (TGraphAsymmErrors*)fileTight->Get("Efficiency_Eta");
  TGraphAsymmErrors* effNpvVeto = (TGraphAsymmErrors*)fileVeto->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNpvLoose = (TGraphAsymmErrors*)fileLoose->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNpvTight = (TGraphAsymmErrors*)fileTight->Get("Efficiency_NPV");
  TGraphAsymmErrors* effPtVetoIsoGivenID = (TGraphAsymmErrors*)fileVetoIsoGivenID->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtCSA14VetoIsoGivenID = (TGraphAsymmErrors*)fileCSA14VetoIsoGivenID->Get("Efficiency_Pt");

  TGraphAsymmErrors* effFakePtVeto = (TGraphAsymmErrors*)fileFakesVeto->Get("Efficiency_Pt");
  TGraphAsymmErrors* effFakePtLoose = (TGraphAsymmErrors*)fileFakesLoose->Get("Efficiency_Pt");
  TGraphAsymmErrors* effFakePtTight = (TGraphAsymmErrors*)fileFakesTight->Get("Efficiency_Pt");
  TGraphAsymmErrors* effFakeEtaVeto = (TGraphAsymmErrors*)fileFakesVeto->Get("Efficiency_Eta");
  TGraphAsymmErrors* effFakeEtaLoose = (TGraphAsymmErrors*)fileFakesLoose->Get("Efficiency_Eta");
  TGraphAsymmErrors* effFakeEtaTight = (TGraphAsymmErrors*)fileFakesTight->Get("Efficiency_Eta");
  TGraphAsymmErrors* effFakeNpvVeto = (TGraphAsymmErrors*)fileFakesVeto->Get("Efficiency_NPV");
  TGraphAsymmErrors* effFakeNpvLoose = (TGraphAsymmErrors*)fileFakesLoose->Get("Efficiency_NPV");
  TGraphAsymmErrors* effFakeNpvTight = (TGraphAsymmErrors*)fileFakesTight->Get("Efficiency_NPV");




  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtVetoIsoGivenID, "Improved Veto Isolation", "LP");
  legend->AddEntry(effPtCSA14VetoIsoGivenID, "CSA14 Veto Isolation", "LP");


  effPtVetoIsoGivenID->SetLineWidth(3);
  effPtVetoIsoGivenID->SetLineColor(kBlue);
  effPtVetoIsoGivenID->GetXaxis()->SetTitle("Electron p_{T} [GeV/c]");
  effPtVetoIsoGivenID->GetYaxis()->SetTitle("Isolation Efficiency");
  effPtVetoIsoGivenID->GetYaxis()->SetTitleOffset(1.2);
  effPtCSA14VetoIsoGivenID->SetLineWidth(3);
  effPtCSA14VetoIsoGivenID->SetLineColor(kRed);

  effPtVetoIsoGivenID->Draw("AP");
  effPtCSA14VetoIsoGivenID->Draw("Psame");
  effPtVetoIsoGivenID->SetMinimum(0.65);

  legend->Draw();  
  cv->SaveAs("ElectronIsolationEfficiency.gif");
  cv->SaveAs("ElectronIsolationEfficiency.pdf");

 

  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtVeto, "Veto", "LP");
  legend->AddEntry(effPtLoose, "Loose", "LP");
  legend->AddEntry(effPtTight, "Tight", "LP");

  effPtVeto->SetLineWidth(3);
  effPtVeto->SetLineColor(kBlack);
  effPtVeto->GetXaxis()->SetTitle("Electron p_{T} [GeV/c]");
  effPtVeto->GetYaxis()->SetTitle("Selection Efficiency");
  effPtVeto->GetYaxis()->SetTitleOffset(1.2);

  effPtLoose->SetLineWidth(3);
  effPtLoose->SetLineColor(kBlue);
  effPtTight->SetLineWidth(3);
  effPtTight->SetLineColor(kRed);

  effPtVeto->Draw("AP");
  effPtLoose->Draw("Psame");
  effPtTight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("ElectronSelectionEfficiencyVsPt.gif");
  cv->SaveAs("ElectronSelectionEfficiencyVsPt.pdf");




  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effEtaVeto, "Veto", "LP");
  legend->AddEntry(effEtaLoose, "Loose", "LP");
  legend->AddEntry(effEtaTight, "Tight", "LP");

  effEtaVeto->SetLineWidth(3);
  effEtaVeto->SetLineColor(kBlack);
  effEtaVeto->GetXaxis()->SetTitle("Electron #eta");
  effEtaVeto->GetYaxis()->SetTitle("Selection Efficiency");
  effEtaVeto->GetYaxis()->SetTitleOffset(1.2);

  effEtaLoose->SetLineWidth(3);
  effEtaLoose->SetLineColor(kBlue);
  effEtaTight->SetLineWidth(3);
  effEtaTight->SetLineColor(kRed);

  effEtaVeto->Draw("AP");
  effEtaLoose->Draw("Psame");
  effEtaTight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("ElectronSelectionEfficiencyVsEta.gif");
  cv->SaveAs("ElectronSelectionEfficiencyVsEta.pdf");



 cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.75,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effNpvVeto, "Veto", "LP");
  legend->AddEntry(effNpvLoose, "Loose", "LP");
  legend->AddEntry(effNpvTight, "Tight", "LP");

  effNpvVeto->SetLineWidth(3);
  effNpvVeto->SetLineColor(kBlack);
  effNpvVeto->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
  effNpvVeto->GetYaxis()->SetTitle("Selection Efficiency");
  effNpvVeto->GetYaxis()->SetTitleOffset(1.2);
  effNpvVeto->GetXaxis()->SetRangeUser(5,35);
  effNpvVeto->GetYaxis()->SetRangeUser(0.5,1.0);

  effNpvLoose->SetLineWidth(3);
  effNpvLoose->SetLineColor(kBlue);
  effNpvTight->SetLineWidth(3);
  effNpvTight->SetLineColor(kRed);

  effNpvVeto->Draw("AP");
  effNpvLoose->Draw("Psame");
  effNpvTight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("ElectronSelectionEfficiencyVsNpv.gif");
  cv->SaveAs("ElectronSelectionEfficiencyVsNpv.pdf");






  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.70,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effFakePtVeto, "Veto", "LP");
  legend->AddEntry(effFakePtLoose, "Loose", "LP");
  legend->AddEntry(effFakePtTight, "Tight", "LP");

  effFakePtVeto->SetLineWidth(3);
  effFakePtVeto->SetLineColor(kBlack);
  effFakePtVeto->GetXaxis()->SetTitle("Electron p_{T} [GeV/c]");
  effFakePtVeto->GetYaxis()->SetTitle("Selection Efficiency");
  effFakePtVeto->GetYaxis()->SetTitleOffset(1.35);
  effFakePtVeto->GetYaxis()->SetRangeUser(0,0.2);

  effFakePtLoose->SetLineWidth(3);
  effFakePtLoose->SetLineColor(kBlue);
  effFakePtTight->SetLineWidth(3);
  effFakePtTight->SetLineColor(kRed);

  effFakePtVeto->Draw("AP");
  effFakePtLoose->Draw("Psame");
  effFakePtTight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("ElectronSelectionFakeRateVsPt.gif");
  cv->SaveAs("ElectronSelectionFakeRateVsPt.pdf");

 
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.70,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effFakeEtaVeto, "Veto", "LP");
  legend->AddEntry(effFakeEtaLoose, "Loose", "LP");
  legend->AddEntry(effFakeEtaTight, "Tight", "LP");

  effFakeEtaVeto->SetLineWidth(3);
  effFakeEtaVeto->SetLineColor(kBlack);
  effFakeEtaVeto->GetXaxis()->SetTitle("Electron #eta");
  effFakeEtaVeto->GetYaxis()->SetTitle("Selection Efficiency");
  effFakeEtaVeto->GetYaxis()->SetTitleOffset(1.35);
  effFakeEtaVeto->GetYaxis()->SetRangeUser(0,0.1);

  effFakeEtaLoose->SetLineWidth(3);
  effFakeEtaLoose->SetLineColor(kBlue);
  effFakeEtaTight->SetLineWidth(3);
  effFakeEtaTight->SetLineColor(kRed);

  effFakeEtaVeto->Draw("AP");
  effFakeEtaLoose->Draw("Psame");
  effFakeEtaTight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("ElectronSelectionFakeRateVsEta.gif");
  cv->SaveAs("ElectronSelectionFakeRateVsEta.pdf");

 


 cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.75,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effFakeNpvVeto, "Veto", "LP");
  legend->AddEntry(effFakeNpvLoose, "Loose", "LP");
  legend->AddEntry(effFakeNpvTight, "Tight", "LP");

  effFakeNpvVeto->SetLineWidth(3);
  effFakeNpvVeto->SetLineColor(kBlack);
  effFakeNpvVeto->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
  effFakeNpvVeto->GetYaxis()->SetTitle("Selection Efficiency");
  effFakeNpvVeto->GetYaxis()->SetTitleOffset(1.35);
  effFakeNpvVeto->GetXaxis()->SetRangeUser(5,35);
  effFakeNpvVeto->GetYaxis()->SetRangeUser(0,0.07);

  effFakeNpvLoose->SetLineWidth(3);
  effFakeNpvLoose->SetLineColor(kBlue);
  effFakeNpvTight->SetLineWidth(3);
  effFakeNpvTight->SetLineColor(kRed);

  effFakeNpvVeto->Draw("AP");
  effFakeNpvLoose->Draw("Psame");
  effFakeNpvTight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("ElectronSelectionFakeRateVsNpv.gif");
  cv->SaveAs("ElectronSelectionFakeRateVsNpv.pdf");




}



void plotElectronTriggerEfficiency() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileTriggerEle27Loose = new TFile("Efficiency_PromptElectron_TTJets_50ns_TriggerEle27Loose.root","READ");
  TFile *fileTriggerEle27Tight = new TFile("Efficiency_PromptElectron_TTJets_50ns_TriggerEle27Tight.root","READ");
  TFile *fileTriggerEle32Tight = new TFile("Efficiency_PromptElectron_TTJets_50ns_TriggerEle32Tight.root","READ");
  TFile *fileTriggerEleCombined = new TFile("Efficiency_PromptElectron_TTJets_50ns_TriggerEleCombined.root","READ");

  

  TGraphAsymmErrors* effPtTriggerEle27Loose = (TGraphAsymmErrors*)fileTriggerEle27Loose->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtTriggerEle27Tight = (TGraphAsymmErrors*)fileTriggerEle27Tight->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtTriggerEle32Tight = (TGraphAsymmErrors*)fileTriggerEle32Tight->Get("Efficiency_Pt");
  TGraphAsymmErrors* effPtTriggerEleCombined = (TGraphAsymmErrors*)fileTriggerEleCombined->Get("Efficiency_Pt");

  TGraphAsymmErrors* effEtaTriggerEle27Loose = (TGraphAsymmErrors*)fileTriggerEle27Loose->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaTriggerEle27Tight = (TGraphAsymmErrors*)fileTriggerEle27Tight->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaTriggerEle32Tight = (TGraphAsymmErrors*)fileTriggerEle32Tight->Get("Efficiency_Eta");
  TGraphAsymmErrors* effEtaTriggerEleCombined = (TGraphAsymmErrors*)fileTriggerEleCombined->Get("Efficiency_Eta");



  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effPtTriggerEle27Loose, "Ele27 Loose", "LP");
  legend->AddEntry(effPtTriggerEle27Tight, "Ele27 Tight", "LP");
  legend->AddEntry(effPtTriggerEle32Tight, "Ele32 Tight", "LP");
  legend->AddEntry(effPtTriggerEleCombined, "OR of above", "LP");

  effPtTriggerEle27Loose->SetLineWidth(3);
  effPtTriggerEle27Loose->SetLineColor(kBlack);
  effPtTriggerEle27Loose->GetXaxis()->SetTitle("Electron p_{T} [GeV/c]");
  effPtTriggerEle27Loose->GetYaxis()->SetTitle("Trigger Efficiency");
  effPtTriggerEle27Loose->GetYaxis()->SetTitleOffset(1.2);

  effPtTriggerEle27Tight->SetLineWidth(3);
  effPtTriggerEle27Tight->SetLineColor(kBlue);
  effPtTriggerEle32Tight->SetLineWidth(3);
  effPtTriggerEle32Tight->SetLineColor(kRed);
  effPtTriggerEleCombined->SetLineWidth(3);
  effPtTriggerEleCombined->SetLineColor(kGreen+2);

  effPtTriggerEle27Loose->Draw("AP");
  effPtTriggerEle27Tight->Draw("Psame");
  effPtTriggerEle32Tight->Draw("Psame");
  effPtTriggerEleCombined->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("SingleEleTriggerEfficiencyVsPt.gif");
  cv->SaveAs("SingleEleTriggerEfficiencyVsPt.pdf");



  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effEtaTriggerEle27Loose, "Ele27 Loose", "LP");
  legend->AddEntry(effEtaTriggerEle27Tight, "Ele27 Tight", "LP");
  legend->AddEntry(effEtaTriggerEle32Tight, "Ele32 Tight", "LP");
  legend->AddEntry(effEtaTriggerEleCombined, "OR of above", "LP");

  effEtaTriggerEle27Loose->SetLineWidth(3);
  effEtaTriggerEle27Loose->SetLineColor(kBlack);
  effEtaTriggerEle27Loose->GetXaxis()->SetTitle("Electron #eta");
  effEtaTriggerEle27Loose->GetYaxis()->SetTitle("Trigger Efficiency");
  effEtaTriggerEle27Loose->GetYaxis()->SetTitleOffset(1.2);

  effEtaTriggerEle27Tight->SetLineWidth(3);
  effEtaTriggerEle27Tight->SetLineColor(kBlue);
  effEtaTriggerEle32Tight->SetLineWidth(3);
  effEtaTriggerEle32Tight->SetLineColor(kRed);
  effEtaTriggerEleCombined->SetLineWidth(3);
  effEtaTriggerEleCombined->SetLineColor(kGreen+2);

  effEtaTriggerEle27Loose->Draw("AP");
  effEtaTriggerEle27Tight->Draw("Psame");
  effEtaTriggerEle32Tight->Draw("Psame");
  effEtaTriggerEleCombined->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("SingleEleTriggerEfficiencyVsEta.gif");
  cv->SaveAs("SingleEleTriggerEfficiencyVsEta.pdf");





}




void plotElectronMiniIsoEfficiency() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *filePromptUnCorr = new TFile("Efficiency_PromptElectron_TTJets_25_MiniIsolationUnCorr.root","READ");
  TFile *filePromptDBCorr = new TFile("Efficiency_PromptElectron_TTJets_25_MiniIsolationDBCorr.root","READ");
  TFile *filePromptEACorr = new TFile("Efficiency_PromptElectron_TTJets_25_MiniIsolationEACorrRhoNeutralCentral.root","READ");
  TFile *filePromptEACorr90 = new TFile("Efficiency_PromptElectron_TTJets_25_MiniIsolationEACorr90RhoNeutralCentral.root","READ");
  TFile *fileFakeUnCorr = new TFile("Efficiency_FakeElectron_TTJets_25_MiniIsolationUnCorr.root","READ");
  TFile *fileFakeDBCorr = new TFile("Efficiency_FakeElectron_TTJets_25_MiniIsolationDBCorr.root","READ");
  TFile *fileFakeEACorr = new TFile("Efficiency_FakeElectron_TTJets_25_MiniIsolationEACorrRhoNeutralCentral.root","READ");
  TFile *fileFakeEACorr90 = new TFile("Efficiency_FakeElectron_TTJets_25_MiniIsolationEACorr90RhoNeutralCentral.root","READ");

  TGraphAsymmErrors* effNPVPromptUnCorr = (TGraphAsymmErrors*)filePromptUnCorr->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNPVPromptDBCorr = (TGraphAsymmErrors*)filePromptDBCorr->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNPVPromptEACorr = (TGraphAsymmErrors*)filePromptEACorr->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNPVPromptEACorr90 = (TGraphAsymmErrors*)filePromptEACorr90->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNPVFakeUnCorr = (TGraphAsymmErrors*)fileFakeUnCorr->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNPVFakeDBCorr = (TGraphAsymmErrors*)fileFakeDBCorr->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNPVFakeEACorr = (TGraphAsymmErrors*)fileFakeEACorr->Get("Efficiency_NPV");
  TGraphAsymmErrors* effNPVFakeEACorr90 = (TGraphAsymmErrors*)fileFakeEACorr90->Get("Efficiency_NPV");

 
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.60,0.90,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  // legend->AddEntry(effNPVPromptUnCorr, "Uncorrected", "LP");
  // legend->AddEntry(effNPVPromptDBCorr, "#Delta#beta Corr", "LP");
  legend->AddEntry(effNPVPromptEACorr, "Eff Area Corr (Mean)", "LP");
  legend->AddEntry(effNPVPromptEACorr90, "Eff Area Corr (90%)", "LP");

  effNPVPromptDBCorr->SetLineWidth(3);
  effNPVPromptDBCorr->SetLineColor(kBlack);
  effNPVPromptDBCorr->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
  effNPVPromptDBCorr->GetYaxis()->SetTitle("MiniIsolation Cut Efficiency");
  effNPVPromptDBCorr->GetYaxis()->SetTitleOffset(1.4);

  effNPVPromptUnCorr->SetLineWidth(3);
  effNPVPromptUnCorr->SetLineColor(kGreen+2);
  effNPVPromptEACorr->SetLineWidth(3);
  effNPVPromptEACorr->SetLineColor(kBlue);
  effNPVPromptEACorr90->SetLineWidth(3);
  effNPVPromptEACorr90->SetLineColor(kRed);

  effNPVPromptEACorr->Draw("AP");
  effNPVPromptEACorr90->Draw("Psame");
  // effNPVPromptDBCorr->Draw("Psame");
  // effNPVPromptUnCorr->Draw("Psame");

  effNPVPromptEACorr->GetYaxis()->SetRangeUser(0.85,0.95);
  effNPVPromptEACorr->GetXaxis()->SetRangeUser(2,30);
  
  legend->Draw();  
  cv->SaveAs("MiniIsoEfficiencyVsNPV_PromptElectrons.gif");


  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.65,0.90,0.85);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  // legend->AddEntry(effNPVFakeUnCorr, "Uncorrected", "LP");
  // legend->AddEntry(effNPVFakeDBCorr, "#Delta#beta Corr", "LP");
  legend->AddEntry(effNPVFakeEACorr, "Eff Area Corr (Mean)", "LP");
  legend->AddEntry(effNPVFakeEACorr90, "Eff Area Corr (90%)", "LP");

  effNPVFakeDBCorr->SetLineWidth(3);
  effNPVFakeDBCorr->SetLineColor(kBlack);
  effNPVFakeDBCorr->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
  effNPVFakeDBCorr->GetYaxis()->SetTitle("MiniIsolation Cut Efficiency");
  effNPVFakeDBCorr->GetYaxis()->SetTitleOffset(1.4);

  effNPVFakeUnCorr->SetLineWidth(3);
  effNPVFakeUnCorr->SetLineColor(kGreen+2);
  effNPVFakeEACorr->SetLineWidth(3);
  effNPVFakeEACorr->SetLineColor(kBlue);
  effNPVFakeEACorr90->SetLineWidth(3);
  effNPVFakeEACorr90->SetLineColor(kRed);

  effNPVFakeEACorr->Draw("AP");
  effNPVFakeEACorr90->Draw("Psame");
  // effNPVFakeDBCorr->Draw("Psame");
  // effNPVFakeUnCorr->Draw("Psame");

  // effNPVFakeEACorr->GetYaxis()->SetRangeUser(0.0,0.30);
  effNPVFakeEACorr->GetYaxis()->SetRangeUser(0.005,0.025);
  effNPVFakeEACorr->GetXaxis()->SetRangeUser(2,30);
  
  legend->Draw();  
  cv->SaveAs("MiniIsoEfficiencyVsNPV_FakeElectrons.gif");

}


void MakeFastsimToFullSimCorrectionFactors() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileFullsimTight = new TFile("Efficiency_PromptElectron_TTJets_25ns_Tight_Fullsim.root","READ");
  TFile *fileFullsimVeto = new TFile("Efficiency_PromptElectron_TTJets_25ns_Veto_Fullsim.root","READ");
  TFile *fileFastsimTight = new TFile("Efficiency_PromptElectron_TTJets_25ns_Tight_Fastsim.root","READ");
  TFile *fileFastsimVeto = new TFile("Efficiency_PromptElectron_TTJets_25ns_Veto_Fastsim.root","READ");

  TH2F* histFullsimTight = (TH2F*)fileFullsimTight->Get("Efficiency_PtEta");
  TH2F* histFullsimVeto = (TH2F*)fileFullsimVeto->Get("Efficiency_PtEta");
  TH2F* histFastsimTight = (TH2F*)fileFastsimTight->Get("Efficiency_PtEta");
  TH2F* histFastsimVeto = (TH2F*)fileFastsimVeto->Get("Efficiency_PtEta");

  TH2F* histSFTight = (TH2F*)histFullsimTight->Clone("ElectronTight_FastsimScaleFactor");
  TH2F* histSFVeto = (TH2F*)histFullsimVeto->Clone("ElectronVeto_FastsimScaleFactor");
  histSFTight->GetXaxis()->SetTitle("Electron p_{T} [GeV/c]");
  histSFTight->GetYaxis()->SetTitle("Electron #eta");
  histSFVeto->GetXaxis()->SetTitle("Electron p_{T} [GeV/c]");
  histSFVeto->GetYaxis()->SetTitle("Electron #eta");

  //Tight WP
  for (int b=1; b<histSFTight->GetXaxis()->GetNbins()+1 ; ++b) {
    for (int c=1; c<histSFTight->GetYaxis()->GetNbins()+1 ; ++c) {
      double sf = histFullsimTight->GetBinContent(b,c) / histFastsimTight->GetBinContent(b,c);
      double sferr = 0;
      if ( histFullsimTight->GetBinContent(b,c) > 0 && histFastsimTight->GetBinContent(b,c) > 0) {
	sferr = sf*sqrt( pow(histFullsimTight->GetBinError(b,c)/histFullsimTight->GetBinContent(b,c),2) +
				pow(histFastsimTight->GetBinError(b,c)/histFastsimTight->GetBinContent(b,c),2) );
      }
      histSFTight->SetBinContent(b,c,sf);
      histSFTight->SetBinError(b,c,sferr);
    }
  }

  //Veto WP
  for (int b=1; b<histSFVeto->GetXaxis()->GetNbins()+1 ; ++b) {
    for (int c=1; c<histSFVeto->GetYaxis()->GetNbins()+1 ; ++c) {
      double sf = histFullsimVeto->GetBinContent(b,c) / histFastsimVeto->GetBinContent(b,c);
      double sferr = 0;
      if ( histFullsimVeto->GetBinContent(b,c) > 0 && histFastsimVeto->GetBinContent(b,c) > 0) {
	sferr = sf*sqrt( pow(histFullsimVeto->GetBinError(b,c)/histFullsimVeto->GetBinContent(b,c),2) +
			 pow(histFastsimVeto->GetBinError(b,c)/histFastsimVeto->GetBinContent(b,c),2) );
      }
      histSFVeto->SetBinContent(b,c,sf);
      histSFVeto->SetBinError(b,c,sferr);
    }
  }

   //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open("ElectronEffFastsimToFullsimCorrectionFactors.root", "UPDATE");
  file->cd();
  file->WriteTObject(histSFTight, "ElectronTight_FastsimScaleFactor", "WriteDelete");  
  file->WriteTObject(histSFVeto, "ElectronVeto_FastsimScaleFactor", "WriteDelete");
  file->WriteTObject(histFullsimTight, "ElectronEff_Tight_Fullsim", "WriteDelete");
  file->WriteTObject(histFullsimVeto, "ElectronEff_Veto_Fullsim", "WriteDelete");
  file->WriteTObject(histFastsimTight, "ElectronEff_Tight_Fastsim", "WriteDelete");
  file->WriteTObject(histFastsimVeto, "ElectronEff_Veto_Fastsim", "WriteDelete");
  file->Close();
  delete file;      

}



//=== MAIN MACRO ================================================================================================= 

void ProduceElectronEfficiencyPlots(const string inputfile, int wp, int option = -1, string label = "") {
 
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;


  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt",";Electron p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorEta = new TH1F ("histDenominatorEta",";Electron Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histNumeratorEta = new TH1F ("histNumeratorEta",";Electron Eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histDenominatorPhi = new TH1F ("histDenominatorPhi",";Electron Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histNumeratorPhi = new TH1F ("histNumeratorPhi",";Electron Phi; Number of Events", 50, 0 , 3.2);
  TH1F *histDenominatorRho = new TH1F ("histDenominatorRho",";Electron Rho; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorRho = new TH1F ("histNumeratorRho",";Electron Rho; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpv = new TH1F ("histDenominatorNpv",";Electron Npv; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpv = new TH1F ("histNumeratorNpv",";Electron Npv; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpu = new TH1F ("histDenominatorNpu",";Electron Npu; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpu = new TH1F ("histNumeratorNpu",";Electron Npu; Number of Events", 50, 0 , 100);

  TH2F *histDenominatorPtEta = new TH2F ("histDenominatorPtEta",";Electron p_{T} [GeV/c] ; Electron #eta; Number of Events", 34, 0 , 170, 60, -3.0, 3.0);
  TH2F *histNumeratorPtEta = new TH2F ("histNumeratorPtEta",";Electron p_{T} [GeV/c] ; Electron #eta; Number of Events", 34, 0 , 170, 60, -3.0, 3.0);

  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  ElectronTree *EleTree = new ElectronTree;
  EleTree->LoadTree(inputfile.c_str());
  EleTree->InitTree(ElectronTree::kEleTreeLight);

  cout << "Total Entries: " << EleTree->tree_->GetEntries() << "\n";
  int nentries = EleTree->tree_->GetEntries();
  for(UInt_t ientry=0; ientry < EleTree->tree_->GetEntries(); ientry++) {       	
    EleTree->tree_->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //Cuts
    if (EleTree->fEleGenPt < 5) continue;
    if (abs(EleTree->fEleGenEta) > 2.5) continue;

    //if (!(EleTree->fElePt > 5)) continue;
    //if (option == 1 && !passPreselection(EleTree)) continue;

    //for isolation efficiency, require that it passes ID first
    if (wp == 11) {
      if (!passIDMVANonTrigVeto(EleTree)) continue;
    }

    //for trigger require pass Tight
    if (wp >= 100) {
      if (!(EleTree->fPassTightSelection)) continue;
    }


    // if (EleTree->fElePt < 5) continue;
    //if (!passImprovedIso(EleTree)) continue;
    // if (!(fabs(EleTree->fEleIP3dSig) < 4)) continue;
    // if(!passIDMVANonTrigVeto(EleTree)) continue;

    // //Temp: require loose ID
    //if (!passEGammaPOGLoose(EleTree)) continue;


    if (option == 0) {
      //**** PT - ETA ****
      histDenominatorPtEta->Fill(EleTree->fEleGenPt,EleTree->fEleGenEta);
      if(PassSelection(EleTree,wp)) {
	histNumeratorPtEta->Fill(EleTree->fEleGenPt,EleTree->fEleGenEta);
      }


      //**** PT ****
      histDenominatorPt->Fill(EleTree->fEleGenPt);

      //Numerator
      if(PassSelection(EleTree,wp)) {
        histNumeratorPt->Fill(EleTree->fEleGenPt);        
      }


      //**** Eta ****
      if (fabs(EleTree->fEleGenPt) > 30) {
	histDenominatorEta->Fill(EleTree->fEleGenEta);

	//Numerator
	if(PassSelection(EleTree,wp)) {
	  histNumeratorEta->Fill(EleTree->fEleGenEta);        
	}

      }

      //**** Phi ****
      if (fabs(EleTree->fEleGenEta) < 2.4) {
	histDenominatorPhi->Fill(EleTree->fEleGenPhi);

	//Numerator
	if(PassSelection(EleTree,wp)) {
	  histNumeratorPhi->Fill(EleTree->fEleGenPhi);        
	}

      }

      //**** Rho ****
      if (fabs(EleTree->fEleGenEta) < 2.4) {
	histDenominatorRho->Fill(EleTree->fRho);

	//Numerator
	if(PassSelection(EleTree,wp)) {
	  histNumeratorRho->Fill(EleTree->fRho);        
	}

      }
      //**** Npv ****
      if (fabs(EleTree->fEleGenEta) < 2.4) {
	histDenominatorNpv->Fill(EleTree->fNVertices);

	//Numerator
	if(PassSelection(EleTree,wp)) {
	  histNumeratorNpv->Fill(EleTree->fNVertices);        
	}

      }

      // //**** Npu ****
      // if (fabs(EleTree->fEleGenEta) < 2.4) {
      //   histDenominatorNpu->Fill(EleTree->);

      //   //Numerator
      //   if(PassSelection(EleTree,wp)) {
      //     histNumeratorNpu->Fill(EleTree->);        
      //   }

      // }
    }

    if (option == 1) {

      if (!(EleTree->fElePt > 0)) continue;

      //**** PT - ETA ****
      histDenominatorPtEta->Fill(EleTree->fElePt,EleTree->fEleEta);
      if(PassSelection(EleTree,wp)) {
	histNumeratorPtEta->Fill(EleTree->fElePt,EleTree->fEleEta);
      }


      //**** PT ****
      histDenominatorPt->Fill(EleTree->fElePt);

      //Numerator
      if(PassSelection(EleTree,wp)) {
        histNumeratorPt->Fill(EleTree->fElePt);        
      }


      //**** Eta ****
      if (fabs(EleTree->fElePt) > 30) {
	histDenominatorEta->Fill(EleTree->fEleEta);

	//Numerator
	if(PassSelection(EleTree,wp)) {
	  histNumeratorEta->Fill(EleTree->fEleEta);        
	}

      }

      //**** Phi ****
      if (fabs(EleTree->fEleEta) < 2.4) {
	histDenominatorPhi->Fill(EleTree->fElePhi);

	//Numerator
	if(PassSelection(EleTree,wp)) {
	  histNumeratorPhi->Fill(EleTree->fElePhi);        
	}

      }

      //**** Rho ****
      if (fabs(EleTree->fEleEta) < 2.4) {
	histDenominatorRho->Fill(EleTree->fRho);

	//Numerator
	if(PassSelection(EleTree,wp)) {
	  histNumeratorRho->Fill(EleTree->fRho);        
	}

      }
      //**** Npv ****
      if (EleTree->fElePt > 30 ) {
	histDenominatorNpv->Fill(EleTree->fNVertices);

	//Numerator
	if(PassSelection(EleTree,wp)) {
	  histNumeratorNpv->Fill(EleTree->fNVertices);        
	}

      }

      // //**** Npu ****
      // if (fabs(EleTree->fEleEta) < 2.4) {
      //   histDenominatorNpu->Fill(EleTree->);

      //   //Numerator
      //   if(PassSelection(EleTree,wp)) {
      //     histNumeratorNpu->Fill(EleTree->);        
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
  efficiency_pt->GetXaxis()->SetTitle("Electron p_{T} [GeV/c]");
  efficiency_pt->GetYaxis()->SetTitle("Efficiency");
  efficiency_pt->GetYaxis()->SetTitleOffset(1.2);
  efficiency_pt->SetLineWidth(3);  
  cv->SaveAs(("Efficiency"+Label+"_Pt.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_eta->Draw("AP");
  efficiency_eta->SetTitle("");
  efficiency_eta->GetYaxis()->SetRangeUser(0.0,1.0);
  efficiency_eta->GetXaxis()->SetTitle("Electron #eta");
  efficiency_eta->GetYaxis()->SetTitle("Efficiency");
  efficiency_eta->GetYaxis()->SetTitleOffset(1.2);
  efficiency_eta->SetLineWidth(3);  
  cv->SaveAs(("Efficiency"+Label+"_Eta.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_phi->Draw("AP");
  efficiency_phi->SetTitle("");
  efficiency_phi->GetYaxis()->SetRangeUser(0.0,1.0);
  efficiency_phi->GetXaxis()->SetTitle("Electron #phi");
  efficiency_phi->GetYaxis()->SetTitle("Efficiency");
  efficiency_phi->GetYaxis()->SetTitleOffset(1.2);
  cv->SaveAs(("Efficiency"+Label+"_Phi.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_rho->Draw("AP");
  efficiency_rho->SetTitle("");
  efficiency_rho->GetYaxis()->SetRangeUser(0.0,1.0);
  efficiency_rho->GetXaxis()->SetTitle("#rho");
  efficiency_rho->GetYaxis()->SetTitle("Efficiency");
  efficiency_rho->GetYaxis()->SetTitleOffset(1.2);
  cv->SaveAs(("Efficiency"+Label+"_Rho.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npv->Draw("AP");
  efficiency_npv->SetTitle("");
  efficiency_npv->GetYaxis()->SetRangeUser(0.0,1.0);
  efficiency_npv->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
  efficiency_npv->GetYaxis()->SetTitle("Efficiency");
  efficiency_npv->GetYaxis()->SetTitleOffset(1.2);
  efficiency_npv->SetLineWidth(3);  
  efficiency_npv->GetXaxis()->SetRangeUser(0,40);
  cv->SaveAs(("Efficiency"+Label+"_Npv.gif").c_str());

  cv = new TCanvas("cv","cv",800,600);
  efficiency_npu->Draw("AP");
  efficiency_npu->SetTitle("");
  efficiency_npu->GetYaxis()->SetRangeUser(0.0,1.0);
  efficiency_npu->GetXaxis()->SetTitle("Number of Pileup Interactions");
  efficiency_npu->GetYaxis()->SetTitle("Efficiency");
  efficiency_npu->GetYaxis()->SetTitleOffset(1.2);
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

void MakeElectronEfficiencyPlots( int option = 0) {

  if (option == 1) {
   
    ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p20/ElectronNtuple_PromptGenLevel_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fullsim.root", 1, 0, "PromptElectron_TTJets_25ns_Veto_Fullsim");
    ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p20/ElectronNtuple_PromptGenLevel_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fullsim.root", 3, 0, "PromptElectron_TTJets_25ns_Tight_Fullsim");    

    
    // ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p19//ElectronNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 1, 1, "FakeElectron_TTJets_25ns_Veto");
    //  ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p19//ElectronNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "FakeElectron_TTJets_25ns_Loose");
    //  ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p19//ElectronNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 3, 1, "FakeElectron_TTJets_25ns_Tight");


  }

  if (option == 2) {   
    ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p20/ElectronNtuple_PromptGenLevel_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fastsim.root", 1, 0, "PromptElectron_TTJets_25ns_Veto_Fastsim");
    ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p20/ElectronNtuple_PromptGenLevel_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns_Fastsim.root", 3, 0, "PromptElectron_TTJets_25ns_Tight_Fastsim");    
  }

  if (option == 3) {
    ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/ElectronNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_50ns.root", 100, 1, "PromptElectron_TTJets_50ns_TriggerEle27Loose");
    ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/ElectronNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_50ns.root", 101, 1, "PromptElectron_TTJets_50ns_TriggerEle27Tight");
    ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/ElectronNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_50ns.root", 102, 1, "PromptElectron_TTJets_50ns_TriggerEle32Tight");
    ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/ElectronNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_50ns.root", 110, 1, "PromptElectron_TTJets_50ns_TriggerEleCombined");
    plotElectronTriggerEfficiency();
    return;
  }

  //test miniisolation efficiency
  if (option == 3) {
     // ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p18/ElectronNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "PromptElectron_TTJets_25_MiniIsolationUnCorr");
     // ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p18/ElectronNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "FakeElectron_TTJets_25_MiniIsolationUnCorr");   
    // ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p18/ElectronNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "PromptElectron_TTJets_25_MiniIsolationDBCorr");
    // ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p18/ElectronNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "FakeElectron_TTJets_25_MiniIsolationDBCorr");   
    // ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p18/ElectronNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "PromptElectron_TTJets_25_MiniIsolationEACorrRhoNeutralCentral");
    // ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p18/ElectronNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "FakeElectron_TTJets_25_MiniIsolationEACorrRhoNeutralCentral");   
     // ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p18/ElectronNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "PromptElectron_TTJets_25_MiniIsolationEACorr90RhoNeutralCentral");
     // ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p18/ElectronNtuple_Fake_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "FakeElectron_TTJets_25_MiniIsolationEACorr90RhoNeutralCentral");   
 

    //ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/V1p18/ElectronNtuple_Prompt_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.root", 2, 1, "PromptElectron_DYJets_25_MiniIsolationEACorr");
 
    plotElectronMiniIsoEfficiency();
    return;
  }


  if (option == 11) {
    ProduceElectronEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/ElectronNtuple_Prompt_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_50ns.root", 11, 1, "ElectronID_TTJets_VetoIso");
  }

  plotElectronEfficiency();

  if (option == 100) {
    MakeFastsimToFullSimCorrectionFactors();
  }

}
