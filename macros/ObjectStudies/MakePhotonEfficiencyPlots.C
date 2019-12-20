//================================================================================================
//
// Simple Example
//
//root -l RazorAnalyzer/macros/ObjectStudies/MakePhotonEfficiencyPlots.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/PhotonNtuple/PhotonNtuple_PromptGenLevel_TTJets_25ns.root",-1,"Photon")'
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
#include <TStyle.h>                
#include <TGraphAsymmErrors.h>                

#include "RazorAnalyzer/macros/ObjectStudies/EfficiencyUtils.hh"
#include "RazorAnalyzer/include/PhotonTree.h"

#endif


bool PassSelection( PhotonTree* PhoTree , int wp = 0 ) {

  bool pass = false;

  //**********************************
  //Tight Selection
  //**********************************
  if (wp == 2 && PhoTree->fPhoIsTight ) {
    pass = true;
  }

  //**********************************
  //Medium Selection
  //**********************************
  if (wp == 1 && PhoTree->fPhoIsMedium ) {
    pass = true;
  }

  // //**********************************
  // //Loose Selection
  // //**********************************
  // if (wp == 0 && PhoTree->fPhoPassLooseID 
  //     && PhoTree->fPhoPassEleVeto 
  //     && PhoTree->fPhoPassLooseIso 
  //     ) {
  //   pass = true;
  // }

  //**********************************
  //Loose Selection with Pixel Seed Veto
  //**********************************
  if (wp == 0 && PhoTree->fPhoPassLooseID 
      && !PhoTree->fPhoHasPixelSeed
      && PhoTree->fPhoPassLooseIso 
      ) {
    pass = true;
  }
 

  return pass;

}


TGraphAsymmErrors* getEffGraph( string filename, string graphname) {
  
  TFile *f = new TFile(filename.c_str(),"READ");
  TGraphAsymmErrors* graph = (TGraphAsymmErrors*)f->Get(graphname.c_str());
  f->Close();
  delete f;
  if (!graph) {
    cout << "Graph " << graphname << " from " << filename << " could not be retrieved\n";
    assert(false);
  }
  return graph;
}


void plotPhotonEfficiency() {

  TGraphAsymmErrors* effVsPt_GJetFlat_50ns_Loose = getEffGraph("Efficiency_GJetFlat_50ns_Loose_NotCloseToParton.root","Efficiency_Pt");
  TGraphAsymmErrors* effVsPt_GJetFlat_50ns_Medium = getEffGraph("Efficiency_GJetFlat_50ns_Medium_NotCloseToParton.root","Efficiency_Pt");
  TGraphAsymmErrors* effVsPt_GJetFlat_50ns_Tight = getEffGraph("Efficiency_GJetFlat_50ns_Tight_NotCloseToParton.root","Efficiency_Pt");

  TGraphAsymmErrors* effVsPt_ttH_50ns_Loose = getEffGraph("Efficiency_ttH_50ns_Loose_NotCloseToParton.root","Efficiency_Pt");
  TGraphAsymmErrors* effVsPt_ttH_50ns_Medium = getEffGraph("Efficiency_ttH_50ns_Medium_NotCloseToParton.root","Efficiency_Pt");
  TGraphAsymmErrors* effVsPt_ttH_50ns_Tight = getEffGraph("Efficiency_ttH_50ns_Tight_NotCloseToParton.root","Efficiency_Pt");

  TGraphAsymmErrors* effVsEta_GJetFlat_50ns_Loose = getEffGraph("Efficiency_GJetFlat_50ns_Loose_NotCloseToParton.root","Efficiency_Eta");
  TGraphAsymmErrors* effVsEta_GJetFlat_50ns_Medium = getEffGraph("Efficiency_GJetFlat_50ns_Medium_NotCloseToParton.root","Efficiency_Eta");
  TGraphAsymmErrors* effVsEta_GJetFlat_50ns_Tight = getEffGraph("Efficiency_GJetFlat_50ns_Tight_NotCloseToParton.root","Efficiency_Eta");

  TGraphAsymmErrors* effVsEta_ttH_50ns_Loose = getEffGraph("Efficiency_ttH_50ns_Loose_NotCloseToParton.root","Efficiency_Eta");
  TGraphAsymmErrors* effVsEta_ttH_50ns_Medium = getEffGraph("Efficiency_ttH_50ns_Medium_NotCloseToParton.root","Efficiency_Eta");
  TGraphAsymmErrors* effVsEta_ttH_50ns_Tight = getEffGraph("Efficiency_ttH_50ns_Tight_NotCloseToParton.root","Efficiency_Eta");

  TGraphAsymmErrors* effVsNPV_GJetFlat_50ns_Loose = getEffGraph("Efficiency_GJetFlat_50ns_Loose_NotCloseToParton.root","Efficiency_NPV");
  TGraphAsymmErrors* effVsNPV_GJetFlat_50ns_Medium = getEffGraph("Efficiency_GJetFlat_50ns_Medium_NotCloseToParton.root","Efficiency_NPV");
  TGraphAsymmErrors* effVsNPV_GJetFlat_50ns_Tight = getEffGraph("Efficiency_GJetFlat_50ns_Tight_NotCloseToParton.root","Efficiency_NPV");

  TGraphAsymmErrors* effVsNPV_ttH_50ns_Loose = getEffGraph("Efficiency_ttH_50ns_Loose_NotCloseToParton.root","Efficiency_NPV");
  TGraphAsymmErrors* effVsNPV_ttH_50ns_Medium = getEffGraph("Efficiency_ttH_50ns_Medium_NotCloseToParton.root","Efficiency_NPV");
  TGraphAsymmErrors* effVsNPV_ttH_50ns_Tight = getEffGraph("Efficiency_ttH_50ns_Tight_NotCloseToParton.root","Efficiency_NPV");


  TGraphAsymmErrors* effVsPt_GJetFlat_50ns_NotCloseToParton_Loose = getEffGraph("Efficiency_GJetFlat_50ns_Loose_NotCloseToParton.root","Efficiency_Pt");
  TGraphAsymmErrors* effVsPt_ttH_50ns_NotCloseToParton_Loose = getEffGraph("Efficiency_ttH_50ns_Loose_NotCloseToParton.root","Efficiency_Pt");




   //****************************************************************************
  //Make Plots
  //****************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;

  //****************************************************************************
  //POG WP Efficiencies
  //****************************************************************************
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effVsPt_GJetFlat_50ns_Loose, "Cut-based Loose", "LP");
  legend->AddEntry(effVsPt_GJetFlat_50ns_Medium, "Cut-based Medium", "LP");
  legend->AddEntry(effVsPt_GJetFlat_50ns_Tight, "Cut-based Tight", "LP");

  effVsPt_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_Loose->SetLineColor(kBlack);
  effVsPt_GJetFlat_50ns_Loose->GetXaxis()->SetTitle("Photon p_{T} [GeV/c]");
  effVsPt_GJetFlat_50ns_Loose->GetYaxis()->SetTitle("Selection Efficiency");
  effVsPt_GJetFlat_50ns_Loose->GetYaxis()->SetTitleOffset(1.2);

  effVsPt_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_Loose->SetLineColor(kRed);
  effVsPt_GJetFlat_50ns_Medium->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_Medium->SetLineColor(kBlue);
  effVsPt_GJetFlat_50ns_Tight->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_Tight->SetLineColor(kGreen+2);

  effVsPt_GJetFlat_50ns_Loose->Draw("AP");
  effVsPt_GJetFlat_50ns_Medium->Draw("Psame");
  effVsPt_GJetFlat_50ns_Tight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("PhotonEfficiencyVsPt.gif");
  cv->SaveAs("PhotonEfficiencyVsPt.pdf");


  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effVsEta_GJetFlat_50ns_Loose, "Cut-based Loose", "LP");
  legend->AddEntry(effVsEta_GJetFlat_50ns_Medium, "Cut-based Medium", "LP");
  legend->AddEntry(effVsEta_GJetFlat_50ns_Tight, "Cut-based Tight", "LP");

  effVsEta_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsEta_GJetFlat_50ns_Loose->SetLineColor(kBlack);
  effVsEta_GJetFlat_50ns_Loose->GetXaxis()->SetTitle("Photon #eta");
  effVsEta_GJetFlat_50ns_Loose->GetYaxis()->SetTitle("Selection Efficiency");
  effVsEta_GJetFlat_50ns_Loose->GetYaxis()->SetTitleOffset(1.2);

  effVsEta_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsEta_GJetFlat_50ns_Loose->SetLineColor(kRed);
  effVsEta_GJetFlat_50ns_Medium->SetLineWidth(3);
  effVsEta_GJetFlat_50ns_Medium->SetLineColor(kBlue);
  effVsEta_GJetFlat_50ns_Tight->SetLineWidth(3);
  effVsEta_GJetFlat_50ns_Tight->SetLineColor(kGreen+2);

  effVsEta_GJetFlat_50ns_Loose->Draw("AP");
  effVsEta_GJetFlat_50ns_Medium->Draw("Psame");
  effVsEta_GJetFlat_50ns_Tight->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("PhotonEfficiencyVsEta.gif");
  cv->SaveAs("PhotonEfficiencyVsEta.pdf");


  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effVsNPV_GJetFlat_50ns_Loose, "Cut-based Loose", "LP");
  legend->AddEntry(effVsNPV_GJetFlat_50ns_Medium, "Cut-based Medium", "LP");
  legend->AddEntry(effVsNPV_GJetFlat_50ns_Tight, "Cut-based Tight", "LP");

  effVsNPV_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsNPV_GJetFlat_50ns_Loose->SetLineColor(kBlack);
  effVsNPV_GJetFlat_50ns_Loose->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
  effVsNPV_GJetFlat_50ns_Loose->GetYaxis()->SetTitle("Selection Efficiency");
  effVsNPV_GJetFlat_50ns_Loose->GetYaxis()->SetTitleOffset(1.2);

  effVsNPV_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsNPV_GJetFlat_50ns_Loose->SetLineColor(kRed);
  effVsNPV_GJetFlat_50ns_Medium->SetLineWidth(3);
  effVsNPV_GJetFlat_50ns_Medium->SetLineColor(kBlue);
  effVsNPV_GJetFlat_50ns_Tight->SetLineWidth(3);
  effVsNPV_GJetFlat_50ns_Tight->SetLineColor(kGreen+2);


  effVsNPV_GJetFlat_50ns_Loose->Draw("AP");
  effVsNPV_GJetFlat_50ns_Medium->Draw("Psame");
  effVsNPV_GJetFlat_50ns_Tight->Draw("Psame");
  effVsNPV_GJetFlat_50ns_Loose->GetXaxis()->SetRangeUser(0,40);
  
  legend->Draw();  
  cv->SaveAs("PhotonEfficiencyVsNPV.gif");
  cv->SaveAs("PhotonEfficiencyVsNPV.pdf");


  //****************************************************************************
  //Compare GJet with ttH
  //****************************************************************************
  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effVsPt_GJetFlat_50ns_Loose, "#gamma+Jet Flat Loose WP", "LP");
  legend->AddEntry(effVsPt_ttH_50ns_Loose, "ttH#rightarrow#gamma#gamma Loose WP", "LP");

  effVsPt_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_Loose->SetLineColor(kBlack);
  effVsPt_GJetFlat_50ns_Loose->GetXaxis()->SetTitle("Photon p_{T} [GeV/c]");
  effVsPt_GJetFlat_50ns_Loose->GetYaxis()->SetTitle("Selection Efficiency");
  effVsPt_GJetFlat_50ns_Loose->GetYaxis()->SetTitleOffset(1.2);

  effVsPt_GJetFlat_50ns_Loose->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_Loose->SetLineColor(kRed);
  effVsPt_ttH_50ns_Loose->SetLineWidth(3);
  effVsPt_ttH_50ns_Loose->SetLineColor(kBlue);

  effVsPt_GJetFlat_50ns_Loose->Draw("AP");
  effVsPt_ttH_50ns_Loose->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("PhotonEfficiencyVsPt_GJetVsTTH_Loose.gif");
  cv->SaveAs("PhotonEfficiencyVsPt_GJetVsTTH_Loose.pdf");



  cv = new TCanvas("cv","cv", 800,600);

  legend = new TLegend(0.50,0.34,0.90,0.54);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(effVsPt_GJetFlat_50ns_NotCloseToParton_Loose, "#gamma+Jet Flat Loose WP", "LP");
  legend->AddEntry(effVsPt_ttH_50ns_NotCloseToParton_Loose, "ttH#rightarrow#gamma#gamma Loose WP", "LP");

  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->SetLineColor(kBlack);
  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->GetXaxis()->SetTitle("Photon p_{T} [GeV/c]");
  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->GetYaxis()->SetTitle("Selection Efficiency");
  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->GetYaxis()->SetTitleOffset(1.2);

  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->SetLineWidth(3);
  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->SetLineColor(kRed);
  effVsPt_ttH_50ns_NotCloseToParton_Loose->SetLineWidth(3);
  effVsPt_ttH_50ns_NotCloseToParton_Loose->SetLineColor(kBlue);

  effVsPt_GJetFlat_50ns_NotCloseToParton_Loose->Draw("AP");
  effVsPt_ttH_50ns_NotCloseToParton_Loose->Draw("Psame");
  
  legend->Draw();  
  cv->SaveAs("PhotonEfficiencyVsPt_GJetVsTTH_NotCloseToParton_Loose.gif");
  cv->SaveAs("PhotonEfficiencyVsPt_GJetVsTTH_NotCloseToParton_Loose.pdf");




 return;


}



void MakeFastsimToFullSimCorrectionFactors() {
  TCanvas *cv =0;
  TLegend *legend =0;

  TFile *fileFullsimLoose = new TFile("Efficiency_GJetFlat_25ns_Loose_Fullsim.root","READ");
  TFile *fileFastsimLoose = new TFile("Efficiency_GJetFlat_25ns_Loose_Fastsim.root","READ");


  TH2F* histFullsimLoose = (TH2F*)fileFullsimLoose->Get("Efficiency_PtEta");
  TH2F* histFastsimLoose = (TH2F*)fileFastsimLoose->Get("Efficiency_PtEta");

  TH2F* histSFLoose = (TH2F*)histFullsimLoose->Clone("ElectronLoose_FastsimScaleFactor");
  histSFLoose->GetXaxis()->SetTitle("Electron p_{T} [GeV/c]");
  histSFLoose->GetYaxis()->SetTitle("Electron #eta");

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
  TFile *file = TFile::Open("PhotonEffFastsimToFullsimCorrectionFactors.root", "UPDATE");
  file->cd();
  file->WriteTObject(histSFLoose, "ElectronLoose_FastsimScaleFactor", "WriteDelete");  
  file->WriteTObject(histFullsimLoose, "ElectronEff_Loose_Fullsim", "WriteDelete");
  file->WriteTObject(histFastsimLoose, "ElectronEff_Loose_Fastsim", "WriteDelete");
  file->Close();
  delete file;      

}




//=== MAIN MACRO ================================================================================================= 

void ProducePhotonEfficiencyPlots(const string inputfile, int wp = 0, int option = -1, bool usePhotonNotNearParton = false, string label = "") {

   // plotPhotonEfficiency();
   // return;

  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  bool printdebug = false;

  TFile *inputFile = new TFile("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonEfficiencies/2016/EGM_PhotonLoose_SF.root","READ");
  TH2F *histSF = (TH2F*)inputFile->Get("EGamma_SF2D");

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************
  TH1F *histDenominatorPt = new TH1F ("histDenominatorPt",";Photon p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 200);
  TH1F *histNumeratorPt = new TH1F ("histNumeratorPt",";Photon p_{T} [GeV/c^{2}]; Number of Events", 50, 0 , 200);
  TH1F *histDenominatorEta = new TH1F ("histDenominatorEta",";Photon #eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histNumeratorEta = new TH1F ("histNumeratorEta",";Photon #eta; Number of Events", 50, -2.5 , 2.5);
  TH1F *histDenominatorPhi = new TH1F ("histDenominatorPhi",";Photon #phi; Number of Events", 50, 0 , 3.2);
  TH1F *histNumeratorPhi = new TH1F ("histNumeratorPhi",";Photon #phi; Number of Events", 50, 0 , 3.2);
  TH1F *histDenominatorRho = new TH1F ("histDenominatorRho",";Rho; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorRho = new TH1F ("histNumeratorRho",";Rho; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpv = new TH1F ("histDenominatorNpv",";Number of Reconstructed Primary Vertices; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpv = new TH1F ("histNumeratorNpv",";Number of Reconstructed Primary Vertices; Number of Events", 50, 0 , 100);
  TH1F *histDenominatorNpu = new TH1F ("histDenominatorNpu",";Number of Pileup Interactions; Number of Events", 50, 0 , 100);
  TH1F *histNumeratorNpu = new TH1F ("histNumeratorNpu",";Number of Pileup Interactions; Number of Events", 50, 0 , 100);

  // TH2F *histDenominatorPtEta = new TH2F ("histDenominatorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 56, 20 , 300, 100, 0.0, 2.5);
  // TH2F *histNumeratorPtEta = new TH2F ("histNumeratorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 56, 20 , 300, 100, 0.0, 2.5);
  //   TH2F *histDenominatorPtEta = new TH2F ("histDenominatorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 28, 20 , 300, 25, 0.0, 2.5);
   // TH2F *histNumeratorPtEta = new TH2F ("histNumeratorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", 28, 20 , 300, 25, 0.0, 2.5);

  //For SUSY/EGM Public Results
  const int NPtBins = 4;
  const int NEtaBins = 10;
  double ptBins[NPtBins+1] = {20, 35, 50, 90, 500};
  double etaBins[NEtaBins+1] = {-2.5, -2.0, -1.566, -1.4442, -0.8, 0, 0.8, 1.4442, 1.566, 2.0, 2.5};

  TH2F *histDenominatorPtEta = new TH2F ("histDenominatorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", NEtaBins, etaBins, NPtBins, ptBins);
  TH2F *histNumeratorPtEta = new TH2F ("histNumeratorPtEta",";Photon p_{T} [GeV/c] ; Photon #eta; Number of Events", NEtaBins, etaBins, NPtBins, ptBins);

  //*******************************************************************************************
  //READ file
  //*******************************************************************************************                
  PhotonTree *PhoTree = new PhotonTree;
  PhoTree->LoadTree(inputfile.c_str());
  PhoTree->InitTree(PhotonTree::kPhotonTreeLight);

  cout << "Total Entries: " << PhoTree->tree_->GetEntries() << "\n";
  int nentries = PhoTree->tree_->GetEntries();
  for(UInt_t ientry=0; ientry < PhoTree->tree_->GetEntries(); ientry++) {       	
    PhoTree->tree_->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //Cuts
    if (PhoTree->fPhoGenPt < 20) continue;
    if (abs(PhoTree->fPhoGenEta) > 2.5) continue;
    //if (abs(PhoTree->fPhoGenEta) <= 1.566) continue;
    //if (abs(PhoTree->fPhoGenEta) >= 1.4442) continue;

    //if (!(PhoTree->fPhoPt > 25)) continue;

    if (usePhotonNotNearParton) {
      if (PhoTree->fDRToClosestParton < 1.0) continue;
    }


    if (option==0) {

      if (PhoTree->fPhoGenPt < 20) continue;      
      

      //**** PT - ETA ****
      //For Internal Plot
      // histDenominatorPtEta->Fill(PhoTree->fPhoGenPt,fabs(PhoTree->fPhoGenEta));
      // if(PassSelection(PhoTree, wp)) {
      // 	histNumeratorPtEta->Fill(PhoTree->fPhoGenPt,fabs(PhoTree->fPhoGenEta));
      // }
      //For Public Plot
      histDenominatorPtEta->Fill(PhoTree->fPhoGenEta,PhoTree->fPhoGenPt);
      if(PassSelection(PhoTree, wp)) {
	histNumeratorPtEta->Fill(PhoTree->fPhoGenEta,PhoTree->fPhoGenPt);
      }									
      
  

      //**** PT ****
      histDenominatorPt->Fill(PhoTree->fPhoGenPt);

      //Numerator
      if(PassSelection(PhoTree, wp)) {
        histNumeratorPt->Fill(PhoTree->fPhoGenPt);        
      }


      //**** Eta ****
      if (fabs(PhoTree->fPhoGenPt) > 30) {
	histDenominatorEta->Fill(PhoTree->fPhoGenEta);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorEta->Fill(PhoTree->fPhoGenEta);        
	}

      }

      //**** Phi ****
      if (fabs(PhoTree->fPhoGenEta) < 2.5) {
	histDenominatorPhi->Fill(PhoTree->fPhoGenPhi);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorPhi->Fill(PhoTree->fPhoGenPhi);        
	}

      }

      //**** Rho ****
      if (fabs(PhoTree->fPhoGenEta) < 2.5) {
	histDenominatorRho->Fill(PhoTree->fRho);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorRho->Fill(PhoTree->fRho);        
	}

      }
      //**** Npv ****
      if (fabs(PhoTree->fPhoGenEta) < 2.5) {
	histDenominatorNpv->Fill(PhoTree->fNVertices);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorNpv->Fill(PhoTree->fNVertices);        
	}

      }

      // //**** Npu ****
      // if (fabs(PhoTree->fPhoGenEta) < 2.5) {
      //   histDenominatorNpu->Fill(PhoTree->);

      //   //Numerator
      //   if(PassSelection(PhoTree, wp)) {
      //     histNumeratorNpu->Fill(PhoTree->);        
      //   }

      // }

    } //end if option  == 0



    if (option==1) {
      //**** PT - ETA ****
      histDenominatorPtEta->Fill(PhoTree->fPhoPt,PhoTree->fPhoEta);
      if(PassSelection(PhoTree, wp)) {
	histNumeratorPtEta->Fill(PhoTree->fPhoPt,PhoTree->fPhoEta);
      }

      //**** PT ****
      histDenominatorPt->Fill(PhoTree->fPhoPt);

      //Numerator
      if(PassSelection(PhoTree, wp)) {
        histNumeratorPt->Fill(PhoTree->fPhoPt);        
      }


      //**** Eta ****
      if (fabs(PhoTree->fPhoPt) > 30) {
	histDenominatorEta->Fill(PhoTree->fPhoEta);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorEta->Fill(PhoTree->fPhoEta);        
	}

      }

      //**** Phi ****
      if (fabs(PhoTree->fPhoEta) < 2.5) {
	histDenominatorPhi->Fill(PhoTree->fPhoPhi);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorPhi->Fill(PhoTree->fPhoPhi);        
	}

      }

      //**** Rho ****
      if (fabs(PhoTree->fPhoEta) < 2.5) {
	histDenominatorRho->Fill(PhoTree->fRho);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorRho->Fill(PhoTree->fRho);        
	}

      }
      //**** Npv ****
      if (fabs(PhoTree->fPhoEta) < 2.5) {
	histDenominatorNpv->Fill(PhoTree->fNVertices);

	//Numerator
	if(PassSelection(PhoTree, wp)) {
	  histNumeratorNpv->Fill(PhoTree->fNVertices);        
	}
      }
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
  // Make Public Efficiency Plot
  // Apply Data to MC Corrections
  //==============================================================================================================
  TH2F *corrEffHist = new TH2F ("PhotonEfficiency",";Photon |#eta|; Photon E_{T} [GeV/c] ; Efficiency", 5, 0 , 5, 4, 0, 4);
  corrEffHist->GetXaxis()->SetTitle( "Photon |#eta|");
  corrEffHist->GetXaxis()->SetTitleSize( 0.05);
  corrEffHist->GetXaxis()->SetLabelSize( 0.05);
  corrEffHist->GetXaxis()->SetBinLabel(1, "(0.0, 0.8)");
  corrEffHist->GetXaxis()->SetBinLabel(2, "(0.8,1.4442)");
  corrEffHist->GetXaxis()->SetBinLabel(3, "(1.4442,1.566)");
  corrEffHist->GetXaxis()->SetBinLabel(4, "(1.566,2.0)");
  corrEffHist->GetXaxis()->SetBinLabel(5, "(2.0,2.5)");

  corrEffHist->GetYaxis()->SetTitle( "Photon E_{T} (GeV)");
  corrEffHist->GetYaxis()->SetTitleSize( 0.05);
  corrEffHist->GetYaxis()->SetLabelSize( 0.05);
  corrEffHist->GetYaxis()->SetBinLabel(1, "20 - 35");
  corrEffHist->GetYaxis()->SetBinLabel(2, "35 - 50");
  corrEffHist->GetYaxis()->SetBinLabel(3, "50 - 90");
  corrEffHist->GetYaxis()->SetBinLabel(4, "90 - 500");

  for (int i=1; i< corrEffHist->GetXaxis()->GetNbins()+1; i++) {
    for (int j=1; j< corrEffHist->GetYaxis()->GetNbins()+1; j++) {

      double tmpPt = efficiency_pteta->GetYaxis()->GetBinCenter(j);
      double tmpEta = 0.5;
      if (i==2) tmpEta = 1.221;
      if (i==3) tmpEta = 1.5;
      if (i==4) tmpEta = 1.75;
      if (i==5) tmpEta = 2.25;

      double EleVetoCorr = 1.0;
      double EleVetoCorrErr = 0.0;
      //CSEV Scale factors
      // if (fabs(tmpEta) < 1.4442) {
      // 	EleVetoCorr = 0.9983;
      // 	EleVetoCorrErr = 0.0119;
      // } else {
      // 	EleVetoCorr = 0.9875;
      // 	EleVetoCorrErr = 0.0044;
      // }

      //Pixel Veto Scale factors
      if (fabs(tmpEta) < 1.4442) {
      	EleVetoCorr = 0.9978;
      	EleVetoCorrErr = 0.0134;
      } else {
      	EleVetoCorr = 0.9931;
      	EleVetoCorrErr = 0.0245;
      }

      double corrEff = 0.5*( efficiency_pteta->GetBinContent(efficiency_pteta->GetXaxis()->FindFixBin(tmpEta), efficiency_pteta->GetYaxis()->FindFixBin(tmpPt))*
  	histSF->GetBinContent( histSF->GetXaxis()->FindFixBin(tmpEta), histSF->GetYaxis()->FindFixBin(tmpPt)) 
  			 +
  			 efficiency_pteta->GetBinContent(efficiency_pteta->GetXaxis()->FindFixBin(-1*tmpEta), efficiency_pteta->GetYaxis()->FindFixBin(tmpPt))*
  			 histSF->GetBinContent( histSF->GetXaxis()->FindFixBin(-1*tmpEta), histSF->GetYaxis()->FindFixBin(tmpPt)) 
  			 )
  	*EleVetoCorr;
      double corrEffErr = corrEff*sqrt( pow( efficiency_pteta->GetBinError(efficiency_pteta->GetXaxis()->FindFixBin(tmpEta), efficiency_pteta->GetYaxis()->FindFixBin(tmpPt))/efficiency_pteta->GetBinContent(efficiency_pteta->GetXaxis()->FindFixBin(tmpEta), efficiency_pteta->GetYaxis()->FindFixBin(tmpPt)),2) +
  					pow( histSF->GetBinError( histSF->GetXaxis()->FindFixBin(tmpEta), histSF->GetYaxis()->FindFixBin(tmpPt)) / histSF->GetBinContent( histSF->GetXaxis()->FindFixBin(tmpEta),  histSF->GetYaxis()->FindFixBin(tmpPt)), 2) + 
  					pow(EleVetoCorrErr / EleVetoCorr, 2));

      corrEffHist->SetBinContent(i,j,corrEff);
      corrEffHist->SetBinError(i,j,corrEffErr);

      //zero out the gap
      if (i==3) {
	corrEffHist->SetBinContent(i,j,0);
	corrEffHist->SetBinError(i,j,0);
      }

      cout << "Bin " << i << " " << j << " : " 
  	   << tmpEta << " " 
  	   << tmpPt << " | "
  	   << efficiency_pteta->GetXaxis()->FindFixBin(tmpEta) << " , " << efficiency_pteta->GetYaxis()->FindFixBin(tmpPt) << " | "
  	   << efficiency_pteta->GetBinContent(efficiency_pteta->GetXaxis()->FindFixBin(tmpEta),efficiency_pteta->GetYaxis()->FindFixBin(tmpPt)) << " +/- " << efficiency_pteta->GetBinError(efficiency_pteta->GetXaxis()->FindFixBin(tmpEta),efficiency_pteta->GetYaxis()->FindFixBin(tmpPt)) << " | " 
  	   << efficiency_pteta->GetBinContent(efficiency_pteta->GetXaxis()->FindFixBin(-1*tmpEta),efficiency_pteta->GetYaxis()->FindFixBin(tmpPt)) << " +/- " << efficiency_pteta->GetBinError(efficiency_pteta->GetXaxis()->FindFixBin(-1*tmpEta),efficiency_pteta->GetYaxis()->FindFixBin(tmpPt)) << " | " 
  	   << histSF->GetXaxis()->FindFixBin(efficiency_pteta->GetXaxis()->GetBinCenter(i)) << ","
  	   << histSF->GetYaxis()->FindFixBin(efficiency_pteta->GetYaxis()->GetBinCenter(j)) << " | " 
  	   << histSF->GetBinContent( histSF->GetXaxis()->FindFixBin(efficiency_pteta->GetXaxis()->GetBinCenter(i)),
  	    			     histSF->GetYaxis()->FindFixBin(efficiency_pteta->GetYaxis()->GetBinCenter(j))) 
  	   << " +/- " 
  	   << histSF->GetBinError( histSF->GetXaxis()->FindFixBin(efficiency_pteta->GetXaxis()->GetBinCenter(i)),
  				   histSF->GetYaxis()->FindFixBin(efficiency_pteta->GetYaxis()->GetBinCenter(j))) 
  	   << " | " 
  	   << EleVetoCorr << " +/- " << EleVetoCorrErr << " | "
  	   << corrEff << " +/- " << corrEffErr << " "
  	   << "\n";
    }
  }

 

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

  //TStyle myStyle;
  cv = new TCanvas("cv","cv",800,600);
  gStyle->SetPaintTextFormat(".3f");
  cv->SetBottomMargin(0.12);
  cv->SetLeftMargin(0.15);
  cv->SetRightMargin(0.12);
  corrEffHist->SetTitle("Efficiency");
  corrEffHist->SetMarkerSize(1.5);
  corrEffHist->SetStats(0);
  corrEffHist->GetYaxis()->SetTitleOffset(1.5);
  corrEffHist->GetZaxis()->SetTitleOffset(1.0);
  corrEffHist->GetZaxis()->SetTitleSize(0.0);
  corrEffHist->SetMinimum(0.0);
  corrEffHist->SetMaximum(1.0);
  corrEffHist->Draw("colztexte1");
  cv->SaveAs("PhotonEfficiency.pdf");

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
  // file->WriteTObject(corrEffHist, "PhotonEfficiency", "WriteDelete");

  file->Close();
  delete file;       

  file = TFile::Open("PhotonEfficiencyPublic.root","RECREATE");
  file->cd();
  file->WriteTObject(corrEffHist, "PhotonEfficiency", "WriteDelete");
  file->Close();
  delete file;

}




void MakePhotonEfficiencyPlots( int option = 0) {

  if (option == 1) {
      
    //Fullsim
    //ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/V3p8/PhotonNtuple_PromptGenLevel_Fullsim_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_25ns.root", 0, 0, true, "GJetFlat_25ns_Loose_Fullsim");
    // ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/V1p17/PhotonNtuple_PromptGenLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_25ns.root", 1, 0, true, "GJetFlat_25ns_Medium_NotCloseToParton");
    // ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/V1p17/PhotonNtuple_PromptGenLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_25ns.root", 2, 0, true, "GJetFlat_25ns_Tight_NotCloseToParton");

    //Fastsim
    //ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/V3p8/PhotonNtuple_PromptGenLevel_Fastsim_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_25ns.root", 0, 0, true, "GJetFlat_25ns_Loose_Fastsim");

    //Fake Photons
    // ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/V1p17/PhotonNtuple_Fake_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_25ns.root", 0, 0, false, "GJetFlat_Fake_25ns_Loose");
    // ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/V1p17/PhotonNtuple_Fake_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_25ns.root", 1, 0, false, "GJetFlat_Fake_25ns_Medium");
    // ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/V1p17/PhotonNtuple_Fake_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_25ns.root", 2, 0, false, "GJetFlat_Fake_25ns_Tight");

    
    //Used to Produce Public Photon Efficiency Plot for ICHEP dataset
    // ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/V3p4/PhotonNtuple_PromptGenLevel_HiggsCombined.root", 0, 0, false, "Photon_LooseID");


    //Used to Produce Public Photon Efficiency Plot for Moriond dataset
    ProducePhotonEfficiencyPlots("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/V3p8/PhotonNtuple_PromptGenLevel_Fullsim_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_25ns.root", 0, 0, false, "Photon_LooseID");
  
  }	



  // plotPhotonEfficiency();
  //MakeFastsimToFullSimCorrectionFactors();

}
