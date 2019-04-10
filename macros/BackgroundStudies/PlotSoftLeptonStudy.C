
// root -l /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/macros/BackgroundStudies/PlotSoftLeptonStudy.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonStudy/VetoLeptonStudy_Bkg_25ns_weighted.root","Bkg_SoftLeptonMultiJet",7)'
// root -l /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/macros/BackgroundStudies/PlotSoftLeptonStudy.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonStudy/VetoLeptonStudy_Bkg_25ns_weighted.root","Bkg_MultiJet",8)'
// root -l /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/macros/BackgroundStudies/PlotSoftLeptonStudy.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonStudy/VetoLeptonStudy_TTJets_25ns_weighted.root","TTJets_MultiJet",7)'
// root -l /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/macros/BackgroundStudies/PlotSoftLeptonStudy.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonStudy/VetoLeptonStudy_ZJetsToNuNu_25ns_weighted.root","ZToNuNu_MultiJet",7)'
// root -l /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/macros/BackgroundStudies/PlotSoftLeptonStudy.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonStudy/VetoLeptonStudy_DYJetsToLL_25ns_weighted.root","DYToLL_MultiJet",7)'
// root -l /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/macros/BackgroundStudies/PlotSoftLeptonStudy.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonStudy/VetoLeptonStudy_WJetsToLNu_25ns_weighted.root","WToLNu_MultiJet",7)'
// root -l /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/macros/BackgroundStudies/PlotSoftLeptonStudy.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonStudy/VetoLeptonStudy_SMS-T1tttt_2J_mGl-1500_mLSP-100_25ns_weighted.root","T1tttt_MultiJet",7)'
// root -l /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/macros/BackgroundStudies/PlotSoftLeptonStudy.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/VetoLeptonStudy/VetoLeptonStudy_SMS-T1bbbb_2J_mGl-1500_mLSP-100_25ns_weighted.root","T1bbbb_MultiJet",7)'

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <vector>
#include <map>
#include <iostream>


 const Int_t NComponents = 3;
 Int_t Colors[NComponents] = { kRed, kGreen+2, kBlue };
string LegendLabels[NComponents] = { "0", "1", "2"}; 

//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
TH1F* NormalizeHist(TH1F *originalHist) {
  TH1F* hist = (TH1F*)originalHist->Clone((string(originalHist->GetName())+"_normalized").c_str());
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return hist;
}


//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void PlotSoftLeptonStudy ( string inputfilemame, string Label, Int_t Option)
{

  double intLumi = 5000; //in units of pb^-1
  string label = "";
  if (Label != "") label = "_" + Label;


  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  vector<TH1F*> MR_nGenTaus;
  vector<TH1F*> Rsq_nGenTaus;
  vector<TH1F*> MR_nGenElectrons;
  vector<TH1F*> Rsq_nGenElectrons;
  vector<TH1F*> MR_nGenMuons;
  vector<TH1F*> Rsq_nGenMuons;
  vector<TH1F*> MR_nGenElectronsAndMuons;
  vector<TH1F*> Rsq_nGenElectronsAndMuons;
  vector<TH1F*> MR_nGenElectronsAndMuonsIncludingTauDecays;
  vector<TH1F*> Rsq_nGenElectronsAndMuonsIncludingTauDecays;
 

  for (UInt_t i = 1; i <= NComponents; ++i) {
    MR_nGenTaus.push_back( new TH1F( Form("MR_nGenTaus_%d",i), ";M_{R} [GeV/c^{2}];Number of Events", 40, 0, 4000));
    MR_nGenTaus[i-1]->SetFillColor(Colors[i-1]);
    MR_nGenTaus[i-1]->SetLineColor(Colors[i-1]);
    MR_nGenTaus[i-1]->SetStats(false);

    Rsq_nGenTaus.push_back( new TH1F( Form("Rsq_nGenTaus_%d",i), ";R^{2} ;Number of Events", 20, 0, 1.5));
    Rsq_nGenTaus[i-1]->SetFillColor(Colors[i-1]);
    Rsq_nGenTaus[i-1]->SetLineColor(Colors[i-1]);
    Rsq_nGenTaus[i-1]->SetStats(false);

    MR_nGenElectrons.push_back( new TH1F( Form("MR_nGenElectrons_%d",i), ";M_{R} [GeV/c^{2}];Number of Events", 40, 0, 4000));
    MR_nGenElectrons[i-1]->SetFillColor(Colors[i-1]);
    MR_nGenElectrons[i-1]->SetLineColor(Colors[i-1]);
    MR_nGenElectrons[i-1]->SetStats(false);

    Rsq_nGenElectrons.push_back( new TH1F( Form("Rsq_nGenElectrons_%d",i), ";R^{2} ;Number of Events", 20, 0, 1.5));
    Rsq_nGenElectrons[i-1]->SetFillColor(Colors[i-1]);
    Rsq_nGenElectrons[i-1]->SetLineColor(Colors[i-1]);
    Rsq_nGenElectrons[i-1]->SetStats(false);

    MR_nGenMuons.push_back( new TH1F( Form("MR_nGenMuons_%d",i), ";M_{R} [GeV/c^{2}];Number of Events", 40, 0, 4000));
    MR_nGenMuons[i-1]->SetFillColor(Colors[i-1]);
    MR_nGenMuons[i-1]->SetLineColor(Colors[i-1]);
    MR_nGenMuons[i-1]->SetStats(false);

    Rsq_nGenMuons.push_back( new TH1F( Form("Rsq_nGenMuons_%d",i), ";R^{2} ;Number of Events", 20, 0, 1.5));
    Rsq_nGenMuons[i-1]->SetFillColor(Colors[i-1]);
    Rsq_nGenMuons[i-1]->SetLineColor(Colors[i-1]);
    Rsq_nGenMuons[i-1]->SetStats(false);

    MR_nGenElectronsAndMuons.push_back( new TH1F( Form("MR_nGenElectronsAndMuons_%d",i), ";M_{R} [GeV/c^{2}];Number of Events", 40, 0, 4000));
    MR_nGenElectronsAndMuons[i-1]->SetFillColor(Colors[i-1]);
    MR_nGenElectronsAndMuons[i-1]->SetLineColor(Colors[i-1]);
    MR_nGenElectronsAndMuons[i-1]->SetStats(false);

    Rsq_nGenElectronsAndMuons.push_back( new TH1F( Form("Rsq_nGenElectronsAndMuons_%d",i), ";R^{2} ;Number of Events", 20, 0, 1.5));
    Rsq_nGenElectronsAndMuons[i-1]->SetFillColor(Colors[i-1]);
    Rsq_nGenElectronsAndMuons[i-1]->SetLineColor(Colors[i-1]);
    Rsq_nGenElectronsAndMuons[i-1]->SetStats(false);

    MR_nGenElectronsAndMuonsIncludingTauDecays.push_back( new TH1F( Form("MR_nGenElectronsAndMuonsIncludingTauDecays_%d",i), ";M_{R} [GeV/c^{2}];Number of Events", 40, 0, 4000));
    MR_nGenElectronsAndMuonsIncludingTauDecays[i-1]->SetFillColor(Colors[i-1]);
    MR_nGenElectronsAndMuonsIncludingTauDecays[i-1]->SetLineColor(Colors[i-1]);
    MR_nGenElectronsAndMuonsIncludingTauDecays[i-1]->SetStats(false);

    Rsq_nGenElectronsAndMuonsIncludingTauDecays.push_back( new TH1F( Form("Rsq_nGenElectronsAndMuonsIncludingTauDecays_%d",i), ";R^{2} ;Number of Events", 20, 0, 1.5));
    Rsq_nGenElectronsAndMuonsIncludingTauDecays[i-1]->SetFillColor(Colors[i-1]);
    Rsq_nGenElectronsAndMuonsIncludingTauDecays[i-1]->SetLineColor(Colors[i-1]);
    Rsq_nGenElectronsAndMuonsIncludingTauDecays[i-1]->SetStats(false);


 }

  //*******************************************************************************************
  //Define Counts
  //*******************************************************************************************
  double NEvents = 0;
  double NEvents_ZeroGenLeptons = 0;
  double NEvents_ZeroGenTaus = 0;
  double NEvents_ZeroGenMuons = 0;
  double NEvents_ZeroGenElectrons = 0;
  double NEvents_ZeroGenElectronsAndMuons = 0;
  double NEvents_ZeroGenElectronsAndMuonsIncludingTauDecays = 0;


  //*******************************************************************************************
  //Read Ntuple
  //*******************************************************************************************
  TFile* inputFile = new TFile(inputfilemame.c_str(),"read");
  assert(inputFile);
  TTree* tree = (TTree*)inputFile->Get("RazorInclusive");
 
  double weight = 0;
  int box = -1;
  int nGenElectrons = 0;
  int nGenMuons = 0;
  int nGenTaus = 0;
  int nGenTauElectrons = 0;
  int nGenTauMuons = 0;
  float MR = 0;
  float Rsq = 0;

  tree->SetBranchAddress("weight",&weight);
  tree->SetBranchAddress("box",&box);
  tree->SetBranchAddress("nGenElectrons",&nGenElectrons);
  tree->SetBranchAddress("nGenMuons",&nGenMuons);
  tree->SetBranchAddress("nGenTaus",&nGenTaus);
  tree->SetBranchAddress("nGenTauElectrons",&nGenTauElectrons);
  tree->SetBranchAddress("nGenTauMuons",&nGenTauMuons);
  tree->SetBranchAddress("MR",&MR);
  tree->SetBranchAddress("Rsq",&Rsq);

  cout << "Bkg Total Events: " << tree->GetEntries() << "\n";
  for (int n=0;n<tree->GetEntries();n++) { 
    
    tree->GetEntry(n);
    if (n % 1000000 == 0) cout << "Processing Event " << n << "\n";       

    //Box Options
    if (Option >= 0 && Option < 10 && box != Option) continue;

    if (Option == 135 && !(box == 3 || box == 5)) continue;
    if (Option == 178 && !(box == 7 || box == 8)) continue;

    //apply baseline cuts
    if (!(MR > 300 && Rsq > 0.1)) continue;

    //counts
    if (Rsq > 0.25 ) {
    //if (MR > 1000) {
    //if ( 0 == 0) { 
      NEvents += intLumi*weight;
      if (nGenTaus == 0)  NEvents_ZeroGenTaus += intLumi*weight;
      if (nGenElectrons == 0)  NEvents_ZeroGenElectrons += intLumi*weight;
      if (nGenMuons == 0)  NEvents_ZeroGenMuons += intLumi*weight;
      if (nGenElectrons + nGenMuons + nGenTaus == 0)  NEvents_ZeroGenLeptons += intLumi*weight;
      if (nGenElectrons + nGenMuons == 0)  NEvents_ZeroGenElectronsAndMuons += intLumi*weight;
      if (nGenElectrons + nGenMuons + nGenTauElectrons + nGenTauMuons == 0)  NEvents_ZeroGenElectronsAndMuonsIncludingTauDecays += intLumi*weight;
    }

    //Fill MR plots
    if (Rsq > 0.25 && MR > 1000) {
      if (nGenTaus == 0) MR_nGenTaus[0]->Fill(MR, intLumi*weight);
      else if (nGenTaus == 1) MR_nGenTaus[1]->Fill(MR, intLumi*weight);
      else if (nGenTaus >= 2) MR_nGenTaus[2]->Fill(MR, intLumi*weight);

      if (nGenElectrons == 0) MR_nGenElectrons[0]->Fill(MR, intLumi*weight);
      else if (nGenElectrons == 1) MR_nGenElectrons[1]->Fill(MR, intLumi*weight);
      else if (nGenElectrons >= 2) MR_nGenElectrons[2]->Fill(MR, intLumi*weight);

      if (nGenMuons == 0) MR_nGenMuons[0]->Fill(MR, intLumi*weight);
      else if (nGenMuons == 1) MR_nGenMuons[1]->Fill(MR, intLumi*weight);
      else if (nGenMuons >= 2) MR_nGenMuons[2]->Fill(MR, intLumi*weight);

      if (nGenElectrons + nGenMuons == 0) MR_nGenElectronsAndMuons[0]->Fill(MR, intLumi*weight);
      else if (nGenElectrons + nGenMuons == 1) MR_nGenElectronsAndMuons[1]->Fill(MR, intLumi*weight);
      else if (nGenElectrons + nGenMuons >= 2) MR_nGenElectronsAndMuons[2]->Fill(MR, intLumi*weight);

      if (nGenElectrons + nGenMuons + nGenTauElectrons + nGenTauMuons == 0) MR_nGenElectronsAndMuonsIncludingTauDecays[0]->Fill(MR, intLumi*weight);
      else if (nGenElectrons + nGenMuons + nGenTauElectrons + nGenTauMuons == 1) MR_nGenElectronsAndMuonsIncludingTauDecays[1]->Fill(MR, intLumi*weight);
      else if (nGenElectrons + nGenMuons + nGenTauElectrons + nGenTauMuons >= 2) MR_nGenElectronsAndMuonsIncludingTauDecays[2]->Fill(MR, intLumi*weight);   

    } 

    //Fill Rsq plots
    if (MR > 1000 && Rsq > 0.25) {
      if (nGenTaus == 0) Rsq_nGenTaus[0]->Fill(Rsq, intLumi*weight);
      else if (nGenTaus == 1) Rsq_nGenTaus[1]->Fill(Rsq, intLumi*weight);
      else if (nGenTaus >= 2) Rsq_nGenTaus[2]->Fill(Rsq, intLumi*weight);

      if (nGenElectrons == 0) Rsq_nGenElectrons[0]->Fill(Rsq, intLumi*weight);
      else if (nGenElectrons == 1) Rsq_nGenElectrons[1]->Fill(Rsq, intLumi*weight);
      else if (nGenElectrons >= 2) Rsq_nGenElectrons[2]->Fill(Rsq, intLumi*weight);

      if (nGenMuons == 0) Rsq_nGenMuons[0]->Fill(Rsq, intLumi*weight);
      else if (nGenMuons == 1) Rsq_nGenMuons[1]->Fill(Rsq, intLumi*weight);
      else if (nGenMuons >= 2) Rsq_nGenMuons[2]->Fill(Rsq, intLumi*weight);

      if (nGenElectrons + nGenMuons == 0) Rsq_nGenElectronsAndMuons[0]->Fill(Rsq, intLumi*weight);
      else if (nGenElectrons + nGenMuons == 1) Rsq_nGenElectronsAndMuons[1]->Fill(Rsq, intLumi*weight);
      else if (nGenElectrons + nGenMuons >= 2) Rsq_nGenElectronsAndMuons[2]->Fill(Rsq, intLumi*weight);

      if (nGenElectrons + nGenMuons + nGenTauElectrons + nGenTauMuons  == 0) Rsq_nGenElectronsAndMuonsIncludingTauDecays[0]->Fill(Rsq, intLumi*weight);
      else if (nGenElectrons + nGenMuons + nGenTauElectrons + nGenTauMuons  == 1) Rsq_nGenElectronsAndMuonsIncludingTauDecays[1]->Fill(Rsq, intLumi*weight);
      else if (nGenElectrons + nGenMuons + nGenTauElectrons + nGenTauMuons  >= 2) Rsq_nGenElectronsAndMuonsIncludingTauDecays[2]->Fill(Rsq, intLumi*weight);   
    }   
  }
  
  


  //*******************************************************************************************
  //Draw Plots
  //*******************************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;


  //*******************************************************************************************
  //MR_nGenTaus
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMR_nGenTaus = new THStack();
  for (Int_t i = MR_nGenTaus.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << MR_nGenTaus[i]->GetSumOfWeights() << "\n";
    if ( MR_nGenTaus[i]->Integral() > 0) {
      stackMR_nGenTaus->Add(MR_nGenTaus[i]);
    }
  }
  for (Int_t i = 0 ; i < MR_nGenTaus.size(); ++i) {
    legend->AddEntry(MR_nGenTaus[i],(LegendLabels[i]+" GenTaus").c_str(), "F");
  }

  stackMR_nGenTaus->Draw();
  stackMR_nGenTaus->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR_nGenTaus->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackMR_nGenTaus->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR_nGenTaus->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs(Form("MR_nGenTaus%s.gif",label.c_str()));

  //*******************************************************************************************
  //Rsq_nGenTaus
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackRsq_nGenTaus = new THStack();
  for (Int_t i = Rsq_nGenTaus.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << Rsq_nGenTaus[i]->GetSumOfWeights() << "\n";
    if ( Rsq_nGenTaus[i]->Integral() > 0) {
      stackRsq_nGenTaus->Add(Rsq_nGenTaus[i]);
    }
  }
  for (Int_t i = 0 ; i < Rsq_nGenTaus.size(); ++i) {
    legend->AddEntry(Rsq_nGenTaus[i],(LegendLabels[i]+" GenTaus").c_str(), "F");
  }

  stackRsq_nGenTaus->Draw();
  stackRsq_nGenTaus->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackRsq_nGenTaus->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackRsq_nGenTaus->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackRsq_nGenTaus->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs(Form("Rsq_nGenTaus%s.gif",label.c_str()));

  //*******************************************************************************************
  //MR_nGenElectrons
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMR_nGenElectrons = new THStack();
  for (Int_t i = MR_nGenElectrons.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << MR_nGenElectrons[i]->GetSumOfWeights() << "\n";
    if ( MR_nGenElectrons[i]->Integral() > 0) {
      stackMR_nGenElectrons->Add(MR_nGenElectrons[i]);
    }
  }
  for (Int_t i = 0 ; i < MR_nGenElectrons.size(); ++i) {
    legend->AddEntry(MR_nGenElectrons[i],(LegendLabels[i]+" GenElectrons").c_str(), "F");
  }

  stackMR_nGenElectrons->Draw();
  stackMR_nGenElectrons->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR_nGenElectrons->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackMR_nGenElectrons->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR_nGenElectrons->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs(Form("MR_nGenElectrons%s.gif",label.c_str()));

  //*******************************************************************************************
  //Rsq_nGenElectrons
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackRsq_nGenElectrons = new THStack();
  for (Int_t i = Rsq_nGenElectrons.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << Rsq_nGenElectrons[i]->GetSumOfWeights() << "\n";
    if ( Rsq_nGenElectrons[i]->Integral() > 0) {
      stackRsq_nGenElectrons->Add(Rsq_nGenElectrons[i]);
    }
  }
  for (Int_t i = 0 ; i < Rsq_nGenElectrons.size(); ++i) {
    legend->AddEntry(Rsq_nGenElectrons[i],(LegendLabels[i]+" GenElectrons").c_str(), "F");
  }

  stackRsq_nGenElectrons->Draw();
  stackRsq_nGenElectrons->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackRsq_nGenElectrons->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackRsq_nGenElectrons->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackRsq_nGenElectrons->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs(Form("Rsq_nGenElectrons%s.gif",label.c_str()));

  //*******************************************************************************************
  //MR_nGenMuons
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMR_nGenMuons = new THStack();
  for (Int_t i = MR_nGenMuons.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << MR_nGenMuons[i]->GetSumOfWeights() << "\n";
    if ( MR_nGenMuons[i]->Integral() > 0) {
      stackMR_nGenMuons->Add(MR_nGenMuons[i]);
    }
  }
  for (Int_t i = 0 ; i < MR_nGenMuons.size(); ++i) {
    legend->AddEntry(MR_nGenMuons[i],(LegendLabels[i]+" GenMuons").c_str(), "F");
  }

  stackMR_nGenMuons->Draw();
  stackMR_nGenMuons->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR_nGenMuons->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackMR_nGenMuons->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR_nGenMuons->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs(Form("MR_nGenMuons%s.gif",label.c_str()));

  //*******************************************************************************************
  //Rsq_nGenMuons
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackRsq_nGenMuons = new THStack();
  for (Int_t i = Rsq_nGenMuons.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << Rsq_nGenMuons[i]->GetSumOfWeights() << "\n";
    if ( Rsq_nGenMuons[i]->Integral() > 0) {
      stackRsq_nGenMuons->Add(Rsq_nGenMuons[i]);
    }
  }
  for (Int_t i = 0 ; i < Rsq_nGenMuons.size(); ++i) {
    legend->AddEntry(Rsq_nGenMuons[i],(LegendLabels[i]+" GenMuons").c_str(), "F");
  }

  stackRsq_nGenMuons->Draw();
  stackRsq_nGenMuons->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackRsq_nGenMuons->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackRsq_nGenMuons->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackRsq_nGenMuons->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs(Form("Rsq_nGenMuons%s.gif",label.c_str()));



  //*******************************************************************************************
  //MR_nGenElectronsAndMuons
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMR_nGenElectronsAndMuons = new THStack();
  for (Int_t i = MR_nGenElectronsAndMuons.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << MR_nGenElectronsAndMuons[i]->GetSumOfWeights() << "\n";
    if ( MR_nGenElectronsAndMuons[i]->Integral() > 0) {
      stackMR_nGenElectronsAndMuons->Add(MR_nGenElectronsAndMuons[i]);
    }
  }
  for (Int_t i = 0 ; i < MR_nGenElectronsAndMuons.size(); ++i) {
    legend->AddEntry(MR_nGenElectronsAndMuons[i],(LegendLabels[i]+" GenElectronsAndMuons").c_str(), "F");
  }

  stackMR_nGenElectronsAndMuons->Draw();
  stackMR_nGenElectronsAndMuons->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR_nGenElectronsAndMuons->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackMR_nGenElectronsAndMuons->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR_nGenElectronsAndMuons->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs(Form("MR_nGenElectronsAndMuons%s.gif",label.c_str()));

  //*******************************************************************************************
  //Rsq_nGenElectronsAndMuons
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackRsq_nGenElectronsAndMuons = new THStack();
  for (Int_t i = Rsq_nGenElectronsAndMuons.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << Rsq_nGenElectronsAndMuons[i]->GetSumOfWeights() << "\n";
    if ( Rsq_nGenElectronsAndMuons[i]->Integral() > 0) {
      stackRsq_nGenElectronsAndMuons->Add(Rsq_nGenElectronsAndMuons[i]);
    }
  }
  for (Int_t i = 0 ; i < Rsq_nGenElectronsAndMuons.size(); ++i) {
    legend->AddEntry(Rsq_nGenElectronsAndMuons[i],(LegendLabels[i]+" GenElectronsAndMuons").c_str(), "F");
  }

  stackRsq_nGenElectronsAndMuons->Draw();
  stackRsq_nGenElectronsAndMuons->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackRsq_nGenElectronsAndMuons->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackRsq_nGenElectronsAndMuons->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackRsq_nGenElectronsAndMuons->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs(Form("Rsq_nGenElectronsAndMuons%s.gif",label.c_str()));




  //*******************************************************************************************
  //MR_nGenElectronsAndMuonsIncludingTauDecays
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMR_nGenElectronsAndMuonsIncludingTauDecays = new THStack();
  for (Int_t i = MR_nGenElectronsAndMuonsIncludingTauDecays.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << MR_nGenElectronsAndMuonsIncludingTauDecays[i]->GetSumOfWeights() << "\n";
    if ( MR_nGenElectronsAndMuonsIncludingTauDecays[i]->Integral() > 0) {
      stackMR_nGenElectronsAndMuonsIncludingTauDecays->Add(MR_nGenElectronsAndMuonsIncludingTauDecays[i]);
    }
  }
  for (Int_t i = 0 ; i < MR_nGenElectronsAndMuonsIncludingTauDecays.size(); ++i) {
    legend->AddEntry(MR_nGenElectronsAndMuonsIncludingTauDecays[i],(LegendLabels[i]+" GenElectronsAndMuonsIncludingTauDecays").c_str(), "F");
  }

  stackMR_nGenElectronsAndMuonsIncludingTauDecays->Draw();
  stackMR_nGenElectronsAndMuonsIncludingTauDecays->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR_nGenElectronsAndMuonsIncludingTauDecays->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackMR_nGenElectronsAndMuonsIncludingTauDecays->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR_nGenElectronsAndMuonsIncludingTauDecays->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs(Form("MR_nGenElectronsAndMuonsIncludingTauDecays%s.gif",label.c_str()));

  //*******************************************************************************************
  //Rsq_nGenElectronsAndMuonsIncludingTauDecays
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackRsq_nGenElectronsAndMuonsIncludingTauDecays = new THStack();
  for (Int_t i = Rsq_nGenElectronsAndMuonsIncludingTauDecays.size()-1; i >= 0; i--) {
    cout << "sample " << i << " : " << Rsq_nGenElectronsAndMuonsIncludingTauDecays[i]->GetSumOfWeights() << "\n";
    if ( Rsq_nGenElectronsAndMuonsIncludingTauDecays[i]->Integral() > 0) {
      stackRsq_nGenElectronsAndMuonsIncludingTauDecays->Add(Rsq_nGenElectronsAndMuonsIncludingTauDecays[i]);
    }
  }
  for (Int_t i = 0 ; i < Rsq_nGenElectronsAndMuonsIncludingTauDecays.size(); ++i) {
    legend->AddEntry(Rsq_nGenElectronsAndMuonsIncludingTauDecays[i],(LegendLabels[i]+" GenElectronsAndMuonsIncludingTauDecays").c_str(), "F");
  }

  stackRsq_nGenElectronsAndMuonsIncludingTauDecays->Draw();
  stackRsq_nGenElectronsAndMuonsIncludingTauDecays->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackRsq_nGenElectronsAndMuonsIncludingTauDecays->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackRsq_nGenElectronsAndMuonsIncludingTauDecays->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackRsq_nGenElectronsAndMuonsIncludingTauDecays->GetHists()->At(0)))->GetYaxis()->GetTitle());
  legend->Draw();
  cv->SaveAs(Form("Rsq_nGenElectronsAndMuonsIncludingTauDecays%s.gif",label.c_str()));



  //*******************************************************************************************
  //Summarize Counts
  //*******************************************************************************************
  cout << "Fraction of Events with 0 GenElectrons : " << NEvents_ZeroGenElectrons / NEvents << "\n";
  cout << "Fraction of Events with 0 GenMuons : " << NEvents_ZeroGenMuons / NEvents << "\n";
  cout << "Fraction of Events with 0 GenTaus : " << NEvents_ZeroGenTaus / NEvents << "\n";
  cout << "Fraction of Events with 0 GenElectronsAndMuons : " << NEvents_ZeroGenElectronsAndMuons / NEvents << "\n";
  cout << "Fraction of Events with 0 GenElectronsAndMuonsIncludingTauDecays : " << NEvents_ZeroGenElectronsAndMuonsIncludingTauDecays / NEvents << "\n";
  cout << "Fraction of Events with 0 GenLeptons : " << NEvents_ZeroGenLeptons / NEvents << "\n";

  

}
