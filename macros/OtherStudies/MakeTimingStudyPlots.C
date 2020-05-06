
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraph.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <vector>
#include <map>
#include <iostream>

const Int_t NComponents = 10;
int color[NComponents] = {kRed, kGreen+2, kBlue, kViolet, kAzure+10, kGray, kOrange+1, kGray+3, kBlack, kBlack};


TH1F* makeNewHist( string title, string axisLabel, int nbins, double rangeLow, double rangeHigh, int color, bool isData) {

  TH1F *newhist = new TH1F( title.c_str(), axisLabel.c_str(), nbins, rangeLow, rangeHigh);
  newhist->SetFillStyle(0);
  newhist->SetLineWidth(2); 
  newhist->SetLineColor(color);    
  newhist->SetStats(false);    
  newhist->Sumw2();
  return newhist;

}

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

double deltaPhi(double phi1, double phi2) {
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
    dphi -= TMath::TwoPi();
  while (dphi <= -TMath::Pi())
    dphi += TMath::TwoPi();
  return dphi;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1,phi2);
  double deta = eta1 - eta2;
  return sqrt( dphi*dphi + deta*deta);
}

// double getMuonEff( double pt, double eta, TH2F* muonEffHist, bool moreSeverePUCondition, bool WithBarrelTiming = false, bool WithEndcapTiming = false) {
//   double eff = 1.0;

//   eff = muonEffHist->GetBinContent( muonEffHist->GetXaxis()->FindFixBin( fmin( pt, 199.9) ), 
// 				    muonEffHist->GetYaxis()->FindFixBin( fmin( fabs(eta), 2.39) ) );
  
//   double pileupDegradationFactor = 0.82/0.95;
//   if (WithBarrelTiming && fabs(eta) < 1.5) pileupDegradationFactor = 0.92/0.95;
//   if (WithEndcapTiming && fabs(eta) > 1.5 && fabs(eta) < 2.4) pileupDegradationFactor = 0.92/0.95;

//   if (moreSeverePUCondition) {
//     pileupDegradationFactor = 0.76/0.95;
//     if (WithBarrelTiming && fabs(eta) < 1.5) pileupDegradationFactor = 0.90/0.95;
//     if (WithEndcapTiming && fabs(eta) > 1.5 && fabs(eta) < 2.4) pileupDegradationFactor = 0.90/0.95;    
//   }

//   eff = eff * pileupDegradationFactor;
//   return eff;
// }

double computeAvergeMuonEff( double pt, double eta, TH2F* muonEffHist, double linearDensity, bool WithBarrelTiming = false, bool WithEndcapTiming = false, bool WithHGCTiming = false) {

  double eff = 1.0;

  eff = muonEffHist->GetBinContent( muonEffHist->GetXaxis()->FindFixBin( fmin( pt, 199.9) ), 
				    muonEffHist->GetYaxis()->FindFixBin( fmin( fabs(eta), 2.39) ) );
  
  // double pileupDegradationFactorWithoutTiming = 0.95 - (linearDensity - 0.3)*(0.95-0.79)/(1.5-0.3);
  // double pileupDegradationFactorWithTiming = 0.95 - (linearDensity - 0.3)*(0.95-0.91)/(1.5-0.3);
  
  //Use a flat line for new estimates which is averaged for the 200PU scenario 
  double pileupDegradationFactorWithoutTiming = 0.90;
  double pileupDegradationFactorWithTiming = 0.95; 
  double pileupDegradationFactor = pileupDegradationFactorWithoutTiming;

  if (fabs(eta) < 1.5) {
    if (WithBarrelTiming) pileupDegradationFactor = 0.95;
    else pileupDegradationFactor = 0.90;    
  } else if (fabs(eta) >= 1.5 && fabs(eta) < 2.4) {
    if (WithEndcapTiming) pileupDegradationFactor = 0.97;  
    else if (WithHGCTiming) pileupDegradationFactor = 0.95;  
    else pileupDegradationFactor = 0.90;    
  }
  
  // if (WithHGCTiming) {
  //   // if (WithEndcapTiming && fabs(eta) > 1.5 && fabs(eta) < 2.4 ) {
  //   cout << "test: " << eta << " " << eff << " " << linearDensity << " " << pileupDegradationFactor << " " << pileupDegradationFactorWithoutTiming << " " << pileupDegradationFactorWithTiming << "\n";
  //   cout << "options: " << WithBarrelTiming << " " << WithEndcapTiming << " " << WithHGCTiming << "\n";
  // }
  
  eff = eff * pileupDegradationFactor;
  return eff;
}

double getMuonEff( double pt, double eta, TH2F* muonEffHist, TH1F *linearDensityHist, bool WithBarrelTiming = false, bool WithEndcapTiming = false, bool WithHGCTiming = false) {

  double eff = 0.0;
  double countBins = 0;

  for (int i=1; i<linearDensityHist->GetXaxis()->GetNbins()+1; i++) {
    double density = linearDensityHist->GetXaxis()->GetBinCenter(i);
    eff += linearDensityHist->GetBinContent(i) * computeAvergeMuonEff( pt, eta, muonEffHist, density, WithBarrelTiming, WithEndcapTiming, WithHGCTiming);
    //cout << i << " " << density << " : " << linearDensityHist->GetBinContent(i) << " " << computeAvergeMuonEff( pt, eta, muonEffHist, density, WithBarrelTiming, WithEndcapTiming) << "\n";
  }
  //cout << "int eff: " << eff << "\n";

  return eff;
}


double getPhoEffFactor( double pt, double eta, TH1F *linearDensityHist, bool WithBarrelTiming = false, bool WithEndcapTiming = false, bool WithHGCTiming = false) {

  double sf = 0.0;
  double countBins = 0;

  for (int i=1; i<linearDensityHist->GetXaxis()->GetNbins()+1; i++) {
    double density = linearDensityHist->GetXaxis()->GetBinCenter(i);

    double tmpsf = 0;
    // double pileupDegradationFactorWithoutTiming = (0.95 - (density - 0.3)*(0.95-0.79)/(1.5-0.3)) / 0.95;
    // double pileupDegradationFactorWithTiming = (0.95 - (density - 0.3)*(0.95-0.91)/(1.5-0.3)) / 0.95;
  
    //Use a flat line for new estimates which is averaged for the 200PU scenario 
    double pileupDegradationFactorWithoutTiming = 0.90;
    double pileupDegradationFactorWithTiming = 0.95; 
    double pileupDegradationFactor = pileupDegradationFactorWithoutTiming;

    if (fabs(eta) < 1.5) {
      if (WithBarrelTiming) pileupDegradationFactor = 0.95;
      else pileupDegradationFactor = 0.90;    
    } else if (fabs(eta) >= 1.5 && fabs(eta) < 2.4) {
      if (WithEndcapTiming) pileupDegradationFactor = 0.97;  
      else if (WithHGCTiming) pileupDegradationFactor = 0.95;  
      else pileupDegradationFactor = 0.90;    
    }
  
    tmpsf = pileupDegradationFactor;

    sf += linearDensityHist->GetBinContent(i) * tmpsf;
    //cout << i << " " << density << " : " << linearDensityHist->GetBinContent(i) << " " << computeAvergeMuonEff( pt, eta, muonEffHist, density, WithBarrelTiming, WithEndcapTiming) << "\n";
  }
  //cout << "int sf: " << sf << "\n";

  return sf;
}

double getBTagEffFactor( double pt, double eta, TH1F *linearDensityHist, bool WithBarrelTiming = false, bool WithEndcapTiming = false, bool WithHGCTiming = false) {

  double sf = 0.0;
  double countBins = 0;

  for (int i=1; i<linearDensityHist->GetXaxis()->GetNbins()+1; i++) {
    double density = linearDensityHist->GetXaxis()->GetBinCenter(i);

    double tmpsf = 0;
    double pileupDegradationFactor = 1;

    if (fabs(eta) < 1.5) {
      double pileupDegradationFactorWithoutTiming = 1.0;
      double pileupDegradationFactorWithTiming = 1.0;
      // if (density < 0.85) {
      // 	pileupDegradationFactorWithoutTiming = (0.78 - (density - 0.3)*(0.78-0.76)/(0.85-0.3)) / 0.78;
      // } else {
      // 	pileupDegradationFactorWithoutTiming = (0.76 - (density - 0.85)*(0.76-0.75)/(1.34-0.85)) / 0.78;
      // }
      //Use a flat line for new estimates which is averaged for the 200PU scenario
      pileupDegradationFactorWithoutTiming = 0.76 / 0.79;

      pileupDegradationFactor = 0.76/0.79;
      if (WithBarrelTiming && fabs(eta) < 1.5) pileupDegradationFactor = 0.79/0.79;    

    } else if (fabs(eta) < 2.5) {
      double pileupDegradationFactorWithoutTiming = 1.0;
      double pileupDegradationFactorWithTiming = 1.0;
      // if (density < 0.85) {
      // 	pileupDegradationFactorWithoutTiming = (0.62 - (density - 0.3)*(0.62-0.59)/(0.85-0.3)) / 0.62;
      // } else {
      // 	pileupDegradationFactorWithoutTiming = (0.59 - (density - 0.85)*(0.59-0.55)/(1.34-0.85)) / 0.62;
      // }
      //Use a flat line for new estimates which is averaged for the 200PU scenario
      pileupDegradationFactorWithoutTiming = 0.66 / 0.70;

      pileupDegradationFactor = 0.66/0.70;
      if (!WithEndcapTiming && WithHGCTiming ) pileupDegradationFactor = 0.68/0.70;
      if (WithEndcapTiming ) pileupDegradationFactor = 0.70/0.70;
    } else if (fabs(eta) < 3.0) {
      double pileupDegradationFactorWithoutTiming = 1.0;
      double pileupDegradationFactorWithTiming = 1.0;
      // if (density < 0.85) {
      // 	pileupDegradationFactorWithoutTiming = (0.62 - (density - 0.3)*(0.62-0.59)/(0.85-0.3)) / 0.62;
      // } else {
      // 	pileupDegradationFactorWithoutTiming = (0.59 - (density - 0.85)*(0.59-0.55)/(1.34-0.85)) / 0.62;
      // }
      //Use a flat line for new estimates which is averaged for the 200PU scenario
      pileupDegradationFactor = 0.53 / 0.55;
      if (!WithEndcapTiming && WithHGCTiming) pileupDegradationFactor = 0.54/0.55;
      if (WithEndcapTiming ) pileupDegradationFactor = 0.55 / 0.55;
    }

    tmpsf = pileupDegradationFactor;

    sf += linearDensityHist->GetBinContent(i) * tmpsf;
    //cout << i << " " << density << " : " << linearDensityHist->GetBinContent(i) << " " <<  tmpsf << "\n";
  }
  //cout << "int sf: " << sf << "\n";

  return sf;
}





void PlotData( TH1F* hist , string dataLabel, string varName, string label, string latexlabel, bool setLogy = false ) {

  TCanvas *cv =0;
  TLegend *legend = 0;

  cv = new TCanvas("cv","cv", 800,700);
  cv->SetHighLightColor(2);
  cv->SetFillColor(0);
  cv->SetBorderMode(0);
  cv->SetBorderSize(2);
  // cv->SetLeftMargin(0.16);
  // cv->SetRightMargin(0.3);
  // cv->SetTopMargin(0.07);
  // cv->SetBottomMargin(0.12);
  // cv->SetFrameBorderMode(0);  

  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(hist,(dataLabel).c_str(), "LP");

  // hist->SetFillColor(kBLack);
  // hist->SetFillStyle(1001);
      
  hist->SetLineWidth(2);
  hist->SetLineColor(kBlack);
  hist->Draw("e1same");

  legend->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Data #sqrt{s} = 8 TeV, #int L = %s fb^{-1}, %s","19.7", latexlabel.c_str()));
  //tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 8 TeV, %s", latexlabel.c_str()));
  tex->Draw();
  

  if(setLogy) {
    cv->SetLogy(true);
    cv->SaveAs(Form("HggRazor_%s%s_Logy.gif",varName.c_str(),label.c_str()));
  } else {
    cv->SetLogy(false);
    cv->SaveAs(Form("HggRazor_%s%s.gif",varName.c_str(), label.c_str()));
  }

}



//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void MakeHZZPlots ( string datafile, string dataLabel,  vector<string> bkgfiles,vector<string> bkgLabels, vector<int> bkgColors, int boxOption = 0, int option = -1, string label = "", string latexlabel = "") {

  string Label  = "";
  if(label != "") Label = "_"+label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TFile *muonEffFile = TFile::Open("/eos/cms/store/group/phys_susy/razor/MTDUpgradePerformance/Efficiencies/Efficiency_PromptMuon_TTJets_25ns_Loose_Fullsim.root","READ");
  TH2F *muonEffHist = (TH2F*)muonEffFile->Get("Efficiency_PtEta");

  TFile *pileupDensityFile = TFile::Open("/eos/cms/store/group/phys_susy/razor/MTDUpgradePerformance/lineDensity.root","READ");
  TH1F *pileupDensityHist140 = (TH1F*)pileupDensityFile->Get("hL140");
  TH1F *pileupDensityHist200 = (TH1F*)pileupDensityFile->Get("hL200");

  vector<string> inputfiles;
  vector<string> processLabels;
  vector<int> color;

  bool hasData = false;
  if (datafile != "") {
    hasData = true;
    inputfiles.push_back(datafile);
    processLabels.push_back(dataLabel);
    color.push_back(kBlack);
  } else {
    hasData = false;
    // inputfiles.push_back("");
    // processLabels.push_back("");    
    // color.push_back(kBlack);
  }
  assert(bkgfiles.size() == bkgLabels.size());
  for (int i=0; i < bkgfiles.size(); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processLabels.push_back(bkgLabels[i]);
     color.push_back(bkgColors[i]);
  }


  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  vector<double> EventCount;
  vector<double> EventCountErrSqr;
  vector<TH1F*> histM4l;
  vector<TH1F*> histMR;
  vector<TH1F*> histRsq;

  TH1F *hPt_140PU = makeNewHist( "hPt_140PU", ";p_{T 4l} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtBarrelEndcapTiming_140PU = makeNewHist( "hPtBarrelEndcapTiming_140PU", ";p_{T 4l} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtBarrelTiming_140PU = makeNewHist( "hPtBarrelTiming_140PU", ";p_{T 4l} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtEndcapTiming_140PU = makeNewHist( "hPtEndcapTiming_140PU", ";p_{T 4l} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPt_200PU = makeNewHist( "hPt_200PU", ";p_{T 4l} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtBarrelEndcapTiming_200PU = makeNewHist( "hPtBarrelEndcapTiming_200PU", ";p_{T 4l} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtBarrelTiming_200PU = makeNewHist( "hPtBarrelTiming_200PU", ";p_{T 4l} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtEndcapTiming_200PU = makeNewHist( "hPtEndcapTiming_200PU", ";p_{T 4l} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtHGCTiming_200PU = makeNewHist( "hPtHGCTiming_200PU", ";p_{T 4l} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );


  TH1F *hRapidity_140PU = makeNewHist( "hRapidity_140PU", ";y_{Higgs};Number of Events", 40, -3, 3, kBlue, false );
  TH1F *hRapidityBarrelEndcapTiming_140PU = makeNewHist( "hRapidityBarrelEndcapTiming_140PU", ";y_{Higgs};Number of Events", 40, -3, 3, kBlue, false );
  TH1F *hRapidityBarrelTiming_140PU = makeNewHist( "hRapidityBarrelTiming_140PU", ";y_{Higgs};Number of Events", 40, -3, 3, kBlue, false );
  TH1F *hRapidityEndcapTiming_140PU = makeNewHist( "hRapidityEndcapTiming_140PU", ";y_{Higgs};Number of Events", 40, -3, 3, kBlue, false );
  TH1F *hRapidity_200PU = makeNewHist( "hRapidity_200PU", ";y_{Higgs};Number of Events", 40, -3, 3, kBlue, false );
  TH1F *hRapidityBarrelEndcapTiming_200PU = makeNewHist( "hRapidityBarrelEndcapTiming_200PU", ";y_{Higgs};Number of Events", 40, -3, 3, kBlue, false );
  TH1F *hRapidityBarrelTiming_200PU = makeNewHist( "hRapidityBarrelTiming_200PU", ";y_{Higgs};Number of Events", 40, -3, 3, kBlue, false );
  TH1F *hRapidityEndcapTiming_200PU = makeNewHist( "hRapidityEndcapTiming_200PU", ";y_{Higgs};Number of Events", 40, -3, 3, kBlue, false );
  TH1F *hRapidityHGCTiming_200PU = makeNewHist( "hRapidityHGCTiming_200PU", ";y_{Higgs};Number of Events", 40, -3, 3, kBlue, false );


  assert (inputfiles.size() == processLabels.size());
  for (int i=0; i < inputfiles.size(); ++i) {        
    EventCount.push_back(0);
    EventCountErrSqr.push_back(0);
    histM4l.push_back( makeNewHist( Form("M4l_%s",processLabels[i].c_str()), ";M_{4l} [GeV/c^{2}];Number of Events", 20, 100, 160, color[i],  (hasData && i==0) ));  
    histMR.push_back( makeNewHist( Form("MR_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}];Number of Events", 25, 0, 500, color[i],  (hasData && i==0) ));  
    histRsq.push_back( makeNewHist( Form("Rsq_%s",processLabels[i].c_str()), ";R^{2} ;Number of Events", 50, 0, 0.4, color[i],  (hasData && i==0) ));  
  }

  //*******************************************************************************************
  //Define Counts
  //*******************************************************************************************
  double TotalCounts_140PU = 0;
  double TotalCounts_140PU_BarrelTiming = 0;
  double TotalCounts_140PU_EndcapTiming = 0;
  double TotalCounts_140PU_BarrelEndcapTiming = 0;
  double TotalCounts_200PU = 0;
  double TotalCounts_200PU_BarrelTiming = 0;
  double TotalCounts_200PU_EndcapTiming = 0;
  double TotalCounts_200PU_BarrelEndcapTiming = 0;
  double TotalCounts_200PU_HGCTiming = 0;


  //*******************************************************************************************
  //Read files
  //*******************************************************************************************
  for (uint i=0; i < inputfiles.size(); ++i) {

    TFile* inputFile = new TFile(inputfiles[i].c_str(),"READ");
    assert(inputFile);
    TTree* tree = 0;
    tree = (TTree*)inputFile->Get("HZZRazor");    
    assert(tree);

    int genlep1Id = 0;
    int genlep2Id = 0;
    int genlep3Id = 0;
    int genlep4Id = 0;
    float genlep1pt = 0;
    float genlep1eta = 0;
    float genlep1phi = 0;
    float genlep2pt = 0;
    float genlep2eta = 0;
    float genlep2phi = 0;
    float genlep3pt = 0;
    float genlep3eta = 0;
    float genlep3phi = 0;
    float genlep4pt = 0;
    float genlep4eta = 0;
    float genlep4phi = 0;
  

    tree->SetBranchAddress("genlep1Id",&genlep1Id);
    tree->SetBranchAddress("genlep1Pt",&genlep1pt);
    tree->SetBranchAddress("genlep1Eta",&genlep1eta);
    tree->SetBranchAddress("genlep1Phi",&genlep1phi);
    tree->SetBranchAddress("genlep2Id",&genlep2Id);
    tree->SetBranchAddress("genlep2Pt",&genlep2pt);
    tree->SetBranchAddress("genlep2Eta",&genlep2eta);
    tree->SetBranchAddress("genlep2Phi",&genlep2phi);
    tree->SetBranchAddress("genlep3Id",&genlep3Id);
    tree->SetBranchAddress("genlep3Pt",&genlep3pt);
    tree->SetBranchAddress("genlep3Eta",&genlep3eta);
    tree->SetBranchAddress("genlep3Phi",&genlep3phi);
    tree->SetBranchAddress("genlep4Id",&genlep4Id);
    tree->SetBranchAddress("genlep4Pt",&genlep4pt);
    tree->SetBranchAddress("genlep4Eta",&genlep4eta);
    tree->SetBranchAddress("genlep4Phi",&genlep4phi);

    cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";

    for (int n=0;n<tree->GetEntries();n++) { 
    
      tree->GetEntry(n);
      if (n % 100000 == 0) cout << "Processing Event " << n << "\n";       
      
      double weight_140PU = 1.0;
      double weightBarrelTiming_140PU = 1.0;
      double weightEndcapTiming_140PU = 1.0;      
      double weightBarrelEndcapTiming_140PU = 1.0;
      double weight_200PU = 1.0;
      double weightBarrelTiming_200PU = 1.0;
      double weightEndcapTiming_200PU = 1.0;
      double weightBarrelEndcapTiming_200PU = 1.0;
      double weightHGCTiming_200PU = 1.0;

      if (!(abs(genlep1Id) == 13 && abs(genlep2Id) == 13 && abs(genlep3Id) == 13 && abs(genlep4Id) == 13)) continue;
      if (!(genlep1pt > 15 && genlep2pt > 15 && genlep3pt > 5 && genlep4pt > 5)) continue;
      if (!(abs(genlep1eta) < 2.4 && abs(genlep2eta) < 2.4 && abs(genlep3eta) < 2.4 && abs(genlep4eta) < 2.4 )) continue;

      // weight_140PU = weight_140PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist140, false, false);
      // weight_140PU = weight_140PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist140, false, false);
      // weight_140PU = weight_140PU * getMuonEff( genlep3pt, genlep3eta, muonEffHist, pileupDensityHist140, false, false);
      // weight_140PU = weight_140PU * getMuonEff( genlep4pt, genlep4eta, muonEffHist, pileupDensityHist140, false, false);
      // weightBarrelTiming_140PU = weightBarrelTiming_140PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist140, true, false);
      // weightBarrelTiming_140PU = weightBarrelTiming_140PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist140, true, false);
      // weightBarrelTiming_140PU = weightBarrelTiming_140PU * getMuonEff( genlep3pt, genlep3eta, muonEffHist, pileupDensityHist140, true, false);
      // weightBarrelTiming_140PU = weightBarrelTiming_140PU * getMuonEff( genlep4pt, genlep4eta, muonEffHist, pileupDensityHist140, true, false);
      // weightEndcapTiming_140PU = weightEndcapTiming_140PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist140, false, true);
      // weightEndcapTiming_140PU = weightEndcapTiming_140PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist140, false, true);
      // weightEndcapTiming_140PU = weightEndcapTiming_140PU * getMuonEff( genlep3pt, genlep3eta, muonEffHist, pileupDensityHist140, false, true);
      // weightEndcapTiming_140PU = weightEndcapTiming_140PU * getMuonEff( genlep4pt, genlep4eta, muonEffHist, pileupDensityHist140, false, true);
      // weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist140, true, true);
      // weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist140, true, true);
      // weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getMuonEff( genlep3pt, genlep3eta, muonEffHist, pileupDensityHist140, true, true);
      // weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getMuonEff( genlep4pt, genlep4eta, muonEffHist, pileupDensityHist140, true, true);

      TotalCounts_140PU += weight_140PU;
      TotalCounts_140PU_BarrelTiming += weightBarrelTiming_140PU;
      TotalCounts_140PU_EndcapTiming += weightEndcapTiming_140PU;
      TotalCounts_140PU_BarrelEndcapTiming += weightBarrelEndcapTiming_140PU;
      //cout << weight_140PU << "\n";

      // if (weight_140PU > 1) {
      // 	cout << weight_140PU << "\n";
      // 	cout << genlep1pt << " " << genlep2pt << " " << genlep3pt << " " << genlep4pt <<  "\n";
      // 	cout << genlep1eta << " " << genlep2eta << " " << genlep3eta << " " << genlep4eta << "\n";
      // 	cout << weight_140PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist140, false, false) << " "
      // 	     << weight_140PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist140, false, false) << " "
      // 	     << weight_140PU * getMuonEff( genlep3pt, genlep3eta, muonEffHist, pileupDensityHist140, false, false) << " "
      // 	     << weight_140PU * getMuonEff( genlep4pt, genlep4eta, muonEffHist, pileupDensityHist140, false, false) << " "
      // 	     << "\n";
      // }


      weight_200PU = weight_200PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist200, false, false);
      weight_200PU = weight_200PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist200, false, false);
      weight_200PU = weight_200PU * getMuonEff( genlep3pt, genlep3eta, muonEffHist, pileupDensityHist200, false, false);
      weight_200PU = weight_200PU * getMuonEff( genlep4pt, genlep4eta, muonEffHist, pileupDensityHist200, false, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist200, true, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist200, true, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getMuonEff( genlep3pt, genlep3eta, muonEffHist, pileupDensityHist200, true, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getMuonEff( genlep4pt, genlep4eta, muonEffHist, pileupDensityHist200, true, false);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist200, false, true);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist200, false, true);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getMuonEff( genlep3pt, genlep3eta, muonEffHist, pileupDensityHist200, false, true);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getMuonEff( genlep4pt, genlep4eta, muonEffHist, pileupDensityHist200, false, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist200, true, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist200, true, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getMuonEff( genlep3pt, genlep3eta, muonEffHist, pileupDensityHist200, true, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getMuonEff( genlep4pt, genlep4eta, muonEffHist, pileupDensityHist200, true, true);
      weightHGCTiming_200PU = weightHGCTiming_200PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist200, false, false, true);
      weightHGCTiming_200PU = weightHGCTiming_200PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist200, false, false, true);
      weightHGCTiming_200PU = weightHGCTiming_200PU * getMuonEff( genlep3pt, genlep3eta, muonEffHist, pileupDensityHist200, false, false, true);
      weightHGCTiming_200PU = weightHGCTiming_200PU * getMuonEff( genlep4pt, genlep4eta, muonEffHist, pileupDensityHist200, false, false, true);


      TotalCounts_200PU += weight_200PU;
      TotalCounts_200PU_BarrelTiming += weightBarrelTiming_200PU;
      TotalCounts_200PU_EndcapTiming += weightEndcapTiming_200PU;
      TotalCounts_200PU_BarrelEndcapTiming += weightBarrelEndcapTiming_200PU;
      TotalCounts_200PU_HGCTiming += weightHGCTiming_200PU;

      // cout << genlep1pt << " " << genlep2pt << " " << genlep3pt << " " << genlep4pt << " " << weight << "\n";
      //cout << genlep1eta << " " << genlep2eta << " " << genlep3eta << " " << genlep4eta << " " << weight << "\n";

      TLorentzVector v1; v1.SetPtEtaPhiM(genlep1pt,genlep1eta,genlep1phi,0.1057);
      TLorentzVector v2; v2.SetPtEtaPhiM(genlep2pt,genlep2eta,genlep2phi,0.1057);
      TLorentzVector v3; v3.SetPtEtaPhiM(genlep3pt,genlep3eta,genlep3phi,0.1057);
      TLorentzVector v4; v4.SetPtEtaPhiM(genlep4pt,genlep4eta,genlep4phi,0.1057);
      
      double pt = (v1+v2+v3+v4).Pt();
      double y = (v1+v2+v3+v4).Rapidity();

      //cout << m << " \n";

      hPt_140PU->Fill(pt, weight_140PU); 
      hPtBarrelTiming_140PU->Fill(pt, weightBarrelTiming_140PU); 
      hPtEndcapTiming_140PU->Fill(pt, weightEndcapTiming_140PU); 
      hPtBarrelEndcapTiming_140PU->Fill(pt, weightBarrelEndcapTiming_140PU); 
      hPt_200PU->Fill(pt, weight_200PU); 
      hPtBarrelTiming_200PU->Fill(pt, weightBarrelTiming_200PU); 
      hPtEndcapTiming_200PU->Fill(pt, weightEndcapTiming_200PU); 
      hPtBarrelEndcapTiming_200PU->Fill(pt, weightBarrelEndcapTiming_200PU); 
      hPtHGCTiming_200PU->Fill(pt, weightHGCTiming_200PU); 

      hRapidity_140PU->Fill(y, weight_140PU); 
      hRapidityBarrelTiming_140PU->Fill(y, weightBarrelTiming_140PU); 
      hRapidityEndcapTiming_140PU->Fill(y, weightEndcapTiming_140PU); 
      hRapidityBarrelEndcapTiming_140PU->Fill(y, weightBarrelEndcapTiming_140PU); 
      hRapidity_200PU->Fill(y, weight_200PU); 
      hRapidityBarrelTiming_200PU->Fill(y, weightBarrelTiming_200PU); 
      hRapidityEndcapTiming_200PU->Fill(y, weightEndcapTiming_200PU); 
      hRapidityBarrelEndcapTiming_200PU->Fill(y, weightBarrelEndcapTiming_200PU); 
      hRapidityHGCTiming_200PU->Fill(y, weightHGCTiming_200PU); 
    }

    inputFile->Close();
    delete inputFile;
  
  }
  
  cout << "Total: " << "\n";
  cout << TotalCounts_140PU << " : " << TotalCounts_140PU_BarrelTiming << " " << TotalCounts_140PU_EndcapTiming << " " << TotalCounts_140PU_BarrelEndcapTiming<< "\n";
  cout << TotalCounts_200PU << " : " << TotalCounts_200PU_HGCTiming << " " << TotalCounts_200PU_BarrelTiming << " " << TotalCounts_200PU_EndcapTiming << " " << TotalCounts_200PU_BarrelEndcapTiming<< "\n";

  //*******************************************************************************************
  //Normalize Hists
  //*******************************************************************************************
  hPt_140PU->Scale(1/TotalCounts_140PU);
  hPtBarrelTiming_140PU->Scale(1/TotalCounts_140PU);
  hPtEndcapTiming_140PU->Scale(1/TotalCounts_140PU);
  hPtBarrelEndcapTiming_140PU->Scale(1/TotalCounts_140PU);

  hPt_200PU->Scale(1/TotalCounts_200PU);
  hPtBarrelTiming_200PU->Scale(1/TotalCounts_200PU);
  hPtEndcapTiming_200PU->Scale(1/TotalCounts_200PU);
  hPtBarrelEndcapTiming_200PU->Scale(1/TotalCounts_200PU);
  hPtHGCTiming_200PU->Scale(1/TotalCounts_200PU);

  hRapidity_140PU->Scale(1/TotalCounts_140PU);
  hRapidityBarrelTiming_140PU->Scale(1/TotalCounts_140PU);
  hRapidityEndcapTiming_140PU->Scale(1/TotalCounts_140PU);
  hRapidityBarrelEndcapTiming_140PU->Scale(1/TotalCounts_140PU);

  hRapidity_200PU->Scale(1/TotalCounts_200PU);
  hRapidityBarrelTiming_200PU->Scale(1/TotalCounts_200PU);
  hRapidityEndcapTiming_200PU->Scale(1/TotalCounts_200PU);
  hRapidityBarrelEndcapTiming_200PU->Scale(1/TotalCounts_200PU);
  hRapidityHGCTiming_200PU->Scale(1/TotalCounts_200PU);

  //*******************************************************************************************
  //Draw Plots
  //*******************************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;
  TLatex *tex = 0;



  //*******************************************************************************************
  //Higgs pT
  //*******************************************************************************************
 
  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend = new TLegend(0.45,0.70,0.88,0.88);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  legend->AddEntry(hPt_140PU, "No Timing");
  legend->AddEntry(hPtBarrelTiming_140PU, "Barrel Timing Only");
  legend->AddEntry(hPtEndcapTiming_140PU, "Endcap Timing Only");
  legend->AddEntry(hPtBarrelEndcapTiming_140PU, "Barrel+Endcap Timing");

 
  hPtBarrelEndcapTiming_140PU->Draw("hist");
  hPtBarrelTiming_140PU->Draw("samehist");
  hPtEndcapTiming_140PU->Draw("samehist");
  hPt_140PU->Draw("samehist");

  hPtBarrelEndcapTiming_140PU->GetYaxis()->SetTitle("Fraction of Events");
  hPtBarrelEndcapTiming_140PU->GetYaxis()->SetRangeUser(0,0.4);
  hPtBarrelEndcapTiming_140PU->GetYaxis()->SetTitleOffset(1.8);
  hPtBarrelEndcapTiming_140PU->GetXaxis()->SetTitleOffset(1.3);
  hPt_140PU->SetLineColor(kBlack);
  hPtBarrelTiming_140PU->SetLineColor(kRed);
  hPtEndcapTiming_140PU->SetLineColor(kGreen+2);
  hPtBarrelEndcapTiming_140PU->SetLineColor(kBlue);

  cv->SetLogx();
  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.5, 0.5, "test");
  tex->DrawLatex(0.45, 0.70, "Increase in Higgs#rightarrow ZZ#rightarrow 4l Yield");
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.45, 0.65, (string(Form("Barrel Timing Only : %.0f", 100*(TotalCounts_140PU_BarrelTiming/TotalCounts_140PU - 1))) + "%").c_str());
  tex->SetTextColor(kGreen+2);
  tex->DrawLatex(0.45, 0.60, (string(Form("Endcap Timing Only : %.0f", 100*(TotalCounts_140PU_EndcapTiming/TotalCounts_140PU - 1))) + "%").c_str());
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.45, 0.55, (string(Form("Barrel+Endcap Timing : %.0f", 100*(TotalCounts_140PU_BarrelEndcapTiming/TotalCounts_140PU - 1)))+"%").c_str());

  tex->SetTextSize(0.040);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.15, 0.94, "Higgs#rightarrow ZZ#rightarrow 4l ( 140 Pileup Distribution )");

  //tex->DrawLatex(0.20, 0.20, "Normalized to");
  //tex->DrawLatex(0.20, 0.15, "\"No Timing\" distribution");

  //tex->Draw();

  cv->SaveAs("HZZPt_TimingStudy_140PU.gif");
  cv->SaveAs("HZZPt_TimingStudy_140PU.pdf");






  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend = new TLegend(0.45,0.70,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  legend->AddEntry(hPt_200PU, "No Timing");
  legend->AddEntry(hPtHGCTiming_200PU, "Timing (p_{T} > 2 GeV) Only");
  legend->AddEntry(hPtBarrelTiming_200PU, "With Barrel Timing");
  legend->AddEntry(hPtEndcapTiming_200PU, "With Endcap Timing");
  legend->AddEntry(hPtBarrelEndcapTiming_200PU, "With Barrel+Endcap Timing");

 
  hPtBarrelEndcapTiming_200PU->Draw("hist");
  hPtHGCTiming_200PU->Draw("samehist");
  hPtBarrelTiming_200PU->Draw("samehist");
  hPtEndcapTiming_200PU->Draw("samehist");
  hPt_200PU->Draw("samehist");

  hPtBarrelEndcapTiming_200PU->GetYaxis()->SetTitle("Fraction of Events");
  hPtBarrelEndcapTiming_200PU->GetYaxis()->SetRangeUser(0,0.4);
  hPtBarrelEndcapTiming_200PU->GetYaxis()->SetTitleOffset(1.8);
  hPtBarrelEndcapTiming_200PU->GetXaxis()->SetTitleOffset(1.3);
  hPt_200PU->SetLineColor(kBlack);
  hPtHGCTiming_200PU->SetLineColor(kViolet);
  hPtBarrelTiming_200PU->SetLineColor(kRed);
  hPtEndcapTiming_200PU->SetLineColor(kGreen+2);
  hPtBarrelEndcapTiming_200PU->SetLineColor(kBlue);

  cv->SetLogx();
  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.5, 0.5, "test");
  tex->DrawLatex(0.45, 0.65, "Increase in Higgs#rightarrow ZZ#rightarrow 4l Yield");
  tex->SetTextColor(kViolet);
  tex->DrawLatex(0.45, 0.60, (string(Form("Timing (p_{T} > 2 GeV) Only : %.0f", 100*(TotalCounts_200PU_HGCTiming/TotalCounts_200PU - 1)))+"%").c_str());
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.45, 0.55, (string(Form("Barrel Timing Only : %.0f", 100*(TotalCounts_200PU_BarrelTiming/TotalCounts_200PU - 1))) + "%").c_str());
  tex->SetTextColor(kGreen+2);
  tex->DrawLatex(0.45, 0.50, (string(Form("Endcap Timing Only : %.0f", 100*(TotalCounts_200PU_EndcapTiming/TotalCounts_200PU - 1))) + "%").c_str());
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.45, 0.45, (string(Form("Barrel+Endcap Timing : %.0f", 100*(TotalCounts_200PU_BarrelEndcapTiming/TotalCounts_200PU - 1)))+"%").c_str());

  tex->SetTextSize(0.040);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.15, 0.94, "Higgs#rightarrow ZZ#rightarrow 4l ( 200 Pileup Distribution )");

  //tex->DrawLatex(0.20, 0.20, "Normalized to");
  //tex->DrawLatex(0.20, 0.15, "\"No Timing\" distribution");

  //tex->Draw();

  cv->SaveAs("HZZPt_TimingStudy_200PU.gif");
  cv->SaveAs("HZZPt_TimingStudy_200PU.pdf");



  //*******************************************************************************************
  //Higgs Y
  //*******************************************************************************************
 
  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend = new TLegend(0.45,0.75,0.85,0.88);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  legend->AddEntry(hRapidity_140PU, "No Timing");
  legend->AddEntry(hRapidityBarrelTiming_140PU, "Barrel Timing Only");
  //legend->AddEntry(hRapidityEndcapTiming_140PU, "Endcap Timing Only");
  legend->AddEntry(hRapidityBarrelEndcapTiming_140PU, "Barrel+Endcap Timing");

  hRapidityBarrelEndcapTiming_140PU->Draw("hist");
  hRapidityBarrelTiming_140PU->Draw("samehist");
  //hRapidityEndcapTiming_140PU->Draw("samehist");
  hRapidity_140PU->Draw("samehist");

  hRapidityBarrelEndcapTiming_140PU->GetYaxis()->SetTitle("Fraction of Events");
  hRapidityBarrelEndcapTiming_140PU->GetYaxis()->SetRangeUser(0,0.4);
  hRapidityBarrelEndcapTiming_140PU->GetYaxis()->SetTitleOffset(1.8);
  hRapidityBarrelEndcapTiming_140PU->GetXaxis()->SetTitleOffset(1.3);
  hRapidity_140PU->SetLineColor(kBlack);
  hRapidityBarrelTiming_140PU->SetLineColor(kRed);
  hRapidityEndcapTiming_140PU->SetLineColor(kGreen+2);
  hRapidityBarrelEndcapTiming_140PU->SetLineColor(kBlue);

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.5, 0.5, "test");
  tex->DrawLatex(0.45, 0.70, "Increase in Higgs#rightarrow ZZ#rightarrow 4l Yield");
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.45, 0.65, (string(Form("Barrel Timing Only : %.0f", 100*(TotalCounts_140PU_BarrelTiming/TotalCounts_140PU - 1))) + "%").c_str());
  // tex->SetTextColor(kGreen+2);
  // tex->DrawLatex(0.45, 0.60, (string(Form("Endcap Timing Only : %.0f", 100*(TotalCounts_140PU_EndcapTiming/TotalCounts_140PU - 1))) + "%").c_str());
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.45, 0.60, (string(Form("Barrel+Endcap Timing : %.0f", 100*(TotalCounts_140PU_BarrelEndcapTiming/TotalCounts_140PU - 1)))+"%").c_str());

  tex->SetTextSize(0.040);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.15, 0.94, "Higgs#rightarrow ZZ#rightarrow 4l ( 140 Pileup Distribution )");

  tex->SetTextSize(0.030);
  //tex->DrawLatex(0.20, 0.03, "Normalized to \"No Timing\" distribution");

  //tex->Draw();

  cv->SaveAs("HZZRapidity_TimingStudy_140PU.gif");
  cv->SaveAs("HZZRapidity_TimingStudy_140PU.pdf");




  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend = new TLegend(0.45,0.70,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  legend->AddEntry(hRapidity_200PU, "No Timing");
  legend->AddEntry(hRapidityHGCTiming_200PU, "Timing (p_{T} > 2 GeV) Only");
  legend->AddEntry(hRapidityBarrelTiming_200PU, "With Barrel Timing");
  legend->AddEntry(hRapidityEndcapTiming_200PU, "With Endcap Timing");
  legend->AddEntry(hRapidityBarrelEndcapTiming_200PU, "With Barrel+Endcap Timing");

 
  hRapidityBarrelEndcapTiming_200PU->Draw("hist");
  hRapidityBarrelTiming_200PU->Draw("samehist");
  hRapidityEndcapTiming_200PU->Draw("samehist");
  hRapidityHGCTiming_200PU->Draw("samehist");
  hRapidity_200PU->Draw("samehist");

  hRapidityBarrelEndcapTiming_200PU->GetYaxis()->SetTitle("Fraction of Events");
  hRapidityBarrelEndcapTiming_200PU->GetYaxis()->SetRangeUser(0,0.4);
  hRapidityBarrelEndcapTiming_200PU->GetYaxis()->SetTitleOffset(1.8);
  hRapidityBarrelEndcapTiming_200PU->GetXaxis()->SetTitleOffset(1.3);
  hRapidity_200PU->SetLineColor(kBlack);
  hRapidityHGCTiming_200PU->SetLineColor(kViolet);
  hRapidityBarrelTiming_200PU->SetLineColor(kRed);
  hRapidityEndcapTiming_200PU->SetLineColor(kGreen+2);
  hRapidityBarrelEndcapTiming_200PU->SetLineColor(kBlue);

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.5, 0.5, "test");
  tex->DrawLatex(0.45, 0.65, "Increase in Higgs#rightarrow ZZ#rightarrow 4l Yield");
  tex->SetTextColor(kViolet);
  tex->DrawLatex(0.45, 0.60, (string(Form("Timing (p_{T} > 2 GeV) Only : %.0f", 100*(TotalCounts_200PU_HGCTiming/TotalCounts_200PU - 1))) + "%").c_str());
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.45, 0.55, (string(Form("Barrel Timing Only : %.0f", 100*(TotalCounts_200PU_BarrelTiming/TotalCounts_200PU - 1))) + "%").c_str());
  tex->SetTextColor(kGreen+2);
  tex->DrawLatex(0.45, 0.50, (string(Form("Endcap Timing Only : %.0f", 100*(TotalCounts_200PU_EndcapTiming/TotalCounts_200PU - 1))) + "%").c_str());
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.45, 0.45, (string(Form("Barrel+Endcap Timing : %.0f", 100*(TotalCounts_200PU_BarrelEndcapTiming/TotalCounts_200PU - 1)))+"%").c_str());

  tex->SetTextSize(0.040);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.15, 0.94, "Higgs#rightarrow ZZ#rightarrow 4l ( 200 Pileup Distribution )");

  tex->SetTextSize(0.030);
  //tex->DrawLatex(0.20, 0.03, "Normalized to \"No Timing\" distribution");

  //tex->Draw();

  cv->SaveAs("HZZRapidity_TimingStudy_200PU.gif");
  cv->SaveAs("HZZRapidity_TimingStudy_200PU.pdf");


 //*******************************************************************************************
  //Higgs Rapidity Ratio
  //*******************************************************************************************
 
  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  TH1F *hRapidityBarrelEndcapTimingRatio_200PU = (TH1F*)hRapidityBarrelEndcapTiming_200PU->Clone("hRapidityBarrelEndcapTimingRatio_200PU");
  TH1F *hRapidityBarrelTimingRatio_200PU = (TH1F*)hRapidityBarrelTiming_200PU->Clone("hRapidityBarrelTimingRatio_200PU");
  TH1F *hRapidityEndcapTimingRatio_200PU = (TH1F*)hRapidityEndcapTiming_200PU->Clone("hRapidityEndcapTimingRatio_200PU");
  TH1F *hRapidityHGCTimingRatio_200PU = (TH1F*)hRapidityHGCTiming_200PU->Clone("hRapidityHGCTimingRatio_200PU");

  for(int i=1; i<hRapidityBarrelEndcapTimingRatio_200PU->GetXaxis()->GetNbins(); i++) {    
    if (  hRapidity_200PU->GetBinContent(i) > 0 ) {
      hRapidityBarrelEndcapTimingRatio_200PU->SetBinContent( i, (hRapidityBarrelEndcapTiming_200PU->GetBinContent(i) / hRapidity_200PU->GetBinContent(i) - 1)*100);
      hRapidityBarrelTimingRatio_200PU->SetBinContent( i, (hRapidityBarrelTiming_200PU->GetBinContent(i) / hRapidity_200PU->GetBinContent(i) - 1)*100);
      hRapidityEndcapTimingRatio_200PU->SetBinContent( i, (hRapidityEndcapTiming_200PU->GetBinContent(i) / hRapidity_200PU->GetBinContent(i) - 1)*100);
      hRapidityHGCTimingRatio_200PU->SetBinContent( i, (hRapidityHGCTiming_200PU->GetBinContent(i) / hRapidity_200PU->GetBinContent(i) - 1)*100);
      //cout << i << " : " << hRapidityBarrelEndcapTiming_200PU->GetBinContent(i) << " " << hRapidityBarrelTiming_200PU->GetBinContent(i) << " , " << hRapidity_200PU->GetBinContent(i) << " : " << hRapidityBarrelEndcapTimingRatio_200PU->GetBinContent(i) << " : " << hRapidityBarrelTimingRatio_200PU->GetBinContent(i) << "\n";

      double IntegratedNoTiming = 0;
      double IntegratedHGCTiming = 0;
      double IntegratedBarrelTiming = 0;
      double IntegratedEndcapTiming = 0;
      double IntegratedBarrelEndcapTiming = 0;
      for(int k=i; k<hRapidityBarrelEndcapTimingRatio_200PU->GetXaxis()->GetNbins(); k++) { 
	IntegratedNoTiming += hRapidity_200PU->GetBinContent(i);
	IntegratedHGCTiming += hRapidityHGCTiming_200PU->GetBinContent(i);
	IntegratedBarrelTiming += hRapidityBarrelTiming_200PU->GetBinContent(i);
	IntegratedEndcapTiming += hRapidityEndcapTiming_200PU->GetBinContent(i);
	IntegratedBarrelEndcapTiming += hRapidityBarrelEndcapTiming_200PU->GetBinContent(i);
      }
      cout << i << " : " << hRapidity_200PU->GetXaxis()->GetBinCenter(i) << " : " 
	//<< IntegratedNoTiming <<  " | " 
	   << IntegratedHGCTiming/IntegratedNoTiming << " " 
	   << IntegratedBarrelTiming/IntegratedNoTiming << " "
	   << IntegratedEndcapTiming/IntegratedNoTiming << " "
	   << IntegratedBarrelEndcapTiming/IntegratedNoTiming << " "
	   << "\n";
    } else {
      hRapidityBarrelEndcapTimingRatio_200PU->SetBinContent( i, 0);
      hRapidityBarrelTimingRatio_200PU->SetBinContent( i, 0);
      hRapidityEndcapTimingRatio_200PU->SetBinContent( i, 0);
      hRapidityHGCTimingRatio_200PU->SetBinContent( i, 0);
    }
  }
  
  legend = new TLegend(0.45,0.70,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  legend->AddEntry(hRapidityHGCTiming_200PU, "Timing (p_{T} > 2 GeV) Only");
  legend->AddEntry(hRapidityBarrelTiming_200PU, "Barrel Timing Only");
  legend->AddEntry(hRapidityEndcapTiming_200PU, "Endcap Timing Only");
  legend->AddEntry(hRapidityBarrelEndcapTiming_200PU, "Barrel+Endcap Timing");

  hRapidityBarrelEndcapTimingRatio_200PU->Draw("hist");
  hRapidityBarrelTimingRatio_200PU->Draw("samehist");
  hRapidityEndcapTimingRatio_200PU->Draw("samehist");
  hRapidityHGCTimingRatio_200PU->Draw("samehist");

  hRapidityBarrelEndcapTimingRatio_200PU->GetYaxis()->SetTitle("Improvement in Efficiency over No-Timing Scenario (%)");
  hRapidityBarrelEndcapTimingRatio_200PU->GetYaxis()->SetRangeUser(0,60);
  hRapidityBarrelEndcapTimingRatio_200PU->GetYaxis()->SetTitleOffset(1.8);
  hRapidityBarrelEndcapTimingRatio_200PU->GetXaxis()->SetTitleOffset(1.3);
  hRapidityBarrelEndcapTimingRatio_200PU->GetXaxis()->SetRangeUser(0,2.3);
  hRapidityHGCTimingRatio_200PU->SetLineColor(kViolet);
  hRapidityBarrelTimingRatio_200PU->SetLineColor(kRed);
  hRapidityEndcapTimingRatio_200PU->SetLineColor(kGreen+2);
  hRapidityBarrelEndcapTimingRatio_200PU->SetLineColor(kBlue);
  hRapidityHGCTimingRatio_200PU->SetLineWidth(4);
  hRapidityBarrelTimingRatio_200PU->SetLineWidth(4);
  hRapidityEndcapTimingRatio_200PU->SetLineWidth(4);
  hRapidityBarrelEndcapTimingRatio_200PU->SetLineWidth(4);
 
  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
 
  tex->SetTextSize(0.040);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.15, 0.94, "Higgs#rightarrow ZZ#rightarrow 4l ( 200 Pileup Distribution )");


  cv->SaveAs("HZZRapidityRatio_TimingStudy_200PU.gif");
  cv->SaveAs("HZZRapidityRatio_TimingStudy_200PU.pdf");



  //*******************************************************************************************
  //Summarize Counts
  //*******************************************************************************************
  // for(int i=0; i<int(processLabels.size()); i++) {
  //   cout << processLabels[i] << " : " << EventCount[i] << " +/- " << sqrt(EventCountErrSqr[i]) << "\n";
  // }
  
   //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("TimingStudyPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();
  
  file->WriteTObject(hPt_200PU, Form("hPt%s",Label.c_str()), "WriteDelete");
 
 
 }



//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void MakeHMMPlots ( string datafile, string dataLabel,  vector<string> bkgfiles,vector<string> bkgLabels, vector<int> bkgColors, int boxOption = 0, int option = -1, string label = "", string latexlabel = "") {

  string Label  = "";
  if(label != "") Label = "_"+label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TFile *muonEffFile = TFile::Open("/eos/cms/store/group/phys_susy/razor/MTDUpgradePerformance/Efficiencies/Efficiency_PromptMuon_TTJets_25ns_Loose_Fullsim.root","READ");
  TH2F *muonEffHist = (TH2F*)muonEffFile->Get("Efficiency_PtEta");

  TFile *pileupDensityFile = TFile::Open("/eos/cms/store/group/phys_susy/razor/MTDUpgradePerformance/lineDensity.root","READ");
  TH1F *pileupDensityHist140 = (TH1F*)pileupDensityFile->Get("hL140");
  TH1F *pileupDensityHist200 = (TH1F*)pileupDensityFile->Get("hL200");

  vector<string> inputfiles;
  vector<string> processLabels;
  vector<int> color;

  bool hasData = false;
  if (datafile != "") {
    hasData = true;
    inputfiles.push_back(datafile);
    processLabels.push_back(dataLabel);
    color.push_back(kBlack);
  } else {
    hasData = false;
    // inputfiles.push_back("");
    // processLabels.push_back("");    
    // color.push_back(kBlack);
  }
  assert(bkgfiles.size() == bkgLabels.size());
  for (int i=0; i < bkgfiles.size(); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processLabels.push_back(bkgLabels[i]);
     color.push_back(bkgColors[i]);
  }


  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  vector<double> EventCount;
  vector<double> EventCountErrSqr;
  vector<TH1F*> histM4l;
  vector<TH1F*> histMR;
  vector<TH1F*> histRsq;

  TH1F *hPt_140PU = makeNewHist( "hPt_140PU", ";p_{T #mu#mu} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtBarrelEndcapTiming_140PU = makeNewHist( "hPtBarrelEndcapTiming_140PU", ";p_{T #mu#mu} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtBarrelTiming_140PU = makeNewHist( "hPtBarrelTiming_140PU", ";p_{T #mu#mu} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtEndcapTiming_140PU = makeNewHist( "hPtEndcapTiming_140PU", ";p_{T #mu#mu} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPt_200PU = makeNewHist( "hPt_200PU", ";p_{T #mu#mu} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtBarrelEndcapTiming_200PU = makeNewHist( "hPtBarrelEndcapTiming_200PU", ";p_{T #mu#mu} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtBarrelTiming_200PU = makeNewHist( "hPtBarrelTiming_200PU", ";p_{T #mu#mu} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );
  TH1F *hPtEndcapTiming_200PU = makeNewHist( "hPtEndcapTiming_200PU", ";p_{T #mu#mu} [GeV/c^{2}];Number of Events", 20, 0, 200, kBlue, false );

  assert (inputfiles.size() == processLabels.size());
  for (int i=0; i < inputfiles.size(); ++i) {        
    EventCount.push_back(0);
    EventCountErrSqr.push_back(0);
    histM4l.push_back( makeNewHist( Form("M4l_%s",processLabels[i].c_str()), ";M_{#mu#mu} [GeV/c^{2}];Number of Events", 20, 100, 160, color[i],  (hasData && i==0) ));  
    histMR.push_back( makeNewHist( Form("MR_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}];Number of Events", 25, 0, 500, color[i],  (hasData && i==0) ));  
    histRsq.push_back( makeNewHist( Form("Rsq_%s",processLabels[i].c_str()), ";R^{2} ;Number of Events", 50, 0, 0.4, color[i],  (hasData && i==0) ));  
  }

  //*******************************************************************************************
  //Define Counts
  //*******************************************************************************************
  double TotalCounts_140PU = 0;
  double TotalCounts_140PU_BarrelTiming = 0;
  double TotalCounts_140PU_EndcapTiming = 0;
  double TotalCounts_140PU_BarrelEndcapTiming = 0;
  double TotalCounts_200PU = 0;
  double TotalCounts_200PU_BarrelTiming = 0;
  double TotalCounts_200PU_EndcapTiming = 0;
  double TotalCounts_200PU_BarrelEndcapTiming = 0;


  //*******************************************************************************************
  //Read files
  //*******************************************************************************************
  for (uint i=0; i < inputfiles.size(); ++i) {

    TFile* inputFile = new TFile(inputfiles[i].c_str(),"READ");
    assert(inputFile);
    TTree* tree = 0;
    tree = (TTree*)inputFile->Get("HZZRazor");    
    assert(tree);

    int genlep1Id = 0;
    int genlep2Id = 0;
    int genlep3Id = 0;
    int genlep4Id = 0;
    float genlep1pt = 0;
    float genlep1eta = 0;
    float genlep1phi = 0;
    float genlep2pt = 0;
    float genlep2eta = 0;
    float genlep2phi = 0;
    float genlep3pt = 0;
    float genlep3eta = 0;
    float genlep3phi = 0;
    float genlep4pt = 0;
    float genlep4eta = 0;
    float genlep4phi = 0;
  

    tree->SetBranchAddress("genlep1Id",&genlep1Id);
    tree->SetBranchAddress("genlep1Pt",&genlep1pt);
    tree->SetBranchAddress("genlep1Eta",&genlep1eta);
    tree->SetBranchAddress("genlep1Phi",&genlep1phi);
    tree->SetBranchAddress("genlep2Id",&genlep2Id);
    tree->SetBranchAddress("genlep2Pt",&genlep2pt);
    tree->SetBranchAddress("genlep2Eta",&genlep2eta);
    tree->SetBranchAddress("genlep2Phi",&genlep2phi);
    tree->SetBranchAddress("genlep3Id",&genlep3Id);
    tree->SetBranchAddress("genlep3Pt",&genlep3pt);
    tree->SetBranchAddress("genlep3Eta",&genlep3eta);
    tree->SetBranchAddress("genlep3Phi",&genlep3phi);
    tree->SetBranchAddress("genlep4Id",&genlep4Id);
    tree->SetBranchAddress("genlep4Pt",&genlep4pt);
    tree->SetBranchAddress("genlep4Eta",&genlep4eta);
    tree->SetBranchAddress("genlep4Phi",&genlep4phi);

    cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";

    for (int n=0;n<tree->GetEntries();n++) { 
    
      tree->GetEntry(n);
      if (n % 100000 == 0) cout << "Processing Event " << n << "\n";       
      
      double weight_140PU = 1.0;
      double weightBarrelTiming_140PU = 1.0;
      double weightEndcapTiming_140PU = 1.0;
      double weightBarrelEndcapTiming_140PU = 1.0;
      double weight_200PU = 1.0;
      double weightBarrelTiming_200PU = 1.0;
      double weightEndcapTiming_200PU = 1.0;
      double weightBarrelEndcapTiming_200PU = 1.0;

      if (!(abs(genlep1Id) == 13 && abs(genlep2Id) == 13 )) continue;
      if (!(genlep1pt > 30 && genlep2pt > 30)) continue;
      if (!(abs(genlep1eta) < 2.4 && abs(genlep2eta) < 2.4)) continue;

      weight_140PU = weight_140PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist140, false, false);
      weight_140PU = weight_140PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist140, false, false);
      weightBarrelTiming_140PU = weightBarrelTiming_140PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist140, true, false);
      weightBarrelTiming_140PU = weightBarrelTiming_140PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist140, true, false);
      weightEndcapTiming_140PU = weightEndcapTiming_140PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist140, false, true);
      weightEndcapTiming_140PU = weightEndcapTiming_140PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist140, false, true);
      weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist140, true, true);
      weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist140, true, true);

      TotalCounts_140PU += weight_140PU;
      TotalCounts_140PU_BarrelTiming += weightBarrelTiming_140PU;
      TotalCounts_140PU_EndcapTiming += weightEndcapTiming_140PU;
      TotalCounts_140PU_BarrelEndcapTiming += weightBarrelEndcapTiming_140PU;
      //cout << weight_140PU << "\n";

      if (weight_140PU > 1) {
	cout << weight_140PU << "\n";
	cout << genlep1pt << " " << genlep2pt << " " << genlep3pt << " " << genlep4pt <<  "\n";
	cout << genlep1eta << " " << genlep2eta << " " << genlep3eta << " " << genlep4eta << "\n";
	cout << weight_140PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist140, false, false) << " "
	     << weight_140PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist140, false, false) << " "
	     << "\n";
      }


      weight_200PU = weight_200PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist200, false, false);
      weight_200PU = weight_200PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist200, false, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist200, true, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist200, true, false);
       weightEndcapTiming_200PU = weightEndcapTiming_200PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist200, false, true);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist200, false, true);
     weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getMuonEff( genlep1pt, genlep1eta, muonEffHist, pileupDensityHist200, true, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getMuonEff( genlep2pt, genlep2eta, muonEffHist, pileupDensityHist200, true, true);
 
      TotalCounts_200PU += weight_200PU;
      TotalCounts_200PU_BarrelTiming += weightBarrelTiming_200PU;
      TotalCounts_200PU_EndcapTiming += weightEndcapTiming_200PU;
      TotalCounts_200PU_BarrelEndcapTiming += weightBarrelEndcapTiming_200PU;

      // cout << genlep1pt << " " << genlep2pt << " " << genlep3pt << " " << genlep4pt << " " << weight << "\n";
      //cout << genlep1eta << " " << genlep2eta << " " << genlep3eta << " " << genlep4eta << " " << weight << "\n";

      TLorentzVector v1; v1.SetPtEtaPhiM(genlep1pt,genlep1eta,genlep1phi,0.1057);
      TLorentzVector v2; v2.SetPtEtaPhiM(genlep2pt,genlep2eta,genlep2phi,0.1057);
      
      double pt = (v1+v2).Pt();

      //cout << m << " \n";

      hPt_140PU->Fill(pt, weight_140PU); 
      hPtBarrelTiming_140PU->Fill(pt, weightBarrelTiming_140PU); 
      hPtEndcapTiming_140PU->Fill(pt, weightEndcapTiming_140PU); 
      hPtBarrelEndcapTiming_140PU->Fill(pt, weightBarrelEndcapTiming_140PU); 
      hPt_200PU->Fill(pt, weight_200PU); 
      hPtBarrelTiming_200PU->Fill(pt, weightBarrelTiming_200PU); 
      hPtEndcapTiming_200PU->Fill(pt, weightEndcapTiming_200PU); 
      hPtBarrelEndcapTiming_200PU->Fill(pt, weightBarrelEndcapTiming_200PU); 

    }

    inputFile->Close();
    delete inputFile;
  
  }
  
  cout << "Total:\n";
  cout << "140PU: " << TotalCounts_140PU << " : " << TotalCounts_140PU_BarrelTiming << " " << TotalCounts_140PU_EndcapTiming << " " << TotalCounts_140PU_BarrelEndcapTiming<< "\n";
  cout << "200PU: " << TotalCounts_200PU << " : " << TotalCounts_200PU_BarrelTiming << " " << TotalCounts_200PU_EndcapTiming << " " << TotalCounts_200PU_BarrelEndcapTiming<< "\n";

  //*******************************************************************************************
  //Normalize Hists
  //*******************************************************************************************
  hPt_140PU->Scale(1/TotalCounts_140PU);
  hPtBarrelTiming_140PU->Scale(1/TotalCounts_140PU);
  hPtEndcapTiming_140PU->Scale(1/TotalCounts_140PU);
  hPtBarrelEndcapTiming_140PU->Scale(1/TotalCounts_140PU);

  hPt_200PU->Scale(1/TotalCounts_200PU);
  hPtBarrelTiming_200PU->Scale(1/TotalCounts_200PU);
  hPtEndcapTiming_200PU->Scale(1/TotalCounts_200PU);
  hPtBarrelEndcapTiming_200PU->Scale(1/TotalCounts_200PU);

  //*******************************************************************************************
  //Draw Plots
  //*******************************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;
  TLatex *tex = 0;



  //*******************************************************************************************
  //M#mu#mu
  //*******************************************************************************************
 
  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend = new TLegend(0.45,0.75,0.85,0.88);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  legend->AddEntry(hPt_140PU, "No Timing");
  legend->AddEntry(hPtBarrelTiming_140PU, "Barrel Timing Only");
  //legend->AddEntry(hPtEndcapTiming_140PU, "Endcap Timing Only");
  legend->AddEntry(hPtBarrelEndcapTiming_140PU, "Barrel+Endcap Timing");

 
  hPtBarrelEndcapTiming_140PU->Draw("hist");
  hPtBarrelTiming_140PU->Draw("samehist");
  //hPtEndcapTiming_140PU->Draw("samehist");
  hPt_140PU->Draw("samehist");

  hPtBarrelEndcapTiming_140PU->GetYaxis()->SetTitle("Fraction of Events");
  hPtBarrelEndcapTiming_140PU->GetYaxis()->SetRangeUser(0,0.4);
  hPtBarrelEndcapTiming_140PU->GetYaxis()->SetTitleOffset(1.8);
  hPtBarrelEndcapTiming_140PU->GetXaxis()->SetTitleOffset(1.3);
  hPt_140PU->SetLineColor(kBlack);
  hPtBarrelTiming_140PU->SetLineColor(kRed);
  hPtEndcapTiming_140PU->SetLineColor(kGreen+2);
  hPtBarrelEndcapTiming_140PU->SetLineColor(kBlue);

  cv->SetLogx();
  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.5, 0.5, "test");
  tex->DrawLatex(0.45, 0.70, "Increase in Higgs#rightarrow #mu#mu Yield");
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.45, 0.65, (string(Form("Barrel Timing Only : %.0f", 100*(TotalCounts_140PU_BarrelTiming/TotalCounts_140PU - 1))) + "%").c_str());
  // tex->SetTextColor(kGreen+2);
  // tex->DrawLatex(0.45, 0.60, (string(Form("Endcap Timing Only : %.0f", 100*(TotalCounts_140PU_EndcapTiming/TotalCounts_140PU - 1))) + "%").c_str());
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.45, 0.60, (string(Form("Barrel+Endcap Timing : %.0f", 100*(TotalCounts_140PU_BarrelEndcapTiming/TotalCounts_140PU - 1)))+"%").c_str());

  tex->SetTextSize(0.040);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.15, 0.94, "Higgs#rightarrow #mu#mu ( 140 Pileup Distribution )");

  // tex->DrawLatex(0.20, 0.20, "Normalized to");
  // tex->DrawLatex(0.20, 0.15, "\"No Timing\" distribution");

  //tex->Draw();

  cv->SaveAs("HMMPt_TimingStudy_140PU.gif");
  cv->SaveAs("HMMPt_TimingStudy_140PU.pdf");




  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend = new TLegend(0.45,0.75,0.85,0.88);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  legend->AddEntry(hPt_200PU, "No Timing");
  legend->AddEntry(hPtBarrelTiming_200PU, "Barrel Timing Only");
  //legend->AddEntry(hPtEndcapTiming_200PU, "Endcap Timing Only");
  legend->AddEntry(hPtBarrelEndcapTiming_200PU, "Barrel+Endcap Timing");

 
  hPtBarrelEndcapTiming_200PU->Draw("hist");
  hPtBarrelTiming_200PU->Draw("samehist");
  //hPtEndcapTiming_200PU->Draw("samehist");
  hPt_200PU->Draw("samehist");

  hPtBarrelEndcapTiming_200PU->GetYaxis()->SetTitle("Fraction of Events");
  hPtBarrelEndcapTiming_200PU->GetYaxis()->SetRangeUser(0,0.4);
  hPtBarrelEndcapTiming_200PU->GetYaxis()->SetTitleOffset(1.8);
  hPtBarrelEndcapTiming_200PU->GetXaxis()->SetTitleOffset(1.3);
  hPt_200PU->SetLineColor(kBlack);
  hPtBarrelTiming_200PU->SetLineColor(kRed);
  hPtEndcapTiming_200PU->SetLineColor(kGreen+2);
  hPtBarrelEndcapTiming_200PU->SetLineColor(kBlue);

  cv->SetLogx();
  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.5, 0.5, "test");
  tex->DrawLatex(0.45, 0.70, "Increase in Higgs#rightarrow #mu#mu Yield");
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.45, 0.65, (string(Form("Barrel Timing Only : %.0f", 100*(TotalCounts_200PU_BarrelTiming/TotalCounts_200PU - 1))) + "%").c_str());
  // tex->SetTextColor(kGreen+2);
  // tex->DrawLatex(0.45, 0.60, (string(Form("Endcap Timing Only : %.0f", 100*(TotalCounts_200PU_EndcapTiming/TotalCounts_200PU - 1))) + "%").c_str());
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.45, 0.60, (string(Form("Barrel+Endcap Timing : %.0f", 100*(TotalCounts_200PU_BarrelEndcapTiming/TotalCounts_200PU - 1)))+"%").c_str());

  tex->SetTextSize(0.040);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.15, 0.94, "Higgs#rightarrow #mu#mu ( 200 Pileup Distribution )");

  // tex->DrawLatex(0.20, 0.20, "Normalized to");
  // tex->DrawLatex(0.20, 0.15, "\"No Timing\" distribution");

  //tex->Draw();

  cv->SaveAs("HMMPt_TimingStudy_200PU.gif");
  cv->SaveAs("HMMPt_TimingStudy_200PU.pdf");





  //*******************************************************************************************
  //Summarize Counts
  //*******************************************************************************************
  // for(int i=0; i<int(processLabels.size()); i++) {
  //   cout << processLabels[i] << " : " << EventCount[i] << " +/- " << sqrt(EventCountErrSqr[i]) << "\n";
  // }
  
   //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("TimingStudyPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();
  
  file->WriteTObject(hPt_200PU, Form("hPt%s",Label.c_str()), "WriteDelete");
 
 
 }


//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void MakeHHbbggPlots ( string datafile, string dataLabel,  vector<string> bkgfiles,vector<string> bkgLabels, vector<int> bkgColors, int boxOption = 0, int option = -1, string label = "", string latexlabel = "") {

  string Label  = "";
  if(label != "") Label = "_"+label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TFile *pileupDensityFile = TFile::Open("/eos/cms/store/group/phys_susy/razor/MTDUpgradePerformance/lineDensity.root","READ");
  TH1F *pileupDensityHist140 = (TH1F*)pileupDensityFile->Get("hL140");
  TH1F *pileupDensityHist200 = (TH1F*)pileupDensityFile->Get("hL200");

  vector<string> inputfiles;
  vector<string> processLabels;
  vector<int> color;

  bool hasData = false;
  if (datafile != "") {
    hasData = true;
    inputfiles.push_back(datafile);
    processLabels.push_back(dataLabel);
    color.push_back(kBlack);
  } else {
    hasData = false;
    // inputfiles.push_back("");
    // processLabels.push_back("");    
    // color.push_back(kBlack);
  }
  assert(bkgfiles.size() == bkgLabels.size());
  for (int i=0; i < bkgfiles.size(); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processLabels.push_back(bkgLabels[i]);
     color.push_back(bkgColors[i]);
  }


  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  vector<double> EventCount;
  vector<double> EventCountErrSqr;
  vector<TH1F*> histM4l;
  vector<TH1F*> histMR;
  vector<TH1F*> histRsq;

  TH1F *hybbgg_140PU = makeNewHist( "hybbgg_140PU", ";y_{HH} [GeV/c^{2}];Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggBarrelEndcapTiming_140PU = makeNewHist( "hybbggBarrelEndcapTiming_140PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggBarrelTiming_140PU = makeNewHist( "hybbggBarrelTiming_140PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggEndcapTiming_140PU = makeNewHist( "hybbggEndcapTiming_140PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbgg_200PU = makeNewHist( "hybbgg_200PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggBarrelEndcapTiming_200PU = makeNewHist( "hybbggBarrelEndcapTiming_200PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggBarrelTiming_200PU = makeNewHist( "hybbggBarrelTiming_200PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggEndcapTiming_200PU = makeNewHist( "hybbggEndcapTiming_200PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggHGCTiming_200PU = makeNewHist( "hybbggHGCTiming_200PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );


  //*******************************************************************************************
  //Define Counts
  //*******************************************************************************************
  double TotalCounts_140PU = 0;
  double TotalCounts_140PU_BarrelTiming = 0;
  double TotalCounts_140PU_EndcapTiming = 0;
  double TotalCounts_140PU_BarrelEndcapTiming = 0;
  double TotalCounts_200PU = 0;
  double TotalCounts_200PU_BarrelTiming = 0;
  double TotalCounts_200PU_EndcapTiming = 0;
  double TotalCounts_200PU_BarrelEndcapTiming = 0;
  double TotalCounts_200PU_HGCTiming = 0;


  //*******************************************************************************************
  //Read files
  //*******************************************************************************************
  for (uint i=0; i < inputfiles.size(); ++i) {

    TFile* inputFile = new TFile(inputfiles[i].c_str(),"READ");
    assert(inputFile);
    TTree* tree = 0;
    tree = (TTree*)inputFile->Get("HHbbgg");    
    assert(tree);

    float mbb = 0;
    float mgg = 0;
    float mbbgg = 0;
    float ybbgg = 0;
    float ptbbgg = 0;
    float pho1Pt = 0;
    float pho1Eta = 0;
    float pho2Pt = 0;
    float pho2Eta = 0;
    float bjet1Pt = 0;
    float bjet1Eta = 0;
    float bjet2Pt = 0;
    float bjet2Eta = 0;

    tree->SetBranchAddress("mbb",&mbb);
    tree->SetBranchAddress("mgg",&mgg);
    tree->SetBranchAddress("mbbgg",&mbbgg);
    tree->SetBranchAddress("ybbgg",&ybbgg);
    tree->SetBranchAddress("ptbbgg",&ptbbgg);
    tree->SetBranchAddress("pho1Pt",&pho1Pt);
    tree->SetBranchAddress("pho1Eta",&pho1Eta);
    tree->SetBranchAddress("pho2Pt",&pho2Pt);
    tree->SetBranchAddress("pho2Eta",&pho2Eta);
    tree->SetBranchAddress("bjet1Pt",&bjet1Pt);
    tree->SetBranchAddress("bjet1Eta",&bjet1Eta);
    tree->SetBranchAddress("bjet2Pt",&bjet2Pt);
    tree->SetBranchAddress("bjet2Eta",&bjet2Eta);

    cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";

    for (int n=0;n<tree->GetEntries();n++) { 
    
      tree->GetEntry(n);
      if (n % 100000 == 0) cout << "Processing Event " << n << "\n";       
      
      if (!(mbb>0)) continue;

      double weight_140PU = 1.0;
      double weightBarrelTiming_140PU = 1.0;      
      double weightEndcapTiming_140PU = 1.0;      
      double weightBarrelEndcapTiming_140PU = 1.0;
      double weight_200PU = 1.0;
      double weightBarrelTiming_200PU = 1.0;
      double weightEndcapTiming_200PU = 1.0;
      double weightBarrelEndcapTiming_200PU = 1.0;
      double weightHGCTiming_200PU = 1.0;

      //cout <<  pho1Eta << " " << pho2Eta << " " << bjet1Eta << " " << bjet2Eta << "\n";

      weight_140PU = weight_140PU * getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist140, false, false);
      weight_140PU = weight_140PU * getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist140, false, false);
      weight_140PU = weight_140PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist140, false, false);
      weight_140PU = weight_140PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist140, false, false);

      //cout << "no timing: " << getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist140, false, false) << " " << getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist140, false, false) << " " << getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist140, false, false) << " " << getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist140, false, false) << "\n";

      weightBarrelTiming_140PU = weightBarrelTiming_140PU * getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist140, true, false);
      weightBarrelTiming_140PU = weightBarrelTiming_140PU * getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist140, true, false);
      weightBarrelTiming_140PU = weightBarrelTiming_140PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist140, true, false);
      weightBarrelTiming_140PU = weightBarrelTiming_140PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist140, true, false);
      weightEndcapTiming_140PU = weightEndcapTiming_140PU * getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist140, false, true);
      weightEndcapTiming_140PU = weightEndcapTiming_140PU * getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist140, false, true);
      weightEndcapTiming_140PU = weightEndcapTiming_140PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist140, false, true);
      weightEndcapTiming_140PU = weightEndcapTiming_140PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist140, false, true);
      weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist140, true, true);
      weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist140, true, true);
      weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist140, true, true);
      weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist140, true, true);

      //cout << "timing: " << getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist140, true, true) << " " << getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist140, true, true) << " " << getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist140, true, true) << " " << getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist140, true, true) << "\n";

      TotalCounts_140PU += weight_140PU;
      TotalCounts_140PU_BarrelTiming += weightBarrelTiming_140PU;
      TotalCounts_140PU_EndcapTiming += weightEndcapTiming_140PU;
      TotalCounts_140PU_BarrelEndcapTiming += weightBarrelEndcapTiming_140PU;
      //cout << weight_140PU << "\n";
     
      weight_200PU = weight_200PU * getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist200, false, false);
      weight_200PU = weight_200PU * getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist200, false, false);
      weight_200PU = weight_200PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist200, false, false);
      weight_200PU = weight_200PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist200, false, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist200, true, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist200, true, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist200, true, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist200, true, false);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist200, false, true);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist200, false, true);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist200, false, true);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist200, false, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist200, true, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist200, true, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist200, true, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist200, true, true);
      weightHGCTiming_200PU = weightHGCTiming_200PU * getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist200, false, false, true);
      weightHGCTiming_200PU = weightHGCTiming_200PU * getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist200, false, false, true);
      weightHGCTiming_200PU = weightHGCTiming_200PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist200, false, false, true);
      weightHGCTiming_200PU = weightHGCTiming_200PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist200, false, false, true);

      TotalCounts_200PU += weight_200PU;
      TotalCounts_200PU_BarrelTiming += weightBarrelTiming_200PU;
      TotalCounts_200PU_EndcapTiming += weightEndcapTiming_200PU;
      TotalCounts_200PU_BarrelEndcapTiming += weightBarrelEndcapTiming_200PU;
      TotalCounts_200PU_HGCTiming += weightHGCTiming_200PU;
      //cout << weight_140PU << "\n";

      hybbgg_140PU->Fill(ybbgg, weight_140PU); 
      hybbggBarrelTiming_140PU->Fill(ybbgg, weightBarrelTiming_140PU); 
      hybbggEndcapTiming_140PU->Fill(ybbgg, weightEndcapTiming_140PU); 
      hybbggBarrelEndcapTiming_140PU->Fill(ybbgg, weightBarrelEndcapTiming_140PU); 
      hybbgg_200PU->Fill(ybbgg, weight_200PU); 
      hybbggBarrelTiming_200PU->Fill(ybbgg, weightBarrelTiming_200PU); 
      hybbggEndcapTiming_200PU->Fill(ybbgg, weightEndcapTiming_200PU); 
      hybbggBarrelEndcapTiming_200PU->Fill(ybbgg, weightBarrelEndcapTiming_200PU);    
      hybbggHGCTiming_200PU->Fill(ybbgg, weightHGCTiming_200PU);    

    }

    inputFile->Close();
    delete inputFile;
  
  }
  
  cout << "Total: " << "\n";
  cout << TotalCounts_140PU << " : " << TotalCounts_140PU_BarrelTiming << " " << TotalCounts_140PU_EndcapTiming << " " << TotalCounts_140PU_BarrelEndcapTiming<< "\n";
  cout << TotalCounts_200PU << " : " << TotalCounts_200PU_HGCTiming << " " << TotalCounts_200PU_BarrelTiming << " " << TotalCounts_200PU_EndcapTiming << " " << TotalCounts_200PU_BarrelEndcapTiming<< "\n";

  //*******************************************************************************************
  //Normalize Hists
  //*******************************************************************************************
  hybbgg_140PU->Scale(1/TotalCounts_140PU);
  hybbggBarrelTiming_140PU->Scale(1/TotalCounts_140PU);
  hybbggEndcapTiming_140PU->Scale(1/TotalCounts_140PU);
  hybbggBarrelEndcapTiming_140PU->Scale(1/TotalCounts_140PU);

  hybbgg_200PU->Scale(1/TotalCounts_200PU);
  hybbggBarrelTiming_200PU->Scale(1/TotalCounts_200PU);
  hybbggEndcapTiming_200PU->Scale(1/TotalCounts_200PU);
  hybbggBarrelEndcapTiming_200PU->Scale(1/TotalCounts_200PU);
  hybbggHGCTiming_200PU->Scale(1/TotalCounts_200PU);

  //*******************************************************************************************
  //Draw Plots
  //*******************************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;
  TLatex *tex = 0;



  //*******************************************************************************************
  //Higgs->bb mass
  //*******************************************************************************************
 
  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend = new TLegend(0.18,0.70,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  legend->AddEntry(hybbgg_140PU, "No Timing");
  legend->AddEntry(hybbggBarrelTiming_140PU, "Barrel Timing Only");
  //legend->AddEntry(hybbggEndcapTiming_140PU, "Endcap Timing Only");
  legend->AddEntry(hybbggBarrelEndcapTiming_140PU, "Barrel+Endcap Timing");

 
  hybbggBarrelEndcapTiming_140PU->Draw("hist");
  hybbggBarrelTiming_140PU->Draw("samehist");
  //hybbggEndcapTiming_140PU->Draw("samehist");
  hybbgg_140PU->Draw("samehist");

  hybbggBarrelEndcapTiming_140PU->GetYaxis()->SetTitle("Fraction of Events");
  hybbggBarrelEndcapTiming_140PU->GetYaxis()->SetRangeUser(0,0.4);
  hybbggBarrelEndcapTiming_140PU->GetYaxis()->SetTitleOffset(1.8);
  hybbggBarrelEndcapTiming_140PU->GetXaxis()->SetTitleOffset(1.3);
  hybbgg_140PU->SetLineColor(kBlack);
  hybbggBarrelTiming_140PU->SetLineColor(kRed);
  hybbggEndcapTiming_140PU->SetLineColor(kGreen+2);
  hybbggBarrelEndcapTiming_140PU->SetLineColor(kBlue);

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.028);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.5, 0.5, "test");
  tex->DrawLatex(0.18, 0.70, "Increase in HH#rightarrow bb#gamma#gamma Yield");
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.18, 0.66, (string(Form("Barrel Timing Only : %.0f", 100*(TotalCounts_140PU_BarrelTiming/TotalCounts_140PU - 1))) + "%").c_str());
  // tex->SetTextColor(kGreen+2);
  // tex->DrawLatex(0.18, 0.62, (string(Form("Endcap Timing Only : %.0f", 100*(TotalCounts_140PU_EndcapTiming/TotalCounts_140PU - 1))) + "%").c_str());
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.18, 0.62, (string(Form("Barrel+Endcap Timing : %.0f", 100*(TotalCounts_140PU_BarrelEndcapTiming/TotalCounts_140PU - 1)))+"%").c_str());

  tex->SetTextSize(0.040);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.18, 0.94, "HH #rightarrow bb#gamma#gamma ( 140 Pileup Distribution )");

  tex->SetTextSize(0.026);
  tex->DrawLatex(0.18, 0.55, "Normalized to");
  tex->DrawLatex(0.18, 0.52, "\"No Timing\" distribution");

  //tex->Draw();

  cv->SaveAs("HHbbgg_ybbgg_TimingStudy_140PU.gif");
  cv->SaveAs("HHbbgg_ybbgg_TimingStudy_140PU.pdf");




  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend = new TLegend(0.18,0.70,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  legend->AddEntry(hybbgg_200PU, "No Timing");
  legend->AddEntry(hybbggHGCTiming_200PU, "Timing (p_{T} > 2 GeV) Only");
  legend->AddEntry(hybbggBarrelTiming_200PU, "Barrel Timing Only");
  legend->AddEntry(hybbggEndcapTiming_200PU, "Endcap Timing Only");
  legend->AddEntry(hybbggBarrelEndcapTiming_200PU, "Barrel+Endcap Timing");

 
  hybbggBarrelEndcapTiming_200PU->Draw("hist");
  hybbggBarrelTiming_200PU->Draw("samehist");
  hybbggEndcapTiming_200PU->Draw("samehist");
  hybbgg_200PU->Draw("samehist");
  hybbggHGCTiming_200PU->Draw("samehist");

  hybbggBarrelEndcapTiming_200PU->GetYaxis()->SetTitle("Fraction of Events");
  hybbggBarrelEndcapTiming_200PU->GetYaxis()->SetRangeUser(0,0.4);
  hybbggBarrelEndcapTiming_200PU->GetYaxis()->SetTitleOffset(1.8);
  hybbggBarrelEndcapTiming_200PU->GetXaxis()->SetTitleOffset(1.3);
  hybbgg_200PU->SetLineColor(kBlack);
  hybbggHGCTiming_200PU->SetLineColor(kViolet);
  hybbggBarrelTiming_200PU->SetLineColor(kRed);
  hybbggEndcapTiming_200PU->SetLineColor(kGreen+2);
  hybbggBarrelEndcapTiming_200PU->SetLineColor(kBlue);

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.028);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.5, 0.5, "test");
  tex->DrawLatex(0.18, 0.66, "Increase in HH#rightarrow bb#gamma#gamma Yield");
  tex->SetTextColor(kViolet);
  tex->DrawLatex(0.18, 0.62, (string(Form("Timing (p_{T} > 2 GeV) Only : %.0f", 100*(TotalCounts_200PU_HGCTiming/TotalCounts_200PU - 1))) + "%").c_str());
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.18, 0.58, (string(Form("Barrel Timing Only : %.0f", 100*(TotalCounts_200PU_BarrelTiming/TotalCounts_200PU - 1))) + "%").c_str());
  tex->SetTextColor(kGreen+2);
  tex->DrawLatex(0.18, 0.54, (string(Form("Endcap Timing Only : %.0f", 100*(TotalCounts_200PU_EndcapTiming/TotalCounts_200PU - 1))) + "%").c_str());
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.18, 0.50, (string(Form("Barrel+Endcap Timing : %.0f", 100*(TotalCounts_200PU_BarrelEndcapTiming/TotalCounts_200PU - 1)))+"%").c_str());

  tex->SetTextSize(0.040);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.18, 0.94, "HH #rightarrow bb#gamma#gamma ( 200 Pileup Distribution )");

  // tex->SetTextSize(0.026);
  // tex->DrawLatex(0.18, 0.55, "Normalized to");
  // tex->DrawLatex(0.18, 0.52, "\"No Timing\" distribution");

  //tex->Draw();

  cv->SaveAs("HHbbgg_ybbgg_TimingStudy_200PU.gif");
  cv->SaveAs("HHbbgg_ybbgg_TimingStudy_200PU.pdf");

  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("TimingStudyPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();
  
  file->WriteTObject(hybbgg_200PU, Form("hybbgg%s",Label.c_str()), "WriteDelete");
 
 
 }


//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void MakeHHbbbbPlots ( string datafile, string dataLabel,  vector<string> bkgfiles,vector<string> bkgLabels, vector<int> bkgColors, int boxOption = 0, int option = -1, string label = "", string latexlabel = "") {

  string Label  = "";
  if(label != "") Label = "_"+label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TFile *pileupDensityFile = TFile::Open("/eos/cms/store/group/phys_susy/razor/MTDUpgradePerformance/lineDensity.root","READ");
  TH1F *pileupDensityHist140 = (TH1F*)pileupDensityFile->Get("hL140");
  TH1F *pileupDensityHist200 = (TH1F*)pileupDensityFile->Get("hL200");

  vector<string> inputfiles;
  vector<string> processLabels;
  vector<int> color;

  bool hasData = false;
  if (datafile != "") {
    hasData = true;
    inputfiles.push_back(datafile);
    processLabels.push_back(dataLabel);
    color.push_back(kBlack);
  } else {
    hasData = false;
    // inputfiles.push_back("");
    // processLabels.push_back("");    
    // color.push_back(kBlack);
  }
  assert(bkgfiles.size() == bkgLabels.size());
  for (int i=0; i < bkgfiles.size(); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processLabels.push_back(bkgLabels[i]);
     color.push_back(bkgColors[i]);
  }


  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  vector<double> EventCount;
  vector<double> EventCountErrSqr;
  vector<TH1F*> histM4l;
  vector<TH1F*> histMR;
  vector<TH1F*> histRsq;

  TH1F *hybbgg_140PU = makeNewHist( "hybbgg_140PU", ";y_{HH} [GeV/c^{2}];Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggBarrelEndcapTiming_140PU = makeNewHist( "hybbggBarrelEndcapTiming_140PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggBarrelTiming_140PU = makeNewHist( "hybbggBarrelTiming_140PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggEndcapTiming_140PU = makeNewHist( "hybbggEndcapTiming_140PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbgg_200PU = makeNewHist( "hybbgg_200PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggBarrelEndcapTiming_200PU = makeNewHist( "hybbggBarrelEndcapTiming_200PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggBarrelTiming_200PU = makeNewHist( "hybbggBarrelTiming_200PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggEndcapTiming_200PU = makeNewHist( "hybbggEndcapTiming_200PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );
  TH1F *hybbggHGCTiming_200PU = makeNewHist( "hybbggHGCTiming_200PU", ";y_{HH} ;Number of Events", 20, -3, 3, kBlue, false );


  //*******************************************************************************************
  //Define Counts
  //*******************************************************************************************
  double TotalCounts_140PU = 0;
  double TotalCounts_140PU_BarrelTiming = 0;
  double TotalCounts_140PU_EndcapTiming = 0;
  double TotalCounts_140PU_BarrelEndcapTiming = 0;
  double TotalCounts_200PU = 0;
  double TotalCounts_200PU_BarrelTiming = 0;
  double TotalCounts_200PU_EndcapTiming = 0;
  double TotalCounts_200PU_BarrelEndcapTiming = 0;
  double TotalCounts_200PU_HGCTiming = 0;


  //*******************************************************************************************
  //Read files
  //*******************************************************************************************
  for (uint i=0; i < inputfiles.size(); ++i) {

    TFile* inputFile = new TFile(inputfiles[i].c_str(),"READ");
    assert(inputFile);
    TTree* tree = 0;
    tree = (TTree*)inputFile->Get("HHbbgg");    
    assert(tree);

    float mbb = 0;
    float mgg = 0;
    float mbbgg = 0;
    float ybbgg = 0;
    float ptbbgg = 0;
    float pho1Pt = 0;
    float pho1Eta = 0;
    float pho2Pt = 0;
    float pho2Eta = 0;
    float bjet1Pt = 0;
    float bjet1Eta = 0;
    float bjet2Pt = 0;
    float bjet2Eta = 0;

    tree->SetBranchAddress("mbb",&mbb);
    tree->SetBranchAddress("mgg",&mgg);
    tree->SetBranchAddress("mbbgg",&mbbgg);
    tree->SetBranchAddress("ybbgg",&ybbgg);
    tree->SetBranchAddress("ptbbgg",&ptbbgg);
    tree->SetBranchAddress("pho1Pt",&pho1Pt);
    tree->SetBranchAddress("pho1Eta",&pho1Eta);
    tree->SetBranchAddress("pho2Pt",&pho2Pt);
    tree->SetBranchAddress("pho2Eta",&pho2Eta);
    tree->SetBranchAddress("bjet1Pt",&bjet1Pt);
    tree->SetBranchAddress("bjet1Eta",&bjet1Eta);
    tree->SetBranchAddress("bjet2Pt",&bjet2Pt);
    tree->SetBranchAddress("bjet2Eta",&bjet2Eta);

    cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";

    for (int n=0;n<tree->GetEntries();n++) { 
    
      tree->GetEntry(n);
      if (n % 100000 == 0) cout << "Processing Event " << n << "\n";       
      
      if (!(mbb>0)) continue;

      double weight_140PU = 1.0;
      double weightBarrelTiming_140PU = 1.0;      
      double weightEndcapTiming_140PU = 1.0;      
      double weightBarrelEndcapTiming_140PU = 1.0;
      double weight_200PU = 1.0;
      double weightBarrelTiming_200PU = 1.0;
      double weightEndcapTiming_200PU = 1.0;
      double weightBarrelEndcapTiming_200PU = 1.0;
      double weightHGCTiming_200PU = 1.0;

      //cout <<  pho1Eta << " " << pho2Eta << " " << bjet1Eta << " " << bjet2Eta << "\n";

      weight_140PU = weight_140PU * getBTagEffFactor( pho1Pt, pho1Eta, pileupDensityHist140, false, false);
      weight_140PU = weight_140PU * getBTagEffFactor( pho2Pt, pho2Eta, pileupDensityHist140, false, false);
      weight_140PU = weight_140PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist140, false, false);
      weight_140PU = weight_140PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist140, false, false);

      //cout << "no timing: " << getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist140, false, false) << " " << getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist140, false, false) << " " << getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist140, false, false) << " " << getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist140, false, false) << "\n";

      weightBarrelTiming_140PU = weightBarrelTiming_140PU * getBTagEffFactor( pho1Pt, pho1Eta, pileupDensityHist140, true, false);
      weightBarrelTiming_140PU = weightBarrelTiming_140PU * getBTagEffFactor( pho2Pt, pho2Eta, pileupDensityHist140, true, false);
      weightBarrelTiming_140PU = weightBarrelTiming_140PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist140, true, false);
      weightBarrelTiming_140PU = weightBarrelTiming_140PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist140, true, false);
      weightEndcapTiming_140PU = weightEndcapTiming_140PU * getBTagEffFactor( pho1Pt, pho1Eta, pileupDensityHist140, false, true);
      weightEndcapTiming_140PU = weightEndcapTiming_140PU * getBTagEffFactor( pho2Pt, pho2Eta, pileupDensityHist140, false, true);
      weightEndcapTiming_140PU = weightEndcapTiming_140PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist140, false, true);
      weightEndcapTiming_140PU = weightEndcapTiming_140PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist140, false, true);
      weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getBTagEffFactor( pho1Pt, pho1Eta, pileupDensityHist140, true, true);
      weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getBTagEffFactor( pho2Pt, pho2Eta, pileupDensityHist140, true, true);
      weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist140, true, true);
      weightBarrelEndcapTiming_140PU = weightBarrelEndcapTiming_140PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist140, true, true);

      //cout << "timing: " << getPhoEffFactor( pho1Pt, pho1Eta, pileupDensityHist140, true, true) << " " << getPhoEffFactor( pho2Pt, pho2Eta, pileupDensityHist140, true, true) << " " << getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist140, true, true) << " " << getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist140, true, true) << "\n";

      TotalCounts_140PU += weight_140PU;
      TotalCounts_140PU_BarrelTiming += weightBarrelTiming_140PU;
      TotalCounts_140PU_EndcapTiming += weightEndcapTiming_140PU;
      TotalCounts_140PU_BarrelEndcapTiming += weightBarrelEndcapTiming_140PU;
      //cout << weight_140PU << "\n";
     
      weight_200PU = weight_200PU * getBTagEffFactor( pho1Pt, pho1Eta, pileupDensityHist200, false, false);
      weight_200PU = weight_200PU * getBTagEffFactor( pho2Pt, pho2Eta, pileupDensityHist200, false, false);
      weight_200PU = weight_200PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist200, false, false);
      weight_200PU = weight_200PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist200, false, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getBTagEffFactor( pho1Pt, pho1Eta, pileupDensityHist200, true, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getBTagEffFactor( pho2Pt, pho2Eta, pileupDensityHist200, true, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist200, true, false);
      weightBarrelTiming_200PU = weightBarrelTiming_200PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist200, true, false);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getBTagEffFactor( pho1Pt, pho1Eta, pileupDensityHist200, false, true);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getBTagEffFactor( pho2Pt, pho2Eta, pileupDensityHist200, false, true);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist200, false, true);
      weightEndcapTiming_200PU = weightEndcapTiming_200PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist200, false, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getBTagEffFactor( pho1Pt, pho1Eta, pileupDensityHist200, true, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getBTagEffFactor( pho2Pt, pho2Eta, pileupDensityHist200, true, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist200, true, true);
      weightBarrelEndcapTiming_200PU = weightBarrelEndcapTiming_200PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist200, true, true);
      weightHGCTiming_200PU = weightHGCTiming_200PU * getBTagEffFactor( pho1Pt, pho1Eta, pileupDensityHist200, false, false, true);
      weightHGCTiming_200PU = weightHGCTiming_200PU * getBTagEffFactor( pho2Pt, pho2Eta, pileupDensityHist200, false, false, true);
      weightHGCTiming_200PU = weightHGCTiming_200PU * getBTagEffFactor( bjet1Pt, bjet1Eta, pileupDensityHist200, false, false, true);
      weightHGCTiming_200PU = weightHGCTiming_200PU * getBTagEffFactor( bjet2Pt, bjet2Eta, pileupDensityHist200, false, false, true);

      TotalCounts_200PU += weight_200PU;
      TotalCounts_200PU_BarrelTiming += weightBarrelTiming_200PU;
      TotalCounts_200PU_EndcapTiming += weightEndcapTiming_200PU;
      TotalCounts_200PU_BarrelEndcapTiming += weightBarrelEndcapTiming_200PU;
      TotalCounts_200PU_HGCTiming += weightHGCTiming_200PU;
      //cout << weight_140PU << "\n";

      hybbgg_140PU->Fill(ybbgg, weight_140PU); 
      hybbggBarrelTiming_140PU->Fill(ybbgg, weightBarrelTiming_140PU); 
      hybbggEndcapTiming_140PU->Fill(ybbgg, weightEndcapTiming_140PU); 
      hybbggBarrelEndcapTiming_140PU->Fill(ybbgg, weightBarrelEndcapTiming_140PU); 
      hybbgg_200PU->Fill(ybbgg, weight_200PU); 
      hybbggBarrelTiming_200PU->Fill(ybbgg, weightBarrelTiming_200PU); 
      hybbggEndcapTiming_200PU->Fill(ybbgg, weightEndcapTiming_200PU); 
      hybbggBarrelEndcapTiming_200PU->Fill(ybbgg, weightBarrelEndcapTiming_200PU);    
      hybbggHGCTiming_200PU->Fill(ybbgg, weightHGCTiming_200PU);    

    }

    inputFile->Close();
    delete inputFile;
  
  }
  
  cout << "Total: " << "\n";
  cout << TotalCounts_140PU << " : " << TotalCounts_140PU_BarrelTiming << " " << TotalCounts_140PU_EndcapTiming << " " << TotalCounts_140PU_BarrelEndcapTiming<< "\n";
  cout << TotalCounts_200PU << " : " << TotalCounts_200PU_HGCTiming << " " << TotalCounts_200PU_BarrelTiming << " " << TotalCounts_200PU_EndcapTiming << " " << TotalCounts_200PU_BarrelEndcapTiming<< "\n";

  //*******************************************************************************************
  //Normalize Hists
  //*******************************************************************************************
  hybbgg_140PU->Scale(1/TotalCounts_140PU);
  hybbggBarrelTiming_140PU->Scale(1/TotalCounts_140PU);
  hybbggEndcapTiming_140PU->Scale(1/TotalCounts_140PU);
  hybbggBarrelEndcapTiming_140PU->Scale(1/TotalCounts_140PU);

  hybbgg_200PU->Scale(1/TotalCounts_200PU);
  hybbggBarrelTiming_200PU->Scale(1/TotalCounts_200PU);
  hybbggEndcapTiming_200PU->Scale(1/TotalCounts_200PU);
  hybbggBarrelEndcapTiming_200PU->Scale(1/TotalCounts_200PU);
  hybbggHGCTiming_200PU->Scale(1/TotalCounts_200PU);

  //*******************************************************************************************
  //Draw Plots
  //*******************************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;
  TLatex *tex = 0;



  //*******************************************************************************************
  //Higgs->bb mass
  //*******************************************************************************************
 
  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend = new TLegend(0.18,0.75,0.58,0.88);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  legend->AddEntry(hybbgg_140PU, "No Timing");
  legend->AddEntry(hybbggBarrelTiming_140PU, "Barrel Timing Only");
  // legend->AddEntry(hybbggEndcapTiming_140PU, "Endcap Timing Only");
  legend->AddEntry(hybbggBarrelEndcapTiming_140PU, "Barrel+Endcap Timing");

 
  hybbggBarrelEndcapTiming_140PU->Draw("hist");
  hybbggBarrelTiming_140PU->Draw("samehist");
  // hybbggEndcapTiming_140PU->Draw("samehist");
  hybbgg_140PU->Draw("samehist");

  hybbggBarrelEndcapTiming_140PU->GetYaxis()->SetTitle("Fraction of Events");
  hybbggBarrelEndcapTiming_140PU->GetYaxis()->SetRangeUser(0,0.4);
  hybbggBarrelEndcapTiming_140PU->GetYaxis()->SetTitleOffset(1.8);
  hybbggBarrelEndcapTiming_140PU->GetXaxis()->SetTitleOffset(1.3);
  hybbgg_140PU->SetLineColor(kBlack);
  hybbggBarrelTiming_140PU->SetLineColor(kRed);
  hybbggEndcapTiming_140PU->SetLineColor(kGreen+2);
  hybbggBarrelEndcapTiming_140PU->SetLineColor(kBlue);

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.028);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.5, 0.5, "test");
  tex->DrawLatex(0.18, 0.70, "Increase in HH#rightarrow bbbb Yield");
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.18, 0.66, (string(Form("Barrel Timing Only : %.0f", 100*(TotalCounts_140PU_BarrelTiming/TotalCounts_140PU - 1))) + "%").c_str());
  // tex->SetTextColor(kGreen+2);
  // tex->DrawLatex(0.18, 0.62, (string(Form("Endcap Timing Only : %.0f", 100*(TotalCounts_140PU_EndcapTiming/TotalCounts_140PU - 1))) + "%").c_str());
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.18, 0.62, (string(Form("Barrel+Endcap Timing : %.0f", 100*(TotalCounts_140PU_BarrelEndcapTiming/TotalCounts_140PU - 1)))+"%").c_str());

  tex->SetTextSize(0.040);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.18, 0.94, "HH #rightarrow bbbb ( 140 Pileup Distribution )");

  // tex->SetTextSize(0.026);
  // tex->DrawLatex(0.18, 0.55, "Normalized to");
  // tex->DrawLatex(0.18, 0.52, "\"No Timing\" distribution");

  //tex->Draw();

  cv->SaveAs("HHbbbb_ybbbb_TimingStudy_140PU.gif");
  cv->SaveAs("HHbbbb_ybbbb_TimingStudy_140PU.pdf");




  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend = new TLegend(0.18,0.70,0.60,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->SetFillStyle(0);
  legend->AddEntry(hybbgg_200PU, "No Timing");
  legend->AddEntry(hybbggHGCTiming_200PU, "Timing (p_{T} > 2 GeV) Only");
  legend->AddEntry(hybbggBarrelTiming_200PU, "Barrel Timing Only");
  legend->AddEntry(hybbggEndcapTiming_200PU, "Endcap Timing Only");
  legend->AddEntry(hybbggBarrelEndcapTiming_200PU, "Barrel+Endcap Timing");

 
  hybbggBarrelEndcapTiming_200PU->Draw("hist");
  hybbggBarrelTiming_200PU->Draw("samehist");
  hybbggEndcapTiming_200PU->Draw("samehist");
  hybbgg_200PU->Draw("samehist");
  hybbggHGCTiming_200PU->Draw("samehist");

  hybbggBarrelEndcapTiming_200PU->GetYaxis()->SetTitle("Fraction of Events");
  hybbggBarrelEndcapTiming_200PU->GetYaxis()->SetRangeUser(0,0.4);
  hybbggBarrelEndcapTiming_200PU->GetYaxis()->SetTitleOffset(1.8);
  hybbggBarrelEndcapTiming_200PU->GetXaxis()->SetTitleOffset(1.3);
  hybbgg_200PU->SetLineColor(kBlack);
  hybbggHGCTiming_200PU->SetLineColor(kViolet);
  hybbggBarrelTiming_200PU->SetLineColor(kRed);
  hybbggEndcapTiming_200PU->SetLineColor(kGreen+2);
  hybbggBarrelEndcapTiming_200PU->SetLineColor(kBlue);

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.028);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.5, 0.5, "test");
  tex->DrawLatex(0.18, 0.66, "Increase in HH#rightarrow bbbb Yield");
  tex->SetTextColor(kViolet);
  tex->DrawLatex(0.18, 0.62, (string(Form("Timing (p_{T} > 2 GeV) Only : %.0f", 100*(TotalCounts_200PU_HGCTiming/TotalCounts_200PU - 1))) + "%").c_str());
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.18, 0.58, (string(Form("Barrel Timing Only : %.0f", 100*(TotalCounts_200PU_BarrelTiming/TotalCounts_200PU - 1))) + "%").c_str());
  tex->SetTextColor(kGreen+2);
  tex->DrawLatex(0.18, 0.54, (string(Form("Endcap Timing Only : %.0f", 100*(TotalCounts_200PU_EndcapTiming/TotalCounts_200PU - 1))) + "%").c_str());
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.18, 0.50, (string(Form("Barrel+Endcap Timing : %.0f", 100*(TotalCounts_200PU_BarrelEndcapTiming/TotalCounts_200PU - 1)))+"%").c_str());

  tex->SetTextSize(0.040);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.18, 0.94, "HH #rightarrow bbbb ( 200 Pileup Distribution )");

  // tex->SetTextSize(0.026);
  // tex->DrawLatex(0.18, 0.55, "Normalized to");
  // tex->DrawLatex(0.18, 0.52, "\"No Timing\" distribution");

  //tex->Draw();

  cv->SaveAs("HHbbbb_ybbbb_TimingStudy_200PU.gif");
  cv->SaveAs("HHbbbb_ybbbb_TimingStudy_200PU.pdf");

  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("TimingStudyPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();
  
  file->WriteTObject(hybbgg_200PU, Form("hybbbb%s",Label.c_str()), "WriteDelete");
 
 
 }



void MakeHiggsImprovementVsPileupPlot() {

  const int nPoints = 20;
  float linDensity[nPoints];
  float effImprovementHZZ[nPoints];
  float effImprovementHMM[nPoints];

  for (int i=0; i<nPoints; i++) {

    linDensity[i] = 0.3 + i * (2.0 - 0.3)/20;
    double pileupDegradationFactor = 0.95 - (i*(2.0-0.3)/20)*(0.95-0.79)/(1.5-0.3);
    double pileupDegradationFactorWithTiming = 0.95 - (i*(2.0-0.3)/20)*(0.95-0.91)/(1.5-0.3);
 
    effImprovementHZZ[i] = 100* (pow( pileupDegradationFactorWithTiming / pileupDegradationFactor , 4) - 1);
    effImprovementHMM[i] = 100 * (pow( pileupDegradationFactorWithTiming / pileupDegradationFactor , 2) - 1);
  }

  TGraph *graphHZZ = new TGraph( nPoints, linDensity, effImprovementHZZ);
  TGraph *graphHMM = new TGraph( nPoints, linDensity, effImprovementHMM);



  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;
  TLatex *tex = 0;
  TLine *line = 0;

  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend = new TLegend(0.20,0.75,0.60,0.88);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(graphHZZ, "Higgs #rightarrow ZZ #rightarrow 4l");
  legend->AddEntry(graphHMM, "Higgs #rightarrow #mu#mu");
 
  graphHZZ->Draw("AP");
  graphHMM->Draw("P");
  legend->Draw();

  graphHZZ->SetFillStyle(0);
  graphHZZ->SetLineColor(kRed);
  graphHZZ->SetLineWidth(0);
  graphHZZ->SetMarkerColor(kRed);
  graphHZZ->SetMarkerStyle(20);
  graphHZZ->SetMarkerSize(2);
  graphHMM->SetFillStyle(0);
  graphHMM->SetLineColor(kBlue);
  graphHMM->SetLineWidth(0);
  graphHMM->SetMarkerColor(kBlue);
  graphHMM->SetMarkerStyle(21);
  graphHMM->SetMarkerSize(2);

  graphHZZ->SetTitle("");
  graphHZZ->GetXaxis()->SetTitle("Linear Pileup Density (events / mm)");
  graphHZZ->GetXaxis()->SetTitleSize(0.045);
  graphHZZ->GetXaxis()->SetTitleOffset(1.1);
  graphHZZ->GetYaxis()->SetTitle("Relative Increase in Effective Luminosity (%)");
  graphHZZ->GetYaxis()->SetTitleOffset(1.2);
  graphHZZ->GetYaxis()->SetTitleSize(0.045);

  // line = new TLine(0.14, 27.3, 2.08, 27.3);
  // line->SetLineColor(kRed);
  // line->SetLineWidth(4);
  // line->SetLineStyle(2);
  // line->Draw();
  line = new TLine(0.14, 62, 2.08, 62);
  line->SetLineColor(kRed);
  line->SetLineWidth(4);
  line->SetLineStyle(2);
  line->Draw();
  line = new TLine(0.14, 27, 2.08, 27);
  line->SetLineColor(kBlue);
  line->SetLineWidth(4);
  line->SetLineStyle(2);
  line->Draw();
  // line = new TLine(0.14, 13, 2.08, 13);
  // line->SetLineColor(kBlue);
  // line->SetLineWidth(4);
  // line->SetLineStyle(2);
  // line->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kRed);
  // tex->DrawLatex(0.18, 0.32, "140 Pileup Average");
  // tex->DrawLatex(0.18, 0.29, "(Higgs #rightarrow ZZ #rightarrow 4l)");
  tex->DrawLatex(0.18, 0.53, "200 Pileup Average ");
  tex->DrawLatex(0.18, 0.50, "(Higgs #rightarrow ZZ #rightarrow 4l)");
  tex->SetTextColor(kBlue);
  // tex->DrawLatex(0.65, 0.17, "140 Pileup Average");
  // tex->DrawLatex(0.65, 0.14, "(Higgs #rightarrow #mu#mu)");
  tex->DrawLatex(0.65, 0.25, "200 Pileup Average ");
  tex->DrawLatex(0.65, 0.22, "(Higgs #rightarrow #mu#mu)");

  tex->SetTextSize(0.040);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.15, 0.94, "Improvement in Higgs yield with Timing");


  cv->SaveAs("HiggsYieldIncreaseVsLinearDensity.gif");
  cv->SaveAs("HiggsYieldIncreaseVsLinearDensity.pdf");


}





 void RunMakeHZZPlots() {

   string datafile = "";
   string dataLabel = "";  

   vector<string> bkgfiles;
   vector<string> bkgLabels;
   vector<int> bkgColors;

   bkgfiles.push_back("/eos/cms/store/group/phys_susy/razor/Run2Analysis/HZZRazor/HZZRazor_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root");
   bkgLabels.push_back("Higgs");
   bkgColors.push_back(kBlue);

   MakeHZZPlots(datafile,dataLabel,bkgfiles,bkgLabels,bkgColors,0,1,"HZZ","HZZ");
 
  
 
 }
 
void RunMakeHMMPlots() {

   string datafile = "";
   string dataLabel = "";  

   vector<string> bkgfiles;
   vector<string> bkgLabels;
   vector<int> bkgColors;

   bkgfiles.push_back("/eos/cms/store/group/phys_susy/razor/Run2Analysis/HZZRazor/HZZRazor_GluGlu_HToMuMu_M125_13TeV_powheg_pythia8.root");
   bkgLabels.push_back("Higgs");
   bkgColors.push_back(kBlue);

   MakeHMMPlots(datafile,dataLabel,bkgfiles,bkgLabels,bkgColors,0,1,"HMM","HMM");
 
  
 
 }
 
void RunMakeHHbbggPlots() {

   string datafile = "";
   string dataLabel = "";  

   vector<string> bkgfiles;
   vector<string> bkgLabels;
   vector<int> bkgColors;

   bkgfiles.push_back("/eos/cms/store/group/phys_susy/razor/Run2Analysis/HHbbgg/HHbbgg_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root");
   bkgLabels.push_back("HH");
   bkgColors.push_back(kBlue);

   MakeHHbbggPlots(datafile,dataLabel,bkgfiles,bkgLabels,bkgColors,0,1,"HH","HH");
 
 }
 
void RunMakeHHbbbbPlots() {

   string datafile = "";
   string dataLabel = "";  

   vector<string> bkgfiles;
   vector<string> bkgLabels;
   vector<int> bkgColors;

   bkgfiles.push_back("/eos/cms/store/group/phys_susy/razor/Run2Analysis/HHbbgg/HHbbgg_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root");
   bkgLabels.push_back("HH");
   bkgColors.push_back(kBlue);

   MakeHHbbbbPlots(datafile,dataLabel,bkgfiles,bkgLabels,bkgColors,0,1,"HH","HH");
 
 }
 

void MakeTimingStudyPlots() {
  RunMakeHZZPlots();
  //RunMakeHMMPlots();  

  //RunMakeHHbbggPlots();
  //RunMakeHHbbbbPlots();

  // MakeHiggsImprovementVsPileupPlot();

}
