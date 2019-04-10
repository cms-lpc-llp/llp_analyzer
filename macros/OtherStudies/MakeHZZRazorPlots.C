
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <vector>
#include <map>
#include <iostream>

const Int_t NComponents = 10;
int color[NComponents] = {kRed, kGreen+2, kBlue, kViolet, kAzure+10, kGray, kOrange+1, kGray+3, kBlack, kBlack};


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

TH1F* makeNewHist( string title, string axisLabel, int nbins, double rangeLow, double rangeHigh, int color, bool isData) {

  TH1F *newhist = new TH1F( title.c_str(), axisLabel.c_str(), nbins, rangeLow, rangeHigh);
  if (!isData) newhist->SetFillColor(color);
  if (isData) newhist->SetLineWidth(3); 
  newhist->SetLineColor(color);    
  newhist->SetStats(false);    
  newhist->Sumw2();
  return newhist;

}

void PlotDataAndStackedBkg( vector<TH1F*> hist , vector<string> processLabels, vector<int> color,  bool hasData, string varName, string label ) {

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

  legend = new TLegend(0.60,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stack = new THStack();
  TH1F *histDataOverMC = (TH1F*)hist[0]->Clone("histDataOverMC");

  if (hasData) {
    for (int i = hist.size()-1; i >= 1; --i) {
      hist[i]->SetFillColor(color[i]);
      hist[i]->SetFillStyle(1001);
      double intError = 0;
      for(int j=1; j < hist[i]->GetNbinsX()+1; j++) {
	intError += pow(hist[i]->GetBinError(j),2);
      }
      intError = sqrt(intError);
      //cout << processLabels[i] << " : " << hist[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( hist[i]->Integral() > 0) {
  	stack->Add(hist[i]);
      }
    }
  } else {
    for (int i = hist.size()-1; i >= 0; --i) {
      hist[i]->SetFillColor(color[i]);
      hist[i]->SetFillStyle(1001);
      double intError = 0;
      for(int j=1; j < hist[i]->GetNbinsX()+1; j++) {
	intError += pow(hist[i]->GetBinError(j),2);
      }
      intError = sqrt(intError);
      //cout << processLabels[i] << " : " << hist[i]->GetSumOfWeights() << " +/- " << intError << "\n";
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
    stack->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stack->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stack->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stack->GetHists()->At(0)))->GetYaxis()->GetTitle());
    stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    stack->SetMaximum( 2* fmax( stack->GetMaximum(), hist[0]->GetMaximum()) );
    stack->SetMinimum( 0.001 );

    if (hasData) {
      hist[0]->SetLineWidth(2);
      hist[0]->SetLineColor(color[0]);
      hist[0]->Draw("e1same");
    }
    legend->Draw();
  }
  cv->cd();
  cv->Update();


  // if (hasData) {
  //   cout << processLabels[0] << " : " << hist[0]->GetSumOfWeights() << "\n";
  // }



  TPad *pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
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
  histDataOverMC->GetXaxis()->SetTitleSize(0.125);
  histDataOverMC->GetXaxis()->SetTitleOffset(1.2);
  histDataOverMC->SetStats(false);
  histDataOverMC->Draw("e1");
  
  pad1->SetLogy(false);
  cv->SaveAs(Form("HggRazor_%s%s.gif",varName.c_str(), label.c_str()));
  
  pad1->SetLogy(true);
  cv->SaveAs(Form("HggRazor_%s%s_Logy.gif",varName.c_str(),label.c_str()));

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
void RunMakeRazorPlots ( string datafile, string dataLabel,  vector<string> bkgfiles,vector<string> bkgLabels, vector<int> bkgColors, int boxOption = 0, int option = -1, string label = "", string latexlabel = "") {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TRandom3 *random = new TRandom3(0);
  double intLumi = 19700; //in units of pb^-1
  string Label = "";
  if (label != "") Label = "_" + label;

  TFile *pileupWeightFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Run1PileupWeights.root", "READ");
  TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("PUWeight_Run1");
  assert(pileupWeightHist);

  TFile *eleEffSFFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/ScaleFactors/ElectronSelection_Run2012ReReco_53X.root","READ");
  TH2D *eleEffSFHist = (TH2D*)eleEffSFFile->Get("sfLOOSE");
  assert(eleEffSFHist);


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
    hasData = true;
    inputfiles.push_back("");
    processLabels.push_back("Hypothetical Data");    
    color.push_back(kBlack);
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

  // vector<TH1F*> histNJets80;
  // vector<TH1F*> histNJets60;
  // vector<TH1F*> histNJets40;
  // vector<TH1F*> histNJets30;
  // vector<TH1F*> histPhoton1Pt;
  // vector<TH1F*> histPhoton2Pt;
  // vector<TH1F*> histJet1Pt;
  // vector<TH1F*> histJet2Pt;
  // vector<TH1F*> histDPhiRazor;
  // vector<TH1F*> histDPhiHiggsMET;  
  // vector<TH1F*> histPtggPeakRegion;
  // vector<TH1F*> histPtggSidebandRegion;
  // vector<TH1F*> histMinDRJetsToPhotons;
  // vector<TH1F*> histMHiggsClosestJet;
  // vector<TH1F*> histMHiggsLeadingJet;
  // vector<TH1F*> histMTLeadJet;
  // vector<TH1F*> histMTMinDPhiMetJet;
  // vector<TH1F*> histDPhiLeadJetMET;
  // vector<TH1F*> histDRLeadJetHiggs;
  // vector<TH1F*> histMjjMaxPtPair;
  // vector<TH1F*> histMjjMaxPtPairExcludeJetClosestToHiggs;
  // vector<TH1F*> histMTJetMaxPt;
  // vector<TH1F*> histMTJJMaxPt;

  // vector<TH1F*> histMassOppositeHiggs;
  // vector<TH1F*> histMTOppositeHiggs;
  // vector<TH1F*> histMassOppositeLeadJet;
  // vector<TH1F*> histMTOppositeLeadJet;

  // vector<TH1F*> histCosThetaStar;



  assert (inputfiles.size() == processLabels.size());
  for (int i=0; i < inputfiles.size(); ++i) {        
    EventCount.push_back(0);
    EventCountErrSqr.push_back(0);
    histM4l.push_back( makeNewHist( Form("M4l_%s",processLabels[i].c_str()), ";M_{4l} [GeV/c^{2}];Number of Events", 20, 100, 160, color[i],  (hasData && i==0) ));  
    histMR.push_back( makeNewHist( Form("MR_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}];Number of Events", 25, 0, 500, color[i],  (hasData && i==0) ));  
    histRsq.push_back( makeNewHist( Form("Rsq_%s",processLabels[i].c_str()), ";R^{2} ;Number of Events", 50, 0, 0.4, color[i],  (hasData && i==0) ));  

    // histMjjMaxPtPair.push_back( makeNewHist( Form("MjjMaxPtPair_%s",processLabels[i].c_str()), ";M_{jj}^{max pTjj} [GeV/c^{2}];Number of Events", 100, 0, 1000, color[i],  (hasData && i==0) ));  
    // histMjjMaxPtPairExcludeJetClosestToHiggs.push_back( makeNewHist( Form("MjjMaxPtPairExcludeJetClosestToHiggs_%s",processLabels[i].c_str()), ";M_{jj}^{max pTjj excl jet nearest Higgs} [GeV/c^{2}];Number of Events", 100, 0, 1000, color[i],  (hasData && i==0) ));  
    // histMassOppositeHiggs.push_back( makeNewHist( Form("histMassOppositeHiggs_%s",processLabels[i].c_str()), ";M_{all objects except #gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000, color[i],  (hasData && i==0) ));  
    // histMTOppositeHiggs.push_back( makeNewHist( Form("histMTOppositeHiggs_%s",processLabels[i].c_str()), ";M_{T}^{all objects except #gamma#gamma} [GeV/c^{2}];Number of Events", 100, 0, 1000, color[i],  (hasData && i==0) ));  
    // histMassOppositeLeadJet.push_back( makeNewHist( Form("histMassOppositeLeadJet_%s",processLabels[i].c_str()), ";M_{all objects except lead jet} [GeV/c^{2}];Number of Events", 100, 0, 1000, color[i],  (hasData && i==0) ));  
    // histMTOppositeLeadJet.push_back( makeNewHist( Form("histMTOppositeLeadJet_%s",processLabels[i].c_str()), ";M_{T}^{all objects except lead jet} [GeV/c^{2}];Number of Events", 100, 0, 1000, color[i],  (hasData && i==0) ));  
    // histMTJetMaxPt.push_back( makeNewHist( Form("MTJetMaxPt_%s",processLabels[i].c_str()), ";M_{jj}^{max pTjj} [GeV/c^{2}];Number of Events", 100, 0, 1000, color[i],  (hasData && i==0) ));  
    // histMTJJMaxPt.push_back( makeNewHist( Form("MTJJMaxPt_%s",processLabels[i].c_str()), ";M_{jj}^{max pTjj} [GeV/c^{2}];Number of Events", 100, 0, 1000, color[i],  (hasData && i==0) ));  

    // histMinDRJetsToPhotons.push_back( makeNewHist( Form("MinDRJetsToPhotons_%s",processLabels[i].c_str()), ";min#Delta R(jet,#gamma);Number of Events", 100, 0, 4.0, color[i], (hasData && i==0) ));  
    // histDRLeadJetHiggs.push_back( makeNewHist( Form("DRLeadJetHiggs_%s",processLabels[i].c_str()), ";#Delta R(Leading jet, Higgs);Number of Events", 100, 0, 10, color[i], (hasData && i==0) ));  
    // histMTMinDPhiMetJet.push_back( makeNewHist( Form("MTMinDPhiMetJet_%s",processLabels[i].c_str()), ";M_{T}^{jet nearest MET} [GeV/c^{2}];Number of Events", 100, 0, 1000, color[i],  (hasData && i==0) ));  
    // histDPhiLeadJetMET.push_back( makeNewHist( Form("DPhiLeadJetMET_%s",processLabels[i].c_str()), ";#Delta#phi(Leading jet, MET);Number of Events", 100, 0, 3.15, color[i],  (hasData && i==0) ));  
    // histMHiggsClosestJet.push_back( makeNewHist( Form("MHiggsClosestJet_%s",processLabels[i].c_str()), ";M_{#gamma#gamma nearest jet} [GeV/c^{2}];Number of Events", 100, 0, 1000, color[i],  (hasData && i==0) ));  
    // histMHiggsLeadingJet.push_back( makeNewHist( Form("MHiggsLeadingJet_%s",processLabels[i].c_str()), ";M_{#gamma#gamma lead jet} [GeV/c^{2}];Number of Events", 100, 0, 1000, color[i],  (hasData && i==0) ));  
    // histMTLeadJet.push_back( makeNewHist( Form("MTLeadJet_%s",processLabels[i].c_str()), ";M_{T}^{leading jet} [GeV/c^{2}];Number of Events", 100, 0, 1000, color[i],  (hasData && i==0) ));  

    // histNJets80.push_back( makeNewHist( Form("NJets80_%s",processLabels[i].c_str()), ";Number of Jets with p_{T} > 80 GeV/c;Number of Events", 10, -0.5, 9.5, color[i],  (hasData && i==0) ));  
    // histNJets60.push_back( makeNewHist( Form("NJets60_%s",processLabels[i].c_str()), ";Number of Jets with p_{T} > 80 GeV/c;Number of Events", 10, -0.5, 9.5, color[i],  (hasData && i==0) ));  
    // histNJets40.push_back( makeNewHist( Form("NJets40_%s",processLabels[i].c_str()), ";Number of Jets with p_{T} > 80 GeV/c;Number of Events", 10, -0.5, 9.5, color[i],  (hasData && i==0) ));   
    // histNJets30.push_back( makeNewHist( Form("NJets30_%s",processLabels[i].c_str()), ";Number of Jets with p_{T} > 80 GeV/c;Number of Events", 10, -0.5, 9.5, color[i],  (hasData && i==0) ));  
    // histPhoton1Pt.push_back( makeNewHist( Form("Photon1Pt_%s",processLabels[i].c_str()), ";Leading Photon p_{T} [GeV/c];Number of Events", 50, 0, 200, color[i],  (hasData && i==0) ));  
    // histPhoton2Pt.push_back( makeNewHist( Form("Photon2Pt_%s",processLabels[i].c_str()), ";Second Photon p_{T} [GeV/c];Number of Events", 50, 0, 200, color[i],  (hasData && i==0) ));  
    // histJet1Pt.push_back( makeNewHist( Form("Jet1Pt_%s",processLabels[i].c_str()), ";Leading Jet p_{T} [GeV/c];Number of Events", 125, 0, 500, color[i],  (hasData && i==0) ));  
    // histJet2Pt.push_back( makeNewHist( Form("Jet2Pt_%s",processLabels[i].c_str()), ";Second Jet p_{T} [GeV/c];Number of Events", 125, 0, 500, color[i],  (hasData && i==0) ));  
    // histDPhiHiggsMET.push_back( makeNewHist( Form("DPhiHiggsMET_%s",processLabels[i].c_str()), ";#Delta#phi(#gamma#gamma,MET)];Number of Events", 20, 0, 3.15, color[i],  (hasData && i==0) ));  
    // histPtggPeakRegion.push_back( makeNewHist( Form("PtggPeakRegion_%s",processLabels[i].c_str()), ";p_{T #gamma#gamma} [GeV/c];Number of Events", 20, 0, 500, color[i],  (hasData && i==0) ));  
    // histPtggSidebandRegion.push_back( makeNewHist( Form("PtggSidebandRegion_%s",processLabels[i].c_str()), ";p_{T #gamma#gamma} [GeV/c];Number of Events", 20, 0, 500, color[i],  (hasData && i==0) ));  
    // histDPhiRazor.push_back( makeNewHist( Form("DPhiRazor_%s",processLabels[i].c_str()), ";#Delta#phi Hemispheres ;Number of Events", 50, 0, 3.14, color[i],  (hasData && i==0) ));  

    // histCosThetaStar.push_back( makeNewHist( Form("CosThetaStar_%s",processLabels[i].c_str()), ";Cos(#theta^{*}_{CS}) ;Number of Events", 50, 0, 1, color[i],  (hasData && i==0) ));  
  }
  //*******************************************************************************************
  //Define Counts
  //*******************************************************************************************



  //*******************************************************************************************
  //Read files
  //*******************************************************************************************
  for (uint i=0; i < inputfiles.size(); ++i) {

    TFile* inputFile = new TFile(inputfiles[i].c_str(),"READ");
    assert(inputFile);
    TTree* tree = 0;
    tree = (TTree*)inputFile->Get("HZZRazor");    
    assert(tree);

    float w = 0;
    int njets = 0;
    float jetE[40]; //[njets]
    float jetpt[40]; //[njets]
    float jeteta[40]; //[njets]
    float jetphi[40]; //[njets]
    float lep1pt = 0;
    float lep1eta = 0;
    float lep1phi = 0;
    float lep2pt = 0;
    float lep2eta = 0;
    float lep2phi = 0;
    float lep3pt = 0;
    float lep3eta = 0;
    float lep3phi = 0;
    float lep4pt = 0;
    float lep4eta = 0;
    float lep4phi = 0;
    int nBTaggedJets = 0;
    float dPhiRazor = 0;
    float MR = 0;
    float m4l = 0;
    float pt4l = 0;
    float eta4l = 0;
    float phi4l = 0;
    float metphi = 0;
    float met = 0;
    float Rsq = 0;
    uint run = 0;
    uint lumi = 0;
    uint event = 0;
    uint npu = 0;

    tree->SetBranchAddress("run",&run);
    tree->SetBranchAddress("lumi",&lumi);
    tree->SetBranchAddress("event",&event);
    tree->SetBranchAddress("npu",&npu);
    tree->SetBranchAddress("weight",&w);
    tree->SetBranchAddress("n_Jets",&njets);
    tree->SetBranchAddress("jet_E",&jetE);
    tree->SetBranchAddress("jet_Pt",&jetpt);
    tree->SetBranchAddress("jet_Eta",&jeteta);
    tree->SetBranchAddress("jet_Phi",&jetphi);
    tree->SetBranchAddress("lep1Pt",&lep1pt);
    tree->SetBranchAddress("lep1Eta",&lep1eta);
    tree->SetBranchAddress("lep1Phi",&lep1phi);
    tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
    tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
    tree->SetBranchAddress("MR",&MR);
    tree->SetBranchAddress("m4l",&m4l);
    tree->SetBranchAddress("pt4l",&pt4l);
    tree->SetBranchAddress("eta4l",&eta4l);
    tree->SetBranchAddress("phi4l",&phi4l);
    tree->SetBranchAddress("t1Rsq",&Rsq);
    tree->SetBranchAddress("t1MET",&met);
    tree->SetBranchAddress("t1METPhi",&metphi);

    cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";

    for (int n=0;n<tree->GetEntries();n++) { 
    
      double weight = 1.0;
      double puWeight = 1.0;
      puWeight = pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(npu));
 
      if (!(hasData && i==0)) weight = w*puWeight*intLumi;
      //if (i >= 3 && i<=6) weight *= 2;
      //if (!(hasData && i==0)) weight *= 1.6;

      //cout << "weight: " << weight << "\n";

      tree->GetEntry(n);
      if (n % 1000000 == 0) cout << "Processing Event " << n << "\n";       

      //Options
      if (option == 0 ) {
      }

      //apply selection cuts
      if (!(m4l > 100 && m4l < 160)) continue;
      //if (!(MR> 0 && Rsq > 0.0 )) continue;
      //if (!(MR> 0 && Rsq > 0.035 )) continue; 
      if (!(MR> 0 && Rsq > 0.1 )) continue; 


      //count jets
      int njets80 = 0;
      int njets60 = 0;
      int njets40 = 0;
      double leadjetpt = 0;
      double secondjetpt = 0;

      for(int j=0; j < njets; ++j) {

	if (fabs(jeteta[j]) >=  3.0) continue;

	if (jetpt[j] > leadjetpt) {
	  secondjetpt = leadjetpt;
	  leadjetpt = jetpt[j];

	  TLorentzVector v1; v1.SetPtEtaPhiE(jetpt[j],jeteta[j],jetphi[j],jetE[j]);
	  TLorentzVector v2; v2.SetPtEtaPhiM(pt4l,eta4l,phi4l,m4l);
	} else if (jetpt[j] > secondjetpt) {
	  secondjetpt = jetpt[j];
	}

      	if (jetpt[j] > 80 && fabs(jeteta[j]) < 3.0) {
      	  njets80++;
      	}
      	if (jetpt[j] > 60 && fabs(jeteta[j]) < 3.0) {
      	  njets60++;
      	}
      	if (jetpt[j] > 40 && fabs(jeteta[j]) < 3.0) {
      	  njets40++;
      	}
      }

      
      //cout << "lead jets: " << leadjetpt << " " << secondjetpt << "\n";
      if (i==0) {
	cout << run << " " << lumi << " " <<  event << " : " << m4l << " " << " " << MR << " " << Rsq << " \n";
      }

      histM4l[i]->Fill(m4l, weight);

      if (m4l > 120 && m4l < 130) {
	EventCount[i] += weight;
	EventCountErrSqr[i] += weight*weight;	
	histMR[i]->Fill(MR,weight);
	if (Rsq != Rsq) Rsq = 0;
	histRsq[i]->Fill(Rsq,weight);
      }
  
    }

    inputFile->Close();
    delete inputFile;
  
  }
  

  //*******************************************************************************************
  //Draw Plots
  //*******************************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;
  TLatex *tex = 0;



  //*******************************************************************************************
  //M4l
  //*******************************************************************************************
  PlotDataAndStackedBkg( histM4l , processLabels, color, true, "M4l", Label);
  PlotDataAndStackedBkg( histMR , processLabels, color, true, "MR", Label);
  PlotDataAndStackedBkg( histRsq , processLabels, color, true, "Rsq", Label);

 
  //*******************************************************************************************
  //Summarize Counts
  //*******************************************************************************************
  for(int i=0; i<int(processLabels.size()); i++) {
    cout << processLabels[i] << " : " << EventCount[i] << " +/- " << sqrt(EventCountErrSqr[i]) << "\n";
  }
  
   //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("RazorPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histM4l[i], Form("histM4l_%s",processLabels[i].c_str()), "WriteDelete");
    // file->WriteTObject(histNJets80[i], Form("histNJets80_%s",processLabels[i].c_str()), "WriteDelete");
    // file->WriteTObject(histNJets60[i], Form("histNJets60_%s",processLabels[i].c_str()), "WriteDelete");
    // file->WriteTObject(histNJets40[i], Form("histNJets40_%s",processLabels[i].c_str()), "WriteDelete");
    // file->WriteTObject(histNJets30[i], Form("histNJets30_%s",processLabels[i].c_str()), "WriteDelete");
    // file->WriteTObject(histPhoton1Pt[i], Form("histPhoton1Pt_%s",processLabels[i].c_str()), "WriteDelete");
    // file->WriteTObject(histPhoton2Pt[i], Form("histPhoton2Pt_%s",processLabels[i].c_str()), "WriteDelete");
    // file->WriteTObject(histJet1Pt[i], Form("histJet1Pt_%s",processLabels[i].c_str()), "WriteDelete");
    // file->WriteTObject(histJet2Pt[i], Form("histJet2Pt_%s",processLabels[i].c_str()), "WriteDelete");
    // file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    // file->WriteTObject(histDPhiRazor[i], Form("histDPhiRazor_%s",processLabels[i].c_str()), "WriteDelete");
    // file->WriteTObject(histPtggPeakRegion[i], Form("histPtggPeakRegion_%s",processLabels[i].c_str()), "WriteDelete");
    // file->WriteTObject(histPtggSidebandRegion[i], Form("histPtggSidebandRegion_%s",processLabels[i].c_str()), "WriteDelete");
  }
  
  // file->WriteTObject(stackMR, "stackMR", "WriteDelete");
  // file->WriteTObject(stackRsq, "stackRsq", "WriteDelete");  
  // file->WriteTObject(stackDPhiRazor, "stackDPhiRazor", "WriteDelete");  

 }


 void MakeHZZRazorPlots() {

   string datafile = "";
   string dataLabel = "";

   // string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/MiniIso/RazorInclusive_SMS-T1qqqq_2J_mGl-1400_mLSP-100_1pb_weighted.root";  
   // string dataLabel = "T1qqqq m_{G}=1400 m_{LSP}=100";
   // string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/MiniIso/RazorInclusive_SMS-T1bbbb_2J_mGl-1500_mLSP-100_1pb_weighted.root";  
   // string dataLabel = "T1bbbb m_{G}=1500 m_{LSP}=100";
   //string datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/MiniIso/RazorInclusive_SMS-T1tttt_2J_mGl-1500_mLSP-100_1pb_weighted.root";  
   //string dataLabel = "T1tttt m_{G}=1500 m_{LSP}=100";   

   dataLabel = "data";
   datafile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorHZZ/HZZRazor_Data_DoubleLeptonTriggers_Run2012_GoodLumi_DuplicateRemoved.root";

   vector<string> bkgfiles;
   vector<string> bkgLabels;
   vector<int> bkgColors;

   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorHZZ/HZZRazor_ZZ_1pb_weighted.root");
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorHZZ/HZZRazor_Top_1pb_weighted.root");
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorHZZ/HZZRazor_DYJetsToLL_HTBinned_1pb_weighted.root");
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorHZZ/HZZRazor_GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6_1pb_weighted.root");
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorHZZ/HZZRazor_VBF_HToZZTo4L_M-125_8TeV-powheg-pythia6_1pb_weighted.root");
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorHZZ/HZZRazor_VH_1pb_weighted.root");
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorHZZ/HZZRazor_TTbarH_HToZZTo4L_M-125_8TeV-pythia6_1pb_weighted.root");
   

   bkgLabels.push_back("ZZ");
   bkgLabels.push_back("t#bar{t}+X");
   bkgLabels.push_back("Z+Jets");
   bkgLabels.push_back("GluonFusion Higgs");
   bkgLabels.push_back("VBF Higgs");
   bkgLabels.push_back("W/Z + Higgs");
   bkgLabels.push_back("t#bar{t} Higgs");

   bkgColors.push_back(kRed);
   bkgColors.push_back(kGreen+2);
   bkgColors.push_back(kBlue);
   bkgColors.push_back(kAzure+10);
   bkgColors.push_back(kGray);
   bkgColors.push_back(kOrange+1);  
   bkgColors.push_back(kViolet);  

   //RunMakeRazorPlots(datafile,dataLabel,bkgfiles,bkgLabels,bkgColors,0,1,"MR250Rsq0p035","H#rightarrowZZ Razor");
   //RunMakeRazorPlots(datafile,dataLabel,bkgfiles,bkgLabels,bkgColors,0,1,"MR250Rsq0p035_EverythingTimes1p6","H#rightarrowZZ Razor");
   
   //RunMakeRazorPlots(datafile,dataLabel,bkgfiles,bkgLabels,bkgColors,0,1,"HZZInclusive","H#rightarrowZZ Razor");
   //RunMakeRazorPlots(datafile,dataLabel,bkgfiles,bkgLabels,bkgColors,0,1,"MR0Rsq0","H#rightarrowZZ Razor");
   //RunMakeRazorPlots(datafile,dataLabel,bkgfiles,bkgLabels,bkgColors,0,1,"MR0Rsq0p035","H#rightarrowZZ Razor");
   RunMakeRazorPlots(datafile,dataLabel,bkgfiles,bkgLabels,bkgColors,0,1,"MR250","H#rightarrowZZ Razor");
   
 
 }
 
