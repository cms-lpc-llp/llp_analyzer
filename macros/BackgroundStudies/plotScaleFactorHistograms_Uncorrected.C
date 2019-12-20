
#include "macros/tdrstyle.C"
#include "macros/CMS_lumi.C"

void plotScaleFactor() {

  //TFile *inf = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_CorrectedToMultiJet.root","READ");
  TFile *inf = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_Uncorrected.root","READ");
  // TFile *inf = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_GJetsInv.root","READ");

  TH2Poly *ttbarNominal = (TH2Poly*)inf->Get("TTJetsScaleFactors");
  TH2Poly *ttbarUp = (TH2Poly*)inf->Get("TTJetsScaleFactorsUp");
  TH2Poly *ttbarDown = (TH2Poly*)inf->Get("TTJetsScaleFactorsDown");
  TH2Poly *wNominal = (TH2Poly*)inf->Get("WJetsScaleFactors");
  TH2Poly *wUp = (TH2Poly*)inf->Get("WJetsScaleFactorsUp");
  TH2Poly *wDown = (TH2Poly*)inf->Get("WJetsScaleFactorsDown");
  TH2Poly *wInvNominal = (TH2Poly*)inf->Get("WJetsInvScaleFactors");
  TH2Poly *wInvUp = (TH2Poly*)inf->Get("WJetsInvScaleFactorsUp");
  TH2Poly *wInvDown = (TH2Poly*)inf->Get("WJetsInvScaleFactorsDown");
  TH2Poly *GJetInvNominal = (TH2Poly*)inf->Get("GJetsInvScaleFactors");


  TCanvas *cv = 0;
  gStyle->SetPaintTextFormat("4.2f");

  //****************************************************
  //Plot GJetsInv Scale Factors
  //****************************************************
  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(53);
  GJetInvNominal->Draw("colztexte1");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  GJetInvNominal->GetXaxis()->SetRangeUser(400,4000);
  GJetInvNominal->GetYaxis()->SetRangeUser(0.25,1.5);
  GJetInvNominal->GetZaxis()->SetTitle("Data to MC Correction Factor");
  GJetInvNominal->GetXaxis()->SetTitle("M_{R} [GeV/c^{2}]");
  GJetInvNominal->GetYaxis()->SetTitle("R^{2}");
  GJetInvNominal->SetTitle("");
  GJetInvNominal->GetZaxis()->SetLabelSize(0.05);
  GJetInvNominal->GetZaxis()->SetTitleSize(0.05);
  GJetInvNominal->GetXaxis()->SetLabelSize(0.05);
  GJetInvNominal->GetXaxis()->SetTitleSize(0.05);
  GJetInvNominal->GetXaxis()->SetTitleOffset(0.8);
  GJetInvNominal->GetYaxis()->SetLabelSize(0.05);
  GJetInvNominal->GetYaxis()->SetTitleSize(0.05);
  GJetInvNominal->GetYaxis()->SetTitleOffset(0.8);
  GJetInvNominal->SetStats(false);
  GJetInvNominal->SetMaximum(1.8);
  GJetInvNominal->SetMinimum(0.35);


  lumi_13TeV = "2.3 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("GJetsInvScaleFactor_Uncorrected.png");
  cv->SaveAs("GJetsInvScaleFactor_Uncorrected.pdf");

  TH2Poly *GJetInvUncertainties = (TH2Poly*)GJetInvNominal->Clone("GJetInvUncertainties");
  for (int i=1; i<GJetInvUncertainties->GetNumberOfBins()+1; ++i) {
    GJetInvUncertainties->SetBinContent(i,100*GJetInvNominal->GetBinError(i) / GJetInvNominal->GetBinContent(i));
    cout << i << " : " << GJetInvNominal->GetBinError(i) << " " << GJetInvNominal->GetBinContent(i) << " : " << GJetInvNominal->GetBinError(i) / GJetInvNominal->GetBinContent(i) << "\n";
  }

  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(1);
  gStyle->SetPaintTextFormat("4.2f");
  GJetInvUncertainties->Draw("colztext");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  GJetInvUncertainties->SetMarkerSize(2.0);
  GJetInvUncertainties->SetTitle("");
  GJetInvUncertainties->GetXaxis()->SetTitle("M_{R} [GeV/c^{2}]");
  GJetInvUncertainties->GetYaxis()->SetTitle("R^{2}");
  GJetInvUncertainties->GetXaxis()->SetRangeUser(400,4000);
  GJetInvUncertainties->GetYaxis()->SetRangeUser(0.25,1.5);
  GJetInvUncertainties->GetZaxis()->SetTitle("Systematic Uncertainty (%)");
  GJetInvUncertainties->GetZaxis()->SetLabelSize(0.05);
  GJetInvUncertainties->GetZaxis()->SetTitleSize(0.05);
  GJetInvUncertainties->GetXaxis()->SetLabelSize(0.05);
  GJetInvUncertainties->GetXaxis()->SetTitleSize(0.05);
  GJetInvUncertainties->GetXaxis()->SetTitleOffset(0.8);
  GJetInvUncertainties->GetYaxis()->SetLabelSize(0.05);
  GJetInvUncertainties->GetYaxis()->SetTitleSize(0.05);
  GJetInvUncertainties->GetYaxis()->SetTitleOffset(0.8);
  GJetInvUncertainties->SetStats(false);
  GJetInvUncertainties->SetMaximum(50);
  GJetInvUncertainties->SetMinimum(0.0);
  lumi_13TeV = "2.3 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("GJetsInvScaleFactorUncertainty.png");
  cv->SaveAs("GJetsInvScaleFactorUncertainty.pdf");



  //****************************************************
  //Plot WJetsInv Scale Factors
  //****************************************************
  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(53);
  wInvNominal->Draw("colztexte1");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  wInvNominal->GetXaxis()->SetRangeUser(300,4000);
  wInvNominal->GetYaxis()->SetRangeUser(0.15,1.5);
  wInvNominal->GetZaxis()->SetTitle("Data to MC Correction Factor");
  wInvNominal->GetZaxis()->SetLabelSize(0.05);
  wInvNominal->GetZaxis()->SetTitleSize(0.05);
  wInvNominal->GetXaxis()->SetLabelSize(0.05);
  wInvNominal->GetXaxis()->SetTitleSize(0.05);
  wInvNominal->GetXaxis()->SetTitleOffset(0.8);
  wInvNominal->GetYaxis()->SetLabelSize(0.05);
  wInvNominal->GetYaxis()->SetTitleSize(0.05);
  wInvNominal->GetYaxis()->SetTitleOffset(0.8);
  wInvNominal->SetStats(false);
  wInvNominal->SetMaximum(1.8);
  wInvNominal->SetMinimum(0.0);


  lumi_13TeV = "2.3 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("WJetsInvScaleFactor_Uncorrected.png");
  cv->SaveAs("WJetsInvScaleFactor_Uncorrected.pdf");




 //****************************************************
  //Plot WJets Scale Factors
  //****************************************************
  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(53);
  wNominal->Draw("colztexte1");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  wNominal->GetXaxis()->SetRangeUser(300,4000);
  wNominal->GetYaxis()->SetRangeUser(0.15,1.5);
  wNominal->GetZaxis()->SetTitle("Data to MC Correction Factor");
  wNominal->GetZaxis()->SetLabelSize(0.05);
  wNominal->GetZaxis()->SetTitleSize(0.05);
  wNominal->GetXaxis()->SetLabelSize(0.05);
  wNominal->GetXaxis()->SetTitleSize(0.05);
  wNominal->GetXaxis()->SetTitleOffset(0.8);
  wNominal->GetYaxis()->SetLabelSize(0.05);
  wNominal->GetYaxis()->SetTitleSize(0.05);
  wNominal->GetYaxis()->SetTitleOffset(0.8);
  wNominal->SetStats(false);
  wNominal->SetMaximum(1.8);
  wNominal->SetMinimum(0.35);

  lumi_13TeV = "2.3 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("WJetsScaleFactor_Uncorrected.png");
  cv->SaveAs("WJetsScaleFactor_Uncorrected.pdf");




  TH2Poly *WJetsUncertainties = (TH2Poly*)wNominal->Clone("WJetsUncertainties");
  for (int i=1; i<WJetsUncertainties->GetNumberOfBins()+1; ++i) {
    WJetsUncertainties->SetBinContent(i,100*wNominal->GetBinError(i) / wNominal->GetBinContent(i));
  }

  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(1);
  gStyle->SetPaintTextFormat("4.2f");
  WJetsUncertainties->Draw("colztext");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  WJetsUncertainties->SetMarkerSize(2.0);
  WJetsUncertainties->SetTitle("");
  WJetsUncertainties->GetXaxis()->SetTitle("M_{R} [GeV/c^{2}]");
  WJetsUncertainties->GetYaxis()->SetTitle("R^{2}");
  WJetsUncertainties->GetXaxis()->SetRangeUser(300,4000);
  WJetsUncertainties->GetYaxis()->SetRangeUser(0.15,1.5);
  WJetsUncertainties->GetZaxis()->SetTitle("Systematic Uncertainty (%)");
  WJetsUncertainties->GetZaxis()->SetLabelSize(0.05);
  WJetsUncertainties->GetZaxis()->SetTitleSize(0.05);
  WJetsUncertainties->GetXaxis()->SetLabelSize(0.05);
  WJetsUncertainties->GetXaxis()->SetTitleSize(0.05);
  WJetsUncertainties->GetXaxis()->SetTitleOffset(0.8);
  WJetsUncertainties->GetYaxis()->SetLabelSize(0.05);
  WJetsUncertainties->GetYaxis()->SetTitleSize(0.05);
  WJetsUncertainties->GetYaxis()->SetTitleOffset(0.8);
  WJetsUncertainties->SetStats(false);
  WJetsUncertainties->SetMaximum(50);
  WJetsUncertainties->SetMinimum(0.0);
  lumi_13TeV = "2.3 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("WJetsScaleFactorUncertainty.png");
  cv->SaveAs("WJetsScaleFactorUncertainty.pdf");




 //****************************************************
  //Plot TTBar Scale Factors
  //****************************************************
  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(53);
  ttbarNominal->Draw("colztexte1");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  ttbarNominal->GetXaxis()->SetRangeUser(300,4000);
  ttbarNominal->GetYaxis()->SetRangeUser(0.15,1.5);
  ttbarNominal->GetZaxis()->SetTitle("Data to MC Correction Factor");
  ttbarNominal->GetZaxis()->SetLabelSize(0.05);
  ttbarNominal->GetZaxis()->SetTitleSize(0.05);
  ttbarNominal->GetXaxis()->SetLabelSize(0.05);
  ttbarNominal->GetXaxis()->SetTitleSize(0.05);
  ttbarNominal->GetXaxis()->SetTitleOffset(0.8);
  ttbarNominal->GetYaxis()->SetLabelSize(0.05);
  ttbarNominal->GetYaxis()->SetTitleSize(0.05);
  ttbarNominal->GetYaxis()->SetTitleOffset(0.8);
  ttbarNominal->SetStats(false);
  // ttbarNominal->SetMaximum(10000);
  ttbarNominal->SetMaximum(1.8);
  ttbarNominal->SetMinimum(0.);

  lumi_13TeV = "2.3 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("TTBarScaleFactor_Uncorrected.png");
  cv->SaveAs("TTBarScaleFactor_Uncorrected.pdf");





  TH2Poly *TTBarUncertainties = (TH2Poly*)ttbarNominal->Clone("TTBarUncertainties");
  for (int i=1; i<TTBarUncertainties->GetNumberOfBins()+1; ++i) {
    TTBarUncertainties->SetBinContent(i,100*ttbarNominal->GetBinError(i) / ttbarNominal->GetBinContent(i));
  }

  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(1);
  gStyle->SetPaintTextFormat("4.2f");
  TTBarUncertainties->Draw("colztext");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  TTBarUncertainties->SetMarkerSize(2.0);
  TTBarUncertainties->SetTitle("");
  TTBarUncertainties->GetXaxis()->SetTitle("M_{R} [GeV/c^{2}]");
  TTBarUncertainties->GetYaxis()->SetTitle("R^{2}");
  TTBarUncertainties->GetXaxis()->SetRangeUser(300,4000);
  TTBarUncertainties->GetYaxis()->SetRangeUser(0.15,1.5);
  TTBarUncertainties->GetZaxis()->SetTitle("Systematic Uncertainty (%)");
  TTBarUncertainties->GetZaxis()->SetLabelSize(0.05);
  TTBarUncertainties->GetZaxis()->SetTitleSize(0.05);
  TTBarUncertainties->GetXaxis()->SetLabelSize(0.05);
  TTBarUncertainties->GetXaxis()->SetTitleSize(0.05);
  TTBarUncertainties->GetXaxis()->SetTitleOffset(0.8);
  TTBarUncertainties->GetYaxis()->SetLabelSize(0.05);
  TTBarUncertainties->GetYaxis()->SetTitleSize(0.05);
  TTBarUncertainties->GetYaxis()->SetTitleOffset(0.8);
  TTBarUncertainties->SetStats(false);
  TTBarUncertainties->SetMaximum(50);
  TTBarUncertainties->SetMinimum(0.0);
  lumi_13TeV = "2.3 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("TTBarScaleFactorUncertainty.png");
  cv->SaveAs("TTBarScaleFactorUncertainty.pdf");



  


}




void plotGJetsScaleFactorSystematics() {

  TFile *inf = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_CorrectedToMultiJet.root","READ");

  TH2Poly *ttbarNominal = (TH2Poly*)inf->Get("TTJetsScaleFactors");
  TH2Poly *ttbarUp = (TH2Poly*)inf->Get("TTJetsScaleFactorsUp");
  TH2Poly *ttbarDown = (TH2Poly*)inf->Get("TTJetsScaleFactorsDown");
  TH2Poly *wNominal = (TH2Poly*)inf->Get("WJetsScaleFactors");
  TH2Poly *wUp = (TH2Poly*)inf->Get("WJetsScaleFactorsUp");
  TH2Poly *wDown = (TH2Poly*)inf->Get("WJetsScaleFactorsDown");
  TH2Poly *wInvNominal = (TH2Poly*)inf->Get("WJetsInvScaleFactors");
  TH2Poly *wInvUp = (TH2Poly*)inf->Get("WJetsInvScaleFactorsUp");
  TH2Poly *wInvDown = (TH2Poly*)inf->Get("WJetsInvScaleFactorsDown");
  TH2Poly *GJetInvNominal = (TH2Poly*)inf->Get("GJetsInvScaleFactors");


  TCanvas *cv = 0;
  gStyle->SetPaintTextFormat("4.2f");

  //****************************************************
  //Systematic Uncertainty and GJets Down SF Histogram
  //****************************************************
  TH2Poly *GJetsSystematicUnc = (TH2Poly*)GJetInvNominal->Clone("GJetsSystematicUnc");
  TH2Poly *GJetsScaleFactor_Down = (TH2Poly*)GJetInvNominal->Clone("GJetsInvScaleFactors_Down");

  //Get bins of each histogram
  TList *wInvBins = wInvNominal->GetBins();
  TList *gInvBins = GJetInvNominal->GetBins();

  //Loop over GJets bins
  TH2PolyBin *gBin, *wBin; //temp variables to hold bin info
  for (int i = 1; i < GJetsSystematicUnc->GetNumberOfBins()+1; ++i) {

      //Get GJets bin
      gBin = (TH2PolyBin*)gInvBins->At(i-1);

      cout << "In bin " << i << " of GJets histogram" << endl;
      //cout << gBin->GetXMin() << " " << gBin->GetXMax() << " " << gBin->GetYMin() << " " << gBin->GetYMax() << endl;

      //Find out which WJets bin we are in
      int wBinNum = -1;
      for (int j = 1; j < wInvNominal->GetNumberOfBins()+1; ++j) {

          //Get WJets bin
          wBin = (TH2PolyBin*)wInvBins->At(j-1);

          //cout << "In bin " << j << " of WJets histogram" << endl;
          //cout << wBin->GetXMin() << " " << wBin->GetXMax() << " " << wBin->GetYMin() << " " << wBin->GetYMax() << endl;

          //Check if this GJets bin is inside this WJets bin
          if ( gBin->GetXMin() >= wBin->GetXMin() &&
               gBin->GetXMax() <= wBin->GetXMax() &&
               gBin->GetYMin() >= wBin->GetYMin() &&
               gBin->GetYMax() <= wBin->GetYMax() ) {
              cout << "This GJets bin is inside bin " << j << " of WJets histogram" << endl;
              wBinNum = j;
              break;
          }
      }

      double gjet = GJetInvNominal->GetBinContent(i);
      double wjet = wInvNominal->GetBinContent(wBinNum);

      //Set bin content of each histogram
      GJetsSystematicUnc->SetBinContent(i, fabs(gjet - wjet)/gjet );
      GJetsScaleFactor_Down->SetBinContent(i, gjet - (wjet - gjet) );

      cout << "Bin " << i << " : " << gjet << " , " << wjet << " , " <<  gjet - (wjet - gjet) << "\n";
  }

  cv = new TCanvas("cv","cv", 800,600);
  gStyle->SetPalette(1);
  GJetsSystematicUnc->Draw("colztext");
  cv->SetLogx();
  cv->SetLogy();
  cv->SetRightMargin(0.175);
  cv->SetBottomMargin(0.12);
  GJetsSystematicUnc->GetXaxis()->SetRangeUser(400,4000);
  GJetsSystematicUnc->GetYaxis()->SetRangeUser(0.25,1.5);
  GJetsSystematicUnc->GetZaxis()->SetTitle("Systematic Uncertainty");
  GJetsSystematicUnc->GetXaxis()->SetTitle("M_{R} [GeV/c^{2}]");
  GJetsSystematicUnc->GetYaxis()->SetTitle("R^{2}");
  GJetsSystematicUnc->SetTitle("");
  GJetsSystematicUnc->GetZaxis()->SetLabelSize(0.05);
  GJetsSystematicUnc->GetZaxis()->SetTitleSize(0.05);
  GJetsSystematicUnc->GetXaxis()->SetLabelSize(0.05);
  GJetsSystematicUnc->GetXaxis()->SetTitleSize(0.05);
  GJetsSystematicUnc->GetXaxis()->SetTitleOffset(0.8);
  GJetsSystematicUnc->GetYaxis()->SetLabelSize(0.05);
  GJetsSystematicUnc->GetYaxis()->SetTitleSize(0.05);
  GJetsSystematicUnc->GetYaxis()->SetTitleOffset(0.8);
  GJetsSystematicUnc->SetStats(false);
  GJetsSystematicUnc->SetMaximum(1.0);
  GJetsSystematicUnc->SetMinimum(0.0);

  lumi_13TeV = "2.3 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  lumiTextSize = 0.5;
  cmsTextSize = 0.6;
  extraOverCmsTextSize = 0.85;
  CMS_lumi(cv,4,0);
  cv->SaveAs("GJetsVsWJetsSystematic.png");
  cv->SaveAs("GJetsVsWJetsSystematic.pdf");

  TFile *outf = new TFile("RazorScaleFactors_Inclusive_CorrectedToMultiJet.root","UPDATE");
  outf->WriteTObject(ttbarNominal);
  outf->WriteTObject(ttbarUp);
  outf->WriteTObject(ttbarDown);
  outf->WriteTObject(wNominal);
  outf->WriteTObject(wUp);
  outf->WriteTObject(wDown);
  outf->WriteTObject(wInvNominal);
  outf->WriteTObject(wInvUp);
  outf->WriteTObject(wInvDown);
  outf->WriteTObject(GJetInvNominal);
  outf->WriteTObject(GJetsScaleFactor_Down);
  outf->Close();

}



void plotScaleFactorHistograms_Uncorrected() {
  gROOT->SetBatch();
  plotScaleFactor();
  //plotGJetsScaleFactorSystematics();

}
