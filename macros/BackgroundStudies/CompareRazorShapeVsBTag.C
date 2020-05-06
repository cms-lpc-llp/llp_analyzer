
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <vector>
#include <map>
#include <iostream>
#include "RazorAnalyzer/macros/tdrstyle.C"
#include "RazorAnalyzer/macros/CMS_lumi.C"

const Int_t NComponents = 10;
int color[NComponents] = {kRed, kGreen+2, kBlue, kViolet, kAzure+10, kBlack, kOrange+1, kGray, kBlack, kBlack};


//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
TH1D* NormalizeHist(TH1D *originalHist) {
  TH1D* hist = (TH1D*)originalHist->Clone((string(originalHist->GetName())+"_normalized").c_str());
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

void AddLargestWeightUncertainty(TH1D* hist) {
 
  for (int i=0; i<hist->GetXaxis()->GetNbins()+2;++i) {
    //hist->SetBinError(i, sqrt( pow(hist->GetBinError(i),2) + pow( 0.0054 , 2)));
    hist->SetBinError(i, sqrt( pow(hist->GetBinError(i),2) + pow( 0.01 , 2)));
  }

}


// TH1D* MakeCorrectedTwoBTagShape(TH1D* hist, string category) {
 
//   TH1D* result = (TH1D*)hist->Clone( (hist->GetName()+"Corr").c_str());

//   for (int i=0; i<result->GetXaxis()->GetNbins()+2;++i) {
//     double slope = 1;
//     if (category == "MultiJet") slope = 0.00053133;
//     if (category == "MuonMultiJet") slope = 0.00085458;
//     if (category == "ElectronMultiJet") slope = 0.00078621;

//     double value = corr * result->GetBinContent(i);
//     double err = corr * result->GetBinError(i);
//     result->SetBinContent(i, value);
//     result->SetBinError(i, err);
//   }

// }





void PlotTwoBTagVsThreeBTagShapes(string binLabel = "", bool useCorrected = false) {

  string BinLabel = "";
  if (binLabel != "") BinLabel = "_" + binLabel;

  TFile *fileMultijet_TwoBTag = new TFile(Form("RazorPlots_Multijet_TwoBTags%s.root",BinLabel.c_str()),"READ");
  TFile *fileMultijet_ThreeBTag = new TFile(Form("RazorPlots_Multijet_ThreeBTags%s.root",BinLabel.c_str()),"READ");
  // TFile *fileLooseLeptonMultijet_TwoBTag = new TFile(Form("RazorPlots_LooseLeptonMultijet_TwoBTags%s.root",BinLabel.c_str()),"READ");
  // TFile *fileLooseLeptonMultijet_ThreeBTag = new TFile(Form("RazorPlots_LooseLeptonMultijet_ThreeBTags%s.root",BinLabel.c_str()),"READ");
  TFile *fileMuonMultijet_TwoBTag = new TFile(Form("RazorPlots_MuonMultijet_TwoBTags%s.root",BinLabel.c_str()),"READ");
  TFile *fileMuonMultijet_ThreeBTag = new TFile(Form("RazorPlots_MuonMultijet_ThreeBTags%s.root",BinLabel.c_str()),"READ");
  TFile *fileElectronMultijet_TwoBTag = new TFile(Form("RazorPlots_ElectronMultijet_TwoBTags%s.root",BinLabel.c_str()),"READ");
  TFile *fileElectronMultijet_ThreeBTag = new TFile(Form("RazorPlots_ElectronMultijet_ThreeBTags%s.root",BinLabel.c_str()),"READ");

  TH1D *histMR_Multijet_TwoBTag = (TH1D*)fileMultijet_TwoBTag->Get("histMRAllBkg");
  TH1D *histRsq_Multijet_TwoBTag = (TH1D*)fileMultijet_TwoBTag->Get("histRsqAllBkg");
  TH1D *histMR_Multijet_ThreeBTag = (TH1D*)fileMultijet_ThreeBTag->Get("histMRAllBkg");
  TH1D *histRsq_Multijet_ThreeBTag = (TH1D*)fileMultijet_ThreeBTag->Get("histRsqAllBkg");
  // TH1D *histMR_LooseLeptonMultijet_TwoBTag = (TH1D*)fileLooseLeptonMultijet_TwoBTag->Get("histMRAllBkg");
  // TH1D *histRsq_LooseLeptonMultijet_TwoBTag = (TH1D*)fileLooseLeptonMultijet_TwoBTag->Get("histRsqAllBkg");
  // TH1D *histMR_LooseLeptonMultijet_ThreeBTag = (TH1D*)fileLooseLeptonMultijet_ThreeBTag->Get("histMRAllBkg");
  // TH1D *histRsq_LooseLeptonMultijet_ThreeBTag = (TH1D*)fileLooseLeptonMultijet_ThreeBTag->Get("histRsqAllBkg");
  TH1D *histMR_MuonMultijet_TwoBTag = (TH1D*)fileMuonMultijet_TwoBTag->Get("histMRAllBkg");
  TH1D *histRsq_MuonMultijet_TwoBTag = (TH1D*)fileMuonMultijet_TwoBTag->Get("histRsqAllBkg");
  TH1D *histMR_MuonMultijet_ThreeBTag = (TH1D*)fileMuonMultijet_ThreeBTag->Get("histMRAllBkg");
  TH1D *histRsq_MuonMultijet_ThreeBTag = (TH1D*)fileMuonMultijet_ThreeBTag->Get("histRsqAllBkg");
  TH1D *histMR_ElectronMultijet_TwoBTag = (TH1D*)fileElectronMultijet_TwoBTag->Get("histMRAllBkg");
  TH1D *histRsq_ElectronMultijet_TwoBTag = (TH1D*)fileElectronMultijet_TwoBTag->Get("histRsqAllBkg");
  TH1D *histMR_ElectronMultijet_ThreeBTag = (TH1D*)fileElectronMultijet_ThreeBTag->Get("histMRAllBkg");
  TH1D *histRsq_ElectronMultijet_ThreeBTag = (TH1D*)fileElectronMultijet_ThreeBTag->Get("histRsqAllBkg");

  // //Make 2 b-tag shape plus linear correction
  TH1D *histMR_Multijet_TwoBTagCorr = (TH1D*)fileMultijet_TwoBTag->Get("histMRAllBkgCorr");
  TH1D *histMR_MuonMultijet_TwoBTagCorr = (TH1D*)fileMuonMultijet_TwoBTag->Get("histMRAllBkgCorr");
  TH1D *histMR_ElectronMultijet_TwoBTagCorr = (TH1D*)fileElectronMultijet_TwoBTag->Get("histMRAllBkgCorr");


  TCanvas *cv = 0;
  TLegend *legend = 0;
  TPad *pad1 = 0;
  TH1D *ratioHist = 0;
  TH1D *ratioHistCorr = 0;
  TPad *pad2 = 0;
  double chisq = 0;
  int ndof = 0;

  AddLargestWeightUncertainty(histMR_Multijet_ThreeBTag);
  AddLargestWeightUncertainty(histRsq_Multijet_ThreeBTag);
  // AddLargestWeightUncertainty(histMR_LooseLeptonMultijet_ThreeBTag);
  // AddLargestWeightUncertainty(histRsq_LooseLeptonMultijet_ThreeBTag);
  AddLargestWeightUncertainty(histMR_MuonMultijet_ThreeBTag);
  AddLargestWeightUncertainty(histRsq_MuonMultijet_ThreeBTag);
  AddLargestWeightUncertainty(histMR_ElectronMultijet_ThreeBTag);
  AddLargestWeightUncertainty(histRsq_ElectronMultijet_ThreeBTag);
  AddLargestWeightUncertainty(histMR_Multijet_TwoBTagCorr);
  AddLargestWeightUncertainty(histMR_MuonMultijet_TwoBTagCorr);
  AddLargestWeightUncertainty(histMR_ElectronMultijet_TwoBTagCorr);


  histMR_Multijet_TwoBTag = NormalizeHist(histMR_Multijet_TwoBTag);
  histRsq_Multijet_TwoBTag = NormalizeHist(histRsq_Multijet_TwoBTag);
  histMR_Multijet_ThreeBTag = NormalizeHist(histMR_Multijet_ThreeBTag);
  histRsq_Multijet_ThreeBTag = NormalizeHist(histRsq_Multijet_ThreeBTag);
  // histMR_LooseLeptonMultijet_TwoBTag = NormalizeHist(histMR_LooseLeptonMultijet_TwoBTag);
  // histRsq_LooseLeptonMultijet_TwoBTag = NormalizeHist(histRsq_LooseLeptonMultijet_TwoBTag);
  // histMR_LooseLeptonMultijet_ThreeBTag = NormalizeHist(histMR_LooseLeptonMultijet_ThreeBTag);
  // histRsq_LooseLeptonMultijet_ThreeBTag = NormalizeHist(histRsq_LooseLeptonMultijet_ThreeBTag);
  histMR_MuonMultijet_TwoBTag = NormalizeHist(histMR_MuonMultijet_TwoBTag);
  histRsq_MuonMultijet_TwoBTag = NormalizeHist(histRsq_MuonMultijet_TwoBTag);
  histMR_MuonMultijet_ThreeBTag = NormalizeHist(histMR_MuonMultijet_ThreeBTag);
  histRsq_MuonMultijet_ThreeBTag = NormalizeHist(histRsq_MuonMultijet_ThreeBTag);
  histMR_ElectronMultijet_TwoBTag = NormalizeHist(histMR_ElectronMultijet_TwoBTag);
  histRsq_ElectronMultijet_TwoBTag = NormalizeHist(histRsq_ElectronMultijet_TwoBTag);
  histMR_ElectronMultijet_ThreeBTag = NormalizeHist(histMR_ElectronMultijet_ThreeBTag);
  histRsq_ElectronMultijet_ThreeBTag = NormalizeHist(histRsq_ElectronMultijet_ThreeBTag);
  histMR_Multijet_TwoBTagCorr = NormalizeHist(histMR_Multijet_TwoBTagCorr);
  histMR_MuonMultijet_TwoBTagCorr = NormalizeHist(histMR_MuonMultijet_TwoBTagCorr);
  histMR_ElectronMultijet_TwoBTagCorr = NormalizeHist(histMR_ElectronMultijet_TwoBTagCorr);

  //*****************************************************
  //Multijet MR Plot
  //*****************************************************

  TFile *file = TFile::Open("BTagShapeCorrection.root", "UPDATE");
  file->cd();
  // file->WriteTObject(ratioHist, "ratio", "WriteDelete");   
  // file->Close();
  

  cv = new TCanvas("cv","cv", 800,600);
  cv->SetHighLightColor(2);
  cv->SetFillColor(0);
  cv->SetBorderMode(0);
  cv->SetBorderSize(2);
  cv->SetLeftMargin(0.16);
  cv->SetRightMargin(0.3);
  cv->SetTopMargin(0.07);
  cv->SetBottomMargin(0.12);
  cv->SetFrameBorderMode(0);  

  pad1 = new TPad("pad1","pad1", 0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(0.04);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.60,0.70,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  if (useCorrected) {
    legend->AddEntry(histMR_Multijet_TwoBTagCorr, "2 b-tags Corrected", "PL");
  } else {
    legend->AddEntry(histMR_Multijet_TwoBTag, "2 b-tags", "PL");
  }
  legend->AddEntry(histMR_Multijet_ThreeBTag, "3 or more b-tags", "PL");

  ratioHist = (TH1D*)histMR_Multijet_TwoBTag->Clone("histRatioMultijet");
  ratioHistCorr = (TH1D*)histMR_Multijet_TwoBTagCorr->Clone("histRatioCorrMultijet");
  chisq = 0;
  ndof = 0;

  for (int b=0; b<ratioHist->GetXaxis()->GetNbins()+2; ++b) {
    double n = histMR_Multijet_ThreeBTag->GetBinContent(b);
    double nErr = histMR_Multijet_ThreeBTag->GetBinError(b);
    double d = histMR_Multijet_TwoBTag->GetBinContent(b);    
    double dErr = histMR_Multijet_TwoBTag->GetBinError(b);    
    if ( d > 0 ) {
      ratioHist->SetBinContent(b, n/d);
      if ( n > 0) {
	ratioHist->SetBinError(b, (n/d)*sqrt(pow(nErr/n,2) + pow(dErr/d,2) ) );
      } else {
	ratioHist->SetBinError(b, 1.0 );
      }
      cout << "Ratio : " << b << " : " << ratioHist->GetBinContent(b) << " +/- " << ratioHist->GetBinError(b) << "\n";
      chisq += fabs(ratioHist->GetBinContent(b) - 1) / ratioHist->GetBinError(b);
      ndof++;
    } else {
      ratioHist->SetBinContent(b,0);
      ratioHist->SetBinError(b,0);
    }
  }
  cout << "chisq / dof = " << chisq << "/" << ndof << " = " << chisq / ndof << "\n";
  for (int b=0; b<ratioHistCorr->GetXaxis()->GetNbins()+2; ++b) {
    double n = histMR_Multijet_ThreeBTag->GetBinContent(b);
    double nErr = histMR_Multijet_ThreeBTag->GetBinError(b);
    double d = histMR_Multijet_TwoBTagCorr->GetBinContent(b);    
    double dErr = histMR_Multijet_TwoBTagCorr->GetBinError(b);    
    if ( d > 0 ) {
      ratioHistCorr->SetBinContent(b, n/d);
      if ( n > 0) {
	ratioHistCorr->SetBinError(b, (n/d)*sqrt(pow(nErr/n,2) + pow(dErr/d,2) ) );
      } else {
	ratioHistCorr->SetBinError(b, 1.0 );
      }
      cout << "Corr Ratio : " << b << " : " << ratioHistCorr->GetBinContent(b) << " +/- " << ratioHistCorr->GetBinError(b) << "\n";
      chisq += fabs(ratioHistCorr->GetBinContent(b) - 1) / ratioHistCorr->GetBinError(b);
      ndof++;
    } else {
      ratioHistCorr->SetBinContent(b,0);
      ratioHistCorr->SetBinError(b,0);
    }
  }
  // cout << "chisq / dof = " << chisq << "/" << ndof << " = " << chisq / ndof << "\n";


  histMR_Multijet_TwoBTag->SetLineWidth(2);
  histMR_Multijet_TwoBTagCorr->SetLineWidth(2);
  histMR_Multijet_ThreeBTag->SetLineWidth(2);
  histMR_Multijet_TwoBTag->SetLineColor(kBlue);
  histMR_Multijet_TwoBTagCorr->SetLineColor(kBlue);
  histMR_Multijet_ThreeBTag->SetLineColor(kRed);

 histMR_Multijet_TwoBTagCorr->SetFillStyle(3002);
  histMR_Multijet_TwoBTagCorr->SetFillColor(kBlue);

  histMR_Multijet_ThreeBTag->Draw("e1same");
  if (useCorrected) {
    histMR_Multijet_TwoBTagCorr->Draw("e1same");
  } else {
    histMR_Multijet_TwoBTag->Draw("e1same");
  }

  histMR_Multijet_ThreeBTag->GetYaxis()->SetTitleOffset(1.2);
  histMR_Multijet_ThreeBTag->GetYaxis()->SetTitle("Fraction of Events");

  legend->Draw();

  cv->cd();
  cv->Update();

  pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
    
  ratioHist->GetYaxis()->SetTitle("Ratio");
  ratioHist->GetYaxis()->SetNdivisions(306);
  ratioHist->GetYaxis()->SetTitleSize(0.10);
  ratioHist->GetYaxis()->SetTitleOffset(0.3);
  ratioHist->GetYaxis()->SetRangeUser(0.0,2.5);
  ratioHist->GetYaxis()->SetLabelSize(0.10);
  ratioHist->GetXaxis()->SetLabelSize(0.125);
  ratioHist->GetXaxis()->SetTitleSize(0.15);
  ratioHist->GetXaxis()->SetTitleOffset(1.0);
  ratioHist->SetLineColor(kBlack);
  ratioHist->SetMarkerStyle(20);      
  ratioHist->SetMarkerSize(1);
  ratioHist->SetStats(false);
  ratioHist->Draw("pe");

  if (useCorrected) {
    ratioHistCorr->GetYaxis()->SetTitle("Ratio");
    ratioHistCorr->GetYaxis()->SetNdivisions(306);
    ratioHistCorr->GetYaxis()->SetTitleSize(0.10);
    ratioHistCorr->GetYaxis()->SetTitleOffset(0.3);
    ratioHistCorr->GetYaxis()->SetRangeUser(0.0,2.5);
    ratioHistCorr->GetYaxis()->SetLabelSize(0.10);
    ratioHistCorr->GetXaxis()->SetLabelSize(0.125);
    ratioHistCorr->GetXaxis()->SetTitleSize(0.15);
    ratioHistCorr->GetXaxis()->SetTitleOffset(1.0);
    ratioHistCorr->SetLineColor(kBlack);
    ratioHistCorr->SetMarkerStyle(20);      
    ratioHistCorr->SetMarkerSize(1);
    ratioHistCorr->SetStats(false);
    ratioHistCorr->Draw("pe");
  } else {
    ratioHist->GetYaxis()->SetTitle("Ratio");
    ratioHist->GetYaxis()->SetNdivisions(306);
    ratioHist->GetYaxis()->SetTitleSize(0.10);
    ratioHist->GetYaxis()->SetTitleOffset(0.3);
    ratioHist->GetYaxis()->SetRangeUser(0.0,2.5);
    ratioHist->GetYaxis()->SetLabelSize(0.10);
    ratioHist->GetXaxis()->SetLabelSize(0.125);
    ratioHist->GetXaxis()->SetTitleSize(0.15);
    ratioHist->GetXaxis()->SetTitleOffset(1.0);
    ratioHist->SetLineColor(kBlack);
    ratioHist->SetMarkerStyle(20);      
    ratioHist->SetMarkerSize(1);
    ratioHist->SetStats(false);
    ratioHist->Draw("pe");
  }
  
  
  // pad1->SetLogy(false);
  // cv->SaveAs("ShapeComparison_Multijet_MR_TwoBTagVsThreeBTag.gif");
  // cv->SaveAs("ShapeComparison_Multijet_MR_TwoBTagVsThreeBTag.pdf");
  pad1->SetLogy(true);
  if (useCorrected) {
    cv->SaveAs(Form("ShapeComparison_Multijet_MR_TwoBTagCorrVsThreeBTag%s_Logy.gif",BinLabel.c_str()));
    cv->SaveAs(Form("ShapeComparison_Multijet_MR_TwoBTagCorrVsThreeBTag%s_Logy.pdf",BinLabel.c_str()));
  } else {
    cv->SaveAs(Form("ShapeComparison_Multijet_MR_TwoBTagVsThreeBTag%s_Logy.gif",BinLabel.c_str()));
    cv->SaveAs(Form("ShapeComparison_Multijet_MR_TwoBTagVsThreeBTag%s_Logy.pdf",BinLabel.c_str()));
  }

  file->WriteTObject((TH1D*)ratioHist->Clone("ThreeBTagToTwoBTagRatio_MR_MultiJet"), "ThreeBTagToTwoBTagRatio_MR_MultiJet", "WriteDelete"); 


  //*****************************************************
  //Multijet Rsq Plot
  //*****************************************************

  cv = new TCanvas("cv","cv", 800,600);
  cv->SetHighLightColor(2);
  cv->SetFillColor(0);
  cv->SetBorderMode(0);
  cv->SetBorderSize(2);
  cv->SetLeftMargin(0.16);
  cv->SetRightMargin(0.3);
  cv->SetTopMargin(0.07);
  cv->SetBottomMargin(0.12);
  cv->SetFrameBorderMode(0);  

  pad1 = new TPad("pad1","pad1", 0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(0.04);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.60,0.70,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histRsq_Multijet_TwoBTag, "2 b-tags", "PL");
  legend->AddEntry(histRsq_Multijet_ThreeBTag, "3 or more b-tags", "PL");

  ratioHist = (TH1D*)histRsq_Multijet_TwoBTag->Clone("histRatioMultijet");
  chisq = 0;
  ndof = 0;
  for (int b=0; b<ratioHist->GetXaxis()->GetNbins()+2; ++b) {
    double n = histRsq_Multijet_ThreeBTag->GetBinContent(b);
    double nErr = histRsq_Multijet_ThreeBTag->GetBinError(b);
    double d = histRsq_Multijet_TwoBTag->GetBinContent(b);    
    double dErr = histRsq_Multijet_TwoBTag->GetBinError(b);    
    if ( d > 0 ) {
      ratioHist->SetBinContent(b, n/d);
      if ( n > 0) {
	ratioHist->SetBinError(b, (n/d)*sqrt(pow(nErr/n,2) + pow(dErr/d,2) ) );
      } else {
	ratioHist->SetBinError(b, 1.0 );
      }
      cout << "Ratio : " << b << " : " << ratioHist->GetBinContent(b) << " +/- " << ratioHist->GetBinError(b) << "\n";
      if (b<16) {
	chisq += fabs(ratioHist->GetBinContent(b) - 1) / ratioHist->GetBinError(b);
	ndof++;
      }
    } else {
      ratioHist->SetBinContent(b,0);
      ratioHist->SetBinError(b,0);
    }
  }
  cout << "chisq / dof = " << chisq << "/" << ndof << " = " << chisq / ndof << "\n";

  histRsq_Multijet_TwoBTag->SetLineWidth(2);
  histRsq_Multijet_ThreeBTag->SetLineWidth(2);
  histRsq_Multijet_TwoBTag->SetLineColor(kBlue);
  histRsq_Multijet_ThreeBTag->SetLineColor(kRed);
  histRsq_Multijet_TwoBTag->Draw("e1");
  histRsq_Multijet_ThreeBTag->Draw("e1same");

  histRsq_Multijet_TwoBTag->GetYaxis()->SetTitleOffset(1.2);
  histRsq_Multijet_TwoBTag->GetYaxis()->SetTitle("Fraction of Events");
  histRsq_Multijet_TwoBTag->GetXaxis()->SetTitle("R^{2}");

  legend->Draw();

  cv->cd();
  cv->Update();

  pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
    
  ratioHist->GetXaxis()->SetTitle("R^{2}");
  ratioHist->GetYaxis()->SetTitle("Ratio");
  ratioHist->GetYaxis()->SetNdivisions(306);
  ratioHist->GetYaxis()->SetTitleSize(0.10);
  ratioHist->GetYaxis()->SetTitleOffset(0.3);
  ratioHist->GetYaxis()->SetRangeUser(0.0,2.5);
  ratioHist->GetYaxis()->SetLabelSize(0.10);
  ratioHist->GetXaxis()->SetLabelSize(0.125);
  ratioHist->GetXaxis()->SetTitleSize(0.15);
  ratioHist->GetXaxis()->SetTitleOffset(1.0);
  ratioHist->SetLineColor(kBlack);
  ratioHist->SetMarkerStyle(20);      
  ratioHist->SetMarkerSize(1);
  ratioHist->SetStats(false);
  ratioHist->Draw("pe");
  
  // pad1->SetLogy(false);
  // cv->SaveAs("ShapeComparison_Multijet_Rsq_TwoBTagVsThreeBTag.gif");
  // cv->SaveAs("ShapeComparison_Multijet_Rsq_TwoBTagVsThreeBTag.pdf");
  pad1->SetLogy(true);
  cv->SaveAs("ShapeComparison_Multijet_Rsq_TwoBTagVsThreeBTag_Logy.gif");
  cv->SaveAs("ShapeComparison_Multijet_Rsq_TwoBTagVsThreeBTag_Logy.pdf");


  file->WriteTObject((TH1D*)ratioHist->Clone("ThreeBTagToTwoBTagRatio_Rsq_MultiJet"), "ThreeBTagToTwoBTagRatio_Rsq_MultiJet", "WriteDelete"); 





  //*****************************************************
  //MuonMultijet MR Plot
  //*****************************************************

  cv = new TCanvas("cv","cv", 800,600);
  cv->SetHighLightColor(2);
  cv->SetFillColor(0);
  cv->SetBorderMode(0);
  cv->SetBorderSize(2);
  cv->SetLeftMargin(0.16);
  cv->SetRightMargin(0.3);
  cv->SetTopMargin(0.07);
  cv->SetBottomMargin(0.12);
  cv->SetFrameBorderMode(0);  

  pad1 = new TPad("pad1","pad1", 0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(0.04);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.60,0.70,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  if (useCorrected) {
    legend->AddEntry(histMR_MuonMultijet_TwoBTagCorr, "2 b-tags Corrected", "PL");
  } else {
    legend->AddEntry(histMR_MuonMultijet_TwoBTag, "2 b-tags", "PL");
  }
  legend->AddEntry(histMR_MuonMultijet_ThreeBTag, "3 or more b-tags", "PL");

  ratioHist = (TH1D*)histMR_MuonMultijet_TwoBTag->Clone("histRatioMuonMultijet");
  ratioHistCorr = (TH1D*)histMR_MuonMultijet_TwoBTagCorr->Clone("histRatioCorrMuonMultijet");
  chisq = 0;
  ndof = 0;

  for (int b=0; b<ratioHist->GetXaxis()->GetNbins()+2; ++b) {
    double n = histMR_MuonMultijet_ThreeBTag->GetBinContent(b);
    double nErr = histMR_MuonMultijet_ThreeBTag->GetBinError(b);
    double d = histMR_MuonMultijet_TwoBTag->GetBinContent(b);    
    double dErr = histMR_MuonMultijet_TwoBTag->GetBinError(b);    
    if ( d > 0 ) {
      ratioHist->SetBinContent(b, n/d);
      if ( n > 0) {
	ratioHist->SetBinError(b, (n/d)*sqrt(pow(nErr/n,2) + pow(dErr/d,2) ) );
      } else {
	ratioHist->SetBinError(b, 1.0 );
      }
      cout << "Ratio : " << b << " : " << ratioHist->GetBinContent(b) << " +/- " << ratioHist->GetBinError(b) << "\n";
      chisq += fabs(ratioHist->GetBinContent(b) - 1) / ratioHist->GetBinError(b);
      ndof++;
    } else {
      ratioHist->SetBinContent(b,0);
      ratioHist->SetBinError(b,0);
    }
  }
  cout << "chisq / dof = " << chisq << "/" << ndof << " = " << chisq / ndof << "\n";
  for (int b=0; b<ratioHistCorr->GetXaxis()->GetNbins()+2; ++b) {
    double n = histMR_MuonMultijet_ThreeBTag->GetBinContent(b);
    double nErr = histMR_MuonMultijet_ThreeBTag->GetBinError(b);
    double d = histMR_MuonMultijet_TwoBTagCorr->GetBinContent(b);    
    double dErr = histMR_MuonMultijet_TwoBTagCorr->GetBinError(b);    
    if ( d > 0 ) {
      ratioHistCorr->SetBinContent(b, n/d);
      if ( n > 0) {
	ratioHistCorr->SetBinError(b, (n/d)*sqrt(pow(nErr/n,2) + pow(dErr/d,2) ) );
      } else {
	ratioHistCorr->SetBinError(b, 1.0 );
      }
      cout << "Corr Ratio : " << b << " : " << ratioHistCorr->GetBinContent(b) << " +/- " << ratioHistCorr->GetBinError(b) << "\n";
      chisq += fabs(ratioHistCorr->GetBinContent(b) - 1) / ratioHistCorr->GetBinError(b);
      ndof++;
    } else {
      ratioHistCorr->SetBinContent(b,0);
      ratioHistCorr->SetBinError(b,0);
    }
  }
  // cout << "chisq / dof = " << chisq << "/" << ndof << " = " << chisq / ndof << "\n";


  histMR_MuonMultijet_TwoBTag->SetLineWidth(2);
  histMR_MuonMultijet_TwoBTagCorr->SetLineWidth(2);
  histMR_MuonMultijet_ThreeBTag->SetLineWidth(2);
  histMR_MuonMultijet_TwoBTag->SetLineColor(kBlue);
  histMR_MuonMultijet_TwoBTagCorr->SetLineColor(kBlue);
  histMR_MuonMultijet_ThreeBTag->SetLineColor(kRed);

 histMR_MuonMultijet_TwoBTagCorr->SetFillStyle(3002);
  histMR_MuonMultijet_TwoBTagCorr->SetFillColor(kBlue);

  histMR_MuonMultijet_ThreeBTag->Draw("e1same");
  if (useCorrected) {
    histMR_MuonMultijet_TwoBTagCorr->Draw("e1same");
  } else {
    histMR_MuonMultijet_TwoBTag->Draw("e1same");
  }

  histMR_MuonMultijet_ThreeBTag->GetYaxis()->SetTitleOffset(1.2);
  histMR_MuonMultijet_ThreeBTag->GetYaxis()->SetTitle("Fraction of Events");

  legend->Draw();

  cv->cd();
  cv->Update();

  pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
    
  ratioHist->GetYaxis()->SetTitle("Ratio");
  ratioHist->GetYaxis()->SetNdivisions(306);
  ratioHist->GetYaxis()->SetTitleSize(0.10);
  ratioHist->GetYaxis()->SetTitleOffset(0.3);
  ratioHist->GetYaxis()->SetRangeUser(0.0,2.5);
  ratioHist->GetYaxis()->SetLabelSize(0.10);
  ratioHist->GetXaxis()->SetLabelSize(0.125);
  ratioHist->GetXaxis()->SetTitleSize(0.15);
  ratioHist->GetXaxis()->SetTitleOffset(1.0);
  ratioHist->SetLineColor(kBlack);
  ratioHist->SetMarkerStyle(20);      
  ratioHist->SetMarkerSize(1);
  ratioHist->SetStats(false);
  ratioHist->Draw("pe");

  if (useCorrected) {
    ratioHistCorr->GetYaxis()->SetTitle("Ratio");
    ratioHistCorr->GetYaxis()->SetNdivisions(306);
    ratioHistCorr->GetYaxis()->SetTitleSize(0.10);
    ratioHistCorr->GetYaxis()->SetTitleOffset(0.3);
    ratioHistCorr->GetYaxis()->SetRangeUser(0.0,2.5);
    ratioHistCorr->GetYaxis()->SetLabelSize(0.10);
    ratioHistCorr->GetXaxis()->SetLabelSize(0.125);
    ratioHistCorr->GetXaxis()->SetTitleSize(0.15);
    ratioHistCorr->GetXaxis()->SetTitleOffset(1.0);
    ratioHistCorr->SetLineColor(kBlack);
    ratioHistCorr->SetMarkerStyle(20);      
    ratioHistCorr->SetMarkerSize(1);
    ratioHistCorr->SetStats(false);
    ratioHistCorr->Draw("pe");
  } else {
    ratioHist->GetYaxis()->SetTitle("Ratio");
    ratioHist->GetYaxis()->SetNdivisions(306);
    ratioHist->GetYaxis()->SetTitleSize(0.10);
    ratioHist->GetYaxis()->SetTitleOffset(0.3);
    ratioHist->GetYaxis()->SetRangeUser(0.0,2.5);
    ratioHist->GetYaxis()->SetLabelSize(0.10);
    ratioHist->GetXaxis()->SetLabelSize(0.125);
    ratioHist->GetXaxis()->SetTitleSize(0.15);
    ratioHist->GetXaxis()->SetTitleOffset(1.0);
    ratioHist->SetLineColor(kBlack);
    ratioHist->SetMarkerStyle(20);      
    ratioHist->SetMarkerSize(1);
    ratioHist->SetStats(false);
    ratioHist->Draw("pe");
  }
  
  
  // pad1->SetLogy(false);
  // cv->SaveAs("ShapeComparison_MuonMultijet_MR_TwoBTagVsThreeBTag.gif");
  // cv->SaveAs("ShapeComparison_MuonMultijet_MR_TwoBTagVsThreeBTag.pdf");
  pad1->SetLogy(true);
  if (useCorrected) {
    cv->SaveAs(Form("ShapeComparison_MuonMultijet_MR_TwoBTagCorrVsThreeBTag%s_Logy.gif",BinLabel.c_str()));
    cv->SaveAs(Form("ShapeComparison_MuonMultijet_MR_TwoBTagCorrVsThreeBTag%s_Logy.pdf",BinLabel.c_str()));
  } else {
    cv->SaveAs(Form("ShapeComparison_MuonMultijet_MR_TwoBTagVsThreeBTag%s_Logy.gif",BinLabel.c_str()));
    cv->SaveAs(Form("ShapeComparison_MuonMultijet_MR_TwoBTagVsThreeBTag%s_Logy.pdf",BinLabel.c_str()));
  }

  file->WriteTObject((TH1D*)ratioHist->Clone("ThreeBTagToTwoBTagRatio_MR_MuMultiJet"), "ThreeBTagToTwoBTagRatio_MR_MuMultiJet", "WriteDelete"); 



  //*****************************************************
  //MuonMultijet Rsq Plot
  //*****************************************************

  cv = new TCanvas("cv","cv", 800,600);
  cv->SetHighLightColor(2);
  cv->SetFillColor(0);
  cv->SetBorderMode(0);
  cv->SetBorderSize(2);
  cv->SetLeftMargin(0.16);
  cv->SetRightMargin(0.3);
  cv->SetTopMargin(0.07);
  cv->SetBottomMargin(0.12);
  cv->SetFrameBorderMode(0);  

  pad1 = new TPad("pad1","pad1", 0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(0.04);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.60,0.70,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histRsq_MuonMultijet_TwoBTag, "2 b-tags", "PL");
  legend->AddEntry(histRsq_MuonMultijet_ThreeBTag, "3 or more b-tags", "PL");

  ratioHist = (TH1D*)histRsq_MuonMultijet_TwoBTag->Clone("histRatioMuonMultijet");
  chisq = 0;
  ndof = 0;
  for (int b=0; b<ratioHist->GetXaxis()->GetNbins()+2; ++b) {
    double n = histRsq_MuonMultijet_ThreeBTag->GetBinContent(b);
    double nErr = histRsq_MuonMultijet_ThreeBTag->GetBinError(b);
    double d = histRsq_MuonMultijet_TwoBTag->GetBinContent(b);    
    double dErr = histRsq_MuonMultijet_TwoBTag->GetBinError(b);    
    if ( d > 0 ) {
      ratioHist->SetBinContent(b, n/d);
      if ( n > 0) {
	ratioHist->SetBinError(b, (n/d)*sqrt(pow(nErr/n,2) + pow(dErr/d,2) ) );
      } else {
	ratioHist->SetBinError(b, 1.0 );
      }
      cout << "Ratio : " << b << " : " << ratioHist->GetBinContent(b) << " +/- " << ratioHist->GetBinError(b) << "\n";
      if (b<16){
	chisq += fabs(ratioHist->GetBinContent(b) - 1) / ratioHist->GetBinError(b);
	ndof++;
      }
    } else {
      ratioHist->SetBinContent(b,0);
      ratioHist->SetBinError(b,0);
    }
  }
  cout << "chisq / dof = " << chisq << "/" << ndof << " = " << chisq / ndof << "\n";

  histRsq_MuonMultijet_TwoBTag->SetLineWidth(2);
  histRsq_MuonMultijet_ThreeBTag->SetLineWidth(2);
  histRsq_MuonMultijet_TwoBTag->SetLineColor(kBlue);
  histRsq_MuonMultijet_ThreeBTag->SetLineColor(kRed);
  histRsq_MuonMultijet_TwoBTag->Draw("e1");
  histRsq_MuonMultijet_ThreeBTag->Draw("e1same");

  histRsq_MuonMultijet_TwoBTag->GetYaxis()->SetTitleOffset(1.2);
  histRsq_MuonMultijet_TwoBTag->GetYaxis()->SetTitle("Fraction of Events");
  histRsq_MuonMultijet_TwoBTag->GetXaxis()->SetTitle("R^{2}");

  legend->Draw();

  cv->cd();
  cv->Update();

  pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
    
  ratioHist->GetXaxis()->SetTitle("R^{2}");
  ratioHist->GetYaxis()->SetTitle("Ratio");
  ratioHist->GetYaxis()->SetNdivisions(306);
  ratioHist->GetYaxis()->SetTitleSize(0.10);
  ratioHist->GetYaxis()->SetTitleOffset(0.3);
  ratioHist->GetYaxis()->SetRangeUser(0.0,2.5);
  ratioHist->GetYaxis()->SetLabelSize(0.10);
  ratioHist->GetXaxis()->SetLabelSize(0.125);
  ratioHist->GetXaxis()->SetTitleSize(0.15);
  ratioHist->GetXaxis()->SetTitleOffset(1.0);
  ratioHist->SetLineColor(kBlack);
  ratioHist->SetMarkerStyle(20);      
  ratioHist->SetMarkerSize(1);
  ratioHist->SetStats(false);
  ratioHist->Draw("pe");
  
  // pad1->SetLogy(false);
  // cv->SaveAs("ShapeComparison_MuonMultijet_Rsq_TwoBTagVsThreeBTag.gif");
  // cv->SaveAs("ShapeComparison_MuonMultijet_Rsq_TwoBTagVsThreeBTag.pdf");
  pad1->SetLogy(true);
  cv->SaveAs("ShapeComparison_MuonMultijet_Rsq_TwoBTagVsThreeBTag_Logy.gif");
  cv->SaveAs("ShapeComparison_MuonMultijet_Rsq_TwoBTagVsThreeBTag_Logy.pdf");

  file->WriteTObject((TH1D*)ratioHist->Clone("ThreeBTagToTwoBTagRatio_Rsq_MuonMultiJet"), "ThreeBTagToTwoBTagRatio_Rsq_MuonMultiJet", "WriteDelete"); 





  //*****************************************************
  //ElectronMultijet MR Plot
  //*****************************************************

  cv = new TCanvas("cv","cv", 800,600);
  cv->SetHighLightColor(2);
  cv->SetFillColor(0);
  cv->SetBorderMode(0);
  cv->SetBorderSize(2);
  cv->SetLeftMargin(0.16);
  cv->SetRightMargin(0.3);
  cv->SetTopMargin(0.07);
  cv->SetBottomMargin(0.12);
  cv->SetFrameBorderMode(0);  

  pad1 = new TPad("pad1","pad1", 0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(0.04);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.60,0.70,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  if (useCorrected) {
    legend->AddEntry(histMR_ElectronMultijet_TwoBTagCorr, "2 b-tags Corrected", "PL");
  } else {
    legend->AddEntry(histMR_ElectronMultijet_TwoBTag, "2 b-tags", "PL");
  }
  legend->AddEntry(histMR_ElectronMultijet_ThreeBTag, "3 or more b-tags", "PL");

  ratioHist = (TH1D*)histMR_ElectronMultijet_TwoBTag->Clone("histRatioElectronMultijet");
  ratioHistCorr = (TH1D*)histMR_ElectronMultijet_TwoBTagCorr->Clone("histRatioCorrElectronMultijet");
  chisq = 0;
  ndof = 0;

  for (int b=0; b<ratioHist->GetXaxis()->GetNbins()+2; ++b) {
    double n = histMR_ElectronMultijet_ThreeBTag->GetBinContent(b);
    double nErr = histMR_ElectronMultijet_ThreeBTag->GetBinError(b);
    double d = histMR_ElectronMultijet_TwoBTag->GetBinContent(b);    
    double dErr = histMR_ElectronMultijet_TwoBTag->GetBinError(b);    
    if ( d > 0 ) {
      ratioHist->SetBinContent(b, n/d);
      if ( n > 0) {
	ratioHist->SetBinError(b, (n/d)*sqrt(pow(nErr/n,2) + pow(dErr/d,2) ) );
      } else {
	ratioHist->SetBinError(b, 1.0 );
      }
      cout << "Ratio : " << b << " : " << ratioHist->GetBinContent(b) << " +/- " << ratioHist->GetBinError(b) << "\n";
      chisq += fabs(ratioHist->GetBinContent(b) - 1) / ratioHist->GetBinError(b);
      ndof++;
    } else {
      ratioHist->SetBinContent(b,0);
      ratioHist->SetBinError(b,0);
    }
  }
  cout << "chisq / dof = " << chisq << "/" << ndof << " = " << chisq / ndof << "\n";
  for (int b=0; b<ratioHistCorr->GetXaxis()->GetNbins()+2; ++b) {
    double n = histMR_ElectronMultijet_ThreeBTag->GetBinContent(b);
    double nErr = histMR_ElectronMultijet_ThreeBTag->GetBinError(b);
    double d = histMR_ElectronMultijet_TwoBTagCorr->GetBinContent(b);    
    double dErr = histMR_ElectronMultijet_TwoBTagCorr->GetBinError(b);    
    if ( d > 0 ) {
      ratioHistCorr->SetBinContent(b, n/d);
      if ( n > 0) {
	ratioHistCorr->SetBinError(b, (n/d)*sqrt(pow(nErr/n,2) + pow(dErr/d,2) ) );
      } else {
	ratioHistCorr->SetBinError(b, 1.0 );
      }
      cout << "Corr Ratio : " << b << " : " << ratioHistCorr->GetBinContent(b) << " +/- " << ratioHistCorr->GetBinError(b) << "\n";
      chisq += fabs(ratioHistCorr->GetBinContent(b) - 1) / ratioHistCorr->GetBinError(b);
      ndof++;
    } else {
      ratioHistCorr->SetBinContent(b,0);
      ratioHistCorr->SetBinError(b,0);
    }
  }
  // cout << "chisq / dof = " << chisq << "/" << ndof << " = " << chisq / ndof << "\n";


  histMR_ElectronMultijet_TwoBTag->SetLineWidth(2);
  histMR_ElectronMultijet_TwoBTagCorr->SetLineWidth(2);
  histMR_ElectronMultijet_ThreeBTag->SetLineWidth(2);
  histMR_ElectronMultijet_TwoBTag->SetLineColor(kBlue);
  histMR_ElectronMultijet_TwoBTagCorr->SetLineColor(kBlue);
  histMR_ElectronMultijet_ThreeBTag->SetLineColor(kRed);

 histMR_ElectronMultijet_TwoBTagCorr->SetFillStyle(3002);
  histMR_ElectronMultijet_TwoBTagCorr->SetFillColor(kBlue);

  histMR_ElectronMultijet_ThreeBTag->Draw("e1same");
  if (useCorrected) {
    histMR_ElectronMultijet_TwoBTagCorr->Draw("e1same");
  } else {
    histMR_ElectronMultijet_TwoBTag->Draw("e1same");
  }

  histMR_ElectronMultijet_ThreeBTag->GetYaxis()->SetTitleOffset(1.2);
  histMR_ElectronMultijet_ThreeBTag->GetYaxis()->SetTitle("Fraction of Events");

  legend->Draw();

  cv->cd();
  cv->Update();

  pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
    
  ratioHist->GetYaxis()->SetTitle("Ratio");
  ratioHist->GetYaxis()->SetNdivisions(306);
  ratioHist->GetYaxis()->SetTitleSize(0.10);
  ratioHist->GetYaxis()->SetTitleOffset(0.3);
  ratioHist->GetYaxis()->SetRangeUser(0.0,2.5);
  ratioHist->GetYaxis()->SetLabelSize(0.10);
  ratioHist->GetXaxis()->SetLabelSize(0.125);
  ratioHist->GetXaxis()->SetTitleSize(0.15);
  ratioHist->GetXaxis()->SetTitleOffset(1.0);
  ratioHist->SetLineColor(kBlack);
  ratioHist->SetMarkerStyle(20);      
  ratioHist->SetMarkerSize(1);
  ratioHist->SetStats(false);
  ratioHist->Draw("pe");

  if (useCorrected) {
    ratioHistCorr->GetYaxis()->SetTitle("Ratio");
    ratioHistCorr->GetYaxis()->SetNdivisions(306);
    ratioHistCorr->GetYaxis()->SetTitleSize(0.10);
    ratioHistCorr->GetYaxis()->SetTitleOffset(0.3);
    ratioHistCorr->GetYaxis()->SetRangeUser(0.0,2.5);
    ratioHistCorr->GetYaxis()->SetLabelSize(0.10);
    ratioHistCorr->GetXaxis()->SetLabelSize(0.125);
    ratioHistCorr->GetXaxis()->SetTitleSize(0.15);
    ratioHistCorr->GetXaxis()->SetTitleOffset(1.0);
    ratioHistCorr->SetLineColor(kBlack);
    ratioHistCorr->SetMarkerStyle(20);      
    ratioHistCorr->SetMarkerSize(1);
    ratioHistCorr->SetStats(false);
    ratioHistCorr->Draw("pe");
  } else {
    ratioHist->GetYaxis()->SetTitle("Ratio");
    ratioHist->GetYaxis()->SetNdivisions(306);
    ratioHist->GetYaxis()->SetTitleSize(0.10);
    ratioHist->GetYaxis()->SetTitleOffset(0.3);
    ratioHist->GetYaxis()->SetRangeUser(0.0,2.5);
    ratioHist->GetYaxis()->SetLabelSize(0.10);
    ratioHist->GetXaxis()->SetLabelSize(0.125);
    ratioHist->GetXaxis()->SetTitleSize(0.15);
    ratioHist->GetXaxis()->SetTitleOffset(1.0);
    ratioHist->SetLineColor(kBlack);
    ratioHist->SetMarkerStyle(20);      
    ratioHist->SetMarkerSize(1);
    ratioHist->SetStats(false);
    ratioHist->Draw("pe");
  }
  
  
  // pad1->SetLogy(false);
  // cv->SaveAs("ShapeComparison_ElectronMultijet_MR_TwoBTagVsThreeBTag.gif");
  // cv->SaveAs("ShapeComparison_ElectronMultijet_MR_TwoBTagVsThreeBTag.pdf");
  pad1->SetLogy(true);
  if (useCorrected) {
    cv->SaveAs(Form("ShapeComparison_ElectronMultijet_MR_TwoBTagCorrVsThreeBTag%s_Logy.gif",BinLabel.c_str()));
    cv->SaveAs(Form("ShapeComparison_ElectronMultijet_MR_TwoBTagCorrVsThreeBTag%s_Logy.pdf",BinLabel.c_str()));
  } else {
    cv->SaveAs(Form("ShapeComparison_ElectronMultijet_MR_TwoBTagVsThreeBTag%s_Logy.gif",BinLabel.c_str()));
    cv->SaveAs(Form("ShapeComparison_ElectronMultijet_MR_TwoBTagVsThreeBTag%s_Logy.pdf",BinLabel.c_str()));
  }

  file->WriteTObject((TH1D*)ratioHist->Clone("ThreeBTagToTwoBTagRatio_MR_MuMultiJet"), "ThreeBTagToTwoBTagRatio_MR_MuMultiJet", "WriteDelete"); 


  //*****************************************************
  //ElectronMultijet Rsq Plot
  //*****************************************************

  cv = new TCanvas("cv","cv", 800,600);
  cv->SetHighLightColor(2);
  cv->SetFillColor(0);
  cv->SetBorderMode(0);
  cv->SetBorderSize(2);
  cv->SetLeftMargin(0.16);
  cv->SetRightMargin(0.3);
  cv->SetTopMargin(0.07);
  cv->SetBottomMargin(0.12);
  cv->SetFrameBorderMode(0);  

  pad1 = new TPad("pad1","pad1", 0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(0.04);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.60,0.70,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histRsq_ElectronMultijet_TwoBTag, "2 b-tags", "PL");
  legend->AddEntry(histRsq_ElectronMultijet_ThreeBTag, "3 or more b-tags", "PL");

  ratioHist = (TH1D*)histRsq_ElectronMultijet_TwoBTag->Clone("histRatioElectronMultijet");
  chisq = 0;
  ndof = 0;
  for (int b=0; b<ratioHist->GetXaxis()->GetNbins()+2; ++b) {
    double n = histRsq_ElectronMultijet_ThreeBTag->GetBinContent(b);
    double nErr = histRsq_ElectronMultijet_ThreeBTag->GetBinError(b);
    double d = histRsq_ElectronMultijet_TwoBTag->GetBinContent(b);    
    double dErr = histRsq_ElectronMultijet_TwoBTag->GetBinError(b);    
    if ( d > 0 ) {
      ratioHist->SetBinContent(b, n/d);
      if ( n > 0) {
	ratioHist->SetBinError(b, (n/d)*sqrt(pow(nErr/n,2) + pow(dErr/d,2) ) );
      } else {
	ratioHist->SetBinError(b, 1.0 );
      }
      cout << "Ratio : " << b << " : " << ratioHist->GetBinContent(b) << " +/- " << ratioHist->GetBinError(b) << "\n";
      if (b<16){
	chisq += fabs(ratioHist->GetBinContent(b) - 1) / ratioHist->GetBinError(b);
	ndof++;
      }
    } else {
      ratioHist->SetBinContent(b,0);
      ratioHist->SetBinError(b,0);
    }
  }
  cout << "chisq / dof = " << chisq << "/" << ndof << " = " << chisq / ndof << "\n";

  histRsq_ElectronMultijet_TwoBTag->SetLineWidth(2);
  histRsq_ElectronMultijet_ThreeBTag->SetLineWidth(2);
  histRsq_ElectronMultijet_TwoBTag->SetLineColor(kBlue);
  histRsq_ElectronMultijet_ThreeBTag->SetLineColor(kRed);
  histRsq_ElectronMultijet_TwoBTag->Draw("e1");
  histRsq_ElectronMultijet_ThreeBTag->Draw("e1same");

  histRsq_ElectronMultijet_TwoBTag->GetYaxis()->SetTitleOffset(1.2);
  histRsq_ElectronMultijet_TwoBTag->GetYaxis()->SetTitle("Fraction of Events");
  histRsq_ElectronMultijet_TwoBTag->GetXaxis()->SetTitle("R^{2}");

  legend->Draw();

  cv->cd();
  cv->Update();

  pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
    
  ratioHist->GetXaxis()->SetTitle("R^{2}");
  ratioHist->GetYaxis()->SetTitle("Ratio");
  ratioHist->GetYaxis()->SetNdivisions(306);
  ratioHist->GetYaxis()->SetTitleSize(0.10);
  ratioHist->GetYaxis()->SetTitleOffset(0.3);
  ratioHist->GetYaxis()->SetRangeUser(0.0,2.5);
  ratioHist->GetYaxis()->SetLabelSize(0.10);
  ratioHist->GetXaxis()->SetLabelSize(0.125);
  ratioHist->GetXaxis()->SetTitleSize(0.15);
  ratioHist->GetXaxis()->SetTitleOffset(1.0);
  ratioHist->SetLineColor(kBlack);
  ratioHist->SetMarkerStyle(20);      
  ratioHist->SetMarkerSize(1);
  ratioHist->SetStats(false);
  ratioHist->Draw("pe");
  
  // pad1->SetLogy(false);
  // cv->SaveAs("ShapeComparison_ElectronMultijet_Rsq_TwoBTagVsThreeBTag.gif");
  // cv->SaveAs("ShapeComparison_ElectronMultijet_Rsq_TwoBTagVsThreeBTag.pdf");
  pad1->SetLogy(true);
  cv->SaveAs("ShapeComparison_ElectronMultijet_Rsq_TwoBTagVsThreeBTag_Logy.gif");
  cv->SaveAs("ShapeComparison_ElectronMultijet_Rsq_TwoBTagVsThreeBTag_Logy.pdf");

  file->WriteTObject((TH1D*)ratioHist->Clone("ThreeBTagToTwoBTagRatio_Rsq_ElectronMultiJet"), "ThreeBTagToTwoBTagRatio_Rsq_ElectronMultiJet", "WriteDelete"); 

  file->Close();



}

void PlotDataAndStackedBkg( vector<TH1D*> hist , vector<string> processLabels, vector<int> color,  bool hasData, string varName, string label ) {

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

  legend = new TLegend(0.60,0.50,0.90,0.84);
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
    stack->SetMaximum( 1.2* fmax( stack->GetMaximum(), hist[0]->GetMaximum()) );
    stack->SetMinimum( 0.1 );

    if (hasData) {
      hist[0]->SetMarkerStyle(20);      
      hist[0]->SetMarkerSize(1);
      hist[0]->SetLineWidth(1);
      hist[0]->SetLineColor(color[0]);
      hist[0]->Draw("pesame");
      // DataMean = hist[0]->GetMean();
      // DataRMS = hist[0]->GetRMS();
    }
    legend->Draw();

    // for (uint i = 0 ; i < hist.size(); ++i) {
    //   if (processLabels[i] == "Data") {
    // 	DataMean = hist[i]->GetMean();
    // 	DataRMS = hist[i]->GetRMS();
    //   } 
    //   if (processLabels[i] == "DY") {
    // 	MCMean = hist[i]->GetMean();
    // 	MCRMS = hist[i]->GetRMS();
    //   } 
    // }

  }

  //****************************
  //Add CMS and Lumi Labels
  //****************************
  lumi_13TeV = "1.25 fb^{-1}";
  //lumi_13TeV = "Run257396-257400";
  writeExtraText = true;
  relPosX = 0.13;
  CMS_lumi(pad1,4,0);

  // TLatex *StatLabels  = new TLatex;
  // StatLabels->SetTextSize(0.03);
  // StatLabels->DrawLatexNDC(0.8,0.8, Form("Data Mean: %.1f",DataMean));
  // StatLabels->DrawLatexNDC(0.8,0.76, Form("Data RMS: %.1f",DataRMS));
  // StatLabels->DrawLatexNDC(0.8,0.72, Form("MC Mean: %.1f",MCMean));
  // StatLabels->DrawLatexNDC(0.8,0.68, Form("MC RMS: %.1f",MCRMS));


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
  histDataOverMC->GetYaxis()->SetRangeUser(0.5,1.5);
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
  cv->SaveAs(Form("RazorPlots_%s%s.png",varName.c_str(), label.c_str()));
  cv->SaveAs(Form("RazorPlots_%s%s.pdf",varName.c_str(), label.c_str()));
  
  pad1->SetLogy(true);
  cv->SaveAs(Form("RazorPlots_%s%s_Logy.png",varName.c_str(),label.c_str()));
  cv->SaveAs(Form("RazorPlots_%s%s_Logy.pdf",varName.c_str(),label.c_str()));


 

}


//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void RunMakeRazorPlots ( vector<vector<string> > bkgfiles, vector<string> bkgLabels, vector<int> bkgColors, double lumi, int boxOption = 0, int btagOption = -1, int RsqBinOption = -1, string label = "", string latexlabel = "") {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  double intLumi = lumi; //in units of pb^-1
  string Label = "";
  if (label != "") Label = "_" + label;

  vector<vector<string> > inputfiles;
  vector<string> processLabels;
  vector<int> color;

  assert(bkgfiles.size() == bkgLabels.size());
  assert(bkgfiles.size() == bkgColors.size());
  for (int i=0; i < int(bkgfiles.size()); ++i) {
    inputfiles.push_back(bkgfiles[i]);
    processLabels.push_back(bkgLabels[i]);
    color.push_back(bkgColors[i]);
  }

  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  float MRBinLowEdges_MultiJet[] = {500, 600, 700, 900, 1200, 1600, 2500, 4000};
  float RsqBinLowEdges_MultiJet[] = {0.25, 0.30, 0.41, 0.52, 0.64, 1.5};
  const int nMRBins_MultiJet = sizeof(MRBinLowEdges_MultiJet)/sizeof(float)-1;
  const int nRsqBins_MultiJet = sizeof(RsqBinLowEdges_MultiJet)/sizeof(float)-1;
  float MRBinLowEdges_LeptonMultiJet[] = {400, 450, 500, 600, 700, 900, 1200, 1600, 2500, 4000};
  float RsqBinLowEdges_LeptonMultiJet[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 1.5};
  const int nMRBins_LeptonMultiJet = sizeof(MRBinLowEdges_LeptonMultiJet)/sizeof(float)-1;
  const int nRsqBins_LeptonMultiJet = sizeof(RsqBinLowEdges_LeptonMultiJet)/sizeof(float)-1;

 
  TH1D* histMRAllBkg = 0;
  TH1D* histRsqAllBkg = 0;
  TH1D* histMRAllBkgCorr = 0;
  TH1D* histRsqAllBkgCorr = 0;
  // if (boxOption == 100 || boxOption == 101) {
  //   histMRAllBkg =  new TH1D( "MRAllBkg",";M_{R} [GeV/c^{2}];Number of Events", 25, 400, 2400);
  //   histRsqAllBkg =  new TH1D( "RsqAllBkg", ";M_{R} [GeV/c^{2}];Number of Events", 25, 0.25, 1.5);
  // } else {
  //   histMRAllBkg =  new TH1D( "MRAllBkg",";M_{R} [GeV/c^{2}];Number of Events", 25, 300, 2300);
  //   histRsqAllBkg =  new TH1D( "RsqAllBkg", ";M_{R} [GeV/c^{2}];Number of Events", 27, 0.15, 1.5);
  // }
  if (boxOption == 100 || boxOption == 101) {
    histMRAllBkg =  new TH1D( "MRAllBkg",";M_{R} [GeV/c^{2}];Number of Events", nMRBins_MultiJet, MRBinLowEdges_MultiJet);
    histRsqAllBkg =  new TH1D( "RsqAllBkg", ";M_{R} [GeV/c^{2}];Number of Events", nRsqBins_MultiJet, RsqBinLowEdges_MultiJet);
    histMRAllBkgCorr =  new TH1D( "MRAllBkgCorr",";M_{R} [GeV/c^{2}];Number of Events", nMRBins_MultiJet, MRBinLowEdges_MultiJet);
    histRsqAllBkgCorr =  new TH1D( "RsqAllBkgCorr", ";M_{R} [GeV/c^{2}];Number of Events", nRsqBins_MultiJet, RsqBinLowEdges_MultiJet);
  } else {
    histMRAllBkg =  new TH1D( "MRAllBkg",";M_{R} [GeV/c^{2}];Number of Events", nMRBins_LeptonMultiJet, MRBinLowEdges_LeptonMultiJet);
    histRsqAllBkg =  new TH1D( "RsqAllBkg", ";M_{R} [GeV/c^{2}];Number of Events", nRsqBins_LeptonMultiJet, RsqBinLowEdges_LeptonMultiJet);
    histMRAllBkgCorr =  new TH1D( "MRAllBkgCorr",";M_{R} [GeV/c^{2}];Number of Events", nMRBins_LeptonMultiJet, MRBinLowEdges_LeptonMultiJet);
    histRsqAllBkgCorr =  new TH1D( "RsqAllBkgCorr", ";M_{R} [GeV/c^{2}];Number of Events", nRsqBins_LeptonMultiJet, RsqBinLowEdges_LeptonMultiJet);
  }
  histMRAllBkg->SetStats(false);
  histRsqAllBkg->SetStats(false);
  histMRAllBkgCorr->SetStats(false);
  histRsqAllBkgCorr->SetStats(false);
  histMRAllBkg->Sumw2();
  histRsqAllBkg->Sumw2();
  histMRAllBkgCorr->Sumw2();
  histRsqAllBkgCorr->Sumw2();
 

  vector<TH1D*> histMR;
  vector<TH1D*> histRsq; 
  vector<TH1D*> histDPhiRazor;

  assert (inputfiles.size() == processLabels.size());
  for (int i=0; i < int(inputfiles.size()); ++i) {    
    histMR.push_back( new TH1D( Form("MR_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}];Number of Events", 50, 400, 2400));
    histMR[i]->SetStats(false);    
    histMR[i]->Sumw2();

    histRsq.push_back( new TH1D( Form("Rsq_%s",processLabels[i].c_str()), ";R^{2} ;Number of Events", 25, 0.25, 1.5));
    histRsq[i]->SetLineColor(color[i]);
    histRsq[i]->SetStats(false);     

    histDPhiRazor.push_back( new TH1D( Form("DPhiRazor_%s",processLabels[i].c_str()), ";#Delta#phi Hemispheres ;Number of Events", 50, 0, 3.14));
    histDPhiRazor[i]->SetLineColor(color[i]);
    histDPhiRazor[i]->SetStats(false); 
 }

  //*******************************************************************************************
  //Define Counts
  //*******************************************************************************************

  //*******************************************************************************************
  //Read files
  //*******************************************************************************************
  for (uint i=0; i < inputfiles.size(); ++i) {
    for (uint j=0; j < inputfiles[i].size(); ++j) {
      
      bool isData = false; 
      bool isSignal = false;
      if ( processLabels[i] == "Data") isData = true;
      if ( processLabels[i] == "Signal") isSignal = true;
      
      TFile* inputFile = new TFile(inputfiles[i][j].c_str(),"READ");
      assert(inputFile);
      TTree* tree = 0;
      tree = (TTree*)inputFile->Get("RazorInclusive");
      assert(tree);    
      cout << "process: " << processLabels[i] << " | file " << inputfiles[i][j] << " | Total Entries: " << tree->GetEntries() << "\n";

      float weight = 0;
      int box = -1;
      int nBTaggedJets = 0;
      float dPhiRazor = 0;
      float MR = 0;
      float Rsq = 0;
      float mT = 0;
      float mTLoose = 0;

      tree->SetBranchAddress("weight",&weight);
      tree->SetBranchAddress("box",&box);
      tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
      tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
      tree->SetBranchAddress("MR",&MR);
      tree->SetBranchAddress("Rsq",&Rsq);
      tree->SetBranchAddress("mT",&mT);    
      tree->SetBranchAddress("mTLoose",&mTLoose);    

      cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";
      for (int n=0;n<tree->GetEntries();n++) { 
    
	tree->GetEntry(n);
	if (n % 1000000 == 0) cout << "Processing Event " << n << "\n";       


	if (intLumi*weight > 100) continue;

	//Box Options
	if (btagOption == 0 ) {
	  if (!(nBTaggedJets == 0)) continue;
	}
	if (btagOption == 1 ) {
	  if (!(nBTaggedJets == 1)) continue;
	}
	if (btagOption == 2 ) {
	  if (!(nBTaggedJets == 2)) continue;
	}
	if (btagOption == 3 ) {
	  if (!(nBTaggedJets == 3)) continue;
	}
	
	if (boxOption == 11) {
	  if (box != 11) continue;
	} else if (boxOption == 12) {
	  if (box != 12) continue;
	} else if (boxOption == 9) {
	  if (box != 9) continue;
	} else if (boxOption == 10) {
	  if (box != 10) continue;
	} else if (boxOption == 3) {
	  if (box != 3) continue;
	} else if (boxOption == 4) {
	  if (box != 4) continue;
	} else if (boxOption == 6) {
	  if (box != 6) continue;
	} else if (boxOption == 7) {
	  if (box != 7) continue;
	} else if (boxOption == 100) { //multijet
	  if (!(box == 11 || box == 12)) continue;
	} else if (boxOption == 101) { //loose lepton multijet
	  if (!(box == 9 || box == 10)) continue;
	} else if (boxOption == 102) { //muon multijet
	  if (!(box == 3 || box == 4)) continue;
	} else if (boxOption == 103) { //ele multijet
	  if (!(box == 6 || box == 7)) continue;
	} else if (boxOption == 105) { //ele multijet
	  if (!(box == 3 || box == 4 || box == 6 || box == 7)) continue;
	}


	//MT cut for leptonic boxes
	if (box == 3 || box == 4 || box == 5 || box == 6 || box == 7 || box == 8 ) {
	  if (!(mT > 100)) continue;
	}
	if (box == 9 || box == 10 || box == 13) {
	  if (!(mTLoose > 100)) continue;
	}

	//dPhi cut for multijet boxes
	if(box == 11 || box == 12) {
	  if (fabs(dPhiRazor) >= 2.8) continue;
	}

	//Make MR and Rsq baseline cuts
	if (box == 11 || box == 12 || box == 9 || box == 10 || box == 13) {
	  if (!(MR>500 && Rsq > 0.25)) continue;
	} else if (box == 3 || box == 4 || box == 5 || box == 6 || box == 7 || box == 8 ) {
	  if (!(MR>400 && Rsq > 0.15)) continue;
	}

	//***************************************************************
	//Fill Histograms
	//***************************************************************

	//Plots for particular Rsq bins
	if (RsqBinOption == 0) 
	  if (!(Rsq < 0.30)) continue;
	if (RsqBinOption == 1) 
	  if (!(Rsq > 0.30 && Rsq < 0.40)) continue;
	if (RsqBinOption == 2) 
	  if (!(Rsq > 0.40)) continue;

	float slope = 0;
	float offset = 0;
	if (boxOption == 100) { slope = 0.00053133; offset = 500; }
	if (boxOption == 102) { slope = 0.00085458; offset = 400; }
	if (boxOption == 103) { slope = 0.00078621; offset = 400; }
	float corr = 1 + slope*(MR - offset);

	histMRAllBkg->Fill(MR, intLumi*weight);
	histRsqAllBkg->Fill(Rsq, intLumi*weight);
	histMRAllBkgCorr->Fill(MR, intLumi*weight*corr);
	histRsqAllBkgCorr->Fill(Rsq, intLumi*weight*corr);
	histDPhiRazor[i]->Fill(dPhiRazor, intLumi*weight);		
	histMR[i]->Fill(MR, intLumi*weight);
	histRsq[i]->Fill(Rsq, intLumi*weight);	  
	//cout << "Fill: " << n << " " << MR << " " << weight << " " << box << " " << dPhiRazor << "\n";
	
      }

      inputFile->Close();
      delete inputFile;
  
    }
  }
	


  //*******************************************************************************************
  //Draw Plots
  //*******************************************************************************************
  PlotDataAndStackedBkg( histMR, processLabels, color, false, "MR", Label);
  PlotDataAndStackedBkg( histRsq, processLabels, color, false, "Rsq", Label);

  //*******************************************************************************************
  //Summarize Counts
  //*******************************************************************************************
 
   //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("RazorPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDPhiRazor[i], Form("histDPhiRazor_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histMRAllBkg, "histMRAllBkg", "WriteDelete");
    file->WriteTObject(histRsqAllBkg, "histRsqAllBkg", "WriteDelete");
    file->WriteTObject(histMRAllBkgCorr, "histMRAllBkgCorr", "WriteDelete");
    file->WriteTObject(histRsqAllBkgCorr, "histRsqAllBkgCorr", "WriteDelete");
  }
  
  // file->WriteTObject(stackMR, "stackMR", "WriteDelete");
  // file->WriteTObject(stackRsq, "stackRsq", "WriteDelete");  
  // file->WriteTObject(stackDPhiRazor, "stackDPhiRazor", "WriteDelete");  

 }


void CompareRazorShapeVsBTag(int option = -1) {

  vector<string> datafiles;
  vector<vector<string> > bkgfiles;
  vector<string> processLabels;
  vector<int> colors;
   
   
  vector<string> bkgfiles_qcd;
  bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_QCD_HTBinned_1pb_weighted_RazorSkim.root");
  vector<string> bkgfiles_dy;
  bkgfiles_dy.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root");
  // vector<string> bkgfiles_ttbar;
  // bkgfiles_ttbar.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted_RazorSkim.root");
  vector<string> bkgfiles_ttbar;
  bkgfiles_ttbar.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_RazorSkim.root");
  bkgfiles_ttbar.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_RazorSkim.root");
  vector<string> bkgfiles_znunu;
  bkgfiles_znunu.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root");
  vector<string> bkgfiles_singletop;
  bkgfiles_singletop.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_ST_1pb_weighted_RazorSkim.root");
  vector<string> bkgfiles_wjets;
  bkgfiles_wjets.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root");
  vector<string> bkgfiles_other;
  bkgfiles_other.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_Other_1pb_weighted_RazorSkim.root");

   // bkgfiles.push_back(bkgfiles_qcd);
  bkgfiles.push_back(bkgfiles_dy);
    bkgfiles.push_back(bkgfiles_ttbar);
  bkgfiles.push_back(bkgfiles_singletop);
  bkgfiles.push_back(bkgfiles_wjets);
  bkgfiles.push_back(bkgfiles_other);

   // processLabels.push_back("QCD");
  processLabels.push_back("DY");
    processLabels.push_back("TTJets");  
  processLabels.push_back("SingleTop");
  processLabels.push_back("WJets");
  processLabels.push_back("Other");

   // colors.push_back(kMagenta);
  colors.push_back(kGreen+2);
    colors.push_back(kAzure+10);
  colors.push_back(kBlue);
  colors.push_back(kRed);
  colors.push_back(kOrange+1);

  double lumi = 133;
  lumi = 1264;

  if (option == 0 || option == 1) {
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 100,2,0,"Multijet_TwoBTags_Rsq0p25To0p3","Multijet 2 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 100,3,0,"Multijet_ThreeBTags_Rsq0p25To0p3","Multijet 3 b-tags");
    // RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 101,2,0,"LooseLeptonMultijet_TwoBTags_Rsq0p25To0p3","LooseLeptonMultijet 2 b-tags");
    // RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 101,3,0,"LooseLeptonMultijet_ThreeBTags_Rsq0p25To0p3","LooseLeptonMultijet 3 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 102,2,0,"MuonMultijet_TwoBTags_Rsq0p25To0p3","MuonMultijet 2 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 102,3,0,"MuonMultijet_ThreeBTags_Rsq0p25To0p3","MuonMultijet 3 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 103,2,0,"ElectronMultijet_TwoBTags_Rsq0p25To0p3","ElectronMultijet 2 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 103,3,0,"ElectronMultijet_ThreeBTags_Rsq0p25To0p3","ElectronMultijet 3 b-tags");
  }

  if (option == 0 || option == 2) {
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 100,2,1,"Multijet_TwoBTags_Rsq0p3To0p4","Multijet 2 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 100,3,1,"Multijet_ThreeBTags_Rsq0p3To0p4","Multijet 3 b-tags");
    // RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 101,2,1,"LooseLeptonMultijet_TwoBTags_Rsq0p3To0p4","LooseLeptonMultijet 2 b-tags");
    // RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 101,3,1,"LooseLeptonMultijet_ThreeBTags_Rsq0p3To0p4","LooseLeptonMultijet 3 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 102,2,1,"MuonMultijet_TwoBTags_Rsq0p3To0p4","MuonMultijet 2 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 102,3,1,"MuonMultijet_ThreeBTags_Rsq0p3To0p4","MuonMultijet 3 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 103,2,1,"ElectronMultijet_TwoBTags_Rsq0p3To0p4","ElectronMultijet 2 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 103,3,1,"ElectronMultijet_ThreeBTags_Rsq0p3To0p4","ElectronMultijet 3 b-tags");
  }


  if (option == 0 || option == 3) {
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 100,2,2,"Multijet_TwoBTags_Rsq0p4ToInf","Multijet 2 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 100,3,2,"Multijet_ThreeBTags_Rsq0p4ToInf","Multijet 3 b-tags");
    // RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 101,2,2,"LooseLeptonMultijet_TwoBTags_Rsq0p4ToInf","LooseLeptonMultijet 2 b-tags");
    // RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 101,3,2,"LooseLeptonMultijet_ThreeBTags_Rsq0p4ToInf","LooseLeptonMultijet 3 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 102,2,2,"MuonMultijet_TwoBTags_Rsq0p4ToInf","MuonMultijet 2 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 102,3,2,"MuonMultijet_ThreeBTags_Rsq0p4ToInf","MuonMultijet 3 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 103,2,2,"ElectronMultijet_TwoBTags_Rsq0p4ToInf","ElectronMultijet 2 b-tags");
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 103,3,2,"ElectronMultijet_ThreeBTags_Rsq0p4ToInf","ElectronMultijet 3 b-tags");
  }


  if (option == -1) {
    PlotTwoBTagVsThreeBTagShapes("Rsq0p25To0p3", true);
    PlotTwoBTagVsThreeBTagShapes("Rsq0p25To0p3", false);
    PlotTwoBTagVsThreeBTagShapes("Rsq0p3To0p4", true);
    PlotTwoBTagVsThreeBTagShapes("Rsq0p3To0p4", false);
    PlotTwoBTagVsThreeBTagShapes("Rsq0p4ToInf", true);
    PlotTwoBTagVsThreeBTagShapes("Rsq0p4ToInf", false);
  }

}
 
