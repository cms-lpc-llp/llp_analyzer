
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
  for (UInt_t b=1; int(b)<hist->GetXaxis()->GetNbins()+1; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=1; int(b)<hist->GetXaxis()->GetNbins()+1; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return hist;
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


void MakeQCDComparisonPlots(string label) {
  
  string Label = "";
  if (label != "") {
    Label = "_"+label;
  }
  TFile *fileMultijet = fileMultijet = new TFile(Form("RazorPlots_Multijet%s.root",Label.c_str()),"READ");

  TH1D *histMR_Multijet_AllBkg= (TH1D*)fileMultijet->Get("histMRAllBkg");
  TH1D *histRsq_Multijet_AllBkg = (TH1D*)fileMultijet->Get("histRsqAllBkg");
  TH1D *histMR_Multijet_QCD = (TH1D*)fileMultijet->Get("histMR_QCD");
  TH1D *histRsq_Multijet_QCD = (TH1D*)fileMultijet->Get("histRsq_QCD");

  TCanvas *cv = 0;
  TLegend *legend = 0;
 
  histMR_Multijet_AllBkg = NormalizeHist(histMR_Multijet_AllBkg);
  histRsq_Multijet_AllBkg = NormalizeHist(histRsq_Multijet_AllBkg);
  histMR_Multijet_QCD = NormalizeHist(histMR_Multijet_QCD);
  histRsq_Multijet_QCD = NormalizeHist(histRsq_Multijet_QCD);

  TH1D* tmp = 0;

  //*****************************************************
  // MR Plot
  //*****************************************************

  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.70,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histMR_Multijet_AllBkg, "All Bkgs", "FL");
  legend->AddEntry(histMR_Multijet_QCD, "QCD", "L");

  histMR_Multijet_AllBkg->SetLineWidth(2);
  histMR_Multijet_AllBkg->SetLineColor(kBlue);
  histMR_Multijet_AllBkg->Draw("e2");
  histMR_Multijet_AllBkg->SetFillStyle(3002);
  histMR_Multijet_AllBkg->SetFillColor(kBlue);
  tmp = (TH1D*)histMR_Multijet_AllBkg->Clone();
  tmp->SetFillStyle(0);
  tmp->Draw("histsame");
  histMR_Multijet_QCD->SetLineWidth(2);
  histMR_Multijet_QCD->SetLineColor(kRed);
  histMR_Multijet_QCD->Draw("e1same");

  histMR_Multijet_AllBkg->GetYaxis()->SetTitleOffset(0.9);
  histMR_Multijet_AllBkg->GetYaxis()->SetTitleSize(0.05);
  histMR_Multijet_AllBkg->GetYaxis()->SetTitle("Fraction of Events");
  histMR_Multijet_AllBkg->GetXaxis()->SetTitleSize(0.05);
  histMR_Multijet_AllBkg->GetXaxis()->SetTitleOffset(0.8);
  histMR_Multijet_AllBkg->GetXaxis()->SetTitle("M_{R} [GeV/c^{2}]");

  legend->Draw();
  cv->SaveAs(Form("QCDShapeComparison_MR%s.gif",Label.c_str()));
  cv->SaveAs(Form("QCDShapeComparison_MR%s.pdf",Label.c_str()));
  cv->SetLogy();
  cv->SaveAs(Form("QCDShapeComparison_MR%s_Logy.gif",Label.c_str()));
  cv->SaveAs(Form("QCDShapeComparison_MR%s_Logy.pdf",Label.c_str()));


  //*****************************************************
  // Rsq Plot
  //*****************************************************

  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.60,0.70,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(histRsq_Multijet_AllBkg, "All Bkgs", "FL");
  legend->AddEntry(histRsq_Multijet_QCD, "QCD", "L");

  histRsq_Multijet_AllBkg->SetLineWidth(2);
  histRsq_Multijet_AllBkg->SetLineColor(kBlue);
  histRsq_Multijet_AllBkg->Draw("e2");
  histRsq_Multijet_AllBkg->SetFillStyle(3002);
  histRsq_Multijet_AllBkg->SetFillColor(kBlue);
  tmp = (TH1D*)histRsq_Multijet_AllBkg->Clone();
  tmp->SetFillStyle(0);
  tmp->Draw("histsame");
  histRsq_Multijet_QCD->SetLineWidth(2);
  histRsq_Multijet_QCD->SetLineColor(kRed);
  histRsq_Multijet_QCD->Draw("e1same");

  histRsq_Multijet_AllBkg->GetYaxis()->SetTitleOffset(0.9);
  histRsq_Multijet_AllBkg->GetYaxis()->SetTitleSize(0.05);
  histRsq_Multijet_AllBkg->GetYaxis()->SetTitle("Fraction of Events");
  histRsq_Multijet_AllBkg->GetXaxis()->SetTitleSize(0.05);
  histRsq_Multijet_AllBkg->GetXaxis()->SetTitleOffset(0.8);
  histRsq_Multijet_AllBkg->GetXaxis()->SetTitle("R^{2}");

  legend->Draw();
  cv->SaveAs(Form("QCDShapeComparison_Rsq%s.gif",Label.c_str()));
  cv->SaveAs(Form("QCDShapeComparison_Rsq%s.pdf",Label.c_str()));
  cv->SetLogy();
  cv->SaveAs(Form("QCDShapeComparison_Rsq%s_Logy.gif",Label.c_str()));
  cv->SaveAs(Form("QCDShapeComparison_Rsq%s_Logy.pdf",Label.c_str()));




}



//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void RunMakeRazorPlots ( vector<vector<string> > bkgfiles, vector<string> bkgLabels, vector<int> bkgColors, double lumi, int boxOption = 0, int btagOption = -1, string label = "", string latexlabel = "") {

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
  TH1D* histMRAllBkg = 0;
  TH1D* histRsqAllBkg = 0;
  TH1D* histMTAllBkg = 0;
  if (boxOption == 100 || boxOption == 101) {
    histMRAllBkg =  new TH1D( "MRAllBkg",";M_{R} [GeV/c^{2}];Number of Events", 20, 800, 2800);
    histRsqAllBkg =  new TH1D( "RsqAllBkg", ";M_{R} [GeV/c^{2}];Number of Events", 10, 0.25, 1.25);
  } else {
    histMRAllBkg =  new TH1D( "MRAllBkg",";M_{R} [GeV/c^{2}];Number of Events", 20, 400, 2400);
    histRsqAllBkg =  new TH1D( "RsqAllBkg", ";M_{R} [GeV/c^{2}];Number of Events", 10, 0.25, 1.25);
  }
  histMTAllBkg =  new TH1D( "MTAllBkg", ";M_{T} [GeV/c^{2}];Number of Events", 100, 0, 400);
  histMRAllBkg->SetStats(false);
  histRsqAllBkg->SetStats(false);
  histMTAllBkg->SetStats(false);
  histMRAllBkg->Sumw2();
  histRsqAllBkg->Sumw2();
  histMTAllBkg->Sumw2();


  vector<TH1D*> histMR;
  vector<TH1D*> histRsq; 
  vector<TH1D*> histDPhiRazor;
  vector<TH1D*> histMT;

  assert (inputfiles.size() == processLabels.size());
  for (int i=0; i < int(inputfiles.size()); ++i) {    
    histMR.push_back( new TH1D( Form("MR_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}];Number of Events", 20, 800, 2800));
    histMR[i]->SetStats(false);    
    histMR[i]->Sumw2();

    histRsq.push_back( new TH1D( Form("Rsq_%s",processLabels[i].c_str()), ";R^{2} ;Number of Events", 10, 0.25, 1.25));
    histRsq[i]->SetLineColor(color[i]);
    histRsq[i]->SetStats(false);     

    histDPhiRazor.push_back( new TH1D( Form("DPhiRazor_%s",processLabels[i].c_str()), ";#Delta#phi Hemispheres ;Number of Events", 50, 0, 3.14));
    histDPhiRazor[i]->SetLineColor(color[i]);
    histDPhiRazor[i]->SetStats(false); 

    histMT.push_back( new TH1D( Form("MT_%s",processLabels[i].c_str()), ";m_{T} [GeV/c^{2}];Number of Events", 100, 000, 400));
    histMT[i]->SetStats(false);    
    histMT[i]->Sumw2();
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


	if (weight > 0.015) continue;
	if (weight > 0.0004) continue;

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
	  //if (!(mT > 100)) continue;
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
	  if (!(MR>400 && Rsq > 0.25)) continue;
	} else if (box == 3 || box == 4 || box == 5 || box == 6 || box == 7 || box == 8 ) {
	  if (!(MR>1000 && Rsq > 0.15 && Rsq < 0.20)) continue;
	  //if (!(MR>300 && Rsq > 0.15)) continue;
	}

	//***************************************************************
	//Fill Histograms
	//***************************************************************
	// if (!(MR>800)) continue;

	if (processLabels[i] != "QCD") {
	  histMRAllBkg->Fill(MR, intLumi*weight);
	  histRsqAllBkg->Fill(Rsq, intLumi*weight);
	  histMTAllBkg->Fill(mT, intLumi*weight);
	}
	histDPhiRazor[i]->Fill(dPhiRazor, intLumi*weight);		
	histMR[i]->Fill(MR, intLumi*weight);
	histRsq[i]->Fill(Rsq, intLumi*weight);	  
	histMT[i]->Fill(mT, intLumi*weight);
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
  PlotDataAndStackedBkg( histMT, processLabels, color, false, "MT", Label);

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
    file->WriteTObject(histMT[i], Form("histMT_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histMRAllBkg, "histMRAllBkg", "WriteDelete");
    file->WriteTObject(histRsqAllBkg, "histRsqAllBkg", "WriteDelete");
    file->WriteTObject(histMTAllBkg, "histMTAllBkg", "WriteDelete");
  }
  
  // file->WriteTObject(stackMR, "stackMR", "WriteDelete");
  // file->WriteTObject(stackRsq, "stackRsq", "WriteDelete");  
  // file->WriteTObject(stackDPhiRazor, "stackDPhiRazor", "WriteDelete");  

 }


void CompareQCDShape(int option = 0) {

  vector<string> datafiles;
  vector<vector<string> > bkgfiles;
  vector<string> processLabels;
  vector<int> colors;
   
   
  vector<string> bkgfiles_qcd;
  bkgfiles_qcd.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_QCD_HTBinned_1pb_weighted_RazorSkim.root");
  vector<string> bkgfiles_dy;
  bkgfiles_dy.push_back("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root");
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

  bkgfiles.push_back(bkgfiles_qcd);
  bkgfiles.push_back(bkgfiles_dy);
  bkgfiles.push_back(bkgfiles_ttbar);
   bkgfiles.push_back(bkgfiles_znunu);
  bkgfiles.push_back(bkgfiles_singletop);
  bkgfiles.push_back(bkgfiles_wjets);
  bkgfiles.push_back(bkgfiles_other);

  processLabels.push_back("QCD");
  processLabels.push_back("DY");
  processLabels.push_back("TTJets");  
   processLabels.push_back("ZNuNu");  
  processLabels.push_back("SingleTop");
  processLabels.push_back("WJets");
  processLabels.push_back("Other");

  colors.push_back(kMagenta);
  colors.push_back(kBlue);
  colors.push_back(kGreen+2);
   colors.push_back(kAzure+10);
  colors.push_back(kOrange+1);
  colors.push_back(kRed);
  colors.push_back(kGray);

  double lumi = 133;
  lumi = 2093;

  if (option == 1) {
    RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 100,-1,"Multijet","Multijet");
    //RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 100,0,"Multijet_0BTag","Multijet_0BTag");
    //RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 100,1,"Multijet_1BTag","Multijet_1BTag");
    //RunMakeRazorPlots(bkgfiles, processLabels,colors,lumi, 100,2,"Multijet_2BTag","Multijet_2BTag");
  }

  if (option == 0) {
    MakeQCDComparisonPlots("");
    //MakeQCDComparisonPlots("NoZNuNu");
  }

}
 
