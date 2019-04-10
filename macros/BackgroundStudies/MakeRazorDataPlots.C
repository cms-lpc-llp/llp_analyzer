
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
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



void PlotDataAndStackedBkg( vector<TH1D*> hist , vector<string> processType, vector<string> processLabels, vector<int> color,  bool hasData, int intLumi, string varName, string label ) {

  //*******************************************************************************************
  //Draw Plots
  //*******************************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;
  TLatex *tex = 0;

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

  legend = new TLegend(0.50,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  //*******************************************************************************************
  //Top Pad
  //*******************************************************************************************
  TPad *pad1 = new TPad("pad1","pad1", 0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(0.04);
  pad1->Draw();
  pad1->cd();

  THStack *stack = new THStack();
  TH1D *dataHist = 0;
  
  for (Int_t i = hist.size()-1; i >= 0; i--) {
    if (processType[i] == "bkg") {
      double intError = 0;
      for(int j=1; j < hist[i]->GetNbinsX()+1; j++) {
	intError += hist[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << hist[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( hist[i]->Integral() > 0) {
	stack->Add(hist[i]);
      }
    }    
  }

  for (Int_t i = 0 ; i < int(hist.size()); ++i) {
    if (processType[i] == "signal") {
      legend->AddEntry(hist[i],processLabels[i].c_str(), "L");
    } else if (processType[i] == "data")  {
      dataHist = hist[i];
      legend->AddEntry(hist[i],processLabels[i].c_str(), "LP");
    } else {
      legend->AddEntry(hist[i],processLabels[i].c_str(), "F");
    }
  }
  
  stack->Draw("hist");
  stack->GetHistogram()->GetXaxis()->SetTitle(((TH1D*)(stack->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stack->GetHistogram()->GetYaxis()->SetTitle(((TH1D*)(stack->GetHists()->At(0)))->GetYaxis()->GetTitle());
  stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.0);
  stack->GetHistogram()->GetYaxis()->SetTitleSize(0.05);

  double max = fmax(stack->GetHistogram()->GetMaximum() , dataHist->GetMaximum());
  stack->SetMaximum(1.25*max);

  for (Int_t i = 0 ; i < int(hist.size()); ++i) {
    if (processType[i] == "data") {
      hist[i]->SetMarkerStyle(20);      
      hist[i]->SetMarkerSize(1);
      hist[i]->Draw("pesame");
      cout << processLabels[i] << " : " << hist[i]->GetSumOfWeights() << "\n";
    } else if (processType[i] == "signal") {
      hist[i]->Draw("histsame");
      cout << processLabels[i] << " : " << hist[i]->GetSumOfWeights() << "\n";
    }
  }

  legend->Draw();

  //****************************
  //Add CMS and Lumi Labels
  //****************************
  lumi_13TeV = "42 pb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  CMS_lumi(pad1,4,0);

  cv->cd();
  cv->Update();



  //*******************************************************************************************
  //Bottom Pad
  //*******************************************************************************************
  TPad *pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  TH1F *histDataOverMC = (TH1F*)dataHist->Clone("histDataOverMC");
  for (int b=0; b<histDataOverMC->GetXaxis()->GetNbins()+2; ++b) {
    double data = 0;
    if (hasData) {
      data = dataHist->GetBinContent(b);
    }
    double MC = 0;
    double MCErrSqr = 0;
    if (hasData) {
      for (Int_t i = 0 ; i < int(hist.size()); ++i) {
	if (processType[i] == "bkg") {
	  MC += hist[i]->GetBinContent(b);
	  MCErrSqr += pow(hist[i]->GetBinError(b),2);
	}
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
  histDataOverMC->SetStats(false);
  histDataOverMC->Draw("pe");


  pad1->SetLogy(false);
  cv->SaveAs(Form("RazorPlot_%s%s.png",varName.c_str(), label.c_str()));
  cv->SaveAs(Form("RazorPlot_%s%s.pdf",varName.c_str(), label.c_str()));
  
  pad1->SetLogy(true);
  stack->SetMinimum(0.1);  
  cv->SaveAs(Form("RazorPlot_%s%s_Logy.png",varName.c_str(),label.c_str()));
  cv->SaveAs(Form("RazorPlot_%s%s_Logy.pdf",varName.c_str(),label.c_str()));

}

//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void RunMakeRazorDataPlots ( string datafile, string dataLabel,  
			     vector<string> bkgfiles,vector<string> bkgLabels,  vector<int> bkgColors, 
			     string signalfile, string signalLabel, 
			     int boxOption = 0, int option = -1, string label = "", string latexlabel = "") {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TFile *NVtxWeightFile = new TFile("/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/data/NVtxReweight_DiJet.root", "READ");
  TH1D *NVtxWeightHist = (TH1D*)NVtxWeightFile->Get("NVtxReweight");
  assert(NVtxWeightHist);


  double intLumi = 40.0 * (1.0 / (8.126350e-01*2.473122e+02)); //in units of pb^-1
  string Label = "";
  if (label != "") Label = "_" + label;

  vector<string> inputfiles;
  vector<string> processType;
  vector<string> processLabels;
  vector<int> color;

  bool hasData = false;
  if (datafile != "") {
    hasData = true;
    inputfiles.push_back(datafile);
    processType.push_back("data");
    processLabels.push_back(dataLabel);
    color.push_back(kBlack);
  }

  bool hasSignal = false;
  if (signalfile != "") {
    hasSignal = true;
    inputfiles.push_back(signalfile);
    processType.push_back("signal");
    processLabels.push_back(signalLabel);
    color.push_back(kRed);
  }
  assert(bkgfiles.size() == bkgLabels.size());
  for (int i=0; i < bkgfiles.size(); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processType.push_back("bkg");
     processLabels.push_back(bkgLabels[i]);
     color.push_back(bkgColors[i]);
  }

  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  TH1D* histNVtxData = new TH1D( "histNVtxData", ";Number of Reconstructed Primary Vertices;Number of Events", 40, -0.5,39.5);
  TH1D* histNVtxAllBkg = new TH1D( "histNVtxAllBkg", ";Number of Reconstructed Primary Vertices;Number of Events", 40, -0.5,39.5);

  vector<TH1D*> histNVtx;
  vector<TH1D*> histMR;
  vector<TH1D*> histRsq; 
  vector<TH1D*> histDPhiRazor;
  vector<TH1D*> histMET; 
  vector<TH1D*> histNBTags;
  vector<TH1D*> histNJets;

  assert (inputfiles.size() == processLabels.size());
  for (int i=0; i < inputfiles.size(); ++i) {    
    histNVtx.push_back( new TH1D( Form("NVtx_%s",processLabels[i].c_str()), ";Number of Reconstructed Primary Vertices;Number of Events", 40, -0.5,39.5));
    if (processType[i] == "bkg") histNVtx[i]->SetFillColor(color[i]);
    if (processType[i] == "data") histNVtx[i]->SetLineWidth(1);
    histNVtx[i]->SetLineColor(color[i]);    
    histNVtx[i]->SetStats(false);    
    histNVtx[i]->Sumw2();

    histMR.push_back( new TH1D( Form("MR_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}];Number of Events", 25, 0, 3000));
    if (processType[i] == "bkg") histMR[i]->SetFillColor(color[i]);
    if (processType[i] == "data") histMR[i]->SetLineWidth(1);
    histMR[i]->SetLineColor(color[i]);    
    histMR[i]->SetStats(false);    
    histMR[i]->Sumw2();

    histRsq.push_back( new TH1D( Form("Rsq_%s",processLabels[i].c_str()), ";R^{2} ;Number of Events", 50, 0, 0.5));
    if (processType[i] == "bkg") histRsq[i]->SetFillColor(color[i]);
    if (processType[i] == "data") histRsq[i]->SetLineWidth(1);
    histRsq[i]->SetLineColor(color[i]);
    histRsq[i]->SetStats(false);     

    histDPhiRazor.push_back( new TH1D( Form("DPhiRazor_%s",processLabels[i].c_str()), ";#Delta#phi Hemispheres ;Number of Events", 50, 0, 3.14));
    if (processType[i] == "bkg") histDPhiRazor[i]->SetFillColor(color[i]);
    if (processType[i] == "data") histDPhiRazor[i]->SetLineWidth(1);
    histDPhiRazor[i]->SetLineColor(color[i]);
    histDPhiRazor[i]->SetStats(false); 

    histMET.push_back( new TH1D( Form("MET_%s",processLabels[i].c_str()), ";MET [GeV/c^{2}];Number of Events", 100, 0, 1000));
    if (processType[i] == "bkg") histMET[i]->SetFillColor(color[i]);
    if (processType[i] == "data") histMET[i]->SetLineWidth(1);
    histMET[i]->SetLineColor(color[i]);    
    histMET[i]->SetStats(false);    
    histMET[i]->Sumw2();

   histNBTags.push_back( new TH1D( Form("NBTags_%s",processLabels[i].c_str()), "; Number of B-Tagged Jets ;Number of Events", 8, -0.5, 7.5));
    if (processType[i] == "bkg") histNBTags[i]->SetFillColor(color[i]);
    if (processType[i] == "data") histNBTags[i]->SetLineWidth(1);
    histNBTags[i]->SetLineColor(color[i]);
    histNBTags[i]->SetStats(false); 

   histNJets.push_back( new TH1D( Form("NJets_%s",processLabels[i].c_str()), "; Number of Selected Jets ;Number of Events", 15, -0.5, 14.5));
    if (processType[i] == "bkg") histNJets[i]->SetFillColor(color[i]);
    if (processType[i] == "data") histNJets[i]->SetLineWidth(1);
    histNJets[i]->SetLineColor(color[i]);
    histNJets[i]->SetStats(false); 
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
    tree = (TTree*)inputFile->Get("RazorInclusive");
 
    float weight = 0;
    int nvtx = 0;
    int box = -1;
    int nBTaggedJets = 0;
    int nSelectedJets = 0;
    int nJets80 = 0;
    float dPhiRazor = 0;
    float MR = 0;
    float Rsq = 0;
    float met = 0;
    bool  HLTDecision[150];

    tree->SetBranchAddress("weight",&weight);
    tree->SetBranchAddress("box",&box);
    tree->SetBranchAddress("nVtx",&nvtx);
    tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
    tree->SetBranchAddress("nSelectedJets",&nSelectedJets);
    tree->SetBranchAddress("nJets80",&nJets80);
    tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
    tree->SetBranchAddress("MR",&MR);
    tree->SetBranchAddress("Rsq",&Rsq);
    tree->SetBranchAddress("met",&met);
    tree->SetBranchAddress("HLTDecision",&HLTDecision);

    cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";
    for (int n=0;n<tree->GetEntries();n++) { 
    
      tree->GetEntry(n);
      if (n % 1000000 == 0) cout << "Processing Event " << n << "\n";       
      
      if (processType[i]=="data") {
	//cout << MR << " : " << HLTDecision[75] << "\n";
	if (!HLTDecision[77]) continue;
      }
      
      if (!(nJets80 >= 2)) continue;
      if (processType[i]=="data") {
	histNVtxData->Fill(nvtx);
	histNVtx[i]->Fill(nvtx);
	histMR[i]->Fill(MR);
	histRsq[i]->Fill(Rsq);
	histDPhiRazor[i]->Fill(dPhiRazor);
	histMET[i]->Fill(met);
	histNBTags[i]->Fill(nBTaggedJets);
	histNJets[i]->Fill(nSelectedJets);
      } else {

	double NVtxWeight = NVtxWeightHist->GetBinContent(NVtxWeightHist->GetXaxis()->FindFixBin(nvtx));
	double w = intLumi* weight * NVtxWeight;

	histNVtxAllBkg->Fill(nvtx, w);
	histNVtx[i]->Fill(nvtx, w);
	histMR[i]->Fill(MR, w);
	histRsq[i]->Fill(Rsq, w);
	histDPhiRazor[i]->Fill(dPhiRazor, w);
	histMET[i]->Fill(met,w);
	histNBTags[i]->Fill(nBTaggedJets, w);
	histNJets[i]->Fill(nSelectedJets, w);
      }

    }

    inputFile->Close();
    delete inputFile;
  
  }
  
  
  //*******************************************************************************************
  //Make NVtx Reweighting Function
  //*******************************************************************************************
  TH1D *NVtxDataNormalized = NormalizeHist( histNVtxData );
  TH1D *NVtxBkgNormalized = NormalizeHist( histNVtxAllBkg );

  TH1D *NVtxReweight = (TH1D*)NVtxDataNormalized->Clone("NVtxReweight");  
  for (int i=1; i<NVtxReweight->GetXaxis()->GetNbins()+1; i++) {

    double data = 0;
    double bkg = 0;
    if (NVtxBkgNormalized->GetBinContent(i) > 0) {
      NVtxReweight->SetBinContent(i,NVtxDataNormalized->GetBinContent(i)/NVtxBkgNormalized->GetBinContent(i));
    } else if (NVtxDataNormalized->GetBinContent(i) == 0){
      NVtxReweight->SetBinContent(i,0.0);
    } else {
      NVtxReweight->SetBinContent(i,1.0);
    }
  }


  //*******************************************************************************************
  //Make Plots
  //*******************************************************************************************
  PlotDataAndStackedBkg( histMR, processType, processLabels, color, true, intLumi, "MR", Label);
  PlotDataAndStackedBkg( histRsq, processType, processLabels, color, true, intLumi, "Rsq", Label);
  PlotDataAndStackedBkg( histDPhiRazor, processType, processLabels, color, true, intLumi, "DPhiRazor", Label);
  PlotDataAndStackedBkg( histMET, processType, processLabels, color, true, intLumi, "MET", Label);
  PlotDataAndStackedBkg( histNVtx, processType, processLabels, color, true, intLumi, "NVtx", Label);
  PlotDataAndStackedBkg( histNBTags, processType, processLabels, color, true, intLumi, "NBTags", Label);
  PlotDataAndStackedBkg( histNJets, processType, processLabels, color, true, intLumi, "NJets", Label);

  //*******************************************************************************************
  //Summarize Counts
  //*******************************************************************************************
 
   //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open("NVtxReweight_DiJet.root", "UPDATE");
  file->cd();
  file->WriteTObject(NVtxReweight, "NVtxReweight", "WriteDelete");
  file->Close();
  delete file;


  file = TFile::Open(("RazorPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();
  file->WriteTObject(NVtxDataNormalized, "NVtxDataNormalized", "WriteDelete");
  file->WriteTObject(NVtxBkgNormalized, "NVtxBkgNormalized", "WriteDelete");
 
  for(int i=0; i<int(inputfiles.size()); i++) {  
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDPhiRazor[i], Form("histDPhiRazor_%s",processLabels[i].c_str()), "WriteDelete");
  }
  
  // file->WriteTObject(stackMR, "stackMR", "WriteDelete");
  // file->WriteTObject(stackRsq, "stackRsq", "WriteDelete");  
  // file->WriteTObject(stackDPhiRazor, "stackDPhiRazor", "WriteDelete");  

 }


 void MakeRazorDataPlots() {

   string datafile = "/afs/cern.ch/user/s/sixie/eos/cms//store/group/phys_susy/razor/Run2Analysis/RazorInclusive/RazorInclusive_JetHT_Run2015B_GoodLumi.root";
   string dataLabel = "Data";
   string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/NewCategories/Samples/RazorInclusive_SMS-T1tttt_2J_mGl-1500_mLSP-100_1pb_weighted.root";  
   string signalLabel = "T1tttt m_{G}=1500 m_{LSP}=100";
   signalfile = "";
   signalLabel = "";


   vector<string> bkgfiles;
   vector<string> bkgLabels;
   vector<int> bkgColors;
    
   bkgfiles.push_back("/afs/cern.ch/user/s/sixie/eos/cms//store/group/phys_susy/razor/Run2Analysis/RazorInclusive/RazorInclusive_QCD_1pb_weighted.root");
   bkgLabels.push_back("QCD");
   bkgColors.push_back(kViolet);

   
   RunMakeRazorDataPlots(datafile, dataLabel,bkgfiles,bkgLabels,bkgColors,signalfile,signalLabel,11,1,"DiJet","DiJet");
 
 }
 
