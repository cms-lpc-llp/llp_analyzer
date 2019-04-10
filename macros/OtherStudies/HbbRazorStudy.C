
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TRandom3.h>
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


//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void RunMakeRazorPlots ( string signalfile, string signalLabel,  vector<string> bkgfiles,vector<string> bkgLabels, int boxOption = 0, int option = -1, string label = "", string latexlabel = "") {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TRandom3 *random = new TRandom3(0);
  double intLumi = 4000; //in units of pb^-1
  string Label = "";
  if (label != "") Label = "_" + label;

  vector<string> inputfiles;
  vector<string> processLabels;

  bool hasSignal = false;
  if (signalfile != "") {
    hasSignal = true;
    inputfiles.push_back(signalfile);
    processLabels.push_back(signalLabel);
  } else {
    hasSignal = true;
    inputfiles.push_back("");
    processLabels.push_back("Hypothetical Signal");    
  }
  assert(bkgfiles.size() == bkgLabels.size());
  for (int i=0; i < bkgfiles.size(); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processLabels.push_back(bkgLabels[i]);
  }


  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  vector<double> EventCount;
  vector<double> EventCountErrSqr;

  vector<TH1F*> histMR;
  vector<TH1F*> histRsq;
  vector<TH1F*> histMbb; 
  vector<TH1F*> histDPhiRazor;

  assert (inputfiles.size() == processLabels.size());
  for (int i=0; i < inputfiles.size(); ++i) {    
    EventCount.push_back(0);
    EventCountErrSqr.push_back(0);

    histMR.push_back( new TH1F( Form("MR_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}];Number of Events", 25, 0, 3000));
    if (!hasSignal || i != 0) histMR[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histMR[i]->SetLineWidth(3);
    histMR[i]->SetLineColor(color[i]);    
    histMR[i]->SetStats(false);    
    histMR[i]->Sumw2();

    histRsq.push_back( new TH1F( Form("Rsq_%s",processLabels[i].c_str()), ";R^{2} ;Number of Events", 50, 0, 1.5));
    if (!hasSignal || i != 0) histRsq[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histRsq[i]->SetLineWidth(3);
    histRsq[i]->SetLineColor(color[i]);
    histRsq[i]->SetStats(false);     

    histMbb.push_back( new TH1F( Form("Mbb_%s",processLabels[i].c_str()), ";M_{bb} [GeV/c^{2}];Number of Events", 15, 60, 360));
    if (!hasSignal || i != 0) histMbb[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histMbb[i]->SetLineWidth(3);
    histMbb[i]->SetLineColor(color[i]);    
    histMbb[i]->SetStats(false);    
    histMbb[i]->Sumw2();

    histDPhiRazor.push_back( new TH1F( Form("DPhiRazor_%s",processLabels[i].c_str()), ";#Delta#phi Hemispheres ;Number of Events", 50, 0, 3.14));
    if (!hasSignal || i != 0) histDPhiRazor[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histDPhiRazor[i]->SetLineWidth(3);
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

    if (hasSignal && i == 0 && inputfiles[i] == "") continue;

    TFile* inputFile = new TFile(inputfiles[i].c_str(),"READ");
    assert(inputFile);
    TTree* tree = 0;
    tree = (TTree*)inputFile->Get("HbbRazor");
    assert(tree);

    float w = 0;
    int box = -1;
    int nBTaggedJets = 0;
    float dPhiRazor = 0;
    float MR = 0;
    float mbb = 0;
    float ptbb = 0;
    float Rsq = 0;

    tree->SetBranchAddress("weight",&w);
    tree->SetBranchAddress("box",&box);
    tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
    tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
    tree->SetBranchAddress("MR",&MR);
    tree->SetBranchAddress("mbb",&mbb);
    tree->SetBranchAddress("ptbb",&ptbb);
    tree->SetBranchAddress("Rsq",&Rsq);

  cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";
    for (int n=0;n<tree->GetEntries();n++) { 
    
      tree->GetEntry(n);
      if (n % 1000000 == 0) cout << "Processing Event " << n << "\n";       



      if (intLumi*w > 100) continue;
      double weight = w;
      if (processLabels[i] == "QCD") weight = weight * 4;            

      //Box Options
      if (option == 0 ) {
      }

      //apply baseline cuts
      if (!(MR > 200 && Rsq > 0.02)) continue;


      if ( Rsq > 0.035 && MR > 350 
	   && ptbb < 110  
	   ) {
	histDPhiRazor[i]->Fill(dPhiRazor, intLumi*weight);		
	histMbb[i]->Fill(mbb, intLumi*weight);
	if (mbb > 115 && mbb < 135) {
	  EventCount[i] += intLumi*weight;
	  EventCountErrSqr[i] += pow(intLumi*weight,2);	
	  histRsq[i]->Fill(Rsq, intLumi*weight);	  
	  histMR[i]->Fill(MR, intLumi*weight);
	}
      }
  
    }

    inputFile->Close();
    delete inputFile;
  
  }
  
  
  //*******************************************************************************************
  //Fill Hypothetical Signal
  //*******************************************************************************************
  double NSignal = 1900 ;//* (0.6/0.5) * (0.6/0.5);
  if (hasSignal && inputfiles[0] == "") {
    for (int i=0; i<10000; ++i) {
      double m = random->Gaus(125,10);
      double mr = random->Landau(400,50);
      histMbb[0]->Fill(m, NSignal/10000.0);
      if (m > 115 && m < 135) {
	EventCount[0] += NSignal/10000.0;
	histMR[0]->Fill(mr, NSignal/10000.0);
      }
    }
  }


  //*******************************************************************************************
  //Draw Plots
  //*******************************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;
  bool firstdrawn = false;
  TLatex *tex = 0;

  //*******************************************************************************************
  //MR
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.50,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMR = new THStack();

  if (hasSignal) {
    for (Int_t i = histMR.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histMR[i]->GetNbinsX()+1; j++) {
	intError += histMR[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histMR[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histMR[i]->Integral() > 0) {
	stackMR->Add(histMR[i]);
      }
    }    
  } else {
    for (Int_t i = histMR.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histMR[i]->GetNbinsX()+1; j++) {
	intError += histMR[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histMR[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histMR[i]->Integral() > 0) {
	stackMR->Add(histMR[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histMR.size()); ++i) {
    // if (hasSignal && i==0) {
    //   legend->AddEntry(histMR[i],processLabels[i].c_str(), "L");
    // } else {
      legend->AddEntry(histMR[i],processLabels[i].c_str(), "F");
    // }
  }
  

  stackMR->Draw("hist");
  stackMR->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackMR->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetYaxis()->GetTitle());

  if (hasSignal) {
    histMR[0]->Draw("histsame");
    cout << processLabels[0] << " : " << histMR[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 13 TeV, #int L = %d fb^{-1}, %s",int(intLumi/1000), latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("MR%s.gif",Label.c_str()));

 


  //*******************************************************************************************
  //Mbb
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.50,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stackMbb = new THStack();

  if (hasSignal) {
    for (Int_t i = histMbb.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histMbb[i]->GetNbinsX()+1; j++) {
  	intError += histMbb[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histMbb[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histMbb[i]->Integral() > 0) {
	if (i!=0) histMbb[i]->Smooth();
  	stackMbb->Add(histMbb[i]);
      }
    }    
  } else {
    for (Int_t i = histMbb.size()-1; i >= 0; i--) {
      double intError = 0;
      for(int j=1; j < histMbb[i]->GetNbinsX()+1; j++) {
  	intError += histMbb[i]->GetBinError(j);
      }
      cout << processLabels[i] << " : " << histMbb[i]->GetSumOfWeights() << " +/- " << intError << "\n";
      if ( histMbb[i]->Integral() > 0) {
  	stackMbb->Add(histMbb[i]);
      }
    }
  }
  for (Int_t i = 0 ; i < int(histMbb.size()); ++i) {
    // if (hasSignal && i==0) {
    //   legend->AddEntry(histMbb[i],processLabels[i].c_str(), "L");
    // } else {
      legend->AddEntry(histMbb[i],processLabels[i].c_str(), "F");
    // }
  }
  

  stackMbb->Draw("hist");
  stackMbb->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMbb->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackMbb->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMbb->GetHists()->At(0)))->GetYaxis()->GetTitle());
  stackMbb->GetHistogram()->SetLineColor(kBlack);

  if (hasSignal) {
    histMbb[0]->Draw("histsame");
    cout << processLabels[0] << " : " << histMbb[0]->GetSumOfWeights() << "\n";
  }

  legend->Draw();

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 13 TeV, #int L = %d fb^{-1}, %s",int(intLumi/1000), latexlabel.c_str()));
  tex->Draw();
  
  cv->SaveAs(Form("Mbb%s.gif",Label.c_str()));

 

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
    file->WriteTObject(histMbb[i], Form("histMbb_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histDPhiRazor[i], Form("histDPhiRazor_%s",processLabels[i].c_str()), "WriteDelete");
  }
  
  file->WriteTObject(stackMR, "stackMR", "WriteDelete");
  // file->WriteTObject(stackRsq, "stackRsq", "WriteDelete");  
  // file->WriteTObject(stackDPhiRazor, "stackDPhiRazor", "WriteDelete");  

 }


 void HbbRazorStudy() {

   string signalfile = "";
   string signalLabel = "";

   // string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/MiniIso/RazorInclusive_SMS-T1qqqq_2J_mGl-1400_mLSP-100_1pb_weighted.root";  
   // string signalLabel = "T1qqqq m_{G}=1400 m_{LSP}=100";
   // string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/MiniIso/RazorInclusive_SMS-T1bbbb_2J_mGl-1500_mLSP-100_1pb_weighted.root";  
   // string signalLabel = "T1bbbb m_{G}=1500 m_{LSP}=100";
   //string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/MiniIso/RazorInclusive_SMS-T1tttt_2J_mGl-1500_mLSP-100_1pb_weighted.root";  
   //string signalLabel = "T1tttt m_{G}=1500 m_{LSP}=100";
   
   vector<string> bkgfiles;
   vector<string> bkgLabels;
   
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVT/HbbRazor_TTJets_1pb_weighted.root");  
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVT/HbbRazor_DYJetsToLL_1pb_weighted.root");
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVT/HbbRazor_WJetsToLNu_1pb_weighted.root");
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVT/HbbRazor_ZJetsToNuNu_1pb_weighted.root");
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVT/HbbRazor_QCD_1pb_weighted.root"); 
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVT/HbbRazor_SingleTop_1pb_weighted.root"); 
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/HbbRazor/CSVT/HbbRazor_Multiboson_1pb_weighted.root"); 
   
   bkgLabels.push_back("TTJets");
   bkgLabels.push_back("DYJetsToLL");
   bkgLabels.push_back("WJetsToLNu");
   bkgLabels.push_back("ZJetsToNuNu");
   bkgLabels.push_back("QCD");
   bkgLabels.push_back("SingleTop");
   bkgLabels.push_back("Other");

   RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,0,1,"MR350Rsq0p035PT110_CSVT","H#rightarrowbb Razor");
   
 
 }
 
