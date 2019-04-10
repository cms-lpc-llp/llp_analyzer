
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

const Int_t NComponents = 10;
int color[NComponents] = {kRed, kGreen+2, kBlue, kViolet, kAzure+10, kBlack, kOrange+1, kGray, kBlack, kBlack};


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
void RunMakeDeltaPhiRazorPlots ( vector<string> files,vector<string> Labels, int boxOption = 0, int option = -1, string label = "", string latexlabel = "") {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  double intLumi = 4000; //in units of pb^-1
  string Label = "";
  if (label != "") Label = "_" + label;

  vector<string> inputfiles;
  vector<string> processLabels;

  assert(files.size() == Labels.size());
  for (int i=0; i < files.size(); ++i) {
     inputfiles.push_back(files[i]);
     processLabels.push_back(Labels[i]);
  }

  //*******************************************************************************************
  //Define Histograms
  //******************************************************************************************* 
  vector<TH1F*> histMR;
  vector<TH1F*> histRsq; 
  vector<TH1F*> histDPhiRazor;

  assert (inputfiles.size() == processLabels.size());
  for (int i=0; i < inputfiles.size(); ++i) {    
    histMR.push_back( new TH1F( Form("MR_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}];Number of Events", 25, 0, 3000));
    histMR[i]->SetLineWidth(3);
    histMR[i]->SetLineColor(color[i]);    
    histMR[i]->SetStats(false);    
    histMR[i]->Sumw2();

    histRsq.push_back( new TH1F( Form("Rsq_%s",processLabels[i].c_str()), ";R^{2} ;Number of Events", 50, 0, 1.5));
    histRsq[i]->SetLineWidth(3);
    histRsq[i]->SetLineColor(color[i]);
    histRsq[i]->SetStats(false);     

    histDPhiRazor.push_back( new TH1F( Form("DPhiRazor_%s",processLabels[i].c_str()), ";#Delta#phi Hemispheres ;Number of Events", 30, 0, 3.15));
    histDPhiRazor[i]->SetLineWidth(3);
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

    TFile* inputFile = new TFile(inputfiles[i].c_str(),"READ");
    assert(inputFile);
    TTree* tree = 0;
    tree = (TTree*)inputFile->Get("RazorInclusive");

 
    float weight = 0;
    int box = -1;
    int nBTaggedJets = 0;
    float dPhiRazor = 0;
    float MR = 0;
    float Rsq = 0;

    tree->SetBranchAddress("weight",&weight);
    tree->SetBranchAddress("box",&box);
    tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
    tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
    tree->SetBranchAddress("MR",&MR);
    tree->SetBranchAddress("Rsq",&Rsq);


    cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";
    for (int n=0;n<tree->GetEntries();n++) { 
    
      tree->GetEntry(n);
      if (n % 1000000 == 0) cout << "Processing Event " << n << "\n";       


      if (intLumi*weight > 20) continue;

      //Box Options


      //apply baseline cuts
      if ( Rsq > 0.25 && MR > 400
	  ) {
	histDPhiRazor[i]->Fill(dPhiRazor, intLumi*weight);	
	
	if (fabs(dPhiRazor) < 2.7) {
	  histMR[i]->Fill(MR, intLumi*weight);
	  histRsq[i]->Fill(Rsq, intLumi*weight);	  
	}
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
  //DPhiRazor
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.15,0.54,0.65,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  //Normalize Histograms
  for (Int_t i = 0 ; i < int(histDPhiRazor.size()); ++i) {
    histDPhiRazor[i] = NormalizeHist(histDPhiRazor[i]);
  }
  
  for (Int_t i = 0 ; i < int(histDPhiRazor.size()); ++i) {
    legend->AddEntry(histDPhiRazor[i],processLabels[i].c_str(), "L");
  }
  
  double max = 0;
  for (Int_t i = 0 ; i < int(histDPhiRazor.size()); ++i) {
    if (i==0) histDPhiRazor[i]->Draw("hist");
    else histDPhiRazor[i]->Draw("hist,same");

    if (histDPhiRazor[i]->GetMaximum() > max) max = histDPhiRazor[i]->GetMaximum();   
  }
  //histDPhiRazor[0]->SetMaximum(1.25*max);

  

  cv->SetLogy(true);
  histDPhiRazor[0]->SetMaximum(0.5);
  histDPhiRazor[0]->SetMinimum(0.001);
  histDPhiRazor[0]->GetYaxis()->SetTitleOffset(1.2);

  legend->Draw();
  cv->SaveAs(Form("DPhiRazor%s.gif",Label.c_str()));

 

 }


 void MakeDeltaPhiRazorPlots() {

   // string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/NewCategories/Samples/RazorInclusive_SMS-T1qqqq_2J_mGl-1400_mLSP-100_1pb_weighted.root";  
   // string signalLabel = "T1qqqq m_{G}=1400 m_{LSP}=100";
   // string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/NewCategories/Samples/RazorInclusive_SMS-T1qqqq_2J_mGl-1000_mLSP-800_1pb_weighted.root";  
   // string signalLabel = "T1qqqq m_{G}=1000 m_{LSP}=800";
   // string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/NewCategories/Samples/RazorInclusive_SMS-T1bbbb_2J_mGl-1500_mLSP-100_1pb_weighted.root";  
   // string signalLabel = "T1bbbb m_{G}=1500 m_{LSP}=100";
    string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/NewCategories/Samples/RazorInclusive_SMS-T1tttt_2J_mGl-1500_mLSP-100_1pb_weighted.root";  
    string signalLabel = "T1tttt m_{G}=1500 m_{LSP}=100";
   // string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorInclusive/NewCategories/Samples/RazorInclusive_SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola_AVE30BX50_ST_PHYS14_20bx25_1pb_weighted.root";
   // string signalLabel = "T1tttt m_{G}=1500 m_{LSP}=100";
   
   vector<string> files;
   vector<string> Labels;
      
   files.push_back("/afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorInclusive/MCReadiness062015/Backgrounds/RazorInclusive_QCD_HTBinned_1pb_weighted.root");
   files.push_back("/afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorInclusive/MCReadiness062015/Signals/RazorInclusive_SMS-T1bbbb_2J_mGl-1500_mLSP-100_1pb_weighted.root");
   files.push_back("/afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorInclusive/MCReadiness062015/Signals/RazorInclusive_SMS-T1bbbb_2J_mGl-1000_mLSP-900_1pb_weighted.root");
 
   Labels.push_back("QCD");
   Labels.push_back("T1bbbb m_{Gl}=1500, m_{LSP}=100");
   Labels.push_back("T1bbbb m_{Gl}=1000, m_{LSP}=900");
   
   RunMakeDeltaPhiRazorPlots(files,Labels,11,1,"All","All");
 
 
 }
 
