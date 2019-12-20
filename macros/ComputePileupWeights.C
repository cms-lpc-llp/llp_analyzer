
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
void NormalizeHist(TH1 *hist) {
  Double_t norm = 0;
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }
}


//*************************************************************************************************
//Divide Hist
//*************************************************************************************************
void DivideHist(TH1 *ratio, TH1 *num, TH1 *den) {  
  for (UInt_t b=0; int(b)<num->GetXaxis()->GetNbins()+2; ++b) {
    if ( den->GetBinContent(b) > 1.0e-4 ) {
      cout << "Bin: " << b << " " << ratio->GetXaxis()->GetBinCenter(b) << " : " << num->GetBinContent(b) << " / " << den->GetBinContent(b) << " = " << num->GetBinContent(b) / den->GetBinContent(b) << "\n";
      ratio->SetBinContent(b,num->GetBinContent(b) / den->GetBinContent(b));    
      ratio->SetBinError(b, (num->GetBinContent(b) / den->GetBinContent(b))*sqrt( pow(num->GetBinError(b)/num->GetBinContent(b),2) + pow(den->GetBinError(b)/den->GetBinContent(b),2)));
    } else {
      ratio->SetBinContent(b,0);
      ratio->SetBinError(b,0);
    }
  }
}




void DoComputePileupWeights( string targetPileupFilename, string MCPileupFilename, string MCLabel, string outputfile , string label) {

  
  TFile *TargetPileupFile = new TFile( targetPileupFilename.c_str() , "READ");
  TH1D *TargetPileupDistribution = (TH1D*)TargetPileupFile->Get("pileup");
  NormalizeHist(TargetPileupDistribution);
  
  TFile *MCPileupFile = new TFile( MCPileupFilename.c_str() , "READ");
  TH1F *MCPileupDistribution = (TH1F*)MCPileupFile->Get(Form("PUMean_%s",MCLabel.c_str()));
  
  TH1F *PileupWeightHist = (TH1F*)MCPileupDistribution->Clone(Form("PUWeight_%s",label.c_str()));
  DivideHist( PileupWeightHist , TargetPileupDistribution, MCPileupDistribution);

  TFile *outputFile = new TFile( outputfile.c_str(), "UPDATE");
  outputFile->WriteTObject(PileupWeightHist, Form("PUWeight_%s",label.c_str()), "WriteDelete");

 }
 
void ComputePileupWeights() {

  DoComputePileupWeights( "/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/Run1DataPileupTarget.root",
			  "/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/data/MCPileupDistribution.root",
			  "DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball",
			  "Run1PileupWeights.root",
			  "Run1");

}
