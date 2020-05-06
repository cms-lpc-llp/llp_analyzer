#include "TH1F.h"

using namespace std;

void RazorVarCombinedHist()
{
  //declaring histograms
  TH1F* MRHist_chi150phi1000 = new TH1F("MRHist_chi150phi1000", "MR distribution, chi=150, phi=1000", 100, 0, 3000);
  TH1F* MRHist_chi150phi500 = new TH1F("MRHist_chi150phi500", "MR distribution, chi=150, phi=500", 100, 0,3000);
  TH1F* MRHist_chi50phi150 = new TH1F("MRHist_chi50phi150", "MR distribution, chi=50, phi=150", 100,0,3000);
  TH1F* MRHist_chi50phi500 = new TH1F("MRHist_chi50phi500", "MR distribution, chi=50, phi=500", 100,0,3000);

  TH1F* MRHist_chi150phi1000_norm = new TH1F("MRHist_chi150phi1000_norm", "MR distribution, chi=150, phi=1000", 100, 0, 3000);
  TH1F* MRHist_chi150phi500_norm = new TH1F("MRHist_chi150phi500_norm", "MR distribution, chi=150, phi=500", 100, 0,3000);
  TH1F* MRHist_chi50phi150_norm = new TH1F("MRHist_chi50phi150_norm", "MR distribution, chi=50, phi=150", 100,0,3000);
  TH1F* MRHist_chi50phi500_norm = new TH1F("MRHist_chi50phi500_norm", "MR distribution, chi=50, phi=500", 100,0,3000);

  TH1F* R2Hist_chi150phi1000 = new TH1F("R2Hist_chi150phi1000", "R^2 distribution, chi=150, phi=1000", 100,0, 1);
  TH1F* R2Hist_chi150phi500 = new TH1F("R2Hist_chi150phi500", "R^2 distribution, chi=150, phi=500", 100, 0, 1);
  TH1F* R2Hist_chi50phi150 = new TH1F("R2Hist_chi50phi150", "R^2 distribution, chi=50, phi=150", 100,0, 1);
  TH1F* R2Hist_chi50phi500 = new TH1F("R2Hist_chi50phi500", "R^2 distribution, chi=50, phi=500", 100,0, 1);

  TH1F* R2Hist_chi150phi1000_norm = new TH1F("R2Hist_chi150phi1000_norm", "R^2 distribution, chi=150, phi=1000", 100,0, 1);
  TH1F* R2Hist_chi150phi500_norm = new TH1F("R2Hist_chi150phi500_norm", "R^2 distribution, chi=150, phi=500", 100, 0, 1);
  TH1F* R2Hist_chi50phi150_norm = new TH1F("R2Hist_chi50phi150_norm", "R^2 distribution, chi=50, phi=150", 100,0, 1);
  TH1F* R2Hist_chi50phi500_norm = new TH1F("R2Hist_chi50phi500_norm", "R^2 distribution, chi=50, phi=500", 100,0, 1);

  TH1F* metPtHist_chi150phi1000 = new TH1F("metPtHist_chi150phi1000", "metPt distribution, chi=150, phi=1000", 100,0, 1000);
  TH1F* metPtHist_chi150phi500 = new TH1F("metPtHist_chi150phi500", "metPt distribution, chi=150, phi=500", 100, 0, 1000);
  TH1F* metPtHist_chi50phi150 = new TH1F("metPtHist_chi50phi150", "metPt distribution, chi=50, phi=150", 100,0, 1000);
  TH1F* metPtHist_chi50phi500 = new TH1F("metPtHist_chi50phi500", "metPt distribution, chi=50, phi=500", 100,0, 1000);

  TH1F* metPtHist_chi150phi1000_norm = new TH1F("metPtHist_chi150phi1000_norm", "metPt distribution, chi=150, phi=1000", 100,0, 1000);
  TH1F* metPtHist_chi150phi500_norm = new TH1F("metPtHist_chi150phi500_norm", "metPt distribution, chi=150, phi=500", 100, 0, 1000);
  TH1F* metPtHist_chi50phi150_norm = new TH1F("metPtHist_chi50phi150_norm", "metPt distribution, chi=50, phi=150", 100,0, 1000);
  TH1F* metPtHist_chi50phi500_norm = new TH1F("metPtHist_chi50phi500_norm", "metPt distribution, chi=50, phi=500", 100,0, 1000);

  //output file
  TFile* rootfile1 = TFile::Open("CombinedHist.root", "recreate");
  
  //input files
  TFile* f1 = new TFile("chi150phi1000.root");
  TFile* f2 = new TFile("chi150phi500.root");
  TFile* f3 = new TFile("chi50phi150.root");
  TFile* f4 = new TFile("chi50phi500.root");

  //input trees
  TTree* t1 = (TTree*)f1->Get("MultiJet");
  TTree* t2 = (TTree*)f2->Get("MultiJet");
  TTree* t3 = (TTree*)f3->Get("MultiJet");
  TTree* t4 = (TTree*)f4->Get("MultiJet");

  //setting necessary tree variables
  float mr1, rsq1, metpt1, mr2, rsq2, metpt2, mr3, rsq3, metpt3, mr4, rsq4, metpt4;
  t1->SetBranchAddress("metPt", &metpt1);
  t1->SetBranchAddress("MR", &mr1);
  t1->SetBranchAddress("Rsq", &rsq1);

  t2->SetBranchAddress("metPt",&metpt2);
  t2->SetBranchAddress("MR", &mr2);
  t2->SetBranchAddress("Rsq", &rsq2);

  t3->SetBranchAddress("metPt",&metpt3);
  t3->SetBranchAddress("MR", &mr3);
  t3->SetBranchAddress("Rsq", &rsq3);

  t4->SetBranchAddress("metPt",&metpt4);
  t4->SetBranchAddress("MR", &mr4);
  t4->SetBranchAddress("Rsq", &rsq4);

  //number of tree entries
  int n1 = t1->GetEntries();
  int n2 = t2->GetEntries();
  int n3 = t3->GetEntries();
  int n4 = t4->GetEntries();

  //fill histograms
  for (int i = 0; i < n1; i++)
    {
      t1->GetEntry(i);
      MRHist_chi150phi1000->Fill(mr1);
      MRHist_chi150phi1000_norm->Fill(mr1);
      R2Hist_chi150phi1000->Fill(rsq1);
      R2Hist_chi150phi1000_norm->Fill(mr1);
      metPtHist_chi150phi1000->Fill(metpt1);
      metPtHist_chi150phi1000_norm->Fill(metpt1);
    }

  for (int i = 0; i < n2; i++)
    {
      t2->GetEntry(i);
      MRHist_chi150phi500->Fill(mr2);
      MRHist_chi150phi500_norm->Fill(mr2);
      R2Hist_chi150phi500->Fill(rsq2);
      R2Hist_chi150phi500_norm->Fill(rsq2);
      metPtHist_chi150phi500->Fill(metpt2);
      metPtHist_chi150phi500_norm->Fill(metpt2);
    }

  for (int i = 0; i < n3; i++)
    {
      t3->GetEntry(i);
      MRHist_chi50phi150->Fill(mr3);
      MRHist_chi50phi150_norm->Fill(mr3);
      R2Hist_chi50phi150->Fill(rsq3);
      R2Hist_chi50phi150_norm->Fill(rsq3);
      metPtHist_chi50phi150->Fill(metpt3);
      metPtHist_chi50phi150_norm->Fill(metpt3);
    }

  for (int i = 0; i < n4; i++)
    {
      t4->GetEntry(i);
      MRHist_chi50phi500->Fill(mr4);
      MRHist_chi50phi500_norm->Fill(mr4);
      R2Hist_chi50phi500->Fill(rsq4);
      R2Hist_chi50phi500_norm->Fill(mr4);
      metPtHist_chi50phi500->Fill(metpt4);
      metPtHist_chi50phi500_norm->Fill(metpt4);

    }

  //move to output file
  rootfile1->cd();
  
  //Normalizing histograms
  R2Hist_chi150phi1000_norm->Scale(1/(R2Hist_chi150phi1000->Integral("width")));
  R2Hist_chi150phi500_norm->Scale(1/(R2Hist_chi150phi500->Integral("width")));
  R2Hist_chi50phi150_norm->Scale(1/(R2Hist_chi50phi150->Integral("width")));
  R2Hist_chi50phi500_norm->Scale(1/(R2Hist_chi50phi500->Integral("width")));
  MRHist_chi150phi1000_norm->Scale(1/(MRHist_chi150phi1000->Integral("width")));
  MRHist_chi150phi500_norm->Scale(1/(MRHist_chi150phi500->Integral("width")));
  MRHist_chi50phi150_norm->Scale(1/(MRHist_chi50phi150->Integral("width")));
  MRHist_chi50phi500_norm->Scale(1/(MRHist_chi50phi500->Integral("width")));
  metPtHist_chi150phi1000_norm->Scale(1/(metPtHist_chi150phi1000->Integral("width")));
  metPtHist_chi150phi500_norm->Scale(1/(metPtHist_chi150phi500->Integral("width")));
  metPtHist_chi50phi150_norm->Scale(1/(metPtHist_chi50phi150->Integral("width")));
  metPtHist_chi50phi500_norm->Scale(1/(metPtHist_chi50phi500->Integral("width")));

  //write histograms to output file
  R2Hist_chi150phi1000->Write();
  R2Hist_chi150phi500->Write();
  R2Hist_chi50phi150->Write();
  R2Hist_chi50phi500->Write();
  MRHist_chi150phi1000->Write();
  MRHist_chi150phi500->Write();
  MRHist_chi50phi150->Write();
  MRHist_chi50phi500->Write();
  metPtHist_chi150phi1000->Write();
  metPtHist_chi150phi500->Write();
  metPtHist_chi50phi150->Write();
  metPtHist_chi50phi500->Write();

  R2Hist_chi150phi1000_norm->Write();
  R2Hist_chi150phi500_norm->Write();
  R2Hist_chi50phi150_norm->Write();
  R2Hist_chi50phi500_norm->Write();
  MRHist_chi150phi1000_norm->Write();
  MRHist_chi150phi500_norm->Write();
  MRHist_chi50phi150_norm->Write();
  MRHist_chi50phi500_norm->Write();
  metPtHist_chi150phi1000_norm->Write();
  metPtHist_chi150phi500_norm->Write();
  metPtHist_chi50phi150_norm->Write();
  metPtHist_chi50phi500_norm->Write();
  
  rootfile1->Close();
}
