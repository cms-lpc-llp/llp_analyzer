//Macro to compare the MR and Rsq distributions of ZNuNu and DYJets MC and save the DY/ZNuNu scale factors in bins of MR and Rsq

#include <iostream>
#include <map>
#include <string>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TTreeFormula.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPad.h"
#include "assert.h"
#include "math.h"

using namespace std;

void DrawDataVsMCRatioPlot(TH1F *dataHist, THStack *mcStack, TLegend *leg, string xaxisTitle, string printString, bool logX);

void ZNuNuWeighting(){
    gROOT->SetBatch();

    //set color palette 
    const Int_t NCont = 100;
    gStyle->SetNumberContours(NCont);
    gStyle->SetPaintTextFormat("1.0f");

    //for plots
    float nMRBins = 8;
    float nRsqBins = 7;
    float MRBinLowEdges[] = {300, 350, 400, 450, 550, 700, 900, 1200, 4000};
    float RsqBinLowEdges[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 1.5};

    map<string, string> suffixes;
    suffixes["DYJets"] = "_noZ";
    suffixes["ZJets"] = "";

    map<string, string> cuts;
    cuts["DYJets"] = "nLooseMuons == 2 && hlt_dimuon && recoZmass > 71 && recoZmass < 111 && MR_noZ > 300 && Rsq_noZ > 0.15 && numJets80_noZ > 1";
    cuts["ZJets"] = "nLooseMuons == 0 && nLooseElectrons == 0 && hlt_razor && MR > 300 && Rsq > 0.15 && numJets80 > 1";

    //get input files -- assumes one TFile for each process, with weights for different HT bins 
    map<string, TFile*> mcfiles;
    mcfiles["DYJets"] = new TFile("/afs/cern.ch/work/a/apresyan/public/DYJetsRun1_19700pb_weighted.root");
    mcfiles["ZJets"] = new TFile("/afs/cern.ch/work/a/apresyan/public/ZJetsRun1_19700pb_weighted.root");

    //get trees and set branches
    map<string, TTree*> mctrees;
    float weight;
    float MR, Rsq;
    for(auto &file : mcfiles){
        mctrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        mctrees[file.first]->SetBranchStatus("*", 0);
        mctrees[file.first]->SetBranchStatus("met_noZ", 1);
        mctrees[file.first]->SetBranchStatus("MR_noZ", 1);
        mctrees[file.first]->SetBranchStatus("Rsq_noZ", 1);
        mctrees[file.first]->SetBranchStatus("numJets80_noZ", 1);
        mctrees[file.first]->SetBranchStatus("met", 1);
        mctrees[file.first]->SetBranchStatus("MR", 1);
        mctrees[file.first]->SetBranchStatus("Rsq", 1);
        mctrees[file.first]->SetBranchStatus("numJets80", 1);
        mctrees[file.first]->SetBranchStatus("weight", 1);
        mctrees[file.first]->SetBranchStatus("recoZmass", 1);
        mctrees[file.first]->SetBranchStatus("nLooseMuons", 1);
        mctrees[file.first]->SetBranchStatus("nLooseElectrons", 1);
        mctrees[file.first]->SetBranchStatus("hlt_dimuon", 1);
        mctrees[file.first]->SetBranchStatus("hlt_razor", 1);
        mctrees[file.first]->SetBranchAddress("weight", &weight);
    }
    mctrees["DYJets"]->SetBranchAddress("MR_noZ", &MR);
    mctrees["DYJets"]->SetBranchAddress("Rsq_noZ", &Rsq);
    mctrees["ZJets"]->SetBranchAddress("MR", &MR);
    mctrees["ZJets"]->SetBranchAddress("Rsq", &Rsq);

    //Get the MR and Rsq distributions for ZJetsNuNu and DYJets
    //map<string, TH1F> metHistosForReweighting;
    map<string, TH2F> razorHistosForReweighting;
    for(auto &tree : mctrees){
        cout << "Filling MC histograms: " << tree.first << endl;
        //metHistosForReweighting[tree.first] = TH1F(Form("metmc%s", tree.first.c_str()), "MET (GeV); MET(GeV)", nMetBins, MetMin, MetMax);
        razorHistosForReweighting[tree.first] = TH2F(Form("razormc%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
        razorHistosForReweighting[tree.first].Sumw2();
        uint nEntries = tree.second->GetEntries();
        //make TTreeFormula for selection cuts
        TTreeFormula cutsFormula(Form("%sCutsFormula", tree.first.c_str()), cuts[tree.first].c_str(), tree.second);
        cutsFormula.GetNdata();
        //loop over entries
        for(uint i = 0; i < nEntries; i++){
            //get entry
            tree.second->GetEntry(i); 

            //apply selection cuts
            bool passesSelection = cutsFormula.EvalInstance();
            if(!passesSelection) continue;

            float eventWeight = weight;

            //fill each quantity
            razorHistosForReweighting[tree.first].Fill(MR, Rsq, eventWeight);
        }
        //normalize
        razorHistosForReweighting[tree.first].Scale(1.0/razorHistosForReweighting[tree.first].Integral());
    }

    TCanvas c("c", "c", 800, 600);
    c.SetLogx();
    c.SetLogz();
    //print the individual razor distributions
    for(auto &hist : razorHistosForReweighting){
        hist.second.SetTitle(Form("MC for %s", hist.first.c_str()));
        hist.second.GetXaxis()->SetTitle("MR");
        hist.second.GetYaxis()->SetTitle("Rsq");
        hist.second.SetStats(0);
        hist.second.Draw("colz");
        hist.second.Draw("same,text");
        c.Print(Form("ZNuNuWeighting%s.pdf", hist.first.c_str()));
    }

    gStyle->SetPaintTextFormat("1.2f");
    c.SetLogy(false);
    c.SetLogz(false);
    c.SetLogx();

    //get the DY/Z scale factors in each bin of MR and Rsq
    TH2F DYOverZScaleFactors = razorHistosForReweighting["DYJets"];
    DYOverZScaleFactors.Divide(&razorHistosForReweighting["ZJets"]);

    //write out the DY/Z scale factors
    c.SetLogx(true);
    c.SetLogy(true);
    c.SetLogz(false);
    TFile sfFile("ZNuNuToDYScaleFactorsRun1.root", "RECREATE");
    sfFile.cd();
    //plot it
    DYOverZScaleFactors.SetStats(0);
    DYOverZScaleFactors.SetTitle("DY/ZNuNu scale factors");
    DYOverZScaleFactors.SetMinimum(0.0);
    DYOverZScaleFactors.SetMaximum(2.0);
    DYOverZScaleFactors.Draw("colz");
    DYOverZScaleFactors.Draw("same,text");
    c.Print("ZNuNuToDYScaleFactors.pdf");
    //write it 
    DYOverZScaleFactors.Write(); 
}

int main(){
    ZNuNuWeighting();
    return 0;
}

void DrawDataVsMCRatioPlot(TH1F *dataHist, THStack *mcStack, TLegend *leg, string xaxisTitle, string printString, bool logX){
    TCanvas c("c", "c", 800, 600);
    c.Clear();
    c.cd();
    TPad pad1("pad1","pad1",0,0.4,1,1);
    pad1.SetBottomMargin(0);
    pad1.SetLogy();
    if(logX) pad1.SetLogx();
    pad1.Draw();
    pad1.cd();
    mcStack->Draw("hist");
    mcStack->GetYaxis()->SetTitle("Number of events in 19.7/fb");
    mcStack->GetYaxis()->SetLabelSize(0.03);
    mcStack->GetYaxis()->SetTitleOffset(0.45);
    mcStack->GetYaxis()->SetTitleSize(0.05);
    mcStack->SetMinimum(0.1);
    dataHist->SetMarkerStyle(20);
    dataHist->SetMarkerSize(1);
    dataHist->GetYaxis()->SetTitle("Number of events in 19.7/fb");
    dataHist->Draw("pesame");
    pad1.Modified();
    gPad->Update();
    //make ratio histogram
    TList * histList = (TList*)mcStack->GetHists();
    TIter next(histList);
    TH1 *mcTotal = (TH1*) histList->First()->Clone();
    //mcTotal->Sumw2();
    TObject *obj;
    while((obj = next())){
        if(obj == histList->First()) continue;
        mcTotal->Add((TH1*)obj);
    }
    TH1F *dataOverMC = (TH1F*)dataHist->Clone();

    string histoName = dataHist->GetName() ;
    if(histoName.find("mr") != std::string::npos  )
      {
	cout<<"Number of events in data: "<<dataHist->Integral()<<" "<<printString<<endl;
	cout<<"Number of events in MC: "<<mcTotal->Integral()<<" "<<endl;
      }
    
    //dataOverMC->Sumw2();
    dataOverMC->Divide(mcTotal);
    dataOverMC->GetXaxis()->SetTitle(xaxisTitle.c_str());
    dataOverMC->GetYaxis()->SetTitle("Data / MC");
    dataOverMC->SetMinimum(0.5);
    dataOverMC->SetMaximum(1.5);
    dataOverMC->GetXaxis()->SetLabelSize(0.1);
    dataOverMC->GetYaxis()->SetLabelSize(0.08);
    dataOverMC->GetYaxis()->SetTitleOffset(0.35);
    dataOverMC->GetXaxis()->SetTitleOffset(1.00);
    dataOverMC->GetYaxis()->SetTitleSize(0.08);
    dataOverMC->GetXaxis()->SetTitleSize(0.08);
    dataOverMC->SetStats(0);
    leg->Draw();
    c.cd();
    TPad pad2("pad2","pad2",0,0.0,1,0.4);
    pad2.SetTopMargin(0);
    pad2.SetTopMargin(0.008);
    pad2.SetBottomMargin(0.25);
    pad2.SetGridy();
    if(logX) pad2.SetLogx();
    pad2.Draw();
    pad2.cd();
    dataOverMC->Draw("pe");
    pad2.Modified();
    gPad->Update();
    c.Print(Form("%s_noISR.pdf", printString.c_str()));
    // c.Print(Form("%s.root", printString.c_str()));
}
