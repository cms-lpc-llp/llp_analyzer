//Macro to predict the MR and Rsq distributions of the Z->nu nu background in the razor search using Z->mu mu, W->mu nu, and Gamma+Jets events

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
#include "TColor.h"
#include "assert.h"
#include "math.h"

using namespace std;

bool debug = false;
//bool debug = true;

//true: find the translation factors from MC to data
//false: find the translation factors from DY, W, G to Z->nunu
bool computeDataOverMCSFs = true;

void DrawDataVsMCRatioPlot(TH1F *dataHist, THStack *mcStack, TLegend *leg, string xaxisTitle, string printString, bool logX);

void HbbRazorPlots(){
    gROOT->SetBatch();

    //set color palette 
    const Int_t NCont = 100;
    gStyle->SetNumberContours(NCont);
    gStyle->SetPaintTextFormat("1.0f");

    //upper bounds of reweighing histograms
    float maxPhotonPt = 999; 

    //map<string, string> suffixes;
    //suffixes["EMQCD"] = "_noPho";
    //suffixes["TTJets"] = "_noPho";

    //get input files -- assumes one TFile for each process, with weights for different HT bins 
    map<string, TFile*> mcfiles;
    map<string, TFile*> datafiles;
    mcfiles["QCD"] = new TFile("./WK_7_11/BTAGCSV_2016B_PRv2_DAY_7_11/Normalized/6000pb-1/COMPLETE_QCD_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_DAY7_11_6000pb.root");
    mcfiles["TTJets"] = new TFile("./WK_7_11/BTAGCSV_2016B_PRv2_DAY_7_11/Normalized/6000pb-1/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_DAY7_11_6000pb_weighted.root");
    //mcfiles ["test"]= new TFile("./wk_6_27/Hbb/scripts/HbbRazor_ZprimeToA0hToA0chichihbb_2HDM_MZp-1000_MA0-300_13TeV.root");
     
    //datafiles["Hbb"] = new TFile("./wk_6_27/Hbb/scripts/HbbRazor_ZprimeToA0hToA0chichihbb_2HDM_MZp-1000_MA0-300_13TeV.root");
    datafiles["Hbb"] = new TFile("./WK_7_11/BTAGCSV_2016B_PRv2_DAY_7_11/DATATRUE_BTagCSV_2016B_PRv2_Complete.root");
    //get trees and set branches
    map<string, TTree*> mctrees;
    map<string, TTree*> datatrees;
    float weight;
    bool passed_DiPFJet80_DiPFJet30_BTagCSVd07d05, passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05, passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03, passed_DiJet80Eta2p6_BTagIP3DFastPVLoose, passed_QuadJet45, passed_QuadJet50, passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200, passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200;
    int nPU_mean, nVtx, nSelectedPhotons;
    float MR, Rsq, HT, met, mbb;
    int nSelectedJets;
    float jet1pt, b1pt, b2pt, ptbb, etabb, b1phi, b2phi;
    int nBTaggedJets;
    
    for(auto &file : mcfiles){
        mctrees[file.first] = (TTree*)file.second->Get("HbbRazor");

        mctrees[file.first]->SetBranchStatus("*", 0); // disable all
        mctrees[file.first]->SetBranchStatus("weight", 1);
        mctrees[file.first]->SetBranchStatus("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05", 1);
        mctrees[file.first]->SetBranchStatus("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05", 1);
        mctrees[file.first]->SetBranchStatus("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03", 1);
        mctrees[file.first]->SetBranchStatus("passed_DiJet80Eta2p6_BTagIP3DFastPVLoose", 1);
        mctrees[file.first]->SetBranchStatus("passed_QuadJet45", 1);
        mctrees[file.first]->SetBranchStatus("passed_QuadJet50", 1);
        mctrees[file.first]->SetBranchStatus("nVtx", 1); // enable 
        mctrees[file.first]->SetBranchStatus("MR", 1);
        mctrees[file.first]->SetBranchStatus("Rsq", 1);
        mctrees[file.first]->SetBranchStatus("HT", 1);
        mctrees[file.first]->SetBranchStatus("met", 1);
        mctrees[file.first]->SetBranchStatus("nSelectedJets", 1);
        mctrees[file.first]->SetBranchStatus("nBTaggedJets", 1);
        mctrees[file.first]->SetBranchStatus("jet1pt", 1);
        mctrees[file.first]->SetBranchStatus("b1pt", 1);
        mctrees[file.first]->SetBranchStatus("b2pt", 1);
        mctrees[file.first]->SetBranchStatus("b1phi", 1);
        mctrees[file.first]->SetBranchStatus("b2phi", 1);
        mctrees[file.first]->SetBranchStatus("mbb", 1);
        mctrees[file.first]->SetBranchStatus("ptbb", 1);
        mctrees[file.first]->SetBranchStatus("etabb", 1);

        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05", &passed_DiPFJet80_DiPFJet30_BTagCSVd07d05);
        mctrees[file.first]->SetBranchAddress("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05", &passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05);
        mctrees[file.first]->SetBranchAddress("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03", &passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03);
        mctrees[file.first]->SetBranchAddress("passed_DiJet80Eta2p6_BTagIP3DFastPVLoose", &passed_DiJet80Eta2p6_BTagIP3DFastPVLoose);
        mctrees[file.first]->SetBranchAddress("passed_QuadJet45", &passed_QuadJet45);
        mctrees[file.first]->SetBranchAddress("passed_QuadJet50", &passed_QuadJet50);
	mctrees[file.first]->SetBranchAddress("MR", &MR);
        mctrees[file.first]->SetBranchAddress("Rsq", &Rsq);
        mctrees[file.first]->SetBranchAddress("nVtx", &nVtx); // enable 
        mctrees[file.first]->SetBranchAddress("HT", &HT); // enable 
        mctrees[file.first]->SetBranchAddress("nSelectedJets", &nSelectedJets); // enable 
        mctrees[file.first]->SetBranchAddress("nBTaggedJets", &nBTaggedJets); // enable 
        mctrees[file.first]->SetBranchAddress("met", &met); // enable 
        mctrees[file.first]->SetBranchAddress("jet1pt", &jet1pt); // enable 
        mctrees[file.first]->SetBranchAddress("b1pt", &b1pt); // enable 
        mctrees[file.first]->SetBranchAddress("b2pt", &b2pt); // enable 
        mctrees[file.first]->SetBranchAddress("b1phi", &b1phi); // enable 
        mctrees[file.first]->SetBranchAddress("b2phi", &b2phi); // enable 
        mctrees[file.first]->SetBranchAddress("mbb", &mbb); // enable 
        mctrees[file.first]->SetBranchAddress("ptbb", &ptbb); // enable 
        mctrees[file.first]->SetBranchAddress("etabb", &etabb); // enable 
    }

    for(auto &file : datafiles){
        datatrees[file.first] = (TTree*)file.second->Get("HbbRazor");
        
        datatrees[file.first]->SetBranchStatus("*", 0); // disable all
        datatrees[file.first]->SetBranchStatus("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05", 1);
        datatrees[file.first]->SetBranchStatus("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05", 1);
        datatrees[file.first]->SetBranchStatus("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03", 1);
        datatrees[file.first]->SetBranchStatus("passed_DiJet80Eta2p6_BTagIP3DFastPVLoose", 1);
        datatrees[file.first]->SetBranchStatus("passed_QuadJet45", 1);
        datatrees[file.first]->SetBranchStatus("passed_QuadJet50", 1);
//new triggers
        datatrees[file.first]->SetBranchStatus("passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200", 1);
        datatrees[file.first]->SetBranchStatus("passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", 1);

        datatrees[file.first]->SetBranchStatus("nVtx", 1); // enable 
        datatrees[file.first]->SetBranchStatus("MR", 1);
        datatrees[file.first]->SetBranchStatus("Rsq", 1);
        datatrees[file.first]->SetBranchStatus("HT", 1);
        datatrees[file.first]->SetBranchStatus("nBTaggedJets", 1);
        datatrees[file.first]->SetBranchStatus("nSelectedJets", 1);
        datatrees[file.first]->SetBranchStatus("met", 1);
        datatrees[file.first]->SetBranchStatus("jet1pt", 1);
        datatrees[file.first]->SetBranchStatus("b1pt", 1);
        datatrees[file.first]->SetBranchStatus("b2pt", 1);
        datatrees[file.first]->SetBranchStatus("b1phi", 1);
        datatrees[file.first]->SetBranchStatus("b2phi", 1);
        datatrees[file.first]->SetBranchStatus("mbb", 1);
        datatrees[file.first]->SetBranchStatus("ptbb", 1);
        datatrees[file.first]->SetBranchStatus("etabb", 1);

        datatrees[file.first]->SetBranchAddress("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05", &passed_DiPFJet80_DiPFJet30_BTagCSVd07d05);
        datatrees[file.first]->SetBranchAddress("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05", &passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d05);
        datatrees[file.first]->SetBranchAddress("passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03", &passed_DiPFJet80_DiPFJet30_BTagCSVd07d05d03);
        datatrees[file.first]->SetBranchAddress("passed_DiJet80Eta2p6_BTagIP3DFastPVLoose", &passed_DiJet80Eta2p6_BTagIP3DFastPVLoose);
        datatrees[file.first]->SetBranchAddress("passed_QuadJet45", &passed_QuadJet45);
        datatrees[file.first]->SetBranchAddress("passed_QuadJet50", &passed_QuadJet50);
//new triggers
        datatrees[file.first]->SetBranchAddress("passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200", &passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200);
        datatrees[file.first]->SetBranchAddress("passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200);

	datatrees[file.first]->SetBranchAddress("MR", &MR);
        datatrees[file.first]->SetBranchAddress("Rsq", &Rsq);
        datatrees[file.first]->SetBranchAddress("nVtx", &nVtx); // enable 
        datatrees[file.first]->SetBranchAddress("HT", &HT); // enable 
        datatrees[file.first]->SetBranchAddress("nBTaggedJets", &nBTaggedJets); // enable 
        datatrees[file.first]->SetBranchAddress("nSelectedJets", &nSelectedJets); // enable 
        datatrees[file.first]->SetBranchAddress("met", &met); // enable 
        datatrees[file.first]->SetBranchAddress("jet1pt", &jet1pt); // enable 
        datatrees[file.first]->SetBranchAddress("b1pt", &b1pt); // enable 
        datatrees[file.first]->SetBranchAddress("b2pt", &b2pt); // enable 
        datatrees[file.first]->SetBranchAddress("b1phi", &b1phi); // enable 
        datatrees[file.first]->SetBranchAddress("b2phi", &b2phi); // enable 
        datatrees[file.first]->SetBranchAddress("mbb", &mbb); // enable 
        datatrees[file.first]->SetBranchAddress("ptbb", &ptbb); // enable 
        datatrees[file.first]->SetBranchAddress("etabb", &etabb); // enable 
    }
    //define cuts and histograms
    float nMRBins = 10;
    float nRsqBins = 8;
    float MRBinLowEdges[] = {300, 350, 400, 450, 550, 700, 900, 1200, 1600, 2500, 4000};
    float RsqBinLowEdges[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 0.80, 1.5};
    vector<string> cutSequence;
    vector<string> cutName;
    vector<string> cutSequenceData;
    vector<string> cutNameData;

//applying new data cuts
 //   cutSequence.push_back( "MR>350&&Rsq>0.035&&ptbb<110&&b1pt>60&&b2pt>40&&mbb>70&&mbb<190" ); 
   // cutName.push_back("MR>350&&Rsq>0.035&&ptbb<110&&b1pt>60&&b2pt>40&&mbb>70&&mbb<190" );

    //cutSequenceData.push_back( "MR>350&&Rsq>0.035&&ptbb<110&&b1pt>60&&b2pt>40&&mbb>70&&mbb<190" );
    //cutNameData.push_back("MR>350&&Rsq>0.035&&ptbb<110&&b1pt>60&&b2pt>40&&mbb>70&&mbb<190" );

    cutSequence.push_back( "MR>350&&Rsq>0.035&&ptbb<110&&b1pt>60&&b2pt>40&&mbb>70&&mbb<190" );
    cutName.push_back("MR>350&&Rsq>0.035&&ptbb<110&&b1pt>60&&b2pt>40&&mbb>70&&mbb<190" );

    cutSequenceData.push_back( "(passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200||passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200)&&MR>350&&Rsq>0.035&&ptbb<110&&b1pt>60&&b2pt>40&&mbb>70&&mbb<190" );
    cutNameData.push_back("(passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200||passed_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200)&&MR>350&&Rsq>0.035&&ptbb<110&&b1pt>60&&b2pt>40&&mbb>70&&mbb<190" );

    //cutSequence.push_back( "passed_QuadJet50&&MR>350&&Rsq>0.035&&ptbb<110&&(mbb<110||mbb>140)" );
    //cutName.push_back( "passed_QuadJet50&&MR>350&&Rsq>0.035&&ptbb<110&&(mbb<110||mbb>140)" );

    //cutSequence.push_back( "passed_QuadJet50&&MR>400&&Rsq>0.05&&ptbb<50&&(mbb<110||mbb>140)" );
    //cutName.push_back( "passed_QuadJet50&&MR>400&&Rsq>0.05&&ptbb<50&&(mbb<110||mbb>140)" );

    // cutSequence.push_back( "passed_DiPFJet80_DiPFJet30_BTagCSVd07d05&&MR>350&&Rsq>0.035&&ptbb<110&&(mbb<110||mbb>140)" );
    // cutName.push_back( "passed_DiPFJet80_DiPFJet30_BTagCSVd07d05&&MR>350&&Rsq>0.035&&ptbb<110&&(mbb<110||mbb>140)" );

    // cutSequence.push_back( "passed_DiPFJet80_DiPFJet30_BTagCSVd07d05&&MR>400&&Rsq>0.05&&ptbb<50&&(mbb<110||mbb>140)" );
    // cutName.push_back( "passed_DiPFJet80_DiPFJet30_BTagCSVd07d05&&MR>400&&Rsq>0.05&&ptbb<50&&(mbb<110||mbb>140)" );

    // cutSequence.push_back( "jet1pt>80&&b1pt>40&&b2pt>40&&ptbb<120&&HT>400&&MR>400&&Rsq>0.035&&(mbb<110||mbb>140)" );
    // cutName.push_back( "jet1pt>80&&b1pt>40&&b2pt>40&&ptbb<120&&HT>400&&MR>400&&Rsq>0.035&&(mbb<110||mbb>140)" );


    map<string, vector<TH1F *> > mcNJets, mcMR, mcRsq,  mcMet, mcNvtx, mcHT, mcj1pt, mcb1pt, mcb2pt, mcmbb, mcptbb, mcetabb, mcb1b2DPhi;
    vector<TH1F *> dataNJets, dataMR, dataRsq, dataMet, dataNvtx, dataHT, dataj1pt, datab1pt, datab2pt, datambb, dataptbb, dataetabb, datab1b2DPhi;
    for(auto &tree : mctrees){
        mcNJets[tree.first] = vector<TH1F *>();
        mcMR[tree.first] = vector<TH1F *>();
        mcRsq[tree.first] = vector<TH1F *>();
        mcMet[tree.first] = vector<TH1F *>();
        mcNvtx[tree.first] = vector<TH1F *>();
        mcHT[tree.first] = vector<TH1F *>();
        mcj1pt[tree.first] = vector<TH1F *>();
        mcb1pt[tree.first] = vector<TH1F *>();
        mcb2pt[tree.first] = vector<TH1F *>();
        mcmbb[tree.first] = vector<TH1F *>();
        mcptbb[tree.first] = vector<TH1F *>();
        mcetabb[tree.first] = vector<TH1F *>();
    }
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        for(auto &tree : mctrees){
            mcNJets[tree.first].push_back(new TH1F(Form("mcNJets%s%d", tree.first.c_str(), cut), Form("%s; Number of jets 40 GeV", cutName[cut].c_str()), 10, 0, 10));
            mcNvtx[tree.first].push_back(new TH1F(Form("mcNvtx%s%d", tree.first.c_str(), cut), Form("%s; NVtx (GeV)", cutName[cut].c_str()), 50, 0, 50));
            mcMR[tree.first].push_back(new TH1F(Form("mcMR%s%d", tree.first.c_str(), cut), Form("%s; MR (GeV)", cutName[cut].c_str()), 20, 300, 1000));//20, 300, 4000
            mcRsq[tree.first].push_back(new TH1F(Form("mcRsq%s%d", tree.first.c_str(), cut), Form("%s; Rsq (GeV)", cutName[cut].c_str()), 15, .035, .2));// nRsqBins, RsqBinLowEdges)
            mcMet[tree.first].push_back(new TH1F(Form("mcMet%s%d", tree.first.c_str(), cut), Form("%s; MET (GeV)", cutName[cut].c_str()), 200, 0, 1000));
            mcHT[tree.first].push_back(new TH1F(Form("mcHT%s%d", tree.first.c_str(), cut), Form("%s; HT (GeV)", cutName[cut].c_str()), 50, 0, 1000));
            mcj1pt[tree.first].push_back(new TH1F(Form("mcj1pt%s%d", tree.first.c_str(), cut), Form("%s; j1pt (GeV)", cutName[cut].c_str()), 100, 0, 500));
		//changed x axis dim
            mcb1pt[tree.first].push_back(new TH1F(Form("mcb1pt%s%d", tree.first.c_str(), cut), Form("%s; b1pt (GeV)", cutName[cut].c_str()), 30, 50, 200));//100,0 500

            mcb2pt[tree.first].push_back(new TH1F(Form("mcb2pt%s%d", tree.first.c_str(), cut), Form("%s; b2pt (GeV)", cutName[cut].c_str()), 30, 50, 200));//100,0 500
            mcmbb[tree.first].push_back(new TH1F(Form("mcmbb%s%d", tree.first.c_str(), cut), Form("%s; Mbb (GeV)", cutName[cut].c_str()), 12, 70, 190)); //40 0 400
            mcptbb[tree.first].push_back(new TH1F(Form("mcptbb%s%d", tree.first.c_str(), cut), Form("%s; Pt BB (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcetabb[tree.first].push_back(new TH1F(Form("mcetabb%s%d", tree.first.c_str(), cut), Form("%s; Eta BB (GeV)", cutName[cut].c_str()), 50, -3, 3));
            mcb1b2DPhi[tree.first].push_back(new TH1F(Form("mcb1b2DPhi%s%d", tree.first.c_str(), cut), Form("%s; DPhi BB", cutName[cut].c_str()), 50, 0, 3.2));

            mcNJets[tree.first][cut]->Sumw2();
            mcMR[tree.first][cut]->Sumw2();
            mcRsq[tree.first][cut]->Sumw2();
            mcMet[tree.first][cut]->Sumw2();
            mcNvtx[tree.first][cut]->Sumw2();
	    mcHT[tree.first][cut]->Sumw2();
	    mcj1pt[tree.first][cut]->Sumw2();
	    mcb1pt[tree.first][cut]->Sumw2();
	    mcb2pt[tree.first][cut]->Sumw2();
	    mcmbb[tree.first][cut]->Sumw2();
	    mcptbb[tree.first][cut]->Sumw2();	    
	    mcetabb[tree.first][cut]->Sumw2();
	    mcb1b2DPhi[tree.first][cut]->Sumw2();	    
        }
        dataNJets.push_back(new TH1F(Form("dataNJets%d", cut), Form("%s; Number of jets 40 GeV", cutNameData[cut].c_str()), 10, 0, 10));
        dataMR.push_back(new TH1F(Form("dataMR%d", cut), Form("%s; MR (GeV)", cutNameData[cut].c_str()), 20, 300, 1000));// ... 40000
        dataRsq.push_back(new TH1F(Form("dataRsq%d", cut), Form("%s; Rsq (GeV)", cutNameData[cut].c_str()), 15, .035, .2));// nRsqBins, RsqBinLowEdges));
        dataMet.push_back(new TH1F(Form("dataMet%d", cut), Form("%s; mcMet (GeV)", cutNameData[cut].c_str()), 200, 0, 1000));
        dataNvtx.push_back(new TH1F(Form("dataNvtx%d", cut), Form("%s; NVtx (GeV)", cutNameData[cut].c_str()), 50, 0, 50));
        dataHT.push_back(new TH1F(Form("dataHT%d", cut), Form("%s; HT (GeV)", cutNameData[cut].c_str()), 50, 0, 1000));
        dataj1pt.push_back(new TH1F(Form("dataj1pt%d", cut), Form("%s; j1pt (GeV)", cutNameData[cut].c_str()), 100, 0, 500));

        datab1pt.push_back(new TH1F(Form("datab1pt%d", cut), Form("%s; b1pt (GeV)", cutNameData[cut].c_str()), 16, 60, 140));//100,0 500

        datab2pt.push_back(new TH1F(Form("datab2pt%d", cut), Form("%s; b2pt (GeV)", cutNameData[cut].c_str()), 16, 40, 120));//100,0 500
        datambb.push_back(new TH1F(Form("datambb%d", cut), Form("%s; Mbb (GeV)", cutNameData[cut].c_str()), 12, 70, 190)); //40 0 400
        dataptbb.push_back(new TH1F(Form("dataptbb%d", cut), Form("%s; Pt BB (GeV)", cutNameData[cut].c_str()), 50, 0, 500));
        dataetabb.push_back(new TH1F(Form("dataetabb%d", cut), Form("%s; Eta BB (GeV)", cutNameData[cut].c_str()), 50, -3, 3));
        datab1b2DPhi.push_back(new TH1F(Form("datab1b2DPhi%d", cut), Form("%s; DPhi BB", cutNameData[cut].c_str()), 50, 0, 3.2));

        dataNJets[cut]->Sumw2();
        dataMR[cut]->Sumw2();
        dataRsq[cut]->Sumw2();
        dataMet[cut]->Sumw2();
	dataNvtx[cut]->Sumw2();
	dataHT[cut]->Sumw2();
	dataj1pt[cut]->Sumw2();
	datab1pt[cut]->Sumw2();
	datab2pt[cut]->Sumw2();
	datambb[cut]->Sumw2();
	dataptbb[cut]->Sumw2();
	dataetabb[cut]->Sumw2();
	datab1b2DPhi[cut]->Sumw2();
    }
    for(auto &tree : mctrees){
        cout << "Filling MC histograms: " << tree.first << endl;
        uint nEntries = tree.second->GetEntries();
        if(debug) nEntries = 100000;
        //make TTreeFormulas for selection cuts
        vector<TTreeFormula *> cuts;
        for(uint cut = 0; cut < cutSequence.size(); cut++){
            cuts.push_back(new TTreeFormula(Form("%sCutsFormula%d", tree.first.c_str(), cut), cutSequence[cut].c_str(), tree.second));
            cuts[cut]->GetNdata();
        }
        //loop over entries
        for(uint i = 0; i < nEntries; i++){
	  //get entry
	  tree.second->GetEntry(i); 
	  
	  //get event weight
	  float eventWeight = weight;
	  
	  //apply selection cuts and fill the appropriate histograms
	  for(uint cut = 0; cut < cutSequence.size(); cut++){
	    bool passesCut = cuts[cut]->EvalInstance();
	    if(!passesCut) continue;

	    (mcNJets[tree.first])[cut]->Fill(nSelectedJets, eventWeight);
	    (mcMR[tree.first])[cut]->Fill(MR, eventWeight);
	    (mcRsq[tree.first])[cut]->Fill(Rsq, eventWeight);
	    (mcMet[tree.first])[cut]->Fill(met, eventWeight);
	    (mcNvtx[tree.first])[cut]->Fill(nVtx, eventWeight);
	    (mcHT[tree.first])[cut]->Fill(HT, eventWeight);
	    (mcj1pt[tree.first])[cut]->Fill(jet1pt, eventWeight);
	    (mcb1pt[tree.first])[cut]->Fill(b1pt, eventWeight);
	    (mcb2pt[tree.first])[cut]->Fill(b2pt, eventWeight);
	    (mcmbb[tree.first])[cut]->Fill(mbb, eventWeight);
	    (mcptbb[tree.first])[cut]->Fill(ptbb, eventWeight);
	    (mcetabb[tree.first])[cut]->Fill(etabb, eventWeight);
	    (mcb1b2DPhi[tree.first])[cut]->Fill(acos(cos(b1phi-b2phi)), eventWeight);
	  }
        }

        for(uint cut = 0; cut < cutSequence.size(); cut++){
            delete cuts[cut];
        }
    }

    for(auto &tree : datatrees){
        cout << "Filling data histograms: " << tree.first << endl;
        uint nEntries = tree.second->GetEntries();
        if(debug) nEntries = 100000;
        //make TTreeFormulas for selection cuts
        vector<TTreeFormula *> cuts;
        for(uint cut = 0; cut < cutSequenceData.size(); cut++){
            cuts.push_back(new TTreeFormula(Form("%sCutsFormula%d", tree.first.c_str(), cut), cutSequenceData[cut].c_str(), tree.second));
            cuts[cut]->GetNdata();
        }
        //loop over entries
        for(uint i = 0; i < nEntries; i++){
            //get entry
            tree.second->GetEntry(i); 

            //get event weight
            float eventWeight = 1.0;
	    
            //apply selection cuts and fill the appropriate histograms
            for(uint cut = 0; cut < cutSequenceData.size(); cut++){
                bool passesCut = cuts[cut]->EvalInstance();
                if(!passesCut) continue;
		
		dataNJets[cut]->Fill(nSelectedJets, eventWeight);
                dataMR[cut]->Fill(MR, eventWeight);
                dataRsq[cut]->Fill(Rsq, eventWeight);
                dataMet[cut]->Fill(met, eventWeight);
		dataNvtx[cut]->Fill(nVtx, eventWeight);
 		dataHT[cut]->Fill(HT, eventWeight);
 		dataj1pt[cut]->Fill(jet1pt, eventWeight);
 		datab1pt[cut]->Fill(b1pt, eventWeight);
 		datab2pt[cut]->Fill(b2pt, eventWeight);
 		datambb[cut]->Fill(mbb, eventWeight);
 		dataptbb[cut]->Fill(ptbb, eventWeight);
 		dataetabb[cut]->Fill(etabb, eventWeight);
 		datab1b2DPhi[cut]->Fill(acos(cos(b1phi-b2phi)), eventWeight);
           }
        }
        for(uint cut = 0; cut < cutSequenceData.size(); cut++){
            delete cuts[cut];
        }
    }

    //print out plots
    TCanvas c("c", "c", 800, 600);
    //c.SetLogy();

    //colors and legend
    map<string, int> colors;
    colors["QCD"] = 8;
    colors["TTJets"] = 7;
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(dataNJets[0], "2016B Data, BTagCSV");
    legend->AddEntry(mcNJets["QCD"][0], "QCD MC");
    legend->AddEntry(mcNJets["TTJets"][0], "TTJets MC");
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        //create histogram stacks for MC
        THStack NumJetsMC(Form("NumJetsStack%d", cut), cutName[cut].c_str());
        THStack MRMC(Form("MRStack%d", cut), cutName[cut].c_str());
        THStack RsqMC(Form("RsqStack%d", cut), cutName[cut].c_str());
        THStack MetMC(Form("MetStack%d", cut), cutName[cut].c_str());
        THStack NVtxMC(Form("NVtxStack%d", cut), cutName[cut].c_str());
        THStack HTMC(Form("HTStack%d", cut), cutName[cut].c_str());
        THStack J1PTMC(Form("J1PTStack%d", cut), cutName[cut].c_str());
        THStack B1PTMC(Form("B1PTStack%d", cut), cutName[cut].c_str());
        THStack B2PTMC(Form("B2PTStack%d", cut), cutName[cut].c_str());
        THStack MBBMC(Form("MBBStack%d", cut), cutName[cut].c_str());
        THStack PTBBMC(Form("PTBBStack%d", cut), cutName[cut].c_str());
        THStack ETABBMC(Form("ETABBStack%d", cut), cutName[cut].c_str());
        THStack DPHIBBMC(Form("DPHIBBStack%d", cut), cutName[cut].c_str());
        
        //add the histograms to the stack in order
        vector<string> orderedtrees {"QCD", "TTJets"};
        for(auto &tree : orderedtrees){
	    mcNJets[tree][cut]->SetFillColor(colors[tree]);
            mcMR[tree][cut]->SetFillColor(colors[tree]);
            mcRsq[tree][cut]->SetFillColor(colors[tree]);
            mcMet[tree][cut]->SetFillColor(colors[tree]);
            mcNvtx[tree][cut]->SetFillColor(colors[tree]);
	    mcHT[tree][cut]->SetFillColor(colors[tree]);
	    mcj1pt[tree][cut]->SetFillColor(colors[tree]);
	    mcb1pt[tree][cut]->SetFillColor(colors[tree]);
	    mcb2pt[tree][cut]->SetFillColor(colors[tree]);
	    mcmbb[tree][cut]->SetFillColor(colors[tree]);
	    mcptbb[tree][cut]->SetFillColor(colors[tree]);
	    mcetabb[tree][cut]->SetFillColor(colors[tree]);
	    mcb1b2DPhi[tree][cut]->SetFillColor(colors[tree]);

            NumJetsMC.Add(mcNJets[tree][cut]);
            MRMC.Add(mcMR[tree][cut]);
            RsqMC.Add(mcRsq[tree][cut]);
	    MetMC.Add(mcMet[tree][cut]);
	    NVtxMC.Add(mcNvtx[tree][cut]);
	    HTMC.Add(mcHT[tree][cut]);
	    J1PTMC.Add(mcj1pt[tree][cut]);
	    B1PTMC.Add(mcb1pt[tree][cut]);
	    B2PTMC.Add(mcb2pt[tree][cut]);
	    MBBMC.Add(mcmbb[tree][cut]);
	    PTBBMC.Add(mcptbb[tree][cut]);
	    ETABBMC.Add(mcetabb[tree][cut]);
	    DPHIBBMC.Add(mcb1b2DPhi[tree][cut]);
        }
        DrawDataVsMCRatioPlot(dataNJets[cut], &NumJetsMC, legend, "Number of jets 40 GeV", "HbbRazor_newQCD_NumJets"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataMR[cut], &MRMC, legend, "MR (GeV)", "HbbRazor_newQCD_MR"+to_string(cut), true);
        DrawDataVsMCRatioPlot(dataRsq[cut], &RsqMC, legend, "Rsq", "HbbRazor_newQCD_Rsq"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataMet[cut], &MetMC, legend, "Met", "HbbRazor_newQCD_MET"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataNvtx[cut], &NVtxMC, legend, "NVtx", "HbbRazor_newQCD_NVtx"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataHT[cut], &HTMC, legend, "HT", "HbbRazor_newQCD_HT"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataj1pt[cut], &J1PTMC, legend, "J1PT", "HbbRazor_newQCD_J1PT"+to_string(cut), false);
        DrawDataVsMCRatioPlot(datab1pt[cut], &B1PTMC, legend, "B1PT", "HbbRazor_newQCD_B1PT"+to_string(cut), false);
        DrawDataVsMCRatioPlot(datab2pt[cut], &B2PTMC, legend, "B2PT", "HbbRazor_newQCD_B2PT"+to_string(cut), false);
        DrawDataVsMCRatioPlot(datambb[cut], &MBBMC, legend, "MBB", "HbbRazor_newQCD_MBB"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataptbb[cut], &PTBBMC, legend, "PTBB", "HbbRazor_newQCD_PTBB"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataetabb[cut], &ETABBMC, legend, "ETABB", "HbbRazor_newQCD_ETABB"+to_string(cut), false);
        DrawDataVsMCRatioPlot(datab1b2DPhi[cut], &DPHIBBMC, legend, "DPHIBB", "HbbRazor_newQCD_DPHIBB"+to_string(cut), false);
    }

    delete legend;
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        delete dataNJets[cut];
        delete dataMR[cut];
        delete dataRsq[cut];
 	delete dataMet[cut];
	delete dataNvtx[cut];
	delete dataHT[cut];
	delete dataj1pt[cut];
	delete datab1pt[cut];
	delete datab2pt[cut];
	delete datambb[cut];
	delete dataptbb[cut];
	delete dataetabb[cut];
	delete datab1b2DPhi[cut];
	for(auto &tree : mctrees){
	  delete mcNJets[tree.first][cut];
	  delete mcMR[tree.first][cut];
	  delete mcRsq[tree.first][cut];
	  delete mcMet[tree.first][cut];
	  delete mcNvtx[tree.first][cut];
	  delete mcHT[tree.first][cut];
	  delete mcj1pt[tree.first][cut];
	  delete mcb1pt[tree.first][cut];
	  delete mcb2pt[tree.first][cut];
	  delete mcmbb[tree.first][cut];
	  delete mcptbb[tree.first][cut];
	  delete mcetabb[tree.first][cut];
	  delete mcb1b2DPhi[tree.first][cut];
        }
    }

    mcfiles["QCD"]->Close();
    mcfiles["TTJets"]->Close();
    datafiles["Hbb"]->Close();
    delete mcfiles["QCD"];
    delete mcfiles["TTJets"];
    delete datafiles["Hbb"];
}

int main(){
    HbbRazorPlots();
    return 0;
}

void DrawDataVsMCRatioPlot(TH1F *dataHist, THStack *mcStack, TLegend *leg, string xaxisTitle, string printString, bool logX){
    TCanvas c("c", "c", 800, 600);
    c.Clear();
    c.cd();
    TPad pad1("pad1","pad1",0,0.4,1,1);
    pad1.SetBottomMargin(0);
    //pad1.SetLogy();
    if(logX) pad1.SetLogx();
    pad1.Draw();
    pad1.cd();
    dataHist->SetStats(0);
    dataHist->SetMinimum(0.5);
    mcStack->Draw("hist");
    dataHist->Draw("pesame");
    // dataHist->SetMinimum(0.5);
    // mcStack->Draw("hist");
    mcStack->GetYaxis()->SetTitle("Number of events in 6/fb");
    mcStack->GetYaxis()->SetLabelSize(0.03);
    mcStack->GetYaxis()->SetTitleOffset(0.45);
    mcStack->GetYaxis()->SetTitleSize(0.05);
    mcStack->SetMinimum(0.1);
    dataHist->SetMarkerStyle(20);
    dataHist->SetMarkerSize(1);
    dataHist->GetYaxis()->SetTitle("Number of events in 2.6/fb");
    dataHist->Draw("pesame");
    pad1.Modified();
    gPad->Update();
    //make ratio histogram
    TList * histList = (TList*)mcStack->GetHists();
    TIter next(histList);
    TH1 *mcTotal = (TH1*) histList->First()->Clone();
    TObject *obj;
    while((obj = next())){
        if(obj == histList->First()) continue;
        mcTotal->Add((TH1*)obj);
    }
    TH1F *dataOverMC = (TH1F*)dataHist->Clone();
    dataOverMC->SetTitle("");
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
    string histoName = dataHist->GetName() ;
    if(histoName.find("MR") != std::string::npos  )
      {
	cout<<"Number of events in data: "<<dataHist->Integral()<<" "<<printString<<endl;
	cout<<"Number of events in MC: "<<mcTotal->Integral()<<" "<<endl;
      }
    if(histoName.find("datadeltaPhi") != std::string::npos
       ||
       histoName.find("b1b2DPhi") != std::string::npos
       )
      {
	leg->SetX1NDC(0.1); leg->SetX2NDC(0.3); leg->SetY1NDC(0.7); leg->SetY2NDC(0.9);
      }
    else 
      {
	leg->SetX1NDC(0.7); leg->SetX2NDC(0.9); leg->SetY1NDC(0.7); leg->SetY2NDC(0.9);
      }
    if(histoName.find("mbb") != std::string::npos
       ||
       histoName.find("ptbb") != std::string::npos
       ||
       histoName.find("b1b2DPhi") != std::string::npos
       )
      {
	//pad1.SetLogy(0);
      }
    pad1.Modified();
    gPad->Update();
    leg->Draw();
    c.cd();
    TPad pad2("pad2","pad2",0,0.0,1,0.4);
    pad2.SetTopMargin(0);
    pad2.SetTopMargin(0.008);
    pad2.SetBottomMargin(0.25);
    pad2.SetGridy();
    if(logX) pad2.SetLogx();
    if(dataOverMC->GetMaximum() > 1.5) dataOverMC->SetMaximum(1.5);
    pad2.Draw();
    pad2.cd();
    dataOverMC->Draw("pe");
    pad2.Modified();
    gPad->Update();
    c.Print(Form("BINTEST1_%s.pdf", printString.c_str()));
    c.Print(Form("BINTEST1_%s.png", printString.c_str()));
    // c.Print(Form("%s.root", printString.c_str()));
}
