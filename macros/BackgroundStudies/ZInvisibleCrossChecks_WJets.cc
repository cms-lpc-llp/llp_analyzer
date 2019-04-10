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

void ZInvisibleCrossChecks(){
    gROOT->SetBatch();

    float maxMuonPt = 999;

    //set color palette 
    const Int_t NCont = 100;
    gStyle->SetNumberContours(NCont);
    gStyle->SetPaintTextFormat("1.0f");
    
    map<string, string> suffixes;
    suffixes["WJets"] = "_noW";
    suffixes["TTJets"] = "_noW";
    suffixes["Top"] = "_noW";
    suffixes["TTW"] = "_noW";
    suffixes["TTZ"] = "_noW";

    //get input files -- assumes one TFile for each process, with weights for different HT bins 
    map<string, TFile*> mcfiles;
    map<string, TFile*> datafiles;
    mcfiles["TTJets"] = new TFile("./TTJetsRun1_19700pb_weighted.root");
    mcfiles["Top"] = new TFile("./SingleTopRun1_19700pb_weighted.root");
    mcfiles["TTW"] = new TFile("./TTWJetsRun1_19700pb_weighted.root");
    mcfiles["TTZ"] = new TFile("./TTZJetsRun1_19700pb_weighted.root");
    mcfiles["WJets"] = new TFile("./WJetsRun1_19700pb_weighted.root");

    datafiles["WJets"] = new TFile("./SingleMuRun1_goodlumi.root");
    //get trees and set branches
    map<string, TTree*> mctrees;
    map<string, TTree*> datatrees;
    map<string, float> mets;
    map<string, float> mrs;
    map<string, float> rsqs;
    map<string, int> njets;
    map<string, int> njets80;
    map<string, float> hts;

    
    float weight;
    int nPU_mean, nVtx, nTightMuons, nLooseMuons, nBTaggedJets;
    float leadingMuonPt, leadingMuonEta, recoZpt, recoZeta, recoZmass, subleadingMuonPt, subleadingMuonEta, mTLepMet;
    bool Flag_HBHENoiseFilter, Flag_CSCTightHaloFilter, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_eeBadScFilter, Flag_ecalLaserCorrFilter;
    float bjet1Pt, bjet2Pt;
    bool bjet1PassMedium, bjet2PassMedium;
    float genZPt, genWPt;
    bool hlt_singlemu, hlt_dimuon;
    float leadingTightMuonPt, leadingTightMuonEta;
    float MR, Rsq;

    for(auto &file : mcfiles){
        mets[file.first] = 0.;
        mrs[file.first] = 0.;
        rsqs[file.first] = 0.;
        njets[file.first] = 0.;
        njets80[file.first] = 0.;
        hts[file.first] = 0.;

        mctrees[file.first] = (TTree*)file.second->Get("RazorInclusive");

        mctrees[file.first]->SetBranchStatus("*", 0); // disable all
        mctrees[file.first]->SetBranchStatus("weight", 1);
        mctrees[file.first]->SetBranchStatus("nPU_mean", 1);
        mctrees[file.first]->SetBranchStatus("nVtx", 1); // enable 
        mctrees[file.first]->SetBranchStatus("nBTaggedJets", 1);
        mctrees[file.first]->SetBranchStatus("leadingTightMuonEta", 1);
        mctrees[file.first]->SetBranchStatus("subleadingMuonPt", 1);
        mctrees[file.first]->SetBranchStatus("nTightMuons", 1);
        mctrees[file.first]->SetBranchStatus("nLooseMuons", 1);
        mctrees[file.first]->SetBranchStatus("genZpt", 1);
        mctrees[file.first]->SetBranchStatus("genWpt", 1);
        mctrees[file.first]->SetBranchStatus("bjet1PassMedium", 1);
        mctrees[file.first]->SetBranchStatus("bjet2PassMedium", 1);
        mctrees[file.first]->SetBranchStatus("bjet1Pt", 1);
        mctrees[file.first]->SetBranchStatus("bjet2Pt", 1);
        mctrees[file.first]->SetBranchStatus("hlt_singlemu", 1);
        mctrees[file.first]->SetBranchStatus("hlt_dimuon", 1);
        mctrees[file.first]->SetBranchStatus("mTLepMet", 1);
        mctrees[file.first]->SetBranchStatus("met_noW", 1);
        mctrees[file.first]->SetBranchStatus("MR_noW", 1);
        mctrees[file.first]->SetBranchStatus("Rsq_noW", 1);
        mctrees[file.first]->SetBranchStatus("HT_noW", 1);
        mctrees[file.first]->SetBranchStatus("numJets80_noW", 1);
        mctrees[file.first]->SetBranchStatus("numJets_noW", 1);
        mctrees[file.first]->SetBranchStatus("MR", 1);
        mctrees[file.first]->SetBranchStatus("Rsq", 1);

        mctrees[file.first]->SetBranchAddress(Form("met%s", suffixes[file.first].c_str()), &mets[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("MR%s", suffixes[file.first].c_str()), &mrs[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("Rsq%s", suffixes[file.first].c_str()), &rsqs[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("numJets%s", suffixes[file.first].c_str()), &njets[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("numJets80%s", suffixes[file.first].c_str()), &njets80[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("HT%s", suffixes[file.first].c_str()), &hts[file.first]);
        mctrees[file.first]->SetBranchAddress("MR", &MR);
        mctrees[file.first]->SetBranchAddress("Rsq", &Rsq);
        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("nPU_mean", &nPU_mean);
        mctrees[file.first]->SetBranchAddress("nVtx", &nVtx); // enable 
        mctrees[file.first]->SetBranchAddress("genZpt", &genZPt);
        mctrees[file.first]->SetBranchAddress("genWpt", &genWPt);
        mctrees[file.first]->SetBranchAddress("bjet1PassMedium", &bjet1PassMedium);
        mctrees[file.first]->SetBranchAddress("bjet2PassMedium", &bjet2PassMedium);
        mctrees[file.first]->SetBranchAddress("bjet1Pt", &bjet1Pt);
        mctrees[file.first]->SetBranchAddress("bjet2Pt", &bjet2Pt);
        mctrees[file.first]->SetBranchAddress("nBTaggedJets", &nBTaggedJets);
        mctrees[file.first]->SetBranchAddress("nTightMuons", &nTightMuons);
        mctrees[file.first]->SetBranchAddress("nLooseMuons", &nLooseMuons);
        mctrees[file.first]->SetBranchAddress("leadingMuonPt", &leadingMuonPt);
        mctrees[file.first]->SetBranchAddress("leadingMuonEta", &leadingMuonEta);
        mctrees[file.first]->SetBranchAddress("hlt_singlemu", &hlt_singlemu);
        mctrees[file.first]->SetBranchAddress("hlt_dimuon", &hlt_dimuon);
        mctrees[file.first]->SetBranchAddress("mTLepMet", &mTLepMet);
    }
    for(auto &file : datafiles){
        datatrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        
        datatrees[file.first]->SetBranchStatus("*", 0); // disable all
        datatrees[file.first]->SetBranchStatus("nVtx", 1); // enable 
        datatrees[file.first]->SetBranchStatus("Flag_HBHENoiseFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_CSCTightHaloFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_eeBadScFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_ecalLaserCorrFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("nBTaggedJets", 1);
        datatrees[file.first]->SetBranchStatus("nTightMuons", 1);
        datatrees[file.first]->SetBranchStatus("nLooseMuons", 1);
        datatrees[file.first]->SetBranchStatus("hlt_dimuon", 1);
        datatrees[file.first]->SetBranchStatus("hlt_singlemu", 1);
        datatrees[file.first]->SetBranchStatus("leadingTightMuonPt", 1);
        datatrees[file.first]->SetBranchStatus("leadingTightMuonEta", 1);
        datatrees[file.first]->SetBranchStatus("hlt_razor", 1);
        datatrees[file.first]->SetBranchStatus("bjet1PassMedium", 1);
        datatrees[file.first]->SetBranchStatus("bjet2PassMedium", 1);
        datatrees[file.first]->SetBranchStatus("bjet1Pt", 1);
        datatrees[file.first]->SetBranchStatus("bjet2Pt", 1);
        datatrees[file.first]->SetBranchStatus("hlt_singlemu", 1);
        datatrees[file.first]->SetBranchStatus("hlt_dimuon", 1);
        datatrees[file.first]->SetBranchStatus("mTLepMet", 1);
        datatrees[file.first]->SetBranchStatus("met_noW", 1);
        datatrees[file.first]->SetBranchStatus("MR_noW", 1);
        datatrees[file.first]->SetBranchStatus("Rsq_noW", 1);
        datatrees[file.first]->SetBranchStatus("HT_noW", 1);
        datatrees[file.first]->SetBranchStatus("numJets80_noW", 1);
        datatrees[file.first]->SetBranchStatus("numJets_noW", 1);

        datatrees[file.first]->SetBranchAddress(Form("met%s", suffixes[file.first].c_str()), &mets[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("MR%s", suffixes[file.first].c_str()), &mrs[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("Rsq%s", suffixes[file.first].c_str()), &rsqs[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("numJets%s", suffixes[file.first].c_str()), &njets[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("numJets80%s", suffixes[file.first].c_str()), &njets80[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("HT%s", suffixes[file.first].c_str()), &hts[file.first]);
        datatrees[file.first]->SetBranchAddress("nVtx", &nVtx); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter); // enable 
        datatrees[file.first]->SetBranchAddress("mTLepMet", &mTLepMet);
        datatrees[file.first]->SetBranchAddress("bjet1PassMedium", &bjet1PassMedium);
        datatrees[file.first]->SetBranchAddress("bjet2PassMedium", &bjet2PassMedium);
        datatrees[file.first]->SetBranchAddress("bjet1Pt", &bjet1Pt);
        datatrees[file.first]->SetBranchAddress("bjet2Pt", &bjet2Pt);
        datatrees[file.first]->SetBranchAddress("nBTaggedJets", &nBTaggedJets);
        datatrees[file.first]->SetBranchAddress("nTightMuons", &nTightMuons);
        datatrees[file.first]->SetBranchAddress("hlt_singlemu", &hlt_singlemu);
        datatrees[file.first]->SetBranchAddress("leadingTightMuonPt", &leadingTightMuonPt);
        datatrees[file.first]->SetBranchAddress("leadingTightMuonEta", &leadingTightMuonEta);
        datatrees[file.first]->SetBranchAddress("hlt_dimuon", &hlt_dimuon);
    }

    //load efficiency/acceptance histograms
    TFile effFile("Run1LeptonPhotonEfficiency.root");
    TH2F muonLooseEffHisto = *(TH2F *)effFile.Get("MuonEfficiency");
    TH2F muonTightEffHisto = *(TH2F *)effFile.Get("MuonEfficiencyTight");
    TH2F zAccHisto = *(TH2F *)effFile.Get("MuonAcceptance");

    //load muon efficiency scale factor histogram
    TFile muIdSFFile("data/ScaleFactors/MuonEfficiencies_ID_Run2012ReReco_53X.root");
    TFile muIsoSFFile("data/ScaleFactors/MuonEfficiencies_ISO_Run_2012ReReco_53X.root");
    //TODO: add muon efficiency scale factors
    float singleMuTriggerSF = 0.97;
    float doubleMuTriggerSF = 0.97;
    float doubleMuNormalizationSF = 0.97;

    //load pileup reweighting histogram
    TFile *pileupWeightFile = new TFile("data/Run1PileupWeights.root", "READ");
    TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("PUWeight_Run1");
    assert(pileupWeightHist);

    //load TTbar scale factor histograms
    TFile *TTBarDileptonScaleFactorsFile = new TFile("data/ScaleFactors/Run1/TTBarDileptonScaleFactors.root");
    TH2F *TTBarDileptonScaleFactor = (TH2F*)TTBarDileptonScaleFactorsFile->Get("TTBarDileptonScaleFactor");
    TFile *TTBarSingleLeptonScaleFactorsFile = new TFile("data/ScaleFactors/Run1/TTBarSingleLeptonScaleFactors.root");
    TH2F *TTBarSingleLeptonScaleFactor = (TH2F*)TTBarSingleLeptonScaleFactorsFile->Get("TTBarSingleLeptonScaleFactor");
    float maxMRForTTJetsSF = 650;
    float maxRsqForTTJetsSF = 1.0;
    
    //define cuts and histograms
    float nMRBins = 10;
    float nRsqBins = 8;
    float MRBinLowEdges[] = {300, 350, 400, 450, 550, 700, 900, 1200, 1600, 2500, 4000};
    float RsqBinLowEdges[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 0.80, 1.5};
    vector<string> cutSequence;
    vector<string> cutName;

    cutSequence.push_back( "nBTaggedJets == 0 && nTightMuons == 1 && nLooseMuons == 1 && hlt_singlemu && MR_noW > 300 && Rsq_noW > 0.15 && numJets80_noW > 1 && mTLepMet > 30 && mTLepMet < 100" );
    cutName.push_back( "W+Jets Control Sample" );

    map<string, vector<TH1F *> > mcNJets, mcNJets80, mcMR, mcRsq,  mcMet, mcNvtx,  mcHT;
    vector<TH1F *>  dataNJets, dataNJets80, dataMR, dataRsq, dataMet, dataNvtx, dataHT;
    for(auto &tree : mctrees){
        mcNJets[tree.first] = vector<TH1F *>();
        mcNJets80[tree.first] = vector<TH1F *>();
        mcMR[tree.first] = vector<TH1F *>();
        mcRsq[tree.first] = vector<TH1F *>();
        mcMet[tree.first] = vector<TH1F *>();
        mcNvtx[tree.first] = vector<TH1F *>();
        mcHT[tree.first] = vector<TH1F *>();
    }
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        for(auto &tree : mctrees){
            mcNJets[tree.first].push_back(new TH1F(Form("mcNJets%s%d", tree.first.c_str(), cut), Form("%s; Number of jets 40 GeV", cutName[cut].c_str()), 10, 0, 10));
            mcNvtx[tree.first].push_back(new TH1F(Form("mcNvtx%s%d", tree.first.c_str(), cut), Form("%s; NVtx (GeV)", cutName[cut].c_str()), 50, 0, 50));
            mcNJets80[tree.first].push_back(new TH1F(Form("mcNJets80%s%d", tree.first.c_str(), cut), Form("%s; Number of jets 80 GeV", cutName[cut].c_str()), 10, 0, 10));
            mcMR[tree.first].push_back(new TH1F(Form("mcMR%s%d", tree.first.c_str(), cut), Form("%s; MR (GeV)", cutName[cut].c_str()), 20, 300, 4000));
            // mcMR[tree.first].push_back(new TH1F(Form("mcMR%s%d", tree.first.c_str(), cut), Form("%s; MR (GeV)", cutName[cut].c_str()), nMRBins, MRBinLowEdges));
            mcRsq[tree.first].push_back(new TH1F(Form("mcRsq%s%d", tree.first.c_str(), cut), Form("%s; Rsq (GeV)", cutName[cut].c_str()), nRsqBins, RsqBinLowEdges));
            mcMet[tree.first].push_back(new TH1F(Form("mcMet%s%d", tree.first.c_str(), cut), Form("%s; MET (GeV)", cutName[cut].c_str()), 200, 0, 1000));
            mcHT[tree.first].push_back(new TH1F(Form("mcHT%s%d", tree.first.c_str(), cut), Form("%s; HT (GeV)", cutName[cut].c_str()), 50, 0, 1000));

            mcNJets[tree.first][cut]->Sumw2();
            mcNJets80[tree.first][cut]->Sumw2();
            mcMR[tree.first][cut]->Sumw2();
            mcRsq[tree.first][cut]->Sumw2();
            mcMet[tree.first][cut]->Sumw2();
            mcNvtx[tree.first][cut]->Sumw2();
	    mcHT[tree.first][cut]->Sumw2();	    
        }
        dataNJets.push_back(new TH1F(Form("dataNJets%d", cut), Form("%s; Number of jets 40 GeV", cutName[cut].c_str()), 10, 0, 10));
        dataNJets80.push_back(new TH1F(Form("dataNJets80%d", cut), Form("%s; Number of jets 80 GeV", cutName[cut].c_str()), 10, 0, 10));
        // dataMR.push_back(new TH1F(Form("dataMR%d", cut), Form("%s; MR (GeV)", cutName[cut].c_str()), nMRBins, MRBinLowEdges));
        dataMR.push_back(new TH1F(Form("dataMR%d", cut), Form("%s; MR (GeV)", cutName[cut].c_str()), 20, 300, 4000));
        dataRsq.push_back(new TH1F(Form("dataRsq%d", cut), Form("%s; Rsq (GeV)", cutName[cut].c_str()), nRsqBins, RsqBinLowEdges));
        dataMet.push_back(new TH1F(Form("dataMet%d", cut), Form("%s; mcMet (GeV)", cutName[cut].c_str()), 200, 0, 1000));
        dataNvtx.push_back(new TH1F(Form("dataNvtx%d", cut), Form("%s; NVtx (GeV)", cutName[cut].c_str()), 50, 0, 50));
        dataHT.push_back(new TH1F(Form("dataHT%d", cut), Form("%s; HT (GeV)", cutName[cut].c_str()), 50, 0, 1000));

        dataNJets[cut]->Sumw2();
        dataNJets80[cut]->Sumw2();
        dataMR[cut]->Sumw2();
        dataRsq[cut]->Sumw2();
        dataMet[cut]->Sumw2();
	dataNvtx[cut]->Sumw2();
	dataHT[cut]->Sumw2();
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
            eventWeight *= pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(nPU_mean));

            //reweigh according to selection efficiency and acceptance
	    if(!computeDataOverMCSFs){
	      double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindFixBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));
		    
	      if(effFactor > 1e-5) eventWeight /= effFactor;
	      else{ 
		eventWeight = 0;
	      }
	    }

	    //trigger and normalization scale factors
	    eventWeight *= singleMuTriggerSF;
	    
	    //b-tagging scale factors
	    double btagScaleFactor = 1.0;
	    double bjet1EventScaleFactor = 1.0;
	    double bjet2EventScaleFactor = 1.0;
	    if (bjet1Pt > 20) {
	      double bjet1ScaleFactor = 0.938887 + 0.00017124 * bjet1Pt + (-2.76366e-07) * bjet1Pt * bjet1Pt ;
	      double MCEff = 1.0;
	      if (bjet1Pt < 50) MCEff = 0.65;
	      else if (bjet1Pt < 80) MCEff = 0.70;
	      else if (bjet1Pt < 120) MCEff = 0.73;
	      else if (bjet1Pt < 210) MCEff = 0.73;
	      else MCEff = 0.66;				 
	      if (bjet1PassMedium) bjet1EventScaleFactor = bjet1ScaleFactor;
	      else bjet1EventScaleFactor = ( 1/MCEff - bjet1ScaleFactor) / ( 1/MCEff - 1);
	    }
	    if (bjet2Pt > 20) {
	      double bjet2ScaleFactor = 0.938887 + 0.00017124 * bjet2Pt + (-2.76366e-07) * bjet2Pt * bjet2Pt ;
	      double MCEff = 1.0;
	      if (bjet2Pt < 50) MCEff = 0.65;
	      else if (bjet2Pt < 80) MCEff = 0.70;
	      else if (bjet2Pt < 120) MCEff = 0.73;
	      else if (bjet2Pt < 210) MCEff = 0.73;
	      else MCEff = 0.66;		 
	      if (bjet2PassMedium) bjet2EventScaleFactor = bjet2ScaleFactor;
	      else bjet2EventScaleFactor = ( 1/MCEff - bjet2ScaleFactor) / ( 1/MCEff - 1);
	    }
	    btagScaleFactor = bjet1EventScaleFactor * bjet1EventScaleFactor;
	    
	    eventWeight *= btagScaleFactor;
	    //TTJets SF
	    if(tree.first == "TTJets"){
	      double ttjetsSF = TTBarDileptonScaleFactor->GetBinContent(TTBarDileptonScaleFactor->FindFixBin(min(MR, maxMRForTTJetsSF), min(Rsq, maxRsqForTTJetsSF)));
	      if(ttjetsSF > 1e-5) eventWeight *= ttjetsSF;
	    }

            //apply selection cuts and fill the appropriate histograms
            for(uint cut = 0; cut < cutSequence.size(); cut++){
                bool passesCut = cuts[cut]->EvalInstance();
                if(!passesCut) continue;


                (mcNJets[tree.first])[cut]->Fill(njets[tree.first], eventWeight);
                (mcNJets80[tree.first])[cut]->Fill(njets80[tree.first], eventWeight);
                (mcMR[tree.first])[cut]->Fill(mrs[tree.first], eventWeight);
                (mcRsq[tree.first])[cut]->Fill(rsqs[tree.first], eventWeight);
                (mcMet[tree.first])[cut]->Fill(mets[tree.first], eventWeight);
                (mcNvtx[tree.first])[cut]->Fill(nVtx, eventWeight);
                (mcHT[tree.first])[cut]->Fill(hts[tree.first], eventWeight);
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
        for(uint cut = 0; cut < cutSequence.size(); cut++){
            cuts.push_back(new TTreeFormula(Form("%sCutsFormula%d", tree.first.c_str(), cut), cutSequence[cut].c_str(), tree.second));
            cuts[cut]->GetNdata();
        }
        //loop over entries
        for(uint i = 0; i < nEntries; i++){
            //get entry
            tree.second->GetEntry(i); 

            //get event weight
            float eventWeight = 1.0;

	    if(!Flag_HBHENoiseFilter || !Flag_CSCTightHaloFilter || !Flag_eeBadScFilter ) continue;
	    
            //reweigh according to selection efficiency and acceptance
	    if(!computeDataOverMCSFs){
	      double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindFixBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));

	      if(effFactor > 1e-5) eventWeight /= effFactor;
	      else{ 
		eventWeight = 0;
	      }
	    }

            //apply selection cuts and fill the appropriate histograms
            for(uint cut = 0; cut < cutSequence.size(); cut++){
                bool passesCut = cuts[cut]->EvalInstance();
                if(!passesCut) continue;

		dataNJets[cut]->Fill(njets[tree.first], eventWeight);
                dataNJets80[cut]->Fill(njets80[tree.first], eventWeight);
                dataMR[cut]->Fill(mrs[tree.first], eventWeight);
                dataRsq[cut]->Fill(rsqs[tree.first], eventWeight);
                dataMet[cut]->Fill(mets[tree.first], eventWeight);
		dataNvtx[cut]->Fill(nVtx, eventWeight);
 		dataHT[cut]->Fill(hts[tree.first], eventWeight);
           }
        }
        for(uint cut = 0; cut < cutSequence.size(); cut++){
            delete cuts[cut];
        }
    }

    //print out plots
    TCanvas c("c", "c", 800, 600);
    c.SetLogy();

    //colors and legend
    map<string, int> colors;
    colors["WJets"] = kOrange+10;
    colors["Top"] = kViolet-5;
    colors["TTJets"] = kViolet-6;
    colors["TTW"] = kRed+2;
    colors["TTZ"] = kOrange-3;
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(dataNJets[0], "2012 Data, WJets CS");
    legend->AddEntry(mcNJets["WJets"][0], "W+Jets MC");
    legend->AddEntry(mcNJets["TTJets"][0], "TTJets MC MC");
    legend->AddEntry(mcNJets["Top"][0], "Single Top MC");
    legend->AddEntry(mcNJets["TTW"][0], "TTW MC");
    legend->AddEntry(mcNJets["TTZ"][0], "TTZ MC");
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        //create histogram stacks for MC
        THStack NumJetsMC(Form("NumJetsStack%d", cut), cutName[cut].c_str());
        THStack NumJets80MC(Form("NumJets80Stack%d", cut), cutName[cut].c_str());
        THStack MRMC(Form("MRStack%d", cut), cutName[cut].c_str());
        THStack RsqMC(Form("RsqStack%d", cut), cutName[cut].c_str());
        THStack MetMC(Form("MetStack%d", cut), cutName[cut].c_str());
        THStack NVtxMC(Form("NVtxStack%d", cut), cutName[cut].c_str());
        THStack HTMC(Form("HTStack%d", cut), cutName[cut].c_str());
        THStack DPhiMC(Form("DPhiMC%d", cut), cutName[cut].c_str());

        //add the histograms to the stack in order
        vector<string> orderedtrees {"TTZ", "TTW", "Top", "TTJets", "WJets"};
        for(auto &tree : orderedtrees){
	    mcNJets[tree][cut]->SetFillColor(colors[tree]);
            mcNJets80[tree][cut]->SetFillColor(colors[tree]);
            mcMR[tree][cut]->SetFillColor(colors[tree]);
            mcRsq[tree][cut]->SetFillColor(colors[tree]);
            mcMet[tree][cut]->SetFillColor(colors[tree]);
            mcNvtx[tree][cut]->SetFillColor(colors[tree]);
	    mcHT[tree][cut]->SetFillColor(colors[tree]);

            NumJetsMC.Add(mcNJets[tree][cut]);
            NumJets80MC.Add(mcNJets80[tree][cut]);
            MRMC.Add(mcMR[tree][cut]);
            RsqMC.Add(mcRsq[tree][cut]);
	    MetMC.Add(mcMet[tree][cut]);
	    NVtxMC.Add(mcNvtx[tree][cut]);
	    HTMC.Add(mcHT[tree][cut]);
        }
        DrawDataVsMCRatioPlot(dataNJets[cut], &NumJetsMC, legend, "Number of jets 40 GeV", "WJetsCrossChecks_NumJets"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataNJets80[cut], &NumJets80MC, legend, "Number of jets 80 GeV", "WJetsCrossChecks_NumJets80"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataMR[cut], &MRMC, legend, "MR (GeV)", "WJetsCrossChecks_MR"+to_string(cut), true);
        DrawDataVsMCRatioPlot(dataRsq[cut], &RsqMC, legend, "Rsq", "WJetsCrossChecks_Rsq"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataMet[cut], &MetMC, legend, "Met", "WJetsCrossChecks_MET"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataNvtx[cut], &NVtxMC, legend, "NVtx", "WJetsCrossChecks_NVtx"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataHT[cut], &HTMC, legend, "HT", "WJetsCrossChecks_HT"+to_string(cut), false);
    }

    delete legend;
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        delete dataNJets[cut];
        delete dataNJets80[cut];
        delete dataMR[cut];
        delete dataRsq[cut];
 	delete dataMet[cut];
	delete dataNvtx[cut];
	delete dataHT[cut];
	for(auto &tree : mctrees){
	  delete mcNJets[tree.first][cut];
	  delete mcNJets80[tree.first][cut];
	  delete mcMR[tree.first][cut];
	  delete mcRsq[tree.first][cut];
	  delete mcMet[tree.first][cut];
	  delete mcNvtx[tree.first][cut];
        }
    }
    
    mcfiles["WJets"]->Close();
    mcfiles["Top"]->Close();
    mcfiles["TTJets"]->Close();
    mcfiles["TTW"]->Close();
    mcfiles["TTZ"]->Close();
    datafiles["WJets"]->Close();
    pileupWeightFile->Close();
    delete mcfiles["WJets"];
    delete mcfiles["Top"];
    delete mcfiles["TTJets"];
    delete mcfiles["TTW"];
    delete mcfiles["TTZ"];
    delete datafiles["WJets"];
    delete pileupWeightFile;
}

int main(){
    ZInvisibleCrossChecks();
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
    if(histoName.find("datadeltaPhi") != std::string::npos  )
      {
	leg->SetX1NDC(0.1); leg->SetX2NDC(0.3); leg->SetY1NDC(0.7); leg->SetY2NDC(0.9);
      }
    else 
      {
	leg->SetX1NDC(0.7); leg->SetX2NDC(0.9); leg->SetY1NDC(0.7); leg->SetY2NDC(0.9);
      }
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
    c.Print(Form("%s.pdf", printString.c_str()));
    // c.Print(Form("%s.root", printString.c_str()));
}
