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
#include "assert.h"
#include "math.h"

using namespace std;

//true: find the translation factors from MC to data
//false: find the translation factors from DY, W, G to Z->nunu
bool computeDataOverMCSFs = true;

void DrawDataVsMCRatioPlot(TH1F *dataHist, THStack *mcStack, TLegend *leg, string xaxisTitle, string printString, bool logX);

void ZInvisibleControlSamples(){
    gROOT->SetBatch();

    //set color palette 
    const Int_t NCont = 100;
    gStyle->SetNumberContours(NCont);
    gStyle->SetPaintTextFormat("1.0f");

    //choose which sample to normalize to
    //string normalizeTo = "DYJets";
    string normalizeTo = "ZJets";

    //for plots
    float MetMin = 0.;
    float MetMax = 1000;
    float nMetBins = 20;
    float nMRBins = 8;
    float nRsqBins = 7;
    float MRBinLowEdges[] = {300, 350, 400, 450, 550, 700, 900, 1200, 4000};
    float RsqBinLowEdges[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 1.5};

    //upper bounds of reweighting histograms
    float maxPhotonPt = 999; 
    float maxMuonPt = 999;
    float maxZPt = 2999;

    //decide to reweight by MET or by MR and Rsq
    //(reweightByRazor = false to reweight by MET, true to reweight by MR and Rsq)
    bool reweightByRazor = true; 

    map<string, string> suffixes;
    suffixes["DYJets"] = "_noZ";
    suffixes["TTJetsDY"] = "_noZ";
    suffixes["TopDY"] = "_noZ";
    suffixes["TTWDY"] = "_noZ";
    suffixes["TTZDY"] = "_noZ";
    suffixes["WJets"] = "_noW";
    suffixes["TTJets"] = "_noW";
    suffixes["Top"] = "_noW";
    suffixes["TTW"] = "_noW";
    suffixes["TTZ"] = "_noW";
    suffixes["GJets"] = "_noPho";
    suffixes["QCD"] = "_noPho";
    suffixes["WG"] = "_noPho";
    suffixes["ZG"] = "_noPho";
    suffixes["TTG"] = "_noPho";
    suffixes["ZJets"] = "";

    map<string, string> cuts;
    cuts["DYJets"] = "nLooseMuons == 2 && hlt_dimuon && recoZmass > 71 && recoZmass < 111 && MR_noZ > 300 && Rsq_noZ > 0.15 && numJets80_noZ > 1";
    cuts["WJets"] = "nBTaggedJets == 0 && nTightMuons == 1 && nLooseMuons == 1 && hlt_singlemu && MR_noW > 300 && Rsq_noW > 0.15 && numJets80_noW > 1 && mTLepMet > 30 && mTLepMet < 100";
    cuts["GJets"] = "hlt_photon && MR_noPho > 300 && Rsq_noPho > 0.15 && numJets80_noPho > 1 && leadingPhotonPt > 80";
    cuts["ZJets"] = "nLooseMuons == 0 && nLooseElectrons == 0 && hlt_razor && MR > 300 && Rsq > 0.15 && numJets80 > 1";

    cuts["TTJets"] = cuts["WJets"];
    cuts["Top"] = cuts["WJets"];
    cuts["TTZ"] = cuts["WJets"];
    cuts["TTW"] = cuts["WJets"];
    cuts["TTJetsDY"] = cuts["DYJets"];
    cuts["TopDY"] = cuts["DYJets"];
    cuts["TTWDY"] = cuts["DYJets"];
    cuts["TTZDY"] = cuts["DYJets"];
    cuts["QCD"] = cuts["GJets"];
    cuts["WG"] = cuts["GJets"];
    cuts["ZG"] = cuts["GJets"];
    cuts["TTG"] = cuts["GJets"];

    //get input files -- assumes one TFile for each process, with weights for different HT bins 
    map<string, TFile*> mcfiles;
    map<string, TFile*> datafiles;
    //signal processes
    mcfiles["DYJets"] = new TFile("./DYJetsRun1_19700pb_weighted.root");
    mcfiles["WJets"] = new TFile("./WJetsRun1_19700pb_weighted.root");
    mcfiles["GJets"] = new TFile("./GJetsRun1_19700pb_weighted.root");
    mcfiles["ZJets"] = new TFile("./ZJetsRun1_19700pb_weighted.root");
    //backgrounds for WJets
    mcfiles["TTJets"] = new TFile("./TTJetsRun1_19700pb_weighted.root");
    mcfiles["Top"] = new TFile("./SingleTopRun1_19700pb_weighted.root");
    mcfiles["TTW"] = new TFile("./TTWJetsRun1_19700pb_weighted.root");
    mcfiles["TTZ"] = new TFile("./TTZJetsRun1_19700pb_weighted.root");
    //backgrounds for DYJets
    mcfiles["TTJetsDY"] = mcfiles["TTJets"];
    mcfiles["TopDY"] = mcfiles["Top"];
    mcfiles["TTWDY"] = mcfiles["TTW"];
    mcfiles["TTZDY"] = mcfiles["TTZ"];
    //backgrounds for GJets
    mcfiles["QCD"] = new TFile("./QCDRun1_19700pb_weighted.root");
    mcfiles["WG"] = new TFile("./WGJetsRun1_19700pb_weighted.root");
    mcfiles["ZG"] = new TFile("./ZGJetsRun1_19700pb_weighted.root");
    mcfiles["TTG"] = new TFile("./TTGJetsRun1_19700pb_weighted.root");
    //data
    datafiles["DYJets"] = new TFile("./DoubleMuRun1_goodlumi.root");
    datafiles["WJets"] = new TFile("./SingleMuRun1_goodlumi.root");
    datafiles["GJets"] = new TFile("./PhotonRun1_goodlumi.root");

    //get trees and set branches
    map<string, TTree*> mctrees;
    map<string, TTree*> datatrees;
    map<string, float> mets;
    map<string, float> mrs;
    map<string, float> rsqs;
    float weight;
    float leadingMuonPt, leadingMuonEta, leadingPhotonPt, leadingPhotonEta, recoZpt, recoZeta, recoZmass, subleadingMuonPt, subleadingMuonEta, mTLepMet;
    float leadingTightMuonPt, leadingTightMuonEta;
    float bjet1Pt, bjet2Pt;
    bool bjet1PassMedium, bjet2PassMedium;
    bool passedHLTPhoton50, passedHLTPhoton75, passedHLTPhoton90, passedHLTPhoton135, passedHLTPhoton150;
    bool Flag_HBHENoiseFilter, Flag_CSCTightHaloFilter, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_eeBadScFilter, Flag_ecalLaserCorrFilter;
    float genZPt, genWPt;
    int nPU_mean;
    float MR, Rsq;
    for(auto &file : mcfiles){
        mets[file.first] = 0.;
        mrs[file.first] = 0.;
        rsqs[file.first] = 0.;
        mctrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        mctrees[file.first]->SetBranchStatus("*", 0);
        mctrees[file.first]->SetBranchStatus("met_noW", 1);
        mctrees[file.first]->SetBranchStatus("MR_noW", 1);
        mctrees[file.first]->SetBranchStatus("Rsq_noW", 1);
        mctrees[file.first]->SetBranchStatus("numJets80_noW", 1);
        mctrees[file.first]->SetBranchStatus("met_noZ", 1);
        mctrees[file.first]->SetBranchStatus("MR_noZ", 1);
        mctrees[file.first]->SetBranchStatus("Rsq_noZ", 1);
        mctrees[file.first]->SetBranchStatus("numJets80_noZ", 1);
        mctrees[file.first]->SetBranchStatus("met_noPho", 1);
        mctrees[file.first]->SetBranchStatus("MR_noPho", 1);
        mctrees[file.first]->SetBranchStatus("Rsq_noPho", 1);
        mctrees[file.first]->SetBranchStatus("numJets80_noPho", 1);
        mctrees[file.first]->SetBranchStatus("met", 1);
        mctrees[file.first]->SetBranchStatus("MR", 1);
        mctrees[file.first]->SetBranchStatus("Rsq", 1);
        mctrees[file.first]->SetBranchStatus("numJets80", 1);
        mctrees[file.first]->SetBranchStatus("weight", 1);
        mctrees[file.first]->SetBranchStatus("leadingMuonPt", 1);
        mctrees[file.first]->SetBranchStatus("leadingMuonEta", 1);
        mctrees[file.first]->SetBranchStatus("leadingTightMuonPt", 1);
        mctrees[file.first]->SetBranchStatus("leadingTightMuonEta", 1);
        mctrees[file.first]->SetBranchStatus("subleadingMuonPt", 1);
        mctrees[file.first]->SetBranchStatus("subleadingMuonEta", 1);
        mctrees[file.first]->SetBranchStatus("leadingPhotonPt", 1);
        mctrees[file.first]->SetBranchStatus("leadingPhotonEta", 1);
        mctrees[file.first]->SetBranchStatus("recoZpt", 1);
        mctrees[file.first]->SetBranchStatus("recoZeta", 1);
        mctrees[file.first]->SetBranchStatus("recoZmass", 1);
        mctrees[file.first]->SetBranchStatus("mTLepMet", 1);
        mctrees[file.first]->SetBranchStatus("nPU_mean", 1);
        mctrees[file.first]->SetBranchStatus("nLooseMuons", 1);
        mctrees[file.first]->SetBranchStatus("nLooseElectrons", 1);
        mctrees[file.first]->SetBranchStatus("nTightMuons", 1);
        mctrees[file.first]->SetBranchStatus("hlt_photon", 1);
        mctrees[file.first]->SetBranchStatus("hlt_dimuon", 1);
        mctrees[file.first]->SetBranchStatus("hlt_singlemu", 1);
        mctrees[file.first]->SetBranchStatus("hlt_razor", 1);
        mctrees[file.first]->SetBranchStatus("nBTaggedJets", 1);
        mctrees[file.first]->SetBranchStatus("genZpt", 1);
        mctrees[file.first]->SetBranchStatus("genWpt", 1);
        mctrees[file.first]->SetBranchStatus("bjet1PassMedium", 1);
        mctrees[file.first]->SetBranchStatus("bjet2PassMedium", 1);
        mctrees[file.first]->SetBranchStatus("bjet1Pt", 1);
        mctrees[file.first]->SetBranchStatus("bjet2Pt", 1);

        mctrees[file.first]->SetBranchAddress(Form("met%s", suffixes[file.first].c_str()), &mets[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("MR%s", suffixes[file.first].c_str()), &mrs[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("Rsq%s", suffixes[file.first].c_str()), &rsqs[file.first]);
        mctrees[file.first]->SetBranchAddress("MR", &MR);
        mctrees[file.first]->SetBranchAddress("Rsq", &Rsq);
        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("leadingMuonPt", &leadingMuonPt);
        mctrees[file.first]->SetBranchAddress("leadingMuonEta", &leadingMuonEta);
        mctrees[file.first]->SetBranchAddress("leadingTightMuonPt", &leadingTightMuonPt);
        mctrees[file.first]->SetBranchAddress("leadingTightMuonEta", &leadingTightMuonEta);
        mctrees[file.first]->SetBranchAddress("subleadingMuonPt", &subleadingMuonPt);
        mctrees[file.first]->SetBranchAddress("subleadingMuonEta", &subleadingMuonEta);
        mctrees[file.first]->SetBranchAddress("leadingPhotonPt", &leadingPhotonPt);
        mctrees[file.first]->SetBranchAddress("leadingPhotonEta", &leadingPhotonEta);
        mctrees[file.first]->SetBranchAddress("recoZpt", &recoZpt);
        mctrees[file.first]->SetBranchAddress("recoZeta", &recoZeta);
        mctrees[file.first]->SetBranchAddress("recoZmass", &recoZmass);
        mctrees[file.first]->SetBranchAddress("mTLepMet", &mTLepMet);
        mctrees[file.first]->SetBranchAddress("nPU_mean", &nPU_mean);
        mctrees[file.first]->SetBranchAddress("genZpt", &genZPt);
        mctrees[file.first]->SetBranchAddress("genWpt", &genWPt);
        mctrees[file.first]->SetBranchAddress("bjet1PassMedium", &bjet1PassMedium);
        mctrees[file.first]->SetBranchAddress("bjet2PassMedium", &bjet2PassMedium);
        mctrees[file.first]->SetBranchAddress("bjet1Pt", &bjet1Pt);
        mctrees[file.first]->SetBranchAddress("bjet2Pt", &bjet2Pt);
    }
    for(auto &file : datafiles){
        datatrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        datatrees[file.first]->SetBranchStatus("*", 0);

        datatrees[file.first]->SetBranchStatus("met_noW", 1);
        datatrees[file.first]->SetBranchStatus("MR_noW", 1);
        datatrees[file.first]->SetBranchStatus("Rsq_noW", 1);
        datatrees[file.first]->SetBranchStatus("numJets80_noW", 1);
        datatrees[file.first]->SetBranchStatus("met_noZ", 1);
        datatrees[file.first]->SetBranchStatus("MR_noZ", 1);
        datatrees[file.first]->SetBranchStatus("Rsq_noZ", 1);
        datatrees[file.first]->SetBranchStatus("numJets80_noZ", 1);
        datatrees[file.first]->SetBranchStatus("met_noPho", 1);
        datatrees[file.first]->SetBranchStatus("MR_noPho", 1);
        datatrees[file.first]->SetBranchStatus("Rsq_noPho", 1);
        datatrees[file.first]->SetBranchStatus("numJets80_noPho", 1);
        datatrees[file.first]->SetBranchStatus("met", 1);
        datatrees[file.first]->SetBranchStatus("MR", 1);
        datatrees[file.first]->SetBranchStatus("Rsq", 1);
        datatrees[file.first]->SetBranchStatus("numJets80", 1);
        datatrees[file.first]->SetBranchStatus("leadingMuonPt", 1);
        datatrees[file.first]->SetBranchStatus("leadingMuonEta", 1);
        datatrees[file.first]->SetBranchStatus("leadingTightMuonPt", 1);
        datatrees[file.first]->SetBranchStatus("leadingTightMuonEta", 1);
        datatrees[file.first]->SetBranchStatus("subleadingMuonPt", 1);
        datatrees[file.first]->SetBranchStatus("subleadingMuonEta", 1);
        datatrees[file.first]->SetBranchStatus("leadingPhotonPt", 1);
        datatrees[file.first]->SetBranchStatus("leadingPhotonEta", 1);
        datatrees[file.first]->SetBranchStatus("recoZpt", 1);
        datatrees[file.first]->SetBranchStatus("recoZeta", 1);
        datatrees[file.first]->SetBranchStatus("recoZmass", 1);
        datatrees[file.first]->SetBranchStatus("mTLepMet", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton50", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton75", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton90", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton135", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton150", 1);
        datatrees[file.first]->SetBranchStatus("nLooseMuons", 1);
        datatrees[file.first]->SetBranchStatus("nLooseElectrons", 1);
        datatrees[file.first]->SetBranchStatus("nTightMuons", 1);
        datatrees[file.first]->SetBranchStatus("hlt_photon", 1);
        datatrees[file.first]->SetBranchStatus("nBTaggedJets", 1);
        datatrees[file.first]->SetBranchStatus("hlt_dimuon", 1);
        datatrees[file.first]->SetBranchStatus("hlt_singlemu", 1);
        datatrees[file.first]->SetBranchStatus("hlt_razor", 1);
        datatrees[file.first]->SetBranchStatus("bjet1PassMedium", 1);
        datatrees[file.first]->SetBranchStatus("bjet2PassMedium", 1);
        datatrees[file.first]->SetBranchStatus("bjet1Pt", 1);
        datatrees[file.first]->SetBranchStatus("bjet2Pt", 1);
        datatrees[file.first]->SetBranchStatus("Flag_HBHENoiseFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_CSCTightHaloFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_eeBadScFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_ecalLaserCorrFilter", 1); // enable

        datatrees[file.first]->SetBranchAddress(Form("met%s", suffixes[file.first].c_str()), &mets[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("MR%s", suffixes[file.first].c_str()), &mrs[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("Rsq%s", suffixes[file.first].c_str()), &rsqs[file.first]);
        datatrees[file.first]->SetBranchAddress("leadingMuonPt", &leadingMuonPt);
        datatrees[file.first]->SetBranchAddress("leadingMuonEta", &leadingMuonEta);
        datatrees[file.first]->SetBranchAddress("leadingTightMuonPt", &leadingTightMuonPt);
        datatrees[file.first]->SetBranchAddress("leadingTightMuonEta", &leadingTightMuonEta);
        datatrees[file.first]->SetBranchAddress("subleadingMuonPt", &subleadingMuonPt);
        datatrees[file.first]->SetBranchAddress("subleadingMuonEta", &subleadingMuonEta);
        datatrees[file.first]->SetBranchAddress("leadingPhotonPt", &leadingPhotonPt);
        datatrees[file.first]->SetBranchAddress("leadingPhotonEta", &leadingPhotonEta);
        datatrees[file.first]->SetBranchAddress("recoZpt", &recoZpt);
        datatrees[file.first]->SetBranchAddress("recoZeta", &recoZeta);
        datatrees[file.first]->SetBranchAddress("recoZmass", &recoZmass);
        datatrees[file.first]->SetBranchAddress("mTLepMet", &mTLepMet);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton50", &passedHLTPhoton50);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton75", &passedHLTPhoton75);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton90", &passedHLTPhoton90);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton135", &passedHLTPhoton135);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton150", &passedHLTPhoton150);
        datatrees[file.first]->SetBranchAddress("bjet1PassMedium", &bjet1PassMedium);
        datatrees[file.first]->SetBranchAddress("bjet2PassMedium", &bjet2PassMedium);
        datatrees[file.first]->SetBranchAddress("bjet1Pt", &bjet1Pt);
        datatrees[file.first]->SetBranchAddress("bjet2Pt", &bjet2Pt);
        datatrees[file.first]->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter); // enable 
    }
    //luminosities collected by the various photon triggers
    float lumi_HLTPhoton50  = 1.353e0 + 4.921e0 + 7.947e0 + 8.131e0;
    float lumi_HLTPhoton75  = 8.111e0 + 2.953e1 + 4.768e1 + 4.879e1;
    float lumi_HLTPhoton90  = 1.622e1 + 6.408e1 + 1.010e2 + 9.948e1;
    float lumi_HLTPhoton135 = 8.893e2 + 1.476e2 + 5.429e3 + 7.318e3;
    float lumi_HLTPhoton150 = 8.893e2 + 4.429e3 + 7.152e3 + 7.318e3;

    //load efficiency/acceptance histograms
    TFile effFile("Run1LeptonPhotonEfficiency.root");
    TH2F muonLooseEffHisto = *(TH2F *)effFile.Get("MuonEfficiency");
    TH2F muonTightEffHisto = *(TH2F *)effFile.Get("MuonEfficiencyTight");
    TH2F photonEffHisto = *(TH2F *)effFile.Get("PhotonEfficiency");
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

    //Step 1: Get the distributions to reweight by: MET, MR, Rsq
    map<string, TH1F> metHistosForReweighting;
    map<string, TH2F> razorHistosForReweighting;
    map<string, TH1F> MRHistosBeforeReweighting;
    map<string, TH1F> RsqHistosBeforeReweighting;
    TH1F mcPhotonPt("mcPhotonPt", "Photon Pt; photon pt", 100, 0, 1000);
    mcPhotonPt.Sumw2();
    for(auto &tree : mctrees){
        cout << "Filling MC histograms: " << tree.first << endl;
        metHistosForReweighting[tree.first] = TH1F(Form("metmc%s", tree.first.c_str()), "MET (GeV); MET(GeV)", nMetBins, MetMin, MetMax);
        razorHistosForReweighting[tree.first] = TH2F(Form("razormc%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
        // MRHistosBeforeReweighting[tree.first] = TH1F(Form("mrmc%s", tree.first.c_str()), "; MR (GeV)", nMRBins, MRBinLowEdges);
        // RsqHistosBeforeReweighting[tree.first] = TH1F(Form("rsqmc%s", tree.first.c_str()), "; Rsq", nRsqBins, RsqBinLowEdges);
        MRHistosBeforeReweighting[tree.first] = TH1F(Form("mrmc%s", tree.first.c_str()), "; MR (GeV)", 20, 300, 4000);
        RsqHistosBeforeReweighting[tree.first] = TH1F(Form("rsqmc%s", tree.first.c_str()), "; Rsq",  10, 0.15, 1.5);
        MRHistosBeforeReweighting[tree.first].Sumw2();
        RsqHistosBeforeReweighting[tree.first].Sumw2();
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
            eventWeight *= pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(nPU_mean));

            //reweight according to selection efficiency and acceptance
            if(tree.first == "GJets" || tree.first == "QCD" || tree.first == "WG" || tree.first == "ZG" || tree.first == "TTG"){
                if(!computeDataOverMCSFs){
                    double effFactor = photonEffHisto.GetBinContent(photonEffHisto.FindFixBin(min(leadingPhotonPt, maxPhotonPt), fabs(leadingPhotonEta)));

                    if(effFactor > 1e-5) eventWeight /= effFactor;
                    else{ 
                        eventWeight = 0;
                        //cout << "Warning: efficiency histogram gives 0 (pt " << leadingPhotonPt << ", eta " << leadingPhotonEta << "); setting event weight to 0" << endl;
                    }
                }
                mcPhotonPt.Fill(leadingPhotonPt, eventWeight);
            }
            else if(tree.first == "WJets" || tree.first == "Top" || tree.first == "TTJets" || tree.first == "TTW" || tree.first == "TTZ"){
                //ID efficiency
                if(!computeDataOverMCSFs){
                    double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindFixBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));

                    if(effFactor > 1e-5) eventWeight /= effFactor;
                    else{ 
                        eventWeight = 0;
                        //cout << "Warning: efficiency histogram gives 0 (pt " << leadingTightMuonPt << ", eta " << leadingTightMuonEta << "); setting event weight to 0" << endl;
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
		
		//ISR reweighting for WJets
                // if(tree.first == "WJets"){
                //     double isrWeight = 1.0;
                //     if(genWPt > 120 && genWPt < 150){
                //         isrWeight = 0.95;
                //     }
                //     else if(genWPt > 150 && genWPt < 250){
                //         isrWeight = 0.90;
                //     }
                //     else if(genWPt > 250){
                //         isrWeight = 0.80;
                //     }
                //     eventWeight *= isrWeight;
                // }


                //TTJets SF
                if(tree.first == "TTJets"){
                    double ttjetsSF = TTBarDileptonScaleFactor->GetBinContent(TTBarDileptonScaleFactor->FindFixBin(min(MR, maxMRForTTJetsSF), min(Rsq, maxRsqForTTJetsSF)));
                    if(ttjetsSF > 1e-5) eventWeight *= ttjetsSF;
                }
            }
            else if(tree.first == "DYJets" || tree.first == "TopDY" || tree.first == "TTJetsDY" || tree.first == "TTWDY" || tree.first == "TTZDY"){
                if(!computeDataOverMCSFs){
                    double effFactor = muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindFixBin(min(leadingMuonPt, maxMuonPt), fabs(leadingMuonEta)));
                    effFactor *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindFixBin(min(subleadingMuonPt, maxMuonPt), fabs(subleadingMuonEta)));
                    effFactor *= zAccHisto.GetBinContent(zAccHisto.FindFixBin(min(recoZpt, maxZPt), fabs(recoZeta)));

                    if(effFactor > 1e-5) eventWeight /= effFactor;
                    else{
                        //cout << "Warning: efficiency histogram gives 0; (lead pt " << leadingMuonPt << ", leading eta " << leadingMuonEta << ", subleading pt " << subleadingMuonPt << ", subleading eta " << subleadingMuonEta << ", z pt " << recoZpt << ", z eta " << recoZeta << "); setting event weight to 0" << endl;
                        eventWeight = 0;
                    }
                }
                //trigger and normalization scale factors
                eventWeight *= doubleMuTriggerSF;
                eventWeight *= doubleMuNormalizationSF;

                //ISR reweighting for DYJets
                // if(tree.first == "DYJets"){
                //     double isrWeight = 1.0;
                //     if(genZPt > 120 && genZPt < 150){
                //         isrWeight = 0.95;
                //     }
                //     else if(genZPt > 150 && genZPt < 250){
                //         isrWeight = 0.90;
                //     }
                //     else if(genZPt > 250){
                //         isrWeight = 0.80;
                //     }
                //     eventWeight *= isrWeight;
                // }

                //TTJets SF
                if(tree.first == "TTJetsDY"){
                    double ttjetsSF = TTBarDileptonScaleFactor->GetBinContent(TTBarDileptonScaleFactor->FindFixBin(min(MR, maxMRForTTJetsSF), min(Rsq, maxRsqForTTJetsSF)));
                    if(ttjetsSF > 1e-5) eventWeight *= ttjetsSF;
                }
            }
            else if(tree.first != "ZJets"){
                cout << "Warning: unexpected dataset " << tree.first << " encountered!  Check the code." << endl;
            }
            //fill each quantity
            metHistosForReweighting[tree.first].Fill(mets[tree.first], eventWeight);
            razorHistosForReweighting[tree.first].Fill(mrs[tree.first], rsqs[tree.first], eventWeight);
            MRHistosBeforeReweighting[tree.first].Fill(mrs[tree.first], eventWeight);
            RsqHistosBeforeReweighting[tree.first].Fill(rsqs[tree.first], eventWeight);
        }
    }

    //Step 2: Sanity check: apply the reweighting factors to MC
    map<string, TH2F> razorHistosMC;
    if(!computeDataOverMCSFs){
        for(auto &tree : mctrees){
            cout << "Filling reweighted MC histograms: " << tree.first << endl;
            razorHistosMC[tree.first] = TH2F(Form("razorMCReweighted%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
            razorHistosMC[tree.first].Sumw2();
            uint nEntries = tree.second->GetEntries();
            TTreeFormula cutsFormula(Form("%sCutsFormula", tree.first.c_str()), cuts[tree.first].c_str(), tree.second);
            cutsFormula.GetNdata();
            for(uint i = 0; i < nEntries; i++){
                //get entry
                tree.second->GetEntry(i);

                //apply selection cuts
                bool passesSelection = cutsFormula.EvalInstance();
                if(!passesSelection) continue;

                float eventWeight = weight;
                eventWeight *= pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(nPU_mean));

                //reweight by efficiency and acceptance
                if(tree.first == "GJets" || tree.first == "QCD" || tree.first == "TTG" || tree.first == "WG" || tree.first == "ZG"){
                    if(!computeDataOverMCSFs){
                        double effFactor = photonEffHisto.GetBinContent(photonEffHisto.FindFixBin(min(leadingPhotonPt, maxPhotonPt), fabs(leadingPhotonEta)));
                        if(effFactor > 1e-5) eventWeight /= effFactor;
                        else{
                            eventWeight = 0;
                            //cout << "Warning: efficiency histogram gives 0 (pt " << leadingPhotonPt << ", eta " << leadingPhotonEta << "); setting event weight to 0" << endl;
                        }
                    }
                }
                else if(tree.first == "WJets" || tree.first == "Top" || tree.first == "TTJets" || tree.first == "TTW" || tree.first == "TTZ"){
                    if(!computeDataOverMCSFs){
                        double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindFixBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));

                        if(effFactor > 1e-5) eventWeight /= effFactor;
                        else{ 
                            eventWeight = 0;
                            //cout << "Warning: efficiency histogram gives 0 (pt " << leadingTightMuonPt << ", eta " << leadingTightMuonEta << "); setting event weight to 0" << endl;
                        }
                    }
                    //trigger scale factor
                    eventWeight *= singleMuTriggerSF;
                }
                else if(tree.first == "DYJets" || tree.first == "TopDY" || tree.first == "TTJetsDY" || tree.first == "TTWDY" || tree.first == "TTZDY"){
                    if(!computeDataOverMCSFs){
                        double effFactor = muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindFixBin(min(leadingMuonPt, maxMuonPt), fabs(leadingMuonEta)));
                        effFactor *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindFixBin(min(subleadingMuonPt, maxMuonPt), fabs(subleadingMuonEta)));
                        effFactor *= zAccHisto.GetBinContent(zAccHisto.FindFixBin(min(recoZpt, maxZPt), fabs(recoZeta)));

                        if(effFactor > 1e-5) eventWeight /= effFactor;
                        else{
                            //cout << "Warning: efficiency histogram gives 0; (lead pt " << leadingMuonPt << ", leading eta " << leadingMuonEta << ", subleading pt " << subleadingMuonPt << ", subleading eta " << subleadingMuonEta << ", z pt " << recoZpt << ", z eta " << recoZeta << "); setting event weight to 0" << endl;
                            eventWeight = 0;
                        }
                    }
                    //trigger and normalization scale factors
                    eventWeight *= doubleMuTriggerSF;
                    eventWeight *= doubleMuNormalizationSF;
                }
                else if(tree.first != "ZJets"){
                    cout << "Warning: unexpected dataset " << tree.first << " encountered!  Check the code." << endl;
                }

                if(reweightByRazor){ //reweight by MR and Rsq
                    //get the factor to reweight by
                    float denominator = razorHistosForReweighting[tree.first].GetBinContent(razorHistosForReweighting[tree.first].FindFixBin(mrs[tree.first], rsqs[tree.first]));
                    float numerator = razorHistosForReweighting[normalizeTo].GetBinContent(razorHistosForReweighting[normalizeTo].FindFixBin(mrs[tree.first], rsqs[tree.first]));
                    if(denominator > 0){
                        eventWeight *= numerator / denominator;
                    }
                } 
                else{ //reweight by MET
                    //get the factor to reweight by
                    float denominator = metHistosForReweighting[tree.first].GetBinContent(metHistosForReweighting[tree.first].FindFixBin(mets[tree.first]));
                    float numerator = metHistosForReweighting[normalizeTo].GetBinContent(metHistosForReweighting[normalizeTo].FindFixBin(mets[tree.first]));
                    if(denominator > 0){
                        eventWeight *= numerator / denominator;    
                    }
                }
                razorHistosMC[tree.first].Fill(mrs[tree.first], rsqs[tree.first], eventWeight);
            }
        }
    }

    //Step 3: make data distributions
    map<string, TH2F> razorHistosDataBeforeReweighting; //apply only efficiency and acceptance corrections
    map<string, TH1F> MRHistosDataBeforeReweighting;
    map<string, TH1F> RsqHistosDataBeforeReweighting;
    TH1F dataPhotonPt("dataPhotonPt", "Photon Pt; photon pt", 100, 0, 1000);
    dataPhotonPt.Sumw2();
    for(auto &tree : datatrees){
        cout << "Filling data histograms: " << tree.first << endl;
        razorHistosDataBeforeReweighting[tree.first] = TH2F(Form("razordatabeforereweighting%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
        // MRHistosDataBeforeReweighting[tree.first] = TH1F(Form("mrdata%s", tree.first.c_str()), "; MR (GeV)", nMRBins, MRBinLowEdges);
        // RsqHistosDataBeforeReweighting[tree.first] = TH1F(Form("rsqdata%s", tree.first.c_str()), "; Rsq (GeV)", nRsqBins, RsqBinLowEdges);
        MRHistosDataBeforeReweighting[tree.first] = TH1F(Form("mrdata%s", tree.first.c_str()), "; MR (GeV)", 20, 300, 4000);
        RsqHistosDataBeforeReweighting[tree.first] = TH1F(Form("rsqdata%s", tree.first.c_str()), "; Rsq (GeV)", 10, 0.15, 1.5);
        razorHistosDataBeforeReweighting[tree.first].Sumw2();
        MRHistosDataBeforeReweighting[tree.first].Sumw2();
        RsqHistosDataBeforeReweighting[tree.first].Sumw2();
        uint nEntries = tree.second->GetEntries();
        TTreeFormula cutsFormula(Form("%sCutsFormula", tree.first.c_str()), cuts[tree.first].c_str(), tree.second);
        cutsFormula.GetNdata();
        for(uint i = 0; i < nEntries; i++){
            //get entry
            tree.second->GetEntry(i);

            //apply selection cuts
            bool passesSelection = cutsFormula.EvalInstance();
            if(!passesSelection) continue;

	    if(!Flag_HBHENoiseFilter || !Flag_CSCTightHaloFilter || !Flag_eeBadScFilter ) continue;

            float eventWeight = 1.0;
            //reweight by efficiency and acceptance
            if(tree.first == "GJets"){
                if(!computeDataOverMCSFs){
                    double effFactor = photonEffHisto.GetBinContent(photonEffHisto.FindFixBin(min(leadingPhotonPt, maxPhotonPt), fabs(leadingPhotonEta)));
                    if(effFactor > 1e-5) eventWeight /= effFactor;
                    else{ 
                        eventWeight = 0;
                        //cout << "Warning: efficiency histogram gives 0 (pt " << leadingPhotonPt << ", eta " << leadingPhotonEta << "); setting event weight to 0" << endl;
                    }
                }
                double triggerWeightRestricted = 0.0;
                //get weight if associate each photon trigger with a particular pt range
                if(passedHLTPhoton150 && leadingPhotonPt > 165){ 
                    triggerWeightRestricted = 1.0;
                }
                else if(passedHLTPhoton135 && leadingPhotonPt > 150 && leadingPhotonPt < 165){
                    triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton135;
                }
                else if(passedHLTPhoton90 && leadingPhotonPt > 100 && leadingPhotonPt < 150){
                    triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton90;
                }
                else if(passedHLTPhoton75 && leadingPhotonPt > 90 && leadingPhotonPt < 100){
                    triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton75; 
                }
                else if(passedHLTPhoton50 && leadingPhotonPt < 90){
                    triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton50;
                }
                eventWeight *= triggerWeightRestricted;

		if(leadingPhotonPt>5000) continue; // reject noise

                dataPhotonPt.Fill(leadingPhotonPt, eventWeight);
            }
            else if(tree.first == "WJets"){
                if(!computeDataOverMCSFs){
                    double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindFixBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));
                    if(effFactor > 1e-5) eventWeight /= effFactor;
                    else{ 
                        eventWeight = 0;
                        //cout << "Warning: efficiency histogram gives 0 (pt " << leadingTightMuonPt << ", eta " << leadingTightMuonEta << "); setting event weight to 0" << endl;
                    }
                }
            }
            else if(tree.first == "DYJets"){
                if(!computeDataOverMCSFs){
                    double effFactor = muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindFixBin(min(leadingMuonPt, maxMuonPt), fabs(leadingMuonEta)));
                    effFactor *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindFixBin(min(subleadingMuonPt, maxMuonPt), fabs(subleadingMuonEta)));
                    effFactor *= zAccHisto.GetBinContent(zAccHisto.FindFixBin(min(recoZpt, maxZPt), fabs(recoZeta)));
                    if(effFactor > 1e-5) eventWeight /= effFactor;
                    else{
                        //cout << "Warning: efficiency histogram gives 0; (lead pt " << leadingMuonPt << ", leading eta " << leadingMuonEta << ", subleading pt " << subleadingMuonPt << ", subleading eta " << subleadingMuonEta << ", z pt " << recoZpt << ", z eta " << recoZeta << "); setting event weight to 0" << endl;
                        eventWeight = 0;
                    }
                }
            }
            else if(tree.first != "ZJets"){
                cout << "Warning: unexpected dataset " << tree.first << " encountered!  Check the code." << endl;
            }
            razorHistosDataBeforeReweighting[tree.first].Fill(mrs[tree.first], rsqs[tree.first], eventWeight);
            MRHistosDataBeforeReweighting[tree.first].Fill(mrs[tree.first], eventWeight);
            RsqHistosDataBeforeReweighting[tree.first].Fill(rsqs[tree.first], eventWeight);
        }
    }
    //for rare background processes, include a 20% uncertainty on the total yield in each bin, summed in quadrature with the statistical uncertainty
    double sysErrorFrac = 0.2;
    //for QCD, assign a 100% uncertainty
    double qcdErrorFrac = 1.0;
    for(auto &tree : mctrees){
        //only do this for rare processes (not signal processes or TTJets)
        if(tree.first == "DYJets" || tree.first == "WJets" || tree.first == "GJets" || tree.first == "ZJets" || tree.first == "TTJets") continue; 
        for(int i = 0; i < razorHistosForReweighting[tree.first].GetNbinsX()+1; i++){
            for(int j = 0; j < razorHistosForReweighting[tree.first].GetNbinsY()+1; j++){
                double error = 0.0;
                if(tree.first == "QCD"){
                    error = qcdErrorFrac*razorHistosForReweighting[tree.first].GetBinContent(i, j);
                }
                else{
                    error = sysErrorFrac*razorHistosForReweighting[tree.first].GetBinContent(i, j);
                }
                razorHistosForReweighting[tree.first].SetBinError(i, j, sqrt(pow(razorHistosForReweighting[tree.first].GetBinError(i, j), 2) + error*error));
            }
        }
    }
    //subtract the top background from the WJets histogram
    razorHistosDataBeforeReweighting["WJets"] = razorHistosDataBeforeReweighting["WJets"] - razorHistosForReweighting["Top"];
    razorHistosDataBeforeReweighting["WJets"] = razorHistosDataBeforeReweighting["WJets"] - razorHistosForReweighting["TTJets"];
    razorHistosDataBeforeReweighting["WJets"] = razorHistosDataBeforeReweighting["WJets"] - razorHistosForReweighting["TTW"];
    razorHistosDataBeforeReweighting["WJets"] = razorHistosDataBeforeReweighting["WJets"] - razorHistosForReweighting["TTZ"];
    //subtract the top background from the DYJets histogram
    razorHistosDataBeforeReweighting["DYJets"] = razorHistosDataBeforeReweighting["DYJets"] - razorHistosForReweighting["TTJetsDY"];
    razorHistosDataBeforeReweighting["DYJets"] = razorHistosDataBeforeReweighting["DYJets"] - razorHistosForReweighting["TopDY"];
    razorHistosDataBeforeReweighting["DYJets"] = razorHistosDataBeforeReweighting["DYJets"] - razorHistosForReweighting["TTWDY"];
    razorHistosDataBeforeReweighting["DYJets"] = razorHistosDataBeforeReweighting["DYJets"] - razorHistosForReweighting["TTZDY"];
    //subtract the QCD background from the photon+jets histogram
    razorHistosDataBeforeReweighting["GJets"] = razorHistosDataBeforeReweighting["GJets"] - razorHistosForReweighting["QCD"];
    razorHistosDataBeforeReweighting["GJets"] = razorHistosDataBeforeReweighting["GJets"] - razorHistosForReweighting["TTG"];
    razorHistosDataBeforeReweighting["GJets"] = razorHistosDataBeforeReweighting["GJets"] - razorHistosForReweighting["WG"];
    razorHistosDataBeforeReweighting["GJets"] = razorHistosDataBeforeReweighting["GJets"] - razorHistosForReweighting["ZG"];

    //Step 4: apply reweighting factors to data
    map<string, TH2F> razorHistosData;
    if(!computeDataOverMCSFs){
        if(reweightByRazor){
            for(auto &tree : datatrees){
                cout << "Making weighted data histograms: " << tree.first << endl;
                razorHistosData[tree.first] = TH2F(Form("razordata%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
                razorHistosData[tree.first].Sumw2();
                for(int i = 0; i < razorHistosData[tree.first].GetNbinsX()+1; i++){
                    for(int j = 0; j < razorHistosData[tree.first].GetNbinsY()+1; j++){
                        float numerator = razorHistosForReweighting[normalizeTo].GetBinContent(i, j);
                        float denominator = razorHistosForReweighting[tree.first].GetBinContent(i, j);
                        float numeratorError = razorHistosForReweighting[normalizeTo].GetBinError(i, j);
                        float denominatorError = razorHistosForReweighting[tree.first].GetBinError(i, j);
                        if(denominator > 0){
                            razorHistosData[tree.first].SetBinContent(i, j, razorHistosDataBeforeReweighting[tree.first].GetBinContent(i, j)*numerator/denominator);
                            //compute the uncertainty on the bin, folding in uncertainties on the scale factor numerator/denominator
                            razorHistosData[tree.first].SetBinError(i, j, sqrt(pow(razorHistosDataBeforeReweighting[tree.first].GetBinError(i, j)*numerator/denominator, 2) + pow(razorHistosDataBeforeReweighting[tree.first].GetBinContent(i, j)*numeratorError/denominator, 2) + pow(razorHistosDataBeforeReweighting[tree.first].GetBinContent(i, j)*numerator*denominatorError/(denominator*denominator), 2)));
                        }
                        else{
                            razorHistosData[tree.first].SetBinContent(i, j, 0.);
                            razorHistosData[tree.first].SetBinError(i, j, 0.);
                        }
                    }
                }   
            }
        }
        else{ 
            for(auto &tree : datatrees){
                cout << "Making weighted data histograms: " << tree.first << endl;
                razorHistosData[tree.first] = TH2F(Form("razordata%s", tree.first.c_str()), "; MR (GeV); Rsq", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges);
                razorHistosData[tree.first].Sumw2();

                uint nEntries = tree.second->GetEntries();
                TTreeFormula cutsFormula(Form("%sCutsFormula", tree.first.c_str()), cuts[tree.first].c_str(), tree.second);
                cutsFormula.GetNdata();
                for(uint i = 0; i < nEntries; i++){
                    //get entry
                    tree.second->GetEntry(i);

                    //apply selection cuts
                    bool passesSelection = cutsFormula.EvalInstance();
                    if(!passesSelection) continue;

                    float reweightFactor = 1.0;
                    //reweight by efficiency and acceptance
                    if(tree.first == "GJets"){
                        double effFactor = photonEffHisto.GetBinContent(photonEffHisto.FindFixBin(min(leadingPhotonPt, maxPhotonPt), fabs(leadingPhotonEta)));
                        if(effFactor > 1e-5) reweightFactor /= effFactor;
                        else{ 
                            reweightFactor = 0;
                            //cout << "Warning: efficiency histogram gives 0 (pt " << leadingPhotonPt << ", eta " << leadingPhotonEta << "); setting event weight to 0" << endl;
                        }
                        //get weight if associate each photon trigger with a particular pt range
                        double triggerWeightRestricted = 0.0;
                        if(passedHLTPhoton150 && leadingPhotonPt > 165){ 
                            triggerWeightRestricted = 1.0;
                        }
                        else if(passedHLTPhoton135 && leadingPhotonPt > 150 && leadingPhotonPt < 165){
                            triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton135;
                        }
                        else if(passedHLTPhoton90 && leadingPhotonPt > 100 && leadingPhotonPt < 150){
                            triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton90;
                        }
                        else if(passedHLTPhoton75 && leadingPhotonPt > 90 && leadingPhotonPt < 100){
                            triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton75; 
                        }
                        else if(passedHLTPhoton50 && leadingPhotonPt < 90){
                            triggerWeightRestricted = lumi_HLTPhoton150/lumi_HLTPhoton50;
                        }
                        reweightFactor *= triggerWeightRestricted;
                    }
                    else if(tree.first == "WJets"){
                        double effFactor = muonTightEffHisto.GetBinContent(muonTightEffHisto.FindFixBin(min(leadingTightMuonPt, maxMuonPt), fabs(leadingTightMuonEta)));
                        if(effFactor > 1e-5) reweightFactor /= effFactor;
                        else{ 
                            reweightFactor = 0;
                            //cout << "Warning: efficiency histogram gives 0 (pt " << leadingTightMuonPt << ", eta " << leadingTightMuonEta << "); setting event weight to 0" << endl;
                        }
                    }
                    else if(tree.first == "DYJets"){
                        double effFactor = muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindFixBin(min(leadingMuonPt, maxMuonPt), fabs(leadingMuonEta)));
                        effFactor *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindFixBin(min(subleadingMuonPt, maxMuonPt), fabs(subleadingMuonEta)));
                        effFactor *= zAccHisto.GetBinContent(zAccHisto.FindFixBin(min(recoZpt, maxZPt), fabs(recoZeta)));
                        if(effFactor > 1e-5) reweightFactor /= effFactor;
                        else{
                            //cout << "Warning: efficiency histogram gives 0; (lead pt " << leadingMuonPt << ", leading eta " << leadingMuonEta << ", subleading pt " << subleadingMuonPt << ", subleading eta " << subleadingMuonEta << ", z pt " << recoZpt << ", z eta " << recoZeta << "); setting event weight to 0" << endl;
                            reweightFactor = 0;
                        }
                    }
                    //reweight by MET
                    //get the factor to reweight by
                    float denominator = metHistosForReweighting[tree.first].GetBinContent(metHistosForReweighting[tree.first].FindFixBin(mets[tree.first]));
                    float numerator = metHistosForReweighting[normalizeTo].GetBinContent(metHistosForReweighting[normalizeTo].FindFixBin(mets[tree.first]));
                    if(denominator > 0){
                        reweightFactor *= numerator / denominator;    
                    }

                    razorHistosData[tree.first].Fill(mrs[tree.first], rsqs[tree.first], reweightFactor);
                }
            }
        }
    }

    TCanvas c("c", "c", 800, 600);
    c.SetLogx();
    //print "step 1" histograms used for reweighting
    c.SetLogz();
    for(auto &hist : razorHistosForReweighting){
        if(computeDataOverMCSFs) hist.second.SetTitle(Form("MC for %s", hist.first.c_str()));
        else hist.second.SetTitle(Form("MC reweighted by efficiency and acceptance, for %s", hist.first.c_str()));
        hist.second.GetXaxis()->SetTitle("MR");
        hist.second.GetYaxis()->SetTitle("Rsq");
        hist.second.SetStats(0);
        hist.second.Draw("colz");
        hist.second.Draw("same,text");
        c.Print(Form("controlSampleMCHistogram%s_noISR.pdf", hist.first.c_str()));
        // c.Print(Form("controlSampleMCHistogram%s.root", hist.first.c_str()));
    }
    //print "step 2" histograms
    if(!computeDataOverMCSFs){
        for(auto &hist : razorHistosMC){
            hist.second.SetTitle(Form("MC with all corrections applied, for %s", hist.first.c_str()));
            hist.second.GetXaxis()->SetTitle("MR");
            hist.second.GetYaxis()->SetTitle("Rsq");
            hist.second.SetStats(0);
            hist.second.Draw("colz");
            hist.second.Draw("same,text");
            c.Print(Form("controlSampleReweightedMCHistogram%s_noISR.pdf", hist.first.c_str()));
            // c.Print(Form("controlSampleReweightedMCHistogram%s.root", hist.first.c_str()));
        }
    }
    //print "step 4" razor histograms from data
    if(!computeDataOverMCSFs){
        for(auto &hist : razorHistosData){
            hist.second.SetTitle(Form("Prediction for %s", hist.first.c_str()));
            hist.second.GetXaxis()->SetTitle("MR");
            hist.second.GetYaxis()->SetTitle("Rsq");
            hist.second.SetStats(0);
            hist.second.Draw("colz");
            hist.second.Draw("same,text");
            c.Print(Form("controlSampleHistogram%s_noISR.pdf", hist.first.c_str()));
            // c.Print(Form("controlSampleHistogram%s.root", hist.first.c_str()));
        }
    }
    //print razor histograms from data before reweighting by MR/Rsq/MET
    for(auto &hist : razorHistosDataBeforeReweighting){
        hist.second.SetTitle(Form("Data (background subtracted) for %s", hist.first.c_str()));
        hist.second.GetXaxis()->SetTitle("MR");
        hist.second.GetYaxis()->SetTitle("Rsq");
        hist.second.SetStats(0);
        hist.second.Draw("colz");
        hist.second.Draw("same,text");
        c.Print(Form("controlSampleHistogramBeforeReweighting%s_noISR.pdf", hist.first.c_str()));
        // c.Print(Form("controlSampleHistogramBeforeReweighting%s.root", hist.first.c_str()));
    }
    //print MR histograms, comparing data to MC
    c.SetLogy();
    //WJets
    THStack SingleMuonMC("SingleMuonMC", "MR in 1-muon control sample");
    MRHistosBeforeReweighting["WJets"].SetFillColor(kOrange+10);
    MRHistosBeforeReweighting["Top"].SetFillColor(kViolet-5);
    MRHistosBeforeReweighting["TTJets"].SetFillColor(kViolet-6);
    MRHistosBeforeReweighting["TTW"].SetFillColor(kRed+2);
    MRHistosBeforeReweighting["TTZ"].SetFillColor(kOrange-3);
    SingleMuonMC.Add(&MRHistosBeforeReweighting["TTZ"]);
    SingleMuonMC.Add(&MRHistosBeforeReweighting["TTW"]);
    SingleMuonMC.Add(&MRHistosBeforeReweighting["Top"]);
    SingleMuonMC.Add(&MRHistosBeforeReweighting["TTJets"]);
    SingleMuonMC.Add(&MRHistosBeforeReweighting["WJets"]);
    MRHistosDataBeforeReweighting["WJets"].SetMarkerStyle(20);
    MRHistosDataBeforeReweighting["WJets"].SetMarkerSize(1);
    TLegend *SingleMuonLegend = new TLegend(0.7, 0.7, 0.9, 0.9);
    SingleMuonLegend->AddEntry(&MRHistosBeforeReweighting["WJets"], "WJets MC");
    SingleMuonLegend->AddEntry(&MRHistosBeforeReweighting["TTJets"], "TTJets MC");
    SingleMuonLegend->AddEntry(&MRHistosBeforeReweighting["Top"], "Single Top MC");
    SingleMuonLegend->AddEntry(&MRHistosBeforeReweighting["TTW"], "TTW MC");
    SingleMuonLegend->AddEntry(&MRHistosBeforeReweighting["TTZ"], "TTZ MC");
    SingleMuonLegend->AddEntry(&MRHistosDataBeforeReweighting["WJets"], "2012 Data, Single Muon CS");
    DrawDataVsMCRatioPlot(&MRHistosDataBeforeReweighting["WJets"], &SingleMuonMC, SingleMuonLegend, "MR (GeV)", "controlSampleMRBackgroundSingleMuon", true);
    //DYJets
    THStack DoubleMuonMC("DoubleMuonMC", "MR in 2-muon control sample");
    MRHistosBeforeReweighting["DYJets"].SetFillColor(kAzure);
    MRHistosBeforeReweighting["TopDY"].SetFillColor(kViolet-5);
    MRHistosBeforeReweighting["TTJetsDY"].SetFillColor(kViolet-6);
    MRHistosBeforeReweighting["TTWDY"].SetFillColor(kRed+2);
    MRHistosBeforeReweighting["TTZDY"].SetFillColor(kOrange-3);
    DoubleMuonMC.Add(&MRHistosBeforeReweighting["TTZDY"]);
    DoubleMuonMC.Add(&MRHistosBeforeReweighting["TTWDY"]);
    DoubleMuonMC.Add(&MRHistosBeforeReweighting["TopDY"]);
    DoubleMuonMC.Add(&MRHistosBeforeReweighting["TTJetsDY"]);
    DoubleMuonMC.Add(&MRHistosBeforeReweighting["DYJets"]);
    MRHistosDataBeforeReweighting["DYJets"].SetMarkerStyle(20);
    MRHistosDataBeforeReweighting["DYJets"].SetMarkerSize(1);
    TLegend *DoubleMuonLegend = new TLegend(0.7, 0.7, 0.9, 0.9);
    DoubleMuonLegend->AddEntry(&MRHistosBeforeReweighting["DYJets"], "DYJets MC");
    DoubleMuonLegend->AddEntry(&MRHistosBeforeReweighting["TTJetsDY"], "TTJets MC");
    DoubleMuonLegend->AddEntry(&MRHistosBeforeReweighting["TopDY"], "Single Top MC");
    DoubleMuonLegend->AddEntry(&MRHistosBeforeReweighting["TTWDY"], "TTW MC");
    DoubleMuonLegend->AddEntry(&MRHistosBeforeReweighting["TTZDY"], "TTZ MC");    
    DoubleMuonLegend->AddEntry(&MRHistosDataBeforeReweighting["DYJets"], "2012 Data, Double Muon CS");
    DrawDataVsMCRatioPlot(&MRHistosDataBeforeReweighting["DYJets"], &DoubleMuonMC, DoubleMuonLegend, "MR (GeV)", "controlSampleMRBackgroundDoubleMuon", true);
    //Gamma+Jets
    THStack PhotonMC("PhotonMC", "MR in photon control sample");
    MRHistosBeforeReweighting["GJets"].SetFillColor(9);
    MRHistosBeforeReweighting["QCD"].SetFillColor(8);
    MRHistosBeforeReweighting["TTG"].SetFillColor(7);
    MRHistosBeforeReweighting["WG"].SetFillColor(38);
    MRHistosBeforeReweighting["ZG"].SetFillColor(4);
    PhotonMC.Add(&MRHistosBeforeReweighting["ZG"]);
    PhotonMC.Add(&MRHistosBeforeReweighting["WG"]);
    PhotonMC.Add(&MRHistosBeforeReweighting["TTG"]);
    PhotonMC.Add(&MRHistosBeforeReweighting["QCD"]);
    PhotonMC.Add(&MRHistosBeforeReweighting["GJets"]);
    MRHistosDataBeforeReweighting["GJets"].SetMarkerStyle(20);
    MRHistosDataBeforeReweighting["GJets"].SetMarkerSize(1);
    TLegend *PhotonLegend = new TLegend(0.7, 0.7, 0.9, 0.9);
    PhotonLegend->AddEntry(&MRHistosBeforeReweighting["GJets"], "GJets MC");
    PhotonLegend->AddEntry(&MRHistosBeforeReweighting["QCD"], "QCD MC");
    PhotonLegend->AddEntry(&MRHistosBeforeReweighting["TTG"], "TTG MC");
    PhotonLegend->AddEntry(&MRHistosBeforeReweighting["WG"], "WG MC");
    PhotonLegend->AddEntry(&MRHistosBeforeReweighting["ZG"], "ZG MC");
    PhotonLegend->AddEntry(&MRHistosDataBeforeReweighting["GJets"], "2012 Data, Photon CS");
    PhotonLegend->Draw();
    DrawDataVsMCRatioPlot(&MRHistosDataBeforeReweighting["GJets"], &PhotonMC, PhotonLegend, "MR (GeV)", "controlSampleMRBackgroundPhoton", true);

    //print Rsq histograms, comparing data to MC
    c.SetLogy();
    c.SetLogx(kFALSE);
    //WJets
    THStack SingleMuonRsqMC("SingleMuonMC", "Rsq in 1-muon control sample");
    RsqHistosBeforeReweighting["WJets"].SetFillColor(kOrange+10);
    RsqHistosBeforeReweighting["Top"].SetFillColor(kViolet-5);
    RsqHistosBeforeReweighting["TTJets"].SetFillColor(kViolet-6);
    RsqHistosBeforeReweighting["TTW"].SetFillColor(kRed+2);
    RsqHistosBeforeReweighting["TTZ"].SetFillColor(kOrange-3);
    SingleMuonRsqMC.Add(&RsqHistosBeforeReweighting["TTZ"]);
    SingleMuonRsqMC.Add(&RsqHistosBeforeReweighting["TTW"]);
    SingleMuonRsqMC.Add(&RsqHistosBeforeReweighting["Top"]);
    SingleMuonRsqMC.Add(&RsqHistosBeforeReweighting["TTJets"]);
    SingleMuonRsqMC.Add(&RsqHistosBeforeReweighting["WJets"]);
    RsqHistosDataBeforeReweighting["WJets"].SetMarkerStyle(20);
    RsqHistosDataBeforeReweighting["WJets"].SetMarkerSize(1);
    DrawDataVsMCRatioPlot(&RsqHistosDataBeforeReweighting["WJets"], &SingleMuonRsqMC, SingleMuonLegend, "Rsq", "controlSampleRsqBackgroundSingleMuon", false);
    //DYJets
    THStack DoubleMuonRsqMC("DoubleMuonMC", "Rsq in 2-muon control sample");
    RsqHistosBeforeReweighting["DYJets"].SetFillColor(kAzure);
    RsqHistosBeforeReweighting["TopDY"].SetFillColor(kViolet-5);
    RsqHistosBeforeReweighting["TTJetsDY"].SetFillColor(kViolet-6);
    RsqHistosBeforeReweighting["TTWDY"].SetFillColor(kRed+2);
    RsqHistosBeforeReweighting["TTZDY"].SetFillColor(kOrange-3);
    DoubleMuonRsqMC.Add(&RsqHistosBeforeReweighting["TTZDY"]);
    DoubleMuonRsqMC.Add(&RsqHistosBeforeReweighting["TTWDY"]);
    DoubleMuonRsqMC.Add(&RsqHistosBeforeReweighting["TopDY"]);
    DoubleMuonRsqMC.Add(&RsqHistosBeforeReweighting["TTJetsDY"]);
    DoubleMuonRsqMC.Add(&RsqHistosBeforeReweighting["DYJets"]);
    RsqHistosDataBeforeReweighting["DYJets"].SetMarkerStyle(20);
    RsqHistosDataBeforeReweighting["DYJets"].SetMarkerSize(1);
    DrawDataVsMCRatioPlot(&RsqHistosDataBeforeReweighting["DYJets"], &DoubleMuonRsqMC, DoubleMuonLegend, "Rsq", "controlSampleRsqBackgroundDoubleMuon", false);
    //Gamma+Jets
    THStack PhotonRsqMC("PhotonMC", "Rsq in photon control sample");
    RsqHistosBeforeReweighting["GJets"].SetFillColor(9);
    RsqHistosBeforeReweighting["QCD"].SetFillColor(8);
    RsqHistosBeforeReweighting["TTG"].SetFillColor(7);
    RsqHistosBeforeReweighting["WG"].SetFillColor(38);
    RsqHistosBeforeReweighting["ZG"].SetFillColor(4);
    PhotonRsqMC.Add(&RsqHistosBeforeReweighting["ZG"]);
    PhotonRsqMC.Add(&RsqHistosBeforeReweighting["WG"]);
    PhotonRsqMC.Add(&RsqHistosBeforeReweighting["TTG"]);
    PhotonRsqMC.Add(&RsqHistosBeforeReweighting["QCD"]);
    PhotonRsqMC.Add(&RsqHistosBeforeReweighting["GJets"]);
    RsqHistosDataBeforeReweighting["GJets"].SetMarkerStyle(20);
    RsqHistosDataBeforeReweighting["GJets"].SetMarkerSize(1);
    DrawDataVsMCRatioPlot(&RsqHistosDataBeforeReweighting["GJets"], &PhotonRsqMC, PhotonLegend, "Rsq", "controlSampleRsqBackgroundPhoton", false);

    gStyle->SetPaintTextFormat("1.2f");
    c.SetLogy(false);
    c.SetLogz(false);
    c.SetLogx();
    if(!computeDataOverMCSFs){
        //quantify agreement between DYJets and WJets predictions
        TH2F *DYWComparisonHist = (TH2F*)razorHistosData["DYJets"].Clone("DYWComparisonHist");
        TH2F *DYWSigmaComparisonHist = (TH2F*)razorHistosData["DYJets"].Clone("DYWSigmaComparisonHist");
        for(int i = 0; i < DYWComparisonHist->GetNbinsX()+1; i++){
            for(int j = 0; j < DYWComparisonHist->GetNbinsY()+1; j++){
                //set bin content to (WJets - DYJets)/DYJets
                DYWComparisonHist->SetBinContent(i, j, (razorHistosData["WJets"].GetBinContent(i, j) - razorHistosData["DYJets"].GetBinContent(i, j))/razorHistosData["DYJets"].GetBinContent(i, j));
                //set bin content to (WJets - DYJets)/(error on difference)
                float sigma1 = razorHistosData["WJets"].GetBinError(i, j);
                float sigma2 = razorHistosData["DYJets"].GetBinError(i, j);
                DYWSigmaComparisonHist->SetBinContent(i, j, (razorHistosData["WJets"].GetBinContent(i, j) - razorHistosData["DYJets"].GetBinContent(i, j))/sqrt(sigma1*sigma1+sigma2*sigma2));
            }
        }
        DYWComparisonHist->SetTitle("(WJets Prediction - DYJets Prediction)/DYJets Prediction");
        DYWComparisonHist->GetXaxis()->SetTitle("MR");
        DYWComparisonHist->GetYaxis()->SetTitle("Rsq");
        DYWComparisonHist->SetStats(0);
        DYWComparisonHist->SetMinimum(-1.0);
        DYWComparisonHist->SetMaximum(1.0);
        DYWComparisonHist->Draw("colz");
        DYWComparisonHist->Draw("same,text");
        c.Print("controlSampleHistogramComparisonDYW_noISR.pdf");
        // c.Print("controlSampleHistogramComparisonDYW.root");
        DYWSigmaComparisonHist->SetTitle("(WJets Prediction - DYJets Prediction)/#sigma_{W - DY}");
        DYWSigmaComparisonHist->GetXaxis()->SetTitle("MR");
        DYWSigmaComparisonHist->GetYaxis()->SetTitle("Rsq");
        DYWSigmaComparisonHist->SetStats(0);
        DYWSigmaComparisonHist->SetMinimum(-3);
        DYWSigmaComparisonHist->SetMaximum(3);
        DYWSigmaComparisonHist->Draw("colz");
        DYWSigmaComparisonHist->Draw("same,text");
        c.Print("controlSampleSigmaHistogramComparisonDYW_noISR.pdf");
        // c.Print("controlSampleSigmaHistogramComparisonDYW.root");

        //do the same for DYJets vs GJets
        TH2F *DYGComparisonHist = (TH2F*)razorHistosData["DYJets"].Clone("DYGComparisonHist");
        TH2F *DYGSigmaComparisonHist = (TH2F*)razorHistosData["DYJets"].Clone("DYGSigmaComparisonHist");
        for(int i = 0; i < DYGComparisonHist->GetNbinsX()+1; i++){
            for(int j = 0; j < DYGComparisonHist->GetNbinsY()+1; j++){
                //set bin content to (GJets - DYJets)/DYJets
                DYGComparisonHist->SetBinContent(i, j, (razorHistosData["GJets"].GetBinContent(i, j) - razorHistosData["DYJets"].GetBinContent(i, j))/razorHistosData["DYJets"].GetBinContent(i, j));
                //set bin content to (GJets - DYJets)/(error on difference)
                float sigma1 = razorHistosData["GJets"].GetBinError(i, j);
                float sigma2 = razorHistosData["DYJets"].GetBinError(i, j);
                DYGSigmaComparisonHist->SetBinContent(i, j, (razorHistosData["GJets"].GetBinContent(i, j) - razorHistosData["DYJets"].GetBinContent(i, j))/sqrt(sigma1*sigma1+sigma2*sigma2));
            }
        }
        DYGComparisonHist->SetTitle("(GJets Prediction - DYJets Prediction)/DYJets Prediction");
        DYGComparisonHist->GetXaxis()->SetTitle("MR");
        DYGComparisonHist->GetYaxis()->SetTitle("Rsq");
        DYGComparisonHist->SetStats(0);
        DYGComparisonHist->SetMinimum(-1.0);
        DYGComparisonHist->SetMaximum(1.0);
        DYGComparisonHist->Draw("colz");
        DYGComparisonHist->Draw("same,text");
        c.Print("controlSampleHistogramComparisonDYG_noISR.pdf");
        // c.Print("controlSampleHistogramComparisonDYG.root");
        DYGSigmaComparisonHist->SetTitle("(GJets Prediction - DYJets Prediction)/#sigma_{G - DY}");
        DYGSigmaComparisonHist->GetXaxis()->SetTitle("MR");
        DYGSigmaComparisonHist->GetYaxis()->SetTitle("Rsq");
        DYGSigmaComparisonHist->SetStats(0);
        DYGSigmaComparisonHist->SetMinimum(-3);
        DYGSigmaComparisonHist->SetMaximum(3);
        DYGSigmaComparisonHist->Draw("colz");
        DYGSigmaComparisonHist->Draw("same,text");
        c.Print("controlSampleSigmaHistogramComparisonDYG_noISR.pdf");
        // c.Print("controlSampleSigmaHistogramComparisonDYG.root");

        //and for WJets vs GJets
        TH2F *WGComparisonHist = (TH2F*)razorHistosData["WJets"].Clone("WGComparisonHist");
        TH2F *WGSigmaComparisonHist = (TH2F*)razorHistosData["WJets"].Clone("WGSigmaComparisonHist");
        for(int i = 0; i < WGComparisonHist->GetNbinsX()+1; i++){
            for(int j = 0; j < WGComparisonHist->GetNbinsY()+1; j++){
                //set bin content to (GJets - WJets)/WJets
                WGComparisonHist->SetBinContent(i, j, (razorHistosData["GJets"].GetBinContent(i, j) - razorHistosData["WJets"].GetBinContent(i, j))/razorHistosData["WJets"].GetBinContent(i, j));
                //set bin content to (GJets - WJets)/(error on difference)
                float sigma1 = razorHistosData["GJets"].GetBinError(i, j);
                float sigma2 = razorHistosData["WJets"].GetBinError(i, j);
                WGSigmaComparisonHist->SetBinContent(i, j, (razorHistosData["GJets"].GetBinContent(i, j) - razorHistosData["WJets"].GetBinContent(i, j))/sqrt(sigma1*sigma1+sigma2*sigma2));
            }
        }
        WGComparisonHist->SetTitle("(GJets Prediction - WJets Prediction)/WJets Prediction");
        WGComparisonHist->GetXaxis()->SetTitle("MR");
        WGComparisonHist->GetYaxis()->SetTitle("Rsq");
        WGComparisonHist->SetStats(0);
        WGComparisonHist->SetMinimum(-1.0);
        WGComparisonHist->SetMaximum(1.0);
        WGComparisonHist->Draw("colz");
        WGComparisonHist->Draw("same,text");
        c.Print("controlSampleHistogramComparisonWG_noISR.pdf");
        // c.Print("controlSampleHistogramComparisonWG.root");
        WGSigmaComparisonHist->SetTitle("(GJets Prediction - WJets Prediction)/#sigma_{G - W}");
        WGSigmaComparisonHist->GetXaxis()->SetTitle("MR");
        WGSigmaComparisonHist->GetYaxis()->SetTitle("Rsq");
        WGSigmaComparisonHist->SetStats(0);
        WGSigmaComparisonHist->SetMinimum(-3);
        WGSigmaComparisonHist->SetMaximum(3);
        WGSigmaComparisonHist->Draw("colz");
        WGSigmaComparisonHist->Draw("same,text");
        c.Print("controlSampleSigmaHistogramComparisonWG_noISR.pdf");
        // c.Print("controlSampleSigmaHistogramComparisonWG.root");
    }

    //plot the photon pt distribution in data and MC
    mcPhotonPt.SetStats(0);
    mcPhotonPt.SetLineColor(kViolet);
    dataPhotonPt.SetStats(0);
    mcPhotonPt.Draw();
    c.SetLogy(false);
    dataPhotonPt.Draw("pesame");
    c.Print("controlSamplePhotonPt_noISR.pdf");
    // c.Print("controlSamplePhotonPt.root");
    c.SetLogy();
    c.Print("controlSamplePhotonPtLog_noISR.pdf");
    // c.Print("controlSamplePhotonPtLog.root");

    //get the data/MC scale factors in each bin of MR and Rsq
    map<string, TH2F> dataOverMCScaleFactors;
    for(auto &tree : datatrees){
        dataOverMCScaleFactors[tree.first] = *((TH2F*)razorHistosDataBeforeReweighting[tree.first].Clone(Form("%sScaleFactors", tree.first.c_str())));
        dataOverMCScaleFactors[tree.first].Divide(&razorHistosForReweighting[tree.first]);
    }

    //write out the data/MC scale factors
    c.SetLogx(true);
    c.SetLogy(true);
    c.SetLogz(false);
    TFile sfFile("ZInvisibleScaleFactorsRun1.root", "RECREATE");
    sfFile.cd();
    for(auto &hist : dataOverMCScaleFactors){
        //plot it
        hist.second.SetStats(0);
        hist.second.SetTitle(Form("Data/MC scale factors, %s", hist.first.c_str()));
        hist.second.SetMinimum(0.0);
        hist.second.SetMaximum(2.0);
        hist.second.Draw("colz");
        hist.second.Draw("same,text");
        c.Print(Form("DataOverMC%s_noISR.pdf", hist.first.c_str()));
        // c.Print(Form("DataOverMC%s.root", hist.first.c_str()));
        //write it 
        hist.second.Write(); 
    }

    //plot the difference between the WJets and DYJets scale factors, in #sigmas and in %
    TH2F *DYWScaleHist = (TH2F*)razorHistosDataBeforeReweighting["DYJets"].Clone("DYWScaleHist");
    TH2F *DYWScaleHistPerc = (TH2F*)razorHistosDataBeforeReweighting["DYJets"].Clone("DYWScaleHistPerc");
    for(int i = 0; i < DYWScaleHist->GetNbinsX()+1; i++){
        for(int j = 0; j < DYWScaleHist->GetNbinsY()+1; j++){
            //set bin content to (WJets - DYJets)/(error on difference)
            float sigma1 = dataOverMCScaleFactors["WJets"].GetBinError(i, j);
            float sigma2 = dataOverMCScaleFactors["DYJets"].GetBinError(i, j);
            DYWScaleHist->SetBinContent(i, j, (dataOverMCScaleFactors["WJets"].GetBinContent(i, j) - dataOverMCScaleFactors["DYJets"].GetBinContent(i, j))/sqrt(sigma1*sigma1+sigma2*sigma2));
            //set bin content to (WJets - DYJets)/DYJets
            DYWScaleHistPerc->SetBinContent(i, j, (dataOverMCScaleFactors["WJets"].GetBinContent(i, j) - dataOverMCScaleFactors["DYJets"].GetBinContent(i, j))/dataOverMCScaleFactors["DYJets"].GetBinContent(i,j));
        }
    }
    DYWScaleHist->SetTitle("(WJets SF - DYJets SF)/#sigma_{W-DY}");
    DYWScaleHist->GetXaxis()->SetTitle("MR");
    DYWScaleHist->GetYaxis()->SetTitle("Rsq");
    DYWScaleHist->SetStats(0);
    DYWScaleHist->SetMinimum(-3.0);
    DYWScaleHist->SetMaximum(3.0);
    DYWScaleHist->Draw("colz");
    DYWScaleHist->Draw("same,text");
    c.Print("ScaleDYWnSigma_noISR.pdf");
    // c.Print("ScaleDYWnSigma.root");
    DYWScaleHistPerc->SetTitle("(WJets SF - DYJets SF)/DYJets SF");
    DYWScaleHistPerc->GetXaxis()->SetTitle("MR");
    DYWScaleHistPerc->GetYaxis()->SetTitle("Rsq");
    DYWScaleHistPerc->SetStats(0);
    DYWScaleHistPerc->SetMinimum(-1.0);
    DYWScaleHistPerc->SetMaximum(1.0);
    DYWScaleHistPerc->Draw("colz");
    DYWScaleHistPerc->Draw("same,text");
    c.Print("ScaleDYWPercent_noISR.pdf");
    // c.Print("ScaleDYWPercent.root");

    //do the same for WJets and GJets
    TH2F *WGScaleHist = (TH2F*)razorHistosDataBeforeReweighting["WJets"].Clone("WGScaleHist");
    TH2F *WGScaleHistPerc = (TH2F*)razorHistosDataBeforeReweighting["WJets"].Clone("WGScaleHistPerc");
    for(int i = 0; i < WGScaleHist->GetNbinsX()+1; i++){
        for(int j = 0; j < WGScaleHist->GetNbinsY()+1; j++){
            //set bin content to (GJets - WJets)/(error on difference)
            float sigma1 = dataOverMCScaleFactors["WJets"].GetBinError(i, j);
            float sigma2 = dataOverMCScaleFactors["GJets"].GetBinError(i, j);
            WGScaleHist->SetBinContent(i, j, (dataOverMCScaleFactors["GJets"].GetBinContent(i, j) - dataOverMCScaleFactors["WJets"].GetBinContent(i, j))/sqrt(sigma1*sigma1+sigma2*sigma2));
            //set bin content to (GJets - WJets)/WJets
            WGScaleHistPerc->SetBinContent(i, j, (dataOverMCScaleFactors["GJets"].GetBinContent(i, j) - dataOverMCScaleFactors["WJets"].GetBinContent(i, j))/dataOverMCScaleFactors["WJets"].GetBinContent(i,j));
        }
    }
    WGScaleHist->SetTitle("(GJets SF - WJets SF)/#sigma_{G-W}");
    WGScaleHist->GetXaxis()->SetTitle("MR");
    WGScaleHist->GetYaxis()->SetTitle("Rsq");
    WGScaleHist->SetStats(0);
    WGScaleHist->SetMinimum(-3.0);
    WGScaleHist->SetMaximum(3.0);
    WGScaleHist->Draw("colz");
    WGScaleHist->Draw("same,text");
    c.Print("ScaleWGnSigma_noISR.pdf");
    // c.Print("ScaleWGnSigma.root");
    WGScaleHistPerc->SetTitle("(GJets SF - WJets SF)/WJets SF");
    WGScaleHistPerc->GetXaxis()->SetTitle("MR");
    WGScaleHistPerc->GetYaxis()->SetTitle("Rsq");
    WGScaleHistPerc->SetStats(0);
    WGScaleHistPerc->SetMinimum(-1.0);
    WGScaleHistPerc->SetMaximum(1.0);
    WGScaleHistPerc->Draw("colz");
    WGScaleHistPerc->Draw("same,text");
    c.Print("ScaleWGPercent_noISR.pdf");
    // c.Print("ScaleWGPercent.root");

    //do the same for DYJets and GJets
    TH2F *DYGScaleHist = (TH2F*)razorHistosDataBeforeReweighting["DYJets"].Clone("DYGScaleHist");
    TH2F *DYGScaleHistPerc = (TH2F*)razorHistosDataBeforeReweighting["DYJets"].Clone("DYGScaleHistPerc");
    for(int i = 0; i < DYGScaleHist->GetNbinsX()+1; i++){
        for(int j = 0; j < DYGScaleHist->GetNbinsY()+1; j++){
            //set bin content to (GJets - DYJets)/(error on difference)
            float sigma1 = dataOverMCScaleFactors["GJets"].GetBinError(i, j);
            float sigma2 = dataOverMCScaleFactors["DYJets"].GetBinError(i, j);
            DYGScaleHist->SetBinContent(i, j, (dataOverMCScaleFactors["GJets"].GetBinContent(i, j) - dataOverMCScaleFactors["DYJets"].GetBinContent(i, j))/sqrt(sigma1*sigma1+sigma2*sigma2));
            //set bin content to (GJets - DYJets)/DYJets
            DYGScaleHistPerc->SetBinContent(i, j, (dataOverMCScaleFactors["GJets"].GetBinContent(i, j) - dataOverMCScaleFactors["DYJets"].GetBinContent(i, j))/dataOverMCScaleFactors["DYJets"].GetBinContent(i, j));
        }
    }
    DYGScaleHist->SetTitle("(GJets SF - DYJets SF)/#sigma_{G-DY}");
    DYGScaleHist->GetXaxis()->SetTitle("MR");
    DYGScaleHist->GetYaxis()->SetTitle("Rsq");
    DYGScaleHist->SetStats(0);
    DYGScaleHist->SetMinimum(-3.0);
    DYGScaleHist->SetMaximum(3.0);
    DYGScaleHist->Draw("colz");
    DYGScaleHist->Draw("same,text");
    c.Print("ScaleDYGnSigma_noISR.pdf");
    // c.Print("ScaleDYGnSigma.root");
    DYGScaleHistPerc->SetTitle("(GJets SF - DYJets SF)/DYJets SF");
    DYGScaleHistPerc->GetXaxis()->SetTitle("MR");
    DYGScaleHistPerc->GetYaxis()->SetTitle("Rsq");
    DYGScaleHistPerc->SetStats(0);
    DYGScaleHistPerc->SetMinimum(-1.0);
    DYGScaleHistPerc->SetMaximum(1.0);
    DYGScaleHistPerc->Draw("colz");
    DYGScaleHistPerc->Draw("same,text");
    c.Print("ScaleDYGPercent_noISR.pdf");
    // c.Print("ScaleDYGPercent.root");
}

int main(){
    ZInvisibleControlSamples();
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
