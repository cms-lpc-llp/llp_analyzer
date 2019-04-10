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
#include "TLatex.h"
#include "assert.h"
#include "math.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "MacroHelper.h"

using namespace std;

bool debug = false;
// bool debug = true;

//lepton Pt, Yields, NJets 40, Electrons

//true: find the translation factors from MC to data
//false: find the translation factors from DY, W, G to Z->nunu
bool computeDataOverMCSFs = false;

void ZInvisibleCrossChecks_WJetsRun2_muele(){
    gROOT->SetBatch();

    float maxMuonPt = 999;

    //set color palette 
    const Int_t NCont = 100;
    gStyle->SetNumberContours(NCont);
    gStyle->SetPaintTextFormat("1.0f");

    const int lumi = 20.3844; 
    
    map<string, string> suffixes;
    suffixes["WJets"] = "";
    suffixes["TTJets"] = "";
    suffixes["Top"] = "";
    // suffixes["TTW"] = "_noW";
    // suffixes["TTZ"] = "_noW";

    // // load file with fit results
    // TFile *fitResultFile = new TFile("./ControlSampleFits_v3.root", "READ");
    // vector<TH1F*> fitHists;
    // if(fitResultFile){
    //     fitHists.push_back((TH1F*)fitResultFile->Get("ControlSampleFits/WSingleLepton/Sideband/h_MR"));
    //     fitHists.push_back((TH1F*)fitResultFile->Get("ControlSampleFits/TTBarSingleLepton/Sideband/h_MR"));
    // }

    //get input files -- assumes one TFile for each process, with weights for different HT bins 
    map<string, TFile*> mcfiles;
    map<string, TFile*> datafiles;

    // Reduced
    mcfiles["TTJets"] = new TFile("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p17/OneLeptonReduced/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
    mcfiles["Top"] = new TFile("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p17/OneLeptonReduced/SingleTop_1pb_weighted.root");
    mcfiles["WJets"] = new TFile("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p17/OneLeptonReduced/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root");

    datafiles["WJets"] = new TFile("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p17/OneLeptonReduced_v2/SingleMuon_Run2015C-GOLDEN.root");
    //datafiles["WJets"] = new TFile("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p17/OneLeptonReduced/SingleElectron_Run2015C-GOLDEN.root");

    //get trees and set branches
    map<string, TTree*> mctrees;
    map<string, TTree*> datatrees;
    map<string, Float_t> mets;
    map<string, Float_t> mrs;
    map<string, Float_t> rsqs;
    map<string, UInt_t> njets40;
    map<string, UInt_t> njets80;
    map<string, Float_t> hts;
    map<string, Float_t> mhts;
    map<string, Float_t> mhtnohfs;
    
    float weight;
    int nPU_mean, nTightMuons, nLooseMuons;
    float leadingMuonPt, leadingMuonEta, recoZpt, recoZeta, recoZmass, subleadingMuonPt, subleadingMuonEta, mTLepMet, mTLepMetNoHF;
    bool Flag_HBHENoiseFilter, Flag_CSCTightHaloFilter, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_eeBadScFilter, Flag_ecalLaserCorrFilter;
    float bjet1Pt, bjet2Pt;
    bool bjet1PassMedium, bjet2PassMedium;
    float genZPt, genWPt;
    bool hlt_singlemu, hlt_dimuon;
    float lep1Pt, lep1Eta;
    float MR, Rsq;
    UInt_t nVtx, nBTaggedJets;
    Bool_t HLTNames[156];
    Bool_t lep1passTight;
    Int_t lep1Type;

    for(auto &file : mcfiles){
        mets[file.first] = 0.;
        mrs[file.first] = 0.;
        rsqs[file.first] = 0.;
        njets40[file.first] = 0.;
        njets80[file.first] = 0.;
        hts[file.first] = 0.;
        mhts[file.first] = 0.;
        mhtnohfs[file.first] = 0.;

       mctrees[file.first] = (TTree*)file.second->Get("ControlSampleEvent");

        mctrees[file.first]->SetBranchStatus("*", 0); // disable all
        mctrees[file.first]->SetBranchStatus("weight", 1);
        mctrees[file.first]->SetBranchStatus("NPV", 1); // enable 
        mctrees[file.first]->SetBranchStatus("NJets40", 1);
        mctrees[file.first]->SetBranchStatus("NJets80", 1);
        mctrees[file.first]->SetBranchStatus("MR", 1);
        mctrees[file.first]->SetBranchStatus("Rsq", 1);
        mctrees[file.first]->SetBranchStatus("MET", 1);
        mctrees[file.first]->SetBranchStatus("lep1MT", 1);
        mctrees[file.first]->SetBranchStatus("lep1MTnoHF", 1);
        mctrees[file.first]->SetBranchStatus("NBJetsMedium", 1);
        mctrees[file.first]->SetBranchStatus("lep1PassTight", 1);
        mctrees[file.first]->SetBranchStatus("HLTDecision", 1);
        mctrees[file.first]->SetBranchStatus("HT", 1);
        mctrees[file.first]->SetBranchStatus("METnoHF", 1);
        mctrees[file.first]->SetBranchStatus("lep1Pt", 1);
        mctrees[file.first]->SetBranchStatus("lep1Eta", 1);
        mctrees[file.first]->SetBranchStatus("lep1Type", 1);
	mctrees[file.first]->SetBranchStatus("HLTDecision", 1);

	mctrees[file.first]->SetBranchAddress("HLTDecision", HLTNames);

        mctrees[file.first]->SetBranchAddress(Form("MET%s", suffixes[file.first].c_str()), &mets[file.first]);
	mctrees[file.first]->SetBranchAddress(Form("MR%s", suffixes[file.first].c_str()), &mrs[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("Rsq%s", suffixes[file.first].c_str()), &rsqs[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("NJets40%s", suffixes[file.first].c_str()), &njets40[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("NJets80%s", suffixes[file.first].c_str()), &njets80[file.first]);
	mctrees[file.first]->SetBranchAddress(Form("HT%s", suffixes[file.first].c_str()), &hts[file.first]);
	mctrees[file.first]->SetBranchAddress(Form("METnoHF%s", suffixes[file.first].c_str()), &mhtnohfs[file.first]);
	       
        mctrees[file.first]->SetBranchAddress("NBJetsMedium", &nBTaggedJets);
        mctrees[file.first]->SetBranchAddress("lep1MT", &mTLepMet);
        mctrees[file.first]->SetBranchAddress("lep1MTnoHF", &mTLepMetNoHF);
        mctrees[file.first]->SetBranchAddress("lep1PassTight", &lep1passTight);
        mctrees[file.first]->SetBranchAddress("lep1Pt", &lep1Pt);
        mctrees[file.first]->SetBranchAddress("lep1Eta", &lep1Eta);
        mctrees[file.first]->SetBranchAddress("lep1Type", &lep1Type);

        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("NPV", &nVtx); // enable 
    }
    for(auto &file : datafiles){
        datatrees[file.first] = (TTree*)file.second->Get("ControlSampleEvent");
        
        datatrees[file.first]->SetBranchStatus("*", 0); // disable all
        datatrees[file.first]->SetBranchStatus("NPV", 1); // enable 
        datatrees[file.first]->SetBranchStatus("Flag_HBHENoiseFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_CSCTightHaloFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_eeBadScFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_ecalLaserCorrFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("NBJetsMedium", 1);
        datatrees[file.first]->SetBranchStatus("MET", 1);
        datatrees[file.first]->SetBranchStatus("MR", 1);
        datatrees[file.first]->SetBranchStatus("Rsq", 1);
        datatrees[file.first]->SetBranchStatus("NJets40", 1);
        datatrees[file.first]->SetBranchStatus("NJets80", 1);
        datatrees[file.first]->SetBranchStatus("lep1MT", 1);
        datatrees[file.first]->SetBranchStatus("lep1MTnoHF", 1);
        datatrees[file.first]->SetBranchStatus("lep1PassTight", 1);
        datatrees[file.first]->SetBranchStatus("lep1Pt", 1);
        datatrees[file.first]->SetBranchStatus("HLTDecision", 1);
        datatrees[file.first]->SetBranchStatus("HT", 1);
        datatrees[file.first]->SetBranchStatus("METnoHF", 1);
        datatrees[file.first]->SetBranchStatus("lep1Pt", 1);
        datatrees[file.first]->SetBranchStatus("lep1Eta", 1);
        datatrees[file.first]->SetBranchStatus("lep1Type", 1);
	datatrees[file.first]->SetBranchStatus("HLTDecision", 1);

	datatrees[file.first]->SetBranchAddress("HLTDecision", HLTNames);
        datatrees[file.first]->SetBranchAddress("lep1MT", &mTLepMet);
        datatrees[file.first]->SetBranchAddress("lep1MTnoHF", &mTLepMetNoHF);
        datatrees[file.first]->SetBranchAddress(Form("MET%s", suffixes[file.first].c_str()), &mets[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("MR%s", suffixes[file.first].c_str()), &mrs[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("Rsq%s", suffixes[file.first].c_str()), &rsqs[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("NJets40%s", suffixes[file.first].c_str()), &njets40[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("NJets80%s", suffixes[file.first].c_str()), &njets80[file.first]);
	datatrees[file.first]->SetBranchAddress(Form("HT%s", suffixes[file.first].c_str()), &hts[file.first]);
	datatrees[file.first]->SetBranchAddress(Form("METnoHF%s", suffixes[file.first].c_str()), &mhtnohfs[file.first]);
        datatrees[file.first]->SetBranchAddress("NPV", &nVtx); // enable 
        datatrees[file.first]->SetBranchAddress("lep1Pt", &lep1Pt); // enable 
        datatrees[file.first]->SetBranchAddress("lep1Eta", &lep1Eta); // enable 
        datatrees[file.first]->SetBranchAddress("lep1Type", &lep1Type);

        datatrees[file.first]->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter); // enable 
        datatrees[file.first]->SetBranchAddress("lep1PassTight", &lep1passTight);
        datatrees[file.first]->SetBranchAddress("NBJetsMedium", &nBTaggedJets);
    }

    //load efficiency/acceptance histograms
    //TFile *effFile = new TFile("./Run1LeptonPhotonEfficiency.root", "READ");
    //TH2F muonLooseEffHisto = *(TH2F *)effFile->Get("MuonEfficiency");
    //TH2F muonTightEffHisto = *(TH2F *)effFile->Get("MuonEfficiencyTight");
    //TH2F zAccHisto = *(TH2F *)effFile->Get("MuonAcceptance");

    //load muon efficiency scale factor histogram
    TFile muIdSFFile("data/ScaleFactors/MuonEfficiencies_ID_Run2012ReReco_53X.root");
    TFile muIsoSFFile("data/ScaleFactors/MuonEfficiencies_ISO_Run_2012ReReco_53X.root");
    //TODO: add muon efficiency scale factors
    float singleMuTriggerSF = 0.97;
    float doubleMuTriggerSF = 0.97;
    float doubleMuNormalizationSF = 0.97;

    //load pileup reweighting histogram
    // TFile *pileupWeightFile = new TFile("data/Run1PileupWeights.root", "READ");
    // TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("PUWeight_Run1");
    TFile *pileupWeightFile = new TFile("./NVtx_Run2015C_SingleMuon.root", "READ");
    TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("dataNvtx0");
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

    cutSequence.push_back( "NBJetsMedium == 0 && lep1PassTight == 1 && lep1MTnoHF > 30 && lep1MTnoHF < 100 && lep1Pt > 25 && TMath::Abs(lep1Eta)<2.5 && METnoHF > 50" );
    cutName.push_back( "" );

    cutSequence.push_back( "NBJetsMedium > 0 && lep1PassTight == 1 && lep1MTnoHF > 30 && lep1MTnoHF < 100 && lep1Pt > 25 && TMath::Abs(lep1Eta)<2.5 && METnoHF > 30" );
    cutName.push_back( "" );

    cutSequence.push_back( "lep1PassTight == 1 && lep1Pt > 25 && TMath::Abs(lep1Eta)<2.5" );
    cutName.push_back( "" );

    map<string, vector<TH1F *> > mcNJets40, mcNJets80, mcMR, mcRsq,  mcMet, mcNvtx,  mcHT, mcMT, mcMTnoHF, mcLep1Pt, mcLep1Eta, mcMETnoHF;
    vector<TH1F *>  dataNJets40, dataNJets80, dataMR, dataRsq, dataMet, dataNvtx, dataHT, dataMT, dataMTnoHF, dataLep1Pt, dataLep1Eta, dataMETnoHF;
    for(auto &tree : mctrees){
        mcNvtx[tree.first] = vector<TH1F *>();
        mcNJets40[tree.first] = vector<TH1F *>();
        mcNJets80[tree.first] = vector<TH1F *>();
        mcMR[tree.first] = vector<TH1F *>();
        mcRsq[tree.first] = vector<TH1F *>();
        mcMet[tree.first] = vector<TH1F *>();
        mcHT[tree.first] = vector<TH1F *>();
        mcMT[tree.first] = vector<TH1F *>();
        mcMTnoHF[tree.first] = vector<TH1F *>();
        mcLep1Pt[tree.first] = vector<TH1F *>();
        mcLep1Eta[tree.first] = vector<TH1F *>();
        mcMETnoHF[tree.first] = vector<TH1F *>();
    }
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        for(auto &tree : mctrees){
            mcNvtx[tree.first].push_back(new TH1F(Form("mcNvtx%s%d", tree.first.c_str(), cut), Form("%s; NVtx (GeV)", cutName[cut].c_str()), 50, 0, 50));
            mcNJets40[tree.first].push_back(new TH1F(Form("mcNJets40%s%d", tree.first.c_str(), cut), Form("%s; Number of jets 40 GeV", cutName[cut].c_str()), 10, 0, 10));
            mcNJets80[tree.first].push_back(new TH1F(Form("mcNJets80%s%d", tree.first.c_str(), cut), Form("%s; Number of jets 80 GeV", cutName[cut].c_str()), 10, 0, 10));
            mcMR[tree.first].push_back(new TH1F(Form("mcMR%s%d", tree.first.c_str(), cut), Form("%s; MR [GeV/c^2]", cutName[cut].c_str()), 20, 300, 4300));
            mcRsq[tree.first].push_back(new TH1F(Form("mcRsq%s%d", tree.first.c_str(), cut), Form("%s; Rsq (GeV)", cutName[cut].c_str()), nRsqBins, RsqBinLowEdges));
            mcMet[tree.first].push_back(new TH1F(Form("mcMet%s%d", tree.first.c_str(), cut), Form("%s; MET (GeV)", cutName[cut].c_str()), 200, 0, 1000));
            mcMT[tree.first].push_back(new TH1F(Form("mcMT%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTnoHF[tree.first].push_back(new TH1F(Form("mcMTnoHF%s%d", tree.first.c_str(), cut), Form("%s; MT noHF (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcHT[tree.first].push_back(new TH1F(Form("mcHT%s%d", tree.first.c_str(), cut), Form("%s; HT (GeV)", cutName[cut].c_str()), 50, 0, 1000));
            mcMETnoHF[tree.first].push_back(new TH1F(Form("mcMETnoHF%s%d", tree.first.c_str(), cut), Form("%s; MET noHF (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcLep1Pt[tree.first].push_back(new TH1F(Form("mcLep1Pt%s%d", tree.first.c_str(), cut), Form("%s; Lepton Pt (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcLep1Eta[tree.first].push_back(new TH1F(Form("mcLep1Eta%s%d", tree.first.c_str(), cut), Form("%s; Lepton Eta", cutName[cut].c_str()), 50, -4, 4));

            mcNJets40[tree.first][cut]->Sumw2();
            mcNJets80[tree.first][cut]->Sumw2();
            mcMR[tree.first][cut]->Sumw2();
            mcRsq[tree.first][cut]->Sumw2();
            mcMet[tree.first][cut]->Sumw2();
            mcNvtx[tree.first][cut]->Sumw2();
	    mcHT[tree.first][cut]->Sumw2();	    
	    mcMT[tree.first][cut]->Sumw2();	    
	    mcMTnoHF[tree.first][cut]->Sumw2();	    
	    mcLep1Pt[tree.first][cut]->Sumw2();	    
	    mcLep1Eta[tree.first][cut]->Sumw2();	    
	    mcMETnoHF[tree.first][cut]->Sumw2();	    
        }
        dataNJets40.push_back(new TH1F(Form("dataNJets40%d", cut), Form("%s; Number of jets 40 GeV", cutName[cut].c_str()), 10, 0, 10));
        dataNJets80.push_back(new TH1F(Form("dataNJets80%d", cut), Form("%s; Number of jets 80 GeV", cutName[cut].c_str()), 10, 0, 10));
        dataMR.push_back(new TH1F(Form("dataMR%d", cut), Form("%s; MR [GeV/c^2]", cutName[cut].c_str()), 20, 300, 4300));
        dataRsq.push_back(new TH1F(Form("dataRsq%d", cut), Form("%s; Rsq (GeV)", cutName[cut].c_str()), nRsqBins, RsqBinLowEdges));
        dataMet.push_back(new TH1F(Form("dataMet%d", cut), Form("%s; mcMet (GeV)", cutName[cut].c_str()), 200, 0, 1000));
        dataNvtx.push_back(new TH1F(Form("dataNvtx%d", cut), Form("%s; NVtx (GeV)", cutName[cut].c_str()), 50, 0, 50));
        dataHT.push_back(new TH1F(Form("dataHT%d", cut), Form("%s; HT (GeV)", cutName[cut].c_str()), 50, 0, 1000));
        dataMT.push_back(new TH1F(Form("dataMT%d", cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
        dataMTnoHF.push_back(new TH1F(Form("dataMTnoHF%d", cut), Form("%s; MT noHF (GeV)", cutName[cut].c_str()), 50, 0, 500));
        dataLep1Pt.push_back(new TH1F(Form("datalep1Pt%d", cut), Form("%s; Lepton Pt (GeV)", cutName[cut].c_str()), 50, 0, 500));
        dataLep1Eta.push_back(new TH1F(Form("datalep1Eta%d", cut), Form("%s; Lepton Eta", cutName[cut].c_str()), 50, -4, 4));
        dataMETnoHF.push_back(new TH1F(Form("dataMETnoHF%d", cut), Form("%s; MET noHF (GeV)", cutName[cut].c_str()), 50, 0, 500));

        dataNJets40[cut]->Sumw2();
        dataNJets80[cut]->Sumw2();
        dataMR[cut]->Sumw2();
        dataRsq[cut]->Sumw2();
        dataMet[cut]->Sumw2();
	dataNvtx[cut]->Sumw2();
	dataHT[cut]->Sumw2();
	dataMT[cut]->Sumw2();
	dataMTnoHF[cut]->Sumw2();
	dataLep1Pt[cut]->Sumw2();
	dataMETnoHF[cut]->Sumw2();
    }
    for(auto &tree : mctrees){
        cout << "Filling MC histograms: " << tree.first << endl;
        uint nEntries = tree.second->GetEntries();
        if(debug) nEntries = 1000;
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

	    if(i % 100000 == 0) cout << "Processing entry " << i << " of "<<tree.first<<endl;
            //get event weight
            float eventWeight = weight;
            //eventWeight *= pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(nVtx));

	    bool trigger_passed = false;
	    if (HLTNames[8] && TMath::Abs(lep1Type) == 13 && TMath::Abs(lep1Eta)<2.4) trigger_passed = true; //muon triggers
	    //if (TMath::Abs(lep1Type) == 11 && lep1Pt > 35) trigger_passed = true; //electron triggers

	    if (trigger_passed == false) continue;
	    
	    eventWeight *= lumi;

            //apply selection cuts and fill the appropriate histograms
            for(uint cut = 0; cut < cutSequence.size(); cut++){
                bool passesCut = cuts[cut]->EvalInstance();
                if(!passesCut) continue;


                (mcNJets40[tree.first])[cut]->Fill(njets40[tree.first], eventWeight);
                (mcNJets80[tree.first])[cut]->Fill(njets80[tree.first], eventWeight);
                (mcMR[tree.first])[cut]->Fill(mrs[tree.first], eventWeight); 
                (mcRsq[tree.first])[cut]->Fill(rsqs[tree.first], eventWeight);
                (mcMet[tree.first])[cut]->Fill(mets[tree.first], eventWeight);
                (mcNvtx[tree.first])[cut]->Fill(nVtx, eventWeight);
		(mcHT[tree.first])[cut]->Fill(hts[tree.first], eventWeight);
		(mcMT[tree.first])[cut]->Fill(mTLepMet, eventWeight);
		(mcMTnoHF[tree.first])[cut]->Fill(mTLepMetNoHF, eventWeight);
		(mcLep1Pt[tree.first])[cut]->Fill(lep1Pt, eventWeight);
		(mcLep1Eta[tree.first])[cut]->Fill(lep1Eta, eventWeight);
		(mcMETnoHF[tree.first])[cut]->Fill(mhtnohfs[tree.first], eventWeight);
            }
        }

        for(uint cut = 0; cut < cutSequence.size(); cut++){
            delete cuts[cut];
        }
    }

    for(auto &tree : datatrees){
        cout << "Filling data histograms: " << tree.first << endl;
        uint nEntries = tree.second->GetEntries();
        if(debug) nEntries = 1000;
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
	    if(i % 10000 == 0) cout << "Processing entry " << i << " of "<<tree.first<<endl;

            //get event weight
            float eventWeight = 1.0;
	    
	    bool trigger_passed = false;
	    if (HLTNames[8] && TMath::Abs(lep1Type) == 13 && TMath::Abs(lep1Eta)<2.4) trigger_passed = true; //muon triggers
	    //if (TMath::Abs(lep1Type) == 11 && lep1Pt > 35) trigger_passed = true; //electron triggers

	    if (trigger_passed == false) continue;

            //apply selection cuts and fill the appropriate histograms
            for(uint cut = 0; cut < cutSequence.size(); cut++){
                bool passesCut = cuts[cut]->EvalInstance();
                if(!passesCut) continue;

		dataNJets40[cut]->Fill(njets40[tree.first], eventWeight);
                dataNJets80[cut]->Fill(njets80[tree.first], eventWeight);
                dataMR[cut]->Fill(mrs[tree.first], eventWeight);
                dataRsq[cut]->Fill(rsqs[tree.first], eventWeight);
                dataMet[cut]->Fill(mets[tree.first], eventWeight);
		dataNvtx[cut]->Fill(nVtx, eventWeight);
		dataHT[cut]->Fill(hts[tree.first], eventWeight);
		dataMT[cut]->Fill(mTLepMet, eventWeight);
		dataMTnoHF[cut]->Fill(mTLepMetNoHF, eventWeight);
		dataLep1Pt[cut]->Fill(lep1Pt, eventWeight);
		dataLep1Eta[cut]->Fill(lep1Eta, eventWeight);
		dataMETnoHF[cut]->Fill(mhtnohfs[tree.first], eventWeight);
           }
        }
        for(uint cut = 0; cut < cutSequence.size(); cut++){
            delete cuts[cut];
        }
    }

    //print out plots
    TCanvas c("c", "c", 800, 700);
    c.SetLogy();

    //colors and legend
    map<string, int> colors;
    colors["WJets"] = kRed;
    colors["Top"] = kBlue;
    colors["TTJets"] = kAzure+10;

    TLegend *legend = new TLegend(0.6, 0.54, 0.9, 0.84);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);


    // TLegend *legend = new TLegend(0.5, 0.5, 0.9, 0.9);
    legend->AddEntry(dataNJets40[0], "Data");
    legend->AddEntry(mcNJets40["WJets"][0], "W+Jets MC");
    legend->AddEntry(mcNJets40["TTJets"][0], "TTJets MC");
    legend->AddEntry(mcNJets40["Top"][0], "Single Top MC");
    // TLegend *legendWithFit = (TLegend*)legend->Clone();
    // legend->AddEntry(fitHists[0], "Fit");

    for(uint cut = 0; cut < cutSequence.size(); cut++){
        //create histogram stacks for MC
        THStack NumJets40MC(Form("NumJets40Stack%d", cut), cutName[cut].c_str());
        THStack NumJets80MC(Form("NumJets80Stack%d", cut), cutName[cut].c_str());
        THStack MRMC(Form("MRStack%d", cut), cutName[cut].c_str());
        THStack RsqMC(Form("RsqStack%d", cut), cutName[cut].c_str());
        THStack MetMC(Form("MetStack%d", cut), cutName[cut].c_str());
        THStack NVtxMC(Form("NVtxStack%d", cut), cutName[cut].c_str());
	THStack HTMC(Form("HTStack%d", cut), cutName[cut].c_str());
	THStack MTMC(Form("MTMC%d", cut), cutName[cut].c_str());
	THStack MTnoHFMC(Form("MTnoHFMC%d", cut), cutName[cut].c_str());
	THStack LEP1PTMC(Form("LEP1PTMC%d", cut), cutName[cut].c_str());
	THStack LEP1ETAMC(Form("LEP1ETAMC%d", cut), cutName[cut].c_str());
	THStack METNOHFMC(Form("METNOHFMC%d", cut), cutName[cut].c_str());
	
        //add the histograms to the stack in order
        // vector<string> orderedtrees {"TTZ", "TTW", "Top", "TTJets", "WJets"};
        vector<string> orderedtrees {"Top", "TTJets", "WJets"};
        for(auto &tree : orderedtrees){
	    mcNJets40[tree][cut]->SetFillColor(colors[tree]);
            mcNJets80[tree][cut]->SetFillColor(colors[tree]);
            mcMR[tree][cut]->SetFillColor(colors[tree]);
            mcRsq[tree][cut]->SetFillColor(colors[tree]);
            mcMet[tree][cut]->SetFillColor(colors[tree]);
            mcNvtx[tree][cut]->SetFillColor(colors[tree]);
	    mcHT[tree][cut]->SetFillColor(colors[tree]);
	    mcMT[tree][cut]->SetFillColor(colors[tree]);
	    mcMTnoHF[tree][cut]->SetFillColor(colors[tree]);
	    mcLep1Pt[tree][cut]->SetFillColor(colors[tree]);
	    mcLep1Eta[tree][cut]->SetFillColor(colors[tree]);
	    mcMETnoHF[tree][cut]->SetFillColor(colors[tree]);
	    
	    mcNJets40[tree][cut]->SetLineColor(colors[tree]);
            mcNJets80[tree][cut]->SetLineColor(colors[tree]);
            mcMR[tree][cut]->SetLineColor(colors[tree]);
            mcRsq[tree][cut]->SetLineColor(colors[tree]);
            mcMet[tree][cut]->SetLineColor(colors[tree]);
            mcNvtx[tree][cut]->SetLineColor(colors[tree]);
	    mcHT[tree][cut]->SetLineColor(colors[tree]);
	    mcMT[tree][cut]->SetLineColor(colors[tree]);
	    mcMTnoHF[tree][cut]->SetLineColor(colors[tree]);
	    mcLep1Pt[tree][cut]->SetLineColor(colors[tree]);
	    mcLep1Eta[tree][cut]->SetLineColor(colors[tree]);
	    mcMETnoHF[tree][cut]->SetLineColor(colors[tree]);

	    mcNJets40[tree][cut]->SetMarkerColor(colors[tree]);
            mcNJets80[tree][cut]->SetMarkerColor(colors[tree]);
            mcMR[tree][cut]->SetMarkerColor(colors[tree]);
            mcRsq[tree][cut]->SetMarkerColor(colors[tree]);
            mcMet[tree][cut]->SetMarkerColor(colors[tree]);
            mcNvtx[tree][cut]->SetMarkerColor(colors[tree]);
	    mcHT[tree][cut]->SetMarkerColor(colors[tree]);
	    mcMT[tree][cut]->SetMarkerColor(colors[tree]);
	    mcMTnoHF[tree][cut]->SetMarkerColor(colors[tree]);
	    mcLep1Pt[tree][cut]->SetMarkerColor(colors[tree]);
	    mcLep1Eta[tree][cut]->SetMarkerColor(colors[tree]);
	    mcMETnoHF[tree][cut]->SetMarkerColor(colors[tree]);

            NumJets40MC.Add(mcNJets40[tree][cut]);
            NumJets80MC.Add(mcNJets80[tree][cut]);
            MRMC.Add(mcMR[tree][cut]);
            RsqMC.Add(mcRsq[tree][cut]);
	    MetMC.Add(mcMet[tree][cut]);
	    NVtxMC.Add(mcNvtx[tree][cut]);
	    HTMC.Add(mcHT[tree][cut]);
	    MTMC.Add(mcMT[tree][cut]);
	    MTnoHFMC.Add(mcMTnoHF[tree][cut]);
	    LEP1PTMC.Add(mcLep1Pt[tree][cut]);
	    LEP1ETAMC.Add(mcLep1Eta[tree][cut]);
	    METNOHFMC.Add(mcMETnoHF[tree][cut]);
        }
	DrawDataVsMCRatioPlot(dataNJets40[cut], &NumJets40MC, legend, "Number of jets 40 GeV", "ControlRegionPlots_Mu_NumJets40"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataNJets80[cut], &NumJets80MC, legend, "Number of jets 80 GeV", "ControlRegionPlots_Mu_NumJets80"+to_string(cut), false);
        // DrawDataVsMCRatioPlot(dataMR[cut], &MRMC, legend, "M_{R} [GeV]", "ControlRegionPlots_Mu_MR"+to_string(cut), false, "40 pb^{-1}", fitHists[cut]);
        DrawDataVsMCRatioPlot(dataMR[cut], &MRMC, legend, "MR (GeV)", "ControlRegionPlots_Mu_MR"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataRsq[cut], &RsqMC, legend, "Rsq", "ControlRegionPlots_Mu_Rsq"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataMet[cut], &MetMC, legend, "Met", "ControlRegionPlots_Mu_MET"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataNvtx[cut], &NVtxMC, legend, "NVtx", "ControlRegionPlots_Mu_NVtx"+to_string(cut), false);
	DrawDataVsMCRatioPlot(dataHT[cut], &HTMC, legend, "HT", "ControlRegionPlots_Mu_HT"+to_string(cut), false);
	DrawDataVsMCRatioPlot(dataMT[cut], &MTMC, legend, "MT", "ControlRegionPlots_Mu_MT"+to_string(cut), false);
	DrawDataVsMCRatioPlot(dataMTnoHF[cut], &MTnoHFMC, legend, "MTnoHF", "ControlRegionPlots_Mu_MTnoHF"+to_string(cut), false);
	DrawDataVsMCRatioPlot(dataLep1Pt[cut], &LEP1PTMC, legend, "Lep1Pt", "ControlRegionPlots_Mu_Lep1Pt"+to_string(cut), false);
	DrawDataVsMCRatioPlot(dataLep1Eta[cut], &LEP1ETAMC, legend, "Lep1Eta", "ControlRegionPlots_Mu_Lep1Eta"+to_string(cut), false);
	DrawDataVsMCRatioPlot(dataMETnoHF[cut], &METNOHFMC, legend, "METnoHF", "ControlRegionPlots_Mu_METnoHF"+to_string(cut), false);
    }

    delete legend;
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        delete dataNJets40[cut];
        delete dataNJets80[cut];
        delete dataMR[cut];
        delete dataRsq[cut];
 	delete dataMet[cut];
	delete dataNvtx[cut];
	delete dataHT[cut];
	delete dataMT[cut];
	delete dataMTnoHF[cut];
	delete dataLep1Pt[cut];
	delete dataLep1Eta[cut];
	delete dataMETnoHF[cut];
	for(auto &tree : mctrees){
	  delete mcNJets40[tree.first][cut];
	  delete mcNJets80[tree.first][cut];
	  delete mcMR[tree.first][cut];
	  delete mcRsq[tree.first][cut];
	  delete mcMet[tree.first][cut];
	  delete mcNvtx[tree.first][cut];
	  delete mcHT[tree.first][cut];
	  delete mcMT[tree.first][cut];
	  delete mcMTnoHF[tree.first][cut];
	  delete mcLep1Pt[tree.first][cut];
	  delete mcLep1Eta[tree.first][cut];
	  delete mcMETnoHF[tree.first][cut];
        }
    }

    mcfiles["WJets"]->Close();
    mcfiles["Top"]->Close();
    mcfiles["TTJets"]->Close();
    // mcfiles["TTW"]->Close();
    // mcfiles["TTZ"]->Close();
    datafiles["WJets"]->Close();
    pileupWeightFile->Close();
    //effFile->Close();
    delete mcfiles["WJets"];
    delete mcfiles["Top"];
    delete mcfiles["TTJets"];
    // delete mcfiles["TTW"];
    // delete mcfiles["TTZ"];
    delete datafiles["WJets"];
    delete pileupWeightFile;
    //delete effFile;
    cout<<"HERRE 5"<<endl;

}

int main(){
    ZInvisibleCrossChecks_WJetsRun2_muele();
    return 0;
}
