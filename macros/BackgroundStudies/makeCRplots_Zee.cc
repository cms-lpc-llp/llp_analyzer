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
#include "TLorentzVector.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "MacroHelper.h"

using namespace std;

bool debug = false;
//bool debug = true;

//lepton Pt, Yields, NJets 40, Electrons

//true: find the translation factors from MC to data
//false: find the translation factors from DY, W, G to Z->nunu
bool computeDataOverMCSFs = false;

float GetElectronScaleCorrection( double pt, double eta ) {
  double scaleCorr = 1.0;
  if ( pt > 0 && fabs(eta) < 1.5) {
    scaleCorr = 1.015;
  } else {
    scaleCorr = 1.05;
  }
  return scaleCorr;
}

void makeCRplots_GJets(){
    gROOT->SetBatch();

    float maxMuonPt = 999;

    //set color palette 
    const Int_t NCont = 100;
    gStyle->SetNumberContours(NCont);
    gStyle->SetPaintTextFormat("1.0f");

    const int lumi = 16.1; 
    
    map<string, string> suffixes;
    suffixes["ZJets"] = "_NoPho";

    //get input files -- assumes one TFile for each process, with weights for different HT bins 
    map<string, TFile*> mcfiles;
    map<string, TFile*> datafiles;

    // Reduced
    //    mcfiles["ZJets"] = new TFile("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p17/PhotonJets/DYJetsToLL_HT-100toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");
     mcfiles["ZJets"] = new TFile("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p17/PhotonJets/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root");

    datafiles["ZJets"] = new TFile("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p17/PhotonJets/DoubleEG_Run2015C-GOLDEN.root");

    //get trees and set branches
    map<string, TTree*> mctrees;
    map<string, TTree*> datatrees;
    map<string, Float_t> mets;
    map<string, Float_t> metphis;
    map<string, Float_t> mrs;
    map<string, Float_t> rsqs;
    map<string, UInt_t> njets40;
    map<string, UInt_t> njets80;
    map<string, Float_t> hts;
    map<string, Float_t> mhts;
    map<string, Float_t> mhtnohfs;
    
    float weight;
    int nPU_mean, nTightMuons, nLooseMuons;
    bool Flag_HBHENoiseFilter, Flag_CSCTightHaloFilter, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_eeBadScFilter, Flag_ecalLaserCorrFilter;
    bool bjet1PassMedium, bjet2PassMedium;
    float MR, Rsq;
    UInt_t nVtx, nBTaggedJets, nSelectedPhotons;
    TLorentzVector* pho1 = NULL; 
    TLorentzVector* pho2 = NULL; 
    TLorentzVector* jet1 = NULL; 
    TLorentzVector* jet2 = NULL; 
    Bool_t HLTDecision[156];
    Int_t pho1_motherID, pho2_motherID;
    Float_t pho1_sigmaietaieta, pho2_sigmaietaieta;

    for(auto &file : mcfiles){
        mets[file.first] = 0.;
	metphis[file.first] = 0.;
	mrs[file.first] = 0.;
        rsqs[file.first] = 0.;
        njets40[file.first] = 0.;
        njets80[file.first] = 0.;
        hts[file.first] = 0.;
        mhts[file.first] = 0.;
        mhtnohfs[file.first] = 0.;

	mctrees[file.first] = (TTree*)file.second->Get("ControlSampleEvent");
       
	mctrees[file.first]->SetBranchAddress("HLTDecision", HLTDecision);

        mctrees[file.first]->SetBranchAddress(Form("MET%s", suffixes[file.first].c_str()), &mets[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("METPhi%s", suffixes[file.first].c_str()), &metphis[file.first]);
	mctrees[file.first]->SetBranchAddress(Form("MR%s", suffixes[file.first].c_str()), &mrs[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("Rsq%s", suffixes[file.first].c_str()), &rsqs[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("NJets%s", suffixes[file.first].c_str()), &njets40[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("NJets80%s", suffixes[file.first].c_str()), &njets80[file.first]);
	mctrees[file.first]->SetBranchAddress(Form("HT%s", suffixes[file.first].c_str()), &hts[file.first]);
	       
        mctrees[file.first]->SetBranchAddress("nSelectedPhotons", &nSelectedPhotons);
        mctrees[file.first]->SetBranchAddress("pho1", &pho1);
        mctrees[file.first]->SetBranchAddress("pho2", &pho2);
        mctrees[file.first]->SetBranchAddress("jet1", &jet1);
        mctrees[file.first]->SetBranchAddress("jet2", &jet2);

        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("pho1_motherID", &pho1_motherID);
        mctrees[file.first]->SetBranchAddress("pho1_sigmaietaieta", &pho1_sigmaietaieta);
        mctrees[file.first]->SetBranchAddress("NPV", &nVtx); // enable 
    }
    for(auto &file : datafiles){
        datatrees[file.first] = (TTree*)file.second->Get("ControlSampleEvent");
        
	datatrees[file.first]->SetBranchAddress("HLTDecision", HLTDecision);
        datatrees[file.first]->SetBranchAddress(Form("MET%s", suffixes[file.first].c_str()), &mets[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("METPhi%s", suffixes[file.first].c_str()), &metphis[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("MR%s", suffixes[file.first].c_str()), &mrs[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("Rsq%s", suffixes[file.first].c_str()), &rsqs[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("NJets%s", suffixes[file.first].c_str()), &njets40[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("NJets80%s", suffixes[file.first].c_str()), &njets80[file.first]);
	datatrees[file.first]->SetBranchAddress(Form("HT%s", suffixes[file.first].c_str()), &hts[file.first]);
        datatrees[file.first]->SetBranchAddress("NPV", &nVtx); // enable 

        datatrees[file.first]->SetBranchAddress("pho1", &pho1);
        datatrees[file.first]->SetBranchAddress("pho2", &pho2);
        datatrees[file.first]->SetBranchAddress("jet1", &jet1);
        datatrees[file.first]->SetBranchAddress("jet2", &jet2);
        datatrees[file.first]->SetBranchAddress("pho1_sigmaietaieta", &pho1_sigmaietaieta);
  }

     //lumis for Photon paths
    float lumi_HLTPhoton50  = lumi/380.; // prescale 2850
    float lumi_HLTPhoton75  = lumi/40.; // prescale 600
    float lumi_HLTPhoton90  = lumi/20.; // prescale 300
    float lumi_HLTPhoton120 = lumi/7.; // prescale 100
    float lumi_HLTPhoton175 = lumi/1.; // prescale 1
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
    TFile *pileupWeightFile = new TFile("./NVtx_Run2015C_DoubleEG.root", "READ");
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

    cutSequence.push_back( "pho1->Pt() > 30 && pho2->Pt() > 15" );
    cutName.push_back( "" );

    map<string, vector<TH1F *> > mcNJets40, mcNJets80, mcMR, mcRsq,  mcMet, mcNvtx,  mcHT, mcPho1Pt, mcnPho, mcPhoSiSi_EB, 
      mcPhoSiSi_EE, mcjet1Pt, mcjet2Pt, mcjet1Eta, mcjet2Eta, mcu1, mcu2, mcZmassEB, mcZmassEE;
    vector<TH1F *>  dataNJets40, dataNJets80, dataMR, dataRsq, dataMet, dataNvtx, dataHT,  dataPho1Pt, datanPho, dataPhoSiSi_EB, 
      dataPhoSiSi_EE, datajet1Pt, datajet2Pt, datajet1Eta, datajet2Eta, datau1, datau2, dataZmassEB, dataZmassEE;

   for(auto &tree : mctrees){
        mcNvtx[tree.first] = vector<TH1F *>();
        mcNJets40[tree.first] = vector<TH1F *>();
        mcNJets80[tree.first] = vector<TH1F *>();
        mcMR[tree.first] = vector<TH1F *>();
        mcRsq[tree.first] = vector<TH1F *>();
        mcMet[tree.first] = vector<TH1F *>();
        mcHT[tree.first] = vector<TH1F *>();
        mcPho1Pt[tree.first] = vector<TH1F *>();
        mcjet1Pt[tree.first] = vector<TH1F *>();
        mcjet2Pt[tree.first] = vector<TH1F *>();
        mcjet1Eta[tree.first] = vector<TH1F *>();
        mcjet2Eta[tree.first] = vector<TH1F *>();
        mcnPho[tree.first] = vector<TH1F *>();
        mcPhoSiSi_EB[tree.first] = vector<TH1F *>();
        mcPhoSiSi_EE[tree.first] = vector<TH1F *>();
        mcu1[tree.first] = vector<TH1F *>();
	mcu2[tree.first] = vector<TH1F *>();
	mcZmassEB[tree.first] = vector<TH1F *>();
	mcZmassEE[tree.first] = vector<TH1F *>();
   }
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        for(auto &tree : mctrees){
            mcNvtx[tree.first].push_back(new TH1F(Form("mcNvtx%s%d", tree.first.c_str(), cut), Form("%s; NVtx (GeV)", cutName[cut].c_str()), 50, 0, 50));
            mcNJets40[tree.first].push_back(new TH1F(Form("mcNJets40%s%d", tree.first.c_str(), cut), Form("%s; Number of jets 40 GeV", cutName[cut].c_str()), 10, 0, 10));
            mcNJets80[tree.first].push_back(new TH1F(Form("mcNJets80%s%d", tree.first.c_str(), cut), Form("%s; Number of jets 80 GeV", cutName[cut].c_str()), 10, 0, 10));
            mcMR[tree.first].push_back(new TH1F(Form("mcMR%s%d", tree.first.c_str(), cut), Form("%s; MR [GeV/c^2]", cutName[cut].c_str()), 100, 0, 2000));
            mcRsq[tree.first].push_back(new TH1F(Form("mcRsq%s%d", tree.first.c_str(), cut), Form("%s; Rsq (GeV)", cutName[cut].c_str()), nRsqBins, RsqBinLowEdges));
            mcMet[tree.first].push_back(new TH1F(Form("mcMet%s%d", tree.first.c_str(), cut), Form("%s; MET (GeV)", cutName[cut].c_str()), 200, 0, 500));
            mcHT[tree.first].push_back(new TH1F(Form("mcHT%s%d", tree.first.c_str(), cut), Form("%s; HT (GeV)", cutName[cut].c_str()), 100, 0, 2000));
            mcPho1Pt[tree.first].push_back(new TH1F(Form("mcPho1Pt%s%d", tree.first.c_str(), cut), Form("%s; Photon Pt (GeV)", cutName[cut].c_str()), 100, 0, 500));
            mcjet1Pt[tree.first].push_back(new TH1F(Form("mcJet1Pt%s%d", tree.first.c_str(), cut), Form("%s; Jet 1 Pt (GeV)", cutName[cut].c_str()), 100, 0, 1000));
            mcjet2Pt[tree.first].push_back(new TH1F(Form("mcJet2Pt%s%d", tree.first.c_str(), cut), Form("%s; Jet 2 Pt (GeV)", cutName[cut].c_str()), 100, 0, 1000));
            mcjet1Eta[tree.first].push_back(new TH1F(Form("mcJet1Eta%s%d", tree.first.c_str(), cut), Form("%s; Jet 1 Eta", cutName[cut].c_str()), 60, -5, 5));
            mcjet2Eta[tree.first].push_back(new TH1F(Form("mcJet2Eta%s%d", tree.first.c_str(), cut), Form("%s; Jet 2 Eta", cutName[cut].c_str()), 60, -5, 5));
            mcnPho[tree.first].push_back(new TH1F(Form("mcnPho%s%d", tree.first.c_str(), cut), Form("%s; NPhoton", cutName[cut].c_str()), 10, 0, 10));
            mcPhoSiSi_EB[tree.first].push_back(new TH1F(Form("mcPhoSiSi_EB%s%d", tree.first.c_str(), cut), Form("%s; #sigma_{i#etai#eta}", cutName[cut].c_str()), 100, 0, 0.04));
            mcPhoSiSi_EE[tree.first].push_back(new TH1F(Form("mcPhoSiSi_EE%s%d", tree.first.c_str(), cut), Form("%s; #sigma_{i#etai#eta}", cutName[cut].c_str()), 100, 0, 0.04));
            mcu1[tree.first].push_back(new TH1F(Form("mcu1%s%d", tree.first.c_str(), cut), Form("%s; u1", cutName[cut].c_str()), 100, -600, 200));
            mcu2[tree.first].push_back(new TH1F(Form("mcu2%s%d", tree.first.c_str(), cut), Form("%s; u2", cutName[cut].c_str()), 100, -200, 200));
            mcZmassEB[tree.first].push_back(new TH1F(Form("mcZmassEB%s%d", tree.first.c_str(), cut), Form("%s; Mass", cutName[cut].c_str()), 60, 60, 120));
            mcZmassEE[tree.first].push_back(new TH1F(Form("mcZmassEE%s%d", tree.first.c_str(), cut), Form("%s; Mass", cutName[cut].c_str()), 60, 60, 120));

            mcNJets40[tree.first][cut]->Sumw2();
            mcNJets80[tree.first][cut]->Sumw2();
            mcMR[tree.first][cut]->Sumw2();
            mcRsq[tree.first][cut]->Sumw2();
            mcMet[tree.first][cut]->Sumw2();
            mcNvtx[tree.first][cut]->Sumw2();
	    mcHT[tree.first][cut]->Sumw2();	    
	    mcPho1Pt[tree.first][cut]->Sumw2();	    
	    mcjet1Pt[tree.first][cut]->Sumw2();	    
	    mcjet2Pt[tree.first][cut]->Sumw2();	    
	    mcjet1Eta[tree.first][cut]->Sumw2();	    
	    mcjet2Eta[tree.first][cut]->Sumw2();	    
	    mcnPho[tree.first][cut]->Sumw2();
	    mcPhoSiSi_EB[tree.first][cut]->Sumw2();	    
	    mcPhoSiSi_EE[tree.first][cut]->Sumw2();	    
 	    mcu1[tree.first][cut]->Sumw2();	    
 	    mcu2[tree.first][cut]->Sumw2();	    
	    mcZmassEB[tree.first][cut]->Sumw2(); 
 	    mcZmassEE[tree.first][cut]->Sumw2(); 
      }
        dataNJets40.push_back(new TH1F(Form("dataNJets40%d", cut), Form("%s; Number of jets 40 GeV", cutName[cut].c_str()), 10, 0, 10));
        dataNJets80.push_back(new TH1F(Form("dataNJets80%d", cut), Form("%s; Number of jets 80 GeV", cutName[cut].c_str()), 10, 0, 10));
        dataMR.push_back(new TH1F(Form("dataMR%d", cut), Form("%s; MR [GeV/c^2]", cutName[cut].c_str()), 100, 0, 2000));
        dataRsq.push_back(new TH1F(Form("dataRsq%d", cut), Form("%s; Rsq (GeV)", cutName[cut].c_str()), nRsqBins, RsqBinLowEdges));
        dataMet.push_back(new TH1F(Form("dataMet%d", cut), Form("%s; mcMet (GeV)", cutName[cut].c_str()), 200, 0, 500));
        dataNvtx.push_back(new TH1F(Form("dataNvtx%d", cut), Form("%s; NVtx (GeV)", cutName[cut].c_str()), 50, 0, 50));
        dataHT.push_back(new TH1F(Form("dataHT%d", cut), Form("%s; HT (GeV)", cutName[cut].c_str()), 100, 0, 2000));
        dataPho1Pt.push_back(new TH1F(Form("datapho1Pt%d", cut), Form("%s; Photon Pt (GeV)", cutName[cut].c_str()), 100, 0, 500));
        datajet1Pt.push_back(new TH1F(Form("datajet1Pt%d", cut), Form("%s; Jet 1 Pt (GeV)", cutName[cut].c_str()), 100, 0, 1000));
        datajet2Pt.push_back(new TH1F(Form("datajet2Pt%d", cut), Form("%s; Jet 2 Pt (GeV)", cutName[cut].c_str()), 100, 0, 1000));
        datajet1Eta.push_back(new TH1F(Form("datajet1Eta%d", cut), Form("%s; Jet 1 Eta", cutName[cut].c_str()), 100, -5, 5));
        datajet2Eta.push_back(new TH1F(Form("datajet2Eta%d", cut), Form("%s; Jet 2 Eta", cutName[cut].c_str()), 100, -5, 5));
        datanPho.push_back(new TH1F(Form("datanPho%d", cut), Form("%s; NPhoton", cutName[cut].c_str()), 10, 0, 10));
        dataPhoSiSi_EB.push_back(new TH1F(Form("dataPhoSiSi_EB%d", cut), Form("%s; #sigma_{i#etai#eta}", cutName[cut].c_str()), 100, 0, 0.04));
        dataPhoSiSi_EE.push_back(new TH1F(Form("dataPhoSiSi_EE%d", cut), Form("%s; #sigma_{i#etai#eta}", cutName[cut].c_str()), 100, 0, 0.04));
        datau1.push_back(new TH1F(Form("datau1%d", cut), Form("%s; u1", cutName[cut].c_str()), 100, -600, 200));
        datau2.push_back(new TH1F(Form("datau2%d", cut), Form("%s; u2", cutName[cut].c_str()), 100, -200, 200));
        dataZmassEB.push_back(new TH1F(Form("dataZmassEB%d", cut), Form("%s; Mass", cutName[cut].c_str()), 60,  60, 120));
        dataZmassEE.push_back(new TH1F(Form("dataZmassEE%d", cut), Form("%s; Mass", cutName[cut].c_str()), 60,  60, 120));

        dataNJets40[cut]->Sumw2();
        dataNJets80[cut]->Sumw2();
        dataMR[cut]->Sumw2();
        dataRsq[cut]->Sumw2();
        dataMet[cut]->Sumw2();
	dataNvtx[cut]->Sumw2();
	dataHT[cut]->Sumw2();
	dataPho1Pt[cut]->Sumw2();
	datajet1Pt[cut]->Sumw2();
	datajet2Pt[cut]->Sumw2();
	datajet1Eta[cut]->Sumw2();
	datajet2Eta[cut]->Sumw2();
	datanPho[cut]->Sumw2(); 
	dataPhoSiSi_EB[cut]->Sumw2(); 
	dataPhoSiSi_EE[cut]->Sumw2(); 
	datau1[cut]->Sumw2(); 
	datau2[cut]->Sumw2(); 
 	dataZmassEB[cut]->Sumw2(); 
 	dataZmassEE[cut]->Sumw2(); 
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
	    eventWeight *= pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(nVtx));

	    bool trigger_passed = false;
	    if (HLTDecision[28] || HLTDecision[29] ) trigger_passed = true; //electron triggers
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
		(mcPho1Pt[tree.first])[cut]->Fill(pho1->Pt(), eventWeight);
 		(mcnPho[tree.first])[cut]->Fill(nSelectedPhotons, eventWeight);
 		if( fabs(pho1->Eta()) < 1.44 )(mcPhoSiSi_EB[tree.first])[cut]->Fill(pho1_sigmaietaieta, eventWeight);
 		if( fabs(pho1->Eta()) > 1.566 &&  fabs(pho1->Eta()) <2.5 )(mcPhoSiSi_EE[tree.first])[cut]->Fill(pho1_sigmaietaieta, eventWeight);

		(mcjet1Pt[tree.first])[cut]->Fill(jet1->Pt(), eventWeight);
		(mcjet2Pt[tree.first])[cut]->Fill(jet2->Pt(), eventWeight);
		(mcjet1Eta[tree.first])[cut]->Fill(jet1->Eta(), eventWeight);
		(mcjet2Eta[tree.first])[cut]->Fill(jet2->Eta(), eventWeight);

		// Compute the recoil for Gamma
		TVector2 vGPt(pho1->Pt()*cos(pho1->Phi()),pho1->Pt()*sin(pho1->Phi())); // Gamma boson pT
		TVector2 vMet(mets[tree.first]*cos(metphis[tree.first]), mets[tree.first]*sin(metphis[tree.first])); // MET vector
		TVector2 vU = -1.0*(vGPt+vMet); // recoil vector
		float u1 = ((pho1->Px())*(vU.Px()) + (pho1->Py())*(vU.Py()))/(pho1->Pt());  // u1 = (pT . u)/|pT|
		float u2 = ((pho1->Px())*(vU.Py()) - (pho1->Py())*(vU.Px()))/(pho1->Pt());  // u2 = (pT x u)/|pT|

		(mcu1[tree.first])[cut]->Fill(u1+pho1->Pt(), eventWeight);
		(mcu2[tree.first])[cut]->Fill(u2, eventWeight);
		
		if(fabs(pho1->Eta())<1.4442 && fabs(pho2->Eta())<1.4442)
		   (mcZmassEB[tree.first])[cut]->Fill( (*pho1 + *pho2).M(), eventWeight);
		
		if( 
		   (fabs(pho1->Eta())>1.566 && fabs(pho1->Eta())<2.5) 
		   ||
		   (fabs(pho2->Eta())>1.566 && fabs(pho2->Eta())<2.5)
		    )
		  (mcZmassEE[tree.first])[cut]->Fill( (*pho1 + *pho2).M(), eventWeight);
	    }
        }

        for(uint cut = 0; cut < cutSequence.size(); cut++){
            delete cuts[cut];
        }
    }

    for(auto &tree : datatrees){
        cout << "Filling data histograms: " << tree.first << endl;
        uint nEntries = tree.second->GetEntries();
        if(debug) nEntries = 10000;
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
	    if (HLTDecision[28] || HLTDecision[29]) trigger_passed = true; //electron triggers
	    if (trigger_passed == false) continue;

            //apply selection cuts and fill the appropriate histograms
            for(uint cut = 0; cut < cutSequence.size(); cut++){
                bool passesCut = cuts[cut]->EvalInstance();
                if(!passesCut) continue;

		TLorentzVector photon1;
		TLorentzVector photon2;
		photon1.SetPtEtaPhiM(pho1->Pt()*GetElectronScaleCorrection(pho1->Pt(),pho1->Eta()),
				     pho1->Eta(),
				     pho1->Phi(),
				     pho1->M());
		photon2.SetPtEtaPhiM(pho2->Pt()*GetElectronScaleCorrection(pho2->Pt(),pho2->Eta()),
				     pho2->Eta(),
				     pho2->Phi(),
				     pho2->M());

		dataNJets40[cut]->Fill(njets40[tree.first], eventWeight);
                dataNJets80[cut]->Fill(njets80[tree.first], eventWeight);
                dataMR[cut]->Fill(mrs[tree.first], eventWeight);
                dataRsq[cut]->Fill(rsqs[tree.first], eventWeight);
                dataMet[cut]->Fill(mets[tree.first], eventWeight);
		dataNvtx[cut]->Fill(nVtx, eventWeight);
		dataHT[cut]->Fill(hts[tree.first], eventWeight);
		dataPho1Pt[cut]->Fill(photon1.Pt(), eventWeight);
 		datanPho[cut]->Fill(nSelectedPhotons, eventWeight);
 		if( fabs(photon1.Eta()) < 1.44 ) dataPhoSiSi_EB[cut]->Fill(pho1_sigmaietaieta, eventWeight);
 		if( fabs(photon1.Eta()) > 1.566 &&  fabs(photon1.Eta()) <2.5 ) dataPhoSiSi_EE[cut]->Fill(pho1_sigmaietaieta, eventWeight);

		datajet1Pt[cut]->Fill(jet1->Pt(), eventWeight);
		datajet2Pt[cut]->Fill(jet2->Pt(), eventWeight);
		datajet1Eta[cut]->Fill(jet1->Eta(), eventWeight);
		datajet2Eta[cut]->Fill(jet2->Eta(), eventWeight);
		
		if(fabs(pho1->Eta())<1.4442 && fabs(pho2->Eta())<1.4442)
		  dataZmassEB[cut]->Fill((photon1+photon2).M(), eventWeight);
		
		if( 
		   (fabs(pho1->Eta())>1.566 && fabs(pho1->Eta())<2.5) 
		   ||
		   (fabs(pho2->Eta())>1.566 && fabs(pho2->Eta())<2.5)
		    )
		  dataZmassEE[cut]->Fill((photon1+photon2).M(), eventWeight);
	    }
        }
        for(uint cut = 0; cut < cutSequence.size(); cut++){
            delete cuts[cut];
        }
    }

    //print out plots
    TCanvas c("c", "c", 800, 700);
    c.SetLogy();
    
    //    dataNvtx[0]->Scale(4917.77/4518.);

    //colors and legend
    map<string, int> colors;
    colors["ZJets"] = kCyan+2;
    TLegend *legend = new TLegend(0.7, 0.8, 0.95, 0.95);
    legend->AddEntry(dataNJets40[0], "Data");
    legend->AddEntry(mcNJets40["ZJets"][0], "DY + Jets MC");
    cout<<"HERRE 2"<<endl;
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        //create histogram stacks for MC
        THStack NumJets40MC(Form("NumJets40Stack%d", cut), cutName[cut].c_str());
        THStack NumJets80MC(Form("NumJets80Stack%d", cut), cutName[cut].c_str());
        THStack MRMC(Form("MRStack%d", cut), cutName[cut].c_str());
        THStack RsqMC(Form("RsqStack%d", cut), cutName[cut].c_str());
        THStack MetMC(Form("MetStack%d", cut), cutName[cut].c_str());
        THStack NVtxMC(Form("NVtxStack%d", cut), cutName[cut].c_str());
	THStack HTMC(Form("HTStack%d", cut), cutName[cut].c_str());
	THStack PHO1PTMC(Form("PHO1PTMC%d", cut), cutName[cut].c_str());
	THStack NPHOMC(Form("NPHOMC%d", cut), cutName[cut].c_str());
	THStack PHO1SiSi_EB(Form("PHO1SiSi_EB%d", cut), cutName[cut].c_str());
	THStack PHO1SiSi_EE(Form("PHO1SiSi_EE%d", cut), cutName[cut].c_str());
	THStack JET1PTMC(Form("JET1PTMC%d", cut), cutName[cut].c_str());
	THStack JET2PTMC(Form("JET2PTMC%d", cut), cutName[cut].c_str());
	THStack JET1ETAMC(Form("JET1ETAMC%d", cut), cutName[cut].c_str());
	THStack JET2ETAMC(Form("JET2ETAMC%d", cut), cutName[cut].c_str());
	THStack U1MC(Form("U1MC%d", cut), cutName[cut].c_str());
	THStack U2MC(Form("U2MC%d", cut), cutName[cut].c_str());
	THStack ZMMCEB(Form("ZMMCEB%d", cut), cutName[cut].c_str());
	THStack ZMMCEE(Form("ZMMCEE%d", cut), cutName[cut].c_str());

        //add the histograms to the stack in order
        vector<string> orderedtrees {"ZJets"};
        for(auto &tree : orderedtrees){
	    mcNJets40[tree][cut]->SetFillColor(colors[tree]);
            mcNJets80[tree][cut]->SetFillColor(colors[tree]);
            mcMR[tree][cut]->SetFillColor(colors[tree]);
            mcRsq[tree][cut]->SetFillColor(colors[tree]);
            mcMet[tree][cut]->SetFillColor(colors[tree]);
            mcNvtx[tree][cut]->SetFillColor(colors[tree]);
	    mcHT[tree][cut]->SetFillColor(colors[tree]);
	    mcPho1Pt[tree][cut]->SetFillColor(colors[tree]);
	    mcnPho[tree][cut]->SetFillColor(colors[tree]);
	    mcPhoSiSi_EB[tree][cut]->SetFillColor(colors[tree]);
	    mcPhoSiSi_EE[tree][cut]->SetFillColor(colors[tree]);
	    mcjet1Pt[tree][cut]->SetFillColor(colors[tree]);
	    mcjet2Pt[tree][cut]->SetFillColor(colors[tree]);
	    mcjet1Eta[tree][cut]->SetFillColor(colors[tree]);
	    mcjet2Eta[tree][cut]->SetFillColor(colors[tree]);
	    mcu1[tree][cut]->SetFillColor(colors[tree]);
	    mcu2[tree][cut]->SetFillColor(colors[tree]);
	    mcZmassEB[tree][cut]->SetFillColor(colors[tree]);
	    mcZmassEE[tree][cut]->SetFillColor(colors[tree]);

            NumJets40MC.Add(mcNJets40[tree][cut]);
            NumJets80MC.Add(mcNJets80[tree][cut]);
            MRMC.Add(mcMR[tree][cut]);
            RsqMC.Add(mcRsq[tree][cut]);
	    MetMC.Add(mcMet[tree][cut]);
	    NVtxMC.Add(mcNvtx[tree][cut]);
	    HTMC.Add(mcHT[tree][cut]);
	    PHO1PTMC.Add(mcPho1Pt[tree][cut]);
 	    NPHOMC.Add(mcnPho[tree][cut]);
 	    PHO1SiSi_EB.Add(mcPhoSiSi_EB[tree][cut]);
 	    PHO1SiSi_EE.Add(mcPhoSiSi_EE[tree][cut]);
	    JET1PTMC.Add(mcjet1Pt[tree][cut]);
	    JET2PTMC.Add(mcjet2Pt[tree][cut]);
	    JET1ETAMC.Add(mcjet1Eta[tree][cut]);
	    JET2ETAMC.Add(mcjet2Eta[tree][cut]);
	    U1MC.Add(mcu1[tree][cut]);
 	    U2MC.Add(mcu2[tree][cut]);
 	    ZMMCEB.Add(mcZmassEB[tree][cut]);
 	    ZMMCEE.Add(mcZmassEE[tree][cut]);
      }
	DrawDataVsMCRatioPlot(dataNJets40[cut], &NumJets40MC, legend, "Number of jets 40 GeV", "ControlRegionPlots_Zee__NumJets40"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataNJets80[cut], &NumJets80MC, legend, "Number of jets 80 GeV", "ControlRegionPlots_Zee__NumJets80"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataMR[cut], &MRMC, legend, "MR (GeV)", "ControlRegionPlots_Zee__MR"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataRsq[cut], &RsqMC, legend, "Rsq", "ControlRegionPlots_Zee__Rsq"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataMet[cut], &MetMC, legend, "Met", "ControlRegionPlots_Zee__MET"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataNvtx[cut], &NVtxMC, legend, "NVtx", "ControlRegionPlots_Zee__NVtx"+to_string(cut), false);
	DrawDataVsMCRatioPlot(dataHT[cut], &HTMC, legend, "HT", "ControlRegionPlots_Zee__HT"+to_string(cut), false);
	DrawDataVsMCRatioPlot(dataPho1Pt[cut], &PHO1PTMC, legend, "Pho1Pt", "ControlRegionPlots_Zee__Pho1Pt"+to_string(cut), false);
	DrawDataVsMCRatioPlot(datanPho[cut], &NPHOMC, legend, "NPho", "ControlRegionPlots_Zee__nPho"+to_string(cut), false);
	DrawDataVsMCRatioPlot(dataPhoSiSi_EB[cut], &PHO1SiSi_EB, legend, "Pho1SiSi_EB", "ControlRegionPlots_Zee__pho1SiSi_EB"+to_string(cut), false);
	DrawDataVsMCRatioPlot(dataPhoSiSi_EE[cut], &PHO1SiSi_EE, legend, "Pho1SiSi_EE", "ControlRegionPlots_Zee__pho1SiSi_EE"+to_string(cut), false);

	DrawDataVsMCRatioPlot(datajet1Pt[cut], &JET1PTMC, legend, "Jet1Pt", "ControlRegionPlots_Zee__Jet1Pt"+to_string(cut), false);
	DrawDataVsMCRatioPlot(datajet2Pt[cut], &JET2PTMC, legend, "Jet2Pt", "ControlRegionPlots_Zee__Jet2Pt"+to_string(cut), false);
	DrawDataVsMCRatioPlot(datajet1Eta[cut], &JET1ETAMC, legend, "Jet1Eta", "ControlRegionPlots_Zee__Jet1Eta"+to_string(cut), false);
	DrawDataVsMCRatioPlot(datajet2Eta[cut], &JET2ETAMC, legend, "Jet2Eta", "ControlRegionPlots_Zee__Jet2Eta"+to_string(cut), false);
 	DrawDataVsMCRatioPlot(datau1[cut], &U1MC, legend, "U1", "ControlRegionPlots_Zee__U1"+to_string(cut), false);
  	DrawDataVsMCRatioPlot(datau2[cut], &U2MC, legend, "U2", "ControlRegionPlots_Zee__U2"+to_string(cut), false);
  	DrawDataVsMCRatioPlot(dataZmassEB[cut], &ZMMCEB, legend, "ZMEB", "ControlRegionPlots_Zee__MassEB"+to_string(cut), false, false);
  	DrawDataVsMCRatioPlot(dataZmassEE[cut], &ZMMCEE, legend, "ZMEE", "ControlRegionPlots_Zee__MassEE"+to_string(cut), false, false);
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
	delete dataPho1Pt[cut];
	delete datajet1Pt[cut];
	delete datajet2Pt[cut];
	delete datajet1Eta[cut];
	delete datajet2Eta[cut];
	delete datanPho[cut];
	delete dataPhoSiSi_EB[cut];
	delete dataPhoSiSi_EE[cut];
	delete datau1[cut];
	delete datau2[cut];
	delete dataZmassEB[cut];
	delete dataZmassEE[cut];
	for(auto &tree : mctrees){
	  delete mcNJets40[tree.first][cut];
	  delete mcNJets80[tree.first][cut];
	  delete mcMR[tree.first][cut];
	  delete mcRsq[tree.first][cut];
	  delete mcMet[tree.first][cut];
	  delete mcNvtx[tree.first][cut];
	  delete mcHT[tree.first][cut];
	  delete mcPho1Pt[tree.first][cut];
	  delete mcjet1Pt[tree.first][cut];
	  delete mcjet2Pt[tree.first][cut];
	  delete mcjet1Eta[tree.first][cut];
	  delete mcjet2Eta[tree.first][cut];
	  delete mcnPho[tree.first][cut];
	  delete mcPhoSiSi_EB[tree.first][cut];
	  delete mcPhoSiSi_EE[tree.first][cut];
 	  delete mcu1[tree.first][cut];
 	  delete mcu2[tree.first][cut];
 	  delete mcZmassEB[tree.first][cut];
 	  delete mcZmassEE[tree.first][cut];
       }
    }

    mcfiles["ZJets"]->Close();
    datafiles["ZJets"]->Close();
    pileupWeightFile->Close();
    // effFile->Close();
    delete mcfiles["ZJets"];
    delete datafiles["ZJets"];
    delete pileupWeightFile;
    // delete effFile;
    cout<<"HERRE 5"<<endl;

}

int main(){
    makeCRplots_GJets();
    return 0;
}
