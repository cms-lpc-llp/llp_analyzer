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
// bool debug = true;

//lepton Pt, Yields, NJets 40, Electronsxs

//true: find the translation factors from MC to data
//false: find the translation factors from DY, W, G to Z->nunu
bool computeDataOverMCSFs = false;
double deltaPhi(double phi1, double phi2);

void MTplots(){
    gROOT->SetBatch();

    float maxMuonPt = 999;

    //set color palette 
    const Int_t NCont = 100;
    gStyle->SetNumberContours(NCont);
    gStyle->SetPaintTextFormat("1.0f");

    const int lumi = 1264.;
    
    map<string, string> suffixes;
    suffixes["WJets"] = "";

    //get input files -- assumes one TFile for each process, with weights for different HT bins 
    map<string, TFile*> mcfiles;

    // Reduced
    // mcfiles["WJets"] = new TFile("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p21/OneLeptonMT/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
    mcfiles["WJets"] = new TFile("eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p21/OneLeptonMT/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");

    //get trees and set branches
    map<string, TTree*> mctrees;
    map<string, Float_t> mets;
    map<string, Float_t> metphis;
    
    float weight;
    UInt_t nPU_mean;
    float leadingMuonPt, leadingMuonEta, recoZpt, recoZeta, recoZmass, subleadingMuonPt, subleadingMuonEta, mTLepMet, mTLepMetNoHF;
    UInt_t nVtx, nBTaggedJets, nSelectedPhotons;
    Bool_t HLTNames[156];
    Bool_t lep1passTight;
    UInt_t event;
    TLorentzVector* lep1 = NULL; 
    Int_t lep1Type;

  float metType1PtJetResUp;
  float metType1PtJetResDown;
  float metType1PtJetEnUp;
  float metType1PtJetEnDown;
  float metType1PtMuonEnUp;
  float metType1PtMuonEnDown;
  float metType1PtElectronEnUp;
  float metType1PtElectronEnDown;
  float metType1PtTauEnUp;
  float metType1PtTauEnDown;
  float metType1PtUnclusteredEnUp;
  float metType1PtUnclusteredEnDown;
  float metType1PtPhotonEnUp;
  float metType1PtPhotonEnDown;       
  
  float metType1PhiJetResUp;
  float metType1PhiJetResDown;
  float metType1PhiJetEnUp;
  float metType1PhiJetEnDown;
  float metType1PhiMuonEnUp;
  float metType1PhiMuonEnDown;
  float metType1PhiElectronEnUp;
  float metType1PhiElectronEnDown;
  float metType1PhiTauEnUp;
  float metType1PhiTauEnDown;
  float metType1PhiUnclusteredEnUp;
  float metType1PhiUnclusteredEnDown;
  float metType1PhiPhotonEnUp;
  float metType1PhiPhotonEnDown;

    for(auto &file : mcfiles){
        mets[file.first] = 0.;
        metphis[file.first] = 0.;

	mctrees[file.first] = (TTree*)file.second->Get("ControlSampleEvent");

	mctrees[file.first]->SetBranchAddress("lep1", &lep1);
	mctrees[file.first]->SetBranchAddress("HLTDecision", HLTNames);

        mctrees[file.first]->SetBranchAddress(Form("MET%s", suffixes[file.first].c_str()), &mets[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("METPhi%s", suffixes[file.first].c_str()), &metphis[file.first]);
        mctrees[file.first]->SetBranchAddress("event", &event);
        mctrees[file.first]->SetBranchAddress("NBJetsMedium", &nBTaggedJets);
        mctrees[file.first]->SetBranchAddress("lep1MT", &mTLepMet);
        mctrees[file.first]->SetBranchAddress("lep1PassTight", &lep1passTight);
        mctrees[file.first]->SetBranchAddress("lep1Type", &lep1Type);

	//        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("NPV", &nVtx); // enable 
        mctrees[file.first]->SetBranchAddress("NPU_0", &nPU_mean); // enable 

	mctrees[file.first]->SetBranchAddress("metType1PtJetResUp", &metType1PtJetResUp);
	mctrees[file.first]->SetBranchAddress("metType1PtJetResDown", &metType1PtJetResDown);
	mctrees[file.first]->SetBranchAddress("metType1PtJetEnUp", &metType1PtJetEnUp);
	mctrees[file.first]->SetBranchAddress("metType1PtJetEnDown", &metType1PtJetEnDown);
	mctrees[file.first]->SetBranchAddress("metType1PtMuonEnUp", &metType1PtMuonEnUp);
	mctrees[file.first]->SetBranchAddress("metType1PtMuonEnDown", &metType1PtMuonEnDown);
	mctrees[file.first]->SetBranchAddress("metType1PtElectronEnUp", &metType1PtElectronEnUp);
	mctrees[file.first]->SetBranchAddress("metType1PtElectronEnDown", &metType1PtElectronEnDown);
	mctrees[file.first]->SetBranchAddress("metType1PtTauEnUp", &metType1PtTauEnUp);
	mctrees[file.first]->SetBranchAddress("metType1PtTauEnDown", &metType1PtTauEnDown);
	mctrees[file.first]->SetBranchAddress("metType1PtUnclusteredEnUp", &metType1PtUnclusteredEnUp);
	mctrees[file.first]->SetBranchAddress("metType1PtUnclusteredEnDown", &metType1PtUnclusteredEnDown);
	mctrees[file.first]->SetBranchAddress("metType1PtPhotonEnUp", &metType1PtPhotonEnUp);
	mctrees[file.first]->SetBranchAddress("metType1PtPhotonEnDown", &metType1PtPhotonEnDown);
	
	mctrees[file.first]->SetBranchAddress("metType1PhiJetResUp", &metType1PhiJetResUp);
	mctrees[file.first]->SetBranchAddress("metType1PhiJetResDown", &metType1PhiJetResDown);
	mctrees[file.first]->SetBranchAddress("metType1PhiJetEnUp", &metType1PhiJetEnUp);
	mctrees[file.first]->SetBranchAddress("metType1PhiJetEnDown", &metType1PhiJetEnDown);
	mctrees[file.first]->SetBranchAddress("metType1PhiMuonEnUp", &metType1PhiMuonEnUp);
	mctrees[file.first]->SetBranchAddress("metType1PhiMuonEnDown", &metType1PhiMuonEnDown);
	mctrees[file.first]->SetBranchAddress("metType1PhiElectronEnUp", &metType1PhiElectronEnUp);
	mctrees[file.first]->SetBranchAddress("metType1PhiElectronEnDown", &metType1PhiElectronEnDown);
	mctrees[file.first]->SetBranchAddress("metType1PhiTauEnUp", &metType1PhiTauEnUp);
	mctrees[file.first]->SetBranchAddress("metType1PhiTauEnDown", &metType1PhiTauEnDown);
	mctrees[file.first]->SetBranchAddress("metType1PhiUnclusteredEnUp", &metType1PhiUnclusteredEnUp);
	mctrees[file.first]->SetBranchAddress("metType1PhiUnclusteredEnDown", &metType1PhiUnclusteredEnDown);
	mctrees[file.first]->SetBranchAddress("metType1PhiPhotonEnUp", &metType1PhiPhotonEnUp);
	mctrees[file.first]->SetBranchAddress("metType1PhiPhotonEnDown", &metType1PhiPhotonEnDown);
   }

    //load pileup reweighting histogram
    TFile *pileupWeightFile = new TFile("eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PileupWeights/NVtxReweight_ZToMuMu_2015D_911ipb.root", "READ");
    TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("NVtxReweight");
    assert(pileupWeightHist);

    //define cuts and histograms
    float nMRBins = 10;
    float nRsqBins = 8;
    float MRBinLowEdges[] = {300, 350, 400, 450, 550, 700, 900, 1200, 1600, 2500, 4000};
    float RsqBinLowEdges[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 0.80, 1.5};
    vector<string> cutSequence;
    vector<string> cutName;

    cutSequence.push_back( "NBJetsMedium == 0 && lep1PassTight == 1 && lep1->Pt() > 25 && TMath::Abs(lep1->Eta())<2.5 " );
    cutName.push_back( "" );

    cutSequence.push_back( "NBJetsMedium == 0 && lep1PassTight == 1 && lep1->Pt() > 25 && TMath::Abs(lep1->Eta())<2.5 && MR > 400 && MR < 600 " );
    cutName.push_back( "" );
    
    cutSequence.push_back( "NBJetsMedium == 0 && lep1PassTight == 1 && lep1->Pt() > 25 && TMath::Abs(lep1->Eta())<2.5 && MR > 600 && MR < 800 " );
    cutName.push_back( "" );

    cutSequence.push_back( "NBJetsMedium == 0 && lep1PassTight == 1 && lep1->Pt() > 25 && TMath::Abs(lep1->Eta())<2.5 && MR > 800 && MR < 1000 " );
    cutName.push_back( "" );

    cutSequence.push_back( "NBJetsMedium == 0 && lep1PassTight == 1 && lep1->Pt() > 25 && TMath::Abs(lep1->Eta())<2.5 && MR > 1000 " );
    cutName.push_back( "" );

   map<string, vector<TH1F *> > mcMT, mcMTLepScaleUp, mcMTLepScaleDown, mcMTJetResUp, mcMTJetResDown, mcMTJetEnUp, mcMTJetEnDown,mcMTMuonEnUp,mcMTMuonEnDown,mcMTElectronEnUp,mcMTElectronEnDown,mcMTTauEnUp,mcMTTauEnDown,mcMTUnclusteredEnUp,mcMTUnclusteredEnDown,mcMTPhotonEnUp,mcMTPhotonEnDown;

    for(auto &tree : mctrees){
        mcMT[tree.first] = vector<TH1F *>();

        mcMTLepScaleUp[tree.first] = vector<TH1F *>();
        mcMTLepScaleDown[tree.first] = vector<TH1F *>();

        mcMTJetResUp[tree.first] = vector<TH1F *>();
        mcMTJetResDown[tree.first] = vector<TH1F *>();

        mcMTJetEnUp[tree.first] = vector<TH1F *>();
        mcMTJetEnDown[tree.first] = vector<TH1F *>();

        mcMTMuonEnUp[tree.first] = vector<TH1F *>();
        mcMTMuonEnDown[tree.first] = vector<TH1F *>();

        mcMTElectronEnUp[tree.first] = vector<TH1F *>();
        mcMTElectronEnDown[tree.first] = vector<TH1F *>();

        mcMTTauEnUp[tree.first] = vector<TH1F *>();
        mcMTTauEnDown[tree.first] = vector<TH1F *>();

        mcMTUnclusteredEnUp[tree.first] = vector<TH1F *>();
        mcMTUnclusteredEnDown[tree.first] = vector<TH1F *>();

        mcMTPhotonEnUp[tree.first] = vector<TH1F *>();
        mcMTPhotonEnDown[tree.first] = vector<TH1F *>();

    }
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        for(auto &tree : mctrees){
            mcMT[tree.first].push_back(new TH1F(Form("mcMT%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTLepScaleUp[tree.first].push_back(new TH1F(Form("mcMTLepScaleUp%s%d", tree.first.c_str(), cut), Form("%s; MT Lep Scale Up (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTLepScaleDown[tree.first].push_back(new TH1F(Form("mcMTLepScaleDown%s%d", tree.first.c_str(), cut), Form("%s; MT Lep Scale Down (GeV)", cutName[cut].c_str()), 50, 0, 500));
	    mcMTJetResUp [tree.first].push_back(new TH1F(Form("mcMTJetResUp%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
	    mcMTJetResDown[tree.first].push_back(new TH1F(Form("mcMTJetResDown%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTJetEnUp[tree.first].push_back(new TH1F(Form("mcMTJetEnUp%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTJetEnDown[tree.first].push_back(new TH1F(Form("mcMTJetEnDown%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTMuonEnUp[tree.first].push_back(new TH1F(Form("mcMTMuonEnUp%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTMuonEnDown[tree.first].push_back(new TH1F(Form("mcMTMuonEnDown%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTElectronEnUp[tree.first].push_back(new TH1F(Form("mcMTElectronEnUp%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTElectronEnDown[tree.first].push_back(new TH1F(Form("mcMTElectronEnDown%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTTauEnUp[tree.first].push_back(new TH1F(Form("mcMTTauEnUp%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTTauEnDown[tree.first].push_back(new TH1F(Form("mcMTTauEnDown%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTUnclusteredEnUp[tree.first].push_back(new TH1F(Form("mcMTUnclusteredEnUp%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTUnclusteredEnDown[tree.first].push_back(new TH1F(Form("mcMTUnclusteredEnDown%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTPhotonEnUp[tree.first].push_back(new TH1F(Form("mcMTPhotonEnUp%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));
            mcMTPhotonEnDown[tree.first].push_back(new TH1F(Form("mcMTPhotonEnDown%s%d", tree.first.c_str(), cut), Form("%s; MT (GeV)", cutName[cut].c_str()), 50, 0, 500));


	    mcMT[tree.first][cut]->Sumw2();
	    mcMTLepScaleUp[tree.first][cut]->Sumw2();
	    mcMTLepScaleDown[tree.first][cut]->Sumw2();
        }
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
        for(uint i = 0; i < 2000000; i++){
            //get entry
            tree.second->GetEntry(i); 

	    if(i % 100000 == 0) cout << "Processing entry " << i << " of "<<tree.first<<endl;
            //get event weight
            //float eventWeight = weight;
            float eventWeight = 1.0;
            eventWeight *= pileupWeightHist->GetBinContent(pileupWeightHist->GetXaxis()->FindFixBin(nPU_mean));

	    bool trigger_passed = false;
	    bool passedMuonLeptonTrigger = bool( HLTNames[2] || HLTNames[7] || HLTNames[12] 
						     || HLTNames[11] || HLTNames[15] );
	    
	    bool passedElectronLeptonTrigger = bool( HLTNames[18] || HLTNames[19] || HLTNames[20] 
						 || HLTNames[21] || HLTNames[28] || HLTNames[29]);

	    if (passedMuonLeptonTrigger && TMath::Abs(lep1Type) == 13 && TMath::Abs(lep1->Eta())<2.4) trigger_passed = true; //muon triggers
	    // if (passedElectronLeptonTrigger && TMath::Abs(lep1Type) == 11 && TMath::Abs(lep1->Eta())<2.4 && lep1->Pt() > 35) trigger_passed = true; //muon triggers
	    
	    if (trigger_passed == false) continue;
	    
	    eventWeight *= lumi;

            //apply selection cuts and fill the appropriate histograms
            for(uint cut = 0; cut < cutSequence.size(); cut++){
                bool passesCut = cuts[cut]->EvalInstance();
                if(!passesCut) continue;

		float lepscaleUPshift = 1.0;
		float lepscaleDOWNshift = 1.0;

		if(TMath::Abs(lep1Type) == 11 && fabs(lep1->Eta())<1.5) { lepscaleUPshift = 1.006; lepscaleDOWNshift = 0.994; }
		if(TMath::Abs(lep1Type) == 11 && fabs(lep1->Eta())>1.5) { lepscaleUPshift = 1.015; lepscaleDOWNshift = 0.985; }

		if(TMath::Abs(lep1Type) == 13 && lep1->Pt()<100) { lepscaleUPshift = 1.002; lepscaleDOWNshift = 0.998; }
		if(TMath::Abs(lep1Type) == 13 && lep1->Pt()>=100) { lepscaleUPshift = 1.05; lepscaleDOWNshift = 0.95; }

		float mcMTLepUp = sqrt(lep1->M2() + 2*mets[tree.first]*lep1->Pt()*1.002*(1 - cos(deltaPhi(metphis[tree.first],lep1->Phi()))));
		float mcMTLepDown = sqrt(lep1->M2() + 2*mets[tree.first]*lep1->Pt()*0.998*(1 - cos(deltaPhi(metphis[tree.first],lep1->Phi()))));

		float MTJetResUp = sqrt(lep1->M2() + 2*metType1PtJetResUp*lep1->Pt()*(1 - cos(deltaPhi(metType1PhiJetResUp,lep1->Phi()))));
		float MTJetResDown = sqrt(lep1->M2() + 2*metType1PtJetResDown*lep1->Pt()*(1 - cos(deltaPhi(metType1PhiJetResDown,lep1->Phi()))));

		float MTJetEnUp = sqrt(lep1->M2() + 2*metType1PtJetEnUp*lep1->Pt()*(1 - cos(deltaPhi(metType1PhiJetEnUp,lep1->Phi()))));
		float MTJetEnDown = sqrt(lep1->M2() + 2*metType1PtJetEnDown*lep1->Pt()*(1 - cos(deltaPhi(metType1PhiJetEnDown,lep1->Phi()))));
		
		float MTMuonEnUp = sqrt(lep1->M2() + 2*metType1PtMuonEnUp*lep1->Pt()*lepscaleUPshift*(1 - cos(deltaPhi(metType1PhiMuonEnUp,lep1->Phi()))));
		float MTMuonEnDown = sqrt(lep1->M2() + 2*metType1PtMuonEnDown*lep1->Pt()*lepscaleDOWNshift*(1 - cos(deltaPhi(metType1PhiMuonEnDown,lep1->Phi()))));

		float MTElectronEnUp = sqrt(lep1->M2() + 2*metType1PtElectronEnUp*lep1->Pt()*lepscaleUPshift*(1 - cos(deltaPhi(metType1PhiElectronEnUp,lep1->Phi()))));
		float MTElectronEnDown = sqrt(lep1->M2() + 2*metType1PtElectronEnDown*lep1->Pt()*lepscaleDOWNshift*(1 - cos(deltaPhi(metType1PtElectronEnDown,lep1->Phi()))));

		float MTTauEnUp = sqrt(lep1->M2() + 2*metType1PtTauEnUp*lep1->Pt()*(1 - cos(deltaPhi(metType1PhiTauEnUp,lep1->Phi()))));
		float MTTauEnDown = sqrt(lep1->M2() + 2*metType1PtTauEnDown*lep1->Pt()*(1 - cos(deltaPhi(metType1PhiTauEnDown,lep1->Phi()))));

		float MTUnclusteredEnUp = sqrt(lep1->M2() + 2*metType1PtUnclusteredEnUp*lep1->Pt()*(1 - cos(deltaPhi(metType1PhiUnclusteredEnUp,lep1->Phi()))));
		float MTUnclusteredEnDown = sqrt(lep1->M2() + 2*metType1PtUnclusteredEnDown*lep1->Pt()*(1 - cos(deltaPhi(metType1PhiUnclusteredEnDown,lep1->Phi()))));

		float MTPhotonEnUp = sqrt(lep1->M2() + 2*metType1PtPhotonEnUp*lep1->Pt()*(1 - cos(deltaPhi(metType1PhiPhotonEnUp,lep1->Phi()))));
		float MTPhotonEnDown = sqrt(lep1->M2() + 2*metType1PtPhotonEnDown*lep1->Pt()*(1 - cos(deltaPhi(metType1PhiPhotonEnDown,lep1->Phi()))));

				
		(mcMT[tree.first])[cut]->Fill(mTLepMet, eventWeight);
		(mcMTLepScaleUp[tree.first])[cut]->Fill(mcMTLepUp, eventWeight);
		(mcMTLepScaleDown[tree.first])[cut]->Fill(mcMTLepDown, eventWeight);

		(mcMTJetResUp[tree.first])[cut]->Fill(MTJetResUp, eventWeight);
		(mcMTJetResDown[tree.first])[cut]->Fill(MTJetResDown, eventWeight);

		(mcMTJetEnUp[tree.first])[cut]->Fill(MTJetEnUp, eventWeight);
		(mcMTJetEnDown[tree.first])[cut]->Fill(MTJetEnDown, eventWeight);

		(mcMTMuonEnUp[tree.first])[cut]->Fill(MTMuonEnUp, eventWeight);
		(mcMTMuonEnDown[tree.first])[cut]->Fill(MTMuonEnDown, eventWeight);

		(mcMTElectronEnUp[tree.first])[cut]->Fill(MTElectronEnUp, eventWeight);
		(mcMTElectronEnDown[tree.first])[cut]->Fill(MTElectronEnDown, eventWeight);

		(mcMTTauEnUp[tree.first])[cut]->Fill(MTTauEnUp, eventWeight);
		(mcMTTauEnDown[tree.first])[cut]->Fill(MTTauEnDown, eventWeight);

		(mcMTUnclusteredEnUp[tree.first])[cut]->Fill(MTUnclusteredEnUp, eventWeight);
		(mcMTUnclusteredEnDown[tree.first])[cut]->Fill(MTUnclusteredEnDown, eventWeight);

		(mcMTPhotonEnUp[tree.first])[cut]->Fill(MTPhotonEnUp, eventWeight);
		(mcMTPhotonEnDown[tree.first])[cut]->Fill(MTPhotonEnDown, eventWeight);
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

    TLegend *legend = new TLegend(0.6, 0.54, 0.9, 0.84);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    
    for(uint cut = 0; cut < cutSequence.size(); cut++){
      vector<string> orderedtrees {"WJets"};
      for(auto &tree : orderedtrees){
	mcMT[tree][cut]->SetLineColor(1);
	mcMTJetResUp[tree][cut]->SetLineColor(4);
	mcMTJetResDown[tree][cut]->SetLineColor(5);
	mcMTJetEnUp[tree][cut]->SetLineColor(6);
	mcMTJetEnDown[tree][cut]->SetLineColor(7);
	mcMTMuonEnUp[tree][cut]->SetLineColor(8);
	mcMTMuonEnDown[tree][cut]->SetLineColor(9);
	mcMTElectronEnUp[tree][cut]->SetLineColor(10);
	mcMTElectronEnDown[tree][cut]->SetLineColor(11);
	mcMTTauEnUp[tree][cut]->SetLineColor(12);
	mcMTTauEnDown[tree][cut]->SetLineColor(13);
	mcMTUnclusteredEnUp[tree][cut]->SetLineColor(14);
	mcMTUnclusteredEnDown[tree][cut]->SetLineColor(14);
	mcMTPhotonEnUp[tree][cut]->SetLineColor(15);
	mcMTPhotonEnDown[tree][cut]->SetLineColor(16);
      }     
    }

    // TLegend *legend = new TLegend(0.5, 0.5, 0.9, 0.9);
    legend->AddEntry(mcMT["WJets"][0], "W+Jets MC");
    legend->AddEntry(mcMTLepScaleUp["WJets"][0], mcMTLepScaleUp["WJets"][0]->GetName());
    legend->AddEntry(mcMTLepScaleDown["WJets"][0], mcMTLepScaleDown["WJets"][0]->GetName());
    legend->AddEntry(mcMTJetResUp["WJets"][0], mcMTJetResUp["WJets"][0]->GetName());
    legend->AddEntry(mcMTJetResDown["WJets"][0], mcMTJetResDown["WJets"][0]->GetName());
    legend->AddEntry(mcMTJetEnUp["WJets"][0], mcMTJetEnUp["WJets"][0]->GetName());
    legend->AddEntry(mcMTJetEnDown["WJets"][0], mcMTJetEnDown["WJets"][0]->GetName());
    legend->AddEntry(mcMTMuonEnUp["WJets"][0], mcMTMuonEnUp["WJets"][0]->GetName());
    legend->AddEntry(mcMTMuonEnDown["WJets"][0], mcMTMuonEnDown["WJets"][0]->GetName());
    legend->AddEntry(mcMTElectronEnUp["WJets"][0], mcMTElectronEnUp["WJets"][0]->GetName());
    legend->AddEntry(mcMTElectronEnDown["WJets"][0], mcMTElectronEnDown["WJets"][0]->GetName());
    legend->AddEntry(mcMTTauEnUp["WJets"][0], mcMTTauEnUp["WJets"][0]->GetName());
    legend->AddEntry(mcMTTauEnDown["WJets"][0], mcMTTauEnDown["WJets"][0]->GetName());
    legend->AddEntry(mcMTUnclusteredEnUp["WJets"][0], mcMTUnclusteredEnUp["WJets"][0]->GetName());
    legend->AddEntry(mcMTUnclusteredEnDown["WJets"][0], mcMTUnclusteredEnDown["WJets"][0]->GetName());
    legend->AddEntry(mcMTPhotonEnUp["WJets"][0], mcMTPhotonEnUp["WJets"][0]->GetName());
    legend->AddEntry(mcMTPhotonEnDown["WJets"][0], mcMTPhotonEnDown["WJets"][0]->GetName());


    for(uint cut = 0; cut < cutSequence.size(); cut++){
      //create histogram stacks for MC
      mcMT["WJets"][cut]->Draw();
      mcMTLepScaleUp["WJets"][cut]->Draw("same");
      mcMTLepScaleDown["WJets"][cut]->Draw("same");
      mcMTJetResUp["WJets"][cut]->Draw("same");
      mcMTJetResDown["WJets"][cut]->Draw("same");
      mcMTJetEnUp["WJets"][cut]->Draw("same");
      mcMTJetEnDown["WJets"][cut]->Draw("same");
      mcMTMuonEnUp["WJets"][cut]->Draw("same");
      mcMTMuonEnDown["WJets"][cut]->Draw("same");
      mcMTElectronEnUp["WJets"][cut]->Draw("same");
      mcMTElectronEnDown["WJets"][cut]->Draw("same");
      mcMTTauEnUp["WJets"][cut]->Draw("same");
      mcMTTauEnDown["WJets"][cut]->Draw("same");
      mcMTUnclusteredEnUp["WJets"][cut]->Draw("same");
      mcMTUnclusteredEnDown["WJets"][cut]->Draw("same");
      mcMTPhotonEnUp["WJets"][cut]->Draw("same");
      mcMTPhotonEnDown["WJets"][cut]->Draw("same");

      float AllMT = mcMT["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float JetResUpMT = mcMTJetResUp["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float JetResDownMT = mcMTJetResDown["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float JetEnUpMT = mcMTJetEnUp["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float JetEnDownMT = mcMTJetEnDown["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float MuonEnUpMT = mcMTMuonEnUp["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float MuonEnDownMT = mcMTMuonEnDown["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float ElectronEnUpMT = mcMTElectronEnUp["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float ElectronEnDownMT = mcMTElectronEnDown["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float TauEnUpMT = mcMTTauEnUp["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float TauEnDownMT = mcMTTauEnDown["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float UnclusteredEnUpMT = mcMTUnclusteredEnUp["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float UnclusteredEnDownMT = mcMTUnclusteredEnDown["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float PhotonEnUpMT = mcMTPhotonEnUp["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);
      float PhotonEnDownMT = mcMTPhotonEnDown["WJets"][cut]->Integral(mcMT["WJets"][cut]->GetXaxis()->FindFixBin(100.), mcMT["WJets"][cut]->GetNbinsX()+1);


      std::cout<<"Event yields: "<<mcMT["WJets"][cut]->GetName()<<" "<<AllMT<<" "<<endl;
      std::cout<<"Event yields in : "<<mcMTJetResUp["WJets"][cut]->GetName()<<" "<<(JetResUpMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTJetResDown["WJets"][cut]->GetName()<<" "<<(JetResDownMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTJetEnUp["WJets"][cut]->GetName()<<" "<<(JetEnUpMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTJetEnDown["WJets"][cut]->GetName()<<" "<<(JetEnDownMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTMuonEnUp["WJets"][cut]->GetName()<<" "<<(MuonEnUpMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTMuonEnDown["WJets"][cut]->GetName()<<" "<<(MuonEnDownMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTElectronEnUp["WJets"][cut]->GetName()<<" "<<(ElectronEnUpMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTElectronEnDown["WJets"][cut]->GetName()<<" "<<(ElectronEnDownMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTTauEnUp["WJets"][cut]->GetName()<<" "<<(TauEnUpMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTTauEnDown["WJets"][cut]->GetName()<<" "<<(TauEnDownMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTUnclusteredEnUp["WJets"][cut]->GetName()<<" "<<(UnclusteredEnUpMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTUnclusteredEnDown["WJets"][cut]->GetName()<<" "<<(UnclusteredEnDownMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTPhotonEnUp["WJets"][cut]->GetName()<<" "<<(PhotonEnUpMT-AllMT)/AllMT*100<<" %"<<endl;
      std::cout<<"Event yields in : "<<mcMTPhotonEnDown["WJets"][cut]->GetName()<<" "<<(PhotonEnDownMT-AllMT)/AllMT*100<<" %"<<endl;

      float JESunc = TMath::Sqrt(fabs((JetEnUpMT-AllMT)/AllMT) * fabs( (JetEnDownMT-AllMT)/AllMT ));
      float MuonEnunc = TMath::Sqrt(fabs((MuonEnUpMT-AllMT)/AllMT) * fabs( (MuonEnDownMT-AllMT)/AllMT ));
      float UnclEnunc = TMath::Sqrt(fabs((UnclusteredEnUpMT-AllMT)/AllMT) * fabs( (UnclusteredEnDownMT-AllMT)/AllMT ));
      float EleEnunc  = (ElectronEnUpMT-AllMT)/AllMT;
      cout<<"JES "<<JESunc*100<<endl;
      cout<<"MuonEn "<<MuonEnunc*100<<endl;
      cout<<"ElectronEn "<<(ElectronEnUpMT-AllMT)/AllMT*100<<endl;
      cout<<"UnclusteredEn "<<UnclEnunc*100<<endl;

      cout<<"Muon events Total uncertainty: "<<TMath::Sqrt(JESunc*JESunc+MuonEnunc*MuonEnunc+UnclEnunc*UnclEnunc)*100<<endl;
      cout<<"Ele events Total uncertainty: " <<TMath::Sqrt(JESunc*JESunc+EleEnunc*EleEnunc+UnclEnunc*UnclEnunc)*100<<endl;
 }

    legend->Draw();
    c.SaveAs("test.root");

    delete legend;
    for(uint cut = 0; cut < cutSequence.size(); cut++){
	for(auto &tree : mctrees){
	  delete mcMT[tree.first][cut];
	  delete mcMTLepScaleUp[tree.first][cut];
	  delete mcMTLepScaleDown[tree.first][cut];
	  delete mcMTJetResUp[tree.first][cut];
	  delete mcMTJetResDown[tree.first][cut];
	  delete mcMTJetEnUp[tree.first][cut];
	  delete mcMTJetEnDown[tree.first][cut];
	  delete mcMTMuonEnUp[tree.first][cut];
	  delete mcMTMuonEnDown[tree.first][cut];
	  delete mcMTElectronEnUp[tree.first][cut];
	  delete mcMTElectronEnDown[tree.first][cut];
	  delete mcMTTauEnUp[tree.first][cut];
	  delete mcMTTauEnDown[tree.first][cut];
	  delete mcMTUnclusteredEnUp[tree.first][cut];
	  delete mcMTUnclusteredEnDown[tree.first][cut];
	  delete mcMTPhotonEnUp[tree.first][cut];
	  delete mcMTPhotonEnDown[tree.first][cut];
        }
    }

    mcfiles["WJets"]->Close();
    pileupWeightFile->Close();
    delete mcfiles["WJets"];
    delete pileupWeightFile;
}

int main(){
    MTplots();
    return 0;
}

double deltaPhi(double phi1, double phi2) {
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
    dphi -= TMath::TwoPi();
  while (dphi <= -TMath::Pi())
    dphi += TMath::TwoPi();
  return dphi;
}
