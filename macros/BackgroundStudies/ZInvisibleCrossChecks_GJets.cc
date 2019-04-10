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

    //set color palette 
    const Int_t NCont = 100;
    gStyle->SetNumberContours(NCont);
    gStyle->SetPaintTextFormat("1.0f");

    //upper bounds of reweighing histograms
    float maxPhotonPt = 999; 

    map<string, string> suffixes;
    suffixes["GJets"] = "_noPho";
    suffixes["EMQCD"] = "_noPho";
    suffixes["TTG"] = "_noPho";
    suffixes["WG"] = "_noPho";
    suffixes["ZG"] = "_noPho";

    //get input files -- assumes one TFile for each process, with weights for different HT bins 
    map<string, TFile*> mcfiles;
    map<string, TFile*> datafiles;
    mcfiles["GJets"] = new TFile("./GJetsRun1_19700pb_weighted.root");
    mcfiles["EMQCD"] = new TFile("./QCDRun1_19700pb_weighted.root");
    mcfiles["TTG"] = new TFile("./TTGJetsRun1_19700pb_weighted.root");
    mcfiles["WG"] = new TFile("./WGJetsRun1_19700pb_weighted.root");
    mcfiles["ZG"] = new TFile("./ZGJetsRun1_19700pb_weighted.root");

    datafiles["GJets"] = new TFile("./PhotonRun1_goodlumi.root");
    //get trees and set branches
    map<string, TTree*> mctrees;
    map<string, TTree*> datatrees;
    float weight;
    float leadingPhotonPt, leadingPhotonEta;
    bool hlt_photon;
    bool passedHLTPhoton50, passedHLTPhoton75, passedHLTPhoton90, passedHLTPhoton135, passedHLTPhoton150;
    int nPU_mean, nVtx, nSelectedPhotons;
    float MR_noPho, Rsq_noPho, HT_noPho, met_noPho, deltaPhi_noPho;
    int numJets_noPho, numJets80_noPho;
    bool Flag_HBHENoiseFilter, Flag_CSCTightHaloFilter, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_eeBadScFilter, Flag_ecalLaserCorrFilter;
      
    for(auto &file : mcfiles){
        mctrees[file.first] = (TTree*)file.second->Get("RazorInclusive");

        mctrees[file.first]->SetBranchStatus("*", 0); // disable all
        mctrees[file.first]->SetBranchStatus("weight", 1);
        mctrees[file.first]->SetBranchStatus("leadingPhotonPt", 1);
        mctrees[file.first]->SetBranchStatus("leadingPhotonEta", 1);
        mctrees[file.first]->SetBranchStatus("hlt_photon", 1);
        mctrees[file.first]->SetBranchStatus("nPU_mean", 1);
        mctrees[file.first]->SetBranchStatus("passedHLTPhoton50", 1);
        mctrees[file.first]->SetBranchStatus("passedHLTPhoton75", 1);
        mctrees[file.first]->SetBranchStatus("passedHLTPhoton90", 1);
        mctrees[file.first]->SetBranchStatus("passedHLTPhoton135", 1);
        mctrees[file.first]->SetBranchStatus("passedHLTPhoton150", 1);
        mctrees[file.first]->SetBranchStatus("numJets_noPho", 1);
        mctrees[file.first]->SetBranchStatus("numJets80_noPho", 1);
        mctrees[file.first]->SetBranchStatus("nVtx", 1); // enable 
        mctrees[file.first]->SetBranchStatus("MR_noPho", 1);
        mctrees[file.first]->SetBranchStatus("Rsq_noPho", 1);
        mctrees[file.first]->SetBranchStatus("deltaPhi_noPho", 1);
        mctrees[file.first]->SetBranchStatus("nSelectedPhotons", 1); // enable 
        mctrees[file.first]->SetBranchStatus("met_noPho", 1); // enable 
        mctrees[file.first]->SetBranchStatus("HT_noPho", 1); // enable 

        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("leadingPhotonPt", &leadingPhotonPt);
        mctrees[file.first]->SetBranchAddress("leadingPhotonEta", &leadingPhotonEta);
        mctrees[file.first]->SetBranchAddress("hlt_photon", &hlt_photon);
        mctrees[file.first]->SetBranchAddress("nPU_mean", &nPU_mean);
        mctrees[file.first]->SetBranchAddress("passedHLTPhoton50", &passedHLTPhoton50);
        mctrees[file.first]->SetBranchAddress("passedHLTPhoton75", &passedHLTPhoton75);
        mctrees[file.first]->SetBranchAddress("passedHLTPhoton90", &passedHLTPhoton90);
        mctrees[file.first]->SetBranchAddress("passedHLTPhoton135", &passedHLTPhoton135);
        mctrees[file.first]->SetBranchAddress("passedHLTPhoton150", &passedHLTPhoton150);
        mctrees[file.first]->SetBranchAddress("numJets_noPho", &numJets_noPho);
        mctrees[file.first]->SetBranchAddress("numJets80_noPho", &numJets80_noPho);
        mctrees[file.first]->SetBranchAddress("MR_noPho", &MR_noPho);
        mctrees[file.first]->SetBranchAddress("Rsq_noPho", &Rsq_noPho);
        mctrees[file.first]->SetBranchAddress("nSelectedPhotons", &nSelectedPhotons);
        mctrees[file.first]->SetBranchAddress("HT_noPho", &HT_noPho);
        mctrees[file.first]->SetBranchAddress("met_noPho", &met_noPho);
        mctrees[file.first]->SetBranchAddress("nVtx", &nVtx); // enable 
        mctrees[file.first]->SetBranchAddress("deltaPhi_noPho", &deltaPhi_noPho); // enable 
    }
    for(auto &file : datafiles){
        datatrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        
        datatrees[file.first]->SetBranchStatus("*", 0); // disable all
        datatrees[file.first]->SetBranchStatus("leadingPhotonPt", 1);
        datatrees[file.first]->SetBranchStatus("leadingPhotonEta", 1);
        datatrees[file.first]->SetBranchStatus("hlt_photon", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton50", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton75", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton90", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton135", 1);
        datatrees[file.first]->SetBranchStatus("passedHLTPhoton150", 1);
        datatrees[file.first]->SetBranchStatus("numJets_noPho", 1);
        datatrees[file.first]->SetBranchStatus("numJets80_noPho", 1);
        datatrees[file.first]->SetBranchStatus("MR_noPho", 1);
        datatrees[file.first]->SetBranchStatus("Rsq_noPho", 1);
        datatrees[file.first]->SetBranchStatus("deltaPhi_noPho", 1);
        datatrees[file.first]->SetBranchStatus("nSelectedPhotons", 1); // enable 
        datatrees[file.first]->SetBranchStatus("HT_noPho", 1); // enable 
        datatrees[file.first]->SetBranchStatus("met_noPho", 1); // enable 
        datatrees[file.first]->SetBranchStatus("nVtx", 1); // enable 
        datatrees[file.first]->SetBranchStatus("Flag_HBHENoiseFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_CSCTightHaloFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_eeBadScFilter", 1); // enable
        datatrees[file.first]->SetBranchStatus("Flag_ecalLaserCorrFilter", 1); // enable

        datatrees[file.first]->SetBranchAddress("leadingPhotonPt", &leadingPhotonPt);
        datatrees[file.first]->SetBranchAddress("leadingPhotonEta", &leadingPhotonEta);
        datatrees[file.first]->SetBranchAddress("hlt_photon", &hlt_photon);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton50", &passedHLTPhoton50);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton75", &passedHLTPhoton75);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton90", &passedHLTPhoton90);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton135", &passedHLTPhoton135);
        datatrees[file.first]->SetBranchAddress("passedHLTPhoton150", &passedHLTPhoton150);
        datatrees[file.first]->SetBranchAddress("numJets_noPho", &numJets_noPho);
        datatrees[file.first]->SetBranchAddress("numJets80_noPho", &numJets80_noPho);
        datatrees[file.first]->SetBranchAddress("MR_noPho", &MR_noPho);
        datatrees[file.first]->SetBranchAddress("Rsq_noPho", &Rsq_noPho);
        datatrees[file.first]->SetBranchAddress("nSelectedPhotons", &nSelectedPhotons);
        datatrees[file.first]->SetBranchAddress("HT_noPho", &HT_noPho);
        datatrees[file.first]->SetBranchAddress("met_noPho", &met_noPho);
        datatrees[file.first]->SetBranchAddress("nVtx", &nVtx); // enable 
        datatrees[file.first]->SetBranchAddress("deltaPhi_noPho", &deltaPhi_noPho); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter); // enable 
        datatrees[file.first]->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter); // enable 
    }

    //load efficiency/acceptance histograms
    TFile effFile("Run1LeptonPhotonEfficiency.root");
    TH2F photonEffHisto = *(TH2F *)effFile.Get("PhotonEfficiency");

    //load pileup reweighting histogram
    TFile *pileupWeightFile = new TFile("data/Run1PileupWeights.root", "READ");
    TH1F *pileupWeightHist = (TH1F*)pileupWeightFile->Get("PUWeight_Run1");
    assert(pileupWeightHist);
    
    //define cuts and histograms
    float lumi_HLTPhoton50  = 1.353e0 + 4.921e0 + 7.947e0 + 8.131e0;
    float lumi_HLTPhoton75  = 8.111e0 + 2.953e1 + 4.768e1 + 4.879e1;
    float lumi_HLTPhoton90  = 1.622e1 + 6.408e1 + 1.010e2 + 9.948e1;
    float lumi_HLTPhoton135 = 8.893e2 + 1.476e2 + 5.429e3 + 7.318e3;
    float lumi_HLTPhoton150 = 8.893e2 + 4.429e3 + 7.152e3 + 7.318e3;
    float nMRBins = 10;
    float nRsqBins = 8;
    float MRBinLowEdges[] = {300, 350, 400, 450, 550, 700, 900, 1200, 1600, 2500, 4000};
    float RsqBinLowEdges[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 0.80, 1.5};
    vector<string> cutSequence;
    vector<string> cutName;

    cutSequence.push_back( "hlt_photon && MR_noPho > 300 && Rsq_noPho > 0.15 && numJets80_noPho > 1 && leadingPhotonPt > 80" );
    // cutSequence.push_back( "leadingPhotonPt > 80 && hlt_photon" );
    cutName.push_back( "Photon Control Region" );

    // cutSequence.push_back( "leadingPhotonPt > 80 && hlt_photon && numJets_noPho > 1" );
    // cutName.push_back( "Require two 40-GeV Jets" );

    // cutSequence.push_back( "leadingPhotonPt > 80 && hlt_photon && HT_noPho > 160 " );
    // cutName.push_back( "Require HT_noPho > 160 GeV" );

    // cutSequence.push_back( "leadingPhotonPt > 80 && hlt_photon && numJets_noPho > 1 && deltaPhi_noPho < 2.7" );
    // cutName.push_back( "Require two 40-GeV Jets and #Delta #phi < 2.7" );

    // cutSequence.push_back( "leadingPhotonPt > 80 && hlt_photon && deltaPhi_noPho < 2.7 && numJets80_noPho > 1" );
    // cutName.push_back( "Require #Delta #phi < 2.7 and two 80-GeV Jets" );

    // cutSequence.push_back( "leadingPhotonPt > 80 && hlt_photon && deltaPhi_noPho < 2.7 && numJets80_noPho > 1 && MR_noPho > 300 && Rsq_noPho > 0.15" );
    // cutName.push_back( "Require #Delta #phi < 2.7, two 80-GeV Jets, MR > 300 GeV, Rsq > 0.15" );

    // cutSequence.push_back( "met_noPho > 100 && hlt_photon" );
    // cutName.push_back( "MET > 100 GeV" );

    map<string, vector<TH1F *> > mcPhotonPt, mcPhotonEta, mcNJets, mcNJets80, mcMR, mcRsq,  mcMet, mcNvtx, mcnSelectedPhotons, mcHT_noPho, mcdeltaPhi_noPho;
    vector<TH1F *> dataPhotonPt, dataPhotonEta, dataNJets, dataNJets80, dataMR, dataRsq, dataMet, dataNvtx, datanSelectedPhotons, dataHT_noPho, datadeltaPhi_noPho;
    for(auto &tree : mctrees){
        mcPhotonPt[tree.first] = vector<TH1F *>();
        mcPhotonEta[tree.first] = vector<TH1F *>();
        mcNJets[tree.first] = vector<TH1F *>();
        mcNJets80[tree.first] = vector<TH1F *>();
        mcMR[tree.first] = vector<TH1F *>();
        mcRsq[tree.first] = vector<TH1F *>();
        mcMet[tree.first] = vector<TH1F *>();
        mcNvtx[tree.first] = vector<TH1F *>();
        mcnSelectedPhotons[tree.first] = vector<TH1F *>();
        mcHT_noPho[tree.first] = vector<TH1F *>();
        mcdeltaPhi_noPho[tree.first] = vector<TH1F *>();
    }
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        for(auto &tree : mctrees){
            mcPhotonPt[tree.first].push_back(new TH1F(Form("mcPhotonPt%s%d", tree.first.c_str(), cut), Form("%s; photon pt (GeV)", cutName[cut].c_str()), 50, 20., 1000));
            mcPhotonEta[tree.first].push_back(new TH1F(Form("mcPhotonEta%s%d", tree.first.c_str(), cut), Form("%s; photon Eta ", cutName[cut].c_str()), 50, -3, 3));
            mcNJets[tree.first].push_back(new TH1F(Form("mcNJets%s%d", tree.first.c_str(), cut), Form("%s; Number of jets 40 GeV", cutName[cut].c_str()), 10, 0, 10));
            mcNvtx[tree.first].push_back(new TH1F(Form("mcNvtx%s%d", tree.first.c_str(), cut), Form("%s; NVtx (GeV)", cutName[cut].c_str()), 50, 0, 50));
            mcNJets80[tree.first].push_back(new TH1F(Form("mcNJets80%s%d", tree.first.c_str(), cut), Form("%s; Number of jets 80 GeV", cutName[cut].c_str()), 10, 0, 10));
            mcMR[tree.first].push_back(new TH1F(Form("mcMR%s%d", tree.first.c_str(), cut), Form("%s; MR (GeV)", cutName[cut].c_str()), 20, 300, 4000));
            // mcMR[tree.first].push_back(new TH1F(Form("mcMR%s%d", tree.first.c_str(), cut), Form("%s; MR (GeV)", cutName[cut].c_str()), nMRBins, MRBinLowEdges));
            mcRsq[tree.first].push_back(new TH1F(Form("mcRsq%s%d", tree.first.c_str(), cut), Form("%s; Rsq (GeV)", cutName[cut].c_str()), nRsqBins, RsqBinLowEdges));
            mcMet[tree.first].push_back(new TH1F(Form("mcMet%s%d", tree.first.c_str(), cut), Form("%s; MET (GeV)", cutName[cut].c_str()), 200, 0, 1000));
            mcnSelectedPhotons[tree.first].push_back(new TH1F(Form("mcnSelectedPhotons%s%d", tree.first.c_str(), cut), Form("%s; NSelected Photons", cutName[cut].c_str()), 10, 0, 10));
            mcHT_noPho[tree.first].push_back(new TH1F(Form("mcHT_noPho%s%d", tree.first.c_str(), cut), Form("%s; HT NoPho", cutName[cut].c_str()), 50, 0, 1000));
            mcdeltaPhi_noPho[tree.first].push_back(new TH1F(Form("mcdeltaPhi_noPho%s%d", tree.first.c_str(), cut), Form("%s; Delta Phi", cutName[cut].c_str()), 50, 0, 3.2));

            mcPhotonPt[tree.first][cut]->Sumw2();
            mcNJets[tree.first][cut]->Sumw2();
            mcNJets80[tree.first][cut]->Sumw2();
            mcMR[tree.first][cut]->Sumw2();
            mcRsq[tree.first][cut]->Sumw2();
            mcPhotonEta[tree.first][cut]->Sumw2();
            mcMet[tree.first][cut]->Sumw2();
            mcNvtx[tree.first][cut]->Sumw2();
	    mcnSelectedPhotons[tree.first][cut]->Sumw2();
	    mcHT_noPho[tree.first][cut]->Sumw2();
	    mcdeltaPhi_noPho[tree.first][cut]->Sumw2();
	    
        }
        dataPhotonPt.push_back(new TH1F(Form("dataPhotonPt%d", cut), Form("%s; photon pt (GeV)", cutName[cut].c_str()), 50, 20., 1000));
        dataPhotonEta.push_back(new TH1F(Form("dataPhotonEta%d", cut), Form("%s; photon Eta", cutName[cut].c_str()), 50, -3., 3));
        dataNJets.push_back(new TH1F(Form("dataNJets%d", cut), Form("%s; Number of jets 40 GeV", cutName[cut].c_str()), 10, 0, 10));
        dataNJets80.push_back(new TH1F(Form("dataNJets80%d", cut), Form("%s; Number of jets 80 GeV", cutName[cut].c_str()), 10, 0, 10));
        // dataMR.push_back(new TH1F(Form("dataMR%d", cut), Form("%s; MR (GeV)", cutName[cut].c_str()), nMRBins, MRBinLowEdges));
        dataMR.push_back(new TH1F(Form("dataMR%d", cut), Form("%s; MR (GeV)", cutName[cut].c_str()), 20, 300, 4000));
        dataRsq.push_back(new TH1F(Form("dataRsq%d", cut), Form("%s; Rsq (GeV)", cutName[cut].c_str()), nRsqBins, RsqBinLowEdges));
        dataMet.push_back(new TH1F(Form("dataMet%d", cut), Form("%s; mcMet (GeV)", cutName[cut].c_str()), 200, 0, 1000));
        dataNvtx.push_back(new TH1F(Form("dataNvtx%d", cut), Form("%s; NVtx (GeV)", cutName[cut].c_str()), 50, 0, 50));
        datanSelectedPhotons.push_back(new TH1F(Form("datanSelectedPhotons%d", cut), Form("%s; NSelected Photons", cutName[cut].c_str()), 10, 0, 10));
        dataHT_noPho.push_back(new TH1F(Form("dataHT_noPho%d", cut), Form("%s; HT NoPho", cutName[cut].c_str()), 50, 0, 1000));
        datadeltaPhi_noPho.push_back(new TH1F(Form("datadeltaPhi_noPho%d", cut), Form("%s; Delta Phi", cutName[cut].c_str()), 50, 0, 3.2));

        dataPhotonPt[cut]->Sumw2();
        dataPhotonEta[cut]->Sumw2();
        dataNJets[cut]->Sumw2();
        dataNJets80[cut]->Sumw2();
        dataMR[cut]->Sumw2();
        dataRsq[cut]->Sumw2();
        dataMet[cut]->Sumw2();
	dataNvtx[cut]->Sumw2();
	datanSelectedPhotons[cut]->Sumw2();
	dataHT_noPho[cut]->Sumw2();
	datadeltaPhi_noPho[cut]->Sumw2();
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
	      double effFactor = photonEffHisto.GetBinContent(photonEffHisto.FindFixBin(min(leadingPhotonPt, maxPhotonPt), fabs(leadingPhotonEta)));
	      
	      if(effFactor > 1e-5) eventWeight /= effFactor;
	      else{ 
		eventWeight = 0;
		//cout << "Warning: efficiency histogram gives 0 (pt " << leadingPhotonPt << ", eta " << leadingPhotonEta << "); setting event weight to 0" << endl;
	      }
	    }
            //apply selection cuts and fill the appropriate histograms
            for(uint cut = 0; cut < cutSequence.size(); cut++){
                bool passesCut = cuts[cut]->EvalInstance();
                if(!passesCut) continue;

                //check that the MC event passes the correct trigger
                bool passesCorrectTrigger = false;
                if(leadingPhotonPt > 165 && passedHLTPhoton150){
                    passesCorrectTrigger = true;
                }
                else if(leadingPhotonPt > 150 && passedHLTPhoton135){
                    passesCorrectTrigger = true;
                }
                else if(leadingPhotonPt > 100 && passedHLTPhoton90){
                    passesCorrectTrigger = true;
                }
                else if(leadingPhotonPt > 90 && passedHLTPhoton75){
                    passesCorrectTrigger = true;
                }
                else if(passedHLTPhoton50){
                    passesCorrectTrigger = true;
                }
                if(!passesCorrectTrigger) continue;

                (mcPhotonPt[tree.first])[cut]->Fill(leadingPhotonPt, eventWeight);
                (mcPhotonEta[tree.first])[cut]->Fill(leadingPhotonEta, eventWeight);
                (mcNJets[tree.first])[cut]->Fill(numJets_noPho, eventWeight);
                (mcNJets80[tree.first])[cut]->Fill(numJets80_noPho, eventWeight);
                (mcMR[tree.first])[cut]->Fill(MR_noPho, eventWeight);
                (mcRsq[tree.first])[cut]->Fill(Rsq_noPho, eventWeight);
                (mcMet[tree.first])[cut]->Fill(met_noPho, eventWeight);
                (mcNvtx[tree.first])[cut]->Fill(nVtx, eventWeight);
                (mcnSelectedPhotons[tree.first])[cut]->Fill(nSelectedPhotons, eventWeight);
                (mcHT_noPho[tree.first])[cut]->Fill(HT_noPho, eventWeight);
                (mcdeltaPhi_noPho[tree.first])[cut]->Fill(deltaPhi_noPho, eventWeight);
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
	    if(leadingPhotonPt>5000) continue;
	    
            //reweigh according to selection efficiency and acceptance
	    if(!computeDataOverMCSFs){
	      double effFactor = photonEffHisto.GetBinContent(photonEffHisto.FindBin(min(leadingPhotonPt, maxPhotonPt), fabs(leadingPhotonEta)));

	      if(effFactor > 1e-5) eventWeight /= effFactor;
	      else{ 
		eventWeight = 0;
	      }
	    }
            //get weight if associate each trigger with a particular pt range
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
            eventWeight *= triggerWeightRestricted;

            //apply selection cuts and fill the appropriate histograms
            for(uint cut = 0; cut < cutSequence.size(); cut++){
                bool passesCut = cuts[cut]->EvalInstance();
                if(!passesCut) continue;

                dataPhotonPt[cut]->Fill(leadingPhotonPt, eventWeight);
		dataPhotonEta[cut]->Fill(leadingPhotonEta, eventWeight);
		dataNJets[cut]->Fill(numJets_noPho, eventWeight);
                dataNJets80[cut]->Fill(numJets80_noPho, eventWeight);
                dataMR[cut]->Fill(MR_noPho, eventWeight);
                dataRsq[cut]->Fill(Rsq_noPho, eventWeight);
                dataMet[cut]->Fill(met_noPho, eventWeight);
		dataNvtx[cut]->Fill(nVtx, eventWeight);
		datanSelectedPhotons[cut]->Fill(nSelectedPhotons, eventWeight);
 		dataHT_noPho[cut]->Fill(HT_noPho, eventWeight);
		datadeltaPhi_noPho[cut]->Fill(deltaPhi_noPho, eventWeight);
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
    colors["GJets"] = 9;
    colors["EMQCD"] = 8;
    colors["TTG"] = 7;
    colors["WG"] = 38;
    colors["ZG"] = 4;
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(dataPhotonPt[0], "2012 Data, Photon CS");
    legend->AddEntry(mcPhotonPt["GJets"][0], "Gamma+Jets MC");
    legend->AddEntry(mcPhotonPt["EMQCD"][0], "QCD MC");
    legend->AddEntry(mcPhotonPt["TTG"][0], "TT+Gamma MC");
    legend->AddEntry(mcPhotonPt["WG"][0], "W+Gamma MC");
    legend->AddEntry(mcPhotonPt["ZG"][0], "Z+Gamma MC");
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        //create histogram stacks for MC
        THStack PhotonPtMC(Form("PhotonPtStack%d", cut), cutName[cut].c_str());
        THStack PhotonEtaMC(Form("PhotonEtaStack%d", cut), cutName[cut].c_str());
        THStack NumJetsMC(Form("NumJetsStack%d", cut), cutName[cut].c_str());
        THStack NumJets80MC(Form("NumJets80Stack%d", cut), cutName[cut].c_str());
        THStack MRMC(Form("MRStack%d", cut), cutName[cut].c_str());
        THStack RsqMC(Form("RsqStack%d", cut), cutName[cut].c_str());
        THStack MetMC(Form("MetStack%d", cut), cutName[cut].c_str());
        THStack NVtxMC(Form("NVtxStack%d", cut), cutName[cut].c_str());
        THStack NPhotonsMC(Form("NPhotonStack%d", cut), cutName[cut].c_str());
        THStack HTMC(Form("HTStack%d", cut), cutName[cut].c_str());
        THStack DPhiMC(Form("DPhiMC%d", cut), cutName[cut].c_str());
        
        //add the histograms to the stack in order
        vector<string> orderedtrees {"ZG", "WG", "TTG", "EMQCD", "GJets"};
        for(auto &tree : orderedtrees){
            mcPhotonPt[tree][cut]->SetFillColor(colors[tree]);
	    mcPhotonEta[tree][cut]->SetFillColor(colors[tree]);
	    mcNJets[tree][cut]->SetFillColor(colors[tree]);
            mcNJets80[tree][cut]->SetFillColor(colors[tree]);
            mcMR[tree][cut]->SetFillColor(colors[tree]);
            mcRsq[tree][cut]->SetFillColor(colors[tree]);
            mcMet[tree][cut]->SetFillColor(colors[tree]);
            mcNvtx[tree][cut]->SetFillColor(colors[tree]);
	    mcnSelectedPhotons[tree][cut]->SetFillColor(colors[tree]);
	    mcHT_noPho[tree][cut]->SetFillColor(colors[tree]);
	    mcdeltaPhi_noPho[tree][cut]->SetFillColor(colors[tree]);

            PhotonPtMC.Add(mcPhotonPt[tree][cut]);
            PhotonEtaMC.Add(mcPhotonEta[tree][cut]);
            NumJetsMC.Add(mcNJets[tree][cut]);
            NumJets80MC.Add(mcNJets80[tree][cut]);
            MRMC.Add(mcMR[tree][cut]);
            RsqMC.Add(mcRsq[tree][cut]);
	    MetMC.Add(mcMet[tree][cut]);
	    NVtxMC.Add(mcNvtx[tree][cut]);
	    NPhotonsMC.Add(mcnSelectedPhotons[tree][cut]);
	    HTMC.Add(mcHT_noPho[tree][cut]);
	    DPhiMC.Add(mcdeltaPhi_noPho[tree][cut]);
        }
        DrawDataVsMCRatioPlot(dataPhotonPt[cut], &PhotonPtMC, legend, "Photon pt (GeV)", "MediumPhotonCrossChecks_newQCD_PhotonPt"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataPhotonEta[cut], &PhotonEtaMC, legend, "Photon Eta", "MediumPhotonCrossChecks_newQCD_PhotonEta"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataNJets[cut], &NumJetsMC, legend, "Number of jets 40 GeV", "MediumPhotonCrossChecks_newQCD_NumJets"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataNJets80[cut], &NumJets80MC, legend, "Number of jets 80 GeV", "MediumPhotonCrossChecks_newQCD_NumJets80"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataMR[cut], &MRMC, legend, "MR (GeV)", "MediumPhotonCrossChecks_newQCD_MR"+to_string(cut), true);
        DrawDataVsMCRatioPlot(dataRsq[cut], &RsqMC, legend, "Rsq", "MediumPhotonCrossChecks_newQCD_Rsq"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataMet[cut], &MetMC, legend, "Met", "MediumPhotonCrossChecks_newQCD_MET"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataNvtx[cut], &NVtxMC, legend, "NVtx", "MediumPhotonCrossChecks_newQCD_NVtx"+to_string(cut), false);
        DrawDataVsMCRatioPlot(datanSelectedPhotons[cut], &NPhotonsMC, legend, "nPhotons", "MediumPhotonCrossChecksN_newQCD_Pho"+to_string(cut), false);
        DrawDataVsMCRatioPlot(dataHT_noPho[cut], &HTMC, legend, "HT", "MediumPhotonCrossChecks_newQCD_HT"+to_string(cut), false);
        DrawDataVsMCRatioPlot(datadeltaPhi_noPho[cut], &DPhiMC, legend, "DPhi", "MediumPhotonCrossChecks_newQCD_DPhi"+to_string(cut), false);
    }

    delete legend;
    for(uint cut = 0; cut < cutSequence.size(); cut++){
        delete dataPhotonPt[cut];
        delete dataPhotonEta[cut];
        delete dataNJets[cut];
        delete dataNJets80[cut];
        delete dataMR[cut];
        delete dataRsq[cut];
 	delete dataMet[cut];
	delete dataNvtx[cut];
	delete datanSelectedPhotons[cut];
	delete dataHT_noPho[cut];
	delete datadeltaPhi_noPho[cut];
	for(auto &tree : mctrees){
	  delete mcPhotonPt[tree.first][cut];
	  delete mcPhotonEta[tree.first][cut];
	  delete mcNJets[tree.first][cut];
	  delete mcNJets80[tree.first][cut];
	  delete mcMR[tree.first][cut];
	  delete mcRsq[tree.first][cut];
	  delete mcMet[tree.first][cut];
	  delete mcNvtx[tree.first][cut];
	  delete mcnSelectedPhotons[tree.first][cut];
	  delete mcdeltaPhi_noPho[tree.first][cut];
        }
    }

    mcfiles["GJets"]->Close();
    mcfiles["EMQCD"]->Close();
    mcfiles["TTG"]->Close();
    mcfiles["WG"]->Close();
    mcfiles["ZG"]->Close();
    datafiles["GJets"]->Close();
    pileupWeightFile->Close();
    delete mcfiles["GJets"];
    delete mcfiles["EMQCD"];
    delete mcfiles["TTG"];
    delete mcfiles["WG"];
    delete mcfiles["ZG"];
    delete datafiles["GJets"];
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
