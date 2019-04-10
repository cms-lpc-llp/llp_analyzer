#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

void makeDummyScaleFactorFiles(){

    //make dummy pileup reweighting histogram
    TFile *pileupFile = new TFile("DummyRun2PileupWeights.root", "recreate");
    TH1F *pileupHist = new TH1F("PUWeight_Run2", "PUWeight_Run2", 1, 0, 100);
    pileupHist->SetBinContent(1, 1.0);
    pileupHist->Write();
    pileupFile->Close();

    //make dummy muon scale factor histogram
    TFile *muonFile = new TFile("DummyRun2MuonWeights.root", "recreate");
    TH2F *muonLooseHist = new TH2F("MuonWeight_Run2_Loose", "MuonWeight_Run2_Loose", 1, -3, 3, 1, 0, 200);
    TH2F *muonTightHist = new TH2F("MuonWeight_Run2_Tight", "MuonWeight_Run2_Tight", 1, -3, 3, 1, 0, 200);
    muonLooseHist->SetBinContent(1, 1, 1.0);
    muonTightHist->SetBinContent(1, 1, 1.0);
    muonLooseHist->SetBinError(1, 1, .01);
    muonTightHist->SetBinError(1, 1, .01);
    muonLooseHist->Write();
    muonTightHist->Write();
    muonFile->Close();

    //make dummy electron scale factor histogram
    TFile *electronFile = new TFile("DummyRun2EleWeights.root", "recreate");
    TH2F *electronLooseHist = new TH2F("EleWeight_Run2_Loose", "EleWeight_Run2_Loose", 1, -3, 3, 1, 0, 200);
    TH2F *electronTightHist = new TH2F("EleWeight_Run2_Tight", "EleWeight_Run2_Tight", 1, -3, 3, 1, 0, 200);
    electronLooseHist->SetBinContent(1, 1, 1.0);
    electronTightHist->SetBinContent(1, 1, 1.0);
    electronLooseHist->SetBinError(1, 1, .01);
    electronTightHist->SetBinError(1, 1, .01);
    electronLooseHist->Write();
    electronTightHist->Write();
    electronFile->Close();

    //make dummy ttjets scale factor histogram
    TFile *ttjetsFile = new TFile("DummyRun2TTJetsSF.root", "recreate");
    TH2F *ttjetsHist = new TH2F("TTJetsSingleLepton", "TTJetsSingleLepton", 4, 0., 4000, 1, 0., 1.5);
    TH2F *ttjetsHistUp = new TH2F("TTJetsSingleLeptonUp", "TTJetsSingleLeptonUp", 4, 0., 4000, 1, 0., 1.5);
    TH2F *ttjetsHistDown = new TH2F("TTJetsSingleLeptonDown", "TTJetsSingleLeptonDown", 4, 0., 4000, 1, 0., 1.5);
    ttjetsHist->SetBinContent(1, 1, 1.0);
    ttjetsHist->SetBinError(1, 1, .1);
    ttjetsHist->SetBinContent(2, 1, 1.0);
    ttjetsHist->SetBinError(2, 1, .3);
    ttjetsHist->SetBinContent(3, 1, 1.0);
    ttjetsHist->SetBinError(3, 1, .1);
    ttjetsHist->SetBinContent(4, 1, 1.0);
    ttjetsHist->SetBinError(4, 1, .3);
    ttjetsHistUp->SetBinContent(1, 1, 1.05);
    ttjetsHistDown->SetBinContent(1, 1, 0.95);
    ttjetsHistUp->SetBinContent(2, 1, 1.05);
    ttjetsHistDown->SetBinContent(2, 1, 0.95);
    ttjetsHistUp->SetBinContent(3, 1, 1.05);
    ttjetsHistDown->SetBinContent(3, 1, 0.95);
    ttjetsHistUp->SetBinContent(4, 1, 1.05);
    ttjetsHistDown->SetBinContent(4, 1, 0.95);
    ttjetsHist->Write();
    ttjetsHistUp->Write();
    ttjetsHistDown->Write();
    ttjetsFile->Close();

    //make dummy wjets scale factor histogram
    TFile *wjetsFile = new TFile("DummyRun2WJetsSF.root", "recreate");
    TH2F *wjetsHist = new TH2F("WJetsSingleLepton", "WJetsSingleLepton", 1, 0., 4000, 1, 0., 1.5);
    TH2F *wjetsHistUp = new TH2F("WJetsSingleLeptonUp", "WJetsSingleLeptonUp", 1, 0., 4000, 1, 0., 1.5);
    TH2F *wjetsHistDown = new TH2F("WJetsSingleLeptonDown", "WJetsSingleLeptonDown", 1, 0., 4000, 1, 0., 1.5);
    wjetsHist->SetBinContent(1, 1, 1.0);
    wjetsHist->SetBinError(1, 1, .1);
    wjetsHistUp->SetBinContent(1, 1, 1.05);
    wjetsHistDown->SetBinContent(1, 1, 0.95);
    wjetsHist->Write();
    wjetsHistUp->Write();
    wjetsHistDown->Write();
    wjetsFile->Close();
    
    //make dummy dyjets scale factor histogram
    TFile *dyjetsFile = new TFile("DummyRun2DYJetsSF.root", "recreate");
    TH2F *dyjetsHist = new TH2F("DYJetsDilepton", "DYJetsDilepton", 1, 0., 4000, 1, 0., 1.5);
    TH2F *dyjetsHistUp = new TH2F("DYJetsDileptonUp", "DYJetsDileptonUp", 1, 0., 4000, 1, 0., 1.5);
    TH2F *dyjetsHistDown = new TH2F("DYJetsDileptonDown", "DYJetsDileptonDown", 1, 0., 4000, 1, 0., 1.5);
    dyjetsHist->SetBinContent(1, 1, 1.0);
    dyjetsHist->SetBinError(1, 1, .1);
    dyjetsHistUp->SetBinContent(1, 1, 1.05);
    dyjetsHistDown->SetBinContent(1, 1, .95);
    dyjetsHist->Write();
    dyjetsHistUp->Write();
    dyjetsHistDown->Write();
    dyjetsFile->Close();

    //make dummy zjets scale factor histogram
    TFile *zjetsFile = new TFile("DummyRun2ZNuNuSF.root", "recreate");
    TH2F *zjetsHist = new TH2F("ZNuNuGJets", "ZNuNuGJets", 1, 0., 4000, 1, 0., 1.5);
    TH2F *zjetsHistUp = new TH2F("ZNuNuGJetsUp", "ZNuNuGJetsUp", 1, 0., 4000, 1, 0., 1.5);
    TH2F *zjetsHistDown = new TH2F("ZNuNuGJetsDown", "ZNuNuGJetsDown", 1, 0., 4000, 1, 0., 1.5);
    TH2F *zjetsHistDY = new TH2F("ZNuNuDilepton", "ZNuNuDilepton", 1, 0., 4000, 1, 0., 1.5);
    zjetsHist->SetBinContent(1, 1, 1.0);
    zjetsHist->SetBinError(1, 1, .1);
    zjetsHistDY->SetBinContent(1, 1, 0.9);
    zjetsHistUp->SetBinContent(1, 1, 1.05);
    zjetsHistDown->SetBinContent(1, 1, .95);
    zjetsHist->Write();
    zjetsHistUp->Write();
    zjetsHistDown->Write();
    zjetsHistDY->Write();
    zjetsFile->Close();
}
