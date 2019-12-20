

//================================================================================================
//
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TChain.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TH1D.h>                
#include <TH1F.h>                
#include <TH2F.h>                
#include <TCanvas.h>                
#include <TLegend.h> 
#include <TPaveText.h>
#include <THStack.h> 
#include <TFractionFitter.h> 


#endif


void PhotonControlSample_MakeTemplatesPrompt(int option) {

  // The template for prompts in barrel was made with this snippet on command line
  
  TChain * chain = new TChain("ControlSampleEvent");
  chain->Add("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDCuts/RunTwoRazorControlRegions_PhotonFull_GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  chain->Add("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDCuts/RunTwoRazorControlRegions_PhotonFull_GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  chain->Add("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDCuts/RunTwoRazorControlRegions_PhotonFull_GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  chain->Add("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDCuts/RunTwoRazorControlRegions_PhotonFull_GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  chain->Add("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDCuts/RunTwoRazorControlRegions_PhotonFull_GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");


  double sigmaietaieta_bins_EB[100] = {0.};
  double sigmaietaieta_bins_EE[100] = {0.};
  
  sigmaietaieta_bins_EB[21] = 0.0103;
  sigmaietaieta_bins_EE[55] = 0.0271;
  const int NBinsSigmaietaieta = sizeof(sigmaietaieta_bins_EB)/sizeof(double)-1;
  
  for(int a = 1; a < 100; a++) { if(a!=21) sigmaietaieta_bins_EB[a] = sigmaietaieta_bins_EB[a-1] + 0.0005; }
  for(int a = 1; a < 100; a++) { if(a!=55) sigmaietaieta_bins_EE[a] = sigmaietaieta_bins_EE[a-1] + 0.0005; }
  
  
  TH1F* hBarrel = new TH1F("PhotonSigmaIEtaIEtaTemplate_Prompt_Barrel","PhotonSigmaIEtaIEtaTemplate_Prompt_Barrel",NBinsSigmaietaieta,sigmaietaieta_bins_EB);
  hBarrel->Sumw2();
  TH1F* hEndcap = new TH1F("PhotonSigmaIEtaIEtaTemplate_Prompt_Endcap","PhotonSigmaIEtaIEtaTemplate_Prompt_Endcap",NBinsSigmaietaieta,sigmaietaieta_bins_EE);
  hEndcap->Sumw2();

  string label = "";

  
  if (option == 0) {
    label = "_MR300Rsq0p15";
    chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Prompt_Barrel","weight*(pho1.Pt() > 185 && MR_NoPho>300 && Rsq_NoPho>0.15 && fabs(pho1.Eta())<1.479&&pho1_chargediso<2.5&&pho1_sigmaietaieta<0.015)");
    chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Prompt_Endcap","weight*(pho1.Pt() > 185 && MR_NoPho>300 && Rsq_NoPho>0.15 && fabs(pho1.Eta())>1.479 && fabs(pho1.Eta()) < 2.5&&pho1_chargediso<2.5&&pho1_sigmaietaieta>0.015)");
  } 
  if (option == 5) {
    label = "_Inclusive";
  chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Prompt_Barrel","weight*(fabs(pho1.Eta())<1.479&&pho1_chargediso<2.5&&pho1_sigmaietaieta<0.015)");
  chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Prompt_Endcap","weight*(fabs(pho1.Eta())>1.479 && fabs(pho1.Eta()) < 2.5&&pho1_chargediso<2.5&&pho1_sigmaietaieta>0.015)");
  } 


  TFile *outputFile = new TFile ( Form("PhotonTemplatesPrompt%s.root",label.c_str()), "UPDATE");
  outputFile->WriteTObject( hBarrel, "PhotonSigmaIEtaIEtaTemplate_Prompt_Barrel", "WriteDelete");
  outputFile->WriteTObject( hEndcap, "PhotonSigmaIEtaIEtaTemplate_Prompt_Endcap", "WriteDelete");
  outputFile->Close();
    
}

void PhotonControlSample_MakeTemplatesFake(int option) {

  // The template for prompts in barrel was made with this snippet on command line
  
  TChain * chain = new TChain("ControlSampleEvent");
  chain->Add("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDIsoCuts/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_2016B_PRv2_GoodLumiGolden.root");
  chain->Add("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDIsoCuts/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_2016C_PRv2_GoodLumiGolden.root");
  chain->Add("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDIsoCuts/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_2016D_PRv2_GoodLumiGolden.root");
  chain->Add("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDIsoCuts/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_2016E_PRv2_GoodLumiGolden.root");
  chain->Add("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDIsoCuts/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_2016F_PRv1_GoodLumiGolden.root");
  chain->Add("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p6_25October2016_CustomType1MET/PhotonAddToMETNoIDIsoCuts/RunTwoRazorControlRegions_PhotonFull_SinglePhoton_2016G_PRv1_GoodLumiGolden.root");

 

  double sigmaietaieta_bins_EB[100] = {0.};
  double sigmaietaieta_bins_EE[100] = {0.};
  
  sigmaietaieta_bins_EB[21] = 0.0103;
  sigmaietaieta_bins_EE[55] = 0.0271;
  const int NBinsSigmaietaieta = sizeof(sigmaietaieta_bins_EB)/sizeof(double)-1;
  
  for(int a = 1; a < 100; a++) { if(a!=21) sigmaietaieta_bins_EB[a] = sigmaietaieta_bins_EB[a-1] + 0.0005; }
  for(int a = 1; a < 100; a++) { if(a!=55) sigmaietaieta_bins_EE[a] = sigmaietaieta_bins_EE[a-1] + 0.0005; }
  
  
  TH1F* hBarrel = new TH1F("PhotonSigmaIEtaIEtaTemplate_Fake_Barrel","PhotonSigmaIEtaIEtaTemplate_Fake_Barrel",NBinsSigmaietaieta,sigmaietaieta_bins_EB);
  hBarrel->Sumw2();
  TH1F* hEndcap = new TH1F("PhotonSigmaIEtaIEtaTemplate_Fake_Endcap","PhotonSigmaIEtaIEtaTemplate_Fake_Endcap",NBinsSigmaietaieta,sigmaietaieta_bins_EE);
  hEndcap->Sumw2();

  string label = "";

  //Use non-R9 triggers
  if (option == 0) {
    label = "_NoR9Triggers_Pt185MR300Rsq0p15";
    chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Barrel","weight*(pho1.Pt() > 185 && MR_NoPho>300 && Rsq_NoPho>0.15 && fabs(pho1.Eta())<1.479 &&pho1_chargediso>2.5&&pho1_sigmaietaieta<0.015&&HLTDecision[102]&& pho1HLTFilter[36])");
    chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Endcap","weight*(pho1.Pt() > 185 && MR_NoPho>300 && Rsq_NoPho>0.15 && fabs(pho1.Eta())>1.479 && fabs(pho1.Eta()) < 2.5&&pho1_chargediso>2.5&&pho1_sigmaietaieta>0.015&&HLTDecision[102]&& pho1HLTFilter[36])");
  } 
  if (option == 1) {
    label = "_NoR9Triggers_Pt185";
    chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Barrel","weight*(pho1.Pt() > 185 && fabs(pho1.Eta())<1.479 &&pho1_chargediso>2.5&&pho1_sigmaietaieta<0.015&&HLTDecision[102]&& pho1HLTFilter[36])");
    chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Endcap","weight*(pho1.Pt() > 185 && fabs(pho1.Eta())>1.479 && fabs(pho1.Eta()) < 2.5&&pho1_chargediso>2.5&&pho1_sigmaietaieta>0.015&&HLTDecision[102]&& pho1HLTFilter[36])");
  } 
  if (option == 2) {
    label = "_NoR9Triggers_Inclusive";
    chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Barrel","weight*(fabs(pho1.Eta())<1.479 &&pho1_chargediso>2.5&&pho1_sigmaietaieta<0.015&&HLTDecision[102]&& pho1HLTFilter[36])");
    chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Endcap","weight*(fabs(pho1.Eta())>1.479 && fabs(pho1.Eta()) < 2.5&&pho1_chargediso>2.5&&pho1_sigmaietaieta>0.015&&HLTDecision[102]&& pho1HLTFilter[36])");
  } 

  if (option == 3) {
    label = "_NoR9Triggers_Pt185MR700Rsq0p15";
    chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Barrel","weight*(pho1.Pt() > 185 && MR_NoPho>700 && Rsq_NoPho>0.15 && fabs(pho1.Eta())<1.479 &&pho1_chargediso>2.5&&pho1_sigmaietaieta<0.015&&HLTDecision[102]&& pho1HLTFilter[36])");
    chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Endcap","weight*(pho1.Pt() > 185 && MR_NoPho>700 && Rsq_NoPho>0.15 && fabs(pho1.Eta())>1.479 && fabs(pho1.Eta()) < 2.5&&pho1_chargediso>2.5&&pho1_sigmaietaieta>0.015&&HLTDecision[102]&& pho1HLTFilter[36])");
  } 
  //Use R9 triggers
  if (option == 10) {
    label = "_R9Triggers_Pt185MR300Rsq0p15";
     chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Barrel","weight*(pho1.Pt() > 185 && MR_NoPho>300 && Rsq_NoPho>0.15 && fabs(pho1.Eta())<1.479&&pho1_chargediso>2.5&&pho1_sigmaietaieta<0.015&&HLTDecision[112]&& pho1HLTFilter[27])");
     chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Endcap","weight*(pho1.Pt() > 185 && MR_NoPho>300 && Rsq_NoPho>0.15 && fabs(pho1.Eta())>1.479 && fabs(pho1.Eta()) < 2.5&&pho1_chargediso>2.5&&pho1_sigmaietaieta>0.015&&HLTDecision[112]&& pho1HLTFilter[27])");
  }
  if (option == 11) {
    label = "_R9Triggers_Pt185";
     chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Barrel","weight*(pho1.Pt() > 185 && fabs(pho1.Eta())<1.479&&pho1_chargediso>2.5&&pho1_sigmaietaieta<0.015&&HLTDecision[112]&& pho1HLTFilter[27])");
     chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Endcap","weight*(pho1.Pt() > 185 && fabs(pho1.Eta())>1.479 && fabs(pho1.Eta()) < 2.5&&pho1_chargediso>2.5&&pho1_sigmaietaieta>0.015&&HLTDecision[112]&& pho1HLTFilter[27])");
  }
  if (option == 12) {
    label = "_R9Triggers_Inclusive";
     chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Barrel","weight*(fabs(pho1.Eta())<1.479&&pho1_chargediso>2.5&&pho1_sigmaietaieta<0.015&&HLTDecision[112]&& pho1HLTFilter[27])");
     chain->Draw("pho1_sigmaietaieta>>PhotonSigmaIEtaIEtaTemplate_Fake_Endcap","weight*(fabs(pho1.Eta())>1.479 && fabs(pho1.Eta()) < 2.5&&pho1_chargediso>2.5&&pho1_sigmaietaieta>0.015&&HLTDecision[112]&& pho1HLTFilter[27])");
  }

  TFile *outputFile = new TFile ( Form("PhotonTemplates%s.root",label.c_str()), "UPDATE");
  outputFile->WriteTObject( hBarrel, "PhotonSigmaIEtaIEtaTemplate_Fake_Barrel", "WriteDelete");
  outputFile->WriteTObject( hEndcap, "PhotonSigmaIEtaIEtaTemplate_Fake_Endcap", "WriteDelete");
  outputFile->Close();
    
}

void PhotonControlSample_MakeTemplates(int option) {

  PhotonControlSample_MakeTemplatesPrompt(option);
  //PhotonControlSample_MakeTemplatesFake(option);

}
