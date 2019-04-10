//================================================================================================
//
// Simple Example
//
//root -l RazorAnalyzer/macros/ObjectStudies/MakeElectronEfficiencyPlots.C+'("/afs/cern.ch/work/s/sixie/public/Run2SUSY/ElectronNtuple/ElectronNtuple_PromptGenLevel_TTJets_25ns.root",-1,"Electron")'
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
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
#include <THStack.h> 

#include "RazorAnalyzer/include/ControlSampleEvents.h"
#include "RazorAnalyzer/macros/tdrstyle.C"
#include "RazorAnalyzer/macros/CMS_lumi.C"

#endif

void PlotDataAndStackedBkg( vector<TH1D*> hist , vector<string> processLabels, vector<int> color,  bool hasData, string varName, string label ) {

  TCanvas *cv =0;
  TLegend *legend = 0;

  cv = new TCanvas("cv","cv", 800,700);
  cv->SetHighLightColor(2);
  cv->SetFillColor(0);
  cv->SetBorderMode(0);
  cv->SetBorderSize(2);
  cv->SetLeftMargin(0.16);
  cv->SetRightMargin(0.3);
  cv->SetTopMargin(0.07);
  cv->SetBottomMargin(0.12);
  cv->SetFrameBorderMode(0);  

  TPad *pad1 = new TPad("pad1","pad1", 0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->SetRightMargin(0.04);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.60,0.50,0.90,0.84);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  THStack *stack = new THStack();
  TH1D *histDataOverMC = (TH1D*)hist[0]->Clone("histDataOverMC");

  if (hasData) {
    for (int i = hist.size()-1; i >= 1; --i) {
      hist[i]->SetFillColor(color[i]);
      hist[i]->SetFillStyle(1001);
      
      if ( hist[i]->Integral() > 0) {
  	stack->Add(hist[i]);
      }
    }
  } else {
    for (int i = hist.size()-1; i >= 0; --i) {
      hist[i]->SetFillColor(color[i]);
      hist[i]->SetFillStyle(1001);
      
      if ( hist[i]->Integral() > 0) {
  	stack->Add(hist[i]);
      }
    }
  }

  for (uint i = 0 ; i < hist.size(); ++i) {
    if (hasData && i==0) {
      legend->AddEntry(hist[i],(processLabels[i]).c_str(), "LP");
    } else {
      legend->AddEntry(hist[i],(processLabels[i]).c_str(), "F");
    }
  }

  if (stack->GetHists()->GetEntries() > 0) {
    stack->Draw("hist");
    stack->GetHistogram()->GetXaxis()->SetTitle(((TH1D*)(stack->GetHists()->At(0)))->GetXaxis()->GetTitle());
    stack->GetHistogram()->GetYaxis()->SetTitle(((TH1D*)(stack->GetHists()->At(0)))->GetYaxis()->GetTitle());    
    stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.0);
    stack->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    stack->GetHistogram()->GetXaxis()->SetTitleSize(0.15);
    stack->SetMaximum( 1.2* fmax( stack->GetMaximum(), hist[0]->GetMaximum()) );
    stack->SetMinimum( 0.1 );

    if (hasData) {
      hist[0]->SetMarkerStyle(20);      
      hist[0]->SetMarkerSize(1);
      hist[0]->SetLineWidth(1);
      hist[0]->SetLineColor(color[0]);
      hist[0]->Draw("pesame");
    }
    legend->Draw();
  }

  //****************************
  //Add CMS and Lumi Labels
  //****************************
  // lumi_13TeV = "42 pb^{-1}";
  lumi_13TeV = "35.9 fb^{-1}";
  writeExtraText = true;
  relPosX = 0.13;
  CMS_lumi(pad1,4,0);

  cv->cd();
  cv->Update();


  TPad *pad2 = new TPad("pad2","pad2", 0,0,1,0.25);
  pad2->SetTopMargin(0.01);
  pad2->SetBottomMargin(0.37);
  pad2->SetRightMargin(0.04);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
    
  for (int b=0; b<histDataOverMC->GetXaxis()->GetNbins()+2; ++b) {
    double data = 0;
    if (hasData) {
      data = hist[0]->GetBinContent(b);
    }
    double MC = 0;
    double MCErrSqr = 0;
    if (hasData) {
      for (uint i = 1 ; i < hist.size(); ++i) {
	MC += hist[i]->GetBinContent(b);
	MCErrSqr += pow(hist[i]->GetBinError(b),2);
      }
    } else {
      MC = 1;
    }
      
    if (MC > 0) {
      histDataOverMC->SetBinContent(b, data / MC);
      histDataOverMC->SetBinError(b, (data / MC)*sqrt(1/data + MCErrSqr/pow(MC,2) ));
    } else {
      histDataOverMC->SetBinContent(b, 0);
      histDataOverMC->SetBinError(b, 0);
    }
    //cout << "bin " << b << " : " << histDataOverMC->GetBinContent(b) << " " << histDataOverMC->GetBinError(b) << "\n";
  }

  histDataOverMC->GetYaxis()->SetTitle("Data/MC");
  histDataOverMC->GetYaxis()->SetNdivisions(306);
  histDataOverMC->GetYaxis()->SetTitleSize(0.10);
  histDataOverMC->GetYaxis()->SetTitleOffset(0.3);
  histDataOverMC->GetYaxis()->SetRangeUser(0.5,1.5);
  histDataOverMC->GetYaxis()->SetLabelSize(0.10);
  histDataOverMC->GetXaxis()->SetLabelSize(0.125);
  histDataOverMC->GetXaxis()->SetTitleSize(0.15);
  histDataOverMC->GetXaxis()->SetTitleOffset(1.0);
  histDataOverMC->SetLineColor(kBlack);
  histDataOverMC->SetMarkerStyle(20);      
  histDataOverMC->SetMarkerSize(1);
  histDataOverMC->SetStats(false);
  histDataOverMC->Draw("pe");
  
  pad1->SetLogy(false);
  cv->SaveAs(Form("WWZ_%s%s.png",varName.c_str(), label.c_str()));
  cv->SaveAs(Form("WWZ_%s%s.pdf",varName.c_str(), label.c_str()));
  
  pad1->SetLogy(true);
  cv->SaveAs(Form("WWZ_%s%s_Logy.png",varName.c_str(),label.c_str()));
  cv->SaveAs(Form("WWZ_%s%s_Logy.pdf",varName.c_str(),label.c_str()));


 

}




//=== MAIN MACRO ================================================================================================= 


void RunSelectWWZTo4L(  vector<string> datafiles, vector<vector<string> > bkgfiles, vector<string> bkgLabels, vector<int> bkgColors, double lumi, string option, int channelOption = -1, string label = "") {
  
  string Label = "";
  if (label != "") Label = "_" + label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  float MR = 0;
  float Rsq = 0;
  float mll = 0;


  bool printdebug = false;

  //*****************************************************************************************
  //Make some histograms
  //*****************************************************************************************

  vector<vector<string> > inputfiles;
  vector<string> processLabels;
  vector<int> color;

  inputfiles.push_back(datafiles);
  processLabels.push_back("Data");
  color.push_back(kBlack);
  
  assert(bkgfiles.size() == bkgLabels.size());
  assert(bkgfiles.size() == bkgColors.size());
  for (int i=0; i < int(bkgfiles.size()); ++i) {
    inputfiles.push_back(bkgfiles[i]);
    processLabels.push_back(bkgLabels[i]);
    color.push_back(bkgColors[i]);
  }

  vector<double> EventYield;
  vector<double> EventYieldErrSqr;
  vector<TH1D*> histZMass;
  vector<TH1D*> histZPt;
  vector<TH1D*> histLep3Pt;
  vector<TH1D*> histLep4Pt;
  vector<TH1D*> histMET;
  vector<TH1D*> histLep3MT;
  vector<TH1D*> histLep4MT;
  vector<TH1D*> histLep34MT;
  vector<TH1D*> histLep1234MT;
  vector<TH1D*> histMass4L;
  vector<TH1D*> histLep34Mass;
  vector<TH1D*> histNBJet30;
  vector<TH1D*> histNBJet20;
  vector<TH1D*> histNJet30;
  vector<TH1D*> histNJet20;
  vector<TH1D*> histPhi0;
  vector<TH1D*> histTheta0;
  vector<TH1D*> histPhi;
  vector<TH1D*> histTheta1;
  vector<TH1D*> histTheta2;

  assert (inputfiles.size() == processLabels.size());
  for (uint i=0; i < inputfiles.size(); ++i) {
    EventYield.push_back(0);
    EventYieldErrSqr.push_back(0);

    histZMass.push_back(new TH1D(Form("histZMass_%s",processLabels[i].c_str()), ";  Z Candidate Mass[GeV/c^{2}]; Number of Events", 30, 60, 120));   
    histZMass[i]->Sumw2();
  
    histZPt.push_back(new TH1D(Form("histZPt_%s",processLabels[i].c_str()), ";  [GeV/c^{2}]; Number of Events", 25, 0, 500));   
    histZPt[i]->Sumw2();
  
    histLep3Pt.push_back(new TH1D(Form("histLep3Pt_%s",processLabels[i].c_str()), "; Lepton 3 p_{T} [GeV/c]; Number of Events", 15, 0, 300));   
    histLep3Pt[i]->Sumw2();
  
    histLep4Pt.push_back(new TH1D(Form("histLep4Pt_%s",processLabels[i].c_str()), ";  [GeV/c^{2}]; Number of Events", 15, 0, 300));   
    histLep4Pt[i]->Sumw2();
  
    histMET.push_back(new TH1D(Form("histMET_%s",processLabels[i].c_str()), ";  MET [GeV]; Number of Events", 10, 0, 200));   
    histMET[i]->Sumw2();
  
    histLep3MT.push_back(new TH1D(Form("histLep3MT_%s",processLabels[i].c_str()), ";  [GeV/c^{2}]; Number of Events", 15, 0, 300));   
    histLep3MT[i]->Sumw2();
  
    histLep4MT.push_back(new TH1D(Form("histLep4MT_%s",processLabels[i].c_str()), ";  [GeV/c^{2}]; Number of Events", 15, 0, 300));   
    histLep4MT[i]->Sumw2();
  
    histLep34MT.push_back(new TH1D(Form("histLep34MT_%s",processLabels[i].c_str()), ";  [GeV/c^{2}]; Number of Events", 15, 0, 300));   
    histLep34MT[i]->Sumw2();
  
    histLep1234MT.push_back(new TH1D(Form("histLep1234MT_%s",processLabels[i].c_str()), ";  m_{T} (4L,MET) [GeV/c^{2}]; Number of Events", 25, 0, 500));   
    histLep1234MT[i]->Sumw2();
  
    histMass4L.push_back(new TH1D(Form("histMass4L_%s",processLabels[i].c_str()), ";  m_{4L} [GeV/c^{2}]; Number of Events", 20, 0, 1000));   
    histMass4L[i]->Sumw2();
  
    histLep34Mass.push_back(new TH1D(Form("histLep34Mass_%s",processLabels[i].c_str()), ";  [GeV/c^{2}]; Number of Events", 30, 0, 300));   
    histLep34Mass[i]->Sumw2();
  
    histNBJet20.push_back(new TH1D(Form("histNBJet20_%s",processLabels[i].c_str()), ";  [GeV/c^{2}]; Number of Events", 6, -0.5, 5.5));   
    histNBJet20[i]->Sumw2();
  
    histNBJet30.push_back(new TH1D(Form("histNBJet30_%s",processLabels[i].c_str()), ";  [GeV/c^{2}]; Number of Events", 6, -0.5, 5.5));   
    histNBJet30[i]->Sumw2();
  
    histNJet20.push_back(new TH1D(Form("histNJet20_%s",processLabels[i].c_str()), ";  [GeV/c^{2}]; Number of Events", 6, -0.5, 5.5));   
    histNJet20[i]->Sumw2();
  
    histNJet30.push_back(new TH1D(Form("histNJet30_%s",processLabels[i].c_str()), ";  [GeV/c^{2}]; Number of Events", 6, -0.5, 5.5));   
    histNJet30[i]->Sumw2();    

    histPhi0.push_back(new TH1D(Form("histPhi0_%s",processLabels[i].c_str()), ";  #Phi_{0}; Number of Events", 6, 0, 2*3.142));   
    histPhi0[i]->Sumw2();    

    histTheta0.push_back(new TH1D(Form("histTheta0_%s",processLabels[i].c_str()), ";  #Theta_{0}; Number of Events", 6, 0, 3.142));   
    histTheta0[i]->Sumw2();    

    histPhi.push_back(new TH1D(Form("histPhi_%s",processLabels[i].c_str()), ";  #Phi; Number of Events", 6, 0, 2*3.142));   
    histPhi[i]->Sumw2();    

    histTheta1.push_back(new TH1D(Form("histTheta1_%s",processLabels[i].c_str()), ";  #Theta_{1}; Number of Events", 6, 0, 3.142));   
    histTheta1[i]->Sumw2();    

    histTheta2.push_back(new TH1D(Form("histTheta2_%s",processLabels[i].c_str()), ";  #Theta_{2}; Number of Events", 6, 0, 3.142));   
    histTheta2[i]->Sumw2();    

  }
 
  double dataYield = 0;


  //*******************************************************************************************
  //Read file
  //*******************************************************************************************                
  for (uint i=0; i < inputfiles.size(); ++i) {

    //for duplicate event checking
    map<pair<uint,uint>, bool > processedRunEvents;

    for (uint j=0; j < inputfiles[i].size(); ++j) {
      if (inputfiles[i][j] == "") continue;
      TFile* inputFile = new TFile(inputfiles[i][j].c_str(),"READ");
      assert(inputFile);
      TTree* tree = 0;
      tree = (TTree*)inputFile->Get("WWZAnalysis");    
      assert(tree);
      

      float eventweight;
      float pileupWeight, pileupWeightUp, pileupWeightDown;
      float triggerEffWeight;
      float triggerEffSFWeight;
      uint NPU, nPV;
      int lep1Id;
      float lep1Pt, lep1Eta, lep1Phi;
      bool lep1PassLooseMVAID;
      int lep2Id;
      float lep2Pt, lep2Eta, lep2Phi;
      bool lep2PassLooseMVAID;
      int lep3Id;
      float lep3Pt, lep3Eta, lep3Phi;
      bool lep3PassLooseMVAID;
      int lep4Id;
      float lep4Pt, lep4Eta, lep4Phi;
      bool lep4PassLooseMVAID;
      float ZMass, ZPt;
      float lep3MT, lep4MT;
      float lep34MT;
      int NJet20;
      int NJet30;
      int NBJet20;
      int NBJet30;
      float minDRJetToLep3;
      float minDRJetToLep4;      
      float MET, MET_JESUp, MET_JESDown; 
      float METPhi;
      unsigned int run, lumiSec, event;
      float phi0, theta0, phi, theta1, theta2, phiH;
      bool HLTDecision[300];

      tree->SetBranchAddress("weight", &eventweight);
      tree->SetBranchAddress("pileupWeight", &pileupWeight);
      tree->SetBranchAddress("pileupWeightUp", &pileupWeightUp);
      tree->SetBranchAddress("pileupWeightDown", &pileupWeightDown);
      tree->SetBranchAddress("triggerEffWeight", &triggerEffWeight);
      tree->SetBranchAddress("triggerEffSFWeight", &triggerEffSFWeight);
      tree->SetBranchAddress("run", &run);
      tree->SetBranchAddress("lumi", &lumiSec);
      tree->SetBranchAddress("event", &event);
      tree->SetBranchAddress("NPU", &NPU);
      tree->SetBranchAddress("nPV", &nPV);
      tree->SetBranchAddress("MET", &MET);
      tree->SetBranchAddress("MET_JESUp", &MET_JESUp);
      tree->SetBranchAddress("MET_JESDown", &MET_JESDown);
      tree->SetBranchAddress("METPhi", &METPhi);
      tree->SetBranchAddress("NJet20", &NJet20);
      tree->SetBranchAddress("NJet30", &NJet30);
      tree->SetBranchAddress("NBJet20", &NBJet20);
      tree->SetBranchAddress("NBJet30", &NBJet30);
      tree->SetBranchAddress("lep1Id", &lep1Id);
      tree->SetBranchAddress("lep1Pt", &lep1Pt);
      tree->SetBranchAddress("lep1Eta", &lep1Eta);
      tree->SetBranchAddress("lep1Phi", &lep1Phi);
      tree->SetBranchAddress("lep1PassLooseMVAID", &lep1PassLooseMVAID);
      tree->SetBranchAddress("lep2Id", &lep2Id);
      tree->SetBranchAddress("lep2Pt", &lep2Pt);
      tree->SetBranchAddress("lep2Eta", &lep2Eta);
      tree->SetBranchAddress("lep2Phi", &lep2Phi);
      tree->SetBranchAddress("lep2PassLooseMVAID", &lep2PassLooseMVAID);
      tree->SetBranchAddress("lep3Id", &lep3Id);
      tree->SetBranchAddress("lep3Pt", &lep3Pt);
      tree->SetBranchAddress("lep3Eta", &lep3Eta);
      tree->SetBranchAddress("lep3Phi", &lep3Phi);
      tree->SetBranchAddress("lep3PassLooseMVAID", &lep3PassLooseMVAID);
      tree->SetBranchAddress("lep4Id", &lep4Id);
      tree->SetBranchAddress("lep4Pt", &lep4Pt);
      tree->SetBranchAddress("lep4Eta", &lep4Eta);
      tree->SetBranchAddress("lep4Phi", &lep4Phi);
      tree->SetBranchAddress("lep4PassLooseMVAID", &lep4PassLooseMVAID);
      tree->SetBranchAddress("ZMass", &ZMass);
      tree->SetBranchAddress("ZPt", &ZPt);
      tree->SetBranchAddress("lep3MT", &lep3MT);
      tree->SetBranchAddress("lep4MT", &lep4MT);
      tree->SetBranchAddress("lep34MT", &lep34MT);
      // tree->SetBranchAddress("minDRJetToLep3", &minDRJetToLep3);
      // tree->SetBranchAddress("minDRJetToLep4", &minDRJetToLep4);
      tree->SetBranchAddress("phi0", &phi0);
      tree->SetBranchAddress("theta0", &theta0);
      tree->SetBranchAddress("phi", &phi);
      tree->SetBranchAddress("theta1", &theta1);
      tree->SetBranchAddress("theta2", &theta2);
      tree->SetBranchAddress("phiH", &phiH);
      tree->SetBranchAddress("HLTDecision", &HLTDecision);


      bool isData = false;
      if ( processLabels[i] == "Data") isData = true;
    
      cout << "process: " << processLabels[i] << " | file " << inputfiles[i][j] << " | Total Entries: " << tree->GetEntries() << "\n";
      for(UInt_t ientry=0; ientry < tree->GetEntries(); ientry++) {       	
	tree->GetEntry(ientry);
      
	if (ientry % 100000 == 0) cout << "Event " << ientry << endl;      

	double puWeight = 1;      
	double weight = 1;
	if (!isData) {
	  weight = lumi * eventweight;
	}

	if (isnan(eventweight) || isinf(eventweight)) {
	  cout << "...bad event: " << eventweight << "\n";
	  continue;
	}


	//******************************
	//Trigger Selection
	//******************************
	bool passTrigger = false;

	// //Use Single Lepton Triggers
	// if ( events->HLTDecision[2] || events->HLTDecision[7] || events->HLTDecision[12] || events->HLTDecision[11] || events->HLTDecision[15])  
	//   passTrigger = true;

	// if (isData) {
	//   if ( events->HLTDecision[22] || events->HLTDecision[23] || events->HLTDecision[24] || 
	//        events->HLTDecision[25] || events->HLTDecision[26] ||
	//        events->HLTDecision[27] || events->HLTDecision[28] || events->HLTDecision[29]	  
	//        ) passTrigger = true;
	// } else {
	//   if ( events->HLTDecision[18] || events->HLTDecision[19] || events->HLTDecision[20] || 
	//        events->HLTDecision[21] ||
	//        events->HLTDecision[28] || events->HLTDecision[29]	  
	//        ) passTrigger = true;
	// }

	passTrigger = true;
	if (!passTrigger) continue;


	//******************************
	//Lorentz Vectors
	//******************************
	TLorentzVector vLep1;
	double lep1Mass = 0;
	if (abs(lep1Id) == 11) lep1Mass = 0.000511;
	else if (abs(lep1Id) == 13) lep1Mass = 0.1057;
	vLep1.SetPtEtaPhiM(lep1Pt, lep1Eta, lep1Phi,lep1Mass);

	TLorentzVector vLep2;
	double lep2Mass = 0;
	if (abs(lep2Id) == 11) lep2Mass = 0.000511;
	else if (abs(lep2Id) == 13) lep2Mass = 0.1057;
	vLep2.SetPtEtaPhiM(lep2Pt, lep2Eta, lep2Phi,lep2Mass);

	TLorentzVector vLep3;
	double lep3Mass = 0;
	if (abs(lep3Id) == 11) lep3Mass = 0.000511;
	else if (abs(lep3Id) == 13) lep3Mass = 0.1057;
	vLep3.SetPtEtaPhiM(lep3Pt, lep3Eta, lep3Phi,lep3Mass);

	TLorentzVector vLep4;
	double lep4Mass = 0;
	if (abs(lep4Id) == 11) lep4Mass = 0.000511;
	else if (abs(lep4Id) == 13) lep4Mass = 0.1057;
	vLep4.SetPtEtaPhiM(lep4Pt, lep4Eta, lep4Phi,lep4Mass);

	TLorentzVector ZCandidate = vLep1+vLep2;
	TLorentzVector v4L = vLep1+vLep2+vLep3+vLep4;
	TLorentzVector vLep34 = vLep3+vLep4;
	
	//******************************
	//Selection Cuts 
	//******************************
	//4Lepton
	if (!(lep1Pt>10 && lep2Pt>10 && lep3Pt>10 && lep4Pt>10 )) continue;

	//leading lepton pt selection
	double leadLeptonPt = 0;
	double subleadLeptonPt = 0;
	vector<double> leptonPtVector; 
	leptonPtVector.push_back(lep1Pt);
	leptonPtVector.push_back(lep2Pt);
	leptonPtVector.push_back(lep3Pt);
	leptonPtVector.push_back(lep4Pt);
	std::sort(leptonPtVector.begin(), leptonPtVector.end());   
	//cout << "Check: " << leptonPtVector[0] << " " << leptonPtVector[1] << " " << leptonPtVector[2] << " " << leptonPtVector[3] << "\n";
	leadLeptonPt = leptonPtVector[3];
	subleadLeptonPt = leptonPtVector[2];
	if (!(leadLeptonPt > 25 && subleadLeptonPt > 15)) continue;

	//ZMass Window
	if (!(ZMass > 76 && ZMass < 106)) continue;

	//Opposite Charge on Lep3 and Lep4
	if ( abs(lep3Id)/lep3Id ==  abs(lep4Id)/lep4Id ) continue;

	//2nd Z Veto
	if ( fabs(vLep34.M() - 91) < 15 ) continue;

	//MET 
	if (!(MET > 50)) continue;

	//BJet Veto
	if (!(NBJet20 == 0)) continue;
 
	//Jet Veto
	//if (!(NJet30 == 0)) continue;
 
	// // //require lep3 and lep4 to pass a tighter electron WP
	// if (! (lep3PassLooseMVAID && lep4PassLooseMVAID)) continue;

	// // //if lep4 is an electron, require lep4 pt > 15
	// if ( abs(lep4Id) == 11 && lep4Pt < 15) continue;

	//******************************
	//Categories
	//******************************
	//Different Flavor
	//if ( abs(lep3Id) == abs(lep4Id) ) continue;

	//Same Flavor
	//if ( abs(lep3Id) != abs(lep4Id) ) continue;


	//******************************
	//Apply Scale Factors
	//******************************
	if (!isData) {
	  double triggerEffScaleFactor = 1.0;
	  double leptonEffScaleFactor = 1.0;
	}



	//******************************
	//Fill histograms
	//******************************
	if (isData) {	 
	} else {
	  EventYield[i] += weight;
	  EventYieldErrSqr[i] += weight*weight;

	  histZMass[i]->Fill( ZCandidate.M(), weight );      
	  histZPt[i]->Fill( ZCandidate.Pt(), weight );
	  histLep3Pt[i]->Fill(lep3Pt, weight );
	  histLep4Pt[i]->Fill(lep4Pt, weight );
	  histMET[i]->Fill(MET, weight );
	  histLep3MT[i]->Fill(lep3MT, weight );
	  histLep4MT[i]->Fill(lep4MT, weight );
	  histLep34MT[i]->Fill(lep34MT, weight );
	  histMass4L[i]->Fill(v4L.M(),weight);
	  histLep34Mass[i]->Fill(vLep34.M(),weight);
	  histNBJet20[i]->Fill(NBJet20,weight);
	  histNBJet30[i]->Fill(NBJet30,weight);
	  histNJet20[i]->Fill(NJet20,weight);
	  histNJet30[i]->Fill(NJet30,weight);
	  histPhi0[i]->Fill(phi0,weight);
	  histTheta0[i]->Fill(theta0,weight);
	  histPhi[i]->Fill(phi,weight);
	  histTheta1[i]->Fill(theta1,weight);
	  histTheta2[i]->Fill(theta2,weight);

	}
      }
    }
  }
 
  //--------------------------------------------------------------------------------------------------------------
  // Make Plots
  //==============================================================================================================


  //--------------------------------------------------------------------------------------------------------------
  // Draw
  //==============================================================================================================
  TCanvas *cv =0;
  TLegend *legend = 0;

  //*******************************************************************************************
  //MR
  //*******************************************************************************************
  PlotDataAndStackedBkg( histZMass, processLabels, color, false, "ZMass", Label);
  PlotDataAndStackedBkg( histZPt, processLabels, color, false, "ZPt", Label);
  PlotDataAndStackedBkg( histLep3Pt, processLabels, color, false, "Lep3Pt", Label);
  PlotDataAndStackedBkg( histLep4Pt, processLabels, color, false, "Lep4Pt", Label);
  PlotDataAndStackedBkg( histMET, processLabels, color, false, "MET", Label);
  PlotDataAndStackedBkg( histLep3MT, processLabels, color, false, "Lep3MT", Label);
  PlotDataAndStackedBkg( histLep4MT, processLabels, color, false, "Lep4MT", Label);
  PlotDataAndStackedBkg( histLep34MT, processLabels, color, false, "Lep34MT", Label);
  PlotDataAndStackedBkg( histMass4L, processLabels, color, false, "Mass4L", Label);
  PlotDataAndStackedBkg( histLep34Mass, processLabels, color, false, "Lep34Mass", Label);
  PlotDataAndStackedBkg( histNBJet20, processLabels, color, false, "NBJet20", Label);
  PlotDataAndStackedBkg( histNBJet30, processLabels, color, false, "NBJet30", Label);
  PlotDataAndStackedBkg( histNJet20, processLabels, color, false, "NJet20", Label);
  PlotDataAndStackedBkg( histNJet30, processLabels, color, false, "NJet30", Label);
  PlotDataAndStackedBkg( histPhi0, processLabels, color, false, "Phi0", Label);
  PlotDataAndStackedBkg( histTheta0, processLabels, color, false, "Theta0", Label);
  PlotDataAndStackedBkg( histPhi, processLabels, color, false, "Phi", Label);
  PlotDataAndStackedBkg( histTheta1, processLabels, color, false, "Theta1", Label);
  PlotDataAndStackedBkg( histTheta2, processLabels, color, false, "Theta2", Label);

  //--------------------------------------------------------------------------------------------------------------
  // Tables
  //==============================================================================================================
  cout << "For Luminosity = " << lumi << " pb^-1\n";
  // cout << "Yields : MR > 300 && Rsq > 0.1\n";
  //cout << "TTJets: " << 

  cout << "Selected Event Yield \n";
  cout << "Data: " << dataYield << "\n";
  double totalMCEventYield = 0;
  double totalMCEventYieldErrSqr = 0;
  for (uint i=0; i < inputfiles.size(); ++i) {
    cout << "Process " << processLabels[i] << " : " << EventYield[i] << " +/- " <<  sqrt(EventYieldErrSqr[i]) << "\n";
    totalMCEventYield += EventYield[i];
    totalMCEventYieldErrSqr += EventYieldErrSqr[i];
  }
  cout << "Total MC Event Yield: " << totalMCEventYield << " +/- " << sqrt(totalMCEventYieldErrSqr) << "\n";

  // //--------------------------------------------------------------------------------------------------------------
  // // Output
  // //==============================================================================================================
  TFile *file = TFile::Open(("TTBarDileptonControlRegionPlots"+Label+".root").c_str(), "UPDATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histZMass[i], Form("histZMass_%s",processLabels[i].c_str()), "WriteDelete");
   }


 
  file->Close();
  delete file;       

}







void SelectWWZTo4L( int option = 0) {

  vector<string> datafiles;
  vector<vector<string> > bkgfiles;
  vector<string> processLabels;
  vector<int> colors;

  string datafile = "";


  //No Skims  
  if (option == 0 ) {
    datafiles.push_back("");  
  }
 

  vector<string> bkgfiles_WWZ;
  vector<string> bkgfiles_ttZ;
  vector<string> bkgfiles_ZZ;
  vector<string> bkgfiles_WZ;  
  vector<string> bkgfiles_ttW; 

  bkgfiles_WWZ.push_back("/eos/cms/store/group/phys_susy/razor/Run2Analysis/WWZAnalysis/2016/11012017/WWZAnalysis_WWZJetsTo4L2Nu_4f_TuneCUETP8M1_13TeV_aMCatNLOFxFx_pythia8_1pb_weighted.root");
  // bkgfiles_ttZ.push_back("/eos/cms/store/group/phys_susy/razor/Run2Analysis/WWZAnalysis/2016/04242017/WWZAnalysis_ttZJets_13TeV_madgraphMLM_1pb_weighted.root");
  // bkgfiles_ZZ.push_back("/eos/cms/store/group/phys_susy/razor/Run2Analysis/WWZAnalysis/2016/04242017/WWZAnalysis_ZZTo4L_13TeV_powheg_pythia8_1pb_weighted.root");
  bkgfiles_WZ.push_back("/eos/cms/store/group/phys_susy/razor/Run2Analysis/WWZAnalysis/2016/11012017/WWZAnalysis_WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root");
  // bkgfiles_ttW.push_back("/eos/cms/store/group/phys_susy/razor/Run2Analysis/WWZAnalysis/2016/04242017/WWZAnalysis_ttWJets_13TeV_madgraphMLM_1pb_weighted.root");


  bkgfiles.push_back(bkgfiles_WWZ);
  // bkgfiles.push_back(bkgfiles_ttZ);
  //  bkgfiles.push_back(bkgfiles_ZZ);
    bkgfiles.push_back(bkgfiles_WZ);
  // bkgfiles.push_back(bkgfiles_ttW);


  processLabels.push_back("WWZ");  
  // processLabels.push_back("t#bar{t}Z");
  //  processLabels.push_back("ZZ");
    processLabels.push_back("WZ");
  // processLabels.push_back("t#bar{t}W");

  colors.push_back(kAzure+10);
  // colors.push_back(kGreen+2);
  //  colors.push_back(kGray);
   colors.push_back(kBlue);
  // colors.push_back(kRed);

  double lumi = 35900;

  //*********************************************************************
  //E-Mu Control Region
  //*********************************************************************
  
  if (option == 0) RunSelectWWZTo4L(datafiles, bkgfiles,processLabels, colors, lumi,"WWZ",-1,"WWZ");



}


// ***********************************************
//Baseline
// Process Data : 0 +/- 0
// Process WWZ : 0.835752 +/- 0.00322842
// Process t#bar{t}Z : 0.544556 +/- 0.0362233
// Process ZZ : 1.29373 +/- 0.0266649
// Process WZ : 0.691432 +/- 0.202694
// Process t#bar{t}W : 0.00523152 +/- 0.00369924
// Total MC Event Yield: 3.3707 +/- 0.207683
// ***********************************************

// ***********************************************
// With LooseMVAID on Lep3 and Lep4
// Process Data : 0 +/- 0
// Process WWZ : 0.811107 +/- 0.00317999
// Process t#bar{t}Z : 0.525279 +/- 0.0355764
// Process ZZ : 1.15688 +/- 0.0252152
// Process WZ : 0.432145 +/- 0.140031
// Process t#bar{t}W : 0.00261576 +/- 0.00261576
// Total MC Event Yield: 2.92803 +/- 0.146721
// ***********************************************

// ***********************************************
// Require Lep4Pt>15 if Lep4 is an electron
// Process Data : 0 +/- 0
// Process WWZ : 0.795551 +/- 0.00314966
// Process t#bar{t}Z : 0.501184 +/- 0.0347509
// Process ZZ : 1.01674 +/- 0.0236386
// Process WZ : 0.38893 +/- 0.129643
// Process t#bar{t}W : 0.00261576 +/- 0.00261576
// Total MC Event Yield: 2.70502 +/- 0.136347
// ***********************************************


// Process WWZ : 1.42379 +/- 0.00421209
// Process WZ : 0.842682 +/- 0.165968


// Process WWZ : 1.49669 +/- 0.00431831
// Process WZ : 1.44768 +/- 0.249187
