#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <vector>
#include <map>
#include <iostream>

const bool density = true;

// int color[] = {kAzure+4, kMagenta, kBlue+1, kCyan+1, kOrange-3, kRed+1, kGreen+2}; // for multijet box
// int color[] = {kBlack, kAzure+4, kMagenta, kBlue+1, kCyan+1, kOrange-3, kRed+1, kGreen+2}; // for multijet box
int color[] = {kAzure+4, kMagenta, kBlue+1, kCyan+1, kOrange-3, kRed+1, kGreen+2, kGreen+4, kGreen-4}; // for lepton box
// int color[] = {kAzure+4, kMagenta, kBlue+1, kCyan+1, kOrange-3, kRed+1, kBlack, kGreen+2}; // for multijet 0L vs 1L
// int color[] = {10, 21, 32, 43, 54, 65, 76, 87, 28, 39, 40, kAzure+4, kMagenta, kBlue+1, kCyan+1, kOrange-3, kRed+1, kGreen+2, kGreen+4, kGreen-4}; // for lepton box

//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
TH1F* NormalizeHist(TH1F *originalHist) {
  TH1F* hist = (TH1F*)originalHist->Clone((string(originalHist->GetName())+"_normalized").c_str());
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return hist;
}


//------------------------------------------------------------------------------
// PlotHiggsRes_LP
//------------------------------------------------------------------------------
void RunMakeRazorPlots ( string signalfile, string signalLabel,  vector<string> bkgfiles,vector<string> bkgLabels, int boxOption = 0, int option = -1, string label = "", string latexlabel = "") {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  double intLumi = 2000; //in units of pb^-1
  string Label = "";
  if (label != "") Label = "_" + label;

  vector<string> inputfiles;
  vector<string> processLabels;

  bool hasSignal = false;
  if (signalfile != "") {
    hasSignal = true;
    inputfiles.push_back(signalfile);
    processLabels.push_back(signalLabel);
  }
  assert(bkgfiles.size() == bkgLabels.size());
  for (int i=0; i < bkgfiles.size(); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processLabels.push_back(bkgLabels[i]);
  }

  //*******************************************************************************************
  //Define Histograms
  //*******************************************************************************************
  TH1F* histMRAllBkg =  new TH1F( "MRAllBkg",";M_{R} [GeV/c^{2}];Number of Events", 100, 400, 2400);
  TH1F* histRsqAllBkg =  new TH1F( "RsqAllBkg", ";R^{2};Number of Events", 24, 0.25, 1.45);
  histMRAllBkg->SetStats(false);
  histRsqAllBkg->SetStats(false);
  
  vector<TH1F*> histMR;
  vector<TH1F*> histRsq; 

  vector<TH2F*> histMRRsq;
  vector<TH1F*> histUnrolled; 
  vector<TH1F*> histUnrolled2bins; 
  vector<TH1F*> histUnrolledPercentage; 
  vector<TH1F*> histUnrolledPercentage2bins; 

  //  float MRBinLowEdges[] = {500, 600, 700, 900, 1200, 1600, 2500, 4000};
  //  float RsqBinLowEdges[] = {0.25, 0.30, 0.41, 0.52, 0.64, 1.5};

  // float MRBinLowEdges[] = {500, 600, 700, 900, 1200, 1600, 2500, 4000}; // Multijet Bins
  // float RsqBinLowEdges[] = {0.25, 0.30, 0.41, 0.52, 0.64, 1.5}; // Multijet Bins

  float MRBinLowEdges[] = {400, 500, 600, 700, 900, 1200, 1600, 2500, 4000}; // Lepton boxes
  float RsqBinLowEdges[] = {0.15, 0.20, 0.25, 0.30, 0.41, 0.52, 0.64, 1.5};  // Lepton boxes


  const int nMRBins = sizeof(MRBinLowEdges)/sizeof(float)-1;
  const int nRsqBins = sizeof(RsqBinLowEdges)/sizeof(float)-1;

  std::cout<<"AAAAAAA "<<nMRBins<<" "<<nRsqBins<<std::endl;

  assert (inputfiles.size() == processLabels.size());

  for (int i=0; i < inputfiles.size(); ++i) {    
    // histMR.push_back( new TH1F( Form("MR_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}];Number of Events", 100, 400, 2400));
    histMR.push_back( new TH1F( Form("MR_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}];Number of Events", nMRBins, MRBinLowEdges));
    if (!hasSignal || i != 0) histMR[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histMR[i]->SetLineWidth(3);
    histMR[i]->SetLineColor(color[i]);    
    histMR[i]->SetStats(false);    
    histMR[i]->Sumw2();

    // histRsq.push_back( new TH1F( Form("Rsq_%s",processLabels[i].c_str()), ";R^{2} ;Number of Events", 24, 0.25, 1.45));
    histRsq.push_back( new TH1F( Form("Rsq_%s",processLabels[i].c_str()), ";R^{2} ;Number of Events", nRsqBins, RsqBinLowEdges));
    if (!hasSignal || i != 0) histRsq[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histRsq[i]->SetLineWidth(3);
    histRsq[i]->SetLineColor(color[i]);
    histRsq[i]->SetStats(false);

    histMRRsq.push_back( new TH2F( Form("MRRsq_%s",processLabels[i].c_str()), ";M_{R} [GeV/c^{2}]; R^{2}", nMRBins, MRBinLowEdges, nRsqBins, RsqBinLowEdges));
    if (!hasSignal || i != 0) histMRRsq[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histMRRsq[i]->SetLineWidth(3);
    histMRRsq[i]->SetLineColor(color[i]);
    histMRRsq[i]->SetStats(false);
    histMRRsq[i]->Sumw2();

    histUnrolled.push_back( new TH1F( Form("Unrolled_%s",processLabels[i].c_str()), ";Bin Number ;Number of Events", nMRBins*nRsqBins, 0, nMRBins*nRsqBins));
    if (!hasSignal || i != 0) histUnrolled[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histUnrolled[i]->SetLineWidth(3);
    histUnrolled[i]->SetLineColor(color[i]);
    histUnrolled[i]->SetStats(false);     

    histUnrolled2bins.push_back( new TH1F( Form("Unrolled2bins_%s",processLabels[i].c_str()), ";Bin Number ;Event Density", 3, 0, 3));
    if (!hasSignal || i != 0) histUnrolled2bins[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histUnrolled2bins[i]->SetLineWidth(3);
    histUnrolled2bins[i]->SetLineColor(color[i]);
    histUnrolled2bins[i]->SetStats(false);     

    histUnrolledPercentage2bins.push_back( new TH1F( Form("UnrolledPercentage2bins_%s",processLabels[i].c_str()), ";;Event Density", 3, 0, 3));
    if (!hasSignal || i != 0) histUnrolledPercentage2bins[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histUnrolledPercentage2bins[i]->SetLineWidth(3);
    histUnrolledPercentage2bins[i]->SetLineColor(color[i]);
    histUnrolled2bins[i]->SetStats(false);     

    histUnrolledPercentage.push_back( new TH1F( Form("UnrolledPercentage_%s",processLabels[i].c_str()), ";Bin Number ; Fraction of total", nMRBins*nRsqBins, 0, nMRBins*nRsqBins));
    if (!hasSignal || i != 0) histUnrolledPercentage[i]->SetFillColor(color[i]);
    if (hasSignal && i==0) histUnrolledPercentage[i]->SetLineWidth(3);
    histUnrolledPercentage[i]->SetLineColor(color[i]);
    histUnrolledPercentage[i]->SetStats(false);     
}

  //*******************************************************************************************
  //Define Counts
  //*******************************************************************************************


  //*******************************************************************************************
  //Read files
  //*******************************************************************************************
  for (uint i=0; i < inputfiles.size(); ++i) {

    TFile* inputFile = new TFile(inputfiles[i].c_str(),"READ");
    assert(inputFile);
    TTree* tree = 0;
    tree = (TTree*)inputFile->Get("RazorInclusive");
 
    float weight = 0;
    int box = -1;
    int nBTaggedJets = 0;
    float dPhiRazor = 0;
    float MR = 0;
    float Rsq = 0;
    float mT = 0;
    int nGenMuons = 0;
    int nGenElectrons = 0;
    int nGenTaus = 0;

    // bool Flag_HBHENoiseFilter = false;
    // bool Flag_goodVertices    = false;
    // bool Flag_eeBadScFilter   = false;
    // bool Flag_EcalDeadCellTriggerPrimitiveFilter = false;


    tree->SetBranchAddress("weight",&weight);
    tree->SetBranchAddress("box",&box);
    tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
    tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
    tree->SetBranchAddress("MR",&MR);
    tree->SetBranchAddress("Rsq",&Rsq);
    tree->SetBranchAddress("mT",&mT);
    tree->SetBranchAddress("nGenMuons",&nGenMuons);
    tree->SetBranchAddress("nGenElectrons",&nGenElectrons);
    tree->SetBranchAddress("nGenTaus",&nGenTaus);
    // tree->SetBranchAddress("Flag_HBHENoiseFilter",&Flag_HBHENoiseFilter);
    // tree->SetBranchAddress("Flag_goodVertices",&Flag_goodVertices);
    // tree->SetBranchAddress("Flag_eeBadScFilter",&Flag_eeBadScFilter);
    // tree->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter",&Flag_EcalDeadCellTriggerPrimitiveFilter);
 

    cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";
    for (int n=0;n<tree->GetEntries();n++) { 
    // for (int n=0;n<10000;n++) { 
    
      tree->GetEntry(n);
      if (n % 1000000 == 0) cout << "Processing Event " << n << "\n";       

      // if (intLumi*weight > 100) continue;

      //Box Options
      if (option == 0 ) {
	if (nBTaggedJets != 0) continue; 
      }
      if (option == 1 ) {
	if (nBTaggedJets != 1) continue;
      }
      if (option == 2 ) {
	if (nBTaggedJets != 2) continue;
      }
      if (option == 3 ) {
	if (nBTaggedJets < 3) continue;
      }
      if (option == 4 ) {
	if (nBTaggedJets < 0) continue; // all b-tag categories combined
      }

     //Box Options
      if (boxOption == 0) { // Multijet Box for Jamboree
	if( !(box == 11 || box == 12) ) continue;
      }
      if (boxOption == 1) { // MuonMultijet Box for Jamboree
	if( !(box == 3 || box == 4) ) continue;
      } 
      if (boxOption == 2) { // EleMultijet Box for Jamboree
	if( !(box == 6 || box == 7) ) continue;
      }

      // LeptonMultijet Box for Jamboree
      if(boxOption == 1 || boxOption == 2)
	if(mT<120) continue;

      // Multijet Box for Jamboree
      if (boxOption == 0)
	if(fabs(dPhiRazor) > 2.8) continue;

      //apply baseline cuts
      if(boxOption == 1 || boxOption == 2) 
	if (!(MR > 400 && Rsq > 0.15)) continue;
      
      if(boxOption == 0) 
	if (!(MR > 500 && Rsq > 0.25)) continue;

      // if(!Flag_HBHENoiseFilter) continue;
      // if(!Flag_goodVertices) continue;
      // if(!Flag_EcalDeadCellTriggerPrimitiveFilter) continue;
      // if(!Flag_eeBadScFilter) continue;
      
      // fill the histos
      if (!hasSignal || i>0) {
	histMRAllBkg->Fill(MR, intLumi*weight);
	histRsqAllBkg->Fill(Rsq, intLumi*weight);
      }
            
      if(strstr(processLabels[i].c_str(), "TTJets")==NULL && strstr(processLabels[i].c_str(), "T1bbbb")==NULL)
      	histMRRsq[i]->Fill(MR, Rsq, intLumi*weight);
      
      // if(strstr(processLabels[i].c_str(), "T1bbbb")!=NULL)
      // 	histMRRsq[i]->Fill(MR, Rsq, intLumi*weight*2.69506e-07);
      
      histMR[i]->Fill(MR, intLumi*weight);
      histRsq[i]->Fill(Rsq, intLumi*weight);
      
      // separate by number of gen leptons for lepton boxes
      if(boxOption==1 || boxOption==2)
	{
	  if(strstr(processLabels[i].c_str(), "TTJets")!=NULL && strstr(bkgLabels[i].c_str(), "2L")!=NULL) {
	    if(nGenMuons+nGenElectrons>=2)
	      histMRRsq[i]->Fill(MR, Rsq, intLumi*weight);
	  }
	  else if(strstr(processLabels[i].c_str(), "TTJets")!=NULL && strstr(bkgLabels[i].c_str(), "Tau")!=NULL) {
	    if((nGenMuons+nGenElectrons+nGenTaus>=2) && !(nGenMuons+nGenElectrons>=2))
	      histMRRsq[i]->Fill(MR, Rsq, intLumi*weight);
	  }
	  else if(strstr(processLabels[i].c_str(), "TTJets")!=NULL && strstr(bkgLabels[i].c_str(), "1L")!=NULL) {
	    if(!(nGenMuons+nGenElectrons>=2) && !(nGenMuons+nGenElectrons+nGenTaus>=2))
	      histMRRsq[i]->Fill(MR, Rsq, intLumi*weight);
	  }
	}
      
	// Multijet box top
      	if(boxOption==0)
	  if(strstr(processLabels[i].c_str(), "TTJets")!=NULL) {
	    histMRRsq[i]->Fill(MR, Rsq, intLumi*weight);
	  }
      }

    inputFile->Close();
    delete inputFile;
  }
  
  //*******************************************************************************************
  //Draw Plots
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.7,0.53,0.90,0.88);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  for (Int_t i = histMRRsq.size()-1 ; i >= 0; --i) {
    if (hasSignal && i==0) {
      legend->AddEntry(histMRRsq[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histMRRsq[i],processLabels[i].c_str(), "F");
    }
  }

  THStack *stackUnrolled = new THStack();
  THStack *stackUnrolled2bins = new THStack();
  THStack *stackUnrolledPercentage = new THStack();
  THStack *stackUnrolledPercentage2bins = new THStack();

  float bintotal[nMRBins*nRsqBins] = {0.};

  // fill out the unrolled histograms 
  for (uint i=0; i < histMRRsq.size(); ++i) {
    
    int binN = 0;
    float total_SB = 0.;
    float total_SR = 0.;

    for(int ii = 0; ii<nMRBins; ii++)
      for (int jj = 0; jj<nRsqBins; jj++)      
  	{      
  	  float value = (histMRRsq[i]->GetBinContent(ii+1, jj+1) > 0) ? histMRRsq[i]->GetBinContent(ii+1, jj+1) : 0. ;
	  
	  float Xrange = histMRRsq[i]->GetXaxis()->GetBinLowEdge(ii+2) - histMRRsq[i]->GetXaxis()->GetBinLowEdge(ii+1);
	  float Yrange = histMRRsq[i]->GetYaxis()->GetBinLowEdge(jj+2) - histMRRsq[i]->GetYaxis()->GetBinLowEdge(jj+1);

	  float area =1.;
	  
	  if(density) area = Xrange*Yrange; //normalize each bin by its area

	  histUnrolled[i]->SetBinContent(binN+1, value/area);

	  if(!hasSignal || i>0)
	    bintotal[binN+1] += value/area;	 
	  
	  if(ii<1 || jj<1)
	    total_SB += value/area;
	  else
	    total_SR += value/area;
	  	  
  	  binN++;
  	}

    histUnrolled2bins[i]->SetBinContent(1, total_SB);
    histUnrolled2bins[i]->SetBinContent(2, total_SR);

    histUnrolled[i]->SetMinimum(0.00001);

    if ( histUnrolled[i]->Integral() > 0) {
      if( !hasSignal || i > 0 )
	stackUnrolled->Add(histUnrolled[i]);
    }

    if ( histUnrolled[i]->Integral() > 0) {
      if( !hasSignal || i > 0 )
	stackUnrolled2bins->Add(histUnrolled2bins[i]);
    }

    cout << "Process : " << processLabels[i] << "\n";	  
  }


  // Unroll into two bins for fractions
  float AllBkg_SB = 0;
  float AllBkg_SR = 0;

  for (uint i=0; i < histMRRsq.size(); ++i) {
    if( !hasSignal || i > 0 ){
      AllBkg_SB += histUnrolled2bins[i]->GetBinContent(1);
      AllBkg_SR += histUnrolled2bins[i]->GetBinContent(2);
    }
  }

  for (uint i=0; i < histMRRsq.size(); ++i) {
    if( !hasSignal || i > 0 ){
      histUnrolledPercentage2bins[i]->SetBinContent(1, histUnrolled2bins[i]->GetBinContent(1)/AllBkg_SB);
      histUnrolledPercentage2bins[i]->SetBinContent(2, histUnrolled2bins[i]->GetBinContent(2)/AllBkg_SR);
    }
    
    if ( histUnrolled2bins[i]->Integral() > 0) {
      if( !hasSignal || i > 0 )
	stackUnrolledPercentage2bins->Add(histUnrolledPercentage2bins[i]);
    }
  }
  
  ///
  // fill out the unrolled percentage histograms 
  for (uint i=0; i < histMRRsq.size(); ++i) {
    if( hasSignal && i == 0 ) continue;
      
    int binN = 0;

    for(int ii = 0; ii<nMRBins; ii++)
      for (int jj = 0; jj<nRsqBins; jj++)      
  	{      
  	  float value = (histMRRsq[i]->GetBinContent(ii+1, jj+1) > 0) ? histMRRsq[i]->GetBinContent(ii+1, jj+1) : 0. ;
	  
	  float Xrange = histMRRsq[i]->GetXaxis()->GetBinLowEdge(ii+2) - histMRRsq[i]->GetXaxis()->GetBinLowEdge(ii+1);
	  float Yrange = histMRRsq[i]->GetYaxis()->GetBinLowEdge(jj+2) - histMRRsq[i]->GetYaxis()->GetBinLowEdge(jj+1);

	  float area =1.;
	  
	  if(density) area = Xrange*Yrange; //normalize each bin by its area

	  if(bintotal[binN+1]>0)
	    histUnrolledPercentage[i]->SetBinContent(binN+1, (value/area)/bintotal[binN+1]);

	  binN++;
	}
		
    if ( histUnrolled[i]->Integral() > 0) {
	stackUnrolledPercentage->Add(histUnrolledPercentage[i]);
    }

    cout << "Unrolling Percentage for Process : " << processLabels[i] << "\n";	  
  }

  /// Unrolled plots in bins of R&MR
  TLatex t1(0.1,0.92, "CMS Preliminary");
  TLatex t2(0.6,0.92, "#sqrt{s}=13 TeV, L = 2 fb^{-1}");
  TLatex t3(0.4,0.92, Form("%s",latexlabel.c_str()) );
  t1.SetNDC();
  t2.SetNDC();
  t3.SetNDC();
  t1.SetTextSize(0.05);
  t2.SetTextSize(0.05);
  t3.SetTextSize(0.02);
  t1.SetTextFont(42);
  t2.SetTextFont(42);
  t3.SetTextFont(42);

  stackUnrolled->Draw();
  stackUnrolled->SetMinimum(0.0001);
  // stackUnrolled->SetMaximum(1000);
  cv->SetLogy();
  stackUnrolled->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackUnrolled->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackUnrolled->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackUnrolled->GetHists()->At(0)))->GetYaxis()->GetTitle());
  stackUnrolled->Draw();
  if(hasSignal) histUnrolled[0]->Draw("same hist");
  legend->Draw();
  t1.Draw();
  t2.Draw();
  t3.Draw();
  cv->SaveAs(Form("Unrolled%s.pdf",Label.c_str()));

  // Unrolled plots in percentages
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.85,0.20,0.95,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);

  for (Int_t i = histMRRsq.size()-1 ; i >= 0; --i) {
    if (hasSignal && i==0) {
      legend->AddEntry(histMRRsq[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histMRRsq[i],processLabels[i].c_str(), "F");
    }
  }
  stackUnrolledPercentage->Draw();
  stackUnrolledPercentage->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackUnrolledPercentage->GetHists()->At(0)))->GetXaxis()->GetTitle());
  // stackUnrolledPercentage->GetHistogram()->GetXaxis()->SetRangeUser(0, 35);
  stackUnrolledPercentage->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackUnrolledPercentage->GetHists()->At(0)))->GetYaxis()->GetTitle());
  if(hasSignal) histUnrolledPercentage[0]->Draw("same hist"); 
  legend->Draw();
  t1.Draw();
  t2.Draw();
  t3.Draw();
  cv->SaveAs(Form("UnrolledPercentage%s.pdf",Label.c_str()));

  // Unrolled plots in sideband vs signal box
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.85,0.20,0.95,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);

  for (Int_t i = histMRRsq.size()-1 ; i >= 0; --i) {
    if (hasSignal && i==0) {
      legend->AddEntry(histMRRsq[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histMRRsq[i],processLabels[i].c_str(), "F");
    }
  }
  stackUnrolled2bins->Draw();
  stackUnrolled2bins->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackUnrolled2bins->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackUnrolled2bins->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackUnrolled2bins->GetHists()->At(0)))->GetYaxis()->GetTitle());
  if(hasSignal) histUnrolled2bins[0]->Draw("same hist");
  legend->Draw();
  t1.Draw();
  t2.Draw();
  t3.Draw();
  cv->SaveAs(Form("Unrolled2bins%s.pdf",Label.c_str()));

  // Unrolled plots in sideband vs signal box in fractions
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.7,0.23,0.90,0.88);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);

  for (Int_t i = histMRRsq.size()-1 ; i >= 0; --i) {
    if (hasSignal && i==0) {
      legend->AddEntry(histMRRsq[i],processLabels[i].c_str(), "L");
    } else {
      legend->AddEntry(histMRRsq[i],processLabels[i].c_str(), "F");
    }
  }

  stackUnrolledPercentage2bins->Draw();
  stackUnrolledPercentage2bins->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackUnrolledPercentage2bins->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackUnrolledPercentage2bins->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackUnrolledPercentage2bins->GetHists()->At(0)))->GetYaxis()->GetTitle());
  stackUnrolledPercentage2bins->GetHistogram()->GetXaxis()->SetBinLabel(1, "Sideband");
  stackUnrolledPercentage2bins->GetHistogram()->GetXaxis()->SetBinLabel(2, "Signal Sensitive Region");
  if(hasSignal) histUnrolledPercentage2bins[0]->Draw("same hist");
  legend->Draw();
  t1.Draw();
  t2.Draw();
  t3.Draw();
  cv->SaveAs(Form("UnrolledPercentage2bins%s.pdf",Label.c_str()));


   //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("RazorPlots"+Label+".root").c_str(), "RECREATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMR[i], Form("histMR_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histRsq[i], Form("histRsq_%s",processLabels[i].c_str()), "WriteDelete");
    file->WriteTObject(histMRRsq[i], Form("histMRRsq_%s",processLabels[i].c_str()), "WriteDelete");
    histUnrolled[i]->Write();  
    histUnrolled2bins[i]->Write();  
    histUnrolledPercentage[i]->Write();  
    histUnrolledPercentage2bins[i]->Write();  
  }
  
  stackUnrolled->Write();
  stackUnrolled2bins->Write();
  stackUnrolledPercentage->Write();
  stackUnrolledPercentage2bins->Write();
 }


 void MakeRazorPlots_Unrolled() {
   
   vector<string> bkgfiles;
   vector<string> bkgLabels;      

   //   string signalfile = "eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_Other_1pb_weighted_RazorSkim.root";
   // string signalfile = "eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForPreapp20151122/jobs/combined/SMS-T1bbbb_1500_100.root";
   // string signalLabel = "T1bbbb m_{G}=1500 m_{LSP}=100";

   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_Other_1pb_weighted_RazorSkim.root");
   // bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorInclusive_WW_TuneCUETP8M1_13TeV-pythia8_1pb_weighted.root");
   // bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorInclusive_WZ_TuneCUETP8M1_13TeV-pythia8_1pb_weighted.root");
   // bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorInclusive_ZZ_TuneCUETP8M1_13TeV-pythia8_1pb_weighted.root");
   // bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorInclusive_WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_1pb_weighted.root");
   // bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorInclusive_WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_1pb_weighted.root");
   // bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorInclusive_ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_1pb_weighted.root");
   // bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorInclusive_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_1pb_weighted.root");
   // bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorInclusive_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_1pb_weighted.root");
   // bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorInclusive_TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_1pb_weighted.root");
   // bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorInclusive_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_1pb_weighted.root");
   // bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorInclusive_TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8_1pb_weighted.root");



   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_QCD_HTBinned_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_ST_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_TTJets_Madgraph_Leptonic_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_TTJets_Madgraph_Leptonic_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_TTJets_Madgraph_Leptonic_1pb_weighted_RazorSkim.root");

   bkgLabels.push_back("Other");
   // bkgLabels.push_back("WW");
   // bkgLabels.push_back("WZ");
   // bkgLabels.push_back("ZZ");
   // bkgLabels.push_back("WWZ");
   // bkgLabels.push_back("WZZ");
   // bkgLabels.push_back("ZZZ");
   // bkgLabels.push_back("TTZToLLNuNu");
   // bkgLabels.push_back("TTZToQQ");
   // bkgLabels.push_back("TTWToQQ");
   // bkgLabels.push_back("TTWToLnu");
   // bkgLabels.push_back("TTTT");

   bkgLabels.push_back("QCD");
   bkgLabels.push_back("DYJetsToLL");
   bkgLabels.push_back("ZJetsToNuNu");
   bkgLabels.push_back("SingleTop");
   bkgLabels.push_back("WJetsToLNu");
   // bkgLabels.push_back("TTJets");
   bkgLabels.push_back("TTJets 1L");
   bkgLabels.push_back("TTJets 2L");
   bkgLabels.push_back("TTJets L+Tau");

   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,1,0,"MuonMultijet_0BTag_Others", "MuonMultijet Box 0 b-tag");
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,1,1,"MuonMultijet_1BTag_Others", "MuonMultijet Box 1 b-tag");
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,1,2,"MuonMultijet_2BTag_Others", "MuonMultijet Box 2 b-tag");
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,1,3,"MuonMultijet_3BTag_Others", "MuonMultijet Box 3 b-tag");
 
  // RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,0,0,"T1bbbb_1500_100_MultiJet_0BTag", "MultiJet Box 0 b-tag");
  // RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,0,1,"T1bbbb_1500_100_MultiJet_1BTag", "MultiJet Box 1 b-tag");
  // RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,0,2,"T1bbbb_1500_100_MultiJet_2BTag", "MultiJet Box 2 b-tag");
  // RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,0,3,"T1bbbb_1500_100_MultiJet_3BTag", "MultiJet Box 3 b-tag");
 
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,0,0,"MultiJet_0BTag", "MultiJet Box 0 b-tag");
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,0,1,"MultiJet_1BTag", "MultiJet Box 1 b-tag");
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,0,2,"MultiJet_2BTag", "MultiJet Box 2 b-tag");
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,0,3,"MultiJet_3BTag", "MultiJet Box 3 b-tag");
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,0,4,"MultiJet_CombinedBTag", "MultiJet Box All b-tag");

   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,0,0,"MultiJet_0BTag_byTopLeps", "MultiJet Box 0 b-tag");
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,0,1,"MultiJet_1BTag_byTopLeps", "MultiJet Box 1 b-tag");
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,0,2,"MultiJet_2BTag_byTopLeps", "MultiJet Box 2 b-tag");
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,0,3,"MultiJet_3BTag_byTopLeps", "MultiJet Box 3 b-tag");

   RunMakeRazorPlots("","",bkgfiles,bkgLabels,2,0,"EleMultijet_0BTag", "EleMultijet Box 0 b-tag");
   RunMakeRazorPlots("","",bkgfiles,bkgLabels,2,1,"EleMultijet_1BTag", "EleMultijet Box 1 b-tag");
   RunMakeRazorPlots("","",bkgfiles,bkgLabels,2,2,"EleMultijet_2BTag", "EleMultijet Box 2 b-tag");
   RunMakeRazorPlots("","",bkgfiles,bkgLabels,2,3,"EleMultijet_3BTag", "EleMultijet Box 3 b-tag");
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,2,4,"EleMultijet_CombinedBTag", "EleMultijet Box All b-tag");

   RunMakeRazorPlots("","",bkgfiles,bkgLabels,1,0,"MuonMultijet_0BTag", "MuonMultijet Box 0 b-tag");
   RunMakeRazorPlots("","",bkgfiles,bkgLabels,1,1,"MuonMultijet_1BTag", "MuonMultijet Box 1 b-tag");
   RunMakeRazorPlots("","",bkgfiles,bkgLabels,1,2,"MuonMultijet_2BTag", "MuonMultijet Box 2 b-tag");
   RunMakeRazorPlots("","",bkgfiles,bkgLabels,1,3,"MuonMultijet_3BTag", "MuonMultijet Box 3 b-tag");
   // RunMakeRazorPlots("","",bkgfiles,bkgLabels,1,4,"MuonMultijet_CombinedBTag", "MuonMultijet Box All b-tag");




 }
 
