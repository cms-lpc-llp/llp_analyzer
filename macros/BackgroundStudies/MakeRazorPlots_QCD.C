
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <vector>
#include <map>
#include <iostream>

const bool density = false;

const Int_t NComponents = 10;
int color[] = {kBlack, kAzure+4, kMagenta, kBlue+1, kCyan+1, kOrange-3, kRed+1, kGreen+2}; // for multijet box

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
  double intLumi = 2100; //in units of pb^-1
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
  vector<TH1F*> histUnrolled; 
  float MRBinLowEdges[] = {500, 600, 700, 900, 1200, 1600, 2500, 4000}; // Multijet Bins
  float RsqBinLowEdges[] = {0.25, 0.30, 0.41, 0.52, 0.64, 1.5}; // Multijet Bins

  const int nMRBins = 7;
  const int nRsqBins = 5;

  TH1F* histMRAllBkg =  new TH1F( "MRAllBkg",";M_{R} [GeV/c^{2}];Number of Events", nMRBins, MRBinLowEdges);
  TH1F* histRsqAllBkg =  new TH1F( "RsqAllBkg", ";R^{2};Number of Events", nRsqBins, RsqBinLowEdges);
  histMRAllBkg->SetStats(false);
  histRsqAllBkg->SetStats(false);  
  histRsqAllBkg->Sumw2();
  histMRAllBkg->Sumw2();

  TH1F* histMRQCD =  new TH1F( "MRQCD",";M_{R} [GeV/c^{2}];Number of Events", nMRBins, MRBinLowEdges);
  TH1F* histRsqQCD =  new TH1F( "RsqQCD", ";R^{2};Number of Events", nRsqBins, RsqBinLowEdges);
  histMRQCD->SetStats(false);
  histRsqQCD->SetStats(false);  
  histRsqQCD->Sumw2();
  histMRQCD->Sumw2();

  TH1F* histMRData =  new TH1F( "MRData",";M_{R} [GeV/c^{2}];Number of Events", nMRBins, MRBinLowEdges);
  TH1F* histRsqData =  new TH1F( "RsqData", ";R^{2};Number of Events", nRsqBins, RsqBinLowEdges);

  vector<TH1F*> histMR;
  vector<TH1F*> histRsq; 
  vector<TH2F*> histMRRsq;

  histMRQCD->SetFillColor(kAzure+4);
  histMRAllBkg->SetFillColor(kMagenta);
  histMRQCD->SetFillStyle(1001);
  histMRAllBkg->SetFillStyle(1001);

  histMRQCD->SetLineColor(kAzure+4);
  histMRAllBkg->SetLineColor(kMagenta);

  assert (inputfiles.size() == processLabels.size());
  for (int i=0; i < inputfiles.size(); ++i) {    
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
  }
  THStack *stackUnrolled = new THStack();

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
  // if (box == 0) {
    //   tree = (TTree*)inputFile->Get("MultiJet");
    // } else if (box == 1) {
    //   tree = (TTree*)inputFile->Get("LooseLeptonMultiJet");
    // } else if (box == 2) {
    //   tree = (TTree*)inputFile->Get("MuMultiJet");
    // } else if (box == 3) {
    //   tree = (TTree*)inputFile->Get("EleMultiJet");
    // }
 
    float weight = 0;
    int box = -1;
    int nBTaggedJets = 0;
    float dPhiRazor = 0;
    float MR = 0;
    float Rsq = 0;
    float mT = 0;

    tree->SetBranchAddress("weight",&weight);
    tree->SetBranchAddress("box",&box);
    tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
    tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
    tree->SetBranchAddress("MR",&MR);
    tree->SetBranchAddress("Rsq",&Rsq);
    tree->SetBranchAddress("mT",&mT);

    cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";
    for (int n=0;n<tree->GetEntries();n++) { 
     // for (int n=0;n<1000;n++) { 
   
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

      if (boxOption == 0) { // Multijet Box for Jamboree
	if( !(box == 11 || box == 12) ) continue;
      } 
      if (boxOption == 1) { // LeptonJet Box for Jamboree
	if( !(box == 3 || box == 4 || box == 6 || box == 7) ) continue;
      } 
      if (boxOption == 2) { // Multijet Box for Jamboree
	if( !(box == 14) ) continue;
      } 

      //apply baseline cuts
      if (!(MR > 400 && Rsq > 0.25)) continue;

      // if (!(MR < 500 )) continue;
      // if (!(Rsq < 0.3)) continue;

      if (!(fabs(dPhiRazor) > 2.8)) continue;

      if (!hasSignal || i>1) {
	histMRAllBkg->Fill(MR, intLumi*weight);
	histRsqAllBkg->Fill(Rsq, intLumi*weight);
	histMRRsq[i]->Fill(MR, Rsq, intLumi*weight);
      }

      if(i==1){
	if (intLumi*weight > 30) continue;
	float qcdweight = 1.56841;
	histMRQCD->Fill(MR, intLumi*weight*qcdweight);
	histRsqQCD->Fill(Rsq, intLumi*weight*qcdweight);	
  	histMRRsq[i]->Fill(MR, Rsq, intLumi*weight*qcdweight);
    }

      if (hasSignal && i==0) {
	histMRData->Fill(MR);
	histRsqData->Fill(Rsq);	
	histMRRsq[i]->Fill(MR, Rsq);
      }
    }

    inputFile->Close();
    delete inputFile;
  
  }
  
  std::cout<<"Data: "<<histMRData->Integral()<<", All Backgrounds: "<<histMRAllBkg->Integral()<<" , QCD: "<<histMRQCD->Integral()<<", Data2 "<< histMRRsq[0]->Integral() <<" QCD2 "<<histMRRsq[1]->Integral() <<std::endl;


  //*******************************************************************************************
  //Draw Plots
  //*******************************************************************************************
 // fill out the unrolled histograms 
  std::cout<<"Rsq bins: "<<nRsqBins<<" "<<nMRBins<<std::endl;
  for (uint i=0; i < histMRRsq.size(); ++i) {
    
    int binN = 0;

    for(int ii = 0; ii<nMRBins; ii++)
      for (int jj = 0; jj<nRsqBins; jj++)      
  	{      
  	  float value = (histMRRsq[i]->GetBinContent(ii+1, jj+1) > 0) ? histMRRsq[i]->GetBinContent(ii+1, jj+1) : 0. ;
	  
	  float Xrange = histMRRsq[i]->GetXaxis()->GetBinLowEdge(jj+2) - histMRRsq[i]->GetXaxis()->GetBinLowEdge(jj+1);
	  float Yrange = histMRRsq[i]->GetYaxis()->GetBinLowEdge(ii+2) - histMRRsq[i]->GetYaxis()->GetBinLowEdge(ii+1);

	  float area =1.;
	  
	  if(density) area = Xrange*Yrange; //normalize each bin by its area

	  histUnrolled[i]->SetBinContent(binN+1, value/area);
  	  binN++;
  	}

    histUnrolled[i]->SetMinimum(0.00001);

    if ( histUnrolled[i]->Integral() > 0) {
      if( !hasSignal || i > 0 )
	stackUnrolled->Add(histUnrolled[i]);
    }

    cout << "Process : " << processLabels[i] << "\n";	  
  }

  TCanvas *cv = 0;
  TLegend *legend = 0;
  TLatex *tex = 0;
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

  /// Unrolled plots in bins of R&MR
  TLatex t1(0.1,0.92, "CMS Preliminary");
  TLatex t2(0.6,0.92, "#sqrt{s}=13 TeV, L = 2.1 fb^{-1}");
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
  stackUnrolled->SetMinimum(0.01);
  stackUnrolled->SetMaximum(1000);
  cv->SetLogy();
  stackUnrolled->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackUnrolled->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackUnrolled->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackUnrolled->GetHists()->At(0)))->GetYaxis()->GetTitle());
  stackUnrolled->Draw();
  if(hasSignal) histUnrolled[0]->Draw("same PE");
  legend->Draw();
  t1.Draw();
  t2.Draw();
  t3.Draw();
  cv->SaveAs(Form("Unrolled_QCD%s.root",Label.c_str()));
  ////

  //*******************************************************************************************
  //MR
  //*******************************************************************************************
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.50,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.2, 0.92, Form("CMS Simulation #sqrt{s} = 13 TeV, #int L = %d fb^{-1}, %s",int(intLumi/1000), latexlabel.c_str()));

  THStack *stackMR = new THStack("stackMR", "");
  THStack *stackRsq = new THStack();

  //*******************************************************************************************
  //MR Before and After DPhi Cut
  //*******************************************************************************************
  //////////////////
  stackMR->Add(histMRAllBkg);
  stackMR->Add(histMRQCD);
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.50,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  stackMR->Draw();
  stackMR->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackMR->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackMR->GetHists()->At(0)))->GetYaxis()->GetTitle());
  stackMR->Draw("");
  histMRData->Draw("same PE");
  legend->Draw();
  cv->SetLogy();
  cv->SaveAs(Form("MRStack_QCD_%s.pdf",Label.c_str()));

  ///////////////////////
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.50,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  histMRAllBkg->SetLineColor(kRed);
  histMRAllBkg->GetYaxis()->SetTitle("Number of Events");
  histMRAllBkg->GetYaxis()->SetTitleOffset(1.2);
  histMRData->SetMarkerStyle(8);

  legend->AddEntry(histMRAllBkg, "All Backgrounds", "L");
  legend->AddEntry(histMRData, "Data", "L");

  histMRAllBkg->Add(histMRQCD, 1.0);

  histMRAllBkg->Draw("hist");
  histMRData->Draw("PE same");

  legend->Draw();
  cv->SetLogy();
  cv->SaveAs(Form("MR_QCD_%s.pdf",Label.c_str()));


  //////

  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.50,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  histRsqAllBkg->SetLineColor(kRed);
  histRsqAllBkg->GetYaxis()->SetTitle("Number of Events");
  histRsqAllBkg->GetYaxis()->SetTitleOffset(1.2);
  histRsqData->SetMarkerStyle(8);

  legend->AddEntry(histRsqAllBkg, "All Backgrounds", "L");
  legend->AddEntry(histRsqData, "Data", "L");

  histRsqAllBkg->Add(histRsqQCD, 1.0);

  histRsqAllBkg->Draw("hist");
  histRsqData->Draw("PE same");

  legend->Draw();
  cv->SetLogy();
  cv->SaveAs(Form("Rsq_QCD_%s.pdf",Label.c_str()));

  //////////////////
  histRsqQCD->SetFillColor(kAzure+4);
  histRsqAllBkg->SetFillColor(kMagenta);

  stackRsq->Add(histRsqAllBkg);
  stackRsq->Add(histRsqQCD);
  cv = new TCanvas("cv","cv", 800,600);
  legend = new TLegend(0.50,0.54,0.90,0.84);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  stackRsq->Draw();
  stackRsq->GetHistogram()->GetXaxis()->SetTitle(((TH1F*)(stackRsq->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackRsq->GetHistogram()->GetYaxis()->SetTitle(((TH1F*)(stackRsq->GetHists()->At(0)))->GetYaxis()->GetTitle());
  stackRsq->Draw();
  histRsqData->Draw("same PE");
  legend->Draw();
  cv->SetLogy();
  cv->SaveAs(Form("RsqStack_QCD_%s.pdf",Label.c_str()));
 
   //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  TFile *file = TFile::Open(("RazorPlots"+Label+".root").c_str(), "RECREATE");
  file->cd();

  for(int i=0; i<int(inputfiles.size()); i++) {
    file->WriteTObject(histMRRsq[i], Form("histMRRsq_%s",processLabels[i].c_str()), "WriteDelete");
    histUnrolled[i]->Write();  
  }
  
  stackUnrolled->Write();
 }


 void MakeRazorPlots_QCD() {

   string signalfile = "eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorInclusive_HTMHT_Run2015D_2093pb_GoodLumiGolden_RazorSkim_FilteredNew.root";
   string signalLabel = "Data";

   vector<string> bkgfiles;
   vector<string> bkgLabels;      

   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_QCD_HTBinned_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_Other_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_ST_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root");
   bkgfiles.push_back("eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/RazorSkim/RazorInclusive_TTJets_Madgraph_Leptonic_1pb_weighted_RazorSkim.root");
    
   bkgLabels.push_back("QCD");
   bkgLabels.push_back("Other");
   bkgLabels.push_back("DYJetsToLL");
   bkgLabels.push_back("ZJetsToNuNu");
   bkgLabels.push_back("SingleTop");
   bkgLabels.push_back("WJetsToLNu");
   bkgLabels.push_back("TTJets");

   RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,2,4,"test_CombinedBTag", ""); 
 }
 
