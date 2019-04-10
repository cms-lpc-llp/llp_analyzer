#define ZeeTiming_cxx
//#include "ZeeTiming.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TF1.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

void ZeeTiming()
{
	Float_t         weight;
	Float_t         pileupWeight;
	Float_t         pileupWeightUp;
	Float_t         pileupWeightDown;
	UInt_t          run;
	UInt_t          lumi;
	UInt_t          event;
	UInt_t          NPU;
	UInt_t          nPV;
	Float_t         mass;
	Float_t         t1;
	Float_t         t2;
	Float_t         t1_seed;
	Float_t         t2_seed;
	Float_t         t1calib_seed;
	Float_t         t2calib_seed;
	Float_t         t1raw_seed;
	Float_t         t2raw_seed;
	Float_t         ele1Pt;
	Bool_t          ele1IsEB;
	Int_t           ele1SeedIEta;
	Int_t           ele1SeedIPhi;
	Int_t           ele1SeedIX;
	Int_t           ele1SeedIY;
	Float_t         ele2Pt;
	Bool_t          ele2IsEB;
	Int_t           ele2SeedIEta;
	Int_t           ele2SeedIPhi;
	Int_t           ele2SeedIX;
	Int_t           ele2SeedIY;

	float deltaT = 0;
	float sigma = 0;

	//how many runs there are in the sample, lowR is the lowest run number, highR is the highest
	int lowR = 278800;
	int highR = 280400;
	int R = 0;
	R = highR - lowR;

	TFile *f = new TFile("ZeeTiming_SingleElectron_2016H_03Feb2017v3_25ns.Job0Of12.root");
	TTree *fChain = (TTree*)f->Get("ZeeTiming");   //!pointer to the analyzed TTree or TChain

	fChain->SetBranchAddress("t1calib_seed", &t1calib_seed);
	fChain->SetBranchAddress("t2calib_seed", &t2calib_seed);
	fChain->SetBranchAddress("ele1SeedIEta", &ele1SeedIEta);
	fChain->SetBranchAddress("ele2SeedIEta", &ele2SeedIEta);
	fChain->SetBranchAddress("ele1SeedIPhi", &ele1SeedIPhi);
	fChain->SetBranchAddress("ele2SeedIPhi", &ele2SeedIPhi);
	fChain->SetBranchAddress("ele1IsEB", &ele1IsEB);
	fChain->SetBranchAddress("ele2IsEB", &ele2IsEB);
	fChain->SetBranchAddress("run", &run);
	fChain->SetBranchAddress("ele1Pt", &ele1Pt);
	fChain->SetBranchAddress("ele2Pt", &ele2Pt);

 	//makes histogram
 	TH1F *histTRC_12;
	histTRC_12 = new TH1F("histTRC_12","; Time [ns];Number of Events", 100, -3,3);

	TH1F *histE_12;
	histE_12 = new TH1F("histE_12","; #Delta #eta;Number of Events", 100, -200,200);

	TH1F *histP_12;
	histP_12 = new TH1F("histP_12","; #Delta #phi;Number of Events", 100, 0, 450);

	TH2F *histRunT;
	histRunT = new TH2F("histRunT","; Run Number ;#Deltat ;Number of Events", 100, 284038, 284045, 100, -3 , 3 );
	
	TH2F *histDT_PT1;
	histDT_PT1 = new TH2F("histDT_PT1","; Pt1 Value ;#Deltat ;PT value", 100, 0, 3000, 100, -3 , 3 );

	TH2F *histDT_PT2;
	histDT_PT2 = new TH2F("histDT_PT2","; Pt2 Value ;#Deltat ;PT value", 100, 0, 3000, 100, -3 , 3 );

	
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
		fChain->GetEntry(jentry);

		// for the t1calib_seed - t2calib_seed time resolution determination
		histTRC_12->Fill( t1calib_seed - t2calib_seed ); 

		deltaT = t1calib_seed - t2calib_seed;

		histRunT->Fill( run, deltaT );

		histDT_PT1->Fill( ele1Pt, deltaT );
		histDT_PT2->Fill( ele2Pt, deltaT );

		// for the ele1SeediEta - ele2SeediEta time resolution determination
		if (  !(ele1IsEB==1 && ele2IsEB==1) ) continue;
		histE_12->Fill( ele1SeedIEta - ele2SeedIEta );  
		histP_12->Fill( abs(ele1SeedIPhi - ele2SeedIPhi) );
	}

	histRunT->FitSlicesY(0, lowR, highR, 1);


	// for the t1calib_seed - t2calib_seed time resolution determination
	TFile *file = new TFile("output_Zee_test.root", "RECREATE");
	file->cd();
	file->WriteTObject(histTRC_12,"histTRC_12", "WriteDelete");

	TCanvas *c1 = new TCanvas ("c1","c1",800, 600); 
	TLatex *tex1 = new TLatex();
	c1->cd();
	TF1* f1_g1;

	double mean1 = histTRC_12->GetMean();
	double rms1 = histTRC_12->GetRMS();
	double xmin1 = mean1 - 2.0*rms1;
	double xmax1 = mean1 + 2.0*rms1;
	f1_g1 = new TF1( Form("g_fit"), "gaus(0)", xmin1, xmax1);
	histTRC_12->Fit(Form("g_fit"),"QMLES","",xmin1,xmax1);

	histTRC_12->GetXaxis()->SetTitleSize( 0.045 );
	histTRC_12->GetXaxis()->SetTitleOffset( 1 );
	histTRC_12->GetYaxis()->SetTitleSize( 0.045 );
	histTRC_12->GetYaxis()->SetTitleOffset( 1.1 ); 

	histTRC_12->Draw();
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	tex1->DrawLatex(0.8, 150, Form("#sigma = %.1f #pm %.1f ps", 1000*f1_g1->GetParameter(2), 1000*f1_g1->GetParError(2)));
	c1->SaveAs( Form("TR_12.C") );
	c1->SaveAs( Form("TR_12.pdf") ); 

	// for the ele1SeedIEta - ele2SeedIEta difference
	file->WriteTObject(histE_12,"histE_12", "WriteDelete");

	TCanvas *c2 = new TCanvas ("c2","c2",800, 600); 
	TLatex *tex2 = new TLatex();
	c2->cd();
	TF1* f1_g2;

	double mean2 = histE_12->GetMean();
	double rms2 = histE_12->GetRMS();
	double xmin2 = mean2 - 2.0*rms2;
	double xmax2 = mean2 + 2.0*rms2;
	f1_g2 = new TF1( Form("g_fit"), "gaus(0)", xmin2, xmax2);
	histE_12->Fit(Form("g_fit"),"QMLES","",xmin2,xmax2);

	histE_12->GetXaxis()->SetTitleSize( 0.045 );
	histE_12->GetXaxis()->SetTitleOffset( 1 );
	histE_12->GetYaxis()->SetTitleSize( 0.045 );
	histE_12->GetYaxis()->SetTitleOffset( 1.1 ); 

	histE_12->Draw();
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	tex2->DrawLatex(2, 105, Form("#sigma = %.1f #pm %.1f ps", 1000*f1_g2->GetParameter(2), 1000*f1_g2->GetParError(2)));
	c2->SaveAs( Form("E_12.C") );
	c2->SaveAs( Form("E_12.pdf") ); 

	// for the ele1SeedIPhi - ele2SeedIPhi difference
	file->WriteTObject(histP_12,"histP_12", "WriteDelete");

	TCanvas *c3 = new TCanvas ("c3","c3",800, 600); 
	TLatex *tex3 = new TLatex();
	c3->cd();
	TF1* f1_g3;

	double mean3 = histP_12->GetMean();
	double rms3 = histP_12->GetRMS();
	double xmin3 = mean3 - 2.0*rms3;
	double xmax3 = mean3 + 2.0*rms3;
	f1_g3 = new TF1( Form("g_fit"), "gaus(0)", xmin3, xmax3);
	histP_12->Fit(Form("g_fit"),"QMLES","",xmin3,xmax3);

	histP_12->GetXaxis()->SetTitleSize( 0.045 );
	histP_12->GetXaxis()->SetTitleOffset( 1 );
	histP_12->GetYaxis()->SetTitleSize( 0.045 );
	histP_12->GetYaxis()->SetTitleOffset( 1.1 ); 

	histP_12->Draw();
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

	tex3->DrawLatex(225, 250, Form("#sigma = %.1f #pm %.1f ps", 1000*f1_g3->GetParameter(2), 1000*f1_g3->GetParError(2)));
	c3->SaveAs( Form("P_12.C") );
	c3->SaveAs( Form("P_12.pdf") ); 

	// for the fitSlices for time vs run number
	file->WriteTObject(histRunT,"histRunT", "WriteDelete");

	TCanvas *c4 = new TCanvas ("c4","c4",800, 600); 
	TLatex *tex4 = new TLatex();
	c4->cd();

	histRunT->GetXaxis()->SetTitleSize( 0.045 );
	histRunT->GetXaxis()->SetTitleOffset( 1 );
	histRunT->GetYaxis()->SetTitleSize( 0.045 );
	histRunT->GetYaxis()->SetTitleOffset( 1.1 ); 

	histRunT->Draw();

	c4->SaveAs( Form("histRunT.C") );
	c4->SaveAs( Form("histRunT.pdf") );

}
