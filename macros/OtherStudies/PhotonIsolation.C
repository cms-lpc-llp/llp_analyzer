#include <TGraph.h>



void PhotonIsolation()

{
  const float chargeHadronIsoCut_LWP_EB = 1.3;
  const float chargeHadronIsoCut_MWP_EB = 0.44;
  const float chargeHadronIsoCut_TWP_EB = 0.2;

  const float chargeHadronIsoCut_LWP_EE = 1.36;
  const float chargeHadronIsoCut_MWP_EE = 0.82;
  const float chargeHadronIsoCut_TWP_EE = 0.27;

  const float HoverECut_LWP_EB = 0.06;
  const float HoverECut_MWP_EB = 0.04;
  const float HoverECut_TWP_EB = 0.027;

  const float HoverECut_LWP_EE = 0.027;
  const float HoverECut_MWP_EE = 0.021;
  const float HoverECut_TWP_EE = 0.020;

  const float SigIEtaIEtaCut_LWP_EB = 0.01031;
  const float SigIEtaIEtaCut_MWP_EB = 0.01023;
  const float SigIEtaIEtaCut_TWP_EB = 0.00994;

  const float SigIEtaIEtaCut_LWP_EE = 0.0269;
  const float SigIEtaIEtaCut_MWP_EE = 0.0259;
  const float SigIEtaIEtaCut_TWP_EE = 0.0258;

  const float beamSpotZ = 0.282329;
  const float beamSpotSigmaZ = 43.131299;
  const float NPU = 140.0;

  const int n_density = 12;

  float IsoEff_NoTiming[n_density] = {0.0};// for n_density pu density bins, ranging from 0 to 2.0, bin width 0.1
  float IsoEff_Timing50_TrkVtx[n_density] = {0.0}; 
  float IsoEff_Timing80_TrkVtx[n_density] = {0.0}; 
  float IsoEff_Timing100_TrkVtx[n_density] = {0.0}; 

  float IsoEff_Timing50_TrkPho[n_density] = {0.0}; 
  float IsoEff_Timing80_TrkPho[n_density] = {0.0}; 
  float IsoEff_Timing100_TrkPho[n_density] = {0.0}; 


  int N_Pho_Total_NoBin = 0;
  int N_Pho_Total[n_density] = {0};
  int N_Pho_PassIso_NoTiming[n_density] = {0}; 
  int N_Pho_PassIso_Timing50_TrkVtx[n_density] = {0}; 
  int N_Pho_PassIso_Timing80_TrkVtx[n_density] = {0}; 
  int N_Pho_PassIso_Timing100_TrkVtx[n_density] = {0}; 
 
  int N_Pho_PassIso_Timing50_TrkPho[n_density] = {0}; 
  int N_Pho_PassIso_Timing80_TrkPho[n_density] = {0}; 
  int N_Pho_PassIso_Timing100_TrkPho[n_density] = {0}; 
 

  TTree * tree = 0;
  TFile *f  =new TFile("/afs/cern.ch/work/z/zhicaiz/public/release/CMSSW_8_1_0_pre15/src/SUSYBSMAnalysis/RazorTuplizer/python/razorNtuple_PU140_Timing_Iso.root");
  TDirectory * dir = (TDirectory*)f->Get("/afs/cern.ch/work/z/zhicaiz/public/release/CMSSW_8_1_0_pre15/src/SUSYBSMAnalysis/RazorTuplizer/python/razorNtuple_PU140_Timing_Iso.root:/ntuples");
  dir->GetObject("RazorEvents",tree);
  //cout<<tree->GetEntries()<<endl;  
  int NEntries = tree->GetEntries();

  Int_t           nPhotons;

  Float_t         pvZ;
  Float_t         phoSigmaIetaIeta[700];
  Float_t         pho_superClusterEta[700];
  Float_t         pho_superClusterEnergy[700];
  Float_t         pho_HoverE[700];
  Bool_t          pho_passEleVeto[700];
  
  Float_t         pho_sumChargedHadronPt_NewPV_NoTiming[700];
  Float_t         pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx[700];
  Float_t         pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx[700];
  Float_t         pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx[700];
 
  Float_t         pho_sumChargedHadronPt_NewPV_Timing50_TrkPho[700];
  Float_t         pho_sumChargedHadronPt_NewPV_Timing80_TrkPho[700];
  Float_t         pho_sumChargedHadronPt_NewPV_Timing100_TrkPho[700];
  
  tree->SetBranchAddress("nPhotons", &nPhotons);
  tree->SetBranchAddress("genVertexZ", &pvZ);
  tree->SetBranchAddress("phoSigmaIetaIeta", phoSigmaIetaIeta);
  tree->SetBranchAddress("pho_superClusterEta", pho_superClusterEta);
  tree->SetBranchAddress("pho_superClusterEnergy", pho_superClusterEnergy);
  tree->SetBranchAddress("pho_HoverE", pho_HoverE);
  tree->SetBranchAddress("pho_passEleVeto", pho_passEleVeto);
  
  tree->SetBranchAddress("pho_sumChargedHadronPt_NewPV_NoTiming", pho_sumChargedHadronPt_NewPV_NoTiming);
  tree->SetBranchAddress("pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx", pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx);
  tree->SetBranchAddress("pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx", pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx);
  tree->SetBranchAddress("pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx", pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx);

  tree->SetBranchAddress("pho_sumChargedHadronPt_NewPV_Timing50_TrkPho", pho_sumChargedHadronPt_NewPV_Timing50_TrkPho);
  tree->SetBranchAddress("pho_sumChargedHadronPt_NewPV_Timing80_TrkPho", pho_sumChargedHadronPt_NewPV_Timing80_TrkPho);
  tree->SetBranchAddress("pho_sumChargedHadronPt_NewPV_Timing100_TrkPho", pho_sumChargedHadronPt_NewPV_Timing100_TrkPho);


  TH1F * h_pho_sumChargedHadronPt_NewPV_NoTiming = new TH1F("h_pho_sumChargedHadronPt_NewPV_NoTiming","h_pho_sumChargedHadronPt_NewPV_NoTiming", 100,0,5);
  TH1F * h_pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx = new TH1F("h_pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx","h_pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx", 100,0,5);
  TH1F * h_pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx = new TH1F("h_pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx","h_pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx", 100,0,5);
  TH1F * h_pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx = new TH1F("h_pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx","h_pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx", 100,0,5);

  TH1F * h_pho_sumChargedHadronPt_NewPV_Timing50_TrkPho = new TH1F("h_pho_sumChargedHadronPt_NewPV_Timing50_TrkPho","h_pho_sumChargedHadronPt_NewPV_Timing50_TrkPho", 100,0,5);
  TH1F * h_pho_sumChargedHadronPt_NewPV_Timing80_TrkPho = new TH1F("h_pho_sumChargedHadronPt_NewPV_Timing80_TrkPho","h_pho_sumChargedHadronPt_NewPV_Timing80_TrkPho", 100,0,5);
  TH1F * h_pho_sumChargedHadronPt_NewPV_Timing100_TrkPho = new TH1F("h_pho_sumChargedHadronPt_NewPV_Timing100_TrkPho","h_pho_sumChargedHadronPt_NewPV_Timing100_TrkPho", 100,0,5);


  for(int i=0;i<NEntries;i++)
 { 
	tree->GetEntry(i);
	for(int j=0;j<nPhotons;j++)
	{
        bool passCut = false;
	if(abs(pho_superClusterEta[j])<1.47 && phoSigmaIetaIeta[j]<SigIEtaIEtaCut_LWP_EB && pho_HoverE[j]<HoverECut_LWP_EB && pho_superClusterEnergy[j]/cosh(pho_superClusterEta[j])> 40.0 && pho_passEleVeto[j]==1) passCut = true;
	else if(abs(pho_superClusterEta[j])<2.5 && phoSigmaIetaIeta[j]<SigIEtaIEtaCut_LWP_EE && pho_HoverE[j]<HoverECut_LWP_EE && pho_superClusterEnergy[j]/cosh(pho_superClusterEta[j])> 40.0 && pho_passEleVeto[j]==1) passCut = true;
	if(!passCut) continue;

        N_Pho_Total_NoBin ++;		

	float pu_density =  NPU*TMath::Gaus(pvZ*10.0,beamSpotZ,beamSpotSigmaZ,1);
	//cout<<"event "<<i<<"  photon "<<j<<" pu_density "<<pu_density<<endl;	
	int pu_density_bin = int(pu_density*n_density/1.2);
 	if(pu_density_bin >= n_density) continue;

  	N_Pho_Total[pu_density_bin] ++;

        if(pho_sumChargedHadronPt_NewPV_NoTiming[j]>0.0) h_pho_sumChargedHadronPt_NewPV_NoTiming->Fill(pho_sumChargedHadronPt_NewPV_NoTiming[j]);	
        if(pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx[j]>0.0) h_pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx->Fill(pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx[j]);	
        if(pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx[j]>0.0) h_pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx->Fill(pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx[j]);	
        if(pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx[j]>0.0) h_pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx->Fill(pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx[j]);	
	
 	if(pho_sumChargedHadronPt_NewPV_Timing50_TrkPho[j]>0.0) h_pho_sumChargedHadronPt_NewPV_Timing50_TrkPho->Fill(pho_sumChargedHadronPt_NewPV_Timing50_TrkPho[j]);	
        if(pho_sumChargedHadronPt_NewPV_Timing80_TrkPho[j]>0.0) h_pho_sumChargedHadronPt_NewPV_Timing80_TrkPho->Fill(pho_sumChargedHadronPt_NewPV_Timing80_TrkPho[j]);	
        if(pho_sumChargedHadronPt_NewPV_Timing100_TrkPho[j]>0.0) h_pho_sumChargedHadronPt_NewPV_Timing100_TrkPho->Fill(pho_sumChargedHadronPt_NewPV_Timing100_TrkPho[j]);	
	

	if(abs(pho_superClusterEta[j])<1.47)
	{
	if( pho_sumChargedHadronPt_NewPV_NoTiming[j] < chargeHadronIsoCut_LWP_EB) N_Pho_PassIso_NoTiming[pu_density_bin]++;
	if( pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx[j] < chargeHadronIsoCut_LWP_EB) N_Pho_PassIso_Timing50_TrkVtx[pu_density_bin]++;
	if( pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx[j] < chargeHadronIsoCut_LWP_EB) N_Pho_PassIso_Timing80_TrkVtx[pu_density_bin]++;
	if( pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx[j] < chargeHadronIsoCut_LWP_EB) N_Pho_PassIso_Timing100_TrkVtx[pu_density_bin]++;

	if( pho_sumChargedHadronPt_NewPV_Timing50_TrkPho[j] < chargeHadronIsoCut_LWP_EB) N_Pho_PassIso_Timing50_TrkPho[pu_density_bin]++;
	if( pho_sumChargedHadronPt_NewPV_Timing80_TrkPho[j] < chargeHadronIsoCut_LWP_EB) N_Pho_PassIso_Timing80_TrkPho[pu_density_bin]++;
	if( pho_sumChargedHadronPt_NewPV_Timing100_TrkPho[j] < chargeHadronIsoCut_LWP_EB) N_Pho_PassIso_Timing100_TrkPho[pu_density_bin]++;

	}
	else if(abs(pho_superClusterEta[j])<2.5)
	{
	if( pho_sumChargedHadronPt_NewPV_NoTiming[j] < chargeHadronIsoCut_LWP_EE) N_Pho_PassIso_NoTiming[pu_density_bin]++;
	if( pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx[j] < chargeHadronIsoCut_LWP_EE) N_Pho_PassIso_Timing50_TrkVtx[pu_density_bin]++;
	if( pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx[j] < chargeHadronIsoCut_LWP_EE) N_Pho_PassIso_Timing80_TrkVtx[pu_density_bin]++;
	if( pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx[j] < chargeHadronIsoCut_LWP_EE) N_Pho_PassIso_Timing100_TrkVtx[pu_density_bin]++;

	if( pho_sumChargedHadronPt_NewPV_Timing50_TrkPho[j] < chargeHadronIsoCut_LWP_EE) N_Pho_PassIso_Timing50_TrkPho[pu_density_bin]++;
	if( pho_sumChargedHadronPt_NewPV_Timing80_TrkPho[j] < chargeHadronIsoCut_LWP_EE) N_Pho_PassIso_Timing80_TrkPho[pu_density_bin]++;
	if( pho_sumChargedHadronPt_NewPV_Timing100_TrkPho[j] < chargeHadronIsoCut_LWP_EE) N_Pho_PassIso_Timing100_TrkPho[pu_density_bin]++;
	}
	}
 }

  cout<<"Total photons in sample that pass cuts: "<<N_Pho_Total_NoBin<<endl;


  float pu_density[n_density];

  for(int i=0;i<n_density;i++)
  {
	pu_density[i] = 1.2*i/n_density;
	if(N_Pho_Total[i]>0)
	{
	IsoEff_NoTiming[i] =100.0* (N_Pho_PassIso_NoTiming[i]*1.0)/(N_Pho_Total[i]);
	IsoEff_Timing50_TrkVtx[i] = 100.0* (N_Pho_PassIso_Timing50_TrkVtx[i]*1.0)/(N_Pho_Total[i]);
	IsoEff_Timing80_TrkVtx[i] = 100.0* (N_Pho_PassIso_Timing80_TrkVtx[i]*1.0)/(N_Pho_Total[i]);
	IsoEff_Timing100_TrkVtx[i] = 100.0* (N_Pho_PassIso_Timing100_TrkVtx[i]*1.0)/(N_Pho_Total[i]);

	IsoEff_Timing50_TrkPho[i] = 100.0* (N_Pho_PassIso_Timing50_TrkPho[i]*1.0)/(N_Pho_Total[i]);
	IsoEff_Timing80_TrkPho[i] = 100.0* (N_Pho_PassIso_Timing80_TrkPho[i]*1.0)/(N_Pho_Total[i]);
	IsoEff_Timing100_TrkPho[i] = 100.0* (N_Pho_PassIso_Timing100_TrkPho[i]*1.0)/(N_Pho_Total[i]);
	}
  }


  float MaxY = 0.0;
  float MinY = 999.9;

  TGraph *gr_NoTiming = new TGraph(n_density, pu_density, IsoEff_NoTiming);
  TGraph *gr_Timing100_TrkVtx = new TGraph(n_density, pu_density, IsoEff_Timing100_TrkVtx);
  TGraph *gr_Timing80_TrkVtx = new TGraph(n_density, pu_density, IsoEff_Timing80_TrkVtx);

  TGraph *gr_Timing100_TrkPho = new TGraph(n_density, pu_density, IsoEff_Timing100_TrkPho);
  TGraph *gr_Timing80_TrkPho = new TGraph(n_density, pu_density, IsoEff_Timing80_TrkPho);


  gStyle->SetOptStat(0);


  TCanvas *cv = 0;
  TLegend *legend = 0;
  TLegend *legend_h = 0;
  bool firstdrawn = false;
  TLatex *tex = 0;
  TLine *line = 0;


  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend_h = new TLegend(0.20,0.6,0.60,0.88, NULL,"brNDC");
  legend_h->SetTextSize(0.03);
  legend_h->SetBorderSize(0);
  legend_h->SetLineColor(1);
  legend_h->SetLineStyle(1);
  legend_h->SetLineWidth(1);
  legend_h->SetFillColor(0);
  legend_h->SetFillStyle(1001);

  legend_h->AddEntry(h_pho_sumChargedHadronPt_NewPV_NoTiming, "No Timing cut (PU140)");
  legend_h->AddEntry(h_pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx, "Track-Vertex Timing 80ps cut (PU140)");
  legend_h->AddEntry(h_pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx, "Track-Vertex Timing 100ps cut (PU140)");
  legend_h->AddEntry(h_pho_sumChargedHadronPt_NewPV_Timing80_TrkPho, "Track-Photon Timing 80ps cut (PU140)");
  legend_h->AddEntry(h_pho_sumChargedHadronPt_NewPV_Timing100_TrkPho, "Track-Photon Timing 100ps cut (PU140)");

 
  h_pho_sumChargedHadronPt_NewPV_NoTiming->Draw();
  h_pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx->Draw("same");
  h_pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx->Draw("same");
  h_pho_sumChargedHadronPt_NewPV_Timing80_TrkPho->Draw("same");
  h_pho_sumChargedHadronPt_NewPV_Timing100_TrkPho->Draw("same");
  legend_h->Draw();

  h_pho_sumChargedHadronPt_NewPV_NoTiming->SetLineColor(2);
  h_pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx->SetLineColor(1);
  h_pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx->SetLineColor(3);
  h_pho_sumChargedHadronPt_NewPV_Timing80_TrkPho->SetLineColor(4);
  h_pho_sumChargedHadronPt_NewPV_Timing100_TrkPho->SetLineColor(5);

  h_pho_sumChargedHadronPt_NewPV_NoTiming->SetTitle("");
  h_pho_sumChargedHadronPt_NewPV_NoTiming->GetXaxis()->SetTitle("sum charged hadron Pt");
  h_pho_sumChargedHadronPt_NewPV_NoTiming->GetXaxis()->SetTitleSize(0.045);
  h_pho_sumChargedHadronPt_NewPV_NoTiming->GetXaxis()->SetTitleOffset(1.1);
  h_pho_sumChargedHadronPt_NewPV_NoTiming->GetYaxis()->SetTitle("Events");
  h_pho_sumChargedHadronPt_NewPV_NoTiming->GetYaxis()->SetTitleOffset(1.2);
  h_pho_sumChargedHadronPt_NewPV_NoTiming->GetYaxis()->SetTitleSize(0.045);
  h_pho_sumChargedHadronPt_NewPV_NoTiming->GetYaxis()->SetRangeUser(0,200);


  cv->SaveAs("~/www/sharebox/tomyself/Timing/h_chargedHadronPt.pdf");
  cv->SaveAs("~/www/sharebox/tomyself/Timing/h_chargedHadronPt.png");
  cv->SaveAs("~/www/sharebox/tomyself/Timing/h_chargedHadronPt.C");

  legend = new TLegend(0.20,0.6,0.60,0.88, NULL,"brNDC");
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(1001);

  legend->AddEntry(gr_NoTiming, "No Timing cut (PU140)");
  legend->AddEntry(gr_Timing80_TrkVtx, "Track-Vertex Timing 80ps cut (PU140)");
  legend->AddEntry(gr_Timing100_TrkVtx, "Track-Vertex Timing 100ps cut (PU140)");
  legend->AddEntry(gr_Timing80_TrkPho, "Track-Photon Timing 80ps cut (PU140)");
  legend->AddEntry(gr_Timing100_TrkPho, "Track-Photon Timing 100ps cut (PU140)");


  gr_NoTiming->Draw("AP");
  gr_Timing80_TrkVtx->Draw("P");
  gr_Timing100_TrkVtx->Draw("P");

  gr_Timing80_TrkPho->Draw("P");
  gr_Timing100_TrkPho->Draw("P");
 
  legend->Draw();

  gr_NoTiming->SetFillStyle(0);
  gr_NoTiming->SetLineColor(kRed);
  gr_NoTiming->SetLineWidth(0);
  gr_NoTiming->SetMarkerColor(kRed);
  gr_NoTiming->SetMarkerStyle(20);
  gr_NoTiming->SetMarkerSize(2);
  gr_Timing100_TrkVtx->SetFillStyle(0);
  gr_Timing100_TrkVtx->SetLineColor(kBlue);
  gr_Timing100_TrkVtx->SetLineWidth(0);
  gr_Timing100_TrkVtx->SetMarkerColor(kBlue);
  gr_Timing100_TrkVtx->SetMarkerStyle(21);
  gr_Timing100_TrkVtx->SetMarkerSize(2);

  gr_Timing80_TrkVtx->SetFillStyle(0);
  gr_Timing80_TrkVtx->SetLineColor(kBlue);
  gr_Timing80_TrkVtx->SetLineWidth(0);
  gr_Timing80_TrkVtx->SetMarkerColor(kBlue);
  gr_Timing80_TrkVtx->SetMarkerStyle(22);
  gr_Timing80_TrkVtx->SetMarkerSize(2);


  gr_Timing100_TrkPho->SetFillStyle(0);
  gr_Timing100_TrkPho->SetLineColor(kGreen);
  gr_Timing100_TrkPho->SetLineWidth(0);
  gr_Timing100_TrkPho->SetMarkerColor(kGreen);
  gr_Timing100_TrkPho->SetMarkerStyle(21);
  gr_Timing100_TrkPho->SetMarkerSize(2);

  gr_Timing80_TrkPho->SetFillStyle(0);
  gr_Timing80_TrkPho->SetLineColor(kGreen);
  gr_Timing80_TrkPho->SetLineWidth(0);
  gr_Timing80_TrkPho->SetMarkerColor(kGreen);
  gr_Timing80_TrkPho->SetMarkerStyle(22);
  gr_Timing80_TrkPho->SetMarkerSize(2);

 
  gr_NoTiming->SetTitle("");
  gr_NoTiming->GetXaxis()->SetTitle("Linear Pileup Density (events / mm)");
  gr_NoTiming->GetXaxis()->SetTitleSize(0.045);
  gr_NoTiming->GetXaxis()->SetTitleOffset(1.1);
  gr_NoTiming->GetYaxis()->SetTitle("Photon Charged Isolation Eff (%)");
  gr_NoTiming->GetYaxis()->SetTitleOffset(1.2);
  gr_NoTiming->GetYaxis()->SetTitleSize(0.045);
  gr_NoTiming->GetYaxis()->SetRangeUser(60,120);

  cv->SaveAs("~/www/sharebox/tomyself/Timing/PhotonIsolationEff.pdf");
  cv->SaveAs("~/www/sharebox/tomyself/Timing/PhotonIsolationEff.C");
  cv->SaveAs("~/www/sharebox/tomyself/Timing/PhotonIsolationEff.png");

}
