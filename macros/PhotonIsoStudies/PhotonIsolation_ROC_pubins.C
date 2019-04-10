#include <TGraph.h>


void PhotonIsolation_ROC_pubins()
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
  const float NPU = 200.0;


  const float bin1_low = 0.2;
  const float bin1_high = 0.6;

  const float bin2_low = 1.1;
  const float bin2_high = 1.4;



  TFile f("../../PhotonNtupler_Hgg_UpgradeTiming_PU200_Timing_Iso.root");
  TTree *tree = (TTree*)f.Get("Photons");
  int NEntries = tree->GetEntries();

  TFile f_TTBar("../../PhotonNtupler_TTBar_UpgradeTiming_PU200_Timing_Iso.root");
  TTree *tree_TTBar = (TTree*)f_TTBar.Get("Photons");
  int NEntries_TTBar = tree_TTBar->GetEntries();


cout<<"NEntries "<<NEntries<<"   "<<NEntries_TTBar<<endl;

  Float_t         pvZ;
  Float_t         PhoSigmaIetaIeta;
  Float_t         PhoEta;
  Float_t         PhoPt;
  Float_t         PhoHOverE;
  Bool_t          PhoPassEleVeto;
  
  Float_t         PhoChargedHadronIso_NewPV;
 
  
  tree->SetBranchAddress("pvZ_New", &pvZ);
  tree->SetBranchAddress("PhoSigmaIetaIeta", &PhoSigmaIetaIeta);
  tree->SetBranchAddress("PhoEta", &PhoEta);
  tree->SetBranchAddress("PhoPt", &PhoPt);
  tree->SetBranchAddress("PhoHOverE", &PhoHOverE);
  tree->SetBranchAddress("PhoPassEleVeto", &PhoPassEleVeto);
  
  tree->SetBranchAddress("PhoChargedHadronIso_NewPV_Timing80_TrkVtx", &PhoChargedHadronIso_NewPV);


  tree_TTBar->SetBranchAddress("pvZ_New", &pvZ);
  tree_TTBar->SetBranchAddress("PhoSigmaIetaIeta", &PhoSigmaIetaIeta);
  tree_TTBar->SetBranchAddress("PhoEta", &PhoEta);
  tree_TTBar->SetBranchAddress("PhoPt", &PhoPt);
  tree_TTBar->SetBranchAddress("PhoHOverE", &PhoHOverE);
  tree_TTBar->SetBranchAddress("PhoPassEleVeto", &PhoPassEleVeto);
  
  tree_TTBar->SetBranchAddress("PhoChargedHadronIso_NewPV_Timing80_TrkVtx", &PhoChargedHadronIso_NewPV);




  const int N_IsoCut = 300;

  const float IsoCut_Low = 0.0;
  const float IsoCut_High = 4.5;
  float Eff_Sig_bin2[N_IsoCut] = {0.0};
  float Eff_Bkg_bin2[N_IsoCut] = {0.0};


  float Eff_Sig_bin1[N_IsoCut] = {0.0};
  float Eff_Bkg_bin1[N_IsoCut] = {0.0};


  int N_Pho_Total_NoBin_bin1 = 0;
  int N_Pho_Total_NoBin_bin2 = 0;
  int N_Pho_PassIso_bin2 = 0;
  int N_Pho_PassIso_bin1 = 0;


 for(int idx_iso = 0; idx_iso < N_IsoCut; idx_iso++)
{

	N_Pho_Total_NoBin_bin1 = 0;
	N_Pho_Total_NoBin_bin2 = 0;
  	N_Pho_PassIso_bin2 = 0;
  	N_Pho_PassIso_bin1 = 0;
	float IsoCut_thisBin = IsoCut_Low + idx_iso*1.0*(IsoCut_High-IsoCut_Low)*1.0*(idx_iso)*0.01/(N_IsoCut);
	//float IsoCut_thisBin = IsoCut_Low + idx_iso*1.0*(IsoCut_High-IsoCut_Low)*1.0/(N_IsoCut);

  for(int i=0;i<NEntries;i++)
 { 

	tree->GetEntry(i);

	float pu_density =  200*TMath::Gaus(pvZ*10.0,beamSpotZ,beamSpotSigmaZ,1);

        bool passCut = false;
	//if(/*pu_density>1.2 && pu_density <1.5 && */abs(PhoEta)<1.47 && PhoSigmaIetaIeta<SigIEtaIEtaCut_LWP_EB && PhoHOverE<HoverECut_LWP_EB && PhoPt> 30.0 && PhoPassEleVeto==1) passCut = true;
	if(/*pu_density>1.2 && pu_density <1.5 && */abs(PhoEta)<1.47 && PhoPt>30.0 ) passCut = true;
	//else if(abs(PhoEta)<2.5 && PhoSigmaIetaIeta<SigIEtaIEtaCut_LWP_EE && PhoHOverE<HoverECut_LWP_EE && PhoPt> 40.0 && PhoPassEleVeto==1) passCut = true;
	if(!passCut) continue;

	if(pu_density<bin1_high && pu_density>bin1_low) N_Pho_Total_NoBin_bin1 ++;		
	if(pu_density<bin2_high && pu_density>bin2_low) N_Pho_Total_NoBin_bin2 ++;		

	if( PhoChargedHadronIso_NewPV < IsoCut_thisBin && pu_density<bin1_high && pu_density>bin1_low) N_Pho_PassIso_bin1 ++;
	if( PhoChargedHadronIso_NewPV < IsoCut_thisBin && pu_density<bin2_high && pu_density>bin2_low) N_Pho_PassIso_bin2 ++;
 }
	
	if(N_Pho_Total_NoBin_bin1>0) 
	{
	Eff_Sig_bin1[idx_iso] = 100.*(N_Pho_PassIso_bin1*1.0)/(N_Pho_Total_NoBin_bin1*1.0);
	}
	
	if(N_Pho_Total_NoBin_bin2>0) 
	{
	Eff_Sig_bin2[idx_iso] = 100.*(N_Pho_PassIso_bin2*1.0)/(N_Pho_Total_NoBin_bin2*1.0);
	}
	
	N_Pho_Total_NoBin_bin1 = 0;
	N_Pho_Total_NoBin_bin2 = 0;
  	N_Pho_PassIso_bin2 = 0;
  	N_Pho_PassIso_bin1 = 0;

  for(int i=0;i<NEntries_TTBar;i++)
 { 

	tree_TTBar->GetEntry(i);
	float pu_density =  200*TMath::Gaus(pvZ*10.0,beamSpotZ,beamSpotSigmaZ,1);
        bool passCut = false;
	//if( /*pu_density>1.2 && pu_density <1.5  && */abs(PhoEta)<1.47 && PhoSigmaIetaIeta<SigIEtaIEtaCut_LWP_EB && PhoHOverE<HoverECut_LWP_EB && PhoPt> 30.0 && PhoPassEleVeto==1) passCut = true;
	if(/*pu_density>1.2 && pu_density <1.5 && */abs(PhoEta)<1.47 && PhoPt>30.0 ) passCut = true;
	//else if(abs(PhoEta)<2.5 && PhoSigmaIetaIeta<SigIEtaIEtaCut_LWP_EE && PhoHOverE<HoverECut_LWP_EE && PhoPt> 40.0 && PhoPassEleVeto==1) passCut = true;
	if(!passCut) continue;

	if(pu_density<bin1_high && pu_density>bin1_low) N_Pho_Total_NoBin_bin1 ++;		
	if(pu_density<bin2_high && pu_density>bin2_low) N_Pho_Total_NoBin_bin2 ++;		

	if( PhoChargedHadronIso_NewPV < IsoCut_thisBin && pu_density<bin1_high && pu_density>bin1_low) N_Pho_PassIso_bin1 ++;
	if( PhoChargedHadronIso_NewPV < IsoCut_thisBin && pu_density<bin2_high && pu_density>bin2_low) N_Pho_PassIso_bin2 ++;

 }
	if(N_Pho_Total_NoBin_bin1>0) 
	{
	Eff_Bkg_bin1[idx_iso] = 100.*(N_Pho_PassIso_bin1*1.0)/(N_Pho_Total_NoBin_bin1*1.0);
	}
	
	if(N_Pho_Total_NoBin_bin2>0) 
	{
	Eff_Bkg_bin2[idx_iso] = 100.*(N_Pho_PassIso_bin2*1.0)/(N_Pho_Total_NoBin_bin2*1.0);
	}
	

	cout<<"idx_iso: "<<idx_iso <<"  iso_cut = "<<IsoCut_thisBin<<"  Bkg_Eff (1/2) = "<<Eff_Bkg_bin2[idx_iso]<<"  "<<Eff_Bkg_bin1[idx_iso]<<"  Sig_Eff (1/2) = "<<Eff_Sig_bin2[idx_iso]<<"  "<<Eff_Sig_bin1[idx_iso]<<endl;
//	cout<<"Iso cut: "<< IsoCut_thisBin<<"  Eff_Sig_bin2 = "<<Eff_Sig_bin2[idx_iso]<<"  Eff_Bkg_bin2 = "<<Eff_Bkg_bin2[idx_iso]<<endl;
}


	TCanvas * myC = new TCanvas("c1","c1",100,100,700,700);
        myC->cd();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);

        TGraph *gr_1 = new TGraph(N_IsoCut,Eff_Bkg_bin2,Eff_Sig_bin2);
        gr_1->SetLineWidth(3);
        gr_1->SetLineColor( 4 );
        gr_1->GetXaxis()->SetTitle("Fake photon iso Eff (%)");
        gr_1->GetYaxis()->SetTitle("Hgg photon iso Eff (%)");
        gr_1->GetXaxis()->SetRangeUser(0.,10.);
        gr_1->GetYaxis()->SetRangeUser(50.0,110.);
        gr_1->SetTitle("");
        gr_1->GetYaxis()->SetTitleOffset(1.0);
        gr_1->GetYaxis()->SetTitleSize(0.06);
        gr_1->GetYaxis()->SetLabelSize(0.042);
        gr_1->GetYaxis()->SetLabelOffset(0.0);
        gr_1->GetXaxis()->SetTitleOffset(1.0);
        gr_1->GetXaxis()->SetTitleSize(0.06);
        gr_1->GetXaxis()->SetLabelSize(0.042);
        gr_1->GetXaxis()->SetLabelOffset(0.0);
        gr_1->Draw("AL");

        TGraph *gr_2 = new TGraph(N_IsoCut,Eff_Bkg_bin1,Eff_Sig_bin1);
        gr_2->SetLineColor( 2 );
        gr_2->SetLineWidth( 3 );
        gr_2->Draw("same");


        TLegend *leg = new TLegend(0.20,0.78,0.60,0.88, NULL,"brNDC"); 
        leg->SetBorderSize(0);
        leg->SetTextSize(0.028);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->SetFillStyle(1001);
        leg->AddEntry(gr_2, "0.2 < pu density < 0.6 (Timing80_TrkVtx)" ,"l");
        leg->AddEntry(gr_1, "1.1 < pu density < 1.4 (Timing80_TrkVtx)" ,"l");
        leg->Draw();

 	myC->SaveAs("~/www/sharebox/tomyself/Timing/PhotonIsolationEff_ROC_pu_bins_Timing80_TrkVtx.pdf");
        myC->SaveAs("~/www/sharebox/tomyself/Timing/PhotonIsolationEff_ROC_pu_bins_Timing80_TrkVtx.C");
        myC->SaveAs("~/www/sharebox/tomyself/Timing/PhotonIsolationEff_ROC_pu_bins_Timing80_TrkVtx.png");

}
