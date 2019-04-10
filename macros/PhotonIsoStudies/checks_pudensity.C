#include <TGraph.h>



void checks_pudensity()

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

  const int n_density = 7;
  float max_pu_density = 200.0*TMath::Gaus(beamSpotZ,beamSpotZ,beamSpotSigmaZ,1);


  

  TFile f200("../../PhotonNtupler_TTBar_UpgradeTiming_PU200_Timing_Iso.root");
  TTree *tree200 = (TTree*)f200.Get("Photons");


  TH1F *pu_density_lowIso_NoTiming = new TH1F("pu_density_lowIso_NoTiming","pu_density_losIso_NoTiming",100,0.,2.0);
  TH1F *pu_density_allIso_NoTiming = new TH1F("pu_density_allIso_NoTiming","pu_density_allIso_NoTiming",100,0.,2.0);

  TH1F *pu_density_lowIso_Timing80_TrkVtx = new TH1F("pu_density_lowIso_Timing80_TrkVtx","pu_density_losIso_Timing80_TrkVtx",100,0.,2.0);
  TH1F *pu_density_allIso_Timing80_TrkVtx = new TH1F("pu_density_allIso_Timing80_TrkVtx","pu_density_allIso_Timing80_TrkVtx",100,0.,2.0);


  tree200->Draw("(200.0*TMath::Gaus(pvZ_New*10.0,0.282329,43.131299,1))>>pu_density_lowIso_NoTiming","PhoChargedHadronIso_NewPV_NoTiming<1.3 && PhoPt>30.0 && abs(PhoEta)<1.47 ");
  tree200->Draw("(200.0*TMath::Gaus(pvZ_New*10.0,0.282329,43.131299,1))>>pu_density_allIso_NoTiming","PhoChargedHadronIso_NewPV_NoTiming>-1. && PhoPt>30.0 && abs(PhoEta)<1.47 ");


  tree200->Draw("(200.0*TMath::Gaus(pvZ_New*10.0,0.282329,43.131299,1))>>pu_density_lowIso_Timing80_TrkVtx","PhoChargedHadronIso_NewPV_Timing80_TrkVtx<1.3 && PhoPt>30.0 && abs(PhoEta)<1.47 ");
  tree200->Draw("(200.0*TMath::Gaus(pvZ_New*10.0,0.282329,43.131299,1))>>pu_density_allIso_Timing80_TrkVtx","PhoChargedHadronIso_NewPV_Timing80_TrkVtx>-1. && PhoPt>30.0 && abs(PhoEta)<1.47 ");

  gStyle->SetOptStat(0);



  TCanvas *cv = 0;
  TLegend *legend_h = 0;
  TLegend *legend_h_Timing = 0;


  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

  legend_h = new TLegend(0.20,0.75,0.60,0.88, NULL,"brNDC");
  legend_h->SetTextSize(0.03);
  legend_h->SetBorderSize(0);
  legend_h->SetLineColor(1);
  legend_h->SetLineStyle(1);
  legend_h->SetLineWidth(1);
  legend_h->SetFillColor(0);
  legend_h->SetFillStyle(1001);

  legend_h->AddEntry(pu_density_lowIso_NoTiming, "iso sum pT < 1.3GeV");
  legend_h->AddEntry(pu_density_allIso_NoTiming, "all iso sum pT");

  pu_density_lowIso_NoTiming->Draw();
  pu_density_allIso_NoTiming->Draw("same");
  legend_h->Draw();
  pu_density_lowIso_NoTiming->SetLineColor(2);
  pu_density_allIso_NoTiming->SetLineColor(4);

  pu_density_lowIso_NoTiming->SetTitle("");
  pu_density_lowIso_NoTiming->GetXaxis()->SetTitle("pu density (events/mm)");
  pu_density_lowIso_NoTiming->GetXaxis()->SetTitleSize(0.045);
  pu_density_lowIso_NoTiming->GetXaxis()->SetTitleOffset(1.1);
  pu_density_lowIso_NoTiming->GetYaxis()->SetTitle("Events");
  pu_density_lowIso_NoTiming->GetYaxis()->SetTitleOffset(1.2);
  pu_density_lowIso_NoTiming->GetYaxis()->SetTitleSize(0.045);
  pu_density_lowIso_NoTiming->GetYaxis()->SetRangeUser(0,500);

  cv->SaveAs("~/www/sharebox/tomyself/Timing/h_pu_density_LowHighIso_NoTiming_TTBar.pdf");
  cv->SaveAs("~/www/sharebox/tomyself/Timing/h_pu_density_LowHighIso_NoTiming_TTBar.png");
  cv->SaveAs("~/www/sharebox/tomyself/Timing/h_pu_density_LowHighIso_NoTiming_TTBar.C");

  legend_h_Timing = new TLegend(0.20,0.75,0.60,0.88, NULL,"brNDC");
  legend_h_Timing->SetTextSize(0.03);
  legend_h_Timing->SetBorderSize(0);
  legend_h_Timing->SetLineColor(1);
  legend_h_Timing->SetLineStyle(1);
  legend_h_Timing->SetLineWidth(1);
  legend_h_Timing->SetFillColor(0);
  legend_h_Timing->SetFillStyle(1001);

  legend_h_Timing->AddEntry(pu_density_lowIso_Timing80_TrkVtx, "iso sum pT < 1.3GeV");
  legend_h_Timing->AddEntry(pu_density_allIso_Timing80_TrkVtx, "all iso sum pT");

  pu_density_lowIso_Timing80_TrkVtx->Draw();
  pu_density_allIso_Timing80_TrkVtx->Draw("same");
  legend_h_Timing->Draw();
  pu_density_lowIso_Timing80_TrkVtx->SetLineColor(2);
  pu_density_allIso_Timing80_TrkVtx->SetLineColor(4);

  pu_density_lowIso_Timing80_TrkVtx->SetTitle("");
  pu_density_lowIso_Timing80_TrkVtx->GetXaxis()->SetTitle("pu density (events/mm)");
  pu_density_lowIso_Timing80_TrkVtx->GetXaxis()->SetTitleSize(0.045);
  pu_density_lowIso_Timing80_TrkVtx->GetXaxis()->SetTitleOffset(1.1);
  pu_density_lowIso_Timing80_TrkVtx->GetYaxis()->SetTitle("Events");
  pu_density_lowIso_Timing80_TrkVtx->GetYaxis()->SetTitleOffset(1.2);
  pu_density_lowIso_Timing80_TrkVtx->GetYaxis()->SetTitleSize(0.045);
  pu_density_lowIso_Timing80_TrkVtx->GetYaxis()->SetRangeUser(0,500);

  cv->SaveAs("~/www/sharebox/tomyself/Timing/h_pu_density_LowHighIso_Timing80_TrkVtx_TTBar.pdf");
  cv->SaveAs("~/www/sharebox/tomyself/Timing/h_pu_density_LowHighIso_Timing80_TrkVtx_TTBar.png");
  cv->SaveAs("~/www/sharebox/tomyself/Timing/h_pu_density_LowHighIso_Timing80_TrkVtx_TTBar.C");


}
