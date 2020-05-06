#include <TGraph.h>



void checks_trackdT_PV()

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

  const int n_density = 7;
  float max_pu_density = 200.0*TMath::Gaus(beamSpotZ,beamSpotZ,beamSpotSigmaZ,1);


  

  TFile f200("root://eoscms//eos/cms//store/group/phys_susy/razor/run2/Run2RazorNtupleFromAODV4.0/MC/Upgrade/zhicaiz/RelValTTbar_13/16Mar2017/razorNtuple_PU200_Timing_TTBar_Iso.root");
  TTree *tree200 = (TTree*)f200.Get("ntuples/RazorEvents");




  TH1F *trackdT_lowPUdensity = new TH1F("trackdT_lowPUdensity","pvZ_losIso_NoTiming",100,-0.1,0.1);
  TH1F *trackdT_highPUdensity = new TH1F("trackdT_highPUdensity","pvZ_losIso_NoTiming",100,-0.1,0.1);

  tree200->Draw("allTrackdT_GenMatchPV>>trackdT_lowPUdensity","abs(pvZ_New_dzt)>6.5 && !isFakePV");
  tree200->Draw("allTrackdT_GenMatchPV>>trackdT_highPUdensity","abs(pvZ_New_dzt)<1.0 && !isFakePV");


  trackdT_lowPUdensity->Scale(1.0/trackdT_lowPUdensity->GetEntries());
  trackdT_highPUdensity->Scale(1.0/trackdT_highPUdensity->GetEntries());

  gStyle->SetOptStat(0);


  TCanvas *cv = 0;
  TLegend *legend_h = 0;

  cv = new TCanvas("cv","cv", 800,800);
  cv->SetLeftMargin(0.15);
  cv->SetBottomMargin(0.12);

//  cv->SetLogy();
	
  legend_h = new TLegend(0.20,0.75,0.60,0.88, NULL,"brNDC");
  legend_h->SetTextSize(0.03);
  legend_h->SetBorderSize(0);
  legend_h->SetLineColor(1);
  legend_h->SetLineStyle(1);
  legend_h->SetLineWidth(1);
  legend_h->SetFillColor(0);
  legend_h->SetFillStyle(1001);

  legend_h->AddEntry(trackdT_lowPUdensity, "low pu density - abs(pvZ) > 6.5");
  legend_h->AddEntry(trackdT_highPUdensity, "high pu density - abs(pvZ) < 1.0");

  trackdT_lowPUdensity->Draw();
  trackdT_highPUdensity->Draw("same");
  legend_h->Draw();

  trackdT_lowPUdensity->SetLineColor(2);
  trackdT_highPUdensity->SetLineColor(4);

  trackdT_lowPUdensity->SetTitle("");
  trackdT_lowPUdensity->GetXaxis()->SetTitle("dT (track, vertex) / ps");
  trackdT_lowPUdensity->GetXaxis()->SetTitleSize(0.045);
  trackdT_lowPUdensity->GetXaxis()->SetTitleOffset(1.1);
  trackdT_lowPUdensity->GetYaxis()->SetTitle("Events");
  trackdT_lowPUdensity->GetYaxis()->SetTitleOffset(1.2);
  trackdT_lowPUdensity->GetYaxis()->SetTitleSize(0.045);
  trackdT_lowPUdensity->GetYaxis()->SetRangeUser(0.,0.05);

  cv->SaveAs("~/www/sharebox/tomyself/Timing/h_trackdT_differentPUdensity_withValidTiming.pdf");
  cv->SaveAs("~/www/sharebox/tomyself/Timing/h_trackdT_differentPUdensity_withValidTiming.png");
  cv->SaveAs("~/www/sharebox/tomyself/Timing/h_trackdT_differentPUdensity_withValidTiming.C");


}
