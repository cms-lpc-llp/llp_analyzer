void doSplit( string process ) {

  TFile *inf = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_Uncorrected.root","READ");
  string outfilename = "";
  if (process == "TTJets") outfilename = "RazorScaleFactors_Inclusive_TTJets.root";
  if (process == "WJets") outfilename = "RazorScaleFactors_Inclusive_WJets.root";
  if (process == "WJetsInv") outfilename = "RazorScaleFactors_Inclusive_WJetsInv.root";
  TFile *outf = new TFile(outfilename.c_str(),"RECREATE");

  TH2Poly *ttbarNominal = (TH2Poly*)inf->Get("TTJetsScaleFactors");
  TH2Poly *ttbarUp = (TH2Poly*)inf->Get("TTJetsScaleFactorsUp");
  TH2Poly *ttbarDown = (TH2Poly*)inf->Get("TTJetsScaleFactorsDown");
  TH2Poly *wNominal = (TH2Poly*)inf->Get("WJetsScaleFactors");
  TH2Poly *wUp = (TH2Poly*)inf->Get("WJetsScaleFactorsUp");
  TH2Poly *wDown = (TH2Poly*)inf->Get("WJetsScaleFactorsDown");
  TH2Poly *wInvNominal = (TH2Poly*)inf->Get("WJetsInvScaleFactors");
  TH2Poly *wInvUp = (TH2Poly*)inf->Get("WJetsInvScaleFactorsUp");
  TH2Poly *wInvDown = (TH2Poly*)inf->Get("WJetsInvScaleFactorsDown");

  if (process == "WJets") {
    outf->WriteTObject(wNominal,"WJetsScaleFactors");
    outf->WriteTObject(wUp,"WJetsScaleFactorsUp");
    outf->WriteTObject(wDown,"WJetsScaleFactorsDown");
  }
  if (process == "TTJets") {
    outf->WriteTObject(ttbarNominal,"TTJetsScaleFactors");
    outf->WriteTObject(ttbarUp,"TTJetsScaleFactorsUp");
    outf->WriteTObject(ttbarDown,"TTJetsScaleFactorsDown");
  }
  if (process == "WJetsInv") {  
    outf->WriteTObject(wInvNominal,"WJetsInvScaleFactors");
    outf->WriteTObject(wInvUp,"WJetsInvScaleFactorsUp");
    outf->WriteTObject(wInvDown,"WJetsInvScaleFactorsDown");
  }
  
  outf->Close();
  inf->Close();
}


void doCopy() {

  TFile *inf = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_TTJets.root","READ");
  TFile *inf2 = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_WJets.root","READ");
  TFile *inf3 = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_WJetsInv.root","READ");
  TFile *inf4 = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_GJetsInv.root","READ");
  TFile *outf = new TFile("RazorScaleFactors_Inclusive_Uncorrected.root","RECREATE");

  TH2Poly *ttbarNominal = (TH2Poly*)inf->Get("TTJetsScaleFactors");
  TH2Poly *ttbarUp = (TH2Poly*)inf->Get("TTJetsScaleFactorsUp");
  TH2Poly *ttbarDown = (TH2Poly*)inf->Get("TTJetsScaleFactorsDown");
  TH2Poly *wNominal = (TH2Poly*)inf2->Get("WJetsScaleFactors");
  TH2Poly *wUp = (TH2Poly*)inf2->Get("WJetsScaleFactorsUp");
  TH2Poly *wDown = (TH2Poly*)inf2->Get("WJetsScaleFactorsDown");
  TH2Poly *wInvNominal = (TH2Poly*)inf3->Get("WJetsInvScaleFactors");
  TH2Poly *wInvUp = (TH2Poly*)inf3->Get("WJetsInvScaleFactorsUp");
  TH2Poly *wInvDown = (TH2Poly*)inf3->Get("WJetsInvScaleFactorsDown");
  TH2Poly *GJetInvNominal = (TH2Poly*)inf4->Get("GJetsInvScaleFactors");

  outf->WriteTObject(wNominal,"WJetsScaleFactors");
  outf->WriteTObject(wUp,"WJetsScaleFactorsUp");
  outf->WriteTObject(wDown,"WJetsScaleFactorsDown");
  outf->WriteTObject(ttbarNominal,"TTJetsScaleFactors");
  outf->WriteTObject(ttbarUp,"TTJetsScaleFactorsUp");
  outf->WriteTObject(ttbarDown,"TTJetsScaleFactorsDown");
  outf->WriteTObject(wInvNominal,"WJetsInvScaleFactors");
  outf->WriteTObject(wInvUp,"WJetsInvScaleFactorsUp");
  outf->WriteTObject(wInvDown,"WJetsInvScaleFactorsDown");
  outf->WriteTObject(GJetInvNominal,"GJetsInvScaleFactors");
  outf->Close();
  inf->Close();

}

void doModify( TH2Poly *f, double correction, double error) {
  
    for (int nb = 1; nb < f->GetNumberOfBins()+1; ++nb) {
      double val = correction*f->GetBinContent(nb);
      double err = sqrt( pow( correction*f->GetBinError(nb), 2) + pow( error*f->GetBinContent(nb), 2) );
      f->SetBinContent(nb, val);
      f->SetBinError(nb-1, err); //work around TH2Poly bug (sets error for the wrong bin)
    }
}

void addMultiJetCorrection() {

  TFile *inf = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_Uncorrected.root","READ");
  assert(inf);
  TFile *outf = new TFile("RazorScaleFactors_Inclusive_CorrectedToMultiJet.root","RECREATE");

  TH2Poly *ttbarNominal = (TH2Poly*)inf->Get("TTJetsScaleFactors");
  TH2Poly *ttbarUp = (TH2Poly*)inf->Get("TTJetsScaleFactorsUp");
  TH2Poly *ttbarDown = (TH2Poly*)inf->Get("TTJetsScaleFactorsDown");
  TH2Poly *wNominal = (TH2Poly*)inf->Get("WJetsScaleFactors");
  TH2Poly *wUp = (TH2Poly*)inf->Get("WJetsScaleFactorsUp");
  TH2Poly *wDown = (TH2Poly*)inf->Get("WJetsScaleFactorsDown");
  TH2Poly *wInvNominal = (TH2Poly*)inf->Get("WJetsInvScaleFactors");
  TH2Poly *wInvUp = (TH2Poly*)inf->Get("WJetsInvScaleFactorsUp");
  TH2Poly *wInvDown = (TH2Poly*)inf->Get("WJetsInvScaleFactorsDown");
  TH2Poly *GJetsInvNominal = (TH2Poly*)inf->Get("GJetsInvScaleFactors");

  ttbarNominal->SetName("TTJetsScaleFactors");
  ttbarUp->SetName("TTJetsScaleFactorsUp");
  ttbarDown->SetName("TTJetsScaleFactorsDown");
  wNominal->SetName("WJetsScaleFactors");
  wUp->SetName("WJetsScaleFactorsUp");
  wDown->SetName("WJetsScaleFactorsDown");
  wInvNominal->SetName("WJetsInvScaleFactors");
  wInvUp->SetName("WJetsScaleInvFactorsUp");
  wInvDown->SetName("WJetsInvScaleFactorsDown");

  doModify(ttbarNominal, 0.9, 0.1);
  doModify(ttbarUp, 0.9, 0.1);
  doModify(ttbarDown, 0.9, 0.1);
  doModify(wNominal, 0.9, 0.1);
  doModify(wUp, 0.9, 0.1);
  doModify(wDown, 0.9, 0.1);
  doModify(wInvNominal, 0.9, 0.1);
  doModify(wInvUp, 0.9, 0.1);
  doModify(wInvDown, 0.9, 0.1);
  doModify(GJetsInvNominal, 0.87, 0.05);

  outf->WriteTObject(ttbarNominal);
  outf->WriteTObject(ttbarUp);
  outf->WriteTObject(ttbarDown);
  outf->WriteTObject(wNominal);
  outf->WriteTObject(wUp);
  outf->WriteTObject(wDown);
  outf->WriteTObject(wInvNominal);
  outf->WriteTObject(wInvUp);
  outf->WriteTObject(wInvDown);
  outf->WriteTObject(GJetsInvNominal);
  outf->Close();
  inf->Close();

}

void addSevenJetCorrection() {

  TFile *inf = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_Uncorrected.root","READ");
  TFile *inf2 = new TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_GJetsInv.root","READ");
  assert(inf);
  assert(inf2);
  TFile *outf = new TFile("RazorScaleFactors_Inclusive_CorrectedTo7Jet.root","RECREATE");

  TH2Poly *ttbarNominal = (TH2Poly*)inf->Get("TTJetsScaleFactors");
  TH2Poly *ttbarUp = (TH2Poly*)inf->Get("TTJetsScaleFactorsUp");
  TH2Poly *ttbarDown = (TH2Poly*)inf->Get("TTJetsScaleFactorsDown");
  TH2Poly *wNominal = (TH2Poly*)inf->Get("WJetsScaleFactors");
  TH2Poly *wUp = (TH2Poly*)inf->Get("WJetsScaleFactorsUp");
  TH2Poly *wDown = (TH2Poly*)inf->Get("WJetsScaleFactorsDown");
  TH2Poly *wInvNominal = (TH2Poly*)inf->Get("WJetsInvScaleFactors");
  TH2Poly *wInvUp = (TH2Poly*)inf->Get("WJetsInvScaleFactorsUp");
  TH2Poly *wInvDown = (TH2Poly*)inf->Get("WJetsInvScaleFactorsDown");
  TH2Poly *GJetsInvNominal = (TH2Poly*)inf2->Get("GJetsInvScaleFactors");

  ttbarNominal->SetName("TTJetsScaleFactors");
  ttbarUp->SetName("TTJetsScaleFactorsUp");
  ttbarDown->SetName("TTJetsScaleFactorsDown");
  wNominal->SetName("WJetsScaleFactors");
  wUp->SetName("WJetsScaleFactorsUp");
  wDown->SetName("WJetsScaleFactorsDown");
  wInvNominal->SetName("WJetsInvScaleFactors");
  wInvUp->SetName("WJetsScaleInvFactorsUp");
  wInvDown->SetName("WJetsInvScaleFactorsDown");

  doModify(ttbarNominal, 0.55, 0.1);
  doModify(ttbarUp, 0.55, 0.1);
  doModify(ttbarDown, 0.55, 0.1);
  doModify(wNominal, 0.55, 0.1);
  doModify(wUp, 0.55, 0.1);
  doModify(wDown, 0.55, 0.1);
  doModify(wInvNominal, 0.70, 0.1);
  doModify(wInvUp, 0.70, 0.1);
  doModify(wInvDown, 0.70, 0.1);
  doModify(GJetsInvNominal, 1.07, 0.05);

  outf->WriteTObject(ttbarNominal);
  outf->WriteTObject(ttbarUp);
  outf->WriteTObject(ttbarDown);
  outf->WriteTObject(wNominal);
  outf->WriteTObject(wUp);
  outf->WriteTObject(wDown);
  outf->WriteTObject(wInvNominal);
  outf->WriteTObject(wInvUp);
  outf->WriteTObject(wInvDown);
  outf->WriteTObject(GJetsInvNominal);
  outf->Close();
  inf->Close();
  inf2->Close();

}


void copyScaleFactorHistograms( int option = 0) {

  if (option == 0) {
    doSplit("TTJets");
    doSplit("WJets");
    doSplit("WJetsInv");
  }

  if (option == 1) {
    doCopy();
  }

  if (option == 2) {
    addMultiJetCorrection();
    addSevenJetCorrection();
  }

}
