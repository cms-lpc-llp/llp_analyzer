#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"

using namespace std;


//I haven't updated the output stuff yet.

void TriggerTurnOnHist()
{
  //For plotting stuff
  TCanvas* C = new TCanvas("C", "I hope my dimensions are right", 400, 500);
  TLegend* leg;


  //setting bin bounds for variable width histograms - bins and nBins correspond to JetPt charts
  double R2Bins[11] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 1.2};
  int nR2Bins = 10;

  double bins[9] = {80, 100, 120, 140, 160, 180, 200, 240, 1000};
  int nBins = 8;

  double MRBins[11] = {0, 100, 200, 300, 400, 600, 800, 1000, 2000, 3000, 5000};
  int nMRBins = 10;

  //declaring histograms

  TH1F* jetPt0Pass_rsq36 = new TH1F("jetPt0Pass_rsq36", "jetPt values of elements that passed the trigger, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300)>260", nBins, bins);
  TH1F* jetPt0All = new TH1F("jetPt0All", "jetPt values of elements, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300)>260", nBins, bins);
  TH1F* jetPt1Pass_rsq36 = new TH1F("jetPt1Pass_rsq36", "jetPt values of elements that passed the trigger, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300)>260", nBins, bins);
  TH1F* jetPt1All = new TH1F("jetPt1All", "jetPt values of elements, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300)>260", nBins, bins);
  TH1F* jetPt0Pass_trigger2 = new TH1F("jetPt0Pass_trigger2", "jetPt values of elements that passed the trigger, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300)>260", nBins, bins);
  TH1F* jetPt1Pass_trigger2 = new TH1F("jetPt1Pass_trigger2", "jetPt values of elements that passed the trigger, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300)>260", nBins, bins);

  TH1F* jetPt0Pass_rsq36_rsq40 = new TH1F("jetPt0Pass_rsq36_rsq40", "jetPt values of elements that passed the trigger, R^2 > 0.4", nBins, bins);
  TH1F* jetPt0All_rsq40 = new TH1F("jetPt0All_rsq40", "All jetPt values of elements, R^2 > 0.4", nBins, bins);
  TH1F* jetPt0Pass_rsq36_rsq50 = new TH1F("jetPt0Pass_rsq36_rsq50", "jetPt values of elements that passed the trigger, R^2 > 0.5", nBins, bins);
  TH1F* jetPt0All_rsq50 = new TH1F("jetPt0All_rsq50", "All jetPt values of elements, R^2 > 0.5", nBins,bins);

  TH1F* jetPt1Pass_rsq36_rsq40 = new TH1F("jetPt1Pass_rsq36_rsq40", "jetPt values of elements that passed the trigger, R^2 > 0.4", nBins, bins);
  TH1F* jetPt1All_rsq40 = new TH1F("jetPt1All_rsq40", "All jetPt values of elements, R^2 > 0.4", nBins,bins);
  TH1F* jetPt1Pass_rsq36_rsq50 = new TH1F("jetPt1Pass_rsq36_rsq50", "jetPt values of elements that passed the trigger, R^2 > 0.5", nBins, bins);
  TH1F* jetPt1All_rsq50 = new TH1F("jetPt1All_rsq50", "All jetPt values of elements, R^2 > 0.5", nBins,bins);

  TH1F* jetPt0Pass_trigger2_3const = new TH1F("jetPt0Pass_trigger2_3const", "jetPt values of elements that passed the trigger, MR > 400, R^2 > 0.25, (R^2+0.25)(MR+300)>260", nBins, bins);
  TH1F* jetPt0All_3const = new TH1F("jetPt0All_3const", "jetPt values of all elements, MR > 400, R^2 > 0.25, (R^2+0.25)(MR+300)>260", nBins, bins);
  TH1F* jetPt1Pass_trigger2_3const = new TH1F("jetPt1Pass_trigger2_3const", "jetPt values of elements that passed the trigger, MR > 400, R^2 > 0.25, (R^2+0.25)(MR+300)>260", nBins, bins);
  TH1F* jetPt1All_3const = new TH1F("jetPt1All_3const", "jetPt values of all elements, MR > 400, R^2 > 0.25, (R^2+0.25)(MR+300)>260", nBins, bins);



  TH1F* MRPass_rsq36_3const = new TH1F("MRPass_rsq36_3const", "MR values of elements that passed the trigger, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300)>260", nMRBins, MRBins);
  TH1F* MRAll_3const = new TH1F("MRAll_3const", "MR values of elements, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300)>260", nMRBins, MRBins);
  TH1F* MRPass_rsq36_2const =  new TH1F("MRPass_rsq36_2const", "MR values of elements that passed the trigger, R^2 > 0.25, (R^2+.25)(MR+300)>260", nMRBins, MRBins);
  TH1F* MRAll_2const = new TH1F("MRAll_2const", "MR values of elements, R^2 > 0.25, (R^2+.25)(MR+300)>260", nMRBins, MRBins);
  TH1F* MRPass_trigger2_3const = new TH1F("MRPass_trigger2_3const", "MR values of elements that passed the trigger, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300)>260", nMRBins, MRBins);
  TH1F* MRPass_trigger2_2const = new TH1F("MRPass_trigger2_2const", "MR values of elements that passed the trigger, R^2 > 0.25, (R^2+.25)(MR+300)>260", nMRBins, MRBins);

  
  TH1F* R2Pass_rsq36_3const = new TH1F("R2Pass_rsq36_3const", "R^2 values of elements that passed the trigger, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300)>260", nR2Bins, R2Bins);
  TH1F* R2All_3const = new TH1F("R2All_3const", "R^2 values of elements, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300)>260", nR2Bins, R2Bins);
  TH1F* R2Pass_rsq36_2const = new TH1F("R2Pass_rsq36_2const", "R^2 values of elements, MR > 200, (R^2+.25)(MR+300)>260", nR2Bins, R2Bins);
  TH1F* R2All_2const = new TH1F("R2All_2const", "R^2 values of elements, MR > 200, (R^2+.25)(MR+300)>260", nR2Bins, R2Bins);
  TH1F* R2Pass_trigger2_3const = new TH1F("R2Pass_trigger2_3const", "R^2 values of elements that passed the trigger, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300)>260", nR2Bins, R2Bins);
  TH1F* R2Pass_trigger2_2const = new TH1F("R2Pass_trigger2_2const", "R^2 values of elements, MR > 200, (R^2+.25)(MR+300)>260", nR2Bins, R2Bins);


  TH1F* jetPt_MRconstPassUneven = new TH1F("jetPt_MRconstPassUneven", "jetPt values of elements that passed the trigger, MR > 200", nBins, bins);
  TH1F* jetPt_MRconstAllUneven = new TH1F("jetPt_MRconstAllUneven", "jetPt values of elements that passed the trigger, MR > 200", nBins, bins);
  TH1F* jetPt_R2constPassUneven = new TH1F("jetPt_R2constPassUneven", "jetPt values of elements that passed the trigger, R^2 > 0.25", nBins, bins);
  TH1F* jetPt_R2constAllUneven = new TH1F("jetPtconst_R2constAllUneven", "jetPt values of elements that passed the trigger, R^2 > 0.25", nBins, bins);
  TH1F* TriggerTurnonJetPt_MRconstUneven = new TH1F("TriggerTurnonJetPt_MRconstUneven", "Trigger turnon curve for JetPt, MR > 200", nBins, bins);
  TH1F* TriggerTurnonJetPt_R2constUneven = new TH1F("TriggerTurnonJetPt_R2constUneven", "Trigger turnon curve for JetPt, R^2 > 0.25", nBins, bins);

  TH1F* R2Pass_rsq36_noconst = new TH1F("R2Pass_rsq36_noconst", "R^2 values of elements that passed the trigger", nR2Bins, R2Bins);
  TH1F* R2All_noconst = new TH1F("R2All_noconst", "R^2 values of all elements", nR2Bins, R2Bins);
  TH1F* R2Pass_rsq36_mr200 = new TH1F("R2Pass_rsq36_mr200", "R^2 values of elements that passed the trigger, MR > 200", nR2Bins, R2Bins);
  TH1F* R2All_mr200 = new TH1F("R2All_mr200", "R^2 values of all elements, MR > 200", nR2Bins, R2Bins);
  TH1F* R2Pass_rsq36_mr400 = new TH1F("R2Pass_rsq36_mr400", "R^2 values of elements that passed the trigger, MR > 400", nR2Bins, R2Bins);
  TH1F* R2All_mr400 = new TH1F("R2All_mr400", "R^2 values of all elements, MR > 400", nR2Bins, R2Bins);

  TH1F* MRPass_rsq36_rsq36 = new TH1F("MRPass_rsq36_rsq36", "MR values of elements that passed the trigger, R^2 > 0.36", nMRBins, MRBins);
  TH1F* MRAll_rsq36 = new TH1F("MRAll_rsq36", "MR values of all elements, R^2 > 0.36", nMRBins, MRBins);
  TH1F* MRPass_rsq36_rsq40 = new TH1F("MRPass_rsq36_rsq40", "MR values of elements that passed the trigger, R^2 > 0.40", nMRBins, MRBins);
  TH1F* MRAll_rsq40 = new TH1F("MRAll_rsq40", "MR values of all elements, R^2 > 0.40", nMRBins, MRBins);
  TH1F* MRPass_rsq36_rsq50 = new TH1F("MRPass_rsq36_rsq50", "MR values of elements that passed the trigger, R^2 > 0.50", nMRBins, MRBins);
  TH1F* MRAll_rsq50 = new TH1F("MRAll_rsq50", "MR values of all elements, R^2 > 0.50", nMRBins,MRBins);

  TH1F* R2Pass_trigger2_mr200hyper = new TH1F("R2Pass_trigger2_mr200hyper", "R^2 values of elements that passed the trigger, MR>200, (R^2+.25)(MR+300)>260", nR2Bins, R2Bins); 
  TH1F* R2All_mr200hyper = new TH1F("R2All_mr200hyper", "R^2 values of all elements, MR > 200, (R^2+.25)(MR+300)>260", nR2Bins, R2Bins);
  TH1F* R2Pass_trigger2_mr400hyper = new TH1F("R2Pass_trigger2_mr400hyper", "R^2 values of elements that passed the trigger, MR>400, (R^2+.25)(MR+300)>260", nR2Bins, R2Bins);
  TH1F* R2All_mr400hyper = new TH1F("R2All_mr400hyper", "R^2 values of all elements, MR > 400, (R^2+.25)(MR+300)>260", nR2Bins,R2Bins);

  TH1F* MRPass_trigger2_rsq09hyper = new TH1F("MRPass_trigger2_rsq09hyper", "MR values of elements that passed the trigger, R^2>0.09, (R^2+.25)(MR+300)>260", nMRBins, MRBins);
  TH1F* MRAll_rsq09hyper = new TH1F("MRAll_rsq09hyper", "MR values of all elements, R^2>0.09, (R^2+.25)(MR+300)>260", nMRBins,MRBins);
  TH1F* MRPass_trigger2_rsq25hyper = new TH1F("MRPass_trigger2_rsq25hyper", "MR values of elements that passed the trigger, R^2>0.25, (R^2+.25)(MR+300)>260", nMRBins, MRBins);
  TH1F* MRAll_rsq25hyper = new TH1F("MRAll_rsq25hyper", "MR values of all elements, R^2>0.25, (R^2+.25)(MR+300)>260", nMRBins,MRBins);




  TH1F* R2All = new TH1F("R2All", "All R^2 values", 100, 0, 1.2);
  TH1F* R2Pass = new TH1F("R2Pass", "R^2 values of events that passed the trigger", 100, 0, 1.2);
  TH1F* R2_MRconstAll = new TH1F("R2_MRconstAll", "R^2 Values with MR > 400", 100, 0, 1.2);
  TH1F* R2_MRconstPass = new TH1F("R2_MRconstPass", "R^2 Values of events that passed the trigger with MR > 400", 100, 0, 1.2);
  TH1F* MR_R2constAll = new TH1F("MR_R2constAll", "MR Values with R^2 > 0.25", 100, 0,5000);
  TH1F* MR_R2constPass = new TH1F("MR_R2constPass", "MR Values of events that passed the trigger with R^2 > 0.25", 100, 0, 5000);
  TH1F* combPass = new TH1F("combPass", "(MR+300)(R^2+0.25) for events that passed the trigger", 100, 0, 5000);
  TH1F* combAll = new TH1F("combAll", "(MR+300)(R^2+0.25) for all events", 100, 0, 1200);
  TH2F* bothPass = new TH2F("bothPass", "MR and R^2 values of elements that passed the trigger", 50, 0, 5000, 50, 0, 1.2);
  TH2F* bothAll = new TH2F("bothAll", "MR and R^2 values of all elements", 50, 0, 5000, 50, 0, 1.2);
  TH1F* doubleR2Pass = new TH1F("doubleR2Pass", "R^2 values of elements that passed both triggers", 100, 0, 1.2);
  TH1F* jetPtPass = new TH1F("jetPtPass", "jetPt values of elements that passed the trigger", 100, 0, 1000);
  TH1F* jetPtAll = new TH1F("jetPtAll", "jetPt values of elements", 100, 0, 1000);
  TH1F* jetPtconstPass = new TH1F("jetPtconstPass", "jetPt values of elements that passed the trigger, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300) > 260", 50, 0, 1000);
  TH1F* jetPtconstAll = new TH1F("jetPtconstAll", "jetPt values of elemnts, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300) > 260", 50, 0, 1000);
  TH1F* jetPt_MRconstPass = new TH1F("jetPt_MRconstPass", "jetPt values of elements that passed the trigger, MR > 200", 50, 0, 1000);
  TH1F* jetPt_MRconstAll = new TH1F("jetPt_MRconstAll", "jetPt values of elements, MR > 200", 50, 0, 1000);
  TH1F* jetPt_R2constPass = new TH1F("jetPt_R2constPass", "jetPt values of elements that passed the trigger, R^2 > 0.25", 50, 0, 1000);
  TH1F* jetPt_R2constAll = new TH1F("jetPt_R2constAll", "jetPt values of elements, R^2 > 0.25", 50, 0, 1000);

  TH1F* TriggerTurnonR2 = new TH1F("TriggerTurnonR2", "Trigger turnon curve for R^2, R^2 > 0.36 trigger", 100, 0, 1.2);
  TH1F* TriggerTurnonMR_R2const = new TH1F("TriggerTurnonMR_R2constraint", "Trigger turnon curve for MR, R^2 > 0.25", 100, 0, 5000);
  TH1F* TriggerTurnonR2_MRconst = new TH1F("TriggerTurnonR2_MRconstraint", "Trigger turnon curve for R^2, MR > 400", 100, 0, 1.2);
  TH1F* TriggerTurnonComb = new TH1F("TriggerTurnonComb", "Trigger turnon curve for (MR+300)*(R^2+25), R^2 > 0.25, MR > 200", 100, 0, 5000);
  TH2F* TriggerTurnonBoth = new TH2F("TriggerTurnonBoth", "Trigger turnon curve for both MR and R^2", 50, 0, 5000, 50, 0, 1.2);
  TH1F* TriggerTurnonDoubleR2 = new TH1F("TriggerTurnonDoubleR2", "Trigger turnon curve for events that passed either of 2 triggers", 100, 0, 1.2);
  TH1F* TriggerTurnonJetPt = new TH1F("TriggerTurnonJetPt", "Trigger turnon curve for JetPt", 100, 0, 1000);
  TH1F* TriggerTurnonJetPtconst = new TH1F("TriggerTurnonJetPtconst", "Trigger turnon curve for JetPt, R^2 > 0.25, MR > 200, (R^2+.25)(MR+300) > 260", 50, 0, 1000); 
  TH1F* TriggerTurnonJetPt_MRconst = new TH1F("TriggerTurnonJetPt_MRconst", "Trigger turnon curve for JetPt, MR > 200", 50, 0, 1000);
  TH1F* TriggerTurnonJetPt_R2const = new TH1F("TriggerTurnonJetPt_R2const", "Trigger turnon curve for JetPt, R^2 > 0.25", 50, 0, 1000);  

  TH1F* R2AllUneven = new TH1F("R2AllUneven", "R^2 values of all elements", nR2Bins, R2Bins);
  TH1F* R2PassUneven = new TH1F("R2AllUneven", "R^2 values of all elements that passed the trigger", nR2Bins, R2Bins);
  TH1F* TriggerTurnonR2Uneven = new TH1F("TriggerTurnonR2Uneven", "Trigger turnon curve for R^2, R^2 > 0.36 trigger", nR2Bins, R2Bins);



  //Using Sumw2() on all of the histograms
  R2All->TH1F::Sumw2();
  R2Pass->TH1F::Sumw2();
  R2_MRconstAll->TH1F::Sumw2();
  R2_MRconstPass->TH1F::Sumw2();
  MR_R2constAll->TH1F::Sumw2();
  MR_R2constPass->TH1F::Sumw2();
  combPass->TH1F::Sumw2();
  combAll->TH1F::Sumw2();
  bothPass->TH2F::Sumw2();
  bothAll->TH2F::Sumw2();
  doubleR2Pass->TH1F::Sumw2();
  jetPtPass->TH1F::Sumw2();
  jetPtAll->TH1F::Sumw2();
  jetPtconstPass->TH1F::Sumw2();
  jetPtconstAll->TH1F::Sumw2();
  jetPt_MRconstPass->TH1F::Sumw2();
  jetPt_MRconstAll->TH1F::Sumw2();
  jetPt_R2constPass->TH1F::Sumw2();
  jetPt_R2constAll->TH1F::Sumw2();
  jetPt_MRconstPassUneven->TH1F::Sumw2();
  jetPt_MRconstAllUneven->TH1F::Sumw2();
  jetPt_R2constPassUneven->TH1F::Sumw2();
  jetPt_R2constAllUneven->TH1F::Sumw2();
  R2AllUneven->TH1F::Sumw2();
  R2PassUneven->TH1F::Sumw2();

  R2Pass_rsq36_noconst->TH1F::Sumw2();
  R2All_noconst->TH1F::Sumw2();
  R2Pass_rsq36_mr200->TH1F::Sumw2();
  R2All_mr200->TH1F::Sumw2();
  R2Pass_rsq36_mr400->TH1F::Sumw2();
  R2All_mr400->TH1F::Sumw2();
  MRPass_rsq36_rsq36->TH1F::Sumw2();
  MRAll_rsq36->TH1F::Sumw2();
  MRPass_rsq36_rsq40->TH1F::Sumw2();
  MRAll_rsq40->TH1F::Sumw2();
  MRPass_rsq36_rsq50->TH1F::Sumw2();
  MRAll_rsq50->TH1F::Sumw2();

  TriggerTurnonR2Uneven->TH1F::Sumw2();
  TriggerTurnonR2->TH1F::Sumw2();
  TriggerTurnonMR_R2const->TH1F::Sumw2();
  TriggerTurnonR2_MRconst->TH1F::Sumw2();
  TriggerTurnonComb->TH1F::Sumw2();
  TriggerTurnonBoth->TH2F::Sumw2();
  TriggerTurnonDoubleR2->TH1F::Sumw2();
  TriggerTurnonJetPt->TH1F::Sumw2();
  TriggerTurnonJetPtconst->TH1F::Sumw2();
  TriggerTurnonJetPt_MRconst->TH1F::Sumw2();
  TriggerTurnonJetPt_R2const->TH1F::Sumw2();
  TriggerTurnonJetPt_MRconstUneven->TH1F::Sumw2();
  TriggerTurnonJetPt_R2constUneven->TH1F::Sumw2();

  jetPt0Pass_rsq36_rsq40->TH1F::Sumw2();
  jetPt0All_rsq40->TH1F::Sumw2();
  jetPt0Pass_rsq36_rsq50->TH1F::Sumw2();
  jetPt0All_rsq50->TH1F::Sumw2();
  jetPt1Pass_rsq36_rsq40->TH1F::Sumw2();
  jetPt1All_rsq40->TH1F::Sumw2();
  jetPt1Pass_rsq36_rsq50->TH1F::Sumw2();
  jetPt1All_rsq50->TH1F::Sumw2();

  jetPt0Pass_trigger2_3const->TH1F::Sumw2();
  jetPt0All_3const->TH1F::Sumw2();
  jetPt1Pass_trigger2_3const->TH1F::Sumw2();
  jetPt0All_3const->TH1F::Sumw2();

  //output file
  TFile* rootfile1 = TFile::Open("TriggerTurnon.root", "recreate");
  
  //input files
  TFile* f1 = new TFile("Output.root");
 

  //input trees
  TTree* t1 = (TTree*)f1->Get("MultiJet");
 
  //setting necessary tree variables
  float rsq;
  float mr;
  bool hlt[100];
  float combVar;
  float jetptArr[100];
  float jetpt, jetpt1;

  t1->SetBranchAddress("MR", &mr);
  t1->SetBranchAddress("Rsq", &rsq);
  t1->SetBranchAddress("HLTDecision", hlt);
  t1->SetBranchAddress("JetPt", jetptArr);

  //get number of entries
  int n1 = t1->GetEntries();

  //fill histograms
  for (int i = 0; i < n1; i++)
    {
      t1->GetEntry(i);    
      jetpt = jetptArr[0];
      jetpt = jetptArr[1];

      //setting too high values within acceptable region to register in rightmost bin
      if (jetpt >= 1000)
	jetpt = 999;
      if (jetpt1 >= 1000)
	jetpt1 = 999;
      if (rsq >= 1.2)
	rsq = 1;
      if (mr > 5000)
	mr = 4999;

      bothAll->Fill(mr, rsq);
      if (hlt[42] == 1)
	{
	  bothPass->Fill(mr, rsq);
	  TriggerTurnonBoth->Fill(mr, rsq);
	}
      if (mr > 400)
	{

	  R2_MRconstAll->Fill(rsq);
	  if (hlt[38] == 1)
	    {

	      R2_MRconstPass->Fill(rsq);
	      TriggerTurnonR2_MRconst->Fill(rsq);
	    }

	  
	}
      if (rsq > 0.25)
	{
	  MR_R2constAll->Fill(mr);
	  if (hlt[38] == 1)
	    {
	      MR_R2constPass->Fill(mr);
	      TriggerTurnonMR_R2const->Fill(mr);
	    }
	}
      R2All->Fill(rsq);
      R2AllUneven->Fill(rsq);
      if (hlt[42] == 1)
	{
	  R2PassUneven->Fill(rsq);
	  R2Pass->Fill(rsq);
	  TriggerTurnonR2->Fill(rsq);
	  TriggerTurnonR2Uneven->Fill(rsq);
	}
      if ((hlt[40] == 1) || (hlt[42] == 1))
	{
	  doubleR2Pass->Fill(rsq);
	  TriggerTurnonDoubleR2->Fill(rsq);
	}
      if ((mr > 200) && (rsq > 0.25))
	{ 
	  combVar = (mr+300)*(rsq+0.25);
	  combAll->Fill(combVar);
	  if (hlt[38] == 1)
	    {
	      combPass->Fill(combVar);
	      TriggerTurnonComb->Fill(combVar);
	    }
	}
      jetPtAll->Fill(jetpt);
      if (hlt[38] == 1)
	{
	  jetPtPass->Fill(jetpt);
	  TriggerTurnonJetPt->Fill(jetpt);
	}
      if (rsq > 0.25 && mr > 200 && (mr+300)*(rsq+0.25) > 260)
	{
	  jetPtconstAll->Fill(jetpt);
	  if (hlt[38] == 1)
	    {
	      jetPtconstPass->Fill(jetpt);
	      TriggerTurnonJetPtconst->Fill(jetpt);
	    }
	}
      if (rsq > 0.25)
	{
	  jetPt_R2constAllUneven->Fill(jetpt);
	  jetPt_R2constAll->Fill(jetpt);
	  if (hlt[38] == 1)
	    {
	      jetPt_R2constPassUneven->Fill(jetpt);
	      TriggerTurnonJetPt_R2constUneven->Fill(jetpt);
	      jetPt_R2constPass->Fill(jetpt);
	      TriggerTurnonJetPt_R2const->Fill(jetpt);
	    }
	}
      if (mr > 200)
	{
	  jetPt_MRconstAll->Fill(jetpt);
	  jetPt_MRconstAllUneven->Fill(jetpt);
	  if (hlt[38] == 1)
	    {
	      jetPt_MRconstPassUneven->Fill(jetpt);
	      TriggerTurnonJetPt_MRconstUneven->Fill(jetpt);
	      jetPt_MRconstPass->Fill(jetpt);
	      TriggerTurnonJetPt_MRconst->Fill(jetpt);
	    }
	}
      if (mr > 200 && rsq > 0.25 && (mr+300)*(rsq+0.25)>260)
	{
	  jetPt0All->Fill(jetpt);
	  jetPt1All->Fill(jetpt1);
	  MRAll_3const->Fill(mr);
	  R2All_3const->Fill(rsq);
	  if (hlt[42] == 1)
	    {
	      jetPt0Pass_rsq36->Fill(jetpt);
	      jetPt1Pass_rsq36->Fill(jetpt1);
	      MRPass_rsq36_3const->Fill(mr);
	      R2Pass_rsq36_3const->Fill(rsq);
	    }
	  if (hlt[38] == 1)
	    {
	      jetPt0Pass_trigger2->Fill(jetpt);
	      jetPt1Pass_trigger2->Fill(jetpt1);
	      MRPass_trigger2_3const->Fill(mr);
	      R2Pass_trigger2_3const->Fill(rsq);
	    }

	}
      if (mr > 200 && (mr+300)*(rsq+0.25)>260)
     	{
	  R2All_2const->Fill(rsq);
	  if (hlt[42] == 1)
	    R2Pass_rsq36_2const->Fill(rsq);
	  if (hlt[38] == 1)
	    R2Pass_trigger2_2const->Fill(rsq);
	}
      if (rsq > 0.25 && (mr+300)*(rsq+0.25)>260)
	{
	  MRAll_2const->Fill(mr);
	  if (hlt[42] == 1)
	    MRPass_rsq36_2const->Fill(mr);
	  if (hlt[38] == 1)
	    MRPass_trigger2_2const->Fill(mr);
	}
      R2All_noconst->Fill(rsq);
      if (hlt[42] == 1)
	R2Pass_rsq36_noconst->Fill(rsq);
      if (mr > 200)
	{
	  R2All_mr200->Fill(rsq);
	  if (hlt[42] == 1)
	    R2Pass_rsq36_mr200->Fill(rsq);
	  if (mr > 400)
	    {
	      R2All_mr400->Fill(rsq);
	      if (hlt[42] == 1)
		R2Pass_rsq36_mr400->Fill(rsq);
	    }
	}
      if (rsq > 0.36)
	{
	  MRAll_rsq36->Fill(mr);
	  if (hlt[42] == 1)
	    MRPass_rsq36_rsq36->Fill(mr);
	  if (rsq > 0.4)
	    {
	      MRAll_rsq40->Fill(mr);
	      if (hlt[42] == 1)
		MRPass_rsq36_rsq40->Fill(mr);
	      if (rsq > 0.5)
		{
		  MRAll_rsq50->Fill(mr);
		  if (hlt[42] == 1)
		    MRPass_rsq36_rsq50->Fill(mr);
		}
	    }
	}
      if (rsq > 0.09 && (rsq+0.25)*(mr+300)>260)
	{
	  MRAll_rsq09hyper->Fill(mr);
	  if (hlt[38] == 1)
	      MRPass_trigger2_rsq09hyper->Fill(mr);
	  if (rsq > 0.25)
	    {
	      MRAll_rsq25hyper->Fill(mr);
	      if (hlt[38] == 1)
		MRPass_trigger2_rsq25hyper->Fill(mr);
	    }
	}

      if (mr > 200 && (rsq+0.25)*(mr+300)>260)
	{
	  R2All_mr200hyper->Fill(rsq);
	  if (hlt[38] == 1)
	    R2Pass_trigger2_mr200hyper->Fill(rsq);
	  if (mr > 400)
	    {
	      R2All_mr400hyper->Fill(rsq);
	      if (hlt[38] == 1)
		R2Pass_trigger2_mr400hyper->Fill(rsq);
	    }
	}

      if (rsq > 0.40)
	{
	  jetPt0All_rsq40->Fill(jetpt);
	  jetPt1All_rsq40->Fill(jetpt1);
	  if (hlt[42] == 1)
	    {
	      jetPt0Pass_rsq36_rsq40->Fill(jetpt);
	      jetPt1Pass_rsq36_rsq40->Fill(jetpt1);
	    }

	  if (rsq > 0.50)
	    {
	      jetPt0All_rsq50->Fill(jetpt);
	      jetPt1All_rsq50->Fill(jetpt1);
	      if (hlt[42] == 1)
		{
		  jetPt0Pass_rsq36_rsq50->Fill(jetpt);
		  jetPt1Pass_rsq36_rsq50->Fill(jetpt1);
		}

	    }
	}
      if (rsq > 0.25 && mr > 400 && (rsq+0.25)(mr+300)>260)
	{
	  jetPt0All_3const->Fill(jetpt);
	  jetPt1All_3const->Fill(jetpt1);
	  if (hlt[38] == 1)
	    {
	      jetPt0Pass_trigger2_3const->Fill(jetpt);
	      jetPt1Pass_trigger2_3const->Fill(jetpt1);
	    }

	}
    }

  TEfficiency* Eff_JetPt0_rsq36 = new TEfficiency(*jetPt0Pass_rsq36, *jetPt0All);
  TEfficiency* Eff_JetPt1_rsq36 = new TEfficiency(*jetPt1Pass_rsq36, *jetPt1All);
  TEfficiency* Eff_JetPt0_trigger2 = new TEfficiency(*jetPt0Pass_trigger2, *jetPt0All);
  TEfficiency* Eff_JetPt1_trigger2 = new TEfficiency(*jetPt1Pass_trigger2, *jetPt1All);
  TEfficiency* Eff_MR_rsq36_3const = new TEfficiency(*MRPass_rsq36_3const, *MRAll_3const);
  TEfficiency* Eff_MR_rsq36_2const = new TEfficiency(*MRPass_rsq36_2const, *MRAll_2const);
  TEfficiency* Eff_MR_trigger2_3const = new TEfficiency(*MRPass_trigger2_3const, *MRAll_3const);
  TEfficiency* Eff_MR_trigger2_2const = new TEfficiency(*MRPass_trigger2_2const, *MRAll_2const);
  TEfficiency* Eff_R2_rsq36_3const = new TEfficiency(*R2Pass_rsq36_3const, *R2All_3const);
  TEfficiency* Eff_R2_rsq36_2const = new TEfficiency(*R2Pass_rsq36_2const, *R2All_2const);
  TEfficiency* Eff_R2_trigger2_3const = new TEfficiency(*R2Pass_trigger2_3const, *R2All_3const);
  TEfficiency* Eff_R2_trigger2_2const = new TEfficiency(*R2Pass_trigger2_2const, *R2All_2const);

  TEfficiency* Eff_R2_rsq36_noconst = new TEfficiency(*R2Pass_rsq36_noconst, *R2All_noconst);
  TEfficiency* Eff_R2_rsq36_mr200 = new TEfficiency(*R2Pass_rsq36_mr200, *R2All_mr200);
  TEfficiency* Eff_R2_rsq36_mr400 = new TEfficiency(*R2Pass_rsq36_mr400, *R2All_mr400);
  TEfficiency* Eff_MR_rsq36_rsq36 = new TEfficiency(*MRPass_rsq36_rsq36, *MRAll_rsq36);
  TEfficiency* Eff_MR_rsq36_rsq40 = new TEfficiency(*MRPass_rsq36_rsq40, *MRAll_rsq40);
  TEfficiency* Eff_MR_rsq36_rsq50 = new TEfficiency(*MRPass_rsq36_rsq50, *MRAll_rsq50);

  TEfficiency* Eff_R2_trigger2_mr200hyper = new TEfficiency(*R2Pass_trigger2_mr200hyper, *R2All_mr200hyper);
  TEfficiency* Eff_R2_trigger2_mr400hyper = new TEfficiency(*R2Pass_trigger2_mr400hyper, *R2All_mr400hyper);
  TEfficiency* Eff_MR_trigger2_rsq09hyper = new TEfficiency(*MRPass_trigger2_rsq09hyper, *MRAll_rsq09hyper);
  TEfficiency* Eff_MR_trigger2_rsq25hyper = new TEfficiency(*MRPass_trigger2_rsq25hyper, *MRAll_rsq25hyper);

  TEfficiency* Eff_jetPt0_rsq36_rsq40 = new TEfficiency(*jetPt0Pass_rsq36_rsq40, *jetPt0All_rsq40);
  TEfficiency* Eff_jetPt0_rsq36_rsq50 = new TEfficiency(*jetPt0Pass_rsq36_rsq50, *jetPt0All_rsq50);
  TEfficiency* Eff_jetPt1_rsq36_rsq40 = new TEfficiency(*jetPt1Pass_rsq36_rsq40, *jetPt1All_rsq40);
  TEfficiency* Eff_jetPt1_rsq36_rsq50 = new TEfficiency(*jetPt1Pass_rsq36_rsq50, *jetPt1All_rsq50);

  TEfficiency* Eff_jetPt0_trigger2_3const = new TEfficiency(*jetPt0Pass_trigger2_3const, *jetPt0All_3const);
  TEfficiency* Eff_jetPt1_trigger2_3const = new TEfficiency(*jetPt1Pass_trigger2_3const, *jetPt1All_3const);

  TriggerTurnonComb->Divide(combAll);
  TriggerTurnonR2->Divide(R2All);
  TriggerTurnonMR_R2const->Divide(MR_R2constAll);
  TriggerTurnonR2_MRconst->Divide(R2_MRconstAll);
  TriggerTurnonBoth->Divide(bothAll);
  TriggerTurnonDoubleR2->Divide(R2All);
  TriggerTurnonJetPt->Divide(jetPtAll);
  TriggerTurnonJetPtconst->Divide(jetPtconstAll);
  TriggerTurnonJetPt_MRconst->Divide(jetPt_MRconstAll);
  TriggerTurnonJetPt_R2const->Divide(jetPt_R2constAll);
  TriggerTurnonJetPt_MRconstUneven->Divide(jetPt_MRconstAllUneven);
  TriggerTurnonJetPt_R2constUneven->Divide(jetPt_R2constAllUneven);
  TriggerTurnonR2Uneven->Divide(R2AllUneven);

  //move to output file
  rootfile1->cd();
  
  //write histograms to output file
  bothPass->Write();
  bothAll->Write();
  R2Pass->Write();
  R2All->Write();
  R2_MRconstPass->Write();
  R2_MRconstAll->Write();
  MR_R2constPass->Write();
  MR_R2constAll->Write();
  combPass->Write();
  combAll->Write();
  doubleR2Pass->Write();
  jetPtPass->Write();
  jetPtAll->Write();
  jetPtconstPass->Write();
  jetPtconstAll->Write();
  jetPt_MRconstPass->Write();
  jetPt_MRconstAll->Write();
  jetPt_R2constPass->Write();
  jetPt_R2constAll->Write();
  jetPt_MRconstPassUneven->Write();
  jetPt_MRconstAllUneven->Write();
  jetPt_R2constPassUneven->Write();
  jetPt_R2constAllUneven->Write();
  R2PassUneven->Write();
  R2AllUneven->Write();

  TriggerTurnonBoth->Write();
  TriggerTurnonMR_R2const->Write();
  TriggerTurnonR2_MRconst->Write();
  TriggerTurnonR2->Write();
  TriggerTurnonComb->Write();
  TriggerTurnonDoubleR2->Write();
  TriggerTurnonJetPt->Write();
  TriggerTurnonJetPtconst->Write();
  TriggerTurnonJetPt_MRconst->Write();
  TriggerTurnonJetPt_R2const->Write();
  TriggerTurnonJetPt_MRconstUneven->Write();
  TriggerTurnonJetPt_R2constUneven->Write();
  TriggerTurnonR2Uneven->Write();

  Eff_JetPt0_rsq36->Write("Eff_JetPt0_rsq36");
  Eff_JetPt1_rsq36->Write("Eff_JetPt1_rsq36");
  Eff_JetPt0_trigger2->Write("Eff_JetPt0_trigger2");
  Eff_JetPt1_trigger2->Write("Eff_JetPt1_trigger2");
  Eff_MR_rsq36_3const->Write("Eff_MR_rsq36_3const");
  Eff_MR_rsq36_2const->Write("Eff_MR_rsq36_2const");
  Eff_MR_trigger2_3const->Write("Eff_MR_trigger2_3const");
  Eff_MR_trigger2_2const->Write("Eff_MR_trigger2_2const");
  Eff_R2_rsq36_3const->Write("Eff_R2_rsq36_3const");
  Eff_R2_rsq36_2const->Write("Eff_R2_rsq36_2const");
  Eff_R2_trigger2_3const->Write("Eff_R2_trigger2_3const");
  Eff_R2_trigger2_2const->Write("Eff_R2_trigger2_2const");

  Eff_R2_rsq36_noconst->Write("Eff_R2_rsq36_noconst");
  Eff_R2_rsq36_mr200->Write("Eff_R2_rsq36_mr200");
  Eff_R2_rsq36_mr400->Write("Eff_R2_rsq36_mr400");
  Eff_MR_rsq36_rsq36->Write("Eff_MR_rsq36_rsq36");
  Eff_MR_rsq36_rsq40->Write("Eff_MR_rsq36_rsq40");
  Eff_MR_rsq36_rsq50->Write("Eff_MR_rsq36_rsq50");

  Eff_R2_trigger2_mr200hyper->Write("Eff_R2_trigger2_mr200hyper");
  Eff_R2_trigger2_mr400hyper->Write("Eff_R2_trigger2_mr400hyper");
  Eff_MR_trigger2_rsq09hyper->Write("Eff_MR_trigger2_rsq09hyper");
  Eff_MR_trigger2_rsq25hyper->Write("Eff_MR_trigger2_rsq25hyper");

  Eff_jetPt0_rsq36_rsq40->Write("Eff_jetPt0_rsq36_rsq40");
  Eff_jetPt0_rsq36_rsq50->Write("Eff_jetPt0_rsq36_rsq50");
  Eff_jetPt1_rsq36_rsq40->Write("Eff_jetPt1_rsq36_rsq40");
  Eff_jetPt1_rsq36_rsq50->Write("Eff_jetPt1_rsq36_rsq50");

  Eff_jetPt0_trigger2_3const->Write("Eff_jetPt0_trigger2_3const");
  Eff_jetPt1_trigger2_3const->Write("Eff_jetPt1_trigger2_3const");

  rootfile1->Close();


}
