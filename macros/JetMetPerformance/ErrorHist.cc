#include "TH1F.h"

using namespace std;

void ErrorHist()
{
  //Define boundaries on intervalRatios histogram: nBins = size of binBounds - 1
  float ptBinBounds[] = {30,40,50,60,70,80,90,100,200,500,1000};
  int ptNBins = 10;

  float etaBinBounds[] = {-6, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 6};
  int etaNBins = 12;

  float phiBinBounds[] = {-3.15, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3.15};
  int phiNBins = 12;

  //Getting the values from the tree in RazorDM.root
  const char* filenames[] = {"chi150phi1000.root", "chi150phi500.root", "chi50phi150.root", "chi50phi500.root"};
  int N = 4;

  for (int entry = 0; entry < N; entry++)
    {
      std::string namestr(filenames[entry]);
      string outstr = "Hists_" + namestr;
      outName = outstr.c_str();
      cout << filenames[entry] << ":" << outName << endl;
      TFile* rootfile1 = TFile::Open(outName, "recreate");
      TFile *f = new TFile(filenames[entry]);
      TTree *t = (TTree*)f->Get("MultiJet");
      
      //Define histograms

      TH1F* ptErrRatio30_40 = new TH1F("ptErrRatio30_40", "(Pt_reco - Pt_gen)/Pt_gen with 30 < Pt_gen < 40", 100, -3, 3);
      TH1F* ptErrRatio40_50 = new TH1F("ptErrRatio40_50", "(Pt_reco - Pt_gen)/Pt_gen with 40 < Pt_gen < 50", 100, -3, 3);
      TH1F* ptErrRatio50_60 = new TH1F("ptErrRatio50_60", "(Pt_reco - Pt_gen)/Pt_gen with 50 < Pt_gen < 60", 100, -3, 3);
      TH1F* ptErrRatio60_70 = new TH1F("ptErrRatio60_70", "(Pt_reco - Pt_gen)/Pt_gen with 60 < Pt_gen < 70", 100, -3, 3);
      TH1F* ptErrRatio70_80 = new TH1F("ptErrRatio70_80", "(Pt_reco - Pt_gen)/Pt_gen with 70 < Pt_gen < 80", 100, -3, 3);
      TH1F* ptErrRatio80_90 = new TH1F("ptErrRatio80_90", "(Pt_reco - Pt_gen)/Pt_gen with 80 < Pt_gen < 90", 100, -3, 3);
      TH1F* ptErrRatio90_100 = new TH1F("ptErrRatio90_100", "(Pt_reco - Pt_gen)/Pt_gen with 90 < Pt_gen < 100", 100, -3, 3);
      TH1F* ptErrRatio100_200 = new TH1F("ptErrRatio100_200", "(Pt_reco - Pt_gen)/Pt_gen with 100 < Pt_gen < 200", 100, -3, 3);
      TH1F* ptErrRatio200_500 = new TH1F("ptErrRatio200_500", "(Pt_reco - Pt_gen)/Pt_gen with 200 < Pt_gen < 500", 100, -3, 3);
      TH1F* ptErrRatio500_1000 = new TH1F("ptErrRatio500_1000", "(Pt_reco - Pt_gen)/Pt_gen with 500 < Pt_gen < 1000", 100, -3, 3);

      TH1F* ptErrIntervalRatios = new TH1F("ptErrIntervalRatios", "Mean value of (Pt_reco - Pt_gen)/Pt_gen for each pt interval", ptNBins, ptBinBounds);

   
      short fTitleColor;      

      TH1F* phiRatio1 = new TH1F("phiRatio_nPi_n5_2", "Pt_reco/Pt_gen with -Pi < phi_gen < -2.5", 100, 0, 5);
      TH1F* phiRatio2 = new TH1F("phiRatio_n5_2_n2", "Pt_reco/Pt_gen with -2.5 < phi_gen < -2", 100, 0, 5);
      TH1F* phiRatio3 = new TH1F("phiRatio_n2_n3_2", "Pt_reco/Pt_gen with -2 < phi_gen < -1.5", 100, 0, 5);
      TH1F* phiRatio4 = new TH1F("phiRatio_n3_2_n1", "Pt_reco/Pt_gen with -1.5 < phi_gen < -1", 100, 0, 5);
      TH1F* phiRatio5 = new TH1F("phiRatio_n1_n1_2", "Pt_reco/Pt_gen with -1 < phi_gen < -0.5", 100, 0, 5);
      TH1F* phiRatio6 = new TH1F("phiRatio_n1_2_0", "Pt_reco/Pt_gen with -0.5 < phi_gen < 0", 100, 0, 5);
      TH1F* phiRatio7 = new TH1F("phiRatio_0_1_2", "Pt_reco/Pt_gen with 0 < phi_gen < 0.5", 100, 0, 5);
      TH1F* phiRatio8 = new TH1F("phiRatio_1_2_1", "Pt_reco/Pt_gen with 0.5 < phi_gen < 1", 100, 0, 5);
      TH1F* phiRatio9 = new TH1F("phiRatio_1_3_2", "Pt_reco/Pt_gen with 1 < phi_gen < 1.5", 100, 0, 5);
      TH1F* phiRatio10 = new TH1F("phiRatio_3_2_2", "Pt_reco/Pt_gen with 1.5 < phi_gen < 2", 100, 0, 5);
      TH1F* phiRatio11 = new TH1F("phiRatio_2_5_2", "Pt_reco/Pt_gen with 2 < phi_gen < 2.5", 100, 0, 5);
      TH1F* phiRatio12 = new TH1F("phiRatio_5_2_Pi", "Pt_reco/Pt_gen with 2.5 < phi_gen < Pi", 100, 0, 5);

      TH1F* phiIntervalRatios = new TH1F("phiIntervalRatios", "Mean ratio of Pt_reco/Pt_gen for each phi interval", phiNBins, phiBinBounds);


      TH1F* etaRatio1 = new TH1F("etaRatio_n6_n3", "Pt_reco/Pt_gen with -6 < eta_gen < -3", 100, 0, 5);
      TH1F* etaRatio2 = new TH1F("etaRatio_n3_n2", "Pt_reco/Pt_gen with -3 < eta_gen < -2", 100, 0, 5);
      TH1F* etaRatio3 = new TH1F("etaRatio_n2_n3_2", "Pt_reco/Pt_gen with -2 < eta_gen < -1.5", 100, 0, 5);
      TH1F* etaRatio4 = new TH1F("etaRatio_n3_2_n1", "Pt_reco/Pt_gen with -1.5 < eta_gen < -1", 100, 0, 5);
      TH1F* etaRatio5 = new TH1F("etaRatio_n1_n1_2", "Pt_reco/Pt_gen with -1 < eta_gen < -0.5", 100, 0, 5);
      TH1F* etaRatio6 = new TH1F("etaRatio_n1_2_0", "Pt_reco/Pt_gen with -0.5 < eta_gen < -0", 100, 0, 5);
      TH1F* etaRatio7 = new TH1F("etaRatio_0_1_2", "Pt_reco/Pt_gen with 0 < eta_gen < 0.5", 100, 0, 5);
      TH1F* etaRatio8 = new TH1F("etaRatio_1_2_1", "Pt_reco/Pt_gen with 0.5 < eta_gen < 1", 100, 0, 5);
      TH1F* etaRatio9 = new TH1F("etaRatio_1_3_2", "Pt_reco/Pt_gen with 1 < eta_gen < 1.5", 100, 0, 5);
      TH1F* etaRatio10 = new TH1F("etaRatio_3_2_2", "Pt_reco/Pt_gen with 1.5 < eta_gen < 2", 100, 0, 5);
      TH1F* etaRatio11 = new TH1F("etaRatio_2_3", "Pt_reco/Pt_gen with 2 < eta_gen < 3", 100, 0, 5);
      TH1F* etaRatio12 = new TH1F("etaRatio_3_6", "Pt_reco/Pt_gen with 3 < eta_gen < 6", 100, 0, 5);
     
      TH1F* etaIntervalRatios = new TH1F("etaIntervalRatios", "Mean ratio of Pt_reco/Pt_gen for each eta interval", etaNBins, etaBinBounds);

   
      TH1F* ptRatio30_40 = new TH1F("ptRatio30_40", "Pt_reco/Pt_gen with 30 < Pt_gen < 40", 100, 0, 5);
      TH1F* ptRatio40_50 = new TH1F("ptRatio40_50", "Pt_reco/Pt_gen with 40 < Pt_gen < 50", 100, 0, 5);
      TH1F* ptRatio50_60 = new TH1F("ptRatio50_60", "Pt_reco/Pt_gen with 50 < Pt_gen < 60", 100, 0, 5);
      TH1F* ptRatio60_70 = new TH1F("ptRatio60_70", "Pt_reco/Pt_gen with 60 < Pt_gen < 70", 100, 0, 5);
      TH1F* ptRatio70_80 = new TH1F("ptRatio70_80", "Pt_reco/Pt_gen with 70 < Pt_gen < 80", 100, 0, 5);
      TH1F* ptRatio80_90 = new TH1F("ptRatio80_90", "Pt_reco/Pt_gen with 80 < Pt_gen < 90", 100, 0, 5);
      TH1F* ptRatio90_100 = new TH1F("ptRatio90_100", "Pt_reco/Pt_gen with 90 < Pt_gen < 100", 100, 0, 5);
      TH1F* ptRatio100_200 = new TH1F("ptRatio100_200", "Pt_reco/Pt_gen with 100 < Pt_gen < 200", 100, 0, 5);
      TH1F* ptRatio200_500 = new TH1F("ptRatio200_500", "Pt_reco/Pt_gen with 200 < Pt_gen < 500", 100, 0, 5);
      TH1F* ptRatio500_1000 = new TH1F("ptRatio500_1000", "Pt_reco/Pt_gen with 500 < Pt_gen < 1000", 100, 0, 5);
      TH1F* ptRatio1000_up = new TH1F("ptRatio1000_up", "Pt_reco/Pt_gen with 1000 < Pt_gen", 100, 0, 5);

      TH1F* ptIntervalRatios = new TH1F("ptIntervalRatios", "Mean ratio of Pt_reco/Pt_gen for each pt interval", ptNBins, ptBinBounds);
      
      TH1F* phiErr = new TH1F("phiErr", "phi comparison", 100, -0.3, 0.3);
      TH1F* etaErr = new TH1F("etaErr", "eta comparison", 100, -0.5, 0.5);
      TH1F* EErr = new TH1F("EErr", "E comparison", 100, -100, 100);
      TH1F* ptErr = new TH1F("ptErr", "pt comparison", 100, -100, 100);
      TH1F* ptErr_uncorr = new TH1F("ptErr_uncorr", "uncorrected pt comparison", 100, -100, 100);
      TH1F* metPtErr = new TH1F("metPtErr", "MET Pt comparison", 100, -200, 600);
      TH1F* metPhiErr = new TH1F("metPhiErr", "MET Phi comparison", 100, -10, 10);
      TH1F* metT0PtErr = new TH1F("metT0PtErr", "MET Type 0 Pt comparison", 100, -200,600);
      TH1F* metT0PhiErr = new TH1F("metT0PhiErr", "MET Type 0 Phi comparison", 100, -10,10);
      TH1F* metT1PtErr = new TH1F("metT1PtErr", "MET Type 1 Pt comparison", 100, -200,600);
      TH1F* metT1PhiErr = new TH1F("metT1PhiErr", "MET Type 1 Phi comparison", 100, -10,10);
      TH1F* metT0T1PtErr = new TH1F("metT0T1PtErr", "MET Type 0+1 Pt comparison", 100, -200,600);
      TH1F* metT0T1PhiErr = new TH1F("metT0T1PhiErr", "MET Type 0+1 Phi comparison", 100, -10,10);
      TH1F* metPt = new TH1F("metPt", "MET Pt values", 100,0, 1000);
      TH1F* metT0Pt = new TH1F("metT0Pt", "MET Type 0 Pt values",100, -400, 400);
      TH1F* metT1Pt = new TH1F("metT1Pt", "MET Type 1 Pt values",100, -400, 400);
      TH1F* metT0T1Pt = new TH1F("metT0T1Pt", "MET Type 0+1 Pt values",100, -400, 400);      
      

      //Assign variables to necessary tree branches of MultiJet tree
      float uncorrPt[30], corrE[30], corrPt[30], eta[30], phi[30], genE[100], genPt[100], genPhi[100], genEta[100];
      bool hasMatch[30];
      int selected, matchIndex[30];
      t->SetBranchAddress("genJetPhi", genPhi);
      t->SetBranchAddress("genJetPt", genPt);
      t->SetBranchAddress("genJetEta", genEta);
      t->SetBranchAddress("genJetE", genE);
      t->SetBranchAddress("JetPhi", phi);
      t->SetBranchAddress("JetEta", eta);
      t->SetBranchAddress("JetE", corrE);
      t->SetBranchAddress("JetPt", corrPt);
      t->SetBranchAddress("JetPt_uncorr", uncorrPt);
      t->SetBranchAddress("nSelectedJets", &selected);
      t->SetBranchAddress("hasMatchingGenJet", hasMatch);
      t->SetBranchAddress("matchingGenJetIndex", matchIndex);

      float mPt, mPhi, mT0Pt, mT0Phi, mT1Pt, mT1Phi, mT0T1Pt, mT0T1Phi, gmPt, gmPhi;
      t->SetBranchAddress("genMetPt", &gmPt);
      t->SetBranchAddress("genMetPhi", &gmPhi);
      t->SetBranchAddress("metPt", &mPt);
      t->SetBranchAddress("metPhi", &mPhi);
      t->SetBranchAddress("metType0Pt", &mT0Pt);
      t->SetBranchAddress("metType0Phi", &mT0Phi);
      t->SetBranchAddress("metType1Pt", &mT1Pt);
      t->SetBranchAddress("metType1Phi", &mT1Phi);
      t->SetBranchAddress("metType0Plus1Pt", &mT0T1Pt);
      t->SetBranchAddress("metType0Plus1Phi", &mT0T1Phi);


      int nentries = t->GetEntries();
 
      //Writing values into histogram
      for (int entryIndex = 0; entryIndex < nentries; entryIndex++) {
	t->GetEntry(entryIndex);
	
	//Fill MET histograms
	metPtErr->Fill(mPt - gmPt);
	metPhiErr->Fill(mPhi - gmPhi);
	metT0PtErr->Fill(mT0Pt - gmPt);
	metT0PhiErr->Fill(mT0Phi - gmPhi);
	metT1PtErr->Fill(mT1Pt - gmPt);
	metT1PhiErr->Fill(mT1Phi - gmPhi);
	metT0T1PtErr->Fill(mT0T1Pt - gmPt);
	metT0T1PhiErr->Fill(mT0T1Phi - gmPhi);
	metPt->Fill(mPt);
	metT0Pt->Fill(mT0Pt);
	metT1Pt->Fill(mT1Pt);
	metT0T1Pt->Fill(mT0T1Pt);
	//Fill reconstructed jet histograms and ratio histograms
	for (int i = 0; i < selected; i++){
	  if (hasMatch[i]){
	    
	    phiErr->Fill(phi[i] - genPhi[matchIndex[i]]);
	    etaErr->Fill(eta[i] - genEta[matchIndex[i]]);
	    EErr->Fill(corrE[i] - genE[matchIndex[i]]);
	    ptErr->Fill(corrPt[i] - genPt[matchIndex[i]]);
	    ptErr_uncorr->Fill(uncorrPt[i] - genPt[matchIndex[i]]);
	   
	    //Filling up the ptRatio plots
	    if (genPt[matchIndex[i]] > 1000){
	      ptRatio1000_up->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }
	    
	    if (genPt[matchIndex[i]] > 500 && genPt[matchIndex[i]] < 1000){
	      ptRatio500_1000->Fill(corrPt[i]/genPt[matchIndex[i]]);
	      ptErrRatio500_1000->Fill((corrPt[i] - genPt[matchIndex[i]])/genPt[matchIndex[i]]);
	    }
	   
	    if (genPt[matchIndex[i]] > 200 && genPt[matchIndex[i]] < 500){
	      ptRatio200_500->Fill(corrPt[i]/genPt[matchIndex[i]]);
	      ptErrRatio200_500->Fill((corrPt[i] - genPt[matchIndex[i]])/genPt[matchIndex[i]]);
	    }

	    if (genPt[matchIndex[i]] > 100 && genPt[matchIndex[i]] < 200){
	      ptRatio100_200->Fill(corrPt[i]/genPt[matchIndex[i]]);
	      ptErrRatio100_200->Fill((corrPt[i] - genPt[matchIndex[i]])/genPt[matchIndex[i]]);
	    }

	    if (genPt[matchIndex[i]] > 90 && genPt[matchIndex[i]] < 100){
	      ptRatio90_100->Fill(corrPt[i]/genPt[matchIndex[i]]);
	      ptErrRatio90_100->Fill((corrPt[i] - genPt[matchIndex[i]])/genPt[matchIndex[i]]);
	    }

	    if (genPt[matchIndex[i]] > 80 && genPt[matchIndex[i]] < 90){
	      ptRatio80_90->Fill(corrPt[i]/genPt[matchIndex[i]]);
	      ptErrRatio80_90->Fill((corrPt[i] - genPt[matchIndex[i]])/genPt[matchIndex[i]]);
	    }

	    if (genPt[matchIndex[i]] > 70 && genPt[matchIndex[i]] < 80){
	      ptRatio70_80->Fill(corrPt[i]/genPt[matchIndex[i]]);
	      ptErrRatio70_80->Fill((corrPt[i] - genPt[matchIndex[i]])/genPt[matchIndex[i]]);
	    }

	    if (genPt[matchIndex[i]] > 60 && genPt[matchIndex[i]] < 70){
	      ptRatio60_70->Fill(corrPt[i]/genPt[matchIndex[i]]);
	      ptErrRatio60_70->Fill((corrPt[i] - genPt[matchIndex[i]])/genPt[matchIndex[i]]);
	    }

	    if (genPt[matchIndex[i]] > 50 && genPt[matchIndex[i]] < 60){
	      ptRatio50_60->Fill(corrPt[i]/genPt[matchIndex[i]]);
	      ptErrRatio50_60->Fill((corrPt[i] - genPt[matchIndex[i]])/genPt[matchIndex[i]]);
	    }

	    if (genPt[matchIndex[i]] > 40 && genPt[matchIndex[i]] < 50){
	      ptRatio40_50->Fill(corrPt[i]/genPt[matchIndex[i]]);
	      ptErrRatio40_50->Fill((corrPt[i] - genPt[matchIndex[i]])/genPt[matchIndex[i]]);
	    }

	    if (genPt[matchIndex[i]] > 30 && genPt[matchIndex[i]] < 40){
	      ptRatio30_40->Fill(corrPt[i]/genPt[matchIndex[i]]);
	      ptErrRatio30_40->Fill((corrPt[i] - genPt[matchIndex[i]])/genPt[matchIndex[i]]);
	    }

	    //Filling up the etaRatio plots


	    if (genEta[matchIndex[i]] > -6 && genEta[matchIndex[i]] < -3){
	      etaRatio1->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genEta[matchIndex[i]] > -3 && genEta[matchIndex[i]] < -2){
	      etaRatio2->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genEta[matchIndex[i]] > -2 && genEta[matchIndex[i]] < -1.5){
	      etaRatio3->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genEta[matchIndex[i]] > -1.5 && genEta[matchIndex[i]] < -1){
	      etaRatio4->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genEta[matchIndex[i]] > -1 && genEta[matchIndex[i]] < -0.5){
	      etaRatio5->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genEta[matchIndex[i]] > -0.5 && genEta[matchIndex[i]] < 0){
	      etaRatio6->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genEta[matchIndex[i]] > 0 && genEta[matchIndex[i]] < 0.5){
	      etaRatio7->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genEta[matchIndex[i]] > 0.5 && genEta[matchIndex[i]] < 1){
	      etaRatio8->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genEta[matchIndex[i]] > 1 && genEta[matchIndex[i]] < 1.5){
	      etaRatio9->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genEta[matchIndex[i]] > 1.5 && genEta[matchIndex[i]] < 2){
	      etaRatio10->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genEta[matchIndex[i]] > 2 && genEta[matchIndex[i]] < 6){
	      etaRatio11->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genEta[matchIndex[i]] > 3 && genEta[matchIndex[i]] < 6){
	      etaRatio12->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    //Filling up the phiRatio plots

	    if (genPhi[matchIndex[i]] > -3.15 && genEta[matchIndex[i]] < -2.5){
	      phiRatio1->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genPhi[matchIndex[i]] > -2.5 && genEta[matchIndex[i]] < -2){
	      phiRatio2->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genPhi[matchIndex[i]] > -2 && genEta[matchIndex[i]] < -1.5){
	      phiRatio3->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genPhi[matchIndex[i]] > -1.5 && genEta[matchIndex[i]] < -1){
	      phiRatio4->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genPhi[matchIndex[i]] > -1 && genEta[matchIndex[i]] < -0.5){
	      phiRatio5->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genPhi[matchIndex[i]] > -0.5 && genEta[matchIndex[i]] < 0){
	      phiRatio6->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genPhi[matchIndex[i]] > 0 && genEta[matchIndex[i]] < 0.5){
	      phiRatio7->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genPhi[matchIndex[i]] > 0.5 && genEta[matchIndex[i]] < 1){
	      phiRatio8->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genPhi[matchIndex[i]] > 1 && genEta[matchIndex[i]] < 1.5){
	      phiRatio9->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genPhi[matchIndex[i]] > 1.5 && genEta[matchIndex[i]] < 2){
	      phiRatio10->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genPhi[matchIndex[i]] > 2 && genEta[matchIndex[i]] < 2.5){
	      phiRatio11->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }

	    if (genPhi[matchIndex[i]] > 2.5 && genEta[matchIndex[i]] < 3.15){
	      phiRatio12->Fill(corrPt[i]/genPt[matchIndex[i]]);
	    }
	  }
	  
	}
	      
      }

      //Filling up ptIntervalRatios
      ptIntervalRatios->SetBinContent(1, ptRatio30_40->GetMean());
      ptIntervalRatios->SetBinError(1, ptRatio30_40->GetRMS());
      ptIntervalRatios->SetBinContent(2, ptRatio40_50->GetMean());
      ptIntervalRatios->SetBinError(2, ptRatio40_50->GetRMS());
      ptIntervalRatios->SetBinContent(3, ptRatio50_60->GetMean());
      ptIntervalRatios->SetBinError(3, ptRatio50_60->GetRMS());
      ptIntervalRatios->SetBinContent(4, ptRatio60_70->GetMean());
      ptIntervalRatios->SetBinError(4, ptRatio60_70->GetRMS());
      ptIntervalRatios->SetBinContent(5, ptRatio70_80->GetMean());
      ptIntervalRatios->SetBinError(5, ptRatio70_80->GetRMS());
      ptIntervalRatios->SetBinContent(6, ptRatio80_90->GetMean());
      ptIntervalRatios->SetBinError(6, ptRatio80_90->GetRMS());
      ptIntervalRatios->SetBinContent(7, ptRatio90_100->GetMean());
      ptIntervalRatios->SetBinError(7, ptRatio90_100->GetRMS());
      ptIntervalRatios->SetBinContent(8, ptRatio100_200->GetMean());
      ptIntervalRatios->SetBinError(8, ptRatio100_200->GetRMS());
      ptIntervalRatios->SetBinContent(9, ptRatio200_500->GetMean());
      ptIntervalRatios->SetBinError(9, ptRatio200_500->GetRMS());
      ptIntervalRatios->SetBinContent(10, ptRatio500_1000->GetMean());
      ptIntervalRatios->SetBinError(10, ptRatio500_1000->GetRMS());

      //Filling up ptErrIntervalRatios
      ptErrIntervalRatios->SetBinContent(1, ptErrRatio30_40->GetMean());
      ptErrIntervalRatios->SetBinError(1, ptErrRatio30_40->GetRMS());
      ptErrIntervalRatios->SetBinContent(2, ptErrRatio40_50->GetMean());
      ptErrIntervalRatios->SetBinError(2, ptErrRatio40_50->GetRMS());
      ptErrIntervalRatios->SetBinContent(3, ptErrRatio50_60->GetMean());
      ptErrIntervalRatios->SetBinError(3, ptErrRatio50_60->GetRMS());
      ptErrIntervalRatios->SetBinContent(4, ptErrRatio60_70->GetMean());
      ptErrIntervalRatios->SetBinError(4, ptErrRatio60_70->GetRMS());
      ptErrIntervalRatios->SetBinContent(5, ptErrRatio70_80->GetMean());
      ptErrIntervalRatios->SetBinError(5, ptErrRatio70_80->GetRMS());
      ptErrIntervalRatios->SetBinContent(6, ptErrRatio80_90->GetMean());
      ptErrIntervalRatios->SetBinError(6, ptErrRatio80_90->GetRMS());
      ptErrIntervalRatios->SetBinContent(7, ptErrRatio90_100->GetMean());
      ptErrIntervalRatios->SetBinError(7, ptErrRatio90_100->GetRMS());
      ptErrIntervalRatios->SetBinContent(8, ptErrRatio100_200->GetMean());
      ptErrIntervalRatios->SetBinError(8, ptErrRatio100_200->GetRMS());
      ptErrIntervalRatios->SetBinContent(9, ptErrRatio200_500->GetMean());
      ptErrIntervalRatios->SetBinError(9, ptErrRatio200_500->GetRMS());
      ptErrIntervalRatios->SetBinContent(10, ptErrRatio500_1000->GetMean());
      ptErrIntervalRatios->SetBinError(10, ptErrRatio500_1000->GetRMS());

      //Filling up etaIntervalRatios
      etaIntervalRatios->SetBinContent(1, etaRatio1->GetMean());
      etaIntervalRatios->SetBinContent(2, etaRatio2->GetMean());
      etaIntervalRatios->SetBinContent(3, etaRatio3->GetMean());
      etaIntervalRatios->SetBinContent(4, etaRatio4->GetMean());
      etaIntervalRatios->SetBinContent(5, etaRatio5->GetMean());
      etaIntervalRatios->SetBinContent(6, etaRatio6->GetMean());
      etaIntervalRatios->SetBinContent(7, etaRatio7->GetMean());
      etaIntervalRatios->SetBinContent(8, etaRatio8->GetMean());
      etaIntervalRatios->SetBinContent(9, etaRatio9->GetMean());
      etaIntervalRatios->SetBinContent(10, etaRatio10->GetMean());
      etaIntervalRatios->SetBinContent(11, etaRatio11->GetMean());
      etaIntervalRatios->SetBinContent(12, etaRatio12->GetMean());
      etaIntervalRatios->SetBinError(1, etaRatio1->GetRMS());
      etaIntervalRatios->SetBinError(2, etaRatio2->GetRMS());
      etaIntervalRatios->SetBinError(3, etaRatio3->GetRMS());
      etaIntervalRatios->SetBinError(4, etaRatio4->GetRMS());
      etaIntervalRatios->SetBinError(5, etaRatio5->GetRMS());
      etaIntervalRatios->SetBinError(6, etaRatio6->GetRMS());
      etaIntervalRatios->SetBinError(7, etaRatio7->GetRMS());
      etaIntervalRatios->SetBinError(8, etaRatio8->GetRMS());
      etaIntervalRatios->SetBinError(9, etaRatio9->GetRMS());
      etaIntervalRatios->SetBinError(10, etaRatio10->GetRMS());
      etaIntervalRatios->SetBinError(11, etaRatio11->GetRMS());
      etaIntervalRatios->SetBinError(12, etaRatio12->GetRMS());


      //Filling up phiIntervalRatios
      phiIntervalRatios->SetBinContent(1, phiRatio1->GetMean());
      phiIntervalRatios->SetBinContent(2, phiRatio2->GetMean());
      phiIntervalRatios->SetBinContent(3, phiRatio3->GetMean());
      phiIntervalRatios->SetBinContent(4, phiRatio4->GetMean());
      phiIntervalRatios->SetBinContent(5, phiRatio5->GetMean());
      phiIntervalRatios->SetBinContent(6, phiRatio6->GetMean());
      phiIntervalRatios->SetBinContent(7, phiRatio7->GetMean());
      phiIntervalRatios->SetBinContent(8, phiRatio8->GetMean());
      phiIntervalRatios->SetBinContent(9, phiRatio9->GetMean());
      phiIntervalRatios->SetBinContent(10, phiRatio10->GetMean());
      phiIntervalRatios->SetBinContent(11, phiRatio11->GetMean());
      phiIntervalRatios->SetBinContent(12, phiRatio12->GetMean());
      phiIntervalRatios->SetBinError(1, phiRatio1->GetRMS());
      phiIntervalRatios->SetBinError(2, phiRatio2->GetRMS());
      phiIntervalRatios->SetBinError(3, phiRatio3->GetRMS());
      phiIntervalRatios->SetBinError(4, phiRatio4->GetRMS());
      phiIntervalRatios->SetBinError(5, phiRatio5->GetRMS());
      phiIntervalRatios->SetBinError(6, phiRatio6->GetRMS());
      phiIntervalRatios->SetBinError(7, phiRatio7->GetRMS());
      phiIntervalRatios->SetBinError(8, phiRatio8->GetRMS());
      phiIntervalRatios->SetBinError(9, phiRatio9->GetRMS());
      phiIntervalRatios->SetBinError(10, phiRatio10->GetRMS());
      phiIntervalRatios->SetBinError(11, phiRatio11->GetRMS());
      phiIntervalRatios->SetBinError(12, phiRatio12->GetRMS());


      //Writing output
      rootfile1->cd();

      //Main interval histograms
      ptIntervalRatios->Write();
      etaIntervalRatios->Write();
      phiIntervalRatios->Write();
      ptErrIntervalRatios->Write();
      //All the individual ratio plots
      ptErrRatio30_40->Write();
      ptErrRatio40_50->Write();
      ptErrRatio50_60->Write();
      ptErrRatio60_70->Write();
      ptErrRatio70_80->Write();
      ptErrRatio80_90->Write();
      ptErrRatio90_100->Write();
      ptErrRatio100_200->Write();
      ptErrRatio200_500->Write();
      ptErrRatio500_1000->Write();
      ptRatio30_40->Write();
      ptRatio40_50->Write();
      ptRatio50_60->Write();
      ptRatio60_70->Write();
      ptRatio70_80->Write();
      ptRatio80_90->Write();
      ptRatio90_100->Write();
      ptRatio100_200->Write();
      ptRatio200_500->Write();
      ptRatio500_1000->Write();
      ptRatio1000_up->Write();
      etaRatio1->Write();
      etaRatio2->Write();
      etaRatio3->Write();
      etaRatio4->Write();
      etaRatio5->Write();
      etaRatio6->Write();
      etaRatio7->Write();
      etaRatio8->Write();
      etaRatio9->Write();
      etaRatio10->Write();
      etaRatio11->Write();
      etaRatio12->Write();
      phiRatio1->Write();
      phiRatio2->Write();
      phiRatio3->Write();
      phiRatio4->Write();
      phiRatio5->Write();
      phiRatio6->Write();
      phiRatio7->Write();
      phiRatio8->Write();
      phiRatio9->Write();
      phiRatio10->Write();
      phiRatio11->Write();
      phiRatio12->Write();


      //MET error histograms
      metPtErr->Write();
      metPhiErr->Write();
      metT0PtErr->Write();
      metT0PhiErr->Write();
      metT1PtErr->Write();
      metT1PhiErr->Write();
      metT0T1PtErr->Write();
      metT0T1PhiErr->Write();

      //MET histograms
      metPt->Write();
      metT0Pt->Write();
      metT1Pt->Write();
      metT0T1Pt->Write();
    
      //Error histograms
      phiErr->Write();
      etaErr->Write();
      EErr->Write();
      ptErr->Write();
      ptErr_uncorr->Write();

      rootfile1->Close();

      cout << "Finished item " << entry << " of " << N << endl;
      }



    }
