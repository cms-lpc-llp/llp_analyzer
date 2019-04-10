
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

const Int_t NComponents = 10;
int color[NComponents] = {kRed, kGreen+2, kBlue, kViolet, kAzure+10, kBlack, kOrange+1, kGray, kBlack, kBlack};


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
void MakeBestBinSummary ( string signalfile, string signalLabel,  vector<string> bkgfiles,vector<string> bkgLabels, vector<double> bkgKFactors, int boxOption = 0, int option = -1, string label = "", string latexlabel = "") {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  double intLumi = 4000; //in units of pb^-1
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
  for (int i=0; i < int(bkgfiles.size()); ++i) {
     inputfiles.push_back(bkgfiles[i]);
     processLabels.push_back(bkgLabels[i]);
  }

  //*******************************************************************************************
  //Bins
  //Assume underflow bin is rejected, but overflow bin is kept
  //*******************************************************************************************
  // const int nBinsMR = 7;
  // const int nBinsRsq = 5;
  // const int nBinsNBTag = 4;
  // double MRBins[nBinsMR] = {400, 500, 600, 700, 900, 1200, 1600};
  // double RsqBins[nBinsRsq] = {0.25,0.30,0.41,0.52,0.64};
  // double NBTagBins[nBinsNBTag] = {0, 1, 2, 3};

  const int nBinsMR = 5;
  const int nBinsRsq = 3;
  const int nBinsNBTag = 4;
  double MRBins[nBinsMR] = {400, 500, 600, 700, 1000};
  double RsqBins[nBinsRsq] = {0.25,0.30,0.40};
  double NBTagBins[nBinsNBTag] = {0, 1, 2, 3};

  //int tmpBin = tmpMRBin + (nBinsMR)*(tmpRsqBin + (nBinsRsq)*tmpNBTagBin);
  int NTotalBins = nBinsMR * nBinsRsq * nBinsNBTag;

  //*******************************************************************************************
  //Define Counts
  //*******************************************************************************************
  vector<vector<double> > BkgYieldByProcess;
  vector<double> TotalBkgYield;
  vector<double> SignalYield;
  vector<vector<double> > BkgYieldErrSqrByProcess;
  vector<double> TotalBkgYieldErrSqr;
  vector<double> SignalYieldErrSqr;


  for (int j=0;j<int(bkgfiles.size()); ++j) {
    vector<double> tmp;
    vector<double> tmpErrSqr;
    for (int i=0; i<NTotalBins;++i) {
      tmp.push_back(0.0);      
      tmpErrSqr.push_back(0.0);      
    }
    BkgYieldByProcess.push_back(tmp);
    BkgYieldErrSqrByProcess.push_back(tmpErrSqr);
  }
  for (int i=0; i<NTotalBins;++i) {
    TotalBkgYield.push_back(0);
    SignalYield.push_back(0);
    TotalBkgYieldErrSqr.push_back(0); 
    SignalYieldErrSqr.push_back(0);
  }
  



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
    float minMTBJet = 0;

    tree->SetBranchAddress("weight",&weight);
    tree->SetBranchAddress("box",&box);
    tree->SetBranchAddress("nBTaggedJets",&nBTaggedJets);
    tree->SetBranchAddress("dPhiRazor",&dPhiRazor);
    tree->SetBranchAddress("MR",&MR);
    tree->SetBranchAddress("Rsq",&Rsq);
    tree->SetBranchAddress("minMTbjet",&minMTBJet);


    cout << "Process : " << processLabels[i] << " : Total Events: " << tree->GetEntries() << "\n";
    for (int n=0;n<tree->GetEntries();n++) { 
    
      tree->GetEntry(n);
      if (n % 1000000 == 0) cout << "Processing Event " << n << "\n";       


      if (intLumi*weight > 100) continue;


      double tmpWeight = intLumi*weight;

      int bkgProcessIndex = -1;
      if (hasSignal) {
	bkgProcessIndex = i-1;
      } else {
	bkgProcessIndex = i;
      }      
      if (!hasSignal || i != 0) {
	tmpWeight = intLumi*weight*bkgKFactors[bkgProcessIndex];
      }

      

      //Box Options
      // if (option == 0 ) {
      // 	if (!(nBTaggedJets == 0)) continue;
      // }
      // if (option == 1 ) {
      // 	if (!(nBTaggedJets >= 1)) continue;
      // }

      if (boxOption == 0) {
	if (box != 8) continue;
      } else if (boxOption == 1) {
	if (box != 7) continue;
      } else if (boxOption == 2) {
	if (box != 3) continue;
      } else if (boxOption == 3) {
	if (box != 5) continue;
      } else if (boxOption == 11) {
	if (!(box == 7 || box == 3 || box == 5)) continue;
      }

      //apply baseline cuts
      if (!(MR > 300 && Rsq > 0.10)) continue;

 
      if (
	  Rsq > 0.25 && MR > 400 
	  && fabs(dPhiRazor) < 2.7
	  && minMTBJet>200
	  ) {

	//determine the bin
	int tmpMRBin = -1;
	for (int k=0; k<nBinsMR; ++k) {
	  if (k<nBinsMR-1) {
	    if (MR >= MRBins[k] && MR < MRBins[k+1]) {
	      tmpMRBin = k; break;
	    }
	  } else {
	    if (MR >= MRBins[k]) {
	      tmpMRBin = k; break;
	    }
	  }
	}
	int tmpRsqBin = -1;
	for (int k=0; k<nBinsRsq; ++k) {
	  if (k<nBinsRsq-1) {
	    if (Rsq >= RsqBins[k] && Rsq < RsqBins[k+1]) {
	      tmpRsqBin = k; break;
	    }
	  } else {
	    if (Rsq >= RsqBins[k]) {
	      tmpRsqBin = k; break;
	    }
	  }
	}
	int tmpNBTagBin = -1;
	for (int k=0; k<nBinsNBTag; ++k) {
	  if (k<nBinsNBTag-1) {
	    if (nBTaggedJets >= NBTagBins[k] && nBTaggedJets < NBTagBins[k+1]) {
	      tmpNBTagBin = k; break;
	    }
	  } else {
	    if (nBTaggedJets >= NBTagBins[k]) {
	      tmpNBTagBin = k; break;
	    }
	  }
	}

	int tmpBin = tmpMRBin + (nBinsMR)*(tmpRsqBin + (nBinsRsq)*tmpNBTagBin);
	
	if (hasSignal && i==0) {
	  SignalYield[tmpBin] += tmpWeight;
	  SignalYieldErrSqr[tmpBin] += tmpWeight*tmpWeight;
	 } else {
	  
	  if (!(tmpBin >= 0 && tmpBin < NTotalBins) ) {
	    cout << "BAD : " << tmpBin << " " << tmpMRBin << " " << tmpRsqBin << " " << tmpNBTagBin << "\n";
	    cout << MR << " " << Rsq << " " << nBTaggedJets << "\n";
	  }

	  if (bkgProcessIndex >= 0 && bkgProcessIndex < int(bkgfiles.size())) {
	    BkgYieldByProcess[bkgProcessIndex][tmpBin] += tmpWeight;
	    BkgYieldErrSqrByProcess[bkgProcessIndex][tmpBin] += tmpWeight*tmpWeight;
	  }
	  TotalBkgYield[tmpBin] += tmpWeight;
	  TotalBkgYieldErrSqr[tmpBin] += tmpWeight*tmpWeight;
	}
	
      } //pass cuts
  

    } //loop over events

    inputFile->Close();
    delete inputFile;

  } //loop over input files
  

  //*******************************************************************************************
  //Sort Bins By S/B or S/sqrt(B)
  //*******************************************************************************************
  vector<int> binIndexMap;
  for (int i=0; i<NTotalBins;++i) binIndexMap.push_back(i);

  const int sortOption = 0;
  for(int i = 0; i < NTotalBins ; i++) {
    for (int j=0; j < NTotalBins-1; j++) {

      double disc_j = 0;
      double disc_jplusone = 0;
      if (sortOption == 0) {
	disc_j = SignalYield[j] / TotalBkgYield[j];
	disc_jplusone = SignalYield[j+1] / TotalBkgYield[j+1];
      } else if (sortOption == 1) {
	disc_j = SignalYield[j] / sqrt(TotalBkgYield[j]);
	disc_jplusone = SignalYield[j+1] / sqrt(TotalBkgYield[j+1]);
      }

      if (disc_jplusone > disc_j) { 

	// swap elements j and j+1
	int tmpBin;
	double tmpSignalYield;
	double tmpSignalYieldErrSqr;
	double tmpTotalBkgYield;
	double tmpTotalBkgYieldErrSqr;
	double tmpBkgYieldByProcess[bkgfiles.size()];
	double tmpBkgYieldErrSqrByProcess[bkgfiles.size()];

	tmpBin = binIndexMap[j] ;
	binIndexMap[j] = binIndexMap[j+1];
	binIndexMap[j+1] = tmpBin;

	tmpSignalYield = SignalYield[j];
	tmpSignalYieldErrSqr = SignalYieldErrSqr[j];
	tmpTotalBkgYield = TotalBkgYield[j];
	tmpTotalBkgYieldErrSqr = TotalBkgYieldErrSqr[j];
	for (int k=0;k<int(bkgfiles.size());++k) {
	  tmpBkgYieldByProcess[k] = BkgYieldByProcess[k][j];
	  tmpBkgYieldErrSqrByProcess[k] = BkgYieldErrSqrByProcess[k][j];
	}

	SignalYield[j] = SignalYield[j+1];
	SignalYieldErrSqr[j] = SignalYieldErrSqr[j+1];
	TotalBkgYield[j] = TotalBkgYield[j+1];
	TotalBkgYieldErrSqr[j] = TotalBkgYieldErrSqr[j+1];
	for (int k=0;k<int(bkgfiles.size());++k) {
	  BkgYieldByProcess[k][j] = BkgYieldByProcess[k][j+1];
	  BkgYieldErrSqrByProcess[k][j] = BkgYieldErrSqrByProcess[k][j+1];
	}

	SignalYield[j+1] = tmpSignalYield;
	SignalYieldErrSqr[j+1] = tmpSignalYieldErrSqr;
	TotalBkgYield[j+1] = tmpTotalBkgYield;
	TotalBkgYieldErrSqr[j+1] = tmpTotalBkgYieldErrSqr;
	for (int k=0;k<int(bkgfiles.size());++k) {
	  BkgYieldByProcess[k][j+1] = tmpBkgYieldByProcess[k];
	  BkgYieldErrSqrByProcess[k][j+1] = tmpBkgYieldErrSqrByProcess[k];
	}

      }
    }
  }
  


  //*******************************************************************************************
  //Summarize Yields
  //*******************************************************************************************
  cout << "Yield Summary \n";
  for( int i=0;i<NTotalBins;i++) {

    int binMR = binIndexMap[i] % nBinsMR;
    int binRsq = int(floor(binIndexMap[i] / nBinsMR)) % nBinsRsq;
    int binNBTag = int(floor(floor(binIndexMap[i] / nBinsMR)/nBinsRsq));

    string BinLabel = "";
    if (binMR < nBinsMR-1) BinLabel += string(Form("MR [%d,%d],",int(MRBins[binMR]),int(MRBins[binMR+1])));
    else BinLabel += string(Form("MR [%d,Inf],",int(MRBins[binMR])));

    if (binRsq < nBinsRsq-1) BinLabel += string(Form("Rsq [%.2f,%.2f],",RsqBins[binRsq],RsqBins[binRsq+1]));
    else BinLabel += string(Form("Rsq [%.2f,Inf],",RsqBins[binRsq]));

    if (binNBTag < nBinsNBTag-1) BinLabel += string(Form("NBtag = %d,",int(NBTagBins[binNBTag])));
    else BinLabel += " NBtag >= 3";


    if (SignalYield[i] > 0.25) {
      //cout << "Bin " << binIndexMap[i] << " : " << binMR << " " << binRsq << " " << binNBTag << " |  " << BinLabel << " |  Signal = " << SignalYield[i] << " +/- " << sqrt(SignalYieldErrSqr[i]) << " , Bkg = " << TotalBkgYield[i] << " +/- " << sqrt(TotalBkgYieldErrSqr[i]) << " | " << SignalYield[i] / TotalBkgYield[i] << " | " << SignalYield[i] / sqrt(TotalBkgYield[i]) << " \n";
      cout << "Bin " << BinLabel << " |  " << Form (" S = %.2f +/- %.2f , B = %.2f +/- %.2f | S/B = %.2f | S/sqrt(B) = %.2f", SignalYield[i],sqrt(SignalYieldErrSqr[i]),  TotalBkgYield[i] ,sqrt(TotalBkgYieldErrSqr[i]),SignalYield[i] / TotalBkgYield[i] , SignalYield[i] / sqrt(TotalBkgYield[i])) << "\n";

      //Signal = " << SignalYield[i] << " +/- " << sqrt(SignalYieldErrSqr[i]) << " , Bkg = " << TotalBkgYield[i] << " +/- " << sqrt(TotalBkgYieldErrSqr[i]) << " | " << SignalYield[i] / TotalBkgYield[i] << " | " << SignalYield[i] / sqrt(TotalBkgYield[i]) << " \n";
    }
  }


 }


 void BestBinSummary() {

   // string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/Inclusive/RazorAnalysis_SMS-T1qqqq_2J_mGl-1400_mLSP-100_25ns_1pb_weighted.root";  
   // string signalLabel = "T1qqqq m_{G}=1400 m_{LSP}=100";
   // string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/new/RazorAnalysis_SMS-T1bbbb_2J_mGl-1500_mLSP-100_25ns_1pb_weighted.root";  
   // string signalLabel = "T1bbbb m_{G}=1500 m_{LSP}=100";
   // string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/newWithElectronD0Cut/Inclusive/RazorAnalysis_SMS-T1tttt_2J_mGl-1500_mLSP-100_25ns_1pb_weighted.root";  
   //  string signalLabel = "T1tttt m_{G}=1500 m_{LSP}=100";
   string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/new/RazorAnalysis_SMS-T2tt_2J_mStop-850_mLSP-100_25ns_1pb_weighted.root";  
   string signalLabel = "T2tt m_{#tilde{t}}=850 m_{LSP}=100";
   //  string signalfile = "/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/new/RazorAnalysis_SMS-T2tt_2J_mStop-425_mLSP-325_25ns_1pb_weighted.root";  
   // string signalLabel = "T2tt m_{#tilde{t}}=425 m_{LSP}=325";

   vector<string> bkgfiles;
   vector<string> bkgLabels;
   vector<double> bkgKFactors;
   
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/new/RazorAnalysis_TTJets_25ns_1pb_weighted.root");  
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/new/RazorAnalysis_DYJetsToLL_HT100ToInf_25ns_1pb_weighted.root");
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/new/RazorAnalysis_WJetsToLNu_HT100ToInf_25ns_1pb_weighted.root");
   bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/new/RazorAnalysis_ZJetsToNuNu_HT100ToInf_25ns_1pb_weighted.root");
   //bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/new/RazorAnalysis_QCD_HT100ToInf_25ns_1pb_weighted.root"); 
   //bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/new/RazorAnalysis_SingleTop_25ns_1pb_weighted.root"); 
   //bkgfiles.push_back("/afs/cern.ch/work/s/sixie/public/Run2SUSY/RazorAnalysis/new/RazorAnalysis_Multiboson_25ns_1pb_weighted.root"); 
   
   bkgLabels.push_back("TTJets");
   bkgLabels.push_back("DYJetsToLL");
   bkgLabels.push_back("WJetsToLNu");
   bkgLabels.push_back("ZJetsToNuNu");
   //bkgLabels.push_back("QCD");
   //bkgLabels.push_back("SingleTop");
   //bkgLabels.push_back("Other");
   
   bkgKFactors.push_back(1.1*689.1/424.5);
   bkgKFactors.push_back(1.1*3.*2008.4/5482.);
   bkgKFactors.push_back(1.1*3.*20508.9/50100.0);
   bkgKFactors.push_back(1.1*3.*2008.4/5482.);
   //bkgKFactors.push_back(1.0);
   //bkgKFactors.push_back(1.1);
   //bkgKFactors.push_back(1.1);
   
   //RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,0,0,"T1qqqq_MultiJet_ZeroBTags", "MultiJet Box 0 b-tag");
   // RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,0,1,"T1bbbb_MultiJet_OneOrMoreBTags","MultiJet Box #geq 1 b-tag");
   //RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,0,1,"T1tttt_MultiJet_OneOrMoreBTags","MultiJet Box #geq 1 b-tag");
   //RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,1,1,"T1tttt_LooseLeptonMultiJet_OneOrMoreBTags","LooseLeptonMultiJet Box #geq 1 b-tag" );
   //RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,2,1,"T1tttt_MuMultiJet_OneOrMoreBTags","MuMultiJet Box #geq 1 b-tag");
   //RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,3,1,"T1tttt_EleMultiJet_OneOrMoreBTags","EleMultiJet Box #geq 1 b-tag");
   //RunMakeRazorPlots(signalfile,signalLabel,bkgfiles,bkgLabels,11,1,"T1tttt_LeptonMultiJet_OneOrMoreBTags","Lepton+MultiJet Box #geq 1 b-tag");

   MakeBestBinSummary(signalfile,signalLabel,bkgfiles,bkgLabels,bkgKFactors,0,1,"T2tt_MultiJet_OneOrMoreBTags","MultiJet Box #geq 1 b-tag");

 
 }
 



//*******************************************************************************
//Preliminary results from PHYS14 samples below
//MR>400, Rsq > 0.25, intL = 4fb^-1 
//*******************************************************************************

//*******************************************************************************
//No minMTbjet cut
//*******************************************************************************
// Bin MR [700,1000],Rsq [0.40,Inf], NBtag >= 3 |   S = 0.31 +/- 0.02 , B = 1.12 +/- 0.29 | S/B = 0.28 | S/sqrt(B) = 0.29
// Bin MR [1000,Inf],Rsq [0.40,Inf],NBtag = 2, |   S = 0.80 +/- 0.02 , B = 3.01 +/- 0.48 | S/B = 0.27 | S/sqrt(B) = 0.46
// Bin MR [1000,Inf],Rsq [0.30,0.40],NBtag = 2, |   S = 0.69 +/- 0.02 , B = 3.48 +/- 0.52 | S/B = 0.20 | S/sqrt(B) = 0.37
// Bin MR [1000,Inf],Rsq [0.25,0.30],NBtag = 2, |   S = 0.49 +/- 0.02 , B = 3.74 +/- 0.55 | S/B = 0.13 | S/sqrt(B) = 0.25
// Bin MR [700,1000],Rsq [0.40,Inf],NBtag = 2, |   S = 1.58 +/- 0.03 , B = 12.21 +/- 1.02 | S/B = 0.13 | S/sqrt(B) = 0.45
// Bin MR [1000,Inf],Rsq [0.40,Inf],NBtag = 1, |   S = 0.96 +/- 0.03 , B = 9.41 +/- 0.59 | S/B = 0.10 | S/sqrt(B) = 0.31
// Bin MR [1000,Inf],Rsq [0.30,0.40],NBtag = 1, |   S = 0.89 +/- 0.03 , B = 11.27 +/- 0.75 | S/B = 0.08 | S/sqrt(B) = 0.27
// Bin MR [700,1000],Rsq [0.40,Inf],NBtag = 1, |   S = 1.90 +/- 0.04 , B = 32.50 +/- 1.18 | S/B = 0.06 | S/sqrt(B) = 0.33
// Bin MR [1000,Inf],Rsq [0.25,0.30],NBtag = 1, |   S = 0.59 +/- 0.02 , B = 10.58 +/- 0.71 | S/B = 0.06 | S/sqrt(B) = 0.18
// Bin MR [600,700],Rsq [0.40,Inf],NBtag = 2, |   S = 0.46 +/- 0.02 , B = 14.32 +/- 1.12 | S/B = 0.03 | S/sqrt(B) = 0.12
// Bin MR [700,1000],Rsq [0.30,0.40],NBtag = 2, |   S = 0.50 +/- 0.02 , B = 17.06 +/- 1.24 | S/B = 0.03 | S/sqrt(B) = 0.12
// Bin MR [600,700],Rsq [0.40,Inf],NBtag = 1, |   S = 0.58 +/- 0.02 , B = 36.87 +/- 1.53 | S/B = 0.02 | S/sqrt(B) = 0.10
// Bin MR [700,1000],Rsq [0.30,0.40],NBtag = 1, |   S = 0.59 +/- 0.02 , B = 41.07 +/- 1.59 | S/B = 0.01 | S/sqrt(B) = 0.09
// Bin MR [1000,Inf],Rsq [0.40,Inf],NBtag = 0, |   S = 0.34 +/- 0.02 , B = 27.04 +/- 0.62 | S/B = 0.01 | S/sqrt(B) = 0.06
// Bin MR [1000,Inf],Rsq [0.30,0.40],NBtag = 0, |   S = 0.28 +/- 0.01 , B = 29.93 +/- 0.92 | S/B = 0.01 | S/sqrt(B) = 0.05
// Bin MR [500,600],Rsq [0.40,Inf],NBtag = 2, |   S = 0.30 +/- 0.01 , B = 35.01 +/- 1.87 | S/B = 0.01 | S/sqrt(B) = 0.05
// Bin MR [700,1000],Rsq [0.25,0.30],NBtag = 1, |   S = 0.31 +/- 0.02 , B = 43.53 +/- 1.72 | S/B = 0.01 | S/sqrt(B) = 0.05
// Bin MR [700,1000],Rsq [0.40,Inf],NBtag = 0, |   S = 0.66 +/- 0.02 , B = 99.62 +/- 1.32 | S/B = 0.01 | S/sqrt(B) = 0.07
// Bin MR [500,600],Rsq [0.40,Inf],NBtag = 1, |   S = 0.40 +/- 0.02 , B = 78.49 +/- 2.49 | S/B = 0.01 | S/sqrt(B) = 0.05

//*******************************************************************************
//minMTbjet > 200 cut
//*******************************************************************************
// Bin MR [700,1000],Rsq [0.40,Inf], NBtag >= 3 |   S = 0.26 +/- 0.01 , B = 0.41 +/- 0.14 | S/B = 0.65 | S/sqrt(B) = 0.41
// Bin MR [1000,Inf],Rsq [0.40,Inf],NBtag = 2, |   S = 0.74 +/- 0.02 , B = 1.87 +/- 0.34 | S/B = 0.40 | S/sqrt(B) = 0.54
// Bin MR [1000,Inf],Rsq [0.30,0.40],NBtag = 2, |   S = 0.63 +/- 0.02 , B = 1.75 +/- 0.32 | S/B = 0.36 | S/sqrt(B) = 0.48
// Bin MR [700,1000],Rsq [0.40,Inf],NBtag = 2, |   S = 1.48 +/- 0.03 , B = 6.42 +/- 0.63 | S/B = 0.23 | S/sqrt(B) = 0.58
// Bin MR [1000,Inf],Rsq [0.25,0.30],NBtag = 2, |   S = 0.43 +/- 0.02 , B = 1.88 +/- 0.36 | S/B = 0.23 | S/sqrt(B) = 0.31
// Bin MR [1000,Inf],Rsq [0.40,Inf],NBtag = 1, |   S = 0.92 +/- 0.03 , B = 7.94 +/- 0.50 | S/B = 0.12 | S/sqrt(B) = 0.33
// Bin MR [1000,Inf],Rsq [0.30,0.40],NBtag = 1, |   S = 0.85 +/- 0.03 , B = 9.07 +/- 0.64 | S/B = 0.09 | S/sqrt(B) = 0.28
// Bin MR [1000,Inf],Rsq [0.25,0.30],NBtag = 1, |   S = 0.55 +/- 0.02 , B = 7.58 +/- 0.51 | S/B = 0.07 | S/sqrt(B) = 0.20
// Bin MR [600,700],Rsq [0.40,Inf],NBtag = 2, |   S = 0.43 +/- 0.02 , B = 5.87 +/- 0.58 | S/B = 0.07 | S/sqrt(B) = 0.18
// Bin MR [700,1000],Rsq [0.40,Inf],NBtag = 1, |   S = 1.84 +/- 0.04 , B = 25.48 +/- 0.95 | S/B = 0.07 | S/sqrt(B) = 0.36
// Bin MR [700,1000],Rsq [0.30,0.40],NBtag = 2, |   S = 0.42 +/- 0.02 , B = 6.20 +/- 0.62 | S/B = 0.07 | S/sqrt(B) = 0.17
// Bin MR [500,600],Rsq [0.40,Inf],NBtag = 2, |   S = 0.25 +/- 0.01 , B = 12.03 +/- 0.96 | S/B = 0.02 | S/sqrt(B) = 0.07
// Bin MR [600,700],Rsq [0.40,Inf],NBtag = 1, |   S = 0.56 +/- 0.02 , B = 28.91 +/- 1.28 | S/B = 0.02 | S/sqrt(B) = 0.10
// Bin MR [700,1000],Rsq [0.30,0.40],NBtag = 1, |   S = 0.54 +/- 0.02 , B = 29.68 +/- 1.24 | S/B = 0.02 | S/sqrt(B) = 0.10
// Bin MR [1000,Inf],Rsq [0.40,Inf],NBtag = 0, |   S = 0.34 +/- 0.02 , B = 27.04 +/- 0.62 | S/B = 0.01 | S/sqrt(B) = 0.06
// Bin MR [1000,Inf],Rsq [0.30,0.40],NBtag = 0, |   S = 0.28 +/- 0.01 , B = 29.93 +/- 0.92 | S/B = 0.01 | S/sqrt(B) = 0.05
// Bin MR [700,1000],Rsq [0.25,0.30],NBtag = 1, |   S = 0.27 +/- 0.01 , B = 31.26 +/- 1.38 | S/B = 0.01 | S/sqrt(B) = 0.05
// Bin MR [700,1000],Rsq [0.40,Inf],NBtag = 0, |   S = 0.66 +/- 0.02 , B = 99.62 +/- 1.32 | S/B = 0.01 | S/sqrt(B) = 0.07
// Bin MR [500,600],Rsq [0.40,Inf],NBtag = 1, |   S = 0.38 +/- 0.02 , B = 58.84 +/- 2.09 | S/B = 0.01 | S/sqrt(B) = 0.05

//*******************************************************************************
//minMTbjet > 200 cut, No Single Top , no ttZ
//*******************************************************************************
// Bin MR [700,1000],Rsq [0.40,Inf], NBtag >= 3 |   S = 0.26 +/- 0.01 , B = 0.29 +/- 0.13 | S/B = 0.91 | S/sqrt(B) = 0.49
// Bin MR [1000,Inf],Rsq [0.30,0.40],NBtag = 2, |   S = 0.63 +/- 0.02 , B = 1.10 +/- 0.19 | S/B = 0.58 | S/sqrt(B) = 0.60
// Bin MR [1000,Inf],Rsq [0.40,Inf],NBtag = 2, |   S = 0.74 +/- 0.02 , B = 1.32 +/- 0.28 | S/B = 0.57 | S/sqrt(B) = 0.65
// Bin MR [1000,Inf],Rsq [0.25,0.30],NBtag = 2, |   S = 0.43 +/- 0.02 , B = 0.90 +/- 0.10 | S/B = 0.48 | S/sqrt(B) = 0.45
// Bin MR [700,1000],Rsq [0.40,Inf],NBtag = 2, |   S = 1.48 +/- 0.03 , B = 4.66 +/- 0.49 | S/B = 0.32 | S/sqrt(B) = 0.68
// Bin MR [1000,Inf],Rsq [0.40,Inf],NBtag = 1, |   S = 0.92 +/- 0.03 , B = 7.23 +/- 0.43 | S/B = 0.13 | S/sqrt(B) = 0.34
// Bin MR [1000,Inf],Rsq [0.30,0.40],NBtag = 1, |   S = 0.85 +/- 0.03 , B = 7.75 +/- 0.52 | S/B = 0.11 | S/sqrt(B) = 0.31
// Bin MR [600,700],Rsq [0.40,Inf],NBtag = 2, |   S = 0.43 +/- 0.02 , B = 4.16 +/- 0.46 | S/B = 0.10 | S/sqrt(B) = 0.21
// Bin MR [700,1000],Rsq [0.30,0.40],NBtag = 2, |   S = 0.42 +/- 0.02 , B = 4.35 +/- 0.47 | S/B = 0.10 | S/sqrt(B) = 0.20
// Bin MR [1000,Inf],Rsq [0.25,0.30],NBtag = 1, |   S = 0.55 +/- 0.02 , B = 6.75 +/- 0.45 | S/B = 0.08 | S/sqrt(B) = 0.21
// Bin MR [700,1000],Rsq [0.40,Inf],NBtag = 1, |   S = 1.84 +/- 0.04 , B = 23.39 +/- 0.85 | S/B = 0.08 | S/sqrt(B) = 0.38
// Bin MR [500,600],Rsq [0.40,Inf],NBtag = 2, |   S = 0.25 +/- 0.01 , B = 9.14 +/- 0.81 | S/B = 0.03 | S/sqrt(B) = 0.08
// Bin MR [600,700],Rsq [0.40,Inf],NBtag = 1, |   S = 0.56 +/- 0.02 , B = 26.29 +/- 1.20 | S/B = 0.02 | S/sqrt(B) = 0.11
// Bin MR [700,1000],Rsq [0.30,0.40],NBtag = 1, |   S = 0.54 +/- 0.02 , B = 26.93 +/- 1.15 | S/B = 0.02 | S/sqrt(B) = 0.10
// Bin MR [1000,Inf],Rsq [0.40,Inf],NBtag = 0, |   S = 0.34 +/- 0.02 , B = 26.53 +/- 0.57 | S/B = 0.01 | S/sqrt(B) = 0.07
// Bin MR [1000,Inf],Rsq [0.30,0.40],NBtag = 0, |   S = 0.28 +/- 0.01 , B = 29.42 +/- 0.88 | S/B = 0.01 | S/sqrt(B) = 0.05
// Bin MR [700,1000],Rsq [0.25,0.30],NBtag = 1, |   S = 0.27 +/- 0.01 , B = 28.50 +/- 1.28 | S/B = 0.01 | S/sqrt(B) = 0.05
// Bin MR [500,600],Rsq [0.40,Inf],NBtag = 1, |   S = 0.38 +/- 0.02 , B = 53.88 +/- 1.98 | S/B = 0.01 | S/sqrt(B) = 0.05
// Bin MR [700,1000],Rsq [0.40,Inf],NBtag = 0, |   S = 0.66 +/- 0.02 , B = 98.09 +/- 1.26 | S/B = 0.01 | S/sqrt(B) = 0.07
