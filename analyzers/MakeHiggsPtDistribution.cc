#include "MakeHiggsPtDistribution.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>

//ROOT includes
#include "TH1F.h"

using namespace std;


void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }
}


void MakeHiggsPtDistribution::Analyze(bool isData, int option, string outputFileName, string label)
{
  //initialization: create one TTree for each analysis box
  std::cout << "Initializing..." << std::endl;
  if (outputFileName.empty()){
    cout << "MakeHiggsPtDistribution: Output filename not specified!" << endl << "Using default output name MakeHiggsPtDistribution.root" << endl;
    outputFileName = "MakeHiggsPtDistribution.root";
  }
  //---------------------------

  TFile outFile(outputFileName.c_str(), "RECREATE");

  string Label = label;
  if (label != "") Label = "_"+label;
  // TH1F* higgsPt =  new TH1F( ("higgsPt"+Label).c_str(),";higgsPt;Number of Events", 200,0,1600);

  //variable binning
  Float_t bins [] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230,
  240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460,
  470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 625, 650, 675, 700, 750, 800, 1000, 1600};
  int size = sizeof(bins)/sizeof(bins[0])-1;

  TH1F* higgsPt =  new TH1F( ("higgsPt"+Label).c_str(),";higgsPt;Number of Events", size, bins);


  //begin loop
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;


  //cout << "Total Events: " << fChain->GetEntries() << "\n";
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //begin event
    if(jentry % 1000000 == 0) cout << "Processing entry " << jentry << endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //find in-time pileup
    float gHiggsPt = 0.0;
    for (int i=0; i < nGenParticle; i++)
    {
      if (abs(gParticleId[i])== 25)
      {
        gHiggsPt = gParticlePt[i];

      }
    }

    higgsPt->Fill(gHiggsPt, genWeight);

  }//end of event loop

  //Normalize the histograms
  // NormalizeHist(higgsPt);

  cout << "Writing output ..." << endl;
  outFile.WriteTObject(higgsPt, ("higgsPt"+Label).c_str(), "WriteDelete");
  outFile.Close();
}
