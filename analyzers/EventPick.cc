#include "EventPick.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>

//ROOT includes
#include "TH1F.h"
#include "TH2F.h"

#include "SimpleTable.h"

using namespace std;

bool getFileContent(std::string fileName, std::vector<unsigned long int> & vecOfStrs)
{

	// Open the File
	std::ifstream in(fileName.c_str());

	// Check if object is valid
	if(!in)
	{
		std::cerr << "Cannot open the File : "<<fileName<<std::endl;
		return false;
	}

	std::string str;
	// Read the next line from File untill it reaches the end.
	while (std::getline(in, str))
	{
		// Line contains string of length > 0 then save it in vector
		if(str.size() > 0)
			vecOfStrs.push_back(std::stoul(str));
	}
	//Close The File
	in.close();
	return true;
}


void EventPick::Analyze(bool isData, int option, string outputFileName, string label)
{
  //initialization: create one TTree for each analysis box
  std::cout << "Initializing..." << std::endl;
  if (outputFileName.empty()){
    cout << "EventPick: Output filename not specified!" << endl << "Using default output name EventPick.root" << endl;
    outputFileName = "EventPick.root";
  }
  //---------------------------

  TFile outFile(outputFileName.c_str(), "RECREATE");
  TTree *newtree = fChain->CloneTree(0);




  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

  fChain->SetBranchAddress("cscRechitClusterMuonVetoPt",cscRechitClusterMuonVetoPt);
  fChain->SetBranchAddress("cscRechitClusterJetVetoPt",cscRechitClusterJetVetoPt);
  fChain->SetBranchAddress("cscRechitClusterEta",cscRechitClusterEta);
  fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus11",cscRechitClusterNRechitChamberPlus11);
  fChain->SetBranchAddress("cscRechitClusterNRechitChamberPlus12",cscRechitClusterNRechitChamberPlus12);
  fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus11",cscRechitClusterNRechitChamberMinus11);
  fChain->SetBranchAddress("cscRechitClusterNRechitChamberMinus12",cscRechitClusterNRechitChamberMinus12);
  fChain->SetBranchAddress("cscRechitClusterTime",cscRechitClusterTime);
  fChain->SetBranchAddress("gLLP_decay_vertex_x",gLLP_decay_vertex_x);
  fChain->SetBranchAddress("gLLP_decay_vertex_y",gLLP_decay_vertex_y);

  fChain->SetBranchAddress("gLLP_decay_vertex_z",gLLP_decay_vertex_z);
  fChain->SetBranchAddress("nCscRechitClusters",&nCscRechitClusters);
  fChain->SetBranchAddress("eventNum",&eventNum);

  //load events to pick
  std::vector<unsigned long int> eventList;
  string pathname = getenv("CMSSW_BASE");
  pathname = Form("%s/src/llp_analyzer/%s", pathname.c_str(), label.c_str());
  bool result = getFileContent(pathname, eventList);
  if(!result) cout<<"eventPick file not found: "<<pathname<<endl;



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
    bool found = false;
    for(unsigned int i = 0;i<eventList.size();i++)
    {
      if (eventList[i] == eventNum) found = true;
    }
    // if (isData)
    // {
    if (!found) continue;
    // }
    // else
    // {
    //   float gLLP_decay_vertex_r= sqrt(gLLP_decay_vertex_x[0]*gLLP_decay_vertex_x[0]+gLLP_decay_vertex_y[0]*gLLP_decay_vertex_y[0]);
    //   float gLLP_decay_vertex_d= sqrt(gLLP_decay_vertex_r*gLLP_decay_vertex_r+gLLP_decay_vertex_z[0]*gLLP_decay_vertex_z[0]);
    //   if(!(gLLP_decay_vertex_r>800 || abs(gLLP_decay_vertex_z[0])>1200 || gLLP_decay_vertex_d <200)) continue;
		//
    //   gLLP_decay_vertex_r = sqrt(gLLP_decay_vertex_x[1]*gLLP_decay_vertex_x[1]+gLLP_decay_vertex_y[1]*gLLP_decay_vertex_y[1]);
    //   gLLP_decay_vertex_d= sqrt(gLLP_decay_vertex_r*gLLP_decay_vertex_r+gLLP_decay_vertex_z[1]*gLLP_decay_vertex_z[1]);
    //   if(!(gLLP_decay_vertex_r>800 || abs(gLLP_decay_vertex_z[1])>1200 ||gLLP_decay_vertex_d<200)) continue;
		//
    //   int nCluster = 0;
    //   for (int i = 0;i<nCscRechitClusters;i++)
    //   {
    //     if(cscRechitClusterNRechitChamberPlus11[i]!=0) continue;
    //     if(cscRechitClusterNRechitChamberPlus12[i]!=0) continue;
    //     if(cscRechitClusterNRechitChamberMinus11[i]!=0) continue;
    //     if(cscRechitClusterNRechitChamberMinus12[i]!=0) continue;
    //     if(cscRechitClusterMuonVetoPt[i]>=20) continue;
    //     if(cscRechitClusterJetVetoPt[i]>=10) continue;
    //     if(abs(cscRechitClusterEta[i])>=2.1) continue;
    //     if(cscRechitClusterTime[i]>=12.5) continue;
    //     if(cscRechitClusterTime[i]<=-5) continue;
    //     nCluster++;
    //   }
    //   if (nCluster==0)continue;
    // }


    newtree->Fill();
    if (isData) NEvents->Fill(1);
    else NEvents->Fill(1,genWeight);
  }//end of event loop

  cout << "Writing output ..." << endl;
  // cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
  outFile.cd();
  newtree->Write();
  NEvents->Write();
  outFile.Close();

}
