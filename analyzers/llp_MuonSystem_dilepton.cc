#include "llp_MuonSystem_dilepton.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

void llp_MuonSystem_dilepton::Analyze(bool isData, int Option, string outputfilename, string label)
{

    cout << "Initializing..." << endl;
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "HHTo4BNtuple.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");

    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

    //output TTree
    TTree *outputTree = new TTree("tree", "");

    //------------------------
    //declare branch variables
    //------------------------
    float weight = 0;
    int lep1Id = 0;
    float lep1Pt = -1;
    float lep1Eta = -1;
    float lep1Phi = -1;
    int lep2Id = 0;
    float lep2Pt = -1;
    float lep2Eta = -1;
    float lep2Phi = -1;
    float dileptonMass = -1;
    float dileptonPt = -1.;
    float MET = -1;
    int NJets = 0;
    int NBTags = 0;
    float gLepPt[2];
    float gLepPhi[2];
    float gLepE[2];
    float gLepEta[2];
    float gZmass = -1.;
    float gZPt = -1.;


    //------------------------
    //set branches on big tree
    //------------------------
    outputTree->Branch("weight", &weight, "weight/F");
    outputTree->Branch("lep1Id", &lep1Id, "lep1Id/I");
    outputTree->Branch("lep1Pt", &lep1Pt, "lep1Pt/F");
    outputTree->Branch("lep1Eta", &lep1Eta, "lep1Eta/F");
    outputTree->Branch("lep1Phi", &lep1Phi, "lep1Phi/F");
    outputTree->Branch("lep2Id", &lep2Id, "lep2Id/I");
    outputTree->Branch("lep2Pt", &lep2Pt, "lep2Pt/F");
    outputTree->Branch("lep2Eta", &lep2Eta, "lep2Eta/F");
    outputTree->Branch("lep2Phi", &lep2Phi, "lep2Phi/F");
    outputTree->Branch("dileptonMass", &dileptonMass, "dileptonMass/F");
    outputTree->Branch("dileptonPt", &dileptonPt, "dileptonPt/F");

    outputTree->Branch("MET", &MET, "MET/F");
    outputTree->Branch("NJets", &NJets, "NJets/I");
    outputTree->Branch("NBTags", &NBTags, "NBTags/I");

    outputTree->Branch("gLepPt", &gLepPt, "gLepPt[2]/F");
    outputTree->Branch("gLepPhi", &gLepPhi, "gLepPhi[2]/F");
    outputTree->Branch("gLepE", &gLepE, "gLepE[2]/F");
    outputTree->Branch("gLepEta", &gLepEta, "gLepEta[2]/F");
    outputTree->Branch("gZmass", &gZmass, "gZmass/F");
    outputTree->Branch("gZPt", &gZPt, "gZPt/F");

  //   outputTree->Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",                                        &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,                                       "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ/O");
  // outputTree->Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",                                        &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,                                       "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL/O");
  // outputTree->Branch("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",                                        &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,                                       "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL/O");
  // outputTree->Branch("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",                                        &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,                                       "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ/O");
  // outputTree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",                                        &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,                                       "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ/O");
  // outputTree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",                                        &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,                                       "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL/O");
    cout << "Run With Option = " << Option << "\n";

    UInt_t NEventsFilled = 0;

    //begin loop
    if (fChain == 0) return;
    UInt_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;

    cout << "nentries = " << nentries << "\n";
    for (UInt_t jentry=0; jentry<nentries;jentry++) {
      //begin event
      if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //fill normalization histogram
      weight = genWeight / fabs(genWeight);
      NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + weight);


      //reset tree variables
      lep1Id = 0;
      lep1Pt = -1;
      lep1Eta = -1;
      lep1Phi = -1;
      lep2Id = 0;
      lep2Pt = -1;
      lep2Eta = -1;
      lep2Phi = -1;
      dileptonMass = -1;
      dileptonPt = -1.;
      MET = -1;
      NJets = 0;
      NBTags = 0;
      for(int i = 0; i <2;i++)
      {
        gLepPt[i] = 0.0;
        gLepPhi[i] = 0.0;
        gLepE[i] = 0.0;
        gLepEta[i] = 0.0;
      }
      gZPt = -999.;
      gZmass = -999.;


      int count = 0;
      for (int i=0; i < nGenParticle; i++)
      {
        if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == 23)
        {
          gLepPt[count] = gParticlePt[i];
          gLepPhi[count] = gParticlePhi[i];
          gLepE[count] = gParticleE[i];
          gLepEta[count] = gParticleEta[i];
          count++;


        }
        else if (abs(gParticleId[i]) == 23)
        {
          TLorentzVector genZ; genZ.SetPtEtaPhiE(gParticlePt[i], gParticleEta[i], gParticlePhi[i], gParticleE[i]);
          gZmass = genZ.M();
          gZPt = gParticlePt[i];

        }

      }





      //***********************
      //Select Leptons
      //***********************
      vector<TLorentzVector> GoodLeptons;
      vector<int> GoodLeptonId;
      vector<float> GoodLeptonPt;
      vector<float> GoodLeptonEta;
      vector<float> GoodLeptonPhi;

      for(int i = 0; i < nMuons; i++){
      	if(muonPt[i] < 20) continue;
      	if(abs(muonEta[i]) > 2.4) continue;
        if(!isMuonPOGLooseMuon(i)) continue;
        if(!isMuonPOGTightMuon(i, true, true)) continue;

      	// if (!( muonIsTight[i] && muonpfRelIso03_all[i] < 0.25)) continue;



      	GoodLeptonId.push_back(13 * muonCharge[i]);
      	GoodLeptonPt.push_back(muonPt[i]);
      	GoodLeptonEta.push_back(muonEta[i]);
      	GoodLeptonPhi.push_back(muonPhi[i]);
      	TLorentzVector thisMuon; thisMuon.SetPtEtaPhiM(muonPt[i], muonEta[i], muonPhi[i], 0.106);
      	GoodLeptons.push_back(thisMuon);
      }

      for(int i = 0; i < nElectrons;  i++){
      	if(elePt[i] < 20) continue;
      	if(abs(eleEta[i]) > 2.5) continue;
      	// if (!( elemvaFall17V2Iso_WP80[i] )) continue;
        if (!isEGammaPOGLooseElectron(i, true, true, true, "Summer16")) continue;

      	GoodLeptonId.push_back(11 * eleCharge[i]);
      	GoodLeptonPt.push_back(elePt[i]);
      	GoodLeptonEta.push_back(eleEta[i]);
      	GoodLeptonPhi.push_back(elePhi[i]);
      	TLorentzVector thisElectron; thisElectron.SetPtEtaPhiM(elePt[i], eleEta[i], elePhi[i], 0.0005);
      	GoodLeptons.push_back(thisElectron);
      }

      TLorentzVector lep1;
      TLorentzVector lep2;
      for (int i=0; i<GoodLeptonId.size();i++) {

      	if (i==0) {
      	  lep1Id = GoodLeptonId[i];
      	  lep1Pt = GoodLeptonPt[i];
      	  lep1Eta = GoodLeptonEta[i];
      	  lep1Phi = GoodLeptonPhi[i];
      	  lep1 = GoodLeptons[i];
      	}
      	if (i==1) {
      	  lep2Id = GoodLeptonId[i];
      	  lep2Pt = GoodLeptonPt[i];
      	  lep2Eta = GoodLeptonEta[i];
      	  lep2Phi = GoodLeptonPhi[i];
      	  lep2 = GoodLeptons[i];
      	}

      }

      dileptonMass = (lep1+lep2).M();
      dileptonPt = (lep1+lep2).Pt();
      MET = metType1Pt;


      //***********************
      //Fill Event
      //***********************
      NEventsFilled++;
      outputTree->Fill();
    }//end of event loop

    cout << "Filled Total of " << NEventsFilled << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}
