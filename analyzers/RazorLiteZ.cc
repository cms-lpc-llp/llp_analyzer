#include "RazorLiteZ.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"

using namespace std;

struct greater_than_pt{
    inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){
        return p1.Pt() > p2.Pt();
    }
};

class RazorLiteTree {
  
public:
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  matchGen;
  UInt_t  category;
  UInt_t  npv, npu;
  Float_t rho, weight;
  Float_t met, metPhi;
  Int_t   q1, q2;

  TLorentzVector *dilep=0, *lep1=0, *lep2=0;

  Float_t activity1, activity2;
  Float_t pfChIso1, pfGamIso1, pfNeuIso1, pfCombIso1, pfChIso2, pfGamIso2, pfNeuIso2, pfCombIso2;
  Float_t chMiniIso1, chMiniIso2, neuMiniIso1, neuMiniIso2, puChMiniIso1, puChMiniIso2;

  Float_t d01, dz1, d02, dz2;
  UInt_t typeBits1, typeBits2;
  UInt_t isTight2, isMedium2, isLoose2, isVeto2;

  ///// electron specific /////
  Float_t sigieie1, hovere1, eoverp1, sigieie2, hovere2, eoverp2;
  Float_t dphi1, deta1, dphi2, deta2;
  UInt_t  isConv1, nexphits1, isConv2, nexphits2; 
  //TLorentzVector *sc1=0, *sc2=0;

  ///// muon specific /////
  UInt_t muontype1, muontype2, muonquality1, muonquality2;
  Float_t muon_ip3d1, muon_ip3dSignificance1;
  Float_t muon_ip3d2, muon_ip3dSignificance2;
  //TLorentzVector *sta1=0, *sta2=0;

  TTree *tree_;
  TFile *f_;

  RazorLiteTree() {
    InitVariables();
  };

  ~RazorLiteTree() {
    if (f_) f_->Close();
  };

  void InitVariables() {
    runNum=0; lumiSec=0; evtNum=0; matchGen=0; category=0;
    npv=0; npu=0; rho=-1; weight=-1;
    met=-1; metPhi=-1;
    q1=0; q2=0;

    dilep=0; lep1=0; lep2=0;

    activity1=-1; activity2=-1;
    pfChIso1=-1; pfGamIso1=-1; pfNeuIso1=-1; pfCombIso1=-1; pfChIso2=-1; pfGamIso2=-1; pfNeuIso2=-1; pfCombIso2=-1;
    chMiniIso1=-1; chMiniIso2=-1; neuMiniIso1=-1; neuMiniIso2=-1; puChMiniIso1=-1; puChMiniIso2=-1;
    d01=-1; dz1=-1; d02=-1; dz2=-1;

    isTight2=0; isMedium2=0; isLoose2=0; isVeto2=0;

    sigieie1=-1; hovere1=-1; eoverp1=-1; sigieie2=-1; hovere2=-1; eoverp2=-1;
    dphi1=-1; deta1=-1; dphi2=-1; deta2=-1;
    isConv1=0; nexphits1=0; isConv2=0; nexphits2=0;

    muontype1=0; muontype2=0; muonquality1=0; muonquality2=0;
    muon_ip3d1=0; muon_ip3dSignificance1=0;
    muon_ip3d2=0; muon_ip3dSignificance2=0;
    
  };

  void LoadTree(const char* file) {
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->Get("TP"));
    InitTree();
    assert(tree_);

  }

  void CreateTree() {
    tree_ = new TTree("TP","TP");
    f_ = 0;
    
    tree_->Branch("runNum",      &runNum,     "runNum/i");      // event run number
    tree_->Branch("lumiSec",     &lumiSec,    "lumiSec/i");     // event lumi section
    tree_->Branch("evtNum",      &evtNum,     "evtNum/i");      // event number
    tree_->Branch("matchGen",    &matchGen,   "matchGen/i");    // event has both leptons matched to MC Z->ll
    tree_->Branch("category",    &category,   "category/i");    // dilepton category
    tree_->Branch("npv",         &npv,        "npv/i");         // number of primary vertices
    tree_->Branch("npu",         &npu,        "npu/i");         // number of in-time PU events (MC)
    tree_->Branch("weight",      &weight,     "weight/F");
    tree_->Branch("rho",         &rho,        "rho/F");
    tree_->Branch("met",         &met,        "met/F");         // MET
    tree_->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)
    tree_->Branch("q1",          &q1,         "q1/I");          // charge of tag lepton
    tree_->Branch("q2",          &q2,         "q2/I");          // charge of probe lepton
    tree_->Branch("dilep",       "TLorentzVector", &dilep);     // di-lepton 4-vector
    tree_->Branch("lep1",        "TLorentzVector", &lep1);      // tag lepton 4-vector
    tree_->Branch("lep2",        "TLorentzVector", &lep2);      // probe lepton 4-vector
    tree_->Branch("activity1",     &activity1,     "activity1/F");
    tree_->Branch("activity2",     &activity2,     "activity2/F");
    tree_->Branch("pfChIso1",    &pfChIso1,    "pfChIso1/F");      // PF charged hadron isolation of tag lepton
    tree_->Branch("pfChIso2",    &pfChIso2,    "pfChIso2/F");      // PF charged hadron isolation of probe lepton
    tree_->Branch("pfGamIso1",   &pfGamIso1,   "pfGamIso1/F");     // PF photon isolation of tag lepton
    tree_->Branch("pfGamIso2",   &pfGamIso2,   "pfGamIso2/F");     // PF photon isolation of probe lepton
    tree_->Branch("pfNeuIso1",   &pfNeuIso1,   "pfNeuIso1/F");     // PF neutral hadron isolation of tag lepton
    tree_->Branch("pfNeuIso2",   &pfNeuIso2,   "pfNeuIso2/F");     // PF neutral hadron isolation of probe lepton
    tree_->Branch("pfCombIso1",  &pfCombIso1,  "pfCombIso1/F");    // PF combined isolation of tag lepton
    tree_->Branch("pfCombIso2",  &pfCombIso2,  "pfCombIso2/F");    // PF combined isolation of probe lepton    
    tree_->Branch("chMiniIso1",  &chMiniIso1,  "chMiniIso1/F");
    tree_->Branch("chMiniIso2",  &chMiniIso2,  "chMiniIso2/F");
    tree_->Branch("neuMiniIso1",  &neuMiniIso1,  "neuMiniIso1/F");
    tree_->Branch("neuMiniIso2",  &neuMiniIso2,  "neuMiniIso2/F");
    tree_->Branch("puChMiniIso1",  &puChMiniIso1,  "puChMiniIso1/F");
    tree_->Branch("puChMiniIso2",  &puChMiniIso2,  "puChMiniIso2/F");

    tree_->Branch("isTight2",         &isTight2,         "isTight2/i");
    tree_->Branch("isMedium2",         &isMedium2,         "isMedium2/i");
    tree_->Branch("isLoose2",         &isLoose2,         "isLoose2/i");
    tree_->Branch("isVeto2",         &isVeto2,         "isVeto2/i");
    tree_->Branch("d01",         &d01,         "d01/F");           // transverse impact parameter of tag lepton
    tree_->Branch("d02",         &d02,         "d02/F");           // transverse impact parameter of probe lepton 
    tree_->Branch("dz1",         &dz1,         "dz1/F");           // longitudinal impact parameter of tag lepton
    tree_->Branch("dz2",         &dz2,         "dz2/F");           // longitudinal impact parameter of probe lepton 

    tree_->Branch("muontype1",         &muontype1,         "muontype1/i");
    tree_->Branch("muontype2",         &muontype2,         "muontype2/i");
    tree_->Branch("muonquality1",         &muonquality1,         "muonquality1/i");
    tree_->Branch("muonquality2",         &muonquality2,         "muonquality2/i");

    tree_->Branch("muon_ip3d1", &muon_ip3d1, "muon_ip3d1/F");   
    tree_->Branch("muon_ip3d2", &muon_ip3d2, "muon_ip3d2/F");   
    tree_->Branch("muon_ip3dSignificance1", &muon_ip3dSignificance1, "muon_ip3dSignificance1/F");   
    tree_->Branch("muon_ip3dSignificance2", &muon_ip3dSignificance2, "muon_ip3dSignificance2/F");   

    tree_->Branch("sigieie1",   &sigieie1,   "sigieie1/F");    // sigma-ieta-ieta of tag
    tree_->Branch("sigieie2",   &sigieie2,   "sigieie2/F");    // sigma-ieta-ieta of probe
    tree_->Branch("hovere1",    &hovere1,    "hovere1/F");     // H/E of tag
    tree_->Branch("hovere2",    &hovere2,    "hovere2/F");     // H/E of probe
    tree_->Branch("eoverp1",    &eoverp1,    "eoverp1/F");     // E/p of tag
    tree_->Branch("eoverp2",    &eoverp2,    "eoverp2/F");     // E/p of probe 
    tree_->Branch("dphi1",      &dphi1,      "dphi1/F");       // GSF track - ECAL dphi of tag
    tree_->Branch("dphi2",      &dphi2,      "dphi2/F");       // GSF track - ECAL dphi of probe 
    tree_->Branch("deta1",      &deta1,      "deta1/F");       // GSF track - ECAL deta of tag
    tree_->Branch("deta2",      &deta2,      "deta2/F");       // GSF track - ECAL deta of probe
    tree_->Branch("isConv1",    &isConv1,    "isConv1/i");     // conversion filter flag of tag lepton
    tree_->Branch("isConv2",    &isConv2,    "isConv2/i");     // conversion filter flag of probe lepton
    tree_->Branch("nexphits1",  &nexphits1,  "nexphits1/i");   // number of missing expected inner hits of tag lepton
    tree_->Branch("nexphits2",  &nexphits2,  "nexphits2/i");   // number of missing expected inner hits of probe lepton

  };
  
  void InitTree() {

    assert(tree_);
    InitVariables();

    tree_->SetBranchAddress("runNum",      &runNum);
    tree_->SetBranchAddress("lumiSec",     &lumiSec);
    tree_->SetBranchAddress("evtNum",      &evtNum);
    tree_->SetBranchAddress("matchGen",    &matchGen);
    tree_->SetBranchAddress("category",    &category);
    tree_->SetBranchAddress("npv",         &npv);
    tree_->SetBranchAddress("npu",         &npu);
    tree_->SetBranchAddress("weight",      &weight);
    tree_->SetBranchAddress("rho",         &rho);
    tree_->SetBranchAddress("met",         &met);
    tree_->SetBranchAddress("metPhi",      &metPhi);
    tree_->SetBranchAddress("q1",          &q1);
    tree_->SetBranchAddress("q2",          &q2);
    tree_->SetBranchAddress("dilep",       &dilep);     // di-lepton 4-vector
    tree_->SetBranchAddress("lep1",        &lep1);      // tag lepton 4-vector
    tree_->SetBranchAddress("lep2",        &lep2);      // probe lepton 4-vector
    tree_->SetBranchAddress("activity1",     &activity1);
    tree_->SetBranchAddress("activity2",     &activity2);
    tree_->SetBranchAddress("pfChIso1",    &pfChIso1);
    tree_->SetBranchAddress("pfChIso2",    &pfChIso2);
    tree_->SetBranchAddress("pfGamIso1",   &pfGamIso1);
    tree_->SetBranchAddress("pfGamIso2",   &pfGamIso2);
    tree_->SetBranchAddress("pfNeuIso1",   &pfNeuIso1);
    tree_->SetBranchAddress("pfNeuIso2",   &pfNeuIso2);
    tree_->SetBranchAddress("pfCombIso1",  &pfCombIso1);
    tree_->SetBranchAddress("pfCombIso2",  &pfCombIso2);
    tree_->SetBranchAddress("chMiniIso1",  &chMiniIso1);
    tree_->SetBranchAddress("chMiniIso2",  &chMiniIso2);
    tree_->SetBranchAddress("neuMiniIso1",  &neuMiniIso1);
    tree_->SetBranchAddress("neuMiniIso2",  &neuMiniIso2);
    tree_->SetBranchAddress("puChMiniIso1",  &puChMiniIso1);
    tree_->SetBranchAddress("puChMiniIso2",  &puChMiniIso2);

    tree_->SetBranchAddress("isTight2",         &isTight2);
    tree_->SetBranchAddress("isMedium2",         &isMedium2);
    tree_->SetBranchAddress("isLoose2",         &isLoose2);
    tree_->SetBranchAddress("isVeto2",         &isVeto2);
    tree_->SetBranchAddress("d01",         &d01);
    tree_->SetBranchAddress("d02",         &d02);
    tree_->SetBranchAddress("dz1",         &dz1);
    tree_->SetBranchAddress("dz2",         &dz2);

    tree_->Branch("muontype1",         &muontype1);
    tree_->Branch("muontype2",         &muontype2);
    tree_->Branch("muonquality1",         &muonquality1);
    tree_->Branch("muonquality2",         &muonquality2);
    tree_->SetBranchAddress("muon_ip3d1", &muon_ip3d1);
    tree_->SetBranchAddress("muon_ip3d2", &muon_ip3d2);
    tree_->SetBranchAddress("muon_ip3dSignificance1", &muon_ip3dSignificance1);
    tree_->SetBranchAddress("muon_ip3dSignificance2", &muon_ip3dSignificance2);

    tree_->SetBranchAddress("sigieie1",   &sigieie1);
    tree_->SetBranchAddress("sigieie2",   &sigieie2);
    tree_->SetBranchAddress("hovere1",    &hovere1);
    tree_->SetBranchAddress("hovere2",    &hovere2);
    tree_->SetBranchAddress("eoverp1",    &eoverp1);
    tree_->SetBranchAddress("eoverp2",    &eoverp2);
    tree_->SetBranchAddress("dphi1",      &dphi1);
    tree_->SetBranchAddress("dphi2",      &dphi2);
    tree_->SetBranchAddress("deta1",      &deta1);
    tree_->SetBranchAddress("deta2",      &deta2);
    tree_->SetBranchAddress("isConv1",    &isConv1);
    tree_->SetBranchAddress("isConv2",    &isConv2);
    tree_->SetBranchAddress("nexphits1",  &nexphits1);
    tree_->SetBranchAddress("nexphits2",  &nexphits2);
    
  };
};

void RazorLiteZ::Analyze(bool isData, int option, string outputfilename, string label)
{
    //initialization: create one TTree for each analysis box 
    cout << "Initializing..." << endl;    
    cout << "IsData = " << isData << "\n";

    Float_t ELE_MASS = 0.000511;
    Float_t MU_MASS  = 0.105658;

    //TRandom3 *random = new TRandom3(33333); //Artur wants this number 33333

    //*************************************************************************
    //Set up Output File
    //*************************************************************************
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "RazorLite.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    RazorLiteTree *TPPair = new RazorLiteTree;
    TPPair->CreateTree();  
    TPPair->tree_->SetAutoFlush(0);

    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  
    //*************************************************************************
    //Look over Input File Events
    //*************************************************************************
    if (fChain == 0) return;
    cout << "Total Events: " << fChain->GetEntries() << "\n";
    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {
      
      //begin event
      if(jentry % 1000 == 0) cout << "Processing entry " << jentry << endl;
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      //fill normalization histogram

      if (isData) {
	NEvents->Fill(1);
	TPPair->weight = 1;
      }
      else {
	NEvents->Fill(genWeight);
	TPPair->weight = genWeight;     
      }
      //event info
      TPPair->runNum = runNum;
      TPPair->lumiSec = lumiNum;
      TPPair->evtNum = eventNum;
	
      //get NPU
      for (int i=0; i < nBunchXing; ++i) {
	if (BunchXing[i] == 0) {
	  TPPair->npu = nPUmean[i];
	}
      }
      TPPair->npv = nPV;
      TPPair->rho = fixedGridRhoFastjetAll;
      TPPair->met = metType1Pt;
      TPPair->metPhi = metType1Phi;

      //******************************************
      //Find Generated leptons
      //******************************************
      vector<int> genLeptonIndex;      
      //find gen electrons
      if ( option == 1 ) {
	for(int j = 0; j < nGenParticle; j++){

	  //look for electrons
	  if (abs(gParticleId[j]) == 11 && gParticleStatus[j] == 1 	      
	      && abs(gParticleEta[j]) < 3.0 && gParticlePt[j] > 3
	      ) {
	    if ( abs(gParticleMotherId[j]) == 23 )  {
	      genLeptonIndex.push_back(j);
	    }
	  }
	} //loop over gen particles
      }
      
      //look for muons
      if (option == 2) {
	for(int j = 0; j < nGenParticle; j++){
	  if (abs(gParticleId[j]) == 13 && gParticleStatus[j] == 1  
	      && abs(gParticleEta[j]) < 3.0 && gParticlePt[j] > 3
	      ) {
	    if ( abs(gParticleMotherId[j]) == 23 ) {	     
	      genLeptonIndex.push_back(j);
	    }	    
	  }
	} //loop over gen particles		
      }
	

      //*********************************************************
      //Electrons
      //*********************************************************
      if (option == 1) {
	//*******************************************************
	//Loop over Tag electrons
	//*******************************************************
	for(int indexTag = 0; indexTag < nElectrons; indexTag++){
	
	  if(elePt[indexTag] < 30) continue;
	  if(fabs(eleEta[indexTag]) > 2.5) continue;

	  //For MC, Match to Gen level electron
	  if (!isData) {
	    bool genmatch = false;
	    for (int q=0;q<int(genLeptonIndex.size()); q++) {
	      if ( deltaR(eleEta[indexTag],elePhi[indexTag],
			  gParticleEta[genLeptonIndex[q]], 
			  gParticlePhi[genLeptonIndex[q]]) < 0.1) {
		genmatch = true;
	      }
	    }
	    if (!genmatch) continue;
	  }

	  //tag must pass tight cuts
	  if ( !isEGammaPOGMediumElectron(indexTag) ) continue;

	  //Tag must match single electron HLT Filters OR tag leg filter of the dedicated T&P trigger
          //Trigger: HLT_Ele27_WPTight_Gsf_v* (filter 6)
	  if ( !(ele_passHLTFilter[indexTag][6] && HLTDecision[35]) ) continue;

	  TLorentzVector vtag;
	  vtag.SetPtEtaPhiM(elePt[indexTag], eleEta[indexTag], elePhi[indexTag], ELE_MASS);
	  //*******************************************************
	  //Loop over Probe electrons
	  //*******************************************************
	  for(int indexProbe = 0; indexProbe < nElectrons; indexProbe++){
	    
	    if(elePt[indexProbe] < 5) continue;
	    if(fabs(eleEta[indexProbe]) > 2.5) continue;
	  
	    //skip the tag
	    if (indexTag == indexProbe) continue;
	  
	    //For MC, Match to Gen level electron
	    if (!isData) {
	      for (int q=0;q<int(genLeptonIndex.size()); q++) {
		if ( deltaR(eleEta[indexProbe],elePhi[indexProbe],
			    gParticleEta[genLeptonIndex[q]], 
			    gParticlePhi[genLeptonIndex[q]]) < 0.1) {
		  TPPair->matchGen = 1;
		}
	      }
	    }
	    
	    TLorentzVector vprobe;
	    vprobe.SetPtEtaPhiM(elePt[indexProbe], eleEta[indexProbe], elePhi[indexProbe], ELE_MASS);
	    TLorentzVector vdilep = vtag+vprobe;

	    if (vdilep.M() < 60 || vdilep.M() > 120) continue;

	    TPPair->lep1 = &vtag;
	    TPPair->lep2 = &vprobe;
	    TPPair->dilep = &vdilep;
	    TPPair->q1 = eleCharge[indexTag];
	    TPPair->q2 = eleCharge[indexProbe];

	    TPPair->d01 = ele_d0[indexTag];
	    TPPair->d02 = ele_d0[indexProbe];
	    TPPair->dz1 = ele_dZ[indexTag];
	    TPPair->dz2 = ele_dZ[indexProbe];
	    TPPair->sigieie1 = eleFull5x5SigmaIetaIeta[indexTag];
	    TPPair->sigieie2 = eleFull5x5SigmaIetaIeta[indexProbe];
	    TPPair->hovere1 = ele_HoverE[indexTag];
	    TPPair->hovere2 = ele_HoverE[indexProbe];
	    TPPair->eoverp1 = ele_OneOverEminusOneOverP[indexTag];
	    TPPair->eoverp2 = ele_OneOverEminusOneOverP[indexProbe];
	    TPPair->dphi1 = ele_dPhi[indexTag];
	    TPPair->dphi2 = ele_dPhi[indexProbe];
	    TPPair->deta1 = ele_dEta[indexTag];
	    TPPair->deta2 = ele_dEta[indexProbe];
	    TPPair->isConv1 = ( ele_PassConvVeto[indexTag] ? 0 : 1);
	    TPPair->isConv2 = ( ele_PassConvVeto[indexProbe] ? 0 : 1);
	    TPPair->nexphits1 = ele_MissHits[indexTag];
	    TPPair->nexphits2 = ele_MissHits[indexProbe];

	    TPPair->isTight2 = ( isEGammaPOGTightElectron(indexProbe,true,false,true) ? 1 : 0);
	    TPPair->isMedium2 = ( isEGammaPOGMediumElectron(indexProbe,true,false,true) ? 1 : 0);
	    TPPair->isLoose2 = ( isEGammaPOGLooseElectron(indexProbe,true,false,true) ? 1 : 0);
	    TPPair->isVeto2 = ( isMVANonTrigVetoElectron(indexProbe,true,false) ? 1 : 0);

	    TPPair->activity1 = ele_activityMiniIsoAnnulus[indexTag];
	    TPPair->pfChIso1 = ele_pileupIso[indexTag] + ele_chargedIso[indexTag];
	    TPPair->pfGamIso1 = ele_photonIso[indexTag];
	    TPPair->pfNeuIso1 = ele_neutralHadIso[indexTag];
	    TPPair->pfCombIso1 =  ele_pileupIso[indexTag] + ele_chargedIso[indexTag] + ele_photonIso[indexTag] + ele_neutralHadIso[indexTag];
	    TPPair->chMiniIso1 = ele_chargedMiniIso[indexTag];
	    TPPair->neuMiniIso1 = ele_photonAndNeutralHadronMiniIso[indexTag];
	    TPPair->puChMiniIso1 = ele_chargedPileupMiniIso[indexTag];
	    TPPair->activity2 = ele_activityMiniIsoAnnulus[indexProbe];
	    TPPair->pfChIso2 = ele_pileupIso[indexProbe] + ele_chargedIso[indexProbe];
	    TPPair->pfGamIso2 = ele_photonIso[indexProbe];
	    TPPair->pfNeuIso2 = ele_neutralHadIso[indexProbe];
	    TPPair->pfCombIso2 =  ele_pileupIso[indexProbe] + ele_chargedIso[indexProbe] + ele_photonIso[indexProbe] + ele_neutralHadIso[indexProbe];
	    TPPair->chMiniIso2 = ele_chargedMiniIso[indexProbe];
	    TPPair->neuMiniIso2 = ele_photonAndNeutralHadronMiniIso[indexProbe];
	    TPPair->puChMiniIso2 = ele_chargedPileupMiniIso[indexProbe];

	    TPPair->tree_->Fill();
	    
	  } //loop over probe electrons
	  
	} // loop over tag electrons
	
      } //if objects are electrons
      
      //*********************************************************
      //Muons
      //*********************************************************
      if (option == 2) {

	//*******************************************************
	//Loop over Tag muons
	//*******************************************************
	for(int indexTag = 0; indexTag < nMuons; indexTag++){

	  if(muonPt[indexTag] < 25) continue;
	  if(fabs(muonEta[indexTag]) > 2.4) continue;

	  //For MC, Match to Gen level electron
	  if (!isData) {
	    bool genmatch = false;
	    for (int q=0;q<int(genLeptonIndex.size()); q++) {
	      if ( deltaR(muonEta[indexTag],muonPhi[indexTag],
			  gParticleEta[genLeptonIndex[q]], 
			  gParticlePhi[genLeptonIndex[q]]) < 0.1) {
		genmatch = true;
	      }
	    }
	    if (!genmatch) continue;
	  }

	  //tag must pass tight cuts
	  if ( !isMuonPOGTightMuon(indexTag) ) continue;
	  
	  //Tag must match single electron HLT Filters OR tag leg filter of the dedicated T&P trigger
	  if ( !(muon_passHLTFilter[indexTag][11] && HLTDecision[16]) ) continue;

	  TLorentzVector vtag;
	  vtag.SetPtEtaPhiM(muonPt[indexTag], muonEta[indexTag], muonPhi[indexTag], MU_MASS);

	  //*******************************************************
	  //Loop over Probe muons
	  //*******************************************************
	  for(int indexProbe = 0; indexProbe < nMuons; indexProbe++){
	    
	    if(muonPt[indexProbe] < 5) continue;
	    if(fabs(muonEta[indexProbe]) > 2.4) continue;
	  
	    //skip the tag
	    if (indexTag == indexProbe) continue;
	  
	    //For MC, Match to Gen level electron
	    if (!isData) {
	      for (int q=0;q<int(genLeptonIndex.size()); q++) {
		if ( deltaR(muonEta[indexProbe],muonPhi[indexProbe],
			    gParticleEta[genLeptonIndex[q]], 
			    gParticlePhi[genLeptonIndex[q]]) < 0.1) {
		  TPPair->matchGen = true;
		}
	      }
	    }
	    
	    TLorentzVector vprobe;
	    vprobe.SetPtEtaPhiM(muonPt[indexProbe], muonEta[indexProbe], muonPhi[indexProbe], MU_MASS);

	    TLorentzVector vdilep = vtag+vprobe;

	    if (vdilep.M() < 60 || vdilep.M() > 120) continue;

            TPPair->lep1 = &vtag;
            TPPair->lep2 = &vprobe;
            TPPair->dilep = &vdilep;
            TPPair->q1 = muonCharge[indexTag];
            TPPair->q2 = muonCharge[indexProbe];

	    TPPair->isTight2 = ( muonIsLoose[indexProbe] ? 1:0);
            TPPair->isMedium2 = ( muonIsMedium[indexProbe] ? 1:0);
            TPPair->isLoose2 = ( muonIsLoose[indexProbe] ? 1:0);
            TPPair->isVeto2 = 0;

	    TPPair->muontype1 = muonType[indexTag];
	    TPPair->muontype2 = muonType[indexProbe];
	    TPPair->muonquality1 = muonQuality[indexTag];
	    TPPair->muonquality2 = muonQuality[indexProbe];
	    TPPair->muon_ip3d1 = muon_ip3d[indexTag];
	    TPPair->muon_ip3d2 = muon_ip3d[indexProbe];
	    TPPair->muon_ip3dSignificance1 = muon_ip3dSignificance[indexTag];
	    TPPair->muon_ip3dSignificance2 = muon_ip3dSignificance[indexProbe];

            TPPair->activity1 = muon_activityMiniIsoAnnulus[indexTag];
            TPPair->pfChIso1 = muon_pileupIso[indexTag] + muon_chargedIso[indexTag];
            TPPair->pfGamIso1 = muon_photonIso[indexTag];
            TPPair->pfNeuIso1 = muon_neutralHadIso[indexTag];
            TPPair->pfCombIso1 =  muon_pileupIso[indexTag] + muon_chargedIso[indexTag] + muon_photonIso[indexTag] + muon_neutralHadIso[indexTag];
            TPPair->chMiniIso1 = muon_chargedMiniIso[indexTag];
            TPPair->neuMiniIso1 = muon_photonAndNeutralHadronMiniIso[indexTag];
            TPPair->puChMiniIso1 = muon_chargedPileupMiniIso[indexTag];
            TPPair->activity2 = muon_activityMiniIsoAnnulus[indexProbe];
            TPPair->pfChIso2 = muon_pileupIso[indexProbe] + muon_chargedIso[indexProbe];
            TPPair->pfGamIso2 = muon_photonIso[indexProbe];
            TPPair->pfNeuIso2 = muon_neutralHadIso[indexProbe];
            TPPair->pfCombIso2 =  muon_pileupIso[indexProbe] + muon_chargedIso[indexProbe] + muon_photonIso[indexProbe] + muon_neutralHadIso[indexProbe];
            TPPair->chMiniIso2 = muon_chargedMiniIso[indexProbe];
            TPPair->neuMiniIso2 = muon_photonAndNeutralHadronMiniIso[indexProbe];
            TPPair->puChMiniIso2 = muon_chargedPileupMiniIso[indexProbe];

	    TPPair->tree_->Fill();
	    
	  } //loop over probe muons
	  
	} // loop over tag muons

      } //if objects are muons

    }//end of event loop


    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}

