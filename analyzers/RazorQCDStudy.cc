#include "RazorQCDStudy.h"
#include "RazorHelper.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "ControlSampleEvents.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

struct greater_than_pt{
    inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){
        return p1.Pt() > p2.Pt();
    }
};
 
void RazorQCDStudy::Analyze(bool isData, int option, string outputfilename, string label)
{

  //initialization: create one TTree for each analysis box 
  //cout << "Initializing..." << endl;
  bool isRunOne = (option == 1);
  //cout << "IsData = " << isData << "\n";
  bool isFastsimSMS = (option == 1 || option == 11);
  
  TRandom3 *random = new TRandom3(33333); //Artur wants this number 33333
  
  bool printSyncDebug = false;
  
  RazorHelper *helper = 0;
  string analysisTag = "Razor2016_80X";
  helper = new RazorHelper(analysisTag, isData, isFastsimSMS);
  
  // Get jet corrector
  std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
  std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();
  //JetCorrectorParameters *JetResolutionParameters = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt",pathname.c_str()));
  //SimpleJetResolution *JetResolutionCalculator = new SimpleJetResolution(*JetResolutionParameters);
  
  //*************************************************************************
  //Set up Output File
  //*************************************************************************
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "RazorControlRegions.root";
  TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
  TTree *outTree = new TTree("QCDTree", "Info on selected razor inclusive events");
  //tree variables
  Float_t                 weight;
  Float_t                 pileupWeight;
  UInt_t                  run;
  UInt_t                  lumi;
  UInt_t                  event;
  UInt_t                  NPU;
  UInt_t                  NPV;
  Float_t                 Rho;
  Bool_t                  hltDecision[150];
  Int_t                   hltPrescale[150];
  Float_t                 MR;
  Float_t                 Rsq;
  Float_t                 minDPhi;
  Float_t                 minDPhiN;
  Float_t                 dPhiRazor;
  Float_t                 MET;
  Float_t                 HT;
  Float_t                 MHT;
  Float_t                 MHTPhi;
  Float_t                 diffMetMht;
  UInt_t                  NJets40;
  UInt_t                  NJets80;
  UInt_t                  NBJetsLoose;
  UInt_t                  NBJetsMedium;
  UInt_t                  NBJetsTight;
  Float_t                 genJetMR;
  Float_t                 genJetRsq;
  Float_t                 genJetDPhiRazor;    
  Float_t                 genJetHT;
  Float_t                 maxJetGenDiff;
  int NJets;
  int box;
  int nLooseMuons, nTightMuons, nLooseElectrons, nTightElectrons, nTightTaus;
  int nVetoMuons, nVetoElectrons, nLooseTaus;
  float leadingJetPt, subLeadingJetPt;
  float leadingTightMuPt, leadingTightElePt;
  float leadingMuPt, leadingElePt;
  float zPt, zEta, zPhi, zMass;
  bool passedDileptonTrigger;
  bool passedSingleLeptonTrigger;
  bool passedHadronicTrigger;
  bool passedDijetTrigger;
  float JetE[99];
  float JetPt[99];
  float JetEta[99];
  float JetPhi[99];
  float GenJetE[99];
  float GenJetPt[99];
  float GenJetEta[99];
  float GenJetPhi[99];
  bool  JetIDTight[99];
  float PtNeutrinoClosestToJet[99]; 
  float DRNeutrinoClosestToJet[99];
  float JetChargedEMEnergyFraction[99];
  float JetNeutralEMEnergyFraction[99];
  float JetChargedHadEnergyFraction[99];
  float JetNeutralHadEnergyFraction[99];
  float JetPartonFlavor[99];
  float JetPileupID[99];
  float JetGenDiff[99];
  
  //book the branches that go in all types of trees
  outTree->Branch("weight",&weight,"weight/F");
  outTree->Branch("pileupWeight",&pileupWeight,"pileupWeight/F");
  outTree->Branch("run",&run,"run/i");
  outTree->Branch("lumi",&lumi,"lumi/i");
  outTree->Branch("event",&event,"event/i");
  outTree->Branch("NPU",&NPU,"NPU/i");
  //outTree->Branch("NPU_Minus1",&NPU_Minus1,"NPU_Minus1/i");
  //outTree->Branch("NPU_Plus1",&NPU_Plus1,"NPU_Plus1/i");
  outTree->Branch("NPV",&NPV,"NPV/i");
  outTree->Branch("Rho",&Rho,"Rho/F");
  outTree->Branch("HLTDecision",&hltDecision,"HLTDecision[150]/O");
  outTree->Branch("HLTPrescale",&hltPrescale,"HLTPrescale[150]/I");
  outTree->Branch("MR",&MR,"MR/F");
  outTree->Branch("Rsq",&Rsq,"Rsq/F");
  outTree->Branch("minDPhi",&minDPhi,"minDPhi/F"); 
  outTree->Branch("minDPhiN",&minDPhiN,"minDPhiN/F"); 
  outTree->Branch("dPhiRazor",&dPhiRazor,"dPhiRazor/F");
  outTree->Branch("MET",&MET,"MET/F");
  outTree->Branch("HT",&HT,"HT/F");
  outTree->Branch("MHT",&MHT,"MHT/F");
  outTree->Branch("MHTPhi",&MHTPhi,"MHTPhi/F");
  outTree->Branch("diffMetMht",&diffMetMht,"diffMetMht/F");
  outTree->Branch("NJets40",&NJets40,"NJets40/i");
  outTree->Branch("NJets80",&NJets80,"NJets80/i");
  outTree->Branch("NBJetsLoose",&NBJetsLoose,"NBJetsLoose/i");
  outTree->Branch("NBJetsMedium",&NBJetsMedium,"NBJetsMedium/i");
  outTree->Branch("NBJetsTight",&NBJetsTight,"NBJetsTight/i");
  outTree->Branch("leadingJetPt",&leadingJetPt,"leadingJetPt/F");
  outTree->Branch("subLeadingJetPt",&subLeadingJetPt,"subLeadingJetPt/F");
  outTree->Branch("leadingMuPt",&leadingMuPt,"leadingMuPt/F");
  outTree->Branch("leadingElePt",&leadingElePt,"leadingElePt/F");
  outTree->Branch("genJetMR",&genJetMR,"genJetMR/F");
  outTree->Branch("genJetRsq",&genJetRsq,"genJetRsq/F");
  outTree->Branch("genJetDPhiRazor",&genJetDPhiRazor,"genJetDPhiRazor/F");
  outTree->Branch("genJetHT",&genJetHT,"genJetHT/F");	  
  outTree->Branch("genMET",&genMetPt,"genMET/F");	         
  outTree->Branch("genMETPhi",&genMetPhi,"genMETPhi/F");	         
  outTree->Branch("maxJetGenDiff",&maxJetGenDiff,"maxJetGenDiff/F");
  outTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter,"Flag_HBHENoiseFilter/O");
  outTree->Branch("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, "Flag_HBHEIsoNoiseFilter/O");
  outTree->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter,"Flag_CSCTightHaloFilter/O");
  outTree->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter,"Flag_hcalLaserEventFilter/O");
  outTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter,"Flag_EcalDeadCellTriggerPrimitiveFilter/O");
  outTree->Branch("Flag_goodVertices", &Flag_goodVertices,"Flag_goodVertices/O");
  outTree->Branch("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter,"Flag_EcalDeadCellBoundaryEnergyFilter/O");
  outTree->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter,"Flag_trackingFailureFilter/O");
  outTree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter,"Flag_eeBadScFilter/O");
  outTree->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter,"Flag_ecalLaserCorrFilter/O");
  outTree->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters,"Flag_trkPOGFilters/O");
  outTree->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X,"Flag_trkPOG_manystripclus53X/O");
  outTree->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X,"Flag_trkPOG_toomanystripclus53X/O");
  outTree->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters,"Flag_trkPOG_logErrorTooManyClusters/O");
  outTree->Branch("Flag_METFilters", &Flag_METFilters,"Flag_METFilters/O");	
  outTree->Branch("passedDileptonTrigger", &passedDileptonTrigger,"passedDileptonTrigger/O");	
  outTree->Branch("passedSingleLeptonTrigger", &passedSingleLeptonTrigger,"passedSingleLeptonTrigger/O");	
  outTree->Branch("passedHadronicTrigger", &passedHadronicTrigger,"passedHadronicTrigger/O");	
  outTree->Branch("passedDijetTrigger", &passedDijetTrigger,"passedDijetTrigger/O");	
  outTree->Branch("nJets", &NJets,"nJets/I");
  outTree->Branch("box", &box, "box/I");
  outTree->Branch("nVetoMuons", &nVetoMuons, "nVetoMuons/I");
  outTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
  outTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
  outTree->Branch("nVetoElectrons", &nVetoElectrons, "nVetoElectrons/I");
  outTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
  outTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
  outTree->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
  outTree->Branch("zPt",&zPt,"zPt/F");
  outTree->Branch("zEta",&zEta,"zEta/F");
  outTree->Branch("zPhi",&zPhi,"zPhi/F");
  outTree->Branch("zMass",&zMass,"zMass/F");
  outTree->Branch("JetE", JetE,"JetE[nJets]/F");
  outTree->Branch("JetPt", JetPt,"JetPt[nJets]/F");
  outTree->Branch("JetEta", JetEta,"JetEta[nJets]/F");
  outTree->Branch("JetPhi", JetPhi,"JetPhi[nJets]/F");
  outTree->Branch("GenJetE", GenJetE,"GenJetE[nJets]/F");
  outTree->Branch("GenJetPt", GenJetPt,"GenJetPt[nJets]/F");
  outTree->Branch("GenJetEta", GenJetEta,"GenJetEta[nJets]/F");
  outTree->Branch("GenJetPhi", GenJetPhi,"GenJetPhi[nJets]/F");
  outTree->Branch("JetGenDiff", JetGenDiff,"JetGenDiff[nJets]/F");
  outTree->Branch("JetIDTight", JetIDTight,"JetIDTight[nJets]/O");
  outTree->Branch("PtNeutrinoClosestToJet", PtNeutrinoClosestToJet,"PtNeutrinoClosestToJet[nJets]/F");
  outTree->Branch("DRNeutrinoClosestToJet", DRNeutrinoClosestToJet,"DRNeutrinoClosestToJet[nJets]/F");
  outTree->Branch("JetChargedEMEnergyFraction", JetChargedEMEnergyFraction,"JetChargedEMEnergyFraction[nJets]/F");
  outTree->Branch("JetNeutralEMEnergyFraction", JetNeutralEMEnergyFraction,"JetNeutralEMEnergyFraction[nJets]/F");
  outTree->Branch("JetChargedHadEnergyFraction", JetChargedHadEnergyFraction,"JetChargedHadEnergyFraction[nJets]/F");
  outTree->Branch("JetNeutralHadEnergyFraction", JetNeutralHadEnergyFraction,"JetNeutralHadEnergyFraction[nJets]/F");
  outTree->Branch("JetPartonFlavor", JetPartonFlavor,"JetPartonFlavor[nJets]/F");
  outTree->Branch("JetPileupID", JetPileupID,"JetPileupID[nJets]/F");
  
  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 0.5, 1.5);
  TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
  
  //*************************************************************************
  //Look over Input File Events
  //*************************************************************************
  if (fChain == 0) return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  //for (Long64_t jentry=0; jentry<100;jentry++) {
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
      //begin event
      if(jentry % 1000000 == 0) cout << "Processing entry " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      printSyncDebug = false;

      //Fill normalization histogram
      NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
      SumWeights->Fill(1.0, weight);

      box=0;
      nVetoMuons = 0;
      nLooseMuons = 0;
      nTightMuons = 0;
      nVetoElectrons = 0;
      nLooseElectrons = 0;
      nTightElectrons = 0;
      nLooseTaus = 0;
      nTightTaus = 0;
      //mT = -1;
      //mTLoose = -1;
      leadingJetPt = -1;
      subLeadingJetPt = -1;
      leadingTightMuPt = -1;
      leadingTightElePt = -1;
      leadingMuPt = -1;
      leadingElePt = -1;
      zPt=-1;
      zEta=-1; 
      zPhi=-1;
      zMass=-1;

      //event info
      weight = genWeight;
      run = runNum;
      lumi = lumiNum;
      event = eventNum;

      pileupWeight = 1.0;
      if(!isData){
	//Get number of PU interactions
	for (int i = 0; i < nBunchXing; i++) {
	  if (BunchXing[i] == 0) {
	    NPU= nPUmean[i];
	  }
	}
	pileupWeight = helper->getPileupWeight(NPU);
	//pileupWeightUp = helper->getPileupWeightUp(NPU) / pileupWeight;
	//pileupWeightDown = helper->getPileupWeightDown(NPU) / pileupWeight;
      }
      
      NPV = nPV;

      passedDileptonTrigger = false;
      passedSingleLeptonTrigger = false;
      passedHadronicTrigger= false;

      vector<int> dileptonTriggerNums = helper->getDileptonTriggerNums();
      vector<int> singleLeptonTriggerNums = helper->getSingleLeptonTriggerNums();
      vector<int> hadronicTriggerNums = helper->getHadronicTriggerNums();
      for( unsigned int itrig = 0; itrig < dileptonTriggerNums.size(); itrig++ ) {
	if (HLTDecision[dileptonTriggerNums[itrig]]) {
	  passedDileptonTrigger = true;
	  break;
	}
      }
      for( unsigned int itrig = 0; itrig < singleLeptonTriggerNums.size(); itrig++ ) {
	if (HLTDecision[singleLeptonTriggerNums[itrig]]) {
	  passedSingleLeptonTrigger = true;
	  break;
	}
      }
      for( unsigned int itrig = 0; itrig < hadronicTriggerNums.size(); itrig++ ) {
	if (HLTDecision[hadronicTriggerNums[itrig]]) {
	  passedHadronicTrigger = true;
	  break;
	}
      }
      if (HLTDecision[129]) passedDijetTrigger = true ;

      //ignore trigger for Fastsim, and for 80X MC
      if(analysisTag == "Razor2016_80X" && !isData){
	passedDileptonTrigger = true;
	passedSingleLeptonTrigger = true;
	passedHadronicTrigger = true;
      }

      vector<TLorentzVector> GoodLeptons; //leptons used to compute hemispheres
      TLorentzVector leadingTightMu, leadingTightEle; //used for mT calculation
      TLorentzVector mu1, mu2; //used for Z control region
      TLorentzVector ele1, ele2; //used for Z control region
      mu1.SetPtEtaPhiM(0,0,0,0); mu2.SetPtEtaPhiM(0,0,0,0);
      ele1.SetPtEtaPhiM(0,0,0,0); ele2.SetPtEtaPhiM(0,0,0,0);

      for(int i=0; i<nMuons; i++) {
	if (muonPt[i] < 20 || abs(muonEta[i])>2.4) continue;
        TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
	if(!isVetoMuon(i)) continue;

	if (muonPt[i]>mu1.Pt()) {
          mu2=mu1;
          mu1=thisMuon;
          leadingMuPt=muonPt[i];
        }
        else if (muonPt[i]>mu2.Pt()) mu2=thisMuon;
      }

      for(int i = 0; i < nMuons; i++){

        if(muonPt[i] < 5) continue;
        if(abs(muonEta[i]) > 2.4) continue;

        //Calculate MC->Data Scale Factors
        if (RazorAnalyzer::matchesGenMuon(muonEta[i],muonPhi[i])) {
	  //apply muon efficiency scale factors
        }

        TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 

        if(isVetoMuon(i)) nVetoMuons++;
        if(isLooseMuon(i) && muonPt[i] >= 20 ) nLooseMuons++;
        if(isTightMuon(i) && muonPt[i] >= 20){
	  nTightMuons++;
	  if(muonPt[i] > leadingTightMuPt){
	    leadingTightMu = thisMuon;
	    leadingTightMuPt = muonPt[i];
	  }
        }

        if(!isVetoMuon(i)) continue;  
	if (thisMuon==mu1||thisMuon==mu2) continue;
        GoodLeptons.push_back(thisMuon); 
      }

      for(int i=0; i<nElectrons; i++) {
	if(elePt[i]<25 || abs(eleEta[i])>2.5) continue;
	TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
	if(!isVetoElectron(i)) continue;
	if (elePt[i]>ele1.Pt()) {
          ele2=ele1;
          ele1=thisElectron;
          leadingElePt=elePt[i];
        }
	else if (elePt[i]>ele2.Pt()) ele2=thisElectron;
      }

      for(int i = 0; i < nElectrons; i++){
        if(elePt[i] < 5) continue;
        if(fabs(eleEta[i]) > 2.5) continue;

        //if (RazorAnalyzer::matchesGenElectron(eleEta[i],elePhi[i])) {
	  //No Efficiency Scale Factors for Run2 Yet
	//}

        TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);

        if(isVetoElectron(i)) nVetoElectrons++;
        if( isLooseElectron(i) && elePt[i] > 25 ) nLooseElectrons++;
        if( isTightElectron(i) && elePt[i] > 25 ){
	  nTightElectrons++;
	  if(elePt[i] > leadingTightElePt){
	    leadingTightEle = thisElectron;
	    leadingTightElePt = elePt[i];
	  }
        }

        //remove overlaps
        bool overlap = false;
        for(auto& lep : GoodLeptons){
	  if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
        }
        if(overlap) continue;

        if(!isVetoElectron(i)) continue; 
	if(thisElectron==ele1||thisElectron==ele2) continue;
        GoodLeptons.push_back(thisElectron);            
      }

      TLorentzVector z; z.SetPtEtaPhiM(0,100,0,0);
      if (mu1.Pt()>20 && mu2.Pt()>20) z=mu1+mu2;
      else if (ele1.Pt()>25 && ele2.Pt()>25) z=ele1+ele2;

      zPt=z.Pt();
      zEta=z.Eta();
      zPhi=z.Phi();
      zMass=z.M();

      if(zMass<60 || zMass>120) {
	z.SetPtEtaPhiM(0,100,0,0);
	if (mu1.Pt()>0) GoodLeptons.push_back(mu1);
	if (mu2.Pt()>0) GoodLeptons.push_back(mu2);
	if (ele1.Pt()>0) GoodLeptons.push_back(ele1);
	if (ele2.Pt()>0) GoodLeptons.push_back(ele2);
      }

      //******************************
      //Only Do Taus for Run2
      //******************************
      for(int i = 0; i < nTaus; i++){ 
        if (tauPt[i] < 20) continue;
        if (fabs(tauEta[i]) > 2.4) continue;

        if(isLooseTau(i)) nLooseTaus++;
        if(isTightTau(i)) nTightTaus++;

        //remove overlaps
        bool overlap = false;
        for(auto& lep : GoodLeptons){
	  if (RazorAnalyzer::deltaR(tauEta[i],tauPhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
        }
        if(overlap) continue;

        if (!isLooseTau(i)) continue;
        TLorentzVector thisTau; thisTau.SetPtEtaPhiM(tauPt[i], tauEta[i], tauPhi[i], 1.777);
        GoodLeptons.push_back(thisTau);  
      }

      //************************************************************************
      //Find all Relevant Jets	
      //************************************************************************
      vector<TLorentzVector> GoodJets;
      int numJetsAbove80GeV = 0;
      int numJetsAbove40GeV = 0;
      double MetX_Type1Corr = 0;
      double MetY_Type1Corr = 0;
      int nBJetsLoose20GeV = 0;
      int nBJetsMedium20GeV = 0;
      int nBJetsTight20GeV = 0;
      minDPhi = 9999;
      minDPhiN = 9999;
      maxJetGenDiff=0;
      if (printSyncDebug) cout << "NJets: " << nJets << "\n";
      NJets = nJets;
      for(int i = 0; i < nJets; i++){

	//*******************************************************
	//apply jet iD
	//*******************************************************
	if (!jetPassIDTight[i]) continue;

	//*****************************************************************
        //exclude selected muons and electrons from the jet collection
        //*****************************************************************
        double deltaR = -1;
        for(auto& lep : GoodLeptons){
	  double thisDR = RazorAnalyzer::deltaR(jetEta[i],jetPhi[i],lep.Eta(),lep.Phi());  
	  if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
        }
        if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton


	//*******************************************************
	//Correct Jet Energy Scale and Resolution
	//*******************************************************
	double tmpRho = fixedGridRhoFastjetAll;
	if (isRunOne) tmpRho = fixedGridRhoAll;
	double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
					       tmpRho, jetJetArea[i], runNum, JetCorrectorIOV,JetCorrector);

	//double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
	//tmpRho, jetJetArea[i], runNum,
	//JetCorrectorIOV,JetCorrector);   
	
	double JECLevel1 = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
						     tmpRho, jetJetArea[i], runNum, JetCorrectorIOV,JetCorrector, 0);  
	
	//double JECLevel1 = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
	//tmpRho, jetJetArea[i], runNum,
	//JetCorrectorIOV,JetCorrector, 0); 
	
	Rho = tmpRho;
	
	double jetEnergySmearFactor = 1.0;
	if (!isData) {
	  std::vector<float> fJetEta, fJetPtNPU;
	  fJetEta.push_back(jetEta[i]);  
	  fJetPtNPU.push_back(jetPt[i]*JEC); 
	  fJetPtNPU.push_back(NPU); 
	  if (printSyncDebug) {
	    cout << "Jet: " << jetPt[i] << " " << jetEta[i] << " " << jetPhi[i] << "\n";
	    cout << "Jet Resolution : " << jetPt[i]*JEC << " " << jetEta[i] << " " << jetPhi[i] << " \n ";
	    //<< JetResolutionCalculator->resolution(fJetEta,fJetPtNPU) << "\n";
	  }
	  //jetEnergySmearFactor = JetEnergySmearingFactor( jetPt[i]*JEC, jetEta[i], NPU_0, JetResolutionCalculator, random);
	  jetEnergySmearFactor = 1.0;
	}
	if (printSyncDebug) {
	  cout << "Jet Smearing Factor " << jetEnergySmearFactor << "\n";
	}
	
	TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);
	TLorentzVector L1CorrJet = makeTLorentzVector(jetPt[i]*JECLevel1, jetEta[i], jetPhi[i], jetE[i]*JECLevel1);
	TLorentzVector UnCorrJet = makeTLorentzVector(jetPt[i], jetEta[i], jetPhi[i], jetE[i]);
	double jetCorrPt = jetPt[i]*JEC*jetEnergySmearFactor;
	
	//*******************************
	//Add to Type1 Met Correction
	//*******************************
	if (jetCorrPt > 15 && 
	    jetChargedEMEnergyFraction[i] + jetNeutralEMEnergyFraction[i] <= 0.9
	    ) {
	  MetX_Type1Corr += -1 * ( thisJet.Px() - L1CorrJet.Px()  );
	  MetY_Type1Corr += -1 * ( thisJet.Py() - L1CorrJet.Py()  );
	  if (printSyncDebug) cout << "Met Type1 Corr: " << thisJet.Px() - L1CorrJet.Px() << " " << thisJet.Py() - L1CorrJet.Py() << "\n";
	}

	//*******************************************************
	//Jet Variables
	//*******************************************************
	JetE[i] = jetE[i]*JEC*jetEnergySmearFactor;
	JetPt[i] = jetPt[i]*JEC*jetEnergySmearFactor;
	JetEta[i] = jetEta[i];
	JetPhi[i] = jetPhi[i];
	JetIDTight[i] = jetPassIDTight[i];
	JetChargedEMEnergyFraction[i] = jetChargedEMEnergyFraction[i];
	JetNeutralEMEnergyFraction[i] = jetNeutralEMEnergyFraction[i];
	JetChargedHadEnergyFraction[i] = jetChargedHadronEnergyFraction[i];
	JetNeutralHadEnergyFraction[i] = jetNeutralHadronEnergyFraction[i];
	JetPartonFlavor[i] = jetPartonFlavor[i];
	JetPileupID[i] = jetPileupId[i];

	if (JetPt[i]>leadingJetPt && JetPt[i]<14000) {
	  subLeadingJetPt=leadingJetPt;
	  leadingJetPt=JetPt[i];
	}
	else if (JetPt[i]>subLeadingJetPt && JetPt[i]<14000) {
	  subLeadingJetPt=JetPt[i];
	}

	//Match To GenJet
	GenJetE[i] = -999;
	GenJetPt[i] = -999;
	GenJetEta[i] = -999;
	GenJetPhi[i] = -999;
	double minDRToGenJet = 9999;	
	for(int j = 0; j < nGenJets; j++){
	  double DR = RazorAnalyzer::deltaR( genJetEta[j], genJetPhi[j], jetEta[i], jetPhi[i]);
	  if (DR > 0.4) continue;
	
	  if (DR < minDRToGenJet) {
	    minDRToGenJet = DR;
	    GenJetE[i] = genJetE[j];
	    GenJetPt[i] = genJetPt[j];
	    GenJetEta[i] = genJetEta[j];
	    GenJetPhi[i] = genJetPhi[j];
	    JetGenDiff[i] = JetPt[i]-GenJetPt[i];
	    if (fabs(JetGenDiff[i])>fabs(maxJetGenDiff)) maxJetGenDiff=JetGenDiff[i];
	  }     
	}
	//./include/RazorAnalyzer.h:double deltaR(double eta1, double phi1, double eta2, double phi2);
	
	//Check for neutrinos near jet
	PtNeutrinoClosestToJet[i] = -999;
	double minDRNeutrino = 9999;
	for(int j = 0; j < nGenParticle; j++){
	  if (!(abs(gParticleId[j]) == 12 || abs(gParticleId[j]) == 14 || abs(gParticleId[j]) == 16)) continue;
	  double DR = RazorAnalyzer::deltaR( jetEta[i], jetPhi[i], gParticleEta[j], gParticlePhi[j] );	  
	  if ( DR < minDRNeutrino ) {
	    minDRNeutrino = DR;
	    PtNeutrinoClosestToJet[i] = gParticlePt[j];
	  }
	}
	DRNeutrinoClosestToJet[i] = minDRNeutrino;

	if(jetCorrPt < 40) continue;
	if(fabs(jetEta[i]) > 3.0) continue;

	if (jetPt[i]*JEC > 20 && isCSVL(i)) nBJetsLoose20GeV++;
	if (jetPt[i]*JEC > 20 && isCSVM(i)) nBJetsMedium20GeV++;
	if (jetPt[i]*JEC > 20 && isCSVT(i)) nBJetsTight20GeV++;
	
	numJetsAbove40GeV++;

	if(jetCorrPt > 80) {
	  //cout << jentry << ": " << jetCorrPt << ", " << jetEta[i] << endl;
	  numJetsAbove80GeV++;	  
	}

	GoodJets.push_back(thisJet);	  

      } //loop over jets

      ///sort good jets
      TLorentzVector tmpjet;
      for (int i=0;i<int(GoodJets.size());i++) {
	for (int j=0;j<int(GoodJets.size()-1);j++) {
	  if (GoodJets[j+1].Pt() > GoodJets[j].Pt()) {
	    tmpjet = GoodJets[j]; //swap elements
	    GoodJets[j] = GoodJets[j+1];
	    GoodJets[j+1] = tmpjet;
	  }
	}
      }
  
      //compute minDPhi variables
      double JER = 0.1; //average jet energy resolution
      for (int i=0;i<int(GoodJets.size());i++) {
	if (i>2) continue; //only consider first 3 jets for minDPhi 	  
	double dPhi = fmin( fabs(deltaPhi(metPhi,GoodJets[i].Phi())) , fabs( 3.1415926 - fabs(deltaPhi(metPhi,GoodJets[i].Phi()))));
	if (dPhi < minDPhi) minDPhi = dPhi;
	double deltaT = 0;
	for(auto& jet2 : GoodJets) {
	  deltaT += pow(JER*jet2.Pt()*sin(fabs(deltaPhi(GoodJets[i].Phi(),jet2.Phi()))),2);
	}
	double dPhiN = dPhi / atan2( sqrt(deltaT) , metPt);
	if (dPhiN < minDPhiN) minDPhiN = dPhiN;
      }

      //Compute the razor variables using the selected jets and possibly leptons
      vector<TLorentzVector> GoodPFObjects;
      for(auto& jet : GoodJets) GoodPFObjects.push_back(jet);
      for(auto& lep : GoodLeptons) GoodPFObjects.push_back(lep);
      if (z.Eta()!=100) GoodPFObjects.push_back(z);

      double PFMetCustomType1CorrectedX = metPt*cos(metPhi) + MetX_Type1Corr;
      double PFMetCustomType1CorrectedY = metPt*sin(metPhi) + MetY_Type1Corr;
      TLorentzVector PFMETCustomType1Corrected;
      PFMETCustomType1Corrected.SetPxPyPzE(PFMetCustomType1CorrectedX, PFMetCustomType1CorrectedY, 0,
					   sqrt( pow(PFMetCustomType1CorrectedX,2) + pow(PFMetCustomType1CorrectedY,2)));
      TLorentzVector PFMETUnCorr = makeTLorentzVectorPtEtaPhiM(metPt, 0, metPhi, 0);
      TLorentzVector PFMETType1 = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);
      TLorentzVector PFMETType0Plus1 = makeTLorentzVectorPtEtaPhiM(metType0Plus1Pt, 0, metType0Plus1Phi, 0);
      TLorentzVector PFMETNoHF = makeTLorentzVectorPtEtaPhiM(metNoHFPt, 0, metNoHFPhi, 0);
      TLorentzVector MyMET = PFMETCustomType1Corrected; //This is the MET that will be used below.

      TLorentzVector PFMET = MyMET; //PFMET.SetPxPyPzE(PFMetX, PFMetY, 0, sqrt(PFMetX*PFMetX + PFMetY*PFMetY));

      if (printSyncDebug) {
	cout << "UnCorrectedMET: " << PFMETUnCorr.Pt() << " " << PFMETUnCorr.Phi() << "\n";
	cout << "Corrected PFMET: " << PFMET.Pt() << " " << PFMET.Phi() << " | X,Y Correction :  " << MetX_Type1Corr << " " << MetY_Type1Corr << "\n";
      }

      //*************************************************************************
      //Make hemispheres
      //*************************************************************************	
      MR = 0;
      Rsq = 0;

      float MhtX=0, MhtY=0;
      HT = 0;
      for(auto& obj : GoodPFObjects) { HT += obj.Pt(); MhtX += obj.Px(); MhtY += obj.Py(); }

      vector<TLorentzVector> hemispheres = getHemispheres(GoodPFObjects);
      MR = computeMR(hemispheres[0], hemispheres[1]); 
      Rsq = computeRsq(hemispheres[0], hemispheres[1], MyMET);
      dPhiRazor = deltaPhi(hemispheres[0].Phi(),hemispheres[1].Phi());
      MET = MyMET.Pt();

      TLorentzVector MyMHT;
      MyMHT.SetPxPyPzE(- MhtX, -MhtY, 0, sqrt( pow(MhtX,2) + pow(MhtY,2)));

      MHT = MyMHT.Pt();
      MHTPhi = MyMHT.Phi();

      TLorentzVector diff = MyMHT - MyMET;
      diffMetMht = abs(diff.Pt())/MyMET.Pt();

      NJets40 = numJetsAbove40GeV;
      NJets80 = numJetsAbove80GeV;
      NBJetsLoose = nBJetsLoose20GeV;
      NBJetsMedium = nBJetsMedium20GeV;
      NBJetsTight = nBJetsTight20GeV;
	
      ////*************************************************************************
      ////Make GenJet vector for hemispheres
      ////*************************************************************************
      double GenMetX = genMetPt*cos(genMetPhi);
      double GenMetY = genMetPt*sin(genMetPhi);
      TLorentzVector GenMET; GenMET.SetPxPyPzE(GenMetX, GenMetY, 0, sqrt(GenMetX*GenMetX + GenMetY*GenMetY));
      
      genJetMR = 0;
      genJetRsq = 0;
      genJetDPhiRazor = 0;
      genJetHT = 0;
      vector<TLorentzVector> GenJetObjects;
      for(int j = 0; j < nGenJets; j++){
      	if (genJetPt[j] > 40 && fabs(genJetEta[j]) < 3) {
      	  TLorentzVector thisGenJet = makeTLorentzVector(genJetPt[j], genJetEta[j], genJetPhi[j], genJetE[j]);
      	  GenJetObjects.push_back(thisGenJet);
      	  genJetHT += genJetPt[j];
      	}
      }
      if (GenJetObjects.size() >= 2 ) {
      	vector<TLorentzVector> tmpHemispheres = getHemispheres(GenJetObjects);
      	genJetMR = computeMR(tmpHemispheres[0], tmpHemispheres[1]); 
      	genJetRsq = computeRsq(tmpHemispheres[0], tmpHemispheres[1], GenMET);
      	genJetDPhiRazor = deltaPhi(tmpHemispheres[0].Phi(),tmpHemispheres[1].Phi());
      }      
     
      //*************************************************************************
      //save HLT Decisions
      //*************************************************************************
      for(int k=0; k<150; ++k) {
      	hltDecision[k] = HLTDecision[k];
	hltPrescale[k] = HLTPrescale[k];
      }

      //if (HLTDecision[106]) cout << jentry << ", " << hltDecision[106] << endl;
      
      //*************************************************************************
      //Skimming
      //*************************************************************************
      //bool passSkim = true;
      //if (option == 1) {
      //if (!(MR>1000 && Rsq>0.25)) {
      //passSkim = false;
      //}
      //}

      //*************************************************************************
      //Fill Tree
      //*************************************************************************
      //bool combineTrees=true;
      
      //MuMu Box
      if((passedDileptonTrigger||passedSingleLeptonTrigger) && nTightMuons > 0 && nLooseMuons > 1){// && NJets40 > 0 ){
	box = MuMu;
	if (zMass<60 || zMass>120) continue;
      }
      //EleEle Box
      else if((passedDileptonTrigger||passedSingleLeptonTrigger) && nTightElectrons > 0 && nLooseElectrons > 1){// && NJets40 > 0 ){
	box = EleEle;
	if (zMass<60 || zMass>120) continue;
      }
      else if (nVetoMuons+nVetoElectrons>0) {
      	continue;
      }
      //MultiJet Box                                
      else if(passedHadronicTrigger && NJets80 >= 2 && NJets40 > 3){
      	box = FourJet;
      }
      else if (NJets80>1 && nVetoElectrons==0 && nVetoMuons==0) {
      	box = DiJet;
      }
      if (!(box==MuMu||box==EleEle||box==FourJet||box==DiJet)) continue;
      //if (!(box==MuMu||box==EleEle)) continue;
      outTree->Fill();
      
    }
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}

