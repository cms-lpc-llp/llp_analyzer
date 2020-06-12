#include "SusyLLP.h"
#include "RazorHelper.h"
#include "SusyLLPTree.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"

//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"

#define _debug 0
#define _debug_npu 0
#define _debug_sync 0
#define _debug_met 0
#define _debug_jet 0
#define _debug_trk 0
#define _debug_lep 0
#define _debug_match 0
#define _debug_avgH 0
#define _debug_trg 0

#define N_MAX_LLP_DAUGHTERS 4
#define N_MAX_LLP_GRAND_DAUGHTERS 4
#define N_MAX_LLPS 2
#define N_MAX_LEPTONS 20
#define N_MAX_JETS 20
#define NTriggersMAX 602 //Number of trigger in the .dat file
using namespace std;

struct leptons
{
	TLorentzVector lepton;
	int pdgId;
	// float dZ;
	// bool passLooseId;
	// bool passMediumId;
	// bool passTightId;
	// bool passId;
};

struct taus
{
	TLorentzVector tau;
};

struct photons
{
	TLorentzVector photon;
};


struct jets
{
	TLorentzVector jet;
	//pat::Jet jet;
	float time;
	int jetNeutralHadronMultiplicity;
	int jetChargedHadronMultiplicity;
	int jetMuonMultiplicity;
	int jetElectronMultiplicity;
	int jetPhotonMultiplicity;
	float jetNeutralHadronEnergyFraction;
	float jetChargedHadronEnergyFraction;
	float jetMuonEnergyFraction;
	float jetElectronEnergyFraction;
	float jetPhotonEnergyFraction;
	float jetCSV;
	int ecalNRechits;
	float ecalRechitE;

	float jetGammaMax_ET;
	float jetMinDeltaRPVTracks;

	float jetChargedEMEnergyFraction;
	float jetNeutralEMEnergyFraction;

	int jetChargedMultiplicity;
	//int jetNHits;
	//int jetNPixelHits;
	float jetNPixelHitsMedian;
	float jetNHitsMedian;

	//int   jetNSV;
	//int   jetNSVCand;
	int   jetNVertexTracks;
	int   jetNSelectedTracks;
	float jetDRSVJet;
	//float jetFlightDist2D;
	//float jetFlightDist2DError;
	//float jetFlightDist3D;
	//float jetFlightDist3DError;
	//float jetSV_x;
	//float jetSV_y;
	//float jetSV_z;
	//int   jetSVNTracks;
	float jetSVMass;
};

//pt comparison
//not used so far
struct greater_than_pt
{
	inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
};

//lepton highest pt comparator
struct largest_pt_lep
{
	inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt_lep;

//jet highest pt comparator
struct largest_pt_jet
{
	inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;

//Analyze
void SusyLLP::Analyze(bool isData, int options, string outputfilename, string analysisTag, string process)
{
	//initialization: create one TTree for each analysis box
	cout << "Initializing..." << endl;
	cout << "IsData = " << isData << "\n";
	cout << "options = " << options << "\n";

	//---------------------------
	//-----------option----------
	//---------------------------
	int option;
	std::string label;
	bool pf;

	//HUNDRED'S DIGIT
	//option of run condor or locally
	if (options < 200){
		option = 1; // used when running condor
	}
	else{
		option = 0;// used when running locally
		//options need to be larger than 200, if do local test run
		cout << "option = 0, running locally, load aux locally \n";
	}

	//TEN'S DIGIT
	// label of signal / bkg
	if ((options/10)%10 == 1){
		label = "wH";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 2){
		label = "zH";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 3){
		label = "bkg_wH";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 4){
		label = "bkg_zH";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 5){
		label = "HH";
		cout << "process label: " << label << "\n";
	}
	else if ((options/10) % 10 == 6){
		label = "bkg_HH";
		cout << "process label: " << label << "\n";
	}
	else{
		cout << "What process it is? Label not defined. \n";
	}

	//UNIT'S DIGIT
	// pf option
	if(options%10==1){
		pf = true;
	}
	else{
		pf = false;
	}

	// DATA or MC
	if( isData )
	{
		std::cout << "[INFO]: running on data with label: " << label << " and option: " << option << " and pfjet is " << pf << std::endl;
	}
	else
	{
		std::cout << "[INFO]: running on MC with label: " << label << " and option: " << option << " and pfjet is " << pf << std::endl;
	}

	const float ELE_MASS = 0.000511;
	const float MU_MASS  = 0.105658;
	const float TAU_MASS  = 1.77686;
	const float Z_MASS   = 91.2;
	const float H_MASS   = 125.0;

	string dataset = "94X";
	//Analysis Tag
	//reference in RazorHelper
	if (analysisTag == ""){
		analysisTag = "Razor2016_80X";
	}

	int wzId;
	int NTrigger;//Number of trigger in trigger paths
	int elePt_cut = 0;
	int muonPt_cut = 0;
	uint nLepton_cut = 0;



	if (label == "HH" || label == "bkg_HH" ){
		NTrigger = 1;
		muonPt_cut = 15;
		elePt_cut = 15;
		nLepton_cut = 0;
	}
	if (label == "zH" || label == "bkg_zH" ){
		NTrigger = 4;
		muonPt_cut = 15;
		elePt_cut = 15;
		nLepton_cut = 2;
	}
	//else{}
	if (label == "wH" || label == "bkg_wH" ){
		NTrigger = 2;
		muonPt_cut = 27;
		elePt_cut = 32;
		nLepton_cut = 1;
	}

	int trigger_paths[NTrigger];
	if (label == "wH" || label == "bkg_wH"){
		wzId = 24;
		trigger_paths[0] = 87;
		trigger_paths[1] = 135;
		// trigger_paths[2] = 310;
	}
	else if (label == "zH" || label == "bkg_zH"){
		wzId = 23;
		trigger_paths[0] = 177;
		trigger_paths[1] = 362;
		// trigger_paths[2] = 310;
		trigger_paths[2] = 87;
		trigger_paths[3] = 135;
	}
	else if (label == "HH" || label == "bkg_HH"){
		wzId = 25;
		trigger_paths[0] = 310;
	}

	//--------------------------------
	//Initialize helper
	//--------------------------------
	RazorHelper *helper = 0;
	if (process == ""){
		process = "ZJetsToNuNu_HT-100ToInf_13TeV-madgraph_Summer16_2016";
	}

	helper = new RazorHelper(analysisTag, isData, false, process.c_str());

	//-----------------------------------------------
	//Set up Output File
	//-----------------------------------------------
	string outfilename = outputfilename;
	if (outfilename == "") outfilename = "SusyLLPTree.root";
	TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
	//RazorLiteTree *llp_tree = new RazorLiteTree;
	SusyLLPTree *llp_tree = new SusyLLPTree;
	llp_tree->CreateTree();
	llp_tree->tree_->SetAutoFlush(0);
	llp_tree->InitTree();
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
		if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		if(_debug) std::cout << "passed break " << jentry << std::endl;
		//GetEntry(ientry);
		nb = fChain->GetEntry(jentry); 
		if(_debug) std::cout << "passed GetEntry " << jentry << std::endl;
		nbytes += nb;
		if(_debug) std::cout << "passed nbytes " << jentry << std::endl;

		//fill normalization histogram
		if(_debug) std::cout << "deb0 " << jentry << std::endl;
		llp_tree->InitVariables();
		if(_debug) std::cout << "deb1 " << jentry << std::endl;
		if (label =="bkg_wH"|| label == "bkg_zH" || label == "bkg_HH"){
			if (isData)
			{
				NEvents->Fill(1);
				llp_tree->weight = 1;
			}
			else
			{
				//NEvents->Fill(genWeight);
				//llp_tree->weight = genWeight;
				NEvents->Fill(1);
				llp_tree->weight = 1;
			}

		}
		else{
			//generatedEvents->Fill(1);
			llp_tree->weight = 1;
		}
		if(_debug) std::cout << "deb2 " << jentry << std::endl;
		//event info
		llp_tree->runNum = runNum;
		llp_tree->lumiSec = lumiNum;
		llp_tree->evtNum = eventNum;
		if(_debug) std::cout << "deb3 " << jentry << std::endl;
		if(_debug) std::cout << "nBunchXing " << nBunchXing << std::endl;
		if (label == "zH" || label == "wH" ){
			NEvents->Fill(1);
			bool wzFlag = false;
			for (int i=0; i < nGenParticle; ++i)
			{
				// if (abs(gParticleId[i]) == wzId && gParticleStatus[i] == 22)
				//if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == wzId)
				if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] && abs(gParticleMotherId[i]) == wzId)
				{
					wzFlag = true;
				}

			}
			if ( wzFlag == false ) continue;

		}
		else if (label == "HH"){
			NEvents->Fill(1);
			bool wzFlag = false;
			for (int i=0; i < nGenParticle; ++i)
			{
				if ((abs(gParticleId[i]) == 5) && gParticleStatus[i] && abs(gParticleMotherId[i]) == wzId)
				{
					wzFlag = true;
				}

			}
			//if(_debug) std::cout << "wzFlag " << wzFlag << std::endl;
			//if ( wzFlag == false ) continue;

		}

		if(!isData)
		{
			
			//   for(int i=0; i < nGenJets; i++)
			//   {
			//   llp_tree->genJetE[llp_tree->nGenJets] = genJetE[i];
			//   llp_tree->genJetPt[llp_tree->nGenJets] = genJetPt[i];
			//   llp_tree->genJetEta[llp_tree->nGenJets] = genJetEta[i];
			//   llp_tree->genJetPhi[llp_tree->nGenJets] = genJetPhi[i];
			//// llp_tree->nGenJets++;
			//}
			
			llp_tree->genMetPtTrue = genMetPtTrue;
			llp_tree->genMetPhiTrue = genMetPhiTrue;
			llp_tree->genMetPtCalo = genMetPtCalo;
			llp_tree->genMetPhiCalo = genMetPhiCalo;

			if(_debug) std::cout << "nBunchXing " << nBunchXing << std::endl;
			for (int i=0; i < nBunchXing; i++)
			{
				if (BunchXing[i] == 0)
				{
					llp_tree->npu = int(nPUmean[i]);
					if(_debug_npu) 
					{
						std::cout << "npu " << llp_tree->npu << std::endl;
						std::cout << "nPUmean[i] " << nPUmean[i] << std::endl;
					}
				}
			}

			llp_tree->pileupWeight = 1;
			llp_tree->pileupWeight = helper->getPileupWeight(llp_tree->npu);
			llp_tree->pileupWeightUp = helper->getPileupWeightUp(llp_tree->npu) / llp_tree->pileupWeight;
			llp_tree->pileupWeightDown = helper->getPileupWeightDown(llp_tree->npu) / llp_tree->pileupWeight;
			if(_debug_npu) 
			{
				std::cout << "pileupWeightUp " << llp_tree->pileupWeightUp << std::endl;
				std::cout << "pileupWeightDown " << llp_tree->pileupWeightDown << std::endl;
			}
		}

		if(_debug_met) std::cout << "npu " << llp_tree->npu << std::endl;
		if(_debug && llp_tree->npu != 0 ) std::cout << "npu " << llp_tree->npu << std::endl;
		//get NPU
		llp_tree->npv = nPV;
		llp_tree->rho = fixedGridRhoFastjetAll;
		llp_tree->met = metType1Pt;
		if(_debug_met) std::cout << "met " << llp_tree->met << std::endl;
		//if( llp_tree->met < 120. ) continue;
		//if( llp_tree->met < 150. ) continue;
		if( llp_tree->met < 200. ) continue;
		if(_debug_met) std::cout << "metType1Pt passed" << metType1Pt << std::endl;
		llp_tree->metPhi = metType1Phi;
		if(_debug) std::cout << "npv " << llp_tree->npv << std::endl;
		TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM( metType1Pt, 0, metType1Phi, 0 );

		//met filters
		llp_tree->Flag2_globalSuperTightHalo2016Filter          = Flag2_globalSuperTightHalo2016Filter;
		llp_tree->Flag2_globalTightHalo2016Filter               = Flag2_globalTightHalo2016Filter;
		llp_tree->Flag2_goodVertices                            = Flag2_goodVertices;
		llp_tree->Flag2_BadChargedCandidateFilter               = Flag2_BadChargedCandidateFilter;
		llp_tree->Flag2_BadPFMuonFilter                         = Flag2_BadPFMuonFilter;
		llp_tree->Flag2_EcalDeadCellTriggerPrimitiveFilter      = Flag2_EcalDeadCellTriggerPrimitiveFilter;
		llp_tree->Flag2_HBHENoiseFilter                         = Flag2_HBHENoiseFilter;
		llp_tree->Flag2_HBHEIsoNoiseFilter                      = Flag2_HBHEIsoNoiseFilter;
		llp_tree->Flag2_ecalBadCalibFilter                      = Flag2_ecalBadCalibFilter;
		llp_tree->Flag2_eeBadScFilter                           = Flag2_eeBadScFilter;

		//Triggers
		for(int i = 0; i < NTriggersMAX; i++){
			llp_tree->HLTDecision[i] = HLTDecision[i];
		}
		if(_debug_trg) std::cout << "begin: 310 " << HLTDecision[310] << std::endl;
		if(_debug_trg) std::cout << "begin: 310 " << llp_tree->HLTDecision[310] << std::endl;
		if(!HLTDecision[310]) continue; 
		bool triggered = false;
		for(int i = 0; i < NTrigger; i++)
		{
			if(_debug_trg) std::cout << "i " << i << ", NTrigger "<< NTrigger << std::endl;

			int trigger_temp = trigger_paths[i];
			if(_debug_trg) std::cout << "temp  " << trigger_paths[i] << ", triggered "<< triggered << std::endl;

			triggered = triggered || HLTDecision[trigger_temp];
			if(_debug_trg) std::cout << "is triggered ?"<< triggered << std::endl;

		}
		//if (triggered) trig->Fill(1);
		if(_debug) std::cout << "triggered " << triggered << std::endl;
		//*************************************************************************
		//Start Object Selection
		//*************************************************************************
		//sync 
		if(_debug_sync)
		{
			std::cout << "nMuons " << nMuons << std::endl;
			std::cout << "nElectron " << nElectrons << std::endl;
			std::cout << "nPhoton " << nPhotons << std::endl;
			std::cout << "nTau " << nTaus << std::endl;
		}


		std::vector<leptons> Leptons;
		//-------------------------------
		//Muons
		//-------------------------------
		//twiki (https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation)
		for( int i = 0; i < nMuons; i++ )
		{
			if(_debug_sync)
			{
				std::cout << "nMuons " << nMuons << std::endl;
				std::cout << "iMuon " << i << ", Pt " << muonPt[i] << ", Eta " << muonEta[i]<< ", Phi " << muonPhi[i]<< ", E " << muonE[i] << std::endl;
				std::cout << "iMuon " << i << ", muon_chargedIso[i] " << muon_chargedIso[i] << ", muon_neutralHadIso[i] " << muon_neutralHadIso[i]<< ", muon_photonIso[i] " << muon_photonIso[i]<< ", muon_pileupIso[i] " << muon_pileupIso[i] << std::endl;
				std::cout << "iMuon " << i << ", muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i] " << muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i] << ", std::max( stuff, 0)  " << std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) << ", muon_chargedIso[i] + std::max( stuff, 0) " << (muon_chargedIso[i] + std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) ) << ", (muon_chargedIso[i] + std::max( stuff, 0))/muonPt " << (muon_chargedIso[i] + std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) )/muonPt[i] << std::endl;
			}

			if(muonPt[i] <= 10 ) continue;
			if(fabs(muonEta[i]) > 2.4) continue;

			if(!muonIsLoose[i]) continue;
			if( (muon_chargedIso[i] + std::max(muon_neutralHadIso[i] + muon_photonIso[i] - 0.5*muon_pileupIso[i], 0.) )/muonPt[i] >= 0.25) continue;


			//remove overlaps
			bool overlap = false;
			for(auto& lep : Leptons)
			{
				if (RazorAnalyzerLLP::deltaR(muonEta[i],muonPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
			}
			if(overlap) continue;

			leptons tmpMuon;
			tmpMuon.lepton.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i], MU_MASS);
			tmpMuon.pdgId = 13 * -1 * muonCharge[i];

			Leptons.push_back(tmpMuon);

			llp_tree->muon_pileupIso[llp_tree->nMuons] = muon_pileupIso[i];
			llp_tree->muon_chargedIso[llp_tree->nMuons] = muon_chargedIso[i];
			llp_tree->muon_photonIso[llp_tree->nMuons] = muon_photonIso[i];
			llp_tree->muon_neutralHadIso[llp_tree->nMuons] = muon_neutralHadIso[i];

			llp_tree->muonIsLoose[llp_tree->nMuons] = muonIsLoose[i];

			llp_tree->muonPt[llp_tree->nMuons] = muonPt[i];
			llp_tree->muonEta[llp_tree->nMuons] = muonEta[i];
			llp_tree->muonE[llp_tree->nMuons] = muonE[i];
			llp_tree->muonPhi[llp_tree->nMuons] = muonPhi[i];

			llp_tree->nMuons++;
		}
		if(_debug) std::cout << "nMuons " << llp_tree->nMuons << std::endl;

		//-------------------------------
		//Electrons
		//-------------------------------
		for( int i = 0; i < nElectrons; i++ )
		{
			if(_debug_sync)
			{
				std::cout << "nElectrons " << nElectrons << std::endl;
				std::cout << "iElectron " << i << ", Pt " << elePt[i] << ", Eta " << eleEta[i]<< ", Phi " << elePhi[i]<< ", E " << eleE[i] << std::endl;
				std::cout << "iElectron " << i << ", ele_passCutBasedIDVeto[i] " << ele_passCutBasedIDVeto[i]<< std::endl;
			}

			if(elePt[i] <= 10 ) continue;
			if(fabs(eleEta[i]) > 2.5) continue;

			if(!ele_passCutBasedIDVeto[i]) continue;

			//remove overlaps
			bool overlap = false;
			for(auto& lep : Leptons)
			{
				if (RazorAnalyzerLLP::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
			}
			if(overlap) continue;

			leptons tmpElectron;
			tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);

			Leptons.push_back(tmpElectron);

			llp_tree->ele_passCutBasedIDVeto[llp_tree->nElectrons] = ele_passCutBasedIDVeto[i];

			llp_tree->elePt[llp_tree->nElectrons] = elePt[i];
			llp_tree->eleEta[llp_tree->nElectrons] = eleEta[i];
			llp_tree->eleE[llp_tree->nElectrons] = eleE[i];
			llp_tree->elePhi[llp_tree->nElectrons] = elePhi[i];

			llp_tree->nElectrons++;
		}
		if(_debug) std::cout << "nElectrons " << llp_tree->nElectrons << std::endl;

		std::vector<taus> Taus;
		//-------------------------------
		//Taus
		//-------------------------------
		for( int i = 0; i < nTaus; i++ )
		{
			if(_debug_sync)
			{
				std::cout << "nTaus " << nTaus << std::endl;
				std::cout << "iTau " << i << ", Pt " << tauPt[i] << ", Eta " << tauEta[i]<< ", Phi " << tauPhi[i]<< ", E " << tauE[i] << std::endl;
				std::cout << "iTau " << i << ", tau_IsLoose[i] " << tau_IsLoose[i]<< std::endl;
			}

			if(tauPt[i] <= 18 ) continue;
			if(fabs(tauEta[i]) > 2.3) continue;

			if(!tau_IsLoose[i]) continue;

			//remove overlaps
			bool overlap = false;
			for(auto& lep : Leptons)
			{
				if (RazorAnalyzerLLP::deltaR(tauEta[i],tauPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
			}
			if(overlap) continue;

			bool overlap_tau = false;
			for(auto& tau : Taus)
			{
				if (RazorAnalyzerLLP::deltaR(tauEta[i],tauPhi[i],tau.tau.Eta(),tau.tau.Phi()) < 0.3) overlap_tau = true;
			}
			if(overlap_tau) continue;

			taus tmpTau;
			tmpTau.tau.SetPtEtaPhiM(tauPt[i],tauEta[i], tauPhi[i], TAU_MASS);

			Taus.push_back(tmpTau);

			llp_tree->tau_IsLoose[llp_tree->nTaus] = tau_IsLoose[i];

			llp_tree->tauPt[llp_tree->nTaus] = tauPt[i];
			llp_tree->tauEta[llp_tree->nTaus] = tauEta[i];
			llp_tree->tauE[llp_tree->nTaus] = tauE[i];
			llp_tree->tauPhi[llp_tree->nTaus] = tauPhi[i];

			llp_tree->nTaus++;
		}
		if(_debug) std::cout << "nTaus " << llp_tree->nTaus << std::endl;

		std::vector<photons> Photons;
		//-------------------------------
		//Photons
		//-------------------------------
		for( int i = 0; i < nPhotons; i++ )
		{
			if(_debug_sync)
			{
				std::cout << "nPhotons " << nPhotons << std::endl;
				std::cout << "iPhoton " << i << ", Pt " << phoPt[i] << ", Eta " << phoEta[i]<< ", Phi " << phoPhi[i]<< ", E " << phoE[i] << std::endl;
				std::cout << "iPhoton " << i << ", pho_passCutBasedIDLoose[i] " << pho_passCutBasedIDLoose[i]<< std::endl;
			}

			if(phoPt[i] <= 15 ) continue;
			if(fabs(phoEta[i]) > 2.5) continue;

			if(!pho_passCutBasedIDLoose[i]) continue;

			//remove overlaps
			bool overlap = false;
			for(auto& lep : Leptons)
			{
				if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
			}
			if(overlap) continue;

			bool overlap_tau = false;
			for(auto& tau : Taus)
			{
				if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],tau.tau.Eta(),tau.tau.Phi()) < 0.3) overlap_tau = true;
			}
			if(overlap_tau) continue;

			bool overlap_pho = false;
			for(auto& pho : Photons)
			{
				if (RazorAnalyzerLLP::deltaR(phoEta[i],phoPhi[i],pho.photon.Eta(),pho.photon.Phi()) < 0.3) overlap_pho = true;
			}
			if(overlap_pho) continue;

			photons tmpPhoton;
			tmpPhoton.photon.SetPtEtaPhiM(phoPt[i],phoEta[i], phoPhi[i], 0.);

			Photons.push_back(tmpPhoton);

			llp_tree->pho_passCutBasedIDLoose[llp_tree->nPhotons] = pho_passCutBasedIDLoose[i];

			llp_tree->phoPt[llp_tree->nPhotons] = phoPt[i];
			llp_tree->phoEta[llp_tree->nPhotons] = phoEta[i];
			llp_tree->phoE[llp_tree->nPhotons] = phoE[i];
			llp_tree->phoPhi[llp_tree->nPhotons] = phoPhi[i];

			llp_tree->nPhotons++;
		}
		if(_debug) std::cout << "nPhotons " << llp_tree->nPhotons << std::endl;

		//-------------------------------
		//Leptons
		//-------------------------------
  	 	TLorentzVector lepp4;
		sort(Leptons.begin(), Leptons.end(), my_largest_pt_lep);
		for ( auto &tmp : Leptons )
		{
			llp_tree->lepE[llp_tree->nLeptons]      = tmp.lepton.E();
			llp_tree->lepPt[llp_tree->nLeptons]     = tmp.lepton.Pt();
			llp_tree->lepEta[llp_tree->nLeptons]    = tmp.lepton.Eta();
			llp_tree->lepPhi[llp_tree->nLeptons]    = tmp.lepton.Phi();
			llp_tree->lepPdgId[llp_tree->nLeptons]  = tmp.pdgId;
			//llp_tree->lepDZ[llp_tree->nLeptons]     = tmp.dZ;
			//llp_tree->lepPassId[llp_tree->nLeptons] = tmp.passId;
			if(_debug) std::cout << "lepE " << tmp.lepton.E() << std::endl;

			// std::cout << "lepton pdg " << llp_tree->lepPdgId[llp_tree->nLeptons] << std::endl;
			llp_tree->nLeptons++;

			lepp4 += tmp.lepton;
		}
		float dPhi = RazorAnalyzerLLP::deltaPhi(llp_tree->metPhi, lepp4.Phi());
		llp_tree->MT = sqrt(2*(llp_tree->met)*lepp4.Pt()*(1-cos(dPhi)));
		if(_debug_lep) std::cout<<"MT "<<llp_tree->MT<<std::endl;
		//if (triggered) trig_lepId->Fill(1);
		if(_debug) std::cout << "nLeptons " << llp_tree->nLeptons << std::endl;


		//-----------------------------------------------
		//Select Jets
		//-----------------------------------------------
		//std::vector<double> jetPtVector;
		//std::vector<double> jetCISVVector;
		std::vector<jets> Jets;
		//auto highest = [](auto a, auto b) { return a > b; };
		//cout <<"nJets :" << nJets << std::endl;

		if(_debug) std::cout << "nJets " << nJets << std::endl;
		if(_debug_jet) std::cout << "nJets " << nJets << std::endl;
		if(_debug_jet) std::cout << "jetE 0 " << jetE[0] << std::endl;
		//if(_debug_jet) std::cout << "jetChargedEMEnergyFraction 0 " << jetChargedEMEnergyFraction[0] << std::endl;
		//if(_debug_jet) std::cout << "jetNeutralEMEnergyFraction 0 " << jetNeutralEMEnergyFraction[0] << std::endl;
		if(_debug_jet) std::cout << "jetGammaMax_ET 0 " << jetGammaMax_ET[0] << std::endl;

		float ht = 0.;

		for(int i = 0; i < nJets; i++)
		{

			ht += jetPt[i];

			//------------------------------------------------------------
			//exclude selected muons and electrons from the jet collection
			//------------------------------------------------------------
			double deltaR = -1;
			for(auto& lep : Leptons){
				double thisDR = RazorAnalyzerLLP::deltaR(jetEta[i],jetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
				if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
			}
			if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

			//------------------------------------------------------------
			//exclude selected taus from the jet collection
			//------------------------------------------------------------
			double deltaR_tau = -1;
			for(auto& tau : Taus){
				double thisDR_tau = RazorAnalyzerLLP::deltaR(jetEta[i],jetPhi[i],tau.tau.Eta(),tau.tau.Phi());
				if(deltaR_tau < 0 || thisDR_tau < deltaR_tau) deltaR_tau = thisDR_tau;
			}
			if(deltaR_tau > 0 && deltaR_tau < 0.4) continue; //jet matches a selected tau

			//------------------------------------------------------------
			//exclude selected photons from the jet collection
			//------------------------------------------------------------
			double deltaR_pho = -1;
			for(auto& pho : Photons){
				double thisDR_pho = RazorAnalyzerLLP::deltaR(jetEta[i],jetPhi[i],pho.photon.Eta(),pho.photon.Phi());
				if(deltaR_pho < 0 || thisDR_pho < deltaR_pho) deltaR_pho = thisDR_pho;
			}
			if(deltaR_pho > 0 && deltaR_pho < 0.4) continue; //jet matches a selected photon

			//------------------------------------------------------------
			//Apply Jet Energy and Resolution Corrections
			//------------------------------------------------------------

			TLorentzVector thisJet = makeTLorentzVector( jetPt[i], jetEta[i], jetPhi[i], jetE[i] );
			//double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
			//		fixedGridRhoFastjetAll, jetJetArea[i] , JetCorrector);
		
			//double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
			//		fixedGridRhoFastjetAll, jetJetArea[i],
			//		runNum,
			//		JetCorrectorIOV,JetCorrector);
			////cout <<"JEC :" << JEC << std::endl;

			//TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );
		
			if( thisJet.Pt() < 30 ) continue;//According to the April 1st 2015 AN
			//if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
			if( fabs( thisJet.Eta() ) >= 2.4 ) continue;
			//if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
			// if ( !jetPassIDLoose[i] ) continue;
			// if (!(jetRechitE[i] > 0.0)) continue;
			// if(jetNRechits[i]<10) continue;
			// if(jetRechitT[i] < 0.0) continue;
			// if ((jetChargedHadronEnergyFraction[i]+jetChargedEMEnergyFraction[i]) > 0.4) continue;
			// if ((jetChargedHadronEnergyFraction[i]+jetNeutralHadronEnergyFraction[i])/(jetChargedEMEnergyFraction[i]+jetNeutralEMEnergyFraction[i]) < 0.2) continue;

			// std::cout <<jetRechitT[i] << "," << jetRechitE[i] <<  "," << jetNRechits[i] << std::endl;

			jets tmpJet;
			tmpJet.jet    = thisJet;
			tmpJet.time   = jetRechitT[i];
			tmpJet.ecalNRechits = jetNRechits[i];
			tmpJet.ecalRechitE = jetRechitE[i];
			tmpJet.jetNeutralHadronMultiplicity = jetNeutralHadronMultiplicity[i];
			tmpJet.jetChargedHadronMultiplicity = jetChargedHadronMultiplicity[i];
			tmpJet.jetMuonMultiplicity = jetMuonMultiplicity[i];
			tmpJet.jetElectronMultiplicity = jetElectronMultiplicity[i];
			tmpJet.jetPhotonMultiplicity = jetPhotonMultiplicity[i];
			tmpJet.jetNeutralHadronEnergyFraction = jetNeutralHadronEnergyFraction[i];
			tmpJet.jetChargedHadronEnergyFraction = jetChargedHadronEnergyFraction[i];
			tmpJet.jetMuonEnergyFraction = jetMuonEnergyFraction[i];
			tmpJet.jetElectronEnergyFraction = jetElectronEnergyFraction[i];
			tmpJet.jetPhotonEnergyFraction = jetPhotonEnergyFraction[i];
			tmpJet.jetCSV = jetCISV[i];

			tmpJet.jetGammaMax_ET = jetGammaMax_ET[i];
			tmpJet.jetMinDeltaRPVTracks = jetMinDeltaRPVTracks[i];

			tmpJet.jetChargedEMEnergyFraction = jetChargedEMEnergyFraction[i];
			tmpJet.jetNeutralEMEnergyFraction = jetNeutralEMEnergyFraction[i];

			tmpJet.jetNVertexTracks = jetNVertexTracks[i];
			tmpJet.jetNSelectedTracks = jetNSelectedTracks[i];
			tmpJet.jetDRSVJet = jetDRSVJet[i];
			tmpJet.jetSVMass = jetSVMass[i];

			tmpJet.jetChargedMultiplicity = jetChargedHadronMultiplicity[i]+jetElectronMultiplicity[i]+jetMuonMultiplicity[i];

			if(_debug_trk) std::cout << "nTracks" << nTracks << std::endl;
			std::vector<float> nPixelHits;
			std::vector<float> nHits;
			//int nChargedMultiplicity=0;
			for(int iTrack=0; iTrack<nTracks; iTrack++)
			{
				double thisDR_trk = RazorAnalyzerLLP::deltaR(jetEta[i],jetPhi[i],track_Eta[iTrack], track_Phi[iTrack]);
				if(thisDR_trk<0.4)
				{
					//nChargedMultiplicity++;
					nPixelHits.push_back(track_nPixelHits[iTrack]);
					nHits.push_back(track_nHits[iTrack]);
				}

			}

			std::sort(nPixelHits.begin(), nPixelHits.end());
			float nPixelHitsMedian = 0;
			if (nPixelHits.size() > 0) {
				if (nPixelHits.size() % 2 ==0) nPixelHitsMedian = ((nPixelHits[nPixelHits.size()/2 -1] + nPixelHits[nPixelHits.size()/2]) /2);
				else nPixelHitsMedian = nPixelHits[nPixelHits.size()/2];
			}    
			std::sort(nHits.begin(), nHits.end());
			float nHitsMedian = 0;
			if (nHits.size() > 0) {
				if (nHits.size() % 2 ==0) nHitsMedian = ((nHits[nHits.size()/2 -1] + nHits[nHits.size()/2]) /2);
				else nHitsMedian = nHits[nHits.size()/2];
			}

			//tmpJet.jetChargedMultiplicity = nChargedMultiplicity;
			tmpJet.jetNPixelHitsMedian = nPixelHitsMedian;
			tmpJet.jetNHitsMedian = nHitsMedian;

			Jets.push_back(tmpJet);

		}
		llp_tree->HT = ht;


		//-----------------------------
		//Require at least 2 jets
		//-----------------------------
		///
		//   if(pf)
		//   {
		//   if( Jets.size() < 1 ) continue;

		//   }
		//   /
		//if (triggered) trig_lepId_dijet->Fill(1);
		sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

		if (Jets.size()>0)
		{
			llp_tree->jetMet_dPhi = RazorAnalyzerLLP::deltaPhi(jetPhi[0],metType1Phi);
			//TLorentzVector t1PFMET = makeTLorentzVectorPtEtaPhiM( metType1Pt, 0, metType1Phi, 0 );
			TLorentzVector jet0 = makeTLorentzVectorPtEtaPhiM( jetPt[0], 0, jetPhi[0], 0 );
			llp_tree->jetMet_dPhiStar = RazorAnalyzerLLP::deltaPhi(jetPhi[0],  (t1PFMET+jet0).Phi() );
		}
		else{
			llp_tree->jetMet_dPhi = -999.;
			llp_tree->jetMet_dPhiStar = -999.;
		}
		float jetMet_dPhiMin_temp = 999 ; 
		float jetMet_dPhiStarMin_temp = 999 ; 
		float jetMet_dPhiMin4_temp = 999 ; 

		for ( auto &tmp : Jets )
		{
			llp_tree->jetNeutralHadronMultiplicity[llp_tree->nJets] = tmp.jetNeutralHadronMultiplicity;
			llp_tree->jetChargedHadronMultiplicity[llp_tree->nJets] = tmp.jetChargedHadronMultiplicity;
			llp_tree->jetMuonMultiplicity[llp_tree->nJets] = tmp.jetMuonMultiplicity;
			llp_tree->jetElectronMultiplicity[llp_tree->nJets] = tmp.jetElectronMultiplicity;
			llp_tree->jetPhotonMultiplicity[llp_tree->nJets] = tmp.jetPhotonMultiplicity;
			llp_tree->jetNeutralHadronEnergyFraction[llp_tree->nJets] = tmp.jetNeutralHadronEnergyFraction;
			llp_tree->jetChargedHadronEnergyFraction[llp_tree->nJets] = tmp.jetChargedHadronEnergyFraction;
			llp_tree->jetMuonEnergyFraction[llp_tree->nJets] = tmp.jetMuonEnergyFraction;
			llp_tree->jetElectronEnergyFraction[llp_tree->nJets] = tmp.jetElectronEnergyFraction;
			llp_tree->jetPhotonEnergyFraction[llp_tree->nJets] = tmp.jetPhotonEnergyFraction;
			llp_tree->jetCSV[llp_tree->nJets] = tmp.jetCSV;

			llp_tree->jetPt[llp_tree->nJets] = tmp.jet.Pt();
			llp_tree->jetEta[llp_tree->nJets] = tmp.jet.Eta();
			llp_tree->jetE[llp_tree->nJets] = tmp.jet.E();
			llp_tree->jetPhi[llp_tree->nJets] = tmp.jet.Phi();
			llp_tree->jetTime[llp_tree->nJets] = tmp.time;
			llp_tree->ecalNRechits[llp_tree->nJets] = tmp.ecalNRechits;
			llp_tree->ecalRechitE[llp_tree->nJets] = tmp.ecalRechitE;

			llp_tree->jetGammaMax_ET[llp_tree->nJets] = tmp.jetGammaMax_ET;
			llp_tree->jetMinDeltaRPVTracks[llp_tree->nJets] = tmp.jetMinDeltaRPVTracks;

			llp_tree->jetChargedEMEnergyFraction[llp_tree->nJets] = tmp.jetChargedEMEnergyFraction;
			llp_tree->jetNeutralEMEnergyFraction[llp_tree->nJets] = tmp.jetNeutralEMEnergyFraction;

			llp_tree->jetEcalE[llp_tree->nJets] = tmp.jet.E()*(tmp.jetElectronEnergyFraction+tmp.jetPhotonEnergyFraction);
			llp_tree->jetHcalE[llp_tree->nJets] = tmp.jet.E()*(tmp.jetNeutralHadronEnergyFraction+tmp.jetChargedHadronEnergyFraction);

			llp_tree->jetNVertexTracks[llp_tree->nJets] = tmp.jetNVertexTracks;
			llp_tree->jetNSelectedTracks[llp_tree->nJets] = tmp.jetNSelectedTracks;
			llp_tree->jetDRSVJet[llp_tree->nJets] = tmp.jetDRSVJet;
			llp_tree->jetSVMass[llp_tree->nJets] = tmp.jetSVMass;

			llp_tree->jetChargedMultiplicity[llp_tree->nJets] = tmp.jetChargedMultiplicity;
			llp_tree->jetNPixelHitsMedian[llp_tree->nJets] = tmp.jetNPixelHitsMedian;
			llp_tree->jetNHitsMedian[llp_tree->nJets] = tmp.jetNHitsMedian;

			//std::cout << "jetEta " << tmp.jet.Eta() << std::endl;
			//std::cout << "jetEta " << llp_tree->jetEta[llp_tree->nJets] << std::endl;

			if(jetMet_dPhiMin4_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi)) && llp_tree->nJets < 4)
			{
				jetMet_dPhiMin4_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi));
			}
			if (jetMet_dPhiMin_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi)))
			{
				jetMet_dPhiMin_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(),metType1Phi));
			}     
			TLorentzVector jet_temp = makeTLorentzVectorPtEtaPhiM( tmp.jet.Pt(), 0, tmp.jet.Phi(), 0 );
			if (jetMet_dPhiStarMin_temp > abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(), (t1PFMET+jet_temp).Phi() )))
			{
				jetMet_dPhiStarMin_temp = abs(RazorAnalyzerLLP::deltaPhi(tmp.jet.Phi(), (t1PFMET+jet_temp).Phi() ));
			}     

			llp_tree->nJets++;
		}
		if(_debug) std::cout << "nJets in tree " << llp_tree->nJets << std::endl;
		llp_tree->jetMet_dPhiMin = jetMet_dPhiMin_temp;
		llp_tree->jetMet_dPhiStarMin = jetMet_dPhiStarMin_temp;
		llp_tree->jetMet_dPhiMin4 = jetMet_dPhiMin4_temp;



	llp_tree->tree_->Fill();
	if(_debug) std::cout << "nJets in tree " << llp_tree->nJets << std::endl;
	}

	cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
	outFile->cd();
	cout << "Writing output trees..." << endl;
	llp_tree->tree_->Write();
	cout << "Writing output NEvents..." << endl;
	NEvents->Write();
	cout << "Closing output trees..." << endl;
	outFile->Write();
	outFile->Close();

	delete helper;
}
