//Class for analyzing ntuples produced by the llp_ntupler framework
//
//Author: Cristian Pena & Si Xie

#ifndef RazorAnalyzer_h
#define RazorAnalyzer_h

#include "llp_event.h" //This is a MakeClass of the llp tree in the ntuple to be analyzed
#include "FactorizedJetCorrector.h"
#include "SimpleJetResolution.h"

//ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include "TLorentzVector.h"
#include "TRandom3.h"

//C++ includes
#include <map>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

class RazorAnalyzer: public llp_event {
    public :
        RazorAnalyzer(TTree *tree=0);
        virtual ~RazorAnalyzer();

        void EnableEventInfo();
        void EnablePVAll();
        void EnableMuons();
        void EnableElectrons();
        void EnableTaus();
        void EnableIsoPFCandidates();
        void EnablePhotons();
        void EnableJets();
        void EnableCaloJets();
        void EnableFatJets();
        void EnableMet();
        void EnablePileup();
        void EnableMC();
        void EnableLLP();
        void EnableGenParticles();
        void EnableRazor();
        void EnableCSC();
        void EnableDT();

        void EnableEcalRechits();
        void EnableAll();
        void EnableAllWithEcalRechits();

        //------ LIST OF ANALYSES ------//
        virtual void Analyze(bool isData, int option, string outputFileName, string label);

        //functions in RazorAuxMuon.cc
	float GetMuonEffectiveAreaMean(int i, string type );
	float GetMuonEffectiveArea90(int i, string EraName );
	bool isMuonPOGLooseMuon(int i, bool applyID = true, bool applyIso = true);
	bool isMuonPOGMediumMuon(int i, bool applyID = true, bool applyIso = true);
        bool isMuonPOGTightMuon(int i, bool applyID = true, bool applyIso = true);
	bool isVetoMuon(int i, bool applyID = true, bool applyIso = true);
	bool isLooseMuon(int i, bool applyID = true, bool applyIso = true);
        bool isTightMuon(int i, bool applyID = true, bool applyIso = true);
        bool passHZZMuonPreselection(int i);
        bool isHZZMuon(int i);
	bool matchMuonHLTFilters( int i, string HLTFilter);
	bool matchTagMuonHLTFilters( int i);

        //functions in RazorAuxElectron.cc
	float GetElectronScaleCorrection( double pt, double eta );
	float GetElectronEffectiveAreaMean(int i, bool use25nsCuts = true);
	float GetElectronEffectiveArea90(int i, string EraName = "Spring15");
        bool isEGammaPOGVetoElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = true, string EraName = "Spring15");
        bool isEGammaPOGLooseElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = true, string EraName = "Spring15");
        bool isEGammaPOGMediumElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = true, string EraName = "Spring15");
        bool isEGammaPOGTightElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = true, string EraName = "Spring15");
        bool isVetoElectron(int i, bool applyID = true, bool applyIso = true, string EraName = "Spring15");
        bool isLooseElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = true, string EraName = "Spring15");
        bool isMediumElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = true, string EraName = "Spring15");
        bool isTightElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = true, string EraName = "Spring15");
	bool isMVANonTrigVetoElectron(int i, bool applyID = true, bool applyIso = true, string EraName = "Spring15");
        bool passEGammaPOGVetoElectronID(int i, bool use25nsCuts = true, string EraName = "Spring15");
        bool passEGammaPOGLooseElectronID(int i, bool use25nsCuts = true, string EraName = "Spring15");
        bool passEGammaPOGMediumElectronID(int i, bool use25nsCuts = true, string EraName = "Spring15");
        bool passEGammaPOGTightElectronID(int i, bool use25nsCuts = true, string EraName = "Spring15");
	bool passMVANonTrigVetoElectronID(int i, string EraName = "Spring15");
	bool passMVAVetoElectronID(int i, string EraName = "Spring15");
	bool passMVALooseElectronID(int i, string EraName = "Spring15");
        bool passEGammaPOGVetoElectronIso(int i, bool use25nsCuts = true);
        bool passEGammaPOGLooseElectronIso(int i, bool use25nsCuts = true);
        bool passEGammaPOGMediumElectronIso(int i, bool use25nsCuts = true);
        bool passEGammaPOGTightElectronIso(int i, bool use25nsCuts = true);
	bool passMVANonTrigVetoElectronIso(int i);
	bool passHZZElectronIso(int i);
	bool passHZZElectronPreselection(int i);
	bool isHZZElectron(int i);
	bool matchElectronHLTFilters( int i, string HLTFilter, string analysisTag);
	bool matchElectronHLTFilters2015( int i, string HLTFilter);
	bool matchElectronHLTFilters2016( int i, string HLTFilter);
	bool matchProbeElectronHLTFilters( int i);
	bool matchProbeSCHLTFilters( int i);
	bool matchTagElectronHLTFilters( int i);

        //functions in RazorAuxTau.cc
        bool isLooseTau(int i);
        bool isMediumTau(int i);
        bool isTightTau(int i);

        //functions in RazorAuxPhoton.cc
        bool photonPassesElectronVeto(int i);
	double getPhotonSminorSmajor(int ind_pho, bool isSminor = true);
	void getPhotonEffAreaRun2( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho );
	void getPhotonEffArea90( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho );
	void getPhotonEffAreaPFClusterIso( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho );
        bool photonPassesIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut, bool useEffectiveArea90, bool usePrivatePF = false, bool usePFClusterIso = false);
	bool photonPassLooseIDWithoutEleVeto(int i, bool use25nsCuts = true);
	bool photonPassMediumIDWithoutEleVeto(int i, bool use25nsCuts = true);
	bool photonPassTightIDWithoutEleVeto(int i, bool use25nsCuts = true);
	bool photonPassLooseDelayedIDWithoutEleVeto(int i, bool use25nsCuts = true);
	bool photonPassMediumDelayedIDWithoutEleVeto(int i, bool use25nsCuts = true);
	bool photonPassTightDelayedIDWithoutEleVeto(int i, bool use25nsCuts = true);
	bool photonPassLooseID(int i, bool use25nsCuts = true);
	bool photonPassMediumID(int i, bool use25nsCuts = true);
	bool photonPassTightID(int i, bool use25nsCuts = true);
	bool photonPassLooseIso(int i, bool use25nsCuts = true, bool usePrivatePF = false, bool usePFClusterIso = false);
	//bool photonPassMediumIso(int i, bool use25nsCuts = true, bool usePrivatePF = false, bool usePFClusterIso = false);
  bool photonPassMediumIso(int i, bool use25nsCuts = true, bool usePrivatePF = false);
	bool photonPassTightIso(int i, bool use25nsCuts = true, bool usePrivatePF = false, bool usePFClusterIso = false);
        bool isLoosePhoton(int i, bool use25nsCuts = true);
        bool isMediumPhoton(int i, bool use25nsCuts = true);
        bool isTightPhoton(int i, bool use25nsCuts = true);
        bool isLoosePhotonWithoutEleVeto(int i, bool use25nsCuts = true);
        bool isMediumPhotonWithoutEleVeto(int i, bool use25nsCuts = true);
        bool isTightPhotonWithoutEleVeto(int i, bool use25nsCuts = true);
        bool isLooseDelayedPhotonWithoutEleVeto(int i, bool use25nsCuts = true);
        bool isMediumDelayedPhotonWithoutEleVeto(int i, bool use25nsCuts = true);
        bool isTightDelayedPhotonWithoutEleVeto(int i, bool use25nsCuts = true);
	bool matchPhotonHLTFilters( int i, string HLTFilter);
	void getPhotonEffAreaExo15004( float eta, double& effAreaPho );
	bool photonPassLooseIDWithoutEleVetoExo15004(int i);
	bool photonPassesIsolationExo15004(int i, double PFChHadIsoCut, double PFPhotIsoCut );
	bool photonPassLooseIsoExo15004(int i);
	TLorentzVector GetCorrectedMomentum( TVector3 vtx, TVector3 phoPos, double phoE );

        //functions for 2017 photon ID 92X
	void getPhotonEffAreaRun2_2017( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho );
	void getPhotonEffArea90_2017( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho );
        bool photonPassesIsolation_2017(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut, bool useEffectiveArea90, bool usePrivatePF = false);
	bool photonPassLooseIDWithoutEleVeto_2017(int i, bool use25nsCuts = true);
	bool photonPassMediumIDWithoutEleVeto_2017(int i, bool use25nsCuts = true);
	bool photonPassTightIDWithoutEleVeto_2017(int i, bool use25nsCuts = true);
	bool photonPassLooseID_2017(int i, bool use25nsCuts = true);
	bool photonPassMediumID_2017(int i, bool use25nsCuts = true);
	bool photonPassTightID_2017(int i, bool use25nsCuts = true);
	bool photonPassLooseIso_2017(int i, bool use25nsCuts = true, bool usePrivatePF = false);
	bool photonPassMediumIso_2017(int i, bool use25nsCuts = true, bool usePrivatePF = false);
	bool photonPassTightIso_2017(int i, bool use25nsCuts = true, bool usePrivatePF = false);
        bool isLoosePhoton_2017(int i, bool use25nsCuts = true);
        bool isMediumPhoton_2017(int i, bool use25nsCuts = true);
        bool isTightPhoton_2017(int i, bool use25nsCuts = true);
        bool isLoosePhotonWithoutEleVeto_2017(int i, bool use25nsCuts = true);
        bool isMediumPhotonWithoutEleVeto_2017(int i, bool use25nsCuts = true);
        bool isTightPhotonWithoutEleVeto_2017(int i, bool use25nsCuts = true);

        //functions for 2017 photon ID 94X
	void getPhotonEffAreaRun2_94X( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho );
	void getPhotonEffArea90_94X( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho );
        bool photonPassesIsolation_94X(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut, bool useEffectiveArea90, bool usePrivatePF = false);
	bool photonPassLooseIDWithoutEleVeto_94X(int i, bool use25nsCuts = true);
	bool photonPassMediumIDWithoutEleVeto_94X(int i, bool use25nsCuts = true);
	bool photonPassTightIDWithoutEleVeto_94X(int i, bool use25nsCuts = true);
	bool photonPassLooseID_94X(int i, bool use25nsCuts = true);
	bool photonPassMediumID_94X(int i, bool use25nsCuts = true);
	bool photonPassTightID_94X(int i, bool use25nsCuts = true);
	bool photonPassLooseIso_94X(int i, bool use25nsCuts = true, bool usePrivatePF = false);
	bool photonPassMediumIso_94X(int i, bool use25nsCuts = true, bool usePrivatePF = false);
	bool photonPassTightIso_94X(int i, bool use25nsCuts = true, bool usePrivatePF = false);
        bool isLoosePhoton_94X(int i, bool use25nsCuts = true);
        bool isMediumPhoton_94X(int i, bool use25nsCuts = true);
        bool isTightPhoton_94X(int i, bool use25nsCuts = true);
        bool isLoosePhotonWithoutEleVeto_94X(int i, bool use25nsCuts = true);
        bool isMediumPhotonWithoutEleVeto_94X(int i, bool use25nsCuts = true);
        bool isTightPhotonWithoutEleVeto_94X(int i, bool use25nsCuts = true);
	/* //function in HggRazorAuxPhoton.cc */
	/* // R u n 1   C u t   B a s e d   I D */
	/* //---------------------------------- */
	/* bool isGoodPhotonRun1( int i, bool _iso, bool _debug ); */
	/* bool photonPassIsoRun1( int i, bool _debug ); */
	/* bool photonPassIsoRun1( int i , WP wp, bool _debug ); */
	/* void getPhotonEffAreaRun1( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho ); */
	/* bool passEleVetoRun1( int i ); */
	/* // R u n 2   C u t   B a s e d   I D */
	/* //---------------------------------- */
	/* bool isGoodPhotonRun2( int i, bool _iso, WP wp, bool _debug ); */
	/* bool photonPassIsoRun2( int i, WP wp ,bool _debug ); */
	/* void getPhotonEffAreaRun2( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho ); */


        //functions in RazorAuxJet.cc
        bool isCSVL(int i, string dataset = "80X");
        bool isCSVM(int i, string dataset = "80X");
        bool isCSVT(int i, string dataset = "80X");
	double JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
					  double rho, double jetArea,
					  int run,
					  std::vector<std::pair<int,int> > JetCorrectionsIOV,
					  std::vector<FactorizedJetCorrector*> jetcorrector,
					  int jetCorrectionLevel = -1,
					  bool printDebug = false);
	double JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
					  double rho, double jetArea,
					  FactorizedJetCorrector* jetcorrector,
					  int jetCorrectionLevel = -1,
					  bool printDebug = false);
	double JetEnergySmearingFactor( double jetPt, double jetEta, double NPU,
  					SimpleJetResolution *JetResolutionCalculator,
                                        TRandom3 *random);
        double UpDownJetEnergySmearingFactor(double unsmearedPt, double jetEta, double NPU,
                                             SimpleJetResolution *JetResolutionCalculator,
                                             double smearedPt, string option);
        double BTagScaleFactor( double jetPt, bool CSVM, string option="");

        //functions in RazorAuxMisc.cc
	double deltaPhi(double phi1, double phi2);
	double deltaR(double eta1, double phi1, double eta2, double phi2);
        TLorentzVector makeTLorentzVector(double pt, double eta, double phi, double energy);
	TLorentzVector makeTLorentzVectorPtEtaPhiM(double pt, double eta, double phi, double mass);
	vector<TLorentzVector> getHemispheres(vector<TLorentzVector> jets);
	std::vector< std::vector<int> > getHemispheresV2( std::vector<TLorentzVector> jets);

    double computeMR(TLorentzVector hem1, TLorentzVector hem2);
        double computeRsq(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector met);
	double GetMT( TLorentzVector visible, TVector3 met );
	double GetMTEnergy( TLorentzVector visible, TVector3 met );
	double GetMT( TLorentzVector visible, TLorentzVector met );
	double GetMTEnergy( TLorentzVector visible, TLorentzVector met );
	double GetDphi( TLorentzVector visible, TVector3 met );
	double GetDphi( TLorentzVector visible, TLorentzVector met );
        bool matchesVetoLepton(float eta, float phi, float dR=0.4);

    double GetAlphaT(vector<TLorentzVector> jets) ;
    double GetDPhiMin(vector<TLorentzVector> jets);

        bool passesHadronicRazorBaseline(double MR, double Rsq);
        bool passesLeptonicRazorBaseline(double MR, double Rsq);
        int SubtractParticleFromCollection(TLorentzVector ToSubtract, vector<TLorentzVector>& Collection, float deltaRMatch=0.4);

	double calcMT2(float testMass, bool massive, std::vector<TLorentzVector> jets, TLorentzVector MET, int hemi_seed, int hemi_association);

	//functions in src/RazorAuxGenLevel.cc
	bool matchesGenMuon(double eta, double phi);
	bool matchesGenElectron(double eta, double phi);
        bool isHadronicDecay(int index, int daughterStatus=23);
        int getMatchingHardProcessParticleIndex(double eta, double phi,
                int id, int status=22, double r=0.8);
        int getMatchingGenWIndex(double eta, double phi, double r=0.8);
        int getMatchingGenTopIndex(double eta, double phi, double r=0.8);
	bool isGenTau(int index);
	bool isGenLeptonicTau(int index);
	int findClosestGenElectron(double eta, double phi);
	int findClosestGenMuon(double eta, double phi);
	int findClosestGenTau(double eta, double phi);
	int findClosestRecoTau(double eta, double phi);
	int GetTauMatchedID(double eta, double phi);
	int findClosestParton(float eta, float phi);
	double getGenHT();
        int getNISR( std::vector<FactorizedJetCorrector*> &JetCorrector, std::vector<std::pair<int,int> > &JetCorrectorIOV ); //count number of gen-level ISR jets

	//Added to src/RazorAuxGenLevel.cc
	int findClosestGenJet(double eta, double phi);

	//conversion between DetId <-> ieta/ix/iphi/iy

	int detID_from_iEtaiPhi(int iEta_or_iX, int iPhi_or_iY, bool isEB, bool isEEMinus);
	int iEta_or_iX_from_detID(int detID, bool isEB);
	int iPhi_or_iY_from_detID(int detID, bool isEB);

        //enums
	// OLD Categories without 6jet category
        /* enum RazorBox { //boxes for razor inclusive analysis */
	/*   MuEle = 0,  */
	/*   MuMu = 1, */
	/*   EleEle = 2, */
	/*   MuMultiJet = 3, */
	/*   MuJet = 4, */
	/*   EleMultiJet = 5, */
	/*   EleJet = 6, */
	/*   LooseLeptonMultiJet = 7, */
	/*   MultiJet = 8, */
	/*   LooseLeptonDiJet = 9, */
	/*   DiJet = 10, */
	/*   TwoBJet = 10, */
	/*   OneBJet = 11, */
	/*   ZeroBJet = 12, */
	/*   NONE = 999 */
        /* }; */
        enum RazorBox { //boxes for razor inclusive analysis
	  MuEle = 0,
	  MuMu = 1,
	  EleEle = 2,
	  MuSixJet = 3,
	  MuFourJet = 4,
	  MuJet = 5,
	  EleSixJet = 6,
	  EleFourJet = 7,
	  EleJet = 8,
	  LooseLeptonSixJet = 9,
	  LooseLeptonFourJet = 10,
	  SixJet = 11,
	  FourJet = 12,
	  LooseLeptonDiJet = 13,
	  DiJet = 14,
	  TwoBJet = 15,
	  OneBJet = 16,
	  ZeroBJet = 17,
	  MuMultiJet = 18,
	  EleMultiJet = 19,
	  LooseLeptonMultiJet = 20,
	  MultiJet = 21,
	  NONE = 999
        };
};

#endif
