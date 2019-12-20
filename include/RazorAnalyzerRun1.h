//Class for analyzing ntuples produced by the RazorTuplizer framework
//
//Author: Caltech Razor team

#ifndef RazorAnalyzerRun1_h
#define RazorAnalyzerRun1_h

#include "RazorEventsRun1.h" //This is a MakeClass of the RazorEventsRun1 tree in the ntuple to be analyzed
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



constexpr bool  EB_EleVeto           = true;
constexpr float EB_HoverECut         = 0.05;
constexpr float EB_SigmaIetaIetaCut  = 0.012;
constexpr float EB_PFchHadIsoCut     = 2.6;
constexpr float EB_PFnHadIsoConst    = 3.5;
constexpr float EB_PFnHadIsoSlope    = 0.04;
constexpr float EB_PFphoIsoConst     = 1.3;
constexpr float EB_PFphoIsoSlope     = 0.005;

constexpr bool  EE_EleVeto           = true;
constexpr float EE_HoverECut         = 0.05;
constexpr float EE_SigmaIetaIetaCut  = 0.034;
constexpr float EE_PFchHadIsoCut     = 2.3;
constexpr float EE_PFnHadIsoConst    = 2.9;
constexpr float EE_PFnHadIsoSlope    = 0.04;

enum class WP { VeryLoose, Loose, Medium, Tight };

float EB_GetPFchHadIsoCut( );
float EB_GetPFchHadIsoCut_VL( );
float EB_GetPFchHadIsoCut_M( );
float EB_GetPFchHadIsoCut_T( );

float EB_GetPFnHadIsoCut( float pt );
float EB_GetPFnHadIsoCut_VL( float pt );
float EB_GetPFnHadIsoCut_M( float pt );
float EB_GetPFnHadIsoCut_T( float pt );

float EB_GetPFphoIsoCut( float pt );
float EB_GetPFphoIsoCut_VL( float pt );
float EB_GetPFphoIsoCut_M( float pt );
float EB_GetPFphoIsoCut_T( float pt );

float EE_GetPFchHadIsoCut( );
float EE_GetPFchHadIsoCut_VL( );
float EE_GetPFchHadIsoCut_M( );
float EE_GetPFchHadIsoCut_T( );

float EE_GetPFnHadIsoCut( float pt );
float EE_GetPFnHadIsoCut_VL( float pt );
float EE_GetPFnHadIsoCut_M( float pt );
float EE_GetPFnHadIsoCut_T( float pt );

float EE_GetPFphoIsoCut_M( float pt );
float EE_GetPFphoIsoCut_T( float pt );



class RazorAnalyzerRun1: public RazorEventsRun1 {
    public :
        RazorAnalyzerRun1(TTree *tree=0);
        virtual ~RazorAnalyzerRun1();

        void EnableEventInfo();
//        void EnablePVAll();
        void EnableMuons();
        void EnableElectrons();
        void EnableTaus();
        void EnableIsoPFCandidates();
        void EnablePhotons();
        void EnableJets();
        void EnableFatJets();
        void EnableMet();
        void EnablePileup();
        void EnableMC();
        void EnableGenParticles();
        void EnableRazor();

        void EnableAll();

        //------ LIST OF ANALYSES ------//
        virtual void Analyze(bool isData, int option, string outputFileName, string label);

        //functions in RazorAuxMuon.cc
	float GetMuonEffectiveAreaMean(int i, string type );
	bool isMuonPOGLooseMuon(int i, bool applyID = true, bool applyIso = true);
//	bool isMuonPOGMediumMuon(int i, bool applyID = true, bool applyIso = true);
        bool isMuonPOGTightMuon(int i, bool applyID = true, bool applyIso = true);
//	bool isVetoMuon(int i, bool applyID = true, bool applyIso = true);
//	bool isLooseMuon(int i, bool applyID = true, bool applyIso = true);
//      bool isTightMuon(int i, bool applyID = true, bool applyIso = true);
        bool passHZZMuonPreselection(int i);
        bool isHZZMuon(int i);
//	bool matchMuonHLTFilters( int i, string HLTFilter);
//	bool matchTagMuonHLTFilters( int i);

        //functions in RazorAuxElectron.cc
	float GetElectronScaleCorrection( double pt, double eta );
	float GetElectronEffectiveAreaMean(int i, bool use25nsCuts = false);
	float GetElectronEffectiveArea90(int i);
        bool isEGammaPOGVetoElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = false);
        bool isEGammaPOGLooseElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = false);
        bool isEGammaPOGMediumElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = false);
        bool isEGammaPOGTightElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = false);
//      bool isVetoElectron(int i, bool applyID = true, bool applyIso = true);
//      bool isLooseElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = false);
//      bool isMediumElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = false);
//      bool isTightElectron(int i, bool applyID = true, bool applyIso = true, bool use25nsCuts = false);
//	bool isMVANonTrigVetoElectron(int i, bool applyID = true, bool applyIso = true);
        bool passEGammaPOGVetoElectronID(int i, bool use25nsCuts = false);
        bool passEGammaPOGLooseElectronID(int i, bool use25nsCuts = false);
        bool passEGammaPOGMediumElectronID(int i, bool use25nsCuts = false);
        bool passEGammaPOGTightElectronID(int i, bool use25nsCuts = false);
//	bool passMVANonTrigVetoElectronID(int i);
        bool passEGammaPOGVetoElectronIso(int i, bool use25nsCuts = false);
        bool passEGammaPOGLooseElectronIso(int i, bool use25nsCuts = false);
        bool passEGammaPOGMediumElectronIso(int i, bool use25nsCuts = false);
        bool passEGammaPOGTightElectronIso(int i, bool use25nsCuts = false);
//	bool passMVANonTrigVetoElectronIso(int i);
	bool passHZZElectronIso(int i);
	bool passHZZElectronPreselection(int i);
	bool isHZZElectron(int i);
//	bool matchElectronHLTFilters( int i, string HLTFilter, string analysisTag);
//	bool matchElectronHLTFilters2015( int i, string HLTFilter);
//	bool matchElectronHLTFilters2016( int i, string HLTFilter);
//	bool matchProbeElectronHLTFilters( int i);
//	bool matchProbeSCHLTFilters( int i);
//	bool matchTagElectronHLTFilters( int i);

        //functions in RazorAuxTau.cc
        bool isLooseTau(int i);
        bool isMediumTau(int i);
        bool isTightTau(int i);

        //functions in RazorAuxPhoton.cc
        bool photonPassesElectronVeto(int i);
	void getPhotonEffAreaRun2( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho );
	void getPhotonEffArea90( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho );
      bool photonPassesIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut, bool useEffectiveArea90);
	bool photonPassLooseIDWithoutEleVeto(int i, bool use25nsCuts = false);
	bool photonPassMediumIDWithoutEleVeto(int i, bool use25nsCuts = false);
	bool photonPassTightIDWithoutEleVeto(int i, bool use25nsCuts = false);
	bool photonPassLooseID(int i, bool use25nsCuts = false);
	bool photonPassMediumID(int i, bool use25nsCuts = false);
	bool photonPassTightID(int i, bool use25nsCuts = false);
	bool photonPassLooseIso(int i, bool use25nsCuts = false);
	bool photonPassMediumIso(int i, bool use25nsCuts = false);
	bool photonPassTightIso(int i, bool use25nsCuts = false);
        bool isLoosePhoton(int i, bool use25nsCuts = false);
        bool isMediumPhoton(int i, bool use25nsCuts = false);
        bool isTightPhoton(int i, bool use25nsCuts = false);
        bool isLoosePhotonWithoutEleVeto(int i, bool use25nsCuts = false);
        bool isMediumPhotonWithoutEleVeto(int i, bool use25nsCuts = false);
        bool isTightPhotonWithoutEleVeto(int i, bool use25nsCuts = false);
	bool matchPhotonHLTFilters( int i, string HLTFilter);
	void getPhotonEffAreaExo15004( float eta, double& effAreaPho );
	bool photonPassLooseIDWithoutEleVetoExo15004(int i);
//	bool photonPassesIsolationExo15004(int i, double PFChHadIsoCut, double PFPhotIsoCut );
//	bool photonPassLooseIsoExo15004(int i);
	TLorentzVector GetCorrectedMomentum( TVector3 vtx, TVector3 phoPos, double phoE );

	/* //function in HggRazorAuxPhoton.cc */
	/* // R u n 1   C u t   B a s e d   I D */
	/* //---------------------------------- */
	 bool isGoodPhotonRun1( int i, bool _iso, bool _debug); 
	 bool photonPassIsoRun1( int i, bool _debug =false); 
	 bool photonPassIsoRun1( int i , WP wp, bool _debug=false );
	 void getPhotonEffAreaRun1( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho ); 
	 bool passEleVetoRun1( int i ); 
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
					  FactorizedJetCorrector *jetcorrector,  
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
	
    double GetAlphaT(vector<TLorentzVector> jets) ;
    double GetDPhiMin(vector<TLorentzVector> jets);

        bool passesHadronicRazorBaseline(double MR, double Rsq);
        bool passesLeptonicRazorBaseline(double MR, double Rsq);
        int SubtractParticleFromCollection(TLorentzVector ToSubtract, vector<TLorentzVector>& Collection, float deltaRMatch=0.4);
	
	double calcMT2(float testMass, bool massive, std::vector<TLorentzVector> jets, TLorentzVector MET, int hemi_seed, int hemi_association);
	
	//functions in src/RazorAuxGenLevel.cc
	bool matchesGenMuon(double eta, double phi);
	bool matchesGenElectron(double eta, double phi);
	bool isGenTau(int index);
	bool isGenLeptonicTau(int index);
	int findClosestGenElectron(double eta, double phi);
	int findClosestGenMuon(double eta, double phi);
	int findClosestGenTau(double eta, double phi);
	int findClosestRecoTau(double eta, double phi);
	int GetTauMatchedID(double eta, double phi);
	int findClosestParton(float eta, float phi);

	//Added to src/RazorAuxGenLevel.cc
	int findClosestGenJet(double eta, double phi);
	
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
