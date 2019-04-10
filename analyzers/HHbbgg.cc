//LOCAL INCLUDES
#include "HHbbgg.h"
#include "RazorHelper.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
//C++ INCLUDES
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h> 
//ROOT INCLUDES
#include <TH1F.h>
#include <TH2D.h>
#include "TRandom3.h"
#include "TMVA/Reader.h"

using namespace std;


struct PhotonCandidate
{                                                  
  int   Index;
  TLorentzVector photon;
  TLorentzVector photonSC;
  float scE;
  float scPt;
  float scEta;
  float scPhi;
  float scX;
  float scY;
  float scZ;
  float SigmaIetaIeta;                                                                        
  float R9;                                                                                  
  float HoverE;                                                                        
  float sumChargedHadronPt;                                                                
  float sumNeutralHadronEt;                                                     
  float sumPhotonEt;                                            
  float sigmaEOverE;
  bool  _passEleVeto;
  bool  _passIso;
  int convType;
  float convTrkZ;
  float convTrkClusZ;
  float vtxSumPx[200];
  float vtxSumPy[200];
};

struct evt
{
  std::string run;
  std::string event;
};

#define _phodebug 0
#define _debug    0
#define _info     1

const double EB_R = 129.0;
const double EE_Z = 317.0;

const double JET_CUT = 30.;
const int NUM_PDF_WEIGHTS = 60;

const float SIGMATRKZ[12] = { 0.0125255,
                              0.0178107,
                              0.0581667,
                              0.152157,
                              0.716301,
                              1.3188,
                              0.38521,
                              0.702755,
                              3.17615,
                              2.23662,
                              1.67937,
                              2.46599
                            };
                            
const float SIGMATRKCLUSZ[12] = { 0.0298574,
                                  0.0935307,
                                  0.180419,
                                  0.577081,
                                  0.414393,
                                  0.756568,
                                  0.494722,
                                  0.892751,
                                  1.06805,
                                  0.62143,
                                  1.21941,
                                  1.56638
                                };

//Testing branching and merging
void HHbbgg::Analyze(bool isData, int option, string outFileName, string label)
{
  //*****************************************************************************
  //Settings
  //*****************************************************************************
  TRandom3 random(3003);
  bool doPhotonScaleCorrection = true;
  bool isFastsimSMS = false;
  bool doMVAVertex = true;

  string analysisTag = "Razor2016_80X";
  if ( label != "") analysisTag = label;

  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  if (analysisTag == "Razor2015_76X") helper = new RazorHelper("Razor2015_76X", isData, isFastsimSMS);
  else if (analysisTag == "Razor2016_80X") helper = new RazorHelper("Razor2016_80X", isData, isFastsimSMS);
  else helper = new RazorHelper(analysisTag, isData, isFastsimSMS);
  

  //initialization: create one TTree for each analysis box 
  if ( _info ) std::cout << "Initializing..." << std::endl;
  
  if ( outFileName.empty() )
    {
      if ( _info ) std::cout << "HHbbgg: Output filename not specified!" << endl << "Using default output name HHbbgg.root" << std::endl;
      outFileName = "HHbbgg.root";
    }
  TFile* outFile = new TFile( outFileName.c_str(), "RECREATE" );
  //---------------------------
  //one tree to hold all events
  //---------------------------
  TTree *razorTree = new TTree("HHbbgg", "Info on selected razor inclusive events");
  
  //Get CMSSW Directory
  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");


  //--------------------------------
  //Photon Energy Scale and Resolution Corrections
  //--------------------------------
  std::string photonCorrectionPath = "./";
  // if ( cmsswPath != NULL ) photonCorrectionPath = string(cmsswPath) + "/src/RazorAnalyzer/data/PhotonCorrections/";

  EnergyScaleCorrection_class *photonCorrector = 0;
  if (analysisTag == "Razor2015_76X") {
    photonCorrector = new EnergyScaleCorrection_class(Form("%s/76X_16DecRereco_2015", photonCorrectionPath.c_str()));
  } else if (analysisTag == "Razor2016_80X") {
    photonCorrector = new EnergyScaleCorrection_class(Form("%s/80X_2016", photonCorrectionPath.c_str()));
  } else if (analysisTag == "Razor2016_MoriondRereco") {
    photonCorrector = new EnergyScaleCorrection_class(Form("%s/Winter_2016_reReco_v1_ele", photonCorrectionPath.c_str()));
  }

  if(!isData) {
    photonCorrector->doScale = false; 
    photonCorrector->doSmearings = true;
  } else {
    photonCorrector->doScale = true; 
    photonCorrector->doSmearings = false;
  }


  //--------------------------------
  //MVA Vertex Selection
  //--------------------------------
  float ptasym = 0.;
  float ptbal = 0.;
  float logsumpt2 = 0.;
  float pull_conv = 0.;
  float nConv = 0.;
  
  
  TMVA::Reader *vtxmvareader = 0;
  if (doMVAVertex) {
    vtxmvareader = new TMVA::Reader( "!Color:Silent" );
    vtxmvareader->AddVariable("ptasym", &ptasym );
    vtxmvareader->AddVariable("ptbal", &ptbal );
    vtxmvareader->AddVariable("logsumpt2", &logsumpt2 );
    vtxmvareader->AddVariable("limPullToConv", &pull_conv );
    vtxmvareader->AddVariable("nConv", &nConv );
    
    std::string vtxpathname;
    if ( cmsswPath != NULL ) vtxpathname = string(cmsswPath) + "/src/RazorAnalyzer/data/";
    vtxmvareader->BookMVA("BDT",Form("%s/TMVAClassification_BDTVtxId_SL_2015.xml",vtxpathname.c_str()));
  }

  //--------------------------------
  //Including Jet Energy Corrections
  //--------------------------------
  std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
  std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();


  //----------
  //pu histo
  //----------
  TH1D* puhisto = new TH1D("pileup", "", 50, 0, 50);
  

  //histogram containing total number of processed events (for normalization)
  TH1F *histNISRJets = new TH1F("NISRJets", "NISRJets", 7, -0.5, 6.5);
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
  TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
  TH1F *SumPdfWeights = new TH1F("SumPdfWeights", "SumPdfWeights", NUM_PDF_WEIGHTS, -0.5, NUM_PDF_WEIGHTS-0.5);

  //--------------
  //tree variables
  //--------------
  float weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;  
  int nSelectedPhotons;
  float pho1Pt, pho1Eta, pho1Phi;
  float pho2Pt, pho2Eta, pho2Phi;
  float bjet1Pt, bjet1Eta, bjet1Phi, bjet1E;
  float bjet2Pt, bjet2Eta, bjet2Phi, bjet2E;
  float mGammaGamma, pTGammaGamma;  
  float mbb, mbbgg;
  float ybbgg, ptbbgg;
  int NPU;
  unsigned int run, lumi, event;

  //------------------------
  //set branches on big tree
  //------------------------

  razorTree->Branch("weight", &weight, "weight/F");
  razorTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
  razorTree->Branch("pileupWeightUp", &pileupWeightUp, "pileupWeightUp/F");
  razorTree->Branch("pileupWeightDown", &pileupWeightDown, "pileupWeightDown/F");
  razorTree->Branch("mGammaGamma", &mGammaGamma, "mGammaGamma/F");
  razorTree->Branch("pTGammaGamma", &pTGammaGamma, "pTGammaGamma/F");      
  razorTree->Branch("mbb", &mbb, "mbb/F");
  razorTree->Branch("mbbgg", &mbbgg, "mbbgg/F");
  razorTree->Branch("ybbgg", &ybbgg, "ybbgg/F");
  razorTree->Branch("ptbbgg", &ptbbgg, "ptbbgg/F");
  razorTree->Branch("NPU", &NPU, "npu/i");
  razorTree->Branch("run", &run, "run/i");
  razorTree->Branch("lumi", &lumi, "lumi/i");
  razorTree->Branch("event", &event, "event/i");

  razorTree->Branch("pho1Pt", &pho1Pt, "pho1Pt/F");
  razorTree->Branch("pho1Eta", &pho1Eta, "pho1Eta/F");
  razorTree->Branch("pho1Phi", &pho1Phi, "pho1Phi/F");
  razorTree->Branch("pho2Pt", &pho2Pt, "pho2Pt/F");
  razorTree->Branch("pho2Eta", &pho2Eta, "pho2Eta/F");
  razorTree->Branch("pho2Phi", &pho2Phi, "pho2Phi/F");
  razorTree->Branch("bjet1E", &bjet1E, "bjet1E/F");
  razorTree->Branch("bjet1Pt", &bjet1Pt, "bjet1Pt/F");
  razorTree->Branch("bjet1Eta", &bjet1Eta, "bjet1Eta/F");
  razorTree->Branch("bjet1Phi", &bjet1Phi, "bjet1Phi/F");
  razorTree->Branch("bjet2E", &bjet2E, "bjet2E/F");
  razorTree->Branch("bjet2Pt", &bjet2Pt, "bjet2Pt/F");
  razorTree->Branch("bjet2Eta", &bjet2Eta, "bjet2Eta/F");
  razorTree->Branch("bjet2Phi", &bjet2Phi, "bjet2Phi/F");


  //begin loop
  if ( fChain == 0 ) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  std::cout << "[INFO]: Total Entries = " << fChain->GetEntries() << "\n";
  for ( Long64_t jentry=0; jentry < nentries; jentry++ )
    {
      //begin event
      if( _info && (jentry % 10000 == 0) ) std::cout << "[INFO]: Processing entry " << jentry << std::endl;
      Long64_t ientry = LoadTree( jentry );
      if ( ientry < 0 ) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
    
      //fill normalization histogram    
      NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
      weight = genWeight;
      SumWeights->Fill(1.0, weight);
      
      //reset tree variables     
      pileupWeight      = 1.0;
      pileupWeightUp    = 1.0;
      pileupWeightDown  = 1.0;
      run = runNum;
      lumi = lumiNum; 
      event = eventNum;

      mGammaGamma    = -1;
      pTGammaGamma   = -1;
      mbb    = -1;
      mbbgg    = -1;

      pho1Pt = 0;
      pho1Eta = 0;
      pho1Phi = 0;
      pho2Pt = 0;
      pho2Eta = 0;
      pho2Phi = 0;
      bjet1Pt = 0;
      bjet1Eta = 0;
      bjet1Phi = 0;
      bjet2Pt = 0;
      bjet2Eta = 0;
      bjet2Phi = 0;

      //------------------
      //Pileup reweighting
      //------------------
      pileupWeight = 1.0;
      if( !isData ) {
	//Get number of PU interactions
	for (int i = 0; i < nBunchXing; i++) {
	  if (BunchXing[i] == 0) {
	    NPU = nPUmean[i];
	  }
	}
	puhisto->Fill(NPU);
	pileupWeight = helper->getPileupWeight(NPU);
	pileupWeightUp = helper->getPileupWeightUp(NPU) / pileupWeight;
	pileupWeightDown = helper->getPileupWeightDown(NPU) / pileupWeight;	
      }
      
      
      //photon selection
      vector<TLorentzVector> GoodPhotons;
      vector<double> GoodPhotonSigmaE; // energy uncertainties of selected photons
      vector<bool> GoodPhotonPassesIso; //store whether each photon is isolated
      std::vector< PhotonCandidate> phoCand;//PhotonCandidate defined in RazorAuxPhoton.hh
      int nPhotonsAbove40GeV = 0;
      for(int i = 0; i < nPhotons; i++) {
	  
	  double scale = photonCorrector->ScaleCorrection(run, (fabs(pho_superClusterEta[i]) < 1.5), phoR9[i], pho_superClusterEta[i], phoE[i]/cosh(pho_superClusterEta[i]));
	  double smear = photonCorrector->getSmearingSigma(run, (fabs(pho_superClusterEta[i]) < 1.5), phoR9[i], pho_superClusterEta[i], phoE[i]/cosh(pho_superClusterEta[i]), 0., 0.); 

	  //Defining Corrected Photon momentum
	  float pho_pt_corr = phoPt[i];
	  if (isData) {
	    pho_pt_corr = phoPt[i]*scale; 
	    if (_phodebug) std::cout << "[DEBUG] : Photon Energy Scale Corrections: " << phoPt[i] << " * " << scale << " --> " << pho_pt_corr << "\n";
	  } else {
	    pho_pt_corr = phoPt[i]*(1+smear*random.Gaus());
	  }
	  TVector3 vec;
	  vec.SetPtEtaPhi( pho_pt_corr, phoEta[i], phoPhi[i] );
	

	  //**********************************************************
	  //Isolation, electron veto, and Barrel requirements are introduced here 
	  //if we want to use the "regular" selection sequence
	  //**********************************************************
	  if ( !photonPassLooseIDWithoutEleVeto(i) ) continue;
	  if (!(pho_passEleVeto[i])) continue;
	  if (!(photonPassLooseIso(i))) continue;
	
	  if ( fabs(pho_superClusterEta[i]) > 1.4442 && fabs(pho_superClusterEta[i]) < 1.566 ) continue;
	  if ( fabs(pho_superClusterEta[i]) > 2.5 ) continue;
	  if ( phoPt[i] < 20.0 ) continue;
		
	  //setting up photon 4-momentum with zero mass
	  TLorentzVector thisPhoton;
	  thisPhoton.SetVectM( vec, .0 );

	  //-----------------------------
	  //uncorrected photon 4-momentum
	  //-----------------------------
	  TVector3 vtx( pvX, pvY, pvZ );
	  TVector3 phoPos;
	  if ( fabs( pho_superClusterEta[i] ) < 1.479 )
	    {
	      phoPos.SetXYZ( EB_R*cos( pho_superClusterPhi[i]), EB_R*sin( pho_superClusterPhi[i] ), EB_R*sinh( pho_superClusterEta[i] ) );
	    }
	  else
	    {
	      double R = fabs( EE_Z/sinh( pho_superClusterEta[i] ) );
	    
	      if ( pho_superClusterEta[i] > .0 )
		{
		  phoPos.SetXYZ( R*cos( pho_superClusterPhi[i] ), R*sin( pho_superClusterPhi[i] ), EE_Z);
		}
	      else
		{
		  phoPos.SetXYZ( R*cos( pho_superClusterPhi[i] ), R*sin( pho_superClusterPhi[i] ), -EE_Z);
		}
	    
	    }
	
	  TLorentzVector phoSC = GetCorrectedMomentum( vtx, phoPos, pho_RegressionE[i] );
	
	  //Filling Photon Candidate
	  PhotonCandidate tmp_phoCand;
	  tmp_phoCand.Index = i;
	  tmp_phoCand.photon = thisPhoton;
	  tmp_phoCand.photonSC = phoSC;
	  tmp_phoCand.scEta = pho_superClusterEta[i];
	  tmp_phoCand.scPhi = pho_superClusterPhi[i];
	  tmp_phoCand.SigmaIetaIeta = phoFull5x5SigmaIetaIeta[i];
	  tmp_phoCand.R9 = phoR9[i];
	  tmp_phoCand.HoverE = pho_HoverE[i];
	  tmp_phoCand.sumChargedHadronPt = pho_pfIsoChargedHadronIso[i];
	  tmp_phoCand.sumNeutralHadronEt = pho_pfIsoNeutralHadronIso[i];
	  tmp_phoCand.sumPhotonEt = pho_pfIsoPhotonIso[i];
	  tmp_phoCand.sigmaEOverE = pho_RegressionEUncertainty[i]/pho_RegressionE[i];
	  tmp_phoCand._passEleVeto = pho_passEleVeto[i];
	  tmp_phoCand._passIso = photonPassLooseIso(i);
	  phoCand.push_back( tmp_phoCand );
	
	  nSelectedPhotons++;
	}
    

      //--------------------------------------
      //Require at least two photon candidates
      //--------------------------------------
      if ( phoCand.size() < 2 ) {
      	if ( _debug ) std::cout << "[INFO]: not enough photon, nphotons: " 
      				<< phoCand.size() << std::endl;
      	for(int i = 0; i < nPhotons; i++) {
      	  if ( _debug ) std::cout << "pho# " << i << " phopt1: " << phoPt[i] 
      				  << " pho_eta: " << phoEta[i] 
      				  << " SIetaIeta: " << phoFull5x5SigmaIetaIeta[i] << std::endl;
      	}
      	continue;
      }
      
      
      if ( _debug ) std::cout << "[DEBUG]: nphotons--> " << phoCand.size() 
			      << " " << nSelectedPhotons << std::endl;
    
      //----------------------------------------
      //find the "best" photon pair, highest Pt!
      //----------------------------------------
      TLorentzVector HiggsCandidate(0,0,0,0);
      TLorentzVector HiggsCandidateSC(0,0,0,0);
      int HiggsPhoIndex1 = -1;
      int HiggsPhoIndex2 = -1;
      double bestSumPt = -99.;
      std::vector< PhotonCandidate > phoSelectedCand;
      PhotonCandidate bestCand[2];
      for ( size_t i = 0; i < phoCand.size(); i++ )
	{
	  for ( size_t j = i+1; j < phoCand.size(); j++ )
	    {
	      PhotonCandidate pho1 = phoCand[i];
	      PhotonCandidate pho2 = phoCand[j];


	      if (doMVAVertex) {

		float pho1E = pho1.photon.E();
		float pho2E = pho2.photon.E();
              
		float maxbdtval = -99.;
		int ipvmax = 0;
                            
		for (int ipv=0; ipv<nPVAll; ++ipv) {
		  float vtxX = pvAllX[ipv];
		  float vtxY = pvAllY[ipv];
		  float vtxZ = pvAllZ[ipv];
                
		  TVector3 pho1dir(pho1.scX-vtxX,pho1.scY-vtxY,pho1.scZ-vtxZ);
		  TVector3 pho2dir(pho2.scX-vtxX,pho2.scY-vtxY,pho2.scZ-vtxZ);
		  TVector3 diphomom = pho1E*pho1dir.Unit() + pho2E*pho2dir.Unit();
		  TVector3 diphoPt(diphomom.x(), diphomom.y(), 0.);
                
		  TVector3 vtxSumPt(pvAllSumPx[ipv]-pho1.vtxSumPx[ipv]-pho2.vtxSumPx[ipv], pvAllSumPy[ipv]-pho1.vtxSumPy[ipv]-pho2.vtxSumPy[ipv],0.);
                

		  if (_debug) {
		    std::cout << " Pho Debug: " << pho1E << " " << pho2E << " : " 
			      << diphoPt.Eta() << " " << diphoPt.Phi() << " "
			      << "\n";
		  }

		  logsumpt2 = pvAllLogSumPtSq[ipv];
		  ptbal = -vtxSumPt.Dot(diphoPt.Unit()); 
		  ptasym = (vtxSumPt.Mag() - diphoPt.Mag())/(vtxSumPt.Mag() + diphoPt.Mag());
                
		  bool hasconv1 = pho1.convType>=0;
		  bool hasconv2 = pho2.convType>=0;
                
		  float convz1 = -99.;
		  float convsz1 = -99.;
		  float convz2 = -99.;
		  float convsz2 = -99.;
                                
		  nConv = 0.;
		  float convz = -99.;
		  float convsz = -99.;
		  if (hasconv1) {
		    convz1 = SIGMATRKZ[pho1.convType] < SIGMATRKCLUSZ[pho1.convType] ? pho1.convTrkZ : pho1.convTrkClusZ;
		    convsz1 = std::min(SIGMATRKZ[pho1.convType],SIGMATRKCLUSZ[pho1.convType]);                  
                  
		    convz = convz1;
		    convsz = convsz1;
		    nConv = 1;
		  }
		  if (hasconv2) {
		    convz2 = SIGMATRKZ[pho2.convType] < SIGMATRKCLUSZ[pho2.convType] ? pho2.convTrkZ : pho2.convTrkClusZ;
		    convsz2 = std::min(SIGMATRKZ[pho2.convType],SIGMATRKCLUSZ[pho2.convType]);
                  
		    convz = convz2;
		    convsz = convsz2;
		    nConv = 1;
		  }
		  if (hasconv1 && hasconv2) {
		    double w1 = 1./convsz1/convsz1;
		    double w2 = 1./convsz2/convsz2;
		    convz = (w1*convz1 + w2*convz2)/(w1+w2);
		    convsz = sqrt(1./(w1+w2));
		    nConv = 2;
		  }
                
		  if (hasconv1 || hasconv2) {
		    pull_conv = std::abs(convz-vtxZ)/convsz;
		    pull_conv = std::min(10.f,pull_conv); 
		  }
		  else {
		    pull_conv = 10.;
		  }
                
		  float bdtval = vtxmvareader->EvaluateMVA("BDT");

		  if (_debug) std::cout << "Vertex: " << ipv << " " << pvAllZ[ipv] << " : " 
					<< ptasym << " " << ptbal << " " << logsumpt2 << " " << pull_conv << " " << nConv << " : "
					<< bdtval << "\n";
                              
		  if (bdtval > maxbdtval) {
		    maxbdtval = bdtval;
		    ipvmax = ipv;
		  }
                
		}

		if (_debug) std::cout << "Vertex Selected: " << ipvmax << "\n\n";
	      
		//vtxIndex = ipvmax;
		//vtxZCoordinate =  pvAllZ[ipvmax];

		float vtxX = pvAllX[ipvmax];
		float vtxY = pvAllY[ipvmax];
		float vtxZ = pvAllZ[ipvmax];
              
		TVector3 pho1dir(pho1.scX-vtxX,pho1.scY-vtxY,pho1.scZ-vtxZ);
		pho1dir = pho1dir.Unit();
		pho1.photon.SetPxPyPzE(pho1dir.x()*pho1E,pho1dir.y()*pho1E,pho1dir.z()*pho1E,pho1E);
              
		TVector3 pho2dir(pho2.scX-vtxX,pho2.scY-vtxY,pho2.scZ-vtxZ);
		pho2dir = pho2dir.Unit();
		pho2.photon.SetPxPyPzE(pho2dir.x()*pho2E,pho2dir.y()*pho2E,pho2dir.z()*pho2E,pho2E);
              
	      } //end if do MVA Vertex
            

	      if ( _debug )
		{
		  std::cout << "[DEBUG]: pho1-> " << pho1.photon.Pt()
			    << "\n[DEBUG]: pho2->" << pho2.photon.Pt() 
			    << std::endl;
		}
	      //need one photon in the pair to have pt > 40 GeV
	      if ( pho1.photon.Pt() < 40.0 && pho2.photon.Pt() < 40.0 )
		{
		  if ( _debug ) std::cout << "[DEBUG]: both photons failed PT > 40 GeV" << std::endl; 
		  //continue;
		}
	      //need diphoton mass between > 100 GeV as in AN (April 1st)
	      double diphotonMass = (pho1.photon + pho2.photon).M();
	      if ( _debug )
		{
		  std::cout << "[DEBUG] Diphoton Sum pT: " << pho1.photon.Pt() + pho2.photon.Pt() << std::endl;
		}
	    
	      if( diphotonMass < 100 )
		{
		  if ( _debug ) std::cout << "[DEBUG]: Diphoton mass < 50 GeV: mgg-> " << diphotonMass << std::endl;
		  if ( _debug ) std::cout << "... pho1Pt: " << pho1.photon.Pt()  << " pho2Pt: " << pho2.photon.Pt()  << std::endl;
		  continue;
		}
	      //---------------------------------------------
	      //if the sum of the photon pT's is larger than 
	      //that of the current Higgs candidate, 
	      //make this the Higgs candidate
	      //---------------------------------------------
	      if( pho1.photon.Pt() + pho2.photon.Pt() > bestSumPt )
		{
		  bestSumPt = pho1.photon.Pt() + pho2.photon.Pt();
		  HiggsCandidate = pho1.photon + pho2.photon;
		  HiggsCandidateSC = pho1.photonSC + pho2.photonSC;
		  if ( pho1.photon.Pt() >= pho2.photon.Pt() )
		    {
		      if ( _debug ) std::cout << "assign photon candidate, pho1Pt > pho2Pt" << std::endl;
		      bestCand[0] = pho1;
		      bestCand[1] = pho2;
		      HiggsPhoIndex1 = pho1.Index;
		      HiggsPhoIndex2 = pho2.Index;  
		    }
		  else
		    {
		      if ( _debug ) std::cout << "assign photon candidate, pho2Pt > pho1Pt" << std::endl;
		      bestCand[0] = pho2;
		      bestCand[1] = pho1;
		      HiggsPhoIndex1 = pho2.Index;
		      HiggsPhoIndex2 = pho1.Index;
		    }
		}//best pt if
	    } //second photon loop
	}//first photon loop
    
    
      //---------------------------------------
      //just use this container for convenience
      //to parse the data into TTree
      //---------------------------------------
      phoSelectedCand.push_back(bestCand[0]);
      phoSelectedCand.push_back(bestCand[1]);
    
      //-----------------------------------
      //Filling Selected Photon Information
      //-----------------------------------
      TLorentzVector pho_cand_vec[2];
      int _pho_index = 0;
      for ( auto& tmpPho : phoSelectedCand )
	{
	  if ( !( tmpPho.Index == HiggsPhoIndex1 || tmpPho.Index == HiggsPhoIndex2 ) ) continue;
	  if( _pho_index > 1 ) std::cerr << "[ERROR]: Photon index larger than 1!" << std::endl;

	  if (_pho_index ==0) {
	    pho1Pt = tmpPho.photon.Pt();
	    pho1Eta = tmpPho.photon.Eta();
	    pho1Phi = tmpPho.photon.Phi();
	  } else {
	    pho2Pt = tmpPho.photon.Pt();
	    pho2Eta = tmpPho.photon.Eta();
	    pho2Phi = tmpPho.photon.Phi();
	  }

	  pho_cand_vec[_pho_index]           = tmpPho.photon;
	  _pho_index++;
	}
    
      //removing events with less than two good photon candidates
      if ( _pho_index < 2 ) continue;
    
     

      //record higgs candidate info
      mGammaGamma    = HiggsCandidate.M();
      pTGammaGamma   = HiggsCandidate.Pt();



      //Find BJets
      //Choose two bjets with highest CSV score as the higgs candidate
      int nbjets = 0;
      double bjet1CSV = -999;
      double bjet2CSV = -999;
      TLorentzVector bjet1;
      TLorentzVector bjet2;
      for(int i = 0; i < nJets; i++) {
	  //Jet Corrections                                                                      
	  double JEC = JetEnergyCorrectionFactor( jetPt[i], jetEta[i], jetPhi[i], jetE[i],
						  fixedGridRhoAll, jetJetArea[i], runNum,
						  JetCorrectorIOV, JetCorrector );
      
	  TLorentzVector thisJet = makeTLorentzVector( jetPt[i]*JEC, jetEta[i], jetPhi[i], jetE[i]*JEC );
	
	  if( thisJet.Pt() < 30 ) continue;//According to the April 1st 2015 AN
	  if( fabs( thisJet.Eta() ) >= 2.4 ) continue;
	  if ( !jetPassIDLoose[i] ) continue;
	
	  //exclude selected photons from the jet collection
	  double deltaRJetPhoton = min( thisJet.DeltaR( pho_cand_vec[0] ), thisJet.DeltaR( pho_cand_vec[1] ) );
	  if ( deltaRJetPhoton <= 0.5 ) continue;//According to the April 1st 2015 AN
      
	  if (!isCSVM(i)) continue;

	  nbjets++;

	  if (jetCISV[i] > bjet1CSV) {
	    bjet2CSV = bjet1CSV;
	    bjet2Pt = bjet1Pt;
	    bjet2Eta = bjet1Eta;
	    bjet2Phi = bjet1Phi;
	    bjet2E = bjet1E;
	    bjet2 = bjet1;
	    bjet1CSV = jetCISV[i];
	    bjet1Pt = thisJet.Pt();
	    bjet1Eta = thisJet.Eta();
	    bjet1Phi = thisJet.Phi();
	    bjet1E = thisJet.E();
	    bjet1 = thisJet;
	  } else if (jetCISV[i] > bjet2CSV) {
	    bjet2CSV = jetCISV[i];
	    bjet2Pt = thisJet.Pt();
	    bjet2Eta = thisJet.Eta();
	    bjet2Phi = thisJet.Phi();
	    bjet2E = thisJet.E();
	    bjet2 = thisJet;
	  }
      } //loop over jets

      
      if (nbjets >= 2) {
	mbb = (bjet1+bjet2).M();
	mbbgg = (bjet1+bjet2+HiggsCandidate).M();
	ybbgg = (bjet1+bjet2+HiggsCandidate).Rapidity();
	ptbbgg = (bjet1+bjet2+HiggsCandidate).Pt();
      }
   
      razorTree->Fill();
    

      //end of event loop
    }
  
  if ( _info ) std::cout << "[INFO]: Number of events processed: " << NEvents->Integral() << std::endl;

  if ( _info ) std::cout << "[INFO]: Writing output trees..." << std::endl;    
  outFile->cd();
  razorTree->Write();
  NEvents->Write();
  SumWeights->Write();
  SumScaleWeights->Write();
  SumPdfWeights->Write();
  histNISRJets->Write();
  puhisto->Write();
  
  outFile->Close();
  delete photonCorrector;
  delete helper;

}
