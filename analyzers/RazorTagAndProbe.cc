
#include "RazorTagAndProbe.h"
#include "JetCorrectorParameters.h"
#include "TagAndProbePair.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

struct greater_than_pt{
    inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){
        return p1.Pt() > p2.Pt();
    }
};

void RazorTagAndProbe::Analyze(bool isData, int option, string outputfilename, string label)
{
    //initialization: create one TTree for each analysis box
    cout << "Initializing..." << endl;
    cout << "IsData = " << isData << "\n";
    cout << "Option = " << option << "\n";

    //---------------------------
    if( isData )
    {
      std::cout << "[INFO]: running on data with label: " << label << " and option: " << option << std::endl;
    }
    else
    {
      std::cout << "[INFO]: running on MC with label: " << label << " and option: " << option << std::endl;
    }

    //***************************************************
    //ten thousands digit refers to object type.
    //***************************************************
    //1: electrons
    //2: muons
    //3: photons
    //4: taus
    //***************************************************
    //thousands and hundreds digit refers to denominator
    //***************************************************
    //0: Basic Object (Track or Supercluster)
    //1: reco object
    //2: pass ID+Iso veto
    //3: pass ID+Iso loose
    //4: pass ID+Iso medium
    //5: pass ID+Iso tight
    //6: pass ID veto
    //7: pass ID loose
    //8: pass ID medium
    //9: pass ID tight
    //10: pass Iso veto
    //11: pass Iso loose
    //12: pass Iso medium
    //13: pass Iso tight
    //14: pass delayed photon ID + Iso loose
    //15: pass delayed photon ID + Iso medium
    //16: pass delayed photon ID + Iso tight

    //***************************************************
    //ones and tens digit refers to numerator cuts
    //***************************************************
    //1: reco object
    //2: pass ID+Iso veto
    //3: pass ID+Iso loose
    //4: pass ID+Iso medium
    //5: pass ID+Iso tight
    //6: pass ID veto
    //7: pass ID loose
    //8: pass ID medium
    //9: pass ID tight
    //10: pass Iso veto
    //11: pass Iso loose
    //12: pass Iso medium
    //13: pass Iso tight
    //14: pass delayed photon ID + Iso loose
    //15: pass delayed photon ID + Iso medium
    //16: pass delayed photon ID + Iso tight
    //50 - 99: pass HLT Filter ( see specific mapping in the code below )

    //Define 25ns for photons
    bool _is25ns = true;

    int objectTypeOption = floor(float(option) / 10000);
    int denominatorType = floor( float(option - objectTypeOption*10000) / 100);
    int numeratorType = option - objectTypeOption*10000 - denominatorType*100;
    if (objectTypeOption == 1) cout << "Object Type : Electrons\n";
    else if (objectTypeOption == 2) cout << "Object Type : Muons\n";
    else if (objectTypeOption == 3) cout << "Object Type : Photons\n";
    else if (objectTypeOption == 4) cout << "Object Type : Taus\n";

    if (denominatorType == 0) cout << "Denominator Type : Basic Object (Track or Supercluster)\n";
    else if (denominatorType == 1) cout << "Denominator Type : Reco Object\n";
    else if (denominatorType == 2) cout << "Denominator Type : ID+Iso Veto\n";
    else if (denominatorType == 3) cout << "Denominator Type : ID+Iso Loose\n";
    else if (denominatorType == 4) cout << "Denominator Type : ID+Iso Medium\n";
    else if (denominatorType == 5) cout << "Denominator Type : ID+Iso Tight\n";
    else if (denominatorType == 6) cout << "Denominator Type : ID Veto\n";
    else if (denominatorType == 7) cout << "Denominator Type : ID Loose\n";
    else if (denominatorType == 8) cout << "Denominator Type : ID Medium\n";
    else if (denominatorType == 9) cout << "Denominator Type : ID Tight\n";
    else if (denominatorType == 10) cout << "Denominator Type : Iso Veto\n";
    else if (denominatorType == 11) cout << "Denominator Type : Iso Loose\n";
    else if (denominatorType == 12) cout << "Denominator Type : Iso Medium\n";
    else if (denominatorType == 13) cout << "Denominator Type : Iso Tight\n";
    else if (denominatorType == 14) cout << "Denominator Type : delayed photon ID + Iso Loose\n";
    else if (denominatorType == 15) cout << "Denominator Type : delayed photon ID + Iso Medium\n";
    else if (denominatorType == 16) cout << "Denominator Type : delayed photon ID + Iso Tight\n";
    else if (denominatorType == 23) cout << "Denominator Type : EGamma ID+Iso Loose\n";
    else if (denominatorType == 24) cout << "Denominator Type : EGamma ID+Iso Medium\n";
    else if (denominatorType == 25) cout << "Denominator Type : EGamma ID+Iso Tight\n";
    else if (denominatorType >= 50) cout << "Denominator Type : pass HLT Filters\n";

    if (numeratorType == 1) cout << "Numerator Type : Reco Object\n";
    else if (numeratorType == 2) cout << "Numerator Type : ID+Iso Veto\n";
    else if (numeratorType == 3) cout << "Numerator Type : ID+Iso Loose\n";
    else if (numeratorType == 4) cout << "Numerator Type : ID+Iso Medium\n";
    else if (numeratorType == 5) cout << "Numerator Type : ID+Iso Tight\n";
    else if (numeratorType == 6) cout << "Numerator Type : ID Veto\n";
    else if (numeratorType == 7) cout << "Numerator Type : ID Loose\n";
    else if (numeratorType == 8) cout << "Numerator Type : ID Medium\n";
    else if (numeratorType == 9) cout << "Numerator Type : ID Tight\n";
    else if (numeratorType == 10) cout << "Numerator Type : Iso Veto\n";
    else if (numeratorType == 11) cout << "Numerator Type : Iso Loose\n";
    else if (numeratorType == 12) cout << "Numerator Type : Iso Medium\n";
    else if (numeratorType == 13) cout << "Numerator Type : Iso Tight\n";
    else if (numeratorType == 14) cout << "Numerator Type : delayed photon ID + Iso Loose\n";
    else if (numeratorType == 15) cout << "Numerator Type : delayed photon ID + Iso Medium\n";
    else if (numeratorType == 16) cout << "Numerator Type : delayed photon ID + Iso Tight\n";
    else if (numeratorType >= 50) cout << "Numerator Type : pass HLT Filters\n";

    Float_t ELE_MASS = 0.000511;
    Float_t MU_MASS  = 0.105658;

    //TRandom3 *random = new TRandom3(33333); //Artur wants this number 33333

    bool printSyncDebug = false;

    //*************************************************************************
    //Set up Output File
    //*************************************************************************
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "RazorControlRegions.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    TagAndProbePair *TPPair = new TagAndProbePair;
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
      if(jentry % 100000 == 0) cout << "Processing entry " << jentry << endl;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      printSyncDebug = false;
      if (printSyncDebug) {
	cout << "\n****************************************************************\n";
	cout << "Debug Event : " << runNum << " " << lumiNum << " " << eventNum << "\n";
      }

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
      TPPair->run = runNum;
      TPPair->lumi = lumiNum;
      TPPair->event = eventNum;

      //get NPU
      for (int i=0; i < nBunchXing; ++i) {
	if (BunchXing[i] == 0) {
	  TPPair->NPU_0 = nPUmean[i];
	}
	if (BunchXing[i] == -1) {
	  TPPair->NPU_Minus1 = nPUmean[i];
	}
	if (BunchXing[i] == 1) {
	  TPPair->NPU_Plus1 = nPUmean[i];
	}
      }
      TPPair->NPV = nPV;
      TPPair->Rho = fixedGridRhoFastjetAll;
      TPPair->met = metType1Pt;

      //******************************************
      //Find Generated leptons
      //******************************************
      vector<int> genLeptonIndex;
      //find gen electrons
      if (objectTypeOption == 1 || objectTypeOption == 3 ) {
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
      if (objectTypeOption == 2) {
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
      if (objectTypeOption == 1) {

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
	  if ( !isEGammaPOGTightElectron(indexTag) ) continue;

	  //Tag must match single electron HLT Filters OR tag leg filter of the dedicated T&P trigger
	  if ( !matchTagElectronHLTFilters(indexTag)) continue;

	  TLorentzVector vtag;
	  vtag.SetPtEtaPhiM(elePt[indexTag], eleEta[indexTag], elePhi[indexTag], ELE_MASS);

	  TPPair->passTighterTag = false;
	  if (
	      ( (ele_chargedIso[indexTag] + fmax(0.0,  ele_photonIso[indexTag] + ele_neutralHadIso[indexTag] - GetElectronEffectiveAreaMean(indexTag)*fixedGridRhoFastjetAll)) / elePt[indexTag] < 0.02)
	      && fabs(ele_d0[indexTag]) < 0.01
	      ) {
	    TPPair->passTighterTag = true;
	  }

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
	      bool genmatch = false;
	      for (int q=0;q<int(genLeptonIndex.size()); q++) {
		if ( deltaR(eleEta[indexProbe],elePhi[indexProbe],
			    gParticleEta[genLeptonIndex[q]],
			    gParticlePhi[genLeptonIndex[q]]) < 0.1) {
		  genmatch = true;
		}
	      }
	      if (!genmatch) continue;
	    }

	    // //Probe must match probe leg filter of the dedicated T&P trigger
	    // if ( !matchProbeElectronHLTFilters(indexTag)) continue;
	    // if ( !matchProbeSCHLTFilters(indexTag)) continue;

	    //cout << "Probe: " << elePt[indexProbe] << " " << eleEta[indexProbe] << " " <<  elePhi[indexProbe] << " : " << isTightElectron(indexProbe) << "\n";

	    //*******************************************************
	    //denominator selection
	    //*******************************************************
	    if (denominatorType == 1) {
	      // reco object doesn't require any additional cuts
	    }
	    if (denominatorType == 2) {
	      if ( !isVetoElectron(indexProbe) ) continue;
	    }
	    if (denominatorType == 3) {
	      if ( !isLooseElectron(indexProbe) ) continue;
	    }
	    if (denominatorType == 5) {
	      if ( !isTightElectron(indexProbe) ) continue;
	    }
	    if (denominatorType == 6) {
	      if ( !isVetoElectron(indexProbe, true, false) ) continue;
	    }
	    if (denominatorType == 7) {
	      if ( !isLooseElectron(indexProbe, true, false) ) continue;
	    }
	    if (denominatorType == 9) {
	      if ( !isTightElectron(indexProbe,true, false) ) continue;
	    }
	    if (denominatorType == 10) {
	      if ( !isVetoElectron(indexProbe, false, true) ) continue;
	    }

	    TLorentzVector vprobe;
	    vprobe.SetPtEtaPhiM(elePt[indexProbe], eleEta[indexProbe], elePhi[indexProbe], ELE_MASS);

	    TPPair->mass = (vtag+vprobe).M();
	    TPPair->pt = elePt[indexProbe];
	    TPPair->eta = eleEta[indexProbe];
	    TPPair->phi = elePhi[indexProbe];
	    TPPair->charge = eleCharge[indexProbe];
	    TPPair->Activity = ele_activityMiniIsoAnnulus[indexProbe];

	    //****************************************
	    //PASS OR FAIL
	    //****************************************
	    bool pass = false;
	    if (numeratorType == 2) {
	      pass = isVetoElectron(indexProbe);
	    }
	    if (numeratorType == 3) {
	      pass = isLooseElectron(indexProbe);
	    }
	    if (numeratorType == 5) {
	      pass = isTightElectron(indexProbe);
	    }
	    if (numeratorType == 6) {
	      pass = isVetoElectron(indexProbe, true, false);
	    }
	    if (numeratorType == 7) {
	      pass = isLooseElectron(indexProbe, true, false);
	    }
	    if (numeratorType == 9) {
	      pass = isTightElectron(indexProbe, true, false);
	    }
	    if (numeratorType == 10) {
	      pass = isVetoElectron(indexProbe, false, true);
	    }
	    if (numeratorType == 11) {
	      pass = isLooseElectron(indexProbe, false, true);
	    }
	    if (numeratorType == 13) {
	      pass = isTightElectron(indexProbe, false, true);
	    }
	    if (numeratorType == 23) {
	      pass = isEGammaPOGLooseElectron(indexProbe);
	    }
	    if (numeratorType == 25) {
	      pass = isEGammaPOGTightElectron(indexProbe);
	    }
	    if (numeratorType == 50) {
	      pass = matchElectronHLTFilters(indexProbe, "SingleElectron", "2016");
	    }
	    if (numeratorType == 51) {
	      pass = matchElectronHLTFilters(indexProbe, "Ele23Loose", "2016");
	    }
	    if (numeratorType == 52) {
	      pass = matchElectronHLTFilters(indexProbe, "Ele27Loose", "2016");
	    }
	    if (numeratorType == 53) {
	      pass = matchElectronHLTFilters(indexProbe, "Ele27Tight", "2016");
	    }
	    if (numeratorType == 54) {
	      pass = matchElectronHLTFilters(indexProbe, "Ele32Tight", "2016");
	    }
	    if (numeratorType == 60) {
	      pass = bool( matchElectronHLTFilters(indexProbe, "Ele23Loose", "2016") ||
			   matchElectronHLTFilters(indexProbe, "Ele27Loose", "2016") ||
			   matchElectronHLTFilters(indexProbe, "Ele27Tight", "2016") ||
			   matchElectronHLTFilters(indexProbe, "Ele32Tight", "2016"));
	    }
	    if (numeratorType == 61) {
	      pass = bool( matchElectronHLTFilters(indexProbe, "Ele23Loose", "2016") ||
			   matchElectronHLTFilters(indexProbe, "Ele27Loose", "2016") ||
			   matchElectronHLTFilters(indexProbe, "Ele27Tight", "2016") ||
			   matchElectronHLTFilters(indexProbe, "Ele32Loose", "2016") ||
			   matchElectronHLTFilters(indexProbe, "Ele32Tight", "2016") ||
			   matchElectronHLTFilters(indexProbe, "Ele105", "2016") ||
			   matchElectronHLTFilters(indexProbe, "Ele115", "2016")
			   );
	    }
	    TPPair->pass = pass;
	    // cout << " TP Pass: " << TPPair->pass << "\n";

	    //****************************************
	    //Fill Output Tree
	    //****************************************
	    TPPair->tree_->Fill();

	  } //loop over probe electrons

	} // loop over tag electrons

      } //if objects are electrons

      //*********************************************************
      //Muons
      //*********************************************************
      if (objectTypeOption == 2) {

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
	  if ( !matchTagMuonHLTFilters(indexTag)) continue;

	  TLorentzVector vtag;
	  vtag.SetPtEtaPhiM(muonPt[indexTag], muonEta[indexTag], muonPhi[indexTag], MU_MASS);

	  TPPair->passTighterTag = false;
	  if (
	      ( (muon_chargedIso[indexTag] + fmax(0.0,  muon_photonIso[indexTag] + muon_neutralHadIso[indexTag] - 0.5*muon_pileupIso[indexTag])) / muonPt[indexTag] < 0.07)
	      && fabs(muon_ip3dSignificance[indexTag]) < 3
	      ) {
	    TPPair->passTighterTag = true;
	  }


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
	      bool genmatch = false;
	      for (int q=0;q<int(genLeptonIndex.size()); q++) {
		if ( deltaR(muonEta[indexProbe],muonPhi[indexProbe],
			    gParticleEta[genLeptonIndex[q]],
			    gParticlePhi[genLeptonIndex[q]]) < 0.1) {
		  genmatch = true;
		}
	      }
	      if (!genmatch) continue;
	    }

	    //*******************************************************
	    //denominator selection
	    //*******************************************************
	    if (denominatorType == 1) {
	      // reco object doesn't require any additional cuts
	    }
	    if (denominatorType == 2) {
	      if ( !isVetoMuon(indexProbe) ) continue;
	    }
	    if (denominatorType == 3) {
	      if ( !isLooseMuon(indexProbe) ) continue;
	    }
	    if (denominatorType == 5) {
	      if ( !isTightMuon(indexProbe) ) continue;
	    }
	    if (denominatorType == 6) {
	      if ( !isVetoMuon(indexProbe, true, false) ) continue;
	    }
	    if (denominatorType == 7) {
	      if ( !isLooseMuon(indexProbe, true, false) ) continue;
	    }
	    if (denominatorType == 9) {
	      if ( !isTightMuon(indexProbe, true, false) ) continue;
	    }
	    if (denominatorType == 10) {
	      if ( !isVetoMuon(indexProbe, false, true) ) continue;
	    }

	    TLorentzVector vprobe;
	    vprobe.SetPtEtaPhiM(muonPt[indexProbe], muonEta[indexProbe], muonPhi[indexProbe], MU_MASS);

	    TPPair->mass = (vtag+vprobe).M();
	    TPPair->pt = muonPt[indexProbe];
	    TPPair->eta = muonEta[indexProbe];
	    TPPair->phi = muonPhi[indexProbe];
	    TPPair->charge = muonCharge[indexProbe];
	    TPPair->Activity = muon_activityMiniIsoAnnulus[indexProbe];

	    //****************************************
	    //PASS OR FAIL
	    //****************************************
	    bool pass = false;
	    if (numeratorType == 2) {
	      pass = isVetoMuon(indexProbe);
	    }
	    if (numeratorType == 3) {
	      pass = isLooseMuon(indexProbe);
	    }
	    if (numeratorType == 5) {
	      pass = isTightMuon(indexProbe);
	    }
	    if (numeratorType == 6) {
	      pass = isVetoMuon(indexProbe, true, false);
	    }
	    if (numeratorType == 7) {
	      pass = isLooseMuon(indexProbe, true, false);
	    }
	    if (numeratorType == 9) {
	      pass = isTightMuon(indexProbe, true, false);
	    }
	    if (numeratorType == 10) {
	      pass = isVetoMuon(indexProbe, false, true);
	    }
	    if (numeratorType == 11) {
	      pass = isLooseMuon(indexProbe, false, true);
	    }
	    if (numeratorType == 13) {
	      pass = isTightMuon(indexProbe, false, true);
	    }
	    if (numeratorType == 50) {
	      pass = matchMuonHLTFilters(indexProbe, "SingleMuon");
	    }
	    if (numeratorType == 51) {
	      pass = bool(matchMuonHLTFilters(indexProbe, "IsoMu20") || matchMuonHLTFilters(indexProbe, "IsoTkMu20"));
	    }
	    if (numeratorType == 52) {
	      pass = bool(matchMuonHLTFilters(indexProbe, "IsoMu27") || matchMuonHLTFilters(indexProbe, "IsoTkMu27"));
	    }
	    if (numeratorType == 53) {
	      pass = matchMuonHLTFilters(indexProbe, "Mu50");
	    }
	    if (numeratorType == 60) {
	      pass = bool(matchMuonHLTFilters(indexProbe, "IsoMu20")
			  || matchMuonHLTFilters(indexProbe, "IsoTkMu20")
			  || matchMuonHLTFilters(indexProbe, "IsoMu22")
			  || matchMuonHLTFilters(indexProbe, "IsoTkMu22")
			  || matchMuonHLTFilters(indexProbe, "IsoMu24")
			  || matchMuonHLTFilters(indexProbe, "IsoTkMu24")
			  || matchMuonHLTFilters(indexProbe, "IsoMu27")
			  || matchMuonHLTFilters(indexProbe, "IsoTkMu27")
			  || matchMuonHLTFilters(indexProbe, "Mu50")
			  || matchMuonHLTFilters(indexProbe, "TkMu50")
			  );
	    }
	    if (numeratorType == 61) {
	      pass = bool(matchMuonHLTFilters(indexProbe, "IsoMu20")
			  || matchMuonHLTFilters(indexProbe, "IsoTkMu20")
			  || matchMuonHLTFilters(indexProbe, "Mu50"));
	    }
	    TPPair->pass = pass;
	    //cout << " TP Pass: " << TPPair->pass << "\n";

	    //****************************************
	    //Fill Output Tree
	    //****************************************
	    TPPair->tree_->Fill();

	  } //loop over probe muons

	} // loop over tag muons

      } //if objects are muons

      //*********************************************************
      //Photons
      //*********************************************************
      if (objectTypeOption == 3) {
	//std::cout << "event: \n";
	//*******************************************************
	//Loop over Tag electrons
	//*******************************************************
	for(int indexTag = 0; indexTag < nElectrons; indexTag++){
	  //std::cout << "Tag : " << elePt[indexTag] << " " << eleEta[indexTag] << "\n";
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
	  if ( !isEGammaPOGTightElectron(indexTag) ) continue;

	  //Tag must match single electron HLT Filters OR tag leg filter of the dedicated T&P trigger
	  if ( !matchTagElectronHLTFilters(indexTag)) continue;

	  TLorentzVector vtag;
	  vtag.SetPtEtaPhiM(elePt[indexTag], eleEta[indexTag], elePhi[indexTag], ELE_MASS);

	  //*******************************************************
	  //Loop over Probe electrons with photon ID
	  //*******************************************************

	  for(int indexProbe = 0; indexProbe < nPhotons; indexProbe++){

	    if(phoPt[indexProbe] < 5) continue;
	    if(fabs(phoEta[indexProbe]) > 2.5) continue;

	    //don't overlap with tag
	    if ( deltaR(eleEta[indexTag],elePhi[indexTag], phoEta[indexProbe], phoPhi[indexProbe]) < 0.4 ) continue;

	    //For MC, Match to Gen level electron
	    if (!isData) {
	      bool genmatch = false;
	      for (int q=0;q<int(genLeptonIndex.size()); q++) {
		if ( deltaR(phoEta[indexProbe],elePhi[indexProbe],
			    gParticleEta[genLeptonIndex[q]],
			    gParticlePhi[genLeptonIndex[q]]) < 0.1) {
		  genmatch = true;
		}
	      }
	      if (!genmatch) continue;
	    }

	    // //Probe must match probe leg filter of the dedicated T&P trigger
	    // if ( !matchProbeElectronHLTFilters(indexTag)) continue;
	    // if ( !matchProbeSCHLTFilters(indexTag)) continue;

	    //cout << "Probe: " << phoPt[indexProbe] << " " << phoEta[indexProbe] << " " <<  phoPhi[indexProbe] << " : " << isLoosePhotonWithoutEleVeto(indexProbe) << "\n";

	    //*******************************************************
	    //denominator selection
	    //*******************************************************
	    if (denominatorType == 1) {
	      // reco object doesn't require any additional cuts
	    }
	    if (denominatorType == 3) {
	      if ( !isLoosePhotonWithoutEleVeto(indexProbe, _is25ns) ) continue;
	    }
	    if (denominatorType == 4) {
	      if ( !isMediumPhotonWithoutEleVeto(indexProbe, _is25ns) ) continue;
	    }
	    if (denominatorType == 5) {
	      if ( !isTightPhotonWithoutEleVeto(indexProbe, _is25ns) ) continue;
	    }
	    if (denominatorType == 14) {
	      if ( !isLooseDelayedPhotonWithoutEleVeto(indexProbe, _is25ns) ) continue;
	    }
	    if (denominatorType == 15) {
	      if ( !isMediumDelayedPhotonWithoutEleVeto(indexProbe, _is25ns) ) continue;
	    }
	    if (denominatorType == 16) {
	      if ( !isTightDelayedPhotonWithoutEleVeto(indexProbe, _is25ns) ) continue;
	    }

	    TLorentzVector vprobe;
	    vprobe.SetPtEtaPhiM(phoPt[indexProbe], phoEta[indexProbe], phoPhi[indexProbe], ELE_MASS);

	    TPPair->mass = (vtag+vprobe).M();
	    TPPair->pt = phoPt[indexProbe];
	    TPPair->eta = phoEta[indexProbe];
	    TPPair->phi = phoPhi[indexProbe];
	    //TPPair->charge = phoCharge[indexProbe];
	    TPPair->charge = 0;

	    //****************************************
	    //PASS OR FAIL
	    //****************************************
	    bool pass = false;
	    if (numeratorType == 3) {
	      pass = isLoosePhotonWithoutEleVeto(indexProbe, _is25ns);
	    }
	    if (numeratorType == 4) {
	      pass = isMediumPhotonWithoutEleVeto(indexProbe, _is25ns);
	    }
	    if (numeratorType == 5) {
	      pass = isTightPhotonWithoutEleVeto(indexProbe, _is25ns);
	    }
	    if (numeratorType == 14) {
	      pass = isLooseDelayedPhotonWithoutEleVeto(indexProbe, _is25ns);
	    }
	    if (numeratorType == 15) {
	      pass = isMediumDelayedPhotonWithoutEleVeto(indexProbe, _is25ns);
	    }
	    if (numeratorType == 16) {
	      pass = isTightDelayedPhotonWithoutEleVeto(indexProbe, _is25ns);
	    }
	    if (numeratorType == 50) {
	      pass = matchPhotonHLTFilters(indexProbe, "DiPhoton30_18_WithPixMatch_Leg1");
	    }
	    if (numeratorType == 51) {
	      pass = matchPhotonHLTFilters(indexProbe, "DiPhoton30_18_WithPixMatch_Leg2");
	    }
	    if (numeratorType == 52) {
	      pass = matchPhotonHLTFilters(indexProbe, "Photon42_Photon25_Mass15_Leg1");
	    }
	    if (numeratorType == 53) {
	      pass = matchPhotonHLTFilters(indexProbe, "Photon42_Photon25_Mass15_Leg2");
	    }


	    TPPair->pass = pass;
	    // cout << " TP Pass: " << TPPair->pass << "\n";

	    //****************************************
	    //Fill Output Tree
	    //****************************************
	    TPPair->tree_->Fill();

	  } //loop over probe electrons

	} // loop over tag electrons

      } //end if photon object


    }//end of event loop


    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}
