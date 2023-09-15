#include "llp_MuonSystem_bparking_rechits.h"
#include "RazorHelper.h"
#include "TreeMuonSystemBParking_rechits.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
#include "DBSCAN.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include <iostream>
#include <random>
//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"


using namespace std::chrono;
using namespace std;
using namespace ROOT::Math;

struct greater_than_pt
{
  inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
};

struct leptons
{
  TLorentzVector lepton;
  int pdgId;
  float dZ;
  float dXY;
  float dXYErr;
  // bool passLooseId;
  // bool passMediumId;
  bool tightId;
  bool looseId;
  bool passVetoId;
  bool passLooseIso;
  bool passTightIso;
  bool passVTightIso;
  bool passVVTightIso;
  //variables for muons only
  unsigned int muonType;
  unsigned int muonQuality;
  bool muon_passHLTFilter[MAX_MuonHLTFilters];

};


struct jets
{
  TLorentzVector jet;
  float time;
  bool passId;
  // bool passLooseId;
  // bool passMediumId;
  // bool passTightId;
  bool isCSVL;
  int ecalNRechits;
  float ecalRechitE;
  float jetChargedEMEnergyFraction;
  float jetNeutralEMEnergyFraction;
  float jetChargedHadronEnergyFraction;
  float jetNeutralHadronEnergyFraction;
  bool jetPassMuFrac;
  float jetPtJESUp;
  float jetPtJESDown;
  float jetEJESUp;
  float jetEJESDown;
  float JecUnc;

  float electronEnergyFraction;
  float neutralEmEnergyFraction;
  float chargedHadronEnergyFraction;
  float neutralHadronEnergyFraction;
  float muonEnergyFraction;

};

//lepton highest pt comparator
struct largest_pt
{
  inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt;

//jet highest pt comparator
struct largest_pt_jet
{
  inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;


void llp_MuonSystem_bparking_rechits::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  int option = options%10;

  if( isData )
  {
    std::cout << "[INFO]: running on data with option: " << option << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with option: " << option << std::endl;
  }

  const float ELE_MASS = 0.000511;
  const float MU_MASS  = 0.105658;
  const float Z_MASS   = 91.2;

  if (analysisTag == ""){
    // analysisTag = "Razor2018_17SeptEarlyReReco";
    analysisTag = "BParking_Source2018";
  }

  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile;
  outFile = new TFile(outfilename.c_str(), "RECREATE");

  TreeMuonSystemBParking_rechits *MuonSystem = new TreeMuonSystemBParking_rechits;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *MuTrig = new TH1F("MuTrig", "MuTrig", 1, 1, 2);
  TH1F *SingleMu = new TH1F("SingleMu", "SingleMu", 1, 1, 2);

//----  char* cmsswPath;
//----  cmsswPath = getenv("CMSSW_BASE");
//----  string pathname;
//----  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
//----  if(cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1

  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  helper = new RazorHelper(analysisTag, isData, false);
  helper->load_BParking_SF();

  //*************************************************************************
  //Look over Input File Events
  //*************************************************************************
  if (fChain == 0) return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nbytes = 0, nb = 0;
  clock_t start, end;
  start = clock();
  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {
    //begin event
    //if (jentry>20) break;
    if(jentry % 10000 == 0)
    {
      end = clock();
      double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
      cout << "Processing entry " << jentry << endl;
      cout << "Time taken by program is : " << time_taken << endl;
      start = clock();
    }
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //GetEntry(ientry);
    nb = fChain->GetEntry(jentry); nbytes += nb;

    //fill normalization histogram
    MuonSystem->InitVariables();

    //event info
    if (isData)
    {
      NEvents->Fill(1);
    }
    else
    {
      NEvents->Fill(1, genWeight);
    }

    MuonSystem->runNum = runNum;
    MuonSystem->lumiSec = lumiNum;
    MuonSystem->evtNum = eventNum;

    if (!isData)
    {
        if (analysisTag=="Razor2016_07Aug2017Rereco") MuonSystem->MC_condition = 2016;
        else if (analysisTag=="Razor2017_17Nov2017Rereco") MuonSystem->MC_condition = 2017;
        else if (analysisTag=="Razor2018_17SeptEarlyReReco") MuonSystem->MC_condition = 2018;
    }

    //get NPU
    MuonSystem->npv = nPV;
    MuonSystem->rho = fixedGridRhoFastjetAll;
    MuonSystem->met = metType1Pt;
    MuonSystem->metPhi = metType1Phi;

    for (int i=0; i < nBunchXing; i++)
    {
        if (BunchXing[i] == 0)
        {
            MuonSystem->npu = nPUmean[i];
        }
    }

    // cout<<"nBunchXing: "<<nBunchXing<<endl;
    // cout<<"nPUmean[1]: "<<nPUmean[1]<<endl;
    // MuonSystem->npu = 0;

    MuonSystem->pileupWeight = helper->getPileupWeight(MuonSystem->npu);
    MuonSystem->pileupWeightUp = helper->getPileupWeightUp(MuonSystem->npu) / MuonSystem->pileupWeight;
    MuonSystem->pileupWeightDown = helper->getPileupWeightDown(MuonSystem->npu) / MuonSystem->pileupWeight;

    std::pair<double,double> corrected_met;
    if (analysisTag=="Razor2016_07Aug2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2016, !isData, nPV);
    else if (analysisTag=="Razor2017_17Nov2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2017, !isData, nPV);
    else if (analysisTag=="Razor2018_17SeptEarlyReReco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2018, !isData, nPV);

    //Can keep as corrected MEt
    //----MuonSystem->metEENoise = corrected_met.first;
    //----MuonSystem->metPhiEENoise = corrected_met.second;


    //Triggers
    for(int i = 0; i < NTriggersMAX; i++){
      MuonSystem->HLTDecision[i] = HLTDecision[i];
      //std::cout<<HLTDecision[i];
    }
    //std::cout<<std::endl;
    bool passBParkingTrig = false;
    for(int i = 1157; i <= 1196; i++){
      if (MuonSystem->HLTDecision[i]) passBParkingTrig = true;
    }
    if (passBParkingTrig)
    {
      if (isData)
      {
        MuTrig->Fill(1);
      }
      else
      {
        MuTrig->Fill(1, genWeight);
      }
    }

    MuonSystem->Flag2_all = Flag2_HBHENoiseFilter && Flag2_HBHEIsoNoiseFilter && Flag2_BadPFMuonFilter && Flag2_globalSuperTightHalo2016Filter && Flag2_EcalDeadCellTriggerPrimitiveFilter;
    if (isData) MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_eeBadScFilter;

    int llp_mother = 0;
    for (int i = 0; i < nGenParticle; i++)
    {
      if (gParticleId[i] == 9900015) {
        llp_mother = gParticleMotherIndex[i];
        break;
      }
    }

    // cout<<nGenParticle<<endl;
    bool genMuonFlag = false;
    for(int i = 0; i < nGenParticle;i++)
    {
      if(gParticleId[i] == 9900015)
      {
        TLorentzVector genParticle = makeTLorentzVector( gParticlePt[i], gParticleEta[i], gParticlePhi[i], gParticleE[i] );

        MuonSystem->gLLP_eta = gParticleEta[i];
        MuonSystem->gLLP_phi = gParticlePhi[i];
        MuonSystem->gLLP_e = gParticleE[i];
        MuonSystem->gLLP_pt = gParticlePt[i];
        MuonSystem->gLLP_beta = genParticle.Beta();

      }
      else if (gParticleMotherId[i] == 9900015)
      {
        MuonSystem->gLLP_decay_vertex_x = gParticleProdVertexX[i];
        MuonSystem->gLLP_decay_vertex_y = gParticleProdVertexY[i];
        MuonSystem->gLLP_decay_vertex_z = gParticleProdVertexZ[i];
        MuonSystem->gLLP_decay_vertex_r = sqrt(gParticleProdVertexX[i]*gParticleProdVertexX[i]+gParticleProdVertexY[i]*gParticleProdVertexY[i]);
        float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r, 2) + pow(MuonSystem->gLLP_decay_vertex_z,2));
        float gamma = 1.0/sqrt(1-MuonSystem->gLLP_beta*MuonSystem->gLLP_beta);
        MuonSystem->gLLP_ctau = gLLP_decay_vertex/(MuonSystem->gLLP_beta * gamma);

      }
      if(abs(gParticleId[i]) == 13 && abs(gParticleEta[i])<1.5 && gParticlePt[i]>7)genMuonFlag = true;
      if (gParticleMotherIndex[i] == llp_mother) {
        MuonSystem->gParticlePt[MuonSystem->nGenParticles] = gParticlePt[i];
        MuonSystem->gParticleId[MuonSystem->nGenParticles] = gParticleId[i];
        MuonSystem->gParticleE[MuonSystem->nGenParticles] = gParticleE[i];
        MuonSystem->gParticlePhi[MuonSystem->nGenParticles] = gParticlePhi[i];
        MuonSystem->gParticleEta[MuonSystem->nGenParticles] = gParticleEta[i];
        MuonSystem->gParticleMotherId[MuonSystem->nGenParticles] = gParticleMotherId[i];
        MuonSystem->gParticleMotherIndex[MuonSystem->nGenParticles] = gParticleMotherIndex[i];
        MuonSystem->nGenParticles++;
      }
    }
    //if (!genMuonFlag)continue;
    if (abs(MuonSystem->gLLP_eta) < 2.4
        && abs(MuonSystem->gLLP_decay_vertex_z)<1100 && abs(MuonSystem->gLLP_decay_vertex_z)>400
        && MuonSystem->gLLP_decay_vertex_r < 695.5
       ) MuonSystem->gLLP_csc = true;
    if (abs(MuonSystem->gLLP_decay_vertex_z)< 661.0
        && MuonSystem->gLLP_decay_vertex_r < 800
        && MuonSystem->gLLP_decay_vertex_r > 200.0
       ) MuonSystem->gLLP_dt = true;

    //*************************************************************************
    //Start Object Selection
    //*************************************************************************
    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
    bool SingleMu_flag = false;
    for( int i = 0; i < nMuons; i++ )
    {
      if(fabs(muonEta[i]) > 2.4) continue;
      leptons tmpMuon;
      tmpMuon.lepton.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i], MU_MASS);
      tmpMuon.pdgId = 13 * -1 * muonCharge[i];
      tmpMuon.dZ = muon_dZ[i];
      tmpMuon.dXY = muon_d0[i];
      tmpMuon.dXYErr = muon_d0Err[i];
      tmpMuon.tightId = isMuonPOGTightMuon(i);
      tmpMuon.looseId = isMuonPOGLooseMuon(i);
      float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
      tmpMuon.passLooseIso = muonIso<0.25;
      tmpMuon.passTightIso = muonIso<0.15;
      tmpMuon.passVTightIso = muonIso<0.10;
      tmpMuon.passVVTightIso = muonIso<0.05;
      tmpMuon.passVetoId = false;
     tmpMuon.muonType = muonType[i];
     tmpMuon.muonQuality = muonQuality[i];
     for(int j=0; j<MAX_MuonHLTFilters; j++) tmpMuon.muon_passHLTFilter[j] = muon_passHLTFilter[i][j];
      if (muonPt[i] > 8 && tmpMuon.looseId) SingleMu_flag = true;
     Leptons.push_back(tmpMuon);
    }
    if (SingleMu_flag)
    {
      if (isData)
      {
        SingleMu->Fill(1);
      }
      else
      {
        SingleMu->Fill(1, genWeight);
      }
    }
    //-------------------------------
    //Electrons
    //-------------------------------
    for( int i = 0; i < nElectrons; i++ )
    {
      if(fabs(eleEta[i]) > 2.5) continue;
      //remove overlaps
      bool overlap = false;
      for(auto& lep : Leptons)
      {
        if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
      }
      if(overlap) continue;
      leptons tmpElectron;
      tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);
      tmpElectron.pdgId = 11 * -1 * eleCharge[i];
      tmpElectron.dZ = ele_dZ[i];
      tmpElectron.tightId = isEGammaPOGTightElectron(i, true, true, true, "Summer16");
      tmpElectron.passVetoId = isEGammaPOGVetoElectron(i, true, true, true, "Summer16");
      Leptons.push_back(tmpElectron);
    }
    sort(Leptons.begin(), Leptons.end(), my_largest_pt);

    for ( auto &tmp : Leptons )
    {
      MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
      MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
      MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
      MuonSystem->lepDXY[MuonSystem->nLeptons]     = tmp.dXY;
      MuonSystem->lepDXYErr[MuonSystem->nLeptons]     = tmp.dXYErr;
      MuonSystem->lepSF[MuonSystem->nLeptons]     = helper->getBParkingTriggerSF(tmp.lepton.Pt(), tmp.dXY/tmp.dXYErr);
      MuonSystem->lepLooseId[MuonSystem->nLeptons] = tmp.looseId;
      MuonSystem->lepTightId[MuonSystem->nLeptons] = tmp.tightId;
      MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
      MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
      MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
      MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;

	    if (abs(tmp.pdgId)==13)
	    {
		     MuonSystem->lepMuonType[MuonSystem->nLeptons] = tmp.muonType;//only assigned for muons
    		 MuonSystem->lepMuonQuality[MuonSystem->nLeptons] = tmp.muonQuality;//only assigned for muons
    		 //MuonSystem->lepMuon_passHLTFilter[MuonSystem->nLeptons] = tmp.muon_passHLTFilter;//only assigned for muons
	       for(int j=0; j<MAX_MuonHLTFilters; j++) MuonSystem->lepMuon_passHLTFilter[MuonSystem->nLeptons][j] = tmp.muon_passHLTFilter[j];
	    }
            MuonSystem->nLeptons++;
      }

        //-----------------------------------------------
        //Select Jets
        //-----------------------------------------------

        std::vector<jets> Jets;
        float MetXCorr_JESUp = 0.;
        float MetYCorr_JESUp = 0.;
        float MetXCorr_JESDown = 0.;
        float MetYCorr_JESDown = 0.;
        float MetXCorr_HEM = 0.;
        float MetYCorr_HEM = 0.;
        float MetXCorr_EENoise = 0.;
        float MetYCorr_EENoise = 0.;

        for(int i = 0; i < nJets; i++)
        {
          //------------------------------------------------------------
          //exclude selected muons and electrons from the jet collection
          //------------------------------------------------------------
          double deltaR = -1;
          for(auto& lep : Leptons){
            double thisDR = RazorAnalyzer::deltaR(jetEta[i],jetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
            if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
          }
          if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

          double jetCorrPt = jetPt[i];
          double jetCorrE = jetE[i];
          TLorentzVector thisJet = makeTLorentzVector( jetCorrPt, jetEta[i], jetPhi[i], jetCorrE );

          if (fabs(thisJet.Eta()) >= 3.0)continue;

          jets tmpJet;
          tmpJet.jet    = thisJet;
          tmpJet.passId = isPFTightJet(i, true,analysisTag);

          Jets.push_back(tmpJet);

          }

          sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

          for ( auto &tmp : Jets )
          {
            if(tmp.jet.Pt()<30)continue;
            if(abs(tmp.jet.Eta())>2.5)continue;
            MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
            MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
            MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
            MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
            MuonSystem->jetTightPassId[MuonSystem->nJets] = tmp.passId;
            MuonSystem->nJets++;
          }

      MuonSystem->nDtRechits  = nDtRechits;
      int nDTRechitsChamberMinus12 = 0, nDTRechitsChamberMinus11 = 0, nDTRechitsChamber10 = 0, nDTRechitsChamberPlus11 = 0, nDTRechitsChamberPlus12 = 0,
      nDTRechitsChamberMinus22 = 0, nDTRechitsChamberMinus21 = 0, nDTRechitsChamber20 = 0, nDTRechitsChamberPlus21 = 0, nDTRechitsChamberPlus22 = 0,
      nDTRechitsChamberMinus32 = 0, nDTRechitsChamberMinus31 = 0, nDTRechitsChamber30 = 0, nDTRechitsChamberPlus31 = 0, nDTRechitsChamberPlus32 = 0,
      nDTRechitsChamberMinus42 = 0, nDTRechitsChamberMinus41 = 0, nDTRechitsChamber40 = 0, nDTRechitsChamberPlus41 = 0, nDTRechitsChamberPlus42 = 0;
      int nDTRechitsStation1 = 0,
      nDTRechitsStation2 = 0, nDTRechitsStation3 = 0, nDTRechitsStation4 = 0,
      nDTRechitsWheelMinus2 = 0, nDTRechitsWheelMinus1 = 0, nDTRechitsWheel0 = 0, nDTRechitsWheelPlus1 = 0, nDTRechitsWheelPlus2 = 0;

      for (int i = 0; i < nDtRechits; i++) {

        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus12++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 0) nDTRechitsChamber10++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus12++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus22++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 0) nDTRechitsChamber20++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus22++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus32++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 0) nDTRechitsChamber30++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus32++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus42++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 0) nDTRechitsChamber40++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus42++;

        if (dtRechitStation[i] == 1) nDTRechitsStation1++;
        if (dtRechitStation[i] == 2) nDTRechitsStation2++;
        if (dtRechitStation[i] == 3) nDTRechitsStation3++;
        if (dtRechitStation[i] == 4) nDTRechitsStation4++;

        if (dtRechitWheel[i] == -2) nDTRechitsWheelMinus2++;
        if (dtRechitWheel[i] == -1) nDTRechitsWheelMinus1++;
        if (dtRechitWheel[i] == 0) nDTRechitsWheel0++;
        if (dtRechitWheel[i] == 1) nDTRechitsWheelPlus1++;
        if (dtRechitWheel[i] == 2) nDTRechitsWheelPlus2++;
      }

      if ( nDTRechitsChamberMinus12 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus11 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber10 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus11 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus12 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus22 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus21 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber20 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus21 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus22 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus32 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus31 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber30 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus31 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus32 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus42 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus41 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber40 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus41 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus42 > 50) MuonSystem->nDtRings++;

      if (nDTRechitsStation1 > 25) MuonSystem->nDtStations25++;
      if (nDTRechitsStation2 > 25) MuonSystem->nDtStations25++;
      if (nDTRechitsStation3 > 25) MuonSystem->nDtStations25++;
      if (nDTRechitsStation4 > 25) MuonSystem->nDtStations25++;
      if (nDTRechitsWheelMinus2 > 25) MuonSystem->nDtWheels25++;
      if (nDTRechitsWheelMinus1 > 25) MuonSystem->nDtWheels25++;
      if (nDTRechitsWheel0 > 25) MuonSystem->nDtWheels25++;
      if (nDTRechitsWheelPlus1 > 25) MuonSystem->nDtWheels25++;
      if (nDTRechitsWheelPlus2 > 25) MuonSystem->nDtWheels25++;


      vector<Point> points;
      vector<int> cscRechitsClusterId;
      points.clear();
      MuonSystem->nCscRechits  = 0;

      int nCscRechitsChamberPlus11 = 0, nCscRechitsChamberPlus12 = 0, nCscRechitsChamberPlus13 = 0, nCscRechitsChamberPlus21 = 0, nCscRechitsChamberPlus22 = 0,
      nCscRechitsChamberPlus31 = 0, nCscRechitsChamberPlus32 = 0, nCscRechitsChamberPlus41 = 0, nCscRechitsChamberPlus42 = 0;
      int nCscRechitsChamberMinus11 = 0, nCscRechitsChamberMinus12 = 0, nCscRechitsChamberMinus13 = 0, nCscRechitsChamberMinus21 = 0, nCscRechitsChamberMinus22 = 0,
      nCscRechitsChamberMinus31 = 0, nCscRechitsChamberMinus32 = 0, nCscRechitsChamberMinus41 = 0, nCscRechitsChamberMinus42 = 0;
      for (int i = 0; i < ncscRechits; i++) {
        //pick out the right bits for chamber
        int chamber = ((cscRechitsDetId[i] >> 3) & 077); //https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonDetId/interface/CSCDetId.h#L147
        int layer = (cscRechitsDetId[i] & 07);

        MuonSystem->point_clusterID[i] = 100;
        // std::cout << i << std::endl;

        Point p;
        p.phi = cscRechitsPhi[i];
        p.eta = cscRechitsEta[i];
        p.x = cscRechitsX[i];
        p.y = cscRechitsY[i];
        p.z = cscRechitsZ[i];
        p.t = cscRechitsTpeak[i];
        p.twire = cscRechitsTwire[i];
        p.station = cscRechitsStation[i];
        p.chamber = cscRechitsChamber[i];
        p.layer = layer;
        p.superlayer = 0;
        p.clusterID = UNCLASSIFIED;
        points.push_back(p);
        cscRechitsClusterId.push_back(-1);

        
          
        if (cscRechitsChamber[i] == 11) nCscRechitsChamberPlus11++;
        if (cscRechitsChamber[i] == 12) nCscRechitsChamberPlus12++;
        if (cscRechitsChamber[i] == 13) nCscRechitsChamberPlus13++;
        if (cscRechitsChamber[i] == 21) nCscRechitsChamberPlus21++;
        if (cscRechitsChamber[i] == 22) nCscRechitsChamberPlus22++;
        if (cscRechitsChamber[i] == 31) nCscRechitsChamberPlus31++;
        if (cscRechitsChamber[i] == 32) nCscRechitsChamberPlus32++;
        if (cscRechitsChamber[i] == 41) nCscRechitsChamberPlus41++;
        if (cscRechitsChamber[i] == 42) nCscRechitsChamberPlus42++;
        if (cscRechitsChamber[i] == -11) nCscRechitsChamberMinus11++;
        if (cscRechitsChamber[i] == -12) nCscRechitsChamberMinus12++;
        if (cscRechitsChamber[i] == -13) nCscRechitsChamberMinus13++;
        if (cscRechitsChamber[i] == -21) nCscRechitsChamberMinus21++;
        if (cscRechitsChamber[i] == -22) nCscRechitsChamberMinus22++;
        if (cscRechitsChamber[i] == -31) nCscRechitsChamberMinus31++;
        if (cscRechitsChamber[i] == -32) nCscRechitsChamberMinus32++;
        if (cscRechitsChamber[i] == -41) nCscRechitsChamberMinus41++;
        if (cscRechitsChamber[i] == -42) nCscRechitsChamberMinus42++;
      }
      MuonSystem->nCscRings = 0;
      if ( nCscRechitsChamberPlus11 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus12 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus13 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus21 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus22 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus31 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus32 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus41 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus42 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus11 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus12 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus13 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus21 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus22 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus31 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus32 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus41 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus42 > 50) MuonSystem->nCscRings++;

      //Do DBSCAN Clustering
      int min_point = 50;  //minimum number of segments to call it a cluster
      float epsilon = 0.2; //cluster radius parameter
      DBSCAN ds(min_point, epsilon, points);
      ds.run();
      ds.result();

      ds.clusterMoments();
      ds.sort_clusters();
      ds.merge_clusters();
      ds.result();
      ds.clusterMoments();
      ds.sort_clusters();
          
      MuonSystem->nCscRechitClusters = 0;
      for ( auto &tmp : ds.CscCluster ) {
        MuonSystem->cscRechitClusterX[MuonSystem->nCscRechitClusters] =tmp.x;
          
        std::cout << MuonSystem->nCscRechitClusters << " " << tmp.x << std::endl;
          
        MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] =tmp.y;
        MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] =tmp.z;
        MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] = tmp.t;
        MuonSystem->cscRechitClusterTimeTotal[MuonSystem->nCscRechitClusters] = tmp.tTotal;
        MuonSystem->cscRechitClusterTimeWeighted[MuonSystem->nCscRechitClusters] = tmp.tWeighted;
        MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] =tmp.eta;
        MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters] = tmp.phi;
	  for (int j = 0; j < nTracks; j++) {
	    //printf("Tracks: %10f %10f %10f \n", track_Pt [j], track_Eta[j], track_Phi[j]);
            // ----------------------------------------------------------------------------------------------------------------
	    float dr = RazorAnalyzer::deltaR(track_Eta[j], track_Phi[j], tmp.eta, tmp.phi);
	    if (dr < 0.2 && track_dzToPV[j] < 0.5 && track_dxyToBS[j] < 0.2) {
		MuonSystem->cscRechitClusterMatchedTrackSumPt_trk_pos_0p2[MuonSystem->nCscRechitClusters] += track_Pt[j];
		MuonSystem->cscRechitClusterMatchedTrackSize_trk_pos_0p2[MuonSystem->nCscRechitClusters]++;
		if (track_Pt[j] > MuonSystem->cscRechitClusterMatchedTrackLeadPt_trk_pos_0p2[MuonSystem->nCscRechitClusters]) {
		    MuonSystem->cscRechitClusterMatchedTrackLeadPt_trk_pos_0p2[MuonSystem->nCscRechitClusters] = track_Pt[j];
		}
	    }
	    if (dr < 0.2) {
		MuonSystem->cscRechitClusterMatchedTrackSumPt_0p2[MuonSystem->nCscRechitClusters] += track_Pt[j];
		MuonSystem->cscRechitClusterMatchedTrackSize_0p2[MuonSystem->nCscRechitClusters]++;
		if (track_Pt[j] > MuonSystem->cscRechitClusterMatchedTrackLeadPt_0p2[MuonSystem->nCscRechitClusters]) {
		    MuonSystem->cscRechitClusterMatchedTrackLeadPt_0p2[MuonSystem->nCscRechitClusters] = track_Pt[j];
		}
	    }
	    // ----------------------------------------------------------------------------------------------------------------
	    if (dr < 0.3 && track_dzToPV[j] < 0.5 && track_dxyToBS[j] < 0.2) {
		MuonSystem->cscRechitClusterMatchedTrackSumPt_trk_pos_0p3[MuonSystem->nCscRechitClusters] += track_Pt[j];
		MuonSystem->cscRechitClusterMatchedTrackSize_trk_pos_0p3[MuonSystem->nCscRechitClusters]++;
		if (track_Pt[j] > MuonSystem->cscRechitClusterMatchedTrackLeadPt_trk_pos_0p3[MuonSystem->nCscRechitClusters]) {
		    MuonSystem->cscRechitClusterMatchedTrackLeadPt_trk_pos_0p3[MuonSystem->nCscRechitClusters] = track_Pt[j];
		}
	    }
	    if (dr < 0.3) {
		MuonSystem->cscRechitClusterMatchedTrackSumPt_0p3[MuonSystem->nCscRechitClusters] += track_Pt[j];
		MuonSystem->cscRechitClusterMatchedTrackSize_0p3[MuonSystem->nCscRechitClusters]++;
		if (track_Pt[j] > MuonSystem->cscRechitClusterMatchedTrackLeadPt_0p3[MuonSystem->nCscRechitClusters]) {
		    MuonSystem->cscRechitClusterMatchedTrackLeadPt_0p3[MuonSystem->nCscRechitClusters] = track_Pt[j];
		}
	    }
	    // ----------------------------------------------------------------------------------------------------------------
	    if (dr < 0.4 && track_dzToPV[j] < 0.5 && track_dxyToBS[j] < 0.2) {
		MuonSystem->cscRechitClusterMatchedTrackSumPt_trk_pos_0p4[MuonSystem->nCscRechitClusters] += track_Pt[j];
		MuonSystem->cscRechitClusterMatchedTrackSize_trk_pos_0p4[MuonSystem->nCscRechitClusters]++;
		if (track_Pt[j] > MuonSystem->cscRechitClusterMatchedTrackLeadPt_trk_pos_0p4[MuonSystem->nCscRechitClusters]) {
		    MuonSystem->cscRechitClusterMatchedTrackLeadPt_trk_pos_0p4[MuonSystem->nCscRechitClusters] = track_Pt[j];
		}
	    }
	    if (dr < 0.4) {
		MuonSystem->cscRechitClusterMatchedTrackSumPt_0p4[MuonSystem->nCscRechitClusters] += track_Pt[j];
		MuonSystem->cscRechitClusterMatchedTrackSize_0p4[MuonSystem->nCscRechitClusters]++;
		if (track_Pt[j] > MuonSystem->cscRechitClusterMatchedTrackLeadPt_0p4[MuonSystem->nCscRechitClusters]) {
		    MuonSystem->cscRechitClusterMatchedTrackLeadPt_0p4[MuonSystem->nCscRechitClusters] = track_Pt[j];
		}
	    }
	    // ----------------------------------------------------------------------------------------------------------------
	    if (dr < 0.5 && track_dzToPV[j] < 0.5 && track_dxyToBS[j] < 0.2) {
		MuonSystem->cscRechitClusterMatchedTrackSumPt_trk_pos_0p5[MuonSystem->nCscRechitClusters] += track_Pt[j];
		MuonSystem->cscRechitClusterMatchedTrackSize_trk_pos_0p5[MuonSystem->nCscRechitClusters]++;
		if (track_Pt[j] > MuonSystem->cscRechitClusterMatchedTrackLeadPt_trk_pos_0p5[MuonSystem->nCscRechitClusters]) {
		    MuonSystem->cscRechitClusterMatchedTrackLeadPt_trk_pos_0p5[MuonSystem->nCscRechitClusters] = track_Pt[j];
		}
	    }
	    if (dr < 0.5) {
		MuonSystem->cscRechitClusterMatchedTrackSumPt_0p5[MuonSystem->nCscRechitClusters] += track_Pt[j];
		MuonSystem->cscRechitClusterMatchedTrackSize_0p5[MuonSystem->nCscRechitClusters]++;
		if (track_Pt[j] > MuonSystem->cscRechitClusterMatchedTrackLeadPt_0p5[MuonSystem->nCscRechitClusters]) {
		    MuonSystem->cscRechitClusterMatchedTrackLeadPt_0p5[MuonSystem->nCscRechitClusters] = track_Pt[j];
		}
	    }
	  }
          MuonSystem->cscRechitClusterTimeSpread[MuonSystem->nCscRechitClusters] = tmp.TSpread;
          MuonSystem->cscRechitClusterTimeSpreadWeighted[MuonSystem->nCscRechitClusters] = tmp.TSpreadWeighted;
          MuonSystem->cscRechitClusterTimeSpreadWeightedAll[MuonSystem->nCscRechitClusters] = tmp.TSpreadWeightedAll;
          MuonSystem->cscRechitClusterSize[MuonSystem->nCscRechitClusters] = tmp.nCscSegments;

          MuonSystem->cscRechitClusterNRechitChamberPlus11[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus11;
          MuonSystem->cscRechitClusterNRechitChamberPlus12[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus12;
          MuonSystem->cscRechitClusterNRechitChamberPlus13[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus13;
          MuonSystem->cscRechitClusterNRechitChamberPlus21[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus21;
          MuonSystem->cscRechitClusterNRechitChamberPlus22[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus22;
          MuonSystem->cscRechitClusterNRechitChamberPlus31[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus31;
          MuonSystem->cscRechitClusterNRechitChamberPlus32[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus32;
          MuonSystem->cscRechitClusterNRechitChamberPlus41[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus41;
          MuonSystem->cscRechitClusterNRechitChamberPlus42[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberPlus42;
          MuonSystem->cscRechitClusterNRechitChamberMinus11[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus11;
          MuonSystem->cscRechitClusterNRechitChamberMinus12[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus12;
          MuonSystem->cscRechitClusterNRechitChamberMinus13[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus13;
          MuonSystem->cscRechitClusterNRechitChamberMinus21[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus21;
          MuonSystem->cscRechitClusterNRechitChamberMinus22[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus22;
          MuonSystem->cscRechitClusterNRechitChamberMinus31[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus31;
          MuonSystem->cscRechitClusterNRechitChamberMinus32[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus32;
          MuonSystem->cscRechitClusterNRechitChamberMinus41[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus41;
          MuonSystem->cscRechitClusterNRechitChamberMinus42[MuonSystem->nCscRechitClusters] = tmp.nCscSegmentChamberMinus42;
          MuonSystem->cscRechitClusterMaxChamber[MuonSystem->nCscRechitClusters] = tmp.maxChamber;
          MuonSystem->cscRechitClusterMaxChamberRatio[MuonSystem->nCscRechitClusters] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
          MuonSystem->cscRechitClusterNChamber[MuonSystem->nCscRechitClusters] = tmp.nChamber;
          MuonSystem->cscRechitClusterMaxStation[MuonSystem->nCscRechitClusters] = tmp.maxStation;
          MuonSystem->cscRechitClusterMaxStationRatio[MuonSystem->nCscRechitClusters] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;

          MuonSystem->cscRechitClusterNStation10[MuonSystem->nCscRechitClusters] = tmp.nStation10;
          MuonSystem->cscRechitClusterAvgStation10[MuonSystem->nCscRechitClusters] = tmp.avgStation10;

          // MuonSystem->cscRechitClusterMe11Ratio[MuonSystem->nCscRechitClusters] = tmp.Me11Ratio;
          // MuonSystem->cscRechitClusterMe12Ratio[MuonSystem->nCscRechitClusters] = tmp.Me12Ratio;

          //Jet veto/ muon veto
          MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;

          
            
          if (!isData)
          {

            for(int i = 0; i < nGenParticle; i++)
            {
	      if (fabs(gParticleId[i])!=13)continue;
              if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && gParticlePt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] ) {
                MuonSystem->cscRechitClusterGenMuonVetoPt[MuonSystem->nCscRechitClusters]  = gParticlePt[i];
              }
              if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.8 && gParticlePt[i] > MuonSystem->cscRechitClusterGenMuonVetoPt_dR0p8[MuonSystem->nCscRechitClusters] ) {
                MuonSystem->cscRechitClusterGenMuonVetoPt_dR0p8[MuonSystem->nCscRechitClusters]  = gParticlePt[i];
              }
            }
          }
          // jet veto
          for(int i = 0; i < nJets; i++)
          {
            if (fabs(jetEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]  = jetPt[i];              // MuonSystem->cscRechitClusterJetVetoPhi[MuonSystem->nCscRechitClusters]  = jetPhi[i];
            }

          }
          for(int i = 0; i < nMuons; i++)
          {
            if (fabs(muonEta[i]>3.0)) continue;
            float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
            if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]  = muonPt[i];
            }
          }


          //match to MB1 DT segments
          for (int i = 0; i < nDtSeg; i++) {
            if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
            {
              MuonSystem->cscRechitCluster_match_dtSeg_0p4[MuonSystem->nCscRechitClusters]++;
              if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters]++;
            }

          }
          //match to RPC hits in RE1/2
          for (int i = 0; i < nRpc; i++) {
            float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
            if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
            {
              if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)
              {
                MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters]++;
              }
              if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
              {
                MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters]++;
              }

            }
          }




          MuonSystem->cscRechitCluster_match_gLLP_deltaR[MuonSystem->nCscRechitClusters] = RazorAnalyzer::deltaR(MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->gLLP_eta, MuonSystem->gLLP_phi);
// can make alternate with MET and corr MET
//          MuonSystem->cscRechitClusterMetEENoise_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhiEENoise);

          //// =========================== CUTS ===========================
          //
          //if (MuonSystem->cscRechitClusterSize[MuonSystem->nCscRechitClusters] < 50) continue;
          //if (MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] > 20) continue;
          //if (MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters] != 0) continue;
          //if (MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters] != 0) continue;
          //if (MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters] != 0) continue;

          MuonSystem->nCscRechitClusters++;
      }
      // DT cluster

      points.clear();
      for (int i = 0; i < nDtRechits; i++) {
        Point p;
        // p.phi = dtRechitPhi[i];
        // p.eta = dtRechitEta[i];
        // p.x = dtRechitX[i];
        // p.y = dtRechitY[i];
        // p.z = dtRechitZ[i];
        p.t = dtRechitTime[i];
        p.twire = dtRechitTime[i];
        p.station = dtRechitStation[i];
        p.chamber = dtRechitWheel[i];
        // p.superlayer = 0;

        p.phi = dtRechitCorrectPhi[i];
        p.eta = dtRechitCorrectEta[i];
        p.x = dtRechitCorrectX[i];
        p.y = dtRechitCorrectY[i];
        p.z = dtRechitCorrectZ[i];
        p.superlayer = dtRechitSuperLayer[i];

        p.clusterID = UNCLASSIFIED;
        points.push_back(p);

      }
      //Do DBSCAN Clustering
      int min_point_dt = isData?50:30;  //minimum number of segments to call it a cluster
      float epsilon_dt = 0.2; //cluster radius parameter
      DBSCAN ds_dtRechit(min_point_dt, epsilon_dt, points);
      ds_dtRechit.run();
      ds_dtRechit.result();
      ds_dtRechit.clusterMoments();
      ds_dtRechit.sort_clusters();
      ds_dtRechit.merge_clusters();
      ds_dtRechit.result();
      ds_dtRechit.clusterMoments();
      ds_dtRechit.sort_clusters();


      MuonSystem->nDtRechitClusters = 0;

      for ( auto &tmp : ds_dtRechit.CscCluster ) {
        //remove overlaps
        bool overlap = false;
        for(int i = 0; i < MuonSystem->nCscRechitClusters; i++)
        {
          if (RazorAnalyzer::deltaR(MuonSystem->cscRechitClusterEta[i],MuonSystem->cscRechitClusterPhi[i],tmp.eta, tmp.phi)<0.4) overlap = true;
        }
        if(overlap) MuonSystem->dtRechitClusterOverlap[MuonSystem->nDtRechitClusters] = true;
        //
        // for (unsigned int j = 0; j < tmp.segment_id.size(); j++)
        // {
        //  MuonSystem->dtRechitsClusterId[tmp.segment_id[j]] = MuonSystem->nDtRechitClusters;
        // }
          MuonSystem->dtRechitClusterX[MuonSystem->nDtRechitClusters] =tmp.x;
          MuonSystem->dtRechitClusterY[MuonSystem->nDtRechitClusters] =tmp.y;
          MuonSystem->dtRechitClusterZ[MuonSystem->nDtRechitClusters] =tmp.z;
	  for (int j = 0; j < nTracks; j++) {
	    float dr = RazorAnalyzer::deltaR(track_Eta[j], track_Phi[j], tmp.eta, tmp.phi);
	    if (dr < 0.2 && track_dzToPV[j] < 0.5 && track_dxyToBS[j] < 0.2) {
		    MuonSystem->dtRechitClusterMatchedTrackSumPt_trk_pos_0p2[MuonSystem->nDtRechitClusters] += track_Pt[j];
		    MuonSystem->dtRechitClusterMatchedTrackSize_trk_pos_0p2[MuonSystem->nDtRechitClusters]++;
		    if (track_Pt[j] > MuonSystem->dtRechitClusterMatchedTrackLeadPt_trk_pos_0p2[MuonSystem->nDtRechitClusters]) {
		      MuonSystem->dtRechitClusterMatchedTrackLeadPt_trk_pos_0p2[MuonSystem->nDtRechitClusters] = track_Pt[j];
		    }
	    }
	    if (dr < 0.2) {
		    MuonSystem->dtRechitClusterMatchedTrackSumPt_0p2[MuonSystem->nDtRechitClusters] += track_Pt[j];
		    MuonSystem->dtRechitClusterMatchedTrackSize_0p2[MuonSystem->nDtRechitClusters]++;
		    if (track_Pt[j] > MuonSystem->dtRechitClusterMatchedTrackLeadPt_0p2[MuonSystem->nDtRechitClusters]) {
		      MuonSystem->dtRechitClusterMatchedTrackLeadPt_0p2[MuonSystem->nDtRechitClusters] = track_Pt[j];
		    }
	    }
	    // ----------------------------------------------------------------------------------------------------------------
	    if (dr < 0.3 && track_dzToPV[j] < 0.5 && track_dxyToBS[j] < 0.2) {
		MuonSystem->dtRechitClusterMatchedTrackSumPt_trk_pos_0p3[MuonSystem->nDtRechitClusters] += track_Pt[j];
		MuonSystem->dtRechitClusterMatchedTrackSize_trk_pos_0p3[MuonSystem->nDtRechitClusters]++;
		if (track_Pt[j] > MuonSystem->dtRechitClusterMatchedTrackLeadPt_trk_pos_0p3[MuonSystem->nDtRechitClusters]) {
		    MuonSystem->dtRechitClusterMatchedTrackLeadPt_trk_pos_0p3[MuonSystem->nDtRechitClusters] = track_Pt[j];
		}
	    }
	    if (dr < 0.3) {
		MuonSystem->dtRechitClusterMatchedTrackSumPt_0p3[MuonSystem->nDtRechitClusters] += track_Pt[j];
		MuonSystem->dtRechitClusterMatchedTrackSize_0p3[MuonSystem->nDtRechitClusters]++;
		if (track_Pt[j] > MuonSystem->dtRechitClusterMatchedTrackLeadPt_0p3[MuonSystem->nDtRechitClusters]) {
		    MuonSystem->dtRechitClusterMatchedTrackLeadPt_0p3[MuonSystem->nDtRechitClusters] = track_Pt[j];
		}
	    }
	    // ----------------------------------------------------------------------------------------------------------------
	    if (dr < 0.4 && track_dzToPV[j] < 0.5 && track_dxyToBS[j] < 0.2) {
		MuonSystem->dtRechitClusterMatchedTrackSumPt_trk_pos_0p4[MuonSystem->nDtRechitClusters] += track_Pt[j];
		MuonSystem->dtRechitClusterMatchedTrackSize_trk_pos_0p4[MuonSystem->nDtRechitClusters]++;
		if (track_Pt[j] > MuonSystem->dtRechitClusterMatchedTrackLeadPt_trk_pos_0p4[MuonSystem->nDtRechitClusters]) {
		    MuonSystem->dtRechitClusterMatchedTrackLeadPt_trk_pos_0p4[MuonSystem->nDtRechitClusters] = track_Pt[j];
		}
	    }
	    if (dr < 0.4) {
		MuonSystem->dtRechitClusterMatchedTrackSumPt_0p4[MuonSystem->nDtRechitClusters] += track_Pt[j];
		MuonSystem->dtRechitClusterMatchedTrackSize_0p4[MuonSystem->nDtRechitClusters]++;
		if (track_Pt[j] > MuonSystem->dtRechitClusterMatchedTrackLeadPt_0p4[MuonSystem->nDtRechitClusters]) {
		    MuonSystem->dtRechitClusterMatchedTrackLeadPt_0p4[MuonSystem->nDtRechitClusters] = track_Pt[j];
		}
	    }
	    // ----------------------------------------------------------------------------------------------------------------
	    if (dr < 0.5 && track_dzToPV[j] < 0.5 && track_dxyToBS[j] < 0.2) {
		MuonSystem->dtRechitClusterMatchedTrackSumPt_trk_pos_0p5[MuonSystem->nDtRechitClusters] += track_Pt[j];
		MuonSystem->dtRechitClusterMatchedTrackSize_trk_pos_0p5[MuonSystem->nDtRechitClusters]++;
		if (track_Pt[j] > MuonSystem->dtRechitClusterMatchedTrackLeadPt_trk_pos_0p5[MuonSystem->nDtRechitClusters]) {
		    MuonSystem->dtRechitClusterMatchedTrackLeadPt_trk_pos_0p5[MuonSystem->nDtRechitClusters] = track_Pt[j];
		}
	    }
	    if (dr < 0.5) {
		MuonSystem->dtRechitClusterMatchedTrackSumPt_0p5[MuonSystem->nDtRechitClusters] += track_Pt[j];
		MuonSystem->dtRechitClusterMatchedTrackSize_0p5[MuonSystem->nDtRechitClusters]++;
		if (track_Pt[j] > MuonSystem->dtRechitClusterMatchedTrackLeadPt_0p5[MuonSystem->nDtRechitClusters]) {
		    MuonSystem->dtRechitClusterMatchedTrackLeadPt_0p5[MuonSystem->nDtRechitClusters] = track_Pt[j];
		}
	    }
	  }
	  
          if (abs(tmp.z) < 126.8) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 0;
          else if (tmp.z > 126.8 && tmp.z < 395.4) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 1;
          else if (tmp.z < -126.8 && tmp.z > -395.4)MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -1;
          else if (tmp.z<0) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -2;
          else MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 2;


          MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters] =tmp.eta;
          MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters] =tmp.phi;

          MuonSystem->dtRechitClusterSize[MuonSystem->nDtRechitClusters] = tmp.nCscSegments;

          unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
          default_random_engine generator (seed);

	        // default_random_engine generator;
    	    uniform_real_distribution<double> distribution(0.0,1.0);
	         float prob = 0.03;
         for (int i=0; i<12; ++i) {
           if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters]++;
         }
         for (int i=0; i<12; ++i) {
           if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters]++;
      	 }
         for (int i=0; i<12; ++i) {
           if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters]++;
      	 }
         for (int i=0; i<8; ++i) {
           if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters]++;
      	 }

         MuonSystem->dtRechitClusterNoiseHit[MuonSystem->nDtRechitClusters] = MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters] +
                                                                              MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters] +
                                                                              MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters] +
                                                                              MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters];


          MuonSystem->dtRechitClusterNSegmentStation1[MuonSystem->nDtRechitClusters] = tmp.nDtSegmentStation1;
        	MuonSystem->dtRechitClusterNSegmentStation2[MuonSystem->nDtRechitClusters] = tmp.nDtSegmentStation2;
        	MuonSystem->dtRechitClusterNSegmentStation3[MuonSystem->nDtRechitClusters] = tmp.nDtSegmentStation3;
        	MuonSystem->dtRechitClusterNSegmentStation4[MuonSystem->nDtRechitClusters] = tmp.nDtSegmentStation4;

        	MuonSystem->dtRechitClusterMaxChamber[MuonSystem->nDtRechitClusters] = tmp.maxChamber;
        	MuonSystem->dtRechitClusterMaxChamberRatio[MuonSystem->nDtRechitClusters] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
        	MuonSystem->dtRechitClusterNChamber[MuonSystem->nDtRechitClusters] = tmp.nChamber;
        	MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters] = tmp.maxStation;
        	MuonSystem->dtRechitClusterMaxStationRatio[MuonSystem->nDtRechitClusters] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
        	MuonSystem->dtRechitClusterNStation[MuonSystem->nDtRechitClusters] = tmp.nStation;
          MuonSystem->dtRechitClusterNStation10[MuonSystem->nDtRechitClusters] = tmp.nStation10;
          MuonSystem->dtRechitClusterAvgStation10[MuonSystem->nDtRechitClusters] = tmp.avgStation10;

          //Jet veto/ muon veto
          MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = 0.0;

          if (!isData)
          {
            for(int i = 0; i < nGenParticle; i++)
            {
	      if (fabs(gParticleId[i])!=13)continue;
              if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && gParticlePt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] ) {
                MuonSystem->dtRechitClusterGenMuonVetoPt[MuonSystem->nDtRechitClusters]  = gParticlePt[i];
              }
              if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.8 && gParticlePt[i] > MuonSystem->dtRechitClusterGenMuonVetoPt_dR0p8[MuonSystem->nDtRechitClusters] ) {
                MuonSystem->dtRechitClusterGenMuonVetoPt_dR0p8[MuonSystem->nDtRechitClusters]  = gParticlePt[i];
              }
            }
          }
          // jet veto
          for(int i = 0; i < nJets; i++)
          {
            if (fabs(jetEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters]  = jetPt[i];
              // MuonSystem->dtRechitClusterJetVetoEta[MuonSystem->nDtRechitClusters]  = jetEta[i];
              // MuonSystem->dtRechitClusterJetVetoPhi[MuonSystem->nDtRechitClusters]  = jetPhi[i];

            }


          }


          for(int i = 0; i < nMuons; i++)
          {
            if (fabs(muonEta[i]>3.0)) continue;
            float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
            if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters]  = muonPt[i];


            }
          }



          for (int i = 0; i < nDtSeg; i++) {
              if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
                if (dtSegStation[i] == 1) MuonSystem->dtRechitClusterNSegStation1[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 2) MuonSystem->dtRechitClusterNSegStation2[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 3) MuonSystem->dtRechitClusterNSegStation3[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 4) MuonSystem->dtRechitClusterNSegStation4[MuonSystem->nDtRechitClusters]  +=1;
              }
              if (abs(RazorAnalyzer::deltaPhi(dtSegPhi[i],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]))>2) {
                if (dtSegStation[i] == 1) MuonSystem->dtRechitClusterNOppositeSegStation1[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 2) MuonSystem->dtRechitClusterNOppositeSegStation2[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 3) MuonSystem->dtRechitClusterNOppositeSegStation3[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 4) MuonSystem->dtRechitClusterNOppositeSegStation4[MuonSystem->nDtRechitClusters]  +=1;
              }
         }

          //match to MB1 DT segments
          MuonSystem->nCscRechits = ncscRechits;

          for (int i = 0; i < nDtRechits; i++) {
            // cout<<"RECHIT: "<<i<<","<<dtRechitEta[i]<<", "<<dtRechitPhi[i]<<", "<<dtRechitStation[i]<<endl;
            // cout<<"Cluster: "<<MuonSystem->nDtRechitClusters<<", "<<MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters]<<", "<<MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]<<endl;
            // cout<<MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters]<<", "<<RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])<<endl;
            if (RazorAnalyzer::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5 )
            {
              if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters]++;
            }
            if (RazorAnalyzer::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
            {
              if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p4[MuonSystem->nDtRechitClusters]++;
            }
            if(abs(dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])==1 && dtRechitStation[i] == 1)
            {
              if (abs(RazorAnalyzer::deltaPhi(dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < TMath::Pi()/4.0 )
              {
                if (dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] == 1)
                {
                  MuonSystem->dtRechitCluster_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters]++;
                }
                else
                {
                  MuonSystem->dtRechitCluster_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters]++;
                }
              }
            }

          }

         std::vector<int> dtRechitCluster_match_rpcBx;

         //match to RPC hits with dPhi<0.5 and same wheel in DT
         for (int i = 0; i < nRpc; i++) {
           float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
           if (rpcRegion[i]!=0) continue;
           if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5 )
           {
             if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
             {
               dtRechitCluster_match_rpcBx.push_back(rpcBx[i]);
               MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]++;
               if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters]++;

             }
           }
           if(RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
           {
             if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_0p4[MuonSystem->nDtRechitClusters]++;
           }
         }


         int max_occurence = 0;
         int max_bx = -999;
         for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++)
         {
           int counter = 0;
           for(unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j++)
           {
             if (dtRechitCluster_match_rpcBx[j] == dtRechitCluster_match_rpcBx[l]) counter++;
           }
           if (counter>max_occurence)
           {
             max_occurence = counter;
             max_bx = dtRechitCluster_match_rpcBx[l];
           }
           // cout<<dtRechitCluster_match_rpcBx[l]<<endl;
         }
         // cout<<max_occurence<<","<<max_bx<<endl;
         // cout<<"###################################################"<<endl;
           MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters] = max_bx;
           MuonSystem->dtRechitClusterMetEENoise_dPhi[MuonSystem->nDtRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhiEENoise);
           MuonSystem->dtRechitCluster_match_gLLP_deltaR[MuonSystem->nDtRechitClusters] = RazorAnalyzer::deltaR(MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters], MuonSystem->gLLP_eta, MuonSystem->gLLP_phi);

           //// ========================== CUTS ===========================

           //if (MuonSystem->dtRechitClusterSize[MuonSystem->nDtRechitClusters] < 50) continue;
           //if (MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters] <= 0) continue;
           //if (MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] > 20) continue;

           MuonSystem->nDtRechitClusters++;
        }

      for(int i = 0; i < MuonSystem->nDtRechitClusters; i++)
      {
        double max_deltaPhi = 0.;
        int index = 999;
        for(int j = 0; j < MuonSystem->nDtRechitClusters; j++)
        {
          // check if there is another dtrechitcluster on the opposite end


          double current_delta_phi = abs(RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[i], MuonSystem->dtRechitClusterPhi[j]));

          if (current_delta_phi > max_deltaPhi)
          {
            max_deltaPhi = current_delta_phi;
            index = j;
          }
        }
        MuonSystem->dtRechitClusterMaxDPhi[i] = max_deltaPhi;
        MuonSystem->dtRechitClusterMaxDPhi_index[i] = index;

      }
      // cout<<MuonSystem->nDtRechitClusters + MuonSystem->nCscRechitClusters<<endl;
      //if ( (MuonSystem->nDtRechitClusters + MuonSystem->nCscRechitClusters) == 0.0)continue;
       if ( isData && MuonSystem->nDtRechitClusters + MuonSystem->nCscRechitClusters < 1)continue;
      //if ( MuonSystem->nDtRechitClusters + MuonSystem->nCscRechitClusters < 1)continue;
      // if ( isData &&  MuonSystem->nCscRechitClusters < 1)continue;
      //if (!isData && !(MuonSystem->gLLP_csc || MuonSystem->gLLP_dt))continue;
      //if ( MuonSystem->nCscRechitClusters < 1)continue;
      // if ( MuonSystem->nDtRechitClusters +  MuonSystem->nDtRechitClusters< 1)continue;

      MuonSystem->tree_->Fill();
    }

    
          
      if (!isData)
      {
         cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
         cout << "Writing output trees..." << endl;
         outFile->cd();
         MuonSystem->tree_->Write();
         NEvents->Write();
         SingleMu->Write();
         MuTrig->Write();
         outFile->Close();
      }

      else
      {
        cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
        cout << "Writing output trees..." << endl;
        outFile->cd();
        MuonSystem->tree_->Write();

        NEvents->Write();
        SingleMu->Write();
        MuTrig->Write();
        // outFile->Write();
        outFile->Close();
      }
}


