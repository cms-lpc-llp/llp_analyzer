#include "llp_MuonSystem_bparking.h"
#include "RazorHelper.h"
#include "TreeMuonSystemBParking.h"
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
  // bool passLooseId;
  // bool passMediumId;
  bool tightId;
  bool looseId;
  bool passVetoId;
  bool passLooseIso;
  bool passTightIso;
  bool passVTightIso;
  bool passVVTightIso;

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


void llp_MuonSystem_bparking::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";



  bool signalScan = int(options/10) == 1;
  int option = options%10;

  if( isData )
  {
    std::cout << "[INFO]: running on data with option: " << option << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with option: " << option << std::endl;
  }
  if( signalScan )
  {
    std::cout << "[INFO]: running with Signal scan" << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running without Signal scan " << option << std::endl;
  }



  const float ELE_MASS = 0.000511;
  const float MU_MASS  = 0.105658;
  const float Z_MASS   = 91.2;

  if (analysisTag == ""){
    analysisTag = "Razor2016_80X";

  }
  int wzId;


  const int zh_lepton0_cut = 15;
  const int zh_lepton1_cut = 15;

  const int wh_muonPt_cut = 25;
  const int wh_elePt_cut = 35;


  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile;
  if (isData || !signalScan) outFile = new TFile(outfilename.c_str(), "RECREATE");


  TreeMuonSystemBParking *MuonSystem = new TreeMuonSystemBParking;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  // for signals, need one output file for each signal point
  map<pair<int,int>, TFile*> Files2D;
  map<pair<int,int>, TTree*>Trees2D;
  map<pair<int,int>, TH1F*> NEvents2D;
  map<pair<int,int>, TH1F*> accep2D;
  map<pair<int,int>, TH1F*> accep_met2D;
  map<pair<int,int>, TH1F*> Total2D;



  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *Total = new TH1F("Total", "Total", 1, 1, 2);

  TH1F *accep = new TH1F("accep", "acceptance", 1, 1, 2);
  TH1F *accep_met = new TH1F("accep_met", "acceptance_met", 1, 1, 2);

  TH1F *Nmet200 = new TH1F("Nmet200", "Nmet200", 1, 1, 2);
  TH1F *NmetFilter = new TH1F("NmetFilter", "NmetFilter", 1, 1, 2);
  TH1F *Nlep0 = new TH1F("Nlep0", "Nlep0", 1, 1, 2);
  TH1F *Njet1 = new TH1F("Njet1", "Njet1", 1, 1, 2);
  TH1F *NcosmicVeto = new TH1F("NcosmicVeto", "NcosmicVeto", 1, 1, 2);



  // TH1F *NEvents_genweight = new TH1F("NEvents_genweight", "NEvents_genweight", 1, 1, 2);

  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
  if(cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1

  // cout << "Getting JEC parameters from " << pathname << endl;
  //
  // std::vector<JetCorrectorParameters> correctionParameters;
  // correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
  // correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
  // correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));
  //
  // FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);




  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  helper = new RazorHelper(analysisTag, isData, false);


  std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
  std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();


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
    // if (jentry>10000) continue;
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
    //std::cout << "deb0 " << jentry << std::endl;
    MuonSystem->InitVariables();
    // std::cout << "deb1 " << jentry << std::endl;




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

    // if (eventNum!=476483009 && eventNum!=346074531 && eventNum!=228964571 && eventNum!=501389779)continue;
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

    std::pair<double,double> corrected_met;
    if (analysisTag=="Razor2016_07Aug2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2016, !isData, nPV);
    else if (analysisTag=="Razor2017_17Nov2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2017, !isData, nPV);
    else if (analysisTag=="Razor2018_17SeptEarlyReReco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2018, !isData, nPV);



    MuonSystem->metEENoise = corrected_met.first;
    MuonSystem->metPhiEENoise = corrected_met.second;


    //Triggers
    for(int i = 0; i < NTriggersMAX; i++){
      MuonSystem->HLTDecision[i] = HLTDecision[i];
    }


    MuonSystem->Flag2_HBHENoiseFilter = Flag2_HBHENoiseFilter;
    MuonSystem->Flag2_HBHEIsoNoiseFilter = Flag2_HBHEIsoNoiseFilter;
    MuonSystem->Flag2_BadPFMuonFilter = Flag2_BadPFMuonFilter;
    MuonSystem->Flag2_globalSuperTightHalo2016Filter = Flag2_globalSuperTightHalo2016Filter;
    MuonSystem->Flag2_globalTightHalo2016Filter = Flag2_globalTightHalo2016Filter;
    MuonSystem->Flag2_BadChargedCandidateFilter =Flag2_BadChargedCandidateFilter;
    // Flag2_goodVertices = Flag2_goodVertices;
    MuonSystem->Flag2_EcalDeadCellTriggerPrimitiveFilter = Flag2_EcalDeadCellTriggerPrimitiveFilter;
    MuonSystem->Flag2_ecalBadCalibFilter = Flag2_ecalBadCalibFilter;
    MuonSystem->Flag2_eeBadScFilter = Flag2_eeBadScFilter;
    MuonSystem->Flag2_all = Flag2_HBHENoiseFilter && Flag2_HBHEIsoNoiseFilter && Flag2_BadPFMuonFilter && Flag2_globalSuperTightHalo2016Filter && Flag2_EcalDeadCellTriggerPrimitiveFilter;
    if (isData) MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_eeBadScFilter;

    if (analysisTag!="Razor2016_07Aug2017Rereco")MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_ecalBadCalibFilter;

    // cout<<nGenParticle<<endl;
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
    }

    if (abs(MuonSystem->gLLP_eta) < 2.4
      && abs(MuonSystem->gLLP_decay_vertex_z)<1100 && abs(MuonSystem->gLLP_decay_vertex_z)>400
      && MuonSystem->gLLP_decay_vertex_r < 695.5) MuonSystem->gLLP_csc = true;
    if (abs(MuonSystem->gLLP_decay_vertex_z)< 661.0
      && MuonSystem->gLLP_decay_vertex_r < 800
       && MuonSystem->gLLP_decay_vertex_r > 200.0) MuonSystem->gLLP_dt = true;


          //*************************************************************************
          //Start Object Selection
          //*************************************************************************

          std::vector<leptons> Leptons;
          //-------------------------------
          //Muons
          //-------------------------------
          for( int i = 0; i < nMuons; i++ )
          {

            // if(!isMuonPOGLooseMuon(i)) continue;
            // if(muonPt[i] < 25) continue;
            if(fabs(muonEta[i]) > 2.4) continue;

            //remove overlaps
            bool overlap = false;
            for(auto& lep : Leptons)
            {
              if (RazorAnalyzer::deltaR(muonEta[i],muonPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
            }
            if(overlap) continue;

            leptons tmpMuon;
            tmpMuon.lepton.SetPtEtaPhiM(muonPt[i],muonEta[i], muonPhi[i], MU_MASS);
            tmpMuon.pdgId = 13 * -1 * muonCharge[i];
            tmpMuon.dZ = muon_dZ[i];
            tmpMuon.looseId = isMuonPOGTightMuon(i);
            tmpMuon.tightId = isMuonPOGLooseMuon(i);


            float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];

            tmpMuon.passLooseIso = muonIso<0.25;
            tmpMuon.passTightIso = muonIso<0.15;
            tmpMuon.passVTightIso = muonIso<0.10;
            tmpMuon.passVVTightIso = muonIso<0.05;

            tmpMuon.passVetoId = false;
            Leptons.push_back(tmpMuon);
          }

          //-------------------------------
          //Electrons
          //-------------------------------
          for( int i = 0; i < nElectrons; i++ )
          {

            // if (!isEGammaPOGLooseElectron(i, true, true, true, "Summer16")) continue;
            // if (!isEGammaPOGLooseElectron(i, true, false, true, "vid")) continue;
            // if(!isEGammaPOGLooseElectron(i, true, true, true, "2017_94X"))continue;
            // if (!isEGammaPOGLooseElectron(i, true, false, true, "vid")) continue;

            // if(elePt[i] < 35) continue;
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
            MuonSystem->lepLooseId[MuonSystem->nLeptons] = tmp.looseId;
            MuonSystem->lepTightId[MuonSystem->nLeptons] = tmp.tightId;
            MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
            MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
            MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
            MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;


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

          if (fabs(thisJet.Eta())> 2.65 && fabs(thisJet.Eta())<3.139 && thisJet.Pt() < 50  && analysisTag == "Razor2017_17Nov2017Rereco")
          {
            MetXCorr_EENoise += thisJet.Px();
            MetYCorr_EENoise += thisJet.Py();
          }
          if (fabs(thisJet.Eta()) >= 3.0)continue;
          // if( thisJet.Pt() < 30 ) continue;//According to the April 1st 2015 AN

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




      TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(MuonSystem->metEENoise, 0, MuonSystem->metPhiEENoise, 0);

      //EENoise
      float PFMetXEENoise   = PFMET.Px() + MetXCorr_EENoise;
      float PFMetYEENoise   = PFMET.Py() + MetYCorr_EENoise;
      MuonSystem->metEENoise    = sqrt( pow(PFMetXEENoise,2) + pow(PFMetYEENoise,2) );
      MuonSystem->metPhiEENoise    = atan(PFMetYEENoise/PFMetXEENoise);
      if  (PFMetXEENoise < 0.0) MuonSystem->metPhiEENoise = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiEENoise,0.0);



      MuonSystem->nDTRechits  = 0;
      for (int i = 0; i < nDtRechits; i++) {

        if (dtRechitY[i]>=0.0) MuonSystem->nDTPositiveYRechits++;
        else MuonSystem->nDTNegativeYRechits++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -2) MuonSystem->nDTRechitsChamberMinus12++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -1) MuonSystem->nDTRechitsChamberMinus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 0) MuonSystem->nDTRechitsChamber10++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 1) MuonSystem->nDTRechitsChamberPlus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 2) MuonSystem->nDTRechitsChamberPlus12++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -2) MuonSystem->nDTRechitsChamberMinus22++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -1) MuonSystem->nDTRechitsChamberMinus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 0) MuonSystem->nDTRechitsChamber20++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 1) MuonSystem->nDTRechitsChamberPlus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 2) MuonSystem->nDTRechitsChamberPlus22++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -2) MuonSystem->nDTRechitsChamberMinus32++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -1) MuonSystem->nDTRechitsChamberMinus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 0) MuonSystem->nDTRechitsChamber30++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 1) MuonSystem->nDTRechitsChamberPlus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 2) MuonSystem->nDTRechitsChamberPlus32++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -2) MuonSystem->nDTRechitsChamberMinus42++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -1) MuonSystem->nDTRechitsChamberMinus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 0) MuonSystem->nDTRechitsChamber40++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 1) MuonSystem->nDTRechitsChamberPlus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 2) MuonSystem->nDTRechitsChamberPlus42++;



        if (dtRechitStation[i] == 1) MuonSystem->nDTRechitsStation1++;
        if (dtRechitStation[i] == 2) MuonSystem->nDTRechitsStation2++;
        if (dtRechitStation[i] == 3) MuonSystem->nDTRechitsStation3++;
        if (dtRechitStation[i] == 4) MuonSystem->nDTRechitsStation4++;

        if (dtRechitWheel[i] == -2) MuonSystem->nDTRechitsWheelMinus2++;
        if (dtRechitWheel[i] == -1) MuonSystem->nDTRechitsWheelMinus1++;
        if (dtRechitWheel[i] == 0) MuonSystem->nDTRechitsWheel0++;
        if (dtRechitWheel[i] == 1) MuonSystem->nDTRechitsWheelPlus1++;
        if (dtRechitWheel[i] == 2) MuonSystem->nDTRechitsWheelPlus2++;


      }


      if ( MuonSystem->nDTRechitsChamberMinus12 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus11 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamber10 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus11 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus12 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus22 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus21 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamber20 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus21 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus22 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus32 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus31 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamber30 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus31 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus32 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus42 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberMinus41 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamber40 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus41 > 50) MuonSystem->nDtRings++;
      if ( MuonSystem->nDTRechitsChamberPlus42 > 50) MuonSystem->nDtRings++;


      if (MuonSystem->nDTRechitsStation1 > 25) MuonSystem->nDtStations25++;
      if (MuonSystem->nDTRechitsStation2 > 25) MuonSystem->nDtStations25++;
      if (MuonSystem->nDTRechitsStation3 > 25) MuonSystem->nDtStations25++;
      if (MuonSystem->nDTRechitsStation4 > 25) MuonSystem->nDtStations25++;

      if (MuonSystem->nDTRechitsWheelMinus2 > 25) MuonSystem->nDtWheels25++;
      if (MuonSystem->nDTRechitsWheelMinus1 > 25) MuonSystem->nDtWheels25++;
      if (MuonSystem->nDTRechitsWheel0 > 25) MuonSystem->nDtWheels25++;
      if (MuonSystem->nDTRechitsWheelPlus1 > 25) MuonSystem->nDtWheels25++;
      if (MuonSystem->nDTRechitsWheelPlus2 > 25) MuonSystem->nDtWheels25++;


      vector<Point> points;
      vector<int> cscRechitsClusterId;
      points.clear();
      MuonSystem->nCscRechits  = 0;

      for (int i = 0; i < ncscRechits; i++) {
        // if (cscRechitsQuality[i]!=1) continue;
        // if (cscRechitsQuality[i]>2) cout<<cscRechitsQuality[i]<<endl;
        // MuonSystem->cscRechitsPhi[MuonSystem->nCscRechits]           = cscRechitsPhi[i];   //[nCsc]
        // MuonSystem->cscRechitsEta[MuonSystem->nCscRechits]           = cscRechitsEta[i];   //[nCsc]
        // MuonSystem->cscRechitsX[MuonSystem->nCscRechits]             = cscRechitsX[i];   //[nCsc]
        // MuonSystem->cscRechitsY[MuonSystem->nCscRechits]             = cscRechitsY[i];   //[nCsc]
        // MuonSystem->cscRechitsZ[MuonSystem->nCscRechits]             = cscRechitsZ[i];   //[nCsc]
        // MuonSystem->cscRechitsTpeak[MuonSystem->nCscRechits] = cscRechitsTpeak[i];
        // MuonSystem->cscRechitsTwire[MuonSystem->nCscRechits] = cscRechitsTwire[i];
        // MuonSystem->cscRechitsQuality[MuonSystem->nCscRechits] = cscRechitsQuality[i];
        // MuonSystem->cscRechitsStation[MuonSystem->nCscRechits] = cscRechitsStation[i];
        // MuonSystem->cscRechitsChamber[MuonSystem->nCscRechits] = cscRechitsChamber[i];
        // static int station(int index) { return ((index >> START_STATION) & MASK_STATION); }
        // cout<<cscRechitsStation[i]<<", " << ((cscRechitsDetId[i] >> 12) & 07)<<","<<((cscRechitsDetId[i] >> 3) & 077)<<endl;

        // MuonSystem->cscRechitsChamber[MuonSystem->nCscRechits] = ((index >> START_CHAMBER) & MASK_CHAMBER);
        //pick out the right bits for chamber
        int chamber = ((cscRechitsDetId[i] >> 3) & 077); //https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonDetId/interface/CSCDetId.h#L147

        int layer = (cscRechitsDetId[i] & 07);

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
        if (cscRechitsY[i]>=0.0)
        {
          MuonSystem->nCscPositiveYRechits++;
          MuonSystem->cscPosTpeak = MuonSystem->cscPosTpeak + cscRechitsTpeak[i];
        }
        else
        {
          MuonSystem->nCscNegativeYRechits++;
          MuonSystem->cscNegTpeak = MuonSystem->cscNegTpeak + cscRechitsTpeak[i];
        }
        if (cscRechitsTpeak[i]<-12.5)MuonSystem->nEarlyCscRechits++;
        if (cscRechitsTpeak[i]>12.5)MuonSystem->nLateCscRechits++;
        if (cscRechitsTpeak[i]<-25)MuonSystem->nEarly2CscRechits++;
        if (cscRechitsTpeak[i]>25)MuonSystem->nLate2CscRechits++;
        if (cscRechitsChamber[i] == 11) MuonSystem->nCscRechitsChamberPlus11++;
        if (cscRechitsChamber[i] == 12) MuonSystem->nCscRechitsChamberPlus12++;
        if (cscRechitsChamber[i] == 13) MuonSystem->nCscRechitsChamberPlus13++;
        if (cscRechitsChamber[i] == 21) MuonSystem->nCscRechitsChamberPlus21++;
        if (cscRechitsChamber[i] == 22) MuonSystem->nCscRechitsChamberPlus22++;
        if (cscRechitsChamber[i] == 31) MuonSystem->nCscRechitsChamberPlus31++;
        if (cscRechitsChamber[i] == 32) MuonSystem->nCscRechitsChamberPlus32++;
        if (cscRechitsChamber[i] == 41) MuonSystem->nCscRechitsChamberPlus41++;
        if (cscRechitsChamber[i] == 42) MuonSystem->nCscRechitsChamberPlus42++;
        if (cscRechitsChamber[i] == -11) MuonSystem->nCscRechitsChamberMinus11++;
        if (cscRechitsChamber[i] == -12) MuonSystem->nCscRechitsChamberMinus12++;
        if (cscRechitsChamber[i] == -13) MuonSystem->nCscRechitsChamberMinus13++;
        if (cscRechitsChamber[i] == -21) MuonSystem->nCscRechitsChamberMinus21++;
        if (cscRechitsChamber[i] == -22) MuonSystem->nCscRechitsChamberMinus22++;
        if (cscRechitsChamber[i] == -31) MuonSystem->nCscRechitsChamberMinus31++;
        if (cscRechitsChamber[i] == -32) MuonSystem->nCscRechitsChamberMinus32++;
        if (cscRechitsChamber[i] == -41) MuonSystem->nCscRechitsChamberMinus41++;
        if (cscRechitsChamber[i] == -42) MuonSystem->nCscRechitsChamberMinus42++;
        // MuonSystem->nCscRechits++;
      }
      MuonSystem->nCscRings = 0;
      if ( MuonSystem->nCscRechitsChamberPlus11 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberPlus12 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberPlus13 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberPlus21 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberPlus22 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberPlus31 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberPlus32 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberPlus41 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberPlus42 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberMinus11 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberMinus12 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberMinus13 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberMinus21 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberMinus22 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberMinus31 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberMinus32 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberMinus41 > 50) MuonSystem->nCscRings++;
      if ( MuonSystem->nCscRechitsChamberMinus42 > 50) MuonSystem->nCscRings++;
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
          MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] =tmp.y;
          MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] =tmp.z;
          MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] = tmp.t;
          MuonSystem->cscRechitClusterTimeTotal[MuonSystem->nCscRechitClusters] = tmp.tTotal;
          MuonSystem->cscRechitClusterTimeWeighted[MuonSystem->nCscRechitClusters] = tmp.tWeighted;
          MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] =tmp.eta;
          MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters] = tmp.phi;

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


          MuonSystem->cscRechitClusterMe11Ratio[MuonSystem->nCscRechitClusters] = tmp.Me11Ratio;
          MuonSystem->cscRechitClusterMe12Ratio[MuonSystem->nCscRechitClusters] = tmp.Me12Ratio;


          // MuonSystem->cscRechitClusterNLayersChamberPlus11[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus11;
          // MuonSystem->cscRechitClusterNLayersChamberPlus12[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus12;
          // MuonSystem->cscRechitClusterNLayersChamberPlus13[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus13;
          // MuonSystem->cscRechitClusterNLayersChamberPlus21[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus21;
          // MuonSystem->cscRechitClusterNLayersChamberPlus22[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus22;
          // MuonSystem->cscRechitClusterNLayersChamberPlus31[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus31;
          // MuonSystem->cscRechitClusterNLayersChamberPlus32[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus32;
          // MuonSystem->cscRechitClusterNLayersChamberPlus41[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus41;
          // MuonSystem->cscRechitClusterNLayersChamberPlus42[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus42;
          // MuonSystem->cscRechitClusterNLayersChamberMinus11[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus11;
          // MuonSystem->cscRechitClusterNLayersChamberMinus12[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus12;
          // MuonSystem->cscRechitClusterNLayersChamberMinus13[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus13;
          // MuonSystem->cscRechitClusterNLayersChamberMinus21[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus21;
          // MuonSystem->cscRechitClusterNLayersChamberMinus22[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus22;
          // MuonSystem->cscRechitClusterNLayersChamberMinus31[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus31;
          // MuonSystem->cscRechitClusterNLayersChamberMinus32[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus32;
          // MuonSystem->cscRechitClusterNLayersChamberMinus41[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus41;
          // MuonSystem->cscRechitClusterNLayersChamberMinus42[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus42;



          //Jet veto/ muon veto
          MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;

          // jet veto
          for(int i = 0; i < nJets; i++)
          {
            if (fabs(jetEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]  = jetPt[i];
              MuonSystem->cscRechitClusterJetVetoEta[MuonSystem->nCscRechitClusters]  = jetEta[i];
              MuonSystem->cscRechitClusterJetVetoPhi[MuonSystem->nCscRechitClusters]  = jetPhi[i];
            }
            if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetE[i] > MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters]  = jetE[i];
            }

          }
          for(int i = 0; i < nMuons; i++)
          {
            if (fabs(muonEta[i]>3.0)) continue;
            float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
            if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]  = muonPt[i];
              MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters]  = muonE[i];
              MuonSystem->cscRechitClusterMuonVetoPhi[MuonSystem->nCscRechitClusters]  = muonPhi[i];
              MuonSystem->cscRechitClusterMuonVetoEta[MuonSystem->nCscRechitClusters]  = muonEta[i];
              MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters]  = muon_isGlobal[i];
              MuonSystem->cscRechitClusterMuonVetoLooseId[MuonSystem->nCscRechitClusters]  = muonIsLoose[i];


            }
          }


          //match to MB1 DT segments
          for (int i = 0; i < nDtSeg; i++) {
            if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
            {
              MuonSystem->cscRechitCluster_match_dtSeg_0p4[MuonSystem->nCscRechitClusters] ++;
              if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters] ++;
            }



          }
          //match to RPC hits in RE1/2
          for (int i = 0; i < nRpc; i++) {
            float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
            if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
            {
              if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)
              {
                MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters] ++;
              }
              if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)
              {
                MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters] ++;
              }

            }
          }


          MuonSystem->cscRechitCluster_match_gLLP_deltaR[MuonSystem->nCscRechitClusters] = RazorAnalyzer::deltaR(MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->gLLP_eta, MuonSystem->gLLP_phi);





          MuonSystem->cscRechitClusterMetEENoise_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhiEENoise);


          MuonSystem->nCscRechitClusters++;
      }



      // DT cluster

      points.clear();
      for (int i = 0; i < nDtRechits; i++) {
        Point p;
        p.phi = dtRechitPhi[i];
        p.eta = dtRechitEta[i];
        p.x = dtRechitX[i];
        p.y = dtRechitY[i];
        p.z = dtRechitZ[i];
        p.t = dtRechitTime[i];
        p.twire = dtRechitTime[i];
        p.station = dtRechitStation[i];
        p.chamber = dtRechitWheel[i];
        p.superlayer = 0;

        p.phi = dtRechitCorrectPhi[i];
        p.eta = dtRechitCorrectEta[i];
        p.x = dtRechitCorrectX[i];
        p.y = dtRechitCorrectY[i];
        p.z = dtRechitCorrectZ[i];
        p.superlayer = dtRechitSuperLayer[i];
        // Int_t           dtRechitLayer[20000];   //[nDtRechits]
        // Float_t         dtRechitCorrectX[20000];   //[nDtRechits]
        // Float_t         dtRechitCorrectY[20000];   //[nDtRechits]
        // Float_t         dtRechitCorrectZ[20000];   //[nDtRechits]
        // Float_t         dtRechitCorrectEta[20000];   //[nDtRechits]
        // Float_t         dtRechitCorrectPhi[20000];   //[nDtRechits]

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
          if (abs(tmp.z) < 126.8) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 0;
          else if (tmp.z > 126.8 && tmp.z < 395.4) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 1;
          else if (tmp.z < -126.8 && tmp.z > -395.4)MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -1;
          else if (tmp.z<0) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -2;
          else MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 2;


          MuonSystem->dtRechitClusterTime[MuonSystem->nDtRechitClusters] =tmp.t;
          MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters] =tmp.eta;
          MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters] =tmp.phi;
          MuonSystem->dtRechitClusterTimeSpread[MuonSystem->nDtRechitClusters] =tmp.TSpread;


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
          MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = 0.0;


          // jet veto
          for(int i = 0; i < nJets; i++)
          {
            if (fabs(jetEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters]  = jetPt[i];
              MuonSystem->dtRechitClusterJetVetoEta[MuonSystem->nDtRechitClusters]  = jetEta[i];
              MuonSystem->dtRechitClusterJetVetoPhi[MuonSystem->nDtRechitClusters]  = jetPhi[i];

            }
            if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetE[i] > MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters]  = jetE[i];
            }


          }


          for(int i = 0; i < nMuons; i++)
          {
            if (fabs(muonEta[i]>3.0)) continue;
            float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
            if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters]  = muonPt[i];
              MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters]  = muonE[i];
              MuonSystem->dtRechitClusterMuonVetoPhi[MuonSystem->nDtRechitClusters]  = muonPhi[i];
              MuonSystem->dtRechitClusterMuonVetoEta[MuonSystem->nDtRechitClusters]  = muonEta[i];
              MuonSystem->dtRechitClusterMuonVetoGlobal[MuonSystem->nDtRechitClusters]  = muon_isGlobal[i];
              MuonSystem->dtRechitClusterMuonVetoLooseId[MuonSystem->nDtRechitClusters]  = muonIsLoose[i];

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
            if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5 )
            {
              if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters] ++;
            }
            if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
            {
              if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p4[MuonSystem->nDtRechitClusters] ++;
            }
            if(abs(dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])==1 && dtRechitStation[i] == 1)
            {
              if (abs(RazorAnalyzer::deltaPhi(dtRechitPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < TMath::Pi()/4.0 )
              {
                if (dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] == 1)
                {
                  MuonSystem->dtRechitCluster_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters] ++;
                }
                else
                {
                  MuonSystem->dtRechitCluster_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters] ++;
                }
              }
            }



          }

         //  MuonSystem->nRpc = nRpc;
         //  //match to RPC hits with dPhi<0.5 and same wheel in DT
         //  for (int i = 0; i < nRpc; i++) {
         //    float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
         //    if (rpcRegion[i]!=0) continue;
         //    if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5 )
         //    {
         //      if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
         //      {
         //        MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters] ++;
         //        MuonSystem->dtRechitCluster_match_RPCTime_dPhi0p5[MuonSystem->nDtRechitClusters] += rpcBx[i];
         //
         //        MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcBx[i]);
         //        MuonSystem->dtRechitCluster_match_RPCX_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcX[i]);
         //        MuonSystem->dtRechitCluster_match_RPCY_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcY[i]);
         //        MuonSystem->dtRechitCluster_match_RPCZ_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcZ[i]);
         //        MuonSystem->dtRechitCluster_match_RPCPhi_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcPhi[i]);
         //        MuonSystem->dtRechitCluster_match_RPCEta_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcEta[i]);
         //        MuonSystem->dtRechitCluster_match_RPCRing_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcRing[i]);
         //        MuonSystem->dtRechitCluster_match_RPCLayer_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcLayer[i]);
         //        MuonSystem->dtRechitCluster_match_RPCSector_dPhi0p5[MuonSystem->nDtRechitClusters].push_back(rpcSector[i]);
         //
         //      }
         //
         //    }
         //    if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
         //    {
         //        MuonSystem->dtRechitCluster_match_RPChits_dR0p4[MuonSystem->nDtRechitClusters] ++;
         //        MuonSystem->dtRechitCluster_match_RPCTime_dR0p4[MuonSystem->nDtRechitClusters]  += rpcBx[i];
         //        MuonSystem->dtRechitCluster_match_RPCBx_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcBx[i]);
         //        MuonSystem->dtRechitCluster_match_RPCX_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcX[i]);
         //        MuonSystem->dtRechitCluster_match_RPCY_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcY[i]);
         //        MuonSystem->dtRechitCluster_match_RPCZ_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcZ[i]);
         //        MuonSystem->dtRechitCluster_match_RPCPhi_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcPhi[i]);
         //        MuonSystem->dtRechitCluster_match_RPCEta_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcEta[i]);
         //        MuonSystem->dtRechitCluster_match_RPCRing_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcRing[i]);
         //        MuonSystem->dtRechitCluster_match_RPCLayer_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcLayer[i]);
         //        MuonSystem->dtRechitCluster_match_RPCSector_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcSector[i]);
         //        if (rpcStation[i] == MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters])
         //        {
         //          MuonSystem->dtRechitCluster_match_RPChits_sameStation_dR0p4[MuonSystem->nDtRechitClusters] ++;
         //          MuonSystem->dtRechitCluster_match_RPCTime_sameStation_dR0p4[MuonSystem->nDtRechitClusters]  += rpcBx[i];
         //          MuonSystem->dtRechitCluster_match_RPCBx_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcBx[i]);
         //          MuonSystem->dtRechitCluster_match_RPCX_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcX[i]);
         //          MuonSystem->dtRechitCluster_match_RPCY_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcY[i]);
         //          MuonSystem->dtRechitCluster_match_RPCZ_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcZ[i]);
         //          MuonSystem->dtRechitCluster_match_RPCPhi_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcPhi[i]);
         //          MuonSystem->dtRechitCluster_match_RPCEta_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcEta[i]);
         //          MuonSystem->dtRechitCluster_match_RPCRing_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcRing[i]);
         //          MuonSystem->dtRechitCluster_match_RPCLayer_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcLayer[i]);
         //          MuonSystem->dtRechitCluster_match_RPCSector_sameStation_dR0p4[MuonSystem->nDtRechitClusters].push_back(rpcSector[i]);
         //        }
         //    }
         //
         //  }
         // MuonSystem->dtRechitCluster_match_RPCTime_dPhi0p5[MuonSystem->nDtRechitClusters] /= 1.0* MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters];
         // MuonSystem->dtRechitCluster_match_RPCTime_dR0p4[MuonSystem->nDtRechitClusters] /= 1.0* MuonSystem->dtRechitCluster_match_RPChits_dR0p4[MuonSystem->nDtRechitClusters];
         // MuonSystem->dtRechitCluster_match_RPCTime_sameStation_dR0p4[MuonSystem->nDtRechitClusters] /= 1.0* MuonSystem->dtRechitCluster_match_RPChits_sameStation_dR0p4[MuonSystem->nDtRechitClusters];
         // for (int i = 0; i < nRpc; i++) {
         //   if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5 )
         //   {
         //     if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
         //     {
         //       MuonSystem->dtRechitCluster_match_RPCTimeSpread_dPhi0p5[MuonSystem->nDtRechitClusters] += pow(rpcBx[i] - MuonSystem->dtRechitCluster_match_RPCTime_dPhi0p5[MuonSystem->nDtRechitClusters],2);
         //     }
         //
         //   }
         //   if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
         //   {
         //       MuonSystem->dtRechitCluster_match_RPCTimeSpread_dR0p4[MuonSystem->nDtRechitClusters]  += pow(rpcBx[i] - MuonSystem->dtRechitCluster_match_RPCTime_dR0p4[MuonSystem->nDtRechitClusters],2);
         //       if (rpcStation[i] == MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters])
         //       {
         //         MuonSystem->dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[MuonSystem->nDtRechitClusters]  += pow(rpcBx[i] - MuonSystem->dtRechitCluster_match_RPCTime_sameStation_dR0p4[MuonSystem->nDtRechitClusters],2);
         //       }
         //   }
         //
         // }
         //
         // MuonSystem->dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[MuonSystem->nDtRechitClusters] = sqrt(MuonSystem->dtRechitCluster_match_RPCTimeSpread_sameStation_dR0p4[MuonSystem->nDtRechitClusters]/MuonSystem->dtRechitCluster_match_RPChits_sameStation_dR0p4[MuonSystem->nDtRechitClusters]);
         // MuonSystem->dtRechitCluster_match_RPCTimeSpread_dR0p4[MuonSystem->nDtRechitClusters] = sqrt(MuonSystem->dtRechitCluster_match_RPCTimeSpread_dR0p4[MuonSystem->nDtRechitClusters]/MuonSystem->dtRechitCluster_match_RPChits_dR0p4[MuonSystem->nDtRechitClusters]);
         // MuonSystem->dtRechitCluster_match_RPCTimeSpread_dPhi0p5[MuonSystem->nDtRechitClusters] = sqrt(MuonSystem->dtRechitCluster_match_RPCTimeSpread_dPhi0p5[MuonSystem->nDtRechitClusters]/MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]);
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
               if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters] ++;

             }
           }
           if(RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
           {
             if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_0p4[MuonSystem->nDtRechitClusters] ++;
           }
         }
         int max_occurence = 0;
         int max_bx = -999;
         for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++)
         {
           int counter = 0;
           for(unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j ++)
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
      if ( isData && MuonSystem->nDtRechitClusters + MuonSystem->nCscRechitClusters < 1)continue;
      //if (!isData && !(MuonSystem->gLLP_csc || MuonSystem->gLLP_dt))continue;
      //if ( MuonSystem->nCscRechitClusters < 1)continue;
      // if ( MuonSystem->nDtRechitClusters +  MuonSystem->nDtRechitClusters< 1)continue;



      MuonSystem->tree_->Fill();



    }
      if(!isData && signalScan)
      {
        for(auto &filePtr : Files2D)
         {
           cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
           filePtr.second->cd();
           Trees2D[filePtr.first]->Write();
           NEvents2D[filePtr.first]->Write("NEvents");
           Total2D[filePtr.first]->Write("Total");
           accep2D[filePtr.first]->Write("acceptance");
           accep_met2D[filePtr.first]->Write("acceptance_met");
           filePtr.second->Close();

         }
      }
      else if (!isData)
      {
         cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
         cout << "Writing output trees..." << endl;
         outFile->cd();
         MuonSystem->tree_->Write();
         NEvents->Write();
         accep->Write("acceptance");
         accep_met->Write("acceptance_met");
         outFile->Close();
      }


      else
      {
        cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
        cout << "Writing output trees..." << endl;
        outFile->cd();
        MuonSystem->tree_->Write();
        Nmet200->Write();
        NmetFilter->Write();
        Nlep0->Write();
        Njet1->Write();
        NcosmicVeto->Write();
        NEvents->Write();
        // outFile->Write();
        outFile->Close();
      }
}
