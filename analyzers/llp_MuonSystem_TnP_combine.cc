#include "llp_MuonSystem_TnP_combine.h"
#include "RazorHelper.h"
#include "TreeMuonSystemCombination.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
#include "DBSCAN.h"
//C++ includes
#include "assert.h"
#include <random>

//ROOT includes
#include "TH1F.h"


using namespace std::chrono;
using namespace std;

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
  bool isGlobal;
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



void llp_MuonSystem_TnP_combine::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  //options format: MH/MX/ctau/condor: 1000/300/0/1
  // mh can be 3-4 digits, mx is always 3 digits, ctau is one digit(number of zeros), last digit is condor option
  // mh can be 3-4 digits, mx is always 3 digits, ctau is 2 digit(number of zeros), last digit is condor option
  //
  //
  // int mx = int(options/1000)%1000;
  // int mh = options/1000000;
  // int ctau = pow(10, int(options/10)%10) * int(int(options/100)%10);
  //
  // cout<<"mh "<<mh<<", mx "<<mx<<", ctau "<<ctau<<endl;

  bool signalScan = int(options/10) == 1;
  int option = options%10;
  // if (options % 1){
  //   option = 1; // used when running condor
  // }
  // else{
  //   option = 0;// used when running locally
  // }

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


  TreeMuonSystemCombination *MuonSystem = new TreeMuonSystemCombination;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  // for signals, need one output file for each signal point
  map<pair<int,int>, TFile*> Files2D;
  map<pair<int,int>, TTree*> Trees2D;
  map<pair<int,int>, TH1F*> NEvents2D;


  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  // TH1F *NEvents_genweight = new TH1F("NEvents_genweight", "NEvents_genweight", 1, 1, 2);

  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
  if(cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1



  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  helper = new RazorHelper(analysisTag, isData, false);



  std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
  std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();
  //----------
  //pu histo
  //----------

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
    nb = fChain->GetEntry(jentry); nbytes += nb;

    MuonSystem->InitVariables();


    if (!isData && signalScan)
    {

      string mh_substring = lheComments->substr(lheComments->find("MH-")+3);
      int mh = stoi(mh_substring.substr(0,mh_substring.find('_')));
      string mx_substring = lheComments->substr(lheComments->find("MS-")+3);
      int mx = stoi(mx_substring.substr(0,mx_substring.find('_')));
      string ctau_substring = lheComments->substr(lheComments->find("ctauS-")+6);
      int ctau = stoi(ctau_substring.substr(0,ctau_substring.find('_')));
      MuonSystem->mH = mh;
      MuonSystem->mX = mx;
      MuonSystem->ctau = ctau;


      pair<int,int> signalPair = make_pair(mx, ctau);

      if (Files2D.count(signalPair) == 0){ //create file and tree
        //format file name
        string thisFileName = outfilename;
        thisFileName.erase(thisFileName.end()-5, thisFileName.end());
        thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";

        Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
        Trees2D[signalPair] =  MuonSystem->tree_->CloneTree(0);
        NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1,0.5,1.5);


        cout << "Created new output file " << thisFileName << endl;
      }
      //Fill NEvents hist
      NEvents2D[signalPair]->Fill(1.0, genWeight);



    }

    //event info
    if (isData)
    {
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    }
    else
    {
      // cout<<*lheComments<<endl;

      MuonSystem->weight = genWeight;
      NEvents->Fill(1, genWeight);
      // NEvents_genweight->Fill(1);
    }
    MuonSystem->runNum = runNum;
    MuonSystem->lumiSec = lumiNum;
    MuonSystem->evtNum = eventNum;

    bool wzFlag = false;
    // cout<<"ngenparticles: "<<nGenParticle<<endl;
    if (!isData)
    {
      for (int i=0; i < nGenParticle; i++)
      {

        if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == wzId)
        { // choosing only the W->munu events
          wzFlag = true;
          MuonSystem->gLepId = gParticleId[i];
          MuonSystem->gLepPt = gParticlePt[i];
          MuonSystem->gLepEta = gParticleEta[i];
          MuonSystem->gLepE = gParticleE[i];
          MuonSystem->gLepPhi = gParticlePhi[i];
        }
        else if (abs(gParticleId[i]) == 15 && gParticleStatus[i] == 2 && abs(gParticleMotherId[i]) == wzId){
          wzFlag = true;
          MuonSystem->gLepId = gParticleId[i];
          MuonSystem->gLepPt = gParticlePt[i];
          MuonSystem->gLepEta = gParticleEta[i];
          MuonSystem->gLepE = gParticleE[i];
          MuonSystem->gLepPhi = gParticlePhi[i];
        }
        if (abs(gParticleId[i])== 25 || abs(gParticleId[i] == 35))
        {
	        MuonSystem->gHiggsPt = gParticlePt[i];
          MuonSystem->gHiggsEta = gParticleEta[i];
          MuonSystem->gHiggsPhi = gParticlePhi[i];
          MuonSystem->gHiggsE = gParticleE[i];

        }
        MuonSystem->gParticleStatus[MuonSystem->nGenParticle] = gParticleStatus[i];
        MuonSystem->gParticleId[MuonSystem->nGenParticle]  = gParticleId[i];
        MuonSystem->gParticleMotherId[MuonSystem->nGenParticle]  = gParticleMotherId[i];
        MuonSystem->gParticlePt[MuonSystem->nGenParticle]  = gParticlePt[i];
        MuonSystem->gParticleEta[MuonSystem->nGenParticle]  = gParticleEta[i];
        MuonSystem->gParticlePhi[MuonSystem->nGenParticle]  = gParticlePhi[i];
        MuonSystem->gParticleE[MuonSystem->nGenParticle]  = gParticleE[i];
        // MuonSystem->nGenParticle++;
        // cout<<"genparticles: "<<MuonSystem->nGenParticle<<endl;


      }



      MuonSystem->genMetPtTrue = genMetPtTrue;
      MuonSystem->genMetPhiTrue = genMetPhiTrue;
      MuonSystem->genMetPtCalo = genMetPtCalo;
      MuonSystem->genMetPhiCalo = genMetPhiCalo;

      for (int i=0; i < nBunchXing; i++)
      {
        if (BunchXing[i] == 0)
        {
          MuonSystem->npu = nPUmean[i];
        }
      }
      MuonSystem->pileupWeight = helper->getPileupWeight(MuonSystem->npu);
      MuonSystem->pileupWeightUp = helper->getPileupWeightUp(MuonSystem->npu) / MuonSystem->pileupWeight;
      MuonSystem->pileupWeightDown = helper->getPileupWeightDown(MuonSystem->npu) / MuonSystem->pileupWeight;
    }

    //get NPU
    MuonSystem->npv = nPV;
    MuonSystem->rho = fixedGridRhoFastjetAll;
    MuonSystem->met = metType1Pt;
    MuonSystem->metPhi = metType1Phi;
    MuonSystem->metJESUp = MuonSystem->met;
    MuonSystem->metJESDown = MuonSystem->met;
    std::pair<double,double> corrected_met;
    if (analysisTag=="Razor2016_07Aug2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2016, !isData, nPV);
    else if (analysisTag=="Razor2017_17Nov2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2017, !isData, nPV);
    else if (analysisTag=="Razor2018_17SeptEarlyReReco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2018, !isData, nPV);

    MuonSystem->metXYCorr = corrected_met.first;
    MuonSystem->metPhiXYCorr = corrected_met.second;


    // if (MuonSystem->met < 200) continue;
    // if (nCscRechitClusters==0) continue;


    //Triggers
    for(int i = 0; i < NTriggersMAX; i++){
      MuonSystem->HLTDecision[i] = HLTDecision[i];
    }
    // flags
    MuonSystem->Flag2_all = Flag2_HBHENoiseFilter && Flag2_HBHEIsoNoiseFilter && Flag2_BadPFMuonFilter && Flag2_globalSuperTightHalo2016Filter && Flag2_EcalDeadCellTriggerPrimitiveFilter;
    if (isData) MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_eeBadScFilter;

    if (analysisTag!="Razor2016_07Aug2017Rereco")
    {
      MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_ecalBadCalibFilter;
    }

    //*************************************************************************
    //Start Object Selection
    //*************************************************************************

    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
    for( int i = 0; i < nMuons; i++ )
    {
      // if(!muon_isGlobal[i])continue;
      // MuonSystem->GlobalMuonPt[MuonSystem->nGlobalMuons] = muonPt[i];
      // MuonSystem->GlobalMuonEta[MuonSystem->nGlobalMuons] = muonEta[i];
      // MuonSystem->GlobalMuonPhi[MuonSystem->nGlobalMuons] = muonPhi[i];
      // MuonSystem->GlobalMuonLooseId[MuonSystem->nGlobalMuons] = muonIsLoose[i];
      //
      // MuonSystem->nGlobalMuons++;


      // if(!isMuonPOGLooseMuon(i)) continue;
      // if(muonPt[i] < 50) continue;
      // if(fabs(muonEta[i]) > 2.4) continue;
      float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];

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
      tmpMuon.tightId = isMuonPOGTightMuon(i, true, false);
      tmpMuon.looseId = isMuonPOGLooseMuon(i);
      tmpMuon.isGlobal = muon_isGlobal[i];
      tmpMuon.passLooseIso = muonIso<0.25;
      tmpMuon.passTightIso = muonIso<0.15;
      tmpMuon.passVTightIso = muonIso<0.10;
      tmpMuon.passVVTightIso = muonIso<0.05;


      // tmpMuon.passId = isLooseMuon(i);
      tmpMuon.passVetoId = false;
      // tmpMuon.passVetoPOGId = false;

      Leptons.push_back(tmpMuon);
    }

    //-------------------------------
    //Electrons
    //-------------------------------
    for( int i = 0; i < nElectrons; i++ )
    {

      if (!isEGammaPOGLooseElectron(i, true, true, true, "Summer16")) continue;
      if(elePt[i] < zh_lepton1_cut) continue;
      if(fabs(eleEta[i]) > 2.4) continue;

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
      // tmpElectron.passId = isEGammaPOGTightElectron(i, true, true, true, "Summer16");
      // tmpElectron.passVetoId = isEGammaPOGVetoElectron(i, true, true, true, "Summer16");
      Leptons.push_back(tmpElectron);
    }

    sort(Leptons.begin(), Leptons.end(), my_largest_pt);


    //----------------
    //Find Z Candidate
    //----------------

    double ZMass = -999;
    double ZPt = -999;
    double tmpDistToZPole = 9999;
    pair<uint,uint> ZCandidateLeptonIndex;
    bool foundZ = false;
    TLorentzVector ZCandidate;
    double leadingLepPt = 0.0;
    for( uint i = 0; i < Leptons.size(); i++ )
    {
      if(!Leptons[i].looseId) continue;
      if(Leptons[i].lepton.Pt() < 50) continue;
      if(fabs(Leptons[i].lepton.Eta()) > 2.4) continue;

      for( uint j = i+1; j < Leptons.size(); j++ )
      {
        if(!Leptons[j].looseId) continue;
        if(Leptons[j].lepton.Pt() < 50) continue;
        if(fabs(Leptons[j].lepton.Eta()) > 2.4) continue;

        if (!( Leptons[i].pdgId == -1*Leptons[j].pdgId )) continue;// same flavor opposite charge
        double tmpMass = (Leptons[i].lepton+Leptons[j].lepton).M();
        //select the pair closest to Z pole mass
        if ( fabs( tmpMass - Z_MASS) < tmpDistToZPole)
        {
          tmpDistToZPole = tmpMass;
          if (Leptons[i].pdgId > 0)
          {
            ZCandidateLeptonIndex = pair<int,int>(i,j);
          }
          else
          {
            ZCandidateLeptonIndex = pair<int,int>(j,i);
          }
          ZMass = tmpMass;
          ZPt = (Leptons[i].lepton+Leptons[j].lepton).Pt();
          ZCandidate = Leptons[i].lepton+Leptons[j].lepton;
          leadingLepPt = max(Leptons[i].lepton.Pt(),Leptons[j].lepton.Pt());
          foundZ = true;
        }
      }
    }

    // if (foundZ  && Leptons.size() == 2 && leadingLepPt > zh_lepton0_cut)
    if (foundZ)
    {
      MuonSystem->ZMass = ZMass;
      MuonSystem->ZPt   = ZPt;
      MuonSystem->ZEta  = ZCandidate.Eta();
      MuonSystem->ZPhi  = ZCandidate.Phi();
      MuonSystem->ZleptonIndex1 = ZCandidateLeptonIndex.first;
      MuonSystem->ZleptonIndex2 = ZCandidateLeptonIndex.second;
      MuonSystem->category = 2;
    } // endif foundZ
    else{
      continue;
    }
    bool tag = false;
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
      MuonSystem->lepGlobal[MuonSystem->nLeptons] = tmp.isGlobal;
      MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
      MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
      MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
      MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
      if (MuonSystem->nLeptons == MuonSystem->ZleptonIndex1 || MuonSystem->nLeptons == MuonSystem->ZleptonIndex2) MuonSystem->lepFromZ[MuonSystem->nLeptons] = true;
      else MuonSystem->lepFromZ[MuonSystem->nLeptons] = false;
      // if (!isData)
      // {
      //   MuonSystem->lepTriggerSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), true, true, true);
      //   MuonSystem->lepTightIdSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, false);
      //   MuonSystem->lepLooseIdSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, true);
      //   MuonSystem->lepTightIsoSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, false);
      //   MuonSystem->lepLooseIsoSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, true);
      //   MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), true, true, true);
      //   MuonSystem->lepTightIdMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, false);
      //   MuonSystem->lepLooseIdMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, true);
      //   MuonSystem->lepTightIsoMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, false);
      //   MuonSystem->lepLooseIsoMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, true);
      // }
      //
      // if (MuonSystem->lepPassTightIso[MuonSystem->nLeptons] && MuonSystem->lepPassId[MuonSystem->nLeptons] && !tag)
      // {
      //   MuonSystem->lepTag[MuonSystem->nLeptons] = true;
      //   if(!isData)
      //   {
      //     MuonSystem->lepSF[MuonSystem->nLeptons] = MuonSystem->lepTriggerSF[MuonSystem->nLeptons] * MuonSystem->lepTightIdSF[MuonSystem->nLeptons] * MuonSystem->lepTightIsoSF[MuonSystem->nLeptons];
      //     MuonSystem->lepEff[MuonSystem->nLeptons] = MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepTightIdMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepTightIsoMCEfficiency[MuonSystem->nLeptons];
      //   }
      //   else{
      //     MuonSystem->lepSF[MuonSystem->nLeptons] = 1.0;
      //     MuonSystem->lepEff[MuonSystem->nLeptons] = 1.0;
      //   }
      //   tag = true;
      // }
      // else
      // {
      //   MuonSystem->lepTag[MuonSystem->nLeptons] = false;
      //   if(!isData)
      //   {
      //     MuonSystem->lepSF[MuonSystem->nLeptons] = MuonSystem->lepTriggerSF[MuonSystem->nLeptons] * MuonSystem->lepLooseIdSF[MuonSystem->nLeptons] * MuonSystem->lepLooseIsoSF[MuonSystem->nLeptons];
      //     MuonSystem->lepEff[MuonSystem->nLeptons] = MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepLooseIdMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepLooseIsoMCEfficiency[MuonSystem->nLeptons];
      //   }
      //   else{
      //     MuonSystem->lepSF[MuonSystem->nLeptons] = 1.0;
      //     MuonSystem->lepEff[MuonSystem->nLeptons] = 1.0;
      //   }
      // }



      if(!isData)
      {
        for (int i=0; i < nGenParticle; i++)
        { float tmpDR = deltaR(gParticleEta[i],gParticlePhi[i],MuonSystem->lepEta[MuonSystem->nLeptons],MuonSystem->lepPhi[MuonSystem->nLeptons]);
          if ((abs(gParticleId[i]) == 13) && abs(gParticleMotherId[i]) == 23 && tmpDR<0.4) MuonSystem->lepFromZ[MuonSystem->nLeptons] = true;
        }
      }
      MuonSystem->nLeptons++;
    }

    if(MuonSystem->nLeptons<2)continue;
    // if(MuonSystem->category!=2)continue;
    // if (abs(MuonSystem->lepPdgId[0])!=13)continue;
    if (abs(MuonSystem->ZMass)<50)continue;
    // if (abs(MuonSystem->lepPt[0])<50)continue;
    // if (abs(MuonSystem->lepPt[1])<50)continue;

    if(MuonSystem->ZMass<120)continue;
    // if(MuonSystem->lepTag[0] == MuonSystem->lepTag[1]) continue; //require one tag one probe

    if(!isData)
    {
      MuonSystem->lepOverallSF = 1.0 - (1.0 - MuonSystem->lepSF[0] * MuonSystem->lepEff[0]) * (1 - MuonSystem->lepSF[1] * MuonSystem->lepEff[1]);
      MuonSystem->lepOverallSF = MuonSystem->lepOverallSF / (1.0 - (1.0 - MuonSystem->lepEff[0]) * (1 - MuonSystem->lepEff[1]));
    }
    else{
      MuonSystem->lepOverallSF = 1.0;
    }


    TLorentzVector met;
    met.SetPtEtaPhiE(metType1Pt,0,metType1Phi,metType1Pt);
    if ( Leptons.size() > 0 )
    {
      TLorentzVector visible = Leptons[0].lepton;
      MuonSystem->MT = GetMT(visible,met);
    }

    //-----------------------------------------------
    //Select Jets
    //-----------------------------------------------
    std::vector<jets> Jets;
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

      //------------------------------------------------------------
      //Apply Jet Energy and Resolution Corrections
      //------------------------------------------------------------
      // double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
         // fixedGridRhoFastjetAll, jetJetArea[i] , JetCorrector);
     // cout<<"before JEC"<<endl;
     double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
                 fixedGridRhoFastjetAll, jetJetArea[i],
                 runNum,
                 JetCorrectorIOV,JetCorrector);
      JEC = 1.0;
      double jetCorrPt = jetPt[i]*JEC;
      double jetCorrE = jetE[i]*JEC;
      // cout<<JEC<<endl;
      // cout<<"corrected pt: "<<jetCorrPt<<", "<<"eta: "<<jetEta[i]<<endl;
      TLorentzVector thisJet = makeTLorentzVector( jetCorrPt, jetEta[i], jetPhi[i], jetCorrE );

      // if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
      if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
      // if ( !jetPassIDLoose[i] ) continue;

      jets tmpJet;
      tmpJet.jet    = thisJet;
      Jets.push_back(tmpJet);

      }

      //-----------------------------
      //Require at least 2 jets
      //-----------------------------
      // if( Jets.size() < 2 ) continue;
      // if (triggered) trig_lepId_dijet->Fill(1);
      sort(Jets.begin(), Jets.end(), my_largest_pt_jet);

      for ( auto &tmp : Jets )
      {
        MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
        MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
        MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
        MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
        MuonSystem->nJets++;
      }



    //-----------------------------
    // CSC INFO
    //-----------------------------
    //
    // MuonSystem->nRpc  = 0;
    // for (int i = 0; i < nRpc; i++) {
    //   MuonSystem->rpcPhi[MuonSystem->nRpc]           = rpcPhi[i];   //[ndt]
    //   MuonSystem->rpcEta[MuonSystem->nRpc]           = rpcEta[i];   //[ndt]
    //   float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
    //   if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)MuonSystem->rpc_RE12[MuonSystem->nRpc] = true;
    //   if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->rpc_RB1[MuonSystem->nRpc] = true;
    //
    //   MuonSystem->nRpc++;
    // }
    //
    //
    // MuonSystem->nDtSeg  = 0;
    // for (int i = 0; i < nDtRechits; i++) {
    //   MuonSystem->dtSegPhi[MuonSystem->nDtSeg]           = dtSegPhi[i];   //[ndt]
    //   MuonSystem->dtSegEta[MuonSystem->nDtSeg]           = dtSegEta[i];   //[ndt]
    //   MuonSystem->dtSegStation[MuonSystem->nDtSeg] = dtSegStation[i];
    //   MuonSystem->dtSegWheel[MuonSystem->nDtSeg] = dtSegWheel[i];
    //   MuonSystem->nDtSeg++;
    // }

    // count dt rechits
    MuonSystem->nDTRechits  = 0;
    for (int i = 0; i < nDtRechits; i++) {
      MuonSystem->dtRechitsPhi[MuonSystem->nDTRechits]           = dtRechitPhi[i];   //[ndt]
      MuonSystem->dtRechitsEta[MuonSystem->nDTRechits]           = dtRechitEta[i];   //[ndt]
      MuonSystem->dtRechitsStation[MuonSystem->nDTRechits] = dtRechitStation[i];
      MuonSystem->dtRechitsWheel[MuonSystem->nDTRechits] = dtRechitWheel[i];
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

      double half_sector = TMath::Pi()/12.0; //12 sector of DT in 360 degree

      if (dtRechitPhi[i] < 1*half_sector && dtRechitPhi[i] >= 1*half_sector )MuonSystem->nDTRechitsSector[dtRechitStation[i]-1][dtRechitWheel[i]+2][0]++;
      if (dtRechitPhi[i] < 3*half_sector && dtRechitPhi[i] >= 1*half_sector )MuonSystem->nDTRechitsSector[dtRechitStation[i]-1][dtRechitWheel[i]+2][1]++;
      if (dtRechitPhi[i] < 5*half_sector && dtRechitPhi[i] >= 3*half_sector )MuonSystem->nDTRechitsSector[dtRechitStation[i]-1][dtRechitWheel[i]+2][2]++;
      if (dtRechitPhi[i] < 7*half_sector && dtRechitPhi[i] >= 5*half_sector )MuonSystem->nDTRechitsSector[dtRechitStation[i]-1][dtRechitWheel[i]+2][3]++;
      if (dtRechitPhi[i] < 9*half_sector && dtRechitPhi[i] >= 7*half_sector )MuonSystem->nDTRechitsSector[dtRechitStation[i]-1][dtRechitWheel[i]+2][4]++;
      if (dtRechitPhi[i] < 11*half_sector && dtRechitPhi[i] >= 9*half_sector )MuonSystem->nDTRechitsSector[dtRechitStation[i]-1][dtRechitWheel[i]+2][5]++;
      if (dtRechitPhi[i] < -11*half_sector && dtRechitPhi[i] >= 11*half_sector )MuonSystem->nDTRechitsSector[dtRechitStation[i]-1][dtRechitWheel[i]+2][6]++;
      if (dtRechitPhi[i] < -9*half_sector && dtRechitPhi[i] >= -11*half_sector )MuonSystem->nDTRechitsSector[dtRechitStation[i]-1][dtRechitWheel[i]+2][7]++;
      if (dtRechitPhi[i] < -7*half_sector && dtRechitPhi[i] >= -9*half_sector )MuonSystem->nDTRechitsSector[dtRechitStation[i]-1][dtRechitWheel[i]+2][8]++;
      if (dtRechitPhi[i] < -5*half_sector && dtRechitPhi[i] >= -7*half_sector )MuonSystem->nDTRechitsSector[dtRechitStation[i]-1][dtRechitWheel[i]+2][9]++;
      if (dtRechitPhi[i] < -3*half_sector && dtRechitPhi[i] >= -5*half_sector )MuonSystem->nDTRechitsSector[dtRechitStation[i]-1][dtRechitWheel[i]+2][10]++;
      if (dtRechitPhi[i] < -1*half_sector && dtRechitPhi[i] >= -3*half_sector )MuonSystem->nDTRechitsSector[dtRechitStation[i]-1][dtRechitWheel[i]+2][11]++;



      if (dtRechitStation[i] == 1) MuonSystem->nDTRechitsStation1++;
      if (dtRechitStation[i] == 2) MuonSystem->nDTRechitsStation2++;
      if (dtRechitStation[i] == 3) MuonSystem->nDTRechitsStation3++;
      if (dtRechitStation[i] == 4) MuonSystem->nDTRechitsStation4++;

      if (dtRechitWheel[i] == -2) MuonSystem->nDTRechitsWheelMinus2++;
      if (dtRechitWheel[i] == -1) MuonSystem->nDTRechitsWheelMinus1++;
      if (dtRechitWheel[i] == 0) MuonSystem->nDTRechitsWheel0++;
      if (dtRechitWheel[i] == 1) MuonSystem->nDTRechitsWheelPlus1++;
      if (dtRechitWheel[i] == 2) MuonSystem->nDTRechitsWheelPlus2++;


      MuonSystem->nDTRechits++;
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



    // cout << "Number of rec hits: "<<ncscRechits<<endl;
    vector<Point> points;
    vector<int> cscRechitsClusterId;
    points.clear();
    MuonSystem->nCscRechits  = 0;

    for (int i = 0; i < ncscRechits; i++) {

      //pick out the right bits for chamber
      int chamber = ((cscRechitsDetId[i] >> 3) & 077); //https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonDetId/interface/CSCDetId.h#L147
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
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);
      cscRechitsClusterId.push_back(-1);


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
        // if (abs(tmp.nCscSegments) < min_point)continue;

        MuonSystem->cscRechitClusterX[MuonSystem->nCscRechitClusters] =tmp.x;
        MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] =tmp.y;
        MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] =tmp.z;
        MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] = tmp.t;
        MuonSystem->cscRechitClusterTimeTotal[MuonSystem->nCscRechitClusters] = tmp.tTotal;
        MuonSystem->cscRechitClusterTimeWeighted[MuonSystem->nCscRechitClusters] = tmp.tWeighted;

        // MuonSystem->cscRechitClusterTimeWire[MuonSystem->nCscRechitClusters] = tmp.tWire;
        // MuonSystem->cscRechitClusterTimeWirePruned[MuonSystem->nCscRechitClusters] = tmp.tWirePruned;
        MuonSystem->cscRechitClusterTimeSpreadWeighted[MuonSystem->nCscRechitClusters] = tmp.TSpreadWeighted;
        MuonSystem->cscRechitClusterTimeSpreadWeightedAll[MuonSystem->nCscRechitClusters] = tmp.TSpreadWeightedAll;

        MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] =tmp.eta;
        MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters] = tmp.phi;
        MuonSystem->cscRechitClusterMajorAxis[MuonSystem->nCscRechitClusters] =tmp.MajorAxis;
        MuonSystem->cscRechitClusterMinorAxis[MuonSystem->nCscRechitClusters] =tmp.MinorAxis;
        MuonSystem->cscRechitClusterXSpread[MuonSystem->nCscRechitClusters] =tmp.XSpread;
        MuonSystem->cscRechitClusterYSpread[MuonSystem->nCscRechitClusters] =tmp.YSpread;
        MuonSystem->cscRechitClusterZSpread[MuonSystem->nCscRechitClusters] =tmp.ZSpread;
        MuonSystem->cscRechitClusterEtaPhiSpread[MuonSystem->nCscRechitClusters] =tmp.EtaPhiSpread;
        MuonSystem->cscRechitClusterXYSpread[MuonSystem->nCscRechitClusters] =tmp.XYSpread;
        MuonSystem->cscRechitClusterRSpread[MuonSystem->nCscRechitClusters] =tmp.RSpread;


        MuonSystem->cscRechitClusterEtaSpread[MuonSystem->nCscRechitClusters] =tmp.EtaSpread;
        MuonSystem->cscRechitClusterPhiSpread[MuonSystem->nCscRechitClusters] = tmp.PhiSpread;
        MuonSystem->cscRechitClusterTimeSpread[MuonSystem->nCscRechitClusters] = tmp.TSpread;
        // MuonSystem->cscRechitClusterTimeWireSpread[MuonSystem->nCscRechitClusters] = tmp.TWireSpread;

        // MuonSystem->cscRechitClusterTimeTotalSpread[MuonSystem->nCscRechitClusters] = tmp.TTotalSpread;
        // MuonSystem->cscRechitClusterTimeTotalSpreadPruned[MuonSystem->nCscRechitClusters] = tmp.TTotalSpreadPruned;

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
        MuonSystem->cscRechitClusterNStation[MuonSystem->nCscRechitClusters] = tmp.nStation;
        MuonSystem->cscRechitClusterNStation5[MuonSystem->nCscRechitClusters] = tmp.nStation5;
        MuonSystem->cscRechitClusterNStation10[MuonSystem->nCscRechitClusters] = tmp.nStation10;
        MuonSystem->cscRechitClusterNStation10perc[MuonSystem->nCscRechitClusters] = tmp.nStation10perc;
        MuonSystem->cscRechitClusterAvgStation[MuonSystem->nCscRechitClusters] = tmp.avgStation;
        MuonSystem->cscRechitClusterAvgStation5[MuonSystem->nCscRechitClusters] = tmp.avgStation5;
        MuonSystem->cscRechitClusterAvgStation10[MuonSystem->nCscRechitClusters] = tmp.avgStation10;
        MuonSystem->cscRechitClusterAvgStation10perc[MuonSystem->nCscRechitClusters] = tmp.avgStation10perc;
        MuonSystem->cscRechitClusterMe11Ratio[MuonSystem->nCscRechitClusters] = tmp.Me11Ratio;
        MuonSystem->cscRechitClusterMe12Ratio[MuonSystem->nCscRechitClusters] = tmp.Me12Ratio;
        if(MuonSystem->category == 2)
        {
          if (RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex1], MuonSystem->lepPhi[MuonSystem->ZleptonIndex1], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4) {
            MuonSystem->cscRechitClusterZLep1[MuonSystem->nCscRechitClusters] = true;
            MuonSystem->cscRechitClusterZLep1Id[MuonSystem->nCscRechitClusters] = MuonSystem->lepPdgId[MuonSystem->ZleptonIndex1];
            // MuonSystem->cscRechitClusterZLep1LooseIso[MuonSystem->nCscRechitClusters]  = MuonSystem->lepPassLooseIso[MuonSystem->ZleptonIndex1];
            // MuonSystem->cscRechitClusterZLep1TightIso[MuonSystem->nCscRechitClusters]  = MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
            // MuonSystem->cscRechitClusterZLep1VTightIso[MuonSystem->nCscRechitClusters]  = MuonSystem->lepPassVTightIso[MuonSystem->ZleptonIndex1];
            // MuonSystem->cscRechitClusterZLep1VVTightIso[MuonSystem->nCscRechitClusters]  = MuonSystem->lepPassVVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitClusterZLep1TightId[MuonSystem->nCscRechitClusters]  = MuonSystem->lepTightId[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitClusterZLep1Tag[MuonSystem->nCscRechitClusters] = MuonSystem->lepTag[MuonSystem->ZleptonIndex1];

          }
          if (RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex2], MuonSystem->lepPhi[MuonSystem->ZleptonIndex2], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4) {
            MuonSystem->cscRechitClusterZLep2[MuonSystem->nCscRechitClusters] = true;
            MuonSystem->cscRechitClusterZLep2Id[MuonSystem->nCscRechitClusters] = MuonSystem->lepPdgId[MuonSystem->ZleptonIndex2];
            // MuonSystem->cscRechitClusterZLep2LooseIso[MuonSystem->nCscRechitClusters]  = MuonSystem->lepPassLooseIso[MuonSystem->ZleptonIndex1];
            // MuonSystem->cscRechitClusterZLep2TightIso[MuonSystem->nCscRechitClusters]  = MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
            // MuonSystem->cscRechitClusterZLep2VTightIso[MuonSystem->nCscRechitClusters]  = MuonSystem->lepPassVTightIso[MuonSystem->ZleptonIndex1];
            // MuonSystem->cscRechitClusterZLep2VVTightIso[MuonSystem->nCscRechitClusters]  = MuonSystem->lepPassVVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitClusterZLep2TightId[MuonSystem->nCscRechitClusters]  = MuonSystem->lepTightId[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitClusterZLep2Tag[MuonSystem->nCscRechitClusters] = MuonSystem->lepTag[MuonSystem->ZleptonIndex1];


          }



        }
        //Jet veto/ muon veto
        MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterGenMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterGenMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterIsoMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;


        float min_deltaR = 15.;
        int index = 999;
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
            MuonSystem->cscRechitClusterMuonVetoLooseIso[MuonSystem->nCscRechitClusters]  = muonIso<0.25;
            MuonSystem->cscRechitClusterMuonVetoTightIso[MuonSystem->nCscRechitClusters]  = muonIso<0.15;
            MuonSystem->cscRechitClusterMuonVetoVTightIso[MuonSystem->nCscRechitClusters]  = muonIso<0.10;
            MuonSystem->cscRechitClusterMuonVetoVVTightIso[MuonSystem->nCscRechitClusters]  = muonIso<0.05;
            MuonSystem->cscRechitClusterMuonVetoTightId[MuonSystem->nCscRechitClusters]  = isMuonPOGTightMuon(i);

          }
        }

        // match to gen-level muon
        if(!isData)
        {
          for(int i = 0; i < nGenParticle; i++)
          {
            if (abs(gParticleId[i])!=13) continue;
            if (abs(gParticleMotherId[i])>24 || abs(gParticleMotherId[i])<23) continue;
            if (abs(gParticleStatus[i])!=1)continue;
            // if (fabs(muonEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 ) {
              MuonSystem->cscRechitClusterGenMuonVetoPt[MuonSystem->nCscRechitClusters]  = gParticlePt[i];
              MuonSystem->cscRechitClusterGenMuonVetoE[MuonSystem->nCscRechitClusters]  = gParticleE[i];

            }

          }
        }

        // cout<<MuonSystem->nCscRechitClusters<<endl;
        MuonSystem->nCscRechitClusters++;
    }


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
        p.clusterID = UNCLASSIFIED;
        points.push_back(p);
      }
    //Do DBSCAN Clustering for DT
    // int min_point_dt = 30;  //minimum number of segments to call it a cluster
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

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator (seed);

        // default_random_engine generator;
        uniform_real_distribution<double> distribution(0.0,1.0);

       for (int i=0; i<12; ++i) {
         if ( distribution(generator) < 0.05) MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters]++;
       }
       for (int i=0; i<12; ++i) {
         if ( distribution(generator) < 0.05) MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters]++;
       }
       for (int i=0; i<12; ++i) {
         if ( distribution(generator) < 0.05) MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters]++;
       }
       for (int i=0; i<8; ++i) {
         if ( distribution(generator) < 0.05) MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters]++;
       }

       MuonSystem->dtRechitClusterNoiseHit[MuonSystem->nDtRechitClusters] = MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters] +
                                                                            MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters] +
                                                                            MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters] +
                                                                            MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters];


        MuonSystem->dtRechitClusterSize[MuonSystem->nDtRechitClusters] = tmp.nCscSegments;


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

        if(MuonSystem->category == 2)
        {
          if (RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex1], MuonSystem->lepPhi[MuonSystem->ZleptonIndex1], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
            MuonSystem->dtRechitClusterZLep1[MuonSystem->nDtRechitClusters] = true;
            MuonSystem->dtRechitClusterZLep1Id[MuonSystem->nDtRechitClusters] = MuonSystem->lepPdgId[MuonSystem->ZleptonIndex1];
            MuonSystem->dtRechitClusterZLep1LooseIso[MuonSystem->nDtRechitClusters]  = MuonSystem->lepPassLooseIso[MuonSystem->ZleptonIndex1];
            MuonSystem->dtRechitClusterZLep1TightIso[MuonSystem->nDtRechitClusters]  = MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->dtRechitClusterZLep1VTightIso[MuonSystem->nDtRechitClusters]  = MuonSystem->lepPassVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->dtRechitClusterZLep1VVTightIso[MuonSystem->nDtRechitClusters]  = MuonSystem->lepPassVVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->dtRechitClusterZLep1TightId[MuonSystem->nDtRechitClusters]  = MuonSystem->lepTightId[MuonSystem->ZleptonIndex1];
            MuonSystem->dtRechitClusterZLep1Tag[MuonSystem->nDtRechitClusters] = MuonSystem->lepTag[MuonSystem->ZleptonIndex1];

          }
          if (RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex2], MuonSystem->lepPhi[MuonSystem->ZleptonIndex2], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
            MuonSystem->dtRechitClusterZLep2[MuonSystem->nDtRechitClusters] = true;
            MuonSystem->dtRechitClusterZLep2Id[MuonSystem->nDtRechitClusters] = MuonSystem->lepPdgId[MuonSystem->ZleptonIndex2];
            MuonSystem->dtRechitClusterZLep2LooseIso[MuonSystem->nDtRechitClusters]  = MuonSystem->lepPassLooseIso[MuonSystem->ZleptonIndex1];
            MuonSystem->dtRechitClusterZLep2TightIso[MuonSystem->nDtRechitClusters]  = MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->dtRechitClusterZLep2VTightIso[MuonSystem->nDtRechitClusters]  = MuonSystem->lepPassVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->dtRechitClusterZLep2VVTightIso[MuonSystem->nDtRechitClusters]  = MuonSystem->lepPassVVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->dtRechitClusterZLep2TightId[MuonSystem->nDtRechitClusters]  = MuonSystem->lepTightId[MuonSystem->ZleptonIndex1];
            MuonSystem->dtRechitClusterZLep2Tag[MuonSystem->nDtRechitClusters] = MuonSystem->lepTag[MuonSystem->ZleptonIndex1];


          }



        }
        //Jet veto/ muon veto


        float min_deltaR = 15.;
        int index = 999;


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
            MuonSystem->dtRechitClusterMuonVetoLooseId[MuonSystem->nDtRechitClusters]  = isMuonPOGLooseMuon(i);



          }
        }
        // match to gen-level muon
        if(!isData)
        {
          for(int i = 0; i < nGenParticle; i++)
          {
            if (abs(gParticleId[i])!=13) continue;
            if (abs(gParticleMotherId[i])>24 || abs(gParticleMotherId[i])<23) continue;
            if (abs(gParticleStatus[i])!=1)continue;
            // if (fabs(muonEta[i]>3.0)) continue;
            if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 ) {
              MuonSystem->dtRechitClusterGenMuonVetoPt[MuonSystem->nDtRechitClusters]  = gParticlePt[i];
              MuonSystem->dtRechitClusterGenMuonVetoE[MuonSystem->nDtRechitClusters]  = gParticleE[i];

            }

          }
        }


        //match to RPC hits with dPhi<0.5 and same wheel in DT
        for (int i = 0; i < nRpc; i++) {
          float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
          if (rpcRegion[i]!=0) continue;
          if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5 )
          {
            if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
            {
              MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]++;

            }
          }
        }

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
         MuonSystem->dtRechitClusterMetEENoise_dPhi[MuonSystem->nDtRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhiEENoise);



        MuonSystem->nDtRechitClusters++;
      }




    if(!isData && signalScan)
    {
      pair<int,int> smsPair = make_pair(MuonSystem->mX, MuonSystem->ctau);
      Trees2D[smsPair]->Fill();
    }
    else
    {
      MuonSystem->tree_->Fill();
    }



  }
  if(!isData && signalScan)
  {
    for(auto &filePtr : Files2D)
     {
       cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
       filePtr.second->cd();
       Trees2D[filePtr.first]->Write();
       NEvents2D[filePtr.first]->Write("NEvents");
       filePtr.second->Close();

     }
  }
  else
  {
    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();
    NEvents->Write();
    // outFile->Write();
    outFile->Close();
  }

}
