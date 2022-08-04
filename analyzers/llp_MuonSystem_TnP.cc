#include "llp_MuonSystem_TnP.h"
#include "RazorHelper.h"
#include "LiteTreeMuonSystem.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
#include "DBSCAN.h"
//C++ includes
#include "assert.h"

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
  bool passId;
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



void llp_MuonSystem_TnP::Analyze(bool isData, int options, string outputfilename, string analysisTag)
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


  LiteTreeMuonSystem *MuonSystem = new LiteTreeMuonSystem;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  // for signals, need one output file for each signal point
  map<pair<int,int>, TFile*> Files2D;
  map<pair<int,int>, TTree*> Trees2D;
  map<pair<int,int>, TH1F*> NEvents2D;
  // map<pair<int,int>, TH1F*> smsSumWeights2D;
  // map<pair<int,int>, TH1F*> smsSumScaleWeights2D;
  // map<pair<int,int>, TH1F*> smsSumPdfWeights2D;
  // map<pair<int,int>, TH1F*> smsNISRJets2D;
  // map<pair<int,int>, TH1F*> smsPtISR2D;
  // map<pair<int,int>, TH1F*> smsNPV2D;




  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
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
    //GetEntry(ientry);
    nb = fChain->GetEntry(jentry); nbytes += nb;

    //fill normalization histogram
    //std::cout << "deb0 " << jentry << std::endl;
    MuonSystem->InitVariables();
    //std::cout << "deb1 " << jentry << std::endl;


    // cout<<*lheComments<<endl;
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

      // if (mh2 != mh || mx2!=mx || ctau2!=ctau) continue;
      // cout<<*lheComments<<endl;

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
      MuonSystem->higgsPtWeight = helper->getHiggsPtWeight(MuonSystem->gHiggsPt);
      for (int i = 0; i < 9; i++)
      {
        MuonSystem->higgsPtWeightSys[i] = helper->getHiggsPtWeightSys(MuonSystem->gHiggsPt, i) / MuonSystem->higgsPtWeight;

      }
      MuonSystem->sf_facScaleUp = MuonSystem->higgsPtWeightSys[5];
      MuonSystem->sf_facScaleDown = MuonSystem->higgsPtWeightSys[3];
      MuonSystem->sf_renScaleUp = MuonSystem->higgsPtWeightSys[7];
      MuonSystem->sf_renScaleDown = MuonSystem->higgsPtWeightSys[1];
      MuonSystem->sf_facRenScaleUp = MuonSystem->higgsPtWeightSys[8];
      MuonSystem->sf_facRenScaleDown = MuonSystem->higgsPtWeightSys[0];

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

    // if (!MuonSystem->Flag2_all) continue;
    //*************************************************************************
    //Start Object Selection
    //*************************************************************************

    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
    MuonSystem->nMuons = 0;
    for( int i = 0; i < nMuons; i++ )
    {
      if(muonPt[i] < 20) continue;
      MuonSystem->muonPt[MuonSystem->nMuons] = muonPt[i];
      MuonSystem->muonEta[MuonSystem->nMuons] = muonEta[i];
      MuonSystem->muonPhi[MuonSystem->nMuons] = muonPhi[i];
      MuonSystem->nMuons++;
      float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];


      if(!isMuonPOGLooseMuon(i)) continue;
      if(muonPt[i] < 50) continue;
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
      tmpMuon.passId = isMuonPOGTightMuon(i, true, false);
      tmpMuon.passLooseIso = muonIso<0.25;
      tmpMuon.passTightIso = muonIso<0.15;
      tmpMuon.passVTightIso = muonIso<0.10;
      tmpMuon.passVVTightIso = muonIso<0.05;


      // tmpMuon.passId = isLooseMuon(i);
      tmpMuon.passVetoId = false;
      // tmpMuon.passVetoPOGId = false;

      Leptons.push_back(tmpMuon);
    }
    // cout<<"jere"<<endl;

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
      tmpElectron.passId = isEGammaPOGTightElectron(i, true, true, true, "Summer16");
      tmpElectron.passVetoId = isEGammaPOGVetoElectron(i, true, true, true, "Summer16");
      Leptons.push_back(tmpElectron);
    }

    sort(Leptons.begin(), Leptons.end(), my_largest_pt);
    //std::cout << "deb7 " << jentry << std::endl;


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
      for( uint j = i+1; j < Leptons.size(); j++ )
      {
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
    if (foundZ  && Leptons.size() == 2 )
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
      for ( unsigned int i = Leptons.size(); i>0; --i )
      {
        int index = i-1;
        if (abs(Leptons[index].pdgId) == 13 && Leptons[index].lepton.Pt() < wh_muonPt_cut)
        {
           Leptons.erase(Leptons.begin() + index);
        }
        else if (abs(Leptons[index].pdgId) == 11 && Leptons[index].lepton.Pt() < wh_elePt_cut)
        {
          Leptons.erase(Leptons.begin() + index);
        }
      }
      if (Leptons.size() == 1) MuonSystem->category = 1;
      else MuonSystem->category = 0;
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
      MuonSystem->lepPassId[MuonSystem->nLeptons] = tmp.passId;
      MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
      MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
      MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
      MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
      if (!isData)
      {
        MuonSystem->lepTriggerSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), true, true, true);
        MuonSystem->lepTightIdSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, false);
        MuonSystem->lepLooseIdSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, true);
        MuonSystem->lepTightIsoSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, false);
        MuonSystem->lepLooseIsoSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, true);
        MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), true, true, true);
        MuonSystem->lepTightIdMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, false);
        MuonSystem->lepLooseIdMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, true);
        MuonSystem->lepTightIsoMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, false);
        MuonSystem->lepLooseIsoMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, true);
      }



      if (MuonSystem->lepPassTightIso[MuonSystem->nLeptons] && MuonSystem->lepPassId[MuonSystem->nLeptons] && !tag)
      {
        MuonSystem->lepTag[MuonSystem->nLeptons] = true;
        if(!isData)
        {
          MuonSystem->lepSF[MuonSystem->nLeptons] = MuonSystem->lepTriggerSF[MuonSystem->nLeptons] * MuonSystem->lepTightIdSF[MuonSystem->nLeptons] * MuonSystem->lepTightIsoSF[MuonSystem->nLeptons];
          MuonSystem->lepEff[MuonSystem->nLeptons] = MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepTightIdMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepTightIsoMCEfficiency[MuonSystem->nLeptons];
        }
        else{
          MuonSystem->lepSF[MuonSystem->nLeptons] = 1.0;
          MuonSystem->lepEff[MuonSystem->nLeptons] = 1.0;
        }
        tag = true;
      }
      else
      {
        MuonSystem->lepTag[MuonSystem->nLeptons] = false;
        if(!isData)
        {
          MuonSystem->lepSF[MuonSystem->nLeptons] = MuonSystem->lepTriggerSF[MuonSystem->nLeptons] * MuonSystem->lepLooseIdSF[MuonSystem->nLeptons] * MuonSystem->lepLooseIsoSF[MuonSystem->nLeptons];
          MuonSystem->lepEff[MuonSystem->nLeptons] = MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepLooseIdMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepLooseIsoMCEfficiency[MuonSystem->nLeptons];
        }
        else{
          MuonSystem->lepSF[MuonSystem->nLeptons] = 1.0;
          MuonSystem->lepEff[MuonSystem->nLeptons] = 1.0;
        }
      }



      if(!isData)
      {
        for (int i=0; i < nGenParticle; i++)
        { float tmpDR = deltaR(gParticleEta[i],gParticlePhi[i],MuonSystem->lepEta[MuonSystem->nLeptons],MuonSystem->lepPhi[MuonSystem->nLeptons]);
          if ((abs(gParticleId[i]) == 13) && abs(gParticleMotherId[i]) == 23 && tmpDR<0.4) MuonSystem->lepFromZ[MuonSystem->nLeptons] = true;
        }
      }
      MuonSystem->nLeptons++;
    }

    if(MuonSystem->nLeptons!=2)continue;
    if(MuonSystem->category!=2)continue;
    if (abs(MuonSystem->lepPdgId[0])!=13)continue;
    if (abs(MuonSystem->ZMass)<50)continue;
    if (abs(MuonSystem->lepPt[0])<50)continue;
    if (abs(MuonSystem->lepPt[1])<50)continue;
    if(MuonSystem->ZMass<120)continue;
    {
      if(MuonSystem->lepTag[0] == MuonSystem->lepTag[1]) continue; //require one tag one probe
      if(!isData)
      {
        MuonSystem->lepOverallSF = 1.0 - (1.0 - MuonSystem->lepSF[0] * MuonSystem->lepEff[0]) * (1 - MuonSystem->lepSF[1] * MuonSystem->lepEff[1]);
        MuonSystem->lepOverallSF = MuonSystem->lepOverallSF / (1.0 - (1.0 - MuonSystem->lepEff[0]) * (1 - MuonSystem->lepEff[1]));
      }
      else{
        MuonSystem->lepOverallSF = 1.0;
      }


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

    MuonSystem->nRpc  = 0;
    for (int i = 0; i < nRpc; i++) {
      MuonSystem->rpcPhi[MuonSystem->nRpc]           = rpcPhi[i];   //[ndt]
      MuonSystem->rpcEta[MuonSystem->nRpc]           = rpcEta[i];   //[ndt]
      float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
      if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730)MuonSystem->rpc_RE12[MuonSystem->nRpc] = true;
      if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->rpc_RB1[MuonSystem->nRpc] = true;

      MuonSystem->nRpc++;
    }


    MuonSystem->nDtSeg  = 0;
    for (int i = 0; i < nDtRechits; i++) {
      MuonSystem->dtSegPhi[MuonSystem->nDtSeg]           = dtSegPhi[i];   //[ndt]
      MuonSystem->dtSegEta[MuonSystem->nDtSeg]           = dtSegEta[i];   //[ndt]
      MuonSystem->dtSegStation[MuonSystem->nDtSeg] = dtSegStation[i];
      MuonSystem->dtSegWheel[MuonSystem->nDtSeg] = dtSegWheel[i];
      MuonSystem->nDtSeg++;
    }

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







    // cout << "Number of rec hits: "<<ncscRechits<<endl;
    vector<Point> points;
    vector<int> cscRechitsClusterId;
    points.clear();
    MuonSystem->nCscRechits  = 0;

    for (int i = 0; i < ncscRechits; i++) {
      MuonSystem->cscRechitsPhi[MuonSystem->nCscRechits]           = cscRechitsPhi[i];   //[nCsc]
      MuonSystem->cscRechitsEta[MuonSystem->nCscRechits]           = cscRechitsEta[i];   //[nCsc]
      MuonSystem->cscRechitsTpeak[MuonSystem->nCscRechits] = cscRechitsTpeak[i];
      MuonSystem->cscRechitsTwire[MuonSystem->nCscRechits] = cscRechitsTwire[i];
      MuonSystem->cscRechitsStation[MuonSystem->nCscRechits] = cscRechitsStation[i];
      MuonSystem->cscRechitsChamber[MuonSystem->nCscRechits] = cscRechitsChamber[i];
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
      MuonSystem->nCscRechits++;
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
    int min_point = 130;  //minimum number of segments to call it a cluster
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




    MuonSystem->nCscRechitClusters3 = 0;
    for ( auto &tmp : ds.CscCluster ) {
        // if (abs(tmp.tTotal) > 12.5)continue;
        if(abs(tmp.eta)>=2.0) continue;
        // if(tmp.TSpread > 20)continue;
        if (abs(tmp.maxChamber)<=12)continue;
        if (abs(tmp.nCscSegments) < min_point)continue;

        MuonSystem->cscRechitCluster3X[MuonSystem->nCscRechitClusters3] =tmp.x;
        MuonSystem->cscRechitCluster3Y[MuonSystem->nCscRechitClusters3] =tmp.y;
        MuonSystem->cscRechitCluster3Z[MuonSystem->nCscRechitClusters3] =tmp.z;
        MuonSystem->cscRechitCluster3Time[MuonSystem->nCscRechitClusters3] = tmp.t;
        MuonSystem->cscRechitCluster3TimeTotal[MuonSystem->nCscRechitClusters3] = tmp.tTotal;
        // MuonSystem->cscRechitCluster3TimeWire[MuonSystem->nCscRechitClusters3] = tmp.tWire;
        // MuonSystem->cscRechitCluster3TimeWirePruned[MuonSystem->nCscRechitClusters3] = tmp.tWirePruned;

        MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3] =tmp.eta;
        MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3] = tmp.phi;
        MuonSystem->cscRechitCluster3MajorAxis[MuonSystem->nCscRechitClusters3] =tmp.MajorAxis;
        MuonSystem->cscRechitCluster3MinorAxis[MuonSystem->nCscRechitClusters3] =tmp.MinorAxis;
        MuonSystem->cscRechitCluster3XSpread[MuonSystem->nCscRechitClusters3] =tmp.XSpread;
        MuonSystem->cscRechitCluster3YSpread[MuonSystem->nCscRechitClusters3] =tmp.YSpread;
        MuonSystem->cscRechitCluster3ZSpread[MuonSystem->nCscRechitClusters3] =tmp.ZSpread;
        MuonSystem->cscRechitCluster3EtaPhiSpread[MuonSystem->nCscRechitClusters3] =tmp.EtaPhiSpread;
        MuonSystem->cscRechitCluster3XYSpread[MuonSystem->nCscRechitClusters3] =tmp.XYSpread;
        MuonSystem->cscRechitCluster3RSpread[MuonSystem->nCscRechitClusters3] =tmp.RSpread;


        MuonSystem->cscRechitCluster3EtaSpread[MuonSystem->nCscRechitClusters3] =tmp.EtaSpread;
        MuonSystem->cscRechitCluster3PhiSpread[MuonSystem->nCscRechitClusters3] = tmp.PhiSpread;
        MuonSystem->cscRechitCluster3TimeSpread[MuonSystem->nCscRechitClusters3] = tmp.TSpread;
        // MuonSystem->cscRechitCluster3TimeWireSpread[MuonSystem->nCscRechitClusters3] = tmp.TWireSpread;

        // MuonSystem->cscRechitCluster3TimeTotalSpread[MuonSystem->nCscRechitClusters3] = tmp.TTotalSpread;
        // MuonSystem->cscRechitCluster3TimeTotalSpreadPruned[MuonSystem->nCscRechitClusters3] = tmp.TTotalSpreadPruned;

        MuonSystem->cscRechitCluster3Size[MuonSystem->nCscRechitClusters3] = tmp.nCscSegments;

        MuonSystem->cscRechitCluster3NRechitChamberPlus11[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus11;
        MuonSystem->cscRechitCluster3NRechitChamberPlus12[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus12;
        MuonSystem->cscRechitCluster3NRechitChamberPlus13[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus13;
        MuonSystem->cscRechitCluster3NRechitChamberPlus21[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus21;
        MuonSystem->cscRechitCluster3NRechitChamberPlus22[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus22;
        MuonSystem->cscRechitCluster3NRechitChamberPlus31[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus31;
        MuonSystem->cscRechitCluster3NRechitChamberPlus32[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus32;
        MuonSystem->cscRechitCluster3NRechitChamberPlus41[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus41;
        MuonSystem->cscRechitCluster3NRechitChamberPlus42[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberPlus42;
        MuonSystem->cscRechitCluster3NRechitChamberMinus11[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus11;
        MuonSystem->cscRechitCluster3NRechitChamberMinus12[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus12;
        MuonSystem->cscRechitCluster3NRechitChamberMinus13[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus13;
        MuonSystem->cscRechitCluster3NRechitChamberMinus21[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus21;
        MuonSystem->cscRechitCluster3NRechitChamberMinus22[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus22;
        MuonSystem->cscRechitCluster3NRechitChamberMinus31[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus31;
        MuonSystem->cscRechitCluster3NRechitChamberMinus32[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus32;
        MuonSystem->cscRechitCluster3NRechitChamberMinus41[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus41;
        MuonSystem->cscRechitCluster3NRechitChamberMinus42[MuonSystem->nCscRechitClusters3] = tmp.nCscSegmentChamberMinus42;
        MuonSystem->cscRechitCluster3MaxChamber[MuonSystem->nCscRechitClusters3] = tmp.maxChamber;
        MuonSystem->cscRechitCluster3MaxChamberRatio[MuonSystem->nCscRechitClusters3] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
        MuonSystem->cscRechitCluster3NChamber[MuonSystem->nCscRechitClusters3] = tmp.nChamber;
        MuonSystem->cscRechitCluster3MaxStation[MuonSystem->nCscRechitClusters3] = tmp.maxStation;
        MuonSystem->cscRechitCluster3MaxStationRatio[MuonSystem->nCscRechitClusters3] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
        MuonSystem->cscRechitCluster3NStation[MuonSystem->nCscRechitClusters3] = tmp.nStation;
        MuonSystem->cscRechitCluster3NStation5[MuonSystem->nCscRechitClusters3] = tmp.nStation5;
        MuonSystem->cscRechitCluster3NStation10[MuonSystem->nCscRechitClusters3] = tmp.nStation10;
        MuonSystem->cscRechitCluster3NStation10perc[MuonSystem->nCscRechitClusters3] = tmp.nStation10perc;
        MuonSystem->cscRechitCluster3AvgStation[MuonSystem->nCscRechitClusters3] = tmp.avgStation;
        MuonSystem->cscRechitCluster3AvgStation5[MuonSystem->nCscRechitClusters3] = tmp.avgStation5;
        MuonSystem->cscRechitCluster3AvgStation10[MuonSystem->nCscRechitClusters3] = tmp.avgStation10;
        MuonSystem->cscRechitCluster3AvgStation10perc[MuonSystem->nCscRechitClusters3] = tmp.avgStation10perc;
        MuonSystem->cscRechitCluster3Me11Ratio[MuonSystem->nCscRechitClusters3] = tmp.Me11Ratio;
        MuonSystem->cscRechitCluster3Me12Ratio[MuonSystem->nCscRechitClusters3] = tmp.Me12Ratio;
        if(MuonSystem->category == 2)
        {
          if (RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex1], MuonSystem->lepPhi[MuonSystem->ZleptonIndex1], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4) {
            MuonSystem->cscRechitCluster3ZLep1[MuonSystem->nCscRechitClusters3] = true;
            MuonSystem->cscRechitCluster3ZLep1Id[MuonSystem->nCscRechitClusters3] = MuonSystem->lepPdgId[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep1LooseIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassLooseIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep1TightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep1VTightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep1VVTightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassVVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep1TightId[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassId[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep1Tag[MuonSystem->nCscRechitClusters3] = MuonSystem->lepTag[MuonSystem->ZleptonIndex1];

          }
          if (RazorAnalyzer::deltaR(MuonSystem->lepEta[MuonSystem->ZleptonIndex2], MuonSystem->lepPhi[MuonSystem->ZleptonIndex2], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4) {
            MuonSystem->cscRechitCluster3ZLep2[MuonSystem->nCscRechitClusters3] = true;
            MuonSystem->cscRechitCluster3ZLep2Id[MuonSystem->nCscRechitClusters3] = MuonSystem->lepPdgId[MuonSystem->ZleptonIndex2];
            MuonSystem->cscRechitCluster3ZLep2LooseIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassLooseIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep2TightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep2VTightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep2VVTightIso[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassVVTightIso[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep2TightId[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassId[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep2TightId[MuonSystem->nCscRechitClusters3]  = MuonSystem->lepPassId[MuonSystem->ZleptonIndex1];
            MuonSystem->cscRechitCluster3ZLep2Tag[MuonSystem->nCscRechitClusters3] = MuonSystem->lepTag[MuonSystem->ZleptonIndex1];


          }



        }
        //Jet veto/ muon veto
        MuonSystem->cscRechitCluster3JetVetoPt[MuonSystem->nCscRechitClusters3] = 0.0;
        MuonSystem->cscRechitCluster3JetVetoE[MuonSystem->nCscRechitClusters3] = 0.0;
        MuonSystem->cscRechitCluster3MuonVetoPt[MuonSystem->nCscRechitClusters3] = 0.0;
        MuonSystem->cscRechitCluster3MuonVetoE[MuonSystem->nCscRechitClusters3] = 0.0;
        MuonSystem->cscRechitCluster3GenMuonVetoPt[MuonSystem->nCscRechitClusters3] = 0.0;
        MuonSystem->cscRechitCluster3GenMuonVetoE[MuonSystem->nCscRechitClusters3] = 0.0;
        MuonSystem->cscRechitCluster3IsoMuonVetoPt[MuonSystem->nCscRechitClusters3] = 0.0;


        float min_deltaR = 15.;
        int index = 999;

        for(int i = 0; i < nMuons; i++)
        {
          if (fabs(muonEta[i]>3.0)) continue;
          float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / muonPt[i];
          if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 && muonPt[i] > MuonSystem->cscRechitCluster3MuonVetoPt[MuonSystem->nCscRechitClusters3] ) {
            MuonSystem->cscRechitCluster3MuonVetoPt[MuonSystem->nCscRechitClusters3]  = muonPt[i];
            MuonSystem->cscRechitCluster3MuonVetoE[MuonSystem->nCscRechitClusters3]  = muonE[i];
            MuonSystem->cscRechitCluster3MuonVetoPhi[MuonSystem->nCscRechitClusters3]  = muonPhi[i];
            MuonSystem->cscRechitCluster3MuonVetoEta[MuonSystem->nCscRechitClusters3]  = muonEta[i];
            MuonSystem->cscRechitCluster3MuonVetoLooseIso[MuonSystem->nCscRechitClusters3]  = muonIso<0.25;
            MuonSystem->cscRechitCluster3MuonVetoTightIso[MuonSystem->nCscRechitClusters3]  = muonIso<0.15;
            MuonSystem->cscRechitCluster3MuonVetoVTightIso[MuonSystem->nCscRechitClusters3]  = muonIso<0.10;
            MuonSystem->cscRechitCluster3MuonVetoVVTightIso[MuonSystem->nCscRechitClusters3]  = muonIso<0.05;
            MuonSystem->cscRechitCluster3MuonVetoTightId[MuonSystem->nCscRechitClusters3]  = isMuonPOGTightMuon(i);

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
            if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem->cscRechitCluster3Eta[MuonSystem->nCscRechitClusters3],MuonSystem->cscRechitCluster3Phi[MuonSystem->nCscRechitClusters3]) < 0.4 ) {
              MuonSystem->cscRechitCluster3GenMuonVetoPt[MuonSystem->nCscRechitClusters3]  = gParticlePt[i];
              MuonSystem->cscRechitCluster3GenMuonVetoE[MuonSystem->nCscRechitClusters3]  = gParticleE[i];

            }

          }
        }

        // cout<<MuonSystem->nCscRechitClusters3<<endl;
        MuonSystem->nCscRechitClusters3++;
    }


    //  for (unsigned int j = 0; j < tmp.segment_id.size(); j++)
    // {
    //   MuonSystem->cscRechitsCluster2Id[tmp.segment_id[j]] = MuonSystem->nCscRechitClusters2;
    // }
    // if(MuonSystem->nCscRechitClusters3==0 && MuonSystem->nCscRechitClusters2==0 && MuonSystem->nCscRechitClusters==0) continue;
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
