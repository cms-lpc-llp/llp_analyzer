#include "llp_MuonSystem_bdt.h"
#include "RazorHelper.h"
#include "LiteTreeMuonSystem.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"
#include "EnergyScaleCorrection_class.hh"
// #include "DBSCAN.h"
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

int cscChamber(double x, double y, double z)
{
  double r = sqrt(x*x+y*y);
  int sign_z = TMath::Sign(1.0, z);
  // if (r > 80 && r < 283 && abs(z) > 568  && abs(z) < 632) return sign_z*11;
  // if (r > 255 && r < 465.0 && abs(z) > 668.3 && abs(z) < 724) return sign_z*12;
  // if (r > 485.5 && r < 695.5 && abs(z) > 686 && abs(z) < 724) return sign_z*13;
  // if (r > 118.5 && r < 345 && abs(z) > 791 && abs(z) < 849.5) return sign_z*21;
  // if (r > 337.5 && r < 695.5 && abs(z) > 791 && abs(z) < 849.5) return sign_z*22;
  // if (r > 140.5 && r < 345 && abs(z) > 911.5 && abs(z) < 970) return sign_z*31;
  // if (r > 337.5 && r < 695.5 && abs(z) > 911.5 && abs(z) < 970) return sign_z*32;
  // if (r > 157.5 && r < 345 && abs(z) > 1002 && abs(z) < 1060.5) return sign_z*41;
  // if (r > 337.5 && r < 695.5 && abs(z) > 1002 && abs(z) < 1060.5) return sign_z*42;
  if (r < 283 && abs(z) > 568  && abs(z) < 632) return sign_z*11;
  if (r < 470.0 && abs(z) > 668.3 && abs(z) < 724) return sign_z*12;
  if (r > 480.0 && abs(z) > 686 && abs(z) < 724) return sign_z*13;
  if (r < 345 && abs(z) > 791 && abs(z) < 849.5) return sign_z*21;
  if (r > 337.5 && abs(z) > 791 && abs(z) < 849.5) return sign_z*22;
  if (r < 345 && abs(z) > 911.5 && abs(z) < 970) return sign_z*31;
  if (r > 337.5 && abs(z) > 911.5 && abs(z) < 970) return sign_z*32;
  if (r < 345 && abs(z) > 1002 && abs(z) < 1060.5) return sign_z*41;
  if (r > 337.5 && abs(z) > 1002 && abs(z) < 1060.5) return sign_z*42;
  return -999;
};
int cscStation(double x, double y, double z)
{
  double r = sqrt(x*x+y*y);
  // z = abs(z);
  int sign_z = TMath::Sign(1.0, z);
  if (r < 283 && abs(z) > 568  && abs(z) < 632) return sign_z*1;
  if (r < 470.0 && abs(z) > 668.3 && abs(z) < 724) return sign_z*1;
  if (r > 480.0 && abs(z) > 686 && abs(z) < 724) return sign_z*1;
  if (r < 345 && abs(z) > 791 && abs(z) < 849.5) return sign_z*2;
  if (r > 337.5 && abs(z) > 791 && abs(z) < 849.5) return sign_z*2;
  if (r < 345 && abs(z) > 911.5 && abs(z) < 970) return sign_z*3;
  if (r > 337.5 && abs(z) > 911.5 && abs(z) < 970) return sign_z*3;
  if (r < 345 && abs(z) > 1002 && abs(z) < 1060.5) return sign_z*4;
  if (r > 337.5 && abs(z) > 1002 && abs(z) < 1060.5) return sign_z*4;
  return -999;
};


void llp_MuonSystem_bdt::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  int option;
  if (options == 1){
    option = 1; // used when running condor
  }
  else{
    option = 0;// used when running locally
  }


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
    analysisTag = "Razor2016_80X";

  }
  int wzId;


  const int zh_lepton0_cut = 15;
  const int zh_lepton1_cut = 15;

  const int wh_muonPt_cut = 25;
  const int wh_elePt_cut = 35;
  //
  // if (label == "zH" || label == "bkg_zH" ){
  //   NTrigger = 4;
  //   muonPt_cut = 15;
  //   elePt_cut = 15;
  //   nLepton_cut = 2;
  // }
  // else{
  //   NTrigger = 2;
  //   muonPt_cut = 25;
  //   elePt_cut = 35;
  //   nLepton_cut = 1;
  // }


  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
  LiteTreeMuonSystem *MuonSystem = new LiteTreeMuonSystem;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();
  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *NEvents_genweight = new TH1F("NEvents_genweight", "NEvents_genweight", 1, 1, 2);

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
  if (analysisTag == "Razor2018_17SeptEarlyReReco") helper = new RazorHelper("Razor2018_17SeptEarlyReReco", isData, false);
  else helper = new RazorHelper(analysisTag, isData, false);



  std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
  std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();
  //----------
  //pu histo
  //----------
  //TH1D* puhisto = new TH1D("pileup", "", 50, 0, 50);
  //histogram containing total number of processed events (for normalization)
  //TH1F *histNPV = new TH1F("NPV", "NPV", 2, -0.5, 1.5);
  //TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  //TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
  //TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
  //TH1F *SumPdfWeights = new TH1F("SumPdfWeights", "SumPdfWeights", NUM_PDF_WEIGHTS, -0.5, NUM_PDF_WEIGHTS-0.5);



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

    if (isData)
    {
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    }
    else
    {
      //NEvents->Fill(genWeight);
      MuonSystem->weight = genWeight;
      NEvents->Fill(1);
      NEvents_genweight->Fill(1, genWeight);
    }


    //event info
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
        if (abs(gParticleId[i])== 25)
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
        MuonSystem->nGenParticle++;
        // cout<<"genparticles: "<<MuonSystem->nGenParticle<<endl;


      }
      // cout<<nGenJets<<endl;
      // MuonSystem->nGenJets = 0;
      for(int i=0; i < nGenJets; i++)
      {
        // cout<<genJetE[i]<<","<<MuonSystem->nGenJets<<","<<nGenJets<<endl;
        MuonSystem->genJetE[MuonSystem->nGenJets] = genJetE[i];
        MuonSystem->genJetPt[MuonSystem->nGenJets] = genJetPt[i];
        MuonSystem->genJetEta[MuonSystem->nGenJets] = genJetEta[i];
        MuonSystem->genJetPhi[MuonSystem->nGenJets] = genJetPhi[i];
        MuonSystem->genJetMET[MuonSystem->nGenJets] = genJetMET[i];
        MuonSystem->nGenJets++;
      }
      MuonSystem->genMetPtTrue = genMetPtTrue;
      MuonSystem->genMetPhiTrue = genMetPhiTrue;
      MuonSystem->genMetPtCalo = genMetPtCalo;
      MuonSystem->genMetPhiCalo = genMetPhiCalo;
      for(int i = 0; i < 2;i++)
      {
        MuonSystem->gLLP_eta[i] = gLLP_eta[i];
        MuonSystem->gLLP_phi[i] = gLLP_phi[i];
        MuonSystem->gLLP_decay_vertex_r[i] = sqrt(gLLP_decay_vertex_x[i]*gLLP_decay_vertex_x[i]+gLLP_decay_vertex_y[i]*gLLP_decay_vertex_y[i]);
        MuonSystem->gLLP_decay_vertex_x[i] = gLLP_decay_vertex_x[i];
        MuonSystem->gLLP_decay_vertex_y[i] = gLLP_decay_vertex_y[i];
        MuonSystem->gLLP_decay_vertex_z[i] = gLLP_decay_vertex_z[i];
        float beta = gLLP_beta[i];
        float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[i], 2) + pow(MuonSystem->gLLP_decay_vertex_z[i],2));
        float gamma = 1.0/sqrt(1-beta*beta);
        MuonSystem->gLLP_ctau[i] = gLLP_decay_vertex/(beta * gamma);
        MuonSystem->gLLP_beta[i] = gLLP_beta[i];

        if (abs(MuonSystem->gLLP_eta[i]) < 2.4 && abs(MuonSystem->gLLP_eta[i]) > 0.9
          && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>568
          && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) MuonSystem->gLLP_csc[i] = true;


      }
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
    // if (MuonSystem->met < 200) continue;
    //Triggers
    for(int i = 0; i < NTriggersMAX; i++){
      MuonSystem->HLTDecision[i] = HLTDecision[i];
    }
    // flags
    MuonSystem->Flag_HBHENoiseFilter = Flag_HBHENoiseFilter;
    MuonSystem->Flag_HBHEIsoNoiseFilter = Flag_HBHEIsoNoiseFilter;
    MuonSystem->Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
    MuonSystem->Flag_CSCTightHaloFilter = Flag_CSCTightHaloFilter;
    MuonSystem->Flag_goodVertices = Flag_goodVertices;
    MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;
    MuonSystem->Flag_all = Flag_goodVertices && Flag_HBHEIsoNoiseFilter && Flag_BadPFMuonFilter && Flag_CSCTightHaloFilter && Flag_goodVertices && Flag_ecalBadCalibFilter;
    if (isData && ! MuonSystem->Flag_all) continue;
    //*************************************************************************
    //Start Object Selection
    //*************************************************************************

    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
    for( int i = 0; i < nMuons; i++ )
    {

      if(!isMuonPOGLooseMuon(i)) continue;
      if(muonPt[i] < zh_lepton1_cut) continue;
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
      tmpMuon.passId = isMuonPOGTightMuon(i);
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
      for( uint j = 0; j < Leptons.size(); j++ )
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

    if (foundZ && fabs(ZMass-Z_MASS) < 30.0 && Leptons.size() == 2 && leadingLepPt > zh_lepton0_cut)
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
    for ( auto &tmp : Leptons )
    {
      MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
      MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
      MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
      MuonSystem->lepPassId[MuonSystem->nLeptons] = tmp.passId;
      MuonSystem->nLeptons++;
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
    double jetCorrPt = jetPt[i]*JEC;
    double jetCorrE = jetE[i]*JEC;
    // cout<<JEC<<endl;
    // cout<<"corrected pt: "<<jetCorrPt<<", "<<"eta: "<<jetEta[i]<<endl;
    TLorentzVector thisJet = makeTLorentzVector( jetCorrPt, jetEta[i], jetPhi[i], jetCorrE );

    if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
    if( fabs( thisJet.Eta() ) >= 3.0 ) continue;
    if ( !jetPassIDLoose[i] ) continue;

    jets tmpJet;
    tmpJet.jet    = thisJet;
    tmpJet.passId = jetPassIDTight[i];

    tmpJet.jetPassMuFrac = jetPassMuFrac[i];
    tmpJet.jetNeutralHadronEnergyFraction = jetNeutralHadronEnergyFraction[i];
    tmpJet.jetNeutralEMEnergyFraction = jetNeutralEMEnergyFraction[i];
    tmpJet.jetChargedEMEnergyFraction = jetChargedEMEnergyFraction[i];
    tmpJet.jetChargedHadronEnergyFraction = jetChargedHadronEnergyFraction[i];

    // calculate jet energy scale uncertainty
    double unc = helper->getJecUnc( jetCorrPt, jetEta[i], runNum ); //use run=999 as default
    double jetPtJESUp = jetCorrPt*(1+unc);
    double jetPtJESDown = jetCorrPt/(1+unc);
    double jetEJESUp = jetCorrE*(1+unc);
    double jetEJESDown = jetCorrE/(1+unc);
    TLorentzVector thisJetJESUp = makeTLorentzVector(jetPtJESUp, jetEta[i], jetPhi[i], jetEJESUp);
    TLorentzVector thisJetJESDown = makeTLorentzVector(jetPtJESDown, jetEta[i], jetPhi[i], jetEJESDown);
    if (jetPtJESUp > 20) MuonSystem->metJESUp +=  -1 * (thisJetJESUp.Pt() - thisJet.Pt());
    if (jetPtJESDown > 20) MuonSystem->metJESDown +=  -1 * (thisJetJESDown.Pt() - thisJet.Pt());


    Jets.push_back(tmpJet);

    }

    //-----------------------------
    //Require at least 2 jets
    //-----------------------------
    // if( Jets.size() < 2 ) continue;
    // if (triggered) trig_lepId_dijet->Fill(1);
    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);
    if (Jets.size()>0)
    {
      MuonSystem->jetMet_dPhi = RazorAnalyzer::deltaPhi(jetPhi[0],metType1Phi);
    }
    else{
      MuonSystem->jetMet_dPhi = -999.;
    }
    double jetMet_dPhiMin_temp = 999.;
    double jetMet_dPhiMin4_temp = 999.;
    for ( auto &tmp : Jets )
    {
      MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
      MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
      MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
      MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
      MuonSystem->jetTime[MuonSystem->nJets] = tmp.time;
      if (jetMet_dPhiMin4_temp > abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi)) && MuonSystem->nJets < 4)
      {
        jetMet_dPhiMin4_temp = abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi));

      }
      if (jetMet_dPhiMin_temp > abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi)))
      {
        jetMet_dPhiMin_temp = abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi));

      }
      MuonSystem->jetPassId[MuonSystem->nJets] = tmp.passId;
      // MuonSystem->ecalNRechits[MuonSystem->nJets] = tmp.ecalNRechits;
      // MuonSystem->ecalNRechits[MuonSystem->nJets] = tmp.ecalRechitE;
      MuonSystem->jetChargedEMEnergyFraction[MuonSystem->nJets] = tmp.jetChargedEMEnergyFraction;
      MuonSystem->jetNeutralEMEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralEMEnergyFraction;
      MuonSystem->jetChargedHadronEnergyFraction[MuonSystem->nJets] = tmp.jetChargedHadronEnergyFraction;
      MuonSystem->jetNeutralHadronEnergyFraction[MuonSystem->nJets] = tmp.jetNeutralHadronEnergyFraction;
      MuonSystem->jetPassMuFrac[MuonSystem->nJets] = tmp.jetPassMuFrac;
      float min_deltaR = 15.;
      int index = 999;
      // cout<<"nGenJets"<<nGenJets<<endl;
      for(int i=0; i < nGenJets; i++)
      {

        double current_delta_r = RazorAnalyzer::deltaPhi(genJetPhi[i],jetPhi[MuonSystem->nJets]);
        if (current_delta_r < min_deltaR)
        {
          min_deltaR = current_delta_r;
          index = i;
        }
      }
      if (min_deltaR < 0.4)
      {
        MuonSystem->jet_match_genJet_minDeltaR[MuonSystem->nJets] = min_deltaR;
        MuonSystem->jet_match_genJet_index[MuonSystem->nJets] = index;
        MuonSystem->jet_match_genJet_pt[MuonSystem->nJets] = jetPt[index];
      }
      MuonSystem->nJets++;
    }


    MuonSystem-> jetMet_dPhiMin = jetMet_dPhiMin_temp;
    MuonSystem-> jetMet_dPhiMin4 = jetMet_dPhiMin4_temp;

    //-----------------------------
    // CSC INFO
    //-----------------------------
    //if(nCscClusters == 0) continue;

    /*if(nCscRechitClusters == 0) continue;

    for(int i = 0; i < nCscRechitClusters; i++)
    {
      if (isData && cscRechitClusterTime[i] > -12.5)continue;
      if (abs(cscRechitClusterMaxChamber[i]) <= 12) continue;

      MuonSystem->cscRechitClusterX[MuonSystem->nCscRechitClusters] =cscRechitClusterX[i];
      MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] =cscRechitClusterY[i];
      MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] =cscRechitClusterZ[i];
      MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] =cscRechitClusterTime[i];
      MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] =cscRechitClusterEta[i];
      MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters] =cscRechitClusterPhi[i];
      MuonSystem->cscRechitClusterMajorAxis[MuonSystem->nCscRechitClusters] =cscRechitClusterMajorAxis[i];
      MuonSystem->cscRechitClusterMinorAxis[MuonSystem->nCscRechitClusters] =cscRechitClusterMinorAxis[i];
      MuonSystem->cscRechitClusterXSpread[MuonSystem->nCscRechitClusters] =cscRechitClusterXSpread[i];
      MuonSystem->cscRechitClusterYSpread[MuonSystem->nCscRechitClusters] =cscRechitClusterYSpread[i];
      MuonSystem->cscRechitClusterZSpread[MuonSystem->nCscRechitClusters] =cscRechitClusterZSpread[i];
      MuonSystem->cscRechitClusterEtaPhiSpread[MuonSystem->nCscRechitClusters] =cscRechitClusterEtaPhiSpread[i];
      MuonSystem->cscRechitClusterEtaSpread[MuonSystem->nCscRechitClusters] =cscRechitClusterEtaSpread[i];
      MuonSystem->cscRechitClusterPhiSpread[MuonSystem->nCscRechitClusters] = cscRechitClusterPhiSpread[i];
      MuonSystem->cscRechitClusterTimeSpread[MuonSystem->nCscRechitClusters] =cscRechitClusterTimeSpread[i];

      MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] =cscRechitClusterJetVetoPt[i];
      MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] =cscRechitClusterJetVetoE[i];
      MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] =cscRechitClusterMuonVetoPt[i];
      MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] =cscRechitClusterMuonVetoE[i];

      MuonSystem->cscRechitClusterSize[MuonSystem->nCscRechitClusters] = cscRechitClusterSize[i];
      MuonSystem->cscRechitClusterMaxChamber[MuonSystem->nCscRechitClusters] =cscRechitClusterMaxChamber[i];
      MuonSystem->cscRechitClusterMaxChamberRatio[MuonSystem->nCscRechitClusters] = cscRechitClusterMaxChamberRatio[i];
      MuonSystem->cscRechitClusterNChamber[MuonSystem->nCscRechitClusters] =cscRechitClusterNChamber[i];
      MuonSystem->cscRechitClusterMaxStation[MuonSystem->nCscRechitClusters] = cscRechitClusterMaxStation[i];
      MuonSystem->cscRechitClusterMaxStationRatio[MuonSystem->nCscRechitClusters] = cscRechitClusterMaxStationRatio[i];
      MuonSystem->cscRechitClusterNStation[MuonSystem->nCscRechitClusters] = cscRechitClusterNStation[i];
      MuonSystem->cscRechitClusterMe11Ratio[MuonSystem->nCscRechitClusters] =cscRechitClusterMe11Ratio[i];
      MuonSystem->cscRechitClusterMe12Ratio[MuonSystem->nCscRechitClusters] = cscRechitClusterMe12Ratio[i];

      MuonSystem->cscRechitClusterNRechitChamberPlus11[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus11[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus12[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus12[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus13[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus13[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus21[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus21[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus22[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus22[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus31[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus31[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus32[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus32[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus41[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus41[i];
      MuonSystem->cscRechitClusterNRechitChamberPlus42[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberPlus42[i];

      MuonSystem->cscRechitClusterNRechitChamberMinus11[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus11[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus12[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus12[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus13[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus13[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus21[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus21[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus22[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus22[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus31[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus31[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus32[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus32[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus41[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus41[i];
      MuonSystem->cscRechitClusterNRechitChamberMinus42[MuonSystem->nCscRechitClusters] = cscRechitClusterNRechitChamberMinus42[i];
      bool me1112_veto = abs(MuonSystem->cscRechitClusterMaxChamber[MuonSystem->nCscRechitClusters]) > 12;
      // if (MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] < JET_PT_CUT) MuonSystem->nCsc_JetVetoCluster0p4 += cscRechitClusterSize[i];
      // if (MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] < JET_PT_CUT && MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] < MUON_PT_CUT) MuonSystem->nCsc_JetMuonVetoCluster0p4 += cscRechitClusterSize[i];
      // if (MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] < JET_PT_CUT && me1112_veto) MuonSystem->nCsc_JetVetoCluster0p4_Me1112Veto+= cscRechitClusterSize[i];
      if (MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] < JET_PT_CUT && MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] < MUON_PT_CUT && me1112_veto) MuonSystem->nCsc_JetMuonVetoRechitCluster0p4_Me1112Veto+= cscRechitClusterSize[i];
      //
      float min_deltaR = 15.;
      int index = 999;
      for(int j = 0; j < 2;j++)
      {

        double current_delta_r = RazorAnalyzer::deltaR(cscRechitClusterEta[MuonSystem->nCscRechitClusters], cscRechitClusterPhi[MuonSystem->nCscRechitClusters], gLLP_eta[j], gLLP_phi[j]);
        if (current_delta_r < min_deltaR)
        {
          min_deltaR = current_delta_r;
          index = j;
        }
      }
      if (min_deltaR < 0.4)
      {
        MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = true;
        MuonSystem->cscRechitCluster_match_gLLP_minDeltaR[MuonSystem->nCscRechitClusters] = min_deltaR;
        MuonSystem->cscRechitCluster_match_gLLP_index[MuonSystem->nCscRechitClusters] = index;
      }
      else{
        if (!isData) continue;
      }
      MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhi);

      MuonSystem->nCscRechitClusters++;
    }*/

    if (MuonSystem->nCscRechitClusters == 0) continue;
    MuonSystem->tree_->Fill();
  }

    cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();
}
