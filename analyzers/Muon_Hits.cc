#include "Muon_Hits.h"

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

#include "TNtuple.h"

using namespace std::chrono;
using namespace std;

struct greater_than_pt {
   inline bool operator()(const TLorentzVector & p1,
      const TLorentzVector & p2) {
      return p1.Pt() > p2.Pt();
   }
};

struct leptons {
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

struct jets {
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

};

//lepton highest pt comparator
struct largest_pt {
   inline bool operator()(const leptons & p1,
      const leptons & p2) {
      return p1.lepton.Pt() > p2.lepton.Pt();
   }
}
my_largest_pt;

//jet highest pt comparator
struct largest_pt_jet {
   inline bool operator()(const jets & p1,
      const jets & p2) {
      return p1.jet.Pt() > p2.jet.Pt();
   }
}
my_largest_pt_jet;
int N_Good_Events = 0;
void Muon_Hits::Analyze(bool isData, int options, string outputfilename, string analysisTag) {
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

   bool signalScan = int(options / 10) == 1;
   int option = options % 10;
   // if (options % 1){
   //   option = 1; // used when running condor
   // }
   // else{
   //   option = 0;// used when running locally
   // }

   if (isData) {
      std::cout << "[INFO]: running on data with option: " << option << std::endl;
   } else {
      std::cout << "[INFO]: running on MC with option: " << option << std::endl;
   }
   if (signalScan) {
      std::cout << "[INFO]: running with Signal scan" << std::endl;
   } else {
      std::cout << "[INFO]: running without Signal scan " << option << std::endl;
   }

   const float ELE_MASS = 0.000511;
   const float MU_MASS = 0.105658;
   const float Z_MASS = 91.2;

   if (analysisTag == "") {
      analysisTag = "Razor2016_80X";
   }
   int wzId;

   const int zh_lepton0_cut = 15;
   const int zh_lepton1_cut = 15;

   const int wh_muonPt_cut = 25;
   const int wh_elePt_cut = 18; //5 35;

   //-----------------------------------------------
   //Set up Output File
   //-----------------------------------------------
   string outfilename = outputfilename;
   // if (outfilename == "") outfilename = "/eos/user/a/arhayrap/OUTPUT_DATA/Muon_Hits.root";
   if (outfilename == "") outfilename = "./Muon_Hits.root";
   TFile * outFile;
   if (isData || !signalScan) outFile = new TFile(outfilename.c_str(), "RECREATE");

   LiteTreeMuonSystem * MuonSystem = new LiteTreeMuonSystem;
   MuonSystem -> CreateTree();
   MuonSystem -> tree_ -> SetAutoFlush(0);
   MuonSystem -> InitTree();

   // for signals, need one output file for each signal point
   map < pair < int, int > , TFile * > Files2D;
   map < pair < int, int > , TTree * > Trees2D;
   map < pair < int, int > , TH1F * > NEvents2D;
   map < pair < int, int > , TH1F * > accep2D;
   map < pair < int, int > , TH1F * > accep_met2D;
   map < pair < int, int > , TH1F * > Total2D;

   // map<pair<int,int>, TH1F*> smsSumScaleWeights2D;
   // map<pair<int,int>, TH1F*> smsSumPdfWeights2D;
   // map<pair<int,int>, TH1F*> smsNISRJets2D;
   // map<pair<int,int>, TH1F*> smsPtISR2D;
   // map<pair<int,int>, TH1F*> smsNPV2D;

   //histogram containing total number of processed events (for normalization)
   TH1F * NTriggers = new TH1F("NTriggers", "", 1001, 0, 1001);
   TH1F * Matched_PT = new TH1F("Matched_PT", "", 800, 0, 400);
   TH1F * tau_pt = new TH1F("tau_pt", "", 800, 0, 400);
   TH1F * NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
   TH1F * lepton_and_hadron = new TH1F("lepton_and_hadron", "lepton_and_hadron", 10, -5, 5);
   TH1F * gllp_custer_dr = new TH1F("gllp_custer_dr", "gllp_custer_dr", 1000, 0, 10);
   TH1F * ds_csc_cluster = new TH1F("ds_csc_cluster", "ds_csc_cluster", 100, 0, 5);
   TH1F * csc_cluster = new TH1F("csc_cluster", "csc_cluster", 100, 0, 5);
   TH1F * Muon_track_match = new TH1F("Muon_track_match", "Muon_track_match", 50, -5, 5);
   TH1F * accep = new TH1F("accep", "acceptance", 1, 1, 2);
   TH1F * accep_met = new TH1F("accep_met", "acceptance_met", 1, 1, 2);

   TNtuple * EVENTS         = new TNtuple("EVENTS"         , "EVENTS"         , "id");
   TNtuple * RECHITS        = new TNtuple("RECHITS"        , "RECHITS"        , "id:eta:phi:x:y:z:station:which");
   TNtuple * SEGMENTS       = new TNtuple("SEGMENTS"       , "SEGMENTS"       , "id:eta:phi:x:y:z:station:which");
   TNtuple * RECHITS_CL     = new TNtuple("RECHITS_CL"     , "RECHITS_CL"     , "id:eta:phi:x:y:z:time:size:stations:avg_stations:which");
   TNtuple * RECHITS_CL_ADD = new TNtuple("RECHITS_CL_ADD" , "RECHITS_CL_ADD" , "cscJetVetoPt:cscMuonVetoPt:cscNRechitChamberPlus11:cscNRechitChamberPlus12:cscNRechitChamberMinus11:cscNRechitChamberMinus12");
   TNtuple * SEGMENTS_CL    = new TNtuple("SEGMENTS_CL"    , "SEGMENTS_CL"    , "id:eta:phi:x:y:z:time:size:stations:avg_stations:which");

   TNtuple * Neutrino_tau = new TNtuple("Neutrino_tau", "Neutrino_tau", "pt:eta:phi:e");
   TNtuple * Leptons_other = new TNtuple("Leptons_other", "Leptons_other", "pt:eta:phi:e");
   TNtuple * Neutrino_other = new TNtuple("Neutrino_other", "Neutrino_other", "pt:eta:phi:e");
   TNtuple * Hadrons = new TNtuple("Hadrons", "Hadrons", "pt:eta:phi:e"); // hadrons
   TNtuple * Leptons = new TNtuple("Leptons", "Leptons", "pt:eta:phi:e"); // leptons

   TH1F * Nmet200 = new TH1F("Nmet200", "Nmet200", 1, 1, 2);
   TH1F * NmetFilter = new TH1F("NmetFilter", "NmetFilter", 1, 1, 2);
   TH1F * Nlep0 = new TH1F("Nlep0", "Nlep0", 1, 1, 2);
   TH1F * Njet1 = new TH1F("Njet1", "Njet1", 1, 1, 2);
   TH1F * NcosmicVeto = new TH1F("NcosmicVeto", "NcosmicVeto", 1, 1, 2);
   TH1F * NEvents_genweight = new TH1F("NEvents_genweight", "NEvents_genweight", 1, 1, 2);

   /*
     //histogram containing total number of processed events (for normalization)
     TH1F * NTriggers = new TH1F("NTriggers", "", 982, 0, 982);
     TH1F * NMatched_Triggers = new TH1F("NMatched_Triggers", "", 982, 0, 982);
     TH1F * NMatched = new TH1F("NMatched", "", 100, 0, 20);
     TH1F * Matched_PT = new TH1F("Matched_PT", "", 800, 0, 400);
     TH1F * tau_pt = new TH1F("tau_pt", "", 800, 0, 400);
     TH1F * NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
     TH1F * accep = new TH1F("accep", "acceptance", 1, 1, 2);
     TH1F * accep_met = new TH1F("accep_met", "acceptance_met", 1, 1, 2);
     //TH1F * N_selected_events = new TH1F("N_selected_events", "N_selected_events", 1, 1, 2);
     //TH1F * Ntau_motherID     = new TH1F("Ntau_motherID", "Ntau_motherID", 250000, -500000, 500000);

     TH1F * Nmet200 = new TH1F("Nmet200", "Nmet200", 1, 1, 2);
     TH1F * NmetFilter = new TH1F("NmetFilter", "NmetFilter", 1, 1, 2);
     TH1F * Nlep0 = new TH1F("Nlep0", "Nlep0", 1, 1, 2);
     TH1F * Njet1 = new TH1F("Njet1", "Njet1", 1, 1, 2);
     TH1F * NcosmicVeto = new TH1F("NcosmicVeto", "NcosmicVeto", 1, 1, 2);
   */

   double drmc = 0.4, drmt = 0.3, iso = 0.3, rcone = 0.3, dzpv = 0.5, dxybs = 0.5, min_pt_track = 50, min_pt_muon = 50, min_clustering_rechits = 50, min_clustering_segments = 10;

   char * cmsswPath;
   cmsswPath = getenv("CMSSW_BASE");
   string pathname;
   if (cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
   if (cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1

   //--------------------------------
   //Initialize helper
   //--------------------------------
   RazorHelper * helper = 0;
   helper = new RazorHelper(analysisTag, isData, false);

   std::vector < FactorizedJetCorrector * > JetCorrector = helper -> getJetCorrector();
   std::vector < std::pair < int, int > > JetCorrectorIOV = helper -> getJetCorrectionsIOV();
   //----------//
   //pu histo  //
   //----------//

   int count_gllp = 0;
   int count_stau = 0;

   //*************************************************************************
   //Look over Input File Events
   //*************************************************************************
   if (fChain == 0) return;
   cout << "Total Events: " << fChain -> GetEntries() << "\n";
   Long64_t nbytes = 0, nb = 0;
   clock_t start, end;
   start = clock();
   for (Long64_t jentry = 0; jentry < fChain -> GetEntries(); jentry++) {
      //begin event
      cout << "-------------------------------------- //|event began//| --------------------------------------" << endl;

      if (jentry % 10000 == 0) {
         end = clock();
         double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
         start = clock();
      }
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      //GetEntry(ientry);
      nb = fChain -> GetEntry(jentry);
      nbytes += nb;

      //fill normalization histogram
      MuonSystem -> InitVariables();

      if (!isData && signalScan) {

         string mh_substring = lheComments -> substr(lheComments -> find("MH-") + 3);
         int mh = stoi(mh_substring.substr(0, mh_substring.find('_')));
         string mx_substring = lheComments -> substr(lheComments -> find("MS-") + 3);
         int mx = stoi(mx_substring.substr(0, mx_substring.find('_')));
         string ctau_substring = lheComments -> substr(lheComments -> find("ctauS-") + 6);
         int ctau = stoi(ctau_substring.substr(0, ctau_substring.find('_')));
         // MuonSystem -> mH = mh;
         // MuonSystem -> mX = mx;
         // MuonSystem -> ctau = ctau;
         // if (mh2 != mh || mx2!=mx || ctau2!=ctau) continue;
         cout << * lheComments << endl;

         pair < int, int > signalPair = make_pair(mx, ctau);

         if (Files2D.count(signalPair) == 0) { //create file and tree
            //format file name
            string thisFileName = outfilename;
            thisFileName.erase(thisFileName.end() - 5, thisFileName.end());
            thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";

            Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
            Trees2D[signalPair] = MuonSystem -> tree_ -> CloneTree(0);
            NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1, 0.5, 1.5);
            Total2D[signalPair] = new TH1F(Form("Total%d%d", mx, ctau), "Total", 1, 0.5, 1.5);
            accep2D[signalPair] = new TH1F(Form("accep2D%d%d", mx, ctau), "acceptance", 1, 0.5, 1.5);
            accep_met2D[signalPair] = new TH1F(Form("accep_met2D%d%d", mx, ctau), "acceptance_met", 1, 0.5, 1.5);

            cout << "Created new output file " << thisFileName << endl;
         }
         //Fill NEvents hist
         NEvents2D[signalPair] -> Fill(1.0, genWeight);

      }

      //event info
      if (isData) {
         NEvents -> Fill(1);
         // MuonSystem -> weight = 1;
      } else {
         // cout<<*lheComments<<endl;
         // MuonSystem -> weight = genWeight;
         NEvents -> Fill(1, genWeight);
         NEvents_genweight -> Fill(1);
      }
      // MuonSystem -> runNum = runNum;
      // MuonSystem -> lumiSec = lumiNum;
      // MuonSystem -> evtNum = eventNum;
      //tau_
      bool hnl_model = true;
      bool wzFlag = false;
      if (!isData) {
         cout << "                  *************/////*************                                 " << nGenParticle << endl;
         for (int i = 0; i < nGenParticle; i++) {
            if (abs(gParticleId[i]) == 9900012 || abs(gParticleId[i]) == 9900014 || abs(gParticleId[i]) == 9900015 || abs(gParticleId[i]) == 9900016) hnl_model = true;
            if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && gParticleStatus[i] == 1 && abs(gParticleMotherId[i]) == 24) { // choosing only the W->munu events
               wzFlag = true;
               // MuonSystem -> gLepId = gParticleId[i];
               // MuonSystem -> gLepPt = gParticlePt[i];
               // MuonSystem -> gLepEta = gParticleEta[i];
               // MuonSystem -> gLepE = gParticleE[i];
               // MuonSystem -> gLepPhi = gParticlePhi[i];
            } else if (abs(gParticleId[i]) == 15 && gParticleStatus[i] == 2 && abs(gParticleMotherId[i]) == 24) {
               wzFlag = true;
               // MuonSystem -> gLepId = gParticleId[i];
               // MuonSystem -> gLepPt = gParticlePt[i];
               // MuonSystem -> gLepEta = gParticleEta[i];
               // MuonSystem -> gLepE = gParticleE[i];
               // MuonSystem -> gLepPhi = gParticlePhi[i];
            }
            /*
            if (abs(gParticleId[i]) == 24 && gParticleStatus[i] == 62) {
               MuonSystem -> gWPt = gParticlePt[i];
            }
            if (abs(gParticleId[i]) == 25 || abs(gParticleId[i] == 35)) {
               MuonSystem -> gHiggsPt = gParticlePt[i];
               MuonSystem -> gHiggsEta = gParticleEta[i];
               MuonSystem -> gHiggsPhi = gParticlePhi[i];
               MuonSystem -> gHiggsE = gParticleE[i];
            }
            */
            /*
            MuonSystem -> gParticleStatus[MuonSystem -> nGenParticle] = gParticleStatus[i];
            MuonSystem -> gParticleId[MuonSystem -> nGenParticle] = gParticleId[i];
            MuonSystem -> gParticleMotherId[MuonSystem -> nGenParticle] = gParticleMotherId[i];
            MuonSystem -> gParticlePt[MuonSystem -> nGenParticle] = gParticlePt[i];
            MuonSystem -> gParticleEta[MuonSystem -> nGenParticle] = gParticleEta[i];
            MuonSystem -> gParticlePhi[MuonSystem -> nGenParticle] = gParticlePhi[i];
            MuonSystem -> gParticleE[MuonSystem -> nGenParticle] = gParticleE[i];
            MuonSystem -> nGenParticle++;
            */
         }
         MuonSystem -> higgsPtWeight = helper -> getHiggsPtWeight(MuonSystem -> gHiggsPt);
         //for (unsigned int i = 0; i < 9; i++) {
         // MuonSystem -> higgsPtWeightSys[i] = helper -> getHiggsPtWeightSys(MuonSystem -> gHiggsPt, i) / MuonSystem -> higgsPtWeight;
         //MuonSystem -> scaleWeights[i] = ( * scaleWeights)[i] / genWeight;
         //}
         /*
         MuonSystem -> sf_facScaleUp = MuonSystem -> higgsPtWeightSys[5];
         MuonSystem -> sf_facScaleDown = MuonSystem -> higgsPtWeightSys[3];
         MuonSystem -> sf_renScaleUp = MuonSystem -> higgsPtWeightSys[7];
         MuonSystem -> sf_renScaleDown = MuonSystem -> higgsPtWeightSys[1];
         MuonSystem -> sf_facRenScaleUp = MuonSystem -> higgsPtWeightSys[8];
         MuonSystem -> sf_facRenScaleDown = MuonSystem -> higgsPtWeightSys[0];

         MuonSystem -> genMetPtTrue = genMetPtTrue;
         MuonSystem -> genMetPhiTrue = genMetPhiTrue;
         MuonSystem -> genMetPtCalo = genMetPtCalo;
         MuonSystem -> genMetPhiCalo = genMetPhiCalo;
         for (int i = 0; i < 2; i++) {
            cout << "GLLP ->    " << gLLP_pt[i] << "    " << gLLP_eta[i] << "    " << gLLP_phi[i] << endl;
            MuonSystem -> gLLP_eta[i] = gLLP_eta[i];
            MuonSystem -> gLLP_phi[i] = gLLP_phi[i];
            MuonSystem -> gLLP_decay_vertex_r[i] = sqrt(gLLP_decay_vertex_x[i] * gLLP_decay_vertex_x[i] + gLLP_decay_vertex_y[i] * gLLP_decay_vertex_y[i]);
            MuonSystem -> gLLP_decay_vertex_x[i] = gLLP_decay_vertex_x[i];
            MuonSystem -> gLLP_decay_vertex_y[i] = gLLP_decay_vertex_y[i];
            MuonSystem -> gLLP_decay_vertex_z[i] = gLLP_decay_vertex_z[i];
            float beta = gLLP_beta[i];
            float gLLP_decay_vertex = sqrt(pow(MuonSystem -> gLLP_decay_vertex_r[i], 2) + pow(MuonSystem -> gLLP_decay_vertex_z[i], 2));
            float gamma = 1.0 / sqrt(1 - beta * beta);
            MuonSystem -> gLLP_ctau[i] = gLLP_decay_vertex / (beta * gamma);
            MuonSystem -> gLLP_beta[i] = gLLP_beta[i];
            MuonSystem -> gLLP_e[i] = gLLP_e[i];
            MuonSystem -> gLLP_pt[i] = gLLP_pt[i];
            MuonSystem -> gLLP_lepdPhi[i] = deltaPhi(MuonSystem -> gLepPhi, MuonSystem -> gLLP_phi[i]);
            if (hnl_model) {
               MuonSystem -> gLLP_EMFracP[i] = gLLP_daughter_pt[0] * cosh(gLLP_daughter_eta[0]) / (gLLP_pt[i] * cosh(gLLP_eta[i]));
               MuonSystem -> gLLP_EMFracPz[i] = gLLP_daughter_pt[0] * sinh(gLLP_daughter_eta[0]) / (gLLP_pt[i] * sinh(gLLP_eta[i]));
               MuonSystem -> gLLP_EMFracE[i] = gLLP_daughter_e[0] / (gLLP_e[i]);
               MuonSystem -> gLLP_EMFracEz[i] = gLLP_daughter_e[0] / cosh(gLLP_daughter_eta[0]) * sinh(gLLP_daughter_eta[0]) / (gLLP_e[i] / cosh(gLLP_eta[i]) * sinh(gLLP_eta[i]));

               MuonSystem -> gLLP_visE[i] = gLLP_e[i];
               MuonSystem -> gLLP_visEz[i] = gLLP_e[i] / cosh(gLLP_eta[i]) * sinh(gLLP_eta[i]);
               MuonSystem -> gLLP_visP[i] = gLLP_pt[i] * cosh(gLLP_eta[i]);
               MuonSystem -> gLLP_visPz[i] = gLLP_pt[i] * sinh(gLLP_eta[i]);
            }

            // cout<<MuonSystem->gLLP_e[i]/gamma<<endl;
            
            if (abs(MuonSystem -> gLLP_eta[i]) < 2.4 && abs(MuonSystem -> gLLP_eta[i]) > 0.9 &&
               abs(MuonSystem -> gLLP_decay_vertex_z[i]) < 1100 && abs(MuonSystem -> gLLP_decay_vertex_z[i]) > 568 &&
               MuonSystem -> gLLP_decay_vertex_r[i] < 695.5) MuonSystem -> gLLP_csc[i] = true;
            else MuonSystem -> gLLP_csc[i] = false;
            */

            /*
            if (abs(MuonSystem -> gLLP_eta[i]) < 2.4 &&
               abs(MuonSystem -> gLLP_decay_vertex_z[i]) < 1100 && abs(MuonSystem -> gLLP_decay_vertex_z[i]) > 400 &&
               MuonSystem -> gLLP_decay_vertex_r[i] < 695.5) MuonSystem -> gLLP_csc[i] = true;
            */
            
        /*
         }
         MuonSystem -> gLLP_deltaR = deltaR(gLLP_eta[0], gLLP_phi[0], gLLP_eta[1], gLLP_phi[1]);

         for (int i = 0; i < 4; i++) {
            MuonSystem -> gLLP_daughter_id[i] = gLLP_daughter_id[i];
            MuonSystem -> gLLP_daughter_pt[i] = gLLP_daughter_pt[i];
            MuonSystem -> gLLP_daughter_eta[i] = gLLP_daughter_eta[i];
            MuonSystem -> gLLP_daughter_phi[i] = gLLP_daughter_phi[i];
            MuonSystem -> gLLP_daughter_e[i] = gLLP_daughter_e[i];
            MuonSystem -> gLLP_daughter_mass[i] = gLLP_daughter_mass[i];
         }
         */
         /*
               HLTDecision 136 IsoMu_24

               HLTDecision 137 IsoMu_26

               HLTDecision 646 IsoMu_30
         */

         // RECHIT CLUSTERING FOR DT ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
         vector < Point > points;
         points.clear();

         for (int i = 0; i < nDtRechits; i++) {
            Point p;
            p.phi = dtRechitPhi[i];
            p.eta = dtRechitEta[i];
            p.x = dtRechitX[i];
            p.y = dtRechitY[i];
            p.z = dtRechitZ[i];
            p.t = dtRechitTime[i];
            // p.twire = dtRechitTime[i];
            p.chamber = dtRechitWheel[i];
            p.station = dtRechitStation[i];
            p.clusterID = UNCLASSIFIED;
            points.push_back(p);
         }
         //Do DBSCAN Clustering
         int min_point_dt = min_clustering_rechits; //minimum number of segments to call it a cluster
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

         points.clear();
         MuonSystem -> nDtRechitClusters = 0;

         for (auto & tmp: ds_dtRechit.CscCluster) {
            cout << "DT clustering -> " << tmp.eta << "  " << tmp.phi << "  " << tmp.x << "  " << tmp.y << "  " << tmp.z << "  " << tmp.t << endl;
         }

        // ------------------------------------------------------- DT RECHIT CLUSTERS TO RPC MATCHING -------------------------------------------------------
        std::vector<int> dtRechitCluster_match_rpcBx;
        for (auto & tmp: ds_dtRechit.CscCluster) {
          cout<<"nRpc "<<nRpc<<endl;
          for (int i = 0; i < nRpc; i++) {
            cout<<"CHECKING: "<<rpcPhi[i]<<"   "<<rpcRegion[i]<<"  "<<rpcRing[i]<<" <---> "<< tmp.nChamber<<endl;
            float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
            if (rpcRegion[i] != 0) continue;
            if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], tmp.phi)) < 0.5)
            {
              if (rpcRing[i] == tmp.nChamber)
              {
                dtRechitCluster_match_rpcBx.push_back(rpcBx[i]);
                cout<<"Rechits rpcBx[i]: "<<rpcBx[i]<<endl;
              }
            }
          }
        }
        cout<<endl<<endl;

         // RECHIT CLUSTERING FOR CSC ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

         for (int i = 0; i < ncscRechits; i++) {
            Point p;
            p.phi = cscRechitsPhi[i];
            p.eta = cscRechitsEta[i];
            p.x = cscRechitsX[i];
            p.y = cscRechitsY[i];
            p.z = cscRechitsZ[i];
            p.t = cscRechitsTpeak[i];
            p.twire = cscRechitsTwire[i];
            p.chamber = cscRechitsChamber[i];
            p.station = cscRechitsStation[i];
            p.clusterID = UNCLASSIFIED;
            points.push_back(p);
         }

         //Do DBSCAN Clustering
         int min_point_csc = min_clustering_rechits; //minimum number of segments to call it a cluster
         float epsilon_csc = 0.2; //cluster radius parameter
         DBSCAN ds_cscRechit(min_point_csc, epsilon_csc, points);
         ds_cscRechit.run();
         ds_cscRechit.result();
         ds_cscRechit.clusterMoments();
         ds_cscRechit.sort_clusters();
         ds_cscRechit.merge_clusters();
         ds_cscRechit.result();
         ds_cscRechit.clusterMoments();
         ds_cscRechit.sort_clusters();

         points.clear();

         for (auto & tmp: ds_cscRechit.CscCluster) {
            cout << "CSC clustering -> " << tmp.eta << "  " << tmp.phi << "  " << tmp.x << "  " << tmp.y << "  " << tmp.z << "  " << tmp.t << endl;
         }

         // SEGMENT CLUSTERING FOR DT ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
         // vector<Point> points;
         points.clear();

         for (int i = 0; i < nDtSeg; i++) {
            Point p;
            p.phi = dtSegPhi[i];
            p.eta = dtSegEta[i];
            p.x = dtSegX[i];
            p.y = dtSegY[i];
            p.z = dtSegZ[i];
            p.t = dtSegTime[i];
            // p.twire = dtSegTime[i];
            p.chamber = dtSegWheel[i];
            p.station = dtSegStation[i];
            p.clusterID = UNCLASSIFIED;
            points.push_back(p);
         }
         //Do DBSCAN Clustering
         min_point_dt = min_clustering_segments; //minimum number of segments to call it a cluster
         epsilon_dt = 0.2; //cluster radius parameter
         DBSCAN ds_dtSeg(min_point_dt, epsilon_dt, points);
         ds_dtSeg.run();
         ds_dtSeg.result();
         ds_dtSeg.clusterMoments();
         ds_dtSeg.sort_clusters();
         ds_dtSeg.merge_clusters();
         ds_dtSeg.result();
         ds_dtSeg.clusterMoments();
         ds_dtSeg.sort_clusters();

         points.clear();
         // MuonSystem->nDtSegClusters = 0;

         for (auto & tmp: ds_dtSeg.CscCluster) {
            cout << "DT SEGMENT clustering -> " << tmp.eta << "  " << tmp.phi << "  " << tmp.x << "  " << tmp.y << "  " << tmp.z << "  " << tmp.t << endl;
         }

        // ------------------------------------------------------- DT SEGMENT CLUSTERS TO RPC MATCHING -------------------------------------------------------
        std::vector<int> dtSegmentCluster_match_rpcBx;
        for (auto & tmp: ds_dtSeg.CscCluster) {
          for (int i = 0; i < nRpc; i++) {
            float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
            if (rpcRegion[i] != 0) continue;
            if (abs(RazorAnalyzer::deltaPhi(rpcPhi[i], tmp.phi)) < 0.5)
            {
              if (rpcRing[i] == tmp.nChamber)
              {
                dtSegmentCluster_match_rpcBx.push_back(rpcBx[i]);
                cout<<"Segments rpcBx[i]: "<<rpcBx[i]<<endl;
              }
            }
          }
        }

         // SEGMENT CLUSTERING FOR CSC ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
         for (int i = 0; i < nCscSeg; i++) {
            Point p;
            p.phi = cscSegPhi[i];
            p.eta = cscSegEta[i];
            p.x = cscSegX[i];
            p.y = cscSegY[i];
            p.z = cscSegZ[i];
            p.t = cscSegT[i];
            // p.twire = cscSegT[i];
            p.chamber = cscSegChamber[i];
            p.station = cscSegStation[i];
            p.clusterID = UNCLASSIFIED;
            points.push_back(p);
         }
         //Do DBSCAN Clustering
         min_point_csc = min_clustering_segments; //minimum number of segments to call it a cluster
         epsilon_csc = 0.2; //cluster radius parameter
         DBSCAN ds_cscSeg(min_point_csc, epsilon_csc, points);
         ds_cscSeg.run();
         ds_cscSeg.result();
         ds_cscSeg.clusterMoments();
         ds_cscSeg.sort_clusters();
         ds_cscSeg.merge_clusters();
         ds_cscSeg.result();
         ds_cscSeg.clusterMoments();
         ds_cscSeg.sort_clusters();

         points.clear();

         for (auto & tmp: ds_cscSeg.CscCluster) {
            cout << "CSC SEGMENT clustering -> " << tmp.eta << "  " << tmp.phi << "  " << tmp.x << "  " << tmp.y << "  " << tmp.z << "  " << tmp.t << endl;
         }

         cout<<"======================> "<<ds_cscRechit.CscCluster.size()<<"   "<<ds_dtRechit.CscCluster.size()<<" ---------- "<<nCscRechitClusters<<"   "<<nDtRechitClusters<<endl;

         if (true) //twin higgs model && !hnl_model
         {
            int temp_id = 999;
            float temp_energy = 999.;
            int temp_status = 999;
            bool has_kaon = false;
            int nstau = 0;
            int nMother_stau = 0;
            int istau = 0;
            float mass = 0, P = 0, px = 0, py = 0, pz = 0;

            for (int i = 0; i < NTriggersMAX; i++) {
               if (HLTDecision[i] == 1) {
                  MuonSystem -> HLTDecision[i] = HLTDecision[i];
                  NTriggers -> Fill(i);
               }
            }

            //***********************************************************************************************************************************************************|
            //***********************************************************************************************************************************************************|
            //**************************************************************** NUMBER OF RECHITS ************************************************************************|
            //***********************************************************************************************************************************************************|
            //***********************************************************************************************************************************************************|
            for (int i = 0; i < nDtRechits; i++) {
               RECHITS  -> Fill(N_Good_Events, dtRechitEta[i], dtRechitPhi[i], dtRechitX[i], dtRechitY[i], dtRechitZ[i], dtRechitStation[i], 0);
            }
            for (int i = 0; i < ncscRechits; i++) {
               RECHITS  -> Fill(N_Good_Events, cscRechitsEta[i], cscRechitsPhi[i], cscRechitsX[i], cscRechitsY[i], cscRechitsZ[i], cscRechitsStation[i], 1);
            }
            //***********************************************************************************************************************************************************|
            //***********************************************************************************************************************************************************|
            //**************************************************************** NUMBER OF SEGMENTS ***********************************************************************|
            //***********************************************************************************************************************************************************|
            //***********************************************************************************************************************************************************|
            for (int i = 0; i < nDtSeg; i++) {
               SEGMENTS -> Fill(N_Good_Events, dtRechitEta[i], dtRechitPhi[i], dtRechitX[i], dtRechitY[i], dtRechitZ[i], dtRechitStation[i], 0);
            }
            for (int i = 0; i < nCscSeg; i++) {
               SEGMENTS -> Fill(N_Good_Events, cscRechitsEta[i], cscRechitsPhi[i], cscRechitsX[i], cscRechitsY[i], cscRechitsZ[i], cscRechitsStation[i], 1);
            }
            //***********************************************************************************************************************************************************|
            //***********************************************************************************************************************************************************|
            //**************************************************************** NUMBER OF RECHITS CL *********************************************************************|
            //***********************************************************************************************************************************************************|
            //***********************************************************************************************************************************************************|
            int i = 0;
            for (auto & tmp: ds_dtRechit.CscCluster) {
              cout<<"tmp.t "<<tmp.t<<endl;
              RECHITS_CL     -> Fill(N_Good_Events, tmp.eta, tmp.phi, tmp.x, tmp.y, tmp.z, tmp.t, tmp.nCscSegments, tmp.nStation10, tmp.avgStation10, 0);
              RECHITS_CL_ADD -> Fill(100000, 100000, 100000, 100000, 100000, 100000);
              i++;
            }
            i = 0;
            for (auto & tmp: ds_cscRechit.CscCluster) {
              cout<<"tmp.t "<<tmp.t<<endl;
              cout<<"VETOES: "<<cscRechitClusterJetVetoPt[i]<<"   "<<cscRechitClusterMuonVetoPt[i]<<endl;
              RECHITS_CL     -> Fill(N_Good_Events, tmp.eta, tmp.phi, tmp.x, tmp.y, tmp.z, tmp.t, tmp.nCscSegments, tmp.nStation10, tmp.avgStation10, 1);
              RECHITS_CL_ADD -> Fill(cscRechitClusterJetVetoPt[i], cscRechitClusterMuonVetoPt[i], cscRechitClusterNRechitChamberPlus11[i], cscRechitClusterNRechitChamberPlus12[i], cscRechitClusterNRechitChamberMinus11[i], cscRechitClusterNRechitChamberMinus12[i]);
              i++;
            }
            //***********************************************************************************************************************************************************|
            //***********************************************************************************************************************************************************|
            //**************************************************************** NUMBER OF SEGMENTS CL ********************************************************************|
            //***********************************************************************************************************************************************************|
            //***********************************************************************************************************************************************************|
            i = 0;
            for (auto & tmp: ds_dtSeg.CscCluster) {
              SEGMENTS_CL -> Fill(N_Good_Events, tmp.eta, tmp.phi, tmp.x, tmp.y, tmp.z, tmp.t, tmp.nCscSegments, tmp.nStation10, tmp.avgStation10, 0);
              i++;
            }
            i = 0;
            for (auto & tmp: ds_cscSeg.CscCluster) {
              SEGMENTS_CL -> Fill(N_Good_Events, tmp.eta, tmp.phi, tmp.x, tmp.y, tmp.z, tmp.t, tmp.nCscSegments, tmp.nStation10, tmp.avgStation10, 1);
              i++;
            }
            //***********************************************************************************************************************************************************|
            //***********************************************************************************************************************************************************|
            //***********************************************************************************************************************************************************|
            //***********************************************************************************************************************************************************|
            //***********************************************************************************************************************************************************|
            // EVENTS -> Fill(N_Good_Events, N_Rechits_dt, N_Segments_dt, N_Rechits_csc, N_Segments_csc, N_Rechits_CL_dt, N_Segments_CL_dt, N_Rechits_CL_csc, N_Segments_CL_csc);
            EVENTS -> Fill(N_Good_Events);

            N_Good_Events++;
            // }
            // ******************************************************************************** My Changes END

         } {
            for (int i = 0; i < 4; i++) {
               int llp_index = 1;
               if (i < 2) llp_index = 0;

               MuonSystem -> gLLP_visE[llp_index] += gLLP_daughter_e[i];
               MuonSystem -> gLLP_visEz[llp_index] += gLLP_daughter_e[i] / cosh(gLLP_daughter_eta[i]) * sinh(gLLP_daughter_eta[i]);
               MuonSystem -> gLLP_visP[llp_index] += gLLP_daughter_pt[i] * cosh(gLLP_daughter_eta[i]);
               MuonSystem -> gLLP_visPz[llp_index] += gLLP_daughter_pt[i] * sinh(gLLP_daughter_eta[i]);
               if (abs(gLLP_daughter_id[i]) == 11 || abs(gLLP_daughter_id[i]) == 22 || abs(gLLP_daughter_id[i]) == 111) //EM don't count muon energy abs(gParticleId[i])==13 ||
               {
                  MuonSystem -> gLLP_EMFracE[llp_index] += gLLP_daughter_e[i];
                  MuonSystem -> gLLP_EMFracEz[llp_index] += gLLP_daughter_e[i] / cosh(gLLP_daughter_eta[i]) * sinh(gLLP_daughter_eta[i]);
                  MuonSystem -> gLLP_EMFracP[llp_index] += gLLP_daughter_pt[i] * cosh(gLLP_daughter_eta[i]);
                  MuonSystem -> gLLP_EMFracPz[llp_index] += gLLP_daughter_pt[i] * sinh(gLLP_daughter_eta[i]);
               }
            }
         }
         /*
         for (int i = 0; i < 2; i++) {
            if (MuonSystem -> gLLP_visE[i] > 0) {
               MuonSystem -> gLLP_EMFracE[i] /= MuonSystem -> gLLP_visE[i];
               MuonSystem -> gLLP_EMFracEz[i] /= MuonSystem -> gLLP_visEz[i];
               MuonSystem -> gLLP_EMFracP[i] /= MuonSystem -> gLLP_visP[i];
               MuonSystem -> gLLP_EMFracPz[i] /= MuonSystem -> gLLP_visPz[i];
            }
         }
         // else MuonSystem->gLLP_daughterKaon[llp_index] = false;

      // MuonSystem->gLLP_daughter_deltaR[0] = deltaR(gLLP_eta[0],gLLP_phi[0],gLLP_eta[1],gLLP_phi[1]);
      // MuonSystem->gLLP_daughter_deltaR[1] = deltaR(gLLP_eta[2],gLLP_phi[2],gLLP_eta[3],gLLP_phi[3]);
      // cout<<"LLP"<<endl;

      for (int i = 0; i < nBunchXing; i++) {
         if (BunchXing[i] == 0) {
            MuonSystem -> npu = nPUmean[i];
         }
      }
      MuonSystem -> pileupWeight = helper -> getPileupWeight(MuonSystem -> npu);
      MuonSystem -> pileupWeightUp = helper -> getPileupWeightUp(MuonSystem -> npu) / MuonSystem -> pileupWeight;
      MuonSystem -> pileupWeightDown = helper -> getPileupWeightDown(MuonSystem -> npu) / MuonSystem -> pileupWeight;
      */
   }

   //get NPU
   /*
   MuonSystem -> npv = nPV;
   MuonSystem -> rho = fixedGridRhoFastjetAll;
   MuonSystem -> met = metType1Pt;
   MuonSystem -> metPhi = metType1Phi;
   MuonSystem -> metJESUp = MuonSystem -> met;
   MuonSystem -> metJESDown = MuonSystem -> met;
   MuonSystem -> metSF = helper -> getMetTriggerSF(MuonSystem -> met);

   if (signalScan && !isData) Total2D[make_pair(MuonSystem -> mX, MuonSystem -> ctau)] -> Fill(1.0, genWeight * MuonSystem -> higgsPtWeight * MuonSystem -> pileupWeight);

   std::pair < double, double > corrected_met;
   if (analysisTag == "Razor2016_07Aug2017Rereco") corrected_met = helper -> METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2016, !isData, nPV);
   else if (analysisTag == "Razor2017_17Nov2017Rereco") corrected_met = helper -> METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2017, !isData, nPV);
   else if (analysisTag == "Razor2018_17SeptEarlyReReco") corrected_met = helper -> METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2018, !isData, nPV);

   MuonSystem -> metXYCorr = corrected_met.first;
   MuonSystem -> metPhiXYCorr = corrected_met.second;

   if (!isData && !hnl_model) {
      if (MuonSystem -> gLLP_csc[0] == false && MuonSystem -> gLLP_csc[1] == false) continue;
   }
   // CHANGE FOR CSC ACCEPTANCE
   //if (MuonSystem -> gLLP_csc[0] == 0 && MuonSystem -> gLLP_csc[1] == 0) continue;
   //cout << MuonSystem -> gLLP_csc[0] << " ==== " << MuonSystem -> gLLP_csc[1] << endl;

   if (signalScan && !isData) accep2D[make_pair(MuonSystem -> mX, MuonSystem -> ctau)] -> Fill(1.0, genWeight * MuonSystem -> higgsPtWeight * MuonSystem -> pileupWeight);
   else if (!isData) accep -> Fill(1.0, genWeight * MuonSystem -> higgsPtWeight * MuonSystem -> pileupWeight);

   if (analysisTag == "Razor2016_07Aug2017Rereco" || analysisTag == "Razor2016_Source2018") {
      MuonSystem -> METTrigger = HLTDecision[310] || HLTDecision[467];
      MuonSystem -> METNoMuTrigger = HLTDecision[467];
   } else {
      MuonSystem -> METTrigger = HLTDecision[310] || HLTDecision[467] || HLTDecision[703] || HLTDecision[717] || HLTDecision[710] || HLTDecision[709];
      MuonSystem -> METNoMuTrigger = HLTDecision[467] || HLTDecision[717] || HLTDecision[710];
   }

   // if (MuonSystem -> met < 200 || !MuonSystem -> METNoMuTrigger) continue; // CHANGE FOR GLLP MATCHING

   if (signalScan && !isData) accep_met2D[make_pair(MuonSystem -> mX, MuonSystem -> ctau)] -> Fill(1.0, genWeight * MuonSystem -> higgsPtWeight * MuonSystem -> pileupWeight * MuonSystem -> metSF);
   else if (!isData) accep_met -> Fill(1.0, genWeight * MuonSystem -> higgsPtWeight * MuonSystem -> pileupWeight);
   else Nmet200 -> Fill(1.0);

   // flags
   MuonSystem -> Flag_HBHENoiseFilter = Flag_HBHENoiseFilter;
   MuonSystem -> Flag_HBHEIsoNoiseFilter = Flag_HBHEIsoNoiseFilter;
   MuonSystem -> Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
   // MuonSystem->Flag_CSCTightHaloFilter = Flag_CSCTightHaloFilter;
   MuonSystem -> Flag_globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
   MuonSystem -> Flag_goodVertices = Flag_goodVertices;
   MuonSystem -> Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;
   MuonSystem -> Flag_BadChargedCandidateFilter = Flag_BadChargedCandidateFilter;
   MuonSystem -> Flag_eeBadScFilter = Flag_eeBadScFilter;
   // MuonSystem->Flag_all = Flag_goodVertices && Flag_HBHEIsoNoiseFilter && Flag_BadPFMuonFilter && Flag_CSCTightHaloFilter && Flag_goodVertices && Flag_ecalBadCalibFilter;
   // cout<<"HBHE noise filter"<<endl;
   MuonSystem -> Flag2_HBHENoiseFilter = Flag2_HBHENoiseFilter;
   MuonSystem -> Flag2_HBHEIsoNoiseFilter = Flag2_HBHEIsoNoiseFilter;
   MuonSystem -> Flag2_BadPFMuonFilter = Flag2_BadPFMuonFilter;
   MuonSystem -> Flag2_globalSuperTightHalo2016Filter = Flag2_globalSuperTightHalo2016Filter;
   MuonSystem -> Flag2_globalTightHalo2016Filter = Flag2_globalTightHalo2016Filter;
   MuonSystem -> Flag2_BadChargedCandidateFilter = Flag2_BadChargedCandidateFilter;
   // Flag2_goodVertices = Flag2_goodVertices;
   MuonSystem -> Flag2_EcalDeadCellTriggerPrimitiveFilter = Flag2_EcalDeadCellTriggerPrimitiveFilter;
   MuonSystem -> Flag2_ecalBadCalibFilter = Flag2_ecalBadCalibFilter;
   MuonSystem -> Flag2_eeBadScFilter = Flag2_eeBadScFilter;
   MuonSystem -> Flag2_all = Flag2_HBHENoiseFilter && Flag2_HBHEIsoNoiseFilter && Flag2_BadPFMuonFilter && Flag2_globalSuperTightHalo2016Filter && Flag2_EcalDeadCellTriggerPrimitiveFilter;
   if (isData) MuonSystem -> Flag2_all = MuonSystem -> Flag2_all && Flag2_eeBadScFilter;

   if (analysisTag != "Razor2016_07Aug2017Rereco" && analysisTag != "Razor2016_Source2018") {
      MuonSystem -> Flag2_all = MuonSystem -> Flag2_all && Flag2_ecalBadCalibFilter;
   }

   if (isData) {
      if (!MuonSystem -> Flag2_all) continue;
      NmetFilter -> Fill(1.0);
   }

   //*************************************************************************
   //Start Object Selection
   //*************************************************************************
   // cout<<"reached muons"<<endl;

   std::vector < leptons > Leptons;

   //-------------------------------
   //Muons
   //-------------------------------
   for (int i = 0; i < nMuons; i++) {
      float muonIso = (muon_chargedIso[i] + fmax(0.0, muon_photonIso[i] + muon_neutralHadIso[i] - 0.5 * muon_pileupIso[i])) / muonPt[i];
      //      if (!isMuonPOGLooseMuon(i)) continue;
      //      if (muonPt[i] < 20) continue;
      //if (fabs(muonEta[i]) > 2.4) continue;
      //remove overlaps

      bool overlap = false;
      for (auto & lep: Leptons) {
         if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], lep.lepton.Eta(), lep.lepton.Phi()) < 0.3) overlap = true;
      }
      if (overlap) continue;

      leptons tmpMuon;
      tmpMuon.lepton.SetPtEtaPhiM(muonPt[i], muonEta[i], muonPhi[i], MU_MASS);
      tmpMuon.pdgId = 13 * -1 * muonCharge[i];
      tmpMuon.dZ = muon_dZ[i];
      tmpMuon.passId = isMuonPOGTightMuon(i, true, false);
      tmpMuon.passLooseIso = muonIso < 0.25;
      tmpMuon.passTightIso = muonIso < 0.15;
      tmpMuon.passVTightIso = muonIso < 0.10;
      tmpMuon.passVVTightIso = muonIso < 0.05;

      tmpMuon.passId = isLooseMuon(i);
      tmpMuon.passVetoId = false;
      //tmpMuon.passVetoPOGId = false;

      Leptons.push_back(tmpMuon);
   }
   // cout<<"jere"<<endl;

   //-------------------------------
   //Electrons
   //-------------------------------

   for (int i = 0; i < nElectrons; i++) {

      // if (!isEGammaPOGLooseElectron(i, true, true, true, "Summer16")) continue;
      // if (!isEGammaPOGLooseElectron(i, true, true, true, "2017_94X")) continue;

      if (!isEGammaPOGLooseElectron(i, true, false, true, "vid")) continue;

      if (elePt[i] < zh_lepton1_cut) continue;
      if (fabs(eleEta[i]) > 2.4) continue;

      //remove overlaps
      bool overlap = false;
      for (auto & lep: Leptons) {
         if (RazorAnalyzer::deltaR(eleEta[i], elePhi[i], lep.lepton.Eta(), lep.lepton.Phi()) < 0.3) overlap = true;
      }
      if (overlap) continue;
      leptons tmpElectron;
      tmpElectron.lepton.SetPtEtaPhiM(elePt[i], eleEta[i], elePhi[i], ELE_MASS);
      tmpElectron.pdgId = 11 * -1 * eleCharge[i];
      tmpElectron.dZ = ele_dZ[i];
      tmpElectron.passId = isEGammaPOGLooseElectron(i, true, true, true, "2017_94X");

      tmpElectron.passVetoId = isEGammaPOGLooseElectron(i, true, true, true, "vid");
      Leptons.push_back(tmpElectron);
   }

   sort(Leptons.begin(), Leptons.end(), my_largest_pt);

   //----------------
   //Find Z Candidate
   //----------------

   double ZMass = -999;
   double ZPt = -999;
   double tmpDistToZPole = 9999;
   pair < uint, uint > ZCandidateLeptonIndex;
   bool foundZ = false;
   TLorentzVector ZCandidate;
   double leadingLepPt = 0.0;
   for (uint i = 0; i < Leptons.size(); i++) {
      for (uint j = i + 1; j < Leptons.size(); j++) {
         if (!(Leptons[i].pdgId == -1 * Leptons[j].pdgId)) continue; // same flavor opposite charge
         double tmpMass = (Leptons[i].lepton + Leptons[j].lepton).M();
         //select the pair closest to Z pole mass
         if (fabs(tmpMass - Z_MASS) < tmpDistToZPole) {
            tmpDistToZPole = tmpMass;
            if (Leptons[i].pdgId > 0) {
               ZCandidateLeptonIndex = pair < int, int > (i, j);
            } else {
               ZCandidateLeptonIndex = pair < int, int > (j, i);
            }
            ZMass = tmpMass;
            ZPt = (Leptons[i].lepton + Leptons[j].lepton).Pt();
            ZCandidate = Leptons[i].lepton + Leptons[j].lepton;
            leadingLepPt = max(Leptons[i].lepton.Pt(), Leptons[j].lepton.Pt());
            foundZ = true;
         }
      }
   }

   // if (foundZ  && Leptons.size() == 2 && leadingLepPt > zh_lepton0_cut)
   if (foundZ && Leptons.size() == 2) {
      MuonSystem -> ZMass = ZMass;
      MuonSystem -> ZPt = ZPt;
      MuonSystem -> ZEta = ZCandidate.Eta();
      MuonSystem -> ZPhi = ZCandidate.Phi();
      MuonSystem -> ZleptonIndex1 = ZCandidateLeptonIndex.first;
      MuonSystem -> ZleptonIndex2 = ZCandidateLeptonIndex.second;
      MuonSystem -> category = 2;
   } // endif foundZ
   else {
      for (unsigned int i = Leptons.size(); i > 0; --i) {
         int index = i - 1;
         if (abs(Leptons[index].pdgId) == 13 && Leptons[index].lepton.Pt() < wh_muonPt_cut) {
            Leptons.erase(Leptons.begin() + index);
         } else if (abs(Leptons[index].pdgId) == 11 && Leptons[index].lepton.Pt() < wh_elePt_cut) {
            Leptons.erase(Leptons.begin() + index);
         }
      }
      if (Leptons.size() == 1) MuonSystem -> category = 1;
      else MuonSystem -> category = 0;
   }
   bool tag = false;
   */
   /*
   for (auto & tmp: Leptons) {
     MuonSystem -> lepE[MuonSystem -> nLeptons] = tmp.lepton.E();
     MuonSystem -> lepPt[MuonSystem -> nLeptons] = tmp.lepton.Pt();
     MuonSystem -> lepEta[MuonSystem -> nLeptons] = tmp.lepton.Eta();
     MuonSystem -> lepPhi[MuonSystem -> nLeptons] = tmp.lepton.Phi();
     MuonSystem -> lepPdgId[MuonSystem -> nLeptons] = tmp.pdgId;
     MuonSystem -> lepDZ[MuonSystem -> nLeptons] = tmp.dZ;
     MuonSystem -> lepPassId[MuonSystem -> nLeptons] = tmp.passId;
     MuonSystem -> lepPassVetoId[MuonSystem -> nLeptons] = tmp.passVetoId;

     MuonSystem -> lepPassLooseIso[MuonSystem -> nLeptons] = tmp.passLooseIso;
     MuonSystem -> lepPassTightIso[MuonSystem -> nLeptons] = tmp.passTightIso;
     MuonSystem -> lepPassVTightIso[MuonSystem -> nLeptons] = tmp.passVTightIso;
     MuonSystem -> lepPassVVTightIso[MuonSystem -> nLeptons] = tmp.passVVTightIso;
     if (!isData) {
       MuonSystem -> lepTriggerSF[MuonSystem -> nLeptons] = helper -> getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), true, true, true);
       MuonSystem -> lepTightIdSF[MuonSystem -> nLeptons] = helper -> getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, false);
       MuonSystem -> lepLooseIdSF[MuonSystem -> nLeptons] = helper -> getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, true);
       MuonSystem -> lepTightIsoSF[MuonSystem -> nLeptons] = helper -> getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, false);
       MuonSystem -> lepLooseIsoSF[MuonSystem -> nLeptons] = helper -> getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, true);
       MuonSystem -> lepTriggerMCEfficiency[MuonSystem -> nLeptons] = helper -> getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), true, true, true);
       MuonSystem -> lepTightIdMCEfficiency[MuonSystem -> nLeptons] = helper -> getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, false);
       MuonSystem -> lepLooseIdMCEfficiency[MuonSystem -> nLeptons] = helper -> getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, true);
       MuonSystem -> lepTightIsoMCEfficiency[MuonSystem -> nLeptons] = helper -> getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, false);
       MuonSystem -> lepLooseIsoMCEfficiency[MuonSystem -> nLeptons] = helper -> getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, true);
     }

     if (MuonSystem -> lepPassTightIso[MuonSystem -> nLeptons] && MuonSystem -> lepPassId[MuonSystem -> nLeptons] && !tag) {
       MuonSystem -> lepTag[MuonSystem -> nLeptons] = true;
       if (!isData) {
         MuonSystem -> lepSF[MuonSystem -> nLeptons] = MuonSystem -> lepTriggerSF[MuonSystem -> nLeptons] * MuonSystem -> lepTightIdSF[MuonSystem -> nLeptons] * MuonSystem -> lepTightIsoSF[MuonSystem -> nLeptons];
         MuonSystem -> lepEff[MuonSystem -> nLeptons] = MuonSystem -> lepTriggerMCEfficiency[MuonSystem -> nLeptons] * MuonSystem -> lepTightIdMCEfficiency[MuonSystem -> nLeptons] * MuonSystem -> lepTightIsoMCEfficiency[MuonSystem -> nLeptons];
       } else {
         MuonSystem -> lepSF[MuonSystem -> nLeptons] = 1.0;
         MuonSystem -> lepEff[MuonSystem -> nLeptons] = 1.0;
       }
       tag = true;
     } else {
       MuonSystem -> lepTag[MuonSystem -> nLeptons] = false;
       if (!isData) {
         MuonSystem -> lepSF[MuonSystem -> nLeptons] = MuonSystem -> lepTriggerSF[MuonSystem -> nLeptons] * MuonSystem -> lepLooseIdSF[MuonSystem -> nLeptons] * MuonSystem -> lepLooseIsoSF[MuonSystem -> nLeptons];
         MuonSystem -> lepEff[MuonSystem -> nLeptons] = MuonSystem -> lepTriggerMCEfficiency[MuonSystem -> nLeptons] * MuonSystem -> lepLooseIdMCEfficiency[MuonSystem -> nLeptons] * MuonSystem -> lepLooseIsoMCEfficiency[MuonSystem -> nLeptons];
       } else {
         MuonSystem -> lepSF[MuonSystem -> nLeptons] = 1.0;
         MuonSystem -> lepEff[MuonSystem -> nLeptons] = 1.0;
       }
     }

     if (!isData) {
       for (int i = 0; i < nGenParticle; i++) {
         if ((abs(gParticleId[i]) == 13 || abs(gParticleId[i]) == 11) && abs(gParticleMotherId[i]) == 23) MuonSystem -> lepFromZ[MuonSystem -> nLeptons] = true;
       }
     }
     MuonSystem -> nLeptons++;
   }
   */
   
   /*
   TLorentzVector met;
   met.SetPtEtaPhiE(metType1Pt, 0, metType1Phi, metType1Pt);
   if (Leptons.size() > 0) {
      TLorentzVector visible = Leptons[0].lepton;
      MuonSystem -> MT = GetMT(visible, met);
   }

   //cout<<"njets:"<<nJets<<endl;
   //-----------------------------------------------
   //Select Jets
   //-----------------------------------------------
   std::vector < jets > Jets;
   float MetXCorr_JESUp = 0;
   float MetYCorr_JESUp = 0;
   float MetXCorr_JESDown = 0;
   float MetYCorr_JESDown = 0;
   float MetXCorr_EENoise = 0;
   float MetYCorr_EENoise = 0;
   float MetXCorr_HEM = 0;
   float MetYCorr_HEM = 0;
   for (int i = 0; i < nJets; i++) {
      //------------------------------------------------------------
      //exclude selected muons and electrons from the jet collection
      //------------------------------------------------------------

      double deltaR = -1;
      for (auto & lep: Leptons) {
         double thisDR = RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], lep.lepton.Eta(), lep.lepton.Phi());
         if (deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      }
      if (deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

      //------------------------------------------------------------
      //Apply Jet Energy and Resolution Corrections
      //------------------------------------------------------------
      //double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
      //fixedGridRhoFastjetAll, jetJetArea[i] , JetCorrector);
      // cout<<"before JEC"<<endl;
      //double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i],
      //fixedGridRhoFastjetAll, jetJetArea[i], runNum, JetCorrectorIOV, JetCorrector);
      double JEC = 1.0;
      double jetCorrPt = jetPt[i] * JEC;
      double jetCorrE = jetE[i] * JEC;
      // cout<<JEC<<endl;
      // cout<<"corrected pt: "<<jetCorrPt<<", "<<"eta: "<<jetEta[i]<<endl;
      TLorentzVector thisJet = makeTLorentzVector(jetCorrPt, jetEta[i], jetPhi[i], jetCorrE);

      if (thisJet.Eta() > -3.0 && thisJet.Eta() < -1.3 && thisJet.Phi() > -1.57 && thisJet.Phi() < -0.87 && analysisTag == "Razor2018_17SeptEarlyReReco") {
         MetXCorr_HEM += thisJet.Px();
         MetYCorr_HEM += thisJet.Py();
      }
      if (fabs(thisJet.Eta()) > 2.65 && fabs(thisJet.Eta()) < 3.139 && thisJet.Pt() < 50 && analysisTag == "Razor2017_17Nov2017Rereco") {
         MetXCorr_EENoise += thisJet.Px();
         MetYCorr_EENoise += thisJet.Py();
         //if (eventNum==123969624)
         //{
         // cout<<i<<", "<<thisJet.Eta()<<","<<jetEta[i]<<endl;
         // cout<<MetXCorr_EENoise<<", "<<MetYCorr_EENoise<<", "<<thisJet.Px()<<", "<<thisJet.Py()<<","<<thisJet.Pt()<<","<<jetPt[i]<<endl;
         //}
      }
      if (fabs(thisJet.Eta()) > 2.25 && fabs(thisJet.Eta()) < 3.0 && thisJet.Pt() > 100 && (analysisTag == "Razor2016_07Aug2017Rereco" || analysisTag == "Razor2017_17Nov2017Rereco")) {
         MuonSystem -> EE_prefiring = false;
      }

      // if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
      if (fabs(thisJet.Eta()) >= 3.0) continue;
      jets tmpJet;
      tmpJet.jet = thisJet;
      tmpJet.passId = jetPassIDTight[i];
      // calculate jet energy scale uncertainty

      //double unc = helper -> getJecUnc(jetCorrPt, jetEta[i], runNum); //use run=999 as default
      //tmpJet.jetPtJESUp   = jetCorrPt * (1 + unc);
      //tmpJet.jetPtJESDown = jetCorrPt * (1 - unc);
      //tmpJet.jetEJESUp    = jetCorrE * (1 + unc);
      //tmpJet.jetEJESDown  = jetCorrE * (1 - unc);
      //tmpJet.JecUnc = unc;

      TLorentzVector thisJetJESUp = makeTLorentzVector(tmpJet.jetPtJESUp, jetEta[i], jetPhi[i], tmpJet.jetEJESUp);
      TLorentzVector thisJetJESDown = makeTLorentzVector(tmpJet.jetPtJESDown, jetEta[i], jetPhi[i], tmpJet.jetEJESDown);
      if (tmpJet.jetPtJESUp > 10) {
         MetXCorr_JESUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
         MetYCorr_JESUp += -1 * (thisJetJESUp.Py() - thisJet.Py());
      }
      if (tmpJet.jetPtJESDown > 10) {
         MetXCorr_JESDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
         MetYCorr_JESDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
      }
      if (thisJet.Pt() < 20) continue; //According to the April 1st 2015 AN
      //if (!jetPassIDLoose[i]) continue;
      Jets.push_back(tmpJet);
   }

   // cout<<"here"<<nJets<<endl;

   sort(Jets.begin(), Jets.end(), my_largest_pt_jet);
   // cout<<"here"<<nJets<<endl;
   if (Jets.size() > 0) {
      MuonSystem -> jetMet_dPhi = RazorAnalyzer::deltaPhi(jetPhi[0], metType1Phi);
   } else {
      MuonSystem -> jetMet_dPhi = -999.;
   }
   double jetMet_dPhiMin_temp = 999.;
   double jetMet_dPhiMin4_temp = 999.;
   */
   /*
   for (auto & tmp: Jets) {
     if (tmp.jet.Pt() < 50) continue;
     if (abs(tmp.jet.Eta()) > 2.4) continue;
     MuonSystem -> jetE[MuonSystem -> nJets] = tmp.jet.E();
     MuonSystem -> jetPt[MuonSystem -> nJets] = tmp.jet.Pt();
     MuonSystem -> jetEta[MuonSystem -> nJets] = tmp.jet.Eta();
     MuonSystem -> jetPhi[MuonSystem -> nJets] = tmp.jet.Phi();
     MuonSystem -> jetTime[MuonSystem -> nJets] = tmp.time;

     MuonSystem -> jetPtJESUp[MuonSystem -> nJets] = tmp.jetPtJESUp;
     MuonSystem -> jetPtJESDown[MuonSystem -> nJets] = tmp.jetPtJESDown;
     MuonSystem -> jetEJESUp[MuonSystem -> nJets] = tmp.jetEJESUp;
     MuonSystem -> jetEJESDown[MuonSystem -> nJets] = tmp.jetEJESDown;
     MuonSystem -> JecUnc[MuonSystem -> nJets] = tmp.JecUnc;

     if (jetMet_dPhiMin4_temp > abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(), metType1Phi)) && MuonSystem -> nJets < 4) {
       jetMet_dPhiMin4_temp = abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(), metType1Phi));

     }
     if (jetMet_dPhiMin_temp > abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(), metType1Phi))) {
       if (tmp.jet.Pt() > 30 && abs(tmp.jet.Eta()) < 2.4) {
         jetMet_dPhiMin_temp = abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(), metType1Phi));

       }
     }
     MuonSystem -> jetTightPassId[MuonSystem -> nJets] = tmp.passId;
     MuonSystem -> nJets++;
   }
   */
   // if (isData && MuonSystem->nJets==0)continue;
   // cout<<"here"<<MuonSystem->nJets<<endl;
   /*
   Njet1 -> Fill(1.0);
   MuonSystem -> jetMet_dPhiMin = jetMet_dPhiMin_temp;
   MuonSystem -> jetMet_dPhiMin4 = jetMet_dPhiMin4_temp;
   TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);

   //JES up
   float PFMetXJESUp = PFMET.Px() + MetXCorr_JESUp;
   float PFMetYJESUp = PFMET.Py() + MetYCorr_JESUp;
   MuonSystem -> metJESUp = sqrt(pow(PFMetXJESUp, 2) + pow(PFMetYJESUp, 2));
   MuonSystem -> metJESUpSF = helper -> getMetTriggerSF(MuonSystem -> metJESUp) / MuonSystem -> metSF;
   cout << PFMetXJESUp << "    " << PFMetYJESUp << endl;

   //JES down
   float PFMetXJESDown = PFMET.Px() + MetXCorr_JESDown;
   float PFMetYJESDown = PFMET.Py() + MetYCorr_JESDown;
   MuonSystem -> metJESDown = sqrt(pow(PFMetXJESDown, 2) + pow(PFMetYJESDown, 2));
   MuonSystem -> metJESDownSF = helper -> getMetTriggerSF(MuonSystem -> metJESDown) / MuonSystem -> metSF;
   cout << PFMetXJESDown << "    " << PFMetYJESDown << endl;

   //EENoise
   float PFMetXEENoise = PFMET.Px() + MetXCorr_EENoise;
   float PFMetYEENoise = PFMET.Py() + MetYCorr_EENoise;
   MuonSystem -> metEENoise = sqrt(pow(PFMetXEENoise, 2) + pow(PFMetYEENoise, 2));
   MuonSystem -> metPhiEENoise = atan(PFMetYEENoise / PFMetXEENoise);
   cout << PFMetXEENoise << "    " << PFMetYEENoise << endl;

   cout << "met" << endl;

   // count dt rechits
   MuonSystem -> nDTRechits = 0;
   for (int i = 0; i < nDtRechits; i++) {

      if (dtRechitY[i] >= 0.0) MuonSystem -> nDTPositiveYRechits++;
      else MuonSystem -> nDTNegativeYRechits++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -2) MuonSystem -> nDTRechitsChamberMinus12++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -1) MuonSystem -> nDTRechitsChamberMinus11++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 0) MuonSystem -> nDTRechitsChamber10++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 1) MuonSystem -> nDTRechitsChamberPlus11++;
      if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 2) MuonSystem -> nDTRechitsChamberPlus12++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -2) MuonSystem -> nDTRechitsChamberMinus22++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -1) MuonSystem -> nDTRechitsChamberMinus21++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 0) MuonSystem -> nDTRechitsChamber20++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 1) MuonSystem -> nDTRechitsChamberPlus21++;
      if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 2) MuonSystem -> nDTRechitsChamberPlus22++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -2) MuonSystem -> nDTRechitsChamberMinus32++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -1) MuonSystem -> nDTRechitsChamberMinus31++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 0) MuonSystem -> nDTRechitsChamber30++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 1) MuonSystem -> nDTRechitsChamberPlus31++;
      if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 2) MuonSystem -> nDTRechitsChamberPlus32++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -2) MuonSystem -> nDTRechitsChamberMinus42++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -1) MuonSystem -> nDTRechitsChamberMinus41++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 0) MuonSystem -> nDTRechitsChamber40++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 1) MuonSystem -> nDTRechitsChamberPlus41++;
      if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 2) MuonSystem -> nDTRechitsChamberPlus42++;

      MuonSystem -> nDTRechits++;
   }
   if (MuonSystem -> nDTRechitsChamberMinus12 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberMinus11 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamber10 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberPlus11 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberPlus12 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberMinus22 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberMinus21 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamber20 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberPlus21 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberPlus22 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberMinus32 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberMinus31 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamber30 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberPlus31 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberPlus32 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberMinus42 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberMinus41 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamber40 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberPlus41 > 50) MuonSystem -> nDtRings++;
   if (MuonSystem -> nDTRechitsChamberPlus42 > 50) MuonSystem -> nDtRings++;

   // cout<<"jets"<<endl;
   cout << "Number of rec hits: " << ncscRechits << endl;

   vector < Point > points;
   vector < int > cscRechitsClusterId;
   points.clear();
   MuonSystem -> nCscRechits = 0;

   for (int i = 0; i < ncscRechits; i++) {
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
      if (cscRechitsY[i] >= 0.0) {
         MuonSystem -> nCscPositiveYRechits++;
         MuonSystem -> cscPosTpeak = MuonSystem -> cscPosTpeak + cscRechitsTpeak[i];
      } else {
         MuonSystem -> nCscNegativeYRechits++;
         MuonSystem -> cscNegTpeak = MuonSystem -> cscNegTpeak + cscRechitsTpeak[i];
      }
      if (cscRechitsTpeak[i] < -12.5) MuonSystem -> nEarlyCscRechits++;
      if (cscRechitsTpeak[i] > 12.5) MuonSystem -> nLateCscRechits++;
      if (cscRechitsTpeak[i] < -25) MuonSystem -> nEarly2CscRechits++;
      if (cscRechitsTpeak[i] > 25) MuonSystem -> nLate2CscRechits++;
      if (cscRechitsChamber[i] == 11) MuonSystem -> nCscRechitsChamberPlus11++;
      if (cscRechitsChamber[i] == 12) MuonSystem -> nCscRechitsChamberPlus12++;
      if (cscRechitsChamber[i] == 13) MuonSystem -> nCscRechitsChamberPlus13++;
      if (cscRechitsChamber[i] == 21) MuonSystem -> nCscRechitsChamberPlus21++;
      if (cscRechitsChamber[i] == 22) MuonSystem -> nCscRechitsChamberPlus22++;
      if (cscRechitsChamber[i] == 31) MuonSystem -> nCscRechitsChamberPlus31++;
      if (cscRechitsChamber[i] == 32) MuonSystem -> nCscRechitsChamberPlus32++;
      if (cscRechitsChamber[i] == 41) MuonSystem -> nCscRechitsChamberPlus41++;
      if (cscRechitsChamber[i] == 42) MuonSystem -> nCscRechitsChamberPlus42++;
      if (cscRechitsChamber[i] == -11) MuonSystem -> nCscRechitsChamberMinus11++;
      if (cscRechitsChamber[i] == -12) MuonSystem -> nCscRechitsChamberMinus12++;
      if (cscRechitsChamber[i] == -13) MuonSystem -> nCscRechitsChamberMinus13++;
      if (cscRechitsChamber[i] == -21) MuonSystem -> nCscRechitsChamberMinus21++;
      if (cscRechitsChamber[i] == -22) MuonSystem -> nCscRechitsChamberMinus22++;
      if (cscRechitsChamber[i] == -31) MuonSystem -> nCscRechitsChamberMinus31++;
      if (cscRechitsChamber[i] == -32) MuonSystem -> nCscRechitsChamberMinus32++;
      if (cscRechitsChamber[i] == -41) MuonSystem -> nCscRechitsChamberMinus41++;
      if (cscRechitsChamber[i] == -42) MuonSystem -> nCscRechitsChamberMinus42++;
      // MuonSystem->nCscRechits++;
   }
   MuonSystem -> nCscRings = 0;
   if (MuonSystem -> nCscRechitsChamberPlus11 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberPlus12 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberPlus13 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberPlus21 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberPlus22 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberPlus31 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberPlus32 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberPlus41 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberPlus42 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberMinus11 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberMinus12 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberMinus13 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberMinus21 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberMinus22 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberMinus31 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberMinus32 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberMinus41 > 50) MuonSystem -> nCscRings++;
   if (MuonSystem -> nCscRechitsChamberMinus42 > 50) MuonSystem -> nCscRings++;
   // if(isData && (MuonSystem->nDtRings+MuonSystem->nCscRings)>=10)continue;
   // NcosmicVeto -> Fill(1.0);

   //Do DBSCAN Clustering
   int min_point = min_clustering_rechits; //minimum number of segments to call it a cluster
   float epsilon = 0.2; //cluster radius parameter
   DBSCAN ds(min_point, epsilon, points);
   ds.run();
   ds.result();
   ds.clusterMoments();
   ds.sort_clusters();
   */

   //Save cluster information
   /**/
   // cluster merging

   //cout << "done clustering before: " << ds.CscCluster.size() << " , " << nCscRechitClusters << endl;

   //ds.merge_clusters();
   //ds.result();
   //ds.clusterMoments();
   //ds.sort_clusters();

   //cout << "done clustering after: " << ds.CscCluster.size() << " , " << nCscRechitClusters << endl;
   
   /*
   ds_csc_cluster -> Fill(ds.CscCluster.size());
   csc_cluster -> Fill(nCscRechitClusters);

   MuonSystem -> nCscRechitClusters3 = 0;
   for (auto & tmp: ds.CscCluster) {

      // if (tmp.tTotal>-12.5)continue;
      // if( tmp.nStation10 <2)continue;
      MuonSystem -> cscRechitCluster3X[MuonSystem -> nCscRechitClusters3] = tmp.x;
      MuonSystem -> cscRechitCluster3Y[MuonSystem -> nCscRechitClusters3] = tmp.y;
      MuonSystem -> cscRechitCluster3Z[MuonSystem -> nCscRechitClusters3] = tmp.z;
      MuonSystem -> cscRechitCluster3Time[MuonSystem -> nCscRechitClusters3] = tmp.t;
      MuonSystem -> cscRechitCluster3TimeTotal[MuonSystem -> nCscRechitClusters3] = tmp.tTotal;
      MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3] = tmp.eta;
      MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3] = tmp.phi;
      MuonSystem -> cscRechitCluster3MajorAxis[MuonSystem -> nCscRechitClusters3] = tmp.MajorAxis;
      MuonSystem -> cscRechitCluster3MinorAxis[MuonSystem -> nCscRechitClusters3] = tmp.MinorAxis;
      MuonSystem -> cscRechitCluster3XSpread[MuonSystem -> nCscRechitClusters3] = tmp.XSpread;
      MuonSystem -> cscRechitCluster3YSpread[MuonSystem -> nCscRechitClusters3] = tmp.YSpread;
      MuonSystem -> cscRechitCluster3ZSpread[MuonSystem -> nCscRechitClusters3] = tmp.ZSpread;
      MuonSystem -> cscRechitCluster3EtaPhiSpread[MuonSystem -> nCscRechitClusters3] = tmp.EtaPhiSpread;
      MuonSystem -> cscRechitCluster3XYSpread[MuonSystem -> nCscRechitClusters3] = tmp.XYSpread;
      MuonSystem -> cscRechitCluster3RSpread[MuonSystem -> nCscRechitClusters3] = tmp.RSpread;
      MuonSystem -> cscRechitCluster3DeltaRSpread[MuonSystem -> nCscRechitClusters3] = tmp.DeltaRSpread;

      MuonSystem -> cscRechitCluster3EtaSpread[MuonSystem -> nCscRechitClusters3] = tmp.EtaSpread;
      MuonSystem -> cscRechitCluster3PhiSpread[MuonSystem -> nCscRechitClusters3] = tmp.PhiSpread;
      MuonSystem -> cscRechitCluster3TimeSpread[MuonSystem -> nCscRechitClusters3] = tmp.TSpread;
      MuonSystem -> cscRechitCluster3Size[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegments;

      MuonSystem -> cscRechitCluster3NRechitChamberPlus11[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberPlus11;
      MuonSystem -> cscRechitCluster3NRechitChamberPlus12[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberPlus12;
      MuonSystem -> cscRechitCluster3NRechitChamberPlus13[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberPlus13;
      MuonSystem -> cscRechitCluster3NRechitChamberPlus21[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberPlus21;
      MuonSystem -> cscRechitCluster3NRechitChamberPlus22[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberPlus22;
      MuonSystem -> cscRechitCluster3NRechitChamberPlus31[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberPlus31;
      MuonSystem -> cscRechitCluster3NRechitChamberPlus32[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberPlus32;
      MuonSystem -> cscRechitCluster3NRechitChamberPlus41[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberPlus41;
      MuonSystem -> cscRechitCluster3NRechitChamberPlus42[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberPlus42;
      MuonSystem -> cscRechitCluster3NRechitChamberMinus11[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberMinus11;
      MuonSystem -> cscRechitCluster3NRechitChamberMinus12[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberMinus12;
      MuonSystem -> cscRechitCluster3NRechitChamberMinus13[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberMinus13;
      MuonSystem -> cscRechitCluster3NRechitChamberMinus21[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberMinus21;
      MuonSystem -> cscRechitCluster3NRechitChamberMinus22[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberMinus22;
      MuonSystem -> cscRechitCluster3NRechitChamberMinus31[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberMinus31;
      MuonSystem -> cscRechitCluster3NRechitChamberMinus32[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberMinus32;
      MuonSystem -> cscRechitCluster3NRechitChamberMinus41[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberMinus41;
      MuonSystem -> cscRechitCluster3NRechitChamberMinus42[MuonSystem -> nCscRechitClusters3] = tmp.nCscSegmentChamberMinus42;
      MuonSystem -> cscRechitCluster3MaxChamber[MuonSystem -> nCscRechitClusters3] = tmp.maxChamber;
      MuonSystem -> cscRechitCluster3MaxChamberRatio[MuonSystem -> nCscRechitClusters3] = 1.0 * tmp.maxChamberSegment / tmp.nCscSegments;
      MuonSystem -> cscRechitCluster3NChamber[MuonSystem -> nCscRechitClusters3] = tmp.nChamber;
      MuonSystem -> cscRechitCluster3MaxStation[MuonSystem -> nCscRechitClusters3] = tmp.maxStation;
      MuonSystem -> cscRechitCluster3MaxStationRatio[MuonSystem -> nCscRechitClusters3] = 1.0 * tmp.maxStationSegment / tmp.nCscSegments;
      MuonSystem -> cscRechitCluster3NStation[MuonSystem -> nCscRechitClusters3] = tmp.nStation;
      MuonSystem -> cscRechitCluster3NStation5[MuonSystem -> nCscRechitClusters3] = tmp.nStation5;
      MuonSystem -> cscRechitCluster3NStation10[MuonSystem -> nCscRechitClusters3] = tmp.nStation10;
      MuonSystem -> cscRechitCluster3NStation10perc[MuonSystem -> nCscRechitClusters3] = tmp.nStation10perc;
      MuonSystem -> cscRechitCluster3AvgStation[MuonSystem -> nCscRechitClusters3] = tmp.avgStation;
      MuonSystem -> cscRechitCluster3AvgStation5[MuonSystem -> nCscRechitClusters3] = tmp.avgStation5;
      MuonSystem -> cscRechitCluster3AvgStation10[MuonSystem -> nCscRechitClusters3] = tmp.avgStation10;
      MuonSystem -> cscRechitCluster3AvgStation10perc[MuonSystem -> nCscRechitClusters3] = tmp.avgStation10perc;
      MuonSystem -> cscRechitCluster3Me11Ratio[MuonSystem -> nCscRechitClusters3] = tmp.Me11Ratio;
      MuonSystem -> cscRechitCluster3Me12Ratio[MuonSystem -> nCscRechitClusters3] = tmp.Me12Ratio;
      //Jet veto/ muon veto
      MuonSystem -> cscRechitCluster3JetVetoPt[MuonSystem -> nCscRechitClusters3] = 0.0;
      MuonSystem -> cscRechitCluster3JetVetoE[MuonSystem -> nCscRechitClusters3] = 0.0;
      MuonSystem -> cscRechitCluster3MuonVetoPt[MuonSystem -> nCscRechitClusters3] = 0.0;
      MuonSystem -> cscRechitCluster3MuonVetoE[MuonSystem -> nCscRechitClusters3] = 0.0;
      MuonSystem -> cscRechitCluster3GenMuonVetoPt[MuonSystem -> nCscRechitClusters3] = 0.0;
      MuonSystem -> cscRechitCluster3GenMuonVetoE[MuonSystem -> nCscRechitClusters3] = 0.0;
      MuonSystem -> cscRechitCluster3IsoMuonVetoPt[MuonSystem -> nCscRechitClusters3] = 0.0;

      cout << "CSC cluster data for gllp: " << MuonSystem -> nCscRechitClusters3;
      printf("%20.3f %10.3f\n", MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]);

      // jet veto
      for (int i = 0; i < nJets; i++) {
         if (fabs(jetEta[i] > 3.0)) continue;
         if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.4 && jetPt[i] > MuonSystem -> cscRechitCluster3JetVetoPt[MuonSystem -> nCscRechitClusters3]) {
            MuonSystem -> cscRechitCluster3JetVetoPt[MuonSystem -> nCscRechitClusters3] = jetPt[i];
         }
         if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.4 && jetE[i] > MuonSystem -> cscRechitCluster3JetVetoE[MuonSystem -> nCscRechitClusters3]) {
            MuonSystem -> cscRechitCluster3JetVetoE[MuonSystem -> nCscRechitClusters3] = jetE[i];
         }
      }
      // genjet veto
      for (int i = 0; i < nGenJets; i++) {
         if (fabs(genJetEta[i] > 3.0)) continue;
         if (RazorAnalyzer::deltaR(genJetEta[i], genJetPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.4 && genJetPt[i] > MuonSystem -> cscRechitCluster3GenJetVetoPt[MuonSystem -> nCscRechitClusters3]) {
            MuonSystem -> cscRechitCluster3GenJetVetoPt[MuonSystem -> nCscRechitClusters3] = genJetPt[i];
         }
         if (RazorAnalyzer::deltaR(genJetEta[i], genJetPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.4 && genJetE[i] > MuonSystem -> cscRechitCluster3GenJetVetoE[MuonSystem -> nCscRechitClusters3]) {
            MuonSystem -> cscRechitCluster3GenJetVetoE[MuonSystem -> nCscRechitClusters3] = genJetE[i];
         }
      }
      float min_deltaR = 15.;
      int index = 999;

      for (int i = 0; i < nMuons; i++) {
         if (fabs(muonEta[i] > 3.0)) continue;
         float muonIso = (muon_chargedIso[i] + fmax(0.0, muon_photonIso[i] + muon_neutralHadIso[i] - 0.5 * muon_pileupIso[i])) / muonPt[i];
         if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.4 && muonPt[i] > MuonSystem -> cscRechitCluster3MuonVetoPt[MuonSystem -> nCscRechitClusters3]) {
            MuonSystem -> cscRechitCluster3MuonVetoPt[MuonSystem -> nCscRechitClusters3] = muonPt[i];
            MuonSystem -> cscRechitCluster3MuonVetoE[MuonSystem -> nCscRechitClusters3] = muonE[i];
            MuonSystem -> cscRechitCluster3MuonVetoPhi[MuonSystem -> nCscRechitClusters3] = muonPhi[i];
            MuonSystem -> cscRechitCluster3MuonVetoEta[MuonSystem -> nCscRechitClusters3] = muonEta[i];
            MuonSystem -> cscRechitCluster3MuonVetoLooseIso[MuonSystem -> nCscRechitClusters3] = muonIso < 0.25;
            MuonSystem -> cscRechitCluster3MuonVetoTightIso[MuonSystem -> nCscRechitClusters3] = muonIso < 0.15;
            MuonSystem -> cscRechitCluster3MuonVetoVTightIso[MuonSystem -> nCscRechitClusters3] = muonIso < 0.10;
            MuonSystem -> cscRechitCluster3MuonVetoVVTightIso[MuonSystem -> nCscRechitClusters3] = muonIso < 0.05;
            MuonSystem -> cscRechitCluster3MuonVetoTightId[MuonSystem -> nCscRechitClusters3] = isMuonPOGTightMuon(i);

         }
         //check if muon is isolated

      }

      // match to gen-level muon
      //if (!isData) {
      //if (true) {
      for (int i = 0; i < nGenParticle; i++) {
         if (abs(gParticleId[i]) != 13) continue;
         if (abs(gParticleMotherId[i]) > 24 || abs(gParticleMotherId[i]) < 23) continue;
         if (abs(gParticleStatus[i]) != 1) continue;
         // if (fabs(muonEta[i]>3.0)) continue;
         if (RazorAnalyzer::deltaR(gParticleEta[i], gParticlePhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.4) {
            MuonSystem -> cscRechitCluster3GenMuonVetoPt[MuonSystem -> nCscRechitClusters3] = gParticlePt[i];
            MuonSystem -> cscRechitCluster3GenMuonVetoE[MuonSystem -> nCscRechitClusters3] = gParticleE[i];

         }

      }

      // match to gen level LLP
      min_deltaR = 15.;
      index = 999;
      for (int j = 0; j < 2; j++) {
         printf("%10.3f %10.3f %10.3f %20.3f %10.3f %10.3f\n", gLLP_pt[j], gLLP_eta[j], gLLP_phi[j], gLLP_decay_vertex_x[j], gLLP_decay_vertex_y[j], gLLP_decay_vertex_z[j]);
         printf("%20.3f %10.3f\n", MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]);

         double current_delta_r = RazorAnalyzer::deltaR(MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3], gLLP_eta[j], gLLP_phi[j]);
         if (current_delta_r < min_deltaR) {
            min_deltaR = current_delta_r;
            index = j;
         }
      }
      // if (min_deltaR < 0.4) {
      // if (min_deltaR < drmc && MuonSystem -> gLLP_csc[index] == 1) {
      if (min_deltaR < drmc && MuonSystem -> gLLP_csc[index] == 1) {
         gllp_custer_dr -> Fill(min_deltaR);
         cout << "GLLP MATCH: "<< min_deltaR << endl;
         count_gllp++;
         printf("%10.3f %10.3f %10.3f %20.3f %10.3f %10.3f\n", gLLP_pt[index], gLLP_eta[index], gLLP_phi[index], gLLP_decay_vertex_x[index], gLLP_decay_vertex_y[index], gLLP_decay_vertex_z[index]);
         MuonSystem -> cscRechitCluster3_match_gLLP[MuonSystem -> nCscRechitClusters3] = true;
         MuonSystem -> cscRechitCluster3_match_gLLP_minDeltaR[MuonSystem -> nCscRechitClusters3] = min_deltaR;
         MuonSystem -> cscRechitCluster3_match_gLLP_index[MuonSystem -> nCscRechitClusters3] = index;
         MuonSystem -> cscRechitCluster3_match_gLLP_eta[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_eta[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_phi[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_phi[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_decay_r[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_decay_vertex_r[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_decay_x[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_decay_vertex_x[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_decay_y[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_decay_vertex_y[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_decay_z[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_decay_vertex_z[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_ctau[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_ctau[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_beta[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_beta[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_csc[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_csc[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_e[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_e[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_pt[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_pt[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_lepdPhi[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_lepdPhi[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_EMFracP[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_EMFracP[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_EMFracPz[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_EMFracPz[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_EMFracE[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_EMFracE[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_EMFracEz[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_EMFracEz[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_visP[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_visP[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_visPz[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_visPz[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_visE[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_visE[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_visEz[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_visEz[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_daughterKaon[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_daughterKaon[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_multiplicity[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_multiplicity[index];
         MuonSystem -> cscRechitCluster3_match_gLLP_EM_multiplicity[MuonSystem -> nCscRechitClusters3] = MuonSystem -> gLLP_EM_multiplicity[index];
      }
      */


      /*
      bool match = (count_gllp == count_stau);
      cout << "GLLP AND STAU MATCH: " << match << endl;
      cout << count_gllp << " {}{}{} " << count_stau << endl;
      count_gllp = 0;
      count_stau = 0;
      for (int i = 0; i < ncscRechits; i++) {
         if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.4) {
            MuonSystem -> cscRechitCluster3_match_cscRechits_0p4[MuonSystem -> nCscRechitClusters3]++;
         }
         if (!(abs(cscRechitsChamber[i]) == 11 || abs(cscRechitsChamber[i]) == 12)) continue;
         if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.4) {
            MuonSystem -> cscRechitCluster3_match_Me1112_0p4[MuonSystem -> nCscRechitClusters3]++;
            if (abs(cscRechitsChamber[i]) == 11) MuonSystem -> cscRechitCluster3_match_Me11_0p4[MuonSystem -> nCscRechitClusters3]++;
            if (abs(cscRechitsChamber[i]) == 12) MuonSystem -> cscRechitCluster3_match_Me12_0p4[MuonSystem -> nCscRechitClusters3]++;

         }
         if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.6) {
            MuonSystem -> cscRechitCluster3_match_Me1112_0p6[MuonSystem -> nCscRechitClusters3]++;
            if (abs(cscRechitsChamber[i]) == 11) MuonSystem -> cscRechitCluster3_match_Me11_0p6[MuonSystem -> nCscRechitClusters3]++;
            if (abs(cscRechitsChamber[i]) == 12) MuonSystem -> cscRechitCluster3_match_Me12_0p6[MuonSystem -> nCscRechitClusters3]++;
         }
         if (RazorAnalyzer::deltaR(cscRechitsEta[i], cscRechitsPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.8) {
            MuonSystem -> cscRechitCluster3_match_Me1112_0p8[MuonSystem -> nCscRechitClusters3]++;
            if (abs(cscRechitsChamber[i]) == 11) MuonSystem -> cscRechitCluster3_match_Me11_0p8[MuonSystem -> nCscRechitClusters3]++;
            if (abs(cscRechitsChamber[i]) == 12) MuonSystem -> cscRechitCluster3_match_Me12_0p8[MuonSystem -> nCscRechitClusters3]++;
         }

      }
      for (int i = 0; i < nCscSeg; i++) {
         if (RazorAnalyzer::deltaR(cscSegEta[i], cscSegPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.4) {
            MuonSystem -> cscRechitCluster3_match_cscSeg_0p4[MuonSystem -> nCscRechitClusters3]++;
            if (abs(cscSegChamber[i]) == 11) MuonSystem -> cscRechitCluster3_match_ME11Seg_0p4[MuonSystem -> nCscRechitClusters3]++;
            if (abs(cscSegChamber[i]) == 12) MuonSystem -> cscRechitCluster3_match_ME12Seg_0p4[MuonSystem -> nCscRechitClusters3]++;

         }

         if (RazorAnalyzer::deltaR(cscSegEta[i], cscSegPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.6) {
            MuonSystem -> cscRechitCluster3_match_cscSeg_0p6[MuonSystem -> nCscRechitClusters3]++;
            if (abs(cscSegChamber[i]) == 11) MuonSystem -> cscRechitCluster3_match_ME11Seg_0p6[MuonSystem -> nCscRechitClusters3]++;
            if (abs(cscSegChamber[i]) == 12) MuonSystem -> cscRechitCluster3_match_ME12Seg_0p6[MuonSystem -> nCscRechitClusters3]++;
         }
      }
      //match to MB1 DT hits
      for (int i = 0; i < nDtRechits; i++) {
         if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.4) {
            MuonSystem -> cscRechitCluster3_match_dtRechits_0p4[MuonSystem -> nCscRechitClusters3]++;
            if (dtRechitStation[i] == 1) MuonSystem -> cscRechitCluster3_match_MB1_0p4[MuonSystem -> nCscRechitClusters3]++;
         }

         if (RazorAnalyzer::deltaR(dtRechitEta[i], dtRechitPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.6) {
            MuonSystem -> cscRechitCluster3_match_dtRechits_0p6[MuonSystem -> nCscRechitClusters3]++;
            if (dtRechitStation[i] == 1) MuonSystem -> cscRechitCluster3_match_MB1_0p6[MuonSystem -> nCscRechitClusters3]++;
         }
         //match to check for beam halo
         if (RazorAnalyzer::deltaPhi(dtRechitPhi[i], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.2) {
            MuonSystem -> cscRechitCluster3_match_dtRechits_phi0p2[MuonSystem -> nCscRechitClusters3]++;
         }

      }
      //match to MB1 DT segments
      for (int i = 0; i < nDtSeg; i++) {
         if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.4) {
            MuonSystem -> cscRechitCluster3_match_dtSeg_0p4[MuonSystem -> nCscRechitClusters3]++;
            if (dtSegStation[i] == 1) MuonSystem -> cscRechitCluster3_match_MB1Seg_0p4[MuonSystem -> nCscRechitClusters3]++;
         }

         if (RazorAnalyzer::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.6) {
            MuonSystem -> cscRechitCluster3_match_dtSeg_0p6[MuonSystem -> nCscRechitClusters3]++;
            if (dtSegStation[i] == 1) MuonSystem -> cscRechitCluster3_match_MB1Seg_0p6[MuonSystem -> nCscRechitClusters3]++;
         }

      }
      //match to RPC hits in RE1/2
      for (int i = 0; i < nRpc; i++) {
         float rpcR = sqrt(rpcX[i] * rpcX[i] + rpcY[i] * rpcY[i]);
         if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.4) {
            if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730) {
               MuonSystem -> cscRechitCluster3_match_RE12_0p4[MuonSystem -> nCscRechitClusters3]++;
            }
            if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661) {
               MuonSystem -> cscRechitCluster3_match_RB1_0p4[MuonSystem -> nCscRechitClusters3]++;
            }

         }
         if (RazorAnalyzer::deltaR(rpcEta[i], rpcPhi[i], MuonSystem -> cscRechitCluster3Eta[MuonSystem -> nCscRechitClusters3], MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3]) < 0.6) {
            if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730) {
               MuonSystem -> cscRechitCluster3_match_RE12_0p6[MuonSystem -> nCscRechitClusters3]++;
            }
            if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661) {
               MuonSystem -> cscRechitCluster3_match_RB1_0p6[MuonSystem -> nCscRechitClusters3]++;
            }

         }

      }

      MuonSystem -> cscRechitCluster3Met_dPhi[MuonSystem -> nCscRechitClusters3] = RazorAnalyzer::deltaPhi(MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3], MuonSystem -> metPhi);
      MuonSystem -> cscRechitCluster3MetXYCorr_dPhi[MuonSystem -> nCscRechitClusters3] = RazorAnalyzer::deltaPhi(MuonSystem -> cscRechitCluster3Phi[MuonSystem -> nCscRechitClusters3], MuonSystem -> metPhiXYCorr);
      MuonSystem -> nCscRechitClusters3++;
   }
   */
   cout << "+++++++++++++++++++++++++++++[[[[[[[ FILLING ]]]]]]]+++++++++++++++++++++++++++++" << endl;
   if (!isData && signalScan) {
      pair < int, int > smsPair = make_pair(MuonSystem -> mX, MuonSystem -> ctau);
      Trees2D[smsPair] -> Fill();
   } else {
      MuonSystem -> tree_ -> Fill();
   }

}
if (!isData && signalScan) {
   for (auto & filePtr: Files2D) {
      cout << "Writing output tree (" << filePtr.second -> GetName() << ")" << endl;
      filePtr.second -> cd();
      Trees2D[filePtr.first] -> Write();
      NEvents2D[filePtr.first] -> Write("NEvents");
      Total2D[filePtr.first] -> Write("Total");
      accep2D[filePtr.first] -> Write("acceptance");
      accep_met2D[filePtr.first] -> Write("acceptance_met");
      filePtr.second -> Close();

   }
} else if (!isData) {
   cout << "Filled Total of " << NEvents -> GetBinContent(1) << " Events\n";
   cout << "Writing output trees..." << endl;
   outFile -> cd();
   MuonSystem -> tree_ -> Write();
   NEvents -> Write();
   NTriggers -> Write();
   Matched_PT -> Write();

   EVENTS -> Write();
   RECHITS -> Write();
   SEGMENTS -> Write();
   RECHITS_CL -> Write();
   RECHITS_CL_ADD -> Write();
   SEGMENTS_CL -> Write();

   gllp_custer_dr -> Write();
   // Neutrino_tau -> Write();
   tau_pt -> Write();
   // Hadrons -> Write();
   // Leptons -> Write();
   // Neutrino_other -> Write();
   accep -> Write("acceptance");
   accep_met -> Write("acceptance_met");
   outFile -> Close();
} else {
   cout << "Filled Total of " << NEvents -> GetBinContent(1) << " Events\n";
   cout << "Writing output trees..." << endl;
   outFile -> cd();
   MuonSystem -> tree_ -> Write();
   Nmet200 -> Write();
   NmetFilter -> Write();
   Nlep0 -> Write();
   Njet1 -> Write();
   NcosmicVeto -> Write();
   NEvents -> Write();
   outFile -> Write();
   outFile -> Close();
}
}
