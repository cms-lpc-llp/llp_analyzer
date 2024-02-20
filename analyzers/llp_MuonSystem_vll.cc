#include "llp_MuonSystem_vll.h"
#include "RazorHelper.h"
#include "TreeMuonSystemVLL.h"
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


void llp_MuonSystem_vll::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData  = " << isData << "\n";
  cout << "options = " << options << "\n";
  cout << "analisisTag = " << analysisTag << "\n";
  if (analysisTag == ""){ analysisTag = "Razor2018_17SeptEarlyReReco"; }

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


  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");

  TreeMuonSystemVLL *MuonSystem = new TreeMuonSystemVLL;
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  //histogram containing total number of processed events (for normalizationS)
  TH1F *CutFlow       = new TH1F("CutFlow"   , "CutFlow"         , 10, 0, 10);
  TH1F *Acceptance    = new TH1F("Acceptance", "Acceptance"      , 10, 0, 10);
  TH1F *AcceptanceCsc = new TH1F("AcceptanceCsc", "AcceptanceCsc", 10, 0, 10);
  TH1F *AcceptanceDt  = new TH1F("AcceptanceDt", "AcceptanceDt"  , 10, 0, 10);

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

  //*************************************************************************
  //Look over Input File Events
  //*************************************************************************
  if (fChain == 0) return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {
  //for (Long64_t jentry=0; jentry<100000;jentry++) { 
    //begin event
    if(jentry % 10000 == 0) { cout << "Processing entry " << jentry << endl;}
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //GetEntry(ientry);
    nb = fChain->GetEntry(jentry); nbytes += nb;

    //fill normalization histogram
    MuonSystem->InitVariables();

    //GenFilter?
    bool GenVeto          = false;
    bool isVLLsignal      = false;
    bool isWtausignal     = false;

    //Fill gen-level information (gParticle, genMet, PUweights)
    if (!isData)
    {
      //gParticle Collection Matching
      int pid, pmotherid, pgrandmid;
      std::vector<int> all_vlls_idxs,all_llps_idxs,all_taus_idxs,all_gams_idxs,all_elcs_idxs,all_tauele_idxs,all_taumu_idxs;
      //Get VLLs, Taus, LLPs, Gammas
      for (int i=0; i < nGenParticle; i++)
      {
        pid        = abs(gParticleId[i]);
        pmotherid  = abs(gParticleMotherId[i]);
        pgrandmid  = abs(gParticleMotherId[gParticleMotherIndex[i]]);
        if (pid == 5000002) {all_vlls_idxs.push_back(i); isVLLsignal=true;}
        if (pid == 5000001 && pmotherid==5000002) all_llps_idxs.push_back(i);
        if (pid == 15      && pmotherid==5000002) all_taus_idxs.push_back(i);
        if (pid == 22      && pmotherid==5000001) all_gams_idxs.push_back(i);
        if (pid != 22      && pmotherid==22 && pgrandmid==5000001) all_elcs_idxs.push_back(i);
        //Muons/eles from tau-leptonic decays  
        if (pid == 11  && pmotherid==15) all_tauele_idxs.push_back(i);   
        if (pid == 13  && pmotherid==15) all_taumu_idxs.push_back(i);
        //is there a tau from W decay?
        if (pid == 15 && pmotherid==24) {isWtausignal=true;}
      }
      if (isVLLsignal) { //VLL gen info
          //Pair photons to their LLP
          std::vector<std::vector<int>> llps_gammas_idxs;
          std::vector<pair<int,int>> llps_gammas_tmp;
          for (u_int i=0; i < all_gams_idxs.size(); i++)
          {
            int g = all_gams_idxs[i];
            std::pair<int,int> tmp_pair;
            if (abs( gParticleStatus[g] ) == 1) //Stable photon
            {
                   if (gParticleMotherId[g] == gParticleMotherId[gParticleMotherIndex[g]]) tmp_pair = make_pair(gParticleMotherIndex[gParticleMotherIndex[g]],gParticleMotherIndex[g]);
                   else tmp_pair = make_pair(gParticleMotherIndex[g],g); 
                   llps_gammas_tmp.push_back(tmp_pair);
            }
            else{ //Unstable photon
                   bool found=false;
                   for (u_int k=0; k < all_elcs_idxs.size(); k++)
                   {
                       int e = all_elcs_idxs[k];
                       if (found) continue;
                       if (g==abs(gParticleMotherIndex[e]))
                       {
                        tmp_pair = make_pair(gParticleMotherIndex[g],g); llps_gammas_tmp.push_back(tmp_pair);found=true;
                       }
                   }
            }
          }
          if (llps_gammas_tmp.size()!=4) cout<<"ERROR: The number of photons is "<<llps_gammas_tmp.size()<<endl;
          llps_gammas_idxs.push_back( {llps_gammas_tmp[0].first,llps_gammas_tmp[0].second,llps_gammas_tmp[1].second} );
          llps_gammas_idxs.push_back( {llps_gammas_tmp[2].first,llps_gammas_tmp[2].second,llps_gammas_tmp[3].second} );   
          //Pair taus with llps and gammas
          std::vector<int> taus_idxs;
          std::vector<std::vector<int>> tausllpsgammas_idxs;
          for (u_int k=0; k < all_taus_idxs.size(); k++)
          {
              int t = all_taus_idxs[k];
              if (abs(gParticleStatus[t]) == 2)
              {
                     if ( gParticleMotherId[t] == gParticleMotherId[gParticleMotherIndex[t]] ) taus_idxs.push_back(gParticleMotherIndex[t]); 
                     else taus_idxs.push_back(t);
              }
          }
          for (u_int i=0; i < taus_idxs.size(); i++)
          {
                 int g = taus_idxs[i];  
                 for (u_int j=0; j < llps_gammas_idxs.size(); j++)
                    {
                       std::vector<int> k = llps_gammas_idxs[j];
                       if ( abs( gParticleMotherIndex[k[0]] ) == abs(gParticleMotherIndex[g])) tausllpsgammas_idxs.push_back({g,k[0],k[1],k[2]});
                       if ( abs( gParticleMotherIndex[k[0]] ) == abs(gParticleMotherIndex[gParticleMotherIndex[g]])) tausllpsgammas_idxs.push_back({g,k[0],k[1],k[2]}) ;                   
                    }
          }
          if (tausllpsgammas_idxs.size()!=2) cout<<"ERROR: The number of taus/llps is "<<llps_gammas_tmp.size()<<endl;
          //Pair taus,llps,photons with last VLL copy
          std::vector<std::vector<int>> idxs;
          int maxstatus=-1;
          for (u_int i=0; i < all_vlls_idxs.size(); i++)
          { 
              int k = all_vlls_idxs[i]; 
              if( abs(gParticleStatus[k])>maxstatus) maxstatus = abs(gParticleStatus[k]);
          }
          for (u_int i=0; i < all_vlls_idxs.size(); i++)
          {
              int g = all_vlls_idxs[i];
              if ( abs(gParticleStatus[g])==maxstatus)
              {
                  for (u_int i=0; i < tausllpsgammas_idxs.size(); i++)
                  { 
                      std::vector<int> k = tausllpsgammas_idxs[i]; 
                      if ( abs(gParticleMotherIndex[k[1]]) == abs(gParticleMotherIndex[g]) ) idxs.push_back({g,k[0],k[1],k[2],k[3]});
                      if ( abs(gParticleMotherIndex[k[1]]) == abs(gParticleMotherIndex[gParticleMotherIndex[g]]) ) idxs.push_back({g,k[0],k[1],k[2],k[3]});
                  } 
              }
          }
          if (idxs.size()!=2) cout<<"ERROR: The number of VLLs is "<<llps_gammas_tmp.size()<<endl;
          //shuffle VLLs to remove any matching bias
          std::random_shuffle(std::begin(idxs),std::end(idxs)); //fix random number

          //Find muon from taus (if exists)
          int taumu1=-1,taumu2=-1;
          for (u_int i=0; i < all_taumu_idxs.size(); i++)
          { 
             int u = all_taumu_idxs[i];
             if ( abs(gParticleStatus[u]) != 1 ) continue;//stable muon
             if ( abs(gParticleMotherIndex[u]) == idxs[0][1]) taumu1=u;
             if ( abs(gParticleMotherIndex[u]) == idxs[1][1]) taumu2=u;
          }
          std::vector<int> taumus = {taumu1 ,taumu2 };
          //Find ele from taus (if exists)
          int tauele1=-1,tauele2=-1;
          for (u_int i=0; i < all_tauele_idxs.size(); i++)
          { 
             int e = all_tauele_idxs[i];
             if ( abs(gParticleStatus[e]) != 1 ) continue;//stable muon
             if ( abs(gParticleMotherIndex[e]) == idxs[0][1]) tauele1=e;
             if ( abs(gParticleMotherIndex[e]) == idxs[1][1]) tauele2=e;
          }
          std::vector<int> taueles= {tauele1,tauele2};
          //Make the VLL,Tau,LLP p4 
          TLorentzVector VLL1_p4,Tau1_p4,LLP1_p4,VLL2_p4,Tau2_p4,LLP2_p4;
          VLL1_p4.SetPtEtaPhiE(gParticlePt[idxs[0][0]],gParticleEta[idxs[0][0]],gParticlePhi[idxs[0][0]],gParticleE[idxs[0][0]]);
          Tau1_p4.SetPtEtaPhiE(gParticlePt[idxs[0][1]],gParticleEta[idxs[0][1]],gParticlePhi[idxs[0][1]],gParticleE[idxs[0][1]]);
          LLP1_p4.SetPtEtaPhiE(gParticlePt[idxs[0][2]],gParticleEta[idxs[0][2]],gParticlePhi[idxs[0][2]],gParticleE[idxs[0][2]]);
          VLL2_p4.SetPtEtaPhiE(gParticlePt[idxs[1][0]],gParticleEta[idxs[1][0]],gParticlePhi[idxs[1][0]],gParticleE[idxs[1][0]]);
          Tau2_p4.SetPtEtaPhiE(gParticlePt[idxs[1][1]],gParticleEta[idxs[1][1]],gParticlePhi[idxs[1][1]],gParticleE[idxs[1][1]]);
          LLP2_p4.SetPtEtaPhiE(gParticlePt[idxs[1][2]],gParticleEta[idxs[1][2]],gParticlePhi[idxs[1][2]],gParticleE[idxs[1][2]]);
          std::vector<TLorentzVector> VLL_p4 = {VLL1_p4,VLL2_p4};
          std::vector<TLorentzVector> Tau_p4 = {Tau1_p4,Tau2_p4};
          std::vector<TLorentzVector> LLP_p4 = {LLP1_p4,LLP2_p4};
          //gVLL Collection
          for(int i = 0; i < 2;i++)
          {
             MuonSystem->gVLL_eta[i] = VLL_p4[i].Eta(); 
             MuonSystem->gVLL_phi[i] = VLL_p4[i].Phi();
             MuonSystem->gVLL_e[i]   = VLL_p4[i].E();  
             MuonSystem->gVLL_pt[i]  = VLL_p4[i].Pt(); 
             MuonSystem->gTau_eta[i] = Tau_p4[i].Eta(); 
             MuonSystem->gTau_phi[i] = Tau_p4[i].Phi();
             MuonSystem->gTau_e[i]   = Tau_p4[i].E();  
             MuonSystem->gTau_pt[i]  = Tau_p4[i].Pt(); 
             MuonSystem->gLLP_eta[i] = LLP_p4[i].Eta(); 
             MuonSystem->gLLP_phi[i] = LLP_p4[i].Phi();
             MuonSystem->gLLP_e[i]   = LLP_p4[i].E();  
             MuonSystem->gLLP_pt[i]  = LLP_p4[i].Pt(); 
          }  
          //gTau Collection variables
          for(int i = 0; i < 2;i++)
          {
              if (taumus[i]<0  && taueles[i]<0)
              {
                   MuonSystem->gTau_decaymode[i] = 0; //Tau hadronic
              } 
              else if (taumus[i]<0  && taueles[i]>=0)
              {
                MuonSystem->gTau_decaymode[i] = 1; //Tau electron
                MuonSystem->gTau_ele_e[i]     = gParticleE[taueles[i]];
                MuonSystem->gTau_ele_pt[i]    = gParticlePt[taueles[i]];
                MuonSystem->gTau_ele_eta[i]   = gParticleEta[taueles[i]];
                MuonSystem->gTau_ele_phi[i]   = gParticlePhi[taueles[i]];
              } 
              else if (taumus[i]>=0 && taueles[i]<0) 
              {
                MuonSystem->gTau_decaymode[i] = 2; //Tau muon
                MuonSystem->gTau_mu_e[i]      = gParticleE[taumus[i]];
                MuonSystem->gTau_mu_pt[i]     = gParticlePt[taumus[i]];
                MuonSystem->gTau_mu_eta[i]    = gParticleEta[taumus[i]];
                MuonSystem->gTau_mu_phi[i]    = gParticlePhi[taumus[i]];
              } 
              else
              {  
              //nothing
              }
          }

          //gLLP Collection additional variables
          MuonSystem->gLLP_decay_vertex_x[0] = gParticleProdVertexX[idxs[0][3]]; MuonSystem->gLLP_decay_vertex_x[1] = gParticleProdVertexX[idxs[1][3]];
          MuonSystem->gLLP_decay_vertex_y[0] = gParticleProdVertexY[idxs[0][3]]; MuonSystem->gLLP_decay_vertex_y[1] = gParticleProdVertexY[idxs[1][3]];
          MuonSystem->gLLP_decay_vertex_z[0] = gParticleProdVertexZ[idxs[0][3]]; MuonSystem->gLLP_decay_vertex_z[1] = gParticleProdVertexZ[idxs[1][3]];
          MuonSystem->gLLP_decay_vertex_r[0] = sqrt(MuonSystem->gLLP_decay_vertex_x[0]*MuonSystem->gLLP_decay_vertex_x[0]+MuonSystem->gLLP_decay_vertex_y[0]*MuonSystem->gLLP_decay_vertex_y[0]);
          MuonSystem->gLLP_decay_vertex_r[1] = sqrt(MuonSystem->gLLP_decay_vertex_x[1]*MuonSystem->gLLP_decay_vertex_x[1]+MuonSystem->gLLP_decay_vertex_y[1]*MuonSystem->gLLP_decay_vertex_y[1]);      //    float beta = gLLP_beta[i];
          float gLLP1_decay_length = sqrt( pow(MuonSystem->gLLP_decay_vertex_x[0]-gParticleProdVertexX[idxs[0][2]], 2) + pow(MuonSystem->gLLP_decay_vertex_y[0]-gParticleProdVertexY[idxs[0][2]], 2) + pow(MuonSystem->gLLP_decay_vertex_z[0]-gParticleProdVertexZ[idxs[0][2]], 2) );
          float gLLP2_decay_length = sqrt( pow(MuonSystem->gLLP_decay_vertex_x[1]-gParticleProdVertexX[idxs[1][2]], 2) + pow(MuonSystem->gLLP_decay_vertex_y[1]-gParticleProdVertexY[idxs[1][2]], 2) + pow(MuonSystem->gLLP_decay_vertex_z[1]-gParticleProdVertexZ[idxs[1][2]], 2) );
          float gLLP1_beta = LLP1_p4.Beta(); float gLLP1_gamma = LLP1_p4.Gamma();  
          float gLLP2_beta = LLP2_p4.Beta(); float gLLP2_gamma = LLP2_p4.Gamma();
          MuonSystem->gLLP_decaylength[0] = gLLP1_decay_length;   
          MuonSystem->gLLP_decaylength[1] = gLLP2_decay_length;  
          MuonSystem->gLLP_decaylength_rest[0] = gLLP1_decay_length/(gLLP1_beta*gLLP1_gamma);   
          MuonSystem->gLLP_decaylength_rest[1] = gLLP2_decay_length/(gLLP2_beta*gLLP2_gamma);  
          MuonSystem->gLLP_beta[0]  = gLLP1_beta;   
          MuonSystem->gLLP_beta[1]  = gLLP2_beta;  
          MuonSystem->gLLP_gamma[0] = gLLP1_gamma;
          MuonSystem->gLLP_gamma[1] = gLLP2_gamma;
          if ( MuonSystem->gLLP_decaylength_rest[0] < 2.0 || MuonSystem->gLLP_decaylength_rest[1] < 2.0) GenVeto=true;
          for(int i = 0; i < 2;i++)
          {
              //LLP detector tags
              if (abs(MuonSystem->gLLP_eta[i]) < 2.4
                  && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>400
                  && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) MuonSystem->gLLP_csc[i] = true;
              if (abs(MuonSystem->gLLP_decay_vertex_z[i])< 661.0
                  && MuonSystem->gLLP_decay_vertex_r[i] < 800
                   && MuonSystem->gLLP_decay_vertex_r[i] > 200.0) MuonSystem->gLLP_dt[i] = true;
          }
      }//end vllsignal info

      //gen MET variables
      MuonSystem->genMetPtTrue = genMetPtTrue;
      MuonSystem->genMetPhiTrue = genMetPhiTrue;
      //MC PU info
      for (int i=0; i < nBunchXing; i++)
      {
          if (BunchXing[i] == 0) MuonSystem->npu = nPUmean[i];
      }
      MuonSystem->pileupWeight     = helper->getPileupWeight(MuonSystem->npu);
      MuonSystem->pileupWeightUp   = helper->getPileupWeightUp(MuonSystem->npu) / MuonSystem->pileupWeight;
      MuonSystem->pileupWeightDown = helper->getPileupWeightDown(MuonSystem->npu) / MuonSystem->pileupWeight;      
    }//end of MC genlevel infos
 

    //Veto events using genveto 
    if (GenVeto) continue;
    //Veto events without a tau from W decay
    //if (!isData && !isVLLsignal && !isWtausignal) continue;

    MuonSystem->weight = 1;//unweighted for now //FIX ME
    MuonSystem->runNum = runNum;
    MuonSystem->lumiSec = lumiNum;
    MuonSystem->evtNum = eventNum;

    //Fill MC conditions
    if (!isData)
    {
        if (analysisTag=="Razor2016_07Aug2017Rereco" || "Razor2016_UL" || "Razor2016APV_UL") MuonSystem->MC_condition = 2016;
        else if (analysisTag=="Razor2017_17Nov2017Rereco" || "Razor2017_UL") MuonSystem->MC_condition = 2017;
        else if (analysisTag=="Razor2018_17SeptEarlyReReco" || "Razor2018_UL") MuonSystem->MC_condition = 2018;
    }

    //Fill NPU info
    MuonSystem->npv = nPV;

    //Fill MET info
    MuonSystem->met        = metType1Pt;
    MuonSystem->metPhi     = metType1Phi;
    MuonSystem->metJESUp   = MuonSystem->met;
    MuonSystem->metJESDown = MuonSystem->met;
    MuonSystem->metSF      = helper->getMetTriggerSF(MuonSystem->met);

    std::pair<double,double> corrected_met;
    if (analysisTag=="Razor2016_07Aug2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2016, !isData, nPV);
    else if (analysisTag=="Razor2017_17Nov2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2017, !isData, nPV);
    else if (analysisTag=="Razor2018_17SeptEarlyReReco") corrected_met = helper->METXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, 2018, !isData, nPV);
    else if (analysisTag=="Razor2016_UL") corrected_met = helper->ULMETXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, "2016nonAPV", !isData, nPV, true);
    else if (analysisTag=="Razor2016APV_UL") corrected_met = helper->ULMETXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, "2016APV", !isData, nPV, true);
    else if (analysisTag=="Razor2017_UL") corrected_met = helper->ULMETXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, "2017", !isData, nPV, true);
    else if (analysisTag=="Razor2018_UL") corrected_met = helper->ULMETXYCorr_Met_MetPhi(metType1Pt, metType1Phi, runNum, "2018", !isData, nPV, true);
    MuonSystem->metXYCorr     = corrected_met.first;
    MuonSystem->metEENoise    = MuonSystem->metXYCorr;
    MuonSystem->metHEM        = MuonSystem->metXYCorr;
    MuonSystem->metPhiXYCorr  = corrected_met.second;
    MuonSystem->metPhiEENoise = MuonSystem->metPhiXYCorr;
    MuonSystem->metPhiHEM     = MuonSystem->metPhiXYCorr;

    //Trigger bits
    for(int i = 0; i < NTriggersMAX; i++){
      MuonSystem->HLTDecision[i] = HLTDecision[i];
    }
    if (analysisTag=="Razor2016_07Aug2017Rereco" || analysisTag=="Razor2016_UL" || analysisTag=="Razor2016APV_UL" )
    {
      MuonSystem->METTrigger = HLTDecision[310] || HLTDecision[467];
      MuonSystem->METNoMuTrigger = HLTDecision[467];
    }
    else
    {
      MuonSystem->METTrigger = HLTDecision[310] || HLTDecision[467] || HLTDecision[703] || HLTDecision[717] || HLTDecision[710] || HLTDecision[709];
      MuonSystem->METNoMuTrigger = HLTDecision[467] ||  HLTDecision[717] || HLTDecision[710];
    }

    //Fill Noise Filter Flag
    //NOTE: REMOVE Flag2_HBHEIsoNoiseFilter
    MuonSystem->Flag2_all = Flag2_HBHENoiseFilter && Flag2_BadPFMuonFilter && Flag2_globalSuperTightHalo2016Filter && Flag2_EcalDeadCellTriggerPrimitiveFilter;
    MuonSystem->Flag2_HBHEIsoNoiseFilter= Flag2_HBHEIsoNoiseFilter;
    if (isData) MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_eeBadScFilter;
    if (analysisTag!="Razor2016_07Aug2017Rereco" || analysisTag!="Razor2016_UL" || analysisTag!="Razor2016APV_UL") MuonSystem->Flag2_all = MuonSystem->Flag2_all && Flag2_ecalBadCalibFilter;

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
        if(muonPt[i] < 25) continue;
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
        if (!isEGammaPOGLooseElectron(i, true, false, true, "vid")) continue;
        if(elePt[i] < 35) continue;
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
    for ( auto &tmp : Leptons )
    {
        MuonSystem->lepE[MuonSystem->nLeptons]            = tmp.lepton.E();
        MuonSystem->lepPt[MuonSystem->nLeptons]           = tmp.lepton.Pt();
        MuonSystem->lepEta[MuonSystem->nLeptons]          = tmp.lepton.Eta();
        MuonSystem->lepPhi[MuonSystem->nLeptons]          = tmp.lepton.Phi();
        MuonSystem->lepPdgId[MuonSystem->nLeptons]        = tmp.pdgId;
        MuonSystem->lepPassId[MuonSystem->nLeptons]       = tmp.passId; 
        MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
        MuonSystem->lepDZ[MuonSystem->nLeptons]           = tmp.dZ;
        MuonSystem->nLeptons++;
    }

    //-----------------------------------------------
    //Jets
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

        //Era dependent special corrections to MET (HEM=Remove jets in affected region, EENoise:Remove noise jets, EEprefiring: Tag events with jets in prefiring region)
        if (thisJet.Eta()>-3.0 && thisJet.Eta()<-1.3 && thisJet.Phi() >-1.57 && thisJet.Phi() <-0.87 && (analysisTag == "Razor2018_17SeptEarlyReReco" || analysisTag == "Razor2018_UL") )
        {
          MetXCorr_HEM += thisJet.Px();
          MetYCorr_HEM += thisJet.Py();
        }
        if (fabs(thisJet.Eta())> 2.65 && fabs(thisJet.Eta())<3.139 && thisJet.Pt() < 50  && (analysisTag == "Razor2017_17Nov2017Rereco" || analysisTag == "Razor2017_UL" )  )
        {
          MetXCorr_EENoise += thisJet.Px();
          MetYCorr_EENoise += thisJet.Py();
        }
        if (fabs(thisJet.Eta())> 2.25 && fabs(thisJet.Eta())<3.0 && thisJet.Pt() > 100 )
        {
          MuonSystem->EE_prefiring = false;
        }

        if (fabs(thisJet.Eta()) >= 3.0) continue;
        jets tmpJet;
        tmpJet.jet    = thisJet;
        tmpJet.passId = isPFTightJet(i, true,analysisTag);

        //JES variations  
        double unc = helper->getJecUnc( jetCorrPt, jetEta[i], runNum ); //use run=999 as default
        tmpJet.jetPtJESUp = jetCorrPt*(1+unc);
        tmpJet.jetPtJESDown = jetCorrPt*(1-unc);
        tmpJet.jetEJESUp = jetCorrE*(1+unc);
        tmpJet.jetEJESDown = jetCorrE*(1-unc);
        tmpJet.JecUnc = unc;
        TLorentzVector thisJetJESUp = makeTLorentzVector(tmpJet.jetPtJESUp, jetEta[i], jetPhi[i], tmpJet.jetEJESUp);
        TLorentzVector thisJetJESDown = makeTLorentzVector(tmpJet.jetPtJESDown, jetEta[i], jetPhi[i], tmpJet.jetEJESDown);
        if (tmpJet.jetPtJESUp > 10)
        {
          MetXCorr_JESUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
          MetYCorr_JESUp += -1 * (thisJetJESUp.Py() - thisJet.Py());
        }
        if (tmpJet.jetPtJESDown > 10)
        {
          MetXCorr_JESDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
          MetYCorr_JESDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
        }
        if( thisJet.Pt() < 30 ) continue; //According to the April 1st 2015 AN
        Jets.push_back(tmpJet);
    }


    //Fill Jet information
    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);
    double jetMet_dPhiMin_temp = 999.;
    for ( auto &tmp : Jets )
    {
      if(tmp.jet.Pt()<30) continue;
      if(abs(tmp.jet.Eta())>2.4) continue;
      MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
      MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
      MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
      MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
      MuonSystem->jetPtJESUp[MuonSystem->nJets] = tmp.jetPtJESUp;
      MuonSystem->jetPtJESDown[MuonSystem->nJets] = tmp.jetPtJESDown;
      MuonSystem->jetEJESUp[MuonSystem->nJets] = tmp.jetEJESUp;
      MuonSystem->jetEJESDown[MuonSystem->nJets] = tmp.jetEJESDown;
      if (jetMet_dPhiMin_temp > abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi)))
      {
        if (tmp.jet.Pt()>30 && abs(tmp.jet.Eta())<2.4)
        {
          jetMet_dPhiMin_temp = abs(RazorAnalyzer::deltaPhi(tmp.jet.Phi(),metType1Phi));
        }
      }
      MuonSystem->jetTightPassId[MuonSystem->nJets] = tmp.passId;
      MuonSystem->nJets++;
    }
    MuonSystem-> jetMet_dPhiMin = jetMet_dPhiMin_temp;

    //Fill era dependent special corrections to MET + JES variation effects to MET: Nominal branch is called EENoise (Christina's)
    TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(MuonSystem->metXYCorr, 0, MuonSystem->metPhiXYCorr, 0);
    //JES up
    float PFMetXJESUp          = PFMET.Px() + MetXCorr_JESUp;
    float PFMetYJESUp          = PFMET.Py() + MetYCorr_JESUp;
    MuonSystem->metJESUp       = sqrt( pow(PFMetXJESUp,2) + pow(PFMetYJESUp,2) );
    MuonSystem->metPhiJESUp    = atan(PFMetYJESUp/PFMetXJESUp);
    if  (PFMetXJESUp < 0.0) MuonSystem->metPhiJESUp = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiJESUp,0.0);
    MuonSystem->metJESUpSF     = helper->getMetTriggerSF(MuonSystem->metJESUp)/MuonSystem->metSF;
    //JES down
    float PFMetXJESDown        = PFMET.Px() + MetXCorr_JESDown;
    float PFMetYJESDown        = PFMET.Py() + MetYCorr_JESDown;
    MuonSystem->metJESDown     =  sqrt( pow(PFMetXJESDown,2) + pow(PFMetYJESDown,2) );
    MuonSystem->metPhiJESDown  = atan(PFMetYJESDown/PFMetXJESDown);
    if  (PFMetXJESUp < 0.0) MuonSystem->metPhiJESDown = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiJESDown,0.0);
    MuonSystem->metJESDownSF   = helper->getMetTriggerSF(MuonSystem->metJESDown)/MuonSystem->metSF;
    //HEM
    float PFMetXHEM            = PFMET.Px() + MetXCorr_HEM;
    float PFMetYHEM            = PFMET.Py() + MetYCorr_HEM;
    MuonSystem->metHEM         = sqrt( pow(PFMetXHEM,2) + pow(PFMetYHEM,2) );
    MuonSystem->metPhiHEM      = atan(PFMetYHEM/PFMetXHEM);
    if  (PFMetXHEM < 0.0) MuonSystem->metPhiHEM = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiHEM,0.0);
    //EENoise
    float PFMetXEENoise        = PFMET.Px() + MetXCorr_EENoise;
    float PFMetYEENoise        = PFMET.Py() + MetYCorr_EENoise;
    MuonSystem->metEENoise     = sqrt( pow(PFMetXEENoise,2) + pow(PFMetYEENoise,2) );
    MuonSystem->metPhiEENoise  = atan(PFMetYEENoise/PFMetXEENoise);
    if  (PFMetXEENoise < 0.0) MuonSystem->metPhiEENoise = RazorAnalyzer::deltaPhi(TMath::Pi() + MuonSystem->metPhiEENoise,0.0);

    //-----------------------------------------------
    //DT RecHits
    //-----------------------------------------------
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
        if (dtRechitWheel[i]   ==-2) MuonSystem->nDTRechitsWheelMinus2++;
        if (dtRechitWheel[i]   ==-1) MuonSystem->nDTRechitsWheelMinus1++;
        if (dtRechitWheel[i]   == 0) MuonSystem->nDTRechitsWheel0++;
        if (dtRechitWheel[i]   == 1) MuonSystem->nDTRechitsWheelPlus1++;
        if (dtRechitWheel[i]   == 2) MuonSystem->nDTRechitsWheelPlus2++;
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

    //-----------------------------------------------
    //CSC RecHits
    //-----------------------------------------------
    vector<Point> points;
    vector<int> cscRechitsClusterId;
    points.clear();
    MuonSystem->nCscRechits  = 0;
    for (int i = 0; i < ncscRechits; i++) {
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

    //-----------------------------------------------
    //CSC DBSCAN Clustering
    //-----------------------------------------------
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
        MuonSystem->cscRechitClusterNLayersChamberPlus11[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus11;
        MuonSystem->cscRechitClusterNLayersChamberPlus12[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus12;
        MuonSystem->cscRechitClusterNLayersChamberPlus13[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus13;
        MuonSystem->cscRechitClusterNLayersChamberPlus21[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus21;
        MuonSystem->cscRechitClusterNLayersChamberPlus22[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus22;
        MuonSystem->cscRechitClusterNLayersChamberPlus31[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus31;
        MuonSystem->cscRechitClusterNLayersChamberPlus32[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus32;
        MuonSystem->cscRechitClusterNLayersChamberPlus41[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus41;
        MuonSystem->cscRechitClusterNLayersChamberPlus42[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberPlus42;
        MuonSystem->cscRechitClusterNLayersChamberMinus11[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus11;
        MuonSystem->cscRechitClusterNLayersChamberMinus12[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus12;
        MuonSystem->cscRechitClusterNLayersChamberMinus13[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus13;
        MuonSystem->cscRechitClusterNLayersChamberMinus21[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus21;
        MuonSystem->cscRechitClusterNLayersChamberMinus22[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus22;
        MuonSystem->cscRechitClusterNLayersChamberMinus31[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus31;
        MuonSystem->cscRechitClusterNLayersChamberMinus32[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus32;
        MuonSystem->cscRechitClusterNLayersChamberMinus41[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus41;
        MuonSystem->cscRechitClusterNLayersChamberMinus42[MuonSystem->nCscRechitClusters] = tmp.nLayersChamberMinus42;
        //Jet veto veto
        MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterTightJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
        MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;
        //Jet veto
        for(int i = 0; i < nJets; i++)
        {
          if (fabs(jetEta[i]>3.0)) continue;
          if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] ) {
            MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]   = jetPt[i];
            MuonSystem->cscRechitClusterJetVetoEta[MuonSystem->nCscRechitClusters]  = jetEta[i];
            MuonSystem->cscRechitClusterJetVetoPhi[MuonSystem->nCscRechitClusters]  = jetPhi[i];
            double unc = helper->getJecUnc( jetPt[i], jetEta[i], runNum ); //use run=999 as default
            MuonSystem->cscRechitClusterJetVetoPtJESUp[MuonSystem->nCscRechitClusters]  = jetPt[i]*(1+unc);
            MuonSystem->cscRechitClusterJetVetoPtJESDown[MuonSystem->nCscRechitClusters]  = jetPt[i]*(1-unc);    
          }
          //if (!isPFTightJet(i, true,analysisTag)) continue;
          //if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterTightJetVetoPt[MuonSystem->nCscRechitClusters] ) {
          //  MuonSystem->cscRechitClusterTightJetVetoPt[MuonSystem->nCscRechitClusters]  = jetPt[i];
          //  double unc = helper->getJecUnc( jetPt[i], jetEta[i], runNum ); //use run=999 as default
          //  MuonSystem->cscRechitClusterTightJetVetoPtJESUp[MuonSystem->nCscRechitClusters]  = jetPt[i]*(1+unc);
          //  MuonSystem->cscRechitClusterTightJetVetoPtJESDown[MuonSystem->nCscRechitClusters]  = jetPt[i]*(1-unc);
          //}
      }
      //MuonVeto
      for(int i = 0; i < nMuons; i++)
      {
           if (fabs(muonEta[i]>3.0)) continue;
           if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] ) {
             MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]       = muonPt[i];
             MuonSystem->cscRechitClusterMuonVetoPhi[MuonSystem->nCscRechitClusters]      = muonPhi[i];
             MuonSystem->cscRechitClusterMuonVetoEta[MuonSystem->nCscRechitClusters]      = muonEta[i];
             MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters]   = muon_isGlobal[i];
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
      MuonSystem->cscRechitClusterMetEENoise_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhiEENoise);
      MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters]        =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhi);
      MuonSystem->cscRechitClusterMetJESUp_dPhi[MuonSystem->nCscRechitClusters]   =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhiJESUp);
      MuonSystem->cscRechitClusterMetJESDown_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhiJESDown);
      MuonSystem->cscRechitClusterMetXYCorr_dPhi[MuonSystem->nCscRechitClusters]  =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhiXYCorr);
      MuonSystem->cscRechitClusterMetHEM_dPhi[MuonSystem->nCscRechitClusters]     =  RazorAnalyzer::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhiHEM);
      MuonSystem->nCscRechitClusters++;
    }
    //-----------------------------------------------
    //DT DBSCAN Clustering
    //-----------------------------------------------
    MuonSystem->nDtSeg = 0;
    for(int d=0; d<nDtSeg; d++)
    {
      MuonSystem->dtSegPhi[MuonSystem->nDtSeg] = dtSegPhi[d];
      MuonSystem->dtSegEta[MuonSystem->nDtSeg] = dtSegEta[d];
      MuonSystem->dtSegX[MuonSystem->nDtSeg] = dtSegX[d];
      MuonSystem->dtSegY[MuonSystem->nDtSeg] = dtSegY[d];
      MuonSystem->dtSegZ[MuonSystem->nDtSeg] = dtSegZ[d];
      MuonSystem->dtSegStation[MuonSystem->nDtSeg] = dtSegStation[d];
      MuonSystem->nDtSeg++;
    }
    points.clear();
    for (int i = 0; i < nDtRechits; i++) {
      Point p;
//Old samples
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
//New samples
//      p.phi = dtRechitCorrectPhi[i];
//      p.eta = dtRechitCorrectEta[i];
//      p.x = dtRechitCorrectX[i];
//      p.y = dtRechitCorrectY[i];
//      p.z = dtRechitCorrectZ[i];
//      p.t = dtRechitTime[i];
//      p.twire = dtRechitTime[i];
//      p.station = dtRechitStation[i];
//      p.chamber = dtRechitWheel[i];
//      // p.superlayer = 999;
//      p.superlayer = dtRechitSuperLayer[i];
//      p.superlayer = 0;

      p.clusterID = UNCLASSIFIED;
      points.push_back(p);
    }
    int min_point_dt = 50;//isData?50:30;  //minimum number of segments to call it a cluster
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
          for(int i = 0; i < MuonSystem->nCscRechitClusters; i++){
              if (RazorAnalyzer::deltaR(MuonSystem->cscRechitClusterEta[i],MuonSystem->cscRechitClusterPhi[i],tmp.eta, tmp.phi)<0.4) overlap = true;
          }
	        if (overlap) continue;
          for (unsigned int j = 0; j < tmp.segment_id.size(); j++){
              MuonSystem->dtRechitsClusterId[tmp.segment_id[j]] = MuonSystem->nDtRechitClusters;
          }
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
          // dt noise 2016 remove
          if (isData && runNum >= 275750 && runNum <= 275850 )
          {
            if (tmp.phi < TMath::Pi()/4 && tmp.phi > TMath::Pi()/12 && tmp.maxStation == 2 && MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] ==1)
            {
              MuonSystem->dtRechitClusterNoiseVeto[MuonSystem->nDtRechitClusters] = true;
            }
          }
          //Jet veto/ muon veto
          MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterTightJetVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          // jet veto
          for(int i = 0; i < nJets; i++) {
              if (fabs(jetEta[i]>3.0)) continue;
              if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] ) {
                MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters]   = jetPt[i];
                MuonSystem->dtRechitClusterJetVetoEta[MuonSystem->nDtRechitClusters]  = jetEta[i];
                MuonSystem->dtRechitClusterJetVetoPhi[MuonSystem->nDtRechitClusters]  = jetPhi[i];
                double unc = helper->getJecUnc( jetPt[i], jetEta[i], runNum ); //use run=999 as default
                MuonSystem->dtRechitClusterJetVetoPtJESUp[MuonSystem->nDtRechitClusters]  = jetPt[i]*(1+unc);
                MuonSystem->dtRechitClusterJetVetoPtJESDown[MuonSystem->nDtRechitClusters]  = jetPt[i]*(1-unc);
              }
              //if (!isPFTightJet(i, true,analysisTag)) continue;
              //if (RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->dtRechitClusterTightJetVetoPt[MuonSystem->nDtRechitClusters] ) {
              //  MuonSystem->dtRechitClusterTightJetVetoPt[MuonSystem->nDtRechitClusters]  = jetPt[i];
              //  double unc = helper->getJecUnc( jetPt[i], jetEta[i], runNum ); //use run=999 as default
              //  MuonSystem->dtRechitClusterTightJetVetoPtJESUp[MuonSystem->nDtRechitClusters]  = jetPt[i]*(1+unc);
              //  MuonSystem->dtRechitClusterTightJetVetoPtJESDown[MuonSystem->nDtRechitClusters]  = jetPt[i]*(1-unc);
              //}
          }
          //Muon Veto
          for(int i = 0; i < nMuons; i++){
              if (fabs(muonEta[i]>3.0)) continue;
              if (RazorAnalyzer::deltaR(muonEta[i], muonPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && muonPt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] ) {
                MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters]       = muonPt[i];
                MuonSystem->dtRechitClusterMuonVetoPhi[MuonSystem->nDtRechitClusters]      = muonPhi[i];
                MuonSystem->dtRechitClusterMuonVetoEta[MuonSystem->nDtRechitClusters]      = muonEta[i];
                MuonSystem->dtRechitClusterMuonVetoGlobal[MuonSystem->nDtRechitClusters]   = muon_isGlobal[i];
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
          MuonSystem->nDtSeg=nDtSeg;
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
          for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++){
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
          }
          MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters] = max_bx;
          MuonSystem->dtRechitClusterMetEENoise_dPhi[MuonSystem->nDtRechitClusters]      =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhiEENoise);
          MuonSystem->dtRechitClusterMet_dPhi[MuonSystem->nDtRechitClusters]             =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhi);
          MuonSystem->dtRechitClusterMetJESUp_dPhi[MuonSystem->nDtRechitClusters]        =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhiJESUp);
          MuonSystem->dtRechitClusterMetJESDown_dPhi[MuonSystem->nDtRechitClusters]      =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhiJESDown);
          MuonSystem->dtRechitClusterMetXYCorr_dPhi[MuonSystem->nDtRechitClusters]       =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhiXYCorr);
          MuonSystem->dtRechitClusterMetHEM_dPhi[MuonSystem->nDtRechitClusters]          =  RazorAnalyzer::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhiHEM);
          MuonSystem->nDtRechitClusters++;
    }
    for(int i = 0; i < MuonSystem->nDtRechitClusters; i++){
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
    //Match LLPs to the closest CSC or DT clusters
    if (!isData)
    {
        //CSC loop
        for(int j = 0; j < 2; j++)
        {
          float minDelta_csc = 0.4;
          int index_csc=-999;
          for(int i = 0; i < MuonSystem->nCscRechitClusters; i++)
          {
            double current_delta_r = abs(RazorAnalyzer::deltaR(MuonSystem->gLLP_eta[j],MuonSystem->gLLP_phi[j], MuonSystem->cscRechitClusterEta[i], MuonSystem->cscRechitClusterPhi[i]));
            if (current_delta_r < minDelta_csc ) {  minDelta_csc = current_delta_r; index_csc = i;}
          }
          if (index_csc >= 0) {MuonSystem->gLLP_match_cscRechitCluster[j] = 1;}
          else{MuonSystem->gLLP_match_cscRechitCluster[j] = 0;}
          MuonSystem->gLLP_match_cscRechitClusterSize[j] = MuonSystem->cscRechitClusterSize[index_csc];
          MuonSystem->gLLP_match_cscRechitClusterIdx[j] = index_csc;
        }
        //DT loop
        for(int j = 0; j < 2; j++)
        {
          float minDelta_dt = 0.4;
          int index_dt=-999;
          for(int i = 0; i < MuonSystem->nDtRechitClusters; i++)
          {
            double current_delta_r = abs(RazorAnalyzer::deltaR(MuonSystem->gLLP_eta[j],MuonSystem->gLLP_phi[j], MuonSystem->dtRechitClusterEta[i], MuonSystem->dtRechitClusterPhi[i]));
            if (current_delta_r < minDelta_dt ) {  minDelta_dt = current_delta_r; index_dt = i;}
          }
          if (index_dt >= 0) {MuonSystem->gLLP_match_dtRechitCluster[j] = 1;}
          else{MuonSystem->gLLP_match_dtRechitCluster[j] = 0;}
          MuonSystem->gLLP_match_dtRechitClusterSize[j]  = MuonSystem->dtRechitClusterSize[index_dt];
          MuonSystem->gLLP_match_dtRechitClusterIdx[j]   = index_dt;
        }
    }

    //-----------------------------------------------
    //Event selection and cutflow histogram 
    //-----------------------------------------------
    //Bin 0: Total number of events
    CutFlow->Fill(0);
    
    //Bin 0 Acceptance or NoCuts
    bool acceptance=false;
    bool acceptancecsc=false;    
    bool acceptancedt=false;    
    if (!isData)
    {
      if (MuonSystem->gLLP_csc[0]==true || MuonSystem->gLLP_csc[1]==true) acceptancecsc=true;
      if (MuonSystem->gLLP_dt[0]==true  || MuonSystem->gLLP_dt[1]==true)  acceptancedt=true;
      if (acceptancecsc==true || acceptancedt==true) acceptance=true;
    }
    if(acceptance) Acceptance->Fill(0);
    if(acceptancecsc) AcceptanceCsc->Fill(0);
    if(acceptancedt) AcceptanceDt->Fill(0);

    //Fill Trees
    MuonSystem->tree_->Fill();
  }
  if (!isData)
  {
     cout << "Total events processed " << CutFlow->GetBinContent(1) << " Events\n";
     cout << "Writing output trees..." << endl;
     outFile->cd();
     MuonSystem->tree_->Write();
     CutFlow->Write();
     Acceptance->Write();
     AcceptanceCsc->Write();
     AcceptanceDt->Write();
     outFile->Close();
  }
  else
  {
    cout << "Total events processed " << CutFlow->GetBinContent(1) << " Events\n";
    cout << "Writing output trees..." << endl;
    outFile->cd();
    MuonSystem->tree_->Write();
    CutFlow->Write();
    outFile->Close();
  }
}
