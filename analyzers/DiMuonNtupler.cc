#include "DiMuonNtupler.h"
#include "DiMuonTree.h"
#include "JetCorrectorParameters.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

struct MuonCandidate
{
  int   Index;
  TLorentzVector muon;
  int muonCharge;
  int isTightMuon;
  float mu_chargedMiniIso;
  float mu_photonAndNeutralHadronMiniIso;
  float mu_d0;
  float mu_dZ;
  float mu_ip3d;
  float mu_ip3dSignificance;
};

bool _debug = false;

void DiMuonNtupler::Analyze(bool isData, int Option, string outputfilename, string label)
{
    //initialization: create one TTree for each analysis box
  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/RazorAnalyzer/data/JEC/";
  cout << "Getting JEC parameters from " << pathname << endl;

  std::vector<JetCorrectorParameters> correctionParameters;
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
  correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));

  FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(correctionParameters);

    cout << "Initializing..." << endl;
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "DiMuonNtuple.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    DiMuonTree* mumuTree = new DiMuonTree;
    mumuTree->CreateTree();
    mumuTree->tree_->SetAutoFlush(0);

    std::cout << "Run With Option = " << Option << "\n";

    UInt_t NDiMuonsFilled = 0;

    //begin loop
    if (fChain == 0) return;
    UInt_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;

    cout << "nentries = " << nentries << "\n";
    for (UInt_t jentry=0; jentry<nentries;jentry++) {
      //begin event
      if(jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //****************************************
      //Tree entries based on reco objects
      //****************************************
      if (Option < 10 ) {
        //--------------
        //good muon selection
        //--------------
        vector<TLorentzVector> GoodMuons;
        std::vector< MuonCandidate > muCand;
        unsigned int nLooseMuons = 0;
        unsigned int nTightMuons = 0;
        for( int i = 0; i < nMuons; i++ )
        {
          //TLorentzVector for this muon
          TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]);
          if(!isVetoMuon(i)) continue;
          if(muonPt[i] < 20) continue;
          if(abs(muonEta[i]) > 2.4) continue;
          nLooseMuons++;
          GoodMuons.push_back(thisMuon);
          if( isTightMuon(i) ) nTightMuons++;
          //Filling Muon Candidate
          MuonCandidate tmp_muCand;
          tmp_muCand.Index = i;
          tmp_muCand.muon = thisMuon;
          tmp_muCand.muonCharge = muonCharge[i];
          tmp_muCand.isTightMuon = isTightMuon(i);
          tmp_muCand.mu_chargedMiniIso = muon_chargedMiniIso[i];
          tmp_muCand.mu_photonAndNeutralHadronMiniIso = muon_photonAndNeutralHadronMiniIso[i];
          tmp_muCand.mu_d0 = muon_d0[i];
          tmp_muCand.mu_dZ = muon_dZ[i];
          tmp_muCand.mu_ip3d = muon_ip3d[i];
          tmp_muCand.mu_ip3dSignificance = muon_ip3dSignificance[i];
          muCand.push_back( tmp_muCand );
        }


        //-------------------------------
        //1) Look for Zmm Candidate
        //-------------------------------

        //Find two muons with the highest pt
        TLorentzVector ZCandidate(0,0,0,0);
        int ZMuIndex1 = -1;
        int ZMuIndex2 = -1;
        double bestDimuonPt = -1;
        std::vector< MuonCandidate > muSelectedCand;
        MuonCandidate bestMuCand[2];

        if( muCand.size() > 1 )
        {
          for ( size_t i = 0; i < muCand.size(); i++ )
          {
            for ( size_t j = i+1; j < muCand.size(); j++ )
            {
              MuonCandidate mu1 = muCand[i];
              MuonCandidate mu2 = muCand[j];
              //need dimuon mass between [76, 106] GeV
              double dimuonMass = (mu1.muon + mu2.muon).M();
              if( mu1.muon.Pt() + mu2.muon.Pt() > bestDimuonPt )
              {
                bestDimuonPt = mu1.muon.Pt() + mu2.muon.Pt();
                ZCandidate = mu1.muon + mu2.muon;
                if ( mu1.muon.Pt() >= mu2.muon.Pt() )
                {
                  if ( _debug ) std::cout << "assign muon candidate, mu1Pt > mu2Pt" << std::endl;
                  bestMuCand[0] = mu1;
                  bestMuCand[1] = mu2;
                  ZMuIndex1 = mu1.Index;
                  ZMuIndex2 = mu2.Index;
                }
                else
                {
                  if ( _debug ) std::cout << "assign muon candidate, mu2Pt > mu1Pt" << std::endl;
                  bestMuCand[0] = mu2;
                  bestMuCand[1] = mu1;
                  ZMuIndex1 = mu2.Index;
                  ZMuIndex2 = mu1.Index;
                }
              }//best pt if
            }
          }
          //---------------------------------------
          //just use this container for convenience
          //to parse the data into TTree
          //---------------------------------------
          muSelectedCand.push_back(bestMuCand[0]);
          muSelectedCand.push_back(bestMuCand[1]);

          mumuTree->f_mu_type[0]     = 13 * -1 * bestMuCand[0].muonCharge;
          mumuTree->f_mu_pt[0]       = bestMuCand[0].muon.Pt();
          mumuTree->f_mu_eta[0]      = bestMuCand[0].muon.Eta();
          mumuTree->f_mu_phi[0]      = bestMuCand[0].muon.Phi();
          mumuTree->f_mu_pass_sel[0] = 1 + 2 * bestMuCand[0].isTightMuon;
          mumuTree->f_mu_chargedMiniIso[0] = bestMuCand[0].mu_chargedMiniIso;
          mumuTree->f_mu_photonAndNeutralHadronMiniIso[0] = bestMuCand[0].mu_photonAndNeutralHadronMiniIso;
          mumuTree->f_mu_d0[0] = bestMuCand[0].mu_d0;
          mumuTree->f_mu_dZ[0] = bestMuCand[0].mu_dZ;
          mumuTree->f_mu_ip3d[0] = bestMuCand[0].mu_ip3d;
          mumuTree->f_mu_ip3dSignificance[0] = bestMuCand[0].mu_ip3dSignificance;
          //
          mumuTree->f_mu_type[1]     = 13 * -1 * bestMuCand[1].muonCharge;
          mumuTree->f_mu_pt[1]       = bestMuCand[1].muon.Pt();
          mumuTree->f_mu_eta[1]      = bestMuCand[1].muon.Eta();
          mumuTree->f_mu_phi[1]      = bestMuCand[1].muon.Phi();
          mumuTree->f_mu_pass_sel[1] = 1 + 2 * bestMuCand[1].isTightMuon;
          mumuTree->f_mu_chargedMiniIso[1] = bestMuCand[1].mu_chargedMiniIso;
          mumuTree->f_mu_photonAndNeutralHadronMiniIso[1] = bestMuCand[1].mu_photonAndNeutralHadronMiniIso;
          mumuTree->f_mu_d0[1] = bestMuCand[1].mu_d0;
          mumuTree->f_mu_dZ[1] = bestMuCand[1].mu_dZ;
          mumuTree->f_mu_ip3d[1] = bestMuCand[1].mu_ip3d;
          mumuTree->f_mu_ip3dSignificance[1] = bestMuCand[1].mu_ip3dSignificance;

          //met
          mumuTree->f_met_pt = metType1Pt;
          mumuTree->f_met_phi = metType1Phi;
          //for MC apply lepton eff scale factor
          if (!isData )
          {
            //if ( matchesGenMuon(lep1Eta,lep1Phi)) leptonEffSF *=  helper->getVetoMuonScaleFactor( lep1Pt, lep1Eta, true);
            //if ( matchesGenMuon(lep2Eta,lep2Phi)) leptonEffSF *=  helper->getVetoMuonScaleFactor( lep2Pt, lep2Eta, true);
          }
          //record Z candidate inf
          mumuTree->f_mumu_mass   = ZCandidate.M();
          mumuTree->f_mumu_pt     = ZCandidate.Pt();
          mumuTree->f_mumu_eta     = ZCandidate.Eta();
          mumuTree->f_mumu_phi     = ZCandidate.Phi();
          //if ( _debug ) std::cout << "[DEBUG]: dimuon mass-> " << dileptonMass << " dimuon pT->" << bestDimuonPt << std::endl;
          //}
        }//end if muCand.size() > 1
        //***********************
        //Fill Jet
        //***********************
        //NJetsFilled++;

        if( ZCandidate.Pt() > 0 )
        {
          mumuTree->tree_->Fill();
          mumuTree->f_i_evt++;
        }

      }
    }
    //end of event loop

    //cout << "Filled Total of " << NJetsFilled << " Jets\n";
    cout << "Writing output trees..." << endl;
    outFile->Write();
    outFile->Close();

}
