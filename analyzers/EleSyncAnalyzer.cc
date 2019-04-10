// This analyzer produces a ROOT tree containing electron ID info for each event.
// If an option number is given, only events with runNum == option will be filled.
// This was created for the electron synchronization exercise done on 21 Sept. 2016.

// C++ includes
#include <map>

// ROOT includes
#include "TLorentzVector.h"

// Local includes
#include "EleSyncAnalyzer.h"

float ELE_PT_CUT = 30;
float ELE_ETA_CUT = 2.5;
float ECAL_GAP_MIN = 1.4442;
float ECAL_GAP_MAX = 1.566;

class EleVarCollection {
    // Collection of electron ID variables.  
    public:
        // Contructor.
        // tag: label used for tree branch names
        // tree: TTree on which to make branches
        EleVarCollection(string tag_, TTree *tree) {
            tag = tag_; 
            resetEvent();
            elePtr = &ele;
            setBranches(tree); 
        }
        // Destructor
        ~EleVarCollection() {}

        // Member functions
        void resetEvent() {
            // Set all member variables to default values.
            // This function should be called once per event.
            ele.SetPtEtaPhiE(0,0,0,0);
            passID = false;
            conv = false;
            chHadIso = 0;
            neuHadIso = 0;
            gammaIso = 0;
            ea = 0;
            iso = 0;
            sieie = 0;
            dPhiIn = 0;
            dEtaIn = 0;
            hovere = 0;
            ooemoop = 0;
            d0 = 0;
            dz = 0;
            nMiss = 0;

            currentIndex = -1;
        }
            
        void setBranches(TTree *t) {
            // Add each variable as a tree branch
            t->Branch(("P"+tag).c_str(), "TLorentzVector", &elePtr);
            t->Branch(("passID"+tag).c_str(), &passID, "passID/O");
            t->Branch(("conv"+tag).c_str(), &conv, "conv/O");
            t->Branch(("chHadIso"+tag).c_str(), &chHadIso, "chHadIso/F");
            t->Branch(("neuHadIso"+tag).c_str(), &neuHadIso, "neuHadIso/F");
            t->Branch(("gammaIso"+tag).c_str(), &gammaIso, "gammaIso/F");
            t->Branch(("ea"+tag).c_str(), &ea, "ea/F");
            t->Branch(("iso"+tag).c_str(), &iso, "iso/F");
            t->Branch(("sieie"+tag).c_str(), &sieie, "sieie/F");
            t->Branch(("dPhiIn"+tag).c_str(), &dPhiIn, "dPhiIn/F");
            t->Branch(("dEtaIn"+tag).c_str(), &dEtaIn, "dEtaIn/F");
            t->Branch(("hovere"+tag).c_str(), &hovere, "hovere/F");
            t->Branch(("thatEP"+tag).c_str(), &ooemoop, "thatEP/F");
            t->Branch(("d0"+tag).c_str(), &d0, "d0/F");
            t->Branch(("dz"+tag).c_str(), &dz, "dz/F");
            t->Branch(("nMiss"+tag).c_str(), &nMiss, "nMiss/I");
        }

        void populate(RazorAnalyzer *event, int index) {
            // Look up the electron information in the RazorAnalyzer object and 
            // use it to populate all data members.  
            // Returns immediately if called with a negative index number.
            // 
            //   index: index of the electron in the RazorAnalyzer electron collection
            if ( index < 0 ) {
                return;
            }
            ele.SetPtEtaPhiE( event->elePt[index], event->eleEta[index], 
                    event->elePhi[index], event->eleE[index] );
            passID = event->isEGammaPOGMediumElectron( index );
            conv = !event->ele_PassConvVeto[index];
            chHadIso = event->ele_chargedIso[index];
            neuHadIso = event->ele_neutralHadIso[index];
            gammaIso = event->ele_photonIso[index];
            ea = event->GetElectronEffectiveAreaMean( index );
            iso = ( event->ele_chargedIso[index] + fmax( 0.0, event->ele_photonIso[index] + 
                    event->ele_neutralHadIso[index] - ea*event->fixedGridRhoFastjetAll ) );
            sieie = event->eleFull5x5SigmaIetaIeta[index];
            dPhiIn = event->ele_dPhi[index];
            dEtaIn = event->ele_dEta[index];
            hovere = event->ele_HoverE[index];
            ooemoop = event->ele_OneOverEminusOneOverP[index];
            d0 = event->ele_d0[index];
            dz = event->ele_dZ[index];
            nMiss = event->ele_MissHits[index];

            currentIndex = index;
        }

        // Variables
        string tag;
        int currentIndex;
        TLorentzVector ele;
        TLorentzVector *elePtr;
        bool passID, conv;
        float chHadIso, neuHadIso, gammaIso, ea, iso;
        float sieie, dPhiIn, dEtaIn, hovere, ooemoop, d0, dz;
        int nMiss;
};

bool passGoodLumi( int lumi ) {
    // For Run 278345, electrons
    if ( ( lumi >= 84 && lumi <= 125 ) ||
         ( lumi >= 127 && lumi <= 138 ) ||
         ( lumi >= 140 && lumi <= 198 ) ) return true;
    return false;
}

void EleSyncAnalyzer::Analyze(bool isData, int option, string outputFileName, string label)
{
    if (fChain == 0) return;

    // Disable all unused branches to maximize speed
    fChain->SetBranchStatus("*",0);
    EnableEventInfo();
    EnableElectrons();

    // Keep track of which files we filled events from
    std::map<std::string, bool> fileMap;
    
    Long64_t nentries = fChain->GetEntriesFast();

    if ( outputFileName == "" ) {
        if ( option >= 0 ) {
            outputFileName = "EventList"+std::to_string(option)+".root";
        }
        else {
            outputFileName = "EventList.root";
        }
    }

    // Set up output tree
    TFile *outFile = new TFile(outputFileName.c_str(), "RECREATE");
    TTree *outTree = new TTree("Events", "Events");
    int nEle = 0;
    outTree->Branch("runNum", &runNum, "runNum/I");
    outTree->Branch("lumiNum", &lumiNum, "lumiNum/I");
    outTree->Branch("eventNum", &eventNum, "eventNum/I");
    outTree->Branch("rho", &fixedGridRhoFastjetAll, "rho/F");
    outTree->Branch("nEle", &nEle, "nEle/I");

    // All electron variables are handled by the EleVarCollection class
    EleVarCollection *vars_e1 = new EleVarCollection("_e1", outTree);
    EleVarCollection *vars_e2 = new EleVarCollection("_e2", outTree);

    std::cout << "Writing event list tree to " << outputFileName << std::endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        if(jentry % 100000 == 0) cout << "Processing entry " << jentry << endl;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        fChain->GetEntry(jentry);  

        // Apply good lumi selection
        if ( !passGoodLumi( lumiNum ) ) continue;

        // Reset event
        nEle = 0;
        vars_e1->resetEvent();
        vars_e2->resetEvent();

        // Select electrons
        for ( int i = 0; i < nElectrons; i++ ) {

            // Cuts here
            if ( elePt[i] < ELE_PT_CUT ) continue;
            if ( fabs(eleEta[i]) > ELE_ETA_CUT ) continue;
            if ( fabs(eleEta[i]) > ECAL_GAP_MIN && fabs(eleEta[i]) < ECAL_GAP_MAX ) continue;
            if ( fabs(eleEta_SC[i]) > ECAL_GAP_MIN && fabs(eleEta_SC[i]) < ECAL_GAP_MAX ) continue;
            if ( !isEGammaPOGMediumElectron(i) ) continue;
            // End cuts

            nEle++;

            // Find the two highest-pt electrons
            // and store their information in the tree variables
            if ( elePt[i] > vars_e1->ele.Pt() ) {
                vars_e2->populate( this, vars_e1->currentIndex );
                vars_e1->populate( this, i );
            }
            else if ( elePt[i] > vars_e2->ele.Pt() ) {
                vars_e2->populate( this, i );
            }
        }

        // Apply extra cuts according to the option label
        if ( label == "trigger" ) {
            if ( !HLTDecision[35] ) continue;
        }
        else if ( label == "ee" ) {
            if ( !HLTDecision[35] || nEle < 2 ) continue;
        }
        else if ( label == "zee" ) {
            float mll = (vars_e1->ele + vars_e2->ele).M();
            if ( !HLTDecision[35] || nEle < 2 || mll < 60 || mll > 120 ) continue;
        }

        if( option < 0 || runNum == static_cast<unsigned int>(option) ) {
            outTree->Fill();
            std::string fileName = fChain->GetCurrentFile()->GetName();
            fileMap[fileName] = true;
        }
    } // End event loop

    outTree->Write();
    outFile->Close();

    std::cout << "Filled events from these files:" << std::endl;
    for ( auto str : fileMap ) {
        std::cout << str.first << std::endl;
    }
}
