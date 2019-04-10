#include "FullRazorDM.h"
#include "RazorHelper.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "BTagCalibrationStandalone.h"

//C++ includes
#include <sys/stat.h>
#include <assert.h>
#include <sstream>

//ROOT includes
#include "TH1F.h"
#include "TH2D.h"

using namespace std;

const int NUM_PDF_WEIGHTS = 60;

// Collection to hold main analysis variables
class RazorVarCollection {

    public:
        // Constructor
        RazorVarCollection(string tag_) { tag = tag_; resetVars(); }
        // Destructor
        ~RazorVarCollection() { }

        // Member functions 
        void resetVars() { //call for each event
            MR = -1; Rsq = -1; dPhiRazor = -9;
            leadingJetPt = -1; subleadingJetPt = -1; 
            leadingTightMuPt = -1; leadingTightElePt = -1;
            mT = -1; mTLoose = -1;
            nSelectedJets = 0; nBTaggedJets = 0; nJets80 = 0;
            nVetoMuons = 0; nTightMuons = 0; nVetoElectrons = 0; nTightElectrons = 0;
            box = RazorAnalyzer::NONE;

            MetXCorr = 0; MetYCorr = 0;
            metOverCaloMet = 0;
            leadingTightMu = TLorentzVector(); leadingTightEle = TLorentzVector();
            GoodJets = vector<TLorentzVector>();
            GoodLeptons = vector<TLorentzVector>();
        }   
        void setBranches(TTree *t) { //add branches to tree
            string conn = "_";
            if (tag == "") { conn = ""; } // remove underscore if not needed
            t->Branch(("MR"+conn+tag).c_str(), &MR, ("MR"+conn+tag+"/F").c_str());
            t->Branch(("Rsq"+conn+tag).c_str(), &Rsq, ("Rsq"+conn+tag+"/F").c_str());
            t->Branch(("metOverCaloMet"+conn+tag).c_str(), &metOverCaloMet, ("metOverCaloMet"+conn+tag+"/F").c_str());
            t->Branch(("dPhiRazor"+conn+tag).c_str(), &dPhiRazor, ("dPhiRazor"+conn+tag+"/F").c_str());
            t->Branch(("leadingJetPt"+conn+tag).c_str(), &leadingJetPt, 
                    ("leadingJetPt"+conn+tag+"/F").c_str());
            t->Branch(("subleadingJetPt"+conn+tag).c_str(), &subleadingJetPt, 
                    ("subleadingJetPt"+conn+tag+"/F").c_str());
            t->Branch(("leadingTightMuPt"+conn+tag).c_str(), &leadingTightMuPt, 
                    ("leadingTightMuPt"+conn+tag+"/F").c_str());
            t->Branch(("leadingTightElePt"+conn+tag).c_str(), &leadingTightElePt, 
                    ("leadingTightElePt"+conn+tag+"/F").c_str());
            t->Branch(("mT"+conn+tag).c_str(), &mT, ("mT"+conn+tag+"/F").c_str());
            t->Branch(("mTLoose"+conn+tag).c_str(), &mTLoose, ("mTLoose"+conn+tag+"/F").c_str());
            t->Branch(("nSelectedJets"+conn+tag).c_str(), &nSelectedJets, 
                    ("nSelectedJets"+conn+tag+"/I").c_str());
            t->Branch(("nBTaggedJets"+conn+tag).c_str(), &nBTaggedJets, 
                    ("nBTaggedJets"+conn+tag+"/I").c_str());
            t->Branch(("nJets80"+conn+tag).c_str(), &nJets80, ("nJets80"+conn+tag+"/I").c_str());
            t->Branch(("box"+conn+tag).c_str(), &box, ("box"+conn+tag+"/I").c_str());
            if (tag == "") {
                t->Branch(("nVetoMuons"+conn+tag).c_str(), &nVetoMuons, 
                        ("nVetoMuons"+conn+tag+"/I").c_str());
                t->Branch(("nTightMuons"+conn+tag).c_str(), &nTightMuons, 
                        ("nTightMuons"+conn+tag+"/I").c_str());
                t->Branch(("nVetoElectrons"+conn+tag).c_str(), &nVetoElectrons, 
                        ("nVetoElectrons"+conn+tag+"/I").c_str());
                t->Branch(("nTightElectrons"+conn+tag).c_str(), &nTightElectrons, 
                        ("nTightElectrons"+conn+tag+"/I").c_str());
            }
        }

        // List of variables
        float MR,Rsq,dPhiRazor,leadingJetPt,subleadingJetPt,leadingTightMuPt,leadingTightElePt,mT,mTLoose;
        int nSelectedJets,nBTaggedJets,nJets80;
        int nVetoMuons, nTightMuons, nVetoElectrons, nTightElectrons;
        float metOverCaloMet;
        RazorAnalyzer::RazorBox box;
        // Non-tree variables
        float MetXCorr, MetYCorr;
        TLorentzVector leadingTightMu, leadingTightEle;
        vector<TLorentzVector> GoodJets, GoodLeptons;
        string tag;
};

void FullRazorDM::Analyze(bool isData, int option, string outFileName, string label)
{
    /////////////////////////////////
    //Basic setup
    /////////////////////////////////

    cout << "Initializing..." << endl;
    bool isFastsimSMS = (option == 1 || option == 11);

    //Output file
    if (outFileName.empty()){
        cout << "FullRazorDM: Output filename not specified!" << endl << "Using default output name FullRazorDM.root" << endl;
        outFileName = "FullRazorDM.root";
    }
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");

    //Output tree
    TTree *razorTree = new TTree("RazorDM", "Info on selected razor DM events");

    //For signal samples, create one output file and tree per signal mass point
    map<pair<int,int>, TFile*> smsFiles;
    map<pair<int,int>, TTree*> smsTrees;
    map<pair<int,int>, TH1F*> smsNEvents;
    map<pair<int,int>, TH1F*> smsSumWeights;
    map<pair<int,int>, TH1F*> smsSumScaleWeights;
    map<pair<int,int>, TH1F*> smsSumPdfWeights;

    //Histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 0.5, 1.5);
    TH1F *SumTopPtWeights = new TH1F("SumTopPtWeights", "SumTopPtWeights", 1, 0.5, 1.5);
    TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
    TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
    TH1F *SumPdfWeights = new TH1F("SumPdfWeights", "SumPdfWeights", NUM_PDF_WEIGHTS, -0.5, NUM_PDF_WEIGHTS-0.5);

    //Initialize helper
    RazorHelper *helper = 0;
    string analysisTag = "Razor2016_80X";
    if ( label != "") analysisTag = label;
    helper = new RazorHelper(analysisTag, isData, isFastsimSMS);

    // Get jet corrector
    std::vector<FactorizedJetCorrector*> JetCorrector = helper->getJetCorrector();
    std::vector<std::pair<int,int> > JetCorrectorIOV = helper->getJetCorrectionsIOV();

    /////////////////////////////////
    //Tree Initialization
    /////////////////////////////////

    //To hold main variables
    map<string,RazorVarCollection*> mainVars;
    vector<string> varCollectionNames;
    if(isData) {
        varCollectionNames = { "" };
    }
    else {
        varCollectionNames = { "", "JESUp", "JESDown", "MESUp", "MESDown", 
            "EESUp", "EESDown", "JERUp", "JERDown" };
    }
    for (auto &str : varCollectionNames) {
        mainVars[str] = new RazorVarCollection(str);
    }

    //Basic tree variables
    int nVtx, NPU; 
    float weight = 1.0;
    float btagCorrFactor, topPtWeight;
    //For signal ISR systematic uncertainty
    float ISRSystWeightUp, ISRSystWeightDown;
    //For pileup systematic uncertainty
    float pileupWeight, pileupWeightUp, pileupWeightDown;
    //For lepton efficiency scale factor uncertainty
    float sf_muonEffUp, sf_muonEffDown;
    float sf_eleEffUp, sf_eleEffDown;
    float sf_muonTrigUp, sf_muonTrigDown;
    float sf_eleTrigUp, sf_eleTrigDown;
    float sf_tauEffUp, sf_tauEffDown;
    float sf_vetoMuonEffUp, sf_vetoMuonEffDown;
    float sf_vetoEleEffUp, sf_vetoEleEffDown;
    //For btag scale factor uncertainty
    float sf_btagUp, sf_btagDown;
    float sf_bmistagUp, sf_bmistagDown;
    //For scale variation uncertainties
    float sf_facScaleUp, sf_facScaleDown;
    float sf_renScaleUp, sf_renScaleDown;
    float sf_facRenScaleUp, sf_facRenScaleDown;
    //For Fastsim scale factor uncertainties
    float sf_muonEffFastsimSFUp, sf_muonEffFastsimSFDown;
    float sf_eleEffFastsimSFUp, sf_eleEffFastsimSFDown;
    float sf_vetoMuonEffFastsimSFUp, sf_vetoMuonEffFastsimSFDown;
    float sf_vetoEleEffFastsimSFUp, sf_vetoEleEffFastsimSFDown;
    float sf_btagFastsimSFUp, sf_btagFastsimSFDown;
    //For jet uncertainties
    int nLooseTaus;
    float met, HT;
    float mjj_leadingJets, mjj_hemispheres;
    float leadingGenLeptonPt;
    float leadingGenLeptonEta;
    int   leadingGenLeptonType;
    float subLeadingGenLeptonEta;
    float subLeadingGenLeptonPt;
    int   subLeadingGenLeptonType;
    int NGenBJets;
    float genHT;
    int NISRJets;
    //SMS parameters 
    int mGluino, mLSP;
    int nCharginoFromGluino, ntFromGluino;
    bool Flag_hasEcalGainSwitch;

    //Set branches
    for (auto &vars : mainVars) {
        vars.second->setBranches(razorTree);
    }

    razorTree->Branch("nVtx", &nVtx, "nVtx/I");
    razorTree->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
    razorTree->Branch("HT", &HT, "HT/F");
    razorTree->Branch("met", &met, "met/F");
    razorTree->Branch("mjj_leadingJets", &mjj_leadingJets, "mjj_leadingJets/F");
    razorTree->Branch("mjj_hemispheres", &mjj_hemispheres, "mjj_hemispheres/F");
    razorTree->Branch("HLTDecision", &HLTDecision, "HLTDecision[300]/O");
    //MET filters
    razorTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
    razorTree->Branch("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, "Flag_HBHEIsoNoiseFilter/O");
    razorTree->Branch("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter, "Flag_badChargedCandidateFilter/O");
    razorTree->Branch("Flag_badMuonFilter", &Flag_badMuonFilter, "Flag_badMuonFilter/O");
    razorTree->Branch("Flag_badGlobalMuonFilter", &Flag_badGlobalMuonFilter, "Flag_badGlobalMuonFilter/O");
    razorTree->Branch("Flag_duplicateMuonFilter", &Flag_duplicateMuonFilter, "Flag_duplicateMuonFilter/O");
    razorTree->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, "Flag_CSCTightHaloFilter/O");
    razorTree->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, "Flag_hcalLaserEventFilter/O");
    razorTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/O");
    razorTree->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/O");
    razorTree->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, "Flag_trackingFailureFilter/O");
    razorTree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/O");
    razorTree->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, "Flag_ecalLaserCorrFilter/O");
    razorTree->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters, "Flag_trkPOGFilters/O");
    razorTree->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, "Flag_trkPOG_manystripclus53X/O");
    razorTree->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, "Flag_trkPOG_toomanystripclus53X/O");
    razorTree->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, "Flag_trkPOG_logErrorTooManyClusters/O");
    razorTree->Branch("Flag_METFilters", &Flag_METFilters, "Flag_METFilters/O");
    razorTree->Branch("Flag_hasEcalGainSwitch", &Flag_hasEcalGainSwitch, "Flag_hasEcalGainSwitch/O");

    if (!isData) {    
        razorTree->Branch("genWeight", &genWeight, "genWeight/F");
        razorTree->Branch("weight", &weight, "weight/F");
        razorTree->Branch("ISRSystWeightUp", &ISRSystWeightUp, "ISRSystWeightUp/F");
        razorTree->Branch("ISRSystWeightDown", &ISRSystWeightDown, "ISRSystWeightDown/F");
        razorTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/F");
        razorTree->Branch("pileupWeightUp", &pileupWeightUp, "pileupWeightUp/F");
        razorTree->Branch("pileupWeightDown", &pileupWeightDown, "pileupWeightDown/F");
        razorTree->Branch("btagCorrFactor", &btagCorrFactor, "btagCorrFactor/F");
        razorTree->Branch("topPtWeight", &topPtWeight, "topPtWeight/F");
        razorTree->Branch("NPU", &NPU, "NPU/I");
        razorTree->Branch("leadingGenLeptonPt", &leadingGenLeptonPt, "leadingGenLeptonPt/F");
        razorTree->Branch("leadingGenLeptonEta", &leadingGenLeptonEta, "leadingGenLeptonEta/F");
        razorTree->Branch("leadingGenLeptonType", &leadingGenLeptonType, "leadingGenLeptonType/I");
        razorTree->Branch("subLeadingGenLeptonPt", &subLeadingGenLeptonPt, "subLeadingGenLeptonPt/F");
        razorTree->Branch("subLeadingGenLeptonEta", &subLeadingGenLeptonEta, "subLeadingGenLeptonEta/F");
        razorTree->Branch("subLeadingGenLeptonType", &subLeadingGenLeptonType, "subLeadingGenLeptonType/I");
        razorTree->Branch("NGenBJets", &NGenBJets, "NGenBJets/I");
        razorTree->Branch("genHT", &genHT, "genHT/F");
        razorTree->Branch("NISRJets", &NISRJets, "NISRJets/I");
        razorTree->Branch("sf_muonEffUp", &sf_muonEffUp, "sf_muonEffUp/F");
        razorTree->Branch("sf_muonEffDown", &sf_muonEffDown, "sf_muonEffDown/F");
        razorTree->Branch("sf_vetoMuonEffUp", &sf_vetoMuonEffUp, "sf_vetoMuonEffUp/F");
        razorTree->Branch("sf_vetoMuonEffDown", &sf_vetoMuonEffDown, "sf_vetoMuonEffDown/F");
        razorTree->Branch("sf_eleEffUp", &sf_eleEffUp, "sf_eleEffUp/F");
        razorTree->Branch("sf_eleEffDown", &sf_eleEffDown, "sf_eleEffDown/F");
        razorTree->Branch("sf_vetoEleEffUp", &sf_vetoEleEffUp, "sf_vetoEleEffUp/F");
        razorTree->Branch("sf_vetoEleEffDown", &sf_vetoEleEffDown, "sf_vetoEleEffDown/F");
        razorTree->Branch("sf_tauEffUp", &sf_tauEffUp, "sf_tauEffUp/F");
        razorTree->Branch("sf_tauEffDown", &sf_tauEffDown, "sf_tauEffDown/F");
        razorTree->Branch("sf_muonTrigUp", &sf_muonTrigUp, "sf_muonTrigUp/F");
        razorTree->Branch("sf_muonTrigDown", &sf_muonTrigDown, "sf_muonTrigDown/F");
        razorTree->Branch("sf_eleTrigUp", &sf_eleTrigUp, "sf_eleTrigUp/F");
        razorTree->Branch("sf_eleTrigDown", &sf_eleTrigDown, "sf_eleTrigDown/F");
        razorTree->Branch("sf_btagUp", &sf_btagUp, "sf_btagUp/F");
        razorTree->Branch("sf_btagDown", &sf_btagDown, "sf_btagDown/F");
        razorTree->Branch("sf_bmistagUp", &sf_bmistagUp, "sf_bmistagUp/F");
        razorTree->Branch("sf_bmistagDown", &sf_bmistagDown, "sf_bmistagDown/F");
        razorTree->Branch("sf_muonEffFastsimSFUp", &sf_muonEffFastsimSFUp, "sf_muonEffFastsimSFUp/F");
        razorTree->Branch("sf_muonEffFastsimSFDown", &sf_muonEffFastsimSFDown, "sf_muonEffFastsimSFDown/F");
        razorTree->Branch("sf_eleEffFastsimSFUp", &sf_eleEffFastsimSFUp, "sf_eleEffFastsimSFUp/F");
        razorTree->Branch("sf_eleEffFastsimSFDown", &sf_eleEffFastsimSFDown, "sf_eleEffFastsimSFDown/F");
        razorTree->Branch("sf_btagFastsimSFUp", &sf_btagFastsimSFUp, "sf_btagFastsimSFUp/F");
        razorTree->Branch("sf_btagFastsimSFDown", &sf_btagFastsimSFDown, "sf_btagFastsimSFDown/F");
        razorTree->Branch("sf_vetoMuonEffFastsimSFUp", &sf_vetoMuonEffFastsimSFUp, "sf_vetoMuonEffFastsimSFUp/F");
        razorTree->Branch("sf_vetoMuonEffFastsimSFDown", &sf_vetoMuonEffFastsimSFDown, "sf_vetoMuonEffFastsimSFDown/F");
        razorTree->Branch("sf_vetoEleEffFastsimSFUp", &sf_vetoEleEffFastsimSFUp, "sf_vetoEleEffFastsimSFUp/F");
        razorTree->Branch("sf_vetoEleEffFastsimSFDown", &sf_vetoEleEffFastsimSFDown, "sf_vetoEleEffFastsimSFDown/F");
        razorTree->Branch("sf_facScaleUp", &sf_facScaleUp, "sf_facScaleUp/F");
        razorTree->Branch("sf_facScaleDown", &sf_facScaleDown, "sf_facScaleDown/F");
        razorTree->Branch("sf_renScaleUp", &sf_renScaleUp, "sf_renScaleUp/F");
        razorTree->Branch("sf_renScaleDown", &sf_renScaleDown, "sf_renScaleDown/F");
        razorTree->Branch("sf_facRenScaleUp", &sf_facRenScaleUp, "sf_facRenScaleUp/F");
        razorTree->Branch("sf_facRenScaleDown", &sf_facRenScaleDown, "sf_facRenScaleDown/F");
        razorTree->Branch("pdfWeights", "std::vector<float>",&pdfWeights); //get PDF weights directly from RazorEvents
        if(isFastsimSMS){
            razorTree->Branch("mGluino", &mGluino, "mGluino/I");
            razorTree->Branch("mLSP", &mLSP, "mLSP/I");
            razorTree->Branch("nCharginoFromGluino", &nCharginoFromGluino, "nCharginoFromGluino/I");
            razorTree->Branch("ntFromGluino", &ntFromGluino, "ntFromGluino/I");
        }
    } 
    else {
        razorTree->Branch("run", &runNum, "run/i");
        razorTree->Branch("lumi", &lumiNum, "lumi/i");
        razorTree->Branch("event", &eventNum, "event/i");
    }

    /////////////////////////////////
    //Event loop
    /////////////////////////////////

    if (fChain == 0) return;
    cout << "Total Events: " << fChain->GetEntries() << "\n";
    Long64_t nbytes = 0, nb = 0;    
    Long64_t nentries = fChain->GetEntriesFast();
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        /////////////////////////////////
        //Begin event
        /////////////////////////////////

        //Initialize
        if (jentry % 10000 == 0) cout << "Processing entry " << jentry << endl;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;

        nb = fChain->GetEntry(jentry);

        //Reset tree variables
        for (auto &vars : mainVars) {
            vars.second->resetVars();
        }
        nVtx = nPV;
        mjj_leadingJets = -1;
        mjj_hemispheres = -1;
        weight = genWeight;
        nLooseTaus = 0;
        if(!isData){
            NPU = -1;
            leadingGenLeptonPt = -9;
            leadingGenLeptonEta = -9;
            leadingGenLeptonType = 0;
            subLeadingGenLeptonPt = -9;
            subLeadingGenLeptonEta = -9;
            subLeadingGenLeptonType = 0;
            NGenBJets = 0;
            genHT = 0;
            NISRJets = 0;
            ISRSystWeightUp = 1.0;
            ISRSystWeightDown = 1.0;
            pileupWeight = 1.0;
            pileupWeightUp = 1.0;
            pileupWeightDown = 1.0;
            btagCorrFactor = 1.0;
            topPtWeight = 1.0;
            sf_muonEffUp = 1.0;
            sf_muonEffDown = 1.0;
            sf_vetoMuonEffUp = 1.0;
            sf_vetoMuonEffDown = 1.0;
            sf_eleEffUp = 1.0;
            sf_eleEffDown = 1.0;
            sf_vetoEleEffUp = 1.0;
            sf_vetoEleEffDown = 1.0;
            sf_muonTrigUp = 1.0;
            sf_muonTrigDown = 1.0;
            sf_eleTrigUp = 1.0;
            sf_eleTrigDown = 1.0;
            sf_tauEffUp = 1.0;
            sf_tauEffDown = 1.0;
            sf_btagUp = 1.0;
            sf_btagDown = 1.0;
            sf_bmistagUp = 1.0;
            sf_bmistagDown = 1.0;
            sf_facScaleUp = 1.0;
            sf_facScaleDown = 1.0;
            sf_renScaleUp = 1.0;
            sf_renScaleDown = 1.0;
            sf_facRenScaleUp = 1.0;
            sf_facRenScaleDown = 1.0;
            sf_muonEffFastsimSFUp = 1.0;
            sf_muonEffFastsimSFDown = 1.0;
            sf_eleEffFastsimSFUp = 1.0;
            sf_eleEffFastsimSFDown = 1.0;
            sf_vetoMuonEffFastsimSFUp = 1.0;
            sf_vetoMuonEffFastsimSFDown = 1.0;
            sf_vetoEleEffFastsimSFUp = 1.0;
            sf_vetoEleEffFastsimSFDown = 1.0;
            sf_btagFastsimSFUp = 1.0;
            sf_btagFastsimSFDown = 1.0;
            if(isFastsimSMS){
                mGluino = -1;
                mLSP = -1;
                nCharginoFromGluino = 0;
                ntFromGluino = 0;
            }
        }

        /////////////////////////////////
        //MC particles
        /////////////////////////////////
        genHT = getGenHT();
        NISRJets = getNISR( JetCorrector, JetCorrectorIOV );
        float ptTop = -1;
        float ptAntitop = -1;
        if(!isData) {
            for(int j = 0; j < nGenParticle; j++){

                //electron or muon
                if ( gParticleStatus[j] == 1 &&
                        (abs(gParticleId[j]) == 11 || abs(gParticleId[j]) == 13)
                        && (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23 || abs(gParticleMotherId[j]) == 15)
                   ) {
                    if (gParticlePt[j] > leadingGenLeptonPt) {
                        subLeadingGenLeptonPt = leadingGenLeptonPt;
                        subLeadingGenLeptonEta = leadingGenLeptonEta;
                        subLeadingGenLeptonType = leadingGenLeptonType;
                        leadingGenLeptonPt = gParticlePt[j];
                        leadingGenLeptonEta = gParticleEta[j];
                        leadingGenLeptonType = gParticleId[j];
                    } else if ( gParticlePt[j] > subLeadingGenLeptonPt) {
                        subLeadingGenLeptonPt = gParticlePt[j];
                        subLeadingGenLeptonEta = gParticleEta[j];
                        subLeadingGenLeptonType = gParticleId[j];
                    }
                }

                //hadronic tau
                if ( gParticleStatus[j] == 2 && abs(gParticleId[j]) == 15
                        && (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23)
                   ) {
                    bool isLeptonicTau = false;
                    for(int k = 0; k < nGenParticle; k++){
                        if ( (abs(gParticleId[k]) == 11 || abs(gParticleId[k]) == 13) && gParticleMotherIndex[k] == j) {
                            isLeptonicTau = true;
                            break;
                        }
                    }
                    if (isLeptonicTau) continue;

                    if (gParticlePt[j] > leadingGenLeptonPt) {
                        subLeadingGenLeptonPt = leadingGenLeptonPt;
                        subLeadingGenLeptonEta = leadingGenLeptonEta;
                        subLeadingGenLeptonType = leadingGenLeptonType;
                        leadingGenLeptonPt = gParticlePt[j];
                        leadingGenLeptonEta = gParticleEta[j];
                        leadingGenLeptonType = gParticleId[j];
                    } else if ( gParticlePt[j] > subLeadingGenLeptonPt) {
                        subLeadingGenLeptonPt = gParticlePt[j];
                        subLeadingGenLeptonEta = gParticleEta[j];
                        subLeadingGenLeptonType = gParticleId[j];
                    }
                }

                //top
                if ( gParticleStatus[j] == 22 && gParticleId[j] == 6  && ptTop < 0 ) {
                    ptTop = gParticlePt[j];
                }
                //antitop
                if ( gParticleStatus[j] == 22 && gParticleId[j] == -6 && ptAntitop < 0 ) {
                    ptAntitop = gParticlePt[j];
                }

            } //loop over gen particles
            // get top pt weight
            if ( ptTop > 0 && ptAntitop > 0 ) {
                topPtWeight = helper->getTopPtWeight( ptTop, ptAntitop );
            }
        } //if !isData

        if(isFastsimSMS){

            //Count gluino to b quark decays and gluino to top quark decays
            int tmp_nbFromGluino = 0;
            int tmp_ntopFromGluino = 0;
            int tmp_nCharginoFromGluino = 0;
            for(int j = 0; j < nGenParticle; j++){
                //cout << "Particle " << j << " : " << gParticleId[j] << " " << gParticleStatus[j] << " : " << gParticlePt[j] << " " << gParticleEta[j] << " " << gParticlePhi[j] << " : " << gParticleMotherIndex[j] << " " << gParticleMotherId[j] << "\n";

                if ( abs(gParticleId[j]) == 5 && gParticleMotherIndex[j] >= 0 
                        && gParticleId[gParticleMotherIndex[j]] == 1000021 
                        && gParticleStatus[gParticleMotherIndex[j]] == 22
                   ) tmp_nbFromGluino++;

                if ( abs(gParticleId[j]) == 6 && gParticleMotherIndex[j] >= 0 
                        && gParticleId[gParticleMotherIndex[j]] == 1000021 
                        && gParticleStatus[gParticleMotherIndex[j]] == 22
                   ) tmp_ntopFromGluino++;														    

                if ( abs(gParticleId[j]) == 1000024 && gParticleMotherIndex[j] >= 0 
                        && gParticleId[gParticleMotherIndex[j]] == 1000021 
                        && gParticleStatus[gParticleMotherIndex[j]] == 22
                   ) tmp_nCharginoFromGluino++;														    
            }
            ntFromGluino = tmp_ntopFromGluino;
            nCharginoFromGluino = tmp_nCharginoFromGluino;	  


            //Get Gen level info for ISR systematics
            TLorentzVector *gluino1PreShowering = 0;
            TLorentzVector *gluino2PreShowering = 0;
            TLorentzVector *gluino1PostShowering = 0;
            TLorentzVector *gluino2PostShowering = 0;
            for(int j = 0; j < nGenParticle; j++){

                if (gParticleId[j] == 1000021 && gParticleStatus[j] == 22) {
                    if (!gluino1PreShowering) {
                        gluino1PreShowering = new TLorentzVector;
                        gluino1PreShowering->SetPtEtaPhiE(gParticlePt[j],gParticleEta[j],gParticlePhi[j],gParticleE[j]);
                    } else if (!gluino2PreShowering) {
                        gluino2PreShowering = new TLorentzVector;
                        gluino2PreShowering->SetPtEtaPhiE(gParticlePt[j],gParticleEta[j],gParticlePhi[j],gParticleE[j]);
                    } else {
                        cout << "Warning More than 2 status 22 gluinos\n";
                    }
                }

                if (gParticleId[j] == 1000021 && gParticleStatus[j] == 62) {
                    if (!gluino1PostShowering) {
                        gluino1PostShowering = new TLorentzVector;
                        gluino1PostShowering->SetPtEtaPhiE(gParticlePt[j],gParticleEta[j],gParticlePhi[j],gParticleE[j]);
                    } else if (!gluino2PostShowering) {
                        gluino2PostShowering = new TLorentzVector;
                        gluino2PostShowering->SetPtEtaPhiE(gParticlePt[j],gParticleEta[j],gParticlePhi[j],gParticleE[j]);
                    } else {
                        cout << "Warning More than 2 status 62 gluinos\n";
                    }
                }

                if (gluino1PreShowering && gluino2PreShowering && gluino1PostShowering && gluino2PostShowering) break;
            }

            if (gluino1PreShowering && gluino2PreShowering) {
                // cout << "PreShowering Gluino Pair System: " 
                // 	 << ((*gluino1PreShowering) + (*gluino2PreShowering)).Pt() << " "
                // 	 << ((*gluino1PreShowering) + (*gluino2PreShowering)).Eta() << " "
                // 	 << ((*gluino1PreShowering) + (*gluino2PreShowering)).Phi() << " "
                // 	 << ((*gluino1PreShowering) + (*gluino2PreShowering)).M() << " "
                // 	 << " \n";	    
            }
            if (gluino1PostShowering && gluino2PostShowering) {
                // cout << "PostShowering Gluino Pair System: " 
                // 	 << ((*gluino1PostShowering) + (*gluino2PostShowering)).Pt() << " "
                // 	 << ((*gluino1PostShowering) + (*gluino2PostShowering)).Eta() << " "
                // 	 << ((*gluino1PostShowering) + (*gluino2PostShowering)).Phi() << " "
                // 	 << ((*gluino1PostShowering) + (*gluino2PostShowering)).M() << " "
                // 	 << " \n";	 

                ISRSystWeightUp = 1.0; 
                if ( ((*gluino1PostShowering) + (*gluino2PostShowering)).Pt() > 400 && 
                        ((*gluino1PostShowering) + (*gluino2PostShowering)).Pt() <= 600 ) {
                    ISRSystWeightUp = 1.15;
                    ISRSystWeightDown = 0.85;
                } else if ( ((*gluino1PostShowering) + (*gluino2PostShowering)).Pt() > 600) {
                    ISRSystWeightUp = 1.30;
                    ISRSystWeightDown = 0.70;
                }		 
            }

            if (gluino1PreShowering) delete gluino1PreShowering;
            if (gluino2PreShowering) delete gluino2PreShowering;
            if (gluino1PostShowering) delete gluino1PostShowering;
            if (gluino2PostShowering) delete gluino2PostShowering;	 
        } // end if fastsim signals

        /////////////////////////////////
        //Trigger
        /////////////////////////////////

        bool passedDileptonTrigger = false;
        bool passedSingleLeptonTrigger = false;
        bool passedHadronicTrigger= false;

        vector<int> dileptonTriggerNums = helper->getDileptonTriggerNums();
        vector<int> singleLeptonTriggerNums = helper->getSingleLeptonTriggerNums();
        vector<int> hadronicTriggerNums = helper->getHadronicTriggerNums();
        for( unsigned int itrig = 0; itrig < dileptonTriggerNums.size(); itrig++ ) {
            if (HLTDecision[dileptonTriggerNums[itrig]]) {
                passedDileptonTrigger = true;
                break;
            }
        }
        for( unsigned int itrig = 0; itrig < singleLeptonTriggerNums.size(); itrig++ ) {
            if (HLTDecision[singleLeptonTriggerNums[itrig]]) {
                passedSingleLeptonTrigger = true;
                break;
            }
        }
        for( unsigned int itrig = 0; itrig < hadronicTriggerNums.size(); itrig++ ) {
            if (HLTDecision[hadronicTriggerNums[itrig]]) {
                passedHadronicTrigger = true;
                break;
            }
        }
        //ignore trigger for Fastsim, and for 80X MC
        if(isFastsimSMS || 
                ((analysisTag == "Razor2016_MoriondRereco" 
                  || analysisTag == "Razor2016_80X" 
                  || analysisTag == "Razor2016G_80X" 
                  || analysisTag == "Razor2016G_SUSYUnblind_80X") && !isData)
          ){
            passedDileptonTrigger = true;
            passedSingleLeptonTrigger = true;
            passedHadronicTrigger = true;
        }

        /////////////////////////////////
        //Pileup reweighting
        /////////////////////////////////

        pileupWeight = 1.0;
        if(!isData){
            //Get number of PU interactions
            for (int i = 0; i < nBunchXing; i++) {
                if (BunchXing[i] == 0) {
                    NPU = nPUmean[i];
                }
            }
            pileupWeight = helper->getPileupWeight(NPU);
            pileupWeightUp = helper->getPileupWeightUp(NPU) / pileupWeight;
            pileupWeightDown = helper->getPileupWeightDown(NPU) / pileupWeight;
        }










        /////////////////////////////////
        //Muon selection
        /////////////////////////////////

        //Scale factor
        float muonEffCorrFactor = 1.0;
        float muonTrigCorrFactor = 1.0;
        float vetoMuonEffCorrFactor = 1.0;
        //Cut parameters
        const float MUON_VETO_CUT = 5;
        const float MUON_TIGHT_CUT = 25;
        //Loop muons
        for (int i = 0; i < nMuons; i++){
            //TLorentzVector for this muon
            TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 
            //acceptance cut
            if (abs(muonEta[i]) > 2.4) continue;
            //lepton scale uncertainty
            if (!isData) {
                float muonPtUp = muonPt[i];
                float muonPtDown = muonPt[i];
                float muonEUp = muonE[i];
                float muonEDown = muonE[i];
                //get up/down muon momenta
                float sfUp = 1.0;
                float sfDown = 1.0;
                if (muonPt[i] < 100) {
                    sfUp = 1.002;
                    sfDown = 0.998;
                }
                else {
                    sfUp = 1.05;
                    sfDown = 0.95;
                }
                muonPtUp *= sfUp;
                muonPtDown *= sfDown;
                muonEUp *= sfUp;
                muonEDown *= sfDown;
                TLorentzVector thisMuonUp = makeTLorentzVector(muonPtUp, muonEta[i], muonPhi[i], muonEUp); 
                TLorentzVector thisMuonDown = makeTLorentzVector(muonPtDown, muonEta[i], muonPhi[i], muonEDown); 
                //Veto selection
                if (isVetoMuon(i)) {
                    if (muonPtUp > MUON_VETO_CUT) {
                        mainVars["MESUp"]->nVetoMuons++;
                        mainVars["MESUp"]->GoodLeptons.push_back(thisMuonUp);
                        mainVars["MESUp"]->MetXCorr += thisMuon.Px() - thisMuonUp.Px();
                        mainVars["MESUp"]->MetYCorr += thisMuon.Py() - thisMuonUp.Py();
                    }
                    if (muonPtDown > MUON_VETO_CUT) {
                        mainVars["MESDown"]->nVetoMuons++;
                        mainVars["MESDown"]->GoodLeptons.push_back(thisMuonDown);
                        mainVars["MESDown"]->MetXCorr += thisMuon.Px() - thisMuonDown.Px();
                        mainVars["MESDown"]->MetYCorr += thisMuon.Py() - thisMuonDown.Py();
                    }
                }
                //Tight selection
                if (isTightMuon(i)) {
                    if (muonPtUp >= MUON_TIGHT_CUT) {
                        mainVars["MESUp"]->nTightMuons++;
                        if (muonPtUp > mainVars["MESUp"]->leadingTightMuPt){
                            mainVars["MESUp"]->leadingTightMu = thisMuonUp;
                            mainVars["MESUp"]->leadingTightMuPt = muonPtUp;
                        }
                    }
                    if (muonPtDown >= MUON_TIGHT_CUT) {
                        mainVars["MESDown"]->nTightMuons++;
                        if (muonPtDown > mainVars["MESDown"]->leadingTightMuPt){
                            mainVars["MESDown"]->leadingTightMu = thisMuonDown;
                            mainVars["MESDown"]->leadingTightMuPt = muonPtDown;
                        }
                    }
                }
            }
            //baseline pt cut
            if (muonPt[i] < MUON_VETO_CUT) continue;

            //Trigger scale factor
            //Note: The current scheme for applying the trigger scale factor is only correct for the 1-lepton box
            //For the 2-lepton boxes, the way we apply these scale factors now is wrong.
            if(!isData && muonPt[i] >= MUON_TIGHT_CUT){
                helper->updateSingleMuTriggerScaleFactors( muonPt[i], muonEta[i], isTightMuon(i), 
                        passedSingleLeptonTrigger, muonTrigCorrFactor, sf_muonTrigUp, sf_muonTrigDown );
            }
            //Veto selection
            if (isVetoMuon(i)){
                for (auto &vars : mainVars) {
                    if (vars.first != "MESUp" && vars.first != "MESDown") {
                        vars.second->nVetoMuons++;
                        vars.second->GoodLeptons.push_back(thisMuon);
                    }
                }
            }
            //Tight selection
            if (isTightMuon(i) && muonPt[i] >= MUON_TIGHT_CUT) {
                for (auto &vars : mainVars) {
                    if (vars.first != "MESUp" && vars.first != "MESDown") {
                        vars.second->nTightMuons++;
                        if (muonPt[i] > vars.second->leadingTightMuPt){
                            vars.second->leadingTightMu = thisMuon;
                            vars.second->leadingTightMuPt = muonPt[i];
                        }
                    }
                }
            }
        }




        /////////////////////////////////
        //Electron selection
        /////////////////////////////////
        float eleEffCorrFactor = 1.0;
        float vetoEleEffCorrFactor = 1.0;
        float eleTrigCorrFactor = 1.0;
        //Cut parameters
        const float ELE_VETO_CUT = 5;
        const float ELE_TIGHT_CUT = 30;
        //Loop electrons
        for (int i = 0; i < nElectrons; i++){
            //Remove overlaps
            bool overlap = false;
            for (auto& lep : mainVars[""]->GoodLeptons){
                if (RazorAnalyzer::deltaR(eleEta[i],elePhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
            }
            if (overlap) continue;
            //TLorentzVector for this electron
            TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
            //acceptance cut
            if (fabs(eleEta[i]) > 2.5) continue;
            //lepton scale uncertainty
            if (!isData) {
                float elePtUp = elePt[i];
                float elePtDown = elePt[i];
                float eleEUp = eleE[i];
                float eleEDown = eleE[i];
                //get up/down ele momenta
                float sfUp = 1.0;
                float sfDown = 1.0;
                if (fabs(eleEta[i]) < 1.5) {
                    sfUp = 1.006;
                    sfDown = 0.994;
                }
                else {
                    sfUp = 1.015;
                    sfDown = 0.985;
                }
                elePtUp *= sfUp;
                elePtDown *= sfDown;
                eleEUp *= sfUp;
                eleEDown *= sfDown;
                TLorentzVector thisElectronUp = makeTLorentzVector(elePtUp, eleEta[i], elePhi[i], eleEUp); 
                TLorentzVector thisElectronDown = makeTLorentzVector(elePtDown, eleEta[i], elePhi[i], eleEDown); 
                //Veto selection
                if (isVetoElectron(i)) {
                    if (elePtUp > ELE_VETO_CUT) {
                        mainVars["EESUp"]->nVetoElectrons++;
                        mainVars["EESUp"]->GoodLeptons.push_back(thisElectronUp); 
                        mainVars["EESUp"]->MetXCorr += thisElectron.Px() - thisElectronUp.Px();
                        mainVars["EESUp"]->MetYCorr += thisElectron.Py() - thisElectronUp.Py();
                    }
                    if (elePtDown > ELE_VETO_CUT) {
                        mainVars["EESDown"]->nVetoElectrons++;
                        mainVars["EESDown"]->GoodLeptons.push_back(thisElectronDown);
                        mainVars["EESDown"]->MetXCorr += thisElectron.Px() - thisElectronDown.Px();
                        mainVars["EESDown"]->MetYCorr += thisElectron.Py() - thisElectronDown.Py();
                    }
                }
                //Tight selection
                if (isTightElectron(i)) {
                    if (elePtUp >= ELE_TIGHT_CUT) {
                        mainVars["EESUp"]->nTightElectrons++;
                        if (elePtUp > mainVars["EESUp"]->leadingTightElePt){
                            mainVars["EESUp"]->leadingTightEle = thisElectronUp;
                            mainVars["EESUp"]->leadingTightElePt = elePtUp;
                        }
                    }
                    if (elePtDown >= ELE_TIGHT_CUT) {
                        mainVars["EESDown"]->nTightElectrons++;
                        if (elePtDown > mainVars["EESDown"]->leadingTightElePt){
                            mainVars["EESDown"]->leadingTightEle = thisElectronDown;
                            mainVars["EESDown"]->leadingTightElePt = elePtDown;
                        }
                    }
                }
            }
            //baseline pt cut
            if (elePt[i] < ELE_VETO_CUT) continue;

            //Trigger scale factor
            //Note: The current scheme for applying the trigger scale factor is only correct for the 1-lepton box
            //For the 2-lepton boxes, the way we apply these scale factors now is wrong.
            if(!isData && elePt[i] > ELE_TIGHT_CUT){
                helper->updateSingleEleTriggerScaleFactors( elePt[i], eleEta[i], isTightElectron(i), 
                        passedSingleLeptonTrigger, eleTrigCorrFactor, sf_eleTrigUp, sf_eleTrigDown );
            }
            //Veto selection
            if (isVetoElectron(i)){
                for (auto &vars : mainVars) {
                    if (vars.first != "EESUp" && vars.first != "EESDown") {
                        vars.second->nVetoElectrons++;
                        vars.second->GoodLeptons.push_back(thisElectron);            
                    }
                }
            }
            if (isTightElectron(i) && elePt[i] > ELE_TIGHT_CUT){
                for (auto &vars : mainVars) {
                    if (vars.first != "EESUp" && vars.first != "EESDown") {
                        vars.second->nTightElectrons++;
                        if (elePt[i] > vars.second->leadingTightElePt){
                            vars.second->leadingTightEle = thisElectron;
                            vars.second->leadingTightElePt = elePt[i];
                        }
                    }
                }
            }
        }


        /////////////////////////////////////////////////////////////////////////////////////
        //Compute Electron and Muon Selection Efficiency Correction Factors
        //Use Gen-Level Leptons as denominator
        /////////////////////////////////////////////////////////////////////////////////////
        if (!isData) {
            for(int j = 0; j < nGenParticle; j++){

                //look for electrons or muons
                if ( (abs(gParticleId[j]) == 11 ||abs(gParticleId[j]) == 13) 
                        && gParticleStatus[j] == 1 	      
                        && ( (abs(gParticleId[j]) == 11 && abs(gParticleEta[j]) < 2.5) ||  
                            (abs(gParticleId[j]) == 13 && abs(gParticleEta[j]) < 2.4))
                        && gParticlePt[j] > 5
                        &&
                        ( abs(gParticleMotherId[j]) == 24 
                          || abs(gParticleMotherId[j]) == 23 
                          || ( (abs(gParticleMotherId[j]) == 15 || abs(gParticleMotherId[j]) == 13 || abs(gParticleMotherId[j]) == 11 )
                              && gParticleMotherIndex[j] >= 0 
                              && (abs(gParticleMotherId[gParticleMotherIndex[j]]) == 24 || 
                                  abs(gParticleMotherId[gParticleMotherIndex[j]]) == 23)
                             )
                        )
                   )  {	      

                    //match to tight or veto ele
                    if (abs(gParticleId[j]) == 11) {
                        bool isSelectedTight = false;
                        bool isSelectedVeto = false;
                        for (int i = 0; i < nElectrons; i++){
                            if (fabs(eleEta[i]) > 2.5) continue;
                            if (elePt[i] < ELE_VETO_CUT) continue;
                            if (deltaR( eleEta[i], elePhi[i], gParticleEta[j], gParticlePhi[j]) > 0.1) continue;
                            if (elePt[i] > ELE_TIGHT_CUT && isTightElectron(i)) isSelectedTight = true;
                            if (isVetoElectron(i)) isSelectedVeto = true;
                        }

                        //if event passes single lepton trigger, apply correction factor for tight
                        if (passedSingleLeptonTrigger && gParticlePt[j] > ELE_TIGHT_CUT) {
                            helper->updateTightElectronScaleFactors(gParticlePt[j], gParticleEta[j], isSelectedTight,
                                    eleEffCorrFactor, sf_eleEffUp, sf_eleEffDown, 
                                    sf_eleEffFastsimSFUp, sf_eleEffFastsimSFDown);
                        }

                        //if event passes hadronic Trigger, apply correction factor for veto
                        if (passedHadronicTrigger) {
                            //for pT below 10 GeV, use the correction for 10 GeV
                            helper->updateVetoElectronScaleFactors( fmax( gParticlePt[j], 10.01) , gParticleEta[j], isSelectedVeto,
                                    vetoEleEffCorrFactor, sf_vetoEleEffUp, sf_vetoEleEffDown, 
                                    sf_vetoEleEffFastsimSFUp, sf_vetoEleEffFastsimSFDown);
                        }

                    } //end if electrons

                    //match to tight or veto muon
                    if (abs(gParticleId[j]) == 13) {
                        bool isSelectedTight = false;
                        bool isSelectedVeto = false;
                        for (int i = 0; i < nMuons; i++){
                            if (fabs(muonEta[i]) > 2.4) continue;
                            if (muonPt[i] < MUON_VETO_CUT) continue;
                            if (deltaR( muonEta[i], muonPhi[i], gParticleEta[j], gParticlePhi[j]) > 0.1) continue;
                            if (muonPt[i] > MUON_TIGHT_CUT && isTightMuon(i)) isSelectedTight = true;
                            if (isVetoMuon(i)) isSelectedVeto = true;
                        }		

                        //if event passes single lepton trigger, apply correction factor for tight
                        if (passedSingleLeptonTrigger && gParticlePt[j] > MUON_TIGHT_CUT) {
                            helper->updateTightMuonScaleFactors( gParticlePt[j], gParticleEta[j], isSelectedTight,
                                    muonEffCorrFactor, sf_muonEffUp, sf_muonEffDown, 
                                    sf_muonEffFastsimSFUp, sf_muonEffFastsimSFDown);
                        }

                        //if event passes hadronic Trigger, apply correction factor for veto
                        if (passedHadronicTrigger) {
                            //for pT below 10 GeV, use the correction for 10 GeV
                            helper->updateVetoMuonScaleFactors( fmax( gParticlePt[j], 10.01) , gParticleEta[j], isSelectedVeto,
                                    vetoMuonEffCorrFactor, sf_vetoMuonEffUp, sf_vetoMuonEffDown, 
                                    sf_vetoMuonEffFastsimSFUp, sf_vetoMuonEffFastsimSFDown );	
                        }

                    } //end if muons

                }//match gen leptons
            }//loop over gen particles
        }//end if data



        /////////////////////////////////
        //Tau selection
        /////////////////////////////////

        const float TAU_LOOSE_CUT = 20;
        float tauEffCorrFactor = 1.0;
        //Loop taus
        for (int i = 0; i < nTaus; i++){	 
            //Baseline cuts
            if (tauPt[i] < TAU_LOOSE_CUT) continue;
            if (fabs(tauEta[i]) > 2.4) continue;

            //remove overlaps
            bool overlap = false;
            for (auto& lep : mainVars[""]->GoodLeptons){
                if (RazorAnalyzer::deltaR(tauEta[i],tauPhi[i],lep.Eta(),lep.Phi()) < 0.4) overlap = true;
            }
            if (overlap) continue;

            //Loose selection
            if (isLooseTau(i)){
                nLooseTaus++;
                TLorentzVector thisTau = makeTLorentzVectorPtEtaPhiM(tauPt[i], tauEta[i], tauPhi[i], 1.777);
                for (auto &vars : mainVars) {
                    vars.second->GoodLeptons.push_back(thisTau);
                }
            }
        }

        /////////////////////////////////
        //Jet selection
        /////////////////////////////////

        //Hadronic trigger efficiency scale factor
        float hadronicTrigCorrFactor = 1.0; //flat trigger scale factor
        if (isFastsimSMS) {
            hadronicTrigCorrFactor *= 0.975;
        }
        //Jet cuts
        const int JET_CUT = 40;
        const int BJET_CUT = 40;
        //Loop jets
        for (int i = 0; i < nJets; i++){
            //Apply Jet ID only on fullsim. fastsim jet ID is broken. 
            if (!isFastsimSMS) {
                if (!jetPassIDTight[i]) continue;
            }
            //Get jet energy correction
            double tmpRho = fixedGridRhoFastjetAll;
            double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
                    tmpRho, jetJetArea[i], 
                    runNum,
                    JetCorrectorIOV,JetCorrector);   
            //Get jet energy resolution correction, with up/down variants
            double jetEnergySmearFactor = 1.0;
            double jetEnergySmearFactorUp = 1.0;
            double jetEnergySmearFactorDown = 1.0;
            //Get L1-only jet energy correction
            double JECLevel1 = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
                    tmpRho, jetJetArea[i], runNum, 
                    JetCorrectorIOV,JetCorrector, 0);   
            //TLorentzVector for this jet
            double jetCorrPt = jetPt[i]*JEC*jetEnergySmearFactor;
            double jetCorrE = jetE[i]*JEC*jetEnergySmearFactor;
            TLorentzVector thisJet = makeTLorentzVector(jetCorrPt, jetEta[i], jetPhi[i], jetCorrE);
            TLorentzVector L1CorrJet = makeTLorentzVector(jetPt[i]*JECLevel1, jetEta[i], jetPhi[i], 
                    jetE[i]*JECLevel1);
            //for MET and jet variables
            bool matchesLepton = false;
            for (auto &vars : mainVars) {
                // loop over leptons to find the closest
                double deltaR = -1;
                for (auto &lep : vars.second->GoodLeptons) {
                    double thisDR = RazorAnalyzer::deltaR(jetEta[i], jetPhi[i], lep.Eta(), lep.Phi());  
                    if (deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
                }
                if (deltaR < 0 || deltaR > 0.4) { // jet does not match a lepton
                    if (jetCorrPt > 15 && 
                            jetChargedEMEnergyFraction[i] + jetNeutralEMEnergyFraction[i] <= 0.9) {
                        // correct the MET
                        vars.second->MetXCorr += -1 * (thisJet.Px() - L1CorrJet.Px());
                        vars.second->MetYCorr += -1 * (thisJet.Py() - L1CorrJet.Py());
                    }
                    if (vars.first != "JESUp" && vars.first != "JESDown" && 
                            vars.first != "JERUp" && vars.first != "JERDown") { //these ones are handled below
                        if (jetCorrPt > BJET_CUT && fabs(jetEta[i]) < 3.0 && isCSVM(i)){
                            // count it as a b-jet
                            vars.second->nBTaggedJets++;
                        }
                        if (jetCorrPt > JET_CUT && fabs(jetEta[i]) < 3.0) {
                            // add to good jets list
                            vars.second->GoodJets.push_back(thisJet);
                            vars.second->nSelectedJets++;
                            if (jetCorrPt > 80) vars.second->nJets80++;
                        }
                    }
                }
                else if (vars.first == "") { //jet matches a lepton
                    matchesLepton = true;
                }
            }
            //Remove overlaps
            if (matchesLepton) continue;

            //Count Number of Gen-Level Matched BJets
            if (abs(jetPartonFlavor[i]) == 5 && jetCorrPt > 40 && fabs(jetEta[i]) < 2.4) NGenBJets++;

            //Apply b-tagging correction factor 
            if (!isData && abs(jetEta[i]) < 2.4 && jetCorrPt > BJET_CUT) { 
                helper->updateBTagScaleFactors( jetCorrPt, jetEta[i], jetPartonFlavor[i], isCSVM(i),
                        btagCorrFactor, sf_btagUp, sf_btagDown, sf_btagFastsimSFUp, sf_btagFastsimSFDown,
                        sf_bmistagUp, sf_bmistagDown );
            }
            //Cut on jet eta
            if (fabs(jetEta[i]) > 3.0) continue;
            //Get uncertainty on JEC and JER
            if(!isData){
                double unc = helper->getJecUnc( jetCorrPt, jetEta[i], 999 ); //use run=999 as default
                double jetPtJESUp = jetCorrPt*(1+unc);
                double jetPtJESDown = jetCorrPt/(1+unc);
                double jetPtJERUp = jetPt[i]*JEC*jetEnergySmearFactorUp;
                double jetPtJERDown = jetPt[i]*JEC*jetEnergySmearFactorDown;
                double jetEJESUp = jetCorrE*(1+unc);
                double jetEJESDown = jetCorrE/(1+unc);
                double jetEJERUp = jetE[i]*JEC*jetEnergySmearFactorUp;
                double jetEJERDown = jetE[i]*JEC*jetEnergySmearFactorDown;
                TLorentzVector thisJetJESUp = makeTLorentzVector(jetPtJESUp, jetEta[i], jetPhi[i], jetEJESUp);
                TLorentzVector thisJetJESDown = makeTLorentzVector(jetPtJESDown, jetEta[i], jetPhi[i], jetEJESDown);
                TLorentzVector thisJetJERUp = makeTLorentzVector(jetPtJERUp, jetEta[i], jetPhi[i], jetEJERUp);
                TLorentzVector thisJetJERDown = makeTLorentzVector(jetPtJERDown, jetEta[i], jetPhi[i], jetEJERDown);
                //Propagate uncertainties to the MET
                if (jetPtJESUp > 20) {
                    mainVars["JESUp"]->MetXCorr += -1 * (thisJetJESUp.Px() - thisJet.Px());
                    mainVars["JESUp"]->MetYCorr += -1 * (thisJetJESUp.Py() - thisJet.Py());
                }
                if (jetPtJESDown > 20) {
                    mainVars["JESDown"]->MetXCorr += -1 * (thisJetJESDown.Px() - thisJet.Px());
                    mainVars["JESDown"]->MetYCorr += -1 * (thisJetJESDown.Py() - thisJet.Py());
                }
                if (jetPtJERUp > 20) {
                    mainVars["JERUp"]->MetXCorr += -1 * (thisJetJERUp.Px() - thisJet.Px());
                    mainVars["JERUp"]->MetYCorr += -1 * (thisJetJERUp.Py() - thisJet.Py());
                }
                if (jetPtJERDown > 20) {
                    mainVars["JERDown"]->MetXCorr += -1 * (thisJetJERDown.Px() - thisJet.Px());
                    mainVars["JERDown"]->MetYCorr += -1 * (thisJetJERDown.Py() - thisJet.Py());
                }
                //Record jets that pass the cut
                if(jetPtJESUp > BJET_CUT && isCSVM(i)) mainVars["JESUp"]->nBTaggedJets++;
                if(jetPtJESDown > BJET_CUT && isCSVM(i)) mainVars["JESDown"]->nBTaggedJets++;
                if(jetPtJERUp > BJET_CUT && isCSVM(i)) mainVars["JERUp"]->nBTaggedJets++;
                if(jetPtJERDown > BJET_CUT && isCSVM(i)) mainVars["JERDown"]->nBTaggedJets++;
                if(jetPtJESUp > JET_CUT){
                    mainVars["JESUp"]->GoodJets.push_back(thisJetJESUp);
                    mainVars["JESUp"]->nSelectedJets++;
                    if (thisJetJESUp.Pt() > 80) mainVars["JESUp"]->nJets80++;
                }
                if(jetPtJESDown > JET_CUT){
                    mainVars["JESDown"]->GoodJets.push_back(thisJetJESDown);
                    mainVars["JESDown"]->nSelectedJets++;
                    if (thisJetJESDown.Pt() > 80) mainVars["JESDown"]->nJets80++;
                }
                if(jetPtJERUp > JET_CUT){
                    mainVars["JERUp"]->GoodJets.push_back(thisJetJERUp);
                    mainVars["JERUp"]->nSelectedJets++;
                    if (thisJetJERUp.Pt() > 80) mainVars["JERUp"]->nJets80++;
                }
                if(jetPtJERDown > JET_CUT){
                    mainVars["JERDown"]->GoodJets.push_back(thisJetJERDown);
                    mainVars["JERDown"]->nSelectedJets++;
                    if (thisJetJERDown.Pt() > 80) mainVars["JERDown"]->nJets80++;
                }
            }
        }

        //Get leading and subleading jet pt
        for (auto &vars : mainVars) {
            TLorentzVector leadingJet;
            TLorentzVector subleadingJet;
            for (auto &jet : vars.second->GoodJets) {
                if (jet.Pt() > vars.second->leadingJetPt){
                    vars.second->subleadingJetPt = vars.second->leadingJetPt;
                    subleadingJet = leadingJet;
                    vars.second->leadingJetPt = jet.Pt();
                    leadingJet = jet;

                }
                else if (jet.Pt() > vars.second->subleadingJetPt){
                    vars.second->subleadingJetPt = jet.Pt();
                    subleadingJet = jet;
                }
            }
            if (vars.first == "") mjj_leadingJets = (leadingJet + subleadingJet).M();
        }

        /////////////////////////////////
        //Compute razor variables and mT
        /////////////////////////////////

        //Combine jet and lepton collections
        for (auto &vars : mainVars) {
            for (auto &lep : vars.second->GoodLeptons) {
                vars.second->GoodJets.push_back(lep);
            }
        }

        //Get HT
        HT = 0;
        for (auto& jet : mainVars[""]->GoodJets) HT += jet.Pt();

        for (auto &vars : mainVars) {
            // Make MET
            double PFMetX = metPt*cos(metPhi) + vars.second->MetXCorr;
            double PFMetY = metPt*sin(metPhi) + vars.second->MetYCorr;
            TLorentzVector MyMET;
            MyMET.SetPxPyPzE( PFMetX, PFMetY, 0, sqrt( PFMetX*PFMetX + PFMetY*PFMetY ) );
            // Compute MR, Rsq, dPhiRazor
            vector<TLorentzVector> hemispheres = getHemispheres(vars.second->GoodJets);
            vars.second->MR = computeMR(hemispheres[0], hemispheres[1]); 
            vars.second->Rsq = computeRsq(hemispheres[0], hemispheres[1], MyMET);
            vars.second->dPhiRazor = deltaPhi(hemispheres[0].Phi(),hemispheres[1].Phi());
            vars.second->metOverCaloMet = MyMET.Pt()/metCaloPt;
            // Compute transverse mass 
            if (vars.second->nTightMuons + vars.second->nTightElectrons > 0){
                TLorentzVector leadingLepton;
                if (vars.second->leadingTightMuPt > vars.second->leadingTightElePt) {
                    leadingLepton = vars.second->leadingTightMu;
                }
                else {
                    leadingLepton = vars.second->leadingTightEle;
                }
                float deltaPhiLepMet = leadingLepton.DeltaPhi(MyMET);
                vars.second->mT = sqrt(2*leadingLepton.Pt()*MyMET.Pt()*(1.0 - cos(deltaPhiLepMet))); 
            }
            // Transverse mass with leading lepton, regardless of quality
            if (vars.second->GoodLeptons.size() > 0){
                //get the highest pt lepton
                float maxLepPt = -1;
                int maxLepIndex = -1;
                for (uint i = 0; i < vars.second->GoodLeptons.size(); i++){
                    if (vars.second->GoodLeptons[i].Pt() > maxLepPt){
                        maxLepPt = vars.second->GoodLeptons[i].Pt();
                        maxLepIndex = i;
                    }
                }
                if (maxLepIndex >= 0){
                    float deltaPhiLepMet = vars.second->GoodLeptons[maxLepIndex].DeltaPhi(MyMET);
                    vars.second->mTLoose = sqrt(2*vars.second->GoodLeptons[maxLepIndex].Pt()*MyMET.Pt()*(1.0 
                                - cos(deltaPhiLepMet)));
                }
            }
            // Additional quantities
            if (vars.first == "") {
                met = MyMET.Pt();
                mjj_hemispheres = (hemispheres[0] + hemispheres[1]).M();
            }
        }

        //////////////////////////////////////////////////////
        //Check for any photons with Ecal Gain Switch
        //////////////////////////////////////////////////////
        Flag_hasEcalGainSwitch = false;
        for (int i = 0; i < nPhotons; i++){
            if (phoPt[i] > 50
                    && 
                    (pho_seedRecHitSwitchToGain6[i] || 
                     pho_seedRecHitSwitchToGain1[i] || 
                     pho_anyRecHitSwitchToGain6[i] || 
                     pho_anyRecHitSwitchToGain1[i]  
                    )
               ) {
                Flag_hasEcalGainSwitch = true;
            }
        }

        /////////////////////////////////
        //Categorize into boxes
        /////////////////////////////////

        for (auto &vars : mainVars) {
            if (passedDileptonTrigger && vars.second->nTightElectrons > 0 && vars.second->nTightMuons > 0){
                vars.second->box = MuEle;
            }
            else if (passedDileptonTrigger && vars.second->nTightMuons > 1){
                vars.second->box = MuMu;
            }
            else if (passedDileptonTrigger && vars.second->nTightElectrons > 1){
                vars.second->box = EleEle;
            }
            else if (passedSingleLeptonTrigger && vars.second->nTightMuons > 0){
                if (vars.second->nSelectedJets > 5) vars.second->box = MuSixJet;
                else if (vars.second->nSelectedJets > 3) vars.second->box = MuFourJet;
                else vars.second->box = MuJet;
            }
            else if (passedSingleLeptonTrigger && vars.second->nTightElectrons > 0){
                if (vars.second->nSelectedJets > 5) vars.second->box = EleSixJet;
                else if (vars.second->nSelectedJets > 3) vars.second->box = EleFourJet;
                else vars.second->box = EleJet;
            }
            else if (passedHadronicTrigger && nLooseTaus + vars.second->nVetoElectrons + 
                    vars.second->nVetoMuons > 0 && vars.second->nJets80 >= 2){
                if (vars.second->nSelectedJets > 5) vars.second->box = LooseLeptonSixJet;
                else if (vars.second->nSelectedJets > 3) vars.second->box = LooseLeptonFourJet;
                else vars.second->box = LooseLeptonDiJet;
            }
            else if (passedHadronicTrigger && vars.second->nJets80 >= 2 && nLooseTaus
                    + vars.second->nVetoElectrons + vars.second->nVetoMuons 
                    + vars.second->nTightElectrons + vars.second->nTightMuons == 0){
                if (vars.second->nSelectedJets > 5) vars.second->box = SixJet;
                else if (vars.second->nSelectedJets > 3) vars.second->box = FourJet;
                else vars.second->box = DiJet;
            }
        }

        //***********************************
        //Filter out the Pathological Events
        //***********************************
        if (isFastsimSMS) {
            bool isPathologicalFastsimEvent = false;		  

            for (int i = 0; i < nJets; i++){
                double JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], 
                        fixedGridRhoFastjetAll, jetJetArea[i], runNum, JetCorrectorIOV,JetCorrector);   	  
                double jetCorrPt = jetPt[i]*JEC;
                if (jetCorrPt < 20) continue;
                if (fabs(jetEta[i]) > 2.5) continue;

                //Match to Gen Jet
                bool isMatch = false;
                for(int j = 0; j < nGenJets; j++){
                    double tmpDR = deltaR( genJetEta[j],genJetPhi[j], jetEta[i],jetPhi[i]);
                    if ( tmpDR < 0.4
                       ) {	
                        isMatch = true;
                    }
                }

                // these are the pathological fastsim jets
                if (!isMatch && jetChargedHadronEnergyFraction[i] < 0.1 ) {
                    isPathologicalFastsimEvent = true;
                }
            }
            //reject event if it's pathological
            if (isPathologicalFastsimEvent) continue;	  
        }

        /////////////////////////////////
        //Scale and PDF variations
        /////////////////////////////////

        if ((*scaleWeights).size() >= 9) {
            sf_facScaleUp = (*scaleWeights)[1]/genWeight;
            sf_facScaleDown = (*scaleWeights)[2]/genWeight;
            sf_renScaleUp = (*scaleWeights)[3]/genWeight;
            sf_renScaleDown = (*scaleWeights)[6]/genWeight;
            sf_facRenScaleUp = (*scaleWeights)[4]/genWeight;
            sf_facRenScaleDown = (*scaleWeights)[8]/genWeight;
        }

        SumScaleWeights->Fill(0.0, sf_facScaleUp);
        SumScaleWeights->Fill(1.0, sf_facScaleDown);
        SumScaleWeights->Fill(2.0, sf_renScaleUp);
        SumScaleWeights->Fill(3.0, sf_renScaleDown);
        SumScaleWeights->Fill(4.0, sf_facRenScaleUp);
        SumScaleWeights->Fill(5.0, sf_facRenScaleDown);

        for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) {
            SumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
        }

        /////////////////////////////////
        //Apply scale factors
        /////////////////////////////////

        //Nominal event weight
        if(!isData){
            weight *= pileupWeight; 

            if (passedSingleLeptonTrigger && 
                    (mainVars[""]->box == MuEle || mainVars[""]->box == MuMu || mainVars[""]->box == EleEle || 
                     mainVars[""]->box == MuSixJet || mainVars[""]->box == MuFourJet || mainVars[""]->box == MuJet 
                     || mainVars[""]->box == EleSixJet || mainVars[""]->box == EleFourJet || mainVars[""]->box == EleJet
                    )
               ) {
                weight *= muonEffCorrFactor;
                weight *= muonTrigCorrFactor;
                weight *= eleEffCorrFactor;
                weight *= eleTrigCorrFactor;
            }
            else if(passedHadronicTrigger &&
                    (mainVars[""]->box == LooseLeptonSixJet || mainVars[""]->box == LooseLeptonFourJet 
                     || mainVars[""]->box == SixJet || mainVars[""]->box == FourJet
                     ||  mainVars[""]->box == LooseLeptonDiJet || mainVars[""]->box == DiJet
                    )
                   ) {
                weight *= vetoMuonEffCorrFactor;
                weight *= vetoEleEffCorrFactor;
                weight *= tauEffCorrFactor;
                weight *= hadronicTrigCorrFactor;
            }
            weight *= btagCorrFactor;   

            // if (weight < 0.1) {
            // cout << "weight: " << weight << " | "
            // 	   << pileupWeight << " " << NPU << " | "
            // 	   << passedSingleLeptonTrigger << " / " << passedHadronicTrigger << " | "
            // 	   << muonEffCorrFactor << " "
            // 	   << muonTrigCorrFactor << " "
            // 	   << eleEffCorrFactor << " " 
            // 	   << eleTrigCorrFactor << " "
            // 	   << vetoMuonEffCorrFactor << " "
            // 	   << vetoEleEffCorrFactor << " "
            // 	   << tauEffCorrFactor << " "
            // 	   << hadronicTrigCorrFactor << " "
            // 	   << btagCorrFactor << " "
            // 	   << "\n";
            // }

        }

        //Fill normalization histogram
        NEvents->SetBinContent( 1, NEvents->GetBinContent(1) + genWeight);
        SumTopPtWeights->SetBinContent( 1, SumTopPtWeights->GetBinContent(1) + topPtWeight );
        SumWeights->Fill(1.0, weight);

        /////////////////////////////////
        //SMS information
        /////////////////////////////////

        bool parsedLHE = false;
        if(isFastsimSMS && lheComments){
            //parse lhe comment string to get gluino and LSP masses
            stringstream parser(*lheComments);
            string item;
            getline(parser, item, '_'); //prefix
            if(getline(parser, item, '_')){ //gluino mass 
                mGluino = atoi(item.c_str());
                if(mGluino == 0) { //fix for the case where the model name contains an underscore
                    if(getline(parser, item, '_')){
                        mGluino = atoi(item.c_str());
                    }
                }
                if(getline(parser, item, '_')){ //LSP mass 
                    mLSP = atoi(item.c_str());
                    pair<int,int> smsPair = make_pair(mGluino, mLSP);
                    parsedLHE = true;
                    if (smsFiles.count(smsPair) == 0){ //create file and tree
                        //format file name
                        string thisFileName = outFileName;
                        thisFileName.erase(thisFileName.end()-5, thisFileName.end());
                        thisFileName += "_" + to_string(mGluino) + "_" + to_string(mLSP) + ".root";

                        smsFiles[smsPair] = new TFile(thisFileName.c_str(), "recreate");
                        smsTrees[smsPair] = razorTree->CloneTree(0);
                        smsNEvents[smsPair] = new TH1F(Form("NEvents%d%d", mGluino, mLSP), "NEvents", 1,0.5,1.5);
                        smsSumWeights[smsPair] = new TH1F(Form("SumWeights%d%d", mGluino, mLSP), "SumWeights", 1,0.5,1.5);
                        smsSumScaleWeights[smsPair] = new TH1F(Form("SumScaleWeights%d%d", mGluino, mLSP), "SumScaleWeights", 6,-0.5,5.5);
                        smsSumPdfWeights[smsPair] = new TH1F(Form("SumPdfWeights%d%d", mGluino, mLSP), "SumPdfWeights", NUM_PDF_WEIGHTS,-0.5,NUM_PDF_WEIGHTS-0.5);
                        cout << "Created new output file " << thisFileName << endl;
                    }
                    //Fill NEvents hist 
                    smsNEvents[smsPair]->Fill(1.0, genWeight);
                    smsSumWeights[smsPair]->Fill(1.0, weight);

                    smsSumScaleWeights[smsPair]->Fill(0.0, sf_facScaleUp);
                    smsSumScaleWeights[smsPair]->Fill(1.0, sf_facScaleDown);
                    smsSumScaleWeights[smsPair]->Fill(2.0, sf_renScaleUp);
                    smsSumScaleWeights[smsPair]->Fill(3.0, sf_renScaleDown);
                    smsSumScaleWeights[smsPair]->Fill(4.0, sf_facRenScaleUp);
                    smsSumScaleWeights[smsPair]->Fill(5.0, sf_facRenScaleDown);

                    for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) {
                        smsSumPdfWeights[smsPair]->Fill(double(iwgt),(*pdfWeights)[iwgt]);
                    }
                }
            }
        }

        /////////////////////////////////
        //Baseline cuts
        /////////////////////////////////

        //Razor
        bool passCuts = false;
        for (auto &vars : mainVars) {
            if (vars.second->MR > 150 && vars.second->Rsq > 0.15 && vars.second->box != NONE) passCuts = true;
        }
        if (!passCuts) continue;

        //Trigger
        if(!passedDileptonTrigger && !passedSingleLeptonTrigger && !passedHadronicTrigger) {
            continue;
        }

        /////////////////////////////////
        //Noise filters
        /////////////////////////////////

        if(!isFastsimSMS){
            if(!Flag_HBHENoiseFilter) continue;
            if(!Flag_HBHEIsoNoiseFilter) continue;
            if(!Flag_goodVertices) continue;
            if(!Flag_eeBadScFilter) continue;
        }

        //Fill tree
        if(!isFastsimSMS){
            razorTree->Fill();
        }
        else if(parsedLHE){
            pair<int,int> smsPair = make_pair(mGluino, mLSP);
            smsTrees[smsPair]->Fill();
        }

    }//end of event loop

    cout << "NBytes: " << nbytes << "\n";

    if(!isFastsimSMS){
        cout << "Writing output tree..." << endl;
        outFile->cd();
        razorTree->Write();
        NEvents->Write();
        SumTopPtWeights->Write();
        SumWeights->Write();
        SumScaleWeights->Write();
        SumPdfWeights->Write();
    }
    else{
        for(auto &filePtr : smsFiles){
            cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
            filePtr.second->cd();
            smsTrees[filePtr.first]->Write();
            smsNEvents[filePtr.first]->Write("NEvents");
            smsSumWeights[filePtr.first]->Write("SumWeights");
            smsSumScaleWeights[filePtr.first]->Write("SumScaleWeights");
            smsSumPdfWeights[filePtr.first]->Write("SumPdfWeights");
        }
    }

    outFile->Close();

    delete helper;

}
