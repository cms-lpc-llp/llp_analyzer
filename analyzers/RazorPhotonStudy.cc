//Produces trees for estimating lepton selection efficiency (using gen-level information) and for studying the use of DY+Jets, W+Jets, and G+Jets control regions to model the distribution of the Z->invisible background in the razor analysis.

#include "RazorPhotonStudy.h"
#include "JetCorrectorParameters.h"

//C++ includes

//ROOT includes
#include "TH1F.h"

using namespace std;

struct greater_than_pt{
    inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){
        return p1.Pt() > p2.Pt();
    }
};

void RazorPhotonStudy::Analyze(bool isData, int option, string outputfilename, string label)
{
    //****************************************************//
    //            Initialization of the tree              //
    //****************************************************//

    cout << "Initializing..." << endl;
    bool isRunOne = true;
    bool filterEvents = !(option == 1);

    //random number generator for jet smearing
    TRandom3 *random = new TRandom3();
    random->SetSeed(33333);

    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "RazorPhotonStudy.root";
    TFile outFile(outfilename.c_str(), "RECREATE");

    //one tree to hold all events
    TTree *razorTree = new TTree("RazorInclusive", "Info on selected razor inclusive events");

    //histogram containing total number of processed events (for normalization)
    TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);

    //tree variables
    int nVtx, nPU_mean;
    UInt_t run, lumi, event;
    bool hlt_dimuon, hlt_singlemu, hlt_photon, hlt_razor;
    float hlt_photon_weight;
    int nSelectedJets, nBTaggedJets;
    int nVetoMuons, nLooseMuons, nTightMuons;
    int nVetoElectrons, nLooseElectrons, nTightElectrons;
    int nLooseTaus, nMediumTaus, nTightTaus;
    int nGenMuons, nGenElectrons, nGenTauMuons, nGenTauElectrons, nGenTaus, nGenPhotons;
    float theMR, MR_noPho, MR_noZ, MR_noW, MR_noGenZ;
    float theRsq, Rsq_noPho, Rsq_noZ, Rsq_noW, Rsq_noGenZ; 
    float deltaPhi, deltaPhi_noPho, deltaPhi_noZ, deltaPhi_noW, deltaPhi_noGenZ; 
    float met,genmet,met_noPho, met_noZ, met_noW, genZmass, met_noGenZ;
    float leadingGenMuonPt, leadingGenPhotonPt;
    float leadingGenMuonEta, leadingGenPhotonEta;
    float leadingGenMuonPhi, leadingGenPhotonPhi;
    float leadingGenMuonE, leadingGenPhotonE;
    float subleadingGenMuonPt, subleadingGenPhotonPt;
    float subleadingGenMuonEta, subleadingGenPhotonEta;
    float subleadingGenMuonPhi, subleadingGenPhotonPhi;
    float subleadingGenMuonE, subleadingGenPhotonE;
    float leadingMuonPt, leadingTightMuonPt, leadingPhotonPt;
    float leadingMuonEta, leadingTightMuonEta, leadingPhotonEta;
    float leadingMuonPhi, leadingTightMuonPhi, leadingPhotonPhi;
    float leadingMuonE, leadingTightMuonE, leadingPhotonE;
    float subleadingMuonPt, subleadingPhotonPt;
    float subleadingMuonEta, subleadingPhotonEta;
    float subleadingMuonPhi, subleadingPhotonPhi;
    float subleadingMuonE, subleadingPhotonE;
    float metphi, genmetphi, metphi_noZ, metphi_noW, metphi_noPho, metphi_noGenZ;
    float HT, HT_noZ, HT_noW, HT_noPho, HT_noGenZ;
    int numJets, numJets_noZ, numJets_noW, numJets_noPho, numJets_noGenZ; 
    int numJets80, numJets80_noZ, numJets80_noW, numJets80_noPho, numJets80_noGenZ; 
    float genZpt, recoZpt, genZeta, recoZeta, genZphi, recoZphi, recoZmass, genWpt, recoWpt, genWeta, genWphi, recoWphi;
    bool leadGenMuonIsFound, leadGenPhotonIsFound;
    bool leadGenMuonIsFoundTight;
    float ptMatchingLeadGenMuon, ptMatchingLeadGenPhoton; //pt of the object matching the gen particle
    int nSelectedPhotons;    
    float mTLepMet;
    bool passedHLTPhoton50, passedHLTPhoton75, passedHLTPhoton90, passedHLTPhoton135, passedHLTPhoton150, passedHLTPhoton160;
    int nBJetsLoose20GeV, nBJetsMedium20GeV, nBJetsTight20GeV;
    bool bjet1PassLoose, bjet1PassMedium, bjet1PassTight, bjet2PassLoose, bjet2PassMedium, bjet2PassTight;
    float bjet1Pt, bjet2Pt;

    //set branches on big tree
    if(!isData){
        razorTree->Branch("nGenMuons", &nGenMuons, "nGenMuons/I");
        razorTree->Branch("nGenPhotons", &nGenPhotons, "nGenPhotons/I");
        razorTree->Branch("leadingGenMuonPt", &leadingGenMuonPt, "leadingGenMuonPt/F");
        razorTree->Branch("leadingGenMuonEta", &leadingGenMuonEta, "leadingGenMuonEta/F");
        razorTree->Branch("leadingGenMuonPhi", &leadingGenMuonPhi, "leadingGenMuonPhi/F");
        razorTree->Branch("leadingGenMuonE", &leadingGenMuonE, "leadingGenMuonE/F");
        razorTree->Branch("leadingGenPhotonPt", &leadingGenPhotonPt, "leadingGenPhotonPt/F");
        razorTree->Branch("leadingGenPhotonEta", &leadingGenPhotonEta, "leadingGenPhotonEta/F");
        razorTree->Branch("leadingGenPhotonPhi", &leadingGenPhotonPhi, "leadingGenPhotonPhi/F");
        razorTree->Branch("leadingGenPhotonE", &leadingGenPhotonE, "leadingGenPhotonE/F");
        razorTree->Branch("subleadingGenMuonPt", &subleadingGenMuonPt, "subleadingGenMuonPt/F");
        razorTree->Branch("subleadingGenMuonEta", &subleadingGenMuonEta, "subleadingGenMuonEta/F");
        razorTree->Branch("subleadingGenMuonPhi", &subleadingGenMuonPhi, "subleadingGenMuonPhi/F");
        razorTree->Branch("subleadingGenMuonE", &subleadingGenMuonE, "subleadingGenMuonE/F");
        razorTree->Branch("genZpt", &genZpt, "genZpt/F");
        razorTree->Branch("genZeta", &genZeta, "genZeta/F");
        razorTree->Branch("genZphi", &genZphi, "genZphi/F");
        razorTree->Branch("genZmass", &genZmass, "genZmass/F");
        razorTree->Branch("leadGenMuonIsFound", &leadGenMuonIsFound, "leadGenMuonIsFound/O");
        razorTree->Branch("leadGenPhotonIsFound", &leadGenPhotonIsFound, "leadGenPhotonIsFound/O");
        razorTree->Branch("leadGenMuonIsFoundTight", &leadGenMuonIsFoundTight, "leadGenMuonIsFoundTight/O");
        if(filterEvents){
            razorTree->Branch("nGenElectrons", &nGenElectrons, "nGenElectrons/I");
            razorTree->Branch("nGenTauMuons", &nGenTauMuons, "nGenTauMuons/I");
            razorTree->Branch("nGenTauElectrons", &nGenTauElectrons, "nGenTauElectrons/I");
            razorTree->Branch("nGenTaus", &nGenTaus, "nGenTaus/I");
            razorTree->Branch("subleadingGenPhotonPt", &subleadingGenPhotonPt, "subleadingGenPhotonPt/F");
            razorTree->Branch("subleadingGenPhotonEta", &subleadingGenPhotonEta, "subleadingGenPhotonEta/F");
            razorTree->Branch("subleadingGenPhotonPhi", &subleadingGenPhotonPhi, "subleadingGenPhotonPhi/F");
            razorTree->Branch("subleadingGenPhotonE", &subleadingGenPhotonE, "subleadingGenPhotonE/F");
            razorTree->Branch("MR_noGenZ", &MR_noGenZ, "MR_noGenZ/F");
            razorTree->Branch("Rsq_noGenZ", &Rsq_noGenZ, "Rsq_noGenZ/F");
            razorTree->Branch("met_noGenZ", &met_noGenZ, "met_noGenZ/F");
            razorTree->Branch("deltaPhi_noGenZ", &deltaPhi_noGenZ, "deltaPhi_noGenZ/F");
            razorTree->Branch("metphi_noGenZ", &metphi_noGenZ, "metphi_noGenZ/F");
            razorTree->Branch("genmet", &genmet, "genmet/F");
            razorTree->Branch("genmetphi", &genmetphi, "genmetphi/F");
            razorTree->Branch("HT_noGenZ", &HT_noGenZ, "HT_noGenZ/F");
            razorTree->Branch("numJets_noGenZ", &numJets_noGenZ, "numJets_noGenZ/I");
            razorTree->Branch("numJets80_noGenZ", &numJets80_noGenZ, "numJets80_noGenZ/I");
            razorTree->Branch("genWpt", &genWpt, "genWpt/F");
            razorTree->Branch("genWeta", &genWeta, "genWeta/F");
            razorTree->Branch("genWphi", &genWphi, "genWphi/F");
            razorTree->Branch("ptMatchingLeadGenMuon", &ptMatchingLeadGenMuon, "ptMatchingLeadGenMuon/F");
            razorTree->Branch("ptMatchingLeadGenPhoton", &ptMatchingLeadGenPhoton, "ptMatchingLeadGenPhoton/F");
        }
    }
    razorTree->Branch("run", &run, "run/i");
    razorTree->Branch("lumi", &lumi, "lumi/i");
    razorTree->Branch("event", &event, "event/i");
    razorTree->Branch("nVtx", &nVtx, "nVtx/I");
    razorTree->Branch("nPU_mean", &nPU_mean, "nPU_mean/I");
    razorTree->Branch("hlt_dimuon", &hlt_dimuon, "hlt_dimuon/O");
    razorTree->Branch("hlt_singlemu", &hlt_singlemu, "hlt_singlemu/O");
    razorTree->Branch("hlt_photon", &hlt_photon, "hlt_photon/O");
    razorTree->Branch("hlt_photon_weight", &hlt_photon_weight, "hlt_photon_weight/F");
    razorTree->Branch("passedHLTPhoton50", &passedHLTPhoton50, "passedHLTPhoton50/O");
    razorTree->Branch("passedHLTPhoton75", &passedHLTPhoton75, "passedHLTPhoton75/O");
    razorTree->Branch("passedHLTPhoton90", &passedHLTPhoton90, "passedHLTPhoton90/O");
    razorTree->Branch("passedHLTPhoton135", &passedHLTPhoton135, "passedHLTPhoton135/O");
    razorTree->Branch("passedHLTPhoton150", &passedHLTPhoton150, "passedHLTPhoton150/O");
    razorTree->Branch("passedHLTPhoton160", &passedHLTPhoton160, "passedHLTPhoton160/O");
    razorTree->Branch("hlt_razor", &hlt_razor, "hlt_razor/O");
    razorTree->Branch("nLooseMuons", &nLooseMuons, "nLooseMuons/I");
    razorTree->Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
    razorTree->Branch("nSelectedPhotons", &nSelectedPhotons, "nSelectedPhotons/I");
    razorTree->Branch("leadingMuonPt", &leadingMuonPt, "leadingMuonPt/F");
    razorTree->Branch("leadingMuonEta", &leadingMuonEta, "leadingMuonEta/F");
    razorTree->Branch("leadingMuonPhi", &leadingMuonPhi, "leadingMuonPhi/F");
    razorTree->Branch("leadingMuonE", &leadingMuonE, "leadingMuonE/F");
    razorTree->Branch("leadingTightMuonPt", &leadingTightMuonPt, "leadingTightMuonPt/F");
    razorTree->Branch("leadingTightMuonEta", &leadingTightMuonEta, "leadingTightMuonEta/F");
    razorTree->Branch("leadingTightMuonPhi", &leadingTightMuonPhi, "leadingTightMuonPhi/F");
    razorTree->Branch("leadingTightMuonE", &leadingTightMuonE, "leadingTightMuonE/F");
    razorTree->Branch("leadingPhotonPt", &leadingPhotonPt, "leadingPhotonPt/F");
    razorTree->Branch("leadingPhotonEta", &leadingPhotonEta, "leadingPhotonEta/F");
    razorTree->Branch("leadingPhotonPhi", &leadingPhotonPhi, "leadingPhotonPhi/F");
    razorTree->Branch("leadingPhotonE", &leadingPhotonE, "leadingPhotonE/F");
    razorTree->Branch("subleadingMuonPt", &subleadingMuonPt, "subleadingMuonPt/F");
    razorTree->Branch("subleadingMuonEta", &subleadingMuonEta, "subleadingMuonEta/F");
    razorTree->Branch("subleadingMuonPhi", &subleadingMuonPhi, "subleadingMuonPhi/F");
    razorTree->Branch("subleadingMuonE", &subleadingMuonE, "subleadingMuonE/F");
    razorTree->Branch("recoZpt", &recoZpt, "recoZpt/F");
    razorTree->Branch("recoZeta", &recoZeta, "recoZeta/F");
    razorTree->Branch("recoZphi", &recoZphi, "recoZphi/F");
    razorTree->Branch("recoZmass", &recoZmass, "recoZmass/F");
    razorTree->Branch("MR_noPho", &MR_noPho, "MR_noPho/F");
    razorTree->Branch("Rsq_noPho", &Rsq_noPho, "Rsq_noPho/F");
    razorTree->Branch("deltaPhi_noPho", &deltaPhi_noPho, "deltaPhi_noPho/F");
    razorTree->Branch("numJets_noPho", &numJets_noPho, "numJets_noPho/I");
    razorTree->Branch("numJets80_noPho", &numJets80_noPho, "numJets80_noPho/I");
    razorTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
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
    if(filterEvents){
        razorTree->Branch("met_noPho", &met_noPho, "met_noPho/F");
        razorTree->Branch("metphi_noPho", &metphi_noPho, "metphi_noPho/F");
        razorTree->Branch("HT_noPho", &HT_noPho, "HT_noPho/F");
        razorTree->Branch("nSelectedJets", &nSelectedJets, "nSelectedJets/I");
        razorTree->Branch("nBTaggedJets", &nBTaggedJets, "nBTaggedJets/I");
        razorTree->Branch("nVetoMuons", &nVetoMuons, "nVetoMuons/I");
        razorTree->Branch("nVetoElectrons", &nVetoElectrons, "nVetoElectrons/I");
        razorTree->Branch("nLooseElectrons", &nLooseElectrons, "nLooseElectrons/I");
        razorTree->Branch("nTightElectrons", &nTightElectrons, "nTightElectrons/I");
        razorTree->Branch("nLooseTaus", &nLooseTaus, "nLooseTaus/I");
        razorTree->Branch("nMediumTaus", &nMediumTaus, "nMediumTaus/I");
        razorTree->Branch("nTightTaus", &nTightTaus, "nTightTaus/I");
        razorTree->Branch("subleadingPhotonPt", &subleadingPhotonPt, "subleadingPhotonPt/F");
        razorTree->Branch("subleadingPhotonEta", &subleadingPhotonEta, "subleadingPhotonEta/F");
        razorTree->Branch("subleadingPhotonPhi", &subleadingPhotonPhi, "subleadingPhotonPhi/F");
        razorTree->Branch("subleadingPhotonE", &subleadingPhotonE, "subleadingPhotonE/F");
        razorTree->Branch("recoWpt", &recoWpt, "recoWpt/F");
        razorTree->Branch("recoWphi", &recoWphi, "recoWphi/F");
        razorTree->Branch("MR", &theMR, "MR/F");
        razorTree->Branch("MR_noZ", &MR_noZ, "MR_noZ/F");
        razorTree->Branch("MR_noW", &MR_noW, "MR_noW/F");
        razorTree->Branch("Rsq", &theRsq, "Rsq/F");
        razorTree->Branch("Rsq_noZ", &Rsq_noZ, "Rsq_noZ/F");
        razorTree->Branch("Rsq_noW", &Rsq_noW, "Rsq_noW/F");
        razorTree->Branch("deltaPhi", &deltaPhi, "deltaPhi/F");
        razorTree->Branch("deltaPhi_noZ", &deltaPhi_noZ, "deltaPhi_noZ/F");
        razorTree->Branch("deltaPhi_noW", &deltaPhi_noW, "deltaPhi_noW/F");
        razorTree->Branch("met", &met, "met/F");
        razorTree->Branch("metphi", &metphi, "metphi/F");
        razorTree->Branch("met_noZ", &met_noZ, "met_noZ/F");
        razorTree->Branch("met_noW", &met_noW, "met_noW/F");
        razorTree->Branch("metphi_noZ", &metphi_noZ, "metphi_noZ/F");
        razorTree->Branch("metphi_noW", &metphi_noW, "metphi_noW/F");
        razorTree->Branch("HT", &HT, "HT/F");
        razorTree->Branch("HT_noZ", &HT_noZ, "HT_noZ/F");
        razorTree->Branch("HT_noW", &HT_noW, "HT_noW/F");
        razorTree->Branch("numJets", &numJets, "numJets/I");
        razorTree->Branch("numJets_noZ", &numJets_noZ, "numJets_noZ/I");
        razorTree->Branch("numJets_noW", &numJets_noW, "numJets_noW/I");
        razorTree->Branch("numJets80", &numJets80, "numJets80/I");
        razorTree->Branch("numJets80_noZ", &numJets80_noZ, "numJets80_noZ/I");
        razorTree->Branch("numJets80_noW", &numJets80_noW, "numJets80_noW/I");
        razorTree->Branch("mTLepMet", &mTLepMet, "mTLepMet/F");
        razorTree->Branch("nBJetsLoose20GeV", &nBJetsLoose20GeV, "nBJetsLoose20GeV/I");
        razorTree->Branch("nBJetsMedium20GeV", &nBJetsMedium20GeV, "nBJetsMedium20GeV/I");
        razorTree->Branch("nBJetsTight20GeV", &nBJetsTight20GeV, "nBJetsTight20GeV/I");
        razorTree->Branch("bjet1PassLoose", &bjet1PassLoose, "bjet1PassLoose/O");
        razorTree->Branch("bjet1PassMedium", &bjet1PassMedium, "bjet1PassMedium/O");
        razorTree->Branch("bjet1PassTight", &bjet1PassTight, "bjet1PassTight/O");
        razorTree->Branch("bjet2PassLoose", &bjet2PassLoose, "bjet2PassLoose/O");
        razorTree->Branch("bjet2PassMedium", &bjet2PassMedium, "bjet2PassMedium/O");
        razorTree->Branch("bjet2PassTight", &bjet2PassTight, "bjet2PassTight/O");
        razorTree->Branch("bjet1Pt", &bjet1Pt, "bjet1Pt/F");
        razorTree->Branch("bjet2Pt", &bjet2Pt, "bjet2Pt/F");
    }
    
    //****************************************************//
    //                  Set up JEC                        //
    //****************************************************//
   

    //get the jet correction parameters
    std::vector<JetCorrectorParameters> correctionParameters;
    char* cmsswPath;
    FactorizedJetCorrector *JetCorrector;
    JetCorrectorParameters *JetResolutionParameters;
    SimpleJetResolution *JetResolutionCalculator;
    cmsswPath = getenv("CMSSW_BASE");
    if(cmsswPath != NULL){
        string pathname(cmsswPath);
        pathname = pathname+"/src/RazorAnalyzer/data/";
        cout << "Getting JEC parameters from " << pathname << endl;
        if (isRunOne) {
            if(isData){
                correctionParameters.push_back(JetCorrectorParameters(Form("%s/Winter14_V8_DATA_L1FastJet_AK5PF.txt", pathname.c_str())));
                correctionParameters.push_back(JetCorrectorParameters(Form("%s/Winter14_V8_DATA_L2Relative_AK5PF.txt", pathname.c_str())));
                correctionParameters.push_back(JetCorrectorParameters(Form("%s/Winter14_V8_DATA_L3Absolute_AK5PF.txt", pathname.c_str()))); 
                correctionParameters.push_back(JetCorrectorParameters(Form("%s/Winter14_V8_DATA_L2L3Residual_AK5PF.txt", pathname.c_str()))); 
            }
            else{
                correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer13_V4_MC_L1FastJet_AK5PF.txt", pathname.c_str())));
                correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer13_V4_MC_L2Relative_AK5PF.txt", pathname.c_str())));
                correctionParameters.push_back(JetCorrectorParameters(Form("%s/Summer13_V4_MC_L3Absolute_AK5PF.txt", pathname.c_str()))); 
            }
        }
        else{ //Run 2
            correctionParameters.push_back(JetCorrectorParameters(Form("%s/PHYS14_V2_MC_L1FastJet_AK4PFchs.txt", pathname.c_str())));
            correctionParameters.push_back(JetCorrectorParameters(Form("%s/PHYS14_V2_MC_L2Relative_AK4PFchs.txt", pathname.c_str())));
            correctionParameters.push_back(JetCorrectorParameters(Form("%s/PHYS14_V2_MC_L3Absolute_AK4PFchs.txt", pathname.c_str())));    
        }
        JetCorrector = new FactorizedJetCorrector(correctionParameters);
        JetResolutionParameters = new JetCorrectorParameters(Form("%s/JetResolutionInputAK5PF.txt", pathname.c_str()));
        JetResolutionCalculator = new SimpleJetResolution(*JetResolutionParameters);
    } 
    else{
        cout << "Error: CMSSW_BASE is not defined!  Please set up CMSSW." << endl <<  "Exiting..." << endl;
        return;
    }

    //****************************************************//
    //            Begin the event loop                    //
    //****************************************************//

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {

        //****************************************************//
        //               Initialize the event                 //
        //****************************************************//
        if(jentry % 1000 == 0) cout << "Processing entry " << jentry << endl;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //fill normalization histogram
        NEvents->Fill(1.0);

        //reset tree variables
        if(!isData){
            nGenMuons = 0;
            nGenElectrons = 0;
            nGenTauMuons = 0;
            nGenTauElectrons = 0;
            nGenTaus = 0;
            nGenPhotons = 0;
            MR_noGenZ = -1;
            genZpt = -1;
            genZeta = -99;
            genZphi = -99;
            genZmass = -1;
            genWpt = -1;
            genWeta = -999;
            genWphi = -999;
            met_noGenZ = -1.;
            metphi_noGenZ = -99.;
            genmet = genMetPt;
            genmetphi = genMetPhi;
            Rsq_noGenZ = -1;
            deltaPhi_noGenZ = -1;
            HT_noGenZ = 0;
            numJets_noGenZ = 0;
            numJets80_noGenZ = 0;
            leadGenMuonIsFound = false;
            leadGenPhotonIsFound = false;
            leadGenMuonIsFoundTight = false;
            ptMatchingLeadGenMuon = -1;
            ptMatchingLeadGenPhoton = -1;
            leadingGenMuonPt = 0;
            leadingGenPhotonPt = 0;
            leadingGenMuonEta = -999;
            leadingGenPhotonEta = -999;
            leadingGenMuonPhi = -999;
            leadingGenPhotonPhi = -999;
            leadingGenMuonE = 0;
            leadingGenPhotonE = 0;
            subleadingGenMuonPt = 0;
            subleadingGenPhotonPt = 0;
            subleadingGenMuonEta = -999;
            subleadingGenPhotonEta = -999;
            subleadingGenMuonPhi = -999;
            subleadingGenPhotonPhi = -999;
            subleadingGenMuonE = 0;
            subleadingGenPhotonE = 0;
        }
        run = 0;
        lumi = 0;
        event = 0;
        nVtx = 0;
        nPU_mean = 0;
        hlt_dimuon = false;
        hlt_singlemu = false;
        hlt_photon = false;
        hlt_razor = false;
        hlt_photon_weight = 0.;
        passedHLTPhoton50 = false;
        passedHLTPhoton75 = false;
        passedHLTPhoton90 = false;
        passedHLTPhoton135 = false;
        passedHLTPhoton150 = false;
        passedHLTPhoton160 = false;
        nSelectedJets = 0;
        nBTaggedJets = 0;
        nVetoMuons = 0;
        nLooseMuons = 0;
        nTightMuons = 0;
        nVetoElectrons = 0;
        nLooseElectrons = 0;
        nTightElectrons = 0;
        nLooseTaus = 0;
        nMediumTaus = 0;
        nTightTaus = 0;
        nSelectedPhotons = 0;
        theMR = -1;
        MR_noZ = -1;
        MR_noW = -1;
        MR_noPho = -1;
        recoZpt = -1;
        recoZeta = -999;
        recoZphi = -999;
        recoZmass = -1;
        recoWpt = -1;
        recoWphi = -999;
        met = metPt;
        met_noPho = -1.;
        met_noZ = -1.;
        met_noW = -1.;
        metphi_noZ = -99.;
        metphi_noW = -99.;
        metphi_noPho = -99.;
        metphi = -99.;
        leadingMuonPt = -1;
        leadingMuonEta = -999;
        leadingMuonPhi = -999;
        leadingMuonE = -999;
        leadingTightMuonPt = -1;
        leadingTightMuonEta = -999;
        leadingTightMuonPhi = -999;
        leadingTightMuonE = -999;
        leadingPhotonPt = -1;
        leadingPhotonEta = -999;
        leadingPhotonPhi = -999;
        leadingPhotonE = -999;
        subleadingMuonPt = -1;
        subleadingMuonEta = -999;
        subleadingMuonPhi = -999;
        subleadingMuonE = -999;
        subleadingPhotonPt = -1;
        subleadingPhotonEta = -999;
        subleadingPhotonPhi = -999;
        subleadingPhotonE = -999;
        theRsq = -1;
        Rsq_noPho = -1;
        Rsq_noZ = -1;
        Rsq_noW = -1;
        deltaPhi = -1;
        deltaPhi_noPho = -1;
        deltaPhi_noZ = -1;
        deltaPhi_noW = -1;
        HT = 0;
        HT_noPho = 0;
        HT_noZ = 0;
        HT_noW = 0;
        numJets = 0;
        numJets_noPho = 0;
        numJets_noZ = 0;
        numJets_noW = 0;
        numJets80 = 0;
        numJets80_noPho = 0;
        numJets80_noZ = 0;
        numJets80_noW = 0;
        mTLepMet = -1;
	nBJetsLoose20GeV = nBJetsMedium20GeV = nBJetsTight20GeV = 0;
	bjet1PassLoose = bjet1PassMedium = bjet1PassTight = false;
	bjet2PassLoose = bjet2PassMedium = bjet2PassTight = false;
	bjet1Pt = bjet2Pt = -999.;
    
        //****************************************************//
        //               Select PU and trigger                //
        //****************************************************//
        nVtx = nPV;
        if(!isData)
            for(int i=0; i<nBunchXing; i++)
                if(BunchXing[i]==0) nPU_mean = nPUmean[i];

        //dimuon trigger
        //3 HLT_Mu17_Mu8
        //4 HLT_Mu17_TkMu8
        if(HLTDecision[3] == 1 || HLTDecision[4] == 1 )
            hlt_dimuon = true;

        //single muon trigger
        //0 HLT_IsoMu24
        //1 HLT_IsoMu24_eta2p1
        if(HLTDecision[0] == 1 || HLTDecision[1] == 1)
            hlt_singlemu = true;

        //photon trigger: 
        //29 HLT_Photon50_CaloIdVL
        //30 HLT_Photon75_CaloIdVL
        //31 HLT_Photon90_CaloIdVL
        //32 HLT_Photon135
        //33 HLT_Photon150
        //34 HLT_Photon160
        for(int i = 29; i <= 34; i++){
            if(HLTDecision[i] == 1) hlt_photon = true;
        }
        //see AN_2013_373 for estimates of the luminosities collected by each trigger
        float lumi_HLTPhoton50  = 1.353e0 + 4.921e0 + 7.947e0 + 8.131e0;
        float lumi_HLTPhoton75  = 8.111e0 + 2.953e1 + 4.768e1 + 4.879e1;
        float lumi_HLTPhoton90  = 1.622e1 + 6.408e1 + 1.010e2 + 9.948e1;
        float lumi_HLTPhoton135 = 8.893e2 + 1.476e2 + 5.429e3 + 7.318e3;
        float lumi_HLTPhoton150 = 8.893e2 + 4.429e3 + 7.152e3 + 7.318e3;
        //data -- scale each event up according to the prescale of the tightest trigger it passed
        //(apply trigger weights EITHER to data OR MC -- not both!)
        if(isData && hlt_photon){
            if(HLTDecision[33] == 1 || HLTDecision[34] == 1){ 
                hlt_photon_weight = 1.0;
            }
            else if(HLTDecision[32] == 1){
                hlt_photon_weight = lumi_HLTPhoton150/lumi_HLTPhoton135;
            }
            else if(HLTDecision[31] == 1){
                hlt_photon_weight = lumi_HLTPhoton150/lumi_HLTPhoton90;
            }
            else if(HLTDecision[30] == 1){
                hlt_photon_weight = lumi_HLTPhoton150/lumi_HLTPhoton75; 
            }
            else if(HLTDecision[29] == 1){
                hlt_photon_weight = lumi_HLTPhoton150/lumi_HLTPhoton50;
            }
            //save the trigger bits
            if(HLTDecision[34] == 1){
                passedHLTPhoton160 = true;
            }
            if(HLTDecision[33] == 1){
                passedHLTPhoton150 = true;
            }
            if(HLTDecision[32] == 1){
                passedHLTPhoton135 = true;
            }
            if(HLTDecision[31] == 1){
                passedHLTPhoton90 = true;
            }
            if(HLTDecision[30] == 1){
                passedHLTPhoton75 = true;
            }
            if(HLTDecision[29] == 1){
                passedHLTPhoton50 = true;
            }
        }
        //MC -- scale each event down according to the prescale of the tightest trigger it passed
        //(apply trigger weights EITHER to data OR  MC -- not both!)
        else if(hlt_photon){
            if(HLTDecision[33] == 1 || HLTDecision[34] == 1){ 
                hlt_photon_weight = 1.0;
            }
            else if(HLTDecision[32] == 1){
                hlt_photon_weight = lumi_HLTPhoton135/lumi_HLTPhoton150;
            }
            else if(HLTDecision[31] == 1){
                hlt_photon_weight = lumi_HLTPhoton90/lumi_HLTPhoton150;
            }
            else if(HLTDecision[30] == 1){
                hlt_photon_weight = lumi_HLTPhoton75/lumi_HLTPhoton150; 
            }
            else if(HLTDecision[29] == 1){
                hlt_photon_weight = lumi_HLTPhoton50/lumi_HLTPhoton150;
            }
            //save the trigger bits
            if(HLTDecision[34] == 1){
                passedHLTPhoton160 = true;
            }
            if(HLTDecision[33] == 1){
                passedHLTPhoton150 = true;
            }
            if(HLTDecision[32] == 1){
                passedHLTPhoton135 = true;
            }
            if(HLTDecision[31] == 1){
                passedHLTPhoton90 = true;
            }
            if(HLTDecision[30] == 1){
                passedHLTPhoton75 = true;
            }
            if(HLTDecision[29] == 1){
                passedHLTPhoton50 = true;
            }
        }

        for(int i = 46; i <= 50; i++){
            if(HLTDecision[i] == 1) hlt_razor = true;
        }

        if(filterEvents && !hlt_dimuon && !hlt_singlemu && !hlt_photon && !hlt_razor) continue;

        run = runNum;
        lumi = lumiNum;
        event = eventNum;

        //****************************************************//
        //               Select gen particles                 //
        //****************************************************//
        vector<TLorentzVector> GoodGenMuons; //for removing gen muons from jet collection later
        if(!isData){
            for(int j = 0; j < nGenParticle; j++){
                //electrons
                if (abs(gParticleId[j]) == 11 && gParticleStatus[j] == 1) {
                    if (  (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23) ) {
                        nGenElectrons++;
                    }
                    if ( abs(gParticleMotherId[j]) == 15) {
                        nGenTauElectrons++;
                    }
                }
                //muons
                if (abs(gParticleId[j]) == 13 && gParticleStatus[j] == 1) {
                    if ( (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23)) {
                        nGenMuons++;
                        TLorentzVector thisGenMuon = makeTLorentzVector(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]);
                        GoodGenMuons.push_back(thisGenMuon);
                        if (gParticlePt[j] > leadingGenMuonPt) {
                            //make leading gen muon into subleading
                            subleadingGenMuonPt = leadingGenMuonPt;
                            subleadingGenMuonEta = leadingGenMuonEta;
                            subleadingGenMuonPhi = leadingGenMuonPhi;
                            subleadingGenMuonE = leadingGenMuonE; 
                            //make this the leading gen muon
                            leadingGenMuonPt = gParticlePt[j];
                            leadingGenMuonEta = gParticleEta[j];
                            leadingGenMuonPhi = gParticlePhi[j];
                            leadingGenMuonE = gParticleE[j];
                        }
                        else if(gParticlePt[j] > subleadingGenMuonPt){
                            //make this the subleading gen muon
                            subleadingGenMuonPt = gParticlePt[j];
                            subleadingGenMuonEta = gParticleEta[j];
                            subleadingGenMuonPhi = gParticlePhi[j];
                            subleadingGenMuonE = gParticleE[j];
                        }
                    }
                    if ( abs(gParticleMotherId[j]) == 15) {
                        nGenTauMuons++;
                    }
                }
                //taus
                if (abs(gParticleId[j]) == 15 && gParticleStatus[j] == 2 
                        && (abs(gParticleMotherId[j]) == 24 || abs(gParticleMotherId[j]) == 23)
                   ) nGenTaus++;
                //photons
                if (abs(gParticleId[j]) == 22 && 
                        ( (abs(gParticleMotherId[j]) >= 1 && abs(gParticleMotherId[j]) <= 5) || 
                          (abs(gParticleMotherId[j]) == 21) || (abs(gParticleMotherId[j]) == 2212) ) && 
                        abs(gParticleStatus[j]) == 1){	     
                    nGenPhotons++; 
                    if(gParticlePt[j] > leadingGenPhotonPt){
                        //make leading gen photon into subleading
                        subleadingGenPhotonPt = leadingGenPhotonPt;
                        subleadingGenPhotonEta = leadingGenPhotonEta;
                        subleadingGenPhotonPhi = leadingGenPhotonPhi;
                        subleadingGenPhotonE = leadingGenPhotonE;
                        //make this the leading gen photon
                        leadingGenPhotonPt = gParticlePt[j];
                        leadingGenPhotonEta = gParticleEta[j];
                        leadingGenPhotonPhi = gParticlePhi[j];
                        leadingGenPhotonE = gParticleE[j];
                    }
                    else if(gParticlePt[j] > subleadingGenPhotonPt){
                        //make this the subleading photon
                        subleadingGenPhotonPt = gParticlePt[j];
                        subleadingGenPhotonEta = gParticleEta[j];
                        subleadingGenPhotonPhi = gParticlePhi[j];
                        subleadingGenPhotonE = gParticleE[j];
                    }
                }
            }

            // gen level Z pt
            for(int j = 0; j < nGenParticle; j++){
                if(gParticleStatus[j] != 22 && gParticleStatus[j] != 3) continue; //gen-level Z and W have pythia8 status 22 in Run 2 ntuples, status 3 in Run 1 ntuples
                TLorentzVector boson = makeTLorentzVector(gParticlePt[j], gParticleEta[j], gParticlePhi[j], gParticleE[j]); 
                if(abs(gParticleId[j]) == 23){ //Z boson
                    genZpt = gParticlePt[j];
                    genZeta = gParticleEta[j];
                    genZphi = gParticlePhi[j];
                    genZmass = boson.M();
                }
                else if(abs(gParticleId[j]) == 24){ //W boson
                    genWpt = gParticlePt[j];
                    genWeta = gParticleEta[j];
                    genWphi = gParticlePhi[j];
                }
            }
        }

        //****************************************************//
        //               Select muons                         //
        //****************************************************//
        vector<TLorentzVector> GoodMuons; 
        vector<TLorentzVector> GoodMuonsTight;
        for(int i = 0; i < nMuons; i++){

            if(!isLooseMuon(i)) continue;
            if(muonPt[i] < 10) continue;
            if(abs(muonEta[i]) > 2.4) continue;
            TLorentzVector thisMuon = makeTLorentzVector(muonPt[i], muonEta[i], muonPhi[i], muonE[i]); 

            if(isVetoMuon(i)) nVetoMuons++;
            if(isTightMuon(i)){
                nTightMuons++;
                GoodMuonsTight.push_back(thisMuon);

                //check if leading tight muon
                if(muonPt[i] > leadingTightMuonPt){
                    leadingTightMuonPt = muonPt[i];
                    leadingTightMuonEta = muonEta[i];
                    leadingTightMuonPhi = muonPhi[i];
                    leadingTightMuonE = muonE[i];
                }
            }
            nLooseMuons++;

            GoodMuons.push_back(thisMuon);

            //check if this muon is leading or subleading
            if(muonPt[i] > leadingMuonPt){
                //make leading muon into subleading
                subleadingMuonPt = leadingMuonPt;
                subleadingMuonEta = leadingMuonEta;
                subleadingMuonPhi = leadingMuonPhi;
                subleadingMuonE = leadingMuonE;
                //make this the leading muon
                leadingMuonPt = muonPt[i];
                leadingMuonEta = muonEta[i];
                leadingMuonPhi = muonPhi[i];
                leadingMuonE = muonE[i];
            }
            else if(muonPt[i] > subleadingMuonPt){
                //make this the subleading muon
                subleadingMuonPt = muonPt[i];
                subleadingMuonEta = muonEta[i];
                subleadingMuonPhi = muonPhi[i];
                subleadingMuonE = muonE[i];
            }

        }

        //****************************************************//
        //               Select electrons                     //
        //****************************************************//
        for(int i = 0; i < nElectrons; i++){
            if(elePt[i] < 10) continue;
            if(fabs(eleEta[i]) > 2.5) continue;
            if(isMVANonTrigVetoElectron(i)) nVetoElectrons++;

	    if(!isLooseElectron(i)) continue;
	    if(isTightElectron(i)) nTightElectrons++;

            TLorentzVector thisElectron = makeTLorentzVector(elePt[i], eleEta[i], elePhi[i], eleE[i]);
            nLooseElectrons++;
        }

        //****************************************************//
        //               Select taus                          //
        //****************************************************//
        for(int i = 0; i < nTaus; i++){

            if(isLooseTau(i)){
                nLooseTaus++;
            }
            if(isMediumTau(i)){
                nMediumTaus++;
            }
            if(isTightTau(i)){
                nTightTaus++;
            }

        }

        //****************************************************//
        //               Select jets                          //
        //****************************************************//
	bool bjet1Found = false;
	bool bjet2Found = false;

        vector<TLorentzVector> GoodJets; //will contain leptons above 40 GeV in addition to jets
        TVector3 metCorrection; //contains p_T - p_T_JEC summed over all jets
        for(int i = 0; i < nJets; i++){
            //apply JEC
            double JEC = 1.0;
	    JEC = JetEnergyCorrectionFactor(jetPt[i], jetEta[i], jetPhi[i], jetE[i], fixedGridRhoFastjetAll, jetJetArea[i], JetCorrector);   

            //apply jet resolution smearing to MC
            double jetEnergySmearFactor = 1.0;
            if(!isData){
                jetEnergySmearFactor = JetEnergySmearingFactor(jetPt[i]*JEC, jetEta[i], nPU_mean, JetResolutionCalculator, random);
            }
            //check for noise
	    if(!jetPassIDTight[i]) continue;
            
            //jet TLorentzVector, with corrections included
            TLorentzVector thisJet = makeTLorentzVector(jetPt[i]*JEC*jetEnergySmearFactor, jetEta[i], jetPhi[i], jetE[i]*JEC*jetEnergySmearFactor);

            //if the corrected pt is above 20 GeV, this jet is used for the MET correction
            TVector3 unCorrJetPerp, corrJetPerp;
            unCorrJetPerp.SetPtEtaPhi(jetPt[i], 0, jetPhi[i]);
            corrJetPerp.SetPtEtaPhi(jetPt[i]*JEC*jetEnergySmearFactor, 0, jetPhi[i]);
            //propagate the correction to the MET
            if(jetPt[i]*JEC*jetEnergySmearFactor > 20){
                metCorrection = metCorrection + unCorrJetPerp - corrJetPerp;
            }

            //pt and eta cuts
            if(jetPt[i]*JEC*jetEnergySmearFactor < 40) continue;
            if(fabs(jetEta[i]) > 3.0) continue;

            //apply jet PU ID
	    if(isCSVM(i)){ 
	      nBTaggedJets++;
	    }
 
	    // bjets in MC
	    if(!isData)
	      if (abs(jetPartonFlavor[i]) == 5) {
		if (!bjet1Found) {
		  bjet1Found = true;
		  bjet1Pt = jetPt[i];
		  bjet1PassLoose  = false;
		  bjet1PassMedium = false;
		  bjet1PassTight  = false;
		  if(isCSVL(i)) bjet1PassLoose = true;	      
		  if(isCSVM(i)) bjet1PassMedium = true;
		  if(isCSVT(i)) bjet1PassTight = true;	      
		} else if ( jetPt[i] > bjet1Pt ) {		  
		  bjet2Pt = bjet1Pt;
		  bjet2PassLoose  = bjet1PassLoose;
		  bjet2PassMedium = bjet1PassMedium;
		  bjet2PassTight  = bjet1PassTight;
		  
		  bjet1Pt = jetPt[i];
		  bjet1PassLoose  = false;
		  bjet1PassMedium = false;
		  bjet1PassTight  = false;
		  if(isCSVL(i)) bjet1PassLoose = true;	      
		  if(isCSVM(i)) bjet1PassMedium = true;
		  if(isCSVT(i)) bjet1PassTight = true;	      
		} else {
		  if (!bjet2Found || jetPt[i] > bjet2Pt ) {
		    bjet2Found = true;
		    bjet2Pt = jetPt[i];
		    bjet2PassLoose  = false;
		    bjet2PassMedium = false;
		    bjet2PassTight  = false;
		    if(isCSVL(i)) bjet2PassLoose = true;	      
		    if(isCSVM(i)) bjet2PassMedium = true;
		    if(isCSVT(i)) bjet2PassTight = true;	   	    
		  }
		}
	      } //if it's a bjet		
	    

            if(jetPt[i]*JEC*jetEnergySmearFactor > 80) numJets80++;
            GoodJets.push_back(thisJet);
            nSelectedJets++;

        }
        sort(GoodJets.begin(), GoodJets.end(), greater_than_pt());

        //****************************************************//
        //     Compute the razor variables and HT, nJets      //
        //****************************************************//

        //get the type-1 corrected PFMET
        TVector3 metUncorr;
        metUncorr.SetPtEtaPhi(metPt, 0, metPhi);
        TVector3 metCorr = metUncorr + metCorrection;

        TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metCorr.Pt(), 0, metCorr.Phi(), 0);
        metphi = PFMET.Phi();

        //count jets and compute HT
        numJets = GoodJets.size();
        for(auto& pf : GoodJets) HT += pf.Pt();

        // compute R and MR
        if(GoodJets.size() >= 2 && GoodJets.size() < 20){
            vector<TLorentzVector> hemispheres = getHemispheres(GoodJets);
            theMR = computeMR(hemispheres[0], hemispheres[1]); 
            theRsq = computeRsq(hemispheres[0], hemispheres[1], PFMET);
            deltaPhi = fabs(hemispheres[0].DeltaPhi(hemispheres[1]));
        }

        //****************************************************//
        //             Select photons                         //
        //****************************************************//
        vector<TLorentzVector> GoodPhotons;
        int nPhotonsAbove40GeV = 0;
        for(int i = 0; i < nPhotons; i++){
            if(phoPt[i] < 10) continue;
            if(fabs(phoEta[i]) > 2.5) continue;
	    if(!isTightPhoton(i)) continue;

            if(phoPt[i] > 40) nPhotonsAbove40GeV++;
            TLorentzVector thisPhoton = makeTLorentzVector(phoPt[i], phoEta[i], phoPhi[i], pho_RegressionE[i]);
            GoodPhotons.push_back(thisPhoton);
            nSelectedPhotons++;

            //check if this photon is leading or subleading
            if(phoPt[i] > leadingPhotonPt){
                //make leading photon into subleading
                subleadingPhotonPt = leadingPhotonPt;
                subleadingPhotonEta = leadingPhotonEta;
                subleadingPhotonPhi = leadingPhotonPhi;
                subleadingPhotonE = leadingPhotonE;
                //make this the leading photon
                leadingPhotonPt = phoPt[i];
                leadingPhotonEta = phoEta[i];
                leadingPhotonPhi = phoPhi[i];
                leadingPhotonE = phoE[i];
            }
            else if(phoPt[i] > subleadingPhotonPt){
                //make this the subleading photon
                subleadingPhotonPt = phoPt[i];
                subleadingPhotonEta = phoEta[i];
                subleadingPhotonPhi = phoPhi[i];
                subleadingPhotonE = phoE[i];
            }
        }

        //****************************************************//
        //        Match gen-level and reco objects            //
        //****************************************************//
        if(!isData){
            if(nGenMuons > 0){
                //see if we selected a muon matching the leading gen muon
                TLorentzVector leadGenMuon = makeTLorentzVector(leadingGenMuonPt, leadingGenMuonEta, leadingGenMuonPhi, leadingGenMuonE);
                for(auto& mu : GoodMuons){
                    float thisDeltaR = mu.DeltaR(leadGenMuon);
                    if(thisDeltaR < 0.1){ //muon matches leading gen muon
                        leadGenMuonIsFound = true;
                        ptMatchingLeadGenMuon = mu.Pt(); //fill pt of the muon matching the gen muon
                        break;
                    }
                }
                for(auto& mu : GoodMuonsTight){
                    float thisDeltaR = mu.DeltaR(leadGenMuon);
                    if(thisDeltaR < 0.1){ //muon matches leading gen muon
                        leadGenMuonIsFoundTight = true;
                        break;
                    }
                }
            }
            if(nGenPhotons > 0){
                //see if we selected a photon matching the leading gen photon
                TLorentzVector leadGenPhoton = makeTLorentzVector(leadingGenPhotonPt, leadingGenPhotonEta, leadingGenPhotonPhi, leadingGenPhotonE);
                for(auto& pho : GoodPhotons){
                    float thisDeltaR = pho.DeltaR(leadGenPhoton);
                    if(thisDeltaR < 0.1){ //photon matches leading gen photon
                        leadGenPhotonIsFound = true;
                        ptMatchingLeadGenPhoton = pho.Pt(); //fill pt of the photon matching the gen photon
                        break;
                    }
                }
            }
        }

        //****************************************************//
        //    Compute razor vars for DY, W, Gamma samples     //
        //****************************************************//
        //photons
        if(GoodPhotons.size()>0){
            sort(GoodPhotons.begin(), GoodPhotons.end(), greater_than_pt());

            //compute MET with leading photon added
            TLorentzVector m1 = GoodPhotons[0];
            TLorentzVector m2 = PFMET;
            TLorentzVector photonPlusMet_perp = makeTLorentzVectorPtEtaPhiM((m1 + m2).Pt(), 0., (m1 + m2).Phi(), 0.0);

            met_noPho = photonPlusMet_perp.Pt();
            metphi_noPho = photonPlusMet_perp.Phi();

            //remove leading photon from collection of selected jets
            vector<TLorentzVector> GoodJetsNoLeadPhoton = GoodJets;
            int subtractedIndex = SubtractParticleFromCollection(GoodPhotons[0], GoodJetsNoLeadPhoton);
            if(subtractedIndex >= 0){
                if(GoodJetsNoLeadPhoton[subtractedIndex].Pt() < 40){ //erase this jet
                    GoodJetsNoLeadPhoton.erase(GoodJetsNoLeadPhoton.begin()+subtractedIndex);
                }
            }
            //count the number of jets above 80 GeV now
            for(auto& jet : GoodJetsNoLeadPhoton){
                if(jet.Pt() > 80) numJets80_noPho++;
            }

            //count jets and compute HT
            numJets_noPho = GoodJetsNoLeadPhoton.size();
            for(auto& pf : GoodJetsNoLeadPhoton) HT_noPho += pf.Pt();

            if(GoodJetsNoLeadPhoton.size() >= 2 && GoodJetsNoLeadPhoton.size() <20){
                //remake the hemispheres using the new jet collection
                vector<TLorentzVector> hemispheresNoLeadPhoton = getHemispheres(GoodJetsNoLeadPhoton);
                TLorentzVector PFMET_NOPHO = makeTLorentzVectorPtEtaPhiM(met_noPho, 0, metphi_noPho, 0);
                MR_noPho = computeMR(hemispheresNoLeadPhoton[0], hemispheresNoLeadPhoton[1]); 
                Rsq_noPho = computeRsq(hemispheresNoLeadPhoton[0], hemispheresNoLeadPhoton[1], PFMET_NOPHO);
                deltaPhi_noPho = fabs(hemispheresNoLeadPhoton[0].DeltaPhi(hemispheresNoLeadPhoton[1]));
            }
        } 
        else{ //save some info even if no photons are found
            numJets_noPho = numJets;
            numJets80_noPho = numJets80;
            met_noPho = met;
            metphi_noPho = metphi;
            HT_noPho = HT;
        }

        // Muons for Z

        //remove selected muons from collection of selected jets and add them to the MET
        vector<TLorentzVector> GoodJetsNoMuons = GoodJets;
        TLorentzVector TotalMuonVec;
        for(auto& mu : GoodMuons){
            TotalMuonVec = TotalMuonVec + mu; //add this muon's momentum to the sum
            int subtractedIndex = SubtractParticleFromCollection(mu, GoodJetsNoMuons);
            if(subtractedIndex >= 0){
                if(GoodJetsNoMuons[subtractedIndex].Pt() < 40){ //erase this jet
                    GoodJetsNoMuons.erase(GoodJetsNoMuons.begin()+subtractedIndex);
                }
            }
        }
        //remove selected TIGHT muons from collection of selected jets and add them to the MET
        vector<TLorentzVector> GoodJetsNoTightMuons = GoodJets;
        TLorentzVector TotalTightMuonVec;
        for(auto& mu : GoodMuonsTight){
            TotalTightMuonVec = TotalTightMuonVec + mu; //add this muon's momentum to the sum
            int subtractedIndex = SubtractParticleFromCollection(mu, GoodJetsNoTightMuons);
            if(subtractedIndex >= 0){
                if(GoodJetsNoTightMuons[subtractedIndex].Pt() < 40){ //erase this jet
                    GoodJetsNoTightMuons.erase(GoodJetsNoTightMuons.begin()+subtractedIndex);
                }
            }
        }

        //do the same for GEN muons
        vector<TLorentzVector> GoodJetsNoGenMuons = GoodJets;
        TLorentzVector TotalGenMuonVec;
        if(!isData){
            for(auto& mu : GoodGenMuons){
                TotalGenMuonVec = TotalGenMuonVec + mu;
                int subtractedIndex = SubtractParticleFromCollection(mu, GoodJetsNoGenMuons);
                if(subtractedIndex >= 0){
                    if(GoodJetsNoGenMuons[subtractedIndex].Pt() < 40){ //erase this jet
                        GoodJetsNoGenMuons.erase(GoodJetsNoGenMuons.begin()+subtractedIndex);
                    }
                }
            }
        }

        //make the MET vector with the muons (or gen muons) added
        TLorentzVector ZPlusMet_perp = makeTLorentzVector((TotalMuonVec + PFMET).Pt(), 0., (TotalMuonVec + PFMET).Phi(), 0.);
        met_noZ = ZPlusMet_perp.Pt();
        metphi_noZ = ZPlusMet_perp.Phi();

        TLorentzVector WPlusMet_perp = makeTLorentzVector((TotalTightMuonVec + PFMET).Pt(), 0., (TotalTightMuonVec + PFMET).Phi(), 0.);
        met_noW = WPlusMet_perp.Pt();
        metphi_noW = WPlusMet_perp.Phi(); 

        TLorentzVector ZPlusMetGen_perp;
        if(!isData){
            ZPlusMetGen_perp = makeTLorentzVectorPtEtaPhiM((TotalGenMuonVec + PFMET).Pt(), 0., (TotalGenMuonVec + PFMET).Phi(), 0.);
            met_noGenZ = ZPlusMetGen_perp.Pt();
            metphi_noGenZ = ZPlusMetGen_perp.Phi();
        }

        //count jets and compute HT
        //Z
        numJets_noZ = GoodJetsNoMuons.size();
        for(auto& jet : GoodJetsNoMuons){
            HT_noZ += jet.Pt();
            if(jet.Pt() > 80) numJets80_noZ++;
        }
        //W
        numJets_noW = GoodJetsNoTightMuons.size();
        for(auto& jet : GoodJetsNoTightMuons){
            HT_noW += jet.Pt();
            if(jet.Pt() > 80) numJets80_noW++;
        }

        //get reco Z information
        if(GoodMuons.size() >= 1){
            recoZpt = TotalMuonVec.Pt();
            recoZeta = TotalMuonVec.Eta();
            recoZphi = TotalMuonVec.Phi();
            recoZmass = TotalMuonVec.M();
        }

        //compute reco Z information and razor variables for DY
        if(numJets_noZ > 1 && GoodJets.size()<20)
        {
            vector<TLorentzVector> hemispheresNoZ = getHemispheres(GoodJetsNoMuons);
            Rsq_noZ = computeRsq(hemispheresNoZ[0], hemispheresNoZ[1], ZPlusMet_perp);
            MR_noZ = computeMR(hemispheresNoZ[0], hemispheresNoZ[1]); 
            deltaPhi_noZ = fabs(hemispheresNoZ[0].DeltaPhi(hemispheresNoZ[1])); 
        }
        if(!isData){
            //Gen Z
            numJets_noGenZ = GoodJetsNoGenMuons.size();
            for(auto& jet : GoodJetsNoGenMuons){
                HT_noGenZ += jet.Pt();
                if(jet.Pt() > 80) numJets80_noGenZ++;
            }
            //razor variables using GEN muons
            if(numJets_noGenZ > 1 && GoodJets.size()<20)
            {
                vector<TLorentzVector> hemispheresNoGenZ = getHemispheres(GoodJetsNoGenMuons);
                Rsq_noGenZ = computeRsq(hemispheresNoGenZ[0], hemispheresNoGenZ[1], ZPlusMetGen_perp);
                MR_noGenZ = computeMR(hemispheresNoGenZ[0], hemispheresNoGenZ[1]); 
                deltaPhi_noGenZ = fabs(hemispheresNoGenZ[0].DeltaPhi(hemispheresNoGenZ[1])); 
            }
        }
        //razor variables using tight muons (for W)
        if(numJets_noW > 1 && GoodJets.size()<20){
            vector<TLorentzVector> hemispheresNoW = getHemispheres(GoodJetsNoTightMuons);
            Rsq_noW = computeRsq(hemispheresNoW[0], hemispheresNoW[1], WPlusMet_perp);
            MR_noW = computeMR(hemispheresNoW[0], hemispheresNoW[1]); 
            deltaPhi_noW = fabs(hemispheresNoW[0].DeltaPhi(hemispheresNoW[1])); 
        }

        //for W, also get the transverse mass of the first tight muon and the MET
        if(GoodMuonsTight.size() > 0) 
        {
            TLorentzVector m1 = GoodMuonsTight[0];
            TLorentzVector m2 = PFMET;
            double deltaPhiLepMet = m1.DeltaPhi(m2);
            mTLepMet = sqrt(2*m2.Pt()*m1.Pt()*( 1.0 - cos( deltaPhiLepMet ) ) ); //transverse mass calculation

            //store reco W information
            recoWpt = (m1+m2).Pt();
            recoWphi = (m1+m2).Phi();
        }

        //************************//
        //*****Filter events******//
        //************************//
        if(filterEvents){
	  // if(nSelectedPhotons < 1 ) continue;
            if(numJets80 < 2) continue; //event fails to have two 80 GeV jets
            // if(GoodMuons.size() == 0 && GoodPhotons.size() == 0) continue; //don't save event if no muons or photons
            if(theMR < 300 && MR_noZ < 300 && MR_noW < 300 && MR_noPho < 300) continue;
            if(theRsq < 0.15 && Rsq_noZ < 0.15 && Rsq_noW < 0.15 && Rsq_noPho < 0.15) continue;
        }

        razorTree->Fill();
    }//end of event loop

    cout << "Writing output tree..." << endl;
    razorTree->Write();
    NEvents->Write();

    outFile.Close();
}


