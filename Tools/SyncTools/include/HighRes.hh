//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May  8 17:13:27 2015 by ROOT version 6.03/03
// from TTree HighRes/HighRes
// found on file: RazorHggCat_Total.root
//////////////////////////////////////////////////////////

#ifndef HighRes_hh
#define HighRes_hh 1

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class HighRes {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           nLooseBTaggedJets;
   Int_t           nMediumBTaggedJets;
   Int_t           nLooseMuons;
   Int_t           nTightMuons;
   Int_t           nLooseElectrons;
   Int_t           nTightElectrons;
   Int_t           nTightTaus;
   Float_t         MR;
   Float_t         Rsq;
   Float_t         t1Rsq;
   Float_t         MET;
   Float_t         t1MET;
   Int_t           nSelectedPhotons;
   Float_t         mGammaGamma;
   Float_t         pTGammaGamma;
   Float_t         pho1E;
   Float_t         pho1Pt;
   Float_t         Pho1Eta;
   Float_t         pho1Phi;
   Float_t         pho1SigmaIetaIeta;
   Float_t         pho1R9;
   Float_t         pho1HoverE;
   Float_t         pho1sumChargedHadronPt;
   Float_t         pho1sumNeutralHadronEt;
   Float_t         pho1sumPhotonEt;
   Float_t         pho1sigmaEOverE;
   Bool_t          pho1passEleVeto;
   Bool_t          pho1passIso;
   Float_t         pho2E;
   Float_t         pho2Pt;
   Float_t         Pho2Eta;
   Float_t         pho2Phi;
   Float_t         pho2SigmaIetaIeta;
   Float_t         pho2R9;
   Float_t         pho2HoverE;
   Float_t         pho2sumChargedHadronPt;
   Float_t         pho2sumNeutralHadronEt;
   Float_t         pho2sumPhotonEt;
   Float_t         pho2sigmaEOverE;
   Bool_t          pho2passEleVeto;
   Bool_t          pho2passIso;
   Float_t         mbbZ;
   Float_t         mbbH;
   Int_t           n_Jets;
   Float_t         jet_E[9];   //[n_Jets]
   Float_t         jet_Pt[9];   //[n_Jets]
   Float_t         jet_Eta[9];   //[n_Jets]
   Float_t         jet_Phi[9];   //[n_Jets]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nLooseBTaggedJets;   //!
   TBranch        *b_nMediumBTaggedJets;   //!
   TBranch        *b_nLooseMuons;   //!
   TBranch        *b_nTightMuons;   //!
   TBranch        *b_nLooseElectrons;   //!
   TBranch        *b_nTightElectrons;   //!
   TBranch        *b_nTightTaus;   //!
   TBranch        *b_MR;   //!
   TBranch        *b_Rsq;   //!
   TBranch        *b_t1Rsq;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_t1MET;   //!
   TBranch        *b_nSelectedPhotons;   //!
   TBranch        *b_mGammaGamma;   //!
   TBranch        *b_pTGammaGamma;   //!
   TBranch        *b_pho1E;   //!
   TBranch        *b_pho1Pt;   //!
   TBranch        *b_pho1Eta;   //!
   TBranch        *b_pho1Phi;   //!
   TBranch        *b_pho1SigmaIetaIeta;   //!
   TBranch        *b_pho1R9;   //!
   TBranch        *b_pho1HoverE;   //!
   TBranch        *b_pho1sumChargedHadronPt;   //!
   TBranch        *b_pho1sumNeutralHadronEt;   //!
   TBranch        *b_pho1sumPhotonEt;   //!
   TBranch        *b_pho1sigmaEOverE;   //!
   TBranch        *b_pho1passEleVeto;   //!
   TBranch        *b_pho1passIso;   //!
   TBranch        *b_pho2E;   //!
   TBranch        *b_pho2Pt;   //!
   TBranch        *b_pho2Eta;   //!
   TBranch        *b_pho2Phi;   //!
   TBranch        *b_pho2SigmaIetaIeta;   //!
   TBranch        *b_pho2R9;   //!
   TBranch        *b_pho2HoverE;   //!
   TBranch        *b_pho2sumChargedHadronPt;   //!
   TBranch        *b_pho2sumNeutralHadronEt;   //!
   TBranch        *b_pho2sumPhotonEt;   //!
   TBranch        *b_pho2sigmaEOverE;   //!
   TBranch        *b_pho2passEleVeto;   //!
   TBranch        *b_pho2passIso;   //!
   TBranch        *b_mbbZ;   //!
   TBranch        *b_mbbH;   //!
   TBranch        *b_n_Jets;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_Pt;   //!
   TBranch        *b_jet_Eta;   //!
   TBranch        *b_jet_Phi;   //!

   HighRes(TTree *tree=0);
   virtual ~HighRes();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HighRes_cxx
HighRes::HighRes(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RazorHggCat_Total.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RazorHggCat_Total.root");
      }
      f->GetObject("HighRes",tree);

   }
   Init(tree);
}

HighRes::~HighRes()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HighRes::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HighRes::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HighRes::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nLooseBTaggedJets", &nLooseBTaggedJets, &b_nLooseBTaggedJets);
   fChain->SetBranchAddress("nMediumBTaggedJets", &nMediumBTaggedJets, &b_nMediumBTaggedJets);
   fChain->SetBranchAddress("nLooseMuons", &nLooseMuons, &b_nLooseMuons);
   fChain->SetBranchAddress("nTightMuons", &nTightMuons, &b_nTightMuons);
   fChain->SetBranchAddress("nLooseElectrons", &nLooseElectrons, &b_nLooseElectrons);
   fChain->SetBranchAddress("nTightElectrons", &nTightElectrons, &b_nTightElectrons);
   fChain->SetBranchAddress("nTightTaus", &nTightTaus, &b_nTightTaus);
   fChain->SetBranchAddress("MR", &MR, &b_MR);
   fChain->SetBranchAddress("Rsq", &Rsq, &b_Rsq);
   fChain->SetBranchAddress("t1Rsq", &t1Rsq, &b_t1Rsq);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("t1MET", &t1MET, &b_t1MET);
   fChain->SetBranchAddress("nSelectedPhotons", &nSelectedPhotons, &b_nSelectedPhotons);
   fChain->SetBranchAddress("mGammaGamma", &mGammaGamma, &b_mGammaGamma);
   fChain->SetBranchAddress("pTGammaGamma", &pTGammaGamma, &b_pTGammaGamma);
   fChain->SetBranchAddress("pho1E", &pho1E, &b_pho1E);
   fChain->SetBranchAddress("pho1Pt", &pho1Pt, &b_pho1Pt);
   fChain->SetBranchAddress("Pho1Eta", &Pho1Eta, &b_pho1Eta);
   fChain->SetBranchAddress("pho1Phi", &pho1Phi, &b_pho1Phi);
   fChain->SetBranchAddress("pho1SigmaIetaIeta", &pho1SigmaIetaIeta, &b_pho1SigmaIetaIeta);
   fChain->SetBranchAddress("pho1R9", &pho1R9, &b_pho1R9);
   fChain->SetBranchAddress("pho1HoverE", &pho1HoverE, &b_pho1HoverE);
   fChain->SetBranchAddress("pho1sumChargedHadronPt", &pho1sumChargedHadronPt, &b_pho1sumChargedHadronPt);
   fChain->SetBranchAddress("pho1sumNeutralHadronEt", &pho1sumNeutralHadronEt, &b_pho1sumNeutralHadronEt);
   fChain->SetBranchAddress("pho1sumPhotonEt", &pho1sumPhotonEt, &b_pho1sumPhotonEt);
   fChain->SetBranchAddress("pho1sigmaEOverE", &pho1sigmaEOverE, &b_pho1sigmaEOverE);
   fChain->SetBranchAddress("pho1passEleVeto", &pho1passEleVeto, &b_pho1passEleVeto);
   fChain->SetBranchAddress("pho1passIso", &pho1passIso, &b_pho1passIso);
   fChain->SetBranchAddress("pho2E", &pho2E, &b_pho2E);
   fChain->SetBranchAddress("pho2Pt", &pho2Pt, &b_pho2Pt);
   fChain->SetBranchAddress("Pho2Eta", &Pho2Eta, &b_pho2Eta);
   fChain->SetBranchAddress("pho2Phi", &pho2Phi, &b_pho2Phi);
   fChain->SetBranchAddress("pho2SigmaIetaIeta", &pho2SigmaIetaIeta, &b_pho2SigmaIetaIeta);
   fChain->SetBranchAddress("pho2R9", &pho2R9, &b_pho2R9);
   fChain->SetBranchAddress("pho2HoverE", &pho2HoverE, &b_pho2HoverE);
   fChain->SetBranchAddress("pho2sumChargedHadronPt", &pho2sumChargedHadronPt, &b_pho2sumChargedHadronPt);
   fChain->SetBranchAddress("pho2sumNeutralHadronEt", &pho2sumNeutralHadronEt, &b_pho2sumNeutralHadronEt);
   fChain->SetBranchAddress("pho2sumPhotonEt", &pho2sumPhotonEt, &b_pho2sumPhotonEt);
   fChain->SetBranchAddress("pho2sigmaEOverE", &pho2sigmaEOverE, &b_pho2sigmaEOverE);
   fChain->SetBranchAddress("pho2passEleVeto", &pho2passEleVeto, &b_pho2passEleVeto);
   fChain->SetBranchAddress("pho2passIso", &pho2passIso, &b_pho2passIso);
   fChain->SetBranchAddress("mbbZ", &mbbZ, &b_mbbZ);
   fChain->SetBranchAddress("mbbH", &mbbH, &b_mbbH);
   fChain->SetBranchAddress("n_Jets", &n_Jets, &b_n_Jets);
   fChain->SetBranchAddress("jet_E", jet_E, &b_jet_E);
   fChain->SetBranchAddress("jet_Pt", jet_Pt, &b_jet_Pt);
   fChain->SetBranchAddress("jet_Eta", jet_Eta, &b_jet_Eta);
   fChain->SetBranchAddress("jet_Phi", jet_Phi, &b_jet_Phi);
   Notify();
}

Bool_t HighRes::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HighRes::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HighRes::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HighRes_cxx
