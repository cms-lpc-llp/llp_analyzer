#ifndef TagAndProbePair_H
#define TagAndProbePair_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"
#include <Rtypes.h>
#include "TLorentzVector.h"
#include "TH1F.h"

class TagAndProbePair {
  
 public:
  
  /// bit map
  /// DON'T CHANGE ORDER
  enum TreeType { kTPType_Default = 0,	                    // default		 
		  kTPType_Full = 1,                         // full
  };
  

  /// variables
  Bool_t                  pass;
  Float_t                 weight;
  UInt_t                  run;
  UInt_t                  lumi;
  UInt_t                  event;
  UInt_t                  NPU_0;
  UInt_t                  NPU_Minus1;
  UInt_t                  NPU_Plus1;
  UInt_t                  NPV;
  Float_t                 Rho;
  Float_t                 pt;
  Float_t                 eta;
  Float_t                 phi;
  Float_t                 mass;
  Int_t                   charge;
  Float_t                 miniIso;
  Float_t                 met;
  Bool_t                  passTighterTag;
  Float_t                 Activity;
  Float_t                 HT;
  Float_t                 MinDRToParton;
  UInt_t                  NJets;


 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  
  /// default constructor  
  TagAndProbePair()  {
    InitVariables();
  };

  /// default destructor
  ~TagAndProbePair(){ 
    if (f_) f_->Close();  
  };
    
  /// initialize varibles and fill list of available variables
  void InitVariables() {
    pass                  = false;
    weight                = 0.0;
    run                   = 0;
    lumi                  = 0;
    event                 = 0;
    NPU_0                 = 0;
    NPU_Minus1            = 0;
    NPU_Plus1             = 0;
    NPV                   = 0;
    Rho                   = 0;
    pt                    = 0;
    eta                   = 0;
    phi                   = 0;
    mass                  = 0;
    charge                = 0;
    miniIso               = 0;
    met                   = 0;
    passTighterTag        = false;
    Activity              = 0;
    HT                    = 0;
    MinDRToParton         = 0;
    NJets                 = 0;
  }
    
  /// load a TagAndProbePair
  void LoadTree(const char* file, int treeType = kTPType_Default){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->Get("TP"));
    InitTree(treeType);
    assert(tree_);
  }
    
  /// create a TagAndProbePair
  void CreateTree(int treeType = kTPType_Default){
    tree_ = new TTree("TP","TP");
    f_ = 0;

    //book the branches that go in all types of trees
    tree_->Branch("pass",&pass,"pass/O");
    tree_->Branch("weight",&weight,"weight/F");   
    tree_->Branch("run",&run,"run/i");           
    tree_->Branch("lumi",&lumi,"lumi/i");          
    tree_->Branch("event",&event,"event/i");         
    tree_->Branch("NPU_0",&NPU_0,"NPU_0/i");         
    tree_->Branch("NPU_Minus1",&NPU_Minus1,"NPU_Minus1/i");    
    tree_->Branch("NPU_Plus1",&NPU_Plus1,"NPU_Plus1/i");     
    tree_->Branch("NPV",&NPV,"NPV/i");           
    tree_->Branch("Rho",&Rho,"Rho/F");           
    tree_->Branch("pt",&pt,"pt/F");            
    tree_->Branch("eta",&eta,"eta/F");           
    tree_->Branch("phi",&phi,"phi/F");           
    tree_->Branch("mass",&mass,"mass/F");    
    tree_->Branch("charge",&charge,"charge/I");
    tree_->Branch("miniIso",&miniIso,"miniIso/F");
    tree_->Branch("met",&met,"met/F");
    tree_->Branch("passTighterTag",&passTighterTag,"passTighterTag/O");
    if (treeType == kTPType_Full) {
      tree_->Branch("Activity",&Activity,"Activity/F");
      tree_->Branch("HT",&HT,"HT/F");
      tree_->Branch("MinDRToParton",&MinDRToParton,"MinDRToParton/F");
      tree_->Branch("NJets",&NJets,"NJets/i");
    }    
  }
  

  // initialze a TagAndProbePair
  void InitTree(int treeType = kTPType_Default){
    assert(tree_);
    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;

    tree_->SetBranchAddress("pass",&pass);
    tree_->SetBranchAddress("weight",&weight);
    tree_->SetBranchAddress("run",&run);
    tree_->SetBranchAddress("lumi",&lumi);
    tree_->SetBranchAddress("event",&event);
    tree_->SetBranchAddress("NPU_0",&NPU_0);
    tree_->SetBranchAddress("NPU_Minus1",&NPU_Minus1);
    tree_->SetBranchAddress("NPU_Plus1",&NPU_Plus1);
    tree_->SetBranchAddress("NPV",&NPV);
    tree_->SetBranchAddress("Rho",&Rho);
    tree_->SetBranchAddress("pt",&pt);
    tree_->SetBranchAddress("eta",&eta);
    tree_->SetBranchAddress("phi",&phi);
    tree_->SetBranchAddress("mass",&mass);
    tree_->SetBranchAddress("charge",&charge);
    tree_->SetBranchAddress("miniIso",&miniIso);
    tree_->SetBranchAddress("met",&met);
    tree_->SetBranchAddress("passTighterTag",&passTighterTag);
    if (treeType == kTPType_Full) {
      tree_->SetBranchAddress("Activity",&Activity);
      tree_->SetBranchAddress("HT",&HT);
      tree_->SetBranchAddress("MinDRToParton",&MinDRToParton);
      tree_->SetBranchAddress("NJets",&NJets);
    }    

    gErrorIgnoreLevel = currentState;
  }

 private:
  
      
}; 


#endif

