#ifndef DEF_llp_MuonSystem_bparking_rechits
#define DEF_llp_MuonSystem_bparking_rechits

#include "RazorAnalyzer.h"

class llp_MuonSystem_bparking_rechits: public RazorAnalyzer {
    public: 
        llp_MuonSystem_bparking_rechits(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
