#ifndef DEF_llp_MuonSystem_bparking
#define DEF_llp_MuonSystem_bparking

#include "RazorAnalyzer.h"

class llp_MuonSystem_bparking: public RazorAnalyzer {
    public: 
        llp_MuonSystem_bparking(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
