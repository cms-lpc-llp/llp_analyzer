#ifndef DEF_llp_vH_MuonSystem
#define DEF_llp_vH_MuonSystem

#include "RazorAnalyzer.h"

class llp_vH_MuonSystem: public RazorAnalyzer {
    public: 
        llp_vH_MuonSystem(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
