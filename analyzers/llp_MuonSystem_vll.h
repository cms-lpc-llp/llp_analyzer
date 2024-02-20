#ifndef DEF_llp_MuonSystem_vll
#define DEF_llp_MuonSystem_vll

#include "RazorAnalyzer.h"

class llp_MuonSystem_vll: public RazorAnalyzer {
    public: 
        llp_MuonSystem_vll(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
