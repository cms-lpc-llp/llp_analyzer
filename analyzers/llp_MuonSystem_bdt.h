#ifndef DEF_llp_MuonSystem_bdt
#define DEF_llp_MuonSystem_bdt

#include "RazorAnalyzer.h"

class llp_MuonSystem_bdt: public RazorAnalyzer {
    public:
        llp_MuonSystem_bdt(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
