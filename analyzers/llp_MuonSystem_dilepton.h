#ifndef DEF_llp_MuonSystem_dilepton
#define DEF_llp_MuonSystem_dilepton

#include "RazorAnalyzer.h"

class llp_MuonSystem_dilepton: public RazorAnalyzer {
    public:
        llp_MuonSystem_dilepton(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData,  int option, string outputFileName, string label);
};

#endif
