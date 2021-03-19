#ifndef DEF_llp_MuonSystem_hnl
#define DEF_llp_MuonSystem_hnl

#include "RazorAnalyzer.h"

class llp_MuonSystem_hnl: public RazorAnalyzer {
    public: 
        llp_MuonSystem_hnl(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
