#ifndef DEF_llp_MuonSystem_TnP
#define DEF_llp_MuonSystem_TnP

#include "RazorAnalyzer.h"

class llp_MuonSystem_TnP: public RazorAnalyzer {
    public: 
        llp_MuonSystem_TnP(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
