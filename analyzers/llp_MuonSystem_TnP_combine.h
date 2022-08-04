#ifndef DEF_llp_MuonSystem_TnP_combine
#define DEF_llp_MuonSystem_TnP_combine

#include "RazorAnalyzer.h"

class llp_MuonSystem_TnP_combine: public RazorAnalyzer {
    public: 
        llp_MuonSystem_TnP_combine(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
