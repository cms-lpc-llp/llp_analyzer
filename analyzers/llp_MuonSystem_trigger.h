#ifndef DEF_llp_MuonSystem_trigger
#define DEF_llp_MuonSystem_trigger

#include "RazorAnalyzer.h"

class llp_MuonSystem_trigger: public RazorAnalyzer {
    public: 
        llp_MuonSystem_trigger(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
