#ifndef DEF_llp_MuonSystem_darkshower
#define DEF_llp_MuonSystem_darkshower

#include "RazorAnalyzer.h"

class llp_MuonSystem_darkshower: public RazorAnalyzer {
    public: 
        llp_MuonSystem_darkshower(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
