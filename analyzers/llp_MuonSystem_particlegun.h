#ifndef DEF_llp_MuonSystem_particlegun
#define DEF_llp_MuonSystem_particlegun

#include "RazorAnalyzer.h"

class llp_MuonSystem_particlegun: public RazorAnalyzer {
    public: 
        llp_MuonSystem_particlegun(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
