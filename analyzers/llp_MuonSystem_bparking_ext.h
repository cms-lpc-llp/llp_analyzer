#ifndef DEF_llp_MuonSystem_bparking_ext
#define DEF_llp_MuonSystem_bparking_ext

#include "RazorAnalyzer.h"

class llp_MuonSystem_bparking_ext: public RazorAnalyzer {
    public: 
        llp_MuonSystem_bparking_ext(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif