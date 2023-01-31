#ifndef DEF_llp_MuonSystem_bparking_with_SF
#define DEF_llp_MuonSystem_bparking_with_SF

#include "RazorAnalyzer.h"

class llp_MuonSystem_bparking_with_SF: public RazorAnalyzer {
    public: 
        llp_MuonSystem_bparking_with_SF(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
