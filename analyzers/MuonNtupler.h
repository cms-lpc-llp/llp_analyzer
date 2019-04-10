#ifndef DEF_MuonNtupler
#define DEF_MuonNtupler

#include "RazorAnalyzer.h"

class MuonNtupler: public RazorAnalyzer {
    public: 
        MuonNtupler(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
