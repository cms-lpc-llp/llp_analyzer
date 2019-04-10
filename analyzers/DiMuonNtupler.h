#ifndef DEF_DiMuonNtupler
#define DEF_DiMuonNtupler

#include "RazorAnalyzer.h"

class DiMuonNtupler: public RazorAnalyzer {
    public: 
        DiMuonNtupler(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
