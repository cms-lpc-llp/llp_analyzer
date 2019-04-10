#ifndef DEF_MatchedRazorInclusive
#define DEF_MatchedRazorInclusive

#include "RazorAnalyzer.h"

class MatchedRazorInclusive: public RazorAnalyzer {
    public: 
        MatchedRazorInclusive(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
