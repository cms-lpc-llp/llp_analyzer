#ifndef DEF_FullRazorInclusive
#define DEF_FullRazorInclusive

#include "RazorAnalyzer.h"

class FullRazorInclusive: public RazorAnalyzer {
    public: 
        FullRazorInclusive(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
