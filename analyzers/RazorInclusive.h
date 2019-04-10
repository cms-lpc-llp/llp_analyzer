#ifndef DEF_RazorInclusive
#define DEF_RazorInclusive

#include "RazorAnalyzer.h"

class RazorInclusive: public RazorAnalyzer {
    public: 
        RazorInclusive(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
