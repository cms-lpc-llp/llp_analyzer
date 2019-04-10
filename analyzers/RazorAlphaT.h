#ifndef DEF_RazorAlphaT
#define DEF_RazorAlphaT

#include "RazorAnalyzer.h"

class RazorAlphaT: public RazorAnalyzer {
    public: 
        RazorAlphaT(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
