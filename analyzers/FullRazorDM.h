#ifndef DEF_FullRazorDM
#define DEF_FullRazorDM

#include "RazorAnalyzer.h"

class FullRazorDM: public RazorAnalyzer {
    public: 
        FullRazorDM(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
