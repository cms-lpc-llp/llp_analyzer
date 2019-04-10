#ifndef DEF_RazorDM
#define DEF_RazorDM

#include "RazorAnalyzer.h"

class RazorDM: public RazorAnalyzer {
    public: 
        RazorDM(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
