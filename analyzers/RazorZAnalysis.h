#ifndef DEF_RazorZAnalysis
#define DEF_RazorZAnalysis

#include "RazorAnalyzer.h"

class RazorZAnalysis: public RazorAnalyzer {
    public: 
        RazorZAnalysis(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
