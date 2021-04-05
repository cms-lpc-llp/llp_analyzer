#ifndef DEF_SlimSusyLLP
#define DEF_SlimSusyLLP

#include "LLPAnalysis/llpAnalyzer/interface/RazorAnalyzerLLP.h"

class SlimSusyLLP: public RazorAnalyzerLLP {
    public: 
        SlimSusyLLP(TTree *tree=0): RazorAnalyzerLLP(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label, string process);
};

#endif
