#ifndef DEF_SlimSusyLLPPF
#define DEF_SlimSusyLLPPF

#include "LLPAnalysis/llpAnalyzer/interface/RazorAnalyzerLLP.h"

class SlimSusyLLPPF: public RazorAnalyzerLLP {
    public: 
        SlimSusyLLPPF(TTree *tree=0): RazorAnalyzerLLP(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label, string process);
};

#endif
