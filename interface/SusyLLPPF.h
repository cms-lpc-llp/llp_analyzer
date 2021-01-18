#ifndef DEF_SusyLLPPF
#define DEF_SusyLLPPF

#include "LLPAnalysis/llpAnalyzer/interface/RazorAnalyzerLLP.h"

class SusyLLPPF: public RazorAnalyzerLLP {
    public: 
        SusyLLPPF(TTree *tree=0): RazorAnalyzerLLP(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label, string process);
};

#endif
