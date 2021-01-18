#ifndef DEF_SusyLLPPassFail
#define DEF_SusyLLPPassFail

#include "RazorAnalyzerLLP.h"

class SusyLLPPassFail: public RazorAnalyzerLLP {
    public: 
        SusyLLPPassFail(TTree *tree=0): RazorAnalyzerLLP(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label, string process);
};

#endif
