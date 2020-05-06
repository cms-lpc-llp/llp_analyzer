#ifndef DEF_SusyLLP
#define DEF_SusyLLP

#include "RazorAnalyzerLLP.h"

class SusyLLP: public RazorAnalyzerLLP {
    public: 
        SusyLLP(TTree *tree=0): RazorAnalyzerLLP(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label, string process);
};

#endif
