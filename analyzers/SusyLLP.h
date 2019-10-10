#ifndef DEF_SusyLLP
#define DEF_SusyLLP

#include "RazorAnalyzer.h"

class SusyLLP: public RazorAnalyzer {
    public: 
        SusyLLP(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
