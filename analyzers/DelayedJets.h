#ifndef DEF_DelayedJets
#define DEF_DelayedJets

#include "RazorAnalyzer.h"

class DelayedJets: public RazorAnalyzer {
    public: 
        DelayedJets(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
