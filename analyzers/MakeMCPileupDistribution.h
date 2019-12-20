#ifndef DEF_MakeMCPileupDistribution
#define DEF_MakeMCPileupDistribution

#include "RazorAnalyzer.h"

class MakeMCPileupDistribution: public RazorAnalyzer {
    public: 
        MakeMCPileupDistribution(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
