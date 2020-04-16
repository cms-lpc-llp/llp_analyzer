#ifndef DEF_MakeHiggsPtDistribution
#define DEF_MakeHiggsPtDistribution

#include "RazorAnalyzer.h"

class MakeHiggsPtDistribution: public RazorAnalyzer {
    public: 
        MakeHiggsPtDistribution(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
