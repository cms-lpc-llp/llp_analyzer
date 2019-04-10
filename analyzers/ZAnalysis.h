#ifndef DEF_ZAnalysis
#define DEF_ZAnalysis

#include "RazorAnalyzer.h"

class ZAnalysis: public RazorAnalyzer {
    public: 
        ZAnalysis(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
