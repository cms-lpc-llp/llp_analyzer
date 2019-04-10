#ifndef DEF_WWZAnalysis
#define DEF_WWZAnalysis

#include "RazorAnalyzer.h"

class WWZAnalysis: public RazorAnalyzer {
    public: 
        WWZAnalysis(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
