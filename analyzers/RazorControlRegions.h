#ifndef DEF_RazorControlRegions
#define DEF_RazorControlRegions

#include "RazorAnalyzer.h"

class RazorControlRegions: public RazorAnalyzer {
    public: 
        RazorControlRegions(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
