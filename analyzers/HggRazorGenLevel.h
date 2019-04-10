#ifndef DEF_HggRazorGenLevel
#define DEF_HggRazorGenLevel

#include "RazorAnalyzer.h"

class HggRazorGenLevel: public RazorAnalyzer {
    public: 
        HggRazorGenLevel(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
