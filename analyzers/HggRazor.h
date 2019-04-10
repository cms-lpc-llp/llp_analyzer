#ifndef DEF_HggRazor
#define DEF_HggRazor

#include "RazorAnalyzer.h"

class HggRazor: public RazorAnalyzer {
    public: 
        HggRazor(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
