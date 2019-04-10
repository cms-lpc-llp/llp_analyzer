#ifndef DEF_HZZRazor
#define DEF_HZZRazor

#include "RazorAnalyzer.h"

class HZZRazor: public RazorAnalyzer {
    public: 
        HZZRazor(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
