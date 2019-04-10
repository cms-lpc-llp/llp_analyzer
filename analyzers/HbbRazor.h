#ifndef DEF_HbbRazor
#define DEF_HbbRazor

#include "RazorAnalyzer.h"

class HbbRazor: public RazorAnalyzer {
    public: 
        HbbRazor(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
