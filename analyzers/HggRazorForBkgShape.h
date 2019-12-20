#ifndef DEF_HggRazorForBkgShape
#define DEF_HggRazorForBkgShape

#include "RazorAnalyzer.h"

class HggRazorForBkgShape: public RazorAnalyzer {
    public: 
        HggRazorForBkgShape(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
