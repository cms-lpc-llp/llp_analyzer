#ifndef DEF_HggRazorRun1WithRun2Bin
#define DEF_HggRazorRun1WithRun2Bin

#include "RazorAnalyzerRun1.h"

class HggRazorRun1WithRun2Bin: public RazorAnalyzerRun1 {
    public: 
        HggRazorRun1WithRun2Bin(TTree *tree=0): RazorAnalyzerRun1(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
