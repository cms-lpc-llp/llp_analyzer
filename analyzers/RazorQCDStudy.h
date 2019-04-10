#ifndef DEF_RazorQCDStudy
#define DEF_RazorQCDStudy

#include "RazorAnalyzer.h"

class RazorQCDStudy: public RazorAnalyzer {
    public: 
        RazorQCDStudy(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
