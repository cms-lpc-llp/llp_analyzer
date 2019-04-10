#ifndef DEF_RazorPhotonStudy
#define DEF_RazorPhotonStudy

#include "RazorAnalyzer.h"

class RazorPhotonStudy: public RazorAnalyzer {
    public: 
        RazorPhotonStudy(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
