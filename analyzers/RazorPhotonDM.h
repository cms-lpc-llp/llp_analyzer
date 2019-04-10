#ifndef DEF_RazorPhotonDM
#define DEF_RazorPhotonDM

#include "RazorAnalyzer.h"

class RazorPhotonDM: public RazorAnalyzer {
    public: 
        RazorPhotonDM(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
