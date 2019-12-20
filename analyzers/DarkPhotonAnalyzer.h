#ifndef DEF_DarkPhotonAnalyzer
#define DEF_DarkPhotonAnalyzer

#include "RazorAnalyzer.h"

class DarkPhotonAnalyzer: public RazorAnalyzer {
    public: 
        DarkPhotonAnalyzer(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
