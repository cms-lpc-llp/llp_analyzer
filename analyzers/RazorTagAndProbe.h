#ifndef DEF_RazorTagAndProbe
#define DEF_RazorTagAndProbe

#include "RazorAnalyzer.h"

class RazorTagAndProbe: public RazorAnalyzer {
    public: 
        RazorTagAndProbe(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
