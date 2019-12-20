#ifndef DEF_RazorLiteZ
#define DEF_RazorLiteZ

#include "RazorAnalyzer.h"

class RazorLiteZ: public RazorAnalyzer {
    public: 
        RazorLiteZ(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
