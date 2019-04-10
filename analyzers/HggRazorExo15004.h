#ifndef DEF_HggRazorExo15004
#define DEF_HggRazorExo15004

#include "RazorAnalyzer.h"

class HggRazorExo15004: public RazorAnalyzer {
    public: 
        HggRazorExo15004(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
