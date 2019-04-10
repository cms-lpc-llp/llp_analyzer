#ifndef DEF_HggRazorLeptons
#define DEF_HggRazorLeptons

#include "RazorAnalyzer.h"

class HggRazorLeptons: public RazorAnalyzer {
    public: 
        HggRazorLeptons(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
