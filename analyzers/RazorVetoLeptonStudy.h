#ifndef DEF_RazorVetoLeptonStudy
#define DEF_RazorVetoLeptonStudy

#include "RazorAnalyzer.h"

class RazorVetoLeptonStudy: public RazorAnalyzer {
    public: 
        RazorVetoLeptonStudy(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
