#ifndef DEF_SusyEwkHgg
#define DEF_SusyEwkHgg

#include "RazorAnalyzer.h"

class SusyEwkHgg: public RazorAnalyzer {
    public: 
        SusyEwkHgg(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
