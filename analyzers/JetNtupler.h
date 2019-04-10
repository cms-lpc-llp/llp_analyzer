#ifndef DEF_JetNtupler
#define DEF_JetNtupler

#include "RazorAnalyzer.h"

class JetNtupler: public RazorAnalyzer {
    public: 
        JetNtupler(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
