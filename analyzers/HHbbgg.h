#ifndef DEF_HHbbgg
#define DEF_HHbbgg

#include "RazorAnalyzer.h"

class HHbbgg: public RazorAnalyzer {
    public: 
        HHbbgg(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
