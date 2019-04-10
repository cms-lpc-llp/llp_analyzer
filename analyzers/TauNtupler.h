#ifndef DEF_TauNtupler
#define DEF_TauNtupler

#include "RazorAnalyzer.h"

class TauNtupler: public RazorAnalyzer {
    public: 
        TauNtupler(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
