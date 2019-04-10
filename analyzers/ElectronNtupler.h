#ifndef DEF_ElectronNtupler
#define DEF_ElectronNtupler

#include "RazorAnalyzer.h"

class ElectronNtupler: public RazorAnalyzer {
    public: 
        ElectronNtupler(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
