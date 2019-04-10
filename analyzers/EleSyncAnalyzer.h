#ifndef DEF_EleSyncAnalyzer
#define DEF_EleSyncAnalyzer

#include "RazorAnalyzer.h"

class EleSyncAnalyzer: public RazorAnalyzer {
    public: 
        EleSyncAnalyzer(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
