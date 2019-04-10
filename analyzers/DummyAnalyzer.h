#ifndef DEF_DummyAnalyzer
#define DEF_DummyAnalyzer

#include "RazorAnalyzer.h"

class DummyAnalyzer: public RazorAnalyzer {
    public: 
        DummyAnalyzer(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
