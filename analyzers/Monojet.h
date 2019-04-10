#ifndef DEF_Monojet
#define DEF_Monojet

#include "RazorAnalyzer.h"

class Monojet: public RazorAnalyzer {
    public: 
        Monojet(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
