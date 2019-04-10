#ifndef DEF_RazorMetAna
#define DEF_RazorMetAna

#include "RazorAnalyzer.h"

class RazorMetAna: public RazorAnalyzer {
    public: 
        RazorMetAna(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
