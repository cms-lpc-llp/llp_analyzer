#ifndef DEF_Muon_Hits
#define DEF_Muon_Hits

#include "RazorAnalyzer.h"

class Muon_Hits: public RazorAnalyzer {
    public: 
        Muon_Hits(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
