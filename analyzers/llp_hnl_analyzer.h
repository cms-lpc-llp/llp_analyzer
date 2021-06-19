#ifndef DEF_llp_hnl_analyzer
#define DEF_llp_hnl_analyzer

#include "RazorAnalyzer.h"

class llp_hnl_analyzer: public RazorAnalyzer {
    public: 
        llp_hnl_analyzer(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
