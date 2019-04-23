#ifndef DEF_llp_vH
#define DEF_llp_vH

#include "RazorAnalyzer.h"

class llp_vH: public RazorAnalyzer {
    public:
        llp_vH(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
