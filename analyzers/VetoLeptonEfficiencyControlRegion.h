#ifndef DEF_VetoLeptonEfficiencyControlRegion
#define DEF_VetoLeptonEfficiencyControlRegion

#include "RazorAnalyzer.h"

class VetoLeptonEfficiencyControlRegion: public RazorAnalyzer {
    public: 
        VetoLeptonEfficiencyControlRegion(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
