#ifndef DEF_HggRazorUpgradeTiming
#define DEF_HggRazorUpgradeTiming

#include "RazorAnalyzerUpgradeTiming.h"

class HggRazorUpgradeTiming: public RazorAnalyzerUpgradeTiming {
    public: 
        HggRazorUpgradeTiming(TTree *tree=0): RazorAnalyzerUpgradeTiming(tree) { }
        void Analyze(bool isData, bool useTiming, bool usePhoChi2, bool useOddEvent, int option, string outputFileName, string label);
};

#endif
