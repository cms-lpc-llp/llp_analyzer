#ifndef DEF_EventPick
#define DEF_EventPick

#include "RazorAnalyzer.h"

class EventPick: public RazorAnalyzer {
    public: 
        EventPick(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
