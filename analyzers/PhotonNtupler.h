#ifndef DEF_PhotonNtupler
#define DEF_PhotonNtupler

#include "RazorAnalyzer.h"

class PhotonNtupler: public RazorAnalyzer {
    public: 
        PhotonNtupler(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
