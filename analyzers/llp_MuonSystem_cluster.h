#ifndef DEF_llp_MuonSystem_cluster
#define DEF_llp_MuonSystem_cluster

#include "RazorAnalyzer.h"

class llp_MuonSystem_cluster: public RazorAnalyzer {
    public: 
        llp_MuonSystem_cluster(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
