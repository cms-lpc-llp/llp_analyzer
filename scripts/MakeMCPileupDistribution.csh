#!/bin/tcsh

./RazorRun /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p5-Run1/DYToMuMu_powheg.cern.txt MakeMCPileupDistribution  MCPileupDistribution.root -1 DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6
./RazorRun /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p5-Run1/DYToEE_powheg.cern.txt MakeMCPileupDistribution  MCPileupDistribution.root -1 DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6
./RazorRun /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p5-Run1/DYJetsToLL_MG.cern.txt MakeMCPileupDistribution  MCPileupDistribution.root -1 DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball

