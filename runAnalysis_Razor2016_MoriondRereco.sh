#!/bin/sh
set -e
python python/ComputeScaleFactorsMacro.py 
python python/ComputeScaleFactorsNJetCorrectionByProcess.py 
python python/OneLeptonClosureTestMacro.py
python python/BTagClosureTestMacro.py
python python/TTJetsDileptonCrossCheck.py 
python python/GJetsBTagClosureTestMacro.py
python python/DYJetsInvCrossCheck.py 
root -l 'macros/BackgroundStudies/plotScaleFactorHistograms.C("Razor2016_MoriondRereco")'
python python/DYJetsInvCrossCheck.py --closure
python python/DYJetsInvCrossCheck_Nb.py 
python python/VetoLeptonCrossCheck.py 
python python/SignalRegionMacro.py --unblind --fine-grained
for nb in 0 1 2 3; do
    python python/SignalRegionPlotMacro.py --box MultiJet --btags $nb --fine-grained --unblind
    python python/SignalRegionPlotMacro.py --box LeptonMultiJet --btags $nb --fine-grained --unblind
done
for nb in 0 1 2; do
    python python/SignalRegionPlotMacro.py --box DiJet --btags $nb --fine-grained --unblind
    python python/SignalRegionPlotMacro.py --box LeptonJet --btags $nb --fine-grained  --unblind
done
