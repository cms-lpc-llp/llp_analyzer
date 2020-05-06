#!/bin/sh
# This is the order to run the scripts in this folder
# in order to create the QCD scale factor files
# needed to make the signal region prediction.
root -l calculateMCScaleFactors.C++
python getInclusiveTransferFactors.py
python QCDFits2D.py
python EvaluateFits.py
python getBTagFractions.py
cp qcd_data_btags.root ../../../data/ScaleFactors/RazorMADD2015/RazorQCDBTagScaleFactors_Razor2016_MoriondRereco.root
