#!/bin/bash

#python python/DustinTuple2RooDataSet.py -c config/run2_2016.config -d Datasets -b MultiJet --data -l 4400 Run2016G_Unblinded/FullRazorInclusive_HTMHT_2016G_GoodLumiGolden_SUSYUnblind.root

#python python/DustinTuple2RooDataSet.py -c config/run2_2016.config -d Datasets -b DiJet --data -l 4400 Run2016G_Unblinded/FullRazorInclusive_HTMHT_2016G_GoodLumiGolden_SUSYUnblind.root

#python python/DustinTuple2RooDataSet.py -c config/run2_2016.config -d Datasets -b LeptonMultiJet --data -l 4400 Run2016G_Unblinded/FullRazorInclusive_SingleLepton_2016G_GoodLumiGolden_SUSYUnblind.root

#python python/DustinTuple2RooDataSet.py -c config/run2_2016.config -d Datasets -b LeptonJet --data -l 4400 Run2016G_Unblinded/FullRazorInclusive_SingleLepton_2016G_GoodLumiGolden_SUSYUnblind.root

#MultiJet Full

#python python/BinnedFit.py -b MultiJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/MultiJet/Full/ Datasets/FullRazorInclusive_HTMHT_2016G_GoodLumiGolden_SUSYUnblind_lumi-4.400_0-3btag_MultiJet.root --data -l 4400 

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/MultiJet/Full/ -b MultiJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/MultiJet/Full/BinnedFitResults_MultiJet.root --no-stat -t 1000

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/MultiJet/Full/ -b MultiJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/MultiJet/Full/BinnedFitResults_MultiJet.root  -t 1000

#python python/PlotFit.py -b MultiJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/MultiJet/Full/ -i fits_2016_11_14/Run2016G_Unblinded/MultiJet/Full/BinnedFitResults_MultiJet.root --data -l 4400 -t fits_2016_11_14/Run2016G_Unblinded/MultiJet/Full/toys_Bayes_MultiJet.root -s fits_2016_11_14/Run2016G_Unblinded/MultiJet/Full/toys_Bayes_noStat_MultiJet.root --no-stat --print-errors


#MultiJet Sideband

#python python/BinnedFit.py -b MultiJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/MultiJet/Sideband/ Datasets/FullRazorInclusive_HTMHT_2016G_GoodLumiGolden_SUSYUnblind_lumi-4.400_0-3btag_MultiJet.root --data -l 4400 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/MultiJet/Sideband/ -b MultiJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/MultiJet/Sideband/BinnedFitResults_MultiJet.root --no-stat -t 1000 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/MultiJet/Sideband/ -b MultiJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/MultiJet/Sideband/BinnedFitResults_MultiJet.root  -t 1000 --fit-region LowMR,LowRsq

#python python/PlotFit.py -b MultiJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/MultiJet/Sideband/ -i fits_2016_11_14/Run2016G_Unblinded/MultiJet/Sideband/BinnedFitResults_MultiJet.root --data -l 4400 -t fits_2016_11_14/Run2016G_Unblinded/MultiJet/Sideband/toys_Bayes_MultiJet.root -s fits_2016_11_14/Run2016G_Unblinded/MultiJet/Sideband/toys_Bayes_noStat_MultiJet.root --fit-region LowMR,LowRsq --no-stat --print-errors

#DiJet Sideband

#python python/BinnedFit.py -b DiJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/DiJet/Sideband/ Datasets/FullRazorInclusive_HTMHT_2016G_GoodLumiGolden_SUSYUnblind_lumi-4.400_0-2btag_DiJet.root --data -l 4400 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/DiJet/Sideband/ -b DiJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/DiJet/Sideband/BinnedFitResults_DiJet.root --no-stat -t 1000 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/DiJet/Sideband/ -b DiJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/DiJet/Sideband/BinnedFitResults_DiJet.root  -t 1000 --fit-region LowMR,LowRsq

#python python/PlotFit.py -b DiJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/DiJet/Sideband/ -i fits_2016_11_14/Run2016G_Unblinded/DiJet/Sideband/BinnedFitResults_DiJet.root --data -l 4400 -t fits_2016_11_14/Run2016G_Unblinded/DiJet/Sideband/toys_Bayes_DiJet.root -s fits_2016_11_14/Run2016G_Unblinded/DiJet/Sideband/toys_Bayes_noStat_DiJet.root --fit-region LowMR,LowRsq --no-stat --print-errors

#DiJet Full

#python python/BinnedFit.py -b DiJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/DiJet/Full/ Datasets/FullRazorInclusive_HTMHT_2016G_GoodLumiGolden_SUSYUnblind_lumi-4.400_0-2btag_DiJet.root --data -l 4400 

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/DiJet/Full/ -b DiJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/DiJet/Full/BinnedFitResults_DiJet.root --no-stat -t 1000

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/DiJet/Full/ -b DiJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/DiJet/Full/BinnedFitResults_DiJet.root  -t 1000

#python python/PlotFit.py -b DiJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/DiJet/Full/ -i fits_2016_11_14/Run2016G_Unblinded/DiJet/Full/BinnedFitResults_DiJet.root --data -l 4400 -t fits_2016_11_14/Run2016G_Unblinded/DiJet/Full/toys_Bayes_DiJet.root -s fits_2016_11_14/Run2016G_Unblinded/DiJet/Full/toys_Bayes_noStat_DiJet.root --no-stat --print-errors

#LeptonMultiJet Sideband

#python python/BinnedFit.py -b LeptonMultiJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Sideband/ Datasets/FullRazorInclusive_SingleLepton_2016G_GoodLumiGolden_SUSYUnblind_lumi-4.400_0-3btag_LeptonMultiJet.root --data -l 4400 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Sideband/ -b LeptonMultiJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Sideband/BinnedFitResults_LeptonMultiJet.root --no-stat -t 1000 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Sideband/ -b LeptonMultiJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Sideband/BinnedFitResults_LeptonMultiJet.root  -t 1000 --fit-region LowMR,LowRsq

#python python/PlotFit.py -b LeptonMultiJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Sideband/ -i fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Sideband/BinnedFitResults_LeptonMultiJet.root --data -l 4400 -t fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Sideband/toys_Bayes_LeptonMultiJet.root -s fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Sideband/toys_Bayes_noStat_LeptonMultiJet.root --fit-region LowMR,LowRsq --no-stat --print-errors

#LeptonMultiJet Full

#python python/BinnedFit.py -b LeptonMultiJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Full/ Datasets/FullRazorInclusive_SingleLepton_2016G_GoodLumiGolden_SUSYUnblind_lumi-4.400_0-3btag_LeptonMultiJet.root --data -l 4400 

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Full/ -b LeptonMultiJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Full/BinnedFitResults_LeptonMultiJet.root --no-stat -t 1000

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Full/ -b LeptonMultiJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Full/BinnedFitResults_LeptonMultiJet.root  -t 1000

#python python/PlotFit.py -b LeptonMultiJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Full/ -i fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Full/BinnedFitResults_LeptonMultiJet.root --data -l 4400 -t fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Full/toys_Bayes_LeptonMultiJet.root -s fits_2016_11_14/Run2016G_Unblinded/LeptonMultiJet/Full/toys_Bayes_noStat_LeptonMultiJet.root --no-stat --print-errors

#LeptonJet Sideband

#python python/BinnedFit.py -b LeptonJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Sideband/ Datasets/FullRazorInclusive_SingleLepton_2016G_GoodLumiGolden_SUSYUnblind_lumi-4.400_0-3btag_LeptonJet.root --data -l 4400 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Sideband/ -b LeptonJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Sideband/BinnedFitResults_LeptonJet.root --no-stat -t 1000 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Sideband/ -b LeptonJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Sideband/BinnedFitResults_LeptonJet.root  -t 1000 --fit-region LowMR,LowRsq

#python python/PlotFit.py -b LeptonJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Sideband/ -i fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Sideband/BinnedFitResults_LeptonJet.root --data -l 4400 -t fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Sideband/toys_Bayes_LeptonJet.root -s fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Sideband/toys_Bayes_noStat_LeptonJet.root --fit-region LowMR,LowRsq --no-stat --print-errors

#LeptonJet Full

#python python/BinnedFit.py -b LeptonJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Full/ Datasets/FullRazorInclusive_SingleLepton_2016G_GoodLumiGolden_SUSYUnblind_lumi-4.400_0-2btag_LeptonJet.root --data -l 4400 

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Full/ -b LeptonJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Full/BinnedFitResults_LeptonJet.root --no-stat -t 1000

#python python/RunToys.py -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Full/ -b LeptonJet -l 4400 -i fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Full/BinnedFitResults_LeptonJet.root  -t 1000

#python python/PlotFit.py -b LeptonJet -c config/run2_2016.config -d fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Full/ -i fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Full/BinnedFitResults_LeptonJet.root --data -l 4400 -t fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Full/toys_Bayes_LeptonJet.root -s fits_2016_11_14/Run2016G_Unblinded/LeptonJet/Full/toys_Bayes_noStat_LeptonJet.root --no-stat --print-errors
