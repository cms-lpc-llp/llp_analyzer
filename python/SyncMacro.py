##### Created 26 September 2015
##### Macro to reproduce control region plots using Run2015D data

import sys
import copy
import argparse
import ROOT as rt
import numpy as np

#local imports
from macro import macro
from macro.razorAnalysis import *

LUMI = 157 #in /pb
MCLUMI = 1 

SAMPLES_DYJ2L = ["SingleTop", "WJets", "VV", "TTJets", "DYJets"]
SAMPLES_TTJ2L = ["SingleTop", "WJets", "VV", "DYJets", "TTJets"]

weightfilenames = {
        "muon": "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/20150924_PR_2015D/efficiency_results_TightMuonSelectionEffDenominatorReco_2015D.root",
        "ele": "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/LeptonEfficiencies/20150924_PR_2015D/efficiency_results_TightElectronSelectionEffDenominatorReco_2015D.root",
        "pileup": "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PileupWeights/NVtxReweight_ZToMuMu_2015D_25ns_20150923.root",
        }
weighthistnames = {
        "muon": "ScaleFactor_TightMuonSelectionEffDenominatorReco",
        "ele": "ScaleFactor_TightElectronSelectionEffDenominatorReco",
        "pileup": "NVtxReweight",
        }

dir_1L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p18_OLD/"
FILENAMES_1L = {
        "DYJets"   :dir_1L+"RunTwoRazorControlRegions_DileptonFull_DileptonSkim_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root",
        "TTJets"   :dir_1L+"RunTwoRazorControlRegions_DileptonFull_DileptonSkim_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root",
        "WJets"    :dir_1L+"RunTwoRazorControlRegions_DileptonFull_DileptonSkim_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root",
        "VV"       :dir_1L+"RunTwoRazorControlRegions_DileptonFull_DileptonSkim_VV_1pb_weighted.root",
        "SingleTop":dir_1L+"RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleTop_1pb_weighted.root",
        "DataEMu"  :dir_1L+"RunTwoRazorControlRegions_DileptonFull_DileptonSkim_SingleLepton_Run2015D_GoodLumi_NoDuplicates.root",
            }

dyjetsDielectronCuts = "abs(lep1Type) == 11 && abs(lep2Type) == 11 && lep1PassTight && lep2PassTight && lep1.Pt() > 30 && lep2.Pt() > 20 && mll > 60 && mll < 120"
dyjetsDielectronCutsData = appendTriggerCuts(dyjetsDielectronCuts, singleLeptonTriggerNumsData)
dyjetsDielectronCutsMC = appendTriggerCuts(dyjetsDielectronCuts, singleLeptonTriggerNumsMC)

dyjetsDimuonCuts = "abs(lep1Type) == 13 && abs(lep2Type) == 13 && lep1PassTight && lep2PassTight && lep1.Pt() > 30 && lep2.Pt() > 20 && mll > 60 && mll < 120"
dyjetsDimuonCutsData = appendTriggerCuts(dyjetsDimuonCuts, singleLeptonTriggerNumsData)
dyjetsDimuonCutsMC = appendTriggerCuts(dyjetsDimuonCuts, singleLeptonTriggerNumsMC)

dyjetsDileptonCuts = "((abs(lep1Type) == 11 && abs(lep2Type) == 11) || (abs(lep1Type) == 13 && abs(lep2Type) == 13)) && lep1PassTight && lep2PassTight && lep1.Pt() > 30 && lep2.Pt() > 20 && mll > 60 && mll < 120"
dyjetsDileptonCutsData = appendTriggerCuts(dyjetsDileptonCuts, singleLeptonTriggerNumsData)
dyjetsDileptonCutsMC = appendTriggerCuts(dyjetsDileptonCuts, singleLeptonTriggerNumsMC)

dyjetsDileptonFullCutsMC = dyjetsDileptonCutsMC+" && NJets80 >= 2 && NBJetsMedium == 0 && MR > 300"
dyjetsDileptonFullCutsData = dyjetsDileptonCutsData+" && NJets80 >= 2 && NBJetsMedium == 0 && MR > 300"

dyjetsDileptonBins = {
        "MR"     : np.arange(200, 2700, 50),
        "Rsq"    : np.arange(0., 1., 0.05),
        "mll"    : np.arange(60, 120, 1.0),
        "MET"    : np.arange(0, 200, 4.0), 
        "NJets40": np.arange(-0.5, 9.5, 1),
        "NJets80": np.arange(-0.5, 9.5, 1),
        }

ttjetsEMuCuts = "((abs(lep1Type) == 11 && abs(lep2Type) == 13) || (abs(lep1Type) == 13 && abs(lep2Type) == 11)) && lep1PassTight && lep2PassTight && lep1.Pt() > 30 && lep2.Pt() > 20 && mll > 20"
ttjetsEMuCutsMC = appendTriggerCuts(ttjetsEMuCuts, singleLeptonTriggerNumsMC)
ttjetsEMuCutsData = appendTriggerCuts(ttjetsEMuCuts, singleLeptonTriggerNumsData)

ttjetsEMuFullCuts = "((abs(lep1Type) == 11 && abs(lep2Type) == 13) || (abs(lep1Type) == 13 && abs(lep2Type) == 11)) && lep1PassTight && lep2PassTight && lep1.Pt() > 30 && lep2.Pt() > 20 && mll > 50 && NBJetsMedium > 0"
ttjetsEMuFullCutsMC = appendTriggerCuts(ttjetsEMuFullCuts, singleLeptonTriggerNumsMC)
ttjetsEMuFullCutsData = appendTriggerCuts(ttjetsEMuFullCuts, singleLeptonTriggerNumsData)

ttjetsDileptonBins = {
        "MR"  : np.arange(0, 1500, 100),
        "Rsq" : np.arange(0, 1.5, 0.01),
        "MET" : np.arange(0, 500, 20),
        "mll" : np.arange(0, 500, 12.5)
        }

weightOpts = ["doPileupWeights", "doLep1Weights", "doLep2Weights"]

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug

    #initialize
    weightHists = loadWeightHists(filenames=weightfilenames, histnames=weighthistnames, debugLevel=debugLevel)
    sfHists = {}

    #DYJets control sample
    #dyjetsDielectronHists = makeControlSampleHists("DYJetsDielectron", filenames=FILENAMES_1L, samples=SAMPLES_DYJ2L, 
    #            cutsMC=dyjetsDielectronCutsMC, cutsData=dyjetsDielectronCutsData, 
    #            bins=dyjetsDileptonBins, logX=False, lumiMC=MCLUMI, lumiData=LUMI, 
    #            weightHists=weightHists, sfHists=sfHists, dataName="DataEMu", weightOpts=weightOpts, debugLevel=debugLevel)

    #dyjetsDimuonHists = makeControlSampleHists("DYJetsDimuon", filenames=FILENAMES_1L, samples=SAMPLES_DYJ2L, 
    #            cutsMC=dyjetsDimuonCutsMC, cutsData=dyjetsDimuonCutsData, 
    #            bins=dyjetsDileptonBins, logX=False, lumiMC=MCLUMI, lumiData=LUMI, 
    #            weightHists=weightHists, sfHists=sfHists, dataName="DataEMu", weightOpts=weightOpts, debugLevel=debugLevel)

    #dyjetsDileptonHists = makeControlSampleHists("DYJetsDilepton", filenames=FILENAMES_1L, samples=SAMPLES_DYJ2L, 
    #            cutsMC=dyjetsDileptonCutsMC, cutsData=dyjetsDileptonCutsData, 
    #            bins=dyjetsDileptonBins, logX=False, lumiMC=MCLUMI, lumiData=LUMI, 
    #            weightHists=weightHists, sfHists=sfHists, dataName="DataEMu", weightOpts=weightOpts, debugLevel=debugLevel)
    #
    #dyjetsDileptonFullHists = makeControlSampleHists("DYJetsDileptonFull", filenames=FILENAMES_1L, samples=SAMPLES_DYJ2L, 
    #            cutsMC=dyjetsDileptonFullCutsMC, cutsData=dyjetsDileptonFullCutsData, 
    #            bins=dyjetsDileptonBins, logX=False, lumiMC=MCLUMI, lumiData=LUMI, 
    #            weightHists=weightHists, sfHists=sfHists, dataName="DataEMu", weightOpts=weightOpts, debugLevel=debugLevel)

    #ttjetsEMuHists = makeControlSampleHists("TTJetsEMu", filenames=FILENAMES_1L, samples=SAMPLES_TTJ2L, 
    #            cutsMC=ttjetsEMuCutsMC, cutsData=ttjetsEMuCutsData, 
    #            bins=ttjetsDileptonBins, logX=False, lumiMC=MCLUMI, lumiData=LUMI, 
    #            weightHists=weightHists, sfHists=sfHists, dataName="DataEMu", weightOpts=weightOpts, debugLevel=debugLevel)

    #ttjetsEMuFullHists = makeControlSampleHists("TTJetsEMuFull", filenames=FILENAMES_1L, samples=SAMPLES_TTJ2L, 
    #            cutsMC=ttjetsEMuFullCutsMC, cutsData=ttjetsEMuFullCutsData, 
    #            bins=ttjetsDileptonBins, logX=False, lumiMC=MCLUMI, lumiData=LUMI, 
    #            weightHists=weightHists, sfHists=sfHists, dataName="DataEMu", weightOpts=weightOpts, debugLevel=debugLevel)
