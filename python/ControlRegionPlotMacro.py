import sys,os
import argparse
import copy
import ROOT as rt

#local imports
from framework import Config
from macro.razorAnalysis import xbinsSignal, colsSignal
from macro.razorMacros import plotControlSampleHists
import macro.macro as macro
from SidebandMacro import LUMI, MCLUMI
from ComputeScaleFactorsMacro import SAMPLES_TTJ1L, SAMPLES_WJ1L, SAMPLES_WJ1L_INV, ControlRegionBinning
from CrossCheckRegionMacro import SAMPLES_VetoLepton, SAMPLES_VetoTau, SAMPLES_TTJ2L

SAMPLES_GJ = ["Other", "QCD", "GJetsFrag", "GJets"]
SAMPLES = { "WJetControlRegion":SAMPLES_WJ1L, "WJetInvControlRegion":SAMPLES_WJ1L_INV, 
    "TTJetsSingleLeptonControlRegion":SAMPLES_TTJ1L, "VetoLeptonControlRegion":SAMPLES_VetoLepton, 
    "VetoTauControlRegion":SAMPLES_VetoTau, "GJetsInvControlRegion":SAMPLES_GJ,
    "TTJetsDileptonControlRegion":SAMPLES_TTJ2L }

SAMPLES_TTJ_REDUCED = ["Other", "WJets", "TTJets"]
SAMPLES_WJ_REDUCED = ["Other", "QCD", "TTJets", "WJets"]
SAMPLES_DY_REDUCED = ["Other", "WJets", "TTJets", "DYJets"]
SAMPLES_VETO_REDUCED = ["Other", "ZInv", "QCD", "WJets", "TTJets"]
SAMPLES_WJI_REDUCED = ["Other", "ZInv", "QCD", "TTJets", "WJetsInv"]
SAMPLES_REDUCED = { "WJetControlRegion":SAMPLES_WJ_REDUCED, "WJetInvControlRegion":SAMPLES_WJI_REDUCED, 
    "TTJetsSingleLeptonControlRegion":SAMPLES_TTJ_REDUCED, "VetoLeptonControlRegion":SAMPLES_VETO_REDUCED, 
    "VetoTauControlRegion":SAMPLES_VETO_REDUCED, "GJetsInvControlRegion":SAMPLES_GJ,
    "TTJetsDileptonControlRegion":SAMPLES_TTJ_REDUCED }

boxMapping = { "WJetControlRegion":"WJetsSingleLepton", "WJetInvControlRegion":"WJetsSingleLeptonInv", 
    "TTJetsSingleLeptonControlRegion":"TTJetsSingleLepton", "VetoLeptonControlRegion":"VetoLeptonControlRegion", 
    "VetoTauControlRegion":"VetoTauControlRegion", "GJetsInvControlRegion":"GJetsInvControlRegion",
    "TTJetsDileptonControlRegion":"TTJetsDileptonControlRegion" }

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument('--box', help="choose a box")
    parser.add_argument('--dir', help="output directory (should contain the ROOT files with the razor histograms)", default="ControlSamplePlots", dest='dirName')
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    dirName = args.dirName

    plotOpts = { "ymin":1e-3, 'comment':False, 'SUS15004CR':True}
    plotOpts["combineBackgrounds"] = {
            "Other":["SingleTop","DYJets","Other"] }

    boxesToUse = ["WJetControlRegion","TTJetsSingleLeptonControlRegion","WJetInvControlRegion",
            "VetoLeptonControlRegion","VetoTauControlRegion","GJetsInvControlRegion","TTJetsDileptonControlRegion"]
    if args.box is not None:
        boxesToUse = [args.box]

    #make output directory
    os.system('mkdir -p '+dirName)

    #estimate yields in signal region
    for boxName in boxesToUse:

        #apply options
        samplesToUse = SAMPLES[boxName]
        plotOpts["combineSamples"] = SAMPLES_REDUCED[boxName]

        print "\n---",boxName,"---"

        if boxName in xbinsSignal:
            unrollBins = (xbinsSignal[boxName]['0B'], colsSignal[boxName]['0B'])
        else:
            unrollBins = (None,None)
        inFile = dirName+'/controlHistograms'+boxMapping[boxName]+'.root'

        plotControlSampleHists(boxMapping[boxName], inFile, samples=samplesToUse, plotOpts=plotOpts, lumiMC=MCLUMI, lumiData=LUMI, boxName=boxName, debugLevel=debugLevel, printdir=dirName, unrollBins=unrollBins)
