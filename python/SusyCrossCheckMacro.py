import sys,os
import argparse
import ROOT as rt

#local imports
import macro.macro as macro
from macro.razorAnalysis import *
from macro.razorWeights import *
from macro.razorMacros import *

LUMI = 1264 #in /pb
MCLUMI = 1 

SAMPLES = ["TTV", "VV", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets"]

DIR_MC= "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC"
DIR_DATA= "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/Data"
PREFIX = "RazorInclusive"
DATA_NAME='SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates'
FILENAMES = {
        "TTJets"    : DIR_MC+"/"+PREFIX+"_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root",
        "WJets"     : DIR_MC+"/"+PREFIX+"_WJetsToLNu_HTBinned_1pb_weighted.root",
        "SingleTop" : DIR_MC+"/"+PREFIX+"_SingleTop_1pb_weighted.root",
        "VV" : DIR_MC+"/"+PREFIX+"_VV_1pb_weighted.root",
        "TTV" : DIR_MC+"/"+PREFIX+"_TTV_1pb_weighted.root",
        "DYJets"     : DIR_MC+"/"+PREFIX+"_DYJetsToLL_M-5to50_HTBinned_1pb_weighted.root",
        "ZInv"     : DIR_MC+"/"+PREFIX+"_ZJetsToNuNu_HTBinned_1pb_weighted.root",
        "Data"      : DIR_DATA+'/'+PREFIX+'_'+DATA_NAME+'.root',
        }

weightOpts = ["doNVtxWeights"]
shapeErrors = []
miscErrors = []

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
    weightHists = loadWeightHists(weightfilenames_DEFAULT, weighthistnames_DEFAULT, debugLevel)
    sfHists = {}

    #get scale factor histograms
    #sfHists = loadScaleFactorHists(processNames=SAMPLES, debugLevel=debugLevel)

    cutsForSusyCrossCheck = "HT > 500 && nSelectedJets >= 2 && met > 250 && leadingTightMuPt > 25 && mT < 120"
    cutsForSusyCrossCheckMC = appendTriggerCuts(cutsForSusyCrossCheck, singleLeptonTriggerNumsMC)
    cutsForSusyCrossCheckData = appendTriggerCuts(cutsForSusyCrossCheck, singleLeptonTriggerNumsData)

    binsForSusyCrossCheck = {
            "nSelectedJets":range(2,11),
            "nBTaggedJets":range(6),
            "met":range(250,800,10),
            "leadingTightMuPt":range(25,400,5),
            "HT":range(500,2000,30),
            }

    makeControlSampleHists("CrossCheck", 
            filenames=FILENAMES, samples=SAMPLES, 
            cutsMC=cutsForSusyCrossCheckMC, cutsData=cutsForSusyCrossCheckData, 
            bins=binsForSusyCrossCheck, lumiMC=MCLUMI, lumiData=LUMI, 
            weightHists=weightHists, sfHists=sfHists, treeName="RazorInclusive", 
            weightOpts=weightOpts, shapeErrors=shapeErrors, miscErrors=miscErrors,
            plotOpts={"logx":False, "ymin":0.1, "comment":False}, debugLevel=debugLevel)
