### fit MC cocktail to validate the fit function

import sys,os
import argparse
import ROOT as rt

#local imports
from framework import Config
import macro.macro as macro
from macro.razorAnalysis import *
from macro.razorWeights import *
from macro.razorMacros import *

LUMI = 2000 #in /pb
MCLUMI = 1 

SAMPLES = ["Other", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets"]

DIR_MC = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/forfit"
PREFIX = "RazorInclusive"
FILENAMES_HADRONIC = {
        "TTJets"    : DIR_MC+"/"+PREFIX+"_TTJets_Madgraph_Leptonic_1pb_weighted_RazorSkim.root",
        "WJets"     : DIR_MC+"/"+PREFIX+"_WJetsToLNu_AlternativeHTBinned_1pb_weighted_RazorSkim.root",
        "SingleTop" : DIR_MC+"/"+PREFIX+"_ST_1pb_weighted_RazorSkim.root",
        "Other" : DIR_MC+"/"+PREFIX+"_Other_1pb_weighted_RazorSkim.root",
        "DYJets"     : DIR_MC+"/"+PREFIX+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
        "ZInv"     : DIR_MC+"/"+PREFIX+"_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root",
        }

FILENAMES_LEPTONIC = {
        "TTJets"    : DIR_MC+"/"+PREFIX+"_TTJets_Madgraph_Leptonic_1pb_weighted_RazorSkim.root",
        "WJets"     : DIR_MC+"/"+PREFIX+"_WJetsToLNu_AlternativeHTBinned_1pb_weighted_RazorSkim.root",
        "SingleTop" : DIR_MC+"/"+PREFIX+"_ST_1pb_weighted_RazorSkim.root",
        "Other" : DIR_MC+"/"+PREFIX+"_Other_1pb_weighted_RazorSkim.root",
        "DYJets"     : DIR_MC+"/"+PREFIX+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
        "ZInv"     : DIR_MC+"/"+PREFIX+"_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root",
        }

FILENAMES={
        "MultiJet":FILENAMES_HADRONIC,
        "MuMultiJet":FILENAMES_LEPTONIC,
        "EleMultiJet":FILENAMES_LEPTONIC,
        }

config = "config/run2_20151108_Preapproval.config"
FIT_DIR = "FitPlots_2fb"
TOYS_FILES = {
        "MultiJet":FIT_DIR+"/toys_Bayes_noStat_MultiJet.root",
        "MuMultiJet":FIT_DIR+"/toys_Bayes_noStat_MuMultiJet.root",
        "EleMultiJet":FIT_DIR+"/toys_Bayes_noStat_EleMultiJet.root",
        }

cfg = Config.Config(config)
binsMRHad = cfg.getBinning("MultiJet")[0]
binsRsqHad = cfg.getBinning("MultiJet")[1]
hadronicBinning = { "MR":binsMRHad, "Rsq":binsRsqHad }
binsMRLep = cfg.getBinning("MuMultiJet")[0]
binsRsqLep = cfg.getBinning("MuMultiJet")[1]
leptonicBinning = { "MR":binsMRLep, "Rsq":binsRsqLep }
binning = { "MultiJet":hadronicBinning, "MuMultiJet":leptonicBinning, "EleMultiJet":leptonicBinning}

weightOpts = []
shapeErrors = []
miscErrors = []

dirName = "ForJamboree_2fb"

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

    #make output directory
    os.system('mkdir -p '+dirName)

    #estimate yields in leptonic signal region
    #for lepType in ["Ele"]:
    #for lepType in ["Mu"]:
    for lepType in [""]:
    #for lepType in ["", "Mu", "Ele"]:
        for jets in ["MultiJet"]:
            boxName = lepType+jets
            #btaglist = [0]
            btaglist = [0,1,2,3]
            for btags in btaglist:
                print "\n---",boxName,"Box,",btags,"B-tags ---"
                #get correct cuts string
                thisBoxCuts = razorCuts[boxName]
                if btags < len(btaglist)-1:
                    thisBoxCuts += " && nBTaggedJets == "+str(btags)
                else:
                    thisBoxCuts += " && nBTaggedJets >= "+str(btags)

                if len(btaglist) > 1:
                    extboxName = boxName+str(btags)+"BTag"
                    nBtags = btags
                else:
                    extboxName = boxName
                    nBtags = -1
                #check fit file and create if necessary
                if not os.path.isfile(TOYS_FILES[boxName]):
                    print "Fit file",TOYS_FILES[boxName],"not found, trying to recreate it"
                    runFitAndToysMC(FIT_DIR, boxName, LUMI, [FILENAMES[boxName][x] for x in FILENAMES[boxName]], DIR_MC, config=config, sideband=True, numToys=4000, noStat=True)
                    #check
                    if not os.path.isfile(TOYS_FILES[boxName]):
                        print "Error creating fit file",TOYS_FILES[boxName]
                        sys.exit()
                makeControlSampleHists(extboxName, 
                        filenames=FILENAMES[boxName], samples=SAMPLES, 
                        cutsMC=thisBoxCuts, cutsData=thisBoxCuts, 
                        bins=binning[boxName], lumiMC=MCLUMI, lumiData=LUMI, 
                        weightHists=weightHists, sfHists=sfHists, treeName="RazorInclusive", 
                        weightOpts=weightOpts, shapeErrors=shapeErrors, miscErrors=miscErrors,
                        fitToyFiles=TOYS_FILES, boxName=boxName, 
                        btags=nBtags, debugLevel=debugLevel, printdir=dirName)

