import sys, os, argparse
import ROOT as rt

from macro import macro
from macro.razorAnalysis import Analysis
from macro.razorMacros import makeControlSampleHistsForAnalysis

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument("--muons", help="require muons", action='store_true')
    parser.add_argument("--electrons", help="require electrons", action='store_true')
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = "Razor2016"

    #initialize
    plotOpts = { 'comment':False }

    regionsOrder = ["0BSusySync", "1BSusySync", "2BSusySync"]
    regions = {
            "0BSusySync":Analysis("SusySync", tag=tag, nbMax=0),
            "1BSusySync":Analysis("SusySync", tag=tag, nbMin=1),
            "2BSusySync":Analysis("SusySync", tag=tag, nbMin=2),
            }

    for region in regionsOrder:
        analysis = regions[region]
        name = region
        #apply extra lepton cuts
        if args.muons:
            analysis.cutsData += " && abs(lep1Type) == 13"
            analysis.cutsMC += " && abs(lep1Type) == 13"
            name += "_Muons"
        if args.electrons:
            analysis.cutsData += " && abs(lep1Type) == 11"
            analysis.cutsMC += " && abs(lep1Type) == 11"
            name += "_Electrons"
        #make output directory
        outdir = 'Plots/'+tag+'/'+name
        os.system('mkdir -p '+outdir)
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, 
                printdir=outdir, debugLevel=debugLevel ) 
        #export histograms
        macro.exportHists(hists, outFileName='controlHistograms'+name+'.root', 
            outDir=outdir, debugLevel=debugLevel)
