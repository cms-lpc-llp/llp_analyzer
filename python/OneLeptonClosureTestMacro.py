import sys, os
import ROOT as rt

import macro.razorAnalysis as razor
from macro import macro, razorWeights
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors

def isSevenJet(analysis):
    return analysis.njetsMin >= 7 or 'SevenJet' in analysis.region

def adjustForRegion(analysis, gjets=False):
    # Reduce the SevenJet binning due to low stats
    if isSevenJet(analysis):
        if gjets:
            if 'MR_NoPho' in analysis.binning:
                print "Using reduced binning for seven jet category"
                analysis.binning['MR_NoPho'] = [400, 600, 900, 4000]
                analysis.binning['Rsq_NoPho'] = [0.25, 0.30, 0.41, 1.5]
        elif analysis.jetVar == 'NJets40': # the condition is to avoid doing this for signal region
            print "Using reduced binning for seven jet category"
            analysis.binning['MR'] = [300, 500, 700, 900, 4000]
            analysis.binning['Rsq'] = [0.15, 0.20, 0.25, 0.30, 0.41, 1.5]
            analysis.unrollBins = (analysis.binning['MR'],
                    [[0.15, 0.20, 0.25, 0.30, 0.41, 1.5],
                     [0.15, 0.20, 0.25, 0.30, 0.41, 1.5],
                     [0.15, 0.20, 0.25, 0.30, 1.5],
                     [0.15, 0.20, 1.5]])


if __name__ == "__main__":
    rt.gROOT.SetBatch()

    parser = razor.make_parser()
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    boostCuts = not args.noBoostCuts

    plotOpts = { "comment":False, 'SUS15004CR':True } 
    regions = {}
    regionsOrder = []

    #define all tests
    jetsOrder = ["SevenJet", "MultiJet", "DiJet"]
    jetsLimit = [(7,-1),(4,6),(2,3)]
    for name,jets in zip(jetsOrder, jetsLimit):
        regionName = 'OneLepton'+name
        regions[regionName] = razor.Analysis("SingleLepton",tag=tag,
                njetsMin=jets[0], njetsMax=jets[1], boostCuts=boostCuts)
        regionsOrder.append(regionName)

    sfHists = macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag), 
            processNames=regions.itervalues().next().samples, debugLevel=debugLevel)
    sfVars = ("MR","Rsq")
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    sfHists['NJetsTTJets'] = sfNJetsFile.Get("TTJetsScaleFactors")
    sfHists['NJetsWJets'] = sfNJetsFile.Get("WJetsScaleFactors")

    for region in regionsOrder:
        print "\nRegion:",region,"\n"
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)
        analysis = regions[region]
        auxSFs = razorWeights.getNJetsSFs(analysis)
        adjustForRegion(analysis)
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHists,
                sfVars=sfVars, printdir=outdir, auxSFs=auxSFs, 
                debugLevel=debugLevel, noFill=args.noFill )
        if not args.noSave:
            macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                    outDir=outdir, debugLevel=debugLevel )
