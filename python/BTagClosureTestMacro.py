import sys, os
import ROOT as rt

import macro.razorAnalysis as razor
from macro import macro, razorWeights
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors

def getOutputFilename(tag='Razor2016_MoriondRereco'):
    return "data/ScaleFactors/RazorMADD2015/RazorBTagScaleFactors_%s.root"%(tag)

def loadScaleFactors(sfHists={}, tag='Razor2016_MoriondRereco', gjets=False):
    f = rt.TFile.Open(getOutputFilename(tag))
    sfKey = 'MR'
    if gjets:
        sfKey = 'MRInv'
    for jets in ['DiJet', 'MultiJet', 'SevenJet']:
        for btags in range(4):
            if (jets == 'DiJet' or gjets) and btags > 2:
                continue
            name = jets+str(btags)+'B'
            h = f.Get(sfKey+name)
            assert(h)
            h.SetDirectory(0)
            sfHists[sfKey+jets+str(btags)+'B'] = h
    return sfHists

def isMultiJet(analysis):
    return (analysis.njetsMin >= 4 and analysis.njetsMin < 7) or 'MultiJet' in analysis.region

def isSevenJet(analysis):
    return analysis.njetsMin >= 7 or 'SevenJet' in analysis.region

def getSFHistName(analysis, gjets=False):
    jetName = 'DiJet'
    if isMultiJet(analysis):
        jetName = 'MultiJet'
    if isSevenJet(analysis):
        jetName = 'SevenJet'
    nbtags = max(0, analysis.nbMin)
    sfKey = 'MR'
    if gjets:
        sfKey = 'MRInv'
        if nbtags > 2:
            nbtags = 2
    return sfKey+jetName+str(nbtags)+'B'

def adjustForRegion(analysis, sfHists, auxSFs, gjets=False):
    """
    Adjusts scale factor dictionaries so that the correct
    corrections are applied for the current region
    """
    sfKey = 'MR'
    if gjets:
        sfKey = 'MRInv'
    if sfKey in sfHists:
        # This avoids using last region's scale factors
        del sfHists[sfKey]
    sfHistName = getSFHistName(analysis, gjets=gjets)
    if sfHistName in sfHists:
        # This sets up the correct reweighting histogram
        sfHists[sfKey] = sfHists[sfHistName]
        print "Using histogram {} for {}".format(sfHistName, sfKey)
    else:
        # This avoids trying to apply SFs we don't have
        for proc, sfs in auxSFs.iteritems():
            if sfKey in sfs:
                print "Omitting {} scale factors for {}".format(sfKey, proc)
                del sfs[sfKey]
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

def adjustForRegionBInclusive(analysis, sfHists, auxSFs, gjets=False):
    """
    Same as adjustForRegion, but sets up all b-tag 
    histograms simulataneously.  Use to process a sample
    that contains events of different b-tag multiplicities.
    """
    sfKey = 'MR'
    if gjets:
        sfKey = 'MRInv'
    nbMax = 2
    jets = 'DiJet'
    if isMultiJet(analysis):
        if not gjets:
            nbMax = 3
        jets = 'MultiJet'
    if isSevenJet(analysis):
        if not gjets:
            nbMax = 3
        jets = 'SevenJet'
    for nb in range(nbMax + 1):
        thisSFKey = sfKey+str(nb)+'B'
        if thisSFKey in sfHists:
            # This avoids using last region's scale factors
            del sfHists[thisSFKey]
        sfHistName = sfKey+jets+str(nb)+'B'
        if sfHistName in sfHists:
            # This sets up the correct reweighting histogram
            sfHists[thisSFKey] = sfHists[sfHistName]
            print "Using histogram {} for {}".format(sfHistName, thisSFKey)
        else:
            # This avoids trying to apply SFs we don't have
            for proc, sfs in auxSFs.iteritems():
                if thisSFKey in sfs:
                    print "Omitting {} scale factors for {}".format(thisSFKey, proc)
                    del sfs[thisSFKey]
    if isSevenJet(analysis):
        print "Using reduced binning for seven jet category"
        if "MR_NoZ" in analysis.binning:
            analysis.binning["MR_NoZ"] = [400, 4000]
            analysis.binning["Rsq_NoZ"] = [0.25, 1.5]
        if "MR" in analysis.binning:
            analysis.binning["MR"] = [300, 4000]
            analysis.binning["Rsq"] = [0.15, 1.5]


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
        for nb in reversed(range(4)):
            if nb >= 3 and name == 'DiJet': 
                continue
            regionName = 'OneLepton'+name+str(nb)+'B'
            regions[regionName] = razor.Analysis("SingleLepton",tag=tag,
                    njetsMin=jets[0], njetsMax=jets[1], nbMin=nb, nbMax=nb, 
                    boostCuts=boostCuts)
            regionsOrder.append(regionName)
            regions[regionName+'MRCorr'] = razor.Analysis("SingleLepton",tag=tag,
                    njetsMin=jets[0], njetsMax=jets[1], nbMin=nb, nbMax=nb, 
                    boostCuts=boostCuts)
            regionsOrder.append(regionName+'MRCorr')

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
        auxSFs = razorWeights.addBTagSFs(analysis, auxSFs)
        adjustForRegion(analysis, sfHists, auxSFs)
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHists,
                sfVars=sfVars, printdir=outdir, auxSFs=auxSFs, btags=analysis.nbMin,
                debugLevel=debugLevel, noFill=args.noFill )

        sfHistName = getSFHistName(analysis)
        sfProcs = ['TTJets', 'WJets']
        sfHistsTmp = sfHists.copy()
        if 'MRCorr' in region:
            sfHistName = sfHistName.replace('MR', 'Rsq')
            appendScaleFactors( sfProcs, hists, sfHistsTmp, 
                    lumiData=analysis.lumi, var="Rsq", useUncertainty=True,
                    debugLevel=debugLevel, printdir=outdir )
        else:
            appendScaleFactors( sfProcs, hists, sfHistsTmp, 
                    lumiData=analysis.lumi, var="MR",
                    debugLevel=debugLevel, printdir=outdir )
        sfHists[sfHistName] = sfHistsTmp['_'.join(sfProcs)]
        sfHists[sfHistName].SetName(sfHistName)
        if not args.noSave:
            macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                    outDir=outdir, debugLevel=debugLevel )
            outfile = rt.TFile(getOutputFilename(tag), "UPDATE")
            histToWrite = sfHists[sfHistName]
            print "Writing scale factor histogram",histToWrite.GetName(),"to file"
            outfile.cd()
            histToWrite.Write( histToWrite.GetName() )
            outfile.Close()
