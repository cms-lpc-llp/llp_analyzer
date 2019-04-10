import sys, os
import ROOT as rt

from macro import macro, razorWeights
from macro.razorAnalysis import Analysis, make_parser
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    parser = make_parser()
    parser.add_argument("--no-corr", dest="noCorr", action='store_true',
                                help="Don't apply scale factor correction")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    noCorr = args.noCorr
    boostCuts = not args.noBoostCuts

    #initialize
    plotOpts = { "comment":False, 'SUS15004CR':True } 
    regions = {}
    regionsOrder = []
    #define all tests
    jetsOrder = ["Inclusive","OneJet","DiJet","MultiJet"]
    jetsLimit = [(-1,-1),(1,1),(2,3),(4,-1)]
    for name,jets in zip(jetsOrder, jetsLimit):
        regionName = "OneLepton"+name+"ClosureTest"
        if noCorr:
            regionName += "_NoCorr"
        regions[regionName] = Analysis("SingleLepton",tag=tag,
                njetsMin=jets[0], njetsMax=jets[1], boostCuts=boostCuts)
        regionsOrder.append(regionName)
        for nb in range(3):
            regions[regionName+str(nb)+"B"] = Analysis("SingleLepton",tag=tag,
                    njetsMin=jets[0], njetsMax=jets[1], nbMin=nb, nbMax=nb, 
                    boostCuts=boostCuts)
            regionsOrder.append(regionName+str(nb)+"B")
    #add 3B test for MultiJet 
    regionName = "OneLeptonMultiJetClosureTest"
    if noCorr:
        regionName += "_NoCorr"
    regions[regionName+"3B"] = Analysis("SingleLepton",tag=tag,
            njetsMin=4, nbMin=3, nbMax=3, boostCuts=boostCuts)
    regionsOrder.append(regionName+"3B")

    sfHists = macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag), 
            processNames=regions.itervalues().next().samples, debugLevel=debugLevel)
    sfVars = ("MR","Rsq")
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    sfHists['NJetsTTJets'] = sfNJetsFile.Get("TTJetsScaleFactors")
    sfHists['NJetsWJets'] = sfNJetsFile.Get("WJetsScaleFactors")
    if noCorr:
        sfHists = {}

    for region in regionsOrder:
        analysis = regions[region]
        print "\nRegion:",region,"\n"
        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)
        #set up analysis
        (xbins,cols) = analysis.unrollBins
        auxSFs = razorWeights.getNJetsSFs(analysis)
        if noCorr:
            auxSFs = {}
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHists,
                sfVars=sfVars, printdir=outdir, auxSFs=auxSFs, btags=analysis.nbMin,
                debugLevel=debugLevel, noFill=args.noFill )
        #compute scale factors
        appendScaleFactors( region+"MR", hists, sfHists, lumiData=analysis.lumi, var="MR",
                debugLevel=debugLevel, signifThreshold=1.0, printdir=outdir )
        appendScaleFactors( region+"Rsq", hists, sfHists, lumiData=analysis.lumi, var="Rsq",
                debugLevel=debugLevel, signifThreshold=1.0, printdir=outdir )
        if not args.noSave:
            #export histograms
            macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                    outDir=outdir, debugLevel=debugLevel )
            #write out scale factors
            if not noCorr:
                outfile = rt.TFile("data/ScaleFactors/RazorMADD2015/RazorBTagClosureTests_%s.root"%(tag),
                        "UPDATE")
                print "Writing scale factor histogram",sfHists[region+"MR"].GetName(),"to file"
                outfile.cd()
                sfHists[region+"MR"].Write( sfHists[region+"MR"].GetName() )
                sfHists[region+"Rsq"].Write( sfHists[region+"Rsq"].GetName() )
                outfile.Close()
