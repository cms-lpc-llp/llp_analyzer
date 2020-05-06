import sys, os
import ROOT as rt

from macro import macro, razorWeights
from macro.razorAnalysis import Analysis, make_parser
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    parser = make_parser()
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    boostCuts = not args.noBoostCuts

    #initialize
    plotOpts = { 'comment':False, "SUS15004CR":True }
    regions = {}
    regionsOrder = []
    #define all tests
    for name,jets in {"DiJet":(2,3),"MultiJet":(4,-1)}.iteritems():
        regionName = "GJetsInv"+name+"ClosureTest"
        regions[regionName] = Analysis("GJetsInv",tag=tag,
                njetsMin=jets[0], njetsMax=jets[1],boostCuts=boostCuts)
        regionsOrder.append(regionName)
        maxB = 3
        if name == 'DiJet':
            maxB = 2
        for nb in range(maxB+1):
            nbMax = nb
            if nb == maxB:
                nbMax = -1
            regionName = "GJetsInv"+name+"ClosureTest"+str(nb)+"B"
            regionsOrder.append(regionName)
            regions[regionName] = Analysis("GJetsInv",tag=tag,
                    njetsMin=jets[0], njetsMax=jets[1], nbMin=nb,
                    nbMax=nbMax, boostCuts=boostCuts)

    sfHists = macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s_Uncorr.root"%(tag), 
            processNames=regions.itervalues().next().samples, debugLevel=debugLevel)
    sfVars = ("MR_NoPho","Rsq_NoPho")
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    sfHists['NJetsInv'] = sfNJetsFile.Get("GJetsInvScaleFactors")
    razorWeights.loadPhotonPurityHists(sfHists, tag, debugLevel)
    for region in regionsOrder:
        analysis = regions[region]
        print "\nRegion:",region,"\n"
        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)
        #set up analysis
        (xbins,cols) = analysis.unrollBins
        auxSFs = razorWeights.getNJetsSFs(analysis, jetName=analysis.jetVar)
        razorWeights.getPhotonPuritySFs(auxSFs)
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHists,
                sfVars=sfVars, printdir=outdir, auxSFs=auxSFs, btags=analysis.nbMin,
                dataDrivenQCD=True, debugLevel=debugLevel, noFill=args.noFill )
        #compute scale factors
        appendScaleFactors( region+"MR", hists, sfHists, lumiData=analysis.lumi, 
                var="MR_NoPho", debugLevel=debugLevel, signifThreshold=1.0, printdir=outdir )
        appendScaleFactors( region+"Rsq", hists, sfHists, lumiData=analysis.lumi, 
                var="Rsq_NoPho", debugLevel=debugLevel, signifThreshold=1.0, printdir=outdir )
        #export histograms
        macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                outDir=outdir, debugLevel=debugLevel )
        if not args.noSave:
            #write out scale factors
            outfile = rt.TFile(
                    "data/ScaleFactors/RazorMADD2015/RazorGJetsBTagClosureTests_%s.root"%(tag),
                    "UPDATE")
            print "Writing scale factor histogram",sfHists[region+"MR"].GetName(),"to file"
            print "Writing scale factor histogram",sfHists[region+"Rsq"].GetName(),"to file"
            outfile.cd()
            sfHists[region+"MR"].Write( sfHists[region+"MR"].GetName() )
            sfHists[region+"Rsq"].Write( sfHists[region+"Rsq"].GetName() )
            outfile.Close()
