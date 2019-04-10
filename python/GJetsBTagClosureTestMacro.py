import sys, os
import copy
import ROOT as rt

from macro import macro, razorWeights
from macro.razorAnalysis import Analysis, make_parser
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors
import BTagClosureTestMacro as bclosure


if __name__ == "__main__":
    rt.gROOT.SetBatch()

    parser = make_parser()
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    boostCuts = not args.noBoostCuts

    plotOpts = { 'comment':False, "SUS15004CR":True }
    regions = {}
    regionsOrder = []

    #define all tests
    for name,jets in {
            "DiJet":(2,3),
            "MultiJet":(4,6),
            "SevenJet":(7,-1)
            }.iteritems():
        for nb in reversed(range(4)):
            regionName = 'GJetsInv'+name+str(nb)+'B'
            nbMax = nb
            if nb >= 2:
                nbMax = -1
            if nb > 2 and name == 'DiJet': 
                continue
            if nb < 3: # use 2b correction for 3b category
                regionsOrder.append(regionName)
                regions[regionName] = Analysis("GJetsInv", tag=tag,
                        njetsMin=jets[0], njetsMax=jets[1], nbMin=nb,
                        nbMax=nbMax, boostCuts=boostCuts)
                regionsOrder.append(regionName+'MRCorr')
                regions[regionName+'MRCorr'] = Analysis("GJetsInv", tag=tag,
                        njetsMin=jets[0], njetsMax=jets[1], nbMin=nb,
                        nbMax=nbMax, boostCuts=boostCuts)

    sfHists = macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag), 
            processNames=regions.itervalues().next().samples, debugLevel=debugLevel)
    sfVars = ("MR_NoPho","Rsq_NoPho")
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    sfHists['NJetsInv'] = sfNJetsFile.Get("GJetsInvScaleFactors")
    razorWeights.loadPhotonPurityHists(sfHists, tag, debugLevel)
    for region in regionsOrder:
        print "\nRegion:",region,"\n"
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)
        analysis = regions[region]
        auxSFs = razorWeights.getNJetsSFs(analysis, jetName=analysis.jetVar)
        auxSFs = razorWeights.addBTagSFs(analysis, auxSFs, var='MR_NoPho', gjets=True)
        bclosure.adjustForRegion(analysis, sfHists, auxSFs, gjets=True)
        razorWeights.getPhotonPuritySFs(auxSFs)
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHists,
                sfVars=sfVars, printdir=outdir, auxSFs=auxSFs, btags=analysis.nbMin,
                dataDrivenQCD=True, debugLevel=debugLevel, noFill=args.noFill )

        sfHistName = bclosure.getSFHistName(analysis, gjets=True)
        sfHistsTmp = copy.copy(sfHists)
        if 'MRCorr' in region:
            sfHistName = sfHistName.replace('MR', 'Rsq')
            appendScaleFactors( 'GJetsInv', hists, sfHistsTmp, lumiData=analysis.lumi, 
                    var="Rsq_NoPho", debugLevel=debugLevel, printdir=outdir )
        else:
            appendScaleFactors( 'GJetsInv', hists, sfHistsTmp, lumiData=analysis.lumi, 
                    var="MR_NoPho", debugLevel=debugLevel, printdir=outdir )
        sfHists[sfHistName] = sfHistsTmp['GJetsInv']
        if not args.noSave:
            macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                    outDir=outdir, debugLevel=debugLevel )
            #write out scale factors
            outfile = rt.TFile(bclosure.getOutputFilename(tag), 'UPDATE')
            histToWrite = sfHists[sfHistName]
            histToWrite.SetName(sfHistName)
            print "Writing scale factor histogram",histToWrite.GetName(),"to file"
            outfile.cd()
            histToWrite.Write( histToWrite.GetName() )
            outfile.Close()
