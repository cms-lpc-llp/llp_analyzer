import sys, os, copy
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

    plotOpts = { "comment":False, 'SUS15004CR':True } 
    wjetsRegionNames = ["WJetsSingleLeptonInvNLO", 
            "WJetsSingleLeptonInvNLOForNJets"] 
    dyjetsRegionNames = ["DYJetsDileptonInvNLO", "DYJetsDileptonInvNLODiJet", 
            "DYJetsDileptonInvNLOMultiJet"]
    regionsOrder = wjetsRegionNames + dyjetsRegionNames
    regions = {
            "WJetsSingleLeptonInvNLO":Analysis("WJetsSingleLeptonInv", 
                nbMax=0, tag=tag, boostCuts=boostCuts),
            "WJetsSingleLeptonInvNLOForNJets":Analysis("WJetsSingleLeptonInv", 
                nbMax=0, tag=tag, boostCuts=boostCuts),
            "DYJetsDileptonInvNLO":Analysis("DYJetsDileptonInv",
                tag=tag,boostCuts=boostCuts),
            "DYJetsDileptonInvNLODiJet":Analysis("DYJetsDileptonInv",
                tag=tag, njetsMin=2, njetsMax=3, boostCuts=boostCuts),
            "DYJetsDileptonInvNLOMultiJet":Analysis("DYJetsDileptonInvMultiJet",
                tag=tag, njetsMin=4, boostCuts=boostCuts),
            }
    sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag)
    sfHists = macro.loadScaleFactorHists( sfFilename=sfFilename,
            processNames=regions["DYJetsDileptonInvNLO"].samples, 
            debugLevel=debugLevel )
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    sfHists['NJetsTTJets'] = sfNJetsFile.Get("TTJetsScaleFactors")
    sfHists['NJetsWJets'] = sfNJetsFile.Get("WJetsScaleFactors")
    sfVars = { "WJets":("MR","Rsq"), "TTJets":("MR","Rsq"), 
            "DYJetsInv":("MR_NoZ","Rsq_NoZ"), "WJetsInv":("MR_NoW", "Rsq_NoW") }

    # We need to recompute the W+jets invisible scale factors 
    # using the NLO sample.
    for reg in wjetsRegionNames:
        regions[reg].filenames['WJetsInv'] = regions[reg].filenames[
                'WJetsInv'].replace('WJets_', 'WJetsPtBinned_')
    for reg in dyjetsRegionNames:
        regions[reg].filenames['DYJetsInv'] = regions[reg].filenames[
                'DYJetsInv'].replace('DYJets_', 'DYJetsPtBinned_')

    if not args.noSave:
        outfile = rt.TFile(
            "data/ScaleFactors/RazorMADD2015/RazorDYJetsDileptonInvNLOCrossCheck_%s.root"%(tag), "RECREATE")

    for region in regionsOrder:
        analysis = regions[region]
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)

        if region in wjetsRegionNames:
            jetName = 'NJets_NoW'
        else:
            jetName = 'NJets_NoZ'
        auxSFs = razorWeights.getNJetsSFs(analysis,jetName=jetName)
        if 'WJetsInv' in auxSFs:
            # W+jets invisible njets SFs are not computed yet
            auxSFs['WJetsInv'] = {}
        sfHistsToUse = sfHists

        hists = makeControlSampleHistsForAnalysis(analysis, plotOpts=plotOpts, 
                sfHists=sfHistsToUse, sfVars=sfVars, printdir=outdir, 
                auxSFs=auxSFs, debugLevel=debugLevel, noFill=args.noFill)

        if region == 'WJetsSingleLeptonInvNLO':
            xbins, cols = analysis.unrollBins
            appendScaleFactors('WJetsInv', hists, sfHists, lumiData=analysis.lumi,
                    th2PolyXBins=xbins, th2PolyCols=cols, debugLevel=debugLevel,
                    var=sfVars['WJetsInv'], printdir=outdir)
            sfHists['DYJetsInv'] = sfHists['WJetsInv']
            histToWrite = sfHists['WJetsInv']
        elif region == 'WJetsSingleLeptonInvNLOForNJets':
            sfHistsCopy = sfHists.copy()
            appendScaleFactors('WJetsInv', hists, sfHistsCopy, 
                    lumiData=analysis.lumi, debugLevel=debugLevel,
                    var=jetName, printdir=outdir)
            sfHists['NJetsInv'] = sfHistsCopy['WJetsInv']
            histToWrite = sfHistsCopy['WJetsInv']
        else:
            sfHistsCopy = sfHists.copy()
            del sfHistsCopy['DYJetsInv']
            appendScaleFactors('DYJetsInv', hists, sfHistsCopy,
                    lumiData=analysis.lumi, debugLevel=debugLevel,
                    var=sfVars['DYJetsInv'], signifThreshold=1.0,
                    printdir=outdir)
            histToWrite = sfHistsCopy['DYJetsInv']
        if not args.noSave:
            print "Writing histogram {} to file".format(
                    histToWrite.GetName())
            outfile.cd()
            histToWrite.Write(region+"ScaleFactors")
            macro.exportHists(hists, outDir=outdir, debugLevel=debugLevel,
                    outFileName='controlHistograms'+region+'.root')

    if not args.noSave:
        outfile.Close()
