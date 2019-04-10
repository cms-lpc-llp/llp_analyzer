import sys, os, copy
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

    #initialize
    plotOpts = { "comment":False, 'SUS15004CR':True } 
    regionsOrder = ["DYJetsDileptonInvBCheck"] 
    regions = {
            "DYJetsDileptonInvBCheck":Analysis("DYJetsDileptonInv",
                tag=tag, boostCuts=boostCuts) }
    sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag)
    sfHists = macro.loadScaleFactorHists( sfFilename=sfFilename,
            processNames=regions["DYJetsDileptonInvBCheck"].samples, 
            scaleFactorNames={ "DYJetsInv":"GJetsInv" }, debugLevel=debugLevel )
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    for d in [sfHists]:
        d['NJetsTTJets'] = sfNJetsFile.Get("TTJetsScaleFactors")
        d['NJetsWJets'] = sfNJetsFile.Get("WJetsScaleFactors")
        d['NJetsInv'] = sfNJetsFile.Get("GJetsInvScaleFactors")
        d['NJetsWJetsInv'] = sfNJetsFile.Get("WJetsInvScaleFactors")
        bclosure.loadScaleFactors(d, tag=tag)
        bclosure.loadScaleFactors(d, tag=tag, gjets=True)
    sfVars = { "WJets":("MR","Rsq"), "TTJets":("MR","Rsq"), "DYJetsInv":("MR_NoZ","Rsq_NoZ") }
    for region in regionsOrder:
        analysis = regions[region]
        analysis.cutsMC = analysis.cutsMC.replace(' && NBJetsMedium == 0', '')
        analysis.cutsData = analysis.cutsData.replace(' && NBJetsMedium == 0', '')
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)
        auxSFs = razorWeights.getNJetsSFs(analysis,jetName='NJets_NoZ')
        auxSFs = razorWeights.addAllBTagSFs(analysis, auxSFs)
        auxSFs = razorWeights.addAllBTagSFs(analysis, auxSFs, var='MR_NoZ', gjets=True)
        bclosure.adjustForRegionBInclusive(analysis, sfHists, auxSFs)
        bclosure.adjustForRegionBInclusive(analysis, sfHists, auxSFs, gjets=True)
        sfHistsToUse = sfHists
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHistsToUse,
            sfVars = sfVars, printdir=outdir, auxSFs=auxSFs, debugLevel=debugLevel, noFill=args.noFill )
        sfHistsTmp = copy.copy(sfHists)
        appendScaleFactors("DYJetsInv", hists, sfHistsTmp, lumiData=analysis.lumi,
                debugLevel=debugLevel, var='NBJetsMedium', printdir=outdir)
        if not args.noSave:
            macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                    outDir=outdir, debugLevel=debugLevel )
            outfile = rt.TFile(
                    "data/ScaleFactors/RazorMADD2015/RazorDYJetsDileptonInvBCheck_{}.root".format(
                        tag), "RECREATE")
            sfHistsTmp["DYJetsInv"].Write(region+"ScaleFactors")
