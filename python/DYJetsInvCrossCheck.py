import sys, os, copy
import ROOT as rt

from macro import macro, razorWeights
from macro.razorAnalysis import Analysis, make_parser
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors
import BTagClosureTestMacro as bclosure

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    parser = make_parser()
    parser.add_argument("--closure", action="store_true", help="include uncertainties from scale factor cross check")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    closure = args.closure
    boostCuts = not args.noBoostCuts

    #initialize
    plotOpts = { "comment":False, 'SUS15004CR':True } 
    #Process inclusive sample twice; the first pass will compute the overall normalization 
    #and the second pass will be a rerun with the corrected normalization
    if closure:
        regionsOrder = []
        for jets in ['SevenJet', 'MultiJet', 'DiJet', '']:
            regionsOrder.append('DYJetsDileptonInv'+jets)
    else:
        regionsOrder = [
                "DYJetsDileptonInvUncorr", 
                "DYJetsDileptonInv", 
                "DYJetsDileptonInvDiJet", 
                "DYJetsDileptonInvMultiJet", 
                "DYJetsDileptonInvSevenJet"
                ]
    regions = {
            "DYJetsDileptonInvUncorr":Analysis("DYJetsDileptonInv",
                tag=tag,boostCuts=boostCuts),
            "DYJetsDileptonInv":Analysis("DYJetsDileptonInv",
                tag=tag,boostCuts=boostCuts),
            "DYJetsDileptonInvDiJet":Analysis("DYJetsDileptonInv",
                tag=tag,njetsMin=2,njetsMax=3, boostCuts=boostCuts),
            "DYJetsDileptonInvMultiJet":Analysis("DYJetsDileptonInvMultiJet",
                tag=tag,njetsMin=4,njetsMax=6, boostCuts=boostCuts),
            "DYJetsDileptonInvSevenJet":Analysis("DYJetsDileptonInvMultiJet",
                tag=tag,njetsMin=7, boostCuts=boostCuts),
            "DYJetsDileptonInvNoSFs":Analysis("DYJetsDileptonInv",
                tag=tag,boostCuts=boostCuts),
            }
    sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag)
    #make two dictionaries of scale factor histograms, one with GJets and one with WJets corrections
    sfHists = macro.loadScaleFactorHists( sfFilename=sfFilename,
            processNames=regions["DYJetsDileptonInvDiJet"].samples, 
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
    if not args.noSave:
        outfile = rt.TFile(
            "data/ScaleFactors/RazorMADD2015/RazorDYJetsDileptonInvCrossCheck_%s.root"%(tag), "RECREATE")
    #optionally inflate scale factor uncertainties to cover difference between G+jets and W+jets SFs
    if args.closure:
        sfFile = rt.TFile.Open(sfFilename)
        downHist = sfFile.Get("GJetsInvScaleFactors_Down")
        for bn in range(sfHists["DYJetsInv"].GetNumberOfBins()+1):
            err = sfHists["DYJetsInv"].GetBinError(bn)
            sysErr = sfHists["DYJetsInv"].GetBinContent(bn) - downHist.GetBinContent(bn)
            newErr = ( err*err + sysErr*sysErr )**(0.5)
            sfHists["DYJetsInv"].SetBinError(bn-1, newErr) # adjust bin number by 1 to account for bug in ROOT
            print "Increasing error on DYJets scale factor bin",bn,"from",err,"to",sfHists["DYJetsInv"].GetBinError(bn)
    
    updateNorm = (not closure)
    for region in regionsOrder:
        analysis = regions[region]
        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        if closure:
            outdir += '_Closure'
        os.system('mkdir -p '+outdir)
        #prepare analysis
        auxSFs = razorWeights.getNJetsSFs(analysis,jetName='NJets_NoZ')
        auxSFs = razorWeights.addAllBTagSFs(analysis, auxSFs)
        auxSFs = razorWeights.addAllBTagSFs(analysis, auxSFs, var='MR_NoZ', 
                gjets=True)
        bclosure.adjustForRegionBInclusive(analysis, sfHists, auxSFs)
        bclosure.adjustForRegionBInclusive(analysis, sfHists, auxSFs, gjets=True)
        #use the correct set of scale factors
        if 'NoSFs' in region:
            sfHistsToUse = {}
            auxSFs = { proc:{} for proc in auxSFs }
        else:
            sfHistsToUse = sfHists
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHistsToUse,
            sfVars=sfVars, printdir=outdir, auxSFs=auxSFs, debugLevel=debugLevel, noFill=args.noFill )
        for v in ['MR', 'Rsq']:
            tmpSFHists = copy.copy(sfHists)
            del tmpSFHists["DYJetsInv"]
            appendScaleFactors("DYJetsInv", hists, tmpSFHists, lumiData=analysis.lumi, 
                debugLevel=debugLevel, var=v+'_NoZ', printdir=outdir)
            histName = region+v+'ScaleFactors'
            tmpSFHists['DYJetsInv'].SetName(histName)
            if not args.noSave:
                print "Writing histogram",tmpSFHists["DYJetsInv"].GetName(),"to file"
                outfile.cd()
                tmpSFHists["DYJetsInv"].Write(histName)

        #in the first pass, update the normalization of the G+jets scale factor histogram
        if updateNorm:
            dataNorm = hists["Data"]["1"].GetBinContent(1)
            dataNormErr = hists["Data"]["1"].GetBinError(1)
            mcNorm = 0
            mcNormErr = 0
            for proc in analysis.samples:
                mcNorm += hists[proc]["1"].GetBinContent(1)
                mcNormErr += (hists[proc]["1"].GetBinError(1))**2
            mcNormErr = mcNormErr**(0.5)
            normUpdate = dataNorm / mcNorm
            normUpdateErr = ( (dataNormErr/mcNorm)**2 + (dataNorm*mcNormErr/(mcNorm*mcNorm))**2 )**(0.5)

            if not args.noSave:
                sfFilenameUncorr="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s_Uncorr.root"%(tag)
                print "Writing old G+jets scale factor histogram to",sfFilenameUncorr
                sfFileUncorr = rt.TFile.Open(sfFilenameUncorr, "RECREATE")
                sfHistUncorr = sfHists["DYJetsInv"].Clone()
                sfFileUncorr.WriteTObject(sfHistUncorr, "GJetsInvScaleFactors", "WriteDelete")
                sfFileUncorr.Close()
            print "Scaling G+jets scale factor histogram by",normUpdate,"( = %.3f / %.3f )"%(dataNorm,mcNorm)
            print "Propagating error on scale factor:",normUpdateErr
            for bn in range(sfHists["DYJetsInv"].GetNumberOfBins()+1):
                sfHists["DYJetsInv"].SetBinContent( bn, sfHists["DYJetsInv"].GetBinContent(bn)*normUpdate )
                sfHists["DYJetsInv"].SetBinError( bn-1, sfHists["DYJetsInv"].GetBinError(bn)*normUpdate )
                sfHists["DYJetsInv"].SetBinError( bn-1, 
                        ( ( sfHists["DYJetsInv"].GetBinError(bn))**2 + normUpdateErr*normUpdateErr )**(0.5) )
            sfHistCorr = sfHists["DYJetsInv"].Clone()
            if not args.noSave:
                #Write the corrected G+jets scale factor file over the old one
                sfFile = rt.TFile.Open(sfFilename, "UPDATE")
                sfFile.WriteTObject(sfHistCorr, "GJetsInvScaleFactors", "WriteDelete")
                sfFile.Close()
            updateNorm = False

        if not args.noSave:
            #export histograms
            macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                    outDir=outdir, debugLevel=debugLevel )

    if not args.noSave:
        outfile.Close()

