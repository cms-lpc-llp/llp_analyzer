### This script is like ComputeScaleFactorsNJetsCorrection.py except that it 
### computes the NJets correction separately for TT+jets and W+jets
import sys, os
import ROOT as rt

from macro import macro, razorWeights
from macro.razorAnalysis import Analysis, make_parser
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    parser = make_parser()
    parser.add_argument("--tight-cuts", action='store_true', dest='tightCuts',
            help="Cut on delta phi and number of 80 GeV jets, and cut tighter on MR and Rsq")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    boostCuts = not args.noBoostCuts

    plotOpts = { 'comment':False, "SUS15004CR":True }
    regionsOrder = ["GJetsInvForNJets", "TTJetsForNJets", 
            "WJetsForNJets", "WJetsInvForNJets"]
    regions = {
            "GJetsInvForNJets":Analysis("GJetsInv",tag=tag,boostCuts=boostCuts),
            "TTJetsForNJets":Analysis("TTJetsSingleLepton",tag=tag,nbMin=1,boostCuts=boostCuts),
            "WJetsForNJets":Analysis("WJetsSingleLepton",tag=tag,nbMax=0,boostCuts=boostCuts),
            "WJetsInvForNJets":Analysis("WJetsSingleLeptonInv",tag=tag,nbMax=0,boostCuts=boostCuts),
            }
    sfVars = {
            "GJetsInvForNJets":("MR_NoPho","Rsq_NoPho"),
            "TTJetsForNJets":("MR","Rsq"),
            "WJetsForNJets":("MR","Rsq"),
            "WJetsInvForNJets":("MR_NoW","Rsq_NoW"),
            }
    sfHists = macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag), 
            processNames=regions["TTJetsForNJets"].samples+["GJetsInv",'WJetsInv'], 
            debugLevel=debugLevel)
    sfNames = { 
            "GJetsInvForNJets":"GJetsInv", 
            "TTJetsForNJets":"TTJets", 
            "WJetsForNJets":"WJets", 
            "WJetsInvForNJets":"WJetsInv", 
            }
    njetsNames = {
            "GJetsInvForNJets":"NJets_NoPho",
            "TTJetsForNJets":"NJets40",
            "WJetsForNJets":"NJets40",
            "WJetsInvForNJets":"NJets_NoW",
            }
    outfile_name = "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag)
    if args.tightCuts:
        outfile_name = outfile_name.replace(".root", "_TightCuts.root")
    outfile = rt.TFile(outfile_name, "UPDATE")

    for region in regionsOrder:
        analysis = regions[region]
        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        if args.tightCuts:
            outdir += "_TightCuts"
        os.system('mkdir -p '+outdir)
        #set up analysis
        process = sfNames[region]
        (xbins,cols) = analysis.unrollBins
        sfVarsToUse = sfVars[region]
        njetsName = njetsNames[region]
        #use proper set of NJets correction factors.
        #do not correct WJets until WJets scale factors have been derived
        auxSFs = razorWeights.getNJetsSFs(analysis,
                jetName=njetsNames[region])
        if region == "TTJetsForNJets" or region == "WJetsForNJets":
            auxSFs["WJets"] = {}
        if region == "TTJetsForNJets":
            auxSFs["TTJets"] = {}
        if region == "GJetsInvForNJets":
            auxSFs["GJetsInv"] = {}
        if region == "WJetsInvForNJets":
            auxSFs["WJetsInv"] = {}
        if region == "GJetsInvForNJets":
            dataDrivenQCD = True
            razorWeights.getPhotonPuritySFs(auxSFs)
            razorWeights.loadPhotonPurityHists(sfHists, tag, debugLevel)
        else:
            dataDrivenQCD = False
        if args.tightCuts:
            if (region == "GJetsInvForNJets" or region == "WJetsInvForNJets"):
                continue
            analysis.cutsData += " && MR > 500 && Rsq > 0.25 && NJets80 >= 2"
            analysis.cutsMC += " && MR > 500 && Rsq > 0.25 && NJets80 >= 2"
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, 
                sfHists=sfHists, sfVars=sfVarsToUse, 
                printdir=outdir, auxSFs=auxSFs, noFill=args.noFill,
                dataDrivenQCD=dataDrivenQCD, debugLevel=debugLevel )
        #compute scale factors
        sfHistsCopy = sfHists.copy()
        appendScaleFactors( process, hists, sfHistsCopy, lumiData=analysis.lumi, 
                debugLevel=debugLevel, var=njetsName, printdir=outdir )
        if region == "TTJetsForNJets":
            sfHists["NJetsTTJets"] = sfHistsCopy["TTJets"]
            tmpTTJets = sfHists["NJetsTTJets"].Clone()
        elif region == "WJetsForNJets":
            sfHists["NJetsWJets"] = sfHistsCopy["WJets"]
        if not args.noSave:
            #export histograms
            macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                    outDir=outdir, debugLevel=debugLevel )
            #write out scale factors
            print "Writing scale factor histogram",sfHistsCopy[process].GetName(),"to file"
            outfile.cd()
            sfHistsCopy[process].Write( sfHistsCopy[process].GetName() )

    outfile.Close()
