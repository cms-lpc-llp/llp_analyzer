#!/usr/bin/env python
import sys, os
import ROOT as rt

import macro.razorAnalysis as razor
from macro import macro, razorWeights
from macro.razorAnalysis import Analysis, make_parser
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors
from ntupling.ControlRegionNtuples2016_V3p15 import SAMPLES

def getNormErrFractions():
    """
    Returns a dictionary of relative cross section uncertainties
    for rare processes
    """
    return {
        'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1': 0.37, # ATLAS TOPQ-2015-16 (8 TeV)
        'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1': 0.13, # CMS TOP-16-003
        'ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1': 0.13,
        'ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1': 0.32, # ATLAS TOPQ-2015-16
        'ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1': 0.32,
        'WWTo2L2Nu_13TeV-powheg': 0.1, # ATLAS STDM-2015-20
        'WWTo4Q_13TeV-powheg': 0.1,
        'WWToLNuQQ_13TeV-powheg': 0.1,
        'WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8': 0.07, # ATLAS STDM-2015-19
        'WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8': 0.07,
        'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8': 0.07,
        'WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8': 0.07,
        'ZZTo2L2Nu_13TeV_powheg_pythia8': 0.06, # CMS SMP-16-017
        'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8': 0.06,
        'ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8': 0.06,
        'ZZTo4L_13TeV_powheg_pythia8': 0.06,
        'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8': 0.06,
        'ttZJets_13TeV_madgraphMLM': 0.15, # CMS TOP-17-005
        'ttWJets_13TeV_madgraphMLM': 0.23,
        'WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8': 0.05, # ATLAS https://arxiv.org/pdf/1610.05088.pdf
        'WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8': 0.05, # https://arxiv.org/abs/1307.7403
        'WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8': 0.08, # https://arxiv.org/pdf/1507.03693
        'ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8': 0.12, # https://arxiv.org/pdf/1610.05876
        }

def getSubprocs(region, proc):
    regions = {
            'TTJetsSingleLepton':'1L',
            'WJetsSingleLepton':'1L',
            'WJetsSingleLeptonInv':'1LInv',
            'GJetsInv':'Photon',
            }
    if proc == 'QCD':
        return ['QCD']
    return SAMPLES[regions[region]][proc]

def getProcFilename(proc, analysis):
    """
    Returns the path to the MC sample file for the given process.
    """
    if proc == 'QCD':
        # QCD is data driven and should not be split into subprocesses
        return analysis.filenames[proc]
    dirs = {
            'TTJetsSingleLepton':razor.dir1L2016,
            'WJetsSingleLepton':razor.dir1L2016,
            'WJetsSingleLeptonInv':razor.dir1LInv2016,
            'GJetsInv':razor.dirPhoton2016,
            }
    prefixes = {
            'TTJetsSingleLepton':razor.prefixes2016[
                analysis.tag]['SingleLepton'],
            'WJetsSingleLepton':razor.prefixes2016[
                analysis.tag]['SingleLepton'],
            'WJetsSingleLeptonInv':razor.prefixes2016[
                analysis.tag]['SingleLeptonInv'],
            'GJetsInv':razor.prefixes2016[
                analysis.tag]['Photon'],
            }
    return "{}/{}_{}_1pb_weighted{}.root".format(dirs[analysis.region],
            prefixes[analysis.region], proc, razor.skimstr)

def adjustForFineGrainedMCPred(analysis, plotOpts,
        combProcs=['SingleTop', 'Other']):
    """
    Modifies the analysis to run on each rare background MC sample
    individually, rather than processing aggregated groups of samples.
    analysis: Analysis object to be modified
    plotOpts: plot option dict
    combProcs: list of processes that should be split into their component parts
    """
    newSamples = []
    newFilenames = {'Data':analysis.filenames['Data']}
    newExtraWeightOpts = {}
    for proc in analysis.samples:
        if proc not in combProcs:
            newSamples.append(proc)
            newFilenames[proc] = analysis.filenames[proc]
            if proc in analysis.extraWeightOpts:
                newExtraWeightOpts[proc] = analysis.extraWeightOpts[proc]
        else:
            for subproc in getSubprocs(analysis.region, proc):
                newSamples.append(subproc)
                newFilenames[subproc] = getProcFilename(subproc, analysis)
                if proc in analysis.extraWeightOpts:
                    newExtraWeightOpts[subproc] = analysis.extraWeightOpts[proc]
    if 'combineBackgrounds' in plotOpts:
        combBkgs = plotOpts['combineBackgrounds']
        newCombBkgs = {}
        for proc in combBkgs:
            bkgs = []
            for p in combBkgs[proc]:
                if p not in combProcs:
                    bkgs.append(p)
                else:
                    try:
                        bkgs += getSubprocs(analysis.region, p)
                    except KeyError:
                        # this happens when the process is not used
                        continue
            newCombBkgs[proc] = bkgs
        for proc in analysis.samples:
            alreadyDone = False
            for _, procs in combBkgs.iteritems():
                if proc in procs: 
                    alreadyDone = True
                    break
            if not alreadyDone:
                try:
                    newCombBkgs[proc] = getSubprocs(analysis.region, proc)
                except KeyError: # process has no subprocs
                    newCombBkgs[proc] = [proc]
        plotOpts['combineBackgrounds'] = newCombBkgs
    analysis.samples = newSamples
    analysis.filenames = newFilenames
    analysis.extraWeightOpts = newExtraWeightOpts
    return analysis, plotOpts


if __name__ == "__main__":
    rt.gROOT.SetBatch()

    parser = make_parser()
    parser.add_argument('--delta-phi-cut', help='cut on delta phi variable', action='store_true',
            dest='deltaPhiCut')
    parser.add_argument('--njets80-cut', help='require two jets of 80 GeV', action='store_true',
            dest='njets80Cut')
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    boostCuts = not args.noBoostCuts

    #initialize
    sfHists = {}
    regionsOrder = ["TTJetsSingleLepton", "WJetsSingleLepton", "WJetsSingleLeptonInv"]
    if tag != "Razor2015":
        regionsOrder.insert(0, "GJetsInv")
    regions = {
            "TTJetsSingleLepton":Analysis("TTJetsSingleLepton",tag=tag,nbMin=1,boostCuts=boostCuts),
            "WJetsSingleLepton":Analysis("WJetsSingleLepton",tag=tag,nbMax=0,boostCuts=boostCuts),
            "WJetsSingleLeptonInv":Analysis("WJetsSingleLeptonInv",tag=tag,nbMax=0,boostCuts=boostCuts),
            "GJetsInv":Analysis("GJetsInv",tag=tag,boostCuts=boostCuts)
            }

    for region in regionsOrder:
        analysis = regions[region]
        process = region.replace('SingleLepton','')
        (xbins,cols) = analysis.unrollBins
        plotOpts = { 'comment':False, 
                "SUS15004CR":True,
                "combineSamples":analysis.samplesReduced,
                "combineBackgrounds":{ 
                    "Other":["SingleTop","DYJets","Other"]
                    },
                }

        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        if args.noSave:
            outdir += '_Test'
        if args.deltaPhiCut:
            outdir += '_DPhiCut'
        if args.njets80Cut:
            outdir += '_NJets80Cut'
        os.system('mkdir -p '+outdir)
        #get correct variable names
        sfVars = ("MR","Rsq")
        if region == "GJetsInv": sfVars = ("MR_NoPho","Rsq_NoPho")
        elif "Inv" in region: sfVars = ("MR_NoW", "Rsq_NoW")
        #QCD estimate for photon region
        if region == 'GJetsInv':
            dataDrivenQCD = True
            auxSFs = razorWeights.getPhotonPuritySFs()
            razorWeights.loadPhotonPurityHists(sfHists, tag, debugLevel)
        else:
            dataDrivenQCD = False
            auxSFs = {}
        #optionally cut on dPhi
        if args.deltaPhiCut:
            dPhiVar = "dPhiRazor"
            if region == 'GJetsInv':
                dPhiVar += '_NoPho'
            analysis.cutsData += " && abs(%s) < 2.8"%(dPhiVar)
            analysis.cutsMC += " && abs(%s) < 2.8"%(dPhiVar)
        if args.njets80Cut:
            njets80Var = "NJets80"
            if region == 'GJetsInv':
                njets80Var += '_NoPho'
            analysis.cutsData += " && %s >= 2"%(njets80Var)
            analysis.cutsMC += " && %s >= 2"%(njets80Var)
        analysis, plotOpts = adjustForFineGrainedMCPred(analysis, plotOpts=plotOpts)
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, 
                sfHists=sfHists, sfVars=sfVars, printdir=outdir, debugLevel=debugLevel, 
                auxSFs=auxSFs, noFill=args.noFill, dataDrivenQCD=dataDrivenQCD) 
        #compute scale factors
        normErrFractions = getNormErrFractions()
        normErrFractions['QCD'] = 0.02
        appendScaleFactors(process, hists, sfHists, lumiData=analysis.lumi, th2PolyXBins=xbins, 
                th2PolyCols=cols, debugLevel=debugLevel, var=sfVars, printdir=outdir,
                normErrFraction=normErrFractions)
        #export histograms
        if not args.noSave:
            macro.exportHists(hists, outFileName='controlHistograms'+region+'.root', 
                outDir=outdir, debugLevel=debugLevel)

    #write scale factors
    if not args.noSave:
        outname = "data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag)
        if args.deltaPhiCut:
            outname = outname.replace('.root','_DPhiCut.root')
        if args.njets80Cut:
            outname = outname.replace('.root','_NJets80Cut.root')
        outfile = rt.TFile(outname, "RECREATE")
        for name in sfHists:
            print "Writing scale factor histogram",sfHists[name].GetName(),"to file"
            sfHists[name].Write(sfHists[name].GetName().replace("Poly",""))
        outfile.Close()
