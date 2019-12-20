import os
from sys import exit
import argparse
import copy
import ROOT as rt

#local imports
from macro import macro
from macro.razorAnalysis import Analysis, razorSignalDirs, signalConfig, controlConfig
from macro.razorMacros import unrollAndStitch, unrollAndStitchFromFiles, getMaxBtags
from RunCombine import exec_me
from SMSTemplates import makeSMSTemplates, signalShapeUncerts, SMSOpts
from DustinTuples2DataCard import uncorrelate, uncorrelateSFs, uncorrelateSFs1D, writeDataCard_th1
from framework import Config
from SignalRegionMacro import adjustForFineGrainedMCPred, getSubprocs
from SignalRegionPlotMacro import adjustShapes
import BTagClosureTestMacro as bclosure
import CheckSignalContamination as contam

BACKGROUND_DIR = "/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorMADD2016_08Apr2018"

def getModelName(model, mass1, mass2):
    return "SMS-%s_%d_%d"%(model, mass1, mass2)

def getBranchingFracsFromModelName(model):
    """Returns (x branching ratio, y branching ratio)"""
    xBR = float(model[model.find('x')+1:
        model.find('y')].replace('p','.'))
    yBR = float(model[model.find('y')+1:].replace(
        'p','.'))
    return xBR, yBR

def getCardName(modelName, box, outDir):
    return outDir+'/RazorInclusiveMADD_%s_%s.txt'%(
                modelName, box)

def consolidateBackgroundHists(hists, origProcs):
    """
    Input: a dictionary of histograms, and a list of sample names
    Output: a dictionary of histograms in which the histograms
        belonging to a common MC process are aggregated together.
        Ex: histograms from all HT bins of W+jets will be combined
    """
    # reverse lookup table to speed up find/replaces
    lookup = {}
    for proc in origProcs:
        for subproc in getSubprocs(proc):
            lookup[subproc] = proc
    newHists = {'data_obs': hists['data_obs']}
    for name, hist in hists.iteritems():
        if '_stat' in name:
            continue
        for subproc, proc in lookup.iteritems():
            # note: assumes that no sample name is a substring of any other
            if subproc in name:
                newName = name.replace(subproc, proc)
                if newName not in newHists:
                    newHists[newName] = hist.Clone(newName)
                else:
                    newHists[newName].Add(hist)
                break
    # stat uncertainties need to be treated separately
    # because they add in quadrature
    for name, hist in hists.iteritems():
        if '_stat' not in name:
            continue
        for subproc, proc in lookup.iteritems():
            if subproc in name:
                newName = name.replace(subproc, proc)
                central = hists[subproc]
                combinedCentral = newHists[proc]
                if newName not in newHists:
                    newHists[newName] = combinedCentral.Clone(newName)
                for bx in range(1, central.GetNbinsX()+1):
                    base = combinedCentral.GetBinContent(bx)
                    oldErr = newHists[newName].GetBinContent(bx) - base
                    curBase = central.GetBinContent(bx)
                    curErr = hist.GetBinContent(bx) - curBase
                    newErr = (oldErr*oldErr + curErr*curErr)**(0.5)
                    if name.endswith('Up'):
                        newHists[newName].SetBinContent(bx, base + newErr)
                    elif name.endswith('Down'):
                        newHists[newName].SetBinContent(bx, 
                                max(0, base - newErr))
                    else:
                        raise ValueError("Histogram name should end with Up or Down")
                break
    return newHists

def addHistSuffix(hists, unc, suffix):
    """
    Renames histograms with name *_<unc>
    to *_<unc><suffix>. 
    """
    newHists = {}
    for name, hist in hists.iteritems():
        newName = name
        if name.endswith(unc+'Up') or name.endswith(unc+'Down'):
            newName = name.replace(unc, unc+suffix)
            hist.SetName(hist.GetName().replace(unc, unc+suffix))
        newHists[newName] = hist
    return newHists


if __name__ == "__main__":
    rt.gROOT.SetBatch()

    # Verbosity
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true",
            help="display detailed output messages")
    parser.add_argument("-d", "--debug", action="store_true",
            help="display excruciatingly detailed output messages")
    # Basic configuration
    parser.add_argument('--tag', default='Razor2016_MoriondRereco')
    parser.add_argument('--box', help="choose a box")
    parser.add_argument('--dir', help="output directory",
            default="SignalRegionPlots", dest='outDir')
    parser.add_argument('--bkg-dir', default=BACKGROUND_DIR, dest='bkgDir',
            help='name of directory containing background histograms')
    # Customization
    parser.add_argument('--fine-grained', dest='fineGrained',
            action='store_true', help='Use un-aggregated MC process list')
    parser.add_argument('--no-limit', dest='noCombine', 
            action='store_true', 
            help='do not call combine, make template histograms only')
    parser.add_argument('--signif', action='store_true', 
            help='compute significance rather than limit')
    parser.add_argument('--no-sys',dest="noSys",default=False,
            action='store_true', help="no systematic templates")
    parser.add_argument('--no-stat',dest="noStat",
            default=False,action='store_true', 
            help="no MC statistics uncertainty")
    parser.add_argument('--no-signal-contam', dest='noSignalContam',
            action='store_true', help='ignore signal contamination')
    parser.add_argument('--fit-sys', dest="addMCVsFit", 
            action='store_true', help="add MC vs fit systematic")
    parser.add_argument('--save-workspace', dest='saveWorkspace', action='store_true',
            help='save combine workspace in output file')
    parser.add_argument('--no-boost-cuts', dest='noBoostCuts',
            action='store_true')
    # Signal model
    parser.add_argument('-m','--model', default="T1bbbb", 
            help="signal model name")
    parser.add_argument('--mGluino',default=-1,type=int, 
            help="mass of gluino")
    parser.add_argument('--mStop',default=-1,type=int, 
            help="mass of stop")
    parser.add_argument('--mLSP',default=-1,type=int, 
            help="mass of LSP")

    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    outDir = args.outDir
    cfg = Config.Config(signalConfig)
    boostCuts = not args.noBoostCuts

    if args.box is None:
        print "Please choose an analysis box with --box"
        exit()
    box = args.box
    boxList = box.split('_') # list of boxes to combine, separated by _

    # make output directory
    os.system('mkdir -p '+outDir)

    for curBox in boxList:
        btagsMax = getMaxBtags(curBox)

        # retrieve binning and other info
        analyses = []
        unrollBins = []
        for nb in range(btagsMax + 1):
            analyses.append(Analysis(curBox, args.tag, nbMin=nb))
            unrollBins.append(analyses[-1].unrollBins)
        lumi = analyses[0].lumi
        sfNames = { "TTJets1L":"TTJets", "TTJets2L":"TTJets", 
                "ZInv":"GJetsInv" }
        sfHists = macro.loadScaleFactorHists(
            "data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(args.tag), 
            processNames=["TTJets1L","TTJets2L","WJets","ZInv"], 
            scaleFactorNames=sfNames, debugLevel=debugLevel)
        sfFileNameBClosure = 'data/ScaleFactors/RazorMADD2015/RazorBTagScaleFactors_%s.root'%(args.tag)
        sfFileBClosure = rt.TFile.Open(sfFileNameBClosure)
        sfHistsForUncorrSFs1D = {}
        sfHistsForUncorrSFs1DMR = {}
        jets = 'MultiJet'
        if curBox in ['DiJet', 'LeptonJet']:
            jets = 'DiJet'
        elif curBox in ['SevenJet', 'LeptonSevenJet']:
            jets = 'SevenJet'
        for name in ['TTJets1L', 'TTJets2L', 'WJets']:
            sfHistsForUncorrSFs1D[name] = sfFileBClosure.Get("Rsq{}0B".format(jets))
            assert(sfHistsForUncorrSFs1D[name])
            sfHistsForUncorrSFs1DMR[name] = sfFileBClosure.Get("MR{}0B".format(jets))
            assert(sfHistsForUncorrSFs1DMR[name])
        sfHistsForUncorrSFs1D['ZInv'] = sfFileBClosure.Get("RsqInv{}0B".format(jets))
        assert(sfHistsForUncorrSFs1D['ZInv'])
        sfHistsForUncorrSFs1DMR['ZInv'] = sfFileBClosure.Get("MRInv{}0B".format(jets))
        assert(sfHistsForUncorrSFs1DMR['ZInv'])

        # assess signal contamination in control regions
        contamHists = None
        if not args.noSignalContam:
            print "Computing signal contamination level in control regions"
            ttContamHist = contam.checkSignalContamination(controlConfig, 
                   outDir=outDir, lumi=lumi, debugLevel=debugLevel,
                   box="TTJetsSingleLeptonForSignal", 
                   model=args.model, mLSP=args.mLSP, 
                   mGluino=args.mGluino, mStop=args.mStop, 
                   mergeBins=True)
            wContamHist = contam.checkSignalContamination(controlConfig,
                   outDir=outDir, lumi=lumi, debugLevel=debugLevel,
                   box="WJetsSingleLeptonForSignal", model=args.model,
                   mLSP=args.mLSP, mGluino=args.mGluino, 
                   mStop=args.mStop, mergeBins=True)
            contamHists = { "TTJets1L":ttContamHist, 
                    "TTJets2L":ttContamHist, "WJets":wContamHist }

        if args.fineGrained:
            origProcs = analyses[0].samples
            origSFHists = sfHists
            analyses[0], sfHists, contamHists, _, _ = adjustForFineGrainedMCPred(
                    analyses[0], sfHists, contamHists, {}, {})
        samples = analyses[0].samples

        # make combined unrolled histograms for background
        print "Retrieving background histograms from files"
        doEmptyBinErrs = args.fineGrained
        backgroundHists = unrollAndStitchFromFiles(curBox, 
                samples=samples, inDir=args.bkgDir,
                outDir=outDir, unrollBins=unrollBins, noSys=args.noSys, 
                addStatUnc=(not args.noStat), doEmptyBinErrs=doEmptyBinErrs,
                addMCVsFit=args.addMCVsFit, debugLevel=debugLevel,
                weightHists=analyses[0].weightHists)

        # get file name for signal input
        signalDir = razorSignalDirs[args.tag]
        isGluinos = ('T1' in args.model or 'T5' in args.model)
        if isGluinos:
            modelName = getModelName(args.model, args.mGluino, args.mLSP)
        else:
            modelName = getModelName(args.model, args.mStop, args.mLSP)
        signalFilename=signalDir+'/'+modelName+'.root'

        # to modify branching ratios in T1ttbb sample
        xBR = yBR = -1
        if 'T1x' in args.model:
            xBR, yBR = getBranchingFracsFromModelName(args.model)
            signalFilename = signalDir+'/'+getModelName('T1ttbb', 
                    args.mGluino, args.mLSP)+'.root'

        # get correct list of uncertainties for this box
        uncerts = copy.copy(signalShapeUncerts)
        if curBox in ['DiJet','MultiJet','SevenJet']:
            uncertsToRemove = ['tightmuoneff','tighteleeff',
                    'muontrig','eletrig','tightmuonfastsim',
                    'tightelefastsim']
        else:
            uncertsToRemove = ['vetomuoneff','vetoeleeff',
                    'vetomuonfastsim','vetoelefastsim']
        for u in uncertsToRemove:
            if u in uncerts:
                uncerts.remove(u)

        # call SMS template maker
        smsOpts = SMSOpts(xBR=xBR, yBR=yBR, doNPVExtrap=True,
                doGenMetVsPFMet=True)
        signalHists = makeSMSTemplates(curBox, signalFilename,
                uncertainties=uncerts, debugLevel=debugLevel,
                tag=args.tag, opts=smsOpts, boostCuts=boostCuts)
    
        # reduced efficiency method -- corrects for signal contamination
        if not args.noSignalContam:
            beforeHist = signalHists['Signal'].Clone()
            beforeHist.SetLineColor(rt.kRed)
            macro.doDeltaBForReducedEfficiencyMethod(backgroundHists, 
                    signalHists, contamHists, sfHists, 
                    unrollBins=unrollBins, debugLevel=debugLevel)

        # combine signal and background dictionaries
        if args.fineGrained:
            backgroundHists = consolidateBackgroundHists(
                    backgroundHists, origProcs)
            samples = origProcs
            sfHists = origSFHists
        hists = backgroundHists.copy()
        hists.update(signalHists)

        # do not correlate closure test uncertainties between
        # boxes with different numbers of jets
        jet_closure_uncs = ['btagcrosscheckrsq', 'btaginvcrosscheck',
                'qcdnorm', 'qcdbtag', 'qcdbtagsys', 'qcdaltfunction',
                'ttcrosscheck', 'vetolepetacrosscheck',
                'vetolepptcrosscheck', 'vetotauetacrosscheck',
                'vetotauptcrosscheck', 'zllcrosscheckmr', 'zllcrosscheckrsq',
                'sfstatttjetsNJetsTTJets', 'sfstatwjetsNJetsWJets',
                'sfstatzinvNJetsInv']
        for unc in jet_closure_uncs:
            hists = addHistSuffix(hists, unc, jets)

        # treat appropriate uncertainties as uncorrelated bin to bin
        toUncorrelate = ['stat'+curBox+sample for sample in samples
                if sample != 'QCD']
        for sys in toUncorrelate:
            if 'stat' in sys:
                if args.noStat: continue
                suppressLevel = 0.1
            else:
                if args.noSys: continue
                suppressLevel = 0.0
            uncorrelate(hists, sys, 
                    suppressLevel=suppressLevel)

        # scale factor uncertainties are correlated according to 
        # the bin they are in
        toUncorrelateSF = ['sfstatttjets','sfstatwjets','sfstatzinv']
        for sys in toUncorrelateSF:
            uncorrelateSFs(hists, sys, sfHists, cfg, 
                    curBox, unrollBins=unrollBins)
        toUncorrelateSF1D = ['btaginvcrosscheck', 'btagcrosscheckrsq']
        for sys in toUncorrelateSF1D:
            uncorrelateSFs1D(hists, sys, sfHistsForUncorrSFs1D, unrollBins)
        toUncorrelateSF1DMR = ['sfstatttjetsMR', 'sfstatwjetsMR', 
                'sfstatzinvMRInv']
        for sys in toUncorrelateSF1DMR:
            uncorrelateSFs1D(hists, sys, sfHistsForUncorrSFs1DMR, unrollBins,
                    useRsq=False)
        # two uncertainties that are only binned in nbtags
        dummyHist = rt.TH1F("tmp", "tmp", 1, 0, 4000)
        dummyHistDict = {'QCD': dummyHist, 'ZInv': dummyHist}
        toUncorrelateSF1DBtag = ['qcdbtag', 'sfstatzinvNBTagsInv']
        for sys in toUncorrelateSF1DBtag:
            uncorrelateSFs1D(hists, sys, dummyHistDict, 
                    unrollBins, xInclusive=True)

        # write histograms to ROOT file
        cardName = getCardName(modelName, curBox, outDir)
        outFileName = cardName.replace('.txt', '.root')
        outFile = rt.TFile(outFileName, 'recreate')
        sortedKeys = sorted(hists.keys())
        for key in sortedKeys:
            print "Writing", key
            hists[key].Write()
        outFile.Close()

        writeDataCard_th1(curBox,cardName,hists,samples)

    if args.noCombine: 
        exit()

    # get combine parameters
    if args.signif:
        combineMethod = 'ProfileLikelihood'
        combineFlags = '--signif -t -1 --toysFreq'
    else:
        combineMethod = 'Asymptotic'
        combineFlags = ''
    if args.saveWorkspace:
        combineFlags += ' --saveWorkspace'

    # run combine 
    if len(boxList) == 1:
        cardName = getCardName(modelName, boxList[0], outDir)
        combineName = 'MADD_'+boxList[0]+'_'+modelName
        exec_me('combine -M %s %s -n %s %s'%(
            combineMethod,cardName,combineName,combineFlags), False)
        exec_me('mv higgsCombine%s.%s.mH120.root %s/'%(
            combineName,combineMethod,outDir),False)
    elif len(boxList) > 1:
        cardNames = []
        combineName = 'MADD_'+('_'.join(boxList))+'_'+modelName
        for curBox in boxList:
            cardName = getCardName(modelName, curBox, outDir)
            cardNames.append(os.path.basename(cardName))
        combinedCardName = ('RazorInclusiveMADD_%s_%s.txt'%(
            modelName, '_'.join(boxList)))
        exec_me('cd '+outDir+'; combineCards.py '+(' '.join(cardNames))
                +' > '+combinedCardName+'; cd ..', False)
        exec_me('combine -M '+combineMethod+' '+outDir+'/'
                +combinedCardName+' -n '+combineName+' '+combineFlags, 
                False)
        exec_me('mv higgsCombine%s.%s.mH120.root %s/'%(
            combineName,combineMethod,outDir),False)
