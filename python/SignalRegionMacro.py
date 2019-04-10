import sys,os,copy
import ROOT as rt

from framework import Config
from macro import macro, razorWeights
import macro.razorAnalysis as razor
from macro.razorMacros import makeControlSampleHistsForAnalysis
from ntupling.ControlRegionNtuples2016_V3p15 import SAMPLES
import BTagClosureTestMacro as bclosure

commonShapeErrors = [
        ('singletopnorm',"SingleTop"),
        ('othernorm',"Other"),
        ('qcdnorm','QCD'),
        ('qcdbtag','QCD'),
        ('qcdaltfunction', 'QCD'),
        'btag', 'bmistag', 'pileup', 'facscale', 'renscale', 'facrenscale',
        ('btaginvcrosscheck',['ZInv']),
        ('btagcrosscheckrsq',['TTJets1L','TTJets2L','WJets']),
        ('sfstatzinv',['ZInv']),
        ('sfsyszinv',['ZInv']), 'jes',
        ('ttcrosscheck',['TTJets2L']),
        ('zllcrosscheckmr',['ZInv']),
        ('zllcrosscheckrsq',['ZInv']),
        ('sfstatttjets',['TTJets1L','TTJets2L']),
        ('sfsysttjets',['TTJets1L','TTJets2L']),
        ('sfstatwjets',['WJets']),
        ('sfsyswjets',['WJets'])
        ]
lepShapeErrors = commonShapeErrors+['tightmuoneff','tighteleeff','muontrig','eletrig']
hadShapeErrors = commonShapeErrors+['wtag', 'vetolepptcrosscheck','vetotauptcrosscheck',
        'vetolepetacrosscheck','vetotauetacrosscheck','vetomuoneff','vetoeleeff']
shapes = { 'MultiJet':hadShapeErrors, 'LeptonMultiJet':lepShapeErrors, 
           'DiJet':hadShapeErrors, 'LeptonJet':lepShapeErrors, 
           'SevenJet':hadShapeErrors, 'LeptonSevenJet':lepShapeErrors,
           }

def makeSignalRegionParser():
    parser = razor.make_parser()
    parser.add_argument("--unblind", help="do not blind signal sensitive region", action='store_true')
    parser.add_argument('--no-mc', help="do not process MC, do data and fit only", 
            action='store_true', dest="noMC")
    parser.add_argument('--no-fit', help="do not load fit results, process data and MC only", 
            action='store_true', dest='noFit')
    parser.add_argument('--no-data', help="do not process data, do fit and MC only", 
            action='store_true', dest='noData')
    parser.add_argument('--no-sys', help="no shape unncertainties or cross check systematics", 
            action="store_true", dest='noSys')
    parser.add_argument('--no-qcd', help="do not include QCD prediction", action="store_true", 
            dest='noQCD')
    parser.add_argument('--no-sfs', help="ignore MC scale factors", action="store_true", 
            dest="noSFs")
    parser.add_argument('--box', help="choose a box")
    parser.add_argument('--btags', type=int, help="choose a number of btags")
    parser.add_argument('--b-inclusive', help='do not bin in btags', action='store_true',
            dest='bInclusive')
    parser.add_argument('--qcd-mc', dest='qcdMC', help='make qcd prediction using MC',
            action='store_true')
    parser.add_argument('--nlo-zinv', dest='nloZInv', action='store_true',
            help='use NLO sample for Z->nu nu prediction')
    parser.add_argument('--fine-grained', dest='fineGrained', action='store_true',
            help='Process each MC sample separately')
    parser.add_argument('--sideband', action='store_true',
            help='Process sideband instead of extrapolation region')
    parser.add_argument('--zoom', action='store_true',
            help='Draw restricted range on ratio plot')
    parser.add_argument('--super-region', action='store_true',
            help='Run on aggregated super analysis regions')
    return parser

def getDirSuffix(args):
    dirSuffix = ""
    if not args.unblind:
        dirSuffix += 'Blinded'
    if args.noSFs:
        dirSuffix += 'NoSFs'
    if args.noSys:
        dirSuffix += 'NoSys'
    if args.bInclusive:
        dirSuffix += 'BInclusive'
    if args.qcdMC:
        dirSuffix += 'QCDMC'
    if args.nloZInv:
        dirSuffix += 'NLOZInv'
    if args.fineGrained:
        dirSuffix += 'FineGrained'
    if args.sideband:
        dirSuffix += 'Sideband'
    if args.noBoostCuts:
        dirSuffix += 'NoBoostCuts'
    if args.super_region:
        dirSuffix += 'Super'

    return dirSuffix

def getPlotOpts(args, analysis):
    plotOpts = {
            'SUS15004':True,
            "combineBackgrounds":{ 
                "Other":["SingleTop","DYJets","Other"], 
                "TTJets":["TTJets1L","TTJets2L","TTJets"] },
            "combineSamples":analysis.samplesReduced,
            'sideband':True,
            }
    if args.noFit:
        del plotOpts['sideband']
    if args.zoom:
        plotOpts['zoom'] = True
    return plotOpts

def getBoxesAndBtags(args):
    """Returns list of boxes to process
        and list of b-tag categories"""
    boxlist = ["SevenJet", "LeptonSevenJet", "MultiJet", "LeptonMultiJet", "DiJet"]
    if args.box is not None:
        boxlist = [args.box]
    if args.btags is not None:
        btaglist = [args.btags]
    elif args.bInclusive:
        btaglist = [0]
    else:
        btaglist = [0,1,2,3]
    return boxlist, btaglist

def getNBMinMax(box, btags, args=None):
    """Returns appropriate min and max number of b-tags"""
    if args is not None and args.bInclusive:
        return 0, -1

    nbMin = btags

    btagCap = 3
    if box in ["DiJet", "MuJet", "EleJet", "LeptonJet"]:
        btagCap = 2
    
    if btags > btagCap: 
        raise ValueError("Too many b-tags for box {}".format(box))
    elif btags == btagCap:
        nbMax = -1 #no upper limit
    else:
        nbMax = btags #exclusive

    return nbMin, nbMax

def loadAllScaleFactorHists(tag, args, processNames, debugLevel=0):
    #scale factor file names
    sfdir = "data/ScaleFactors/RazorMADD2015/"
    sfFile = sfdir+'/RazorScaleFactors_%s.root'%(tag)
    sfFile_nJets = sfdir+'/RazorNJetsScaleFactors_%s.root'%(tag)
    vetolepFile = sfdir+'/RazorVetoLeptonClosureTests_%s.root'%(tag)
    ttFile = sfdir+'/RazorTTJetsDileptonCrossCheck_%s.root'%(tag)
    dyFile = sfdir+'/RazorDYJetsDileptonInvCrossCheck_%s.root'%(tag)
    dybFile = sfdir+'/RazorDYJetsDileptonInvBCheck_%s.root'%(tag)
    btagFile = sfdir+'/RazorBTagScaleFactors_%s.root'%(tag)
    dynloFile = sfdir+'/RazorDYJetsDileptonInvNLOCrossCheck_%s.root'%(tag)

    #get MR-Rsq scale factor histograms
    sfNames={
            "ZInv":"GJetsInv",
            "TTJets1L":"TTJets",
            "TTJets2L":"TTJets",
            }
    sfHists = macro.loadScaleFactorHists(sfFilename=sfFile, processNames=processNames, 
            scaleFactorNames=sfNames, debugLevel=debugLevel)
    bclosure.loadScaleFactors(sfHists, tag=tag)
    bclosure.loadScaleFactors(sfHists, tag=tag, gjets=True)
    #reopen the file and grab the ZNuNu up/down histograms
    #down scale factors are (gjets - (wjets-gjets))
    sfTFile = rt.TFile.Open(sfFile)
    sfHists['ZInvUp'] = sfTFile.Get('WJetsInvScaleFactors')
    sfHists['ZInvDown'] = sfTFile.Get('GJetsInvScaleFactors_Down') 
    #get njets scale factor histogram
    sfNJetsFile = rt.TFile.Open(sfFile_nJets)
    sfHists['NJetsTTJets'] = sfNJetsFile.Get("TTJetsScaleFactors")
    sfHists['NJetsWJets'] = sfNJetsFile.Get("WJetsScaleFactors")
    sfHists['NJetsInv'] = sfNJetsFile.Get("GJetsInvScaleFactors")
    #nb double ratio scale factors
    dybTFile = rt.TFile.Open(dybFile)
    sfHists['NBTagsInv'] = dybTFile.Get("DYJetsDileptonInvBCheckScaleFactors")
    #get veto lepton/tau, DYJets, and TTBar Dilepton cross check scale factor histograms
    #and b-tag closure results
    vlFile = rt.TFile.Open(vetolepFile)
    ttTFile = rt.TFile.Open(ttFile)
    dyTFile = rt.TFile.Open(dyFile)
    btagTFile = rt.TFile.Open(btagFile)
    dynloTFile = rt.TFile.Open(dynloFile)
    if args.nloZInv:
        sfHists['ZInv'] = dynloTFile.Get("WJetsSingleLeptonInvNLOScaleFactors")
        sfHists['NJetsInv'] = dynloTFile.Get("WJetsSingleLeptonInvNLOForNJetsScaleFactors")
    for jtype in ['DiJet','MultiJet','SevenJet']:
        sfHists['TTJetsDilepton'+jtype+'Up'] = ttTFile.Get("TTJetsDilepton"+jtype+"ScaleFactors")
        sfHists['TTJetsDilepton'+jtype+'Down'] = macro.invertHistogram(
                sfHists['TTJetsDilepton'+jtype+'Up'])
        sfHists['DYJetsInvMR'+jtype+'Up'] = dyTFile.Get('DYJetsDileptonInv'+jtype+'MRScaleFactors')
        sfHists['DYJetsInvRsq'+jtype+'Up'] = dyTFile.Get('DYJetsDileptonInv'+jtype+'RsqScaleFactors')
        if args.nloZInv:
            sfHists['DYJetsInvMR'+jtype+'Up'] = dynloTFile.Get("DYJetsDileptonInvNLO"+jtype+"MRScaleFactors")
            sfHists['DYJetsInvRsq'+jtype+'Up'] = dynloTFile.Get("DYJetsDileptonInvNLO"+jtype+"RsqScaleFactors")
        sfHists['DYJetsInvMR'+jtype+'Down'] = macro.invertHistogram(
                sfHists['DYJetsInvMR'+jtype+'Up'])
        sfHists['DYJetsInvRsq'+jtype+'Down'] = macro.invertHistogram(
                sfHists['DYJetsInvRsq'+jtype+'Up'])
        for ltype in ['VetoLepton','VetoTau']:
            name = jtype+'For'+ltype
            sfHists[ltype+jtype+'PtUp'] = vlFile.Get(name+'ScaleFactors')
            sfHists[ltype+jtype+'PtDown'] = macro.invertHistogram(sfHists[ltype+jtype+'PtUp'])
            sfHists[ltype+jtype+'EtaUp'] = vlFile.Get(name+'PtCorrScaleFactors')
            sfHists[ltype+jtype+'EtaDown'] = macro.invertHistogram(sfHists[ltype+jtype+'EtaUp'])
        for b in range(4):
            if jtype == 'DiJet' and b > 2: continue
            bs = str(b)
            sfHists['Rsq'+jtype+bs+'BUp'] = btagTFile.Get(
                    'Rsq'+jtype+bs+'B')
            sfHists['Rsq'+jtype+bs+'BDown'] = macro.invertHistogram(sfHists['Rsq'+jtype+bs+'BUp'])
            if b > 2: continue
            sfHists['ZInv'+jtype+bs+'BUp'] = btagTFile.Get(
                    'RsqInv'+jtype+bs+'B')
            sfHists['ZInv'+jtype+bs+'BDown'] = macro.invertHistogram(
                    sfHists['ZInv'+jtype+bs+'BUp'])

    #check that everything came out correctly
    for h,hist in sfHists.iteritems():
        if debugLevel > 0:
            print "Checking scale factor histogram:",h
        #assert hist
        #hist.SetDirectory(0)
        if hist:
            hist.SetDirectory(0)
        elif 'MultiJet' not in h and 'DiJet' not in h:
            raise ValueError("Hist {} is null".format(h))
    return sfHists

def adjustCuts(analysis, boxName, sideband=False):
    if sideband:
        print "Adjust baseline cuts to select only sideband"
        if boxName in ['DiJet', 'MultiJet', 'SevenJet']:
            minMR = 500
            maxMR = 650
            minRsq = 0.25
            maxRsq = 0.30
        else:
            minMR = 400
            maxMR = 550
            minRsq = 0.15
            maxRsq = 0.20
        # quick hack, not very stable
        replaceStr = 'MR > {} && MR < 4000 && Rsq > {} && Rsq < 1.5'.format(minMR, minRsq)
        newStr = 'MR > {} && Rsq > {} && (MR < {} || Rsq < {})'.format(
                    minMR, minRsq, maxMR, maxRsq)
        analysis.cutsData = analysis.cutsData.replace(replaceStr, newStr)
        analysis.cutsMC = analysis.cutsMC.replace(replaceStr, newStr)
        analysis.unrollBins = (
                razor.xbinsSignalSideband[boxName]["{}B".format(analysis.nbMin)],
                razor.colsSignalSideband[boxName]["{}B".format(analysis.nbMin)])
        # need binning that extends down into sideband region
        cfg = Config.Config("config/run2_2017_03_13_SeparateBtagFits_forToys.config")
        analysis.binning['MR'] = cfg.getBinning(boxName)[0]
        analysis.binning['Rsq'] = cfg.getBinning(boxName)[1]

def applyAnalysisOptions(analysis, args, boxName=None):
    """
    Changes Analysis object attributes according to specified options"""
    print "Applying options to {}".format(boxName)
    if args.noQCD and 'QCD' in analysis.samples:
        analysis.samples.remove('QCD')
        analysis.samplesReduced.remove('QCD')
    elif args.qcdMC:
        analysis.filenames['QCD'] = "Backgrounds/Signal/FullRazorInclusive_%s_QCD_1pb_weighted.root"%tag
    if args.nloZInv:
        analysis.filenames['ZInv'] = analysis.filenames['ZInv'].replace('ZInv_', 'ZInvPtBinned_')
    if args.noMC: analysis.samples = []
    if analysis.samples is None or len(analysis.samples) == 0:
        analysis.filenames = {"Data":analysis.filenames["Data"]}
    if args.noData: 
        del analysis.filenames['Data']
    if boxName is not None:
        adjustCuts(analysis, boxName, sideband=args.sideband)
        if args.super_region:
            analysis.binning = razor.razorBinning['Super'+boxName]
            analysis.unrollBins = (razor.xbinsSignal['Super'+boxName],
                    razor.colsSignal['Super'+boxName])

def getScaleFactorHistsForBox(sfHists, boxName, btags):
    sfHistsToUse = sfHists.copy()
    if 'SevenJet' in boxName:
        jtype = 'SevenJet'
    elif 'MultiJet' in boxName:
        jtype = 'MultiJet'
    else:
        jtype = 'DiJet'
    #veto lepton/tau
    for ltype in ['VetoLepton','VetoTau']:
        for pteta in ['Pt','Eta']:
            for updown in ['Up','Down']:
                sfHistsToUse[ltype+pteta+updown] = sfHistsToUse[ltype+jtype+pteta+updown]
    ##ttbar dilepton, dyjets dilepton, and zinv b-tag
    for name in ['TTJetsDilepton','DYJetsInvMR', 'DYJetsInvRsq']:
        for updown in ['Up','Down']:
            sfHistsToUse[name+updown] = sfHistsToUse[name+jtype+updown]
    #b-tag closure tests
    for updown in ['BUp','BDown']:
        sfHistsToUse['Rsq'+updown] = sfHistsToUse[
                'Rsq'+jtype+str(btags)+updown]
        sfHistsToUse['ZInv'+updown] = sfHistsToUse[
                'ZInv'+jtype+str(min(btags, 2))+updown]
    return sfHistsToUse

def removeSFShapes(shapes):
    toRemove = ['btaginvcrosscheck','btagcrosscheckrsq',
            'btagcrosscheckmr','sfsyszinv','ttcrosscheck',
            'zllcrosscheckmr','zllcrosscheckrsq','sfsysttjets','sfsyswjets',
            'vetolepptcrosscheck','vetotauptcrosscheck',
            'vetolepetacrosscheck','vetotauetacrosscheck']
    #remove scale factor cross check uncertainties
    shapes = [s for s in shapes if s not in toRemove]
    #this removes scale factor uncertainties that are listed as tuples
    shapes = [s for s in shapes if not (hasattr(s, '__getitem__') and s[0] in toRemove)] 
    return shapes

def getSubprocs(proc):
    """
    Returns a list of all subprocesses associated with 
    the given physics process.
    """
    if proc == 'QCD':
        # QCD is data driven and should not be split into processes
        return ['QCD']
    return SAMPLES['Signal'][proc]

def makeFineGrainedShapeErrors(shapes):
    """
    Breaks each MC process into its component processes
    and assigns appropriate uncertainties to each component.
    Input: shape uncertainty dictionary of the form
    { box:[unc1, unc2, ...], ... }
    where each uncertainty is either 1) a string,
    2) a pair of strings (unc name, process name), or
    3) a pair (unc name, [process1, process2, ...])
    Output: a shape uncertainty dictionary in which each
    process name has been replaced by a list of its component processes.
    """
    outShapes = []
    for unc in shapes:
        try:
            uncName, uncProcs = unc
            outProcs = []
            if isinstance(uncProcs, basestring):
                # there is only one process for this uncertainty
                outProcs += getSubprocs(uncProcs)
            else:
                # there are multiple top-level processes
                for proc in uncProcs:
                    outProcs += getSubprocs(proc)
            outShapes.append((uncName, outProcs))
        except ValueError:
            # this uncertainty applies to all processes
            outShapes.append(unc)
    return outShapes

def getProcFilename(proc, analysis):
    """
    Returns the path to the MC sample file for the given process.
    """
    if proc == 'QCD':
        # QCD is data driven and should not be split into subprocesses
        return analysis.filenames[proc]
    return "{}/{}_{}_1pb_weighted{}.root".format(razor.dirSignal2016,
            razor.prefixes2016[analysis.tag]['Signal'], proc, razor.skimstr)

def adjustForFineGrainedMCPred(analysis, sfHists, auxSFs, shapes, plotOpts):
    """
    Modifies the analysis to run on each MC sample individually
    rather than processing aggregated groups of samples.
    analysis: Analysis object
    sfHists: dictionary of scale factor histograms
    auxSFs: dictionary of auxiliary scale factor specifications
    shapes: dictionary of uncertainty options
    plotOpts: dictionary of plotting options
    """
    newSamples = []
    newFilenames = {'Data':analysis.filenames['Data']}
    newExtraWeightOpts = {}
    newSFHists = {}
    newAuxSFs = {}
    for proc in analysis.samples:
        for subproc in getSubprocs(proc):
            newSamples.append(subproc)
            newFilenames[subproc] = getProcFilename(subproc, analysis)
            if proc in auxSFs:
                newAuxSFs[subproc] = auxSFs[proc]
            if proc in analysis.extraWeightOpts:
                newExtraWeightOpts[subproc] = analysis.extraWeightOpts[proc]
    for name, hist in sfHists.iteritems():
        isProcessSFHist = False
        for proc in analysis.samples:
            if name == proc:
                isProcessSFHist = True
                for subproc in getSubprocs(proc):
                    newSFHists[subproc] = hist
            elif name == proc+'Up':
                isProcessSFHist = True
                for subproc in getSubprocs(proc):
                    newSFHists[subproc+'Up'] = hist
            elif name == proc+'Down':
                isProcessSFHist = True
                for subproc in getSubprocs(proc):
                    newSFHists[subproc+'Down'] = hist
        if not isProcessSFHist:
            newSFHists[name] = hist
    if 'combineBackgrounds' in plotOpts:
        combBkgs = plotOpts['combineBackgrounds']
        newCombBkgs = {}
        for proc in combBkgs:
            bkgs = []
            for p in combBkgs[proc]:
                try:
                    bkgs += getSubprocs(p)
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
                newCombBkgs[proc] = getSubprocs(proc)
        plotOpts['combineBackgrounds'] = newCombBkgs
    analysis.samples = newSamples
    analysis.filenames = newFilenames
    analysis.extraWeightOpts = newExtraWeightOpts
    newShapes = makeFineGrainedShapeErrors(shapes)
    return analysis, newSFHists, newAuxSFs, newShapes, plotOpts

def getAllAuxSFs(analysis):
    auxSFs = razorWeights.getNJetsSFs(analysis, jetName='nSelectedJets')
    auxSFs = razorWeights.addBTagSFs(analysis, auxSFs)
    auxSFs = razorWeights.addBTagSFs(analysis, auxSFs, gjets=True)
    auxSFs = razorWeights.addBTagDoubleRatioSFs(analysis, auxSFs)
    return auxSFs


if __name__ == "__main__":
    rt.gROOT.SetBatch()
    parser = makeSignalRegionParser()
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    boostCuts = not args.noBoostCuts
    dirSuffix = getDirSuffix(args)
    if not args.noFit:
        print "Disabling fit option until needed again"
        args.noFit = True
    if args.noFit: 
        toysToUse = {}
    else:
        toysToUse = razor.razorFitFiles[tag]
    boxesToUse, btaglist = getBoxesAndBtags(args)

    regionsOrder = []
    regions = {}
    for box in boxesToUse:
        for btags in btaglist:
            if box in ['DiJet', 'LeptonJet'] and btags > 2: continue
            nbMin, nbMax = getNBMinMax(box, btags, args)
            extBox = '%s%dB%s'%(box,btags,dirSuffix)
            regionsOrder.append(extBox)
            regions[extBox] = razor.Analysis(box, tag=tag, boostCuts=boostCuts,
                    nbMin=nbMin, nbMax=nbMax)
    processNames = regions.itervalues().next().samples
    sfHists = loadAllScaleFactorHists(tag, args, 
            processNames, debugLevel=debugLevel)

    for region in regionsOrder:
        analysis = regions[region]
        plotOpts = getPlotOpts(args, analysis)
        boxName = region.replace(dirSuffix,'')[:-2]
        btags = int(region.replace(dirSuffix,'')[-2])
        shapesToUse = copy.copy(shapes[boxName])
        auxSFs = getAllAuxSFs(analysis)
        dataDrivenQCD = True
        blindBins = [(x,y) for x in range(1,len(analysis.binning["MR"])+1) 
                for y in range(1,len(analysis.binning["Rsq"])+1)]
        bclosure.adjustForRegion(analysis, sfHists, auxSFs)
        bclosure.adjustForRegion(analysis, sfHists, auxSFs, gjets=True)
        sfHistsToUse = getScaleFactorHistsForBox(sfHists, boxName, btags)
        auxSFsToUse = auxSFs.copy()

        # modify analysis according to options
        applyAnalysisOptions(analysis, args, boxName)
        if args.unblind: 
            blindBins = None
        if args.qcdMC:
            dataDrivenQCD = False
        if args.noSys:
            shapesToUse = []
        if args.noSFs:
            print "Ignoring all Data/MC scale factors"
            sfHistsToUse = {}
            auxSFsToUse = {}
            shapesToUse = removeSFShapes(shapesToUse)
        if args.fineGrained:
            (analysis, sfHistsToUse, auxSFsToUse, shapesToUse, 
                    plotOpts) = adjustForFineGrainedMCPred(
                            analysis, sfHistsToUse, auxSFsToUse, shapesToUse, plotOpts)

        print "\nBox:",region,"("+boxName,str(btags),"B-tag)"

        #make output directory
        outdir = "Plots/"+tag+"/"+region
        os.system('mkdir -p '+outdir)

        #run analysis
        hists = makeControlSampleHistsForAnalysis(analysis,
                sfHists=sfHistsToUse, treeName="RazorInclusive", 
                shapeErrors=shapesToUse, fitToyFiles=toysToUse, boxName=boxName, 
                blindBins=blindBins, btags=btags, debugLevel=debugLevel, 
                auxSFs=auxSFsToUse, dataDrivenQCD=dataDrivenQCD, printdir=outdir, 
                plotOpts=plotOpts, noFill=args.noFill, exportShapeErrs=True, 
                propagateScaleFactorErrs=False, makePlots=False)
        if not args.noSave:
            #export histograms
            macro.exportHists(hists, outFileName='razorHistograms'+region+'.root', 
                    outDir=outdir, debugLevel=debugLevel)
