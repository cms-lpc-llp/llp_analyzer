import sys
import argparse
import copy
import ROOT as rt
from array import array

#local imports
from DustinTuples2DataCard import uncorrelate, getTheoryCrossSectionAndError
from macro.razorAnalysis import Analysis
from macro.razorMacros import makeControlSampleHistsForAnalysis, getMaxBtags, unrollAndStitch, makeRazor2DTable
from macro.razorWeights import getNPVHist, getNISRScaleFactor

signalShapeUncerts = ['tightmuoneff','tighteleeff','vetomuoneff',
        'vetoeleeff','jes','muontrig','eletrig','btag','bmistag',
        'tightmuonfastsim','tightelefastsim','vetomuonfastsim',
        'vetoelefastsim','btagfastsim','facscale','renscale',
        'facrenscale','isr']

class SMSOpts(object):
    """
    Structure holding SMS processing options.
    Attributes:
    -xBR, yBR (floats): branching ratios for event reweighting
    -doNPVExtrap (bool): whether to perform NPV reweighting procedure
    """
    defaults = {
            "xBR":-1,
            "yBR":-1,
            "doNPVExtrap":False,
            "doGenMetVsPFMet":False,
            }
    def __init__(self, **kwargs):
        for param, default in self.defaults.iteritems():
            setattr(self, param, kwargs.get(param, default))

def getModelInfoFromFilename(f, xBR=-1, yBR=-1):
    """f: filename
       xBR, yBR: branching fractions
       returns: model_name, mass1, mass2"""
    modelString = '_'.join(f.split('/')[-1].split(
        '.root')[0].split('_')[:-2])
    if xBR>-1 and yBR>-1:
        modelString = modelString.replace('T1ttbb',
                ('T1x%.2fy%.2f'%(xBR,yBR)).replace('.','p'))
    modelName = modelString.split('-')[-1]
    massPoint = '_'.join(f.split('/')[-1].split(
        '.root')[0].split('_')[1:])
    m1 = massPoint.split('_')[-2]
    m2 = massPoint.split('_')[-1]
    return modelName, m1, m2

def getNPVLowHighHist(analysis):
    """From Analysis object, gets two-bin NPV histogram"""
    signalFile = rt.TFile.Open(analysis.filenames['Signal'])
    npvLowHighHist = signalFile.Get('NPV')
    assert npvLowHighHist
    npvLowHighHist.SetDirectory(0)
    return npvLowHighHist

def doNPVExtrapolation(hists, npvHist, npvLowHighHist, scale):
    """
    Implements the NPV extrapolation procedure described here:
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsMoriond17
    hists: dictionary of histograms, should include 
            the keys 'Signal_npvextrapDown/Up'
    npvHist: 1D histogram with the full NPV distribution in data
    npvLowHighHist: two-bin 1D histogram storing the number of 
            signal events with less/more than 20 primary vertices
    scale: event weight, used to assign uncertainties on empty bins
    """
    print "Performing NPV extrapolation procedure for signal MC"
    histLowNPV = hists['Signal_npvextrapDown']
    histHighNPV = hists['Signal_npvextrapUp']

    x = array('d', [14.68, 24.26]) # bin centers of npvLowHighHist
    ex = array('d', [0,0])
    for ibin in range(1, histLowNPV.GetNbinsX()+1):
        lowYield = histLowNPV.GetBinContent(ibin)
        lowError = histLowNPV.GetBinError(ibin)
        highYield = histHighNPV.GetBinContent(ibin)
        highError = histHighNPV.GetBinError(ibin)
        # Bins with zero yield still need to have nonzero uncertainty,
        # otherwise chi squared cannot be computed. Assign an error
        # based on the appropriate Poisson confidence interval.
        if not lowError:
            lowError = 1.83 * scale
        if not highError:
            highError = 1.83 * scale

        npvNorm = (npvLowHighHist.GetBinContent(1) 
                + npvLowHighHist.GetBinContent(2))
        lowNPVFrac = npvLowHighHist.GetBinContent(1)/npvNorm
        highNPVFrac = npvLowHighHist.GetBinContent(2)/npvNorm

        # Perform linear fit to 2-bin NPV distribution
        y = array('d', [lowYield/lowNPVFrac, highYield/highNPVFrac])
        ey = array('d', [lowError/lowNPVFrac, highError/highNPVFrac])
        graph = rt.TGraphErrors(2, x, y, ex, ey)
        fitResult = graph.Fit("pol1", "SMF")
        p0 = fitResult.Parameter(0)
        p1 = fitResult.Parameter(1)
        fitCov = fitResult.GetCovarianceMatrix()
     
        # Convolve linear fit with NPV distribution in data
        averageYield = 0
        averageYieldErr = 0
        for npvBin in range(1, npvHist.GetNbinsX()+1):
            npv = npvHist.GetXaxis().GetBinCenter(npvBin)
            npvWeight = npvHist.GetBinContent(npvBin)

            fitPred = p0 + npv * p1
            fitError = ( npv*npv * fitCov(1,1) + 
                    fitCov(0,0) + 2*npv * fitCov(0,1) )**(0.5)

            averageYield += npvWeight * fitPred
            averageYieldErr += npvWeight * fitError

        averageYield = max(0., averageYield)
        # Scale all signal histograms based on these weights
        nominal = hists['Signal'].GetBinContent(ibin)
        poisson = hists['Signal'].GetBinError(ibin) # nominal error
        if nominal > 0:
            weightedOverNominal = (averageYield / nominal)
            print "In bin %d, weighted yield is %.3f of nominal"%(
                    ibin,weightedOverNominal)
            for hName, hist in hists.iteritems():
                if 'npvextrap' in hName:
                    continue
                hist.SetBinContent(ibin, 
                        hist.GetBinContent(ibin) * weightedOverNominal)

        # Put the up/down errors from this procedure 
        # into the npvextrap histograms.
        # If the uncertainty obtained is larger than
        # 1.6 times the statistical uncertainty from
        # Poisson statistics, it is capped at that level.
        # This prevents instability in limit setting
        # due to large uncertainties from this procedure.
        errThreshold = 1.6 * poisson
        if averageYieldErr > errThreshold:
            print "Reducing uncertainty on bin {} ({:.2f} events) from {:.2f} to {:.2f}".format(
                    ibin, nominal, averageYieldErr, errThreshold)
            averageYieldErr = errThreshold
        hists['Signal_npvextrapUp'].SetBinContent(ibin,
                averageYield + averageYieldErr)
        hists['Signal_npvextrapDown'].SetBinContent(ibin,
                max(0.,averageYield - averageYieldErr))

    return hists

def doGenMetVsPFMetSystematic(hists):
    """Implements the procedure described here:
       https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsMoriond17.
       The central yield in each bin is the average of the yields
       from PF MET and gen-MET.  All up/down shape histograms are
       modified accordingly. The uncertainty on the procedure is
       taken to be half of the difference between the acceptances."""
    print "Correcting yields using gen-MET vs PF MET comparison"
    # The 'down' histogram may not exist yet
    hists['Signal_genmetvspfmetDown'] = (
            hists['Signal_genmetvspfmetUp'].Clone(
                'Signal_genmetvspfmetDown'))
    for ibin in range(1, hists['Signal'].GetNbinsX()+1):
        pfmetYield = hists['Signal'].GetBinContent(ibin)
        genmetYield = hists['Signal_genmetvspfmetUp'].GetBinContent(ibin)
        central = (pfmetYield+genmetYield)/2.0
        unc = (pfmetYield-genmetYield)/2.0

        hists['Signal_genmetvspfmetUp'].SetBinContent(
                ibin, central + unc)
        hists['Signal_genmetvspfmetDown'].SetBinContent(
                ibin, central - unc)

        if pfmetYield <= 0:
            continue
        print "Multiplying bin content by",central/pfmetYield
        for hName,hist in hists.iteritems():
            if 'genmetvspfmet' in hName:
                continue
            hist.SetBinContent(ibin, hist.GetBinContent(ibin)
                    * central/pfmetYield)

def getGlobalNISRScaleFactor(f):
    """
    Gets the scale factor needed to avoid
    the NISR weights changing the overall
    signal cross section
    """
    nisr = f.Get("NISRJets") 
    numUnweighted = 0
    numWeighted = 0
    for nisrJets in range(nisr.GetNbinsX()):
        nisrEvents = nisr.GetBinContent(nisrJets+1)
        nisrWeight = getNISRScaleFactor(nisrJets)
        numUnweighted += nisrEvents
        numWeighted += nisrEvents * nisrWeight
    return numUnweighted / numWeighted

def makeSMSTemplates(box, inFile, uncertainties=[], debugLevel=0,
        tag="Razor2016_MoriondRereco", opts=None, boostCuts=True,
        doUncorrelate=True):
    """Returns a dictionary of histograms representing predicted
        yields for the indicated signal sample.
        'opts' should be an SMSOpts instance containing desired
        values for configurable parameters."""
    print '\nInput file is %s' % inFile
    uncerts = copy.copy(uncertainties)

    # process configuration options
    if opts is None:
        opts = SMSOpts()
    xBR = opts.xBR
    yBR = opts.yBR
    doMCStat = True
    if opts.doNPVExtrap:
        uncerts.append('npvextrap')
        doMCStat = False
    if opts.doGenMetVsPFMet:
        uncerts.append('genmetvspfmet')

    # get mass point information
    model, m1, m2 = getModelInfoFromFilename(inFile, xBR, yBR)
    isGluinos = ("T1" in model or "T5" in model)
    if isGluinos:
        thyXsec, _ = getTheoryCrossSectionAndError(
                mGluino=m1)
    else:
        thyXsec, _ = getTheoryCrossSectionAndError(
                mStop=m1)
    print "\n--- %s %s %s ---\n"%(model, m1, m2)
    print "Theory cross section: %.3f pb"%thyXsec

    minBtags = 0
    maxBtags = getMaxBtags(box)
    # special case: 1L control regions
    if box == "TTJetsSingleLeptonForSignal":
        minBtags = maxBtags = 1
    elif box == "WJetsSingleLeptonForSignal":
        minBtags = maxBtags = 0
    hists = []
    unrollBins = []
    for nb in range(minBtags, maxBtags+1):
        # get analysis info
        nbMax = nb
        if nb == maxBtags:
            nbMax = -1
        # special case: 1L control regions
        if box == "WJetsSingleLeptonForSignal":
            nbMax = 0
        analysis = Analysis(box, tag, nbMin=nb, nbMax=nbMax,
                boostCuts=boostCuts)
        unrollBins.append(analysis.unrollBins)

        # modify for signal sample
        analysis.filenames = { 'Signal':inFile }
        analysis.samples = ['Signal']
        analysis.samplesReduced = analysis.samples
            
        # scale according to cross section
        f = rt.TFile.Open(inFile)
        nEvents = f.Get('NEvents').Integral()
        globalScaleFactor = thyXsec/nEvents 
        nisrFactor = getGlobalNISRScaleFactor(f)
        globalScaleFactor *= nisrFactor
        print "Number of events: %d"%nEvents
        print "Integrated luminosity: %d /pb"%analysis.lumi
        print "Overall NISR scale factor: {}".format(nisrFactor)
        print "Overall scale factor: %.3f"%(
                analysis.lumi * globalScaleFactor)

        # fill the histograms
        hists.append(makeControlSampleHistsForAnalysis(analysis, 
                treeName="RazorInclusive", shapeErrors=uncerts,
                boxName=box, btags=nb, makePlots=False, 
                exportShapeErrs=True, propagateScaleFactorErrs=False,
                lumiMC=1./globalScaleFactor, debugLevel=debugLevel))

        makeRazor2DTable(pred=None, 
                boxName='{}-{}-{}-{}-'.format(box, model, m1, m2), 
                mcNames=['{} {} {}'.format(model, m1, m2)], 
                mcHists=[hists[-1]['Signal'][('MR','Rsq')]], 
                btags=nb, unrollBins=unrollBins[-1], listAllMC=True)
            
    signalHists = unrollAndStitch(hists, box, samples=analysis.samples,
            debugLevel=debugLevel, unrollBins=unrollBins, 
            addStatUnc=doMCStat)

    # apply pileup weight extrapolation procedure
    if 'Signal_npvextrapUp' in signalHists:
        npvLowHighHist = getNPVLowHighHist(analysis)
        npvHist = getNPVHist(tag)
        doNPVExtrapolation(signalHists, npvHist, npvLowHighHist,
                scale=analysis.lumi*globalScaleFactor)
        # need this extrapolation to be uncorrelated across bins & boxes
        for sysName in ['Signal_npvextrapUp', 'Signal_npvextrapDown']:
            signalHists[sysName.replace(
                'npvextrap','npvextrap'+box)] = signalHists[sysName]
            del signalHists[sysName]
        if doUncorrelate:
            uncorrelate(signalHists, 'npvextrap'+box, suppressLevel=0.1)

    # apply gen-MET vs PF MET systematic
    if 'Signal_genmetvspfmetUp' in signalHists:
        doGenMetVsPFMetSystematic(signalHists)

    return signalHists
