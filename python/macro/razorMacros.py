## Machinery for inclusive razor analysis

import sys,os
import copy
from array import array
import numpy as np
import ROOT as rt
from collections import namedtuple

#local imports
from RunCombine import exec_me
from framework import Config
import macro
import plotting
from razorAnalysis import *
from razorWeights import *

###########################################
### RAZOR FIT
###########################################

def runFitAndToys(fitDir, boxName, lumi, dataName, dataDir='./', config='config/run2.config', sideband=False, numToys=4000, noStat=False):
    #make folder
    if not os.path.isdir(fitDir):
        exec_me('mkdir -p '+fitDir, False)
    #make RooDataSet
    exec_me('python python/DustinTuple2RooDataSet.py -b '+boxName+' -c '+config+' -l '+str(lumi)+' --data -d '+fitDir+' '+dataDir+'/'+dataName+'.root', False)
    #do fit
    if not sideband:
        exec_me('python python/BinnedFit.py -c '+config+' -d '+fitDir+' -l '+str(lumi)+' -b '+boxName+' --data '+fitDir+'/'+dataName+'_lumi-'+('%1.3f' % (lumi*1.0/1000))+'_0-3btag_'+boxName+'.root', False)
    else:
        exec_me('python python/BinnedFit.py -c '+config+' -d '+fitDir+' -l '+str(lumi)+' -b '+boxName+' --data --fit-region LowMR,LowRsq '+fitDir+'/'+dataName+'_lumi-'+('%1.3f' % (lumi*1.0/1000))+'_0-3btag_'+boxName+'.root', False)
    #run toys
    exec_me('python python/RunToys.py -b '+boxName+' -c '+config+' -i '+fitDir+'/BinnedFitResults_'+boxName+'.root -d '+fitDir+' -t '+str(numToys)+((noStat)*" --no-stat"), False)

def runFitAndToysMC(fitDir, boxName, lumi, fileNames, mcDir='./', config='config/run2.config', sideband=False, numToys=4000, noStat=False):
    #make folder
    if not os.path.isdir(fitDir):
        exec_me('mkdir -p '+fitDir, False)
    #make RooDataSet
    exec_me('python python/DustinTuple2RooDataSet.py -w -b '+boxName+' -c '+config+' -l '+str(lumi)+' -d '+fitDir+' '+' '.join(fileNames), False)
    #do fit
    if not sideband:
        exec_me('python python/BinnedFit.py -c '+config+' -d '+fitDir+' -l '+str(lumi)+' -b '+boxName+' '+fitDir+'/RazorInclusive_SMCocktail_weighted_lumi-'+('%1.3f' % (lumi*1.0/1000))+'_0-3btag_'+boxName+'.root', False)
    else:
        exec_me('python python/BinnedFit.py -c '+config+' -d '+fitDir+' -l '+str(lumi)+' -b '+boxName+' --fit-region LowMR,LowRsq '+fitDir+'/RazorInclusive_SMCocktail_weighted_lumi-'+('%1.3f' % (lumi*1.0/1000))+'_0-3btag_'+boxName+'.root', False)
    #run toys
    exec_me('python python/RunToys.py -b '+boxName+' -c '+config+' -i '+fitDir+'/BinnedFitResults_'+boxName+'.root -d '+fitDir+' -t '+str(numToys)+((noStat)*" --no-stat"), False)

def get2DNSigmaHistogram(data, bins, fitToyFiles, boxName, btags=-1, btagsMax=3, debugLevel=0):
    print "Making Nsigma histogram using fit information"
    #set up histogram
    nsigma = data.Clone(data.GetName()+"Nsigma")
    nsigma.Reset()

    #load fit information, including toys
    toyFile = rt.TFile.Open(fitToyFiles[boxName])
    assert toyFile
    if debugLevel > 0: print "Opened file",fitToyFiles[boxName],"to get Bayesian toy results"
    toyTree = toyFile.Get("myTree")
    assert toyTree
    if debugLevel > 0: print "Got tree myTree"

    #make dummy options tuple
    Opt = namedtuple("Opt", ["printErrors","noStat"])
    options = Opt(False,False)

    z = range(btagsMax+2) #b-tag binning
    if btags < 0: #be inclusive in b-tags
        zmin = 0
        zmax = len(z)-1
    else:
        zmin = btags+1
        zmax = btags+1

    from PlotFit import getNsigma2D
    nsigma = getNsigma2D(nsigma, data, toyTree, options, 'yx', 0, len(bins['MR'])-1, 0, len(bins['Rsq'])-1, zmin, zmax, bins['MR'], bins['Rsq'], z)
    return nsigma

def getFitCorrelationMatrix(config, box, fitToyFile, debugLevel=0):
    print "Creating fit correlation matrix"

    cfg = Config.Config(config)
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning

    #load fit information, including toys
    toyFile = rt.TFile.Open(fitToyFile)
    assert toyFile
    if debugLevel > 0: print "Opened file",fitToyFile,"to get Bayesian toy results"
    toyTree = toyFile.Get("myTree")
    assert toyTree
    if debugLevel > 0: print "Got tree myTree"

    from PlotFit import getCorrelationMatrix
    corrMatrix = getCorrelationMatrix(toyTree,'yx',0,len(x)-1,0,len(y)-1,0,len(z)-1,x,y,z)
    return corrMatrix

def import2DRazorFitHistograms(hists, bins, fitToyFile, c, dataName="Data", btags=-1, btagsMax=3, debugLevel=0, noStat=False,
        force=False):
    FIT_SCALE_FACTOR = 1.0
    print "Loading fit histograms"
    if noStat: 
        print "Using only systematic errors on fit points"
        if 'noStat' not in fitToyFile:
            fitToyFile = fitToyFile.replace('Bayes','Bayes_varyN_noStat')
    #sanity check
    if "Fit" in hists and not force:
        print "Error in import2DFitHistograms: fit histogram already exists!"
        return
    #make histogram for fit result
    hists["Fit"] = {}
    for v in ["MR","Rsq",("MR","Rsq")]:
        hists["Fit"][v] = next(hists.itervalues())[v].Clone(next(hists.itervalues())[v].GetName()+"Fit")
        hists["Fit"][v].Reset()
        hists['Fit'][v].SetDirectory(0)

    #load fit information, including toys
    toyFile = rt.TFile.Open(fitToyFile)
    print "Opening fit file",fitToyFile
    assert toyFile
    if debugLevel > 0: print "Opened file",fitToyFile,"to get Bayesian toy results"
    toyTree = toyFile.Get("myTree")
    assert toyTree
    if debugLevel > 0: print "Got tree myTree"

    #get uncertainties on fit prediction
    #using code imported from Javier's plotting script
    from PlotFit import getBinSumDicts, getBestFitRms, getErrors1D
    z = range(btagsMax+2) # b-tag binning
    if btags < 0: #be inclusive in b-tags
        zmin = 0
        zmax = len(z)-1
    else:
        zmin = btags+1
        zmax = btags+1

    #make dummy options tuple
    Opt = namedtuple("Opt", ["printErrors","noStat"])
    options = Opt(False, noStat)

    #store best fit and uncertainty in each bin
    #1D
    if debugLevel > 0:
        print bins
    binSumDict = getBinSumDicts('x', 0,len(bins['MR'])-1,0,len(bins['Rsq'])-1,zmin,zmax,bins['MR'],bins['Rsq'],z)
    if debugLevel > 0:
        print binSumDict
    for i,sumName in binSumDict.iteritems():
        if dataName in hists: nObs = hists[dataName]["MR"].GetBinContent(i)
        else: nObs = 0
        bestFit, rms, pvalue, nsigma, c = getBestFitRms(toyTree,sumName,nObs,c,options,"")
        hists["Fit"]["MR"].SetBinContent(i,bestFit * FIT_SCALE_FACTOR)
        hists["Fit"]["MR"].SetBinError(i,rms * FIT_SCALE_FACTOR)
    binSumDict = getBinSumDicts('y', 0, len(bins['MR'])-1,0,len(bins['Rsq'])-1,zmin,zmax,bins['MR'],bins['Rsq'],z)
    for i,sumName in binSumDict.iteritems():
        if dataName in hists: nObs = hists[dataName]["Rsq"].GetBinContent(i)
        else: nObs = 0
        bestFit, rms, pvalue, nsigma, c = getBestFitRms(toyTree,sumName,nObs,c,options,"")
        hists["Fit"]["Rsq"].SetBinContent(i,bestFit * FIT_SCALE_FACTOR)
        hists["Fit"]["Rsq"].SetBinError(i,rms * FIT_SCALE_FACTOR)
    #2D
    binSumDict = getBinSumDicts('yx', 0, len(bins['MR'])-1,0,len(bins['Rsq'])-1,zmin,zmax,bins['MR'],bins['Rsq'],z)
    for (i,j),sumName in binSumDict.iteritems():
        if dataName in hists: nObs = hists[dataName][("MR","Rsq")].GetBinContent(i,j)
        else: nObs = 0
        bestFit, rms, pvalue, nsigma, c = getBestFitRms(toyTree,sumName,nObs,c,options=options,plotName="")
        hists["Fit"][("MR","Rsq")].SetBinContent(i,j,bestFit * FIT_SCALE_FACTOR)
        hists["Fit"][("MR","Rsq")].SetBinError(i,j,rms * FIT_SCALE_FACTOR)

def get3DRazorFitHistogram(configFile, fitToyFile, boxName, debugLevel=0):
    print "Getting 3D fit histogram from",fitToyFile
    #load binning
    cfg = Config.Config(configFile)
    binsX = array('d', cfg.getBinning(boxName)[0])
    binsY = array('d', cfg.getBinning(boxName)[1])
    binsZ = array('d', cfg.getBinning(boxName)[2])
    nBins = (len(binsX)-1)*(len(binsY)-1)*(len(binsZ)-1)

    #load functions from Javier
    from PlotFit import getBinSumDicts, getBestFitRms

    #make fit histogram 
    fitHist3D = rt.TH3F("fit"+boxName, "fit"+boxName, len(binsX)-1, binsX, len(binsY)-1, binsY, len(binsZ)-1, binsZ)

    #load fit information, including toys
    toyFile = rt.TFile.Open(fitToyFile)
    assert toyFile
    if debugLevel > 0: print "Opened file",fitToyFile,"to get Bayesian toy results"
    toyTree = toyFile.Get("myTree")
    assert toyTree
    if debugLevel > 0: print "Got tree myTree"

    #make dummy options tuple
    Opt = namedtuple("Opt", ["printErrors","noStat"])
    options = Opt(False, False)

    #store best fit and uncertainty in each bin
    binSumDict = getBinSumDicts('zyx', 0,len(binsX)-1,0,len(binsY)-1,0,len(binsZ)-1,binsX,binsY,binsZ)
    c = rt.TCanvas("placeholder","placeholder",800,600)
    for (i,j,k),sumName in binSumDict.iteritems():
        bestFit, rms, pvalue, nsigma, c = getBestFitRms(toyTree,sumName,0,c,options,"")
        fitHist3D.SetBinContent(i,j,k,bestFit)
        fitHist3D.SetBinError(i,j,k,rms)

    return fitHist3D

def makeRazor3DTable(hist, boxName, signalHist=None, signalName="T1bbbb"):
    """Print latex table with prediction and uncertainty in each bin"""
    xbinLowEdges = []
    xbinUpEdges = []
    ybinLowEdges = []
    ybinUpEdges = []
    zbinLowEdges = []
    zbinUpEdges = []
    predictions = []
    uncerts = []
    signal = []
    #for each bin, get values for all table columns
    for bx in range(1, hist.GetNbinsX()+1):
        for by in range(1, hist.GetNbinsY()+1):
            for bz in range(1, hist.GetNbinsZ()+1):
                xbinLowEdges.append('%.0f' % (hist.GetXaxis().GetBinLowEdge(bx)))
                xbinUpEdges.append('%.0f' % (hist.GetXaxis().GetBinUpEdge(bx)))
                ybinLowEdges.append(str(hist.GetYaxis().GetBinLowEdge(by)))
                ybinUpEdges.append(str(hist.GetYaxis().GetBinUpEdge(by)))
                zbinLowEdges.append('%.0f' % (hist.GetZaxis().GetBinLowEdge(bz)))
                zbinUpEdges.append('%.0f' % (hist.GetZaxis().GetBinUpEdge(bz)))
                pred = hist.GetBinContent(bx,by,bz)
                predictions.append('%.2f' % (pred))
                unc = hist.GetBinError(bx,by,bz)
                if pred != 0:
                    uncerts.append('%.1f\\%%' % (100*(unc/pred)))
                else:
                    uncerts.append('%.2f' % (unc))
                if signalHist is not None:
                    signal.append('%.2f' % (signalHist.GetBinContent(bx,by,bz)))
    xRanges = [low+'-'+high for (low, high) in zip(xbinLowEdges, xbinUpEdges)]
    yRanges = [low+'-'+high for (low, high) in zip(ybinLowEdges, ybinUpEdges)]
    zRanges = copy.copy(zbinLowEdges)
    headers=["$M_R$", "$R^2$", "B-tags", "Prediction", "Uncertainty"]
    cols = [xRanges, yRanges, zRanges, predictions, uncerts]
    if signalHist is not None:
        headers.append(signalName)
        cols.append(signal)
    plotting.table_basic(headers, cols, caption="Fit prediction for the "+boxName+" box", printstr="razorFitTable"+boxName)

def makeRazor2DTable(pred, boxName, nsigma=None, obs=None, mcNames=[], mcHists=[], btags=-1, unrollBins=None, useMCFitSys=False, printdir='.', emptyBinErrs=None, listAllMC=False):
    """Print latex table with prediction and uncertainty in each bin"""
    xbinLowEdges = []
    xbinUpEdges = []
    ybinLowEdges = []
    ybinUpEdges = []
    zbinLowEdges = []
    preds = []
    obses = []
    nsigmas = []
    mcs = []
    totalMCs = []
    for m in mcNames:
        mcs.append([])
    #for each bin, get values for all table columns
    if pred is not None:
        refHist = pred
    elif len(mcHists) > 0:
        refHist = mcHists[0]
    else:
        refHist = obs
    if unrollBins is None: #TH2F case
        for bx in range(1, refHist.GetNbinsX()+1):
            for by in range(1, refHist.GetNbinsY()+1):
                xbinLowEdges.append('%.0f' % (refHist.GetXaxis().GetBinLowEdge(bx)))
                xbinUpEdges.append('%.0f' % (refHist.GetXaxis().GetBinUpEdge(bx)))
                ybinLowEdges.append(str(refHist.GetYaxis().GetBinLowEdge(by)))
                ybinUpEdges.append(str(refHist.GetYaxis().GetBinUpEdge(by)))
                zbinLowEdges.append('%.0f' % (btags))
                if pred is not None:
                    prediction = pred.GetBinContent(bx,by)
                    uncert = pred.GetBinError(bx,by)
                    preds.append('%.2f $\\pm$ %.2f' % (prediction, uncert))
                if obs is not None: 
                    observed = obs.GetBinContent(bx,by)
                    obses.append('%.2f' % (observed))
                if nsigma is not None:
                    nsig = nsigma.GetBinContent(bx,by)
                    nsigmas.append('%.2f' % (nsig))
                totalMC = 0.0
                totalMCErr = 0.0
                for i in range(len(mcNames)):
                    mcs[i].append('%.3f $\\pm$ %.3f' % (max(0,mcHists[i].GetBinContent(bx,by)),mcHists[i].GetBinError(bx,by)))
                    totalMC += mcHists[i].GetBinContent(bx,by)
                    totalMCErr = ( totalMCErr**2 + (mcHists[i].GetBinError(bx,by))**2 )**(0.5)
                if len(mcNames) > 0:
                    if useMCFitSys and pred is not None: #add (MC-fit) in quadrature with MC error
                        totalMCErr = ( totalMCErr**2 + (totalMC-prediction)**2 )**(0.5)
                    totalMCs.append('%.2f $\\pm$ %.2f' % (max(0,totalMC), totalMCErr))
    else: #some bins are merged TH2Poly style
        print "Merging bins according to unrolled binning provided"
        if nsigma is not None:
            print "Warning: nsigma histogram not supported for TeX table with merged bins!"
        xbinLowEdges, xbinUpEdges, ybinLowEdges, ybinUpEdges = macro.getBinBoundariesFromColumns(unrollBins[0], unrollBins[1]) #get the bin low/high edges
        xbinLowEdges, xbinUpEdges, ybinLowEdges, ybinUpEdges = [map(str, bins) for bins in [xbinLowEdges, xbinUpEdges, ybinLowEdges, ybinUpEdges]] 
        #apply unrolled binning
        mergedFit = plotting.unroll2DHistograms([pred], unrollBins[0], unrollBins[1])[0]
        mergedObs = plotting.unroll2DHistograms([obs], unrollBins[0], unrollBins[1])[0] 
        mergedMCs = plotting.unroll2DHistograms(mcHists, unrollBins[0], unrollBins[1]) 
        if pred is not None:
            refHist = mergedFit
        elif len(mcNames) > 0:
            refHist = mergedMCs[0]
        else:
            refHist = mergedObs
        for bx in range(1, refHist.GetNbinsX()+1):
                zbinLowEdges.append('%.0f' % (btags))
                if pred is not None:
                    prediction = mergedFit.GetBinContent(bx)
                    uncert = mergedFit.GetBinError(bx)
                    preds.append('%.2f $\\pm$ %.2f' % (prediction, uncert))
                if obs is not None: 
                    observed = mergedObs.GetBinContent(bx)
                    obses.append('%.2f' % (observed))
                totalMC = 0.0
                totalMCErr = 0.0
                for i in range(len(mcNames)):
                    mcErr = mergedMCs[i].GetBinError(bx)
                    if emptyBinErrs is not None:
                        emptyErr = emptyBinErrs[mcNames[i]][bx]
                        mcErr = (mcErr*mcErr + emptyErr*emptyErr)**(0.5)
                    mcs[i].append('%.3f $\\pm$ %.3f' % (max(0,mergedMCs[i].GetBinContent(bx)), mcErr))
                    totalMC += mergedMCs[i].GetBinContent(bx)
                    totalMCErr = ( totalMCErr**2 + mcErr**2 )**(0.5)
                if len(mcNames) > 0:
                    if useMCFitSys and pred is not None: #add (MC-fit) in quadrature with MC error
                        totalMCErr = ( totalMCErr**2 + (totalMC-prediction)**2 )**(0.5)
                    totalMCs.append('%.2f $\\pm$ %.2f' % (max(0,totalMC), totalMCErr))
        
    xRanges = [low+'-'+high for (low, high) in zip(xbinLowEdges, xbinUpEdges)]
    yRanges = [low+'-'+high for (low, high) in zip(ybinLowEdges, ybinUpEdges)]
    zRanges = copy.copy(zbinLowEdges)
    caption = "Comparison of event yields for the "+boxName+" box"
    label = 'yields'+boxName
    if btags >= 0:
        headers=["$M_R$", "$R^2$", "B-tags"]
        cols = [xRanges, yRanges, zRanges]
        if obs is not None: 
            cols.append(obses)
            headers.append("Observed")
        if not useMCFitSys and pred is not None:
            headers = headers + ["Fit Prediction"]
            cols = cols + [preds]
            if nsigma is not None and len(nsigmas) > 0: 
                cols.append(nsigmas)
                headers.append("Number of sigmas")
        caption += " ("+str(btags)+" b-tags)"
        label += str(btags)+'B'
    else:
        headers=["$M_R$", "$R^2$"]
        cols = [xRanges, yRanges]
        if obs is not None:
            cols.append(obses)
            headers.append("Observed")
        if not useMCFitSys and pred is not None:
            headers = headers + ["Fit Prediction"]
            cols = cols + [preds]
            if nsigma is not None and len(nsigmas) > 0: 
                cols.append(nsigmas)
                headers.append("Number of sigmas")
    if len(mcNames) > 0:
        headers.extend(["MC Prediction"])
        cols.extend([totalMCs])
        if listAllMC:
            for i,n in enumerate(mcNames):
                headers.append(n)
                cols.extend([mcs[i]])
    if listAllMC:
        printstr='razor2DFitTableFull'+boxName+str(btags)+'btag'
    else:
        printstr='razor2DFitTable'+boxName+str(btags)+'btag'
    plotting.table_basic(headers, cols, caption=caption, label=label, printstr=printstr, printdir=printdir, landscape=True, size='normalsize')

###########################################
### BASIC HISTOGRAM FILLING/PLOTTING MACRO
###########################################

def makeControlSampleHists(regionName="TTJetsSingleLepton", filenames={}, samples=[], cutsMC="", cutsData="", bins={}, plotOpts={}, lumiMC=1, lumiData=3000, weightHists={}, sfHists={}, treeName="ControlSampleEvent",dataName="Data", weightOpts=[], shapeErrors=[], miscErrors=[], fitToyFiles=None, boxName=None, btags=-1, blindBins=None, makePlots=True, debugLevel=0, printdir=".", plotDensity=True, sfVars = ("MR","Rsq"), auxSFs={}, dataDrivenQCD=False, unrollBins=(None,None), noFill=False, exportShapeErrs=False, propagateScaleFactorErrs=True, extraWeightOpts={}, extraCuts={}, dataWeightOpts=[]):
    """Basic function for filling histograms and making plots.
        NOTE: for most purposes please call makeControlSampleHistsForAnalysis.
        Arguments:
            regionName: name of the box/bin/control region (used for plot labels)
            filenames: dictionary of process:filename pairs for loading ntuples
            samples: list of samples to process, in the order that they should 
                appear in stacked plots, legends, etc
            cutsMC, cutsData: strings, to be used with TTreeFormula to 
                make selection cuts
            bins: dictionary formatted like { "variable1":[bin0,bin1,bin2,...], 
                "variable2":[bin0,bin1,bin2,...]}
            lumiMC, lumiData: in /pb
            weightHists: dictionary of weight histograms, like that produced by 
                razorWeights.loadWeightHists
            weightOpts: list of strings with directives for applying weights to the MC
            miscErrors: currently unused
            fitToyFiles: optional -- dictionary of boxName:toyFile pairs for loading razor fit results
            doDataShapeUnc: create shape uncertainty histograms for data
    """
    # this is used to customize histogram titles 
    # for certain kinematic variables
    titles = {
        "MR": "M_{R} [GeV]", 
        "Rsq": "R^{2}",
        "mll": "m_{ll} [GeV]",
        "NBJetsMedium" : "Number of B-tagged Jets",
        "NJets80" : "Number of Jets with p_{T} > 80 GeV",
        "NJets40" : "Number of Jets",        
        "lep1.Pt()": "lepton p_{T} [GeV]",
        "lep2.Pt()": "lepton p_{T} [GeV]",
        "lep1.Eta()": "lepton #eta",
        "lep2.Eta()": "lepton #eta",
        }

    #Getting data-driven QCD prediction from high deltaPhi region
    if dataDrivenQCD and not 'QCD' in samples:
        print "Note: ignoring dataDrivenQCD option because QCD is not in the list of samples to process!"
        dataDrivenQCD = False
    if dataDrivenQCD:
        print "\nWill first process the data to obtain QCD prediction"
        cutsForQCDBkg = cutsMC.replace('abs(dPhiRazor) <','abs(dPhiRazor) >')
        cutsForQCDData = cutsData.replace('abs(dPhiRazor) <','abs(dPhiRazor) >')
        samplesForQCD = copy.copy(samples)
        samplesForQCD.remove('QCD')
        weightOptsForQCD = copy.copy(weightOpts)
        if 'GJets' in regionName:
            weightOptsForQCD.append('qcdphoton')
        else:
            cutsForQCDBkg = macro.removeVarCuts(cutsForQCDBkg, 'nBTaggedJets')
            cutsForQCDData = macro.removeVarCuts(cutsForQCDData, 'nBTaggedJets')
        #recursion
        histsForQCD = makeControlSampleHists(
                regionName=regionName+"QCDControlRegion", 
                filenames=filenames, samples=samplesForQCD, 
                cutsMC=cutsForQCDBkg, cutsData=cutsForQCDData, 
                bins=bins, lumiMC=lumiMC, lumiData=lumiData, 
                weightHists=weightHists, sfHists=sfHists, 
                treeName=treeName, dataName="QCD", 
                weightOpts=weightOptsForQCD, boxName=boxName, 
                btags=btags, debugLevel=debugLevel, 
                sfVars=sfVars, auxSFs=auxSFs, makePlots=False, 
                noFill=noFill, extraCuts=extraCuts, 
                extraWeightOpts=extraWeightOpts, exportShapeErrs=True,
                dataWeightOpts=copy.copy(dataWeightOpts))
        print "Now back to our signal region..."
    else:
        histsForQCD = None

    ##Get plotting options (for customizing plot behavior)
    special = ""
    #set log scale
    if "logx" in plotOpts: logx = plotOpts["logx"]
    else: logx = True
    if "ymin" in plotOpts: ymin = plotOpts["ymin"]
    else: ymin = 0.1
    #allow disabling comment string (normally written at the top of each plot)
    if "comment" in plotOpts: comment = plotOpts["comment"]
    else: comment = True
    #use sideband fit result 
    if "sideband" in plotOpts:
        if plotOpts['sideband']:
            special += 'sideband'
        else:
            special += 'full'
    #SUS-15-004 style for unrolled plots
    if "SUS15004" in plotOpts and plotOpts["SUS15004"]:
        special += "SUS15004"
    elif "SUS15004CR" in plotOpts and plotOpts["SUS15004CR"]:
        special += "CR15004"

    #set up files and trees
    if debugLevel > 0: print ""
    inputs = filenames
    files = {name:rt.TFile.Open(inputs[name]) for name in inputs} #get input files
    for name in inputs: 
        print "Checking",name,"..."
        assert files[name] #fails if ROOT could not open the file
        if debugLevel > 0: print "Opened file",inputs[name]
    trees = macro.makeTreeDict(files, treeName, debugLevel)

    #extra scale factors can be defined per-process or globally for all processes
    #(the auxSFs dictionary will have an extra layer in that case)
    auxSFsPerProcess = False
    for key,item in auxSFs.iteritems():
        if isinstance(item, dict):
            auxSFsPerProcess = True
            break

    #split histograms into those that can be applied via per-event weights, and those that require further processing
    sfShapes, otherShapes = splitShapeErrorsByType(shapeErrors)
    #get scale factor options for each shape uncertainty
    shapeAuxSFs = {}
    for shape in sfShapes:
        if not isinstance(shape,basestring): #tuple (shape, [list of processes])
            curShape = shape[0]
        else:
            curShape = shape
        shapeAuxSFs[curShape+'Up'] = copy.deepcopy(auxSFs)
        getAuxSFsForErrorOpt(auxSFs=shapeAuxSFs[curShape+'Up'], errorOpt=curShape+'Up',
                auxSFsPerProcess=auxSFsPerProcess)
        shapeAuxSFs[curShape+'Down'] = copy.deepcopy(auxSFs)
        getAuxSFsForErrorOpt(auxSFs=shapeAuxSFs[curShape+'Down'], errorOpt=curShape+'Down',
                auxSFsPerProcess=auxSFsPerProcess)
    print "\nThese shape uncertainties will be applied via event-level weights:"
    print sfShapes
    print "\nOther shape uncertainties:"
    print otherShapes
    print "\n These event reweighting options will be used:"
    print weightOpts

    #define histograms to fill
    hists,shapeHists = macro.setupHistograms(regionName, inputs, 
            samples, bins, titles, shapeErrors, dataName)
    listOfVars = hists.itervalues().next().keys() #list of the variable names
    
    samplesToUse = copy.copy(samples) #MC samples to process
    if dataDrivenQCD:
        #do not process QCD as a MC sample
        samplesToUse.remove('QCD')
        #store QCD result from above 
        postprocessQCDHists(hists, shapeHists,
                histsForQCD, weightHists, regionName, btags, debugLevel)
    #use QCD extrapolation on data
    auxSFsData = {}
    auxSFHistsData={}
    if 'qcdphoton' in map(str.lower, weightOpts):
        auxSFsData = auxSFs["QCD"]
        samplesToUse = [] #skip running MC samples
        del auxSFs["QCD"] 
        auxSFHistsData = {name:sfHists[name] for name in auxSFsData}
            
    #sum of weights histograms are needed to compute 
    #factorization/renormalization scale uncertainties
    if ('facscale' in shapeErrors or 'renscale' in shapeErrors
            or 'facrenscale' in shapeErrors):
        print "Loading sum of weights histograms for all processes"
        weightHistsPerProcess = makeWeightHistDict(files, debugLevel)
    else:
        weightHistsPerProcess = {}
    print "\n These scale factors will be used for data:"
    print auxSFsData

    #loop over data separately from MC
    if dataName in trees:
        print("\nData:")
        macro.loopTree(trees[dataName], weightF=weight_data, cuts=cutsData, 
                hists=hists[dataName], weightHists=weightHists, 
                auxSFs=auxSFsData, auxSFHists=auxSFHistsData, noFill=noFill, 
                debugLevel=debugLevel, weightOpts=dataWeightOpts) 

    print("\nMC:")
    if debugLevel > 0:
        print "\nMisc SF hists to use:"
        print auxSFs
    macro.loopTrees(trees, weightF=weight_mc, cuts=cutsMC, hists={name:hists[name] for name in samplesToUse}, weightHists=weightHists, sfHists=sfHists, scale=lumiData*1.0/lumiMC, weightOpts=weightOpts, sfVars=sfVars, statErrOnly=False, auxSFs=auxSFs, shapeHists=shapeHists, shapeNames=sfShapes, shapeAuxSFs=shapeAuxSFs, noFill=noFill, propagateScaleFactorErrs=propagateScaleFactorErrs, debugLevel=debugLevel, extraCuts=extraCuts, extraWeightOpts=extraWeightOpts, weightHistsPerProcess=weightHistsPerProcess) 

    #Some shape uncertainties cannot easily be computed in the same pass as
    #the central histogram values (usually because the selection cuts are different).
    #We compute these separately by looping over the MC again.
    for shape in otherShapes:
        if not isinstance(shape,basestring): #tuple (shape, [list of processes])
            shapeSamplesToUse = filter(lambda n: n in shape[1], samplesToUse)
            curShape = shape[0]
        else:
            shapeSamplesToUse = samplesToUse
            curShape = shape
        print "\n"+curShape,"Up:"
        #get any scale factor histograms needed to apply this up variation
        auxSFsToUse = copy.deepcopy(auxSFs)
        getAuxSFsForErrorOpt(auxSFs=auxSFsToUse, errorOpt=curShape+"Up",
                auxSFsPerProcess=auxSFsPerProcess)
        if debugLevel > 0:
            print "Auxiliary SF hists to use:"
            print auxSFsToUse
        macro.loopTrees(trees, weightF=weight_mc, cuts=cutsMC, hists={name:shapeHists[name][curShape+"Up"] for name in shapeSamplesToUse}, weightHists=weightHists, sfHists=sfHists, scale=lumiData*1.0/lumiMC, weightOpts=weightOpts, errorOpt=curShape+"Up", boxName=boxName, sfVars=sfVars, statErrOnly=True, auxSFs=auxSFsToUse, noFill=noFill, propagateScaleFactorErrs=False, debugLevel=debugLevel, extraCuts=extraCuts, extraWeightOpts=extraWeightOpts)
        print "\n"+curShape,"Down:"
        #get any scale factor histograms needed to apply this down variation
        auxSFsToUse = copy.deepcopy(auxSFs)
        getAuxSFsForErrorOpt(auxSFs=auxSFsToUse, errorOpt=curShape+"Down",
                auxSFsPerProcess=auxSFsPerProcess)
        if debugLevel > 0:
            print "Auxiliary SF hists to use:"
            print auxSFsToUse
        macro.loopTrees(trees, weightF=weight_mc, cuts=cutsMC, hists={name:shapeHists[name][curShape+"Down"] for name in shapeSamplesToUse}, weightHists=weightHists, sfHists=sfHists, scale=lumiData*1.0/lumiMC, weightOpts=weightOpts, errorOpt=curShape+"Down", boxName=boxName, sfVars=sfVars, statErrOnly=True, auxSFs=auxSFsToUse, noFill=noFill, propagateScaleFactorErrs=False, debugLevel=debugLevel, extraCuts=extraCuts, extraWeightOpts=extraWeightOpts)

    #we either save the shape histograms separately, or propagate the uncertainties
    #directly to the central histograms.
    if exportShapeErrs:
        hists['Sys'] = shapeHists
        for proc in hists['Sys']:
            print proc
            for sys in hists['Sys'][proc]:
                print "   ",sys
    else:
        macro.propagateShapeSystematics(hists, samples, bins, shapeHists, shapeErrors, miscErrors, boxName, debugLevel=debugLevel)

    c = rt.TCanvas(regionName+"c", regionName+"c", 800, 600)

    #optionally compare with fit
    nsigmaFitData = None
    nsigmaFitMC = None
    dataForTable = None
    if fitToyFiles is not None and boxName in fitToyFiles and fitToyFiles[boxName] is not None and ("MR","Rsq") in bins:
        noFitStat=True
        if boxName in ['DiJet', 'LeptonJet']:
            btagsMax = 2
        else:
            btagsMax = 3
        print "Ignoring statistical uncertainty on fit prediction (except for nsigma plot)."
        import2DRazorFitHistograms(hists, bins, fitToyFiles[boxName], c, dataName, btags, btagsMax, debugLevel, noStat=noFitStat)
        print "Making TeX table with predictions in each analysis bin"
        if dataName in hists: 
            print "Including observed data yields in table"
            dataForTable=hists[dataName][("MR","Rsq")]
            nsigmaFitData = get2DNSigmaHistogram(hists[dataName][("MR","Rsq")], bins, fitToyFiles, boxName, btags, btagsMax, debugLevel)
        makeRazor2DTable(pred=hists["Fit"][("MR","Rsq")], obs=dataForTable,
                nsigma=nsigmaFitData, mcNames=samples, mcHists=[hists[s][("MR","Rsq")] for s in samples], boxName=boxName, btags=btags, unrollBins=unrollBins, useMCFitSys=False, printdir=printdir)
        makeRazor2DTable(pred=hists["Fit"][("MR","Rsq")], obs=dataForTable,
                nsigma=nsigmaFitData, mcNames=samples, mcHists=[hists[s][("MR","Rsq")] for s in samples], boxName=boxName, btags=btags, unrollBins=unrollBins, printdir=printdir, listAllMC=True)

    #print histograms
    rt.SetOwnership(c, False)
    if makePlots:
        histsCopy = { proc:hs.copy() for proc, hs in hists.iteritems() }
        if 'combineBackgrounds' in plotOpts:
            print "Combining backgrounds according to dictionary",plotOpts['combineBackgrounds']
            macro.combineBackgroundHists(histsCopy, plotOpts['combineBackgrounds'], 
                    listOfVars, debugLevel=debugLevel)
            if 'combineSamples' in plotOpts:
                samples = plotOpts['combineSamples']
            else: #get new list of samples, possibly not keeping the same ordering
                samples = macro.combineBackgroundNames(samples, plotOpts["combineBackgrounds"])
        macro.basicPrint(histsCopy, mcNames=samples, varList=listOfVars, c=c, printName=regionName, logx=logx, dataName=dataName, ymin=ymin, comment=comment, lumistr=('%.1f' % (lumiData*1.0/1000))+" fb^{-1}", boxName=boxName, btags=btags, blindBins=blindBins, nsigmaFitData=nsigmaFitData, printdir=printdir, doDensity=plotDensity, special=special, unrollBins=unrollBins, vartitles=titles, debugLevel=debugLevel)

    #close files and return
    for f in files: files[f].Close()
    return hists

def makeControlSampleHistsForAnalysis(analysis, plotOpts={}, sfHists={}, sfVars=("MR","Rsq"), 
        treeName="ControlSampleEvent", shapeErrors=[], fitToyFiles=None, boxName=None, btags=-1, 
        blindBins=None, makePlots=True, printdir=".", auxSFs={}, dataDrivenQCD=False, noFill=False, 
        exportShapeErrs=False, propagateScaleFactorErrs=True, lumiMC=1, debugLevel=0):
    """This is the function that should be used to fill histograms for an 
        analysis control/signal region.
        Arguments:
            analysis: macro.razorAnalysis.Analysis object containing 
                control/signal region configuration 
            plotOpts: dictionary containing special plotting options, if any.
                TODO: document the supported options
            sfHists: dictionary of the form
                    { "name1":hist1, ... }
                where the values hist1,... are scale factor histograms
            sfVars: string or tuple indicating which variable(s) should be used 
                when applying MC scale factors
            treeName: name of the input ROOT tree
            shapeErrors: list of uncertainties for which up/down shapes should be computed.
                There are three acceptable formats for uncertainties:
                1) string -- name of the uncertainty 
                    (for uncertainties affecting all physics processes)
                2) tuple of strings ("unc","proc")
                    (for uncertainties affecting a single physics process)
                3) tuple ("unc", ["proc1", "proc2", ...])
                    (for uncertainties affecting several physics processes)
            fitToyFiles: paths to files containing fit results computed with toys
            boxName: name of the analysis box (razor signal region only)
            btags: number of b-tags (razor signal region only)
            blindBins: list of 2-tuples [(x,y),...] indicating bins that should be blinded
            makePlots: bool.  If True, results will be plotted.
            printdir: directory where plots will be stored.
            auxSFs: dictionary specifying additional scale factor options.
                Two formats are supported:
                1) { "scaleFactorName":("varName","cut"), ... }
                    where
                    -"scaleFactorName" is the desired key in the sfHists dictionary
                    -"varName" is the variable to be reweighted
                    -"cut" is a string a la TTree::Draw specifying which events to reweight
                2) { "proc1":{ ... }, "proc2":{ ... }, ... }
                    where "proc1", etc, are physics process names
                    and { ... } are dictionaries formatted like in (1)
            dataDrivenQCD: bool.  if True, data-driven QCD prediction will be obtained
            noFill: if True, set up histograms but do not fill them
                (note: may cause plotting code to crash)
            exportShapeErrs: if True, up/down shape uncertainties are saved as separate histograms.
                If False, they will be propagated to the central histogram.
            propagateScaleFactorErrs: if True, uncertainties from applied scale factors
                will be propagated to the filled histograms
            debugLevel: 0 for standard mode, 1 for verbose mode, 2 for debug mode
        Returns: dictionary of the form 
            { "process1":{"var1":hist1, "var2":hist2, ...}, "process2":{...} }
        containing the desired set of histograms, filled for each physics process."""
            

    return makeControlSampleHists( regionName=analysis.region, filenames=analysis.filenames, 
            samples=analysis.samples, cutsMC=analysis.cutsMC, cutsData=analysis.cutsData, 
            bins=analysis.binning, plotOpts=plotOpts, lumiData=analysis.lumi, 
            weightHists=analysis.weightHists, sfHists=sfHists, treeName=treeName, 
            weightOpts=analysis.weightOpts, shapeErrors=shapeErrors, fitToyFiles=fitToyFiles, 
            boxName=boxName, btags=btags, blindBins=blindBins, makePlots=makePlots, 
            printdir=printdir, auxSFs=auxSFs, dataDrivenQCD=dataDrivenQCD, 
            unrollBins=analysis.unrollBins, noFill=noFill, exportShapeErrs=exportShapeErrs, 
            extraCuts=analysis.extraCuts, extraWeightOpts=analysis.extraWeightOpts, 
            propagateScaleFactorErrs=propagateScaleFactorErrs, 
            sfVars=sfVars, lumiMC=lumiMC, debugLevel=debugLevel,
            dataWeightOpts=analysis.dataWeightOpts)

def postprocessQCDHists(hists, shapeHists,
        qcdHists, wHists, region, btags, debugLevel=0):
    """Performs extrapolation of QCD prediction to signal region
       and stores results in histogram dictionary.
       -hists: dictionary in which to store histograms
       -shapeHists: dictionary to store up/down histograms
       -qcdHists: histograms to process and store
       -wHists: dictionary of histograms for weighting
       -region: name of the analysis region"""
    # We subtract the backgrounds in data in order to isolate
    # the QCD component. This is not relevant for the gamma+jets
    # case because we already apply the photon purity calculation
    # event by event. 
    isSignalRegion = not ('GJets' in region)
    macro.exportHists(qcdHists, 'qcdHists.root', delete=False)
    if isSignalRegion:
        macro.subtractBkgsInData(process='QCD', hists=qcdHists, 
                dataName='QCD', debugLevel=debugLevel)

    # Merge the two sets of histograms
    for var in hists['QCD']:
        if debugLevel > 0:
            print "\nStoring QCD histogram for",var
        hists['QCD'][var] = qcdHists['QCD'][var].Clone()
        hists['QCD'][var].SetDirectory(0)

    # Note: for signal regions, the histograms produced 
    # here are only guaranteed to be correct for the 2D MR-Rsq
    # distribution. The extrapolation is not done on 1D hists.
    if isSignalRegion:
        var = ('MR','Rsq')
        qcd2DHist = qcdHists['QCD'][var]
        hists['QCD'][var] = (
                makeQCDExtrapolation(qcd2DHist, 
                wHists, region, btags, debugLevel))
        if 'QCD' in shapeHists and 'qcdnormUp' in shapeHists['QCD']:
            for errType in ['qcdnormUp', 'qcdnormDown', 
                    'qcdbtagUp', 'qcdbtagDown',
                    'qcdaltfunctionUp', 'qcdaltfunctionDown']:
                shapeHists['QCD'][errType][var] = (
                        makeQCDExtrapolation(qcd2DHist, wHists, region, 
                            btags, debugLevel, errOpt=errType))

def plotControlSampleHists(regionName="TTJetsSingleLepton", inFile="test.root", samples=[], plotOpts={}, lumiMC=1, lumiData=3000, dataName="Data", boxName=None, btags=-1, blindBins=None, debugLevel=0, printdir=".", plotDensity=True, unrollBins=(None,None), shapeErrors=[], doEmptyBinErrs=False, weightHists=None):
    """Loads the output of makeControlSampleHists from a file and creates plots"""

    titles = {
        "MR": "M_{R} [GeV]", 
        "Rsq": "R^{2}",
        "mll": "m_{ll} [GeV]",
        "NBJetsMedium" : "Number of B-tagged Jets",
        "NJets80" : "Number of Jets with p_{T} > 80 GeV",
        "NJets40" : "Number of Jets",        
        "lep1.Pt()": "Lepton p_{T} [GeV]",
        "lep2.Pt()": "Lepton p_{T} [GeV]",
        "lep1.Eta()": "Lepton #eta",
        "lep2.Eta()": "Lepton #eta",
        }
    if 'Tau' in regionName: 
        titles['lep1.Pt()'] = 'Tau p_{T} [GeV]'

    #load the histograms
    hists = macro.importHists(inFile, debugLevel)

    ##Get plotting options (for customizing plot behavior)
    special = ""
    #set log scale
    if "logx" in plotOpts: logx = plotOpts["logx"]
    else: logx = True
    if "ymin" in plotOpts: ymin = plotOpts["ymin"]
    else: ymin = 0.1
    #allow disabling comment string (normally written at the top of each plot)
    if "comment" in plotOpts: comment = plotOpts["comment"]
    else: comment = True
    #use sideband fit result 
    if "sideband" in plotOpts:
        if plotOpts['sideband']:
            special += 'sideband'
        else:
            special += 'full'
    #SUS-15-004 style for unrolled plots
    if "SUS15004" in plotOpts and plotOpts["SUS15004"]:
        special += "SUS15004"
    elif "SUS15004CR" in plotOpts and plotOpts["SUS15004CR"]:
        special += "CR15004"
    if 'zoom' in plotOpts:
        special += 'zoomratio'

    listOfVars = hists.itervalues().next().keys() #list of the variable names
    
    if 'Sys' in hists: 
        #Shape errors have not been applied; propagate them to histograms them before plotting
        if 'QCD' in samples:
            # add additional b-tag uncertainty on QCD
            makeQCDBtagSysHists(hists, weightHists, boxName, btags)  
        print "Propagating shape uncertainties found in Sys directory"
        shapeHists = hists['Sys']
        miscErrors = []
        del hists['Sys']
        macro.propagateShapeSystematics(hists, samples, listOfVars, shapeHists, shapeErrors, miscErrors, boxName, debugLevel=debugLevel)
    emptyBinErrs = None
    emptyBinErrsByProcess = None
    if doEmptyBinErrs:
        if unrollBins is None or unrollBins[0] is None or unrollBins[1] is None:
            raise ValueError("Need to specify binning in order to compute empty bin errors")
        emptyBinErrs = macro.computeEmptyBinErrs(hists, unrollBins)
        emptyBinErrsByProcess = macro.computeEmptyBinErrs(hists, unrollBins, 
                aggregate=False)
    #optionally combine some background processes
    if 'combineBackgrounds' in plotOpts:
        macro.combineBackgroundHists(hists, plotOpts['combineBackgrounds'], listOfVars, debugLevel=debugLevel)
        if emptyBinErrsByProcess is not None:
            emptyBinErrsByProcess = macro.combineEmptyBinErrs(
                    emptyBinErrsByProcess, plotOpts['combineBackgrounds'])
        if 'combineSamples' in plotOpts:
            samples = plotOpts['combineSamples']
        else: #get new list of samples, possibly not keeping the same ordering
            samples = macro.combineBackgroundNames(samples, plotOpts["combineBackgrounds"])

    c = rt.TCanvas(regionName+"c", regionName+"c", 800, 600)

    #make tex table
    if ("MR","Rsq") in listOfVars:
        dataForTable = None
        fitForTable = None
        print "Making TeX table with predictions in each analysis bin"
        if dataName in hists: 
            dataForTable=hists[dataName][("MR","Rsq")]
            print "Including observed data yields in table"
        if 'Fit' in hists:
            fitForTable=hists['Fit'][('MR','Rsq')]
            print "Including fit predictions in table"
        makeRazor2DTable(pred=fitForTable,obs=dataForTable,
                mcNames=samples, mcHists=[hists[s][("MR","Rsq")] for s in samples], 
                boxName=boxName, btags=btags, unrollBins=unrollBins, 
                useMCFitSys=False, printdir=printdir, emptyBinErrs=emptyBinErrsByProcess)
        makeRazor2DTable(pred=fitForTable,obs=dataForTable,
                mcNames=samples, mcHists=[hists[s][("MR","Rsq")] for s in samples], 
                boxName=boxName, btags=btags, unrollBins=unrollBins, printdir=printdir, 
                emptyBinErrs=emptyBinErrsByProcess, listAllMC=True)

    #print histograms
    rt.SetOwnership(c, False)
    macro.basicPrint(hists, mcNames=samples, varList=listOfVars, c=c, 
            printName=regionName, logx=logx, dataName=dataName, ymin=ymin, 
            comment=comment, lumistr=('%.1f' % (lumiData*1.0/1000))+" fb^{-1}", 
            boxName=boxName, btags=btags, blindBins=blindBins, nsigmaFitData=None, 
            nsigmaFitMC=None, printdir=printdir, doDensity=plotDensity, 
            special=special, unrollBins=unrollBins, vartitles=titles, 
            emptyBinErrs=emptyBinErrs, debugLevel=debugLevel)


#######################################
### MAKE SCALE FACTORS FROM HISTOGRAMS
#######################################

def appendScaleFactors(process="TTJets", hists={}, sfHists={}, 
        var=("MR","Rsq"), dataName="Data", normErrFraction=0.2, 
        lumiData=0, signifThreshold=0., th2PolyXBins=None, useUncertainty=False,
        th2PolyCols=None, debugLevel=0, printdir=".",makePlots=True):
    """Subtract backgrounds and make the data/MC histogram for the given process.
    Also makes up/down histograms corresponding to uncertainty on the background normalization (controlled by the normErrFraction argument).
    normErrFraction can be provided as a dictionary to implement process-dependent
    normalization uncertainties.
    
    process: MC physics process for which scale factors should be computed.  
    if process is not in the histogram collection, will compute Data/Total MC scale factors instead. 
    if process is a list, will compute scale factors for the sum of the processes in the list.
    hists: dictionary of data and MC histograms like that produced by the macro.loopTrees function
    sfHists: dictionary of existing scale factor histograms. the new scale factor histograms will be inserted into this dictionary.
    var: variable or tuple of variables in which scale factors should be computed (usually ("MR","Rsq") is used)
    signifThreshold: scale factors that are within N sigma of 1.0, where N=signifThreshold, are set to 1.
    useUncertainty: if True, ignores the data/MC ratio value and instead 
        fills the output histogram with 1 plus the error on the data/MC ratio.
    th2PolyXBins, th2PolyCols: if non-grid binning is desired, provide here the binning for a TH2Poly with rectangular bins.  th2PolyXBins should be a list of bin edges in the x-direction.  th2PolyCols should be a list of lists containing bin low edges in the y-direction for each x-bin.
    """

    ##Sanity checks
    if debugLevel > 0: 
        print "Scale factor histograms so far:"
        print sfHists
    #warn if this scale factor histogram already exists
    try:
        if process in sfHists:
            print "Warning in appendScaleFactors: ",process," scale factor histogram already exists!  Will overwrite..."
    except TypeError: # this deals with errors when process is a list
        pass
    #warn if the needed histograms are not found
    if dataName not in hists:
        print "Error in appendScaleFactors: target data histogram (",dataName,") was not found!"
        return
    if var not in hists[dataName]:
        print "Error in appendScaleFactors: could not find ",var," in hists[",dataName,"]!"
        return
    #if process is not in the input collection, we will do Data/Total MC scale factors
    doTotalMC = False
    doMultipleProcs = False
    if not isinstance(process, basestring):
        process = '_'.join(process)
        print "Will compute scale factors for",process
        doMultipleProcs = True
    elif process not in hists:
        print "The requested process",process,"was not found in the input collection.  Will compute Data/Total MC scale factors."
        doTotalMC=True
    elif var not in hists[process]:
        print "Error in appendScaleFactors: could not find ",var," in hists[",process,"]!"
        return

    #make the scale factor histogram (clone data hist; later subtract backgrounds and divide by MC hist)
    print "Making scale factor histogram for",process
    sfHists[process] = hists[dataName][var].Clone(process+"ScaleFactors")
    sfHists[process].SetDirectory(0)

    #save yields for debugging
    dataDebug = sfHists[process].Clone()

    #make histos to store uncertainties from MC normalization
    mcSysUncUp = sfHists[process].Clone()
    mcSysUncUp.Reset()
    mcSysUncDown = sfHists[process].Clone()
    mcSysUncDown.Reset()
    if not doTotalMC:
        #subtract backgrounds in data
        bgProcesses = [mcProcess for mcProcess in hists 
                if mcProcess != process and mcProcess not in process.split('_')
                and mcProcess != dataName and mcProcess != "Fit"] 
        for mcProcess in bgProcesses:

            #make sure relevant background histogram exists
            if var not in hists[mcProcess]:
                print "Error in appendScaleFactors: could not find",var," in hists[",mcProcess,"]!  Returning from appendScaleFactors..."
                return
            #subtract it
            print "Subtracting",mcProcess,"from",dataName,"distribution"
            sfHists[process].Add(hists[mcProcess][var], -1) 
            if mcProcess not in sfHists: #if we have not computed scale factors for this process, apply a flat normalization uncertainty to its yield
                try:
                    errFrac = normErrFraction[mcProcess]
                except TypeError: # normErrFraction is not a dictionary
                    errFrac = normErrFraction
                except KeyError: # no norm err specified for this process
                    errFrac = 0.
                print "Process",mcProcess,"has no associated scale factors.  Its normalization will be given a ",int(errFrac*100),"% uncertainty"
                for bx in range(1,mcSysUncUp.GetSize()+1):
                    mcSysUncUp.SetBinContent(bx, ( (mcSysUncUp.GetBinContent(bx))**2 + (errFrac*hists[mcProcess][var].GetBinContent(bx))**2 )**(0.5))
                    mcSysUncDown.SetBinContent(bx, ( (mcSysUncDown.GetBinContent(bx))**2 + ((1/(1+errFrac)-1)*hists[mcProcess][var].GetBinContent(bx))**2 )**(0.5))
        #make up/down systematic error histograms (this is the systematic due to MC background normalization)
        sfHists[process+"NormUp"] = sfHists[process].Clone(process+"ScaleFactorsUp")
        sfHists[process+"NormUp"].SetDirectory(0)
        sfHists[process+"NormUp"].Add(mcSysUncUp)
        sfHists[process+"NormDown"] = sfHists[process].Clone(process+"ScaleFactorsDown")
        sfHists[process+"NormDown"].SetDirectory(0)
        sfHists[process+"NormDown"].Add(mcSysUncDown, -1)
    else: 
        #make total MC histogram
        bgProcesses = [mcProcess for mcProcess in hists if mcProcess != dataName and mcProcess != "Fit"] #get list of non-data, non-fit samples
        mcTotal = hists[bgProcesses[0]][var].Clone()
        mcTotal.Reset()
        for mcProcess in bgProcesses:
            if var not in hists[mcProcess]:
                print "Error in appendScaleFactors: could not find",var," in hists[",mcProcess,"]!  Returning from appendScaleFactors..."
                return
            if debugLevel > 0: print "Adding",mcProcess,"to total MC histogram"
            mcTotal.Add(hists[mcProcess][var])
            if mcProcess not in sfHists: #if we have not computed scale factors for this process, apply a flat normalization uncertainty to its yield
                if debugLevel > 0: print "Process",mcProcess,"has no associated scale factors.  Its normalization will be given a ",int(normErrFraction*100),"% uncertainty"
                for bx in range(1,mcSysUncUp.GetSize()+1):
                    mcSysUncUp.SetBinContent(bx, ( (mcSysUncUp.GetBinContent(bx))**2 + (normErrFraction*hists[mcProcess][var].GetBinContent(bx))**2 )**(0.5))
                    mcSysUncDown.SetBinContent(bx, ( (mcSysUncDown.GetBinContent(bx))**2 + ((1/(1+normErrFraction)-1)*hists[mcProcess][var].GetBinContent(bx))**2 )**(0.5))
        #propagate uncertainties to up/down histograms
        mcTotalUp = mcTotal.Clone(mcTotal.GetName()+'Up')
        mcTotalUp.Add(mcSysUncUp)
        mcTotalDown = mcTotal.Clone(mcTotal.GetName()+'Down')
        mcTotalDown.Add(mcSysUncDown, -1)
        sfHists[process+"NormUp"] = sfHists[process].Clone(process+"ScaleFactorsUp")
        sfHists[process+"NormUp"].SetDirectory(0)
        sfHists[process+"NormDown"] = sfHists[process].Clone(process+"ScaleFactorsDown")
        sfHists[process+"NormDown"].SetDirectory(0)

    #save yields for debugging
    dataDebugSubtr = sfHists[process].Clone()

    ###################### TH2Poly Case
    if th2PolyXBins is not None and th2PolyCols is not None: 
        if doMultipleProcs or useUncertainty:
            raise NotImplementedError("This combination of options isn't implemented")
        print "Converting TH2 to TH2Poly"
        sfHists[process] = macro.makeTH2PolyRatioHist(sfHists[process], hists[process][var], th2PolyXBins, th2PolyCols)
        sfHists[process+"NormUp"] = macro.makeTH2PolyRatioHist(sfHists[process+"NormUp"], hists[process][var], th2PolyXBins, th2PolyCols)
        sfHists[process+"NormDown"] = macro.makeTH2PolyRatioHist(sfHists[process+"NormDown"], hists[process][var], th2PolyXBins, th2PolyCols)
        #zero any negative scale factors
        for bn in range(1, sfHists[process].GetNumberOfBins()+1):
            sfHists[process].SetBinContent(bn, max(0., sfHists[process].GetBinContent(bn)))
            sfHists[process+"NormUp"].SetBinContent(bn, max(0., sfHists[process+"NormUp"].GetBinContent(bn)))
            sfHists[process+"NormDown"].SetBinContent(bn, max(0., sfHists[process+"NormDown"].GetBinContent(bn)))
        #Suppress scale factors consistent with 1
        if signifThreshold > 0:
            print "Ignoring scale factors compatible with 1.0 (",signifThreshold,"sigma significance )"
            for bn in range(1, sfHists[process].GetNumberOfBins()+1):
                thisSysErr = abs(sfHists[process+'NormUp'].GetBinContent(bn)-sfHists[process+'NormDown'].GetBinContent(bn))/2.0
                thisStatErr = sfHists[process].GetBinError(bn)
                thisBinErr = (thisSysErr*thisSysErr + thisStatErr*thisStatErr)**(0.5)
                if thisBinErr == 0: continue
                nsigma = abs(sfHists[process].GetBinContent(bn)-1.0)/thisBinErr
                if nsigma < signifThreshold: 
                    sfHists[process].SetBinContent(bn,1.0)
                    #NOTE: TH2Poly has a bug that causes SetBinError(x) to set the error of bin x+1. Beware!!!
                    sfHists[process].SetBinError(bn-1,0.0)

    ###################### Ordinary TH2 Case
    else:
        #divide data/MC
        if doMultipleProcs:
            procs = process.split('_')
            denom = hists[procs[0]][var].Clone()
            for p in procs[1:]:
                denom.Add(hists[p][var])
            sfHists[process].Divide(denom)
            sfHists[process+"NormUp"].Divide(denom)
            sfHists[process+"NormDown"].Divide(denom)
        elif not doTotalMC:
            sfHists[process].Divide(hists[process][var])
            sfHists[process+"NormUp"].Divide(hists[process][var])
            sfHists[process+"NormDown"].Divide(hists[process][var])
        else:
            sfHists[process].Divide(mcTotal)
            sfHists[process+"NormUp"].Divide(mcTotalUp)
            sfHists[process+"NormDown"].Divide(mcTotalDown)
        #save yields for debugging 
        dataDebugRatio = sfHists[process].Clone()

        if useUncertainty:
            for bx in range(sfHists[process].GetSize()+1):
                print "Setting scale factor bin content {} to {}".format(bx, 1+sfHists[process].GetBinError(bx))
                sfHists[process].SetBinContent(bx, 
                        1 + sfHists[process].GetBinError(bx))

        #zero any negative scale factors
        for bx in range(sfHists[process].GetSize()+1):
            sfHists[process].SetBinContent(bx,max(0., sfHists[process].GetBinContent(bx)))
            sfHists[process+"NormUp"].SetBinContent(bx,max(0., sfHists[process+"NormUp"].GetBinContent(bx)))
            sfHists[process+"NormDown"].SetBinContent(bx,max(0., sfHists[process+"NormDown"].GetBinContent(bx)))

        #suppress scale factors consistent with 1
        if signifThreshold > 0:
            print "Ignoring scale factors compatible with 1.0 (",signifThreshold,"sigma significance )"
            for bx in range(sfHists[process].GetSize()+1):
                thisSysErr = abs(sfHists[process+'NormUp'].GetBinContent(bx)-sfHists[process+'NormDown'].GetBinContent(bx))/2.0
                thisStatErr = sfHists[process].GetBinError(bx)
                thisBinErr = (thisSysErr*thisSysErr + thisStatErr*thisStatErr)**(0.5)
                if thisBinErr == 0: continue
                nsigma = abs(sfHists[process].GetBinContent(bx)-1.0)/thisBinErr
                if nsigma < signifThreshold: 
                    sfHists[process].SetBinContent(bx,1.0)
                    sfHists[process].SetBinError(bx,0.0)

        if debugLevel > 0:
            print "Printing scale factor calculation in each bin:"
            for bx in range(sfHists[process].GetSize()+1):
                print "Bin",bx,":"
                print "   Data yield:     %.0f +/- %.2f"%(dataDebug.GetBinContent(bx),dataDebug.GetBinError(bx))
                print "   Bkg subtracted: %.2f"%(dataDebug.GetBinContent(bx)-dataDebugSubtr.GetBinContent(bx))
                print "   Error on BG-subtracted data: %.2f"%(dataDebugSubtr.GetBinError(bx))
                if process in hists:
                    print "   Scale factor:   %.2f = %.2f / %.2f"%(dataDebugRatio.GetBinContent(bx),dataDebugSubtr.GetBinContent(bx),hists[process][var].GetBinContent(bx))
                print "   Error on scale factor:   %.2f"%(dataDebugRatio.GetBinError(bx))
            print "\nScale factor histograms after adding",process,":"
            print sfHists

    #plot scale factors in 2D (not yet implemented for 1 or 3 dimensions)
    c = rt.TCanvas("c"+process+"SFs", "c", 800, 600)
    if not isinstance(var, basestring) and len(var) == 2 and makePlots==True: 
        plotting.draw2DHist(c, sfHists[process], xtitle=var[0], ytitle=var[1], zmin=0.0, zmax=1.8, printstr=process+"ScaleFactors", lumistr=('%.1f' % (lumiData*1.0/1000))+" fb^{-1}", commentstr="", drawErrs=True, logz=False, numDigits=2, printdir=printdir)

        if not doTotalMC and not th2PolyXBins:
            #print table
            xbinLowEdges = []
            xbinUpEdges = []
            ybinLowEdges = []
            ybinUpEdges = []
            sysUncerts = {mcProcess:[] for mcProcess in bgProcesses}
            statUncerts = []
            sfs = []
            #for each bin, get values for all table columns
            for bx in range(1, sfHists[process].GetNbinsX()+1):
                for by in range(1, sfHists[process].GetNbinsY()+1):
                    xbinLowEdges.append('%.0f' % (sfHists[process].GetXaxis().GetBinLowEdge(bx)))
                    xbinUpEdges.append('%.0f' % (sfHists[process].GetXaxis().GetBinUpEdge(bx)))
                    ybinLowEdges.append(str(sfHists[process].GetYaxis().GetBinLowEdge(by)))
                    ybinUpEdges.append(str(sfHists[process].GetYaxis().GetBinUpEdge(by)))
                    scaleFactor = sfHists[process].GetBinContent(bx,by)
                    sfs.append('%.3f' % (scaleFactor))
                    if scaleFactor > 0:
                        statUncerts.append('%.1f\\%%' % (100*(sfHists[process].GetBinError(bx,by)/scaleFactor)))
                    else: 
                        statUncerts.append('--')
                    for mcProcess in bgProcesses: 
                        try:
                            errFrac = normErrFraction[mcProcess]
                        except TypeError: # normErrFraction is not a dictionary
                            errFrac = normErrFraction
                        except KeyError: # no norm err specified for this process
                            errFrac = 0.
                        dataYield = hists[dataName][var].GetBinContent(bx,by)
                        if dataYield > 0 and scaleFactor > 0:
                            sysUncerts[mcProcess].append('%.1f\\%%' % (100*abs(hists[mcProcess][var].GetBinContent(bx,by)*errFrac*1.0/dataYield/scaleFactor)))
                        else: 
                            sysUncerts[mcProcess].append('--')
            xRanges = [low+'-'+high for (low, high) in zip(xbinLowEdges, xbinUpEdges)]
            yRanges = [low+'-'+high for (low, high) in zip(ybinLowEdges, ybinUpEdges)]
            headers=[var[0], var[1], "Scale factor", "Stat.\\ unc."]
            cols = [xRanges, yRanges, sfs, statUncerts]
            for mcProcess in bgProcesses: 
                headers.extend(["Unc.\\ from "+mcProcess])
                cols.extend([sysUncerts[mcProcess]])
            plotting.table_basic(headers, cols, caption="Scale factors for "+process+" background", printstr="scaleFactorTable"+process, printdir=printdir)

def makeVetoLeptonCorrectionHist(hists={}, var=("MR","Rsq"), dataName="Data", lumiData=0, 
        signifThreshold=0., debugLevel=0, regionName="Veto Lepton", 
        normErrFraction=0.2, sfHists={}, histsToCorrect=None, 
        signalRegionVar=None, doDataOverMC=False, mtEfficiencyHist=None, 
        dPhiEfficiencyHist=None, printdir=".", useUncertainty=False):
    """Compare data and MC in veto lepton control region.  Makes ratio histogram and returns it.  
    Arguments are similar to those for appendScaleFactors()

    doDataOverMC: if True, will save Data/MC scale factors.  if False, will save MC-Data. 
    histsToCorrect: the additive Data-MC corrections will be applied to these histograms, and Corrected/Uncorrected scale factors will be derived from the corrected histogram and saved.  (do not use with doDataOverMC option) 
    signalRegionVar: histsToCorrect (signal region) variable name
    mtEfficiencyHist: histogram of MT cut efficiency, used to correct the control region yields
    dPhiEfficiencyHist: histogram of dPhi cut efficiency, used to correct the control region yields
    useUncertainty: output histogram bins contain 1 + uncertainty on data/MC ratio
    """
    regionNameReduced = regionName.replace(' ','')

    #warn if the needed histograms are not found
    if dataName not in hists:
        print "Error in makeVetoLeptonCorrectionHist: target data histogram (",dataName,") was not found!"
        return
    if var not in hists[dataName]:
        print "Error in makeVetoLeptonCorrectionHist: could not find ",var," in hists[",dataName,"]!"

    #make total MC histograms (up, down, central)
    bgProcesses = [mcProcess for mcProcess in hists if mcProcess != dataName and mcProcess != "Fit"] #get list of non-data, non-fit samples
    mcTotal = hists[dataName][var].Clone()
    mcTotal.Reset()
    #make histos to store uncertainties from MC normalization
    mcSysUncUp = mcTotal.Clone()
    mcSysUncDown = mcTotal.Clone()
    for p in bgProcesses:
        if var not in hists[p]:
            print "Error in makeVetoLeptonCorrectionHist: histogram for",p,"not found!"
            return
        mcTotal.Add(hists[p][var])
        if p not in sfHists: #vary normalization of processes not controlled by scale factors
            for bx in range(1,mcSysUncUp.GetSize()+1):
                if debugLevel > 0:
                    print "Error from",p,"normalization:",normErrFraction*hists[p][var].GetBinContent(bx),"up,",(1/(1+normErrFraction)-1)*hists[p][var].GetBinContent(bx),"down"
                mcSysUncUp.SetBinContent(bx, ( (mcSysUncUp.GetBinContent(bx))**2 + (normErrFraction*hists[p][var].GetBinContent(bx))**2 )**(0.5))
                mcSysUncDown.SetBinContent(bx, ( (mcSysUncDown.GetBinContent(bx))**2 + ((1/(1+normErrFraction)-1)*hists[p][var].GetBinContent(bx))**2 )**(0.5))
    mcTotalUp = mcTotal.Clone(mcTotal.GetName()+'Up')
    mcTotalUp.Add(mcSysUncUp)
    mcTotalDown = mcTotal.Clone(mcTotal.GetName()+'Down')
    mcTotalDown.Add(mcSysUncDown, -1)

    #make signal region histogram that will receive the correction
    if histsToCorrect is not None:
        bgProcessesSR = [mcProcess for mcProcess in histsToCorrect if mcProcess != dataName and mcProcess != "Fit"] #get list of non-data, non-fit samples
        histToCorrect = hists[dataName][var].Clone("histToCorrect"+regionNameReduced) #set up correct binning
        histToCorrect.Reset()
        for p in bgProcessesSR:
            if signalRegionVar not in histsToCorrect[p]:
                print "Error in makeVetoLeptonCorrectionHist: signal region histogram for",p,"not found!"
                return
            histToCorrect.Add(histsToCorrect[p][signalRegionVar])
    else:
        histToCorrect = None

    #output histograms
    vlHists = {}
    if not doDataOverMC: #do MC minus Data
        #subtract MC-data
        print "Making correction histogram for veto lepton control region (MC minus Data)"
        vlHists['Central'] = mcTotal.Clone(regionNameReduced+"Correction")
        vlHists['Up'] = mcTotalUp.Clone(regionNameReduced+"CorrectionUp")
        vlHists['Down'] = mcTotalDown.Clone(regionNameReduced+"CorrectionDown")
        for n,h in vlHists.iteritems():
            h.Add(hists[dataName][var],-1.0)
    else: #do Data/MC
        #divide data by MC
        print "Making correction histogram for veto lepton control region (Data/MC)"
        vlHists['Central'] = hists[dataName][var].Clone(regionNameReduced+"Correction")
        vlHists['Central'].Divide(mcTotal)
        #up histogram
        vlHists['Up'] = hists[dataName][var].Clone(regionNameReduced+"CorrectionUp")
        vlHists['Up'].Divide(mcTotalUp)
        #down histogram
        vlHists['Down'] = hists[dataName][var].Clone(regionNameReduced+"CorrectionDown")
        vlHists['Down'].Divide(mcTotalDown)
    for n,h in vlHists.iteritems():
        h.SetDirectory(0)

    #suppress corrections consistent with 0
    targetVal = 0.0 #check for consistency with this number
    if doDataOverMC:
        targetVal = 1.0
    if signifThreshold > 0:
        print "Ignoring corrections compatible with 0 at",signifThreshold,"sigma significance"
        for bx in range(vlHists['Central'].GetSize()+1):
            thisSysErr = abs(vlHists['Up'].GetBinContent(bx)-vlHists['Down'].GetBinContent(bx))/2.0
            thisStatErr = vlHists['Central'].GetBinError(bx)
            thisBinErr = (thisSysErr*thisSysErr + thisStatErr*thisStatErr)**(0.5)
            if thisBinErr == 0: continue
            nsigma = abs(vlHists['Central'].GetBinContent(bx)-targetVal)/thisBinErr
            if nsigma < signifThreshold: 
                #reset central value to 0
                vlHists['Central'].SetBinContent(bx,targetVal)
                vlHists['Central'].SetBinError(bx,0.0)
                #reset up/down histograms 
                vlHists['Up'].SetBinContent(bx,targetVal)
                vlHists['Up'].SetBinError(bx,0.0)
                vlHists['Down'].SetBinContent(bx,targetVal)
                vlHists['Down'].SetBinError(bx,0.0)

    if debugLevel > 0:
        print "\nCorrections before MT and dPhi efficiencies:"
        vlHists['Central'].Print("all")

    vetoLeptonOutfile = rt.TFile("Razor"+regionNameReduced+"CrossCheck.root", "RECREATE")
    c = rt.TCanvas("cVetoLeptonCorrs", "c", 800, 600)
    #optionally, apply corrections to signal region histograms and derive Corrected/Uncorrected scale factors
    #also correct for MT cut efficiency for extrapolation into the signal region
    if histToCorrect is not None and not doDataOverMC:
        print "Correcting histogram",vlHists['Central']," and its up/down variants using MT cut efficiencies taken from",mtEfficiencyHist.GetName(),"and using dPhi cut efficiencies from",dPhiEfficiencyHist.GetName()
        #create MT up/down versions of the correction histogram
        vlHists['MTUp'] = vlHists['Central'].Clone()
        vlHists['MTDown'] = vlHists['Central'].Clone()
        #create dPhi up/down versions of the correction histogram
        vlHists['DPhiUp'] = vlHists['Central'].Clone()
        vlHists['DPhiDown'] = vlHists['Central'].Clone()
        for bx in range(1, vlHists['Central'].GetSize()+1):
            #divide all bin contents by the efficiency of the MT cut
            for n,h in vlHists.iteritems():
                if 'MT' not in n:
                    h.SetBinContent(bx, h.GetBinContent(bx)/mtEfficiencyHist.GetBinContent(bx))
                    h.SetBinError(bx, h.GetBinError(bx)/mtEfficiencyHist.GetBinContent(bx))
                elif 'Up' in n:
                    h.SetBinContent(bx, h.GetBinContent(bx)/(mtEfficiencyHist.GetBinContent(bx)+mtEfficiencyHist.GetBinError(bx)))
                    h.SetBinError(bx, h.GetBinError(bx)/(mtEfficiencyHist.GetBinContent(bx)+mtEfficiencyHist.GetBinError(bx)))
                elif 'Down' in n:
                    h.SetBinContent(bx, h.GetBinContent(bx)/(mtEfficiencyHist.GetBinContent(bx)-mtEfficiencyHist.GetBinError(bx)))
                    h.SetBinError(bx, h.GetBinError(bx)/(mtEfficiencyHist.GetBinContent(bx)-mtEfficiencyHist.GetBinError(bx)))
            #multiply all bin contents by the efficiency of the dPhi cut
                if 'DPhi' not in n:
                    h.SetBinContent(bx, h.GetBinContent(bx)*dPhiEfficiencyHist.GetBinContent(bx))
                    h.SetBinError(bx, h.GetBinError(bx)*dPhiEfficiencyHist.GetBinContent(bx))
                elif 'Up' in n:
                    h.SetBinContent(bx, h.GetBinContent(bx)*(dPhiEfficiencyHist.GetBinContent(bx)+dPhiEfficiencyHist.GetBinError(bx)))
                    h.SetBinError(bx, h.GetBinError(bx)*(dPhiEfficiencyHist.GetBinContent(bx)+dPhiEfficiencyHist.GetBinError(bx)))
                elif 'Down' in n:
                    h.SetBinContent(bx, h.GetBinContent(bx)*(dPhiEfficiencyHist.GetBinContent(bx)-dPhiEfficiencyHist.GetBinError(bx)))
                    h.SetBinError(bx, h.GetBinError(bx)*(dPhiEfficiencyHist.GetBinContent(bx)-dPhiEfficiencyHist.GetBinError(bx)))

        print "Correcting histogram",histToCorrect.GetName(),"using the additive corrections just derived."
        #set all histToCorrect errors to 0
        for bx in range(1, histToCorrect.GetSize()+1):
            histToCorrect.SetBinError(bx, 0.0)
        #corrected hist is histToCorrect + correction
        correctedHists = { n:histToCorrect.Clone(histToCorrect.GetName()+"VetoLepton"+n) for n in vlHists }
        for n in vlHists:
            correctedHists[n].Add(vlHists[n])
        #signal region scale factors are correctedHist / histToCorrect
        signalRegionSFs = { n:h.Clone(regionNameReduced+"ScaleFactors"+n.replace('Central','')) for n,h in correctedHists.iteritems() }

        #divide and write out
        for n,h in signalRegionSFs.iteritems():
            h.SetDirectory(0)
            h.Divide(histToCorrect)
            #zero any negative scale factors
            for bx in range(h.GetSize()+1):
                if useUncertainty:
                    h.SetBinContent(bx, 1 + h.GetBinError(bx))
                else:
                    h.SetBinContent(bx,max(0., h.GetBinContent(bx)))
            print "Writing histogram",h.GetName(),"to file"
            h.Write()

        if debugLevel > 0:
            print "\nCorrections:"
            vlHists['Central'].Print("all")
            print "\nBefore correction:"
            histToCorrect.Print("all")
            print "\nAfter correction:"
            correctedHists['Central'].Print("all")
            print "\nScale factors:"
            for n,h in signalRegionSFs.iteritems():
                print "\n"+n+":"
                h.Print("all")

    #otherwise, write corrections to file
    else:
        for n,h in vlHists.iteritems():
            #zero any negative scale factors
            for bx in range(h.GetSize()+1):
                h.SetBinContent(bx,max(0., h.GetBinContent(bx)))
            print "Writing histogram",h.GetName(),"to file"
            h.Write(regionNameReduced+"ScaleFactors"+n)

    #plot correction factors in 2D (not yet implemented for 1 or 3 dimensions)
    comment = "MC-Data"
    if doDataOverMC:
        comment = "Data/MC"
    if not isinstance(var, basestring) and len(var) == 2: 
        for n,h in vlHists.iteritems():
            commentstr=comment+" "+n+", "+regionName+" Control Region"
            plotting.draw2DHist(c, h, xtitle=var[0], ytitle=var[1], zmin=-200, zmax=200, printstr=regionNameReduced+"Correction"+n, lumistr=('%.1f' % (lumiData*1.0/1000))+" fb^{-1}", commentstr="", drawErrs=True, logz=False, numDigits=2, printdir=printdir)
        if histToCorrect is not None:
            for n,h in signalRegionSFs.iteritems():
                commentstr="Signal Region Scale Factors "+n+", "+regionName+" Control Region"
                plotting.draw2DHist(c, h, xtitle=var[0], ytitle=var[1], zmin=-200, zmax=200, printstr=regionNameReduced+"SignalRegionSFs"+n, lumistr=('%.1f' % (lumiData*1.0/1000))+" fb^{-1}", commentstr="", drawErrs=True, logz=False, numDigits=2, printdir=printdir)

        xbinLowEdges = []
        xbinUpEdges = []
        ybinLowEdges = []
        ybinUpEdges = []
        statUncerts = []
        sfs = []
        #for each bin, get values for all table columns
        for bx in range(1, vlHists['Central'].GetNbinsX()+1):
            for by in range(1, vlHists['Central'].GetNbinsY()+1):
                xbinLowEdges.append('%.0f' % (vlHists['Central'].GetXaxis().GetBinLowEdge(bx)))
                xbinUpEdges.append('%.0f' % (vlHists['Central'].GetXaxis().GetBinUpEdge(bx)))
                ybinLowEdges.append(str(vlHists['Central'].GetYaxis().GetBinLowEdge(by)))
                ybinUpEdges.append(str(vlHists['Central'].GetYaxis().GetBinUpEdge(by)))
                corrFactor = vlHists['Central'].GetBinContent(bx,by)
                sfs.append('%.3f' % (corrFactor))
                if corrFactor != 0:
                    statUncerts.append('%.1f\\%%' % (100*(vlHists['Central'].GetBinError(bx,by)/corrFactor)))
                else: 
                    statUncerts.append('--')
                for mcProcess in bgProcesses: 
                    dataYield = hists[dataName][var].GetBinContent(bx,by)
        xRanges = [low+'-'+high for (low, high) in zip(xbinLowEdges, xbinUpEdges)]
        yRanges = [low+'-'+high for (low, high) in zip(ybinLowEdges, ybinUpEdges)]
        headers=[var[0], var[1], comment, "Stat.\\ unc."]
        cols = [xRanges, yRanges, sfs, statUncerts]
        plotting.table_basic(headers, cols, caption="Difference between MC and data yields in veto lepton control region", printstr="corrFactorTable"+regionNameReduced, printdir=printdir)

    if histToCorrect is not None and not doDataOverMC:
        return signalRegionSFs['Central']
    else:
        return vlHists['Central']

#########################################
### PREPARE HISTOGRAMS FOR LIMIT SETTING
#########################################

def unrollAndStitchOne(hists, unrollBins, suffix=""):
    """
    One-physics-process version of unrollAndStitch.
    """
    # Note that each histogram is assumed to have its 
    # own set of unroll bins
    unrolled = [plotting.unroll2DHistograms(
            [h], xbins, cols, suffix=suffix+str(i))[0]
            for i, (h, (xbins, cols)) in enumerate(
                zip(hists, unrollBins))]
    return macro.stitch(unrolled)

def doQCDStatUnc(qcdHists, unrollBins, boxName):
    """
    The QCD prediction is made by applying transfer factors to the
    high dPhi event yields, which are counted inclusively in
    the number of b-tags.  Thus, QCD statistical uncertainties
    in different b-tag bins in the same box are correlated.
    This function returns a dictionary containing appropriately
    correlated QCD stat uncertainty histograms.
    """
    def getHistBaseName(box, bx, by):
        return 'QCD_stat{}QCD{}x{}y'.format(
                box, bx, by)
    qcdStatHists = {}
    notUsed = {}
    refHist = qcdHists[0] # use this to count bins
    for refbx in range(1, refHist.GetNbinsX()+1):
        for refby in range(1, refHist.GetNbinsY()+1):
            histBaseName = getHistBaseName(boxName, refbx, refby)
            notUsed[histBaseName] = True
            for updown in ['Up', 'Down']:
                histName = histBaseName + updown
                qcdStatHists[histName] = unrollAndStitchOne(
                    qcdHists, unrollBins, 
                    suffix='qcdstat{}x{}y'.format(refbx, refby))
                qcdStatHists[histName].Reset()
                # Now we have an empty unrolled hist
    ix = 0 # cumulative bin number
    for hist, (xbins, cols) in zip(qcdHists, unrollBins):
        for i, xdo in enumerate(xbins[:-1]):
            xup = xbins[i+1]
            for j, ydo in enumerate(cols[i][:-1]):
                yup = cols[i][j+1]
                ix += 1
                for xbin in range(1, hist.GetNbinsX()+1):
                    for ybin in range(1, hist.GetNbinsY()+1):
                        centerX = hist.GetXaxis().GetBinCenter(xbin)
                        centerY = hist.GetYaxis().GetBinCenter(ybin)
                        if not (centerX > xdo and centerX < xup
                                and centerY > ydo and centerY < yup):
                            continue
                        histBaseName = getHistBaseName(
                                boxName, xbin, ybin)
                        # If reference bin is current bin,
                        # update with up/down variation.
                        # Otherwise, use the nominal counts
                        for name, h in qcdStatHists.iteritems():
                            curContent = h.GetBinContent(ix)
                            newContent = (curContent 
                                    + hist.GetBinContent(xbin, ybin))
                            if name.startswith(histBaseName):
                                if histBaseName in notUsed:
                                    del notUsed[histBaseName]
                                err = hist.GetBinError(xbin, ybin)
                                if name.endswith('Up'):
                                    newContent += err
                                else:
                                    newContent -= err
                                    newContent = max(newContent, 0.)
                            h.SetBinContent(ix, newContent)
    for name in notUsed:
        print "Not using histogram {}".format(name)
        del qcdStatHists[name+'Up']
        del qcdStatHists[name+'Down']
    return qcdStatHists

def getMaxBtags(boxName):
    maxBtags = 3
    if boxName in ['DiJet','LeptonJet']:
        maxBtags = 2
    return maxBtags

def unrollHistsForDataCard(hists, samples=[], noSys=False,
        var=('MR','Rsq'), debugLevel=0, unrollBins=None):
    """
    Unrolls each histogram.
    hists: list of dictionaries, 
        each the output of, eg, makeControlSampleHists.
    samples: list of MC sample names
    Returns: unrolled MC, data, fit, and shape histograms.
        Data and fit histograms are returned as lists,
        while MC and shape uncertainty hists are returned
        as dictionaries of lists, indexed by sample name.
    """
    unrolledMC = {s:[] for s in samples}
    unrolledData = []
    unrolledFit = []
    unrolledShapeHists = {s:{} for s in samples}
    for i, curHists in enumerate(hists):
        unrollRows = unrollCols = None
        if unrollBins is not None:
            unrollRows, unrollCols = unrollBins[i]

        #get shape histograms
        if 'Sys' in curHists: 
            shapeHists = curHists['Sys']
            del curHists['Sys']
        elif not noSys:
            raise KeyError("No shape unc. histograms were found!")
        if noSys: 
            shapeHists = {s:{} for s in samples}

        #unroll each MC histogram
        unrolledMCs = plotting.unroll2DHistograms(
                [curHists[s][var] for s in samples], 
                unrollRows, unrollCols, suffix="_%d"%(i))
        for n,s in enumerate(samples):
            unrolledMC[s].append(unrolledMCs[n])
            for shape in shapeHists[s]:
                #omit some uncertainties
                if 'singletopnorm' in shape or 'othernorm' in shape: 
                    continue
                if shape not in unrolledShapeHists[s]:
                    unrolledShapeHists[s][shape] = []
                unrolledShapeHists[s][shape].append(
                        plotting.unroll2DHistograms(
                            [shapeHists[s][shape][var]], 
                            unrollRows, unrollCols, suffix="_%d"%(i))[0])
        #unroll data and fit
        if 'Data' in curHists:
            unrolledData.append(plotting.unroll2DHistograms(
                [curHists['Data'][var]], unrollRows, unrollCols, 
                suffix="_%d"%(i))[0])
        if 'Fit' in curHists:
            unrolledFit.append(plotting.unroll2DHistograms(
                [curHists['Fit'][var]], unrollRows, unrollCols, 
                suffix="_%d"%(i))[0])

    return unrolledMC, unrolledData, unrolledFit, unrolledShapeHists

def stitchHistsForDataCard(unrolledMC, unrolledData, 
    unrolledShapeHists, samples=[]): 
    #piece together histograms from different b-tag bins
    histsForDataCard = {}
    for s in samples:
        histsForDataCard[s] = macro.stitch(unrolledMC[s])
        for shape in unrolledShapeHists[s]:
            #protect against empty histograms (for QCD)
            isEmpty = True
            needsTempPatch = False
            for hist in unrolledShapeHists[s][shape]:
                if hist.Integral() > 0:
                    isEmpty = False
                    break
            # qcdnormDown histogram is occasionally all zeros due to large uncertainty
            if isEmpty and shape != 'qcdnormDown': 
                print "Warning: empty shape histograms for",s,shape
                if s == 'SingleTop' and (shape == 'renscaleUp' or shape == 'facrenscaleUp'):
                    print "Ignoring empty shape histogram for now -- PLEASE FIX THOUGH"
                    needsTempPatch = True
                else:
                    continue
            #create histogram
            histsForDataCard[s+'_'+shape] = macro.stitch(
                    unrolledShapeHists[s][shape])
            if needsTempPatch:
                print "WARNING: USING DOWN HISTOGRAM INSTEAD OF UP, DUE TO ISSUE WITH WEIGHTS NORMALIZATION."
                histsForDataCard[s+'_'+shape] = macro.stitch(
                        unrolledShapeHists[s][shape.replace('Up','Down')])
    if len(unrolledData) > 0:
        histsForDataCard['data_obs'] = macro.stitch(unrolledData)
    return histsForDataCard

def addMCVsFitUncertainty(histsForDataCard, unrolledFit, samples=[]):
    fitTotal = macro.stitch(unrolledFit)

    #make total MC histogram to compare with fit
    mcTotal = histsForDataCard[samples[0]].Clone('mcTotal')
    mcTotal.Reset()
    for s in samples:
        mcTotal.Add(histsForDataCard[s])
        histsForDataCard[s+'_fitmccrosscheckUp'] = histsForDataCard[s].Clone()
        histsForDataCard[s+'_fitmccrosscheckDown'] = histsForDataCard[s].Clone()
        histsForDataCard[s+'_fitmccrosscheckUp'].SetName(s+'_fitmccrosscheckUp')
        histsForDataCard[s+'_fitmccrosscheckDown'].SetName(s+'_fitmccrosscheckDown')
        histsForDataCard[s+'_fitmccrosscheckUp'].SetTitle(s+'_fitmccrosscheckUp')
        histsForDataCard[s+'_fitmccrosscheckDown'].SetTitle(s+'_fitmccrosscheckDown')
    
    #get fit/MC and propagate to individual MC processes
    for bx in range(1, histsForDataCard[samples[0]].GetNbinsX()+1):
        if mcTotal.GetBinContent(bx) > 0:
            fitOverMC = fitTotal.GetBinContent(bx) / mcTotal.GetBinContent(bx)

            for s in samples:
                histsForDataCard[s+'_fitmccrosscheckUp'].SetBinContent(bx, 
                        histsForDataCard[s].GetBinContent(bx) * fitOverMC)
                histsForDataCard[s+'_fitmccrosscheckDown'].SetBinContent(bx, 
                        histsForDataCard[s].GetBinContent(bx) / fitOverMC)

def stitchEmptyBinErrs(emptyBinErrs):
    """
    Converts a list of empty bin dictionaries into a single one 
    in which the bin numbers are indexed according to the
    corresponding stitched-together histogram 
    """
    combinedDict = {}
    globalCounters = {}
    for curBinErrs in emptyBinErrs:
        for proc, procBinErrs in curBinErrs.iteritems():
            if proc not in globalCounters:
                globalCounters[proc] = 0
            if proc not in combinedDict:
                combinedDict[proc] = {}
            localCounter = 0
            for ibin in sorted(procBinErrs):
                globalCounters[proc] += 1
                localCounter += 1
                if ibin != localCounter:
                    raise ValueError(
                            "Missing bin {} in empty bin dictionary".format(localCounter))
                combinedDict[proc][globalCounters[proc]] = procBinErrs[ibin]
    return combinedDict

def addStatUncertainty(histsForDataCard, boxName, samples, emptyBinErrs=None):
    for sample in samples:
        if sample == 'QCD':
            continue # QCD stat errors are handled by doQCDStatUnc()
        print "Adding stat uncertainty for {}".format(sample)
        histsForDataCard[sample+'_stat'+boxName+sample+'Up'] = (
                histsForDataCard[sample].Clone())
        histsForDataCard[sample+'_stat'+boxName+sample+'Down'] = (
                histsForDataCard[sample].Clone())
        for bx in range(1, histsForDataCard[sample].GetNbinsX()+1):
            histsForDataCard[
                sample+'_stat'+boxName+sample+'Up'].SetBinContent(
                    bx, histsForDataCard[sample].GetBinContent(bx) 
                    + histsForDataCard[sample].GetBinError(bx))
            histsForDataCard[
                sample+'_stat'+boxName+sample+'Down'].SetBinContent(
                    bx, histsForDataCard[sample].GetBinContent(bx) 
                    - histsForDataCard[sample].GetBinError(bx))
            if emptyBinErrs is not None:
                # Increase up uncertainty equal to twice the error assigned
                # by computeEmptyBinErrs, because that function assumes a 
                # symmetrical error while we want to add the whole thing to the 
                # 'up' uncertainty
                histsForDataCard[sample+'_stat'+boxName+sample+'Up'].SetBinContent(
                        bx, histsForDataCard[
                            sample+'_stat'+boxName+sample+'Up'].GetBinContent(
                                bx) + 2 * emptyBinErrs[sample][bx])
                print "Increasing error on bin {} by {}".format(
                        bx, emptyBinErrs[sample][bx])
        histsForDataCard[sample+'_stat'+boxName+sample+'Up'].SetName(
            sample+'_stat'+boxName+sample+'Up')
        histsForDataCard[sample+'_stat'+boxName+sample+'Down'].SetName(
            sample+'_stat'+boxName+sample+'Down')
        histsForDataCard[sample+'_stat'+boxName+sample+'Up'].SetTitle(
            sample+'_stat'+boxName+sample+'Up')
        histsForDataCard[sample+'_stat'+boxName+sample+'Down'].SetTitle(
            sample+'_stat'+boxName+sample+'Down')

def unrollAndStitch(hists, boxName, samples=[], var=('MR','Rsq'),
        debugLevel=0, unrollBins=None, noSys=False, addStatUnc=True,
        addMCVsFit=False, emptyBinErrs=None):
    """
    Unrolls each histogram and pieces together the different
    b-tag bins.
    """
    unrolledMC, unrolledData, unrolledFit, unrolledShapeHists = (
            unrollHistsForDataCard(hists, samples, noSys,
            var, debugLevel, unrollBins))

    histsForDataCard = stitchHistsForDataCard(unrolledMC,
            unrolledData, unrolledShapeHists, samples)

    if emptyBinErrs is not None:
        emptyBinErrs = stitchEmptyBinErrs(emptyBinErrs)
    if addMCVsFit:
        addMCVsFitUncertainty(histsForDataCard, unrolledFit, samples)
    if addStatUnc:
        addStatUncertainty(histsForDataCard, boxName, samples,
                emptyBinErrs=emptyBinErrs)
        if 'QCD' in samples:
            qcdStatHists = doQCDStatUnc(
                    [h['QCD'][('MR', 'Rsq')] for h in hists], 
                    unrollBins, boxName)
            histsForDataCard.update(qcdStatHists)

    #set names on histograms
    if 'data_obs' in histsForDataCard:
        histsForDataCard['data_obs'].SetName('data_obs')
        histsForDataCard['data_obs'].SetTitle('data_obs')
    for sample in samples:
        histsForDataCard[sample].SetName(sample)
        histsForDataCard[sample].SetTitle(sample)
        for shape in unrolledShapeHists[sample]:
            if sample+'_'+shape in histsForDataCard:
                histsForDataCard[sample+'_'+shape].SetName(
                    sample+'_'+shape)
                histsForDataCard[sample+'_'+shape].SetTitle(
                    sample+'_'+shape)

    return histsForDataCard
         
def unrollAndStitchFromFiles(boxName, samples=[], inDir=".", outDir=".", 
        var=('MR','Rsq'), debugLevel=0, unrollBins=None, export=True, 
        noSys=False, addStatUnc=True, addMCVsFit=False, doEmptyBinErrs=False, weightHists=None):
    """
    Loads the output of makeControlSampleHists, unrolls each 
    histogram, and pieces together the different b-tag bins 
    to get the histograms used for limit setting.
    """

    maxBtags = getMaxBtags(boxName)
    filenames = [inDir+"/razorHistograms"+boxName+str(b)+"BFineGrained.root" 
            for b in range(maxBtags+1)]

    hists = []
    emptyBinErrs = None
    for ifile, f in enumerate(filenames):
        hists.append(macro.importHists(f, debugLevel))
        if doEmptyBinErrs:
            if emptyBinErrs is None:
                emptyBinErrs = []
            emptyBinErrs.append(macro.computeEmptyBinErrs(hists[ifile],
                unrollBins[ifile], aggregate=False))
        if 'QCD' in samples:
            makeQCDBtagSysHists(hists[ifile], weightHists, boxName, 
                    ifile)

    histsForDataCard = unrollAndStitch(hists, boxName, samples, 
            var, debugLevel, unrollBins, noSys, addStatUnc, 
            addMCVsFit, emptyBinErrs=emptyBinErrs)

    if export:
        macro.exportHists(histsForDataCard, 
            outFileName='razorBackgroundHists'+boxName+'.root', 
            outDir=outDir, useDirectoryStructure=False, delete=False, 
            debugLevel=debugLevel)

    return histsForDataCard

def printRazorEvidenceTable(evidenceHist, obsEvidenceHist, 
        unrolledMCTotal, unrolledData, out_f):
    nbins = evidenceHist.GetNbinsX()
    evidences = np.zeros(nbins)
    obsEvidences = np.zeros(nbins)
    mcs  = np.zeros(nbins)
    mcerrs  = np.zeros(nbins)
    datas = np.zeros(nbins)
    for ix in range(nbins):
        evidences[ix] = max(0, evidenceHist.GetBinContent(ix+1))
        obsEvidences[ix] = obsEvidenceHist.GetBinContent(ix+1)
        mcs[ix] = unrolledMCTotal.GetBinContent(ix+1)
        mcerrs[ix] = unrolledMCTotal.GetBinError(ix+1)
        datas[ix] = unrolledData.GetBinContent(ix+1)
    indices = np.argsort(-evidences)
    out_f.write("\n")
    for i, (ev, mc, err, obs, obsev) in enumerate(np.dstack((evidences[indices], 
        mcs[indices], mcerrs[indices], datas[indices], obsEvidences[indices]))[0]):
        if ev == 0: continue
        out_f.write("{}: Evidence {:.3f}, Predicted {:.2f} +/- {:.2f}, Observed {:.0f}, Observed Evidence {:.3f}\n".format(
                i, ev, mc, err, obs, obsev))
    out_f.write("Totals: Evidence {:.3f}, Predicted {:.2f}, Observed {:.0f}, Observed Evidence {:.3f}\n".format(
        sum(evidences), sum(mcs), sum(datas), sum(obsEvidences)))

def makeRazorBinEvidencePlots(boxName, samples, inDir='.', signalHist=None, outDir='.', unrollBins=None, 
        zmin=1e-3, unblind=False, debugLevel=0):
    nbMax = 3
    if boxName == 'DiJet':
        nbMax = 2
    filenames = [inDir+"/razorHistograms"+boxName+str(b)+"BFineGrained.root" for b in range(nbMax+1)]
    c = rt.TCanvas("d", "d", 800, 600)
    
    signalBin = 0
    if unblind:
        out_txt = open(outDir+"/yields.txt", 'w')
    for i,f in enumerate(filenames):
        unrollRows = unrollBins[i][0]
        unrollCols = unrollBins[i][1]
        unrolledMCTotal = macro.makeRazorMCTotalUnrolledHist(f, samples, unrollRows, unrollCols, debugLevel)
        
        #loop over bins
        evidenceHist = unrolledMCTotal.Clone("evidence"+str(i))
        evidenceHist.Reset()
        if unblind:
            unrolledData = macro.makeRazorMCTotalUnrolledHist(f, ['Data'], unrollRows, unrollCols, debugLevel)
            obsEvidenceHist = unrolledData.Clone("obsevidence"+str(i))
            obsEvidenceHist.Reset()
        for bx in range(1, evidenceHist.GetNbinsX()+1):
            signalBin += 1
            s = signalHist.GetBinContent(signalBin)
            b = unrolledMCTotal.GetBinContent(bx)
            evidence = macro.getBinEvidence(b, b, s) #expected contribution to the likelihood
            if unblind:
                obs = unrolledData.GetBinContent(bx)
                obsEvidenceHist.SetBinContent(bx, 
                        macro.getBinEvidence(obs, b, s))
            if evidence > zmin:
                evidenceHist.SetBinContent(bx, evidence)
            else:
                evidenceHist.SetBinContent(bx, -999)
            evidenceHist.SetBinError(bx, 0)
        if unblind:
            printRazorEvidenceTable(evidenceHist, obsEvidenceHist,
                    unrolledMCTotal, unrolledData, out_txt)

        #plot the evidence
        plotting.plotEvidenceHist(c, evidenceHist, printstr="evidence"+boxName+str(i)+"BTag"+signalHist.GetName(), obs=obsEvidenceHist,
                printdir=outDir, unrollBins=unrollBins[i], zmin=zmin)

        evidenceHist2D = macro.makeTH2PolyFromColumns(evidenceHist.GetName()+"2D", evidenceHist.GetName()+"2D", 
                unrollBins[i][0], unrollBins[i][1])
        macro.fillTH2PolyFromTH1(evidenceHist, evidenceHist2D)
        plotting.draw2DHist(c, evidenceHist2D, xtitle="MR", ytitle="Rsq", zmin=zmin, zmax=3.84,
                printstr="evidence2D"+boxName+str(i)+"BTag"+signalHist.GetName(), lumistr="2.3 fb^{-1}", dotext=True, palette=56,
                numDigits=3, textSize=1.0, printdir=outDir)

def makeRazorMCTotalPlusSignalPlot(boxName, samples, inDir='.', 
        signalHist=None, outDir='.', unrollBins=None, unblind=False,
        signalString="Signal", modelName="Signal", debugLevel=0):
    nbmax = 3
    if boxName == 'DiJet':
        nbmax = 2
    filenames = [inDir+"/razorHistograms"+boxName+str(b)+"BTag.root" for b in range(nbmax+1)]
    hists = macro.makeRazorMCTotalUnrolledHists(
            boxName, samples, inDir, unrollBins, debugLevel, doData=unblind)
    if unblind:
        mcTotalHists, dataHists = hists
    else:
        mcTotalHists = hists

    c = rt.TCanvas("d", "d", 800, 600)
    rt.SetOwnership(c, False)
    
    unrolledSignals = macro.splitByUnrollBins(signalHist, unrollBins)
    for i,f in enumerate(filenames):
        unrolledRowsCols = unrollBins[i]
        unrolledMCTotal = mcTotalHists[i]
        unrolledSignal = unrolledSignals[i]
        if unblind:
            unrolledData = dataHists[i]
        else:
            unrolledData = None
        
        comment = '#it{'+boxName+' '+str(i)+' b-tag}'
        plotting.plot_SUS15004_MCTotalWithSignal(c, mcTotalUnrolled=unrolledMCTotal, signalUnrolled=unrolledSignal, 
                printstr="MRRsq"+boxName+str(i)+"BTagUnrolled"+modelName, lumistr="2.3 fb^{-1}", 
                commentstr=comment, unrollBins=unrolledRowsCols, printdir=outDir, signalString=signalString, dataUnrolled=unrolledData,
                ratiomin=0.5, ratiomax=1.5)
