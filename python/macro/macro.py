import os
import ROOT as rt
import copy
import math
from array import array
from ast import literal_eval

#local imports
from PlotFit import setFFColors
from RunCombine import exec_me
from razorWeights import applyMTUncertainty1D, applyMTUncertainty2D, getSFHistNameForErrorOpt
from printJson import walk
import plotting

def walkDictionary(d, path=""):
    """
    Recursively iterate over a nested dictionary.
    Returns ordered pairs (p, x) where x is a value stored at some level of the dictionary and 
    p is the concatenation of all keys (joined with '/') needed to reach x from the top level.
    """
    for key,val in d.iteritems():
        #add key to path
        newpath = path+'/'+str(key)
        if newpath.startswith('/'): 
            newpath = newpath[1:]
        #yield or recurse
        if isinstance(val, dict):
            for x in walkDictionary(val, newpath):
                yield x
        else:
            yield (newpath, val)

def exportHists(hists, outFileName='hists.root', outDir='.', useDirectoryStructure=True, varName=None, delete=True, debugLevel=0):
    """
    Writes histograms from the given dictionary to a ROOT file.  
    If useDirectoryStructure is True, histograms will be put in folders according to the structure of the dictionary.
    If varName is provided, only histograms of the provided variable will be saved.
    """

    print "Storing histograms in file",outFileName
    outFile = rt.TFile(outDir+'/'+outFileName, 'recreate')
    for pair in walkDictionary(hists): #get (path_to_histogram, histogram) pairs
        if useDirectoryStructure:
            path = pair[0]
            if varName is not None:
                if str(varName) not in path: #skip histograms of other variables
                    continue
                else: #remove variable name from path
                    path = path.replace(str(varName),'').replace('//','/')
            tdir = outFile.GetDirectory(path)
            if tdir==None:
                if debugLevel > 0:
                    print "Making directory",path
                outFile.mkdir(path)
                tdir = outFile.GetDirectory(path)
                tdir.cd()
                if debugLevel > 0:
                    print "Writing histogram",pair[1].GetName(),"to directory",path
                pair[1].Write()
                if delete: pair[1].Delete()
        else:
            if debugLevel > 0:
                print "Writing histogram",pair[1].GetName()
            pair[1].Write()
            if delete: pair[1].Delete()
    outFile.Close()

def importHists(inFileName='hists.root', debugLevel=0):
    """
    Creates a dictionary according to the directory structure of the input file, 
    and populates the dictionary with the histograms in the file.
    """

    hists = {}
    print "\nGetting histograms from file",inFileName
    inFile = rt.TFile.Open(inFileName)
    for dirpath, dirnames, filenames, tdirectory in walk(inFile):
        if len(filenames) > 0: #there are objects to retrieve
            dirInFile = dirpath.split(':')[-1].split('/')[1:] #get path to histogram as a list of subdirectories
            currentLayer = hists #keep track of current layer of the dictionary
            storedHist = False
            if dirInFile != ['']: #dirInFile will be [''] if there are no subdirectories
                for i,subdir in enumerate(dirInFile): #descend along path to histogram
                    if subdir not in currentLayer: 
                        if len(filenames) == 1 and i == len(dirInFile)-1:
                            #if there is only one histogram, store it now
                            try:
                                lastSubdir = literal_eval(subdir)
                                if debugLevel > 0:
                                    print "Converted",subdir,"to tuple"
                            except ValueError:
                                lastSubdir = subdir
                            currentLayer[lastSubdir] = tdirectory.Get(filenames[0])
                            currentLayer[lastSubdir].SetDirectory(0)
                            if debugLevel > 0: 
                                print "Retrieved histogram",filenames[0]
                            storedHist = True
                            break
                        currentLayer[subdir] = {} #create dict layer if needed
                    currentLayer = currentLayer[subdir] #now we go one level lower and repeat
            #now currentLayer is the lowest-level dictionary, where the remaining histograms should be inserted
            if not storedHist:
                for tkey in filenames:
                    currentLayer[tkey] = tdirectory.Get(tkey)
                    currentLayer[tkey].SetDirectory(0)
                    if debugLevel > 0: 
                        print "Retrieved histogram",tkey
    inFile.Close()
    return hists

def stitch(th1s = []):
    """Join the input histograms end-to-end to form one long 1D histogram.
    Bins in the output histogram have width 1 and the x-axis starts at 0."""
    if len(th1s) == 0:
        return
    #get total number of bins across all inputs
    totalNBins = sum([th1.GetNbinsX() for th1 in th1s])
    #concatenate hist names
    name = 'and'.join([th1.GetName() for th1 in th1s])
    #make output hist
    out = rt.TH1F(name+"stitched", name, totalNBins, 0, totalNBins)
    out.SetDirectory(0)
    bn = 1
    for i,th1 in enumerate(th1s):
        for bx in range(1, th1.GetNbinsX()+1):
            out.SetBinContent(bn, th1.GetBinContent(bx))
            out.SetBinError(bn, th1.GetBinError(bx))
            bn += 1
    return out

def getBinBoundariesFromColumns(xbins, cols):
    """
    Get lower and upper bin edges in the x and y directions
    (Useful for making yield tables)
    """
    xLow = []
    xHigh = []
    yLow = []
    yHigh = []
    for i in range(len(xbins)-1):
        for j in range(len(cols[i])-1):
            xLow.append(xbins[i])
            xHigh.append(xbins[i+1])
            yLow.append(cols[i][j])
            yHigh.append(cols[i][j+1])
    return xLow, xHigh, yLow, yHigh

def makeTH2PolyFromColumns(name, title, xbins, cols):
    """
    Makes a TH2Poly histogram with rectangular bins. 
    xbins: list of bin edges in the x-direction.
    cols: 2D list indicating bin edges in the y-direction for each column
    """
    #construct TH2Poly
    minY = min([col[0] for col in cols])
    poly = rt.TH2Poly(name, title, xbins[0], xbins[-1], minY, cols[0][-1])
    #add bins in each column
    for i in range(len(xbins)-1):
        for j in range(len(cols[i])-1):
            poly.AddBin( xbins[i], cols[i][j], xbins[i+1], cols[i][j+1] )

    poly.SetDirectory(0)
    poly.Sumw2()
    return poly

def fillTH2PolyFromTH1(th1, poly):
    """Fills the given TH2Poly histogram using the bin contents of the given TH1, with TH2Poly bins being filled in 
    the order they were created."""
    for bx in range(1, th1.GetNbinsX()+1):
        poly.SetBinContent(bx, th1.GetBinContent(bx))
        poly.SetBinError(bx-1, th1.GetBinError(bx))

def fillTH2PolyFromTH2(th2, poly):
    """Fills the given TH2Poly histogram using the bin contents of the given TH2.  If multiple bins in the TH2 map to the same bin in the TH2Poly, the bin contents will be added and the errors will be summed in quadrature."""
    for bx in range(1, th2.GetNbinsX()+1):
        for by in range(1, th2.GetNbinsY()+1):
            #get bin number in TH2Poly
            centerX = th2.GetXaxis().GetBinCenter(bx)
            centerY = th2.GetYaxis().GetBinCenter(by)
            bn = poly.FindBin(centerX, centerY)
            poly.SetBinContent(bn, poly.GetBinContent(bn) + th2.GetBinContent(bx,by))
            #NOTE: TH2Poly has a bug that causes SetBinError(x) to set the error of bin x+1. Beware!!!
            poly.SetBinError(bn-1, ( poly.GetBinError(bn)**2 + th2.GetBinError(bx,by)**2 )**(0.5) )

def makeTH2PolyRatioHist(num, denom, xbins, cols):
    """Makes a TH2Poly using makeTH2PolyFromColumns on the provided numerator and denominator histograms; 
    then divides the two histograms and propagates uncertainty assuming uncorrelated num. and denom."""

    #make TH2Poly with custom binning
    numPoly = makeTH2PolyFromColumns( num.GetName()+"Poly", num.GetTitle()+"Poly", xbins, cols)
    denomPoly = makeTH2PolyFromColumns( denom.GetName()+"Poly", denom.GetTitle()+"Poly", xbins, cols)

    #populate with numerator and denominator histogram bin contents
    fillTH2PolyFromTH2(num, numPoly) 
    fillTH2PolyFromTH2(denom, denomPoly)

    #need to divide by hand
    for bn in range(1, numPoly.GetNumberOfBins()+1):
        if denomPoly.GetBinContent(bn) != 0:
            numPoly.SetBinContent(bn, numPoly.GetBinContent(bn)/denomPoly.GetBinContent(bn))
            #NOTE: TH2Poly has a bug that causes SetBinError(x) to set the error of bin x+1. Beware!!!
            numPoly.SetBinError(bn-1, ( ( numPoly.GetBinError(bn)/denomPoly.GetBinContent(bn) )**2 + ( denomPoly.GetBinError(bn)*numPoly.GetBinContent(bn)/(denomPoly.GetBinContent(bn))**2 )**2 )**(0.5) )

    return numPoly

def blindHistograms(histList, blindBins):
    """blindBins should be a list of ordered pairs corresponding to bin coordinates.
    Sets blinded bin contents to -999."""
    for hist in histList:
        for (x,y) in blindBins:
            hist.SetBinContent(x,y,-999)

def setupHistograms(regionName, inputs, samples, bins, titles, shapeErrors, dataName):
    """Creates dictionary of histograms with specified binning.

    regionName: used only in histogram names and titles
    inputs: dictionary of process:filename pairs
    samples: list of processes for which histograms should be made
    bins: dictionary.  key = name of observable, value = list of bins for that observable
    titles: dictionary.  key = name of observable, value = title for histograms of that quantity 
    shapeErrors: list of shape uncertainties to apply
    """

    hists = {name:{} for name in inputs}
    shapeHists = {name:{} for name in inputs}
    if inputs is None: return
    for name in inputs:
        for var in bins:
            if isinstance(var, basestring): 
                #1D histograms
                if var in titles: title=titles[var]
                else: title = var
                hists[name][var] = rt.TH1F(regionName+var+name, title+';'+title+';', len(bins[var])-1, array('d',bins[var]))
                #add up/down histograms for each systematic uncertainty
                if samples is not None and name in samples:
                    for shape in shapeErrors:
                        if not isinstance(shape,basestring): #tuple (shape, [list of processes])
                            if name not in shape[1]: continue
                            curShape = shape[0]
                        else:
                            curShape = shape
                        if curShape+"Down" not in shapeHists[name]: shapeHists[name][curShape+"Down"] = {}
                        if curShape+"Up" not in shapeHists[name]: shapeHists[name][curShape+"Up"] = {}
                        shapeHists[name][curShape+"Down"][var] = hists[name][var].Clone(hists[name][var].GetName()+curShape+"Down")
                        shapeHists[name][curShape+"Up"][var] = hists[name][var].Clone(hists[name][var].GetName()+curShape+"Up")
            elif len(var) == 2: 
                #2D histograms
                title = [titles[v] if v in titles else v for v in var]
                hists[name][var] = rt.TH2F(regionName+var[0]+var[1]+name, ';'+title[0]+';'+title[1], len(bins[var[0]])-1, array('d',bins[var[0]]), len(bins[var[1]])-1, array('d',bins[var[1]]))
                if samples is not None and name in samples:
                    for shape in shapeErrors:
                        if not isinstance(shape,basestring): #tuple (shape, [list of processes])
                            if name not in shape[1]: continue
                            curShape = shape[0]
                        else:
                            curShape = shape
                        if curShape+"Down" not in shapeHists[name]: shapeHists[name][curShape+"Down"] = {}
                        if curShape+"Up" not in shapeHists[name]: shapeHists[name][curShape+"Up"] = {}
                        shapeHists[name][curShape+"Down"][var] = hists[name][var].Clone(hists[name][var].GetName()+curShape+"Down")
                        shapeHists[name][curShape+"Up"][var] = hists[name][var].Clone(hists[name][var].GetName()+curShape+"Up")
            elif len(var) == 3:
                #3D histograms
                title = [titles[v] if v in titles else v for v in var]
                hists[name][var] = rt.TH3F(regionName+var[0]+var[1]+var[2]+name, ';'+title[0]+';'+title[1]+';'+title[2], len(bins[var[0]])-1, array('d',bins[var[0]]), len(bins[var[1]])-1, array('d',bins[var[1]]), len(bins[var[2]])-1, array('d',bins[var[2]]))
                if samples is not None and name in samples:
                    for shape in shapeErrors:
                        if not isinstance(shape,basestring): #tuple (shape, [list of processes])
                            if name not in shape[1]: continue
                            curShape = shape[0]
                        else:
                            curShape = shape
                        if curShape+"Down" not in shapeHists[name]: shapeHists[name][curShape+"Down"] = {}
                        if curShape+"Up" not in shapeHists[name]: shapeHists[name][curShape+"Up"] = {}
                        shapeHists[name][curShape+"Down"][var] = hists[name][var].Clone(hists[name][var].GetName()+curShape+"Down")
                        shapeHists[name][curShape+"Up"][var] = hists[name][var].Clone(hists[name][var].GetName()+curShape+"Up")
        for var in hists[name]: 
            hists[name][var].Sumw2()
            hists[name][var].SetDirectory(0)
            if samples is not None and name in samples:
                for shape in shapeErrors:
                    if not isinstance(shape,basestring):
                        if name not in shape[1]: continue
                        curShape = shape[0]
                    else:
                        curShape = shape
                    shapeHists[name][curShape+"Down"][var].Sumw2()
                    shapeHists[name][curShape+"Up"][var].Sumw2()
                    shapeHists[name][curShape+"Down"][var].SetDirectory(0)
                    shapeHists[name][curShape+"Up"][var].SetDirectory(0)
    #deal correctly with Poisson errors on data
    if dataName in hists:
        for var in hists[dataName]:
            hists[dataName][var].SetBinErrorOption(rt.TH1.kPoisson)
    return hists,shapeHists

def propagateShapeSystematics(hists, samples, varList, shapeHists, shapeErrors, miscErrors=[], boxName="", debugLevel=0):
    """For each bin of the central histogram, add the appropriate uncertainties in quadrature with the statistical uncertainty.
    List of arguments is similar to razorMacros.makeControlSampleHists
    """
    for var in varList:
        for name in samples:
            for shape in shapeErrors:
                if not isinstance(shape,basestring): #tuple (shape, [list of processes])
                    if name not in shape[1]: continue
                    curShape = shape[0]
                else:
                    curShape = shape
                if debugLevel > 0: print "Adding",curShape,"uncertainty in quadrature with",name,"errors for",var
                #loop over histogram bins
                for bx in range(hists[name][var].GetSize()+1):
                    #use difference between Up and Down histograms as uncertainty
                    sysErr = abs(shapeHists[name][curShape+'Up'][var].GetBinContent(bx) - shapeHists[name][curShape+'Down'][var].GetBinContent(bx))/2.0
                    #add in quadrature with existing error
                    oldErr = hists[name][var].GetBinError(bx)
                    if not rt.TMath.IsNaN(sysErr):
                        hists[name][var].SetBinError(bx, (oldErr**2 + sysErr**2)**(0.5))
                    if debugLevel > 0 and sysErr > 0: print curShape,": Error on bin ",bx,"increases from",oldErr,"to",hists[name][var].GetBinError(bx),"after adding",sysErr,"in quadrature"
            for source in miscErrors:
                #MT uncertainty (deprecated)
                if source.lower() == "mt" and var == "MR":
                    if isinstance(var, basestring): #1D histogram
                        applyMTUncertainty1D(hists[name][var], process=name+"_"+boxName, debugLevel=debugLevel)
                    else: #2D histogram
                        applyMTUncertainty2D(hists[name][var], process=name+"_"+boxName, debugLevel=debugLevel)

def subtractBkgsInData(process, hists={}, dataName="Data", debugLevel=0):
    """
    Subtracts all MC backgrounds (except 'process') from each data yield histogram.
    The data histogram is modified in-place in the hists dictionary
    """
    #warn if the needed histograms are not found
    if dataName not in hists:
        print "Error in subtractBkgsInData: target data histogram (",dataName,") was not found!"
        return
    bgProcesses = [mcProcess for mcProcess in hists if mcProcess != process and mcProcess != dataName and mcProcess != "Fit"] #get list of non-data, non-fit, non-signal samples
    if debugLevel > 0:
        print "Will isolate",process,"process by subtracting these backgrounds from the data yield:"
        print bgProcesses
    #subtract backgrounds in each data histogram
    for var in hists[dataName]:
        for p in bgProcesses:
            #make sure relevant background histogram exists
            if var not in hists[p]:
                print "Warning in subtractBkgsInData: could not find",var," in hists[",p,"]."
                continue
            #subtract it
            if debugLevel > 0: 
                print "Subtracting",p,"from",dataName,"distribution for",var
            hists[dataName][var].Add(hists[p][var], -1) 
        #zero any negative bins
        for bx in range(1, hists[dataName][var].GetSize()+1):
            if hists[dataName][var].GetBinContent(bx) < 0:
                print "Zeroing negative prediction",hists[dataName][var].GetBinContent(bx),"for",dataName
                hists[dataName][var].SetBinContent(bx, 0)

def cleanVarName(var):
    """Remove unwanted characters from string for use in filenames"""
    #replace dots and slashes with underscores
    out = var.replace('.','_').replace('/','_')
    #remove other characters
    toremove = '(){}<>#$%+-=*&|[]'
    return out.translate(None,toremove)

def computeEmptyBinErrs(hists, unrollBins, aggregate=True):
    """
    For each MC histogram, assign an uncertainty on each empty bin equal to
    1.83/2 (symmetrized Poisson confidence interval)
    times the average event weight in the histogram.
    If aggregate is True, returns a dictionary {binNum:uncertainty}
    If it's False, returns a dictionary {process:{binNum:uncertainty}}
    """
    toExclude = [ # processes with low stats but not expected to contribute
            "WJetsToLNu_Wpt-0To50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
            "WJetsToLNu_Wpt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
            'QCD'
            ]
    print "Computing uncertainties to apply to empty bins"
    errs = {}
    for proc, hs in hists.iteritems():
        if proc in ['Data', 'Fit', 'Sys']: continue
        curDict = errs
        if not aggregate:
            errs[proc] = {}
            curDict = errs[proc]
        hist = hs[('MR','Rsq')]
        if hist.GetEntries() > 0: 
            avgWeight = hist.Integral() / hist.GetEntries()
        else:
            avgWeight = 0.0
        err = 1.83/2 * avgWeight
        unrolled = plotting.unroll2DHistograms([hist], 
                unrollBins[0], unrollBins[1])[0]
        for ix in range(1, unrolled.GetSize()-1):
            if ix not in curDict:
                curDict[ix] = 0.0
            if unrolled.GetBinContent(ix) == 0 and proc not in toExclude and err > 0:
                curDict[ix] = (curDict[ix]*curDict[ix] + err*err)**(0.5)
    if aggregate:
        print "Overall uncertainties on empty bins:"
        for b in sorted(errs.keys()):
            print "{} : {}".format(b, errs[b])
    return errs

def basicPrint(histDict, mcNames, varList, c, printName="Hist", dataName="Data", logx=False, ymin=0.1, lumistr="40 pb^{-1}", boxName=None, btags=None, comment=True, blindBins=None, nsigmaFitData=None, nsigmaFitMC=None, doDensity=False, printdir=".", special="", unrollBins=(None,None), vartitles={}, debugLevel=0, emptyBinErrs=None):
    """Make stacked plots of quantities of interest, with data overlaid"""
    print "Preparing histograms for plotting..."
    #format MC histograms
    for name in mcNames: 
        for var in histDict[name]: plotting.setHistColor(histDict[name][var], name)
    mcTitles = {'WJets':'W+Jets', 'DYJets':'Z #rightarrow ll','TTJets':'t#bar{t}+Jets', 'TTJets2L':'2l t#bar{t}+Jets', 'TTJets1L':'1l t#bar{t}+Jets', 'SingleTop':'Single top', 'QCD':'QCD', 'ZInv':'Z #rightarrow #nu #nu', 'GJets':'#gamma+Jets', 'GJetsInv':'#gamma+Jets', 'GJetsFrag':'#gamma+Jets (frag.)','Other':'Other'}
    titles = {name:(mcTitles[name] if name in mcTitles else name) for name in mcNames}

    #get data histograms
    if dataName in histDict: dataHists = histDict[dataName]
    else: dataHists = None

    #make correct comment string
    if comment:
        if boxName is None or boxName == "":
            commentstr = printName
        else:
            commentstr = '#it{'+boxName
            if btags is not None and btags >= 0:
                commentstr += " "+str(btags)+" b-tag"
            #if 'sideband' in special:
                #commentstr += " Sideband Fit"
            #elif 'full' in special:
                #commentstr += ' Full Fit'
            commentstr += '}'
    else:
        commentstr = ""

    legend=None
    plotFit = ("Fit" in histDict)
    for i,var in enumerate(varList): 
        print "Variable:",var
        if not isinstance(var, basestring) and len(var) == 2: #2D plots
            unrollBinsToUse = (None,None)
            if 'MR' in var[0] and 'Rsq' in var[1]:
                unrollBinsToUse = unrollBins
            mcDict = None 
            if len(mcNames) > 0:
                mcDict = {} #for stacked unrolled plots
                mcPrediction = histDict[mcNames[0]][var].Clone(histDict[mcNames[0]][var].GetName()+"mcPrediction")
                mcPrediction.Reset()
                for name in mcNames: 
                    mcPrediction.Add(histDict[name][var])
                    mcDict[name] = histDict[name][var] #for stacked unrolled plots
            else:
                mcPrediction = 0
            #copy data and fit histograms
            if dataHists is not None: obsData = dataHists[var].Clone("obsData")
            else: obsData = 0
            if not plotFit: 
                fitPrediction = 0
            else: 
                fitPrediction = histDict["Fit"][var]
            #blind signal region if necessary
            if blindBins is not None and dataHists is not None:
                blindHistograms([obsData], blindBins)
                if nsigmaFitData is not None:
                    blindHistograms([nsigmaFitData], blindBins)
            #get axis titles
            if var[0] in vartitles:
                xtitle = vartitles[var[0]]
            else:
                xtitle = var[0]
            if var[1] in vartitles:
                ytitle = vartitles[var[1]]
            else:
                ytitle = var[1]
            printstr = "%s_%s"%(cleanVarName(var[0]+var[1]),printName)
            #print histograms
            if debugLevel > 0:
                if dataHists is not None:
                    print "Data histogram:"
                    obsData.Print("all")
                if len(mcNames) > 0:
                    print "MC histogram:"
                    mcPrediction.Print("all")
                    for name in mcNames:
                        print name,"histogram:"
                        mcDict[name].Print("all")
                
            #make plots
            ratioMin = 0.0
            ratioMax = 5.0
            if 'zoomratio' in special:
                ratioMin = 0.5
                ratioMax = 1.5
            if 'SUS15004' in special:
                plotting.plot_SUS15004(c, data=obsData, fit=fitPrediction, printstr=printstr, 
                        lumistr=lumistr, commentstr=commentstr, mcDict=mcDict, mcSamples=mcNames, 
                        unrollBins=unrollBinsToUse, printdir=printdir, emptyBinErrs=emptyBinErrs, ratiomax=ratioMax, ratiomin=ratioMin)
                #do MC total (no stack)
                plotting.plot_SUS15004_FitVsMCTotal(c, mcTotal=mcPrediction, fit=fitPrediction, 
                        printstr=printstr+'MCTotal', lumistr=lumistr, commentstr=commentstr, 
                        unrollBins=unrollBinsToUse, printdir=printdir)
            elif 'CR15004' in special:
                plotting.plot_SUS15004(c, data=obsData, fit=fitPrediction, printstr=printstr, 
                        lumistr=lumistr, commentstr=commentstr, mcDict=mcDict, mcSamples=mcNames, 
                        unrollBins=unrollBinsToUse, printdir=printdir, ratiomin=0, ratiomax=2, controlRegion=True)
            else:
                plotting.plot_basic_2D(c, mc=mcPrediction, data=obsData, fit=fitPrediction, xtitle=xtitle, 
                        ytitle=ytitle, printstr=printstr, lumistr=lumistr, commentstr=commentstr, 
                        saveroot=True, savepdf=True, savepng=True, nsigmaFitData=nsigmaFitData, 
                        nsigmaFitMC=nsigmaFitMC, mcDict=mcDict, mcSamples=mcNames, ymin=ymin, 
                        unrollBins=unrollBinsToUse, printdir=printdir)
                #do MC total (no stack)
                plotting.plot_basic_2D(c, mc=mcPrediction, data=obsData, fit=fitPrediction, xtitle=xtitle, 
                        ytitle=ytitle, printstr=printstr+'MCTotal', lumistr=lumistr, 
                        commentstr=commentstr, saveroot=True, savepdf=True, savepng=True, 
                        nsigmaFitData=nsigmaFitData, nsigmaFitMC=nsigmaFitMC, ymin=ymin, unrollBins=unrollBinsToUse, 
                        printdir=printdir)
            #draw bin mapping
            if unrollBinsToUse[0] is not None and unrollBinsToUse[1] is not None:
                plotting.drawUnrolledBinMapping(c, unrollBinsToUse, xtitle=xtitle, ytitle=ytitle, printstr=printstr+"BINNING", printdir=printdir)
        #for other variables make 1D plots
        if not isinstance(var, basestring): continue #only consider strings
        varHists = {name:histDict[name][var].Clone(histDict[name][var].GetName()+"Clone") for name in mcNames}
        if doDensity: #divide each histogram bin by its width
            for name in varHists:
                for bx in range(1,varHists[name].GetNbinsX()+1):
                    varHists[name].SetBinContent(bx, varHists[name].GetBinContent(bx)*1.0/varHists[name].GetXaxis().GetBinWidth(bx))
                    varHists[name].SetBinError(bx, varHists[name].GetBinError(bx)*1.0/varHists[name].GetXaxis().GetBinWidth(bx))
        stack = plotting.makeStack(varHists, mcNames, var)
        if len(mcNames) == 0: stack = None #plotting function can't handle an empty stack
        if dataHists is not None:
            obsData = dataHists[var].Clone(dataHists[var].GetName()+"obsData")
            if doDensity:
                for bx in range(1,obsData.GetNbinsX()+1):
                    obsData.SetBinContent(bx, obsData.GetBinContent(bx)*1.0/obsData.GetXaxis().GetBinWidth(bx))
                    obsData.SetBinError(bx, obsData.GetBinError(bx)*1.0/obsData.GetXaxis().GetBinWidth(bx))
        else:
            obsData = None
        if not plotFit:
            fitPrediction = None
        else:
            fitPrediction = histDict["Fit"][var].Clone(histDict["Fit"][var].GetName()+"FitPrediction")
            if doDensity:
                for bx in range(1,fitPrediction.GetNbinsX()+1):
                    fitPrediction.SetBinContent(bx, fitPrediction.GetBinContent(bx)*1.0/fitPrediction.GetXaxis().GetBinWidth(bx))
                    fitPrediction.SetBinError(bx, fitPrediction.GetBinError(bx)*1.0/fitPrediction.GetXaxis().GetBinWidth(bx))
        if not legend:
            legend = rt.TLegend(0.75,0.6,0.9,0.9)
            rt.SetOwnership(legend, False)
            if blindBins is None and obsData is not None: legend.AddEntry(obsData, dataName)
            if plotFit and fitPrediction is not None: legend.AddEntry(fitPrediction, "Fit")
            for name in reversed(mcNames): legend.AddEntry(varHists[name], titles[name], 'f')
        if doDensity:
            if var in ['MR','MR_NoW','MR_NoZ','MR_NoPho','lep1.Pt()','genlep1.Pt()','leadingGenLeptonPt']:
                ytitle = 'Events / GeV'
            else:
                ytitle = "Events / Bin Width"
            ymin = None
        else:
            ytitle = "Events"
        #set x title
        if var in vartitles:
            xtitle = vartitles[var]
        else:
            xtitle = var
        if var in ['MR','Rsq','MR_NoW',"Rsq_NoW","MR_NoZ","Rsq_NoZ", "lep1.Pt()", "genlep1.Pt()","leadingGenLeptonPt"]: logx = True
        else: logx = False
        printstr = "%s_%s"%(cleanVarName(var),printName)
        if blindBins is None:
            if 'SUS15004' in special:
                plotting.plot_SUS15004_1D(c, mc=stack, data=obsData, leg=legend, xtitle=xtitle, ytitle=ytitle, printstr=printstr, logx=logx, lumistr=lumistr, ymin=ymin, commentstr=commentstr, printdir=printdir)
            elif 'CR15004' in special:
                plotting.plot_SUS15004_1D(c, mc=stack, data=obsData, leg=legend, xtitle=xtitle, ytitle=ytitle, printstr=printstr, logx=logx, lumistr=lumistr, ymin=ymin, commentstr=commentstr, printdir=printdir)
            else:
                plotting.plot_basic(c, mc=stack, data=obsData, fit=fitPrediction, leg=legend, xtitle=xtitle, ytitle=ytitle, printstr=printstr, logx=logx, lumistr=lumistr, ymin=ymin, commentstr=commentstr, saveroot=True, savepdf=True, savepng=True, printdir=printdir)
        else:
            plotting.plot_basic(c, mc=stack, data=None, fit=fitPrediction, leg=legend, xtitle=xtitle, ytitle=ytitle, printstr=printstr, logx=logx, lumistr=lumistr, ymin=ymin, commentstr=commentstr, saveroot=True, savepdf=True, savepng=True, printdir=printdir)
        legend.Delete()

def transformVarsInString(string, varNames, suffix):
    outstring = copy.copy(string)
    for var in varNames:
        outstring = outstring.replace(var, var+suffix) 
    return outstring

JETVARS = ["MR","Rsq","nBTaggedJets","dPhiRazor","leadingJetPt","subleadingJetPt","nSelectedJets","nJets80","mT","box"] #quantities susceptible to jet uncertainties
def transformVarString(event, string, errorOpt, process="", debugLevel=0):
    outstring = string
    if errorOpt == "jesUp":
        outstring = transformVarsInString(string, JETVARS, "_JESUp")
    elif errorOpt == "jesDown":
        outstring = transformVarsInString(string, JETVARS, "_JESDown")
    elif errorOpt == "mesUp":
        outstring = transformVarsInString(string, JETVARS, "_MESUp")
    elif errorOpt == "mesDown":
        outstring = transformVarsInString(string, JETVARS, "_MESDown")
    elif errorOpt == "eesUp":
        outstring = transformVarsInString(string, JETVARS, "_EESUp")
    elif errorOpt == "eesDown":
        outstring = transformVarsInString(string, JETVARS, "_EESDown")
    elif errorOpt == "jerUp":
        outstring = transformVarsInString(string, JETVARS, "_JERUp")
    elif errorOpt == "jerDown":
        outstring = transformVarsInString(string, JETVARS, "_JERDown")
    elif errorOpt == "npvextrapUp":
        if '&&' in string: # only applies to cut string
            outstring = string + " && nVtx >= 20"
    elif errorOpt == "npvextrapDown":
        if '&&' in string: # only applies to cut string
            outstring = string + " && nVtx < 20"
    # use gen-level MET instead of PF MET
    elif errorOpt == "genmetvspfmetUp":
        outstring = transformVarsInString(string, ['Rsq', 'mT'], 'GenMet')
    # 'down' option does not exist for this option; select no events
    elif errorOpt == "genmetvspfmetDown": 
        if '&&' in string:
            print ("This 'down' histogram will be filled manually "
                    "during postprocessing; skipping for now")
            outstring = "0"

    if debugLevel > 1:
        if outstring != string: print "For option",errorOpt,"Replacing string '",string,"' with '",outstring,"'"
    return outstring

def getAdditionalCuts(tree, errorOpt, process, debugLevel=0):
    """Implement here any event-by-event cuts that can't be handled with TTree::Draw"""
    pass

def basicFill(tree, hists={}, weight=1.0, sysErrSquaredHists={}, sysErr=0.0, errorOpt=None, additionalCuts=None, formulas={}, debugLevel=0):
    """Fills each histogram with the corresponding variable in the tree.
    'hists' should be a dictionary of histograms, with keys being the variable names to fill.
    Ex: hists['MR'] should be the histogram you want to fill with MR values.
    A key that is a tuple of variables (ex: ('MR','Rsq')) should be paired with a multidimensional histogram.
    In this case, the given variables will be used to fill the histogram."""
    # if true, fill each systematic uncertainty separately
    splitSysErrs = isinstance(sysErr, dict) 
    #make additional cuts 
    if additionalCuts is not None:
        additionalCutsBool = eval(additionalCuts)
        if debugLevel > 1: 
            print "Additional cuts:",additionalCuts,"\nDecision:",additionalCutsBool
        if not additionalCutsBool: return
    #fill each variable
    for varName, hist in hists.iteritems(): 
        if isinstance(varName, basestring): #if varName is a string
            #transform variable name
            if errorOpt is not None: varName = transformVarString(tree, varName, errorOpt, debugLevel=debugLevel)
            if varName in formulas: #if TTreeFormula
                varValue = formulas[varName].EvalInstance()
            else: #if tree variable
                varValue =  getattr(tree, varName)
            hist.Fill(varValue, weight)
            if debugLevel > 1: print "Filling",varName,"=",varValue,"with weight",weight
            if splitSysErrs:
                # each uncertainty goes into its own histogram
                for sf, err in sysErr.iteritems():
                    if sf not in sysErrSquaredHists:
                        continue
                    sysErrSquared = weight * weight * err * err
                    sysErrSquaredHists[sf][varName].Fill(
                            varValue, sysErrSquared)
            elif varName in sysErrSquaredHists: 
                # there is only one uncertainty to fill
                sysErrSquared = weight*weight*sysErr*sysErr
                sysErrSquaredHists[varName].Fill(varValue, sysErrSquared)
                if debugLevel > 1: print "Sys. Error =",sysErr,"; Filling (w*sysErr)^2 histogram with",sysErrSquared
        else: #treat it as a tuple of variables that should be filled
            #transform each variable
            if errorOpt is not None: varName = tuple([transformVarString(tree, v, errorOpt, debugLevel=debugLevel) for v in varName])
            toFillBase = [formulas[v].EvalInstance() if v in formulas
                    else getattr(tree, v) for v in varName]
            toFill = toFillBase+[weight]
            if debugLevel > 1: print "Filling",varName,":",toFill
            hist.Fill(*toFill)
            if splitSysErrs:
                for sf, err in sysErr.iteritems():
                    if sf not in sysErrSquaredHists:
                        continue
                    sysErrSquared = weight * weight * err * err
                    toFillErr = toFillBase + [sysErrSquared]
                    sysErrSquaredHists[sf][varName].Fill(*toFillErr)
            elif varName in sysErrSquaredHists:
                sysErrSquared = weight*weight*sysErr*sysErr
                toFillErr = toFillBase+[sysErrSquared]
                sysErrSquaredHists[varName].Fill(*toFillErr)
                if debugLevel > 1: print "Sys. Error =",sysErr,"; Filling (w*sysErr)^2 histogram with",sysErrSquared

def makeTreeDict(fileDict, treeName, debugLevel=0):
    """gets a tree called treeName from each file in fileDict, and returns a dict of trees"""
    trees = {}
    for name in fileDict:
        if debugLevel > 0: print("Loading tree "+treeName+" for process "+name)
        trees[name] = fileDict[name].Get(treeName)
        if debugLevel > 0: print("Got tree containing "+str(trees[name].GetEntries())+" entries")
        assert trees[name]
    if debugLevel > 0: 
        print("Trees loaded: ") 
        print trees
    return trees

def getScaleFactorAndError(tree, sfHist, sfVars=("MR","Rsq"), formulas={}, debugLevel=0):
    if debugLevel > 1:
        print "Getting scale factor from histogram",sfHist.GetName()
    if isinstance(sfVars, basestring):
        sfVars = (sfVars,) #cast into tuple 
    #get variables
    var = [formulas[v].EvalInstance() if v in formulas else getattr(tree, v) for v in sfVars]
    #constrain variables to be within the bounds of the histogram
    var[0] = min(var[0], sfHist.GetXaxis().GetXmax()*0.999)
    var[0] = max(var[0], sfHist.GetXaxis().GetXmin()*1.001)
    if len(var) > 1:
        var[1] = min(var[1], sfHist.GetYaxis().GetXmax()*0.999)
        var[1] = max(var[1], sfHist.GetYaxis().GetXmin()*1.001)
    if len(var) > 2:
        var[2] = min(var[2], sfHist.GetZaxis().GetXmax()*0.999)
        var[2] = max(var[2], sfHist.GetZaxis().GetXmin()*1.001)
    #TH2Poly case
    if sfHist.InheritsFrom('TH2Poly'):
        scaleFactor = sfHist.GetBinContent(sfHist.FindBin(*var))
        scaleFactorErr = sfHist.GetBinError(sfHist.FindBin(*var))
    #TH2F case
    else:
        scaleFactor = sfHist.GetBinContent(sfHist.FindFixBin(*var))
        scaleFactorErr = sfHist.GetBinError(sfHist.FindFixBin(*var))

    #add protection against unphysics scale factors
    #if scale factor is 0, then use scale factor of 1 and add 100% systematic uncertainty
    if scaleFactor == 0:
        scaleFactor = 1.0
        scaleFactorErr = math.sqrt( 1 + scaleFactorErr*scaleFactorErr )
    #protect against NaN uncertainties
    if math.isnan(scaleFactorErr):
        scaleFactorErr = 0.0

    if debugLevel > 1: print "Applying scale factor: ",scaleFactor,"with error",scaleFactorErr
    return (scaleFactor, scaleFactorErr)

def addToTH2ErrorsInQuadrature(hists, sysErrSquaredHists, debugLevel=0):
    """For each histogram in hists, look for the corresponding histogram in sysErrSquaredHists.
    Treats the values of sysErrSquaredHists as sums of (weight*error)^2 in each bin, and adds these errors in quadrature with the existing bin errors in hists"""
    for name in hists:
        if name in sysErrSquaredHists:
            if debugLevel > 0: print "Including systematic errors on ",name
            for bx in range(hists[name].GetSize()+1):
                squaredError = sysErrSquaredHists[name].GetBinContent(bx)
                oldErr = hists[name].GetBinError(bx)
                hists[name].SetBinError(bx,(oldErr*oldErr + squaredError)**(0.5))
                if debugLevel > 0 and squaredError > 0: print name,": Error on bin",bx,"increases from",oldErr,"to",hists[name].GetBinError(bx),"after adding",(squaredError**(0.5)),"in quadrature"

def propagateScaleFactorStatErrors(sysErrSquaredHists, upHists, downHists, debugLevel=0):
    """Use scale factor uncertainties to make up and down shape histograms"""
    for var in upHists:
        if var in sysErrSquaredHists:
            if debugLevel > 0: print "Making shape histogram for scale factor errors on",var
            for bx in range(1,sysErrSquaredHists[var].GetSize()+1):
                squaredError = sysErrSquaredHists[var].GetBinContent(bx)
                #increase up histogram bin content by 1 sigma, decrease down histogram bin content
                centralValue = upHists[var].GetBinContent(bx)
                upHists[var].SetBinContent(bx, centralValue + (squaredError)**(0.5))
                if centralValue > 0:
                    percentChange = (upHists[var].GetBinContent(bx) - centralValue)/centralValue
                    downHists[var].SetBinContent(bx, centralValue/(1+percentChange))

def isCustomBadEvent(event):
    """
    Returns True if this event needs to be removed.
    (nothing personal)
    """
    # These are ECAL spike events listed on
    # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    if not hasattr(event, 'event'):
        return False
    if not hasattr(event, 'run'):
        return False
    badEvents = [(597707687, 274244),
                 (1018772012, 274338),
                 (352051316, 274200),
                 (21612513, 274338)]
    r = event.run
    e = event.event
    # Event numbers are stored in the tree as unsigned 32-bit ints.
    # When pyROOT reads the tree it casts the number into a signed int
    # so if the event number is too large it will show as negative.
    # Fix this by adding 2^32 to the event number.
    if e < 0:
        e = e + 2**32
    for badE, badR in badEvents:
        if e == badE and r == badR:
            print "Bad event: {} {}".format(e, r)
            return True
    return False

def loopTree(tree, weightF, cuts="", hists={}, weightHists={}, sfHist=None, scale=1.0, fillF=basicFill, sfVars=("MR","Rsq"), statErrOnly=False, weightOpts=[], errorOpt=None, process="", auxSFs={}, auxSFHists={}, shapeHists={}, shapeNames=[], shapeSFHists={}, shapeAuxSFs={}, shapeAuxSFHists={}, propagateScaleFactorErrs=True, noFill=False, debugLevel=0):
    """Loop over a single tree and fill histograms.
    Returns the sum of the weights of selected events.
    Arguments:
        tree: ROOT tree to process
        weightF: python function that computes event weights
        cuts: string indicating the cuts that should be applied
        hists: dictionary { "var1":hist1, ... }
            where "var1" is the name of a variable and hist1 is a ROOT histogram
        weightHists: dictionary { "name1":hist1, ... }
        sfHist: ROOT histogram containing scale factors
        scale: scale factor to be applied to all event weights
        fillF: python function that fills histograms
        sfVars: string of tuple of strings indicating the reweighting variable
        statErrOnly: if True, do not compute systematic errors
        weightOpts: list of strings indicating weighting options
        errorOpt: uncertainty option that should be applied
        process: name of the physics process
        auxSFs: dictionary containing additional scale factor options
            (see makeControlSampleHistsForAnalysis for description)
        auxSFHists: dictionary { "name1":hist1, ... }
            where "name1" is a reweighting option included as a key in "auxSFs"
        shapeHists: dictionary { "unc1":hists1, ... }
            where "unc1", ... are shape uncertainties and hists1, ... are 
            dictionaries of histograms as for the "hists" argument
        shapeNames: list of shape uncertainty names
        shapeSFHists: scale factor histograms needed to compute shape uncertainties:
            { "unc1":hists1, ... }
            where hists1, ... are dictionaries of scale factor histograms
        shapeAuxSFs: dictionary { "unc1":auxSFs1, ... }
            where auxSFs1, ... are dictionaries of scale factors as the "auxSFs" argument
        shapeAuxSFHists: dictionary { "unc1":hists1, ... }
            where hists1, .. are dictionaries of scale factor histograms
        propagateScaleFactorErrs: if True, uncertainties from scale factors are propagated to
            the filled histograms
        noFill: if True, histograms are not filled
        debugLevel: 0 for standard mode, 1 for verbose mode, 2 for debug mode"""

    if debugLevel > 0: print ("Looping tree "+tree.GetName())
    #make TTreeFormulas for variables not found in the tree
    formulas = {}
    for var in hists:
        if isinstance(var, basestring) and not hasattr(tree, var): #if it's not in the tree
            formulas[var] = rt.TTreeFormula(var, var, tree)
            formulas[var].GetNdata()
    #make TTreeFormulas for auxSFs
    auxSFForms = {}
    for name,pair in auxSFs.iteritems(): #auxSFs pairs look like "ScaleFactor":("Var to reweight","Cuts")
        if debugLevel > 0: print "Making TTreeFormula for",name,"with formula",pair[1],"for reweighting",pair[0]
        auxSFForms[name] = rt.TTreeFormula(name+"Cuts", pair[1], tree)
        auxSFForms[name].GetNdata()
        #add the reweighting variable to the formula list
        if isinstance(pair[0], basestring):
            if pair[0] not in formulas:
                formulas[pair[0]] = rt.TTreeFormula(pair[0], pair[0], tree)
                formulas[pair[0]].GetNdata()
        else: #treat as tuple
            for var in pair[0]:
                if var not in formulas:
                    formulas[var] = rt.TTreeFormula(var, var, tree)
                    formulas[var].GetNdata()

    #make TTreeFormulas for shape histogram scale factors
    #A list of TTreeFormulas is maintained in shapeAuxSFForms, and shapeAuxSFFormsLookup matches each shape uncertainty with the index of a TTreeFormula in shapeAuxSFForms.
    shapeAuxSFForms = []
    shapeAuxSFFormsLookup = {} #stores the number of the TTreeFormula for each shape uncertainty cut
    cutStrings = []
    for n in shapeNames:
        shapeAuxSFFormsLookup[n+'Up'] = {}
        shapeAuxSFFormsLookup[n+'Down'] = {}
        for name,pair in shapeAuxSFs[n+'Up'].iteritems():
            #add reweighting variable to formula list
            if pair[0] not in formulas and isinstance(pair[0], basestring):
                formulas[pair[0]] = rt.TTreeFormula(pair[0], pair[0], tree)
                formulas[pair[0]].GetNdata()
            if pair[1] not in cutStrings:
                #make the TTreeFormula for this cut string
                if debugLevel > 0: print "Making TTreeFormula for",name,"with formula",pair[1],"for reweighting",pair[0],"(error option",n+"Up)"
                shapeAuxSFFormsLookup[n+'Up'][name] = len(shapeAuxSFForms) #index of TTreeFormula to use
                shapeAuxSFForms.append(rt.TTreeFormula("shapeAuxSFCuts"+str(len(shapeAuxSFForms)), pair[1], tree))
                shapeAuxSFForms[-1].GetNdata()
                cutStrings.append(pair[1])
            else:  
                #refer to the matching TTreeFormula that already exists
                index = cutStrings.index(pair[1])
                shapeAuxSFFormsLookup[n+'Up'][name] = index
        for name,pair in shapeAuxSFs[n+'Down'].iteritems():
            #add reweighting variable to formula list
            if pair[0] not in formulas and isinstance(pair[0], basestring):
                formulas[pair[0]] = rt.TTreeFormula(pair[0], pair[0], tree)
                formulas[pair[0]].GetNdata()
            if pair[1] not in cutStrings:
                #make the TTreeFormula for this cut string
                if debugLevel > 0: print "Making TTreeFormula for",name,"with formula",pair[1],"for reweighting",pair[0],"(error option",n+"Down)"
                shapeAuxSFFormsLookup[n+'Down'][name] = len(shapeAuxSFForms) #index of TTreeFormula to use
                shapeAuxSFForms.append(rt.TTreeFormula("shapeAuxSFCuts"+str(len(shapeAuxSFForms)), pair[1], tree))
                shapeAuxSFForms[-1].GetNdata()
                cutStrings.append(pair[1])
            else:  
                #refer to the matching TTreeFormula that already exists
                index = cutStrings.index(pair[1])
                shapeAuxSFFormsLookup[n+'Down'][name] = index

    #transform cuts 
    if errorOpt is not None:
        if debugLevel > 0: print "Error option is:",errorOpt
        cuts = transformVarString(tree, cuts, errorOpt, process=process, debugLevel=debugLevel+1)
    print ("Cuts: "+cuts)
    #create histograms for scale factor systematics
    sysErrSquaredHists = {}
    if not statErrOnly:
        for name in hists: 
            sysErrSquaredHists[name] = hists[name].Clone(hists[name].GetName()+"SFERRORS")
            sysErrSquaredHists[name].Reset()
            if debugLevel > 0: print "Created temp histogram",sysErrSquaredHists[name].GetName(),"to hold systematic errors from",sfVars,"scale factors"
        # If scale factor uncertainties will be saved for further
        # processing, we have to make a separate histogram for each
        # uncertainty source.
        if not propagateScaleFactorErrs:
            sysErrSquaredHists = {'': sysErrSquaredHists}
            for sf in auxSFs:
                sysErrSquaredHists[sf] = {}
                for name in hists:
                    sysErrSquaredHists[sf][name] = hists[name].Clone(
                            hists[name].GetName()+"SFERRORS"+sf)
                    sysErrSquaredHists[sf][name].Reset()
    count = 0
    sumweight = 0.0
    if not noFill:
        #get list of entries passing the cuts
        tree.Draw('>>elist', cuts, 'entrylist')
        elist = rt.gDirectory.Get('elist')
        print "Total entries passing cuts:",elist.GetN()
        while True:
            #load the next entry
            entry = elist.Next()
            if entry == -1: break
            if count > 0 and count % 100000 == 0: print "Processing entry",count
            elif debugLevel > 0 and count % 10000 == 0: print "Processing entry",count
            elif debugLevel > 1: print "Processing entry",count
            tree.GetEntry(entry)
            if isCustomBadEvent(tree):
                continue
            w = weightF(tree, weightHists, scale, weightOpts, errorOpt, debugLevel=debugLevel)

            err = 0.0
            #get scale factor
            if sfHist is not None: 
                sf, err = getScaleFactorAndError(tree, sfHist, sfVars, formulas, debugLevel=debugLevel)
                #apply scale factor to event weight
                w *= sf
                if not propagateScaleFactorErrs:
                    err = {'': err} # keep track of each unc source separately
            for name in auxSFs: #apply misc scale factors (e.g. veto lepton correction)
                if auxSFForms[name].EvalInstance(): #check if this event should be reweighted
                    auxSF, auxErr = getScaleFactorAndError(tree, auxSFHists[name], sfVars=auxSFs[name][0], formulas=formulas, debugLevel=debugLevel)
                    w *= auxSF
                    if propagateScaleFactorErrs:
                        err = (err*err + auxErr*auxErr)**(0.5)
                    else:
                        err[name] = auxErr
            #protection for case of infinite-weight events
            if math.isinf(w):
                if debugLevel > 0: 
                    print "Warning: infinite-weight event encountered!"
                continue
            if propagateScaleFactorErrs and math.isnan(err):
                print "Error: err is nan!"
            #fill with weight
            fillF(tree, hists, w, sysErrSquaredHists, err, errorOpt, additionalCuts=None, formulas=formulas, debugLevel=debugLevel)

            ###PROCESS EXTRA SHAPE HISTOGRAMS
            #get weights for filling additional shape histograms, if provided
            wsUp = [weightF(tree, weightHists, scale, weightOpts, e+'Up', debugLevel=debugLevel) for e in shapeNames]
            wsDown = [weightF(tree, weightHists, scale, weightOpts, e+'Down', debugLevel=debugLevel) for e in shapeNames]
            #evaluate all TTreeFormulas 
            shapeAuxSFFormResults = [form.EvalInstance() for form in shapeAuxSFForms]
            if debugLevel > 1 and len(shapeAuxSFFormResults) > 0:
                print "Cut decisions for systematic uncertainty scale factors:",shapeAuxSFFormResults
            #apply scale factor to shape up/down weights too (if needed)
            for i,n in enumerate(shapeNames):
                if debugLevel > 1: print "Error option:",n
                #main scale factor
                if shapeSFHists[n+'Up'] is not None:
                    shapeSFUp, shapeErrUp = getScaleFactorAndError(tree, shapeSFHists[n+'Up'], sfVars, formulas, debugLevel=debugLevel)
                    wsUp[i] *= shapeSFUp
                if shapeSFHists[n+'Down'] is not None:
                    shapeSFDown, shapeErrDown = getScaleFactorAndError(tree, shapeSFHists[n+'Down'], sfVars, formulas, debugLevel=debugLevel)
                    wsDown[i] *= shapeSFDown
                #misc scale factors
                for name in shapeAuxSFs[n+'Up']:
                    if shapeAuxSFFormResults[ shapeAuxSFFormsLookup[n+'Up'][name] ]:
                        shapeAuxSFUp, shapeAuxErrUp = getScaleFactorAndError(tree, shapeAuxSFHists[n+'Up'][name], sfVars=shapeAuxSFs[n+'Up'][name][0], formulas=formulas, debugLevel=debugLevel)
                        wsUp[i] *= shapeAuxSFUp
                for name in shapeAuxSFs[n+'Down']:
                    if shapeAuxSFFormResults[ shapeAuxSFFormsLookup[n+'Down'][name] ]:
                        shapeAuxSFDown, shapeAuxErrDown = getScaleFactorAndError(tree, shapeAuxSFHists[n+'Down'][name], sfVars=shapeAuxSFs[n+'Down'][name][0], formulas=formulas, debugLevel=debugLevel)
                        wsDown[i] *= shapeAuxSFDown
                #fill up/down shape histograms
                fillF(tree, shapeHists[n+'Up'], wsUp[i], formulas=formulas, debugLevel=debugLevel)
                fillF(tree, shapeHists[n+'Down'], wsDown[i], formulas=formulas, debugLevel=debugLevel)
            #####

            sumweight += w
            count += 1
    if propagateScaleFactorErrs:
        # directly add systematic error in quadrature with statistical
        addToTH2ErrorsInQuadrature(hists, sysErrSquaredHists, debugLevel)
    else:
        # propagate systematics to each histogram separately
        knownSFErrs = ['sfstatttjets', 'sfstatwjets', 'sfstatzinv']
        for unc in knownSFErrs:
            if unc not in shapeNames: 
                continue
            # create all the necessary histograms
            for sf in sysErrSquaredHists:
                if sf == '': 
                    continue # already have this one
                for updown in ['Up', 'Down']:
                    shapeHists[unc+sf+updown] = {n: h.Clone(
                                h.GetName().replace(unc, unc+sf))
                            for n, h in shapeHists[
                                unc+updown].iteritems()}
                    for _, h in shapeHists[unc+sf+updown].iteritems():
                        h.SetDirectory(0)
            # propagate
            for sf, hs in sysErrSquaredHists.iteritems():
                propagateScaleFactorStatErrors(
                        hs, upHists=shapeHists[unc+sf+'Up'], 
                        downHists=shapeHists[unc+sf+'Down'], 
                        debugLevel=debugLevel)
    print "Sum of weights for this sample:",sumweight
    return sumweight

def loopTrees(treeDict, weightF, cuts="", hists={}, weightHists={}, sfHists={}, scale=1.0, weightOpts=[], errorOpt=None, fillF=basicFill, sfVars=("MR","Rsq"), statErrOnly=False, boxName="NONE", auxSFs={}, shapeHists={}, shapeNames=[], shapeAuxSFs={}, noFill=False, propagateScaleFactorErrs=True, extraWeightOpts={}, extraCuts={}, debugLevel=0, weightHistsPerProcess={}):
    """calls loopTree on each tree in the dictionary.  
    Arguments:
        treeDict: dictionary { "process1":tree1, "process2":tree2, ... }
            where tree1, tree2, ... are TTrees
        weightF: python function that should be used to weight the MC
        cuts: string with cuts to be applied
        hists: dictionary { "process1":hists1, "process2":hists2, ... }
            where hists1, hists2, ... are dictionaries of the form
                { "var1":hist1, "var2":hist2, ... }
        weightHists: dictionary { "name1":hist1, ... }
            See the weighting function for details on which weight options are supported
        sfHists: dictionary { "name1":hist1, ... } of scale factor histograms
        scale: overall scale factor to apply to the events
        weightOpts: list of weighting options.
            See the weighting function for more details.
        errorOpt: indicates which uncertainty should be used
        fillF: python function that should be used to fill histograms for the event.
        sfVars: string or tuple of strings: variable to reweight
        statErrOnly: if True, do not compute systematic uncertainties
        boxName: name of analysis box (for razor signal region only)
        auxSFs: auxiliary scale factor options (see description at makeControlSampleHistsForAnalysis)
        shapeHists: dictionary { "unc1":{ ... }, ... }
            where "unc1", ... are the names of uncertainties and { ... } 
            are dictionaries of the form used by the "hists" argument
        shapeNames: list of shape uncertainty names (see makeControlSampleHistsForAnalysis)
        shapeAuxSFs: dictionary { "unc1":{ ... }, ... }
            where "unc1", ... are the names of uncertainties and { ... }
            are dictionaries of the form used by the "auxSFs" argument
        noFill: do not fill histograms
        propagateScaleFactorErrs: if True, propagate uncertainties on scale factors
            to the filled histograms
        extraWeightOpts: additional per-process weight options:
            { "proc1":[opt1, opt2,...], ... }
        extraCuts: additional per-process cuts to apply: 
            { "proc1":"cuts1", ... }
        debugLevel: 0 for standard mode, 1 for verbose mode, 2 for debug mode
        weightHistsPerProcess: dictionary { "proc1":{<weight histograms>}, ...}
            used for NEvents and SumScaleWeights histograms
    """
    sumweights=0.0
    for name in treeDict: 
        if name not in hists: continue
        print("Filling histograms for tree "+name)

        #get correct scale factor histograms
        sfHistToUse = None
        if name in sfHists:
            sfHistToUse = sfHists[name]
            print("Using scale factors from histogram "+sfHistToUse.GetName())
        if name in auxSFs and isinstance(auxSFs[name],dict):
            auxSFsToUse = auxSFs[name]
        else:
            auxSFsToUse = auxSFs
        print "Also using scale factors",auxSFsToUse
        #get appropriate scale factor histograms for misc reweightings
        auxSFHists = {name:sfHists[name] for name in auxSFsToUse} 
        #get correct weight histograms
        weightHistsToUse = copy.copy(weightHists)
        if name in weightHistsPerProcess:
            weightHistsToUse.update(weightHistsPerProcess[name])
        #get correct variables for scale factor reweighting.
        #if sfVars is a dictionary, get the appropriate value from it.  otherwise use sfVars directly.
        sfVarsToUse = sfVars
        if sfHistToUse is not None and isinstance(sfVars,dict):
            sfVarsToUse = sfVars[name]
            print "Reweighting in",sfVarsToUse
        weightOptsToUse = copy.copy(weightOpts)
        #add additional per-process weight options
        if name in extraWeightOpts:
            weightOptsToUse += extraWeightOpts[name]
            print "Using extra weight options",extraWeightOpts[name],"for this process"
        #remove inapplicable option
        if 'nbjets' in weightOptsToUse and (name != 'TTJets' or not hasattr(treeDict[name],'NGenBJets')):
            weightOptsToUse.remove('nbjets')

        #get cuts
        cutsToUse = cuts
        if name in extraCuts:
            cutsToUse += " && "+extraCuts[name]

        #prepare additional shape histograms if needed
        shapeHistsToUse = {}
        shapeNamesToUse = []
        shapeAuxSFsToUse = {}
        shapeAuxSFHists = {}
        shapeSFHists = {}
        for shape in shapeNames:
            if not isinstance(shape,basestring): #tuple (shape, [list of processes])
                if name not in shape[1]: continue
                curShape = shape[0]
            else:
                curShape = shape
            if curShape+'Up' in shapeHists[name]:
                shapeHistsToUse[curShape+'Up'] = shapeHists[name][curShape+'Up']
                shapeHistsToUse[curShape+'Down'] = shapeHists[name][curShape+'Down']
                shapeNamesToUse.append(curShape)
                if name in shapeAuxSFs[curShape+'Up'] and isinstance(shapeAuxSFs[curShape+'Up'],dict):
                    shapeAuxSFsToUse[curShape+'Up'] = shapeAuxSFs[curShape+'Up'][name]
                    shapeAuxSFsToUse[curShape+'Down'] = shapeAuxSFs[curShape+'Down'][name]
                else:
                    shapeAuxSFsToUse[curShape+'Up'] = shapeAuxSFs[curShape+'Up']
                    shapeAuxSFsToUse[curShape+'Down'] = shapeAuxSFs[curShape+'Down']
                shapeAuxSFHists[curShape+'Up'] = {n:sfHists[n] for n in shapeAuxSFsToUse[curShape+'Up']} 
                shapeAuxSFHists[curShape+'Down'] = {n:sfHists[n] for n in shapeAuxSFsToUse[curShape+'Down']} 
                shapeSFHists[curShape+'Up'] = None
                shapeSFHists[curShape+'Down'] = None
                if name in sfHists:
                    shapeSFHists[curShape+'Up'] = sfHists[getSFHistNameForErrorOpt(curShape+'Up', name)]
                    shapeSFHists[curShape+'Down'] = sfHists[getSFHistNameForErrorOpt(curShape+'Down', name)]
        if debugLevel > 0:
            print "Will fill histograms for these shape uncertainties:"
            print shapeNamesToUse

        sumweights += loopTree(treeDict[name], weightF, cutsToUse, hists[name], weightHistsToUse, sfHistToUse, scale, fillF, sfVarsToUse, statErrOnly, weightOptsToUse, errorOpt, process=name+"_"+boxName, auxSFs=auxSFsToUse, auxSFHists=auxSFHists, shapeHists=shapeHistsToUse, shapeNames=shapeNamesToUse, shapeSFHists=shapeSFHists, shapeAuxSFs=shapeAuxSFsToUse, shapeAuxSFHists=shapeAuxSFHists, noFill=noFill, propagateScaleFactorErrs=propagateScaleFactorErrs, debugLevel=debugLevel)
        # deal with the case where loopTree creates new histograms
        for sf in shapeHistsToUse:
            shapeHists[name][sf] = shapeHistsToUse[sf]
    print "Sum of event weights for all processes:",sumweights

def correctScaleFactorUncertaintyForSignalContamination(centralHist, upHist, downHist, sfHist, contamHist, debugLevel=0):
    """
    sigHist should be the histogram that needs to be corrected.
    sfHist should be the histogram of scale factors.
    contamHist should be a histogram of the same binning as sfHist giving % signal contamination in each bin of sfHist.
    Increases the uncertainty on each bin of sfHist to account for the level of signal contamination.
    Assumes that the signal histogram uncertainties reflect MC statistics only!  
    (i.e. no other systematics have been propagated)
    Supports TH2s as well as TH2Polys for the scale factor histogram
    """
    if debugLevel > 0:
        print "Adding uncertainty for signal contamination in",sfHist.GetName()
    if sfHist.InheritsFrom('TH2Poly'):
        for bx in range(1,centralHist.GetNbinsX()+1):
            for by in range(1, centralHist.GetNbinsY()+1):
                #find the scale factor bin corresponding to this signal region bin
                xCoord = centralHist.GetXaxis().GetBinCenter(bx)
                yCoord = centralHist.GetYaxis().GetBinCenter(by)
                sfBin = sfHist.FindBin(xCoord, yCoord)
            
                #get error, scale factor error, and level of signal contamination
                sf = sfHist.GetBinContent(sfBin)
                contam = contamHist.GetBinContent(sfBin)
                contamErr = sf*contam

                curErr = upHist.GetBinContent(bx,by) - centralHist.GetBinContent(bx,by)
                newErr = ( curErr**2 + contamErr**2 )**(0.5)

                upHist.SetBinContent(bx,by, centralHist.GetBinContent(bx,by) + newErr)
                if centralHist.GetBinContent(bx,by) > 0:
                    percentChange = newErr/centralHist.GetBinContent(bx,by)
                    downHist.SetBinContent(bx,by, centralHist.GetBinContent(bx,by)/(1+percentChange))

                if debugLevel > 0:
                    print "Signal contamination in bin",bx,by,"is",contam,"; uncertainty goes from",curErr,"to",newErr
    else:
        print "Error in correctScaleFactorUncertaintyForSignalContamination: function implemented for TH2Poly only!"
            
def getExcludedSignalStrength(dirName, model, mGluino=-1, mStop=-1, mLSP=-1, debugLevel=0): 
    """ Retrieve the expected signal exclusion for the given model, using previously computed limits """ 
    if mLSP < 0:
        print "Error in getExcludedSignalStrength: please specify mLSP!"
        return 0

    #open file and get limit results
    fName = model+'_MultiJet_EleMultiJet_MuMultiJet_results.root'
    hName = "xsecUL_Exp_%s_MultiJet_EleMultiJet_MuMultiJet"%(model)
    resultFile = rt.TFile.Open(dirName+'/'+fName)
    xsecULHist = resultFile.Get(hName)
    if not xsecULHist:
        print "Error in getExcludedSignalStrength: histogram",hName,"not found in",dirName+'/'+fName
        return 0

    #get excluded cross section
    if 'T2' in model:
        xsecUL = xsecULHist.GetBinContent(xsecULHist.FindFixBin(mStop,mLSP))
    else:
        xsecUL = xsecULHist.GetBinContent(xsecULHist.FindFixBin(mGluino,mLSP))

    #get theoretical cross section
    if mGluino!=-1:
        for line in open('data/gluino13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(mGluino))==line.split(',')[0]:
                xsecTheory = float(line.split(',')[1]) #pb
    elif mStop!=-1:
        for line in open('data/stop13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(mStop))==line.split(',')[0]:
                xsecTheory = float(line.split(',')[1]) #pb
    else:
        print "Error in getExcludedSignalStrength: please specify either mStop or mGluino!"
        return 0
    if debugLevel > 0:
        print "Expected UL",xsecUL,"pb"
        print "Theory cross section",xsecTheory,"pb"

    #return the ratio
    return xsecUL/xsecTheory

def doDeltaBForReducedEfficiencyMethod(backgroundHists, signalHists, contamHists, sfHists, unrollBins, debugLevel=0):
    for proc, contamHist in contamHists.iteritems():
        #get change in background histogram due to signal contamination
        backgroundHist = backgroundHists[proc]
        sfHist = sfHists[proc]
        bn = 0 #keep track of unrolled bin number
        for btagUnrollBins in unrollBins: #loop over btags
            for bx in range(len(btagUnrollBins[0])-1): #loop over columns
                xCoord = (btagUnrollBins[0][bx+1] + btagUnrollBins[0][bx])/2.0
                for by in range(len(btagUnrollBins[1][bx])-1): #loop over y-bins in this column
                    bn += 1 
                    yCoord = (btagUnrollBins[1][bx][by] + btagUnrollBins[1][bx][by+1])/2.0
                    #find the scale factor bin corresponding to this signal region bin
                    xCoord = min(xCoord, sfHist.GetXaxis().GetXmax()*0.999)
                    xCoord = max(xCoord, sfHist.GetXaxis().GetXmin()*1.001)
                    yCoord = min(yCoord, sfHist.GetYaxis().GetXmax()*0.999)
                    yCoord = max(yCoord, sfHist.GetYaxis().GetXmin()*1.001)
                    sfBin = sfHist.FindBin(xCoord, yCoord)
                
                    #get level of signal contamination
                    sf = sfHist.GetBinContent(sfBin)
                    contam = contamHist.GetBinContent(sfBin)
                    background = backgroundHist.GetBinContent(bn)
                    if sf > 0:
                        deltaB = (background/sf)*contam #amount the background should be corrected for sig. contam.

                    #subtract deltaB from all of the signal histograms
                    if deltaB > 0:
                        for sname,signal in signalHists.iteritems():
                            if debugLevel > 0: 
                                print "Subtracting",deltaB,"from signal histogram",sname,"due to contamination.  ",
                                print "Bin",bn,"content changes from",signal.GetBinContent(bn),
                            signal.SetBinContent( bn, max(0, signal.GetBinContent(bn) - deltaB) )
                            if debugLevel > 0:
                                print "to",signal.GetBinContent(bn)

def combineEmptyBinErrs(emptyBinErrs, combineBkgs):
    """
    Combine empty bin uncertainties on subprocesses sharing the same umbrella process.
    Inputs and outputs are both dictionaries of the form
    { process:{ bin1:err1, bin2:err2, ... }, ... }
    combineBkgs is of the form
    { proc1:[subproc1, ...], ... }
    """
    newEmptyBinErrs = { proc:{} for proc in combineBkgs }
    for proc, subprocs in combineBkgs.iteritems():
        for subproc in subprocs:
            for bx, err in emptyBinErrs[subproc].iteritems():
                curErr = newEmptyBinErrs[proc].get(bx, 0.0)
                newEmptyBinErrs[proc][bx] = (curErr*curErr + err*err)**(0.5)
    return newEmptyBinErrs

def combineBackgroundHists(hists, combineBackgrounds, listOfVars, debugLevel=0):
    """Combine many background histograms into one, for plotting purposes.
    combineBackgrounds should be a dictionary whose keys are the desired (combined) process names.
    The value for each key should be a list of backgrounds that should be combined."""
    for combProcess,combList in combineBackgrounds.iteritems():
        tmphists = {}
        if debugLevel > 0:
            print "Combining background histograms for",combProcess
        #check which processes are present
        combListPresent = filter(lambda c: c in hists, combList)
        if len(combListPresent) != len(combList):
            print "Warning in combineBackground hists:",(len(combList)-len(combListPresent)),"of",len(combList),"requested background processes are not present in the histogram collection"
            print "(Looking for %s, Found %s)"%(' '.join(combList),' '.join(combListPresent))
        if len(combListPresent) == 0: continue
        for v in listOfVars: #loop over variables
            #make a new histogram for the combined backgrounds
            combHist = hists[combListPresent[0]][v].Clone()
            combHist.SetName( combHist.GetName().replace(combListPresent[0], combProcess) )
            combHist.Reset()
            #add together the backgrounds
            for process in combListPresent:
                combHist.Add(hists[process][v])
                if debugLevel > 0:
                    print " Including",process
                #delete it from the dictionary
                #hists[process][v].Delete()
                del hists[process][v]
            #insert the new histogram
            tmphists[v] = combHist
        #clean up dictionary
        for process in combListPresent:
            del hists[process]
        hists[combProcess] = tmphists

def combineBackgroundNames(oldNames, combineBackgrounds):
    """Get list of sample names remaning after combining the specified backgrounds together"""
    processesToRemove = []
    processesToAdd = []
    for combProcess, combList in combineBackgrounds.iteritems():
        processesToAdd.append(combProcess)
        for proc in combList:
            processesToRemove.append(proc)
    newNames = filter((lambda x: x not in processesToRemove), oldNames)
    newNames = list(set( newNames + processesToAdd ))
    return newNames

def th1ToTGraphAsymmErrors(th1):
    """Convert a TH1 to a TGraphAsymmErrors with appropriate Poisson uncertainties."""
    #get x bin centers
    nbins = th1.GetNbinsX()
    xarray = array('d', [0 for x in range(nbins)])
    th1.GetXaxis().GetCenter(xarray)
    #get bin contents
    ylist = []
    xerrslist = []
    uperrslist = []
    downerrslist = []
    for bx in range(1, nbins+1):
        content = th1.GetBinContent(bx)
        ylist.append(content)
        xerrslist.append(0.5)
        uperrslist.append(getUpPoissonError(content))
        downerrslist.append(getDownPoissonError(content))
    yarray = array('d',ylist)
    xerrsarray = array('d',xerrslist)
    uperrsarray = array('d',uperrslist)
    downerrsarray = array('d',downerrslist)
    graph = rt.TGraphAsymmErrors(nbins, xarray, yarray, xerrsarray, xerrsarray, downerrsarray, uperrsarray)
    return graph

def getUpPoissonError(n):
    """Get upper range of Poisson 1-sigma confidence interval for a bin with n counts"""
    alpha = 1.- 0.682689492
    if n < 0:
        print "Error in getUpPoissonError: cannot get error for bin with negative contents!"
        return 0
    return rt.Math.gamma_quantile_c( alpha/2, n+1, 1 ) - n

def getDownPoissonError(n):
    """Get lower range of Poisson 1-sigma confidence interval for a bin with n counts"""
    alpha = 1.- 0.682689492
    if n < 0:
        print "Error in getDownPoissonError: cannot get error for bin with negative contents!"
        return 0
    if n == 0: 
        return 0
    return n - rt.Math.gamma_quantile( alpha/2, n, 1. )

def getBinEvidence(x, b, s):
    """Compute the contribution of bin with x observed, b background and s signal expected, to the poisson likelihood"""
    return -2 * (x * math.log( (s+b)*1.0/b ) - s)

def makeRazorMCTotalUnrolledHist(infile, samples, unrollRows, unrollCols, debugLevel=0):
    unrolledMC = {s:[] for s in samples}
    hists = importHists(infile, debugLevel)

    #make total background histogram
    unrolledMCs = plotting.unroll2DHistograms([hists[s][('MR','Rsq')] for s in samples], unrollRows, 
            unrollCols, labelBins=True)
    unrolledMCTotal = unrolledMCs[0].Clone()
    unrolledMCTotal.Reset()
    for n,s in enumerate(samples):
        unrolledMCTotal.Add(unrolledMCs[n])
    return unrolledMCTotal

def makeRazorMCTotalUnrolledHists(boxName, samples, inDir='.', 
        unrollBins=None, debugLevel=0, doData=False):
    """Retrieve MC histograms, unroll, and add together"""
    nbmax = 3
    if boxName == 'DiJet':
        nbmax = 2
    filenames = [inDir+"/razorHistograms"+boxName+str(b)+"BFineGrained.root" for b in range(nbmax+1)]
    unrolledMCs = []
    for i,f in enumerate(filenames):
        unrollRows = unrollBins[i][0]
        unrollCols = unrollBins[i][1]
        unrolledMCs.append( makeRazorMCTotalUnrolledHist(f, samples, unrollRows, unrollCols, debugLevel) )
    if doData:
        unrolledData = []
        for i, f in enumerate(filenames):
            unrollRows = unrollBins[i][0]
            unrollCols = unrollBins[i][1]
            unrolledData.append( makeRazorMCTotalUnrolledHist(
                f, ['Data'], unrollRows, unrollCols, debugLevel) )
        return unrolledMCs, unrolledData
    return unrolledMCs

def splitByUnrollBins(hist, unrollBins):
    """For splitting razor signal templates into histograms for individual b-tag bins"""
    outHists = []
    bn = 0
    for i, bins in enumerate(unrollBins):
        rows = bins[0]
        cols = bins[1]
        numBins = sum( [ len(col)-1 for col in cols ] )
        outHist = rt.TH1F(hist.GetName()+str(i), hist.GetTitle(), numBins, 0, numBins)
        outHist.SetDirectory(0)
        
        for bx in range(1, numBins+1):
            bn += 1
            outHist.SetBinContent(bx, hist.GetBinContent(bn))
            outHist.SetBinError(bx, hist.GetBinError(bn))
        outHists.append(outHist)
    return outHists

def loadScaleFactorHists(sfFilename="RazorScaleFactors.root", processNames=[], scaleFactorNames={}, debugLevel=0):
    """Returns a dict with available scale factor histograms"""
    if debugLevel > 0:
        print "Opening scale factor file",sfFilename
    sfFile = rt.TFile.Open(sfFilename)
    assert sfFile
    sfHists = {}
    for pname in processNames:
        histname = pname
        if pname in scaleFactorNames:
            histname = scaleFactorNames[pname]
        if debugLevel > 0: print "Looking for scale factor histogram",histname,"for",pname,"...",        

        tmp = sfFile.Get(histname+"ScaleFactors")
        if tmp: 
            sfHists[pname] = tmp
            sfHists[pname].SetDirectory(0)
            if debugLevel > 0: print "Found!"
        else:
            if debugLevel > 0: print ""

        #get up/down histograms
        tmp = sfFile.Get(histname+"ScaleFactorsUp")
        if tmp: 
            sfHists[pname+"Up"] = tmp
            sfHists[pname+"Up"].SetDirectory(0)
            if debugLevel > 0: print "Up histogram found!"
        tmp = sfFile.Get(histname+"ScaleFactorsDown")
        if tmp: 
            sfHists[pname+"Down"] = tmp
            sfHists[pname+"Down"].SetDirectory(0)
            if debugLevel > 0: print "Down histogram found!"

    sfFile.Close()
    return sfHists

def invertHistogram(hist):
    """Replaces contents of each hist bin with 1/(contents).  Updates bin errors accordingly.
       For bins with no contents, does nothing."""
    if not hist:
        raise ValueError("Hist does not exist")
    ret = hist.Clone(hist.GetName()+"Inverted")
    for b in range(hist.GetSize()+1):
        if hist.GetBinContent(b) != 0:
            ret.SetBinContent( b, 1.0/hist.GetBinContent(b) )
            ret.SetBinError( b, hist.GetBinError(b) / (hist.GetBinContent(b))**2 )
    return ret

def removeVarCuts(cuts, varName):
    """
    Removes all cuts on the given variable 
    from the given string.
    """
    while varName in cuts:
        split = cuts.split(' && ')
        cuts = ' && '.join([cut for cut in split if varName not in cut])
    return cuts

