from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
import os
import random
import sys
import math
from scipy.integrate import quad
from itertools import *
from operator import *

densityCorr = False

def densityCorrect(h):
    if "_densitycorr" in h.GetName():
        return h
    h_densitycorr = h.Clone(h.GetName()+"_densitycorr")
    h_densitycorr.Sumw2()
    for iBinX in range(1,h.GetNbinsX()+1):
        val = h.GetBinContent(iBinX)
        err = h.GetBinError(iBinX)
        width = h.GetBinWidth(iBinX)
        h_densitycorr.SetBinContent(iBinX,val/width)
        h_densitycorr.SetBinError(iBinX,err/width)
    return h_densitycorr

def convertSideband(name,w,x,y,z):
    if name=="Full":
        return "Full"
    names = name.split(',')
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    iBinX = -1
    sidebandBins = []
    for ix in range(1,len(x)):
        for iy in range(1,len(y)):
            for iz in range(1,len(z)):
                iBinX+=1
                w.var('MR').setVal((x[ix]+x[ix-1])/2.)
                w.var('Rsq').setVal((y[iy]+y[iy-1])/2.)
                w.var('nBtag').setVal((z[iz]+z[iz-1])/2.)
                inSideband = 0
                for fitname in names:
                    inSideband += ( w.var('MR').inRange(fitname) * w.var('Rsq').inRange(fitname) * w.var('nBtag').inRange(fitname) )
                if inSideband: sidebandBins.append(iBinX)

    sidebandGroups = []
    for k, g in groupby(enumerate(sidebandBins), lambda (i,x):i-x):
        consecutiveBins = map(itemgetter(1), g)
        sidebandGroups.append([consecutiveBins[0],consecutiveBins[-1]+1])
        
    newsidebands = ''
    nameNoComma = name.replace(',','')
        
    for iSideband, sidebandGroup in enumerate(sidebandGroups):
        if not w.var('th1x').hasRange('%s%i'%(nameNoComma,iSideband)):
            w.var('th1x').setRange("%s%i"%(nameNoComma,iSideband),sidebandGroup[0],sidebandGroup[1])
        newsidebands+='%s%i,'%(nameNoComma,iSideband)
    newsidebands = newsidebands[:-1]
    return newsidebands

def setStyle():
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetPaintTextFormat("1.2g")
    rt.gROOT.SetBatch()
    rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.FATAL)

def findLastBin(h):
    for i in range(1,h.GetXaxis().GetNbins()):
        thisbin = h.GetXaxis().GetNbins()-i
        if h.GetBinContent(thisbin)>=0.1: return thisbin+1
    return h.GetXaxis().GetNbins()

def useThisRho(minX,maxX,htemp):

    return 2.0

    if maxX<10:
        return 4.0
    elif maxX<40:
        return 3.0
    else:
        return 2.0

def getKDE(sumName,myTree,myHist):
    maxX = myHist.GetXaxis().GetBinUpEdge(myHist.GetNbinsX())
    varNames = sumName.split("+")
    nExp = [rt.RooRealVar(iVar,iVar,0,maxX) for iVar in varNames]
    nExpSet = rt.RooArgSet("nExpSet")
    nExpList = rt.RooArgList("nExpList")
    for expVar in nExp:
        nExpSet.add(expVar)
        nExpList.add(expVar)
    dataset = rt.RooDataSet("dataset","dataset",myTree,nExpSet)
    sumExp = rt.RooFormulaVar("sumExp","sumExp",sumName,nExpList)
    sumExpData = dataset.addColumn(sumExp)
    sumExpData.setRange(0,maxX)
    rho = useThisRho(0.,maxX,myHist)
    rkpdf = rt.RooKeysPdf("rkpdf","rkpdf",sumExpData,dataset,rt.RooKeysPdf.NoMirror,rho)
    func = rkpdf.asTF(rt.RooArgList(sumExpData))
    return func, rkpdf, dataset


def find68ProbRangeFromKDEMode(maxX,func, probVal=0.6827):
    mean = func.Mean(0.,maxX)
    funcMax = func.GetMaximum()
    mode = func.GetX(funcMax,0.,maxX)
    totalProb = func.Integral(0,maxX)
    nprobSum = int(1.0)
    probSum = array("d",[0.5])
    q = array("d",[0])
    func.GetQuantiles(nprobSum,q,probSum)
    median = q[0]
    # iterate first with a coarse epsilon
    # THEN iterate again with a smaller epsilon
    probRange = 0.

    if mode>0:
        epsilon=mode/3.0
    else:
        epsilon = 0.1
    #epsilon = max(mode/2, (maxX - mode)/2)
    above68 = False
    numIter = 0
    sigmaMinus = mode
    sigmaPlus = mode
    while abs(probRange - probVal)>0.01 and numIter < 100.:
        numIter += 1
        if probRange < probVal:
            if above68:  epsilon = epsilon/2.0
            above68 = False
            sigmaMinus = sigmaMinus - epsilon
            sigmaPlus = sigmaPlus + epsilon
        else:
            if not above68: epsilon = epsilon/2.0
            above68 = True
            sigmaMinus = sigmaMinus + epsilon
            sigmaPlus = sigmaPlus - epsilon
            
        if sigmaMinus<=0: sigmaMinus = 0.
        if sigmaPlus>=maxX: sigmaPlus = maxX
        if sigmaMinus<sigmaPlus : probRange = func.Integral(sigmaMinus,sigmaPlus)/totalProb
        else: probRange = 0.
        if numIter%5 == 0: print "iteration = %d"%numIter
    print "sigmaPlus = %f"%sigmaPlus
    print "sigmaMinus = %f"%sigmaMinus
    print "Int_[sigmaMinus,sigmaPlus] f(x) dx = %f"%probRange
   
    funcFill68 = func.Clone("funcFill68")
    funcFill68.SetRange(sigmaMinus,sigmaPlus)
    ic68 = rt.TColor(1404, 0.49, 0.60, 0.82,"")
    funcFill68.SetFillColor(ic68.GetColor(0.49,0.60,.82))
    funcFill68.SetFillStyle(3144)
    return mode,sigmaMinus,sigmaPlus,probRange,funcFill68


def getPValueFromKDE(nObs,maxX,func):
    funcObs = func.Eval(nObs)
    funcMax = func.GetMaximum()
    otherRoot = 0.
    rightSide = False
    veryNearMax = False
    epsilon = max(0.2,nObs/100)
    pvalKDE = 0
    totalProb = func.Integral(0,maxX)
    if abs(float(funcObs)-funcMax)/funcMax < 0.003:
        veryNearMax = True
        pvalKDE = 1.
        otherRoot = nObs
    elif func.Derivative(nObs)<0.:
        rightSide = True
        otherRoot = func.GetX(funcObs,0,nObs-epsilon)
        if otherRoot >= nObs-epsilon:
            otherRoot = func.GetX(funcObs*(0.99),0,nObs-epsilon)
        elif otherRoot < nObs-epsilon:
            print "contribution from 0 to otherRoot =", func.Integral(0,otherRoot)
            pvalKDE = func.Integral(0,otherRoot)
        print "contribution from n0bs to maxX =", func.Integral(nObs,maxX)
        pvalKDE += func.Integral(nObs,maxX)
    else:
        otherRoot = func.GetX(funcObs,nObs+epsilon,maxX)
        if otherRoot <= nObs+epsilon:
            otherRoot = func.GetX(funcObs*(0.99),nObs+epsilon,maxX)
        elif otherRoot > nObs+epsilon:
            print "contribution from 0 to nObs =", func.Integral(0,nObs)
            pvalKDE = func.Integral(0,nObs)
        print "contribution from otherRoot to maxX =", func.Integral(otherRoot,maxX)
        pvalKDE += func.Integral(otherRoot,maxX)
    pvalKDE = pvalKDE/totalProb
    if pvalKDE>1.:
        pvalKDE = 1.
    # DRAWING FUNCTION AND FILLS
    ic = rt.TColor(1398, 0.75, 0.92, 0.68,"")
    func.SetLineColor(ic.GetColor(0.1, .85, 0.5))
    funcFillRight = func.Clone("funcFillRight")
    funcFillLeft = func.Clone("funcFillLeft")
    if veryNearMax:
        funcFillRight.SetRange(nObs,maxX)
        funcFillLeft.SetRange(0,nObs)
    elif rightSide:
        funcFillRight.SetRange(nObs,maxX)
        funcFillLeft.SetRange(0,otherRoot)
    else:
        funcFillRight.SetRange(otherRoot,maxX)
        funcFillLeft.SetRange(0,nObs)
    funcFillRight.SetFillColor(ic.GetColor(0.1, .85, 0.5))
    funcFillRight.SetFillStyle(3002)
    funcFillLeft.SetFillColor(ic.GetColor(0.1, .85, 0.5))
    funcFillLeft.SetFillStyle(3002)
    func.Draw("same")
    funcFillRight.Draw("fcsame")
    drawLeft = (rightSide and otherRoot < nObs-epsilon) or not rightSide
    # PRINTING INFORMATION
    print "nObs =  %d"%(nObs)
    print "f(nObs) =  %f"%(funcObs)
    print "fMax = %f"%(funcMax)
    print "percent diff = %f"%(abs(float(funcObs)-funcMax)/funcMax)
    print "other root = %f"%(otherRoot)
    print "total prob =  %f"%(totalProb)
    print "pvalKDE = %f"%pvalKDE
    return pvalKDE,funcFillRight,funcFillLeft,drawLeft


def find68ProbRange(hToy, probVal=0.6827):
    minVal = 0.
    maxVal = 100000.
    mode = hToy.GetBinLowEdge(hToy.GetMaximumBin())
    integral = hToy.Integral()
    if integral<=0: 
        return mode,minVal,maxVal
    
    ### range calculation assumes a unimodal toy distribution

    #print "nbinsX      %i"%(hToy.GetNbinsX()) 
    #print "nbinsXaxis  %i"%(hToy.GetXaxis().GetNbins())
    #get a good number of histogram bins 
    targetNumBins = 100
    rebinBy = int(hToy.GetNbinsX()/targetNumBins)
    if rebinBy > 1: hToy = hToy.Rebin(rebinBy,hToy.GetName()+"rebinned")
    mode = hToy.GetBinLowEdge(hToy.GetMaximumBin()) #recompute mode
    #print "nbinsX      %i"%(hToy.GetNbinsX()) 
    #print "nbinsXaxis  %i"%(hToy.GetXaxis().GetNbins())

    #start at the mode
    leftBin = hToy.GetMaximumBin()
    rightBin = hToy.GetMaximumBin()
    curProb = hToy.GetBinContent(leftBin)/integral

    #edge case: the best fit bin has enough probability already
    if curProb >= probVal:
        binFractionToTake = probVal/curProb
        minVal = hToy.GetXaxis().GetBinCenter(leftBin) - hToy.GetXaxis().GetBinWidth(leftBin)*binFractionToTake/2
        maxVal = hToy.GetXaxis().GetBinCenter(leftBin) + hToy.GetXaxis().GetBinWidth(leftBin)*binFractionToTake/2

    #add bins on either side of the best fit bin until probVal is reached
    #print "before entering while loop"
    while curProb < probVal:
        #print "starting while loop, curProb = %f < %f"%(curProb, probVal)
        if leftBin <= 1 and rightBin >= hToy.GetNbinsX():
            print "Error in find68ProbRange: failed to reach the target probability value",probVal
            return mode, minVal, maxVal

        #get probabilities to the left and right
        if leftBin-1 > 0:
            leftProb = hToy.GetBinContent(leftBin-1)/integral
        else:
            leftProb = 0.0
        if rightBin+1 < hToy.GetNbinsX()+1:
            rightProb = hToy.GetBinContent(rightBin+1)/integral
        else:
            rightProb = 0.0

        #decide which way to go
        if leftProb > rightProb or (rightProb<=0.0 and leftBin-1 > 0): 
            #go left
            #print "going left, leftBin = %i"%leftBin
            leftBin -= 1
            if curProb + leftProb >= probVal:
                #get left edge (correct for overshooting)
                binFractionToTake = (probVal-curProb)/leftProb
                minVal = hToy.GetXaxis().GetBinUpEdge(leftBin) - hToy.GetXaxis().GetBinWidth(leftBin)*binFractionToTake
                maxVal = hToy.GetXaxis().GetBinUpEdge(rightBin)
            #update probability sum
            curProb += leftProb
            #print "update curProb = %f with leftProb = %f" %(curProb,leftProb)
        else: 
            #go right
            #print "going right, rightBin = %i"%rightBin
            rightBin += 1
            if curProb + rightProb >= probVal:
                #get right edge (correct for overshooting)
                binFractionToTake = (probVal-curProb)/rightProb
                minVal = hToy.GetXaxis().GetBinLowEdge(leftBin)
                maxVal = hToy.GetXaxis().GetBinLowEdge(rightBin) + hToy.GetXaxis().GetBinWidth(rightBin)*binFractionToTake
            #update probability sum
            curProb += rightProb
            #print "update curProb = %f with rightProb = %f" %(curProb,rightProb)

    return mode, minVal,maxVal

def findMedian(myHisto):
    prob = 0
    median = 0
    for iBin in range(1, myHisto.GetNbinsX()+1):
        if prob <= 0.5 and prob+myHisto.GetBinContent(iBin) > 0.5:
            median = myHisto.GetBinCenter(iBin)
        prob = prob + myHisto.GetBinContent(iBin)
    return median

def getPValue(n, hToy):
    if hToy.Integral() <= 0.: return 0.
    Prob_n = hToy.GetBinContent(hToy.FindBin(n))
    Prob = 0
    for i in range(1, hToy.GetNbinsX()+1):
        if hToy.GetBinContent(i)<= Prob_n: Prob += hToy.GetBinContent(i)
    Prob = Prob/hToy.Integral()
    return Prob

def getSigmaFromPval(n, bf, hToy, pVal):
    medianVal = findMedian(hToy)
    if n>bf: return rt.TMath.NormQuantile(1 - pVal/2.)
    else: return rt.TMath.NormQuantile( pVal/2. )

        
def getBinSumDicts(sumType, minX, maxX, minY, maxY, minZ, maxZ, x, y, z):

    binSumDict = {}

    if sumType.lower()=="x":
        for i in range(1,len(x)):
            binSumDict[i] = ''
    elif sumType.lower()=="y":
        for j in range(1,len(y)):
            binSumDict[j] = ''
    elif sumType.lower()=="z":
        for k in range(1,len(z)):
            binSumDict[k] = ''
    elif sumType.lower()=="yx":
        for i in range(1,len(x)):
            for j in range(1,len(y)):
                binSumDict[i,j] = ''
    elif sumType.lower()=="zyx":
        for i in range(1,len(x)):
            for j in range(1,len(y)):
                for k in range(1,len(z)):
                    binSumDict[i,j,k] = ''
            
    iBinX = -1
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            for k in range(1,len(z)):
                iBinX += 1
                if i>=minX and i<=maxX and j>=minY and j<=maxY and k>=minZ and k<=maxZ:
                    if sumType.lower()=="x":
                        binSumDict[i] += 'b%i+'%iBinX
                    elif sumType.lower()=="y":
                        binSumDict[j] += 'b%i+'%iBinX
                    elif sumType.lower()=="z":
                        binSumDict[k] += 'b%i+'%iBinX
                    elif sumType.lower()=="yx":
                        binSumDict[i,j] += 'b%i+'%iBinX
                    elif sumType.lower()=="zyx":
                        binSumDict[i,j,k] = 'b%i'%iBinX
                
    if sumType.lower()=="x":
        for i in range(1,len(x)):
            binSumDict[i] = binSumDict[i][0:-1]            
    elif sumType.lower()=="y":
        for j in range(1,len(y)):
            binSumDict[j] = binSumDict[j][0:-1]          
    elif sumType.lower()=="z":
        for k in range(1,len(z)):
            binSumDict[k] = binSumDict[k][0:-1]     
    elif sumType.lower()=="yx":
        for i in range(1,len(x)):
            for j in range(1,len(y)):
                binSumDict[i,j] = binSumDict[i,j][0:-1]
            
    return binSumDict

def getCorrelationCoefficient(myTree, sumName1, sumName2):
    #get mean and standard deviation of each quantity
    c = rt.TCanvas('empty','empty',400,300)
    myTree.Draw('%s>>htest%s'%(sumName1,myTree.GetName()+sumName1.replace('+','')))
    htemp1 = rt.gPad.GetPrimitive("htest%s"%(myTree.GetName()+sumName1.replace('+','')))
    mean1 = htemp1.GetMean()
    sigma1 = htemp1.GetRMS()
    myTree.Draw('%s>>htest%s'%(sumName2,myTree.GetName()+sumName2.replace('+','')))
    htemp2 = rt.gPad.GetPrimitive("htest%s"%(myTree.GetName()+sumName2.replace('+','')))
    mean2 = htemp2.GetMean()
    sigma2 = htemp2.GetRMS()

    #get covariance
    myTree.Draw('(%s - %f)*(%s - %f)>>hcov%s%s'%(sumName1,mean1,sumName2,mean2,myTree.GetName()+sumName1.replace('+',''),myTree.GetName()+sumName2.replace('+','')))
    htempCov = rt.gPad.GetPrimitive('hcov%s%s'%(myTree.GetName()+sumName1.replace('+',''),myTree.GetName()+sumName2.replace('+','')))
    covariance = htempCov.GetMean()

    #compute correlation coefficient
    corrCoeff = covariance/(sigma1*sigma2)
    return corrCoeff

def getCorrelationMatrix(myTree, sumType, minX, maxX, minY, maxY, minZ, maxZ, x, y, z):
    binSumDict = getBinSumDicts(sumType, minX, maxX, minY, maxY, minZ, maxZ, x, y, z)

    #histogram for output
    nbins = (maxX-minX)*(maxY-minY)*(maxZ-minZ)
    h = rt.TH2F("correlationMatrix", "correlationMatrix",nbins,0,nbins,nbins,0,nbins)
    h.SetDirectory(0)

    i = 0
    for ix1 in range(minX, maxX):
        for iy1 in range(minY, maxY):
            for iz1 in range(minZ, maxZ):
                j = 0
                i += 1
                sumName1 = binSumDict[(ix1+1,iy1+1,iz1+1)]
                for ix2 in range(minX, maxX):
                    for iy2 in range(minY, maxY):
                        for iz2 in range(minZ, maxZ):
                            j += 1
                            if j > i: break
                            sumName2 = binSumDict[(ix2+1,iy2+1,iz2+1)]
                            corrCoeff = getCorrelationCoefficient(myTree, sumName1, sumName2)
                            h.SetBinContent(i,j,corrCoeff)
                            h.SetBinContent(j,i,corrCoeff)
                                
    return h

def getBestFitRms(myTree, sumName, nObs, d, options, plotName):
    myTree.GetEntry(0)
    bestFit = eval(sumName.replace('b','myTree.b'))
    myTree.Draw('%s>>htest%s'%(sumName,myTree.GetName()+sumName.replace('+','')))
    htemp = rt.gPad.GetPrimitive("htest%s"%(myTree.GetName()+sumName.replace('+','')))
    #myTree.Draw('%s>>htest'%(sumName))    
    #htemp = rt.gPad.GetPrimitive("htest")
    xmax = int(max(htemp.GetXaxis().GetXmax()+1,nObs+1))
    xmin = int(max(0,htemp.GetXaxis().GetXmin()))
    mean = htemp.GetMean()
    
    # just make 20 bins if (xmax-xmin)>40 and mean>10 (because we're using KDE)
    useKDE = ((xmax-xmin)>40) and (mean>10)
    useKDE = False
    if useKDE:
        #print "useKDE"
        if (xmax-xmin)%20:
            xmax = xmax + 20-(xmax-xmin)%20
        nbins=20
    #no-stat bins should be much smaller (to avoid having too many nsigma=0 bins)
    elif options.noStat and ('nsigma' not in plotName):
        #print "noStat and noNsigma"
        nbins = (xmax-xmin)*50
        #nbins = (xmax-xmin)
    else:
        #print "else"
        nbins = xmax-xmin

    myTree.Draw('%s>>htemp%s(%i,%f,%f)'%(sumName,myTree.GetName()+sumName.replace('+',''),nbins,xmin,xmax))
    htemp = rt.gPad.GetPrimitive("htemp%s"%(myTree.GetName()+sumName.replace('+','')))
    #myTree.Draw('%s>>htemp(%i,%f,%f)'%(sumName,nbins,xmin,xmax))
    #htemp = rt.gPad.GetPrimitive("htemp")
    mean = htemp.GetMean()
    probRange = 0.6827
    mode,rangeMin,rangeMax = find68ProbRange(htemp,probRange)
    range68 = (rangeMax-rangeMin)
    rms = htemp.GetRMS()
    pvalue = getPValue(nObs,htemp)    
    if useKDE:        
        func, rkpdf, dataset = getKDE(sumName,myTree,htemp)
        pvalue,funcFillRight,funcFillLeft,drawLeft = getPValueFromKDE(nObs,xmax,func)
        mode,rangeMin,rangeMax,probRange,funcFill68 = find68ProbRangeFromKDEMode(xmax,func)
        range68 = (rangeMax-rangeMin)
        
    nsigma = getSigmaFromPval(nObs, bestFit, htemp, pvalue)
    
    if pvalue <= 0.:
        print "pvalue = 0 from histogram method: reverting to gaussian approximation"
        #print rms
        if rms<=0:
            print "rms==0, this is sad"
            nsigma = 100
        else:
            nsigma = (nObs-bestFit)/rms
        pvalue = 2.*rt.Math.gaussian_cdf_c(abs(nsigma))
    print '%s, bestFit %f, mean %.1f, mode %.1f, rms %.1f, pvalue %f, nsigma %.1f, range68 %.1f, %.1f'%(sumName, bestFit,mean,mode,rms,pvalue,nsigma,rangeMax,rangeMin)

    if options.printErrors:
        if htemp.Integral()>0.: 
            htemp.Scale(1./htemp.Integral()/htemp.GetBinWidth(1))
        htemp.SetMinimum(0)
        htemp.SetMaximum(1.25*(htemp.GetMaximum()+htemp.GetBinError(htemp.GetMaximumBin())))
        htemp.SetMarkerStyle(20)
        htemp.SetMarkerColor(rt.kBlack)
        htemp.SetLineColor(rt.kBlack)
        tgraph = rt.TGraph(4)
        tgraph.SetPoint(1, rangeMin,0)
        tgraph.SetPoint(2, rangeMin,htemp.GetMaximum())
        tgraph.SetPoint(3, rangeMax,htemp.GetMaximum())
        tgraph.SetPoint(4, rangeMax,0)
        tgraph.SetFillColor(rt.kBlue-10)
        tgraph.SetLineColor(rt.kBlue-10)
        tline = rt.TLine(bestFit,0,bestFit,htemp.GetMaximum())
        tline.SetLineColor(rt.kBlue)
        tline.SetLineWidth(2)
        tlineObs = rt.TLine(nObs,0,nObs,htemp.GetMaximum())
        tlineObs.SetLineColor(rt.kBlack)
        tlineObs.SetLineWidth(2)
        htemp.Draw("pe")
        if useKDE:
            func.Draw("same")
            #funcFill68.Draw("fcsame")
            funcFillRight.Draw("fcsame")
            if drawLeft: funcFillLeft.Draw("fcsame")
        else:
            tgraph.Draw("fsame")
            
        tline.Draw("same")
        tlineObs.Draw("same")
        htemp.Draw("pesame")

        
        tleg = rt.TLegend(0.65,.6,.89,.89)
        tleg.AddEntry(tlineObs,"observed = %.2f"%nObs,"l")
        tleg.AddEntry(tline,"best fit = %.2f"%bestFit,"l")
        #tleg.AddEntry(tgraph,"rms = %.2f"%rms,"f")
        if useKDE:
            tleg.AddEntry(funcFillRight,"p-value = %.4f"%pvalue,"f")
            tleg.AddEntry(funcFill68,"%.1f%% range = [%.1f,%.1f]"%(probRange*100,rangeMin,rangeMax),"f")
        else:
            tleg.AddEntry(tlineObs,"p-value = %.4f"%pvalue,"f")
            tleg.AddEntry(tgraph,"%.1f%% range = [%.1f,%.1f]"%(probRange*100,rangeMin,rangeMax),"f")
            
        tleg.SetFillColor(rt.kWhite)
        tleg.SetLineColor(rt.kWhite)
        tleg.Draw("same")
                
        d.SetLogy(0)
        d.SetLogx(0)
        os.system("mkdir -p "+options.outDir+"/errors")
        d.Print(options.outDir+"/errors/"+plotName)
        d.Print(options.outDir+"/errors/"+plotName.replace(".pdf",".C"))
    return bestFit, range68/2.0, pvalue, nsigma, d

def getPads(c):    
    pad1 = rt.TPad(c.GetName()+"_pad1","pad1",0,0.25,1,1)
    pad2 = rt.TPad(c.GetName()+"_pad2","pad2",0,0,1,0.25)
    pad1.Range(-213.4588,-0.3237935,4222.803,5.412602);
    pad2.Range(-213.4588,-2.206896,4222.803,3.241379);
    pad1.SetLeftMargin(0.15)
    pad2.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad2.SetRightMargin(0.05)
    pad1.SetTopMargin(0.12)
    pad2.SetTopMargin(0.)
    pad1.SetBottomMargin(0.)
    pad2.SetBottomMargin(0.47)
    pad1.Draw()
    pad1.cd()
    rt.gPad.SetLogy(1)
    return pad1, pad2

def getPadsNs(c):    
    pad1 = rt.TPad(c.GetName()+"_pad1","pad1",0,0.25,1,1)
    pad2 = rt.TPad(c.GetName()+"_pad2","pad2",0,0,1,0.25)
    pad2.SetGridy(1)
    pad1.Range(-213.4588,-0.3237935,4222.803,5.412602);
    pad2.Range(-213.4588,-2.206896,4222.803,3.241379);
    pad1.SetLeftMargin(0.15)
    pad2.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad2.SetRightMargin(0.05)
    pad1.SetTopMargin(0.12)
    pad2.SetTopMargin(0.)
    pad1.SetBottomMargin(0.)
    pad2.SetBottomMargin(0.47)
    pad1.Draw()
    pad1.cd()
    rt.gPad.SetLogy(1)
    return pad1, pad2

def setDataHist(h_data,xTitle,yTitle,densityCorr=False,color=rt.kBlack):
    h_data.SetMarkerColor(color)
    h_data.SetMarkerStyle(20)
    h_data.SetLineColor(color)
    h_data.GetXaxis().SetTitle(xTitle)
    h_data.GetYaxis().SetTitle(yTitle)
    h_data.GetXaxis().SetLabelOffset(0.16)
    h_data.GetXaxis().SetLabelSize(0.06)
    h_data.GetYaxis().SetLabelSize(0.06)
    h_data.GetXaxis().SetTitleSize(0.06)
    h_data.GetYaxis().SetTitleSize(0.08)
    h_data.GetXaxis().SetTitleOffset(0.8)
    h_data.GetYaxis().SetTitleOffset(0.7)
    h_data.GetXaxis().SetTicks("+-")
    if "h_Rsq_" in h_data.GetName() and "MR" in h_data.GetName():
        h_data.SetMaximum(max(10,math.pow(h_data.GetBinContent(h_data.GetMaximumBin()),1.5)))
        h_data.SetMinimum(5e-2)
    else:        
        h_data.SetMaximum(max(10,math.pow(h_data.GetBinContent(h_data.GetMaximumBin()),1.25)))
        #h_data.SetMinimum(max(1e-1,1e-1*h_data.GetBinContent(h_data.GetMinimumBin())))
        if densityCorr and "_MR_" in h_data.GetName():
            h_data.SetMinimum(1e-6)
        elif '_Rsq_' in h_data.GetName() or 'nBtagRsqMR_y' in h_data.GetName():
            h_data.SetMinimum(5e-2) 
        elif 'MR' in h_data.GetName() or 'nBtagRsqMR_x' in h_data.GetName():
            h_data.SetMinimum(5e-4)
        else:
            h_data.SetMinimum(5e-4) # for th1x
    return h_data

def getDivideHistos(h,hClone,h_data,xTitle,divTitle):
    h.Sumw2()
    hClone.Sumw2()
    h_data.Sumw2()

    hDivide = h.Clone(h.GetName()+"Divide") 
    hCloneDivide = hClone.Clone(hClone.GetName()+"Divide") 
    hDataDivide = h_data.Clone(h_data.GetName()+"Divide")
    hDivide.Sumw2()
    hCloneDivide.Sumw2()
    hDataDivide.Sumw2()
    hCloneDivide.GetYaxis().SetLabelSize(0.18)
    hCloneDivide.SetTitle("")
    hCloneDivide.SetMaximum(3.5)
    hCloneDivide.SetMinimum(0.)
    hCloneDivide.GetXaxis().SetLabelSize(0.22)
    hCloneDivide.GetXaxis().SetTitleSize(0.22)

    
    for i in range(1, h_data.GetNbinsX()+1):
        tmpVal = hCloneDivide.GetBinContent(i)
        if tmpVal != -0.:
            hDataDivide.SetBinContent(i, hDataDivide.GetBinContent(i)/tmpVal)
            hDataDivide.SetBinError(i, hDataDivide.GetBinError(i)/tmpVal)
            hCloneDivide.SetBinContent(i, hCloneDivide.GetBinContent(i)/tmpVal)
            hCloneDivide.SetBinError(i, hCloneDivide.GetBinError(i)/tmpVal)
            hDivide.SetBinContent(i, hDivide.GetBinContent(i)/tmpVal)
            hDivide.SetBinError(i, hDivide.GetBinError(i)/tmpVal)

            
    hCloneDivide.GetXaxis().SetTitleOffset(0.97)
    hCloneDivide.GetXaxis().SetLabelOffset(0.02)
    hCloneDivide.GetXaxis().SetTitle(xTitle)

    hCloneDivide.GetYaxis().SetNdivisions(504,rt.kTRUE)
    hCloneDivide.GetYaxis().SetTitleOffset(0.2)
    hCloneDivide.GetYaxis().SetTitleSize(0.22)
    hCloneDivide.GetYaxis().SetTitle(divTitle)
    hCloneDivide.GetXaxis().SetTicks("+")
    hCloneDivide.GetXaxis().SetTickLength(0.07)
    hCloneDivide.SetMarkerColor(hCloneDivide.GetFillColor())
    
    return hDivide, hCloneDivide, hDataDivide
    
def print1DProj(c,rootFile,h,h_data,printName,xTitle,yTitle,lumiLabel="",boxLabel="",plotLabel="",isData=False,doSignalInj=False,options=None,tLeg=None,h_components=[],h_colors=[],h_labels=[]):
    
    if densityCorr:
        h_densitycorr = densityCorrect(h)
        h_data_densitycorr = densityCorrect(h_data)        
        h_components_densitycorr = [densityCorrect(h_comp) for h_comp in h_components]
        
        h = h_densitycorr
        h_data = h_data_densitycorr
        h_components = h_components_densitycorr
        
    pad1, pad2 = getPads(c)

    h.SetLineWidth(2)
    h.SetLineColor(rt.kBlue)
    hClone = h.Clone(h.GetName()+"Clone")
    hClone.SetLineColor(rt.kBlue)
    hClone.SetFillColor(rt.kBlue-10)
    
    h_data = setDataHist(h_data,xTitle,yTitle,densityCorr)

    if 'th1x' in h_data.GetName() and '_MultiJet' in h_data.GetName() and h_data.GetNbinsX()>140:
        h_data.GetXaxis().SetRange(0,140)
    
    if doSignalInj:
        h_data.SetMaximum(pow(h_data.GetBinContent(h_data.GetMaximumBin()),1.75))
        if h_data.GetBinContent(h_data.GetMaximumBin()) < 15:
            h_data.SetMaximum(225)
    h_data.Draw("pe")
    hClone.Draw("e2same")
    h.SetFillStyle(0)
    for h_comp, color, label in zip(h_components, h_colors, h_labels):
        h_comp.SetLineColor(color)
        h_comp.SetLineWidth(2)            
        h_comp.Draw("histsame")
        
    if doSignalInj:
        c4 = rt.gROOT.GetColor(rt.kGray+2)
        c4.SetAlpha(1.0)
        #signal is always last component
        h_components[-1].SetLineColor(rt.kBlack)
        h_components[-1].SetFillColor(rt.kGray+2)
        h_components[-1].SetLineStyle(2)
        h_components[-1].SetFillStyle(3005)
        h_components[-1].SetLineWidth(2)
        h_components[-1].Draw("histfsame")
    h.DrawCopy("histsame")
    h_data.Draw("pesame")
    pad1.Draw()
    c.Update()
    c.cd()
    pad2.Draw()
    pad2.cd()
    rt.gPad.SetLogy(0)

    hDivide, hCloneDivide, hDataDivide  = getDivideHistos(h, hClone, h_data, xTitle, "Data/Fit")
    
    if 'th1x' in hCloneDivide.GetName() and '_MultiJet' in hCloneDivide.GetName() and hCloneDivide.GetNbinsX()>140:
        hCloneDivide.GetXaxis().SetRange(0,140)
        
    hCloneDivide.Draw("e2")
    #hDivide.Draw("histsame")
    hDataDivide.Draw('pesame')
    hCloneDivide.Draw("axissame")


    pad2.Update()
    pad1.cd()
    pad1.Update()
    pad1.Draw()

    if tLeg==None:
        if len(h_components)>=7 or (doSignalInj and len(h_components)>=5):
            tLeg = rt.TLegend(0.7,0.3,0.9,0.8)
        elif len(h_components)==6 or (doSignalInj and len(h_components)==4):
            tLeg = rt.TLegend(0.7,0.35,0.9,0.8)
        elif len(h_components)==5 or (doSignalInj and len(h_components)==3):
            tLeg = rt.TLegend(0.7,0.4,0.9,0.8)
        elif len(h_components)==4 or (doSignalInj and len(h_components)==2):
            tLeg = rt.TLegend(0.7,0.45,0.9,0.8)
        elif len(h_components)==3 or (doSignalInj and len(h_components)==1):
            tLeg = rt.TLegend(0.7,0.5,0.9,0.8)
        elif len(h_components)==2 or (doSignalInj and len(h_components)==0):
            tLeg = rt.TLegend(0.7,0.55,0.9,0.8)            
        else:
            tLeg = rt.TLegend(0.7,0.6,0.9,0.8)
        tLeg.SetFillColor(0)
        tLeg.SetTextFont(42)
        tLeg.SetLineColor(0)
        if isData:
            tLeg.AddEntry(h_data,"Data","lep")
        else:
            tLeg.AddEntry(h_data,"Sim. Data","lep")
        tLeg.AddEntry(hClone,"Fit Total","lf")
        for h_comp, color, label in zip(h_components, h_colors, h_labels):                
            tLeg.AddEntry(h_comp,label,"l")
        if doSignalInj:
            if options.model=="T1bbbb":
                tLeg.AddEntry(h_components[-1],'pp#rightarrow#tilde{g}#tilde{g}, #mu = %.1f'%options.r,'lf')
                tLeg.AddEntry(None,'#tilde{g}#rightarrowb#bar{b}#tilde{#chi}_{1}^{0}','')
            
    tLeg.Draw("same")

    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    if isData:
        l.DrawLatex(0.15,0.9,"CMS preliminary")
    else:
        l.DrawLatex(0.15,0.9,"CMS simulation")
    l.DrawLatex(0.75,0.9,"%s"%lumiLabel)
    l.SetTextFont(52)
    l.SetTextSize(0.045)
    l.DrawLatex(0.2,0.82,boxLabel)
    l.DrawLatex(0.3,0.77,plotLabel)
    if doSignalInj:        
        if options.model=="T1bbbb":
            l.SetTextFont(42)
            l.DrawLatex(0.3,0.77,"m_{#tilde{g}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(options.mGluino,options.mLSP))

    c.cd()
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')
    cWrite = c.Clone(os.path.splitext(printName)[0].split('/')[-1])
    rootFile.cd()
    c.Write(os.path.splitext(printName)[0].split('/')[-1])

    
def print1DProjNs(c,rootFile,h,h_data,h_ns,printName,xTitle,yTitle,lumiLabel="",boxLabel="",plotLabel="",isData=False,doSignalInj=False,options=None,tLeg=None,h_components=[],h_colors=[],h_labels=[],cfg=None):
    
    if densityCorr:
        h_densitycorr = densityCorrect(h)
        h_data_densitycorr = densityCorrect(h_data)        
        h_components_densitycorr = [densityCorrect(h_comp) for h_comp in h_components]
        
        h = h_densitycorr
        h_data = h_data_densitycorr
        h_components = h_components_densitycorr
        
    pad1, pad2 = getPadsNs(c)

    h.SetLineWidth(2)
    h.SetLineColor(rt.kBlue)
    hClone = h.Clone(h.GetName()+"Clone")
    hClone.SetLineColor(rt.kBlue)
    hClone.SetFillColor(rt.kBlue-10)
    
    h_data = setDataHist(h_data,xTitle,yTitle,densityCorr)

    if 'th1x' in h_data.GetName() and '_MultiJet' in h_data.GetName() and h_data.GetNbinsX()>140:
        h_data.GetXaxis().SetRange(0,140)

    if doSignalInj:
        h_data.SetMaximum(pow(h_data.GetBinContent(h_data.GetMaximumBin()),1.75))
        if h_data.GetBinContent(h_data.GetMaximumBin()) < 15:
            h_data.SetMaximum(225)
    h_data.Draw("pe")
    hClone.Draw("e2same")

    
    nBinsMR = 8
    #nBinsMR = 6
    #if 'MuMultiJet' in boxLabel or 'EleMultiJet' in boxLabel or 'LeptonMultiJet' in boxLabel or 'EleJet' in boxLabel or 'MuJet' in boxLabel or 'LeptonJet' in boxLabel:
    #    nBinsMR = 7
    
    
    if 'th1x' in h.GetName():
        if options.fitRegion=="LowMR,LowRsq":
            hFit = hClone.Clone("h_fitregion")
            for i in range(1,h_data.GetNbinsX()+1):
                if i>=nBinsMR+1:
                    if (i-1)%(nBinsMR)!=0:
                        hFit.SetBinContent(i, 0)
            hFit.SetFillColor(rt.kGreen-4)
            hFit.SetFillStyle(3144)
            hFit.Draw("e2same")
                
    h.SetFillStyle(0)
    for h_comp, color, label in zip(h_components, h_colors, h_labels):
        h_comp.SetLineColor(color)
        h_comp.SetLineWidth(2)            
        h_comp.Draw("histsame")

    if doSignalInj:
        c4 = rt.gROOT.GetColor(rt.kGray+2)
        c4.SetAlpha(1.0)
        #signal is always last component
        h_components[-1].SetLineColor(rt.kBlack)
        h_components[-1].SetFillColor(rt.kGray+2)
        h_components[-1].SetLineStyle(2)
        h_components[-1].SetFillStyle(3005)
        h_components[-1].SetLineWidth(2)
        h_components[-1].Draw("histfsame")
    h.DrawCopy("histsame")
    h_data.Draw("pesame")

        
    pad1.Draw()
    c.Update()
    
    if 'th1x' in h.GetName():
        tlines = []
        for i in range(1,h_data.GetNbinsX()):
            if i%nBinsMR==0:
                #tlines.append(rt.TLine(i, 0, i, h_data.GetBinContent(h_data.GetMaximumBin())))
                ymax = pad1.GetUymax()
                ymin = pad1.GetUymin()
                newYmax = (ymax-ymin)/3.+ymin
                tlines.append(rt.TLine(i, pow(10,ymin), i, pow(10,newYmax)))
        for tline in tlines:
            #tline.SetNDC()
            tline.SetLineStyle(2)
            tline.Draw()

    c.cd()
    pad2.Draw()
    pad2.cd()
    rt.gPad.SetLogy(0)

    
    hDivide, hCloneDivide, hDataDivide  = getDivideHistos(h, hClone, h_data, xTitle, "Stat.+Sys. n#sigma")    
    if cfg is not None:
        box = boxLabel.split(' ')[1]
        x = array('d', cfg.getBinning(box)[0]) # MR binning
        y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    hDataDivideNs = get1DHistoFrom2D(h_ns,x,y,h_ns.GetName()+'1d')
    hDataDivideNs.SetMarkerStyle(8)
    hDataDivideNs.SetMarkerSize(0.75)
    #hDataDivideNs.SetMarkerStyle(20)
    hCloneDivide.SetMaximum(6)
    hCloneDivide.SetMinimum(-6)    
    hCloneDivide.GetYaxis().SetNdivisions(503,True)
    hCloneDivide.GetYaxis().SetTitleSize(0.16)
    hCloneDivide.GetYaxis().SetTitleOffset(0.25)
    
    hTwoSigma = hCloneDivide.Clone("hTwoSigma")
    hDivide.SetLineColor(rt.kBlue)
    hCloneDivide.SetFillColor(rt.kYellow)
    hTwoSigma.SetFillColor(rt.kGreen)
    for i in range(1,hCloneDivide.GetNbinsX()+1):        
        hCloneDivide.SetBinContent(i,0)
        hCloneDivide.SetBinError(i,1)
        hTwoSigma.SetBinContent(i,0)
        hTwoSigma.SetBinError(i,2)        
        hDivide.SetBinContent(i,0)
        hDivide.SetBinError(i,0)
        hDataDivide.SetBinContent(i,hDataDivideNs.GetBinContent(i))
        hDataDivide.SetBinError(i,0)
        hDataDivideNs.SetBinError(i,0)
    
    if 'th1x' in hCloneDivide.GetName() and '_MultiJet' in hCloneDivide.GetName() and hCloneDivide.GetNbinsX()>140:
        hCloneDivide.GetXaxis().SetRange(0,140)
        
    hTwoSigma.Draw("e2")
    hCloneDivide.Draw("e2same")
    hDivide.Draw("lsame")
    #hDivide.Draw("histsame")
    hDataDivideNs.Draw('pesame')
    hCloneDivide.Draw("axissame")


    pad2.Update()
    pad1.cd()
    pad1.Update()
    pad1.Draw()

    if tLeg==None:
        if len(h_components)>=7 or (doSignalInj and len(h_components)>=5):
            tLeg = rt.TLegend(0.7,0.3,0.9,0.8)
        elif len(h_components)==6 or (doSignalInj and len(h_components)==4):
            tLeg = rt.TLegend(0.7,0.35,0.9,0.8)
        elif len(h_components)==5 or (doSignalInj and len(h_components)==3):
            tLeg = rt.TLegend(0.7,0.4,0.9,0.8)
        elif len(h_components)==4 or (doSignalInj and len(h_components)==2):
            tLeg = rt.TLegend(0.7,0.45,0.9,0.8)
        elif len(h_components)==3 or (doSignalInj and len(h_components)==1):
            tLeg = rt.TLegend(0.7,0.5,0.9,0.8)
        elif len(h_components)==2 or (doSignalInj and len(h_components)==0):
            tLeg = rt.TLegend(0.7,0.55,0.9,0.8)            
        else:
            tLeg = rt.TLegend(0.7,0.6,0.9,0.8)
        tLeg.SetFillColor(0)
        tLeg.SetTextFont(42)
        tLeg.SetLineColor(0)
        if isData:
            tLeg.AddEntry(h_data,"Data","lep")
        else:
            tLeg.AddEntry(h_data,"Sim. Data","lep")
        tLeg.AddEntry(hClone,"Fit Total","lf")
        for h_comp, color, label in zip(h_components, h_colors, h_labels):                
            tLeg.AddEntry(h_comp,label,"l")
        if doSignalInj:
            if options.model=="T1bbbb":
                tLeg.AddEntry(h_components[-1],'pp#rightarrow#tilde{g}#tilde{g}, #mu = %.1f'%options.r,'lf')
                tLeg.AddEntry(None,'#tilde{g}#rightarrowb#bar{b}#tilde{#chi}_{1}^{0}','')
            
    tLeg.Draw("same")

    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    if isData:
        l.DrawLatex(0.15,0.9,"CMS preliminary")
    else:
        l.DrawLatex(0.15,0.9,"CMS simulation")
    l.DrawLatex(0.75,0.9,"%s"%lumiLabel)
    l.SetTextFont(52)
    #l.SetTextSize(0.045)
    l.SetTextSize(0.055)
    l.DrawLatex(0.2,0.82,boxLabel)
    l.DrawLatex(0.3,0.77,plotLabel)
    if doSignalInj:        
        if options.model=="T1bbbb":
            l.SetTextFont(42)
            l.DrawLatex(0.3,0.77,"m_{#tilde{g}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(options.mGluino,options.mLSP))

    c.cd()
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')
    cWrite = c.Clone(os.path.splitext(printName)[0].split('/')[-1])
    rootFile.cd()
    c.Write(os.path.splitext(printName)[0].split('/')[-1])


def print1DSlice(c,rootFile,h_slices,h_data_slices,printName,xTitle,yTitle,lumiLabel="",boxLabel="",plotLabel="",isData=False,tLeg=None,h_colors=[],h_labels=[]):

    if densityCorr:
        h_slices_densitycorr = [densityCorrect(h_slice) for h_slice in h_slices]
        h_data_slices_densitycorr = [densityCorrect(h_data_slice) for h_data_slice in h_data_slices]
        h_slices = h_slices_densitycorr
        h_data_slices = h_data_slices_densitycorr
        
    pad1, pad2 = getPads(c)

    for h,color in zip(h_slices,h_colors):
        h.SetLineWidth(2)
        h.SetLineColor(color)

    for h_data,color in zip(h_data_slices,h_colors):
        h_data = setDataHist(h_data,xTitle,yTitle,densityCorr,color)
    #h_data_slices[0].SetMaximum(1e2*h_data_slices[0].GetMaximum())
            
    first = True
    for h_data,color in zip(h_data_slices,h_colors):
        if first:            
            h_data.Draw("pe")
            first = False
        else:
            h_data.Draw("pesame")
    
    for h,color in zip(h_slices,h_colors): h.Draw("histsame")
        
    for h_data, color in zip(h_data_slices,h_colors): h_data.Draw("pesame")
        
    pad1.Draw()
    c.Update()
    c.cd()
    pad2.Draw()
    pad2.cd()
    rt.gPad.SetLogy(0)

    hDataDivides=[]
    hDivides = []
    hCloneDivides = []
    for h, h_data, color in zip(h_slices, h_data_slices, h_colors):
        hClone = h.Clone(h.GetName()+"Clone")
        hDivide, hCloneDivide, hDataDivide  = getDivideHistos(h, hClone, h_data, xTitle, "Data/Fit")
        hCloneDivide.SetLineColor(rt.kWhite)
        
        hDataDivides.append(hDataDivide)
        hCloneDivides.append(hCloneDivide)
        hDivides.append(hDivide)

    hCloneDivides[0].Draw("axis")
    hDataDivides[0].Draw("pesame")
    for i, hDataDivide in enumerate(hDataDivides):
        if i>0: hDataDivide.Draw("pesame")
    #hCloneDivides[0].Draw("axissame")
    
    pad2.Update()
    pad1.cd()
    pad1.Update()
    pad1.Draw()

    if tLeg==None:
        if len(h_slices)>=8:
            tLeg = rt.TLegend(0.5,0.55,0.7,0.8)
            tLeg2 = rt.TLegend(0.7,0.55,0.9,0.8)
        elif len(h_slices)==7:
            tLeg = rt.TLegend(0.7,0.4,0.9,0.8)
        elif len(h_slices)==6:
            tLeg = rt.TLegend(0.7,0.4,0.9,0.8)
        elif len(h_slices)==5:
            tLeg = rt.TLegend(0.7,0.45,0.9,0.8)
        elif len(h_slices)==4:
            tLeg = rt.TLegend(0.7,0.55,0.9,0.8)
        elif len(h_slices)==3:
            tLeg = rt.TLegend(0.7,0.55,0.9,0.8)
        elif len(h_slices)==2:
            tLeg = rt.TLegend(0.7,0.55,0.9,0.8)            
        else:
            tLeg = rt.TLegend(0.7,0.6,0.9,0.8)
        tLeg.SetFillColor(0)
        tLeg.SetTextFont(42)
        tLeg.SetLineColor(0)
        if len(h_slices)>=8:
            tLeg2.SetFillColor(0)
            tLeg2.SetTextFont(42)
            tLeg2.SetLineColor(0)
            
        legendEntry = 0
        for h, color, label in zip(h_slices, h_colors, h_labels):
            legendEntry+=1     
            if legendEntry>4 and len(h_slices)>=8:
                tLeg2.AddEntry(h,label,"l")
            else:            
                tLeg.AddEntry(h,label,"l")
    tLeg.Draw("same")
    
    if len(h_slices)>=8:
        tLeg2.Draw("same")

    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    if isData:
        l.DrawLatex(0.15,0.9,"CMS preliminary")
    else:
        l.DrawLatex(0.15,0.9,"CMS simulation")
    l.DrawLatex(0.75,0.9,"%s"%lumiLabel)
    l.SetTextFont(52)
    l.SetTextSize(0.045)
    l.DrawLatex(0.2,0.82,boxLabel)
    l.DrawLatex(0.3,0.77,plotLabel)

    c.cd()
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')
    cWrite = c.Clone(os.path.splitext(printName)[0].split('/')[-1])
    rootFile.cd()
    c.Write(os.path.splitext(printName)[0].split('/')[-1])
    
def setFFColors(hNS, minZ=-5.1, maxZ=5.1):
    Red = array('d',  [0.00, 0.70, 0.90, 1.00, 1.00, 1.00, 1.00])
    Green = array('d',[0.00, 0.70, 0.90, 1.00, 0.90, 0.70, 0.00])
    Blue = array('d', [1.00, 1.00, 1.00, 1.00, 0.90, 0.70, 0.00])
    Length =array('d',[0.00, 0.20, 0.35, 0.50, 0.65, 0.8, 1.00]) # colors get darker faster at 4sigma
    rt.TColor.CreateGradientColorTable(7,Length,Red,Green,Blue,999)
    hNS.SetMaximum(maxZ)
    hNS.SetMinimum(minZ) # so the binning is 0 2 4
    hNS.SetContour(999)

    
def setRainbowColors(hNS, minZ=0, maxZ=100):
    rt.gStyle.SetPalette(1)
    Red = array('d',  [0.00, 0.00, 0.87, 1.00, 0.51])
    Green = array('d',[0.00, 0.81, 1.00, 0.20, 0.00])
    Blue = array('d', [0.51, 1.00, 0.12, 0.00, 0.00])
    Length =array('d',[0.00, 0.34, 0.61, 0.84, 1.00]) 
    rt.TColor.CreateGradientColorTable(5,Length,Red,Green,Blue,999)
    hNS.SetMaximum(maxZ)
    hNS.SetMinimum(minZ) # so the binning is 0 2 4
    hNS.SetContour(999)

    
def print2DCanvas(c,rootFile,h,printName):
    c.SetLogx(1)
    c.SetLogy(1)
    h.Draw("textcolz")
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')
    cWrite = c.Clone(os.path.splitext(printName)[0].split('/')[-1])
    rootFile.cd()
    c.Write(os.path.splitext(printName)[0].split('/')[-1])

def set2DHisto(h2D,xTitle,yTitle,zTitle):    
    h2D.GetXaxis().SetMoreLogLabels()
    h2D.GetXaxis().SetNoExponent()
    h2D.GetYaxis().SetMoreLogLabels()
    h2D.GetYaxis().SetNoExponent()
    h2D.SetMarkerSize(1.5)
    h2D.GetXaxis().SetTitle(xTitle)
    h2D.GetYaxis().SetTitle(yTitle)
    h2D.GetZaxis().SetTitle(zTitle)    
    h2D.GetXaxis().SetTitleSize(0.056)
    h2D.GetXaxis().SetLabelSize(0.056)
    h2D.GetYaxis().SetTitleSize(0.056)
    h2D.GetYaxis().SetLabelSize(0.056)
    h2D.GetZaxis().SetLabelSize(0.056)
    h2D.GetZaxis().SetTitleSize(0.056)
    h2D.GetZaxis().SetTitleOffset(1)
    #h2D.GetXaxis().SetTitleOffset(0.8)
    #h2D.GetYaxis().SetTitleOffset(0.7)
    #h2D.GetXaxis().SetTicks("+-")
    return h2D

def getGrayLines(x,y,sidebandFit=None):        
    # the gray lines
    xLines = []
    yLines = []

    lastX = len(x)-1
    lastY = len(y)-1

    for i in range(1,lastY):
        xLines.append(rt.TLine(x[0], y[i], x[lastX], y[i]))
        xLines[i-1].SetLineStyle(2);
        xLines[i-1].SetLineColor(rt.kGray);
        
    for i in range(1,lastX):
        yLines.append(rt.TLine(x[i], y[0], x[i], y[lastY]))
        yLines[i-1].SetLineStyle(2)
        yLines[i-1].SetLineColor(rt.kGray)

        
    if sidebandFit:
        mrSide = sidebandFit[0]
        rsqSide = sidebandFit[1]
        yLines.append(rt.TLine(mrSide, rsqSide, mrSide, y[-1]))
        yLines[-1].SetLineStyle(2)
        yLines[-1].SetLineWidth(2)
        yLines[-1].SetLineColor(rt.kGreen)
        xLines.append(rt.TLine(mrSide, rsqSide, x[-1], rsqSide))
        xLines[-1].SetLineStyle(2)
        xLines[-1].SetLineWidth(2)
        xLines[-1].SetLineColor(rt.kGreen)
        
    return xLines,yLines

def dressFrenchFlag(hNS):    
    fGrayGraphs = []
    tlatexList = []
    col1 = rt.gROOT.GetColor(rt.kGray+1)
    col1.SetAlpha(0.3)
    for iBinX in range(1,hNS.GetNbinsX()+1):
        for iBinY in range(1,hNS.GetNbinsY()+1):
            if hNS.GetBinContent(iBinX,iBinY)!= -999: continue
            xBinLow = hNS.GetXaxis().GetBinLowEdge(iBinX)
            xBinHigh = xBinLow+hNS.GetXaxis().GetBinWidth(iBinX)
            yBinLow = hNS.GetYaxis().GetBinLowEdge(iBinY)
            yBinHigh = yBinLow+hNS.GetYaxis().GetBinWidth(iBinY)
            fGray = rt.TGraph(5)
            fGray.SetPoint(0,xBinLow,yBinLow)
            fGray.SetPoint(1,xBinLow,yBinHigh)
            fGray.SetPoint(2,xBinHigh,yBinHigh)
            fGray.SetPoint(3,xBinHigh,yBinLow)
            fGray.SetPoint(4,xBinLow,yBinLow)
            fGray.SetFillColor(rt.kGray+1)
            fGrayGraphs.append(fGray)
    for iBinX in range(1,hNS.GetNbinsX()+1):
        for iBinY in range(1,hNS.GetNbinsY()+1):
            binCont = hNS.GetBinContent(iBinX,iBinY)
            if binCont == -999: continue
            if abs(binCont)==0: continue
            yBin = hNS.GetYaxis().GetBinLowEdge(iBinY) + .3*hNS.GetYaxis().GetBinWidth(iBinY) # bottom of TLatex 30% across the binwidth in X
            if iBinX==1 and binCont>=0.:
                xBin  = hNS.GetXaxis().GetBinLowEdge(iBinX) + .1*hNS.GetXaxis().GetBinWidth(iBinX)
            elif iBinX==1 and binCont<0.:
                xBin  = hNS.GetXaxis().GetBinLowEdge(iBinX) + .05*hNS.GetXaxis().GetBinWidth(iBinX) # left side of TLatex 10% across the binwidth in X
            elif binCont>=0.:
                xBin  = hNS.GetXaxis().GetBinLowEdge(iBinX) + .25*hNS.GetXaxis().GetBinWidth(iBinX)
            elif binCont<0.:
                xBin  = hNS.GetXaxis().GetBinLowEdge(iBinX) + .1*hNS.GetXaxis().GetBinWidth(iBinX) # left side of TLatex 10% across the binwidth in X

            tlatex = rt.TLatex(xBin,yBin,"%2.1f"%binCont)
            tlatex.SetTextSize(0.04)
            tlatex.SetTextFont(42)
            tlatexList.append(tlatex)
    return fGrayGraphs, tlatexList

def set2DCanvas(c):
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    c.SetRightMargin(0.17)
    c.SetLogx(1)
    c.SetLogy(1)
    c.SetLogz(0)
    return c

def print2DScatter(c,rootFile,h,printName,xTitle,yTitle,zTitle,lumiLabel,boxLabel,plotLabel,x,y,zMin,zMax,isData=False,sidebandFit=None,doSignalInj=False,options=None,drawOpt="colz"):

    c = set2DCanvas(c)
    c.SetLogz(1)
    
    h = set2DHisto(h,xTitle,yTitle,zTitle)
    setRainbowColors(h,zMin,zMax)
    
    h.Draw(drawOpt)
    xLines, yLines = getGrayLines(x,y,sidebandFit)
    
    [xLine.Draw("l") for xLine in xLines]
    [yLine.Draw("l") for yLine in yLines]

    
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    if isData:
        l.DrawLatex(0.15,0.91,"CMS preliminary")
    else:
        l.DrawLatex(0.15,0.91,"CMS simulation")
    l.DrawLatex(0.65,0.91,"%s"%lumiLabel)
    l.SetTextFont(52)
    l.SetTextSize(0.04)
    l.DrawLatex(0.2,0.85,boxLabel)
    l.DrawLatex(0.3,0.8,plotLabel)
    if doSignalInj:        
        if options.model=="T1bbbb":
            l.SetTextFont(42)
            l.DrawLatex(0.3,0.8,'pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrowb#bar{b}#tilde{#chi}_{1}^{0}')
            l.DrawLatex(0.3,0.75,"m_{#tilde{g}} = %i GeV, m_{#tilde{#chi}} = %i GeV, #mu = %.1f"%(options.mGluino,options.mLSP,options.r))
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')    
    cWrite = c.Clone(os.path.splitext(printName)[0].split('/')[-1])
    rootFile.cd()
    c.Write(os.path.splitext(printName)[0].split('/')[-1])
    c.SetLogy(0)
    c.SetLogx(0)
    #c.Print(printName.replace('log','lin'))
    #c.Print(os.path.splitext(printName.replace('log','lin'))[0]+'.C')
    #cWrite = c.Clone(os.path.splitext(printName.replace('log','lin'))[0].split('/')[-1])
    #rootFile.cd()
    #c.Write(os.path.splitext(printName.replace('log','lin'))[0].split('/')[-1])
    
def print2DResiduals(c,rootFile,h_resi,printName,xTitle,yTitle,zTitle,lumiLabel,boxLabel,plotLabel,x,y,isData=False,sidebandFit=None,doSignalInj=False,options=None,drawOpt="colz"):
    
    c = set2DCanvas(c)
    absMax = max(abs(h_resi.GetMinimum()),abs(h_resi.GetMaximum()))    
    h_resi = set2DHisto(h_resi,xTitle,yTitle,zTitle)
    
    if "nsigma" in h_resi.GetName():
        setFFColors(h_resi,-5.1,5.1)
    else:
        setFFColors(h_resi,-1.5*absMax,1.5*absMax)
    h_resi.Draw(drawOpt)
    xLines, yLines = getGrayLines(x,y,sidebandFit)

    fGrayGraphs, tlatexList = dressFrenchFlag(h_resi)
    
    [xLine.Draw("l") for xLine in xLines]
    [yLine.Draw("l") for yLine in yLines]
    #[fGray.Draw("f") for fGray in fGrayGraphs]    
    [tlatex.Draw() for tlatex in tlatexList]
    
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()    
    if isData:
        l.DrawLatex(0.15,0.91,"CMS preliminary")
    else:
        l.DrawLatex(0.15,0.91,"CMS simulation")
    l.DrawLatex(0.65,0.91,"%s"%lumiLabel)
    l.SetTextFont(52)
    l.SetTextSize(0.04)
    l.DrawLatex(0.2,0.85,boxLabel)
    l.DrawLatex(0.3,0.8,plotLabel)
    if doSignalInj:        
        if options.model=="T1bbbb":
            l.SetTextFont(42)
            l.DrawLatex(0.3,0.8,'pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrowb#bar{b}#tilde{#chi}_{1}^{0}')
            l.DrawLatex(0.3,0.75,"m_{#tilde{g}} = %i GeV, m_{#tilde{#chi}} = %i GeV, #mu = %.1f"%(options.mGluino,options.mLSP,options.r))
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')
    cWrite = c.Clone(os.path.splitext(printName)[0].split('/')[-1])
    rootFile.cd()
    c.Write(os.path.splitext(printName)[0].split('/')[-1])
    c.SetLogy(0)
    c.SetLogx(0)
    #c.Print(printName.replace('log','lin'))
    #c.Print(os.path.splitext(printName.replace('log','lin'))[0]+'.C')
    #cWrite = c.Clone(os.path.splitext(printName.replace('log','lin'))[0].split('/')[-1])
    #rootFile.cd()
    #c.Write(os.path.splitext(printName.replace('log','lin'))[0].split('/')[-1])
        

def get3DHistoFrom1D(h1D,x,y,z,name):   
    h3D = rt.TH3D(name,name,len(x)-1,x,len(y)-1,y,len(z)-1,z)

    iBinX=-1
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            for k in range(1,len(z)):
                iBinX += 1
                h3D.SetBinContent(i,j,k,h1D.GetBinContent(iBinX+1))
                h3D.SetBinError(i,j,k,h1D.GetBinError(iBinX+1))
    return h3D


def get1DHistoFrom2D(h2D,x,y,name):
    nBins = (len(x)-1)*(len(y)-1)
    h1D = rt.TH1D(name,name,nBins,0,nBins)

    iBinX=-1
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            iBinX += 1
            h1D.SetBinContent(iBinX+1,h2D.GetBinContent(i,j))
            h1D.SetBinError(iBinX+1,h2D.GetBinError(i,j))       
    return h1D


def get1DHistoFrom3D(h3D,x,y,z,name):
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    h1D = rt.TH1D(name,name,nBins,0,nBins)

    iBinX=-1
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            for k in range(1,len(z)):
                iBinX += 1
                h1D.SetBinContent(iBinX+1,h3D.GetBinContent(i,j,k))
                h1D.SetBinError(iBinX+1,h3D.GetBinError(i,j,k))
    return h1D

def Gamma(a, x):
    return rt.TMath.Gamma(a) * rt.Math.inc_gamma_c(a,x)

def Gfun(x, y, X0, Y0, B, N):
    return Gamma(N,B*N*rt.TMath.Power((x-X0)*(y-Y0),1/N))

def integrand(Y,X0,Y0,Y1,Y2,B,N,xmin,xmax):
    integral = ( (xmin-X0)*rt.TMath.Exp(B*N*rt.TMath.Power(xmax-X0,1/N)*rt.TMath.Power(Y-Y0,1/N)) - (xmax-X0)*rt.TMath.Exp(B*N*rt.TMath.Power(xmin-X0,1/N)*rt.TMath.Power(Y-Y0,1/N)) )*rt.TMath.Exp(-B*N*(rt.TMath.Power(xmin-X0,1/N)+rt.TMath.Power(xmax-X0,1/N))*rt.TMath.Power(Y-Y0,1/N))
    suppress = 0.5*rt.TMath.Erfc((Y-Y1)/Y2)
    
    return integral*suppress

def getBinEvents(i, j, k, x, y, z, workspace,box):
    errorFlag = False
    
    bkg = "TTj%ib"%z[k-1]

    pdfs = workspace.allPdfs()

    B = -99999999
    N = -99999999
    X0 = -99999999
    Y0 = -99999999
    Y1 = -99999999
    Y2 = -99999999
    NTOT = -99999999
    for pdf in rootTools.RootIterator.RootIterator(pdfs):
        pdfParams = pdf.getParameters(workspace.data("obs_data"))        
        for pdfParam in rootTools.RootIterator.RootIterator(pdfParams):
            if bkg in pdf.GetName():
                if pdfParam.GetName().find("b_")==0:
                    B = pdfParam.getVal()                    
                elif pdfParam.GetName().find("n_")==0:
                    N = pdfParam.getVal()                    
                elif pdfParam.GetName().find("MR0_")==0:
                    X0 = pdfParam.getVal()
                elif pdfParam.GetName().find("R0_")==0:
                    Y0 = pdfParam.getVal()
                elif pdfParam.GetName().find("R1_")==0:
                    Y1 = pdfParam.getVal()
                elif pdfParam.GetName().find("R2_")==0:
                    Y2 = pdfParam.getVal()
            if "Ntot_%s"%bkg in pdfParam.GetName():
                NTOT = pdfParam.getVal()
    
    xmin  = x[0]
    xmax  = x[-1]
    ymin  = y[0]
    ymax  = y[-1]
    if Y1>-99999999 and Y2>-99999999:
        total_numerical_integral = quad(integrand, ymin, ymax, args=(X0,Y0,Y1,Y2,B,N,xmin,xmax))
        total_integral = total_numerical_integral[0]
    else:
        #total_integral = (N/rt.TMath.Power(B*N,N))*(Gfun(xmin,ymin,X0,Y0,B,N)-Gfun(xmin,ymax,X0,Y0,B,N)-Gfun(xmax,ymin,X0,Y0,B,N)+Gfun(xmax,ymax,X0,Y0,B,N))
        total_integral = Gfun(xmin,ymin,X0,Y0,B,N)-Gfun(xmin,ymax,X0,Y0,B,N)-Gfun(xmax,ymin,X0,Y0,B,N)+Gfun(xmax,ymax,X0,Y0,B,N)
        
    xmin  = x[i-1]
    xmax  = x[i]
    ymin  = y[j-1]
    ymax  = y[j]
    if Y1>-99999999 and Y2>-99999999:
        numerical_integral = quad(integrand, ymin, ymax, args=(X0,Y0,Y1,Y2,B,N,xmin,xmax))
        integral = numerical_integral[0]
    else:
        #integral = (N/rt.TMath.Power(B*N,N))*(Gfun(xmin,ymin,X0,Y0,B,N)-Gfun(xmin,ymax,X0,Y0,B,N)-Gfun(xmax,ymin,X0,Y0,B,N)+Gfun(xmax,ymax,X0,Y0,B,N))
        integral = Gfun(xmin,ymin,X0,Y0,B,N)-Gfun(xmin,ymax,X0,Y0,B,N)-Gfun(xmax,ymin,X0,Y0,B,N)+Gfun(xmax,ymax,X0,Y0,B,N)

    if integral > 0. and total_integral > 0.:        
        bin_events =  NTOT*integral/total_integral

    else: 
        errorFlag = True
        #print "ERROR: bin razor pdf integral =", integral
        #print "ERROR: total razor pdf integral =", total_integral
        return 0., errorFlag
    
    return bin_events, errorFlag

def getErrors1D(h, h_data,  myTree, options, sumType,minX, maxX, minY, maxY, minZ, maxZ, x, y, z):
    binSumDict = getBinSumDicts(sumType, minX, maxX, minY, maxY, minZ, maxZ, x, y, z)
    
    d = rt.TCanvas('d','d',500,400)
    for i, sumName in binSumDict.iteritems():
        nObs = h_data.GetBinContent(i)
        #print "---- geterrors1d"
        bestFit, rms, pvalue, nsigma, d = getBestFitRms(myTree,sumName,nObs,d,options,"%s_error_%i.pdf"%(h.GetName(),i))
        h.SetBinError(i,rms)
        #print "i,rms=", i, rms
    return h

def getErrors2D(h, h_data,  myTree, options, sumType,minX, maxX, minY, maxY, minZ, maxZ, x, y, z):
    binSumDict = getBinSumDicts(sumType, minX, maxX, minY, maxY, minZ, maxZ, x, y, z)
    
    d = rt.TCanvas('d','d',500,400)
    for (i,j), sumName in binSumDict.iteritems():
        nObs = h_data.GetBinContent(i, j)
        #print "---- geterrors2d"
        bestFit, rms, pvalue, nsigma, d = getBestFitRms(myTree,sumName,nObs,d,options,"%s_error_%i_%i.pdf"%(h.GetName(),i,j))
        #print "i, j, rms=", i, j, rms
        h.SetBinError(i,j,rms)
    return h

def getErrors3D(h, h_data,  myTree, options, sumType,minX, maxX, minY, maxY, minZ, maxZ, x, y, z):
    binSumDict = getBinSumDicts(sumType, minX, maxX, minY, maxY, minZ, maxZ, x, y, z)
    
    d = rt.TCanvas('d','d',500,400)
    for (i,j, k), sumName in binSumDict.iteritems():
        nObs = h_data.GetBinContent(i, j, k)
        #print "---- geterrors3d"
        bestFit, rms, pvalue, nsigma, d = getBestFitRms(myTree,sumName,nObs,d,options,"%s_error_%i_%i_%i.pdf"%(h.GetName(),i,j,k))
        #print "i, j, k, rms=", i, j, k, rms
        h.SetBinError(i,j,k,rms)
    return h

def getNsigma2D(h, h_data,  myTree, options, sumType,minX, maxX, minY, maxY, minZ, maxZ, x, y, z):
    binSumDict = getBinSumDicts(sumType, minX, maxX, minY, maxY, minZ, maxZ, x, y, z)
    
    d = rt.TCanvas('d','d',500,400)
    for (i,j), sumName in binSumDict.iteritems():
        nObs = h_data.GetBinContent(i,j)
        #print "---- getnsigma2d"
        bestFit, rms, pvalue, nsigma, d = getBestFitRms(myTree,sumName,nObs,d,options,"%s_error_%i_%i.pdf"%(h.GetName(),i,j))
        h.SetBinContent(i,j,nsigma)
        #print "i, j, rms=", i, j, rms
        #if abs(nsigma)==0 and bestFit<1 and nObs==0:
        #    h.SetBinContent(i,j,-999)
    return h


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store results")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default="BinnedFitResults_MultiJet.root",type="string",
                  help="input fit file")
    parser.add_option('-t','--input-toy-file',dest="inputToyFile", default=None,type="string",
                  help="input toy file with stat+sys uncertainties")
    parser.add_option('-s','--input-sys-file',dest="inputSysFile", default=None,type="string",
                  help="input toy file with sys-only uncertainties")
    parser.add_option('--print-errors',dest="printErrors", default=False,action='store_true',
                  help="print plots of individual error calculation")
    parser.add_option('--data',dest="isData", default=False,action='store_true',
                  help="changes plots for data")
    parser.add_option('--fit-region',dest="fitRegion",default="Full",type="string",
                  help="Fit region")
    parser.add_option('--plot-region',dest="plotRegion",default="Full",type="string",
                  help="Plot region")
    parser.add_option('-w','--weight',dest="useWeight",default=False,action='store_true',
                  help="use weight")
    parser.add_option('--no-stat', dest='noStat', default=False, action='store_true',
                  help='make 1d plots with toys thrown with systematic uncertainties only')
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model")
    parser.add_option('--mGluino',dest="mGluino", default=1500,type="float",
                  help="mgluino")
    parser.add_option('--mLSP',dest="mLSP", default=100,type="float",
                  help="mass of LSP")
    parser.add_option('-r','--signal-strength',dest="r", default=-1,type="float",
                  help="signal strength")
    
    (options,args) = parser.parse_args()
     
    box = options.box
    lumi = options.lumi
    cfg = Config.Config(options.config)
    fitRegion = options.fitRegion
    plotRegion = options.plotRegion
    doSignalInj = (options.r > -1)
    
    print "plotting in range:", plotRegion
    print "fit in range:", fitRegion 

    inputFitFile = rt.TFile.Open(options.inputFitFile,"read")

            
    toyTree = None
    if options.inputToyFile is not None:
        toyFiles = options.inputToyFile.split(',')
        toyTree = rt.TChain("myTree")
        for toyFile in toyFiles:
            toyTree.Add(toyFile)
        toyTree.SetName("toyTree")
        
    computeErrors = (toyTree is not None)

    w = inputFitFile.Get("w"+box)
        
    th1x = w.var('th1x')
    dataHist = w.data("data_obs")

    setStyle()
    
    extRazorPdf = w.pdf('extRazorPdf')
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)

    sysTree = toyTree
    if (options.inputSysFile is not None) and options.noStat:        
        sysFiles = options.inputSysFile.split(',')
        sysTree = rt.TChain("myTree")
        for sysFile in sysFiles:
            sysTree.Add(sysFile)
        sysTree.SetName("sysTree")

    xFine = array('d', [x[0]+i*(x[-1]-x[0])/100. for i in range(0,101)]) # MR binning fine
    yFine = array('d', [y[0]+i*(y[-1]-y[0])/100. for i in range(0,101)]) # Rsq binning fine
    zFine = array('d', cfg.getBinning(box)[2]) # nBtag binning fine
    nBinsFine = (len(xFine)-1)*(len(yFine)-1)*(len(zFine)-1)    
    
    th1x.setBins(nBins)

    asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('central'),rt.RooFit.Asimov())
    
    
    rt.TH1D.SetDefaultSumw2()
    rt.TH2D.SetDefaultSumw2()
    rt.TH3D.SetDefaultSumw2()
    
    # start writing output
    c = rt.TCanvas('c','c',500,400)    
    rootFile = rt.TFile.Open(options.outDir + '/' + 'Plots_%s'%box + '.root','recreate')
    tdirectory = rootFile.GetDirectory(options.outDir)
    if tdirectory==None:
        rootFile.mkdir(options.outDir)
        tdirectory = rootFile.GetDirectory(options.outDir)

        
    plotband = convertSideband(plotRegion,w,x,y,z)
    
    opt = [rt.RooFit.CutRange(myRange) for myRange in plotband.split(',')]
    
    asimov_reduce = asimov.reduce(opt[0])
    dataHist_reduce = dataHist.reduce(opt[0])
    for iOpt in range(1,len(opt)):
        asimov_reduce.add(asimov.reduce(opt[iOpt]))
        dataHist_reduce.add(dataHist.reduce(opt[iOpt]))

        
    h_th1x = asimov_reduce.createHistogram('h_th1x',th1x)
    h_data_th1x = dataHist_reduce.createHistogram('h_data_th1x',th1x)
    
    h_data_nBtagRsqMR_fine = rt.TH3D("h_data_nBtagRsqMR_fine","h_data_nBtagRsqMR_fine",len(xFine)-1,xFine,len(yFine)-1,yFine,len(zFine)-1,zFine)
    h_nBtagRsqMR_fine = rt.TH3D("h_nBtagRsqMR_fine","h_nBtagRsqMR_fine",len(xFine)-1,xFine,len(yFine)-1,yFine,len(zFine)-1,zFine)
    
    opt = [rt.RooFit.CutRange(myRange) for myRange in plotRegion.split(',')]
    data_reduce = w.data("RMRTree").reduce(opt[0])
    for iOpt in range(1,len(opt)):
        data_reduce.append(w.data("RMRTree").reduce(opt[iOpt]))
    data_reduce.fillHistogram(h_data_nBtagRsqMR_fine,rt.RooArgList(w.var("MR"),w.var("Rsq"),w.var("nBtag")))

    if doSignalInj:        
        extSignalPdf = w.pdf('extSignalPdf')
        signal = extSignalPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('signal'),rt.RooFit.Asimov())
        signal_reduce = signal.reduce(*opt)
        h_sig_th1x = signal_reduce.createHistogram('h_sig_th1x',th1x)
        h_sig_nBtagRsqMR = get3DHistoFrom1D(h_sig_th1x,x,y,z,"h_sig_nBtagRsqMR")
        h_sig_RsqMR = h_sig_nBtagRsqMR.Project3D("yxe")
        h_sig_MR = h_sig_nBtagRsqMR.Project3D("xe")
        h_sig_Rsq = h_sig_nBtagRsqMR.Project3D("ye")
        
    for i in range(1,len(xFine)):
        for j in range(1,len(yFine)):
            for k in range(1,len(zFine)):
                w.var('MR').setVal((xFine[i]+xFine[i-1])/2.)
                w.var('Rsq').setVal((yFine[j]+yFine[j-1])/2.)
                w.var('nBtag').setVal((zFine[k]+zFine[k-1])/2.)
                inSideband = 0
                for myRange in plotRegion.split(','):
                    inSideband += ( w.var('MR').inRange(myRange) * w.var('Rsq').inRange(myRange) * w.var('nBtag').inRange(myRange) )
                if not inSideband: continue  
                value, errorFlag = getBinEvents(i,j,k,xFine,yFine,zFine,w,box)
                if not errorFlag:
                    h_nBtagRsqMR_fine.SetBinContent(i,j,k,value)

                    
    h_data_RsqMR_fine = h_data_nBtagRsqMR_fine.Project3D("yxe")
    h_RsqMR_fine = h_nBtagRsqMR_fine.Project3D("yxe")
    
    h_data_nBtagRsqMR = get3DHistoFrom1D(h_data_th1x,x,y,z,"h_data_nBtagRsqMR")
    h_nBtagRsqMR = get3DHistoFrom1D(h_th1x,x,y,z,"h_nBtagRsqMR")    
        
    h_data_RsqMR = h_data_nBtagRsqMR.Project3D("yxe")
    h_data_RsqMR.SetName("h_data_RsqMR")
    h_data_MR = h_data_nBtagRsqMR.Project3D("xe")
    h_data_MR.SetName("h_data_MR")
    h_data_Rsq = h_data_nBtagRsqMR.Project3D("ye")
    h_data_Rsq.SetName("h_data_Rsq")
    h_RsqMR = h_nBtagRsqMR.Project3D("yxe")
    h_RsqMR.SetName("h_RsqMR")
    h_MR = h_nBtagRsqMR.Project3D("xe")
    h_MR.SetName("h_MR")
    h_Rsq = h_nBtagRsqMR.Project3D("ye")
    h_Rsq.SetName("h_Rsq")

    if computeErrors:
        h_MR = getErrors1D(h_MR,h_data_MR,sysTree,options,"x",0,len(x)-1,0,len(y)-1,0,len(z)-1,x,y,z)
        h_Rsq = getErrors1D(h_Rsq,h_data_Rsq,sysTree,options,"y",0,len(x)-1,0,len(y)-1,0,len(z)-1,x,y,z)
        h_RsqMR = getErrors2D(h_RsqMR,h_data_RsqMR,sysTree,options,"yx",0,len(x)-1,0,len(y)-1,0,len(z)-1,x,y,z)
        h_nBtagRsqMR = getErrors3D(h_nBtagRsqMR,h_data_nBtagRsqMR,sysTree,options,"zyx",0,len(x)-1,0,len(y)-1,0,len(z)-1,x,y,z)
        
        # convert after getting the 3D errors if provided
        h_th1x = get1DHistoFrom3D(h_nBtagRsqMR,x,y,z,'h_th1x_wErrors')
        
    h_MR_slices = []
    h_MR_integrals = []
    h_data_MR_slices = []
    h_data_MR_integrals = []
    h_MR_slice_labels = []
    h_MR_integral_labels = []

    h_MR_slice_components = []
    h_data_MR_slice_components = []
    h_MR_slice_component_labels = []
    h_MR_integral_components = []
    h_data_MR_integral_components = []
    h_MR_integral_component_labels = []
    for j in range(1,len(y)):
        h_MR_slices.append(h_nBtagRsqMR.ProjectionX(("h_MR_%.2fRsq%.2f"%(y[j-1],y[j])).replace('.','p'),j,j,0,-1,""))
        h_MR_integrals.append(h_nBtagRsqMR.ProjectionX(("h_MR_Rsq%.2f"%(y[j-1])).replace('.','p'),j,len(y)-1,0,-1,""))
        h_data_MR_slices.append(h_data_nBtagRsqMR.ProjectionX(("h_MR_data_%.2fRsq%.2f"%(y[j-1],y[j])).replace('.','p'),j,j,0,-1,""))
        h_data_MR_integrals.append(h_data_nBtagRsqMR.ProjectionX(("h_MR_data_Rsq%.2f"%(y[j-1])).replace('.','p'),j,len(y)-1,0,-1,""))
        h_MR_slice_labels.append("%.2f #leq R^{2} < %.2f"%(y[j-1],y[j]))
        h_MR_integral_labels.append("R^{2} #geq %.2f"%(y[j-1]))
        if computeErrors:
            h_MR_slices[-1] = getErrors1D(h_MR_slices[-1],h_data_MR_slices[-1],sysTree,options,"x",0,len(x)-1,j,j,0,len(z)-1,x,y,z)
            h_MR_integrals[-1] = getErrors1D(h_MR_integrals[-1],h_data_MR_integrals[-1],sysTree,options,"x",0,len(x)-1,j,len(y)-1,0,len(z)-1,x,y,z)            
        if len(z)>2:
            h_MR_comp = []
            h_data_MR_comp = []
            h_label_comp = []
            h_MR_int_comp = []
            h_data_MR_int_comp = []
            h_label_int_comp = []
            for k in range(1,len(z)):
                h_MR_comp.append(h_nBtagRsqMR.ProjectionX(("h_MR_%ibtag_%.2fRsq%.2f"%(z[k-1],y[j-1],y[j])).replace('.','p'),j,j,k,k,""))
                h_data_MR_comp.append(h_data_nBtagRsqMR.ProjectionX(("h_MR_data_%ibtag_%.2fRsq%.2f"%(z[k-1],y[j-1],y[j])).replace('.','p'),j,j,k,k,"")) 
                h_MR_int_comp.append(h_nBtagRsqMR.ProjectionX(("h_MR_%ibtag_Rsq%.2f"%(z[k-1],y[j-1])).replace('.','p'),j,len(y)-1,k,k,""))               
                h_data_MR_int_comp.append(h_data_nBtagRsqMR.ProjectionX(("h_MR_data_%ibtag_Rsq%.2f"%(z[k-1],y[j-1])).replace('.','p'),j,len(y)-1,k,k,"")) 
                if computeErrors:
                    h_MR_comp[-1] = getErrors1D(h_MR_comp[-1],h_data_MR_comp[-1],sysTree,options,"x",0,len(x)-1,j,j,k,k,x,y,z)
                    h_MR_int_comp[-1] = getErrors1D(h_MR_int_comp[-1],h_data_MR_int_comp[-1],sysTree,options,"x",0,len(x)-1,j,len(y)-1,k,k,x,y,z)
                if z[k-1]==3 and z[-1]==4:
                    h_label_comp.append("#geq %i b-tag, %.2f #leq R^{2} < %.2f " % (z[k-1],y[j-1],y[j]) )
                    h_label_int_comp.append("#geq %i b-tag, R^{2} #geq %.2f " % (z[k-1],y[j-1]) )
                if z[k-1]==1 and z[-1]==4 and box in ['MuEle','MuMu','EleEle']:                
                    h_label_comp.append("#geq %i b-tag, %.2f #leq R^{2} < %.2f " % (z[k-1],y[j-1],y[j]) )
                    h_label_int_comp.append("#geq %i b-tag, R^{2} #geq %.2f " % (z[k-1],y[j-1]) )
                else:            
                    h_label_comp.append("%i b-tag, %.2f #leq R^{2} < %.2f " % (z[k-1],y[j-1],y[j]) )
                    h_label_int_comp.append("%i b-tag, R^{2} #geq %.2f " % (z[k-1],y[j-1]) )
            h_MR_slice_components.append(h_MR_comp)
            h_MR_integral_components.append(h_MR_int_comp)
            h_data_MR_slice_components.append(h_data_MR_comp)
            h_data_MR_integral_components.append(h_data_MR_int_comp)
            h_MR_slice_component_labels.append(h_label_comp)
            h_MR_integral_component_labels.append(h_label_int_comp)
        
    h_Rsq_slices = []
    h_Rsq_integrals = []
    h_data_Rsq_slices = []
    h_data_Rsq_integrals = []
    h_Rsq_slice_labels = []
    h_Rsq_integral_labels = []
    
    h_Rsq_slice_components = []
    h_data_Rsq_slice_components = []
    h_Rsq_slice_component_labels = []
    h_Rsq_integral_components = []
    h_data_Rsq_integral_components = []
    h_Rsq_integral_component_labels = []
    for i in range(1,len(x)):
        h_Rsq_slices.append(h_nBtagRsqMR.ProjectionY("h_Rsq_%iMR%i"%(x[i-1],x[i]),i,i,0,-1,""))
        h_Rsq_integrals.append(h_nBtagRsqMR.ProjectionY("h_Rsq_MR%i"%(x[i-1]),i,len(x)-1,0,-1,""))
        h_data_Rsq_slices.append(h_data_nBtagRsqMR.ProjectionY("h_Rsq_data_%iMR%i"%(x[i-1],x[i]),i,i,0,-1,""))
        h_data_Rsq_integrals.append(h_data_nBtagRsqMR.ProjectionY("h_Rsq_data_MR%i"%(x[i-1]),i,len(x)-1,0,-1,""))
        h_Rsq_slice_labels.append("%i #leq M_{R} < %i"%(x[i-1],x[i]))
        h_Rsq_integral_labels.append("M_{R} #geq %i"%(x[i-1]))
        
        if computeErrors:
            h_Rsq_slices[-1] = getErrors1D(h_Rsq_slices[-1],h_data_Rsq_slices[-1],sysTree,options,"y",i,i,0,len(y)-1,0,len(z)-1,x,y,z)
            h_Rsq_integrals[-1] = getErrors1D(h_Rsq_integrals[-1],h_data_Rsq_integrals[-1],sysTree,options,"y",i,i,0,len(y)-1,0,len(z)-1,x,y,z)
        if len(z)>2:
            h_Rsq_comp = []
            h_data_Rsq_comp = []
            h_label_comp = []
            h_Rsq_int_comp = []
            h_data_Rsq_int_comp = []
            h_label_int_comp = []
            for k in range(1,len(z)):
                h_Rsq_comp.append(h_nBtagRsqMR.ProjectionY("h_Rsq_%ibtag_%iMR%i"%(z[k-1],x[i-1],x[i]),i,i,k,k,""))
                h_Rsq_int_comp.append(h_nBtagRsqMR.ProjectionY("h_Rsq_%ibtag_MR%i"%(z[k-1],x[i-1]),i,len(x)-1,k,k,""))
                h_data_Rsq_comp.append(h_data_nBtagRsqMR.ProjectionY("h_Rsq_data_%ibtag_%iMR%i"%(z[k-1],x[i-1],x[i]),i,i,k,k,""))     
                h_data_Rsq_int_comp.append(h_data_nBtagRsqMR.ProjectionY("h_Rsq_data_%ibtag_MR%i"%(z[k-1],x[i-1]),i,len(x)-1,k,k,""))                
                if computeErrors:
                    h_Rsq_comp[-1] = getErrors1D(h_Rsq_comp[-1],h_data_Rsq_comp[-1],sysTree,options,"y",i,i,0,len(y)-1,k,k,x,y,z)    
                    h_Rsq_int_comp[-1] = getErrors1D(h_Rsq_int_comp[-1],h_data_Rsq_int_comp[-1],sysTree,options,"y",i,len(x)-1,0,len(y)-1,k,k,x,y,z)           
                if z[k-1]==3 and z[-1]==4:
                    h_label_comp.append("#geq %i b-tag, %i #leq M_{R} < %i " % (z[k-1],x[i-1],x[i]) )
                    h_label_int_comp.append("#geq %i b-tag, M_{R} #geq %i " % (z[k-1],x[i-1]) )
                if z[k-1]==1 and z[-1]==4 and box in ['MuEle','MuMu','EleEle']:                
                    h_label_comp.append("#geq %i b-tag, %i #leq M_{R} < %i " % (z[k-1],x[i-1],x[i]) )
                    h_label_int_comp.append("#geq %i b-tag, M_{R} #geq %i " % (z[k-1],x[i-1]) )
                else:            
                    h_label_comp.append("%i b-tag, %i #leq M_{R} < %i " % (z[k-1],x[i-1],x[i]) )
                    h_label_int_comp.append("%i b-tag, M_{R} #geq %i " % (z[k-1],x[i-1]) )
            h_Rsq_slice_components.append(h_Rsq_comp)
            h_Rsq_integral_components.append(h_Rsq_int_comp)
            h_data_Rsq_slice_components.append(h_data_Rsq_comp)
            h_data_Rsq_integral_components.append(h_data_Rsq_int_comp)
            h_Rsq_slice_component_labels.append(h_label_comp)
            h_Rsq_integral_component_labels.append(h_label_int_comp)
                        
    h_MR_components = []
    h_Rsq_components = []
    h_RsqMR_components = []
    h_RsqMR_fine_components = []
    h_data_MR_components = []
    h_data_Rsq_components = []
    h_data_RsqMR_components = []
    h_data_RsqMR_fine_components = []
    h_th1x_components = []
    h_data_th1x_components = []
    h_sig_th1x_components = []
    h_sig_MR_components = []
    h_sig_Rsq_components = []
    h_sig_RsqMR_components = []
    
    h_labels = []        
    h_colors = [rt.kOrange,rt.kViolet,rt.kRed,rt.kGreen]
    if len(z)>2:
        for k in range(1,len(z)):
            h_MR_components.append(h_nBtagRsqMR.ProjectionX("h_MR_%ibtag"%z[k-1],0,-1,k,k,""))
            h_Rsq_components.append(h_nBtagRsqMR.ProjectionY("h_Rsq_%ibtag"%z[k-1],0,-1,k,k,""))
            h_data_MR_components.append(h_data_nBtagRsqMR.ProjectionX("h_MR_data_%ibtag"%z[k-1],0,-1,k,k,""))
            h_data_Rsq_components.append(h_data_nBtagRsqMR.ProjectionY("h_Rsq_data_%ibtag"%z[k-1],0,-1,k,k,""))
            if doSignalInj:
                h_sig_MR_components.append(h_sig_nBtagRsqMR.ProjectionX("h_MR_sig_%ibtag"%z[k-1],0,-1,k,k,""))
                h_sig_Rsq_components.append(h_sig_nBtagRsqMR.ProjectionY("h_Rsq_sig_%ibtag"%z[k-1],0,-1,k,k,""))
            
            h_nBtagRsqMR.GetZaxis().SetRange(k,k)
            h_RsqMR_components.append(h_nBtagRsqMR.Project3D("%ibtag_yx"%z[k-1]))
            h_nBtagRsqMR_fine.GetZaxis().SetRange(k,k)
            h_RsqMR_fine_components.append(h_nBtagRsqMR_fine.Project3D("%ibtag_yx"%z[k-1]))
            if doSignalInj:
                h_sig_nBtagRsqMR.GetZaxis().SetRange(k,k)
                h_sig_RsqMR_components.append(h_sig_nBtagRsqMR.Project3D("%ibtag_yx"%z[k-1]))
            
            h_data_nBtagRsqMR.GetZaxis().SetRange(k,k)
            h_data_RsqMR_components.append(h_data_nBtagRsqMR.Project3D("%ibtag_yx"%z[k-1]))
            h_data_nBtagRsqMR_fine.GetZaxis().SetRange(k,k)
            h_data_RsqMR_fine_components.append(h_data_nBtagRsqMR_fine.Project3D("%ibtag_yx"%z[k-1]))
            
            if computeErrors:
                h_MR_components[-1] = getErrors1D(h_MR_components[-1],h_data_MR_components[-1],sysTree,options,"x",0,len(x)-1,0,len(y)-1,k,k,x,y,z)
                h_Rsq_components[-1] = getErrors1D(h_Rsq_components[-1],h_data_Rsq_components[-1],sysTree,options,"y",0,len(x)-1,0,len(y)-1,k,k,x,y,z)
                h_RsqMR_components[-1] = getErrors2D(h_RsqMR_components[-1],h_data_RsqMR_components[-1],sysTree,options,"yx",0,len(x)-1,0,len(y)-1,k,k,x,y,z)

            # convert after getting the 2D errors
            h_th1x_components.append(get1DHistoFrom2D(h_RsqMR_components[-1],x,y,'h_th1x_%ibtag'%(z[k-1])))
            h_data_th1x_components.append(get1DHistoFrom2D(h_data_RsqMR_components[-1],x,y,'h_th1x_data_%ibtag'%(z[k-1])))
            if doSignalInj:
                h_sig_th1x_components.append(get1DHistoFrom2D(h_sig_RsqMR_components[-1],x,y,'h_th1x_sig_%ibtag'%(z[k-1])))
            
            if z[k-1]==3 and z[-1]==4:
                h_labels.append("#geq %i b-tag" % z[k-1] )
            if z[k-1]==1 and z[-1]==4 and box in ['MuEle','MuMu','EleEle']:                
                h_labels.append("#geq %i b-tag" % z[k-1] )
            else:            
                h_labels.append("%i b-tag" % z[k-1] )
                
    h_RsqMR_residuals = h_data_RsqMR.Clone("h_RsqMR_residuals")
    h_RsqMR_residuals.Add(h_RsqMR,-1.)
    #h_RsqMR_percentdiff = h_data_RsqMR.Clone("h_RsqMR_percentdiff")
    #h_RsqMR_percentdiff.Add(h_RsqMR,-1.)
    #h_RsqMR_percentdiff.Divide(h_RsqMR)
    
    h_RsqMR_statnsigma = h_data_RsqMR.Clone("h_RsqMR_statnsigma")
    for i in range(1,h_RsqMR_statnsigma.GetNbinsX()+1):
        for j in range(1,h_RsqMR_statnsigma.GetNbinsY()+1):
            fit = h_RsqMR.GetBinContent(i,j)
            dat = h_data_RsqMR.GetBinContent(i,j)
            if fit > 0.0:
                h_RsqMR_statnsigma.SetBinContent(i,j,(dat-fit)/rt.TMath.Sqrt(fit))
            else:
                print "ERROR FIT = 0, SKIPPING BIN"
    
    h_RsqMR_nsigma = h_RsqMR_statnsigma.Clone("h_RsqMR_nsigma")
    if computeErrors:
        h_RsqMR_nsigma = getNsigma2D(h_RsqMR_nsigma,h_data_RsqMR,toyTree,options,"yx",0,len(x)-1,0,len(y)-1,0,len(z)-1,x,y,z)
        if options.plotRegion!='Full':            
            for i in range(1,h_RsqMR_nsigma.GetNbinsX()+1):
                for j in range(1,h_RsqMR_nsigma.GetNbinsY()+1):
                    w.var('MR').setVal(h_RsqMR_nsigma.GetXaxis().GetBinCenter(i))
                    w.var('Rsq').setVal(h_RsqMR_nsigma.GetYaxis().GetBinCenter(j))
                    w.var('nBtag').setVal(0.)
                    inSideband = 0
                    for fitname in options.plotRegion.split(','):
                        inSideband += ( w.var('MR').inRange(fitname) * w.var('Rsq').inRange(fitname) * w.var('nBtag').inRange(fitname) )
                    if not inSideband:
                        h_RsqMR_nsigma.SetBinContent(i,j,0)


    h_RsqMR_nsigma_components = []
    h_RsqMR_statnsigma_components = []
    #h_RsqMR_percentdiff_components = []
    h_RsqMR_residuals_components = []
    if len(z)>2:
        for k in range(1,len(z)):
            h_RsqMR_nsigma_btag = h_RsqMR_nsigma.Clone("h_RsqMR_nsigma_%ibtag"%z[k-1])
            if computeErrors:
                h_RsqMR_nsigma_btag = getNsigma2D(h_RsqMR_nsigma_btag,h_data_RsqMR_components[k-1],toyTree,options,"yx",0,len(x)-1,0,len(y)-1,k,k,x,y,z)
                if options.plotRegion!='Full':            
                    for i in range(1,h_RsqMR_nsigma_btag.GetNbinsX()+1):
                        for j in range(1,h_RsqMR_nsigma_btag.GetNbinsY()+1):
                            w.var('MR').setVal(h_RsqMR_nsigma_btag.GetXaxis().GetBinCenter(i))
                            w.var('Rsq').setVal(h_RsqMR_nsigma_btag.GetYaxis().GetBinCenter(j))
                            w.var('nBtag').setVal(z[k-1])
                            inSideband = 0
                            for fitname in options.plotRegion.split(','):
                                inSideband += ( w.var('MR').inRange(fitname) * w.var('Rsq').inRange(fitname) * w.var('nBtag').inRange(fitname) )
                            if not inSideband:
                                h_RsqMR_nsigma_btag.SetBinContent(i,j,0)                                
                h_RsqMR_nsigma_components.append(h_RsqMR_nsigma_btag)
                            
            h_RsqMR_statnsigma_btag = h_RsqMR_statnsigma.Clone("h_RsqMR_statnsigma_%ibtag"%z[k-1])
            for i in range(1,h_RsqMR_statnsigma_btag.GetNbinsX()+1):
                for j in range(1,h_RsqMR_statnsigma_btag.GetNbinsY()+1):
                    fit = h_RsqMR_components[k-1].GetBinContent(i,j)
                    dat = h_data_RsqMR_components[k-1].GetBinContent(i,j)
                    if fit > 0.0:
                        h_RsqMR_statnsigma_btag.SetBinContent(i,j,(dat-fit)/rt.TMath.Sqrt(fit))                        
                    else:
                        print "ERROR FIT = 0, SKIPPING BIN"
            h_RsqMR_statnsigma_components.append(h_RsqMR_statnsigma_btag)
            #h_RsqMR_percentdiff_btag = h_data_RsqMR_components[k-1].Clone("h_RsqMR_percentdiff_%ibtag"%z[k-1])
            #h_RsqMR_percentdiff_btag.Add(h_RsqMR_components[k-1],-1.)
            #h_RsqMR_percentdiff_btag.Divide(h_RsqMR_components[k-1])
            #h_RsqMR_percentdiff_components.append(h_RsqMR_percentdiff_btag)            
            h_RsqMR_residuals_btag = h_data_RsqMR_components[k-1].Clone("h_RsqMR_residuals_%ibtag"%z[k-1])
            h_RsqMR_residuals_btag.Add(h_RsqMR_components[k-1],-1.)
            h_RsqMR_residuals_components.append(h_RsqMR_residuals_btag)


    btagLabel = ""
    if z[-1] == z[0]+1 and z[-1]==4:
        btagLabel = "#geq %i b-tag" % z[0]
    elif z[-1] == z[0]+1:
        btagLabel = "%i b-tag" % z[0]
    elif z[-1]==4:
        btagLabel = "#geq %i b-tag" % z[0]    
    else:
        btagLabel = "%i-%i b-tag" % (z[0],z[-2])

    lumiLabel = "%.1f fb^{-1} (13 TeV)" % (lumi/1000)
    boxLabel = "razor %s %s %s Fit" % (box,btagLabel,fitRegion.replace('LowMR,LowRsq','Sideband'))
    #plotLabel = "%s Projection" % (plotRegion)
    plotLabel = ""

    sidebandFit = None
    if fitRegion=="LowMR,LowRsq":
        mrSide = w.var('MR').getMax('LowMR')
        rsqSide = w.var('Rsq').getMax('LowRsq')
        sidebandFit = [mrSide, rsqSide]
        
    if options.isData:
        dataString = "Data"
    else:
        dataString = "Sim. Data"
    

    for h in [h_nBtagRsqMR,h_data_nBtagRsqMR,h_RsqMR,h_data_RsqMR,h_MR,h_data_MR,h_Rsq,h_data_Rsq]:
        tdirectory.cd()
        h.Write()

    for  h in h_RsqMR_components:
        tdirectory.cd()
        h.Write()
        
    eventsLabel = "Events"
    if densityCorr:
        eventsLabel = "Events/Bin Width"

    h_sig_total_th1x_components = []
    if doSignalInj:
        h_MR_components.append(h_sig_MR)
        h_Rsq_components.append(h_sig_Rsq)
        h_sig_total_th1x_components.append(h_sig_th1x)
        
    print1DProj(c,tdirectory,h_th1x,h_data_th1x,options.outDir+"/h_th1x_%s.pdf"%box,"Bin Number",eventsLabel,lumiLabel,boxLabel,plotLabel,options.isData,doSignalInj,options,None,h_sig_total_th1x_components)
    print1DProj(c,tdirectory,h_MR,h_data_MR,options.outDir+"/h_MR_%s.pdf"%box,"M_{R} [GeV]",eventsLabel,lumiLabel,boxLabel,plotLabel,options.isData,doSignalInj,options,None,h_MR_components,h_colors,h_labels)
    print1DProj(c,tdirectory,h_Rsq,h_data_Rsq,options.outDir+"/h_Rsq_%s.pdf"%box,"R^{2}",eventsLabel,lumiLabel,boxLabel,plotLabel,options.isData,doSignalInj,options,None,h_Rsq_components,h_colors,h_labels)
    
    more_colors = [rt.kBlack,rt.kBlue]
    more_colors.extend(h_colors)
    more_colors.extend([rt.kMagenta,rt.kGray,rt.kCyan,rt.kYellow])
    #print1DSlice(c,tdirectory,h_MR_integrals,h_data_MR_integrals,options.outDir+"/h_MR_slicesRsq_%s.pdf"%(box),"M_{R} [GeV]",eventsLabel,lumiLabel,boxLabel,plotLabel,options.isData,None,more_colors,h_MR_integral_labels)
    #print1DSlice(c,tdirectory,h_Rsq_integrals,h_data_Rsq_integrals,options.outDir+"/h_Rsq_slicesMR_%s.pdf"%(box),"R^{2}",eventsLabel,lumiLabel,boxLabel,plotLabel,options.isData,None,more_colors,h_Rsq_integral_labels)
    
    if len(z)>2:
        for k in range(1,len(z)):            
            newBoxLabel = "razor %s %s %s Fit"%(box,h_labels[k-1],fitRegion.replace('LowMR,LowRsq','Sideband'))
            #print1DSlice(c,tdirectory,zip(*h_MR_integral_components)[k-1],zip(*h_data_MR_integral_components)[k-1],options.outDir+"/h_MR_slicesRsq_%ibtag_%s.pdf"%(z[k-1],box),"M_{R} [GeV]",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData,None,more_colors,zip(*h_MR_integral_component_labels)[k-1])
            #print1DSlice(c,tdirectory,zip(*h_Rsq_integral_components)[k-1],zip(*h_data_Rsq_integral_components)[k-1],options.outDir+"/h_Rsq_slicesMR_%ibtag_%s.pdf"%(z[k-1],box),"R^{2}",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData,None,more_colors,zip(*h_Rsq_integral_component_labels)[k-1])

    #print2DResiduals(c,tdirectory,h_RsqMR_residuals,options.outDir+"/h_RsqMR_residuals_log_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Residuals (%s - Fit)"%dataString,lumiLabel,boxLabel,plotLabel,x,y,options.isData,sidebandFit,doSignalInj,options)
    #print2DResiduals(c,tdirectory,h_RsqMR_percentdiff,options.outDir+"/h_RsqMR_percentdiff_log_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Percent Diff. (%s - Fit)/Fit"%dataString,lumiLabel,boxLabel,plotLabel,x,y,options.isData,sidebandFit,doSignalInj,options)
    print2DResiduals(c,tdirectory,h_RsqMR_statnsigma,options.outDir+"/h_RsqMR_statnsigma_log_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Stat. n#sigma (%s - Fit)/sqrt(Fit)"%dataString,lumiLabel,boxLabel,plotLabel,x,y,options.isData,sidebandFit,doSignalInj,options)
    if computeErrors:
        print2DResiduals(c,tdirectory,h_RsqMR_nsigma,options.outDir+"/h_RsqMR_nsigma_log_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Stat.+Sys. n#sigma",lumiLabel,boxLabel,plotLabel,x,y,options.isData,sidebandFit,doSignalInj,options)
    #print2DScatter(c,tdirectory,h_RsqMR_fine,options.outDir+"/h_RsqMR_scatter_log_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Fit",lumiLabel,boxLabel,plotLabel,x,y,h_data_RsqMR_fine.GetMinimum(),h_data_RsqMR_fine.GetMaximum(),options.isData,sidebandFit,doSignalInj,options)
    print2DScatter(c,tdirectory,h_data_RsqMR_fine,options.outDir+"/h_RsqMR_scatter_data_log_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", dataString,lumiLabel,boxLabel,plotLabel,x,y,h_data_RsqMR_fine.GetMinimum(),h_data_RsqMR_fine.GetMaximum(),options.isData,sidebandFit,doSignalInj,options)

    if len(z)>2:
        for k in range(0,len(z)-1):
            newBoxLabel = "razor %s %s %s Fit"%(box,h_labels[k],fitRegion.replace('LowMR,LowRsq','Sideband'))
            if doSignalInj:
                #print1DProj(c,tdirectory,h_th1x_components[k],h_data_th1x_components[k],options.outDir+"/h_th1x_%ibtag_%s.pdf"%(z[k],box),"Bin Number",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData,doSignalInj,options,None,[h_sig_th1x_components[k]])
                print1DProj(c,tdirectory,h_MR_components[k],h_data_MR_components[k],options.outDir+"/h_MR_%ibtag_%s.pdf"%(z[k],box),"M_{R} [GeV]",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData,doSignalInj,options,None,[h_sig_MR_components[k]])
                print1DProj(c,tdirectory,h_Rsq_components[k],h_data_Rsq_components[k],options.outDir+"/h_Rsq_%ibtag_%s.pdf"%(z[k],box),"R^{2}",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData,doSignalInj,options,None,[h_sig_Rsq_components[k]])
            else:
                #print1DProj(c,tdirectory,h_th1x_components[k],h_data_th1x_components[k],options.outDir+"/h_th1x_%ibtag_%s.pdf"%(z[k],box),"Bin Number",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData)
                print1DProj(c,tdirectory,h_MR_components[k],h_data_MR_components[k],options.outDir+"/h_MR_%ibtag_%s.pdf"%(z[k],box),"M_{R} [GeV]",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData)
                print1DProj(c,tdirectory,h_Rsq_components[k],h_data_Rsq_components[k],options.outDir+"/h_Rsq_%ibtag_%s.pdf"%(z[k],box),"R^{2}",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData)
            if computeErrors:
                if doSignalInj:
                    print1DProjNs(c,tdirectory,h_th1x_components[k],h_data_th1x_components[k],h_RsqMR_nsigma_components[k],options.outDir+"/h_th1x_ns_%ibtag_%s.pdf"%(z[k],box),"Bin Number",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData,doSignalInj,options,None,[h_sig_th1x_components[k]], cfg=cfg)
                else:
                    print1DProjNs(c,tdirectory,h_th1x_components[k],h_data_th1x_components[k],h_RsqMR_nsigma_components[k],options.outDir+"/h_th1x_ns_%ibtag_%s.pdf"%(z[k],box),"Bin Number",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData,doSignalInj,options, cfg=cfg)
                print2DResiduals(c,tdirectory,h_RsqMR_nsigma_components[k],options.outDir+"/h_RsqMR_nsigma_log_%ibtag_%s.pdf"%(z[k],box),"M_{R} [GeV]", "R^{2}", "Stat.+Sys. n#sigma",lumiLabel,newBoxLabel,plotLabel,x,y,options.isData,sidebandFit,doSignalInj,options)
            else:                
                if doSignalInj:
                    print1DProj(c,tdirectory,h_th1x_components[k],h_data_th1x_components[k],options.outDir+"/h_th1x_%ibtag_%s.pdf"%(z[k],box),"Bin Number",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData,doSignalInj,options,None,[h_sig_th1x_components[k]])
                else:
                    print1DProj(c,tdirectory,h_th1x_components[k],h_data_th1x_components[k],options.outDir+"/h_th1x_%ibtag_%s.pdf"%(z[k],box),"Bin Number",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData,doSignalInj,options)
            #print2DResiduals(c,tdirectory,h_RsqMR_residuals_components[k],options.outDir+"/h_RsqMR_residuals_log_%ibtag_%s.pdf"%(z[k],box),"M_{R} [GeV]", "R^{2}", "Residuals (%s - Fit)"%dataString,lumiLabel,newBoxLabel,plotLabel,x,y,options.isData,sidebandFit,doSignalInj,options)
            #print2DResiduals(c,tdirectory,h_RsqMR_percentdiff_components[k],options.outDir+"/h_RsqMR_percentdiff_log_%ibtag_%s.pdf"%(z[k],box),"M_{R} [GeV]", "R^{2}", "Percent Diff. (%s - Fit)/Fit"%dataString,lumiLabel,newBoxLabel,plotLabel,x,y,options.isData)
            print2DResiduals(c,tdirectory,h_RsqMR_statnsigma_components[k],options.outDir+"/h_RsqMR_statnsigma_log_%ibtag_%s.pdf"%(z[k],box),"M_{R} [GeV]", "R^{2}", "Stat. n#sigma (%s - Fit)/sqrt(Fit)"%dataString,lumiLabel,newBoxLabel,plotLabel,x,y,options.isData,sidebandFit,doSignalInj,options)
            #print2DScatter(c,tdirectory,h_RsqMR_fine_components[k],options.outDir+"/h_RsqMR_scatter_log_%ibtag_%s.pdf"%(z[k],box),"M_{R} [GeV]", "R^{2}", "Fit",lumiLabel,newBoxLabel,plotLabel,x,y,h_data_RsqMR_fine_components[k].GetMinimum(),h_data_RsqMR_fine_components[k].GetMaximum(),options.isData,sidebandFit,doSignalInj,options)
            print2DScatter(c,tdirectory,h_data_RsqMR_fine_components[k],options.outDir+"/h_RsqMR_scatter_data_log_%ibtag_%s.pdf"%(z[k],box),"M_{R} [GeV]", "R^{2}", dataString,lumiLabel,newBoxLabel,plotLabel,x,y,h_data_RsqMR_fine_components[k].GetMinimum(),h_data_RsqMR_fine_components[k].GetMaximum(),options.isData,sidebandFit,doSignalInj,options)


        
    #for j in range(0,len(y)-1):
    #    newBoxLabel = "razor %s %s %s Fit"%(box,h_MR_slice_labels[j],fitRegion.replace('LowMR,LowRsq','Sideband'))
    #    print1DProj(c,tdirectory,h_MR_slices[j],h_data_MR_slices[j],options.outDir+"/h_MR_%.2fRsq%.2f_%s.pdf"%(y[j],y[j+1],box),"M_{R} [GeV]",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData,None,h_MR_slice_components[j],h_colors,h_labels)
    #    for k in range(0,len(z)-1):
    #        newBoxLabel = "razor %s %s %s Fit"%(box,h_MR_slice_component_labels[j][k],fitRegion.replace('LowMR,LowRsq','Sideband'))
    #        print1DProj(c,tdirectory,h_MR_slice_components[j][k],h_data_MR_slice_components[j][k],options.outDir+"/h_MR_%ibtag_%.2fRsq%.2f_%s.pdf"%(z[k],y[j],y[j+1],box),"M_{R} [GeV]",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData)
    #    
    #for i in range(0,len(x)-1):
    #    newBoxLabel = "razor %s %s %s Fit"%(box,h_Rsq_slice_labels[i],fitRegion.replace('LowMR,LowRsq','Sideband'))
    #    print1DProj(c,tdirectory,h_Rsq_slices[i],h_data_Rsq_slices[i],options.outDir+"/h_Rsq_%iMR%i_%s.pdf"%(x[i],x[i+1],box),"R^{2}",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData,None,h_Rsq_slice_components[i],h_colors,h_labels)
    #    for k in range(0,len(z)-1):
    #        newBoxLabel = "razor %s %s %s Fit"%(box,h_Rsq_slice_component_labels[i][k],fitRegion.replace('LowMR,LowRsq','Sideband'))
    #        print1DProj(c,tdirectory,h_Rsq_slice_components[i][k],h_data_Rsq_slice_components[i][k],options.outDir+"/h_Rsq_%ibtag_%iMR%i_%s.pdf"%(z[k],x[i],x[i+1],box),"R^{2}",eventsLabel,lumiLabel,newBoxLabel,plotLabel,options.isData)

            
