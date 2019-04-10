import os
import ROOT as rt
import copy
import array
import random

#local imports
from PlotFit import setFFColors
from RunCombine import exec_me
from razorWeights import applyMTUncertainty1D, applyMTUncertainty2D
import macro as macro

def makeLegend(hists, titles, ordering, x1=0.75, y1=0.6, x2=0.9, y2=0.9):
    """Takes a dict of histograms, a dict of histogram titles, and an ordered list of names, and returns a legend with the histograms in the desired order"""
    leg = rt.TLegend(x1, y1, x2, y2)
    rt.SetOwnership(leg, False)
    for name in ordering: 
        leg.AddEntry(hists[name], titles[name], "f")
    return leg

def makeStack(hists, ordering, title="Stack"):
    """Takes a dict of histograms and an ordered list of names, and returns a THStack containing the histograms stacked in the desired order"""
    stack = rt.THStack("thstack"+title.replace(" ",""), title)
    rt.SetOwnership(stack, False)
    for name in ordering: 
        if name in hists:
            stack.Add(hists[name])
        else: 
            print("Warning in makeStack: histogram "+name+" not found in histogram dictionary!")
    return stack

def makeGrayGraphs(hNS):    
    if hNS is None or hNS == 0: return
    fGrayGraphs = []
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
    return fGrayGraphs

def drawUnrolledBinMapping(c, unrollBins, xtitle="", ytitle="", printstr="", printdir="."):
    """Draw a 2D histogram showing the bin boundaries for the indicated TH2Poly-style binning
    unrollBins should be a tuple (xbins, columns)"""

    #make temporary TH2Poly
    poly = macro.makeTH2PolyFromColumns("Temp", "", unrollBins[0], unrollBins[1])
    #set bin content = bin number
    for bn in range(1, poly.GetNumberOfBins()+1):
        poly.SetBinContent(bn, bn-1)
    poly.SetBinContent(1, 0.1) #hack to get text to display

    #draw it
    draw2DHist(c, poly, xtitle=xtitle, ytitle=ytitle, printstr=printstr, logz=False, numDigits=0, textSize=2, printdir=printdir, drawCMSPreliminary=False, drawColz=False)

    #delete temp histogram
    poly.Delete()
        

def getLinesForUnrolled(hist):        
    # the gray lines
    lines = []
    nBinsX = hist.GetNbinsX()
    nBinsY = hist.GetNbinsY()
    for i in range(1,nBinsX):
        lineX = i*nBinsY
        lines.append(rt.TLine(lineX, 0, lineX, hist.GetYaxis().GetXmax()))
        lines[-1].SetLineStyle(2)
        lines[-1].SetLineColor(rt.kBlack)
    return lines

def setHistColor(hist, name):
    colors = {"WJets":rt.kRed+1, "WJetsInv":rt.kRed+1, "DYJets":rt.kBlue+1, "DYJetsInv":rt.kBlue+1, "TTJets":rt.kGreen+2, "TTJets1L":rt.kGreen+2, "TTJets2L":rt.kGreen+3, "ZInv":rt.kCyan+1, "QCD":rt.kMagenta, "SingleTop":rt.kOrange-3, "VV":rt.kViolet+3, "TTV":rt.kGreen-7, "DYJetsLow":rt.kBlue+1, "GJets":rt.kOrange, "GJetsInv":rt.kOrange, "GJetsFrag":(rt.kOrange+4), "Other":rt.kAzure+4}
    """Sets histogram color"""
    if name in colors: hist.SetFillColor(colors[name])
    else: 
        hist.SetFillColor(int(random.random()*rt.gStyle.GetNumberOfColors()))


def plot_several(c, hists=0, leg=0, colors=[], xtitle="", ytitle="Events", ymin=None, ymax=None, printstr="hist", logx=False, logy=True, lumistr="40 pb^{-1}", commentstr="", ratiomin=0.5, ratiomax=1.5, saveroot=True, savepdf=True, savepng=True, savec=True, printdir='.'):
    """Draw several histograms as colored lines"""
    #setup
    c.Clear()
    c.cd()
    if len(hists) < 1: 
        print "Error in plot_several: no histograms provided!"
        return
    elif len(hists) > 1: pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.4, 1, 1)
    else: pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.1, 1, 1)
    pad1.SetBottomMargin(0)
    pad1.SetLogx(logx)
    pad1.SetLogy(logy)
    pad1.Draw()
    pad1.cd()

    #draw 
    for i in range(len(hists)-1,-1,-1): #draw in reverse order
        hist = hists[i]
        if i < len(colors):
            hist.SetMarkerColor(colors[i])
        hist.SetLineWidth(2)
        hist.SetTitle("")
        hist.GetYaxis().SetTitle(ytitle)
        hist.GetYaxis().SetLabelSize(0.03)
        hist.SetStats(0)
        if logy: 
            hist.GetXaxis().SetMoreLogLabels()
            hist.GetXaxis().SetNoExponent()
        if len(hists) > 1:
            hist.GetYaxis().SetTitleOffset(0.50)
        else:
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitleOffset(0.45)
        hist.GetYaxis().SetTitleSize(0.05)
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(0.5)
        if ymin is not None: hist.SetMinimum(ymin)
        if ymax is not None: hist.SetMaximum(ymax)
        if i == len(hists)-1:
            hist.Draw('p')
        else:
            hist.Draw("psame")
    pad1.Modified()
    rt.gPad.Update()

    #add legend and LaTeX 
    leg.Draw()
    t1 = rt.TLatex(0.1,0.94, "CMS preliminary")
    t2 = rt.TLatex(0.55,0.94, ((lumistr != "")*((lumistr)+' ('))+'13 TeV'+((lumistr != "")*(')')))
    t1.SetNDC()
    t2.SetNDC()
    t1.SetTextSize(0.06)
    t2.SetTextSize(0.06)
    t1.Draw()
    t2.Draw()
    if commentstr != "":
        t3 = rt.TLatex(0.40, 0.92, commentstr)
        t3.SetNDC()
        t3.SetTextSize(0.05)
        t3.Draw()

    #ratio histograms
    c.cd()
    pad2 = rt.TPad(printstr+"pad2",printstr+"pad2",0,0.0,1,0.4)
    pad2.SetTopMargin(0)
    pad2.SetTopMargin(0.008)
    pad2.SetBottomMargin(0.25)
    pad2.SetGridy()
    pad2.SetLogx(logx)
    lowerHists = []
    for i in range(len(hists)-1,-1,-1): #draw in reverse order
        if i == 0: continue
        hist = hists[i]
        lowerHists.append(make1DRatioHistogram(hist, hists[0], xtitle, ratiomin, ratiomax, logx))
        lowerHists[-1].SetMarkerSize(0.5)
        lowerHists[-1].SetMarkerStyle(20)
        lowerHists[-1].GetYaxis().SetTitle("Ratio")

        if i == len(hists)-1:
            pad2.Draw()
            pad2.cd()
            lowerHists[-1].Draw('p') 
        else:
            lowerHists[-1].Draw('psame')
        pad2.Modified()
        rt.gPad.Update()

    #save
    if savepng: c.Print(printdir+'/'+printstr+".png")
    if savepdf: c.Print(printdir+'/'+printstr+".pdf")
    if saveroot: c.Print(printdir+'/'+printstr+".root")
    if savec: c.Print(printdir+'/'+printstr+".C")

def plot_basic(c, mc=0, data=0, fit=0, leg=0, xtitle="", ytitle="Events", ymin=None, ymax=None, printstr="hist", logx=False, logy=True, lumistr="40 pb^{-1}", commentstr="", ratiomin=0.5, ratiomax=1.5, pad2Opt="Ratio", fitColor=rt.kBlue, mcErrColor=rt.kBlack, customPad2Hist=None, saveroot=True, savepdf=True, savepng=True, savec=True, printdir='.', grayLines=None):
    """Plotting macro with options for data, MC, and fit histograms.  Creates data/MC ratio if able."""
    #setup
    c.Clear()
    c.cd()
    if (data and mc) or (data and fit) or (mc and fit): pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.4, 1, 1)
    else: pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.1, 1, 1)
    rt.SetOwnership(pad1, False)
    pad1.SetBottomMargin(0)
    pad1.SetLogx(logx)
    pad1.SetLogy(logy)
    pad1.Draw()
    pad1.cd()

    #draw MC
    if mc:
        mc.SetTitle("")
        #make total MC histogram
        histList = mc.GetHists()
        mcTotal = histList.First().Clone()
        mcTotal.SetTitle("")
        mcTotal.SetStats(0)
        mcTotal.Reset()
        numMCHists = 0
        for h in histList:
            mcTotal.Add(h)
            numMCHists += 1
        mcTotal.SetFillColor(mcErrColor)
        mcTotal.GetYaxis().SetTitle(ytitle)
        mcTotal.GetYaxis().SetTitleSize(0.1)
        mcTotal.GetYaxis().SetLabelSize(0.06)
        if logx: 
            mcTotal.GetXaxis().SetMoreLogLabels()
            mcTotal.GetXaxis().SetNoExponent()
        if not data: mcTotal.GetXaxis().SetTitle(xtitle)
        mcTotal.GetYaxis().SetTitleOffset(0.6)
        mcTotal.GetYaxis().SetTitleSize(0.05)
        if ymin is not None: mcTotal.SetMinimum(ymin)
        if ymax is not None: mcTotal.SetMaximum(ymax)
        if data and data.GetMaximum() > mcTotal.GetMaximum() and ymax is None: 
            mcTotal.SetMaximum(data.GetMaximum())
        mc.Draw('hist')
        if data and data.GetMaximum() > mc.GetMaximum() and ymax is None: 
            mc.SetMaximum(data.GetMaximum())
        if ymin is not None: mc.SetMinimum(ymin)
        if ymax is not None: mc.SetMaximum(ymax)
        if logx: 
            mc.GetXaxis().SetMoreLogLabels()
            mc.GetXaxis().SetNoExponent()
        mc.GetXaxis().SetTitle(xtitle)
        mc.GetXaxis().SetTitleSize(0.06)
        mc.GetYaxis().SetTitle(ytitle)
        mc.GetYaxis().SetTitleOffset(0.60)
        mc.GetYaxis().SetTitleSize(0.06)
        mc.GetYaxis().SetLabelSize(0.06)
        if fit and numMCHists == 1: 
            mcTotal.SetLineWidth(2)
            mcTotal.SetFillStyle(0)
            mcTotal.SetLineColor(rt.kBlue)
            mcTotal.SetFillColor(rt.kBlue)

            #Draw MC as points
            #mcTotal.SetLineColor(rt.kBlack)
            #mcTotal.SetMarkerStyle(21)
            #mcTotal.SetMarkerSize(1)
            #mcTotal.SetMarkerColor(rt.kBlack)
            #mcTotal.GetXaxis().SetTitle(xtitle)
            #mcTotal.GetYaxis().SetTitle(ytitle)
            #mcTotal.GetYaxis().SetTitleOffset(0.60)
            #mcTotal.GetYaxis().SetTitleSize(0.06)
            #mcTotal.GetYaxis().SetLabelSize(0.06)
            #if ymin is not None: 
            #    mcTotal.SetMinimum(ymin)
            #if ymax is not None: 
            #    mcTotal.SetMaximum(ymax)
            #mcTotal.Draw()
            #if logx: mcTotal.GetXaxis().SetMoreLogLabels()

            #Draw MC as blue line with lines for up/down errors
            mcTotalCopy = mcTotal.Clone()
            mcTotalCopy.SetFillStyle(0)
            mcTotalCopy.SetLineWidth(2)
            mcLower = mcTotal.Clone()
            mcLower.SetLineColor(rt.kBlue+2)
            mcTotal.SetFillColor(rt.kAzure+6)
            mcTotal.SetFillStyle(3001)
            mcLower.SetLineWidth(1)
            mcUpper = mcLower.Clone()
            for bx in range(1, mcTotal.GetNbinsX()+1):
                mcLower.SetBinContent(bx, mcTotal.GetBinContent(bx)-mcTotal.GetBinError(bx))
                mcUpper.SetBinContent(bx, mcTotal.GetBinContent(bx)+mcTotal.GetBinError(bx))
                mcTotalBin = mcTotal.GetBinContent(bx)
                mcTotalErr = mcTotal.GetBinError(bx)
                if mcTotalBin < ymin and mcTotalBin+mcTotalErr > ymin: #push it into the visible range
                    mcTotal.SetBinContent(bx, ymin)
                    mcTotal.SetBinError(bx, mcTotalErr + mcTotal.GetBinContent(bx) - mcTotalBin)
            mcLower.Draw('hist')
            mcUpper.Draw('histsame')
            mcTotal.Draw('e2same')
            mcTotalCopy.Draw("histsame") 

            #for legend
            mcHist = histList.First()
            mcHist.SetLineWidth( mcTotal.GetLineWidth() )
            mcHist.SetLineColor( mcTotal.GetLineColor() )
            mcHist.SetFillColor( mcTotal.GetFillColor() )
            mcHist.SetFillStyle( mcTotal.GetFillStyle() )
            
        if numMCHists > 1:
            mcTotal.SetFillStyle(3001)
            mcTotal.SetFillColor(rt.kGray+2)
            mcTotal.Draw("e2same")

    #draw fit
    if fit:
        fit.SetStats(0)
        fit.SetMarkerStyle(20)
        fit.SetLineWidth(2)
        fit.SetMarkerSize(0)
        fit.SetLineColor(fitColor)
        fit.SetMarkerColor(fitColor)
        fit.GetYaxis().SetTitle(ytitle)
        fit.SetTitle("")
        if ymin is not None: fit.SetMinimum(ymin)
        if ymax is not None: fit.SetMaximum(ymax)
        #blue line for fit
        fitCopy = fit.Clone()
        fitCopy.SetFillStyle(0)        
        fitCopy.SetLineWidth(2)

        if data:
            #Fit = blue line with shading for error band
            fit.SetFillStyle(3001)
            fit.SetFillColor(rt.kAzure+6)
            fit.Draw("e2same")
            fitCopy.Draw("histsame") 
        elif mc:
            #Fit = blue line with lines for up/down errors
            #fitLower = fit.Clone()
            #fitLower.SetLineColor(rt.kBlue+2)
            #fit.SetFillColor(rt.kAzure+6)
            #fit.SetFillStyle(3001)
            #fitLower.SetLineWidth(1)
            #fitUpper = fitLower.Clone()
            #for bx in range(1, fit.GetNbinsX()+1):
            #    fitLower.SetBinContent(bx, fit.GetBinContent(bx)-fit.GetBinError(bx))
            #    fitUpper.SetBinContent(bx, fit.GetBinContent(bx)+fit.GetBinError(bx))
            #fitLower.Draw('histsame')
            #fitUpper.Draw('histsame')
            #if numMCHists == 1:
            #    #fitCopy.SetMarkerSize(1)
            #    #fitCopy.SetLineWidth(1)
            #    fit.Draw('e2same')
            #    mcHist.Draw('same')
            #fitCopy.Draw("histsame") 
            
            #Fit = black points
            fit.SetMarkerColor(rt.kBlack)
            fit.SetLineColor(rt.kBlack)
            fit.SetMarkerSize(1)
            fit.Draw('pe1same')
        else:
            fitCopy.Draw("histsame")
    #draw data
    if data:
        data.SetStats(0)
        data.SetMarkerStyle(20)
        data.SetMarkerSize(1)
        data.SetLineWidth(2)
        data.SetLineColor(rt.kBlack)
        data.GetYaxis().SetTitle(ytitle)
        data.SetTitle("")
        if ymin is not None: data.SetMinimum(ymin)
        if ymax is not None: data.SetMaximum(ymax)
        data.SetBinErrorOption(rt.TH1.kPoisson)
        data.Draw("pe1same")
    pad1.Modified()
    rt.gPad.Update()

    #gray lines
    if grayLines is not None:
        for line in grayLines: 
            line.Draw()

    #add legend and LaTeX 
    leg.Draw()
    CMS_lumi(pad1, lumistr=lumistr)

    if commentstr != "":
        t3 = rt.TLatex(0.40, 0.92, commentstr)
        t3.SetNDC()
        t3.SetTextSize(0.05)
        t3.Draw()

    #lower pad plot
    lowerPadHist = None
    lowerPadHist2 = None
    lowerPadHist2Central = None
    lowerPadHist2Lower = None
    lowerPadHist2Upper = None

    #set up custom histogram if provided
    if customPad2Hist is not None:
        customPad2Hist.SetTitle("")
        customPad2Hist.GetXaxis().SetTitle(xtitle)
        customPad2Hist.SetMinimum(ratiomin)
        customPad2Hist.SetMaximum(ratiomax)
        customPad2Hist.SetStats(0)
        customPad2Hist.GetXaxis().SetLabelSize(0.06)
        customPad2Hist.GetYaxis().SetLabelSize(0.06)
        customPad2Hist.GetYaxis().SetTitleOffset(0.35)
        customPad2Hist.GetXaxis().SetTitleOffset(1.50)
        customPad2Hist.GetYaxis().SetTitleSize(0.08)
        customPad2Hist.GetXaxis().SetTitleSize(0.08)
        customPad2Hist.SetTitle("")
        if logx: 
            customPad2Hist.GetXaxis().SetMoreLogLabels()
            customPad2Hist.GetXaxis().SetNoExponent()

    #make ratio data/MC
    if data and mc:
        if customPad2Hist is not None:
            lowerPadHist = customPad2Hist.Clone("lowerPadHist")
            lowerPadHist.GetYaxis().SetTitle("Nsigma")
            lowerPadHist.SetMarkerStyle(20)
            lowerPadHist.SetMarkerSize(1)
        elif pad2Opt.lower() == "ratio":
            #draw data/MC with poisson errors from data
            lowerPadHist = make1DRatioHistogram(data, mcTotal, xtitle, ratiomin, ratiomax, logx, ignoreDenominatorErrs=True)
            lowerPadHist.GetYaxis().SetTitle("Data / MC")
            #draw MC uncertainties on the ratio
            #get plot style from MC histogram
            dataWithMCStyle = mcTotal.Clone()
            for bx in range(1, dataWithMCStyle.GetNbinsX()+1):
                dataWithMCStyle.SetBinContent(bx, data.GetBinContent(bx))
                dataWithMCStyle.SetBinError(bx, data.GetBinError(bx))
            lowerPadHist2 = make1DRatioHistogram(dataWithMCStyle, mcTotal, xtitle, ratiomin, ratiomax, logx, ignoreNumeratorErrs=True)
            lowerPadHist2.GetYaxis().SetTitle(lowerPadHist.GetYaxis().GetTitle())
            for bx in range(1, lowerPadHist2.GetNbinsX()+1):
                lowerPadHist2.SetBinContent(bx,1)
            lowerPadHist2.SetFillStyle(3001)
            lowerPadHist2.SetFillColor(rt.kGray+2)
        elif pad2Opt.lower() == "nsigma" or pad2Opt.lower() == "pulls" or pad2Opt.lower() == "ff":
            lowerPadHist = make1DPullHistogram(data, mcTotal, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("(Data - MC)/#sigma")
        elif pad2Opt.lower() == "perc" or pad2Opt.lower() == "percent":
            lowerPadHist = make1DPercentDiffHistogram(data, mcTotal, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("(Data - MC)/MC")
    #make ratio data/fit
    elif data and fit:
        if customPad2Hist is not None:
            lowerPadHist = customPad2Hist.Clone("lowerPadHist")
            lowerPadHist.GetYaxis().SetTitle("Nsigma (Stat+Sys)")
            lowerPadHist.SetMarkerStyle(20)
            lowerPadHist.SetMarkerSize(1)
        elif pad2Opt.lower() == "ratio":
            lowerPadHist = make1DRatioHistogram(data, fit, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("Data / Fit")
        elif pad2Opt.lower() == "nsigma" or pad2Opt.lower() == "pulls" or pad2Opt.lower() == "ff":
            lowerPadHist = make1DPullHistogram(data, fit, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("(Data - Fit)/#sigma")
        elif pad2Opt.lower() == "perc" or pad2Opt.lower() == "percent":
            lowerPadHist = make1DPercentDiffHistogram(data, fit, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("(Data - Fit)/Fit")
    #make ratio mc/fit
    elif mc and fit:
        if customPad2Hist is not None:
            lowerPadHist = customPad2Hist.Clone("lowerPadHist")
            lowerPadHist.GetYaxis().SetTitle("Nsigma")
            lowerPadHist.SetMarkerStyle(fit.GetMarkerStyle())
            lowerPadHist.SetMarkerSize(fit.GetMarkerSize())
        elif pad2Opt.lower() == "ratio":
            if numMCHists == 1: #total MC plot
                lowerPadHist = make1DRatioHistogram(fit, mcTotal, xtitle, ratiomin, ratiomax, logx, ignoreDenominatorErrs=True)
                lowerPadHist.GetYaxis().SetTitle("Fit / MC")
                #draw relative MC uncertainties 
                lowerPadHist2 = make1DRatioHistogram(mcTotal, mcTotal, xtitle, ratiomin, ratiomax, logx, ignoreDenominatorErrs=True)
                lowerPadHist2.GetYaxis().SetTitle(lowerPadHist.GetYaxis().GetTitle())
                lowerPadHist2Upper = lowerPadHist2.Clone()
                lowerPadHist2Upper.SetFillStyle(0)
                lowerPadHist2Upper.SetLineWidth(1)
                lowerPadHist2Lower = lowerPadHist2Upper.Clone()
                lowerPadHist2Central = lowerPadHist2Upper.Clone()
                lowerPadHist2Central.SetLineWidth(2)
                for bx in range(1, lowerPadHist2.GetNbinsX()+1):
                    lowerPadHist2Upper.SetBinContent(bx, 1+lowerPadHist2.GetBinError(bx))
                    lowerPadHist2Lower.SetBinContent(bx, 1-lowerPadHist2.GetBinError(bx))
            else: #MC stacked plot
                lowerPadHist = make1DRatioHistogram(fit, mcTotal, xtitle, ratiomin, ratiomax, logx)
                lowerPadHist.GetYaxis().SetTitle("Fit / MC")
                lowerPadHist.SetMarkerSize(1)
                lowerPadHist.SetMarkerColor(rt.kBlack)
        elif pad2Opt.lower() == "nsigma" or pad2Opt.lower() == "pulls" or pad2Opt.lower() == "ff":
            lowerPadHist = make1DPullHistogram(fit, mcTotal, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("(Fit - MC)/#sigma")
        elif pad2Opt.lower() == "perc" or pad2Opt.lower() == "percent":
            lowerPadHist = make1DPercentDiffHistogram(fit, mcTotal, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("(Fit - MC)/Fit")
    c.cd()
    if (data and mc) or (data and fit) or (mc and fit): 
        pad2 = rt.TPad(printstr+"pad2",printstr+"pad2",0,0.0,1,0.4)
        pad2.SetTopMargin(0)
        pad2.SetTopMargin(0.008)
        pad2.SetBottomMargin(0.25)
        pad2.SetGridy()
        pad2.SetLogx(logx)
        if lowerPadHist2 is not None:
            pad2.Draw()
            pad2.cd()
            lowerPadHist2.Draw("e2same")
            pad2.Modified()
            rt.gPad.Update()
        if lowerPadHist is not None:
            pad2.Draw()
            pad2.cd()
            if mc and fit:
                lowerPadHist.Draw("pe1same")
            else:
                lowerPadHist.Draw("pesame")
            pad2.Modified()
            rt.gPad.Update()
        if lowerPadHist2Central is not None:
            lowerPadHist2Upper.Draw("histsame")
            lowerPadHist2Lower.Draw("histsame")
            lowerPadHist2.Draw('e2same')
            lowerPadHist2Central.Draw("histsame")
            if mc and fit:
                lowerPadHist.Draw("pe1same")
            else:
                lowerPadHist.Draw("pesame")

    #save
    if savepng: c.Print(printdir+'/'+printstr+".png")
    if savepdf: c.Print(printdir+'/'+printstr+".pdf")
    if saveroot: c.Print(printdir+'/'+printstr+".root")
    if savec: c.Print(printdir+'/'+printstr+".C")

    pad1.Delete()

def draw2DHist(c, hist, xtitle="", ytitle="", ztitle="", zmin=None, zmax=None, printstr="hist", logx=True, logy=True, logz=True, lumistr="", commentstr="", dotext=True, drawErrs=False, palette='RAINBOW', grayGraphs=None, saveroot=True, savepdf=True, savepng=True, savec=True, numDigits=1, textSize=2.0, printdir='.', drawCMSPreliminary=True, drawColz=True, commentX=0.40, commentY=0.92):
    """Draw a single 2D histogram and print to file"""
    rt.gStyle.SetNumberContours(99)
    if palette == "FF":
        if zmin is None or zmax is None:
            setFFColors(hist, -5.1, 5.1)
        else:
            setFFColors(hist, zmin, zmax)
    elif palette == "RAINBOW":
        pass
    else:
        rt.gStyle.SetPalette(palette)
    c.Clear()
    c.cd()
    c.SetLogx(logx)
    c.SetLogy(logy)
    c.SetLogz(logz)
    c.Draw()
    c.cd()
    hist.SetTitle("")
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetYaxis().SetTitle(ytitle)
    #hist.GetZaxis().SetTitle(ztitle)
    hist.GetZaxis().SetTitle("") #until we can get the z-axis to display correctly
    hist.GetYaxis().SetLabelSize(0.03)
    hist.GetYaxis().SetTitleOffset(0.50)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.SetStats(0)
    if zmin is not None: hist.SetMinimum(zmin)
    elif hist.GetMinimum() < -10: hist.SetMinimum(0.1) #avoid drawing -999's etc
    if zmax is not None: hist.SetMaximum(zmax)
    if drawColz:
        hist.Draw("colz")
    if grayGraphs is not None: 
        for g in grayGraphs: g.Draw("f")
    if dotext:
        rt.gStyle.SetPaintTextFormat('4.%df' % numDigits)
        hist.SetMarkerSize(textSize)
        if not drawErrs: hist.Draw('text'+((drawColz)*('same')))
        else: 
            hist.SetMarkerSize(textSize-1)
            hist.Draw('texte'+((drawColz)*('same')))
    #add LaTeX 
    if drawCMSPreliminary:
        CMS_lumi(c, .15, lumistr=lumistr)
        #t1 = rt.TLatex(0.1,0.94, "CMS preliminary")
        #t2 = rt.TLatex(0.55,0.94, ((lumistr != "")*((lumistr)+' ('))+'13 TeV'+((lumistr != "")*(')')))
        #t1.SetNDC()
        #t2.SetNDC()
        #t1.SetTextSize(0.05)
        #t2.SetTextSize(0.05)
        #t1.Draw()
        #t2.Draw()
    if commentstr != "":
        t3 = rt.TLatex(commentX, commentY, commentstr)
        t3.SetNDC()
        t3.SetTextSize(0.05)
        t3.Draw()
    #save
    if savepng: c.Print(printdir+'/'+printstr+".png")
    if savepdf: c.Print(printdir+'/'+printstr+".pdf")
    if saveroot: c.Print(printdir+'/'+printstr+".root")
    if savec: c.Print(printdir+'/'+printstr+".C")

def make1DRatioHistogram(num, denom, xtitle="", ratiomin=0.25, ratiomax=2.0, logx=False, forPad2=True, ignoreNumeratorErrs=False, ignoreDenominatorErrs=False):
    ratio = num.Clone()
    denomClone = denom.Clone()
    if ignoreNumeratorErrs:
        for bx in range(1, ratio.GetNbinsX()+1):
            ratio.SetBinError(bx, 0.0)
    if ignoreDenominatorErrs:
        for bx in range(1, denomClone.GetNbinsX()+1):
            denomClone.SetBinError(bx, 0.0)
    ratio.SetTitle("")
    ratio.Divide(denomClone)
    if forPad2:
        applyPad2RatioStyle(ratio, xtitle, ratiomin, ratiomax, logx)
    return ratio

def applyPad2RatioStyle(ratio, xtitle, ratiomin, ratiomax, logx):
    ratio.GetXaxis().SetTitle(xtitle)
    ratio.SetMinimum(ratiomin)
    ratio.SetMaximum(ratiomax)
    if hasattr(ratio, "SetStats"): ratio.SetStats(0)
    ratio.GetXaxis().SetLabelSize(0.1)
    ratio.GetYaxis().SetLabelSize(0.08)
    ratio.GetYaxis().SetTitleOffset(0.35)
    ratio.GetXaxis().SetTitleOffset(1.50)
    ratio.GetYaxis().SetTitleSize(0.08)
    ratio.GetXaxis().SetTitleSize(0.08)
    ratio.SetTitle("")

def make1DPullHistogram(h1, h2, xtitle="", ymin=-5.0, ymax=5.0, logx=False, forPad2=True, suppress=True, suppressLevel=0.1):
    """Makes (h1 - h2)/sigma histogram, where sigma is the error on the difference"""
    ret = h1.Clone(h1.GetName()+h2.GetName()+"Pulls")
    ret.Add(h2, -1)
    for bx in range(1, h1.GetNbinsX()+1):
        content = ret.GetBinContent(bx)
        err1 = h1.GetBinError(bx)
        err2 = h2.GetBinError(bx)
        err = (err1*err1 + err2*err2)**(0.5)
        if err > 0:
            ret.SetBinContent(bx,content*1.0/err)
            ret.SetBinError(bx, 0.0)
        else:
            ret.SetBinContent(bx,-9999)
        if suppress:
            if h1.GetBinContent(bx) < suppressLevel and h2.GetBinContent(bx) < suppressLevel:
                ret.SetBinContent(bx,-9999)
    ret.SetMinimum(ymin)
    ret.SetMaximum(ymax)
    if forPad2:
        ret.GetXaxis().SetLabelSize(0.1)
        ret.GetYaxis().SetLabelSize(0.08)
        ret.GetYaxis().SetTitleOffset(0.35)
        ret.GetXaxis().SetTitleOffset(1.50)
        ret.GetYaxis().SetTitleSize(0.08)
        ret.GetXaxis().SetTitleSize(0.08)
        ret.SetTitle("")
        if logx: 
            ret.GetXaxis().SetMoreLogLabels()
            ret.GetXaxis().SetNoExponent()

    return ret

def make1DPercentDiffHistogram(h1, h2, xtitle="", ymin=-1.0, ymax=1.0, logx=False, forPad2=True):
    """Makes (h1 - h2)/h2 histogram"""
    ret = h1.Clone(h1.GetName()+h2.GetName()+"PercDiff")
    ret.Add(h2, -1)
    ret.Divide(h2)
    ret.SetMinimum(ymin)
    ret.SetMaximum(ymax)
    if forPad2:
        ret.GetXaxis().SetLabelSize(0.1)
        ret.GetYaxis().SetLabelSize(0.08)
        ret.GetYaxis().SetTitleOffset(0.35)
        ret.GetXaxis().SetTitleOffset(1.50)
        ret.GetYaxis().SetTitleSize(0.08)
        ret.GetXaxis().SetTitleSize(0.08)
        ret.SetTitle("")
        if logx: 
            ret.GetXaxis().SetMoreLogLabels()
            ret.GetXaxis().SetNoExponent()
    return ret

def make2DPullHistogram(h1, h2, suppress=True):
    """Makes (h1 - h2)/sigma histogram, where sigma is the error on the difference"""
    ret = h1.Clone(h1.GetName()+h2.GetName()+"Pulls")
    ret.Add(h2, -1)
    for bx in range(1, h1.GetNbinsX()+1):
        for by in range(1, h1.GetNbinsY()+1):
            content = ret.GetBinContent(bx,by)
            err1 = h1.GetBinError(bx,by)
            err2 = h2.GetBinError(bx,by)
            err = (err1*err1 + err2*err2)**(0.5)
            if err > 0:
                ret.SetBinContent(bx,by,content*1.0/err)
            else:
                ret.SetBinContent(bx,by,0)
            #suppress bins that have < 1 entry in both histograms
            if suppress and h1.GetBinContent(bx,by) < 1 and h2.GetBinContent(bx,by) < 1:
                ret.SetBinContent(bx,by,-9999)
    return ret

def make2DPercentDiffHistogram(h1, h2, suppress=True):
    """Makes (h1 - h2)/h2 histogram"""
    ret = h1.Clone(h1.GetName()+h2.GetName()+"PercentDiff")
    ret.Add(h2, -1)
    ret.Divide(h2)
    if suppress:
        for bx in range(1, h1.GetNbinsX()+1):
            for by in range(1, h1.GetNbinsY()+1):
                #suppress bins that have < 1 entry in both histograms
                if h1.GetBinContent(bx,by) < 1 and h2.GetBinContent(bx,by) < 1:
                    ret.SetBinContent(bx,by,-9999)
    return ret

def make2DRelativeUncertaintyHistogram(h, suppress=True, suppressLevel=10.0):
    ret = h.Clone(h.GetName()+"RelUnc")
    for bx in range(1, h.GetNbinsX()+1):
        for by in range(1, h.GetNbinsY()+1):
            if h.GetBinContent(bx,by) != 0:
                ret.SetBinContent(bx,by,h.GetBinError(bx,by)*1.0/h.GetBinContent(bx,by))
            else:
                ret.SetBinContent(bx,by,-9999)
            if suppress and ret.GetBinContent(bx,by) > suppressLevel:
                ret.SetBinContent(bx,by,-9999)
    return ret

def unroll2DHistograms(hists, xbins=None, cols=None, labelBins=False, suffix=""):
    """Convert a 2D histogram into a 1D histogram.  
    Bins of the input histogram can be selectively merged. Provide a list (xbins) of bin edges in
    the x-direction and a list of lists (cols) giving the bin edges in the y-direction for each
    bin in x.  
    """
    if xbins is None or cols is None: labelBins = False
    out = [] 
    for hist in hists:
        if hist is None or hist == 0: 
            out.append(None)
            continue
        if xbins is None or cols is None: #unroll using existing 2D binning
            numbins = hist.GetNbinsX()*hist.GetNbinsY()
            outHist = rt.TH1F(hist.GetName()+"Unroll"+suffix, 
                    hist.GetTitle(), numbins, 0, numbins)
            ibin = 0
            for bx in range(1, hist.GetNbinsX()+1):
                for by in range(1, hist.GetNbinsY()+1):
                    ibin += 1
                    outHist.SetBinContent(ibin, hist.GetBinContent(bx,by))
                    outHist.SetBinError(ibin, hist.GetBinError(bx,by))
        else: #make temporary TH2Poly from the binning provided, and unroll that
            poly = macro.makeTH2PolyFromColumns(hist.GetName()+suffix, 
                    hist.GetTitle(), xbins, cols)
            macro.fillTH2PolyFromTH2(hist, poly)
            numbins = poly.GetNumberOfBins()
            outHist = rt.TH1F(hist.GetName()+"Unroll"+suffix, 
                    hist.GetTitle(), numbins, 0, numbins)
            for bn in range(1, numbins+1):
                outHist.SetBinContent(bn, poly.GetBinContent(bn))
                outHist.SetBinError(bn, poly.GetBinError(bn))
        if labelBins:
            labelUnrollBins(outHist, xbins, cols)
        out.append(outHist)
    return out

def plot_basic_2D(c, mc=0, data=0, fit=0, xtitle="", ytitle="", ztitle="Events", zmin=None, zmax=None, printstr="hist", logx=True, logy=True, logz=True, lumistr="", commentstr="", dotext=True, saveroot=True, savepdf=True, savec=True, savepng=True, nsigmaFitData=None, nsigmaFitMC=None, mcDict=None, mcSamples=None, ymin=0.1, unrollBins=(None,None), printdir="."):
    """Plotting macro for data, MC, and/or fit yields.  Creates french flag plots comparing data/MC/fit if able.
    mc, data, fit: 2d histograms for MC prediction, data, and fit prediction (all are optional)
    lumistr: tex-formatted string indicating the integrated luminosity
    commentstr: additional string that will be displayed at the top of the plot 
    nsigmaFitData, nsigmaFitMC (optional): externally provided nsigma histograms to use instead of (yield-prediction)/sigma 
    mcDict (optional): dictionary of MC histograms, for making stacked unrolled plots (provide a list of sample names, mcSamples, to enforce an ordering on the MC histograms in the stack)
    unrollBins (optional): custom binning to use when unrolling 2D histograms into 1D (see unroll2DHistograms function)
    """
    #make a gray square for each -999 bin
    grayGraphs = [makeGrayGraphs(hist) for hist in [mc,fit,data]]
    #unroll 2D hists to plot in 1D
    unrolled = unroll2DHistograms([mc, data, fit], unrollBins[0], unrollBins[1])
    #if individual MC histograms are available, unroll all of them
    if mcDict is not None:
        if mcSamples is not None:
            mcUnrolledList = unroll2DHistograms([mcDict[s] for s in mcSamples], unrollBins[0], unrollBins[1])
            mcUnrolledDict = {mcSamples[n]:mcUnrolledList[n] for n in range(len(mcSamples))}
        else:
            mcSamples = [s for s in mcDict]
            mcUnrolledList = unroll2DHistograms([mcDict[s] for s in mcSamples], unrollBins[0], unrollBins[1])
            mcUnrolledDict = {mcSamples[n]:mcUnrolledList[n] for n in range(len(mcSamples))}
        for s in mcSamples: setHistColor(mcUnrolledDict[s], s)
    if data: unrolled[1].SetBinErrorOption(rt.TH1.kPoisson) #get correct error bars on data
    #synchronize z-axes
    if zmin is None:
        smallestZ = 0.1
        if mc: smallestZ = min(smallestZ, max(0.1,mc.GetMinimum()/1.2))
        if data: smallestZ = min(smallestZ, max(0.1,data.GetMinimum()/1.2))
        if fit: smallestZ = min(smallestZ, max(0.1,fit.GetMinimum()/1.2))
        if mc: mc.SetMinimum(smallestZ)
        if data: data.SetMinimum(smallestZ)
        if fit: fit.SetMinimum(smallestZ)
    if zmax is None:
        largestZ = -float('inf')
        if mc: largestZ = max(largestZ, mc.GetMaximum()*1.2)
        if data: largestZ = max(largestZ, data.GetMaximum()*1.2)
        if fit: largestZ = max(largestZ, fit.GetMaximum()*1.2)
        if mc: mc.SetMaximum(largestZ)
        if data: data.SetMaximum(largestZ)
        if fit: fit.SetMaximum(largestZ)
    #draw each histogram and all relevant combinations
    if mc:
        #draw 2D bin mapping for unrolled plot
        unrolled[0].SetLineColor(rt.kBlack)
        unrolled[0].SetFillColor(38)
        if mcDict is not None: #draw full stack of MC on unrolled plot
            mcStack = makeStack(mcUnrolledDict, mcSamples, "MC")
        else: #draw sum of MC only
            mcStack = makeStack({"MC":unrolled[0]}, ["MC"], "MC")
        draw2DHist(c, mc, xtitle, ytitle, ztitle, zmin, zmax, printstr+'MC', lumistr=lumistr, commentstr=commentstr+", MC prediction", dotext=dotext, grayGraphs=grayGraphs[0], drawErrs=True, saveroot=saveroot, savepdf=savepdf, savepng=savepng, savec=savec, printdir=printdir)
        if data:
            #unroll and compare
            #(first remove any blinded bins)
            if mcDict is not None:
                legDataMC = rt.TLegend(0.75, 0.6, 0.9, 0.9)
                blindMCUnrolledDict = {}
                #blind and create stack
                for s in mcSamples: 
                    blindMCUnrolledDict[s] = mcUnrolledDict[s].Clone()
                    for bx in range(1,blindMCUnrolledDict[s].GetNbinsX()+1):
                        if unrolled[1].GetBinContent(bx) < 0:
                            blindMCUnrolledDict[s].SetBinContent(bx, -999)
                blindStack = makeStack(blindMCUnrolledDict, mcSamples, "MC")
                #add each sample to legend
                for n in range(len(mcSamples)-1,-1,-1):
                    s = mcSamples[n]
                    legDataMC.AddEntry(blindMCUnrolledDict[s], s)
            else:
                legDataMC = rt.TLegend(0.75, 0.7, 0.9, 0.9)
                blindMC = unrolled[0].Clone("blindMC")
                for bx in range(1,blindMC.GetNbinsX()+1):
                    if unrolled[1].GetBinContent(bx) < 0:
                        blindMC.SetBinContent(bx,-999)
                blindStack = makeStack({"MC":blindMC}, ["MC"], "MC")
                legDataMC.AddEntry(blindMC, "MC Prediction")
            legDataMC.AddEntry(unrolled[1], "Data")
            rt.SetOwnership(legDataMC, False)
            plot_basic(c, blindStack, unrolled[1], None, legDataMC, xtitle="Bin", ymin=ymin, printstr=printstr+"UnrolledDataMC", lumistr=lumistr, mcErrColor=rt.kBlack, commentstr=commentstr.replace('Sideband Fit',''), ratiomin=0.0,ratiomax=5, pad2Opt="ratio", saveroot=True, savec=savec, printdir=printdir)
            legDataMC.Delete()
            blindStack.Delete()
        if fit: 
            #do (fit - mc)/unc
            if nsigmaFitMC is None:
                #do (MC - fit)/unc
                mcFitPulls = make2DPullHistogram(fit,mc)
                note="(MC - Fit)/#sigma"
                printnote="MCFitPulls"
            else:
                #make nsigma plot
                note="Nsigmas"
                printnote="MCFitNSigma"
                mcFitPulls = nsigmaFitMC
            draw2DHist(c, mcFitPulls, xtitle, ytitle, ztitle, None, None, printstr+printnote, lumistr=lumistr, commentstr=commentstr+note, palette="FF", logz=False, dotext=dotext, grayGraphs=grayGraphs[0], saveroot=saveroot, savepdf=savepdf, savepng=savepng, savec=savec, printdir=printdir)
            nsigmaUnrolledFitMC = None
            if nsigmaFitMC is not None:
                nsigmaUnrolledFitMC = unroll2DHistograms([nsigmaFitMC], unrollBins[0], unrollBins[1])[0]
                for bx in range(1, nsigmaUnrolledFitMC.GetNbinsX()+1):
                    nsigmaUnrolledFitMC.SetBinError(bx,0.0)

            if mcDict is None:
                legMCFit = rt.TLegend(0.75, 0.7, 0.9, 0.9)
            else:
                legMCFit = rt.TLegend(0.75, 0.6, 0.9, 0.9)
            rt.SetOwnership(legMCFit, False)
            legMCFit.AddEntry(unrolled[2], "Fit Prediction")
            if mcDict is None:
                #add one legend entry for MC
                legMCFit.AddEntry(unrolled[0], "MC Prediction")
            else:
                #add each sample to legend
                for n in range(len(mcSamples)-1,-1,-1): #iterate backwards
                    s = mcSamples[n]
                    legMCFit.AddEntry(mcDict[s], s)
            if nsigmaUnrolledFitMC is not None:
                pad2OptToUse = 'ff'
                ratiominToUse = -5
                ratiomaxToUse = 5
            else:
                pad2OptToUse = 'ratio'
                ratiominToUse= 0.0
                ratiomaxToUse=5.0
            plot_basic(c, mcStack, None, unrolled[2], legMCFit, xtitle="Bin", ymin=ymin, printstr=printstr+"UnrolledMCFit", lumistr=lumistr, commentstr=commentstr, ratiomin=ratiominToUse, ratiomax=ratiomaxToUse, pad2Opt=pad2OptToUse, saveroot=True, savec=savec, mcErrColor=rt.kBlack, customPad2Hist=nsigmaUnrolledFitMC, printdir=printdir)
            mcStack.Delete()
            legMCFit.Delete()
    if data:
        draw2DHist(c, data, xtitle, ytitle, ztitle, zmin=max(0.1,zmin), printstr=printstr+'Data', lumistr=lumistr, commentstr=commentstr+", Data", dotext=dotext, grayGraphs=grayGraphs[2], saveroot=saveroot, savepdf=savepdf, savepng=savepng, savec=savec, printdir=printdir)
        if fit: 
            if nsigmaFitData is None:
                #do (data - fit)/unc
                dataFitPulls = make2DPullHistogram(data,fit)
                note=", (Data - Fit)/#sigma"
                printnote="DataFitPulls"
            else:
                #make nsigma plot
                note=", Nsigmas"
                printnote="DataFitNSigma"
                dataFitPulls = nsigmaFitData
            draw2DHist(c, dataFitPulls, xtitle, ytitle, ztitle, None, None, printstr+printnote, lumistr=lumistr, commentstr=commentstr+note, dotext=dotext, palette="FF", logz=False, grayGraphs=grayGraphs[2], saveroot=saveroot, savepdf=savepdf, savepng=savepng, savec=savec, printdir=printdir)
            #unroll and compare
            blindFit = unrolled[2].Clone("blindFit")

            #unroll and compare
            grayLines = getLinesForUnrolled(fit)

            nsigmaUnrolled = None
            if nsigmaFitData is not None:
                nsigmaUnrolled = unroll2DHistograms([nsigmaFitData], unrollBins[0], unrollBins[1])[0]
                for bx in range(1,nsigmaUnrolled.GetNbinsX()+1):
                    nsigmaUnrolled.SetBinError(bx,0.0)

            for bx in range(1,blindFit.GetNbinsX()+1):
                if unrolled[1].GetBinContent(bx) < 0:
                    blindFit.SetBinContent(bx, -999)

            legDataFit = rt.TLegend(0.7, 0.7, 0.9, 0.9)
            rt.SetOwnership(legDataFit, False)
            legDataFit.AddEntry(blindFit, "Fit Prediction")
            legDataFit.AddEntry(unrolled[1], "Data")

            plot_basic(c, None, unrolled[1], blindFit, legDataFit, xtitle="Bin", ymin=ymin, printstr=printstr+"UnrolledDataFit", lumistr=lumistr, commentstr=commentstr, ratiomin=-5., ratiomax=5.0, pad2Opt="ff", fitColor=rt.kBlue, saveroot=True, savec=savec, customPad2Hist=nsigmaUnrolled, printdir=printdir, grayLines=grayLines)
            legDataFit.Delete()
    if fit:
        draw2DHist(c, fit, xtitle, ytitle, ztitle, zmin, zmax, printstr+'Fit', lumistr=lumistr, commentstr=commentstr+", Fit prediction", grayGraphs=grayGraphs[1], dotext=dotext, drawErrs=True, saveroot=saveroot, savepdf=savepdf, savepng=savepng, savec=savec, printdir=printdir)

def makeStackAndPlot(canvas, mcHists={}, dataHist=None, dataName="Data", mcOrdering=[], titles=[], mcTitle="Stack", xtitle="", ytitle="Events", printstr="hist", logx=False, logy=True, lumistr="40 pb^{-1}", saveroot=True, savepdf=True, savepng=True, savec=True, ymin=None, ymax=None):
    #make stack
    stack = makeStack(mcHists, mcOrdering, mcTitle)
    #make legend
    hists = copy.copy(mcHists)
    hists[dataName] = dataHist
    ordering = copy.copy(mcOrdering)
    ordering.append(dataName)
    leg = makeLegend(hists, titles, ordering)
    #plot
    plot_basic(canvas, stack, dataHist, leg=leg, xtitle=xtitle, ytitle=ytitle, printstr=printstr, logx=logx, logy=logy, lumistr=lumistr, saveroot=saveroot, savepdf=savepdf, savepng=savepng, savec=savec, ymin=ymin)
    stack.Delete()  
    leg.Delete()

def table_basic(headers=[], cols=[], caption="", label="", printstr='table', landscape=False, printdir='.', size='footnotesize'):
    #check for input
    if len(cols) == 0:
        print "table_basic: no columns provided.  doing nothing."
        return
    #check that all columns have the same length
    for col in cols:
        if len(col) != len(cols[0]):
            print "Error in table_basic: columns do not have equal lengths!"
            print "(",[len(c) for c in cols],")"
            return
    #check that there is a header for each column
    if len(headers) != len(cols):
        print "Error in table_basic: number of headers does not equal number of columns!"
        return

    with open(printdir+'/'+printstr+'.tex', 'w') as f:
        f.write('\\newgeometry{margin=0.2cm}\n')
        if landscape: f.write('\\begin{landscape}\n')
        f.write('\\begin{center}\n\\'+size+'\n\\begin{longtable}{|'+('|'.join(['c' for c in cols]))+'|}\n')
        f.write('\\hline\n')
        f.write(' & '.join(headers)+' \\\\\n\\hline\n')
        for row in range(len(cols[0])):
            f.write((' & '.join([col[row] for col in cols]))+' \\\\\n\\hline\n')
        f.write('\\caption{'+caption+'}\n\\label{tab:'+label+'}\n')
        f.write('\\end{longtable}\n\\end{center}\n')
        if landscape: f.write('\\end{landscape}\n')
        f.write('\\restoregeometry\n')
        print "Created LaTeX scale factor table",(printstr+".tex")

def setTDRStyle():
    tdrStyle = rt.TStyle("tdrStyle","Style for P-TDR")

    # For the canvas:
    tdrStyle.SetCanvasBorderMode(0)
    tdrStyle.SetCanvasColor(kWhite)
    tdrStyle.SetCanvasDefH(600) #Height of canvas
    tdrStyle.SetCanvasDefW(600) #Width of canvas
    tdrStyle.SetCanvasDefX(0)   #POsition on screen
    tdrStyle.SetCanvasDefY(0)

    # For the Pad:
    tdrStyle.SetPadBorderMode(0)
    tdrStyle.SetPadColor(kWhite)
    tdrStyle.SetPadGridX(False)
    tdrStyle.SetPadGridY(False)
    tdrStyle.SetGridColor(0)
    tdrStyle.SetGridStyle(3)
    tdrStyle.SetGridWidth(1)
  
    # For the frame:
    tdrStyle.SetFrameBorderMode(0)
    tdrStyle.SetFrameBorderSize(1)
    tdrStyle.SetFrameFillColor(0)
    tdrStyle.SetFrameFillStyle(0)
    tdrStyle.SetFrameLineColor(1)
    tdrStyle.SetFrameLineStyle(1)
    tdrStyle.SetFrameLineWidth(1)
    
    # For the histo:
    tdrStyle.SetHistLineColor(1)
    tdrStyle.SetHistLineStyle(0)
    tdrStyle.SetHistLineWidth(1)
  
    tdrStyle.SetEndErrorSize(2)
    
    tdrStyle.SetMarkerStyle(20)
    
    # For the fit/function:
    tdrStyle.SetOptFit(1)
    tdrStyle.SetFitFormat("5.4g")
    tdrStyle.SetFuncColor(2)
    tdrStyle.SetFuncStyle(1)
    tdrStyle.SetFuncWidth(1)
  
    # For the date:
    tdrStyle.SetOptDate(0)
  
    # For the statistics box:
    tdrStyle.SetOptFile(0)
    tdrStyle.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")
    tdrStyle.SetStatColor(kWhite)
    tdrStyle.SetStatFont(42)
    tdrStyle.SetStatFontSize(0.025)
    tdrStyle.SetStatTextColor(1)
    tdrStyle.SetStatFormat("6.4g")
    tdrStyle.SetStatBorderSize(1)
    tdrStyle.SetStatH(0.1)
    tdrStyle.SetStatW(0.15)
  
    # Margins:
    tdrStyle.SetPadTopMargin(0.10)
    tdrStyle.SetPadBottomMargin(0.13)
    tdrStyle.SetPadLeftMargin(0.10)
    tdrStyle.SetPadRightMargin(0.02)
  
    # For the Global title:
    tdrStyle.SetOptTitle(0)
    tdrStyle.SetTitleFont(42)
    tdrStyle.SetTitleColor(1)
    tdrStyle.SetTitleTextColor(1)
    tdrStyle.SetTitleFillColor(10)
    tdrStyle.SetTitleFontSize(0.05)
  
    # For the axis titles:
    tdrStyle.SetTitleColor(1, "XYZ")
    tdrStyle.SetTitleFont(42, "XYZ")
    tdrStyle.SetTitleSize(0.06, "XYZ")
    tdrStyle.SetTitleXOffset(0.9)
    tdrStyle.SetTitleYOffset(1.25)
  
    # For the axis labels:
    tdrStyle.SetLabelColor(1, "XYZ")
    tdrStyle.SetLabelFont(42, "XYZ")
    tdrStyle.SetLabelOffset(0.007, "XYZ")
    tdrStyle.SetLabelSize(0.05, "XYZ")
  
    # For the axis:
    tdrStyle.SetAxisColor(1, "XYZ")
    tdrStyle.SetStripDecimals(kTRUE)
    tdrStyle.SetTickLength(0.03, "XYZ")
    tdrStyle.SetNdivisions(510, "XYZ")
  
    # Change for log plots:
    tdrStyle.SetOptLogx(0)
    tdrStyle.SetOptLogy(0)
    tdrStyle.SetOptLogz(0)
  
    # Postscript options:
    tdrStyle.SetPaperSize(20.,20.)
  
    tdrStyle.SetHatchesLineWidth(5)
    tdrStyle.SetHatchesSpacing(0.05)
  
def CMS_lumi(pad, relPosX=0.1, iPeriod=4, iPosX=0, lumistr="2.3 fb^{-1}", writeExtraText=True):
    lumi_13TeV = lumistr
    cmsText = "CMS"
    cmsTextFont = 61 # default is helvetica-bold
    extraText = " Preliminary"
    extraTextFont = 52 # default is helvetica-italics
    lumiTextSize = 0.6
    lumiTextOffset = 0.2
    cmsTextSize = 0.75
    relPosY = 0.035
    relExtraDY = 1.2
    extraOverCmsTextSize = 0.76 # ratio of "CMS" and extra text size
    drawLogo = False

    outOfFrame = False
    if iPosX/10 == 0: 
        outOfFrame = True

    alignY_=3
    alignX_=2
    if iPosX/10==0: alignX_=1
    if iPosX==0: alignX_=1
    if iPosX==0: alignY_=1
    if iPosX/10==1: alignX_=1
    if iPosX/10==2: alignX_=2
    if iPosX/10==3: alignX_=3

    align_ = 10*alignX_ + alignY_
    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()

    pad.cd()

    lumiText = ""
    if iPeriod==4:
        lumiText += lumi_13TeV
        lumiText += " (13 TeV)"
    elif iPeriod==0:
        lumiText += lumi_sqrtS
     
    latex = rt.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(rt.kBlack)    

    extraTextSize = extraOverCmsTextSize*cmsTextSize

    latex.SetTextFont(42)
    latex.SetTextAlign(31) 
    latex.SetTextSize(lumiTextSize*t)    
    latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText)

    if outOfFrame:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11) 
        latex.SetTextSize(cmsTextSize*t)    
        latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText)
    
    pad.cd()

    posX_=0
    if iPosX%10<=1:
        posX_ = l + relPosX*(1-l-r)
    elif iPosX%10==2:
        posX_ = l + 0.5*(1-l-r)
    elif iPosX%10==3:
        posX_ =  1-r - relPosX*(1-l-r)
    posY_ = 1-t - relPosY*(1-t-b)
    if not outOfFrame:
        if drawLogo:
            posX_ = l + 0.045*(1-l-r)*W/H
            posY_ = 1-t - 0.045*(1-t-b)
            xl_0 = posX_
            yl_0 = posY_ - 0.15
            xl_1 = posX_ + 0.15*H/W
            yl_1 = posY_
            pad_logo = rt.TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 )
            pad_logo.Draw()
            pad_logo.cd()
            pad_logo.Modified()
            pad.cd()
        else:
            latex.SetTextFont(cmsTextFont)
            latex.SetTextSize(cmsTextSize*t)
            latex.SetTextAlign(align_)
            latex.DrawLatex(posX_, posY_, cmsText)
            if writeExtraText:
                latex.SetTextFont(extraTextFont)
                latex.SetTextAlign(align_)
                latex.SetTextSize(extraTextSize*t)
                latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText)
    elif writeExtraText:
        if iPosX==0:
            posX_ =   l +  relPosX*(1-l-r)
            posY_ =   1-t+lumiTextOffset*t
        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize*t)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, extraText)      

def getUnrollBinsFromHistogram(hist):
    """Get binning from histogram, in a format suitable for the unrollBins function"""
    xbins = []
    cols  = []
    for bx in range(1,hist.GetNbinsX()+1):
        xbins.append(hist.GetXaxis().GetBinLowEdge(bx))
        if bx == hist.GetNbinsX(): xbins.append(hist.GetXaxis().GetBinUpEdge(bx)) #last bin edge
        cols.append([])
        for by in range(1, hist.GetNbinsY()+1):
            cols[-1].append(hist.GetYaxis().GetBinLowEdge(by))
            if by == hist.GetNbinsY(): cols[-1].append(hist.GetYaxis().GetBinUpEdge(by)) #last bin edge
    return (xbins,cols)

def plot_SUS15004(c, data=0, fit=0, printstr="hist", lumistr="", 
        commentstr="", mcDict=None, mcSamples=None, unrollBins=(None,None), 
        printdir=".", ratiomin=0, ratiomax=5, controlRegion=False,
        emptyBinErrs=None):
    if mcDict is None or mcSamples is None:
        print "Error in plot_SUS15004: please provide list of MC samples and associated histograms!"
        return
    ymin=5e-3
    mcTitles = {'WJets':'W+Jets', 'DYJets':'Z #rightarrow ll','TTJets':'t#bar{t}+Jets', 'TTJets2L':'2l t#bar{t}+Jets', 'TTJets1L':'1l t#bar{t}+Jets', 'SingleTop':'Single top', 'QCD':'QCD', 'ZInv':'Z #rightarrow #nu #nu', 'GJets':'#gamma+Jets', 'GJetsInv':'#gamma+Jets', 'GJetsFrag':'#gamma+Jets (frag.)', 'Other':'Other'}
    mcOrdering = ['GJets','GJetsInv','GJetsFrag','WJets','ZInv','TTJets','TTJets1L','TTJets2L','DYJets','SingleTop','QCD','Other']
    #unroll 2D hists to plot in 1D
    if unrollBins[0] is None or unrollBins[1] is None: #if unrolledBins not provided, use histogram binning
        unrollBins = getUnrollBinsFromHistogram(data)
    unrolled = unroll2DHistograms([data, fit], unrollBins[0], unrollBins[1], labelBins=True)
    mcUnrolledList = unroll2DHistograms([mcDict[s] for s in mcSamples], unrollBins[0], unrollBins[1], labelBins=True)
    mcUnrolledDict = {mcSamples[n]:mcUnrolledList[n] for n in range(len(mcSamples))}
    for s in mcSamples: setHistColor(mcUnrolledDict[s], s)
    #draw full stack of MC on unrolled plot
    mcStack = makeStack(mcUnrolledDict, mcSamples, "MC")
    leg = rt.TLegend(0.15, 0.21, 0.45, 0.34)
    rt.SetOwnership(leg, False)
    if data:
        leg.AddEntry(unrolled[0], "Data")
    elif fit: 
        leg.AddEntry(unrolled[1], "Fit Prediction")
    #add each sample to legend
    for s in mcOrdering:
        if s in mcSamples:
            leg.AddEntry(mcUnrolledDict[s], mcTitles[s], "f")
    #add any samples not in the predefined list
    for s in mcSamples:
        if s not in mcOrdering:
            leg.AddEntry(mcUnrolledDict[s], s, "f")
    plot_SUS15004_Unrolled(c, mcStack, unrolled[0], None, leg, ymin=ymin, 
            printstr=printstr+"UnrolledDataMC", lumistr=lumistr, 
            commentstr=commentstr, ratiomin=ratiomin, ratiomax=ratiomax, 
            printdir=printdir, unrollBins=unrollBins, 
            controlRegion=controlRegion, emptyBinErrs=emptyBinErrs)
    mcStack.Delete()
    leg.Delete()

def plot_SUS15004_FitVsMCTotal(c, mcTotal=0, fit=0, printstr="hist", lumistr="", commentstr="", unrollBins=(None,None), printdir=".", ratiomin=0, ratiomax=5):
    ymin=5e-3
    #unroll 2D hists to plot in 1D
    if unrollBins[0] is None or unrollBins[1] is None: #if unrolledBins not provided, use histogram binning
        unrollBins = getUnrollBinsFromHistogram(mcTotal)
    unrolled = unroll2DHistograms([mcTotal, fit], unrollBins[0], unrollBins[1], labelBins=True)
    unrolled[0].SetFillColor(rt.kAzure-9)
    mcStack = makeStack({"MC":unrolled[0]}, ["MC"], "MC")
    leg = rt.TLegend(0.15, 0.21, 0.5, 0.27)
    rt.SetOwnership(leg, False)
    leg.AddEntry(unrolled[0], "Method A Pred.","f")
    leg.AddEntry(unrolled[1], "Method B Pred.")
    plot_SUS15004_Unrolled(c, mcStack, None, unrolled[1], leg, ymin=ymin, printstr=printstr+"UnrolledMCFit", 
            lumistr=lumistr, commentstr=commentstr, ratiomin=ratiomin, ratiomax=ratiomax, printdir=printdir, 
            unrollBins=unrollBins, isPreliminary=True)
    mcStack.Delete()
    leg.Delete()

def plot_SUS15004_MCTotalWithSignal(c, mcTotalUnrolled=0, signalUnrolled=0, 
        printstr="hist", lumistr="", commentstr="", unrollBins=(None,None), 
        printdir=".", ratiomin=0, ratiomax=5, signalString="Signal", 
        dataUnrolled=None):
    ymin=5e-3
    #unroll 2D hists to plot in 1D
    mcTotalUnrolled.SetFillColor(rt.kAzure-9)
    mcStack = makeStack({"MC":mcTotalUnrolled}, ["MC"], "MC")
    leg = rt.TLegend(0.15, 0.21, 0.55, 0.27)
    rt.SetOwnership(leg, False)
    leg.AddEntry(mcTotalUnrolled, "Background Pred.","f")
    leg.AddEntry(signalUnrolled, signalString, "lf")
    plot_SUS15004_Unrolled(c, mcStack, dataUnrolled, None, leg, 
            ymin=ymin, printstr=printstr, lumistr=lumistr, 
            commentstr=commentstr, ratiomin=ratiomin, ratiomax=ratiomax, 
            printdir=printdir, unrollBins=unrollBins, 
            signalHist=signalUnrolled)
    mcStack.Delete()
    leg.Delete()

def labelUnrollBins(hist, xbins=None, cols=None):
    """Label each bin in the unrolled plot with its corresponding y-axis range"""
    bn = 0
    for i in range(len(xbins)-1):
        for j in range(len(cols[i])-1):
            bn += 1
            label = '[%.2f, %.2f]'%(cols[i][j],cols[i][j+1])
            hist.GetXaxis().SetBinLabel(bn, label)

def drawMRLabels(pad, xbins, cols):
    """Draw the labels for the MR bins on unrolled plots"""
    pad.cd()
    latex = rt.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(rt.kBlack)    

    latex.SetTextFont(62)
    latex.SetTextAlign(22) 
    latex.SetTextSize(.035)

    #determine what fraction of the pad is taken up by each MR bin
    plotLeftEdge=0.10
    plotRightEdge=0.90
    binWidths=[]
    for i in range(len(xbins)-1):
        binWidths.append(len(cols[i])-1)
    nbins = sum(binWidths)
    #get left edge of each MR bin
    offsets=[plotLeftEdge]
    for i in range(len(xbins)-1):
        offsets.append(offsets[-1]+binWidths[i]*(plotRightEdge-plotLeftEdge)*1.0/nbins)
    #get center of each MR bin
    binCenters = [(offsets[i]+offsets[i+1])/2.0 for i in range(len(offsets)-1)]
    smallBin = [offsets[i+1]-offsets[i] < 0.09 for i in range(len(offsets)-1)]

    #label each MR bin
    latex.DrawLatex(binCenters[0], 0.85, 'M_{R} [GeV]')
    latex.DrawLatex(plotLeftEdge-0.04, 0.12, 'R^{2}')
    for i in range(len(xbins)-1):
        mrString = "[%d, %d]"%(xbins[i],xbins[i+1])
        if smallBin[i]: 
            latex.SetTextAngle(90)
            latex.DrawLatex(binCenters[i], 0.72, mrString)
        else:
            latex.SetTextAngle(0)
            latex.DrawLatex(binCenters[i], 0.8, mrString)

def getMRLines(pad, xbins, cols, hist, ymin=None, ymax=None):
    """Draws lines at MR borders"""
    if ymin is None:
        ymin = 0
    if ymax is None:
        ymax = hist.GetMaximum()*2
    lines = []
    ibin = 0
    for i in range(0,len(xbins)-2):
        ibin += len(cols[i])-1
        lines.append(rt.TLine(ibin,ymin,ibin,ymax))
        lines[-1].SetLineStyle(2)
        lines[-1].SetLineWidth(1)
        lines[-1].SetLineColor(rt.kBlack)
        lines[-1].Draw()
    return lines

def makeRatioTGraphAsymmErrorsTH1(num, denom, xtitle="", ratiomin=0.25, ratiomax=2.0, logx=False, forPad2=True):
    """Divide a TGraphAsymmErrors by a TH1, ignoring uncertainties on the TH1"""
    npoints = num.GetN()
    ratio = num.Clone()
    xs = num.GetX()
    ys = num.GetY()
    uperrs = num.GetEYhigh()
    downerrs = num.GetEYlow()
    xerrs = num.GetEXlow()
    for bx in range(npoints):
        denomContent = denom.GetBinContent(bx+1)
        if denomContent > 0:
            ratio.SetPoint(bx, xs[bx], ys[bx]/denomContent)
            ratio.SetPointError(bx, xerrs[bx], xerrs[bx], downerrs[bx]/denomContent, uperrs[bx]/denomContent)
        else:
            ratio.SetPoint(bx, ratio.GetXaxis().GetXmin()-1, -999) #move it off the plot
            ratio.SetPointError(bx, 0, 0, 0, 0)
    if forPad2: 
        applyPad2RatioStyle(ratio, xtitle, ratiomin, ratiomax, logx)
    return ratio

def plot_SUS15004_Unrolled(c, mc=0, data=0, fit=0, leg=0, ymin=None, ymax=None, printstr="hist", lumistr="", commentstr="", ratiomin=0.5, ratiomax=1.5, pad2Opt="Ratio", printdir='.', unrollBins=(None,None), controlRegion=False, signalHist=None, isPreliminary=False, emptyBinErrs=None):
    #setup
    c.Clear()
    c.cd()
    pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.28, 1, 1)
    rt.SetOwnership(pad1, False)
    pad1.SetBottomMargin(0.2)
    pad1.SetLogy()
    pad1.Draw()
    pad1.cd()

    ### draw MC

    mc.SetTitle("")

    #make total MC histogram
    histList = mc.GetHists()
    mcTotal = histList.First().Clone()
    mcTotal.SetTitle("")
    mcTotal.SetStats(0)
    mcTotal.Reset()
    numMCHists = 0
    for h in histList:
        mcTotal.Add(h)
        numMCHists += 1
        #for legend
        h.SetMarkerSize(0)
        h.SetLineWidth(0)
        h.SetMarkerColor(rt.kBlack)
        h.SetLineColor(rt.kBlack)
    if emptyBinErrs is not None:
        for ix, err in emptyBinErrs.iteritems():
            oldErr = mcTotal.GetBinError(ix)
            mcTotal.SetBinError(ix,
                    (oldErr*oldErr + err*err)**(0.5))
            print "Error on bin {} increases from {} to {}".format(ix, oldErr, mcTotal.GetBinError(ix))
    mcTotal.SetFillColor(rt.kBlack)
    mc.Draw('hist')
    if ymin is not None: mc.SetMinimum(ymin)
    if ymax is not None: mc.SetMaximum(ymax)
    else: mc.SetMaximum(10*mc.GetMaximum())
    mc.GetXaxis().SetTitle("")
    mc.GetYaxis().SetTitle("Events")
    mc.GetYaxis().SetTitleOffset(0.70)
    mc.GetYaxis().SetTitleSize(0.06)
    mc.GetXaxis().SetLabelSize(0.045)
    mc.GetXaxis().LabelsOption("v")
    mc.GetXaxis().SetLabelFont(62)
    mc.GetYaxis().SetLabelSize(0.06)
    mcTotal.SetFillStyle(3001)
    mcTotal.SetFillColor(rt.kGray+2)
    mcTotal.Draw("e2same")

    #draw fit
    if fit:
        fit.SetStats(0)
        fit.SetMarkerStyle(20)
        fit.SetLineWidth(2)
        fit.SetMarkerSize(0)
        fit.SetLineColor(rt.kBlue)
        fit.SetMarkerColor(rt.kBlue)
        #blue line for fit
        fitCopy = fit.Clone()
        fitCopy.SetFillStyle(0)        
        fitCopy.SetLineWidth(2)
        #Fit = black points
        fit.SetMarkerColor(rt.kBlack)
        fit.SetLineColor(rt.kBlack)
        fit.SetMarkerSize(1)
        fit.Draw('pesame')

    #draw data
    if data:
        #for legend
        data.SetMarkerStyle(20)
        data.SetMarkerSize(1)
        data.SetLineWidth(1)
        data.SetLineColor(rt.kBlack)
        #make TGraphAsymmErrors with appropriate poisson uncertainties
        data = macro.th1ToTGraphAsymmErrors(data)
        data.SetMarkerStyle(20)
        data.SetMarkerSize(1)
        data.SetLineWidth(1)
        data.SetLineColor(rt.kBlack)
        data.Draw("pZ0same")
        #make a copy with no markers to fix error bar issue
        dataNoMarkers = data.Clone()
        dataNoMarkers.SetMarkerSize(0)
        dataNoMarkers.SetMarkerStyle(1)
        dataNoMarkers.Draw("pZ0same")

    #draw signal
    if signalHist:
        signalHist.SetLineWidth(2)
        signalHist.SetLineColor(rt.kAzure-1)
        signalHist.SetFillStyle(3004)
        signalHist.SetFillColor(rt.kAzure-2)
        signalHist.Draw("histsame")

    pad1.Modified()
    rt.gPad.Update()

    #add legend and LaTeX 
    CMS_lumi(pad1, lumistr=lumistr, writeExtraText=isPreliminary)
    if unrollBins[0] is not None and unrollBins[1] is not None:
        drawMRLabels(pad1, unrollBins[0], unrollBins[1])
        lines = getMRLines(pad1, unrollBins[0], unrollBins[1], mcTotal)
        for l in lines: l.Draw()
    leg.SetNColumns(3)
    leg.SetFillColor(rt.kWhite)
    leg.SetTextSize(0.03)
    leg.Draw()

    if commentstr != "":
        commentstr = commentstr.replace("LeptonMultiJet","Lepton Multijet")
        commentstr = commentstr.replace("LeptonJet","Lepton Jet")
        commentstr = commentstr.replace("MultiJet","Multijet")
        commentstr = commentstr.replace("DiJet","Dijet")
        t3 = rt.TLatex(0.40, 0.85, commentstr)
        t3.SetNDC()
        t3.SetTextSize(0.05)
        t3.SetTextFont(42)
        t3.Draw()

    #lower pad plot
    lowerPadHist = None
    lowerPadHist2 = None
    lowerPadHist2Central = None
    lowerPadHist2Lower = None
    lowerPadHist2Upper = None

    #make ratio data/MC
    if data:
        #draw data/MC with poisson errors from data
        lowerPadHist = makeRatioTGraphAsymmErrorsTH1(data, mcTotal, "", ratiomin, ratiomax)
        lowerPadHist.GetYaxis().SetTitle("Data / pred.")
    elif fit:
        lowerPadHist = make1DRatioHistogram(fit, mcTotal, "", ratiomin, ratiomax, ignoreDenominatorErrs=True)
        lowerPadHist.GetYaxis().SetTitle("Method B / Method A")
    #draw relative MC uncertainties 
    lowerPadHist2 = make1DRatioHistogram(mcTotal, mcTotal, "", ratiomin, ratiomax, ignoreDenominatorErrs=True)
    if data or fit:
        lowerPadHist2.GetYaxis().SetTitle(lowerPadHist.GetYaxis().GetTitle())
    else:
        lowerPadHist2.GetYaxis().SetTitle("MC rel. err.")
        lowerPadHist2.SetMaximum(3)
    lowerPadHist2.GetYaxis().SetNdivisions(5, 5, 0)
    lowerPadHist2.GetYaxis().SetLabelSize(0.1)
    lowerPadHist2.GetYaxis().SetTitleSize(0.090)
    lowerPadHist2Upper = lowerPadHist2.Clone()
    lowerPadHist2Upper.SetFillStyle(0)
    lowerPadHist2Upper.SetLineWidth(1)
    lowerPadHist2Lower = lowerPadHist2Upper.Clone()
    lowerPadHist2Central = lowerPadHist2Upper.Clone()
    lowerPadHist2Central.SetLineWidth(2)
    for bx in range(1, lowerPadHist2.GetNbinsX()+1):
        lowerPadHist2Upper.SetBinContent(bx, 1+lowerPadHist2.GetBinError(bx))
        lowerPadHist2Lower.SetBinContent(bx, 1-lowerPadHist2.GetBinError(bx))
        lowerPadHist2.GetXaxis().SetBinLabel(bx, "")

    #draw lower pad
    c.cd()
    pad2 = rt.TPad(printstr+"pad2",printstr+"pad2",0,0.0,1,0.3)
    pad2.SetTopMargin(0.035)
    pad2.SetBottomMargin(0.10)
    pad2.SetGridy()
    pad2.Draw()
    pad2.cd()

    #lower pad hist 2
    lowerPadHist2.Draw("e2same")
    pad2.Modified()
    rt.gPad.Update()
    pad2.Draw()
    pad2.cd()

    #lower pad hist 1
    if data or fit:
        lowerPadHist.Draw("pZ0same")
        pad2.Modified()
        rt.gPad.Update()

    #lower pad hist central/upper/lower
    lowerPadHist2Upper.Draw("histsame")
    lowerPadHist2Lower.Draw("histsame")
    lowerPadHist2.Draw('e2same')
    lowerPadHist2Central.Draw("histsame")
    if data or fit:
        lowerPadHist.Draw("pZ0same")

    #save
    c.Print(printdir+'/'+printstr+".png")
    c.Print(printdir+'/'+printstr+".pdf")
    c.Print(printdir+'/'+printstr+".root")
    c.Print(printdir+'/'+printstr+".C")

    pad1.Delete()

def plot_SUS15004_1D(c, mc=0, data=0, leg=0, xtitle="", ytitle="Events", ymin=None, ymax=None, printstr="hist", logx=False, lumistr="40 pb^{-1}", commentstr="", ratiomin=0.5, ratiomax=1.5, printdir='.'):
    #setup
    c.Clear()
    c.cd()
    pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.4, 1, 1)
    rt.SetOwnership(pad1, False)
    pad1.SetBottomMargin(0.0130)
    pad1.SetLogx(logx)
    pad1.SetLogy()
    pad1.Draw()
    pad1.cd()

    ### draw MC
    mc.SetTitle("")

    #make total MC histogram
    histList = mc.GetHists()
    mcTotal = histList.First().Clone()
    mcTotal.SetStats(0)
    mcTotal.Reset()
    numMCHists = 0
    for h in histList:
        mcTotal.Add(h)
        numMCHists += 1
    mc.Draw('hist')
    if ymin is not None: mc.SetMinimum(ymin)
    if ymax is not None: mc.SetMaximum(ymax)
    mc.GetXaxis().SetMoreLogLabels()
    mc.GetXaxis().SetNoExponent()
    mc.GetXaxis().SetTitle(xtitle)
    mc.GetXaxis().SetTitleSize(0.06)
    mc.GetYaxis().SetTitle(ytitle)
    mc.GetYaxis().SetTitleOffset(0.60)
    mc.GetYaxis().SetTitleSize(0.06)
    mc.GetYaxis().SetLabelSize(0.06)
    mcTotal.SetFillStyle(3001)
    mcTotal.SetFillColor(rt.kGray+2)
    mcTotal.GetXaxis().SetMoreLogLabels()
    mcTotal.GetXaxis().SetNoExponent()
    mcTotal.Draw("e2same")

    #draw data
    if data:
        data.SetStats(0)
        data.SetMarkerStyle(20)
        data.SetMarkerSize(1)
        data.SetLineWidth(1)
        data.SetLineColor(rt.kBlack)
        data.SetTitle("")
        data.Draw("pesame")
    pad1.Modified()
    rt.gPad.Update()

    #add legend and LaTeX 
    leg.SetX1(0.7)
    leg.SetY1(0.4)
    leg.SetX2(0.89)
    leg.SetY2(0.85)
    leg.SetBorderSize(0)
    leg.SetFillColor(rt.kWhite)
    leg.Draw()
    CMS_lumi(pad1, lumistr=lumistr, writeExtraText=False)

    if commentstr != "":
        t3 = rt.TLatex(0.40, 0.92, commentstr)
        t3.SetNDC()
        t3.SetTextSize(0.05)
        t3.Draw()

    #lower pad plot
    lowerPadHist = None
    lowerPadHist2 = None
    lowerPadHist2Central = None
    lowerPadHist2Lower = None
    lowerPadHist2Upper = None

    #make ratio data/MC
    #draw data/MC with poisson errors from data
    if data:
        lowerPadHist = make1DRatioHistogram(data, mcTotal, xtitle, ratiomin, ratiomax, logx, ignoreDenominatorErrs=True)
        lowerPadHist.GetYaxis().SetTitle("Data / pred.")
    #draw relative MC uncertainties
    lowerPadHist2 = make1DRatioHistogram(mcTotal, mcTotal, xtitle, ratiomin, ratiomax, logx, ignoreNumeratorErrs=True)
    if data:
        lowerPadHist2.GetYaxis().SetTitle(lowerPadHist.GetYaxis().GetTitle())
    lowerPadHist2.GetXaxis().SetTitleSize(0.1)
    lowerPadHist2.GetYaxis().SetTitleSize(0.08)
    lowerPadHist2.GetXaxis().SetTitleOffset(1.2)
    lowerPadHist2.GetYaxis().SetTitleOffset(.4)
    lowerPadHist2.GetYaxis().SetNdivisions(5, 5, 0)
    lowerPadHist2.GetYaxis().SetLabelSize(0.08)
    lowerPadHist2.SetFillStyle(3001)
    lowerPadHist2.SetFillColor(rt.kGray+2)

    #draw the pad
    c.cd()
    pad2 = rt.TPad(printstr+"pad2",printstr+"pad2",0,0.0,1,0.4)
    pad2.SetTopMargin(0.1)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.25)
    pad2.SetGridy()
    pad2.SetLogx(logx)
    pad2.Draw()
    pad2.cd()

    lowerPadHist2.Draw("e2same")
    if data:
        lowerPadHist.Draw("pesame")
    pad2.Modified()
    rt.gPad.Update()

    if lowerPadHist2Central is not None:
        lowerPadHist2Upper.Draw("histsame")
        lowerPadHist2Lower.Draw("histsame")
        lowerPadHist2.Draw('e2same')
        lowerPadHist2Central.Draw("histsame")
        lowerPadHist.Draw("pesame")

    #save
    c.Print(printdir+'/'+printstr+".png")
    c.Print(printdir+'/'+printstr+".pdf")
    c.Print(printdir+'/'+printstr+".root")
    c.Print(printdir+'/'+printstr+".C")
    pad1.Delete()

def plotEvidenceHist(c, ev, printstr="hist", printdir='.', lumistr="2.3 fb^{-1}", unrollBins=(None,None), zmin=1e-3, obs=None):
    #setup
    c.Clear()
    c.cd()
    pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.28, 1, 1)
    rt.SetOwnership(pad1, False)
    pad1.SetBottomMargin(0.2)
    pad1.SetLogy(False)
    pad1.Draw()
    pad1.cd()

    ### draw evidence plot
    ev.SetTitle("")
    ev.SetStats(0)
    ev.GetXaxis().SetTitle("")
    ev.GetYaxis().SetTitle("Contribution to -2DLL")
    ev.GetYaxis().SetTitleOffset(1.00)
    ev.GetYaxis().SetTitleSize(0.04)
    ev.GetXaxis().SetLabelSize(0.045)
    ev.GetXaxis().LabelsOption("v")
    ev.GetXaxis().SetLabelFont(62)
    ev.GetYaxis().SetLabelSize(0.04)
    ev.SetMarkerStyle(20)
    ev.SetMarkerColor(9)
    ev.SetMarkerSize(1)
    ev.SetLineWidth(1)
    ev.SetLineColor(rt.kBlack)
    #ev.SetMinimum(zmin)
    ymin = -5
    ymax = 5
    ev.SetMaximum(ymax)
    ev.SetMinimum(ymin)
    ev.Draw("pe")
    if obs is not None:
        obs.SetMarkerStyle(20)
        obs.SetMarkerSize(1)
        obs.Draw('pesame')
    pad1.Modified()
    rt.gPad.Update()

    #add LaTeX 
    CMS_lumi(pad1, lumistr=lumistr)
    if unrollBins[0] is not None and unrollBins[1] is not None:
        drawMRLabels(pad1, unrollBins[0], unrollBins[1])
        lines = getMRLines(pad1, unrollBins[0], unrollBins[1], ev,
                ymin=ymin, ymax=ymax)
        for l in lines: l.Draw()

    #save
    c.Print(printdir+'/'+printstr+".png")
    c.Print(printdir+'/'+printstr+".pdf")
    c.Print(printdir+'/'+printstr+".root")
    c.Print(printdir+'/'+printstr+".C")

    pad1.Delete()
