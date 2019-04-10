import ROOT as rt
import sys
import rootTools
import glob
from math import *
import os
from array import *
from scipy.interpolate import Rbf
import numpy as np
from Get2DContour import getModelSettings,interpolate2D


def plotT1Limits(boxModel, model, sigHist, obsGraphModel, expGraphModel, doHybridNew):

    mgMin, mgMax, mchiMin, mchiMax, binWidth, nRebins, xsecMin, xsecMax, diagonalOffset, smoothing, fixLSP0 = getModelSettings(model)

    #mgMin = mgMin+binWidth/2.
    #mgMax = mgMax+binWidth/2.
    #mchiMin = mchiMin+binWidth/2.
    #mchiMax = mchiMax+binWidth/2.

    mgMin = 600
    mgMax = 1950
    mchiMin = 0
    mchiMax = 1800

    print mgMin, mgMax
    emptyHist = rt.TH2D("emptyHist", "", int((mgMax-mgMin)/25.), mgMin, mgMax, int((mchiMax-mchiMin)/25.), mchiMin, mchiMax)

    modelLabel = "pp #rightarrow #tilde{g}#tilde{g}"
    moreModelLabel = {"T1bbbb":"100% #tilde{g} #rightarrow b#bar{b}#tilde{#chi}^{0}_{1}",
                      "T1tttt":"100% #tilde{g} #rightarrow t#bar{t}#tilde{#chi}^{0}_{1}",
                      "T1x0p50y0p00":"50% #tilde{g} #rightarrow tb#tilde{#chi}^{#pm}_{1}, 50% #tilde{g} #rightarrow bb#tilde{#chi}^{0}_{1}",
                      "T1x0p50y0p50":"50% #tilde{g} #rightarrow tt#tilde{#chi}^{0}_{1}, 50% #tilde{g} #rightarrow bb#tilde{#chi}^{0}_{1}",
                      "T1x0p50y0p25":"25% #tilde{g} #rightarrow tb#tilde{#chi}^{#pm}_{1}, 25% #tilde{g} #rightarrow tt#tilde{#chi}^{0}_{1}, 50% #tilde{g} #rightarrow bb#tilde{#chi}^{0}_{1}",
                      "T1x0p25y0p50":"25% #tilde{g} #rightarrow tb#tilde{#chi}^{#pm}_{1}, 50% #tilde{g} #rightarrow tt#tilde{#chi}^{0}_{1}, 25% #tilde{g} #rightarrow bb#tilde{#chi}^{0}_{1}",
                      "T1x0p00y0p00":"100% #tilde{g} #rightarrow tb#tilde{#chi}^{#pm}_{1}",
                      "T1x0p25y0p25":"50% #tilde{g} #rightarrow tb#tilde{#chi}^{#pm}_{1}, 25% #tilde{g} #rightarrow tt#tilde{#chi}^{0}_{1}, 25% #tilde{g} #rightarrow bb#tilde{#chi}^{0}_{1}",
                      "T1x0p00y0p50":"50% #tilde{g} #rightarrow tb#tilde{#chi}^{#pm}_{1}, 50% #tilde{g} #rightarrow tt#tilde{#chi}^{0}_{1}",
                      "T1tttt":"100% #tilde{g} #rightarrow t#bar{t}#tilde{#chi}^{0}_{1}",
                      "T1qqqq":"100% #tilde{g} #rightarrow q#bar{q}#tilde{#chi}^{0}_{1}",
                      "T1bri":""}
    sparticleLabel = "m_{#tilde{g}} [GeV]"
    massLabel = "m_{#tilde{#chi}^{#pm}_{1}}-m_{#tilde{#chi}^{0}_{1}} = 5 GeV"
    #massLabel = ''
    asymptoticLabel = "Asymptotic"
    if doHybridNew: asymptoticLabel = ""
        
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetPadTopMargin(0.08)
    rt.gStyle.SetPadRightMargin(0.19)
    rt.gStyle.SetPadBottomMargin(0.14)
    rt.gStyle.SetPadLeftMargin(0.165)



    
        #from PlotsSMS/python/smsPlotXSEC.py
    stops = [0.00, 0.20, 0.70, 0.90, 1.00]
    
    red   = [0.00, 0.00, 0.60, 1.00, 1.00, 0.00]
    green = [0.00, 0.50, 0.90, 0.60, 0.00, 1.00]
    blue  = [1.00, 0.60, 0.18, 0.00, 0.00, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)
    
    npoints = len(s)
    ncontours = 5
    rt.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    rt.gStyle.SetNumberContours(ncontours)


    rt.gStyle.SetPaintTextFormat("2.1f")
    
    c = rt.TCanvas("c","c",600,600)
    
    emptyHist.GetXaxis().SetRangeUser(mgMin, mgMax)
    emptyHist.GetYaxis().SetRangeUser(mchiMin, mchiMax)
    
    emptyHist.GetZaxis().SetLabelFont(42)
    emptyHist.GetZaxis().SetTitleFont(42)
    emptyHist.GetZaxis().SetLabelSize(0.035)
    emptyHist.GetZaxis().SetTitleSize(0.035)


    emptyHist.GetXaxis().SetLabelFont(42)
    emptyHist.GetXaxis().SetLabelSize(0.04)
    emptyHist.GetXaxis().SetNdivisions(407,True)
    emptyHist.GetXaxis().SetTitleFont(42)
    emptyHist.GetXaxis().SetTitleSize(0.05)
    emptyHist.GetXaxis().SetTitleOffset(1.2)
    emptyHist.GetXaxis().SetTitle(sparticleLabel)

    # set y axis
    emptyHist.GetYaxis().SetLabelFont(42)
    emptyHist.GetYaxis().SetLabelSize(0.04)
    #emptyHist.GetYaxis().SetNdivisions(407,True)
    emptyHist.GetYaxis().SetTitleFont(42)
    emptyHist.GetYaxis().SetTitleSize(0.05)
    emptyHist.GetYaxis().SetTitleOffset(1.6)
    emptyHist.GetYaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} [GeV]")
    
    emptyHist.Draw()
    

    rt.gPad.Update()
    

    #textCMS = rt.TLatex(0.15,0.98,"CMS %s                 %.1f fb^{-1} (%s TeV)" %('preliminary', 2.300, 13))
    textCMS = rt.TLatex(0.15,0.98,"CMS")
    textCMS.SetNDC()
    textCMS.SetTextAlign(13)
    textCMS.SetTextFont(62)
    textCMS.SetTextSize(0.05)
    textCMS.Draw()
    textCMS1 = rt.TLatex(0.57,0.98,"%.1f fb^{-1} (%s TeV)" %(2.3, 13))
    textCMS1.SetNDC()
    textCMS1.SetTextAlign(13)
    textCMS1.SetTextFont(42)
    textCMS1.SetTextSize(0.038)
    textCMS1.Draw()

    
    textCOLZ = rt.TLatex(0.98,0.45,"95% C.L. exclusion contour")
    textCOLZ.SetNDC()
    #textCOLZ.SetTextAlign(13)
    textCOLZ.SetTextFont(42)
    textCOLZ.SetTextSize(0.045)
    textCOLZ.SetTextAngle(90)
    #textCOLZ.Draw()

    
    xRange = mgMax - mgMin
    yRange = mchiMax - mchiMin

    
    
    # MORE MODEL LABELS
    obsColor = {"T1bbbb":rt.TColor.GetColor(red[0],green[0],blue[0]),
                "T1x0p50y0p50":rt.TColor.GetColor(red[1],green[1],blue[1]),
                "T1x0p50y0p00":rt.TColor.GetColor(red[1],green[1],blue[1]),
                "T1x0p25y0p50":rt.TColor.GetColor(red[1],green[1],blue[1]),
                "T1x0p50y0p25":rt.TColor.GetColor(red[1],green[1],blue[1]),
                "T1x0p00y0p00":rt.TColor.GetColor(red[2],green[2],blue[2]),
                "T1x0p25y0p25":rt.kGray+1,
                "T1x0p00y0p50":rt.TColor.GetColor(red[3],green[3],blue[3]),
                "T1tttt":rt.TColor.GetColor(red[4],green[4],blue[4]),
                "T1qqqq":rt.TColor.GetColor(red[5],green[5],blue[5]),
                "T1bri":rt.kBlack
                }
        
    #red   = [0.00, 0.00, 0.60, 1.00, 1.00, 0.00]
    #green = [0.00, 0.50, 0.90, 0.60, 0.00, 1.00]
    #blue  = [1.00, 0.60, 0.18, 0.00, 0.00, 0.00]
    
    expColor = obsColor
        
    obsDist = {"T1bbbb":1.0-0.1,
               "T1x0p50y0p00":1.6-0.1,
               "T1x0p25y0p50":1.6-0.1,
               "T1x0p50y0p25":1.6-0.1,
               "T1x0p50y0p50":1.6-0.1,
               "T1x0p00y0p00":2.2-0.1,
               "T1x0p25y0p25":2.8-0.1,
               "T1x0p00y0p50":3.4-0.1,
               "T1tttt":4.0-0.1,
               "T1bri":0}
        
    expDist = {"T1bbbb":1.2-0.1,
               "T1x0p50y0p00":1.8-0.1,
               "T1x0p50y0p25":1.8-0.1,
               "T1x0p50y0p50":1.8-0.1,
               "T1x0p25y0p00":1.8-0.1,
               "T1x0p00y0p00":2.4-0.1,
               "T1x0p25y0p25":3.0-0.1,
               "T1x0p00y0p50":3.6-0.1,
               "T1tttt":4.2-0.1,
               "T1bri":0}
        
    diagX = {"T1bbbb":array('d',[0,20000,mgMin]),
             "T1x0p50y0p00":array('d',[0,20000,mgMin]),
             "T1x0p50y0p25":array('d',[0,20000,mgMin]),
             "T1x0p25y0p50":array('d',[0,20000,mgMin]),
             "T1x0p50y0p50":array('d',[0,20000,mgMin]),
             "T1x0p00y0p00":array('d',[0,20000,mgMin]),
             "T1x0p25y0p25":array('d',[0,20000,mgMin]),
             "T1x0p00y0p50":array('d',[0,20000,mgMin]),
             "T1tttt":array('d',[0,20000,mgMin]),
             "T1bri":array('d',[0,20000,mgMin]),
             "T1qqqq":array('d',[0,20000,mgMin])}
    
    diagY = {"T1bbbb":array('d',[-25, 20000-25,mgMax]),
             "T1x0p50y0p00":array('d',[-225, 20000-225,mgMax]),
             "T1x0p50y0p50":array('d',[-225, 20000-225,mgMax]),
             "T1x0p50y0p25":array('d',[-225, 20000-225,mgMax]),
             "T1x0p25y0p50":array('d',[-225, 20000-225,mgMax]),
             "T1x0p00y0p00":array('d',[-225, 20000-225,mgMax]),     
             "T1x0p25y0p25":array('d',[-225, 20000-225,mgMax]),             
             "T1x0p00y0p50":array('d',[-225, 20000-225,mgMax]),
             "T1tttt":array('d',[-225, 20000-225,mgMax]),
             "T1bri":array('d',[-225, 20000-225,mgMax]),
             "T1qqqq":array('d',[-25, 20000-25,mgMax])}

    LObs = {}
    LExp = {}
    textObs = {}
    diagonal = {}
    for model in models:
        diagonal[model] = rt.TGraph(3, diagX[model], diagY[model])
        diagonal[model].SetName("diagonal%s"%model)
        #diagonal[model].SetFillColor(rt.kWhite)
        diagonal[model].SetLineColor(rt.kGray)
        diagonal[model].SetLineStyle(2)
        #diagonal[model].Draw("FSAME")
        diagonal[model].Draw("LSAME")
        
    for model in models:
        expGraphModel[model].SetLineColor(expColor[model])
        expGraphModel[model].SetLineStyle(2)
        expGraphModel[model].Draw("csame")
        obsGraphModel[model].SetLineColor(obsColor[model])
        obsGraphModel[model].Draw("csame")

        
    c.RedrawAxis()
    # white background
    graphWhite = rt.TGraph(5)
    graphWhite.SetName("white")
    graphWhite.SetTitle("white")
    graphWhite.SetFillColor(rt.kWhite)
    graphWhite.SetFillStyle(1001)
    graphWhite.SetLineColor(rt.kBlack)
    graphWhite.SetLineStyle(1)
    graphWhite.SetLineWidth(3)
    graphWhite.SetPoint(0,mgMin, mchiMax)
    graphWhite.SetPoint(1,mgMax, mchiMax)
    graphWhite.SetPoint(2,mgMax, mchiMax*0.565)
    graphWhite.SetPoint(3,mgMin, mchiMax*0.565)
    graphWhite.SetPoint(4,mgMin, mchiMax)
    graphWhite.Draw("FSAME")
    graphWhite.Draw("LSAME")

    # MODEL LABEL
    textModelLabel= rt.TLatex(0.185,0.91,"%s        95%% C.L. NLO+NLL exclusion" %modelLabel)
    textModelLabel.SetNDC()
    textModelLabel.SetTextAlign(13)
    textModelLabel.SetTextFont(42)
    textModelLabel.SetTextSize(0.036)
    textModelLabel.Draw()
    
    # MASS LABEL
    textMassLabel= rt.TLatex(0.6,0.86,"%s"%massLabel)
    textMassLabel.SetNDC()
    textMassLabel.SetTextAlign(13)
    textMassLabel.SetTextFont(42)
    textMassLabel.SetTextSize(0.032)
    textMassLabel.Draw()
        
    for model in models:
        LObs[model] = rt.TGraph(2)
        LObs[model].SetName("LObs%s"%model)
        LObs[model].SetTitle("LObs%s"%model)
        LObs[model].SetLineColor(obsColor[model])
        LObs[model].SetLineStyle(1)
        LObs[model].SetLineWidth(3)
        LObs[model].SetMarkerStyle(20)
        LObs[model].SetPoint(0,mgMin+3*xRange/100, mchiMax-obsDist[model]*yRange/100*10)
        LObs[model].SetPoint(1,mgMin+10*xRange/100, mchiMax-obsDist[model]*yRange/100*10)        
        LExp[model] = rt.TGraph(2)
        LExp[model].SetName("LExp%s"%model)
        LExp[model].SetTitle("LExp%s"%model)
        LExp[model].SetLineColor(expColor[model])
        LExp[model].SetLineStyle(7)
        LExp[model].SetLineWidth(3)
        LExp[model].SetPoint(0,mgMin+3*xRange/100, mchiMax-expDist[model]*yRange/100*10)
        LExp[model].SetPoint(1,mgMin+10*xRange/100, mchiMax-expDist[model]*yRange/100*10)
        textObs[model] = rt.TLatex(mgMin+11*xRange/100, mchiMax-expDist[model]*yRange/100*10, moreModelLabel[model])
        textObs[model].SetTextFont(42)
        textObs[model].SetTextSize(0.032)
        LObs[model].Draw("lsame")
        LExp[model].Draw("lsame")
        textObs[model].Draw()



    model = "test"
    textExp = {}
    LObs[model] = rt.TGraph(2)
    LObs[model].SetName("LObs%s"%model)
    LObs[model].SetTitle("LObs%s"%model)
    LObs[model].SetLineColor(rt.kBlack)
    LObs[model].SetLineStyle(1)
    LObs[model].SetLineWidth(3)
    LObs[model].SetMarkerStyle(20)
    LObs[model].SetPoint(0,mgMin+72*xRange/100-20, mchiMax-3.5*yRange/100*10-20)
    LObs[model].SetPoint(1,mgMin+79*xRange/100-20, mchiMax-3.5*yRange/100*10-20)        
    LExp[model] = rt.TGraph(2)
    LExp[model].SetName("LExp%s"%model)
    LExp[model].SetTitle("LExp%s"%model)
    LExp[model].SetLineColor(rt.kBlack)
    LExp[model].SetLineStyle(7)
    LExp[model].SetLineWidth(3)
    LExp[model].SetPoint(0,mgMin+72*xRange/100-20, mchiMax-3.9*yRange/100*10-20)
    LExp[model].SetPoint(1,mgMin+79*xRange/100-20, mchiMax-3.9*yRange/100*10-20)
    textObs[model] = rt.TLatex(mgMin+81*xRange/100-20, mchiMax-3.6*yRange/100*10-20, "Observed")
    textObs[model].SetTextFont(42)
    textObs[model].SetTextSize(0.032)
    textExp[model] = rt.TLatex(mgMin+81*xRange/100-20, mchiMax-4.0*yRange/100*10-20, "Expected")
    textExp[model].SetTextFont(42)
    textExp[model].SetTextSize(0.032)
    LObs[model].Draw("lsame")
    LExp[model].Draw("lsame")
    textObs[model].Draw()
    textExp[model].Draw()
    
    if doHybridNew:
        c.Print("%s%s%s%s.pdf"%("T1","HybridNew","","BARE"))
        c.Print("%s%s%s%s.C"%("T1","HybridNew","","BARE"))
    else:
        c.Print("%s%s%s%s.pdf"%("T1","Asymptotic","","BARE"))
        c.Print("%s%s%s%s.C"%("T1","Asymptotic","","BARE"))

if __name__ == '__main__':
    directory = sys.argv[1]
    #rt.gROOT.SetBatch()

    models = ["T1bbbb","T1x0p50y0p00","T1x0p00y0p00","T1x0p25y0p25","T1x0p00y0p50","T1tttt"]

    boxModel = {"T1bbbb":"MultiJet",
                "T1tttt":"MuMultiJet_EleMultiJet_MultiJet",
                "T1x0p50y0p00":"MuMultiJet_EleMultiJet_MultiJet",
                "T1x0p50y0p50":"MuMultiJet_EleMultiJet_MultiJet",
                "T1x0p50y0p25":"MuMultiJet_EleMultiJet_MultiJet",
                "T1x0p25y0p50":"MuMultiJet_EleMultiJet_MultiJet",
                "T1x0p00y0p50":"MuMultiJet_EleMultiJet_MultiJet",
                "T1x0p25y0p25":"MuMultiJet_EleMultiJet_MultiJet",
                "T1x0p00y0p00":"MuMultiJet_EleMultiJet_MultiJet",
                "T1bri": "MuMultiJet_EleMultiJet_MultiJet",
                "T1qqqq":"MultiJet"}
    obsGraphModel = {}
    expGraphModel = {}
    for model in models:
        resultsFile = rt.TFile.Open("%s/%s/%s_%s_results.root"%(directory,model,model,boxModel[model]))
        sigHist = resultsFile.Get("xsecUL_Obs_%s_%s"%(model,boxModel[model]))
        obsGraphModel[model] = resultsFile.Get("Obs_%s_%s"%(model,boxModel[model]))
        expGraphModel[model] = resultsFile.Get("Exp_%s_%s"%(model,boxModel[model]))
                              
    doHybridNew = False
    plotT1Limits(boxModel, "T1bbbb", sigHist, obsGraphModel, expGraphModel, doHybridNew)
