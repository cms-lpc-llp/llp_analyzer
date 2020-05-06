from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
from WriteDataCard import *
import os
import random
import sys
import math
from scipy.integrate import quad
import csv

def setStyle():
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptFit(0000)
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetPaintTextFormat("1.2g")
    rt.gROOT.SetBatch()
    rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.FATAL)
    
    rt.gStyle.SetStatY(1.9)
    rt.gStyle.SetStatX(1.9)
    rt.gStyle.SetStatW(0.1)
    rt.gStyle.SetStatH(0.1)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store results")
    parser.add_option('--numerator',dest="numerator",default="HLT_RsqMR240_Rsq0p09_MR200",type="string", 
                  help="numerator trigger")
    parser.add_option('--denominator',dest="denominator",default="HLT_Ele27_eta2p1_WPLoose_Gsf",type="string",
                  help="denominator trigger")
    parser.add_option('--path-names',dest="pathNames",default="data/RazorHLTPathnames.dat",type="string",
                  help="text file containing mapping between array index and path name")
    parser.add_option('--tree-name',dest="treeName",default="RazorInclusive",type="string",
                  help="tree name to use")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    
    
    (options,args) = parser.parse_args()
     

    x = array('d', [0,25,50,100,125,150,175,200,225,250,275,300,325,350,400,500,600,700,900,1200])# MR binning
    y = array('d', [0,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.62,0.74,0.86,1,1.5]) # Rsq binning


    xCuts = [0,300,400,500]
    yCuts = [0.15,0.25]

    setStyle()

    f = open(options.pathNames)
    csvin = csv.reader(f,delimiter=' ')
    triggerDict = {}
    for row in csvin:
        triggerDict[row[-1]] = int(row[0])
    
    
    lumi_in = 0.
    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            tree = rootFile.Get(options.treeName)
                
                
    dPhiCut = 999999
    #dPhiCut = 2.8
    mRmin = 0
    mRmax = 999999
    rsqMin = 0
    rsqMax = 999999
    
    listDenom = ['HLTDecision[%i]'%triggerDict[denom] for denom in options.denominator.split(',')]
    stringDenom = '||'.join(listDenom)

    metFlags = ['Flag_HBHENoiseFilter','Flag_goodVertices','Flag_eeBadScFilter']
    metString = '&&'.join(metFlags)
    
    boxString = 'box == 6 || box == 7 || box ==  8'
    
    muonString = 'leadingMuonPt < 100 || allMuonPt < 100'    
    #muonString = 'nVetoMuons == 0'
    #muonString = '1'
        
    tree.Draw('>>elist',
              'MR > %f && MR < %f && Rsq > %f && Rsq < %f && abs(dPhiRazor) < %f && (%s) && (%s) && (%s) && (%s)' % (mRmin,mRmax,rsqMin,rsqMax,dPhiCut,stringDenom,metString,boxString,muonString),
              'entrylist')
    
        
    elist = rt.gDirectory.Get('elist')

    pNum2D =  [rt.TH2D("num2D","R2-MR numerator;M_{R} [GeV];R^{2};numerator",len(x)-1,x,len(y)-1,y)]
    pDenom2D =  [rt.TH2D("denom2D","R2-MR denominator;M_{R} [GeV];R^{2};denominator",len(x)-1,x,len(y)-1,y)]
    pEff2D =  [rt.TEfficiency("eff2D","R2-MR efficiency;M_{R} [GeV];R^{2};efficiency",len(x)-1,x,len(y)-1,y)]
    [pEff.SetStatisticOption(rt.TEfficiency.kFCP) for pEff in pEff2D]

    
    pNumMR =  [rt.TH1D("numMR_Rsq%.2f"%yCut,"MR numerator;M_{R} [GeV];numerator",len(x)-1,x) for yCut in yCuts]
    pDenomMR =  [rt.TH1D("denomMR_Rsq%.2f"%yCut,"MR denominator;M_{R} [GeV];denominator",len(x)-1,x) for yCut in yCuts]
    pEffMR =  [rt.TEfficiency("effMR_Rsq%.2f"%yCut,"MR efficiency;M_{R} [GeV];efficiency",len(x)-1,x) for yCut in yCuts]
    [pEff.SetStatisticOption(rt.TEfficiency.kFCP) for pEff in pEffMR]
    
    pNumRsq =  [rt.TH1D("numRsq_MR%i"%xCut,"Rsq numerator;R^{2};numerator",len(y)-1,y) for xCut in xCuts]
    pDenomRsq =  [rt.TH1D("denomRsq_MR%i"%xCut,"Rsq denominator;R^{2};denominator",len(y)-1,y) for xCut in xCuts]
    pEffRsq =  [rt.TEfficiency("effRsq_MR%i"%xCut,"Rsq efficiency;R^{2};efficiency",len(y)-1,y) for xCut in xCuts]
    [pEff.SetStatisticOption(rt.TEfficiency.kFCP) for pEff in pEffRsq]
    
    entry = -1
    while True:
        entry = elist.Next()
        #if entry>10000: break
        if entry == -1: break
        tree.GetEntry(entry)
        bNum = any([tree.HLTDecision[triggerDict[num]] for num in options.numerator.split(',')]) and any([tree.HLTDecision[triggerDict[denom]] for denom in options.denominator.split(',') ])
        bDenom = any([tree.HLTDecision[triggerDict[denom]] for denom in options.denominator.split(',') ])
        if bDenom:            
            for pDenom in pDenom2D:
                pDenom.Fill(tree.MR,tree.Rsq)
            for pDenom,yCut in zip(pDenomMR,yCuts):
                if tree.Rsq>yCut: pDenom.Fill(tree.MR)
            for pDenom,xCut in zip(pDenomRsq,xCuts):
                if tree.MR>xCut: pDenom.Fill(tree.Rsq)
            if bNum:                
                for pNum in pNum2D:
                    pNum.Fill(tree.MR,tree.Rsq)
                for pNum,yCut in zip(pNumMR,yCuts):
                    if tree.Rsq>yCut: pNum.Fill(tree.MR)
                for pNum,xCut in zip(pNumRsq,xCuts):
                    if tree.MR>xCut: pNum.Fill(tree.Rsq)
                        
            for pEff in pEff2D:
                pEff.Fill(bNum,tree.MR,tree.Rsq)
            for pEff,yCut in zip(pEffMR,yCuts):
                if tree.Rsq>yCut: pEff.Fill(bNum,tree.MR)
            for pEff,xCut in zip(pEffRsq,xCuts):
                if tree.MR>xCut: pEff.Fill(bNum,tree.Rsq)


    setStyle()
    c = rt.TCanvas("c","c",500,400)
    c.SetRightMargin(0.15)
    
    colors = [rt.kViolet,rt.kRed,rt.kGreen,rt.kBlue,rt.kBlack,rt.kGray,rt.kOrange]
    for pEff in pEff2D:     
        pEff.Draw("colztext")
        rt.gPad.Update()
        pEff.GetPaintedHistogram().GetZaxis().SetTitle("efficiency")       
        pEff.GetPaintedHistogram().SetMaximum(1)
        pEff.GetPaintedHistogram().SetMinimum(0)   
        l = rt.TLatex()
        l.SetTextAlign(11)
        l.SetTextSize(0.045)
        l.SetTextFont(42)
        l.SetNDC()
        l.DrawLatex(0.12,0.92,"CMS preliminary")
        l.DrawLatex(0.6,0.92,"13 TeV (%.0f pb^{-1})"%options.lumi)
        l.SetTextSize(0.02)
        l.SetTextFont(42)
        l.DrawLatex(0.12,0.85,"signal:       %s"%(' || '.join(options.numerator.split(','))))
        l.DrawLatex(0.12,0.82,"reference:  %s"%(' || '.join(options.denominator.split(','))))
        c.Print(options.outDir+"/"+'_'.join(options.numerator.split(','))+"_"+'_'.join(options.denominator.split(','))+"_"+pEff.GetName().replace(".","p")+".pdf")
        c.Print(options.outDir+"/"+'_'.join(options.numerator.split(','))+"_"+'_'.join(options.denominator.split(','))+"_"+pEff.GetName().replace(".","p")+".C")

        
    sigmoid = rt.TF1("sigmoidMR","[0]/(1.0+exp(-(x-[1])/[2]))",0.,1200.)
    sigmoid.SetParameter(0,1)
    sigmoid.SetParLimits(0,0,1)
    sigmoid.SetParameter(1,200)
    sigmoid.SetParameter(2,30)
    first = True
    for pEff,color,yCut in zip(pEffMR,colors,yCuts):
        #pEff.SetLineColor(color)
        pEff.SetMarkerSize(0.8)
        pEff.SetMarkerStyle(20)
        pEff.Draw("apez")
        pEff.Fit(sigmoid,"I")
        rt.gPad.Update()        
        #pEff.GetPaintedHistogram().GetXaxis().SetRangeUser(0,1200)
        pEff.GetPaintedGraph().SetMarkerStyle(8)
        pEff.GetPaintedGraph().SetMarkerSize(20)        
        pEff.GetPaintedGraph().SetMinimum(0)
        pEff.GetPaintedGraph().SetMaximum(1.3)
        rt.gPad.Update()
        #if first:
        #    pEff.Draw("pe")
        #    first = False
        #else:
        #    pEff.Draw("pesame")        
        l = rt.TLatex()
        l.SetTextAlign(11)
        l.SetTextSize(0.045)
        l.SetTextFont(42)
        l.SetNDC()
        l.DrawLatex(0.12,0.92,"CMS preliminary")
        l.DrawLatex(0.6,0.92,"13 TeV (%.0f pb^{-1})"%options.lumi)
        l.SetTextFont(52)
        l.DrawLatex(0.7,0.75,"R^{2} > %.2f"%yCut)
        l.SetTextSize(0.02)
        l.SetTextFont(42)        
        l.DrawLatex(0.12,0.85,"signal:       %s"%(' || '.join(options.numerator.split(','))))
        l.DrawLatex(0.12,0.82,"reference:  %s"%(' || '.join(options.denominator.split(','))))
        c.Print(options.outDir+"/"+'_'.join(options.numerator.split(','))+"_"+'_'.join(options.denominator.split(','))+"_"+pEff.GetName().replace(".","p")+".pdf")
        c.Print(options.outDir+"/"+'_'.join(options.numerator.split(','))+"_"+'_'.join(options.denominator.split(','))+"_"+pEff.GetName().replace(".","p")+".C")


    sigmoid = rt.TF1("sigmoidRsq","[0]/(1.0+exp(-(x-[1])/[2]))",0.,1.5)
    sigmoid.SetParameter(0,1)
    sigmoid.SetParLimits(0,0,1)
    sigmoid.SetParameter(1,0.25)
    sigmoid.SetParameter(2,0.05)
    first = True
    for pEff,color,xCut in zip(pEffRsq,colors,xCuts):
        rt.gPad.Update()
        pEff.Fit(sigmoid,"I")
        #pEff.SetLineColor(color)
        pEff.SetMarkerSize(0.8)
        pEff.SetMarkerStyle(20)
        pEff.Draw("apez")
        rt.gPad.Update()
        #pEff.GetPaintedHistogram().GetXaxis().SetRangeUser(0,1.5)
        pEff.GetPaintedGraph().SetMinimum(0)
        pEff.GetPaintedGraph().SetMaximum(1.3)
        #if first:
        #    pEff.Draw("pe")
        #    first = False
        #else:
        #    pEff.Draw("pesame")
        l = rt.TLatex()
        l.SetTextAlign(11)
        l.SetTextSize(0.045)
        l.SetTextFont(42)
        l.SetNDC()
        l.DrawLatex(0.12,0.92,"CMS preliminary")
        l.DrawLatex(0.6,0.92,"13 TeV (%.0f pb^{-1})"%options.lumi)
        l.SetTextFont(52)
        l.DrawLatex(0.7,0.75,"M_{R} > %i"%xCut)
        l.SetTextSize(0.02)
        l.SetTextFont(42)        
        l.DrawLatex(0.12,0.85,"signal:       %s"%(' || '.join(options.numerator.split(','))))
        l.DrawLatex(0.12,0.82,"reference:  %s"%(' || '.join(options.denominator.split(','))))
        c.Print(options.outDir+"/"+'_'.join(options.numerator.split(','))+"_"+'_'.join(options.denominator.split(','))+"_"+pEff.GetName().replace(".","p")+".pdf")
        c.Print(options.outDir+"/"+'_'.join(options.numerator.split(','))+"_"+'_'.join(options.denominator.split(','))+"_"+pEff.GetName().replace(".","p")+".C")
