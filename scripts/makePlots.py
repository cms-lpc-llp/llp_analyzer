#! /usr/bin/env python
import ROOT as rt
import os.path
import sys, glob, re
from array import *
import math

def setstyle():
    # For the canvas:
    rt.gStyle.SetCanvasBorderMode(0)
    rt.gStyle.SetCanvasColor(rt.kWhite)
    rt.gStyle.SetCanvasDefH(400) #Height of canvas
    rt.gStyle.SetCanvasDefW(600) #Width of canvas
    rt.gStyle.SetCanvasDefX(0)   #POsition on screen
    rt.gStyle.SetCanvasDefY(0)
    
    # For the Pad:
    rt.gStyle.SetPadBorderMode(0)
    # rt.gStyle.SetPadBorderSize(Width_t size = 1)
    rt.gStyle.SetPadColor(rt.kWhite)
    rt.gStyle.SetPadGridX(False)
    rt.gStyle.SetPadGridY(False)
    rt.gStyle.SetGridColor(0)
    rt.gStyle.SetGridStyle(3)
    rt.gStyle.SetGridWidth(1)
    
    # For the frame:
    rt.gStyle.SetFrameBorderMode(0)
    rt.gStyle.SetFrameBorderSize(1)
    rt.gStyle.SetFrameFillColor(0)
    rt.gStyle.SetFrameFillStyle(0)
    rt.gStyle.SetFrameLineColor(1)
    rt.gStyle.SetFrameLineStyle(1)
    rt.gStyle.SetFrameLineWidth(1)
    
    # set the paper & margin sizes
    rt.gStyle.SetPaperSize(20,26)
    rt.gStyle.SetPadTopMargin(0.09)
    #rt.gStyle.SetPadRightMargin(0.065)
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.17)
    
    # use large Times-Roman fonts
    rt.gStyle.SetTitleFont(42,"xyz")  # set the all 3 axes title font
    rt.gStyle.SetTitleFont(42," ")    # set the pad title font
    rt.gStyle.SetTitleSize(0.065,"xyz") # set the 3 axes title size
    rt.gStyle.SetTitleSize(0.055," ")   # set the pad title size
    rt.gStyle.SetLabelFont(42,"xyz")
    rt.gStyle.SetLabelSize(0.065,"xyz")
    rt.gStyle.SetLabelColor(1,"xyz")
    rt.gStyle.SetTextFont(42)
    rt.gStyle.SetTextSize(0.08)
    rt.gStyle.SetStatFont(42)
    
    # use bold lines and markers
    #rt.gStyle.SetMarkerStyle(20)
    rt.gStyle.SetLineStyleString(2,"[12 12]") # postscript dashes
    
    #..Get rid of X error bars
    #rt.gStyle.SetErrorX(0.001)
    
    # do not display any of the standard histogram decorations
    #rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptFit(1111)
    rt.gStyle.SetStatY(0.85)        
    rt.gStyle.SetStatX(0.92)                
    rt.gStyle.SetStatW(0.15)                
    rt.gStyle.SetStatH(0.15)                
    
    # put tick marks on top and RHS of plots
    rt.gStyle.SetPadTickX(1)
    rt.gStyle.SetPadTickY(1)
    
    ncontours = 999
    
    stops = [ 0.00, 0.1, 0.25, 0.65, 1.00 ]
    #stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
    red =   [ 1.0,   0.95,  0.95,  0.65,   0.15 ]
    green = [ 1.0,  0.85, 0.7, 0.5,  0.3 ]
    blue =  [ 0.95, 0.6 , 0.3,  0.45, 0.65 ]
    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)
        
    npoints = len(s)
    #rt.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    rt.gStyle.SetNumberContours(ncontours)
   
    rt.gStyle.cd()

def makeHistos(fileName):
    
    rt.gStyle.SetOptStat(0)
    c  = rt.TCanvas("c","c",1,1,650,376)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    c.Update()
    tleg = rt.TLegend(0.68, 0.75, 0.9, 0.9)
    
    tfile = rt.TFile(fileName)
    tree = tfile.Get('HbbRazor')

    histoData = {('mbb',''):
                 ['mbb','BTagCSV data',100,0,1000,"m_{bb} [GeV]","Events",""],
                 ('MR',''):
                 ['MR','BTagCSV data',100,0,4000,"M_{R} [GeV]","Events",""],
                 ('Rsq',''):
                 ['Rsq','BTagCSV data',100,0,1.2,"R^{2}","Events",""],
                 ('Rsq:MR',''):
                 ['RsqMR','BTagCSV data',100,0,4000,100,0,1.2,"M_{R} [GeV]","R^{2}","colz"],
                }

    histo = {}
    
    for (var,wgt), histoDatum in histoData.iteritems():
        if len(var.split(":"))==1:
            histo[var,wgt] = rt.TH1D(histoDatum[0],histoDatum[1],histoDatum[2],histoDatum[3],histoDatum[4])
            histo[var,wgt].SetXTitle(histoDatum[5])
            histo[var,wgt].SetYTitle(histoDatum[6])
            tree.Project(histo[var,wgt].GetName(),var,wgt)
            #histo[var,wgt].Scale(1./histo[var,wgt].Integral())
            histo[var,wgt].Draw(histoDatum[7])
        elif len(var.split(":"))==2:
            histo[var,wgt] = rt.TH2D(histoDatum[0],histoDatum[1],histoDatum[2],histoDatum[3],histoDatum[4],histoDatum[5],histoDatum[6],histoDatum[7])
            histo[var,wgt].SetXTitle(histoDatum[8])
            histo[var,wgt].SetYTitle(histoDatum[9])
            tree.Project(histo[var,wgt].GetName(),var,wgt)
            #histo[var,wgt].Scale(1./histo[var,wgt].Integral())
            histo[var,wgt].Draw(histoDatum[10])
        

        c.Print(histo[var,wgt].GetName()+".png")
        c.Print(histo[var,wgt].GetName()+".pdf")
        c.Print(histo[var,wgt].GetName()+".C")
        
    outputFile = rt.TFile.Open("%s_plots.root"%fileName.split('.root')[0],"recreate")
    for var,histogram in histo.iteritems():
        histogram.Write()
    

    
if __name__ == '__main__':
    setstyle()
    rt.gStyle.SetTitleFont(42,"xyz")  # set the all 3 axes title font
    rt.gStyle.SetTitleFont(42," ")    # set the pad title font
    rt.gStyle.SetLabelFont(42,"xyz")
    rt.gStyle.SetTextFont(42)
    rt.gStyle.SetNumberContours(999)
    fileName = sys.argv[1]
    makeHistos(fileName)
