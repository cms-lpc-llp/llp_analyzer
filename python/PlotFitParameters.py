#! /usr/bin/env python
from optparse import OptionParser

import ROOT as rt
import rootTools
from framework import Config
import os.path
from array import *
import sys

boxes = ["MuMultiJet","EleMultiJet","MultiJet"]

boxDict = {"MuMultiJet":0,"EleMultiJet":1,"MultiJet":2}

parDict = {"MR0_TTj0b":"M_{R,0b-tag}^{0}","R0_TTj0b":"R^{2}_{0,0b-tag}",
           "MR0_TTj1b":"M_{R,1b-tag}^{0}","R0_TTj1b":"R^{2}_{0,1b-tag}",
           "MR0_TTj2b":"M_{R,2b-tag}^{0}","R0_TTj2b":"R^{2}_{0,2b-tag}",
           "MR1_TTj3b":"M_{R,#geq3b-tag}^{1}",
           "n_TTj0b":"n_{0b-tag}","b_TTj0b":"b_{0b-tag}",
           "n_TTj1b":"n_{1b-tag}","b_TTj1b":"b_{1b-tag}",
           "n_TTj2b":"n_{2b-tag}","b_TTj2b":"b_{2b-tag}",
           "Ntot_TTj0b":"N_{0b-tag}","Ntot_TTj1b":"N_{1b-tag}","Ntot_TTj2b":"N_{2b-tag}","Ntot_TTj3b":"N_{#geq3b-tag}"}
    
texParDict = {"MR0_TTj0b":"M_{R,0\\textrm{b-tag}}^{0}","R0_TTj0b":"R^{2}_{0,0\\textrm{b-tag}}",
           "MR0_TTj1b":"M_{R,1\\textrm{b-tag}}^{0}","R0_TTj1b":"R^{2}_{0,1\\textrm{b-tag}}",
           "MR0_TTj2b":"M_{R,2\\textrm{b-tag}}^{0}","R0_TTj2b":"R^{2}_{0,2\\textrm{b-tag}}",
           "MR1_TTj3b":"M_{R,\geq3\\textrm{b-tag}}^{1}",
           "n_TTj0b":"n_{0\\textrm{b-tag}}","b_TTj0b":"b_{0\\textrm{b-tag}}",
           "n_TTj1b":"n_{1\\textrm{b-tag}}","b_TTj1b":"b_{1\\textrm{b-tag}}",
           "n_TTj2b":"n_{2\\textrm{b-tag}}","b_TTj2b":"b_{2\\textrm{b-tag}}",
           "Ntot_TTj0b":"N_{0\\textrm{b-tag}}","Ntot_TTj1b":"N_{1\\textrm{b-tag}}","Ntot_TTj2b":"N_{2\\textrm{b-tag}}","Ntot_TTj3b":"N_{\geq3\\textrm{b-tag}}"}
    
def box_sort_key(parFile):
    label = parFile[0]
    box, fit = label
    boxNum = boxDict[box]
    boxNum*=2
    if fit == "Full": boxNum+=1
    return boxNum

def param_sort_key(parName):    
    keyDict = {"Ntot_TTj0b":0,
               "MR0_TTj0b":1,
               "R0_TTj0b":2,
               "b_TTj0b":3,
               "n_TTj0b":4,
               "Ntot_TTj1b":5,
               "MR0_TTj1b":6,
               "R0_TTj1b":7,
               "b_TTj1b":8,
               "n_TTj1b":9,
               "Ntot_TTj2b":10,
               "MR0_TTj2b":11,
               "R0_TTj2b":12,
               "b_TTj2b":13,
               "n_TTj2b":14,
               "Ntot_TTj3b":15,
               "MR1_TTj3b":16}
        
    return keyDict[parName]


def readFitResult(label, fileName):
    box, fit = label
    rootFile = rt.TFile.Open(fileName)
    #read the variables from the workspace
    w = rootFile.Get('w'+box)    
    if w.obj("fitresult_extRazorPdf_data_obs") != None:
        fr = w.obj("fitresult_extRazorPdf_data_obs")
    elif w.obj("nll_extRazorPdf_data_obs") != None:
        fr = w.obj("nll_extRazorPdf_data_obs")
    elif w.obj("fitresult_extRazorPdf_data_obs_with_constr") != None:
        fr = w.obj("fitresult_extRazorPdf_data_obs_with_constr")
    elif w.obj("nll_extRazorPdf_data_obs_with_constr") != None:
        fr = w.obj("nll_extRazorPdf_data_obs_with_constr")
    fr.Print("v")
    #get the final values from the fit
    parList = fr.floatParsFinal()
    fitPars = {}
    for p in rootTools.RootIterator.RootIterator(parList):
        fitPars[p.GetName().replace('_%s'%box,'')] = [p.getVal(),p.getError(),p.getMin(),p.getMax()]
    return fitPars


def getTexTable(parFiles,parNameList):        
    texTable = ' Parameter'
    
    for label, parValErr in sorted(parFiles.iteritems(),key=box_sort_key):     
        box, fit = label
        texTable += ' & %s %s Fit' %(box,fit)
    texTable += ' \\\ \n'
    for parName in sorted(parNameList,key=param_sort_key):
        texTable += ' $%s$'% texParDict[parName]
        
        for label, parValErr in sorted(parFiles.iteritems(),key=box_sort_key):
            box, fit = label
            parVal = parValErr[parName][0]
            parErr = parValErr[parName][1]
            parMin = parValErr[parName][2]
            parMax = parValErr[parName][3]
            if 'MR1' in parName:
                addToTable = ' & $%.2e\pm%.2e$'%(parVal,parErr)
                addToTable = addToTable.replace('e-03','\\times10^{-3}')
                addToTable = addToTable.replace('e-04','\\times10^{-4}')
                texTable += addToTable
            else:
                texTable += ' & $%.2f\pm%.2f$'%(parVal,parErr)
        texTable += ' \\\ \n'
    return texTable

def getParHisto(parName,parFiles):
    numBins = 0
    
    for label, parValErr in sorted(parFiles.iteritems(),key=box_sort_key):
        try:
            parVal = parValErr[parName][0]
        except KeyError:
            continue
        numBins+=1
    parHisto = rt.TH1D(parName,parDict[parName],numBins,0,numBins)
    binNum = 0
    setAxis = parHisto.GetXaxis()
    setAxis.SetTitle("")
    for label, parValErr in sorted(parFiles.iteritems(),key=box_sort_key):
        box, fit = label
        try:
            parVal = parValErr[parName][0]
            parErr = parValErr[parName][1]
            parMin = parValErr[parName][2]
            parMax = parValErr[parName][3]
            binNum +=1
            parHisto.SetBinContent(binNum,parVal)
            parHisto.SetBinError(binNum,parErr)
        except KeyError:
            continue
            
        parHisto.SetMarkerStyle(8)
        #parHisto.SetMarkerColor(rt.kViolet)
        #parHisto.SetLineColor(rt.kAzure)
        parHisto.SetMarkerColor(rt.kCyan+2)
        parHisto.SetLineColor(rt.kCyan+2)
        parHisto.SetMarkerSize(1.2)
        setAxis.SetBinLabel(binNum,"%s %s Fit"%(box,fit))
    return parHisto

def getBoxHisto(box,parFiles):
    numPars = len(parFiles[box,'Full'].keys())
    boxHisto = rt.TH1D(box,box,2*numPars,0,2*numPars)
    setAxis = boxHisto.GetXaxis()
    
    binNum=0
    
    for label, parValErr in sorted(parFiles.iteritems(),key=box_sort_key):        
        if label[0]!=box: continue
        box, fit = label

        for parName in parValErr.keys():
            parVal = parValErr[parName][0]
            parErr = parValErr[parName][1]
            parMin = parValErr[parName][2]
            parMax = parValErr[parName][3]
            binNum +=1
            boxHisto.SetBinContent(binNum,parVal)
            boxHisto.SetBinError(binNum,parErr) 
            setAxis.SetBinLabel(binNum,"%s %s"%(parName, fit))
        boxHisto.SetMarkerStyle(8)
        boxHisto.SetMarkerColor(rt.kViolet)
        boxHisto.SetLineColor(rt.kAzure)
        boxHisto.SetMarkerSize(1.2)
            
    return boxHisto

def setstyle():
    rt.gStyle.SetCanvasBorderMode(0)
    rt.gStyle.SetCanvasColor(rt.kWhite)
    rt.gStyle.SetCanvasDefH(400) 
    rt.gStyle.SetCanvasDefW(600) 
    rt.gStyle.SetCanvasDefX(0)  
    rt.gStyle.SetCanvasDefY(0)

    rt.gStyle.SetPadBorderMode(0)
    rt.gStyle.SetPadColor(rt.kWhite)
    rt.gStyle.SetPadGridX(False)
    rt.gStyle.SetPadGridY(False)
    rt.gStyle.SetGridColor(0)
    rt.gStyle.SetGridStyle(3)
    rt.gStyle.SetGridWidth(1)
    
    rt.gStyle.SetFrameBorderMode(0)
    rt.gStyle.SetFrameBorderSize(1)
    rt.gStyle.SetFrameFillColor(0)
    rt.gStyle.SetFrameFillStyle(0)
    rt.gStyle.SetFrameLineColor(1)
    rt.gStyle.SetFrameLineStyle(1)
    rt.gStyle.SetFrameLineWidth(1)
    
    rt.gStyle.SetPaperSize(20,26)
    rt.gStyle.SetPadTopMargin(0.1)
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.2)
    rt.gStyle.SetPadLeftMargin(0.15)
    
    rt.gStyle.SetTitleFont(42,"xyz") 
    rt.gStyle.SetTitleFont(42," ")   
    rt.gStyle.SetTitleSize(0.06,"xyz")
    rt.gStyle.SetTitleSize(0.06," ")  
    rt.gStyle.SetLabelFont(42,"xyz")
    rt.gStyle.SetLabelSize(0.06,"xyz")
    rt.gStyle.SetLabelColor(1,"xyz")
    rt.gStyle.SetLabelOffset(0.015,"x")
    rt.gStyle.SetTextFont(42)
    rt.gStyle.SetTextSize(0.08)
    rt.gStyle.SetStatFont(132)
    rt.gStyle.SetMarkerStyle(8)
    #rt.gStyle.SetHistLineWidth((rt.Width_t) 1.85)
    #rt.gStyle.SetLineStyleString(2,"[12 12]")
    rt.gStyle.SetErrorX(0.2)
    rt.gStyle.SetOptTitle(1)
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptFit(11111111)
    rt.gStyle.SetPadTickX(1)
    rt.gStyle.SetPadTickY(1)

    def set_palette(name="default", ncontours=255):
        """Set a color palette from a given RGB list
        stops, red, green and blue should all be lists of the same length
        see set_decent_colors for an example"""

        if name == "gray" or name == "grayscale":
            stops = [0.00, 0.34, 0.61, 0.84, 1.00]
            red   = [1.00, 0.95, 0.95, 0.65, 0.15]
            green = [1.00, 0.85, 0.7, 0.5, 0.3]
            blue  = [0.95, 0.6, 0.3, 0.45, 0.65]
            # elif name == "whatever":
            # (define more palettes)
        elif name == "chris":
            stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
            red =   [ 1.0,   0.95,  0.95,  0.65,   0.15 ]
            green = [ 1.0,  0.85, 0.7, 0.5,  0.3 ]
            blue =  [ 0.95, 0.6 , 0.3,  0.45, 0.65 ]
        else:
            # default palette, looks cool
            stops = [0.00, 0.34, 0.61, 0.84, 1.00]
            red   = [0.00, 0.00, 0.87, 1.00, 0.51]
            green = [0.00, 0.81, 1.00, 0.20, 0.00]
            blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

        s = array('d', stops)
        r = array('d', red)
        g = array('d', green)
        b = array('d', blue)

        npoints = len(s)
        rt.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
        rt.gStyle.SetNumberContours(ncontours)

    set_palette("chris")

    rt.gStyle.cd()

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default=None,
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outdir",default="./",type="string",
                  help="Output directory to store datasets")
    
    (options,args) = parser.parse_args()

    cfg = Config.Config(options.config)
    print 'Input files are %s' % ', '.join(args)

    setstyle()
    
    labels = {}
    for f in args:
        if f.lower().endswith('.root'):
            if 'Sideband' in f:
                fit = 'Sideband'
            else:
                fit = 'Full'
            for boxTest in boxes:
                if '/%s/'%boxTest in f: box = boxTest
            if not labels.has_key((box,fit)):
                labels[box,fit] = f            
        else:
            "File '%s' of unknown type. Looking for .root files only" % f

    parFiles = {}
    for label, files in labels.iteritems():
        parFiles[label] = readFitResult(label, files)

    parNameSet = set([])
    for label, files in labels.iteritems():
        parNameSet = parNameSet.union(set(parFiles[label].keys()))

    
    rt.gStyle.SetPadTopMargin(0.11)
    c = rt.TCanvas("c","c",600,400)
    c.SetLogy(0)
    
    # for box in boxes:
    #     boxHisto = getBoxHisto(box,parFiles)
    #     boxHisto.Draw('E1')
    #     c.Print("%s/%s.pdf"%(options.outdir,box))


    for parName in parNameSet:
        parHisto = getParHisto(parName,parFiles)
        tlines = []

        if parName.find("n_")!=-1 or parName.find("b_")!=-1:
            rt.gPad.SetLogy()
            parHisto.SetMinimum(.3*parHisto.GetBinContent(parHisto.GetMinimumBin()))
            parHisto.SetMaximum(3*parHisto.GetBinContent(parHisto.GetMaximumBin()))
        else:
            rt.gPad.SetLogy(0)
        parHisto.Draw('E1')
        rt.gPad.Update()
        
        c.Update()
        for i in xrange(2, parHisto.GetNbinsX()+1,2):
            
            if parName.find("n_")!=-1 or parName.find("b_")!=-1:
                tlines.append(rt.TLine(i, parHisto.GetMinimum(), i, parHisto.GetMaximum()))
            else:
                tlines.append(rt.TLine(i, rt.gPad.GetFrame().GetY1(), i, rt.gPad.GetFrame().GetY2()))

        for tline in tlines:
            tline.Draw("same")
        c.Print("%s/%s.pdf"%(options.outdir,parName))
        c.Print("%s/%s.C"%(options.outdir,parName))


    parNameList = list(parNameSet)
    texTable = getTexTable(parFiles,parNameList)
    print texTable
    
    
