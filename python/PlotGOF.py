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
from PlotFit import *

def getGOFHistos(varName,toyTree):
    toyTree.GetEntry(0)
    var_data = eval('toyTree.%s_%s'%(varName,box))
        
    toyTree.Draw('%s>>htest%s'%('%s_%s'%(varName,box),'%s_%s'%(varName,box)))
    htemp = rt.gPad.GetPrimitive("htest%s"%('%s_%s'%(varName,box)))
    rms = htemp.GetRMS()
    mean = htemp.GetMean()
    
    xmax = max(mean+3.*rms,var_data+1)
    
    xmin = int(max(0,htemp.GetXaxis().GetXmin()))
    
    h = rt.TH1D('h_%s'%varName,'h_%s'%varName,70,xmin,xmax)
    h_cut = rt.TH1D('h_%s_cut'%varName,'h_%s_cut'%varName,70,xmin,xmax)
    toyTree.Project('h_%s'%varName,'%s_%s'%(varName,box))
    toyTree.Project('h_%s_cut'%varName,'%s_%s'%(varName,box),'%s_%s>%f'%(varName,box,var_data))

    return h, h_cut, var_data

def setHist(h_data,xTitle,yTitle,color=rt.kBlack):
    h_data.SetMarkerColor(color)
    h_data.SetLineColor(color)
    h_data.SetMarkerStyle(20)
    h_data.SetLineColor(color)
    h_data.GetXaxis().SetTitle(xTitle)
    h_data.GetYaxis().SetTitle(yTitle)
    #h_data.GetXaxis().SetLabelOffset(0.16)
    h_data.GetXaxis().SetLabelSize(0.05)
    h_data.GetYaxis().SetLabelSize(0.05)
    h_data.GetXaxis().SetTitleSize(0.05)
    h_data.GetYaxis().SetTitleSize(0.05)
    #h_data.GetXaxis().SetTitleOffset(0.8)
    #h_data.GetYaxis().SetTitleOffset(0.7)
    h_data.SetMaximum(pow(h_data.GetMaximum(),1.7))
    h_data.SetMinimum(2e-1)
    return h_data
    
def print1DGOF(c,rootFile,h,h_cut,func,chi2_data,printName,xTitle,yTitle,lumiLabel="",boxLabel="",plotLabel="",isData=False,tLeg=None):

    c.SetLogy(1)
    setHist(h,xTitle,yTitle)
    
    #h.SetLineWidth(2)    
    h.Draw("pe")
    h_cut.SetLineColor(rt.kViolet-10)
    h_cut.SetFillColor(rt.kViolet-10)
    h_cut.Draw("histsame")
    h.Draw("pesame")
    if func!=None:
        npoints = 1000
        listx = []
        listy = []
        xmin = h.GetXaxis().GetXmin()
        xmax = h.GetXaxis().GetXmax()
        for ix in range(0,npoints+1):
            x = xmin+ix*(xmax-xmin)/npoints
            listx.append(x)
            #listy.append(h.GetEntries()*func.Eval(x)/func.Integral(xmin,xmax))     
        #graph = rt.TGraph(npoints, array('d',listx), array('d',listy))
        #graph.SetLineColor(rt.kRed)
        #graph.SetLineWidth(2)
        #graph.Draw("same")
    
    pvalue = h_cut.Integral(0,h_cut.GetNbinsX()+1)/h.Integral(0,h.GetNbinsX()+1)
    for i in range(1,h_cut.GetNbinsX()+1):
        if h_cut.GetXaxis().GetBinLowEdge(i) < chi2_data:
            h_cut.SetBinContent(i,0)
    
    tlineObs = rt.TArrow(chi2_data,0,chi2_data,h.GetBinContent(h.GetMaximumBin()),0.04,"<")
    tlineObs.SetLineColor(rt.kBlack)
    tlineObs.SetLineWidth(4)
    tlineObs.Draw()

        
    tLeg = rt.TLegend(0.69,0.5,0.89,0.8)
    tLeg.SetLineColor(rt.kWhite)
    tLeg.AddEntry(h,"Toy Data","lep")
    tLeg.AddEntry(tlineObs,"Observed = %.1f"%chi2_data,"l")
    tLeg.AddEntry(h_cut,"p-value = %.2f"%(pvalue),"f")
            
    tLeg.Draw("same")

    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    if isData:
        l.DrawLatex(0.12,0.91,"CMS preliminary")
    else:
        l.DrawLatex(0.12,0.91,"CMS simulation")
    l.DrawLatex(0.7,0.91,"%s"%lumiLabel)
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
    parser.add_option('-t','--input-toy-file',dest="inputToyFile", default=None,type="string",
                  help="input toy file")
    parser.add_option('--data',dest="isData", default=False,action='store_true',
                  help="changes plots for data")
    parser.add_option('--fit-region',dest="fitRegion",default="Full",type="string",
                  help="Fit region")
    
    (options,args) = parser.parse_args()
     
    setStyle()
    
    box = options.box
    lumi = options.lumi
    cfg = Config.Config(options.config)
    fitRegion = options.fitRegion

    
    toyTree = None
    if options.inputToyFile is not None:
        toyFiles = options.inputToyFile.split(',')
        toyTree = rt.TChain("myTree")
        for toyFile in toyFiles:
            toyTree.Add(toyFile)
    
    c = rt.TCanvas('c','c',500,400)    
    c.SetLeftMargin(0.12) 
    c.SetBottomMargin(0.12)
    rootFile = rt.TFile.Open(options.outDir + '/' + 'Plots_%s'%box + '.root','recreate')
    tdirectory = rootFile.GetDirectory(options.outDir)
    if tdirectory==None:
        rootFile.mkdir(options.outDir)
        tdirectory = rootFile.GetDirectory(options.outDir)

    if options.isData:
        dataString = "Data"
    else:
        dataString = "Sim. Data"

    eventsLabel = "Toy Datasets"
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
        
    h_n2llr, h_n2llr_cut, n2llr_data = getGOFHistos('n2llr',toyTree)
    h_chi2, h_chi2_cut, chi2_data = getGOFHistos('chi2',toyTree)

    h_n2llr_btag = [None for k in range(0,len(z)-1)]
    h_n2llr_cut_btag = [None for k in range(0,len(z)-1)]
    n2llr_data_btag = [None for k in range(0,len(z)-1)]
    h_chi2_btag = [None for k in range(0,len(z)-1)]
    h_chi2_cut_btag = [None for k in range(0,len(z)-1)]
    chi2_data_btag = [None for k in range(0,len(z)-1)]
    
    if len(z)>2:
        for k in range(0,len(z)-1):
            ibtag = z[k]
            h_n2llr_btag[k], h_n2llr_cut_btag[k], n2llr_data_btag[k] = getGOFHistos('n2llr_%ibtag'%(ibtag),toyTree)
            h_chi2_btag[k], h_chi2_cut_btag[k], chi2_data_btag[k] = getGOFHistos('chi2_%ibtag'%(ibtag),toyTree)
            
    
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
    boxLabel = "razor %s %s %s Fit" % (box,btagLabel,fitRegion)

    chi2 = rt.RooRealVar('chi2','chi2',0,0,10000)
    ndof = rt.RooRealVar('ndof','ndof',1,0,200)
    chi2_pdf = rt.RooChiSquarePdf('chi2_pdf','chi2_pdf',chi2,ndof)
    ndof.setVal(140-17) #roughly ndof
    chi2_func = chi2_pdf.asTF(rt.RooArgList(chi2))
    
    print1DGOF(c,tdirectory,h_n2llr,h_n2llr_cut,chi2_func,n2llr_data,options.outDir+"/gof_n2llr_%s.pdf"%box,"-2 log #lambda",eventsLabel,lumiLabel,boxLabel,'',options.isData,None)
    print1DGOF(c,tdirectory,h_chi2,h_chi2_cut,chi2_func,chi2_data,options.outDir+"/gof_chi2_%s.pdf"%box,"#chi^{2}",eventsLabel,lumiLabel,boxLabel,'',options.isData,None)
    
    if len(z)>2:
        for k in range(0,len(z)-1):             
            newBtagLabel = '%i b-tag' %z[k]
            if z[k]==3 and z[-1]==4:
                newBtagLabel = '#geq %i b-tag' %z[k]                
            newBoxLabel = "razor %s %s %s Fit"%(box,newBtagLabel,fitRegion)
            ndof.setVal(34-5) #roughly ndof
            chi2_func = chi2_pdf.asTF(rt.RooArgList(chi2))
            print1DGOF(c,tdirectory,h_n2llr_btag[k],h_n2llr_cut_btag[k],chi2_func,n2llr_data_btag[k],options.outDir+"/gof_n2llr_%ibtag_%s.pdf"%(z[k],box),"-2 log #lambda",eventsLabel,lumiLabel,newBoxLabel,'',options.isData,None)
            print1DGOF(c,tdirectory,h_chi2_btag[k],h_chi2_cut_btag[k],chi2_func,chi2_data_btag[k],options.outDir+"/gof_chi2_%ibtag_%s.pdf"%(z[k],box),"#chi^{2}",eventsLabel,lumiLabel,newBoxLabel,'',options.isData,None)

 
