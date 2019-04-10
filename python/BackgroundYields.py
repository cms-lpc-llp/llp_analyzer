#Plotting script to visualize the yields from each background MC sample for a number of kinematic distributions
import sys, os
import string
import numpy as np
import itertools
import ROOT as rt
import rootTools
from array import * 

def draw1D(boxName, box, files, titleString, plotString, cutString, binning, outpath = ""):
    print("Producing yield plot for box "+boxName)

    #signal SMS info
    sigXSec = 0.002 #in picobarns
    nSig = 105964 #number of processed events in signal sample
    sigWeight = sigXSec / nSig
    
    #create histograms to stack
    histNames = ["TTJets","WJets","ZJetsToNuNu","DYJetsToLL","QCD"]
    sigName = "T1bbbb_2J_mGl-1500_mLSP-100"
    colors = [8,2,7,4,1]
    hists = [rt.TH1F(name, titleString, len(binning)-1, binning) for name in histNames]
    sigHist = rt.TH1F(sigName, titleString, len(binning) - 1, binning)
    for index, hist in enumerate(hists): 
        hist.SetFillColor(colors[index])
        hist.SetLineColor(1)
    #transparent fill for signal histogram
    sigHist.SetLineStyle(2)
    sigHist.SetLineColor(6)
    sigHist.SetLineWidth(4)

    stack = rt.THStack("stack", titleString)
    leg = rt.TLegend(0.6, 0.6, 0.9, 0.9)
    for index, tree in enumerate(box): 
        #determine which histogram to fill
        for name in histNames:
            if name in files[index]:
                tree.Draw(plotString+">>+"+name, "weight*("+cutString+")")
        if sigName in files[index]:
            tree.Draw(plotString+">>+"+sigName, str(sigWeight)+"*("+cutString+")")
    hists.sort(key = lambda h: h.Integral()) #sort by histogram integral
    #print out histogram yields
    for hist in hists: hist.Print("all")
    for hist in hists: 
        for bin in range(1, len(binning)): hist.SetBinContent(bin, hist.GetBinContent(bin)/hist.GetBinWidth(bin)) #divide each histogram bin by its width
        stack.Add(hist)
        if hist.Integral() > 0: leg.AddEntry(hist, hist.GetName())
    for bin in range(1, len(binning)): sigHist.SetBinContent(bin, sigHist.GetBinContent(bin)/sigHist.GetBinWidth(bin))
    if sigHist.Integral() > 0: leg.AddEntry(sigHist, sigHist.GetName())

    #plot and print
    c = rt.TCanvas("c", titleString, 1000, 1000)
    c.SetLogy()
    stack.SetTitle(titleString+", "+boxName+" Box")
    stack.Draw("h")
    stack.GetHistogram().GetXaxis().SetTitle(titleString)
    sigHist.Draw("same")
    leg.Draw()
    c.Print(outpath+"/Background"+titleString.replace(" ","").replace("{","").replace("}","").replace("^","").replace("_","")+boxName+".pdf")

##initialization

#check that input file was specified
if len(sys.argv) < 2: 
    print 'Usage: BackgroundYields.py filename1.root ...'
    exit()

#switch ROOT to batch mode
rt.gROOT.SetBatch()

#get the path of this script and load the razor library
pathname = os.path.dirname(sys.argv[0]) 
fullpath = os.path.abspath(pathname) #get the full path of this script
outpath = fullpath+"/output"
libpath = fullpath+"/lib"
rt.gSystem.Load(libpath+"/libRazorRun2.so")

print 'Will create yield histograms for standard model background events in the chosen file(s)'

#load the TTrees from the input files
boxNames = [
        'MuEle',
        'MuMu',
        'EleEle',
        'MuMultiJet',
        'MuJet',
        'EleMultiJet',
        'EleJet',
        'MultiJet',
        'TwoBJet',
        'OneBJet',
        'ZeroBJet'
        ]
filenames = [sys.argv[i] for i in range(1, len(sys.argv))]
files = [rt.TFile(filename) for filename in filenames]
razorBoxes = dict((box, [file.Get(box) for file in files]) for box in boxNames)

##MR yields, Rsq > 0.15
#MRBins = [300, 350, 400, 450, 550, 700, 900, 1200, 1600, 2500, 4000]
MRBins = [400, 450, 550, 700, 900, 1200, 1600, 2500, 4000]
MRBinArray = array('d',MRBins)
MRTitle = "M_{R}"
MRPlot = "MR"
#MRCut = "Rsq > 0.15"
MRCut = "Rsq > 0.25"
for key in razorBoxes: draw1D(key, razorBoxes[key], filenames, MRTitle, MRPlot, MRCut, MRBinArray, outpath)

##Rsq yields, MR > 300
#RsqBins = [0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 0.64, 0.80, 1.5]
RsqBins = [0.25, 0.3, 0.41, 0.52, 0.64, 0.80, 1.5]
RsqBinArray = array('d', RsqBins)
RsqTitle = "R^{2}"
RsqPlot = "Rsq"
#RsqCut = "MR > 300"
RsqCut = "MR > 400"
for key in razorBoxes: draw1D(key, razorBoxes[key], filenames, RsqTitle, RsqPlot, RsqCut, RsqBinArray, outpath)
