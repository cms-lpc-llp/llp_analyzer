from optparse import OptionParser
import ROOT as rt
import sys
import rootTools
import glob
from math import *
import os
from array import *
import numpy as np
from scipy.interpolate import Rbf, interp1d
import itertools
from GChiPairs import gchipairs
import operator
from limits.SMSConfig import sms_models

toFix = []
def interpolate2D(hist,epsilon=1,smooth=0,diagonalOffset=0,fixLSP0=False,refHist=None):
    x = array('d',[])
    y = array('d',[])
    z = array('d',[])
    
    binWidth = float(hist.GetXaxis().GetBinWidth(1))
    
    for i in range(1, hist.GetNbinsX()+1):
        for j in range(1, hist.GetNbinsY()+1):
            if hist.GetBinContent(i,j)>0.:
                if refHist!=None and refHist.GetBinContent(i,j) > 0.:
                        x.append(hist.GetXaxis().GetBinCenter(i))
                        y.append(hist.GetYaxis().GetBinCenter(j))
                        z.append(rt.TMath.Log(hist.GetBinContent(i,j)/refHist.GetBinContent(i,j)))
                else:
                    x.append(hist.GetXaxis().GetBinCenter(i))
                    y.append(hist.GetYaxis().GetBinCenter(j))
                    z.append(rt.TMath.Log(hist.GetBinContent(i,j)))

    mgMin = hist.GetXaxis().GetBinCenter(1)
    mgMax = hist.GetXaxis().GetBinCenter(hist.GetNbinsX())
    mchiMin = hist.GetYaxis().GetBinCenter(1)
    mchiMax = hist.GetYaxis().GetBinCenter(hist.GetNbinsY())
    
    myX = np.linspace(mgMin, mgMax,int((mgMax-mgMin)/binWidth+1))
    myY = np.linspace(mchiMin, mchiMax, int((mchiMax-mchiMin)/binWidth+1))
    myXI, myYI = np.meshgrid(myX,myY)

    rbf = Rbf(x, y, z,function='multiquadric', epsilon=epsilon,smooth=smooth)
    myZI = rbf(myXI, myYI)
    
    rbf_nosmooth = Rbf(x, y, z, function='multiquadric',epsilon=epsilon,smooth=10)
    otherY = array('d',[mchiMin])
    lineXI, lineYI = np.meshgrid(myX,otherY)
    lineZI = rbf_nosmooth(lineXI, lineYI)

    for i in range(1, hist.GetNbinsX()+1):
        for j in range(1, hist.GetNbinsY()+1):
            xLow = hist.GetXaxis().GetBinCenter(i)
            yLow = hist.GetYaxis().GetBinCenter(j)
            if j==1 and fixLSP0:
                hist.SetBinContent(i,j,rt.TMath.Exp(lineZI[j-1][i-1]))
                continue
            if xLow >= yLow+diagonalOffset-binWidth/2.:
                if refHist!=None:                    
                    hist.SetBinContent(i,j,refHist.GetBinContent(i,j)*rt.TMath.Exp(myZI[j-1][i-1]))
                else:
                    hist.SetBinContent(i,j,rt.TMath.Exp(myZI[j-1][i-1]))
    return hist


def fix_hist_byhand(hist, model, box, clsType):
    if 'Obs' in clsType:
        for (mg,mchi) in gchipairs(model):            
            obs = hist.GetBinContent(hist.FindBin(mg,mchi))
            exp = xsecUL['Exp'].GetBinContent(xsecUL['Exp'].FindBin(mg,mchi))
            expPlus2 = xsecUL['ExpPlus2'].GetBinContent(xsecUL['ExpPlus2'].FindBin(mg,mchi))
            expPlus = xsecUL['ExpPlus'].GetBinContent(xsecUL['ExpPlus'].FindBin(mg,mchi))
            expMinus = xsecUL['ExpMinus'].GetBinContent(xsecUL['ExpMinus'].FindBin(mg,mchi))
            expMinus2 = xsecUL['ExpMinus2'].GetBinContent(xsecUL['ExpMinus2'].FindBin(mg,mchi))
            if hist.GetBinContent(hist.FindBin(mg,mchi))==0:
                if (mg,mchi) not in toFix:                
                    toFix.append((mg,mchi))
                        
def set_palette(name="default", ncontours=255):
    # For the canvas:
    rt.gStyle.SetCanvasBorderMode(0)
    rt.gStyle.SetCanvasColor(rt.kWhite)
    rt.gStyle.SetCanvasDefH(400) #Height of canvas
    rt.gStyle.SetCanvasDefW(500) #Width of canvas
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
    rt.gStyle.SetPadTopMargin(0.085)
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.17)
    
    # use large Times-Roman fonts
    rt.gStyle.SetTitleFont(42,"xyz")  # set the all 3 axes title font
    rt.gStyle.SetTitleFont(42," ")    # set the pad title font
    rt.gStyle.SetTitleSize(0.06,"xyz") # set the 3 axes title size
    rt.gStyle.SetTitleSize(0.06," ")   # set the pad title size
    rt.gStyle.SetLabelFont(42,"xyz")
    rt.gStyle.SetLabelSize(0.05,"xyz")
    rt.gStyle.SetLabelColor(1,"xyz")
    rt.gStyle.SetTextFont(42)
    rt.gStyle.SetTextSize(0.08)
    rt.gStyle.SetStatFont(42)
    
    # use bold lines and markers
    rt.gStyle.SetMarkerStyle(8)
    #rt.gStyle.SetHistLineWidth(1.85)
    rt.gStyle.SetLineStyleString(2,"[12 12]") # postscript dashes
    
    #..Get rid of X error bars
    rt.gStyle.SetErrorX(0.001)
    
    # do not display any of the standard histogram decorations
    rt.gStyle.SetOptTitle(1)
    rt.gStyle.SetOptStat(0)
    #rt.gStyle.SetOptFit(11111111)
    rt.gStyle.SetOptFit(0)
    
    # put tick marks on top and RHS of plots
    rt.gStyle.SetPadTickX(1)
    rt.gStyle.SetPadTickY(1)
    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.95, 0.95, 0.65, 0.15]
        green = [1.00, 0.85, 0.7, 0.5, 0.3]
        blue  = [0.95, 0.6, 0.3, 0.45, 0.65]
    elif name == "chris":
	stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
	red =   [ 1.0,   0.95,  0.95,  0.65,   0.15 ]
	green = [ 1.0,  0.85, 0.7, 0.5,  0.3 ]
	blue =  [ 0.95, 0.6 , 0.3,  0.45, 0.65 ]
    elif name == "blue":
	stops = [ 0.0, 0.5, 1.0]
	red =   [ 0.0, 1.0, 0.0]
	green = [ 0.0, 0.0, 0.0]
	blue =  [ 1.0, 0.0, 0.0]
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]
        #from PlotsSMS/python/smsPlotXSEC.py
        stops = [0.00, 0.20, 0.70, 0.90, 1.00]
        red   = [0.00, 0.00, 1.00, 1.00, 1.00]
        green = [0.00, 1.00, 1.00, 0.30, 0.00]
        blue  = [1.00, 1.00, 0.00, 0.20, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)
    
    npoints = len(s)
    rt.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    rt.gStyle.SetNumberContours(ncontours)


if __name__ == '__main__':
    
    rt.gROOT.SetBatch()
    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Input/Output directory to store output")
    parser.add_option('--toys',dest="doHybridNew",default=False,action='store_true',
                  help="for toys instead of asymptotic")
    parser.add_option('--xsec-file',dest="refXsecFile",default="./data/gluino13TeV.txt",type="string",
                  help="Input directory")
    parser.add_option('--no-smooth', dest='noSmooth', action='store_true', 
                  help='Draw grid without interpolation')
    parser.add_option('--degen', type=int, help='Degeneracy of SUSY states (for T2qq limits)', default=1)
    
    (options,args) = parser.parse_args()
    
    box = options.box
    model = options.model
    directory = options.outDir

    
    refXsecFile = options.refXsecFile
    doHybridNew = options.doHybridNew
                
    set_palette("rainbow",255)
    rt.gStyle.SetOptStat(0)
    rt.gROOT.ProcessLine(".L macros/swissCrossInterpolate.h+")
    rt.gSystem.Load("macros/swissCrossInterpolate_h.so")

    try:
        sms = sms_models[model]
    except KeyError:
        sys.exit("Model {} is not implemented!".format(model))
    mgMin = sms.mgMin
    mgMax = sms.mgMax
    mchiMin = sms.mchiMin
    mchiMax = sms.mchiMax 
    binWidth = sms.binWidth
    nRebins = sms.nRebins
    xsecMin = sms.xsecMin
    xsecMax = sms.xsecMax
    diagonalOffset = sms.diagonalOffset
    smoothing = sms.smoothing
    epsilon = sms.epsilon
    fixLSP0 = sms.fixLSP0

    if model=="T1bri":
        xsecFile = rt.TFile.Open("%s/smoothXsecUL_%s.root"%(directory,box))
        xsecTree = xsecFile.Get("smoothXsecTree")
    else: 
        if doHybridNew:
            xsecFile = rt.TFile.Open("%s/xsecUL_HybridNew_%s.root"%(directory,box))
        else: 
            xsecFile = rt.TFile.Open("%s/xsecUL_Asymptotic_%s.root"%(directory,box))
        xsecTree = xsecFile.Get("xsecTree")
    xsecGluino =  rt.TH2D("xsecGluino","xsecGluino",int((mgMax-mgMin)/binWidth),mgMin, mgMax,int((mchiMax-mchiMin)/binWidth), mchiMin, mchiMax)
    xsecGluinoPlus =  rt.TH2D("xsecGluinoPlus","xsecGluinoPlus",int((mgMax-mgMin)/binWidth),mgMin, mgMax,int((mchiMax-mchiMin)/binWidth), mchiMin, mchiMax)
    xsecGluinoMinus =  rt.TH2D("xsecGluinoMinus","xsecGluinoMinus",int((mgMax-mgMin)/binWidth),mgMin, mgMax,int((mchiMax-mchiMin)/binWidth), mchiMin, mchiMax)
    
    clsTypes = ["Exp","ExpMinus","ExpMinus2","ExpPlus","ExpPlus2","Obs","ObsMinus","ObsPlus"]

    smooth = {}
    for clsType in clsTypes:
        if 'Obs' in clsType:
            smooth[clsType] = smoothing
        else:
            smooth[clsType] = 50
    
    titleMap = {"Exp":"Expected","ExpMinus":"Expected-1#sigma","ExpPlus":"Expected+1#sigma",
                "ExpMinus2":"Expected-2#sigma","ExpPlus2":"Expected+2#sigma",
                "ObsMinus":"Observed-1#sigma", "ObsPlus":"Observed+1#sigma","Obs":"Observed"}
    whichCLsVar = {"Obs":"xsecULObs_%s"%(box),"ObsPlus":"xsecULObs_%s"%(box),"ObsMinus":"xsecULObs_%s"%(box),
                "Exp":"xsecULExp_%s"%(box),"ExpPlus":"xsecULExpMinus_%s"%(box),"ExpMinus":"xsecULExpPlus_%s"%(box),
                "ExpPlus2":"xsecULExpMinus2_%s"%(box),"ExpMinus2":"xsecULExpPlus2_%s"%(box)}
                   
    xsecUL = {}
    logXsecUL = {}
    rebinXsecUL = {}
    subXsecUL = {}
    contourFinal = {}

    
    subboxes = box.split("_")
    
    thyXsec = {}
    thyXsecErr = {}
    if refXsecFile is not None:
        print "INFO: Input ref xsec file!"
        for mg in range(int(mgMin+12.5),int(mgMax-12.5)+25,25):            
            for mchi in range(int(mchiMin+12.5),mg,25):
                for line in open(refXsecFile,'r'):
                    line = line.replace('\n','')
                    if str(mg)==line.split(',')[0]:
                        thyXsec[mg,mchi] = options.degen * float(line.split(',')[1]) #pb
                        thyXsecErr[mg,mchi] = 0.01*float(line.split(',')[2])               
    else: 
        print "ERROR: no xsec file; exiting"
        sys.exit()  
    
    for i in xrange(1,xsecGluino.GetNbinsX()+1):
        xLow = xsecGluino.GetXaxis().GetBinCenter(i)
        for j in xrange(1,xsecGluino.GetNbinsY()+1):
            yLow = xsecGluino.GetYaxis().GetBinCenter(j)
            if xLow >= yLow+diagonalOffset and xLow <= mgMax-binWidth/2:
                xsecVal = thyXsec[int(xLow),int(yLow)]
                xsecErr =  thyXsecErr[int(xLow),int(yLow)]
                xsecGluino.SetBinContent(i,j,xsecVal)
                xsecGluinoPlus.SetBinContent(i,j,xsecVal*(1+xsecErr))
                xsecGluinoMinus.SetBinContent(i,j,xsecVal*(1-xsecErr))
                
    # now rebin xsecGluino the correct number of times
    for i in xrange(0,nRebins):
        xsecGluino = rt.swissCrossRebin(xsecGluino,"NE")
        xsecGluinoPlus = rt.swissCrossRebin(xsecGluinoPlus,"NE")
        xsecGluinoMinus = rt.swissCrossRebin(xsecGluinoMinus,"NE")
                
    xyPairExp = {}
    for clsType in clsTypes:
        xsecUL[clsType] = rt.TH2D("xsecUL_%s"%clsType,"xsecUL_%s"%clsType,int((mgMax-mgMin)/binWidth),mgMin, mgMax,int((mchiMax-mchiMin)/binWidth), mchiMin, mchiMax)
        xsecTree.Project("xsecUL_%s"%clsType,"mchi:mg",whichCLsVar[clsType])
        if model=="T1bri":
            brValues = [(0.00, 1.00), (0.25, 0.25), (0.50, 0.00), (0.00, 0.50), (0.00, 0.00), (0.25, 0.50), (0.50, 0.50), (0.50, 0.25)]
            tempXsecUL = {}
            for (x, y) in brValues:
                brString = ('x%.2fy%.2f'%(x,y)).replace('.','p')
                tempXsecUL[(x,y)] = rt.TH2D("xsecUL_%s_%s"%(clsType,brString),"xsecUL_%s_%s"%(clsType,brString),int((mgMax-mgMin)/binWidth),mgMin, mgMax,int((mchiMax-mchiMin)/binWidth), mchiMin, mchiMax)
                xsecTree.Project("xsecUL_%s_%s"%(clsType,brString),"mchi:mg","%s*(x==%.2f && y==%.2f)"%(whichCLsVar[clsType],x,y))
            for iBin in range(1,xsecUL[clsType].GetNbinsX()+1):
                for jBin in range(1,xsecUL[clsType].GetNbinsY()+1):
                    allValues = {}
                    for (x, y) in brValues:
                        allValues[(x,y)] = tempXsecUL[(x,y)].GetBinContent(iBin,jBin)
                    if clsType=="Exp":
                        xyPairExp[(iBin,jBin)] = max(allValues.iteritems(), key=operator.itemgetter(1))[0]
                    xsecUL[clsType].SetBinContent(iBin,jBin,allValues[xyPairExp[(iBin,jBin)]])


        fix_hist_byhand(xsecUL[clsType],model,box,clsType)
        
            
        print "INFO: doing interpolation for %s"%(clsType)
        
        # do swiss cross average in real domain
        rebinXsecUL[clsType] = rt.swissCrossInterpolate(xsecUL[clsType],"NE")

        # do scipy multi-quadratic interpolation in log domain
        if 'Obs' in clsType:
            this_eps = epsilon
        else:
            this_eps = 5
        rebinXsecUL[clsType] = interpolate2D(rebinXsecUL[clsType],epsilon=this_eps,smooth=smooth[clsType],diagonalOffset=diagonalOffset,fixLSP0=fixLSP0)

        # do swiss cross rebin + average in real domain (should be log??)
        for i in xrange(0,nRebins):
            rebinXsecUL[clsType] = rt.swissCrossRebin(rebinXsecUL[clsType],"NE")

        # only for display purposes of underlying heat map: do swiss cross average then scipy interpolation 
        if not options.noSmooth:
            xsecUL[clsType] = rt.swissCrossInterpolate(xsecUL[clsType],"NE")
            xsecUL[clsType] = interpolate2D(xsecUL[clsType], epsilon=this_eps,smooth=smooth[clsType],diagonalOffset=diagonalOffset,fixLSP0=fixLSP0)

        # fix axes
        xsecUL[clsType].GetXaxis().SetRangeUser(xsecUL[clsType].GetXaxis().GetBinCenter(1),xsecUL[clsType].GetXaxis().GetBinCenter(xsecUL[clsType].GetNbinsX()))
        xsecUL[clsType].GetYaxis().SetRangeUser(xsecUL[clsType].GetYaxis().GetBinCenter(1),xsecUL[clsType].GetYaxis().GetBinCenter(xsecUL[clsType].GetNbinsY()))
        
    c = rt.TCanvas("c","c",500,500)
    
    for clsType in clsTypes:
        xsecUL[clsType].SetName("xsecUL_%s_%s_%s"%(clsType,model,box))
        xsecUL[clsType].SetTitle("%s %s, %s #sigma #times Branching Fraction"%(model,box.replace("_","+"),titleMap[clsType]))
        xsecUL[clsType].SetTitleOffset(1.5)
        xsecUL[clsType].SetXTitle("gluino Mass [GeV]")
        xsecUL[clsType].GetXaxis().SetTitleOffset(0.95)
        xsecUL[clsType].SetYTitle("LSP Mass [GeV]")
        xsecUL[clsType].GetYaxis().SetTitleOffset(1.25)

        # subtract the predicted xsec
        subXsecUL[clsType] = rebinXsecUL[clsType].Clone()
        if clsType=="ObsMinus":
            subXsecUL[clsType].Add(xsecGluinoMinus,-1)
        elif clsType=="ObsPlus":
            subXsecUL[clsType].Add(xsecGluinoPlus,-1)
        else:
            subXsecUL[clsType].Add(xsecGluino,-1)

        contours = array('d',[0.0])
        subXsecUL[clsType].SetContour(1,contours)
        
        xsecUL[clsType].SetMinimum(xsecMin)
        xsecUL[clsType].SetMaximum(xsecMax)
        rebinXsecUL[clsType].SetMinimum(xsecMin)
        rebinXsecUL[clsType].SetMaximum(xsecMax)
        subXsecUL[clsType].SetMaximum(1.)
        subXsecUL[clsType].SetMinimum(-1.)
        
        c.SetLogz(0)
        subXsecUL[clsType].Draw("CONT Z LIST")
        c.Update()
        
        conts = rt.gROOT.GetListOfSpecials().FindObject("contours")

        xsecUL[clsType].Draw("COLZ")
        
        contour0 = conts.At(0)
        curv = contour0.First()
        finalcurv = rt.TGraph(1)
        try:
            curv.SetLineWidth(3)
            curv.SetLineColor(rt.kBlack)
            curv.Draw("lsame")
            finalcurv = curv.Clone()
            maxN = curv.GetN()
        except AttributeError:
            print "ERROR: no curve drawn for box=%s, clsType=%s -- no limit "%(box, clsType)
        
        for i in xrange(1, contour0.GetSize()):
            curv = contour0.After(curv)
            curv.SetLineWidth(3)
            curv.SetLineColor(rt.kBlack)
            curv.Draw("lsame")
            if curv.GetN()>maxN:
                maxN = curv.GetN()
                finalcurv = curv.Clone()

        contourFinal[clsType] = finalcurv

        contourFinal[clsType].SetName("%s_%s_%s"%(clsType,model,box))
        
        c.SetLogz(1)
        c.Print("%s/%s_INTERP_%s_%s.pdf"%(directory,model,box,clsType))
    outName = "%s/%s_%s_results.root"%(directory,model,box)
    if options.degen != 1:
        outName = outName.replace('.root', '_EXTRA.root')
    outFile = rt.TFile.Open(outName,"recreate")
    for clsType in clsTypes:
        contourFinal[clsType].Write()
        xsecUL[clsType].Write()

    smoothOutFile = rt.TFile.Open("%s/smoothXsecUL_%s_%s.root"%(directory,model,box), "recreate")
    
    smoothXsecTree = rt.TTree("smoothXsecTree", "smoothXsecTree")
    myStructCmd = "struct MyStruct{Double_t mg;Double_t mchi;Double_t x;Double_t y;"
    ixsecUL = 0
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+0)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+1)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+2)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+3)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+4)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+5)
    ixsecUL+=6
    myStructCmd += "}"
    rt.gROOT.ProcessLine(myStructCmd)
    from ROOT import MyStruct

    s = MyStruct()
    smoothXsecTree.Branch("mg", rt.AddressOf(s,"mg"),'mg/D')
    smoothXsecTree.Branch("mchi", rt.AddressOf(s,"mchi"),'mchi/D')
    smoothXsecTree.Branch("x", rt.AddressOf(s,"x"),'x/D')
    smoothXsecTree.Branch("y", rt.AddressOf(s,"y"),'y/D')
    if 'T1x' in model:
        s.x = float(model[model.find('x')+1:model.find('y')].replace('p','.'))
        s.y = float(model[model.find('y')+1:].replace('p','.'))
    elif model == 'T1bbbb':
        s.x = 1
        s.y = 0
    elif model == 'T1tttt':
        s.x = 0
        s.y = 1
    else:
        s.x = -1
        s.y = -1
    

    ixsecUL = 0
    smoothXsecTree.Branch("xsecULObs_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+0)),'xsecUL%i/D'%(ixsecUL+0))
    smoothXsecTree.Branch("xsecULExpPlus2_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+1)),'xsecUL%i/D'%(ixsecUL+1))
    smoothXsecTree.Branch("xsecULExpPlus_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+2)),'xsecUL%i/D'%(ixsecUL+2))
    smoothXsecTree.Branch("xsecULExp_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+3)),'xsecUL%i/D'%(ixsecUL+3))
    smoothXsecTree.Branch("xsecULExpMinus_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+4)),'xsecUL%i/D'%(ixsecUL+4))
    smoothXsecTree.Branch("xsecULExpMinus2_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+5)),'xsecUL%i/D'%(ixsecUL+5))

    for mg in range(int(mgMin-12.5),int(mgMax+binWidth-12.5),int(binWidth)):
        for mchi in range(int(mchiMin-12.5),int(mchiMax+binWidth-12.5),int(binWidth)):
            s.mg = mg
            s.mchi = mchi
            
            xsecULObs = xsecUL["Obs"].GetBinContent(xsecUL["Obs"].FindBin(mg,mchi))
            xsecULExp = xsecUL["Exp"].GetBinContent(xsecUL["Exp"].FindBin(mg,mchi))
            xsecULExpPlus = xsecUL["ExpMinus"].GetBinContent(xsecUL["ExpMinus"].FindBin(mg,mchi))
            xsecULExpMinus = xsecUL["ExpPlus"].GetBinContent(xsecUL["ExpPlus"].FindBin(mg,mchi))
            xsecULExpPlus2 = xsecUL["ExpMinus2"].GetBinContent(xsecUL["ExpMinus2"].FindBin(mg,mchi))
            xsecULExpMinus2 = xsecUL["ExpPlus2"].GetBinContent(xsecUL["ExpPlus2"].FindBin(mg,mchi))
            
            exec 's.xsecUL%i = xsecULObs'%(ixsecUL+0)
            exec 's.xsecUL%i = xsecULExpPlus2'%(ixsecUL+1)
            exec 's.xsecUL%i = xsecULExpPlus'%(ixsecUL+2)
            exec 's.xsecUL%i = xsecULExp'%(ixsecUL+3)
            exec 's.xsecUL%i = xsecULExpMinus'%(ixsecUL+4)
            exec 's.xsecUL%i = xsecULExpMinus2'%(ixsecUL+5)

            if xsecULObs > 0.:
                smoothXsecTree.Fill()

    smoothOutFile.cd()
    smoothXsecTree.Write()
    smoothOutFile.Close()

    print len(gchipairs(model)), "total points"
    print len(toFix), "points to fix"
    print toFix
