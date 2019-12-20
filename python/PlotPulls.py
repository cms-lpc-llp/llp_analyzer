import ROOT as rt
import sys
import rootTools
import glob
from math import *
import os
from array import *
from optparse import OptionParser

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
    rt.gStyle.SetPadRightMargin(0.065)
    rt.gStyle.SetPadBottomMargin(0.12)
    rt.gStyle.SetPadLeftMargin(0.1)
    
    # use large Times-Roman fonts
    rt.gStyle.SetTitleFont(42,"xyz")  # set the all 3 axes title font
    rt.gStyle.SetTitleFont(42," ")    # set the pad title font
    rt.gStyle.SetTitleSize(0.05,"xyz") # set the 3 axes title size
    rt.gStyle.SetTitleSize(0.05," ")   # set the pad title size
    rt.gStyle.SetLabelFont(42,"xyz")
    rt.gStyle.SetLabelSize(0.05,"xyz")
    rt.gStyle.SetLabelColor(1,"xyz")
    rt.gStyle.SetTextFont(42)
    rt.gStyle.SetTextSize(0.08)
    rt.gStyle.SetStatFont(42)
    
    # use bold lines and markers
    #rt.gStyle.SetMarkerStyle(8)
    #rt.gStyle.SetLineStyleString(2,"[12 12]") # postscript dashes
    
    #..Get rid of X error bars
    rt.gStyle.SetErrorX(0.001)
    
    # do not display any of the standard histogram decorations
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptFit(1111)
    rt.gStyle.SetStatY(0.85);                
    rt.gStyle.SetStatX(0.92);                
    rt.gStyle.SetStatW(0.15);                
    rt.gStyle.SetStatH(0.15);                
    
    # put tick marks on top and RHS of plots
    #rt.gStyle.SetPadTickX(1)
    #rt.gStyle.SetPadTickY(1)
    
    ncontours = 999
    
    stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
    red =   [ 1.0,   0.95,  0.95,  0.65,   0.15 ]
    green = [ 1.0,  0.85, 0.7, 0.5,  0.3 ]
    blue =  [ 0.95, 0.6 , 0.3,  0.45, 0.65 ]
    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)
        
    npoints = len(s)
    rt.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    rt.gStyle.SetNumberContours(ncontours)
   
    rt.gStyle.cd()
        
def getFileName(mg, mchi, r, box,directory):
    rString = str('%.3f'%r).replace(".","p")
    fileName = "%s/toys_Freq_r%s_%s.root"%(directory,rString,box)
    print fileName
    return fileName



def getPulls(mg, mchi, r, box, directory,c):
    
    LzCut = "covQual_%s==3"%(box)
    #LzCut = "covQual_%s==3&&migrad_%s==0&&minos_%s==0"%(box,box,box)
    #LzCut = "r_error>0"
    
    fileName = getFileName(mg,mchi,r,box,directory)
            
    hypoTree = rt.TChain("myTree")
    addToChain = fileName.replace("//","/")+"/myTree"
    print "adding to chain: %s"% addToChain
    hypoTree.Add(addToChain)
    
    pullHist = rt.TH1D("pullHist","pullHist",40,-4,4)
    pullHist.SetMarkerStyle(20)
    pullHist.SetLineStyle(20)
    pullHist.SetMarkerSize(1)

    hypoTree.Project('pullHist','(r - %f)/r_error'%(r),LzCut)
    pullHist.SetMaximum(1.75*pullHist.GetBinContent(pullHist.GetMaximumBin()))
    pullHist.Draw("pe")
    gauspull = rt.TF1("gauspull","gaus",-4, 4)
    pullHist.Fit("gauspull","RL","pe")
    pullHist.GetXaxis().SetTitle("Pull (#hat{#mu}-#mu)/#delta#hat#mu")
    pullHist.GetYaxis().SetTitle("Toy Datasets")
    rString = str('%.3f'%r).replace(".","p")    
    l = rt.TLatex()
    l.SetTextAlign(12)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    
    l.DrawLatex(0.1,0.95,"CMS simulation")
    l.DrawLatex(0.65,0.95,"%.1f fb^{-1} (13 TeV)"%(options.lumi/1000.))
    if model=="T1bbbb":
        #l.DrawLatex(0.1,0.95,"m_{#tilde{g}} = %.0f GeV; m_{#tilde{#chi}} = %.0f GeV; %s Box"%(mg,mchi,box))
        l.SetTextFont(52)
        l.DrawLatex(0.12,0.85,"razor %s"%(box))
        l.SetTextFont(42)
        l.DrawLatex(0.12,0.78,"pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrowb#bar{b}#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.12,0.7,"m_{#tilde{g}} = %.0f GeV, m_{#tilde{#chi}} = %.0f GeV"%(mg,mchi))
        l.DrawLatex(0.12,0.65,"#mu = %.1f"%(r))
    c.Print("%s/pulls_%s_%i_%i_%s.pdf"%(directory,model,mg,mchi,rString))
    c.Print("%s/pulls_%s_%i_%i_%s.C"%(directory,model,mg,mchi,rString))


    rHist = rt.TH1D("rHist","rHist",40,r-min(max(r,2),3),r+min(max(r,2),3))
    rHist.SetMarkerStyle(20)
    rHist.SetLineStyle(20)
    rHist.SetMarkerSize(1)
    
    
    hypoTree.Project('rHist','r',LzCut)
    
    rHist.SetMaximum(1.75*rHist.GetBinContent(rHist.GetMaximumBin()))
    rHist.Draw("pe")
    #gausr = rt.TF1("gausr","gaus",r-0.25*max(r,1),r+0.25*max(r,1))
    #rHist.Fit("gausr","RL","pe")
    rHist.GetXaxis().SetTitle("#hat{#mu}")
    rHist.GetYaxis().SetTitle("Toy Datasets")
    l.DrawLatex(0.1,0.95,"CMS simulation")
    l.DrawLatex(0.65,0.95,"%.1f fb^{-1} (13 TeV)"%(options.lumi/1000.))
    if model=="T1bbbb":
        #l.DrawLatex(0.1,0.95,"m_{#tilde{g}} = %.0f GeV; m_{#tilde{#chi}} = %.0f GeV; %s Box"%(mg,mchi,box))
        l.SetTextFont(52)
        l.DrawLatex(0.12,0.85,"razor %s"%(box))
        l.SetTextFont(42)
        l.DrawLatex(0.12,0.78,"pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrowb#bar{b}#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.12,0.7,"m_{#tilde{g}} = %.0f GeV, m_{#tilde{#chi}} = %.0f GeV"%(mg,mchi))
        l.DrawLatex(0.12,0.65,"#mu = %.1f"%(r))
    c.Print("%s/r_%s_%i_%i_%s.pdf"%(directory,model,mg,mchi,rString))
    c.Print("%s/r_%s_%i_%i_%s.C"%(directory,model,mg,mchi,rString))


    mean = rHist.GetMean()
    sigma = rHist.GetRMS()
    
    return mean, sigma

    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store datasets")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="model name")

    (options,args) = parser.parse_args()
    
    gchipairs = [(1500,100)]
    #gchipairs = [(1000,900)]
    
    rRange = [0,0.5,1,2,5,10]
    
    model = options.model
    box = options.box
    directory = options.outDir

    
    setstyle()
    c = rt.TCanvas("c","c",500, 400)

    for mg,mchi in gchipairs:        
        rFit = []
        rErr = []
        for r in rRange:
            mean, sigma = getPulls(mg, mchi, r, box, directory,c)
            rFit.append(mean)
            rErr.append(sigma)

        emptyHisto = rt.TH1D("empty","empty",1,-1,12)        
        emptyHisto.GetXaxis().SetTitleSize(0.045)
        emptyHisto.GetXaxis().SetTitle("Injected signal strength #mu")
        emptyHisto.GetYaxis().SetTitleSize(0.045)
        emptyHisto.GetYaxis().SetTitle("Fit signal strength #hat{#mu}")
        emptyHisto.GetYaxis().SetRangeUser(-1,12)
        emptyHisto.SetLineColor(rt.kWhite)
        emptyHisto.Draw("")
        
        gr = rt.TGraphErrors(len(rRange),array('d',rRange),array('d',rFit),array('d',[0 for r in rRange]), array('d',rErr))
        gr.SetMarkerColor(rt.kBlue)
        gr.SetLineColor(rt.kBlue)
        gr.SetMarkerStyle(20)
        gr.Draw('Psame')
        #c.Update()
        #gr.GetXaxis().SetRangeUser(-1,12)
        #gr.Draw('APsame')
        c.Update()
        
        l = rt.TLatex()
        l.SetTextAlign(12)
        l.SetTextSize(0.05)
        l.SetTextFont(42)
        l.SetNDC()
        
        l.DrawLatex(0.1,0.95,"CMS simulation")
        l.DrawLatex(0.65,0.95,"%.1f fb^{-1} (13 TeV)"%(options.lumi/1000.))
        if model=="T1bbbb":
            #l.DrawLatex(0.1,0.95,"m_{#tilde{g}} = %.0f GeV; m_{#tilde{#chi}} = %.0f GeV; %s Box"%(mg,mchi,box))
            l.SetTextFont(52)
            l.DrawLatex(0.12,0.85,"razor %s"%(box))
            l.SetTextFont(42)
            l.DrawLatex(0.12,0.78,"pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrowb#bar{b}#tilde{#chi}^{0}_{1}")
            l.DrawLatex(0.12,0.7,"m_{#tilde{g}} = %.0f GeV, m_{#tilde{#chi}} = %.0f GeV "%(mg,mchi))

        l = rt.TGraph(2)
        l.SetPoint(0, -1, -1)
        l.SetPoint(1, 12, 12)
        l.SetLineWidth(2)
        l.SetLineColor(rt.kGray)
        l.Draw("lsame")
        c.Update()
        c.Print('%s/bias_%s_%i_%i.pdf'%(directory,model,mg,mchi))
        c.Print('%s/bias_%s_%i_%i.C'%(directory,model,mg,mchi))
