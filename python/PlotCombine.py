from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import re
from framework import Config

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('--mGluino',dest="mGluino", default=-1,type="float",
                  help="mass of gluino")
    parser.add_option('--mStop',dest="mStop", default=-1,type="float",
                  help="mass of stop")
    parser.add_option('--mLSP',dest="mLSP", default=-1,type="float",
                  help="mass of LSP")
    parser.add_option('--lumi-array',dest="lumi_array", default="0.2,4,10",type="string",
                  help="lumi array in fb^-1, e.g.: 0.2,4,10")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store plots")
    parser.add_option('-i','--indir',dest="inDir",default="./",type="string",
                  help="Input directory")
    parser.add_option('--signif',dest="signif",default=False,action='store_true',
                  help="plot significance instead of limit")

    (options,args) = parser.parse_args()

    signif = options.signif

    inDir = options.inDir
    lumiArray = array('d',[float(lumi) for lumi in options.lumi_array.split(',')])
    box = options.box
    model = options.model
    mGluino = options.mGluino
    mLSP = options.mLSP
    mStop = options.mStop
    
    cfg = Config.Config(options.config)    
    
    z = array('d', cfg.getBinning(box.split('_')[0])[2]) # nBtag binning
    btagMin = z[0]
    btagMax = z[-1]        
    if btagMax-1>btagMin:          
        btag = '%i-%ibtag'%(btagMin,btagMax-1)
    else:
        btag = '%ibtag'%(btagMin)    

    btagLabel = ""
    if z[-1]==4:
        btagLabel = "#geq %i b-tag" % z[0]
    elif z[-1] == z[0]+1:
        btagLabel = "%i b-tag" % z[0]
    else:
        btagLabel = "%i-%i b-tag" % (z[0],z[-1]-1)
        
    thyXsec = -1
    thyXsecErr = -1
    
    if mGluino!=-1:
        massPoint = "mGl-%i_mLSP-%i"%(mGluino,mLSP)
        for line in open('data/gluino13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(mGluino))==line.split(',')[0]:
                thyXsec = float(line.split(',')[1])*1000 # convert pb to fb
                thyXsecErr = 0.01*float(line.split(',')[2])
    if mStop!=-1:
        massPoint = "mStop-%i_mLSP-%i"%(mStop,mLSP)
        for line in open('data/stop13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(mStop))==line.split(',')[0]:
                thyXsec = float(line.split(',')[1])*1000 # convert pb to fb
                thyXsecErr = 0.01*float(line.split(',')[2])

    expArray = array('d')
    expp1sigmaArray = array('d')
    expm1sigmaArray = array('d')
    expp2sigmaArray = array('d')
    expm2sigmaArray = array('d')
    zeroArray = array('d',[0 for lumi in lumiArray])
    xsecArray = array('d',[thyXsec for lumi in lumiArray])
    xsecp1sigmaArray = array('d',[thyXsec*(thyXsecErr)/2. for lumi in lumiArray])
    xsecm1sigmaArray = array('d',[thyXsec*(thyXsecErr)/2. for lumi in lumiArray])

    sigArray = array('d')
    
    if signif:
        for lumi in lumiArray:
            tfile = rt.TFile.Open('%s/higgsCombine%s_%s_lumi-%.1f_%s_%s.ProfileLikelihood.mH120.root'%(inDir,model,massPoint,lumi,btag,box))
            
            limit = tfile.Get('limit')
            limit.Draw('>>elist','','entrylist')
            elist = rt.gDirectory.Get('elist')
            entry = elist.Next()
            limits = []
            while entry>=0:
                limit.GetEntry(entry)
                limits.append(limit.limit)
                entry = elist.Next()
            sigArray.append(limits[0])
            
        print sigArray
    else:
        
        for lumi in lumiArray:

            tfile = rt.TFile.Open('%s/higgsCombine%s_%s_lumi-%s_%s_%s.Asymptotic.mH120.root'%(inDir,model,massPoint,lumi,btag,box))
            tfile.Print('v')
            limit = tfile.Get('limit')
            limit.Draw('>>elist','','entrylist')
            elist = rt.gDirectory.Get('elist')
            entry = elist.Next()
            limits = []
            while entry>=0:
                limit.GetEntry(entry)
                limits.append(limit.limit)
                entry = elist.Next()
            expArray.append(limits[2]*thyXsec)
            expp2sigmaArray.append((limits[4]-limits[2])*thyXsec)
            expp1sigmaArray.append((limits[3]-limits[2])*thyXsec)
            expm1sigmaArray.append((limits[2]-limits[1])*thyXsec)
            expm2sigmaArray.append((limits[2]-limits[0])*thyXsec)

        print expp1sigmaArray
        print expm1sigmaArray
        print expp2sigmaArray
        print expm2sigmaArray
    
    
        expGraph = rt.TGraph(len(lumiArray),lumiArray,expArray)
        expGraph.SetLineStyle(2)
    
        exp2sigmaGraph = rt.TGraphAsymmErrors(len(lumiArray), lumiArray, expArray ,zeroArray,zeroArray, expm2sigmaArray,expp2sigmaArray)
        exp2sigmaGraph.SetLineColor(5)
        exp2sigmaGraph.SetFillColor(5)
        exp2sigmaGraph.SetFillStyle(1001)

        exp1sigmaGraph = rt.TGraphAsymmErrors(len(lumiArray), lumiArray, expArray ,zeroArray,zeroArray, expm1sigmaArray,expp1sigmaArray)
        exp1sigmaGraph.SetLineColor(rt.kGreen-7)
        exp1sigmaGraph.SetFillColor(rt.kGreen-7)
        exp1sigmaGraph.SetFillStyle(1001)

        xsecGraphErr = rt.TGraphAsymmErrors(len(lumiArray),lumiArray,xsecArray,zeroArray,zeroArray,xsecm1sigmaArray,xsecp1sigmaArray)
        xsecGraphErr.SetFillStyle(1001)
        xsecGraphErr.SetLineColor(rt.kOrange)
        xsecGraphErr.SetFillColor(rt.kBlue-7)
    
        xsecGraph = rt.TGraph(len(lumiArray),lumiArray,xsecArray)
        xsecGraph.SetMarkerSize(0)
        xsecGraph.SetLineStyle(1)
        xsecGraph.SetLineColor(rt.kOrange)

    if signif:
        i=0
        while i < len(sigArray):
            if sigArray[i]>10000:
                sigArray.pop(i)
                lumiArray.pop(i)
            i+=1

        print lumiArray
        print sigArray
        sigGraph = rt.TGraph(len(sigArray),lumiArray,sigArray)
        sigGraph.SetLineStyle(1)
        threeSigGraph = rt.TGraph(len(sigArray),lumiArray,array('d',[3. for sig in sigArray]))
        threeSigGraph.SetLineStyle(2)
        threeSigGraph.SetLineColor(rt.kRed)
        fiveSigGraph = rt.TGraph(len(sigArray),lumiArray,array('d',[5. for sig in sigArray]))
        fiveSigGraph.SetLineStyle(2)
        fiveSigGraph.SetLineColor(rt.kRed)
    
    h_limit = rt.TMultiGraph()
    c = rt.TCanvas('c','c',600,400)
    if signif:
        c.SetLogy(0)
        h_limit.Add(sigGraph)
        h_limit.Add(threeSigGraph)
        h_limit.Add(fiveSigGraph)
        h_limit.SetMaximum(10)
        h_limit.SetMinimum(0)
        h_limit.Draw("a3")
        h_limit.GetXaxis().SetLimits(lumiArray[0],lumiArray[-1])
        h_limit.GetXaxis().SetTitle("Integrated Luminosity [fb^{-1}]")
        h_limit.GetYaxis().SetTitle("Expected Signal Significance [std. dev.]")
        h_limit.Draw("a3")
        threeSigGraph.Draw("l same")
        fiveSigGraph.Draw("l same")
        sigGraph.Draw("l same")
    else:
        c.SetLogy(1)
        h_limit.Add(exp2sigmaGraph)
        h_limit.Add(exp1sigmaGraph)
        h_limit.Add(xsecGraphErr)
        h_limit.Draw("a3")
        h_limit.GetXaxis().SetLimits(lumiArray[0],lumiArray[-1])
        h_limit.GetXaxis().SetTitle("Integrated Luminosity [fb^{-1}]")
        h_limit.GetYaxis().SetTitle("95% C.L. Upper Limit Cross Section [fb]")
        #h_limit.SetMaximum(1.e3)
        #h_limit.SetMaximum(1.e4)
        h_limit.SetMaximum(100*max(xsecArray))
        h_limit.SetMinimum(0.1*min(xsecArray))
        #h_limit.SetMinimum(2.)
        #h_limit.SetMinimum(50.)
        h_limit.Draw("a3")
    
        expGraph.Draw("l same")
        xsecGraph.Draw("l same")
        demo = exp1sigmaGraph.Clone()
        demo.SetLineColor(rt.kBlack)
        demo.SetLineStyle(2)
        
    rt.gPad.Update()


    
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.045)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.15,0.84,"CMS simulation (13 TeV)")
    l.SetTextFont(52)
    boxes = box.split('_')
    if len(boxes)>1:
        l.DrawLatex(0.15,0.77,"razor %s"%'+'.join(boxes[0:2]))
        for i in range(2,len(boxes)-1):
            l.DrawLatex(0.15,0.84-0.07*i,"+%s"%boxes[i])
        if len(boxes)>2:
            l.DrawLatex(0.15,0.84-0.07*(len(boxes)-1),"+%s, %s"%(boxes[-1],btagLabel))
        else:            
            l.DrawLatex(0.15,0.84-0.07*2,"%s"%(btagLabel))
            
    else:
        l.DrawLatex(0.15,0.77,"razor %s box, %s"%(boxes[0],btagLabel))
    l.SetTextFont(42)
    if model=="T1bbbb":
        l.DrawLatex(0.52,0.84,"pp #rightarrow #tilde{g}#tilde{g},  #tilde{g}#rightarrowb#bar{b}#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.52,0.77,"m_{#tilde{g}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(mGluino,mLSP))
    elif model=="T1tttt":
        l.DrawLatex(0.52,0.84,"pp #rightarrow #tilde{g}#tilde{g},  #tilde{g}#rightarrowt#bar{t}#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.52,0.77,"m_{#tilde{g}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(mGluino,mLSP))
    elif model=="T1qqqq":
        l.DrawLatex(0.52,0.84,"pp #rightarrow #tilde{g}#tilde{g},  #tilde{g}#rightarrowq#bar{q}#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.52,0.77,"m_{#tilde{g}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(mGluino,mLSP))
    elif model=="T2tt":
        l.DrawLatex(0.52,0.84,"pp #rightarrow #tilde{t}#tilde{t},  #tilde{t}#rightarrowt#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.52,0.77,"m_{#tilde{t}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(mStop,mLSP))
    elif model=="T2bb":
        l.DrawLatex(0.52,0.84,"pp #rightarrow #tilde{b}#tilde{b},  #tilde{b}#rightarrowb#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.52,0.77,"m_{#tilde{b}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(mStop,mLSP))
    elif model=="T2qq":
        l.DrawLatex(0.52,0.84,"pp #rightarrow #tilde{q}#tilde{q},  #tilde{b}#rightarrowb#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.52,0.77,"m_{#tilde{b}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(mStop,mLSP))

    if signif:
        leg = rt.TLegend(0.52,0.64,0.89,0.74)
    else:
        leg = rt.TLegend(0.52,0.54,0.89,0.74)
    leg.SetTextFont(42)
    leg.SetFillColor(rt.kWhite)
    leg.SetLineColor(rt.kWhite)

    if signif:
        leg.AddEntry(sigGraph, "signal significance [std. dev.]", "l")
    else:
        
        if model.find("T1")!=-1:
            leg.AddEntry(xsecGraphErr, "#sigma_{NLO+NLL} (#tilde{g}#tilde{g}) #pm 1 #sigma (theory)","lf")
        elif model.find("T2tt")!=-1:
            leg.AddEntry(xsecGraphErr, "#sigma_{NLO+NLL} (#tilde{t}#tilde{t}) #pm 1 #sigma (theory)","lf")
        elif model.find("T2bb")!=-1:
            leg.AddEntry(xsecGraphErr, "#sigma_{NLO+NLL} (#tilde{b}#tilde{b}) #pm 1 #sigma (theory)","lf")
        leg.AddEntry(demo, "expected #pm 1 #sigma (experiment)","lf")

    leg.Draw()

    #rightmax = 10
    #rt.gPad.Update()
    #axis = rt.TGaxis(rt.gPad.GetUxmax(),rt.gPad.GetUymin(),rt.gPad.GetUxmax(),10*rt.gPad.GetUymax(),0,rightmax,510,"+L")
    #axis.SetLineColor(rt.kRed)
    #axis.SetLabelColor(rt.kRed)
    #axis.Draw()
    if model.find("T1")!=-1:
        if signif:
            c.Print("%s/signif_%s_%i_%i_%s_%s.pdf"%(options.outDir,model,mGluino,mLSP,btag,box))
            c.Print("%s/signif_%s_%i_%i_%s_%s.C"%(options.outDir,model,mGluino,mLSP,btag,box))
        else:
            c.Print("%s/limit_%s_%i_%i_%s_%s.pdf"%(options.outDir,model,mGluino,mLSP,btag,box))
            c.Print("%s/limit_%s_%i_%i_%s_%s.C"%(options.outDir,model,mGluino,mLSP,btag,box))
    elif model.find("T2")!=-1:
        if signif:
            c.Print("%s/signif_%s_%i_%i_%s_%s.pdf"%(options.outDir,model,mStop,mLSP,btag,box))
            c.Print("%s/signif_%s_%i_%i_%s_%s.C"%(options.outDir,model,mStop,mLSP,btag,box))
        else:
            c.Print("%s/limit_%s_%i_%i_%s_%s.pdf"%(options.outDir,model,mStop,mLSP,btag,box))
            c.Print("%s/limit_%s_%i_%i_%s_%s.C"%(options.outDir,model,mStop,mLSP,btag,box))
