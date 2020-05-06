import ROOT as rt
import sys
import rootTools
import glob
from math import *
import os
from array import *
import random
from optparse import OptionParser
from framework import Config
from PlotFit import *
import numpy as np

strategy = 2
tol = 1e-5
maxFuncCalls = 100000
maxIter = 100000
npoints = 50

def reset(w,fr):
    for p in rootTools.RootIterator.RootIterator(fr.floatParsFinal()):
        w.var(p.GetName()).setVal(p.getVal())
        w.var(p.GetName()).setError(p.getError())
    return True

def prepareFrame(c,rFrame,options,type='data'):
    tlines = []
    cl = 0.95
    crossing = rt.TMath.Power(rt.Math.normal_quantile(1-0.5*(1-cl), 1.0),2)
    tline = rt.TLine(0,crossing,options.rMax,crossing)
    tline.SetLineColor(rt.kRed)
    tline.SetLineWidth(2)
    tlines.append(tline)
    for tline in tlines:
        rFrame.addObject(tline,"")

        
    rFrame.Draw()
    rFrame.SetMinimum(0)
    rFrame.SetMaximum(6)
    
    rFrame.SetXTitle("#mu (signal strength)")
    rFrame.SetYTitle("-2 #Delta log L(%s)"%type)
    rFrame.SetTitleSize(0.04,"X")
    rFrame.SetTitleOffset(0.85,"X")
    rFrame.SetTitleSize(0.04,"Y")
    rFrame.SetTitleOffset(0.8,"Y")
    rFrame.SetLabelSize(0.04,"X")
    rFrame.SetLabelSize(0.04,"Y")
    rFrame.SetNdivisions(505,"X")
    
    leg = rt.TLegend(0.7,0.15,0.89,0.3)
    leg.SetTextFont(42)
    leg.SetFillColor(rt.kWhite)
    leg.SetLineColor(rt.kWhite)
    leg.AddEntry("p2ll_%s"%type, "stat + syst","l")
    leg.AddEntry("n2ll_%s"%type, "stat only","l")
    leg.Draw()
    
    l = rt.TLatex()
    l.SetTextAlign(12)
    l.SetTextSize(0.04)
    l.SetTextFont(42)
    l.SetNDC()
    l.SetTextSize(0.04)
    if options.isData:
        l.DrawLatex(0.12,0.85,"CMS preliminary, %.1f fb^{-1} (13 TeV)"%(options.lumi/1000))
    else:
        l.DrawLatex(0.12,0.85,"CMS simulation, %.1f fb^{-1} (13 TeV)"%(options.lumi/1000))
    l.DrawLatex(0.12,0.80, "razor %s"%options.box.replace("_","+"))    
    if options.mGluino!=-1:
        for line in open('data/gluino13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(options.mGluino))==line.split(',')[0]:
                thyXsec = float(line.split(',')[1])*1000 # convert pb to fb
    if options.mStop!=-1:
        for line in open('data/stop13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(options.mStop))==line.split(',')[0]:
                thyXsec = float(line.split(',')[1])*1000 # convert pb to fb
    if options.mGluino!=-1:
        if options.model=="T1tttt":
            l.DrawLatex(0.12,0.75, "pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow t#bar{t}#tilde{#chi}^{0}, #sigma = %.1f fb"%thyXsec)
        if options.model=="T1bbbb":
            l.DrawLatex(0.12,0.75, "pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow b#bar{b}#tilde{#chi}^{0}, #sigma = %.1f fb"%thyXsec)
        if options.model=="T1qqqq":
            l.DrawLatex(0.12,0.75, "pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow q#bar{q}#tilde{#chi}^{0}, #sigma = %.1f fb"%thyXsec)
        l.DrawLatex(0.12,0.69,"m_{#tilde{g}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(options.mGluino,options.mLSP))
    if options.mStop!=-1:
        if options.model=="T2tt":
            l.DrawLatex(0.12,0.75, "pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}, #sigma = %.1f fb"%thyXsec)
            l.DrawLatex(0.12,0.69,"m_{#tilde{t}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(options.mStop,options.mLSP))
        if options.model=="T2bb":
            l.DrawLatex(0.12,0.75, "pp#rightarrow#tilde{b}#tilde{b}, #tilde{b}#rightarrow b#tilde{#chi}^{0}, #sigma = %.1f fb"%thyXsec)
            l.DrawLatex(0.12,0.69,"m_{#tilde{b}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(options.mStop,options.mLSP))
        if options.model=="T2qq":
            l.DrawLatex(0.12,0.75, "pp#rightarrow#tilde{q}#tilde{q}, #tilde{q}#rightarrow q#tilde{#chi}^{0}, #sigma = %.1f fb"%thyXsec)
            l.DrawLatex(0.12,0.69,"m_{#tilde{q}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(options.mStop,options.mLSP))

    c.Print(options.outDir+"/deltaLL_%s_%s.pdf"%(type,options.box))
    c.Print(options.outDir+"/deltaLL_%s_%s.C"%(type,options.box))


def testWorkspace(w,outFile,box,options):
    boxes = box.split("_")
    
     
    rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.FATAL)
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0)
    
    CMS_th1x = w.var("th1x")
    CMS_channel = w.cat("CMS_channel")
    
    CMS_set = rt.RooArgSet()
    CMS_set.add(CMS_channel)
    CMS_set.add(CMS_th1x)
    
    data_obs = w.data("data_obs")


    bPdfs = ["%s=pdf_bin%s_bonly"%(singleBox,singleBox) for singleBox in boxes]
    sbPdfs = ["%s=pdf_bin%s"%(singleBox,singleBox) for singleBox in boxes]
    
    w.factory("SIMUL::model_b2(CMS_channel,%s)"%(','.join(bPdfs)))
    w.factory("SIMUL::model_s2(CMS_channel,%s)"%(','.join(sbPdfs)))
    w.Print("v")
    model_b = w.pdf("model_b")
    #model_b = w.pdf("pdf_bin%s_bonly"%singleBox)
    model_sb = w.pdf("model_s")
    #model_s = w.pdf("pdf_bin%s"%singleBox)


    r = w.var('r')
    r.setMin(0.)
    r.setMax(options.rMax)
    poi = w.set('POI')
    
    allParams = model_sb.getParameters(data_obs)
    rt.RooStats.RemoveConstantParameters(allParams)

    opt = rt.RooLinkedList()
    opt.Add(rt.RooFit.CloneData(False))
    opt.Add(rt.RooFit.Constrain(allParams))

    nll_b = model_b.createNLL(data_obs,opt)
    nll_sb = model_sb.createNLL(data_obs,opt)

    r.setVal(0.)
    r.setConstant(True)
    minim_b = rt.RooMinimizer(nll_b)
    minim_b.setPrintLevel(-1)
    minim_b.setPrintEvalErrors(-1)
    minim_b.setStrategy(strategy)
    minim_b.setEps(tol)
    minim_b.optimizeConst(2)    
    minim_b.setMaxFunctionCalls(maxFuncCalls)
    minim_b.setMaxIterations(maxIter)
    status_b = minim_b.minimize('Minuit2','migrad')
    status_b = minim_b.minimize('Minuit2','improve')
    fr_b = minim_b.save()
    fr_b.Print("v")

    r.setVal(options.rAsimov)
    asimov_b = model_sb.generateBinned(CMS_set,rt.RooFit.Asimov())
    asimov_b.Print("V")
    
    r.setConstant(False)
    minim_sb = rt.RooMinimizer(nll_sb)
    minim_sb.setPrintLevel(-1)
    minim_sb.setPrintEvalErrors(-1)
    minim_sb.setStrategy(strategy)
    minim_sb.setEps(tol)
    minim_sb.optimizeConst(2)
    minim_sb.setMaxFunctionCalls(maxFuncCalls)
    minim_sb.setMaxIterations(maxIter)
    status_sb = minim_sb.minimize('Minuit2','migrad')
    status_sb = minim_sb.minimize('Minuit2','improve')
    fr_sb = minim_sb.save()
    fr_sb.Print("v")

    asimov_sb = model_sb.generateBinned(CMS_set,rt.RooFit.Asimov())
    
    minNll_b = fr_b.minNll()
    minNll_sb = fr_sb.minNll()

    rBestFit = r.getVal()
            
    cfg = Config.Config(options.config)

    c = rt.TCanvas('d','d',500,400)
    
    for singleBox in boxes:
        CMS_channel.setLabel(singleBox)
        x = array('d', cfg.getBinning(singleBox)[0]) # MR binning
        y = array('d', cfg.getBinning(singleBox)[1]) # Rsq binning
        z = array('d', cfg.getBinning(singleBox)[2]) # nBtag binning
        nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)


        if w.pdf('shapeSig_%s_%s_%s_morph'%(singleBox,singleBox,options.model)) != None:
            # naming convention with shape systematics
            w.factory('SUM::pdf_bin%s_sonly(n_exp_final_bin%s_proc_%s_%s*shapeSig_%s_%s_%s_morph)'%(singleBox,singleBox,singleBox,options.model,singleBox,singleBox,options.model))
        else:
            # naming convention without shape systematics
            w.factory('SUM::pdf_bin%s_sonly(n_exp_bin%s_proc_%s_%s*shapeSig_%s_%s_%sPdf)'%(singleBox,singleBox,singleBox,options.model,singleBox,options.model,singleBox))
        model_s = w.pdf('pdf_bin%s_sonly'%singleBox)
        asimov_s = model_s.generateBinned(rt.RooArgSet(CMS_th1x),rt.RooFit.Asimov())
        h_sig_th1x = asimov_s.createHistogram('h_sig_th1_%s'%singleBox,CMS_th1x)
        
        data_obs_red = data_obs.reduce(rt.RooFit.Cut("CMS_channel==CMS_channel::%s"%singleBox))
        h_data_th1x = data_obs_red.createHistogram('h_data_th1x_%s'%singleBox,CMS_th1x)
        
        asimov_sb_red = asimov_sb.reduce(rt.RooFit.Cut("CMS_channel==CMS_channel::%s"%singleBox))
        h_sigbkgd_th1x = asimov_sb_red.createHistogram('h_sigbkgd_th1_%s'%singleBox,CMS_th1x)

        asimov_b_red = asimov_b.reduce(rt.RooFit.Cut("CMS_channel==CMS_channel::%s"%singleBox))
        h_bkgd_th1x = asimov_b_red.createHistogram('h_bkgd_th1x_%s'%singleBox,CMS_th1x)
        
        h_data_th1x.SetLineColor(rt.kBlack)
        h_data_th1x.SetMarkerStyle(20)
        h_sigbkgd_th1x.SetLineColor(rt.kBlue)
        h_sigbkgd_th1x.SetLineWidth(2)
        h_sig_th1x.SetLineColor(rt.kRed)
        h_sig_th1x.SetLineWidth(2)
        h_bkgd_th1x.SetLineColor(rt.kGreen)
        h_bkgd_th1x.SetLineWidth(2)
        
        #c.SetLogy(1)
        #frame = CMS_th1x.frame(rt.RooFit.Bins(nBins),rt.RooFit.Range(0,nBins),rt.RooFit.Title("th1x frame"))
        #data_obs.plotOn(frame,rt.RooFit.LineColor(rt.kBlack),rt.RooFit.Name("data"),rt.RooFit.Cut("CMS_channel==CMS_channel::%s"%singleBox),rt.RooFit.Invisible())
        ##model_sb.plotOn(frame,rt.RooFit.LineColor(rt.kBlue),rt.RooFit.Name("sigbkgd"),rt.RooFit.Slice(CMS_channel,singleBox),rt.RooFit.ProjWData(rt.RooArgSet(CMS_channel),data_obs))
        ##model_s.plotOn(frame,rt.RooFit.LineColor(rt.kRed),rt.RooFit.Name("sig"),rt.RooFit.ProjWData(rt.RooArgSet(CMS_channel),data_obs))
        ##model_b.plotOn(frame,rt.RooFit.LineColor(rt.kGreen),rt.RooFit.Name("bkgd"),rt.RooFit.Slice(CMS_channel,singleBox),rt.RooFit.ProjWData(rt.RooArgSet(CMS_channel),data_obs))
        #frame.addTH1(h_data_th1x,'pe')
        #frame.addTH1(h_sigbkgd_th1x,'hist')
        ##frame.addTH1(h_bkgd_th1x,'hist')
        #frame.SetMinimum(1e-1)
        #frame.Draw()
        #c.Print(options.outDir+"/c_%s.pdf"%singleBox)
        #c.Print(options.outDir+"/c_%s.C"%singleBox)
        #c.SetLogy(0)                
    
        h_data_nBtagRsqMR = get3DHistoFrom1D(h_data_th1x,x,y,z,"h_data_nBtagRsqMR_%s"%singleBox)
        h_bkgd_nBtagRsqMR = get3DHistoFrom1D(h_bkgd_th1x,x,y,z,"h_bkgd_nBtagRsqMR_%s"%singleBox)
        h_sig_nBtagRsqMR = get3DHistoFrom1D(h_sig_th1x,x,y,z,"h_sig_nBtagRsqMR_%s"%singleBox)
        h_sigbkgd_nBtagRsqMR = get3DHistoFrom1D(h_sigbkgd_th1x,x,y,z,"h_sigbkgd_nBtagRsqMR_%s"%singleBox)
    
        h_data_MR = h_data_nBtagRsqMR.Project3D("xe")
        h_data_Rsq = h_data_nBtagRsqMR.Project3D("ye")
        h_bkgd_MR = h_bkgd_nBtagRsqMR.Project3D("xe")
        h_bkgd_Rsq = h_bkgd_nBtagRsqMR.Project3D("ye")
        h_sig_MR = h_sig_nBtagRsqMR.Project3D("xe")
        h_sig_Rsq = h_sig_nBtagRsqMR.Project3D("ye")
        h_sigbkgd_MR = h_sigbkgd_nBtagRsqMR.Project3D("xe")
        h_sigbkgd_Rsq = h_sigbkgd_nBtagRsqMR.Project3D("ye")

        h_colors = [rt.kRed]
        h_th1x_components = [h_sig_th1x]
        h_labels = ['signal']
    
        h_colors = [rt.kRed]
        h_MR_components = [h_sig_MR]
        h_labels = ['signal']    
    
        h_colors = [rt.kRed]
        h_Rsq_components = [h_sig_Rsq]
        h_labels = ['signal']

        btagLabel = "#geq %i b-tag" % z[0]
        lumiLabel = "%.0f pb^{-1} (13 TeV)" % (options.lumi)
        boxLabel = "razor %s %s" % (singleBox,btagLabel)
        dataString = "Data"
        plotLabel = "Full Projection"
    
        print1DProj(c,outFile,h_sigbkgd_th1x,h_data_th1x,options.outDir+"/h_th1x_%s.pdf"%singleBox,"Bin Number","Events",lumiLabel,boxLabel,plotLabel,True,False,-1,None,h_th1x_components,h_colors,h_labels)
        print1DProj(c,outFile,h_sigbkgd_MR,h_data_MR,options.outDir+"/h_MR_%s.pdf"%singleBox,"M_{R} [GeV]","Events",lumiLabel,boxLabel,plotLabel,True,False,-1,None,h_MR_components,h_colors,h_labels)
        print1DProj(c,outFile,h_sigbkgd_Rsq,h_data_Rsq,options.outDir+"/h_Rsq_%s.pdf"%singleBox,"R^{2}","Events",lumiLabel,boxLabel,plotLabel,True,False,-1,None,h_Rsq_components,h_colors,h_labels)

        h_sigbkgd_MR_components = []
        h_sigbkgd_Rsq_components = []
        h_sigbkgd_RsqMR_components = []
        h_sigbkgd_th1x_components = []
        h_sig_MR_components = []
        h_sig_Rsq_components = []
        h_sig_RsqMR_components = []
        h_sig_th1x_components = []
        h_data_MR_components = []
        h_data_Rsq_components = []
        h_data_RsqMR_components = []
        h_data_th1x_components = []
        h_labels = []        
        h_colors = [rt.kOrange,rt.kViolet,rt.kRed,rt.kGreen]
        if len(z)>2:
            for k in range(1,len(z)):
                h_sigbkgd_MR_components.append(h_sigbkgd_nBtagRsqMR.ProjectionX("h_sigbkgd_MR_%ibtag_%s"%(z[k-1],singleBox),0,-1,k,k,""))
                h_sigbkgd_Rsq_components.append(h_sigbkgd_nBtagRsqMR.ProjectionY("h_sigbkgd_Rsq_%ibtag_%s"%(z[k-1],singleBox),0,-1,k,k,""))
                h_sig_MR_components.append(h_sig_nBtagRsqMR.ProjectionX("h_sig_MR_%ibtag_%s"%(z[k-1],singleBox),0,-1,k,k,""))
                h_sig_Rsq_components.append(h_sig_nBtagRsqMR.ProjectionY("h_sig_Rsq_%ibtag_%s"%(z[k-1],singleBox),0,-1,k,k,""))
                h_data_MR_components.append(h_data_nBtagRsqMR.ProjectionX("h_MR_data_%ibtag_%s"%(z[k-1],singleBox),0,-1,k,k,""))
                h_data_Rsq_components.append(h_data_nBtagRsqMR.ProjectionY("h_Rsq_data_%ibtag_%s"%(z[k-1],singleBox),0,-1,k,k,""))
            
                h_sigbkgd_nBtagRsqMR.GetZaxis().SetRange(k,k)
                h_sigbkgd_RsqMR_components.append(h_sigbkgd_nBtagRsqMR.Project3D("%ibtag_yx"%z[k]))
                h_sigbkgd_th1x_components.append(get1DHistoFrom2D(h_sigbkgd_RsqMR_components[-1],x,y,'h_sigbgkd_th1x_%ibtag_%s'%(z[k-1],singleBox)))
                
                h_sig_nBtagRsqMR.GetZaxis().SetRange(k,k)
                h_sig_RsqMR_components.append(h_sig_nBtagRsqMR.Project3D("%ibtag_yx"%z[k]))
                h_sig_th1x_components.append(get1DHistoFrom2D(h_sig_RsqMR_components[-1],x,y,'h_sig_th1x_%ibtag_%s'%(z[k-1],singleBox)))            
            
                h_data_nBtagRsqMR.GetZaxis().SetRange(k,k)
                h_data_RsqMR_components.append(h_data_nBtagRsqMR.Project3D("%ibtag_yx"%z[k]))
                h_data_th1x_components.append(get1DHistoFrom2D(h_data_RsqMR_components[-1],x,y,'h_data_th1x_%ibtag_%s'%(z[k-1],singleBox)))
                
                if z[k-1]==3 and z[-1]==4:
                    h_labels.append("#geq %i b-tag" % z[k-1] )
                if z[k-1]==1 and z[-1]==4 and box in ['MuEle','MuMu','EleEle']:                
                    h_labels.append("#geq %i b-tag" % z[k-1] )
                else:            
                    h_labels.append("%i b-tag" % z[k-1] )
                
        if len(z)>2:
            for k in range(0,len(z)-1):
                newBoxLabel = "razor %s %s"%(box,h_labels[k])
                print1DProj(c,outFile,h_sigbkgd_th1x_components[k],h_data_th1x_components[k],options.outDir+"/h_th1x_%ibtag_%s.pdf"%(z[k],box),"Bin Number","Events",lumiLabel,newBoxLabel,plotLabel,True,False,-1,None,[h_sig_th1x_components[k]],[rt.kRed],['signal'])
                print1DProj(c,outFile,h_sigbkgd_MR_components[k],h_data_MR_components[k],options.outDir+"/h_MR_%ibtag_%s.pdf"%(z[k],box),"M_{R} [GeV]","Events",lumiLabel,newBoxLabel,plotLabel,True,False,-1,None,[h_sig_MR_components[k]],[rt.kRed],['signal'])
                print1DProj(c,outFile,h_sigbkgd_Rsq_components[k],h_data_Rsq_components[k],options.outDir+"/h_Rsq_%ibtag_%s.pdf"%(z[k],box),"R^{2}","Events",lumiLabel,newBoxLabel,plotLabel,True,False,-1,None,[h_sig_Rsq_components[k]],[rt.kRed],['signal'])
                
    pll = nll_sb.createProfile(poi)
    n2ll = rt.RooFormulaVar("n2ll","2*@0-2*%f"%minNll_sb,rt.RooArgList(nll_sb))
    p2ll = n2ll.createProfile(poi)
    
    print "signal+background nll = %f on data at r = %f"%(minNll_sb,rBestFit)
    
    xs, ys = np.linspace(0, options.rMax, npoints + 1), []    
    for x in xs:
        r.setVal(x)
        r.setConstant(True)
        ys.append(2.*nll_sb.getVal() - 2.*minNll_sb)
    xp, yp = np.linspace(0, options.rMax, npoints + 1), []
    for x in xp:        
        r.setVal(x)
        r.setConstant(True)
        minim_sb.minimize('Minuit2','migrad')
        minim_sb.minimize('Minuit2','improve')     
        if nll_sb.getVal() - minNll_sb < 0:
            reset(w,fr_sb)
            r.setVal(x)
            r.setConstant(True)
            minim_sb.minimize('Minuit2','migrad')
            minim_sb.minimize('Minuit2','improve')
        yp.append(2.*nll_sb.getVal() - 2.*minNll_sb)
        
    gr_s = rt.TGraph(len(xs), array('f', xs), array('f', ys))
    gr_s.SetLineStyle(2)
    gr_s.SetLineColor(rt.kBlue)
    gr_s.SetLineWidth(3)
    gr_s.SetName("n2ll_data")
    
    gr_p = rt.TGraph(len(xp), array('f', xp), array('f', yp))
    gr_p.SetLineColor(rt.kBlack)
    gr_p.SetLineWidth(3)
    gr_p.SetName("p2ll_data")
    
    rFrame = r.frame(rt.RooFit.Bins(npoints),rt.RooFit.Range(0.0,options.rMax),rt.RooFit.Title("r frame (data)"))
    rFrame.SetMinimum(0)
    rFrame.SetMaximum(6)

    #n2ll.plotOn(rFrame,rt.RooFit.ShiftToZero(),rt.RooFit.LineStyle(2),rt.RooFit.Name("n2ll_data"))
    #p2ll.plotOn(rFrame,rt.RooFit.LineColor(rt.kBlack),rt.RooFit.Name("p2ll_data"),rt.RooFit.Precision(-1))
    rFrame.addObject(gr_s, 'L')
    rFrame.addObject(gr_p, 'L')

    prepareFrame(c,rFrame,options,type='data')

    print "now doing asimov"
    # now do asimov
    
    reset(w,fr_b)
    
    r.setVal(0.)
    r.setConstant(True)
    
    nll_b_asimov = model_b.createNLL(asimov_b,opt)
    
    minim_b_asimov = rt.RooMinimizer(nll_b_asimov)
    minim_b_asimov.setStrategy(strategy)   
    minim_b_asimov.setPrintLevel(-1)
    minim_b_asimov.setPrintEvalErrors(-1)
    minim_b_asimov.setEps(tol)
    minim_b_asimov.optimizeConst(2)
    minim_b_asimov.setMaxFunctionCalls(maxFuncCalls)
    minim_b_asimov.setMaxIterations(maxIter)
    status_b_asimov = minim_b_asimov.minimize('Minuit2','migrad')
    status_b_asimov = minim_b_asimov.minimize('Minuit2','improve')
    fr_b_asimov = minim_b_asimov.save()
    fr_b_asimov.Print("v")
    
    minNll_b_asimov = fr_b_asimov.minNll()
    
    r.setVal(options.rAsimov)
    r.setConstant(False)

    nll_sb_asimov = model_sb.createNLL(asimov_b,opt)
    
    minim_sb_asimov = rt.RooMinimizer(nll_sb_asimov)
    minim_sb_asimov.setStrategy(strategy)   
    minim_sb_asimov.setPrintLevel(-1)
    minim_sb_asimov.setPrintEvalErrors(-1)
    minim_sb_asimov.setEps(tol)
    minim_sb_asimov.optimizeConst(2)
    minim_sb_asimov.setMaxFunctionCalls(maxFuncCalls)
    minim_sb_asimov.setMaxIterations(maxIter)
    status_sb_asimov = minim_sb_asimov.minimize('Minuit2','migrad')
    status_sb_asimov = minim_sb_asimov.minimize('Minuit2','improve')
    fr_sb_asimov = minim_sb_asimov.save()
    fr_sb_asimov.Print("v")
    
    minNll_sb_asimov = fr_sb_asimov.minNll()
    rBestFit_asimov = r.getVal()
    
    pll_asimov = nll_sb_asimov.createProfile(poi)
    n2ll_asimov = rt.RooFormulaVar("n2ll_asimov","2*@0-2*%f"%minNll_sb_asimov,rt.RooArgList(nll_sb_asimov))
    p2ll_asimov = n2ll_asimov.createProfile(poi)
    
    print "signal+background nll = %f on asimov at r = %f"%(minNll_sb_asimov,rBestFit_asimov)

    xs_a, ys_a = np.linspace(0, options.rMax, npoints + 1), []    
    for x in xs_a:
        r.setVal(x)
        r.setConstant(True)
        ys_a.append(2.*nll_sb_asimov.getVal() - 2.*minNll_sb_asimov)
    xp_a, yp_a = np.linspace(0, options.rMax, npoints + 1), []
    for x in xp_a:
        r.setVal(x)
        r.setConstant(True)
        minim_sb_asimov.minimize('Minuit2','migrad')
        minim_sb_asimov.minimize('Minuit2','improve')
        if nll_sb_asimov.getVal() - minNll_sb_asimov < 0:
            reset(w,fr_sb_asimov)
            r.setVal(x)
            r.setConstant(True)
            minim_sb_asimov.minimize('Minuit2','migrad')
            minim_sb_asimov.minimize('Minuit2','improve')
        yp_a.append(2.*nll_sb_asimov.getVal() - 2.*minNll_sb_asimov)
        
    gr_sa = rt.TGraph(len(xs_a), array('f', xs_a), array('f', ys_a))
    gr_sa.SetLineStyle(2)
    gr_sa.SetLineColor(rt.kBlue)
    gr_sa.SetLineWidth(3)
    gr_sa.SetName("n2ll_asimov")
    
    gr_pa = rt.TGraph(len(xp_a), array('f', xp_a), array('f', yp_a))
    gr_pa.SetLineColor(rt.kBlack)
    gr_pa.SetLineWidth(3)
    gr_pa.SetName("p2ll_asimov")
    
    rFrame_asimov = r.frame(rt.RooFit.Bins(npoints),rt.RooFit.Range(0.0,options.rMax),rt.RooFit.Title("r frame (asimov)"))
    rFrame_asimov.SetMinimum(0)
    rFrame_asimov.SetMaximum(6)
    
    #n2ll_asimov.plotOn(rFrame_asimov,rt.RooFit.ShiftToZero(),rt.RooFit.LineStyle(2),rt.RooFit.Name("n2ll_asimov"))
    #p2ll_asimov.plotOn(rFrame_asimov,rt.RooFit.LineColor(rt.kBlack),rt.RooFit.Name("p2ll_asimov"),rt.RooFit.Precision(-1))
    rFrame_asimov.addObject(gr_sa, 'L')
    rFrame_asimov.addObject(gr_pa, 'L')

    prepareFrame(c,rFrame_asimov,options,type='asimov')

    print ''
    print "recap:\n"
    
    print "background-only   nll = %e on data at   r = %e"%(minNll_b,0)
    print "signal+background nll = %e on data at   r = %e"%(minNll_sb,rBestFit)
    print "background-only   nll = %e on asimov at r = %e"%(minNll_b_asimov,0)
    print "signal+background nll = %e on asimov at r = %e"%(minNll_sb_asimov,rBestFit_asimov)
    print ''
    
    print "signal+background n2ll profile scan on data:\n"
    print np.array(zip(xp,yp))

    print "signal+background n2ll scan on asimov:\n"
    print np.array(zip(xp_a,yp_a))
    
    
    outFile.cd()
    outFile.Close()
        

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1tttt",type="string",
                  help="signal model name")
    parser.add_option('--mGluino',dest="mGluino", default=-1,type="float",
                  help="mass of gluino")
    parser.add_option('--mStop',dest="mStop", default=-1,type="float",
                  help="mass of stop")
    parser.add_option('--mLSP',dest="mLSP", default=-1,type="float",
                  help="mass of LSP")
    parser.add_option('--rMax',dest="rMax", default=3,type="float",
                  help="maximum of r (signal strength) in profile likelihood plot")
    parser.add_option('--rAsimov',dest="rAsimov", default=0,type="float",
                  help="r (signal strength) for asimov dataset generation")
    parser.add_option('-d','--outDir',dest="outDir",default="TestWorkspace/",type="string",
                  help="Output file to store results")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('--data',dest="isData", default=False,action='store_true',
                  help="changes plots for data")

    (options,args) = parser.parse_args()

     
    rt.gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so")
    
    for f in args:
        if f.lower().endswith('.root'):
            workspaceFileName = f
    workspaceFile = rt.TFile.Open(workspaceFileName,"READ")
    outFile = rt.TFile.Open(options.outDir+"/TestWorkspace.root","RECREATE")
    w = workspaceFile.Get("w")

    testWorkspace(w,outFile,options.box,options)
