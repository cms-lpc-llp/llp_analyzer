from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
from itertools import *
from operator import *
from WriteDataCard import initializeWorkspace, convertDataset2TH1
import os
import random
import sys
import math
from PlotFit import setStyle,print1DProj,print2DScatter,get3DHistoFrom1D,getBinEvents,convertSideband,densityCorr
from macro.razorFits import binnedFit

densityCorr = False
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-s','--signal',dest="signalFileName", default="None",type="string",
                  help="input dataset file for signal pdf")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model")
    parser.add_option('--mGluino',dest="mGluino", default=1500,type="float",
                  help="mgluino")
    parser.add_option('--mLSP',dest="mLSP", default=100,type="float",
                  help="mass of LSP")
    parser.add_option('-r','--signal-strength',dest="r", default=-1,type="float",
                  help="signal strength")
    parser.add_option('--no-fit',dest="noFit",default=False,action='store_true',
                  help="Turn off fit (useful for visualizing initial parameters)")
    parser.add_option('--fit-region',dest="fitRegion",default="Full",type="string",
                  help="Fit region")
    parser.add_option('--plot-region',dest="plotRegion",default="Full",type="string",
                  help="Plot region")
    parser.add_option('--data',dest="isData", default=False,action='store_true',
                  help="changes plots for data")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default=None,type="string",
                  help="input fit file")
    parser.add_option('-w','--weight',dest="useWeight",default=False,action='store_true',
                  help="use weight")
    parser.add_option('--no-plots',dest="noPlots",default=False,action='store_true',
                  help="no plots")



    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)
    
    box = options.box
    lumi = options.lumi
    noFit = options.noFit
    fitRegion = options.fitRegion
    plotRegion = options.plotRegion
    
    lumi_in = 0.

    data = None
    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            workspace = rootFile.Get('w'+box)
            data = workspace.data('RMRTree')
            lumi_in = 1000.*float([g.replace('lumi-','') for g in f.split('_') if g.find('lumi')!=-1][0])
    if data is None:
        print "give a root file as input"
  
    w = rt.RooWorkspace("w"+box)
    rootTools.Utils.importToWS(w,data)

    if options.isData:        
        paramNames, bkgs = initializeWorkspace(w,cfg,box)
    else:
        paramNames, bkgs = initializeWorkspace(w,cfg,box,lumi/lumi_in)

        
    if options.inputFitFile is not None:
        inputRootFile = rt.TFile.Open(options.inputFitFile,"r")
        wIn = inputRootFile.Get("w"+box).Clone("wIn"+box)            
        if wIn.obj("fitresult_extRazorPdf_data_obs") != None:
            frIn = wIn.obj("fitresult_extRazorPdf_data_obs")
        elif wIn.obj("nll_extRazorPdf_data_obs") != None:
            frIn = wIn.obj("nll_extRazorPdf_data_obs")
        elif wIn.obj("fitresult_extRazorPdf_data_obs_with_constr") != None:
            frIn = wIn.obj("fitresult_extRazorPdf_data_obs_with_constr")
        elif wIn.obj("nll_extRazorPdf_data_obs_with_constr") != None:
            frIn = wIn.obj("nll_extRazorPdf_data_obs_with_constr")
                        
        print "restoring parameters from fit"
        frIn.Print("V")
        for p in rootTools.RootIterator.RootIterator(frIn.floatParsFinal()):
            #if 'Ntot' in p.GetName(): continue
            w.var(p.GetName()).setVal(p.getVal())
            w.var(p.GetName()).setError(p.getError())
            
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    
    xFine = array('d', [x[0]+i*(x[-1]-x[0])/100. for i in range(0,101)]) # MR binning fine
    yFine = array('d', [y[0]+i*(y[-1]-y[0])/100. for i in range(0,101)]) # Rsq binning fine
    zFine = array('d', cfg.getBinning(box)[2]) # nBtag binning fine
    nBinsFine = (len(xFine)-1)*(len(yFine)-1)*(len(zFine)-1)
    
    th1x = w.var('th1x')
    
    sideband = convertSideband(fitRegion,w,x,y,z)
    plotband = convertSideband(plotRegion,w,x,y,z)
    
    myTH1 = convertDataset2TH1(data, cfg, box, w, options.useWeight)

    signalDs = None
    doSignalInj = (options.signalFileName != "None") and (options.r > -1)
    if doSignalInj:
        sigRootFile = rt.TFile(options.signalFileName)
        #sigWorkspace = sigRootFile.Get('w'+box)
        #signalDs = sigWorkspace.data('RMRTree')
        model = options.signalFileName.split('.root')[0].split('-')[1].split('_')[0]
        massPoint = '_'.join(options.signalFileName.split('.root')[0].split('_')[1:3])
        #sigTH1 = convertDataset2TH1(signalDs, cfg, box, w,"signal")
        #sigTH1.Scale(lumi/lumi_in)
        sigTH1 = sigRootFile.Get('%s_%s'%(box,model)).Clone('hist_%s_%s'%(box,model))
        sigDataHist = rt.RooDataHist('%s_%s'%(box,model),'%s_%s'%(box,model),rt.RooArgList(th1x), sigTH1)
        sigPdf = rt.RooHistPdf('%s_Signal'%box,'%s_Signal'%box,rt.RooArgSet(th1x), sigDataHist)
        rootTools.Utils.importToWS(w,sigDataHist)
        rootTools.Utils.importToWS(w,sigPdf)
        w.factory('r[%f]'%options.r)
        w.var('r').setConstant(False)
        w.factory('Ntot_Signal_%s_In[%f]'%(box,sigTH1.Integral()))
        w.factory('expr::Ntot_Signal_%s("r*Ntot_Signal_%s_In",r,Ntot_Signal_%s_In)'%(box,box,box))
        w.factory('SUM::extSpBPdf(Ntot_Signal_%s*%s_Signal,Ntot_TTj0b_%s*%s_TTj0b,Ntot_TTj1b_%s*%s_TTj1b,Ntot_TTj2b_%s*%s_TTj2b,Ntot_TTj3b_%s*%s_TTj3b)'%(box,box,box,box,box,box,box,box,box,box))
        w.factory('SUM::extSignalPdf(Ntot_Signal_%s*%s_Signal)'%(box,box))
        
        w.Print('v')

        
    extRazorPdf = w.pdf('extRazorPdf')
    
    if not options.isData:        
        myTH1.Scale(lumi/lumi_in)       
        if doSignalInj:
            if options.r > 0:
                #gendata = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('forsignalinj'),rt.RooFit.Asimov()) #weighted approach                
                #rt.RooRandom.randomGenerator().SetSeed(1989) # for r = 1, new                       
                rt.RooRandom.randomGenerator().SetSeed(1988) # for r = 1, new          
                #gendata = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('forsignalinj')) # real unweighted dataset
                gendata = w.pdf('extSpBPdf').generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('forsignalinj')) # real unweighted dataset
                myTH1 = gendata.createHistogram('gendata',th1x)
                #Npois = rt.RooRandom.randomGenerator().Poisson(options.r*sigTH1.Integral())       
                #myTH1.Add(sigTH1,options.r) # "weighted approach"      
                #myTH1.FillRandom(sigTH1,Npois) # "unweighted approach" - generating a real dataset
    dataHist = rt.RooDataHist("data_obs","data_obs",rt.RooArgList(th1x), rt.RooFit.Import(myTH1))
    dataHist.Print('v')
    
    
    rootTools.Utils.importToWS(w,dataHist)

    setStyle()
    
    if noFit and options.inputFitFile is not None:
        fr = frIn
        fr.Print('v')    
        rootTools.Utils.importToWS(w,fr)
    elif noFit:
        fr = rt.RooFitResult()
    else:
        fr = binnedFit(extRazorPdf,dataHist,sideband,options.useWeight,box=box)
        total = extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
        nll = 0
        iBinX = -1
        observed_counts = []
        for i in range(1,len(x)):
            for j in range(1,len(y)):
                for k in range(1,len(z)):
                    iBinX += 1
                    th1x.setVal(iBinX+0.5)
                    inSideband = any([th1x.inRange(band) for band in sideband.split(',')])
                    observed = float(dataHist.weight(rt.RooArgSet(th1x)))
                    observed_counts.append(int(observed))
                    if not inSideband: continue            
                    #expected, errorFlag = getBinEvents(i,j,k,x,y,z,w,box)
                    expected = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
                    if expected>0:
                        nll -= observed*rt.TMath.Log(expected) - expected
        print "fr      min -log(L) = ", fr.minNll()
        print "by hand min -log(L) = ", nll
            
        fr.Print('v')    
        rootTools.Utils.importToWS(w,fr)
        
        if doSignalInj:
            extSpBPdf = w.pdf('extSpBPdf')
            #frSpB = binnedFit(extSpBPdf,dataHist,sideband)
            #frSpB.Print('v')    
            #rootTools.Utils.importToWS(w,frSpB)
    
    th1x.setBins(nBins)

    asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('central'),rt.RooFit.Asimov())
        
    opt = [rt.RooFit.CutRange(myRange) for myRange in plotband.split(',')]
    asimov_reduce = asimov.reduce(opt[0])
    dataHist_reduce = dataHist.reduce(opt[0])
    for iOpt in range(1,len(opt)):
        asimov_reduce.add(asimov.reduce(opt[iOpt]))
        dataHist_reduce.add(dataHist.reduce(opt[iOpt]))
    rt.TH1D.SetDefaultSumw2()
    rt.TH2D.SetDefaultSumw2()
    rt.TH3D.SetDefaultSumw2()
    
    # start writing output
    c = rt.TCanvas('c','c',500,400)
    rootFile = rt.TFile.Open(options.outDir + '/' + 'Plots_%s'%box + '.root','recreate')
    tdirectory = rootFile.GetDirectory(options.outDir)
    print tdirectory
    if tdirectory==None:
        print "making directory"
        rootFile.mkdir(options.outDir)
        tdirectory = rootFile.GetDirectory(options.outDir)
        tdirectory.Print('v')
        
    h_th1x = asimov_reduce.createHistogram('h_th1x',th1x)
    h_data_th1x = dataHist_reduce.createHistogram('h_data_th1x',th1x)
    
    h_data_nBtagRsqMR = get3DHistoFrom1D(h_data_th1x,x,y,z,"h_data_nBtagRsqMR")
    h_nBtagRsqMR = get3DHistoFrom1D(h_th1x,x,y,z,"h_nBtagRsqMR")

    h_data_RsqMR = h_data_nBtagRsqMR.Project3D("yxe")
    h_data_MR = h_data_nBtagRsqMR.Project3D("xe")
    h_data_Rsq = h_data_nBtagRsqMR.Project3D("ye")
    h_RsqMR = h_nBtagRsqMR.Project3D("yxe")
    h_MR = h_nBtagRsqMR.Project3D("xe")
    h_Rsq = h_nBtagRsqMR.Project3D("ye")
    
    if doSignalInj:        
        extSignalPdf = w.pdf('extSignalPdf')
        signal = extSignalPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('signal'),rt.RooFit.Asimov())
        signal_reduce = signal.reduce(*opt)
        h_sig_th1x = signal_reduce.createHistogram('h_sig_th1x',th1x)
        h_sig_nBtagRsqMR = get3DHistoFrom1D(h_sig_th1x,x,y,z,"h_sig_nBtagRsqMR")
        h_sig_RsqMR = h_sig_nBtagRsqMR.Project3D("yxe")
        h_sig_MR = h_sig_nBtagRsqMR.Project3D("xe")
        h_sig_Rsq = h_sig_nBtagRsqMR.Project3D("ye")    
    
    if len(z)>1:
        h_MR_components = []
        h_Rsq_components = []
        h_labels = []        
        h_colors = [rt.kOrange,rt.kViolet,rt.kRed,rt.kGreen,rt.kGray+2]
        for k in range(1,len(z)):
            h_MR_components.append(h_nBtagRsqMR.ProjectionX("MR_%ibtag"%z[k-1],0,-1,k,k,""))
            h_Rsq_components.append(h_nBtagRsqMR.ProjectionY("Rsq_%ibtag"%z[k-1],0,-1,k,k,""))
            if z[k-1]==3 and z[-1]==4:
                h_labels.append("#geq %i b-tag" % z[k-1] )
            elif z[k-1]==1 and z[-1]==4 and len(z)==2:
                h_labels.append("#geq %i b-tag" % z[k-1] )
            else:            
                h_labels.append("%i b-tag" % z[k-1] )
                
    ## h_data_nBtagRsqMR_fine = rt.TH3D("h_data_nBtagRsqMR_fine","h_data_nBtagRsqMR_fine",len(xFine)-1,xFine,len(yFine)-1,yFine,len(zFine)-1,zFine)
    ## h_nBtagRsqMR_fine = rt.TH3D("h_nBtagRsqMR_fine","h_nBtagRsqMR_fine",len(xFine)-1,xFine,len(yFine)-1,yFine,len(zFine)-1,zFine)
    ## opt = [rt.RooFit.CutRange(myRange) for myRange in plotRegion.split(',')]
    ## data_reduce = w.data("RMRTree").reduce(opt[0])
    ## for iOpt in range(1,len(opt)):
    ##     data_reduce.append(w.data("RMRTree").reduce(opt[iOpt]))
    ## data_reduce.fillHistogram(h_data_nBtagRsqMR_fine,rt.RooArgList(w.var("MR"),w.var("Rsq"),w.var("nBtag")))
    
    ## for i in range(1,len(xFine)):
    ##     for j in range(1,len(yFine)):
    ##         for k in range(1,len(zFine)):                
    ##             w.var('MR').setVal((xFine[i]+xFine[i-1])/2.)
    ##             w.var('Rsq').setVal((yFine[j]+yFine[j-1])/2.)
    ##             w.var('nBtag').setVal((zFine[k]+zFine[k-1])/2.)
    ##             inSideband = 0
    ##             for myRange in plotRegion.split(','):
    ##                 inSideband += ( w.var('MR').inRange(myRange) * w.var('Rsq').inRange(myRange) * w.var('nBtag').inRange(myRange) )
    ##             if not inSideband: continue  
    ##             value, errorFlag = getBinEvents(i,j,k,xFine,yFine,zFine,w,box)
    ##             if not errorFlag:
    ##                 h_nBtagRsqMR_fine.SetBinContent(i,j,k,value)
            
    ## h_data_RsqMR_fine = h_data_nBtagRsqMR_fine.Project3D("yxe")
    ## h_RsqMR_fine = h_nBtagRsqMR_fine.Project3D("yxe")

    btagLabel = ""
    if z[-1] == z[0]+1 and z[-1]==4:
        btagLabel = "#geq %i b-tag" % z[0]
    elif z[-1] == z[0]+1:
        btagLabel = "%i b-tag" % z[0]
    elif z[-1]==4:
        btagLabel = "#geq %i b-tag" % z[0]
    elif len(z)==2 and z[0]==1 and z[-1]==4:
        btagLabel = "#geq %i b-tag" % z[0]        
    else:
        btagLabel = "%i-%i b-tag" % (z[0],z[-2])

    if options.isData:
        lumiLabel = "%.0f pb^{-1} (13 TeV)" % (lumi)
    else:        
        lumiLabel = "%.0f fb^{-1} (13 TeV)" % (lumi/1000)
    boxLabel = "razor %s %s %s Fit" % (box,btagLabel,fitRegion)
    plotLabel = "%s Projection" % (plotRegion)

    
    if options.isData:
        dataString = "Data"
    else:
        dataString = "Sim. Data"
        
    for h in [h_nBtagRsqMR,h_data_nBtagRsqMR,h_RsqMR,h_data_RsqMR,h_MR,h_data_MR,h_Rsq,h_data_Rsq]:
        tdirectory.cd()
        h.Write()
        
    eventsLabel = "Events"
    if densityCorr:
        eventsLabel = "Events/Bin Width"
        
    h_sig_total_th1x_components = []
    if doSignalInj:
        h_MR_components.append(h_sig_MR)
        h_Rsq_components.append(h_sig_Rsq)
        h_sig_total_th1x_components.append(h_sig_th1x)

    if not options.noPlots:
        
        print1DProj(c,tdirectory,h_th1x,h_data_th1x,options.outDir+"/h_th1x_%s.pdf"%box,"Bin Number",eventsLabel,lumiLabel,boxLabel,plotLabel,options.isData,doSignalInj,options,None,h_sig_total_th1x_components)
        print1DProj(c,tdirectory,h_MR,h_data_MR,options.outDir+"/h_MR_%s.pdf"%box,"M_{R} [GeV]",eventsLabel,lumiLabel,boxLabel,plotLabel,options.isData,doSignalInj,options,None,h_MR_components,h_colors,h_labels)
        print1DProj(c,tdirectory,h_Rsq,h_data_Rsq,options.outDir+"/h_Rsq_%s.pdf"%box,"R^{2}",eventsLabel,lumiLabel,boxLabel,plotLabel,options.isData,doSignalInj,options,None,h_Rsq_components,h_colors,h_labels)
    
        #print2DScatter(c,tdirectory,h_RsqMR_fine,options.outDir+"/h_RsqMR_scatter_log_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Fit",lumiLabel,boxLabel,plotLabel,x,y,1e-1,h_data_RsqMR_fine.GetMaximum(),options.isData)
        #print2DScatter(c,tdirectory,h_data_RsqMR_fine,options.outDir+"/h_RsqMR_scatter_data_log_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Sim. Data",lumiLabel,boxLabel,plotLabel,x,y,1e-1,h_data_RsqMR_fine.GetMaximum(),options.isData)

    outFileName = "BinnedFitResults_%s.root"%(box)
    outFile = rt.TFile.Open(options.outDir+"/"+outFileName,'recreate')
    outFile.cd()
    w.Write()
    outFile.Close()
