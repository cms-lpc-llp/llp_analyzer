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
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, FileTransferSpeed, FormatLabel, Percentage, ProgressBar, ReverseBar, RotatingMarker, SimpleProgress, Timer
        

def getTree(myTree,paramNames,nBins,box,z):
    
    rando = random.randint(1,999999)
    # first structure
    stringMyStruct1 = "void tempMacro_%d(){struct MyStruct1{"%(rando)
    
    stringMyStruct1 = stringMyStruct1+"int toy_num; int migrad_sf; int migrad_ff; int hesse_sf; int hesse_ff; int covQual_sf; int covQual_ff;"
    
    for paramName in paramNames:
        stringMyStruct1 = stringMyStruct1+'float %s_ff; float %s_ff_error; float %s_sf; float %s_sf_error;' %(paramName,paramName,paramName,paramName)
    for iBinX in range(0,nBins):
        stringMyStruct1 = stringMyStruct1+"float b%i_sf;" %(iBinX)        
    for iBinX in range(0,nBins):
        stringMyStruct1 = stringMyStruct1+"float b%i_ff;" %(iBinX)
    for iBinX in range(0,nBins):
        stringMyStruct1 = stringMyStruct1+"float b%i_toy;" %(iBinX)
    tempMacro = open("tempMacro_%d.C"%rando,"w")
    tempMacro.write(stringMyStruct1+"};}")
    tempMacro.close()
    rt.gROOT.Macro("tempMacro_%d.C"%rando)
    #rt.gROOT.ProcessLine(".x tempMacro_%d.C"%rando)
    #rt.gROOT.ProcessLine(stringMyStruct1+"};")
    from ROOT import MyStruct1

    # fills the bins list and associate the bins to the corresponding variables in the structure
    s1 = MyStruct1()
    myTree.Branch("toy_num", rt.AddressOf(s1,"toy_num"),'toy_num/I')
    myTree.Branch("migrad_sf", rt.AddressOf(s1,"migrad_sf"),'migrad_sf/I')
    myTree.Branch("migrad_ff", rt.AddressOf(s1,"migrad_ff"),'migrad_ff/I')
    myTree.Branch("hesse_sf", rt.AddressOf(s1,"hesse_sf"),'hesse_sf/I')
    myTree.Branch("hesse_ff", rt.AddressOf(s1,"hesse_ff"),'hesse_ff/I')
    myTree.Branch("covQual_sf", rt.AddressOf(s1,"covQual_sf"),'covQual_sf/I')
    myTree.Branch("covQual_ff", rt.AddressOf(s1,"covQual_ff"),'covQual_ff/I')
    
    for paramName in paramNames:
        myTree.Branch('%s_ff'%paramName , rt.AddressOf(s1,'%s_ff'%paramName),'%s_ff/F' %paramName)
        myTree.Branch('%s_ff_error'%paramName , rt.AddressOf(s1,'%s_ff_error'%paramName),'%s_ff_error/F' %paramName)
        myTree.Branch('%s_sf'%paramName , rt.AddressOf(s1,'%s_sf'%paramName),'%s_sf/F' %paramName)
        myTree.Branch('%s_sf_error'%paramName , rt.AddressOf(s1,'%s_sf_error'%paramName),'%s_sf_error/F' %paramName)
    
    for ix in range(0, nBins):
        myTree.Branch("b%i_sf" %(ix) , rt.AddressOf(s1,"b%i_sf" %(ix)),'b%i_sf/F' %ix)
        myTree.Branch("b%i_ff" %(ix) , rt.AddressOf(s1,"b%i_ff" %(ix)),'b%i_ff/F' %ix)
        myTree.Branch("b%i_toy" %(ix) , rt.AddressOf(s1,"b%i_toy" %(ix)),'b%i_toy/F' %ix)

    os.system("rm tempMacro_%d.C"%rando)
    return s1


def runToys(w,options,cfg,seed):
    
    setStyle()
    rt.RooRandom.randomGenerator().SetSeed(seed)
    
    extRazorPdf = w.pdf('extRazorPdf')
    dataHist = w.data("data_obs")    
    if w.obj("fitresult_extRazorPdf_data_obs") != None:
        fr = w.obj("fitresult_extRazorPdf_data_obs")
    elif w.obj("nll_extRazorPdf_data_obs") != None:
        fr = w.obj("nll_extRazorPdf_data_obs")
    elif w.obj("fitresult_extRazorPdf_data_obs_with_constr") != None:
        fr = w.obj("fitresult_extRazorPdf_data_obs_with_constr")
    elif w.obj("nll_extRazorPdf_data_obs_with_constr") != None:
        fr = w.obj("nll_extRazorPdf_data_obs_with_constr")

    fr.Print("V")
    if options.r>-1:
        extSpBPdf = w.pdf('extSpBPdf')
                
    #nll = w.function('nll_extRazorPdf_data_obs')
    th1x = w.var("th1x")
    
    params = extRazorPdf.getParameters(dataHist)
    paramsToRemove = []
    for p in rootTools.RootIterator.RootIterator(params):
        if p.isConstant(): paramsToRemove.append(p)

    [params.remove(p) for p in paramsToRemove]
    paramNames = [p.GetName() for p in rootTools.RootIterator.RootIterator(params)]
    paramNames.sort()    
    if options.r>-1: paramNames.append('r')
        
    x = array('d', cfg.getBinning(options.box)[0]) # MR binning
    y = array('d', cfg.getBinning(options.box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(options.box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    
    th1x.setBins(nBins)
        
    if seed>-1:            
        output = rt.TFile.Open(options.outDir+'/genfittoys_Freq_s%i_%s.root'%(seed,options.box),'recreate')
    else:
        output = rt.TFile.Open(options.outDir+'/genfittoys_Freq_%s.root'%(options.box),'recreate')
        
    output.cd()
    myTree = rt.TTree("myTree", "myTree")
    
    s1 = getTree(myTree, paramNames, nBins, options.box, z)
    
    pSet = fr.floatParsFinal()
    for p in rootTools.RootIterator.RootIterator(pSet):
        w.var(p.GetName()).setVal(p.getVal())
        w.var(p.GetName()).setError(p.getError())
        asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('true'),rt.RooFit.Asimov())
    
        value = setattr(s1, 'toy_num', -1) #save the toy number (for differentiating from true pdf)
        iBinX = -1    
        for i in range(1,len(x)):
            for j in range(1,len(y)):
                for k in range(1,len(z)):
                    iBinX += 1
                    th1x.setVal(iBinX+0.5)
                    expected = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))                    
                    toy = float(asimov.weight(rt.RooArgSet(th1x)))
                    value = setattr(s1, 'b%i_ff'%iBinX, expected) #save predicted fit yield
                    value = setattr(s1, 'b%i_sf'%iBinX, expected) #save predicted fit yield
                    value = setattr(s1, 'b%i_toy'%iBinX, toy) #save toy yield
                  
    myTree.Fill()
    iToy = 0
    
    
    widgets = ['Running Freq toys ', Percentage(), ' ', Bar(marker=RotatingMarker()),' ', ETA(), ' ', FileTransferSpeed()]
    pbar = ProgressBar(widgets=widgets, max_value=options.nToys).start()
    while iToy < options.nToys:    
        pSet = fr.floatParsFinal()
        for p in rootTools.RootIterator.RootIterator(pSet):
            w.var(p.GetName()).setVal(p.getVal())
            w.var(p.GetName()).setError(p.getError())
            if 'Ntot' in p.GetName():                
                w.var(p.GetName()).setVal(options.scaleFactor*p.getVal())                
            #print "%s = %f +- %f"%(p.GetName(),p.getVal(),p.getError())


        #print "good pars"                        
        errorCountBefore = rt.RooMsgService.instance().errorCount()
        
        asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy'))

        errorCountAfter = rt.RooMsgService.instance().errorCount()   
        if errorCountAfter > errorCountBefore:
            #print "can't generate toy=%i"%iToy
            continue

        #print "SUCCESS: generated toy=%i"%iToy

        pSetSave = pSet
        migrad_status_sf = -1
        hesse_status_sf = -1
        minos_status_sf = -1
        migrad_status_ff = -1
        hesse_status_ff = -1
        minos_status_ff = -1

        
        sideband = convertSideband('LowMR,LowRsq',w,x,y,z)
    
        nll_func_toy_ff = extRazorPdf.createNLL(asimov,rt.RooFit.Extended(True))     
        nll_func_toy_sf = extRazorPdf.createNLL(asimov,rt.RooFit.Extended(True),rt.RooFit.Range(sideband))                 
        m = rt.RooMinimizer(nll_func_toy_ff)
        m.setStrategy(0)
        m.setPrintLevel(-1)
        m.setPrintEvalErrors(-1)
        migrad_status_ff = m.minimize('Minuit2','migrad')
        #hesse_status_ff = m.minimize('Minuit2','hesse')
        fr_ff = m.save()
        covQual_ff = fr_ff.covQual()

        
        value = setattr(s1, 'toy_num', iToy) #save the toy number (for differentiating from true pdf)
        value = setattr(s1, 'migrad_ff', migrad_status_ff) #save migrad status full fit
        value = setattr(s1, 'hesse_ff', hesse_status_ff) #save hesse status full fit
        value = setattr(s1, 'covQual_ff', covQual_ff) #save cov qual full fit
        iBinX = -1    
        for i in range(1,len(x)):
            for j in range(1,len(y)):
                for k in range(1,len(z)):
                    iBinX += 1
                    th1x.setVal(iBinX+0.5)
                    expected = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))                    
                    toy = float(asimov.weight(rt.RooArgSet(th1x)))
                    value = setattr(s1, 'b%i_ff'%iBinX, expected) #save predicted full fit yield
                    value = setattr(s1, 'b%i_toy'%iBinX, toy) #save toy yield

        # save full fit parameters
        for p in rootTools.RootIterator.RootIterator(fr_ff.floatParsFinal()):
            value = setattr(s1, p.GetName()+"_ff", p.getVal())
            value = setattr(s1, p.GetName()+"_ff_error", p.getError())    
                  
        m = rt.RooMinimizer(nll_func_toy_sf)
        m.setStrategy(0)
        m.setPrintLevel(-1)
        m.setPrintEvalErrors(-1)
        migrad_status_sf = m.minimize('Minuit2','migrad')
        #hesse_status_sf = m.minimize('Minuit2','hesse')
        fr_sf = m.save()
        covQual_sf = fr_sf.covQual()
        
        value = setattr(s1, 'migrad_sf', migrad_status_sf) #save migrad status sideband fit
        value = setattr(s1, 'hesse_sf', hesse_status_sf) #save hesse status sideband fit
        value = setattr(s1, 'covQual_sf', covQual_sf) #save cov qual sideband fit
        iBinX = -1    
        for i in range(1,len(x)):
            for j in range(1,len(y)):
                for k in range(1,len(z)):
                    iBinX += 1
                    th1x.setVal(iBinX+0.5)
                    expected = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
                    value = setattr(s1, 'b%i_sf'%iBinX, expected) #save predicted sideband fit yield
                    

        # save sideband fit parameters
        for p in rootTools.RootIterator.RootIterator(fr_sf.floatParsFinal()):
            value = setattr(s1, p.GetName()+"_sf", p.getVal())
            value = setattr(s1, p.GetName()+"_sf_error", p.getError())    


                    
        pbar.update(iToy)
        myTree.Fill()
        iToy+=1        
    rt.RooMsgService.instance().reset()
    pbar.finish()
    
    w.Print('v')
    output.cd()
    myTree.Write()
    w.Write()
    output.Close()
    return output.GetName()
    
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
    parser.add_option('-s','--seed',dest="seed", default=-1,type="int",
                  help="random seed")
    parser.add_option('-r','--signal-strength',dest="r", default=-1,type="float",
                  help="signal strength")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default=None,type="string",
                  help="input fit file")
    parser.add_option('-t','--toys',dest="nToys", default=3000,type="int",
                  help="number of toys")
    parser.add_option('--scale-factor',dest="scaleFactor", default=1,type="float",
                  help="scale factor for toys")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)    
    lumi = options.lumi
    inputFitFile = options.inputFitFile
    seed = options.seed
    
    lumi_in = 0.

    if inputFitFile is not None:
        rootFile = rt.TFile.Open(inputFitFile,"r")
        w = rootFile.Get("w"+options.box)
        if w.obj("fitresult_extRazorPdf_data_obs") != None:
            fr = w.obj("fitresult_extRazorPdf_data_obs")
        elif w.obj("nll_extRazorPdf_data_obs") != None:
            fr = w.obj("nll_extRazorPdf_data_obs")
        elif w.obj("fitresult_extRazorPdf_data_obs_with_constr") != None:
            fr = w.obj("fitresult_extRazorPdf_data_obs_with_constr")
        elif w.obj("nll_extRazorPdf_data_obs_with_constr") != None:
            fr = w.obj("nll_extRazorPdf_data_obs_with_constr")

        
    outputName = runToys(w,options,cfg,options.seed)

    print "writing tree to %s"%(outputName)
