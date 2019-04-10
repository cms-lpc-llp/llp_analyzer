from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
from WriteDataCard import *
import os
import random
import sys
import math
from PlotFit import *
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, FileTransferSpeed, FormatLabel, Percentage, ProgressBar, ReverseBar, RotatingMarker, SimpleProgress, Timer
        
def getTree(myTree,paramNames,nBins,box,z):
    
    rando = random.randint(1,999999)
    # first structure
    stringMyStruct1 = "void tempMacro_%d(){struct MyStruct1{"%(rando)
    stringMyStruct1 += "float toy_num; float nll_%s; float n2llr_%s; float chi2_%s;"%(box,box,box)
    for k in range(1,len(z)):
        ibtag = z[k-1]      
        stringMyStruct1 += "float nll_%ibtag_%s; float n2llr_%ibtag_%s; float chi2_%ibtag_%s;"%(ibtag,box,ibtag,box,ibtag,box)
    stringMyStruct1 += "int covQual_%s; int migrad_%s; int hesse_%s; int minos_%s;"%(box,box,box,box)
    for paramName in paramNames:
        stringMyStruct1 = stringMyStruct1+"float %s; float %s_error;" %(paramName,paramName)
        if paramName=='r':
            stringMyStruct1 = stringMyStruct1+"float r_errorlo; float r_errorhi;"           
    for iBinX in range(0,nBins):
        stringMyStruct1 = stringMyStruct1+"float b%i;" %(iBinX)
    tempMacro = open("tempMacro_%d.C"%rando,"w")
    tempMacro.write(stringMyStruct1+"};}")
    tempMacro.close()
    rt.gROOT.Macro("tempMacro_%d.C"%rando)
    #rt.gROOT.ProcessLine(".x tempMacro_%d.C"%rando)
    #rt.gROOT.ProcessLine(stringMyStruct1+"};")
    from ROOT import MyStruct1

    # fills the bins list and associate the bins to the corresponding variables in the structure
    s1 = MyStruct1()
    myTree.Branch('toy_num' , rt.AddressOf(s1,'toy_num'),'toy_num/F')
    myTree.Branch('nll_%s'%box , rt.AddressOf(s1,'nll_%s'%box),'nll_%s/F' %box)
    myTree.Branch('n2llr_%s'%box , rt.AddressOf(s1,'n2llr_%s'%box),'n2llr_%s/F' %box)
    myTree.Branch('chi2_%s'%box , rt.AddressOf(s1,'chi2_%s'%box),'chi2_%s/F' %box)    
    for k in range(1,len(z)):
        ibtag = z[k-1]      
        myTree.Branch('nll_%ibtag_%s'%(ibtag,box) , rt.AddressOf(s1,'nll_%ibtag_%s'%(ibtag,box)),'nll_%ibtag_%s/F' %(ibtag,box))
        myTree.Branch('n2llr_%ibtag_%s'%(ibtag,box) , rt.AddressOf(s1,'n2llr_%ibtag_%s'%(ibtag,box)),'n2llr_%ibtag_%s/F' %(ibtag,box))
        myTree.Branch('chi2_%ibtag_%s'%(ibtag,box) , rt.AddressOf(s1,'chi2_%ibtag_%s'%(ibtag,box)),'chi2_%ibtag_%s/F' %(ibtag,box))
    myTree.Branch('covQual_%s'%box , rt.AddressOf(s1,'covQual_%s'%box),'covQual_%s/I' %box)
    myTree.Branch('migrad_%s'%box , rt.AddressOf(s1,'migrad_%s'%box),'migrad_%s/I' %box)
    myTree.Branch('hesse_%s'%box , rt.AddressOf(s1,'hesse_%s'%box),'hesse_%s/I' %box)
    myTree.Branch('minos_%s'%box , rt.AddressOf(s1,'minos_%s'%box),'minos_%s/I' %box)
    for paramName in paramNames:
        myTree.Branch(paramName , rt.AddressOf(s1,paramName),'%s/F' %paramName)
        myTree.Branch('%s_error'%paramName , rt.AddressOf(s1,'%s_error'%paramName),'%s_error/F' %paramName)
        if paramName=='r':            
            myTree.Branch('r_errorlo' , rt.AddressOf(s1,'r_errorlo'),'r_errorlo/F')
            myTree.Branch('r_errorhi' , rt.AddressOf(s1,'r_errorhi'),'r_errorhi/F')
    for ix in range(0, nBins):
        myTree.Branch("b%i" %(ix) , rt.AddressOf(s1,"b%i" %(ix)),'b%i/F' %ix)

    os.system("rm tempMacro_%d.C"%rando)
    return s1


def runToys(w,options,cfg,seed):
    
    setStyle()
    if seed>-1:
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
    
    fitband = convertSideband(options.fitRegion,w,x,y,z)
    sideband = convertSideband('LowMR,LowRsq',w,x,y,z)
    ixMin = 3
    iyMin = 3
    if options.box in ['MuMultiJet','EleMultiJet']:
        if x[2]==500:
            ixMin = 3
        else:
            ixMin = 2
        iyMin = 3

    unc = 'Bayes'
    if options.varyN and options.noStat: unc = "Bayes_varyN_noStat"
    elif options.varyN: unc = "Bayes_varyN"
    elif options.noStat: unc = "Bayes_noStat"
    elif options.noSys: unc = "Bayes_noSys"
    elif options.oneSigma: unc = 'oneSigma'

    if options.freq and options.noStat: unc = 'Freq_noStat'
    elif options.freq: unc = 'Freq_varyN'
        
    if options.r>-1:
        rString = str('%.3f'%options.r).replace(".","p")
        if seed>-1:
            output = rt.TFile.Open(options.outDir+'/toys_%s_r%s_s%i_%s.root'%(unc,rString,seedoptions.box),'recreate')
        else:
            output = rt.TFile.Open(options.outDir+'/toys_%s_r%s_%s.root'%(unc,rString,options.box),'recreate')
    else:
        if seed>-1:            
            output = rt.TFile.Open(options.outDir+'/toys_%s_s%i_%s.root'%(unc,seed,options.box),'recreate')
        else:
            output = rt.TFile.Open(options.outDir+'/toys_%s_%s.root'%(unc,options.box),'recreate')
        
    output.cd()
    myTree = rt.TTree("myTree", "myTree")
    
    s1 = getTree(myTree, paramNames, nBins, options.box, z)
    value =  setattr(s1, 'toy_num', -1) # set toy number to -1
    
    for p in rootTools.RootIterator.RootIterator(fr.floatParsFinal()):
        w.var(p.GetName()).setVal(p.getVal())
        w.var(p.GetName()).setError(p.getError())
        value = setattr(s1, p.GetName(), p.getVal())
        value = setattr(s1, p.GetName()+'_error', p.getError())
    if options.r>-1:
        value =  setattr(s1, 'r', 0)
        value =  setattr(s1, 'r_error', 0)
        value =  setattr(s1, 'r_errorlo', 0)
        value =  setattr(s1, 'r_errorhi', 0)
        
    asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('central'),rt.RooFit.Asimov())

    chi2_data = 0
    n2llr_data = 0
    nll_data = 0
    chi2_data_btag = [0 for k in range(1,len(z))]
    n2llr_data_btag = [0 for k in range(1,len(z))]
    nll_data_btag = [0 for k in range(1,len(z))]
    bestFitByBin = []
    iBinX = -1    
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            for k in range(1,len(z)):
                iBinX += 1
                th1x.setVal(iBinX+0.5)                
                #expected = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))  
                expected = float(asimov.weight(rt.RooArgSet(th1x)))
                bestFitByBin.append(expected)
                observed = float(dataHist.weight(rt.RooArgSet(th1x)))
                value = setattr(s1, 'b%i'%iBinX, expected)
                if expected>0:
                    chi2_data += ( observed - expected ) * ( observed - expected ) / ( expected )
                    chi2_data_btag[k-1] += ( observed - expected ) * ( observed - expected ) / ( expected )
                    #if k==1:                        
                    #    print "bin", i, j, k, "observed:          ", observed
                    #    print "bin", i, j, k, "expected:          ", expected
                    #    print "bin", i, j, k, "chi2 contribution: ", ( observed - expected ) * ( observed - expected ) / ( expected )
                    #    print "bin", k, "chi2 sum so far:       ", chi2_data_btag[k-1]
                if expected>0:            
                    nll_data -= observed*rt.TMath.Log(expected) - expected
                    nll_data_btag[k-1] -= observed*rt.TMath.Log(expected) - expected
                if observed>0 and expected>0:
                    n2llr_data += 2 * ( observed*rt.TMath.Log(observed/expected) - observed )
                    n2llr_data_btag[k-1] += 2 * ( observed*rt.TMath.Log(observed/expected) - observed )
                if expected>0:                
                    n2llr_data += 2 * ( expected )
                    n2llr_data_btag[k-1] += 2 * ( expected )
                    #if k==1:                        
                    #    print "bin", i, j, k, "n2llr contribution:", 2 * ( observed*rt.TMath.Log(observed/expected) - observed ) + 2 * ( expected )
                    #    print "bin", k, "n2llr sum so far:      ", n2llr_data_btag[k-1]
        
    value = setattr(s1, 'nll_%s'%options.box, nll_data)
    value = setattr(s1, 'n2llr_%s'%options.box, n2llr_data)
    value = setattr(s1, 'chi2_%s'%options.box, chi2_data)
    
    for k in range(1,len(z)):
        ibtag = z[k-1]
        value = setattr(s1, 'nll_%ibtag_%s'%(ibtag,options.box), nll_data_btag[k-1])
        value = setattr(s1, 'n2llr_%ibtag_%s'%(ibtag,options.box), n2llr_data_btag[k-1])
        value = setattr(s1, 'chi2_%ibtag_%s'%(ibtag,options.box), chi2_data_btag[k-1])
    
        
    myTree.Fill()
        
    iToy = 0
    
    nBadPars = 0
    pBest = fr.floatParsFinal()
    pBestVal = {}
    pBestErr = {}
    #mu = rt.RooArgList()
    #hesseParams = rt.RooArgList()
    #iArray = 0
    #nIndexArray = []
    #xFactor = {}
    for p in rootTools.RootIterator.RootIterator(pBest):
        pBestVal[p.GetName()] = p.getVal()
        pBestErr[p.GetName()] = p.getError()
        #p.setConstant(True)
        #mu.add(p)
        #hesseParams.add(w.var(p.GetName()))
        #if 'n_TTj' in p.GetName() or 'b_TTj' in p.GetName():
        #if 'MultiJet' in p.GetName():
        #    nIndexArray.append(iArray)
        #    if 'TTj0b' in p.GetName():
        #        xFactor[iArray] = 1.7
        #    elif 'TTj1b' in p.GetName():
        #        xFactor[iArray] = 1.2
        #    elif 'TTj2b' in p.GetName():
        #        xFactor[iArray] = 1.15
        #    elif 'TTj3b' in p.GetName():
        #        xFactor[iArray] = 1.15
        #iArray+=1
    #maxArray = iArray
        
    #covMatrix = fr.covarianceMatrix()
    #corrMatrix = fr.correlationMatrix()
    #covMatrixClone = covMatrix.Clone(covMatrix.GetName()+"_varyN")
    
    # double the uncertainty for each n parameter
    #print nIndexArray
    #print xFactor
    #for nIndex in nIndexArray:
    #    for otherIndex in nIndexArray:        
    #        covMatrixClone[nIndex][otherIndex] = xFactor[nIndex]*xFactor[otherIndex]*covMatrix[nIndex][otherIndex]
    #        covMatrixClone[otherIndex][nIndex] = xFactor[otherIndex]*xFactor[nIndex]*covMatrix[otherIndex][nIndex]

    #hesseParams.Print('v')
    #mu.Print('v')
    #covMatrix.Print('v')
    #covMatrixClone.Print("V")
    #hessePdf = rt.RooMultiVarGaussian('hessePdf','hessePdf',hesseParams,mu,covMatrixClone)
    #hessePdf.Print('v')
    #hesseDs = hessePdf.generate(params,int(100*options.nToys))
    #corrMatrix.Print('v')
    #c = rt.TCanvas('c','c',500,400)
    #varName = 'n_TTj1b_MultiJet'
    #varName2 = 'b_TTj1b_MultiJet'
    #frame = w.var(varName).frame(rt.RooFit.Range(0,10))
    #frame.Print('v')
    #hesseDs.plotOn(frame)
    #hessePdf.plotOn(frame)
    #frame.Draw()
    #c.Print(varName+'.2pdf')
    #w.var(varName).setMin(0)
    #w.var(varName2).setMin(0)
    #w.var(varName).setMax(10)
    #w.var(varName2).setMax(2.5)
    #hist2d = hesseDs.createHistogram(w.var(varName),w.var(varName2),100,100)
    #hist2d.GetXaxis().SetTitle(varName)
    #hist2d.GetYaxis().SetTitle(varName2)    
    #hist2d.Draw("colz")
    #c.Print(varName+varName2+'.2pdf')
    #sys.exit()
                
    if options.box in ['MultiJet','DiJet','MultiJet_0b','MultiJet_1b','MultiJet_2b','DiJet_0b','DiJet_1b','DiJet_2b']:
        xFactor = [1.8, 1.4, 1.4, 1.4] #xFactor for each b-tag bin
    elif options.box in ['MuMultiJet', 'LeptonMultiJet', 'LeptonJet','LeptonJet_0b','LeptonJet_1b','LeptonJet_2b','LeptonMultiJet_0b','LeptonMultiJet_1b','LeptonMultiJet_2b']:
        xFactor = [2.0, 2.0, 2.0, 2.0] #xFactor for each b-tag bin
    elif options.box == 'EleMultiJet':
        xFactor = [2.0, 1.8, 1.2, 1.2] #xFactor for each b-tag bin
    widgets = ['Running %s toys '%unc, Percentage(), ' ', Bar(marker=RotatingMarker()),' ', ETA(), ' ', FileTransferSpeed()]
    pbar = ProgressBar(widgets=widgets, maxval=options.nToys).start()
    iAttempt = -1
    while iToy < options.nToys:
        iAttempt+=1
        if options.freq:
            pSet = fr.floatParsFinal()
        elif options.oneSigma:
            pSet = fr.randomizePars()
        else:
            if options.noSys:                
                pSet = fr.floatParsFinal()
            elif options.varyN:
                #pSet = hesseDs.get(iAttempt)
                pSet = fr.randomizePars()
            else:                
                pSet = fr.randomizePars()
        for p in rootTools.RootIterator.RootIterator(pSet):
            w.var(p.GetName()).setVal(p.getVal())
            w.var(p.GetName()).setError(p.getError())

        badPars = []
        for bkgd in ['TTj0b','TTj1b','TTj2b','TTj3b']:
            if w.var('n_%s_%s'%(bkgd,options.box))!=None:                
                badPars.append(w.var('n_%s_%s'%(bkgd,options.box)).getVal() <= 0)
            if w.var('b_%s_%s'%(bkgd,options.box))!=None:                                
                badPars.append(w.var('b_%s_%s'%(bkgd,options.box)).getVal() <= 0)
            if w.var('MR0_%s_%s'%(bkgd,options.box))!=None:                                
                badPars.append(w.var('MR0_%s_%s'%(bkgd,options.box)).getVal() <= 0)
                #badPars.append(w.var('MR0_%s_%s'%(bkgd,options.box)).getVal() >= w.var('MR').getMin())
            if w.var('R0_%s_%s'%(bkgd,options.box))!=None:                   
                badPars.append(w.var('R0_%s_%s'%(bkgd,options.box)).getVal() <= 0)
                #badPars.append(w.var('R0_%s_%s'%(bkgd,options.box)).getVal() >= w.var('Rsq').getMin())
        if any(badPars):
            nBadPars+=1
            #print "bad pars toy=%i"%iToy
            #print badPars
            continue

        #print "good pars"                        
        errorCountBefore = rt.RooMsgService.instance().errorCount()

        badVal = False
        for iBinX in range(0,nBins):
            th1x.setVal(iBinX+0.5) # check number of events in each bin
            pdfValV = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
            pdfVal0 = extRazorPdf.getValV(0) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
            if bestFitByBin[iBinX] > 0 and pdfValV/bestFitByBin[iBinX] <= 1e-12:
            #if bestFitByBin[iBinX] > 0 and pdfValV <= 0:
                #print "bin = %i"%iBinX
                #print "best fit = %e"%(bestFitByBin[iBinX])
                #print "pdf valv = %e"%(pdfValV)
                #print "pdf val0 = %e"%(pdfVal0)
                badVal = True                
        if badVal:
            #print "bad val"
            continue
        
        errorCountAfter = rt.RooMsgService.instance().errorCount()
        if errorCountAfter > errorCountBefore:            
            print "can't evaulate pdf toy=%i"%iToy
            continue
        
        
        errorCountBefore = rt.RooMsgService.instance().errorCount()        
        #print "start generating toy=%i"%iToy
        if options.freq:
            if options.r>-1:
                w.var('r').setVal(options.r)
                asimov = extSpBPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy'),rt.RooFit.Extended(True))
            else:
                asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy'),rt.RooFit.Extended(True))
        elif options.noStat:         
            if options.r>-1:
                w.var('r').setVal(options.r)            
                asimov = extSpBPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy'),rt.RooFit.Asimov())
            else:
                asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy'),rt.RooFit.Asimov())
        else:
            if options.r>-1:                
                w.var('r').setVal(options.r)                
                asimov = extSpBPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy'),rt.RooFit.Extended(True))
            else:
                asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy'),rt.RooFit.Extended(True))

        #print "toy entries = %.2f"%asimov.sumEntries()
        errorCountAfter = rt.RooMsgService.instance().errorCount()
        if errorCountAfter > errorCountBefore:
            print "can't generate toy=%i"%iToy
            continue

        #print "SUCCESS: generated toy=%i"%iToy
        
        pSetSave = pSet
        migrad_status = -1
        hesse_status = -1
        minos_status = -1
        if options.freq:                      
            if options.r>-1:
                nll_func_toy = extSpBPdf.createNLL(asimov,rt.RooFit.Extended(True))
                m = rt.RooMinimizer(nll_func_toy)
                m.setStrategy(0)
                m.setPrintLevel(-1)
                m.setPrintEvalErrors(-1)
                rSet = rt.RooArgSet(w.var('r'))
                migrad_status = m.minimize('Minuit2','migrad')
                #hesse_status = m.minimize('Minuit2','hesse')
                #minos_status = m.minos(rSet)                
                fr_toy = m.save()
                value = setattr(s1,'migrad_%s'%options.box, migrad_status)   
                value = setattr(s1,'hesse_%s'%options.box, hesse_status)
                value = setattr(s1,'minos_%s'%options.box, minos_status)
            else:                
                #print "yes"
                nll_func_toy = extRazorPdf.createNLL(asimov,rt.RooFit.Extended(True),rt.RooFit.Range(fitband))               
                m = rt.RooMinimizer(nll_func_toy)
                m.setStrategy(0)
                m.setPrintLevel(-1)
                m.setPrintEvalErrors(-1)
                migrad_status = m.minimize('Minuit2','migrad')
                #migrad_status = m.minimize('Minuit2','migrad')
                #migrad_status = m.minimize('Minuit2','migrad')
                fr_toy = m.save()
            value = setattr(s1,'covQual_%s'%options.box, fr_toy.covQual())   
            value = setattr(s1,'migrad_%s'%options.box, migrad_status)   
            value = setattr(s1,'hesse_%s'%options.box, hesse_status)
            value = setattr(s1,'minos_%s'%options.box, minos_status)
            pSetSave = fr_toy.floatParsFinal()
            
                
        for p in rootTools.RootIterator.RootIterator(pSetSave):
            value = setattr(s1, p.GetName(), p.getVal())
            value = setattr(s1, p.GetName()+"_error", p.getError())
            if p.GetName()=='r':
                value = setattr(s1, "r_errorlo", p.getAsymErrorLo())
                value = setattr(s1, "r_errorhi", p.getAsymErrorHi())
                

        chi2_toy = 0
        n2llr_toy = 0
        nll_toy = 0
        chi2_toy_btag = [0 for k in range(1,len(z))]
        n2llr_toy_btag = [0 for k in range(1,len(z))]
        nll_toy_btag = [0 for k in range(1,len(z))]
        # restore best-fit to calculate expected values        
        for p in rootTools.RootIterator.RootIterator(pSetSave):
            w.var(p.GetName()).setVal(p.getVal())
            w.var(p.GetName()).setError(p.getError())

        #if options.varyN:
        #    xGaus = rt.RooRandom.randomGenerator().Gaus(0,1)                
            
        iBinX = -1
        for i in range(1,len(x)):
            for j in range(1,len(y)):
                for k in range(1,len(z)):
                    iBinX += 1
                    th1x.setVal(iBinX+0.5)                    
                    inSideband = False
                    #print "x, y = ", x[i-1], y[j-1]
                    #print "xMin, yMin = ", x[ixMin-1], y[iyMin-1]
                    if x[i-1] < x[ixMin-1]:
                        inSideband = True
                    if y[j-1] < y[iyMin-1]:
                        inSideband = True
                    #print "inSideband = %s"%inSideband
                    expected = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
                    if options.noStat and options.varyN:                        
                        if inSideband:
                            toy = float(asimov.weight(rt.RooArgSet(th1x)))
                        else:
                            xGaus = rt.RooRandom.randomGenerator().Gaus(0,1)                    
                            toy = float(asimov.weight(rt.RooArgSet(th1x)))*rt.TMath.Power(xFactor[k-1],xGaus)
                    elif options.varyN:                        
                        if inSideband:
                            toy = float(asimov.weight(rt.RooArgSet(th1x)))
                        else:
                            xGaus = rt.RooRandom.randomGenerator().Gaus(0,1)
                            central = float(expected*rt.TMath.Power(xFactor[k-1],xGaus))
                            toy = rt.RooRandom.randomGenerator().Poisson(central)
                    elif options.noStat:
                        toy = float(asimov.weight(rt.RooArgSet(th1x)))
                    else:
                        toy = float(asimov.weight(rt.RooArgSet(th1x)))
                    observed = float(dataHist.weight(rt.RooArgSet(th1x)))  
                    if options.oneSigma:
                        toy = observed # to get nll with respect to original dataset
                        value = setattr(s1, 'b%i'%iBinX, expected) #save expected yield
                    elif options.freq and options.noStat:
                        #print "expected = ", expected, toy, observed
                        value = setattr(s1, 'b%i'%iBinX, expected) 
                    elif options.freq:
                        withStat = rt.RooRandom.randomGenerator().Poisson(expected)
                        #print "with stats = ", withStat
                        value = setattr(s1, 'b%i'%iBinX, withStat) 
                    else:          
                        value = setattr(s1, 'b%i'%iBinX, toy)

                        #print "observed, expected, toy = ", observed, expected, toy

                    if expected>0:
                        chi2_toy += ( toy - expected ) * ( toy - expected ) / ( expected )
                        chi2_toy_btag[k-1] += ( toy - expected ) * ( toy - expected ) / ( expected )
                    if toy>0 and expected>0:
                        n2llr_toy += 2 * ( toy*rt.TMath.Log(toy/expected) - toy )
                        n2llr_toy_btag[k-1] += 2 * ( toy*rt.TMath.Log(toy/expected) - toy )
                    if expected>0:
                        n2llr_toy += 2 * ( expected )
                        n2llr_toy_btag[k-1] += 2 * ( expected )
                    if expected>0:
                        nll_toy -= toy*rt.TMath.Log(expected) - expected
                        nll_toy_btag[k-1] -= toy*rt.TMath.Log(expected) - expected

        # to check  nll, chi2 calculation
        #nll_func_toy = extRazorPdf.createNLL(asimov,rt.RooFit.Extended(True))
        #chi2_func_toy = extRazorPdf.createChi2(asimov,rt.RooFit.Extended(True),rt.RooFit.DataError(rt.RooAbsData.Expected))
        
        #print ''
        #print "chi2 func:   ", chi2_func_toy.getVal()
        #print "chi2 by hand  ", chi2_toy
        #print "nll func:    ", nll_func_toy.getVal()
        #print "nll by hand: ", nll_toy
        if options.oneSigma:
            nsigma = 1.0
            nparam = fr.floatParsFinal().getSize()
            if 2*(nll_toy-nll_data) > rt.Math.chisquared_quantile(rt.Math.erf(nsigma/rt.TMath.Sqrt(2.)),nparam): continue
            
        value = setattr(s1, 'nll_%s'%options.box, nll_toy)
        value = setattr(s1, 'n2llr_%s'%options.box, n2llr_toy)
        value = setattr(s1, 'chi2_%s'%options.box, chi2_toy)
    
        for k in range(1,len(z)):
            ibtag = z[k-1]
            value = setattr(s1, 'nll_%ibtag_%s'%(ibtag,options.box), nll_toy_btag[k-1])
            value = setattr(s1, 'n2llr_%ibtag_%s'%(ibtag,options.box), n2llr_toy_btag[k-1])
            value = setattr(s1, 'chi2_%ibtag_%s'%(ibtag,options.box), chi2_toy_btag[k-1])
    

        
        value =  setattr(s1, 'toy_num', iToy) # save toy number
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
    parser.add_option('-t','--toys',dest="nToys", default=3000,type="int",
                  help="number of toys")
    parser.add_option('-s','--seed',dest="seed", default=-1,type="int",
                  help="random seed")
    parser.add_option('--no-stat',dest="noStat",default=False,action='store_true',
                  help="no statistical uncertainty, just systematic uncertainty, default is statstical + systematic uncertainty")
    parser.add_option('--no-sys',dest="noSys",default=False,action='store_true',
                  help="no systematic uncertainty, just statistical uncertainty, default is statstical + systematic uncertainty")
    parser.add_option('--freq',dest="freq",default=False,action='store_true',
                  help="refit each toy with only statistical fluctuations, as in frequentist approach; default is bayeseian")
    parser.add_option('--one-sigma',dest="oneSigma",default=False,action='store_true',
                  help="sample parameters uniformly (within 2-sigma) and save yields if NLL is within 1-sigma of minimum")
    parser.add_option('-r','--signal-strength',dest="r", default=-1,type="float",
                  help="signal strength => do each fit the the SpB pdf")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default=None,type="string",
                  help="input fit file")
    parser.add_option('--fit-region',dest="fitRegion",default="Full",type="string",
                  help="Fit region")
    parser.add_option('--vary-n',dest="varyN",default=False,action='store_true',
                  help="vary n parameters in addition to standard hesse pdf")
    
    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)
    
    lumi = options.lumi

    inputFitFile = options.inputFitFile    
    
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
