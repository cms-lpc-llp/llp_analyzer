from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
import os
import sys

def fixPars(w, label, doFix=True, setVal=None):
    parSet = w.allVars()
    for par in rootTools.RootIterator.RootIterator(parSet):
        if label in par.GetName():
            par.setConstant(doFix)
            if setVal is not None: par.setVal(setVal)

def initializeWorkspace(w,cfg,box,scaleFactor=1.,penalty=False,x=None,y=None,z=None,emptyHist3D=None):
    
    if x is None or y is None or z is None:
        x = array('d', cfg.getBinning(box)[0]) # MR binning
        y = array('d', cfg.getBinning(box)[1]) # Rsq binning
        z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    #maxBins = 224
    maxBins = nBins
    
    parameters = cfg.getVariables(box, "combine_parameters")
    paramNames = []
    for parameter in parameters:
        if penalty and '_norm' in parameter:
            continue
        w.factory(parameter)
        paramName = parameter.split('[')[0]
        if not ("Cut" in paramName):
            paramNames.append(paramName)
            w.var(paramName).setConstant(False)
            
        
        # float normalization parameters
        fixPars(w,"Ntot",False)

        
        # fix Rsq MR cut parameters
        fixPars(w,"Cut")
        
        # fix Gaussian constraint parameters
        fixPars(w,"In")
        fixPars(w,"Mean")
        fixPars(w,"Sigma")
        
        # turn off shape parameters if no events in b-tag bin (not done yet)
        for k in range(0,len(z)-1):
            if "Ntot" in paramName and "%ib"%z[k] in paramName:
                w.var(paramName).setVal(scaleFactor * (w.data("RMRTree").sumEntries("nBtag>=%i && nBtag<%i"% (z[k],z[k+1] )) ))
                if not w.var(paramName).getVal():
                    fixPars(w,"%ib"%z[k])    
        
    w.factory('th1x[0,0,%i]'%maxBins)
    w.var('th1x').setBins(maxBins)
    if emptyHist3D==None:
        emptyHist3D = rt.TH3D("emptyHist3D","emptyHist3D",len(x)-1,x,len(y)-1,y,len(z)-1,z)
        iBinX = -1
        for ix in range(1,len(x)):
            for iy in range(1,len(y)):
                for iz in range(1,len(z)):
                    iBinX+=1
                    emptyHist3D.SetBinContent(ix,iy,iz,1)
                    emptyHist3D.SetBinError(ix,iy,iz,0)
    else:
        iBinX = -1
        for ix in range(1,len(x)):
            for iy in range(1,len(y)):
                for iz in range(1,len(z)):
                    iBinX+=1
                    if emptyHist3D.GetBinContent(ix,iy,iz)<=0:
                        emptyHist3D.SetBinContent(ix,iy,iz,1)
                    emptyHist3D.SetBinError(ix,iy,iz,0)
                    print ix, iy, iz, emptyHist3D.GetBinContent(ix,iy,iz)
        

    w.Print('v')
    commands = cfg.getVariables(box, "combine_pdfs")
    bkgs = []
    for command in commands:
        lower = command.lower()
        if lower.find('sum::')!=-1 or lower.find('prod::')!=-1 or lower.find('expr::')!=-1 or lower.find('roogaussian::')!=-1:
            w.factory(command)
        else:
            myclass = command.split('::')[0]
            remaining = command.split('::')[1]
            name = remaining.split('(')[0]
            altname = '_'.join(reversed(name.split('_')))
            mytuple = remaining.replace(name,'').replace('(','').replace(')','')
            mylist = mytuple.split(',')
            arglist = [name, name]
            for myvar in mylist:
                if w.var(myvar)!=None:
                    arglist.append(w.var(myvar))
                elif w.function(myvar)!=None:
                    arglist.append(w.function(myvar))
                elif 'MADD_pars_' in myvar:
                        parlist = rt.RooArgList(myvar)
                        listdef = ''
                        for iBinX in range(0,maxBins):
                            if w.function('MADD_bin%i_%s'%(iBinX,box))!=None:
                                parlist.add(w.function('MADD_bin%i_%s'%(iBinX,box)))
                            else:
                                for iz in range(0,len(z)-1):
                                    if iBinX%(len(z)-1)==iz:
                                        w.factory("expr::MADD_bin%i_%s('@0',MADD_TTj%ib_%s)"%(iBinX,box,z[iz],box))
                                parlist.add(w.function('MADD_bin%i_%s'%(iBinX,box)))
                        rootTools.Utils.importToWS(w,parlist)
                        arglist.append(parlist)
                        
            args = tuple(arglist)
            pdf = getattr(rt,myclass)(*args)
            if hasattr(pdf,'setTH3Binning'):
                pdf.setTH3Binning(emptyHist3D)
            rootTools.Utils.importToWS(w,pdf)
            bkg = name.split("_")
            if box in bkg: bkg.remove(box)
            bkgs.append("_".join(bkg))
    return paramNames, bkgs


def initializeWorkspace_noFit(w,cfg,box):
    parameters = cfg.getVariables(box, "combine_parameters")
    paramNames = []
    for parameter in parameters:
        w.factory(parameter)
        paramName = parameter.split('[')[0]
        if paramName.find("Cut")==-1:
            paramNames.append(paramName)
            w.var(paramName).setConstant(False)
        else:
            w.var(paramName).setConstant(True)
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    #maxBins = 224
    maxBins = nBins
    
    w.factory('th1x[0,0,%i]'%maxBins)

    return paramNames

def convertDataset2TH1(data, cfg, box, workspace, useWeight=False, th1Name = 'h', x = array('d',[]), y = array('d',[]), z = array('d',[])):
    """Get the cocktail dataset from the file"""
    
    row = data.get()

    MR = row['MR']
    Rsq = row['Rsq']
    nBtag = row['nBtag']
    
    varSet = rt.RooArgSet(MR,Rsq,nBtag)
    varList = rt.RooArgList(MR,Rsq,nBtag)
    varList2D = rt.RooArgList(MR,Rsq)
        
    mRmin = row['MR'].getMin()
    mRmax = row['MR'].getMax()
    rsqMin = row['Rsq'].getMin()
    rsqMax = row['Rsq'].getMax()
    nbtagMin = row['nBtag'].getMin()
    nbtagMax = row['nBtag'].getMax()

    if cfg is not None:    
        x = array('d', cfg.getBinning(box)[0]) # MR binning
        y = array('d', cfg.getBinning(box)[1]) # Rsq binning
        z = array('d', cfg.getBinning(box)[2]) # nBtag binning

    
    myTH3 = rt.TH3D(th1Name+box, th1Name+box, len(x)-1, x, len(y)-1, y, len(z)-1, z)
    myTH2 = rt.TH2D(th1Name+box+"2d", th1Name+box+"2d", len(x)-1, x, len(y)-1, y)
    if useWeight:
        myTH3.Sumw2()
        myTH2.Sumw2()

    # fills automatically with weight
    data.fillHistogram(myTH3, varList,"MR>%f && MR<%f && Rsq>%f && Rsq<%f && nBtag >= %f && nBtag <= %f"%(x[0],x[-1],y[0],y[-1],z[0],z[-1]))
    data.fillHistogram(myTH2, varList2D,"MR>%f && MR<%f && Rsq>%f && Rsq<%f && nBtag >= %f && nBtag <= %f"%(x[0],x[-1],y[0],y[-1],z[0],z[-1]))
    
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    #maxBins = 224
    maxBins = nBins
    
    if maxBins >= nBins:
        myTH1 = rt.TH1D(th1Name+box+"1d",th1Name+box+"1d",maxBins,0,maxBins)
    else:
        myTH1 = rt.TH1D(th1Name+box+"1d",th1Name+box+"1d",nBins,0,nBins)
        
    if useWeight:
        myTH1.Sumw2()
    i = 0
    for ix in range(1,len(x)):
        for iy in range(1,len(y)):
            for iz in range(1,len(z)):
                i+= 1
                myTH1.SetBinContent(i,myTH3.GetBinContent(ix,iy,iz))
                if useWeight: myTH1.SetBinError(i,myTH3.GetBinError(ix,iy,iz))

    return myTH1


def writeDataCard(box,model,txtfileName,bkgs,paramNames,w,penalty,fixed,shapes=[]):
        obsRate = w.data("data_obs").sumEntries()
        nBkgd = len(bkgs)
        rootFileName = txtfileName.replace('.txt','.root')
        signals = len(model.split('p'))
        if signals>1:
                rates = [w.data("%s_%s"%(box,sig)).sumEntries() for sig in model.split('p')]
                processes = ["%s_%s"%(box,sig) for sig in model.split('p')]
                lumiErrs = [1.027 for sig in model.split('p')]
        else:
                rates = [w.data("%s_%s"%(box,model)).sumEntries()]
                processes = ["%s_%s"%(box,model)]
                lumiErrs = [1.027]
        rates.extend([w.var('Ntot_%s_%s'%(bkg,box)).getVal() for bkg in bkgs])
        processes.extend(["%s_%s"%(box,bkg) for bkg in bkgs])
        lumiErrs.extend([1.00 for bkg in bkgs])
        divider = "------------------------------------------------------------\n"
        datacard = "imax 1 number of channels\n" + \
                   "jmax %i number of processes minus 1\n"%(nBkgd+signals-1) + \
                   "kmax * number of nuisance parameters\n" + \
                   divider + \
                   "observation	%.3f\n"%obsRate + \
                   divider + \
                   "shapes * * %s w%s:$PROCESS w%s:$PROCESS_$SYSTEMATIC\n"%(rootFileName,box,box) + \
                   divider
        binString = "bin"
        processString = "process"
        processNumberString = "process"
        rateString = "rate"
        lumiString = "lumi\tlnN"
        for i in range(0,len(bkgs)+signals):
            binString +="\t%s"%box
            processString += "\t%s"%processes[i]
            processNumberString += "\t%i"%(i-signals+1)
            rateString += "\t%.3f" %rates[i]
            lumiString += "\t%.3f"%lumiErrs[i]
        binString+="\n"; processString+="\n"; processNumberString+="\n"; rateString +="\n"; lumiString+="\n"
        datacard+=binString+processString+processNumberString+rateString+divider
        # now nuisances
        datacard+=lumiString
        for shape in shapes:
            shapeString = '%s\tshape\t'%shape
            for sig in range(0,signals):
                shapeString += '\t1.0'
            for i in range(0,len(bkgs)):
                shapeString += '\t-'
            shapeString += '\n'
            datacard+=shapeString
        for paramName in paramNames:
            if fixed:
                fixPars(w,paramName)    
            elif 'Mean' in paramName or 'Sigma' in paramName:
                fixPars(w,paramName)
            elif 'MR1_' in paramName:
                mean = w.var(paramName.replace('MR1','MR1Mean')).getVal()
                sigma = w.var(paramName.replace('MR1','MR1Sigma')).getVal()
                datacard += "%s\tparam\t%e\t%e\n"%(paramName,mean,sigma)
            elif 'MADD_' in paramName:
                mean = 0.
                sigma = 1.
                datacard += "%s\tparam\t%e\t%e\n"%(paramName,mean,sigma)                
            elif penalty:                    
                mean = w.var(paramName).getVal()
                sigma = w.var(paramName).getError()                
                if "Ntot" in paramName:                    
                    effectString = ''
                    for sig in range(0,signals):
                        effectString += "\t1.0"           
                    for bkg in bkgs:
                        if bkg in paramName:
                            effectString += "\t%.3f"%(1.0+sigma/mean)                            
                        else:
                            effectString += "\t1.0"                    
                    datacard += "%s\tlnN%s\n"%(paramName.replace("Ntot","Norm"),effectString)
                else:
                    datacard += "%s\tparam\t%e\t%e\n"%(paramName,mean,sigma)
                         
            else:
                if "Ntot" in paramName:
                    effectString = ""
                    for sig in range(0,signals):
                        effectString += "\t1.0"
                    for bkg in bkgs:
                        if bkg in paramName:
                            effectString += "\t2.0"                         
                        else:
                            effectString += "\t1.0"                    
                    datacard += "%s\tlnU%s\n"%(paramName.replace("Ntot","Norm"),effectString)
                else:
                    datacard += "%s\tflatParam\n"%(paramName)
            
        txtfile = open(txtfileName,"w")
        txtfile.write(datacard)
        txtfile.close()

        
def writeDataCard_noFit(box,model,txtfileName,bkgs,paramNames,w):
        obsRate = w.data("data_obs").sumEntries()
        nBkgd = len(bkgs)
        rootFileName = txtfileName.replace('.txt','.root')
        rates = [w.data("%s_%s"%(box,model)).sumEntries()]
        rates.extend([w.data("%s_%s"%(box,bkg)).sumEntries() for bkg in bkgs])
        processes = ["%s_%s"%(box,model)]
        processes.extend(["%s_%s"%(box,bkg) for bkg in bkgs])
        lumiErrs = [1.05]
        lumiErrs.extend([1.00 for bkg in bkgs])
        mcErrs = {} #dictionary of uncorrelated mc bkgd lnN uncertainties
        for bkg in bkgs:
                mcErrs[bkg] = [1.00]
                mcErrs[bkg].extend([1.00 + 0.10*(bkg==bkg1) for bkg1 in bkgs])
                
        divider = "------------------------------------------------------------\n"
        datacard = "imax 1 number of channels\n" + \
                   "jmax %i number of backgrounds\n"%nBkgd + \
                   "kmax * number of nuisance parameters\n" + \
                   divider + \
                   "observation	%.3f\n"%obsRate + \
                   divider + \
                   "shapes * * %s w%s:$PROCESS w%s:$PROCESS_$SYSTEMATIC\n"%(rootFileName,box,box) + \
                   divider
                   
        binString = "bin"
        processString = "process"
        processNumberString = "process"
        rateString = "rate"
        lumiString = "lumi\tlnN"
        for i in range(0,len(bkgs)+1):
            binString +="\t%s"%box
            processString += "\t%s"%processes[i]
            processNumberString += "\t%i"%i
            rateString += "\t%.3f" %rates[i]
            lumiString += "\t%.3f"%lumiErrs[i]
        binString+="\n"; processString+="\n"; processNumberString+="\n"; rateString +="\n"; lumiString+="\n"
            
        mcErrStrings = {}
        for bkg in bkgs:
                mcErrStrings[bkg] = "%s_%s_norm\tlnN"%(box,bkg)
                for i in range(0,len(bkgs)+1):                
                        mcErrStrings[bkg] += "\t%.3f"%mcErrs[bkg][i]
                mcErrStrings[bkg]+="\n"
                
        datacard+=binString+processString+processNumberString+rateString+divider
        
        # now nuisances
        datacard+=lumiString
        
        for bkg in bkgs:
                datacard+=mcErrStrings[bkg]

        txtfile = open(txtfileName,"w")
        txtfile.write(datacard)
        txtfile.close()
        
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
    parser.add_option('--no-fit',dest="noFit",default=False,action='store_true',
                  help="Turn off fit (use MC directly)")
    parser.add_option('--print-yields',dest="printYields",default=False,action='store_true',
                  help="print yields")
    parser.add_option('--penalty',dest="penalty",default=False,action='store_true',
                  help="penalty terms on background shape + norm parameters from input fit result")
    parser.add_option('--fixed',dest="fixed",default=False,action='store_true',
                  help="fixed background shape + norm parameters")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default=None,type="string",
                  help="input fit file")
    parser.add_option('--no-signal-sys',dest="noSignalSys",default=False,action='store_true',
                  help="no signal systematic shape uncertainties")
    parser.add_option('--no-backgrounds',dest="noBackgrounds",default=False,action='store_true',
                  help="no background components, just create datacard for signal")
    parser.add_option('--histo-file',dest="histoFile", default=None,type="string",
                  help="input histogram file for MADD/fit systematic")
    #pdf uncertainty options.  current prescription is just to take 10% uncorrelated error on each bin
    #parser.add_option('--num-pdf-weights',dest="numPdfWeights",default=0,type='int',
                  #help='number of pdf nuisance parameters to use')
    #parser.add_option('--compute-pdf-envelope',dest="computePdfEnvelope",default=False,action='store_true',
                  #help="Use the SUS pdf reweighting prescription, summing weights in quadrature")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box = options.box
    lumi = options.lumi
    noFit = options.noFit
    printYields = options.printYields

    lumi_in = 0.
    signalFileName = ''
    model = ''
    massPoint = ''
    for f in args:
        if f.lower().endswith('.root'):
            if f.lower().find('t1')!=-1 or f.lower().find('t5')!=-1 or f.lower().find('t2')!=-1:
                signalFileName = f
                #signalDs = workspace.data('RMRTree')
                model = f.split('.root')[0].split('-')[1].split('_')[0]
                massPoint = '_'.join(f.split('.root')[0].split('_')[1:3])
            else:
                rootFile = rt.TFile(f)
                workspace = rootFile.Get('w'+box)
                data = workspace.data('RMRTree')
            lumi_in = 1000.*float([g.replace('lumi-','') for g in f.split('_') if g.find('lumi')!=-1][0])

    w = rt.RooWorkspace("w"+box)
    
    rootTools.Utils.importToWS(w,data)
    
    if noFit:
        paramNames = initializeWorkspace_noFit(w,cfg,box)
    else:
        if options.histoFile is not None:
            histoFile = rt.TFile.Open(options.histoFile,"r")
            emptyHist3D = histoFile.Get('FullFit_MADD_Fix_%s'%box)       
            paramNames, bkgs = initializeWorkspace(w,cfg,box,scaleFactor=lumi/lumi_in,penalty=False,x=None,y=None,z=None,emptyHist3D=emptyHist3D)
        else:
            paramNames, bkgs = initializeWorkspace(w,cfg,box,scaleFactor=lumi/lumi_in,penalty=options.penalty)
    
    
    th1x = w.var('th1x')
    
    myTH1 = convertDataset2TH1(data, cfg, box, w)
    myTH1.Scale(lumi/lumi_in)
    dataHist = rt.RooDataHist("data_obs", "data_obs", rt.RooArgList(th1x), rt.RooFit.Import(myTH1))
    rootTools.Utils.importToWS(w,dataHist)
    
    if noFit:        
        z = array('d', cfg.getBinning(box)[2]) # nBtag binning
        for k in range(0,len(z)-1):
            data_red = data.reduce("nBtag>=%i && nBtag<%i"%(z[k],z[k+1]))
            myTH1_red = convertDataset2TH1(data_red, cfg, box, w)
            myTH1_red.Scale(lumi/lumi_in)
            dataHist_red = rt.RooDataHist("%s_%s"%(box, "TTj%ib"%z[k]), "%s_%s"%(box, "TTj%ib"%z[k]), rt.RooArgList(th1x), rt.RooFit.Import(myTH1_red))
            rootTools.Utils.importToWS(w,dataHist_red)

    elif options.inputFitFile is not None:
        inputRootFile = rt.TFile.Open(options.inputFitFile,"r")
        wIn = inputRootFile.Get("w"+box).Clone("wIn"+box)            
        if wIn.obj("fitresult_extRazorPdf_data_obs") != None:
            frIn = wIn.obj("fitresult_extRazorPdf_data_obs")
        elif wIn.obj("nll_extRazorPdf_data_obs") != None:
            frIn = wIn.obj("nll_extRazorPdf_data_obs")
        elif wIn.obj("fitresult_extRazorPdf_data_obs_with_constr") != None:
            fr = wIn.obj("fitresult_extRazorPdf_data_obs_with_constr")
        elif wIn.obj("nll_extRazorPdf_data_obs_with_constr") != None:
            frIn = wIn.obj("nll_extRazorPdf_data_obs_with_constr")
        print "restoring parameters from fit"
        frIn.Print("V")
        for p in rootTools.RootIterator.RootIterator(frIn.floatParsFinal()):
            w.var(p.GetName()).setVal(p.getVal())
            w.var(p.GetName()).setError(p.getError())
            if "Ntot" in p.GetName():
                normNameList = p.GetName().replace("Ntot_","").split("_")
                normNameList.reverse()
                normNameList.append("norm")
                normName = "_".join(normNameList)
                #w.var(normName).setError(w.var(p.GetName()).getError()/w.var(p.GetName()).getVal())
            
    
    signalHistos = []
    signalFile = rt.TFile.Open(signalFileName)
    names = [k.GetName() for k in signalFile.GetListOfKeys()]
    for name in names:
        d = signalFile.Get(name)
        if isinstance(d, rt.TH1):
            #d.SetDirectory(rt.gROOT)
            signalHistos.append(d)
            sigDataHist = rt.RooDataHist(d.GetName(), d.GetName(), rt.RooArgList(th1x), rt.RooFit.Import(d))
            rootTools.Utils.importToWS(w,sigDataHist)
    
    if options.noSignalSys:
        shapes = []
    else:
        shapes = ['jes','muontrig','eletrig','btag','btagfastsim','facscale','renscale','facrenscale','ees','mes','pileup','isr','mcstat%s'%box.lower(),'pdf%s'%box.lower()]
        if 'mcstat%s'%box.lower() in shapes:
            for hist in signalHistos:
                if 'mcstat%s'%box.lower() in hist.GetName() and 'Up' in hist.GetName():
                    num = hist.GetName().split('mcstat%s'%box.lower())[-1].split('Up')[0]
                    shapes.append('mcstat%s%s'%(box.lower(),num))                   
            shapes.remove('mcstat%s'%box.lower())            
        if 'pdf%s'%box.lower() in shapes:
            for hist in signalHistos:
                if 'pdf%s'%box.lower() in hist.GetName() and 'Up' in hist.GetName():
                    num = hist.GetName().split('pdf%s'%box.lower())[-1].split('Up')[0]
                    shapes.append('pdf%s%s'%(box.lower(),num))                   
            shapes.remove('pdf%s'%box.lower())            
        #if options.computePdfEnvelope:
            #shapes.append('pdfenvelope')
        #else:
            #shapes.extend(['n'+str(n)+'pdf' for n in range(options.numPdfWeights)]) 
        if box in ["MultiJet", "DiJet", "FourJet", "SixJet", "LooseLeptonMultiJet", "LooseLeptonDiJet", "LooseLeptonFourJet", "LooseLeptonSixJet"]:
            shapes.extend(['vetomuoneff','vetoeleeff','vetomuonfastsim','vetoelefastsim'])
        else:
            shapes.extend(['tightmuoneff','tighteleeff','tightmuonfastsim','tightelefastsim'])
        
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    btagMin = z[0]
    btagMax = z[-1]        
    if btagMax>btagMin+1:            
        outFile = 'razor_combine_%s_%s_lumi-%.3f_%i-%ibtag_%s.root'%(model,massPoint,lumi/1000.,btagMin,btagMax-1,box)
    else:
        outFile = 'razor_combine_%s_%s_lumi-%.3f_%ibtag_%s.root'%(model,massPoint,lumi/1000.,btagMin,box)
    
    outputFile = rt.TFile.Open(options.outDir+"/"+outFile,"recreate")
    if noFit:
        writeDataCard_noFit(box,model,options.outDir+"/"+outFile.replace(".root",".txt"),["TTj%ib"%iz for iz in z[:-1]],paramNames,w)
    else:
        writeDataCard(box,model,options.outDir+"/"+outFile.replace(".root",".txt"),bkgs,paramNames,w,options.penalty,options.fixed,shapes=shapes)
    w.Write()
    os.system("cat %s"%options.outDir+"/"+outFile.replace(".root",".txt"))

    if printYields:
        x = array('d', cfg.getBinning(box)[0]) # MR binning
        y = array('d', cfg.getBinning(box)[1]) # Rsq binning
        z = array('d', cfg.getBinning(box)[2]) # nBtag binning
        zBins = len(z)-1
        yBins = len(y)-1
        xBins = len(x)-1
        csvFile = open(options.outDir+"/"+outFile.replace(".root",".csv"),'w')
        csvOutput = "bin number,MR range,Rsq range,b-tags,signal yield (S),background yield from MC (B),background yield from fit (F),S/B,S/F"
        for iBin in range(0,th1x.getBins()):
            signal = sigDataHist.sumEntries("th1x>=%i && th1x<%i+1"%(iBin,iBin))
            mc = w.data("data_obs").sumEntries("th1x>=%i && th1x<%i+1"%(iBin,iBin))
            asimov = w.pdf("extRazorPdf").generateBinned(rt.RooArgSet(th1x),rt.RooFit.Asimov())
            fit = asimov.sumEntries("th1x>=%i && th1x<%i+1"%(iBin,iBin))
            th1x.setVal(iBin+0.5)
            fit = w.pdf("extRazorPdfNoConst").getVal()
            zBin = iBin % zBins
            yBin = ( (iBin - zBin)/(zBins) ) % (yBins)
            xBin =  (iBin - zBin - yBin*zBins ) / (zBins*yBins)
            csvOutput += "\n%i,%i-%i,%.2f-%.2f,%i,%f,%f,%f,=E%i/F%i,=E%i/G%i" % (iBin, x[xBin],x[xBin+1],y[yBin],y[yBin+1],z[zBin],signal,mc,fit,iBin+2,iBin+2,iBin+2,iBin+2)
        csvFile.write(csvOutput)
        csvFile.close()
        print ""
        print "yields written into %s"%(options.outDir+"/"+outFile.replace(".root",".csv"))
