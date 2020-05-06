from optparse import OptionParser
import ROOT as rt
import sys
import rootTools
import glob
from math import *
from framework import Config
import os
from array import *
from GChiPairs import gchipairs

def getFileName(hybridLimit, mg, mchi, box, model, lumi, btag, directory, method,t):
    if hybridLimit == "higgsCombineToys":
        modelPoint = "%i_%i"%(mg,mchi)
        fileName = "%s/%s%s_%s_lumi-%.3f_%s_%s_%i.%s.mH120.root"%(directory,hybridLimit,model,modelPoint,lumi,btag,box,t,method)
    else:
        modelPoint = "%i_%i"%(mg,mchi)
        fileName = "%s/%s%s_%s_lumi-%.3f_%s_%s.%s.mH120.root"%(directory,hybridLimit,model,modelPoint,lumi,btag,box,method)
    return fileName


def writeXsecTree(box, model, directory, mg, mchi, xsecULObs, xsecULExpPlus2, xsecULExpPlus, xsecULExp, xsecULExpMinus, xsecULExpMinus2):
    outputFileName = "%s/xsecUL_mg_%s_mchi_%s_%s.root" %(directory, mg, mchi, box)
    print "INFO: xsec UL values being written to %s"%outputFileName
    fileOut = rt.TFile.Open(outputFileName, "recreate")
    
    xsecTree = rt.TTree("xsecTree", "xsecTree")
    try:
        from ROOT import MyStruct
    except ImportError:
        myStructCmd = "struct MyStruct{Double_t mg;Double_t mchi; Double_t x; Double_t y;"
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
    xsecTree.Branch("mg", rt.AddressOf(s,"mg"),'mg/D')
    xsecTree.Branch("mchi", rt.AddressOf(s,"mchi"),'mchi/D')
    xsecTree.Branch("x", rt.AddressOf(s,"x"),'x/D')
    xsecTree.Branch("y", rt.AddressOf(s,"y"),'y/D')
    
    
    s.mg = mg
    s.mchi = mchi
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
    xsecTree.Branch("xsecULObs_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+0)),'xsecUL%i/D'%(ixsecUL+0))
    xsecTree.Branch("xsecULExpPlus2_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+1)),'xsecUL%i/D'%(ixsecUL+1))
    xsecTree.Branch("xsecULExpPlus_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+2)),'xsecUL%i/D'%(ixsecUL+2))
    xsecTree.Branch("xsecULExp_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+3)),'xsecUL%i/D'%(ixsecUL+3))
    xsecTree.Branch("xsecULExpMinus_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+4)),'xsecUL%i/D'%(ixsecUL+4))
    xsecTree.Branch("xsecULExpMinus2_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+5)),'xsecUL%i/D'%(ixsecUL+5))
    exec 's.xsecUL%i = xsecULObs[ixsecUL]'%(ixsecUL+0)
    exec 's.xsecUL%i = xsecULExpPlus2[ixsecUL]'%(ixsecUL+1)
    exec 's.xsecUL%i = xsecULExpPlus[ixsecUL]'%(ixsecUL+2)
    exec 's.xsecUL%i = xsecULExp[ixsecUL]'%(ixsecUL+3)
    exec 's.xsecUL%i = xsecULExpMinus[ixsecUL]'%(ixsecUL+4)
    exec 's.xsecUL%i = xsecULExpMinus2[ixsecUL]'%(ixsecUL+5)
    ixsecUL += 4

    xsecTree.Fill()

    fileOut.cd()
    xsecTree.Write()
    
    fileOut.Close()
    
    return outputFileName

if __name__ == '__main__':

    
    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('-l','--lumi',dest="lumi", default=0.210,type="float",
                  help="lumi in fb^-1, e.g.: 0.210")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Input/Output directory to store output")
    parser.add_option('--signif',dest="doSignificance",default=False,action='store_true',
                  help="for significance instead of limit")
    parser.add_option('--toys',dest="doHybridNew",default=False,action='store_true',
                  help="for toys instead of asymptotic")
    parser.add_option('--xsec-file',dest="refXsecFile",default="./data/gluino13TeV.txt",type="string",
                  help="Input directory")
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")

    (options,args) = parser.parse_args()

    boxInput = options.box
    model = options.model
    lumi = options.lumi
    directory = options.outDir
    doHybridNew = options.doHybridNew
    doSignificance = options.doSignificance
    refXsecFile = options.refXsecFile
    

    cfg = Config.Config(options.config)

    boxes = boxInput.split('_')

    btag = ''
    for box in boxes:            
        z = array('d', cfg.getBinning(box)[2]) # nBtag binning
        btagMin = z[0]
        btagMax = z[-1]
        if btagMax-1>btagMin:          
            btag = '%i-%ibtag'%(btagMin,btagMax-1)
        else:
            btag = '%ibtag'%(btagMin)            

    #output = rt.TFile.Open("%s/combine_%s_%s.root"%(directory,model,boxInput),"RECREATE")
    #rt.TTree("combine","combine")

    haddOutputs = []

    #mgMin, mgMax, mchiMin, mchiMax, binWidth, nRebins, xsecMin, xsecMax, diagonalOffset, smoothing = getModelSettings(model)
    #sigHist = rt.TH2D("significance","significance",int((mgMax-mgMin)/binWidth),mgMin, mgMax,int((mchiMax-mchiMin)/binWidth), mchiMin, mchiMax)

    thyXsec = {}
    thyXsecErr = {}
    if refXsecFile is not None:
        print "INFO: Input ref xsec file!"
        for mg, mchi in gchipairs(model):
            for line in open(refXsecFile,'r'):
                line = line.replace('\n','')
                if str(mg)==line.split(',')[0]:
                    thyXsec[mg,mchi] = float(line.split(',')[1]) #pb
                    thyXsecErr[mg,mchi] = 0.01*float(line.split(',')[2])                    
            if model=="T5ttttDM175T2tt":                    
                for line in open('data/stop13TeV.txt','r'):
                    line = line.replace('\n','')
                    if str(mchi+175)==line.split(',')[0]:
                        thyXsecStop = float(line.split(',')[1]) #pb
                #thyXsec[mg,mchi]+=thyXsecStop                    
                
    for mg, mchi in gchipairs(model):
        if refXsecFile is not None:
            refXsec = 1.e3*thyXsec[mg,mchi]
            #print "INFO: ref xsec taken to be: %s mass %d, xsec = %f fb"%(gluinoHistName, mg, refXsec)
        
        if doSignificance and doHybridNew:
            if not glob.glob(getFileName("higgsCombineSignif",mg,mchi,boxInput,model,lumi,btag,directory,"HybridNew",0)): continue
            print "INFO: opening %s"%(getFileName("higgsCombineSignif",mg,mchi,boxInput,model,lumi,btag,directory,"HybridNew",0))
            tFile = rt.TFile.Open(getFileName("higgsCombineSignif",mg,mchi,boxInput,model,lumi,btag,directory,"HybridNew",0))
        elif doHybridNew: 
            if not glob.glob(getFileName("higgsCombineToys",mg,mchi,boxInput,model,lumi,btag,directory,"HybridNew",0)): continue
            print "INFO: opening %s"%(getFileName("higgsCombineToys",mg,mchi,boxInput,model,lumi,btag,directory,"HybridNew",0))
            tFile = rt.TFile.Open(getFileName("higgsCombineToys",mg,mchi,boxInput,model,lumi,btag,directory,"HybridNew",0))
        elif doSignificance: 
            if not glob.glob(getFileName("higgsCombine",mg,mchi,boxInput,model,lumi,btag,directory,"ProfileLikelihood",0)): continue
            print "INFO: opening %s"%(getFileName("higgsCombine",mg,mchi,boxInput,model,lumi,btag,directory,"ProfileLikelihood",0))
            tFile = rt.TFile.Open(getFileName("higgsCombine",mg,mchi,boxInput,model,lumi,btag,directory,"ProfileLikelihood",0))
        else:
            if not glob.glob(getFileName("higgsCombine",mg,mchi,boxInput,model,lumi,btag,directory,"Asymptotic",0)): continue
            print "INFO: opening %s"%(getFileName("higgsCombine",mg,mchi,boxInput,model,lumi,btag,directory,"Asymptotic",0))
            tFile = rt.TFile.Open(getFileName("higgsCombine",mg,mchi,boxInput,model,lumi,btag,directory,"Asymptotic",0))

        try:
            if tFile.InheritsFrom("TFile") is False:
                continue
        except:
            continue
            
        limit = tFile.Get("limit")
        try:
            if limit.InheritsFrom("TTree") is False: 
                tFile.cd()
                tFile.Close()
                continue
        except:
            tFile.cd()
            tFile.Close()
            continue
        if doSignificance and limit.GetEntries() < 1: 
            tFile.cd()
            tFile.Close()
            continue
        if (not doSignificance) and limit.GetEntries() < 6: 
            tFile.cd()
            tFile.Close()
            continue
        limit.Draw('>>elist','','entrylist')
        elist = rt.gDirectory.Get('elist')
        entry = elist.Next()
        limit.GetEntry(entry)
        limits = []
        while True:
            if entry == -1: break
            limit.GetEntry(entry)
            if doSignificance:
                limits.append(max(0.0,limit.limit))
            else:
                limits.append(refXsec*(1.e-3)*limit.limit)
            entry = elist.Next()
        tFile.cd()
        tFile.Close()
            
        limits.reverse()
        print mg, mchi
        print limits
        
        #if doSignificance:
        #    sigHist.SetBinContent(sigHist.FindBin(mg,mchi),limits[0])
        #else:

        haddOutput = writeXsecTree(boxInput, model, directory, mg, mchi, [limits[0]],[limits[1]],[limits[2]],[limits[3]],[limits[4]],[limits[5]])
        haddOutputs.append(haddOutput)


    #if doSignificance:
    #    c = plotSignificance(boxInput,model,sigHist,doHybridNew)
    #else:
    if doHybridNew:
        os.system("hadd -f %s/xsecUL_HybridNew_%s.root %s"%(directory,boxInput," ".join(haddOutputs)))
    else:
        os.system("hadd -f %s/xsecUL_Asymptotic_%s.root %s"%(directory,boxInput," ".join(haddOutputs)))
