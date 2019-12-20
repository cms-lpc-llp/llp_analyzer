import time
from optparse import OptionParser
import os
import ROOT as rt
from array import *
from framework import Config
import sys
import glob
from GChiPairs import gchipairs
    
def writeBashScript(box,btag,model,mg,mchi,lumi,config,submitDir,isData,fit,penalty,inputFitFile,noSignalSys,histoFile,min_tol,min_strat,numPdfWeights,computePdfEnvelope,rMax=-1):
    
    massPoint = "%i_%i"%(mg, mchi)
    dataString = ''
    if isData:
        dataString = '--data'
        
    signalSys = ''
    if noSignalSys:
        signalSys = '--no-signal-sys'

    fitString = ''
    if fit:
        fitString = '--fit'
        
    penaltyString = ''
    if penalty:
        penaltyString = '--penalty'
                
    histoString = ''
    if histoFile:
        histoString = '--histo-file %s'%(options.histoFile)
        
    computePdfEnvelopeString = ''
    if computePdfEnvelope:
        computePdfEnvelopeString = '--compute-pdf-envelope'
        
    # prepare the script to run
    outputname = submitDir+"/submit_"+model+"_"+massPoint+"_lumi-%.3f_"%(lumi)+btag+"_"+box+".src"
        
    ffDir = submitDir+"/logs_"+model+"_"+massPoint+"_"+btag+"_"+box
    user = os.environ['USER']
    pwd = os.environ['PWD']
    
    combineDir = "/afs/cern.ch/work/%s/%s/RAZORRUN2/Limits/%s/"%(user[0],user,model) # directory where combine output files will be copied
    cmsswBase = "/afs/cern.ch/work/%s/%s/RAZORRUN2/CMSSW_7_1_5"%(user[0],user) # directory where 'cmsenv' will be run (needs to have combine and RazorAnalyzer setup)

    script =  '#!/usr/bin/env bash -x\n'
    script += 'mkdir -p %s\n'%combineDir
    
    script += 'echo $SHELL\n'
    script += 'pwd\n'
    script += 'cd %s/src/RazorAnalyzer \n'%(cmsswBase)
    script += 'pwd\n'
    script += "export SCRAM_ARCH=slc6_amd64_gcc481\n"
    script += "export CMSSW_BASE=%s\n"%(cmsswBase)
    script += 'eval `scramv1 runtime -sh`\n'
    script += 'cd - \n'
    script += "export TWD=${PWD}/%s_%s_lumi-%.3f_%s_%s\n"%(model,massPoint,lumi,btag,box)
    script += "mkdir -p $TWD\n"
    script += "cd $TWD\n"
    script += 'pwd\n'
    script += 'git clone git@github.com:RazorCMS/RazorAnalyzer\n'
    script += 'cd RazorAnalyzer\n'
    script += 'git checkout -b Limits Limits20160420v2\n'
    script += 'source setup.sh\n'
    script += 'make\n'
    script += 'mkdir -p Datasets\n'
    script += 'mkdir -p %s\n'%submitDir
    if "T1" in model or "T5" in model:
        script += 'python python/RunCombine.py -i %s -m %s --mGluino %i --mLSP %i %s -c %s --lumi-array %f -d %s -b %s %s %s %s --min-tol %e --min-strat %i --rMax %f %s\n'%(inputFitFile,model,mg,mchi,dataString,config,lumi,submitDir,box,fitString,penaltyString,signalSys,min_tol,min_strat,rMax,histoString)
    else:
        script += 'python python/RunCombine.py -i %s -m %s --mStop %i --mLSP %i %s -c %s --lumi-array %f -d %s -b %s %s %s %s --min-tol %e --min-strat %i --rMax %f %s\n'%(inputFitFile,model,mg,mchi,dataString,config,lumi,submitDir,box,fitString,penaltyString,signalSys,min_tol,min_strat,rMax,histoString)
    script += 'cp %s/higgsCombine* %s/\n'%(submitDir,combineDir) 
    script += 'cd ../..\n'
    script += 'rm -rf $TWD\n'
    
    outputfile = open(outputname,'w')
    outputfile.write(script)
    outputfile.close

    return outputname,ffDir



if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('-l','--lumi',dest="lumi", default=0.210,type="float",
                  help="lumi in fb^-1, e.g.: 0.210")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('--fit',dest="fit",default=False,action='store_true',
                  help="Turn on pre-fit")
    parser.add_option('--data',dest="isData", default=False,action='store_true',
                  help="changes for data")
    parser.add_option('--penalty',dest="penalty",default=False,action='store_true',
                  help="penalty terms on background parameters")
    parser.add_option('--no-sub',dest="noSub", default=False,action='store_true',
                  help="no submission")
    parser.add_option('-q','--queue',dest="queue",default="1nh",type="string",
                  help="queue: 1nh, 8nh, 1nd, etc.")
    parser.add_option('--mg-geq',dest="mgMin",default=-1,type="float",
                  help="mgMin ")
    parser.add_option('--mg-lt',dest="mgMax",default=10000,type="float",
                  help="mgMax ")
    parser.add_option('--mchi-geq',dest="mchiMin",default=-1,type="float",
                  help="mchiMin ")
    parser.add_option('--mchi-lt',dest="mchiMax",default=10000,type="float",
                  help="mchiMax ")
    parser.add_option('--done-file',dest="doneFile",default=None,type="string",
                  help="file containing output files")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default='FitResults/BinnedFitResults.root',type="string",
                  help="input fit file")
    parser.add_option('--no-signal-sys',dest="noSignalSys",default=False,action='store_true',
                  help="no signal shape systematic uncertainties")
    parser.add_option('--min-tol',dest="min_tol",default=0.001,type="float",
                  help="minimizer tolerance (default = 0.001)")
    parser.add_option('--min-strat',dest="min_strat",default=2,type="int",
                  help="minimizer strategy (default = 2)")
    parser.add_option('--asymptotic-file',dest="asymptoticFile", default=None,type="string",
                  help="input file with asymptotic limit results (to set rMax dynamically based on expected limit)")
    parser.add_option('--num-pdf-weights',dest="numPdfWeights",default=0,type='int',
                  help='number of pdf nuisance parameters to use')
    parser.add_option('--compute-pdf-envelope',dest="computePdfEnvelope",default=False,action='store_true',
                  help="Use the SUS pdf reweighting prescription, summing weights in quadrature")
    parser.add_option('--histo-file',dest="histoFile", default=None,type="string",
                  help="input histogram file for MADD/fit systematic")
    (options,args) = parser.parse_args()


    cfg = Config.Config(options.config)

    boxes = options.box.split('_')

    btag = ''
    for box in boxes:            
        z = array('d', cfg.getBinning(box)[2]) # nBtag binning
        btagMin = z[0]
        btagMax = z[-1]
        if btagMax-1>btagMin:          
            btag = '%i-%ibtag'%(btagMin,btagMax-1)
        else:
            btag = '%ibtag'%(btagMin)
                
    nJobs = 0
    donePairs = []
    if options.doneFile is not None:
        with open(options.doneFile,'r') as f:            
            allFiles = [ line.replace('\n','') for line in f.readlines()]
            for (mg, mchi) in gchipairs(options.model):
                outputname = 'higgsCombine%s_%i_%i_lumi-%.3f_%s_%s.Asymptotic.mH120.root'%(options.model,mg,mchi,options.lumi,btag,options.box)
                if outputname in allFiles: donePairs.append((mg,mchi))

    thyXsec = {}
    if "T1" in options.model or "T5" in options.model:
        xsecFile = 'data/gluino13TeV.txt'
    elif "T2" in options.model:
        xsecFile = 'data/stop13TeV.txt'
        
    for line in open(xsecFile,'r'):
        for (mg, mchi) in gchipairs(options.model):
            if str(int(mg))==line.split(',')[0]:
                thyXsec[(mg,mchi)] = float(line.split(',')[1]) #pb
    if options.model=="T5ttttDM175T2tt":
        for line in open('data/stop13TeV.txt','r'):
            for (mg, mchi) in gchipairs(options.model):
                if str(int(mchi+175))==line.split(',')[0]:
                    thyXsec[(mg,mchi)] += float(line.split(',')[1]) #pb
            

    if options.asymptoticFile != None:        
        print "INFO: Input ref xsec file!"
        asymptoticRootFile = rt.TFile.Open(options.asymptoticFile,"READ")
        expMinus2 = asymptoticRootFile.Get("xsecUL_ExpMinus2_%s_%s"%(options.model,options.box))
        expPlus2 = asymptoticRootFile.Get("xsecUL_ExpPlus2_%s_%s"%(options.model,options.box))
        expMinus = asymptoticRootFile.Get("xsecUL_ExpMinus_%s_%s"%(options.model,options.box))
        expPlus = asymptoticRootFile.Get("xsecUL_ExpPlus_%s_%s"%(options.model,options.box))
        exp = asymptoticRootFile.Get("xsecUL_Exp_%s_%s"%(options.model,options.box))
        
    for (mg, mchi) in gchipairs(options.model):
        if not (mg >= options.mgMin and mg < options.mgMax): continue
        if not (mchi >= options.mchiMin and mchi < options.mchiMax): continue
        if (mg, mchi) in donePairs: continue
        nJobs+=1

        rMax = -1
        if options.asymptoticFile != None:
            rExpP2 = expPlus2.GetBinContent(expPlus2.FindBin(mg,mchi)) / thyXsec[(mg,mchi)]
            rExpP = expPlus.GetBinContent(expPlus.FindBin(mg,mchi)) / thyXsec[(mg,mchi)]
            rExp = exp.GetBinContent(exp.FindBin(mg,mchi)) / thyXsec[(mg,mchi)]
            rExpM = expMinus.GetBinContent(expMinus.FindBin(mg,mchi)) / thyXsec[(mg,mchi)]
            rExpM2 = expMinus2.GetBinContent(expMinus2.FindBin(mg,mchi)) / thyXsec[(mg,mchi)]
            rMaxThresholds = [50.,20., 10., 5., 2., 1., 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
            for rMaxTest in rMaxThresholds:
                if rExpM2<rMaxTest: rMax = 2*rMaxTest
            print "expected limit (+2sigma) = %f"%rExpP2
            print "expected limit (+1sigma) = %f"%rExpP
            print "expected limit           = %f"%rExp
            print "expected limit (-1sigma) = %f"%rExpM
            print "expected limit (-2sigma) = %f"%rExpM2
            print "=>          setting rMax = %f"%rMax

        outputname,ffDir = writeBashScript(options.box,btag,
                                           options.model,mg,mchi,
                                           options.lumi,options.config,
                                           options.outDir,options.isData,
                                           options.fit,options.penalty,
                                           options.inputFitFile,options.noSignalSys,
                                           options.histoFile,
                                           options.min_tol,options.min_strat,
                                           options.numPdfWeights,options.computePdfEnvelope,
                                           rMax)
        
        pwd = os.environ['PWD']
        os.system("mkdir -p "+pwd+"/"+ffDir)
        os.system("echo bsub -q "+options.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)        
        #os.system("source "+pwd+"/"+outputname)
        if not options.noSub:
            time.sleep(3)
            os.system("bsub -q "+options.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)

    print "nJobs = %i"%nJobs


