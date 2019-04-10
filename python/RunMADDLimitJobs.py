import os
import glob
import argparse
import ROOT as rt
from array import array

#local imports
from GChiPairs import gchipairs
    
def writeBashScript(tag, box, model, mg, mchi, submitDir, 
        noSys, bkgDir=None, saveWorkspace=False, noBoostCuts=False,
        fineGrained=False):
    
    massPoint = "%i_%i"%(mg, mchi)

    optString = ''
    if noSys:
        optString += '--no-sys --no-stat'
    if saveWorkspace:
        optString += ' --save-workspace'
    if noBoostCuts:
        optString += ' --no-boost-cuts'
    if fineGrained:
        optString += ' --fine-grained'

    particleString = '--mStop'
    if 'T1' in model or 'T5' in model:
        particleString = '--mGluino'

    # prepare the script to run
    outputname = ('Limits/'+submitDir+"/"+box+"/submit_"+model+"_"
            +massPoint+box+".src")
        
    ffDir = 'Limits/'+submitDir+"/"+box+"/logs_"+model+"_"+massPoint+"_"+box
    cmsswBase = os.environ['CMSSW_BASE']

    combineDir = "/eos/cms/store/group/phys_susy/razor/Run2Analysis/Limits/RazorInclusive2016/%s/%s/"%(
            submitDir,model) 

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
    script += "export TWD=${PWD}/%s_%s_%s\n"%(model,massPoint,box)
    script += "mkdir -p $TWD\n"
    script += "cd $TWD\n"
    script += 'pwd\n'
    script += 'git clone git://github.com/RazorCMS/RazorAnalyzer.git\n'
    script += 'cd RazorAnalyzer\n'
    script += 'git checkout -b LimitSetting LimitsMADD20180714\n' 
    script += 'make\n'
    script += 'mkdir -p %s\n'%submitDir
    script += 'python python/WriteRazorMADDCard.py'
    script += ' --tag %s --box %s'%(tag, box)
    if bkgDir is not None:
        script += ' --bkg-dir %s'%(bkgDir)
    script += ' --dir %s --model %s %s %i --mLSP %i %s\n'%(
            submitDir, model, particleString, mg, mchi, optString)
    script += 'cp %s/higgsCombine* %s/\n'%(submitDir,combineDir) 
    script += 'cp %s/RazorInclusiveMADD_* %s/\n'%(submitDir,combineDir) 
    script += 'cd ../..\n'
    script += 'rm -rf $TWD\n'
    
    outputfile = open(outputname,'w')
    outputfile.write(script)
    outputfile.close

    return outputname,ffDir


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tag',default='Razor2016_MoriondRereco')
    parser.add_argument('--box',help="box name",required=True)
    parser.add_argument('--model',help="signal model name",required=True)
    parser.add_argument('--dir',dest="outDir",default="./",
              help="Output directory to store cards")
    parser.add_argument('--bkg-dir', dest='bkgDir',
              help="Directory containing background histograms")
    parser.add_argument('--no-sub',dest="noSub",action='store_true',
              help="no submission")
    parser.add_argument('--queue',default="1nh",
              help="queue: 1nh, 8nh, 1nd, etc.")
    parser.add_argument('--mg-geq',dest="mgMin",default=-1,type=float)
    parser.add_argument('--mg-lt',dest="mgMax",default=10000,type=float)
    parser.add_argument('--mchi-geq',dest="mchiMin",
              default=-1,type=float)
    parser.add_argument('--mchi-lt',dest="mchiMax",default=10000,
              type=float)
    parser.add_argument('--done-file',dest="doneFile",
              help="file containing output file names")
    parser.add_argument('--no-sys',dest="noSys",action='store_true',
              help="no shape systematic uncertainties")
    parser.add_argument('--save-workspace',dest='saveWorkspace', 
              action='store_true',help='save workspace in output file')
    parser.add_argument('--no-boost-cuts', dest='noBoostCuts',
            action='store_true')
    parser.add_argument('--fine-grained', dest='fineGrained', action='store_true',
            help='Use fineGrained option for WriteRazorMADDCard')

    args = parser.parse_args()
    boxes = args.box.split('_')

    nJobs = 0
    donePairs = []
    if args.doneFile is not None:
        combineMethod = 'Asymptotic'
        with open(args.doneFile,'r') as f:            
            allFiles = [line.replace('\n','') for line in f.readlines()]
            for (mg, mchi) in gchipairs(args.model):
                outputname = 'higgsCombineMADD_%s_SMS-%s_%i_%i.%s.mH120.root'%(
                        args.box,args.model,mg,mchi,combineMethod)
                if outputname in allFiles: donePairs.append((mg,mchi))

    thyXsec = {}
    if "T1" in args.model or "T5" in args.model:
        xsecFile = 'data/gluino13TeV.txt'
    else:
        xsecFile = 'data/stop13TeV.txt'
        
    for line in open(xsecFile,'r'):
        for (mg, mchi) in gchipairs(args.model):
            if str(int(mg))==line.split(',')[0]:
                thyXsec[(mg,mchi)] = float(line.split(',')[1]) #pb

    for (mg, mchi) in gchipairs(args.model):
        if not (mg >= args.mgMin and mg < args.mgMax): continue
        if not (mchi >= args.mchiMin and mchi < args.mchiMax): continue
        if (mg, mchi) in donePairs: 
            print (mg,mchi),"is already done; skipping"
            continue
        nJobs+=1

        pwd = os.environ['PWD']
        os.system("mkdir -p "+pwd+"/Limits/"+args.outDir+"/"+args.box)
        outputname,ffDir = writeBashScript(args.tag, args.box, 
                args.model, mg, mchi, args.outDir, args.noSys, args.bkgDir, 
                args.saveWorkspace, noBoostCuts=args.noBoostCuts,
                fineGrained=args.fineGrained)
        
        os.system("mkdir -p "+pwd+"/"+ffDir)
        os.system("echo bsub -q "+args.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)        
        if not args.noSub:
            os.system("bsub -q "+args.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)

    print "nJobs = %i"%nJobs
