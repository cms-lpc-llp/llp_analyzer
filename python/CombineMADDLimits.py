import os
import subprocess as sp
import glob
import argparse
import ROOT as rt

from GChiPairs import gchipairs
from limits.SMSConfig import VERSION, sms_models
    
def writeBashScript(boxes, model, mg, mchi, inDir, subDir):
    
    massPoint = "%i_%i"%(mg, mchi)
    name = "SMS-{}_{}_{}".format(model, massPoint, '_'.join(boxes))
    alt_name = "MADD_{}_SMS-{}_{}".format('_'.join(boxes), model, massPoint)
    outputname = subDir+"/submit_combine_{}.src".format(name)
    inputCards = ['RazorInclusiveMADD_SMS-{}_{}_{}.txt'.format(
            model, massPoint, box) for box in boxes]
    for card in inputCards:
        cardFile = os.path.join(inDir, card)
        if not os.path.isfile(cardFile):
            print "Card {} not found!".format(cardFile)
            return None
        print "Input card: {}".format(cardFile)
    inputFiles = [f.replace('.txt', '.root') for f in inputCards]
    cmsswBase = os.environ['CMSSW_BASE']
    combinedName = "RazorInclusiveMADD_{}.txt".format(name)
    combinedResult = "higgsCombine{}.Asymptotic.mH120.root".format(alt_name)

    script =  '#!/usr/bin/env bash -x\n'
    script += 'pwd\n'
    script += 'cd %s/src/RazorAnalyzer \n'%(cmsswBase)
    script += 'pwd\n'
    script += "export SCRAM_ARCH=slc6_amd64_gcc481\n"
    script += 'eval `scramv1 runtime -sh`\n'
    script += 'cd - \n'
    script += 'mkdir %s\n'%(VERSION)
    for f in inputCards:
        script += 'cp %s/%s .\n'%(inDir, f)
    for f in inputFiles:
        script += 'cp %s/%s %s\n'%(inDir, f, VERSION)
    script += 'combineCards.py %s > %s\n'%(' '.join(inputCards), combinedName)
    script += 'combine -M Asymptotic %s -n %s\n'%(combinedName, alt_name)
    script += 'cp %s %s/\n'%(combinedName, inDir) 
    script += 'cp %s %s/\n'%(combinedResult, inDir) 
    
    outputfile = open(outputname,'w')
    outputfile.write(script)
    outputfile.close

    return outputname


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tag',default='Razor2016_MoriondRereco')
    parser.add_argument('--boxes', nargs='+', required=True)
    parser.add_argument('--model',help="signal model name",required=True)
    parser.add_argument('--dir',dest="outDir",default="./",
              help="Output directory to store cards")
    parser.add_argument('--no-sub',dest="noSub",action='store_true',
              help="no submission")
    parser.add_argument('--queue',default="1nd",
              help="queue: 1nh, 8nh, 1nd, etc.")
    parser.add_argument('--done-file',dest="doneFile",
              help="file containing output file names")

    args = parser.parse_args()

    nJobs = 0
    donePairs = []
    if args.doneFile is not None:
        combineMethod = 'Asymptotic'
        with open(args.doneFile,'r') as f:            
            allFiles = [line.replace('\n','') for line in f.readlines()]
            for (mg, mchi) in gchipairs(args.model):
                outputname = 'higgsCombineMADD_%s_SMS-%s_%i_%i.%s.mH120.root'%(
                        '_'.join(args.boxes),args.model,mg,mchi,combineMethod)
                if outputname in allFiles: donePairs.append((mg,mchi))

    submodels = sms_models[args.model].submodels
    if submodels is None:
        pairs = gchipairs(args.model)
    else:
        pairs = []
        for submodel in submodels:
            pairs += gchipairs(submodel)
    for (mg, mchi) in pairs:
        if (mg, mchi) in donePairs: 
            print (mg,mchi),"is already done; skipping"
            continue
        nJobs+=1

        subDir = '{}/Limits/{}'.format(os.environ['PWD'], VERSION)
        sp.call(['mkdir', '-p', subDir])
        outputname = writeBashScript(args.boxes, args.model,
                mg, mchi, args.outDir, subDir)
        if outputname is None:
            continue
        
        cmd = ['bsub', '-q', args.queue, '-o', subDir+'/log_{}_{}_{}.log'.format(
            args.model, mg, mchi), 'source', outputname]
        print ' '.join(cmd)
        if not args.noSub:
            sp.call(cmd)

    print "nJobs = %i"%nJobs
