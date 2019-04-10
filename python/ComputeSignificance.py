### This script takes an existing combine card and 
### computes the signal significance instead of the limit.

import os
import subprocess as sp
import glob
import argparse
import ROOT as rt

from GChiPairs import gchipairs
from limits.SMSConfig import VERSION, sms_models
from limits.SMSConfig import BOOST_LOCAL_DIR, BOOST_LIMIT_DIR, BOOST_BOXES

def writeBashScript(boxes, model, mg, mchi, inDir, subDir, with_boost):
    cmsswBase = os.environ['CMSSW_BASE']
    massPoint = "%i_%i"%(mg, mchi)
    name = "SMS-{}_{}".format(model, massPoint)
    outputname = subDir+"/submit_signif_{}.src".format(name)
    if with_boost:
        outputname = outputname.replace('signif', 'signif_boostcombine')

    card = 'RazorInclusiveMADD_{}_{}.txt'.format(name, '_'.join(boxes))
    inclusiveInputFiles = ['RazorInclusiveMADD_SMS-{}_{}_{}.root'.format(
        model, massPoint, box) for box in boxes]
    nameForCombine = 'MADDSignificance'

    if with_boost:
        # boost datacards use 0 instead of 1 for lowest LSP mass
        boostName = name
        boostMassPoint = massPoint
        if mchi == 1:
            boostMassPoint = "%i_0"%(mg)
            boostName = "SMS-{}_{}".format(model, boostMassPoint)
        boostDir = '{}/{}'.format(BOOST_LIMIT_DIR, model)
        boostInputFiles = ['RazorBoost_SMS-{}_{}_{}.root'.format(
            model, boostMassPoint, box) for box in BOOST_BOXES]

        if len(boxes) == 1:
            # This is to allow combining razor boost with individual 
            # razor inclusive boxes without breaking the naming convention
            # from before (which by default doesn't use the box name).
            name = boxes[0]+'_'+name
        nameForCombine = nameForCombine.replace('MADD', 'RazorInclusiveBoost')
        card = "RazorInclusiveBoostCombined_{}.txt".format(name)
    outResult = (
        "higgsCombineMADDSignificance_{}.ProfileLikelihood.mH120.root".format(name))
    if with_boost:
        outResult = outResult.replace('MADD', 'RazorInclusiveBoost')

    cardFile = os.path.join(inDir, card)
    if not os.path.isfile(cardFile):
        print "Card {} not found!".format(cardFile)
        return None
    print "Input card: {}".format(cardFile)

    script =  '#!/usr/bin/env bash -x\n'
    script += 'pwd\n'
    script += 'cd %s/src/RazorAnalyzer \n'%(cmsswBase)
    script += 'pwd\n'
    script += "export SCRAM_ARCH=slc6_amd64_gcc481\n"
    script += 'eval `scramv1 runtime -sh`\n'
    script += 'cd - \n'

    script += 'mkdir -p %s\n'%(VERSION)
    script += 'cp %s/%s .\n'%(inDir, card)
    for f in inclusiveInputFiles:
        script += 'cp %s/%s %s\n'%(inDir, f, VERSION)
    if with_boost:
        script += 'mkdir -p %s\n'%(BOOST_LOCAL_DIR)
        for f in boostInputFiles:
            script += 'cp %s/%s %s\n'%(boostDir, f, BOOST_LOCAL_DIR)

    script += 'combine -M ProfileLikelihood --signif %s -n %s\n'%(card, nameForCombine+'_'+name)
    script += 'cp %s %s/\n'%(outResult, inDir) 
    
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
    parser.add_argument('--combined-with-boost', action='store_true',
              help='whether to combine with boosted analysis')

    args = parser.parse_args()

    nJobs = 0
    donePairs = []
    if args.doneFile is not None:
        combineMethod = 'ProfileLikelihood'
        with open(args.doneFile,'r') as f:            
            allFiles = [line.replace('\n','') for line in f.readlines()]
            for (mg, mchi) in gchipairs(args.model):
                outputname = 'higgsCombineMADDSignificance_SMS-%s_%i_%i.%s.mH120.root'%(
                        args.model, mg, mchi, combineMethod)
                if args.combined_with_boost:
                    outputname = outputname.replace('MADD', 'RazorInclusiveBoost')
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
                mg, mchi, args.outDir, subDir, args.combined_with_boost)
        if outputname is None:
            continue
        
        cmd = ['bsub', '-q', args.queue, '-o', 
                subDir+'/log_signif_{}_{}_{}.log'.format(
            args.model, mg, mchi), 'source', outputname]
        print ' '.join(cmd)
        if not args.noSub:
            sp.call(cmd)

    print "nJobs = %i"%nJobs
