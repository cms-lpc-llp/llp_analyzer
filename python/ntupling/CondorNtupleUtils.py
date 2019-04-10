#!/bin/env python

### Sequence for data: submit --> hadd --> skim --> hadd-final --> remove-duplicates --> good-lumi --> copy-local
### Sequence for MC: submit --> hadd --> normalize --> hadd-final --> copy-local
### Sequence for fastsim signal samples: submit --> hadd

import math
import os,sys
import glob
import argparse
from subprocess import call, check_output

from ControlRegionNtuples2016_V3p15 import SAMPLES, TREETYPES, TREETYPEEXT, SKIMS, DIRS, OPTIONS, VERSION, DATA, SUFFIXES, ANALYZERS
import NtupleUtils as nt

RAZOR_EOS_DIR = '/eos/cms/store/group/phys_susy/razor/Run2Analysis/Analyzers/'

def getSubDirName(tag, datasetname):
    return os.path.join("condor", tag, datasetname)

def getSubFileName(subdir):
    return subdir+'/condorsubmit.cmd'

def getCondorSubmitFile(subdir, executable, flavor='espresso'):
    """
    Creates submission file and fills in header info common to all jobs.
    Returns the open file object.
    """
    subfile = getSubFileName(subdir)
    f = open(subfile,"w")
    f.write("Universe = vanilla\n")
    f.write("Executable = "+executable+"\n")
    f.write('Transfer_Output_Files = ""\n')
    f.write('Output = /dev/null\n')
    f.write('Error = /dev/null\n')
    f.write('Log = /dev/null\n')
    f.write("Notification = Never\n")
    f.write('+JobFlavour = "'+flavor+'"\n')
    return f

def writeCondorSubmitFragment(f, subdir, jobname, options):
    """
    Generates the part of the config file for the current job
    and appends it to the file f.
    """
    f.write("\nArguments = "+(' '.join(options))+"\n")
    # log files should be omitted to avoid job IO issues
    #f.write("Output = {}/{}.out\n".format(subdir, jobname))
    #f.write("Log = {}/{}.log\n".format(subdir, jobname))
    #f.write("Error = {}/{}.err\n".format(subdir, jobname))
    f.write("Queue\n")
    return f

def submitCondorJobsForSample(f, subdir, submit=False, spool=False,
        verbose=False):
    """
    Submits the condor jobs for the specified dataset
    """
    f.close()
    subfile = getSubFileName(subdir)
    cmd = "condor_submit -spool "+subfile
    if not spool:
        cmd = cmd.replace(' -spool', '')
    if verbose:
        print cmd
    if submit:
        os.system(cmd)

def submitJobs(analyzer,tag,isData=False,submit=False,label='',
        flavor='espresso', spool=False, verbose=False, fastsim=False):
    local_dir = os.environ['CMSSW_BASE']+'/src/RazorAnalyzer/'
    samples = SAMPLES
    listdir = local_dir+'lists/Run2/razorNtupler'+(VERSION.split('_')[0])+'/MC_Summer16'
    eos_list_dir = '{}/lists/{}/'.format(RAZOR_EOS_DIR, VERSION)
    jobssuffix = '/jobs'
    if isData:
        listdir = listdir.replace('/MC_Summer16','/data')
        samples = DATA
    elif fastsim:
        listdir = listdir.replace('/MC_Summer16', '/MCFastsim')
    filesperjob = 6
    script=local_dir+'scripts/runRazorJob_NoAFS.sh'
    call(['mkdir','-p',DIRS[tag]+'/jobs'])
    call(['mkdir','-p',eos_list_dir])
    # transfer needed files to EOS
    for f in [local_dir+'/RazorRun_NoAFS', 
              local_dir+'/bin/Run'+analyzer,
              local_dir+'/RazorRunAuxFiles_Expanded.tar.gz']:
        call(['cp', f, RAZOR_EOS_DIR])
    #samples loop
    for process in samples[tag]:
        for sample in samples[tag][process]:
            subdir = getSubDirName(tag, sample)
            call(['mkdir','-p',subdir])
            inlist = os.path.join(listdir,sample+'.cern.txt')
            if not os.path.isfile(inlist):
                print "Warning: list file",inlist,"not found!"
                continue
            nfiles = sum([1 for line in open(inlist)])
            maxjob = int(math.ceil( nfiles*1.0/filesperjob ))-1
            print "Sample:",sample
            call(['cp', inlist, eos_list_dir])
            inlist = eos_list_dir+'/'+os.path.basename(inlist)
            f = getCondorSubmitFile(subdir, script, flavor)
            njobs = 0
            for ijob in range(maxjob+1):
                outfile = nt.getJobFileName(analyzer,tag,sample,ijob,maxjob,label=label)
                if not os.path.isfile( DIRS[tag]+jobssuffix+'/'+outfile ):
                    njobs += 1
                    jobname = '_'.join([analyzer,sample,label,str(ijob)])
                    options = [analyzer,inlist,str(int(isData)),str(OPTIONS[tag]),
                            str(filesperjob),str(ijob),outfile,
                            DIRS[tag].replace('/eos/cms','')+jobssuffix,'CMSSW_8_0_26',
                            label]
                    if verbose:
                        print "Job %d of %d"%(ijob,maxjob)
                        print "Command: %s %s"%(script, ' '.join(options))
                    # remove absolute paths from argument list
                    options = [opt.replace(local_dir, '') for opt in options]
                    writeCondorSubmitFragment(f, subdir, jobname, options)
            if njobs > 0:
                print "Number of jobs: {} / {}".format(njobs, maxjob+1)
            submitCondorJobsForSample(f, subdir, submit=(submit and njobs), 
                    spool=spool, verbose=verbose)


if __name__ == '__main__':
    parser = nt.makeParser()
    parser.add_argument('--flavor', default='espresso', help=('Determines time allotted to job.'
        'Options are espresso, microcentury, longlunch, workday, tomorrow, testmatch, nextweek'))
    parser.add_argument('--spool', action='store_true', help=('Spool job input files to submit host.'
        '  Note that no log files are available in this case.'))
    args = parser.parse_args()
    tag = args.tag
    if args.fastsim:
        tag += 'Fastsim'
    analyzer = ANALYZERS[tag]
    isData = args.data
    noSub = args.noSub
    force = args.force
    skim = not args.noSkim
    label = args.label

    print "Analyzer:",analyzer
    print "Tag:",tag
    print "Label:",label

    if args.submit:
        print "Submit batch jobs..."
        submitJobs(analyzer,tag,isData,submit=(not noSub),label=label,
                flavor=args.flavor, spool=args.spool, fastsim=args.fastsim,
                verbose=args.verbose)

    if args.findZombies:
        print "Searching for zombie files..."
        nt.findZombies(analyzer,tag,isData,label=label)

    if args.hadd:
        print "Combine ntuples..."
        if args.fastsim:
            nt.haddFastsimJobs(analyzer,tag,label=label,dryRun=noSub)
        else:
            nt.haddFiles(analyzer,tag,isData,force,label=label)

    if args.normalize:
        print "Normalize ntuples..."
        if isData:
            sys.exit("Error: options --data and --normalize do not make sense together!")
        nt.normalizeFiles(analyzer,tag,force,label=label)

    if args.haddFinal:
        print "Combine normalized files..."
        if isData:
            nt.combineData(analyzer,tag,force,skim,label=label)
        else:
            nt.haddNormalizedFiles(analyzer,tag,force,label=label)

    if args.skim:
        print "Skim finished ntuples..."
        nt.skimNtuples(analyzer,tag,isData,label=label)

    if args.removeDuplicates:
        print "Remove duplicate events..."
        if not isData:
            print "--data was not specified, but I assume you want to use data."
        nt.removeDuplicates(analyzer,tag,skim,label=label)

    if args.goodLumi:
        print "Apply good lumi selection..."
        if not isData:
            print "--data was not specified, but I assume you want to use data."
        nt.goodLumi(analyzer,tag,skim,label=label)

    if args.copyLocal:
        print "Copy files locally..."
        nt.copyLocal(analyzer,tag,isData,skim,label=label)
