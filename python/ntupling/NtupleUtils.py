#!/bin/env python

### Sequence for data: submit --> find-zombies --> hadd --> skim --> hadd-final --> remove-duplicates --> good-lumi --> copy-local
### Sequence for MC: submit --> find-zombies --> hadd --> normalize --> hadd-final --> copy-local
### Sequence for fastsim signal samples: submit --> hadd

# Note: finding zombie files is prohibitively slow for fastsim.
# If using lxbatch for submission, make sure to use a long enough queue (e.g. 8nh)
# so that jobs do not time out and create zombie files.
# Even better, find a way to stop these zombies from being created when jobs time out.

import math
import os,sys
import glob
import argparse
from subprocess import call, check_output

sys.path.append(os.path.dirname(__file__)+'/..')
from haddFastsimSMS import makeFileLists, haddFastsimFiles
from ControlRegionNtuples2016_V3p15 import SAMPLES, TREETYPES, TREETYPEEXT, SKIMS, DIRS, OPTIONS, VERSION, DATA, SUFFIXES, ANALYZERS, EXTRASKIMS

RAZOR_EOS_DIR = '/eos/cms/store/group/phys_susy/razor/Run2Analysis/Analyzers/'

def getSamplePrefix(analyzer,tag,reHLT=False,label=''):
    return analyzer.replace('RazorControl','RunTwoRazorControl')+(
            (TREETYPEEXT[tag]!='')*('_'+TREETYPEEXT[tag]))+(
            (SKIMS[tag]!='')*('_'+SKIMS[tag]))+(
            (reHLT==True)*('_reHLT'))+(
            (label != '')*('_'+label))

def getJobFileName(analyzer,tag,sample,ijob,maxjob,reHLT=False,label=''):
    prefix = getSamplePrefix(analyzer,tag,reHLT,label)
    return '%s_%s.Job%dof%d.root'%(prefix,sample,ijob,maxjob)

def getFileName(analyzer,tag,sample,reHLT=False,label=''):
    prefix = getSamplePrefix(analyzer,tag,reHLT,label)
    return '%s_%s.root'%(prefix,sample)

def submitJobs(analyzer,tag,isData=False,submit=False,reHLT=False,label='',
        queue='8nm', filesperjob=6, fastsim=False, 
        jobsLimit=-1, verbose=False):
    # parameters
    local_dir = os.environ['CMSSW_BASE']+'/src/RazorAnalyzer/'
    samples = SAMPLES
    listdir = local_dir+'lists/Run2/razorNtupler'+(VERSION.split('_')[0])+'/MC_Summer16'
    eos_list_dir = '{}/lists/{}/'.format(RAZOR_EOS_DIR, VERSION)
    if reHLT: listdir += 'reHLT'
    jobssuffix = '/jobs'
    if isData:
        listdir = listdir.replace('/MC_Summer16','/data')
        samples = DATA
    elif fastsim:
        listdir = listdir.replace('/MC_Summer16','/MCFastsim')
    script=local_dir+'scripts/runRazorJob_NoAFS.sh'
    call(['mkdir','-p',DIRS[tag]+'/jobs'])
    call(['mkdir','-p',eos_list_dir])
    # transfer needed files to EOS
    for f in [local_dir+'/RazorRun_NoAFS', 
              local_dir+'/bin/Run'+analyzer,
              local_dir+'/RazorRunAuxFiles_Expanded.tar.gz']:
        call(['cp', f, RAZOR_EOS_DIR])
    # LSF specific configuration
    os.environ['LSB_JOB_REPORT_MAIL'] = 'N'
    if not os.path.isdir(local_dir+'/output/'):
        os.mkdir(local_dir+'/output/')
    #samples loop
    totalJobs = 0
    for process in samples[tag]:
        for sample in samples[tag][process]:
            inlist = os.path.join(listdir,sample+'.cern.txt')
            if not os.path.isfile(inlist):
                print "Warning: list file",inlist,"not found!"
                continue
            nfiles = sum([1 for line in open(inlist)])
            maxjob = int(math.ceil( nfiles*1.0/filesperjob ))-1
            print "Sample:",sample
            call(['cp', inlist, eos_list_dir])
            inlist = eos_list_dir+'/'+os.path.basename(inlist)
            njobs = 0
            for ijob in range(maxjob+1):
                outfile = getJobFileName(analyzer,tag,sample,ijob,maxjob,reHLT,label)
                if not os.path.isfile( DIRS[tag]+jobssuffix+'/'+outfile ):
                    njobs += 1
                    totalJobs += 1
                    jobname = '_'.join([analyzer,sample,label,str(ijob)])
                    logfile = '/dev/null'
                    #logfile = os.path.join(local_dir,'output','%s_%s_%s_%d.out'%(
                            #analyzer,sample,label,ijob))
                    cmd = ['bsub','-q',queue,'-o',logfile,'-J',jobname,script,analyzer,inlist,
                            str(int(isData)),str(OPTIONS[tag]),str(filesperjob),str(ijob),outfile,
                        DIRS[tag].replace('/eos/cms','')+jobssuffix, 'CMSSW_8_0_26',
                        label]
                    if verbose:
                        print "Job %d of %d"%(ijob,maxjob)
                        print ' '.join(cmd)
                    if submit:
                        call(cmd)
                    if jobsLimit >= 0 and totalJobs >= jobsLimit:
                        print "Submitted total of {} jobs; stopping".format(totalJobs)
                        return
            if njobs > 0:
                print "Number of jobs: {} / {}".format(njobs, maxjob+1)

def findZombies(analyzer,tag,isData=False,reHLT=False,label='',fastsim=False):
    """Looks through job files and searches for Zombies.  Prints out a list of bad files."""
    import ROOT as rt
    samples = SAMPLES
    if isData:
        samples = DATA
    zombieFileName = "Zombies_%s_%s.txt"%(tag, label)
    if isData:
        zombieFileName = zombieFileName.replace(".txt","_Data.txt")
    # get all file names ahead of time to avoid repeated expensive ls calls
    allFiles = glob.glob(DIRS[tag]+'/jobs/*.Job*.root')
    allFilesDict = {}
    for f in allFiles:
        sample = os.path.basename(f.split('.')[0])
        if sample not in allFilesDict:
            allFilesDict[sample] = []
        allFilesDict[sample].append(f)
    with open(zombieFileName,'w') as zombieFile:
        for process in samples[tag]:
            for sample in samples[tag][process]:
                print "Sample:",sample
                matchStr = getFileName(analyzer, tag, sample, reHLT, label).replace('.root', '')
                try:
                    jobfiles = allFilesDict[matchStr]
                except KeyError:
                    print "Warning: no files found for {}".format(sample)
                    continue
                if len(jobfiles) > 0:
                    for f in jobfiles:
                        if os.path.getsize(f) < 1000:
                            if fastsim and '_' not in f.split('.')[1]:
                                print "skipping",f
                                continue # this is the empty placeholder file for this signal sample
                            print "\nNO KEYS:",f,"\n"
                            zombieFile.write(f+"\n")
                            continue
                        curFile = rt.TFile.Open(f)
                        if not curFile:
                            print "\nZOMBIE:",f,"\n"
                            zombieFile.write(f+"\n")
                else:
                    print "Warning: no files found (",matchStr,")"

def haddFiles(analyzer,tag,isData=False,force=False,reHLT=False,label=''):
    samples = SAMPLES
    if isData:
        samples = DATA
    for process in samples[tag]:
        for sample in samples[tag][process]:
            print "Sample:",sample
            fname = DIRS[tag]+'/'+getFileName(analyzer,tag,sample,reHLT,label)
            query = DIRS[tag]+'/jobs/'+(getFileName(analyzer,tag,sample,reHLT,label).replace('.root','*.Job*.root'))
            jobfiles = glob.glob( query )
            if os.path.isfile( fname ) and not force:
                print "File",fname,"exists; skipping"
            elif len(jobfiles) > 0:
                if sample in EXTRASKIMS:
                    skim_str = EXTRASKIMS[sample]
                    print "Performing additional skim {} before hadd".format(skim_str)
                    skimfilename = 'skim_{}.txt'.format(sample)
                    with open(skimfilename, 'w') as skimfile:
                        for f in jobfiles:
                            skimfile.write(f+'\n')
                    call(['./SkimNtuple', skimfilename, DIRS[tag]+'/jobs', 'ExtraSkim', skim_str])
                    jobfiles = [f.replace('.root', '_ExtraSkim.root') for f in jobfiles]
                if force:
                    call(['hadd','-f',fname]+jobfiles)
                else:
                    call(['hadd',fname]+jobfiles)
            else:
                print "Warning: no files found (",query,")"

def haddFastsimJobs(analyzer,tag,label='',dryRun=True):
    for process in SAMPLES[tag]:
        fileLists = makeFileLists(DIRS[tag]+'/jobs/', process)
        outDir = DIRS[tag]
        haddFastsimFiles(fileLists, process, DIRS[tag]+'/jobs/', outDir, dryRun=dryRun)

def normalizeFiles(analyzer,tag,force=False,reHLT=False,label=''):
    #make list file for normalizing
    paths = glob.glob( DIRS[tag]+'/*.root' )
    with open('ntuples_'+tag+'.txt','w') as normfile:
        for f in paths:
            #check if normalized file exists
            if (not force) and os.path.isfile( f.replace('.root','_1pb_weighted.root') ): continue
            sample = os.path.basename(f).replace('.root','').replace(
                    getSamplePrefix(analyzer,tag,reHLT,label)+'_','')
            #check if we need this sample
            for process in SAMPLES[tag]:
                if sample in SAMPLES[tag][process]:
                    normfile.write(sample+' '+f+'\n')
                    break
    call(['./NormalizeNtuple','ntuples_'+tag+'.txt'])

def haddNormalizedFiles(analyzer,tag,force=False,reHLT=False,label=''):
    for process in SAMPLES[tag]:
        print "Process:",process
        haddList = []
        for sample in SAMPLES[tag][process]:
            print "Sample:",sample,
            fname = DIRS[tag]+'/'+getFileName(analyzer,tag,sample,reHLT,label).replace(
                    '.root','_1pb_weighted.root')
            if os.path.isfile( fname ):
                print "found!"
                haddList.append( fname )
            else:
                print "not found..."
        outName = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag,reHLT,label)+'_'+process+'_1pb_weighted.root'
        if len(haddList) > 0:
            if force:
                call(['hadd','-f',outName]+haddList)
            else:
                call(['hadd',outName]+haddList)
        else:
            print "Skipping process",process,"because no input files were found!"
            print "Failed to create file",outName

def combineData(analyzer,tag,force=False,skim=True,label=''):
    haddList = []
    for process in DATA[tag]:
        print "Dataset:",process
        for sample in DATA[tag][process]:
            fname = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag,label=label)+'_'+sample+'.root'
            if skim:
                fname = fname.replace('.root','_RazorSkim.root')
            print fname,
            if os.path.isfile( fname ):
                print "found!"
                haddList.append( fname )
            else:
                print "not found..."
    outName = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag,label=label)+'_Data.root'
    if skim:
        outName = outName.replace('.root','_RazorSkim.root')
    if len(haddList) > 0:
        if force:
            call(['hadd','-f',outName]+haddList)
        else:
            call(['hadd',outName]+haddList)
    else:
        print "No input files were found!"
        return

def skimNtuples(analyzer,tag,isData=False,label=''):
    samples = SAMPLES
    if isData:
        samples = DATA
    #make list file for skimming
    with open('skim_'+tag+'.txt','w') as skimfile:
        for process in samples[tag]:
            for sample in samples[tag][process]:
                fname = DIRS[tag]+'/'+getFileName(analyzer,tag,sample,label=label)
                if not isData: fname = fname.replace('.root','_1pb_weighted.root')
                if os.path.isfile( fname ):
                    print "Adding",sample,"to skim file"
                    skimfile.write(fname+'\n')
                else:
                    print "Input file for",sample,"not found!"
                    print "( looking for",fname,")"
    skimString = 'MR%s > 300 && Rsq%s > 0.15'%(SUFFIXES[tag],SUFFIXES[tag])
    print "Skimming with",skimString
    call(['./SkimNtuple','skim_'+tag+'.txt',DIRS[tag],'RazorSkim',skimString])

def removeDuplicates(analyzer,tag,skim=True,label=''):
    outName = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag,label=label)+'_Data.root'
    if skim:
        outName = outName.replace('.root','_RazorSkim.root')
    outNameNoDuplicates = outName.replace('_Data','_Data_NoDuplicates')
    dupMacro = 'macros/RemoveDuplicateEvents.C'
    macroCall = "root -l '%s(\"%s\",\"%s\")'"%(dupMacro,outName,outNameNoDuplicates)
    print macroCall
    call(macroCall,shell=True)

def goodLumi(analyzer,tag,skim=True,label=''):
    #get location of good lumi script and python file containing JSON link
    script = check_output(['which', 'FWLiteGoodLumi'])
    pythonFile = os.path.join(os.environ['CMSSW_BASE'],'src','RazorCommon','Tools','python','loadJson.py')
    #copy file locally
    inName = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag,label=label)+'_Data_NoDuplicates.root'
    if skim:
        inName = inName.replace('.root','_RazorSkim.root')
    outName = inName.replace('.root','_GoodLumiGolden.root')
    #call script
    print "FWLiteGoodLumi %s %s %s"%(pythonFile,inName,outName)
    call(['FWLiteGoodLumi',pythonFile,inName,outName])

def copyLocal(analyzer,tag,isData=False,skim=True,label=''):
    samples = SAMPLES
    if isData:
        samples = DATA
    #make directory
    localdir = 'Backgrounds/'+tag
    call(['mkdir','-p',localdir])
    if isData:
        fname = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag,label=label)+'_Data_NoDuplicates_GoodLumiGolden.root'
        if skim:
            fname = fname.replace('GoodLumiGolden.root','RazorSkim_GoodLumiGolden.root')
        print "cp",fname,localdir
        call(['cp',fname,localdir])
    else:
        for process in samples[tag]:
            fname = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag,label=label)+'_'+process+'_1pb_weighted.root'
            print "cp",fname,localdir
            call(['cp',fname,localdir])

def makeParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('tag',help='1L, 2L, ...')
    parser.add_argument('--data',action='store_true',help='Run on data (MC otherwise)')
    parser.add_argument('--submit',action='store_true',help='Submit batch jobs')
    parser.add_argument('--no-sub', dest='noSub',action='store_true', 
            help='Print commands but do not submit')
    parser.add_argument('--force', action='store_true', help='Force HADD to recreate files')
    parser.add_argument('--hadd',action='store_true',help='Combine ntuple files')
    parser.add_argument('--normalize',action='store_true',help='Normalize ntuple files')
    parser.add_argument('--hadd-final',dest='haddFinal',action='store_true',help='Combine normalized ntuple files')
    parser.add_argument('--skim',action='store_true',help='Apply razor skim to final ntuples')
    parser.add_argument('--remove-duplicates',dest='removeDuplicates',action='store_true',help='Remove duplicates')
    parser.add_argument('--good-lumi',dest='goodLumi',action='store_true',help='Apply good lumi selection')
    parser.add_argument('--find-zombies',dest='findZombies',action='store_true',help='Find zombie files')
    parser.add_argument('--copy-local',dest='copyLocal',action='store_true',help='Copy files locally')
    parser.add_argument('--reHLT',action='store_true',help='Process reHLT samples')
    parser.add_argument('--no-skim',dest='noSkim',action='store_true',help='Do not assume skimmed data')
    parser.add_argument('--fastsim', action='store_true', help='process fastsim SMS ntuples')
    parser.add_argument('--label', help='label for RazorRun',
            default='Razor2016_MoriondRereco') 
    parser.add_argument('--verbose', action='store_true', help='Print job commands')
    return parser


if __name__ == '__main__':
    parser = makeParser()
    parser.add_argument('--queue', default='8nm', help='batch queue')
    parser.add_argument('--files-per-job', default=6, type=int)
    parser.add_argument('--limit', default=-1, type=int,
            help='Max number of jobs to submit')
    args = parser.parse_args()
    tag = args.tag
    if args.fastsim:
        tag += 'Fastsim'
    analyzer = ANALYZERS[tag]
    isData = args.data
    noSub = args.noSub
    force = args.force
    reHLT = args.reHLT
    skim = not args.noSkim
    label = args.label

    print "Analyzer:",analyzer
    print "Tag:",tag

    if args.submit:
        print "Submit batch jobs..."
        submitJobs(analyzer,tag,isData,submit=(not noSub),reHLT=reHLT,fastsim=args.fastsim,
                label=label,queue=args.queue,filesperjob=args.files_per_job, 
                jobsLimit=args.limit, verbose=args.verbose)

    if args.findZombies:
        print "Searching for zombie files..."
        findZombies(analyzer,tag,isData,reHLT=reHLT,label=label,fastsim=args.fastsim)

    if args.hadd:
        print "Combine ntuples..."
        if args.fastsim:
            haddFastsimJobs(analyzer,tag,label=label,dryRun=noSub)
        else:
            haddFiles(analyzer,tag,isData,force,reHLT=reHLT,label=label)

    if args.normalize:
        print "Normalize ntuples..."
        if isData:
            sys.exit("Error: options --data and --normalize do not make sense together!")
        normalizeFiles(analyzer,tag,force,reHLT=reHLT,label=label)

    if args.haddFinal:
        print "Combine normalized files..."
        if isData:
            combineData(analyzer,tag,force,skim,label=label)
        else:
            haddNormalizedFiles(analyzer,tag,force,reHLT=reHLT,label=label)

    if args.skim:
        print "Skim finished ntuples..."
        skimNtuples(analyzer,tag,isData,label=label)

    if args.removeDuplicates:
        print "Remove duplicate events..."
        if not isData:
            print "--data was not specified, but I assume you want to use data."
        removeDuplicates(analyzer,tag,skim,label=label)

    if args.goodLumi:
        print "Apply good lumi selection..."
        if not isData:
            print "--data was not specified, but I assume you want to use data."
        goodLumi(analyzer,tag,skim,label=label)

    if args.copyLocal:
        print "Copy files locally..."
        copyLocal(analyzer,tag,isData,skim,label=label)
