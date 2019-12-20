#! /usr/bin/env python
import os
import sys
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) < 10:
    print "usage SubmitT3-higgsCondor.sh ProcessName sampleList.list files_per_job hadoopDirName tarball pick_events.txt treeName runBr\
anch eventBranch"
    print "e.g: python SubmitT3-higgsCondor.sh ProcessName list.txt  files_per_job hadoopDirName tarball  pick_events.txt treeName runBranch eventBranch"
    sys.exit(1)

tarball = "" 

os.system("export PATH=$PATH:/opt/gridengine/bin/lx26-amd64/")

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################

############################
#configuring Submission
############################
process          = sys.argv[1]
inputlist        =  sys.argv[2]
files_per_job    = int(sys.argv[3])
hadoopDir        = sys.argv[4]
tarball          = sys.argv[5]
pick_events      = sys.argv[6]
treeName         = sys.argv[7]
runBranch        = sys.argv[8]
eventBranch      = sys.argv[9]

if pick_events.find( "/home/cmorgoth/" ) == -1:
    pick_events = pwd + "/" + pick_events
#######################################
#Splitting jobs
#######################################
files_in_last_job = 0
num_lines = sum(1 for line in open(inputlist))
if num_lines%files_per_job == 0:
    njobs = num_lines/files_per_job
else:
    njobs = num_lines/files_per_job + 1
    files_in_last_job = num_lines%files_per_job
input = open(inputlist)
######################################
print "[INFO]: Creating " + str( njobs ) + " jobs"
################################################
os.system("mkdir -p submission")
os.system("mkdir -p submission/"+process)
os.system("mkdir -p submission/"+process+"/log/")
os.system("mkdir -p submission/"+process+"/input/")
os.system("mkdir -p submission/"+process+"/src/")
os.system("mkdir -p submission/"+process+"/out/")
#######################################

submitjobs_name = "submission/"+process+"/src/send_allfiles.sh"
submitjobs = open(submitjobs_name, 'w')
submitjobs.write("#!/bin/bash \n\n") 

for ijob in range(njobs):
    # prepare the list file
    inputfilename = pwd+"/submission/"+process+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    # if it is a normal job get filesperjob lines
    if ijob != (njobs-1):
        for line in range( files_per_job ):
            ntpfile = input.readline() 
            inputfile.write( ntpfile )
            continue
    else:
        # if it is the last job get ALL remaining lines
        ntpfile = input.readline()
        while ntpfile != '':
            inputfile.write( ntpfile )
            ntpfile = input.readline()
            continue
    inputfile.close()

    outputName = "submission/"+process+"/src/submit_"+str(ijob)+".sh"
    logName = "submission/"+process+"/log/submit_"+str(ijob)
    submitName = "submission/"+process+"/src/condor_"+str(ijob)+".sub"
    rootOutput = "skim_" + treeName + "_" + str(ijob)+ ".root"
    #print "OUTPUTFILE: ", outputName
    submitFile = open(submitName,'w')
    submitFile.write("Executable = "+pwd+"/"+outputName)
    submitFile.write("\nUniverse = vanilla")
    submitFile.write("\nOutput = "+pwd+"/"+logName+".out.$(Process)")
    submitFile.write("\nLog = "+pwd+"/"+logName+".log")
    submitFile.write("\nError = "+pwd+"/"+logName+".err")
    submitFile.write("\ngetenv = True\n\nQueue")
    
    outputfile = open(outputName,'w')
    outputfile.write("#!/bin/sh\n\n")
    
    outputfile.write("export SCRAM_ARCH=slc5_amd64_gcc472\n")
    outputfile.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
    outputfile.write("mkdir "+process+"_"+str(ijob))
    outputfile.write("\ncd "+process+"_"+str(ijob))
    outputfile.write("\ncmsrel CMSSW_6_2_0;cd CMSSW_6_2_0/src;cmsenv\n")
    outputfile.write("tar -xf "+tarball+"; cd SkimEvents\n")
    outputfile.write("make clean;make -j 4\n")
    outputfile.write("cp " + inputfilename + " input_list.txt\n")
    outputfile.write("cp " + pick_events + " pick_events.txt\n")
    outputfile.write( "\n./CreateSkimmed --event_list==pick_events.txt"
                      + " --list_of_ntuples==input_list.txt"
                      + " --output_name==" + rootOutput
                      + " --tree_name==" + treeName
                      + " --run_branch==" + runBranch
                      + " --event_branch==" + eventBranch
                      )
    outputfile.write("\nmkdir -p /mnt/hadoop/store/user/cmorgoth/"+hadoopDir)
    outputfile.write("\nmkdir -p /mnt/hadoop/store/user/cmorgoth/"+hadoopDir+"/"+process)
    outputfile.write("\ncp " + rootOutput + " /mnt/hadoop/store/user/cmorgoth/"+hadoopDir+"/"+process+"/.")
    
    outputfile.close()
    os.system("sleep .1\n")
    scriptfile=pwd+"/"+outputName
    submitjobs.write("condor_submit "+pwd+"/"+submitName)
    submitjobs.write("\nsleep 0.2\n")
    ijob = ijob+1
    continue

submitjobs.close()

