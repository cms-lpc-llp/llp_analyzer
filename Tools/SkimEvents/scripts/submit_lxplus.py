#! /usr/bin/env python
import os
import sys
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) < 5:
    print "usage Submit-lxplus.sh ProcessName sampleList.list files_per_job pick_events.txt treeName runBranch eventBranch"
    print "e.g: python Submit-lxplus.sh ProcessName list.txt  files_per_job pick_events.txt treeName runBranch eventBranch"
    sys.exit(1)


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
pick_events      = sys.argv[4]
#tree name
if len(sys.argv) > 5 :
    treeName = sys.argv[5]
else:
    treeName = "ntp1"
#run branch name
if len(sys.argv) > 6:
    runBranch = sys.argv[6]
else:
    runBranch = "runNumber"
    
#event branch name
if len(sys.argv) > 7:
    eventBranch = sys.argv[7]
else:
    eventBranch = "eventNumber"

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

    outputName = pwd + "/submission/" + process + "/src/submit_" + str(ijob) + ".sh"
    rootOutput = pwd + "/submission/" + process + "/out/" + process + "_" + str(ijob) + ".root"
    logName = pwd + "/submission/" + process + "/log/log_" + str(ijob)
    errName = pwd + "/submission/" + process + "/log/err_" + str(ijob)
    outputfile = open(outputName,'w')
    outputfile.write("#!/bin/sh\n\n")
    
    outputfile.write( "cd " + pwd )
    outputfile.write( "\neval `scramv1 run -sh`;")
    outputfile.write( "\n./CreateSkimmed --event_list==" + pick_events  
                      + " --list_of_ntuples==" + inputfilename 
                      + " --output_name==" + rootOutput
                      + " --tree_name==" + treeName
                      + " --run_branch==" + runBranch
                      + " --event_branch==" + eventBranch
                      )
        
    outputfile.close
    os.system("bsub -q 1nh source " + outputName + " -o " +logName+ " -e " +errName)
    os.system("sleep .1\n")
    scriptfile=pwd+"/"+outputName
    ijob = ijob+1
    continue


