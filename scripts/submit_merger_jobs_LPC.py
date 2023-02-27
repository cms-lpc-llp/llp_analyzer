#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys
from collections import OrderedDict

queueType = "longlunch"
option = 0
label = "metoptimizationbkgstudiesfull"

analysis = "MergeNtuples"
outputfile = "MergeNtuples" + "_" + label

cmsswReleaseVersion = "CMSSW_10_6_8"
outputDirectoryBase = "/store/user/lpclonglived/guerrero/merger/"+analysis+"/"+label+"/"
taudirectory        = "root://cmsxrootd.fnal.gov//store/group/lpclonglived/guerrero/analyzer/VLLntupler/mergentuplesfull/taufiles/"

datasetListDir = "mergerVLLNtuple/"

datasetList = OrderedDict()

#signal samples (ready)
#datasetList['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.txt']   = [0, 1, "2018","","WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root"] 
datasetList['VLLPair_VLLToTauS_MVLL300_MS10_ctau300.txt']                     = [0, 10, "2018","","VLLPair_VLLToTauS_MVLL300_MS10_ctau300.root"]   #40
datasetList['VLLPair_VLLToTauS_MVLL300_MS10_ctau1000.txt']                    = [0, 10, "2018","","VLLPair_VLLToTauS_MVLL300_MS10_ctau1000.root"]  #12
datasetList['VLLPair_VLLToTauS_MVLL700_MS10_ctau100.txt']                     = [0, 10, "2018","","VLLPair_VLLToTauS_MVLL700_MS10_ctau100.root"]   #40
datasetList['VLLPair_VLLToTauS_MVLL700_MS10_ctau1000.txt']                    = [0, 10, "2018","","VLLPair_VLLToTauS_MVLL700_MS10_ctau1000.root"]  #7
datasetList['VLLPair_VLLToTauS_MVLL1000_MS10_ctau100.txt']                    = [0, 10, "2018","","VLLPair_VLLToTauS_MVLL1000_MS10_ctau100.root"]  #40
datasetList['VLLPair_VLLToTauS_MVLL1000_MS10_ctau1000.txt']                   = [0, 10, "2018","","VLLPair_VLLToTauS_MVLL1000_MS10_ctau1000.root"] #10

#2016 data (ready)
datasetList['Data2016_Run2016B-HighMET-07Aug17_ver1.txt']    = [1, 2, "2016","","MET_Run2016B-ver1-UL2016.root"] #12
datasetList['Data2016_Run2016B-HighMET-07Aug17_ver2.txt']    = [1, 2, "2016","","MET_Run2016B-ver2-UL2016.root"] #100
datasetList['Data2016_Run2016C-HighMET-07Aug17.txt']         = [1, 2, "2016","","MET_Run2016C-UL2016.root"] #32 
datasetList['Data2016_Run2016D-HighMET-07Aug17.txt']         = [1, 2, "2016","","MET_Run2016D-UL2016.root"] #51
datasetList['Data2016_Run2016E-HighMET-07Aug17.txt']         = [1, 2, "2016","","MET_Run2016E-UL2016.root"] #46
datasetList['Data2016_Run2016F-HighMET-07Aug17.txt']         = [1, 2, "2016","","MET_Run2016F-UL2016.root"] #33
datasetList['Data2016_Run2016G-HighMET-07Aug17.txt']         = [1, 2, "2016","","MET_Run2016G-UL2016.root"] #78
datasetList['Data2016_Run2016H-HighMET-07Aug17.txt']         = [1, 2, "2016","","MET_Run2016H-UL2016.root"] #88

#2017 data (ready)
datasetList['Data2017_Run2017B-HighMET-17Nov2017.txt']       = [1, 2, "2017","","MET_Run2017B-UL2017.root"] #45
datasetList['Data2017_Run2017C-HighMET-17Nov2017.txt']       = [1, 2, "2017","","MET_Run2017C-UL2017.root"] #97
datasetList['Data2017_Run2017D-HighMET-17Nov2017.txt']       = [1, 2, "2017","","MET_Run2017D-UL2017.root"] #48
datasetList['Data2017_Run2017E-HighMET-17Nov2017.txt']       = [1, 2, "2017","","MET_Run2017E-UL2017.root"] #76
datasetList['Data2017_Run2017F-HighMET-17Nov2017.txt']       = [1, 2, "2017","","MET_Run2017F-UL2017.root"] #103 

##2018 data (ready)
datasetList['Data2018_Run2018A-HighMET-17Sep2018.txt']              = [1, 2, "2018","","MET_Run2018A-UL2018.root"] #102
datasetList['Data2018_Run2018B-HighMET-17Sep2018.txt']              = [1, 2, "2018","","MET_Run2018B-UL2018.root"] #50
datasetList['Data2018_Run2018C-HighMET-17Sep2018.txt']              = [1, 2, "2018","","MET_Run2018C-UL2018.root"] #47
datasetList['Data2018_Run2018Dv1-PromptReco-HighMET-17Sep2018.txt'] = [1, 2, "2018","","MET_Run2018D-UL2018.root"] #2 
datasetList['Data2018_Run2018Dv2-PromptReco-HighMET-17Sep2018.txt'] = [1, 5, "2018","","MET_Run2018D-UL2018.root"] #95 

CMSSW_BASE_DIR = os.getenv('CMSSW_BASE')
Analyzer_DIR = CMSSW_BASE_DIR+"/src/llp_analyzer/"

#create directory for condor jobs
for listfile in datasetList.keys():
    datasetName = listfile.replace(".txt","")
    print "Preparing merger workflow for dataset :" + datasetName + "\n"
    if not os.path.exists(Analyzer_DIR+"/lists/" + datasetListDir + "/" + listfile):
        print "listfile: " + Analyzer_DIR+"/lists/" + datasetListDir + "/" + listfile + " does not exist. skipping."
        continue
    outputDirectory = outputDirectoryBase + datasetName + "/"
    tmpListFile     = open(Analyzer_DIR + "/lists/" + datasetListDir + "/" + listfile,"r")
    year       = datasetList[listfile][2]
    sampleName = datasetList[listfile][3]
    tauName    = taudirectory+datasetList[listfile][4]

    #####################################
    #Job Splitting
    #####################################
    isData = "no"
    if (datasetList[listfile][0] == 1): 
        isData = "yes"
    filesPerJob = datasetList[listfile][1]
    tmpJobFileCount = 0
    nJobs = 1
    if os.path.exists(Analyzer_DIR+"/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + "/" + "input_list_" + str(nJobs) + ".txt"):
        print "Warning: condor directory " + Analyzer_DIR + "/condor/analyzer_" + analysis + "_" + label + "/" + datasetName + " is not empty. Skipping."
        continue
    #create condor directories
    os.system("mkdir -p " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName )
    os.system("mkdir -p " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/log/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/out/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/err/")

    tmpOutputListFile = open( Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/" + "input_list_" + str(nJobs) + ".txt","w")
    for line in tmpListFile:        
        #open list file for new job
        if tmpJobFileCount >= filesPerJob:
            tmpOutputListFile.close()
            tmpJobFileCount = 0
            nJobs = nJobs + 1           
            tmpOutputListFile = open( Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/" + "input_list_" + str(nJobs) + ".txt","w")
        #write input file into job list file
        tmpOutputListFile.write(line)
        tmpJobFileCount += 1
    tmpOutputListFile.close()    
    os.system("cd " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/; tar czf input_list.tgz input_list_*.txt")

    #####################################
    #Copy run script and executable
    #####################################
    os.system("cp " + Analyzer_DIR + "/scripts/run_mergerjob_LPC.sh "   + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/")
    os.system("cp " + Analyzer_DIR + analysis    + " " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/")
    #####################################
    #Create Condor JDL file
    #####################################
    tmpCondorJDLFile = open(Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/task.jdl","w+")
    tmpCondorJDLFileTemplate = """
Universe  = vanilla
Executable = ./run_mergerjob_LPC.sh
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    tmpCondorJDLFile.write("Arguments = " + analysis + " " + str(isData) + " " + str(option) + " " + "$(I) " + outputfile + " " + outputDirectory + " " + cmsswReleaseVersion + " " + year + " " + sampleName + " " + tauName + "\n")
    tmpCondorJDLFileTemplate = """
Log = log/job.$(Cluster).$(Process).log
Output = out/job.$(Cluster).$(Process).out
Error = err/job.$(Cluster).$(Process).err
x509userproxy = $ENV(X509_USER_PROXY)
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    tmpCondorJDLFile.write("transfer_input_files = " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/run_mergerjob_LPC.sh, " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/input_list.tgz, " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/" + analysis + "\n")
    tmpCondorJDLFileTemplate = """
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
# Resources request
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    tmpCondorJDLFile.write("RequestMemory = 6000 \n") #MB
    tmpCondorJDLFileTemplate = """
# Jobs selection
Queue I from (
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    for i in range(1,nJobs+1):
        tmpCondorJDLFile.write(str(i)+"\n")
    tmpCondorJDLFile.write(")\n")
    tmpCondorJDLFile.close()

