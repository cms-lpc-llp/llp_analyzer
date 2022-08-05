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
label = "option0"

analysis = "RunDataRunEventIndexing"
outputfile = "RunDataRunEventIndexing" + "_" + label


cmsswReleaseVersion = "CMSSW_10_6_8"
outputDirectoryBase = "/store/user/lpclonglived/sixie/analyzer/"+analysis+"/"+label+"/"
#datasetListDir = "displacedJetMuonNtuple/V1p171/Data2018_UL/MuonHitsOnly/"
#datasetListDir = "displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/v5/sixie/SingleMuon/"
#datasetListDir = "displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/"
datasetListDir = "displacedJetMuonNtuple/V1p17/MC_Fall18/v2/sixie/"

datasetList = OrderedDict()

#2016 hqu ntuples

#datasetList['ParkingBPH1_2018A-20Jun2021-v1.txt'] = [1, 1, "2016", ""]
#datasetList['Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018A-17Sep2018.txt'] = [1, 1, "2018", ""]
#datasetList['Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018.txt'] = [1, 1, "2018", ""] 
#datasetList['Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018C-17Sep2018.txt'] = [1, 1, "2018", ""]
#datasetList['SingleMuon_2018B.txt'] = [1, 1, "2018", ""]
datasetList['ParticleGun_K0Lpt10.txt'] = [1, 1, "2018", ""]
datasetList['ParticleGun_K0Lpt5.txt'] = [1, 1, "2018", ""]    
datasetList['ParticleGun_K0Lpt2.txt'] = [1, 1, "2018", ""]    
datasetList['ParticleGun_Photonpt10.txt'] = [1, 1, "2018", ""]
datasetList['ParticleGun_Photonpt5.txt'] = [1, 1, "2018", ""] 
datasetList['ParticleGun_Photonpt2.txt'] = [1, 1, "2018", ""] 
datasetList['ParticleGun_KPluspt10.txt'] = [1, 1, "2018", ""] 


#CMSSW_BASE_DIR = os.getenv('CMSSW_BASE')
CMSSW_BASE_DIR = "/uscms/home/sxie/work/releases/HH/CMSSW_10_6_5/"
Analyzer_DIR = CMSSW_BASE_DIR+"/src/llp_analyzer/"

#create directory for condor jobs

for listfile in datasetList.keys():

    datasetName = listfile.replace(".txt","")
    print "Preparing analyzer workflow for dataset :" + datasetName + "\n"
    if not os.path.exists(Analyzer_DIR+"/lists/" + datasetListDir + "/" + listfile):
        print "listfile: " + Analyzer_DIR+"/lists/" + datasetListDir + "/" + listfile + " does not exist. skipping."
        continue

    outputDirectory = outputDirectoryBase + datasetName + "/"
    tmpListFile = open(Analyzer_DIR + "/lists/" + datasetListDir + "/" + listfile,"r")

    year = datasetList[listfile][2]
    sampleName = datasetList[listfile][3]

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
    os.system("cp " + Analyzer_DIR + "/scripts/run_job_LPC.sh " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/")
    os.system("cp " + Analyzer_DIR + "/" + analysis + " " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/")
    os.system("mkdir -p " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/")
    #os.system("mkdir -p " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/PileupWeights/")
    #os.system("cp " + Analyzer_DIR + "/data/PileupWeights/PileupWeights.root " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/PileupWeights/")
    #os.system("mkdir -p " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/JEC/Summer16_07Aug2017_V11_MC/")
    #os.system("cp " + Analyzer_DIR + "/data/JEC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_Uncertainty_AK8PFPuppi.txt " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/JEC/Summer16_07Aug2017_V11_MC/")
    #os.system("mkdir -p " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/JEC/Fall17_17Nov2017_V32_MC/")
    #os.system("cp " + Analyzer_DIR + "/data/JEC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK8PFPuppi.txt " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/JEC/Fall17_17Nov2017_V32_MC/")
    #os.system("mkdir -p " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/JEC/Autumn18_V19_MC/")
    #os.system("cp " + Analyzer_DIR + "/data/JEC/Autumn18_V19_MC/Autumn18_V19_MC_Uncertainty_AK8PFPuppi.txt " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/HHBoostedAnalyzer/data/JEC/Autumn18_V19_MC/")


    #####################################
    #Create Condor JDL file
    #####################################
    tmpCondorJDLFile = open(Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/task.jdl","w+")
    tmpCondorJDLFileTemplate = """
Universe  = vanilla
Executable = ./run_job_LPC.sh
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    tmpCondorJDLFile.write("Arguments = " + analysis + " " + str(isData) + " " + str(option) + " " + "$(I) " + outputfile + " " + outputDirectory + " " + cmsswReleaseVersion + " " + year + " " + sampleName + "\n")

    tmpCondorJDLFileTemplate = """
Log = log/job.$(Cluster).$(Process).log
Output = out/job.$(Cluster).$(Process).out
Error = err/job.$(Cluster).$(Process).err
x509userproxy = $ENV(X509_USER_PROXY)
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    tmpCondorJDLFile.write("transfer_input_files = " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/run_job_LPC.sh, " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/input_list.tgz, " + Analyzer_DIR + "/condor/" + analysis + "_" + label + "/" + datasetName + "/" + analysis + "\n")

    tmpCondorJDLFileTemplate = """
should_transfer_files = YES
when_to_transfer_output = ON_EXIT

# Resources request
"""
    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    tmpCondorJDLFile.write("RequestMemory = 2000 \n")

    tmpCondorJDLFileTemplate = """

# Jobs selection
Queue I from (
"""

    tmpCondorJDLFile.write(tmpCondorJDLFileTemplate)
    for i in range(1,nJobs+1):
        tmpCondorJDLFile.write(str(i)+"\n")
    tmpCondorJDLFile.write(")\n")
    tmpCondorJDLFile.close()





