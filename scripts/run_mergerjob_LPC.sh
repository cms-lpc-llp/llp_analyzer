#!/bin/bash

#######################
#debugging purposes
#######################
voms-proxy-info --all
ls -l

############################
#define input parameters
############################
analysisType=$1
isData=$2
option=$3
jobnumber=$4
outputfile=$5
outputDirectory=$6
cmsswReleaseVersion=$7
optionLabel=$8
tauName=$9

############################
#define exec and setup cmssw
############################
workDir=`pwd`
executable=${analysisType}
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc820
#tar -zxvf cms_setup.tar.gz
scramv1 project CMSSW $cmsswReleaseVersion

#########################################
#copy input list and exec to cmssw folder
########################################
cp input_list.tgz $cmsswReleaseVersion/src/
cp ${executable} $cmsswReleaseVersion/src/.
cp JEC.tar.gz *.root $cmsswReleaseVersion/src/
###########################
#get cmssw environment
###########################
cd $cmsswReleaseVersion/src/
eval `scram runtime -sh`
tar vxzf input_list.tgz
tar vxzf JEC.tar.gz
inputfilelist=input_list_${jobnumber}.txt

###################################
#copy input files ahead of time
###################################
mkdir inputs/
for i in `cat $inputfilelist`
do
echo "Copying Cluster Input File: " $i
xrdcp $i ./inputs/
done
hadd -f tmp_input_list.root inputs/*.root 
rm inputs -rv
mkdir tauinputs/
echo "Copying Tau Input File: " ${tauName}
xrdcp ${tauName} ./tauinputs/
hadd -f tmp_tau.root tauinputs/*.root
rm tauinputs -rv 
###########################
#run executable
###########################
echo "Executing Analysis executable:"
echo "./${executable}  tmp_input_list.root tmp_tau.root  ${outputfile}_${jobnumber}.root ${isData}"
./${executable}  tmp_input_list.root tmp_tau.root  ${outputfile}_${jobnumber}.root ${isData}
ls -l
##########################################################
#copy outputfile to /eos space -- define in submitter code
##########################################################
eosmkdir -p ${outputDirectory}
xrdcp -f ${outputfile}_${jobnumber}.root root://cmseos.fnal.gov/${outputDirectory}/${outputfile}_${jobnumber}.root 
rm ${outputfile}_${jobnumber}.root
rm tmp_input_list.root
rm tmp_tau.root
cd -
