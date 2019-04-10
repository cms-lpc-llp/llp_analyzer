#!/bin/tcsh

hostname

echo "Arguments: $*"
set eventList=$1
set inputfilelist=$2
set treename=$3
set runBranch=$4
set eventBranch=$5
set filePerJob=$6
set jobnumber=$7
set outputfile=$8
set outputDirectory=$9

echo " "; echo "Initialize CMSSW"; echo " "
#setenv KRB5CCNAME /home/sixie/.krb5/ticket
set workDir=`pwd`

setenv SCRAM_ARCH slc6_amd64_gcc491
cd    /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/
eval `scramv1 runtime -csh`
cd -

pwd
# env

cp $CMSSW_BASE/src/RazorAnalyzer/Tools/SkimEvents/CreateSkimmed ./

echo " "; echo "Show where we are"; echo " "
hostname
pwd
## env

klist

#setenv STAGE_SVCCLASS cmsprod

#Do Job splitting and make input file list
cat $inputfilelist | awk "NR > ($jobnumber*$filePerJob) && NR <= (($jobnumber+1)*$filePerJob)" >! inputfilelistForThisJob_${jobnumber}.txt
echo ""
echo "************************************"
echo "Running on these input files:"
cat inputfilelistForThisJob_${jobnumber}.txt
echo "************************************"
echo ""


# Get ready to run in your home directory
echo " "; echo "Starting razor run job now"; echo " ";
echo "./CreateSkimmed --event_list=${eventList} --list_of_ntuples=inputfilelistForThisJob_${jobnumber}.txt --tree_name=${treename} --run_branch=${runBranch} --event_branch=${eventBranch} --output_name=${outputfile}"
./CreateSkimmed --event_list=${eventList} --list_of_ntuples=inputfilelistForThisJob_${jobnumber}.txt --tree_name=${treename} --run_branch=${runBranch} --event_branch=${eventBranch} --output_name=${outputfile}

ls -ltr 

echo $outputfile 
echo $outputDirectory

#Do below only for output to CERN EOS
cmsMkdir $outputDirectory
cmsStage -f $outputfile $outputDirectory
#cmsStage -f ${outputfile}.log $outputDirectory

set status=`echo $?`
echo "Status: $status"

hostname

exit $status
