#!/bin/bash

#hostname

echo "Arguments: $*"
cmssw=$1
analysisType=$2
inputfilelist=$3
isData=$4
option=$5
filePerJob=$6
jobnumber=$7
outputfile=$8
outputDirectory=$9

echo " "; echo "Initialize CMSSW"; echo " "
workDir=`pwd`

SCRAM_ARCH=slc6_amd64_gcc491
cd  ${cmssw}
eval `scramv1 runtime -csh`
cd -

pwd
# env

cp $CMSSW_BASE/src/RazorAnalyzer/RazorRun ./

echo " "; echo "Show where we are"; echo " "
hostname
pwd
## env

klist

#Do Job splitting and make input file list
ls $inputfilelist
cat $inputfilelist | awk "NR > ($jobnumber*$filePerJob) && NR <= (($jobnumber+1)*$filePerJob)" > inputfilelistForThisJob_${jobnumber}.txt
echo ""
echo "************************************"
echo "Running on these input files:"
cat inputfilelistForThisJob_${jobnumber}.txt
echo "************************************"
echo ""

# Get ready to run in your home directory
echo " "; echo "Starting razor run job now"; echo " ";
if [[ "${isData}" == "false" ]]; then
    echo ./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} --outputFile=${outputfile} --option=${option}
    ./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} --outputFile=${outputfile} --option=${option} > ${outputfile}.log
else
    echo ./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} --outputFile=${outputfile} --option=${option} --isData
    ./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} --outputFile=${outputfile} --option=${option} --isData > ${outputfile}.log
fi

ls -ltr 

#Do below only for output to CERN EOS

cp $outputfile $outputDirectory

set status=`echo $?`
echo "Status: $status"

hostname

exit $status
