#!/bin/tcsh

hostname

echo "Arguments: $*"
set analysisType=$1
set inputfilelist=$2
set isData=$3
set option=$4
set filePerJob=$5
set jobnumber=$6
set outputfile=$7
set outputDirectory=$8

echo " "; echo "Initialize CMSSW"; echo " "
#setenv KRB5CCNAME /home/sixie/.krb5/ticket
set workDir=`pwd`

setenv SCRAM_ARCH slc6_amd64_gcc491
cd    /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/
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
echo ./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} ${isData} --outputFile=${outputfile} --optionNumber=${option}
./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} ${isData} --outputFile=${outputfile} --optionNumber=${option} |& tee ${outputfile}.log

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
