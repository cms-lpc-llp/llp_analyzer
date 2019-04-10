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
set cmsswDir=$9
set label=$10

echo " "; echo "Initialize CMSSW"; echo " "
#setenv KRB5CCNAME /home/sixie/.krb5/ticket
set workDir=`pwd`

setenv SCRAM_ARCH slc6_amd64_gcc491
cd    $cmsswDir
eval `scramv1 runtime -csh`
cd -

pwd
# env

cp $CMSSW_BASE/src/RazorAnalyzer/RazorRun ./
eos cp /eos/cms/store/user/jmao/razorRun2Analysis/RazorRunAuxFiles_Expanded.tar.gz ./
#wget http://cmsdoc.cern.ch/~duanders/RazorRunAuxFiles_Expanded.tar.gz
tar vxzf RazorRunAuxFiles_Expanded.tar.gz
cp RazorRunAuxFiles_Expanded/* .
tar vxzf JEC_Summer16_23Sep2016V3.tgz
tar vxzf Spring16_FastSimV1.tgz

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
set datastring=""
if (${isData} == 1) then
    set datastring="--isData "
endif

echo " "; echo "Starting razor run job now"; echo " ";
echo $datastring
echo ./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} ${datastring} -f=${outputfile} -n=${option}
./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} ${datastring} -f=${outputfile} -n=${option} -l=${label} |& tee ${outputfile}.log
 
ls -ltr 

echo $outputfile 
echo $outputDirectory

#Do below only for output to CERN EOS
mkdir -p /eos/cms/$outputDirectory
cp -v $outputfile /eos/cms/$outputDirectory

set tempOutputfile = `echo $outputfile | sed 's/.root//'`
foreach f ( ${tempOutputfile}_*.root )
   cp -v $f /eos/cms/$outputDirectory
end

set status=`echo $?`
echo "Status: $status"

hostname

exit $status
