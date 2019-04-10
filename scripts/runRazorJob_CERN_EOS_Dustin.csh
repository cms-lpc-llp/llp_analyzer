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
set workDir=`pwd`

setenv SCRAM_ARCH slc6_amd64_gcc481
cd    $cmsswDir
eval `scramv1 runtime -csh`
cd -

pwd

eosdir="/eos/cms/store/group/phys_susy/razor/Run2Analysis/Analyzers/"
eos cp ${eosdir}/RazorRunAuxFiles_Expanded.tar.gz ./
eos cp ${eosdir}/RazorRun_NoAFS ./
eos cp ${eosdir}/Run${analysisType} ./

echo " "; echo "Show where we are"; echo " "
hostname
pwd

klist

#Do Job splitting and make input file list
cat $inputfilelist | awk "NR > ($jobnumber*$filePerJob) && NR <= (($jobnumber+1)*$filePerJob)" >! inputfilelistForThisJob_${jobnumber}.txt
echo ""
echo "************************************"
echo "Running on these input files:"
cat inputfilelistForThisJob_${jobnumber}.txt
echo "************************************"
echo ""

set datastring=""
if (${isData} == 1) then
    set datastring="--isData "
endif

# Get ready to run in your home directory
echo " "; echo "Starting razor run job now"; echo " ";
echo ./RazorRun_NoAFS inputfilelistForThisJob_${jobnumber}.txt ${analysisType} ${datastring}-f=${outputfile} -n=${option}
./RazorRun_NoAFS inputfilelistForThisJob_${jobnumber}.txt ${analysisType} ${datastring}-f=${outputfile} -n=${option} -l=${label} |& tee ${outputfile}.log

ls -ltr 

echo $outputfile 
echo $outputDirectory

#Do below only for output to CERN EOS
eos mkdir -p $outputDirectory
eos cp $outputfile /eos/cms/$outputDirectory/

set tempOutputfile = `echo $outputfile | sed 's/.root//'`
foreach f ( ${tempOutputfile}_*.root )
   eos cp $f /eos/cms/$outputDirectory/
end

set status=`echo $?`
echo "Status: $status"

hostname

exit $status
