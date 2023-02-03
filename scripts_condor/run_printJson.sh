#!/bin/sh

hostname
echo "Job started"
date
start_time=`date +%s`
inputfilelist=$1
filePerJob=$2
jobnumber=$3
maxjob=$4
sample=${inputfilelist##*/}
sample=${sample%.txt}
outputfile=${sample}_Job${jobnumber}_of_${maxjob}.json
outputDirectory=$5
CMSSW_BASE=${6}
homeDir=${7}
currentDir=`pwd`
#user=${homeDir#*/data/}
user=${homeDir#*/storage/user/}
runDir=${currentDir}/${user}_${analyzerTag}/

rm -rf ${runDir}
mkdir -p ${runDir}

if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then
	#setup cmssw
	#cd $CMSSW_BASE/src/
	cd ${CMSSW_BASE}/src/
	echo "$6"
	echo "${CMSSW_BASE}"
	workDir=`pwd`
	echo "entering directory: ${workDir}"
	ulimit -c 0
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	export SCRAM_ARCH=slc7_amd64_gcc700
	eval `scramv1 runtime -sh`
	echo `which root`

	cd ${runDir}
	echo "entering directory: ${runDir}"
	echo "${CMSSW_BASE}/src/llp_analyzer/python/printJson.py"
	if [ -f ${CMSSW_BASE}/src/llp_analyzer/python/printJson.py ]
	then
		cp $CMSSW_BASE/src/llp_analyzer/python/printJson.py ./
		#get grid proxy
		#export X509_USER_PROXY=${homeDir}x509_proxy
		echo "${homeDir}x509_proxy"
		voms-proxy-info


		#run the job
		cat ${CMSSW_BASE}${inputfilelist} | awk "NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})" > inputfilelistForThisJob_${jobnumber}.txt
		echo ""
		echo "************************************"
		echo "Running on these input files:"
		cat inputfilelistForThisJob_${jobnumber}.txt
		echo "************************************"
		echo ""
		echo " "; echo "Starting razor run job now"; echo " ";
			python printJson.py -r runNum -l lumiNum -o ${outputfile} -i inputfilelistForThisJob_${jobnumber}.txt
		echo ${outputfile}
		echo ${outputDirectory}
		echo "printJson.py finished"
		cat ${outputfile}
		date

		sleep 2
		echo "I slept for 2 second"

		##job finished, copy file to T2
		echo "copying output file to ${outputDirectory}"
		eval `scram unsetenv -sh`
		mkdir -p ${outputDirectory}
		cp ${outputfile} ${outputDirectory}/${outputfile}
		if [ -f ${outputDirectory}/${outputfile} ]
		then
			echo ${outputfile} "copied"
		else
			echo ${outputfile} "not copied"
		fi

	else
		echo echo "WWWWYYYY ============= failed to access file RazorRun_T2, job anandoned"
	fi

else
	echo "VVVVYYYY ============= failed to access file /cvmfs/cms.cern.ch/cmsset_default.sh, job abandoned"
fi

cd ${currentDir}
#rm -rf ${runDir}
echo "Job finished"
date
end_time=`date +%s`
runtime=$((end_time-start_time))
echo ${runtime}
