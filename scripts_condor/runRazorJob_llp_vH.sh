#!/bin/sh

hostname
echo "Job started"
date
start_time=`date +%s`
analysisType=$1
inputfilelist=$2
isData=$3
option=$4
filePerJob=$5
jobnumber=$6
maxjob=$7
sample=${inputfilelist##*/}
sample=${sample%.txt}
outputfile=${sample}_Job${jobnumber}_of_${maxjob}.root
outputDirectory=$8
analyzerTag=$9
CMSSW_BASE=${10}
homeDir=${11}
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
	workDir=`pwd`
	echo "entering directory: ${workDir}"
	ulimit -c 0
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	export SCRAM_ARCH=slc7_amd64_gcc630
	eval `scram runtime -sh`
	echo `which root`

	cd ${runDir}
	echo "entering directory: ${runDir}"
	echo "${CMSSW_BASE}/src/llp_analyzer/RazorRun_T2"
	if [ -f ${CMSSW_BASE}/src/llp_analyzer/RazorRun_T2 ]
	then
		cp $CMSSW_BASE/src/llp_analyzer/RazorRun_T2 ./
		#get grid proxy
		export X509_USER_PROXY=${homeDir}x509_proxy
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
		if [ ${analysisType} == "MakeMCPileupDistribution" ]
		then
			echo "./RazorRun_T2 inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -f=${outputfile}"
			./RazorRun_T2 inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -f=${outputfile}
		else
			echo ./RazorRun_T2 inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile} -l=${analyzerTag}
			./RazorRun_T2 inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile} -l=${analyzerTag}
		fi
		
		echo ${outputfile}
		echo ${outputDirectory}
		ls *root > output.txt
		echo "Output ROOT files: "
		cat output.txt
		##^_^##
		echo "RazorRun_T2 finished"
		date

		sleep 2
		echo "I slept for 2 second"

		##job finished, copy file to T2
		echo "copying output file to ${outputDirectory}"
		eval `scram unsetenv -sh`
		#gfal-mkdir -p davs://xrootd-redir-stageout.ultralight.org:1095/${outputDirectory}
		#gfal-mkdir -p gsiftp://transfer-lb.ultralight.org//${outputDirectory}	
		mkdir -p ${outputDirectory}	
		while IFS= read -r line
		do
        		echo $line	
			#gfal-copy -f --checksum-mode=both ${line} davs://xrootd-redir-stageout.ultralight.org:1095/${outputDirectory}/${line}
			cp  ${line} ${outputDirectory}/${line}
			if [ -f ${outputDirectory}/${line} ]
			then
				echo ${line} "copied"
			else
				echo ${line} "not copied"
			fi
		done <"output.txt"

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
