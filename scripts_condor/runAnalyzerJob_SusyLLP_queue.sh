#!/bin/sh

hostname
echo "Job started"
date

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
process=${12}


currentDir=`pwd`
user=${homeDir#*/storage/user/}
runDir=${currentDir}/${user}_${analyzerTag}_${sample}_${option}_Job${jobnumber}/
#runDir=${currentDir}/${user}_${analyzerTag}_${sample}_Job${jobnumber}/
rm -rf ${runDir}
mkdir -p ${runDir}

if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then

	#setup cmssw
	cd ${CMSSW_BASE}/src/
	#cd ${homeDir}/login-1/jmao/CMSSW_9_4_4/src/
	workDir=`pwd`
	echo "entering directory: ${workDir}"
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	export SCRAM_ARCH=slc7_amd64_gcc630
	ulimit -c 0
	eval `scram runtime -sh`
	echo `which root`

	cd ${runDir}
	echo "entering directory: ${runDir}"
	echo "CMSSW base: ${CMSSW_BASE}"
	if [ -f $CMSSW_BASE/src/cms_lpc_llp/llp_analyzer/RazorRun_T2_JM ]
	then
		cp $CMSSW_BASE/src/cms_lpc_llp/llp_analyzer/RazorRun_T2_JM ./


		ls
		pwd

		mkdir -p JEC
		cp -r $CMSSW_BASE/src/cms_lpc_llp/llp_analyzer/data/JEC/Summer16_23Sep2016V3_MC/ JEC/Summer16_23Sep2016V3_MC
		if [ -d JEC/Summer16_23Sep2016V3_MC ]
		then
			echo "copied JEC parameters"
		else
			echo "didn't copy JEC parameters"
		fi	
		#get grid proxy
		export X509_USER_PROXY=${homeDir}/x509_proxy
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
		echo ./RazorRun_T2_JM inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile} -l=${analyzerTag} -p=${process}
		./RazorRun_T2_JM inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile} -l=${analyzerTag} -p=${process}
		#./RazorRun_T2_JM inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile} -l=${analyzerTag} -p=${process}
		echo ${outputfile}
		echo ${outputDirectory}

		##^_^##
		echo "RazorRun_T2_JM finished"
		date

		sleep 2
		echo "I slept for 2 second"

		##job finished, copy file to T2
		echo "copying output file to /mnt/hadoop/${outputDirectory}"
		if [ -f ${outputfile} ]
		then
			eval `scram unsetenv -sh`
			voms-proxy-info
			echo "gfal-mkdir -p gsiftp://transfer.ultralight.org//${outputDirectory}"
			gfal-mkdir -p gsiftp://transfer.ultralight.org//${outputDirectory}
			echo "mkdir done, start copying"
			ls
			pwd
			gfal-copy -f ${outputfile} gsiftp://transfer.ultralight.org//${outputDirectory}/${outputfile}
		else
			echo "output doesn't exist"
		fi
		if [ -f /mnt/hadoop/${outputDirectory}/${outputfile} ]
		then
			echo "ZZZZAAAA ============ good news, job finished successfully "
		else
			echo "YYYYZZZZ ============ somehow job failed, please consider resubmitting"
		fi
	else
		echo echo "WWWWYYYY ============= failed to access file RazorRun_T2, job anandoned"
	fi

else
	echo "VVVVYYYY ============= failed to access file /cvmfs/cms.cern.ch/cmsset_default.sh, job abandoned"
fi

cd ${currentDir}
rm -rf ${runDir}
echo "Job finished"
date