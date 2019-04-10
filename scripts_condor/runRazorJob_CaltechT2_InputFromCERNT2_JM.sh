#!/bin/sh

hostname
echo "Job started"
date

analysisType=$1
echo "analysisType "${analysisType}
inputfilelist=$2
echo "inputfilelist "${inputfilelist}
isData=$3
echo "isData "${isData}
option=$4
echo "option "${option}
filePerJob=$5
echo "filePerJob "${filePerJob}
jobnumber=$6
echo "jobnumber "${jobnumber}
outputfile=$7
echo "outputfile "${outputfile}
outputDirectory=$8
echo "outputDirectory "${outputDirectory}
#code_dir_suffix=$9
#echo "code_dir_suffix "${code_dir_suffix}
analysisTag=$9
echo "analysisTag "${analysisTag}

currentDir=`pwd`
homeDir=/data/jmao/CMSSW_9_2_1/src/RazorAnalyzer/
#runDir=${currentDir}/jmao_${code_dir_suffix}/
#runDir=${currentDir}/jmao_${analysisType}_Job${jobnumber}/
runDir=${currentDir}/jmao_$RANDOM/
rm -rf ${runDir}
mkdir -p ${runDir}

if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then

	#setup cmssw
	cd ${homeDir}
	workDir=`pwd`
	echo "entering directory: ${workDir}"
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	export SCRAM_ARCH=slc7_amd64_gcc630
	ulimit -c 0
	eval `scram runtime -sh`
	echo `which root`

	cd ${runDir}
	echo "entering directory: ${runDir}"
	#cd ${currentDir}
	#echo "entering directory: ${currentDir}"

	if [ -f $CMSSW_BASE/src/RazorAnalyzer/RazorRun_T2_JM ]
	then 
		cp $CMSSW_BASE/src/RazorAnalyzer/RazorRun_T2_JM ./

		#get grid proxy
		export X509_USER_PROXY=${homeDir}x509_proxy

		#run the job
		cat ${CMSSW_BASE}${inputfilelist} | awk "NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})" > inputfilelistForThisJob_${jobnumber}.txt
		echo ""
		echo "************************************"
		echo "Running on these input files:"
		cat inputfilelistForThisJob_${jobnumber}.txt
		echo "************************************"
		echo ""

		nRemote=`cat inputfilelistForThisJob_${jobnumber}.txt | wc -l`
		#copy file to local directory, no need if the input is from T2
		echo "copy file to local directory ======"

		for ifile in `cat inputfilelistForThisJob_${jobnumber}.txt`
		do
			xrdcp ${ifile} ./
			echo "copied to" 
			echo ${pwd}
		done

		ls razorNtuple_*.root > inputfilelistForThisJob_${jobnumber}_local.txt
		nLocal=`cat inputfilelistForThisJob_${jobnumber}_local.txt |wc -l`

		if [ "${nRemote}" -eq "${nLocal}" ]
		then
			echo " "; echo "Starting razor run job now"; echo " ";
			echo ./RazorRun_T2_JM inputfilelistForThisJob_${jobnumber}_local.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile} 

			./RazorRun_T2_JM inputfilelistForThisJob_${jobnumber}_local.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile} -l=${analysisTag}

			echo ${analysisTag}
			echo ${outputfile}
			echo ${outputDirectory}
			mkdir -p /mnt/hadoop/${outputDirectory}

			##^_^##
			echo "RazorRun_T2_JM finished"
			date

			sleep 2
			echo "I slept for 2 second" 

			##job finished, copy file to T2
			echo "copying output file to /mnt/hadoop/${outputDirectory}"
			cp ${outputfile} /mnt/hadoop/${outputDirectory}
			if [ -f /mnt/hadoop/${outputDirectory}/${outputfile} ]
			then
				echo "ZZZZAAAA ============ good news, job finished successfully "
			else
				echo "YYYYZZZZ ============ somehow job failed, please consider resubmitting"
			fi

		else
			echo "XXXXYYYY ============= copy file failed (${nRemote} -> ${nLocal}), job abandoned"
		fi
	else
		echo echo "WWWWYYYY ============= failed to access file RazorRun_T2_JM, job anandoned"
	fi

else
	echo "VVVVYYYY ============= failed to access file /cvmfs/cms.cern.ch/cmsset_default.sh, job anandoned"
fi

cd ${currentDir}
rm -rf ${runDir}
echo "Job finished"
date
