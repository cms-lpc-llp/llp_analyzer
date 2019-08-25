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
outputfile=$7
outputDirectory=$8
code_dir_suffix=$9

currentDir=`pwd`
homeDir=/data/christiw/
runDir=${currentDir}/christiw_${code_dir_suffix}/
rm -rf ${runDir}
mkdir -p ${runDir}

if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then
	#setup cmssw
	cd ${homeDir}LLP/CMSSW_9_4_4/src/
	workDir=`pwd`
	echo "entering directory: ${workDir}"
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	export SCRAM_ARCH=slc7_amd64_gcc630
	ulimit -c 0
	eval `scram runtime -sh`
	echo `which root`

	cd ${runDir}
	echo "entering directory: ${runDir}"
	if [ -f $CMSSW_BASE/src/llp_analyzer/RazorRun_T2 ]
	then
		cp $CMSSW_BASE/src/llp_analyzer/RazorRun_T2 ./
		mkdir -p JEC
		cp -r $CMSSW_BASE/src/llp_analyzer/data/JEC/Summer16_23Sep2016V3_MC/ JEC/Summer16_23Sep2016V3_MC
		if [ -d JEC/Summer16_23Sep2016V3_MC ]
		then
			echo "copied JEC parameters"
		else
			echo "didn't copy JEC parameters"
		fi	
		#get grid proxy
		export X509_USER_PROXY=${homeDir}x509_proxy

		#copy pedestal file
                if [ ${isData} == "yes" ]
                then
			date
			echo "copying pedestal file to local directory ========="
                        #hadoop fs -get /store/group/phys_susy/razor/Run2Analysis/EcalTiming/EcalPedestals_Legacy2016_time_v1/tree_EcalPedestals_Legacy2016_time_v1_G12rmsonly.root ./
                        cp /mnt/hadoop/store/group/phys_susy/razor/Run2Analysis/EcalTiming/EcalPedestals_Legacy2016_time_v1/tree_EcalPedestals_Legacy2016_time_v1_G12rmsonly.root ./
			echo "done with copying pedestal file to local directory ========="
			date
                fi

		#run the job
		cat ${CMSSW_BASE}${inputfilelist} | awk "NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})" > inputfilelistForThisJob_${jobnumber}.txt
		echo ""
		echo "************************************"
		echo "Running on these input files:"
		cat inputfilelistForThisJob_${jobnumber}.txt
		echo "************************************"
		echo ""
		echo " "; echo "Starting razor run job now"; echo " ";
		echo ./RazorRun_T2 inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile}
		./RazorRun_T2 inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile} -l=${label}
		echo ${outputfile}
		echo ${outputDirectory}

		##^_^##
		echo "RazorRun_T2 finished"
		date

		sleep 2
		echo "I slept for 2 second"

		##job finished, copy file to T2
		echo "copying output file to /mnt/hadoop/${outputDirectory}"
		if [ -f ${outputfile} ]
		then
			eval `scram unsetenv -sh`
			gfal-mkdir -p gsiftp://transfer.ultralight.org//${outputDirectory}
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
#rm -rf ${runDir}
echo "Job finished"
date
