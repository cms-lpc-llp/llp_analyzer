#!/bin/sh

hostname
echo "Job started"
ls
pwd
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
base=`pwd`
user=${homeDir#*/storage/user/}
runDir=${currentDir}/${user}_${analyzerTag}/

rm -rf ${runDir}
mkdir -p ${runDir}

if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        export SCRAM_ARCH=sslc7_amd64_gcc630
        tar -xzf CMSSW_9_4_4.tar.gz
        mkdir -p CMSSW_9_4_4/src
        rm CMSSW_9_4_4.tar.gz
        cd CMSSW_9_4_4/src
	workDir=`pwd`
	echo "entering directory: ${workDir}"
	ulimit -c 0
	eval `scramv1 runtime -sh`
	echo `which root`

	#cd ${runDir}
	#echo "entering directory: ${runDir}"
	echo "Preparing to run printJson.py"
	if [ -f ${base}/printJson.py ]
	then
		cp ${base}/printJson.py ./
		cp ${base}/${sample}.txt ./thisList.txt

		#run the job
                echo ${base}
		cat thisList.txt | awk "NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})" > inputfilelistForThisJob_${jobnumber}.txt
		ls
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
		echo "copying output file $outputfile to ${outputDirectory}/"
		eval `scram unsetenv -sh`
		#mkdir -p ${outputDirectory}
                xrdcp -f $outputfile root://cmseos.fnal.gov//${outputDirectory}/
		#cp ${outputfile} ${outputDirectory}/${outputfile}
		#gfal-copy -f ${outputfile} ${outputDirectory}/${outputfile}
                #if [ -f ${outputDirectory}/${outputfile} ]
		#then
		#	echo ${outputfile} "copied"
		#else
		#	echo ${outputfile} "not copied"
		#fi

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
