#!/bin/sh

hostname
echo "Job started"
date

mode=$1
sample=$2
root_dir=$3/${sample}/
root_file=/store/group/phys_exotica/delayedjets/llp_analyzer/V1p0/MC_Summer16/v2/bkg/${sample}/${sample}.root

currentDir=`pwd`
homeDir=/data/jmao/
runDir=${currentDir}/jmao_${mode}_${sample}/
normalize_file=llp_${mode}_${sample}.txt
rm -rf ${runDir}
mkdir -p ${runDir}

if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then

	#setup cmssw
	cd ${homeDir}/CMSSW_9_4_4/src/
	workDir=`pwd`
	echo "entering directory: ${workDir}"
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	export SCRAM_ARCH=slc7_amd64_gcc630
	ulimit -c 0
	eval `scram runtime -sh`
	echo `which root`

	cd ${runDir}
        echo "entering directory: ${runDir}"

	echo "/mnt/hadoop/${root_dir}/${sample}*_Job*.root"
	hadd ${sample}.root /mnt/hadoop/${root_dir}/${sample}*_Job*.root

	#copy file over
#	eval `scram unsetenv -sh`
#	echo "$root_file"
#	echo "gfal-copy gsiftp://transfer.ultralight.org//$root_file ."
#	if [ -f /mnt/hadoop/$root_file ]
#	then
#        	gfal-copy gsiftp://transfer.ultralight.org//$root_file .
#	else
#        	echo "input ROOT file doesn't exist"
#	fi

	eval `scramv1 runtime -sh`
	if [ -f $CMSSW_BASE/src/cms_lpc_llp/llp_analyzer/data/xSections.dat ]
	then
		mkdir -p data
		cp $CMSSW_BASE/src/cms_lpc_llp/llp_analyzer/data/xSections.dat data/xSections.dat 
	else
		echo "data/xSections.dat doesn't exist"

	fi	
	
	#create normalization file
	rm -f $normalize_file
	dataset=${sample}
	echo "$dataset"
	echo "${dataset} ${runDir}/${sample}.root" > $normalize_file
	if [ -f $normalize_file ]
	then
		echo "normalization file created"
	fi

	#normalize
	if [ -f $CMSSW_BASE/src/cms_lpc_llp/llp_analyzer/NormalizeNtuple ]
        then
                cp $CMSSW_BASE/src/cms_lpc_llp/llp_analyzer/NormalizeNtuple ./
	        ./NormalizeNtuple ${normalize_file} 1
	else
		echo "NormalizeNtuple not found"
	fi
	echo "Normalization done"
	
	sleep 2
        echo "I slept for 2 second"
	
	# copy normalized file back to hadoop
	eval `scram unsetenv -sh`
	gfal-copy -f ${runDir}/${sample}_1pb_weighted.root gsiftp://transfer.ultralight.org//${root_dir}/${sample}_1pb_weighted.root
	

	if [ -f /mnt/hadoop/${root_dir}/${sample}_1pb_weighted.root ]
	then
	        echo "copied succeed"
	else
	        echo "copied failed"
	fi


else
	echo "VVVVYYYY ============= failed to access file /cvmfs/cms.cern.ch/cmsset_default.sh, job abandoned"
fi

cd ${currentDir}
#rm -rf ${runDir}
echo "Job finished"
date
