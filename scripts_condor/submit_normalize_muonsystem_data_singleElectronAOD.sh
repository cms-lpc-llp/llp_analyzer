#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

mode=bkg
job_script=${RazorAnalyzerDir}/scripts_condor/normalize.sh
ver=V1p17

listData2016=(
Run2_displacedJetMuonNtupler_${ver}_Data2016_AOD_Run2016B-07Aug17_ver1
Run2_displacedJetMuonNtupler_${ver}_Data2016_AOD_Run2016B-07Aug17_ver2
Run2_displacedJetMuonNtupler_${ver}_Data2016_AOD_Run2016C-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_AOD_Run2016D-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_AOD_Run2016E-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_AOD_Run2016F-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_AOD_Run2016G-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_AOD_Run2016H-07Aug17
)
listData2017=(
Run2_displacedJetMuonNtupler_${ver}_Data2017_AOD_Run2017B-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_AOD_Run2017C-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_AOD_Run2017D-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_AOD_Run2017E-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_AOD_Run2017F-17Nov2017
)
listData2018ABC=(
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_AOD_Run2018A-17Sep2018-v2_part2
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_AOD_Run2018A-17Sep2018-v2
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_AOD_Run2018B-17Sep2018
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_AOD_Run2018C-17Sep2018
)
listData2018D=(
Run2_displacedJetMuonNtupler_${ver}_Data2018D_17Sept2018_AOD_Run2018D-PromptReco
)
for year in \
Data2016 \
Data2017 \
Data2018ABC \
Data2018D
do
        echo ${year}
	inputDir=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/${year}_AOD/v5/v64/
	echo ${inputDir}
	outputDir=${inputDir}normalized
        sampleList=list${year}[@]
       
	for sample in "${!sampleList}"
        do


		echo "Sample " ${sample}
		analyzer=llp_MuonSystem
		rm -f submit/${analyzer}_normalize_${mode}_${sample}*.jdl
		rm -f log/${analyzer}_normalize_${mode}_${sample}*
		jdl_file=submit/${analyzer}_normalize_${mode}_${year}_${sample}.jdl
		echo "Universe = vanilla" > ${jdl_file}
		echo "Executable = ${job_script}" >> ${jdl_file}
		echo "Arguments =${mode} yes ${sample} ${inputDir} ${outputDir} ${CMSSW_BASE} ${HOME}/"  >> ${jdl_file}

		echo "Log = log/${analyzer}_normalize_${mode}_${sample}_PC.log" >> ${jdl_file}
		echo "Output = log/${analyzer}_normalize_${mode}_${sample}_\$(Cluster).\$(Process).out" >> ${jdl_file}
		echo "Error = log/${analyzer}_normalize_${mode}_${sample}_\$(Cluster).\$(Process).err" >> ${jdl_file}

		#echo "Requirements=TARGET.OpSysAndVer==\"CentOS7\"" >> ${jdl_file}
		echo "Requirements=(TARGET.OpSysAndVer==\"CentOS7\" && regexp(\"blade.*\", TARGET.Machine))" >> ${jdl_file}
		echo "RequestMemory = 4000" >> ${jdl_file}
		echo "RequestCpus = 1" >> ${jdl_file}
		echo "RequestDisk = 4" >> ${jdl_file}

		echo "+RunAsOwner = True" >> ${jdl_file}
		echo "+InteractiveUser = true" >> ${jdl_file}
		#echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7\"" >> ${jdl_file}
		echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7-m202006\"" >> ${jdl_file}
		echo '+SingularityBindCVMFS = True' >> ${jdl_file}
		echo "run_as_owner = True" >> ${jdl_file}
		echo "x509userproxy = ${HOME}/x509_proxy" >> ${jdl_file}
		echo "should_transfer_files = YES" >> ${jdl_file}
		echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
		echo "Queue 1" >> ${jdl_file}
		echo "condor_submit ${jdl_file}"
		condor_submit -batch-name ${sample} ${jdl_file}
	done
done
