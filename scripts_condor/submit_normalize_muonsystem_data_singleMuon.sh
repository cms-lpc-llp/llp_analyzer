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
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016B-ZMu-07Aug17_ver1
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016B-ZMu-07Aug17_ver2
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016C-ZMu-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016D-ZMu-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016E-ZMu-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016F-ZMu-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016G-ZMu-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016H-ZMu-07Aug17
)
listData2017=(
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017B-ZMu-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017C-ZMu-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017D-ZMu-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017E-ZMu-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017F-ZMu-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017G-ZMu-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017H-ZMu-17Nov2017
)
listData2018=(
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018A-ZMu-17Sep2018
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018B-ZMu-17Sep2018
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018C-ZMu-17Sep2018
Run2_displacedJetMuonNtupler_${ver}_Data2018D_17Sept2018_Run2018D-ZMu-PromptReco
)

for year in \
Data2016 \
Data2017 \
Data2018
do
        echo ${year}
	inputDir=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${ver}/${year}/v5/v83/
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
