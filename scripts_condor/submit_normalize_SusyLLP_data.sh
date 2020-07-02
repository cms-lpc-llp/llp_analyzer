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
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016D-HighMET-07Aug17-v1_v5_v1
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016F-HighMET-07Aug17-v1_v5_v1
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-HighMET-07Aug17_ver1-v1_v5_v1
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-HighMET-07Aug17_ver2-v1_v5_v1
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016C-HighMET-07Aug17-v1_v5_v1
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016E-HighMET-07Aug17-v1_v5_v1
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016G-HighMET-07Aug17-v1_v5_v1
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016H-HighMET-07Aug17-v1_v5_v1
)


#Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016B-HighMET-07Aug17_ver1
#Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016B-HighMET-07Aug17_ver2
#Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016C-HighMET-07Aug17
#Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016D-HighMET-07Aug17
#Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016E-HighMET-07Aug17
#Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016F-HighMET-07Aug17
#Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016G-HighMET-07Aug17
#Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016H-HighMET-07Aug17

listData2016_AOD=(
Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016E-07Aug17-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016D-07Aug17-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016F-07Aug17-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016C-07Aug17-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver2-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016G-07Aug17-v1_v5_v2 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016H-07Aug17-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver1-v1_v5_v2 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver1-v1_v5_v1 \
)

listData2017=(
Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017D-HighMET-17Nov2017-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017B-HighMET-17Nov2017-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017F-HighMET-17Nov2017-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017E-HighMET-17Nov2017-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2017_Run2017C-HighMET-17Nov2017-v1_v5_v1 \
)

#Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017B-HighMET-17Nov2017
#Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017C-HighMET-17Nov2017
#Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017D-HighMET-17Nov2017
#Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017E-HighMET-17Nov2017
#Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017F-HighMET-17Nov2017

listData2017_AOD=(
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017F-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017D-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017C-17Nov2017-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017B-17Nov2017-v1_v5_v1 \

Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017C-17Nov2017-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017B-17Nov2017-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017D-17Nov2017-v1_v5_v3 \
Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017B-17Nov2017-v1_v5_v3 \
Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017-v1_v5_v2 \
Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017-v1_v5_v3 \
Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017-v1_v5_v4 \
)

listData2018=(
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_Run2018A-HighMET-17Sep2018-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_Run2018D-HighMET-PromptReco-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_Run2018E-HighMET-PromptReco-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_Run2018D-HighMET-PromptReco-v2_v5_v1 \
)

#Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018A-HighMET-17Sep2018
#Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018
#Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018
#Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018D-HighMET-PromptReco

listData2018_AOD=(
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018C-17Sep2018-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018A-17Sep2018-v1_v5_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018-v1_v5_v2 \
Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018-v1_v5_v4 \

#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018C-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018A-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018-v1_v5_v1 \
#Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_AOD_Run2018D-PromptReco-v2_v5_v2 \
)
for year in \
Data2018_AOD \

#Data2016_AOD \
#Data2017_AOD \

#Data2018
#Data2016 \
#Data2017 \
do
        echo ${year}
	#inputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/${ver}/v6/${year}/SinglePhoton/
	#inputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/${ver}/v7/${year}/MuonEG/
	#inputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/${ver}/v8/${year}/Zmumu/
	#inputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/${ver}/v8/${year}/Zee/
	inputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/${ver}/v9/${year}/JetHT/
	#inputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/${ver}/v8/${year}/SingleMuon/
	#inputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/${ver}/v8/${year}/SingleElectron/
	#inputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/${ver}/v5/${year}/
	echo ${inputDir}
	outputDir=${inputDir}normalized
        sampleList=list${year}[@]
       
	for sample in "${!sampleList}"
        do


		echo "Sample " ${sample}
		analyzer=SusyLLP
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
