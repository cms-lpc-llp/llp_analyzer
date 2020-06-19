#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_llp_vH.sh
filesPerJob=200
#Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018 \
#Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018 \
#Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018D-HighMET-PromptReco
ver=V1p15

listData2016=(
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016B-HighMET-07Aug17_ver1
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016B-HighMET-07Aug17_ver2
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016C-HighMET-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016D-HighMET-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016E-HighMET-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016F-HighMET-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016G-HighMET-07Aug17
Run2_displacedJetMuonNtupler_${ver}_Data2016_Run2016H-HighMET-07Aug17
)
listData2017=(
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017B-HighMET-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017C-HighMET-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017D-HighMET-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017E-HighMET-17Nov2017
Run2_displacedJetMuonNtupler_${ver}_Data2017_Run2017F-HighMET-17Nov2017
)
listData2018=(
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018A-HighMET-17Sep2018
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018
Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018D-HighMET-PromptReco
)

for year in \
Data2016 \
Data2017 \
Data2018
do
	echo ${year}
	sampleList=list${year}[@]
        for sample in "${!sampleList}"
        do
		echo "Sample " ${sample}
		#output=/storage/user/christiw/displacedJetMuonAnalyzer/V1p7/MC_Summer16/v3/bkg/wH/${sample}
		version=/${ver}/${year}/v1/
		output=/store/group/phys_exotica/delayedjets/pickEvents/csc/${version}/v2/${sample}
		echo ${output}
		inputfilelist=/src/llp_analyzer/lists/displacedJetMuonNtuple/${ver}/${year}/${sample}.txt
		nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
	        maxjob=`python -c "print int($nfiles.0/$filesPerJob)+1"`
		mod=`python -c "print int($nfiles.0%$filesPerJob)"`
	        if [ ${mod} -eq 0 ]
	        then
	                maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
	        fi
		#event_list=lists/EventPick/event_pick_Run2_displacedJetMuonNtupler_V1p15_Data2016_Data2017_Data2018-HighMET_goodLumi_oot_highBDT_timeSpread20.txt
		#event_list=list/EventPick/event_pick_Run2_displacedJetMuonNtupler_V1p15_Data2016_Data2017_Data2018-HighMET_goodLumi_intime0p20p7_binD.txt
		event_list=lists/EventPick/event_pick_Run2_displacedJetMuonNtupler_V1p15_Data2016_Data2017_Data2018-HighMET_goodLumi_intime_bdtBin2_Nrechit80_dphi0p75.txt
		analyzer=EventPick
		rm -f submit/${analyzer}_${sample}*
		rm -f log/${analyzer}_${sample}*
	
		echo "job " ${maxjob}
	        jdl_file=submit/${analyzer}_${sample}_${maxjob}.jdl
	        echo "Universe = vanilla" > ${jdl_file}
	        echo "Executable = ${job_script}" >> ${jdl_file}
	        echo "Arguments = ${analyzer} ${inputfilelist} yes 01 ${filesPerJob} \$(ProcId) ${maxjob} ${output} ${event_list} ${CMSSW_BASE} ${HOME}/" >> ${jdl_file}
	
	        # option should always be 1, when running condor
	        echo "Log = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_PC.log" >> ${jdl_file}
	        echo "Output = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
	        echo "Error = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}
	
	        #echo "Requirements=TARGET.OpSysAndVer==\"CentOS7\"" >> ${jdl_file}
	        echo "Requirements=(TARGET.OpSysAndVer==\"CentOS7\" && regexp(\"blade.*\", TARGET.Machine))" >> ${jdl_file}
	
	        echo "RequestMemory = 2000" >> ${jdl_file}
	        echo "RequestCpus = 1" >> ${jdl_file}
	        echo "RequestDisk = 4" >> ${jdl_file}
	
	        echo "+RunAsOwner = True" >> ${jdl_file}
	        echo "+InteractiveUser = true" >> ${jdl_file}
	        echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7\"" >> ${jdl_file}
	        echo '+SingularityBindCVMFS = True' >> ${jdl_file}
	        echo "run_as_owner = True" >> ${jdl_file}
	        echo "x509userproxy = ${HOME}/x509_proxy" >> ${jdl_file}
	        echo "should_transfer_files = YES" >> ${jdl_file}
	        echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
	        echo "Queue ${maxjob}" >> ${jdl_file}
	        echo "condor_submit ${jdl_file}"
	        condor_submit ${jdl_file} -batch-name ${sample}


	done
done
