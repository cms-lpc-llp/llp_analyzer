#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_llp_vH.sh
filesPerJob=1
#Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018 \
#Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018 \
#Run2_displacedJetMuonNtupler_${ver}_Data2018_17Sept2018_Run2018D-HighMET-PromptReco

#ParkingBPH4_Run2018A.txt
for sample in \
ParkingBPH4_2018A
do
	echo "Sample " ${sample}
	ver=V1p19
	year=Data2018_UL
	#version=/${ver}/${year}/MuonHitsOnly/
	version=/${ver}/${year}/
	output=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/bparking/${version}/v9/${sample}

	echo ${output}
	inputfilelist=/src/llp_analyzer/lists/displacedJetMuonNtuple/${version}/${sample}.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
        maxjob=`python -c "print int($nfiles.0/$filesPerJob)+1"`
	mod=`python -c "print int($nfiles.0%$filesPerJob)"`
        if [ ${mod} -eq 0 ]
        then
                maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
        fi
	analyzer=llp_MuonSystem_bparking
        echo "Data18 condition"
        analyzerTag=Razor2018_17SeptEarlyReReco
	rm -f submit/${analyzer}_${sample}*
	rm -f log/${analyzer}_${sample}*

	echo "job " ${maxjob}
        jdl_file=submit/${analyzer}_${sample}_${maxjob}.jdl
        echo "Universe = vanilla" > ${jdl_file}
        echo "Executable = ${job_script}" >> ${jdl_file}
        echo "Arguments = ${analyzer} ${inputfilelist} yes 01 ${filesPerJob} \$(ProcId) ${maxjob} ${output} ${analyzerTag} ${CMSSW_BASE} ${HOME}/" >> ${jdl_file}

        # option should always be 1, when running condor
        echo "Log = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_PC.log" >> ${jdl_file}
        echo "Output = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
        echo "Error = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}


        echo "RequestMemory = 4000" >> ${jdl_file}
        echo "RequestCpus = 1" >> ${jdl_file}
        echo "RequestDisk = 4" >> ${jdl_file}
        echo "+JobQueue=\"Short\"" >>  ${jdl_file}
        echo "+RunAsOwner = True" >> ${jdl_file}
        echo "+InteractiveUser = true" >> ${jdl_file}
        #echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7\"" >> ${jdl_file}
	echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7\"" >> ${jdl_file}
        echo '+SingularityBindCVMFS = True' >> ${jdl_file}
        echo "run_as_owner = True" >> ${jdl_file}
        echo "x509userproxy = ${HOME}/x509_proxy" >> ${jdl_file}
        echo "should_transfer_files = YES" >> ${jdl_file}
        echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
        echo "Queue ${maxjob}" >> ${jdl_file}
        echo "condor_submit ${jdl_file}"
        condor_submit ${jdl_file} -batch-name ${sample}


done
