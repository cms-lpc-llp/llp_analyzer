#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_llp_vH.sh
filesPerJob=5

for sample in \
TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
do
	echo "Sample " ${sample}
	#output=/storage/user/christiw/displacedJetMuonAnalyzer/V1p7/MC_Summer16/v3/bkg/wH/${sample}
	version=/V1p13/MC_Summer16/v6/jmao/
	output=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/V1p12/MC_Summer16/v4/v3/${sample}
	echo ${output}
	inputfilelist=/src/llp_analyzer/lists/displacedJetTimingNtuple/${version}/${sample}.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
        maxjob=`python -c "print int($nfiles.0/$filesPerJob)+1"`
        mod=`python -c "print int($nfiles.0%$filesPerJob)"`
        if [ ${mod} -eq 0 ]
        then
                maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
        fi
	analyzer=llp_MuonSystem
	analyzerTag=Razor2016_07Aug2017Rereco
	rm -f submit/${analyzer}_${sample}_Job*.jdl
	rm -f log/${analyzer}_${sample}_Job*

	echo "job " ${maxjob}
	jdl_file=submit/${analyzer}_${sample}_${maxjob}.jdl
	echo "Universe = vanilla" > ${jdl_file}
	echo "Executable = ${job_script}" >> ${jdl_file}
	echo "Arguments = ${analyzer} ${inputfilelist} no 1 ${filesPerJob} \$(ProcId) ${maxjob} ${output} ${analyzerTag} ${CMSSW_BASE} ${HOME}/" >> ${jdl_file}

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
