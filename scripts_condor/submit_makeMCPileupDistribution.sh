#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

#job_script=${RazorAnalyzerDir}/scripts_condor/runMakeMCPileupDistribution.sh

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_llp_vH.sh
filesPerJob=20

for year in \
MC_Fall17
do
	if [ ${year} == "MC_Summer16" ]
        then
                tune=TuneCUETP8M1
        else
                tune=TuneCP5
        fi
        for sample in \
        ZToMuMu_NNPDF31_13TeV-powheg_M_50_120
	do
		echo "Sample " ${sample}
		version=/V1p17/${year}/v1/
		output=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${version}/pileUp
		#output=${CMSSW_BASE}/src/llp_analyzer/data/PileupWeights/
		inputfilelist=/src/llp_analyzer/lists/displacedJetMuonNtuple/${version}/sixie/${sample}.txt
		analyzer=MakeMCPileupDistribution
		nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
                maxjob=`python -c "print int($nfiles.0/$filesPerJob)+1"`
                mod=`python -c "print int($nfiles.0%$filesPerJob)"`
                if [ ${mod} -eq 0 ]
                then
                        maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
                fi
	
		rm -f submit/${analyzer}_${sample}_Job*.jdl
		rm -f log/${analyzer}_${sample}_Job*
		jdl_file=submit/${analyzer}_${sample}.jdl
		echo "Universe = vanilla" > ${jdl_file}
		echo "Executable = ${job_script}" >> ${jdl_file}
		echo "Arguments =  ${analyzer} ${inputfilelist} no 01 ${filesPerJob} \$(ProcId) ${maxjob} ${output} tag ${CMSSW_BASE} ${HOME}/" >> ${jdl_file}
		echo "Log = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_PC.log" >> ${jdl_file}
                echo "Output = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
                echo "Error = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}

		echo "Requirements=(TARGET.OpSysAndVer==\"CentOS7\" && regexp(\"blade.*\", TARGET.Machine))" >> ${jdl_file}
		echo "RequestMemory = 2000" >> ${jdl_file}
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
		echo "Queue ${maxjob}" >> ${jdl_file}
		echo "condor_submit ${jdl_file}"
		condor_submit ${jdl_file} --batch-name ${sample}
	done
done
