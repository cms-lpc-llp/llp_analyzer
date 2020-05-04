#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runMakeMCPileupDistribution.sh

for year in \
MC_Fall18 \
MC_Fall17 \
MC_Summer16
do
	if [ ${year} == "MC_Summer16" ]
        then
                tune=TuneCUETP8M1
        else
                tune=TuneCP5
        fi
        for sample in \
        ggH_HToSSTobbbb_MH-125_${tune}_13TeV-powheg-pythia8
        do
		echo "Sample " ${sample}
		version=/V1p15/${year}/v3/
		output=${CMSSW_BASE}/src/llp_analyzer/data/PileupWeights/
		inputfilelist=/src/llp_analyzer/lists/displacedJetMuonNtuple/${version}/sixie/${sample}.txt
		analyzer=MakeMCPileupDistribution
		rm -f submit/${analyzer}_${sample}_Job*.jdl
		rm -f log/${analyzer}_${sample}_Job*
		jdl_file=submit/${analyzer}_${sample}.jdl
		echo "Universe = vanilla" > ${jdl_file}
		echo "Executable = ${job_script}" >> ${jdl_file}
		echo "Arguments = ${inputfilelist} PileupSource_${year}_${sample}.root ${output} ${CMSSW_BASE} ${HOME}/" >> ${jdl_file}
		echo "Log = log/${analyzer}_${year}_${sample}_PC.log" >> ${jdl_file}
		echo "Output = log/${analyzer}_${year}_${sample}_\$(Cluster).\$(Process).out" >> ${jdl_file}
		echo "Error = log/${analyzer}_${year}_${sample}_\$(Cluster).\$(Process).err" >> ${jdl_file}
	
		#echo "Requirements=TARGET.OpSysAndVer==\"CentOS7\"" >> ${jdl_file}
		echo "Requirements=(TARGET.OpSysAndVer==\"CentOS7\" && regexp(\"blade.*\", TARGET.Machine))" >> ${jdl_file}
		echo "RequestMemory = 4000" >> ${jdl_file}
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
		echo "Queue 1" >> ${jdl_file}
		echo "condor_submit ${jdl_file}"
		condor_submit ${jdl_file}
	done
done
