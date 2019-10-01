#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

mode=bkg
inputDir=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/V1p7/MC_Summer16/v11/v5/bkg/wH/
outputDir=${inputDir}normalized
job_script=${RazorAnalyzerDir}/scripts_condor/normalize.sh

#WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
#WJetsToLNu_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
#WJetsToLNu_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
#WJetsToLNu_Pt-400To600_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
#WJetsToLNu_Pt-600ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
#WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
for sample in \
WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
do
	echo "Sample " ${sample}
	analyzer=llp_MuonSystem
	rm -f submit/${analyzer}_normalize_${mode}_${sample}*.jdl
	rm -f log/${analyzer}_normalize_${mode}_${sample}*
	jdl_file=submit/${analyzer}_normalize_${mode}_${sample}.jdl
	echo "Universe = vanilla" > ${jdl_file}
	echo "Executable = ${job_script}" >> ${jdl_file}
	echo "Arguments =${mode} ${sample} ${inputDir} ${outputDir} ${CMSSW_BASE} ${HOME}/"  >> ${jdl_file}

	echo "Log = log/${analyzer}_normalize_${mode}_${sample}_PC.log" >> ${jdl_file}
	echo "Output = log/${analyzer}_normalize_${mode}_${sample}_\$(Cluster).\$(Process).out" >> ${jdl_file}
	echo "Error = log/${analyzer}_normalize_${mode}_${sample}_\$(Cluster).\$(Process).err" >> ${jdl_file}

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
	echo "x509userproxy = /data/christiw/x509_proxy" >> ${jdl_file}
	echo "should_transfer_files = YES" >> ${jdl_file}
	echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
	echo "Queue 1" >> ${jdl_file}
	echo "condor_submit ${jdl_file}"
	condor_submit ${jdl_file}
done
