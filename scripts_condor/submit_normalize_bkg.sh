#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

mode=bkg
outputDir=/store/group/phys_exotica/delayedjets/llp_analyzer/V1p6/MC_Summer16/v1/bkg/wH/normalized/
inputDir=/store/group/phys_exotica/delayedjets/llp_analyzer/V1p6/MC_Summer16/v1/bkg/wH/
job_script=${RazorAnalyzerDir}/scripts_condor/normalize.sh

for sample in \
ZJetsToNuNu_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8 \
QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1 \
ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1 \
ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1 \
ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1 \
TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
do
	echo "Sample " ${sample}
	analyzer=llp_vH
	rm -f submit/${analyzer}_normalize_${mode}_${sample}*.jdl
	rm -f log/${analyzer}_normalize_${mode}_${sample}*

	jdl_file=submit/${analyzer}_normalize_${mode}_${sample}.jdl
	echo "Universe = vanilla" > ${jdl_file}
	echo "Executable = ${job_script}" >> ${jdl_file}
	echo "Arguments =${mode} ${sample} ${inputDir} ${outputDir}"  >> ${jdl_file}

	echo "Log = log/${analyzer}_normalize_${mode}_${sample}_PC.log" >> ${jdl_file}
	echo "Output = log/${analyzer}_normalize_${mode}_${sample}_\$(Cluster).\$(Process).out" >> ${jdl_file}
	echo "Error = log/${analyzer}_normalize_${mode}_${sample}_\$(Cluster).\$(Process).err" >> ${jdl_file}

	echo "Requirements=TARGET.OpSysAndVer==\"CentOS7\"" >> ${jdl_file}
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
