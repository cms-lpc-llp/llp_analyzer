#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

mode=bkg
isData=no
#inputDir=/store/group/phys_exotica/jmao/susy_llp/llp_analyzer/
inputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p16/v2/MC_Summer16/
outputDir=${inputDir}normalized
job_script=${RazorAnalyzerDir}/scripts_condor/normalize.sh

for sample in \
WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \

#WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \

#QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \

#ZJetsToNuNu_HT-100To200_13TeV-madgraph \
#ZJetsToNuNu_HT-200To400_13TeV-madgraph \
#ZJetsToNuNu_HT-400To600_13TeV-madgraph \
#ZJetsToNuNu_HT-600To800_13TeV-madgraph \
#ZJetsToNuNu_HT-800To1200_13TeV-madgraph \
#ZJetsToNuNu_HT-1200To2500_13TeV-madgraph \
#ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph \
#TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#TTJets_SingleLeptFromTbar_genMET-150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \

#WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \

#QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#ZJetsToNuNu_HT-100To200_13TeV-madgraph \
#ZJetsToNuNu_HT-200To400_13TeV-madgraph \
#ZJetsToNuNu_HT-400To600_13TeV-madgraph \
#ZJetsToNuNu_HT-600To800_13TeV-madgraph \
#ZJetsToNuNu_HT-800To1200_13TeV-madgraph \
#ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph \
#TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \

#ZJetsToNuNu_HT-1200To2500_13TeV-madgraph \

do
	echo "Sample " ${sample}
	analyzer=SusyLLP
	rm -f submit/${analyzer}_normalize_${mode}_${sample}*.jdl
	rm -f log/${analyzer}_normalize_${mode}_${sample}*
	jdl_file=submit/${analyzer}_normalize_${mode}_${sample}.jdl
	echo "Universe = vanilla" > ${jdl_file}
	echo "Executable = ${job_script}" >> ${jdl_file}
	echo "Arguments =${mode} ${isData} ${sample} ${inputDir} ${outputDir} ${CMSSW_BASE} ${HOME}/"  >> ${jdl_file}

	echo "Log = log/${analyzer}_normalize_${mode}_${sample}_PC.log" >> ${jdl_file}
	echo "Output = log/${analyzer}_normalize_${mode}_${sample}_\$(Cluster).\$(Process).out" >> ${jdl_file}
	echo "Error = log/${analyzer}_normalize_${mode}_${sample}_\$(Cluster).\$(Process).err" >> ${jdl_file}

	#echo "Requirements=TARGET.OpSysAndVer==\"CentOS7\"" >> ${jdl_file}
	echo "Requirements=(TARGET.OpSysAndVer==\"CentOS7\" && regexp(\"blade.*\", TARGET.Machine))" >> ${jdl_file}
	echo "RequestMemory = 20000" >> ${jdl_file}
	echo "RequestCpus = 1" >> ${jdl_file}
	echo "RequestDisk = 4" >> ${jdl_file}

	echo "+RunAsOwner = True" >> ${jdl_file}
	echo "+InteractiveUser = true" >> ${jdl_file}
	echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7\"" >> ${jdl_file}
	echo '+SingularityBindCVMFS = True' >> ${jdl_file}
	echo "run_as_owner = True" >> ${jdl_file}
	echo "x509userproxy = /storage/user/jmao/x509_proxy" >> ${jdl_file}
	echo "should_transfer_files = YES" >> ${jdl_file}
	echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
	echo "Queue 1" >> ${jdl_file}
	echo "condor_submit ${jdl_file}"
	condor_submit ${jdl_file} -batch-name ${sample}
done
