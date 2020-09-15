#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

mode=bkg
isData=no
inputDir=/store/group/phys_exotica/delayedjets/displacedJetTimingAnalyzer/V1p17/v14/MC_Fall18/

#outputDir=${inputDir}normalized
outputDir=${inputDir}1pb_normalized
job_script=${RazorAnalyzerDir}/scripts_condor/normalize.sh

for sample in \
GluGluH2_H2ToSSTobbbb_MH-400_MS-100_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
GluGluH2_H2ToSSTobbbb_MH-400_MS-50_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
GluGluH2_H2ToSSTobbbb_MH-600_MS-150_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
GluGluH2_H2ToSSTobbbb_MH-600_MS-50_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
GluGluH2_H2ToSSTobbbb_MH-1000_MS-400_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
GluGluH2_H2ToSSTobbbb_MH-1000_MS-150_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
GluGluH2_H2ToSSTobbbb_MH-1500_MS-500_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
GluGluH2_H2ToSSTobbbb_MH-1500_MS-200_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
GluGluH2_H2ToSSTobbbb_MH-2000_MS-600_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
GluGluH2_H2ToSSTobbbb_MH-2000_MS-250_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \

#GluGluH2_H2ToSSTobbbb_MH-600_MS-150_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-600_MS-150_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-600_MS-150_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-600_MS-150_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-600_MS-150_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-600_MS-50_ctauS-10000_${tune}_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-600_MS-50_ctauS-1000_${tune}_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-600_MS-50_ctauS-2000_${tune}_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-600_MS-50_ctauS-5000_${tune}_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-600_MS-50_ctauS-500_${tune}_13TeV-pythia8_PRIVATE-MC \

#GluGluH2_H2ToSSTobbbb_MH-1000_MS-150_ctauS-10000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1000_MS-150_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1000_MS-150_ctauS-2000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1000_MS-150_ctauS-5000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1000_MS-150_ctauS-500_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1000_MS-400_ctauS-10000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1000_MS-400_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1000_MS-400_ctauS-2000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1000_MS-400_ctauS-5000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1000_MS-400_ctauS-500_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1500_MS-200_ctauS-10000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1500_MS-200_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1500_MS-200_ctauS-2000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1500_MS-200_ctauS-5000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1500_MS-200_ctauS-500_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1500_MS-500_ctauS-10000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1500_MS-500_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1500_MS-500_ctauS-2000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1500_MS-500_ctauS-5000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-1500_MS-500_ctauS-500_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-2000_MS-250_ctauS-10000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-2000_MS-250_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-2000_MS-250_ctauS-2000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-2000_MS-250_ctauS-5000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-2000_MS-250_ctauS-500_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-2000_MS-600_ctauS-10000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-2000_MS-600_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-2000_MS-600_ctauS-2000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-2000_MS-600_ctauS-5000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-2000_MS-600_ctauS-500_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-200_MS-25_ctauS-10000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-200_MS-25_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-200_MS-25_ctauS-2000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-200_MS-25_ctauS-5000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-200_MS-25_ctauS-500_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-200_MS-50_ctauS-10000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-200_MS-50_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-200_MS-50_ctauS-2000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-200_MS-50_ctauS-5000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-200_MS-50_ctauS-500_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-400_MS-100_ctauS-10000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-400_MS-100_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-400_MS-100_ctauS-2000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-400_MS-100_ctauS-5000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-400_MS-100_ctauS-500_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-400_MS-50_ctauS-10000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-400_MS-50_ctauS-1000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-400_MS-50_ctauS-2000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-400_MS-50_ctauS-5000_TuneCP5_13TeV-pythia8_PRIVATE-MC \
#GluGluH2_H2ToSSTobbbb_MH-400_MS-50_ctauS-500_TuneCP5_13TeV-pythia8_PRIVATE-MC \

do
	echo "Sample " ${sample}
	analyzer=SusyLLP
	#analyzer=SusyLLPPassFail
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
	#echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7\"" >> ${jdl_file}
	echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7-m202006\"" >> ${jdl_file}
	echo '+SingularityBindCVMFS = True' >> ${jdl_file}
	echo "run_as_owner = True" >> ${jdl_file}
	echo "x509userproxy = /storage/user/jmao/x509_proxy" >> ${jdl_file}
	echo "should_transfer_files = YES" >> ${jdl_file}
	echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
	echo "Queue 1" >> ${jdl_file}
	echo "condor_submit ${jdl_file}"
	condor_submit ${jdl_file} -batch-name ${sample}
done
