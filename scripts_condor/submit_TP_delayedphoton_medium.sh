#!/bin/sh

mkdir -p log
mkdir -p submit

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_CaltechT2.sh
filesPerJob=50

for sample in \
SingleElectron_2016B_ver1_06Aug2018 \
SingleElectron_2016B_ver2_06Aug2018 \
SingleElectron_2016C_06Aug2018 \
SingleElectron_2016D_06Aug2018 \
SingleElectron_2016E_06Aug2018 \
SingleElectron_2016F_06Aug2018 \
SingleElectron_2016G_06Aug2018 \
SingleElectron_2016H_06Aug2018

do
	echo "Sample " ${sample}
	inputfilelist=/src/RazorAnalyzer/lists/Run2/razorNtuplerV4p1/Data_2016_reMINIAOD/${sample}.cern.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
	maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
	analyzer=RazorTagAndProbe
	TP_tag=PhoMediumEffDenominatorReco

	rm submit/${analyzer}_${TP_tag}_${sample}_Job*.jdl
	rm log/${analyzer}_${TP_tag}_${sample}_Job*

	for jobnumber in `seq 0 1 ${maxjob}`
	do
		echo "job " ${jobnumber} " out of " ${maxjob}
		jdl_file=submit/${analyzer}_${TP_tag}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
		echo "Universe = vanilla" > ${jdl_file}
		echo "Executable = ${job_script}" >> ${jdl_file}
		echo "Arguments = ${analyzer}_${TP_tag}_${sample}_Job${jobnumber}_Of_${maxjob} /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/${TP_tag}/2016/jobs ${analyzer} ${inputfilelist} yes 30115 ${filesPerJob} ${jobnumber} ${sample}_Job${jobnumber}_Of_${maxjob}.root" >> ${jdl_file}
		echo "Log = log/${analyzer}_${TP_tag}_${sample}_Job${jobnumber}_Of_${maxjob}_PC.log" >> ${jdl_file}
		echo "Output = log/${analyzer}_${TP_tag}_${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
		echo "Error = log/${analyzer}_${TP_tag}_${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}
		echo 'Requirements=TARGET.OpSysAndVer=="CentOS7"' >> ${jdl_file}
		echo "should_transfer_files = YES" >> ${jdl_file}
		echo "RequestMemory = 2000" >> ${jdl_file}
		echo "RequestCpus = 1" >> ${jdl_file}
		echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
		echo "Queue 1" >> ${jdl_file}
		echo "condor_submit submit/${analyzer}_${TP_tag}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl"
		condor_submit submit/${analyzer}_${TP_tag}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
	done
done

