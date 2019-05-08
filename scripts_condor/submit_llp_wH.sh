#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_llp_vH.sh
filesPerJob=1

for sample in \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh1000_mx475_pl10000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh1000_mx475_pl1000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh1000_mx475_pl500_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl10000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl500_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh2000_mx975_pl10000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh2000_mx975_pl1000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh2000_mx975_pl500_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh300_mx125_pl10000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh300_mx125_pl1000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh300_mx125_pl500_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh500_mx225_pl10000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh500_mx225_pl1000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh500_mx225_pl500_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh1000_mx475_pl10000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh1000_mx475_pl1000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh1000_mx475_pl500_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh125_mx50_pl10000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh125_mx50_pl1000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh125_mx50_pl500_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh2000_mx975_pl10000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh2000_mx975_pl1000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh2000_mx975_pl500_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh300_mx125_pl10000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh300_mx125_pl1000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh300_mx125_pl500_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh500_mx225_pl10000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh500_mx225_pl1000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh500_mx225_pl500_ev100000

do
	echo "Sample " ${sample}
	inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/llpntuple/V1p0/MC_Summer16/v1/christiw/${sample}.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
	maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
	echo "new sample is ${sample}"
	mod=`python -c "print int($nfiles.0%$filesPerJob)"`
	if [ ${mod} -eq 0 ]
	then
		maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
	fi
	analyzer=llp_vH
	analyzerTag=Razor2016_80X
	rm -f submit/${analyzer}_${sample}_Job*.jdl
	rm -f log/${analyzer}_${sample}_Job*

	for jobnumber in `seq 0 1 ${maxjob}`
	do
		echo "job " ${jobnumber} " out of " ${maxjob}
		jdl_file=submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
		echo "Universe = vanilla" > ${jdl_file}
		echo "Executable = ${job_script}" >> ${jdl_file}
		echo "Arguments = ${analyzer} ${inputfilelist} no 11 ${filesPerJob} ${jobnumber} ${sample}_Job${jobnumber}_Of_${maxjob}.root /store/group/phys_exotica/delayedjets/llp_analyzer/V1p0/MC_Summer16/v3/signals/wH/${sample} ${analyzerTag} " >> ${jdl_file}

		# option should always be 1, when running condor
		echo "Log = log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}_PC.log" >> ${jdl_file}
		echo "Output = log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
		echo "Error = log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}

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
		echo "condor_submit submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl"
		condor_submit submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
	done
done
