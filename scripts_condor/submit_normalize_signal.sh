#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

mode=wH
inputDir=/store/group/phys_exotica/delayedjets/llp_analyzer/V1p0/MC_Summer16/v1/signals/
job_script=${RazorAnalyzerDir}/scripts_condor/normalize_signal.sh

for sample in \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh1000_mx475_pl10000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh1000_mx475_pl1000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh1000_mx475_pl500_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl10000_ev100000 \
ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl500_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh2000_mx975_pl10000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh2000_mx975_pl1000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh2000_mx975_pl500_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh300_mx125_pl10000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh300_mx125_pl1000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh300_mx125_pl500_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh500_mx225_pl10000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh500_mx225_pl1000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh500_mx225_pl500_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh1000_mx475_pl10000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh1000_mx475_pl1000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh1000_mx475_pl500_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh125_mx50_pl10000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh125_mx50_pl1000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh125_mx50_pl500_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh2000_mx975_pl10000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh2000_mx975_pl1000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh2000_mx975_pl500_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh300_mx125_pl10000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh300_mx125_pl1000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh300_mx125_pl500_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh500_mx225_pl10000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh500_mx225_pl1000_ev100000 \
#ppTohToSS1SS2_SS1Tobb_SS2Toveve_vh_withISR_mh500_mx225_pl500_ev100000

do
	echo "Sample " ${sample}
	analyzer=llp_vH
	analyzerTag=Razor2016_80X
	rm -f submit/${analyzer}_normalize_${mode}_${sample}_Job*.jdl
	rm -f log/${analyzer}_${sample}_Job*

	jdl_file=submit/${analyzer}_normalize_${mode}_${sample}.jdl
	echo "Universe = vanilla" > ${jdl_file}
	echo "Executable = ${job_script}" >> ${jdl_file}
	echo "Arguments = wH ${sample} ${inputDir}"  >> ${jdl_file}

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
