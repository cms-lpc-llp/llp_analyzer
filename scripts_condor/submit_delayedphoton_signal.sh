#!/bin/sh

mkdir -p log
mkdir -p submit

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_CaltechT2.sh
filesPerJob=1

for sample in \
GMSB_L100TeV_Ctau0_1cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau0_1cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L150TeV_Ctau800cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau0_1cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L200TeV_Ctau800cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L250TeV_Ctau800cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau0_1cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L300TeV_Ctau800cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau0_1cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L350TeV_Ctau800cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau0_1cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau10cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau1200cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau200cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau400cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau600cm_13TeV-pythia8 \
GMSB_L400TeV_Ctau800cm_13TeV-pythia8

do
	echo "Sample " ${sample}
	inputfilelist=/src/RazorAnalyzer/lists/Run2/razorNtuplerV4p1/MC_Summer16_reMINIAOD/${sample}.cern.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
	maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
	analyzer=DelayedPhotonAnalyzer

	rm submit/${analyzer}_${sample}_Job*.jdl
	rm log/${analyzer}_${sample}_Job*

	for jobnumber in `seq 0 1 ${maxjob}`
	do
		echo "job " ${jobnumber} " out of " ${maxjob}
		jdl_file=submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
		echo "Universe = vanilla" > ${jdl_file}
		echo "Executable = ${job_script}" >> ${jdl_file}
		echo "Arguments = ${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob} /store/group/phys_susy/razor/Run2Analysis/DelayedPhotonAnalysis/2016/orderByPt/jobs/ ${analyzer} ${inputfilelist} no 10 ${filesPerJob} ${jobnumber} ${sample}_Job${jobnumber}_Of_${maxjob}.root" >> ${jdl_file}
		echo "Log = log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}_PC.log" >> ${jdl_file}
		echo "Output = log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
		echo "Error = log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}
		echo 'Requirements=TARGET.OpSysAndVer=="CentOS7"' >> ${jdl_file}
		echo "should_transfer_files = YES" >> ${jdl_file}
		echo "RequestMemory = 2000" >> ${jdl_file}
		echo "RequestCpus = 1" >> ${jdl_file}
		echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
		echo "Queue 1" >> ${jdl_file}
		echo "condor_submit submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl"
		condor_submit submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
	done
done

