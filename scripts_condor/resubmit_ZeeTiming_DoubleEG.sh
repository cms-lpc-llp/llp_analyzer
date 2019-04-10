#!/bin/sh

mkdir -p log
mkdir -p submit

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_CaltechT2.sh
filesPerJob=15

for sample in \
DoubleEG_2016B_ver1_06Aug2018 \
DoubleEG_2016B_ver2_06Aug2018 \
DoubleEG_2016C_06Aug2018 \
DoubleEG_2016D_06Aug2018 \
DoubleEG_2016E_06Aug2018 \
DoubleEG_2016F_06Aug2018 \
DoubleEG_2016G_06Aug2018 \
DoubleEG_2016H_06Aug2018

do
	echo "Sample " ${sample}
	inputfilelist=/src/RazorAnalyzer/lists/Run2/razorNtuplerV4p1/Data_2016_reMINIAOD/${sample}.cern.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
	maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
	analyzer=ZeeTiming
	#rm submit/${sample}_Job*.jdl
	#rm log/${sample}_Job*

	for jobnumber in `seq 0 1 ${maxjob}`
	do
		jdl_file=submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
		#noFail=`grep YYYY log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}*.out`
		outRoot="/mnt/hadoop/store/group/phys_susy/razor/Run2Analysis/EcalTiming/ntuples_V4p1_31Aug2018/jobs/${sample}_Job${jobnumber}_Of_${maxjob}.root"

		minimumsize=50000
                actualsize=0
                if [ -f ${outRoot} ]
                then
                        actualsize=$(wc -c <"${outRoot}")
                fi
                if [ $actualsize -ge $minimumsize ]
                then
                        echo "job ${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob} finished already "
		#elif [ -z "${noFail}" ]
		#then
		#	echo "job ${sample}_Job${jobnumber}_Of_${maxjob} being processed now, be patient"
		else
			echo "job ${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob} failed, now being resubmitted"
			rm ${outRoot}
			#rm log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}*
			#condor_submit submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
		fi
	done
done

