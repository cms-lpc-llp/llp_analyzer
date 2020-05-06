#!/bin/sh

mkdir -p log
mkdir -p submit

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runAnalyzerJob_SusyLLP.sh
filesPerJob=50

for sample in \
WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_ext2-v2_v4_v1 \

do
	echo "Sample " ${sample}
	outputDir=/store/group/phys_exotica/jmao/susy_llp/llp_analyzer/${sample} 
	inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetTimingNtuple/V1p10/MC_Summer16/${sample}.txt

	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
	maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
	mod=`python -c "print int($nfiles.0%$filesPerJob)"`
	analyzer=SusyLLP
	if [ ${mod} -eq 0 ]
        then
                maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
        fi

	for jobnumber in `seq 0 1 ${maxjob}`
	do
		
		jdl_file=submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
                #noFail=`grep YYYY log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}*.out`
		outRoot="/mnt/hadoop/${outputDir}/${sample}_Job${jobnumber}_Of_${maxjob}.root"
		
		minimumsize=1000
                actualsize=0
                if [ -f ${outRoot} ]
                then
                        actualsize=$(wc -c <"${outRoot}")
                fi
                if [ $actualsize -ge $minimumsize ]
		then
			finished=yes
                else
                        echo "job ${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob} failed, now being resubmitted"
                        rm -f log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}*
                        condor_submit submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
                fi
	done
done

