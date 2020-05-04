#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_llp_vH.sh
filesPerJob=20

#ggH_HToSSTobbbb_ms55_pl1000_RunIIFall18
#ggH_HToSSTobbbb_ms55_pl1000 \
#ggH_HToSSTobbbb_ms55_pl10000
#n3n2-n1-hbb-hbb_mh127_pl1000_ev100000 \
#n3n2-n1-hbb-hbb_mh200_pl1000_ev100000 \
#n3n2-n1-hbb-hbb_mh400_pl1000_ev100000 \
#n3n2-n1-hinc-hgg_mh200_pl100_ev100000
for sample in \
ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch1 \
ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch3 \
ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch4 \
ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch1 \
ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch3 \
ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch4
do
#for sample in \
#WplusH_HToSSTobbbb_ms55_pl10000_ev150000 \
#WminusH_HToSSTobbbb_ms55_pl10000_ev150000
#do
#for sample in \
#WH_HToSSTobbbb_CSCDecayFilter_ms55_pl100000_ev150000
#do
	echo "Sample " ${sample}
	year=Summer16
	#year=Fall18
	version=/V1p16/MC_${year}/v1/
	output=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${version}/v8/${sample}
	echo ${output}
	inputfilelist=/src/llp_analyzer/lists/displacedJetMuonNtuple/${version}/sixie/${sample}.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
        maxjob=`python -c "print int($nfiles.0/$filesPerJob)+1"`
        mod=`python -c "print int($nfiles.0%$filesPerJob)"`
        if [ ${mod} -eq 0 ]
        then
                maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
        fi
	analyzer=llp_MuonSystem_cluster
	if [ ${year} == "Fall18" ]
        then
                echo ${year}
                analyzerTag=Razor2018_17SeptEarlyReReco
        elif [ ${year} == "Fall17" ]
        then
                echo ${year}
                analyzerTag=Razor2017_17Nov2017Rereco
        elif [ ${year} == 'Summer16' ]
        then
                echo ${year}
                analyzerTag=Razor2016_07Aug2017Rereco
        else
                echo "ERROR: NEED TO SET CORRECT YEAR"
        fi
	rm -f submit/${analyzer}_${sample}_Job*.jdl
	rm -f log/${analyzer}_${sample}_Job*

	echo "number of jobs: " ${maxjob}
	jdl_file=submit/${analyzer}_${sample}_${maxjob}.jdl
	echo "Universe = vanilla" > ${jdl_file}
	echo "Executable = ${job_script}" >> ${jdl_file}
	echo "Arguments = ${analyzer} ${inputfilelist} no 01 ${filesPerJob} \$(ProcId) ${maxjob} ${output} ${analyzerTag} ${CMSSW_BASE} ${HOME}/" >> ${jdl_file}

	# option should always be 1, when running condor
	echo "Log = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_PC.log" >> ${jdl_file}
	echo "Output = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
	echo "Error = log/${analyzer}_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}

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
	echo "x509userproxy = ${HOME}/x509_proxy" >> ${jdl_file}
	echo "should_transfer_files = YES" >> ${jdl_file}
	echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
	echo "Queue ${maxjob}" >> ${jdl_file}
	echo "condor_submit ${jdl_file}"
	condor_submit ${jdl_file} --batch-name ${sample}
done
