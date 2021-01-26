#!/bin/sh


mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/normalize.sh
#0 1 5 10 25 50 100 500 1000 2000 5000 10000 100000
ctau_VBF=(
0 1 5 10 25 50 100 500 1000 2000 5000 10000 100000
)
# 1 10 100 1000 10000 100000
#Summer16 \
#Fall17 \
#Fall18
ctau_ggH=(1 10 100 1000 10000 100000)
for year in \
Fall18_FullGenParticles
do
        if [ ${year} == "Summer16" ]
        then
                tune=TuneCUETP8M1
        else
                tune=TuneCP5
        fi
        if [ ${year} == "Fall18" ]
        then
                echo ${year}
                xsec=59740
        elif [ ${year} == "Fall17" ]
        then
                echo ${year}
                xsec=41530
        elif [ ${year} == 'Summer16' ]
        then
                echo ${year}
                xsec=35920
        else
                echo "ERROR: NEED TO SET CORRECT YEAR"
        fi
	xsec=137000
	tune=TuneCP5
	#VBFH_HToSSTo4Tau_MH-125_${tune}_13TeV-powheg-pythia8 \
	#VBFH_HToSSTo4b_MH-125_${tune}_13TeV-powheg-pythia8 \
	#ggH_HToSSTodddd_MH-125_${tune}_13TeV-powheg-pythia8
  #        ggH_HToSSTo4Tau_MH-125_${tune}_13TeV-powheg-pythia8 \

        #ggH_HToSSTobbbb_MH-125_${tune}_13TeV-powheg-pythia8
        for sample in \
        ggH_HToSSTodddd_MH-125_${tune}_13TeV-powheg-pythia8
	do

		mass=()
		if [[ ${sample} == "VBFH_HToSSTo4b_MH-125_${tune}_13TeV-powheg-pythia8" || ${sample} == "ggH_HToSSTobbbb_MH-125_${tune}_13TeV-powheg-pythia8" ]]
		then
			mass=( 15 40 55 )
		else
			mass=( 7 15 40 55 )
		fi
		for mx in "${mass[@]}"
		do
			ctau_list=ctau_${sample:0:3}[@]
                        for ctau in "${!ctau_list}"
                        do
				echo "Sample " ${sample} "mx " ${mx} "ctau " ${ctau}
				sample_long=${sample%_MH125*}_MH-125_MS-${mx}_ctau-${ctau}_${tune}_13TeV-powheg-pythia8
				analyzer=llp_MuonSystem_bdt
				inputDir=/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/V1p17/MC_${year}/v1/v102/
				outputDir=${inputDir}normalized
				echo ${inputDir}
				rm -f submit/${analyzer}_normalize_${year}_${sample_long}.jdl
				rm -f log/${analyzer}_normalize_${year}_${sample_long}*
				jdl_file=submit/${analyzer}_normalize_${year}_${sample_long}.jdl

				echo "Universe = vanilla" > ${jdl_file}
				echo "Executable = ${job_script}" >> ${jdl_file}
				echo "Arguments = central no ${sample}_mx${mx}_ctau${ctau} ${inputDir} ${outputDir} ${CMSSW_BASE} ${HOME}/ ${xsec}"  >> ${jdl_file}
				#echo "Arguments = central no ${sample}_mx${mx}_ctau${ctau} ${inputDir} ${outputDir} ${CMSSW_BASE} ${HOME}/ 1"  >> ${jdl_file}

				echo "Log = log/${analyzer}_normalize_${year}_${sample_long}_PC.log" >> ${jdl_file}
				echo "Output = log/${analyzer}_normalize_${year}_${sample_long}_\$(Cluster).\$(Process).out" >> ${jdl_file}
				echo "Error = log/${analyzer}_normalize_${year}_${sample_long}_\$(Cluster).\$(Process).err" >> ${jdl_file}

				#echo "Requirements=TARGET.OpSysAndVer==\"CentOS7\"" >> ${jdl_file}
				echo "Requirements=(TARGET.OpSysAndVer==\"CentOS7\" && regexp(\"blade.*\", TARGET.Machine))" >> ${jdl_file}
				echo "RequestMemory = 4000" >> ${jdl_file}
				echo "RequestCpus = 1" >> ${jdl_file}
				echo "RequestDisk = 4" >> ${jdl_file}

				echo "+RunAsOwner = True" >> ${jdl_file}
				echo "+InteractiveUser = true" >> ${jdl_file}
				#echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7\"" >> ${jdl_file}
				echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7-m202006\"" >> ${jdl_file}
				echo '+SingularityBindCVMFS = True' >> ${jdl_file}
				echo "run_as_owner = True" >> ${jdl_file}
				echo "x509userproxy = ${HOME}/x509_proxy" >> ${jdl_file}
				echo "should_transfer_files = YES" >> ${jdl_file}
				echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
				echo "Queue 1" >> ${jdl_file}
				echo "condor_submit ${jdl_file}"
				condor_submit ${jdl_file}
			done
		done
	done
done
