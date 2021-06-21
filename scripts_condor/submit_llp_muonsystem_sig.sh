#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_llp_vH.sh
filesPerJob=10

#ggH_HToSSTobbbb_ms55_pl1000_RunIIFall18
#ggH_HToSSTobbbb_ms55_pl1000 \
#ggH_HToSSTobbbb_ms55_pl10000
#n3n2-n1-hbb-hbb_mh127_pl1000_ev100000 \
#n3n2-n1-hbb-hbb_mh200_pl1000_ev100000 \
#n3n2-n1-hbb-hbb_mh400_pl1000_ev100000 \
#n3n2-n1-hinc-hgg_mh200_pl100_ev100000
#listSummer16=(
#ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch1
#ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch3
#ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch4
#ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch1
#ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch3
#ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch4
#WH_HToSSTobbbb_CSCDecayFilter_ms55_pl100000_ev150000
#)
listSummer16=(
ZH_HToSSTobbbb_ms55_pl1000_ev150000
ZH_HToSSTobbbb_ms55_pl10000_ev150000
WH_HToSSTobbbb_CSCDecayFilter_ms55_pl100000_ev150000
)
#listFall18=(
#WplusH_HToSSTobbbb_ms55_pl10000_ev150000
#WminusH_HToSSTobbbb_ms55_pl10000_ev150000
#WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
#)
listFall18=(
ggH_HToSS_SToEE_ms0p1_pl100
ggH_HToSS_SToEE_ms0p1_pl500
ggH_HToSS_SToPi0Pi0_ms1_pl100
ggH_HToSS_SToPi0Pi0_ms1_pl500
ggH_HToSSTobbbb_ms1_pl1000
ggH_HToSS_SToPiPlusPiMinus_ms1_pl500
ggH_HToSS_SToKPlusKMinus_ms1p5_pl500
ggH_HToSS_SToEE_ms0p4_pl500
)
#listFall18=(
##ttH_HToSS_SToBB_ms15_pl100
##ttH_HToSS_SToBB_ms15_pl1000
##ttH_HToSS_SToBB_ms15_pl10000
##ttH_HToSS_SToBB_ms15_pl100000
##ttH_HToSS_SToBB_ms40_pl100
##ttH_HToSS_SToBB_ms40_pl1000
##ttH_HToSS_SToBB_ms40_pl10000
##ttH_HToSS_SToBB_ms40_pl100000
##ttH_HToSS_SToBB_ms55_pl100
##ttH_HToSS_SToBB_ms55_pl1000
##ttH_HToSS_SToBB_ms55_pl10000
##ttH_HToSS_SToBB_ms55_pl100000
##WminusHToSS_SToBB_ms15_pl100
##WminusHToSS_SToBB_ms15_pl1000
##WminusHToSS_SToBB_ms15_pl10000
##WminusHToSS_SToBB_ms15_pl100000
##WminusHToSS_SToBB_ms40_pl100
##WminusHToSS_SToBB_ms40_pl1000
##WminusHToSS_SToBB_ms40_pl10000
##WminusHToSS_SToBB_ms40_pl100000
##WminusHToSS_SToBB_ms55_pl100
##WminusHToSS_SToBB_ms55_pl1000
##WminusHToSS_SToBB_ms55_pl10000
##WminusHToSS_SToBB_ms55_pl100000
##WplusHToSS_SToBB_ms15_pl100
##WplusHToSS_SToBB_ms15_pl1000
##WplusHToSS_SToBB_ms15_pl10000
##WplusHToSS_SToBB_ms15_pl100000
##WplusHToSS_SToBB_ms40_pl100
##WplusHToSS_SToBB_ms40_pl1000
##WplusHToSS_SToBB_ms40_pl10000
##WplusHToSS_SToBB_ms40_pl100000
##WplusHToSS_SToBB_ms55_pl100
#WplusHToSS_SToBB_ms55_pl1000
##WplusHToSS_SToBB_ms55_pl10000
##WplusHToSS_SToBB_ms55_pl100000
##ZHToSS_SToBB_ms15_pl100
##ZHToSS_SToBB_ms15_pl1000
##ZHToSS_SToBB_ms15_pl10000
##ZHToSS_SToBB_ms15_pl100000
##ZHToSS_SToBB_ms40_pl100
##ZHToSS_SToBB_ms40_pl1000
##ZHToSS_SToBB_ms40_pl10000
##ZHToSS_SToBB_ms40_pl100000
##ZHToSS_SToBB_ms55_pl100
##ZHToSS_SToBB_ms55_pl1000
##ZHToSS_SToBB_ms55_pl10000
##ZHToSS_SToBB_ms55_pl100000
#)
#listFall18=(
#ggZHToSS_SToBB_ZToQQ_ms15_pl100
#ggZHToSS_SToBB_ZToQQ_ms15_pl1000
#ggZHToSS_SToBB_ZToQQ_ms15_pl10000
#ggZHToSS_SToBB_ZToQQ_ms15_pl100000
#ggZHToSS_SToBB_ZToQQ_ms40_pl100
#ggZHToSS_SToBB_ZToQQ_ms40_pl1000
#ggZHToSS_SToBB_ZToQQ_ms40_pl10000
#ggZHToSS_SToBB_ZToQQ_ms40_pl100000
#ggZHToSS_SToBB_ZToQQ_ms55_pl100
#ggZHToSS_SToBB_ZToQQ_ms55_pl1000
#ggZHToSS_SToBB_ZToQQ_ms55_pl10000
#ggZHToSS_SToBB_ZToQQ_ms55_pl100000
#)
# listFall17=(
# ZToMuMu_NNPDF31_13TeV-powheg_M_120_200
# ZToMuMu_NNPDF31_13TeV-powheg_M_200_400
# ZToMuMu_NNPDF31_13TeV-powheg_M_3500_4500
# ZToMuMu_NNPDF31_13TeV-powheg_M_50_120
# ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300
# ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500
# ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000
# ZToMuMu_NNPDF31_13TeV-powheg_M_6000_Inf
# ZToMuMu_NNPDF31_13TeV-powheg_M_400_800
# )
# #listFall18=(
# #ggH_HToSSTobbbb_ms1_pl1000
# #)
listFall18=(
QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8
QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8
)
#listFall18=(
#TChiHH_mass127_pl1000
#TChiHH_mass150_pl1000
#TChiHH_mass200_pl1000
#TChiHH_mass300_pl1000
#TChiHH_mass400_pl1000
#)
#
#
#listBkg_Fall18=(
#WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
#)
listFall18=(
ggH_HToSS_SToEE_ms0p4_pl10
)
for year in \
Fall18
do
        echo ${year}
        sampleList=list${year}[@]
        #sampleList=listFall18[@]
	for sample in "${!sampleList}"
        do

		echo "Sample " ${sample}
		version=/V1p17/MC_${year}/v2/
		output=/storage/cms/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/${version}/v118/${sample}
		echo ${output}
	  	inputfilelist=/src/llp_analyzer/lists/displacedJetMuonNtuple/${version}/sixie/${sample}.txt
		#inputfilelist=/src/llp_analyzer/lists/displacedJetMuonNtuple//V1p17/MC_Fall18/v1//sixie/${sample}.txt
		nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
        	maxjob=`python -c "print int($nfiles.0/$filesPerJob)+1"`
        	mod=`python -c "print int($nfiles.0%$filesPerJob)"`
        	if [ ${mod} -eq 0 ]
        	then
        	        maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
        	fi
		analyzer=llp_MuonSystem_cluster
		if [[ ${year} == "Fall18" || ${year}=="Bkg_Fall18" ]]
        	then
        	        echo ${year}
        	        analyzerTag=Razor2018_17SeptEarlyReReco
        	elif [ ${year} == "Fall17" ]
        	then
        	        echo ${year}
        	        # analyzerTag=Razor2017_Source2018
                	analyzerTag=Razor2017_17Nov2017Rereco
        	elif [ ${year} == 'Summer16' ]
        	then
        	        echo ${year}
        	        # analyzerTag=Razor2016_Source2018
                	analyzerTag=Razor2016_07Aug2017Rereco
        	else
        	        echo "ERROR: NEED TO SET CORRECT YEAR"
        	fi
		echo ${analyzerTag}
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
		#echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7\"" >> ${jdl_file}
		echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7-m202006\"" >> ${jdl_file}
		echo '+SingularityBindCVMFS = True' >> ${jdl_file}
		echo "run_as_owner = True" >> ${jdl_file}
		echo "x509userproxy = ${HOME}/x509_proxy" >> ${jdl_file}
		echo "should_transfer_files = YES" >> ${jdl_file}
		echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
		echo "Queue ${maxjob}" >> ${jdl_file}
		echo "condor_submit ${jdl_file}"
		condor_submit ${jdl_file} --batch-name ${sample}
	done
done
