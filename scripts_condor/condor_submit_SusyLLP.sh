#!/bin/sh

mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runAnalyzerJob_SusyLLP.sh
#filesPerJob=15  # qcd/znunu/ttjets bkg, remember to change the arguments 
filesPerJob=10  #bkg, remember to change the arguments 
#filesPerJob=1   #sig, remember to change the arguments
#filesPerJob=50  #x1n2 sig, remember to change the arguments

for sample in \
Run2016B-07Aug17 \

#n3n2-n1-hbb-hbb_mh127_pl1000_ev100000 \
#n3n2-n1-hbb-hbb_mh150_pl1000_ev100000 \
#n3n2-n1-hbb-hbb_mh175_pl1000_ev100000 \

#WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \

#ZJetsToNuNu_HT-100To200_13TeV-madgraph \
#ZJetsToNuNu_HT-200To400_13TeV-madgraph \
#ZJetsToNuNu_HT-400To600_13TeV-madgraph \
#ZJetsToNuNu_HT-600To800_13TeV-madgraph \
#ZJetsToNuNu_HT-800To1200_13TeV-madgraph \
#ZJetsToNuNu_HT-1200To2500_13TeV-madgraph \
#ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph \
#TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \

#QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \
#QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1_v2_v1 \

#x1n2-n1-wlv-hbb_mchi200_mlsp150_pl1000_ev100000 \

#WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6-v1_v2_v1 \

#WplusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-1000 \
#WminusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-1000 \
#x1n2-n1-wlv-hbb_mchi200_mlsp150_pl1000_ev100000 \
#n3n2-n1-zll-hbb_mh200_pl1000_ev100000 \

#WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6-v1_v1_v1 \
#WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_ext2-v2_v1_v1 \

#1
#RunIIFall17_x1n2-n1-wlv-hbb_mchi200_mlsp150_pl1000_ev100000_fullsim_signal_aod \
#RunIIFall17_x1n2-n1-wlv-hbb_mh200_pl1000_ev100000_fullsim_signal_aod \
#WminusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WplusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \

# 50

#RunIIFall17_n3n2-n1-hbb-hbb_mh200_pl1000_ev100000_fullsim_signal_aod \
#RunIIFall17_n3n2-n1-zll-hbb_mh200_pl1000_ev100000_fullsim_signal_aod \
#RunIIFall17_n3n2-n1-zll-hbb_mh200_pl100_ev100000_fullsim_signal_aod \

#QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \


do
	echo "Sample " ${sample}
	outputDir=/store/group/phys_exotica/jmao/susy_llp/llp_analyzer/${sample} 
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/llpntuple/susy_llp/${sample}.txt
	inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetTimingNtuple/V1p13/Data_2016/${sample}.txt
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetTimingNtuple/V1p11/MC_Fullsim/${sample}.txt
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetTimingNtuple/V1p10/MC_Summer16/${sample}.txt
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetTimingNtuple/V1p11/MC_Summer16/${sample}.txt
	#inputfilelist=/src/cms_lpc_llp/llp_analyzer/lists/displacedJetTimingNtuple/V1p13/MC_Summer16/${sample}.txt
	nfiles=`cat ${CMSSW_BASE}$inputfilelist | wc | awk '{print $1}' `
        maxjob=`python -c "print int($nfiles.0/$filesPerJob)"`
        mod=`python -c "print int($nfiles.0%$filesPerJob)"`
	echo "maxjob " ${maxjob}
	echo "Mod " ${mod}
        if [ ${mod} -eq 0 ]
        then
                maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
	fi
	echo "maxjob " ${maxjob}
	#lastjobfile=`python -c "print int($nfiles.0%$filesPerJob)"`
	analyzer=SusyLLP
	analyzerTag=Razor2016_80X
	rm -f submit/${analyzer}_${sample}_Job*.jdl
	rm -f log/${analyzer}_${sample}_Job*

	for jobnumber in `seq 0 1 ${maxjob}`
	do
		echo "job " ${jobnumber} " out of " ${maxjob}
		jdl_file=submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
		echo "Universe = vanilla" > ${jdl_file}
		echo "Executable = ${job_script}" >> ${jdl_file}
		#bkg option 1(condor)3(bkg_wh)1(pf=true)
		#echo "Arguments = ${analyzer} ${inputfilelist} no 131 ${filesPerJob} ${jobnumber} ${sample}_Job${jobnumber}_Of_${maxjob}.root ${outputDir} ${analyzerTag} " >> ${jdl_file}
		#signal option 1(condor)1(wh)1(pf=true)
		echo "Arguments = ${analyzer} ${inputfilelist} no 151 ${filesPerJob} ${jobnumber} ${sample}_Job${jobnumber}_Of_${maxjob}.root ${outputDir} ${analyzerTag} " >> ${jdl_file}

		# option should always be 1, when running condor
		echo "Log = log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}_PC.log" >> ${jdl_file}
		echo "Output = log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
		echo "Error = log/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}

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
		echo "x509userproxy = /storage/user/jmao/x509_proxy" >> ${jdl_file}
		echo "should_transfer_files = YES" >> ${jdl_file}
		echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
		echo "Queue 1" >> ${jdl_file}
		echo "condor_submit submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl"
		condor_submit submit/${analyzer}_${sample}_Job${jobnumber}_Of_${maxjob}.jdl
	done
done
