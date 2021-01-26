#!/bin/bash

#WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
#WJetsToLNu_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
#WJetsToLNu_Pt-400To600_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
#WJetsToLNu_Pt-600ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
#WminusH_HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WminusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WminusH_HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WplusH_HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WplusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WplusH_HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WminusH_HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WminusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WminusH_HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WplusH_HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WplusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WplusH_HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WminusH_HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WminusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WminusH_HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WplusH_HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WplusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8 \
#WplusH_HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8
#DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
#QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
#MC_Summer16 \
#MC_Fall17 \
for year in \
MC_Fall18
do
	version=displacedJetMuonNtuple/V1p17/${year}/v2/sixie/
	root_dir=/mnt/hadoop/store/group/phys_exotica/delayedjets/${version}/
	list_dir=$CMSSW_BASE/src/llp_analyzer/lists/${version}
	echo $list_dir
	mkdir -p $list_dir

	#ggH_HToSSTobbbb_ms55_pl1000_RunIIFall18
	#n3n2-n1-hbb-hbb_mh127_pl1000_ev100000 \
	#n3n2-n1-hbb-hbb_mh200_pl1000_ev100000 \
	#n3n2-n1-hbb-hbb_mh400_pl1000_ev100000 \
	#n3n2-n1-hinc-hgg_mh200_pl100_ev100000
	#ggH_HToSSTobbbb_ms1_pl1000
	#WplusH_HToSSTobbbb_ms55_pl10000_ev150000 \
	#WminusH_HToSSTobbbb_ms55_pl10000_ev150000
	#for sample in \
	#ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch1 \
	#ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch3 \
	#ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch4 \
	#ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch1 \
	#ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch3 \
	#ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch4
	#do
	#for sample in \
	#WH_HToSSTobbbb_CSCDecayFilter_ms55_pl100000_ev150000 \
	#ZH_HToSSTobbbb_ms55_pl1000_ev150000 \
	#ZH_HToSSTobbbb_ms55_pl10000_ev150000
	#do
	#ggH_HToSSTobbbb_ms1_pl1000 \
	#WplusH_HToSSTobbbb_ms55_pl10000_ev150000 \
        #WminusH_HToSSTobbbb_ms55_pl10000_ev150000
	#WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
	#WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
	#WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
	#WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
	#WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
	#WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
	#WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
	#WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
	#for sample in \
	#ggH_HToSSTobbbb_ms1_pl1000 \
	#ggH_HToSS_SToEE_ms0p1_pl100 \
	#ggH_HToSS_SToEE_ms0p1_pl500 \
	#ggH_HToSS_SToPi0Pi0_ms1_pl100 \
	#ggH_HToSS_SToPi0Pi0_ms1_pl500 \
	#ggH_HToSS_SToPiPlusPiMinus_ms1_pl500 \
	#ggH_HToSS_SToKPlusKMinus_ms1p5_pl500 \
	#ggH_HToSS_SToEE_ms0p4_pl500
	#do
	#for sample in \
	#ZToMuMu_NNPDF31_13TeV-powheg_M_50_120 \
	#ZToMuMu_NNPDF31_13TeV-powheg_M_120_200 \
	#ZToMuMu_NNPDF31_13TeV-powheg_M_200_400 \
	#ZToMuMu_NNPDF31_13TeV-powheg_M_400_800 \
	#ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300 \
	#ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500 \
	#ZToMuMu_NNPDF31_13TeV-powheg_M_3500_4500 \
	#ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000 \
	#ZToMuMu_NNPDF31_13TeV-powheg_M_6000_Inf
	#do
	#ttH_HToSS_SToBB_ms15_pl100 \
	#ttH_HToSS_SToBB_ms15_pl1000 \
	#ttH_HToSS_SToBB_ms15_pl10000 \
	#ttH_HToSS_SToBB_ms15_pl100000 \
	#ttH_HToSS_SToBB_ms40_pl100 \
	#ttH_HToSS_SToBB_ms40_pl1000 \
	#ttH_HToSS_SToBB_ms40_pl10000 \
	#ttH_HToSS_SToBB_ms40_pl100000 \
	#ttH_HToSS_SToBB_ms55_pl100 \
	#ttH_HToSS_SToBB_ms55_pl1000 \
	#ttH_HToSS_SToBB_ms55_pl10000 \
	#ttH_HToSS_SToBB_ms55_pl100000 \
	#WminusHToSS_SToBB_ms15_pl100 \
	#WminusHToSS_SToBB_ms15_pl1000 \
	#WminusHToSS_SToBB_ms15_pl10000 \
	#WminusHToSS_SToBB_ms15_pl100000 \
	#WminusHToSS_SToBB_ms40_pl100 \
	#WminusHToSS_SToBB_ms40_pl1000 \
	#WminusHToSS_SToBB_ms40_pl10000 \
	#WminusHToSS_SToBB_ms40_pl100000 \
	#WminusHToSS_SToBB_ms55_pl100 \
	#WminusHToSS_SToBB_ms55_pl1000 \
	#WminusHToSS_SToBB_ms55_pl10000 \
	#WminusHToSS_SToBB_ms55_pl100000 \
	#WplusHToSS_SToBB_ms15_pl100 \
	#WplusHToSS_SToBB_ms15_pl1000 \
	#WplusHToSS_SToBB_ms15_pl10000 \
	#WplusHToSS_SToBB_ms15_pl100000 \
	#WplusHToSS_SToBB_ms40_pl100 \
	#WplusHToSS_SToBB_ms40_pl1000 \
	#WplusHToSS_SToBB_ms40_pl10000 \
	#WplusHToSS_SToBB_ms40_pl100000 \
	#WplusHToSS_SToBB_ms55_pl100 \
	#WplusHToSS_SToBB_ms55_pl1000 \
	#WplusHToSS_SToBB_ms55_pl10000 \
	#WplusHToSS_SToBB_ms55_pl100000 \
	#ZHToSS_SToBB_ms15_pl100 \
	#ZHToSS_SToBB_ms15_pl1000 \
	#ZHToSS_SToBB_ms15_pl10000 \
	#ZHToSS_SToBB_ms15_pl100000 \
	#ZHToSS_SToBB_ms40_pl100 \
	#ZHToSS_SToBB_ms40_pl1000 \
	#ZHToSS_SToBB_ms40_pl10000 \
	#ZHToSS_SToBB_ms40_pl100000 \
	#ZHToSS_SToBB_ms55_pl100 \
	#ZHToSS_SToBB_ms55_pl1000 \
	#ZHToSS_SToBB_ms55_pl10000 \
	#ZHToSS_SToBB_ms55_pl100000 \
	##HNL_testpoint1
	#ggZHToSS_SToBB_ZToLL_ms40_pl1000 \
	#ggZHToSS_SToBB_ZToNuNu_ms40_pl1000 \
	#ggZHToSS_SToBB_ZToQQ_ms40_pl1000
	for sample in \
	ggZHToSS_SToBB_ZToQQ_ms15_pl100 \
	ggZHToSS_SToBB_ZToQQ_ms15_pl1000 \
	ggZHToSS_SToBB_ZToQQ_ms15_pl10000 \
	ggZHToSS_SToBB_ZToQQ_ms15_pl100000 \
	ggZHToSS_SToBB_ZToQQ_ms40_pl100 \
	ggZHToSS_SToBB_ZToQQ_ms40_pl1000 \
	ggZHToSS_SToBB_ZToQQ_ms40_pl10000 \
	ggZHToSS_SToBB_ZToQQ_ms40_pl100000 \
	ggZHToSS_SToBB_ZToQQ_ms55_pl100 \
	ggZHToSS_SToBB_ZToQQ_ms55_pl1000 \
	ggZHToSS_SToBB_ZToQQ_ms55_pl10000 \
	ggZHToSS_SToBB_ZToQQ_ms55_pl100000
	do
	        echo "${list_dir}${sample}.txt"
	        rm -f ${list_dir}${sample}.txt
		sample=${sample%.txt}
	        #find ${root_dir}${sample%_batch*}/*${sample##*ev150000_}*/ -name "*.root" -size +1000c >> ${list_dir}${sample}.txt
		find ${root_dir}${sample}/ -name "*.root" -size +1000c >> ${list_dir}${sample}.txt
		sed -i '/failed/d' ${list_dir}${sample}.txt
	        echo "input list created for $sample"
	done
done
