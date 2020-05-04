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
MC_Summer16
do
	version=displacedJetMuonNtuple/V1p16/${year}/v1/sixie/
	root_dir=/mnt/hadoop/store/group/phys_exotica/delayedjets/${version}/
	list_dir=$CMSSW_BASE/src/llp_analyzer/lists/${version}
	echo $list_dir
	mkdir -p $list_dir
	
	#ggH_HToSSTobbbb_ms55_pl1000_RunIIFall18	
	#n3n2-n1-hbb-hbb_mh127_pl1000_ev100000 \
	#n3n2-n1-hbb-hbb_mh200_pl1000_ev100000 \
	#n3n2-n1-hbb-hbb_mh400_pl1000_ev100000 \
	#n3n2-n1-hinc-hgg_mh200_pl100_ev100000
	#WplusH_HToSSTobbbb_ms55_pl10000_ev150000 \
	#WminusH_HToSSTobbbb_ms55_pl10000_ev150000
	#WH_HToSSTobbbb_CSCDecayFilter_ms55_pl100000_ev150000 \
	for sample in \
	ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch1 \
	ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch3 \
	ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch4 \
	ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch1 \
	ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch3 \
	ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch4
	do
	        echo "${list_dir}${sample}.txt"
	        rm -f ${list_dir}${sample}.txt
		sample=${sample%.txt}
	        find ${root_dir}${sample%_batch*}/*${sample##*ev150000_}*/ -name "*.root" -size +1000c >> ${list_dir}${sample}.txt
	        sed -i '/failed/d' ${list_dir}${sample}.txt
	        echo "input list created for $sample"
	done
done

