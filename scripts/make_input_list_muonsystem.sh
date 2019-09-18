#!/bin/bash
version=displacedJetMuonNtuple/V1p7/MC_Summer16/v11/sixie/
root_dir=/mnt/hadoop/store/group/phys_exotica/delayedjets/${version}
list_dir=$CMSSW_BASE/src/llp_analyzer/lists/${version}
echo $list_dir
mkdir -p $list_dir
for sample in \
WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 \
WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
WJetsToLNu_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
WJetsToLNu_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
WJetsToLNu_Pt-400To600_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
WJetsToLNu_Pt-600ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
WminusH_HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
WminusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
WminusH_HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
WplusH_HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
WplusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
WplusH_HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-10000_TuneCUETP8M1_13TeV-powheg-pythia8 \
WminusH_HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
WminusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
WminusH_HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
WplusH_HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
WplusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
WplusH_HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-1000_TuneCUETP8M1_13TeV-powheg-pythia8 \
WminusH_HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8 \
WminusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8 \
WminusH_HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8 \
WplusH_HToSSTobbbb_WToLNu_MH-125_MS-15_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8 \
WplusH_HToSSTobbbb_WToLNu_MH-125_MS-40_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8 \
WplusH_HToSSTobbbb_WToLNu_MH-125_MS-55_ctauS-100_TuneCUETP8M1_13TeV-powheg-pythia8
do
        echo "${list_dir}${sample}.txt"
        rm -f ${list_dir}${sample}.txt
        find ${root_dir}${sample} -name "*.root" >> ${list_dir}${sample}.txt
        sed -i '/failed/d' ${list_dir}${sample}.txt
        echo "input list created for $sample"

done

