#!/bin/bash
#Data2016/    Data2017/    Data2018/    MC_Autumn18/ MC_Fall17/   MC_Summer16/
version=displacedJetMuonNtuple/V1p17//Data2016/v4/sixie/MET/
root_dir=/mnt/hadoop/store/group/phys_exotica/delayedjets/${version}
list_dir=$CMSSW_BASE/src/cms_lpc_llp/llp_analyzer/lists/${version}
echo $list_dir
mkdir -p $list_dir
for sample in \
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-HighMET-07Aug17_ver1-v1_v4_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016B-HighMET-07Aug17_ver2-v1_v4_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016C-HighMET-07Aug17-v1_v4_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016D-HighMET-07Aug17-v1_v4_v2 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016E-HighMET-07Aug17-v1_v4_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016F-HighMET-07Aug17-v1_v4_v2 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016G-HighMET-07Aug17-v1_v4_v1 \
Run2_displacedJetMuonNtupler_V1p17_Data2016_Run2016H-HighMET-07Aug17-v1_v4_v1 \
do
        echo "${list_dir}${sample}.txt"
        rm -f ${list_dir}${sample}.txt
        find ${root_dir}${sample} -name "*.root" >> ${list_dir}${sample}.txt
        sed -i '/failed/d' ${list_dir}${sample}.txt
        echo "input list created for $sample"

done

