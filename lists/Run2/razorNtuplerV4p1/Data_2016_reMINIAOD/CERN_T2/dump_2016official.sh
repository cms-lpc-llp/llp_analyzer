#!/bin/bash

baseDir=/eos/cms/store/group/phys_susy/razor/run2/Run2RazorNtupleV4.1/Data/DoubleEG2016ReMiniAOD_EcalRechits/v1/zhicaiz/DoubleEG/

for run in 2016B_ver1 2016B_ver2 2016C 2016D 2016E 2016F 2016G 2016H
do
	echo 'xxxxxxxxxx' > DoubleEG_${run}_06Aug2018.cern.txt
	subBase1=${baseDir}crab_prod_Run2RazorNtuplerV4p1_DoubleEG_${run}_ReMiniAOD2016Official_v1_06Aug2018/
	echo ${subBase1}
	for sub1 in `ls ${subBase1}`
	do
		subBase2=${subBase1}${sub1}/
		echo ${subBase2}
		for sub2 in `ls ${subBase2}`
		do
			subBase3=${subBase2}${sub2}/
			echo ${subBase3}
			eos find ${subBase3} | grep root >> DoubleEG_${run}_06Aug2018.cern.txt
		done
	done
	sed -i '/xxxxxxxxxx/d' DoubleEG_${run}_06Aug2018.cern.txt
done
sed -i "s/\/eos\/cms/root\:\/\/cms-xrd-global.cern.ch\//g" `grep \/eos\/cms -rl DoubleEG_2016*06Aug2018.cern.txt`
