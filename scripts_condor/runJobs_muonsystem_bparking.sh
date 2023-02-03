#!/bin/bash

hostname
echo make an ls
ls
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=sslc7_amd64_gcc630
tar -xzf CMSSW_9_4_4.tar.gz
mkdir -p CMSSW_9_4_4/src
rm CMSSW_9_4_4.tar.gz
cd CMSSW_9_4_4/src
eval `scramv1 runtime -sh`
cp ../../Runllp_MuonSystem_bparking .
cp ../../*.root .
cp ../../lists.tgz .
tar -xzf lists.tgz

echo pwd and ls
pwd
ls



start_time=`date +%s`


./Runllp_MuonSystem_bparking lists/$1"_"$2.txt -d=$3 -l=$4 -n=$5 -f=$1"_"$2.root

echo done
ls
