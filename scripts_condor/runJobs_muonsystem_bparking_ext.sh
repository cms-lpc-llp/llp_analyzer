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
cp ../../Runllp_MuonSystem_bparking_ext .
cp ../../*.root .
cp ../../lists.tgz .
cp ../../JEC.tar.gz .
tar -xzf lists.tgz
tar -xzf JEC.tar.gz

echo " "
echo pwd and ls
echo " "
pwd
ls

start_time=`date +%s`

file=$1"_"$2.root
eospath=$6
./Runllp_MuonSystem_bparking_ext lists/$1"_"$2.txt -d=$3 -l=$4 -n=$5 -f=$file
if $7
then
  echo the eospath
  echo $6, $eospath
  xrdcp -f $file root://cmseos.fnal.gov//$eospath/ 2>&1
else
mv $file ../../
fi
echo done
echo "1", $1 
echo "2", $2 
echo "3", $3 
echo "4", $4 
echo "5", $5 
echo "6", $6 
echo "7", $7 
ls
