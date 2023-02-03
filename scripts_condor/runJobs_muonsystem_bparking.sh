#!/bin/bash

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

start_time=`date +%s`


./Runllp_MuonSystem_bparking $1"_"$2.txt -d=$3 -l=$4 -n=$5 -f=$1"_"$2.root
