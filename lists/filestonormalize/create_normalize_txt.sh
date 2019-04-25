#!/bin/bash

function dirnametofilename() {
  for f in $*; do
    bn=$(basename "$f")
    ext="${bn##*.}"
    filename=$(basename "$f")
    filepath=$(dirname "$f")
    dirname=$(basename "$filepath")
    dataset=${filename/_mx*/}
    dataset=${dataset/ppTohToSS1SS2_SS1Tobb_SS2Tobb_/}
    dataset=${dataset/ppTohToSS1SS2_SS1Tobb_SS2Toveve_/}
    echo "$dirname$dataset $f"
  done
}
dir=/mnt/hadoop/store/group/phys_exotica/delayedjets/llp_analyzer/V1p0/MC_Summer16/v1/signals/wH
export -f dirnametofilename

find $dir -name "*vh*.root" -exec bash -c 'dirnametofilename "{}"'  \;

