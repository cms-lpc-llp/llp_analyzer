#!/bin/tcsh

foreach sample( \
#DoubleEG_2016B_03Feb2017 \
DoubleEG_2016C_03Feb2017 \
DoubleEG_2016D_03Feb2017 \
DoubleEG_2016E_03Feb2017 \
DoubleEG_2016F_03Feb2017 \
DoubleEG_2016G_03Feb2017 \
)
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/z/zhicaiz/public/release/CMSSW_9_2_5/src/RazorAnalyzer/lists/Run2/razorNtuplerV3p14/data/${sample}.cern.txt"
  set filesPerJob = 2
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  mkdir -p ../submission
  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 1nh -o /afs/cern.ch/work/z/zhicaiz/public/release/CMSSW_9_2_5/src/RazorAnalyzer/submission/ZeeTiming_${sample}_${jobnumber}.out -J ZeeTiming_${sample}_${jobnumber} runRazorJob_CERN_EOS_zhicai.csh  ZeeTiming $inputfilelist yes 10 $filesPerJob $jobnumber ZeeTiming_${sample}.Job${jobnumber}Of${maxjob}.root /eos/cms/store/group/phys_susy/razor/EcalTiming/ntuples_V3p14_04Aug2017/jobs/
    sleep 0.1
  end

end

