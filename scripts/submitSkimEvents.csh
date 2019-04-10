#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================


##########################
# HZZ Skim
##########################

foreach sample( \
Data_DoubleMuParked_Run2012A \
Data_DoubleMuParked_Run2012B \
Data_DoubleMuParked_Run2012C \
Data_DoubleMuParked_Run2012D \
Data_DoubleElectron_Run2012A \
Data_DoubleElectron_Run2012B \
Data_DoubleElectron_Run2012C \
Data_DoubleElectron_Run2012D \
Data_MuEG_Run2012A \
Data_MuEG_Run2012B \
Data_MuEG_Run2012C \
Data_MuEG_Run2012D \
) 
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/razorNtuplerV1p8-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/run2/HZZEventSkim/jobs/HZZEventSkim_razorNtuple_${sample}.Job${jobnumber}Of${maxjob}.root ) then
      echo "job " $jobnumber " out of " $maxjob
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/HZZEventSkim_${jobnumber}.out -J HZZEventSkim_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runSkimEvent_CERN_EOS.csh /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorHZZ/HZZLegacy/HZZEvent.list $inputfilelist RazorEvents runNum eventNum $filesPerJob $jobnumber HZZEventSkim_razorNtuple_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/run2/HZZEventSkim/jobs/
      sleep 0.1
    endif
  end

end


