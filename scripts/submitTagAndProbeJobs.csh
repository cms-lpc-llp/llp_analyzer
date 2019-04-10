#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================



##########################
# Muon Trigger
##########################
foreach sample( \
SingleMuon_Run2015B \
SingleMuon_Run2015C \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p17/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerIsoMu20EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerIsoMu20EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 20551 $filesPerJob $jobnumber TagAndProbe_MuTriggerIsoMu20EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerIsoMu20EffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerIsoMu27EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerIsoMu27EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 20552 $filesPerJob $jobnumber TagAndProbe_MuTriggerIsoMu27EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerIsoMu27EffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerMu50EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerMu50EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 20553 $filesPerJob $jobnumber TagAndProbe_MuTriggerMu50EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerMu50EffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerIsoMu20ORMu50EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerIsoMu20ORMu50EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 20560 $filesPerJob $jobnumber TagAndProbe_MuTriggerIsoMu20ORMu50EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerIsoMu27ORMu50EffDenominatorTight/jobs/


  end
end

foreach sample( \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p17/MC/25ns/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerIsoMu20EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerIsoMu20EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 20551 $filesPerJob $jobnumber TagAndProbe_MuTriggerIsoMu20EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerIsoMu20EffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerIsoMu27EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerIsoMu27EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 20552 $filesPerJob $jobnumber TagAndProbe_MuTriggerIsoMu27EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerIsoMu27EffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerMu50EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerMu50EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 20553 $filesPerJob $jobnumber TagAndProbe_MuTriggerMu50EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerMu50EffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerIsoMu20ORMu50EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerIsoMu20ORMu50EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 20560 $filesPerJob $jobnumber TagAndProbe_MuTriggerIsoMu20ORMu50EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerIsoMu27ORMu50EffDenominatorTight/jobs/


  end
end





##########################
# Electron Trigger
##########################
foreach sample( \
SingleElectron_Run2015B \
SingleElectron_Run2015C \
DoubleEG_Run2015B \
DoubleEG_Run2015C \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p17/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEle27LooseEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEle27LooseEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10551 $filesPerJob $jobnumber TagAndProbe_EleTriggerEle27LooseEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEle27LooseEffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEle27TightEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEle27TightEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10552 $filesPerJob $jobnumber TagAndProbe_EleTriggerEle27TightEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEle27TightEffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEle32TightEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEle32TightEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10553 $filesPerJob $jobnumber TagAndProbe_EleTriggerEle32TightEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEle32TightEffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEleCombinedEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEleCombinedEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10560 $filesPerJob $jobnumber TagAndProbe_EleTriggerEleCombinedEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEleCombinedEffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEleCombinedExtEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEleCombinedExtEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10561 $filesPerJob $jobnumber TagAndProbe_EleTriggerEleCombinedExtEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEleCombinedExtEffDenominatorTight/jobs/



  end
end

foreach sample( \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p17/MC/25ns/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEle27LooseEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEle27LooseEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10551 $filesPerJob $jobnumber TagAndProbe_EleTriggerEle27LooseEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEle27LooseEffDenominatorTight/jobs/ 
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEle27TightEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEle27TightEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10552 $filesPerJob $jobnumber TagAndProbe_EleTriggerEle27TightEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEle27TightEffDenominatorTight/jobs/ 
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEle32TightEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEle32TightEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10553 $filesPerJob $jobnumber TagAndProbe_EleTriggerEle32TightEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEle32TightEffDenominatorTight/jobs/ 
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEleCombinedEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEleCombinedEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10560 $filesPerJob $jobnumber TagAndProbe_EleTriggerEleCombinedEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEleCombinedEffDenominatorTight/jobs/
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEleCombinedExtEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEleCombinedExtEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10561 $filesPerJob $jobnumber TagAndProbe_EleTriggerEleCombinedExtEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEleCombinedExtEffDenominatorTight/jobs/




  end
end






##########################
# Muon Selection
##########################
foreach sample( \
SingleMuon_Run2015B \
SingleMuon_Run2015C \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p17/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_TightMuonSelectionEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_TightMuonSelectionEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 20105 $filesPerJob $jobnumber TagAndProbe_TightMuonSelectionEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/TightMuonSelectionEffDenominatorReco/jobs/
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_LooseMuonSelectionEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_LooseMuonSelectionEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 20103 $filesPerJob $jobnumber TagAndProbe_LooseMuonSelectionEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/LooseMuonSelectionEffDenominatorReco/jobs/


	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_VetoMuonIDEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_VetoMuonIDEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 20106 $filesPerJob $jobnumber TagAndProbe_VetoMuonIDEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/VetoMuonIDEffDenominatorReco/jobs/
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_VetoMuonSelectionEffDenominatorVetoID_${sample}_${jobnumber}.out -J TagAndProbe_VetoMuonSelectionEffDenominatorVetoID_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 20602 $filesPerJob $jobnumber TagAndProbe_VetoMuonSelectionEffDenominatorVetoID_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/VetoMuonSelectionEffDenominatorVetoID/jobs/


  end
end

foreach sample( \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p17/MC/25ns/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_TightMuonSelectionEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_TightMuonSelectionEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 20105 $filesPerJob $jobnumber TagAndProbe_TightMuonSelectionEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/TightMuonSelectionEffDenominatorReco/jobs/
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_LooseMuonSelectionEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_LooseMuonSelectionEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 20103 $filesPerJob $jobnumber TagAndProbe_LooseMuonSelectionEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/LooseMuonSelectionEffDenominatorReco/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_VetoMuonIDEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_VetoMuonIDEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 20106 $filesPerJob $jobnumber TagAndProbe_VetoMuonIDEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/VetoMuonIDEffDenominatorReco/jobs/
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_VetoMuonSelectionEffDenominatorVetoID_${sample}_${jobnumber}.out -J TagAndProbe_VetoMuonSelectionEffDenominatorVetoID_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 20602 $filesPerJob $jobnumber TagAndProbe_VetoMuonSelectionEffDenominatorVetoID_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/VetoMuonSelectionEffDenominatorVetoID/jobs/

  end
end





##########################
# Electron Selection
##########################
foreach sample( \
SingleElectron_Run2015B \
SingleElectron_Run2015C \
DoubleEG_Run2015B \
DoubleEG_Run2015C \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p17/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_TightElectronSelectionEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_TightElectronSelectionEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10105 $filesPerJob $jobnumber TagAndProbe_TightElectronSelectionEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/TightElectronSelectionEffDenominatorReco/jobs/
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_LooseElectronSelectionEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_LooseElectronSelectionEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10103 $filesPerJob $jobnumber TagAndProbe_LooseElectronSelectionEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/LooseElectronSelectionEffDenominatorReco/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_TightElectronSelectionEffDenominatorTightID_${sample}_${jobnumber}.out -J TagAndProbe_TightElectronSelectionEffDenominatorTightID_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10905 $filesPerJob $jobnumber TagAndProbe_TightElectronSelectionEffDenominatorTightID_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/TightElectronSelectionEffDenominatorTightID/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_VetoElectronIDEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_VetoElectronIDEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10106 $filesPerJob $jobnumber TagAndProbe_VetoElectronIDEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/VetoElectronIDEffDenominatorReco/jobs/
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_VetoElectronSelectionEffDenominatorVetoID_${sample}_${jobnumber}.out -J TagAndProbe_VetoElectronSelectionEffDenominatorVetoID_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10602 $filesPerJob $jobnumber TagAndProbe_VetoElectronSelectionEffDenominatorVetoID_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/VetoElectronSelectionEffDenominatorVetoID/jobs/



  end
end

foreach sample( \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8\
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p17/MC/25ns/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_TightElectronSelectionEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_TightElectronSelectionEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10105 $filesPerJob $jobnumber TagAndProbe_TightElectronSelectionEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/TightElectronSelectionEffDenominatorReco/jobs/
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_LooseElectronSelectionEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_LooseElectronSelectionEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10103 $filesPerJob $jobnumber TagAndProbe_LooseElectronSelectionEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/LooseElectronSelectionEffDenominatorReco/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_TightElectronSelectionEffDenominatorTightID_${sample}_${jobnumber}.out -J TagAndProbe_TightElectronSelectionEffDenominatorTightID_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10905 $filesPerJob $jobnumber TagAndProbe_TightElectronSelectionEffDenominatorTightID_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/TightElectronSelectionEffDenominatorTightID/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_VetoElectronIDEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_VetoElectronIDEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10106 $filesPerJob $jobnumber TagAndProbe_VetoElectronIDEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/VetoElectronIDEffDenominatorReco/jobs/
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_VetoElectronSelectionEffDenominatorVetoID_${sample}_${jobnumber}.out -J TagAndProbe_VetoElectronSelectionEffDenominatorVetoID_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10602 $filesPerJob $jobnumber TagAndProbe_VetoElectronSelectionEffDenominatorVetoID_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/VetoElectronSelectionEffDenominatorVetoID/jobs/


  end
end












##########################
# Photons
##########################
foreach sample( \
SingleElectron_Run2015C \
DoubleEG_Run2015C \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p17/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoLooseEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_PhoLooseEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 30103 $filesPerJob $jobnumber TagAndProbe_PhoLooseEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoLooseEffDenominatorReco/jobs/
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoTightEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_PhoTightEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 30105 $filesPerJob $jobnumber TagAndProbe_PhoTightEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoTightEffDenominatorReco/jobs/
  end
end

foreach sample( \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p17/MC/25ns/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoLooseEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_PhoLooseEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 30103 $filesPerJob $jobnumber TagAndProbe_PhoLooseEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoLooseEffDenominatorReco/jobs/
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoTightEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_PhoTightEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 30105 $filesPerJob $jobnumber TagAndProbe_PhoTightEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoTightEffDenominatorReco/jobs/

  end
end


