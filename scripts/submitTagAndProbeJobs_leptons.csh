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
SingleMuon_2016B_PRv2 \
SingleMuon_2016C_PRv2 \
SingleMuon_2016D_PRv2 \
SingleMuon_2016E_PRv2 \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/lists/Run2/razorNtuplerV3p5/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/submission/MuTrigger/TagAndProbe_MuTrigger_${sample}_${jobnumber}.out -J TagAndProbe_MuTrigger_${sample}_${jobnumber} /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_Jiajing.csh RazorTagAndProbe $inputfilelist true 20560 $filesPerJob $jobnumber TagAndProbe_MuTrigger_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe_ICHEP2016/MuTrigger/jobs/


  end
end

#foreach sample( \
#DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
#)
#  setenv LSB_JOB_REPORT_MAIL N
#  set inputfilelist="/afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/lists/Run2/razorNtuplerV3p5/MC/25ns/${sample}.cern.txt"
#  set filesPerJob = 1
#  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
#  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
#  echo "Sample " $sample " maxjob = " $maxjob
#
#  foreach jobnumber(`seq 0 1 $maxjob`)
#	bsub -q 8nm -o /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/submission/TagAndProbe_MuTrigger_${sample}_${jobnumber}.out -J TagAndProbe_MuTrigger_${sample}_${jobnumber} /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_Jiajing.csh RazorTagAndProbe $inputfilelist false 20560 $filesPerJob $jobnumber TagAndProbe_MuTrigger_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe_ICHEP2016/MuTrigger/jobs/
#


#  end
#end





##########################
# Electron Trigger
##########################
foreach sample( \
SingleElectron_2016B_PRv2 \
SingleElectron_2016C_PRv2 \
SingleElectron_2016D_PRv2 \
SingleElectron_2016E_PRv2 \
DoubleEG_2016B_PRv2 \
DoubleEG_2016C_PRv2 \
DoubleEG_2016D_PRv2 \
DoubleEG_2016E_PRv2 \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/lists/Run2/razorNtuplerV3p5/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/submission/EleTrigger/TagAndProbe_EleTrigger_${sample}_${jobnumber}.out -J TagAndProbe_EleTrigger_${sample}_${jobnumber} /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_Jiajing.csh RazorTagAndProbe $inputfilelist true 10561 $filesPerJob $jobnumber TagAndProbe_EleTrigger_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe_ICHEP2016/EleTrigger/jobs/



  end
end

foreach sample( \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/lists/Run2/razorNtuplerV3p2/MC/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	bsub -q 8nm -o /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/submission/EleTrigger/TagAndProbe_EleTrigger_${sample}_${jobnumber}.out -J TagAndProbe_EleTrigger_${sample}_${jobnumber} /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_Jiajing.csh RazorTagAndProbe $inputfilelist false 10561 $filesPerJob $jobnumber TagAndProbe_EleTrigger_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe_ICHEP2016/EleTrigger/jobs/




  end
end






##########################
# Muon Selection
##########################
foreach sample( \
SingleMuon_2016B_Prv2 \
SingleMuon_2016C_Prv2 \
SingleMuon_2016D_Prv2 \
SingleMuon_2016E_Prv2 \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/lists/Run2/razorNtuplerV3p5/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	echo "job " $jobnumber " out of " $maxjob
        bsub -q 8nm -o /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/submission/TightMuonSelection/TagAndProbe_TightMuonSelectionEff_${sample}_${jobnumber}.out -J TagAndProbe_TightMuonSelectionEff_${sample}_${jobnumber} /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_Jiajing.csh RazorTagAndProbe $inputfilelist true 20105 $filesPerJob $jobnumber TagAndProbe_TightMuonSelectionEff_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe_ICHEP2016/TightMuonSelectionEff/jobs/	


        bsub -q 8nm -o /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/submission/VetoMuonSelection/TagAndProbe_VetoMuonIDEff_${sample}_${jobnumber}.out -J TagAndProbe_VetoMuonIDEff_${sample}_${jobnumber} /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_Jiajing.csh RazorTagAndProbe $inputfilelist true 20102 $filesPerJob $jobnumber TagAndProbe_VetoMuonIDEff_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe_ICHEP2016/VetoMuonIDEff/jobs/	


  end
end

foreach sample( \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/lists/Run2/razorNtuplerV3p2/MC/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
        bsub -q 8nm -o /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/submission/TightMuonSelection/TagAndProbe_TightMuonSelectionEff_${sample}_${jobnumber}.out -J TagAndProbe_TightMuonSelectionEff_${sample}_${jobnumber} /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_Jiajing.csh RazorTagAndProbe $inputfilelist false 20105 $filesPerJob $jobnumber TagAndProbe_TightMuonSelectionEff_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe_ICHEP2016/TightMuonSelectionEff/jobs/	


        bsub -q 8nm -o /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/submission/VetoMuonSelection/TagAndProbe_VetoMuonIDEff_${sample}_${jobnumber}.out -J TagAndProbe_VetoMuonIDEff_${sample}_${jobnumber} /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_Jiajing.csh RazorTagAndProbe $inputfilelist false 20102 $filesPerJob $jobnumber TagAndProbe_VetoMuonIDEff_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe_ICHEP2016/VetoMuonIDEff/jobs/	

  end
end





##########################
# Electron Selection
##########################
foreach sample( \
SingleElectron_2016B_PRv2 \
SingleElectron_2016C_PRv2 \
SingleElectron_2016D_PRv2 \
SingleElectron_2016E_PRv2 \
DoubleEG_2016B_PRv2 \
DoubleEG_2016C_PRv2 \
DoubleEG_2016D_PRv2 \
DoubleEG_2016E_PRv2 \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/lists/Run2/razorNtuplerV3p5/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	echo "job " $jobnumber " out of " $maxjob
        bsub -q 8nm -o /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/submission/TightElectronSelection/TagAndProbe_TightElectronSelection_${sample}_${jobnumber}.out -J TagAndProbe_TightElectronSelection_${sample}_${jobnumber} /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_Jiajing.csh RazorTagAndProbe $inputfilelist true 10105 $filesPerJob $jobnumber TagAndProbe_TightElectronSelection_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe_ICHEP2016/TightElectronSelectionEff/jobs/	


        bsub -q 8nm -o /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/submission/VetoElectronSelection/TagAndProbe_VetoElectronSelectionEff_${sample}_${jobnumber}.out -J TagAndProbe_VetoElectronSelectionEff_${sample}_${jobnumber} /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_Jiajing.csh RazorTagAndProbe $inputfilelist true 10102 $filesPerJob $jobnumber TagAndProbe_VetoElectronSelectionEff_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe_ICHEP2016/VetoElectronSelectionEff/jobs/	


  end
end

foreach sample( \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
#DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8\
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/lists/Run2/razorNtuplerV3p2/MC/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
        bsub -q 8nm -o /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/submission/TightElectronSelection/TagAndProbe_TightElectronSelection_${sample}_${jobnumber}.out -J TagAndProbe_TightElectronSelection_${sample}_${jobnumber} /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_Jiajing.csh RazorTagAndProbe $inputfilelist false 10105 $filesPerJob $jobnumber TagAndProbe_TightElectronSelection_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe_ICHEP2016/TightElectronSelectionEff/jobs/	


        bsub -q 8nm -o /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/submission/VetoElectronSelection/TagAndProbe_VetoElectronSelectionEff_${sample}_${jobnumber}.out -J TagAndProbe_VetoElectronSelectionEff_${sample}_${jobnumber} /afs/cern.ch/work/j/jmao/public/releases/CMSSW_7_4_7/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_Jiajing.csh RazorTagAndProbe $inputfilelist false 10102 $filesPerJob $jobnumber TagAndProbe_VetoElectronSelectionEff_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe_ICHEP2016/VetoElectronSelectionEff/jobs/	
	





  end
end












