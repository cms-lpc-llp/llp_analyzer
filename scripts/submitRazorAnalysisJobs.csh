#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================

##########################
# QCD Study
##########################
foreach sample( \
QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Spring15_25ns \
QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Spring15_25ns \
QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Spring15_25ns \
QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Spring15_25ns \
)
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p12/MC/25ns/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorQCDStudy_${sample}_${jobnumber}.out -J RazorQCDStudy_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorQCDStudy $inputfilelist false -1 $filesPerJob $jobnumber RazorQCDStudy_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RazorQCDStudy/Spring15/jobs/
    sleep 0.1
  end

end



##########################
# Baseline Razor Analysis
##########################
foreach sample( \
DoubleEG_Run2015B\
DoubleMuon_Run2015B\
HTMHT_Run2015B\
JetHT_Run2015B\
Jet_Run2015B\
MET_Run2015B\
MuonEG_Run2015B\
SingleElectron_Run2015B\
SingleMuon_Run2015B\
SingleMu_Run2015B\
Tau_Run2015B\
) 
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p16/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorInclusive_${jobnumber}.out -J RazorInclusive_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh razor $inputfilelist true -1 $filesPerJob $jobnumber RazorInclusive_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RazorInclusive/jobs/   
    sleep 0.1
  end

end

foreach sample( \
QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8\
) 
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p16/MC/50ns/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorInclusive_${jobnumber}.out -J RazorInclusive_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh razor $inputfilelist false -1 $filesPerJob $jobnumber RazorInclusive_${sample}_50ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RazorInclusive/jobs/  
    sleep 0.1
  end

end




foreach sample( \
SMS-T1tttt_2J_mGl-1500_mLSP-100_20bx25 \
SMS-T1bbbb_2J_mGl-1500_mLSP-100_20bx25 \
SMS-T1tttt_2J_mGl-1200_mLSP-800_20bx25 \
SMS-T1bbbb_2J_mGl-1000_mLSP-900_20bx25 \
SMS-T1qqqq_2J_mGl-1400_mLSP-100_20bx25 \
SMS-T1qqqq_2J_mGl-1000_mLSP-800_20bx25 \
SMS-T2bb_2J_mStop-600_mLSP-580_20bx25 \
SMS-T2bb_2J_mStop-900_mLSP-100_20bx25 \
SMS-T2qq_2J_mStop-1200_mLSP-100_20bx25 \
SMS-T2qq_2J_mStop-600_mLSP-550_20bx25 \
SMS-T2tt_2J_mStop-425_mLSP-325_20bx25 \
SMS-T2tt_2J_mStop-500_mLSP-325_20bx25 \
SMS-T2tt_2J_mStop-650_mLSP-325_20bx25 \
SMS-T2tt_2J_mStop-850_mLSP-100_20bx25 \
TTJets_20bx25 \
QCD_HT100To250_20bx25 \
QCD_HT250To500_20bx25 \
QCD_HT500To1000_20bx25 \
QCD_HT1000ToInf_20bx25 \
WJetsToLNu_HT100To200_20bx25 \
WJetsToLNu_HT200To400_20bx25 \
WJetsToLNu_HT400To600_20bx25 \
WJetsToLNu_HT600ToInf_20bx25 \
DYJetsToLL_HT100To200_20bx25 \
DYJetsToLL_HT200To400_20bx25 \
DYJetsToLL_HT400To600_20bx25 \
DYJetsToLL_HT600ToInf_20bx25 \
ZJetsToNuNu_HT100To200_20bx25 \
ZJetsToNuNu_HT200To400_20bx25 \
ZJetsToNuNu_HT400To600_20bx25 \
ZJetsToNuNu_HT600ToInf_20bx25 \
T_tW_20bx25 \
Tbar_tW_20bx25 \
TToLeptons_s_20bx25 \
TBarToLeptons_s_20bx25 \
TToLeptons_t_20bx25 \
TBarToLeptons_t_20bx25 \
WZJetsTo3LNu_20bx25 \
TTWJets_20bx25 \
TTZJets_20bx25 \
SMS-T1qqqq_2J_mGl-1000_mLSP-800_Tune4C_13TeV-madgraph-tauola_AVE30BX50_ST_PHYS14 \
SMS-T1bbbb_2J_mGl-1000_mLSP-900_Tune4C_13TeV-madgraph-tauola_AVE30BX50_ST_PHYS14
SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola_AVE30BX50_ST_PHYS14
SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola_AVE30BX50_ST_PHYS14
SMS-T2tt_2J_mStop-425_mLSP-325_Tune4C_13TeV-madgraph-tauola_AVE30BX50_ST_PHYS14
DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_AVE30BX50_ST_PHYS14
DYJetsToLL_M-50_13TeV-madgraph-pythia8_4bx50_PHYS14
TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_30bx50_PHYS14
) 
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/razorNtuplerV1p5-25ns_v1_v7/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorInclusive_${jobnumber}.out -J RazorInclusive_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh razor $inputfilelist false -1 $filesPerJob $jobnumber RazorInclusive_${sample}_20bx25.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RazorInclusive/jobs/
    bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/HbbRazor_${jobnumber}.out -J HbbRazor_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh hbbrazor $inputfilelist false -1 $filesPerJob $jobnumber HbbRazor_${sample}_20bx25.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/HbbRazor/jobs/
    sleep 0.1
  end

end





##########################
# Veto Lepton Study
##########################
foreach sample( \
T1bbbb_1500
T1tttt_1500
TTJets \
) 
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/lists/razorNtuplerV1p4-25ns_v1_v1/${sample}_20bx25.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/VetoLeptonStudy_${jobnumber}.out -J RazorAnalysis_VetoLeptonStudy_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh razorVetoLeptonStudy $inputfilelist 1 $filesPerJob $jobnumber VetoLeptonStudy_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/VetoLeptonStudy/jobs/
    sleep 0.1
  end

end


##########################
# Object Ntuplers
##########################
foreach sample( \
TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
ttHJetToGG_M120_13TeV_amcatnloFXFX_madspin_pythia8 \
ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8 \
ttHJetToGG_M130_13TeV_amcatnloFXFX_madspin_pythia8 \
GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8 \
QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8 \ 
)
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p15/MC/50ns/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/jobs/ElectronNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/ElectronNtupler_${jobnumber}.out -J RazorAnalysis_ElectronNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh electronNtupler $inputfilelist false 1 $filesPerJob $jobnumber ElectronNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/jobs/  
      sleep 0.1
    endif 

    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/jobs/ElectronNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/ElectronNtupler_${jobnumber}.out -J RazorAnalysis_ElectronNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh electronNtupler $inputfilelist false 0 $filesPerJob $jobnumber ElectronNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/jobs/
      sleep 0.1
    endif 

    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/jobs/ElectronNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/ElectronNtupler_${jobnumber}.out -J RazorAnalysis_ElectronNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh electronNtupler $inputfilelist false 11 $filesPerJob $jobnumber ElectronNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/jobs/
      sleep 0.1
    endif 
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/jobs/MuonNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/MuonNtupler_${jobnumber}.out -J RazorAnalysis_MuonNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh muonNtupler $inputfilelist false 1 $filesPerJob $jobnumber MuonNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/MuonNtuple/jobs/
      sleep 0.1
    endif 

    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/jobs/MuonNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/MuonNtupler_${jobnumber}.out -J RazorAnalysis_MuonNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh muonNtupler $inputfilelist false 0 $filesPerJob $jobnumber MuonNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/MuonNtuple/jobs/
      sleep 0.1
    endif 

    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/jobs/MuonNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/MuonNtupler_${jobnumber}.out -J RazorAnalysis_MuonNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh muonNtupler $inputfilelist false 11 $filesPerJob $jobnumber MuonNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/MuonNtuple/jobs/
      sleep 0.1
    endif 

    bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/JetNtupler_${jobnumber}.out -J RazorAnalysis_JetNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh jetNtupler $inputfilelist false -1 $filesPerJob $jobnumber JetNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/JetNtuple/jobs/
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TauNtupler_${jobnumber}.out -J RazorAnalysis_TauNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh tauNtupler $inputfilelist false 1 $filesPerJob $jobnumber TauNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TauNtuple/jobs/
    sleep 0.1
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TauNtupler_${jobnumber}.out -J RazorAnalysis_TauNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh tauNtupler $inputfilelist false 0 $filesPerJob $jobnumber TauNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TauNtuple/jobs/
    sleep 0.1
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TauNtupler_${jobnumber}.out -J RazorAnalysis_TauNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh tauNtupler $inputfilelist false 11 $filesPerJob $jobnumber TauNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TauNtuple/jobs/
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/jobs/PhotonNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/PhotonNtupler_${jobnumber}.out -J RazorAnalysis_PhotonNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh photonNtupler $inputfilelist false 1 $filesPerJob $jobnumber PhotonNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/jobs/
      sleep 0.1
    endif 

    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/jobs/PhotonNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/PhotonNtupler_${jobnumber}.out -J RazorAnalysis_PhotonNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh photonNtupler $inputfilelist false 0 $filesPerJob $jobnumber PhotonNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/jobs/
      sleep 0.1
    endif 

    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/jobs/PhotonNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/PhotonNtupler_${jobnumber}.out -J RazorAnalysis_PhotonNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh photonNtupler $inputfilelist false 11 $filesPerJob $jobnumber PhotonNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/jobs/
      sleep 0.1
    endif 

  end
end



##########################
# Run 2 Control Region Study
##########################
foreach sample( \
SingleMuon_Run2015B \
SingleElectron_Run2015B \
DoubleMuon_Run2015B \
DoubleEG_Run2015B \
MuonEG_Run2015B \
JetHT_Run2015B \
HTMHT_Run2015B \
Jet_Run2015B \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p16/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/jobs/RunTwoRazorControlRegions_OneLeptonFull_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim_${sample}_${jobnumber}.out -J RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorControlRegions $inputfilelist true 101 $filesPerJob $jobnumber RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/jobs/
     endif

     bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RunTwoRazorControlRegions_OneLeptonReduced_SingleLeptonSkim_${sample}_${jobnumber}.out -J RunTwoRazorControlRegions_OneLeptonReduced_SingleLeptonSkim_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorControlRegions $inputfilelist true 111 $filesPerJob $jobnumber RunTwoRazorControlRegions_OneLeptonReduced_SingleLeptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/jobs/
     sleep 0.1

     if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/jobs/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_${sample}_${jobnumber}.out -J RunTwoRazorControlRegions_DileptonFull_DileptonSkim_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorControlRegions $inputfilelist true 203 $filesPerJob $jobnumber RunTwoRazorControlRegions_DileptonFull_DileptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/jobs/
    endif
  end
end

foreach sample( \
TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8\
ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1\
ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1\
ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1\
ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1\
DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8\
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8\
DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
WWTo2L2Nu_13TeV-powheg\
WZJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8\
ZZTo2L2Nu_13TeV_powheg_pythia8\
TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8\
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p16/MC/50ns/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/jobs/RunTwoRazorControlRegions_OneLeptonFull_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RunTwoRazorControlRegions_OneLeptonFull_${sample}_${jobnumber}.out -J RunTwoRazorControlRegions_OneLeptonFull_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorControlRegions $inputfilelist false 1 $filesPerJob $jobnumber RunTwoRazorControlRegions_OneLeptonFull_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/jobs/
	sleep 0.1
    endif
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/jobs/RunTwoRazorControlRegions_OneLeptonReduced_SingleLeptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RunTwoRazorControlRegions_OneLeptonReduced_SingleLeptonSkim_${sample}_${jobnumber}.out -J RunTwoRazorControlRegions_OneLeptonReduced_SingleLeptonSkim_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorControlRegions $inputfilelist false 111 $filesPerJob $jobnumber RunTwoRazorControlRegions_OneLeptonReduced_SingleLeptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/jobs/
	sleep 0.1
    endif
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/jobs/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RunTwoRazorControlRegions_DileptonFull_DileptonSkim_${sample}_${jobnumber}.out -J RunTwoRazorControlRegions_DileptonFull_DileptonSkim_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorControlRegions $inputfilelist false 203 $filesPerJob $jobnumber RunTwoRazorControlRegions_DileptonFull_DileptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/jobs/
	sleep 0.1
    endif
  end
end


##########################
# Tag And Probe
##########################
foreach sample( \
SingleElectron_Run2015B \
DoubleEG_Run2015B \
SingleMuon_Run2015B \
DoubleMuon_Run2015B \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p16/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleSingleEleTriggerEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleSingleEleTriggerEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10550 $filesPerJob $jobnumber TagAndProbe_EleSingleEleTriggerEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleSingleEleTriggerEffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTightEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_EleTightEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10105 $filesPerJob $jobnumber TagAndProbe_EleTightEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTightEffDenominatorReco/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoLooseEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_PhoLooseEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 30103 $filesPerJob $jobnumber TagAndProbe_PhoLooseEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoLooseEffDenominatorReco/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoTightEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_PhoTightEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 30105 $filesPerJob $jobnumber TagAndProbe_PhoTightEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoTightEffDenominatorReco/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoHLTLeg1EffDenominatorLoose_${sample}_${jobnumber}.out -J TagAndProbe_PhoHLTLeg1EffDenominatorLoose_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 30350 $filesPerJob $jobnumber TagAndProbe_PhoHLTLeg1EffDenominatorLoose_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoHLTLeg1EffDenominatorLoose/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoHLTLeg2EffDenominatorLoose_${sample}_${jobnumber}.out -J TagAndProbe_PhoHLTLeg2EffDenominatorLoose_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 30351 $filesPerJob $jobnumber TagAndProbe_PhoHLTLeg2EffDenominatorLoose_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoHLTLeg2EffDenominatorLoose/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoHLTLeg1EffDenominatorMedium_${sample}_${jobnumber}.out -J TagAndProbe_PhoHLTLeg1EffDenominatorMedium_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 30450 $filesPerJob $jobnumber TagAndProbe_PhoHLTLeg1EffDenominatorMedium_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoHLTLeg1EffDenominatorMedium/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoHLTLeg2EffDenominatorMedium_${sample}_${jobnumber}.out -J TagAndProbe_PhoHLTLeg2EffDenominatorMedium_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 30451 $filesPerJob $jobnumber TagAndProbe_PhoHLTLeg2EffDenominatorMedium_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoHLTLeg2EffDenominatorMedium/jobs/

  end
end

foreach sample( \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8\
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p16/MC/50ns/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleSingleEleTriggerEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleSingleEleTriggerEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10550 $filesPerJob $jobnumber TagAndProbe_EleSingleEleTriggerEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleSingleEleTriggerEffDenominatorTight/jobs/ 

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTightEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_EleTightEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10105 $filesPerJob $jobnumber TagAndProbe_EleTightEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTightEffDenominatorReco/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoLooseEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_PhoLooseEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 30103 $filesPerJob $jobnumber TagAndProbe_PhoLooseEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoLooseEffDenominatorReco/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoTightEffDenominatorReco_${sample}_${jobnumber}.out -J TagAndProbe_PhoTightEffDenominatorReco_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 30105 $filesPerJob $jobnumber TagAndProbe_PhoTightEffDenominatorReco_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoTightEffDenominatorReco/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoHLTLeg1EffDenominatorLoose_${sample}_${jobnumber}.out -J TagAndProbe_PhoHLTLeg1EffDenominatorLoose_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 30350 $filesPerJob $jobnumber TagAndProbe_PhoHLTLeg1EffDenominatorLoose_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoHLTLeg1EffDenominatorLoose/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_PhoHLTLeg2EffDenominatorLoose_${sample}_${jobnumber}.out -J TagAndProbe_PhoHLTLeg2EffDenominatorLoose_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 30351 $filesPerJob $jobnumber TagAndProbe_PhoHLTLeg2EffDenominatorLoose_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/PhoHLTLeg2EffDenominatorLoose/jobs/


  end
end





##########################
# Hgg Razor Analysis
##########################
foreach sample( \
DoubleEG_Run2015B \
SinglePhoton_Run2015B \
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p16/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HggRazor/jobs/HggRazor_OneLeptonFull_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/HggRazor_${sample}_${jobnumber}.out -J HggRazor_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh hggrazor $inputfilelist true -1 $filesPerJob $jobnumber HggRazor_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/HggRazor/jobs/
    endif
  end
end

foreach sample( \
QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8\
QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8\
GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8\
GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8\
DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa\
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8\
TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\
QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8\
QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8\
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p16/MC/50ns/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/jobs/RunTwoRazorControlRegions_OneLeptonFull_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/HggRazor_${sample}_${jobnumber}.out -J HggRazor_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh hggrazor $inputfilelist false -1 $filesPerJob $jobnumber HggRazor_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/HggRazor/jobs/
	sleep 0.1
    endif  
  end
end



##########################
# Razor Z Analysis
##########################
foreach sample( \
#DYJetsToLL_MG \
#DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
DYToEE_powheg
DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
Data_DoubleMuParked_Run2012A \
Data_DoubleMuParked_Run2012B \
Data_DoubleMuParked_Run2012C \
Data_DoubleMuParked_Run2012D \
Data_DoubleElectron_Run2012A \
Data_DoubleElectron_Run2012B \
Data_DoubleElectron_Run2012C \
Data_DoubleElectron_Run2012D \
DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball \
TTJets_FullLeptMGDecays_8TeV-madgraph-tauola \
T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
T_s-channel_TuneZ2star_8TeV-powheg-tauola \
T_t-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
WZJetsTo3LNu_8TeV_TuneZ2Star_madgraph_tauola \
ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola \
TTWJets_8TeV-madgraph \
TTWWJets_8TeV-madgraph \
TTZJets_8TeV-madgraph_v2 \
TTTT_TuneZ2star_8TeV-madgraph-tauola \
) 
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p6-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorZAnalysis/jobs/RazorZAnalysis_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorZAnalysis_${sample}_${jobnumber}.out -J RazorAnalysis_RazorZAnalysis_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh razorZAnalysis $inputfilelist 1 $filesPerJob $jobnumber RazorZAnalysis_${sample}.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorZAnalysis/jobs/
	sleep 0.1
    endif
  end

end




  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorZAnalysis/jobs/RazorZAnalysis_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh RazorZAnalysis $inputfilelist 1 $filesPerJob $jobnumber RazorZAnalysis_${sample}.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorZAnalysis/jobs/
	sleep 0.1
    endif
  end

 foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorZAnalysis/jobs/RazorZAnalysis_DielectronSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh RazorZAnalysis $inputfilelist 2 $filesPerJob $jobnumber RazorZAnalysis_DielectronSkim_${sample}.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorZAnalysis/jobs/
	sleep 0.1
    endif
  end





##########################
# HZZ Razor Analyzer
##########################
foreach sample( \
ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball \
TTJets_FullLeptMGDecays_8TeV-madgraph-tauola \
T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
T_s-channel_TuneZ2star_8TeV-powheg-tauola \
T_t-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola\
WZJetsTo3LNu_8TeV_TuneZ2Star_madgraph_tauola \
ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola \
TTWJets_8TeV-madgraph \
TTWWJets_8TeV-madgraph \
TTZJets_8TeV-madgraph_v2 \
TTTT_TuneZ2star_8TeV-madgraph-tauola \
) 
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/razorNtuplerV1p8-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HZZRazor/jobs/HZZRazor_${sample}.Job${jobnumber}Of${maxjob}.root ) then
      echo "job " $jobnumber " out of " $maxjob
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/HZZRazor_${jobnumber}.out -J HZZRazor_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh hzzRunOneRazor $inputfilelist false -1 $filesPerJob $jobnumber HZZRazor_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/HZZRazor/jobs/    
    endif 
    sleep 0.1
  end

end

foreach sample( \
GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6 \
TTbarH_HToZZTo4L_M-125_8TeV-pythia6 \
VBF_HToZZTo4L_M-125_8TeV-powheg-pythia6 \
WH_HToZZTo4L_M-125_8TeV-pythia6 \
ZH_HToZZTo4L_M-125_8TeV-pythia6 \
GluGluToZZTo4L_8TeV-gg2zz-pythia6 \
GluGluToZZTo2L2L_TuneZ2star_8TeV-gg2zz-pythia6 \
ZZTo2e2mu_8TeV-powheg-pythia6 \
ZZTo2e2mu_8TeV_ext-powheg-pythia6 \
ZZTo2e2tau_8TeV-powheg-pythia6 \
ZZTo2e2tau_8TeV_ext-powheg-pythia6 \
ZZTo2mu2tau_8TeV-powheg-pythia6 \
ZZTo2mu2tau_8TeV_ext-powheg-pythia6 \
ZZTo4e_8TeV-powheg-pythia6 \
ZZTo4e_8TeV_ext-powheg-pythia6 \
ZZTo4mu_8TeV-powheg-pythia6 \
ZZTo4mu_8TeV_ext-powheg-pythia6 \
ZZTo4tau_8TeV-powheg-pythia6 \
ZZTo4tau_8TeV_ext-powheg-pythia6 \
GluGluToHToZZTo4L_M-125_8TeV-powheg15-pythia6 \
ZZTo4e_8TeV_mll8_mZZ95-160-powheg15-pythia6 \
ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6 \
ZZTo4tau_8TeV_mll8_mZZ95-160-powheg15-pythia6 \
GluGluToHToZZTo4L_M-125_8TeV-minloHJJ-pythia6-tauola \
) 
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/razorNtuplerV1p9-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HZZRazor/jobs/HZZRazor_${sample}.Job${jobnumber}Of${maxjob}.root ) then
      echo "job " $jobnumber " out of " $maxjob
      bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/HZZRazor_${jobnumber}.out -J HZZRazor_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh hzzRunOneRazor $inputfilelist false -1 $filesPerJob $jobnumber HZZRazor_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/HZZRazor/jobs/    
    endif
    sleep 0.1
  end

end


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
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HZZRazor/jobs/HZZRazor_${sample}.Job${jobnumber}Of${maxjob}.root ) then
      echo "job " $jobnumber " out of " $maxjob
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/HZZRazor_${jobnumber}.out -J HZZRazor_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh hzzRunOneRazor $inputfilelist true -1 $filesPerJob $jobnumber HZZRazor_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/HZZRazor/jobs/    
      sleep 0.1
    endif
  end

end






##########################
# MC for ZInv
##########################
foreach sample( \
ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext \
WJetsToLNu_HT-150To200_8TeV-madgraph \
WJetsToLNu_HT-200To250_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph_v2 \
WJetsToLNu_HT-300To400_8TeV-madgraph \
WJetsToLNu_HT-300To400_8TeV-madgraph_v2 \
WJetsToLNu_HT-400ToInf_8TeV-madgraph \
WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2 \
) 
  set inputfilelist="/afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p6-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 8nm -o /afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/RazorControlRegions_${jobnumber}.out -J RazorAnalysis_RazorPhoton_${jobnumber} /afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh razorPhotonStudy $inputfilelist 0 $filesPerJob $jobnumber RazorControlRegions_${sample}.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/jobs/
    sleep 0.1
  end

end

# DATA for ZInv
foreach sample( \
DoubleMuParked \
SingleMu \
Photon \
Data_SinglePhoton_Run2012C \
Data_SinglePhoton_Run2012B \
Data_SinglePhotonParked_Run2012D \
) 
  set inputfilelist="/afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p6-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 8nm -o /afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/RazorControlRegions_${jobnumber}.out -J RazorAnalysis_RazorPhoton_${jobnumber} /afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh razorPhotonStudy $inputfilelist 1 $filesPerJob $jobnumber RazorControlRegions_${sample}.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/jobs/
    sleep 0.1
  end

end


## CRAB submit    
foreach sample( \
ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext \
WJetsToLNu_HT-150To200_8TeV-madgraph \
WJetsToLNu_HT-200To250_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph_v2 \
WJetsToLNu_HT-300To400_8TeV-madgraph \
WJetsToLNu_HT-300To400_8TeV-madgraph_v2 \
WJetsToLNu_HT-400ToInf_8TeV-madgraph \
WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2 \
QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_170_250_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_250_350_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_350_EMEnriched_TuneZ2star_8TeV_pythia6 \
GJets_HT-40To100_8TeV-madgraph \
GJets_HT-100To200_8TeV-madgraph \
GJets_HT-200To400_8TeV-madgraph_v2 \
GJets_HT-400ToInf_8TeV-madgraph_v3 \
TTGJets_8TeV-madgraph \	
DoubleMuParked \
SingleMu \
Photon \
Data_SinglePhoton_Run2012B \
Data_SinglePhoton_Run2012C \
Data_SinglePhotonParked_Run2012D \
)
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p6-Run1/${sample}.cern.txt"
  set njobs = `cat $inputfilelist | wc | awk '{print $1}' `
  set anatype = razorPhotonStudy
  set option = 0

  sed "s/sampleName/$sample/" crab_runRazorRun.py > crab_tmp.py
  sed -i "s/runRazorCrab/tmp_runRazorCrab/" crab_tmp.py
  sed -i "s/ntupleName.root/RazorAnalysis_$sample.root/" crab_tmp.py 
  sed -i "s|listfile.txt|$inputfilelist|" crab_tmp.py 
  sed -i "s/999/$njobs/" crab_tmp.py
  sed -i "s/999/$njobs/" crab_tmp.py

  sed "s/listfile.txt/${sample}.cern.txt/" runRazorCrab.sh > tmp_runRazorCrab.sh
  sed -i "s/ntupleName.root/RazorAnalysis_$sample.root/" tmp_runRazorCrab.sh
  sed -i "s/Analyses/$anatype/" tmp_runRazorCrab.sh
  sed -i "s/Option/$option/" tmp_runRazorCrab.sh

  crab submit -c crab_tmp.py
  rm crab_tmp.py tmp_runRazorCrab.sh

end

    

foreach sample ( \
Data_SinglePhotonParked_Run2012D \
)
    set inputfilelist = tmp/cms/store/group/phys_susy/razor/Run2Analysis/RazorNtupleV1.6/Run1/Test/MinBias/crab_$sample/
    set stringList = ""
    foreach dir(`find $inputfilelist -type f -name '*root*' | sed -r 's|/[^/]+$||' | sort | uniq`)
    set root = "/root"
    set stringList = "$stringList $dir$root"
    end

echo $stringList | sed  "s/\/root/\/*root/g" > tmp.txt 
set files = `cat tmp.txt`
hadd $sample.root $files && rm tmp.txt
end

##########################
# HggRazor Bkg shape study - Inverse Photon Isolation - 2015 data
##########################
foreach sample( \
DoubleEG_2015C \
DoubleEG_2015D \
)
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/z/zhicaiz/public/RazorAnalyzer/CMSSW_7_6_3/src/RazorAnalyzer/lists/Run2/razorNtuplerV2p4/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 8nm -o /afs/cern.ch/work/z/zhicaiz/public/RazorAnalyzer/CMSSW_7_6_3/src/RazorAnalyzer/submission/HggBkgShapeStudy_${sample}_${jobnumber}.out -J HggBkgShapeStudy_${sample}_${jobnumber} /afs/cern.ch/work/z/zhicaiz/public/RazorAnalyzer/CMSSW_7_6_3/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_zhicai.csh  HggRazorForBkgShape $inputfilelist yes 10 $filesPerJob $jobnumber HggRazorForBkgShape_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/HggRazor/ICHEP2016Combined/V3p4_BkgShapeStudyOnly_PhotonCorrJuly20_InvertedIso_20160810/2015jobs/
    sleep 0.1
  end

end

##########################
# HggRazor Bkg shape study - Inverse Photon Isolation - 2016 data
##########################

foreach sample( \
DoubleEG_2016B_PRv2 \
DoubleEG_2016C_PRv2 \
DoubleEG_2016D_PRv2 \
)
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/z/zhicaiz/public/RazorAnalyzer/CMSSW_7_6_3/src/RazorAnalyzer/lists/Run2/razorNtuplerV3p4/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 8nm -o /afs/cern.ch/work/z/zhicaiz/public/RazorAnalyzer/CMSSW_7_6_3/src/RazorAnalyzer/submission/HggBkgShapeStudy_${sample}_${jobnumber}.out -J HggBkgShapeStudy_${sample}_${jobnumber} /afs/cern.ch/work/z/zhicaiz/public/RazorAnalyzer/CMSSW_7_6_3/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS_zhicai.csh  HggRazorForBkgShape $inputfilelist yes 20 $filesPerJob $jobnumber HggRazorForBkgShape_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/HggRazor/ICHEP2016Combined/V3p4_BkgShapeStudyOnly_PhotonCorrJuly20_InvertedIso_20160810/2016jobs/
    sleep 0.1
  end

end

