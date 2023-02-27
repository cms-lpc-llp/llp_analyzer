#!/bin/sh


mkdir -p log
mkdir -p submit

echo `pwd`

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/normalize.sh
listSummer16=(
ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch1
ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch3
ZH_HToSSTobbbb_ms55_pl10000_ev150000_batch4
ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch1
ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch3
ZH_HToSSTobbbb_ms55_pl1000_ev150000_batch4
WH_HToSSTobbbb_CSCDecayFilter_ms55_pl100000_ev150000
)
listSummer16=(
ZH_HToSSTobbbb_ms55_pl1000_ev150000
ZH_HToSSTobbbb_ms55_pl10000_ev150000
WH_HToSSTobbbb_CSCDecayFilter_ms55_pl100000_ev150000
)
listFall18=(
WplusH_HToSSTobbbb_ms55_pl10000_ev150000
WminusH_HToSSTobbbb_ms55_pl10000_ev150000
WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
)
listFall18=(
ttH_HToSS_SToBB_ms15_pl100
ttH_HToSS_SToBB_ms15_pl1000
ttH_HToSS_SToBB_ms15_pl10000
ttH_HToSS_SToBB_ms15_pl100000
ttH_HToSS_SToBB_ms40_pl100
ttH_HToSS_SToBB_ms40_pl1000
ttH_HToSS_SToBB_ms40_pl10000
ttH_HToSS_SToBB_ms40_pl100000
ttH_HToSS_SToBB_ms55_pl100
ttH_HToSS_SToBB_ms55_pl1000
ttH_HToSS_SToBB_ms55_pl10000
ttH_HToSS_SToBB_ms55_pl100000
WminusHToSS_SToBB_ms15_pl100
WminusHToSS_SToBB_ms15_pl1000
WminusHToSS_SToBB_ms15_pl10000
WminusHToSS_SToBB_ms15_pl100000
WminusHToSS_SToBB_ms40_pl100
WminusHToSS_SToBB_ms40_pl1000
WminusHToSS_SToBB_ms40_pl10000
WminusHToSS_SToBB_ms40_pl100000
WminusHToSS_SToBB_ms55_pl100
WminusHToSS_SToBB_ms55_pl1000
WminusHToSS_SToBB_ms55_pl10000
WminusHToSS_SToBB_ms55_pl100000
WplusHToSS_SToBB_ms15_pl100
WplusHToSS_SToBB_ms15_pl1000
WplusHToSS_SToBB_ms15_pl10000
WplusHToSS_SToBB_ms15_pl100000
WplusHToSS_SToBB_ms40_pl100
WplusHToSS_SToBB_ms40_pl1000
WplusHToSS_SToBB_ms40_pl10000
WplusHToSS_SToBB_ms40_pl100000
WplusHToSS_SToBB_ms55_pl100
WplusHToSS_SToBB_ms55_pl1000
WplusHToSS_SToBB_ms55_pl10000
WplusHToSS_SToBB_ms55_pl100000
ZHToSS_SToBB_ms15_pl100
ZHToSS_SToBB_ms15_pl1000
ZHToSS_SToBB_ms15_pl10000
ZHToSS_SToBB_ms15_pl100000
ZHToSS_SToBB_ms40_pl100
ZHToSS_SToBB_ms40_pl1000
ZHToSS_SToBB_ms40_pl10000
ZHToSS_SToBB_ms40_pl100000
ZHToSS_SToBB_ms55_pl100
ZHToSS_SToBB_ms55_pl1000
ZHToSS_SToBB_ms55_pl10000
ZHToSS_SToBB_ms55_pl100000
)
#listFall18=(
#QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8
#)
listFall17=(
ZToMuMu_NNPDF31_13TeV-powheg_M_50_120
ZToMuMu_NNPDF31_13TeV-powheg_M_120_200
ZToMuMu_NNPDF31_13TeV-powheg_M_200_400
ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300
ZToMuMu_NNPDF31_13TeV-powheg_M_2300_3500
ZToMuMu_NNPDF31_13TeV-powheg_M_3500_4500
ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000
ZToMuMu_NNPDF31_13TeV-powheg_M_6000_Inf
ZToMuMu_NNPDF31_13TeV-powheg_M_400_800
)
# listFall18=(
# HNL_testpoint1
# )
# listFall18=(
# ggZHToSS_SToBB_ZToLL_ms40_pl1000
# ggZHToSS_SToBB_ZToNuNu_ms40_pl1000
# ggZHToSS_SToBB_ZToQQ_ms40_pl1000
# )
#listFall18=(
#ggH_HToSSTo4Tau_MH-125_TuneCP5_13TeV-powheg-pythia8
#)
#listFall18=(
#QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8
#QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8
#)
#listBkg_Fall18=(
#WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
#)
#listFall18=(
#ggH_HToSS_SToEE_ms0p1_pl100
#ggH_HToSS_SToEE_ms0p1_pl500
#ggH_HToSS_SToPi0Pi0_ms1_pl100
#ggH_HToSS_SToPi0Pi0_ms1_pl500
#ggH_HToSSTobbbb_ms1_pl1000
#ggH_HToSS_SToPiPlusPiMinus_ms1_pl500
#ggH_HToSS_SToKPlusKMinus_ms1p5_pl500
#ggH_HToSS_SToEE_ms0p4_pl500
#)
listFall18=(
ggZHToSS_SToBB_ZToQQ_ms15_pl100
ggZHToSS_SToBB_ZToQQ_ms15_pl1000
ggZHToSS_SToBB_ZToQQ_ms15_pl10000
ggZHToSS_SToBB_ZToQQ_ms15_pl100000
ggZHToSS_SToBB_ZToQQ_ms40_pl100
ggZHToSS_SToBB_ZToQQ_ms40_pl1000
ggZHToSS_SToBB_ZToQQ_ms40_pl10000
ggZHToSS_SToBB_ZToQQ_ms40_pl100000
ggZHToSS_SToBB_ZToQQ_ms55_pl100
ggZHToSS_SToBB_ZToQQ_ms55_pl1000
ggZHToSS_SToBB_ZToQQ_ms55_pl10000
ggZHToSS_SToBB_ZToQQ_ms55_pl100000
ggZHToSS_SToBB_ZToLL_ms40_pl1000
ggZHToSS_SToBB_ZToNuNu_ms40_pl1000
)
#listFall18=(
#TChiHH_mass127_pl1000
#TChiHH_mass150_pl1000
#TChiHH_mass200_pl1000
#TChiHH_mass300_pl1000
#TChiHH_mass400_pl1000
#)
#listFall18_FullGenParticles=(
#ggH_HToSSTodddd_MH-125_TuneCP5_13TeV-powheg-pythia8
#ggH_HToSSTo4Tau_MH-125_TuneCP5_13TeV-powheg-pythia8
#)
#listFall18=(
#ggH_HToSS_SToEE_ms0p4_pl10
#)
#listFall18=(
#ggH_HToSS_SToEE_ms0p4_pl10
#ggH_HToSS_SToEE_ms0p4_pl50
#ggH_HToSS_SToEE_ms1p0_pl125
#ggH_HToSS_SToEE_ms1p0_pl25
#ggH_HToSS_SToEE_ms2p0_pl250
#ggH_HToSS_SToEE_ms2p0_pl50
#ggH_HToSS_SToEE_ms4p0_pl100
#ggH_HToSS_SToEE_ms4p0_pl500
#ggH_HToSS_SToGammaGamma_ms0p4_pl10
#ggH_HToSS_SToGammaGamma_ms0p4_pl50
#ggH_HToSS_SToGammaGamma_ms1p0_pl125
#ggH_HToSS_SToGammaGamma_ms1p0_pl25
#ggH_HToSS_SToGammaGamma_ms2p0_pl250
#ggH_HToSS_SToGammaGamma_ms2p0_pl50
#ggH_HToSS_SToGammaGamma_ms4p0_pl100
#ggH_HToSS_SToGammaGamma_ms4p0_pl500
#ggH_HToSS_SToGlueGlue_ms3p0_pl100
#ggH_HToSS_SToGlueGlue_ms3p0_pl500
#ggH_HToSS_SToPi0Pi0_ms0p4_pl10
#ggH_HToSS_SToPi0Pi0_ms0p4_pl50
#ggH_HToSS_SToPi0Pi0_ms1p0_pl125
#ggH_HToSS_SToPi0Pi0_ms1p0_pl25
#ggH_HToSS_SToPiPlusPiMinus_ms0p4_pl10
#ggH_HToSS_SToPiPlusPiMinus_ms0p4_pl50
#ggH_HToSS_SToPiPlusPiMinus_ms1p0_pl125
#ggH_HToSS_SToPiPlusPiMinus_ms1p0_pl25
#ggH_HToSS_SToPiPlusPiMinus_ms2p0_pl250
#ggH_HToSS_SToPiPlusPiMinus_ms2p0_pl50
#ggH_HToSS_SToPiPlusPiMinus_ms4p0_pl100
#ggH_HToSS_SToPiPlusPiMinus_ms4p0_pl500
#ggH_HToSS_STodd_ms3p0_pl100
#ggH_HToSS_STodd_ms3p0_pl500
#)
#listFall18=(
#ggH_HToSS_SToGammaGamma_ms0p4_pl50
#ggH_HToSS_SToGammaGamma_ms4p0_pl500
#ggH_HToSS_SToGlueGlue_ms3p0_pl500
#ggH_HToSS_SToPi0Pi0_ms1p0_pl125
#ggH_HToSS_SToPi0Pi0_ms1p0_pl25
#)
#listFall18=(
#ggH_HToSSTobbbb_ms1_pl1000
#ggH_HToSS_SToEE_ms0p1_pl100
#ggH_HToSS_SToEE_ms0p1_pl500
#ggH_HToSS_SToEE_ms0p4_pl500
#ggH_HToSS_SToKPlusKMinus_ms1p5_pl500
#ggH_HToSS_SToPi0Pi0_ms1_pl100
#ggH_HToSS_SToPi0Pi0_ms1_pl500
#ggH_HToSS_SToPiPlusPiMinus_ms1_pl500
#)
#
#listFall18=(
#ggH_HToSSTobbbb_ms55_pl1000
#)
#listFall18=(
#WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
#)
listFall18=(
MinBias_TuneCP5_13TeV-pythia8
)
listFall18_FullGenParticles=(
WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8
)
#
#
#listFall18=(
#ParticleGun_K0Lpt10
#ParticleGun_K0Lpt2
#ParticleGun_K0Lpt5
#ParticleGun_KPluspt10
#ParticleGun_KPluspt2
#ParticleGun_KPluspt5
#ParticleGun_Photonpt10
#ParticleGun_Photonpt2
#ParticleGun_Photonpt5
#ParticleGun_PiPluspt10
#ParticleGun_PiPluspt2
#ParticleGun_PiPluspt5
#)
listFall18=(
ParticleGun_Muonpt10
ParticleGun_Muonpt20
ParticleGun_Muonpt40
ParticleGun_Muonpt60
ParticleGun_neutronpt10
ParticleGun_neutronpt2
ParticleGun_neutronpt5
ParticleGun_uquarkpt100
ParticleGun_uquarkpt30
ParticleGun_uquarkpt300
ParticleGun_uquarkpt60
ParticleGun_K0Lpt10
ParticleGun_K0Lpt2
ParticleGun_K0Lpt5
ParticleGun_KPluspt10
ParticleGun_KPluspt5
ParticleGun_KPluspt2
ParticleGun_PiPluspt2
ParticleGun_PiPluspt5
ParticleGun_PiPluspt10
)
listFall18=(
DarkShowerHiggs_darkphoton_M15_pl1000_XIOMEGA1_XILAMBDA1
DarkShowerHiggs_darkphoton_M15_pl100_XIOMEGA1_XILAMBDA1
DarkShowerHiggs_darkphoton_M2_pl1000_XIOMEGA1_XILAMBDA1
DarkShowerHiggs_darkphoton_M2_pl100_XIOMEGA1_XILAMBDA1
)

for year in \
Fall17
do
        echo ${year}
        sampleList=list${year}[@]
        for sample in "${!sampleList}"
        do

                inputDir=/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/V1p17/MC_${year}/v1/v180/

		outputDir=${inputDir}normalized
		#inputDir=/storage/af/user/christiw/login-1/christiw/LLP/displacedJetMuonAnalyzer/csc/V1p17/MC_${year}/v1/v129/
		#outputDir=${inputDir}normalized
		echo "Sample " ${sample}
		analyzer=llp_MuonSystem_cluster
		rm -f submit/${analyzer}_normalize_${sample}*.jdl
		rm -f log/${analyzer}_normalize_${sample}*

		jdl_file=submit/${analyzer}_normalize_${sample}.jdl
		echo "Universe = vanilla" > ${jdl_file}
		echo "Executable = ${job_script}" >> ${jdl_file}
		echo "Arguments = bkg no ${sample} ${inputDir} ${outputDir} ${CMSSW_BASE} ${HOME}/" 1 >> ${jdl_file}

		echo "Log = log/${analyzer}_normalize_${sample}_PC.log" >> ${jdl_file}
		echo "Output = log/${analyzer}_normalize_${sample}_\$(Cluster).\$(Process).out" >> ${jdl_file}
		echo "Error = log/${analyzer}_normalize_${sample}_\$(Cluster).\$(Process).err" >> ${jdl_file}

		echo "RequestMemory = 16000" >> ${jdl_file}
		echo "RequestCpus = 4" >> ${jdl_file}
		echo "RequestDisk = 4" >> ${jdl_file}
                echo "+JobQueue=\"Short\"" >>  ${jdl_file}

		echo "+RunAsOwner = True" >> ${jdl_file}
		echo "+InteractiveUser = true" >> ${jdl_file}
                echo "+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7\"" >> ${jdl_file}
		echo '+SingularityBindCVMFS = True' >> ${jdl_file}
		echo "run_as_owner = True" >> ${jdl_file}
		echo "x509userproxy = ${HOME}/x509_proxy" >> ${jdl_file}
		echo "should_transfer_files = YES" >> ${jdl_file}
		echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
		echo "Queue 1" >> ${jdl_file}
		echo "condor_submit ${jdl_file}"
		condor_submit -batch-name ${sample} ${jdl_file}

	done
done
