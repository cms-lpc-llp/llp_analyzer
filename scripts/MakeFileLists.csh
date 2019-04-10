ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.5/Run1/MC/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/crab_DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6__Summer12_DR53X-PU_S10_START53_V7A-v1__06March2015__V2/150306_230801/0000/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! DYToMuMu_powheg.cern.txt

ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.5/Run1/MC/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/crab_DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6__Summer12_DR53X-PU_S10_START53_V7A-v1__06March2015__V1/150306_230826/0000/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! DYToEE_powheg.cern.txt

ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.5/Run1/MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/crab_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball__Summer12_DR53X-PU_S10_START53_V7A-v1__06March2015__V1/150306_230847/0000/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! DYJetsToLL_MG.cern.txt


ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.5/Run1/Data/DoubleMuParked/crab_DoubleMuParked_Run2012A-22Jan2013-v1__06March2015_V1/150306_231907/0000/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! Data_DoubleMuParked_Run2012A.cern.txt

ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.5/Run1/Data/DoubleMuParked/crab_DoubleMuParked_Run2012B-22Jan2013-v1__06March2015_V1/150307_035418/0000/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! Data_DoubleMuParked_Run2012B.cern.txt

ls ~/eos/cms/store/user/sixie/RazorNtupleV1.5/Run1/Data/DoubleMuParked/crab_DoubleMuParked_Run2012C-22Jan2013-v1__06March2015_V1/150306_233630/0000/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! Data_DoubleMuParked_Run2012C.cern.txt

ls ~/eos/cms/store/user/sixie/RazorNtupleV1.5/Run1/Data/DoubleMuParked/crab_DoubleMuParked_Run2012D-22Jan2013-v1__06March2015_V1/150306_231821/*/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! Data_DoubleMuParked_Run2012D.cern.txt




ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.6/Run1/Data/v1/DoubleMuParked/crab_DoubleMuParked_Run2012A-22Jan2013-v1__06March2015_V1/*/*/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! Data_DoubleMuParked_Run2012A.cern.txt
ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.6/Run1/Data/v1/DoubleMuParked/crab_DoubleMuParked_Run2012B-22Jan2013-v1__06March2015_V1/*/*/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! Data_DoubleMuParked_Run2012B.cern.txt
ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.6/Run1/Data/v1/DoubleMuParked/crab_DoubleMuParked_Run2012C-22Jan2013-v1__06March2015_V1/*/*/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! Data_DoubleMuParked_Run2012C.cern.txt
ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.6/Run1/Data/v1/DoubleMuParked/crab_DoubleMuParked_Run2012D-22Jan2013-v1__06March2015_V1/*/*/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! Data_DoubleMuParked_Run2012D.cern.txt


ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.6/Run1/Data/v1/DoubleElectron/crab_DoubleElectron_Run2012A-22Jan2013-v1__06March2015_V1/*/*/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! Data_DoubleElectron_Run2012A.cern.txt
ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.6/Run1/Data/v1/DoubleElectron/crab_DoubleElectron_Run2012B-22Jan2013-v1__06March2015_V1/*/*/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! Data_DoubleElectron_Run2012B.cern.txt
ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.6/Run1/Data/v1/DoubleElectron/crab_DoubleElectron_Run2012C-22Jan2013-v1__06March2015_V1/*/*/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! Data_DoubleElectron_Run2012C.cern.txt
ls ~/eos/cms/store/group/phys_susy/razor/run2/RazorNtupleV1.6/Run1/Data/v1/DoubleElectron/crab_DoubleElectron_Run2012D-22Jan2013-v1__06March2015_V1/*/*/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! Data_DoubleElectron_Run2012D.cern.txt




foreach f(\
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball \
DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
DiPhotonBox_Pt-10To25_8TeV-pythia6 \
DiPhotonBox_Pt-250ToInf_8TeV-pythia6 \
DiPhotonBox_Pt-250ToInf_8TeV_ext-pythia6 \
DiPhotonBox_Pt-25To250_8TeV-pythia6 \
DiPhotonBox_Pt-25To250_8TeV_ext-pythia6 \
DiPhotonJetsBox_M60_8TeV-sherpa \
DiPhotonJets_8TeV-madgraph \
DiPhotonJets_8TeV-madgraph-tarball-v2 \
DiPhotonJets_M0_8TeV-madgraph \
GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6 \
GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV_ext-pythia6 \
GJets_HT-100To200_8TeV-madgraph \
GJets_HT-200To400_8TeV-madgraph_v2 \
GJets_HT-400ToInf_8TeV-madgraph_v3 \
GJets_HT-40To100_8TeV-madgraph \
GluGluToHToGG_M-125_8TeV-powheg-pythia6 \
QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6 \
QCD_Pt_170_250_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_250_350_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_350_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6 \
TBZToLL_4F_TuneZ2star_8TeV-madgraph-tauola \
TTGJets_8TeV-madgraph \
TTH_HToGG_M-125_8TeV-pythia6 \
TTJets_FullLeptMGDecays_8TeV-madgraph-tauola \
TTJets_HadronicMGDecays_8TeV-madgraph \
TTJets_MSDecays_central_TuneZ2star_8TeV-madgraph-tauola \
TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola \
TTTT_TuneZ2star_8TeV-madgraph-tauola \
TTWJets_8TeV-madgraph \
TTWWJets_8TeV-madgraph \
TTZJets_8TeV-madgraph_v2 \
TT_CT10_TuneZ2star_8TeV-powheg-tauola \
T_s-channel_TuneZ2star_8TeV-powheg-tauola \
T_t-channel_TuneZ2star_8TeV-powheg-tauola \
T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
VBF_HToGG_M-125_8TeV-powheg-pythia6 \
WH_ZH_HToGG_M-125_8TeV-pythia6 \
WJetsToLNu_HT-150To200_8TeV-madgraph \
WJetsToLNu_HT-200To250_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph_v2 \
WJetsToLNu_HT-300To400_8TeV-madgraph \
WJetsToLNu_HT-300To400_8TeV-madgraph_v2 \
WJetsToLNu_HT-400ToInf_8TeV-madgraph \
WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2 \
WWGJets_8TeV-madgraph \
WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
WWWJets_8TeV-madgraph \
WWZNoGstarJets_8TeV-madgraph \
WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola \
WZJetsTo3LNu_8TeV_TuneZ2Star_madgraph_tauola \
WZZNoGstarJets_8TeV-madgraph \
ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp \
ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph_ext \
ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola \
ZZJetsTo2Q2Nu_TuneZ2star_8TeV-madgraph-tauloa \
ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola \
ZZZNoGstarJets_8TeV-madgraph \
)
ls ~/eos/cms/store/group/phys_susy/razor/run2/RunOneRazorNtupleV1.6/Run1/MC/v2/sixie/${f}/runOneRazorNtuplerV1p6_MC_v2_v*/*/*/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! $f.cern.txt
end



foreach f(\
DoubleElectron\
DoubleMuParked\
DoublePhoton\
HT\
HTMHTParked\
MuEG\
Photon\
SingleElectron\
SingleMu\
)
ls ~/eos/cms/store/group/phys_susy/razor/run2/RunOneRazorNtupleV1.6/Run1/Data/v4/apresyan/${f}/runOneRazorNtuplerV1p6_Data_v4_v*/*/*/*.root | sed 's/\/afs\/cern.ch\/user\/s\/sixie\/eos\/cms//' | awk '{print "root://eoscms//eos/cms/"$1}' > ! $f.cern.txt
end

 
