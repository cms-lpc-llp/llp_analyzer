VERSION = "V3p2"
TREETYPES = { '1L':'OneLeptonRazorSkimLeptonSkim',
              '1LInv':'OneLeptonInvRazorSkimLeptonSkim',
              '2L':'DileptonRazorSkimDileptonSkim',
              '2LInv':'DileptonInvRazorSkimDileptonSkim',
              'VetoL':'VetoLeptonRazorSkim',
              'VetoTau':'VetoTauRazorSkim',
              'Photon':'PhotonRazorSkimPhotonSkim',
              'Signal':'FullRazorInclusive',
              }
SKIMS = { '1L':'SingleLeptonSkim',
          '1LInv':'SingleLeptonSkim',
          '2L':'DileptonSkim',
          '2LInv':'DileptonSkim',
          'VetoL':'',
          'VetoTau':'',
          'Photon':'',
          'Signal':'',
          }
DIR = 'eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/'+VERSION
DIRS = { tag:DIR+'/'+TREETYPES[tag] for tag in TREETYPES }
DIRS['Signal'] = 'eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/'+VERSION

OPTIONS = { '1L':1101,
            '1LInv':2102,
            '2L':1203,
            '2LInv':3204,
            'VetoL':1007,
            'VetoTau':1009,
            'Photon':5505,
            'Signal':10
          }

SAMPLES = {}
SAMPLES['1L'] = {
        "DYJets":[
                'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'DYJetsToLL_M-5to50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'DYJetsToLL_M-5to50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'DYJetsToLL_M-5to50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'DYJetsToLL_M-5to50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            ],
        "QCD":[
                'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            ],
        "SingleTop":[
                'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1',
                'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1',
                'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
                'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
            ],
        "TTJets":[
                'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            ],
        "Other":[
                'WWTo2L2Nu_13TeV-powheg',
                'WWTo4Q_13TeV-powheg',
                'WWToLNuQQ_13TeV-powheg',
                'WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8',
                'WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8',
                'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',
                'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8',
                'ZZTo2L2Nu_13TeV_powheg_pythia8',
                'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',
                'ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8',
                'ZZTo4L_13TeV_powheg_pythia8',
                'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8',
                'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8',
                'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8',
                'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
                'TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
                'WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
                'WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
                'WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
                'ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
                'TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
            ],
        "WJets":[
                'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            ],
        "ZInv":[
                'ZJetsToNuNu_HT-100To200_13TeV-madgraph',
                'ZJetsToNuNu_HT-200To400_13TeV-madgraph',
                'ZJetsToNuNu_HT-400To600_13TeV-madgraph',
                'ZJetsToNuNu_HT-600To800_13TeV-madgraph',
                'ZJetsToNuNu_HT-800To1200_13TeV-madgraph',
                'ZJetsToNuNu_HT-1200To2500_13TeV-madgraph',
                'ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph',
            ],
        }
SAMPLES['1LInv'] = SAMPLES['1L'].copy()
SAMPLES['2L'] = SAMPLES['1L'].copy()
SAMPLES['2LInv'] = SAMPLES['1L'].copy()
SAMPLES['VetoL'] = SAMPLES['1L'].copy()
SAMPLES['VetoTau'] = SAMPLES['1L'].copy()
SAMPLES['Photon'] = {
        "GJetsDR0p4":[
                'GJets_DR-0p4_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'GJets_DR-0p4_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'GJets_DR-0p4_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'GJets_DR-0p4_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'GJets_DR-0p4_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            ],
        "GJets":[
                'GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            ],
        "QCD":[
                'QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8',
                'QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8',
                'QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8',
                'QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8',
            ],
        "Other":[
                'TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8',
                'TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8',
                'WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
                'WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
                'ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
            ],
        }
SAMPLES['Signal'] = SAMPLES['1L'].copy()

DATA = {}
DATA['1L'] = {
        "SingleMuon":[
                'SingleMuon_2016B_PRv2',
            ],
        "SingleElectron":[
                'SingleElectron_2016B_PRv2',
            ]
        }
DATA['1LInv'] = DATA['1L'].copy()
DATA['2L'] = DATA['1L'].copy()
DATA['2LInv'] = DATA['1L'].copy()
DATA['VetoL'] = { 'HTMHT':['HTMHT_2016B_PRv2'] }
DATA['VetoTau'] = { 'HTMHT':['HTMHT_2016B_PRv2'] }
DATA['Photon'] = {
        "SinglePhoton":[
                'SinglePhoton_2016B_PRv2',
                ]
        }
DATA['Signal'] = DATA['1L'].copy()
DATA['Signal']["HTMHT"] = ['HTMHT_2016B_PRv2']
