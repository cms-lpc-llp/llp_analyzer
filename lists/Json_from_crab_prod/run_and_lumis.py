#!/bin/env python
import os
sampleName = {}
rootDir = {}
version='V1p12'
years = ['Data2016', 'Data2017']
years = ['Data2018']
sampleName['Data2018'] = [   
	'Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_Run2018A-HighMET-17Sep2018',
        'Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018',
        'Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018',
        'Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_Run2018D-HighMET-PromptReco',
]
rootDir['Data2018'] = [
#	'2018/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_v2_0763c7e2a4217afc523aa9f51cc4ae74_v1/local',
	'2018/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_v2_0763c7e2a4217afc523aa9f51cc4ae74_v3/results',
	'2018/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_v2_6b4ff9e50fb24519a20c06244287f193_v1/local',
	'2018/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_v2_af2984c56f8b0f50d251caae20b613c3_v3/results',
	'2018/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_v2_4840db2e9249b02e276eae53d568a5db_v1/local',
        ]

sampleName['Data2017'] = [
        'Run2_displacedJetMuonNtupler_V1p12_Data2017_Run2017B-HighMET-17Nov2017',
        'Run2_displacedJetMuonNtupler_V1p12_Data2017_Run2017C-HighMET-17Nov2017',
        'Run2_displacedJetMuonNtupler_V1p12_Data2017_Run2017D-HighMET-17Nov2017',
        'Run2_displacedJetMuonNtupler_V1p12_Data2017_Run2017E-HighMET-17Nov2017',
        'Run2_displacedJetMuonNtupler_V1p12_Data2017_Run2017F-HighMET-17Nov2017',
        ]
rootDir['Data2017'] = [
'2017/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2017_v2_a450f5a1182126c6e9ea15e054293037_v3/results',
'2017/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2017_v2_cc5eb93ddcf7fe360a1b23070ba25d2c_v1/local', #2017C
'2017/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2017_v2_f456c4e2cc046c5f67a90dc698f7c163_v1/local', #2017D
'2017/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2017_v2_eb3a4a526934421b9c148480a74d080a_v1/local', #2017E
'2017/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2017_v2_cc49ef85acbf8236f3123f79032f6b4e_v3/results',
]
sampleName['Data2016'] = [
        'Run2_displacedJetMuonNtupler_V1p12_Data2016_Run2016B-HighMET-07Aug17_ver1',
        'Run2_displacedJetMuonNtupler_V1p12_Data2016_Run2016B-HighMET-07Aug17_ver2',
        'Run2_displacedJetMuonNtupler_V1p12_Data2016_Run2016C-HighMET-07Aug17',
        'Run2_displacedJetMuonNtupler_V1p12_Data2016_Run2016D-HighMET-07Aug17',
        'Run2_displacedJetMuonNtupler_V1p12_Data2016_Run2016E-HighMET-07Aug17',
        'Run2_displacedJetMuonNtupler_V1p12_Data2016_Run2016F-HighMET-07Aug17',
        'Run2_displacedJetMuonNtupler_V1p12_Data2016_Run2016G-HighMET-07Aug17',
        'Run2_displacedJetMuonNtupler_V1p12_Data2016_Run2016H-HighMET-07Aug17'
        ]
rootDir['Data2016'] = [
'2016/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2016_v2_c79c957e4df865ddd65ad372670bfed2_v2/results/',#2016B-ver1
'2016/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2016_v2_06ebed3df92e4f48bbc5e368b90ae4a9_v1/local', #2016B-ver2
'2016/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2016_v2_b23ef69cd1fe9b93eee84cb6a4ca35ec_v2/results',#C
'2016/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2016_v2_4f8f8c23f2b00b48a26b22316c8e6714_v2/results',#D
'2016/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2016_v2_beaf9cf97e2b9a0f3632ba5381846b7c_v2/results',#E
'2016/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2016_v2_a9a1e67590d329c2c4df56991b033196_v1/local', #2016F
'2016/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2016_v2_47b9264200c899f36b8d0e3cf9e0880c_v2/results',#G
'2016/data/crab_prod/crab_prod_Run2_displacedJetMuonNtupler_V1p12_Data2016_v2_a0ddee24117bb83cfcc4e0a2d428b336_v2/results',#H
]
failed = 0
succeed = 0
for year in years:
	for i, sample in enumerate(sampleName[year]):
		source_file = '/afs/cern.ch/work/s/sixie/public/Production/displacedJetMuon/V1p12/'+rootDir[year][i]+'/run_and_lumis.tar.gz'
		destination_dir = 'displacedJetMuonNtuple/V1p12/'+year+'/'+sample+"/"
		os.system("mkdir -p "+destination_dir)
		os.system('scp christiw@lxplus.cern.ch:'+source_file+' '+destination_dir)
		if os.path.isfile(destination_dir+'run_and_lumis.tar.gz'):
			print sample + " succeed"
			succeed += 1
		else:
			print sample + " failed"
			failed += 1
print succeed + " datasets succeeded"
print failed + " datasets failed"

