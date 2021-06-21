#!/bin/env python
import os
import glob
years=['Data2018']
#samples = ['Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_v2_MET_Run2018D-HighMET-PromptReco-v2_v1',
#samples =['Run2_displacedJetMuonNtupler_V1p12_Data2018_SingleMuon_17Sept2018_Run2018D-ZMu-PromptReco-v2',]
samples = [   'Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_v2_MET_Run2018A-HighMET-17Sep2018-v1_v1',
        'Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_v2_MET_Run2018B-HighMET-17Sep2018-v1_v1',
        'Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018-v1_v2_v3',
        'Run2_displacedJetMuonNtupler_V1p12_Data2018_17Sept2018_v2_MET_Run2018D-HighMET-PromptReco-v2_v1',]
ntupler_version = ['/v1/sixie/',
	'/v1/sixie/',
	'/v2/sixie/MET/',
	'/v1/sixie/',
	]

for year in years:
    version='displacedJetMuonNtuple/V1p12/'+year
    list_dir=os.getenv("CMSSW_BASE")+'/src/llp_analyzer/lists/displacedJetMuonNtuple/V1p12/'+year
    os.system('mkdir -p '+list_dir)
    for j, sample in enumerate(samples):
	inputfile = list_dir+sample+'_bad.txt'
        root_dir='/mnt/hadoop/store/group/phys_exotica/delayedjets/'+version+ntupler_version[j]+sample
        print(root_dir)
        print(inputfile)
        os.system("rm -f "+inputfile)
	import ROOT as rt
        bad_files = []
        allFiles = glob.glob(root_dir+'/*.root')
        print("processing "+str(len(allFiles))+" files")
        for i, f in enumerate(allFiles):
            if i%100 == 0:
                print("processed "+str(i)+" files")
            if os.path.getsize(f) < 1000:
		bad_files.append(f)
                continue
            curFile = rt.TFile.Open(f)
            if not curFile:
		bad_files.append(f)
                continue
        print(str(len(bad_files))+ " bad files")

	#write good files to input list
	with open(inputfile, 'w') as filehandle:
    	    for f in bad_files:
                filehandle.write('%s\n' % f)
        print( "broken file list created for "+sample)
