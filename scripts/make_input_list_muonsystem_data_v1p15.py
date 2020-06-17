#!/bin/env python
import os
import glob
from FWCore.PythonUtilities.LumiList import LumiList as LumiList
import tarfile
sampleName = {}
rootDir = {}
goldenLumi = {}

years=['Data2016','Data2017','Data2018']
years = ['Data2017','Data2018']
goldenLumi['Data2018'] = os.getenv("CMSSW_BASE")+'/src/llp_analyzer/lists/officialJson/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
sampleName['Data2018'] = [   'Run2_displacedJetMuonNtupler_V1p15_Data2018_17Sept2018_Run2018A-HighMET-17Sep2018',
        'Run2_displacedJetMuonNtupler_V1p15_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018',
        'Run2_displacedJetMuonNtupler_V1p15_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018',
        'Run2_displacedJetMuonNtupler_V1p15_Data2018_17Sept2018_Run2018D-HighMET-PromptReco',
]
rootDir['Data2018'] = [
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2018/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2018_17Sept2018_Run2018A-HighMET-17Sep2018-v1_v1_v2/200226_145657/*/',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2018/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2018_17Sept2018_Run2018B-HighMET-17Sep2018-v1_v1_v1/200219_205929/*/',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2018/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2018_17Sept2018_Run2018C-HighMET-17Sep2018-v1_v1_v1/200219_210117/*/',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2018/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2018_17Sept2018_Run2018D-HighMET-PromptReco-v2_v1_v1/200219_210505/*/',
	]
goldenLumi['Data2017'] = os.getenv("CMSSW_BASE")+'/src/llp_analyzer/lists/officialJson/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
sampleName['Data2017'] = [
	'Run2_displacedJetMuonNtupler_V1p15_Data2017_Run2017B-HighMET-17Nov2017',
	'Run2_displacedJetMuonNtupler_V1p15_Data2017_Run2017C-HighMET-17Nov2017',
	'Run2_displacedJetMuonNtupler_V1p15_Data2017_Run2017D-HighMET-17Nov2017',
	'Run2_displacedJetMuonNtupler_V1p15_Data2017_Run2017E-HighMET-17Nov2017',
	'Run2_displacedJetMuonNtupler_V1p15_Data2017_Run2017F-HighMET-17Nov2017',	
	]
rootDir['Data2017'] = [
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2017/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2017_Run2017B-HighMET-17Nov2017-v1_v1_v1/200219_204336/0000/',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2017/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2017_Run2017C-HighMET-17Nov2017-v1_v1_v1/200219_204348/*/',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2017/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2017_Run2017D-HighMET-17Nov2017-v1_v1_v1/200219_204400/*/',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2017/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2017_Run2017E-HighMET-17Nov2017-v1_v1_v1/200219_204412/*/',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2017/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2017_Run2017F-HighMET-17Nov2017-v1_v1_v2/200225_152658/*/',
]

goldenLumi['Data2016'] = os.getenv("CMSSW_BASE")+'/src/llp_analyzer/lists/officialJson/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt'
sampleName['Data2016'] = [
	'Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016B-HighMET-07Aug17_ver1',
	'Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016B-HighMET-07Aug17_ver2',
	'Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016C-HighMET-07Aug17',
	'Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016D-HighMET-07Aug17',
	'Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016E-HighMET-07Aug17',
	'Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016F-HighMET-07Aug17',
	'Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016G-HighMET-07Aug17',
	'Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016H-HighMET-07Aug17'
	]
rootDir['Data2016'] = [
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2016/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016B-HighMET-07Aug17_ver1-v1_v1_v1/200219_203814/0000/',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2016/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016B-HighMET-07Aug17_ver2-v1_v1_v1/200219_203827/*/',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2016/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016C-HighMET-07Aug17-v1_v1_v1/200219_203839/*/',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2016/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016D-HighMET-07Aug17-v1_v1_v1/200219_203851/*',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2016/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016E-HighMET-07Aug17-v1_v1_v1/200219_203909/0000/',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2016/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016F-HighMET-07Aug17-v1_v1_v1/200219_203921/*/',
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2016/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016G-HighMET-07Aug17-v1_v1_v1/200219_203933/*/',	
'/store/group/phys_exotica/delayedjets/displacedJetMuonNtuple/V1p15/Data2016/v1/sixie/MET/Run2_displacedJetMuonNtupler_V1p15_Data2016_Run2016H-HighMET-07Aug17-v1_v1_v2/200225_152005/*/',
]
for year in years:
    version='displacedJetMuonNtuple/V1p15/'+year
    list_dir=os.getenv("CMSSW_BASE")+'/src/llp_analyzer/lists/displacedJetMuonNtuple/V1p15/'+year+'/'
    os.system('mkdir -p '+list_dir)
    assert(len(rootDir[year])==len(sampleName[year]))
    for j, sample in enumerate(sampleName[year]):
	inputfile = list_dir+sample+'.txt'
        badfilelist = list_dir+sample+'_bad.txt'
 	root_dir = '/mnt/hadoop/'+rootDir[year][j]
        print(root_dir)
        print(inputfile)
        os.system("rm -f "+inputfile)
	os.system("rm -f "+badfilelist)
	import ROOT as rt
        good_files_index = []
	good_files = []
        bad_files = []
	allFiles = glob.glob(root_dir+'/*.root')
	if len(allFiles) == 0:
	    allFiles = glob.glob(root_dir+'/*/*/*.root')
	allFiles.sort()
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
            good_files.append(f)
            index = f.split('_')[-1]
            index = int(index.split('.root',1)[0])
	    good_files_index.append(index)
	assert(len(good_files)+len(bad_files) == len(allFiles))
        print(str(len(good_files))+ " good files")
        print(str(len(bad_files))+ " bad files")

	#write good files to input list
	
	with open(inputfile, 'w') as filehandle:
    	    for f in good_files:
                filehandle.write('%s\n' % f)
        with open(badfilelist, 'w') as filehandle:
            for f in bad_files:
                filehandle.write('%s\n' % f)
	if os.path.isfile(inputfile):
	    print( "input list created for "+sample)
	if os.path.isfile(badfilelist):
	    print( "bad file list created for " + sample)

	# write good lumi json file
#	lumiTarball = os.getenv("CMSSW_BASE")+'/src/llp_analyzer/lists/Json_from_crab_prod/displacedJetMuonNtuple/V1p15/'+year+'/'+sample+'/run_and_lumis.tar.gz'
#	outputJson = os.getenv("CMSSW_BASE")+'/src/llp_analyzer/lists/Json_from_crab_prod/displacedJetMuonNtuple/V1p15/'+year+'/'+sample+'/'+sample+'_goodLumi.json'
#        os.system('rm -rf temp_dir')
#	os.mkdir('temp_dir')
#        os.chdir('temp_dir')
#        os.system("rm -f "+outputJson)
#	tf = tarfile.open(lumiTarball)
#        tf.extractall()
#        lumiFiles = glob.glob('*json')
#	finalLumi = LumiList()
#        for i in good_files_index:
#            tempJson = 'job_lumis_'+str(i)+'.json'
#            if tempJson in lumiFiles:
#                if os.path.isfile(tempJson):
#                    tempLumiList = LumiList('job_lumis_'+str(i)+'.json')
#                    finalLumi = finalLumi | tempLumiList
#                else:
#                    print('WARNING: '+ tempJson+' NOT FOUND!!!')
#        finalLumi = finalLumi & LumiList(goldenLumi[year])
#        finalLumi.writeJSON(fileName=outputJson)
#        os.chdir('../')
#        os.system('rm -r temp_dir')
#	if os.path.isfile(outputJson):
#	    print(outputJson)
#	    print 'good json file created for '+ sample


	




