#!/bin/sh
total=0
for sample in \
/MET/Run2018E-HighMET-PromptReco-v1/RAW-RECO \
/MET/Run2018D-HighMET-PromptReco-v2/RAW-RECO \
/MET/Run2018D-HighMET-PromptReco-v1/RAW-RECO \
/MET/Run2018C-HighMET-17Sep2018-v1/RAW-RECO \
/MET/Run2018B-HighMET-17Sep2018-v1/RAW-RECO \
/MET/Run2018A-HighMET-17Sep2018-v1/RAW-RECO \
/MET/Run2017F-HighMET-17Nov2017-v1/RAW-RECO \
/MET/Run2017E-HighMET-17Nov2017-v1/RAW-RECO \
/MET/Run2017D-HighMET-17Nov2017-v1/RAW-RECO \
/MET/Run2017C-HighMET-17Nov2017-v1/RAW-RECO \
/MET/Run2017B-HighMET-17Nov2017-v1/RAW-RECO \
/MET/Run2016H-HighMET-07Aug17-v1/RAW-RECO \
/MET/Run2016G-HighMET-07Aug17-v1/RAW-RECO \
/MET/Run2016F-HighMET-07Aug17-v1/RAW-RECO \
/MET/Run2016E-HighMET-07Aug17-v1/RAW-RECO \
/MET/Run2016D-HighMET-07Aug17-v1/RAW-RECO \
/MET/Run2016C-HighMET-07Aug17-v1/RAW-RECO \
/MET/Run2016B-HighMET-07Aug17_ver2-v1/RAW-RECO \
/MET/Run2016B-HighMET-07Aug17_ver1-v1/RAW-RECO
do
	tmp=$(dasgoclient -query="summary dataset=${sample} | grep summary.nevents")
	total=`python -c "print int($total+$tmp)"`
	echo ${sample} ${tmp} ${total}
done
