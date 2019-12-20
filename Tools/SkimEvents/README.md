# SkimEvents
Skims events on a tree based on run and event number

This program has been comipled in c++11; I have tested this on t3-higgs using CMSSW_6_2_0.
please do cmsenv within your CMSSW_6_2_0 directory and then compile the code.

Installation:

cmsrel CMSSW_6_2_0
cd CMSSW_6_2_0/src
cmsenv
git clone https://github.com/cmorgoth/SkimEvents.git
make -j 4

To run the example you can do:

./CreateSkimmed -event_list=input/input_skim.txt --list_of_ntuples=lists/input_0.list --output_name=test.root


If everything is fine you should see the following:

[INFO]: pick events file = input/input_skim.txt
[INFO]: list of ntuple file = lists/input_0.list
[INFO]: output file name = test.root
[INFO]: TTree name = ntp1
[INFO]: run TBranch name = runNumber
[INFO]: event TBranch name = eventNumber
[INFO:] map size: 25
[INFO]: adding file #1: /mnt/tier2/store/user/cmorgoth/MuEG-Run2012A-22Jan2013-v1/default_data_100_1_srP.root
[INFO]: adding file #2: /mnt/tier2/store/user/cmorgoth/MuEG-Run2012A-22Jan2013-v1/default_data_101_1_cPo.root
[INFO]: adding file #3: /mnt/tier2/store/user/cmorgoth/MuEG-Run2012A-22Jan2013-v1/default_data_10_1_IKJ.root
[INFO]: 41043 Total Entries
[INFO]: 0 entry
******************************************************************************
*Tree    :ntp1      : ntp1                                                   *
*Entries :        0 : Total =          542767 bytes  File  Size =          0 *
*        :          : Tree compression factor =   1.00                       *
******************************************************************************

...


In the end you should get a root file named: test.root, in this casa the default tree is vecbos. In this case there are no selected events.
to change the tree, run and event number branches names used the other options in the program.

