RazorAnalyzer
=============

Class for analyzing the 2015 razor ntuples

Setup
-------------

    cmsrel CMSSW_9_4_4
    cd CMSSW_9_4_4/src
    git clone git@github.com:cms-lpc-llp/llp_analyzer.git
    cd llp_analyzer
    make
  
Defining a new analysis
-------------
1) Copy analyzers/DummyAnalyzer.cc and replace each instance of "DummyAnalyzer" with the name of your desired analyzer.
   Modify the body of the Analyze function to define your analyzer's behavior.
   DO NOT need to write a header file for the analyzer class; the Makefile will generate one for you automatically.  

2) Do `make`.  This will create an executable `bin/Run<name of your analyzer>`. You can execute your analysis using this program directly or by calling it via the `RazorRun` script. 

Running
------------
After compiling, 

    ./RazorRun_T2 <list of input files> <name of your analyzer> <options>
  
Example: to execute a dummy analysis that does nothing,

    ./RazorRun_T2 lists/TTJets_List_Test.txt DummyAnalyzer

The "options" are the following:
    
    -d   --isData
    -f=  --outputFile=<output filename> (optional)
    -n=  --optionNumber=<option number> (optional)
    -l=  --optionLabel=<option Label> (optional)
    -h   --help


## Run the llp_analyzer
    ./RazorRun_T2 <list of input files> llp_vH -d=${isData} -n=${option} -f=${outputfile} -l=${tag}
* ```isData``` is ```yes``` or ```no```
* ```option``` is 1 if run on condor, else ```option``` can be any other number, if run interactively. (It just sets the directory of JEC parameters differently)
* ```tag``` can be ```Razor2017_17Nov2017Rereco_EOY_RERECO```, ```Razor2016_07Aug2017Rereco```, ```Razor2018_17SeptEarlyReReco``` (different tag for each year, defined in ```src/RazorHelper.cc```)
* list of input files are stored in ```lists/llpntuple/V1p0/MC_Summer16/v1/``` , same format as the llpntuple storage space in ```/mnt/hadoop/```




### Create input list
##### Monte Carlo
Run  ```scripts/make_input_list_muonsystem_sig.sh```
 * Check the ```root_dir```, which is where the ntuples are
 * Check ```list_dir``` where you want the input list to be
##### Data
First run ```lists/Json_from_crab_prod/run_and_lumis.py``` to copy all the processed json files from lxplus to t2
Then run ```scripts/make_input_list_muonsystem_data.py```
* The script checks if there are broken files in the ntuple directory and outputs a good file list, bad file list, and a lumi file that contains: good file lumi AND good lumi.
* ```rootDir``` dictionary stores the location of the ntuples for each run period
* ``` sampleName``` stores the sample name of the ntuples for each run period
* ```goldenLumi``` stores the good lumi list for each year

### Submit condor jobs on tier2
Before submitting jobs, make sure proxy and CMSSW environment is setup.

* run the ```llp_MuonSystem``` analyzer for signal, bkg or data:
	* ```scripts_condor/submit_llp_wH_muonsystem_bkg.sh```
	* ```scripts_condor/submit_llp_wH_muonsystem_sig.sh```
	* ```scripts_condor/submit_llp_wH_muonsystem_data*.sh```
	* Check ```inputfilelist``` where the list you created in step 1 is stored, and ```output``` where you want the output do be stored
* Normalize the ntuples
  * ```scripts_condor/submit_normalize_muonsystem_*.sh```
  * Check ```outputDir``` and ```inputDir```
  

Normalizing the processed ntuples
------------
The NormalizeNtuple macro opens a specified set of files and adds a 'weight' branch to each TTree in each file.  The value of 'weight' is the same for all events in a tree and is equal to lumi * CrossSection/NEvents, where NEvents is the total number of events processed for the given dataset, and lumi is the luminosity normalized to.  The cross sections can be found in the file ```data/xSections.dat```.  To run NormalizeNtuple:

    ./NormalizeNtuple <input file list> [lumi]

See lists/filestonormalize/testTTJets.txt for an example input file to be used with NormalizeNtuple.

* Make sure the dataset being processed have xSections in ```data/xSections.dat```
* Create input file list using ```scripts/create_normalize_txt.py```

The script ```hadd_llp_bkg.sh``` automatically hadd and normalize the llp_analyzer ROOT files for the background samples.



### Hadd ntuples
* hadd signals (hadd WminusH and WplusH)
	*  ```scripts/hadd_displaced.sh```
* hadd different bins of QCD/ZJetstoNunu:
	* ```hadd_qcd.sh``` or ```hadd_ZJetToNuNu.sh```
* hadd different run periods in a particular year for Data:
	* ```hadd_data*.sh```
### Filter good lumi events for data

https://github.com/RazorCMS/RazorCommon


Fitting samples and setting limits
-----------

The main controller of the fit and limit setting comes from the
configuration file, which defines the binning used in the binned fit,
the initlial values of the shape parameters, etc.

    config/run2_sideband.config

Setup combine from lxplus:

	mkdir ~/work/RAZORRUN2/
	cd ~/work/RAZORRUN2/
    cmsrel CMSSW_7_1_5
    cd CMSSW_7_1_5/src/
    git clone https://github.com/RazorCMS/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cd HiggsAnalysis/CombinedLimit
    git pull origin razor1dpdf_71X
    scramv1 b clean; scramv1 b
    cd ../..

Proceed with the usual setup of RazorAnalyzer:

    git clone https://github.com/RazorCMS/RazorAnalyzer.git
    cd RazorAnalyzer
    make
  
Make output directories to make workflow easier 

    mkdir Datasets
    mkdir FitResults; mkdir FitProjections; mkdir cards

Now convert the SM MC ntuples into a SM Cocktail RooDataSet (ignoring
the QCD contirbution).

     python python/DustinTuple2RooDataSet.py -c config/run2_sideband.config -b MultiJet -d Datasets/ -w -l 4000 \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_DYJetsToLL_M-50_HTBinned_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_DYJetsToLL_M-5to50_HTBinned_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_SingleTop_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_TTV_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_VV_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_WJetsToLNu_HTBinned_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_ZJetsToNuNu_HTBinned_1pb_weighted.root 

Similarly you can run over the data,

     python python/DustinTuple2RooDataSet.py -b MultiJet -c config/run2_sideband.config -d Datasets/ --data -l 1264 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/Data/RazorInclusive_HTMHT_Run2015D_Oct05ReMiniAOD_PRv4_GoodLumiGolden.root

To perform the fit on the SM Cocktail,

     python python/BinnedFit.py -b MultiJet -c config/run2_sideband.config -d FitResults -l 4000 Datasets/RazorInclusive_SMCocktail_weighted_lumi-4.000_0-3btag_MultiJet.root

To produce the signal templates,

	 python python/SMSTemplates.py -c config/run2_sideband.config -b MultiJet -d Datasets/ -l 4000 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p22_ForPreappFreezing20151106/jobs/combined/SMS-T1bbbb_1500_100.root 

To run Bayesian toys, marginalizing the shape parameters to
estimate the systematic uncertainty in each bin, and produce a large number of other 1D fit projections, and the 2D
projections in MR and Rsq,

     python python/RunToys.py -b MultiJet -c config/run2_sideband.config -i FitResults/BinnedFitResults_MultiJet.root -d FitResults  -t 10000
     python python/PlotFit.py -b MultiJet -c config/run2_sideband.config -d FitResults -i FitResults/BinnedFitResults_MultiJet.root -t FitResults/toys_Bayes_MultiJet.root
	
Next, to produce the datacards and run combine, execute the following:

     python python/DustinTuple2RooDataSet.py -b MultiJet -c config/run2_sideband.config -w -l 4000 -d Datasets Signals/RazorInclusive_SMS-T1bbbb_2J_mGl-1500_mLSP-100_1pb_weighted.root
     python python/RunCombine.py -b MultiJet -c config/run2_sideband.config -d cards --lumi-array 4 -m T1bbbb --mGluino 1500 --mLSP 100

To submit jobs to run the limits on the data:

	 mkdir ~/work/RAZORRUN2/Limits/
     python python/RunCombineJobs.py -i FitResults/BinnedFitResults.root -c config/run2_sideband.config --data -l 1.264 -m T1bbbb -b MultiJet -d cards/ -q 8nh

After the jobs finish, the output files will be in Limits directory we
just created

	 cp ~/work/RAZORRUN2/Limits/cards/higgs*.root cards/
     python python/GetCombine.py -m T1bbbb -b MultiJet -c config/run2_sideband.config -d cards/ -l 1.264
     python python/Get2DContour.py -m T1bbbb -b MultiJet  -d cards/ 

To make the final (expected) limit plot, we need to check out a different repository. 

     git clone git@github.com:RazorCMS/PlotsSMS
     cd PlotsSMS

Note you have to change the smoothed cross section limit file location
in config/SUS15004/T1bbbb\_Exp\_SUS15004.cfg. Then you can run it,

	 python python/makeSMSplots.py config/SUS15004/T1bbbb_Exp_SUS15004.cfg T1bbbbAsymptotic

To make an "unweighted" dataset (not needed for the preceding commands),
execute

     python python/RooDataSet2UnweightedDataSet.py -b MultiJet -c config/run2_sideband.config -d Datasets Datasets/RazorAnalysis_SMCocktail_weighted_lumi-4.000_0-3btag_MultiJet.root

To make a .csv file of the yields in the cards directory, run the following command

     python python/WriteDataCard.py -b MultiJet -c config/run2_sideband.config -d cards Datasets/RazorInclusive_SMS-T1bbbb_2J_mGl-1500_mLSP-100_weighted_lumi-4.000_0-3btag_MultiJet.root Datasets/RazorInclusive_SMCocktail_weighted_lumi-3.0_1-3btag_MultiJet.root --print-yields 

