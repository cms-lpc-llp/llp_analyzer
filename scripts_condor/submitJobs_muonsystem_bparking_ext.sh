#!/bin/bash

#(NOT ADDED YET)If you want to send to an eos directory put the LFN here
sendToEOS=true
eosdir="/store/group/lpclonglived/B-ParkingLLPs/"
eosprefix="/eos/uscms"

doSubmit=true
version="V1p17_8" #"V1p19_10"
isData="yes" #"yes" or "no"
options=-1
#outputfilename=""
#label="BParking_Source2018"
#label="Razor2016_07Aug2017Rereco"
label="Razor2017_17Nov2017Rereco"
#label="Razor2018_17SeptEarlyReReco"
filesPerJob=5   #40 good for signal, 5 good for Data
#-1 to do all
maxJobs=-1

# For Reg Run2 
mclistdir=${CMSSW_BASE}/src/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Run2Merged
datalistdir=${CMSSW_BASE}/src/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Run2Merged
subdir=${CMSSW_BASE}/src/llp_analyzer/scripts_condor/gitignore/$version

# For B-Parking Data
#mclistdir=${CMSSW_BASE}/src/llp_analyzer/lists/displacedJetMuonNtuple/V1p19/MC_Fall18/v1/sixie
#datalistdir=${CMSSW_BASE}/src/llp_analyzer/lists//displacedJetMuonNtuple/V1p19/Data2018_UL
#subdir=${CMSSW_BASE}/src/llp_analyzer/scripts_condor/gitignore/$version

#for assigning indices to files after splitting
padding=7
printf "Version: ${version}\n"

cd ${CMSSW_BASE}/..
tar --exclude-caches-all --exclude-vcs -zcf CMSSW_9_4_4.tar.gz -C CMSSW_9_4_4/.. CMSSW_9_4_4 --exclude=src --exclude=tmp --exclude="*.root"
cd -
cp ${CMSSW_BASE}/../CMSSW_9_4_4.tar.gz .

samples=(  \
 "SingleElectron_2017B" \
 "SingleElectron_2017C" \
 "SingleElectron_2017C_v2" \
 "SingleElectron_2017D" \
 "SingleElectron_2017E" \
 "SingleElectron_2017F" \
 "SingleMuon_2017B" \
 "SingleMuon_2017C" \
 "SingleMuon_2017D" \
 "SingleMuon_2017E" \
 "SingleMuon_2017F" \
 "SingleMuon_2017G" \
 "SingleMuon_2017H" \
)

# "ParkingBPH1_2018A"      \
# "ParkingBPH2_2018A"      \
# "ParkingBPH3_2018A"      \
# "ParkingBPH4_2018A"      \
# "ParkingBPH5_2018A"      \
# "ParkingBPH6_2018A"      \
# "ParkingBPH1_2018B"      \
# "ParkingBPH2_2018B"      \
# "ParkingBPH3_2018B"      \
# "ParkingBPH4_2018B"      \
# "ParkingBPH5_2018B"      \
# "ParkingBPH6_2018B"      \
# "ParkingBPH1_2018C"      \
# "ParkingBPH2_2018C"      \
# "ParkingBPH3_2018C"      \
# "ParkingBPH4_2018C"      \
# "ParkingBPH5_2018C"      \
# "ParkingBPH1_2018D"      \
# "ParkingBPH2_2018D"      \
# "ParkingBPH3_2018D"      \
# "ParkingBPH4_2018D"      \
# "ParkingBPH5_2018D"      \
# "ParkingBPH1_2018D"      \
# "ParkingBPH2_2018D"      \
# "ParkingBPH3_2018D"      \
# "ParkingBPH4_2018D"      \
# "ParkingBPH5_2018D"      \

# "SingleElectron_2016B" \
# "SingleElectron_2016B_v2" \
# "SingleElectron_2016C" \
# "SingleElectron_2016D" \
# "SingleElectron_2016E" \
# "SingleElectron_2016F" \
# "SingleElectron_2016G" \
# "SingleElectron_2016H" \
# "SingleMuon_2016B" \
# "SingleMuon_2016C" \
# "SingleMuon_2016D" \
# "SingleMuon_2016E" \
# "SingleMuon_2016F" \
# "SingleMuon_2016G" \
# "SingleMuon_2016H" \
# 
# "SingleElectron_2017B" \
# "SingleElectron_2017C" \
# "SingleElectron_2017C_v2" \
# "SingleElectron_2017D" \
# "SingleElectron_2017E" \
# "SingleElectron_2017F" \
# "SingleMuon_2017B" \
# "SingleMuon_2017C" \
# "SingleMuon_2017D" \
# "SingleMuon_2017E" \
# "SingleMuon_2017F" \
# "SingleMuon_2017G" \
# "SingleMuon_2017H" \
# 
# "SingleMuon_2018A" \
# "SingleMuon_2018B" \
# "SingleMuon_2018C" \
# "SingleMuon_2018D" \
# "EGamma_2018A" \
# "EGamma_2018B" \
# "EGamma_2018C" \
# "EGamma_2018D" \
 
mkdir -p $subdir
makeasubmitdir () {
 printf "Making submits for $1\n"
 
 # go to the directory
 submitdir=${subdir}/$1
 mkdir -p ${submitdir} 
 pushd    ${submitdir}  > /dev/null
 printf " The directory is %s\n" $(pwd)
 
 mkdir -p logs
 
 # breaking up input file list
 split -l $filesPerJob -d -a $padding --additional-suffix=.txt $3/${sample}.txt ${sample}"_"
 mkdir -p lists
 mv *.txt lists/
 tar -zcf lists.tgz lists
 
 # write base for submit file
 printf "universe = vanilla\n" > submitfile
 printf "Executable = ${CMSSW_BASE}/src/llp_analyzer/scripts_condor/runJobs_muonsystem_bparking_ext.sh\n" >> submitfile
 printf "Should_Transfer_Files = YES \n" >> submitfile
 printf "WhenToTransferOutput = ON_EXIT\n" >> submitfile
 printf "Transfer_Input_Files = ${submitdir}/lists.tgz,${CMSSW_BASE}/src/llp_analyzer/data/SFs.tar.gz, ${CMSSW_BASE}/src/llp_analyzer/data/PileupWeights.tar.gz ,${submitdir}/lists.tgz,${CMSSW_BASE}/src/llp_analyzer/data/JEC.tar.gz,${CMSSW_BASE}/src/llp_analyzer/scripts_condor/CMSSW_9_4_4.tar.gz,${CMSSW_BASE}/src/llp_analyzer/bin/Runllp_MuonSystem_bparking_ext,${CMSSW_BASE}/src/llp_analyzer/data/PileupWeights/PileupReweight_2018_BParking.root,${CMSSW_BASE}/src/llp_analyzer/data/ScaleFactors/BParking_SF.root,${CMSSW_BASE}/src/llp_analyzer/data/ScaleFactors/METTriggers_SF.root,${CMSSW_BASE}/src/llp_analyzer/data/HiggsPtWeights/ggH_HiggsPtReweight_NNLOPS.root\n" >> submitfile
 printf "\n" >> submitfile
 printf "notify_user = $(whoami)@cern.ch\n" >> submitfile
 printf "x509userproxy = $X509_USER_PROXY\n" >> submitfile
 printf "\n" >> submitfile
 printf "Output = logs/runanalyzer_muonsystem_bparking_\$(Cluster)_\$(Process).stdout\n" >> submitfile
 printf "Error  = logs/runanalyzer_muonsystem_bparking_\$(Cluster)_\$(Process).stderr\n" >> submitfile
 printf "Log    = logs/runanalyzer_muonsystem_bparking_\$(Cluster)_\$(Process).log\n" >> submitfile

 TheEosDir=${eosdir}/${version}/$1
 if [ ${sendToEOS} == true ]; then 
   if [ -d ${eosprefix}${TheEosDir} ]; then
     eos root://cmseos.fnal.gov rm -r ${TheEosDir} 
   fi
   eos root://cmseos.fnal.gov mkdir -p $TheEosDir
   echo making eos dir
   echo ${eosprefix}$TheEosDir
 fi

 # make haddfile (make now for merging expected results)
 haddfile="./haddit.sh"
 printf "#!/bin/bash\n\n" > ${haddfile}    
 # make checker
 checkfile="./checker.sh"
 printf "#!/bin/bash\n\n" > ${checkfile}
 # make cleaner
 cleaner="./cleanup.sh"
 printf "#!/bin/bash\n\n" > ${cleaner}

 jobnrlow=0
 jobnr=0
 haddnr=0

 if [ ${sendToEOS} == true ]; then
   # hadd command, name of the first intermediate merged file
   printf "hadd ${eosprefix}${TheEosDir}/$1-hadd${haddnr}.root"     >>       ${haddfile}
 else
   # hadd command, name of the first intermediate merged file
   printf "hadd ${submitdir}/$1-hadd${haddnr}.root"     >>       ${haddfile}
 fi    

 printf "\n" >> submitfile
 until [ ${jobnrlow} -ge $2 ]; do
  index=$(printf "%0*d" ${padding} ${jobnr});
  listFile=$submitdir/${sample}"_"${index}.txt
  printf "Arguments = $1 $index $isData $label $options $TheEosDir $sendToEOS\n" >> submitfile
  printf "Queue\n" >> submitfile
  printf "\n" >> submitfile

  # add files to be produced to haddfiles
  printf "\\"  >> ${haddfile}    

  if [ ${sendToEOS} == true ]; then
    printf "\n eos root://cmseos.fnal.gov rm ${TheEosDir}/$1_${index}.root"     >> ${cleaner}    
    printf "\n ${eosprefix}${TheEosDir}/$1_${index}.root"     >> ${haddfile}    
    # add file to checker, all histos are made at the same time, so only check one
    printf "\n if [ ! -f ${eosprefix}${TheEosDir}/$1_${index}.root ]; then printf \" ${eosprefix}${TheEosDir}/$1_${index}.root \\n\"; fi " >> ${checkfile}
  else 
    printf "\n rm $(pwd)/$1_${index}.root"     >> ${cleaner}    
    printf "\n $(pwd)/$1_${index}.root"     >> ${haddfile}    
    # add file to checker, all histos are made at the same time, so only check one
    printf "\n if [ ! -f $(pwd)/$1_${index}.root ]; then printf \" $(pwd)/$1_${index}.root \\n\"; fi " >> ${checkfile}
  fi

  if [ $((${jobnr}%190)) == 0 ] && [ ${jobnr} -ne 0 ]; then
    haddnr=$(( ${haddnr} + 1 ))
    printf "\\"  >> ${haddfile}    
    if [ ${sendToEOS} == true ]; then
      # hadd command, name of the nth intermediate merged file
      printf "\n\n" >> ${haddfile}
      printf "hadd ${eosprefix}${TheEosDir}/$1-hadd${haddnr}.root"     >>       ${haddfile}
    else
      # hadd command, name of the nth intermediate merged file
      printf "\n\n" >> ${haddfile}
      printf "hadd ${submitdir}/$1-hadd${haddnr}.root"     >>       ${haddfile}
    fi    
  fi
  # increment filenumber counters
  jobnrlow=$(( ${jobnrlow} + 1 ))
  jobnr=$(( ${jobnr} + 1 ))
 done # until jobnrlow > nTotalJobs

 printf "\n\n\n" >> ${haddfile}
 printf "##--Now Merge the final file \n" >> ${haddfile}    
 # hadd command, name of final merged file
 printf "hadd -f ${eosprefix}${TheEosDir}/$1.root"     >>       ${haddfile}
 iter=0    
 until [ $iter -gt $haddnr ]; do
  printf "\\"  >> ${haddfile}    
  if [ ${sendToEOS} == true ]; then
    printf "\n eos root://cmseos.fnal.gov rm ${TheEosDir}/$1-hadd${iter}.root"     >> ${cleaner}    
    printf "\n ${eosprefix}${TheEosDir}/$1-hadd${iter}.root"     >> ${haddfile}    
  else 
    printf "\n rm $(pwd)/$1-hadd${iter}.root"     >> ${cleaner}    
    printf "\n $(pwd)/$1-hadd${iter}.root"     >> ${haddfile}    
  fi
  iter=$(( ${iter} + 1 ))
 done


 if [ ${doSubmit} = true ]; then
  condor_submit submitfile
 fi
 
 popd > /dev/null
}


# actually call the function
for sample in ${samples[@]}; do
  #if [[ ${sample} == "Parking"* ]]; then 
  #  isData="yes"
  #else
  #  isData="no"
  #fi
  isData="yes"
  #data or MC?
  listdir=$datalistdir
  #if [[ "$isData" == "yes" ]];  then
  #  listdir=$datalistdir
  #else
  #  listdir=$mclistdir
  #fi
  echo is data? $isData

  nlines=`cat ${listdir}/${sample}.txt | wc -l`
  njobs1=$((nlines/filesPerJob))
  njobs2=$((nlines%filesPerJob))
  if [ $njobs2 -gt 0 ]; then
    njobs=$((njobs1+1))
  else  
    njobs=$((njobs1))
  fi

  if [ $maxJobs -eq -1 ]; then 
    nTotalJobs=$((njobs))
  elif [ $maxJobs -lt $njobs ]; then
    nTotalJobs=$maxJobs
  else
    nTotalJobs=$((nTotalJobs+njobs))
  fi
  echo Total number of jobs for $sample : $nTotalJobs
  makeasubmitdir ${sample} ${nTotalJobs} ${listdir}
  echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
done
