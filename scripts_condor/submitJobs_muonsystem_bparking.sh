#!/bin/bash

doSubmit=true
version="V1p19_0"
isData=false
options=-1
outputfilename=""
label="Razor2018_17SeptEarlyReReco"
filesPerJob=5
maxJobs=1

mclistdir=${CMSSW_BASE}/src/llp_analyzer/lists/displacedJetMuonNtuple/V1p19/MC_Fall18/v1/sixie
datalistdir=${CMSSW_BASE}/src/llp_analyzer/lists//displacedJetMuonNtuple/V1p19/Data2018_UL
subdir=${CMSSW_BASE}/src/llp_analyzer/scripts_condor/gitignore/$version

padding=5
printf "Version: ${version}\n"

cd ${CMSSW_BASE}/..
tar --exclude-caches-all --exclude-vcs -zcf CMSSW_9_4_4.tar.gz -C CMSSW_9_4_4/.. CMSSW_9_4_4 --exclude=src --exclude=tmp --exclude="*.root"
cd -
cp ${CMSSW_BASE}/CMSSW_9_4_4.tar.gz .

samples=(  \
 "BToKPhi_MuonLLPDecayGenFilter_PhiToPi0Pi0_mPhi0p3_ctau300"      \
# "ParkingBPH4_2018A"      \
# "BToKPhi_MuonLLPDecayGenFilter_PhiToPiPlusPiMinus_mPhi0p3_ctau300"      \
)

 
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
 printf "Executable = ${CMSSW_BASE}/src/llp_analyzer/scripts_condor/runJobs_muonsystem_bparking.sh\n" >> submitfile
 printf "Should_Transfer_Files = YES \n" >> submitfile
 printf "WhenToTransferOutput = ON_EXIT\n" >> submitfile
 printf "Transfer_Input_Files = ${submitdir}/lists.tgz,${CMSSW_BASE}/src/llp_analyzer/scripts_condor/CMSSW_9_4_4.tar.gz,${CMSSW_BASE}/src/llp_analyzer/Runllp_MuonSystem_bparking,${CMSSW_BASE}/src/llp_analyzer/data/PileupWeights/PileupReweight_MC_Fall18_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root,${CMSSW_BASE}/src/llp_analyzer/data/ScaleFactors/BParking_SF.root,${CMSSW_BASE}/src/llp_analyzer/data/ScaleFactors/METTriggers_SF.root,${CMSSW_BASE}/src/llp_analyzer/data/HiggsPtWeights/ggH_HiggsPtReweight_NNLOPS.root\n" >> submitfile
 printf "\n" >> submitfile
 printf "notify_user = $(whoami)@cern.ch\n" >> submitfile
 printf "x509userproxy = $X509_USER_PROXY\n" >> submitfile
 printf "\n" >> submitfile
 printf "Output = logs/runanalyzer_muonsystem_bparking_\$(Cluster)_\$(Process).stdout\n" >> submitfile
 printf "Error  = logs/runanalyzer_muonsystem_bparking_\$(Cluster)_\$(Process).stderr\n" >> submitfile
 printf "Log    = logs/runanalyzer_muonsystem_bparking_\$(Cluster)_\$(Process).log\n" >> submitfile

 # make haddfile (make now for merging expected results)
 haddfile="./haddit.sh"
 
 hadddir="${rootdir}/${aversion}"
 mkdir -p ${hadddir}
 printf "#!/bin/bash\n\n" > ${haddfile}    

 # make checker
 checkfile="./checker.sh"
 printf "#!/bin/bash\n\n" > ${checkfile}

 # hadd command, name of final merged file
 printf "hadd ${hadddir}/$1.root"     >>       ${haddfile}    

 jobnrlow=0
 jobnr=0

 printf "\n" >> submitfile
 until [ ${jobnrlow} -ge $2 ]
 do
  index=$(printf "%0*d" ${padding} ${jobnr});
  listFile=$submitdir/${sample}"_"${index}.txt
  echo $listFile
  printf "Arguments = $1 $index $isData $label $option\n" >> submitfile
  printf "Queue\n" >> submitfile
  printf "\n" >> submitfile

  # add files to be produced to haddfiles
  printf "\\"  >> ${haddfile}    

  printf "\n $(pwd)/$1_${index}.root"     >> ${haddfile}    

  # add file to checker, all histos are made at the same time, so only check one
  printf "\n if [ ! -f $(pwd)/$1_${index}.root ]; then printf \" $(pwd)/$1_${index}.root \\n\"; fi " >> ${checkfile}

  # increment filenumber counters
  #printf "NFILES: %s %s %s\n" $nfilesinlist $filenrlow $jobfilenr
  jobnrlow=$(( ${jobnrlow} + 1 ))
  jobnr=$(( ${jobnr} + 1 ))

 done # until jobnrlow > nTotalJobs

 printf "\n\n" >> ${haddfile}    

 if [ ${doSubmit} = true ]
 then
  condor_submit submitfile
 fi
 
 popd > /dev/null
}


# actually call the function
for sample in ${samples[@]} 
do
  if [[ ${sample} == "Parking"* ]]
  then 
    isData=true
  else
    isData=false
  fi
  #data or MC?
  listdir=./
  if [ $isData == true ]
  then
    listdir=$datalistdir
  else
    listdir=$mclistdir
  fi
  echo is data? $isData
  echo listdir $listdir
  nlines=`cat ${listdir}/${sample}.txt | wc -l`
  njobs1=$((nlines/filesPerJob))
  njobs2=$((nlines%filesPerJob))
  if [ $njobs2 -gt 0 ]
  then
    echo adding one
    njobs=$((njobs1+1))
  else  
    echo not adding one
    njobs=$((njobs1))
  fi
  nTotalJobs=0
  if [ $maxJobs -lt $((nTotalJobs+njobs)) ]
  then 
    nTotalJobs=$maxJobs
  else
    nTotalJobs=$((nTotalJobs+njobs))
  fi
  echo Total numer of jobs for $sample : $nTotalJobs
  makeasubmitdir ${sample} ${nTotalJobs} ${listdir}
  echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
done
