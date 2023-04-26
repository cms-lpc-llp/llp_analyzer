cp ../lists/displacedJetMuonNtuple/V1p19/Data2018_UL/ParkingBPH1_2018D.txt submit/list.txt
file=ParkingBPH1_2018D
ref=ParkingBPH1_2018D_ref
rm -f ${ref}.txt
rm -f ${ref}_dummy.txt
rm -f missing.txt

for num in {0..2073} 
do 
  echo  list_Job${num}_of_2073.json >> ${ref}_dummy.txt
done
sort ${ref}_dummy.txt > ${ref}.txt

grep -v -f  ${file}.txt ${ref}.txt > missing.txt
#sed -e 's/list_Job\(.*\)_of_\(.*\)/\1/'

rm -f ${ref}.txt
rm -f ${ref}_dummy.txt
