
for cat in MultiJet_0b MultiJet_1b MultiJet_2b LeptonJet_0b LeptonJet_1b LeptonJet_2b LeptonMultiJet_0b LeptonMultiJet_1b LeptonMultiJet_2b DiJet_0b DiJet_1b Dijet_2b #DiJet MultiJet LeptonJet LeptonMultiJet
do
    python python/RazorFitMacro.py --config config/run2_2017_03_13_SeparateBtagFits_forToys.config $cat > ${cat}.txt
    python python/RazorFitMacro.py --config config/run2_2017_03_13_SeparateBtagFits_forToys.config --input-fit-file Plots/Razor2016_MoriondRereco/Fits_${cat}/RazorFitInstance_Razor2016_MoriondRereco_${cat}.root --run-toys --no-fit --do-freq $cat 
done
for cat in DiJet LeptonJet LeptonMultiJet MultiJet
do
    python python/CombineBTagToys.py --config config/run2_2017_03_13_SeparateBtagFits_forToys.config $cat
done