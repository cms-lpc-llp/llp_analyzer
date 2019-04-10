rm -rf weights weights_PU0_NoTiming
root -l -q TMVA_vtx.C+(\"HggRazorUpgradeTiming_PU0_NoTiming_vtx.root\",\"TMVA_PU0_NoTiming_vtx.root\")
mv weights weights_PU0_NoTiming
rm -rf weights weights_PU140_NoTiming
root -l -q TMVA_vtx.C+(\"HggRazorUpgradeTiming_PU140_NoTiming_vtx.root\",\"TMVA_PU140_NoTiming_vtx.root\")
mv weights weights_PU140_NoTiming
rm -rf weights weights_PU0_Timing_vtx
root -l -q TMVA_vtx.C+(\"HggRazorUpgradeTiming_PU0_Timing_vtx.root\",\"TMVA_PU0_Timing_vtx.root\")
mv weights weights_PU0_Timing_vtx
rm -rf weights weights_PU140_Timing_vtx
root -l -q TMVA_vtx.C+(\"HggRazorUpgradeTiming_PU140_Timing_vtx.root\",\"TMVA_PU140_Timing_vtx.root\")
mv weights weights_PU140_Timing_vtx
