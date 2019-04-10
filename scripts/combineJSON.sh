#!/bin/bash -f
# usage createGoodList.sh CERT.json

certjson=$1

# combine all processed lumi summary files from all datasets
compareJSON.py --or RazorAnalyzer/lists/razorNtuplerV1p5-Run1/lumiSummary/lumiSummary_DoubleMuParked_Run2012A.json RazorAnalyzer/lists/razorNtuplerV1p5-Run1/lumiSummary/lumiSummary_DoubleMuParked_Run2012B.json tmpAB.json
compareJSON.py --or RazorAnalyzer/lists/razorNtuplerV1p5-Run1/lumiSummary/lumiSummary_DoubleMuParked_Run2012C.json RazorAnalyzer/lists/razorNtuplerV1p5-Run1/lumiSummary/lumiSummary_DoubleMuParked_Run2012D.json tmpCD.json
compareJSON.py --or tmpAB.json  tmpCD.json RazorAnalyzer/lists/razorNtuplerV1p5-Run1/lumiSummary/lumiSummary_DoubleMuParked_Run2012.json


# and do the AND of the lumi summary files of all datasets that we want
compareJSON.py --and DoubleMu.json DoubleElectron.json output1.json
compareJSON.py --and MuEG.json output1.json CombinedPDs.json

# and now do the OR with the certification json and calc the lumi
lumiCalc2.py -c frontier://LumiCalc/CMS_LUMI_PROD -i CombinedPDs.json --nowarning overview >& CombinedPDs.lumi &
compareJSON.py --and $certjson CombinedPDs.json CombinedPDs_Good.json
lumiCalc2.py -c frontier://LumiCalc/CMS_LUMI_PROD -i CombinedPDs_Good.json --nowarning overview >& CombinedPDs_Good.lumi &

echo "END. Final jsons are:"
echo "---> from crab report CombinedPDs.json (lumi = CombinedPDs.lumi)"
echo "---> from crab report in AND with DQM json CombinedPDs_Good.json (lumi = CombinedPDs_Good.lumi)"
