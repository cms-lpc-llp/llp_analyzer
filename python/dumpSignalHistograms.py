### INPUT: ROOT file produced by RazorAnalyzer python macro code

import sys
import math
import argparse
import ROOT as rt

from macro import macro, plotting

parser = argparse.ArgumentParser()
parser.add_argument("fname",help="input ROOT file name")
parser.add_argument("--all",action='store_true',dest="print_all",
        help="Print all MC histograms")
args = parser.parse_args()
hists = macro.importHists(args.fname)
var = ('MR','Rsq')

tmp = plotting.unroll2DHistograms([hists[hists.iterkeys().next()][var]])[0]
tot = tmp.Clone("tot")
tot.Reset()
sys_dict = {}
qcdHist = None
for proc in hists:
    if proc == 'Sys':
        for sys_proc in hists[proc]:
            if sys_proc == 'Data': continue
            for sys in hists[proc][sys_proc]:
                if sys not in sys_dict:
                    tmp = plotting.unroll2DHistograms([hists['Data'][var]])[0]
                    sys_dict[sys] = tmp.Clone(sys)
                    sys_dict[sys].Reset()
                tmp = plotting.unroll2DHistograms([hists['Sys'][sys_proc][sys][var]])[0]
                for ibin in range(1, tmp.GetNbinsX()+1):
                    if math.isnan(tmp.GetBinContent(ibin)):
                        tmp.SetBinContent(ibin, 0)
                sys_dict[sys].Add(tmp)
    elif proc == "Data" or proc == "Fit":
        print proc
        tmp = plotting.unroll2DHistograms([hists[proc][var]])[0]
        tmp.Print("all")
    else:
        tmp = plotting.unroll2DHistograms([hists[proc][var]])[0]
        tot.Add(tmp)
        if args.print_all:
            tmp.Print("all")
        if proc == 'QCD':
            qcdHist = tmp
print "Total MC"
tot.Print("all")

print "QCD fraction"
if qcdHist is not None:
    qcdFraction = qcdHist.Clone("qcdFraction")
    qcdFraction.Divide(tot)
    qcdFraction.Print("all")

# Compute the size of each systematic uncertainty,
# averaged across all analysis bins
for sys,hist in sys_dict.iteritems():
    if 'Up' in sys:
        downHist = sys_dict[ sys.replace('Up','Down') ]
        hist.Add(downHist, -1)
        hist.Divide(tot)
        print "Systematic:",sys.replace('Up',''),
        print hist.Integral()/hist.GetNbinsX()
