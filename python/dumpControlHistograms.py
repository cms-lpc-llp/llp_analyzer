### INPUT: ROOT file produced by RazorAnalyzer python macro code

import sys
import argparse
import ROOT as rt
from ast import literal_eval

from macro import macro

parser = argparse.ArgumentParser()
parser.add_argument("fname",help="input ROOT file name")
parser.add_argument("varname", help="variable to print")
parser.add_argument("--all", dest="print_all", action='store_true',
        help="Print individual MC histograms")
parser.add_argument("--debug", action="store_true", help="print debug info")
args = parser.parse_args()
hists = macro.importHists(args.fname, debugLevel=args.debug)
# parse variable name if tuple
try:
    parsed_var = literal_eval(args.varname)
except ValueError:
    parsed_var = args.varname

tot = hists["Data"][parsed_var].Clone("tot")
tot.Reset()
for proc in hists:
    if args.print_all or proc == "Data":
        print proc
        hists[proc][parsed_var].Print("all")
    else:
        tot.Add(hists[proc][parsed_var])
print "Total MC"
tot.Print("all")
