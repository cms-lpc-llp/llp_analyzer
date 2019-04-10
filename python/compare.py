import argparse
import ROOT as rt
from macro import macro

titles={
        "Fastsim": "Fastsim",
        "Fullsim": "Fullsim",
        }

def compare(c, fi1, fi2, var, mini, maxi, bins, log=False, cut="", name1="Fastsim", name2="Fullsim", scaleHistsToMatch=True, treeName="RazorInclusive"):
    #get
    f1 = rt.TFile(fi1)
    f2 = rt.TFile(fi2)
    t1 = f1.Get(treeName)
    t2 = f2.Get(treeName)

    #draw
    h1 = rt.TH1F("h1", name1, bins, mini, maxi)
    h2 = rt.TH1F("h2", name2, bins, mini, maxi)
    h1.Sumw2()
    h2.Sumw2()
    t1.Draw(var+'>>h1', cut)
    t2.Draw(var+'>>h2', cut)

    #scale
    if scaleHistsToMatch:
        scale = t1.GetEntries()*1.0/t2.GetEntries()
        h2.Scale(scale)

    #min/max
    ymin = None
    if not log: ymin = min(h1.GetMinimum(), h2.GetMinimum())
    ymax = max(h1.GetMaximum(), h2.GetMaximum())

    #plot
    macro.makeStackAndPlot(c, {name2:h2}, h1, dataName=name1, mcOrdering=[name2], titles=titles, mcTitle=var, xtitle=var, printstr="compare_"+(var.replace('$','').replace('(','').replace(')','').replace('/',''))+cut, logy=log, lumistr='', ytitle="A.U.", savepng=False, savepdf=True, saveroot=True) 

if __name__ == '__main__':
    rt.gROOT.SetBatch()

    parser = argparse.ArgumentParser()
    parser.add_argument('listfile')
    parser.add_argument('f1')
    parser.add_argument('f2')
    parser.add_argument('-t', '--tree', default='RazorInclusive', help='name of input tree')
    args = parser.parse_args()

    c = rt.TCanvas('c','c',800,600)
    with open(args.listfile) as f:
        for line in f:
            splitLine = line.split()
            var = splitLine[0]
            mini = float(splitLine[1])
            maxi = float(splitLine[2])
            bins = 50
            log = False
            cut = ""
            if len(splitLine) > 3: bins = int(splitLine[3])
            if len(splitLine) > 4: log = int(splitLine[4])
            if len(splitLine) > 5: cut = splitLine[5]
            compare(c, args.f1, args.f2, var, mini, maxi, bins, log, cut, treeName=args.tree)
