from __future__ import print_function
from ROOT import TH2F, TFile, TCanvas
from subprocess import call
infile = TFile.Open("Efficiency_RazorTrigger_RsqMR270_Rsq0p09_MR200_All_SingleLeptonData.root")
assert(infile)
eff = infile.Get("Efficiency_MRRsq")

with open("EfficiencyTable.tex",'w') as outfile:

    outfile.write('\\documentclass{article}\n')
    outfile.write('\\usepackage{pdflscape}\n')
    outfile.write('\\usepackage{slashbox}\n')
    outfile.write('\\usepackage[margin=1in]{geometry}\n')
    outfile.write('\\usepackage{graphicx}\n')
    outfile.write('\\usepackage{multicol}\n')
    outfile.write('\\begin{document}\n')
    outfile.write('All cuts: $R^2 > 0.15, M_R > 150$. At least 2 jets.\n\n') 
    outfile.write('Denominator: Pass SingleLepton triggers: \n\\begin{multicols}{2}\n\\begin{itemize}\n')
    outfile.write('\\item        \"HLT\\_IsoMu20\"\n \\item        \"HLT\\_IsoTkMu20\"\n \\item        \"HLT\\_IsoMu22\"\n       \\item        \"HLT\\_IsoTkMu22\"\n        \\item        \"HLT\\_IsoMu24\"\n        \\item        \"HLT\\_IsoTkMu24\"\n        \\item        \"HLT\\_IsoMu27\"\n         \\item        \"HLT\\_IsoTkMu27\"\n         \\item        \"HLT\\_Mu50\"\n        \\item        \"HLT\\_Ele23\\_WPLoose\\_Gsf\"\n        \\item        \"HLT\\_Ele27\\_WPLoose\\_Gsf\"\n        \\item        \"HLT\\_Ele27\\_WPTight\\_Gsf\"\n        \\item        \"HLT\\_Ele27\\_eta2p1\\_WPLoose\\_Gsf\"\n        \\item        \"HLT\\_Ele27\\_eta2p1\\_WPTight\\_Gsf\"\n        \\item        \"HLT\\_Ele32\\_eta2p1\\_WPLoose\\_Gsf\"\n        \\item        \"HLT\\_Ele32\\_eta2p1\\_WPTight\\_Gsf\"\n        \\item        \"HLT\\_Ele105\\_CaloIdVT\\_GsfTrkIdT\"\n        \\item        \"HLT\\_Ele115\\_CaloIdVT\\_GsfTrkIdT\"\n')
    outfile.write('\\end{itemize}\n\\end{multicols}\n')
    outfile.write('Numerator: Pass SingleLepton and Razor Trigger: HLT\\_RsqMR270\\_Rsq0p09\\_MR200.\n\n')
    outfile.write('\\begin{center}\n')
    outfile.write('\\includegraphics[width=0.4\\textwidth]{Efficiency_RazorTrigger_RsqMR270_Rsq0p09_MR200_All_SingleLeptonData_MR.png}\n')
    outfile.write('\\includegraphics[width=0.4\\textwidth]{Efficiency_RazorTrigger_RsqMR270_Rsq0p09_MR200_All_SingleLeptonData_Rsq.png}\\\\\n')
    outfile.write('\\includegraphics[width=0.6\\textwidth]{Efficiency_RazorTrigger_RsqMR270_Rsq0p09_MR200_All_SingleLeptonData_MRRsq.png}\n')
    outfile.write('\\end{center}\n')
    outfile.write('\\newpage\n')
    outfile.write('\\begin{landscape}\n')
    outfile.write('\\centering')
    outfile.write('\\begin{tabular} { | ')
    outfile.write(' p{2.1cm} | ')
    for ix in range(1,int(eff.GetNbinsX()/2)+1):
        outfile.write(' p{1.8cm} | ')
    outfile.write('}\n\\hline\n')
    outfile.write ('\\backslashbox{$R^2$}{$M_R$} ')
    for ix in range(1,int(eff.GetNbinsX()/2)+1):
        outfile.write (" & [%.2f, %.2f] " % (eff.GetXaxis().GetBinLowEdge(ix), eff.GetXaxis().GetBinLowEdge(ix+1)))
    outfile.write ('\\\\\\hline\n')
    for iy in range(1,eff.GetNbinsY()+1):
        outfile.write ("{[%.2f, %.2f]}" % (eff.GetYaxis().GetBinLowEdge(iy), eff.GetYaxis().GetBinLowEdge(iy+1)))
        for ix in range(1,int(eff.GetNbinsX()/2)+1):
            outfile.write (" & %.2f $\\pm$ %.2f" % (eff.GetBinContent(ix,iy), eff.GetBinErrorLow(ix,iy)))
        outfile.write ('\\\\\n')
    outfile.write('\\hline\n')
    outfile.write('\\end{tabular}\n\n')
    outfile.write('\\begin{tabular} { | ')
    outfile.write(' p{2.1cm} | ')
    for ix in range(int(eff.GetNbinsX()/2)+1,eff.GetNbinsX()+1):
        outfile.write(' p{1.8cm} | ')
    outfile.write('}\n\\hline\n')
    outfile.write ('\\backslashbox{$R^2$}{$M_R$} ')
    for ix in range(int(eff.GetNbinsX()/2)+1,eff.GetNbinsX()+1):
        outfile.write (" & [%.2f, %.2f] " % (eff.GetXaxis().GetBinLowEdge(ix), eff.GetXaxis().GetBinLowEdge(ix+1)))
    outfile.write ('\\\\\\hline\n')
    for iy in range(1,eff.GetNbinsY()+1):
        outfile.write ("{[%.2f, %.2f]}" % (eff.GetYaxis().GetBinLowEdge(iy), eff.GetYaxis().GetBinLowEdge(iy+1)))
        for ix in range(int(eff.GetNbinsX()/2)+1,eff.GetNbinsX()+1):
            outfile.write (" & %.2f $\\pm$ %.2f" % (eff.GetBinContent(ix,iy), eff.GetBinErrorLow(ix,iy)))
        outfile.write ('\\\\\n')
    outfile.write('\\hline\n')
    outfile.write('\\end{tabular}\n')
    outfile.write('\\end{landscape}\n')
    outfile.write('\\end{document}')

call(['pdflatex','EfficiencyTable.tex'])
