import sys,os
import argparse
import copy
import ROOT as rt

#local imports
from framework import Config
import macro.plotting as plotting

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('--num-pdf-weights',dest="numPdfWeights",default=60,type=int, help="Number of nuisance parameters to use for PDF uncertainties")
    parser.add_argument('--compute-pdf-envelope',dest='computePdfEnvelope',default=False,action='store_true', help='Collapse pdf variations into a single nuisance parameter')
    parser.add_argument("filename", help="Path to input file")
    args = parser.parse_args() 

    #get input file
    infile = rt.TFile(args.filename)
    assert infile

    #specify shape uncertainties
    shapes = ['muoneff','eleeff','jes','muontrig','eletrig','btag','muonfastsim','elefastsim','btagfastsim','facscale','renscale','facrenscale','ees','mes']
    if args.computePdfEnvelope:
        shapes.append('pdfenvelope')
    else:
        shapes.extend([str(n)+'pdf' for n in range(args.numPdfWeights)])
    titles = {'muoneff':'Muon efficiency', 'eleeff':'Electron efficiency', 'jes':'Jet energy scale',
            'muontrig':'Muon trigger eff.', 'eletrig':'Electron trigger eff.',
            'btag':'B-tagging efficiency', 'muonfastsim':'Fastsim muon eff.',
            'elefastsim':'Fastsim electron eff.', 'btagfastsim':'Fastsim b-tagging eff.',
            'facscale':'Factorization scale', 'renscale':'Renormalization scale',
            'facrenscale':'Diag. scale variation', 'ees':'Electron energy scale',
            'mes':'Muon energy scale', 'pdfenvelope':'PDF variations'}
    for n in range(args.numPdfWeights):
        titles[str(n)+'pdf'] = 'Pdf variation '+str(n)

    #parse model and box name from file
    split = args.filename.replace('.root','').split('_')
    model = split[0].replace('SMS-','')
    box = split[-1]

    #make output directory
    dirName = "SMSPlots"
    os.system('mkdir -p '+dirName)

    #load histograms and plot
    hists = {}
    hists['nominal'] = infile.Get('_'.join([box,model]))
    print "Loading histogram",'_'.join([box,model])
    c = rt.TCanvas('c','c',800,600)
    assert hists['nominal']
    for shape in shapes:
        for updown in ['Up','Down']:
            hists[shape+updown] = infile.Get('_'.join([box,model,shape+updown]))
            print "Loading histogram",'_'.join([box,model,shape+updown])
            assert hists[shape+updown]

        #print some statistics about the size of this systematic
        maxDeviation = 0
        avgDeviation = 0
        for bx in range(1, hists['nominal'].GetNbinsX()+1):
            if hists['nominal'].GetBinContent(bx) > 0:
                percentDiffUp =  abs(hists[shape+'Up'].GetBinContent(bx) - hists['nominal'].GetBinContent(bx))/hists['nominal'].GetBinContent(bx)
                percentDiffDown =  abs(hists[shape+'Down'].GetBinContent(bx) - hists['nominal'].GetBinContent(bx))/hists['nominal'].GetBinContent(bx)
                if percentDiffUp > maxDeviation: maxDeviation = percentDiffUp
                if percentDiffDown > maxDeviation: maxDeviation = percentDiffDown
                avgDeviation = (2*bx*avgDeviation + percentDiffDown + percentDiffUp)/(2*(bx+1))
        print ('\n'+shape),": \nMax deviation",maxDeviation,'\nAvg deviation',avgDeviation

        toplot = [hists['nominal'],hists[shape+'Up'],hists[shape+'Down']]
        colors = [rt.kBlack, rt.kGreen, rt.kBlue]
        title = titles[shape]
        d = {0:toplot[0],1:toplot[1],2:toplot[2]}
        t = {0:'Nominal',1:title+' Up',2:title+' Down'}
        ordering = [0,1,2]
        if shape == 'pdfenvelope': 
            toplot.pop()
            colors.pop()
            ordering.pop()
            del t[2]
            del d[2]
        leg = plotting.makeLegend(d, t, ordering, 0.4, 0.0, 0.6, 0.3)

        plotting.plot_several(c, toplot, leg, colors, 'Bin', 'A.U.', printstr=box+model+shape, lumistr='', commentstr=model+' '+box, printdir=dirName)
