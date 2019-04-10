from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
from DustinTuple2RooDataSet import initializeWorkspace

#seed = 1999 # fits for 10-24-2015
#seed = 1989
#seed = 555 # fits for 11-28-2015
seed = 1912
def convertDataset2UnweightedToy(data, cfg, box, workspace, uwName = 'uw'):
    """Get the cocktail dataset from the file"""
    row = data.get()

    MR = row['MR']
    Rsq = row['Rsq']
    nBtag = row['nBtag']
    
    varSet = rt.RooArgSet(MR,Rsq,nBtag)
    varList = rt.RooArgList(MR,Rsq,nBtag)
    varList2D = rt.RooArgList(MR,Rsq)
    uwdata = rt.RooDataSet(uwName+'tree','Unweighted Cocktail',varSet)
        
    mRmin = row['MR'].getMin()
    mRmax = row['MR'].getMax()
    rsqMin = row['Rsq'].getMin()
    rsqMax = row['Rsq'].getMax()
    btagMin = row['nBtag'].getMin()
    btagMax = row['nBtag'].getMax()
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning

    # use fine binning
    rt.RooRandom.randomGenerator().SetSeed(seed)
    #myTH3 = rt.TH3D(uwName+box, uwName+box, 100, mRmin, mRmax, 100, rsqMin, rsqMax, int(btagMax-btagMin), btagMin, btagMax)
    #myTH2 = rt.TH2D(uwName+box+"2d", uwName+box+"2d", 100, mRmin, mRmax, 100, rsqMin, rsqMax)
    #myTH2Toy = rt.TH2D("h", "h", 100, mRmin, mRmax, 100, rsqMin, rsqMax)

    # use binning written in config
    myTH3 = rt.TH3D(uwName+box, uwName+box, len(x)-1, x, len(y)-1, y, len(z)-1, z)
    myTH2 = rt.TH2D(uwName+box+"2d", uwName+box+"2d", len(x)-1, x, len(y)-1, y)
    myTH2Toy = rt.TH2D("h", "h", len(x)-1, x, len(y)-1, y)
    myTH2.Sumw2()
    myTH2Toy.Sumw2()

    # fills automatically with weight
    data.fillHistogram(myTH3, varList,"MR>0")
    data.fillHistogram(myTH2, varList2D,"MR>0")

    # fix for negative weight bins
    for i in range(1,myTH3.GetNbinsX()+1):        
        for j in range(1,myTH3.GetNbinsY()+1):            
            for k in range(1,myTH3.GetNbinsZ()+1):
                print i, j, k, myTH3.GetBinContent(i,j,k)
                if myTH3.GetBinContent(i,j,k) < 0.:
                    myTH3.SetBinContent(i,j,k,0.)
                    
    c = rt.TCanvas("c","c",600,400)
    #rt.gStyle.SetOptStat(1001000011)
    rt.gStyle.SetOptStat(0)
    myTH2.SetTitle("Weighted %s"%box)
    sumW2 = 0
    for i in range(0,wdata.numEntries()):
       wdata.get(i)
       sumW2+=(wdata.weight())*(wdata.weight())

    print "sum (weights)^2 = %.1f" %sumW2
    print "(sum weights)^2 = %.1f" %((wdata.sumEntries())*(wdata.sumEntries()))
    effEntries = (((wdata.sumEntries())*(wdata.sumEntries()))/sumW2)
    print "effective entries = %.1f"%effEntries
    myTH2.GetXaxis().SetTitle("M_{R}")
    myTH2.GetYaxis().SetTitle("R^{2}")
    myTH2.GetXaxis().SetMoreLogLabels()
    myTH2.GetYaxis().SetMoreLogLabels()
    myTH2.GetXaxis().SetNoExponent()
    myTH2.GetYaxis().SetNoExponent()
    myTH2.Draw("colz")
    c.SetLogy()
    c.SetLogx()

    
    #l = rt.TLatex()
    #l.SetTextAlign(11)
    #l.SetTextSize(0.045)
    #l.SetTextFont(42)
    #l.SetNDC()
    #l.DrawLatex(0.18,0.84,"CMS Simulation 5 fb  ^{-1} (13 TeV)")
    #l.DrawLatex(0.18,0.77,"Razor MultiJet Box, QCD")
    #l.DrawLatex(0.11,0.84,"CMS Simulation 5 fb  ^{-1} (13 TeV)")
    #l.DrawLatex(0.11,0.77,"Razor MultiJet Box")
    #l.DrawLatex(0.56,0.84,"pp #rightarrow #tilde{g}#tilde{g},  #tilde{g}#rightarrowb#bar{b}#tilde{#chi}^{0}_{1}")
    #l.DrawLatex(0.5,0.77,"m_{#tilde{g}} = %i GeV, m  _{#tilde{#chi}} = %i GeV"%(1000,900))

    if btagMax>btagMin+1:
        c.Print(options.outDir+"/TH2D_SMCocktail_weighted_%i-%ibtag_%s.pdf"%(btagMin,btagMax-1,box))
        c.Print(options.outDir+"/TH2D_SMCocktail_weighted_%i-%ibtag_%s.C"%(btagMin,btagMax-1,box))
    else:
        c.Print(options.outDir+"/TH2D_SMCocktail_weighted_%ibtag_%s.pdf"%(btagMin,box))
        c.Print(options.outDir+"/TH2D_SMCocktail_weighted_%ibtag_%s.C"%(btagMin,box))

    print wdata.weight()
    Nev = myTH3.Integral()
    Nent = myTH3.GetEntries()
    print "weighted events %.1f"% Nev
    print "entries  %d"% Nent
    Npois = rt.RooRandom.randomGenerator().Poisson(Nev)
    print "Npois = %d "%Npois
    
    #wdata2d = wdata.reduce(rt.RooArgSet(MR,Rsq),"MR>500&&Rsq>0.3&&nBtag==3")
    #rookeys = rt.RooNDKeysPdf("rookeys", "rookeys", rt.RooArgList(MR,Rsq), wdata2d, "am")
    #uwdata = rookeys.generate(rt.RooArgSet(MR,Rsq),Npois)
    
    for i in range(0,Npois):
       myMR = rt.Double()
       myRsq = rt.Double()
       mynBtag = rt.Double()
       myTH3.GetRandom3(myMR,myRsq,mynBtag)
       mynBtag = int(mynBtag)
       varSet.setRealValue('MR',myMR)
       varSet.setRealValue('Rsq',myRsq)
       varSet.setRealValue('nBtag',mynBtag)
       uwdata.add(varSet)
    

    uwdata.fillHistogram(myTH2Toy, varList2D,"MR>0")
    myTH2Toy.SetTitle("Unweighted %s"%box)
    myTH2Toy.GetXaxis().SetTitle("M_{R}")
    myTH2Toy.GetYaxis().SetTitle("R^{2}")
    myTH2Toy.GetXaxis().SetMoreLogLabels()
    myTH2Toy.GetYaxis().SetMoreLogLabels()
    myTH2Toy.GetXaxis().SetNoExponent()
    myTH2Toy.GetYaxis().SetNoExponent()
    myTH2Toy.Draw("colz")
    if btagMax>btagMin+1:
        c.Print(options.outDir+"/TH2D_SMCocktail_unweighted_%i-%ibtag_%s.pdf"%(btagMin,btagMax-1,box))
        c.Print(options.outDir+"/TH2D_SMCocktail_weighted_%i-%ibtag_%s.C"%(btagMin,btagMax-1,box))
    else:
        c.Print(options.outDir+"/TH2D_SMCocktail_unweighted_%ibtag_%s.pdf"%(btagMin,box))
        c.Print(options.outDir+"/TH2D_SMCocktail_unweighted_%ibtag_%s.C"%(btagMin,box))

    return uwdata


    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store datasets")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box =  options.box

    ids = 0
    
    w = rt.RooWorkspace("w"+box)

    initializeWorkspace(w,cfg,box)
    
    ds = []
    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            workspace = rootFile.Get('w'+box)
            ds.append(workspace.data('RMRTree').Clone('RMRTree_%i'%ids))
            ids+=1

    
    wdata = ds[0].Clone('RMRTree')
    for ids in range(1,len(ds)):
        wdata.append(ds[ids])
    
    uwdata = convertDataset2UnweightedToy(wdata, cfg, box, w, uwName = 'uw')

    uwdata.SetName('RMRTree')
    rootTools.Utils.importToWS(w,uwdata)

    
    inFiles = [f for f in args if f.lower().endswith('.root')]
            
    if len(inFiles)==1:
        outFile = inFiles[0].split('/')[-1].replace('weighted','unweighted')
        
    outFile = rt.TFile.Open(options.outDir+"/"+outFile,'recreate')
    outFile.cd()
    w.Write()
    outFile.Close()
    
            
