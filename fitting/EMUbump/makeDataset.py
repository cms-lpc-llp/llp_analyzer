import ROOT as rt
rt.gSystem.Load('libRooFit')
rt.gSystem.Load('lib/libRazor.so')

if __name__ == '__main__':

    RsqCut = "0.1"
    
    inFile = rt.TFile.Open("data.root")
    tree = inFile.Get('tree')
    print "Reading %i events from tree" %tree.GetEntries()
    MR = rt.RooRealVar("MR","MR",450., 400., 40000.)
    Rsq = rt.RooRealVar("Rsq","Rsq",0.4, float(RsqCut), 1.0)
    mll = rt.RooRealVar("mll","mll",200.,150.,1200.)
    varset = rt.RooArgSet(MR, Rsq, mll)
    dataset = rt.RooDataSet('data','dataset',tree, varset, rt.RooFormulaVar("selection","","1>0",rt.RooArgList()))

    # build the pdf
    ws = rt.RooWorkspace("ws","ws")
    getattr(ws,'import')(MR)
    getattr(ws,'import')(Rsq)
    getattr(ws,'import')(mll)
    #ws.factory("RooCruijff::SigPdf(mll,m0s[200.,100.,800.], sLs[1.,0.5,20.], sRs[1.,0.5,20.], aLs[0.0001,0.,10.], aRs[0.0001,0.,10.])")
    ws.factory("RooGaussian::SigPdf(mll,m0s[400.,200.,800.], s[20.,0.5,20.])")
    ws.factory("RooExponential::BkgPdf(mll,k[-0.1 -100., -0.001])")    
    ws.factory("SUM::extLik(Nsig[0.,0,100]*SigPdf,Nbkg[%i,100,300000]*BkgPdf)" %dataset.sumEntries())
    ws.factory("SUM::bkgLik(Nbkg[%i,100,300000]*BkgPdf)" %dataset.sumEntries())
    ws.Print()

    # draw the data
    frame = mll.frame()
    dataset.plotOn(frame)

    # bkg-only fit
    bkgLik =ws.pdf('bkgLik')
    bkgOnlyFitResult = bkgLik.fitTo(dataset,rt.RooFit.Save(),rt.RooFit.Extended(rt.kTRUE))
        
    # fit the data with S+B
    pdf = ws.pdf('extLik')
    fitResult = pdf.fitTo(dataset,rt.RooFit.Save())
    NS = ws.var("Nsig").getVal()
    ws.var("Nsig").setVal(0.)
    pdf.plotOn(frame,rt.RooFit.LineStyle(rt.kDashed))
    ws.var("Nsig").setVal(NS)
    pdf.plotOn(frame)
    likS = fitResult.minNll()

    # fit the data with B
    ws.var("Nsig").setVal(0.)
    ws.var("Nsig").setConstant(rt.kTRUE)
    ws.var("m0s").setConstant(rt.kTRUE)
    ws.var("s").setConstant(rt.kTRUE)
    fitResult2 = pdf.fitTo(dataset,rt.RooFit.Save())
    likB = fitResult2.minNll()

    print "Signal significance %f for R^2 > %s" %(rt.TMath.Sqrt(2.*rt.TMath.fabs(likS-likB)),RsqCut)
    
    # draw result
    c1 =rt.TCanvas("c1,"",600,600")
    #c1.SetLogy()    
    frame.Draw()
    c1.SaveAs("pippo.jpg")
    #c1.SaveAs("pippo.C")

    
    
