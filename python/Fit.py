from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
import sys
from WriteDataCard import *

def initializeWorkspace(w,cfg):
    variables = cfg.getVariablesRange(box,"variables",w)
    parameters = cfg.getVariables(box, "parameters")
    paramNames = []
    for parameter in parameters:
        w.factory(parameter)
        paramName = parameter.split('[')[0]
        if paramName.find("Cut")==-1 and paramName.find("Ntot")==-1:
            paramNames.append(paramName)
            w.var(paramName).setConstant(False)
        else:
            if paramName.find("Ntot")==-1:
                w.var(paramName).setConstant(True)
            else:
                w.var(paramName).setConstant(False)
                
        # turn off shape parameters if no events in b-tag bin:
        for i in [0, 1, 2, 3]:
            if "Ntot_TTj%ib"%i in paramName:
                print "help me!"
                w.var(paramName).setVal(w.data("RMRTree").sumEntries("nBtag==%i"%i) )
                if not w.var(paramName).getVal():
                    fixPars(w,"TTj%ib"%i)
                print paramName, w.var(paramName).getVal()
                    
    pdfs = cfg.getPdfs(box,"pdfs",w)
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store fit results")
    parser.add_option('--fit-region',dest="fitRegion",default="Full",type="string",
                  help="Fit region")
    parser.add_option('--no-fit',dest="noFit",default=False,action='store_true',
                  help="Turn off fit (useful for visualizing initial parameters)")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)
    box = options.box


    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            workspace = rootFile.Get("w"+box)
            data = workspace.data('RMRTree')
            
    w = rt.RooWorkspace("w"+box)
    rootTools.Utils.importToWS(w,data)
    initializeWorkspace(w,cfg)
    
    
    w.Print('v')
  

    pdf = w.pdf('extRazorPdf')
    
    #force numeical integrals and set precision
    #pdf.forceNumInt(True)
    #rt.RooAbsReal.defaultIntegratorConfig().Print("v")
    rt.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-13) 
    rt.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-13) 
    
    if options.noFit:
        fitResult = rt.RooFitResult()
    else:
        rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.FATAL)
        if options.fitRegion == "Full":
            nll = pdf.createNLL(data)
        else:
            nll = pdf.createNLL(data,rt.RooFit.Range(options.fitRegion))
        m = rt.RooMinuit(nll)
        m.migrad()
        m.hesse()
        fitResult = m.save()
        rt.RooMsgService.instance().reset()

        fitResult.Print('v')
        fitResult.correlationMatrix().Print('v')
    
        rootTools.Utils.importToWS(w,fitResult)
    
    mr = w.var('MR')
    rsq = w.var('Rsq')
    nbtag = w.var('nBtag')

    btagMin = nbtag.getMin()
    btagMax = nbtag.getMax()

    btag = ""
    if btagMax>btagMin+1:
        btag = "%i-%ibtag"%(btagMin,btagMax-1)
    else:
        btag = "%ibtag"%(btagMin)

    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    
    c = rt.TCanvas("c","c",600,400)
    c.SetLogy()
    mrFrame = mr.frame(x[0],x[-1],50)
    #mrFrame = mr.frame(x[0],2500,50)
    mrFrame.SetTitle("")
    mrFrame.SetXTitle("M_{R}")
    rsqFrame = rsq.frame(y[0],y[-1],50)
    #rsqFrame = rsq.frame(y[0],1.,50)
    rsqFrame.SetTitle("")
    rsqFrame.SetXTitle("R^{2}")

    
    def plot1d(data,pdf,var,frame,btag,c):
        data.plotOn(frame,rt.RooFit.Name("Data"),rt.RooFit.Invisible())
        #pdf.plotOn(frame,rt.RooFit.VisualizeError(fitResult,0.25),rt.RooFit.FillColor(rt.kBlue-10),rt.RooFit.Range("Full"),rt.RooFit.NormRange("Full"))
        pdf.plotOn(frame,rt.RooFit.Name("Total"),rt.RooFit.FillColor(rt.kBlue-10),rt.RooFit.Range("Full"),rt.RooFit.NormRange("Full"),rt.RooFit.Normalization(pdf.expectedEvents(w.set('variables')),rt.RooAbsReal.NumEvent))
        pdf.plotOn(frame,rt.RooFit.Name("TTj1b"),rt.RooFit.Components('razor3dPdf_TTj1b'),rt.RooFit.LineColor(rt.kViolet),rt.RooFit.LineStyle(rt.kDashed),rt.RooFit.Range("Full"),rt.RooFit.NormRange("Full"),rt.RooFit.Normalization(w.var('Ntot_TTj1b').getVal(),rt.RooAbsReal.NumEvent))
        pdf.plotOn(frame,rt.RooFit.Name("TTj2b"),rt.RooFit.Components('razor3dPdf_TTj2b'),rt.RooFit.LineColor(rt.kRed),rt.RooFit.LineStyle(rt.kDashed),rt.RooFit.Range("Full"),rt.RooFit.NormRange("Full"),rt.RooFit.Normalization(w.var('Ntot_TTj2b').getVal(),rt.RooAbsReal.NumEvent))
        pdf.plotOn(frame,rt.RooFit.Name("TTj3b"),rt.RooFit.Components('razor3dPdf_TTj3b'),rt.RooFit.LineColor(rt.kGreen),rt.RooFit.LineStyle(rt.kDashed),rt.RooFit.Range("Full"),rt.RooFit.NormRange("Full"),rt.RooFit.Normalization(w.var('Ntot_TTj3b').getVal(),rt.RooAbsReal.NumEvent))
        data.plotOn(frame,rt.RooFit.Name("Data"))
        frame.SetMinimum(0.01)
        frame.SetMaximum(data.sumEntries()*2)
        frame.Draw()

        l = rt.TLatex()
        l.SetTextAlign(11)
        l.SetTextSize(0.05)
        l.SetTextFont(42)
        l.SetNDC()
        l.DrawLatex(0.72,0.92,"3 fb^{-1} (13 TeV)")
        l.DrawLatex(0.15,0.85,"CMS simulation")
        l.SetTextFont(52)
        if options.fitRegion=="Full":
            fitRegion = "fit"
        else:
            fitRegion = "fit"
        l.DrawLatex(0.15,0.80,"razor %s %s"%(box,fitRegion))
        leg = rt.TLegend(0.7,0.59,0.89,0.88)
        leg.SetTextFont(42)
        leg.SetFillColor(rt.kWhite)
        leg.SetLineColor(rt.kWhite)
        leg.AddEntry(frame.findObject("Data"),"Sim Data","pe")
        leg.AddEntry(frame.findObject("Total"),"Total","lf")
        leg.AddEntry(frame.findObject("TTj1b"),"1b-tag","l")
        leg.AddEntry(frame.findObject("TTj2b"),"2b-tag","l")
        leg.AddEntry(frame.findObject("TTj3b"),"3b-tag","l")
        leg.Draw()
    
        c.Print(options.outDir+"/RooPlot_"+var.GetName()+"_"+options.fitRegion.replace(',','_')+"_"+btag+"_"+box+".pdf")
        c.Print(options.outDir+"/RooPlot_"+var.GetName()+"_"+options.fitRegion.replace(',','_')+"_"+btag+"_"+box+".C")

    plot1d(data,pdf,mr,mrFrame,btag,c)
    plot1d(data,pdf,rsq,rsqFrame,btag,c)

    inFiles = [f for f in args if f.lower().endswith('.root')]
            
    if len(inFiles)==1:
        outFile = inFiles[0].split('/')[-1].replace('RazorAnalysis','FitResult').replace(box,options.fitRegion.replace(',','_')+"_"+box)

    outFile = rt.TFile.Open(options.outDir+"/"+outFile,'recreate')
    outFile.cd()
    w.Write()
    outFile.Close()
