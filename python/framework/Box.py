import ROOT as rt
import rootTools
from framework import Config

class Box(object):
    
    def __init__(self, name, cfg, workspace = None):
        self.name = name
        self.cfg = cfg
        
        if workspace is None:
            self.workspace = rt.RooWorkspace("w"+name)
        else:
            self.workspace = workspace
            
        self.fitmodel = 'extRazorPdf'
        

    def fixPars(self, label, doFix=True, setVal=None):
        parSet = self.workspace.allVars()
        for par in RootTools.RootIterator.RootIterator(parSet):
            if label in par.GetName():
                par.setConstant(doFix)
                if setVal is not None: par.setVal(setVal)
