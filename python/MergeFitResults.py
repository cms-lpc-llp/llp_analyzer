#!/usr/bin/env python
import ROOT as rt
import rootTools
from optparse import OptionParser

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option('-o','--output',dest="output",type="string",default='BinnedFitResults.root',
                  help="Name of the output file to use",metavar='FILE') 
    (options,args) = parser.parse_args()

    def getBoxContents(box, store):
        #rt.gDirectory.cd(box)
        keys = rt.gDirectory.GetListOfKeys()
        for k in rootTools.RootIterator.RootIterator(keys):
            name = k.GetName()
            o = rt.gDirectory.Get(name)
            if o is None or not o: continue
            store.add(o,name=name)
        rt.gDirectory.cd()
    
    store = rootTools.RootFile.RootFile(options.output)

    rt.gSystem.Load("python/lib/libRazorRun2")
    for f in args:
        if '.root' in f:
            box = f.split('/')[-1].split('.root')[0].split('_')[-1]
            print box            
            tfile = rt.TFile.Open(f)
            getBoxContents(box,store)
    
    store.write()
