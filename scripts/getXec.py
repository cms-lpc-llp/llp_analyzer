import os
import sys

#if not os.path.isfile('rest.py'):
#    print "getting the mcm rest api"
#    os.system('ln -s /afs/cern.ch/cms/PPD/PdmV/tools/McM/rest.py')

from rest import restful

McM = restful(dev=False)

dataset=sys.argv[1]

rs = McM.getA('requests',query='produce=%s'% dataset)
if len(rs)>1:
    print "this cannot really be"
elif len(rs)==0:
    print dataset,"is not produced via mcm"
else:
    print dataset,"produced by",rs[0]['prepid']
    r = rs[0]
    ## pull out the chains
    crs = McM.getA('chained_requests',query='contains=%s'% r['prepid'])
    if len(crs)>1:
        print "unlikely to have more than one chain for an analysis dataset"
    elif len(crs)==0:
        print r['prepid'],"is not in any chain"
    else:
        infos=[]
        for r in reversed(crs[0]['chain']):
            rr = McM.getA('requests',r)
            if len(rr['generator_parameters']):
                gp = rr['generator_parameters'][-1]
                infos.append( (gp['cross_section'], gp['filter_efficiency'], gp['match_efficiency'] ) )
        print ["%s pb"%i[0] for i in infos]
        print ["%s filtering"%i[1] for i in infos]
        print ["%s matching"%i[2] for i in infos]
        
                             
    
