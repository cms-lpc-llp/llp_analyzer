#! /usr/bin/env python

# V1.3
# Changes:
#   o add --begin and --end functionality 

# V1.2
# Changes:
#   o remove spaces from lumibyls option

# Known problems/missing functionality:
#   o minimal checking of input errors
#   o no lumibyXing option

import argparse
import socket
import time

def getLSLumi(args,vals,nLS,deliveredLumi,recordedLumi) :
    beamMode = vals[fields.index('Mode')].lstrip().rstrip()  
    if (args.beamMode == 'All' or args.beamMode == beamMode) and testInterval(vals) :
        nLS += 1
        lumi = float(vals[fields.index(args.source)])
        deliveredLumi += LStime*lumi
        recordedLumi += LStime*lumi*(1. - float(vals[fields.index('deadtime')]))
    else :
        lumi = 0.

    if args.lumibyls :
        #Run:Fill,LS,UTCTime,Beam Status,E(GeV),Delivered(/ub),Recorded(/ub),avgPU
        #149382:1455,150:150,12/30/15 03:33:59,STABLE BEAMS,6500.1,11.383,11.383,0.752
        delLumi = LStime*lumi
        recLumi = LStime*lumi*(1. - float(vals[fields.index('deadtime')]))
        outString = "{0:s}:{1:s},{2:s}:{2:s},{3:s},{4:s},6500.,{5:.5f},{6:.5f},{7:.5f}".format(
            vals[fields.index('Run')].strip(),vals[fields.index('Fill')].strip(),vals[fields.index('LS')].strip(),
            vals[fields.index('time')].strip(),vals[fields.index('Mode')].strip(),delLumi,recLumi,getPU(vals))
        if args.printLevel > 0 : print outString 
        resultOutput(args,outString) 
            
    return nLS, lumi, deliveredLumi, recordedLumi

def getUnit(integratedLumi) :
    if integratedLumi < 1000.   : return '/ub', 1.
    elif integratedLumi < 1.0e6 : return '/nb', .001
    elif integratedLumi < 1.0e9 : return '/pb', 1.0e-6
    else : return '/fb', 1.0e-9

def getPU(vals) :
    sigmaMB = 78400.
    lumi = float(vals[fields.index('lumi')])
    nBX = float(vals[fields.index('Ncol')])
    fLHC = 11246.
    return sigmaMB*lumi/nBX/fLHC

def getArgs() :
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fill", help="Fill number (or 'all')")
    parser.add_argument("-r","--run", type=int, help="Run number.")
    parser.add_argument("-i","--JSON",help="JSON file")
    parser.add_argument("-o","--outputFile",help="Output file")
    parser.add_argument("-l","--lumibyls",action="store_true",help="Lumi by LS mode")
    parser.add_argument("-b","--beamMode",default="STABLE BEAMS",help="Beam mode (default 'STABLE BEAMS').")
    parser.add_argument("-p","--printLevel", type=int, help="Print level")
    defaultDataDir = "./runcsv/"
    if 'lxplus' in socket.gethostname() : defaultDataDir = "/afs/cern.ch/user/m/marlow/public/lcr2/runcsv/"
    parser.add_argument("-d","--dataDirectory", default=defaultDataDir,help="path to runcsv directory")
    parser.add_argument("-s","--source", default='lumi', help="Lumi source (default = best)")
    parser.add_argument("--begin",help="Begin fill, run, or time in mm/dd/yy hh:mm:ss format")
    parser.add_argument("--end",help="End fill, run, or time in mm/dd/yy hh:mm:ss format")
    return parser.parse_args()

def resultOutput(args,str) :
    global started
    global fOut 
    if not started :
        started = True
        if args.outputFile :
            fOut = open(args.outputFile,'w')
            fOut.write("Run:Fill,LS,UTCTime,Beam Status,E(GeV),Delivered(/ub),Recorded(/ub),avgPU\r\n")
        else :
            print "  Run:Fill,   LS,         UTC,            Status,  E(GeV), Del(/ub), Rec(/ub), avgPU"
                
    if args.outputFile :
        fOut.write(str + '\n')
    else :
        print str

def setBeginAndEnd(args) :
    global firstFill, lastFill, firstRun, lastRun, startTime, endTime 
    firstFill, lastFill, firstRun, lastRun, startTime, endTime = 0, 9999, 0, 999999, 0., time.time()

    if args.begin :
        if len(args.begin) == 4 and args.begin.isdigit():
            firstFill = int(args.begin)
        elif len(args.begin) == 6 and args.begin.isdigit():
            firstRun = int(args.begin)
        else :
            sTime = time.strptime(args.begin, "%m/%d/%y %H:%M:%S")
            UTC = time.mktime(sTime)
            startTime = UTC

    if args.end :
        if len(args.end) == 4 and args.end.isdigit():
            lastFill = int(args.end)
        elif len(args.end) == 6 and args.end.isdigit():
            lastRun = int(args.end)
        else :
            sTime = time.strptime(args.end, "%m/%d/%y %H:%M:%S")
            UTC = time.mktime(sTime)
            endTime = UTC

def testInterval(vals) :
    if firstFill > 0 or lastFill < 9999 :
        fill = int(vals[fields.index('Fill')])
        if fill < firstFill or fill > lastFill : return False
    if firstRun > 0 or lastRun < 999999 :
        run = int(vals[fields.index('Run')])
        if run < firstRun or run > lastRun : return False
    if startTime > 0.  or endTime < time.time() :
        sTime = time.strptime(vals[fields.index('time')].strip(), "%m/%d/%y %H:%M:%S")
        utc = time.mktime(sTime)
        if utc < startTime or utc > endTime : return False
    return True
    
        
# execution starts here

LStime = 23.31 
fields = ['Fill', 'Run', 'LS', 'Mode', 'source', 'time', 'lumi','deadtime', 'Ncol', 'HF', 'PLT', 'BCMF', 'I', 'L[I]']
nLS, deliveredLumi , recordedLumi = 0, 0., 0.
args = getArgs()
setBeginAndEnd(args)

started = False

if args.JSON :
    runList = eval(open(args.JSON, 'r').read())

    for run in runList :
        if args.printLevel > 0 : print "Opening: {0:s}run{1:s}.csv".format(args.dataDirectory,run) 
        lines = open("{0:s}run{1:s}.csv".format(args.dataDirectory,run),'r').readlines()
        for line in lines :
            vals = line.split(',')
            for ls_interval in runList[run]:
                if ls_interval[0] <= int(vals[fields.index('LS')]) <= ls_interval[1]:
                    nLS, lumi, deliveredLumi, recordedLumi = getLSLumi(args,vals,nLS,deliveredLumi,recordedLumi) 

    unit, factor = getUnit(deliveredLumi)
    print "JSON mode: \n runList={0:s} \nnLS={1:d} Integrated Lumi: delivered={2:8.3f} ({3:s}) recorded={4:8.3f} ({3:s})".format(
        runList,nLS,factor*deliveredLumi,unit,factor*recordedLumi)

if args.run :
    runFile = "{0:s}run{1:d}.csv".format(args.dataDirectory,args.run)
    if args.printLevel > 0 : print "Opening: {0:s}".format(runFile) 
    lines = open(runFile,'r').readlines()
    for line in lines : 
        vals = line.split(',')
        nLS, lumi, deliveredLumi, recordedLumi = getLSLumi(args,vals,nLS,deliveredLumi,recordedLumi) 

    unit, factor = getUnit(deliveredLumi)
    print "Run mode: Run = {0:d}, nLS={1:d} Integrated Lumi: delivered={2:8.3f} ({3:s}) recorded={4:8.3f} ({3:s})".format(
        args.run,nLS,factor*deliveredLumi,unit,factor*recordedLumi)

if args.fill :
    if args.fill.lower() == 'all' : targetFill = 0
    else : targetFill = int(args.fill) 
    import glob
    runFileList = glob.glob(args.dataDirectory + '*')
    nRun = len(runFileList) - 1
    for runFile in runFileList :
        if args.printLevel > 0 : print "Opening: {0:s}".format(runFile) 
        lines = open(runFile,'r').readlines()
        nextRun = False 
        for line in lines :
            vals = line.split(',')
            if targetFill == 0 or targetFill == int(vals[fields.index('Fill')])  :
                nLS, lumi, deliveredLumi, recordedLumi = getLSLumi(args,vals,nLS,deliveredLumi,recordedLumi) 
            else : # if the fill in this file isn't right, move on to the next file
                nRun -= 1
                nextRun = True
                break
            
        if nextRun : continue

    if nRun > 0 :
        unit, factor = getUnit(deliveredLumi)
        print "Fill mode: Fill = {0:s}, #Runs={1:d} #LS={2:d} Integrated Lumi: delivered={3:8.3f} ({4:s}) recorded={5:8.3f} ({4:s})".format(
            args.fill,nRun,nLS,factor*deliveredLumi,unit,factor*recordedLumi)
    else:
        print "Fill mode: Fill = {0:s} no runs found.".format(args.fill) 
        
if not (args.fill or args.run or args.JSON) :
    print "Must specify fill, run, or JSON file.  For help enter lcr2.py -h"


    


