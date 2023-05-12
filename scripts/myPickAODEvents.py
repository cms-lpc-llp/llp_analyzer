#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys



if len(sys.argv) < 2:
    print "Usage: python myPickEvents.py [inputfile] [datasetName]\n"
    exit()

inputfileName = sys.argv[1]
datasetName = sys.argv[2]

#######################################################
#Run DAS query to get full map of lumi -> file
#######################################################
lumiToFileMap = dict()
if os.path.exists("tmpDASQueryResult.txt"):
    os.system("rm -f tmpDASQueryResult.txt")
#command = "dasgoclient -query=\"file,run,lumi dataset="+datasetName+"\" > tmpDASQueryResult.txt"
#print(command)

command = "dasgoclient -query=\"file,run,lumi dataset="+datasetName+" instance=prod/phys03\" > tmpDASQueryResult.txt"
print command
os.system(command)
tempfile = open("tmpDASQueryResult.txt","r")
templines = tempfile.readlines()
for line in templines:
    tmpSplitLine = line.split()
    tmpFilename = tmpSplitLine[0]
    tmpLumiList = tmpSplitLine[2].strip("[]").split(",")
    for lumi in tmpLumiList:
        if not lumiToFileMap.has_key(lumi):
            lumiToFileMap[lumi] = tmpFilename
        else:
            print "Warning: lumi "+lumi+" appeared twice"

#print lumiToFileMap

#######################################################
#Loop over requested lumi and event numbers and produce
#list of edmCopyPickMerge commands
#######################################################
inputfile = open(inputfileName,"r")
eventList = inputfile.readlines()
fileToEventListMap = dict()

for event in eventList:
    #print "Pick Event for event: "+event
    eventStringList = event.strip("\n").split(":")
    lumiNum = eventStringList[1]
    eventNum = eventStringList[2]
    filename = ""
    if (lumiToFileMap.has_key(lumiNum)):
        filename = lumiToFileMap[lumiNum]
        if (not fileToEventListMap.has_key(filename)):
            fileToEventListMap[filename] = list()
        fileToEventListMap[filename].append("1:"+lumiNum+":"+eventNum)

#print fileToEventListMap
            
keycount = 0
for key in fileToEventListMap:
    outputCommand = "edmCopyPickMerge outputFile=pickevents_"+str(keycount)+".root eventsToProcess="

    eventListString = ""
    for event in fileToEventListMap[key]:
        eventListString = eventListString + event + ","
    outputCommand = outputCommand + eventListString.strip(",") + " inputFiles=" + key
    print outputCommand
    keycount = keycount + 1
    
exit()
