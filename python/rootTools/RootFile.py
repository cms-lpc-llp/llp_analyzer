#$Revision:$

import ROOT as rt
import os.path

class RootFile(object):
    
    def __init__(self, fileName):
        self.fileName = fileName
        self.plots = {}
        
    def add(self, plot, name = None, dir = None):
        if name is None: name = plot.GetName()
        if dir is not None: name = os.path.join(dir,name)

        l = self.plots.get(name,[])
        l.append(plot)
        self.plots[name] = l

    def write(self):
        out = None
        try:
            out = rt.TFile.Open(self.fileName,'RECREATE')
            for name, plots in self.plots.iteritems():
                out.cd()
                
                dir = os.path.dirname(name)
                objname = os.path.basename(name)
                
                #support for directories if required
                if dir:
                    d = out.Get(dir)
                    if d is not None and d and d.InheritsFrom('TDirectory'):
                        pass
                    else:
                        out.mkdir(dir)
                    out.cd(dir)
                
                if not plots:
                    continue
                elif len(plots) == 1:
                    plots[0].Write(objname)
                else:
                    index = 0
                    for i in xrange(len(plots)):
                        p = plots[i]
                        p.Write('%s_%i' % (objname,i))
            #needed so that the object can be deleted 
            self.plots.clear()
        finally:
            if out is not None: out.Close()

#------------------------------------------------------------------------------
# File: Ntuple.py
# Description: a simple ntuple reader
# Created: 22 Sep 2010 Harrison B. Prosper
#          24 Aug 2012 HBP - handle ntuples with lots of variables
#------------------------------------------------------------------------------
import os, sys
from string import *
from ROOT import *	 
#------------------------------------------------------------------------------
class Buffer:

	def __init__(self, buffer=None, buffermap=None, variable=None):
		self.buffermap= buffermap
		self.variable = variable
		self.buffer   = buffer

	def __getattr__(self, variable):
		if self.buffermap.has_key(variable):
			jj = self.buffermap[variable]
			return self.buffer[jj].__getattribute__(variable)
		else:
			raise AttributeError(variable)
		
	def __str__(self):
		s = 'Event variables:\n'
		for tname, name, maxcount in self.variable:
			s += "  %-12s %-24s:\t" % (tname, name)
			s += "%s\n" % self.__getattr__(name)
		return s

	def clone(self):
		ev = Buffer()
		ev.buffermap = self.buffermap
		ev.variable  = self.variable

		# copy data
		
		ev.buffer = []
		for ii in xrange(len(self.buffer)):
			bufferName = "Buffer%d" % ii
			# import struct
			exec("from ROOT import %s" % bufferName)
			ev.buffer.append(eval("%s()" % bufferName))

		for tname, variable, maxcount in self.variable:
			jj = self.buffermap[variable]
			cmd = 'ev.buffer[jj].%s = self.%s' % (variable, variable)
			exec(cmd)
		return ev
#------------------------------------------------------------------------------
class Ntuple:
	# "self" is Python's equivalent of the "this" pointer in C++
	# self points to the memory allocated for the object
	
	def __init__(self, filename, treename):

		# cache inputs
		
		self.filename = filename # Root file name
		self.treename = treename
		
		# check that root file exists
		if not os.path.exists(filename):
			print "** Ntuple *** "\
			      "root file %s not found" % filename
			sys.exit(0)

		# ------------------------------------
		# open root file
		# ------------------------------------
		print "reading root file\t--> %s" % filename
		
		self.file = TFile(self.filename)
		self.tree = self.file.Get(self.treename)
		tree = self.tree
		
		if tree == None:
			print "** Ntuple *** Tree is non-cooperative"\
				  " - perhaps the name is wrong?"
			sys.exit(0)

		# get names of variables from root file
		branches  = tree.GetListOfBranches()
		# get number of variables 
		nbranches = branches.GetEntries()

		# print variables
		self.vars = []
		print "variables:"
		for i in xrange(nbranches):
			# get the ith branch (aka variable)
			bname = branches[i].GetName()			
			
			# assume that the leaf attached to this branch
			# has the same name as the branch
			leafname = bname
			leaf  = branches[i].GetLeaf(leafname)
			if leaf == None:
				# Get all leaves associated with this branch
				leaves = branches[i].GetListOfLeaves()
				if leaves == None:
					print "No leaves found!"
					sys.exit(0)

				# Assume one leaf/branch
				leaf = leaves[0]
				if leaf == None:
					print "No leaf found"
					sys.exit(0)
					
				leafname = leaf.GetName()
				
			# get leaf type (int, float, double, etc.)
			tname = leaf.GetTypeName()

			#check for leaf counter
			flag = Long(0)
			leafcounter = leaf.GetLeafCounter(flag)
			if leafcounter:
				maxcount = leafcounter.GetMaximum()
			else:
				maxcount = leaf.GetLen()
				
			# store type and variable name
			self.vars.append( (tname, bname, maxcount) )
			print "\t%4d\t%-12s\t%-32s\t%d" % (i, tname, bname, maxcount)
		nlen = len(self.vars)

		# create a map of variable name to column number
		self.varmap = {}
		for ind, var in enumerate(self.vars):
			self.varmap[var] = ind

		# initialize row number
		self.row = 0
		self.entries = self.tree.GetEntries()
		# ------------------------------------
		# set up branches as a struct
		# Root has a limit on how long a
		# string it can cope with, so split
		# into multiple strings
		# ------------------------------------

		bufferCount = 0
		newBuffer = True
		rec = ""
		bufferName = ""
		maxlength  = 2000
		self.buffermap  = {}
		self.buffer = []
		self.buffernames = []
		
		for count, (tname, name, maxcount) in enumerate(self.vars):

			# keep track of map from variable name to buffer count
			self.buffermap[name] = bufferCount
			
			if newBuffer:
				newBuffer = False
				bufferName = "Buffer%d" % bufferCount
				self.buffernames.append(bufferName)
				
				rec = "struct %s {" % bufferName

			if maxcount == 1:
				rec += "%s %s;" % (tname, name)
			else:				
				rec += "%s %s[%d];" % (tname, name, maxcount)

			if (len(rec) > maxlength) or \
				   (count >= len(self.vars)-1):
				rec += "};"
				newBuffer = True
				
				# evaluate it
				gROOT.ProcessLine(rec)

				# now import struct
				exec("from ROOT import %s" % bufferName)

				# add to list of buffers
				self.buffer.append(eval("%s()" % bufferName))

				# remember of update buffer count
				bufferCount += 1

		# create a generic reusable event object
		self.event = Buffer(self.buffer, self.buffermap, self.vars)
		
		# Now that addresses are stable, give address of each variable
		for tname, name, maxcount in self.vars:
			jj = self.buffermap[name]
			tree.SetBranchAddress(name, AddressOf(self.buffer[jj], name))

	# destructor
	def __del__(self):
		pass

	def entries(self):
		return self.entries

	def read(self, row):
		self.tree.GetEntry(row)

	def get(self, variable):
		if buffermap.has_key(variable):
			jj = buffermap[variable]
			return self.buffer[jj].__getattribute__(variable)
		else:
			return None
		
	# Implement Python iterator protocol
	
	def __iter__(self):
		return self

	def next(self):
		if self.row > self.entries-1:
			self.row = 0
			raise StopIteration
		else:
			self.read(self.row)
			self.row += 1
			return self.event

#------------------------------------------------------------------------------
# Pick event with probability proportional to its weight.
# Algorithm:
# Represent an event by an interval equal to its weight, pick a number at
# random between 0 and the summed weights and find in which interval (that is,
# event) the number falls and pick that event.
#
# Example:
# 
# Consider the 5 events below.Event 1 has weight 5, event 2 weight 3, event 3
# weight 3, event 4 weight 1, and event 5 has weight 4. The total length of
# the intervals, 16, represents the summed weights. Pick a number at random
# between 0 and 16. Suppose it is 6.3. This lands in event 2 so it is selected.
#
# 1     2   3   4 5
# *****|***|***|*|****
#
# In practice, we find the interval in which the random number lands using a
# modified version of a binary search. See code at
# http://personal.denison.edu/~havill/102/python/search.py
#------------------------------------------------------------------------------
def pickEventAtRandom(wcdf):
	from random import uniform

	# choose a value at random between 0 and total weight
	# and find in which "event" it lands
	
	w = uniform(0, wcdf[-1])
	first = 0
	last  = len(wcdf) - 1   
	found = False

	while (first <= last) and not found:
		mid = (first + last) / 2
		if w <= wcdf[mid]:
			last = mid       
		else:
			first = mid + 1

		if first >= last:
			mid = first
			found = True

		if found: return mid
	return -1
