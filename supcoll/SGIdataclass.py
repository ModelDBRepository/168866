#### SGI data class for analysis ####
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2009-2013 Jan Moren
#
# This file is mostly provided for the file access functions. The analysis
# accessory functions should work but have not been tested with recent
# versions of the model. Take them as inspiration for your own analysis tools
# rather than as finished, reliable components.



import math
import numpy
import json
import time
import glob
import os
import sys
import bz2
import pickle


from copy import deepcopy
import scipy

if scipy.__version__>="0.7.0":
    HAVEKDTree=True
    from scipy.spatial import KDTree
else:
    HAVEKDTree=False



## Support functions

# Give this a simulation file name
#
# Figure out the path, the parameter file and ending ("sim" or "zim") and the
# path and ending ("spikes" or "zpikes") to the individual spike data files.
# Prioritize the compressed versions.

def find_sim(fname):
    fhead, ftail = os.path.split(fname)
    fbase,ext = os.path.splitext(ftail)

    if ext == ".sim":
	pass

    elif ext == ".zim":
	pass

    else:		# We only got a basename, so look for the 
			# right sim file
	fbase = ftail
	if os.path.isfile(os.path.join(fhead,fbase+".zim")):
	    ext=".zim"
	else:
	    ext=".sim"
    
    return (fhead, fbase, ext)


# Load a simulation parameter file

def load_simfile(spath):
    
    pathsimfile = os.path.join(spath[0], spath[1]+spath[2])
    try:
	if spath[2] == ".sim":
	    f = open(pathsimfile, "r")
	elif spath[2] == ".zim":
	    f = bz2.BZ2File(pathsimfile, "r")
	else:
	    print "Unsupported file ending: ", spath[2]
	    sys.exit(-1)
	try:
	    simdata = json.load(f)

	finally:
	    f.close()

    except IOError:
	
	print "Failed to load simulation data from", pathsimfile
	sys.exit(-1)

    return simdata


# load a spike file

def load_spikefile(spath, name):
    
    spiked = {} 
    if spath[2] == ".sim":
        for infile in glob.glob(os.path.join(spath[0], spath[1], name+'*.spikes')):
	    try: 
		for line in open(infile, "r"):
		    tmp = line.split()
		    gid = tmp[0]
		    t = tmp[1]
		    spiked.setdefault(int(gid),[]).append(float(t))
	    
	    except IOError:
		print "File open error reading "+name+".spikes\n"
		sys.exit(-1)
	
    elif spath[2] == ".zim":
	f=0
	try: 
	    f = bz2.BZ2File(os.path.join(spath[0], spath[1], name+".zpikes"),
		"r")
	    spiked = pickle.load(f)
	    f.close()
	except IOError:
	    print "File open error reading "+name+".zpikes\n"
	    sys.exit(-1)

    else:
	print "file format must be given as .sim or .zim"
	sys.exit(-1)

    return spiked

# load a blob file

def load_blobfile(spath, name):
    
    spiked = [] 
    if spath[2] == ".sim":
        for infile in glob.glob(os.path.join(spath[0], spath[1], name+'*.spikes')):
	    try: 
		for line in open(infile, "r"):
		    spiked.append(float(line))
	    
	    except IOError:
		print "File open error reading "+name+".spikes\n"
		sys.exit(-1)
	
    elif spath[2] == ".zim":
	f=0
	try: 
	    f = bz2.BZ2File(os.path.join(spath[0], spath[1], name+".zpikes"),
		"r")
	    spiked = pickle.load(f)
	    f.close()
	except IOError:
	    print "File open error reading "+name+".zpikes\n"
	    sys.exit(-1)

    else:
	print "file format must be given as .sim or .zim"
	sys.exit(-1)

    return spiked



##############################
##
## Class for a neuron layer 
##

class Surface:
    
    def __init__(self, surface, spath, runs = 1.0):

	self.rows = surface["rows"] 
	self.cols = surface["cols"] 
	self.pitch = surface["pitch"] 
	self.nname = surface["name"] 
	self.xcenter = surface["xcenter"] 
	self.ycenter = surface["ycenter"] 
	self.runs = runs

	self.xdim	= (self.cols)*self.pitch
	self.ydim	= (self.rows)*self.pitch
	self.dirname	= spath[0]
	self.lists	= {}
	inset		= []

	# create a container for our data
	self.lists['orig'] = SpikeTrain(surface)

	spiked = load_spikefile(spath, self.nname)
	
	self.lists['orig'].add_set(spiked)
	self.lists['orig'].proc_set()


    # Create new data subset covering only the neurons in gids

    def new_dataset(self, name, gids=[]):

	if name == "orig":
	    print "'orig' is a protected name for data subsets."
	    exit()

	self.lists[name] = deepcopy(self.lists['orig'])
	if gids:
	   self.lists[name].subset(gids)


    def make_subset(self, gids, dataset="orig"):

	   self.lists[dataset].subset(gids)


	
    #######################################################################
    #
    # Accessory methods for the analysis
    #

    # Calculate the average interspike interval and spikes for firing neurons
    # in a given set of events. 
    #
    # small wrinkle: the very first spike has no preceding event, so a neuron
    # with only the very first spike and no more events does not get recorded
    # as having a positive ISI, while still getting counted for the simple
    # summation. 
    #
    # Another wrinkle: we assume that only events actually happening during
    # the time period will count. A neuron that fired before the start and
    # after the end won't get added here even though it technically has a
    # well-defined ISI during the period.
    #
    # Final issue: When we analyze collated runs, a neuron may have multiple
    # spikes at the same instant. This can break the ISI estimation in some
    # cases. It's good enough for plotting only.
    #
    # return (gid, isi rate, #spikes) tuple

    def get_event_rate(self, events):

	tlist={}
	results =[]

	for ev in events: 
	    gid = ev[1]
	    tdiff = ev[2]

	    tmp = tlist.setdefault(gid,(0.0, 0))	# if not set, {gid} 
							# gets default tuple
	    tlist[gid] = (tmp[0] + tdiff, tmp[1] + 1)

	for gid, t in tlist.iteritems():
	    if t[0] > 0.0:
		results.append((gid, 1000.0*t[1]/t[0], t[1]))
	    else:
		results.append((gid, 0.0, t[1]))

	return results



    # Find gid of neuron closest to specific point

    def neuron_at(self, x, y, dataset="orig"):
	
	a = self.lists[dataset].xy_gid(x,y)
	return a

    def neurons_around(self, x, y, r, dataset="orig"):
	
	return self.lists[dataset].xyr_gids(x, y, r)

    def neuron_at_mm(self, mmx, mmy, dataset="orig"):
	
	return self.lists[dataset].mm_gid(mmx,mmy)

    def neurons_around_mm(self, mmx, mmy, mmr, dataset="orig"):
	
	return self.lists[dataset].mmr_gids(mmx, mmy, mmr, dataset)





    #######################################################
    # slice_rate
    #
    # rate of a single neuron or group of neurons by averaged ISI over a timeslice
    #
    # tstart, tend: start and end time 
    # 
    # tperiod:	step size
    #
    # binsize: bin size for averaging
    #
    # avgtype: "bin" for spike binning, "isi" for interstimulus interval 
    #
    # neurons: "True" to return the individual neuron data
    #
    # avruns: "True" to treat the data as collated from "runs" number of 
    # simulations
    #
    # return a tuple (max value, max time, cumulative value, list of (time, rate)
    # tuples)


    def slice_rate(self, tstart, tend, tperiod, binsize = 0, avgtype="bin",
	    neurons=False, avruns=False, vcut = -1, dataset="orig"):

	i = tstart*1.0
	val=[]
	sev=[]
	cum=0.0
	maxval=0.0
	maxtime=0.0
	neurlist = []
	totsp = 0
	cumsp = 0
	cumtime = -1
	r = 1.0

	# NOTE: doesn't work ATM - set runs manually in init.
#	if avruns:
#	    r = float(self.runs)
    
	#print r
	if binsize == 0:
	    binsize = tperiod

	if vcut>0.0:
	    ev = self.lists[dataset].timeslice(tstart*1.0, tend*1.0)
	    binlist = self.get_event_rate(ev)
	    for gid, isi, b in binlist:
		totsp += b


	while i<tend:
	
	    ev = self.lists[dataset].timeslice(i-binsize*0.5,
		    i+binsize*0.5)
	    
	    binlist = self.get_event_rate(ev)
	    

	    avg_isi = 0.0
	    avg_bin = 0.0
	    avg = 0.0

	    if neurons:
		neurlist.append((i,binlist))

	    for gid, isi, b in binlist:
		avg_bin = avg_bin + b
		avg_isi = avg_isi + isi

	    if vcut > 0.0 and cumtime <0.0:
		ev2 = self.lists[dataset].timeslice(i, i+tperiod+0.01)
		bl2 = self.get_event_rate(ev2)
		
		for gid, isi, b in bl2:
		    cumsp += b
#		print "i, cumsp, totspcut", i, cumsp , totsp
		if cumsp >= totsp*vcut:
		    cumtime = i
#		    print "cumsp, totsp, vcut:", cumsp, totsp, vcut
	    
	    if avgtype == "bin":
		avg = (1000.0*avg_bin/(binsize*r))

	    elif avgtype == "isi":
		avg = avg_isi/r
		
	    else: 
		print "slice_rate: %s is an unknown average type" % (avgtype)
		exit (-1)


	    if maxval < avg:
		maxval = avg
		maxtime = i

	    val.append((i, avg))
	    cum=cum + avg_bin
	    i = i + tperiod

	if neurons:
	    return (maxval, maxtime, cum/r, val, neurlist)
	    
	elif vcut>0.0:
	    return (maxval, maxtime, cum/r, val, cumtime)
	else:
	    return (maxval, maxtime, cum/r, val)

    ###################################################################
    # eyepos
    #
    # figure out the eye position change for a given range of spikes

    def eyepos(self, tstart, tend, totspikes, tperiod=1.0, A=3.0, bx=1.4,
	    by=1.8, avruns=True,  dataset="orig"):

	
	i = tstart*1.0
	xpos = 0.0
	ypos = 0.0
	poslist = []
	sev = self.lists[dataset] 
	while i<tend:
	    
	    if (i-tperiod)>=0.0:
		events = sev.timeslice(i-tperiod, i)
	    else:
		events = sev.timeslice(0.0, i)

	    cumx = 0.0
	    cumy = 0.0
	    for ev in events:

		gid = ev[1]
		(x,y) = sev.gid_xy(gid)
		(mx, my) =sev.xy_mm(x,y)

		wx = A*(math.exp(mx/bx)*math.cos(my/by)-1)
		wy = A*math.exp(mx/bx)*math.sin(my/by)
		
		cumx = cumx + wx
		cumy = cumy + wy

	    xpos = xpos + cumx/totspikes
	    ypos = ypos + cumy/totspikes
	
	    poslist.append((i,xpos,ypos))

	    i=i+tperiod

	return poslist



##############################
##
## Class for a neuron blob 
##
## No spatial positions or gids so much simpler
##

class Blob:
    
    def __init__(self, blob, spath):

	self.nname = blob["name"] 
	self.units = blob["units"] 
	self.runs = blob.setdefault("runs", 1)
	self.dirname = spath[0]

	# create a container for our data
	self.lists = BlobTrain(blob)
	spiked = load_blobfile(spath, self.nname)

	self.lists.add_set(spiked)
	self.lists.proc_set()



#########################################
##
## Storage class for the spike data
##

class SpikeTrain:

    def __init__(self, surface):

	self.ev=[]			    # ordered sequence of spike data
					    # #0 is time in ms, #1 is gid, #2
					    # is time from last spike

	self.tidx={}			    # ms scale index into the other sequence
	self.first = 0
	self.coord={}

	self.pitch = surface["pitch"]
	self.xcenter = surface["xcenter"]
	self.ycenter = surface["ycenter"]
	
	self.cols = surface["cols"]
	self.rows = surface["rows"]
	
	self.xzero = self.xcenter-self.cols*self.pitch/2.0
	self.yzero = self.ycenter-self.rows*self.pitch/2.0


	self.mask = numpy.ones([surface["rows"]+0,
	    surface["cols"]+0], bool)

	coord = surface["coords"]
	# Use nice data structure for fast neighbor search
	if HAVEKDTree:
	    self.coordindex=[]
	    self.anticoord=[]

	    for gid in coord.keys():

		self.coord[int(gid)] = [coord[gid][0], coord[gid][1]]  
		self.mask[coord[gid][1]][coord[gid][0]] = False
		
		# We need this to find nearest neighbors
		self.anticoord.append((coord[gid][0], coord[gid][1]))
		self.coordindex.append(int(gid))

	    self.kdtree = KDTree(self.anticoord)	    # self.kdtree.data contains the
						    # list in the same order
	# Brute force variant
	else:

	    self.anticoord=[]

	    for gid in coord.keys():

		self.coord[int(gid)] = [coord[gid][0], coord[gid][1]]  
		self.mask[coord[gid][1]][coord[gid][0]] = False

		# We need this to find nearest neighbors
		self.anticoord.append(((coord[gid][0], coord[gid][1]),
		    int(gid)))


    # add data to the storage

    def add_set(self, spiked):
	sevapp=self.ev.append
	for neuron, spikes in spiked.iteritems():
	    for s in spikes:
		sevapp([float(s), int(neuron), 0])


    def add_line(self, line):		    # a line of the form "gid    time.dec"

	tmp = line.split()
	self.ev.append([float(tmp[1]), int(tmp[0]), 0]) 
	

    # all data added. Sort and preprocess the data

    def proc_set(self):

	# Empty set - probably no events
	if len(self.ev) == 0:
	    return

	self.ev.sort()			    # automagically sorts in order of
					    # key of sublist, leftmost key
					    # first

	self.first = int(math.floor(self.ev[0][0]))


	gididx = {}
	idx = 0
	i = 0

	for ev in self.ev:    

	    t = ev[0]
	    gid = ev[1]
	   
	    # figure out the ISI for each spike and neuron

	    if gididx.has_key(gid):
		ev[2] = t - gididx[gid]
		gididx[gid] = t
	    else:
		ev[2] = 0
		gididx[gid] = t

	    # create a per-ms index into the spike data so we don't have to
	    # search so much when we plot.
	    #
	    # NOTE: the index points to the next spike not earlier than the
	    # index time. It means that some index times (in the beginning
	    # especially) can point to spikes many milliseconds later than the
	    # index time.

	    ms = int(math.floor(t))

	    while idx<=t:
		self.tidx[idx] = i
		idx = idx+1

	    i = i + 1


    # get spikes covering the range [fromt, tot) 

    def timeslice(self, fromt, tot, debug = False):	
	
	res = []
	fr = int(fromt)
	sev = self.ev
	if (tot>fromt):
	    tt = int(tot)
	else:
	    tt= int(sev[-1][0]+1)

	if not self.tidx.has_key(fr):	    # No such key, so there are no
	    return res			    # spikes later than this

	spike = self.tidx[fr]
	l = len(sev)
	resapp = res.append
	if debug: print "spike, l, sev, sev-1:", spike, l, sev[spike][0], sev[spike-1][0]
	
	while spike<l and sev[spike][0]<tt:
	    resapp(sev[spike])
	    spike = spike + 1
	if debug: print "afterspike, sev, sev-1:", spike, sev[spike][0], sev[spike-1][0]
    
	return res

    # coordinate <-> gid mappings

    def xy_mm(self,x,y):

	return (self.xzero+x*self.pitch, self.yzero+y*self.pitch)

    def mm_xy(self, mmx, mmy):
    
	return ((mmx-self.xzero)/self.pitch, (mmy-self.yzero)/self.pitch)



    # Get closest neuron to coordinates x and y

    def xy_gid(self, x, y):		
	
	if len(self.ev) == 0:
	    return
	if HAVEKDTree:
	    p = self.kdtree.query((x,y))
	    return self.coordindex[p[1]]
	
	else:		    # brute force search. Slower, ugly.
	    tmppt = self.anticoord[0]
	    
	    xd = (tmppt[0][0]-x)
	    yd = (tmppt[0][1]-y)
	   
	    dist = xd*xd + yd*yd
	    
	    gid = tmppt[1]

	    for p in self.anticoord:
		xd = (p[0][0]-x)
		yd = (p[0][1]-y)
		tmpdist = xd*xd + yd*yd
		
		if tmpdist < dist:
		    gid = p[1]
		    dist = tmpdist

	    return gid

    def mm_gid(self, mmx, mmy):
	[x,y] = self.mm_xy(mmx, mmy)
	return self.xy_gid(x, y)


    # Get all neurons within r grid units from (x,y)

    def xyr_gids(self, x, y, r):		
	
	if len(self.ev) == 0:
	    return
	if HAVEKDTree:
	    p = self.kdtree.query_ball_point((x,y), r)
	    res=[]
	    for pp in p:
		res.append(self.coordindex[pp])
	    return res
	
	else:		    # brute force search. Slower, ugly.
	    res = []
	    rs = r*r

	    for p in self.anticoord:
		xd = (p[0][0]-x)
		yd = (p[0][1]-y)
		tmpdist = xd*xd + yd*yd
		
		if tmpdist <= rs:
		    res.append(p[1])

	    return res


    def mmr_gids(self, mmx, mmy, mmr):
	
	[x,y] = self.mm_xy(mmx, mmy)
	r = mmr/self.pitch
	return self.xyr_gids(x, y, r)


    def gid_xy(self, gid):
	return [self.coord[gid][0], self.coord[gid][1]]

    def gid_x(self, gid):
	return self.coord[gid][0]

    def gid_y(self, gid):
	return self.coord[gid][1]


    # Purge the events of any neurons other than gids
    
    def subset(self, gids):
   

	# redo the neighbour lists

	tcoord = {}
	# Use nice data structure for fast neighbor search
	if HAVEKDTree:
	    self.coordindex=[]
	    self.anticoord=[]
	    keylist = [gid for gid in self.coord.keys() if gid in gids]
	    for gid in keylist:
		
		tcoord[int(gid)] = [self.coord[gid][0], self.coord[gid][1]]
		
		# We need this to find nearest neighbors
		self.anticoord.append((self.coord[gid][0], self.coord[gid][1]))
		self.coordindex.append(int(gid))

	    self.kdtree = KDTree(self.anticoord)	    # self.kdtree.data contains the
						    # list in the same order
	    self.coord = tcoord

	# Brute force variant
	else:

	    self.anticoord=[]

	    keylist = [gid for gid in self.coord.keys() if gid in gids]
	    for gid in keylist:
		
		tcoord[int(gid)] = [self.coord[gid][0], self.coord[gid][1]]

		# We need this to find nearest neighbors
		self.anticoord.append(((self.coord[gid][0], self.coord[gid][1]),
		    int(gid)))

	    self.coord = tcoord

	# purge the event list itself and redo the index

	tev = []
	tev = [e for e in self.ev if e[1] in gids]
	self.ev = tev

	self.proc_set()


    # print info (debug)

    def print_set(self):
	i = 0
	for ev in self.ev:    
	    print "[ %8.3f, %5d, %6.1f]" % (ev[0], ev[1], ev[2]),
	    i=i+1
	    if i==4:
		print
		i=0
	print self.tidx

#########################################
##
## Storage class for blob spike data
##

class BlobTrain:

    def __init__(self, blob):

	self.ev=[]			    # ordered sequence of spike data
					    # #0 is time in ms
	self.tidx={}			    # ms scale index into the other sequence
	self.first = 0

    # add data to the storage

    def add_set(self, spiked):

	for spikes in spiked:
	    self.ev.append(float(spikes))

    def add_line(self, line):		    # a line of the form "time.dec"

	self.ev.append(float(line)) 
	

    # all data added. Sort and preprocess the data

    def proc_set(self):

	# Empty set - probably no events
	if len(self.ev) == 0:
	    return

	self.ev.sort()
	
	self.first = int(math.floor(self.ev[0]))

	gididx = {}
	idx = 0
	i = 0

	for ev in self.ev:    

	    t = ev

	    # create a per-ms index into the spike data so we don't have to
	    # search so much when we plot.
	    #
	    # NOTE: the index points to the next spike not earlier than the
	    # index time. It means that some index times (in the beginning
	    # especially) can point to spikes many milliseconds later than the
	    # index time.

	    ms = int(math.floor(t))

	    while idx<t:
		self.tidx[idx] = i
		idx = idx+1

	    i = i + 1


    # get spikes covering the range [fromt, tot) 

    def timeslice(self, fromt, tot):	
	
	res = []
	fr = int(fromt)
	
	if (tot>fromt):
	    tt = int(tot)
	else:
	    tt= int(self.ev[-1]+1)

	if not self.tidx.has_key(fr):	    # No such key, so there are no
	    return res			    # spikes later than this

	spike = self.tidx[fr]
	l = len(self.ev)
	resapp=res.append
	while spike<l and self.ev[spike]<tt:
	    resapp(self.ev[spike])    
	    spike = spike + 1

	return res

    # print info (debug)

    def print_set(self):
	i = 0
	for ev in self.ev:    
	    print "[ %8.3f]" % (ev),
	    i=i+1
	    if i==4:
		print
		i=0
	print self.tidx



