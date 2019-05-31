#### SGI model class library ####
#
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

import nest
from nest import topology
from math import sqrt, sin, cos, exp, pi, log, atan2
import time
import sys

##############################
##
## Class for a neuron layer 
##

class Surface:
    
    def __init__(self, row, col, pitch, A=0.0, bx=0.0, by=0.0):
	self.rows   = row
	self.cols   = col
	self.pitch  = pitch
	
	self.xdim	= (col)*pitch
	self.ydim	= (row)*pitch
	self.xcenter	= self.xdim/2.0
	self.ycenter	= 0.0
	self.simname	= ""
	self.nname	= ""

	if A>0:
	    self.trans = True
	    self.A = A
	    self.bx = bx
	    self.by = by
	else:
	    self.trans = False


    def create_layer(self, nname, inparams={}):

	self.nname = nname
	self.net = topology.CreateLayer({
	    "rows": self.rows,
            "columns": self.cols,
            "extent": [self.xdim, self.ydim],
	    "center": [self.xcenter, self.ycenter],
            "elements": self.nname,
	    "edge_wrap": False})


#	lpstart=time.time()
	self.nmdaw = {}
	if inparams!={}:
	    try:
		params=inparams[nname]
	    except:
		print "%s is not defined in parameter file" % (nname)
		exit (-1)
	    
	    self.xy = params["XY"]	    # grid position
	    self.sxy = params["SXY"]	    # surface position
	    self.freezed = params["freeze"]
	    self.frozen= params["frozen"]

	else:
	
	    # node lookup tables	
	    
	    self.xy = {}	    # grid position
	    self.sxy = {}	    # surface position
	    self.freezed = []
	    self.frozen=0
	    isnode = 0


	    for x in range(self.cols):
		for y in range(self.rows):
		    node_id = topology.GetElement(self.net,[x,y])[0]
	
		    islocal = nest.GetStatus([node_id], 'local')
		    if (islocal[0]):
			[node_p] = topology.GetPosition([node_id])
			(retx, rety) = self.to_ret_xy(node_p)
			if retx<=0:
			    inside = False
			else:
			    inside = True

			if ((not self.trans) or (inside)):
			    
			    self.xy[node_id] = (x, y)
			    self.sxy[node_id] = node_p
			if (not inside):
			    self.freezed.append(node_id)
			    self.frozen+=1


#	lpend=time.time()
#	print "Rank: %d: map time: %.1f s frozen: %d" %(nest.Rank(),
#		lpend-lpstart, self.frozen)

#	lpstart=time.time()

	# Add spike detector
	self.sd = nest.Create("spike_detector")
	nest.SetStatus(self.sd, {"withgid": True })

	nest.ConvergentConnect(nest.GetLeaves(self.net)[0], self.sd)
    
#	lpend=time.time()
#	print "Rank: %d: spike detector time: %.1f s" %(nest.Rank(), lpend-lpstart)

    # Create a non-grid layer with a list of specific node positions
    # return a list of gids, guaranteed in the same order as the input list of
    # node positions
    
    def create_free_layer(self, nname, nodelist):

	self.nname = nname

	self.net = topology.CreateLayer({
            "extent": [self.xdim, self.ydim],
	    "center": [self.xcenter, self.ycenter],
            "elements": self.nname,
	    "edge_wrap": False,
	    "positions": nodelist})

#	lpstart=time.time()


	# node lookup tables	

	self.xy = {}	    # grid position - none in this case
	self.sxy = {}	    # surface position
	self.freezed = []   # not used
	self.frozen=0

	for node_id in nest.GetLeaves(self.net)[0]:
	    islocal = nest.GetStatus([node_id], 'local')
	    if (islocal):
		[node_p] = topology.GetPosition([node_id])
		self.sxy[node_id] = node_p
	
	outlist = []
	for nd in nodelist:

	    md=1000000000
	    nodegid = -1
	    for gid, pt in self.sxy.iteritems(): 
		d = ((pt[0]-nd[0])**2+(pt[1]-nd[1])**2)
		if d<md:
		    nodegid = gid
		    md = d

	    outlist.append(nodegid)

	# Add spike detector
	self.sd = nest.Create("spike_detector")
	nest.SetStatus(self.sd, {"withgid": True })

	nest.ConvergentConnect(nest.GetLeaves(self.net)[0], self.sd)
	
	return outlist

    # Optionally cache the configuration values. Speeds up create_layer() by
    # 2/3, from 30 seconds or so, down to 10. For NEST 2.2 and onwards, it is 
    # necessary to save configuration values, then use the saved values for 
    # multi-process simulations.

    def save_params(self):
	
	params = {}
	params["name"] = self.nname
	params["XY"] = self.xy
	params["SXY"] = self.sxy
	params["freeze"] = self.freezed
	params["frozen"] = self.frozen

	return params

    # we prepare data to save on disk

    def save_data(self, savename):

	
	sname = savename+"_"+self.nname
	nest.SetStatus(self.sd, {"withgid": True, "to_file": True, 
	    "label": sname, "file_extension": "spikes",  "fbuffer_size": 8192})

	surface = {}

	surface["rows"] = self.rows
	surface["cols"] = self.cols
	surface["pitch"] = self.pitch
	surface["name"] = self.nname
	surface["xcenter"] = self.xcenter
	surface["ycenter"] = self.ycenter
	surface["filebasename"] = self.nname
	surface["GID"] = str(nest.GetStatus(self.sd)[0]['global_id'])
	surface["coords"] = self.xy

	return surface


    # Get the nodes in the list that are stored locally on this specific
    # process. Some functions in topology especially can only operate on local
    # data.
    def get_local(self, nodes):

	islocal = nest.GetStatus(nodes, 'local')
        l_nodes = [nodes[i] for i in xrange (len (nodes)) if islocal[i]]
	return l_nodes


    # For real-world simulations we need two colliculi, and a way to connect
    # them. This implements our optional interconection.

    # Create a connection between self and the target layer, where self and
    # target are assumed to belong to opposite SCs. Creates a connection
    # according to conndict, but shifted so the two SC seem to make up a
    # single continuous surface.
    #
    # See [Tabareaux et. al.] for why this simplification is actually quite
    # wrong. But the results are close enough to be useful for some types of
    # simulation.

    def intraconnect(self, sources, pyrngs, conndict):

	# A general search for anchor terms to change, or adding a 
	# new anchor term if there should be one.
	# This is not well tested; some legal connection dictionaries may
	# break.

	# For sanity, we check that mask is always circular, kernel and delay
	# is always constant and weights is always gaussian or constant.

	weight = {}
	mask = {}
	kernel = 0.0
	delay = 0.0
	synapse = "static_synapse"
	do_nmda = False
	
	if (conndict.has_key("synapse_model")):
	    synapse = conndict['synapse_model']

	if (synapse == "nmda_synapse"):
	    do_nmda = True

	if (conndict.has_key("mask") and 
		conndict["mask"].has_key("circular") and 
		conndict["mask"]["circular"].has_key("radius")) :	
	
	    mask["anchor"] = conndict["mask"]["circular"].get("anchor", [0.0,0.0])
	    mask["radius"] = conndict["mask"]["circular"]["radius"]
	
	else:
	    print "intraconnect: no valid mask parameter: ", conndict
	    exit (-1)

	if (conndict.has_key("kernel") and 
		type(conndict["kernel"]) == type (1.0)):
	    
	    kernel = conndict["kernel"]
	
	else: 
	    print "intraconnect: No valid kernel parameter: ", conndict
	    exit (-1)
	
	if (conndict.has_key("delays") and
		type(conndict["delays"]) == type(1.0)):

	    delay = conndict["delays"]
	
	else:
	    print "intraconnect: no valid delay parameter: ", conndict
	    exit (-1)

	if (conndict.has_key("weights") and 
		type(conndict["weights"]) ==type (1.0)):
	    
	    weight["wt"] = conndict["weights"]
	    weight["anchor"] = [0.0, 0.0]
	    weight["constant"] = True

	elif (conndict.has_key("weights") and 
		conndict["weights"].has_key("gaussian") and 
		conndict["weights"]["gaussian"].has_key("p_center") and
		conndict["weights"]["gaussian"].has_key("sigma")):

	    weight["wt"] = conndict["weights"]["gaussian"]["p_center"]
	    weight["sigma"] = conndict["weights"]["gaussian"]["sigma"]
	    weight["anchor"] = conndict["weights"]["gaussian"].get("anchor", [0.0,0.0])
	    weight["constant"] = False

	else:
	    print "intraconnect: no valid weight parameter: ", conndict
	    exit(-1)

	
	# For each target neuron:
	#   pick possible source neurons according to mask,
	#   if weight is constant:
	#	connect sources to target with kernel, delay and weight
	#   if gaussian:
	#	figure out weight for each source:
	#	    connect each source to target with params as above


	# Select only local target nodes

        nodes = self.sxy.keys()
	l_nodes = self.get_local(nodes)
        #islocal = nest.GetStatus(nodes, 'local')
        #l_nodes = [nodes[i] for i in xrange (len (nodes)) if islocal[i]]
        l_xy = nest.GetStatus(l_nodes, ['global_id', 'vp'])
	#node_info = nest.GetStatus(self.sxy.keys(), ['global_id', 'vp', 'local'])
	#l_xy = [(gid, virtp) for gid, virtp, islocal in node_info if islocal]
	
	for tgid, virtp in l_xy: 

	    pt = self.sxy[tgid]

	    # Projection mask and weight
	    ptmask = ((pt[0]-mask["anchor"][0]), (pt[1]-mask["anchor"][1]))
	    ptweight = ((pt[0]-weight["anchor"][0]), (pt[1]-weight["anchor"][1]))


	    # SC -> retinal angular, cross Y=0 -> SC
	    (rx, ry) = self.to_ret_xy(pt)	    
	    pt_t = self.retxy_to_xy((-rx,ry))	    # projected neuron center
	    (rx, ry) = self.to_ret_xy(ptmask)
	    ptmask_t = self.retxy_to_xy((-rx,ry))   # projected mask center
	    (rx, ry) = self.to_ret_xy(ptweight)
	    ptwt_t = self.retxy_to_xy((-rx,ry))	    # projected weight center

	    # pick source neurons, filter those that are not on the local node
	    
	    sneurons=sources.xyr_gids(ptmask_t[0], ptmask_t[1], mask['radius'])

	    sxy=[]
	    
	    if kernel==1.0:
		sxy = sneurons
	    else:
		sxy = [gid for gid in sneurons if
			(pyrngs[virtp].uniform()<kernel)]

	    # constant weight
	    if len(sxy)>0:
		if weight["constant"]:
		    nest.ConvergentConnect(sxy, [tgid], weight["wt"], delay)

		# gaussian weight 
		else:
		    weights = []
		    delays = [delay]*len(sxy)
		    if do_nmda == True:
			w = self.nmdaw[tgid]*3.0
		    else:
			w = weight["wt"]*1.0
		    for gid in sxy:

			(x,y) = sources.sxy[gid]
			d2 = ((ptwt_t[0]-x)**2 + (ptwt_t[1]-y)**2)
			gw = w*exp(-(d2)/(2*weight["sigma"]**2))
			weights.append(float(gw))

		    nest.ConvergentConnect(sxy, [tgid], weight=weights,
			    delay=delays, model=synapse)


    # create a distributed NMDA connection from QV to burst neurons (this
    # layer)

    def nmdacon(self, sources, pyrngs, indict):

	self.inhdict = {
	    "point_a": 2.0,
	    "val_a": 1.0,
	    "point_b": 35.0,
	    "val_b": 0.8,
	    "radius": 1.5,
	    "sigma": 0.5,
	    "kernel": 0.25,
	    "delay": 1.0,
	    "synapse": "nmda_synapse"}	

	for x in indict:
	    if not x in self.inhdict:
		print "Error: non-existent key '%s' input to NMDA mapping.\n" % (x)
		sys.exit(-1)
	    self.inhdict[x] = indict[x]

	pa = float(self.inhdict['point_a'])
	va = float(self.inhdict['val_a'])
	pb = float(self.inhdict['point_b'])
	vb = float(self.inhdict['val_b'])

	delay = float(self.inhdict['delay'])
	synapse = self.inhdict['synapse']

	radius = float(self.inhdict['radius'])
	sigma = float(self.inhdict['sigma'])
	kernel = float(self.inhdict['kernel'])

	(dax, day) = self.retxy_to_xy((pa, 0.0))
	(dbx, dby) = self.retxy_to_xy((pb, 0.0))
	 
	# slope and intercept for the linear distribution of weights between
	# angles a and b.
 	m = (va-vb)/(dax-dbx)	
	b = va-m*dax

        # Find only neurons that are stored locally, so we don't connect all 
        # neurons in parallel across all computing nodes 
	tmpnodes = self.sxy.keys()
        l_nodes = self.get_local(tmpnodes)
        l_xy = nest.GetStatus(l_nodes, ['global_id', 'vp'])
	for gid, virtp in l_xy:

	    xy = self.sxy[gid]
    
	    # sc to retinal angular, to polar, to sc distance coords
	    ret = self.to_ret_xy(xy)
	    (r, theta) = self.retxy_to_rth(ret)
	    (px, py) =self.retxy_to_xy((r,0.0))

	    w = float(m*px+b)
	    self.nmdaw[gid] = w

	    # Pick the set of source nodes in sources within the radius,
	    # get the distance and scale according to gaussian with sigma,
	    # then set connections with kernel probability.
	    weights = []
	    delays = []
	    s_list = sources.xyr_gids(xy[0],xy[1], radius)

	    # pick a random subset
	    elems = int(kernel*len(s_list))
	    t_list = []
	    for i in range(elems):

		t_list.append(self.lselect(s_list, pyrngs[virtp]))
	    
	    # Connect the selected nodes 
	    for t in t_list:
	    
		src = sources.sxy[t]
		d = sqrt((xy[0]-src[0])**2 + (xy[1]-src[1])**2)
	    
	    
		gw = w*exp(-(d*d)/(2*sigma*sigma))
		weights.append(float(gw))

	
	    delays=[delay]*len(t_list)
	    nest.ConvergentConnect(t_list, [gid], weight=weights,
		    delay=delays, model=synapse)
	
    
    # Select a random element in the input list, removing the element

    def lselect(self, data, pyrng):
	if data != []:
	    index = int(pyrng.uniform(0, len(data)))
	    elem = data[index]
	    data[index] = data[-1]
	    del data[-1]
	    return elem
	else:
	    return data

    # create a distributed inhibitory connection to burst neurons (this layer). 
    def inhibcon(self, sources, indict):

	self.inhdict = {
	    "point_a": 2.0,
	    "val_a": 1.0,
	    "point_b": 35.0,
	    "val_b": 0.8,
	    "delay": 1.0,
	    "synapse": "static_synapse"}	

	for x in indict:

	    if not x in self.inhdict:

		print "Error: non-existent key '%s' input to inhibitory mapping.\n" % (x)
		sys.exit(-1)

	    self.inhdict[x] = indict[x]

	pa = float(self.inhdict['point_a'])
	va = float(self.inhdict['val_a'])
	pb = float(self.inhdict['point_b'])
	vb = float(self.inhdict['val_b'])

	delay = float(self.inhdict['delay'])
	synapse = self.inhdict['synapse']

	(dax, day) = self.retxy_to_xy((pa, 0.0))
	(dbx, dby) = self.retxy_to_xy((pb, 0.0))
	
 	m = (va-vb)/(dax-dbx)	
	b = va-m*dax
    

	for gid, xy in self.sxy.iteritems():
	    islocal = nest.GetStatus([gid], 'local')
	    if (islocal[0]):
	    
		# sc to retinal angular, to polar, to sc distance coords
		ret = self.to_ret_xy(xy)
		(r, theta) = self.retxy_to_rth(ret)
		(px, py) =self.retxy_to_xy((r,0.0))

		w = float(m*px+b)

		nest.ConvergentConnect(sources.blob, [gid], w, delay, synapse)


    
    # Create a log-polar connection. It's easier (if slower) to do this at
    # this level rather than to patch the topology subsystem.

    def logpolar(self, target, indict):
	

	nest.Rank() == 0
	self.pdict = {
		'delay': 1.0,
		'radius': 1.0,
		'variance': 1.0,
		'weight': 1.0,
		'shift': 0.0}	    # leftward shift of radius for finding
				    # affected nodes
		
	for x in indict:

	    if not x in self.pdict:

		print "Error: non-existent key '%s' input to log-polar mapping.\n" % (x)
		sys.exit(-1)

	    self.pdict[x] = indict[x]

	d = self.pdict['delay']*1.0
	r = self.pdict['radius']*1.0
	v = self.pdict['variance']*1.0
	wt = self.pdict['weight']*1.0
	shift= self.pdict['shift']*1.0
	
	for xy in self.sxy:

	    [x,y] =self.sxy[xy]
	    t_list = target.xyr_gids(x-shift,y, r)
	    t_filtlist = self.get_local(t_list)
	    self.lp_unit(target, xy, t_filtlist, r, v, wt,d)


    # Connect origin with targets according to log-polar mapping, modified
    # by gaussian

    def lp_unit(self, target, origin, targets, radius, variance, weight, delay):

	# all need to be vectors 
	weights = []
	delays = []
	maxd= 0.0	
	o = self.sxy[origin]
	orig = self.to_ret_xy(o)
	
	if (o[1]<0):
	    radius = -radius

	orig2 = self.to_ret_xy((o[0], o[1] - radius))
	retrad = sqrt((orig[0]-orig2[0])**2 + (orig[1]-orig2[1])**2)
	var = variance*(retrad/radius)
	
	for t in targets:
	    
	    trg = self.to_ret_xy(target.sxy[t])
	    d = sqrt((orig[0]-trg[0])**2 + (orig[1]-trg[1])**2)
	    
	    if (d>maxd):
		maxd = d
	    
	    gw = weight*exp(-(d*d)/(2*var*var))
	    weights.append(gw)
	
	delays=[delay]*len(targets)
	nest.DivergentConnect([origin], targets, weights, delays)


    # list all nodes within a radius r around (x,y)

    def xyr_gids(self, x, y, r):		
	
	res = []
	rs = r*r

	for xyp in self.sxy:
	    [p0, p1] = self.sxy[xyp]
	    xd = (p0-x)
	    yd = (p1-y)
	    tmpdist = xd*xd + yd*yd
	    
	    if tmpdist <= rs:
		res.append(xyp)

	return res


    # finds all nodes in layer "source" that connect to node "unit" in this layer
    # used only for debuggging
    def connection_from(self, unit, source):
	
	sgids = nest.GetLeaves(source.net)[0]
	slist = []

	for sg in sgids:

	    t=topology.GetTargetNodes([sg], self.net)[0]
	    if unit in t:

		slist.append(sg)

	return slist

    # Freeze neurons that are outside the actual (non-rectangular) layer. Once
    # frozen they will not influence calculations and not take take any
    # computation time. Note that they will count toward the total connection
    # limit, if any, for randomized connections, retaining the expected
    # connection density near edges.
    
    def freeze(self):
	nest.SetStatus(self.freezed, "frozen", True)
	

    # SC to retinal angular coordinates
    def to_ret_xy(self, (x, y)):
	
	if self.trans:

	    ex = self.A*(exp(x/self.bx)*cos(y/self.by) - 1)
	    ey = self.A*(exp(x/self.bx)*sin(y/self.by))
	    
	    return (ex, ey)

	else:
	    return (x, y)

    # Retinal angular coordinates to SC surface
    def retxy_to_xy(self, (ex, ey)):

	z = complex (ex, ey)
	t1=(z+self.A)/self.A
	rt=t1.real
	ri=t1.imag
	p = log(sqrt(rt*rt+ri*ri))+complex(0, atan2(ri,rt))

	px = p.real*self.bx
	py = p.imag*self.by

	return (px, py)
    
    # retinal angular to retinal polar coordinates

    def retxy_to_rth(self, (xin, yin)):
	x = xin*pi/180.0
	y = yin*pi/180.0

	ra = sqrt((sin(x)*cos(y))**2 + sin(y)**2)
	rb = cos(x)*cos(y)
	
	r = atan2(ra, rb)*180/pi
	theta = atan2(y,x)*180/pi

	return(r,theta)


# Class to keep track of and save non-spatial components in the model

class Blob:

    def __init__(self):
	self.simname	= ""
	self.nname	= ""
	self.nrunits    = 1

    def create_blob(self, nname, nr):

	self.nname = nname
	self.nrunits = nr
	self.blob = nest.Create(self.nname, self.nrunits)

	self.sd = nest.Create("spike_detector")
	nest.SetStatus(self.sd, {"withgid": False })

	nest.ConvergentConnect(self.blob, self.sd)
    
    # we save the data on disk

    def save_data(self, savename):

	sname = savename+"_"+self.nname
	nest.SetStatus(self.sd, {"withgid": False, "to_file": True, 
	    "label": sname, "file_extension": "spikes",  "fbuffer_size": 8192})

	blob = {}

	blob["units"] = self.nrunits
	blob["name"] = self.nname
	blob["filebasename"] = self.nname

	return blob


