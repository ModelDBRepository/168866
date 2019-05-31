#!/usr/bin/python
#### SGI model main file ####
#
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



# Only use this when actually running batch code. Don't want to miss
# important stuff
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


import nest
from nest import topology
import math
import os, os.path
import re
import errno
import sys
import getopt
import time
import pickle
import numpy
import numpy.random
import json
from copy import deepcopy

from SGI2class import *



nest.sli_run("M_ERROR setverbosity")
nest.SetKernelStatus({"resolution": 0.1})
#nest.SetKernelStatus({"local_num_threads": 8})
#nest.pynestkernel.logstdout("/dev/null")

debug = False
#debug = True
#print "Rank: %d" %(nest.Rank())
###############################################
#
# Parameters
#

## Cols = X = horizontal extent, from 0.0 to xdim
## Rows = Y = vertical extent, from ydim/2 to -ydim/2


paramsets = {}

paramset="original"


# These parameter settings reproduce the model as used in [Moren2013]

paramsets['original'] = {
	'xsize':    4.81,	# surface size, mm
	'ysize':    2.77*2,
	'spar':	    1.0,	# density adjustment parameter
	'in_pitch': 0.0375,	# interneuron distance, mm
	'd_red':    2.0,	# density reduction for sparse layers
	'kernprop': 0.25,	# kernel connection density

	'build':    0.1,	# build-burst strength
	'fback':    0.12,	# build-inhibitory strength

	'a_angle':  2.0,	# Angles for the below parameters
	'b_angle':  35.0,
	'a_nmda':   0.203,	# qv-burst nmda strength by angle
	'b_nmda':   0.1144,
	'a_inhib':  0.374,	# integrator-burst inhibitory strengh
	'b_inhib':  0.2359}


# Small-scale model, used for interfacing with a robot. Not adapted for
# current NEST versions; unlikely to work properly atm.

paramsets['smallscale'] = {
	'xsize':    4.81,	# surface size, mm
	'ysize':    2.77*2,
	'spar':	    1.78,	# density adjustment parameter
	'in_pitch': 0.1,	# interneuron distance, mm
	'd_red':    1.5,	# density reduction for sparse layers
	'kernprop': 1.0,	# kernel connection density

	'build':    0.12,	# build-burst strength
	'fback':    0.145,	# build-inhibitory strength

	'a_angle':  2.0,	# Angles for the below parameters
	'b_angle':  35.0,
	'a_nmda':   0.24,	# qv-burst nmda strength by angle
	'b_nmda':   0.21,
	'a_inhib':  0.29,	# integrator-burst inhibitory strengh
	'b_inhib':  0.22}


xsize	    = paramsets[paramset]['xsize']
ysize	    = paramsets[paramset]['ysize']
spar	    = paramsets[paramset]['spar']
in_pitch    = paramsets[paramset]['in_pitch']
outf	    = paramsets[paramset]['d_red']
kernprop    = paramsets[paramset]['kernprop']

buildpar    = paramsets[paramset]['build']
fback	    = paramsets[paramset]['fback']

a_angle	    = paramsets[paramset]['a_angle']
b_angle	    = paramsets[paramset]['b_angle']
a_nmda	    = paramsets[paramset]['a_nmda']
b_nmda	    = paramsets[paramset]['b_nmda']
a_inhib	    = paramsets[paramset]['a_inhib']
b_inhib	    = paramsets[paramset]['b_inhib']


#
# Surface parameters
#


in_cols = int(xsize/in_pitch) 
in_rows = int(ysize/in_pitch)

out_cols = int(in_cols/outf)	    # number of rows and columns for output layers
out_rows = int(in_rows/outf)    
out_pitch = xsize*1.0/out_cols

if nest.Rank() == 0 and debug:
    print "full: (%d, %d) = %d \t Reduced: (%d, %d) = %d \t pitch %f - %f \n" % (in_cols,
    	in_rows, in_cols*in_rows, out_cols, out_rows, out_cols*out_rows,
	in_pitch, out_pitch)


# Parameters for the synthetic inhibitory input
# Really obsolete, as the input is not actually spatial at all. 

syn_cols = 7		    
syn_rows = 7
syn_pitch = xsize*1.0/syn_cols
xin = 4
yin = 3


# these parameters are used only without angle-dependent parameters when
# tuning the model, and then we always set these through the command line.
inhibpar = 0.4 
nmdapar  = 0.21


# Monkey retina-collicular mapping parameters

A = 3.0		   
bx = 1.4
by = 1.8


# basic synaptic conductance values

ampaw =0.72
nmdaw = 1.2
gabaw = -0.04

# For our simulations - fewer connections mean compensating for it

amw = ampaw*2.0
nmw = nmdaw*2.0
gbw = gabaw *2.0
inw = ampaw
syninw = amw

# Default random seed
sseed=1
# Number of OpenMP threads to use
ompnr=1

synbuildinhw = gabaw*.095 # integrator inhibition of inhibitory neurons (0.095)
syninhinhw = gabaw*3.0	# inhibitory signal to the integrator from SNpr and
			# deep inhibitory neuron layer

synconspeed = 5.0	
synconinspeed = 5.0
synconinhibspeed = 10.0

# initial input
spikerate = 200
spikepop = 500
inhspikerate = 100*120

# default angle for input stimulus. Normally set through the command line
retx = 9.0
rety = 0.0

simtime = 400			# milliseconds


###############################################
#
# General initialization 
#


def show_help():

    print """
Usage: %s [-nibs <val>] [-l] [--savepar] [-h] output

  -n, --nmda=VAL	set QV->burst nmda synapse strength factor
  -i, --inhib=VAL	set cMRF->burst inhibitory strength factor

  NOTE: if either of nmda or inhib is not set, the strengths of the parameter
  not set will be distributed according to the preset parameters in the
  beginning of this file. You need to set both to keep both fixed.

  -b, --build=VAL	set build->burst strength factor
  --fback=VAL		set build ->deep inhibitory feedback strength

  --retx=VAL		Retinal angular coordinates (9, 0)			 
  --rety=VAL
  --rseed=VAL		set initial random seed (integer)
  --spikerate=VAL	set input rate (x500 neurons) 

  --stime=VAL		simulation time (milliseconds)

  --music		Use MUSIC integration
  --omp=VAL		How many OpenMP threads to use (1)

  --nosave		Do not save any simulation results on disk
  --noinhib		Do not activate deep layer inhibitory feedback

  -l, --loadpar		load surface parameters from a precalculated 
			parameter file (%s.par)
  --savepar		calculate and save a surface parameter file
			IMPORTANT: do not run on multiple CPUs

  -h, --help		show a brief help text
    """ % (simname, simname)
simname = os.path.basename(os.path.splitext(sys.argv[0])[0])

argv = sys.argv[1:]

savedata = True	    # Whether to save data files for later analysis
loadparams = False
saveparams = False
doinhibition = True
domusic = False
try:
    opts, args = getopt.getopt(argv, "n:i:b:l", ["nmda=", "inhib=",
	"build=","rseed=", "fback=", "retx=", "rety=", "spikerate=",
	"stime=", "omp=", "noinhib", "loadpar", "savepar", "music", "nosave"])
except:
    show_help()
    sys.exit(-1)

for opt, arg in opts:
    if opt == "--retx":
	retx = float(arg)
    if opt == "--rety":
	rety = float(arg)
    if opt == "--rseed":
	sseed = int(arg)
    if opt == "--spikerate":
	spikerate = int(arg)
    if opt == "--stime":
	simtime = int(arg)
    if opt == "--music":
	domusic=True
    if opt == "--omp":
	ompnr = int(arg)
    if opt == "--nosave":
	savedata=False

    if opt == "--noinhib":
	doinhibition=False

    if opt == "--fback":
	fback = float(arg)

    if opt in ("-n", "--nmda"):
	nmdapar = float(arg)
	a_nmda = nmdapar
	b_nmda = nmdapar
    if opt in ("-i", "--inhib"):
	inhibpar = float(arg)
	a_inhib = inhibpar
	b_inhib = inhibpar
    if opt in ("-b", "--build"):
	buildpar = float(arg)

    if opt in ("-l", "--loadpar"):
	loadparams = True


    if opt == "--savepar":
	saveparams = True
	loadparams = False

if len(args)<1 and not (saveparams or (not savedata)):
    show_help()
    print "ERROR: Specify one output name.\n"
    sys.exit(-1)

if len(args) == 1:
    savename = args[0]
elif len(args) > 1:
    show_help()
    print "ERROR: The script takes only one file argument\n"
    sys.exit(-1)
else:
    savename = "save"

nest.SetKernelStatus({"local_num_threads": ompnr})

nproc = nest.GetKernelStatus('total_num_virtual_procs')
seeds = (sseed-1)*(nproc*2+1)+1
srange = range(seeds,seeds+nproc)

nest.SetKernelStatus({'rng_seeds': srange})
nest.SetKernelStatus({'grng_seed' : (seeds+nproc)})
pyrngs = [numpy.random.RandomState(s) for s in range(seeds+nproc+1,seeds+2*nproc+1)]

## Sanity check

if nproc>1:
    if loadparams == False and saveparams == False:
	print """
    ERROR: you must load parameters from a parameter file when using more 
    than one process for the simulation.

    First run once using _one_ process only to generate the parameter file:

	python [yourmodel.py] --savepar
    """
	sys.exit(-1)

    if saveparams == True:
	print """
    ERROR: You must run "--saveparams" using only one process and no threads:

	python [yourmodel.py] --savepar
    """
	sys.exit(-1)



#### Parameters for the neuron models ###

# General dictionary

gen_dict = {
    'V_peak'    :  0.0,
    'V_reset'   :-65.0,
    't_ref'     :  0.0,
    'g_L'       :  4.0,	# 2.0
    'C_m'       : 62.0,
    'E_ex'      :  0.0,
    'E_in'      :-75.0,
    'E_L'       :-65.0,
    'V_m'	:-65.0,
    'Delta_T'   :  2.0,
    'tau_w'     : 20.0, # 20.0
    'a'         :  0.2,
    'b'         :  30.0,
    'V_th'      :-47.0,
    'tau_syn_ex':  0.2,
    'tau_syn_in':  3.0,
    'I_e'       :  0.0,
    'NMDA_V_max':-43.6,
    'NMDA_V_min':-60.0,
    'NMDA_gain'	:  3.0,
    'E_n'       :  0.0,
    'tau_syn_n' :  3.0
    }

# The synthetic integrator pool parameters
intdict = {
    'Smax': 100.0, 
    'Gamma': 1.,
    'Reset': 0.5
    }
# number of integrators in the pool
synnr = 100

##################
#
# SGS and SO
# 

lpstart=time.time()
paramdata={}
fname = simname+".par"

if loadparams:
    try:
	if debug: 
	    print "Loading parameter file"
	f = open(fname, "rb")
	try:
	   #paramdata = json.load(f)
	   paramdata = pickle.load(f)
	finally:
	    f.close()
    except IOError:
	
	print "Failed to load parameter data from", fname,"\n"
	sys.exit(-1)

## wide field neuron

wifi = Surface(in_rows, in_cols, in_pitch, A, bx, by)

wifi.dict = gen_dict
    	    
nest.CopyModel("aeif_cond_nmda_alpha", "wifi_neuron", wifi.dict)
wifi.create_layer("wifi_neuron", paramdata)


## narrow field neuron

nafi = Surface(in_rows, in_cols, in_pitch, A, bx, by)

nafi.dict = gen_dict
    	    
nest.CopyModel("aeif_cond_nmda_alpha", "nafi_neuron", nafi.dict)
nafi.create_layer("nafi_neuron", paramdata)



#################
#
# SGI layer neurons
#

#
## Quasivisual neuron
#
qv = Surface(in_rows, in_cols, in_pitch, A, bx, by)

qv.dict = deepcopy(gen_dict)
    	   # Fast GABAa
qv.dict['tau_syn_in'] = 1.5
 
nest.CopyModel("aeif_cond_nmda_alpha", "qv_neuron", qv.dict)
qv.create_layer("qv_neuron", paramdata)

#
## qv inhibitory neuron
#
qvinh = Surface(in_rows, in_cols, in_pitch, A, bx, by)

qvinh.dict = deepcopy(gen_dict)
# Fast GABAa
qvinh.dict['tau_syn_in'] = 1.5

nest.CopyModel("aeif_cond_nmda_alpha", "qvinhup_neuron", qvinh.dict)
qvinh.create_layer("qvinhup_neuron", paramdata)


#
## buildup neuron
#

build = Surface(in_rows, in_cols, in_pitch, A, bx, by)

build.dict = deepcopy(gen_dict)

# Fast GABAa
build.dict['tau_syn_in'] = 1.5

nest.CopyModel("aeif_cond_nmda_alpha", "buildup_neuron", build.dict)
build.create_layer("buildup_neuron", paramdata)

#
## buildup inhibitory neuron
#
buildinh = Surface(in_rows, in_cols, in_pitch, A, bx, by)

buildinh.dict = deepcopy(gen_dict)
# Fast GABAa
buildinh.dict['tau_syn_in'] = 1.5

nest.CopyModel("aeif_cond_nmda_alpha", "buildinhup_neuron", buildinh.dict)
buildinh.create_layer("buildinhup_neuron", paramdata)


#
## burst neuron
#

burst = Surface(out_rows, out_cols, out_pitch, A, bx, by)

burst.dict = deepcopy(gen_dict)
#burst.dict['I_e'] = 50.0
burst.dict['C_m'] = 40.0
burst.dict['tau_syn_in'] = 1.5

nest.CopyModel("aeif_cond_nmda_alpha", "burst_neuron", burst.dict)
burst.create_layer("burst_neuron", paramdata)


# Input port for the synapses need to be set at this level, rather than in the
# topological connections

nmdasynapsedict = {
	'receptor_type':1
	}
nest.CopyModel("static_synapse", "nmda_synapse", nmdasynapsedict)

#
## Inhibitory neuron
#
inhib = Surface(out_rows, out_cols, out_pitch, A, bx, by)

inhib.dict = deepcopy(gen_dict)

inhib.dict['tau_syn_in'] = 1.5
inhib.dict['C_m'] = 62.0

nest.CopyModel("aeif_cond_nmda_alpha", "inhib_neuron", inhib.dict)
inhib.create_layer("inhib_neuron", paramdata)

#
## Synthetic integrator neuron pool
#
integ = Blob()

integ.dict=deepcopy(intdict)
nest.CopyModel("synth_integrator", "integrate", integ.dict)

integ.create_blob("integrate", synnr)


# 
# Cache parameters for speedier startup
#


if saveparams:
    paramdata[wifi.nname] = wifi.save_params()
    paramdata[nafi.nname] = nafi.save_params()
    paramdata[qv.nname] = qv.save_params()
    paramdata[qvinh.nname] = qvinh.save_params()
    paramdata[build.nname] = build.save_params()
    paramdata[buildinh.nname] = buildinh.save_params()
    paramdata[burst.nname] = burst.save_params()
    paramdata[inhib.nname] = inhib.save_params()


    try:
	f = open(fname, "wb")
	try:
	    pickle.dump(paramdata, f, -1)
	finally:
	    f.close()
	    if not domusic:
		print "Parameters saved. Exiting."
		sys.exit(0)
    except IOError:
	
	print "Failed to save parameter file."
	sys.exit(-1)



lpend=time.time()
if debug:
    print "Rank: %d: Model creation time: %.1f s" %(nest.Rank(), lpend-lpstart)


lpstart=time.time()
#################################################
#
# Connections
#


# Wide-field projection to QV - may be weak.

wifi_qv_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
		"circular": {"radius": 0.20}},
	"delays": 1.0,
	"weights": amw*spar, # should be on
	"kernel": kernprop,
	"allow_autapses": False,
	"allow_multapses": False
	}

# Wide-field projection to build neurons; will be asymmetrical

wifi_build_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
		"circular": {"radius": 0.20}},
	"delays": 1.0,
	"weights": amw*spar,
	"kernel": kernprop,
	"allow_autapses": False,
	"allow_multapses": False
	}

# Narrow-field projection to QV

nafi_qv_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
		"circular": {"radius": 0.20}},
	"delays": 1.0,
	"weights": amw*spar, #should be 2.0
	"kernel": kernprop,
	"allow_autapses": False,
	"allow_multapses": False
	}

# Narrow-field projection to build neurons


nafi_build_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
		"circular": {"radius": 0.20}},
	"delays": 1.0,
	"weights": amw*spar,
	"kernel": kernprop,
	"allow_autapses": False,
	"allow_multapses": False
	}


###########################################
#
# SGI connections
#

qv_qv_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
               "circular": {"radius": 0.5}},
	"delays": 5.0,
	"weights": {"gaussian": 
                       {   "sigma":.6,
			"p_center": amw*spar,
                           "anchor": [0.3, 0.0]}},
        "kernel": (.7*kernprop), #0.25

	"allow_autapses": False,
	"allow_multapses": False
	}

# qv-inhib connection for rate limitation

qv_qvinh_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
               "circular": {"radius": 0.6}},
	"delays": 1.0,
	"weights": {"gaussian": 
                       {   "sigma":.4,
                           "p_center": amw*spar, # 0.5
                           "anchor": [0.0, 0.0]}},
	"kernel": kernprop,

	"allow_autapses": False,
	"allow_multapses": False
	}

# inhib-qv feedback connection for rate limitation

qvinh_qv_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
               "circular": {"radius": 0.6}},
	"delays": 1.0,
	"weights": {"gaussian": 
                       {   "sigma":.4,		# 0.4 cheap way for flat 
                           "p_center": gbw*spar,
                           "anchor": [0.0, 0.0]}},
	"kernel": kernprop,

	"allow_autapses": False,
	"allow_multapses": False
	}

# Normal input connection to qv neurons
qv_build_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
		"circular": {"radius": .5}},
	"delays": 1.0,
	"weights": {"gaussian": 
			{   "sigma":.5,
			    "p_center": amw*0.2*spar, # orig 0.1
			    "anchor": [0.2, 0.0]}},
	"kernel": kernprop,
	"allow_autapses": False,
	"allow_multapses": False
	}


# The burst neuron connection from the input layer

qv_burst_connect = {
	    "point_a": a_angle,
	    "val_a":   nmw*a_nmda*spar, #nmw*nmdapar,
	    "point_b": b_angle,
	    "val_b":   nmw*b_nmda*spar, #nmw*nmdapar,
	    "radius":  1.5,
	    "sigma":   0.5,
	    "kernel":  kernprop,
	    "delay":   1.0,
	    "synapse": "nmda_synapse"
	    }


# Build neuron intraconnection, to create the activation wave (moving hill)
# that in turn shuts down burst neurons through the inhibitory interneurons

build_build_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
               "circular": {"radius": 0.5}},
	"delays": 5.0,
	"weights": {"gaussian": 
                       {   "sigma":.6,
			"p_center": amw*spar,
                           "anchor": [0.3, 0.0]}}, 
        "kernel": (.7*kernprop),
	"allow_autapses": False,
	"allow_multapses": False
	}

# Build neuron intraconnection, to create the activation wave (moving hill)
# that in turn shuts down burst neurons through the inhibitory interneurons

# build-inhib connection for rate limitation

build_buildinh_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
               "circular": {"radius": 0.6}},
	"delays": 1.0,
	"weights": {"gaussian": 
                       {   "sigma":.4,
                           "p_center": amw*.5*spar, # 0.5
                           "anchor": [0.0, 0.0]}},
	"kernel": kernprop,

	"allow_autapses": False,
	"allow_multapses": False
	}

# inhib-build feedback connection for rate limitation

buildinh_build_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
               "circular": {"radius": 0.6}},
	"delays": 1.0,
	"weights": {"gaussian": 
                       {   "sigma":.4,		# 0.4 cheap way for flat 
                           "p_center": gbw*3*spar,
                           "anchor": [0.0, 0.0]}},
	"kernel": kernprop,

	"allow_autapses": False,
	"allow_multapses": False
	}
# Connect build neurons to their burst counterparts

build_burst_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
		"circular": {"radius": .5}},
	"delays": 1.0,
	"weights": {"gaussian": 
			{   "sigma":.4,
			    "p_center": amw*buildpar*spar, #0.10
			    "anchor": [0.0, 0.0]}},
	"kernel": kernprop,
	"allow_autapses": False,
	"allow_multapses": False
	}

# Return connection from burst to build neurons

burst_build_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
		"circular": {"radius": 0.5}},
	"delays": 1.0,
	"weights": {"gaussian": 
			{   "sigma": 0.3, #0.3
			    "p_center": amw*2.0, # 2.0 for a buildup peak
			    "anchor": [0.0, 0.0]}},
	"kernel": kernprop,
	"allow_autapses": False,
	"allow_multapses": False
	}

#Buildup neurons to the inhibitory interlayer
# we give it a wider prjection area than is motivated by recent research.
# It's a simplification, but a necessary one
build_inhib_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
		"circular": {"radius": 2.0}}, # actually more like 0.5
	"delays": 1.0,
	"weights": {"gaussian": 
			{   "sigma":1.0,
			    "p_center": amw*fback*spar, # amw*0.10
			    "anchor": [0.0, 0.0]}},
	"kernel": .20*kernprop,
	"allow_autapses": False,
	"allow_multapses": False
	}

#  Burst neurons to the inhibitory interlayer

burst_inhib_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
		"circular": {"radius": 0.5}},
	"delays": 1.0,
	"weights": {"gaussian": 
			{   "sigma":.1,
			    "p_center": 0.1*amw*spar,
			    "anchor": [0.0, 0.0]}},
	"kernel": 0.1*kernprop,
	"allow_autapses": False,
	"allow_multapses": False
	}

# inhibitory interlayer coactivation

inhib_inhib_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
		"circular": {"radius": 0.5}},
	"delays": 1.0,
	"weights": {"gaussian": 
			{   "sigma":.4,
			    "p_center": amw*2.3*spar,    # was 2.3
			    "anchor": [0.0, 0.0]}},
	"kernel": kernprop,
	"allow_autapses": False,
	"allow_multapses": False
	}

# backprojections from the deep inhibitory layer
inhib_burst_inh_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
		"circular": {"radius": 0.5}},
	"delays": 1.0,
	"weights": {"gaussian": 
			{   "sigma":0.5,
			    "p_center": gbw*4,
 			    "anchor": [0.0, 0.0]}},
	"kernel": 1.0*kernprop,
	"allow_autapses": False,
	"allow_multapses": False
	}
inhib_qv_inh_connect = {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
		"circular": {"radius": 0.5}},
	"delays": 1.0,
	"weights": {"gaussian": 
			{   "sigma":0.5,
			    "p_center": gbw*4,
 			    "anchor": [0.0, 0.0]}},
	"kernel": kernprop,
	"allow_autapses": False,
	"allow_multapses": False
	}


###################################
#
# Create the connections themselves
#


topology.ConnectLayers(wifi.net, qv.net, wifi_qv_connect)
topology.ConnectLayers(wifi.net, build.net, wifi_build_connect)

topology.ConnectLayers(nafi.net, qv.net, nafi_qv_connect)
topology.ConnectLayers(nafi.net, build.net, nafi_build_connect)

topology.ConnectLayers(qv.net, qv.net, qv_qv_connect)
topology.ConnectLayers(qv.net, qvinh.net, qv_qvinh_connect)
topology.ConnectLayers(qvinh.net, qv.net, qvinh_qv_connect)
topology.ConnectLayers(qv.net, build.net, qv_build_connect)

burst.nmdacon(qv, pyrngs, qv_burst_connect)

topology.ConnectLayers(build.net, build.net, build_build_connect)

topology.ConnectLayers(build.net, buildinh.net, build_buildinh_connect)
topology.ConnectLayers(buildinh.net, build.net, buildinh_build_connect)

topology.ConnectLayers(build.net, burst.net, build_burst_connect)


topology.ConnectLayers(burst.net,build.net, burst_build_connect)

#***
topology.ConnectLayers(build.net, inhib.net, build_inhib_connect)
topology.ConnectLayers(burst.net, inhib.net, burst_inhib_connect)
topology.ConnectLayers(inhib.net, inhib.net, inhib_inhib_connect)


if doinhibition:
    topology.ConnectLayers(inhib.net, burst.net, inhib_burst_inh_connect)
    topology.ConnectLayers(inhib.net, qv.net, inhib_qv_inh_connect)


##################################################
##
## Synthetic integration
##

synoutw = gabaw*inhibpar

integ_burst_connect = {
	
	"point_a": a_angle,		    # degrees
	"val_a": gabaw*a_inhib,
	"point_b": b_angle,
	"val_b":  gabaw*b_inhib,
	"delay": synconspeed,
	"synapse": "static_synapse"
	}

burst.inhibcon(integ, integ_burst_connect)

nest.ConvergentConnect(nest.GetNodes(burst.net)[0], integ.blob,
	syninw*kernprop, synconinspeed, "static_synapse")

nest.ConvergentConnect(integ.blob, nest.GetNodes(buildinh.net)[0],
	synbuildinhw*1.0, synconspeed, "static_synapse")

# integration reset
if doinhibition:
    nest.ConvergentConnect(nest.GetNodes(inhib.net)[0], integ.blob, syninhinhw,
	synconinhibspeed*1.0, "static_synapse")


lpend=time.time()
if debug:
    print "Rank: %d: connection time: %.1f s" %(nest.Rank(), lpend-lpstart)

lpstart=time.time()

####################################################
#
### Synthetic input ###
#

# the pitch and stuff just sets the surface size for standalone running, but
# also specifies the layout of inputs if we use MUSIC integration

spikes = Surface(in_rows, in_cols, in_pitch, A, bx, by)


# SCS input layer

inlatency = -1.0

# This is all only when interfacing the model with remote systems. 
if domusic:
     
    nest.CopyModel("music_event_in_proxy", "spin_proxy", 
	    params={'port_name': 'SC_right_input'})
    spikes.create_layer("spin_proxy")
    
    if inlatency>0:
	nest.SetAcceptableLatency("SC_right_input", inlatency)
    sppar = spikes.save_params()
    spins = sppar["SXY"]
    spfreeze = sppar["freeze"]

    spkeys=spins.keys()				# Get UIDs
    print "inp totals #:", len(spkeys)
    spvid = list(set(spkeys)-set(spfreeze))	# exlude frozen ones
    spvalid=[[k, spins[k]] for k in spvid]	# get compound list with only
						# valid UIDs
    
    snr=sorted(spvalid)

    i=0
    inplist=[]
    for sp in snr:
	nest.SetStatus([sp[0]], {'music_channel': i})
	inplist.append([i, sp[1]])
	i+=1

    # Need to still assign a channel to these or they get a #0 default
    for sp in spfreeze:
	nest.SetStatus([sp], {'music_channel': i})
	i+=1

# Without MUSIC. 

else:
    nest.CopyModel("poisson_generator", "syn_generator", params={'rate': 0.0})
    inp = spikes.retxy_to_xy((retx,rety))
    spinput = spikes.create_free_layer("syn_generator", [inp])

# Test code for spatial input 
# 0 angle 9 degree radius, and 45 angle 8 degree radius
#inp0900 = spikes.retxy_to_xy((9.0, 0.0))
#inp0845 = spikes.retxy_to_xy((5.66, 5.66))
#(spinput1, spinput2) = spikes.create_free_layer("syn_generator", [inp0900,
#    inp0845])

# Wide field input

spikes_wifi_logpo_connect = {

	"radius": 1.5,
	"variance": 0.6,
	"weight": amw*.07,
	"shift": 0.5
	}

# Narrow field input
spikes_nafi_connect =  {
	"connection_type": "convergent",
	"synapse_model": "static_synapse",
	"mask": {
	    "circular": {"radius": .2}},
	"delays": 1.0,
	"weights": amw*.15,
	"kernel": 1.0,
	"allow_autapses": False,
	"allow_multapses": False
	}


spikes.logpolar(wifi, spikes_wifi_logpo_connect)

topology.ConnectLayers(spikes.net, nafi.net, spikes_nafi_connect)
if not domusic:
    nest.SetStatus(spinput, {"rate": spikerate*spikepop*1.0})


# synthetic SNpr inhibition

if domusic:
    
    inspikes = nest.Create("music_event_in_proxy", 1)
    nest.SetStatus(inspikes, {"port_name":"SNpr_right"})
    if inlatency>0:
        nest.SetAcceptableLatency("SNpr_right", inlatency)

    nest.ConvergentConnect(inspikes,
	    nest.GetNodes(burst.net)[0], gbw*1.0, synconspeed, "static_synapse")
    nest.ConvergentConnect(inspikes, integ.blob, syninhinhw*1.0, synconspeed, "static_synapse")

else:

    inspikes = Surface(syn_rows, syn_cols, syn_pitch)
    nest.CopyModel("poisson_generator", "insyn_generator", params={'rate': 0.0})

    inspikes.create_layer("insyn_generator")
    nest.ConvergentConnect(nest.GetNodes(inspikes.net)[0],
	    nest.GetNodes(burst.net)[0], gbw*1.0, synconspeed, "static_synapse")
    nest.ConvergentConnect(nest.GetNodes(inspikes.net)[0], integ.blob, syninhinhw*1.0, synconspeed, "static_synapse")


lpend=time.time()
if debug:
    print "Rank: %d: odds and ends: %.1f s" %(nest.Rank(), lpend-lpstart)

lpstart=time.time()
#########################
#
# Save on disk
#

def mkdir_safe(path):
    try:
	os.makedirs(path)
    except OSError, e:
#	print "errno: ", e
	if e.errno == errno.EEXIST:
#	    print "pass"
	    pass
	else:
#	    print "raise"
	    raise

if savedata==True:
    mkdir_safe(simname)

    if not os.path.isdir(simname):
	print "\n '", simname,"' exists as a file. Can't create directory of\
	the same name.\n"
	sys.exit(-1)

    savepath = simname


    nest.SetKernelStatus({'data_path': savepath, 'overwrite_files': True})

    wifidata = wifi.save_data(savename)
    nafidata = nafi.save_data(savename)
    qvdata = qv.save_data(savename)
    builddata = build.save_data(savename)
    buildinhdata = buildinh.save_data(savename)
    burstdata = burst.save_data(savename)
    inhibdata = inhib.save_data(savename)
    integdata = integ.save_data(savename)

    simdata = {}
    simdata["simtime"] = simtime
    simdata["datadir"] = savename
    simdata["simdir"] = simname
    simdata["params"] = {
	    'inhib': inhibpar,
	    'nmda': nmdapar,
	    'build': buildpar,
	    'rand': sseed,
	    'fback': fback,
	    'retx': retx,
	    'rety': rety,
	    'runs': 1}
    simdata["surfaces"] = [wifidata, nafidata, qvdata, builddata,
	    buildinhdata, burstdata, inhibdata]
    simdata["blobs"] =[integdata]
    
    fname = os.path.join(simname, savename+".sim")

    try:
        f = open(fname, "w")
	try:
	    json.dump(simdata, f, separators=(', ',':'))
	finally:
	    f.close()
    except IOError:
	
	print "Failed to save simulation data. Aborting."
	sys.exit(-1)

if domusic:

    # burst output

    itemburst = burst.save_params()["SXY"]	
    aburst = sorted(itemburst.items()) # sorted (id, (xcoord, ycoord))
    burstnr = len(aburst)
    print "burst #:", burstnr
    burst_motor = nest.Create("music_event_out_proxy",1)
    nest.SetStatus(burst_motor, {"port_name":"burst_right_motor"})
    
    proxylist=[]
    i = 0
    for ab in aburst:
	nest.Connect([ab[0]], burst_motor, {'music_channel':i})
	proxylist.append([i, ab[1]])
	i+=1

    
    # deep inhibitory output
    
    iteminhib = inhib.save_params()["SXY"]
    ainhib = iteminhib.keys()
    inhib_motor = nest.Create("music_event_out_proxy", 1)
    nest.SetStatus(inhib_motor, {"port_name":"inhib_right_motor"})

    if doinhibition:
	for inid in ainhib:
	    nest.Connect([inid], inhib_motor, {'music_channel':0})
    else:
	dummy = nest.Create("iaf_neuron", 1)
	nest.Connect(dummy, inhib_motor, {'music_channel':0})



# save parameters for MUSIC-connected modules
    if saveparams:
	
	fname = "M_"+simname+".json"
	try:
	    f = open(fname, "w")
	    try:
		json.dump({'burst': proxylist, 'spinput': inplist}, f, separators=(', ',':'))
	    finally:
		f.close()
	except IOError:
	    
	    print "Failed to save Music parameter data."
	    sys.exit(-1)


# Hacky way to delete unneeded neurons

lpend=time.time()
if debug:
    print "Rank: %d: saving: %.1f s" %(nest.Rank(), lpend-lpstart)

lpstart=time.time()

wifi.freeze()
nafi.freeze()
qv.freeze()
qvinh.freeze()
build.freeze()
buildinh.freeze()
burst.freeze()
inhib.freeze()
spikes.freeze()

lpend=time.time()
if debug:
    print "Rank: %d: freeze: %.1f s" %(nest.Rank(), lpend-lpstart)



simstart =time.time()
# Specific run for the two-target simulation used in [Moren2013]

# inhibit and stabilize
#nest.SetStatus(topology.GetElement(inspikes.net,[xin,yin]), {"rate": inhspikerate*1.0})
#nest.Simulate(100)
# start the saccade
#nest.SetStatus(topology.GetElement(inspikes.net,[xin,yin]), {"rate": inhspikerate*0.0})
#nest.Simulate(140)
# reinhibit, remove target
#nest.SetStatus(topology.GetElement(inspikes.net,[xin,yin]), {"rate": inhspikerate*1.0})
#nest.SetStatus((spinput1), {"rate": spikerate*spikepop*0.0})
#nest.Simulate(100)
# new target
#nest.SetStatus((spinput2), {"rate": spikerate*spikepop*1.0})
#nest.Simulate(50)
# start new saccade
#nest.SetStatus(topology.GetElement(inspikes.net,[xin,yin]), {"rate": inhspikerate*0.0})
#nest.Simulate(300)

#sys.exit(0)
allstart=0
if not domusic:
    nest.SetStatus(topology.GetElement(inspikes.net,[xin,yin]), {"rate": inhspikerate*1.0})
    nest.Simulate(100);
    nest.SetStatus(topology.GetElement(inspikes.net,[xin,yin]), {"rate": inhspikerate*0.0})
    print "Release"
    allstart = time.time()
    
    if True:
	nest.Simulate(simtime-100)

    else:
	for i in range(2, int((simtime)/50)):
	    stepstart = time.time()

	    if nest.Rank() == 0:
		if debug:
		    print "ms: %d" %(i*50)
	    nest.Simulate(50)
	    stepend = time.time()
	    steptime =stepend-stepstart

	    if nest.Rank() == 0:
		print "at %s-%sms: %.2f s\t %.1f times" % (i*50, (i+1)*50, steptime,
		    (steptime*1000.0)/50.0)
else:

    for i in range(0, int(simtime/50)):
	stepstart = time.time()

	nest.Simulate(50)
	stepend = time.time()
	steptime =stepend-stepstart

	if nest.Rank() == 0:
	    print "at %s-%sms: %.2f s\t %.1f times" % (i*50, (i+1)*50, steptime,
		(steptime*1000.0)/50.0)


simend=time.time()


if nest.Rank() == 0:
    print "Simulation time: %.1f s \t %.1f times" % (simend-simstart,
	    (simend-simstart)*1000.0/simtime)
    print "\tactive time: %.1f s \t %.1f times" % (simend-allstart,
	    (simend-allstart)*1000.0/float(simtime-100))

    


