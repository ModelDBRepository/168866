#!/usr/bin/env python
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
# Reformat and compact the spike data.
#
# * Collate data into one file per surface, regardless of the number of nodes
# used in simulation;
#
# * use a more efficient savings format (one list of spike time per neuron,
# rather than neuron-time pairs);
#
# * Use the binary python 'dump' format rather than text for the spike data
#   and compress with bzip2
#
# * sort spike files into subdirectories if they're not already there.
#


import math
import os.path
import sys
import time
import getopt
import re
import pickle
import bz2
import json
from SGIdataclass import *


###############################################
#
# Parameters
#

dearchive = False	    # Don't read a tgz archive
del_orig = False

###############################################
#
# Read input
#

appname = os.path.splitext(sys.argv[0])[0]

def show_help():

    print """
Usage: %s [-dh] input[.sim]

Convert the results files of an SC simulation into a space-efficient packed
binary format. Good for archiving and for writing efficient tools that need to
deal with a lot of simulation results (collating multiple simulations for
instance). 

For a per-simulation analysis it is likely easier to use the "assemble.py"
tool and simple text files instead.

  -h, --help            show this help text
  -d, --delete		delete the original files after conversion
    """ % (appname)


argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "hd", ["help", "delete"])
except:
    show_help()
    sys.exit(-1)

for opt, arg in opts:
    if opt in ("-h", "--help"):
	show_help()
	sys.exit(0)
    if opt in ("-d", "--delete"):
	del_orig = True

if len(args) != 1:
    show_help()
    print "The script takes only one file argument"
    sys.exit(-1)


# Deal with the simulation header file

(fhead, fbase, ext) = simp =  find_sim(args[0])

print "simfile: ", fbase+ext


delfiles = [os.path.join(fhead, fbase+".sim")]
simdata = load_simfile(simp)
# Fix subdirectories

subdir = simdata['datadir']
if not os.path.exists(subdir):
    try:
        os.makedirs(subdir)
	pass
    except OSError, e:
	print "errno: ", e
        if e.errno == errno.EEXIST:
            pass
        else:
            raise
    
    for infile in glob.glob(os.path.join(fhead, subdir+'*.spikes')):
	rp = re.match(".*"+subdir+"_", infile)
	dest = infile[rp.end():]
	os.rename(os.path.join(fhead,infile), os.path.join(fhead, subdir,
	    dest))

# save in compressed form
f = bz2.BZ2File(os.path.join(fhead, fbase+'.zim'), 'w')

try:
    json.dump(simdata,f)
except:
    print "Failed to save simulation parameters in ", fbase+".zim","\n"
    sys.exit(-1)

f.close()

# Load the spike data
for surface in simdata["surfaces"]:
 
    print "    Unit: ", surface["name"]
    spiked = load_spikefile(simp, surface["name"])
    
    if ext == ".sim":
        for infile in glob.glob(os.path.join(fhead, fbase, surface["name"]+'*.spikes')):
	    delfiles.append(infile)
	
    outfile = os.path.join(fhead, fbase, surface["name"]+'.zpikes')
    f = bz2.BZ2File(outfile, 'w')
    try:
	pickle.dump(spiked, f, -1)
    except:
	print "Failed to save simulation data ", outfile,"\n"
	sys.exit(-1)
    f.close()

for blob in simdata["blobs"]:

    print "    Unit: ", blob["name"]
    spiked = load_blobfile(simp, blob["name"])
    
    if ext == ".sim":
        for infile in glob.glob(os.path.join(fhead, fbase, blob["name"]+'*.spikes')):
	    delfiles.append(infile)
	
    outfile = os.path.join(fhead, fbase, blob["name"]+'.zpikes')
    f = bz2.BZ2File(outfile, 'w')
    try:
	pickle.dump(spiked, f, -1)
    except:
	print "Failed to save simulation data ", outfile,"\n"
	sys.exit(-1)
    f.close()
	
if del_orig == True:
    print "Deleting..."
    for f in delfiles:
	os.remove(f)


