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

import os
import sys
import getopt

neurons = ['buildinhup_neuron', 'buildup_neuron', 'burst_neuron',
	'inhib_neuron', 'integrate', 'nafi_neuron', 'qv_neuron',
	'wifi_neuron']


###############################################
#
# Read input
#

appname = os.path.splitext(sys.argv[0])[0]

def show_help():

    print """
Usage: %s [-h] simname

A simple script to collate the per-surface outputs 
from an SC simulation, and sort the result by spike
time-stamp.

Useful when you want to analyze the data with external 
tools such as R, or study the output directly.

  -h, --help            show this help text
    """ % (appname)


argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "h", ["help"])
except:
    show_help()
    sys.exit(-1)

for opt, arg in opts:
    if opt in ("-h", "--help"):
        show_help()
        sys.exit(0)

if len(args) != 1:
    show_help()
    print "The script takes the simulation output name as argument.\n"
    sys.exit(-1)


##################################################
#
# Collate and sort
#

basefile = args[0]
for nt in neurons:
    cmd = "cat %s_%s* | sort -n -k 2 >%s.%s.spk"%(basefile, nt, nt, basefile)
    os.system(cmd)







