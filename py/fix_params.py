#!/usr/bin/env python2.7

# imports 
import argparse
import sys
import numpy as np


# command-line arguments
parser = argparse.ArgumentParser()

parser.add_argument("-temp","--temperature", type=float, help="temparature of the gas (default 1.e4 K)",default="1.e4")

parser.add_argument("-tau","--tau", type=float, help="line centre optical depth (default 1.e5)",default="1.e5")

parser.add_argument("-box","--boxsize", type=float, help="box size in cm (default 1.e18)",default="1.e18")

parser.add_argument("-r","--radius", type=float, help="sphere radius in box unit (default 0.5)",default="0.5")

if len(sys.argv)<2:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

temperature = args.temperature
tau         = args.tau
boxsize     = args.boxsize
radius      = args.radius

print "<--- input parameters <---"
print "T            =",temperature
print "tau          =",tau
print "boxsize [cm] =",boxsize
print "radius       =",radius

print " "
print "---> fix parameters --->"

R_cm        = radius * boxsize
sigma_0     = 5.88e-14 * (temperature/1.e4)**(-0.5) # cm^2 from Dijkstra14
density     = tau/sigma_0/R_cm
vth_new     = 12.9*np.sqrt(temperature/1.e4)*1.e5

print "# overwrite parameters"
print "  gas_overwrite       = T"
print "  fix_nhi             = ",density
print "  fix_vth             = ",vth_new
print "  fix_box_size_cm     = ",boxsize
  
