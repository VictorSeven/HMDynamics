#!/usr/bin/python

# --- Header
import numpy as np 
import sys
from os import system

# --- Check we have a list of networks

if (len(sys.argv) < 2):
    print("Add the list of filenames you want to process. One filename per network")
    exit()

# --- C++ compilation and preparation 

is_proteus = False
ncores=4

#Set an adequate launch command depending on system
if is_proteus:
    launch = "slanzarv --nomail -c={ncores} ".format(ncores=ncores)
else:
    launch = "./"

#For compiling, also defining some important paths
cppfolder = "../cpp/"
cppfile = cppfolder + "network-dynamics.cpp"
cppoutput = cppfolder + "bin/dynamics.exe"
netfolder = "../networks/"
datafolder = "../../data/pd/"

gcc_flags = "-std=c++11 -O3 -fopenmp -DMODE=DIAGRAM -DNUM_THREADS={ncores}".format(ncores)

#Ensure we have folders for the stuff
system("mkdir {output_path}".format(output_path = cppfolder + "bin/"))
system("mkdir {output_path}".format(output_path = datafolder)) 

#Compile latest version of C++ network generator
system("g++ {file} {flags} -o {out}".format(file=cppfile, flags=gcc_flags, out=cppoutput))

print("Compilation successful")

# --- Run dynamics for each network 

params = {"w0":1.0,  "delta":0.5,  "sigma":0.0,  "q":[0.0,2.0,10]}

for network in sys.argv[1:]:
    netpath = netfolder + network
    outpath = datafolder + network
    system("{launch}{exe} {w0} {delta} {sigma} {q[0]} {q[1]} {q[2]} {netpath} {outpath}".format(**params, launch=launch, exe=cppoutput, netpath=netpath, outpath=outpath))

