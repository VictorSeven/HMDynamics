#!/usr/bin/python

# --- Header
import numpy as np 
import sys
import os

# --- Check if we have a list of networks where we want to run, or if we want all

if (len(sys.argv) < 2):
    use_all = True
else:
    use_all = False

# --- C++ compilation and preparation 

is_proteus = True 
ncores=4

#Set an adequate launch command depending on system
if is_proteus:
    launch = "slanzarv --nomail -c={ncores} ".format(ncores=ncores)
    datafolder = "../data/pd/"
else:
    launch = "./"
    datafolder = "../../data/pd/"

#For compiling, also defining some important paths
cppfolder = "../cpp/"
cppfile = cppfolder + "network-dynamics.cpp"
cppoutput = cppfolder + "bin/dynamics.exe"
netfolder = "../networks/"

gcc_flags = "-std=c++11 -O3 -fopenmp -DMODE=DIAGRAM -DNUM_THREADS={0}".format(ncores)

#Ensure we have folders for the stuff
os.system("mkdir {output_path}".format(output_path = cppfolder + "bin/"))
os.system("mkdir {output_path}".format(output_path = datafolder)) 

#Compile latest version of C++ dynamics 
os.system("g++ {file} {flags} -o {out}".format(file=cppfile, flags=gcc_flags, out=cppoutput))

print("Compilation successful")

# --- Run dynamics for each network 

params = {"w0":1.0,  "delta":0.5,  "sigma":0.0,  "q":[0.0,2.0,100]}

#Get list of files
if use_all:
    networks_files = os.listdir(netfolder)
else:
    networks_files = sys.argv[1:]

for network in networks_files:
    if network.endswith(".mtx"):
        netpath = netfolder + network
        outpath = datafolder + network
        os.system("{launch}{exe} {w0} {delta} {sigma} {q[0]} {q[1]} {q[2]} {netpath} {outpath}".format(**params, launch=launch, exe=cppoutput, netpath=netpath, outpath=outpath))

