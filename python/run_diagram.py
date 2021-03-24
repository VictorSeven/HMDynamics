#!/usr/bin/python

# --- Header
import numpy as np 
import sys
import os

# --- Get absolute path

path_2_this = os.path.dirname(os.path.abspath(__file__))
path_2_this = path_2_this if path_2_this != "" else "."

# --- Check if we have a list of networks where we want to run, or if we want all

if (len(sys.argv) < 2):
    use_all = True
else:
    use_all = False

# --- C++ compilation and preparation 

is_proteus = False 

#Set an adequate launch command depending on system
if is_proteus:
    launch = "slanzarv --nomail -c={ncores} ".format(ncores=ncores)
else:
    launch = " "

#For compiling, also defining some important paths
datafolder = path_2_this + "/../data/pd/"
cppfolder = path_2_this + "/../cpp/"
cppfile = cppfolder + "network-dynamics.cpp"
cppoutput = cppfolder + "bin/dynamics.exe"
netfolder = path_2_this + "/../networks/"

gcc_flags = "-std=c++11 -O3 -DMODE=DIAGRAM" 

#Ensure we have folders for the stuff
os.system("mkdir {output_path}".format(output_path = cppfolder + "bin/"))
os.system("mkdir {output_path}".format(output_path = datafolder)) 

#Compile latest version of C++ dynamics 
os.system("g++ {file} {flags} -o {out}".format(file=cppfile, flags=gcc_flags, out=cppoutput))

print("Compilation successful")

# --- Prepare functions to launch all the programs

#Get list of files
if use_all:
    networks_files = os.listdir(netfolder)
else:
    networks_files = sys.argv[1:]

def launch_runs(params):
    for network in networks_files:
        if network.endswith(".mtx"):
            network = network[:-4]
            netpath = netfolder + network
            outpath = datafolder + network
            os.system("{launch}{exe} {w0} {delta} {q} {variable_a} {fixed} {var[0]} {var[1]} {var[2]} {netpath} {outpath}".format(**params, launch=launch, exe=cppoutput, netpath=netpath, outpath=outpath))


# --- Run dynamics for each network 

#Fixed = a
params = {"w0":1.0, "fixed":0.0, "delta":0.0,  "q":1.0,  "var":[0.0,2.0,100], "variable_a": 0}
launch_runs(params)

#Fixed = s
params = {"w0":1.0, "fixed":0.0, "delta":0.0,  "q":1.0,  "var":[0.0,2.0,100], "variable_a": 1}
launch_runs(params)

