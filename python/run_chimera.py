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

#For compiling, also defining some important paths
datafolder = path_2_this + "/../data/chimera/"
cppfolder = path_2_this + "/../cpp/"
cppfile = cppfolder + "network-dynamics.cpp"
cppoutput = cppfolder + "bin/chimera.exe"
netfolder = path_2_this + "/../networks/"

gcc_flags = "-std=c++11 -O3 -DMODE=CHIMERA"

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

#Automatically launch all possible networks in the folder
def launch_runs(params, extension):
    network_moduli = {"erdos-renyi":10, "r4":30, "c4": 30, "rb7":128, "cb7":128}
    for network in networks_files:
        #Filter just the ones I can read
        if network.endswith(".mtx"):
            #Get the number of moduli of this network (known by construction, check generate_networks)
            network = network[:-4]
            params["nmoduli"] = network_moduli[network] 

            #Get all paths and then run program
            netpath = netfolder + network
            outpath = datafolder + network + '-' + extension
            os.system("slanzarv --nomail -J {procname} {exe} {w0} {delta} {s} {a} {q} {nmoduli} {netpath} {outpath}".format(**params, procname="trace_"+network, exe=cppoutput, netpath=netpath, outpath=outpath))

# --- Run dynamics for each network 

a_list = [0.0, 0.0,  0.0,  0.5, 0.5, 0.5,  0.9,   1.1,   1.07,  0.9,   1.1,   1.0]
d_list = [0.3, 0.45, 0.55, 0.2, 0.3, 0.4,  0.125, 0.125, 0.125, 0.005, 0.005, 0.005] 
#s_list = [0.8, 1.0, 1.2, 0.8, 0.9, 1.0, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1]
name_list = ["hopf_sub", "hopf_crit", "hopf_super", "hopf_exc_sub", "hopf_exc_crit", "hopf_exc_super", "hyb_sub", "hyb_crit", "hyb_super", "snic_sub", "snic_crit", "snic_super"]

for a,delta,name in zip(a_list, d_list, name_list):
    params = {"w0":1.0, "a":a, "delta":delta, "s":0.1, "q":1.0}
    launch_runs(params, name)



