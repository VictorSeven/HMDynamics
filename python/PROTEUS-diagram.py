#!/usr/bin/python

# --- Header
import numpy as np 
import sys
import os

# --- Simulation parameters

delta, s = 0.5, 0.4
q_values = [0.0, 7.0, 100]
a_values = [0.0, 0.8]
w = 1.0

networks = ["erdos-renyi", "random", "core"]
n_averages = 50 

fileout = "../data/phase-diagram-coupling/"

# --- Get absolute path

path_2_this = os.path.dirname(__file__)
path_2_this = path_2_this if path_2_this != "" else "."

# --- C++ compilation and preparation 

datafolder = path_2_this + "/" + fileout 

#For compiling, also defining some important paths
cppfolder = path_2_this + "/../cpp/"
cppfile = cppfolder + "network-dynamics.cpp"
cppoutput = cppfolder + "bin/dynamics.exe"
netfolder = path_2_this + "/../networks/"

gcc_flags = "-std=c++11 -O3 -DMODE=DIAGRAM "

#Ensure we have folders for the stuff
os.system("mkdir {output_path}".format(output_path = cppfolder + "bin/"))
os.system("mkdir {output_path}".format(output_path = datafolder)) 

#Compile latest version of C++ dynamics 
os.system("g++ {file} {flags} -o {out}".format(file=cppfile, flags=gcc_flags, out=cppoutput))

print("Compilation finished")

# --- Prepare functions to launch all the programs

def launch_runs(params, extension, var_space):
    params["var"] = var_space 

    for network in networks:
        for j in range(n_averages):
            netpath = netfolder + network 
            outpath = datafolder + network + extension +  "_{0}".format(j)

            os.system("slanzarv --nomail --short -J {procname} {exe} {w0} {delta} {s} {variable_a} {fixed} {var[0]} {var[1]} {var[2]} {netpath} {outpath} {index}".format(**params, procname="pd"+extension+"-"+network, exe=cppoutput, index=j, netpath=netpath, outpath=outpath))

# --- Run dynamics for each network 

for a in a_values:
    params = {"w0":w, "variable_a":0, "fixed":a, "delta":delta, "s":s}
    launch_runs(params, "a{0:.2f}".format(a), q_values)