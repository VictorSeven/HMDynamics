#!/usr/bin/python

# --- Header
import numpy as np 
import sys
import os

# --- Get absolute path

path_2_this = os.path.dirname(__file__)
path_2_this = path_2_this if path_2_this != "" else "."

# --- Check if we have a list of networks where we want to run, or if we want all

if (len(sys.argv) < 2):
    use_all = True
else:
    use_all = False

# --- C++ compilation and preparation 

datafolder = path_2_this + "/../data/pd/"

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

print("Compilation successful")

# --- Prepare functions to launch all the programs

#Get list of files
if use_all:
    networks_files = os.listdir(netfolder)
else:
    networks_files = sys.argv[1:]


def launch_runs(params, extension, var_space, points_per_file=10):
    var0, varf, npoints = var_space 
    n_simulations = int(npoints // points_per_file)
    var = np.linspace(var0, varf, n_simulations)

    #Read all the networks
    for network in networks_files:
        if network.endswith(".mtx"):
            network = network[:-4] #Eliminate .mtx from name
            netpath = netfolder + network

            #Launch simulations
            for i in range(n_simulations-1):
                params["var"] = [var[i], var[i+1], points_per_file] 
                outpath = datafolder + network + extension + "_part{0}".format(i)
                os.system("slanzarv --nomail -J {procname} {exe} {w0} {delta} {q} {variable_a} {fixed} {var[0]} {var[1]} {var[2]} {netpath} {outpath}".format(**params, procname="pd"+extension+"-"+network, exe=cppoutput, netpath=netpath, outpath=outpath))


# --- Run dynamics for each network 

#Stochastic Kuramoto model simulation
stocht_noise = np.array([0.5, 1.5, 100])
determ_noise = 0.0

params = {"w0":1.0, "variable_a":0, "fixed":0.0, "delta":determ_noise, "q":1.0}
launch_runs(params, "kuramoto", stocht_noise)

#Fixed = a
a_list = [0.0, 0.5, 0.8, 0.92, 1.0, 1.1] 

stocht_noise = np.array([0.3, 1.4, 100])
determ_noise = 0.03

filenames = ["a_{0:.2f}".format(a) for a in a_list]
#filenames = ["-non-excitable", "-exc-hopf"]
for a,outname in zip(a_list, filenames):
    params = {"w0":1.0, "variable_a": 0, "fixed":a, "delta":determ_noise,  "q":1.0}
    launch_runs(params, outname, stocht_noise)

#Fixed = s
stocht_noise = [0.1, 0.3, 0.5, 0.6, 0.8] 
determ_noise = 0.03 

a_list = np.array([0.5,1.5,100])

filenames = ["s_{0:.2f}".format(s) for s in stocht_noise]
#filenames = ["-snic", "-hybrid"]
for s,outname in zip(stocht_noise, filenames):
    params = {"w0":1.0, "variable_a":1, "fixed":s, "delta":determ_noise,  "q":1.0}
    launch_runs(params, outname, a_list)
