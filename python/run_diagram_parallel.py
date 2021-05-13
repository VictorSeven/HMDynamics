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
    n_simulations = npoints // points_per_file
    var = np.linspace(var0, varf, n_simulations)
    for network in networks_files:
        if network.endswith(".mtx"):
            network = network[:-4]
            netpath = netfolder + network
            for i in range(n_simulations-1):
                params["var"] = [var[i], var[i+1], points_per_file] 
                outpath = datafolder + network + extension + "_part{0}".format(i)
                os.system("slanzarv --nomail -J {procname} {exe} {w0} {sigma} {q} {variable_a} {fixed} {var[0]} {var[1]} {var[2]} {netpath} {outpath}".format(**params, procname="pd"+extension+"-"+network, exe=cppoutput, netpath=netpath, outpath=outpath))


# --- Run dynamics for each network 

#Kuramoto model simulation
to_detrm_noise = lambda x: 0.5*x*x

total_noise = np.array([0.5, 1.5])
stocht_noise = 0.0
determ_noise = total_noise - stocht_noise 
determ_noise = to_detrm_noise(determ_noise)

params = {"w0":1.0, "variable_a": 0, "fixed":0.0, "sigma":stocht_noise, "q":1.0}
launch_runs(params, "kuramoto", [determ_noise[0], determ_noise[1], 100])

#Fixed = a
a_list = [0.0, 0.5] 

total_noise = np.array([0.5, 1.5])
stocht_noise = 0.1
determ_noise = total_noise - stocht_noise 
determ_noise = to_detrm_noise(determ_noise)

filenames = ["-non-excitable", "-exc_hopf"]
for a,outname in zip(a_list, filenames):
    params = {"w0":1.0, "variable_a": 0, "fixed":a, "sigma":stocht_noise,  "q":1.0}
    launch_runs(params, outname, [determ_noise[0], determ_noise[1], 100])

#Fixed = s
total_noise = np.array([0.2, 0.5])
stocht_noise = 0.1
determ_noise = total_noise - stocht_noise 
determ_noise = to_detrm_noise(determ_noise)

filenames = ["-snic", "-hybrid"]
for delta,outname in zip(determ_noise, filenames):
    params = {"w0":1.0, "variable_a":1, "fixed":delta, "sigma":stocht_noise,  "q":1.0}
    launch_runs(params, outname, [0.5,1.5,100])
