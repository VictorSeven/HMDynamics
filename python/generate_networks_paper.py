#!/usr/bin/python

# --- Header
import numpy as np 
from os import system

# --- Aux functions

def list2str(lista):
    return [str(elem) for elem in lista]

# --- C++ compilation and preparation 

is_proteus = True 

#Set an adequate launch command depending on system
if is_proteus:
    launch = "slanzarv --nomail --short -J create-networks "
else:
    launch = "./"

#For compiling, also defining some important paths
cppfolder = "../cpp/"
cppfile = cppfolder + "create-network.cpp"
cppoutput = cppfolder + "bin/network.exe"
netfolder = "../networks/"

gcc_flags = "-std=c++11 -O3"

#Ensure we have folders for the stuff
system("mkdir {output_path}".format(output_path = cppfolder + "bin/"))
system("mkdir {output_path}".format(output_path = netfolder)) 

#Compile latest version of C++ network generator
system("g++ {file} {flags} -o {out}".format(file=cppfile, flags=gcc_flags, out=cppoutput))

print("Compilation successful")

# --- Network creation

# - Single modulo and Erdos-Renyi 

def generate_er_network(n,k,name):
    n,k = (n, k) #number of neurons and average degree
    netfile = netfolder + name #Where to save the result

    system("{launch}{exe} --random 1 {n} {k} {save}".format(launch=launch, exe=cppoutput, n=n, k=k, save=netfile))

generate_er_network(100,  12.0, "single-module")  #Complete ER with same connectivity as modulus
generate_er_network(6400, 30.0, "erdos-renyi")   #Complete ER with same connectivity as modulus

print("Erdos-Renyi OK")

# - Binary HM Random
depth = 6

n = [100] + [2 for j in range(depth)] #Number of modules in each level, from bottom to top

#Connectivity of levels, from bottom to top
k0, kf = 5, 0.1  
dk = (kf - k0) / depth
k = [12] + [k0 + dk*j for j in range(depth)]

#Convert lists to string
n = list2str(n)  
k = list2str(k) 
netfile = netfolder + "rb{d}-bist".format(d=depth)

#Create network
system("{launch}{exe} --random {nlv} {n} {k} {save}".format(nlv=depth+1, launch=launch, exe=cppoutput, n=" ".join(n), k=" ".join(k), save=netfile))

print("Binary HMRandom OK")

# - Binary HM Core 
depth = 6

#Same procedure as before
n = [100] + [2 for j in range(depth)]
k0, kf = 5, 0.1
dk = (kf - k0) / depth
k = [12] + [k0 + dk*j for j in range(depth)]
g = [0] + [2 for j in range(depth)] #Scale-free exponent of each module (0 means random)

#Convert lists to string
n = list2str(n)  
k = list2str(k) 
g = list2str(g)
netfile = netfolder + "cb{d}-bist".format(d=depth)

#Create network
system("{launch}{exe} --core {nlv} {n} {k} {g} {save}".format(nlv=depth+1, launch=launch, exe=cppoutput, n=" ".join(n), k=" ".join(k), g=" ".join(g), save=netfile))

print("Binary HMCore OK")
