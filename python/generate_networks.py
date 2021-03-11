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
    launch = "slanzarv --nomail "
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

# - Erdos-Renyi 
n,k = (1000, 10.0) #number of neurons and average degree
netfile = netfolder + "erdos-renyi" #Where to save the result

system("{launch}{exe} --random 1 {n} {k} {save}".format(launch=launch, exe=cppoutput, n=n, k=k, save=netfile))

# - HM Random
n = [100, 5, 3, 2]          #Number of elements in each module
k = [10.0, 5.0, 1.0, 0.5]    #Connectivity of each module

#Convert lists to string
n = list2str(n)  
k = list2str(k) 
netfile = netfolder + "r4"

system("{launch}{exe} --random 4 {n} {k} {save}".format(launch=launch, exe=cppoutput, n=" ".join(n), k=" ".join(k), save=netfile))

print("HMRandom OK")

# - Binary HM Random
depth = 7

n = [30] + [2 for j in range(depth)]
k = [20] + [depth-j for j in range(depth)]

#Convert lists to string
n = list2str(n)  
k = list2str(k) 
netfile = netfolder + "rb{d}".format(d=depth)

system("{launch}{exe} --random {nlv} {n} {k} {save}".format(nlv=depth+1, launch=launch, exe=cppoutput, n=" ".join(n), k=" ".join(k), save=netfile))

print("Binary HMRandom OK")

# - HM Core

n = [100, 5, 3, 2]          #Number of elements in each module
k = [10.0, 5.0, 1.0, 0.5]   #Connectivity of each module
g = [0.0, 2.0, 2.0, 2.0]    #Scale free exponents 

#Convert lists to string
n = list2str(n)  
k = list2str(k) 
g = list2str(g)
netfile = netfolder + "c4"

system("{launch}{exe} --core 4 {n} {k} {g} {save}".format(launch=launch, exe=cppoutput, n=" ".join(n), g=" ".join(g), k=" ".join(k), save=netfile))

print("HMCore OK")

# - Binary HM Core 
depth = 7

n = [30] + [2 for j in range(depth)]
k = [15] + [depth-j for j in range(depth)]
g = [0] + [2 for j in range(depth)]

#Convert lists to string
n = list2str(n)  
k = list2str(k) 
g = list2str(g)

netfile = netfolder + "cb{d}".format(d=depth)

system("{launch}{exe} --core {nlv} {n} {k} {g} {save}".format(nlv=depth+1, launch=launch, exe=cppoutput, n=" ".join(n), k=" ".join(k), g=" ".join(g), save=netfile))

print("Binary HMRandom OK")
