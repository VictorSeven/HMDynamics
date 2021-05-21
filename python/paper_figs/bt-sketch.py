#!/usr/bin/python

import numpy as np 
import matplotlib.pyplot as plt 
#from matplotlib import cm
from scipy.integrate import solve_ivp
import figs_header


x = np.linspace(0,7)
y = np.sin(x)

#Strogatz and Childs, 2008, so we have Ott-Antonsen parameters
def bogdanov_takens(t, xt, a, b):
    r,p = xt 
    r2 = r*r

    rdot = 0.5*5*r*((1.0-r2)) - r + 0.5*a*(1.0-r2)*np.cos(p)
    pdot = -(b + 0.5*a*(r+1.0/r)*np.sin(p))

    return (rdot, pdot)

#Get the time for the solutions
t = np.linspace(0.0, 70.0, 1000)

#Parameters for frequency and excitability
param_set = [(2.96, 3.06), (3, 3.06),   (2.96, 3.1), (2.99, 3.09)]
titles =    ["Excit.",   "Asynch.",     "Sync.",     "Hybrid"]

#Initial conditions
x0_set = [(-0.1, -0.8), (0.0, -0.5), (0.1,-0.1), (0.3, 0.2), (-0.2, -0.8), (-0.15, -0.9), (-0.3, 0.3), (-0.5, 0.6), (-0.4, -0.6)]

#Set colormap per graph
colormaps = ["YlGn", "YlOrBr", "OrRd", "PuBu"]
init_color = 0.3
ncolors = len(x0_set)
palettes = {} 
for cm in colormaps:
    colormap = plt.get_cmap(cm)
    palettes[cm] = [colormap(init_color + j/ncolors) for j in range(ncolors)]


#Set fonts and size
figs_header.set_up_figure()
fig = plt.figure(figsize=figs_header.figsize())


#Prepare all polar axes in the indicated positions
axes = []
marginx, marginy = 0.5, 0.1
axsize = 0.32
sepx, sepy = axsize - 0.1, axsize + 0.15
for i in range(2):
    for j in range(2):
        x, y = marginx + i*sepx, marginy + j*sepy
        axes.append(fig.add_axes([x,y,axsize, axsize], projection="polar"))

#Solve for each one of the parameter sets
for i,(a,b) in enumerate(param_set):

    cm = palettes[colormaps[i]]

    #Different initial conditions
    for j in range(ncolors):
        x0 = x0_set[j] 
        x0 = [np.sqrt(x0[0]**2 + x0[1]**2), np.arctan2(x0[1], x0[0])] #To polar

        #Solve problem
        trajectory = solve_ivp(bogdanov_takens, [t[0], t[-1]], x0, args=(a,b), dense_output=True, method="LSODA")
        trajectory = trajectory.sol(t).T
        
        axes[i].plot(trajectory[:,1], trajectory[:,0], color=cm[j])
    

    #Eliminate tick labels
    axes[i].set_xticklabels([])
    axes[i].set_yticklabels([])
    #Set title
    axes[i].set_title(titles[i], loc="left", color=cm[2*len(cm)//3])

#Last steps
plt.savefig("sketch-phase.svg")
plt.show()