#!/usr/bin/python

import numpy as np 
import matplotlib.pyplot as plt 
#from matplotlib import cm
from scipy.integrate import odeint
import figs_header


x = np.linspace(0,7)
y = np.sin(x)

def bogdanov_takens(xt, t, a, b):
    x,y = xt 

    xdot = y 
    ydot = a + b*x  + x*x - x*y 

    return (xdot, ydot)

t = np.linspace(0.0, 100.0, 1000)


param_set = [(0.05, -0.5)]

x0_set = [(0.01, 0.0), (-0.03, 0.0), (0.0, 0.01), (0.0, -0.01)]

colormap = plt.get_cmap("YlGn")
ncolors = 5
colormap = [colormap(0.3 + j/ncolors) for j in range(ncolors)]

figs_header.set_up_figure()

fig, axes = plt.subplots(ncols=4, nrows=2, constrained_layout=True)
for i,(a,b) in enumerate(param_set):
    fp = 0.5 * (-b + np.sqrt(b*b - 4*a))
    print(i)
    for j in range(4):
        x0 = x0_set[j]
        trajectory = odeint(bogdanov_takens, x0, t, args=(a,b))
        axes[i,0].plot(trajectory[:,0], trajectory[:,1], color = colormap[j])

plt.show()