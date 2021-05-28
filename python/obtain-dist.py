import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmaps
import networkx as nx
import pandas as pd
import subprocess
import bct



def functional_dist(filename, savefig=False, titlefig=None, itsperwindow=2000, dt=0.01, activity_function=lambda z: 0.5*np.abs(z)*(1+np.sin(np.angle(z))), cut=0, mode="spectral"):
    """
    filename: where the data is stored
    n_com: array that says how much communities are inside a given one
    ntracies: number of time traces recorded
    dt: to set the correct temporal scale of the X axis, timestep used in simulation
    ordered: if True, performs a communnity detection algorithm, sorting the displayed matrix accordingly
    """
    
    #Load data and get its properties
    data = pd.read_csv(filename, delimiter = " ").values

    
    x = data[cut:, :-1:2]
    y = data[cut:, 1:-1:2]
    
    z = x + 1.0j * y
    del data, x, y
    
    tsteps = np.size(z[:, 0])
    n_modules = np.size(z[0, :])
    
    activity = activity_function(z)
    
    nwindow = tsteps // itsperwindow
    t = np.linspace(0, tsteps*dt, nwindow-1)

    if mode=="spectral":
        #Define time series correlation matrix
        spectrum = np.empty((nwindow-1, n_modules))

        for i in range(nwindow-1):
            spectrum[i,:] = np.sort(np.linalg.eigvalsh(np.corrcoef(activity[i*itsperwindow:(i+1)*itsperwindow, :], rowvar=False)))

        pairs = np.empty((nwindow-1, nwindow-1))
        for i in range(nwindow-1):
            pairs[i,i] = 0.0
            for j in range(i+1,nwindow-1):
                pairs[i,j] = np.linalg.norm(spectrum[i,:] - spectrum[j,:]) 
                pairs[j,i] = pairs[i,j]
                
    elif mode=="frobenius":
        
        pairs = np.empty((nwindow-1, nwindow-1))
        for i in range(nwindow-1):
            matrix_i = np.corrcoef(activity[i*itsperwindow:(i+1)*itsperwindow, :], rowvar=False)
            pairs[i,i] = 0.0
            for j in range(i+1,nwindow-1):
                matrix_j = np.corrcoef(activity[j*itsperwindow:(j+1)*itsperwindow, :], rowvar=False)

                pairs[i,j] = np.sqrt(np.sum((matrix_i - matrix_j)**2)) 
                pairs[j,i] = pairs[i,j]
                
    elif mode=="induced":
        
        pairs = np.empty((nwindow-1, nwindow-1))
        for i in range(nwindow-1):
            matrix_i = np.corrcoef(activity[i*itsperwindow:(i+1)*itsperwindow, :], rowvar=False)
            pairs[i,i] = 0.0
            for j in range(i+1,nwindow-1):
                matrix_j = np.corrcoef(activity[j*itsperwindow:(j+1)*itsperwindow, :], rowvar=False)
                pairs[i,j] = np.max(np.linalg.svd(matrix_j - matrix_i)[1])
                pairs[j,i] = pairs[i,j]           

    return pairs




def shannon_entropy(dist, w):
    dist = dist[dist>0]
    return -np.sum(dist*np.log(dist/w))


path = "./../data/functional-data-pablo/"
phases = ["_super", "_crit", "_sub"]
phasename = ["Async/Exc", "Critical", "Synch"]

networks = ["er-cb6", "rb6", "cb6"]
netsize  = [6400, 6400, 6400]
dynamics = ["-hopf", "-hopf_exc", "-snic", "-hyb"]
dyncontrol = [r"$\Delta$", r"$\Delta$", r"$a$", r"$a$"]

nbins = 300
bins = np.linspace(0, 60, nbins)

rownames = ["Erdos-Renyi", "Random Binary 6", "Core Binary 6"]
colnames = ["Kur. (noise+het)", "Exc. Hopf", "SNIC", "Hybrid"]
fig, axes = plt.subplots(figsize=(20,30), ncols=len(networks), nrows=len(dynamics))


for i,net in enumerate(networks):
    print(net)
    for j,dyn in enumerate(dynamics):
        name = net + dyn
        data = np.empty((0,3))
        
        for k,phase in enumerate(phases):
            pairs = np.ravel(functional_dist(path+name+phase, itsperwindow=2000, cut=0, mode="induced"))
            
            histograma, edges = np.histogram(pairs, bins=bins, density=True)

            axes[j,i].fill_between(edges[:-1], np.zeros(np.size(histograma)), histograma, alpha=0.5,  label=r"{0}, $S(P_{{ij}})={1:.2f}$".format(phasename[k], shannon_entropy(histograma/np.sum(histograma), edges[1]-edges[0])))
   
        axes[j,i].legend(fontsize=8, loc="best")
        axes[j,i].set_xlim(0,60)
        axes[j,i].set_xlabel("Spect. dist. distr.")

        
for ax, col in zip(axes[0], rownames):
    ax.set_title(col, fontsize=20)

for ax, row in zip(axes[:,0], colnames):
    ax.set_ylabel(row, fontsize=20)

fig.tight_layout()
plt.savefig("functional_dist_induced_axes_pablo.pdf")
plt.show()
