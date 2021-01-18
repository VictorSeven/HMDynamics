# HM Dynamics

Here we study the dynamics of hierarchical modular neuronal networks. 

## Directory structure

- `cpp`: contains C++ code to make simulations.
- `python`: contains Python code for launching programs and generating graphs

## Compiling sources

For C++ programs, it is enough to compile as always, using the C++11 standard. Simulations rely on the [CNetwork library](https://github.com/VictorSeven/CNetwork) to perform computations. Just copy the `CNetwork` folder next to your `.cpp` files.

To generate the plots, Anaconda Python 3.x is required.   

## TO-DO list

- [ ] Split the old QIF code to make files just for network generation. Networks will be stored in files and just read for the dynamics.
- [ ] Implement Kuramoto dynamics in the network. Check possible frustrated synchronization.
- [ ] Implement type-II excitability model. Differences?
- [ ] Implement type-I excitability model. Differences?
- [ ] Start implementing measures for dynamical richness. Chimera index (see references of Villegas, Moretti, Muñoz 2014)
- [ ] Next steps??

  

  

