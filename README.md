# HM Dynamics

Here we study the dynamics of hierarchical modular neuronal networks.  This is the code for  our paper [The broad edge of synchronisation: Griffiths effects and collective phenomena in brain networks](https://doi.org/10.48550/arXiv.2109.11783).

## Directory structure

- `cpp`: contains C++ code to make simulations. The code is separated into a file to create the network topology and another one to run the dynamics over a selected topology file. They both rely on [CNetwork](https://github.com/VictorSeven/CNetwork).
- `python`: contains Python code for launching programs and generating graphs. Code here allows to launch the C++ code massively on the cluster, as well as to invoke the network generation algorithm with desired parameters in a systematic way.

## Compiling sources

For C++ programs, it is enough to compile as always, using the C++11 standard. Simulations rely on the [CNetwork library](https://github.com/VictorSeven/CNetwork) to perform computations. Just copy the `CNetwork` folder next to these `.cpp` files.

To generate the plots, Anaconda Python 3.x is required.  





