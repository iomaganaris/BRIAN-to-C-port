Implementation of a simulation of a biological neural network based on the ADEX model with pre- and post-synaptic expression of STDP in C.

The Python code of *prepostSTDP_savings.py* , the script responsible for running the desired simulation in the BRIAN neural network simulator, as well as the corresponding functions of BRIAN's source code were ported into a standalone simulator for the specific models in C++.

Doxygen documentation can be found in /doc/html/index.html

The main function of ported code is found in *src/simulation.cpp*

The original simulation is in https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=184487&file=%2FprepostSTDP_savings%2FprepostSTDP_savings.py#tabs-1
