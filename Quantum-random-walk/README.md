# Quantum ramdom walk
This folder is called this way because it was inspired in the investigation "Quantum and Classical Correlations in Waveguide Lattices" and their representation of the quantum walk

# What do these files do?

The idea of constructing these programs is to simulate the evolution of the quantum state along the WGA and provide visual representations of the output state, a correlation matrix.
It can show the spread of light along the WGA, but this is only of interest when not injecting in all waveguides.

## How to use them?

To use this simulation, you have to run the main_N.ipynb notebook. The parameters are explained as comments in the notebook, however, it is important to consider that the parameter sigma was used to add disorder to the coupling profile of the WGA.
The "disorder" is aded by making the coupling coefficients random variables normally distributed with mean 1 and standard deviation sigma. If you don't want to introduce any random variable to the simulation just use sigma = 0 
