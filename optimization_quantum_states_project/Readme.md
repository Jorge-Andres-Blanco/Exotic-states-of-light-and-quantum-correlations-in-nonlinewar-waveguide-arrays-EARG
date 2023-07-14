# This folder contains all programs used in the project "Experimental parameters optimization for non-classical quantum states of light in nonlinear WGA"

In the Mathematica notebook we constructed an algebraic representation od the quantum state that depended on some experimental parameters as well as the results of the measurements performed in the WGA.
In this notebook, one can choose the target state wanted at the output of one of the three waveguides, for instance:

$$|\Psi_t\rangle = a|0\rangle+b|1\rangle$$

As a result, the notebook will optimize 5 of these experimental parameters to construct the most similar to the target state it can provide.
It will then save the optimized parameters, the amplitudes of the target and generated states and their fidelities (metric to compare states) in a .csv file.

The .csv can be later used by the "main program", which takes automatically grabs the information of the optimized parameters from the .csv file and makes a more accurate calculation of the output state.
With this output state, one can compare the results from the Mathematica notebook and the python simulation to see if the first is accurate enogh and compare the fidelities given by the mathematica and python simulations.
It is expected that the python numerical simulation gives more accurate results, for this reason, it is considered as a reference.

Finally, it is important to mention that "main_program.py" must be runned in the same folder as the .csv file and all the other .py files (functions) are. Otherwise it will output a FileNotFoundError.
