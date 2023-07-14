# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 13:20:00 2022

@author: earoj


QuTiP implementation of the Array of Nonlinear Waveguides (ANW) propagation
problem.
Here, the scripts in the 
"...\Quantum_Optics\ANW\case_N3\qutip_recursive_run" folder are 
used as a basis to create a script (or set of scripts) that allows to perform
calculations assigning a set of values to the input parameters.

It is based on the scripts in

"...\GitHub\Quantum_Optics\ANW\case_N2"


Author: EDGAR A. ROJAS-GONZALEZ
email: edgar.rojasgonzalez@ucr.ac.cr


log:
    -creation date: 2022-01-20
    -creation date: 2022-01-11
    -2022-01-21: Save the array of evolved states directly in this script
        and give the oportunity to save or not the .qu QuTiP result files
    -2022-03-03: Adapt the script for 3 injections with identical intensity.
    
Dependencies:
    -functions_definition_N3.py:
        It contains the definitions of some extra functions used in this 
        script.

--Jorge:
This file is a function to be called by the main program, the parameters are obtained from it
"""

"""
Import libraries
"""

import numpy as np
import pandas as pd
from qutip import *
import importlib
from functions_definition import *
import os

# To update the function definitions in the "functions_definition_N3.py" file
# if it were changed


def main_data(current_path,path_results_folder,C0,flist, g, L, Nmodes, Ndim, z_points, injection_list, run):
    """
    This function calculates the evolution of the state starting from vacuum and passing
    through a N-mode WGA which has a linear coupling constant of C0 for interactions involving mode 3
    and a linear coupling constant of 4*C0 for the interactions involving first and last mode.
    In this case the phases of mode 1 and mode 2 are the same as those of mode 4 and 5, in order to have a symmetric system
    and all modes have the same injection.
    PARAMETERS:
    The angles of injection (phi_1,phi_2 [rad]) (which are the same as phi_4,phi_5)
    and the linear coupling constant (C_0) [m^-1];
    Coupling profile coefficient list (flist), it has the structure flist=[f1,f2,...,f_(N-1)]
    g [m^-1 W^(1/2)]: nonlinear constant
    L [m]: length of active propagation region
    Ndim: Define the number of dimensions of the Hilbert Space
    z_points, defines the number of points of the propagation for which the state is going to be saved
    Injection will be an array with first column representing the mode, second representing alpha, and third representing eta.    

    RETURNS:
    NONE
    Calculates the evolved state and creates a file that stores the information from it.
    Also prints the average photon number at the output of each mode and some probability
    amplitudes of some states of interest.

    DEPENDECIES:
    -functions_definition_N3.py:
        It contains the definitions of some extra functions used in this 
        script.
    """

    #############################################################################

    """
    User-defined quantities
    """

    # Option to save the .qu result file in the results folder.
    #   -save_qu_result_files=True: 
    #        The .qu QuTiP result files ARE saved in the results folder.
    #   -save_qu_result_files=False:
    #        The .qu QuTiP result files are NOT saved in the results folder.
    # It may be convenient not to save these files because they can consume
    # lots of storage space and their related names can get too long.
    save_qu_result_files=True

    # Number of modes. In this case, number of waveguides.;
    #list_to_file(path_results_folder,Nmodes_list,'nmodes_list.txt')
    #save Ndim value into a file
    Ndim_list=[Ndim];
    #list_to_file(path_results_folder,Ndim_list,'ndim_list.txt')

    # save user-defined quantities in a file
    file_calculation_parameters(path_results_folder,Ndim,C0,g,L,flist)

    ##############################################################################
    # The user-defined lists below should be introduced as tuples converted into
    # numpy arrays


    alpha_abs = np.zeros(Nmodes) #[W^(1/2)] list of absolute values of alpha
    eta_phase = np.zeros(Nmodes) #For a coherent state, list of phases
    for injection in injection_list:
        #Injection[0] represents the mode
        #Injection[1] represents the absolute value of alpha
        #Injection[2] is the phase
        alpha_abs[injection[0]-1] = injection[1]
        eta_phase[injection[0]-1] = injection[2]


    # End of user-defined section
    ##############################################################################


    """
    Define the initial quantum states
    """

    Psi_ini=vacuum_state_multimode(Ndim, Nmodes); #total initial state, vacuum state
    # the function 
    # -vacuum_state_multimode(Ndim, Nmodes)
    # is defined in the file "functions_definitions_N3.py"

    """
    Define the quantum operators
    """
    #the function 
    #-destruction_operator_multimode(Ndim, Nmodes, mode_number)
    #is defined in the file "functions_definitions_N3.py"

    #define annihilation operators for each mode
    a_list = [destruction_operator_multimode(Ndim, Nmodes, mode) for mode in range(1,Nmodes+1)]
    #define creation operators for each mode
    a_dag_list=[creation_operator_multimode(Ndim, Nmodes, mode) for mode in range(1,Nmodes+1)]

    """
    Modify options of dynamic solver if necessary
    see https://qutip.org/docs/latest/guide/dynamics/dynamics-options.html
    if required.
    """

    options = Options()
    options.nsteps = 4000

    """
    Perform the calculation of evolution of the quantum states
    """

    #The counters cabs1, cphase1, cabs2, cphase2, cabs3, and cphase3 used below 
    #are for rotulation purposes to identify the input conditions corresponding 
    #to each calculation

    # navigate to the results folder
    os.chdir(path_results_folder) 

    #obtaining the eta array
    eta_list = g*alpha_abs*np.exp(eta_phase*1j)

 
    """
    Define the momentum operator for the case with N waveguides
    """

    Gop_term1_N=0
    for mode in range(1,Nmodes+1):
        if mode !=Nmodes:
            Gop_term1_N+= flist[mode-1]*C0*a_list[mode]*a_dag_list[mode-1]
            Gop_term1_N+= eta_list[mode-1]*(a_dag_list[mode-1])**2
        else:
            Gop_term1_N+= eta_list[mode-1]*(a_dag_list[mode-1])**2

    Gop_N=Gop_term1_N+Gop_term1_N.dag()


    """
    Calculate evolved state
    """
    #list of distances z at which the evolved state will ve calculates
    zlist=[zp*L/z_points for zp in range(0,z_points+1)] # [m]
        
    evol_state = mcsolve(-Gop_N, Psi_ini, zlist,options=options)
    # the minus sign that multiplies Gop is explained in Logbook #1, page 15
    # Note that here we are using the momentum operator instead of the 
    # hamiltonian. That is
    # the hamiltonian H is exchanged by -Gop 
    # and the time t is exchanged by the propagation distance z

        
    if save_qu_result_files:
        # saves the .qu QuTiP files if save_qu_result_files=True is
        # selected.
        name_solution_case=f'anw_n_{Nmodes}'
        qsave(evol_state,name_solution_case)


    """
    Extract the evolved states
    """
    psi_evol=evol_state.states
    # save information from the solver
    # create name according to current initial input conditions



    data = np.zeros((Nmodes,len(zlist)+1))

    data[:,-1] = flist
    
    for i in range(0,len(zlist)):
        
        """
        Save the evolved states in the evolve_states array
        """ 
        
        evolved_states=psi_evol[i];

        # Average output photon number of all modes
        av_photon_number_list = np.array([expect(a_dag_list[mode]*a_list[mode],evolved_states) for mode in range(0,Nmodes)])
        data[:,i] = av_photon_number_list

    
    column_header_list = [str(round(z,5)) for z in zlist]
    column_header_list.append('f')
    df = pd.DataFrame(data, columns=column_header_list)
        


    #Save the dictionary with the data to a .csv file

    df.to_csv(f'results_nrun_{run}.csv', index=False)
    

    # return to the original main folder
    os.chdir(current_path)