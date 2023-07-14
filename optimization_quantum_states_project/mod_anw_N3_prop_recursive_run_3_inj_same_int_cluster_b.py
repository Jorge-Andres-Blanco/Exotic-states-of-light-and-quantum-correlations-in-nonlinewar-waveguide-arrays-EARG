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
from qutip import *
import importlib
from functions_definition_N3_3_injections_same_intensity import *
exec(open("./functions_definition_N3_3_injections_same_intensity.py").read())
import os
import shutil

# To update the function definitions in the "functions_definition_N3.py" file
# if it were changed


def main_data(phi_1, phi_3, x_med,C0_med, ts,projection_list):
    """
    PARAMETERS:

    This function recieves the angles of injection (phi_1,phi_3), the coupling constant (C_0);
    the target state (ts),
    information of the measurement (projection list) (explained in the main program),
    the results of the measurement: eigenvalue asociated with the x cuadrature of the homodyne measurement (x_med) and its phase (theta),
    and filename: the name of the related .csv file
    
    Some of these parameters are not actually used for calculations, but to generate the name of the files involved

    RETURNS:
    NONE
    Calculates the evolved state and creates a file for storing information about that case.

    DEPENDECIES:
    -functions_definition_N3.py:
        It contains the definitions of some extra functions used in this 
        script.
    """
    
    #LOCATE THE RESULTS FOLDER. The user has to introduce the path to the relevant folder
    
    #Get current path
    current_path = os.getcwd()
    
    #Create folder where the calculation results are going to be stored
    results_folder=f'results_phi1_{int(phi_1)}p{int((round(phi_1,2)-int(phi_1))*100)}\
_phi3_{int(phi_3)}p{int((round(phi_3,2)-int(phi_3))*100)}_x_{int(x_med)}p{int((round(x_med,2)-int(x_med))*100)}_C0_{int(C0_med)}p{int((round(C0_med,2)-int(C0_med))*100)}\
_ts_{ts[0]}_{ts[1]}'
    
    path_results_folder = os.path.join(current_path,results_folder)
    try:
        os.mkdir(path_results_folder)

    except FileExistsError:
       shutil.rmtree(path_results_folder)
       os.mkdir(path_results_folder)


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
    save_qu_result_files=False


    Nmodes=3 # Number of modes. In this case, number of waveguides.
    Nmodes_list=[Nmodes];
    list_to_file(path_results_folder,Nmodes_list,'nmodes_list.txt')
    Ndim=20 #Define the number of dimensions of the Hilbert Space
    #save Ndim value into a file
    Ndim_list=[Ndim];
    list_to_file(path_results_folder,Ndim_list,'ndim_list.txt')

    #quantities common to all waveguides
    C0=C0_med #[m^-1]: linear coupling constant
    #C0=0 #[m^-1]: linear coupling constant

    g=10 #[m^-1 W^(1/2)]: nonlinear constant
    #g=0 #[m^-1 W^(1/2)]: nonlinear constant

    L=0.020 #[m]: length of active propagation region

    #coupling profile coefficient list
    # it has the structure
    # flist=[f1,f2,...,f_(N-1)]
    f1=1; f2=1;
    flist=[f1,f2]

    # save user-defined quantities in a file
    file_calculation_parameters(path_results_folder,Ndim,C0,g,L,flist)

    ##############################################################################
    # The user-defined lists below should be introduced as tuples converted into
    # numpy arrays

    """
    Quantities of Mode 1 (waveguide 1)
    alpha1_abs_list and eta1_phase_list defined by the user
    """
    #alpha1_abs_list = np.sqrt((0.0005,0.005,0.050)); #[W^(1/2)] list of absolute values of alpha_p that characterizes the strong 
    alpha1_abs_list = np.sqrt((0.0050,)); #[W^(1/2)] list of absolute values of alpha_p that characterizes the strong *M*
    #coherent state injected in waveguide 1.
    #save list into a file
    list_to_file(path_results_folder,alpha1_abs_list,'alpha1_abs_list.txt')

    eta1_phase_list=(phi_1,); #For a coherent state, list of phases of alpha for mode 1
    #save list into a file
    list_to_file(path_results_folder,eta1_phase_list,'eta1_phase_list.txt')

    """
    Quantities of Mode 2 (waveguide 2)
    alpha2_abs_list and eta2_phase_list defined by the user
    """
    alpha2_abs_list = alpha1_abs_list #[W^(1/2)] list of absolute values of alpha_p that characterizes the strong 
    #coherent state injected in waveguide 2.
    #save list into a file
    list_to_file(path_results_folder,alpha2_abs_list,'alpha2_abs_list.txt')

    eta2_phase_list=np.array((np.pi/2,)); #For a coherent state, list of phases of alpha for mode 2
    #save list into a file
    list_to_file(path_results_folder,eta2_phase_list,'eta2_phase_list.txt')

    """
    Quantities of Mode 3 (waveguide 3)
    alpha2_abs_list and eta2_phase_list defined by the user
    """
    alpha3_abs_list = alpha1_abs_list #[W^(1/2)] list of absolute values of alpha_p that characterizes the strong 
    #coherent state injected in waveguide 3.
    #save list into a file
    list_to_file(path_results_folder,alpha3_abs_list,'alpha3_abs_list.txt')

    eta3_phase_list=(phi_3,); #For a coherent state, list of phases of alpha for mode 3
    #save list into a file
    list_to_file(path_results_folder,eta3_phase_list,'eta3_phase_list.txt')


    # End of user-defined section
    ##############################################################################


    """
    Define the initial quantum states
    """

    # #initial state of mode 1
    # Psi1_ini=fock(Ndim,0) #vacuum state

    # #initial state of mode 2
    # Psi2_ini=fock(Ndim,0) #vacuum state

    # #initial state of mode 3
    # Psi2_ini=fock(Ndim,0) #vacuum state

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

    #define annihilation operator of mode 1
    a1=destruction_operator_multimode(Ndim, Nmodes, 1)
    #define creation operator of mode 1
    a1dag=creation_operator_multimode(Ndim, Nmodes, 1)

    #define annihilation operator of mode 2
    a2=destruction_operator_multimode(Ndim, Nmodes, 2)
    #define creation operator of mode 2
    a2dag=creation_operator_multimode(Ndim, Nmodes, 2)

    #define annihilation operator of mode 3
    a3=destruction_operator_multimode(Ndim, Nmodes, 3)
    #define creation operator of mode 3
    a3dag=creation_operator_multimode(Ndim, Nmodes, 3)


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

    """
    Initialize arrays whe evolved quantum states are going to be stored
    """
    evolved_states=[];
    for index1abs in range (len(alpha1_abs_list)):
        evolved_states.append([])
        for index1phase in range (len(eta1_phase_list)):
            evolved_states[index1abs].append([])
            for index2abs in range (len(alpha2_abs_list)):
                evolved_states[index1abs][index1phase].append([])
                for index2phase in range (len(eta2_phase_list)):
                    evolved_states[index1abs][index1phase][index2abs].append([])
                    for index3abs in range (len(alpha3_abs_list)):
                            evolved_states[index1abs][index1phase][index2abs][index2phase].append([])
                            for index3phase in range (len(eta3_phase_list)):
                                evolved_states[index1abs][index1phase][index2abs][index2phase][index3abs].append([])
    
    """
    According to this, the evolved_states array will present the following 
    structure

        -evolved_states[index1abs][index1phase][index2abs][index2phase][index3abs][index3phase]

    with:
        -index1abs: the index of the related element of the "alpha1_abs" list.
        -index1phase: the index of the related element of the "eta1_phase" list.
        -index2abs: the index of the related element of the "alpha2_abs" list.
        -index2phase: the index of the related element of the "eta2_phase" list.
        -index3abs: the index of the related element of the "alpha3_abs" list.
        -index3phase: the index of the related element of the "eta3_phase" list.

    """


    #The counters cabs1, cphase1, cabs2, cphase2, cabs3, and cphase3 used below 
    #are for rotulation purposes to identify the input conditions corresponding 
    #to each calculation

    # navigate to the results folder
    os.chdir(path_results_folder) 

    #initialize the counter of the number of calculations
    counter_calculation=0

    #initiallize counter for alpha1_abs_list
    ca1=0
    for alpha1_abs in alpha1_abs_list:
        #initiallize counter for eta1_phase_list
        cp1=0
        for eta1_phase in eta1_phase_list:
            #initiallize counter for eta2_phase_list
            cp2=0
            for eta2_phase in eta2_phase_list:
                #initiallize counter for alpha3_abs_list
                cp3=0
                for eta3_phase in eta3_phase_list:
                    eta1_abs=g*alpha1_abs; #absolute value of eta1
                    eta1=eta1_abs*np.exp(eta1_phase*1j); #phase value of eta1

                    eta2_abs=g*alpha1_abs; #absolute value of eta2
                    eta2=eta2_abs*np.exp(eta2_phase*1j); #phase value of eta2

                    eta3_abs=g*alpha1_abs; #absolute value of eta3
                    eta3=eta3_abs*np.exp(eta3_phase*1j); #phase value of eta3

                    """
                    Define the momentum operator for the case N=2
                    """
                
                    Gop_term1_N3=flist[1-1]*C0*a2*a1dag+flist[2-1]*C0*a3*a2dag\
                        +eta1*(a1dag)**2+eta2*(a2dag)**2+eta3*(a3dag)**2
                
                    Gop_N3=Gop_term1_N3+Gop_term1_N3.dag()


                    """
                    Calculate evolved state
                    """
                    #list of distances z at which the evolved state will ve calculates
                    zlist=np.array((0,L)) # [m]
                        
                    evol_state = mcsolve(-Gop_N3, Psi_ini, zlist,options=options)
                    # the minus sign that multiplies Gop is explained in Logbook #1, page 15
                    # Note that here we are using the momentum operator instead of the 
                    # hamiltonian. That is
                    # the hamiltonian H is exchanged by -Gop 
                    # and the time t is exchanged by the propagation distance z
                            
                    """
                    Extract the evolved states
                    """
                    psi_evol=evol_state.states
                    # save information from the solver
                    # create name according to current initial input conditions
                        
                    if save_qu_result_files:
                        # saves the .qu QuTiP files if save_qu_result_files=True is
                        # selected.
                        name_solution_case=\
                        'anw_n_3_ca1_'+str(ca1)+'_cp1_'+str(cp1)+\
                            '_ca2_'+str(ca2)+'_cp2_'+str(cp2)+\
                                '_ca3_'+str(ca3)+'_cp3_'+str(cp3)
                        qsave(evol_state,name_solution_case)
                    """
                    Save the evolved states in the evolve_states array
                    """ 
                    evolved_states[ca1][cp1][ca1][cp2][ca1][cp3]\
                        =psi_evol[1];

                    #add 1 to the counters
                    counter_calculation+=1
                    print(counter_calculation)
                    cp3+=1
                cp2+=1
            cp1+=1
        ca1+=1
        
    """
    Save the "evolved_states" list in a .qu file
    """

    qsave(evolved_states,'file_evolved_state_qutip_array')

    # return to the original main folder
    os.chdir(current_path)

#Parameters: phi_1,phi_3,x,targetstate (targ.state is a list with the prob. amplitude
# similar to the mathematica function. examp: [0,1])