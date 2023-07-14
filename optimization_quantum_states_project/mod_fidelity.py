# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 10:37:08 2022

@author: earoj

Having a "evolved_states" array calculated from the 
"anw_N3_prop_recursive_run_3_injections_same_intensity.py", calculate the wigner distribution function of
the desired states.


Author: EDGAR A. ROJAS-GONZALEZ
email: edgar.rojasgonzalez@ucr.ac.cr

--Jorge:
This file is a function to be called by the main program, the parameters are obtained from it
"""

"""
Import libraries
"""
#matplotlib inline
from distutils import filelist
import matplotlib as mpl
import csv
from matplotlib import cm

import matplotlib.pyplot as plt
from numpy.core.fromnumeric import transpose

from scipy.optimize import curve_fit

import numpy as np
from IPython.display import Image
from qutip import *
from functions_definition_N3_3_injections_same_intensity import *
import os


# To update the function definitions in the "functions_definition_N3.py" file
# if it were changed


def main_fidelity(phi_1,phi_3,x_med,theta,c_0,ts,amp_ts,projection_list,file_name):
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
    Saves the value of fidelity in the csv file and also prints other fidelities (the eigenvalues of the state). The values 
    of the fidelity are obtained from comparing the calculated state with the target state, in this case, it compares with 
    the target state, the eigenstates that compose it. Compares states of the form:

    a*|n>+b*|m> ; with n,m integers. 

    With the calculated state with obtained from the corresponding projection.

    DEPENDECIES:
    -functions_definition_N3.py:
        It contains the definitions of some extra functions used in this 
        script.
    """
    
    #LOCATE THE RESULTS FOLDER. The user has to introduce the path to the relevant folder
    #Get current path

    current_path = os.getcwd()
    
    results_folder=f'results_phi1_{int(phi_1)}p{int((round(phi_1,2)-int(phi_1))*100)}\
_phi3_{int(phi_3)}p{int((round(phi_3,2)-int(phi_3))*100)}_x_{int(x_med)}p{int((round(x_med,2)-int(x_med))*100)}_C0_{int(c_0)}p{int((round(c_0,2)-int(c_0))*100)}\
_ts_{ts[0]}_{ts[1]}'#This is the name of the results folder
    
    #Navigate to the folder where the calculation results are going to be stored
    path_results_folder = os.path.join(current_path,results_folder)

    os.chdir(path_results_folder)

    ##############################################################################

    """
    upload alphaj_abs_list and etaj_phase_list lists
    """

    """
    Define the dimensionality of the Hilbert spaces
    It should be the same as the one used in the solution to the evolution problem
    """
    # Define the number of modes. That is, the number of waveguides of the ANWs.
    Nmodes_list=txt_file_tonumerical_array('nmodes_list.txt')

    Nmodes=int(Nmodes_list[0])


    # Define the number of dimensions of the Hilbert Space
    # Read the value of Ndim from its respective file in the results folder
    Ndim_list=txt_file_tonumerical_array('ndim_list.txt')

    Ndim=int(Ndim_list[0])

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


    ####################
    #Mode 1

    #alpha1_abs_list
    alpha1_abs=np.array(tuple(txt_file_tonumerical_array('alpha1_abs_list.txt')))

    #eta1_phase_list
    eta1_phase=np.array(tuple(txt_file_tonumerical_array('eta1_phase_list.txt')))

    ####################
    #Mode 2

    #alpha2_abs_list
    alpha2_abs=np.array(tuple(txt_file_tonumerical_array('alpha2_abs_list.txt')))

    #eta2_phase_list
    eta2_phase=np.array(tuple(txt_file_tonumerical_array('eta2_phase_list.txt')))

    ####################
    #Mode 3

    #alpha2_abs_list
    alpha3_abs=np.array(tuple(txt_file_tonumerical_array('alpha3_abs_list.txt')))

    #eta2_phase_list
    eta3_phase=np.array(tuple(txt_file_tonumerical_array('eta3_phase_list.txt')))

    """
    The evolved_states array contains the evolved states and
    presents the following structure

        -evolved_states[index1abs][index1phase][index2abs][index2phase]
            [index3abs][index3phase]

    with:
        -index1abs: the index of the related element of the "alpha1_abs" list.
        -index1phase: the index of the related element of the "eta1_phase" list.
        -index2abs: the index of the related element of the "alpha2_abs" list.
        -index2phase: the index of the related element of the "eta2_phase" list.
        -index3abs: the index of the related element of the "alpha3_abs" list.
        -index3phase: the index of the related element of the "eta3_phase" list.

    """
                    

    """
    Load the "evolved_states" list in a .qu file
    """

    evolved_states=qload('file_evolved_state_qutip_array')

    """
    Return to the initial folder
    """
    # navigate back to the initial current folder
    os.chdir(current_path)


    #########################################################

    # Define array that contains the fidelities between the corresponding quantum state

    """ Case:
        -rho_2=<x3|(<1_1|rho_123|1_1>)|x3>, x3=0.1 compared to
        -rho_2=|1><1|
    """
    calc_state_vs_eigen_1= []

    """ Case:
        -rho_2=<x3|(<1_1|rho_123|1_1>)|x3>, x3=0.1 compared to
        -rho_2=|0><0|
    """
    calc_state_vs_eigen_0= []

    """ Case:
        -rho_2=<x3|(<1_1|rho_123|1_1>)|x3>, x3=0.1 compared to
        -rho_2=(1/2)((|1>+|0>)(<0|+<1|) 
    """
    calc_state_vs_eigen_target= []



    for eta1_phase_index in [0]: #VESTIGE
        calc_state_vs_eigen_1.append([])
        calc_state_vs_eigen_0.append([])
        calc_state_vs_eigen_target.append([])


    ##############################################################################
    counter_loop=0 #VESTIGE
    for eta1_phase_index in [0]:   # VESTIGE
        for eta3_phase_index in [0]:   # VESTIGE
            index_alpha=0; #I made index_alpha=0 instead of 1 because the last was out of range 
            psi_state=evolved_states[index_alpha][eta1_phase_index][index_alpha][0][index_alpha][eta3_phase_index]

            # FIRST, define the stated which we'll project to
            # define the eigenvalue of the in-phase quadrature X used here
            xvalue=x_med
            thetavalue=theta
            # define the operator in-phase quadrature in the Fock space
            xstate=np.pi**(1/4)\
                *((1j)*thetavalue*create(N=Ndim)*destroy(N=Ndim)).expm()*(-(xvalue**2/2)+np.sqrt(2)*xvalue*create(N=Ndim)-(create(N=Ndim))**2/2).expm()*vacuum_state_multimode(Ndim, 1)
            try:
            #THE FOLLOWING CORRRESPONDS TO THE PROJECTIONS

                if bool(projection_list[3]) == False: #If 4th element in list is FALSE, the order of the PROJECTIONS has to be REVERSED
                    
                    projection_list.pop(3) #Eliminates the False element so it doesn't end  at the begining of the list
                    projection_list.reverse() #reverse
                
                firstproj=True #This variable indicates that the FOLLOWING projection IS the FIRST one

                for proj in projection_list: #A loop to access elemets of the projection list and make the projections

                    mode = projection_list.index(proj) 
                    
                    if proj == "x": #IF the PROJECTION corresponds to X (cuadrature)
                        
                        if firstproj: #IF this is the FIRST PROJECTION, there are three modes available
                            
                            tens_list=[identity(Ndim),identity(Ndim),identity(Ndim)] #Initialize the operator for the projection
                            tens_list[mode]=xstate #Puts the ket X IN the corresponding MODE
                            psi_state_projected_1 = (tensor(tens_list).dag()*psi_state).unit() #Makes the projection
                            mode_first_proj = mode #Saves the mode in which the first projection was made, this will be used to eliminate this mode 
                            firstproj=False #The following projection is NOT going to be the first one
                        
                        else: #If this corresponds to the second projection
                        
                            tens_list[mode]=xstate #Puts the ket x in the corresponding mode
                            tens_list.pop(mode_first_proj) #Deletes the element of the list of corresponding to the first projection
                            psi_state_projected = (tensor(tens_list).dag()*psi_state_projected_1).unit() #Makes the projection
                    try:

                        photons=int(proj) #IF the element of the list IS a INTEGER, it corresponds to a PROJECTION to a number of PHOTONS (executes the code inside the try)
                        
                        if firstproj: #IF this is the FIRST PROJECTION, there are three modes available
                        
                            tens_list=[identity(Ndim),identity(Ndim),identity(Ndim)] #Initialize the operator for the projection
                            tens_list[mode]=fock(Ndim,photons) #puts the ket of photons IN the corresponding MODE
                            psi_state_projected_1 = (tensor(tens_list).dag()*psi_state).unit() #Makes the projection
                            mode_first_proj = mode #Saves the mode in which the first projection was made, this will be used to eliminate this mode
                            firstproj=False #The following projection is NOT going to be the first one
                        
                        else: #If this corresponds to the second projection
                        
                            tens_list[mode]=fock(Ndim,photons) #The same idea as before
                            tens_list.pop(mode_first_proj)
                            psi_state_projected = (tensor(tens_list).dag()*psi_state_projected_1).unit()
                    except:
                        pass
                #Evaluates fidelity comparing with target state and its composing eigenstates
                #vs |ts[1]> ---> Second eigenstate of target state
                calc_state_vs_eigen_1[counter_loop].append(metrics.fidelity(psi_state_projected,fock(Ndim,ts[1])))
                #vs |ts[0]> ---> First eigenstate of target state
                calc_state_vs_eigen_0[counter_loop].append(metrics.fidelity(psi_state_projected,fock(Ndim,ts[0])))
                #vs (|0>+|1>)(1/sqrt(2)) ---> Target state *()
                calc_state_vs_eigen_target[counter_loop].append(metrics.fidelity(psi_state_projected,(amp_ts[1]\
                    *fock(Ndim,1)+amp_ts[0]*fock(Ndim,0)+amp_ts[2]*fock(Ndim,2)+amp_ts[3]*fock(Ndim,3)).unit()))
            
            except ZeroDivisionError: #Case of error
                print("Division by Zero")
                calc_state_vs_eigen_1[counter_loop].append(np.nan)
                calc_state_vs_eigen_0[counter_loop].append(np.nan)
                calc_state_vs_eigen_target[counter_loop].append(np.nan)

        counter_loop +=1
    
    file_list=[]

    #This part includes the fidelity in the .csv file

    #Open the .csv and saves it in a variable
    with open(f"./{file_name}", 'r') as file:
            csvreader = csv.reader(file, delimiter=',')
            for row in csvreader:
                file_list.append(row)
            file.close()
    
    #Includes fidelity value for the target state
    for element in file_list:
        if element[0]=="fidelity_py":
            element[1] = str(calc_state_vs_eigen_target[0][0])

    #Rewrites all the .csv
    with open(f"./{file_name}", 'w') as file:
            csvwriter = csv.writer(file)
            csvwriter.writerows(file_list)
            file.close()

    #Prints all fidelity values obtained
    print(f"the value of the fidelity when comparig to the state \
        |0> is {calc_state_vs_eigen_0}")
    print(f"the value of the fidelity when comparig to the state \
        |{ts[1]}> is {calc_state_vs_eigen_1}")
    print(f"the value of the fidelity when comparig to the state \
        (1/sqrt{amp_ts[0]**2+amp_ts[1]**2}({amp_ts[0]}|{ts[0]}>+{amp_ts[1]}|{ts[1]}>)) is {calc_state_vs_eigen_target}")