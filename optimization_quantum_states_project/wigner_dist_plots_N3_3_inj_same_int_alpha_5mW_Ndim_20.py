# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 10:37:08 2022

@author: earoj

Having a "evolved_states" array calculated from the 
"anw_N3_prop_recursive_run_3_injections_same_intensity.py", calculate the wigner distribution function of
the desired states.

--Jorge:
This file is a function to be called by the main program, the parameters are obtained from it
--

Author: EDGAR A. ROJAS-GONZALEZ
email: edgar.rojasgonzalez@ucr.ac.cr

"""


"""
Import libraries
"""

import matplotlib as mpl

from matplotlib import cm

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

import numpy as np

from IPython.display import Image

from qutip import *

from functions_definition_N3_3_injections_same_intensity import *

import os



# To update the function definitions in the "functions_definition_N3.py" file
# if it were changed

exec(open("./functions_definition_N3_3_injections_same_intensity.py").read())

def main_wigner(phi_1,phi_3,x_med,theta,C0,ts,projection_list):
    """
    PARAMETERS:

    This function recieves the angles of injection (phi_1,phi_3), the coupling constant (C_0);
    the target state (ts),
    information of the measurement (projection list) (explained in the main program)
    and the results of the measurement: eigenvalue asociated with the x cuadrature of the homodyne measurement (x_med) and its phase (theta).
    
    Some of these parameters are not actually used for calculations, but to generate the name of the files involved

    RETURNS:
    NONE
    This function saves the plots of the Wigner distributions of the calculated state obtained from the evolved state with the measurements performed to it  
    The figures are saved in formats (.eps, .pdf, .svg, jpg)
    
    DEPENDECIES:
    -functions_definition_N3.py:
        It contains the definitions of some extra functions used in this 
        script.
    """

    
    #LOCATE THE RESULTS FOLDER. The user has to introduce the path to the relevant folder
    #Get current path
    
    current_path = os.getcwd()
    
    
    results_folder=f'results_phi1_{int(phi_1)}p{int((round(phi_1,2)-int(phi_1))*100)}\
_phi3_{int(phi_3)}p{int((round(phi_3,2)-int(phi_3))*100)}_x_{int(x_med)}p{int((round(x_med,2)-int(x_med))*100)}_C0_{int(C0)}p{int((round(C0,2)-int(C0))*100)}\
_ts_{ts[0]}_{ts[1]}' #This is the name of the results folder
    
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
    # Navigate back to the initial current folder
    os.chdir(current_path)

    """
    #Calculate an array with the states on which the Wigner functios are calculated
    and plotted in the figure below.
    Case:
        -rho_1=Trace_3(rho_123 projected on the |x> quadrature state of mode 2)
    Wigner functions for mode 1 after tracing modes 2 and 3
    """

    # define array that contains the states from which the Wigner distribution 
    # functions will be calculated

    """
    indexing convention
    calc_state[j] corresponds to the state 
    used in the jth subfigure in the plots. 
    Here,
    a corresponds to 1,
    b to 2, 
    c to 3,
    etc 
    """

    calc_state = []

    for eta1_phase_index in [0]:   # VESTIGE
        for eta3_phase_index in [0]:   # VESTIGE
            calc_state.append([])


    ##############################################################################
    """
    The following instructions correspond to the projections made to the evolved state.
    The idea is to use the projection_list to make the indicated projections, in the proper order, as stated in the Mathematica notebook.
    I introduced the variable firstproj because the code (states and modes) changes after the first projection and other variables for similar reasons

    Each step is commented

    NOTES:
    -The for loops are unnecesary, the variables eta3_phase_index and 1 can be deleted since they're =0
    """

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
            except ZeroDivisionError:
                #If there's an error, the final state is empty, same as the plot
                psi_state_projected=[] 

            psi_state_1=psi_state_projected

            calc_state[counter_loop]=psi_state_1 #VESTIGE
            counter_loop +=1
    #######################################################

    """
    Making the Wigner distribution plot
    """
    fig1, axes = plt.subplots(figsize=(78,60))
    levels_variables=60

    states_wigner=calc_state

    psi_state_1=states_wigner[0]

    density_operator_1=psi_state_1*psi_state_1.dag()

    xvec = np.linspace(-3,3,201)

    W11 = wigner(density_operator_1, xvec, xvec)

    wmap = wigner_cmap(W11)  # Generate Wigner colormap

    plotwigner = axes.contourf(xvec, xvec, W11, levels_variables, cmap=wmap)

    plt.sca(axes)
    plt.xticks(fontsize=150)

    plt.yticks(fontsize=150)

    cb2a = fig1.colorbar(plotwigner, ax=axes)
    cb2a.ax.tick_params(labelsize=150)

    #TITLE
    axes.set_title(r"Wigner state for $\langle X_{\theta 3}|\langle 1_2|\Psi_{123}\rangle$""\n"f"with $\phi_1={round(phi_1,3)}$\
 ; $\phi_3={round(phi_3,3)}$ ; and $x={round(x_med,3)}$", fontsize = 150, pad = 100)

    axes.set_xlabel(r"$x_1$", fontsize=170)
    axes.set_ylabel(r"$p_1$", fontsize=170)

    axes.set_aspect('equal', 'box')

    #######################################################

    fig1.savefig(f'{results_folder}/wigner_dist_n3_3inj_same_int_alpha_5mw_rho_1_proj_2_n1_proj_x3_{int(x_med)}p{int((round(x_med,2)-int(x_med))*100)}_Ndim_100.eps')
    fig1.savefig(f'{results_folder}/wigner_dist_n3_3inj_same_int_alpha_5mw_rho_1_proj_2_n1_proj_x3_{int(x_med)}p{int((round(x_med,2)-int(x_med))*100)}_Ndim_100.jpg')
    fig1.savefig(f'{results_folder}/wigner_dist_n3_3inj_same_int_alpha_5mw_rho_1_proj_2_n1_proj_x3_{int(x_med)}p{int((round(x_med,2)-int(x_med))*100)}_Ndim_100.svg')
    fig1.savefig(f'{results_folder}/wigner_dist_n3_3inj_same_int_alpha_5mw_rho_1_proj_2_n1_proj_x3_{int(x_med)}p{int((round(x_med,2)-int(x_med))*100)}_Ndim_100.pdf')