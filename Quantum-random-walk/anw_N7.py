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
from functions_definition import *
import os
import pandas as pd

# To update the function definitions in the "functions_definition_N3.py" file
# if it were changed


def main_data(current_path,path_results_folder,phi_1, phi_2,C0,flist, g, L, Ndim, z_points,nrun):
    """
    This function calculates the evolution of the state starting from vacuum and passing
    through a 7-mode WGA which has a linear coupling constant of C0 for interactions involving mode 3
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
    save_qu_result_files=False


    Nmodes=7 # Number of modes. In this case, number of waveguides.
    Nmodes_list=[Nmodes];
    #list_to_file(path_results_folder,Nmodes_list,'nmodes_list.txt')
    #save Ndim value into a file
    Ndim_list=[Ndim];
    #list_to_file(path_results_folder,Ndim_list,'ndim_list.txt')

    # save user-defined quantities in a file
    file_calculation_parameters(path_results_folder,Ndim,C0,g,L,flist)

    ##############################################################################
    # The user-defined lists below should be introduced as tuples converted into
    # numpy arrays

    """
    Quantities of Mode 1 (waveguide 1)
    alpha1_abs and eta1_phase defined by the user
    """

    alpha1_abs = 0; #[W^(1/2)] list of absolute values of alpha_p that characterizes the strong *M*
    #coherent state injected in waveguide 1.
    #save list into a file
    #list_to_file(path_results_folder,alpha1_abs,'alpha1_abs.txt')

    eta1_phase=0; #For a coherent state, list of phases of alpha for mode 1
    #save list into a file
    #list_to_file(path_results_folder,eta1_phase,'eta1_phase.txt')

    """
    Quantities of Mode 2 (waveguide 2)
    alpha2_abs and eta2_phase defined by the user
    """
    alpha2_abs = 0 #[W^(1/2)] list of absolute values of alpha_p that characterizes the strong 
    #coherent state injected in waveguide 2.
    #save list into a file
    #list_to_file(path_results_folder,alpha2_abs,'alpha2_abs.txt')

    eta2_phase=0; #For a coherent state, list of phases of alpha for mode 2
    #save list into a file
    #list_to_file(path_results_folder,eta2_phase,'eta2_phase.txt')

    """
    Quantities of Mode 3 (waveguide 3)
    alpha2_abs and eta2_phase defined by the user
    """
    alpha3_abs = 0 #[W^(1/2)] list of absolute values of alpha_p that characterizes the strong 
    #coherent state injected in waveguide 3.
    #save list into a file
    #list_to_file(path_results_folder,alpha3_abs,'alpha3_abs.txt')

    eta3_phase=0; #For a coherent state, list of phases of alpha for mode 3
    #save list into a file
    #list_to_file(path_results_folder,eta3_phase,'eta3_phase.txt')

    """
    Quantities of Mode 4 (waveguide 4)
    alpha4_abs and eta4_phase defined by the user
    """
    
    alpha4_abs = np.sqrt(0.0050); #[W^(1/2)] list of absolute values of alpha_p that characterizes the strong *M*
    #coherent state injected in waveguide 4.
    #save list into a file
    #list_to_file(path_results_folder,alpha4_abs,'alpha4_abs.txt')

    eta4_phase=0; #For a coherent state, list of phases of alpha for mode 4
    #save list into a file
    #list_to_file(path_results_folder,eta4_phase,'eta4_phase.txt')

    """
    Quantities of Mode 5 (waveguide 5)
    alpha5_abs and eta5_phase defined by the user
    """

    alpha5_abs = 0; #[W^(1/2)] list of absolute values of alpha_p that characterizes the strong *M*
    #coherent state injected in waveguide 5.
    #save list into a file
    #list_to_file(path_results_folder,alpha5_abs,'alpha1_abs.txt')

    eta5_phase=0; #For a coherent state, list of phases of alpha for mode 5
    #save list into a file
    #list_to_file(path_results_folder,eta5_phase,'eta5_phase.txt')

    """
    Quantities of Mode 6 (waveguide 6)
    alpha6_abs and eta6_phase defined by the user
    """

    alpha6_abs = 0; #[W^(1/2)] list of absolute values of alpha_p that characterizes the strong *M*
    #coherent state injected in waveguide 6.
    #save list into a file
    #list_to_file(path_results_folder,alpha6_abs,'alpha1_abs.txt')

    eta6_phase=0; #For a coherent state, list of phases of alpha for mode 6
    #save list into a file
    #list_to_file(path_results_folder,eta6_phase,'eta6_phase.txt')

    """
    Quantities of Mode 7 (waveguide 7)
    alpha7_abs and eta7_phase defined by the user
    """

    alpha7_abs = 0; #[W^(1/2)] list of absolute values of alpha_p that characterizes the strong *M*
    #coherent state injected in waveguide 7.
    #save list into a file
    #list_to_file(path_results_folder,alpha7_abs,'alpha1_abs.txt')

    eta7_phase=0; #For a coherent state, list of phases of alpha for mode 7
    #save list into a file
    #list_to_file(path_results_folder,eta7_phase,'eta7_phase.txt')

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

    #define annihilation operator of mode 4
    a4=destruction_operator_multimode(Ndim, Nmodes, 4)
    #define creation operator of mode 4
    a4dag=creation_operator_multimode(Ndim, Nmodes, 4)

    #define annihilation operator of mode 5
    a5=destruction_operator_multimode(Ndim, Nmodes, 5)
    #define creation operator of mode 5
    a5dag=creation_operator_multimode(Ndim, Nmodes, 5)

    #define annihilation operator of mode 6
    a6=destruction_operator_multimode(Ndim, Nmodes, 6)
    #define creation operator of mode 6
    a6dag=creation_operator_multimode(Ndim, Nmodes, 6)

    #define annihilation operator of mode 7
    a7=destruction_operator_multimode(Ndim, Nmodes, 7)
    #define creation operator of mode 7
    a7dag=creation_operator_multimode(Ndim, Nmodes, 7)

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


    eta1_abs=g*alpha1_abs; #absolute value of eta1
    eta1=eta1_abs*np.exp(eta1_phase*1j); #phase value of eta1

    eta2_abs=g*alpha2_abs; #absolute value of eta2
    eta2=eta2_abs*np.exp(eta2_phase*1j); #phase value of eta2

    eta3_abs=g*alpha3_abs; #absolute value of eta3
    eta3=eta3_abs*np.exp(eta3_phase*1j); #phase value of eta3

    eta4_abs=g*alpha4_abs; #absolute value of eta4
    eta4=eta4_abs*np.exp(eta4_phase*1j); #phase value of eta4

    eta5_abs=g*alpha5_abs; #absolute value of eta5
    eta5=eta5_abs*np.exp(eta5_phase*1j); #phase value of eta5

    eta6_abs=g*alpha6_abs; #absolute value of eta6
    eta6=eta6_abs*np.exp(eta6_phase*1j); #phase value of eta6
 
    eta7_abs=g*alpha7_abs; #absolute value of eta7
    eta7=eta7_abs*np.exp(eta7_phase*1j); #phase value of eta7

 
    """
    Define the momentum operator for the case N=7
    """

    Gop_term1_N7=flist[1-1]*C0*a2*a1dag+flist[2-1]*C0*a3*a2dag+flist[3-1]*C0*a4*a3dag+flist[4-1]*C0*a5*a4dag+flist[5-1]*C0*a6*a5dag+flist[6-1]*C0*a7*a6dag\
+eta1*(a1dag)**2+eta2*(a2dag)**2+eta3*(a3dag)**2+eta4*(a4dag)**2+eta5*(a5dag)**2+eta6*(a6dag)**2+eta7*(a7dag)**2

    Gop_N7=Gop_term1_N7+Gop_term1_N7.dag()


    """
    Calculate evolved state
    """
    #list of distances z at which the evolved state will ve calculates
    zlist=np.array([zp*L/z_points for zp in range(0,z_points+1)]) # [m]
        
    evol_state = mcsolve(-Gop_N7, Psi_ini, zlist,options=options)
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
        name_solution_case='anw_n_7'
        qsave(evol_state,name_solution_case)
    
    #Container to save the data from the run

    data_dict = {
        'z':zlist,
        'av_pho_1':[],
        'av_pho_2':[],
        'av_pho_3':[],
        'av_pho_4':[],
        'av_pho_5':[],
        'av_pho_6':[],
        'av_pho_7':[],
        'prob_1':[],
        'prob_2':[],
        'prob_3':[],
        'prob_4':[],
        'prob_5':[],
        'prob_6':[],
        'prob_7':[],
    }
    
    for i in range(0,len(zlist)):
        
        """
        Save the evolved states in the evolve_states array
        """ 
        
        evolved_states=psi_evol[i];

        # Average output photon number of mode 1
        av_photon_number_list_mode1=expect(a1dag*a1,evolved_states)
        # Average output photon number of mode 2
        av_photon_number_list_mode2=expect(a2dag*a2,evolved_states)
        # Average output photon number of mode 3
        av_photon_number_list_mode3=expect(a3dag*a3,evolved_states)
        # Average output photon number of mode 4
        av_photon_number_list_mode4=expect(a4dag*a4,evolved_states)
        # Average output photon number of mode 5
        av_photon_number_list_mode5=expect(a5dag*a5,evolved_states)
        # Average output photon number of mode 6
        av_photon_number_list_mode6=expect(a6dag*a6,evolved_states)
        # Average output photon number of mode 7
        av_photon_number_list_mode7=expect(a7dag*a7,evolved_states)

        av_total_photon_number = av_photon_number_list_mode1+av_photon_number_list_mode2+av_photon_number_list_mode3+\
            av_photon_number_list_mode4+av_photon_number_list_mode5+av_photon_number_list_mode6+av_photon_number_list_mode7
        
        print(f'z = {zlist[i]}')

        print('average photon number in output state',av_photon_number_list_mode1,av_photon_number_list_mode2,\
            av_photon_number_list_mode3,av_photon_number_list_mode4,av_photon_number_list_mode5,av_photon_number_list_mode6,av_photon_number_list_mode7)

        """
        Save the "evolved_states" list in a .qu file
        """
        #qsave(evolved_states,f'file_evolved_state_qutip_array_{zlist[i]}')

        data_dict['av_pho_1'].append(av_photon_number_list_mode1)
        data_dict['av_pho_2'].append(av_photon_number_list_mode2)
        data_dict['av_pho_3'].append(av_photon_number_list_mode3)
        data_dict['av_pho_4'].append(av_photon_number_list_mode4)
        data_dict['av_pho_5'].append(av_photon_number_list_mode5)
        data_dict['av_pho_6'].append(av_photon_number_list_mode6)
        data_dict['av_pho_7'].append(av_photon_number_list_mode7)


        if av_total_photon_number !=0:
            prob_mode1 = av_photon_number_list_mode1/av_total_photon_number
            prob_mode2 = av_photon_number_list_mode2/av_total_photon_number
            prob_mode3 = av_photon_number_list_mode3/av_total_photon_number
            prob_mode4 = av_photon_number_list_mode4/av_total_photon_number
            prob_mode5 = av_photon_number_list_mode5/av_total_photon_number
            prob_mode6 = av_photon_number_list_mode6/av_total_photon_number
            prob_mode7 = av_photon_number_list_mode7/av_total_photon_number

            print('reason of photons',prob_mode1,prob_mode2,prob_mode3,prob_mode4,prob_mode5,prob_mode6,prob_mode7)

            #Calculate dispersion:
            dispersion = prob_mode1*(1**2)+prob_mode2*(2**2)+prob_mode3*(3**2)+prob_mode4*(4**2)+prob_mode5*(5**2)+prob_mode6*(6**2)+prob_mode7*(7**2)

            data_dict['prob_1'].append(prob_mode1)
            data_dict['prob_2'].append(prob_mode2)
            data_dict['prob_3'].append(prob_mode3)
            data_dict['prob_4'].append(prob_mode4)
            data_dict['prob_5'].append(prob_mode5)
            data_dict['prob_6'].append(prob_mode6)
            data_dict['prob_7'].append(prob_mode7)
        else:
            data_dict['prob_1'].append(0)
            data_dict['prob_2'].append(0)
            data_dict['prob_3'].append(0)
            data_dict['prob_4'].append(0)
            data_dict['prob_5'].append(0)
            data_dict['prob_6'].append(0)
            data_dict['prob_7'].append(0)


    #Save the dictionary with the data to a .csv file
    df = pd.DataFrame(data_dict)
    df.set_index('z',inplace=True)
    df.to_csv(f'results_nrun_{nrun}.csv')
    

    # return to the original main folder
    os.chdir(current_path)

    