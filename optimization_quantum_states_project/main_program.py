"""

@author: jblanco

Author: Jorge A. Blanco-Vásquez
email: jorge.blancovasquez@ucr.ac.cr

All these scripts can be improved using libraries like pandas to manage the .csv files
or making the projections just once.
"""

#Import libraries

from cmath import exp
import csv
import numpy as np
from qutip import *
import importlib
import os
import subprocess

from functions_definition_N3_3_injections_same_intensity import *
from mod_anw_N3_prop_recursive_run_3_inj_same_int_cluster_b import *
from mod_fidelity import *
from wigner_dist_plots_N3_3_inj_same_int_alpha_5mW_Ndim_20 import *
from amp_ratios_N3_3_inj_same_int_alpha_5mW import *
from bigger_prob_amplitudes_N3 import *


def main(file_name):
    """Calculates fidelity and wigner distribution for given experimentals parameters
    from a csv file generated from mathematica
    
    PARAMETERS:
    The name of the .csv file that contains all the data from the mathematica computation

    RETURNS:
    NONE
    Calls other functions that fill the .csv file with the results (fidelity and
    relative amplitudes) of th calculations with QuTiP and saves a figure of the Wigner distribution
    
    DEPENDECIES:
    -functions_definition_N3.py:
        It contains the definitions of some extra functions used in this 
        script.
    -mod_anw_N3_prop_recursive_run_3_inj_same_int_cluster_b
        Generates the evolved state
    -mod_fidelity
        Calculates fidelity
    -wigner_dist_plots_N3_3_inj_same_int_alpha_5mW_Ndim_20
        Obtains the Wigner distributions
    -amp_ratios_N3_3_inj_same_int_alpha_5mW
        Calculates relative amplitudes
    """

#The first part is to obtain the experimental optimized parameters of the function from .csv file from Mathematica
    with open(f"./{file_name}", 'r') as file:
        csvreader = csv.reader(file, delimiter=',')
        for row in csvreader:
            if row[0]=="projection":
                #Gets the projection list, which descibes the projections
                proj_list = [row[i] for i in range(1,5)]
            if row[0]=="target_arg":
                #Gets the target-state complex arguments of the amplitudes
                target_args=[]
                for i in range(1,7):
                    #These conditions change the arguments that are written in terms of pi
                    if row[i]=='Pi':
                        target_args.append(np.pi)
                    elif row[i]=='Pi/2':
                        target_args.append(np.pi/2)
                    elif row[i]=='Pi/3':
                        target_args.append(np.pi/3)
                    elif row[i]=='Pi/4':
                        target_args.append(np.pi/4)
                    else:
                        target_args.append(float(row[i]))
            if row[0]=="target_abs":
                #Gets the absolute value of the target-state amplitudes 
                target_abs = [float(row[i]) for i in range(1,7)]
            if row[0] == "par_values":
                #Gets the values of experimental parameters of the measurement of the state
                #The phase angle of mode 1 with respect tho mode 2
                phi_1= float(row[1])
                #The phase angle of mode 1 with respect tho mode 2
                phi_3= float(row[2])
                #The eigenvalue of the homodyne measurement
                x=float(row[3])
                #The phase angle of the cuadrature projection
                theta=float(row[4])
                #The coupling constant
                c_0 = float(row[5])
        file.close()

    #Reconstructs the probability amplitudes of each of the eigenstates
    prob_amp = [target_abs[i]*np.exp(target_args[i]*(1j)) for i in range(0,6)]

    #Creates a list of the eigenstates of the target state
    target_st = [i for i in range(0,6) if (target_abs[i]>0)]

    #Generates the data and the calculated state
    main_data(phi_1,phi_3,x,c_0,target_st,proj_list)
    #Makes the projections and fill the QuTiP fidelity value in the .csv file
    main_fidelity(phi_1,phi_3,x,theta,c_0,target_st,prob_amp,proj_list,file_name)
    #Makes the projections and saves a figure of the Wigner distribution
    main_wigner(phi_1,phi_3,x,theta,c_0,target_st,proj_list) 
    #Makes the projections and fill the QuTiP relative amplitudes values in the .csv file
    amp_ratios(phi_1,phi_3,x,theta,c_0,target_st,proj_list,file_name)
    #Makes a csv file with the states with bigger amplitudes
    
    results_folder=f'results_phi1_{int(phi_1)}p{int((round(phi_1,2)-int(phi_1))*100)}\
_phi3_{int(phi_3)}p{int((round(phi_3,2)-int(phi_3))*100)}_x_{int(x)}p{int((round(x,2)-int(x))*100)}_C0_{int(c_0)}p{int((round(c_0,2)-int(c_0))*100)}\
_ts_{target_st[0]}_{target_st[1]}'
    
    prob_amplitudes(results_folder, 20, threshold=0.0001)

####MAIN PROGRAM#####
#Obtains all the .csv files from Mathematica that are in the same folder and updates them

#List of .csv files

#If Windows doesn't like it, the following should work

#csv_list=[res_file for res_file in os.listdir() if res_file.startswith("results_mathematica")]

csv_list = subprocess.getoutput(f"ls results_mathematica*").split("\n")
for case in csv_list:
    #Fills with QuTiP results
    main(case)