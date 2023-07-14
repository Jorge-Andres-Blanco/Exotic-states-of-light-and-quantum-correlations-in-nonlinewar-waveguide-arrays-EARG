# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 13:08:54 2022

@author: earoj

Author: EDGAR A. ROJAS-GONZALEZ
email: edgar.rojasgonzalez@ucr.ac.cr
"""
import numpy as np
from qutip import *
import os

def list_to_file(path_results_folder,lname,fname):    
    """ function to save list to file
    fname: name of the file to be created (str)
    lname: list name (list)
    """
    #path_results_folder: path to the results folder
    textfile = open(os.path.join(path_results_folder,fname), "w")
    for element in lname:
        textfile.write(str(element) + "\n")
    textfile.close()
    
def file_calculation_parameters(path_results_folder,Ndim,C0,g,L,flist):
    """Create file with the summary of the calculation parameters""" 
    #path_results_folder: path to the results folder
    textfile = open(os.path.join(path_results_folder,'calculation_parameters_summary.txt'), "w")
    textfile.write('"""User-defined quantities"""'+'\n'+'\n')
    textfile.write('Ndim = '+str(Ndim)+' #Define the number of dimensions of the Hilbert Space'+'\n'+'\n')
    textfile.write('#quantities common to all waveguides'+'\n')
    textfile.write('C0= '+str(C0)+' #[m^-1]: linear coupling constant'+'\n')
    textfile.write('g= '+str(g)+' #[m^-1 W^(1/2)]: nonlinear constant'+'\n') 
    textfile.write('L= '+str(L)+' #[m]: length of active propagation region'+'\n')
    textfile.write('flist=[ \n')
    
    for iindex in range(len(flist)):
        textfile.write('f'+str(iindex+1)+'= '+ str(flist[iindex])+', \n')
    textfile.write('] \n') 
    textfile.write('#coupling profile coefficient '+'\n')
    textfile.close()
    
def txt_file_tonumerical_array(text_file_name):
    #text_file_name: String variable that contains the name of text file from 
    # which the numerical values are extracted. 
    # The .txt ending has to be added to this name!
    f=open(text_file_name)
    lines = f.readlines()
    f.close()
    array_return=[];
    for index in range(len(lines)):
        array_return.append(float(lines[index]))
    return array_return
    
def expectation_value_list(alpha1_abs,eta1_phase,alpha2_abs,eta2_phase,alpha3_abs,eta3_phase,operator,qstate_array):
    """Returns a list of expectation values of a given operator. This list
    presents the same indexing structure as the evolved_stated array
    Structure:        
        Variables:
            -alpha1_abs: 1D Tuple with the values of the absolute values of 
                alpha_p that characterizes the strong coherent state injected
                in waveguide 1.
            -eta1_phase: 1D Tuple with the values of the absolute values of the 
                phases of alpha (or eta) for mode 1.
            -alpha2_abs: 1D Tuple with the values of the absolute values of 
                alpha_p that characterizes the strong coherent state injected
                in waveguide 2.
            -eta2_phase: 1D Tuple with the values of the absolute values of the 
                phases of alpha (or eta) for mode 2.
            -alpha3_abs: 1D Tuple with the values of the absolute values of 
                alpha_p that characterizes the strong coherent state injected
                in waveguide 3.
            -eta3_phase: 1D Tuple with the values of the absolute values of the 
                phases of alpha (or eta) for mode 3.
            -operator: Quantum operator of which the expectation values will
                be calculated.
            -qstate_array: List of evolved quantum state used to calculate the
                expectation values.
    """ 
    # allocate the average photon number list
    # av_photon_number_list: list of the average photon numbers with the
    # same index structure as 
    # evolved_states[index1abs][index1phase][index2abs][index2phase][index3abs][index3phase]
    array_expectation_values=[];
    for index1abs in range (len(alpha1_abs)):
        array_expectation_values.append([])
        for index1phase in range (len(eta1_phase)):
            array_expectation_values[index1abs].append([])
            for index2abs in range (len(alpha2_abs)):
                array_expectation_values[index1abs][index1phase].append([])
                for index2phase in range (len(eta2_phase)):
                    array_expectation_values[index1abs][index1phase][index2abs].append([])
                    for index3abs in range (len(alpha3_abs)):
                        array_expectation_values[index1abs][index1phase][index2abs][index2phase].append([])
                        for index3phase in range (len(eta3_phase)):
                            array_expectation_values[index1abs][index1phase][index2abs][index2phase][index3abs].append([])
                            
    # Calculate the average photon numbers
    for index1abs in range (len(alpha1_abs)):
        for index1phase in range (len(eta1_phase)):
            for index2phase in range (len(eta2_phase)):
                for index3phase in range (len(eta3_phase)):
                    # calculate the expectation value 
                    array_expectation_values[index1abs][index1phase][index1abs][index2phase][index1abs][index3phase]\
                        =expect(operator,qstate_array[index1abs][index1phase][index1abs][index2phase][index1abs][index3phase])
    return array_expectation_values
    
    
def vacuum_state_multimode(Ndim, Nmodes):
    """ function that creates a vacuum state of a number Nmodes of modes. That
    is, the tensor product of the vacuum state of each mode. The variables are
        -Ndim: Hilbert state dimension for each mode.
        -Nmodes: Number of modes to be considered in the calculation. In this
            case, it corresponds to the number of waveguides.
    """
    qstate=fock(Ndim,0);
    for iindex in range(Nmodes-1):
        qstate=tensor(qstate,fock(Ndim,0))
    return qstate

def destruction_operator_multimode(Ndim, Nmodes, mode_number):
    """ function that creates an annihilation operator corresponding to the 
    mode "mode_number" of a multimode system of "Nmodes" number of modes, each 
    one described by a Hilbert space of "Ndim" dimensions. The variables are
        -Ndim: Hilbert state dimension for each mode.
        -Nmodes: Number of modes to be considered in the calculation. In this
            case, it corresponds to the number of waveguides.
        -mode_number: number that describes the mode corresponding to the
            annihilation operator.
    """
    for iindex in range(Nmodes):
        if (iindex==0):
            if (iindex+1==mode_number):
                operator=destroy(N=Ndim)
            else:
                operator=identity(Ndim)
        else: 
            if (iindex+1==mode_number):
                operator=tensor(operator,destroy(N=Ndim))
            else:
                operator=tensor(operator,identity(Ndim))
    return operator

def creation_operator_multimode(Ndim, Nmodes, mode_number):
    """ function that creates an creation operator corresponding to the 
    mode "mode_number" of a multimode system of "Nmodes" number of modes, each 
    one described by a Hilbert space of "Ndim" dimensions. The variables are
        -Ndim: Hilbert state dimension for each mode.
        -Nmodes: Number of modes to be considered in the calculation. In this
            case, it corresponds to the number of waveguides.
        -mode_number: number that describes the mode corresponding to the
            creation operator.
    """
    for iindex in range(Nmodes):
        if (iindex==0):
            if (iindex+1==mode_number):
                operator=create(N=Ndim)
            else:
                operator=identity(Ndim)
        else: 
            if (iindex+1==mode_number):
                operator=tensor(operator,create(N=Ndim))
            else:
                operator=tensor(operator,identity(Ndim))
    return operator
    