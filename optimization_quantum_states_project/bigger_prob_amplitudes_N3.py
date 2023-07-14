"""

@author: jblanco

Author: Jorge A. Blanco-Vásquez
email: jorge.blancovasquez@ucr.ac.cr

"""

#Import libraries

import numpy as np
from qutip import *
import pandas as pd
import os

def prob_amplitudes(path_results_folder, Ndim, threshold=0.0001):
    """
    This function saves the eigenstates and its amplitudes that conform the overall output state
    if the probability amplitude is bigger than some threshold value

    PARAMETERS
    path_results folder: The path to the .qu file
    Ndim: Define the number of dimensions of the Hilbert Space
    threshold: minimum value that the amplitude can take to save the state

    RETURN: NONE
    Saves a .csv file with the states which have a probability amplitude bigger than the treshold
    
    """
    # Loads the evolved state from the .qu file
    os.chdir(path_results_folder)
    evolved_state=qload('file_evolved_state_qutip_array')
    
    #Creates an empty dataframe
    df = pd.DataFrame({'state':[], 'amplitude': [], 'phase': []})

    #Makes all posible projections to obtain the their amplitudes and phases
    for i in range(0,6):
        for j in range(0,6):
            for k in range(0,6):
                        
                #Makes the projection
                comparing_state = tensor(fock(Ndim,i),fock(Ndim,j),fock(Ndim,k))
                projection = (comparing_state.dag()*evolved_state)[0][0][0][0][0][0] #No idea why I need to put all these zeros
                print(projection)
                #Transforms the qobj into an array and obtains its value
                proj_scalar = projection.full()[0][0]

                #Condition to save the state
                if np.abs(proj_scalar)>threshold:
                    df_i = pd.DataFrame({'state':[f'\'{i}{j}{k}\''],
                    'amplitude':[np.abs(proj_scalar)],
                    'phase':[np.angle(proj_scalar)]})

                    #Concatenates the dataframe with the old, this is like appending
                    df = pd.concat([df, df_i],ignore_index=True)

    # Saves the dataframe to .csv
    df.to_csv(f'bigger_amplitudes.csv', index=False)