"""

@author: jblanco

Author: Jorge A. Blanco-Vásquez
email: jorge.blancovasquez@ucr.ac.cr

"""

#Import libraries

from cmath import exp
import numpy as np
from qutip import *
import shutil
from functions_definition import *
from anw_N7 import *


#Calculates the evolved state given the phase angles phi_1,phi_2 and the coupling constant C0
"""User-defined quantities"""

phi_1 = 0
phi_2 = 0
C0 = 30
g=70 #[m^-1 W^(1/2)]: nonlinear constant
L=0.020 #[m]: length of active propagation region

Nmodes=7 # Number of modes. In this case, number of waveguides.
Nmodes_list=[Nmodes];
#list_to_file(path_results_folder,Nmodes_list,'nmodes_list.txt')
Ndim=6 #Define the number of dimensions of the Hilbert Space
#save Ndim value into a file
Ndim_list=[Ndim];
#list_to_file(path_results_folder,Ndim_list,'ndim_list.txt')

#quantities common to all waveguides

#coupling profile coefficient list
# it has the structure
# flist=[f1,f2,...,f_(N-1)]
#In this case, we decided to add disorder to the WGA by changing the f_i to random normal distribution
#The std deviation is sigma
sigma = 0.1
#The abs avoids negative values, but the distribution stays almost the same if sigma is small
flist=[abs(np.random.normal(1,sigma)) for i in range(0,Nmodes)] 


#LOCATE THE RESULTS FOLDER. The user has to introduce the path to the relevant folder
    
#Get current path
current_path = os.getcwd()

#Create folder where the calculation results are going to be stored
results_folder=f'results_C0_{int(C0)}p{int((round(C0,2)-int(C0))*100)}_sigma_{int(sigma)}p{int((round(sigma,2)-int(sigma))*100)}_N_{Nmodes}'

path_results_folder = os.path.join(current_path,results_folder)
try:
    os.mkdir(path_results_folder)

except FileExistsError:
   shutil.rmtree(path_results_folder)
   os.mkdir(path_results_folder)

# save user-defined quantities in a file
file_calculation_parameters(path_results_folder,Ndim,C0,g,L,flist)

#The amount of points from the propagation
z_points=10

nrun = 1

main_data(current_path,path_results_folder,phi_1, phi_2,C0,flist, g, L, Ndim,z_points, nrun)

#threshold = 0.0001
