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
from anw_N import *


#Calculates the evolved state given the phase angles phi_1,phi_2 and the coupling constant C0
"""User-defined quantities"""
#The amount of points from the propagation
# Coupling constant C0
C0 = 250
g=70 #[m^-1 W^(1/2)]: nonlinear constant
L=0.020 #[m]: length of active propagation region
Nmodes=7 # Number of modes. In this case, number of waveguides.
Ndim=6 #Define the number of dimensions of the Hilbert Space
z_points=20 #The amount of points from the propagation to analize
nruns = 10 #Times that the program is run
flist = [1 for i in range(0,Nmodes)]
# The coupling profile coefficient list has the structure flist=[f1,f2,...,f_(N-1)]
#In this case, we decided to add disorder to the WGA by getting f_i from a random normal distribution
# with average = 1 and standard deviation of sigma (variable), like a percentage (~0.1 = 10%)
sigma = 0.1  # Standard deviation
run = 1
#This is a list of tuples that describes the injection to the WGA
#Injection[0] represents the mode, going from  to N
#Injection[1] represents the absolute value of alpha
#Injection[2] is the phase
injection_list = [(4,np.sqrt(0.0050), 0)]
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



main_data(current_path,path_results_folder,C0,flist, g, L, Nmodes,Ndim,z_points, injection_list, run)

#threshold = 0.0001
