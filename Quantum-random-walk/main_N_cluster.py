import numpy as np

from qutip import *
import shutil
from functions_definition import *
from anw_N_data import main_data


# Coupling constant C0
C0 = 1000
g=70 #[m^-1 W^(1/2)]: nonlinear constant
L=0.020 #[m]: length of active propagation region
Nmodes=10 # Number of modes. In this case, number of waveguides.
Ndim=5> #Define the number of dimensions of the Hilbert Space
z_points=200 #The amount of points from the propagation to analize
nruns = 1 #Times that the program is run

# The coupling profile coefficient list has the structure flist=[f1,f2,...,f_(N-1)]
#In this case, we decided to add disorder to the WGA by getting f_i from a random normal distribution
# with average = 1 and standard deviation of sigma (variable), like a percentage (~0.1 = 10%)
sigma = 0  # Standard deviation

#This is a list of tuples that describes the injection to the WGA
#Injection[i][0] represents the mode, going from  to N
#Injection[i][1] represents the absolute value of alpha
#Injection[i][2] is the phase

odd_phase = np.pi/4


#injection_list = [(1,0, np.pi),(2,np.sqrt(0.0050), 0),(3,0, np.pi),(4,np.sqrt(0.0050), 0),(5,0, np.pi)]
#injection_list = [(1,np.sqrt(0.005), 0),(2,np.sqrt(0.005), np.pi/4),(3,np.sqrt(0.005), 0),(4,np.sqrt(0.005), np.pi/4),(5,np.sqrt(0.005), 0),(6,np.sqrt(0.005), np.pi/4),(7,np.sqrt(0.005), 0),(8,np.sqrt(0.005),np.pi/4),(9,np.sqrt(0.005),0),(10,np.sqrt(0.005),np.pi/4),(11,np.sqrt(0.005),0)]

injection_list_odd = [(i, np.sqrt(0.005), odd_phase) for i in range(1,Nmodes+1, 2)]
injection_list_even = [(i, np.sqrt(0.005), 0) for i in range(2,Nmodes+1,2)]


injection_list = injection_list_even+injection_list_odd
injection_list


#Get current path
current_path = os.getcwd()

#Create folder where the calculation results are going to be stored
results_folder=f'results_C0_{int(C0)}p{int((round(C0,2)-int(C0))*100)}_sigma_{int(sigma)}p{int((round(sigma,2)-int(sigma))*100)}_N_{Nmodes}_odd_phase_{int(odd_phase)}_p_{int((round(odd_phase,2)-int(odd_phase))*100)}'
path_results_folder = os.path.join(current_path,results_folder)



try:
    os.mkdir(path_results_folder)

except FileExistsError:
   shutil.rmtree(path_results_folder)
   os.mkdir(path_results_folder)



   # The data is generated calculating the evolved state the amount of times that the variable nruns indicates

for run in range(1,nruns+1):

    # Define the coupling profile coefficient list
    # The abs avoids negative values, but the distribution stays almost the same if sigma is small
    flist=[abs(np.random.normal(1,sigma)) for i in range(0,Nmodes)]
    #Parabolic profile array

    #flist=[np.sqrt(j*(Nmodes-j)/2) for j in range(1,Nmodes+1)]

    #Generates the data calculating the evolved state and extracting its information in every step L/zpoints
    main_data(current_path,path_results_folder,C0,flist, g, L, Nmodes,Ndim,z_points, injection_list, run)

os.chdir(current_path)
