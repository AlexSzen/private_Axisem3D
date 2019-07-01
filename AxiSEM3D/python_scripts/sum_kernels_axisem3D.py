'''

inputs : unperturbed model parameters, perturbed model params,
kernels obtained in perturbed model

purpose : compute misfit by integrating kernels*perturbation over domain.

'''





import numpy as np
from netCDF4 import Dataset
import os




##### INPUTS #####
home = os.getenv("HOME")
INPUT_DIR_WVF = home+'/Desktop/phd/private_Axisem3D/AxiSEM3D/build/output/wavefields/'
INPUT_DIR_KER = home+'/Desktop/phd/private_Axisem3D/AxiSEM3D/build/output/kernels/'

list_files_ker = [f for f in sorted(os.listdir(INPUT_DIR_KER)) if '30deg' in f]
list_files_wvf = [f for f in sorted(os.listdir(INPUT_DIR_WVF)) if 'wavefield_db_fwd' in f]

if (len(list_files_ker) != len(list_files_wvf)):
    raise ValueError('Inconsistent number of wvf and kernel files')
num_files = len(list_files_ker)



chi_tot = 0

for i in range(num_files):
    
    ##### LOAD FILES #####
    dset_fwd = Dataset(INPUT_DIR_WVF+list_files_wvf[i])
    dset_ker = Dataset(INPUT_DIR_KER+list_files_ker[i])

    kernels = dset_ker.variables["Kernels"][:]
    nuKer = dset_ker.variables['Nus'][:]
    int_factor = dset_fwd.variables['integral_factor'][:]
    nelem = len(nuKer)
    

    ##### FFT, MULTIPLY, INTEGRATE ####

    nulineKer = 0
    nulineFwd = 0

    for ielem in range(nelem):
        
        nuelem = nuKer[ielem] 
        chi_elem = 0.
        
        for inu in range(nuelem):
            ker_c = kernels[nulineKer+inu, 0, :, :] + 1j * kernels[nulineKer + inu, 1, :, :]
            
            chi_elem +=np.real(np.sum( np.multiply(int_factor[ielem],ker_c)  ))   
        
        chi_tot +=(chi_elem)
        nulineKer += nuKer[ielem]

print(4 * np.pi * chi_tot)

