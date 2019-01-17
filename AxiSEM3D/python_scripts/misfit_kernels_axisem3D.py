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
INPUT_DIR_WVF = home+'/Desktop/phd/private_Axisem3D/AxiSEM3D/build/output/wavefields/bubble_500km_0.06_2200km/'
INPUT_DIR_KER = home+'/Desktop/phd/private_Axisem3D/AxiSEM3D/build/output/kernels/kernel_nu100_70deg_mystf/'

list_files_ker = [f for f in sorted(os.listdir(INPUT_DIR_KER)) if 'kernels_db_' in f]
#list_files_wvf = [f for f in sorted(os.listdir(INPUT_DIR_WVF)) if '3d_wavefield_db' in f ]
list_files_wvf = [f for f in sorted(os.listdir(INPUT_DIR_WVF)) if 'wavefield_db_' in f and 'fwd' not in f and '3d' not in f]

if (len(list_files_ker) != len(list_files_wvf)):
    raise ValueError('Inconsistent number of wvf and kernel files')
num_files = len(list_files_ker)



chi_tot = 0

for i in range(num_files):
    
    ##### LOAD FILES #####
    dset_fwd = Dataset(INPUT_DIR_WVF+list_files_wvf[i])
    dset_ker = Dataset(INPUT_DIR_KER+list_files_ker[i])

    dm = dset_fwd.variables['vp'][:] - dset_fwd.variables['vp1D'][:]
    m = dset_fwd.variables['vp'][:]
    kernels = dset_ker.variables["Kernels"][:]
    nu = dset_fwd.variables['Nus'][:]
    nuKer = dset_ker.variables['Nus'][:]
    int_factor = dset_fwd.variables['integral_factor'][:]
    nelem = len(nu)
    

    ##### FFT, MULTIPLY, INTEGRATE ####

    nulineKer = 0
    nulineFwd = 0

    for ielem in range(nelem):
        
         
        nuelem = min(nu[ielem],nuKer[ielem])
        #test = np.zeros((5,5), dtype = complex)
        chi_elem = 0.
        
        for inu in range(nuelem):
            dm_c = dm[nulineFwd + inu, 0, :, :] + 1j * dm[nulineFwd + inu, 1, :, :]
            ker_c = kernels[nulineKer+inu, 0, :, :] + 1j * kernels[nulineKer + inu, 1, :, :]
            
            chi_elem +=np.real(np.sum( np.multiply(int_factor[ielem],np.multiply(dm_c, ker_c))  ))   
            #test += dm[nulineFwd + inu, 0, :, :] + 1j * dm[nulineFwd + inu, 1, : ,:]    
    #        print(np.max(dm[nulineFwd + inu, 0, :, :] + 1j * dm[nulineFwd + inu, 1, :, :]) )
            #test += kernels[nulineKer + inu, 4, :, :] + 1j * kernels[nulineKer + inu, 5, : ,:]    
        chi_tot +=(chi_elem)
        nulineKer += nuKer[ielem]
        nulineFwd += nu[ielem]

print(4 * np.pi * chi_tot)

