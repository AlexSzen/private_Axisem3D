'''

inputs : unperturbed model parameters, perturbed model params,
kernels obtained in perturbed model

purpose : compute misfit by integrating kernels*perturbation over domain.

'''





import numpy as np
from netCDF4 import Dataset
import sys




##### INPUTS #####

INPUT_WVF = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/wavefields/wavefield_db.nc4'  
INPUT_WVF1D = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/wavefields/wavefield_db_fwd.nc4'  
INPUT_KER = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/kernels/kernels_db.nc4'


##### LOAD FILES #####

dset_fwd = Dataset(INPUT_WVF)
dset_ker = Dataset(INPUT_KER)
dset_fwd1D = Dataset(INPUT_WVF1D) #just to get domain decomposition and integral factor

dm = dset_fwd.variables['vp'][:] - dset_fwd.variables['vp1D'][:]
m = dset_fwd.variables['vp1D'][:]
kernels = dset_ker.variables["Kernels"][0]
nu = dset_fwd.variables['Nus'][:]
nu_sum = np.cumsum(np.concatenate(([0],nu)))[:-1]
nuKer = dset_ker.variables['Nus'][:]
nuKer_sum = np.cumsum(np.concatenate(([0],nuKer)))[:-1]
nr = dset_fwd.variables['Nrs'][:]
nrKer = dset_ker.variables['Nrs'][:]
int_fact = dset_fwd.variables['integral_factor'][:]

fem_mesh_fwd = dset_fwd.variables['element_mesh'][:]
fem_mesh_ker = dset_fwd1D.variables['element_mesh'][:]

nelem = len(nu)

#a = 0.306567

##### FFT, MULTIPLY, INTEGRATE ####

chi_tot = 0

for ielem in range(nelem):
    
    ielem_fwd = np.argwhere(fem_mesh_fwd==ielem)[0,0]    
    ielem_ker = np.argwhere(fem_mesh_ker==ielem)[0,0] 
    nuelem = min(nu[ielem_fwd],nuKer[ielem_ker])
    nrelem = min(nr[ielem_fwd],nrKer[ielem_ker])
    nulineFwd = nu_sum[ielem_fwd]
    nulineKer = nuKer_sum[ielem_ker]
    chi_elem = 0
    ker_c = kernels[nulineKer:nulineKer+nuelem, 0, :, :] + 1j * kernels[nulineKer:nulineKer+nuelem, 1, :, :]
    dm_c = dm[nulineFwd:nulineFwd + nuelem, 0, :, :] + 1j * dm[nulineFwd:nulineFwd+nuelem, 1, :, :]
    int_elem = 2*np.pi*int_fact[ielem_fwd,:,:]
    
    
    for inu in range(nuelem):
        chi_elem += np.sum( 2*np.real(np.multiply(int_elem,np.multiply(ker_c[inu,:,:], dm_c[inu,:,:])))) #because negative nu are complex conjugate of positive, it simplifies to 2*np.real()

    chi_tot +=(chi_elem)

print(chi_tot)

