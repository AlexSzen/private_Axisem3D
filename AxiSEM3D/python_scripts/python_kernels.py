''' this is a test for tabs'''

from netCDF4 import Dataset
import numpy as np 
from mpi4py import MPI 


### INPUTS

INPUT_FILE_FWD = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/wavefields/wavefield_db_fwd.nc4'
INPUT_FILE_ADJ = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/wavefields/wavefield_db.nc4'

START_TSTEP = 0
END_TSTEP = 1
INT_TSTEP = 1


### start up MPI 

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

### load files. Only master will load the wavefields, then scatterv to workers. 

dset_fwd = Dataset(INPUT_FILE_FWD)
dset_adj = Dataset(INPUT_FILE_ADJ)
nu_fwd = dset_fwd.variables['Nus'][:]
nu_sum_fwd = np.cumsum(np.concatenate(([0],nu_fwd)))[:-1]
nu_adj = dset_adj.variables['Nus'][:]
nu_sum_adj = np.cumsum(np.concatenate(([0],nu_adj)))[:-1]
nelem = len(nu_fwd)

nelem_proc = nelem // nprocs
if (rank == nprocs - 1 and nelem%nprocs != 0):
	nelem_proc = nelem%nprocs 


count_elem_allprocs_temp = np.empty(nprocs) 	
count_elem_allprocs_temp = comm.gather(nelem_proc, root=0)
if (rank == 0):
	start_elem_allprocs = np.cumsum(np.concatenate(([0],count_elem_allprocs_temp)))[:-1].astype('int')
	count_elem_allprocs = np.asarray(count_elem_allprocs_temp).astype('int')
else:
	start_elem_allprocs=np.empty(nprocs,dtype='int')
	count_elem_allprocs=np.empty(nprocs, dtype='int')

comm.Bcast(count_elem_allprocs, root = 0)
comm.Bcast(start_elem_allprocs, root = 0)

start_nuelemfwd_allprocs = nu_sum_fwd[start_elem_allprocs]
count_nuelemfwd_allprocs = nu_fwd[start_elem_allprocs+count_elem_allprocs-1]+nu_sum_fwd[start_elem_allprocs + count_elem_allprocs-1] - nu_sum_fwd[start_elem_allprocs]

start_nuelemadj_allprocs = nu_sum_adj[start_elem_allprocs]
count_nuelemadj_allprocs = nu_adj[start_elem_allprocs+count_elem_allprocs-1]+nu_sum_adj[start_elem_allprocs + count_elem_allprocs-1] - nu_sum_adj[start_elem_allprocs]

split_sizes_input_fwd = count_nuelemfwd_allprocs * 6 * 5 * 5




if (rank==0):
    wavefield = dset_fwd.variables["displacement_wavefield"][END_TSTEP:START_TSTEP:-INT_TSTEP] #domain decomposed wavefield.
    wavefield_adj = dset_adj.variables["displacement_wavefield"][START_TSTEP:END_TSTEP:INT_TSTEP]
    sem_mesh = dset_fwd.variables["sem_mesh"][:] #for each point contains global point tag
    s_temp = dset_fwd.variables["mesh_S"][:]
    z_temp = dset_fwd.variables["mesh_Z"][:]

    vp_temp = dset_fwd.variables['vp1D'][:]
    vp = np.expand_dims(vp_temp, axis = 0) # allows to treat it as a wvf
    rho_temp = dset_fwd.variables['rho1D'][:]
    rho = np.expand_dims(rho_temp, axis = 0) # allows to treat it as a wvf
    int_factor = dset_fwd.variables['integral_factor'][:]
    npoints = np.max(sem_mesh) + 1
    s = np.zeros(npoints, dtype = np.float32)
    z = np.zeros(npoints, dtype = np.float32)
    tim = dset_fwd.variables["time"][:]
    dt = tim[1] - tim[0]

else:
	wavefield = None
	wavefield_adj = None 
	

### master sends wavefields to workers 


	
	
	
	
	




















