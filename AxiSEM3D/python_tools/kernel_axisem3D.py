''' simple computation of a kernel on a slice from axisem3d wavefields.
 objective is to compare to same from axisem, see if wavefields are wrong or
 computation is wrong (maybe both...)'''


from wavefield_utils import WavefieldComputer
from sp_utils import *
from params_kernels import *
import numpy as np
import pyvtk
import sys
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#---------- <LOAD VARIABLES> ----------
dset_fwd = Dataset(INPUT_FILE_FWD)
dset_bwd = Dataset(INPUT_FILE_BWD)

wavefield_fwd = dset_fwd.variables["displacement_wavefield"][START_TSTEP:END_TSTEP:INT_TSTEP] #domain decomposed wavefield.
sem_mesh_fwd = dset_fwd.variables["sem_mesh"][:] #for each point contains global point tag
s_temp_fwd = dset_fwd.variables["mesh_S"][:]
z_temp_fwd = dset_fwd.variables["mesh_Z"][:]
nu_fwd = dset_fwd.variables['Nus'][:]
time_temp = dset_fwd.variables['time'][START_TSTEP:END_TSTEP:INT_TSTEP] 
n_pad = 3 * len(time_temp)
n_buff = 100
time = pad(time_temp, n_buff, n_pad)
dt = time_temp[1] - time_temp[0]
time[:n_buff] = time_temp[0] - np.flip(np.asarray(range(n_buff))*dt,0)
time[n_buff:] = time_temp[0] + np.asarray(range( len(time_temp) + n_pad)) * dt
wavefield_bwd = dset_bwd.variables["displacement_wavefield"][START_TSTEP:END_TSTEP:INT_TSTEP] #domain decomposed wavefield.
sem_mesh_bwd = dset_bwd.variables["sem_mesh"][:] #for each point contains global point tag
nu_bwd = dset_bwd.variables['Nus'][:]

npoints = np.max(sem_mesh_fwd) + 1 
nelem = len(nu_fwd)
s = np.zeros(npoints, dtype = np.float32)
z = np.zeros(npoints, dtype = np.float32)
for ielem in range(nelem):
    s[sem_mesh_fwd[ielem]] = s_temp_fwd[ielem]
    z[sem_mesh_fwd[ielem]] = z_temp_fwd[ielem]


num_steps = len(time)

#---------- </LOAD VARIABLES> ----------

 
#---------- <DEFINE SP> ----------
taper = cosine_taper(len(time_temp))
freq = time_to_freq(time)
filt = log_gabor_filter(freq, 2*np.pi/70., 0.7)
#---------- </DEFINE SP> ----------

wc_fwd = WavefieldComputer(wavefield_fwd, nu_fwd, s, z, sem_mesh_fwd)
wc_bwd = WavefieldComputer(wavefield_bwd, nu_bwd, s, z, sem_mesh_bwd)

xdum,ydum,zdum,wvf_fwd_slice = wc_fwd.compute_slice(-1, 7e10, 0., 's')
xdum,ydum,zdum,wvf_bwd_slice = wc_bwd.compute_slice(-1, 7e10, 0., 's')
n_points_slice = len(xdum)
ker_time = np.zeros((num_steps,n_points_slice), dtype = np.float32)
ker_time_sum = np.zeros(n_points_slice, dtype = np.float32)
report = npoints/10
for ip in range(n_points_slice):
    
    if (ip%report == 0):
        print(str(float(ip)/n_points_slice * 100.) + ' % done')
    trace_fwd = taper * wvf_fwd_slice[:,ip]
    trace_bwd = taper * wvf_bwd_slice[:,ip]

    trace_fwd_pad = pad(trace_fwd, n_buff, n_pad)
    trace_bwd_pad = pad(trace_bwd, n_buff, n_pad)

    vel_fwd_filt = freq * filt *1.0j * np.fft.rfft(trace_fwd_pad)
    vel_bwd_filt = freq * 1.0j * np.fft.rfft(trace_bwd_pad)

    
    ker_freq = vel_bwd_filt * vel_fwd_filt
    ker_time[:,ip] = np.fft.irfft(ker_freq)

for it in range(num_steps):
    if (time[it]>150. + n_buff*dt and time[it]<190. + n_buff*dt):
        ker_time_sum += ker_time[it,:]

if not os.path.exists(OUTPUT_DIR+"kernels"):
    os.makedirs(OUTPUT_DIR+"kernels")



points_ker = list(zip(xdum,ydum,zdum))
range_ker = range(n_points_slice)
vtk = pyvtk.VtkData(
        pyvtk.UnstructuredGrid(points_ker,range_ker),
        pyvtk.PointData(pyvtk.Scalars(ker_time_sum,name='ker_test')),
        'animation'
        )
vtk.tofile(OUTPUT_DIR + 'kernels/' +'_kertest_sum' + '.vtk')


#for it in range(100):
#    vtk = pyvtk.VtkData(
#            pyvtk.UnstructuredGrid(points_ker,range_ker),
#            pyvtk.PointData(pyvtk.Scalars(ker_time[it],name='ker_test')),
#            'animation'
#            )
#    vtk.tofile(OUTPUT_DIR + 'kernels/' +'_kertest_' + str(it) + '.vtk')



#---------- </LOAD VARIABLES> ----------



