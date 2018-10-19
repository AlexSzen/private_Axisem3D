'''

plots Z component seismograms at receiver location from axial vertical point force on north pole
plot Z component seismograms at north pole from off axis source at receiver location 
By reciprocity, these are equal. By plotting them together, we can assess if the gaussian does a good job
at approximating a point force source 
'''



import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import sys

####### INPUTS #######

INPUT_DIR = '/home/alex/Desktop/VariousAxiSEM3DVersions/kuangdai_axisem/AxiSEM3D/build/output/stations/'

file_axis = 'axisem3d_synthetics_axis_allrec.nc'
file_offaxis = 'axisem3d_synthetics.nc'
file_offaxis_2 = 'axisem3d_synthetics_allrec_0.02.nc'
rec_axis = 'II.IAK.SPZ'
rec_offaxis = 'II.IAK.SPZ'

######## INPUTS #######

d_axis = Dataset(INPUT_DIR + file_axis)
d_offaxis = Dataset(INPUT_DIR + file_offaxis)
d_offaxis2  = Dataset(INPUT_DIR + file_offaxis_2)
time = d_axis.variables['time_points'][:]
time2 = d_offaxis.variables['time_points'][:]
s_axis = d_axis.variables[rec_axis][:]
s_offaxis = d_offaxis.variables[rec_offaxis][:]
s_offaxis2 = d_offaxis.variables[rec_offaxis][:]
plt.plot(time, s_axis[:,2]/max(s_axis[:,2]), 'r')
plt.plot(time2, s_offaxis[:,2]/max(s_offaxis[:,2]), 'k')
#plt.plot(time, s_offaxis2[:,2]/max(s_offaxis2[:,2]))


plt.xlabel('Time (s)')
plt.ylabel('Normalized amplitude')
plt.title('Comparison of seismograms from axial and non axial sources')
plt.legend(['Axial source','Off axis source, a = 0.015'])
plt.show()

