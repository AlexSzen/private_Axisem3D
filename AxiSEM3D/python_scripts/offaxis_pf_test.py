########
### compare seismograms from axial source and off axis source
### compute mse
########

### Modules
from netCDF4 import Dataset,Variable
import numpy as np
import matplotlib.pyplot as plt
import copy


### Inputs

INPUT_DIR = '/home/alex/Desktop/phd/private_Axisem3D/AxiSEM3D/build/output/stations/'
OUTPUT_DIR = '/home/alex/Desktop/phd/figures/seismograms/'

STATION = 'II.DAK.ENZ'
STATION2 = 'II.DAK.ENZ'


smgr1 = Dataset(INPUT_DIR + 'axisem3d_synthetics_fwd_archer.nc', 'r')
smgr2 = Dataset(INPUT_DIR + 'axisem3d_synthetics_offaxis_nu1782_a0.005_archer.nc', 'r')


test = smgr2.variables[STATION2]
phi = test.longitude
time = smgr1.variables['time_points'][:]
station = smgr1.variables[STATION][:,2]
station2_temp = smgr2.variables[STATION2][:,2]
station2 = copy.deepcopy(station2_temp)


tot_energy = np.sum(np.abs(station)) 
diff_energy = np.abs(np.abs(station)-np.abs(station2))
diff_energy = np.sum(diff_energy)
misfit = diff_energy/tot_energy
print(misfit)
