#!usr/bin/env python
# -*- coding: utf-8 -*-

'''
netcdf_inversion_input.py

Create netcdf input file ADJOINT_SOURCES for AxiSEM3D inversion.
The file contains variables corresponding to each receiver serving as an adjoint source.
The variable is named after the receiver similar to the seismogram output from AxiSEM3D (e.g. II.AAK.SPZ)
Each variable contains all the measurements for this receiver (rec location, time window, filter info, measurement type, measurement value)

The receivers are loaded from file STATIONS. The rest of the info is loaded from inversion_input.py at the moment

'''

import numpy as np
from netCDF4 import Dataset
import sys
import os
from inversion_input import *



######## CREATE NETCDF FILE #########

### create netcdf file 
nc_adjoint = Dataset(OUTPUT_PATH + 'adjoint_stations.nc4', 'w', 'NETCDF4')

### create dimensions

# one dimension for each measurement number
n_measurement_min = min(NUM_MEASUREMENTS)
n_measurement_max = max(NUM_MEASUREMENTS)

for i in range(n_measurement_min, n_measurement_max + 1):
    nc_adjoint.createDimension('n_measurement_' + str(i), i)

# one dimension for each possible number of params (depends on filter and tomography type)
# for now just consider traveltime/amplitude and log gabor : 5 params 
nc_adjoint.createDimension('station_params_3', 3)
nc_adjoint.createDimension('filter_dim_2', 2)

# write global attributes 
nc_adjoint.measurement_type = MEASUREMENT_TYPE
nc_adjoint.num_sources = NUM_SOURCES
nc_adjoint.num_filters = NUM_FILTERS
nc_adjoint.cp = FILTERS[0][1]
nc_adjoint.sig = FILTERS[0][2]
if MEASUREMENT_TYPE == 'waveform_measurement':
    nc_adjoint.dataPath = DATA_PATH

# load stations file for names and locations 
stations = open(STATIONS_PATH, 'r')

######## CREATE NETCDF FILE #########

####### create filters 

for ifilter in range(NUM_FILTERS):

    filt_name = 'band_' + str(ifilter)
    filt = nc_adjoint.createVariable(filt_name, 'f8', ('filter_dim_2'))
    filt[0] = FILTERS[ifilter][1]
    filt[1] = FILTERS[ifilter][2]
    filt.type = FILTERS[ifilter][0]
    filt.cp = FILTERS[ifilter][1]
    filt.sig = FILTERS[ifilter][2]
######## Create source variables 

for isource in range(NUM_SOURCES):
    
    station = stations.readline().strip().split() 
    s_name = station[1] + '.' + station[0]
    
    # create variable
    adj_source = nc_adjoint.createVariable(s_name, 'f8', ('n_measurement_' + str(NUM_MEASUREMENTS[isource]), 'station_params_3'))
    

    # write measurements parameters 
    for imeasure in range(NUM_MEASUREMENTS[isource]):
        
        adj_source[imeasure, 0] = WINDOWS[isource][imeasure][0]
        adj_source[imeasure, 1] = WINDOWS[isource][imeasure][1]
        adj_source[imeasure, 2] = MEASUREMENTS[isource][imeasure] 
        filt_measure = "filter_" + str(imeasure)
        adj_source.setncattr(filt_measure, int(MEASUREMENT_FILTERS[isource][imeasure][-1]))
        
nc_adjoint.close()
stations.close()








