#!/usr/bin/env python


'''

Measurement inputs for netcdf_inversion_input.py
Here we define all the inputs (except name and location) needed
for adjoint sources : time windows, filters, measurement types,...

'''


########## GENERAL  ##########

STATIONS_PATH = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/input/STATIONS' #path of STATIONS file for AxiSEM3D input. Used to get names and location of adjoint sources
OUTPUT_PATH = '/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/input/' # should be path of AxiSEM3D input
NUM_SOURCES = 1 #Number of adjoint sources. If NUM_SOURCES < number of receivers in STATIONS file, only first NUM_SOURCES will be used. 
NUM_FILTERS = 1
########## GENERAL  ##########

########## NUMBER OF MEASUREMENTS ##########
'''
Number of measurements for each receiver : list of ints
will give first dimension for each netcdf adjoint source file
'''

NUM_MEASUREMENTS = [1,1,1,1,1]

########## NUMBER OF MEASUREMENTS ##########


########## TIME WINDOWS ##########

'''
windows (in seconds) are a list of tuples. 
'''

WINDOWS = [[(560.,760.)],[(580., 720.)],[(160., 260.)], [(150., 170.)],[(150., 170.)],[(150., 170.)]]


########## TIME WINDOWS ##########


########## FILTERS  ##########
'''
filters are a list of tuples : they contain filter type (string) and filter params
'''

FILTERS = [('butter_lowpass', 30., 4)]

########## FILTERS  ##########


########## MEASUREMENTS  ##########
'''
tomography type : traveltime_tomography / amplitude_tomography / waveform_tomography 
Depending on the type, measurements will be a list of 
- traveltime delays (s)
- amplitude differences (?)
- name data to construct misfit. in this case needs to specify DATA_PATH
'''

MEASUREMENT_TYPE = 'traveltime_measurement'

DATA_PATH = ''

MEASUREMENTS = [[2.3], [2.3],[2.3],[2.3],[2.3]]

MEASUREMENT_FILTERS = [['band_0'], ['band_0'], ['band_0'], ['band_0'], ['band_0']]
########## MEASUREMENTS  ##########
