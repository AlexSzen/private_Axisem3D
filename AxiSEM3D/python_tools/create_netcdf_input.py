'''
create netcdf input for inversion.
template ascii inputs provided by Kasra :
- ffproc.receivers
- ffproc.ampstt.bandi
- kernel_filter.txt
- ffproc.source

others

structure : 

- filters

- events :
    - stations :
        - real
        - synthetics :
            - iterations :
                - measurements : 
                    - filter
                    - windows
                    - etc

'''

INPUT_DIR = '/home/alex/Desktop/VariousAxiSEM3DVersions/kuangdai_axisem/AxiSEM3D/build/input/measurement_files/'
OUTPUT_DIR = '/home/alex/Desktop/VariousAxiSEM3DVersions/kuangdai_axisem/AxiSEM3D/build/input/'

import numpy as np
from netCDF4 import Dataset
import os
import sys 

f_cmt = open(INPUT_DIR + 'CMTSOLUTION', 'r')
f_event = open(INPUT_DIR + 'ffproc.source', 'r')
f_receivers = open(INPUT_DIR + 'ffproc2.receivers', 'r')
f_measurements = open(INPUT_DIR + 'ffproc2.ampstt.band01', 'r')

root_grp = Dataset(OUTPUT_DIR + 'inversion_input.nc4', 'w', format = 'NETCDF4')

########## Define event group ##############
f_cmt.readline()
cols = f_cmt.readline().strip().split()
event_grp = root_grp.createGroup(cols[2])

cols = f_cmt.readline().strip().split()
event_grp.time_shift = cols[2]

cols = f_cmt.readline().strip().split()
event_grp.half_duration = cols[2]

cols = f_cmt.readline().strip().split()
event_grp.latitude = float(cols[1])

cols = f_cmt.readline().strip().split()
event_grp.longitude = float(cols[1])

cols = f_cmt.readline().strip().split()
event_grp.depth = float(cols[1])

cols = f_cmt.readline().strip().split()
event_grp.Mrr = float(cols[1])

cols = f_cmt.readline().strip().split()
event_grp.Mrr = float(cols[1])

cols = f_cmt.readline().strip().split()
event_grp.Mpp = float(cols[1])

cols = f_cmt.readline().strip().split()
event_grp.Mrt = float(cols[1])

cols = f_cmt.readline().strip().split()
event_grp.Mrp= float(cols[1])

cols = f_cmt.readline().strip().split()
event_grp.Mtp= float(cols[1])


########## Define station groups ##############

f_receivers.readline() # junk
for i in range(5): # junk
    f_measurements.readline()
    
for line_receiver in f_receivers:
    
    cols = line_receiver.strip().split()
    
    if (cols == []): # can have empty lines before eof  
        continue 

    station_grp = event_grp.createGroup('station_' + cols[0])
    station_grp.grp = float(cols[1])
    station_grp.station_name = cols[2]
    station_grp.latitude = float(cols[3])
    station_grp.longitude = float(cols[4])
    station_grp.elevation = float(cols[5])
    station_grp.burial = float(cols[6])

    real_grp = station_grp.createGroup('Real_data')
    synth_grp = station_grp.createGroup('Synthetic_data')
    iter_grp = synth_grp.createGroup('iteration_1')
    
    last_pos_measure = f_measurements.tell() # if we arrive at measurement of next station, go back one line
    line_measure = f_measurements.readline()
    
    while (line_measure != ''): 

        cols = line_measure.strip().split()
        
        if (cols[0] != station_grp.grp):
            break
        
        


    
root_grp.close()
f_receivers.close()
f_cmt.close()
f_event.close()






