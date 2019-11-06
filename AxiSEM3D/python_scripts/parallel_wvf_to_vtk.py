'''
Load wavefield netcdf file and converts to vtk for paraview visualization.
Wavefield has dimensions [time][fouriers and elems][6 : real and imag parts of s p z][ipol][jpol]
Author : Alexandre Szenicer
'''

from netCDF4 import Dataset
import pyvtk
import numpy as np
import time
import sys
import os
from new_wavefield_utils import WavefieldComputer
from wisdom_kernel import wisdom
from params_animation import *
from multiprocessing import Pool


#---------- <INPUT AND OUTPUT> ----------


NUM_CORES = 20

if (SLICES):
    if not os.path.exists(OUTPUT_DIR+"slices"):
        os.makedirs(OUTPUT_DIR+"slices")

if (SHELLS):
    if not os.path.exists(OUTPUT_DIR+"shells"):
        os.makedirs(OUTPUT_DIR+"shells")


#list_files_ker = [f for f in sorted(os.listdir(INPUT_DIR_KER)) if '30deg_kernels_db_' in f]
list_files_ker = [f for f in sorted(os.listdir(INPUT_DIR_KER))]
### this one for 3d perturb
list_files_wvf = [f for f in sorted(os.listdir(INPUT_DIR_WVF)) if 'wavefield_db_' in f and 'fwd' not in f and '3d' not in f]
###this one for other perturb
#list_files_wvf = [f for f in sorted(os.listdir(INPUT_DIR_WVF)) if '3d' in f ]
###this one for kernel i guess 
#list_files_wvf = [f for f in sorted(os.listdir(INPUT_DIR_WVF)) if 'wavefield_db_fwd' in f ]

if (len(list_files_ker) != len(list_files_wvf)):
    raise ValueError('Inconsistent number of wvf and kernel files')
num_files = len(list_files_ker)

list_files_input = list(zip(list_files_ker,list_files_wvf, list(range(num_files)) ))

num_steps = int((END_TSTEP - START_TSTEP)/INT_TSTEP)

shell_name = 'kernel_test'
slice_name = 'kernel_test'
#---------- </INPUT AND OUTPUT> ----------

#---------- <LOAD VARIABLES> ----------


def netcdf2txt(INPUT_FILE):
    INPUT_FILE_KER, INPUT_FILE_WVF, NUM_FILE = INPUT_FILE
    dset_fwd = Dataset(INPUT_DIR_WVF + INPUT_FILE_WVF)
    dset_ker = Dataset(INPUT_DIR_KER + INPUT_FILE_KER)

    wavefield = dset_fwd.variables["displacement_wavefield"][START_TSTEP:END_TSTEP:INT_TSTEP] #domain decomposed wavefield.
    # TODO: need to fix this in C++ code, add a timestep even if integrated kernel, for compatibility.
    kernels = np.expand_dims(dset_ker.variables["Kernels"][:], axis = 0)

    vp_temp = dset_fwd.variables['vp'][:] - dset_fwd.variables['vp1D'][:]
    vp_pert = np.expand_dims(vp_temp, axis = 0) # allows to treat it as a wvf

    #sem_mesh = dset_fwd.variables["sem_mesh"][:] #for each point contains global point tag
    s = dset_fwd.variables["mesh_S"][:]
    z = dset_fwd.variables["mesh_Z"][:]
    nu = dset_fwd.variables['Nus'][:]
    nuKer = dset_ker.variables['Nus'][:]

    #vp_temp = dset_fwd.variables['vp'][:] - dset_fwd.variables['vp1D'][:]
    #vp = np.expand_dims(vp_temp, axis = 0) # allows to treat it as a wvf
    #int_factor = dset_fwd.variables['integral_factor'][:]
    #npoints = np.max(sem_mesh) + 1
    nelem = len(nu)
    #nelem = len(nuKer)
    if (nelem == 0):
        return
    #---------- </LOAD VARIABLES> ----------


    #---------- <CREATE VTK FILE> ----------

    wc = WavefieldComputer(kernels, nuKer, nuKer, s, z)
#    wc = WavefieldComputer(vp_pert, nu, nu, s, z)
#    wc = WavefieldComputer(wavefield,nu,nu,s,z)


    if PLOT_NU:
        x_nu, y_nu, z_nu, nu_slice = wc.compute_nu_slice()
        points_slice = list(zip(x_nu,y_nu,z_nu))
        numpoints_slice = len(x_nu)
    #    f_coords = open(OUTPUT_DIR + 'coords_' + str(NUM_FILE)+'.txt','w')
#        f_coords.write(points_slice)
#        f_coords.close()

        ### Make sure we remove existing files
        if os.path.exists(OUTPUT_DIR +'/slices/'+ 'coords_nu_' + str(NUM_FILE)+'.txt'):
            os.remove(OUTPUT_DIR +'/slices/'+ 'coords_nu_' + str(NUM_FILE)+'.txt')
        if os.path.exists(OUTPUT_DIR +'/slices/'+ 'values_nu_' + str(NUM_FILE)+'.txt'):
            os.remove(OUTPUT_DIR +'/slices/'+ 'values_nu_' + str(NUM_FILE)+'.txt')
            
        np.savetxt(OUTPUT_DIR +'/slices/'+ 'coords_nu_' + str(NUM_FILE)+'.txt', points_slice)
        np.savetxt(OUTPUT_DIR +'/slices/'+ 'values_nu_' + str(NUM_FILE)+'.txt', nu_slice, fmt='%i', newline = ' ')
        
        ### This overwrites the file directly
        with open(OUTPUT_DIR +'/slices/'+ 'num_points_nu_' + str(NUM_FILE)+'.txt','w') as f:
            f.write(str(numpoints_slice))

    for i_slice in range(SLICES):

        x_slice, y_slice, z_slice, wvf_slice = wc.compute_slice(RMIN[i_slice], RMAX[i_slice], PHIS_SLICES[i_slice], COMP_SLICES[i_slice])
        points_slice = list(zip(x_slice,y_slice,z_slice))
        numpoints_slice = len(x_slice)

        ### Even if no values are written (i.e. no points found), 
        ### wvf has a length of 1, so we have to skip it        
        if wvf_slice.shape == (num_steps, 1):
            continue


        for it in range(num_steps):

            ### Make sure we remove existing files
            if os.path.exists(OUTPUT_DIR +'/slices/' + 'coords_slice_' + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
            + str(PHIS_SLICES[i_slice]) + '_'+ str(it) + '_' + str(NUM_FILE) + '.txt'):
                os.remove(OUTPUT_DIR +'/slices/' + 'coords_slice_' + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
                + str(PHIS_SLICES[i_slice]) + '_'+ str(it) + '_' + str(NUM_FILE) + '.txt')
            if os.path.exists(OUTPUT_DIR +'/slices/'+ 'values_slice_' + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
            + str(PHIS_SLICES[i_slice]) + '_'+ str(it) + '_' + str(NUM_FILE) + '.txt'):
                os.remove(OUTPUT_DIR +'/slices/'+ 'values_slice_' + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
                + str(PHIS_SLICES[i_slice]) + '_'+ str(it) + '_' + str(NUM_FILE) + '.txt')   
            if os.path.exists(OUTPUT_DIR +'/slices/'+ 'num_points_slice_' + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
            + str(PHIS_SLICES[i_slice]) + '_'+ str(it) + '_' + str(NUM_FILE) + '.txt'):
                os.remove(OUTPUT_DIR +'/slices/'+ 'num_points_slice_' + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
                + str(PHIS_SLICES[i_slice]) + '_'+ str(it) + '_' + str(NUM_FILE) + '.txt')   
                         
            ### Save values
            np.savetxt(OUTPUT_DIR +'/slices/' + 'coords_slice_' + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
            + str(PHIS_SLICES[i_slice]) + '_'+ str(it) + '_' + str(NUM_FILE) + '.txt', points_slice)
            np.savetxt(OUTPUT_DIR +'/slices/'+ 'values_slice_' + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
            + str(PHIS_SLICES[i_slice]) + '_'+ str(it) + '_' + str(NUM_FILE) + '.txt', wvf_slice, newline = ' ')
            
            ### This overwrites the file directly
            with open(OUTPUT_DIR +'/slices/'+ 'num_points_slice_' + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
            + str(PHIS_SLICES[i_slice]) + '_'+ str(it) + '_' + str(NUM_FILE) + '.txt','w') as f:
                f.write(str(numpoints_slice))

    for i_shell in range(SHELLS):

        x_shell, y_shell, z_shell, wvf_shell = wc.compute_shell(R[i_shell], R_TOLERANCE[i_shell], SAMPLE_DENSITY[i_shell], PHIS_SHELLS[i_shell], INNER_OUTER[i_shell], COMP_SHELLS[i_shell])
        points_shell = list(zip(x_shell,y_shell,z_shell))
        range_shell = range(len(x_shell))
        numpoints_shell = len(x_shell)

        ### Even if no values are written (i.e. no points found), 
        ### wvf has a length of 1, so we have to skip it        
        if wvf_shell.shape == (num_steps, 1):
            continue

        for it in range(num_steps):
            ### Make sure we remove existing files
            if os.path.exists(OUTPUT_DIR + 'shells/' + 'coords_shell_'+ INNER_OUTER[i_shell]+'_' + str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) + '_' + str(NUM_FILE) + '.txt' ):
                os.remove(OUTPUT_DIR + 'shells/' + 'coords_shell_' + INNER_OUTER[i_shell]+'_'+str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) + '_' + str(NUM_FILE) + '.txt' )
            if os.path.exists(OUTPUT_DIR + 'shells/' + 'values_shell_' + INNER_OUTER[i_shell] + '_' + str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) + '_' + str(NUM_FILE) + '.txt' ):
                os.remove(OUTPUT_DIR + 'shells/' + 'values_shell_' + INNER_OUTER[i_shell] + '_'  +str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) + '_' + str(NUM_FILE) + '.txt' )   
            if os.path.exists(OUTPUT_DIR + 'shells/' + 'num_points_shell_' + INNER_OUTER[i_shell] + '_' + str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) + '_' + str(NUM_FILE) + '.txt' ):
                os.remove(OUTPUT_DIR + 'shells/' + 'num_points_shell_' + INNER_OUTER[i_shell] + '_'  +str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) + '_' + str(NUM_FILE) + '.txt' )   
                         
            ### Save values
            np.savetxt(OUTPUT_DIR + 'shells/' + 'coords_shell_' + INNER_OUTER[i_shell] + '_' + str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) + '_' + str(NUM_FILE) + '.txt', points_shell )
            np.savetxt(OUTPUT_DIR + 'shells/' + 'values_shell_' + INNER_OUTER[i_shell] + '_' + str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) + '_' + str(NUM_FILE) + '.txt', wvf_shell, newline = ' ')
            
            ### This overwrites the file directly
            with open(OUTPUT_DIR + 'shells/' + 'num_points_shell_' + INNER_OUTER[i_shell] + '_' + str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) + '_' + str(NUM_FILE) + '.txt', 'w' ) as f:
                f.write(str(numpoints_shell))






    #---------- </CREATE VTK FILE> ----------




if __name__ == "__main__" :

    print("Start")

    #---------- <BEGIN TIMER> ----------
    t0 = time.time()
    #---------- </BEGIN TIMER> ----------

    p = Pool(NUM_CORES)
    p.map(netcdf2txt, list_files_input)



    if PLOT_NU:
        nus_coords_files_txt = [f for f in sorted(os.listdir(OUTPUT_DIR+'/slices/')) if 'coords_nu' in f]
        nus_values_files_txt = [f for f in sorted(os.listdir(OUTPUT_DIR+'/slices/')) if 'values_nu' in f]
        nus_numpoints_files_txt = [f for f in sorted(os.listdir(OUTPUT_DIR+'/slices/')) if 'num_points_nu' in f]

        nu_vtk_path = OUTPUT_DIR+'/slices/' + 'slice_nu_wisdom.vtk'

        if os.path.exists(nu_vtk_path):
            os.remove(nu_vtk_path)

        with open(nu_vtk_path, "a") as f_vtk:

            ### Write VTK specific headers for compatibility
            f_vtk.write("# vtk DataFile Version 2.0 \n")
            f_vtk.write("animation \n")
            f_vtk.write("ASCII \n")
            f_vtk.write("DATASET UNSTRUCTURED_GRID \n")

            ### Get total number of points
            num_points = 0
            for numpoints_file in nus_numpoints_files_txt:
                num_points += np.loadtxt(open(OUTPUT_DIR+'/slices/'+ numpoints_file))

            f_vtk.write("POINTS " + str(int(num_points)) + " double \n")


            for coord_file in nus_coords_files_txt:
                f_vtk.write(open( (OUTPUT_DIR+'/slices/'+ coord_file) ).read())

            ### More VTK headers
            f_vtk.write("CELLS " + str(int(num_points)) + " " + str(int(2*num_points)) + '\n' )

            ### Weird VTK thing with range of points and cell type
            ### I create list of tuples (1, range_pos) and write that
            list_tuples = [(1, range_pos) for range_pos in range(int(num_points))]
            f_vtk.write('\n'.join('%s %s' % x for x in list_tuples))

            ### More VTK headers
            f_vtk.write("\nCELL_TYPES " + str(int(num_points)) + '\n')
            np.savetxt(f_vtk, np.ones(int(num_points)), fmt = '%i', newline = ' ')

            ### More VTK headers
            f_vtk.write("\nPOINT_DATA " + str(int(num_points)) + '\n')
            f_vtk.write("SCALARS nu_slice int 1 \n")
            f_vtk.write("LOOKUP_TABLE default \n")

            for value_file in nus_values_files_txt:
                f_vtk.write(open( (OUTPUT_DIR+'/slices/'+ value_file) ).read())

    for i_slice in range(SLICES):
        for it in range(num_steps):
                slice_coords_files_txt = [f for f in sorted(os.listdir(OUTPUT_DIR+'/slices/')) if 'coords_slice_'
                + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
                + str(PHIS_SLICES[i_slice]) + '_'+ str(it) in f]
                slice_values_files_txt = [f for f in sorted(os.listdir(OUTPUT_DIR+'/slices/')) if 'values_slice_'
                + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
                + str(PHIS_SLICES[i_slice]) + '_'+ str(it) in f]
                slice_numpoints_files_txt = [f for f in sorted(os.listdir(OUTPUT_DIR+'/slices/')) if 'num_points_slice_'
                + str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_'
                + str(PHIS_SLICES[i_slice]) + '_'+ str(it) in f]

                ### If I split the name on several lines I get an error wtf??
                slice_vtk_path = OUTPUT_DIR + 'slices/'+slice_name+ str(int(RMIN[i_slice]*1.e-3)) + '_' + str(int(RMAX[i_slice]*1.e-3)) + '_' + str(PHIS_SLICES[i_slice]) + '.vtk'


                if os.path.exists(slice_vtk_path):
                    os.remove(slice_vtk_path)

                with open(slice_vtk_path, "a") as f_vtk:

                    ### Write VTK specific headers for compatibility
                    f_vtk.write("# vtk DataFile Version 2.0 \n")
                    f_vtk.write("animation \n")
                    f_vtk.write("ASCII \n")
                    f_vtk.write("DATASET UNSTRUCTURED_GRID \n")

                    ### Get total number of points
                    num_points = 0
                    for numpoints_file in slice_numpoints_files_txt:
                        num_points += np.loadtxt(open(OUTPUT_DIR+'/slices/'+ numpoints_file))

                    f_vtk.write("POINTS " + str(int(num_points)) + " double \n")


                    for coord_file in slice_coords_files_txt:
                        f_vtk.write(open( (OUTPUT_DIR+'/slices/'+ coord_file) ).read())

                    ### More VTK headers
                    f_vtk.write("CELLS " + str(int(num_points)) + " " + str(int(2*num_points)) + '\n' )

                    ### Weird VTK thing with range of points and cell type
                    ### I create list of tuples (1, range_pos) and write that
                    list_tuples = [(1, range_pos) for range_pos in range(int(num_points))]
                    f_vtk.write('\n'.join('%s %s' % x for x in list_tuples))

                    ### More VTK headers
                    f_vtk.write("\nCELL_TYPES " + str(int(num_points)) + '\n')
                    np.savetxt(f_vtk, np.ones(int(num_points)), fmt = '%i', newline = ' ')

                    ### More VTK headers
                    f_vtk.write("\nPOINT_DATA " + str(int(num_points)) + '\n')
                    f_vtk.write("SCALARS slice double 1 \n")
                    f_vtk.write("LOOKUP_TABLE default \n")

                    for value_file in slice_values_files_txt:
                        f_vtk.write(open( (OUTPUT_DIR+'/slices/'+ value_file) ).read())

    for i_shell in range(SHELLS):
        for it in range(num_steps):
                shell_coords_files_txt = [f for f in sorted(os.listdir(OUTPUT_DIR+'/shells/')) if 'coords_shell_'+ INNER_OUTER[i_shell]+'_' + str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it)  in f]
                shell_values_files_txt = [f for f in sorted(os.listdir(OUTPUT_DIR+'/shells/')) if 'values_shell_'
                + INNER_OUTER[i_shell]+'_' + str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) in f]
                shell_numpoints_files_txt = [f for f in sorted(os.listdir(OUTPUT_DIR+'/shells/')) if 'num_points_shell_'
                + INNER_OUTER[i_shell]+'_' + str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0])
        + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) in f]
                ### If I split the name on several lines I get an error wtf??
                shell_vtk_path = OUTPUT_DIR + 'shells/'+shell_name+ '_'+ INNER_OUTER[i_shell]+'_' + str(int(R[i_shell]*1.e-3)) + '_' + str(PHIS_SHELLS[i_shell][0]) + '_' + str(PHIS_SHELLS[i_shell][1]) + '_' + str(it) + '.vtk'

                if os.path.exists(shell_vtk_path):
                    os.remove(shell_vtk_path)

                with open(shell_vtk_path, "a") as f_vtk:

                    ### Write VTK specific headers for compatibility
                    f_vtk.write("# vtk DataFile Version 2.0 \n")
                    f_vtk.write("animation \n")
                    f_vtk.write("ASCII \n")
                    f_vtk.write("DATASET UNSTRUCTURED_GRID \n")

                    ### Get total number of points
                    num_points = 0
                    for numpoints_file in shell_numpoints_files_txt:
                        num_points += np.loadtxt(open(OUTPUT_DIR+'/shells/'+ numpoints_file))
                    f_vtk.write("POINTS " + str(int(num_points)) + " double \n")


                    for coord_file in shell_coords_files_txt:
                        f_vtk.write(open( (OUTPUT_DIR+'/shells/'+ coord_file) ).read())

                    ### More VTK headers
                    f_vtk.write("CELLS " + str(int(num_points)) + " " + str(int(2*num_points)) + '\n' )

                    ### Weird VTK thing with range of points and cell type
                    ### I create list of tuples (1, range_pos) and write that
                    list_tuples = [(1, range_pos) for range_pos in range(int(num_points))]
                    f_vtk.write('\n'.join('%s %s' % x for x in list_tuples))

                    ### More VTK headers
                    f_vtk.write("\nCELL_TYPES " + str(int(num_points)) + '\n')
                    np.savetxt(f_vtk, np.ones(int(num_points)), fmt = '%i', newline = ' ')

                    ### More VTK headers
                    f_vtk.write("\nPOINT_DATA " + str(int(num_points)) + '\n')
                    f_vtk.write("SCALARS shell double 1 \n")
                    f_vtk.write("LOOKUP_TABLE default \n")

                    for value_file in shell_values_files_txt:
                        f_vtk.write(open( (OUTPUT_DIR+'/shells/'+ value_file) ).read())



    #---------- <END TIMER> ----------
    t2 = time.time()
    print("Total runtime is " + str(t2-t0) + "s")
    #---------- </END TIMER> ----------
