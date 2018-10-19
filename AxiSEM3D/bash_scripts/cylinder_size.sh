#!/bin/bash

### script to increase size of cylinder and move output to right name 


n_count=20
i_count=0
size=100.
size_new=100.
i_size=30. # size increment
input_dir_run="/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/input/"  #where the input files are
exec_file="/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/axisem3d" #executable to run 
dir_seis="/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/stations/" #where the outputs are
mfile_name="${input_dir_run}inparam.model"

function misfit_iter () {
    
    
    sed -i "s/ref1D\$.05\$${size}/ref1D\$.05\$${size_new}/" $mfile_name    # change value of perturbation in model file
    mpirun -np 20 $exec_file >& "axisem3d_output_${size_new}km.txt" #run on 20 cores and print output to file
    mv "${dir_seis}axisem3d_synthetics.nc" "${dir_seis}axisem3d_synthetics_traveltime_0.05vp_${size_new}km_dt0.01.nc" #move seismogram to correct name
    
}

while [ $i_count -le $n_count ]
do
    size_new=`echo $size + $i_size | bc`
    echo Running AxiSEM3D with Vp perturbation cylinder of radius $size_new km
    misfit_iter
    size=$size_new
    let "i_count+=1"
done



