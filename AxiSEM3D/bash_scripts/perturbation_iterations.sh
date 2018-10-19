#!/bin/bash

# script to change the value of perturbation in inparam model
# launch an axisem3D run, then move synthetics to correct name

n_count=100
i_count=0
perturbation=.00
perturb_new=.00
i_perturb=.001 # perturbation increment
input_dir_run="/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/input/"  #where the input files are
exec_file="/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/axisem3d" #executable to run 
dir_seis="/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/stations/" #where the outputs are
mfile_name="${input_dir_run}inparam.model"

function misfit_iter () {
    
    
    sed -i "s/ref1D\$${perturbation}/ref1D\$${perturb_new}/" $mfile_name    # change value of perturbation in model file
    mpirun -np 20 $exec_file >& "axisem3d_output_${perturb_new}.txt" #run on 20 cores and print output to file
    mv "${dir_seis}axisem3d_synthetics.nc" "${dir_seis}axisem3d_synthetics_traveltime_${perturb_new}.nc" #move seismogram to correct name
    
}

while [ $i_count -le $n_count ]
do
    perturb_new=`echo $perturbation + $i_perturb | bc`
    echo Running AxiSEM3D with Vp perturbation $perturb_new
    misfit_iter
    perturbation=$perturb_new
    let "i_count+=1"
done



