#!/bin/bash

### script to increase size of cylinder and move output to right name 


n_count=20
i_count=0
size=100
size_new=100
i_size=30 # size increment
input_dir_run="/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/snd_build/input/"  #where the input files are
exec_file="/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/snd_build/axisem3d" #executable to run 
dir_wvf="/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/snd_build/output/wavefields/" #where the outputs are
dir_wvf="/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/snd_build/output/wavefields/"
mfile_name="${input_dir_run}inparam.model"

function misfit_iter () {
    
    
    sed -i "s/ref1D\$.05\$${size}/ref1D\$.05\$${size_new}/" $mfile_name    # change value of perturbation in model file
    mpirun -np 20 $exec_file >& "axisem3d_output_model_${size_new}km.txt" #run on 20 cores and print output to file
    
    mv "${dir_wvf}wavefield_db.nc4" "${dir_wvf}model_${size_new}km_nu200_mesh50s.nc4" #move seismogram to correct name
    killall -9 axisem3d 
}

while [ $i_count -le $n_count ]
do
#    size_new=`echo $size + $i_size | bc`
    let "size_new=size+i_size"
    echo Running AxiSEM3D with Vp perturbation cylinder of radius $size_new km
    misfit_iter
    size=$size_new
    let "i_count+=1"
done



