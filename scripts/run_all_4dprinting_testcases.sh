#!/bin/bash

projectdir=${HOME}"/software/growth_SM2018/"

printpath_dir=${projectdir}"/printpaths/"
exec_dir=${projectdir}"/bin/"

function runOneSimulation(){

    tag=$1
    fname1=$2
    fname2=$3
    mu1=$4
    otherargs=$5

    fdir=${tag}

    if [ ! -d "$fdir" ]; then
	mkdir ${fdir}
    else
	echo "SKIPPING "${tag}" because directory exists already"
	exit
    fi
    
    echo "running "${tag}

    cd "$fdir"
    cp ../shell .
    ./shell -sim bilayer_4dfilaments -case ${tag}"_printpath" -fname_bot ${printpath_dir}${fname1} -fname_top ${printpath_dir}${fname2} -botswellfac ${mu1} ${otherargs}
    cd ../
}

# move to the executable directory
cd ${exec_dir}

# run each one of the simulation cases
runOneSimulation helicoid helicoid_vec_1.dat helicoid_vec_2.dat 1.0
runOneSimulation catenoid catenoid_vec_1.dat catenoid_vec_2.dat 1.025
runOneSimulation sombrero sombrero_vec_1.dat sombrero_vec_2.dat 1.2
runOneSimulation folding_flower folding_flower_2_vec_1.dat folding_flower_2_vec_2.dat 1.2
runOneSimulation logspiral logspiral_vec_1.dat logspiral_vec_2.dat 1.0
runOneSimulation orchid orchid_vec_flipshort_1.dat orchid_vec_flipshort_2.dat 1.0 "-simulateplate true -planeflip true -topswellfac 1.2 -res 0.25" 

