#!/bin/bash

projectdir=${HOME}"/software/growth_SM2018/"

exec_dir=${projectdir}"/bin/"

function runOneSimulation(){

    tag=$1
    s1=$2 # by default s1 is the radial growth direction (because angle is defaulted to 0)
    s2=$3 # by default s2 is the circumferential growth direction
    
    otherargs=$4

    if [ ! -d "$tag" ]; then
        mkdir ${tag}
    else
        echo "SKIPPING "${tag}" because directory exists already"
        exit
    fi

    echo "running "${tag}

    cd "$tag"
    cp ../shell .
    ./shell -sim monolayer_growth -case basic_disk -growthcase ortho -res 512 -h 0.01 -innerR 0.05 -s1 $s1 -s2 $s2 $otherargs
    cd ../
}

# move to the executable directory

cd ${exec_dir}
runOneSimulation radial_growth 0.05 0.0
runOneSimulation circum_growth 0.0 0.05
