#!/bin/bash

projectdir=${HOME}"/software/growth_SM2018/"

exec_dir=${projectdir}"/bin/"

function runOneSimulation(){

    tag=$1
    logh=$2 #actually -log(h)
    h=`awk -v arg="$logh" 'BEGIN {print 10**arg}'`
    res=$3
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
    ./shell -sim monolayer_growth -case basic_disk -growthcase validation -h $h -res $res -innerR 0.1 -R 1.1 -E 1e6 $otherargs
    cd ../
}

# move to the executable directory

cd ${exec_dir}

minlogh=-2 # h = 0.01
maxlogh=0  # h = 1.0
nh=21

for (( i=0; i < ${nh}; i++ )); do
    logh=$(echo "$minlogh + $i * ($maxlogh - $minlogh)/($nh - 1)" | bc -l)
    echo `awk -v arg="$logh" 'BEGIN {print 10**arg}'`
    runOneSimulation "validation_h"$i $logh 352
done

