#!/bin/bash -e

base_dir=$PWD
lmp=~/proj/lammps/build/lmp
for system_dir in $base_dir/mw $base_dir/graph_il; do
    [ -d $system_dir ] || { echo "Error"; exit 1; }
    input_dir=${system_dir}/inputs
    [ -d $input_dir ] || { echo "Error"; exit 1; }
    #for dt in 0.05 0.1 0.2 0.4 0.8 1.6 3.2 6.4; do
    for dt in 0.2 0.4 0.8 1.6 3.2 6.4; do
        name=dt_${dt}
        dir=$system_dir/$name
        if [ -d $dir ]; then
            echo "Skipping $dir"
            continue
        fi
        echo $dir
        mkdir -p $dir
        cp $input_dir/* $dir
        echo "variable dt equal ${dt}" >> $dir/settings.lmp
        cd $dir
        mpirun -np 4 $lmp -i run.lmp &
    done
    wait
done
echo "Done"
