#!/bin/bash -e

base_dir=$PWD

declare -A data_dirs
data_dirs[es]="../gen_trj/es/308204"
data_dirs[conE]="../gen_trj/conE/308205"

declare -A lmp_inputs
lmp_inputs[es]=lmp_scripts/es.lmp
lmp_inputs[conE]=lmp_scripts/conE.lmp

for name in es conE; do
    i=0
    mkdir -p $name
    data_dir=${data_dirs[$name]}
    [ -d $data_dir ] || { echo "Error: $data_dir is not a directory"; exit 1; }
    for data_file in ${data_dir}/data_*.lammps; do
        data_name=${data_file##*/}
        number=${data_file##*_}
        number=${number%%.*}
        dir=$name/trj_$number
        [ ! -d $dir ] || continue
        # [ "$i" -lt 100 ] || break
        i=$((i+1))
        echo "Creating $dir"
        mkdir $dir
        cp $data_file $dir
        cp $base_dir/lmp_scripts/settings.lmp $dir
        cp $base_dir/inv.csv $dir
        cp $base_dir/${lmp_inputs[$name]} $dir/run.lmp
        settings_file=$dir/step_settings.lmp
        echo "variable data_file string ${data_name}" > $settings_file
        echo "variable seed equal ${number}" >> $settings_file
    done
done
        #ln -s $base_dir/inv.csv $step_dir/inv.csv 
        #cp lmp_scripts/settings.lmp $step_dir
        #cp ${lmp_inputs[${name}]} "$step_dir/in.lmp"

