#!/bin/bash -e

# Standard output and error
#SBATCH -o ./out.%j
#SBATCH -e ./err.%j
#
# working directory
#SBATCH -D ./
#
# partition
#SBATCH -p draco-small


#SBATCH --ntasks 32
#SBATCH --cpus-per-task 1
#SBATCH --time 1-00:00:00
#SBATCH --nodes=1

echo "Running on $HOSTNAME"

source .venv/bin/activate

base_dir=$PWD
extract_dir=$base_dir/graph-il_msd

i=0
for name in es conE; do
    for dir in $base_dir/$name/trj_*; do
        [ -d $dir ] || continue
        job_name=${dir##*/}
        dump_dir=$extract_dir/$name/$job_name
        [ ! -f $dump_dir/msd.pickle || ! -f $dump_dir/count.pickle ] || continue
        mkdir -p $dump_dir
        for job_dir in $dir/*; do
            [ -d $job_dir ] || continue
            echo "### $job_dir ###"
            python3 $base_dir/extract_msd.py $job_dir $dump_dir &
            i=$((i+1))
            if [ "$i" -ge 32 ]; then
                echo "Waiting"
                wait
                i=0
            fi
        done
    done
    wait
done

echo "DONE"
