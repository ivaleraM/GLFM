#!/bin/bash
module load matlab/R2017a-fasrc01

for Dz in 50; do
    for Pmissing in 90; do
        for mask_seed in {1..1}; do
            for db_id in 4; do
                sbatch indiv_run.sh $Pmissing $mask_seed $db_id $Dz
            done
        done
    done
done
