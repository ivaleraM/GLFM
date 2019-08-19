#!/bin/bash
module load matlab/R2017a-fasrc01

for Dz in 10 20 50; do
    for Pmissing in 10 20 30 40 50 60 70 80 90; do
        for mask_seed in {1..20}; do
            for db_id in 1 2 3 4 5; do
                sbatch indiv_run.sh $Pmissing $mask_seed $db_id $Dz
            done
        done
    done
done
