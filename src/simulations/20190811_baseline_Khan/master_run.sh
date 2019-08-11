module load matlab/R2017a-fasrc01
Pmissing=$1
  mask_seed=$2
  inCluster=1
  db_id=$3

for Pmissing in 10 20 30 40 50 60 70 80 90; do
    for mask_seed in {1..20}; do
        for db_id in 1 2 3 4 5; do
            sbatch indiv_run.sh $Pmissing $mask_seed $db_id
        done
    done
done
