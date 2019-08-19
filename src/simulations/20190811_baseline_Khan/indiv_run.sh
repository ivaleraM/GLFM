#!/bin/bash
#SBATCH -t 0-48:00 # Runtime in D-HH:MM
#SBATCH -p shared # Partition to submit to | serial_requeue | shared |doshi-velez
#SBATCH --mem=16000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o /n/home11/mfernandezpradier/logs/glfm/log_%j.out # File to which STDOUT will be written
#SBATCH -e /n/home11/mfernandezpradier/logs/glfm/log_%j.err # File to which STDERR will be written

Pmissing=$1
mask_seed=$2
inCluster=1
db_id=$3
Dz=$4
matlab -nodisplay -nosplash -nodesktop -r "run_baseline_mixed_FA($Pmissing, $mask_seed, $inCluster, $db_id, $Dz); exit;"
