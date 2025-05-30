#!/bin/bash
#SBATCH --account=def-jdancker
#SBATCH --time=0-00:07
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=1
#SBATCH --array=401-960
#SBATCH --output=integration_batch.out
#path to config
config=/project/def-jdancker/n24taylo/HCP_neuroimaging/code/config2.txt
subnum=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
#load matlab
module load matlab/2021a
#subject ID to run
srun matlab -nodisplay -singleCompThread -r "integration_supercomputer($subnum,1,0.75)"