#!/bin/sh
#Aggregator_Habanero.sh

#Slurm directives
#
#SBATCH -A sscc
#SBATCH -J FG_aggregator
#SBATCH -c 4
#SBATCH -t 6:00:00
#SBATCH --mem-per-cpu=24gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fg2465@columbia.edu
#SBATCH --array=1-12

module load R

i=${SLURM_ARRAY_TASK_ID}

#Command to execute R code
R CMD BATCH --no-save --vanilla AggregatorByC_Habanero.R rout_aggregator$SLURM_ARRAY_TASK_ID

# End of script
