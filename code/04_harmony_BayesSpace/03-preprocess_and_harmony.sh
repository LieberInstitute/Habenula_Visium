#!/bin/bash

#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 4
#$ -N preprocess_and_harmony
#$ -o ../../processed-data/04_harmony_BayesSpace/logs/03-preprocess_and_harmony.log
#$ -e ../../processed-data/04_harmony_BayesSpace/logs/03-preprocess_and_harmony.log

#SBATCH -q shared
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH --job-name=preprocess_and_harmony
#SBATCH -o ../../processed-data/04_harmony_BayesSpace/logs/03-preprocess_and_harmony.log
#SBATCH -e ../../processed-data/04_harmony_BayesSpace/logs/03-preprocess_and_harmony.log

if [[ ! -z $SLURMD_NODENAME ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
    node_name=$SLURMD_NODENAME
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
    node_name=$HOSTNAME
fi

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

## List current modules for reproducibility
module load conda_R/4.3.x
module list

Rscript 03-preprocess_and_harmony.R

echo "**** Job ends ****"
date