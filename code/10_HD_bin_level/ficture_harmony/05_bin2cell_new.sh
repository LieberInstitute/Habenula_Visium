#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=05_bin2cell
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH --array=1-5%5

MY_VAR=${MY_VAR:-undefined}
LOG_DIR=../../../processed-data/10_HD_bin_level/ficture_harmony/logs
mkdir -p "$LOG_DIR"
OUT_LOG="$LOG_DIR/05_bin2cell_loop_s${SLURM_ARRAY_TASK_ID}_k${MY_VAR}.txt"
echo "Redirecting to $OUT_LOG"
exec > "$OUT_LOG" 2>&1

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "K: ${MY_VAR}"

module load visium_hd/1.0
## List current modules for reproducibility
module list

python 05_bin2cell_new.py

echo "**** Job ends ****"
date