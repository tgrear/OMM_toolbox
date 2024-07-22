#!/bin/sh

#SBATCH --partition=Gemini --job-name=OMM --nodes=1 --ntasks-per-node=64
#SBATCH --gres=gpu:1 --mem=192gb --time=2:00:00 --output=%x_%j_out.out --error=%x_%j_err.err

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"

cd ${SLURM_SUBMIT_DIR}
module purge
module load anaconda3
conda activate OMM_env
python -B OMM_lysozyme_tutorial.py

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
