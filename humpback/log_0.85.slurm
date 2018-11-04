#!/bin/bash
#SBATCH -o log_0.85
#SBATCH -n 11
#SBATCH -t 4:00:00
#SBATCH -C intel
echo $MPI_BIN
echo $SLURM_JOB_NODELIST
export PATH=.:~/bin:$PATH
echo $PATH
cd /gpfs/research/beerli/scratch/beerli/migrate-mlf2/mittag-leffler/humpback_oct272018_rerun
module purge;
y=$(which ./migrate-mlf-mpi)
x=$(strings $y | grep mvapich | wc -l)
if [ "$x"  = "1" ]; then
    module load gnu-mvapich2
else
    module load gnu-openmpi
fi
echo $SLURM_JOB_NODELIST
pwd
srun -n 11 ./migrate-mlf-mpi parmfile_0.85 -nomenu
