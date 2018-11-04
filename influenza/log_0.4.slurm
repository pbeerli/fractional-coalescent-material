#!/bin/bash
#SBATCH -o log_0.4
#SBATCH -n 41
#SBATCH -t 4:00:00
#SBATCH -C intel
echo $MPI_BIN
echo $SLURM_JOB_NODELIST
export PATH=.:~/bin:$PATH
echo $PATH
cd /gpfs/research/beerli/scratch/beerli/migrate-mlf2/mittag-leffler/influenza_rerun_oct27four
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
srun -n 41 ./migrate-mlf-mpi parmfile_0.4 -nomenu
