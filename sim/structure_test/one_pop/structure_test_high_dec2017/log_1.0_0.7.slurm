#!/bin/bash
#SBATCH -o log_1.0_0.7
#SBATCH -n 21
#SBATCH -t 40:00:00
#SBATCH -C intel
echo $MPI_BIN
echo $SLURM_JOB_NODELIST
export PATH=.:~/bin:$PATH
echo $PATH
cd /panfs/storage.local/scs/beerli/scratch/beerli/migrate-mlf2/mittag-leffler/structure_test_high_dec2017
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
srun -n 21 ./migrate-mlf-mpi parmfile_1.0_0.7 -nomenu
