#!/bin/bash

# Request resources:
#SBATCH -N 26		# number of compute nodes. 
#SBATCH -n 1024		# number of MPI ranks (1 per CPU core)
#SBATCH --mem=32G	# memory required per node, in units M, G or T
#SBATCH --time=0:15:0	# time limit for job (format:  days-hours:minutes:seconds)

# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared 

# Modules necessary for job:
module purge
module load gcc openmpi

# compile part2.c into part2
mpicc -lm part2.c -o part2

# run part2 with 4 processes (from #SBATCH -n, above)
time mpirun ./part2
