#!/bin/bash

# Request resources:
#SBATCH -N 1		# number of compute nodes. 
#SBATCH -c 4		# number of CPU cores, one per thread, up to 128
#SBATCH --mem=1G	# memory required, up to 250G on standard nodes
#SBATCH --time=0:35:0	# time limit for job (format:  days-hours:minutes:seconds)

# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared

# Modules necessary for job:
module purge
module load gcc

# compile part1.c into part1
gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=1 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=2 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=3 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=4 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=10 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=20 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=30 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=40 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=50 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=60 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=70 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=80 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=90 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=100 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=110 time ./part1

gcc -fopenmp -lm part1.c -o part1
OMP_NUM_THREADS=128 time ./part1


# # run part1 with 4 threads
# OMP_NUM_THREADS=128 time ./part1

# OMP_NUM_THREADS=1 time ./part1
