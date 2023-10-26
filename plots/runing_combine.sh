#!/bin/sh
#SBATCH --job-name ground_pca_combining
#SBATCH --chdir /home/changfenggroup/nrui/works/codes/ground_pca/plots
#SBATCH --partition changfeng 
#SBATCH --nodes 1
#SBATCH --cpus-per-task 40
export OMP_NUM_THREADS=1
mpiexec -n 40 python combine.py > combine.out 2>&1
