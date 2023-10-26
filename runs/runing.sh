#!/bin/sh
#SBATCH --job-name ground_pca_GR
#SBATCH --chdir /home/changfenggroup/nrui/works/codes/ground_pca/runs
#SBATCH --partition changfeng 
#SBATCH --nodes 1
#SBATCH --ntasks 40
export OMP_NUM_THREADS=1

mpiexec -n 40 python GW190513_205428_GR.py > GW190513_205428_GR.out 2>&1
mpiexec -n 40 python GW190519_153544_GR.py > GW190519_153544_GR.out 2>&1
mpiexec -n 40 python GW190602_175927_GR.py > GW190602_175927_GR.out 2>&1
mpiexec -n 40 python GW190706_222641_GR.py > GW190706_222641_GR.out 2>&1
mpiexec -n 40 python GW190828_063405_GR.py > GW190828_063405_GR.out 2>&1

mpiexec -n 40 python _GW170823_dispersion.py > GW170823_dispersion.out 2>&1
mpiexec -n 40 python _GW190408_181802_dispersion.py > GW190408_181802_dispersion.out 2>&1
mpiexec -n 40 python _GW190503_185404_dispersion.py > GW190503_185404_dispersion.out 2>&1
mpiexec -n 40 python _GW190512_180714_dispersion.py > GW190512_180714_dispersion.out 2>&1
mpiexec -n 40 python _GW190915_235702_dispersion.py > GW190915_235702_dispersion.out 2>&1
