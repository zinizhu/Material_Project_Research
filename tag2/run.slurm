#!/bin/bash -l 
 #SBATCH -J tag2 
 #SBATCH -q regular 
 #SBATCH -N 2 
 #SBATCH -t 02:30:00 
#SBATCH -C haswell 

 module load vasp 
srun -n 64 ../vasp_std 
 module load vasp 

 srun -n 64 ../vasp_std