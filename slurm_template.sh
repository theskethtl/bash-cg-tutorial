#!/bin/bash
#SBATCH --export=ALL
#SBATCH --output=PEPTIDE_slurm.log
#SBATCH --job-name=PEPTIDE_eq_coarse
#SBATCH --account=teaching
#SBATCH --partition=teaching
#SBATCH --time="06:00:00"
#SBATCH --ntasks=16
module purge
module load gromacs/intel-2018.2/2016.5-single-wee-archie
mpirun -np ${SLURM_NTASKS} gmx_mpi mdrun -ntomp 1 -deffnm PEPTIDE_eq
