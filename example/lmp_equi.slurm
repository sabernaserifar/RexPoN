#!/bin/csh
#SBATCH -J NaCl
#SBATCH -o lmp_equi.out
#SBATCH -e lmp_equi.err
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 10:00:00
##SBATCH --mail-user=naserifar.saber@gmail.com
##SBATCH --mail-type=END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH -A projectnumber
pwd; hostname; date
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"= $SLURM_JOB_NODELIST
srun /home/naseri/codes/src-REQM/lmp_intel < lammps.in
date
