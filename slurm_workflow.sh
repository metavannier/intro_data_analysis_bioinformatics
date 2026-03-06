#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type END
#SBATCH --mail-user foo.bar@france-bioinformatique.fr
#
#SBATCH --partition fast
#SBATCH -N 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32GB

module load snakemake

snakemake  --use-conda --use-singularity --cores $SLURM_CPUS_PER_TASK
