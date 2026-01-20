#!/bin/bash
#SBATCH --job-name=cluster_analysis_part_1_batch1_ser
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --time=72:00:00
#SBATCH --constraint=LSDF

module load math/matlab/R2024b
srun matlab -nodisplay -batch "ClusterAnalysis_mESC_Oct2025_part1_batch1_ser"
