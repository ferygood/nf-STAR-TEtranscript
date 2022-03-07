#!/bin/bash
#SBATCH -J nextflow
#SBATCH -D .
#SBATCH -o nf.%j.out
#SBATCH --partition=nowick
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=750M
#SBATCH --time=10-00:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=yaochung41@gmail.com

date
hostname

nextflow run nf-STAR.nf
