#!/bin/bash
#SBATCH --job-name=PLINKpca
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --time=120:00:00
#SBATCH --output=PLINKpca.%j.out
#SBATCH --error=PLINKpca.%j.err

cd /work/kylab/mike/PUFA-GWAS/PLINKPCA

genoindir=("/scratch/mf91122/PUFA-GWAS/Prune-improved/4.bedfinal")
outdir=("/scratch/mf91122/PUFA-GWAS/PLINKPCA/outpc")

ml PLINK/2.00-alpha2.3-x86_64-20210920-dev

mkdir -p $outdir

plink2 \
--bfile $genoindir/mergedbed \
--pca approx 20 \
--out $outdir/fullpca
