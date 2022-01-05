#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=BOLT-REML.refine
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=167:00:00
#SBATCH --mem=400G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --constraint=Intel

cd /work/kylab/mike/PUFA-GWAS/BOLT-REML
#https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.3.4_manual.pdf

bgendir=("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/genotypeQC")
bfiledir=("/scratch/mf91122/PUFA-GWAS/Prune/4.bedfinal")
phenodir=("/scratch/mf91122/PUFA-GWAS/pheno")
outdir=("/scratch/mf91122/PUFA-GWAS/BOLT-REML/norefine")
boltdir=("/home/mf91122/BOLT-LMM/BOLT-LMM_v2.3.4")

mkdir -p $outdir/$j

$boltdir/bolt \
--bed=$bfiledir/chr{1..22}.bed \
--bim=$bfiledir/chr{1..22}.bim \
--fam=$bfiledir/chr1.fam \
--geneticMapFile=$boltdir/tables/genetic_map_hg19_withX.txt.gz \
--phenoFile="$phenodir"/PUFA_GWAS_pheno_M1.txt.pc.txt \
--phenoCol=w3FA_NMR_resinv \
--phenoCol=w6FA_NMR_resinv \
--phenoCol=DHA_NMR_resinv \
--phenoCol=LA_NMR_resinv \
--phenoCol=PUFA_NMR_resinv \
--phenoCol=MUFA_NMR_resinv \
--covarFile="$phenodir"/PUFA_GWAS_pheno_M1.txt.pc.txt \
--qCovarCol=PC{1:20} \
--qCovarCol=Age \
--qCovarCol=Age2 \
--qCovarCol=center{1:18} \
--qCovarCol=center{20:21} \
--qCovarCol=Geno_batch \
--covarCol=Sex \
--numThreads=16 \
--LDscoresFile="$boltdir"/tables/LDSCORE.1000G_EUR.tab.gz \
--reml \
--remlNoRefine
