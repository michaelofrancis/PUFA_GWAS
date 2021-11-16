#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=BOLTM2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=167:00:00
#SBATCH --mem=400G
#SBATCH --output=BOLTM2.%j.out
#SBATCH --error=BOLTM2.%j.err
#SBATCH --constraint=Intel

cd /work/kylab/mike/PUFA-GWAS/BOLT-LMM/M2
#https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.3.4_manual.pdf

bgendir=("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/genotypeQC")
bfiledir=("/scratch/mf91122/PUFA-GWAS/Prune/4.bedfinal")
phenodir=("/scratch/mf91122/PUFA-GWAS/pheno")
outdir=("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/M2")
boltdir=("/home/mf91122/BOLT-LMM/BOLT-LMM_v2.3.4")

#phenotypes=("RankTest" "RankbioavailableTest")

phenotypes=("w3FA_NMR_resinv" "w3FA_NMR_TFAP_resinv" "w6FA_NMR_resinv" 
		"w6FA_NMR_TFAP_resinv" "w6_w3_ratio_NMR_resinv" 
		"DHA_NMR_resinv" "DHA_NMR_TFAP_resinv" "LA_NMR_resinv" 
		"LA_NMR_TFAP_resinv" "PUFA_NMR_resinv" "PUFA_NMR_TFAP_resinv" 
		"MUFA_NMR_resinv" "MUFA_NMR_TFAP_resinv" "PUFA_MUFA_ratio_NMR_resinv")


pheno=("w3FA_NMR_resinv")

m=2

for j in ${phenotypes[@]} 
	do


mkdir -p $outdir/$j

$boltdir/bolt \
--bed=$bfiledir/chr{1..22}.bed \
--bim=$bfiledir/chr{1..22}.bim \
--fam=$bfiledir/chr1.fam \
--geneticMapFile=$boltdir/tables/genetic_map_hg19_withX.txt.gz \
--phenoFile="$phenodir"/PUFA_GWAS_pheno_M2.txt.pc.txt \
--phenoCol="$j" \
--covarFile="$phenodir"/PUFA_GWAS_pheno_M2.txt.pc.txt \
--qCovarCol=PC{1:20} \
--qCovarCol=Age \
--qCovarCol=Age2 \
--qCovarCol=center{1:18} \
--qCovarCol=center{20:21} \
--qCovarCol=Geno_batch \
--covarCol=Sex \
--qCovarCol=BMI \
--qCovarCol=Townsend \
--covarCol=statins \
--bgenFile="$bgendir"/chr{1..22}_b-keep.bgen \
--sampleFile="$bgendir"/chr1_b-keep.sample \
--lmm \
--numThreads=16 \
--LDscoresFile="$boltdir"/tables/LDSCORE.1000G_EUR.tab.gz \
--bgenMinINFO=0.3 \
--statsFile="$outdir"/"$j"/BOLT1-statsFile-m"$m" \
--statsFileBgenSnps="$outdir"/"$j"/BOLT1-statsFile-BgenSnps-m"$m"

done
