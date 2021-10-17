#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=BOLTM1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=167:00:00
#SBATCH --mem=600G
#SBATCH --output=BOLTM1.%j.out
#SBATCH --error=BOLTM1.%j.err
#SBATCH --constraint=Intel

cd /work/kylab/mike/PUFA-GWAS/BOLT-LMM/M1
#https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.3.4_manual.pdf

bgendir=("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/genotypeQC")
bfiledir=("/scratch/mf91122/PUFA-GWAS/Prune/4.bedfinal")
phenodir=("/scratch/mf91122/PUFA-GWAS/pheno")
outdir=("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/M1")
boltdir=("/home/mf91122/BOLT-LMM/BOLT-LMM_v2.3.4")

#phenotypes=("RankTest" "RankbioavailableTest")

phenotypes=("w3FA_NMR_resinv" "w3FA_NMR_TFAP_resinv" "w6FA_NMR_resinv" 
		"w6FA_NMR_TFAP_resinv" "w6_w3_ratio_NMR_resinv" 
		"DHA_NMR_resinv" "DHA_NMR_TFAP_resinv" "LA_NMR_resinv" 
		"LA_NMR_TFAP_resinv" "PUFA_NMR_resinv" "PUFA_NMR_TFAP_resinv" 
		"MUFA_NMR_resinv" "MUFA_NMR_TFAP_resinv" "PUFA_MUFA_ratio_NMR_resinv")


#pheno=("w3FA_NMR_resinv")

m=1


for j in ${phenotypes[@]} 
	do


mkdir -p $outdir/$j

$boltdir/bolt \
--bfile=$bfiledir/mergedbed \
--phenoFile="$phenodir"/PUFA_GWAS_pheno_M1.txt.pc.txt \
--phenoCol="$j" \
--covarFile="$phenodir"/PUFA_GWAS_pheno_M1.txt.pc.txt \
--qCovarCol=PC{1:20} \
--qCovarCol=Age \
--qCovarCol=Age2 \
--qCovarCol=center{1:18} \
--qCovarCol=center{20:21} \
--qCovarCol=Geno_batch \
--covarCol=Sex \
--bgenFile="$bgendir"/chr{1..22}.bgen \
--sampleFile="$bgendir"/chr1.sample \
--lmm \
--numThreads=16 \
--LDscoresFile="$boltdir"/tables/LDSCORE.1000G_EUR.tab.gz \
--bgenMinINFO=0.3 \
--statsFile="$outdir"/"$j"/BOLT1-statsFile-m"$m" \
--statsFileBgenSnps="$outdir"/"$j"/BOLT1-statsFile-BgenSnps-m"$m"

done
