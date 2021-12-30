#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=GCTA-fastGWA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=167:00:00
#SBATCH --mem=200000
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --array=1-2

i=$SLURM_ARRAY_TASK_ID

#fastGWA (https://cnsgenomics.com/software/gcta/#fastGWA) is a new method that can 
#quickly finish GWAS using mixed linear regression model for large cohort, while keeping 
#low false positive rate, controlling for population stratification by PC and relatedness 
#by a sparse genetic relationship matrix.

cd /work/kylab/mike/PUFA-GWAS/UKB-multi-ancestry/fastGWA

gctadir=("/home/mf91122/GCTA/gcta_1.93.2beta")
genoindir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/MLMA/genotypeQC")
phenodir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/pheno/phen")
grmdir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/GCTA-GRM/sparseGRM")
outdir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/fastGWA")

#pop=("AFR" "CSA" "EAS")
pop=("AFR")
#phenotypes=("w3FA_NMR_resinv" "w3FA_NMR_TFAP_resinv" "w6FA_NMR_resinv" 
#                "w6FA_NMR_TFAP_resinv" "w6_w3_ratio_NMR_resinv" 
#                "DHA_NMR_resinv" "DHA_NMR_TFAP_resinv" "LA_NMR_resinv" 
#                "LA_NMR_TFAP_resinv" "PUFA_NMR_resinv" "PUFA_NMR_TFAP_resinv" 
#                "MUFA_NMR_resinv" "MUFA_NMR_TFAP_resinv" "PUFA_MUFA_ratio_NMR_resinv")
phenotypes=("w3FA_NMR_resinv")

for p in ${pop[@]} 
        do
        mkdir -p "$outdir"/"$p"

        for j in ${phenotypes[@]}
                do

$gctadir/gcta64 \
--mbfile $genoindir/$p/merge-list.txt \
--grm-sparse $grmdir/$p/spgrm \
--fastGWA-mlm \
--geno 0.05 \
--pheno $phenodir/PUFA_GWAS_pheno_M"$i"."$p".txt.pc.txt."$j".phen \
--qcovar $phenodir/PUFA_GWAS_pheno_M"$i"."$p".txt.pc.txt.qcovar \
--out $outdir/"$p"/"$j"-M"$i"-fastGWA

done
done
