#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=MLMA.phen2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=144:00:00
#SBATCH --mem=50000
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --array=1-22

i=$SLURM_ARRAY_TASK_ID

cd /work/kylab/mike/PUFA-GWAS/UKB-multi-ancestry/MLMA

#Set directories
gctadir=("/home/mf91122/GCTA/gcta_1.94.0beta")
genoindir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/MLMA/genotypeQC2")
phenodir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/pheno/phen2")
grmdir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/GCTA-GRM/GRM-Y")
outdir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/MLMA/MLMA2-03092022")

pop=("AFR" "CSA" "EAS")
#pop=("AFR")


phenotypes=("w3FA_NMR_resinv" "w3FA_NMR_TFAP_resinv" "w6FA_NMR_resinv" 
                "w6FA_NMR_TFAP_resinv" "w6_w3_ratio_NMR_resinv" 
                "DHA_NMR_resinv" "DHA_NMR_TFAP_resinv" "LA_NMR_resinv" 
                "LA_NMR_TFAP_resinv" "PUFA_NMR_resinv" "PUFA_NMR_TFAP_resinv" 
                "MUFA_NMR_resinv" "MUFA_NMR_TFAP_resinv" "PUFA_MUFA_ratio_NMR_resinv")

#phenotypes=("w3FA_NMR_resinv")

for p in ${pop[@]} 
	do

	for j in ${phenotypes[@]}
		do
		mkdir -p "$outdir"/"$p"/"$j"


$gctadir/gcta64 \
--mlma \
--bfile $genoindir/$p/chr"$i" \
--grm $grmdir/"$p"/full \
--pheno $phenodir/PUFA_GWAS_pheno_M2."$p".txt.pc.txt."$j".phen \
--qcovar $phenodir/PUFA_GWAS_pheno_M2."$p".txt.pc.txt.qcovar \
--thread-num 32 \
--out "$outdir"/"$p"/"$j"/chr"$i"-M2-mlm

done
done
