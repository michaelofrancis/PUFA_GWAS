#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=cojo-slct-disc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=48:00:00
#SBATCH --mem=50000
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --array=1-22

#i for chromosome
i=$SLURM_ARRAY_TASK_ID


cd /work/kylab/mike/PUFA-GWAS/GCTA-COJO2/cojo-slct-discovery

gctadir=("/home/mf91122/GCTA/gcta_1.93.2beta")
sumdir=("/scratch/mf91122/PUFA-GWAS/GCTA-COJO2/cojo-slct-discovery/ma")
bfiledir=("/scratch/mf91122/PUFA-GWAS/GCTA-COJO2/1.Bfile20k/bfile3")
outdir=("/scratch/mf91122/PUFA-GWAS/GCTA-COJO2/2.cojo-slct-discovery3")

phenotypes=("w3FA_NMR_resinv" "w3FA_NMR_TFAP_resinv" "w6FA_NMR_resinv" 
		"w6FA_NMR_TFAP_resinv" "w6_w3_ratio_NMR_resinv" 
		"DHA_NMR_resinv" "DHA_NMR_TFAP_resinv" "LA_NMR_resinv" 
		"LA_NMR_TFAP_resinv" "PUFA_NMR_resinv" "PUFA_NMR_TFAP_resinv" 
		"MUFA_NMR_resinv" "MUFA_NMR_TFAP_resinv" "PUFA_MUFA_ratio_NMR_resinv")


###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###RUN GCTA-COJO-SLCT=-=-=-=-=-=-=-=
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


for j in ${phenotypes[@]}
           do

mkdir -p "$outdir"

$gctadir/gcta64 \
--bfile $bfiledir/chr"$i" \
--chr $i \
--maf 0.01 \
--cojo-file "$sumdir"/"$j".ma \
--cojo-slct \
--thread-num 16 \
--cojo-p 1.678e-8 \
--out "$outdir"/"$j"_chr"$i"

done
