#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=cojo-slct-MA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=167:00:00
#SBATCH --mem=50000
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --array=1-22

#i for chromosome
i=$SLURM_ARRAY_TASK_ID

cd /work/kylab/mike/PUFA-GWAS/GCTA-COJO2/2.cojo-slct

gctadir=("/home/mf91122/GCTA/gcta_1.93.2beta")
sumdir=("/scratch/mf91122/PUFA-GWAS/sumstats-final/UKBEURKETMET/maCOJO")
snplist=("/scratch/mf91122/PUFA-GWAS/METAL/M2-N/GCTA-COJO/snplist")

bfiledir=("/scratch/mf91122/PUFA-GWAS/GCTA-COJO2/1.Bfile20k/bfile3")

outdir=("/scratch/mf91122/PUFA-GWAS/GCTA-COJO2/2.cojo-slct-meta-analysis.fixse")


phenotypes=("FAw3" "FAw6" "DHA" "LA" "MUFA")


mkdir -p "$outdir"
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
--cojo-file "$sumdir"/"$j".UKBEURKETMET.b.ma \
--cojo-p 2.439e-8 \
--cojo-slct \
--thread-num 20 \
--out "$outdir"/"$j"_chr"$i"

done
