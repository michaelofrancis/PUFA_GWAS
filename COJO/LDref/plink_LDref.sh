#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=bgenQC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=180000
#SBATCH --output=bgenQC.%j.out
#SBATCH --error=bgenQC.%j.err
#SBATCH --array=1-22

#Make UKB LD reference bfile from 20k random for COJO

i=$SLURM_ARRAY_TASK_ID

cd /work/kylab/mike/PUFA-GWAS/GCTA-COJO2/1.Bfile20k
ml PLINK/2.00-alpha2.3-x86_64-20210920-dev

genoindir=("/scratch/mf91122/bgen_v1.2_UKBsource")
prepdir=("/scratch/mf91122/PUFA-GWAS/GCTA-COJO2/1.Bfile20k")
outdir=("/scratch/mf91122/PUFA-GWAS/GCTA-COJO2/1.Bfile20k/bfile2")

mkdir -p $outdir

#04-07-2022 took out mind, geno, hwe filtering as I was losing snps that were significant
#In GWAS and also the summary stats in the COJO step were already filtered on these 
#parameters before running these.

plink2 \
--bgen $genoindir/ukb_imp_chr"$i"_v3.bgen ref-first \
--sample $genoindir/ukb_imp_v3.sample \
--extract "$prepdir"/keepsnps.txt \
--keep "$prepdir"/LDref_20kQCunrelated.fam \
--make-bed \
--out $outdir/chr$i
