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

i=$SLURM_ARRAY_TASK_ID

cd /work/kylab/mike/PUFA-GWAS/BOLT-LMM/genotypeQC
#ml PLINK/1.9b_5-x86_64
ml PLINK/2.00-alpha2.3-x86_64-20210920-dev

#---------
#Set which
#steps run
#---------
step1=false
step2=true
step3=false
#---------

if [ $step1 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 1. GET INFO >0.3-=-=-=-=-=-=-=-=-=-=
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 1-=-=-=-=-=-=-=-\n\n"

mfiscoredir=("/work/kylab/mike/UKB/quality-scores/mfi")

awk '{if ($8 >=0.3) print $2}' $mfiscoredir/ukb_mfi_chr"$i"_v3.txt > $mfiscoredir/ukb_mfi_keepsnps_chr"$i"_0.3.txt

fi

if [ $step2 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 2. GENOTYPE QC PLINK-=-=-=-=-=-=-=-=
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 2-=-=-=-=-=-=-=-\n\n"

genoindir=("/scratch/mf91122/bgen_v1.2_UKBsource")
mfiscoredir=("/work/kylab/mike/UKB/quality-scores/mfi")
outdir=("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/genotypeQC")
mkdir -p $outdir

plink2 \
--bgen $genoindir/ukb_imp_chr"$i"_v3.bgen ref-first \
--sample $genoindir/ukb_imp_v3.sample \
--extract $mfiscoredir/ukb_mfi_keepsnps_chr"$i"_0.3.txt \
--mind 0.05 \
--geno 0.05 \
--hwe 1e-08 \
--maf 0.001 \
--autosome \
--maj-ref \
--keep /scratch/mf91122/PUFA-GWAS/pheno/PUFA_GWAS_phenoQC_IDS.txt \
--export bgen-1.2 bits=8 \
--out "$outdir"/chr"$i"

fi #end step 2 if

if [ $step3 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 3. KEEP SELECTED PARTICIPANTS-=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

genoindir=("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/genotypeQC")

#Keep list from (Prune-hard-called-imputed.sh step 7)

plink2 \
--bgen $genoindir/chr"$i".bgen ref-first \
--sample $genoindir/chr"$i".sample \
--keep /scratch/mf91122/PUFA-GWAS/Prune/4.bedfinal/mergedbed_mind0.05.fam.keep \
--export bgen-1.2 bits=8 \
--out $genoindir/chr"$i"_b-keep

fi #end step 3 if
