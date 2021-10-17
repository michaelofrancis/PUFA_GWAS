#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=Prune
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=180000
#SBATCH --output=Prune.%j.out
#SBATCH --error=Prune.%j.err

#cancl----#SBATCH --array=1-22

i=$SLURM_ARRAY_TASK_ID

ml PLINK/2.00-alpha2.3-x86_64-20210920-dev
cd /work/kylab/mike/PUFA-GWAS/Prune

#---------
#Set which
#steps run
#---------
step1=false
step2=false
step3=false
step4=false
step5=false #this merging process screwed things up for the bfile in BOLT so recommend stopping at step4.
#Cancel array jobs below here-----
step6=true
step7=true
#---------

if [ $step1 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 1. Convert BGEN to PFILE=-=-=-=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 1-=-=-=-=-=-=-=-\n\n"

#Start with UKB source genotype bgen files
genoindir=("/scratch/mf91122/bgen_v1.2_UKBsource")
mfiscoredir=("/work/kylab/mike/UKB/quality-scores/mfi")
outdir=("/scratch/mf91122/PUFA-GWAS/Prune/1.pfilefromrawbgen")

mkdir -p $outdir

plink2 \
--bgen $genoindir/ukb_imp_chr"$i"_v3.bgen ref-first \
--sample $genoindir/ukb_imp_v3.sample \
--keep /scratch/mf91122/PUFA-GWAS/pheno/PUFA_GWAS_phenoQC_IDS.txt \
--geno 0.01 \
--hwe 1e-08 \
--maf 0.01 \
--make-pfile \
--out $outdir/chr"$i"

fi


if [ $step2 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 2. Extract INFO>0.8 SNPs=-=-=-=-=-=- 
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

echo "-=-=-=-=-=-=-=-STEP 2-=-=-=-=-=-=-=-\n\n"
genoindir=$outdir
outdir=("/scratch/mf91122/PUFA-GWAS/Prune/2.extractINFOsnps")

mkdir -p $outdir

plink2 \
--pfile $genoindir/chr"$i" \
--extract $mfiscoredir/ukb_mfi_keepsnps_chr"$i"_0.8.txt \
--make-pgen \
--out $outdir/chr"$i"

fi


if [ $step3 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 3. Hard call and LD prune-=-=-=-=-=- 
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 3-=-=-=-=-=-=-=-\n\n"

genoindir=$outdir
outdir=("/scratch/mf91122/PUFA-GWAS/Prune/3.Hardcall-PruneSNPs")

mkdir -p $outdir

plink2 \
--pfile $genoindir/chr"$i" \
--exclude range /work/kylab/mike/CCC/C.QC/Prune1excludeset.txt \
--hard-call-threshold 0.1 \
--set-missing-var-ids @:#:\$r:\$a \
--rm-dup force-first \
--new-id-max-allele-len 414 \
--indep-pairwise 50 5 0.2 \
--out $outdir/chr"$i"

#Post-process
cat $outdir/chr*.prune.in > $outdir/Prune.in.txt

fi

if [ $step4 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 4. Write prune lists to bedfile-=-=- 
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 4-=-=-=-=-=-=-=-\n\n"

#genoindir same as step 2 so don't re set
prunelistdir=$outdir
outdir=("/scratch/mf91122/PUFA-GWAS/Prune/4.bedfinal")


mkdir -p $outdir

plink2 \
--pfile $genoindir/chr"$i" \
--exclude range /work/kylab/mike/CCC/C.QC/Prune1excludeset.txt \
--keep /scratch/mf91122/PUFA-GWAS/pheno/PUFA_GWAS_phenoQC_IDS.txt \
--extract $prunelistdir/Prune.in.txt \
--hard-call-threshold 0.1 \
--make-bed \
--out $outdir/chr"$i"

fi

if [ $step5 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 5. Make list to merge bedfiles=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 5-=-=-=-=-=-=-=-\n\n"

outdir=("/scratch/mf91122/PUFA-GWAS/Prune/4.bedfinal")
rm -f "$outdir"/merge.txt

chromosomes=({1..22})


for i in ${chromosomes[*]}
        do

        echo "$outdir"/chr"$i" >> "$outdir"/merge.txt

done


fi #end step 5


if [ $step6 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 6. Merge bedfiles-=-=-=-=-=-=-=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 6-=-=-=-=-=-=-=-\n\n"

outdir=("/scratch/mf91122/PUFA-GWAS/Prune/4.bedfinal")

plink2 \
--pmerge-list "$outdir"/merge.txt bfile \
--make-pgen \
--merge-max-allele-ct 2 \
--out "$outdir"/merged

fi #end step 6

if [ $step7 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 7. Convert merged pfile to bed=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 7-=-=-=-=-=-=-=-\n\n"

outdir=("/scratch/mf91122/PUFA-GWAS/Prune/4.bedfinal")

plink2 \
--pfile "$outdir"/merged \
--make-bed \
--geno 0.01 \
--hwe 1e-08 \
--maf 0.01 \
--mind 0.05 \
--out "$outdir"/mergedbed_mind0.05

fi #end step 7
