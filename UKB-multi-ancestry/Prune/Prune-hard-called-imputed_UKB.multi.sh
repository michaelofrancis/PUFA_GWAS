#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=Prune_multi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=180000
#SBATCH --output=Prune_multi.%j.out
#SBATCH --error=Prune_multi.%j.err
#SBATCH --array=1-22

i=$SLURM_ARRAY_TASK_ID

ml PLINK/2.00-alpha2.3-x86_64-20210920-dev
cd /work/kylab/mike/PUFA-GWAS/UKB-multi-ancestry/Prune

#---------
#Set which
#steps run
#---------
step1=true
step2=true
step3=true
step4=true
step5=true #this merging process screwed things up for the bfile in BOLT so recommend stopping at step4 for BOLT but need them merged for PCA!!
#Cancel array jobs below here-----
step6=false
step7=false
#step6=true
#step7=true
#---------

popin=("AFR" "CSA" "EAS")

FOR LOOP ANCESTRIES
for p in ${popin[@]} 
	do



if [ $step1 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 1. Convert BGEN to PFILE=-=-=-=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 1-=-=-=-=-=-=-=-\n\n"

#Start with UKB source genotype bgen files
genoindir=("/scratch/mf91122/bgen_v1.2_UKBsource")
mfiscoredir=("/work/kylab/mike/UKB/quality-scores/mfi")
outdir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/Prune/1.pfilefromrawbgen")

mkdir -p $outdir/"$p"

plink2 \
--bgen $genoindir/ukb_imp_chr"$i"_v3.bgen ref-first \
--sample $genoindir/ukb_imp_v3.sample \
--keep /scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/pheno/PUFA_GWAS_phenoQC_IDS_M1."$p".txt \
--geno 0.01 \
--hwe 1e-08 \
--maf 0.01 \
--make-pfile \
--out $outdir/"$p"/chr"$i"

fi


if [ $step2 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 2. Extract INFO>0.8 SNPs=-=-=-=-=-=- 
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

echo "-=-=-=-=-=-=-=-STEP 2-=-=-=-=-=-=-=-\n\n"
genoindir=$outdir
outdir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/Prune/2.extractINFOsnps")

mkdir -p $outdir/"$p"

plink2 \
--pfile $genoindir/"$p"/chr"$i" \
--extract $mfiscoredir/ukb_mfi_keepsnps_chr"$i"_0.8.txt \
--make-pgen \
--out $outdir/"$p"/chr"$i"

fi


if [ $step3 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 3. Hard call and LD prune-=-=-=-=-=- 
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 3-=-=-=-=-=-=-=-\n\n"

genoindir=$outdir
outdir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/Prune/3.Hardcall-PruneSNPs")

mkdir -p $outdir/"$p"

plink2 \
--pfile $genoindir/"$p"/chr"$i" \
--exclude range /work/kylab/mike/CCC/C.QC/Prune1excludeset.txt \
--hard-call-threshold 0.1 \
--set-missing-var-ids @:#:\$r:\$a \
--rm-dup force-first \
--new-id-max-allele-len 414 \
--indep-pairwise 50 5 0.2 \
--out $outdir/"$p"/chr"$i"

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
outdir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/Prune/4.bedfinal")


mkdir -p $outdir/"$p"

plink2 \
--pfile $genoindir/"$p"/chr"$i" \
--exclude range /work/kylab/mike/CCC/C.QC/Prune1excludeset.txt \
--keep /scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/pheno/PUFA_GWAS_phenoQC_IDS_M1."$p".txt \
--extract $prunelistdir/"$p"/Prune.in.txt \
--hard-call-threshold 0.1 \
--make-bed \
--out $outdir/"$p"/chr"$i"

fi

if [ $step5 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 5. Make list to merge bedfiles=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 5-=-=-=-=-=-=-=-\n\n"

outdir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/Prune/4.bedfinal")
rm -f "$outdir"/"$p"/merge.txt

chromosomes=({1..22})


for i in ${chromosomes[*]}
        do

        echo "$outdir"/"$p"/chr"$i" >> "$outdir"/"$p"/merge.txt

done


fi #end step 5


if [ $step6 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 6. Merge bedfiles-=-=-=-=-=-=-=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 6-=-=-=-=-=-=-=-\n\n"

outdir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/Prune/4.bedfinal")

plink2 \
--pmerge-list "$outdir"/"$p"/merge.txt bfile \
--make-pgen \
--merge-max-allele-ct 2 \
--out "$outdir"/"$p"/merged

fi #end step 6

if [ $step7 = true ]; then
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 7. Convert merged pfile to bed=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 7-=-=-=-=-=-=-=-\n\n"

outdir=("/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/Prune/4.bedfinal")

plink2 \
--pfile "$outdir"/"$p"/merged \
--make-bed \
--geno 0.01 \
--hwe 1e-08 \
--maf 0.01 \
--mind 0.05 \
--out "$outdir"/"$p"/mergedbed_mind0.05

fi #end step 7


done #end for loop ancestry
