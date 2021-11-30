#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=ldsc1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=180000
#SBATCH --output=ldsc1.%j.out
#SBATCH --error=ldsc1.%j.err

cd /work/kylab/mike/PUFA-GWAS/ldsc

#Load conda env
#It was necessary to create the conda environment after 
#ml Anaconda 2 because the script is in python 2..?
ml Anaconda2/5.3.0
source activate /home/mf91122/.conda/envs/ldsc

filedir=("/scratch/mf91122/PUFA-GWAS/ldsc/tutorial/1kg_eur")
filedir2=("/scratch/mf91122/PUFA-GWAS/ldsc/tutorial")
ldscdir=("/home/mf91122/ldsc/ldsc")

phenotypes=("w3FA_NMR_resinv" "w3FA_NMR_TFAP_resinv" "w6FA_NMR_resinv" 
                "w6FA_NMR_TFAP_resinv" "w6_w3_ratio_NMR_resinv" 
                "DHA_NMR_resinv" "DHA_NMR_TFAP_resinv" "LA_NMR_resinv" 
                "LA_NMR_TFAP_resinv" "PUFA_NMR_resinv" "PUFA_NMR_TFAP_resinv" 
                "MUFA_NMR_resinv" "MUFA_NMR_TFAP_resinv" "PUFA_MUFA_ratio_NMR_resinv")


sumdir=("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results")
outdir=("/scratch/mf91122/PUFA-GWAS/ldsc/out")

for j in ${phenotypes[@]} 
        do

for m in {1..2}
        do



#https://github.com/bulik/ldsc/wiki
python $ldscdir/munge_sumstats.py \
	--sumstats $sumdir/M"$m"/"$j"/BOLT1-statsFile-BgenSnps-m"$m"-ldsccols.transfer2 \
	--N 118000 \
	--out $outdir/w3test_mungesumstats \
	--merge-alleles $filedir2/w_hm3.snplist

done
done
