#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=ldsc-sumstatprep
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=180000
#SBATCH --output=ldsc-sumstatprep.%j.out
#SBATCH --error=ldsc-sumstatprep.%j.err


#Before running LDSC, prep summary stat files
#1. change headers
#2.convert P values that are too low to zero.

phenotypes=("w3FA_NMR_resinv" "w3FA_NMR_TFAP_resinv" "w6FA_NMR_resinv" 
                "w6FA_NMR_TFAP_resinv" "w6_w3_ratio_NMR_resinv" 
                "DHA_NMR_resinv" "DHA_NMR_TFAP_resinv" "LA_NMR_resinv" 
                "LA_NMR_TFAP_resinv" "PUFA_NMR_resinv" "PUFA_NMR_TFAP_resinv" 
                "MUFA_NMR_resinv" "MUFA_NMR_TFAP_resinv" "PUFA_MUFA_ratio_NMR_resinv")


sumdir=("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results")
outdir=("/scratch/mf91122/PUFA-GWAS/ldsc/out")

###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Pt 1 change headers with sed
#Use '1s/P_BOLT_LMM/P/' twice because there are two columns that start with that.
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

for j in ${phenotypes[@]} 
        do

for m in {1..2}
	do

sed -e '1s/ALLELE0/a2/' -e '1s/P_BOLT_LMM/P_INF/' -e '1s/P_BOLT_LMM/P_raw/' \
$sumdir/M"$m"/"$j"/BOLT1-statsFile-BgenSnps-m"$m" > $sumdir/M"$m"/"$j"/BOLT1-statsFile-BgenSnps-m"$m"-ldsccols

done
done


###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Pt 2 convert p values that are too low to zero.
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
for j in ${phenotypes[@]} 
        do

for m in {1..2}
        do


awk 'BEGIN{FS=OFS="\t"}{if(strtonum($12)<4.94066e-324) print $0,strtonum(0);else{print $0,$12}}' \
	$sumdir/M"$m"/"$j"/BOLT1-statsFile-BgenSnps-m"$m"-ldsccols > \
	$sumdir/M"$m"/"$j"/BOLT1-statsFile-BgenSnps-m"$m"-ldsccols.transfer

sed -e '1s/0/P/' $sumdir/M"$m"/"$j"/BOLT1-statsFile-BgenSnps-m"$m"-ldsccols.transfer > \
        $sumdir/M"$m"/"$j"/BOLT1-statsFile-BgenSnps-m"$m"-ldsccols.transfer2

done
done
