phenotypes=("w3FA_NMR_resinv" "w3FA_NMR_TFAP_resinv" "w6FA_NMR_resinv" 
		"w6FA_NMR_TFAP_resinv" "w6_w3_ratio_NMR_resinv" 
		"DHA_NMR_resinv" "DHA_NMR_TFAP_resinv" "LA_NMR_resinv" 
		"LA_NMR_TFAP_resinv" "PUFA_NMR_resinv" "PUFA_NMR_TFAP_resinv" 
		"MUFA_NMR_resinv" "MUFA_NMR_TFAP_resinv" "PUFA_MUFA_ratio_NMR_resinv")


dir="/scratch/mf91122/PUFA-GWAS/GCTA-COJO2/2.cojo-slct-discovery3"
outdir="/scratch/mf91122/PUFA-GWAS/GCTA-COJO2/2.cojo-slct-discovery3/jma"
mkdir $outdir


for j in ${phenotypes[@]} 
	do

cat "$dir"/"$j"_chr*.jma.cojo > "$outdir"/"$j".jma.cojo

done
