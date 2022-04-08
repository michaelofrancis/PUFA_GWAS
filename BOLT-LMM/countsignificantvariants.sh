phenotypes=("w3FA_NMR_resinv" "w3FA_NMR_TFAP_resinv" "w6FA_NMR_resinv" 
		"w6FA_NMR_TFAP_resinv" "w6_w3_ratio_NMR_resinv" 
		"DHA_NMR_resinv" "DHA_NMR_TFAP_resinv" "LA_NMR_resinv" 
		"LA_NMR_TFAP_resinv" "PUFA_NMR_resinv" "PUFA_NMR_TFAP_resinv" 
		"MUFA_NMR_resinv" "MUFA_NMR_TFAP_resinv" "PUFA_MUFA_ratio_NMR_resinv")


M1dir=("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/M1")
M2dir=("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/M2")

for j in ${phenotypes[@]}
           do


	wc -l $M1dir/"$j"/BOLT1-statsFile-BgenSnps-m1
	awk '{if ($12 < 1.678e-08) print $0}' $M1dir/BOLT1-statsFile-BgenSnps-m1 | wc -l
	wc -l $M2dir/"$j"/BOLT1-statsFile-BgenSnps-m2
        awk '{if ($12 < 1.678e-08) print $0}' $M2dir/BOLT1-statsFile-BgenSnps-m2 | wc -l

done


for j in ${phenotypes[@]}
           do

	awk '{if ($12 < 1.678e-08) print $1}' $M1dir/"$j"/BOLT1-statsFile-BgenSnps-m1 >> $M1dir/allsnps.txt
	awk '{if ($12 < 1.678e-08) print $1}' $M2dir/"$j"/BOLT1-statsFile-BgenSnps-m2 >> $M2dir/allsnps.txt

done


#R---------------------
#snps1<-read.table("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/M1/allsnps.txt", header=F)
#length(unique(snps1$V1))
#[1] 29786

#snps2<-read.table("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/M2/allsnps.txt", header=F)
#length(unique(snps2$V1))
#[1] 39565
