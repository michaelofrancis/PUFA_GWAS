library(plyr)
library(tidyverse)


phenotypes<-c("w3FA_NMR_resinv", "w3FA_NMR_TFAP_resinv", "w6FA_NMR_resinv",
                "w6FA_NMR_TFAP_resinv", "w6_w3_ratio_NMR_resinv",
                "DHA_NMR_resinv", "DHA_NMR_TFAP_resinv", "LA_NMR_resinv", 
                "LA_NMR_TFAP_resinv", "PUFA_NMR_resinv", "PUFA_NMR_TFAP_resinv", 
                "MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv", "PUFA_MUFA_ratio_NMR_resinv")

phenotypesb<-str_split(phenotypes, "_resinv", simplify=T)[,1]

phenotypesb<-c("Omega-3 Fatty Acids", "Omega-3 Fatty Acids to Total Fatty Acids percentage",
		"Omega-6 Fatty Acids", "Omega-6 Fatty Acids to Total Fatty Acids percentage",
		"Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio", 
		"Docosahexaenoic Acid", "Docosahexaenoic Acid to Total Fatty Acids percentage",
		"Linoleic Acid", "Linoleic Acid to Total Fatty Acids percentage",
		"Polyunsaturated Fatty Acids","Polyunsaturated Fatty Acids to Total Fatty Acids percentage",
		"Monounsaturated Fatty Acids", "Monounsaturated Fatty Acids to Total Fatty Acids percentage",
		"Polyunsaturated Fatty Acids to Monounsaturated Fatty Acids ratio")

phelabels<-c("Omega-3", "Omega-3_TFAP", "Omega-6", "Omega-6_TFAP", 
	"Omega-6_Omega-3_ratio", "DHA","DHA_TFAP", "LA", "LA_TFAP", 
	"PUFAs", "PUFAs_TFAP", "MUFAs", "MUFAs_TFAP", "PUFA_MUFA_ratio")

dir1<-"/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/MLMA/MLMA2-03092022/AFR"
dir2<-"/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/MLMA/MLMA2-03092022/CSA"
dir3<-"/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/MLMA/MLMA2-03092022/EAS"

outdir="/scratch/mf91122/PUFA-GWAS/sumstats-final/multi-ancestry"

for (p in 1:length(phenotypes)){

	AFR<-as_tibble(read_delim(
	        paste(dir1, "/", phenotypes[p], "/allchrs.mlma", sep="")
                                        ))
	CSA<-as_tibble(read_delim(
        paste(dir2, "/", phenotypes[p], "/allchrs.mlma", sep="")
                                        ))

	EAS<-as_tibble(read_delim(
        paste(dir3, "/", phenotypes[p], "/allchrs.mlma", sep="")
                                        ))

colnames(AFR)[4:9]<-paste(colnames(AFR)[4:9], "_AFR", sep="")
colnames(CSA)[4:9]<-paste(colnames(CSA)[4:9], "_CSA", sep="")
colnames(EAS)[4:9]<-paste(colnames(EAS)[4:9], "_EAS", sep="")


join1<-full_join(AFR, CSA, by=c("Chr", "SNP", "bp"))
join2<-full_join(join1, EAS, by=c("Chr", "SNP", "bp"))

write.csv(join2, paste(outdir, "/UKB-multi-ancestry-", phelabels[p],".csv", sep=""),
	quote=F, row.names=F)

}
