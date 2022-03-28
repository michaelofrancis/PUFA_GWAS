library(plyr)
library(tidyverse)

pheno<-c("w3FA", "w6FA", "DHA", "LA", "MUFA")
pheno2<-c("W3", "w6_METAL", "DHA", "LA_METAL","MUFA_METAL")

phenotypes<-c("w3FA_NMR_resinv", "w3FA_NMR_TFAP_resinv", "w6FA_NMR_resinv",
                "w6FA_NMR_TFAP_resinv", "w6_w3_ratio_NMR_resinv",
                "DHA_NMR_resinv", "DHA_NMR_TFAP_resinv", "LA_NMR_resinv",
                "LA_NMR_TFAP_resinv", "PUFA_NMR_resinv", "PUFA_NMR_TFAP_resinv",
                "MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv", "PUFA_MUFA_ratio_NMR_resinv")



ph<-c("Omega-3 Fatty Acids",
"Omega-3 Fatty Acids to Total Fatty Acids percentage", 
"Omega-6 Fatty Acids",
"Omega-6 Fatty Acids to Total Fatty Acids percentage",
"Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio",
"Docosahexaenoic Acid", "Docosahexaenoic Acid to Total Fatty Acids percentage",
"Linoleic Acid","Linoleic Acid to Total Fatty Acids percentage",
"Polyunsaturated Fatty Acids", 
"Polyunsaturated Fatty Acids to Total Fatty Acids percentage",
"Monounsaturated Fatty Acids","Monounsaturated Fatty Acids to Total Fatty Acids percentage",
"Polyunsaturated Fatty Acids to Monounsaturated Fatty Acids ratio")


m=2

dir<-paste("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/M", m, "/", sep="")

#Read in BOLT-LMM files
file<-list()
for (i in 1:length(phenotypes)){ 
	file[[i]]<-read.table(paste(dir, phenotypes[i], "/BOLT1-statsFile-BgenSnps-m", m, sep=""),
                header=TRUE, stringsAsFactors=FALSE)
	file[[i]]<-as_tibble(file[[i]])

	file[[i]]$Phenotype<-phenotypes[i]
	file[[i]]<-file[[i]]%>%select(Phenotype, everything())
}

#BOLTsum<-do.call(rbind, file) #too much memory

FUMAdisc<-as_tibble(read_csv("/work/kylab/mike/PUFA-GWAS/FUMA-final/FUMAdisc_S6_03272022.csv"))

FUMAdisc$Phenotype<-mapvalues(FUMAdisc$Phenotype, from=ph, to=phenotypes)



file2<-list()
for (i in 1:length(phenotypes)){
file2[[i]]<-file[[i]][
	file[[i]]$BP %in% FUMAdisc$pos,
	]
}


BOLTsum<-do.call(rbind, file2)

FUMAdisc$PheChrPos<-paste(FUMAdisc$Phenotype, FUMAdisc$CHR, FUMAdisc$POS, sep=":")
BOLTsum$PheChrPos<-paste(BOLTsum$Phenotype, BOLTsum$CHR, BOLTsum$BP, sep=":")


join<-left_join(FUMAdisc, BOLTsum)
join$Phenotype<-mapvalues(join$Phenotype, from=phenotypes, to=ph)


write.csv(join, "/work/kylab/mike/PUFA-GWAS/FUMA-final/table-process/alleles.discovery.csv",
	row.names=F, quote=F)
