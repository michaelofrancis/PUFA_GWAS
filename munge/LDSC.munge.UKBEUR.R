library(MungeSumstats)
suppressMessages(library(dplyr))

dir<-"/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication"

#Read UKB
pheno2<-c("w3FA_NMR_resinv","w3FA_NMR_TFAP_resinv",
                "w6FA_NMR_resinv", "w6FA_NMR_TFAP_resinv",
                "w6_w3_ratio_NMR_resinv",
                "DHA_NMR_resinv","DHA_NMR_TFAP_resinv",
                "LA_NMR_resinv","LA_NMR_TFAP_resinv",
                "PUFA_NMR_resinv","PUFA_NMR_TFAP_resinv",
                "MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv",
                "PUFA_MUFA_ratio_NMR_resinv")

UKBdir<-"/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/"


M=c(1,2)

for (i in 1:length(pheno2)){

	for (m in 1:length(M)){
	
	#i=1;m=1
	UKB<-as_tibble(read.table(paste(UKBdir,"M", m, 
		 "/", pheno2[i], "/BOLT1-statsFile-BgenSnps-m",m, sep=""), header=T))

	#Munge UKB
	UKB2<-UKB%>%select(-P_BOLT_LMM_INF, -GENPOS)
	colnames(UKB2)<-c("SNP", "CHR","BP", "EFFECT_ALLELE", "NON_EFFECT_ALLELE",
		"EFFECT_ALLELE_FRQ","INFO","BETA","SE", "P")
	reformattedUKB<-MungeSumstats::format_sumstats(path=UKB2, ref_genome="GRCh37", ldsc_format=T, 
		compute_n=102000, INFO_filter=0.3)

	refoUKB<-as_tibble(read.table(reformattedUKB, header=T))
	#Remove infinite Z values
	refoUKB<-refoUKB[!is.infinite(refoUKB$Z),]


	write.table(refoUKB, paste(dir, "/munge/UKBldsc/", pheno2[i],".M", m, ".txt", sep=""), quote=F, row.names=F)
	}#end model loop
} #end UKB loop
