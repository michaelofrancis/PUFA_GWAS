library(MungeSumstats)
suppressMessages(library(dplyr))

dir<-"/scratch/mf91122/PUFA-GWAS/UKB-multi-ancestry/fastGWA"

pop<-c("AFR", "CSA", "EAS")

pheno2<-c("w3FA_NMR_resinv","w3FA_NMR_TFAP_resinv",
                "w6FA_NMR_resinv", "w6FA_NMR_TFAP_resinv",
                "w6_w3_ratio_NMR_resinv",
                "DHA_NMR_resinv","DHA_NMR_TFAP_resinv",
                "LA_NMR_resinv","LA_NMR_TFAP_resinv",
                "PUFA_NMR_resinv","PUFA_NMR_TFAP_resinv",
                "MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv",
                "PUFA_MUFA_ratio_NMR_resinv")

for (i in 1:length(pheno2)){
	for (p in 1:length(pop)){
		for (m in 1:2){

	#i=1;p=1
	file<-as_tibble(read.table(paste(dir, "/", pop[p], "/", pheno2[i], 
			"-M", m, ".fastGWA", sep=""), header=T))

	colnames(file)<-c("CHR", "SNP", "BP", "A1", "A2", "N", "FRQ", "BETA", "SE", "P")

        reformattedfile<-MungeSumstats::format_sumstats(path=file, ref_genome="GRCh37")

        refofile<-as_tibble(read.table(reformattedfile, header=T))

	write.table(refofile, paste(dir, "/munge/", pop[p], "/", pheno2[i], ".M",m, ".txt", sep=""), 
                quote=F, row.names=F)

		} #end models
	} #end pop
} #end pheno
