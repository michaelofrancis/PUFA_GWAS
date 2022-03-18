suppressMessages(library(qqman))
suppressMessages(library(tidyverse))

phenotypes<-c("w3FA_NMR_resinv", "w3FA_NMR_TFAP_resinv", "w6FA_NMR_resinv", 
		"w6FA_NMR_TFAP_resinv", "w6_w3_ratio_NMR_resinv", 
		"DHA_NMR_resinv", "DHA_NMR_TFAP_resinv", "LA_NMR_resinv", 
		"LA_NMR_TFAP_resinv", "PUFA_NMR_resinv", "PUFA_NMR_TFAP_resinv", 
		"MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv", "PUFA_MUFA_ratio_NMR_resinv")

phenotypenames<-c("Omega-3", "Omega-3 (TFAP)", "Omega-6" , "Omega-6 (TFAP)", 
		"Omega-6:Omega-3 ratio", "DHA", "DHA (TFAP)", "LA", "LA (TFAP)",
		"PUFAs", "PUFAs (TFAP)", "MUFAs", "MUFAs (TFAP)", "PUFAs:MUFAs ratio")

m<-c(1,2)
#m=2

mancolors<-list(c("#96D4C1", "#699487"),
		c("#0CEBA8", "#066C4E"))


for (m in 1:length(m)){ #loop over models

dir<-paste("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/M", m, "/", sep="")


for (i in 1:length(phenotypes)){ #loop over phenotypes
#m=1;i=1


file<-read.table(paste(dir, phenotypes[i], "/BOLT1-statsFile-BgenSnps-m", m, sep=""), 
		header=TRUE, stringsAsFactors=FALSE)
file<-as_tibble(file)

file2<-file%>%select(CHR, SNP, BP, P_BOLT_LMM)
colnames(file2)[4]<-"P"
#file3<-file%>%select(CHR, SNP, BP, P_BOLT_LMM_INF)
#colnames(file3)[4]<-"P"

suppressWarnings(dir.create(paste("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/plot/M", m, sep="")))

plotoutputpath=paste("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/plot/M", m, "/",
	phenotypes[i], "-Manhattan", sep="")

#Make Manhattan plot P_BOLT_LMM
png(filename=plotoutputpath, type="cairo", width=1000,height=500)

print(
manhattan(file2,
	main = paste("UKB-EUR discovery ", "M", m, " ", phenotypenames[i], sep=""),
	ylim = c(0, 300), cex = 0.6,
	annotatePval = FALSE, annotateTop = FALSE,
	cex.axis = 0.9, 
	col = mancolors[[m]],
	suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
	chrlabs = as.character(1:22)
	)#end manhattan
)

dev.off() 


qqoutputpath=paste("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/plot/M", m, "/", 
        phenotypes[i], "-QQ", sep="")

png(filename=qqoutputpath, type="cairo", width=500,height=500)

qq(file2$P,
	main = paste("UKB-EUR discovery ", "M", m, " ", phenotypenames[i], sep="")
	)

dev.off()

}#end loop over models
} #end loop over phenotypes
