suppressMessages(library(qqman))
suppressMessages(library(tidyverse))

phenotypes<-c("w3FA_NMR_resinv", "w3FA_NMR_TFAP_resinv", "w6FA_NMR_resinv", 
		"w6FA_NMR_TFAP_resinv", "w6_w3_ratio_NMR_resinv", 
		"DHA_NMR_resinv", "DHA_NMR_TFAP_resinv", "LA_NMR_resinv", 
		"LA_NMR_TFAP_resinv", "PUFA_NMR_resinv", "PUFA_NMR_TFAP_resinv", 
		"MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv", "PUFA_MUFA_ratio_NMR_resinv")

m<-c(1,2)

for (m in 1:length(m)){ #loop over models

dir<-paste("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/M", m, "/", sep="")


for (i in 1:length(phenotypes)){ #loop over phenotypes


file<-read.table(paste(dir, phenotypes[i], "/BOLT1-statsFile-BgenSnps-m", m, sep=""), 
		header=TRUE, stringsAsFactors=FALSE)
file<-as_tibble(file)

file2<-file%>%select(CHR, SNP, BP, P_BOLT_LMM)
colnames(file2)[4]<-"P"
#file3<-file%>%select(CHR, SNP, BP, P_BOLT_LMM_INF)
#colnames(file3)[4]<-"P"

suppressWarnings(dir.create(paste(dir, phenotypes[i], "/plot", sep="")))

plotoutputpath=paste(dir, phenotypes[i], "/plot/ManhattanBOLT-10242021-m", m, sep="")


#Make Manhattan plot P_BOLT_LMM
png(filename=plotoutputpath, type="cairo")

manhattan(file2,
	main = paste(phenotypes[i], " BOLT-LMM UKB M", m, sep=""),
	ylim = c(0, 30), cex = 0.6,
	annotatePval = FALSE, annotateTop = FALSE,
	cex.axis = 0.9, col = c('red', 'black'),
	suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
	chrlabs = (1:22)
	)#end manhattan

dev.off() 


qqoutputpath=paste(dir, phenotypes[i], "/plot/QQBOLT-10242021-m", m, sep="")

png(filename=qqoutputpath, type="cairo")

qq(file2$P,
	main= paste(phenotypes[i], " BOLT-LMM UKB M", m, sep="")
	)

dev.off()

}#end loop over models
} #end loop over phenotypes
