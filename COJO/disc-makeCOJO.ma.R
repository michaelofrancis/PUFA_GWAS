library(tidyverse)
library(MungeSumstats)

phenotypes<-c("w3FA_NMR_resinv", "w3FA_NMR_TFAP_resinv", "w6FA_NMR_resinv",
                "w6FA_NMR_TFAP_resinv", "w6_w3_ratio_NMR_resinv",
                "DHA_NMR_resinv", "DHA_NMR_TFAP_resinv", "LA_NMR_resinv",
                "LA_NMR_TFAP_resinv", "PUFA_NMR_resinv", "PUFA_NMR_TFAP_resinv",
                "MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv", "PUFA_MUFA_ratio_NMR_resinv")

UKBdir<-"/scratch/mf91122/PUFA-GWAS/BOLT-LMM/results/"
m=1

for (p in 1:length(phenotypes)){

        UKB<-as_tibble(read.table(paste(UKBdir,"M", m,
                 "/", phenotypes[p], "/BOLT1-statsFile-BgenSnps-m",m, sep=""), header=T))


ma<-UKB%>%select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM)
colnames(ma)<-c("SNP", "A1", "A2", "freq", "b", "se", "p")

ma$N<-101729

write.table(ma,
        paste("/scratch/mf91122/PUFA-GWAS/GCTA-COJO2/cojo-slct-discovery/ma/",
                        phenotypes[p], ".ma", sep=""),
        quote=F, row.names=F)


}
