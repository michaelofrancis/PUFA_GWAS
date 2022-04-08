library(tidyverse)

phe<-as_tibble(read.table("/Users/mike/Documents/Research/PUFA-GWAS/pheno/PUFA_GWAS_pheno_M2.txt",header=T))

pheno<-c("w3FA_NMR","w3FA_NMR_TFAP",
         "w6FA_NMR", "w6FA_NMR_TFAP",
         "w6_w3_ratio_NMR",
         "DHA_NMR","DHA_NMR_TFAP",
         "LA_NMR","LA_NMR_TFAP",
         "PUFA_NMR","PUFA_NMR_TFAP",
         "MUFA_NMR", "MUFA_NMR_TFAP",
         "PUFA_MUFA_ratio_NMR")

pheno<-paste(pheno, "_resinv", sep="")

phe2<-phe%>%select(pheno)
phe2<-phe2[complete.cases(phe2),]

A<-cov(phe2)
ev <- eigen(A)
values <- ev$values
values

top<-(sum(values))^2
bottom<-sum(values^2)

top/bottom
#[1] 2.978671



# meta ----------------------------------------------------------


phe3<-phe2%>%select(w3FA_NMR_resinv, w6FA_NMR_resinv, DHA_NMR_resinv, LA_NMR_resinv, MUFA_NMR_resinv)
A<-cov(phe3)
ev <- eigen(A)
values <- ev$values
values

top<-(sum(values))^2
bottom<-sum(values^2)

top/bottom
