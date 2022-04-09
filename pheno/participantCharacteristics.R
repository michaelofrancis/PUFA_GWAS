suppressMessages(library(tidyverse))
library(psych)

dir<-"/scratch/mf91122/PUFA-GWAS/pheno"

file1<-"PUFA_GWAS_pheno_M1.txt"
file2<-"PUFA_GWAS_pheno_M2.txt"

tab1<-read.table(paste(dir, file1, sep="/"), header=T, stringsAsFactors=F)
tab1<-as_tibble(tab1)
tab2<-read.table(paste(dir, file2, sep="/"), header=T, stringsAsFactors=F)
tab2<-as_tibble(tab2)

options(digits=4)
as.data.frame(tab1[c(2,3,5,16:29)]%>%describe())

tab2<-tab2[c(2:6,16:29,52)]
tab2<-tab2[complete.cases(tab2),]
as.data.frame(describe(tab2))
