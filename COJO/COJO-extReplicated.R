library(plyr)
library(tidyverse)

pheno<-c("FAw3", "FAw6", "DHA", "LA", "MUFA")
pheno2<-c("w3", "w6", "DHA", "LA","MUFA")
ph<-c("Omega-3 Fatty Acids",
      "Omega-6 Fatty Acids",
      "Docosahexaenoic Acid",
      "Linoleic Acid",
      "Monounsaturated Fatty Acids")

jma<-as_tibble(read.csv(
        "/work/kylab/mike/PUFA-GWAS/GCTA-COJO2/2.cojo-slct/table-process/meta-analysis.combined.jma.cojo.csv"))

jma$Phenotype<-mapvalues(jma$Phenotype, from =ph, to=pheno)
jma$replicated<-0

#i=1

for (i in 1:5){

join3<-as_tibble(read.table(
        paste("/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication/munge/UKBEURKETMET.mungecombined.metalinput/",
                pheno[i], ".UKBEURKETMET.mungedcombined.metalinput.txt", sep=""),
        header=T))

snp<-jma$SNP[jma$Phenotype==pheno[i]]
extsnp<-join3[join3$SNP %in% snp,] %>%filter(P_KET<0.05 | P_MET <0.05)%>%select(SNP)
extsnp<-unlist(extsnp)

jma$replicated[jma$Phenotype==pheno[i] & jma$SNP %in% extsnp]<-1

}

write.csv(jma,
"/work/kylab/mike/PUFA-GWAS/GCTA-COJO2/2.cojo-slct/table-process/meta-analysis.combined.jma.cojo.rep.csv",
quote=F, row.names=F)
