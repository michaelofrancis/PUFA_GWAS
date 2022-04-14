library(tidyverse)

dir<-"/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication"

pheno<-c("FAw3", "FAw6", "DHA", "LA", "MUFA")

#Replication table of UKB-EUR with external EUR studies

cols<-c(
    "nrow table",

    "Variants in UKB-EUR (munged)",

    "Significant variants UKB-EUR",

    "Variants in FinMetseq (munged)",

    "Overlapping variants UKB-EUR and FinMetSeq",

    "Associations at P<0.05 FinMetSeq (munged)",

    "Significant variants UKB-EUR replicated in FinMetSeq",

    "Variants in Ket (munged)",

    "Overlapping variants UKB-EUR and KET",

    "Associations at P<0.05 Kettunen et al. (munged)",

    "Significant variants UKB-EUR replicated in KET",

"Significant variants UKB-EUR replicated in KET or MET",

"Significant variants UKB-EUR replicated in KET AND MET")

tab<-as.data.frame(matrix(nrow=5,ncol=length(cols)))
colnames(tab)<-cols


for (i in 1:length(pheno)) {
join3<-as_tibble(read.table(paste("/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication/munge/UKBEURKETMET.mungecombined.metalinput/",
                pheno[i],".UKBEURKETMET.mungedcombined.metalinput.txt", sep=""), header=T))

tab[i,]<-c(
    nrow(join3),
    
    nrow(join3[!is.na(join3$P_UKB),]),
    
    nrow(join3%>%filter(P_UKB < 1.678e-08)),
    
    nrow(join3[!is.na(join3$P_MET),]),
    
    nrow(join3[!is.na(join3$P_UKB) & !is.na(join3$P_MET),]),
    
    nrow(join3%>%filter(P_MET <0.05)),

    nrow(join3%>%filter(P_UKB < 1.678e-08 &
                P_MET <0.05)),
    
    nrow(join3[!is.na(join3$P_KET),]),

    nrow(join3[!is.na(join3$P_UKB) & !is.na(join3$P_KET),]),
    
    nrow(join3%>%filter(P_KET <0.05)),

    nrow(join3%>%filter(P_UKB < 1.678e-08 &
                P_KET <0.05)),
    
    nrow(join3%>%filter(P_UKB < 1.678e-08 &(
                P_KET <0.05 |
                P_MET <0.05
                )
            )),
    
    nrow(join3%>%filter(P_UKB < 1.678e-08 &(
                P_KET <0.05 &
                P_MET <0.05
                )
            )))

}

write.csv(tab,"/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication/munge/UKBEURKETMET.mungecombined.metalinput/replication.csv")
