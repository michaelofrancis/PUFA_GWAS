library(tidyverse)
library(MungeSumstats)

dir<-"/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication"

pheno<-c("FAw3", "FAw6", "DHA", "LA", "MUFA")

pheno2<-c("w3FA", "w6FA", "DHA", "LA", "MUFA")

#i=1
for (i in 1:length(pheno)){

join3<-as_tibble(read.table(
        paste("/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication/munge/UKBEURKETMET.mungecombined.metalinput/", pheno[i], ".UKBEURKETMET.mungedcombined.metalinput.txt", sep=""),
        header=T))


METAL<-as_tibble(read.table(
        paste("/scratch/mf91122/PUFA-GWAS/METAL/mungecombined/META_IVW_", pheno[i],"1.txt",sep=""),
        header=T))

#METAL alleles are default lowercase, set to upper
METAL$Allele1<-toupper(METAL$Allele1)
METAL$Allele2<-toupper(METAL$Allele2)

#Tag metal columns before join
colnames(METAL)[2:11]<-paste(colnames(METAL)[2:11], "_METAL", sep="")


joinM<-full_join(join3, METAL, by=c("SNP"="MarkerName"))

#As of here: A2 is effect allele for join3, Allele1 is effect allele for METAL. Let's label them accordingly
colnames(joinM)[4:5]<-c("NEA", "EA")
colnames(joinM)[24:25]<-c("EA_METAL", "NEA_METAL")

#Where EA = EA_METAL and NEA = NEA_METAL, METAL values do not need to be changed, for example:
#as.data.frame(head(joinM%>%filter(EA==EA_METAL, NEA==NEA_METAL)%>%select(BETA_UKB, BETA_KET, BETA_MET, Direction_METAL, FRQ_UKB, FRQ_KET, FRQ_MET, Freq1_METAL), 20))

#Where EA = NEA_METAL and NEA = EA_METAL:
#BETA = - BETA
#FRQ = 1-FRQ

joinM$Freq1_METAL[joinM$NEA==joinM$EA_METAL & joinM$EA==joinM$NEA_METAL]= 1-joinM$Freq1_METAL[joinM$NEA==joinM$EA_METAL & joinM$EA==joinM$NEA_METAL]
joinM$Effect_METAL[joinM$NEA==joinM$EA_METAL & joinM$EA==joinM$NEA_METAL]= -joinM$Effect_METAL[joinM$NEA==joinM$EA_METAL & joinM$EA==joinM$NEA_METAL]

joinM<-joinM%>%select(-EA_METAL, -NEA_METAL, -Direction_METAL)
        
joinM<-joinM%>%rowwise()%>%mutate(Ntotal=sum(N_UKB, N_KET, NSTUDY_MET, na.rm=TRUE))

write.table(joinM, paste("/scratch/mf91122/PUFA-GWAS/sumstats-final/UKBEURKETMET/", pheno[i], ".UKBEURKETMET.txt", sep=""),
                quote=F, row.names=F, fileEncoding = "ascii")

}
