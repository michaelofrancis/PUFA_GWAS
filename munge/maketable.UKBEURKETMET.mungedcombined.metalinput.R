library(tidyverse)
library(MungeSumstats)

dir<-"/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication"

pheno<-c("FAw3", "FAw6", "DHA", "LA", "MUFA")

pheno2<-c("w3FA", "w6FA", "DHA", "LA", "MUFA")

#i=1
for (i in 1:length(pheno)){

#Load BOLT UKB-EUR files
UKBmunge<-as_tibble(read.table(
	paste(dir, "/munge/UKBN/", pheno2[i], "_NMR_resinv.M2.txt",sep=""),
        header=T))

#LOAD METSIM
	metfile<-paste(dir, "/METSIM/", "Locke_", pheno[i], ".txt", sep="")
	MET<-as_tibble(read.table(metfile, header=T))
	MET2<-MET%>%select(CHROM, BEG, NS, MAF, PVALUE, BETA, SEBETA, SNP, NonRefAllele, RefAllele)
	colnames(MET2)<-c("CHR","BP", "NSTUDY", "MAF", "PVALUE", "BETA", "SE", "SNP", "A2", "A1")
	METmunge<-as_tibble(read.table(paste(dir, "/munge/METSIM-EAF/", pheno[i],".txt", sep=""),
			header=T))

#LOAD KETTUNEN
ketfile<-paste(dir, "/Kettunen/", "MAGNETIC_", pheno[i], ".txt", sep="")
KET<-as_tibble(read.table(ketfile, header=T))
#Munge ket
KET2<-KET
colnames(KET2)<-c("CHR","BP","ID","EFFECT_ALLELE","NON_EFFECT_ALLELE","EFFECT_ALLELE_FREQ","BETA","SE","P","NSTUDY","N")

KETmunge<-as_tibble(read.table(paste(dir, "/munge/Kettunen/", pheno[i],".txt", sep=""), 
	header=T))

print(paste("nrow UKBmunge: ", nrow(UKBmunge)))
print(paste("nrow KETmunge: ", nrow(KETmunge)))
print(paste("nrow METmunge: ", nrow(METmunge)))


#Join UKB munge and KET munge
join<-full_join(UKBmunge, KETmunge, by=c("SNP", "CHR", "BP", "A1", "A2"))
colnames(join)<-c("SNP","CHR","BP","A1","A2","FRQ_UKB",  "INFO_UKB",   "BETA_UKB",
 		"SE_UKB",   "P_UKB",  "N_UKB","FRQ_KET","BETA_KET","SE_KET","P_KET","NSTUDY_KET", "N_KET") 

print(paste("nrow join: ", nrow(join)))

join2<-full_join(join, METmunge, by=c("SNP", "CHR", "BP", "A1", "A2"))
colnames(join2)<-c("SNP","CHR","BP","A1","A2","FRQ_UKB",  "INFO_UKB",   "BETA_UKB",
                "SE_UKB",   "P_UKB",  "N_UKB","FRQ_KET","BETA_KET","SE_KET","P_KET","NSTUDY_KET", "N_KET",
		"NSTUDY_MET", "FRQ_MET", "P_MET", "BETA_MET", "SE_MET", "ID_MET")


print(paste("nrow join2: ", nrow(join2)))

join3<-join2%>%filter(!duplicated(SNP))

print(paste("nrow join3: ", nrow(join3)))

write.table(join3, 
	paste("/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication/munge/",pheno[i], 
		".UKBEURKETMET.mungedcombined.metalinput.txt", sep=""),
	row.names=F, quote=F)

}
