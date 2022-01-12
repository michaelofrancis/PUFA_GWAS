suppressMessages(library(tidyverse))

METKETMETALdir<-"/scratch/mf91122/PUFA-GWAS/METAL/MET_KET_only"
METdir="/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication/METSIM"
KETdir<-"/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication/Kettunen"
UKBdiscoveryEURdir<-"/work/kylab/mike/PUFA-GWAS/FUMA-UKBEURdiscovery"

pheno<-c("DHA", "LA", "MUFA", "FAw3", "FAw6")
pheno1<-c("DHA_NMR", "LA_NMR", "MUFA_NMR", "w3FA_NMR", "w6FA_NMR")
pheno2<-c("DHA.txt", "LA.txt", "MUFA.txt", "FAw3.txt", "FAw6.txt")

#Load FUMA UKB European discovery FUMA output
FUMAUKBEUR<-as_tibble(read.csv(paste(UKBdiscoveryEURdir, "/GenomicLoci12112021.csv" ,sep=""),
	header=T))

#Initialize results table
results<-data.frame()


#=-=-=-=-=-=-==-=-=-=-=-=-==-=-=-=-=-=
#FOR LOOP OVER PHENOTYPES=-=-=-=-=-=-=
#=-=-=-=-=-=-==-=-=-=-=-=-==-=-=-=-=-=

for (p in 1:length(pheno1)){
#p=1

#Initialize empty snptable
snptable<-data.frame()

#Get just the phenotype of interest in this FUMA table
FUMAp<-FUMAUKBEUR[FUMAUKBEUR$Phenotype==pheno1[p],]


#------------------------------------------------------------------
#----*Load METSIM, Kettunen, and meta-analyzed MET+KET datasets----
#------------------------------------------------------------------

#Load MET+KET metal output
MKM<-as_tibble(read.table(paste(METKETMETALdir, "/META_IVW_", pheno[p], ".STDERR1.txt", sep=""),
                header=T))
#print(paste("Total number of rows in", pheno2[p]," for KET+MET meta:", nrow(MKM)))
#Get significant vars
MKMsig<-MKM[MKM$P.value<0.05,]
#print(paste("P<0.05 rows in", pheno2[p]," for KET+MET meta:", nrow(MKMsig)))

#Load MET
MET<-as_tibble(read.table(paste(METdir, "/Locke_", pheno2[p], sep=""), header=T))
#print(paste("Total number of rows in", pheno2[p]," for METSIM:", nrow(MET)))
#Get significant vars
METsig<-MET[MET$PVALUE<0.05,]
#print(paste("P<0.05 rows in", pheno2[p]," for METSIM:", nrow(METsig)))

#Load KET
KET<-as_tibble(read.table(paste(KETdir, "/MAGNETIC_", pheno2[p], sep=""), header=T))
#print(paste("Total number of rows in", pheno2[p]," for Kettunen:", nrow(KET)))
KETsig<-KET[KET$p_value<0.05,]
#print(paste("P<0.05 rows in", pheno2[p]," for Kettunen:", nrow(KETsig)))

#------------------------------------------------------------------

#Process FUMA dataset to get list of SNPs into table

#Get primary rsID SNPs from FUMA into table
SNPs<-FUMAp$rsID
tab<-data.frame(pheno1[p], SNPs)
tab$Rank<-"Primary"
colnames(tab)<-c("Pheno", "SNP", "Rank")

#Get secondary SNPs from FUMA into table
SNP2<-FUMAp%>%select(IndSigSNPs)
SNP2<-as.character(SNP2)
#Deal with the ";" separated values
SNP2<-str_split(SNP2, pattern=";|,|\"")
SNP2<-unlist(SNP2)
SNP2<-SNP2[grep("rs", SNP2)]
tab2<-data.frame(pheno1[p], SNP2)
tab2$Rank<-"Secondary"
colnames(tab2)<-c("Pheno", "SNP", "Rank")

#Combine Primary and secondary SNPs into table
snptable<-as_tibble(rbind(tab, tab2))
#Initialize Replicated columns
snptable$ReplicatedMET<-0
snptable$ReplicatedKET<-0 
snptable$ReplicatedMETKET<-0 

#Identify significant snps from MET, KET, and MET+KET/METAL output in the snptable
snptable$ReplicatedMET[snptable$SNP %in% METsig$SNP]<-1
snptable$ReplicatedKET[snptable$SNP %in% KETsig$ID]<-1
snptable$ReplicatedMETKET[snptable$SNP %in% MKMsig$MarkerName]<-1

#May not need below lineV
#snpreplicated<-snptable$SNP[snptable$Replicated==1]

#Transfer that significant info to the appropriate line in the FUMA table

#Initialize FUMA replicated columns
FUMAp$ReplicatedMET<-0
FUMAp$ReplicatedKET<-0
FUMAp$ReplicatedMETKET<-0

#Identify rsID SNPs replicated in FUMA table
FUMAp$ReplicatedMET[FUMAp$rsID %in% snptable$SNP[snptable$ReplicatedMET==1]]<-1
FUMAp$ReplicatedKET[FUMAp$rsID %in% snptable$SNP[snptable$ReplicatedKET==1]]<-1
FUMAp$ReplicatedMETKET[FUMAp$rsID %in% snptable$SNP[snptable$ReplicatedMETKET==1]]<-1


#The part where you look in the ";" separated IndSigSNPs column to match the significant variants
for (i in 1:nrow(snptable[snptable$ReplicatedMET==1,])){
	FUMAp$ReplicatedMET[grep(pattern=snptable$SNP[snptable$ReplicatedMET==1][i], x=FUMAp$IndSigSNPs)]<-1
}

for (i in 1:nrow(snptable[snptable$ReplicatedKET==1,])){
        FUMAp$ReplicatedKET[grep(pattern=snptable$SNP[snptable$ReplicatedKET==1][i], x=FUMAp$IndSigSNPs)]<-1
}

for (i in 1:nrow(snptable[snptable$ReplicatedMETKET==1,])){
        FUMAp$ReplicatedMETKET[grep(pattern=snptable$SNP[snptable$ReplicatedMETKET==1][i], x=FUMAp$IndSigSNPs)]<-1
}

outdir="/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication/replicatedtable"

write.table(FUMAp, paste(outdir, "/FUMA-", pheno1[p], "-", "UKBEUR-MET-KET-replication.txt", sep=""),
                quote=F, row.names=F)


row<-c(nrow(MET), nrow(METsig),	sum(FUMAp$ReplicatedMET), nrow(KET), nrow(KETsig),sum(FUMAp$ReplicatedKET),
		nrow(MKM), nrow(MKMsig), sum(FUMAp$ReplicatedMETKET))

print(row)

}
