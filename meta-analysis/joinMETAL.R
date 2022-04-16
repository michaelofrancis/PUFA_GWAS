library(tidyverse)
library(MungeSumstats)

dir<-"/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication"

pheno<-c("FAw3", "FAw6", "DHA", "LA", "MUFA")

pheno2<-c("w3FA", "w6FA", "DHA", "LA", "MUFA")

#i=1
for (i in 1:length(pheno)){

join3<-as_tibble(read.table( 
	paste("/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication/munge/UKBEURKETMET.mungecombined.metalinput/", 
		pheno[i], ".UKBEURKETMET.mungedcombined.metalinput.txt", sep=""),
	header=T))


METAL<-as_tibble(read.table(
        paste("/scratch/mf91122/PUFA-GWAS/METAL/mungecombined/META_IVW_", pheno[i],"1.txt",sep=""),
        header=T))


join3<-join3[!is.na(join3$P_KET) | !is.na(join3$P_MET),]



#METAL alleles are default lowercase, set to upper
METAL$Allele1<-toupper(METAL$Allele1)
METAL$Allele2<-toupper(METAL$Allele2)

#Tag metal columns before join
colnames(METAL)[2:11]<-paste(colnames(METAL)[2:11], "_METAL", sep="")


joinM<-left_join(join3, METAL, by=c("SNP"="MarkerName"))

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

		


#joinM<-as_tibble(read.table(
 #       paste("/scratch/mf91122/PUFA-GWAS/sumstats-final/UKBEURKETMET/", pheno[i], 
#	".UKBEURKETMET.txt", sep=""),header=T))

#Prepare LDSC format for meta-analysis results
	ldsc<-joinM%>%select(SNP, CHR, BP, EA,NEA, 
		Freq1_METAL, Effect_METAL, StdErr_METAL, P.value_METAL, Ntotal)
	colnames(ldsc)<-c("SNP", "CHR","BP", "EFFECT_ALLELE", "NON_EFFECT_ALLELE",
                "EFFECT_ALLELE_FRQ","BETA","SE", "P", "N")

        reformatted<-MungeSumstats::format_sumstats(path=ldsc, ref_genome="GRCh37", ldsc_format=T,
                INFO_filter=0.3)

        refo<-as_tibble(read.table(reformatted, header=T))
        #Remove infinite Z values
        refo<-refo[!is.infinite(refo$Z),]

	write.table(refo, paste("/scratch/mf91122/PUFA-GWAS/sumstats-final/UKBEURKETMET/ldsc/", pheno[i], 
		".UKBEURKETMET.txt", sep=""),
                quote=F, row.names=F)


#Prepare .ma file for COJO

ma<-joinM%>%select(SNP, EA,NEA, Freq1_METAL, Effect_METAL, StdErr_METAL, P.value_METAL, Ntotal)

colnames(ma)<-c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")

#------------------------
#Line to fix the GCTA-COJO precision bug
#See: https://gcta.freeforums.net/thread/518/cojo-slct-selecting-lowest-value
ma$se<-abs(ma$b/qnorm(ma$p/2, lower.tail=F))
ma$se[is.na(ma$se)]<-0
#range(ma$p[(ma$se2-ma$se>0.05)])
#[1] 0.9991 0.9999
#------------------------


write.table(ma,
	paste("/scratch/mf91122/PUFA-GWAS/sumstats-final/UKBEURKETMET/maCOJO/", pheno[i],
        ".UKBEURKETMET.ma", sep=""),quote=F, row.names=F)


#Prepare FUMA input--------------------------------------------------------------------------------------
#Just METAL columns into the FUMA input since it's too big with the full file

FUMAin<-joinM%>%select(SNP, CHR, BP, NEA,EA,Freq1_METAL,FreqSE_METAL, Effect_METAL, StdErr_METAL,P.value_METAL, Ntotal)
colnames(FUMAin)<-c("SNP", "CHR", "BP", "non_effect_allele", "effect_allele", "FRQ", "FRQ_SE", "Beta", "SE", "P", "N")

FUMAin$BP<-as.integer(FUMAin$BP)

FUMAin$P[FUMAin$P<1e-300]=1e-300

gz<-gzfile(paste(
	"/scratch/mf91122/PUFA-GWAS/sumstats-final/UKBEURKETMET/FUMAin/",
                pheno[i], ".UKBEURKETMET.FUMAin.txt.gz", sep=""))

write.table(FUMAin,gz,quote=F, row.names=F, fileEncoding = "ascii")

}

#INPUT FOR CMPLOT (FIG 2)--------------------

joinM<-list()
for (i in 1:length(pheno)){
joinM[[i]]<-as_tibble(read.table(paste("/scratch/mf91122/PUFA-GWAS/sumstats-final/UKBEURKETMET/", pheno[i], ".UKBEURKETMET.txt", $
        header=T))
}

dat<-joinM[[1]][1:3]
dat$P_FAw3<-joinM[[1]]$P.value_METAL

dat<-joinM[[2]]%>%select(SNP, CHR, BP, P.value_METAL)%>%inner_join(dat)
colnames(dat)[4]<-"P_FAw6"

dat<-joinM[[3]]%>%select(SNP, CHR, BP, P.value_METAL)%>%inner_join(dat)
colnames(dat)[4]<-"P_DHA"

dat<-joinM[[4]]%>%select(SNP, CHR, BP, P.value_METAL)%>%inner_join(dat)
colnames(dat)[4]<-"P_LA"

dat<-joinM[[5]]%>%select(SNP, CHR, BP, P.value_METAL)%>%inner_join(dat)
colnames(dat)[4]<-"P_MUFA"

dat<-dat%>%select(SNP, CHR, BP, paste("P_", pheno, sep=""))

write.csv(dat, "/scratch/mf91122/PUFA-GWAS/sumstats-final/UKBEURKETMET/CMplot-04162022.csv",
        row.names=F, quote=F)


