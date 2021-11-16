library(tidyverse)

setwd("/work/kylab/mike/PUFA-GWAS/PLINKPCA/addPCtoPheno")

#files<-list.files("/scratch/mf91122/PUFA-GWAS/pheno", full.names=TRUE)[1:2] #CHANGE when M2

files<-c("/scratch/mf91122/PUFA-GWAS/pheno/PUFA_GWAS_pheno_M1.txt",
	"/scratch/mf91122/PUFA-GWAS/pheno/PUFA_GWAS_pheno_M2.txt")

pca<-read.table("/scratch/mf91122/PUFA-GWAS/PLINKPCA/outpc/fullpca.eigenvec", 
		header=FALSE, stringsAsFactors=FALSE)

pca<-as_tibble(pca)

colnames(pca)<-c("FID","IID", sprintf("PC%s", 1:20))


#Join files
for (i in 1:length(files)){
	phen<-read.table(files[i], header=TRUE, stringsAsFactors=FALSE)	
	phen<-as_tibble(phen)
	join<-left_join(phen, pca,
                                by=intersect(colnames(phen),
                                colnames(pca)))
        join<-join%>%select(FID, IID, everything())
	write.table(join, paste(files[i], ".pc.txt", sep=""), quote=FALSE, row.names=FALSE)
}
