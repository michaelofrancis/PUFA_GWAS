library(plyr)
library(tidyverse)

pheno<-c("w3FA", "w6FA", "DHA", "LA", "MUFA")
pheno2<-c("W3", "w6_METAL", "DHA", "LA_METAL","MUFA_METAL")


FUMAdisc<-as_tibble(read.csv("/work/kylab/mike/PUFA-GWAS/FUMA-final/GenomicLoci-DISC-02102022.csv"))
FUMAma<-as_tibble(read.csv("/work/kylab/mike/PUFA-GWAS/FUMA-final/GenomicLoci-MA-UKBEUR-MET-KET-02102022.csv"))

FUMAma$Phenotype<-mapvalues(FUMAma$Phenotype, from=pheno2, to=pheno, warn_missing = TRUE)


#filein<-list()
#for (j in 1:length(pheno)){
#j=1
#	filein[[j]]<-as_tibble(read.table(paste("/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication/munge/UKBN/",
#                pheno[j], "_NMR_resinv.M2.txt", sep=""), header=T))
#	filein[[j]]$Phenotype<-pheno[j] 
 #       filein[[j]]<-filein[[j]]%>%select(Phenotype, everything())
#}

#names(filein)<-pheno

files<-c("META_IVW_FAw3.STDERR-N1.txt.gz","META_IVW_FAw6.STDERR-N1.txt.gz",
"META_IVW_DHA.STDERR1.txt.gz","META_IVW_LA.STDERR1.txt.gz","META_IVW_MUFA.STDERR1.txt.gz")

metaldir<-"/scratch/mf91122/PUFA-GWAS/METAL/M2-N/gz"

metalfile<-list()
for (j in 1:length(pheno)){
	metalfile[[j]]<-as_tibble(read_delim(	
			paste(metaldir, "/", files[j], sep="")
				))
	metalfile[[j]]$Phenotype<-pheno[j]
        metalfile[[j]]<-metalfile[[j]]%>%select(Phenotype, everything())
	}

metalfile2<-do.call(rbind, metalfile)

colnames(FUMAma)[4:6]<-c("MarkerName", "CHR", "BP")

join<-left_join(FUMAma, metalfile2)
#Joining, by = c("Phenotype", "CHR", "BP")

write.csv(join, "/work/kylab/mike/PUFA-GWAS/FUMA-final/table-process/alleles.csv",
	row.names=F, quote=F)





files<-c("META_IVW_FAw3.STDERR1.POSTPROCESS.csv",
    "META_IVW_FAw6.STDERR1.POSTPROCESS.csv",
    "META_IVW_DHA.STDERR1.POSTPROCESS.csv",
"META_IVW_LA.STDERR1.POSTPROCESS.csv",
"META_IVW_MUFA.STDERR1.POSTPROCESS.csv")

metaldir<-"/scratch/mf91122/PUFA-GWAS/METAL/M2/STDERR-postprocess"

metalfile<-list()
for (j in 1:length(pheno)){
        metalfile[[j]]<-as_tibble(read_csv(
                        paste(metaldir, "/", files[j], sep="")
                                ))
        metalfile[[j]]$Phenotype<-pheno[j]
        metalfile[[j]]<-metalfile[[j]]%>%select(Phenotype, everything())
        }
metalfile2<-do.call(rbind, metalfile)
