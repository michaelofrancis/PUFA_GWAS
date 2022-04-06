library(tidyverse)
library(MungeSumstats)

dir<-"/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication"

pheno<-c("FAw3", "FAw6", "DHA", "LA", "MUFA")

pheno2<-c("w3FA", "w6FA", "DHA", "LA", "MUFA")

#REDO munge METSIM with correct EAF column
for (i in 1:length(pheno)){
#METSIM-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#i=1
#p=1
#Udir<-"/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication/METSIM/UMich"
#U<-as_tibble(read_delim(paste(Udir,"/phenocode-", pheno[p], "_combined.tsv.gz", sep="")))

metfile<-paste(dir, "/METSIM/Yfiles/Locke_", pheno[i], ".txt", sep="")
MET<-as_tibble(read.table(metfile, header=T))
MET2<-MET%>%select(CHROM, BEG, NS, EAF, PVALUE, BETA, SEBETA, SNP, NonRefAllele, RefAllele)
colnames(MET2)<-c("CHR","BP", "NSTUDY", "EAF", "PVALUE", "BETA", "SE", "SNP", "EFFECT_ALLELE", "NON_EFFECT_ALLELE")

reformattedMET<-MungeSumstats::format_sumstats(path=MET2, ref_genome="GRCh37",
                        INFO_filter=0.29, snp_ids_are_rs_ids=F)

refoMET<-as_tibble(read.table(reformattedMET, header=T))

write.table(refoMET, paste(dir, "/munge/METSIM-EAF/", pheno[i],".txt", sep=""), quote=F, row.names=F)
}
