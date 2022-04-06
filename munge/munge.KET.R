library(MungeSumstats)
suppressMessages(library(dplyr))

dir<-"/scratch/mf91122/PUFA-GWAS/PUFA-GWAS-replication"

pheno<-c("FAw3", "FAw6", "DHA", "LA", "MUFA")

for (i in 1:length(pheno)){

#Kettunen
ketfile<-paste(dir, "/Kettunen/", "MAGNETIC_", pheno[i], ".txt", sep="")
KET<-as_tibble(read.table(ketfile, header=T))

#Munge ket
KET2<-KET
colnames(KET2)<-c("CHR","BP","ID","EFFECT_ALLELE","NON_EFFECT_ALLELE","EFFECT_ALLELE_FREQ","BETA","SE","P","NSTUDY","N")

reformattedKET<-MungeSumstats::format_sumstats(path=KET2, ref_genome="GRCh37", INFO_filter=0.29)
refoKET<-as_tibble(read.table(reformattedKET, header=T))

write.table(refoKET, paste(dir, "/munge/Kettunen/", pheno[i],".txt", sep=""), quote=F, row.names=F)

}
