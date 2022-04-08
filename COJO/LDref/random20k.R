unrelatedEUR<-read.table("/work/kylab/mike/PUFA-GWAS/GCTA-COJO2/1.Bfile20k/ukb48818_rel_s488282_output.dat")

sampleQC<-read.table("/scratch/mf91122/PUFA-GWAS/BOLT-LMM/genotypeQC/chr1.sample",
                header=F,skip=2)

rand<-unrelatedEUR$V1[! unrelatedEUR$V1 %in% sampleQC$V1]

rand2<-rand[sample(length(rand), 20000)]


rand<-sampleQC$V1[! sampleQC$V1 %in% unrelatedEUR$V1]
set.seed(42)
rand2<-rand[sample(length(rand), 20000)]

rand3<-data.frame("FID"=rand2, "IID"=rand2)

write.table(rand3, "/scratch/mf91122/PUFA-GWAS/GCTA-COJO2/1.Bfile20k/LDref_20kQCunrelated.fam",
                quote=F, row.names=F)
