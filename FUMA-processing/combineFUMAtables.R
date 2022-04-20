suppressMessages(library(tidyverse))

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/discovery/resultsFix/FUMA-M2"

discdirs<-list.dirs(dir, recursive = F)
phen<-str_split(discdirs, "/", simplify=T)[,10]

disc<-list()
for (i in 1:length(discdirs)){
    disc[[i]]<-as_tibble(read.table(paste(discdirs[i], "/GenomicRiskLoci.txt",
                               sep=""), 
                         header=T))
    disc[[i]]$Phenotype<-phen[i]
    
}
disc2<-do.call("rbind", disc)
disc2<-disc2%>%select(Phenotype, everything())
disc2

write.csv(disc2, 
          paste(dir, "/GenomicLoci-DISC-04092022.csv", sep=""),
            quote=F, row.names=F)


#META ANALYSIS

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/FUMA-out-04132022"
madirs<-list.dirs(dir, recursive = F)

maphen<-str_split(madirs, "/", simplify=T)[,9]

ma<-list()
for (i in 1:length(madirs)){
    ma[[i]]<-as_tibble(read.table(paste(madirs[i], "/GenomicRiskLoci.txt",
                                          sep=""), 
                                    header=T))
    ma[[i]]$Phenotype<-maphen[i]
    
}

ma2<-do.call("rbind", ma)
ma2<-ma2%>%select(Phenotype, everything())
ma2


write.csv(ma2, 
          paste(dir, "/GenomicLoci-MA-UKBEUR-MET-KET-04132022.csv", sep=""),
          quote=F, row.names=F)
