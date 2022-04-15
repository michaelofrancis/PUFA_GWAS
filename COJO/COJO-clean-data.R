library(plyr)
library(tidyverse)

#dir="/Users/mike/Documents/Research/PUFA-GWAS/COJO/jma.cojo/"
dir="/Users/mike/Documents/Research/PUFA-GWAS/discovery/jma"
list.files(dir)

files<-list()
for (i in 1:length(list.files(dir))){
    
    files[[i]]<-as_tibble(read_delim(
        list.files(dir, full.names = T)[i]), header=T)
    files[[i]]<-files[[i]][-grep("Chr", files[[i]]$Chr),]
    files[[i]]$Phenotype<-str_split(
        list.files(dir), ".jma.cojo", simplify = T)[i,1]
    files[[i]]<-files[[i]]%>%select(Phenotype, everything())
    }
files<-do.call(rbind, files)




files$freq<-as.numeric(files$freq)
files$freq_geno<-as.numeric(files$freq_geno)
dif<-files$freq-files$freq_geno

range(dif)
#[1] -0.048282  0.024550


pheno<-c("w3FA_NMR_resinv", "w3FA_NMR_TFAP_resinv", "w6FA_NMR_resinv",
         "w6FA_NMR_TFAP_resinv", "w6_w3_ratio_NMR_resinv",
         "DHA_NMR_resinv", "DHA_NMR_TFAP_resinv", "LA_NMR_resinv",
         "LA_NMR_TFAP_resinv", "PUFA_NMR_resinv", "PUFA_NMR_TFAP_resinv",
         "MUFA_NMR_resinv", "MUFA_NMR_TFAP_resinv", "PUFA_MUFA_ratio_NMR_resinv")


phenotypenames<-c("Omega-3", "Omega-3 (TFAP)", "Omega-6" , "Omega-6 (TFAP)",
                  "Omega-6:Omega-3 ratio", "DHA", "DHA (TFAP)", "LA", "LA (TFAP)",
                  "PUFAs", "PUFAs (TFAP)", "MUFAs", "MUFAs (TFAP)", "PUFAs:MUFAs ratio")

ph<-c("Omega-3 Fatty Acids",
      "Omega-3 Fatty Acids to Total Fatty Acids percentage", 
      "Omega-6 Fatty Acids",
      "Omega-6 Fatty Acids to Total Fatty Acids percentage",
      "Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio",
      "Docosahexaenoic Acid", "Docosahexaenoic Acid to Total Fatty Acids percentage",
      "Linoleic Acid","Linoleic Acid to Total Fatty Acids percentage",
      "Polyunsaturated Fatty Acids", 
      "Polyunsaturated Fatty Acids to Total Fatty Acids percentage",
      "Monounsaturated Fatty Acids","Monounsaturated Fatty Acids to Total Fatty Acids percentage",
      "Polyunsaturated Fatty Acids to Monounsaturated Fatty Acids ratio")

files$Phenotype<-mapvalues(files$Phenotype,from = pheno, to=ph)


write.csv(files, paste(dir, "/discovery.combined.jma.cojo.csv", sep=""), row.names = F,
          quote=F)

files<-as_tibble(read.csv(paste(dir, "/discovery.combined.jma.cojo.csv", sep="")))
files%>%group_by(Phenotype)%>%dplyr::summarise(n())


# meta-analysis -------------------------------------------------

dir<-"/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/jma"

list.files(dir)

files<-list()
for (i in 1:length(list.files(dir))){
    
    files[[i]]<-as_tibble(read_delim(
        list.files(dir, full.names = T)[i]), header=T)
    files[[i]]<-files[[i]][-grep("Chr", files[[i]]$Chr),]
    files[[i]]$Phenotype<-str_split(
        list.files(dir), ".jma.cojo", simplify = T)[i,1]
    files[[i]]<-files[[i]]%>%select(Phenotype, everything())
}
files<-do.call(rbind, files)




files$freq<-as.numeric(files$freq)
files$freq_geno<-as.numeric(files$freq_geno)
dif<-files$freq-files$freq_geno


range(dif)
#[1] -0.048282  0.024550
files

pheno<-c("FAw3", "FAw6", "DHA", "LA", "MUFA")

ph<-c("Omega-3 Fatty Acids",
      "Omega-6 Fatty Acids",
      "Docosahexaenoic Acid", 
      "Linoleic Acid",
      "Monounsaturated Fatty Acids")

files$Phenotype<-mapvalues(files$Phenotype,from = pheno, to=ph)


write.csv(files, paste(dir, "/meta-analysis.combined.jma.cojo.csv", sep=""), row.names = F,
          quote=F)


files%>%group_by(Phenotype)%>%dplyr::summarise(n())
