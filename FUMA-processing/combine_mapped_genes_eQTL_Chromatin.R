library(plyr)
library(tidyverse)

madir<-"/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/FUMA-out-livermapped"
dirs<-list.dirs(madir, recursive = F)
phenotypes<-c("Docosahexaenoic acid","Linoleic acid","Monounsaturated fatty acids",
              "Omega-3 fatty acids", "Omega-6 fatty acids")

genes<-list()
for (i in 1:length(dirs)){
    genes[[i]]<-as_tibble(read_delim(
            paste(dirs[i],"/genes.txt",sep="")))
    genes[[i]]$Phenotype<-phenotypes[i]
    genes[[i]]<-genes[[i]]%>%select(Phenotype, everything())
}

map<-do.call(rbind, genes)
write.csv(map, "/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/FUMA-out-livermapped/mapped_genes.csv",
            quote=F, row.names=F)


# Pt 2, join to table -------------------------------------------

map<-as_tibble(read.csv("/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/FUMA-out-livermapped/mapped_genes.csv"))

FUMA<-as_tibble(read.csv("/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/alleles.csv"))
FUMA$BP<-as.numeric(gsub(",", "", FUMA$BP))

FUMA$Phenotype<-mapvalues(FUMA$Phenotype, from=unique(FUMA$Phenotype), to=unique(map$Phenotype))
FUMA

#map<-as_tibble(read.csv("/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/S13.Mapped_Genes.csv"))

identical(unique(FUMA$Phenotype), unique(map$Phenotype))
#^^MUST BE TRUE


FUMA$Gene<-""
FUMA$Gene_type<-""
for (i in 1:nrow(FUMA)){
    for (j in 1:nrow(map)){
        #If same phenotype, and same genomic locus
        if (FUMA[[i,"Phenotype"]]==map[j, "Phenotype"] & (FUMA[[i, "GenomicLocus"]] == map[j, "GenomicLocus"])){ 
            
                if (FUMA[[i, "Gene"]]==""){        
                    FUMA[[i, "Gene"]]<-map[[j,"symbol"]]
                    FUMA[[i, "Gene_type"]]<-map[[j,"type"]]
                    
                    
                }else 
                    FUMA[[i, "Gene"]]<-paste(
                        FUMA[[i, "Gene"]],map[[j,"symbol"]],
                        sep="; ") 
                
                FUMA[[i, "Gene_type"]]<-paste(
                    FUMA[[i, "Gene_type"]],map[[j,"type"]],
                    sep="; ") 
            }#end inner if
        }#end inner for
    }#end outer for



FUMA

write.csv(FUMA,"/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/FUMA-out-livermapped/FUMA-maUKBEURMETKET.withmappedgenes.csv",
          row.names=F, quote=F)
