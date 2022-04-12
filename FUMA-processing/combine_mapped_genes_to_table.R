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

genes<-do.call(rbind, genes)
write.csv(genes, "/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/FUMA-out-livermapped/mapped_genes.csv",
            quote=F, row.names=F)


# Pt 2, join to table -------------------------------------------


FUMA<-as_tibble(read.csv("/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/alleles.csv"))
FUMA$BP<-as.numeric(gsub(",", "", FUMA$BP))

#map<-as_tibble(read.csv("/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/S13.Mapped_Genes.csv"))
map<-as_tibble(read.csv("/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/FUMA-out-livermapped/mapped_genes.csv"))

FUMA$Gene<-""
FUMA$Gene_type<-""

for (i in 1:nrow(FUMA)){
    for (j in 1:nrow(map)){
        
        if (FUMA[[i,"Phenotype"]]==map[j, "Phenotype"]){ 
            #If same phenotype, two situations where FUMA pos can be in map pos
            
            if(
                (FUMA[[i, "pos"]] >= map[j, "start"] & FUMA[[i, "pos"]] <=map[j, "end"]) 
            ){
                # cat(
                #     paste(
                #         "i: ", i,
                #         "j: ", j,
                #     "\nFUMA pheno:", FUMA[i,"Phenotype"],
                #     "\nmap pheno:", map[j,"Phenotype"],
                #     "\nFUMA pos:", FUMA[i,"pos"], 
                #     "\nmap start:", map[j,"start"],
                #     "\nmap end:", map[j,"end"],
                #     "\nmap gene:", map[j,"symbol"],
                #     "\n\n"
                #     )
                # )
                
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
        }#end outer if
    }
}

FUMA<-FUMA%>%select(-Nearest.gene)
write.csv(FUMA,"/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/FUMA-out-livermapped/FUMA-maUKBEURMETKET.withmappedgenes.csv",
          row.names=F, quote=F)
