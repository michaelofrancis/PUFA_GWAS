library(plyr)
library(tidyverse)
library(reshape)
library(data.table)

# dat<-as_tibble(read.csv("/Users/mike/Documents/Research/PUFA-GWAS/discovery/SMultiXCan/genecount-summary.csv"))
# 
# colnames(dat)<-c("Phenotype", "Significant_nonnovel", "Novel", "NovelNotGWAS")
# dat
# 
# dat$Significant_nonnovel<-dat$Significant_nonnovel-dat$Novel
# dat$Novel<-dat$Novel-dat$NovelNotGWAS
# 
# dat$Phenotype<-factor(dat$Phenotype, levels = rev(dat$Phenotype))
# 
# 
# datlong<-as_tibble(melt(setDT(dat), 
#                         id.vars = c("Phenotype")))
# datlong$variable<-factor(datlong$variable, levels = c("NovelNotGWAS", "Novel", 
#                                               "Significant_nonnovel"))
# datlong$variable
# 
# datlong$Phenotype
# 
# ggplot(datlong, aes(x=Phenotype,y=value, fill=variable)) + 
#     geom_bar(position="stack", stat="identity")+
#     coord_flip()+
#     theme_bw()+ylab("CountGenes")+
#     ggtitle("SMultiXcan Significant genes summary")
# 
# 


# Redo with loci ------------------------------------------------

dat2<-as_tibble(read.csv("/Users/mike/Documents/Research/PUFA-GWAS/discovery/SMultiXCan/SMul3.csv"))
dat2
dat2%>%group_by(Phenotype)

dat2$InGWAS<-1
dat2$InGWAS[dat2$SMul_in_FUMAd==0 & dat2$SMul_in_FUMAm==0]<-0
table(dat2$InGWAS)

ph<-unique(dat2$Phenotype)
mydat<-as.data.frame(matrix(nrow=14, ncol=4))
for (p in 1:length(ph)){
    x<-dat2[dat2$Phenotype==ph[p],]
    x<-x%>%select(LocusID, Novel1MB, InGWAS)
    numloc<-length(unique(x$LocusID))
    novel<-length(unique(x$LocusID[x$Novel1MB==1]))
    novelnotgwas<-length(unique(x$LocusID[x$Novel1MB==1 & x$InGWAS==0]))
    mydat[p,]<-c(ph[p], numloc, novel, novelnotgwas)
}


mydat[2:4]<-sapply(mydat[2:4], as.numeric)
colnames(mydat)<-c("Phenotype", "Known", "Novel", "Novel + not in GWAS results")
mydat$Known<-mydat$Known-mydat$Novel
mydat$Novel<-mydat$Novel-mydat$`Novel + not in GWAS results`

mydatlong<-as_tibble(melt(setDT(mydat), 
                        id.vars = c("Phenotype")))
mydatlong$variable<-factor(mydatlong$variable, 
                         levels = rev(c("Novel + not in GWAS results", "Novel", 
                                              "Known")))
mydatlong$variable
phenotypenames<-c("Omega-3", "Omega-3 (TFAP)", "Omega-6" , "Omega-6 (TFAP)",
                  "Omega-6:Omega-3 ratio", "DHA", "DHA (TFAP)", "LA", "LA (TFAP)",
                  "PUFAs", "PUFAs (TFAP)", "MUFAs", "MUFAs (TFAP)", "PUFAs:MUFAs ratio")

ph2<-c("n-3", "n-3 pct", "n-6" , "n-6 pct",
       "n-6:n-3", "DHA", "DHA pct", "LA", "LA pct",
       "PUFAs", "PUFAs pct", "MUFAs", "MUFAs pct", "PUFAs:MUFAs")

mydatlong$Phenotype<-mapvalues(mydatlong$Phenotype, from=phenotypenames, to=ph2)

mydatlong$Phenotype<-factor(mydatlong$Phenotype, levels=rev(ph2))

colnames(mydatlong)[2]="Novelty (UKB only)"


tiff("4D.SMultiscan-summary.tiff",
     width = 7, height = 5,
     units = 'in', res = 600)

ggplot(mydatlong, aes(x=Phenotype,y=value, fill=`Novelty (UKB only)`)) + 
    geom_bar(position="stack", stat="identity")+
    coord_flip()+
    theme_bw()+
    scale_fill_manual(values = c("#2F7FBA", "#BA1613", "#5ABA2F"), 
                      labels=function(x) str_wrap(x, width =15))+
    ylab("Count significant S-MultiXcan loci")
dev.off()


