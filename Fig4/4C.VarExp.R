library(plyr)
suppressMessages(library(tidyverse))
library(reshape)
library(data.table)

# 2. Meta-analysis ----------------------------------------------

FUMAmaFile<-"/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/meta-analysis-loci_05092022.csv"
FUMAma<-as_tibble(read.csv(FUMAmaFile, skip=1))
COJOmaFile<-"/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/jma/meta-analysis.combined.jma.cojo.rep.csv"
COJOma<-as_tibble(read.csv(COJOmaFile))

ph<-c("Omega-3 fatty acids","Docosahexaenoic acid" ,"Omega-6 fatty acids",
       "Linoleic acid","Monounsaturated fatty acids")
phc<-c("FAw3","DHA","FAw6","LA","MUFA")

COJOma$Phenotype<-mapvalues(COJOma$Phenotype, from=phc, to=ph)


# FUMA ----------------------------------------------------------

ve.ma<-matrix(data=NA, nrow=5, ncol=7)

for (i in 1:length(ph)){
    ve.ma[i,1]<-(ph[i])
    x<-FUMAma[FUMAma$Phenotype==ph[i],]
    freq<-x$FRQ
    beta<-x$BETA
    var<-(sum(2*freq*(1-freq)*beta^2)*100)
    ve.ma[i,3]<-(var)
    
}
FUMA.count<-FUMAma%>%group_by(Phenotype)%>%summarise(N=n())
FUMA.count<-FUMA.count[c(4,1,5,2,3),]
ve.ma[,2]<-unlist(FUMA.count[,2])

#Known only
FUMA.known<-FUMAma%>%filter(Novel.locus...1Mb.from.any.known.locus==0)

for (i in 1:length(ph)){
    x<-FUMA.known[FUMA.known$Phenotype==ph[i],]
    freq<-x$FRQ
    beta<-x$BETA
    var<-(sum(2*freq*(1-freq)*beta^2)*100)
    ve.ma[i,5]<-(var)
}


FUMA.count.known<-FUMAma%>%group_by(Phenotype)%>%
    filter(Novel.locus...1Mb.from.any.known.locus==0)%>%summarise(N=n())
FUMA.count.known<-FUMA.count.known[c(4,1,5,2,3),]
ve.ma[,4]<-FUMA.count.known$N


# COJO ----------------------------------------------------------

for (i in 1:length(ph)){
    x<-COJOma[COJOma$Phenotype==ph[i],]
    freq<-x$freq
    beta<-x$bJ
    var<-(sum(2*freq*(1-freq)*beta^2)*100)
    ve.ma[i,7]<-(var)
}

COJO.count<-COJOma%>%group_by(Phenotype)%>%summarise(N=n())
COJO.count<-COJO.count[c(4,1,5,2,3),]
ve.ma[,6]<-unlist(COJO.count$N)

ve.ma<-as.data.frame(ve.ma)
ve.ma

#w3,dha,w6,la,mufa
h2<-c(0.1296, 0.1175, 0.1538, 0.1601,0.1636)
h2se<-c(0.0215,0.0223,0.0209,0.0212,0.0227)
ve.ma$h2<-h2*100
h2ci<-h2se*100*1.96 #1.96 for 95% confidence interval
ve.ma

colnames(ve.ma)<-c("Phenotype", "Count genomic loci (FUMA)",
                   "All genomic loci (FUMA)", "FUMA.count.known", 
                   "Known genomic loci (FUMA)", 
                   "Count independent variants (COJO)",
                   "Independent variants (COJO)",
                   "SNP-based heritability (h2; LDSC)")

vars<-c("Known genomic loci (FUMA)", "All genomic loci (FUMA)", 
        "Independent variants (COJO)", "SNP-based heritability (h2; LDSC)")


ph<-c("Omega-3 fatty acids","Docosahexaenoic acid" ,"Omega-6 fatty acids",
      "Linoleic acid","Monounsaturated fatty acids")
ph2<-c("n-3","DHA","n-6","LA","MUFAs")

ve.ma$Phenotype<-mapvalues(ve.ma$Phenotype, from=ph, to=ph2)


plotme<-as_tibble(melt(setDT(ve.ma%>%
            select(Phenotype, vars)), 
            id.vars = c("Phenotype"), 
            variable.name = "var"))
plotme$value<-as.double(plotme$value)
plotme$ci<-NA
plotme$ci[plotme$var=="SNP-based heritability (h2; LDSC)"]<-h2ci
plotme$Phenotype<-factor(plotme$Phenotype, levels=ph2)

tiff("4C.VarExp.tiff", 
     width = 6, height = 5, 
     units = 'in', res = 600)

ggplot(plotme, aes(fill=var, y=value, x=Phenotype)) + 
    geom_bar(stat = "identity", width=0.7, 
             position = position_dodge(width = 0.8))+
    geom_errorbar(aes(ymin = value - ci, ymax = value + ci),
                  position = position_dodge(width = 0.8),
                  na.rm=T, size=0.5,
                  width=0.5)+
    ylab("Explained phenotypic variance (%)")+
    theme_bw()+
    labs(fill = "")+
    scale_fill_manual(values = c("#6DA1C6", "#87130F", "#25BA55", "#9463C6"), 
                      labels=function(x) str_wrap(x, width =15))+
    theme(axis.text.x = element_text(angle = 45, 
                        vjust = 1, hjust=1))

dev.off()
