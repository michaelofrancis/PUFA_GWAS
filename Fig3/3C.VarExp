library(plyr)
suppressMessages(library(tidyverse))
library(reshape)
library(data.table)


# #Variance explained = (2*fre*(1-fre)*beta^2)
# 
# #1. Discovery--COJO
# #2. Meta-analysis --COJO
# 
# # 1. Discovery -----------------------------------------------------
# 
# 
# ph<-c("Omega-3 Fatty Acids",
#       "Omega-3 Fatty Acids to Total Fatty Acids percentage", 
#       "Omega-6 Fatty Acids",
#       "Omega-6 Fatty Acids to Total Fatty Acids percentage",
#       "Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio",
#       "Docosahexaenoic Acid", "Docosahexaenoic Acid to Total Fatty Acids percentage",
#       "Linoleic Acid","Linoleic Acid to Total Fatty Acids percentage",
#       "Polyunsaturated Fatty Acids", 
#       "Polyunsaturated Fatty Acids to Total Fatty Acids percentage",
#       "Monounsaturated Fatty Acids","Monounsaturated Fatty Acids to Total Fatty Acids percentage",
#       "Polyunsaturated Fatty Acids to Monounsaturated Fatty Acids ratio")
# 
# #Load 
# FUMAdFile<-"/Users/mike/Documents/Research/PUFA-GWAS/discovery/FUMA.disc.05182022.csv"
# FUMAd<-as_tibble(read.csv(FUMAdFile))
# FUMAd
# COJOdFile="/Users/mike/Documents/Research/PUFA-GWAS/discovery/jma/discovery.combined.jma.cojo.csv"
# COJOd<-as_tibble(read.csv(COJOdFile))
# #Phenotype names are already standardized
# 
# # FUMA ----------------------------------------------------------
# 
# ve.d<-matrix(data=NA, nrow=14, ncol=7)
# 
# for (i in 1:length(ph)){
#     ve.d[i,1]<-(ph[i])
#     x<-FUMAd[FUMAd$Phenotype==ph[i],]
#     freq<-x$FREQ
#     beta<-x$BETA
#     var<-(sum(2*freq*(1-freq)*beta^2)*100)
#     ve.d[i,3]<-(var)
# }
# 
# ve.d
# FUMA.count<-FUMAd%>%group_by(Phenotype)%>%summarize(N=n())
# ph2<-as.data.frame(ph)
# FUMA.count<-left_join(ph2,FUMA.count,  by=c("ph"="Phenotype"))
# ve.d[,2]<-unlist(FUMA.count[,2])
# ve.d
# 
# #Known only
# FUMA.known<-FUMAd%>%filter(Novel.locus...1Mb.from.any.known.locus==0)
# 
# for (i in 1:length(ph)){
#     x<-FUMA.known[FUMA.known$Phenotype==ph[i],]
#     freq<-x$FREQ
#     beta<-x$BETA
#     var<-(sum(2*freq*(1-freq)*beta^2)*100)
#     ve.d[i,5]<-(var)
# }
# 
# ve.d
# FUMA.count.known<-FUMAd%>%group_by(Phenotype)%>%
#     filter(Novel.locus...1Mb.from.any.known.locus==0)%>%summarize(N=n())
# FUMA.count.known<-left_join(ph2,FUMA.count.known,  by=c("ph"="Phenotype"))
# ve.d[,4]<-FUMA.count.known$N
# 
# 
# # COJO ----------------------------------------------------------
# 
# for (i in 1:length(ph)){
#     x<-COJOd[COJOd$Phenotype==ph[i],]
#     freq<-x$freq
#     beta<-x$bJ
#     var<-(sum(2*freq*(1-freq)*beta^2)*100)
#     ve.d[i,7]<-(var)
# }
# 
# COJO.count<-COJOd%>%group_by(Phenotype)%>%summarize(N=n())
# COJO.count<-left_join(ph2,COJO.count,  by=c("ph"="Phenotype"))
# ve.d[,6]<-unlist(COJO.count$N)
# 
# ve.d<-as.data.frame(ve.d)
# ve.d
# 
# #Copy vals from S4
# h2.d<-c(0.148,0.1253,0.1724,0.1416,0.1269,0.1325,0.1372,0.1722,0.1231,0.1834,0.1492,0.1849,0.2044,0.1882)
# h2se.d<-c(0.0234,0.0232,0.0234,0.0184,0.0225,0.0233,0.0245,0.0234,0.0183,0.0259,0.0196,0.0265,0.0324,0.0267)
# ve.d$h2<-h2.d*100
# ve.d$h2ci<-h2se.d*100*1.96 #1.96 for 95% confidence interval
# ve.d
# 
# colnames(ve.d)<-c("Phenotype", "Count genomic loci (FUMA)",
#                    "All genomic loci (FUMA)", "FUMA.count.known", 
#                    "Known genomic loci (FUMA)", 
#                    "Count independent variants (COJO)",
#                    "Independent variants (COJO)",
#                    "SNP-based heritability (h2; LDSC)", "h2ci")
# 
# vars<-c("Known genomic loci (FUMA)", "All genomic loci (FUMA)", 
#         "Independent variants (COJO)", "SNP-based heritability (h2; LDSC)")
# 
# plotme<-as_tibble(melt(setDT(ve.d%>%
#                                  select(Phenotype, vars)), 
#                        id.vars = c("Phenotype"), 
#                        variable.name = "var"))
# plotme$value<-as.double(plotme$value)
# plotme$ci<-NA
# plotme$ci[plotme$var=="SNP-based heritability (h2; LDSC)"]<-ve.d$h2ci
# plotme$Phenotype<-factor(plotme$Phenotype, levels=ph)
# 
# 
# ggplot(plotme, aes(fill=var, y=value, x=Phenotype)) + 
#     geom_bar(stat = "identity", position = position_dodge(width = 1))+
#     geom_errorbar(aes(ymin = value - ci, ymax = value + ci),
#                   position = position_dodge(width = 1),
#                   na.rm=T, size=0.5,
#                   width=0.5)+
#     ylab("Explained phenotypic variance (%)")+
#     theme_bw()+
#     labs(fill = "Variant grouping")+
#     theme(axis.text.x = element_text(angle = 45, 
#                                      vjust = 1, hjust=1))
# 





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

# comparing h2 distribution of disc and ma ----------------------
h2<-c(0.1296, 0.1175, 0.1538, 0.1601,0.1636)
h2d<-c(0.148, 0.1724, 0.1325, 0.1722, 0.1849)
h2d

h2d-h2

h2.meta<-h2
h2.disc<-h2d

h2.meta
h2.disc
h2.meta-h2.disc


