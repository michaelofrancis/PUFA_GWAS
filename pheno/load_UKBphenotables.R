library(tidyverse)

setwd("/scratch/mf91122/UKB-pheno")

#Load UK Biobank datasets-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
source('ukb34137_loaddata.r') #15 min
bd_add<-read.table("updated_42606/ukb42606.tab",
                   header=TRUE, sep="\t")
bd4<-read.table("43083/ukb43083.tab",
                header=TRUE, sep="\t")

bd_join<-as_tibble(inner_join(bd4, bd_add,
                              by=intersect(colnames(bd4), colnames(bd_add))))

rm(bd_add, bd4)

bd5<-read.table("43770/ukb43770.tab",
                header=TRUE, sep="\t")
bd_join2<-as_tibble(inner_join(bd5, bd_join,
                               by=intersect(colnames(bd5), colnames(bd_join))))
rm(bd5, bd_join)


bd6<-read.table("45916/ukb45916.tab",
                header=TRUE, sep="\t")

bd_join3<-as_tibble(inner_join(bd6, bd_join2,
                               by=intersect(colnames(bd6), colnames(bd_join2))))
rm(bd6, bd_join2)

bd7<-read.table("47434/ukb47434.tab",
                header=TRUE, sep="\t")

bd_join4<-as_tibble(inner_join(bd7, bd_join3,
                               by=intersect(colnames(bd7), colnames(bd_join3))))

rm(bd_join3)

bd8<-read.table("48364/ukb48364.tab",
                header=TRUE, sep="\t") #updated NMR dataset
