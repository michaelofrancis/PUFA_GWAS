library(plyr)
library(tidyverse)
library(gridExtra)
library(ggpattern)
library(cowplot)

setwd("~/Documents/Research/PUFA-GWAS/Fig5")
adjloglabel<-expression("-log"[10]*"(adj"*italic(P)*")")

pheno<-c("w3", "w6", "DHA", "LA", "MUFA")

phenotypesb<-c("Omega-3 fatty acids", 
               "Omega-6 fatty acids", 
               "Docosahexaenoic acid", 
               "Linoleic acid", 
               "Monounsaturated fatty acids")



# Create novel/nonnovel tables----------------------------------------

#prep novel
novel<-list()
for (i in 1:length(pheno)){
novel[[i]]<-read_delim(paste(
    "/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/novelty/send_novel_back_to_FUMA_pathway_analysis/results/", 
    pheno[i], "/GS.txt", sep=""))
    novel[[i]]$Phenotype<-pheno[i]
    novel[[i]]<-novel[[i]]%>%select(Phenotype, everything())
}
novel<-do.call(rbind, novel)
novel$novel_GS<-1

#Prep non-novel
nonnovel<-list()
for (i in 1:length(pheno)){
    nonnovel[[i]]<-read_delim(paste(
        "/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/novelty/send_known_back_to_FUMA_pathway_analysis/results/", 
        pheno[i], "/GS.txt", sep=""))
    nonnovel[[i]]$Phenotype<-pheno[i]
    nonnovel[[i]]<-nonnovel[[i]]%>%select(Phenotype, everything())
}
nonnovel<-do.call(rbind, nonnovel)
nonnovel$novel_GS<-0


GS<-rbind(novel, nonnovel)



# New plots 7-18-2022 --------------------------------------------
col=c("turquoise2",
    "olivedrab3",   #DHA 
    "tan2" , #w6
     "gold",  #LA
    "mediumpurple3"  #MUFA /innermost phenotype alternating colors (MUFA)
    )

#--------------------------------------------------------------
#PLOT1--------------------------------------------------------------
# Pathways unstrat -------------------------------------------------
#--------------------------------------------------------------

path<-GSall%>%filter(Category=="Wikipathways")
path<-path[rev(order(path$logadjP)),]
path$GeneSet<-gsub("_", " ", path$GeneSet)

#Fix capitalization manually
path$GeneSet<-str_to_sentence(path$GeneSet)
path$GeneSet<-gsub("Dna", "DNA", path$GeneSet)
path$GeneSet<-gsub("dna", "DNA", path$GeneSet)
path$GeneSet<-gsub("ldl, hdl and tg, including diseases", "LDL, HDL and TG", path$GeneSet)
path$GeneSet<-gsub("Fbxl10", "FBXL10", path$GeneSet)
path$GeneSet<-gsub("map/erk", "MAP/ERK", path$GeneSet)
path$GeneSet<-gsub("b-cell", "B-cell", path$GeneSet)
path$GeneSet<-gsub("hutchinson-gilford", "Hutchinson-Gilford", path$GeneSet)
path$GeneSet<-gsub("b12", "B12", path$GeneSet)
path$GeneSet<-gsub("Vitamin a", "Vitamin A", path$GeneSet)
path$GeneSet<-gsub("\\(srebp\\)", "", path$GeneSet)
path$GeneSet<-gsub("pcsk9", "PCSK9", path$GeneSet)
path$GeneSet<-gsub("ldl", "LDL", path$GeneSet)
path$GeneSet<-gsub("The effect of progerin on the involved genes in Hutchinson-Gilford progeria syndrome", 
                   "Progerin in Hutchinson-Gilford progeria syndrome", path$GeneSet)
path$GeneSet<-gsub("FBXL10 enhancement of MAP/ERK signaling in diffuse large B-cell lymphoma", 
                   "FBXL10 + MAP/ERK signaling in B-cell lymphoma", path$GeneSet)
path$GeneSet<-gsub("  ", " ", path$GeneSet)



traitorder<-path%>%select(GeneSet, logadjP)%>%group_by(GeneSet)%>%
    top_n(1, abs(logadjP))%>%select(GeneSet)
traitorder<-rev(unlist(traitorder))
path$GeneSet<-factor(path$GeneSet, levels=traitorder)
path$Phenotype<-factor(path$Phenotype, 
                       levels=rev(c("Monounsaturated fatty acids",
                                    "Linoleic acid",
                                    "Omega-6 fatty acids", 
                                    "Docosahexaenoic acid",
                                    "Omega-3 fatty acids")))

path
keep<-path$GeneSet[path$logadjP>1.612]
keep<-keep[!keep %in% c("Proprotein convertase subtilisin/kexin type 9 (PCSK9) mediated LDL receptor degradation",
             "Metabolism of alpha-linolenic acid",
             "Evolocumab mechanism")]

path2<-path[path$GeneSet %in% keep,]


#tiff("WikiPathAll.tiff", width = 9, height = 6, units = 'in', res = 300)
p1<-ggplot(data = path2, aes(x = GeneSet, y=logadjP, fill=Phenotype))+
    geom_point(shape=21, position =position_dodge(width=posdog),
               size=2.5)+
    coord_flip()+
    scale_fill_manual(values=col)+
    theme_bw()+
    ylab(adjloglabel)+
    geom_hline(yintercept=(-log10(0.05/2.05)),
               linetype="dashed", color = "grey")+
    labs(fill="Phenotype color") 
p1  
#dev.off()

#--------------------------------------------------------------
#PLOT2--------------------------------------------------------------
#GWAS Catalog Not stratified by novelty---------------------------------
#--------------------------------------------------------------


GSall<-as_tibble(read.csv("/Users/mike/Documents/Research/PUFA-GWAS/meta-analysis/GENE2FUNC/combined.csv"))
GSall$logadjP=-log10(GSall$adjP)
GSall<-GSall%>%select(-link, -X)




GSall$Phenotype<-mapvalues(x = GSall$Phenotype, from = unique(GSall$Phenotype),
          to=c("Docosahexaenoic acid",
               "Linoleic acid",
               "Monounsaturated fatty acids",
               "Omega-3 fatty acids",      
                                 "Omega-6 fatty acids"))


GSall
gwas<-GSall%>%filter(Category=="GWAScatalog")
gwas<-gwas[rev(order(gwas$logadjP)),]

keepgwas3<-gwas$GeneSet[c(1,2,7,12,23,24,50,76,91,95,116,140,
                          144,150,159,166,173,177)]

keepgwas3
keepgwas3<-gsub("Cholesterol, total", "Cholesterol", keepgwas3)
keepgwas3<-gsub("Sarcoidosis \\(Lofgren's syndrome vs non-Lofgren's syndrome\\)", 
                "Sarcoidosis", keepgwas3)
keepgwas3<-gsub("Cerebrospinal AB1-42 levels in mild cognitive impairment", 
                "AB1-42 levels in cognitive impairment", keepgwas3)

keepgwas3

gwas3<-gwas[gwas$GeneSet %in% keepgwas3,]

traitorder<-gwas3%>%select(GeneSet, logadjP)%>%group_by(GeneSet)%>%
    top_n(1, abs(logadjP))%>%select(GeneSet)
traitorder<-rev(unlist(traitorder))

gwas3$GeneSet<-factor(gwas3$GeneSet, levels=traitorder)
gwas3$Phenotype<-factor(gwas3$Phenotype, 
                        levels=rev(c("Monounsaturated fatty acids",
                                 "Linoleic acid",
                                 "Omega-6 fatty acids", 
                                 "Docosahexaenoic acid",
                                 "Omega-3 fatty acids")))



#tiff("GWASCatAll.tiff", width = 8, height = 6, units = 'in', res = 300)
p2<-ggplot(data = gwas3, aes(x = GeneSet, y=logadjP, fill=Phenotype))+
    geom_point(shape=21, 
               position =position_dodge(width=posdog),
               size=2.5)+
    scale_fill_manual(values=col)+
    ylab(expression("-log"[10]*"(adjP)"))+
    coord_flip()+
    theme_bw()+
       geom_hline(yintercept=(-log10(0.05/2.05)),
           linetype="dashed", color = "grey")
 p2   
#dev.off()    






# GWASCAT NOVEL ----------------------------------------------

#table(GS$Category)

#Stratified by novelty:
GS$logadjP=-log10(GS$adjP)
gwn<-GS%>%filter(Category=="GWAScatalog", novel_GS==1)
gwn$logadjP=-log10(gwn$adjP)
gwn<-gwn[rev(order(gwn$logadjP)),]
gwn<-gwn%>%select(-link)
gwn$Phenotype<-mapvalues(x = gwn$Phenotype, from = sort(unique(gwn$Phenotype)),
                           to=c("Docosahexaenoic acid",
                                "Linoleic acid",
                                "Monounsaturated fatty acids",
                                "Omega-3 fatty acids",      
                                "Omega-6 fatty acids"))


gwn

gwn$GeneSet<-gsub("Total cholesterol levels", 
                  "Total cholesterol", gwn$GeneSet)
gwn$GeneSet<-gsub("Handedness \\(Right-handed vs. non-right-handed\\)", 
                  "Handedness", gwn$GeneSet)
gwn$GeneSet<-gsub("Sarcoidosis \\(Lofgren's syndrome vs non-Lofgren's syndrome\\)", 
                "Sarcoidosis", gwn$GeneSet)

keepgwn3<-gwn$GeneSet[c(2,5,6,11,16,17,18,21,33,35,
                        40,50,54,58,60)]

gwn3<-gwn[gwn$GeneSet %in% keepgwn3,]

traitorder<-gwn3%>%select(GeneSet, logadjP)%>%group_by(GeneSet)%>%
    top_n(1, abs(logadjP))%>%select(GeneSet)
traitorder<-rev(unlist(traitorder))
gwn3$GeneSet<-factor(gwn3$GeneSet, levels=traitorder)
gwn3$Phenotype<-factor(gwn3$Phenotype, 
                        levels=rev(c("Monounsaturated fatty acids",
                                     "Linoleic acid",
                                     "Omega-6 fatty acids", 
                                     "Docosahexaenoic acid",
                                     "Omega-3 fatty acids")))


#tiff("GWASCatalogNovel.tiff", width = 9, height = 6, units = 'in', res = 300)

p3<-ggplot(data = gwn3, aes(x = GeneSet, y=logadjP, fill=Phenotype))+
    geom_point(shape=21, position =position_dodge(width=posdog),
               size=2.5)+
    scale_fill_manual(values=col)+
    ylab(expression("-log"[10]*"(adjP)"))+
    coord_flip()+
    theme_bw()+
    geom_hline(yintercept=(-log10(0.05/2.05)),
               linetype="dashed", color = "grey")
p3
#dev.off()


#--------------------------------------------
# Alcohol use disorder ------------------------------------------

# join[join$GeneSet=="Alcohol use disorder (total score)",][1:5]
# 
# novelgenes<-novel[novel$GeneSet=="Alcohol use disorder (total score)",]
# ng<-novelgenes$genes
# ng2<-unique(unlist(str_split(ng, ":")))
# ng2[1:12]
# 
# #From gene ontology:
# #All of these genes are with intracellular membrane-bounded organelle
# 
# joinb<-full_join(novel2, nonnovel2)%>%select(-novel.link, -nonnovel.link)
# View(joinb)
# 
# write.csv(joinb, "alcohol_GeneSets_novelvsnon_liver.csv", quote=F, row.names=F)
# 
# 
# joinc<-joinb[complete.cases(joinb),]%>%filter(Phenotype!="MUFA")
# joind<-joinc[rep(seq_len(nrow(joinc)), each = 2), ]
# write.csv(joind, "alcohol_GeneSets_novelvsnon_joind_liver.csv", quote=F, row.names=F)


#Do operation in Excel----


###--------------------------------------------------------------
#MAKE ALCOHOL PLOT- P4-------------------------------------------
###--------------------------------------------------------------
joindplot<-read.csv("/Users/mike/Documents/R_files/PUFA_GWAS/alcohol_GeneSets_novelvsnon_joind_plot_v4.csv",header=T)
#joindplot<-joindplot[1:8,] #ALC USE DISORDER ONLY
joindplot$logP<-(-log10(joindplot$P))
joindplot$logadjP<-(-log10(joindplot$adjP))
joindplot$Novelty<-rep(c("Novel", "Known"), 5)
joindplot$Phenotype<-str_split(joindplot$Phenotype, "\\.", simplify=T)[,1]
phenos<-c("Omega-3 fatty acids","Docosahexaenoic acid",
          "Omega-6 fatty acids","Linoleic acid", "Monounsaturated fatty acids")

joindplot$Phenotype<-mapvalues(joindplot$Phenotype, from=c("w3", "DHA", "w6", "LA", "MUFA"),
                               to=phenos)


joindplot$Phenotype<-factor(joindplot$Phenotype, levels=rev(phenos))

joindplot$Phenotype
#joindplot$GeneSet<-wrap.labels(joindplot$GeneSet, 10)


#PLOT ALCOHOL P4-----------------------==========================

# Alcohol by novelty dot plot -------------------------------------
tiff("AUDbyNovelty.tiff", width = 7, height = 5, units = 'in', res = 300)

p4<-ggplot(joindplot, mapping =aes(x=Phenotype,
                                   y=logadjP, shape=Novelty, 
                               fill=Phenotype)) +
    geom_point(size=3, stroke=0.7)+
    scale_shape_manual(values = c(22,24))+
    coord_flip()+
    theme_bw()+
    labs(y = adjloglabel, shape="Novelty of variants mapped to gene set")+
    scale_fill_manual(values=rev(col), guide="none")+
    geom_hline(yintercept=(-log10(0.05/2.05)),
               linetype="dashed", color = "grey")

p4
dev.off()



# dummy plot P5 with novelty legend I need for cowplot ---------------------

joindplot<-read.csv("/Users/mike/Documents/R_files/PUFA_GWAS/alcohol_GeneSets_novelvsnon_joind_plot_v4.csv",header=T)
#joindplot<-joindplot[1:8,] #ALC USE DISORDER ONLY
joindplot$logP<-(-log10(joindplot$P))
joindplot$logadjP<-(-log10(joindplot$adjP))
joindplot$Novelty<-c(rep(c("All variants", "Novel only", "Known only"), 3), "Novel only")

joindplot$Phenotype<-str_split(joindplot$Phenotype, "\\.", simplify=T)[,1]
phenos<-c("Omega-3 fatty acids","Docosahexaenoic acid",
          "Omega-6 fatty acids","Linoleic acid", "Monounsaturated fatty acids")

joindplot$Phenotype<-mapvalues(joindplot$Phenotype, from=c("w3", "DHA", "w6", "LA", "MUFA"),
                               to=phenos)


joindplot$Phenotype<-factor(joindplot$Phenotype, levels=rev(phenos))

joindplot$Phenotype
#joindplot$GeneSet<-wrap.labels(joindplot$GeneSet, 10)


#DUMMY PLOT P5-----------------------==========================


p5<-ggplot(joindplot, mapping =aes(x=Phenotype,
                                   y=logadjP, shape=Novelty, 
                                   fill=Phenotype)) +
    geom_point(size=3, stroke=0.7)+
    scale_shape_manual(values = c(21,22,24))+
    coord_flip()+
    theme_bw()+
    labs(y = adjloglabel, shape="Novelty of variants 
mapped to genes in 
gene set")+
    scale_fill_manual(values=rev(col), guide="none")+
    geom_hline(yintercept=(-log10(0.05/2.05)),
               linetype="dashed", color = "grey")
    


p5

# cowplot -------------------------------------------------------

wrapwidth=50
posdog=0.08

legend <- get_legend(
    # create some space to the left of the legend
    p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

legend2 <- get_legend(
    #this second legend was generated separately and added in PS
    # create some space to the left of the legend
    p5 + theme(legend.box.margin = margin(0, 0, 0, 12))
)


pgrid<-cowplot::plot_grid(p1+ theme(legend.position="none")+
                       scale_x_discrete(labels = function(x) str_wrap(x, width = wrapwidth)),
                   NULL,
                       p2+theme(legend.position="none")+
                       scale_x_discrete(labels = function(x) str_wrap(x, width = wrapwidth)),
                   p3+theme(legend.position="none")+
                       scale_x_discrete(labels = function(x) str_wrap(x, width =wrapwidth)),
                   NULL,
                   p4+theme(legend.position="none"), 
                   labels = c("A", "", "B", "C", "", "D"), 
                   ncol=3, rel_widths = c(1,0.05, 1))

tiff("Cowplotlegend2.tiff", width = 12, height = 8, units = 'in', res = 300)
plot_grid(pgrid, NULL, legend2, NULL,ncol = 4, rel_widths = c(1,0.05, .1,0.05))
dev.off()


stop()

# GO_bp unstrat -------------------------------------------------

bp<-GSall%>%filter(Category=="GO_bp")
bp<-bp[rev(order(bp$logadjP)),]
bp$GeneSet<-gsub("_", " ", bp$GeneSet)
bp$GeneSet<-str_split(bp$GeneSet, pattern = "GO ", simplify=T)[,2]
bp$GeneSet<-str_to_sentence(bp$GeneSet)
bp$GeneSet<-gsub("Dna", "DNA", bp$GeneSet)
bp$GeneSet<-gsub("dna", "DNA", bp$GeneSet)
traitorder<-bp%>%select(GeneSet, logadjP)%>%group_by(GeneSet)%>%
    top_n(1, abs(logadjP))%>%select(GeneSet)
traitorder<-rev(unlist(traitorder))
bp$GeneSet<-factor(bp$GeneSet, levels=traitorder)
bp$Phenotype<-factor(bp$Phenotype, 
                     levels=rev(c("Monounsaturated fatty acids",
                                  "Linoleic acid",
                                  "Omega-6 fatty acids", 
                                  "Docosahexaenoic acid",
                                  "Omega-3 fatty acids")))

keep<-bp$GeneSet[bp$logadjP>10]
bp2<-bp[bp$GeneSet %in% keep,]


ggplot(data = bp2, aes(x = GeneSet, y=logadjP, fill=Phenotype))+
    geom_point(shape=21, position =position_dodge(width=0.1))+
    coord_flip()+
    scale_fill_manual(values=col)+
    theme_bw()+
    ylab(adjloglabel)+
    ggtitle("GO_bp: all significant variants")+
    geom_hline(yintercept=(-log10(0.05/2.05)),
               linetype="dashed", color = "grey")




# GO_cc unstrat -------------------------------------------------

cc<-GSall%>%filter(Category=="GO_cc")
cc<-cc[rev(order(cc$logadjP)),]
cc$GeneSet<-gsub("_", " ", cc$GeneSet)
cc$GeneSet<-str_split(cc$GeneSet, pattern = "GO ", simplify=T)[,2]
cc$GeneSet<-str_to_sentence(cc$GeneSet)
cc$GeneSet<-gsub("Dna", "DNA", cc$GeneSet)
cc$GeneSet<-gsub("dna", "DNA", cc$GeneSet)
traitorder<-cc%>%select(GeneSet, logadjP)%>%group_by(GeneSet)%>%
    top_n(1, abs(logadjP))%>%select(GeneSet)
traitorder<-rev(unlist(traitorder))
cc$GeneSet<-factor(cc$GeneSet, levels=traitorder)
cc$Phenotype<-factor(cc$Phenotype, 
                     levels=rev(c("Monounsaturated fatty acids",
                                  "Linoleic acid",
                                  "Omega-6 fatty acids", 
                                  "Docosahexaenoic acid",
                                  "Omega-3 fatty acids")))





keep<-cc$GeneSet[cc$logadjP>5]
cc2<-cc[cc$GeneSet %in% keep,]
ggplot(data = cc2, aes(x = GeneSet, y=logadjP, fill=Phenotype))+
    geom_point(shape=21, position =position_dodge(width=0.1))+
    coord_flip()+
    scale_fill_manual(values=col)+
    theme_bw()+
    ylab(adjloglabel)+
    ggtitle("GO_cc: all significant variants")



# GO_mf unstrat -------------------------------------------------

mf<-GSall%>%filter(Category=="GO_mf")
mf<-mf[rev(order(mf$logadjP)),]
mf$GeneSet<-gsub("_", " ", mf$GeneSet)
mf$GeneSet<-str_split(mf$GeneSet, pattern = "GO ", simplify=T)[,2]
mf$GeneSet<-str_to_sentence(mf$GeneSet)
mf$GeneSet<-gsub("Dna", "DNA", mf$GeneSet)
mf$GeneSet<-gsub("dna", "DNA", mf$GeneSet)
traitorder<-mf%>%select(GeneSet, logadjP)%>%group_by(GeneSet)%>%
    top_n(1, abs(logadjP))%>%select(GeneSet)
traitorder<-rev(unlist(traitorder))
mf$GeneSet<-factor(mf$GeneSet, levels=traitorder)
mf$Phenotype<-factor(mf$Phenotype, 
                     levels=rev(c("Monounsaturated fatty acids",
                                  "Linoleic acid",
                                  "Omega-6 fatty acids", 
                                  "Docosahexaenoic acid",
                                  "Omega-3 fatty acids")))




mf
keep<-mf$GeneSet[mf$logadjP>2]
mf2<-mf[mf$GeneSet %in% keep,]
ggplot(data = mf2, aes(x = GeneSet, y=logadjP, fill=Phenotype))+
    geom_point(shape=21, position =position_dodge(width=0.1))+
    coord_flip()+
    scale_fill_manual(values=col)+
    theme_bw()+
    ylab(adjloglabel)+
    ggtitle("GO_mf: all significant variants")
